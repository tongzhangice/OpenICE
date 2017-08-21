#include "phg.h"
#include <unistd.h>
#include <string.h>

#define DUMP_MAT_VEC 0
#define TEST_DOF_ORDER 2

#define Bzero(a) bzero(a, sizeof(a))
#define SQUARE(x) ((x)*(x))
#define INNER_PRODUCT(p, q)                     \
    (*(p    ) * *(q    ) +                      \
     *(p + 1) * *(q + 1) +                      \
     *(p + 2) * *(q + 2))

/*
 * Test func laplacian is generated using matlab
 * >> lap = -diff(f, x, 2) - diff(f, y, 2) - diff(f, z, 2)
 *
 *  */
static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
#if TEST_DOF_ORDER == 1    
# warning test dof order 1
    *value = 3.3 * x - 0.2445 * y + 0.922 * z - 0.34;
#elif TEST_DOF_ORDER == 2    
    *value = (3.3 * x - 0.2445 * y + 0.922 * z) * (1.- x + 3.*z) + 2.*(x-y)*(2.*y-3.*z+1.1);
#elif TEST_DOF_ORDER == 3
    *value = (3.3 * x - 0.2445 * y) * (1.8 + 0.922 * z) * (1.- x + 3.*z) + (2.*z+10.)*(x-y)*(2.*y-3.*z+1.1);
#else
    *value = Cos(2. * M_PI * x) * Cos(2. * M_PI * y) * Cos(2. * M_PI * z);
#endif
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
#if TEST_DOF_ORDER == 1    
# warning test dof order 1
    *value = 0.;
#elif TEST_DOF_ORDER == 2 
    *value = 2267./250.;
#elif TEST_DOF_ORDER == 3    
    *value = (35213.*z)/2500. - (5323713.*y)/500000. - (15639.*x)/2500. + 1297./25.;
#else
    func_u(x, y, z, value);
    *value = 12. * M_PI * M_PI * *value;
#endif
}

static int
my_bc_map(int bctype)
{
#if 1 
    return DIRICHLET;
#else
    switch (bctype) {
    case 0:
	return NEUMANN;	/* Neumann */
    case 1:
	return DIRICHLET;	/* Dirichlet */
    default:
	return DIRICHLET;
    }
#endif
}

double
elapsed_time(GRID *g, BOOLEAN flag, double mflops)
/* returns and prints elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;

    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    mem = phgMemoryUsage(g, NULL);

    if (flag) {
        if (mflops <= 0)
            phgPrintf("[%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
        else
            phgPrintf("[%0.4lgMB %0.4lfs %0.4lgGF]\n",
                      mem / (1024.0 * 1024.0), et, mflops * 1e-3);
    }

    return et;
}


int
main(int argc, char *argv[])
{
    GRID *g;
    SIMPLEX *e;
    DOF *u_h, *f_h, *err;
    MAP *map;
    SOLVER *solver;
    char *fn = "./cube.D6.mesh", vtk_file[1000];
    int pre_refines = 0;
    BOOLEAN debug = FALSE, export = FALSE;
    INT mem_max = 400;
    size_t mem, mem_peak;
    int i, j;

    Unused(mem);
    Unused(mem_peak);
    Unused(solver);
    phgOptionsRegisterNoArg("debug", "debug for gdb", &debug);
    phgOptionsRegisterNoArg("export", "export result", &export);
    phgOptionsRegisterFilename("mesh_file", "Mesh file", &fn);
    phgOptionsRegisterInt("mem_max", "Maximum memory (MB)", &mem_max);
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
    phgInit(&argc, &argv);
    phgOptionsShowUsed();

    g = phgNewGrid(-1);
    phgImportSetBdryMapFunc(my_bc_map);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    //phgDumpGrid(g);

    phgRefineAllElements(g, pre_refines);
    if (export) phgExportMedit(g, "refined.mesh");
    if (export) phgExportVTK(g, "refined.vtk", NULL);
    

#if 1
    if (phgBalanceGrid(g, 1.1, 1, NULL, 0.))
	phgPrintf("Repartition mesh, load imbalance: %lg\n",
		  (double)g->lif);
    phgExportVTK(g, "metis.vtk", NULL);
#endif

    if (debug) {
	int _i_ = 0;
	char hostname[256];
	gethostname(hostname, sizeof(hostname));
	printf("#### Lengweee debug "
	       "PID %d on %s ready for attach\n", getpid(), hostname);
	fflush(stdout);
	while (0 == _i_) {
	    sleep(50000);
	}
	MPI_Barrier(g->comm);
	printf("### PID %d, ready!\n", getpid());
    }




    u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", func_u);
    f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h", func_f);
    solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
    solver->mat->handle_bdry_eqns = TRUE;

    phgPrintf("\n%d DOF(%s), %d elements, %d verts, %d submesh%s, load imbalance: %lg\n",
	      DofGetDataCountGlobal(u_h), u_h->type->name, 
	      g->nleaf_global, g->nvert_global,
	      g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);

    
    /****************************/
    /* Build solver mat & RHS   */
    /****************************/
    ForAllElements(g, e) {
	int N = u_h->type->nbas;
	int q, order = 2;
	FLOAT A[N][N], rhs[N], buffer[N];
	INT I[N];
	QUAD *quad;
	FLOAT vol, area;
	const FLOAT *w, *p, *vf;

	Unused(area);
	Bzero(A); Bzero(rhs); Bzero(buffer);

	for (i = 0; i < N; i++)
	    I[i] = phgMapE2L(solver->mat->cmap, 0, e, i);

	vol = phgGeomGetVolume(g, e);
	quad = phgQuadGetQuad3D(order);
	vf = phgQuadGetDofValues(e, f_h, quad);

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    for (i = 0; i < N; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, u_h, i, quad) + q;       /* phi_i */
		const FLOAT *ggi = phgQuadGetBasisGradient(e, u_h, i, quad) + q*Dim;    /* grad phi_i */
		for (j = 0; j < N; j++) {
		    const FLOAT *gj = phgQuadGetBasisValues(e, u_h, j, quad) + q;       /* phi_i */
		    const FLOAT *ggj = phgQuadGetBasisGradient(e, u_h, j, quad) + q*Dim;    /* grad phi_i */
			
		    Unused(gj);
		    A[j][i] += (*w)*vol * INNER_PRODUCT(ggi, ggj);
		}
		rhs[i] += (*w)*vol * (*vf) * (*gi);
	    }
	    vf++;
	    w++; p += Dim + 1;
	}
	    
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, func_u, buffer, rhs+i,
				  DOF_PROJ_NONE)) {
		phgMatAddEntries(solver->mat, 1, I + i, N, I, buffer); 
	    }
	    else {	/* interior node */
		phgMatAddEntries(solver->mat, 1, I + i, N, I, A[i]); 
	    }
	}
	phgVecAddEntries(solver->rhs, 0, N, I, rhs);
    }

    phgMatAssemble(solver->mat);
    phgVecAssemble(solver->rhs);

    phgSolverUpdateRHS(solver);
    phgPrintf("\nBuild solver ");

    if (DUMP_MAT_VEC) {
	phgPrintf("*** Dumping An, bn\n");
	phgMatDumpMATLAB(solver->mat, "An", "An_.m");
	phgVecDumpMATLAB(solver->rhs, "bn", "bn_.m");
    }


    /* Solve */
    phgDofSetDataByValue(u_h, 0.);
    elapsed_time(g, FALSE, 0.);
    phgSolverSolve(solver, FALSE, u_h, NULL);
    phgPrintf("   nits: %d, res: %e ", solver->nits, solver->residual);
    elapsed_time(g, TRUE, 0.);
    if (export) phgExportVTK(g, "mg.vtk", u_h, NULL);

    /* Error check */
    err = phgDofNew(g, DOF_DEFAULT, 1, "error", func_u);
    phgDofAXPBY(-1., u_h, 1., &err);
    phgPrintf("Error L2: %E\n", phgDofNormL2(err));
    if (export) phgExportVTK(g, "err.vtk", u_h, err, NULL);
    phgDofFree(&err);


    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&err);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
