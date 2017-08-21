/*
 *  Benchmark test of anisotropic prism grid.
 *
 *  The grid is defined by a mapping
 *   from  [0, 1] x [0, 1] x [0, 1]
 *   to    [0, 1] x [0, 1] x [h_za(x, y), h_zb(x, y)]
 *
 *  In benchmark[1] A & B, the grid in vertical direction seems really bad in z-dir,
 *  can this ugly grid work at all?
 *
 *  This simple test is trying to solve a poisson problem on this mesh and see what
 *  is happening...
 *
 *
 *  Ref:
 *  [1] Benchmark experiments for higher-order and full-Stokes ice sheet models (ISMIP–HOM).
 *
 *
 *  
 *  */

#include "phg.h"
#include "ins.h"
#include <string.h>
#include <math.h>

#undef OUTPUT_DIR
#define OUTPUT_DIR "./output/"


/* ---------------------------------
 * Test case 1: sinoid
 *           2: P2
 *           3: P2(z)
 *           4: zero
 *           5: periodic
 *
 * */
#undef TEST_CASE
#define TEST_CASE 2


FLOAT _L_ = 5.;
FLOAT _alpha_ = 0.;
static FLOAT a = .5;
void bench_grid(GRID *g);

void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
#if TEST_CASE == 1
    *value = Cos(2. * M_PI * x) * Cos(2. * M_PI * y) * Cos(2. * M_PI * z);
#elif TEST_CASE == 2
    *value = x*(1-x) + y*z;
#elif TEST_CASE == 3
    *value = z*(2-z);
#elif TEST_CASE == 4
    *value = 0.;
#elif TEST_CASE == 5
    const FLOAT L = _L_;
    const FLOAT alpha = _alpha_ / 180. * M_PI;
    const FLOAT h = tan(alpha);
    const FLOAT c = 1.0;
    const FLOAT b = .3;
    *value = c + b/L*x + (b/L/h - b + b*x)*z;
#else
# error Test case error!!!
#endif
}

void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    func_u(x, y, z, value);
#if TEST_CASE == 1
    *value = 12. * M_PI * M_PI * *value + a * *value;
#elif TEST_CASE == 2
    *value = 2 + a * *value;
#elif TEST_CASE == 3
    *value = 2 + a * *value;
#elif TEST_CASE == 4
    *value = 1.;
#elif TEST_CASE == 5
    *value = a * *value;
#else
# error Test case error!!!
#endif
}




/*
 *  Ice slab map:
 *  z(x, y) 
 *
 *
 *  */
static void
func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord)
{
    const FLOAT alpha = _alpha_ / 180. * M_PI;
    const FLOAT L = _Length_;
    FLOAT Ha = -x * tan(alpha);
#if 0
    FLOAT Hb = Ha - (1./L) + (.5/L) * sin(2*M_PI * x) * sin(2*M_PI * y);
#elif 1
    FLOAT Hb = Ha - (1./L); 
#else
    phgError(1, "Wrong ice bench case!!!\n");
#endif
    FLOAT Hz = Ha - Hb;

    FLOAT rz = z;		/* ratio z in [0, 1] */
    z = Hb + Hz * rz;

    coord[0] = x * L;
    coord[1] = y * L;
    coord[2] = z * L;
}


void ice_grid(GRID *g)
{
    SIMPLEX *e;
    INT i, k, s;
    DOF *coord = phgDofNew(g, DOF_P1, Dim, "coord", func_ice_slab);
    FLOAT *v0, *v1, Area[4];
    const FLOAT eps = 1e-8;
    
    Unused(k);

    /* ----------------------------------------
     * 
     *   Boundary type
     * 
     * ---------------------------------------- */
    bzero(Area, sizeof(Area));
    ForAllElements(g, e) {
	for (s = 0; s < NFace; s++) {
	    FLOAT a, n_top[] = {0, 0, 1}, n_bottom[] = {0, 0, -1};
	    const FLOAT *n;

#if 0
	    e->bound_type[s] &= ~DIRICHLET;
	    if (e->neighbours[s] != NULL)
	    	continue;
	    n = phgGeomGetFaceOutNormal(g, e, s);
	    a = phgGeomGetFaceArea(g, e, s);
	    if (fabs(1 - INNER_PRODUCT(n, n_top)) < eps) {
	    	e->bound_type[s] = DIRICHLET;
	    	Area[0] += a;
	    } else if (fabs(1 - INNER_PRODUCT(n, n_bottom)) < eps) {
	    	e->bound_type[s] = NEUMANN;
	    	Area[1] += a;
	    } else {
	    	e->bound_type[s] = NEUMANN;
	    	if (e->bound_type[s] & DIRICHLET)
	    	    Area[2] += a;
	    	else
	    	    Area[3] += a;
	    }
#else
	    if (e->neighbours[s] == NULL)
		e->bound_type[s] |= DIRICHLET;
#endif
	}
    }
    phgUpdateBoundaryTypes(g);

    phgPrintf("--------------------\n");
    phgPrintf("Set up boundary type\n");
    phgPrintf("  Area top   : %lf\n", Area[0]);
    phgPrintf("  Area bottom: %lf\n", Area[1]);
    phgPrintf("  Area dirich: %lf\n\n", Area[2]);
    phgPrintf("  Area other : %lf\n\n", Area[3]);


#if 1
    /* ----------------------------------------
     * 
     *   Map to physical region
     * 
     * ---------------------------------------- */
    for (i = 0; i < g->nvert; i++) {
	v1 = DofVertexData(coord, i);
	v0 = g->verts[i];
	memcpy(v0, v1, Dim*sizeof(*v0));
    }
    phgGeomInit_(g, TRUE);
#endif
    return;
}



static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
    int i, j;
    GRID *g = u_h->g;
    SIMPLEX *e;

    //vtk_verb = 1;
    assert(u_h->dim == 1);
    ForAllElements(g, e) {
	int N = DofGetNBas(u_h, e);	/* number of basises in the element */
	QUAD *quad;
	FLOAT A[N][N], rhs[N], buffer[N];
	INT I[N], q;
	const FLOAT *w, *p;
	FLOAT vol, vf, det;

	Unused(det);
	bzero(A, sizeof(A));
	bzero(rhs, sizeof(rhs));
	bzero(buffer, sizeof(buffer));
	vol = phgGeomGetVolume(g, e);
	quad = phgQuadGetQuad3D(2*DofTypeOrder(u_h, e));
	
	for (i = 0; i < N; i++)
	    I[i] = phgSolverMapE2L(solver, 0, e, i);

	p = quad->points;
 	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    FLOAT x, y, z;
	    
	    //printf("  quad: %d\n", q);
	    phgGeomLambda2XYZ(g, e, p, &x, &y, &z);
	    func_f(x, y, z, &vf);

	    for (i = 0; i < N; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, u_h, i, quad) + q;       /* phi_i */
		const FLOAT *ggi = phgQuadGetBasisGradient(e, u_h, i, quad) + q*Dim;    /* grad phi_i */

		for (j = 0; j < N; j++) {
		    const FLOAT *gj = phgQuadGetBasisValues(e, u_h, j, quad) + q;       /* phi_j */
		    const FLOAT *ggj = phgQuadGetBasisGradient(e, u_h, j, quad) + q*Dim;    /* grad phi_j */
	
		    A[j][i] += vol*(*w) * ( ggi[0] * ggj[0]
		    			   + ggi[1] * ggj[1]
		    			   + ggi[2] * ggj[2]
					   + a * (*gi) * (*gj));
		}
		rhs[i] += vol*(*w) * vf * (*gi);
	    }
	    w++;
	    p += Dim+1;
	}

	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, func_u, buffer, &rhs[i],
				  DOF_PROJ_NONE)) {
		phgMatAddEntries(solver->mat, 1, I + i, N, I, buffer);
	    } else {
		phgMatAddEntries(solver->mat, 1, I + i, N, I, &A[i][0]);
	    }
	}
	
	phgSolverAddRHSEntries(solver, N, I, rhs);
    }

}




int
my_bc_map(int bctype)
{
    switch (bctype) {
#if 0
    case 0:
	return DIRICHLET;
    case 1:
	return NEUMANN;
    case 2:
	return BDRY_USER1;
#endif
    default:
	return DIRICHLET;
    }
}

FLOAT _Length_;
FLOAT _alpha_;


int
main(int argc, char *argv[])
{
    INT periodicity = 0 /* X_MASK | Y_MASK | Z_MASK */;
    INT mem_max = 400, pre_refines = 0;
    FLOAT tol = 5e-1;
    char *fn = "../test/cube4.dat";
    GRID *g;
    SIMPLEX *e;
    DOF *u_h, *f_h, *grad_u, *error, *u;
    SOLVER *solver;
    FLOAT L2error;
    size_t mem, mem_peak;

    
    Unused(grad_u);
    Unused(e);
    phgOptionsRegisterFloat("a", "Coefficient", &a);
    phgOptionsRegisterFloat("length_scale", "Length Scale", &_Length_);
    phgOptionsRegisterFloat("slop_alpha", "Slope angle", &_alpha_);
    phgOptionsRegisterFloat("tol", "Tolerance", &tol);
    phgOptionsRegisterInt("mem_max", "Maximum memory (MB)", &mem_max);
    phgOptionsRegisterInt("periodicity", "Set periodicity", &periodicity);
    phgOptionsRegisterFilename("mesh_file", "Mesh file", &fn);
    phgOptionsRegisterInt("pre_refines", "pre refines", &pre_refines);

    phgOptionsPreset("-solver hypre "
		     "-hypre_solver pcg "
		     "-hypre_pc boomeramg "
		     "-solver_maxit 1000 "
		     "-solver_rtol 1e-12");

    phgInit(&argc, &argv);

    /* ------------------------------------------------------------
     * 
     *    Import grid
     * 
     * ------------------------------------------------------------ */
    g = phgNewGrid(-1);
    phgSetPeriodicity(g, periodicity);
    phgImportSetBdryMapFunc(my_bc_map);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgRefineAllElements(g, pre_refines);


    phgPrintf("---------------------\n");
    phgPrintf("Length scale: %E\n", _L_);
    phgPrintf("Slop angle  : %E\n", _alpha_);
    phgPrintf("---------------------\n");



    /* ------------------------------------------------------------
     * 
     *    change grid coord
     * 
     * ------------------------------------------------------------ */
    ice_grid(g);








    u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
    phgDofSetDataByValue(u_h, 0.0);

    f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h",  func_f);
    u = phgDofNew(g, DOF_ANALYTIC, 1, "u", func_u);

    double t0 = phgGetTime(NULL);
    
    if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	phgPrintf("Repartition mesh, load imbalance: %lg\n",
		  (double)g->lif);
    phgPrintf("%d DOF, %d elements, %d submeshes, load imbalance: %lg\n",
	      DofGetDataCountGlobal(u_h), g->nleaf_global, g->nprocs,
	      (double)g->lif);

    /* ------------------------------------------------------------
     * 
     *    Solver
     * 
     * ------------------------------------------------------------ */

    solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
    phgPrintf("   DOF: %d, unknowns: %d, Dirichlet bdry: %d\n",
	      DofGetDataCountGlobal(u_h), solver->rhs->map->nglobal,
	      solver->rhs->map->bdry_nglobal);
    build_linear_system(solver, u_h, f_h);
    phgSolverSolve(solver, TRUE, u_h, NULL);
    phgPrintf("   nits = %d\n", solver->nits);
    phgSolverDestroy(&solver);



    error = phgDofNew(g, DOF_DEFAULT, 1, "err", func_u);
    phgDofAXPBY(-1., u_h, 1., &error); 
    L2error = phgDofNormL2(error);
    mem = phgMemoryUsage(g, &mem_peak);
    phgPrintf("*  L2 error = %0.3le\n", (double)L2error);
    phgPrintf("   mem = %0.2lfMB\n", (double)mem_peak / (1024.0 * 1024.0));
    phgPrintf("   Wall time: %0.3le\n", phgGetTime(NULL) - t0);

    phgExportEnsight(g, OUTPUT_DIR"test-grid", u_h, error, NULL);

    phgDofFree(&u);
    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&error);



    phgPrintf("\nSuccesfully done.\n===============\n\n");
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
