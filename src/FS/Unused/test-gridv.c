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
#include <string.h>
#include <math.h>
#include "mat_op3.h"
#include "ins.h"



#define SYM_TENSOR 1



#undef OUTPUT_DIR
#define OUTPUT_DIR "./output/"
int dim = 3;
SURF_BAS *surf_bas = NULL;

#define DIM_RETRUN(d) if (d >= dim) return;

/* ---------------------------------
 * Test case 1: sinoid
 *           2: P2
 *           3: P2(z)
 *           4: zero
 *           5: periodic
 *           6: \laplace u
 *           7: e(u) : e(v), P2
 *           8: e(u) : e(v), P3
 *
 * */
#undef TEST_CASE
#define TEST_CASE 8


FLOAT _L_ = 5.;
FLOAT _alpha_ = 0.;
static FLOAT a = .5;
void bench_grid(GRID *g);

void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *u)
{
#if TEST_CASE == 1
    *u = Cos(2. * M_PI * x) * Cos(2. * M_PI * y) * Cos(2. * M_PI * z);
#elif TEST_CASE == 2
    *(u++) = x*(1-x) + y*z; DIM_RETRUN(1);
    *(u++) = x*(1-x) + y*z; DIM_RETRUN(2);
    *(u++) = x*(1-x) + y*z; DIM_RETRUN(3);
#elif TEST_CASE == 3
    *u = z*(2-z);
#elif TEST_CASE == 4
    *u = 0.;
#elif TEST_CASE == 5
    const FLOAT L = _L_;
    const FLOAT alpha = _alpha_ / 180. * M_PI;
    const FLOAT h = tan(alpha);
    const FLOAT c = 1.0;
    const FLOAT b = .3;
    *u = c + b/L*x + (b/L/h - b + b*x)*z;
#elif TEST_CASE == 6
    u[0] = x * x;
    u[1] = y * z;
    u[2] = y * y + x * z;
#elif TEST_CASE == 7
    u[0] = x * x;
    u[1] = y * z;
    u[2] = y * y + x * z;
#elif TEST_CASE == 8
u[0] = x * x * y + pow(z, 0.3e1) + y * y + 1;
u[1] = y * z + 0.1e1 + 1.5;
u[2] = x * y * z + x * z * z + .5;
#else
# error Test case error!!!
#endif
}

void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *f)
{
    func_u(x, y, z, f);
#if TEST_CASE == 1
    *f = 12. * M_PI * M_PI * *f + a * *f;
#elif TEST_CASE == 2
    f[0] = 2 + a * *f; DIM_RETRUN(1);
    f[1] = 2 + a * *f; DIM_RETRUN(2);
    f[2] = 2 + a * *f; DIM_RETRUN(3);
#elif TEST_CASE == 3
    *f = 2 + a * *f;
#elif TEST_CASE == 4
    *f = 1.;
#elif TEST_CASE == 5
    *f = a * *f;
#elif TEST_CASE == 6
    f[0] = -2; 
    f[1] = 0; 
    f[2] = -2; 
#elif TEST_CASE == 7
    f[0] = -5./2.; 
    f[1] = 0; 
    f[2] = -3./2.; 
#elif TEST_CASE == 8
#  warning test case 8
#  if SYM_TENSOR
    f[0] = -0.1e1 - 0.5e1 / 0.2e1 * y - (double) (4 * z);
    f[1] = -0.3e1 / 0.2e1 * x;
    f[2] = -0.1e1 / 0.2e1 - 0.2e1 * x;
#  else
    f[0] = -2 - 2 * y - 6 * z;
    f[1] = 0;
    f[2] = -2 * x;
#  endif
#else
# error Test case error!!!
#endif
}


void func_g1(FLOAT x, FLOAT y, FLOAT z, FLOAT *g1) {
    g1[0] = 2 * x * y;
    g1[1] = x * x + 2 * y;
    g1[2] = 3 * z * z;
}

void func_g2(FLOAT x, FLOAT y, FLOAT z, FLOAT *g2) {
    g2[0] = 0;
    g2[1] = z;
    g2[2] = y;
}

void func_g3(FLOAT x, FLOAT y, FLOAT z, FLOAT *g3) {
    g3[0] = (double) y * (double) z + (double) z * (double) z;
    g3[1] = x * z;
    g3[2] = x * y + 2 * x * z;
}


FLOAT *
get_gbas_product(const FLOAT *gi, const FLOAT *gj)
{
    static FLOAT prod[Dim][Dim];

    prod[0][0] = gi[0] * gj[0] + .5 * (gi[1] * gj[1] + gi[2] * gj[2]);
    prod[0][1] = .5 * gi[0] * gj[1];
    prod[0][2] = .5 * gi[0] * gj[2];

    prod[1][0] = .5 * gi[1] * gj[0];
    prod[1][1] = gi[1] * gj[1] + .5 * (gi[0] * gj[0] + gi[2] * gj[2]);
    prod[1][2] = .5 * gi[1] * gj[2];

    prod[2][0] = .5 * gi[2] * gj[0];
    prod[2][1] = .5 * gi[2] * gj[1];
    prod[2][2] = gi[2] * gj[2] + .5 * (gi[0] * gj[0] + gi[1] * gj[1]);

    return prod[0];
}

#define INNER_PRODUCT(p, q)			\
    (*(p    ) * *(q    ) +			\
     *(p + 1) * *(q + 1) +			\
     *(p + 2) * *(q + 2))


static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h, DOF *gradu)
{
    int i, j, k, l;
    GRID *g = u_h->g;
    SIMPLEX *e;

    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);

    Unused(l);
    ForAllElements(g, e) {
	int N = DofGetNBas(u_h, e);	/* number of basises in the element */
	QUAD *quad;
	FLOAT A[N][dim][N][dim], rhs[N][dim], buffer[N], tmp[Dim];
	INT I[N][dim], J[dim][N], q;
	const FLOAT *w, *p, *gu;
	FLOAT vol, vf[dim], det;

	Unused(det);
	bzero(A, sizeof(A));
	bzero(rhs, sizeof(rhs));
	bzero(buffer, sizeof(buffer));
	vol = phgGeomGetVolume(g, e);
	quad = phgQuadGetQuad3D(2*DofTypeOrder(u_h, e));
	gu = phgQuadGetDofValues(e, gradu, quad);
	
	for (i = 0; i < N; i++)
	    for (k = 0; k < dim; k++)
		I[i][k] = J[k][i] = phgSolverMapE2L(solver, 0, e, i * dim + k);

	p = quad->points;
 	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    FLOAT x, y, z;
	    

	    //printf("  quad: %d\n", q);
	    phgGeomLambda2XYZ(g, e, p, &x, &y, &z);
	    func_f(x, y, z, vf);

	    for (i = 0; i < N; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, u_h, i, quad) + q;       /* phi_i */
		const FLOAT *ggi = phgQuadGetBasisGradient(e, u_h, i, quad) + q*Dim;    /* grad phi_i */

		for (j = 0; j < N; j++) {
		    const FLOAT *gj = phgQuadGetBasisValues(e, u_h, j, quad) + q;       /* phi_j */
		    const FLOAT *ggj = phgQuadGetBasisGradient(e, u_h, j, quad) + q*Dim;    /* grad phi_j */
	
#if SYM_TENSOR
#  warning tensor product
		    const FLOAT *tp = get_gbas_product(ggi, ggj);
		    //SHOW_M(tp, dim, dim);
		    for (k = 0; k < dim; k++) 
			for (l = 0; l < dim; l++)
			    A[j][l][i][k] += vol*(*w) * tp[k + l*dim];
#else
		    for (k = 0; k < dim; k++) 
			A[j][k][i][k] += vol*(*w) * ( ggi[0] * ggj[0]
						      + ggi[1] * ggj[1]
						      + ggi[2] * ggj[2]);
#endif
		}

		/* for (k = 0; k < dim; k++) */
		/*     rhs[i][k] += vol*(*w) * vf[k] * (*gi); */

		FLOAT eu[Dim*Dim];
		MAT3_SYM(gu, eu);
#if SYM_TENSOR
		for (k = 0; k < dim; k++) {
		    rhs[i][k] -= vol*(*w) * (INNER_PRODUCT(eu + k*dim, ggi) - vf[k]* (*gi));
		}
#else
		for (k = 0; k < dim; k++) {
		    rhs[i][k] -= vol*(*w) * (INNER_PRODUCT(gu + k*dim, ggi) - vf[k]* (*gi));
		}
#endif
	    }
	    w++;
	    p += Dim+1;
	    gu += Dim*Dim;
	}



	/****************/
	/* Outflow bdry */
	/****************/
	int s;
	DOF_USER_FUNC func_gn[] = {func_g1, func_g2, func_g3};
	for (s = 0; s < NFace; s++) 
	    if (e->bound_type[s] & NEUMANN) {
		int v0, v1, v2;
		int nbas_face = NbasFace(u_h);
		SHORT bases[nbas_face];
		FLOAT lambda[Dim + 1], gn[DDim], eu[DDim];
		int order = DofTypeOrder(u_h, e) * 3 - 1;

		phgDofGetBasesOnFace(u_h, e, s, bases);
		v0 = GetFaceVertex(s, 0);
		v1 = GetFaceVertex(s, 1);
		v2 = GetFaceVertex(s, 2);
		lambda[s] = 0.;

		FLOAT area = phgGeomGetFaceArea(g, e, s);
		const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);
		QUAD *quad = phgQuadGetQuad2D(order);
		
		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    FLOAT vu[Dim];
		    FLOAT x, y, z;
		    lambda[v0] = *(p++);
		    lambda[v1] = *(p++);
		    lambda[v2] = *(p++);
		
		    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		    phgDofEval(u_h, e, lambda, vu);
		    for (i = 0; i < nbas_face; i++) {
			int ii = bases[i];
			const FLOAT *gi_u  = u_h->type->BasFuncs(u_h, e, ii, ii + 1, lambda);
			const FLOAT *ggi_u = u_h->type->BasGrads(u_h, e, ii, ii + 1, lambda);

			for (k = 0; k < Dim; k++) 
			    func_gn[k](x, y, z, &gn[k*Dim]);
			MAT3_SYM(gn, eu);

#if SYM_TENSOR
			for (k = 0; k < Dim; k++) {
			    rhs[ii][k] += area*(*w) * INNER_PRODUCT(eu + k*Dim, normal) * (*gi_u);
			}
#else
			for (k = 0; k < Dim; k++) {
			    rhs[ii][k] += area*(*w) * INNER_PRODUCT(gn + k*Dim, normal) * (*gi_u);
			}
#endif
		    }     /* end of bas_i */
		    w++;
		}		/* end of quad point */
	    }		/* end of face outflow */


	/* Rotate bases */
	for (i = 0; i < N; i++) {
	    INT id = phgDofMapE2D(surf_dof, e, i * (Dim*Dim)) / (Dim*Dim);
	    if (!rotated[id])
		continue;	
	    const FLOAT *trans = Trans + i*(Dim*Dim);

	    trans_left  (&A[i][0][0][0], Dim*N, Dim*N, trans);
	    trans_rightT(&A[0][0][i][0], Dim*N, Dim*N, trans);
	    
	    trans_left  (&rhs[i][0], 1, 1, trans);
	}



	for (i = 0; i < N; i++) {
	    //SHOW_V(rhs[i], Dim);
	    for (k = 0; k < Dim; k++) {
		if (phgDofDirichletBC_(u_h, e, i*Dim+k, NULL, buffer, tmp,
				       DOF_PROJ_NONE)) {
		    rhs[i][k] = 0.;
		    phgMatAddEntries(solver->mat, 1, I[i] + k, N, J[k], buffer);
		} else {
		    phgMatAddEntries(solver->mat, 1, I[i] + k, N*Dim, I[0], &A[i][k][0][0]);
		}
	    }
	    //SHOW_V(rhs[i], Dim);
	}

	phgSolverAddRHSEntries(solver, N*dim, I[0], rhs[0]);
    }

    phgMatDumpMATLAB(solver->mat, "A", "A_.m");
    phgVecDumpMATLAB(solver->rhs, "b", "b_.m");
}


#if 0
void dof_set_normal_data(DOF *u_h, SURF_BAS *surf_bas)
{
    GRID *g = u_h->g;
    SIMPLEX *e;

    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);

    ForAllElements(g, e) {
	int i, N = surf_dof->type->nbas;
	FLOAT vu[N][Dim], u0[Dim];

	phgDofGetElementDatas(u_h, e, vu[0]);
	for (i = 0; i < N; i++) {
	    INT id = phgDofMapE2D(surf_dof, e, i * (Dim*Dim)) / (Dim*Dim);
	    if (!rotated[id])
		continue;	

	    const FLOAT *trans = Trans + id*(Dim*Dim);
	    FLOAT *coord = phgDofGetElementCoordinates(u_h, e, i*Dim);
	    func_u(coord[0], coord[1], coord[2], u0);
	    trans_left  (vu[i], 1, 1, trans);
	    vu[i][0] = INNER_PRODUCT(u0, Trans);
	    trans_leftT (vu[i], 1, 1, trans);
	}
	phgDofSetElementDatas(u_h, e, vu[0]);
    }	      /* end elem */

    return;
}
#endif


int
my_bc_map(int bctype)
{
    switch (bctype) {
#if 1
    case 1:
	return DIRICHLET;
    case 2:
	return NEUMANN;
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
    DOF *u_h, *f_h, *grad_u, *error, *u, *gradu, *du;
    SOLVER *solver;
    FLOAT L2error;
    size_t mem, mem_peak;

    
    Unused(grad_u);
    Unused(e);
    phgOptionsRegisterFloat("a", "Coefficient", &a);
    phgOptionsRegisterFloat("L", "Length Scale", &_Length_);
    phgOptionsRegisterFloat("alpha", "Slope angle", &_alpha_);
    phgOptionsRegisterFloat("tol", "Tolerance", &tol);
    phgOptionsRegisterInt("mem_max", "Maximum memory (MB)", &mem_max);
    phgOptionsRegisterInt("periodicity", "Set periodicity", &periodicity);
    phgOptionsRegisterInt("dof_dim", "dof dim", &dim);
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


    assert(dim == Dim);
    phgPrintf("---------------------\n");
    phgPrintf("Length scale: %E\n", _L_);
    phgPrintf("Slop angle  : %E\n", _alpha_);
    phgPrintf("---------------------\n");



    /* ------------------------------------------------------------
     * 
     *    change grid coord
     * 
     * ------------------------------------------------------------ */
    //ice_grid(g);



    u_h = phgDofNew(g, DOF_DEFAULT, dim, "u_h", DofInterpolation);
    du = phgDofNew(g, DOF_DEFAULT, dim, "du", DofNoAction);
    phgDofSetDataByValue(du, 0.0);
    {
#if 1
	BTYPE DB_masks[3] = {DIRICHLET|NEUMANN, DIRICHLET, DIRICHLET};
#else
	BTYPE DB_masks[3] = {DIRICHLET, DIRICHLET, DIRICHLET};
#endif
	phgDofSetDirichletBoundaryMasks(u_h, DB_masks);
    }
    //phgDofSetDataByFunction(u_h, func_u);

    f_h = phgDofNew(g, DOF_DEFAULT, dim, "f_h",  func_f);
    u = phgDofNew(g, DOF_ANALYTIC, dim, "u", func_u);

    double t0 = phgGetTime(NULL);
    
    if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	phgPrintf("Repartition mesh, load imbalance: %lg\n",
		  (double)g->lif);
    phgPrintf("%d DOF, %d elements, %d submeshes, load imbalance: %lg\n",
	      DofGetDataCountGlobal(u_h), g->nleaf_global, g->nprocs,
	      (double)g->lif);
    surf_bas = get_surface_bases(g, u_h->type);
#if 1
    /* Case 1. set inflow & outflow normal,
     *   note the order. */
    dof_set_normal_data(u_h, surf_bas);
    phgDofSetBdryDataByFunction(u_h, func_u, DIRICHLET);
#elif 0
    /* Case 2. set inflow & outflow */
    phgDofSetBdryDataByFunction(u_h, func_u, DIRICHLET|NEUMANN);
#else
    /* Case 3. set all */
    phgDofSetDataByFunction(u_h, func_u);
#endif
    gradu = phgDofGradient(u_h, NULL, NULL, NULL);


    /* ------------------------------------------------------------
     * 
     *    Solver
     * 
     * ------------------------------------------------------------ */

    solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
    phgPrintf("   DOF: %d, unknowns: %d, Dirichlet bdry: %d\n",
	      DofGetDataCountGlobal(u_h), solver->rhs->map->nglobal,
	      solver->rhs->map->bdry_nglobal);
    build_linear_system(solver, u_h, f_h, gradu);
    phgSolverSolve(solver, TRUE, du, NULL);
    rotate_dof_bases(du, surf_bas, FALSE);

    phgDofAXPBY(1., du, 1., &u_h); 
    phgPrintf("   nits = %d\n", solver->nits);
    phgSolverDestroy(&solver);



    error = phgDofNew(g, DOF_DEFAULT, dim, "err", func_u);
    phgDofAXPBY(-1., u_h, 1., &error); 
    L2error = phgDofNormL2(error);
    mem = phgMemoryUsage(g, &mem_peak);
    phgPrintf("*  L2 error = %0.3le\n", (double)L2error);
    phgPrintf("   mem = %0.2lfMB\n", (double)mem_peak / (1024.0 * 1024.0));
    phgPrintf("   Wall time: %0.3le\n", phgGetTime(NULL) - t0);

    phgExportEnsight(g, OUTPUT_DIR"test-grid", u_h, error, NULL);
    phgExportVTK(g, OUTPUT_DIR"test-grid.vtk", u_h, error, NULL);

    phgDofFree(&u);
    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&error);



    phgPrintf("\nSuccesfully done.\n===============\n\n");
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
