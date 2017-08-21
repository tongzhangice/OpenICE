/*
 * Build MG solver matrix for convection-diffusion problem
 *   on coarse grid. 
 *
 *  */
#include "ins.h"
#define _nsp (ns->ns_params)
#define _pcd (ns->pcd)
#define TENSOR_STRAIN 1

void 
build_mg_levels(NSSolver *ns, int i_slv)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    MULTI_GRID *mg = ns->mg;
    MG_LEVEL **ml = mg->ml;
    int i, j, k, l, level, max_level, min_level, start_level;
    DOF *u_L, *p_L;
    FLOAT *dt = ns->dt, Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    BOOLEAN use_GSline = mg_params->use_GSline;
    BOOLEAN use_mg_F  = ns_params->use_mg_F;
    BOOLEAN use_mg_Ap = ns_params->use_mg_Ap;
    BOOLEAN use_mg_Qp = ns_params->use_mg_Qp;

    SURF_BAS *surf_bas = ns->surf_bas;
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);

    Unused(dt);
    Thet1 = 1 - Theta;
    max_level = mg_params->max_level;
    min_level = mg_params->min_level;
    start_level = mg->start_level;
    u_L = p_L = NULL;  

    phgPrintf("   MultiGrid Build Mat on each grid level\n");
    phgMultiGridUpdate(mg); /* not need ? */
    elapsed_time(g, FALSE, 0.);



    /* --------------------------------------------------
     *
     *         MG for solver F
     *
     * --------------------------------------------------
     * */
    if (use_mg_F && i_slv == MG_SOLVER_F) {
	static BOOLEAN initialized = FALSE;
	phgPrintf("   Building MG solver F\n");

	set_active_solver(mg, MG_SOLVER_TYPE_U, 
			  MG_SOLVER_F);
	for (level = max_level-1; level >= 0; level--) {
	    phgPrintf("      Building mg level %d ", level);
	    if (level < start_level) {
		elapsed_time( g, TRUE, 0.);
		continue;
	    }

	    u_L = ml[level]->dofs[0];
	    assert(!strcmp(u_L->name, "dof u"));

	    /* Destroy old mat & solver */
	    phgMatDestroy(&ml[level]->mat_[MG_SOLVER_F]);
	    phgSolverDestroy(&ml[level]->solver_[MG_SOLVER_F]);

	    if (TRUE
		&& mg_params->reuse_mat 
		&& level == max_level - 1 
		) {
		//assert(ml[level]->mat->refcount == 0);

		assert(ns->matF != NULL);
		ml[level]->mat_[MG_SOLVER_F] = ns->matF;
		ns->matF->refcount++;

		if (use_GSline) {
		    ml[level]->mat = ns->matF;
		    ml[level]->block_dofs = ml[level]->block_dofs_[MG_SOLVER_F];
		    phgMultiGridFreeFactLineBlock(ml[level]);
		    phgMultiGridFactorizeLineBlock(ml[level]);
		    ml[level]->block_dofs = NULL;
		    ml[level]->mat = NULL;
		    phgMatUnpack(ns->matF);
		}
		elapsed_time( g, TRUE, 0.);
		continue;
	    }

	    ml[level]->mat = phgMapCreateMat(ml[level]->map, ml[level]->map);
	    ml[level]->mat->handle_bdry_eqns = TRUE;

	    /* dim > 1: build mat & rhs */
	    int nsub = 0;
	    if (DOF_TYPE_SUBELEM(u_L->type)) {
		nsub = u_L->type->order;
		phgPrintf("      Using Quad-Sub: %d, ", nsub);
	    }

	    ForAllElements(ml[level]->sub_grid, e) {
		int N = u_L->type->nbas;
		int q, order = 6;
		FLOAT A[N][Dim][N][Dim], buffer[N], rhsu[N][Dim];
		INT I[N][Dim], J[Dim][N];
		QUAD *quad;
		FLOAT vol, det;
		const FLOAT *w, *p, *gu, *vw;

		Bzero(A); Bzero(buffer); Bzero(rhsu);
		Bzero(I); Bzero(J); 

		for (i = 0; i < N; i++)
		    for (k = 0; k < Dim; k++)
			J[k][i] = I[i][k] = 
			    phgMapE2L(ml[level]->mat->cmap, 0, e, i*Dim + k);

		if (nsub > 0) {
		    quad = phgQuadGetQuad3DSub(2, nsub);
		} else {
		    quad = phgQuadGetQuad3D(order);
		}
		    
		//MG_DEBUGn(4, "   vol: %f\n", vol);
		gu = phgQuadGetDofValues(e, ns->gradu[1], quad);      /* grad u^{n+1}
								       *  */
		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    phgGeomGetCurvedJacobianAtLambda(ml[level]->sub_grid, e, p, &det);
		    vol = fabs(det / 6.);
		    for (i = 0; i < N; i++) {
			const FLOAT *gi_u = phgQuadGetBasisValues(e, u_L, i, quad) + q;               /* phi_i */
			const FLOAT *ggi_u = phgQuadGetBasisCurvedGradient(e, u_L, i, quad, q);    /* grad phi_i */
			for (j = 0; j < N; j++) {
			    const FLOAT *gj_u = phgQuadGetBasisValues(e, u_L, j, quad) + q;           /* phi_i */
			    const FLOAT *ggj_u = phgQuadGetBasisCurvedGradient(e, u_L, j, quad, q);    /* grad phi_i */
			    FLOAT mass = (*gj_u) * (*gi_u);
			    FLOAT diffu = INNER_PRODUCT(ggj_u, ggi_u);

			    Unused(diffu);
			    Unused(mass);
			    nu = get_effective_viscosity(gu, initialized);
#  if TENSOR_STRAIN
			    const FLOAT *tp = get_gbas_product(ggi_u, ggj_u);
			    for (k = 0; k < Dim; k++) 
				for (l = 0; l < Dim; l++) 
				    A[j][l][i][k] += vol*(*w) * EQU_SCALING * nu * tp[k+l*Dim];
#  else
			    for (k = 0; k < Dim; k++) 
				A[i][k][j][k] += vol*(*w) * EQU_SCALING * (.5*nu * diffu);
#  endif
			}
		    }
		    vw += Dim;
		    gu += DDim;
		    w++; p += Dim + 1;
		}
	    

#if TEST_CASE == ICE_BENCH_C			\
    || TEST_CASE == ICE_BENCH_D			\
    || TEST_CASE == ICE_BENCH_E			\
    || TEST_CASE == ICE_BENCH_F
		/* Rotate bases */
		for (i = 0; i < N; i++) {
		    INT id = phgDofMapE2D(surf_dof, e, i * (Dim*Dim)) / (Dim*Dim);
		    if (!rotated[id])
			continue;	
		    const FLOAT *trans = Trans + id*(Dim*Dim);
		
		    //SHOW_M(trans, Dim, Dim);
		    trans_left(&A[i][0][0][0], Dim*N, Dim*N, trans);
		    trans_rightT(&A[0][0][i][0], Dim*N, Dim*N, trans);
		}
#else
		Unused(Trans);
		Unused(rotated);
#endif


		for (i = 0; i < N; i++) {
		    for (k = 0; k < Dim; k++) {
			if (phgDofDirichletBC_(u_L, e, i*Dim+k, NULL, buffer, &rhsu[i][0],
					       DOF_PROJ_NONE)) {
			    phgMatAddEntries(ml[level]->mat, 1, I[i] + k, N, J[k], buffer);
			} else {
			    phgMatAddEntries(ml[level]->mat, 1, I[i] + k, N*Dim, I[0],
					     &(A[i][k][0][0]));
			}
		    }
		}
	    } /* end of elements */
	    phgMatAssemble(ml[level]->mat);

	    /* coarsest grid solve */
	    if (mg_params->solve_coarst
		&& level == mg_params->min_level) {
		int verb = phgVerbosity;
		phgOptionsPush();
		phgOptionsSetOptions("-solver superlu ");
		phgOptionsSetOptions(mg_params->coarst_opts);
		phgVerbosity = -1;
		ml[level]->solver = phgMat2Solver(SOLVER_DEFAULT, ml[level]->mat);
		phgVerbosity = verb;
		ml[level]->solver->warn_maxit = FALSE;
		phgOptionsPop();
	    }

	    if (use_GSline //&& level > mg_params->min_level) {
		&& ml[level]->solver->oem_solver != SOLVER_SUPERLU
		&& ml[level]->solver->oem_solver != SOLVER_MUMPS) {
		ml[level]->block_dofs = ml[level]->block_dofs_[MG_SOLVER_F];
		phgMultiGridFreeFactLineBlock(ml[level]);
		phgMultiGridFactorizeLineBlock(ml[level]); 
		ml[level]->block_dofs = NULL;
	    }
	    if (DUMP_MAT_VEC) {
		char mat_name[100], mat_file[100];
		sprintf(mat_name, "F%d", level);
		sprintf(mat_file, "F%d_.m", level);
		phgPrintf("*** Dumping %s\n", mat_name);
		phgMatDumpMATLAB(ml[level]->mat, mat_name, mat_file);
	    }

	    ml[level]->mat_[MG_SOLVER_F] = ml[level]->mat;
	    ml[level]->solver_[MG_SOLVER_F] = ml[level]->solver; 
	    ml[level]->mat = NULL;
	    ml[level]->solver = NULL;
	    elapsed_time( g, TRUE, 0.);
	}	  /* end of level */
	set_active_solver(mg, -1, -1);
	initialized = TRUE;
    }

    /* --------------------------------------------------
     *
     *         MG for solver Ap, Qp
     *
     * --------------------------------------------------
     * */

    if ((use_mg_Ap || use_mg_Qp)
	&& (i_slv == MG_SOLVER_Ap || 
	    i_slv == MG_SOLVER_Qp)) {
	static BOOLEAN initialized = FALSE;

	phgPrintf("   Building MG solver [AQZ]p\n");
	set_active_solver(mg, MG_SOLVER_TYPE_P, 
			  -1);

	for (level = max_level-1; level >= 0; level--) {
	    phgPrintf("      Building mg level %d ", level);
	    if (level < start_level) {
		elapsed_time( g, TRUE, 0.);
		continue;
	    }

	    u_L = ml[level]->dofs[0];
	    assert(!strcmp(u_L->name, "dof u"));
	    p_L = ml[level]->dofs[1];
	    assert(!strcmp(p_L->name, "dof p"));

	    /* Destroy old mat & solver */
	    phgMatDestroy(&ml[level]->mat_[MG_SOLVER_Ap]);
	    phgMatDestroy(&ml[level]->mat_[MG_SOLVER_Qp]);
	    phgSolverDestroy(&ml[level]->solver_[MG_SOLVER_Ap]);
	    phgSolverDestroy(&ml[level]->solver_[MG_SOLVER_Qp]);

	    ml[level]->mat = NULL;	/* not used */
	    ml[level]->solver = NULL;

	    if (TRUE
		&& mg_params->reuse_mat 
		&& level == max_level - 1 
		) {
		//assert(ml[level]->mat->refcount == 0);

		assert(_pcd->matAp != NULL);
		assert(_pcd->matQp != NULL);
		ml[level]->mat_[MG_SOLVER_Ap] = _pcd->matAp;
		ml[level]->mat_[MG_SOLVER_Qp] = _pcd->matQp;
		_pcd->matAp->refcount++;
		_pcd->matQp->refcount++;

		if (use_GSline) {
		    /* Ap */
		    ml[level]->mat = _pcd->matAp;
		    ml[level]->block_dofs = ml[level]->block_dofs_[MG_SOLVER_Ap];
		    phgMultiGridFreeFactLineBlock(ml[level]);
		    phgMultiGridFactorizeLineBlock(ml[level]);
		    phgMatUnpack(_pcd->matAp);

		    /* Qp */
		    ml[level]->mat = _pcd->matQp;
		    ml[level]->block_dofs = ml[level]->block_dofs_[MG_SOLVER_Qp];
		    phgMultiGridFactorizeLineBlock(ml[level]);
		    phgMatUnpack(_pcd->matQp);

		    ml[level]->block_dofs = NULL;
		    ml[level]->mat = NULL;
		}
		elapsed_time( g, TRUE, 0.);
		continue;
	    }

	    ml[level]->mat_[MG_SOLVER_Ap] = phgMapCreateMat(ml[level]->map, ml[level]->map);
	    ml[level]->mat_[MG_SOLVER_Qp] = phgMapCreateMat(ml[level]->map, ml[level]->map);
	    ml[level]->mat_[MG_SOLVER_Ap]->handle_bdry_eqns = TRUE;
	    ml[level]->mat_[MG_SOLVER_Qp]->handle_bdry_eqns = FALSE;

	    ForAllElements(ml[level]->sub_grid, e) {
		int N = p_L->type->nbas;
		int q, order = DofTypeOrder(p_L, e) * 3;
		FLOAT Ap[N][N], Qp[N][N], Zp[N][N], rhs[N], buffer[N];
		INT I[N];
		QUAD *quad;
		FLOAT vol, det;
		const FLOAT *w, *p;

		Bzero(Ap); Bzero(Qp); Bzero(Zp); 
		Bzero(rhs); Bzero(buffer);
	    
		for (i = 0; i < N; i++)
		    I[i] = phgMapE2L(ml[level]->mat_[MG_SOLVER_Ap]->cmap, 0, e, i);

		quad = phgQuadGetQuad3D(order);
		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    phgGeomGetCurvedJacobianAtLambda(ml[level]->sub_grid, e, p, &det);
		    vol = fabs(det / 6.);
		    for (i = 0; i < N; i++) {
			const FLOAT *gi = phgQuadGetBasisValues(e, p_L, i, quad) + q;       /* phi_i */
			const FLOAT *ggi = phgQuadGetBasisCurvedGradient(e, p_L, i, quad, q);    /* grad phi_i */
			for (j = 0; j < N; j++) {
			    const FLOAT *gj = phgQuadGetBasisValues(e, p_L, j, quad) + q;       /* phi_i */
			    const FLOAT *ggj = phgQuadGetBasisCurvedGradient(e, p_L, j, quad, q);    /* grad phi_i */

			    Ap[j][i] += (*w)*vol * INNER_PRODUCT(ggi, ggj);
			    Qp[j][i] += (*w)*vol * (*gj) * (*gi);
			}
		    }
		    w++; p += Dim + 1;
		}

		for (i = 0; i < N; i++) {
		    if (phgDofDirichletBC(p_L, e, i, NULL, buffer, NULL, DOF_PROJ_NONE)
			&& !(phgDofGetElementBoundaryType(p_L, e, i) & INFLOW) ) {
			phgMatAddEntries(ml[level]->mat_[MG_SOLVER_Ap], 1, I+i, N, I, buffer); 
		    }
		    else {	/* interior node */
			phgMatAddEntries(ml[level]->mat_[MG_SOLVER_Ap], 1, I+i, N, I, Ap[i]); 
		    }
		    phgMatAddEntries(ml[level]->mat_[MG_SOLVER_Qp], 1, I+i, N, I, Qp[i]); 
		}
	    } /* end of elements */

	    phgMatAssemble(ml[level]->mat_[MG_SOLVER_Ap]);
	    phgMatAssemble(ml[level]->mat_[MG_SOLVER_Qp]);

	    /* coarsest grid solve */
	    if (mg_params->solve_coarst
		&& level == mg_params->min_level) {
		phgOptionsPush();
		phgOptionsSetOptions("-solver superlu ");
		phgOptionsSetOptions(mg_params->coarst_opts);
		ml[level]->solver_[MG_SOLVER_Ap] = phgMat2Solver(SOLVER_DEFAULT, ml[level]->mat_[MG_SOLVER_Ap]);
		ml[level]->solver_[MG_SOLVER_Qp] = phgMat2Solver(SOLVER_DEFAULT, ml[level]->mat_[MG_SOLVER_Qp]);
		ml[level]->solver_[MG_SOLVER_Ap]->warn_maxit = FALSE;
		ml[level]->solver_[MG_SOLVER_Qp]->warn_maxit = FALSE;
		phgOptionsPop();
	    }

	    if (use_GSline //&& level > mg_params->min_level) {
		&& ml[level]->solver[MG_SOLVER_Qp].oem_solver != SOLVER_SUPERLU
		&& ml[level]->solver[MG_SOLVER_Qp].oem_solver != SOLVER_MUMPS) {

		/* Ap */
		ml[level]->mat = ml[level]->mat_[MG_SOLVER_Ap];
		ml[level]->block_dofs = ml[level]->block_dofs_[MG_SOLVER_Ap];
		phgMultiGridFreeFactLineBlock(ml[level]);
		phgMultiGridFactorizeLineBlock(ml[level]);
		ml[level]->block_dofs = NULL;
		ml[level]->mat = NULL;

		/* Qp */
		ml[level]->mat = ml[level]->mat_[MG_SOLVER_Ap];
		ml[level]->block_dofs = ml[level]->block_dofs_[MG_SOLVER_Qp];
		phgMultiGridFreeFactLineBlock(ml[level]);
		phgMultiGridFactorizeLineBlock(ml[level]);
		ml[level]->block_dofs = NULL;
		ml[level]->mat = NULL;
	    }

	    if (DUMP_MAT_VEC) {
		char mat_name[100], mat_file[100];

		sprintf(mat_name, "mg_Ap%d", level);
		sprintf(mat_file, "mg_Ap%d_.m", level);
		phgPrintf("*** Dumping %s\n", mat_name);
		phgMatDumpMATLAB(ml[level]->mat_[MG_SOLVER_Ap], mat_name, mat_file);

		sprintf(mat_name, "mg_Qp%d", level);
		sprintf(mat_file, "mg_Qp%d_.m", level);
		phgPrintf("*** Dumping %s\n", mat_name);
		phgMatDumpMATLAB(ml[level]->mat_[MG_SOLVER_Qp], mat_name, mat_file);
	    }

	    ml[level]->mat = NULL;
	    ml[level]->solver = NULL;
	    elapsed_time( g, TRUE, 0.);
	}	  /* end of level */
	set_active_solver(mg, -1, -1);
	initialized = TRUE;
    }

    return;
}

static double time0 = 0.;

void 
mg_pc_proc(SOLVER *pc, VEC *b0, VEC **x0)
{
    MULTI_GRID *mg = (MULTI_GRID *) pc->mat->mv_data[0];
    MG_LEVEL **ml = mg->ml;
    GRID *g = mg->grid;
    MAT *A = pc->mat;
    VEC *x;
    int level;

    /* MG cycle */
    level = mg_params->max_level - 1; /* work on finest level */
    phgInfo(0, "\n--- --- MG solve --- ---\n");

    MPI_Barrier(g->comm);
    time0 = phgGetTime(NULL);             
    x = phgMapCreateVec(A->rmap, 1);

    /* MG cycle setup */
    phgVecCopy(b0, &ml[level]->f);

    /* MG correct on x */
    recursive_MG_cycle(mg, level, TRUE);
    
    /* return smoother solution */
    phgVecCopy(ml[level]->x, x0);

    phgVecDestroy(&x);

    MPI_Barrier(g->comm);
    if (g->rank == 0)                               
	phgInfo(0, "   MG PC cost: %0.8lfs\n",  
		phgGetTime(NULL) - time0);
    return;
}


