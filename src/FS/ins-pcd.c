#include "ins.h"
#include "mat_op3.h"
#define _nsp (ns->ns_params)
#define _pcd (ns->pcd)

#define SUB_SOLVER_VERB -1

/***********************************/
/* Init NS solver's Preconditioner */
/***********************************/
void phgNSInitPc(NSSolver *ns)
{
    GRID *g = ns->g;
    MAP *Pmap = ns->Pmap, *Pbcmap;
    BOOLEAN use_Fu = _nsp->use_Fu;
    int verb;

    /* pcd boundary type test */
    _pcd->dof_inflow = phgDofNew(g, _nsp->ptype, 1, "dof inflow", DofNoAction);
    _pcd->dof_outflow = phgDofNew(g, _nsp->ptype, 1, "dof outflow", DofNoAction);
    _pcd->dof_nobdry = phgDofNew(g, _nsp->ptype, 1, "dof nobdry", DofNoAction);

    phgDofSetDirichletBoundaryMask(_pcd->dof_inflow, INFLOW);
    phgDofSetDirichletBoundaryMask(_pcd->dof_outflow, OUTFLOW);
    phgDofSetDirichletBoundaryMask(_pcd->dof_nobdry, 0);

    _pcd->map_inflow = phgMapCreate(_pcd->dof_inflow, NULL);
    _pcd->map_outflow = phgMapCreate(_pcd->dof_outflow, NULL);
    _pcd->map_nobdry = phgMapCreate(_pcd->dof_nobdry, NULL);

    _pcd->Pbcmap = phgMapCreate(_pcd->pbc, NULL);
    Pbcmap = _pcd->Pbcmap;

    Unused(Pmap);
#warning PCD B.C.: step 1. build mat using map... 
    /*
     * PCD boundary setup: should be consistent with code above
     */
    if (_nsp->pin_node) {
	_pcd->matFp = phgMapCreateMat(Pbcmap, Pbcmap);
	_pcd->matAp = phgMapCreateMat(Pbcmap, Pbcmap);
	_pcd->matQp = phgMapCreateMat(Pbcmap, Pbcmap);
    } else {
	//_pcd->matAp = phgMapCreateMat(_pcd->map_inflow, _pcd->map_inflow);
	//_pcd->matFp = phgMapCreateMat(_pcd->map_inflow, _pcd->map_inflow);
	_pcd->matFp = phgMapCreateMat(_pcd->map_outflow, _pcd->map_outflow);
	_pcd->matAp = phgMapCreateMat(_pcd->map_outflow, _pcd->map_outflow);
	//_pcd->matQp = phgMapCreateMat(_pcd->map_outflow, _pcd->map_outflow);
	//_pcd->matQp = phgMapCreateMat(_pcd->map_inflow, _pcd->map_inflow);
	//_pcd->matFp = phgMapCreateMat(_pcd->map_nobdry, _pcd->map_nobdry);
	//_pcd->matAp = phgMapCreateMat(_pcd->map_nobdry, _pcd->map_nobdry);
	//_pcd->matFp = phgMapCreateMat(_pcd->map_nobdry, _pcd->map_nobdry);
	_pcd->matQp = phgMapCreateMat(_pcd->map_nobdry, _pcd->map_nobdry);
    }

    /* stokes problem: get SYMETRIC mat when assemble.
     * Handle_bdry_eqns means mat is composed with row of boundary row
     * and non-bdry row, and eliminating mat columes of dirichlet dof. 
     */
    if (_nsp->use_symetric) {
	_pcd->matFp->handle_bdry_eqns = TRUE;
	_pcd->matAp->handle_bdry_eqns = TRUE;
	_pcd->matQp->handle_bdry_eqns = TRUE;
    } 
    /* genearl NS: no need to eliminate mat columes of dirichlet dof */
    else {
	_pcd->matFp->handle_bdry_eqns = FALSE;
	_pcd->matAp->handle_bdry_eqns = FALSE;
	_pcd->matQp->handle_bdry_eqns = FALSE;
    }	
    
    _pcd->rhsScale = phgMapCreateVec(_pcd->matQp->rmap, 1);
    phgVecDisassemble(_pcd->rhsScale);

    ns->pc = phgMat2Solver(SOLVER_PreOnly, ns->solver_u->mat);
    if (_nsp->use_PCD)
	phgSolverSetPC(ns->solver_u, ns->pc, pc_proc);

    /* solver F */
    phgOptionsPush();
    phgOptionsSetOptions("-solver hypre "
			 "-hypre_solver pcg "
			 "-hypre_pc boomeramg "
			 "-solver_maxit 10 "
			 "-solver_rtol 1e-4");
    phgOptionsSetOptions(_nsp->F_opts);
    /* use matF in the preconditioning matrix */
    _pcd->solver_F = phgMat2Solver(SOLVER_DEFAULT, ns->matF);
    _pcd->solver_F->verb = SUB_SOLVER_VERB; /* Set user options. */
    _pcd->pc_F = NULL;
#if USE_MG 
    if (ns_params->use_mg_F) {
	MAT *matF = ns->matF;

	assert(ns_params->use_PCD && !ns_params->use_Fu);
	_pcd->pc_F = phgMat2Solver(SOLVER_PreOnly, matF);
	phgOptionsSetOptions("-solver petsc ");
	matF->mv_data = phgAlloc(sizeof(*matF->mv_data));
	matF->mv_data[0] = (void *) ns->mg;
	phgSolverSetPC(_pcd->solver_F, _pcd->pc_F, mg_pc_proc);
    }
#endif /* USE_MG */
    _pcd->solver_F->warn_maxit = FALSE;
    phgOptionsPop();

    /* solver Ap */
    phgOptionsPush();
    phgOptionsSetOptions("-solver hypre "
			 "-hypre_solver gmres "
			 "-hypre_pc boomeramg "
			 "-solver_maxit 10 "
			 "-solver_rtol 1e-3");
    phgOptionsSetOptions(_nsp->Ap_opts);
    _pcd->solver_Ap = phgMat2Solver(SOLVER_DEFAULT, _pcd->matAp);
    _pcd->solver_Ap->warn_maxit = FALSE;
    _pcd->solver_Ap->verb = SUB_SOLVER_VERB;
    phgOptionsPop();
    _pcd->pc_Ap = NULL;
#if USE_MG 
    if (ns_params->use_mg_Ap) {
	MAT *matAp = _pcd->matAp;

	assert(ns_params->use_PCD);
	_pcd->pc_Ap = phgMat2Solver(SOLVER_PreOnly, matAp);
	phgOptionsSetOptions("-solver petsc ");
	matAp->mv_data = phgAlloc(sizeof(*matAp->mv_data));
	matAp->mv_data[0] = (void *) ns->mg;
	phgSolverSetPC(_pcd->solver_Ap, _pcd->pc_Ap, mg_pc_proc);
    }
#endif /* USE_MG */


    /* solver Qp */
    phgOptionsPush();
    phgOptionsSetOptions("-solver hypre "
			 "-hypre_solver pcg "
			 "-hypre_pc boomeramg "
			 "-solver_maxit 10 "
			 "-solver_rtol 1e-3");
    phgOptionsSetOptions(_nsp->Qp_opts);
    _pcd->solver_Qp = phgMat2Solver(SOLVER_DEFAULT, _pcd->matQp);
    _pcd->solver_Qp->warn_maxit = FALSE;
    _pcd->solver_Qp->verb = SUB_SOLVER_VERB;
    phgOptionsPop();
    _pcd->pc_Qp = NULL;
#if USE_MG 
    if (ns_params->use_mg_Qp) {
	MAT *matQp = _pcd->matQp;

	assert(ns_params->use_PCD);
	_pcd->pc_Qp = phgMat2Solver(SOLVER_PreOnly, matQp);
	phgOptionsSetOptions("-solver petsc ");
	matQp->mv_data = phgAlloc(sizeof(*matQp->mv_data));
	matQp->mv_data[0] = (void *) ns->mg;
	phgSolverSetPC(_pcd->solver_Qp, _pcd->pc_Qp, mg_pc_proc);
    }
#endif /* USE_MG */

    /* Fu for solve F^{-1} */
    if (use_Fu) { /* _nsp->implicit_centrip &&  */
	DOF *u1;
	MAP *u1map;
	MAT *matFu;
	u1 = _pcd->u1 = 
	    phgDofNew(g, _nsp->utype, 1, "velocity component u", DofNoAction);
	phgDofSetDirichletBoundaryMask(u1, SETFLOW);
	
	u1map = _pcd->u1map
	    = phgMapCreate(_pcd->u1, NULL);
	matFu = _pcd->matFu 
	    = phgMapCreateMat(u1map, u1map);
	if (_nsp->use_symetric)
	    matFu->handle_bdry_eqns = TRUE;

	/* solver Fu */
	phgOptionsPush();
	phgOptionsSetOptions("-solver hypre "
			     "-hypre_solver pcg "
			     "-hypre_pc boomeramg "
			     "-solver_maxit 10 "
			     "-solver_rtol 1e-4");
	phgOptionsSetOptions(_nsp->Fu_opts);
	_pcd->solver_Fu = phgMat2Solver(SOLVER_DEFAULT, _pcd->matFu);
	_pcd->solver_Fu->warn_maxit = FALSE;
	_pcd->solver_Fu->verb = SUB_SOLVER_VERB;
	phgOptionsPop();

    }

    return;
}


#define eu_xx eu[0]
#define eu_xy eu[1]
#define eu_xz eu[2]
#define eu_yx eu[3]
#define eu_yy eu[4]
#define eu_yz eu[5]
#define eu_zx eu[6]
#define eu_zy eu[7]
#define eu_zz eu[8]

static FLOAT * 
get_gbas_product(const FLOAT *gi, const FLOAT *gj,
		 const FLOAT *gu, LTYPE ltype) 
{
    static FLOAT prod[Dim][Dim];
    FLOAT Gi[Dim], Gj[Dim];
    FLOAT eu[DDim];
    FLOAT eps, eps2;
    int k;

    /* Picard term */
    prod[0][0] = gi[0] * gj[0] + .5 * (gi[1] * gj[1] + gi[2] * gj[2]);
    prod[0][1] = .5 * gi[0] * gj[1];
    prod[0][2] = .5 * gi[0] * gj[2];

    prod[1][0] = .5 * gi[1] * gj[0]; 
    prod[1][1] = gi[1] * gj[1] + .5 * (gi[0] * gj[0] + gi[2] * gj[2]);
    prod[1][2] = .5 * gi[1] * gj[2];

    prod[2][0] = .5 * gi[2] * gj[0]; 
    prod[2][1] = .5 * gi[2] * gj[1]; 
    prod[2][2] = gi[2] * gj[2] + .5 * (gi[0] * gj[0] + gi[1] * gj[1]);

    if (ltype == PICARD) {
	return prod[0];
    }

    /* Newton term */
    MAT3_SYM(gu, eu);
    for (k = 0; k < DDim; k++)
	eu[k] /= LEN_SCALING;
    eps = sqrt(.5) * MAT3_NORM2(eu);
    
    if (eps < MIN_EFFECTIVE_STRAIN) 
	eps = MIN_EFFECTIVE_STRAIN;

    eps2 = - (1./3.) / (eps*eps);

    Gi[0] = eu_xx * gi[0] + eu_xy * gi[1] + eu_xz * gi[2];
    Gi[1] = eu_yx * gi[0] + eu_yy * gi[1] + eu_yz * gi[2];
    Gi[2] = eu_zx * gi[0] + eu_zy * gi[1] + eu_zz * gi[2];

    Gj[0] = eu_xx * gj[0] + eu_xy * gj[1] + eu_xz * gj[2];
    Gj[1] = eu_yx * gj[0] + eu_yy * gj[1] + eu_yz * gj[2];
    Gj[2] = eu_zx * gj[0] + eu_zy * gj[1] + eu_zz * gj[2];
    
    prod[0][0] += Gi[0] * Gj[0] * eps2;
    prod[0][1] += Gi[1] * Gj[0] * eps2;
    prod[0][2] += Gi[2] * Gj[0] * eps2;

    prod[1][0] += Gi[0] * Gj[1] * eps2;
    prod[1][1] += Gi[1] * Gj[1] * eps2;
    prod[1][2] += Gi[2] * Gj[1] * eps2;

    prod[2][0] += Gi[0] * Gj[2] * eps2;
    prod[2][1] += Gi[1] * Gj[2] * eps2;
    prod[2][2] += Gi[2] * Gj[2] * eps2;

    return prod[0];
}

#undef eu_xx
#undef eu_xy
#undef eu_xz
#undef eu_yx
#undef eu_yy
#undef eu_yz
#undef eu_zx
#undef eu_zy
#undef eu_zz


#define DGETRF_FOR dgetrf_
#define DGETRI_FOR dgetri_
void DGETRF_FOR(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
void DGETRI_FOR(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);


#define SHOW_M(matN, mat_m, mat_n) {				\
	printf("\n### rank: %d\n", g->rank);			\
	int i, j;						\
	printf(" --- "#matN":(%3d * %3d)\n", mat_m, mat_n);	\
	for (i = 0; i < mat_m; i++) {				\
	    for (j = 0; j < mat_n; j++){			\
		printf("%14.8e, ", *(matN + i * (mat_n) + j));	\
	    }							\
	    printf("\n");					\
	}							\
    }


/***********************************/
/* Build Preconditioner's Matrices */
/***********************************/
void
phgNSBuildPc(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    FLOAT *dt = ns->dt;
    int i, j, q, s, k, l;
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1, nu0 = 0;
    DOF *tmp_u1 = phgDofNew(g, _nsp->utype, Dim, "tmp u1", func_u);
    int viscosity_type = ns->viscosity_type;
    LTYPE ltype = ns->ltype;


#if STEADY_STATE
    assert(fabs(Theta - 1) < 1e-12);
    Thet1 = 0; Unused(Thet1);
    Unused(dt);
#else
    Thet1 = 1 - Theta;
#endif /* STEADY_STATE */


    ForAllElements(g, e) {
	int M = ns->u[1]->type->nbas;	/* num of bases of Velocity */
	int N = ns->p[1]->type->nbas;	/* num of bases of Pressure */
	int order = 2 * DofTypeOrder(ns->p[1], e) + 
	    DofTypeOrder(ns->u[1], e) - 1; 	/* highest order term (u \nabla p, psi)  */
	FLOAT Ap[N][N], Fp[N][N], Qp[N][N], bufp[N], rhs1 = 1;
	FLOAT F[M*Dim][M*Dim], B[N][M*Dim], Bt[M*Dim][N];
	INT Ip[N];
	QUAD *quad;
	FLOAT vol, det;
	const FLOAT *w, *p, *vw, *gu, *vTe;

	quad = phgQuadGetQuad3D(order);
	vw = phgQuadGetDofValues(e, ns->wind, quad);  /* value wind */
	gu = phgQuadGetDofValues(e, ns->gradu[1], quad);        /* grad u^{n+1} */
	if (ns_params->noniter_temp)
	    vTe = phgQuadGetDofValues(e, ns->T[1], quad);  /* value temp */
	else
	    vTe = phgQuadGetDofValues(e, ns->T[0], quad);  /* value temp */
	
	vol = 0;
	Bzero(Ap); Bzero(Fp); Bzero(Qp); 
	Bzero(F); Bzero(Bt); Bzero(B);
	Bzero(bufp); 

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);

	    for (i = 0; i < N; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, ns->p[1], i, quad) + q;       /* phi_i */
		const FLOAT *ggi = phgQuadGetBasisCurvedGradient(e, ns->p[1], i, quad, q);    /* grad phi_i */
		for (j = 0; j < N; j++) {
		    const FLOAT *gj = phgQuadGetBasisValues(e, ns->p[1], j, quad) + q;       /* phi_j */
		    const FLOAT *ggj = phgQuadGetBasisCurvedGradient(e, ns->p[1], j, quad, q);    /* grad phi_i */
		    
		    nu = get_effective_viscosity(gu, *vTe, 0, viscosity_type);
		    if (i == 0 && j == 0)
			nu0 += nu;
#if ICE_BENCH_TEST ||				\
    ESIMINT_TEST ||				\
    HEINO_TEST ||				\
    TEST_CASE == ICE_EXACT	||		\
    TEST_CASE == ICE_GREEN_LAND
		    Unused(dt);
		    /* Note: B Q^-1 Bt ~ Ap(nu=1),
		     *       Fp(nu varies) is very different to Ap */
		    Ap[i][j] += vol*(*w) * INNER_PRODUCT(ggj, ggi);
#  if USE_QP_ONLY

		    //Qp[i][j] += vol*(*w) * LEN_SCALING * PRES_SCALING /(nu) * (*gj) * (*gi);
		    Qp[i][j] += vol*(*w) * 1. /(EQU_SCALING * nu) * (*gj) * (*gi);
		    /* if (i < NVert && j < NVert) { */
		    /* 	Qp[i][j] += vol*(*w) * LEN_SCALING * PRES_SCALING / (nu) * (*gj) * (*gi); */
		    /* } else if (i == NVert && j == NVert) { */
		    /* 	Qp[i][j] += vol*(*w) * LEN_SCALING * PRES_SCALING / (nu) * (*gj) * (*gi); */
		    /* } */

#  else
		    Qp[i][j] += vol*(*w) * (*gj) * (*gi);
#  endif
		    Fp[i][j] += vol*(*w) * (EQU_SCALING * nu * INNER_PRODUCT(ggj, ggi)
					    );
#elif STEADY_STATE 
		    Ap[i][j] += vol*(*w) * INNER_PRODUCT(ggj, ggi);
		    Qp[i][j] += vol*(*w) * (*gj) * (*gi);
		    Fp[i][j] += vol*(*w) * (nu * INNER_PRODUCT(ggj, ggi) * EQU_SCALING
					    );
#elif TIME_DEP_NON
		    Ap[i][j] += vol*(*w) * INNER_PRODUCT(ggj, ggi);
		    Qp[i][j] += vol*(*w) * (*gj) * (*gi);
		    Fp[i][j] += vol*(*w) * ((*gj) * (*gi) / dt[0]
					    + Theta * (nu * INNER_PRODUCT(ggj, ggi)
						       )
					    );
#else
		    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE */
		}
	    }

	    vw += Dim; 
	    gu += DDim;
	    vTe++;
	    w++; p += Dim+1;
	}


	/* Map: Element -> system */
	for (i = 0; i < N; i++) 
	    Ip[i] = phgMapE2L(_pcd->matFp->cmap, 0, e, i);

	/*
	 * PCD boundary setup I:
	 * Automaticly decide inflow boundary condition using wind direction.
	 *
	 * NOT active.
	 * */
	if (FALSE && !_nsp->pin_node) {
	    for (i = 0; i < N; i++) {
		BOOLEAN flag_inflow = FALSE;
		for (s = 0; s < NFace; s++) {
		    SHORT bases[NbasFace(ns->p[1])];
		    FLOAT *coord, vw_[3]; 
		    const FLOAT *lam, *normal;

		    if (!(e->bound_type[s] & BDRY_MASK))
			//if (!(e->bound_type[s] & INFLOW))
			continue;	/* boundary face */

		    phgDofGetBasesOnFace(ns->p[1], e, s, bases);
		    for (j = 0; j < NbasFace(ns->p[1]); j++) 
			if (i == bases[j]) {
			    normal = phgGeomGetFaceOutNormal(g, e, s);
			    coord = phgDofGetElementCoordinates(ns->p[1], e, i);
			    lam = phgGeomXYZ2Lambda(g, e, coord[0], coord[1], coord[2]);
			    phgDofEval(tmp_u1, e, lam, vw_);
			    if (INNER_PRODUCT(vw_, normal) > 1e-8) 
				flag_inflow = TRUE;
			}
		}
		
		if (flag_inflow) {
		    Bzero(bufp); bufp[i] = 1.0;
		    phgMatAddEntries(_pcd->matAp, 1, Ip + i, N, Ip, bufp);
		    phgMatAddEntries(_pcd->matFp, 1, Ip + i, N, Ip, bufp);
		    //phgMatAddEntries(_pcd->matQp, 1, Ip + i, N, Ip, bufp);
		    phgVecAddEntries(_pcd->rhsScale, 0, 1, Ip + i, &rhs1);
		}
		else {
		    /* interior node Or Neumann */
		    phgMatAddEntries(_pcd->matAp, 1, Ip + i, N, Ip, Ap[i]);
		    phgMatAddEntries(_pcd->matFp, 1, Ip + i, N, Ip, Fp[i]);
		    //phgMatAddEntries(_pcd->matQp, 1, Ip + i, N, Ip, Qp[i]);
		}
		phgMatAddEntries(_pcd->matQp, 1, Ip + i, N, Ip, Qp[i]);
	    }
	} 
	/*
	 * PCD boundary setup II:
	 * Enclose flow: use pinnode boundary.
	 *
	 * Qp is pinned, this is different to open flow.
	 * 
	 * */
	else if (_nsp->pin_node) {

	    for (i = 0; i < N; i++) {
		if (phgDofDirichletBC(_pcd->pbc, e, i, NULL, bufp, NULL, DOF_PROJ_NONE)) {
#if PIN_AT_ROOT 
		    if (g->rank != 0)
		    	phgError(1, "Pinned node only on rank 0!\n");
		    if (e->verts[i] != ns->pinned_node_id)
			phgError(1, "pinned node [%d] & [%d] doesn't coincide when build pc!\n",
				 e->verts[i], ns->pinned_node_id);
#else
		    if (GlobalVertex(g, e->verts[i]) != ns->pinned_node_id)
			phgError(1, "pinned node [%d] & [%d] doesn't coincide when build pc!\n",
				 e->verts[i], ns->pinned_node_id);
#endif /* PIN_AT_ROOT */

		    phgMatAddEntries(_pcd->matAp, 1, Ip + i, N, Ip, bufp);
		    phgMatAddEntries(_pcd->matFp, 1, Ip + i, N, Ip, bufp);
		    phgMatAddEntries(_pcd->matQp, 1, Ip + i, N, Ip, bufp);
		    phgVecAddEntries(_pcd->rhsScale, 0, 1, Ip + i, &rhs1);
		} else {
		    /* interior node Or Neumann */
		    phgMatAddEntries(_pcd->matAp, 1, Ip + i, N, Ip, Ap[i]);
		    phgMatAddEntries(_pcd->matFp, 1, Ip + i, N, Ip, Fp[i]);
		    phgMatAddEntries(_pcd->matQp, 1, Ip + i, N, Ip, Qp[i]);
		}
	    }
	}
	/*
	 * PCD boundary setup III:
	 * Open flow: there could be varies kinds of combination on seting up
	 *   boundary conditon, but Inflow:Robin & Outflow:scaled Dirich is
	 *   prefered. See Ref[2].
	 * 
	 * */
	else {
	    for (i = 0; i < N; i++) {

		/*****************/
                /* Inflow	 */
                /*****************/
#warning PCD B.C.: Step 2.1. build mat, all neumann, add dirich entries
		if (FALSE && phgDofDirichletBC(_pcd->dof_inflow, e, i, NULL, bufp, NULL, DOF_PROJ_NONE)) {
		    phgMatAddEntries(_pcd->matAp, 1, Ip + i, N, Ip, bufp);
		    phgMatAddEntries(_pcd->matFp, 1, Ip + i, N, Ip, bufp);
		    phgMatAddEntries(_pcd->matQp, 1, Ip + i, N, Ip, bufp);
		    phgVecAddEntries(_pcd->rhsScale, 0, 1, Ip + i, &rhs1);
		} else if (FALSE && phgDofDirichletBC(_pcd->dof_outflow, e, i, NULL, bufp, NULL, DOF_PROJ_NONE)
			   && !(phgDofGetElementBoundaryType(ns->p[1], e, i) & INFLOW) ) {

		    ERROR_MSG("Fp, Qp");
		    nu = get_effective_viscosity(NULL, 0, 0, viscosity_type);
		    phgMatAddEntries(_pcd->matAp, 1, Ip + i, N, Ip, bufp);
		    bufp[i] *= EQU_SCALING * nu;
		    phgMatAddEntries(_pcd->matFp, 1, Ip + i, N, Ip, bufp);
		    phgVecAddEntries(_pcd->rhsScale, 0, 1, Ip + i, &rhs1);

		    //phgMatAddEntries(_pcd->matQp, 1, Ip + i, N, Ip, bufp);
		} else if (FALSE && phgDofDirichletBC(_pcd->pbc, e, i, NULL, bufp, NULL, DOF_PROJ_NONE)) {
		    phgMatAddEntries(_pcd->matAp, 1, Ip + i, N, Ip, bufp);
		    phgMatAddEntries(_pcd->matFp, 1, Ip + i, N, Ip, bufp);
		    phgMatAddEntries(_pcd->matQp, 1, Ip + i, N, Ip, bufp);
		    phgVecAddEntries(_pcd->rhsScale, 0, 1, Ip + i, &rhs1);
		}
		else if (FALSE) {
		    /* interior node Or Neumann */

		    ERROR_MSG("Fp, Qp");
		    phgMatAddEntries(_pcd->matAp, 1, Ip + i, N, Ip, Ap[i]);
		    phgMatAddEntries(_pcd->matFp, 1, Ip + i, N, Ip, Fp[i]);
		    //phgMatAddEntries(_pcd->matQp, 1, Ip + i, N, Ip, Qp[i]);
		}

		/******************/
                /* No bdry	  */
                /******************/
		//phgMatAddEntries(_pcd->matFp, 1, Ip + i, N, Ip, Fp[i]);
		phgMatAddEntries(_pcd->matQp, 1, Ip + i, N, Ip, Qp[i]);
	    }
	}


	if (0) {
	    /* Special term <[[p_i]], [[p_j]]> */
	    int face;
	    nu0 /= quad->npoints;
	    for (face = 0; face < NFace; face++) {
		FLOAT area =  phgGeomGetFaceArea(g, e, face);
		//FLOAT value = {area, -area};
		FLOAT values[2] = {vol * 1. /(EQU_SCALING * nu0),
				   -vol * 1. /(EQU_SCALING * nu0)};
		SIMPLEX *e_neigh;

		phgMatAddEntries(_pcd->matQp, 1, Ip+NVert, 1, Ip+NVert, values);
		if ((e_neigh = GetNeighbour(e, face)) != NULL) {
		    INT Ip_neigh = phgMapE2L(_pcd->matFp->cmap, 0, e_neigh, NVert);
		    phgMatAddEntries(_pcd->matQp, 1, Ip+NVert, 1, &Ip_neigh, values + 1);
		}
	    }
	}

    }	/* end element */

    phgMatAssemble(_pcd->matFp);
#if 0
    /* 
     * Once there is a bug in building Fp, which takes dim of wind to be Dim+1,
     *   this causes severe and hidden memory leakage problem that is hard to find.
     *
     * This bug is fixed by ZLB, and the following code serves as a safeguard against
     *   similar problems.
     * 
     *  */
    for (i = 0; i < _pcd->matFp->rmap->nlocal; i++) {
	const MAT_ROW *row = phgMatGetRow(_pcd->matFp, i);
	for (j = 0; j < row->ncols; j++) {
	    assert(row->cols[j] >= 0 && row->cols[j] < _pcd->matFp->cmap->nglobal);
	    assert(Fabs(row->data[j]) < 1e+18);
	}
    }
#endif	/* check matFp assemble */
    phgMatAssemble(_pcd->matAp);
    phgMatAssemble(_pcd->matQp);
    phgVecAssemble(_pcd->rhsScale);

    //phgPrintf("rhs scale norm: %E\n", phgVecNorm2(_pcd->rhsScale, 0, NULL));

    /* scale Dirich */
    if (FALSE) {
	MAP *map = _pcd->rhsScale->map;
	FLOAT *vs = _pcd->rhsScale->data;
	MAT_ROW *rowAp = _pcd->matAp->rows;
	MAT_ROW *rowFp = _pcd->matFp->rows;
	MAT_ROW *rowQp = _pcd->matQp->rows;
        INT i;

	ERROR_MSG("Fp, Qp");
	for (i = 0; i < map->nlocal; i++, vs++,
		 rowAp++, rowQp++, rowFp++) {
	    if (*vs > .5) {
		assert(rowAp->ncols == 1);
		assert(rowFp->ncols == 1);
		rowAp->data[0] /= *vs;
		rowFp->data[0] /= *vs;
	    }
	}
    }
    
    if (DUMP_MAT_VEC) {
	phgPrintf("Dumping [AFQ]p\n");
	phgMatDumpMATLAB(_pcd->matAp, "Ap", "Ap_.m");
	phgMatDumpMATLAB(_pcd->matFp, "Fp", "Fp_.m");
	phgMatDumpMATLAB(_pcd->matQp, "Qp", "Qp_.m");
    }

    phgDofFree(&tmp_u1);
    return;
}

/*
 *  Preconditioning procedure:
 *     p = - Qp^-1 Fp Ap^-1 * q 
 *     u = F^-1 ( bu - Bt * p )
 *
 *  PCD Dirichlet boundary scaling:
 *
 *  Note: instead of scaling diag entries of [AQF]p, we scale rhs entries,
 *    this is only because of LIMITATION of implementation.
 *  Todo: use scaling diag entries.
 *  
 *  Ref[2]: Boundary conditions in approximate commutator preconditioners for
 *   the Navier-Stokes equations, Howard C. Elman and Ray S. Tuminaro
 *
 *  */
void
pc_proc(SOLVER *pc_solver, VEC *b0, VEC **x0)
{
    NSSolver *ns = (NSSolver *) pc_solver->mat->mv_data[0];
#if USE_MG 
    MULTI_GRID *mg = ns->mg;
#endif /* USE_MG */
    GRID *g = ns->g;
    VEC *xu, *xp, *xu2, *xp2;
    INT verb = 1;
    int i, k;
    int Nu = ns->matF->rmap->nlocal;
    int Np = _pcd->matAp->rmap->nlocal;
    SOLVER *solver_F = _pcd->solver_F, *solver_Ap = _pcd->solver_Ap, 
	*solver_Qp = _pcd->solver_Qp;
    FLOAT *rhsF, *rhsAp, *rhsQp;
    BOOLEAN use_Fu = _nsp->use_Fu;
    double t;


    /* Note:
     * 1. To get a symetric matrix,
     *    set handle_bdry_eqns = FALSE, then MatAssemble
     * 2. To make the matrix, symetric or NOT,  directly multiply rhs
     *    set rhs_updated = TRUE, then SolverVecSolve
     * TODO: why ???
     */
    solver_F->rhs_updated = TRUE;
    solver_Ap->rhs_updated = TRUE;
    solver_Qp->rhs_updated = TRUE;

    /* save old rhs */
    rhsF = solver_F->rhs->data;
    rhsAp = solver_Ap->rhs->data;
    rhsQp = solver_Qp->rhs->data;

    xu = phgMapCreateVec(ns->matF->rmap, 1);
    xu2 = phgMapCreateVec(ns->matF->rmap, 1);
    xp = phgMapCreateVec(solver_Ap->mat->rmap, 1);
    xp2 = phgMapCreateVec(solver_Ap->mat->rmap, 1);
    memcpy(xp->data, b0->data + Nu, sizeof(*xp->data) * Np);

    
    /***********************/
    /* PCD preconditioning */
    /***********************/

    /* xp = - Ap^-1 Fp Qp^-1 * q */
    t = phgGetTime(NULL);
    solver_Qp->rhs->data = xp->data;
    solver_Qp->rhs->assembled = TRUE;
    bzero(xp2->data, sizeof(*xp2->data) * Np);
#if USE_MG 
    if (ns_params->use_mg_Qp)
        set_active_solver(mg, MG_SOLVER_TYPE_P, MG_SOLVER_Qp);
#endif /* USE_MG */
    phgSolverVecSolve(solver_Qp, FALSE, xp2); 
    memcpy(xp->data, xp2->data, sizeof(*xp->data) * Np);
    if (verb > 0)
	phgPrintf("\t    Qp: nits = %3d, residual = %0.4le [%0.4lgMB %0.4lfs]\n",
		  solver_Qp->nits, (double)solver_Qp->residual,
		  phgMemoryUsage(g, NULL) / (1024.0 * 1024.0),
		  phgGetTime(NULL) - t);

#  if !USE_QP_ONLY
    #include "Fp.c"
#  endif
    phgVecAXPBY(0., NULL, -1.0, &xp);


    /* xu = F^-1 ( bu - Bt * xp ) */
    t = phgGetTime(NULL);
    memcpy(xu->data, b0->data, sizeof(*xu->data) * Nu);
    phgMatVec(MAT_OP_N, -1.0, ns->matBt, xp, 1., &xu); 


    if (!use_Fu) {
	/* use F in the preconditioning matrix */
	solver_F->rhs->data = xu->data;
	solver_F->rhs->assembled = TRUE;
	bzero(xu2->data, sizeof(*xu2->data) * Nu);
#if USE_MG 
        if (ns_params->use_mg_F)
            set_active_solver(mg, MG_SOLVER_TYPE_U, MG_SOLVER_F);
#endif /* USE_MG */
	phgSolverVecSolve(solver_F, FALSE, xu2); 
	phgVecAXPBY(1.0, xu2, 0.0, &xu);
	if (verb > 0)
	    phgPrintf("\t    F : nits = %3d, residual = %0.4le [%0.4lgMB %0.4lfs]\n",
		      solver_F->nits, (double)solver_F->residual,
		      phgMemoryUsage(g, NULL) / (1024.0 * 1024.0),
		      phgGetTime(NULL) - t);
    } else {
	SOLVER *solver_Fu = _pcd->solver_Fu;
	VEC *xu0;
	FLOAT *rhsFu; 
	INT nits[3], Nu0 = _pcd->matFu->rmap->nlocal;
	float res[3];

	assert(Nu0 * Dim == Nu); 
	
	rhsFu = solver_Fu->rhs->data;
	xu0 = phgMapCreateVec(solver_Fu->mat->rmap, 1);

	t = phgGetTime(NULL);
	Bzero(nits);
	Bzero(res);
	for (k = 0; k < Dim; k++) {
	    /* solve Fu * xu0 = xu_k */
	    for (i = 0; i < Nu0; i++)
		rhsFu[i] = xu->data[3 * i + k];
	    solver_Fu->rhs->assembled = TRUE;
	    solver_Fu->rhs_updated = TRUE;
	    bzero(xu0->data, sizeof(*xu->data) * Nu0);
	    phgSolverVecSolve(solver_Fu, FALSE, xu0);
	    nits[k] = solver_Fu->nits;
	    res[k] = (double)solver_Fu->residual;
	    for (i = 0; i < Nu0; i++)
		xu->data[3 * i + k] = xu0->data[i];
	}
	if (verb > 0)
	    phgPrintf("\t    Fu*3: nits = (%d,%d,%d), residual = (%0.4le,%0.4le,%0.4le)"
		      " [%0.4lgMB %0.4lfs]\n",
		      nits[0], nits[1], nits[2], 
		      res[0], res[1], res[2], 
		      phgMemoryUsage(g, NULL) / (1024.0 * 1024.0),
		      phgGetTime(NULL) - t);

	solver_Fu->rhs->data = rhsFu;
	phgVecDestroy(&xu0);
    } 


    /* Copy xu, xp to x */
    memcpy((*x0)->data, xu->data, sizeof(*xu->data) * Nu);
    memcpy((*x0)->data + Nu, xp->data, sizeof(*xp->data) * Np);

    /* restore rhs */
    solver_F->rhs->data = rhsF;
    solver_Ap->rhs->data = rhsAp;
    solver_Qp->rhs->data = rhsQp;

    phgVecDestroy(&xu);
    phgVecDestroy(&xu2);
    phgVecDestroy(&xp);
    phgVecDestroy(&xp2);

    return;
}


/********************************************/
/* Destroy NS solver's Preconditioner	    */
/* Note: call before phgNSDestroySolver(),  */
/*   to free solver->mat->mv_data.	    */
/********************************************/
void phgNSDestroyPc(NSSolver *ns)
{
    BOOLEAN use_Fu = _nsp->use_Fu;
    //phgFree(ns->matNS->mv_data);

    phgMapDestroy(&_pcd->map_inflow);
    phgMapDestroy(&_pcd->map_outflow);
    phgMapDestroy(&_pcd->map_nobdry);

    phgDofFree(&_pcd->dof_nobdry);
    phgDofFree(&_pcd->dof_inflow);
    phgDofFree(&_pcd->dof_outflow);

    phgMatDestroy(&_pcd->matAp);
    phgMatDestroy(&_pcd->matQp);
    phgMatDestroy(&_pcd->matFp);

    phgMapDestroy(&_pcd->Pbcmap);
    phgSolverDestroy(&ns->pc);

    phgSolverDestroy(&_pcd->solver_F);
    phgSolverDestroy(&_pcd->solver_Ap);
    phgSolverDestroy(&_pcd->solver_Qp);

    phgSolverDestroy(&_pcd->pc_F);
    phgSolverDestroy(&_pcd->pc_Ap);
    phgSolverDestroy(&_pcd->pc_Qp);

    phgVecDestroy(&_pcd->rhsScale);

    if (use_Fu) {
	phgMatDestroy(&_pcd->matFu);
	phgSolverDestroy(&_pcd->solver_Fu);
	phgMapDestroy(&_pcd->u1map);
	phgDofFree(&_pcd->u1);
    }
}
