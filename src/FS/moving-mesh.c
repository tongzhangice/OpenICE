#include "ins.h"
#if USE_PETSC
#  include <petscksp.h>
#endif

#define _nsp (ns->ns_params)


BTYPE phgDofGetBoundaryType(DOF *u, INT index);

static void get_moved_coord_FEM(NSSolver *ns, int nstep);
static void get_moved_coord_FV_cell(NSSolver *ns, int nstep);
static void get_moved_coord_FV_point(NSSolver *ns, int nstep);

static double accum = 0.5; // m a-1


/*
 * Scheme:
 * 1. FEM implicit
 * 2. FEM explicit
 * 3. FV cell based
 * 4. FV point based
 *
 *  */
void 
get_moved_coord(NSSolver *ns, int tstep)
{
    if (_nsp->height_scheme == 1) {
	phgPrintf("Height solver scheme: FEM implicit.\n");
	get_moved_coord_FEM(ns, tstep);
    } else if (_nsp->height_scheme == 0) {
	phgPrintf("Height solver scheme: FEM explicit.\n");
	get_moved_coord_FEM(ns, tstep);
    } else if (_nsp->height_scheme == 2) {
	phgPrintf("Height solver scheme: FV cell.\n");
	get_moved_coord_FV_cell(ns, tstep);
    } else if (_nsp->height_scheme == 3) {
	phgPrintf("Height solver scheme: FV point.\n");
	get_moved_coord_FV_point(ns, tstep);
    } else {
	phgError(0, "Wrong scheme!\n");
    }
}

void inv3(FLOAT a[3][3], FLOAT ia[3][3])
{
    int i, j;
    FLOAT det = a[0][0]*a[1][1]*a[2][2] - a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2] 
	+ a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0];
    FLOAT b[3][3] = {
	{ a[1][1]*a[2][2] - a[1][2]*a[2][1], a[0][2]*a[2][1] - a[0][1]*a[2][2], a[0][1]*a[1][2] - a[0][2]*a[1][1]},
	{ a[1][2]*a[2][0] - a[1][0]*a[2][2], a[0][0]*a[2][2] - a[0][2]*a[2][0], a[0][2]*a[1][0] - a[0][0]*a[1][2]},
	{ a[1][0]*a[2][1] - a[1][1]*a[2][0], a[0][1]*a[2][0] - a[0][0]*a[2][1], a[0][0]*a[1][1] - a[0][1]*a[1][0]}
    };
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
	ia[i][j] = b[i][j] / det;
}

static void 
get_grad_bas(const FLOAT *J, const FLOAT *G,
		  FLOAT *p)
{
    *(p++) = G[0] * J[0 * (Dim + 1) + 0] +
	G[1] * J[1 * (Dim + 1) + 0] +
	G[2] * J[2 * (Dim + 1) + 0] +
	G[3] * J[3 * (Dim + 1) + 0];
    *(p++) = G[0] * J[0 * (Dim + 1) + 1] +
	G[1] * J[1 * (Dim + 1) + 1] +
	G[2] * J[2 * (Dim + 1) + 1] +
	G[3] * J[3 * (Dim + 1) + 1];
    *(p++) = G[0] * J[0 * (Dim + 1) + 2] +
	G[1] * J[1 * (Dim + 1) + 2] +
	G[2] * J[2 * (Dim + 1) + 2] +
	G[3] * J[3 * (Dim + 1) + 2];
}

static void
func_z(FLOAT x, FLOAT y, FLOAT z, FLOAT *v) 
{
    *v = z;
}



/* ================================================================================
 *
 *   update height using FEM 
 *   
 * ================================================================================
 */


/*
 *  
 *  Get Surface height change
 *  Output dH |_{bdry surf}
 *
 *  */
#if 0
void
get_surf_dH(NSSolver *ns)
{
    GRID *g = ns->g;
    DOF *u = ns->u[1];
    DOF *dH = ns->dH;
    FLOAT *dt = ns->dt;
    int height_scheme = _nsp->height_scheme;
    SOLVER *solver;
    SIMPLEX *e;
    INT i, j, s;
    int verb;
    
	DOF *avg_n = phgDofNew(g, DOF_P2, 3, "avg n", DofNoAction);
        get_avg_n(g, avg_n);
        
    //phgDofSetDirichletBoundaryMask(dH, BC_BOT);
    //phgOptionsPush();
    //phgOptionsSetOptions("-solver gmres "
	//		 "-gmres_pc_type jacobi "
	//		 "-solver_maxit 10000 "
	//		 "-solver_rtol 1e-12");
    //phgOptionsSetOptions("-solver petsc -oem_options 'pc_type lu -pc_factor_mat_solver_package superlu'");
    //verb = phgVerbosity;
    //phgVerbosity = 1;
    solver = phgSolverCreate(SOLVER_DEFAULT, dH, NULL);
    //phgVerbosity = verb;
    //phgOptionsPop();


    ForAllElements(g, e) {
	int N = dH->type->nbas;	/* number of basises in the element */
	FLOAT A[N][N], rhs[N], buffer[N];
	INT I[N];
	BTYPE btype[N];
	int q, dof_i;
	FLOAT *avg_n_v;
	FLOAT normal1[Dim];

	QUAD *quad;
	const FLOAT *w, *p, *J;

	for (i = 0; i < N; i++) 
	    I[i] = phgSolverMapE2L(solver, 0, e, i);

	bzero(A, sizeof(A));
	bzero(rhs, sizeof(rhs));
	bzero(btype, sizeof(btype));
	J = phgGeomGetJacobian(g, e);

    /*
	for (i = 0; i < N; i++) {
	    int gind, bind;
	    GTYPE gtype;
	    btype[i] = phgDofGetElementBasisInfo(dH, e, i,
						 &gtype, &gind, &bind);
	}
    */
	//SHOW_iV_(0, btype, N);

	for (s = 0; s < NFace; s++) {
	    int v0, v1, v2, ii, jj;
	    int nbas_face = NbasFace(dH);
	    SHORT bases[nbas_face];
	    FLOAT x, y, z, lambda[Dim + 1], area, vu[Dim],
		gi, gj, ggi[Dim], ggj[Dim], qa;
	    const FLOAT *normal, *Gi, *Gj;

	    if ((e->bound_type[s] & BC_TOP) || (e->bound_type[s] & BC_BOTTOM))
        {
	    
	    phgDofGetBasesOnFace(dH, e, s, bases);
	    v0 = GetFaceVertex(s, 0);
	    v1 = GetFaceVertex(s, 1);
	    v2 = GetFaceVertex(s, 2);
	    lambda[s] = 0.;

	    area = phgGeomGetFaceArea(g, e, s);
	    const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);
	    quad = phgQuadGetQuad2D(5);

	    /* Projection to 2D plan */
	    //area *= normal[Z_DIR];

	    p = quad->points;
	    w = quad->weights;
        const FLOAT *u_quad = phgQuadGetDofValues(e, u, quad);
	    for (q = 0; q < quad->npoints; q++) {
		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);
		
		phgDofEval(u, e, lambda, vu);
        const FLOAT *u_quad_q = u_quad + q*Dim;
		phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
        const FLOAT *bas = dH->type->BasFuncs(dH, e, 0, -1, lambda);

        if (e->bound_type[s] & BC_TOP)
            func_smb_top(x, y, z, &qa);
            // smb: surface mass balance
        else if (e->bound_type[s] & BC_BOTTOM)
            func_smb_bot(x, y, z, &qa);
            // smb beneath the grounded region


		if (height_scheme == 1) {

		    /* ---------------------
		     * 
		     *  Surf move implicit
		     *  See notes about scaling.
		     *  Use dH as H^{n+1}		 *
		     * --------------------- */


		    FLOAT diam = phgGeomGetFaceDiameter(g, e, s);
		    FLOAT normU = sqrt(vu[0]*vu[0] + vu[1]*vu[1] + vu[2]*vu[2]);
		    FLOAT tau = 1 / sqrt( SQUARE(4*normU/diam)
					  + SQUARE(2/dt[0]));

            const FLOAT *bas_v = dH->type->BasFuncs(dH, e, 0, -1, lambda);
            const FLOAT *gbas_v = dH->type->BasGrads(dH, e, 0, -1, lambda);

		    for (ii = 0; ii < nbas_face; ii++) {
			i = bases[ii];
            gi = bas_v[i];
			//gi = *dH->type->BasFuncs(dH, e, i, i + 1, lambda);
			//Gi = dH->type->BasGrads(dH, e, i, i + 1, lambda);
            Gi = gbas_v + i*(Dim+1);
			get_grad_bas(J, Gi, ggi);

			FLOAT ti = (vu[0] * ggi[0] + vu[1] * ggi[1]);

			for (jj = 0; jj < nbas_face; jj++) { 
			    j = bases[jj];
			    //gj = *dH->type->BasFuncs(dH, e, j, j + 1, lambda);
			    //Gj = dH->type->BasGrads(dH, e, j, j + 1, lambda);
                gj = bas_v[j];
                Gj = gbas_v + j*(Dim+1);
			    get_grad_bas(J, Gj, ggj);

			    FLOAT mass = (gj)*(gi);
			    FLOAT conv = (vu[0] * ggj[0] + vu[1] * ggj[1]);
			    //FLOAT conv = 5 * (-0.1 * ggj[0] + 0.1 * ggj[1]);
			    FLOAT diff = 10 * (ggj[0]*ggi[0] + ggj[1]*ggi[1]);

#  if 0
			    /* Pure convection */
			    A[i][j] += area*(*w) * (gj*gi/dt[0] 
						    + conv * (gi)
						    ); 
#  elif 0
			    /* Artificial diff */
			    A[i][j] += area*(*w) * (gj*gi * LEN_SCALING / dt[0] 
						    + conv * (gi)
						    + 0.001 * diff
						    ); 
#  else
			    /* SUPG */
			    A[i][j] += area*(*w) * (gj * LEN_SCALING / dt[0] 
						    + conv 
						    ) * (gi + tau * ti);
#  endif
			} 
            //phgPrintf("qa: %e\n", qa+vu[2]);
			rhs[i] += area*(*w) * (z/dt[0] 
					       + vu[2] + qa
					       ) * (gi + tau * ti*0);

		    }
		    w++;
		} 
        else {
		    /* ---------------------
		     * 
		     *  Surf move explicit
		     *
		     * --------------------- */

		    FLOAT diam = phgGeomGetFaceDiameter(g, e, s);
		    FLOAT normU = sqrt(vu[0]*vu[0] + vu[1]*vu[1] + vu[2]*vu[2]);

		    for (ii = 0; ii < nbas_face; ii++) {
			i = bases[ii];
            dof_i = phgDofMapE2D(avg_n, e, i*Dim);
            avg_n_v = DofData(avg_n);
            normal1[0] = avg_n_v[dof_i + 0];
            normal1[1] = avg_n_v[dof_i + 1];
            normal1[2] = avg_n_v[dof_i + 2];
			//gi = *dH->type->BasFuncs(dH, e, i, i + 1, lambda);
			//Gi = dH->type->BasGrads(dH, e, i, i + 1, lambda);
			//get_grad_bas(J, Gi, ggi);
			for (jj = 0; jj < nbas_face; jj++) { 
			    j = bases[jj];
			    //gj = *dH->type->BasFuncs(dH, e, j, j + 1, lambda);
			    //Gj = dH->type->BasGrads(dH, e, j, j + 1, lambda);
			    //get_grad_bas(J, Gj, ggj);

			    //FLOAT mass = (gj)*(gi);
			    A[i][j] += area*(*w) * bas[i]*bas[j] + 1*area*(*w)*bas[j]*diam/(2*normU);
			}
			FLOAT dh = vu[Z_DIR]
			    - vu[X_DIR] * (-normal1[X_DIR]/normal1[Z_DIR])
			    - vu[Y_DIR] * (-normal1[Y_DIR]/normal1[Z_DIR])
			    + qa;
			//SHOW_V_(0, vu, 3);
			//SHOW_V_(0, normal, 3);
			//phgInfo(0, "  d[%d] += %e\n", i, dh);

			dh *= dt[0];


			rhs[i] += area*(*w) * (dh) * bas[i] / LEN_SCALING;
#  if 0
			printf("   U: (%12.8f %12.8f %12.8f)\n",
			       vu[0], vu[1], vu[2]);
			printf("  gH: (%12.8f %12.8f %12.8f)\n",
			       normal[0]/normal[2], normal[1]/normal[2], 1.);
#  endif
			assert(!isnan(rhs[i])); 
		    }
		    w++;
		} /* end scheme */
	    }	  /* end quad pts */
	}	  /* end face */
    }

	/* volume dof is set to zero */
	for (i = 0; i < N; i++) 
    {
        btype[i] = phgDofGetElementBoundaryType(dH, e, i);
	    if (!(btype[i] & BC_TOP) && !(btype[i] & BC_BOTTOM)) 
        {
            //bzero(A[i], N*sizeof(A[0][0]));
            A[i][i] = 1.;
            rhs[i] = 0.;
	    } 
    }

	//SHOW_V_(0, rhs, N);
	for (i = 0; i < N; i++) 
    {
	    if (phgDofDirichletBC_(dH, e, i, NULL, buffer, rhs+i, DOF_PROJ_NONE)) 
        {
            /* No bdry mask */
            //ERROR_MSG("no B.C. for dH");
            printf("no bc !\n");
            phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer); 
	    }
	    else 
        {
            phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
	    }
	}

        phgSolverAddRHSEntries(solver, N, I, rhs);

    }
    // end all elements


    if (DUMP_MAT_VEC) 
    {
        phgPrintf("Dumping dH\n");
        phgMatDumpMATLAB(solver->mat, "Ah", "Ah_.m");
        phgVecDumpMATLAB(solver->rhs, "bh", "bh_.m");
    }

    if (height_scheme == 0) 
    {
        phgDofSetDataByValue(dH, 0.);
    } 
    else 
    {
        phgDofSetDataByValue(dH, 0.);
        //phgDofSetBdryDataByFunction(dH, func_z, BC_TOP);
    }

    phgSolverSolve(solver, FALSE, dH, NULL);
        phgExportVTK(g, "dH.vtk", dH, NULL);
        // This dH is actually the H, not dH, for the implicit method.
        // However, for the explicit  method, dH is dH.
    phgPrintf("H_new min, max: %e %e\n", phgDofMinValVec(dH), phgDofMaxValVec(dH));
    
    phgPrintf("   DH solve: nits %d, res %e\n", 
	      solver->nits, solver->residual);
    phgSolverDestroy(&solver);

    if (height_scheme == 1) 
    {
        //DOF *H_old = phgDofNew(g, DOF_ANALYTIC, 1, "height old", func_z);
        DOF *H_old = phgDofNew(g, DOF_P1, 1, "height old", DofInterpolation);
        DOF *q = phgDofNew(g, DOF_P1, 1, "q", DofInterpolation);
        phgDofSetDataByFunction(H_old, func_z);
        phgDofSetDataByFunction(q, func_smb_top);
        phgPrintf("H_old min, max: %e %e\n", phgDofMinValVec(H_old), phgDofMaxValVec(H_old));
        phgExportVTK(g, "H.vtk", dH, H_old, q, NULL);
        phgDofAXPY(-1., H_old, &dH);
        phgPrintf("dH min, max: %e %e\n", phgDofMinValVec(dH), phgDofMaxValVec(dH));
        phgExportVTK(g, "dH.vtk", dH, H_old, NULL);
        phgDofFree(&H_old);
        phgDofFree(&q);
    }

    DOF_SCALE(dH, "dh");

    //phgDofDump(elev);
    //phgExportVTK(g, "./output/dHs.vtk", dH, NULL);
    return;
}
#endif

# if 1
void
get_surf_dH(NSSolver *ns)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    DOF *u = ns->u[1];
    DOF *dH = ns->dH;
    FLOAT *dt = ns->dt;
    int height_scheme = _nsp->height_scheme;
    SOLVER *solver;
    SIMPLEX *e;
    INT ii, i, j, k, s;
    int verb_core = 0;
    

    static FLOAT *surf_velo, *surf_velo0, *surf_elev, *surf_elev0;
    if (surf_velo == NULL) {
	PHG_CALLOC(surf_velo, 3*gL->nvert);
	PHG_CALLOC(surf_velo0, 3*gL->nvert);
	PHG_CALLOC(surf_elev, gL->nvert);
	PHG_CALLOC(surf_elev0, gL->nvert);
    }
    /* clear */
    for (i = 0; i < gL->nvert; i++) {
	surf_velo[3*i + 0] = -1e30;
	surf_velo[3*i + 1] = -1e30;
	surf_velo[3*i + 2] = -1e30;
	surf_elev[i] = -1e30;
    }

    /* Collect vertex velocity, elevation */
    for (ii = 0; ii < gL->nvert_bot; ii++) {
        i = gL->vert_bot_Lidx[ii];
        assert(gL->vert_local_lists[i] != NULL);

        INT iG = gL->vert_bot_Gidx[ii];
        assert(iG < gL->nvert);

        int nv = gL->vert_local_lists[i][0];
        int *iL = &gL->vert_local_lists[i][1];

        assert(nv > 0);
        //for (j = 0; j < nv; j++) {
#if 0
	j = nv-1;
#else
	j = 0;
#endif
	{
            FLOAT *vu = DofVertexData(ns->u[1], iL[j]);
	    surf_velo[3*iG  ] = vu[0];
	    surf_velo[3*iG+1] = vu[1];
	    surf_velo[3*iG+2] = vu[2];
	    surf_elev[iG] = g->verts[iL[j]][Z_DIR];
        }
    }

    /* Gather to root */
    MPI_Allreduce(surf_velo, surf_velo0, 3 * gL->nvert,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(surf_elev, surf_elev0, gL->nvert,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    if (1) {// || verb_core && phgRank == 0)
	FILE *fp = fopen("vel.dat", "w");
	for (ii = 0; ii < gL->nvert; ii++) {
	    fprintf(fp, "%e %e %e %e %e\n", 
		   gL->verts[ii][0],
		   surf_velo0[3*ii], 
		   surf_velo0[3*ii+1],
		   surf_velo0[3*ii+2], 
		   surf_elev0[ii]);
	}
	fclose(fp);
    }


    /* FEM 2D Solve */
    if (phgRank == 0) {
	MAP *map = phgMapCreateSimpleMap(MPI_COMM_SELF, gL->nvert, gL->nvert);
	MAT *mat = phgMapCreateMat(map, map);
	MAT *mat0 = phgMapCreateMat(map, map);
	SOLVER *solver = NULL;
	SOLVER *solver0 = NULL;

	VEC *rhs = phgMapCreateVec(map, 1);
	VEC *rhs0 = phgMapCreateVec(map, 1);
	VEC *vec = phgMapCreateVec(map, 1);
	VEC *vec0 = phgMapCreateVec(map, 1);
	phgVecDisassemble(rhs);
	phgVecDisassemble(rhs0);
	
	TRIA *t = gL->trias;
	for (ii = 0; ii < gL->ntria; ii++, t++) {
	    QUAD *quad;
	    int q;
	    const FLOAT *w, *p;
	    int I[3] = {t->verts[0],
			t->verts[1],
			t->verts[2]};
	    FLOAT *x[3] = {gL->verts[I[0]], 
			   gL->verts[I[1]], 
			   gL->verts[I[2]]};

	    FLOAT iJ[3][3] = {
		{x[0][0], x[1][0], x[2][0]},
		{x[0][1], x[1][1], x[2][1]},
		{1, 1, 1},
	    }, JJ[3][3];
	    FLOAT area = t->area;
	    FLOAT A[3][3], b[3];
	    FLOAT A0[3][3], b0[3];

	    bzero(A, sizeof(A));
	    bzero(b, sizeof(b));
	    bzero(A0, sizeof(A0));
	    bzero(b0, sizeof(b0));

	    inv3(iJ, JJ);


	    FLOAT diam = sqrt(SQUARE(x[0][0] - x[1][0])
	    		      + SQUARE(x[1][0] - x[2][0])
	    		      + SQUARE(x[2][0] - x[0][0]));

        //FLOAT diam = sqrt(SQUARE(x[0][0] - x[1][0]) + SQUARE(x[0][1] - x[1][1]))
        //            + sqrt(SQUARE(x[1][0] - x[2][0]) + SQUARE(x[1][1] - x[2][1]))
        //            + sqrt(SQUARE(x[2][0] - x[0][0]) + SQUARE(x[2][1] - x[0][1]));

        FLOAT h1 = fabs(x[0][0] - x[1][0]);
        FLOAT h2 = fabs(x[1][0] - x[2][0]);
        FLOAT h3 = fabs(x[2][0] - x[0][0]);

        FLOAT h = MAX(h1, h2); h = MAX(h, h3);


	    quad = phgQuadGetQuad2D(2);
	    p = quad->points;
	    w = quad->weights;

	    FLOAT P[3] = {1./3., 1./3, 1./3.};
	    FLOAT vu[3] = {0, 0, 0};
	    for (i = 0; i < 3; i++) {
	    	FLOAT gi = p[i];
	    	for (k = 0; k < Dim; k++)
	    	    vu[k] += surf_velo0[3*I[i] + k] * gi;
	    }

	    for (q = 0; q < quad->npoints; q++) {

		if (height_scheme == 1) {
		    /*
		     *
		     * Implicit
		     *
		     * */
		    FLOAT vu[3] = {0, 0, 0}, z = 0;
		    for (i = 0; i < 3; i++) {
			FLOAT gi = p[i];
			FLOAT *ggi = JJ[i];
			for (k = 0; k < Dim; k++)
			    vu[k] += surf_velo0[3*I[i] + k] * gi;
			z += surf_elev0[I[i]] * gi;
		    }


            FLOAT normU = sqrt(vu[0]*vu[0] + vu[1]*vu[1] + 0*vu[2]*vu[2])+1e-10;
            FLOAT tau = h/(2*normU);
            //FLOAT tau = 1 / sqrt( SQUARE(4*normU/diam)
            //              + SQUARE(2/dt[0]*LEN_SCALING));

            //printf("tau1 tau: %f %f\n", tau1, tau);

		    for (i = 0; i < 3; i++) {
			FLOAT gi = p[i];
			FLOAT *ggi = JJ[i];
			FLOAT ti = (vu[0] * ggi[0] + vu[1] * ggi[1]);

			for (j = 0; j < 3; j++) {
			    FLOAT gj = p[j];
			    FLOAT *ggj = JJ[j];

			    A[i][j] += area*(*w) * (gj*LEN_SCALING / dt[0]
						    + 1*(vu[0] * ggj[0] + vu[1] * ggj[1] + 0*gj*h/(2.*normU))) * (gi + 1*tau*ti);
			    A0[i][j] += area*(*w) * (gj*LEN_SCALING / dt[0]
						    + 0*(vu[0] * ggj[0] + vu[1] * ggj[1] + 0*gj*h/(2.*normU))) * (gi + 1*tau*ti);
			    
			}

			b[i] += area*(*w) * (1*z * LEN_SCALING / dt[0]+ 1*vu[2]) * (gi + 1*tau*ti);
			b0[i] += area*(*w) * (1*z * LEN_SCALING / dt[0]+ 0*vu[2]) * (gi + 1*tau*ti);
		    }
		} 
		else {
		    /*
		     *
		     * Explicit
		     *
		     * */
		    FLOAT vu[3] = {0, 0, 0}, dSx = 0, dSy = 0, x = 0;
		    for (i = 0; i < 3; i++) {
			FLOAT gi = p[i];
			FLOAT *ggi = JJ[i];

			for (k = 0; k < Dim; k++)
			    vu[k] += surf_velo0[3*I[i] + k] * gi;

			dSx += surf_elev0[I[i]] * ggi[0];
			dSy += surf_elev0[I[i]] * ggi[1];
			x += gL->verts[I[i]][0] * gi;
		    }

            FLOAT normU = sqrt(vu[0]*vu[0] + vu[1]*vu[1] + vu[2]*vu[2]);

		    for (i = 0; i < 3; i++) {
			FLOAT gi = p[i];
			FLOAT *ggi = JJ[i];
			for (j = 0; j < 3; j++) {
			    FLOAT gj = p[j];
			    FLOAT *ggj = JJ[j];

			    A[i][j] += area*(*w) * (gi*gj) + dt[0]*area*(*w)*gj*diam/(2*normU);
			}

			FLOAT dh = vu[Z_DIR] 
			    - vu[X_DIR] * (dSx)
			    - vu[Y_DIR] * (dSy);
			dh *= dt[0];
			 //printf("%e\t%e\t%e\t%e\n",  
			 //       x,  
			 //       vu[Z_DIR],  
			 //       vu[X_DIR],  
			 //       dSx); 
			b[i] += area*(*w) * (dh) * gi / LEN_SCALING;
		    }
		}
		
		p += 3;
		w++;
	    } /* end quad */

	    for (i = 0; i < 3; i++)
		phgMatAddEntries(mat, 1, I+i, 3, I, A[i]);
	    phgVecAddEntries(rhs, 0, 3, I, b);

	    for (i = 0; i < 3; i++)
		phgMatAddEntries(mat0, 1, I+i, 3, I, A0[i]);
	    phgVecAddEntries(rhs0, 0, 3, I, b0);
	}
	phgMatAssemble(mat);
	phgVecAssemble(rhs);

	phgMatAssemble(mat0);
	phgVecAssemble(rhs0);


	//phgOptionsPush();
	//phgOptionsSetOptions("-solver pcg -pcg_pc_type jacobi");
	//phgOptionsSetOptions("-solver petsc -oem_options -ksp_type bcgsl -oem_options -pc_type sor -oem_options -pc_sor_its 10");
    //phgOptionsSetOptions("-solver petsc -oem_options '-ksp_type preonly -pc_type lu -pc_type lu -pc_factor_mat_solver_package superlu'");
    //phgOptionsSetOptions("-solver petsc");
    //phgOptionsSetOptions("-oem_options '-ksp_type preonly'");
    //phgOptionsSetOptions("-oem_options '-pc_type lu'");
    //phgOptionsSetOptions("-oem_options '-pc_factor_mat_solver_package superlu_dist'");
    //phgOptionsSetOptions("-oem_options '-pc_factor_mat_solver_type superlu'");
	solver = phgMat2Solver(SOLVER_SUPERLU, mat);
	solver0 = phgMat2Solver(SOLVER_SUPERLU, mat0);
	//phgOptionsPop();
	phgMatDestroy(&mat);
	phgSolverAssemble(solver);
	phgMatDestroy(&mat0);
	phgSolverAssemble(solver0);


	/* solve */
	solver->rhs->assembled = TRUE;
	phgVecAXPBY(1., rhs, 0, &solver->rhs);
	solver0->rhs->assembled = TRUE;
	phgVecAXPBY(1., rhs0, 0, &solver0->rhs);
	//phgPrintf("* Solver 2D.\n");
	phgSolverVecSolve(solver, TRUE, vec);
	phgSolverVecSolve(solver0, TRUE, vec0);

    phgVecAXPBY(-1., vec0, 1., &vec);

	memcpy(surf_elev,
	       vec->data, gL->nvert * sizeof(*surf_elev));
	//memcpy(surf_elev0,
	//       vec->data, gL->nvert * sizeof(*surf_elev0));

	/* Dump to check */
	if (1 || verb_core) {
	    FILE *fp = fopen("ds.dat", "w");
	    for (i = 0; i < gL->nvert; i++) 
		fprintf(fp, "%e %e %e %e\n", 
			gL->verts[i][0], 
			gL->verts[i][1], 
			0.,
			vec->data[i]
			);
	    fclose(fp);
	}

	phgVecDestroy(&vec);
	phgMapDestroy(&map);
	//exit(0);
    }
    
    /* Broadcast to other */
    MPI_Bcast(surf_elev, gL->nvert,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);
    
     for (ii = 0; ii < gL->nvert_bot; ii++) { 
     	i = gL->vert_bot_Lidx[ii]; 
     	assert(gL->vert_local_lists[i] != NULL); 

     	INT iG = gL->vert_bot_Gidx[ii]; 
     	assert(iG < gL->nvert); 

     	int nv = gL->vert_local_lists[i][0]; 
     	int *iL = &gL->vert_local_lists[i][1]; 
	    
     	assert(nv > 0); 
        FLOAT *vg = DofVertexData(dH, iL[0]); 
        vg[0] = surf_elev[iG]; 
     } 
    /* //phgExportVTK(g, "grad_surf.vtk", dof_gradS, NULL); */

#if 1
    /* clear */
    for (i = 0; i < gL->nvert; i++) {
	surf_velo[3*i + 0] = -1e30;
	surf_velo[3*i + 1] = -1e30;
	surf_velo[3*i + 2] = -1e30;
	surf_elev[i] = -1e30;
    }
    for (ii = 0; ii < gL->nvert_bot; ii++) {
        i = gL->vert_bot_Lidx[ii];
        assert(gL->vert_local_lists[i] != NULL);

        INT iG = gL->vert_bot_Gidx[ii];
        assert(iG < gL->nvert);

        int nv = gL->vert_local_lists[i][0];
        int *iL = &gL->vert_local_lists[i][1];

        assert(nv > 0);
        //for (j = 0; j < nv; j++) {
#if 1
	j = nv-1;
#else
	j = 0;
#endif
	{
            FLOAT *vu = DofVertexData(ns->u[1], iL[j]);
	    surf_velo[3*iG  ] = vu[0];
	    surf_velo[3*iG+1] = vu[1];
	    surf_velo[3*iG+2] = vu[2];
	    surf_elev[iG] = g->verts[iL[j]][Z_DIR];
        }
    }

    /* Gather to root */
    MPI_Allreduce(surf_velo, surf_velo0, 3 * gL->nvert,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(surf_elev, surf_elev0, gL->nvert,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    if (1) {// || verb_core && phgRank == 0)
	FILE *fp = fopen("vel.dat", "w");
	for (ii = 0; ii < gL->nvert; ii++) {
	    fprintf(fp, "%e %e %e %e %e\n", 
		   gL->verts[ii][0],
		   surf_velo0[3*ii], 
		   surf_velo0[3*ii+1],
		   surf_velo0[3*ii+2], 
		   surf_elev0[ii]);
	}
	fclose(fp);
    }


    /* FEM 2D Solve */
    if (phgRank == 0) {
	MAP *map = phgMapCreateSimpleMap(MPI_COMM_SELF, gL->nvert, gL->nvert);
	MAT *mat = phgMapCreateMat(map, map);
	MAT *mat0 = phgMapCreateMat(map, map);
	SOLVER *solver = NULL;
	SOLVER *solver0 = NULL;

	VEC *rhs = phgMapCreateVec(map, 1);
	VEC *rhs0 = phgMapCreateVec(map, 1);
	VEC *vec = phgMapCreateVec(map, 1);
	VEC *vec0 = phgMapCreateVec(map, 1);
	phgVecDisassemble(rhs);
	phgVecDisassemble(rhs0);
	
	TRIA *t = gL->trias;
	for (ii = 0; ii < gL->ntria; ii++, t++) {
	    QUAD *quad;
	    int q;
	    const FLOAT *w, *p;
	    int I[3] = {t->verts[0],
			t->verts[1],
			t->verts[2]};
	    FLOAT *x[3] = {gL->verts[I[0]], 
			   gL->verts[I[1]], 
			   gL->verts[I[2]]};

	    FLOAT iJ[3][3] = {
		{x[0][0], x[1][0], x[2][0]},
		{x[0][1], x[1][1], x[2][1]},
		{1, 1, 1},
	    }, JJ[3][3];
	    FLOAT area = t->area;
	    FLOAT A[3][3], b[3];
	    FLOAT A0[3][3], b0[3];

	    bzero(A, sizeof(A));
	    bzero(b, sizeof(b));
	    bzero(A0, sizeof(A0));
	    bzero(b0, sizeof(b0));

	    inv3(iJ, JJ);


	    FLOAT diam = sqrt(SQUARE(x[0][0] - x[1][0])
	    		      + SQUARE(x[1][0] - x[2][0])
	    		      + SQUARE(x[2][0] - x[0][0]));

        //FLOAT diam = sqrt(SQUARE(x[0][0] - x[1][0]) + SQUARE(x[0][1] - x[1][1]))
        //            + sqrt(SQUARE(x[1][0] - x[2][0]) + SQUARE(x[1][1] - x[2][1]))
        //            + sqrt(SQUARE(x[2][0] - x[0][0]) + SQUARE(x[2][1] - x[0][1]));

        FLOAT h1 = fabs(x[0][0] - x[1][0]);
        FLOAT h2 = fabs(x[1][0] - x[2][0]);
        FLOAT h3 = fabs(x[2][0] - x[0][0]);

        FLOAT h = MAX(h1, h2); h = MAX(h, h3);


	    quad = phgQuadGetQuad2D(2);
	    p = quad->points;
	    w = quad->weights;

	    FLOAT P[3] = {1./3., 1./3, 1./3.};
	    FLOAT vu[3] = {0, 0, 0};
	    for (i = 0; i < 3; i++) {
	    	FLOAT gi = p[i];
	    	for (k = 0; k < Dim; k++)
	    	    vu[k] += surf_velo0[3*I[i] + k] * gi;
	    }

	    for (q = 0; q < quad->npoints; q++) {

		if (height_scheme == 1) {
		    /*
		     *
		     * Implicit
		     *
		     * */
		    FLOAT vu[3] = {0, 0, 0}, z = 0;
		    for (i = 0; i < 3; i++) {
			FLOAT gi = p[i];
			FLOAT *ggi = JJ[i];
			for (k = 0; k < Dim; k++)
			    vu[k] += surf_velo0[3*I[i] + k] * gi;
			z += surf_elev0[I[i]] * gi;
		    }


            FLOAT normU = sqrt(vu[0]*vu[0] + vu[1]*vu[1] + 0*vu[2]*vu[2])+1e-10;
            FLOAT tau = h/(2*normU);
            //FLOAT tau = 1 / sqrt( SQUARE(4*normU/diam)
            //              + SQUARE(2/dt[0]*LEN_SCALING));

            //printf("tau1 tau: %f %f\n", tau1, tau);

		    for (i = 0; i < 3; i++) {
			FLOAT gi = p[i];
			FLOAT *ggi = JJ[i];
			FLOAT ti = (vu[0] * ggi[0] + vu[1] * ggi[1]);

			for (j = 0; j < 3; j++) {
			    FLOAT gj = p[j];
			    FLOAT *ggj = JJ[j];

			    A[i][j] += area*(*w) * (gj*LEN_SCALING / dt[0]
						    + 1*(vu[0] * ggj[0] + vu[1] * ggj[1] + 0*gj*h/(2.*normU))) * (gi + 1*tau*ti);
			    A0[i][j] += area*(*w) * (gj*LEN_SCALING / dt[0]
						    + 0*(vu[0] * ggj[0] + vu[1] * ggj[1] + 0*gj*h/(2.*normU))) * (gi + 1*tau*ti);
			    
			}

			b[i] += area*(*w) * (1*z * LEN_SCALING / dt[0]+ vu[2] + accum*dt[0]) * (gi + 1*tau*ti);
			b0[i] += area*(*w) * (1*z * LEN_SCALING / dt[0]+ 0*vu[2]) * (gi + 1*tau*ti);
		    }
		} 
		else {
		    /*
		     *
		     * Explicit
		     *
		     * */
		    FLOAT vu[3] = {0, 0, 0}, dSx = 0, dSy = 0, x = 0;
		    for (i = 0; i < 3; i++) {
			FLOAT gi = p[i];
			FLOAT *ggi = JJ[i];

			for (k = 0; k < Dim; k++)
			    vu[k] += surf_velo0[3*I[i] + k] * gi;

			dSx += surf_elev0[I[i]] * ggi[0];
			dSy += surf_elev0[I[i]] * ggi[1];
			x += gL->verts[I[i]][0] * gi;
		    }

            FLOAT normU = sqrt(vu[0]*vu[0] + vu[1]*vu[1] + vu[2]*vu[2]);

		    for (i = 0; i < 3; i++) {
			FLOAT gi = p[i];
			FLOAT *ggi = JJ[i];
			for (j = 0; j < 3; j++) {
			    FLOAT gj = p[j];
			    FLOAT *ggj = JJ[j];

			    A[i][j] += area*(*w) * (gi*gj) + dt[0]*area*(*w)*gj*diam/(2*normU);
			}

			FLOAT dh = vu[Z_DIR] 
			    - vu[X_DIR] * (dSx)
			    - vu[Y_DIR] * (dSy);
			dh *= dt[0];
			 //printf("%e\t%e\t%e\t%e\n",  
			 //       x,  
			 //       vu[Z_DIR],  
			 //       vu[X_DIR],  
			 //       dSx); 
			b[i] += area*(*w) * (dh) * gi / LEN_SCALING;
		    }
		}
		
		p += 3;
		w++;
	    } /* end quad */

	    for (i = 0; i < 3; i++)
		phgMatAddEntries(mat, 1, I+i, 3, I, A[i]);
	    phgVecAddEntries(rhs, 0, 3, I, b);

	    for (i = 0; i < 3; i++)
		phgMatAddEntries(mat0, 1, I+i, 3, I, A0[i]);
	    phgVecAddEntries(rhs0, 0, 3, I, b0);
	}
	phgMatAssemble(mat);
	phgVecAssemble(rhs);

	phgMatAssemble(mat0);
	phgVecAssemble(rhs0);


	phgOptionsPush();
	//phgOptionsSetOptions("-solver petsc -oem_options -ksp_type bcgsl -oem_options -pc_type sor -oem_options -pc_sor_its 10");
    //phgOptionsSetOptions("-solver superlu");
			 //"-solver_maxit 10000 "
			 //"-solver_rtol 1e-10");
	phgOptionsSetOptions("-solver superlu");
	//phgOptionsSetOptions("-solver_maxit 10000");
	//phgOptionsSetOptions("-solver_rtol 1e-10");
    //phgOptionsSetOptions("-oem_options '-pc_type lu'");
    //phgOptionsSetOptions("-oem_options '-pc_factor_mat_solver_type umfpack'");
	solver = phgMat2Solver(SOLVER_DEFAULT, mat);
	solver0 = phgMat2Solver(SOLVER_DEFAULT, mat0);
	//phgOptionsSetOptions("-solver pcg -pcg_pc_type jacobi");
	//phgOptionsSetOptions("-solver petsc -oem_options -pc_type lu -oem_options -pc_factor_mat_solver_type umfpack");
	phgOptionsPop();
	phgMatDestroy(&mat);
	phgSolverAssemble(solver);
	phgMatDestroy(&mat0);
	phgSolverAssemble(solver0);


	/* solve */
	solver->rhs->assembled = TRUE;
	phgVecAXPBY(1., rhs, 0, &solver->rhs);
	solver0->rhs->assembled = TRUE;
	phgVecAXPBY(1., rhs0, 0, &solver0->rhs);
	//phgPrintf("* Solver 2D.\n");
	phgSolverVecSolve(solver, TRUE, vec);
	phgSolverVecSolve(solver0, TRUE, vec0);

    phgVecAXPBY(-1., vec0, 1., &vec);

	memcpy(surf_elev,
	       vec->data, gL->nvert * sizeof(*surf_elev));
	//memcpy(surf_elev0,
	//       vec->data, gL->nvert * sizeof(*surf_elev0));

	/* Dump to check */
	if (1 || verb_core) {
	    FILE *fp = fopen("ds.dat", "w");
	    for (i = 0; i < gL->nvert; i++) 
		fprintf(fp, "%e %e %e %e\n", 
			gL->verts[i][0], 
			gL->verts[i][1], 
			0.,
			vec->data[i]
			);
	    fclose(fp);
	}

	phgVecDestroy(&vec);
	phgMapDestroy(&map);
	//exit(0);
    }
    
    /* Broadcast to other */
    MPI_Bcast(surf_elev, gL->nvert,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);
    
     for (ii = 0; ii < gL->nvert_bot; ii++) { 
     	i = gL->vert_bot_Lidx[ii]; 
     	assert(gL->vert_local_lists[i] != NULL); 

     	INT iG = gL->vert_bot_Gidx[ii]; 
     	assert(iG < gL->nvert); 

     	int nv = gL->vert_local_lists[i][0]; 
     	int *iL = &gL->vert_local_lists[i][1]; 
	    
     	assert(nv > 0); 
        FLOAT *vg = DofVertexData(dH, iL[nv-1]); 
        vg[0] = surf_elev[iG]; 
     } 
#endif

     phgExportVTK(ns->g, "dH.vtk", dH, NULL);
    return;
}
#endif


/*
 * Get grid point z coord change
 * Output dH.
 * */
static void 
get_moved_coord_FEM(NSSolver *ns, int nstep)
{
    GRID *g = ns->g;
    DOF *u = ns->u[0];
    DOF *dH = ns->dH;
    LAYERED_MESH *gL = ns->gL;
    INT i, j;
    int ii, rank;
    FLOAT *sbuf, *rbuf, *dHS, *vH;
    FLOAT dH_max = -1e10, dH_min = 1e10;
    int verb;

    /* Use layer info */
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	assert(nv > 0);
	
	FLOAT dh, dh0, dh1;
	dh0 = *DofVertexData(dH, iL[0]);
	dh1 = *DofVertexData(dH, iL[nv-1]);
    // dh1 is the dH value on top of ii

    FLOAT bed_z;
    FLOAT x_coord = g->verts[iL[0]][0];
    FLOAT y_coord = g->verts[iL[0]][1];
    FLOAT z_coord = g->verts[iL[0]][2]+dh0;
    //if (x_coord < 1e6)
    //printf("dh0 dh1 %e %e %e %e\n", dh0, dh1, x_coord, y_coord);
    //phgPrintf("x y %e %e\n", x_coord, y_coord);
    func_bed_z(x_coord, y_coord, z_coord, &bed_z);
    if (z_coord < bed_z)
    {
        dh0 = bed_z - g->verts[iL[0]][2];
        //printf("dh0 x y %e %e %e %d\n", dh0, x_coord, y_coord, nv);
    }
        // if the node is at the grounded part, we set dh0 = 0
        // if the node is at the floating part, we set dh0 = -z_original + z_bed
        // Thus, we make sure the ice bottom doesn't penetrate the bedrock

    //if ((g->types_vert[iL[0]] & BC_BOTTOM_GRD) && !(g->types_vert[iL[0]] & BC_ISHELF))
    if ((g->types_vert[iL[0]] & BC_BOTTOM_GRD))
        //dh0 = 0;
        // make sure the ice bottom at the grounded part is stick to the bedrock

    //if (phgDofGetBoundaryType(dH, iL[nv-1]) & BC_LATERL)
    //    dh1 = 0;

	/* Set zero height to eps,
	 * dH is not used in this case. */
	if (TRUE) {
	    FLOAT *z = &g->verts[iL[nv-1]][Z_DIR];
	    if (*z + dh1  < HEIGHT_EPS) 
		dh1 = HEIGHT_EPS - *z; /* use H_eps */
	}


	dh = (dh1 - dh0) / (nv-1.0);

	dH_max = MAX(dH_max, dh1);
	dH_min = MIN(dH_min, dh1);

	for (j = 0; j < nv; j++) 
    {
	    //*DofVertexData(dH, iL[j]) = j * dh + dh0;
	    g->verts[iL[j]][Z_DIR] += j*dh + dh0;

        if (j == 0)
        if ((g->types_vert[iL[j]] & BC_BOTTOM_GRD))
            g->verts[iL[j]][Z_DIR] = bed_z;


	}

    if (g->verts[iL[0]][2] > bed_z)
    {
        //printf("bedz z dh0 %e %e %e\n", bed_z, g->verts[iL[0]][2], dh0);
    }


    if (0)
    {
        FILE *fp1 = fopen("bed.txt", "a");
        FILE *fp2 = fopen("sur.txt", "a");

        fprintf(fp1, "%f %f %f\n", x_coord, y_coord, g->verts[iL[0]][2]);
        fprintf(fp2, "%f %f %f\n", x_coord, y_coord, g->verts[iL[nv-1]][2]);
        fclose(fp1); fclose(fp2);
    }

    }
    // update dH and g

    PERIODIC_SYNC(dH);	 
    //phgExportVTK(g, "output/dH.vtk", dH, NULL);

    {
	FLOAT a[2] = {dH_max, -dH_min}, b[2]; 
    //phgPrintf("dH_min: %e\n", dH_min);
	MPI_Allreduce(&a, &b, 2, PHG_MPI_FLOAT, MPI_MAX, g->comm);
	dH_max = b[0];
	dH_min = -b[1];
    //phgPrintf("b1: %e\n", b[1]);
	//phgPrintf("   dH: [%8.4e, %8.4e]\n", dH_min, dH_max);
    }

    return;
}






/* ================================================================================
 *
 *   update height using FV (cell based)
 *   
 * ================================================================================
 */
static void 
get_moved_coord_FV_cell(NSSolver *ns, int nstep)
{
    GRID *g = ns->g;
    DOF *u = ns->u[1];
    LAYERED_MESH *gL = ns->gL;
    int i, j, k;
    DOF *dH = ns->dH;
    BOTTOM_FACE *fb;
    FLOAT dt = ns_params->dt0;

    /* ------------------------------------------------------------
     * 
     * Step 1. compute flux
     *
     * ------------------------------------------------------------ */
    FLOAT h_bar[gL->ntria], 
	h_bar_local[gL->nface_bot],
	h_bar_global[gL->ntria],
	height[gL->nvert];
    FLOAT dh_local[gL->nface_bot], max_dht = -1e10;
    
    fb = gL->face_bot;
    //e = NULL;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	FLOAT flux = 0.;
	TRIA *t = gL->trias + fb->tria;
	FLOAT area_tria = t->area;
	int ii;
	
	for (ii = 0; ii < fb->ne; ii++) {
	    SIMPLEX *e = fb->elems[ii];
	    INT *faces = fb->faces + ii*3;

	    for (j = 0; j < 3; j++) { /* tria edges */
		FLOAT *n = t->normals + j*2;
		FLOAT nt[3] = {n[0], n[1], 0.};
		FLOAT lambda[Dim+1];
		int q, v0, v1, v2;
		int s = faces[j];

		if (s >= 0) {
		    assert(s < NFace);
		} else {
		    continue;
		}
		
		QUAD *quad = phgQuadGetQuad2D(1);
		FLOAT area = phgGeomGetFaceArea(g, e, s);
		lambda[s] = 0;
		v0 = GetFaceVertex(s, 0);
		v1 = GetFaceVertex(s, 1);
		v2 = GetFaceVertex(s, 2);
		
		const FLOAT *p, *w;
		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    FLOAT vu[Dim];
		    FLOAT x, y, z;
		    lambda[v0] = *(p++);
		    lambda[v1] = *(p++);
		    lambda[v2] = *(p++);

		    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		    phgDofEval(u, e, lambda, vu);
		    flux += area*(*w) * INNER_PRODUCT(nt, vu);
		    w++;
		} /* end quad */
	    }	  /* end tri edges */
	}	  /* end elements */


	dh_local[i] = flux / LEN_SCALING / area_tria; /* dh per year */
	max_dht = MAX(fabs(dh_local[i]), max_dht);
    }		  /* end bot face */

#   warning CFL limit
    /*
     * 
     * CFL limit: dh < 5e-3
     *
     * */
    FLOAT dht = max_dht;
    MPI_Allreduce(&dht, &max_dht, 1, PHG_MPI_FLOAT, MPI_MAX, g->comm);
    while (dht * dt > 5e-3)
	dt  /= 2.;
    phgPrintf("  time step change to %f\n", dt);
    ns->dt[0] = dt;

    fb = gL->face_bot;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	FLOAT flux = 0.;
	TRIA *t = gL->trias + fb->tria;
	FLOAT area_tria = t->area;
	int ii;
	
	FLOAT vol = gL->volumes[fb->tria] - dh_local[i] * dt * area_tria;
	phgInfo(0, " volume[%3d]: %e, flux: %e\n",
		fb->tria,
		gL->volumes[fb->tria],
		flux * dt / LEN_SCALING);

	if (vol < HEIGHT_EPS * area_tria) {
	    vol = HEIGHT_EPS * area_tria; /* Fix height */
	    phgInfo(0, "*  fix height[%3d] e: %e\n", fb->tria, vol / area_tria);
	}
	gL->volumes[fb->tria] = vol; /* only local tria
				      * volume */
	h_bar_local[i] = vol / area_tria; 
	phgInfo(0, "  h_bar_local[%3d]: %e\n", i, h_bar_local[i]);
    }
    //printf("%3d %e\n", i, flux);





    /* ------------------------------------------------------------
     * 
     * Step 2.  gather h_bar
     *
     * ------------------------------------------------------------ */
    INT   id_tria_global[gL->ntria],
	id_tria_local[gL->nface_bot];
    int rank, *tria_cnts, *tria_dsps, *tria_idxs, ntotal;

    tria_cnts = phgCalloc(2 * g->nprocs, sizeof(*tria_cnts));
    tria_dsps = tria_cnts + g->nprocs;

    MPI_Allgather(&gL->nface_bot, 1, MPI_INT, 
		  tria_cnts, 1, MPI_INT, g->comm);
    ntotal = 0;
    for (rank = 0; rank < g->nprocs; rank++) {
	tria_dsps[rank] = ntotal;
	ntotal += tria_cnts[rank];
    }

    fb = gL->face_bot;
    //e = NULL;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	id_tria_local[i] = fb->tria;
    }

    MPI_Allgatherv(id_tria_local,  gL->nface_bot, MPI_INT,
		   id_tria_global, tria_cnts, tria_dsps, 
		   MPI_INT, g->comm);
    MPI_Allgatherv(h_bar_local,    gL->nface_bot, PHG_MPI_FLOAT,
		   h_bar_global,   tria_cnts, tria_dsps, 
		   PHG_MPI_FLOAT, g->comm);

    for (i = 0; i < gL->ntria; i++) {
	INT id_tria = id_tria_global[i];
	h_bar[id_tria] = h_bar_global[i];
    }

    /* ------------------------------------------------------------
     * 
     * Step 3.  project from trigs to verts
     *
     * ------------------------------------------------------------ */
    
#if USE_PETSC
    /* \sum_T (h, \phi)_T = \sum_T (h^bar, \phi)_T  */
    {
	Vec x, b, u;
	Mat A;
	KSP ksp;
	PC pc;
	PetscReal norm;
	PetscErrorCode ierr;
	PetscInt i, j, k, n = gL->nvert, col[3], its;
	PetscMPIInt size;
	PetscScalar one = 1.0, value[3][3], r[3];
	int rank = 0;

	//PetscInitialize(&argc, &args, (char *) 0, help);
	MPI_Group orig_group, new_group; 
	MPI_Comm new_comm; 
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group); 
	MPI_Group_incl(orig_group, 1, &rank, &new_group);
	MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm); 

	if (g->rank == 0) {

	    MPI_Comm_size(new_comm, &size);
	    if (size != 1)
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1
		SETERRQ(1, "This is a uniprocessor example only!");
#else
		SETERRQ(MPI_COMM_WORLD, 1, "This is a uniprocessor example only!");
#endif

	    /* Create vectors. */
	    ierr = VecCreate(new_comm, &x);	CHKERRQ(ierr);
	    ierr = PetscObjectSetName((PetscObject) x, "Solution");	CHKERRQ(ierr);
	    ierr = VecSetSizes(x, PETSC_DECIDE, n);	CHKERRQ(ierr);
	    ierr = VecSetFromOptions(x);	CHKERRQ(ierr);
	    ierr = VecDuplicate(x, &b);	CHKERRQ(ierr);
	    ierr = VecDuplicate(x, &u);	CHKERRQ(ierr);

	    /* Create matrix. */
	    ierr = MatCreate(new_comm, &A);	CHKERRQ(ierr);
	    ierr = MatSetType(A,MATSEQAIJ);     CHKERRQ(ierr);
	    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);	CHKERRQ(ierr);
	    //ierr = MatSetFromOptions(A);	CHKERRQ(ierr);

	    /* Assemble matrix */
	    TRIA *t = gL->trias;
	    for (k = 0; k < gL->ntria; k++, t++) {
	    	/* nodes of tria */
	    	PetscInt I[3] = {t->verts[0],
	    			 t->verts[1],
	    			 t->verts[2]};
		FLOAT *coord[3];
	    	PetscScalar area = t->area;

		for (i = 0; i < 3; i++) {
		    coord[i] = gL->verts[I[i]];
		}
		FLOAT x, y;
		x = (coord[0][0]+ coord[1][0] + coord[2][0]) / 3.;
		y = (coord[0][1]+ coord[1][1] + coord[2][1]) / 3.;

		/* New height
		 * Vialov
		 * */
		if (0) {
		    x -= 750.;
		    y -= 750.;
		    FLOAT d = sqrt(x*x + y*y);
		    FLOAT h0 = 300.;
		    FLOAT L = 500.;
		    FLOAT n = POWER_N;

		    if (d < L)
			h_bar[k] = h0 * pow(fabs(1 - pow(d/L, (n+1)/n)), 
					    n/(2*n+2));
		    else
			h_bar[k] = HEIGHT_EPS;
		}

	    	for (i = 0; i < 3; i++) {
	    	    for (j = 0; j < 3; j++)
	    		value[i][j] = area * ((i == j) ? (1./6.) : (1./12.));
	    	    r[i] = area * h_bar[k] * (1./3.);
	    	}

	    	for (i = 0; i < 3; i++) {
	    	    ierr = MatSetValues(A, 1, I+i, 3, I,
	    				value[i], ADD_VALUES);  CHKERRQ(ierr);
	    	}

	    	ierr = VecSetValues(b, 3, I, r, ADD_VALUES);  CHKERRQ(ierr);
	    }

	    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);
	    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

	    /* PetscViewer viewer; */
	    /* PetscViewerASCIIOpen(PETSC_COMM_WORLD,"A_.m",&viewer); */
	    /* PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); */
	    /* ierr = MatView(A, viewer); CHKERRQ(ierr); */
	    /* PetscViewerASCIIOpen(PETSC_COMM_WORLD,"b_.m",&viewer); */
	    /* PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); */
	    /* ierr = VecView(b, viewer);CHKERRQ(ierr); */

	    ierr = VecSet(u, 0);	 CHKERRQ(ierr);

	    /* KSP */
	    ierr = KSPCreate(new_comm, &ksp);	CHKERRQ(ierr);
	    //ierr = KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN);	CHKERRQ(ierr);
	    ierr = KSPGetPC(ksp, &pc);	CHKERRQ(ierr);
	    ierr = PCSetType(pc, PCJACOBI);	CHKERRQ(ierr);
	    ierr =
		KSPSetTolerances(ksp, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT,
				 PETSC_DEFAULT);	CHKERRQ(ierr);
	    ierr = KSPSetFromOptions(ksp);	CHKERRQ(ierr);
	    ierr = KSPSolve(ksp, b, x);	CHKERRQ(ierr);

	    PetscScalar *xx;
	    phgInfo(0, "New height\n");
	    ierr = VecGetArray(x, &xx); CHKERRQ(ierr);
	    for (i = 0; i < gL->nvert; i++) {
		height[i] = xx[i];
		if (height[i] < HEIGHT_EPS) {/* Fix height */
		    phgInfo(0, "*   fix height[%3d] v %e\n", i, height[i]);
		    height[i] = HEIGHT_EPS;
		}
		phgInfo(0, "  %3d %e\n", i, height[i]);
	    }

	    /* Free */
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1
	    ierr = VecDestroy(x);	CHKERRQ(ierr);
	    ierr = VecDestroy(u);	CHKERRQ(ierr);
	    ierr = VecDestroy(b);	CHKERRQ(ierr);
	    ierr = MatDestroy(A);	CHKERRQ(ierr);
	    ierr = KSPDestroy(ksp);	CHKERRQ(ierr);
#else
	    ierr = VecDestroy(&x);	CHKERRQ(ierr);
	    ierr = VecDestroy(&u);	CHKERRQ(ierr);
	    ierr = VecDestroy(&b);	CHKERRQ(ierr);
	    ierr = MatDestroy(&A);	CHKERRQ(ierr);
	    ierr = KSPDestroy(&ksp);	CHKERRQ(ierr);
#endif
	}
    }
#else
    phgError(0, "Petsc solver disabled!!!\n");
#endif

    /* ------------------------------------------------------------
     * 
     * Step 4.  update height
     *
     * ------------------------------------------------------------ */
    MPI_Bcast(height, gL->nvert, PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);    

    int ii;
    FLOAT dH_max = -1e20, dH_min = 1e20;
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	int iv = gL->vert_L2S[i];
	assert(nv > 0);
	
	FLOAT dh, dh0, dh1;
	dh0 = 0;
	//dh1 = *DofVertexData(dH, iL[nv-1]);
	dh1 = height[iv] - g->verts[iL[nv-1]][Z_DIR];

# warning set zero height to eps
	/* Set zero height to eps,
	 * dH is not used in this case. */
	if (TRUE) {
	    FLOAT *z = &g->verts[iL[nv-1]][Z_DIR];
	    if (*z + dh1  < HEIGHT_EPS) 
		dh1 = HEIGHT_EPS - *z; /* use H_eps */
	}

	dh = (dh1 - dh0) / nv;

	dH_max = MAX(dH_max, dh1);
	dH_min = MIN(dH_min, dh1);

	for (j = 0; j < nv; j++) {
	    *DofVertexData(dH, iL[j]) = j * dh;
	}
    }
    
    {
        FLOAT a[2] = {dH_max, -dH_min}, b[2];
        MPI_Allreduce(&a, &b, 2, PHG_MPI_FLOAT, MPI_MAX, g->comm);
        dH_max = b[0];
        dH_min = -b[1];
        phgPrintf("   dH: [%8.4e, %8.4e]\n", dH_min, dH_max);
    }
    phgFree(tria_cnts);

}








/* ================================================================================
 *
 *   update height using FV  (point based)
 *   
 * ================================================================================
 */
static void
get_moved_coord_FV_point(NSSolver *ns, int tstep)
{
    GRID *g = ns->g;
    DOF *u = ns->u[1];
    LAYERED_MESH *gL = ns->gL;
    int i, j, k;
    DOF *dH = ns->dH;
    BOTTOM_FACE *fb;
    FLOAT dt = ns_params->dt0;
    INT nL = gL->max_nlayer;

    

    /* ------------------------------------------------------------
     * 
     * Step 1. compute flux
     *
     * ------------------------------------------------------------ */
    FLOAT max_dht = -1e10;
    FLOAT *U_local, *U_vert, *U; 
    U       = phgCalloc(gL->ntria * 6, sizeof(*U)); /* [ne][3][2] */
    U_local = phgCalloc(gL->ntria * 6, sizeof(*U)); /* [ne][3][2] */
    U_vert  = phgCalloc(gL->nvert * 2, sizeof(*U)); /* [nv][2] */

    fb = gL->face_bot;
    //e = NULL;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	TRIA *t = gL->trias + fb->tria;
	FLOAT area_tri = t->area;
	int ii;

	int e2v_tri[3][2] = {
	    {0, 1}, {1, 2}, {2, 0}
	};
	INT vert_tri[3] = {t->verts[0],
			   t->verts[1],
			   t->verts[2]};
	FLOAT *x_tri[3] = {gL->verts[vert_tri[0]], 
			   gL->verts[vert_tri[1]], 
			   gL->verts[vert_tri[2]]};
	FLOAT mid_tet[NEdge][Dim];
	FLOAT lens[3] = {0, 0, 0};
	FLOAT *u_local = U_local + fb->tria * 6; /* [3][2] */

	for (ii = 0; ii < fb->ne; ii++) { /* tets */
	    SIMPLEX *e = fb->elems[ii];
	    const FLOAT *vu[NEdge];
	    int *edge2to3 = fb->edge2to3 + ii*6; /* [3][2] */
	    
	    /* tetra edge */
	    for (j = 0; j < NEdge; j++) {
		vu[j] = DofEdgeData(ns->u[1], e->edges[j]);
		
		int v0 = GetEdgeVertex(j, 0); /* 3D */
		int v1 = GetEdgeVertex(j, 1);
		FLOAT *x0 = g->verts[e->verts[v0]];
		FLOAT *x1 = g->verts[e->verts[v1]];
		
		mid_tet[j][0] = .5 * (x0[0] + x1[0]);
		mid_tet[j][1] = .5 * (x0[1] + x1[1]);
		mid_tet[j][2] = .5 * (x0[2] + x1[2]);
	    }

	    for (k = 0; k < 3; k++) { /* tri edge */
		int e0 = edge2to3[k*2];
		int e1 = edge2to3[k*2+1];
		if (e1 == -1)
		    continue;
		
		/* flux v0 -> v1 */
		FLOAT len = fabs(mid_tet[e0][2] - mid_tet[e1][2]); /* z dir */

		FLOAT u = .5 * (vu[e0][0] + vu[e1][0]) * len; /* \int u dz */
		FLOAT v = .5 * (vu[e0][1] + vu[e1][1]) * len;

		u_local[2*k   ] += u / LEN_SCALING;
		u_local[2*k +1] += v / LEN_SCALING;  /* unit km/a */

		lens[k] += len;
	    } /* end tri edge */
	}     /* end tet */

	for (k = 0; k < 3; k++) { /* velocity mean */
	    u_local[2*k   ] /= lens[k];
	    u_local[2*k +1] /= lens[k];
	}

	phgInfo(3, "U_local[%d * 6]: %e %e %e %e %e %e\n", i,
		u_local[0], u_local[1], u_local[2],
		u_local[3], u_local[4], u_local[5]);
    }	      /* end bot face */


    /* Reduce to root
     * (actually non is sumed)
     * */
    MPI_Reduce(U_local, U, gL->ntria * 6,
	       PHG_MPI_FLOAT, MPI_SUM, 0, g->comm);
    for (i = 0; i < gL->ntria; i++) {
	FLOAT *u = U + i * 6;
	phgInfo(3, "U[%d * 6]: %e %e %e %e %e %e\n", i, 
		u[0], u[1], u[2], u[3], u[4], u[5]);
    }


    /*
     * FV solve
     * */
    if (g->rank == 0) {
	FLOAT *dht, max_dht = 0;
	dht = phgCalloc(gL->nvert, sizeof(*dht)); 
	fv_update(gL->height, U, dht, U_vert);

	/* CFL limit */
	for (i = 0; i < gL->nvert; i++)
	    if (max_dht < fabs(dht[i]))
		max_dht = fabs(dht[i]);
	phgPrintf("   Max dht %e\n", max_dht);

#warning CFL condition 1e-2
	while (max_dht * dt > 1e-1)
	    dt  /= 2.;
	phgPrintf("   Time step change to %f\n", dt);
	ns->dt[0] = dt;

    phgPrintf("is here!\n");
	for (i = 0; i < gL->nvert; i++) {
    phgPrintf("is here 111!\n");
	    FLOAT h = gL->height[i] + dht[i] * dt;

	    if (h < HEIGHT_EPS) {
		h = HEIGHT_EPS; /* Fix height */
		phgInfo(0, "*  fix height[%3d] e: %e\n", i, h);
	    }

	    gL->height[i] = h;
	}
	
	/* debug */
	if (0) {
	    char name[100];
	    sprintf(name, "./output/H_%05d_.m", tstep);
	    FILE *fp = fopen(name, "w");
	    fprintf(fp, "H = [\n");
	    for (i = 0; i < gL->nvert; i++) {
		fprintf(fp, "%e %e %e %e\n",
			gL->height[i], dht[i],
			U_vert[2*i], U_vert[2*i+1]);
	    }
	    fprintf(fp, "];\n");
	    fclose(fp);
	}
	phgFree(dht);
    }
    phgFree(U_local);
    phgFree(U);
    phgFree(U_vert);


    MPI_Bcast(ns->dt, 1,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);    
    MPI_Bcast(gL->height, gL->nvert,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);    

    int ii;
    const FLOAT *ratio = gL->layer_ratio;
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	int iv = gL->vert_L2S[i];
	assert(nv > 0);

	FLOAT h0, h1;
	h0 = 0;
	h1 = gL->height[iv];

	FLOAT H[nv];

	get_layer_height(H, nv, ratio, h0, h1);
    printf("is here!\n");
		
	assert(gL->max_nlayer + 1 == nv);
	for (j = 0; j < nv; j++) {
	    g->verts[iL[j]][Z_DIR] = H[j];
	}
    }
}



void
check_height(GRID *g, LAYERED_MESH *gL)
{
    INT n = gL->nvert, i;
    FLOAT sum_h = 0., de; 
    
    for (i = 0; i < n; i++){
	sum_h += gL->height[i] * gL->height[i];
    }
    sum_h = sqrt(sum_h);

    {
	FLOAT a[2] = {sum_h, -sum_h}, b[2];
	MPI_Allreduce(a, b, 2, 
		      PHG_MPI_FLOAT, MPI_MAX, g->comm);
	de = (fabs(b[0] - sum_h)
	      + fabs(b[1] + sum_h) ) / 2;
    }
    phgPrintf("   *height %e +- %e\n", sum_h, de);
    
}




/* Note: Move mesh coord is divied to two step,
 * 1. Compute dH(dZ)
 * 2. Z += dH
 * This is because Z is not periodic, but dH is.
 *
 * */
void
move_mesh(NSSolver *ns)
{
#if 1
    GRID *g = ns->g;
    SIMPLEX *e;
    FLOAT *vc, *vz;
    LAYERED_MESH *gL = ns->gL;
    INT i;

    /* g->verts changed */
    phgGeomInit_(g, TRUE);

    /* Set inactive element masks */
    ForAllElements(g, e) {
	for (i = 0; i < NVert; i++)
	    if (g->verts[i][Z_DIR] > HEIGHT_EPS)
		break;
	if (i < NVert) {	/* at least one non-zero: active */
	    e->region_mark = 0;
	} else {		/* all zero: inactive */
	    e->region_mark = 1;
	}
    }


    //phgDofDump(dz);
    //phgExportMedit(g, "moved_mesh.mesh");
    phgDofSetDataByFunction(ns->coord, func_xyz_);
    get_height_depth(ns);
    

    /* Volume */
    FLOAT Vol = 0;
    ForAllElements(g, e) {
	Vol += phgGeomGetVolume(g, e);
    }
    FLOAT tmp = Vol;
    MPI_Allreduce(&tmp, &Vol, 1, MPI_DOUBLE, MPI_SUM, g->comm);
    phgPrintf("   Total volume: %24.12e\n", Vol);
#endif

    return;
}    



/*
 * Build Dof Height and Depth
 *  */
void
get_height_depth(NSSolver *ns)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    DOF *depth_P1, *depth_P2; 
    FLOAT z0, z1;
    int i, j, ii;

    depth_P1 = phgDofNew(g, DOF_P1, 1, "depth P1", DofNoAction);
    depth_P2 = phgDofNew(g, ns_params->T_type, 1, "depth P2", DofNoAction);

    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	assert(nv > 0);

	z0 = g->verts[iL[nv-1]][Z_DIR];
	z1 = g->verts[iL[0]][Z_DIR];

	for (j = 0; j < nv; j++) {
	    *DofVertexData(depth_P1, iL[j])
		= z0 - g->verts[iL[j]][Z_DIR];
	}
    }

    /* interp */
    phgDofCopy(depth_P1, &depth_P2, NULL, "depth P2");
    //phgExportVTK(g, "depth.vtk", depth_P1, depth_P2, NULL);

    if (ns->depth_P1 != NULL)
	phgDofFree(&ns->depth_P1);
    if (ns->depth_P2 != NULL)
	phgDofFree(&ns->depth_P2);

    ns->depth_P1 = depth_P1;
    ns->depth_P2 = depth_P2;

    return;
}


/*
 * Update Dof value when mesh coords change.
 * NOT used in the simulation, since the ice-sheet movement is small.
 *
 * Use algorithm in Y. Di, R. Li, T. Tang, and P. Zhang. Moving mesh
 * finite element methods for the incompressible Navier-Stokes
 * equations. SIAM J. Sci. Comput., 26(3):1036–1056, 2005
 *
 *  */
void 
move_dof(GRID *g, DOF *dz, DOF *u)
{
    SIMPLEX *e;
    MAP *map;
    MAT *mat;
    SOLVER *solver;
    int dim = u->dim;
    INT i, j, k;
    DOF *gradu;
    FLOAT msl = 1.;
    BTYPE DB_mask_save[100];
    int verb;

    gradu = phgDofGradient(u, NULL, NULL, NULL);


    /* Note: save DB_mask of dofs, and currently set it to UNDEFINED
     *       to remove the effect of DIRICHLET bdry nodes. */
    if (u->DB_masks != NULL) {
	memcpy(DB_mask_save, u->DB_masks, dim * sizeof(*DB_mask_save));
	memset(u->DB_masks, dim, sizeof(*DB_mask_save));
    } else
	DB_mask_save[0] = u->DB_mask;

    /* Create dof update solver */
    map = phgMapCreate(u, NULL);
    mat = phgMapCreateMat(map, map);
    phgOptionsPush();
    phgOptionsSetOptions("-solver petsc "
			 "-solver_maxit 10000 "
			 "-solver_rtol 1e-10");
    verb = phgVerbosity;
    phgVerbosity = 1;
    solver = phgMat2Solver(SOLVER_DEFAULT, mat);
    phgVerbosity = verb;
    phgVecDisassemble(solver->rhs);
    phgOptionsPop();

    ForAllElements(g, e) {
	int M = DofGetNBas(u, e);	/* number of basises in the element */
	int order = DofTypeOrder(u, e) * 2;
	FLOAT A[M][dim][M][dim], vol, buf[M], rhs[M][dim];
	INT I[M][dim], J[dim][M];
	QUAD *quad;
	int q;
	const FLOAT *w, *p, *gu, *vu, *vz;

	bzero(A, sizeof(A));
	bzero(rhs, sizeof(rhs));
	bzero(buf, sizeof(buf));
	vol = phgGeomGetVolume(g, e);
	quad = phgQuadGetQuad3D(order);
	gu = phgQuadGetDofValues(e, gradu, quad); 
	vu = phgQuadGetDofValues(e, u, quad);        
	vz = phgQuadGetDofValues(e, dz, quad);     
	    
	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    /* Mat */
	    for (i = 0; i < M; i++) {
		for (j = 0; j < M; j++) {
		    const FLOAT *gi = phgQuadGetBasisValues(e, u, i, quad) + q; /* phi_i */
		    const FLOAT *gj = phgQuadGetBasisValues(e, u, j, quad) + q; /* phi_j */
		    FLOAT m = vol*(*w) * (*gi) * (*gj);
		    for (k = 0; k < dim; k++)
			A[i][k][j][k] += m;
		}
	    }

	    for (i = 0; i < M; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, u, i, quad) + q; /* phi_i */
		for (k = 0; k < dim; k++) {
		    rhs[i][k] += vol*(*w) * (*gi) * (vu[k] + msl * gu[k*Dim+Z_DIR] * *vz
						     );
		}
	    }
	    
	    vu += dim; 
	    gu += Dim*dim; 
	    vz++;
	    w++; p += Dim + 1;
	}

	for (i = 0; i < M; i++)
	    for (k = 0; k < dim; k++)
		J[k][i] = I[i][k] =
		    phgMapE2L(solver->rhs->map, 0, e, i * dim + k);

	for (i = 0; i < M; i++) {
	    for (k = 0; k < dim; k++) {
		if (phgDofDirichletBC_(u, e, i*dim+k, NULL, buf, &rhs[i][0],
				       DOF_PROJ_NONE)) {
		    phgMatAddEntries(mat, 1, I[i] + k, M, J[k], buf);
		} else {
		    phgMatAddEntries(mat, 1, I[i] + k, M*dim, I[0],
				     &(A[i][k][0][0]));
		}
	    }
	}

	phgVecAddEntries(solver->rhs, 0, M * dim, I[0], &rhs[0][0]);
    } /* end build rhs */

    phgSolverSolve(solver, FALSE, u, NULL);
    phgPrintf("   Dof %s moving solve: nits %d, res %e\n", 
	      u->name, solver->nits, solver->residual);
    phgDofFree(&gradu);

    phgSolverDestroy(&solver);
    phgMatDestroy(&mat);
    phgMapDestroy(&map);

    /* restore */
    if (u->DB_masks != NULL)
	memcpy(u->DB_masks, DB_mask_save, dim * sizeof(*DB_mask_save));
    else
	u->DB_mask = DB_mask_save[0];
    return;
}

BTYPE
phgDofGetBoundaryType(DOF *u, INT index)
/* returns boundary type of the unknown with local index 'index' */
{
    INT n;
    DOF_TYPE *type = u->type;

    if (type != NULL && type->base_type != NULL) {
	/* Special case: type = DGn */
	SIMPLEX *e = NULL;
	assert(type->np_vert == 0 && type->np_edge == 0 && type->np_face == 0);
	n = index / u->count_elem;	/* the element index */
	if ((e = u->g->elems[n]) == NULL)
	    return UNREFERENCED;
	return phgDofGetElementBoundaryType(u, e, index % u->count_elem);
    }

    if (index < (n = DofGetVertexDataCount(u)))
	return u->g->types_vert[index / u->count_vert];
    index -= n;

    if (index < (n = DofGetEdgeDataCount(u)))
	return u->g->types_edge[index / u->count_edge];
    index -= n;

    if (index < (n = DofGetFaceDataCount(u)))
	return u->g->types_face[index / u->count_face];
    index -= n;

    assert(index < DofGetElementDataCount(u));
    return u->g->types_elem[index / u->count_elem];
}
