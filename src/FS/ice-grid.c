#include "phg.h"
#include "ins.h"
#include "mat_op3.h"
#include <string.h>
#include <math.h>
#include "layers.h"

#define _nsp (ns->ns_params)

/*
 * Init ice grid using ice height.
 * */
void
ice_grid(GRID *g)
{
    SIMPLEX *e;
    INT i, k, s;
    DOF *coord = phgDofNew(g, DOF_P1, Dim, "coord", func_ice_slab);
    DOF *ice_shelf_mask = phgDofNew(g, DOF_P1, 1, "ice_shelf_mask", func_ice_shelf_mask);
    phgExportVTK(g, "ice_shelf_mask.vtk", ice_shelf_mask, NULL);
    FLOAT *v0, *v1, Area[5];
    const FLOAT eps = 1e-8;
    
    Unused(k);
    Unused(eps);
    bzero(Area, sizeof(Area));

    /* ----------------------------------------
     * 
     *   Update boundary type
     * 
     * ---------------------------------------- */

    if (ns_params->update_bdry_type) {
	/* Update boundary type by 2D mesh */
	phgPrintf("Update boundary type by 2D mesh.\n");
	ForAllElements(g, e) {
	    for (s = 0; s < NFace; s++) {
		FLOAT a, n_top[] = {0, 0, 1}, n_bottom[] = {0, 0, -1}, n_terminus[] = {1, 0, 0}, n_divide[] = {-1, 0, 0};
		const FLOAT *n;

		e->bound_type[s] &= ~DIRICHLET;

		if (e->neighbours[s] != NULL) 
		    continue;

		n = phgGeomGetFaceOutNormal(g, e, s);
		a = phgGeomGetFaceArea(g, e, s);
		if (fabs(1 - INNER_PRODUCT(n, n_top)) < eps) {
		    e->bound_type[s] = BC_TOP;    /* clear other */
		    Area[0] += a;
		} 
        else if (fabs(1 - INNER_PRODUCT(n, n_bottom)) < eps) {
		    e->bound_type[s] = BC_BOTTOM; /* clear other */
		    
            if (ns_params->add_ice_shelf){

                    e->bound_type[s] = (BC_ISHELF | BC_BOTTOM);
               }
            else
                    e->bound_type[s] = (BC_BOTTOM_GRD | BC_BOTTOM);
               


            
		    Area[1] += a;
		} 
        else{
		    if (g->period == NULL) {
			e->bound_type[s] = BC_LATERL;

            if (0){
            if (ns_params->add_ice_shelf){

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                INT local_idx0 = e->verts[v0];
                INT local_idx1 = e->verts[v1];
                INT local_idx2 = e->verts[v2];

                if (((g->verts[local_idx0][0] < 2210) && (g->verts[local_idx0][1] > 7.2e2)) || ((g->verts[local_idx0][0] < 2184) && (g->verts[local_idx0][1] <=7.2e2))){

                    /*
                if (*DofVertexData(ice_shelf_mask, local_idx0) > 0 && 
                    *DofVertexData(ice_shelf_mask, local_idx1) > 0 && 
                    *DofVertexData(ice_shelf_mask, local_idx2) > 0) {
                    */
                
                    e->bound_type[s] = (BC_LATERL_GRD | BC_LATERL);
                }
                else
                    e->bound_type[s] = (BC_ISHELF | BC_LATERL);
            }
            }
            
			Area[2] += a;
		    }

		    if (e->bound_type[s] & DIRICHLET)
			Area[3] += a;
		    else 
			Area[4] += a;
		}
	    }
	}
	phgUpdateBoundaryTypes(g);
    }
    else {
	
	ForAllElements(g, e) {
	    for (s = 0; s < NFace; s++) {
		FLOAT a = phgGeomGetFaceArea(g, e, s);
		if (e->bound_type[s] & BC_TOP)
		    Area[0] += a;
		else if (e->bound_type[s] & BC_BOTTOM)
		    Area[1] += a;
		else if (e->bound_type[s] & BC_LATERL)
		    Area[2] += a;
		else if (e->bound_type[s] & DIRICHLET)
		    Area[3] += a;
		else
		    Area[4] += a;
	    }
	}
    }
	    

    phgPrintf("--------------------\n");
    phgPrintf("Set up boundary type\n");
    phgPrintf("  Area top   : %lf\n", Area[0]);
    phgPrintf("  Area bottom: %lf\n", Area[1]);
    phgPrintf("  Area laterl: %lf\n", Area[2]);
    phgPrintf("  Area dirich: %lf\n", Area[3]);
    phgPrintf("  Area other : %lf\n", Area[4]);


    /* ----------------------------------------
     * 
     *   Map to physical region
     * 
     * ---------------------------------------- */
    ForAllElements(g, e) {
	for (i = 0; i < NVert; i++) {
	    v0 = g->verts[e->verts[i]];
        memcpy(v0, DofVertexData(coord, e->verts[i]), Dim*sizeof(*v0));
	}
    }

    phgGeomInit_(g, TRUE);

    phgDofFree(&coord);
    phgDofFree(&ice_shelf_mask);
    return;
}


/* Partitioner */
BOOLEAN 
iceParter(GRID *g, int nprocs, DOF *weights, FLOAT power)
{
    /* e->mark is some where else, e.g. set by ice-grid */
    return TRUE;
}




/* -------------------------------------------------------------------- 
 *
 *
 *  L2 projection of gradient 
 * 
 *                                                                      
 * -------------------------------------------------------------------- */

static MAT *matGu = NULL;
static MAP *Gu_map = NULL;
static DOF *Gradu = NULL;
static SOLVER *solver_Gu = NULL;

static void phgNSInitSolverGu(NSSolver *ns)
{
    int verb = phgVerbosity;
    
    /* dof map */
    Gu_map = phgMapCreate(Gradu, NULL);

    /* matrices */
    matGu = phgMapCreateMat(Gu_map, Gu_map);
    matGu->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;

    /* solver_Gu */
    phgOptionsPush();
    phgOptionsSetOptions("-solver petsc "
			 "-solver_maxit 10000 "
			 "-solver_rtol 1e-10");
    //phgOptionsSetOptions(Gu_opts);
    phgVerbosity = 0;
    solver_Gu = phgMat2Solver(SOLVER_DEFAULT, matGu);
    phgVerbosity = verb;
    phgOptionsPop();
    phgVecDisassemble(solver_Gu->rhs);
    phgDofSetDataByValue(Gradu, 0.);
    return;
}

static void phgNSDestroySolverGu(NSSolver *ns)
{
    phgInfo(2, "   Destroy solver Gu\n");
    phgMatDestroy(&matGu);
    phgSolverDestroy(&solver_Gu);
    solver_Gu = NULL;
    phgMapDestroy(&Gu_map);
}

void 
phgNSBuildSolverGu(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, k, q;

    ForAllElements(g, e) {
	int M = Gradu->type->nbas;	/* num of bases of Velocity */
	int order = DofTypeOrder(Gradu, e) * 2;
	FLOAT A[M][DDim][M][DDim], rhs[M][DDim];
	INT I[M][DDim];
	QUAD *quad;
	FLOAT vol, det;
	const FLOAT *w, *p, *gu;

	Bzero(A); Bzero(rhs);
	quad = phgQuadGetQuad3D(order);
	gu = phgQuadGetDofValues(e, ns->gradu[1], quad); 

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);

	    for (i = 0; i < M; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, Gradu, i, quad) + q;    
		for (j = 0; j < M; j++) {
		    const FLOAT *gj = phgQuadGetBasisValues(e, Gradu, j, quad) + q;       
		    FLOAT qmass = vol*(*w) * (*gj) * (*gi);
		    for (k = 0; k < DDim; k++)
			A[i][k][j][k] += qmass;
		}
		    
		for (k = 0; k < DDim; k++)
		    rhs[i][k] += vol*(*w) * gu[k] * (*gi); 
	    }
	    gu += DDim;
	    w++; p += Dim + 1;
	}

	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    for (k = 0; k < DDim; k++)
		I[i][k] = phgMapE2L(matGu->cmap, 0, e, i * DDim + k);

	/* Global res */
	for (i = 0; i < M; i++)
	    phgMatAddEntries(matGu, DDim, I[i], M * DDim, I[0],
			     &(A[i][0][0][0])); 
	phgSolverAddRHSEntries(solver_Gu, M*DDim, I[0], &rhs[0][0]);
    }				/* end element */
    
    phgVecAssemble(solver_Gu->rhs);
    solver_Gu->rhs_updated = FALSE;

    if (DUMP_MAT_VEC) {
	phgPrintf("Dumping MatGu, rhsGu\n");
	phgMatDumpMATLAB(solver_Gu->mat, "A_gu", "mat_gu_.m");
	phgVecDumpMATLAB(solver_Gu->rhs, "b_gu", "rhs_gu_.m");
    }
 
    phgOptionsPush();
    phgOptionsSetOptions(_nsp->Gu_opts);
    phgSolverAssemble(solver_Gu);
    phgOptionsPop();
 
    return;
}


/*
 * Build layer height
 * 
 *  */
void 
get_layer_height(FLOAT *H, int nv, const FLOAT *ratio, FLOAT h0, FLOAT h1)
{
    int i;
    for (i = 0; i < nv; i++) {
	H[i] = h0 + (h1-h0) * ratio[i];
    }


#if 0
#  warning Build struct boundary layer: h[5] <= 0.6km

    const int im = 5;
    const FLOAT Hm = .6;
    const FLOAT H1 = h1;
    assert(im < nv);

    if (H[im] > Hm) {
	for (i = 0; i < im; i++)
	    H[i] = H[i] * Hm / H[im];
	for (i = im+2; i < nv; i++)
	    H[i] = H1 - (H1 - H[i]) * (H1 - Hm) / (H1 - H[im]);
	H[im] = Hm;
    }
#endif

    return;
}


#if 1

/*
 * ----------------------------------------------------------------------
 * Monitor:
 * 1. On bottom and/or top surface,
 *    velocity, tau, dp
 * 2. Surface velocity at [0, L] x [L/4]
 *
 * Ouput individually to files, leave post process script to handle.
 * 
 * */

#define FORMAT "%16.8f "

#define XX 0
#define XY 1
#define XZ 2
#define YX 3
#define YY 4
#define YZ 5
#define ZX 6
#define ZY 7
#define ZZ 8

#define PACK_DAT(dat, p, x, y, z, vu, tau, dp, vT, dh)	\
    dat[p++] = x;				\
    dat[p++] = y;				\
    dat[p++] = z;				\
    dat[p++] = vu[0];				\
    dat[p++] = vu[1];				\
    dat[p++] = vu[2];				\
    dat[p++] = tau[XZ];				\
    dat[p++] = tau[YZ];				\
    dat[p++] = dp;				\
    dat[p++] = *vT;				\
    dat[p++] = *dh;				


#define UNPACK_DAT(dat, p, x, y, z, vu, tau, dp, vT, dh)	\
    x = dat[p++];					\
    y = dat[p++];					\
    z = dat[p++];					\
    vu[0] = dat[p++];					\
    vu[1] = dat[p++];					\
    vu[2] = dat[p++];					\
    tau[XZ] = dat[p++];					\
    tau[YZ] = dat[p++];					\
    dp = dat[p++];					\
    *vT = dat[p++];					\
    *dh = dat[p++];


void
ice_monitor(NSSolver *ns, int nonstep)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    char fname[1000];
    DOF *u = ns->u[1];
    DOF *p = ns->p[1]; 
    DOF *T = ns->T[1]; 
    DOF *dH = ns->dH;
    //DOF *gradu = ns->gradu[1];
    FLOAT *data, *sbuf, *rbuf;
    INT i, ii, ip, nvalue;
    int viscosity_type = ns->viscosity_type;

#if TEST_CASE == ICE_GREEN_LAND
    ERROR_MSG("TODO");
#endif

    /* L2 projection of gradu */
    Gradu = phgDofNew(g, DOF_P1, DDim, "Grad u", DofNoAction);
    phgNSInitSolverGu(ns);
    phgNSBuildSolverGu(ns);
    phgSolverSolve(solver_Gu, FALSE, Gradu, NULL);
    phgPrintf("      solver_Gu: nits = %d, resid = %0.4lg ",
	      solver_Gu->nits, solver_Gu->residual);
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
    phgNSDestroySolverGu(ns);

    /*
     * Gathter Dof data on vert:
     * Surf: x, y, z, u[X|Y|Z], tau[XZ|YZ], dp, T
     * Base: x, y, z, u[X|Y|Z], tau[XZ|YZ], dp, T
     *   nvalue = 2 * (3 + 3 + 2 + 3) = 22;
     * */
    nvalue = 2*11;
    PHG_CALLOC(sbuf, nvalue * gL->nvert_bot);
    PHG_CALLOC(rbuf, nvalue * gL->vert_bot_total);
    
    ip = 0;
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	const FLOAT *vu, *gu, *vT, *dh;
	FLOAT vp, x, y, z0, z1, nu,
	    eu[DDim], tau[DDim], scale, dp, p0;
	INT nv, iL0, iL1;

	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	nv = gL->vert_local_lists[i][0];
	iL0 = gL->vert_local_lists[i][1];
	iL1 = gL->vert_local_lists[i][nv];
	assert(nv > 0);

	x = g->verts[iL0][0];
	y = g->verts[iL0][1];
	z0 = g->verts[iL0][2];
	z1 = g->verts[iL1][2];

	const FLOAT rho = RHO_ICE;
	const FLOAT grav = GRAVITY;
	const FLOAT a = SEC_PER_YEAR;
	
	/* Base */
	vu = DofVertexData(u, iL0);
	gu = DofVertexData(Gradu, iL0);
	vT = DofVertexData(T, iL0);
	vp = *DofVertexData(p, iL0);
	dh = DofVertexData(dH, iL0);
	MAT3_SYM(gu, eu);
	nu = get_effective_viscosity(gu, *vT, 0, viscosity_type);
	MAT3_ZERO(tau);
	scale = nu / LEN_SCALING; 
	MAT3_AXPBY(scale, eu, 0., tau);
	p0 = rho*grav* (z1 - z0) * LEN_SCALING; 
	dp = vp*PRES_SCALING - p0;

	PACK_DAT(sbuf, ip, x, y, z0, vu, tau, dp, vT, dh);

	/* Surf */
	vu = DofVertexData(u, iL1);
	gu = DofVertexData(Gradu, iL1);
	vT = DofVertexData(T, iL1);
	vp = *DofVertexData(p, iL1);
	dh = DofVertexData(dH, iL1);
	MAT3_SYM(gu, eu);
	nu = get_effective_viscosity(gu, *vT, 0, viscosity_type);
	MAT3_ZERO(tau);
	scale = nu / LEN_SCALING;
	MAT3_AXPBY(scale, eu, 0., tau);
	p0 = 0;
	dp = 0;

	PACK_DAT(sbuf, ip, x, y, z1, vu, tau, dp, vT, dh);
    }
    MPI_Datatype type;
    MPI_Type_contiguous(nvalue, PHG_MPI_FLOAT, &type);
    MPI_Type_commit(&type);
    MPI_Gatherv(sbuf, gL->nvert_bot, type,
		rbuf, gL->vert_bot_cnts, gL->vert_bot_dsps, 
		type, 0, g->comm);
    MPI_Type_free(&type);

    if (g->rank == 0) {
	int rank;
	FILE *fp_base, *fp_surf;
	FLOAT *data;
	PHG_CALLOC(data, nvalue * gL->nvert);

	for (rank = 0; rank < g->nprocs; rank++) {
	    INT *idx = gL->vert_bot_idxs
		+ gL->vert_bot_dsps[rank];
	    FLOAT *v = rbuf
		+ nvalue * gL->vert_bot_dsps[rank];

	    for (i = 0;
		 i < gL->vert_bot_cnts[rank];
		 i++, idx++, v += nvalue) {
		memcpy(data + nvalue * (*idx), 
		       v, nvalue * sizeof(FLOAT));
	    }
	}

	/* output */
	phgPrintf("   Output surf & base data\n");
	sprintf(fname, OUTPUT_DIR "/ins_" NS_PROBLEM "_surf.%05d", 
		nonstep);
	fp_surf = fopen(fname, "w");
	sprintf(fname, OUTPUT_DIR "/ins_" NS_PROBLEM "_base.%05d", 
		nonstep);
	fp_base = fopen(fname, "w");

	ip = 0;
	for (i = 0; i < gL->nvert; i++) {
	    FLOAT vu[Dim], x, y, z, tau[DDim], dp, vT, dh;

	    UNPACK_DAT(data, ip, x, y, z, vu, tau, dp, &vT, &dh);
	    fprintf(fp_base, FORMAT FORMAT FORMAT, x, y, z);
	    fprintf(fp_base, FORMAT FORMAT FORMAT, vu[0], vu[1], vu[2]);
	    fprintf(fp_base, FORMAT FORMAT FORMAT FORMAT FORMAT,
		    tau[XZ], tau[YZ], dp, vT, dh);
	    fprintf(fp_base, "\n");

	    UNPACK_DAT(data, ip, x, y, z, vu, tau, dp, &vT, &dh);
	    fprintf(fp_surf, FORMAT FORMAT FORMAT, x, y, z);
	    fprintf(fp_surf, FORMAT FORMAT FORMAT, vu[0], vu[1], vu[2]);
	    fprintf(fp_surf, FORMAT FORMAT FORMAT FORMAT FORMAT,
		    tau[XZ], tau[YZ], dp, vT, dh);
	    fprintf(fp_surf, "\n");

	}
	fclose(fp_surf);
	fclose(fp_base);
	phgFree(data);
    }

    phgFree(sbuf);
    phgFree(rbuf);
    phgDofFree(&Gradu);

    return;
}





#else

#  define MATCH_EDGE(v0, v1, vL0, vL1)		\
    ((v0 == vL0 && v1 == vL1)			\
     || (v0 == vL1 && v1 == vL0)))

void
ice_monitor(NSSolver *ns, int tstep) 
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    FILE *fp;


    char name[1000];
    sprintf(name, OUTPUT_DIR"Bot_%05d.p%03d", tstep, g->rank);
    fp = fopen(name, "r");

    fb = gL->face_bot;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	TRIA *t = gL->trias + fb->tria;
	SIMPLEX *e = fb->elems[0];

	for (s = 0; s < NFace; s++) {
	    if (e->bound_type[s] & BC_BOTTOM) {
		int nbas_face = NbasFace(ns->du);
		SHORT bases[nbas_face];
		phgDofGetBasesOnFace(ns->u[1], e, s, bases);
		

		for (j = 0; j < nbas_face; j++) {
		    const FLOAT *p =
			phgDofGetElementCoordinates(ns->u[1], e, bases[j] * Dim);
		    const FLAOT *lambda =
			phgGeomXYZ2Lambda(g, e, p[0], p[1], p[2]);
		    const FLOAT *vu, *vp, *vT, *vH;
		    
		    phgDofEval(ns->u[1], e, lambda, vu);
		    phgDofEval(ns->p[1], e, lambda, vp);
		    phgDofEval(ns->T[1], e, lambda, vT);
		    phgDofEval(ns->depth_P1, e, lambda, vH);
		    fprintf(fp, "%5d %5d ", fb->tria, bases[j]);
		    fprintf(fp, "%15.8e %15.8e %15.8e ", vu[0], vu[1], vu[2]);
		    fprintf(fp, "%15.8e %15.8e %15.8e\n", *vp, *vT, *vH);
		}
	    } /* bot */
	}     /* elem face */
    }	      /* bot face */

    fclose(fp);

    //e = NULL;
    /* for (i = 0; i < gL->nface_bot; i++, fb++) { */
    /* 	TRIA *t = gL->trias + fb->tria; */
    /* 	SIMPLEX *e = fb->elems[0]; */
    /* 	INT ev_ele[3][2] = {{0, 1}, {1, 2}, {2, 0}}; */
    /* 	INT vL[3];  */
    /* 	INT ev_tri[3][2]; */
    /* 	int i, j, k, s; */

    /* 	for (i = 0; i < 3; i++) */
    /* 	    vL[i] = gL->vert_local_lists[fb->vert[i]][1]; */
    /* 	for (i = 0; i < 3; i++) */
    /* 	    for (k = 0; k < 2; k++) */
    /* 		ev_tri[i][k] = vL[ev_ele[i][k]]; */
	
    /* 	for (s = 0; s < NFace; s++) { */
    /* 	    if (e->bound_type[s] & BC_BOTTOM) { */
    /* 		int v[3]; */
    /* 		INT ev_tet[3][2]; */
    /* 		int edge2to3[3]; */
		
    /* 		for (k = 0; k < 3; k++) */
    /* 		    v[k] = GetFaceVertex(s, k); */
    /* 		for (i = 0; i < 3; i++) */
    /* 		    for (k = 0; k < 2; k++) */
    /* 			ev_tet[i][k] = e->verts[v[ev_ele[i][k]]]; */

    /* 		for (i = 0; i < 3; i++) { */
    /* 		    for (j = 0; j < 3; j++) { */
    /* 			if (MATCH_EDGE(ev_tri[i][0], ev_tri[i][1], */
    /* 				       ev_tet[j][0], ev_tet[j][1])) { */
    /* 			    edge2to3[i] = GetEdgeNo(v[ev_ele[i][0]], */
    /* 						    v[ev_ele[i][1]]); */
    /* 			    break; */
    /* 			} /\* match *\/ */
    /* 		    }/\* tet edge *\/ */
    /* 		    assert(j < 3); */
    /* 		}    /\* tri edge *\/ */

    /* 	    } /\* bot *\/ */
    /* 	}     /\* elem face *\/ */
    }	      /* bot face */
}

#endif
