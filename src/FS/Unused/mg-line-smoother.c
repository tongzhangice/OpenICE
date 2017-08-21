#include "multi-grid.h"
#include "vtk-draw.h"
#define _p mg_params
#define _mgp (mg->mg_params)
/*
 * Multi Grid line smoother
 *
 * */

#define CHECK_MAT_ENTRIES 0

/* Lapack routines */
#define integer int
#define doublereal double
int dgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	    doublereal *ab, integer *ldab, integer *ipiv, integer *info);

int dgbtrs_(char *trans, integer *n, integer *kl, integer *
	    ku, integer *nrhs, doublereal *ab, integer *ldab, integer *ipiv,
	    doublereal *b, integer *ldb, integer *info);
#undef integer
#undef doublereal

#define MAX_BLOCK 16
#define MAX_BLOCK_DOF 1024
#define EPS_H_DIAM 0.1
#define IDFB_SIZE 1000		/* max # of dof_Ln in line,
				 * TODO: dynamic decide */

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define INNER_PRODUCT3(p, q)			\
    (*(p    ) * *(q    ) +			\
     *(p + 1) * *(q + 1) +			\
     *(p + 2) * *(q + 2))

#define INNER_PRODUCT4(p, q)			\
    (*(p    ) * *(q    ) +			\
     *(p + 1) * *(q + 1) +			\
     *(p + 2) * *(q + 2) +			\
     *(p + 3) * *(q + 3))

#define CROSS_PRODUCT3(c, a, b)	{		\
	c[0] = a[1]*b[2] - a[2]*b[1];		\
	c[1] = a[2]*b[0] - a[0]*b[2];		\
	c[2] = a[0]*b[1] - a[1]*b[0];		\
    }

#define NORMALIZE3(v) {				\
	FLOAT norm_ = INNER_PRODUCT3(v, v);	\
	norm_ = sqrt(norm_);			\
	assert(norm_ > 1E-8);			\
	v[0] /= norm_;				\
	v[1] /= norm_;				\
	v[2] /= norm_;				\
    }

#define NORMALIZE4(v) {				\
	FLOAT norm_ = INNER_PRODUCT4(v, v);	\
	norm_ = sqrt(norm_);			\
	assert(norm_ > 1E-8);			\
	v[0] /= norm_;				\
	v[1] /= norm_;				\
	v[2] /= norm_;				\
	v[3] /= norm_;				\
    }

char vtk_tmp_str[1000];

/*
 * index                        : local index on element
 * index /= dim                 : remove dim mode
 * index -= n:                  : next geo type
 * index / type->np_geo         : geo index 
 * index % type->np_geo         : index on geo
 *
 * m = type->np_geo
 * n = m * NGeo                : geo lam span in dof type
 * i = index / type->np_geo    : geo index
 * j = index % type->np_geo    : index on geo
 * 
 *  */

FLOAT *
phgDofGetElementLambda(DOF *u, SIMPLEX *e, int index)
/* returns world lambda coordinates of the unknown with local 'index' on element 'e' */
{
    static FLOAT lam[Dim+1];
    GRID *g = u->g;
    FLOAT *L;
    int m, n, i, j;
    DOF_TYPE *type = u->type;

    assert(type != NULL);

    index /= u->dim;

    if (type->base_type != NULL)
	type = type->base_type;

    if (type->points == NULL)
	phgError(1, "%s:%d, DOF points for DOF \"%s\" undefined.\n",
		 __FILE__, __LINE__, type->name);

    Unused(g);
    bzero(lam, sizeof(lam));

    /* vert */
    if (index < (n = (m = type->np_vert) * NVert)) {
	i = index / m;
	j = index % m;
	lam[i] = 1.;
	return lam;
    }
    index -= n;

    /* edge */
    if (index < (n = (m = type->np_edge) * NEdge)) {
	i = index / m;
	j = index % m;
	L = type->points + type->np_vert * 1 + j * 2;
	lam[GetEdgeVertex(i, 0)] = L[0];
	lam[GetEdgeVertex(i, 1)] = L[1];
	return lam;
    }
    index -= n;
    L = type->points + type->np_vert * 1 + type->np_edge * 2;

    /* face */
    if (index < (n = (m = type->np_face) * NFace)) {
	i = index / m;
	j = index % m;
	L += j * 3;
	lam[GetFaceVertex(i, 0)] = L[0];
	lam[GetFaceVertex(i, 1)] = L[1];
	lam[GetFaceVertex(i, 2)] = L[2];
	return lam;
    }
    index -= n;
    L += type->np_face * 3 + (index % type->np_elem) * 4;

    /* element */
    memcpy(lam, L, (Dim+1) * sizeof(FLOAT));
    return lam;
}


static FLOAT crdx[1000][Dim];

static int 
coord_comp(const void *pa, const void *pb)
{
    int *a = (int *) pa, *b = (int *) pb;  
    FLOAT *xa = crdx[*a], *xb = crdx[*b], dx, dy, dz;
    
    dx = xa[0] - xb[0];
    if (fabs(dx) > 1e-8) {
	if (dx > 0)
	    return TRUE;
	else 
	    return FALSE;
    }

    dy = xa[1] - xb[1];
    if (fabs(dy) > 1e-8) {
	if (dy > 0)
	    return TRUE;
	else 
	    return FALSE;
    }

    dz = xa[2] - xb[2];
    if (dz > 0)
	return TRUE;
    else 
	return FALSE;
}


#define DEBUG_ROW_CONNECTION(row)					\
    while (1)								\
	vtkPause(0);							\
	phgError(1, "Complicicated line connection!!!\n");		



void 
phgMultiGridInitLineBlock(MG_LEVEL *ml, int islv, DOF *u0)
{
    GRID *g = ml->grid;	
    SIMPLEX *e;
    MAT *mat_cnnt;
    MAT *mat_fill;		/* mat_fill keeps all possible matrix entries,
				 * zero or not.
				 * TODO: remove it. */
    DOF *u;
    DOF_TYPE *utype = u0->type;
    int order = utype->order;
    MAP *map = NULL;
    INT i, j, k, l, nlocal, localsize;
    INT *dof_iLine, *dof_iLngh;
    INT *line_dofs, *line_dof0, *line_ndof, *dof0;
    INT *lngh_dofs, *lngh_dof0, *lngh_ndof, *ndofs;
    INT nd, nd_alloc, nl, nl_alloc, ndof, next_dof0;
    INT ndof_pt, *dof_pts, nds[10];
    MAT_ROW *row;
    VEF_MAP2 *vef; 
    MG_BLOCK_DOFS *bk;
    FLOAT line_dir[Dim] = {0, 0., 1.};
    int verb = 3;

    FLOAT *crd;
    DOF *Coord;
    MAP *crd_map;
    VEC *crd_vec;

    Unused(order);
    phgPrintf("   Init Line for dof: %s\n", u0->name);

    /*
     * Line smoother for vector dof:
     * Using same type of scalar dof to build line connection info.
     *  */
    u = phgDofNew(g, utype, 1, "scalar u", DofNoAction);
    map = phgMapCreate(u, NULL);
    nlocal = map->nlocal;
    localsize = map->localsize;

    /*
     * set coord for debug
     * TODO: remove dof coord
     * */
    {
	Coord = phgDofNew(g, utype, Dim, "coord", DofNoAction);
	ForAllElements(g, e) {
	    int N = utype->nbas;
	    FLOAT value[N][Dim];
	    for (i = 0; i < N; i++) {
		const FLOAT *coord = phgDofGetElementCoordinates(Coord, e, i*Dim);
		memcpy(value[i], coord, Dim * sizeof(FLOAT));
	    }
	    phgDofSetElementDatas(Coord, e, value[0]);
	}
	//phgDofDump(Coord);

	/* dump coord to matlab */
	crd_map = phgMapCreate(Coord, NULL); 
	crd_vec = phgMapCreateVec(crd_map, 1);
	phgMapDofToLocalData(crd_map, 1, &Coord, crd_vec->data);
	crd = crd_vec->data;
	/* phgVecDumpMATLAB(crd_vec, "coord", "coord_.m"); */
	/* phgVecDestroy(&cvec); */
	/* phgMapDestroy(&cmap); */
    }	/* end of vtk debug */

    /*
     * step 1. Build connection matrix:
     * connection relationship is determined by dof direciton.
     * 
     * */
    verb = 3;

#if DIRECTION_CONNECTED
    phgInfo(verb, "Use direction to build connection matrix.\n");
#else
    phgInfo(verb, "Use distance ratio to build connection matrix.\n");
#endif

    mat_cnnt = phgMapCreateMat(map, map);
    /* Set mode means connection relationship is local.  Incase
     * connecition is set several times, it's decided only by the last
     * time.
     * 
     * Note:
     * 
     * If we change mat type to add, than inter-proc relationship will
     * be establised. Be carefull then assertion like idx <
     * map->localsize will fail!  Such assertion only means line neigh
     * should be local, fix it when need.
     *
     * */
    phgMatSetMode(mat_cnnt, PHG_REPLACE); 
    mat_fill = phgMapCreateMat(map, map);

    ForAllElements(g, e) {
	int N = utype->nbas;
	INT I[N];
	FLOAT diam, A[N][N], lam_dir[Dim+1], lam_dofs[N][Dim+1];
	const FLOAT *J = NULL, *lam;

	/* Only in B.L. */
#if 0
	if (e->region_mark != 1)
	    continue;
#endif

	bzero(A, sizeof(A));
	/* vtkSetColor(verb, "yellow"); */
	/* vtkSetTransparent(verb, 0.2); */
	/* vtkDrawElement(verb, e); */
	/* vtkSetTransparent(verb, 1.); */

	phgInfo(verb, "---\nelement: %d\n", e->index);
	diam = phgGeomGetDiameter(g, e);
	J = phgGeomGetJacobian(g, e);
	for (k = 0; k < Dim+1; k++) {
	    lam_dir[k] = J[0] * line_dir[0] 
		+ J[1] * line_dir[1] + J[2] * line_dir[2];
	    J += 4;
	}
	NORMALIZE4(lam_dir);
	//SHOW_V(lam_dir, Dim+1);

	for (i = 0; i < N; i++) {
	    lam = phgDofGetElementLambda(u, e, i*u->dim);
	    memcpy(lam_dofs[i], lam, (Dim+1) * sizeof(FLOAT));
	}

	for (i = 0; i < N; i++) 
	    I[i] = phgMapE2L(map, 0, e, i);

#if DIRECTION_CONNECTED
#  warning Line smoother: using line direction
	/*
	 * Method I: Use direction to determine a Dofs line.
	 *
	 * Test direction between every two bases.
	 * TODO: use BLAS
	 * */
	for (i = 0; i < N; i++) {
	    FLOAT lam_test[Dim], lam_diff[Dim+1];
	    //memcpy(lam_diff, lam_dofs, N * (Dim+1) * sizeof(FLOAT));

	    /* diff basis_i basis_j */
	    for (j = 0; j < i; j++) {
		FLOAT theta;
		phgInfo(verb, "\n test %d %d\n", i, j);

		for (k = 0; k < Dim+1; k++) 
		    lam_diff[k] = lam_dofs[j][k] - lam_dofs[i][k];
		//SHOW_V(lam_diff, Dim+1);
		CROSS_PRODUCT3(lam_test, lam_diff, lam_dir);
		//SHOW_V(lam_test, Dim+1);

		/* **************************************
		 * Set rule of Dofs line here:
		 *
		 * Here the direction of vector between two Dofs is used to determine
		 * whether they lie in a line. A parallel test based on cross product
		 * of normalized vector is used.
		 *
		 * For CR1 element a tolerance (THETA_TOL) need to be determined.
		 *
		 * ************************************** */
#define THETA_TOL 1e-5
		theta = INNER_PRODUCT3(lam_test, lam_test);
		if (theta < THETA_TOL) {
		    const FLOAT *coord;
		    FLOAT c_[2][Dim];
		    
		    {
			/* vtk debug */
			coord = phgDofGetElementCoordinates(Coord, e, i*Dim);
			memcpy(c_[0], coord, Dim * sizeof(FLOAT));
			coord = phgDofGetElementCoordinates(Coord, e, j*Dim);
			memcpy(c_[1], coord, Dim * sizeof(FLOAT));
		    
			//vtkTmpActorsBegin(verb);
			vtkSetColor(verb, "red");
			sprintf(vtk_tmp_str, "  %d", i);
			vtkDrawPoint(verb, c_[0], vtk_tmp_str)
			    sprintf(vtk_tmp_str, "  %d", j);
			vtkDrawPoint(verb, c_[1], vtk_tmp_str)
			    vtkDrawLine(verb, c_[0], c_[1]);

			phgInfo(verb, "   found connection: (%d, %d)\n", i, j);
			vtkPause(verb);
			//vtkTmpActorsClear(verb);
		    } /* end of vtk debug */
		    
		    A[i][j] = A[j][i] = 1.; /* bi-direction */
		}			    
	    } /* end of basis_j */
	}     /* end of basis_i */

#elif 0
	/*
	 * Method II: Use length ratio determine a Dofs line.
	 * ratio = dist(phi_i, phi_j) / diam(e)
	 * RATIO_TOL determine need to be decided.
	 *
	 * Warning: this method could fail when the element have big angle(s),
	 *    maybe Method III: Use length is needed.
	 *
	 * */
	FLOAT crdb[N][Dim];
	FLOAT line_ar_tol = mg_params->line_ar_tol;
	for (i = 0; i < N; i++) {
	    const FLOAT *c = phgDofGetElementCoordinates(u, e, i*u->dim);
	    memcpy(crdb[i], c, Dim * sizeof(*c));
	}

	for (i = 0; i < N; i++) {
	    FLOAT length, ratio;
	    /* diff basis_i basis_j */
	    for (j = 0; j < i; j++) {
		FLOAT vdiff[Dim];

		for (k = 0; k < Dim; k++) 
		    vdiff[k] = crdb[i][k] - crdb[j][k];

		length = sqrt(INNER_PRODUCT3(vdiff, vdiff));
		ratio = length / diam * order;
		phgInfo(verb, "   length(%d, %d): %12.6f, ratio: %12.6f\n", 
			i, j, length, ratio);

		if (ratio < line_ar_tol) {
		    const FLOAT *coord;
		    FLOAT c_[2][Dim];

		    coord = phgDofGetElementCoordinates(Coord, e, i*Dim);
		    memcpy(c_[0], coord, Dim * sizeof(FLOAT));
		    coord = phgDofGetElementCoordinates(Coord, e, j*Dim);
		    memcpy(c_[1], coord, Dim * sizeof(FLOAT));
		    
		    //vtkTmpActorsBegin(verb);
		    vtkSetColor(verb, "red");
		    sprintf(vtk_tmp_str, "  %d", i);
		    vtkDrawPoint(verb, c_[0], vtk_tmp_str)
			sprintf(vtk_tmp_str, "  %d", j);
		    vtkDrawPoint(verb, c_[1], vtk_tmp_str)
			vtkDrawLine(verb, c_[0], c_[1]);

		    phgInfo(verb, "   found connection: (%d, %d)\n", i, j);
		    vtkPause(verb);
		    //vtkTmpActorsClear(verb);
		    
		    A[i][j] = A[j][i] = 1.; /* bi-direction */
		}
	    }
	}
#else
	/*
	 *  Other method to determine a Dofs line.
	 *
	 *  Method III: Use cylinder coordinates.
	 *  
	 *
	 *
	 *  */


#endif

	/*
	 * DOF P2:
	 *   reduce connection of line on edge to v-e-v
	 *  */
	if (utype == DOF_P2 || utype == DOF_SUB2) {
	    for (j = 0; j < NEdge; j++) {
		int jj = j + NVert, ii[NVert], i0, i1, nn;
		nn = 0;
		for (i = 0; i < NVert; i++)
		    if (A[jj][i] > .5)
			ii[nn++] = i;
		if (nn != 2)
		    continue;
		i0 = ii[0]; i1 = ii[1];
		A[i0][i1] = A[i1][i0] = 0.;
	    }
	}
	/*
	 * DOF P3:
	 *   reduce connection of line
	 *   1) on edge to v-e-e-v
	 *   2) on face to e-f-e
	 *  */
	if (utype == DOF_P3 || utype == DOF_SUB3) {
	    //SHOW_M_(verb, A[0], N, N);
	    /* edge */
	    for (j = 0; j < NEdge; j++) {
		int j0, j1;
		int M[2], i0, i1;

		GetEdgeVertices(e, j, i0, i1);
		phgDofMapEdgeData(utype, e, j, M);
		j0 = NVert + 2*j + M[0];
		j1 = NVert + 2*j + M[1];
		/*  Note: take edge as whole, doesn't allow single
		 *    connection on an edge.
		 *  */
#define EDGE_CONNECT(i, k0, k1, k2)					\
		(A[i][k0] > .5 || A[i][k1] > .5 || A[i][k2] > .5)
#define CLEAR_CONNECT(i, j, k, l)		\
		A[i][j] = A[j][i] = 0;		\
		A[i][k] = A[k][i] = 0;		\
		A[i][l] = A[l][i] = 0;		\
		A[j][k] = A[k][j] = 0;		\
		A[j][l] = A[l][j] = 0;		\
		A[l][k] = A[k][l] = 0;		\

		if (EDGE_CONNECT(i0, i1, j0, j1) || 
		    EDGE_CONNECT(i1, j0, j1, i0) || 
		    EDGE_CONNECT(j0, j1, i0, i1) || 
		    EDGE_CONNECT(j1, i0, i1, j0)) {/* edge-vert */
		    /* v-e-e-v */
		    CLEAR_CONNECT(i0, i1, j0, j1)
		    A[i0][j0] = A[j0][i0] = 1.;
		    A[j0][j1] = A[j1][j0] = 1.;
		    A[j1][i1] = A[i1][j1] =1.;
		} else {
		    CLEAR_CONNECT(i0, i1, j0, j1)
		}
	    }
#undef EDGE_CONNECT
	    /* face */
	    for (j = 0; j < NFace; j++) {
		int jj = j + NVert + NEdge * 2,
		    ii[NEdge], i0, i1, nn;
		nn = 0;
		for (i = 0; i < 2*NEdge; i++) 
		    if (A[jj][i+NVert] > .5) /* face-edge */
			ii[nn++] = i+NVert;
		assert(nn <= 2); /* could be e---f---e, or e---f-x-e
				  * leave it alone */
		if (nn != 2)
		    continue;
		i0 = ii[0]; i1 = ii[1];
		A[i0][i1] = A[i1][i0] = 0.;
	    }
	}

	phgMatSetEntries(mat_cnnt, N, I, N, I, A[0]);
	
	for (i = 0; i < N; i++)
	    for (j = 0; j < N; j++)
		A[i][j] = 1.;
	phgMatAddEntries(mat_fill, N, I, N, I, A[0]);

	//vtkPause(verb);
	//vtkTmpActorsClear(verb);
    }
    vtkTmpActorsClear(verb);
    vtkPause(verb);

    phgMatAssemble(mat_cnnt);
    assert(mat_cnnt->type == PHG_UNPACKED);
    if (DUMP_MAT_VEC)
	phgMatDumpMATLAB(mat_cnnt, "line", "line_.m");
    phgMatAssemble(mat_fill);
    assert(mat_fill->type == PHG_UNPACKED);
    if (DUMP_MAT_VEC)
	phgMatDumpMATLAB(mat_fill, "fill", "fill_.m");

    /* 
     * step 2. Build line block:
     * */
    verb = 3;
    vtkTmpActorsBegin(verb);
    next_dof0 = 0;
    nd = nl = 0;
    nl_alloc = 8;
    nd_alloc = 8;
    line_dofs = phgAlloc(nd_alloc * sizeof(*line_dofs));
    line_ndof = phgAlloc(nl_alloc * sizeof(*line_ndof));
    line_dof0 = phgAlloc(nl_alloc * sizeof(*line_dof0));

    dof_iLine = phgAlloc(map->nlocal * sizeof(*dof_iLine));
    dof_iLngh = phgAlloc(map->nlocal * sizeof(*dof_iLine));
    for (i = 0; i < map->nlocal; i++) {
	dof_iLine[i] = dof_iLngh[i] = -1;
    }
    phgInfo(verb, "MAP: nlocal %d, local_size %d\n", 
	    map->nlocal, map->localsize);
    vef = phgDofSetupVEFMap2(g, u, VERT_FLAG | EDGE_FLAG | FACE_FLAG); /* only for DOF_P1 P2 CR1*/

    row = mat_cnnt->rows;
    for (i = 0; i < map->nlocal; i++, row++) {
	INT i0, i1, idf[IDFB_SIZE], idb[IDFB_SIZE]; /* dof idx forward & backward */
	MAT_ROW *row1;
	int nf, nb;
	int verb_cr = 3;

	/*
	 * skip dof_Ln that has been filled in line block.
	 * */
	if (dof_iLine[i] != -1)
	    continue;


	/* skip dof_Pt */
	if (row->ncols == 0) {
	    continue;
	}
	/* skip line that is cut to Pt
	 * TODO: cylic reduction */
	else if (row->ncols == 1) {	
	    if (row->cols[0] >= nlocal)
		continue;
	} 
	else if (row->ncols == 2) {	
	    if (row->cols[0] >= nlocal
		&& row->cols[1] >= nlocal)
		continue;
	}
	/* a dof can only connected to one or two other dof */
	else {
	    DEBUG_ROW_CONNECTION(row);
	}

	/* new line */
	phgInfo(verb, "\n=== === ===\nstart: %d\n", i);
	dof_iLine[i] = nl;
	vtkSetColor(verb, "red");
	//vtkTmpActorsBegin(verb);
	sprintf(vtk_tmp_str, "  %3d: Line %d", i, nl);
	//SHOW_V(crd + i*Dim, Dim);
	vtkDrawPoint(verb, crd + i*Dim, vtk_tmp_str);
	//vtkPause(verb+1);

	/* 1. forward link */
	nf = 0;
	if ((i1 = row->cols[0]) < nlocal) {
	    i0 = i; 
	    vtkSetColor(verb, "blue");
	    while (TRUE) {
		dof_iLine[i1] = nl;
		idf[nf++] = i1;
		{
		    sprintf(vtk_tmp_str, "  F%3d", i1);
		    vtkDrawPoint(verb, crd + i1*Dim, vtk_tmp_str);
		    vtkDrawLine(verb, crd + i1*Dim, crd + i0*Dim);
		    //vtkPause(verb+1);
		} /* end of vtk debug */
		
		/* set next dof */
		row1 = mat_cnnt->rows + i1;

		/* end */
		if (row1->ncols == 1) {
		    {
			assert(row1->cols[0] == i0); /* end pt must connect to prev pt */
			sprintf(vtk_tmp_str, "      End");
			vtkDrawPoint(verb, crd + i1*Dim, vtk_tmp_str);
			//vtkPause(verb+1);
		    } /* end of vtk debug */

		    break;
		}
		/* next */
		else if (row1->ncols == 2) {
		    if (row1->cols[0] == i0) {
			i0 = i1;
			i1 = row1->cols[1];
		    } else if (row1->cols[1] == i0) {
			i0 = i1;
			i1 = row1->cols[0];
		    } else {
			phgError(1, "line connection broken!\n");
		    }

		    /* off proc */
		    if (i1 >= nlocal) {
			assert(i1 < localsize);
			phgInfo(verb+1, "line cut by proc interface!\n");
			break;
		    }
		}
		else {
		    DEBUG_ROW_CONNECTION(row1);
		    phgError(1, "%d: line connection is not well defined. ncols: %d!\n",
			     __LINE__, row1->ncols);
		}
	    }
	    //SHOW_iV(idf, nf);
	}

	/* 2. backward link */
	nb = 0;
	if (row->ncols == 2 
	    && row->cols[1] < nlocal) {
	    i0 = i; 
	    i1 = row->cols[1];
	    assert(i1 < localsize);
	    vtkSetColor(verb, "yellow");
	    while (TRUE) {
		dof_iLine[i1] = nl;
		idb[nb++] = i1;
		{
		    sprintf(vtk_tmp_str, "  B%3d", i1);
		    vtkDrawPoint(verb, crd + i1*Dim, vtk_tmp_str);
		    vtkDrawLine(verb, crd + i1*Dim, crd + i0*Dim);
		    //vtkPause(verb+1);
		} /* end of vtk debug */

		
		/* set next dof */
		row1 = mat_cnnt->rows + i1;

		/* end */
		if (row1->ncols == 1) {
		    {
			assert(row1->cols[0] == i0);
			sprintf(vtk_tmp_str, "      End");
			vtkDrawPoint(verb, crd + i1*Dim, vtk_tmp_str);
			//vtkPause(verb+1);
		    } /* end of vtk debug */

		    break;
		}
		/* next */
		else if (row1->ncols == 2) {
		    if (row1->cols[0] == i0) {
			i0 = i1;
			i1 = row1->cols[1];
		    } else if (row1->cols[1] == i0) {
			i0 = i1;
			i1 = row1->cols[0];
		    } else {
			phgError(1, "line connection broken!\n");
		    }

		    /* off proc */
		    if (i1 >= nlocal) {
			assert(i1 < localsize);
			phgInfo(verb+1, "line cut by proc interface!\n");
			break;
		    }
		}
		else {
		    DEBUG_ROW_CONNECTION(row1);
		    phgError(1, "%d: line connection is not well defined. ncols: %d!\n",
			     __LINE__, row1->ncols);
		}
	    }
	    //SHOW_iV(idb, nb);
	}
	//vtkPause(verb);
	//vtkTmpActorsClear(verb);

	/*
	 * 3. record line:
	 *    line_dofs[nd; nd_alloc];
	 *    line_dof0[nl; nl_alloc];
	 *    line_ndof[nl; nl_alloc]; 
	 * */
	    
	if (nl >= nl_alloc) {
	    line_dof0 = phgRealloc_(line_dof0, (nl_alloc * 2) * sizeof(*line_dof0),
				    nl_alloc * sizeof(*line_dof0));
	    line_ndof = phgRealloc_(line_ndof, (nl_alloc * 2) * sizeof(*line_ndof),
				    nl_alloc * sizeof(*line_ndof));
	    nl_alloc *= 2;
	    phgInfo(verb+1, "* num line block alloc: %d\n", nl_alloc);
	}
	    
	ndof = nb + nf + 1;
	line_ndof[nl] = ndof;
	line_dof0[nl] = next_dof0;
	assert (ndof > 0);

	/* copy line dof indicies */
	REALLOC_VEC(line_dofs, nd, nd_alloc, ndof);
	dof0 = line_dofs + next_dof0;
	for (j = 0; j < nb; j++)
	    *(dof0++) = idb[nb-j-1];
	*(dof0++) = i;
	for (j = 0; j < nf; j++)
	    *(dof0++) = idf[j];

	dof0 = line_dofs + next_dof0; /* go to begin of dofs */
	SHOW_iV_(verb_cr, dof0, ndof);
	    
	dof0 = line_dofs + next_dof0; /* go to begin of dofs */
	SHOW_iV_(verb_cr, dof0, ndof);

	next_dof0 += ndof;
	nd += ndof;
	nl++;
    } /* end of local dofs */

    vtkPause(verb);
    vtkTmpActorsClear(verb);
    nds[0] = nd;
    SHOW_iV_(verb, line_dof0, nl);
    SHOW_iV_(verb, line_ndof, nl);
    SHOW_iV_(verb, line_dofs, nd);
    if (nl == 0)
	phgWarning("!!! Line Dofs not found !!!\n");

    /*
     * step 3. Set Dofs near line block.
     * 
     * */

    /* set DofP which is neighbour of line */
    verb = 3;
    nd = 0;
    ForAllElements(g, e) {
	int N = utype->nbas;
	INT I[N], ngh[N], min_iLine = -1;

	vtkSetColor(verb, "yellow");
	vtkSetTransparent(verb, 0.2);
	vtkDrawElement(verb, e);
	vtkSetTransparent(verb, 1.);

	//vtkTmpActorsBegin(verb);
	bzero(ngh, sizeof(ngh));
	for (i = 0; i < N; i++) {
	    j = phgMapE2L(map, 0, e, i);
	    j = phgMapL2V(map, j); 
	    I[i] = j;
	    
	    if (j >= nlocal) 
		continue;

	    /* Case 1: local Dof in line */
	    if ((l = dof_iLine[j]) != -1) {
		if (min_iLine != -1)
		    min_iLine = MIN(min_iLine, l);
		else
		    min_iLine = l;
	    }
	    /* Case 2: local Dof near line */
	    else
		ngh[i] = TRUE;
	}

	/* Case 3: no local dof is near line*/
	if (min_iLine == -1) {
	    phgInfo(verb, "Not dofL in elem.\n");

	    vtkTmpActorsBegin(verb);
	    vtkSetColor(verb, "yellow");
	    vtkDrawElement(verb, e);
	    vtkPause(verb);
	    vtkTmpActorsClear(verb);
	    continue;
	}

	/* Case 2:
	 * local Dof near line is set to related
	 * to the line of minimum index. */
	for (i = 0; i < N; i++) {
	    if (!ngh[i])
		continue;
	    if ((l = dof_iLngh[I[i]]) == -1) {
		dof_iLngh[I[i]] = min_iLine;
		nd++;
	    } else {
		dof_iLngh[I[i]] = MIN(l, min_iLine);
	    }
	}

	/* vtkSetColor(verb, "green"); */
	/* sprintf(vtk_tmp_str, "  %3d", I[i]); */
	/* vtkDrawPoint(verb, crd + I[i]*Dim, vtk_tmp_str); */
	/* vtkPause(verb); */

	//vtkPause(verb);
	//vtkTmpActorsClear(verb);
    }
    phgInfo(verb, "* num of Dof neigh: %d\n", nd);

    /*
     *  record dof line neighbour:
     *    lngh_dofs[nd; nd_alloc]
     *    lngh_dof0[nl; nl_alloc]
     *    lngh_ndof[nl; nl_alloc]
     *  */
    lngh_dofs = phgCalloc(nd, sizeof(*lngh_dofs));
    lngh_ndof = phgCalloc(nl, sizeof(*lngh_ndof));
    lngh_dof0 = phgCalloc(nl, sizeof(*lngh_dof0));
    for (i = 0; i < nlocal; i++) {
	if ((l = dof_iLngh[i]) == -1)
	    continue;
	lngh_ndof[l]++;
    }
    SHOW_iV_(verb, lngh_ndof, nl);
	
    next_dof0 = 0;
    for (i = 0; i < nl; i++) {
	lngh_dof0[i] = next_dof0;
	next_dof0 += lngh_ndof[i];
    }
    //assert(next_dof0 == nd);
    SHOW_iV_(verb, lngh_dof0, nl);
	
    ndofs = phgCalloc(nl, sizeof(*ndofs));
    for (i = 0; i < nlocal; i++) {
	if ((l = dof_iLngh[i]) == -1)
	    continue;
	lngh_dofs[lngh_dof0[l] + ndofs[l]++] = i;
    }
    phgFree(ndofs);
    nds[1] = nd;
    SHOW_iV_(verb, lngh_dofs, nd);
	
#if 1
    /* step 3'. Sort line neighbours */
    {
	int ib, idn_o[1000], tmp[1000];

	for (ib = 0; ib < nl; ib++) {
	    INT ii, *idx = lngh_dofs + lngh_dof0[ib];
	    int N = lngh_ndof[ib];
	    
	    //memcpy(x0, crd + idl[0]*Dim, Dim*sizeof(FLOAT));
	    for (i = 0; i < N; i++) {
		ii = idx[i];
		idn_o[i] = i;
		tmp[i] = ii;
		memcpy(crdx[i], crd + ii*Dim, Dim*sizeof(FLOAT));
	    }
	    qsort(idn_o, N, sizeof(*idn_o), coord_comp);
	    //SHOW_iV(idn_o, N);
	    for (i = 0; i < N; i++) {
		idx[i] = tmp[idn_o[i]];
	    }
	    //SHOW_iV(idx, N);
	}
    }
#endif


    /*
     * step 3''. Check line-neigh relationship
     *  */
    verb = CHECK_LINE_NEIGH;
    {
	/* V2D map */
	GTYPE dof_type[nlocal];
	int dof_idx[nlocal];

	ForAllElements(g, e) {
	    int N = utype->nbas;
	    GTYPE gtype;
	    int gidx, bidx;

	    for (i = 0; i < N; i++) {
		phgDofGetElementBasisInfo(u, e, i,
					  &gtype, &gidx, &bidx);
		j = phgMapE2L(map, 0, e, i);
		j = phgMapL2V(map, j); 
		    
		if (j >= nlocal)
		    continue;
		dof_type[j] = gtype;
		    
		if (gtype == VERTEX)
		    dof_idx[j] = e->verts[gidx];
		else if (gtype == EDGE)
		    dof_idx[j] = e->edges[gidx];
		else if (gtype == FACE)
		    dof_idx[j] = e->faces[gidx];
		else
		    dof_idx[j] = -1;
	    }
	}
		
	for (i = 0; i < nl; i++) {
	    //#warning Only plot short line
	    /* if (line_ndof[i] > 6) */
	    /* 	continue; */

	    vtkTmpActorsBegin(verb);

	    /* plot line */
	    ndof = line_ndof[i];
	    dof0 = line_dofs + line_dof0[i];
	    for (j = 0; j < ndof; j++) {
		int id, iv = dof0[j], nelem;
		
		vtkSetColor(verb, "red");
		vtkSetTransparent(verb, 1.);
		vtkDrawPoint(verb, crd + iv*Dim, NULL);

		if (j > 0)
		    vtkDrawLine(verb, crd + dof0[j-1]*Dim, 
				crd + dof0[j]*Dim);

		if (dof_type[iv] == VERTEX) {
		    id = dof_idx[iv];
		    nelem = vef->Vsize[id];
		    //phgInfo(verb, "nelem: %d\n", nelem);
		    for (k = 0; k < nelem; k++) {
			vtkSetColor(verb, "yellow");
			vtkSetTransparent(verb, 0.3);
			vtkDrawElement(verb, vef->Vmap[id][k]);

			/* vtkTmpActorsBegin(verb); */
			/* vtkSetTransparent(verb, 1.); */
			/* vtkSetColor(verb, "blue"); */
			/* vtkDrawElement(verb, vef->Vmap[id][k]); */
			/* vtkPause(verb); */
			/* vtkTmpActorsClear(verb); */
		    }
		} else if (dof_type[iv] == EDGE) {
		    id = dof_idx[iv];
		    nelem = vef->Esize[id];
		    //phgInfo(verb, "nelem: %d\n", nelem);
		    for (k = 0; k < nelem; k++) {
			vtkSetColor(verb, "yellow");
			vtkSetTransparent(verb, 0.1);
			vtkDrawElement(verb, vef->Emap[id][k]);
			//vtkPause(verb);
		    }
		} else if (dof_type[iv] == FACE) {
		    id = dof_idx[iv];
		    nelem = vef->Fsize[id];
		    //phgInfo(verb, "nelem: %d\n", nelem);
		    for (k = 0; k < nelem; k++) {
			vtkSetColor(verb, "yellow");
			vtkSetTransparent(verb, 0.1);
			vtkDrawElement(verb, vef->Fmap[id][k]);
			//vtkPause(verb);
		    }
		}
	    }
	    vtkPause(verb);

	    /* plot neigh */
	    ndof = lngh_ndof[i];
	    dof0 = lngh_dofs + lngh_dof0[i];
	    vtkSetColor(verb, "cyan");
	    vtkSetTransparent(verb, 1.);
	    for (j = 0; j < ndof; j++) {
		vtkDrawPoint(verb, crd + dof0[j]*Dim, NULL);
		//vtkPause(verb);
	    }

	    vtkPause(verb);
	    vtkTmpActorsClear(verb);
	}

    } /* end of check neigh */
    phgDofFreeVEFMap2(&vef);

    /*
     * step 3'''. Line length statistic: 
     * 1. average length
     * 2. num of short lines
     *  */
#define LENN 1000
    {
	int len[LENN];
	phgInfo(1, "*****************************\n");
	phgInfo(1, "*** Line length statistic ***\n");
	phgInfo(1, "*****************************\n");
	phgInfo(1, "total line num: %d\n", nl);

	/* line */
	phgInfo(verb, "averge line length: %E\n", 
		nds[0] / (FLOAT) (nl + 1e-5)); /* avoid divided by zero */
	bzero(len, sizeof(len));
	for (i = 0; i < nl; i++) {
	    ndof = line_ndof[i];
	    if (ndof < LENN-1)
		len[ndof]++;
	    else
		len[LENN-1]++;
	}
	for (k = 0; k < LENN-1; k++)
	    if (len[k] > 0)
		phgInfo(1, "num of line len %2d: %5d\n", 
			k, len[k]);
	phgInfo(1, "num of line len >= %2d: %5d\n", 
		LENN-1, len[LENN-1]);

	/* line neighs */
	phgInfo(1, "averge lngh length: %E\n", 
		nds[1] / (FLOAT) (nl + 1e-5));
	bzero(len, sizeof(len));
	for (i = 0; i < nl; i++) {
	    ndof = lngh_ndof[i];
	    if (ndof < LENN-1)
		len[ndof]++;
	    else
		len[LENN-1]++;
	}
	for (k = 0; k < LENN-1; k++)
	    if (len[k] > 0)
		phgInfo(1, "num of lngh len %2d: %5d\n", 
			k, len[k]);
	phgInfo(1, "num of lngh len >= %2d: %5d\n", 
		LENN-1, len[LENN-1]);

    }
#undef LENNN

    /*
     * Step 4. reorder line neigh dofs.
     *
     *  */
    {
	


    }

    /*
     * Step 5. record left dofs
     *
     *  */

    /* check for special case: no left dofs */
    //assert(nlocal == nds[0] + nds[1]);

    {
	INT *pts;

	nd = 0;
	for (i = 0; i < nlocal; i++)
	    if (dof_iLine[i] == -1
		&& dof_iLngh[i] == -1)
		nd++;
	nds[2] = nd;
	ndof_pt = nd;
	dof_pts = phgCalloc(nd, sizeof(*dof_pts));
	pts = dof_pts;
	
	for (i = 0; i < nlocal; i++)
	    if (dof_iLine[i] == -1
		&& dof_iLngh[i] == -1)
		*(pts++) = i;
    }

    phgInfo(verb, "# Line dofs: %d\n", nds[0]);
    phgInfo(verb, "# Lngh dofs: %d\n", nds[1]);
    phgInfo(verb, "# Left dofs: %d\n", nds[2]);

    /* save line block info */
    /* record in ML */
    bk = phgCalloc(1, sizeof(*bk));
    ml->block_dofs_[islv] = bk;
    memcpy(bk->nds, nds, sizeof(nds));
    bk->mat_cnnt = mat_cnnt;
    bk->map_cnnt = map;
    bk->nline = nl;	
    bk->line_dofs = line_dofs;	
    bk->line_dof0 = line_dof0;	
    bk->line_ndof = line_ndof;	
    bk->lngh_dofs = lngh_dofs; 	
    bk->lngh_dof0 = lngh_dof0;	
    bk->lngh_ndof = lngh_ndof;	
    bk->ndof_pt = ndof_pt;
    bk->dof_pts = dof_pts;
    bk->factorized = FALSE;	

    phgDofFree(&u);
    phgDofFree(&Coord);
    return;
}

/*
 * Debug routines.
 * Set mat entries_{i, j} = i.j; for local indicies i, j;
 *         offp indicies is led with 9.
 * Only for packed mat.
 *  */
#define IDX_BAS 10000
void debugMatrix(MAT *A)
{
    INT i, j, n, *pc, *pc_offp;
    FLOAT *pd, *pd_offp;

    assert(A->type != PHG_DESTROYED);
    if (!A->assembled)
	phgMatAssemble(A);
    assert(A->type != PHG_MATRIX_FREE);
    phgMatPack(A);

    pc = A->packed_cols;
    pd = A->packed_data;

    if (A->cmap->nprocs > 1) {
	pc_offp = A->packed_cols + A->rmap->nlocal + A->nnz_d;
	pd_offp = A->packed_data + A->nnz_d;
    } else {
	pc_offp = NULL;
	pd_offp = NULL;
    }


    for (i = 0; i < A->rmap->nlocal;  i++) {
	INT jcol;

	/* local data */
	if ((n = *(pc++)) != 0) {
	    for (j = 0; j < n; j++) {
		jcol = *(pc++);
		*(pd++) = IDX_BAS + i + jcol / (10.*IDX_BAS) + 0.1;
	    }
	}
	    
	/* remote data */
	if (pc_offp != NULL && (n = *(pc_offp++)) != 0) {
	    for (j = 0; j < n; j++) {
		jcol = *(pc_offp++);
		*(pd_offp++) = IDX_BAS + i + jcol / (10.*IDX_BAS) + 0.9;
	    }
	}
    }

    //phgMatDumpMATLAB2(A, "Aij", "Aij_.m");
    return;
}


static int 
int_comp(const void *pa, const void *pb)
{
    int *a = (int *) pa, *b = (int *) pb;  
    return *a - *b;
}

#define MAX_DOF_NEIGHBOUR 1280
#define MAX_MAT_NCOL 1280

/* Save LU Factorization of each line block, and corresponding
 * pivot indices, use packed format storage.
 *
 * TODO: use pointer arithmatic.
 * */
void 
phgMultiGridFactorizeLineBlock(MG_LEVEL *ml)
{
    GRID *g = ml->grid;
    MAT *A = ml->mat;
    DOF_TYPE *utype = ml->u->type;
    int order = utype->order;
    MG_BLOCK_DOFS *bk = ml->block_dofs;
    //MAP *map = A->rmap;
    int dim = ml->u->dim;
    INT i, ib, k, *pc, *pc0, *pocs, nl = bk->nline;
    size_t *ps;
    INT *line_ndof, *line_dof0,  *line_dofs,
	*line_mat0, *line_piv0, *line_pivs,
	*line_poc0, *line_pocs, /* packed off diagnal block colume entries */
	*ipiv, idx_o[MAX_DOF_NEIGHBOUR][2], pc_o[MAX_MAT_NCOL][2];
    INT mat_alloc = 256, mat_size = 0, 
	piv_alloc = 16, piv_size = 0,
	poc_alloc = 16, poc_size = 0;
    BYTE ly[MAX_MAT_NCOL];
    FLOAT *pd, *pd0, *line_matv, *matv;
    int ngh_size, row_size, row_diag, nly_b;
    int verb = 3;
    
    Unused(g);
    assert(bk != NULL);
    if (bk->factorized)
	return;

    phgPrintf("\n         Factorize %s, ", ml->u->name);
    /* TODO: compress mat */
    /* { */
    /* 	int nlocal = map->nlocal,  */
    /* 	    nnz = A->nnz_d + A->nnz_o; */
    /* 	poc0 = phgAlloc((2*nlocal + nnz) * sizeof(*poc)); */
    /* 	bk->packed_offD_col = poc0; */
    /* 	phgInfo(verb, "# poc size %d\n", 2*nlocal + nnz); */
    /* } */

    line_dof0 = bk->line_dof0;
    line_ndof = bk->line_ndof;
    line_dofs = bk->line_dofs;

    line_mat0 = phgCalloc(nl, sizeof(*line_mat0));             /* line block mat start */
    line_matv = phgCalloc(mat_alloc, sizeof(*line_matv)); /* line block mat entries */
    line_piv0 = phgCalloc(nl, sizeof(*line_piv0));;            /* line block pivot start */
    line_pivs = phgCalloc(piv_alloc, sizeof(*line_pivs)); /* line block pivot indices */
    line_poc0 = phgCalloc(nl, sizeof(*line_poc0));		  /* block diag entries start */
    line_pocs = phgCalloc(poc_alloc, sizeof(*line_pocs)); /* block diag entries  */

    /* Neighbour size:
     *   Scale:   P1: 1, P2: 2, P3: 3
     *   Vector:  P1: 5, P2: 8, P3: 11
     *  */
    if (DOF_TYPE_LAGRANGE(utype))
	ngh_size = bk->neigh_size = (order + 1) * dim - 1;
    else if (DOF_TYPE_SUBELEM(utype))
	ngh_size = bk->neigh_size = (1 + 1) * dim - 1;
    else
	phgError(1, "Unsupported DOF type: %s!\n", utype->name);

    row_size = 2*ngh_size + ngh_size + 1; /* LDAB >= 2*KL+KU+1 */
    row_diag = 2*ngh_size;		  /* diagnol position of row */
    poc_size = 0;		/* # of diag block entries */
    assert(2*ngh_size + 1 < MAX_DOF_NEIGHBOUR);
    phgPrintf("neigh size: %d, ", ngh_size);

#if CHECK_MAT_ENTRIES
    debugMatrix(A);			   /* debug */
    phgMatDumpMATLAB2(A, "Aij", "Aij_.m"); /* debug */
#endif /* CHECK_MAT_ENTRIES */
    assert(A->type != PHG_DESTROYED);
    if (!A->assembled)
	phgMatAssemble(A);
    assert(A->type != PHG_MATRIX_FREE);
    /* Note: It's assumed that the matrix may be used in other solver,
     *       so packed format is always needed. */
#if CHECK_MAT_ENTRIES
    phgMatDumpMATLAB(A, "A", "A_.m"); /* debug */
#endif /* CHECK_MAT_ENTRIES */
    phgMatPack(A);

    pc0 = A->packed_cols;
    pd0 = A->packed_data;
    ps = A->packed_ind;
    //ly0 = bk->line_entry;
	
    verb = 3;
    for (ib = 0; ib < nl; ib++) {
	INT *idx0 = line_dofs + line_dof0[ib];
	int k0, N, N0 = line_ndof[ib]; 	/* # of dofs in line block */
	INT idx[N0*dim];

	/* --------------------------------
	 * Extend idx for vector dofs.
	 * Todo: faster code?
	 * -------------------------------- */
	for (i = 0; i < N0; i++)
	    for (k0 = 0; k0 < dim; k0++)
		idx[i*dim + k0] = idx0[i] * dim + k0;
	N = N0*dim;
	nly_b = 0;		/* # of off diag block entries in block_ib */
	
	phgInfo(verb, "\n### Line: %3d\n", ib);
	SHOW_iV_(verb, idx, N);

	/* realloc for mat and piv */
	REALLOC_VEC(line_matv, mat_size, mat_alloc, N * row_size);
	REALLOC_VEC(line_pivs, piv_size, piv_alloc, N);

	matv = line_matv + mat_size;
	ipiv = line_pivs + piv_size;
	line_mat0[ib] = mat_size;
	line_piv0[ib] = piv_size;
	line_poc0[ib] = poc_size;
	pocs = line_pocs + line_poc0[ib];
	bzero(matv, N * row_size * sizeof(*matv));
	bzero(ipiv, N * sizeof(*ipiv));
	/* band matrix storage use LAPACK style, LDAB >= 2*KL+KU+1,
	 * the extra KL row is used to store factorization. */
	
	/* fill band matrix row_i */
	for (i = 0; i < N; i++) {
	    int ii = idx[i];
	    int ncol, n, kmax, nly = 0; /* # of diag block entries in row_i */
	    INT j, j0, m = 0;
	    
	    phgInfo(verb, "\nline entry[%d]: %d\n", i, ii);
	    /* get neighbour dof index */
	    j = 0;
	    kmax = MIN(ngh_size, i);
	    for (k = 1; k <= kmax; k++, j++) {
		idx_o[j][0] = idx[i-(kmax-k)-1];
		idx_o[j][1] = j;
	    }
	    idx_o[j][0] = ii;
	    idx_o[j][1] = j;
	    j0 = j++; 		/* record mid row pos */
	    kmax = MIN(ngh_size, N-1-i);
	    for (k = 1; k <= kmax; k++, j++) {
		idx_o[j][0] = idx[i+k];
		idx_o[j][1] = j;
	    }
	    n = j;		/* # of neighbour dofs */
	    //assert(n == (((i == 0 || i == N-1) ? ngh_size + 1 : 2*ngh_size + 1))); 

	    /* sort neighbour dof index */
	    SHOW_iM_(verb, idx_o[0], n, 2);
	    qsort(idx_o, n, sizeof(*idx_o), int_comp);
	    SHOW_iM_(verb, idx_o[0], n, 2);

	    /* get sub matrix row */
	    pc = pc0 + PACK_COL(ps, ii);
	    pd = pd0 + PACK_DAT(ps, ii);
	    ncol = *(pc++);
	    assert(ncol != 0);
	    assert(ncol < MAX_MAT_NCOL);
	    //printf("[ncol %d] ", ncol);
	    
	    /* dirichlet entries */
	    if (ncol == 1) { 
		assert(*(pc++) == ii); 
		matv[(row_diag) + i*row_size] = *(pd++);

		nly = 1;
		/* realloc */
		REALLOC_VEC(line_pocs, poc_size, poc_alloc, nly+1);

		pocs = line_pocs + line_poc0[ib] + nly_b;
		*(pocs++) = 1;	/* only one entry, that's itself */
		*(pocs++) = 0;	
		nly_b += nly+1;
		poc_size += nly+1;	
		phgInfo(verb+1, "mat row[%5d]: %d\n", ii, 1);
		continue;
	    }

	    for (k = 0; k < ncol; k++) {
		pc_o[k][0] = pc[k];
		pc_o[k][1] = k;
	    }
	    SHOW_iM_(verb, pc_o[0], ncol, 2);
	    qsort(pc_o, ncol, sizeof(*pc_o), int_comp);
	    SHOW_iM_(verb, pc_o[0], ncol, 2);
		    

	    /* find mat entries for neighbour dofs*/
	    bzero(ly, ncol * sizeof(*ly));
	    j = 0;
	    //phgInfo(verb, "\n--- ---\n row: %d[%d]\n", ii, ncol);
	    for (k = 0; k < n; k++) {
		int ofs;
		while (j < ncol &&
		       (m = idx_o[k][0] - pc_o[j][0]) > 0) {
		    //phgInfo(verb, "pc[%d]_%d  ", j, pc_o[j]);
		    j++;
		}	
		if (j == ncol) 
		    break;
		
		/*  LAPACK band mat storage mapping:
		 *   AB(kl+ku+1+i-j,j) = A(i,j)
		 *   for max(1,j-ku)<=i<=min(m,j+kl)
		 *
		 *  e.g. when  M = N = 6, KL = 2, KU = 1
		 *
		 *      *    *    *    +    +    +    
		 *      *    *    +    +    +    +    
		 *      *   a12  a23  a34  a45  a56   
		 *     a11  a22  a33  a44  a55  a66   <---- mid
		 *     a21  a32  a43  a54  a65   *    
		 *     a31  a42  a53  a64   *    *    
		 *
		 *   */
#if 0
		assert(!m); /* Generally mat entries for neighbour should NOT be
			     * zero, however, it may happen in some case. */
#else
		if (0 & m) {
		    phgWarning("*** Dirich entries: row %d: %d\n", ii, idx_o[k][0]);
		    //SHOW_iM(idx_o[0], n, 2);
		    //SHOW_iM(pc_o[0], ncol, 2);
		    /* if mat handle bdry entries & col happens to be Dirich entries
		     * then A_{i, j} is zero. */
		    if (pc0[ PACK_COL(ps, idx_o[k][0])] != 1) {
			phgInfo(verb, "A[%d][%d] got zero entry ?\n", 
			       ii, idx_o[k][0]);
		    }
		} else {
		    if (!m) 	/* found */
			phgInfo(verb, "   %3d -> %3d[%d]\n", ii, idx_o[k][0], 
			       idx_o[k][1] - j0);
		}
#endif
		ofs = idx_o[k][1] - j0;
		assert(ofs <= ngh_size);
		matv[(row_diag - ofs) + (i+ofs)*row_size] = 
		    (m) ? 0. : pd[pc_o[j][1]];
		if (!m) { /* found */
		    phgInfo(verb, "# i: %d j: %d\n", i, ofs);
		    ly[pc_o[j][1]] = TRUE; /* block entries */
		    nly++;
		}
	    } /* end of neigh dofs of dof_i */
	    //phgInfo(verb, "\n");

	    /* Record block diag entries:
	     * the j-th cols of row_i is block diag entries or not.
	     *
	     * Note: this could be different even for same row!
	     *  Because the row may lie in different lines.
	     * 
	     * Either we record off block diag entries,
	     *  or we record block diag entries,
	     *  it depends on the num of off block diag entries in each row.
	     * TODO: dynamicly decide.
	     *
	     * */
	    assert(nly <= ncol);

	    /* realloc */
	    REALLOC_VEC(line_pocs, poc_size, poc_alloc, nly+1);

	    pocs = line_pocs + line_poc0[ib] + nly_b;
	    *(pocs++) = nly;
	    for (j = 0; j < ncol; j++) 
		if (ly[j]) {	/* Block diag entries */
		    *(pocs++) = j;
		    phgInfo(verb, "   i:%d, j:%d, %d\n", ii, j, pc[j]);
		}
	    nly_b += nly+1;
	    poc_size += nly+1;	
	    phgInfo(verb+1, "mat row[%5d]: %d\n", ii, ncol);
    	} /* end of block band matrix row_i */

	pocs = line_pocs + line_poc0[ib];
	SHOW_iV_(verb, pocs, nly_b);

	SHOW_M_(verb, matv, N, row_size);
	vtkPause(verb);
	/* LAPACK: triangular factorization */
	{
	    int M = N, KL = ngh_size, KU = ngh_size, LDAB = row_size, 
		INFO = 0, *IPIV = ipiv;
	    FLOAT *AB = matv;
	    dgbtrf_(&M, &N, &KL, &KU, AB, &LDAB, IPIV, &INFO);
	    if (INFO)
		phgError(1, "LAPACK DGBTRF error!!! code: %d.\n", INFO);
	}
	SHOW_M_(verb, matv, N, row_size);

	mat_size += N * row_size;
	piv_size += N;
    } /* end of line block */
    SHOW_V_(verb, line_matv, mat_size);
    SHOW_iV_(verb, line_pivs, piv_size);

    bk->line_mat0 = line_mat0;
    bk->line_matv = line_matv;
    bk->line_piv0 = line_piv0;
    bk->line_pivs = line_pivs;
    bk->line_poc0 = line_poc0;
    bk->line_pocs = line_pocs;
    bk->factorized = TRUE;
}



void 
phgMultiGridFreeFactLineBlock(MG_LEVEL *ml) 
{
    MG_BLOCK_DOFS *bk = ml->block_dofs;

#define FREE_SET(pt) phgFree(pt); pt = NULL;
    FREE_SET(bk->line_matv);
    FREE_SET(bk->line_piv0);
    FREE_SET(bk->line_pivs);
    FREE_SET(bk->line_poc0);
    FREE_SET(bk->line_pocs);
    bk->factorized = FALSE;
#undef FREE_SET
}



void 
mg_GaussSidel_line(MAT *A, VEC *x, VEC *b, int nsmooth, void *ctx)
{
    MG_LEVEL *ml = (MG_LEVEL *) ctx;
    MG_BLOCK_DOFS *bk = ml->block_dofs;
    //MAP *map= ml->map;
    INT i, j, k, k_smooth, ib, n, jcol, *pc, *pc0, *poc, *pc_offp, nlocal;
    size_t *ps;
    int dim = ml->u->dim, ngh_size, row_size;
    FLOAT *pd, *pd0, *pd_offp, *vx, *vb;
    FLOAT sum, omega = _p->smooth_damp;
    BOOLEAN remote_data = FALSE, factorized;
    INT *line_ndof, *line_dof0, *line_dofs,
	*lngh_ndof, *lngh_dof0, *lngh_dofs,
	*line_pivs, *line_mat0,	*line_piv0,
	*line_pocs, *line_poc0, *dof_pts;
    FLOAT *line_matv;
    int verb = 3;
#if USE_MPI
    FLOAT *offp_data = NULL;
#endif	/* USE_MPI */
#if TEST_LINE_SMOOTHING
    INT *smoothed = phgCalloc(x->map->nlocal, sizeof(*smoothed));
#endif	/* TEST_LINE_SMOOTHING */
    
#if TEST_LINE_SMOOTHING
    phgPrintf("*** TEST Line smoothing: Level[%d].\n", ml->level);
#endif	/* TEST_LINE_SMOOTHING */

    //phgVecDumpMATLAB(x, "xx", "xx_.m");
    MagicCheck(VEC, x);
    MagicCheck(VEC, b);
    assert(x == NULL || x->nvec == 1);
    assert(b == NULL || b->nvec == 1);
    if (x != NULL && !x->assembled)
	phgVecAssemble(x);
    if (b != NULL && !b->assembled)
	phgVecAssemble(b);
    //phgVecDumpMATLAB(x, "xxx", "xxx_.m");
    
    assert(A->type != PHG_DESTROYED);
    if (!A->assembled)
	phgMatAssemble(A);
    assert(A->type != PHG_MATRIX_FREE);
    phgMatPack(A);

#if USE_MPI
    if (A->cmap->nprocs > 1) {
	offp_data = phgAlloc(A->cinfo->rsize * sizeof(*offp_data));
    }
#endif	/* USE_MPI */

    if (A->cmap->nlocal != A->rmap->nlocal ||
	A->cmap->nlocal != x->map->nlocal  ||
	A->cmap->nlocal != b->map->nlocal)
	phgError(1, "%s:%d: inconsistent matrix-vector.", __FILE__,
		 __LINE__);
    nlocal = A->cmap->nlocal;

#if USE_MPI
    if (A->cmap->nprocs > 1) {
	phgMapScatterBegin(A->cinfo, x->nvec, x->data, offp_data);
	phgMapScatterEnd(A->cinfo, x->nvec, x->data, offp_data);
    }
#endif	/* USE_MPI */

    if (A->cmap->nprocs > 1)
	remote_data = TRUE;

    assert(bk != NULL);
    /* Dofs line*/
    line_ndof = bk->line_ndof;
    line_dof0 = bk->line_dof0;
    line_dofs = bk->line_dofs;
    /* Dofs line neighs */
    lngh_ndof = bk->lngh_ndof;
    lngh_dof0 = bk->lngh_dof0;
    lngh_dofs = bk->lngh_dofs;
    /* factorized matrix */
    line_mat0 = bk->line_mat0;
    line_matv = bk->line_matv;
    line_piv0 = bk->line_piv0;
    line_pivs = bk->line_pivs;
    line_poc0 = bk->line_poc0;
    line_pocs = bk->line_pocs;

    factorized = bk->factorized;
    assert(factorized);
    ngh_size = bk->neigh_size;
    row_size = 2*ngh_size + ngh_size + 1; /* LDAB >= 2*KL+KU+1 */

    vx = x->data;
    vb = b->data;
    pc0 = A->packed_cols;
    pd0 = A->packed_data;
    ps = A->packed_ind;;

    verb = 3;
    /* iteration */
    for (k_smooth = 0; k_smooth < nsmooth; k_smooth++) {

	/*
	 * line block:
	 *
	 * line smooth operation:
	 *   1. take line dof entries of mat as local left hand mat, 
	 *      the factorazition has already been done in
	 *      phgMultiGridFactorizeLineBlock.
	 *   2. multiply rest dof entries of mat to local right hand  
	 *   3. solver the local system
	 *
	 * use LAPACK routines to solve, the factorization of each line
	 * block has already been done.
	 * */
	for (ib = 0; ib < bk->nline; ib++) {
	    {
		INT *idx0 = line_dofs + line_dof0[ib]; 
		int ncols, k0, N, N0 = line_ndof[ib];		       /* # of dof in line */
		INT idx[N0*dim];
		FLOAT *matv = line_matv + line_mat0[ib];
		INT *ipiv = line_pivs + line_piv0[ib];
		static FLOAT *rhs = NULL;              /* rhs of line block */
		static int rhs_size = 0;               /* max rhs size */

		/* extend idx for vector Dof */
		for (i = 0; i < N0; i++)
		    for (k0 = 0; k0 < dim; k0++) {
			idx[i*dim + k0] = idx0[i] * dim + k0;
		    }
		N = N0*dim;

		if (N > rhs_size) {
		    rhs_size = N;
		    rhs = phgRealloc0(rhs, rhs_size * sizeof(*rhs));
		    phgFreeAtExit((void **)(void *)&rhs);
		}

		phgInfo(verb, "\n\n### check line: %d\n", ib);
		for (i = 0; i < N; i++) {
		    int ii = idx[i];
		    rhs[i] = vb[ii];
		}
		phgInfo(verb, "   init rhs\n");
		SHOW_V_(verb, rhs, N);
		
		poc = line_pocs + line_poc0[ib];

		/* local data */
		for (i = 0; i < N; i++) { /* for each dof in line */
		    int ii = idx[i], nly;
		    static FLOAT *pbc = NULL;              /* block columes */
		    static int pbc_size = 0;               /*  */

		    pc = pc0 + PACK_COL(ps, ii);
		    pd = pd0 + PACK_DAT(ps, ii);
		    ncols = *(pc++);
		    nly = *(poc++);
		    assert(ncols != 0);

		    if (ncols > pbc_size) {
			pbc_size = ncols + 1000;
			pbc = phgRealloc0(pbc, pbc_size * sizeof(*pbc));
			phgInfo(verb, "pbc realloc: %d\n", pbc_size);
			phgFreeAtExit((void **)(void *)&pbc);
		    }

		    bzero(pbc, ncols * sizeof(*pbc));
		    phgInfo(verb, "row: %d, nly: %d\n", ii, nly);
		    for (j = 0; j < nly; j++, poc++) {
			pbc[*poc] = 1;
			phgInfo(verb, "   poc: %d\n", *poc);
		    }
		    for (j = 0; j < ncols; j++, pd++, pc++) {
			if (pbc[j] == 1)
			    continue;
			jcol = *pc;
			phgInfo(verb, "rhs[%3d] %14.6e -= a[%3d] %14.6e * x[%3d] %14.6e", 
				ii, rhs[i], jcol, *pd, jcol, vx[jcol]);
			rhs[i] -= *pd * vx[jcol];
			phgInfo(verb, " = %14.6e\n", rhs[i]);
		    }
		}
		
		/* remote data */
		if (remote_data) {
		    for (i = 0; i < N; i++) {
			int ii = idx[i];
			pc_offp = pc0 + PACK_COL_OFFP(ps, ii, nlocal);
			pd_offp = pd0 + PACK_DAT_OFFP(ps, ii, nlocal);
			if ((n = *(pc_offp++)) != 0) {
			    for (j = 0; j < n; j++) {
				jcol = *(pc_offp++);
				rhs[i] -= *(pd_offp++) * offp_data[jcol];
			    }
			}
		    }
		}
		
		phgInfo(verb, "   before solve\n");
		SHOW_V_(verb, rhs, N);
		SHOW_M_(verb, matv, N, row_size);
		SHOW_iV_(verb, ipiv, N);
		/* solve */
		{
		    int KL = ngh_size, KU = ngh_size, LDAB = row_size, NRHS = 1,
			INFO = 0, LDB = N, *IPIV = ipiv;
		    FLOAT *AB = matv, *B = rhs;
		    char TRANS = 'N';
		    dgbtrs_(&TRANS, &N, &KL, &KU, &NRHS, AB, &LDAB, IPIV, 
			    B, &LDB, &INFO);
		    if (INFO)
			phgError(1, "LAPACK DGBTRS error!!! code: %d.\n", INFO);
		}
		phgInfo(verb, "   after solve\n");
		SHOW_V_(verb, rhs, N);	/* when Ax=b, dx = 0 */

		/* check smoothing */ 
		{
		    static FLOAT vvx[10000];
		    for (i = 0; i < N; i++) {
			vvx[i] = vx[idx[i]];          /* dx */
			phgInfo(verb+1, "  [%5d]: %12.5E ==> %12.5E\n", 
				idx[i], vx[idx[i]], vvx[i]);//vx[idx[i]]);
		    }
		    SHOW_V_(verb, vvx, N);
		}

		for (i = 0; i < N; i++) {
		    rhs[i] -= vx[idx[i]];          /* dx */
		    vx[idx[i]] += omega * rhs[i];
		}
		phgInfo(verb, "   should be zero\n");
		SHOW_V_(verb, rhs, N);	/* when Ax=b, dx = 0 */
#if TEST_LINE_SMOOTHING
		for (i = 0; i < N; i++) {
		    assert(fabs(rhs[i]) < 1e-12);
		    smoothed[idx[i]]++;
		}
#endif	/* TEST_LINE_SMOOTHING */
	    } /* end of one line block */

	    /*
	     * Line neighbour:
	     * smooth seperately as dofP.
	     *
	     * TODO: better smooth order ???
	     * */
	    //#warning Do not smooth neigh
	    if (1) {
		INT *idx0 = lngh_dofs + lngh_dof0[ib]; 
		int k0, N, N0 = lngh_ndof[ib];		       /* # of dof in line */
		INT idx[N0*dim];
		FLOAT aa = 0., dx;

		/* extend idx for vector Dof */
		for (i = 0; i < N0; i++)
		    for (k0 = 0; k0 < dim; k0++) {
			idx[i*dim + k0] = idx0[i] * dim + k0;
		    }
		N = N0*dim;

		for (i = 0; i < N; i++) {
		    int ii = idx[i];
		
		    sum = vb[ii];
		    /* local data */
		    pc = pc0 + PACK_COL(ps, ii);
		    pd = pd0 + PACK_DAT(ps, ii);
		    if ((n = *(pc++)) != 0) {
			for (j = 0; j < n; j++) {
			    jcol = *(pc++);
			    if (jcol != ii) {
				sum -= *(pd++) * vx[jcol];
			    } else {
				aa = *(pd++);
			    }
			}
		    }

		    /* remote data */
		    if (remote_data) {
			pc_offp = pc0 + PACK_COL_OFFP(ps, ii, nlocal);
			pd_offp = pd0 + PACK_DAT_OFFP(ps, ii, nlocal);
			if ((n = *(pc_offp++)) != 0) {
			    for (j = 0; j < n; j++) {
				jcol = *(pc_offp++);
				sum -= *(pd_offp++) * offp_data[jcol];
			    }
			}
		    }

		    dx = sum / aa - vx[ii];
#if TEST_LINE_SMOOTHING
		    assert(fabs(dx) < 1e-12);
		    smoothed[ii]++;
#endif	/* TEST_LINE_SMOOTHING */
		    vx[ii] += omega * dx;
		} /* end of line neigh_i */
	    }	  /* end of all line neighs */
	}	  /* end of all line blocks */


	/* Dof_pt */
	dof_pts = bk->dof_pts;
	/* 
	 * scalar smooth
	 * */
	if (dim == 1) {
	    for (i = 0; i < bk->ndof_pt; i++, dof_pts++) {
		INT ii = *dof_pts;
		FLOAT aa = 0., dx;

		sum = vb[ii];
		/* local data */
		pc = pc0 + PACK_COL(ps, ii);
		pd = pd0 + PACK_DAT(ps, ii);
		if ((n = *(pc++)) != 0) {
		    for (j = 0; j < n; j++) {
			jcol = *(pc++);
			if (jcol != ii) {
			    sum -= *(pd++) * vx[jcol];
			} else {
			    aa = *(pd++);
			}
		    }
		}

		/* remote data */
		if (remote_data) {
		    pc_offp = pc0 + PACK_COL_OFFP(ps, ii, nlocal);
		    pd_offp = pd0 + PACK_DAT_OFFP(ps, ii, nlocal);
		    if ((n = *(pc_offp++)) != 0) {
			for (j = 0; j < n; j++) {
			    jcol = *(pc_offp++);
			    sum -= *(pd_offp++) * offp_data[jcol];
			}
		    }
		}

		dx = sum / aa - vx[ii];
#if TEST_LINE_SMOOTHING
		assert(fabs(dx) < 1e-12);
		smoothed[ii]++;
#endif	/* TEST_LINE_SMOOTHING */
		vx[ii] += omega * dx;
	    } 
	}
	/*
	 * vec smoother, like GSv
	 * */
	else {
	    assert(dim == Dim);
	    for (i = 0; i < bk->ndof_pt; i++, dof_pts++) {
		INT ii = (*dof_pts) * dim;
		INT jcol;
		FLOAT sum[Dim], aa[Dim][Dim], det, dx[Dim];

		/* TODO: vec smoother for 3 componets */
		phgError(-1, "Vec smoother not implemented!!!\n");

		memset(aa, 0, sizeof(aa));
		sum[0] = vb[ii  ];
		sum[1] = vb[ii+1];
		sum[2] = vb[ii+2];

		/* local data */
		pc = pc0 + PACK_COL(ps, ii);
		pd = pd0 + PACK_DAT(ps, ii);
		for (k = 0; k < Dim; k++) {
		    if ((n = *(pc++)) != 0) {
			for (j = 0; j < n; j++) {
			    jcol = *(pc++);
			    if (jcol < ii || jcol > ii+2) {
				sum[k] -= *(pd++) * vx[jcol];
			    } else { /* offD, jcol = ii,i+1,ii+2 */
				aa[k][jcol - ii] = *(pd++);
			    }
			}
		    }
		}

		/* remote data */
		if (A->cmap->nprocs > 1) {
		    pc_offp = pc0 + PACK_COL_OFFP(ps, ii, nlocal);
		    pd_offp = pd0 + PACK_DAT_OFFP(ps, ii, nlocal);
		    for (k = 0; k < Dim; k++) {
			if ((n = *(pc_offp++)) != 0) {
			    for (j = 0; j < n; j++) {
				jcol = *(pc_offp++);
				sum[k] -= *(pd_offp++) * offp_data[jcol];
			    }
			}
		    }
		}

		/* solve */
		det = (aa[0][0] * aa[1][1] - aa[0][1] * aa[1][0]);
		assert(Fabs(det) > 1e-12);
		det = 1. / det;
		dx[0] = (aa[1][1] * sum[0] - aa[0][1] * sum[1]) * det - vx[ii  ];
		dx[1] = (aa[0][0] * sum[1] - aa[1][0] * sum[1]) * det - vx[ii+1];
		dx[2] = (1./aa[2][2] * sum[2]) - vx[ii+2];
#if TEST_LINE_SMOOTHING
		assert(fabs(dx[0]) + fabs(dx[1]) + fabs(dx[2]) < 1e-12);
		smoothed[ii  ]++;
		smoothed[ii+1]++;
		smoothed[ii+2]++;
#endif	/* TEST_LINE_SMOOTHING */
		vx[ii  ] += omega * dx[0];
		vx[ii+1] += omega * dx[1];
		vx[ii+2] += omega * dx[2];
	    } /* end of pt dof */
	}     /* end of dof_pt */


#if TEST_LINE_SMOOTHING
	/* Each entry is smoothed once and only once.
	 * For Box smoother, smoothing many happen several times.
	 *
	 * */
	for (i = 0; i < x->map->nlocal; i++)
	    assert(smoothed[i] == 1);
	phgPrintf("*** TEST Line smoothing: passed for Level[%d].\n", ml->level);
#endif	/* TEST_LINE_SMOOTHING */

#if USE_MPI
	if (A->cmap->nprocs > 1) {
	    phgMapScatterBegin(A->cinfo, x->nvec, x->data, offp_data);
	    phgMapScatterEnd(A->cinfo, x->nvec, x->data, offp_data);
	}
#endif	/* USE_MPI */
    }


#if USE_MPI
    phgFree(offp_data);
#endif	/* USE_MPI */

    return;
}
