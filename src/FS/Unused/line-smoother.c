#include "layers.h"
#include "vtk-draw.h"
#include <string.h>
#include <math.h>


/*
 * Line smoother for thermal solver
 *
 * */

#define CHECK_MAT_ENTRIES 0
#define DUMP_MAT_VEC 0


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
#define REALLOC_VEC(vecv, vec_size, vec_alloc, add_alloc) {	\
	INT vec_alloc##0 = vec_alloc;				\
	Unused(vec_alloc##0);					\
	while ((vec_size) + (add_alloc) > vec_alloc) {		\
	    vec_alloc *= 2;					\
	    vecv = phgRealloc_(vecv, vec_alloc * sizeof(*vecv),	\
			       vec_alloc##0 * sizeof(*vecv));	\
	}							\
    }

#define PACK_COL(ps, i) (ps[i])
#define PACK_DAT(ps, i) (ps[i] - i)
#define PACK_COL_OFFP(ps, i, nlocal) (ps[i+nlocal])
#define PACK_DAT_OFFP(ps, i, nlocal) (ps[i+nlocal] - i - nlocal)



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



MG_BLOCK_DOFS *
init_line_block(DOF *u0, LAYERED_MESH *gL)
{
    GRID *g = u0->g;
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

    int nlayer = gL->max_nlayer;
    int ndof_layer = 0;


    Unused(order);
    phgPrintf("   Init Line for dof: %s\n", u0->name);


    if (u0->type == DOF_P1)
	ndof_layer = nlayer + 1;
    else if (u0->type == DOF_P2)
	ndof_layer = 2 * nlayer + 1;
    else
	phgError(1, "Unsupported DOF type: %s\n", 
		 u0->type);

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
	    phgDofSetElementData(Coord, e, value[0]);
	}
	//phgDofDump(Coord);

#if 0
	/* dump coord to matlab */
	crd_map = phgMapCreate(Coord, NULL); 
	crd_vec = phgMapCreateVec(crd_map, 1);
	phgMapDofToLocalData(crd_map, 1, &Coord, crd_vec->data);
	crd = crd_vec->data;
	/* phgVecDumpMATLAB(crd_vec, "coord", "coord_.m"); */
	/* phgVecDestroy(&cvec); */
	/* phgMapDestroy(&cmap); */
#endif
    }	/* end of vtk debug */


    /*
     * step 1. Build connection matrix:
     * connection relationship is determined by dof direciton.
     * 
     * */
    verb = 3;
    phgInfo(verb, "Use direction to build connection matrix.\n");
 
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


	/*
	 * DOF P2:
	 *   reduce connection of line on edge to v-e-v
	 *  */
	if (utype == DOF_P2) {
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
	if (utype == DOF_P3) {
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
    dof_iLngh = NULL;
    for (i = 0; i < map->nlocal; i++) {
	dof_iLine[i] = -1;
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
	
	/*
	 * Special case: ice-sheet
	 * Line length if const
	 *
	 * */
	assert(ndof == ndof_layer);

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
     * Special case: ice-sheet
     * Line length if const
     *
     * */
    assert (nlocal == ndof_layer * nl);




    /*
     * step 3. Set Dofs near line block.
     * 
     * */
    phgDofFreeVEFMap2(&vef);

    /*
     * step 3'''. Line length statistic: 
     * 1. average length
     * 2. num of short lines
     *  */
#define LENN 1000
    verb = 1;
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
    }
#undef LENNN

    /*
     * Step 4. reorder line neigh dofs.
     *
     *  */
    nds[1] = 0;

    /*
     * Step 5. record left dofs
     *
     *  */
    nds[2] = 0;

    verb = 1;
    phgInfo(verb, "# Line dofs: %d\n", nds[0]);
    phgInfo(verb, "# Lngh dofs: %d\n", nds[1]);
    phgInfo(verb, "# Left dofs: %d\n", nds[2]);

    /* save line block info */
    /* record in ML */
    bk = phgCalloc(1, sizeof(*bk));
    memcpy(bk->nds, nds, sizeof(nds));
    bk->dof_type = u0->type;
    bk->mat_cnnt = mat_cnnt;
    bk->map_cnnt = map;
    bk->nline = nl;	
    bk->line_dofs = line_dofs;	
    bk->line_dof0 = line_dof0;	
    bk->line_ndof = line_ndof;	

    phgMatDestroy(&mat_fill);
    phgDofFree(&u);
    phgDofFree(&Coord);
    return bk;
}

void destroy_line_block(MG_BLOCK_DOFS **bk_ptr)
{
    MG_BLOCK_DOFS *bk = *bk_ptr;
    if (bk == NULL)
	return;

    /* phgMatDestroy(&bk->mat_cnnt);     */
    /* phgMapDestroy(&bk->map_cnnt);  */
    /* phgFree(bk->nds); */
    /* phgFree(bk->line_dofs); */
    /* phgFree(bk->line_dof0); */
    /* phgFree(bk->line_ndof); */
    /* phgFree(bk); */

    *bk_ptr = NULL;
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
factorize_line_block(GRID *g, MG_BLOCK_DOFS *bk, MAT *A, MAT_FACTOR **mf_ptr) 
{
    DOF_TYPE *utype = bk->dof_type;
    int order = utype->order;
    //MAP *map = A->rmap;
    int dim = 1;
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
    int ngh_size = 0, row_size, row_diag, nly_b;
    MAT_FACTOR *mf;
    int verb = 3;
    
    Unused(g);
    assert(bk != NULL);

    mf = phgCalloc(1, sizeof(*mf));

    phgPrintf("\n         Factorize %s, ", bk->dof_type->name);
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
    assert(dim == 1);
    if (utype->name[0] == 'P')
	ngh_size = bk->neigh_size = (order + 1) * dim - 1;
    else
	phgError(1, "Unsupported DOF type: %s!\n", utype->name);

    row_size = 2*ngh_size + ngh_size + 1; /* LDAB >= 2*KL+KU+1 */
    row_diag = 2*ngh_size;		  /* diagnol position of row */
    poc_size = 0;		/* # of diag block entries */
    assert(2*ngh_size + 1 < MAX_DOF_NEIGHBOUR);
    phgPrintf("neigh size: %d\n", ngh_size);

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

    mf->line_mat0 = line_mat0;
    mf->line_matv = line_matv;
    mf->line_piv0 = line_piv0;
    mf->line_pivs = line_pivs;
    mf->line_poc0 = line_poc0;
    mf->line_pocs = line_pocs;
    mf->factorized = TRUE;
    *mf_ptr = mf;    
}


void 
free_fact_line_block(MAT_FACTOR **mf_ptr) 
{
    MAT_FACTOR *mf = *mf_ptr;

    if (mf != NULL) {
#define FREE_SET(pt) phgFree(pt); pt = NULL;
	FREE_SET(mf->line_mat0);
	FREE_SET(mf->line_matv);
	FREE_SET(mf->line_piv0);
	FREE_SET(mf->line_pivs);
	FREE_SET(mf->line_poc0);
	FREE_SET(mf->line_pocs);
	mf->factorized = FALSE;
#undef FREE_SET
    }

    phgFree(mf);
    mf_ptr = NULL;
}



