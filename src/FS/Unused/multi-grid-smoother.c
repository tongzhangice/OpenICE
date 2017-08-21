#include "multi-grid.h"
#define _p mg_params
#define _mgp (mg->mg_params)
/*
 * Multi Grid smoother
 *
 * TODO:
 *   Implement smoother from BoomerAMG, Hypre.
 *
 * TODO:   
 *   Schwarz type smoother, in the future.
 *
 * Note:
 *   Smoother for symetric positive problems may not work with
 *   convection-diffusion problems, further test is needed !
 *   
 * */


/*
 * Jacobi smoother.
 *
 *   Implement using phgMatvec
 *   ref: phgSolverJacobi, solver-phg.c
 * */
void 
mg_Jacobi(MAT *A, VEC *x, VEC *b, int nsmooth, void *ctx)
{
    VEC *r;
    INT i, k; 
    FLOAT dx, omega = _p->smooth_damp;

    if (A->diag == NULL)
	phgMatSetupDiagonal(A);

    r = phgMapCreateVec(b->map, 1);
    for (k = 0; k < nsmooth; k++) {
	phgVecCopy(b, &r);
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
	for (i = 0; i < A->cmap->nlocal; i++) {
	    assert(fabs(A->diag[i]) > 1e-12);
	    dx = r->data[i] / A->diag[i]; 
	    x->data[i] += omega * dx;
	}
    }
    phgVecDestroy(&r);
    
    return;
}


/*
 * Jacobi smoother2:
 *   Implement using details of matvec
 *   ref: matvec in matvec.c
 * 
 *  */
void 
mg_Jacobi2(MAT *A, VEC *x, VEC *b, int nsmooth, void *ctx)
{
    INT i, j, k, n, *pc, *pc_offp;
    FLOAT *pd, *pd_offp, *vx, *vx0, *vb;
    FLOAT sum, omega = _p->smooth_damp;;
    VEC *x0;
#if USE_MPI
    FLOAT *offp_data = NULL;
#endif	/* USE_MPI */

    MagicCheck(VEC, x);
    MagicCheck(VEC, b);
    assert(x == NULL || x->nvec == 1);
    assert(b == NULL || b->nvec == 1);
    if (x != NULL && !x->assembled)
	phgVecAssemble(x);
    if (b != NULL && !b->assembled)
	phgVecAssemble(b);
    x0 = phgMapCreateVec(x->map, 1);
    
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

#if USE_MPI
    if (A->cmap->nprocs > 1) {
	phgMapScatterBegin(A->cinfo, x->nvec, x->data, offp_data);
	phgMapScatterEnd(A->cinfo, x->nvec, x->data, offp_data);
    }
#endif	/* USE_MPI */

    /* iteration */
    for (k = 0; k < nsmooth; k++) {
	phgVecCopy(x, &x0);

	/* multiply with local data */
	vx = x->data;
	vx0 = x0->data;
	vb = b->data;
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
	    FLOAT aa = 0., dx;
	    /* x_i = (b_i - \sum_{j ~= i} a_ij * x_j) / a_ii */

	    sum = vb[i];
	    /* local data */
	    if ((n = *(pc++)) != 0) {
		for (j = 0; j < n; j++) {
		    jcol = *(pc++);
		    if (jcol != i) {
			sum -= *(pd++) * vx0[jcol];
		    } else {
			aa = *(pd++);
			assert(fabs(aa) > 1e-14);
		    }
		}
	    }
	    
	    /* remote data */
	    if (pc_offp != NULL && (n = *(pc_offp++)) != 0) {
		for (j = 0; j < n; j++) {
		    jcol = *(pc_offp++);
		    sum -= *(pd_offp++) * offp_data[jcol];
		}
	    }
	    
	    dx = sum / aa - vx[i];
	    vx[i] += omega * dx;
	}

#if USE_MPI
	if (A->cmap->nprocs > 1) {
	    phgMapScatterBegin(A->cinfo, x->nvec, x->data, offp_data);
	    phgMapScatterEnd(A->cinfo, x->nvec, x->data, offp_data);
	}
#endif	/* USE_MPI */
    }


    phgVecDestroy(&x0);
#if USE_MPI
    phgFree(offp_data);
#endif	/* USE_MPI */

    return;
}


/*
 * GaussSidel smoother:
 *   1. exchange off proc data
 *   2. smooth local dof
 *   
 *   ref: matvec in matvec.c
 *  */
void 
mg_GaussSidel(MAT *A, VEC *x, VEC *b, int nsmooth, void *ctx)
{
    INT i, j, k, n, *pc, *pc_offp, nlocal;
    size_t *ps;
    FLOAT *pd, *pd_offp, *vx, *vb;
    FLOAT sum, omega = _p->smooth_damp;
#if USE_MPI
    FLOAT *offp_data = NULL;
#endif	/* USE_MPI */

    MagicCheck(VEC, x);
    MagicCheck(VEC, b);
    assert(x == NULL || x->nvec == 1);
    assert(b == NULL || b->nvec == 1);
    if (x != NULL && !x->assembled)
	phgVecAssemble(x);
    if (b != NULL && !b->assembled)
	phgVecAssemble(b);
    
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

    /* iteration */
    for (k = 0; k < nsmooth; k++) {

	/* multiply with local data */
	vx = x->data;
	vb = b->data;
	pc = A->packed_cols;
	pd = A->packed_data;
	ps = A->packed_ind;
	if (A->cmap->nprocs > 1) {
	    pc_offp = A->packed_cols + A->rmap->nlocal + A->nnz_d;
	    pd_offp = A->packed_data + A->nnz_d;
	} else {
	    pc_offp = NULL;
	    pd_offp = NULL;
	}

	for (i = 0; i < A->rmap->nlocal;  i++) {
	    INT jcol;
	    FLOAT aa = 0., dx;
	    /* x_i = (b_i - \sum_{j ~= i} a_ij * x_j) / a_ii */

	    sum = vb[i];
	    /* local data */
	    if ((n = *(pc++)) != 0) {
		assert(A->packed_cols + PACK_COL(ps, i) + 1== pc);
		assert(A->packed_data + PACK_DAT(ps, i) == pd);
		for (j = 0; j < n; j++) {
		    jcol = *(pc++);
		    if (jcol != i) {
			sum -= *(pd++) * vx[jcol];
		    } else {
			aa = *(pd++);
			assert(fabs(aa) > 1e-14);
		    }
		}
	    }
	    
	    /* remote data */
	    if (pc_offp != NULL && (n = *(pc_offp++)) != 0) {
		assert(A->packed_cols + PACK_COL_OFFP(ps, i, nlocal) + 1 == pc_offp);
		assert(A->packed_data + PACK_DAT_OFFP(ps, i, nlocal) == pd_offp);
		for (j = 0; j < n; j++) {
		    jcol = *(pc_offp++);
		    sum -= *(pd_offp++) * offp_data[jcol];
		}
	    }
	    
	    dx = sum / aa - vx[i];
	    vx[i] += omega * dx;
	}

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

/*
 * Gauss-Sidel smoother2:
 *
 *   1. interior dof (!REMOTE) is smoothed, and then offp data is updated; 
 *   2. proc boundary dof (REMOTE) is smoothed, and then offp data is updated.
 *
 * Note: need types_vec.
 *  */
void 
mg_GaussSidel2(MAT *A, VEC *x, VEC *b, int nsmooth, void *ctx)
{
    MG_LEVEL *ml = (MG_LEVEL *)ctx;
    INT i, j, k, n, *pc, *pc_offp;
    FLOAT *pd, *pd_offp, *vx, *vb;
    FLOAT sum, omega = _p->smooth_damp;;
#if USE_MPI
    FLOAT *offp_data = NULL;
#endif	/* USE_MPI */

    MagicCheck(VEC, x);
    MagicCheck(VEC, b);
    assert(x == NULL || x->nvec == 1);
    assert(b == NULL || b->nvec == 1);
    if (x != NULL && !x->assembled)
	phgVecAssemble(x);
    if (b != NULL && !b->assembled)
	phgVecAssemble(b);
    
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

#if USE_MPI
    if (A->cmap->nprocs > 1) {
	phgMapScatterBegin(A->cinfo, x->nvec, x->data, offp_data);
	phgMapScatterEnd(A->cinfo, x->nvec, x->data, offp_data);
    }
#endif	/* USE_MPI */

    /* iteration */
    for (k = 0; k < nsmooth; k++) {

	/* multiply with local data */
	vx = x->data;
	vb = b->data;

	/* First, interior dof (!REMOTE) is smoothed, and then offp data is updated */
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
	    FLOAT aa = 0., dx;

	    if (ml->types_vec[i] & REMOTE) {
		if ((n = *(pc++)) != 0) {
		    pc += n; 
		    pd += n;
		}
		if (pc_offp != NULL && (n = *(pc_offp++)) != 0) {
		    pc_offp += n; 
		    pd_offp += n;
		}
		continue;
	    }

	    sum = vb[i];
	    /* local data */
	    if ((n = *(pc++)) != 0) {
		for (j = 0; j < n; j++) {
		    jcol = *(pc++);
		    if (jcol != i) {
			sum -= *(pd++) * vx[jcol];
		    } else {
			aa = *(pd++);
			assert(fabs(aa) > 1e-14);
		    }
		}
	    }
	    
	    /* remote data */
	    if (pc_offp != NULL && (n = *(pc_offp++)) != 0) {
		for (j = 0; j < n; j++) {
		    jcol = *(pc_offp++);
		    sum -= *(pd_offp++) * offp_data[jcol];
		}
	    }
	    
	    dx = sum / aa - vx[i];
	    vx[i] += omega * dx;
	}

#if USE_MPI
	if (A->cmap->nprocs > 1) {
	    phgMapScatterBegin(A->cinfo, x->nvec, x->data, offp_data);
	    phgMapScatterEnd(A->cinfo, x->nvec, x->data, offp_data);
	}
#endif	/* USE_MPI */

	/* Second, proc boundary dof (!REMOTE) is smoothed, and then offp data is updated */
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
	    FLOAT aa = 0., dx;

	    if (!(ml->types_vec[i] & REMOTE)) {
		if ((n = *(pc++)) != 0) {
		    pc += n; 
		    pd += n;
		}
		if (pc_offp != NULL && (n = *(pc_offp++)) != 0) {
		    pc_offp += n; 
		    pd_offp += n;
		}
		continue;
	    }

	    sum = vb[i];
	    /* local data */
	    if ((n = *(pc++)) != 0) {
		for (j = 0; j < n; j++) {
		    jcol = *(pc++);
		    if (jcol != i) {
			sum -= *(pd++) * vx[jcol];
		    } else {
			aa = *(pd++);
			assert(fabs(aa) > 1e-14);
		    }
		}
	    }
	    
	    /* remote data */
	    if (pc_offp != NULL && (n = *(pc_offp++)) != 0) {
		for (j = 0; j < n; j++) {
		    jcol = *(pc_offp++);
		    sum -= *(pd_offp++) * offp_data[jcol];
		}
	    }
	    
	    dx = sum / aa - vx[i];
	    vx[i] += omega * dx;
	}

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

/*
 * Gauss-Sidel smoother for vector:
 *
 * As GS: 
 *   1. exchange off proc data
 *   2. smooth local dof
 *
 * Assumption:
 *   1. unknow dof is vector of dim 3
 *   2. Matrix block for vector is 
 *      | a1  -b   0 |
 *      |  b  a2   0 |
 *      |  0   0  a3 |
 * 
 *  */
void 
mg_GaussSidel_vec(MAT *A, VEC *x, VEC *b, int nsmooth, void *ctx)
{
    INT i, j, k, l, n, *pc, *pc0, *pc_offp, nlocal;
    FLOAT *pd, *pd0, *pd_offp, *vx, *vb;
    size_t *ps; 
    FLOAT sum[Dim], omega = _p->smooth_damp;
#if USE_MPI
    FLOAT *offp_data = NULL;
#endif	/* USE_MPI */

    MagicCheck(VEC, x);
    MagicCheck(VEC, b);
    assert(x == NULL || x->nvec == 1);
    assert(b == NULL || b->nvec == 1);
    if (x != NULL && !x->assembled)
	phgVecAssemble(x);
    if (b != NULL && !b->assembled)
	phgVecAssemble(b);
    
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

    if (A->cmap->nlocal % Dim != 0)
	phgError(1, "%s: assume vector dof of dim 3!\n", __FUNCTION__);

#if USE_MPI
    if (A->cmap->nprocs > 1) {
	phgMapScatterBegin(A->cinfo, x->nvec, x->data, offp_data);
	phgMapScatterEnd(A->cinfo, x->nvec, x->data, offp_data);
    }
#endif	/* USE_MPI */

    /* iteration */
    for (l = 0; l < nsmooth; l++) {
	INT i_start, i_add;

	/* multiply with local data */
	vx = x->data;
	vb = b->data;
	pc0 = A->packed_cols;
	pd0 = A->packed_data;
	ps = A->packed_ind;
	nlocal = A->rmap->nlocal;

	/*
	 * lexicographic order: low to high
	 * Note: low->high and high->low alternatively does not help.
	 * */
	if (TRUE || l % 2 == 0) {
	    i_start = 0;
	    i_add = Dim;
	} else {
	    i_start = nlocal - Dim;
	    i_add = -Dim;
	}

	//for (i = i_start; i < A->rmap->nlocal && i >= 0;  i += i_add) {
	for (i = 0; i < nlocal ;  i += Dim) {
	    INT jcol;
	    FLOAT aa[Dim][Dim], det, dx[Dim];

	    memset(aa, 0, sizeof(aa));
	    sum[0] = vb[i  ];
	    sum[1] = vb[i+1];
	    sum[2] = vb[i+2];

	    /* local data */
	    pc = pc0 + PACK_COL(ps, i);
	    pd = pd0 + PACK_DAT(ps, i);
	    for (k = 0; k < Dim; k++) {
		if ((n = *(pc++)) != 0) {
		    for (j = 0; j < n; j++) {
			jcol = *(pc++);
			if (jcol < i || jcol > i+2) {
			    sum[k] -= *(pd++) * vx[jcol];
			} else { /* offD, jcol = i,i+1,i+2 */
			    aa[k][jcol - i] = *(pd++);
			}
		    }
		}
	    }

	    /* remote data */
	    if (A->cmap->nprocs > 1) {
		pc_offp = pc0 + PACK_COL_OFFP(ps, i, nlocal);
		pd_offp = pd0 + PACK_DAT_OFFP(ps, i, nlocal);
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
	    dx[0] = (aa[1][1] * sum[0] - aa[0][1] * sum[1]) * det - vx[i  ];
	    dx[1] = (aa[0][0] * sum[1] - aa[1][0] * sum[1]) * det - vx[i+1];
	    dx[2] = (1./aa[2][2] * sum[2]) - vx[i+2];
	    vx[i  ] += omega * dx[0];
	    vx[i+1] += omega * dx[1];
	    vx[i+2] += omega * dx[2];
	}

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

/* -------------------------------------------------------------------------------- */
/* Register smoother */

static const char *smoother_name[] = {
    "Jacobi", "Jacobi2", "GS", "GS1", "GSv", "GSline", NULL};
static MG_SMOOTHER smoother_list[] = {
    mg_Jacobi, mg_Jacobi2, mg_GaussSidel, mg_GaussSidel2, mg_GaussSidel_vec, 
    mg_GaussSidel_line};
static int smoother_index = 0;	/* default Jacobi */

void
phgMultiGridSmootherRegister(void)
{
    phgOptionsRegisterKeyword("mg_smoother", MG_PREFIX"smoother",
			      smoother_name, &smoother_index);
}


/*
 * Get smoother.
 * Smoother could be a solver, in that case, it's stored in struct ml.
 * */
MG_SMOOTHER
phgMultiGridGetSmoother(MULTI_GRID *mg, int level, int type)
{
#if 0    
    switch (type) {
    case DOWN_CYCLE:
	return mg_Jacobi;
    case UP_CYCLE:
	return mg_Jacobi;
    case COARSEST:
	return mg_Jacobi;
    case FINEST:
	return mg_Jacobi;
    default:
	return mg_Jacobi;
    }
#else
    return smoother_list[smoother_index];
#endif
}

