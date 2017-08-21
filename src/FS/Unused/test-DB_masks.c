/* Steady-state incompressible Navier-Stokes equations:
 *
 *   -1/Re \laplace u + (u \cdot \grad) u + \grad p = f,
 *   \div u = 0,
 *   u = w on \partial \Omega_D, 
 *   1/Re \frac{\partial u}{\partial n} - n p = g on \partial \Omega_N.
 *
 *   The nonliner algebraic system of equations resulting from finite element
 *   discretization is solved using Picard or Newton iterations. In each
 *   nonlinear iteration, the corresponding linearized residual system is
 *   solved using GMRES with a pressure convection-diffusion preconditioner
 *   which is written as:
 *
 *   |  F       -Bt       |^{-1}
 *   |                    |
 *   |  0  -Ap Fp^{-1} Qp |
 *
 *   where
 *   	Fp <==>	-1/Re \laplace + u \cdot \grad
 *		the pressure convection-diffusion matrix.
 *   	Ap <==> the discrete laplace operator
 *   	Qp <==> the lumped pressure mass matrix
 *   All three operators above are in the pressure space.
 *
 *   Reference:
 *	D. Kay, D. Loghin, and A. Wathen, A preconditioner for the steady-state
 *	Navier-Stokes equations, SIAM J. Sci. Comput., 24 (2002), pp. 237-256.
 * 
 * $Id: navier-stokes.c,v 1.51 2009/01/13 02:17:09 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <strings.h>	/* bzero() */
#include <math.h>
#include <stdarg.h>

/* ********************** */
/* *** define equations * */
/* ********************** */
static FLOAT nu;
#ifndef TEST_CASE
# define TEST_CASE 3
/* # define TEST_CASE 100		/\* DRIVEN_CAVITY *\/ */
/* # define TEST_CASE 101*/ /* CYLINDER */
#endif
#include "test-DB_masks.h"

FUNCV_DEF(func_u, u_);
FUNCV_DEF(func_p, p_);
FUNCV_DEF(func_gradu, gradu_);
FUNCV_DEF(func_f, f_);
FUNCV_DEF(func_gx, g1_);
FUNCV_DEF(func_gy, g2_);
FUNCV_DEF(func_gz, g3_);

typedef enum {PICARD, NEWTON} LTYPE;	/* linearization type */

/* enclosed_flow
 *   1. system is singular, pin a node. (Unimplemented)
 *   2. Procond Fp, Ap have no bdry.
 * */
/*#define ENCLOSED_FLOW*/
#ifdef ENCLOSED_FLOW
#define PINNED_NODE BDRY_USER1
#endif

#define OUTPUT_VTK 1

/* Gloable viriables */
static SOLVER *solver_Ap = NULL, *solver_Qp = NULL, *solver_F = NULL;
static MAT *matFp, *matBt, *matF, *matFu = NULL;

static double
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

/* Too many dofs, pack them */
void
packDof(DOF **dofs, int ndof, ...)
{
    int i;
    va_list ap;

    va_start(ap, ndof);
    for (i = 0; i < ndof; i++)
	dofs[i] = va_arg(ap, DOF *);
    va_end(ap);
    return;
}

void
unpackDof(DOF **dofs, int ndof, ...)
{
    int i;
    va_list ap;
    DOF **pdof;

    va_start(ap, ndof);
    for (i = 0; i < ndof; i++) {
	pdof = va_arg(ap, DOF **);
	*pdof = dofs[i];
    }
    va_end(ap);
    return;
}

void
packMat(MAT **mats, int nmat, ...)
{
    int i;
    va_list ap;

    va_start(ap, nmat);
    for (i = 0; i < nmat; i++)
	mats[i] = va_arg(ap, MAT *);
    va_end(ap);
    return;
}

void
unpackMat(MAT **mats, int nmat, ...)
{
    int i;
    va_list ap;
    MAT **pmat;

    va_start(ap, nmat);
    for (i = 0; i < nmat; i++) {
	pmat = va_arg(ap, MAT **);
	*pmat = mats[i];
    }
    va_end(ap);
    return;
}

/******************************/
/* *** Quadrature routine *** */
/******************************/
# define BasisOrder(u, e, i) (!DofIsHP(u) ? (u)->type->order : \
	      (u)->hp->info->types[(u)->hp->elem_order[e->index]]->order)
/* *******************************************************
 * ( (u.\grad) phi_j, phi_i )
 * return 1 value
 *  */
FLOAT
phgQuadDofDotGradBasBas(SIMPLEX *e, DOF *u, DOF *v, int m, int n, int order)
{
    int i;
    const FLOAT *g1, *g2, *gu, *w, *lambda;
    FLOAT d, d0;
    QUAD *quad;

    assert(u->dim == 3);

    if (order < 0)
	order = DofTypeOrder(u, e) +
	      BasisOrder(v, e, m) - 1 + BasisOrder(v, e, n);
    if (order < 0)
	order = 0;
    quad = phgQuadGetQuad3D(order);

    gu = phgQuadGetDofValues(e, u, quad);
    g1 = phgQuadGetBasisGradient(e, v, m, quad);
    g2 = phgQuadGetBasisValues(e, v, n, quad);
    d = 0.;
    lambda = quad->points;
    w = quad->weights;
    for (i = 0; i < quad->npoints; i++) {
	d0 = 0.;
	d0 += (*(gu++)) * (*(g1++)) * (*(g2));	/* u dphi_m/dx phi_n */
	d0 += (*(gu++)) * (*(g1++)) * (*(g2));	/* v dphi_m/dy phi_n */
	d0 += (*(gu++)) * (*(g1++)) * (*(g2));	/* w dphi_m/dz phi_n */
	g2++;
	d += d0 * (*(w++));
	lambda += Dim + 1;
    }

    return d * phgGeomGetVolume(u->g, e);
}

/* ****************************************
 * ( d[uvw]/d[xyz]*phi_j, phi_i )
 * return static float[9] managed by user.
 *  */
void
phgQuadGradDofBasDotBas(SIMPLEX *e, DOF *gradu, DOF *u, int m, int n,
			int order, FLOAT *value)
{
    int i;
    const FLOAT *g1, *g2, *gu, *w, *lambda;
    FLOAT d[9], *pd, vol;
    QUAD *quad;

    assert(gradu->dim == 9 && u->dim == 3);

    if (order < 0)
	order = DofTypeOrder(gradu, e) +
	       BasisOrder(u, e, m) + BasisOrder(u, e, n) ;
    if (order < 0)
	order = 0;
    quad = phgQuadGetQuad3D(order);

    gu = phgQuadGetDofValues(e, gradu, quad);
    g1 = phgQuadGetBasisValues(e, u, m, quad);
    g2 = phgQuadGetBasisValues(e, u, n, quad);
    lambda = quad->points;
    w = quad->weights;
    for (i = 0; i < 9; i++)
	d[i] = 0;
    for (i = 0; i < quad->npoints; i++) {
	pd = d;
	(*pd++) += (*(gu++)) * (*g1) * (*g2) * (*w);	/* du/dx phi_m phi_n */
	(*pd++) += (*(gu++)) * (*g1) * (*g2) * (*w);	/* du/dy phi_m phi_n */
	(*pd++) += (*(gu++)) * (*g1) * (*g2) * (*w);	/* du/dz phi_m phi_n */

	(*pd++) += (*(gu++)) * (*g1) * (*g2) * (*w);	/* dv/dx phi_m phi_n */
	(*pd++) += (*(gu++)) * (*g1) * (*g2) * (*w);	/* dv/dy phi_m phi_n */
	(*pd++) += (*(gu++)) * (*g1) * (*g2) * (*w);	/* dv/dz phi_m phi_n */

	(*pd++) += (*(gu++)) * (*g1) * (*g2) * (*w);	/* dw/dx phi_m phi_n */
	(*pd++) += (*(gu++)) * (*g1) * (*g2) * (*w);	/* dw/dy phi_m phi_n */
	(*pd++) += (*(gu++)) * (*g1) * (*g2) * (*w);	/* dw/dz phi_m phi_n */

	g1++;
	g2++;
	w++;
	lambda += Dim + 1;
    }

    vol = phgGeomGetVolume(u->g, e);
    for (i = 0; i < 9; i++)
	value[i] = d[i] * vol;

    return;
}

/* **********************
 * \grad phi_j \times \phi_i
 * return 3 values
 *  */
void
phgQuadGradBasDotBas3(SIMPLEX *e, DOF *s, int m, DOF *v, int n,
		      FLOAT *values, int order)
{
    int i, j, nvalues = DofTypeDim(v);
    const FLOAT *g1, *g2, *w, *lambda;
    FLOAT d1, d2, d3, dd1, dd2, dd3;
    QUAD *quad;

    assert(!SpecialDofType(s->type) && !SpecialDofType(v->type));

    if (nvalues != DofTypeDim(s))
	phgError(1, "%s:%d, dimensions mismatch: grad(%s) <==> (%s)\n",
		 __FILE__, __LINE__, s->name, v->name);

    if (order < 0)
	order = BasisOrder(s, e, m) - 1 + BasisOrder(v, e, n);
    if (order < 0)
	order = 0;
    quad = phgQuadGetQuad3D(order);

    g1 = phgQuadGetBasisGradient(e, s, m, quad);
    g2 = phgQuadGetBasisValues(e, v, n, quad);
    dd1 = 0.;
    dd2 = 0.;
    dd3 = 0.;
    lambda = quad->points;
    w = quad->weights;
    for (i = 0; i < quad->npoints; i++) {
	d1 = d2 = d3 = 0.;
	for (j = 0; j < nvalues; j++) {
	    d1 += *(g1++) * *g2;
	    d2 += *(g1++) * *g2;
	    d3 += *(g1++) * *g2;
	    g2++;
	}
	dd1 += d1 * *w;
	dd2 += d2 * *w;
	dd3 += d3 * *w;
	w++;
	lambda += Dim + 1;
    }
    values[0] = dd1 * phgGeomGetVolume(s->g, e);
    values[1] = dd2 * phgGeomGetVolume(s->g, e);
    values[2] = dd3 * phgGeomGetVolume(s->g, e);
    return;
}

/* ***************************
 * ( (u.\grad) u, phi_i )
 * return 3 values.
 *  */
void
phgQuadDofDotGradDofBas(SIMPLEX *e, DOF *u, DOF *gradu,
			int n, int order, FLOAT *values)
{
    int i;
    const FLOAT *g2, *gu, *ggu, *w, *lambda;
    FLOAT d1, d2, d3, dd1, dd2, dd3, vol;
    QUAD *quad;

    assert(gradu->dim == 9 && u->dim == 3);

    if (order < 0)
	order = DofTypeOrder(u, e) +
	       DofTypeOrder(gradu, e) + BasisOrder(u, e, n);
    if (order < 0)
	order = 0;
    quad = phgQuadGetQuad3D(order);

    gu = phgQuadGetDofValues(e, u, quad);
    ggu = phgQuadGetDofValues(e, gradu, quad);
    g2 = phgQuadGetBasisValues(e, u, n, quad);
    dd1 = dd2 = dd3 = 0.;
    lambda = quad->points;
    w = quad->weights;
    for (i = 0; i < quad->npoints; i++) {
	d1 = d2 = d3 = 0.;

	d1 += *(gu) * (*(ggu++)) * (*(g2));	/* u du/dx phi_n */
	d1 += *(gu + 1) * (*(ggu++)) * (*(g2));	/* v du/dy phi_n */
	d1 += *(gu + 2) * (*(ggu++)) * (*(g2));	/* w du/dz phi_n */

	d2 += *(gu) * (*(ggu++)) * (*(g2));	/* u dv/dx phi_n */
	d2 += *(gu + 1) * (*(ggu++)) * (*(g2));	/* v dv/dy phi_n */
	d2 += *(gu + 2) * (*(ggu++)) * (*(g2));	/* w dv/dz phi_n */

	d3 += *(gu) * (*(ggu++)) * (*(g2));	/* u dw/dx phi_n */
	d3 += *(gu + 1) * (*(ggu++)) * (*(g2));	/* v dw/dy phi_n */
	d3 += *(gu + 2) * (*(ggu++)) * (*(g2));	/* w dw/dz phi_n */

	dd1 += d1 * (*w);
	dd2 += d2 * (*w);
	dd3 += d3 * (*w);

	gu += 3;
	g2++;
	w++;
	lambda += Dim + 1;
    }

    vol = phgGeomGetVolume(u->g, e);
    values[0] = dd1 * vol;
    values[1] = dd2 * vol;
    values[2] = dd3 * vol;
    return;
}

/* ***********************
 * p \times \grad \phi_m
 * return 3 values.
 * */
void
phgQuadDofAGradBas3(SIMPLEX *e, DOF *p, DOF *v, int m, int order, FLOAT *values)
{
    int i, j;
    const FLOAT *g1, *g2, *w;
    FLOAT d1, d2, d3, vol;
    QUAD *quad;

    assert(!SpecialDofType(v->type));
    assert(p->dim == 1);

    if (order < 0) {
	order = ((j = DofTypeOrder(p, e)) >= 0 ? j : BasisOrder(v, e, m));
	order += BasisOrder(v, e, m) - 1;
    }
    quad = phgQuadGetQuad3D(order);

    g1 = phgQuadGetDofValues(e, p, quad);
    g2 = phgQuadGetBasisGradient(e, v, m, quad);
    w = quad->weights;

    d1 = d2 = d3 = 0.;
    for (i = 0; i < quad->npoints; i++) {
	d1 += *(g1) * (*g2++) * (*w);
	d2 += *(g1) * (*g2++) * (*w);
	d3 += *(g1) * (*g2++) * (*w);
	g1++;
	w++;
    }

    vol = phgGeomGetVolume(p->g, e);
    values[0] = d1 * vol;
    values[1] = d2 * vol;
    values[2] = d3 * vol;
    return;
}

/* *****************************
 * (grad(u_3), grad(\phi))
 * return 3 values.
 * */
void
phgQuadGradDof3GradBas(SIMPLEX *e, DOF *gradu, DOF *u, int m, int order,
		       FLOAT *values)
{
    int i;
    const FLOAT *g2, *gu, *w, *lambda;
    FLOAT d1, d2, d3, dd1, dd2, dd3, vol;
    QUAD *quad;

    assert(gradu->dim == 9 && u->dim == 3);

    if (order < 0)
	order = 2 * BasisOrder(u, e, m) - 1;
    if (order < 0)
	order = 0;
    quad = phgQuadGetQuad3D(order);

    gu = phgQuadGetDofValues(e, gradu, quad);
    g2 = phgQuadGetBasisGradient(e, u, m, quad);
    dd1 = dd2 = dd3 = 0.;
    lambda = quad->points;
    w = quad->weights;
    for (i = 0; i < quad->npoints; i++) {
	d1 = d2 = d3 = 0.;

	d1 += (*(gu++)) * (*(g2));
	d1 += (*(gu++)) * (*(g2 + 1));
	d1 += (*(gu++)) * (*(g2 + 2));

	d2 += (*(gu++)) * (*(g2));
	d2 += (*(gu++)) * (*(g2 + 1));
	d2 += (*(gu++)) * (*(g2 + 2));

	d3 += (*(gu++)) * (*(g2));
	d3 += (*(gu++)) * (*(g2 + 1));
	d3 += (*(gu++)) * (*(g2 + 2));

	dd1 += d1 * (*w);
	dd2 += d2 * (*w);
	dd3 += d3 * (*w);

	g2 += 3;
	w++;
	lambda += Dim + 1;
    }

    vol = phgGeomGetVolume(u->g, e);
    values[0] = dd1 * vol;
    values[1] = dd2 * vol;
    values[2] = dd3 * vol;
    return;
}

/* Preconditioning procedure using diag(matFu, matFu, matFu) */
static void
pc_proc1(SOLVER *pc_solver, VEC *b0, VEC **x0)
{
    int Nu = matF->rmap->nlocal;
    int k, nits;
    FLOAT res;
    INT i;
    VEC *tmp = phgMapCreateVec(pc_solver->rhs->map, 1);
    INT Nu0 = tmp->map->nlocal;

    assert(pc_solver->mat == matFu);
    assert(Nu == Nu0 * 3);
    nits = 0;
    res = 0.;
    for (k = 0; k < 3; k++) {
	for (i = 0; i < Nu0; i++)
	    pc_solver->rhs->data[i] = b0->data[3 * i + k];
	pc_solver->rhs->assembled = TRUE;
	pc_solver->rhs_updated = TRUE;
	bzero(tmp->data, sizeof(*tmp->data) * Nu0);
	phgSolverVecSolve(pc_solver, FALSE, tmp);
	for (i = 0; i < Nu0; i++)
	    (*x0)->data[3 * i + k] = tmp->data[i];
	if (nits < pc_solver->nits)
	    nits = pc_solver->nits;
	if (res < pc_solver->residual)
	    res = pc_solver->residual;
    }
    phgVecDestroy(&tmp);
    pc_solver->nits = nits;
    pc_solver->residual = res;

    return;
}

/*
 *  Preconditioning procedure:
 *     p = Qp^-1 Fp Ap^-1 * q
 *     u = F^-1 ( bu - Bt * p )
 *  */
static void
pc_proc(SOLVER *pc_solver, VEC *b0, VEC **x0)
{
    GRID *g = matF->rmap->dofs[0]->g;
    VEC *xu, *xp, *xu2, *xp2;
    int Nu = matF->rmap->nlocal;
    int Np = solver_Ap->mat->rmap->nlocal;
    INT verb = phgVerbosity;
    FLOAT *rhsF, *rhsAp, *rhsQp;
    double t;

    if (phgVerbosity > 0)
	phgVerbosity--;

    /* save old rhs */
    rhsF = solver_F->rhs->data;
    rhsAp = solver_Ap->rhs->data;
    rhsQp = solver_Qp->rhs->data;

    xu = phgMapCreateVec(matF->rmap, 1);
    xu2 = phgMapCreateVec(matF->rmap, 1);
    xp = phgMapCreateVec(solver_Ap->mat->rmap, 1);
    xp2 = phgMapCreateVec(solver_Ap->mat->rmap, 1);

    /* Ap^-1 * q */
    t = phgGetTime(NULL);
    solver_Ap->rhs->data = xp->data;
    solver_Ap->rhs->assembled = TRUE;
    solver_Ap->rhs_updated = TRUE;
    memcpy(xp->data, b0->data + Nu, sizeof(*xp->data) * Np);
    bzero(xp2->data, sizeof(*xp->data) * Np);
    phgSolverVecSolve(solver_Ap, FALSE, xp2);
    memcpy(xp->data, xp2->data, sizeof(*xp->data) * Np);
    if (verb > 0)
	phgPrintf("\t    Ap: nits = %d, residual = %0.4le [%0.4lgMB %0.4lfs]\n",
		  solver_Ap->nits, (double)solver_Ap->residual,
		  phgMemoryUsage(g, NULL) / (1024.0 * 1024.0),
		  phgGetTime(NULL) - t);

    /* Fp Ap^-1 * q */
    phgMatVec(MAT_OP_N, 1.0, matFp, xp, 0., &xp2);
    memcpy(xp->data, xp2->data, sizeof(*xp->data) * Np);

    /* xp = Qp^-1 Fp Ap^-1 * q */
    t = phgGetTime(NULL);
    solver_Qp->rhs->data = xp->data;
    solver_Qp->rhs->assembled = TRUE;
    solver_Qp->rhs_updated = TRUE;
    bzero(xp2->data, sizeof(*xp->data) * Np);
    phgSolverVecSolve(solver_Qp, FALSE, xp2);
    if (verb > 0)
	phgPrintf("\t    Qp: nits = %d, residual = %0.4le [%0.4lgMB %0.4lfs]\n",
		  solver_Qp->nits, (double)solver_Qp->residual,
		  phgMemoryUsage(g, NULL) / (1024.0 * 1024.0),
		  phgGetTime(NULL) - t);

    /* xu = F^-1 ( bu - Bt * xp ) */
    t = phgGetTime(NULL);
    memcpy(xu->data, b0->data, sizeof(*xu->data) * Nu);
    phgMatVec(MAT_OP_N, -1.0, matBt, xp2, 1., &xu);

    if (solver_F->mat->cmap->nlocal == Nu) {
	/* use F in the preconditioning matrix */
	solver_F->rhs->data = xu->data;
	solver_F->rhs->assembled = TRUE;
	solver_F->rhs_updated = TRUE;
	bzero(xu2->data, sizeof(*xu->data) * Nu);
	phgSolverVecSolve(solver_F, FALSE, xu2);
    }
    else {
	pc_proc1(solver_F, xu, &xu2);
    }

    if (verb > 0)
	phgPrintf("\t     F: nits = %d, residual = %0.4le [%0.4lgMB %0.4lfs]\n",
		  solver_F->nits, (double)solver_F->residual,
		  phgMemoryUsage(g, NULL) / (1024.0 * 1024.0),
		  phgGetTime(NULL) - t);

    /* Copy xu, xp to x */
    memcpy((*x0)->data, xu2->data, sizeof(*xu->data) * Nu);
    memcpy((*x0)->data + Nu, xp2->data, sizeof(*xp->data) * Np);

    /* restore rhs */
    solver_F->rhs->data = rhsF;
    solver_Ap->rhs->data = rhsAp;
    solver_Qp->rhs->data = rhsQp;

    phgVecDestroy(&xu);
    phgVecDestroy(&xu2);
    phgVecDestroy(&xp);
    phgVecDestroy(&xp2);

    phgVerbosity = verb;

    return;
}


static char *
Btype2bit(BTYPE n)
{
    int i;
    static char bit_str[100];

    bit_str[0] = '@';
    for (i = 16; i > 0; i--, n >>= 1)
	bit_str[i] = (n % 2 == 0) ? 'o' : '1';
    bit_str[17] = '\0';

    return bit_str;
}

/****************************************************************
 * Build RHS which is the residual of the nonlinear system.
 ***************************************************************/
static void
build_rhs(SOLVER *solver, SOLVER *pc, DOF **dofs, MAT **mats)
{
    DOF *u = dofs[0], *p = dofs[1];
    DOF *f, *pbc, *gn[3], *gradu, *divu;
    int M = u->type->nbas;	/* num of bases of Velocity */
    int N = p->type->nbas;	/* num of bases of Pressure */
    int i, k;
    GRID *g = u->g;
    SIMPLEX *e;
    FLOAT bufu[M], resu[M][Dim], resp[N], tmp[9];
    INT Iu[M][Dim], Ip[N];
    int s;

    /* Unpack Dofs */
    unpackDof(dofs, 9, &u, &p, &gradu, &divu, &f, &pbc, &gn[0], &gn[1], &gn[2]);

    ForAllElements(g, e) {
	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    for (k = 0; k < Dim; k++)
		Iu[i][k] = phgMapE2L(solver->rhs->map, 0, e, i * Dim + k);

	for (i = 0; i < N; i++) {
	    Ip[i] = phgMapE2L(solver->rhs->map, 1, e, i);
	}

	/* Global Matrix */
	bzero(resu, sizeof(resu));
	bzero(resp, sizeof(resp));
	for (i = 0; i < M; i++) {
	    for (k = 0; k < Dim; k++) {
		/* Dirichle Boundary for velocity. */
		if (phgDofDirichletBC_(u, e, i*Dim + k, func_u, bufu, tmp, DOF_PROJ_NONE)) {
		    resu[i][k] = tmp[k];
		}
		else {		/* interior node or Neumann */
		    /* (f, \phi_i) */
		    phgQuadDofTimesBas(e, f, u, i, 10, tmp);
		    resu[i][k] = tmp[k];
		}
	    }
	} /* end of Block (1,1), (1,2) */

	/* Neumann Bdry */
	for (s = 0; s < NFace; s++) {
	    if (e->bound_type[s] & BDRY_MASK) {
		SHORT bases_u[NbasFace(u)];
		phgDofGetBasesOnFace(u, e, s, bases_u);
		/* printf("\n-------\n e[%3d] bound[%d]: %d\n", e->index, s,  */
		/*        e->bound_type[s]); */

		/* velocity outflow */
		for (i = 0; i < NbasFace(u); i++) {
		    INT ii = bases_u[i];
		    for (k = 0; k < Dim; k++) {
			if (phgDofDirichletBC_(u, e, ii*Dim + k, NULL, bufu, NULL, DOF_PROJ_NONE)) {
			    /* Dirichlet bas on Neumann face, do nothing */
			    //printf("*** DirichN [%d][%4d]\n", k, Iu[ii][k]);
			}
			else {
			    //printf("*** Neumann [%d][%4d]\n", k, Iu[ii][k]);
			    resu[ii][k] += phgQuadFaceDofDotBas(e, s, gn[k],
								DOF_PROJ_DOT, u,
								ii, QUAD_DEFAULT);
			}
		    }
		}		/* end of base on face */
		/* pressure outflow */
	    }			/* end of face neu */
	}			/* end of all neumann face in element */

	/* Global res */
	phgSolverAddRHSEntries(solver, M * Dim, Iu[0], &resu[0][0]);
	phgSolverAddRHSEntries(solver, N, Ip, &resp[0]);
    }				/* end element */

    /* solver->rhs_updated = FALSE;	/\* note: can be safely commented out *\/ */

    return;
}

/***************************************************************************
 * Build matrices *F, *Fu, *B, *Bt, *Fp, *Ap, and *Qp used * by the solvers.
 **************************************************************************/
static void
build_matrices(SOLVER *solver, SOLVER *pc, DOF **dofs, MAT **mats, LTYPE type)
{
    DOF *u = dofs[0], *p = dofs[1];
    DOF *f, *pbc, *gn[3], *gradu, *divu;
    int M = u->type->nbas;	/* num of bases of Velocity */
    int N = p->type->nbas;	/* num of bases of Pressure */
    int i, j, k, l;
    GRID *g = u->g;
    SIMPLEX *e;
    MAT *matB, *matC, *matAp, *matQp;
    FLOAT F[M][Dim][M][Dim], Fu[M][M],
	B[N][M][Dim], Bt[M][Dim][N],
	Ap[N][N], Fp[N][N], Qp[N][N],
	bufu[M], bufp[N], tmp[9];
    INT Iu[M][Dim], Ju[Dim][M], Iu0[M], Ip[N];

    /* Unpack Dofs */
    unpackDof(dofs, 9, &u, &p, &gradu, &divu, &f, &pbc, &gn[0], &gn[1], &gn[2]);
    unpackMat(mats, 4, &matB, &matC, &matAp, &matQp);

#if 0
#warning remove me!
SOLVER *s = phgSolverCreate(SOLVER_PCG, p, NULL);
#endif

    ForAllElements(g, e) {
	/* Map: Element -> system */
	for (i = 0; i < M; i++) {
	    for (k = 0; k < Dim; k++)
		Ju[k][i] = Iu[i][k] = phgMapE2L(matF->cmap, 0, e, i * Dim + k);
	    Iu0[i] = phgMapE2L(matFu->cmap, 0, e, i);
	}
	for (i = 0; i < N; i++) {
	    Ip[i] = phgMapE2L(matFp->cmap, 0, e, i);
	}

	bzero(F, sizeof(F));
	bzero(B, sizeof(B));
	bzero(Bt, sizeof(Bt));
	/* A V W */
	for (i = 0; i < M; i++) {
	    for (j = 0; j < M; j++) {
		FLOAT a, v, w[9];
		/* (\grad \phi_j) \cdot (\grad \phi_i) */
		a = nu * phgQuadGradBasDotGradBas(e, u, j, u, i, QUAD_DEFAULT);
		/* (u \cdot \grad) \phi_j \times \phi_i, 1 item */
		v = phgQuadDofDotGradBasBas(e, u, u, j, i, QUAD_DEFAULT);
		Fu[i][j] = a;
		if (type == PICARD) {
		    memset(w, 0, sizeof(w));
		}
		else {
		    /* \phi_j (\grad u) \times \phi_i, 9 items */
		    phgQuadGradDofBasDotBas(e, gradu, u, j, i, QUAD_DEFAULT, w);
		}
		for (k = 0; k < Dim; k++) {
		    F[i][k][j][k] = a;
		    for (l = 0; l < Dim; l++) {
			//F[i][k][j][l] = *(w + k * Dim + l);
		    }
		}
	    }
	}

	/* B Bt */
	for (i = 0; i < M; i++) {
	    for (j = 0; j < N; j++) {
		FLOAT bxyz[3];
		phgQuadGradBasDotBas3(e, u, i, p, j, bxyz, QUAD_DEFAULT);
		for (k = 0; k < Dim; k++)
		    Bt[i][k][j] = B[j][i][k] = -bxyz[k];
	    }
	}

	/* Ap Qp;  Fp(update) */
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		FLOAT ap;
		ap = phgQuadGradBasDotGradBas(e, p, j, p, i, QUAD_DEFAULT);
		Ap[i][j] = nu * ap;
		Fp[i][j] =
		    phgQuadDofDotGradBasBas(e, u, p, j, i, QUAD_DEFAULT) + nu * ap;
		Qp[i][j] = phgQuadBasDotBas(e, p, j, p, i, QUAD_DEFAULT);
	    }
	}

	/* Global Matrix */
	for (i = 0; i < M; i++) {
	    for (k = 0; k < Dim; k++) {
		/* Dirichle Boundary for velocity. */
		if (phgDofDirichletBC_(u, e, i*Dim + k, NULL, bufu, NULL, DOF_PROJ_NONE)) {
		    phgMatAddEntries(matF, 1, Iu[i] + k, M, Ju[k], bufu);
		    if (k == 0) phgMatAddEntries(matFu, 1, Iu0 + i, M, Iu0, bufu);
		}
		else {		/* interior node Or Neumann */
		    /* Matrix F */
		    phgMatAddEntries(matF, 1, Iu[i] + k, M*Dim, Iu[0], &(F[i][k][0][0]));
		    phgMatAddEntries(matBt, 1, &Iu[i][k], N, Ip, &Bt[i][k][0]);
		    /* Matrix Fu */
		    if (k == 0) phgMatAddEntries(matFu, 1, Iu0 + i, M, Iu0, &(Fu[i][0]));
		}
	    }
	}			/* end of Block (1,1), (1,2) */

	for (i = 0; i < N; i++)
	    phgMatAddEntries(matB, 1, Ip + i, M * Dim, Iu[0], &B[i][0][0]);

	/* Matrix Ap Qp Fp */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(pbc, e, i, NULL, bufp, tmp, DOF_PROJ_NONE)) {
		/* Dirichle Boundary for pressure PC. */
		phgMatAddEntries(matAp, 1, Ip + i, N, Ip, bufp);
		phgMatAddEntries(matFp, 1, Ip + i, N, Ip, bufp);
	    }
	    else {
		/* interior node Or Neumann */
		phgMatAddEntries(matAp, 1, Ip + i, N, Ip, Ap[i]);
		phgMatAddEntries(matFp, 1, Ip + i, N, Ip, Fp[i]);
	    }
#if 0
	    /* use only diagonal of the mass matrix in the preconditioner */
	    phgMatAddEntries(matQp, 1, Ip + i, 1, Ip + i, &Qp[i][i]);
#else
	    /* use full mass matrix in the preconditioner */
	    phgMatAddEntries(matQp, 1, Ip + i, N, Ip, Qp[i]);
#endif
#if 0
phgMatAddEntries(s->mat, 1, Ip + i, N, Ip, Qp[i]);
#endif
	}
    }				/* end element */

#if 0
VEC *v = phgMapCreateVec(s->rhs->map, 1);
memcpy(s->rhs->data, solver->rhs->data + DofGetDataCount(u), s->rhs->map->nlocal * sizeof(FLOAT));
s->rhs_updated = TRUE;
s->rhs->assembled = TRUE;
s->rtol = 1e-10;
s->maxit = 10000;
phgSolverVecSolve(s, TRUE, v);
phgPrintf("v = %e\n", phgVecNormInfty(v, 0, NULL));
phgVecDestroy(&v);
phgSolverDestroy(&s);
#endif

    return;
}

static void
estimate_error(DOF *u, DOF *p, DOF *f, DOF *gradu, DOF *divu, DOF *vortex,
	       DOF *e_H1)
/* compute error indicators L_H1(e_u) + L2(e_p). */
{
    GRID *g = u->g;
    SIMPLEX *e;
#if 0
    int i;
    DOF *jump, *residual, *tmp;
    FLOAT eta, d, h;
    FLOAT diam;

    /* RE = [[nu \grad u - p I]] */
    tmp = phgDofCopy(gradu, NULL, NULL, NULL);
    phgDofAFXPBY(-1.0, f_1to3, p, nu, &tmp);
    jump = phgQuadFaceJump(tmp, DOF_PROJ_DOT, "jumps", QUAD_DEFAULT);
    phgDofFree(&tmp);
    /* RT1 = f + nu \laplace u - (u \cdot \grad) u - \grad p
     * RT2 = \div u */
    tmp = phgDofDivergence(gradu, NULL, NULL, NULL);
    residual = phgDofGetSameOrderDG(u, -1, "residual 1");
    phgDofCopy(f, &residual, NULL, NULL);
    phgDofAXPBY(nu, tmp, 1.0, &residual);
    phgDofFree(&tmp);
    tmp = phgDofGradient(p, NULL, NULL, NULL);
    phgDofAXPY(-1.0, tmp, &residual);
    phgDofMM(MAT_OP_N, MAT_OP_T, 1, 3, 3, 1.0, u, -1, gradu, 1., &residual);
    phgDofFree(&tmp);

    ForAllElements(g, e) {
	diam = phgGeomGetDiameter(g, e);
	e->mark = 0;		/* clear refinement mmark */
	eta = 0.0;
	/* for each face F compute [grad_u \cdot n] */
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] & (DIRICHLET | NEUMANN))
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    d = *DofFaceData(jump, e->faces[i]);
	    eta += d * h;
	}

	eta +=
	    diam * diam * phgQuadDofDotDof(e, residual, residual, QUAD_DEFAULT)
			+ phgQuadDofDotDof(e, divu, divu, QUAD_DEFAULT);

	/* add curved boundary errors (FIXME: how to normalize?) */
	eta += phgGetBoundaryError(g, e);

	*DofElementData(e_H1, e->index) = eta;
    }
    phgDofFree(&jump);
    phgDofFree(&residual);
#else
    /* use \int |vortex|^2 as error indicator */
    ForAllElements(g, e) {
	*DofElementData(e_H1, e->index) =
		phgQuadDofDotDof(e, vortex, vortex, QUAD_DEFAULT);
    }
#endif

    return;
}

static int
my_bc_map(int bctype)
{
    switch (bctype) {
#if 1
    case 0:
	return DIRICHLET;	/* Dirichlet */
    case 1:
	return NEUMANN;        /* Neumann */
    case 2:
	return BDRY_USER0;
    case 3:
	return BDRY_USER1;
    case 4:
	return BDRY_USER2;
    case 5:
	return BDRY_USER3;
    case 6:
	return BDRY_USER4;
    case 7:
	return BDRY_USER5;
    case 8:
	return BDRY_USER6;
    case 9:
	return BDRY_USER7;
    case 10:
	return BDRY_USER8;
    case 11:
	return BDRY_USER9;
    default:
	return NEUMANN;
#elif 1
#warning All DIRICH boundary	
    default:
	return DIRICHLET;
#else
#warning All NEUMANN boundary	
    default:
	return NEUMANN;
#endif
    }
}

static void
checkBdry(GRID *g)
{
    SIMPLEX *e;
    int s;
    double a[5];

    a[0] = a[1] = a[2] = a[3] = a[4] = 0.;
    ForAllElements(g, e) {
	for (s = 0; s < NFace; s++) {
	    FLOAT area = phgGeomGetFaceArea(g, e, s);
	    if (e->bound_type[s] & BDRY_MASK) {
		if (e->bound_type[s] & DIRICHLET)
		    a[0] += area;
		if (e->bound_type[s] & NEUMANN)
		    a[1] += area;
		if (e->bound_type[s] & BDRY_USER1)
		    a[2] += area;
		if (e->bound_type[s] & BDRY_USER2)
		    a[3] += area;
	    }
	    else {
		a[4] += area;
	    }
	}
    }

#if USE_MPI
    {
	double b[5];
	MPI_Reduce(a, b, 5, MPI_DOUBLE, MPI_SUM, 0, g->comm);
	memcpy(a, b, sizeof(b));
    }
#endif

    phgPrintf("Boundary types check:\n");
    phgPrintf("    Dirichlet: %g, Neumann: %g, user1: %g, user2: %g, "
	      "other: %g\n", a[0], a[1], a[2], a[3], a[4]);

    return;
}

int
main(int argc, char *argv[])
{
    GRID *g;
    MAT *matNS, *matB, *matC, *matAp, *matQp;
    DOF_TYPE *utype, *ptype;
    DOF *u0, *u, *gradu, *p, *f, *pbc, *gxbc, *gybc, *gzbc, *divu, *vortex,
	*dofs[64], *du, *dp, *u_exact, *p_exact, *gradu_exact, *eu, *ep,
	*egradu, *e_H1;
    SOLVER *solver, *pc = NULL;
    MAP *Vmap, *V1map, *Pmap, *Pmap0;
    MAT *pmat[4], *mats[16];
    FLOAT res, tol = 1e-3;
    INT step = 0, nstep = 50;
    BOOLEAN use_pc = TRUE, use_Fu = TRUE;
    size_t mem, mem_peak;
    INT mem_max = 400;
    char *utype_name = "P3";
    char *ptype_name = "P2";
    char *F_opts = NULL;
    char *Ap_opts = NULL;
    char *Qp_opts = NULL;

#if 0
    phgOptionsPreset("-solver gmres -gmres_restart 0 "
	 "-solver_maxit 500 -solver_rtol 1e-2");
#else
    phgOptionsPreset("-solver petsc -solver_maxit 50 -solver_rtol 1e-2 "
		     "-oem_options \"-ksp_type gmres -ksp_gmres_restart 50\" "
		     /*"-oem_options \"-ksp_monitor\""*/);
#endif

    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterString("utype", "DOF type for velocity", &utype_name);
    phgOptionsRegisterString("ptype", "DOF type for pressure", &ptype_name);
    phgOptionsRegisterFloat("Re", "Reynolds number", &Re);
    phgOptionsRegisterFloat("tol", "Tolerance", &tol);
    phgOptionsRegisterInt("mem_max", "Max memory per process", &mem_max);
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("nstep", "Max newton iter step", &nstep);
    phgOptionsRegisterNoArg("use_pc", "Use preconditioner", &use_pc);
    phgOptionsRegisterNoArg("use_Fu", "Use Fu in the preconditioner", &use_Fu);
    phgOptionsRegisterString("F_opts", "Solver F options", &F_opts);
    phgOptionsRegisterString("Ap_opts", "Solver Ap options", &Ap_opts);
    phgOptionsRegisterString("Qp_opts", "Solver Qp options", &Qp_opts);

    phgInit(&argc, &argv);
    phgOptionsShowUsed();

    /* set utype and ptype */
    {
	char s[128];
	phgOptionsPush();
	sprintf(s, "-dof_type %s", utype_name);
	phgOptionsSetOptions(s);
	utype = DOF_DEFAULT;
	sprintf(s, "-dof_type %s", ptype_name);
	phgOptionsSetOptions(s);
	ptype = DOF_DEFAULT;
	phgOptionsPop();
	assert(utype->fe_space == FE_H1 && ptype->fe_space == FE_H1
		&& utype->order > ptype->order);
    }

    g = phgNewGrid(-1);
    phgImportSetBdryMapFunc(my_bc_map);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    checkBdry(g);

    phgRefineAllElements(g, pre_refines);

    phgPrintf("Reynolds number = %lg\n", (double)(Re));
    nu = 1. / Re;

    u = phgDofNew(g, utype, 3, "u", DofInterpolation);
    u0 = phgDofNew(g, utype, 1, "u0", DofInterpolation);
    p = phgDofNew(g, ptype, 1, "p", DofInterpolation);
    e_H1 = phgDofNew(g, DOF_P0, 1, "estimator", DofNoAction);

    du = phgDofNew(g, utype, 3, "du", DofNoAction);
    dp = phgDofNew(g, ptype, 1, "dp", DofNoAction);
#if 0
    phgDofSetDataByFunction(u, func_u);
    phgDofSetDataByFunction(p, func_p);
#endif /* Exact solution. */

#ifdef ENCLOSED_FLOW
    phgPrintf("Enclosed flow\n");
#endif /* ENCLOSED_FLOW */

    pbc = phgDofNew(g, ptype, 1, "pbc", func_p);
    f = phgDofNew(g, DOF_ANALYTIC, 3, "f", func_f);
    gxbc = phgDofNew(g, DOF_ANALYTIC, 3, "gxbc", func_gx);
    gybc = phgDofNew(g, DOF_ANALYTIC, 3, "gybc", func_gy);
    gzbc = phgDofNew(g, DOF_ANALYTIC, 3, "gzbc", func_gz);
    u_exact = phgDofNew(g, DOF_ANALYTIC, 3, "exact u", func_u);
    p_exact = phgDofNew(g, DOF_ANALYTIC, 1, "exact p", func_p);
    gradu_exact = phgDofNew(g, DOF_ANALYTIC, 9, "exact u", func_gradu);

    /********************************************/
    /* Set velocity bdry mask: scalar or vector */
    /********************************************/
#if 1
    {
	BTYPE DB_masks[100];
# if 1
	/* all D&N, exact u, p, gradu */
	DB_masks[0] = BDRY_USER1|BDRY_USER3;
	DB_masks[1] = BDRY_USER4|NEUMANN;
	DB_masks[2] = DIRICHLET|BDRY_USER1|BDRY_USER2|NEUMANN;
# elif 1
	/* u[1] only N, exact u+C, p, gradu
	 *   singular does not matter. */
	DB_masks[0] = BDRY_USER1|BDRY_USER3;
	DB_masks[1] = 0; 
	DB_masks[2] = DIRICHLET|BDRY_USER1|BDRY_USER2|NEUMANN;
# elif 0
	/* u all D, exact u, p+C, gradu
	 *   singualar polute u, use mat check */
	DB_masks[0] = BDRY_MASK; 
	DB_masks[1] = BDRY_MASK; 
	DB_masks[2] = BDRY_MASK; 
# endif	 /* DB_masks combination */
	phgDofSetDirichletBoundaryMasks(u, DB_masks);
    }
#else
    u->DB_mask = BDRY_USER4|BDRY_USER1;
#endif

    /* Prepssure has no Dirichlet bdry condition. */
    p->DB_mask = 0; //DIRICHLET;

    /* Time advancing. */
    gradu = NULL;
    divu = NULL;
    while (TRUE) {
	static BOOLEAN initialized = FALSE;

	/* set Dirichlet BC */
	phgDofSetBdryDataByFunction(u, func_u, DIRICHLET);
    	phgDofSetBdryDataByFunction(p, func_p, 0);
	phgExportVTK(g, "test_DBs_init.vtk", u, p, NULL);

	/* (re)initialize grad(u) and div(u) */
	phgDofGradient(u, &gradu, NULL, "grad u");
	phgDofSetFunction(gradu, DofInterpolation);
	phgDofDivergence(u, &divu, NULL, "div u");
	phgDofSetFunction(divu, DofInterpolation);

	phgPrintf("\n%d DOF, %d elements, %d submesh%s, load imbalance: %lg\n",
	     DofGetDataCountGlobal(u), g->nleaf_global, g->nprocs,
	     g->nprocs > 1 ? "es" : "", (double)g->lif);

	if (phgBalanceGrid(g, 1.2, -1, NULL, 0.))
	    phgPrintf("Repartition mesh, %d submesh%s, load imbalance: %lg\n",
		      g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);

	if (!initialized) {
	    /* reset mem_peak */
	    phgMemoryUsageReset();
	    initialized = TRUE;
	}

	Vmap = phgMapCreate(u, NULL);
	V1map = phgMapCreate(u0, NULL);
	Pmap = phgMapCreate(p, NULL);

#ifdef ENCLOSED_FLOW
	pbc->DB_mask = 0;
#else
	pbc->DB_mask = DIRICHLET;
#endif
	Pmap0 = phgMapCreate(pbc, NULL);

	/* Newton nonliner iteration. */
	step = 0;
	while (TRUE) {
	    elapsed_time(g, FALSE, 0.);	/* reset timer */

	    /* matrices */
	    matF = phgMapCreateMat(Vmap, Vmap);
	    matFu = phgMapCreateMat(V1map, V1map);
	    matB = phgMapCreateMat(Pmap, Vmap);
	    matBt = phgMapCreateMat(Vmap, Pmap);
	    matC = /*phgMapCreateMat(Pmap, Pmap)*/NULL;

	    matFp = phgMapCreateMat(Pmap0, Pmap0);
	    matAp = phgMapCreateMat(Pmap0, Pmap0);
	    matQp = phgMapCreateMat(Pmap, Pmap);

	    matF->handle_bdry_eqns = TRUE;
	    matFu->handle_bdry_eqns = TRUE;
	    if (matC != NULL)
		matC->handle_bdry_eqns = TRUE;
	    matB->handle_bdry_eqns = FALSE;
	    matBt->handle_bdry_eqns = FALSE;

	    matAp->handle_bdry_eqns = FALSE;
	    matQp->handle_bdry_eqns = FALSE;
	    matFp->handle_bdry_eqns = FALSE;

	    pmat[0] = matF;
	    pmat[1] = matBt;
	    pmat[2] = matB;
	    pmat[3] = matC;

	    matNS = phgMatCreateBlockMatrix(g, 2, 2, pmat, NULL, NULL);
	    /* Note: can't use phgMat2Solver here because build_rhs()
	     * requires solver->rhs->map to map (u,p) into vector indices */
	    solver = phgSolverCreate(SOLVER_DEFAULT, u, p, NULL);
	    phgMatDestroy(&solver->mat);
	    solver->mat = matNS;
	    solver->rhs->mat = solver->mat;

	    /* Build matrices */
	    packDof(dofs, 9, u, p, gradu, divu, f, pbc, gxbc, gybc, gzbc);	/* gradu,divu changed. */
	    packMat(mats, 4, matB, matC, matAp, matQp);

	    build_rhs(solver, pc, dofs, mats);

	    res = phgVecNorm2(solver->rhs, 0, NULL);
	    mem = phgMemoryUsage(g, &mem_peak);
	    phgPrintf("    Newton step %2d, residual = %le, mem = %0.4lgMB "
			"(%0.4lgMB peak)\n", step, (double)res,
			(double)mem / (1024.0 * 1024.0),
			(double)mem_peak / (1024.0 * 1024.0));

	    phgPrintf("\tBuild matrices: ");
	    if (step < 1)
		build_matrices(solver, pc, dofs, mats, PICARD);
	    else
		build_matrices(solver, pc, dofs, mats, NEWTON);
	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

#if 0
#warning REDUNDENT MATCHECK
	    if (TRUE) {
		FILE *fdofs;
		VEC *x = phgMapCreateVec(solver->rhs->map, 1);
		DOF *dofs[2] = {u, p};

		fdofs = fopen("dofs_.m", "w");
		fprintf(fdofs, "m = %d; n=% d; m3 = %d; nu = %lg;",
			DofGetDataCount(u) / 3,
			DofGetDataCount(p), DofGetDataCount(u), nu);

		phgSolverUpdateRHS(solver);
		phgDofSetDataByFunction(u, func_u);
		phgDofSetDataByFunction(p, func_p);
		phgMapDofToLocalData(solver->rhs->map, 2, dofs, x->data);
		phgVecDumpMATLAB(solver->rhs, "b", "b_.m");
		phgVecDumpMATLAB(x, "x", "x_.m");

		phgMatDumpMATLAB(matF, "F", "F_.m");
		phgMatDumpMATLAB(matFu, "Fu", "Fu_.m");
		phgMatDumpMATLAB(matBt, "Bt", "Bt_.m");
		phgMatDumpMATLAB(matB, "B", "B_.m");
		phgVecDumpMATLAB(solver->rhs, "b", "rhs_.m");

		phgMatVec(MAT_OP_N, -1.0, solver->mat, x, 1.0, &solver->rhs);
		phgPrintf("matcheck, res:%e\n", phgVecNorm2(solver->rhs, 0, NULL));
		phgFinalize();
		exit(__LINE__);
	    }
#endif

	    if (use_pc) {
		/* set up preconditioner (the OEM_SOLVER is irrelavant here) */
		pc = phgMat2Solver(SOLVER_GMRES, solver->mat);
		phgSolverSetPC(solver, pc, pc_proc);

		/* solver F */
		phgOptionsPush();
		phgOptionsSetOptions("-solver hypre "
				     "-hypre_solver gmres "
				     "-hypre_pc boomeramg "
				     "-solver_maxit 10 "
				     "-solver_rtol 1e-3");
		phgOptionsSetOptions(F_opts);
		if (use_Fu) {
		    /* use diag(matFu, matFu, matFu) in the precond. matrix */
		    solver_F = phgMat2Solver(SOLVER_DEFAULT, matFu);
		}
		else {
		    /* use matF in the preconditioning matrix */
		    solver_F = phgMat2Solver(SOLVER_DEFAULT, matF);
		}
		solver_F->warn_maxit = FALSE;
		phgOptionsPop();

		/* solver Ap */
		phgOptionsPush();
		phgOptionsSetOptions("-solver hypre "
				     "-hypre_solver gmres "
				     "-hypre_pc boomeramg "
				     "-solver_maxit 10 "
				     "-solver_rtol 1e-2");
		phgOptionsSetOptions(Ap_opts);
		solver_Ap = phgMat2Solver(SOLVER_DEFAULT, matAp);
		solver_Ap->warn_maxit = FALSE;
		phgOptionsPop();

		/* solver Qp */
		phgOptionsPush();
		phgOptionsSetOptions("-solver pcg ");
		phgOptionsSetOptions("-solver_maxit 10 -solver_rtol 1e-2");
		phgOptionsSetOptions(Qp_opts);
		solver_Qp = phgMat2Solver(SOLVER_DEFAULT, matQp);
		solver_Qp->warn_maxit = FALSE;
		phgOptionsPop();
	    }

	    phgPrintf("\tAssemble linear system: ");
	    phgSolverAssemble(solver);
	    if (use_pc) {
		phgSolverAssemble(solver_F);
		phgSolverAssemble(solver_Ap);
		phgSolverAssemble(solver_Qp);
	    }
	    elapsed_time(g, TRUE, 0.);

	    phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	    phgDofSetDataByValue(du, 0.);
	    phgDofSetDataByValue(dp, 0.);
	    phgPrintf("\tSolution: ");
	    phgSolverSolve(solver, TRUE, u, p, NULL);
	    phgPrintf("nits = %d, resid = %0.4le ",
		      solver->nits, (double)solver->residual);
	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	    phgSolverDestroy(&solver_F);
	    phgSolverDestroy(&solver_Ap);
	    phgSolverDestroy(&solver_Qp);

	    phgMatDestroy(&matF);
	    phgMatDestroy(&matFu);
	    phgMatDestroy(&matB);
	    phgMatDestroy(&matBt);
	    phgMatDestroy(&matC);
	    phgMatDestroy(&matAp);
	    phgMatDestroy(&matQp);
	    phgMatDestroy(&matFp);

	    phgSolverDestroy(&solver);
	    phgSolverDestroy(&pc);

	    //phgDofAXPY(1.0, du, &u);
	    //phgDofAXPY(1.0, dp, &p);

	    phgDofGradient(u, &gradu, NULL, "grad u");
	    phgDofDivergence(u, &divu, NULL, "div u");
	    break;
	}			/* end of nonlinear iteration */

	vortex = phgDofCurl(u, NULL, NULL, "vorticity");

#if OUTPUT_VTK
	if (FALSE) {
	    char vtkfile[100];
	    DOF *u_exact2, *p_exact2;
	    u_exact2 = phgDofNew(g, utype, 3, "u", func_u);
	    p_exact2 = phgDofNew(g, ptype, 1, "p", func_p);
	    sprintf(vtkfile, "ins_%03d.vtk", step);
	    phgExportVTK(g, vtkfile, u, p, vortex, u_exact2, p_exact2, eu,
			     ep, NULL);
	    sprintf(vtkfile, "ins_%03d.dx", step);
	    phgExportDX(g, vtkfile, u, p, vortex, u_exact2, p_exact2, eu,
			    ep, NULL);
	    phgDofFree(&u_exact2);
	    phgDofFree(&p_exact2);
	}
#endif /* OUTPUT_VTK */

	phgMapDestroy(&Vmap);
	phgMapDestroy(&V1map);
	phgMapDestroy(&Pmap);
	phgMapDestroy(&Pmap0);

	/* ----Error check------- */
	phgPrintf("    Errors: ");
	eu = phgDofCopy(u, NULL, utype, "err u");
	ep = phgDofCopy(p, NULL, ptype, "err p");
	egradu = phgDofCopy(gradu, NULL, utype, "err grad u");

	phgDofAXPY(-1.0, u_exact, &eu);
	phgDofAXPY(-1.0, gradu_exact, &egradu);
	phgDofAXPY(-1.0, p_exact, &ep);
	phgPrintf("L2 = (%0.6le,%0.6le)\n            H1 = (%0.6le), "
		  "div(u) = %0.6le ", phgDofNormL2(eu), phgDofNormL2(ep),
		  res = phgDofNormL2(egradu), phgDofNormL2(divu));
	elapsed_time(g, TRUE, 0.);

	phgDofFree(&egradu);
	phgDofFree(&eu);
	phgDofFree(&ep);

	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("    Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		      (double)mem / (1024.0 * 1024.0),
		      (double)mem_peak / (1024.0 * 1024.0));

	if (TRUE || res < tol || mem_peak > 1024 * 1024 * (size_t)mem_max)
	    break;

#if 0
	phgRefineAllElements(g, 1);
#else
    {
	SIMPLEX *e;
	estimate_error(u, p, f, gradu, divu, vortex, e_H1);
	phgMarkRefine(MARK_GERS, e_H1, Pow(0.5, 2), NULL, 0.0, 1, 0.0);
	/* mark boundary elements */
	ForAllElements(g, e) {
	    int i;
	    for (i = 0; i < NFace; i++)
		if (e->bound_type[i] & BDRY_MASK)
		    break;
	    if (i < NFace)
		e->mark = 1;
	}
	phgRefineMarkedElements(g);
    }
#endif
	phgDofFree(&vortex);
    }				/* end of time advaning */

    phgDofFree(&gradu);
    phgDofFree(&divu);

#if OUTPUT_VTK
    phgExportVTK(g, "ins.vtk", u, p, vortex, NULL);
    phgExportDX(g, "ins.dx", u, p, vortex, NULL);
    phgDofFree(&vortex);
#endif /* OUTPUT_VTK */

    phgDofFree(&e_H1);
    phgDofFree(&u);
    phgDofFree(&u0);
    phgDofFree(&p);
    phgDofFree(&f);
    phgDofFree(&u_exact);
    phgDofFree(&gradu_exact);
    phgDofFree(&p_exact);
    phgDofFree(&pbc);
    phgDofFree(&gxbc);
    phgDofFree(&gybc);
    phgDofFree(&gzbc);
    phgDofFree(&du);
    phgDofFree(&dp);

    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
