#include "phg.h"
#define _p new_params
#define _nsp (ns->ns_params)

/*
 * Benchmark experiments for higher-order and full-Stokes ice sheet
 * models (ISMIP–HOM) F. Pattyn, L. Perichon, A. Aschwanden,
 * B. Breuer, B. de Smedt, O. Gagliardini, G. H. Gudmundsson,
 * R. C. A. Hindmarsh, A. Hubbard, J. V. Johnson, T. Kleiner,
 * Y. Konovalov, C. Martin, A. J. Payne, D. Pollard, S. Price,
 * M. Rückamp, F. Saito, O. Souček, S. Sugiyama, and
 * T. Zwinger
 *
 *
 * Experiment. A,B,C,D,F
 *   Domain [0, L] x [0, L]
 *   L = 5, 10, 20, 40, 80, 160km 
 *   Height is different for each experiment.
 * Experiment E
 *   Haut Glacier d'Arolla, 2D grid(xz).
 * 
 *
 *  */

/* ------------------------------------------------------------
 *    
 *    Ice grid
 *    
 * ------------------------------------------------------------ */
void 
iceInit(GRID *g, LAYERED_MESH **gL)
{
    SIMPLEX *e;

    /* Build ice grid with given func_ice_slab */
    ice_grid(g);
    checkBdry(g);

    /* Read the partition info from region_mark. */
    ForAllElements(g, e) {
	e->mark = e->region_mark / 10;
	e->region_mark = 0;
    }

    phgPartUserSetFunc(iceParter);
    phgRedistributeGrid(g);
    phgPrintf("\nRepartition mesh, %d submesh%s, load imbalance: %lg",
		  g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

    if (1) phgExportVTK(g, "parted.vtk", NULL);

}






#if ICE_BENCH_TEST
void
func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord)
{
    const FLOAT alpha = _alpha_ / 180. * M_PI;
    const FLOAT L = _Length_;
#  if TEST_CASE == ICE_BENCH_A
    FLOAT Ha = -x * tan(alpha);
    FLOAT Hb = Ha - (1./L) + (.5/L) * sin(2*M_PI * x) * sin(2*M_PI * y);
#  elif TEST_CASE == ICE_BENCH_B
    FLOAT Ha = -x * tan(alpha);
    FLOAT Hb = Ha - (1./L) + (.5/L) * sin(2*M_PI * x);
#  elif TEST_CASE == ICE_BENCH_C
    FLOAT Ha = -x * tan(alpha);
    FLOAT Hb = Ha - (1./L); 
#  elif TEST_CASE == ICE_BENCH_D
    FLOAT Ha = -x * tan(alpha);
    FLOAT Hb = Ha - (1./L); 
#  elif TEST_CASE == ICE_BENCH_F
    FLOAT Ha = 0;
    FLOAT Hb = Ha - (1./L) + (0.1/L) * exp( - (SQUARE(x - .5) + SQUARE(y - .5)) / SQUARE(10./L) ); 
#  else
    phgError("Wrong ice bench case!!!\n");
#  endif
    FLOAT Hz = Ha - Hb;

    FLOAT rz = z;		/* ratio z in [0, 1] */
    z = Hb + Hz * rz;

    coord[0] = x * L * 1000 / LEN_SCALING;
    coord[1] = y * L * 1000 / LEN_SCALING;
    coord[2] = z * L * 1000 / LEN_SCALING;
#  if TEST_CASE == ICE_BENCH_F
    x = coord[0];
    z = coord[2];
    FLOAT c = cos(alpha), s = sin(alpha);
    coord[0] =  x * c + z * s;
    coord[2] = -x * s + z * c;
#  endif
}

#elif TEST_CASE == ICE_BENCH_E
void
func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord)
{
    coord[0] = x / LEN_SCALING;
    coord[1] = y / LEN_SCALING;
    coord[2] = z / LEN_SCALING;
}

#else
#  error Test case wrong!!!
#endif




/* ------------------------------------------------------------
 *    
 *    Effective viscosity
 *    
 * ------------------------------------------------------------ */

FLOAT 
get_effective_viscosity(const FLOAT *gu, FLOAT T, FLOAT p,
			int viscosity_type)
{
    FLOAT A = 1e-16;
    const FLOAT n = POWER_N;
    const FLOAT a = SEC_PER_YEAR;
    const FLOAT L = _Length_;

    FLOAT yeta, eps;
    FLOAT eu[Dim*Dim];
    int k;

    Unused(a);
    T = 0;			/* Temperature unused */
    assert(p == 0.);		/* Unused */
#  if 0
#  warning const vis -------------------------------------------------
    nu_max = nu_min = 1e7;
    return 1e7;
#  endif

#if TEST_CASE == ICE_BENCH_F
    A = 2.140373e-7;
    yeta = 1. / A;
    nu_max = nu_min = yeta;
    return yeta;
#endif

    if (viscosity_type == VIS_CONST) {
	eps = (100.) / (L * 1000.);	/* initial guess */
    } else if (viscosity_type == VIS_STRAIN) {
	MAT3_SYM(gu, eu);
	for (k = 0; k < DDim; k++)
	    eu[k] /= LEN_SCALING;
	eps = sqrt(.5) * MAT3_NORM2(eu);
    } else {
	phgError(1, "Unknown viscosity type\n!!!");
    }

    if (eps < MIN_EFFECTIVE_STRAIN) 
	eps = MIN_EFFECTIVE_STRAIN;

    yeta = pow(A, -1./n) * pow(eps, (1.-n)/n);

    nu_max = MAX(nu_max, yeta);
    nu_min = MIN(nu_min, yeta);

    return yeta;
}



/* ------------------------------------------------------------
 *    
 *    Fraction coef
 *    
 * ------------------------------------------------------------ */


void func_beta(FLOAT x, FLOAT y, FLOAT z, FLOAT *beta) {
    const FLOAT L = _Length_;
    Unused(L);

#if TEST_CASE == ICE_BENCH_C
    beta[0] = 1000. + 1000. * sin(2*M_PI * x / L / 1000 * LEN_SCALING) *
    	sin(2*M_PI * y / L / 1000 * LEN_SCALING);
#elif TEST_CASE == ICE_BENCH_D
    beta[0] = 1000. + 1000. * sin(2*M_PI * x / L / 1000 * LEN_SCALING);
#elif TEST_CASE == ICE_BENCH_E
    beta[0] = 0.;
#elif TEST_CASE == ICE_BENCH_F
#  if USE_SLIDING_BC
    FLOAT A = 2.140373e-7;
    FLOAT H0 = 1e3;
    beta[0] = 1. / (1. * A * H0);
#  else
    beta[0] = 0.;
#  endif
#endif
}




/* ------------------------------------------------------------
 *    
 *    B.C. funcs
 *    
 * ------------------------------------------------------------ */


void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *u) {
    u[0] = 0;
    u[1] = 0;
    u[2] = 0;
}

void func_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *p) {
    p[0] = 0;
}

void func_gradu(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradu) {
    gradu[0] = 0;
    gradu[1] = 0;
    gradu[2] = 0;
    gradu[3] = 0;
    gradu[4] = 0;
    gradu[5] = 0;
    gradu[6] = 0;
    gradu[7] = 0;
    gradu[8] = 0;
}

void func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *f) {
    f[0] = 0;
    f[1] = 0;
    f[2] = 0;
}

void func_g1(FLOAT x, FLOAT y, FLOAT z, FLOAT *g1) {
    g1[0] = 0;
    g1[1] = 0;
    g1[2] = 0;
}

void func_g2(FLOAT x, FLOAT y, FLOAT z, FLOAT *g2) {
    g2[0] = 0;
    g2[1] = 0;
    g2[2] = 0;
}

void func_g3(FLOAT x, FLOAT y, FLOAT z, FLOAT *g3) {
    g3[0] = 0;
    g3[1] = 0;
    g3[2] = 0;
}

void func_T(FLOAT x, FLOAT y, FLOAT z, FLOAT *T) {
    T[0] = 0;
}

void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *fT) {
    fT[0] = 0;
}

void func_q(FLOAT x, FLOAT y, FLOAT z, FLOAT *q) 
{
    *q = 0.;
}

/* ------------------------------------------------------------
 *    
 *    b.C. map
 *    
 * ------------------------------------------------------------ */

int
my_bc_map(int bctype)
{
#if TEST_CASE == ICE_BENCH_E
    switch (bctype) {
    case 0:
	return DIRICHLET;
    case 1:
	return BC_TOP;
    case 2:
	return BC_BOTTOM;
    case 3:
#  if ZERO_TRACT_ZONE
	return BC_BOTTOM2;
#  else
	return BC_BOTTOM;
#  endif
    default:
	return DIRICHLET;;
    }
#else
    return DIRICHLET;    
#endif
}



void 
set_boundary_mask(NSSolver *ns)    
{
    DOF **u = ns->u;
    DOF **p = ns->p;
    DOF **T = ns->T;

#if TEST_CASE == ICE_BENCH_C			\
    || TEST_CASE == ICE_BENCH_D
    BTYPE DB_masks[3] = {INFLOW, 0, 0,};

#elif TEST_CASE == ICE_BENCH_A			\
    || TEST_CASE == ICE_BENCH_B
    BTYPE DB_masks[3] = {BC_BOTTOM,
    			 BC_BOTTOM,
    			 BC_BOTTOM};
    /* BTYPE DB_masks[3] = {BC_BOTTOM | BC_LATERL, */
    /* 			 BC_BOTTOM | BC_LATERL, */
    /* 			 BC_BOTTOM | BC_LATERL}; */

#elif TEST_CASE == ICE_BENCH_E
    BTYPE DB_masks[3] = {BC_BOTTOM | BC_BOTTOM2, 
			 BC_BOTTOM, 
			 BC_BOTTOM};

#elif TEST_CASE == ICE_BENCH_F 
#  if USE_SLIDING_BC
#  warning sliping bdry
    BTYPE DB_masks[3] = {INFLOW, 0, 0};
#  else
#  warning non-sliping bdry
    BTYPE DB_masks[3] = {INFLOW, INFLOW, INFLOW};
#  endif
#endif

    phgDofSetDirichletBoundaryMasks(u[1], DB_masks);
    if (_nsp->enclosed_flow) {
	if (_nsp->pin_node)
	    phgDofSetDirichletBoundaryMask(p[1], PINNEDNODE);
	else
	    phgDofSetDirichletBoundaryMask(p[1], 0);
    } else {
	phgDofSetDirichletBoundaryMask(p[1], 0);
    }
    phgDofSetDirichletBoundaryMask(T[1], 0);

    return;
}


void iceSetBoundaryTypes(NSSolver *ns)
{
    /* Unused. */
}


void func_a(FLOAT x, FLOAT y, FLOAT *q)
{
}

void func_s(FLOAT x, FLOAT y, FLOAT z, FLOAT *q)
{
}

void func_b(FLOAT x, FLOAT y, FLOAT z, FLOAT *q)
{
}
