/*
 * -------------------- 
 * Test-solverT B.C.
 * --------------------
 *
 * 2D region [0: 1,500km] x [0: 1,500km].
 * Height is zero at beginning, to void degenerated elements,
 *   use a eps height instead.
 *
 * Ref:
 * 
 * [1] The EISMINT benchmarks for testing ice-sheet models,
 * Huybrechts, P. Payne, T., Annals, 1996, VOL 23, pages 1-12
 * 
 * [2] Results from the EISMINT model intercomparison: the effects of
 * thermomechanical coupling Payne, A. J. Huybrechts, P. Abe-Ouchi,
 * A. Calov, R. Fastook, J. L.; Greve, R. Marshall, S. J. Marsiat,
 * I. Ritz, C. Tarasov, L.  Journal of Glaciology, 2000, VOL 46; ISSU
 * 153, pages 227-238
 * 
 * */

/* ------------------------------------------------------------
 *    
 *    Ice grid
 *    
 * ------------------------------------------------------------ */
#warning Test solverT B.C.

#  warning Thickness_0 = H0
#define H0 1.5e3

void
func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord)
{
    FLOAT Ha, Hb, Hz;
    FLOAT rz = z;		/* ratio z in [0, 1] */
    FLOAT d;

    x *= 1.5;
    y *= 1.5;
    d = sqrt((x-.75)*(x-.75) + (y-.75)*(y-.75));
    Hb = 0.;

    Hz = H0 / LEN_SCALING; //HEIGHT_EPS;
    //Hz = 5e3 / LEN_SCALING * (.5 + exp(-d/5.)); //HEIGHT_EPS;

    coord[0] = x;
    coord[1] = y;
    coord[2] = Hb + rz * Hz;

    return;
}

/* ------------------------------------------------------------
 *    
 *    B.C. map
 *    
 * ------------------------------------------------------------ */

int
my_bc_map(int bctype)
{
    switch (bctype) {
    case 0:
	return BC_BOTTOM;
    case 1:
	return BC_TOP;
    default:
	return DIRICHLET;
    }
}

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

    FLOAT yeta, eps = 0;
    FLOAT eu[Dim*Dim];
    int k;

    Unused(a);
    assert(p == 0.);		/* unused */

#  if 0
#  warning Fixme: const vis -------------------------------------------------
    nu_max = MAX(nu_max, 1e9);
    nu_min = MIN(nu_min, 1e9);
    return 1e9;
#  endif

    const FLOAT T0 = ARRHENIUS_T;
    const FLOAT Q0 = ARRHENIUS_Q0;
    const FLOAT a0 = ARRHENIUS_A0;
    const FLOAT Q1 = ARRHENIUS_Q1;
    const FLOAT a1 = ARRHENIUS_A1;
    const FLOAT R  = GAS_CONSTANT;

    assert(T <= 300 && T > 150); /* reasonable temp region */

#if 0
    if (T < T0)
        A = a0 * exp(-Q0 / R / T);
    else
        A = a1 * exp(-Q1 / R / T);
#else
#warning Fixme: const vis A = 1e-16
    A = 1e-16;
#endif

    //#warning Fixme: const vis
    if (viscosity_type == VIS_CONST) {
	/* Initial guess:
	 * 1 m/a on 1km surf
	 * */
	eps = 1 / 1000; 
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
 *    B.C. funcs
 *    
 * ------------------------------------------------------------ */


void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *u) {
#warning set velocity to convection
    u[0] = -50.;
    u[1] = 30.;
    u[2] = -0.01;
}

void func_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *p) {
#if 0
#warning Fixme: debug pres
    const FLOAT rho = RHO_ICE;
    const FLOAT grav = GRAVITY;
    p[0] = - rho * grav * (z * LEN_SCALING - 5e3) / PRES_SCALING;
#else
    p[0] = 0.;
#endif
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
    const FLOAT T_min = 238.15;
    const FLOAT x0 = 750e3 / LEN_SCALING;
    const FLOAT y0 = 750e3 / LEN_SCALING;
    const FLOAT ST = 1.67e-2 / 1e3 * LEN_SCALING;

    FLOAT d;
    x -= x0;
    y -= y0;
    d = sqrt(x*x + y*y);
    *T = T_min + ST * d;
    //#warning Set Temp may above melt point
    *T += 0.02 * (H0 - z * LEN_SCALING);
    if (*T > TEMP_WATER)
	*T = TEMP_WATER;
}

void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *fT) {
    fT[0] = 0;
}

void func_q(FLOAT x, FLOAT y, FLOAT z, FLOAT *q) {
    const FLOAT M_max = .5;
    const FLOAT x0 = 750e3 / LEN_SCALING;
    const FLOAT y0 = 750e3 / LEN_SCALING;
    const FLOAT Rel = 450e3 / LEN_SCALING;
    const FLOAT Sb = 1e-2 / 1e3 * LEN_SCALING; 

    FLOAT d;
    x -= x0;
    y -= y0;
    d = sqrt(x*x + y*y);
    *q = Sb * (Rel - d);
    *q = MIN(M_max, *q);	/* unit: meter */
    return;
}


void 
set_boundary_mask(NSSolver *ns)    
{
    BTYPE DB_masks[3] = {BC_BOTTOM | BC_LATERL, 
			 BC_BOTTOM | BC_LATERL, 
			 BC_BOTTOM | BC_LATERL};
    phgDofSetDirichletBoundaryMasks(ns->u[1], DB_masks);


    //phgDofSetDirichletBoundaryMask(ns->T[1], BC_TOP | BC_LATERL);
    phgDofSetDirichletBoundaryMask(ns->T[1], BC_TOP);

    return;
}


void func_beta(FLOAT x, FLOAT y, FLOAT z, FLOAT *beta) {
    *beta = 0.;
}
