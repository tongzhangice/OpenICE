/*
 *
 * Greenland data from
 *   http://websrv.cs.umt.edu/isis/images/a/a5/Greenland_5km_v1.1.nc
 *   in netCDF format.
 *
 * Domain: [-800km, 700km] x [-3400km, -600km]
 *         1500km x 2800km
 * resolution: 5km
 * Scaled grid region: [0, 300] x [0, 560]
 *
 *
 *
 *  */

#define GRID_SCALING 5000	/* Note: input coord range
				 * [0, 301] x [0, 561]
				 * */


/* ------------------------------------------------------------
 *    
 *    Ice grid
 *    
 * ------------------------------------------------------------ */
void 
iceInit(GRID *g, LAYERED_MESH **gL_ptr)
{
    SIMPLEX *e;
    LAYERED_MESH *gL;

    /* Read data */
    read_nc_data(ns_params->nc_file);

    /* Build ice grid with given func_ice_slab */
    ice_grid(g);
    checkBdry(g);

    /* build layered mesh */
    gL = import_layered_mesh(ns_params->tria_file,
			     ns_params->layer_file,
			     ns_params->nodeZ_file,
			     NULL,
			     phgNProcs);
    build_layered_mesh(g, gL);
    part_layered_mesh(g, gL);	/* Partition saved in e->mark. */
    destory_layerd_mesh(&gL);
    phgPartUserSetFunc(iceParter);
    if (phgBalanceGrid(g, 1.1, -1, NULL, 0.)) {
	phgPrintf("\nRepartition mesh, %d submesh%s, load imbalance: %lg",
		  g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
    }
    if (0) {
	phgExportVTK(g, "parted.vtk", NULL);
    }

    gL = import_layered_mesh(ns_params->tria_file,
			     ns_params->layer_file,
			     ns_params->nodeZ_file,
			     NULL,
			     phgNProcs);
    build_layered_mesh(g, gL);
    if (1) {
	phgExportTecplot(g, "parted.plt", NULL);
	//phgExportVTK(g, OUTPUT_DIR "/ins_" NS_PROBLEM "_init.vtk", NULL);
    }


    *gL_ptr = gL;
    return;
}



void
func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord)
{
    FLOAT Ha, Hb, Hz;
    FLOAT rz = z;		/* ratio z in [0, 1] */

    x *= GRID_SCALING / LEN_SCALING;
    y *= GRID_SCALING / LEN_SCALING;

    nc_data_set_active(data_index_thk);
    nc_data_scaling = 1. / LEN_SCALING;
    interp_nc_data(x, y, z, &Hz);

    if (Hz < HEIGHT_EPS * LEN_SCALING) {
        phgInfo(0, "   Thickness At (%e %e) is too small: %e, "
                "     change to %e\n",
               x, y, Hz, HEIGHT_EPS*LEN_SCALING);
        Hz = HEIGHT_EPS*LEN_SCALING;
    }

    nc_data_set_active(data_index_topg);
    nc_data_scaling = 1. / LEN_SCALING;;
    interp_nc_data(x, y, z, &Hb);

    coord[0] = x;
    coord[1] = y;
    coord[2] = Hb + rz * Hz;

    return;
}

FLOAT
func_ice_slabH(FLOAT x, FLOAT y)
{
    FLOAT Hz, z = 0;
    x *= GRID_SCALING / LEN_SCALING;
    y *= GRID_SCALING / LEN_SCALING;

    nc_data_set_active(data_index_thk);
    nc_data_scaling = 1. / LEN_SCALING;;
    interp_nc_data(x, y, z, &Hz);
    
    return Hz;
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

#if 1
#  warning vis - theraml coupled
    if (T < T0)
        A = a0 * exp(-Q0 / R / T);
    else
        A = a1 * exp(-Q1 / R / T);
    A *= SEC_PER_YEAR;
#elif 0
#  warning const vis A(T = -10)
    T = TEMP_WATER - 10.;
    A = a1 * exp(-Q1 / R / T);
    A *= SEC_PER_YEAR;
#else
#  warning const vis A = 1e-16
    A = 1e-16;
#endif

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
 *    B.C. function
 *    
 * ------------------------------------------------------------ */

void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *u) {
#if 0
    /* u */
    nc_data_set_active(data_index_surfvelx);
    nc_data_scaling = 1. / LEN_SCALING;
    interp_nc_data(x, y, z, u);
    /* v */
    nc_data_set_active(data_index_surfvely);
    nc_data_scaling = 1. / LEN_SCALING;
    interp_nc_data(x, y, z, u+1);
    /* w */
    u[2] = 0;
#else
#  warning func_u = zero
    u[0] = u[1] = u[2] = 0.;
#endif
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
    /* Temp */
    x *= LEN_SCALING;
    y *= LEN_SCALING;

    nc_data_set_active(data_index_temp);
    nc_data_scaling = 1.;
    interp_nc_data_3D(x, y, z, 0, T);

    *T += TEMP_WATER;			
}

void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *fT) {
    fT[0] = 0;
}

void func_beta(FLOAT x, FLOAT y, FLOAT z, FLOAT *beta) 
{
    /* Read from 2D file */

    /* nc_data_set_active(data_index_beta); */
    /* nc_data_scaling = 1.; */
    /* interp_nc_data(x, y, 0, beta); */

    return;
}

void func_q(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
    *value = 0;
}

void func_a(FLOAT x, FLOAT y, FLOAT *value) {
    *value = 0;
}

void func_s(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
    *value = 0;
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
	return BC_LATERL;
    case 1:
	return BC_BOTTOM;
    case 2:
	return BC_TOP;
    default:
	return DIRICHLET;
    }
}

void 
set_boundary_mask(NSSolver *ns)    
{
#if !USE_SLIDING_BC
#   warning BD mask set to (x,x,x)
    BTYPE DB_masks[3] = {BC_BOTTOM,
			 BC_BOTTOM,
			 BC_BOTTOM};
#else
#   warning BD mask set to (x,0,0)
    BTYPE DB_masks[3] = {BC_BOTTOM,
			 0, 
			 0};
#endif

    phgDofSetDirichletBoundaryMasks(ns->u[1], DB_masks);
    phgDofSetDirichletBoundaryMask(ns->T[1], BC_TOP);

    return;
}


void iceSetBoundaryTypes(NSSolver *ns)
{
    /* Do nothing */
}

