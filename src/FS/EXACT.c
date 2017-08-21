#include "phg.h"
#define _p new_params
#define _nsp (ns->ns_params)


/* ------------------------------------------------------------
 *    
 *    Ice grid
 *    
 * ------------------------------------------------------------ */
void 
iceInit(GRID *g, LAYERED_MESH **gL_ptr)
{
    SIMPLEX *e;

    /* Init func for exact-solutions.c */
    func_init_params(_Length_, _alpha_ * M_PI / 180.);

    /* Build ice grid with given func_ice_slab */
    ice_grid(g);
    checkBdry(g);

    /* Struct mesh init */
    struct_mesh_init(g);

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
    *gL_ptr = NULL;
}






void
func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord)
/* input: [0,1]x[0,1]x[0,1], or
 *        [m]/L [m]/L [m]/Z */
{
    const FLOAT alpha = _alpha_ / 180. * M_PI;
    const FLOAT L = _Length_;

    /* unit: m */
    x *= L*1000;
    y *= L*1000;

    static FLOAT Ha, Hb;
    func_s(x, y, 0, &Ha);	
    func_b(x, y, 0, &Hb);

    FLOAT Hz = Ha - Hb;
    FLOAT rz = z;		/* ratio z in [0, 1] */
    z = Hb + Hz * rz;

    /* output: [m]/L */
    coord[0] = x / LEN_SCALING;
    coord[1] = y / LEN_SCALING;
    coord[2] = z / LEN_SCALING;
}

static FLOAT
func_ice_slabH(FLOAT x, FLOAT y)
{
    const FLOAT alpha = _alpha_ / 180. * M_PI;
    const FLOAT L = _Length_;

    /* Unit: m */
    x *= LEN_SCALING; 
    y *= LEN_SCALING;

    static FLOAT Ha, Hb;
    FLOAT z = 0.;
    /* Unit: m */
    func_s(x, y, 0, &Ha);	
    func_b(x, y, 0, &Hb);

    /* unit: m/L */
    return (Ha - Hb) / LEN_SCALING;
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

    FLOAT yeta, eps = 0.;
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
    /* Unused */
    *beta = 1e10;
}






/* ------------------------------------------------------------
 *    
 *    b.C. map
 *    
 * ------------------------------------------------------------ */

int
my_bc_map(int bctype)
{
    return DIRICHLET;    	/* clear it! */
}


void 
set_boundary_mask(NSSolver *ns)    
{
    DOF **u = ns->u;
    DOF **p = ns->p;
    DOF **T = ns->T;

    BTYPE DB_masks[3] = {BC_BOTTOM,
    			 BC_BOTTOM,
    			 BC_BOTTOM};
    /* BTYPE DB_masks[3] = {BC_BOTTOM | BC_LATERL, */
    /* 			 BC_BOTTOM | BC_LATERL, */
    /* 			 BC_BOTTOM | BC_LATERL}; */

/*     BTYPE DB_masks[3] = {BC_BOTTOM | BC_TOP, */
/*     			 BC_BOTTOM | BC_TOP, */
/*     			 BC_BOTTOM | BC_TOP}; */


    phgDofSetDirichletBoundaryMasks(u[1], DB_masks);
    phgDofSetDirichletBoundaryMask(T[1], 0);

    return;
}


void iceSetBoundaryTypes(NSSolver *ns)
{
    /* Unused. */
}



/* ------------------------------------------------------------
 *    
 *    B.C. funcs
 *    
 * ------------------------------------------------------------ */
#define t Time
#include "exact-solutions.c"
