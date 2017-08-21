/*
 * ISMIP-HEINO benchmark:
 *
 * 2D region [0: 1,500km] x [0: 1,500km].
 * Height is zero at beginning, to void degenerated elements,
 *   use a eps height instead.
 *
 * Ref:
 * Calov R, Greve R, Abe-Ouchi A, Bueler E, Huybrechts P, Johnson JV,
 * Pattyn F, Pollard D, Ritz C, Saito F, Tarasov L (2010) Results from
 * the ice sheet model intercomparison project - Heinrich event
 * intercomparison (ISMIP HEINO). J Glaciol 56 (197): 371-383 
 * 
 * */


/* ------------------------------------------------------------
 *    
 *    Ice grid
 *    
 * ------------------------------------------------------------ */

void
func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord)
{
    FLOAT Ha, Hb, Hz;
    FLOAT rz = z;		/* ratio z in [0, 1] */

    x *= 1.5;
    y *= 1.5;
    Hb = 0.;
    Hz = eps_height;

    coord[0] = x;
    coord[1] = y;
    coord[2] = Hb + rz * Hz;

    return;
}

FLOAT
func_ice_slabH(FLOAT x, FLOAT y)
{
    FLOAT Hz;

    Hz = eps_height;    
    return Hz;
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
    default:
	return DIRICHLET;
    }
}
