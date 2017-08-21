/*
 * ESIMINT benchmark:
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
void 
iceInit(GRID *g, LAYERED_MESH **gL_ptr)
{
    /* build ice grid with given func_ice_slab */
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
    if (0) phgExportVTK(g, "parted.vtk", NULL);

    gL = import_layered_mesh(ns_params->tria_file,
			     ns_params->layer_file,
			     ns_params->nodeZ_file,
			     NULL,
			     phgNProcs);
    build_layered_mesh(g, gL);
    if (0) phgExportVTK(g, OUTPUT_DIR "/ins_" NS_PROBLEM "_init.vtk", NULL);


    /* init fv solver */
    fv_solver_init(ns_params->vert_file,
		   ns_params->tria_file,
		   func_Q);

    *gL_ptr = gL;
    return;
}

void
func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord)
{
    FLOAT Ha, Hb, Hz;
    FLOAT rz = z;		/* ratio z in [0, 1] */
    FLOAT d;

    x *= 1;
    y *= 1;
    d = sqrt((x-750.)*(x-750.) + (y-750.)*(y-750.));
    Hb = 0.;

#if 0
#  warning Thickness_0 = 1km
    Hz = 1e3 / LEN_SCALING; //HEIGHT_EPS;
    //Hz = 5e3 / LEN_SCALING * (.5 + exp(-d/5.)); //HEIGHT_EPS;
#else
#  warning Vialov set, h0 = 3.2km
    FLOAT xx = x - 750.;
    FLOAT yy = y - 750.;
    d = sqrt(xx*xx + yy*yy);

    FLOAT h0 = 3.2;
    FLOAT L = 579;
    FLOAT n = POWER_N;

    if (d < L)
	Hz = h0 * pow(fabs(1 - pow(d/L, (n+1)/n)),
		      n/(2*n+2));
    else
	Hz = HEIGHT_EPS;
#endif

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

void iceSetBoundaryTypes(NSSolver *ns)
{
    GRID *g = ns->g;
    
#if TEST_CASE == ESIMINT_H
#   define EPS_T 1e-6		/* Obviously bigger than melting point means melting */
    {
	SIMPLEX *e;
	DOF *dof_depth = NULL;
	int nfrozen = 0;

	assert(ns->T[1]->type == ns->u[1]->type);
	dof_depth = ns->depth_P2;
	
	FLOAT *vT = ns->T[1]->data;
	FLOAT *vu = ns->u[1]->data;
	FLOAT *vH = dof_depth->data;
	ForAllElements(g, e) {
	    int i, k, s, N = ns->u[1]->type->nbas;

	    for (s = 0; s < NFace; s++) {
		if (e->bound_type[s] & BC_BOTTOM) {
		    e->bound_type[s] &= ~BC_BOTTOM2; /* clear mark */

#  if 1
		    /* Method 1: Use all dof on that face */
		    int nbas_face = NbasFace(ns->T[1]);
		    SHORT bases[nbas_face];
		    phgDofGetBasesOnFace(ns->T[1], e, s, bases);

		    for (i = 0; i < nbas_face; i++) {
			INT iD = phgDofMapE2D(ns->T[1], e, bases[i]);

			if (vT[iD] > TEMP_WATER - BETA_MELT * vH[iD] * LEN_SCALING - EPS_T) {
			    /* melt */
			    break;
			} else {
			    /* frozen */
			}
		    }
		    if (i == nbas_face) {
			/* all frozen: froze face */
			e->bound_type[s] |= BC_BOTTOM2;

			nfrozen++;
		    } 	
#  elif 0
		    /* Method 2: Use vert dof on that face */
		    /* Vert */
		    for (j = 0; j < 3; j++) {
			int v = GetEdgeVertex(s, j);
			
		    }
#  else
		    /* Method 3: Use face ceter
		     * */
		    FLOAT vT, vH, 
			lambda[Dim + 1] = {1./3., 1./3., 1./3, 1./3. };
		    lambda[s] = 0.;
		    phgDofEval(ns->T[1], e, lambda, &vT);
		    phgDofEval(ns->depth_P1, e, lambda, &vH);

		    if (vT > TEMP_WATER - BETA_MELT * vH *LEN_SCALING - EPS_T) { 
			/* melting */
		    } else {
			/* frozen */
			e->bound_type[s] |= BC_BOTTOM2;
			int nbas_face = NbasFace(ns->T[1]);
			SHORT bases[nbas_face];
			phgDofGetBasesOnFace(ns->T[1], e, s, bases);

			for (i = 0; i < nbas_face; i++) {
			    INT iD = phgDofMapE2D(ns->u[1], e, bases[i]*Dim);
			    for (k = 0; k < Dim; k++)
				ns->u[1]->data[iD + k] = 0.;
			}
			nfrozen++;
		    }
#   endif

		} /* end bot */
	    }	  /* end face */
	}
	phgInfo(0, "Local frozen: %d\n", nfrozen);
    }

    /* Update BC_BOTTOM2 */
    phgUpdateBoundaryTypes(g);
    phgDofSetBdryDataByFunction(ns->u[1], func_u, BC_BOTTOM2);
    phgDofGradient(ns->u[1], &ns->gradu[1], NULL, "gradu_{n+1}");
#endif
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
#warning vis - theraml coupled
    if (T < T0)
        A = a0 * exp(-Q0 / R / T);
    else
        A = a1 * exp(-Q1 / R / T);
    A *= SEC_PER_YEAR;
#elif 0
#warning const vis A(T = -10)
    T = TEMP_WATER - 10.;
    A = a1 * exp(-Q1 / R / T);
    A *= SEC_PER_YEAR;
#else
#warning const vis A = 1e-16
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
    u[0] = 0;
    u[1] = 0;
    u[2] = 0;
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
#if TEST_CASE == ESIMINT_A
#   warning test case ESIMINT_A
    const FLOAT T_min = 238.15;
    const FLOAT ST = 1.67e-2 / 1e3;
#elif TEST_CASE == ESIMINT_B
#   warning test case ESIMINT_B
    const FLOAT T_min = 243.15;	/* temp +5 */
    const FLOAT ST = 1.67e-2 / 1e3;
#elif TEST_CASE == ESIMINT_F
#   warning test case ESIMINT_F
    const FLOAT T_min = 223.15;	/* Temp -15 */
    const FLOAT ST = 1.67e-2 / 1e3;
#elif TEST_CASE == ESIMINT_C ||			\
    TEST_CASE == ESIMINT_D ||			\
    TEST_CASE == ESIMINT_G ||			\
    TEST_CASE == ESIMINT_H
#   warning test case ESIMINT_[CDGH]
    const FLOAT T_min = 238.15;
    const FLOAT ST = 1.67e-2 / 1e3;
#endif
    const FLOAT x0 = 750e3;
    const FLOAT y0 = 750e3;

    FLOAT d;

    x *= LEN_SCALING;
    y *= LEN_SCALING;

    /* orig zero */
    x -= x0;
    y -= y0;
    d = sqrt(x*x + y*y);

    *T = T_min + ST * d;
}

void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *fT) {
    fT[0] = 0;
}

void func_q(FLOAT x, FLOAT y, FLOAT z, FLOAT *q) {
#if TEST_CASE == ESIMINT_A
#   warning test case ESIMINT_A
    const FLOAT M_max = .5;
    const FLOAT Rel = 450e3 / LEN_SCALING;
    const FLOAT Sb = 1e-2 / 1e3 * LEN_SCALING; 
#elif TEST_CASE == ESIMINT_C
#   warning test case ESIMINT_C
    const FLOAT M_max = .25;		    /* .5->0.25 */
    const FLOAT Rel = 425e3 / LEN_SCALING; /* 450->425 */
    const FLOAT Sb = 1e-2 / 1e3 * LEN_SCALING; 
#elif TEST_CASE == ESIMINT_D
#   warning test case ESIMINT_D
    const FLOAT M_max = .5;
    const FLOAT Rel = 425e3 / LEN_SCALING; /* 450->425 */
    const FLOAT Sb = 1e-2 / 1e3 * LEN_SCALING; 
#elif TEST_CASE == ESIMINT_B ||			\
    TEST_CASE == ESIMINT_E ||			\
    TEST_CASE == ESIMINT_F ||			\
    TEST_CASE == ESIMINT_G ||			\
    TEST_CASE == ESIMINT_H
#   warning test case ESIMINT_[BEFGH]
    const FLOAT M_max = .5;
    const FLOAT Rel = 450e3 / LEN_SCALING;
    const FLOAT Sb = 1e-2 / 1e3 * LEN_SCALING; 
#endif
    const FLOAT x0 = 750e3 / LEN_SCALING;
    const FLOAT y0 = 750e3 / LEN_SCALING;

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
#if TEST_CASE == ESIMINT_G ||\
    TEST_CASE == ESIMINT_H 

    BTYPE DB_masks[3] = {BC_BOTTOM | BC_BOTTOM2,
			 BC_BOTTOM2,
			 BC_BOTTOM2};
#else
    BTYPE DB_masks[3] = {BC_BOTTOM | BC_LATERL, 
			 BC_BOTTOM | BC_LATERL, 
			 BC_BOTTOM | BC_LATERL};
#endif
    phgDofSetDirichletBoundaryMasks(ns->u[1], DB_masks);

    phgDofSetDirichletBoundaryMask(ns->T[1], BC_TOP);

    return;
}

void func_beta(FLOAT x, FLOAT y, FLOAT z, FLOAT *beta) {
#if TEST_CASE == ESIMINT_G ||\
    TEST_CASE == ESIMINT_H
#   warning test case ESIMINT_[GH]
    *beta = 1e3;
#endif
}


