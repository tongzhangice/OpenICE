/* for smoothing data along y
 * poly15: dx0 = 100 m, dy = 500 m, nx = 426, ny = 101 
 * poly17: dx0 = 50 m, dy = 1km, nx = 740, ny = 51
 */

//static double *data_sur, *data_thk, **data_x, **data_y;
static int NX = 17, NY = 114;

static double **data_bed, **data_sur, **data_thk, **data_mask, **data_x, **data_y;
static double **data_sur_grad_x, **data_sur_grad_y;
//double** read_txt_data(char *file_name);
//void interp_txt_data(double **data, double x, double y, double z, double *a);

void
iceInit(GRID *g, LAYERED_MESH **gL_ptr)
{
//    SIMPLEX *e;
    int i, j;
    LAYERED_MESH *gL;
    
    // read data
    data_bed = read_txt_data(ns_params->bed_txt_file);
    data_sur = read_txt_data(ns_params->sur_txt_file);
    data_thk = read_txt_data(ns_params->thk_txt_file);
    data_sur_grad_x = read_txt_data(ns_params->sur_grad_x_txt_file);
    data_sur_grad_y = read_txt_data(ns_params->sur_grad_y_txt_file);


    ice_grid(g);
    phgExportVTK(g, "ice_domain.vtk", NULL);
    
    //free(data_bed[0]);
    //free(data_bed);
    //free(data_sur[0]);
    //free(data_sur);
    //free(data_thk[0]);
    //free(data_thk);
    //free(data_sur_grad_x[0]);
    //free(data_sur_grad_x);
    //free(data_sur_grad_y[0]);
    //free(data_sur_grad_y);

    checkBdry(g);

    /* build layered mesh */
    gL = import_layered_mesh(ns_params->tria_file,
		     ns_params->layer_file,
			     ns_params->nodeZ_file,
			     NULL,
			     phgNProcs);
    build_layered_mesh(g, gL);
    //phgExportVTK(g, "TEST1.vtk", NULL);
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
    if (0) {
	phgExportTecplot(g, "parted.plt", NULL);
	phgExportVTK(g, OUTPUT_DIR "/ins_" NS_PROBLEM "_init.vtk", NULL);
    }
    
    
    for (i = 0; i < gL->nvert; i++) {
    gL->verts[i][0] *= 1000; 
    gL->verts[i][1] *= 1000; 
    //gL->verts[ii][2] *= 1000; 
    }

    TRIA *t = gL->trias;
    for (i = 0; i < gL->ntria; i++, t++) {
    t->area *= 1000*1000;
    } 

    *gL_ptr = gL;

    return;

}


void
func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord)
{

	double bed, sur, thk;
    x *= 1000;
    y *= 1000;    
	interp_txt_data(data_bed, x, y, z, &bed);
	interp_txt_data(data_sur, x, y, z, &sur);
	interp_txt_data(data_thk, x, y, z, &thk);

	if (thk < 1e-10 && thk >= 0) {
	  //printf("thk zero, change to 1m\n.");
      //printf("sur: %lf, bed: %lf\n", sur, bed);
	  thk = 1;
	}

	coord[0] = x;
	coord[1] = y; 
	coord[2] = sur - (1-z)*thk;
    
}



void
func_ice_shelf_mask(FLOAT x, FLOAT y, FLOAT z, FLOAT *ice_shelf_mask)
{
    double bed, sur, thk;
    //double **data_bed, **data_sur, **data_thk;
    x *= 1000;
    y *= 1000;    
	interp_txt_data(data_bed, x, y, z, &bed);
	interp_txt_data(data_sur, x, y, z, &sur);
	interp_txt_data(data_thk, x, y, z, &thk); 
    
    if ((sur-thk-bed)>1){
        ice_shelf_mask[0] = 1;
    }
    else
        ice_shelf_mask[0] = 0;
}





/*
void
func_ice_shelf_pres(FLOAT x, FLOAT y, FLOAT z, FLOAT *ice_shelf_pres)
{
    double bed, sur, thk;
    x *= 1000;
    y *= 1000;    
	interp_txt_data(data_bed, x, y, z, &bed);
	interp_txt_data(data_sur, x, y, z, &sur);
	interp_txt_data(data_thk, x, y, z, &thk); 
    
    if ((sur-thk-bed)>0.00001){
        ice_shelf_pres[0] = -1000*9.8*(sur-thk)/10e5;
        if (ice_shelf_pres[0]<0){
            ice_shelf_pres[0]=0;
        }
    }
    else
        ice_shelf_pres[0] = 0; 
}
*/

void
func_sur_grad_x(FLOAT x, FLOAT y, FLOAT z, FLOAT *sur_grad_x)
{

    double sur_grad;
    x *= 1000;
    y *= 1000;    
	interp_txt_data(data_sur_grad_x, x, y, z, &sur_grad); 
    sur_grad_x[0] = sur_grad;
}
void
func_sur_grad_y(FLOAT x, FLOAT y, FLOAT z, FLOAT *sur_grad_y)
{

    double sur_grad;
    x *= 1000;
    y *= 1000;    
	interp_txt_data(data_sur_grad_y, x, y, z, &sur_grad); 
    sur_grad_y[0] = sur_grad;
}
void
func_ice_sur(FLOAT x, FLOAT y, FLOAT z, FLOAT *ice_sur)
{
    double sur;
    x *= 1000;
    y *= 1000;
    interp_txt_data(data_sur, x, y, z, &sur);
    ice_sur[0] = sur;
}
void
func_bed_z(FLOAT x, FLOAT y, FLOAT z, FLOAT *ice_bed)
{
    double bed;
    x *= 1000;
    y *= 1000;
    interp_txt_data(data_bed, x, y, z, &bed);
    ice_bed[0] = bed;
}
void
func_ice_thk(FLOAT x, FLOAT y, FLOAT z, FLOAT *ice_thk)
{
    double thk;
    x *= 1000;
    y *= 1000;
    interp_txt_data(data_thk, x, y, z, &thk);
    ice_thk[0] = thk;
}

/*
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
*/



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
    //const FLOAT L = _Length_;

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
#  warning const vis A = 6.338e-25 mismip+ setting
    A = 0.3*6.338e-25*SEC_PER_YEAR;
#endif

    if (viscosity_type == VIS_CONST) {
	/* Initial guess:
	 * 1 m/a on 1km surf
	 * */
	eps = 1. / 1000.;
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



void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *fT) {
    fT[0] = 0;
}

void func_beta(FLOAT x, FLOAT y, FLOAT z, FLOAT *beta) 
{
    /* Read from 2D file */

    /* nc_data_set_active(data_index_beta); */
    /* nc_data_scaling = 1.; */
    /* interp_nc_data(x, y, 0, beta); */

    *beta = 1e3;
    //return;
}

void func_q(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) {
    *value = 0.3;
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
    BTYPE DB_masks[3] = {BC_BOTTOM_GRD|BC_LATERL_GRD,
			 BC_BOTTOM_GRD|BC_LATERL_GRD,
			 BC_BOTTOM_GRD|BC_LATERL_GRD};
#else
#   warning BD mask set to (x,0,0)
    BTYPE DB_masks[3] = {BC_LATERL_GRD|BC_BOTTOM_GRD,
			 BC_LATERL_GRD, 
			 BC_LATERL_GRD};
#endif

    phgDofSetDirichletBoundaryMasks(ns->u[1], DB_masks);
    phgDofSetDirichletBoundaryMasks(ns->du, DB_masks);
    phgDofSetDirichletBoundaryMask(ns->T[1], BC_TOP);

    return;
}


void iceSetBoundaryTypes(NSSolver *ns)
{
    /* Do nothing */
}


void func_smb_top(FLOAT x, FLOAT y, FLOAT z, FLOAT *value) 
{
    FLOAT Sb = 0.01;
    FLOAT Rel = 450.0;
    FLOAT Mmax = 0.5;

    FLOAT M0 = Sb*(Rel-sqrt(pow(x/1000,2)+pow(y/1000,2)));

    if (Mmax > M0)
        *value = M0;
    else
        *value = Mmax;


    *value = 0.3;
}

void func_smb_bot(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}

void func_smb_slf(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0;
}


void load_dH_from_file(NSSolver *ns, DOF *dof_P1, int up_or_lower)
{
    GRID *g = ns->g;
    DOF *sn = dof_P1;
    LAYERED_MESH *gL = ns->gL;


    int nx = NX, ny = NY;

    int i, ii, j, k;

    static FLOAT *nsv, *nsv0, *x_coord, *x_coord0, *y_coord, *y_coord0;

    PHG_CALLOC(nsv, gL->nvert);
    PHG_CALLOC(nsv0, gL->nvert);
    

    if (phgRank == 0)
    {
    //FLOAT **data_x = read_txt_data("data_x.txt");
    //FLOAT **data_y = read_txt_data("data_y.txt");

    FLOAT snv[nx][ny];
    FLOAT sum_sn[nx];
    FLOAT mean_sn[nx];
    FLOAT m_avg[nx][ny];

    Bzero(snv); Bzero(sum_sn); Bzero(mean_sn); Bzero(m_avg);



    FILE *fp1;

    if (up_or_lower == 0)
    {
        fp1 = fopen("bot_dh_implicit_fdm.txt", "r");
        if (fp1==NULL){
            printf("errors when opening bot dh file!\n");
        }
    }
    
    if (up_or_lower == 1)
    {
        fp1 = fopen("sur_dh_implicit_fdm.txt", "r");
        if (fp1==NULL){
            printf("errors when opening sur dh file!\n");
        }
    }

	FLOAT **data = (double**) calloc(nx, sizeof *data);
	for (j=0;j<nx;j++)
	{
		data[j] = (double*) calloc(ny, sizeof *data[j]);
		for (i=0;i<ny;i++)
		{
			fscanf(fp1, "%lf", &data[j][i]);
            //if (data[j][i] > 0.3)
            //printf("%f\n", data[j][i]);
		}
	}

    fclose(fp1);

    for (k = 0; k < nx; k++)
        for (j = 0; j < ny; j++)
            m_avg[k][j] = data[k][j];


    for (i = 0; i < gL->nvert; i++)
    {
        for (j = 0; j < ny; j++)
        {
            if (fabs(gL->verts[i][1] - data_y[0][j]) < 1e-1)
            {
                for (k = 0; k < nx; k++)
                {
                    if (fabs(gL->verts[i][0]-data_x[k][0]) < 1e-1)
                    {
                            nsv0[i] = m_avg[k][j];
                    }
                }
            }
        }
    }
    }


    MPI_Bcast(nsv0, gL->nvert,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);

     for (ii = 0; ii < gL->nvert_bot; ii++) { 
     	i = gL->vert_bot_Lidx[ii]; 
     	assert(gL->vert_local_lists[i] != NULL); 

     	INT iG = gL->vert_bot_Gidx[ii]; 
     	assert(iG < gL->nvert); 

     	int nv = gL->vert_local_lists[i][0]; 
     	int *iL = &gL->vert_local_lists[i][1]; 
	    
     	assert(nv > 0); 
        FLOAT *vg;

        if (up_or_lower == 0)
            vg = DofVertexData(sn, iL[0]); 
        if (up_or_lower == 1)
            vg = DofVertexData(sn, iL[nv-1]); 

        vg[0] = nsv0[iG]; 
     } 

}

void save_free_surface_velo(NSSolver *ns, int which_dim, int up_or_lower)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;


    int nx = NX, ny = NY;

    int i, ii, j, k;

    static FLOAT *nsv, *nsv0, *x_coord, *x_coord0, *y_coord, *y_coord0;

    PHG_CALLOC(nsv, gL->nvert);
    PHG_CALLOC(nsv0, gL->nvert);

    for (i = 0; i < gL->nvert; i++)
    {
        nsv[i] = -1e30;
    }
    
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

    INT iG = gL->vert_bot_Gidx[ii];
    assert(iG < gL->nvert);
	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	assert(nv > 0);
	
    FLOAT val;

    FLOAT *vu = DofVertexData(ns->u[1], iL[j]);

    if (up_or_lower == 0)
    {
        vu = DofVertexData(ns->u[1], iL[0]);
        val = vu[which_dim];
    }

    if (up_or_lower == 1)
    {
        vu = DofVertexData(ns->u[1], iL[nv-1]);
        val = vu[which_dim];
    }
    

    nsv[iG] = val;

    }

    MPI_Allreduce(nsv, nsv0, gL->nvert,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    //MPI_Allreduce(x_coord, x_coord0, gL->nvert,
   // 		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    //MPI_Allreduce(y_coord, y_coord0, gL->nvert,
   // 		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

    if (1 && phgRank == 0){
	FILE *fp = fopen("nsv.dat", "w");
	for (ii = 0; ii < gL->nvert; ii++) {
	    fprintf(fp, "%e\n", 
		   nsv0[ii]);
	}
	fclose(fp);
    }

    if (phgRank == 0)
    {
    //FLOAT **data_x = read_txt_data("data_x.txt");
    //FLOAT **data_y = read_txt_data("data_y.txt");

    FLOAT snv[nx][ny];

    Bzero(snv);

    for (i = 0; i < gL->nvert; i++)
    {
        for (j = 0; j < ny; j++)
        {
            if (fabs(gL->verts[i][1] - data_y[0][j]) < 1e-1)
            {
                for (k = 0; k < nx; k++)
                {
                    //printf("%f %f\n", gL->verts[i][0], data_x[0][k]);
                    if (fabs(gL->verts[i][0]-data_x[k][0]) < 1e-1)
                    {
                        snv[k][j] = nsv0[i];
                        //printf("nsv: %f %d %d %d\n", nsv0[i], i, j, k);
                        //printf("nsv: %f %d %d %d\n", snv[j][k], i, j, k);
                    }
                }
            }
        }
    }



    FILE *fp0;

    if (up_or_lower == 0)
    {
        if (which_dim == 0)
            fp0 = fopen("u_bot.txt", "w");
        else if (which_dim == 1)
            fp0 = fopen("v_bot.txt", "w");
        else if (which_dim == 2)
            fp0 = fopen("w_bot.txt", "w");
        else
            printf("Wrong dimension !");

    }

    if (up_or_lower == 1)
    {
        if (which_dim == 0)
            fp0 = fopen("u_sur.txt", "w");
        else if (which_dim == 1)
            fp0 = fopen("v_sur.txt", "w");
        else if (which_dim == 2)
            fp0 = fopen("w_sur.txt", "w");
        else
            printf("Wrong dimension !");
    }


    for (k = 0; k < nx; k++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf(fp0, "%f ", snv[k][j]);
        }
        fprintf(fp0, "\n");
    }

    fclose(fp0);

    }

}

void save_free_surface_elev(NSSolver *ns, int up_or_lower)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;

    int avg_style = 1;

    int nx = NX, ny = NY;

    int i, ii, j, k;

    static FLOAT *nsv, *nsv0, *x_coord, *x_coord0, *y_coord, *y_coord0;

    PHG_CALLOC(nsv, gL->nvert);
    PHG_CALLOC(nsv0, gL->nvert);

    for (i = 0; i < gL->nvert; i++)
    {
        nsv[i] = -1e30;
    }
    
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

    INT iG = gL->vert_bot_Gidx[ii];
    assert(iG < gL->nvert);
	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	assert(nv > 0);
	
    FLOAT val;

    if (up_or_lower == 0)
        val = g->verts[iL[0]][2];
        //val = *DofVertexData(sn, iL[0]);
        // smooth the values at the lower surface
    if (up_or_lower == 1)
        val = g->verts[iL[nv-1]][2];
        //val = *DofVertexData(sn, iL[nv-1]);
        // smooth the values of the upper surface
    

    nsv[iG] = val;

    }

    MPI_Allreduce(nsv, nsv0, gL->nvert,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    //MPI_Allreduce(x_coord, x_coord0, gL->nvert,
   // 		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    //MPI_Allreduce(y_coord, y_coord0, gL->nvert,
   // 		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

    if (1 && phgRank == 0){
	FILE *fp = fopen("nsv.dat", "w");
	for (ii = 0; ii < gL->nvert; ii++) {
	    fprintf(fp, "%e\n", 
		   nsv0[ii]);
	}
	fclose(fp);
    }

    if (phgRank == 0)
    {

    FLOAT snv[nx][ny];

    Bzero(snv);

    for (i = 0; i < gL->nvert; i++)
    {
        for (j = 0; j < ny; j++)
        {
            if (fabs(gL->verts[i][1] - data_y[0][j]) < 1e-1)
            {
                for (k = 0; k < nx; k++)
                {
                    //printf("%f %f\n", gL->verts[i][0], data_x[0][k]);
                    if (fabs(gL->verts[i][0]-data_x[k][0]) < 1e-1)
                    {
                        snv[k][j] = nsv0[i];
                        //printf("nsv: %f %d %d %d\n", nsv0[i], i, j, k);
                        //printf("nsv: %f %d %d %d\n", snv[j][k], i, j, k);
                    }
                }
            }
        }
    }



    FILE *fp0;

    if (up_or_lower == 0)
        fp0 = fopen("s_bot.txt", "w");

    if (up_or_lower == 1)
        fp0 = fopen("s_sur.txt", "w");


    for (k = 0; k < nx; k++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf(fp0, "%f ", snv[k][j]);
        }
        fprintf(fp0, "\n");
    }

    phgPrintf("save surface elevation!\n");

    fclose(fp0);

    }

}

void get_smooth_surface_values(NSSolver *ns, DOF *dof_P1, int up_or_lower)
{
    GRID *g = ns->g;
    DOF *sn = dof_P1;
    LAYERED_MESH *gL = ns->gL;

    int avg_style = 1;

    int nx = NX, ny = NY;

    int i, ii, j, k;

    static FLOAT *nsv, *nsv0, *x_coord, *x_coord0, *y_coord, *y_coord0;

    PHG_CALLOC(nsv, gL->nvert);
    PHG_CALLOC(nsv0, gL->nvert);

    for (i = 0; i < gL->nvert; i++)
    {
        nsv[i] = -1e30;
    }
    
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

    INT iG = gL->vert_bot_Gidx[ii];
    assert(iG < gL->nvert);
	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	assert(nv > 0);
	
    FLOAT val;

    if (up_or_lower == 0)
        val = *DofVertexData(sn, iL[0]);
        // smooth the values at the lower surface
    if (up_or_lower == 1)
        val = *DofVertexData(sn, iL[nv-1]);
        // smooth the values of the upper surface
    

    nsv[iG] = val;

    }

    MPI_Allreduce(nsv, nsv0, gL->nvert,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    //MPI_Allreduce(x_coord, x_coord0, gL->nvert,
   // 		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    //MPI_Allreduce(y_coord, y_coord0, gL->nvert,
   // 		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

    if (0 && phgRank == 0){
	FILE *fp = fopen("nsv.dat", "w");
	for (ii = 0; ii < gL->nvert; ii++) {
	    fprintf(fp, "%e\n", 
		   nsv0[ii]);
	}
	fclose(fp);
    }

    if (phgRank == 0)
    {
    //FLOAT **data_x = read_txt_data("data_x.txt");
    //FLOAT **data_y = read_txt_data("data_y.txt");

    FLOAT snv[nx][ny];
    FLOAT sum_sn[nx];
    FLOAT mean_sn[nx];
    FLOAT m_avg[nx][ny];

    FLOAT data4fit[ny];

    Bzero(snv); Bzero(sum_sn); Bzero(mean_sn); Bzero(m_avg);

    Bzero(data4fit);

    for (i = 0; i < gL->nvert; i++)
    {
        for (j = 0; j < ny; j++)
        {
            if (fabs(gL->verts[i][1] - data_y[0][j]) < 1e-1)
            {
                for (k = 0; k < nx; k++)
                {
                    //printf("%f %f\n", gL->verts[i][0], data_x[0][k]);
                    if (fabs(gL->verts[i][0]-data_x[k][0]) < 1e-1)
                    {
                        snv[k][j] = nsv0[i];
                        //printf("nsv: %f %d %d %d\n", nsv0[i], i, j, k);
                        //printf("nsv: %f %d %d %d\n", snv[j][k], i, j, k);
                    }
                }
            }
        }
    }


#if 0
    FILE *fp0 = fopen("data_orig.txt", "w");

    for (k = 0; k < nx; k++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf(fp0, "%f ", snv[k][j]);
        }
        fprintf(fp0, "\n");
    }

    fclose(fp0);


    system("python get_smooth_data_along_y.py");
#endif


    for (k = 0; k < nx; k++)
    {
        //if (k==0)
        //printf("%f %f %f %f %f %f %f %f %f %f %f %d\n", snv[0][k], snv[1][k], snv[2][k], snv[3][k], snv[4][k], snv[5][k],snv[6][k],snv[7][k],snv[8][k],snv[9][k],snv[10][k],k);
        //FILE *fp = fopen("data_polyfit.txt", "w");

        for (j = 0; j < ny; j++)
        {
            sum_sn[k] += snv[k][j];

            if (j == 0)
                m_avg[k][j] = (snv[k][0] + snv[k][1])/2;
            else if (j == ny-1)
                m_avg[k][j] = (snv[k][ny-2] + snv[k][ny-1])/2;
            else
                m_avg[k][j] = (snv[k][j-1] + snv[k][j] + snv[k][j+1])/3;

            //data4fit[j] = snv[k][j];

            //fprintf(fp, "%d %f\n", j, data4fit[j]);
        }

        //fclose(fp);

        //system("python get_smooth_data_along_y.py");



        mean_sn[k] = sum_sn[k]/ny;
        //printf("mean_dH: %f %f %f\n", snv[j][k], sum_sn[k], mean_sn[k]);
    }


#if 0
    FILE *fp1 = fopen("data_fitted.txt", "r");
    if (fp1==NULL){
        printf("errors when opening fitting file!\n");
    }

	FLOAT **data = (double**) calloc(nx, sizeof *data);
	for (j=0;j<nx;j++)
	{
		data[j] = (double*) calloc(ny, sizeof *data[j]);
		for (i=0;i<ny;i++)
		{
			fscanf(fp1, "%f", &data[j][i]);
		}
	}

    fclose(fp1);

    for (k = 0; k < nx; k++)
        for (j = 0; j < ny; j++)
            m_avg[k][j] = data[k][j];

#endif

    for (i = 0; i < gL->nvert; i++)
    {
        for (j = 0; j < ny; j++)
        {
            if (fabs(gL->verts[i][1] - data_y[0][j]) < 1e-1)
            {
                for (k = 0; k < nx; k++)
                {
                    if (fabs(gL->verts[i][0]-data_x[k][0]) < 1e-1)
                    {
                        if (avg_style == 1)
                            nsv0[i] = mean_sn[k];
#if 0
                        if (avg_style == 2)
                            nsv0[i] = m_avg[k][j];
#endif
                    }
                }
            }
        }
    }
    }


    MPI_Bcast(nsv0, gL->nvert,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);

     for (ii = 0; ii < gL->nvert_bot; ii++) { 
     	i = gL->vert_bot_Lidx[ii]; 
     	assert(gL->vert_local_lists[i] != NULL); 

     	INT iG = gL->vert_bot_Gidx[ii]; 
     	assert(iG < gL->nvert); 

     	int nv = gL->vert_local_lists[i][0]; 
     	int *iL = &gL->vert_local_lists[i][1]; 
	    
     	assert(nv > 0); 
        FLOAT *vg;

        if (up_or_lower == 0)
            vg = DofVertexData(sn, iL[0]); 
        if (up_or_lower == 1)
            vg = DofVertexData(sn, iL[nv-1]); 

        vg[0] = nsv0[iG]; 
     } 

}
