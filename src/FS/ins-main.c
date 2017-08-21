
#include "ins.h"
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#if USE_PETSC
#  include <petscksp.h>
#endif 


/***************/
/* GLOBAL vars */
/***************/
NSParams *ns_params = NULL;

//static double **data_sur_T;
//double** read_txt_data(char *file_name);
//void interp_txt_data(double **data, FLOAT x, FLOAT y, FLOAT z, FLOAT *a);
char vtk_file[1000];

void func_T(FLOAT x, FLOAT y, FLOAT z, FLOAT *T) {
    //double sur = -0.018*(x-1e5)+200;
    *T = -5+TEMP_WATER;			
}

void func_u0(FLOAT x, FLOAT y, FLOAT z, FLOAT *u)
{
    u[0] = 10;
    u[1] = 0;
    u[2] = 1;
}

/****************/
/* Main program */
/****************/
int
main(int argc, char *argv[])
{
    GRID *g;
    SIMPLEX *e;
    DOF **u, **p, **T, **gradu,
	*eu, *ep, *egradu, *ediv, *eT,
	*dH = NULL;
    FLOAT Time, *dt, res, non_du, non_dp, non_dT;
    INT tstep = 0, nelem;
    char mesh_file[100], hostname[256],
	data_file[100], data_u[100], data_p[100], data_T[100], data_Crd[100];
    size_t mem, mem_peak;
    int verb;
    double tt[3], tt1[3];
    /* ---------- NS ---------- */
    NSSolver *ns = NULL;
    SURF_BAS *surf_bas = NULL;
    LAYERED_MESH *gL = NULL;
    //GEO_INFO *geo = NULL;
    /* MG_BLOCK_DOFS *bk = NULL; */

    /* ================================================================================
     *
     *         Initialize Grid & parameters
     *
     * ================================================================================
     */
    /* Global (static) options */
    Unused(verb); 
     
    ns_params = phgParametersCreate();     
    
    phgInit(&argc, &argv);
    
    phgOptionsShowUsed();

    g = phgNewGrid(-1);
    //phgSetPeriodicity(g, ns_params->periodicity);
    phgImportSetBdryMapFunc(my_bc_map);
    if (ns_params->resume) {
	phgResumeStage(g, &Time, &tstep, mesh_file, data_file);
	phgPrintf("================================\n\n");
	phgPrintf("* RESUME from time:%E, tstep:%d\n", Time, tstep);
	phgPrintf("*             mesh:%s\n", mesh_file);
	phgPrintf("*             data:%s\n", data_file);
	phgPrintf("================================\n");

	if (!phgImport(g, mesh_file, FALSE))
	    phgError(1, "can't read file \"%s\".\n", ns_params->fn);

    }
    else {
	phgPrintf("Using mesh: %s\n", ns_params->fn);
	if (!phgImport(g, ns_params->fn, FALSE))
	    phgError(1, "can't read file \"%s\".\n", ns_params->fn);
    }
    checkBdry(g);
    elapsed_time(g, FALSE, 0.);	/* reset timer */
    gethostname(hostname, sizeof(hostname));
    printf("#%5d# runing PID %5d on %s \n", phgRank, getpid(), hostname);

    NsSolver_Options();
    phgPrintf("  Pre-refine & repartition ");
    phgRefineAllElements(g, ns_params->pre_refines);

    /* Set Reynolds number */
    Time = ns_params->time_start;	/* default: 0 */
    setFuncTime(Time);
    setFlowParameter(ns_params->Re, ns_params->nu, Time);
    


    /* ================================================================================
     *
     *         build ice grid
     *
     * ================================================================================
     */

    iceInit(g, &gL);
    phgPrintf("Geometry initialization done!\n");
    phgExportVTK(g, "ice_domain.vtk", NULL, NULL);


    /* ================================================================================
     *
     *         Create INS solver
     *
     * ================================================================================
     */

    /* Note: pointers u, p, gradu, dt
     *       DIRECTLY access private member of INS solver.
     */
     
    phgPrintf("  Create INS solver");
    tstep = 1;			/* time step start at 1 */
    setFuncTime(Time);          /* in file ins-bc.c: static */

    //data_sur_T = read_txt_data("./LAS_temp_for_mesh.txt");

    ns = phgNSCreate(g, ns_params);

    //get_mask_bot(ns);
    
    ns->time[1] = Time;
    u = ns->u; 
    p = ns->p;
    T = ns->T;
    gradu = ns->gradu;
    dH = ns->dH;
    dt = ns->dt;		/* direct accses ns  */
    dt[0] = ns_params->dt0;
    ns->gL = gL;
    //ns->bk = bk = NULL; //init_line_block(T[1], gL);   /* Use line block */
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

    /* Init height & depth */
    if (gL != NULL){
	get_height_depth(ns);
    }

    /* surf bases */
    //surf_bas = ns->surf_bas;

#if 0
    DOF *beta = phgDofNew(g, DOF_P1, 1, "beta", func_beta);
    phgExportEnsight(g, "check", beta, NULL);
    phgFinalize();
    exit(1);
#endif



    /* ------------------------------------------------------------
     * 
     * Resume dof data.
     * 
     * ------------------------------------------------------------ */
    if (ns_params->resume) {
	FILE *fp = NULL;
	char fname[1000];
	phgResumeStage(g, &Time, &tstep, mesh_file, data_file);

	/* resume coord */
	{
	    const FLOAT *v = DofData(ns->coord);
	    int i, k;
	    //sprintf(data_Crd, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat.Crd", tstep - 1);
	    sprintf(data_Crd, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat.Crd", 2);
	    assert(ns->coord->type == DOF_P1);
	    load_dof_data3(g, ns->coord, data_Crd, mesh_file);
	    for (i = 0; i < g->nvert; i++) 
		for (k = 0; k < Dim; k++)
		    g->verts[i][k] = *(v++);
	    phgGeomInit_(g, TRUE);
	}

#if 1
	/* resmue u_{n-1} */
	//sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", tstep - 2);
	sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", 1);
	DATA_FILE_SURFIX;
	sprintf(fname, "%s.p%03d", data_u, g->rank);
	if ((fp = fopen(fname, "r")) == NULL) {
	    phgPrintf("*  u_{%d} unavailable.\n", tstep - 2);
	} 
	else {
	    fclose(fp);

	    phgDofCopy(u[1], &u[0], NULL, "u_{n}");
	    phgDofCopy(p[1], &p[0], NULL, "p_{n}");
	    phgDofCopy(T[1], &T[0], NULL, "T_{n}");
	    gradu[0] = NULL;

	    load_dof_data3(g, u[0], data_u, mesh_file);
	    load_dof_data3(g, p[0], data_p, mesh_file);
	    load_dof_data3(g, T[0], data_T, mesh_file);
	    phgPrintf("   Resume u_ {%5d}[%8d]:%24.12E p_ {%5d}[%8d]:%24.12E\n", 
		      tstep - 2, DofGetDataCountGlobal(u[0]), phgDofNormL2(u[0]), 
		      tstep - 2, DofGetDataCountGlobal(p[0]), phgDofNormL2(p[0]));
	    phgPrintf("   Resume T_{%5d}[%8d]:%24.12E\n", 
		      tstep - 2, DofGetDataCountGlobal(T[0]), phgDofNormL2(T[0]));
      
	    phgDofGradient(u[0], &gradu[0], NULL, "gradu_{n}");
	    phgDofSetFunction(u[0], DofInterpolation);
	    phgDofSetFunction(p[0], DofInterpolation);
	    //phgDofSetBdryDataByFunction(u[0], func_u, SETFLOW);
	    DOF_SCALE(u[0], "resume");
	    DOF_SCALE(p[0], "resume");
	    DOF_SCALE(T[0], "resume");
	    DOF_SCALE(gradu[0], "resume");

	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
	}

	/* resmue u_{n} */
	//sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", tstep - 1);
	sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", 2);
	DATA_FILE_SURFIX;
	sprintf(fname, "%s.p%03d", data_u, g->rank);
	if ((fp = fopen(fname, "r")) == NULL) {
	    phgError(1, "read Dof data %s failed!\n", data_file);
	} else {
	    fclose(fp);
	    load_dof_data3(g, u[1], data_u, mesh_file);
	    load_dof_data3(g, p[1], data_p, mesh_file);
	    load_dof_data3(g, T[1], data_T, mesh_file);
	    phgPrintf("   Resume u_ {%5d}[%8d]:%24.12E p_ {%5d}[%8d]:%24.12E\n", 
		      tstep - 1, DofGetDataCountGlobal(u[1]), phgDofNormL2(u[1]), 
		      tstep - 1, DofGetDataCountGlobal(p[1]), phgDofNormL2(p[1]));
	    phgPrintf("   Resume T_{%5d}[%8d]:%24.12E\n", 
		      tstep - 1, DofGetDataCountGlobal(T[1]), phgDofNormL2(T[1])); 
      
	    phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
	    phgDofSetFunction(u[1], DofInterpolation);
	    phgDofSetFunction(p[1], DofInterpolation);
	    //phgDofSetBdryDataByFunction(u[1], func_u, SETFLOW);
	    DOF_SCALE(u[1], "resume");
	    DOF_SCALE(p[1], "resume");
	    DOF_SCALE(T[1], "resume");
	    DOF_SCALE(gradu[1], "resume");

	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
	}
#endif

	/* Re init height & depth */
	if (gL != NULL) {
	    get_height_depth(ns);
	    build_layered_mesh_height(g, gL);
	    check_height(g, gL);
	}


	/* reconstruct last time step */
	//Time -= dt[0];
	Time = tstep - 1;
	ns->time[1] = Time;
	ns->time[0] = Time;
	setFuncTime(Time);
#if 0
	phgExportEnsightT(g, OUTPUT_DIR "/ins_" NS_PROBLEM , tstep, tstep, u[1], p[1], NULL);
	phgFinalize();
	return 0;
#endif	/* debug exit */

    }	/* end of resume */


    /* Init temp field */
    DOF_SCALE(ns->beta, "test");
    if (ns_params->solve_temp && tstep == 1) {
        phgPrintf("Init temp field!\n");
	phgNSTempInit(ns);
    }

	
    /* ================================================================================
     *
     * 
     *    Main loop:
     *       1. Steady state:   adaptive refinement.
     *       2. Time dependent: time advance.
     *
     * ================================================================================
     */
    while (TRUE) {


    FLOAT all_memory_usage = phgMemoryUsage(g, NULL)/(1024.0*1024.0);
    if (all_memory_usage > 10000)
        break;

    get_surf_bot_elev(ns);
    DOF *grad_surf_elev = phgDofGradient(ns->surf_elev_P1, NULL, NULL, "gradu_surf_elev");
    ns->grad_surf_elev = grad_surf_elev;
    



        //phgPrintf("modify mask bot!!\n");
        //modify_mask_bot(ns);
        //modify_mask_bot(ns);

    //phgDofSetDataByValue(u[1], 0.);
    //phgDofSetDataByValue(p[1], 0.);

    ns->surf_bas = get_surface_bases(g, DOF_P2);
    surf_bas = ns->surf_bas;

	static BOOLEAN initialized = FALSE;
	FLOAT time_end = ns_params->time_end;

	elapsed_time(g, FALSE, 0.);	/* reset timer */
	phgGetTime(tt);

	if (Fabs(time_end - Time) < 1e-12) {
	    phgPrintf("\n=======\nTime reach end: %lf, exit.\n", Time);
    phgPrintf("Time End %f\n", time_end);
	    break;
	}

	if (tstep > ns_params->max_tstep) {
	    phgPrintf("\n=======\nTime step reach end: %d, exit.\n", 
		      tstep);
	    break;
	}

#if 0
	/* use time t^{n+1} */
	dt[-1] = dt[0];
	if (Time + dt[0] > time_end)
	    dt[0] = time_end - Time;

	Time += dt[0];
	setFuncTime(Time);
#endif

	phgPrintf("\n==========\ntime: %lf, step:%d\n", (double)Time, tstep);
	phgPrintf("    %d DOF (u:%d, p:%d), %d elements, %d submesh%s, load imbalance: %lg\n",
		  DofGetDataCountGlobal(u[1]) + DofGetDataCountGlobal(p[1]), 
		  DofGetDataCountGlobal(u[1]), DofGetDataCountGlobal(p[1]), 
		  g->nleaf_global, g->nprocs,
		  g->nprocs > 1 ? "es" : "", (double)g->lif);

	/* save mesh */
	if (ns_params->record
	    && tstep % ns_params->step_span == 0) {			
	    phgResumeLogUpdate(g, &Time, &tstep, ns_params->fn, NULL);
	}

	if (!initialized) {
	    /* reset mem_peak */
	    phgMemoryUsageReset();
	    initialized = TRUE;
	}




	/* ------------------------------------------------------------
	 * 
	 *   Time marching 
	 * 
	 * ------------------------------------------------------------ */
	/* update variale,  
	 *   call this routine after time update and/or grid change.*/
	phgNSTimeAdvance(ns, Time, tstep);
	phgPrintf("    update solution");
	elapsed_time(g, TRUE, 0.);
	if (ns_params->pin_node)
	    phgNSPinNode(ns);

	/* --------------------------------------------------------------------------------
	 * 
	 *  Step 3.
	 *
	 *  Momentum Equations.
	 *
	 * -------------------------------------------------------------------------------- */
    int mask_iter = 0;
    INT IF_CHANGE_MASK;

#if 1
    while (TRUE)
    {
        // iteration for ice shelf mask updating

    phgPrintf("----------------------------\n");
    phgPrintf("ice shelf mask iteration: %d\n", mask_iter);
    phgPrintf("----------------------------\n");
	elapsed_time(g, FALSE, 0.);	/* reset timer */
	/*
	 * non-linear iteration.
	 * */
	int max_nonstep = 0, newton_start = 0;
	assert(ns_params->utype == DOF_P2);

	/* For nonlinear iter */
	int nonstep = 0; 
	non_du = non_dp = non_dT = 1e+10;
	DOF *u_last = phgDofCopy(u[1], NULL, NULL, "u_last");
	DOF *p_last = phgDofCopy(p[1], NULL, NULL, "p_last");

	FLOAT non_res_last = 1.;
	LTYPE ltype_last = PICARD;

	/* First step, change max non step.  */
	if (tstep == 1) {	
	    if (ns_params->max_nonstep0 > 0)
		max_nonstep = ns_params->max_nonstep0;
	    else
		max_nonstep = ns_params->max_nonstep;
	    if (ns_params->newton_start0 > 0)
		newton_start = ns_params->newton_start0;
	    else
		newton_start = ns_params->newton_start;
	    phgPrintf("   * Set max nonstep to %d for first step.\n", max_nonstep);
	    phgPrintf("   * Set Newton start to %d for first step.\n", newton_start);
	} else {
	    max_nonstep = ns_params->max_nonstep;
	    newton_start = ns_params->newton_start;
	}


	while (TRUE) {
	    phgPrintf("\n   ==================\n");
	    phgPrintf("   Non-linear interation step: %d\n", nonstep);
	    
	    /* Init const viscosity */
	    if (ns_params->start_const_vis &&
		tstep == 0 && nonstep == 0 && mask_iter == 0) {
		phgPrintf("* vis: const\n");
		ns->viscosity_type = VIS_CONST;
	    } else {		
		phgPrintf("* vis: strain\n");
		ns->viscosity_type = VIS_STRAIN;		
	    }
	    
	    sayHello("Non linear solve begin");
	    phgNSInitSolverU(ns);

	    if (nonstep < newton_start)
		ns->ltype = PICARD;
	    else 
		ns->ltype = NEWTON;

	    phgPrintf("   Build RHS: ");
	    phgNSBuildSolverURHS(ns, 1, nonstep, Time);
	    //phgNSBuildSolverURHS(ns);
	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	    ns->non_res = 
		res = phgVecNorm2(ns->solver_u->rhs, 0, NULL);
	    phgPrintf("   nonlinear residual: %24.12E\n", res);

	    /* Restore Picard if no improvement */
#if 1
	    if (ltype_last == NEWTON &&
		res > non_res_last * .75) {
		phgPrintf("   !!! Newton step failed, use Picard to run again\n");
		ns->ltype = PICARD;
		max_nonstep += 5; /* Add more Picard steps  */

		/* resotre dofs:
		 * Fix me: temprature */
		phgDofCopy(u_last, &u[1], NULL, "u_{n+1}");
		phgDofCopy(p_last, &p[1], NULL, "p_{n+1}");
		phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
		phgNSBuildSolverURHS(ns, 1, nonstep, Time);
		ns->non_res = 
		    res = phgVecNorm2(ns->solver_u->rhs, 0, NULL);
		phgPrintf("   nonlinear residual: %24.12E\n", res);
	    } 
#endif

	    /* save non res */
	    non_res_last = res;
	    ltype_last = ns->ltype;
	    
	    /* build matrices */
	    //if (ns_params->use_PCD)
		//phgNSInitPc(ns);

	    phgPrintf("   Build matrices:\n");
	    phgNSBuildSolverUMat(ns, 1, nonstep, Time);
	    //phgNSBuildSolverUMat(ns);
	    phgPrintf("      done ");

	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

#if 0
	    if (ns_params->use_PCD) {
		phgPrintf("   Build Pc: \n");
		phgNSBuildPc(ns);
		phgPrintf("      done ");
		elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
	    }
#endif


	    //elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	    /*
	     * solve equation and update (u, p)
	     * */
	    phgPrintf("solver tol: %E\n", ns->solver_u->rtol);


	    phgSolverSolve(ns->solver_u, TRUE, ns->du, ns->dp, NULL);


	   // elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
        //get_water_pressure(ns);
#if USE_NODAL_LOADS
        get_nodal_force(ns);
        phgMatDestroy(&ns->matF0);
        phgMatDestroy(&ns->matB0);
        phgMatDestroy(&ns->matBt0);
        phgMatDestroy(&ns->matC0);
        //phgSolverDestroy(&ns->solver_u0);
        ns->solver_u0 = NULL;
#endif
        //get_contact(ns);


#if USE_SLIDING_BC
	    rotate_dof_bases(ns->du, surf_bas, FALSE);
#endif
	    phgPrintf("      solver_u: nits = %d, resid = %0.4lg ",
		      ns->solver_u->nits, ns->solver_u->residual);
	    //elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	    /* save dofs */
	    phgDofCopy(u[1], &u_last, NULL, "u_last");
	    phgDofCopy(p[1], &p_last, NULL, "p_last");
        phgExportVTK(g, "u_test.vtk", ns->u[1], NULL);

	    /* nonlinear correction */
	    phgDofAXPY(1.0, ns->du, &u[1]);
	    phgDofAXPY(1.0, ns->dp, &p[1]);
	    assert(u[1]->type == ns_params->utype);
	    assert(p[1]->type == ns_params->ptype);
#if USE_SLIDING_BC
	    //dof_set_normal_data(u[1], surf_bas);
#else
#endif

	    /* non_du = phgDofNormL2(ns->du); */
	    /* non_dp = phgDofNormL2(ns->dp); */
	    non_du = phgDofNormInftyVec(ns->du);
	    non_dp = phgDofNormInftyVec(ns->dp);
	    phgPrintf(" \n  du: %24.12E dp: %24.12E\n", non_du, non_dp);
	    phgPrintf("   u: [%24.12E, %24.12E]\n", 
		      phgDofMinValVec(u[1]), 
		      phgDofMaxValVec(u[1]));
	    phgPrintf("   p: [%24.12E, %24.12E]\n", 
		      phgDofMinValVec(p[1]), 
		      phgDofMaxValVec(p[1]));
	    phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
        //get_avg_gu(ns);
        //gradu[1] = ns->avg_gu;

	    //elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	    //if (ns_params->use_PCD)
		//phgNSDestroyPc(ns);
        //

	    /* evolution of u */
	    //DOF_SCALE(u[1], "after solve");
	    //DOF_SCALE(p[1], "after solve");
#if 0
#  warning check solution U,p
	    sprintf(vtk_file, OUTPUT_DIR "non_%02d_u.vtk", nonstep);
	    phgExportVTK(g, vtk_file, u[1], p[1], NULL);
	    phgExportEnsightT(g, OUTPUT_DIR "ins_" NS_PROBLEM ,
			      nonstep, nonstep, u[1], p[1], T[1],
			      ns->du, ns->dp, ns->dT, NULL);
#endif



#if 0
	    if (FALSE && nonstep % ns_params->step_span == 0) {
		phgPrintf("   Output solution to ensight ");
		phgExportEnsightT(g, OUTPUT_DIR "/ins_" NS_PROBLEM , nonstep, nonstep,
				  u[1], p[1], T[1], NULL);  /* ensight */
		sprintf(vtk_file, OUTPUT_DIR "non_%02d_T.vtk", nonstep);
		phgExportVTK(g, vtk_file , 
			     u[1], p[1], T[1], ns->du, NULL);
		elapsed_time(g, TRUE, 0.);
		//ice_monitor(ns, nonstep);
	    }
#endif


	    /* Linearized */
#if 0
	    if (!ns_params->non_linear
		&& nonstep >= 0) {
		phgPrintf("   Linearized iteration converges.\n");
		break;
	    }
#endif

	    phgGetTime(tt1);
	    phgPrintf("    time usage of current non step: %lfs\n",
		      (double)(tt1[2] - tt[2]));

	    nonstep++;

	    /*
	     * Nonliner iteration break, 
	     *   converge for characteristic value.
	     * Velocity: 100 m/a
	     * Pressure: 1e8 Pa
	     *
	     *  */

        
        INT u_convergence;

        if (mask_iter < 1)
            u_convergence = check_u_convergence0(ns, u[1], p[1], u_last, p_last, ns_params->u_tol0, ns_params->u_tol0);
        else
            u_convergence = check_u_convergence(ns, u[1], p[1], u_last, p_last, ns_params->u_tol, ns_params->u_tol);


	    //elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	    const FLOAT U0 = 100;
	    const FLOAT P0 = 1e8;
	    if ( nonstep >= ns_params->min_nonstep 
            && ns->viscosity_type != VIS_CONST 
            && u_convergence
		    || nonstep > ns_params->max_nonstep)
        {

            if (nonstep > ns_params->max_nonstep) 
            {
                phgPrintf("   Non-linear iteration reach max step,"
                       " results may be inaccrate!\n");
#if 1 && USE_NODAL_LOADS
                /*
                get_nodal_force(ns);
                phgMatDestroy(&ns->matF0);
                phgMatDestroy(&ns->matB0);
                phgMatDestroy(&ns->matBt0);
                phgMatDestroy(&ns->matC0);
                //phgSolverDestroy(&ns->solver_u0);
                ns->solver_u0 = NULL;
                */
                
                phgMatDestroy(&ns->matF);
                phgMatDestroy(&ns->matB);
                phgMatDestroy(&ns->matBt);
                phgMatDestroy(&ns->matC);
                phgSolverDestroy(&ns->solver_u);
                ns->solver_u = NULL;
                phgMapDestroy(&ns->Vmap);
                phgMapDestroy(&ns->Pmap);
                phgDofFree(&ns->u_shape);
                phgDofFree(&ns->p_shape);
#endif
                break;
            }
            else
            {
                phgPrintf("   Non-linear iteration converges.\n");
#if 1 && USE_NODAL_LOADS
                /*
                get_nodal_force(ns);
                phgMatDestroy(&ns->matF0);
                phgMatDestroy(&ns->matB0);
                phgMatDestroy(&ns->matBt0);
                phgMatDestroy(&ns->matC0);
                //phgSolverDestroy(&ns->solver_u0);
                ns->solver_u0 = NULL;
                */
                phgMatDestroy(&ns->matF);
                phgMatDestroy(&ns->matB);
                phgMatDestroy(&ns->matBt);
                phgMatDestroy(&ns->matC);
                phgSolverDestroy(&ns->solver_u);
                ns->solver_u = NULL;
                phgMapDestroy(&ns->Vmap);
                phgMapDestroy(&ns->Pmap);
                phgDofFree(&ns->u_shape);
                phgDofFree(&ns->p_shape);
#endif
                break;
            }
        }

	    phgNSDestroySolverU(ns);
       
	} /* solve */

	phgDofFree(&u_last);
	phgDofFree(&p_last);


	//phgPrintf("Save Dofs\n");
	//save_dof_data3(g, u[1], OUTPUT_DIR"u.dat");
	//save_dof_data3(g, p[1], OUTPUT_DIR"p.dat");

#elif 0
	/* Set velocity, for debugging */
	ns->viscosity_type = VIS_STRAIN;
	phgDofSetDataByFunction(u[1], func_u0);
	//phgDofSetDataByFunction(p[1], func_p0);
	phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
#else
	phgPrintf("Load Dofs\n");

	load_dof_data3(g, u[1], OUTPUT_DIR"u.dat", NULL);
	load_dof_data3(g, p[1], OUTPUT_DIR"p.dat", NULL);
#endif


	/* Project to continuous gradu */

#if 1
    phgPrintf("\n ======================= \n");
    phgPrintf("project velocity gradient");
    phgPrintf("\n======================== \n");
	proj_gradu(ns, ns->gradu[1]);
#endif


	//DOF_SCALE(gradu[1], "grad u");
	//DOF_SCALE(ns->Gradu, "Grad u");

    phgPrintf("\n----------------------------\n");
    phgPrintf("check ice shelf mask status");
    phgPrintf("\n----------------------------\n");


    get_mask_bot(ns);

    if (!(ns_params->another_run_with_updated_mask))
    {
        mask_iter = 1;
        phgPrintf("manually stops another run with updated mask!\n\n\n");


        IF_CHANGE_MASK = if_update_shelf_mask(ns);
        IF_CHANGE_MASK = 0;
        get_mask_bot(ns);
        if (tstep % ns_params->step_span == 0) { 
            phgPrintf("Save water and nodal forces to VTK \n");
            DOF *water_P1 = phgDofCopy(ns->water_force, NULL, DOF_P1, NULL);
            DOF *nodal_P1 = phgDofCopy(ns->nodal_force, NULL, DOF_P1, NULL);
            sprintf(vtk_file, MASK_OUTPUT_DIR "MASK_%05d.vtk", tstep);
            phgExportVTK(g, vtk_file, water_P1,nodal_P1,ns->contact_force,ns->mask_bot,NULL);
            phgDofFree(&water_P1);phgDofFree(&nodal_P1);
        }
    }

    if (mask_iter < 1)
    {
    IF_CHANGE_MASK = if_update_shelf_mask(ns);
    }



    if ((ns_params->another_run_with_updated_mask)){
    if (IF_CHANGE_MASK == 0)
    {
        phgPrintf("The lower surface mask remains unchanged. Stop the iteration of mask updating !\n");

        if (tstep % ns_params->step_span == 0) { 
            get_mask_bot(ns);
            phgPrintf("Save water and nodal forces to VTK \n");
            DOF *water_P1 = phgDofCopy(ns->water_force, NULL, DOF_P1, NULL);
            DOF *nodal_P1 = phgDofCopy(ns->nodal_force, NULL, DOF_P1, NULL);
            sprintf(vtk_file, MASK_OUTPUT_DIR "MASK_%05d.vtk", tstep);
            phgExportVTK(g, vtk_file, water_P1,nodal_P1,ns->contact_force,ns->mask_bot,NULL);
            phgDofFree(&water_P1);phgDofFree(&nodal_P1);
        }

        phgDofFree(&ns->nodal_force);
        phgDofFree(&ns->water_force);
        phgDofFree(&ns->contact_force);

        break;
    }
    }


    if (mask_iter == 1)
    {
        phgPrintf("\n-------------------------\n");
        phgPrintf("ice shelf mask updated \n");
        phgPrintf("--------------------------\n");
        break;
    }


    if (mask_iter < 1)
    if (tstep % ns_params->step_span == 0) { 
        phgPrintf("Save water and nodal forces to VTK \n");
        DOF *water_P1 = phgDofCopy(ns->water_force, NULL, DOF_P1, NULL);
        DOF *nodal_P1 = phgDofCopy(ns->nodal_force, NULL, DOF_P1, NULL);
        sprintf(vtk_file, MASK_OUTPUT_DIR "MASK_%05d.vtk", tstep);
        phgExportVTK(g, vtk_file, water_P1,nodal_P1,ns->contact_force,ns->mask_bot,NULL);
        phgDofFree(&water_P1);phgDofFree(&nodal_P1);
    }
    

    phgDofFree(&ns->nodal_force);
    phgDofFree(&ns->water_force);
    phgDofFree(&ns->contact_force);


    //phgDofFree(&ns->mask_bot);

    //phgDofFree(&ns->water_pressure);
    //phgDofFree(&ns->stress_nn);
    //    phgPrintf("Free stress_nn !\n");
    //phgDofFree(&ns->stress);

    //phgDofFree(&ns->water_pressure1);
    //phgDofFree(&ns->stress_nn1);
    //phgDofFree(&ns->stress1);
    //phgDofFree(&ns->avg_gu);

    mask_iter++;

    }


        proj_gradu(ns, ns->gradu[1]);

        get_stress(ns, ns->gradu[1], ns->p[1]);



#if 0
        if (1) {
            int i, k;
            DOF *Gu[DDim], *gu[DDim], *guDG0, *stress[DDim];
            guDG0 = phgDofCopy(ns->gradu[1], NULL, DOF_P0, NULL);

            for (k = 0; k < DDim; k++) {
            FLOAT *vGu;
            INT n;
            char name[1000];

            sprintf(name, "Gu%d", k);
            Gu[k] = phgDofNew(g, DOF_P1, 1, name, DofNoAction);
            vGu = ns->Gradu->data; /* DOF_P1 */
            n = DofGetDataCount(Gu[k]);
            for (i = 0; i < n; i++)
                Gu[k]->data[i] = vGu[i * DDim + k];

            sprintf(name, "stress%d", k);
            stress[k] = phgDofNew(g, DOF_P1, 1, name, DofNoAction);
            vGu = ns->stress->data; /* DOF_P1 */
            n = DofGetDataCount(stress[k]);
            for (i = 0; i < n; i++)
                stress[k]->data[i] = vGu[i * DDim + k];
            
            sprintf(name, "gu%d", k);
            gu[k] = phgDofNew(g, DOF_P0, 1, name, DofNoAction);
            vGu = guDG0->data;   /* DOF_P0 */
            n = DofGetDataCount(gu[k]);
            for (i = 0; i < n; i++)
                gu[k]->data[i] = vGu[i * DDim + k];
            }

            
            phgExportVTK(g, "stress.vtk",
                    stress[0], stress[1], stress[2],
                    stress[3], stress[4], stress[5],
                    stress[6], stress[7], stress[8],
                    NULL);

            phgExportVTK(g, "Gu.vtk",
                 Gu[0], Gu[1], Gu[2],
                 Gu[3], Gu[4], Gu[5],
                 Gu[6], Gu[7], Gu[8],
                 gu[0], gu[1], gu[2],
                 gu[3], gu[4], gu[5],
                 gu[6], gu[7], gu[8],
                 NULL);
           // phgFinalize();
        }
#endif






	/* --------------------------------------------------------------------------------
	 *
	 *  Step 4.
	 *
	 *   Solve temperature.
	 *
	 * -------------------------------------------------------------------------------- */

#if 1
	if (ns_params->solve_temp) {
	    phgPrintf("\n   ==================\n");
	    phgPrintf("   Temperature solve \n");
	    phgPrintf("   ==================\n\n");
	    phgPrintf("   T type: %s\n", T[1]->type->name);

	    elapsed_time(g, FALSE, 0.);	/* reset timer */

	    phgNSInitSolverT(ns);

	    phgPrintf("   Build Mat: ");
	    phgNSBuildSolverTMat(ns, FALSE);
	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
	    phgPrintf("   Build RHS: ");
	    phgNSBuildSolverTRHS(ns, FALSE);
	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	    phgNSSolverTBuildConstrain(ns);
	    phgDofCopy(ns->T[1], &ns->dT, NULL, "dT");

	    phgNSSolverTSolve(ns, FALSE);
	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	    phgNSDestroySolverT(ns);

	    //find_melt_region(ns);
	    DOF_SCALE(ns->T[1], "after solve");
	    phgDofAXPY(-1.0, ns->T[1], &ns->dT);
	    non_dT = phgDofNormInftyVec(ns->dT);
	    phgPrintf("   dT: %24.12E\n", non_dT);


	    DOF *temp_diff = phgDofCopy(ns->T[1], NULL, NULL, "Td");
	    {
		FLOAT *vt = temp_diff->data;
		const FLOAT *vh = ns->depth_P2->data;

		INT i, n = DofGetDataCount(temp_diff);
		for (i = 0; i < n; i++, vh++, vt++)
		    *vt = TEMP_WATER - BETA_MELT * (*vh) *LEN_SCALING  - (*vt);
	    }
	} else {
	    phgPrintf("Temp not updated.\n");
	    non_dT = 0.;
	}
#endif




	/* ------------------------------------------------------------
	 * 
	 *   Error check
	 * 
	 * ------------------------------------------------------------ */

#if 0
	if (!ns_params->compute_error) {
	    eu = ep = egradu = eT = NULL;
	    ediv = phgDofDivergence(u[1], NULL, NULL, "err div u");
	    phgPrintf(            "            normL2(u, p) = (%20.12E, %20.12E)\n"
				  "            normH1(u)    = (%20.12E)\n"
				  "            normDiv(u)   = (%20.12E)\n"
				  "            normGadu    = (%20.12E)\n",
				  dofNormL2(u[1]), dofNormL2(p[1]),
				  dofNormL2(gradu[1]), dofNormL2(ediv),
				  dofNormL2(ns->gradu[1]));
	    elapsed_time(g, TRUE, 0.);
	} else {
	    /* ----Error check------- */
	    phgPrintf("    Errors: \n");
	    eu = phgDofCopy(u[1], NULL, ns_params->utype, "erru");
	    ep = phgDofCopy(p[1], NULL, ns_params->ptype, "errp");
	    eT = phgDofCopy(T[1], NULL, NULL, "errT");
	    egradu = phgDofCopy(gradu[1], NULL, NULL, "err grad u");
	    ediv = phgDofDivergence(u[1], NULL, NULL, "err div u");

	    adjust_time(- (1. - ns_params->Theta) * dt[0]);
	    restore_time();

	    phgPrintf("            errL2(u, p) = (%20.12E, %20.12E)\n"
		      "            errH1(u)    = (%20.12E)\n"
		      "            errDiv(u)   = (%20.12E)\n",
		      dofNormL2(eu), dofNormL2(ep),
		      dofNormL2(egradu), dofNormL2(ediv));
	    elapsed_time(g, TRUE, 0.);
	    phgDofFree(&egradu);

	    dof_norm_L2(eT);
	}
#endif


	getPecletNum(g, u[1], ns_params->nu, 6);	
	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("    Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		  (double)mem / (1024.0 * 1024.0),
		  (double)mem_peak / (1024.0 * 1024.0));





	/* ------------------------------------------------------------
	 * 
	 *   Move mesh
	 * 
	 * ------------------------------------------------------------ */
    FLOAT ice_volume_last = get_ice_volume(g);
//	DOF *dH_last = phgDofCopy(ns->dH, NULL, NULL, "dH_last");

	if (ns_params->solve_height) {
	    if (gL != NULL) {
		/* Unstructed layered mesh */
		phgPrintf("Move mesh.\n");
        get_surf_dH(ns);
        //DOF *dH_fem = phgDofCopy(ns->dH, NULL, NULL, NULL);
        phgExportVTK(g, "dH_fem1.vtk", ns->dH, NULL);
        get_smooth_surface_values(ns, ns->dH, 0);
        get_smooth_surface_values(ns, ns->dH, 1);
        phgExportVTK(g, "dH_fem.vtk", ns->dH, NULL);

        /*
        save_free_surface_elev(ns, 0);
        save_free_surface_elev(ns, 1);

        save_free_surface_velo(ns, 0, 0);
        save_free_surface_velo(ns, 1, 0);
        save_free_surface_velo(ns, 2, 0);
        save_free_surface_velo(ns, 0, 1);
        save_free_surface_velo(ns, 1, 1);
        save_free_surface_velo(ns, 2, 1);

        if (phgRank == 0)
        {
            system("python get_upper_ds.py");
            system("python get_lower_ds.py");
        }

        load_dH_from_file(ns, ns->dH, 0);
        load_dH_from_file(ns, ns->dH, 1);
        //DOF *dH_fdm = phgDofCopy(ns->dH, NULL, NULL, NULL);

        //phgDofAXPY(-1, dH_fem, &dH_fdm);

        //phgExportVTK(g, "dH_diff.vtk", dH_fdm, NULL);
        *///phgExportVTK(g, "dH_fdm.vtk", ns->dH, NULL);

		get_moved_coord(ns, tstep);
        phgExportVTK(g, "dH_fem2.vtk", NULL);
		move_mesh(ns);
        //phgExportVTK(g, "moved_geo.vtk", NULL);
		//check_height(g, gL);
	    } else {
		/* Structed layered mesh */
		struct_mesh_update(ns, tstep, Time);
	    }

	    phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
	    phgDofGradient(u[0], &gradu[0], NULL, "gradu_{n}");
	} else {
	    phgPrintf("Mesh not moved.\n");
	}

    phgPrintf("\n-----------------------------------------\n");
    phgPrintf("   dH: [%4.2E, %4.2E]\n", 
          phgDofMinValVec(ns->dH), 
          phgDofMaxValVec(ns->dH));
    phgPrintf("------------------------------------------\n\n");

    FLOAT ice_volume = get_ice_volume(g);

    //DOF *u_P1 = phgDofCopy(u[1], NULL, DOF_P1, NULL);
    //phgExportVTK(g, "u_P1.vtk", u_P1, NULL);

#if 1
    //INT dH_convergence = check_surf_convergence(ns, dH_last);

    
    //if (dH_convergence == 1)
    FLOAT dVdt = fabs(ice_volume-ice_volume_last)/ice_volume/dt[0]; 


    if (Time > 10 && dVdt < ns_params->s_tol)
    {
        phgPrintf("-----------------------------------------------------\n");
        phgPrintf("The ice domain reaches a steady state! Model stops! dVdt %e", dVdt);
        phgPrintf("\n-----------------------------------------------------\n");

        break;
    }
    else
    {
        phgPrintf("\n--------------------------------------------------------------\n");
        phgPrintf("The ice domain is still in an unsteady state! Model continues! dVdt %e\n", dVdt);
        phgPrintf("---------------------------------------------------------------\n\n");
    }
#endif
        


    
    update_grounding_line(ns, tstep);
       // after we update the geometry, we need to check the ice shelf mask again





	/* ------------------------------------------------------------
	 * 
	 *   Output solution
	 * 
	 * ------------------------------------------------------------ */
	if (tstep % ns_params->step_span == 0) { 
	    //ice_monitor(ns, tstep);

#if 1
	    /*  Temp check */
	    phgPrintf("    Output solution to VTK");
/* 	    phgExportEnsightT(g, OUTPUT_DIR "/ins_" NS_PROBLEM , Time, tstep, */
/* 	    		      u[1], p[1], T[1], ns->depth_P1, */
/* 	    		      NULL);  /\* ensight *\/ */
        DOF *u_P1 = phgDofCopy(u[1], NULL, DOF_P1, NULL);
        //phgExportVTK(g,"results.vtk",u_P1,p[1],NULL);
	    sprintf(vtk_file, dH_OUTPUT_DIR "ice_%05d.vtk", tstep);
	    phgExportVTK(g, vtk_file, u_P1, NULL);
        phgDofFree(&u_P1);
	    sprintf(vtk_file, dH_OUTPUT_DIR "dH_%05d.vtk", tstep);
	    phgExportVTK(g, vtk_file, ns->dH, NULL);
	    sprintf(vtk_file, dH_OUTPUT_DIR "stress_%05d.vtk", tstep);
	    phgExportVTK(g, vtk_file, ns->stress, NULL);
	    //sprintf(vtk_file, OUTPUT_DIR "mask_%05d.vtk", tstep);
	    //phgExportVTK(g, vtk_file, ns->mask_bot, NULL);
	    

	    elapsed_time(g, TRUE, 0.);
#endif

	    /* Save coord data */
#if 1
	if (tstep % ns_params->step_span_resume == 0) { 
	    {
		//sprintf(data_Crd, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat.Crd", tstep);
		sprintf(data_Crd, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat.Crd", 2);
		assert(ns->coord->type == DOF_P1);
		save_dof_data3(g, ns->coord, data_Crd);
	    }

	    if (ns_params->record) {			
		/* save dof data for time step {n}, {n+1} */
		//sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", tstep - 1);
		sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", 1);
		DATA_FILE_SURFIX;
		save_dof_data3(g, u[0], data_u);
		save_dof_data3(g, p[0], data_p);
		save_dof_data3(g, T[0], data_T);
		phgPrintf("   Save u_ {%5d}[%8d]:%24.12E p_ {%5d}[%8d]:%24.12E\n", 
			  tstep - 1, DofGetDataCountGlobal(u[0]), phgDofNormL2(u[0]), 
			  tstep - 1, DofGetDataCountGlobal(p[0]), phgDofNormL2(p[0]));
		phgPrintf("   Save T_{%5d}[%8d]:%24.12E\n", 
			  tstep - 1, DofGetDataCountGlobal(T[0]), phgDofNormL2(T[0])); 
		DOF_SCALE(u[0], "save");
		DOF_SCALE(p[0], "save");
		DOF_SCALE(T[0], "save");
		DOF_SCALE(gradu[0], "save");
		

		//sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", tstep);
		sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", 2);
		DATA_FILE_SURFIX;
		save_dof_data3(g, u[1], data_u);
		save_dof_data3(g, p[1], data_p);
		save_dof_data3(g, T[1], data_T);
		phgPrintf("   Save u_ {%5d}[%8d]:%24.12E p_ {%5d}[%8d]:%24.12E\n", 
			  tstep, DofGetDataCountGlobal(u[1]), phgDofNormL2(u[1]), 
			  tstep, DofGetDataCountGlobal(p[1]), phgDofNormL2(p[1]));
		phgPrintf("   Save T_{%5d}[%8d]:%24.12E\n", 
			  tstep, DofGetDataCountGlobal(T[1]), phgDofNormL2(T[1]));
		DOF_SCALE(u[1], "save");
		DOF_SCALE(p[1], "save");
		DOF_SCALE(T[1], "save");
		DOF_SCALE(gradu[1], "save");
		phgResumeLogUpdate(g, NULL, NULL, NULL, data_file);
	    } /* end of record */
	    if (gL != NULL) {
		    //check_height(g, gL);
        }
	    sayHello("After record final solution data");
    }
#endif

	}



    //phgDofFree(&dH_last);
	/* clean up */
	//phgDofFree(&ediv);
	//phgDofFree(&eu);
	//phgDofFree(&ep);
	//phgDofFree(&eT);



	/* ----------------------------------------------------------------------
	 * 
	 * Compute drag force FD
	 *
	 * ---------------------------------------------------------------------- 
	 * */






        phgGetTime(tt1);
        phgPrintf("    total time usage of current time step: %lfs\n",
		  (double)(tt1[2] - tt[2]));

	if (mem_peak > 1024 * (size_t)ns_params->mem_max * 1024) {
	    phgPrintf("\n=======\nMem usage reach max, exit.\n");
	    break;
	}
	tstep++;


#if 1
	/* use time t^{n} */
	dt[-1] = dt[0];
	if (Time + dt[0] > time_end)
	    dt[0] = time_end - Time;

	Time += dt[0];
	setFuncTime(Time);
#endif

    phgDofFree(&surf_bas->dof);
    phgDofFree(&ns->surf_bas->dof);
    }				/* end of time advaning */



    /* destroy line block */
    //destroy_line_block(&ns->bk);

    /* destroy reused solver */
    if (ns->solver_u != NULL) {
	if (ns_params->use_PCD)
	    phgNSDestroyPc(ns);
	phgNSDestroySolverU(ns);
    }


    phgNSFinalize(&ns);
	    
    phgFreeGrid(&g);
    phgFinalize();
    phgFree(ns_params);

    return 0;
}
