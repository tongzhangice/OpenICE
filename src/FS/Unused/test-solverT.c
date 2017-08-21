/*
 * Test constrained solver T.
 * 
 * Constrain on T is Dirich type, T|_base < T_0.  The constraint
 * happens because of both the geothermal effect and the head
 * conduction effect. To solve such a problem, constrained smoothing
 * is used.
 *
 *
 * Two problems arise for this solver:
 * 1. Will the smoother work? 
 * 2. Since the grid is highly anisotopic, will the smoothing be less
 *    effective for this ?
 *
 * This test will answer the two questions.
 *    
 *  */


/*
 *
 * =========================================
 *                                          
 *  Incompressible Navier-Stokes Equations  
 *                                          
 * =========================================
 *
 * ----------
 * Equations 
 * ----------
 * 
 * Steady static:
 *   -1/Re \laplace u + (u \cdot \grad) u + \grad p = f,
 *   \div u = 0,
 *   u = w on \partial \Omega_D, 
 *   1/Re \frac{\partial u}{\partial n} - n p = g on \partial \Omega_N.
 *
 * Time depentdent:
 *   du / dt - 1/Re \laplace u + (u \cdot \grad) u + \grad p = f, ...
 *
 * ----------
 * Problems
 * ----------
 *
 * Benchmark of flow around cylinder.
 * Schafer M, Turek S. The benchmark problem "Flow around a cylinder". In Hirschel EH (ed.), Flow Simulation
 * with High-Performance Computers II, vol. 52. Notes on Numerical Fluid Mechanics, Vieweg.
 * 
 * 
 *   
 * $Id: ins-main.c,v 1.27 2010-02-05 05:06:38 wleng Exp $ */

#include "ins.h"
#include <stdlib.h>
#include <time.h>
#include <unistd.h>


/***************/
/* GLOBAL vars */
/***************/
NSParams *ns_params = NULL;
FLOAT _Length_ = 5.; 
FLOAT _alpha_ = 0.; 
char vtk_file[1000];


void pc_proc(SOLVER *pc_solver, VEC *b0, VEC **x0) 
{}

SURF_BAS *get_surface_bases(GRID *g, DOF_TYPE *u_type)
{}


/****************/
/* Main program */
/****************/
int
main(int argc, char *argv[])
{
    GRID *g;
    SIMPLEX *e;
    DOF **u, **p, **T, **gradu, **gradT,
	*u_exact, *p_exact, *gradu_exact, *T_exact,
	*eu, *ep, *egradu, *ediv, *eT, *err_ind; 
    FLOAT Time, *dt, res, non_du, non_dp, non_dT;
    INT tstep = 0, nonstep = 0, nelem;
    char mesh_file[100], hostname[256],
	data_file[100], data_u[100], data_p[100], data_T[100];
    size_t mem, mem_peak;
    int verb;
    double tt[3], tt1[3];
    BOOLEAN debug = FALSE, use_smooth_solver = FALSE;
    /* ---------- NS --------- */
    NSSolver *ns = NULL;
    SURF_BAS *surf_bas = NULL;
    LAYERED_MESH *gL = NULL;
    MG_BLOCK_DOFS *bk = NULL;

    /* ================================================================================
     *
     *         Initialize Grid & parameters
     *
     * ================================================================================
     */
    /* Global (static) options */
    Unused(verb);  
    ns_params = phgParametersCreate();     
#if USE_MG 
    mg_params = phgMultiGridParametersCreate();
#endif /* USE_MG */
    phgOptionsRegisterNoArg("debug", "mpi debug", &debug);
    phgOptionsRegisterNoArg("use_smooth_solver", "use smooth solver", &use_smooth_solver);
    phgInit(&argc, &argv);
    phgOptionsShowUsed();

    g = phgNewGrid(-1);
    phgSetPeriodicity(g, ns_params->periodicity);
    phgImportSetBdryMapFunc(my_bc_map);
    phgPrintf("Using mesh: %s\n", ns_params->fn);
    if (!phgImport(g, ns_params->fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", ns_params->fn);

    checkBdry(g);
    elapsed_time(g, FALSE, 0.);	/* reset timer */
    gethostname(hostname, sizeof(hostname));
    printf("#%5d# runing PID %5d on %s \n", phgRank, getpid(), hostname);

    if (debug) {
        int _i_ = 0;
        unsigned int t = 5;
        int pid = getpid();

        gethostname(hostname, sizeof(hostname));
        printf("#### Lengweee debug "
               "PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        while (0 == _i_) {
            MPI_Barrier(MPI_COMM_WORLD);
            printf("%d after bar\n", pid);
            fflush(stdout);
            for (t=0; t<50; t++)
                sleep(1);
            printf("%d after sleep\n", pid);
            fflush(stdout);
        }
        printf("### PID %d, ready!\n", getpid());
    }

    NsSolver_Options();
    phgPrintf("  Pre-refine & repartition ");
    phgRefineAllElements(g, ns_params->pre_refines);

    /* set Reynolds number */
    Time = 0.;
    setFlowParameter(ns_params->Re, ns_params->nu, Time);
    
    /* ================================================================================
     *
     *         build ice grid
     *
     * ================================================================================
     */
    if (!strncmp(NS_PROBLEM, "iceA", 3)
	|| !strncmp(NS_PROBLEM, "esimint_", 8)
	|| !strncmp(NS_PROBLEM, "heino_", 6)
	|| !strcmp(NS_PROBLEM, "greenland"))
    	ice_grid(g);
    checkBdry(g);

    gL = import_layered_mesh(ns_params->tria_file,
			     ns_params->layer_file,
			     ns_params->nodeZ_file,
			     phgNProcs);
    build_layered_mesh(g, gL);
    part_layered_mesh(g, gL);
    destory_layerd_mesh(&gL);
    phgOptionsSetOptions("-partitioner metis -metis_use_parted");
    if (phgBalanceGrid(g, 1.1, -1, NULL, 0.)) {
	phgPrintf("\nRepartition mesh, %d submesh%s, load imbalance: %lg",
		  g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
    }
    phgExportVTK(g, "parted.vtk", NULL);

    gL = import_layered_mesh(ns_params->tria_file,
			     ns_params->layer_file,
			     ns_params->nodeZ_file,
			     phgNProcs);
    build_layered_mesh(g, gL);
    phgExportVTK(g, OUTPUT_DIR "/ins_" NS_PROBLEM "_init.vtk", NULL);





    /* ================================================================================
     *
     *         Create INS solver
     *
     * ================================================================================
     */

    /* Note: pointers u, p, gradu, dt
     *       DIRECTLY access private member of INS solver.
     * */
    phgPrintf("  Create INS solver");
    Time = 0;			
    tstep = 1;			/* time step start at 1 */
    setFuncTime(Time);          /* in file ins-bc.c: static */
    ns = phgNSCreate(g, ns_params);
    u = ns->u; 
    p = ns->p;
    T = ns->T;
    gradu = ns->gradu;
    gradT = ns->gradT;
    dt = ns->dt;		/* direct accses ns  */
    dt[0] = ns_params->dt0;
#if USE_MG 
    if (ns_params->use_mg_F)
	ns->mg = mg;
#endif /* USE_MG */
    ns->gL = gL;
    ns->bk = bk = init_line_block(T[1], gL);
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
    get_height_depth(ns);


    /* ================================================================================
     *
     * 
     *    Main loop:
     *       1. Steady state:   adaptive refinement.
     *       2. Time dependent: time advance.
     *
     * ================================================================================
     * */
    while (TRUE) {
	static BOOLEAN initialized = FALSE;
	FLOAT time_end = ns_params->time_end;

	elapsed_time(g, FALSE, 0.);	/* reset timer */
	phgGetTime(tt);

	if (Fabs(time_end - Time) < 1e-12) {
	    phgPrintf("\n=======\nTime reach end: %lf, exit.\n", Time);
	    break;
	}

	if (tstep > ns_params->max_tstep) {
	    phgPrintf("\n=======\nTime step reach end: %d, exit.\n", 
		      tstep);
	    break;
	}

	dt[-1] = dt[0];
	if (Time + dt[0] > time_end)
	    dt[0] = time_end - Time;

	Time += dt[0];
	setFuncTime(Time);

	phgPrintf("\n==========\ntime: %lf, step:%d\n", (double)Time, tstep);
	phgPrintf("    %d DOF (u:%d, p:%d), %d elements, %d submesh%s, load imbalance: %lg\n",
		  DofGetDataCountGlobal(u[1]) + DofGetDataCountGlobal(p[1]), 
		  DofGetDataCountGlobal(u[1]), DofGetDataCountGlobal(p[1]), 
		  g->nleaf_global, g->nprocs,
		  g->nprocs > 1 ? "es" : "", (double)g->lif);



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
	
	elapsed_time(g, FALSE, 0.);	/* reset timer */
	/*
	 * non-linear iteration.
	 * */
	nonstep = 0;
	non_du = non_dp = non_dT = 1e+10;
	while (TRUE) {
	    phgPrintf("\n   ==================\n");
	    phgPrintf("   Non-linear interation step: %d\n", nonstep);
	    
#if ESIMINT_A || HEINO_TEST
	    /* --------------------------------------------------------------------------------
	     * 
	     *  Step 4.
	     *
	     *   Solve temperature.
	     * 
	     * -------------------------------------------------------------------------------- */

	    phgPrintf("\n   ==================\n");
	    phgPrintf("   Temperature solve \n");
	    phgPrintf("   ==================\n\n");
	    phgPrintf("   T type: %s\n", T[1]->type->name);

	    elapsed_time(g, FALSE, 0.);	/* reset timer */
	    phgNSInitSolverT(ns);
	    phgNSBuildSolverTMat(ns); 
	    phgNSBuildSolverTRHS(ns); 
	    phgSolverUpdateRHS(ns->solver_T);
	    phgNSBuildSolverTConstrain(ns);
	    phgDofCopy(ns->T[1], &ns->dT, NULL, "dT");

	    verb = phgVerbosity; 
	    phgVerbosity--;
	    if (use_smooth_solver)
		phgNSSolverTSolve(ns);
	    else
		phgSolverSolve(ns->solver_T, FALSE, ns->T[1], NULL);

	    phgVerbosity = verb; 
	    phgPrintf("      solver_T: nits = %d, resid = %0.4lg ",
		      ns->solver_T->nits, ns->solver_T->residual);
	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
	    find_melt_region(ns);
	    phgNSDestroySolverT(ns);

	    DOF_SCALE(ns->T[1], "after solve");
	    phgDofAXPY(-1.0, ns->T[1], &ns->dT);
	    non_dT = phgDofNormInftyVec(ns->dT);
	    phgPrintf("   dT: %24.12E\n", non_dT);

	    phgDofGradient(T[1], &gradT[1], NULL, "gradT_{n+1}");
#else
	    non_dT = 0.;
#endif


	    /* Linearized */
	    if (!ns_params->non_linear
		&& nonstep >= 0) {
		phgPrintf("   Linearized iteration converges.\n");
		break;
	    }

	    phgGetTime(tt1);
	    phgPrintf("    time usage of current non step: %lfs\n",
		      (double)(tt1[2] - tt[2]));

	    nonstep++;
	    if (//(res < ns_params->non_tol) 
		 (non_du < ns_params->non_tol
		 && non_dp < ns_params->non_tol
		 && non_dT < ns_params->non_tol)
		|| nonstep > ns_params->max_nonstep) {

		if (nonstep > ns_params->max_nonstep) 
		    phgWarning("   Non-linear iteration reach max step,"
			       " results may be inaccrate!\n");
		else
		    phgPrintf("   Non-linear iteration converges.\n");
		break;
	    }
	} /* solve */



	/* ------------------------------------------------------------
	 * 
	 *   Error check
	 * 
	 * ------------------------------------------------------------ */
#if 1
	eu = ep = egradu = eT = NULL;
	ediv = phgDofDivergence(u[1], NULL, NULL, "err div u");
	phgPrintf(            "            normL2(u, p) = (%20.12E, %20.12E)\n"
			      "            normTemp     = (%20.12E)\n",
			      dofNormL2(u[1]), dofNormL2(p[1]),
			      dofNormL2(T[1])
			      );
	elapsed_time(g, TRUE, 0.);
#endif	/* TEST_CASE == 101 */


	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("    Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		  (double)mem / (1024.0 * 1024.0),
		  (double)mem_peak / (1024.0 * 1024.0));




	/* ------------------------------------------------------------
	 * 
	 *   Output solution
	 * 
	 * ------------------------------------------------------------ */
	if (tstep % ns_params->step_span == 0) { 
#if 1
	    phgPrintf("    Output step %d solution to ensight ", tstep);
	    phgExportEnsightT(g, OUTPUT_DIR "/ins_" NS_PROBLEM , 1, 1,
	    		      u[1], T[1], NULL);  /* ensight */
	    elapsed_time(g, TRUE, 0.);
#endif
	}



        phgGetTime(tt1);
        phgPrintf("    total time usage of current time step: %lfs\n",
		  (double)(tt1[2] - tt[2]));

	if (mem_peak > 1024 * (size_t)ns_params->mem_max * 1024) {
	    phgPrintf("\n=======\nMem usage reach max, exit.\n");
	    break;
	}
	tstep++;
    }				/* end of time advaning */



    /* destroy reused solver */
    if (ns->solver_u != NULL) {
	if (ns_params->use_PCD)
	    phgNSDestroyPc(ns);
	phgNSDestroySolverU(ns);
    }
    if (ns->solver_T != NULL) 
	phgNSDestroySolverT(ns);


    phgNSFinalize(&ns);
    phgFreeGrid(&g);
    phgFinalize();
    phgFree(ns_params);

    return 0;
}
