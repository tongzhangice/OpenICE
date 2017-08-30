
#include "ins.h"
#include "io.h"
#include "periodic.h"
#if USE_PETSC
#  include <petscsys.h>
#  include <petscviewer.h>
#endif


#define _p new_params
#define _nsp (ns->ns_params)
#define _pcd (ns->pcd)

#define SUB_SOLVER_VERB -1

/* non linear solver type: picard or newton */
static const char *noniter_name[] = {
  "picard", "newton", NULL};

/*******************************************************/
/* Get default parameters and user defined parameters. */
/*******************************************************/
NSParams *
phgParametersCreate()
{
    NSParams *new_params = phgCalloc(1, sizeof(*new_params));

    /* default settings */
    _p->fn = "../test/cube.dat";

    /* Get user defined parameter from file*/
    /* geom & mesh files */
    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&_p->fn);
    phgOptionsRegisterFilename("tria_file", "Mesh file", (char **)&_p->tria_file);
    phgOptionsRegisterFilename("vert_file", "Mesh file", (char **)&_p->vert_file);
    phgOptionsRegisterFilename("layer_file", "Mesh file", (char **)&_p->layer_file);
    phgOptionsRegisterFilename("dual_file", "Mesh file", (char **)&_p->dual_file);
    phgOptionsRegisterFilename("nodeZ_file", "Mesh file", (char **)&_p->nodeZ_file);
    phgOptionsRegisterFilename("netcdf_file", "Mesh file", (char **)&_p->nc_file);

    phgOptionsRegisterFilename("vx_txt_file", "Mesh file", (char **)&_p->vx_txt_file);
    phgOptionsRegisterFilename("vy_txt_file", "Mesh file", (char **)&_p->vy_txt_file);
    phgOptionsRegisterFilename("bot_txt_file", "Mesh file", (char **)&_p->bot_txt_file);
    phgOptionsRegisterFilename("sur_txt_file", "Mesh file", (char **)&_p->sur_txt_file);
    phgOptionsRegisterFilename("thk_txt_file", "Mesh file", (char **)&_p->thk_txt_file);
    phgOptionsRegisterFilename("mask_txt_file", "Mesh file", (char **)&_p->mask_txt_file);
    phgOptionsRegisterFilename("x_txt_file", "Mesh file", (char **)&_p->x_txt_file);
    phgOptionsRegisterFilename("y_txt_file", "Mesh file", (char **)&_p->y_txt_file);
    phgOptionsRegisterFilename("sur_grad_x_txt_file", "Mesh file", (char **)&_p->sur_grad_x_txt_file);
    phgOptionsRegisterFilename("sur_grad_y_txt_file", "Mesh file", (char **)&_p->sur_grad_y_txt_file);

    phgOptionsRegisterNoArg("layered_mesh", "Layered mesh", &_p->layered_mesh);
    //phgOptionsRegisterFloat("slop_alpha", "Slop alpha", &_alpha_);
    //phgOptionsRegisterFloat("slop_length", "Slop length", &_Length_);
    phgOptionsRegisterInt("periodicity", "Set periodicity", &_p->periodicity);
    phgOptionsRegisterInt("pre_refines", "Pre refines", &_p->pre_refines);
    phgOptionsRegisterInt("n_bdry_layer", "# of boundary layers", &_p->n_bdry_layer);
    phgOptionsRegisterNoArg("curved_bdry", "Use curved boundary", &_p->curved_bdry);
    phgOptionsRegisterNoArg("curved_refine", "curved refinement", &_p->curved_refine);
    phgOptionsRegisterNoArg("update_bdry_type", "update bdry type", &_p->update_bdry_type);

    /* structured mesh */
    phgOptionsRegisterNoArg("struct_mesh", "Structure mesh", &_p->struct_mesh);
    phgOptionsRegisterInt("Nx", "Nx", &_p->Nx);
    phgOptionsRegisterInt("Ny", "Ny", &_p->Ny);
    phgOptionsRegisterInt("Nz", "Nz", &_p->Nz);


    /* FEM discretization */
    phgOptionsRegisterString("utype", "DOF type for velocity", &_p->utype_name);
    phgOptionsRegisterString("ptype", "DOF type for pressure", &_p->ptype_name);
    phgOptionsRegisterString("T_type", "DOF type for Confrom tensor", &_p->T_type_name);
    phgOptionsRegisterFloat("nu", "viscosity number", &_p->nu);
    phgOptionsRegisterNoArg("enclosed_flow", "Enclosed flow", &_p->enclosed_flow);
    phgOptionsRegisterNoArg("pin_node", "Pin a node for enclosed flow", &_p->pin_node);
    phgOptionsRegisterNoArg("use_moc", "Use Method of Characterics", &_p->use_moc);
    phgOptionsRegisterNoArg("implicit_convect", "Implicit scheme for convetion term", 
			    &_p->implicit_convect); 
    phgOptionsRegisterInt("height_scheme", "Height solver scheme", &_p->height_scheme);
    phgOptionsRegisterInt("init_temp_type", "Init temperature type", &_p->init_temp_type);
    phgOptionsRegisterNoArg("solve_temp", "Solve temperature", &_p->solve_temp);
    phgOptionsRegisterNoArg("solve_height", "Solve height", &_p->solve_height);
    phgOptionsRegisterNoArg("start_const_vis", "Start viscosity as constant", &_p->start_const_vis);
    phgOptionsRegisterNoArg("compensate_equ", "Compensate equations", &_p->compensate_equ);
    phgOptionsRegisterNoArg("add_ice_shelf", "add_ice_shelf", &_p->add_ice_shelf);
    phgOptionsRegisterNoArg("lateral_sia_pres_bc", "lateral_sia_pres_bc", &_p->lateral_sia_pres_bc);
    phgOptionsRegisterNoArg("another_run_with_updated_mask", "another_run_with_updated_mask", &_p->another_run_with_updated_mask);
    phgOptionsRegisterNoArg("compute_error", "Compute error", &_p->compute_error);

    phgOptionsRegisterInt("slip_condition", "slip condition", &_p->slip_condition);
    phgOptionsRegisterFloat("slip_beta2", "slip beta2", &_p->slip_beta2);
    phgOptionsRegisterFloat("slip_alpha2", "slip alpha2", &_p->slip_alpha2);
    phgOptionsRegisterFloat("slip_index", "slip index", &_p->slip_index);


    phgOptionsRegisterFloat("u_tol0", "u tol0", &_p->u_tol0);
    phgOptionsRegisterFloat("p_tol0", "p tol0", &_p->p_tol0);
    phgOptionsRegisterFloat("u_tol", "u tol", &_p->u_tol);
    phgOptionsRegisterFloat("p_tol", "p tol", &_p->p_tol);
    phgOptionsRegisterFloat("s_tol", "s tol", &_p->s_tol);

    /* Time discretization */
    phgOptionsRegisterFloat("dt", "Time step", &_p->dt0);
    phgOptionsRegisterFloat("time_start", "Time start", &_p->time_start);
    phgOptionsRegisterFloat("time_end", "Time end", &_p->time_end);
    phgOptionsRegisterFloat("theta", "Time fraction coef", &_p->Theta);
    phgOptionsRegisterInt("step_span", "Step span to output geo file", &_p->step_span);
    phgOptionsRegisterInt("step_span_resume", "Step span to output dof/Crd file", &_p->step_span_resume);
    phgOptionsRegisterInt("max_time_step", "Max time step", &_p->max_tstep);

    /* Nonlinear system */
    phgOptionsRegisterNoArg("non_linear", "non linear problem", &_p->non_linear);
    phgOptionsRegisterFloat("non_tol", "Nonlinear iteration tolerance", &_p->non_tol);
    phgOptionsRegisterFloat("non_sub_tol", 
			    "Nonlinear iteration sub linear solver tolerance", &_p->non_sub_tol);
    phgOptionsRegisterInt("max_non_step", "Max nonlinear iteration step", &_p->max_nonstep);
    phgOptionsRegisterInt("min_non_step", "Min nonlinear iteration step", &_p->min_nonstep);
    phgOptionsRegisterInt("max_non_step0", "Max nonlinear iteration step for 1st step", &_p->max_nonstep0);
    phgOptionsRegisterInt("newton_start", "newton start", &_p->newton_start);
    phgOptionsRegisterInt("newton_start0", "newton start for 1st step", &_p->newton_start0);
    phgOptionsRegisterNoArg("noniter_temp", "non linear Temperature", &_p->noniter_temp);

    /* Linear system */
    phgOptionsRegisterFloat("eps_diagP", "Diagnal entries for pressure mat", &_p->eps_diagP);
    phgOptionsRegisterFloat("fp_scale", "Dirich scale for Fp", &_p->fp_scale);
    phgOptionsRegisterNoArg("use_PCD", "Use PCD preconditioner", &_p->use_PCD);
    phgOptionsRegisterNoArg("use_Fu", "Use Fu in the preconditioner", &_p->use_Fu);
#if USE_MG 
    phgOptionsRegisterNoArg("use_mg_F", "Use MG preconditioner", &_p->use_mg_F);
    phgOptionsRegisterNoArg("use_mg_Ap", "Use MG preconditioner", &_p->use_mg_Ap);
    phgOptionsRegisterNoArg("use_mg_Qp", "Use MG preconditioner", &_p->use_mg_Qp);
#endif /* USE_MG */

    /* Utils */
    phgOptionsRegisterInt("mem_max", "Max memory per process(MB)", &_p->mem_max);
    phgOptionsRegisterNoArg("resume", "resume from previous step", &_p->resume);
    phgOptionsRegisterNoArg("record", "save data to resume", &_p->record);

    phgOptionsRegisterNoArg("use_symetric", "Use symetric Mat and PC", &_p->use_symetric);
    phgOptionsRegisterNoArg("extern_force", "Extern force", &_p->extern_force); 

    /* Solver options */
    phgOptionsRegisterString("Stokes_opts", "Solver Stokes options", &_p->Stokes_opts);
    phgOptionsRegisterString("F_opts", "Solver F options", &_p->F_opts);
    phgOptionsRegisterString("Fu_opts", "Solver Fu options", &_p->Fu_opts);
    phgOptionsRegisterString("Ap_opts", "Solver Ap options", &_p->Ap_opts);
    phgOptionsRegisterString("Qp_opts", "Solver Qp options", &_p->Qp_opts);
    phgOptionsRegisterString("Gu_opts", "Solver Gu options", &_p->Gu_opts);
    phgOptionsRegisterString("T_opts", "temperature solver options", &_p->T_opts);


    phgOptionsPreset("-mesh_file \"../test/cube.dat\" "
		     "-utype P2 -ptype P1 " 
		     "-T_type P2 "
		     "-nu 1.0 "
		     "-dt 1e-5 "
		     "-time_end 1.0 "
		     "-max_time_step 10 "
		     "-eps_diagP 0. "
		     "-theta 0.5 "
		     "-mem_max 100000 "
		     "-step_span 1 "
		     "-pre_refines 0 "
		     "+enclosed_flow "
		     "-non_linear "
		     "-non_sub_tol 1e-3 "
		     "-periodicity 0 "
		     "+resume "
		     "-record "
		     "+pin_node "
		     "+curved_bdry "
		     "-n_bdry_layer 2 "
		     "-use_PCD "
		     "-fp_scale 1 "
		     "+use_moc "
		     "+use_symetric "
		     "+use_Fu "
		     "-implicit_convect "
		     "+extern_force "
		     "-max_non_step 10 " 
		     "-newton_start 100 " 
		     "-max_non_step0 -1 " 
		     "-newton_start0 -1 " 
		     );		     

    return new_params;
}

void phgDofSetName(DOF *dof, const char *name)
{
    if (dof->name != NULL) {
	phgFree(dof->name);
	dof->name = NULL;
    }

    if (name != NULL) {
	dof->name = phgAlloc(strlen(name) + 1);
	strcpy(dof->name, name);
    }

    return;
}

/***********************************/
/* Navier-stokes equations:	   */
/* Set initial values at time = 0. */
/***********************************/
NSSolver *
phgNSCreate(GRID *g, NSParams *ns_params0)
{
    DOF **u , **p, **T, **gradu;
    DOF_TYPE *utype, *ptype;

    NSSolver *ns = (NSSolver *) phgCalloc(1, sizeof(*ns));    
    ns->pcd = (NSPC *) phgCalloc(1, sizeof(*ns->pcd));
    ns->g = g;
    ns->ns_params = ns_params0;

    /* set utype and ptype */
    {
	char s[128];
	phgOptionsPush();

	sprintf(s, "-dof_type %s", _nsp->utype_name);
	phgOptionsSetOptions(s);
	_nsp->utype = DOF_DEFAULT;

	sprintf(s, "-dof_type %s", _nsp->ptype_name);
	phgOptionsSetOptions(s);
	_nsp->ptype = DOF_DEFAULT;

	sprintf(s, "-dof_type %s", _nsp->T_type_name);
	phgOptionsSetOptions(s);
	_nsp->T_type = DOF_DEFAULT;

	phgOptionsPop();

	assert(_nsp->utype->fe_space == FE_H1);
	       // && _nsp->ptype->fe_space == FE_H1 
	       // && _nsp->utype->order >= _nsp->ptype->order);
    }

    utype = _nsp->utype;
    ptype = _nsp->ptype;

    u = ns->u = ns->u_queue + 1;
    p = ns->p = ns->p_queue + 1;
    T = ns->T = ns->T_queue + 1;
    gradu = ns->gradu = ns->gradu_queue + 1;
    
    ns->f = phgDofNew(g, DOF_ANALYTIC, 3, "f_u", func_f);
    ns->gn[0] = phgDofNew(g, DOF_ANALYTIC, 3, "gxbc", func_g1);
    ns->gn[1] = phgDofNew(g, DOF_ANALYTIC, 3, "gybc", func_g2);
    ns->gn[2] = phgDofNew(g, DOF_ANALYTIC, 3, "gzbc", func_g3);
    _pcd->pbc = phgDofNew(g, ptype, 1, "pressure bc for PCD", DofNoAction);

#if 0
    if (_nsp->enclosed_flow) {
	if (_nsp->pin_node)
	    phgDofSetDirichletBoundaryMask(_pcd->pbc, PINNEDNODE);
	else
	    phgDofSetDirichletBoundaryMask(_pcd->pbc, 0);
    } else {
	/* PCD boundary type for open flow:
	 *   now using more delicated method, and _pcd->pbc is not usesd. */
	phgDofSetDirichletBoundaryMask(_pcd->pbc, SETFLOW);
    }
#endif

    /* Init surface bases */
    //ns->surf_bas = get_surface_bases(g, utype);

    /* Set initial value */
    u[1] = phgDofNew(g, utype, Dim, "u_{n+1}", DofInterpolation);
    p[1] = phgDofNew(g, ptype, 1, "p_{n+1}", DofInterpolation);
    T[1] = phgDofNew(g, _nsp->T_type, 1, "T_{n+1}", DofInterpolation);


#if STEADY_STATE
    /* Set initial values only at boundary.
     * For enclosed flow, pinned node has not been decided,
     *   so it is useless to set pinned node value here. */
    if (ns->set_dirichlet_bc)
        phgDofSetBdryDataByFunction(u[1], func_u, BC_LATERL|BC_TOP);
    else
        phgDofSetBdryDataByFunction(u[1], func_u, BC_LATERL);

    
#else
    /* Set initial values in whole domain */
    phgDofSetDataByFunction(u[1], func_u);
    phgDofSetDataByFunction(p[1], func_p);
    phgDofSetDataByFunction(T[1], func_T);
    PERIODIC_SYNC(u[1]);
    PERIODIC_SYNC(p[1]);
    PERIODIC_SYNC(T[1]);
#endif /* STEADY_STATE */
    gradu[1] = phgDofGradient(u[1], NULL, NULL, "gradu_{n+1}");
    phgDofSetDataByFunction(T[1], func_T);
    phgExportVTK(g, "T.vtk", T[1], NULL);
    PERIODIC_SYNC(T[1]);

    phgExportVTK(g, "ns_init.vtk", u[1], p[1], NULL); /* Output to check grid and DOFs. */

    u[0] = NULL;
    p[0] = NULL;
    T[0] = NULL;
    gradu[0] = NULL;

    u[-1] = NULL;
    p[-1] = NULL;
    T[-1] = NULL;
    gradu[-1] = NULL;

    ns->mask_bot = phgDofNew(g, DOF_P1, 1, "mask bot", DofNoAction);
    /* Init height change varible. */
    ns->dH = phgDofNew(g, DOF_P1, 1, "DH", DofNoAction);
    phgDofSetDirichletBoundaryMask(ns->dH, 0);

    /* Time and dt */
    ns->time = ns->time_queue + 1;
    ns->dt = ns->dt_queue + 1;
    Bzero(ns->time_queue);
    Bzero(ns->dt_queue);

    /* Init nonlinear iteration DOF */
#if STEADY_STATE || TIME_DEP_NON
    ns->du = phgDofCopy(u[1], NULL, NULL, "du");
    ns->dp = phgDofCopy(p[1], NULL, NULL, "dp");
    ns->dT = phgDofCopy(T[1], NULL, NULL, "dT");
    set_boundary_mask(ns);	/* Problem depenedent */
    phgDofSetDataByValue(ns->du, 0.);
    phgDofSetDataByValue(ns->dp, 0.);
    phgDofSetDataByValue(ns->dT, 0.);
#else
    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON*/
    ns->non_res = 1e10;
 
    /* Coordinates */
    ns->coord = phgDofNew(g, DOF_P1, Dim, "Coord P1", func_xyz_);

    /* Slip coeficient */
    ns->beta = phgDofNew(g, DOF_P1, 1, "Beta", DofInterpolation);

    return ns;
}

void phgNSFinalize(NSSolver **ns_ptr)
{
    NSSolver *ns = *ns_ptr;
    DOF **u = ns->u, **p = ns->p, **T = ns->T,
	**gradu = ns->gradu;

#if STEADY_STATE || TIME_DEP_NON 
    phgDofFree(&ns->du);
    phgDofFree(&ns->dp);
    phgDofFree(&ns->dT);
#else
    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */
    
    if (u[-1] != NULL) {
	phgDofFree(&u[-1]);
	phgDofFree(&p[-1]);
	phgDofFree(&T[-1]);
	phgDofFree(&gradu[-1]);
    }

    if (u[0] != NULL) {
	phgDofFree(&u[0]);
	phgDofFree(&p[0]);
	phgDofFree(&T[0]); 
	phgDofFree(&gradu[0]);
    }

    phgDofFree(&u[1]);
    phgDofFree(&p[1]);
    phgDofFree(&T[1]);
    phgDofFree(&gradu[1]);

    phgDofFree(&ns->f);
    phgDofFree(&ns->gn[0]);
    phgDofFree(&ns->gn[1]);
    phgDofFree(&ns->gn[2]);
    phgDofFree(&_pcd->pbc);
    phgDofFree(&ns->dH);
    phgDofFree(&ns->mask_bot);
    phgDofFree(&ns->coord);
    phgDofFree(&ns->beta);
    phgDofFree(&ns->depth_P1);
    phgDofFree(&ns->depth_P2);
    //phgDofFree(&ns->stress_nn);
    //phgDofFree(&ns->water_pressure);
    //phgDofFree(&ns->stress);
    phgDofFree(&ns->Gradu);
    //phgFree(&ns->surf_bas);
    phgFree(ns->pcd);
    phgFree(ns);
    *ns_ptr = NULL;

    return;
}



/*******************************************************************************/
/* Navier-stokes equations:						       */
/* Time advance: shift solutions at different time.			       */
/* Note: 								       */
/* 1. this routine MUST be called AFTER time update,			       */
/*     since DOF f is evauled using t_{n+1}.				       */
/* 2. this routine MUST be called AFTER grid change (refine/coercen, balance), */
/*     since gradu is not interpolated during grid change, instead,	       */
/*     they are recaculated using the interpolated u.			       */
/*******************************************************************************/
void phgNSTimeAdvance(NSSolver *ns, FLOAT time, int tstep)
{
    GRID *g = ns->g;
    DOF **u = ns->u, **p = ns->p, **T = ns->T,
	**gradu = ns->gradu;

    elapsed_time(g, FALSE, 0.);	/* reset timer */

    ns->tstep = tstep;
    ns->time[0] = ns->time[1];
    ns->time[1] = time;
    /* phgPrintf("Time advance: "); */

    /* --------------
     * 1. Last level:
     *    u_{n-1} -> discard,
     *    u_{n}   -> u{n-1}.
     *  */
    phgDofFree(&u[-1]);
    phgDofFree(&p[-1]);
    phgDofFree(&T[-1]);
    if (u[0] != NULL) {
	/* phgPrintf("Free u_{n-1}; "); */

	assert (!strcmp(u[0]->name, "u_{n}"));
	assert (!strcmp(p[0]->name, "p_{n}"));
	assert (!strcmp(T[0]->name, "T_{n}"));
	assert (!strcmp(gradu[0]->name, "gradu_{n}"));

	u[-1] = u[0];
	p[-1] = p[0];
	T[-1] = T[0];
	/* grad u_{n-1} unused */
	phgDofFree(&gradu[0]);

	phgDofSetName(u[-1], "u_{n-1}");
	phgDofSetName(p[-1], "p_{n-1}");
	phgDofSetName(T[-1], "T_{n-1}");

    } 
    
    /* --------------
     * 2. mid level:
     *    copy its upper level's value sequentially from 2nd last level to 2nd top level.
     * */

    assert (!strcmp(u[1]->name, "u_{n+1}"));
    assert (!strcmp(p[1]->name, "p_{n+1}"));
    assert (!strcmp(T[1]->name, "T_{n+1}"));
    u[0] = u[1];
    p[0] = p[1];
    T[0] = T[1];
    gradu[0] = gradu[1];
    phgDofSetName(u[0], "u_{n}");
    phgDofSetName(p[0], "p_{n}");
    phgDofSetName(T[0], "T_{n}");
    phgDofSetName(gradu[0], "gradu_{n}");


    /* --------------
     * 3. top level:
     *    created by copying 2nd top level.
     * */
    phgPrintf("      new u");
    u[1] = phgDofCopy(u[0], NULL, NULL, "u_{n+1}");
    p[1] = phgDofCopy(p[0], NULL, NULL, "p_{n+1}");
    T[1] = phgDofCopy(T[0], NULL, NULL, "T_{n+1}");

    /* Set Dirich B.C */
    //phgPrintf("   Set init flow at bdry, ");
    //phgDofSetBdryDataByFunction(u[1], func_u, SETFLOW); /* No inflow */


    /* Pin_node value */
#if 0
    if (_nsp->pin_node) {
	adjust_time(- (1. - _nsp->Theta) * ns->dt[0]);
	phgDofSetBdryDataByFunction(p[1], func_p, PINNEDNODE);
	restore_time();
    }
    PERIODIC_SYNC(u[1]);
    PERIODIC_SYNC(p[1]);
    PERIODIC_SYNC(T[1]);
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
#endif


    /* Recompute grad u_{n+1},
     * grad T_{n+1} not include. */
    phgPrintf("      gradu");
#if STEADY_STATE
    //gradu[1] = phgDofCopy(gradu[0], NULL, NULL, "gradu_{n+1}");
    gradu[1] = phgDofGradient(u[1], NULL, NULL, "gradu_{n+1}");
#elif TIME_DEP_NON
    gradu[1] = phgDofGradient(u[1], NULL, NULL, "gradu_{n+1}");
#else
    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
    ns->non_res = 1e10;

    /* dt change not here*/

    //phgExportVTK(g, "ns_init.vtk", u[1], p[1], NULL); /* Output to check grid and DOFs. */
    return;
}


/* Pin a node to remove the singularity of the system.
 *
 * Note:
 *   1. g->types_vert[pinned_node] will be changed when refine, coarsen, balance grid,
 *      so after calling this routine, grid should be fixed until linear solver finishes.
 *
 *   2. the pinned vertex need to be local to avoid global change,
 *      TODO: fixed this for the case when all verts are remote.
 *      
 * */
INT 
phgNSPinNode(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, id, pinned_node_id = -1;

    if (!_nsp->enclosed_flow)
    	phgError(1, "Pining node is only needed for enclosed flow.\n");

    phgUpdateBoundaryTypes(g); 	/* clean pinned node */
    MPI_Barrier(g->comm);
    pinned_node_id = -1;

#if PIN_AT_ROOT
    if (g->rank == 0) {
	ForAllElements(g, e) {
	    for (j = 0; j < NVert; j++) {
		if (pinned_node_id == -1) {
		    FLOAT vp;
		    i = e->verts[j];
		    if((g->types_vert[i] & REMOTE) || !(g->types_vert[i] & BDRY_MASK))
		    	continue;

		    assert(g->verts != NULL && g->verts[j] != NULL);
		    g->types_vert[i] |= PINNEDNODE;
		    pinned_node_id = i;
		    printf("   Pin a node at[%d] in e[%d]: (%f, %f, %f), at rank: %d\n",
			   i, e->index,
			   g->verts[i][0], g->verts[i][1], g->verts[i][2], 
			   g->rank);
		    func_p(g->verts[i][0], g->verts[i][1], g->verts[i][2], &vp);
		    printf("      value: %16.8E\n", vp);
		}
	    }
	}
	if (pinned_node_id == -1)
	    phgError(-1, "All verts in rank 0 is either remote or interior, "
		     "leaves no one to pin!\n");
    }
#else
    /* pick one node in each proc to be pinned */
    ForAllElements(g, e) {
	for (j = 0; j < NVert; j++) {
	    if (pinned_node_id == -1) {
		//FLOAT vp;
		i = e->verts[j];
		if(!(g->types_vert[i] & BDRY_MASK))
		    continue;
		pinned_node_id = GlobalVertex(g, i);

		assert(g->verts != NULL && g->verts[i] != NULL);
#if 0
		/* check */
		printf("   [%2d] possible pin node at vert local[%d] global[%d]\n",
		       g->rank, i, pinned_node_id);
		printf("       in e[%d]: (%f, %f, %f)\n",
			e->index, g->verts[i][0], g->verts[i][1], g->verts[i][2]);
		func_p(g->verts[i][0], g->verts[i][1], g->verts[i][2], &vp);
		printf("       value: %16.8E\n", vp);
#endif
	    }
	}
    }

    /* Pin the node with MAX global index */
    id = pinned_node_id;
    MPI_Allreduce(&id, &pinned_node_id, 1, MPI_INT, MPI_MAX, g->comm);
    phgPrintf("   Pin node at vert global[%d]\n", pinned_node_id);
    if (pinned_node_id == -1)
	phgError(-1, "Pin verts err!\n");

    /* Set pinned node type */
#if 0
    for (i = 0; i < g->nvert; i++) {
	if (!(g->types_vert[i] & BDRY_MASK))
	    continue;
	if (GlobalVertex(g, i) == pinned_node_id) {
	    g->types_vert[i] |= PINNEDNODE;
	    printf("   [%2d] pin node at vert[%d]\n",
		    g->rank, i);
	    break;
	}
    }
#endif
    MPI_Barrier(g->comm);
#endif	/* PIN_AT_ROOT */

    /* pinned_node_id = -1 on non-root rank */
#if 0
    ns->pinned_node_id = pinned_node_id; 
    if (_nsp->pin_node) {
	adjust_time(- (1. - _nsp->Theta) * ns->dt[0]);
	phgDofSetBdryDataByFunction(ns->p[1], func_p, PINNEDNODE);
	restore_time();
    }
#endif

    MPI_Barrier(g->comm);
    return ns->pinned_node_id;
}


/******************/
/* Init NS solver */
/******************/
void phgNSInitSolverU(NSSolver *ns)
{
    GRID *g = ns->g;
    MAP *Vmap, *Pmap;
    MAT *pmat[9];
    MAT *pmat0[9];
    //DOF **u = ns->u; 
    //FLOAT *dt = ns->dt;
    //INT verb;

    //Unused(dt);  
    //Unused(verb);  
    /* dof copy */
    ns->u_shape = phgDofCopy(ns->du, NULL, NULL, "u shape");
    ns->p_shape = phgDofCopy(ns->dp, NULL, NULL, "p shape");
  
    /* dof map */
    ns->Vmap = phgMapCreate(ns->u_shape, NULL);
    ns->Pmap = phgMapCreate(ns->p_shape, NULL);

    phgPrintf("   solver_u size [%d, %d]\n", ns->Vmap->nglobal, ns->Vmap->bdry_nglobal);
    phgPrintf("   solver_p size [%d, %d]\n", ns->Pmap->nglobal, ns->Pmap->bdry_nglobal);

    Vmap = ns->Vmap;
    Pmap = ns->Pmap;

    /* matrices */
    ns->matF = phgMapCreateMat(Vmap, Vmap);
    ns->matBt = phgMapCreateMat(Vmap, Pmap);
    ns->matB = phgMapCreateMat(Pmap, Vmap);
    ns->matC = phgMapCreateMat(Pmap, Pmap);

#if USE_NODAL_LOADS
    ns->matF0 = phgMapCreateMat(Vmap, Vmap);
    ns->matBt0 = phgMapCreateMat(Vmap, Pmap);
    ns->matB0 = phgMapCreateMat(Pmap, Vmap);
    ns->matC0 = phgMapCreateMat(Pmap, Pmap);
#endif

    ns->matF->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;
    if (ns->matC != NULL)
	ns->matC->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;

    ns->matB->handle_bdry_eqns = FALSE;
    ns->matBt->handle_bdry_eqns = FALSE;

    /*          | F  Bt|
     *  matNS = |      |
     *          | B  C |
     *  */
    pmat[0] = ns->matF;
    pmat[1] = ns->matBt;
    pmat[2] = ns->matB;
    pmat[3] = ns->matC;

#if USE_NODAL_LOADS
    pmat0[0] = ns->matF0;
    pmat0[1] = ns->matBt0;
    pmat0[2] = ns->matB0;
    pmat0[3] = ns->matC0;

    ns->matNS0 = phgMatCreateBlockMatrix(g, 2, 2, pmat0, NULL, NULL);
    ns->matNS0->mv_data = phgAlloc(sizeof(*ns->matNS0->mv_data));
    ns->matNS0->mv_data[0] = (void *) ns;
    ns->matNS0->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;
#endif

    ns->matNS = phgMatCreateBlockMatrix(g, 2, 2, pmat, NULL, NULL);
    ns->matNS->mv_data = phgAlloc(sizeof(*ns->matNS->mv_data));
    ns->matNS->mv_data[0] = (void *) ns;
    ns->matNS->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;

    /* solver_u */
    /* Note: can't use phgMat2Solver here because build_rhs()
     * requires solver->rhs->map to map U into vector indices */
    phgOptionsPush();
    phgOptionsSetOptions(_nsp->Stokes_opts);
/* #if USE_PETSC */
/*     PetscOptionsInsertFile(MPI_COMM_WORLD, "../options/petsc_asm.opts", PETSC_TRUE); */
/*     //PetscOptionsInsertFile(MPI_COMM_WORLD, "../options/petsc_pcd.opts", PETSC_TRUE); */
/* #endif */
    ns->solver_u = phgSolverCreate(SOLVER_PETSC, ns->u_shape, ns->p_shape,
				   NULL);

#if USE_NODAL_LOADS
    ns->solver_u0 = phgSolverCreate(SOLVER_PETSC, ns->u_shape, ns->p_shape,
				   NULL);
#endif
    phgVerbosity = 0;
    phgOptionsPop();

    {
	phgInfo(0, "nlocal1: %d\n", ns->solver_u->rhs->map->nlocal);
	phgInfo(0, "nlocals: %d\n", ns->solver_u->rhs->map->localsize);
    }

    phgMatDestroy(&ns->solver_u->mat);
    ns->solver_u->mat = ns->matNS;
    ns->solver_u->rhs->mat = ns->solver_u->mat;

#if USE_NODAL_LOADS
    phgMatDestroy(&ns->solver_u0->mat);
    ns->solver_u0->mat = ns->matNS0;
    ns->solver_u0->rhs->mat = ns->solver_u0->mat;
#endif

    if (_nsp->non_linear)
	ns->solver_u->rtol = _nsp->non_sub_tol;

#if STEADY_STATE || TIME_DEP_NON
    phgDofSetDataByValue(ns->du, 0.);
    phgDofSetDataByValue(ns->dp, 0.);
    //phgDofFree(&ns->wind);
    //ns->wind = phgDofCopy(u[1], NULL, NULL, "wind");
#else
    TIME_DEP_LINEAR; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */

#if GIVEN_WIND 
    phgDofSetDataByFunction(ns->wind, func_wind);
#endif /* GIVEN_WIND */

    /* Set boundary types of ice-sheet according to temprature. */
    //iceSetBoundaryTypes(ns);

    return;
}

void phgNSReInitSolverU(NSSolver *ns)
{
    DOF **u = ns->u; 

    if (_nsp->non_linear)
	ns->solver_u->rtol = _nsp->non_sub_tol;

#if STEADY_STATE || TIME_DEP_NON
    phgDofSetDataByValue(ns->du, 0.);
    phgDofSetDataByValue(ns->dp, 0.);
    phgDofFree(&ns->wind);
    ns->wind = phgDofCopy(u[1], NULL, NULL, "wind");
#else
    TIME_DEP_LINEAR; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */

#if GIVEN_WIND 
    phgDofSetDataByFunction(ns->wind, func_wind);
#endif /* GIVEN_WIND */

    return;
}


/**********************/
/* Destroy NS solver  */
/**********************/
void phgNSDestroySolverU(NSSolver *ns)
{
    phgInfo(2, "   Destroy solver U\n");
    phgDofFree(&ns->wind);
    ns->wind = NULL;

    phgMatDestroy(&ns->matF);
    phgMatDestroy(&ns->matB);
    phgMatDestroy(&ns->matBt);
    phgMatDestroy(&ns->matC);

#if USE_NODAL_LOADS
    phgPrintf("Free mats and vecs for contact!\n");
    phgMatDestroy(&ns->matF0);
    phgMatDestroy(&ns->matB0);
    phgMatDestroy(&ns->matBt0);
    phgMatDestroy(&ns->matC0);
    //phgVecDestroy(&ns->vec_rhs0);

    phgSolverDestroy(&ns->solver_u0);
    ns->solver_u0 = NULL;
#endif



    phgSolverDestroy(&ns->solver_u);
    ns->solver_u = NULL;
    phgMapDestroy(&ns->Vmap);
    phgMapDestroy(&ns->Pmap);
    phgDofFree(&ns->u_shape);
    phgDofFree(&ns->p_shape);

}










