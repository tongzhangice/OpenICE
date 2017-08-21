#include "multi-grid.h"
#include "ins.h"
#define _p new_params
#define _mgp (mg->mg_params)

MG_PARAMS *mg_params = NULL;

MG_PARAMS *
phgMultiGridParametersCreate()
{

    MG_PARAMS *new_params = (MG_PARAMS *) phgCalloc(1, sizeof(*new_params));

    _p->Ctype = NULL;
    _p->Ctype_name = NULL;
    _p->coarst_opts = NULL;

    /* Register options */
    phgOptionsRegisterInt("mg_max_level", MG_PREFIX "finest grid level", &_p->max_level);
    phgOptionsRegisterInt("mg_min_level", MG_PREFIX"coarsest grid level", &_p->min_level);
    phgOptionsRegisterInt("mg_maxit", MG_PREFIX"max iteration step", &_p->maxit);
    phgOptionsRegisterInt("mg_n_cycle", MG_PREFIX"num of cycles on each level", &_p->n_cycle);
    phgOptionsRegisterInt("mg_n_pre_smooth", MG_PREFIX"pre smooth times", &_p->n_pre_smooth);
    phgOptionsRegisterInt("mg_n_post_smooth", MG_PREFIX"post smooth times", &_p->n_post_smooth);
    phgOptionsRegisterInt("mg_n_coarst_smooth", MG_PREFIX"coarsest smooth times", &_p->n_coarst_smooth);
    phgOptionsRegisterFloat("mg_rtol", MG_PREFIX"convergence relative tolrence", &_p->rtol);
    phgOptionsRegisterFloat("mg_smooth_damp", MG_PREFIX"smooth damping parameter", &_p->smooth_damp);
    phgOptionsRegisterFloat("mg_correct_damp", MG_PREFIX"correct damping parameter", &_p->correct_damp);
    phgOptionsRegisterFloat("mg_line_ar_tol", MG_PREFIX"line smoother AR tol", &_p->line_ar_tol);
    phgOptionsRegisterNoArg("mg_reuse_mat", MG_PREFIX"reuse mat on finest level", &_p->reuse_mat);
    phgOptionsRegisterNoArg("mg_solve_coarst", MG_PREFIX"use coarsest grid solver", &_p->solve_coarst);
    phgOptionsRegisterString("mg_coarst_opts", MG_PREFIX"coarsest grid solver options", &_p->coarst_opts);
    phgOptionsRegisterNoArg("mg_use_low_order", MG_PREFIX"use low order on finest grid", &_p->use_low_order);
    phgOptionsRegisterNoArg("mg_use_GSline", MG_PREFIX"use GS line smoother", &_p->use_GSline);
    phgOptionsRegisterNoArg("mg_use_upwind", MG_PREFIX"use upwind for convetion", &_p->use_upwind);
    phgOptionsRegisterNoArg("mg_export", MG_PREFIX"export level mesh", &_p->export);
    phgOptionsRegisterString("mg_Ctype", MG_PREFIX"coarse Dof type ", &_p->Ctype_name);
    phgOptionsRegisterNoArg("mg_timing", MG_PREFIX"timing", &_p->timing);
    phgMultiGridSmootherRegister();    

    phgOptionsPreset("-mg_max_level 1 "
		     "-mg_min_level 0 "
		     "-mg_maxit 100 "
		     "-mg_rtol 1e-10 "
		     "-mg_n_cycle 1 "
		     "-mg_n_pre_smooth 1 "
		     "-mg_n_post_smooth 1 "
		     "-mg_n_coarst_smooth 1 "
		     "-mg_smooth_damp 1 "
		     "-mg_correct_damp 1 "
		     "-mg_reuse_mat "
		     "+mg_solve_coarst "
		     );

    return new_params;
}


MULTI_GRID *
phgMultiGridCreate(GRID *g, MG_PARAMS *mg_params0)
{
    MULTI_GRID *mg;

    MG_DEBUGn(2, "MultiGrid crete\n");
    mg = phgCalloc(1, sizeof(*mg));
    mg->grid = g;
    mg->ngrid = 1;
    mg->grids[0] = g;
    mg->mg_params = mg_params0;    

    /* set Ctype */
    {
	char s[128];
	DOF_TYPE *type = DOF_DEFAULT;

	phgOptionsPush();
	if (mg_params0->Ctype_name != NULL) {
	    sprintf(s, "-dof_type %s", mg_params0->Ctype_name);
	    phgOptionsSetOptions(s);
	}
	mg_params0->Ctype = DOF_DEFAULT;
	phgOptionsPop();
	DOF_DEFAULT = type;
    }

    return mg;
}

void 
phgMultiGridInit(MULTI_GRID *mg)
{
    MG_DEBUGn(2, "MultiGrid initailize\n");
    mg->ml = phgCalloc(_mgp->max_level, sizeof(*mg->ml));
    mg->level = -1;		/* init level = -1 */

    /* Unused feathers */
    assert(mg_params->reuse_mat == TRUE);
    assert(mg_params->use_low_order == TRUE);
    assert(mg_params->use_upwind == FALSE);
    assert(mg_params->Ctype == DOF_P1);

    return;
}

void 
phgDestroyMultiGrid(MULTI_GRID **mg_ptr)
{
    MULTI_GRID *mg = *mg_ptr;
    MG_LEVEL *ml;
    int i, j;

    MG_DEBUGn(2, "MultiGrid destroy:\n");
    for (i = 0; i < _mgp->max_level; i++) {
	ml = mg->ml[i];
	if (ml != NULL) {
	    MG_DEBUGn(2, "--- clear level: %d\n", i);

	    for (j = 0; j < ml->nslvt; j++) {
		phgMatDestroy(&ml->P_[j]);
		phgMatDestroy(&ml->R_[j]);
		phgVecDestroy(&ml->x_[j]);
		phgVecDestroy(&ml->f_[j]);
		phgVecDestroy(&ml->r_[j]);
		phgMapDestroy(&ml->map_[j]);
		phgFree(ml->types_vec_[j]);
	    }
	    

	    for (j = 0; j < 10; j++) {
		if (ml->mat_[j] != NULL) {
		    phgMatDestroy(&ml->mat_[j]);
		}
		if (ml->solver_[j] != NULL) 
		    phgSolverDestroy(&ml->solver_[j]);
	    }

	    if (ml->transfer != NULL) {
		MG_TRANS_OP *trans = ml->transfer;

		MG_DISABLED("CR");
		//assert(ml->dofs[0]->type == DOF_CR1);
		phgVecDestroy(&trans->vec_P);
		phgMapDestroy(&trans->map_P);
		phgMatDestroy(&trans->mat_CtoP);
		phgMatDestroy(&trans->mat_PtoF);
		
		phgDofClearCache(NULL, trans->dof_P, NULL, NULL, TRUE);
		if (trans->dof_P->data != NULL)
		    phgFree(trans->dof_P->data);
		phgFree(trans->dof_P);
		phgFree(trans);
	    }

	    if (!ml->dummy) {
		phgFree(ml->sub_grid->elems);
		phgFree(ml->sub_grid->geom->data);
		phgFree(ml->sub_grid->geom);
		phgFree(ml->sub_grid);
		for (j = 0; j < ml->ndof; j++) {
		    phgDofClearCache(NULL, ml->dofs[j], NULL, NULL, TRUE);
		    if (ml->dofs[j]->data != NULL)
			phgFree(ml->dofs[j]->data);
		    phgFree(ml->dofs[j]->name);
		    phgFree(ml->dofs[j]);
		}
		phgFree(ml->dofs);
	    }

	    phgFree(ml);
	}
    }

    for (i = 0; i < mg->ngrid - 1; i++) {
	phgFreeGrid(&mg->grids[i]);
    }

    phgFree(mg->ml);
    phgFree(mg);
    *mg_ptr = NULL;

    return;
}

/*
 * Init vector types, used in GaussSidel iteration
 *
 * TODO: since types_vec indicate only remote flag, it should be dim 1.
 * 
 * */
static void
build_types_vec(MG_LEVEL *ml)
{
    GRID *g = ml->grid;
    SIMPLEX *e;
    DOF *u = ml->u;
    INT i, k, n, dim;

#ifdef MG_VTK_DEUBG
    vtkTmpActorsBegin(0);
    vtkSetColor(0, "white");
    vtkSetTransparent(0, 0.3);
    vtkDrawMesh(0);
    vtkSetColor(0, "blue");
    vtkDrawRemoteFace(0);
    vtkSetTransparent(0, 1.0);
#endif /* MG_VTK_DEUBG */

    dim = u->dim;
    assert(u->g != NULL);
    ml->types_vec = phgCalloc(ml->map->nlocal, sizeof(*ml->types_vec));
    ForAllElements(g, e) {
	BTYPE btype;
	GTYPE gtype;
	INT Il, Iv, N = u->type->nbas;
	
	for (i = 0; i < N; i++) {
	    Il = phgMapE2L(ml->map, 0, e, i * dim);
	    Iv = phgMapL2V(ml->map, Il); 
	    btype = phgDofGetElementBasisInfo(u, e, i,
					      &gtype, NULL, NULL);

	    phgInfo(3, " dof[%3d], btype:%3d, gtype:%3d\n", 
		    Iv, btype, gtype);

	    /* TODO:
	     *   What does REMOTE mean on Element DOF ???
	     *
	     * Bug fixed:
	     *   Can NOT use bit operation here!
	     * */
	    if (!(gtype == VERTEX || 
		  gtype == EDGE   || 
		  gtype == FACE))
	    	continue;

#ifdef MG_VTK_DEUBG
	    /* vtk debug */
	    if (btype & REMOTE) {
		FLOAT *coord = phgDofGetElementCoordinates(u, e, i*dim);
		if (btype & OWNER) {
		    vtkSetColor(0, "red");
		} else {
		    vtkSetColor(0, "yellow");
		}
		sprintf(vtk_tmp_str, NULL);
		vtkDrawPoint(0, coord, vtk_tmp_str);
	    }
#endif /* MG_VTK_DEUBG */

	    /* only mark local vec */
	    assert(Iv >=0 && Iv < ml->map->localsize);
	    if (Iv >= ml->map->nlocal)
		continue;

	    for (k = 0; k < dim; k++)
		ml->types_vec[Iv+k] |= btype;
	}
    }
    for (i = 0, n = 0; i < ml->map->nlocal; i++)
	if (ml->types_vec[i] & REMOTE) {
	    n++;
	}
    MG_DEBUGn(1, "   Dofs on proc boundary: %d\n", n);

#ifdef MG_VTK_DEUBG
    vtkPause(0);
    vtkTmpActorsClear(0);
#endif /* MG_VTK_DEUBG */


    return;
}


void 
phgMultiGridBuildLevel_(MULTI_GRID *mg, int level, int ndof, DOF **dof_s, 
			int nslvt, int *solver2dof, BOOLEAN copy_dof_data, 
			BOOLEAN ordinary_level)
{
    GRID *g = mg->grid;
    MG_LEVEL *ml;
    static int level_built = 0;	/* record for check */
    int i;


    /* make sure build levels top down */
    assert(level == level_built 
	   && level == mg->level+1 
	   && mg->ml[level] == NULL);
    assert(level < _mgp->max_level);
    MG_DEBUGn(1, "-------------------------\n", level);
    MG_DEBUGn(1, "MultiGrid build level: %d\n", level);
    mg->level++;
    level_built++;

    ml = mg->ml[level] = 
	phgCalloc(1, sizeof(*ml));
    ml->mg = mg;
    ml->level = level;
    ml->grid = g;		/* tree! */

    if (!ordinary_level) {	/* dummy level */
	/* no Dofs */
	/* no solvers */
	/* no maps */
	//ml->sub_grid = phgCalloc(1, sizeof(*ml->sub_grid));
	/* no vec */
	/* no tranOP */
	/* no types vec */
	/* no redist map */
	//ml->ndof = ndof;
	ml->nslvt = nslvt;	/* needed to free Dist map */
	phgInfo(1, "Dummy level.\n");
	ml->dummy = TRUE;
	return;
    }

    phgInfo(0, "   Local  Tetra: %5d, Verts: %5d, Eedges: %5d, Faces: %5d\n",
	    g->nleaf, g->nvert, g->nedge, g->nface);
    phgInfo(0, "   Global Tetra: %5d, Verts: %5d, Eedges: %5d, Faces: %5d\n",
	    g->nleaf_global, g->nvert_global, g->nedge_global, g->nface_global);


    /* dofs */
    assert(ml->ndof == 0);
    ml->ndof = ndof;
    ml->dofs = phgCalloc(ndof, sizeof(*ml->dofs));
    for (i = 0; i < ndof; i++) {
	DOF *dof = dof_s[i];
	MEM_DUP(ml->dofs[i], dof, 1);
	ml->dofs[i]->name = strdup(dof->name);
	if (copy_dof_data) {
	    DOF_DATA_DUP(ml->dofs[i], dof);
	} else {
	    ml->dofs[i]->data = ml->dofs[i]->data_vert
		= ml->dofs[i]->data_edge 
		= ml->dofs[i]->data_face
		= ml->dofs[i]->data_elem = NULL;
	}
	if (dof->DB_masks != NULL) {
	    MEM_DUP(ml->dofs[i]->DB_masks, dof->DB_masks, dof->dim);
	}
    }

    /* # of solver types */
    for (i = 0; i < nslvt; i++) {
	int i_dof = solver2dof[i];

	ml->solver_dof[i] = ml->dofs[i_dof];
	ml->solver_type[i] = ml->dofs[i_dof]->type;
	phgInfo(0, "   MG solver[%d] for %s(%d)\n", i,
		ml->solver_type[i]->name, 
		ml->solver_dof[i]->dim);

	if (level == mg_params->max_level - 1
	    && i == MG_SOLVER_TYPE_P) {	/* presure redundent level */
	    ml->redunt_level_[i] = TRUE;
	}
    }
    ml->nslvt = nslvt;

    /* map */
    for (i = 0; i < nslvt; i++) 
	ml->map_[i] = phgMapCreate(ml->solver_dof[i], NULL); 

    /* grid dup */
    MEM_DUP(ml->sub_grid, g, 1);
    MEM_DUP(ml->sub_grid->elems, g->elems, g->nelem);
    MEM_DUP(ml->sub_grid->geom, g->geom, 1);
    DOF_DATA_DUP(ml->sub_grid->geom, g->geom);
    for (i = 0; i < ndof; i++) 
	ml->dofs[i]->g = ml->sub_grid;
    ml->sub_grid->curved = g->curved; /* Note: restrict to p-MG only !!! */
    if (g->curved)
	ml->sub_grid->CurvedCoord = g->CurvedCoord;

    /*
     * Mat is not initialized here.
     * */

    /* vec */
    for (i = 0; i < nslvt; i++) {
	ml->x_[i] = phgMapCreateVec(ml->map_[i], 1);
	ml->f_[i] = phgMapCreateVec(ml->map_[i], 1);
	ml->r_[i] = phgMapCreateVec(ml->map_[i], 1);
	DEBUG_TIME_COST("init map-mat-vec");
    }

    phgMultiGridUpdate(mg); 

    for (i = 0; i < nslvt; i++) {
	set_active_solver(mg, i, -1);
	build_types_vec(ml); 
	ml->types_vec_[i] = ml->types_vec;
	ml->types_vec = NULL;
    }
    set_active_solver(mg, -1, -1);
    DEBUG_TIME_COST("types vec");

    return;
}


#if 0
/*
 * phgMultiGridBuildLevel2:
 * Just a wrapper for phgMultiGridBuildLevel, enable using variable list of Dofs.
 * */
void 
phgMultiGridBuildLevel2(MULTI_GRID *mg, int level, DOF *dof, ...)
{
    int ndof;
    DOF **dofs;
    va_list ap;
    
    dofs = phgAlloc(256 * sizeof(*dofs));

    va_start(ap, dof);
    for (ndof = 0; ndof < 256; ndof++) {
	if (dof == NULL)
	    break;
	else
	    MG_DEBUGn(3, "*** build level dof[%d]: %s\n", ndof, dof->name);
	dofs[ndof] = dof;
	dof = va_arg(ap, DOF *);
    }
    va_end(ap);

    phgMultiGridBuildLevel(mg, level, ndof, dofs, TRUE);    
    phgFree(dofs);
}
#endif

void 
phgMultiGridGetDofs(MG_LEVEL *ml, DOF **dof_ptr, ...)
{
    va_list ap;
    int i = 0;

    assert(ml != NULL);

    va_start(ap, dof_ptr);
    while (TRUE) {
	if (dof_ptr == NULL)
	    break;
	*dof_ptr = ml->dofs[i++];
	dof_ptr = va_arg(ap, DOF **);
    }
    va_end(ap);
    assert(i == ml->ndof);

    return;
}


/* Update grid struct in MultiGrid:
 *    MUST call after grid refine !!!
 * */
void 
phgMultiGridUpdate(MULTI_GRID *mg)
{
    GRID *g = mg->grid;
    MG_LEVEL **ml = mg->ml;
    int i;
    phgInfo(1, "   MG current grid  : %X\n", g);

    /* grid info update on coarse grid */
    for (i = mg->start_level; i < mg->level + 1; i++) {
	phgInfo(1, "   MG level:%d, grid: %X\n", 
		i, ml[i]->grid);

	/* Note: only update grids on the same tree. */
	if (mg->grid != ml[i]->grid)
	    continue;

	MG_DEBUGn(1, "MultiGrid update level: %d\n", i);
	ml[i]->sub_grid->verts = g->verts;
	ml[i]->sub_grid->types_vert = g->types_vert;
	ml[i]->sub_grid->types_edge = g->types_edge;
	ml[i]->sub_grid->types_face = g->types_face;
	ml[i]->sub_grid->types_elem = g->types_elem;
    }

    return;
}



void 
set_active_solver(MULTI_GRID *mg, int islvt, int islv)
{
    int level; 

    mg->active_solver_type = islvt;
    mg->active_solver = islv;

    for (level = mg->start_level; level < mg->level+1; level++) {
	MG_LEVEL *ml = mg->ml[level];
    
	if (islvt >= 0) {
	    assert(islvt < ml->nslvt);
	    ml->u = ml->solver_dof[islvt];
	    ml->P = ml->P_[islvt];
	    ml->R = ml->R_[islvt];
	    ml->x = ml->x_[islvt];
	    ml->f = ml->f_[islvt];
	    ml->r = ml->r_[islvt];
	    ml->map = ml->map_[islvt];
	    ml->types_vec = ml->types_vec_[islvt];
	    ml->redunt_level = ml->redunt_level_[islvt];
	} else {
	    ml->u = NULL;
	    ml->P = NULL;
	    ml->R = NULL;
	    ml->x = NULL;
	    ml->f = NULL;
	    ml->r = NULL;
	    ml->map = NULL;
	    ml->types_vec = NULL;
	    ml->redunt_level = FALSE;
	}

	if (islv >= 0) {
	    ml->mat = ml->mat_[islv];
	    ml->solver = ml->solver_[islv];
	    ml->block_dofs = ml->block_dofs_[islv];
	} else {
	    ml->mat = NULL;
	    ml->solver = NULL;
	    ml->block_dofs = NULL;
	}
    }

    return;
}



/******************************/
/* prolongation & restriction */
/******************************/

/* static variable using in traversing */
//static int level_traverse;
static DOF *uC_traverse;
static DOF *uF_traverse;
static MAT *prlg_traverse;
static MAT *rstr_traverse;
static GRID *g_traverse;	

static INT phgTraverseIndex_ = 0; 
static INT phgTraverseDepth_ = 0; /* Depth of the element in the tree.
				     This variable is made global and may
				     be used by the callback functions */
static SIMPLEX *phgTraverseElem0_ = NULL; 


static VEC *vec_scale;		  /* vector record scale of mat entries */
static int prlg_mat_size  = -1;
static FLOAT *prlg_stack;
static INT *Ic = NULL;		  /* dof index of coarse grid */
static INT *If = NULL;		  /* dof index of fine grid */
static INT Isize = 0;


/* Preoder traversal:
 *   traverse_element in utils.c is postorder traversal.
 * */
static BOOLEAN
traverse_element_preorder(SIMPLEX *e, BOOLEAN(*cb) CB_ARGS(e))
{
    BOOLEAN is_leaf = TRUE;

    phgTraverseIndex_++;
    if(!cb CB_ARGS0(g_traverse, e)) 
	return FALSE;

    if (e->children[0] != NULL) {
        is_leaf = FALSE;
	phgTraverseDepth_++;
	if (!traverse_element_preorder(e->children[0], cb))
	    return FALSE;
	phgTraverseDepth_--;
    }
    if (e->children[1] != NULL) {
        is_leaf = FALSE;
	phgTraverseDepth_++;
	if (!traverse_element_preorder(e->children[1], cb))
	    return FALSE;
	phgTraverseDepth_--;
    } 
    
    return TRUE;
}


#define pj_new(i, j) pj_new[(i) * N + j]
#define pj_old(i, j) pj_old[(i) * N + j]
#define pj_new_v(i) pj_new(i, j)
#define pj_new_e(i) pj_new(i+NVert, j)
#define pj_old_v(i) pj_old(i, j)
#define pj_old_e(i) pj_old(i+NVert, j)
/*
 * Local grid transfor operator: Pn -> Pn
 *
 * This routine is faster than prlg_callbackPn. It uses alberta-like
 * refinement information to get the relation of bases of elements from
 * corase and fine grid.
 * 
 * */
static BOOLEAN
prlg_callback_Pn CB_ARGS(e)
{
    SIMPLEX *e0 = e->parent;
    FLOAT *pj_new, *pj_old;
    DOF *dof = uF_traverse;
    int N = dof->type->nbas;
    int i, j, k, Vmap[NVert], Emap[NEdge];
    size_t n_ = N*sizeof(*pj_new);
    
    if (phgTraverseDepth_ == 0) {
    } else {
	assert(e0 != NULL );
	assert(e0->children[0] == e || e0->children[1] == e);
    }

    /* traverse debug */
    MG_DEBUGn(3, "### ");
    for (i = 0; i < phgTraverseDepth_ - 1; i++)
	MG_DEBUGn(3, "### ");
    if (IsLeaf(e)) {
	MG_DEBUGn(3, "##O Traver: %5d\n", e->index);
    } else {
	MG_DEBUGn(3, "### Traver: %5d\n", e->index);
    }
 
    /* get prlg matrix on element */
    if (phgTraverseDepth_ == 0) {
	/* depth 0: use eye mat */
	pj_new = prlg_stack;
	memset(pj_new, 0, prlg_mat_size * sizeof(*pj_new));
	for (i = 0; i < N; i++)
	    pj_new(i, i) = 1.;
    } else {
	/* depth n: inherit from depth n-1 */
	pj_old = prlg_stack + (phgTraverseDepth_ - 1) * prlg_mat_size;
	pj_new = prlg_stack + (phgTraverseDepth_) * prlg_mat_size;

	if (e == e0->children[0]) 
	    k = 0;
	else 
	    k = 1;
	
	phgMapC2P(e0, Vmap, Emap, k);
	
	/* DOF P1 interp */
	if (dof->type == DOF_P1) {
	    for (i = 0; i < NVert - 1; i++) 
		memcpy(pj_new + i*NVert, pj_old + Vmap[i]*NVert, n_);
	    for (j = 0; j < NVert; j++)	
		pj_new_v(3) = .5 * (pj_old_v(0) + pj_old_v(1));	
	} 
	/* DOF P2 interp */
	else if (dof->type == DOF_P2) {
	    /* old vertex */
	    for (i = 0; i < NVert - 1; i++) 
		memcpy(pj_new + i*N, pj_old + Vmap[i]*N, n_);

	    /* new vertex: vnew = eold0 */
	    for (j = 0; j < N; j++)	
		pj_new_v(3) = pj_old_e(0);

	    for (i = 0; i < NEdge; i++) {
		switch (Emap[i]) {
		case -1:	/* new edge */
		    if (Vmap[GetEdgeVertex(i, 0)] == 2) {
			for (j = 0; j < N; j++) 
			    pj_new_e(i) = 1./4. * pj_old_e(0) - 1./8. * pj_old_v(0) - 1./8. * pj_old_v(1) 
				+ 0.5 * (pj_old_e(1) + pj_old_e(3));
		    } else {
			for (j = 0; j < N; j++) 
			    pj_new_e(i) = 1./4. * pj_old_e(0) - 1./8. * pj_old_v(0) - 1./8. * pj_old_v(1) 
				+ 0.5 * (pj_old_e(2) + pj_old_e(4));
		    }
		    break;
		case 0:		/* cut edge */
		    if (k == 0) 
			for (j = 0; j < N; j++)
			    pj_new_e(i) = 3./8. * pj_old_v(0) - 1./8. * pj_old_v(1) + 3./4. * pj_old_e(0);
		    else
			for (j = 0; j < N; j++)
			    pj_new_e(i) = 3./8. * pj_old_v(1) - 1./8. * pj_old_v(0) + 3./4. * pj_old_e(0);
		    break;
		default:	/* old edge */
		    memcpy(pj_new + (i+NVert)*N, pj_old + (NVert+Emap[i])*N, n_);
		    break;
		}
	    }
	} else {
	    phgError(1, "dof type %s has not implement interp Pn proj!!!\n", dof->type->name);
	}
    }
    //SHOW_M(pj_new, NVert, NVert);
    
    /* Add to global prlgection matrix if leaf node */
    if (IsLeaf(e)) {
	int dim = dof->dim;
	if (dim == 1) {
	    for (i = 0; i < N; i++)
		If[i] = phgMapE2L(prlg_traverse->rmap, 0, e, i);

	    REMOVE_ZERO(pj_new, prlg_mat_size);
	    phgMatSetEntries(prlg_traverse, N, If, N, Ic, pj_new); 
	} else {
	    FLOAT pj_dim[N*dim][N*dim];
	    Bzero(pj_dim);
	    REMOVE_ZERO(pj_new, prlg_mat_size);
	    for (i = 0; i < N; i++) {
		for (k = 0; k < dim; k++) {
		    If[i*dim + k] = phgMapE2L(prlg_traverse->rmap, 
					      0, e, i * dim + k);
		    for (j = 0; j < N; j++)
			pj_dim[i*dim+k][j*dim+k] = pj_new[i*N + j];
		}
	    }
	    //SHOW_iV(If, N*dim);
	    phgMatSetEntries(prlg_traverse, N*dim, If, N*dim, Ic, &pj_dim[0][0]); 
	}	
    }

    return TRUE;
}
#undef pj_new_v
#undef pj_new_e
#undef pj_old_v
#undef pj_old_e
#undef pj_new
#undef pj_old


/*
 * Local grid transfor operator: Pm -> Pn or CR1
 *
 * This routine build local grid transfor mat from bases of element of coarse grid,
 * which is of dof type Pm, bases of element of coarse grid, which it Pn.
 *
 * The routine use geometric coordnates as a connection between bases of elements
 * in coarse and fine grid. It's slower than prlg_callback_P1.
 *
 * The transfor matrix of element e to its antecedents E is of size Ne * NE,
 * where Ne is the num of bases in e and NE is the num of bases in E.
 * Because the interpolation relationship reads:
 *    u^e(lambda^e_i) = u^E(x(lambda^e_i))
 *                    = \sum_{j = 1, NE} u^E_j * \phi^E_j(lambda^E(x(lambda^e_i))
 * so the entries of prlgect mat of e to E is,
 *    prlg_mat_{i}{j} = \phi^E_j(lambda^E(x(lambda^e_i))).
 *
 * 1. Prolongation from coarse CONTINUOUS space to fine space:
 *    Use interpolation.
 * 2. Prolongation from coarse DISCONTINUOUS space to fine space:
 *    Use weighted interpolation. For each fine element basis, use ALL coarse element
 *    that contains it.
 * 
 *  */
static BOOLEAN
prlg_callback_Pmn CB_ARGS(e)
{
    SIMPLEX *e0 = e->parent;
    SIMPLEX *e00 = phgTraverseElem0_;
    DOF *uC = uC_traverse;
    DOF *uF = uF_traverse;
    int Nc = uC->type->nbas;
    int Nf = uF->type->nbas;
    int i, j, k, dim = uC->dim;
    BOOLEAN continuous = TRUE;
    FLOAT matP[Nf*dim][Nc*dim], identy[Nf*dim];
    FLOAT *coord;
    const FLOAT *lambda, *v;
    
    continuous = (uC->type->continuity < 0) ? FALSE : TRUE;
    if (phgTraverseDepth_ == 0) {
    } else {
	assert(e0 != NULL );
	assert(e0->children[0] == e || e0->children[1] == e);
    }

    /* traverse debug */
    for (i = 0; i < phgTraverseDepth_; i++)
	MG_DEBUGn(3, "### ");
    if (IsLeaf(e)) {
	MG_DEBUGn(3, "EEE Traver: %5d\n", e->index);
    } else {
	MG_DEBUGn(3, "### Traver: %5d\n", e->index);
    }
    
    if (!continuous)
	for (i = 0; i < Nf*dim; i++)
	    identy[i] = 1.;

    assert(uF->dim == dim);
    /* Add to global prlgection matrix if leaf node */
    if (IsLeaf(e)) {
	Bzero(matP);
	if (dim == 1) {
	    for (i = 0; i < Nf; i++) 
		If[i] = phgMapE2L(prlg_traverse->rmap, 0, e, i);

	    for (i = 0; i < Nf; i++) {
		MG_DEBUGn(4, "_B%d_", i);
		coord = phgDofGetElementCoordinates(uF, e, i);
		lambda = phgGeomXYZ2Lambda(g_traverse, e00, coord[0], coord[1], coord[2]);
		v = uC->type->BasFuncs(uC, e00, 0, -1, lambda);
		memcpy(matP[i], v, Nc * sizeof(*matP[i]));
		REMOVE_ZERO(matP[i], Nc);
	    }
	    //SHOW_M(matP[0], Nf, Nc);
	    if (continuous) { /* continuous coarse space */
		phgMatSetEntries(prlg_traverse, Nf, If, Nc, Ic, &matP[0][0]); 
	    } else {		       /* discontinuous coarse space */
		phgMatAddEntries(prlg_traverse, Nf, If, Nc, Ic, &matP[0][0]); 
	    }
	} else {
	    for (i = 0; i < Nf; i++) 
		for (k = 0; k < dim; k++)
		    If[i*dim + k] = phgMapE2L(prlg_traverse->rmap, 0, e, i*dim + k);

	    //SHOW_iV(If, Nf*dim);
	    for (i = 0; i < Nf; i++) {
		MG_DEBUGn(4, "_B%d_", i);
		coord = phgDofGetElementCoordinates(uF, e, i*dim);
		lambda = phgGeomXYZ2Lambda(g_traverse, e00, coord[0], coord[1], coord[2]);
		v = uC->type->BasFuncs(uC, e00, 0, -1, lambda);
		//SHOW_V(v, Nc);
		for (k = 0; k < dim; k++) 
		    for (j = 0; j < Nc; j++) 
			matP[i*dim + k][j*dim + k] =  v[j];
		REMOVE_ZERO(matP[i*dim], dim*Nc*dim);
		//SHOW_M(matP[i*dim], dim, Nc*dim);
		if (continuous) { /* continuous coarse space */
		    for (k = 0; k < dim; k++) 
			phgMatSetEntries(prlg_traverse, 
					 1, If+i*dim+k, Nc*dim, Ic, &matP[i*dim+k][0]); 
		} else {	           /* discontinuous coarse space */
		    for (k = 0; k < dim; k++) 
			phgMatAddEntries(prlg_traverse,
					 1, If+i*dim+k, Nc*dim, Ic, &matP[i*dim+k][0]); 
		}
	    }
	}
	if (!continuous)
	    phgVecAddEntries(vec_scale, 0, Nf*dim, If, identy);
    }

    return TRUE;
}

/* TODO: DG prolongation */
#if 0
static BOOLEAN
prlg_callback_DGn CB_ARGS(e)
{
    return TRUE;
}
#endif


void 
build_prolong_mat(MAT *prlg, int level, GRID *gF, GRID *gC, 
		  DOF *uF, DOF *uC) 
{
    int i, k, N_, dim = uF->dim;
    SIMPLEX *e;
    PRLG_RSTR_CALLBACK prlg_callback = NULL;
    BOOLEAN continuous = (uC->type->continuity == -1) ? FALSE : TRUE;
    BOOLEAN on_same_grid = (gF == gC) ? TRUE : FALSE;

    phgInfo(1, "*** Prolongation Coarse space: %s --> Fine space: %s\n", 
	    uC->type->name, uF->type->name);
    if (on_same_grid)
	phgInfo(1, "    prolong on same grid.\n");
    continuous = (uC->type->continuity == -1) ? FALSE : TRUE;
    N_ = uF->type->nbas;
    if (uC->type == DOF_P1 && uF->type == DOF_P1) {
	prlg_callback = prlg_callback_Pn;
	prlg_mat_size = N_*N_;
	phgInfo(1, "    using prlg callback Pn\n");
    } else if (uC->type == DOF_P2 && uF->type == DOF_P2) {
	prlg_callback = prlg_callback_Pmn;
	prlg_mat_size = N_*N_;
	phgInfo(1, "    using prlg callback Pmn\n");
    } else {
	prlg_callback = prlg_callback_Pmn;
	prlg_mat_size = 0;
	phgInfo(1, "    using prlg callback Pmn\n");
    }

    N_ = MAX(N_, uC->type->nbas);
    if (Isize == 0 || Isize < N_) {
	dim = uF->dim;
	assert(dim == uC->dim);
	Ic = phgRealloc0(Ic, N_*dim*sizeof(*Ic));
	If = phgRealloc0(If, N_*dim*sizeof(*If));
	prlg_stack = phgRealloc0(prlg_stack, 
				100*prlg_mat_size * sizeof(*prlg_stack));
	phgFreeAtExit((void **)(void *)&Ic);
	phgFreeAtExit((void **)(void *)&If);
	phgFreeAtExit((void **)(void *)&prlg_stack);
    }

    if (continuous) {		/* continuous coarse space:
				 * use interpolation */
	phgMatSetMode(prlg, PHG_REPLACE);
    } else {			/* discontinuous coarse space:
				 * use weighted interpolation */
	vec_scale = phgMapCreateVec(prlg->rmap, 1);
	phgVecDisassemble(vec_scale);
	phgInfo(1, "    discontinuous coarse space: use weight interpolation\n");
    }

    /* static traverse variables */
    g_traverse = gF;		/* static g_traverse */
    prlg_traverse = prlg;	/* prolong matrix */
    uC_traverse = uC;		/* Dof on coarse grid */
    uF_traverse = uF;		/* Dof on fine grid */

    ForAllElements(gC, e) {
	int Nc = uC->type->nbas;

	bzero(Ic, Nc*dim*sizeof(*Ic));
	MG_DEBUGn(3, "--------\n");
	MG_DEBUGn(3, "elem: %d\n", e->index);

	if (phgVerbosity > 2 && FALSE)
	    phgDumpElement(gC, e);

	if (dim == 1) {
	    for (i = 0; i < Nc; i++)
		Ic[i] = phgMapE2L(prlg->cmap, 0, e, i);
	} else {
	    for (i = 0; i < Nc; i++)
		for (k = 0; k < dim; k++)
		    Ic[i*dim + k] = phgMapE2L(prlg->cmap, 0, e, i*dim + k);
	}
	//SHOW_iV(Ic, dim*Nc);
	phgTraverseDepth_ = 0;
	phgTraverseIndex_ = 0;
	phgTraverseElem0_ = e;
	if (on_same_grid)
	    (*prlg_callback) CB_ARGS0(g_traverse, e);
	else
	    traverse_element_preorder(e, prlg_callback);
	assert(phgTraverseDepth_ == 0);
    }

    phgMatAssemble(prlg_traverse);

    /* discontinuous coarse space, use weight to prolong matrix. */
    if (!continuous) {	
	MAT *mat = prlg;
	MAT_ROW *row;
	FLOAT *vs = vec_scale->data;
	int j;

	phgVecAssemble(vec_scale);
	row = mat->rows;
	for (i = 0; i < mat->rmap->nlocal; i++, row++, vs++)
	    if (Fabs(*vs - 1.) > 1e-12)
		for (j = 0; j < row->ncols; j++)
		    row->data[j] /= *vs;
	if (DUMP_MAT_VEC) {
	    char vec_name[100], vec_file[100];
	    sprintf(vec_name, "U%d", level);
	    sprintf(vec_file, "U%d_.m", level);
	    phgPrintf("*** Dumping U%d...\n", level);
	    phgVecDumpMATLAB(vec_scale, vec_name, vec_file);
	}
	phgVecDestroy(&vec_scale);
    }

    return;
}


/* tranfer mat mat-vec multiplication */
static int
trans_mv_func(MAT_OP op, MAT *mat, VEC *x, VEC *y) 
{
    MG_TRANS_OP *trans = (MG_TRANS_OP *) mat->mv_data[0];
    MAT *mat1 = trans->mat_CtoP;
    MAT *mat2 = trans->mat_PtoF;
    VEC *z = trans->vec_P;

    if (op == MAT_OP_D) {
	phgError(1, "transfer mat op diag not defined!\n");
    } 
    else if (op == MAT_OP_N) {
	phgMatVec(MAT_OP_N, 1., mat1, x, 0., &z);
	phgMatVec(MAT_OP_N, 1., mat2, z, 0., &y);
    }
    else if (op == MAT_OP_T) {
	phgMatVec(MAT_OP_T, 1., mat2, x, 0., &z);
	phgMatVec(MAT_OP_T, 1., mat1, z, 0., &y);
    } 

    return 0;
}


/* build prolongation operator matrix */
void 
phgMultiGridBuildProlongationMat(MULTI_GRID *mg, int level)
{
    GRID *g = mg->grid;
    SUB_GRID *gC;
    MG_LEVEL **ml = mg->ml;
    DOF *uC, *uF;
    MAP *mapC, *mapF;
    int i, islvt, nslvt = ml[level]->nslvt;

    /* No Prolong on level 0 */
    if (level == 0) {
	for (i = 0; i < nslvt; i++)
	    ml[level]->P_[i] = NULL;
	return;
    }

    /* Prolongation should only be build on the same tree. */
    assert(g == ml[level-1]->grid);

    for (islvt = 0; islvt < nslvt; islvt++) {
	set_active_solver(mg, islvt, -1);
	uC = ml[level-1]->u;
	uF = ml[level]->u;
	mapC = ml[level-1]->map;
	mapF = ml[level]->map;
	gC = ml[level-1]->sub_grid;

	/* No CR1 */
	{
	    ml[level]->P = phgMapCreateMat(mapF, mapC);
	    build_prolong_mat(ml[level]->P, level, g, gC,
			      uF, uC);

	    if (DUMP_MAT_VEC) {
		char mat_name[100], mat_file[100];
		sprintf(mat_name, "P%d", level);
		sprintf(mat_file, "P%d_%d_.m", islvt, level);
		phgPrintf("*** Dumping %s\n", mat_file);
		phgMatDumpMATLAB(ml[level]->P, mat_name, mat_file);
	    }
	}
	
	ml[level]->P_[islvt] = ml[level]->P;
	ml[level]->P = NULL;
    }
    set_active_solver(mg, -1, -1);

    return;
}

/*
 * Local Restriction implement:
 *
 * Here we use a mat-vec way, which have two advantages:
 *   1. as the multigrid is fixed, the same restriction process occurs every
 *      time at the begining of new time step, a mat-vec way would save a lot
 *      of time compare to bisection-coarsen-interpolation way
 *   2. the level different between antecedents and decendents may be 2 or 3,
 *      I'm not sure how to do it bisection-coarsen-interpolation way
 * 
 * The main target is to locate the bases of antecedents element in which
 * decendent element. We use a list to keep the progress of each basis, when
 * we travel on new element, the list is checked, the basis in the element will
 * find the children it lies in using lambda coordinates, and fresh the list
 * accordingly. In the end, all the bases find the decendent elements it locates.
 * 
 */

typedef struct BASIS_LOCATION_ {
    /* int i;			/\* basis no on coarse elements *\/ */
    /* INT I;			/\* local index of the basis *\/ */
    FLOAT coord[3];		/* coordinates of the basis */
    SIMPLEX *e;			/* current looking up element */
} BASIS_LOCATION; 

static BASIS_LOCATION *bas_loc;	/* dim-less */


/*
 * Restriction from fine CONTINUOUS space to coarse space. 
 * Use interpolation, and locate one coarse space basis in
 * one fine grid element is enough.
 *
 *  */
static BOOLEAN
rstr_callback_continuous CB_ARGS(e)
{
    SIMPLEX *e0 = e->parent;
    DOF *uC = uC_traverse;
    DOF *uF = uF_traverse;
    int Nc = uC->type->nbas;
    int Nf = uF->type->nbas;
    int i, j, k, dim = uC->dim;
    FLOAT matR[Nc*dim][Nf*dim];
    FLOAT *coord;
    const FLOAT *lambda, *v;
    
    assert(uF->type->continuity > -1);
    if (phgTraverseDepth_ == 0) {
    } else {
	assert(e0 != NULL );
	assert(e0->children[0] == e || e0->children[1] == e);
    }

    /* traverse debug */
    for (i = 0; i < phgTraverseDepth_; i++)
	MG_DEBUGn(3, "### ");
    if (IsLeaf(e)) {
	MG_DEBUGn(3, "EEE %d\n", e->index);
    } else {
	MG_DEBUGn(3, "### %d(%d,%d)\n", e->index, 
		  e->children[0]->index, e->children[1]->index);
    }

    assert(dim == uC->dim);
    /* Add to global restriction matrix if leaf node */
    if (IsLeaf(e)) {
	Bzero(matR);
	if (dim == 1) {
	    for (i = 0; i < Nf; i++) 
		If[i] = phgMapE2L(rstr_traverse->cmap, 0, e, i);

	    for (i = 0; i < Nc; i++) {
		if (bas_loc[i].e != e)
		    continue;
	    
		bas_loc[i].e = NULL;  /* locate process end */
		coord = bas_loc[i].coord;
		lambda = phgGeomXYZ2Lambda(g_traverse, e,
					   coord[0], coord[1], coord[2]);
		v = uF->type->BasFuncs(uF, e, 0, -1, lambda);
		memcpy(matR[i], v, Nf * sizeof(*matR[i]));
		REMOVE_ZERO(matR[i], Nf);
		phgMatSetEntries(rstr_traverse,
				 1, Ic+i, Nf, If, &matR[i][0]); 
	    }
	} else {
	    for (i = 0; i < Nf; i++) 
		for (k = 0; k < dim; k++)
		    If[i*dim + k] = phgMapE2L(rstr_traverse->cmap, 0, e, i*dim + k);

	    for (i = 0; i < Nc; i++) {
		if (bas_loc[i].e != e)
		    continue;
	    
		//SHOW_iV(If, Nf*dim);
		bas_loc[i].e = NULL;  /* locate process end */
		coord = bas_loc[i].coord;
		lambda = phgGeomXYZ2Lambda(g_traverse, e,
					   coord[0], coord[1], coord[2]);
		v = uF->type->BasFuncs(uF, e, 0, -1, lambda);
		for (k = 0; k < dim; k++) 
		    for (j = 0; j < Nf; j++) 
			matR[i*dim + k][j*dim + k] = v[j];
		REMOVE_ZERO(matR[i*dim], dim*Nf*dim);
		for (k = 0; k < dim; k++) 
		    phgMatSetEntries(rstr_traverse,
				     1, Ic+i*dim+k, Nf*dim, If, &matR[i*dim+k][0]); 
		//SHOW_M(matR[i*dim], dim, Nf*dim);
	    }
	}
    } else {
	for (i = 0; i < Nc; i++) {
	    if (bas_loc[i].e != e)
		continue;

	    coord = bas_loc[i].coord;
	    /* find the bases in which children of the element */
#if 0
	    /* test both children */
	    for (k = 0; k < 2; k++) {
		lambda = phgGeomXYZ2Lambda(g_traverse, e->children[k], 
					   coord[0], coord[1], coord[2]);
		if (LambdaInElement(lambda)) {
		    bas_loc[i].e = e->children[k];
		    MG_DEBUGn(4, "   bas:%d: %d --> %d\n", i, e->index,
			      e->children[k]->index);
		    break;
		}
	    }
	    assert(k < 2); 
#else
	    /* test one children */
	    MG_DEBUGn(4, "_B%d_", i);
	    lambda = phgGeomXYZ2Lambda(g_traverse, e->children[0], 
				       coord[0], coord[1], coord[2]);
	    if (LambdaInElement(lambda)) {
		k = 0;
		bas_loc[i].e = e->children[0];
	    } else {
		k = 1;
		bas_loc[i].e = e->children[1];
	    }
	    MG_DEBUGn(4, "   bas:%d: %d --> %d\n", i, e->index,
		      e->children[k]->index);
#endif
	}
    }

    return TRUE;
}


/*
 * Restriction from fine DISCONTINUOUS space to coarse space. 
 * Use weighted interpolation, and locate one coarse space basis in
 * ALL fine grid element which contains it. 
 *
 *  */
static BOOLEAN
rstr_callback_discontinuous CB_ARGS(e)
{
    SIMPLEX *e0 = e->parent;
    DOF *uC = uC_traverse;
    DOF *uF = uF_traverse;
    int Nc = uC->type->nbas;
    int Nf = uF->type->nbas;
    int i, j, k, dim = uC->dim;
    FLOAT matR[Nc*dim][Nf*dim];
    FLOAT *coord, identy[dim];
    const FLOAT *lambda, *v;
    
    assert(uF->type->continuity == -1);
    if (phgTraverseDepth_ == 0) {
    } else {
	assert(e0 != NULL );
	assert(e0->children[0] == e || e0->children[1] == e);
    }

    /* traverse debug */
    for (i = 0; i < phgTraverseDepth_; i++)
	MG_DEBUGn(3, "### ");
    if (IsLeaf(e)) {
	MG_DEBUGn(3, "EEE %d\n", e->index);
    } else {
	MG_DEBUGn(3, "### %d(%d,%d)\n", e->index, 
		  e->children[0]->index, e->children[1]->index);
    }

    for (k = 0; k < dim; k++)
	identy[k] = 1.;

    assert(dim == uC->dim);
    /* Add to global restriction matrix if leaf node */
    if (IsLeaf(e)) {
	Bzero(matR);
	if (dim == 1) {
	    for (i = 0; i < Nf; i++) 
		If[i] = phgMapE2L(rstr_traverse->cmap, 0, e, i);

	    for (i = 0; i < Nc; i++) {
		coord = bas_loc[i].coord;
		lambda = phgGeomXYZ2Lambda(g_traverse, e,
					   coord[0], coord[1], coord[2]);
		if (!LambdaInElement(lambda)) 
		    continue;
		
		v = uF->type->BasFuncs(uF, e, 0, -1, lambda);
		memcpy(matR[i], v, Nf * sizeof(*matR[i]));
		REMOVE_ZERO(matR[i], Nf);
		phgMatAddEntries(rstr_traverse,
				  1, Ic+i, Nf, If, &matR[i][0]);
		phgVecAddEntries(vec_scale, 0, 1, Ic+i, identy);
	    }
	} else {
	    for (i = 0; i < Nf; i++) 
		for (k = 0; k < dim; k++)
		    If[i*dim + k] = phgMapE2L(rstr_traverse->cmap, 0, e, i*dim + k);

	    for (i = 0; i < Nc; i++) {
		coord = bas_loc[i].coord;
		lambda = phgGeomXYZ2Lambda(g_traverse, e,
					   coord[0], coord[1], coord[2]);
		if (!LambdaInElement(lambda)) 
		    continue;

		v = uF->type->BasFuncs(uF, e, 0, -1, lambda);
		for (k = 0; k < dim; k++) 
		    for (j = 0; j < Nf; j++) 
			matR[i*dim + k][j*dim + k] = v[j];
		REMOVE_ZERO(matR[i*dim], dim*Nf*dim);
		for (k = 0; k < dim; k++) 
		    phgMatAddEntries(rstr_traverse, 
				     1, Ic+i*dim+k, Nf*dim, If, &matR[i*dim+k][0]); 
		phgVecAddEntries(vec_scale, 0, dim, Ic+i*dim, identy);
		//SHOW_M(matR[i*dim], dim, Nf*dim);
	    }
	}
    }

    return TRUE;
}

/* build restriction operator matrix */
void 
phgMultiGridBuildRestrictionMat(MULTI_GRID *mg, int level)
{
    GRID *g = mg->grid;
    MG_LEVEL **ml = mg->ml;
    SIMPLEX *e;
    DOF *uC, *uF; 
    int i, k, dim, N_, islvt, nslvt = ml[level]->nslvt;
    BOOLEAN continuous = TRUE;
    PRLG_RSTR_CALLBACK rstr_callback = NULL;

    /* No restriction on level 0 */
    if (level == 0) {
	for (i = 0; i < nslvt; i++)
	    ml[level]->R_[i] = NULL;
	return;
    }
    
    /* Restriction should only be build on the same tree. */
    assert(g == ml[level]->grid);

    for (islvt = 0; islvt < nslvt; islvt++) {
	set_active_solver(mg, islvt, -1);
	uC = ml[level-1]->u;
	uF = ml[level]->u;
	phgInfo(1, "*** Restriction Fine space: %s --> Coarse space: %s\n", 
		uF->type->name, uC->type->name);
	continuous = (uF->type->continuity == -1) ? FALSE : TRUE;
	if (continuous) {		/* continuous fine space */
	    rstr_callback = rstr_callback_continuous;
	    phgInfo(1, "   use restriction callback continuous\n");
	} else {			/* discontinuous fine space */
	    rstr_callback = rstr_callback_discontinuous;
	    phgInfo(1, "   use restriction callback discontinuous\n");
	}
	N_ = MAX(uF->type->nbas,
		 uC->type->nbas);
	dim = uF->dim;
	if (Isize == 0 || Isize < N_) {
	    assert(dim == ml[level-1]->u->dim);
	    Ic = phgRealloc0(Ic, N_*dim*sizeof(*Ic));
	    If = phgRealloc0(If, N_*dim*sizeof(*If));
	    bas_loc = phgRealloc0(bas_loc, N_*sizeof(*bas_loc));
	    phgFreeAtExit((void **)(void *)&Ic);
	    phgFreeAtExit((void **)(void *)&If);
	    phgFreeAtExit((void **)(void *)&bas_loc);
	}

	ml[level]->R = phgMapCreateMat(ml[level-1]->map, ml[level]->map);

	if (continuous) {		/* continuous fine space
					 * use interpolation*/
	    phgMatSetMode(ml[level]->R, PHG_REPLACE);
	} else {			/* discontinuous fine space
					 * use weighted interpolation*/ 
	    vec_scale = phgMapCreateVec(ml[level-1]->map, 1);
	    phgVecDisassemble(vec_scale);
	    phgInfo(1, "discontinuous Fine space: use weight interpolation\n");
	}

	/* static traverse variables */
	g_traverse = g;		/* static g_traverse */
	//level_traverse = level;	/* static ml_level */
	uC_traverse = uC;
	uF_traverse = uF;
	rstr_traverse = ml[level]->R;

	ForAllElements(ml[level - 1]->sub_grid, e) {
	    FLOAT *coord;
	    int Nc = uC->type->nbas;

	    bzero(Ic, Nc*dim*sizeof(*Ic));
	    MG_DEBUGn(3, "--------\n");
	    MG_DEBUGn(3, "elem: %d\n", e->index);

	    if (phgVerbosity > 2 && FALSE)
		phgDumpElement(g, e);

	    for (i = 0; i < Nc; i++) {
		for (k = 0; k < dim; k++)
		    Ic[i*dim + k] = phgMapE2L(rstr_traverse->rmap, 0, e, i*dim + k);

		/* bas_loc is dim-less */
		coord = phgDofGetElementCoordinates(uC, e, i*dim);
		memcpy(bas_loc[i].coord, coord, Dim * sizeof(*coord));
		bas_loc[i].e = e;
	    }
	    //SHOW_iV(Ic, Nc*dim);
	
	    phgTraverseDepth_ = 0;
	    phgTraverseIndex_ = 0;
	    phgTraverseElem0_ = e;
	    traverse_element_preorder(e, rstr_callback);

	    if (continuous)
		for (i = 0; i < Nc; i++)
		    assert(bas_loc[i].e == NULL);
	    assert(phgTraverseDepth_ == 0);
	}
	phgMatAssemble(rstr_traverse);


	if (!continuous) {
	    MAT *mat = rstr_traverse;
	    MAT_ROW *row;
	    FLOAT *vs = vec_scale->data;
	    int j;

	    phgVecAssemble(vec_scale);
	    row = mat->rows;
	    for (i = 0; i < mat->rmap->nlocal; i++, row++, vs++)
		if (Fabs(*vs - 1.) > 1e-12)
		    for (j = 0; j < row->ncols; j++)
			row->data[j] /= *vs;
	    if (DUMP_MAT_VEC) {
		char vec_name[100], vec_file[100];
		sprintf(vec_name, "S%d", level);
		sprintf(vec_file, "S%d_.m", level);
		phgPrintf("*** Dumping S%d...\n", level);
		phgVecDumpMATLAB(vec_scale, vec_name, vec_file);
	    }
	    phgVecDestroy(&vec_scale);
	}

	if (DUMP_MAT_VEC) {
	    char mat_name[100], mat_file[100];
	    sprintf(mat_name, "R%d", level);
	    sprintf(mat_file, "R%d_%d_.m", islvt, level);
	    phgPrintf("*** Dumping %s\n", mat_file);
	    phgMatDumpMATLAB(ml[level]->R, mat_name, mat_file);
	}

	ml[level]->R_[islvt] = ml[level]->R;
	ml[level]->R = NULL;
    }
    set_active_solver(mg, -1, -1);

    return;
}


#if 1
# define MG_CYCLE_DEBUG(...) {					\
	int _i;							\
	MULTI_GRID *mg = ml[0]->mg;				\
	phgInfo(0, "Solver: %d[%d]   ", mg->active_solver,	\
		mg->active_solver_type);			\
	for (_i = 0; _i < _mgp->max_level - level; _i++)	\
	    phgInfo(0, "##### ");				\
	phgInfo(0, " L%d ", level);				\
	phgInfo(0, __VA_ARGS__);				\
	phgInfo(0, "\n");					\
    }
#else
# define MG_CYCLE_DEBUG(...)
#endif


static double time0 = 0., time1 = 0.;
#define COARSEST_SMOOTH 1

/*
 * Recursive iteration of cycling.
 * TODO: Non-recursive version of cycling in Boomeramg?
 * 
 * */
void
recursive_MG_cycle(MULTI_GRID *mg, int level, BOOLEAN init_zero)
{
    MG_LEVEL **ml = mg->ml;
    GRID *g = ml[level]->grid;
    int islv = mg->active_solver;
    MG_SMOOTHER smoother;
    MAT *A = ml[level]->mat;
    MAT *P = ml[level]->P;
    VEC *x = ml[level]->x;
    VEC *f = ml[level]->f;
    VEC *r = ml[level]->r;
    int i;

    if (level == _mgp->max_level - 1) {
	phgPrintf("_mg_ ");
	MPI_Barrier(g->comm);
        time1 = phgGetTime(NULL);             
    }

    if (level < mg->start_level) {
	MG_CYCLE_DEBUG("~");
	return;
    }

    if (init_zero)
	ZeroVec(ml[level]->x);

    /*
     * coarsest level solve
     * */
    if (level <= _mgp->min_level) { 
	assert(level >= 0);

	if (_mgp->solve_coarst) {
	    /* use solver on coarsest grid */
	    SOLVER *solver = ml[level]->solver;
	    FLOAT *rhs = ml[level]->solver->rhs->data;
    
#if COARSEST_SMOOTH
	    smoother = phgMultiGridGetSmoother(mg, level, COARSEST);
	    TIMING_BEGIN;
	    smoother(A, x, f, _mgp->n_pre_smooth, (void *)ml[level]);
	    if (smoother == mg_GaussSidel_line) {
		MG_CYCLE_DEBUG("\\L[%d]", _mgp->n_pre_smooth);
		TIMING_END("Pre smoothL");
	    } else {
		MG_CYCLE_DEBUG("\\[%d]", _mgp->n_pre_smooth);
		TIMING_END("Pre smooth");
	    }
#endif

	    solver->rhs_updated = TRUE;
	    solver->rhs->data = f->data;
	    TIMING_BEGIN;
	    phgSolverVecSolve(solver, FALSE, x);
	    TIMING_END("Coarst solve");
	    solver->rhs->data = rhs;
	    phgInfo(1, "   use user solver on coarst level\n");
	    MG_CYCLE_DEBUG("!!!");
	    phgInfo(0, "      coasrt solve nits:%d, res:%e\n",
		    solver->nits, solver->residual);

#if COARSEST_SMOOTH
	    TIMING_BEGIN;
	    smoother(A, x, f, _mgp->n_post_smooth, (void *)ml[level]);
	    if (smoother == mg_GaussSidel_line) {
		MG_CYCLE_DEBUG("\\L[%d]", _mgp->n_post_smooth);
		TIMING_END("Post smoothL");
	    } else {
		MG_CYCLE_DEBUG("\\[%d]", _mgp->n_post_smooth);
		TIMING_END("Post smooth");
	    }
#endif
	} else {
	    /* use smoother on coarsest grid */
	    smoother = phgMultiGridGetSmoother(mg, level, COARSEST);
	    TIMING_BEGIN;
	    smoother(A, x, f, _mgp->n_coarst_smooth, (void *)ml[level]);
	    TIMING_END("Coarst smooth");
	    MG_CYCLE_DEBUG("!");
	}
    }
    /*
     * redundent level 
     * */
    else if (mg->ml[level]->redunt_level) {
	/* just copy */

	MG_CYCLE_DEBUG("-");
	phgVecCopy(f, &ml[level-1]->f);

	recursive_MG_cycle(mg, level-1, init_zero);

	phgVecCopy(ml[level-1]->x, &x);
	MG_CYCLE_DEBUG("-");
    }
    /*
     * intermeidate level solve
     * */
    else {
	/* pre-smooth */
	smoother = phgMultiGridGetSmoother(mg, level, DOWN_CYCLE);
	TIMING_BEGIN;
	smoother(A, x, f, _mgp->n_pre_smooth, (void *)ml[level]);
	if (smoother == mg_GaussSidel_line) {
	    MG_CYCLE_DEBUG("\\L[%d]", _mgp->n_pre_smooth);
	    TIMING_END("Pre smoothL");
	} else {
	    MG_CYCLE_DEBUG("\\[%d]", _mgp->n_pre_smooth);
	    TIMING_END("Pre smooth");
	}

	/* residual:
	 *   r_{k} = b_{k} - A_{k} * x_{k} */
	phgVecCopy(f, &r);
	TIMING_BEGIN;
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
	TIMING_END("Residual");

	/* restriction:
	 *   f_{k-1} = P_{k} r_{k} */
	TIMING_BEGIN;
	phgMatVec(MAT_OP_T, 1.0, P, r, 0., &ml[level-1]->f);
	TIMING_END("restrict");
	
	/* recursive cycle
	 * n_cycle = 1: V-cycle
	 * n_cycle = 2: W-cycle
	 */
	for (i = 0; i < _mgp->n_cycle; i++) {
	    /* Note:
	     *   Not sure in smoothing as done in Alberta is neccessory.
	     *   Without this in smoothing, W(+) cycle may not improve none
	     *   compare to V cycle.
	     * */
	    recursive_MG_cycle(mg, level - 1, (i==0));
	}

	/* prolongation:
	 *   x_{k} += P_{k}^T x_{k-1} */
	TIMING_BEGIN;
	phgMatVec(MAT_OP_N, _mgp->correct_damp, 
		  P, ml[level-1]->x, 1.0, &x);
	TIMING_END("prolong");

	TIMING_BEGIN;
	smoother = phgMultiGridGetSmoother(mg, level, UP_CYCLE);
	smoother(A, x, f, _mgp->n_post_smooth, (void *)ml[level]);
	if (smoother == mg_GaussSidel_line) {
	    MG_CYCLE_DEBUG("/L[%d]", _mgp->n_post_smooth);
	    TIMING_END("Post smoothL");
	} else {
	    MG_CYCLE_DEBUG("/[%d]", _mgp->n_post_smooth);
	    TIMING_END("Post smooth");
	}
    }

    if (level == _mgp->max_level - 1) {
	MPI_Barrier(g->comm);
        if (g->rank == 0)                               
	    phgInfo(0, "   Cycle cost: %0.8lfs\n",  
		    phgGetTime(NULL) - time1);
    }
    return;
}


void
phgMultiGridSolve(MULTI_GRID *mg, SOLVER *solver, DOF *u)
{
    MG_LEVEL **ml = mg->ml;
    FLOAT r0_norm, r_norm, b_norm;
    int nits = 0, level, max_level = _mgp->max_level;
    MAT *A = solver->mat;
    VEC *x, *r, *b0 = solver->rhs;


    /* Show MG info */
    if (phgVerbosity > 0) {
	phgPrintf("\n*********************\n");
	phgPrintf("MultiGrid Info:\n");
	phgPrintf("   levels: %d\n", max_level);
	for (level = mg->start_level; level < max_level; level++) {
	    phgPrintf("      level: %d, DOF: %d, elems: %d\n",
		       level, DofGetDataCountGlobal(ml[level]->u), 
		       ml[level]->sub_grid->nleaf_global); 
	}
	phgPrintf("   damping : %e\n", _mgp->smooth_damp);
	phgPrintf("   # pre  smooth: %d\n", _mgp->n_pre_smooth);
	phgPrintf("   # post smooth: %d\n", _mgp->n_post_smooth);
	phgPrintf("\n");
    }

    /* work on finest level */
    level = _mgp->max_level - 1; 
    x = ml[level]->x;
    r = phgMapCreateVec(solver->rhs->map, 1);
    phgMapDofToLocalData(solver->rhs->map, 1, &u, x->data);


    /* initial residual */
    r_norm = 0.;
    b_norm = phgVecNorm2(b0, 0, NULL);
    phgVecCopy(b0, &r);
    phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r); 	
    r0_norm = phgVecNorm2(r, 0, NULL);
    phgPrintf("   initial err %E\n", r0_norm);
    phgVecDumpMATLAB(r, "r0", "r0_.m");

    /* MG cycle */
    while (nits < _mgp->maxit) {
	/* MG cycle setup */
	phgVecCopy(b0, &ml[level]->f);

	/* MG correct on x */
	recursive_MG_cycle(mg, level, FALSE);
	    
	/* new residual */
	phgVecCopy(b0, &r);
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r); 
	r_norm = phgVecNorm2(r, 0, NULL);
	phgPrintf("   niter: %d, err %E\n", 
		   nits, r_norm / r0_norm);

	nits++;
	if ((r_norm / r0_norm) < _mgp->rtol)
	    break;
    }

    solver->nits = nits;
    solver->residual = r_norm;
    if (nits < _mgp->maxit) {
	phgPrintf("   MG converge at step: %d\n", nits);
    } else
	phgPrintf("WARNING: maxit attained in OEM solver, "
		   "the result may be inaccurate.\n");

    if (0 && DUMP_MAT_VEC) {
	VEC *solu = phgMapCreateVec(solver->rhs->map, 1);
	phgMapDofToLocalData(solver->rhs->map, 1, &u, solu->data);
	phgVecDumpMATLAB(solu, "x0", "x0_.m");
	phgVecDumpMATLAB(x, "x", "x_.m");
	phgVecDestroy(&solu);
    }

    phgMapLocalDataToDof(solver->rhs->map, 1, &u, x->data);
    phgVecDestroy(&r);

    return;
}
