/*
 * Moving mesh routines, get code from AFEPack MovingMesh3D.cpp.
 * Note: this is a constrained boundary version.
 *  */
#include "moving_mesh.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#define MOVINTMESH_DEBUG

#define _p params
#define _m mmesh
#define _mp (mmesh->params)
#define _mb _m
#define _mv (mmesh->info_vert)

static int viz_step = 0;
static char viz_name[100];
static VEC *vec_reduce = NULL;
static MAP *map_reduce = NULL;
static FLOAT _rhs[NVert * Dim];
static INT _I[NVert * Dim];

static BTYPE btype_values[] = {
    DIRICHLET, NEUMANN, 
#ifdef EXTEND_BDRY_TYPE 
    BDRY_USER1, BDRY_USER2, BDRY_USER3, BDRY_USER4, BDRY_USER5,
    BDRY_USER6, BDRY_USER7, BDRY_USER8, BDRY_USER9, BDRY_USER10,
#endif /* EXTEND_BDRY_TYPE */
    UNDEFINED
};

static void getMoveStepLength(MovingMesh *mmesh);
static void getMoveDirection(MovingMesh *mmesh);
static void updateMesh(MovingMesh *mmesh);
static void updateDof(MovingMesh *mmesh);

static double
elapsed_time(GRID *g, BOOLEAN flag, double mflops)
/* returns and prints elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;
 
    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    mem = phgMemoryUsage(g, NULL);

    if (flag) {
	if (mflops <= 0)
	    phgPrintf("[%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
	else
	    phgPrintf("[%0.4lgMB %0.4lfs %0.4lgGF]\n", mem / (1024.0 * 1024.0),
		      et, mflops*1e-3);
    }

    return et;
}

char *
uchar2bit(INT n)
{
    int i;
    static char bit_str[100];

    bit_str[0] = '@';
    for (i = 8; i > 0; i--, n >>= 1)
	bit_str[i] = (n % 2 == 0) ? 'o' : '1';
    bit_str[9] = '\0';

    return bit_str;
}

void
phgDofMaxMin(DOF *u, FLOAT *umax, FLOAT *umin)
{
    GRID *g = u->g;
    int i, n = DofGetDataCount(u);
    FLOAT *p = u->data;
    FLOAT max = -1e10;
    FLOAT min = 1e10;
    FLOAT gmin, gmax;
    
    for (i = 0; i < n; i++) {
	max = MAX(max, *p);
	min = MIN(min, *p);
	p++;
    }
    
#if USE_MPI
    MPI_Reduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
    MPI_Reduce(&min, &gmin, 1, MPI_DOUBLE, MPI_MIN, 0, g->comm);
#endif	/* USE_MPI */

    *umax = gmax;
    *umin = gmin;
    return;
}


PARAMS 
*phgParametersCreate()
{
    PARAMS *params = (PARAMS *) phgAlloc(sizeof(*params));
    
    /* default settings */
    _p->fn = "cube.bdface.mesh";
    _p->fn_bdry = "cube.bdry.dat";

    _p->fix_edge = FALSE;
    _p->fix_bdry = FALSE;

    _p->pre_refines = 0;
    _p->tol = 1e-2;
    _p->n_move_substep = 1;
    _p->n_smooth_monitor = 1;

    _p->move_opts = NULL;
    _p->dof_opts = NULL;

    _p->verb = 0;
    _p->viz_move = FALSE;

    /* get user defined parameter from file*/
    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&_p->fn);
    phgOptionsRegisterFilename("bdry_face_file", "Boundary face file", (char **)&_p->fn_bdry);
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &_p->pre_refines);
    
    phgOptionsRegisterInt("mm_n_move_substep", "Moving Mesh: num of substep in one step", &_p->n_move_substep);
    phgOptionsRegisterInt("mm_n_smooth_monitor", "Moving Mesh: num of smoothing mointor", &_p->n_smooth_monitor);
    phgOptionsRegisterInt("mm_verbosity", "Moving Mesh: Verbosity", &_p->verb);
    phgOptionsRegisterFloat("mm_tol", "Moving Mesh: tolerence of error", &_p->tol);
    phgOptionsRegisterNoArg("mm_fix_edge", "Moving Mesh: fix edge boundary nodes", &_p->fix_edge);
    phgOptionsRegisterNoArg("mm_fix_bdry", "Moving Mesh: fix all boundary nodes", &_p->fix_bdry);
    phgOptionsRegisterNoArg("mm_viz_move", "Moving Mesh: visualize nodes move", &_p->viz_move);
    phgOptionsRegisterString("mm_move_opts", "Moving Mesh Solver of moving options", &_p->move_opts);
    phgOptionsRegisterString("mm_dof_opts", "Moving Mesh Solver of dof updating options", &_p->dof_opts);
    
    _p->fix_edge |= _p->fix_bdry;
    return params;
}


static BOOLEAN
get_token(FILE *fp, char *token)
{
    int c;
    char *p;

    while (TRUE) {
	if (fscanf(fp, "%s", token) != 1)
	    return FALSE;
	if (token[0] != '#')
	    break;
	/* skip to newline */
	do
	    if ((c = fgetc(fp)) == EOF)
		return FALSE;
	while (c != '\n');
    }
    if ((p = strchr(token, '#')) != NULL)
	*p = '\0';
    return TRUE;
}

static void 
getBdryTypes(MovingMesh *mmesh) 
{
    GRID *g = _m->g;
    SIMPLEX *e;
    FILE *fp;
    char token[1024]; 
    int i, n, s;
    int max_bdry_type = sizeof(BTYPE) * 8;
    
    Unused(s);
    Unused(e);
    Unused(g);

    /* Read in boundary faces info from file */
    fp = fopen(_mp->fn_bdry, "r");

    /* Header */
    if (!get_token(fp, token) || strcasecmp(token, "BoundaryFaces")) 
	READ_ERROR;

    if (!get_token(fp, token))
	READ_ERROR;
    n = atoi(token);
    phgInfo(2, "number of bdry faces: %d\n", n);
    _mb->move_bdry_mask = 0;
    _mb->n_bdry = max_bdry_type;
    _mb->normal = phgCalloc(max_bdry_type * Dim, sizeof(*_mb->normal));
    _mb->b = phgCalloc(max_bdry_type, sizeof(*_mb->b));
    
    /* Read in noraml and projection */
    for (i = 0; i < n; i++) {
	int type;
	READ_NUMBER;
	type = atoi(token) - 1;
	if (!get_token(fp, token))
	    READ_ERROR;
	_mb->normal[type*Dim + 0] = atof(token);
	if (!get_token(fp, token))
	    READ_ERROR;
	_mb->normal[type*Dim + 1] = atof(token);
	if (!get_token(fp, token))
	    READ_ERROR;
	_mb->normal[type*Dim + 2] = atof(token);
	if (!get_token(fp, token))
	    READ_ERROR;
	_mb->b[type] = atof(token);
	_mb->move_bdry_mask |= btype_values[type];
    }

    /* End of read */
    if (!get_token(fp, token) || strcasecmp(token, "End"))
	READ_ERROR;
    fclose(fp);

    return;
}

/*
 *  Moving mesh initialization, following jobs are done:
 *  1. create interior vertex map and boundary vertex map,
 *  2. create a logical mesh.
 *  */
MovingMesh *
phgMovingMeshCreate(GRID *g, PARAMS *params)
{
    MovingMesh *mmesh;
    INT i, k;
    FLOAT *v;
    
    _m = phgCalloc(1, sizeof(*_m));
    _m->params = params;
    _m->g = g;

    _m->logical_node = 
	phgDofNew(g, DOF_P1, 3, "logical node coordinates", DofNoAction);
    _m->logical_node->DB_mask[0] = MM_CONSTR0;
    _m->logical_node->DB_mask[1] = MM_CONSTR1;
    _m->logical_node->DB_mask[2] = MM_CONSTR2;

    _m->logical_move = 
	phgDofNew(g, DOF_P1, 3, "logical move direction", DofNoAction);
    _m->move = 
	phgDofNew(g, DOF_P1, 3, "node move direction", DofNoAction);
    _m->monitor = 
	phgDofNew(g, DOF_P0, 1, "monitor", DofNoAction);

    /* P1 bases */
    _m->phi = 
	phgDofNew(g, DOF_P1, 1, "P1 bases", DofNoAction);
    _m->phi->DB_mask[0] = 0;
    _m->map = phgMapCreate(_m->phi, NULL);

    /* Boundary types of verties */
    getBdryTypes(_m);

    /* get logical mesh
     * Note: here logical mesh is identical to original mesh*/
    v = _m->logical_node->data;
    for (i = 0; i < g->nvert; i++) {
	if (!(g->types_vert[i] == UNREFERENCED)) {
	    for (k = 0; k < Dim; k++)
		v[i*Dim + k] = g->verts[i][k];
	}
    } /* Point loop: local */
    DOF_SCALE(_m->logical_node);

    return _m;
}

void 
phgMovingMeshMove(MovingMesh *mmesh, int max_step)
{
    GRID *g = _m->g;
    FLOAT eps = _mp->tol;
    FLOAT min_mv = 2 * eps, max_mv = 0.;
    INT i, j, k, step;
    FLOAT *v;
    char name[100];
    
    Unused(j);
    Unused(k);

    step = 0; 
    while (min_mv > eps) {
	if (max_step > 0 && step >= max_step) {
	    phgPrintf(MESSAGE_HEADER1"Moving mesh stop: reach max step\n");
	    break;
	}
	elapsed_time(g, FALSE, 0.);	/* reset timer */

	if (_mp->verb > 0) {
	    phgPrintf(MESSAGE_HEADER1"Moving mesh substep: %d ", step);
	    elapsed_time(g, TRUE, 0.);
	}

	getMoveDirection(mmesh); 
	max_mv = 0;
	
	v = _m->move->data;
	for (i = 0; i < g->nvert; i++) {
	    FLOAT d = 0.;

	    d += *v * *v; v++;
	    d += *v * *v; v++;
	    d += *v * *v; v++;
	    d = sqrt(d);
	    min_mv = MIN(min_mv, d);
	    max_mv = MAX(max_mv, d);
	}

#if USE_MPI
	{
	    FLOAT min_mv0 = min_mv, max_mv0 = max_mv;
	    MPI_Allreduce(&min_mv0, &min_mv,
			  1, MPI_DOUBLE, MPI_MIN, g->comm);
	    MPI_Allreduce(&max_mv0, &max_mv,
			  1, MPI_DOUBLE, MPI_MAX, g->comm);
	}
#endif /* USE_MPI */

	if (_mp->verb > 0) {
	    phgPrintf(MESSAGE_HEADER1"mesh moving min move = %e\n", min_mv);
	    phgPrintf(MESSAGE_HEADER1"mesh moving max move = %e\n", max_mv);
	}
	getMoveStepLength(_m);

	for (i = 0;i < _mp->n_move_substep;i ++) {
	    updateDof(_m);
	    updateMesh(_m);
#if 1
	    /* use analytic dof to check update dof */
	    {
		DOF **dofs = _m->dofs;
		//char name[100];
		//static int n = 0;
		DOF *dof_err = phgDofCopy(dofs[0], NULL, NULL, "dof err");
		assert(dofs[1]->type == DOF_ANALYTIC);
		phgDofAXPY(-1.0, dofs[1], &dof_err);
		phgPrintf("* update dofs change = %e\n", phgDofNormL2(dof_err));
		//sprintf(name, "dof_error2_%03d.vtk", n++);
		//phgExportVTK(g, name, dof_err, NULL);
		phgDofFree(&dof_err);
	    }
#endif	/* Check dof update */
	}

#if 1
	sprintf(name, "Moving_mesh.dof_%03d.plt", step);
	//phgExportVTK(g, name, _m->dofs[0], NULL);
	//phgExportTecplot(g, name, _m->dofs[0], NULL);

	//phgExportEnsight(g, "Moving", (double)step, _m->monitor, _m->dofs[0], NULL);
#endif	/* export dof to VTK files */
	step++;
    }

    return;
}


char *
ushort2bit(USHORT n)
{
    int i;
    static char bit_str[100];

    bit_str[0] = '@';
    for (i = sizeof(USHORT)*8; i > 0; i--, n >>= 1)
	bit_str[i] = (n % 2 == 0) ? 'o' : '1';
    bit_str[sizeof(USHORT)*8+1] = '\0';

    return bit_str;
}

/* Get local bases for boundary vertex i,
 *  return num of constrained bases, 
 *  and noramlized direction of rotated bases.
 *  */
BYTE 
getBdryVertBases(MovingMesh *mmesh, INT i, const BYTE **bdry, const FLOAT **bases, const FLOAT **bb)
{
    GRID *g = mmesh->g;
    int j, k, m;
    FLOAT norm;
    static BYTE n[Dim];
    static FLOAT H[Dim][Dim], v[Dim], c[Dim], b[Dim],
	bxyz[Dim][Dim] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    bzero(b, sizeof(b));
    bzero(n, sizeof(n));
    *bases = &bxyz[0][0];
    if (_mp->fix_bdry)
	return 3;

    //printf("* vert type: %s\n", ushort2bit(g->types_vert[i])); 
    for (j = 0, k = 0; j < sizeof(btype_values) / sizeof(btype_values[0]) - 4; j++) {
	if (g->types_vert[i] & btype_values[j])
	    n[k++] = j;
    }

    /* More than three constrains, then fixed */
    if (k >= 3) 
	return 3;
	
    /* Use Gram–Schmidt process to get orthogonal bases  */

    /* one constrains */
    if (k == 1) {
	/* fisrt basis */
	memcpy(H[0], _mb->normal + n[0]*Dim, Dim * sizeof(FLOAT));
	b[0] = _mb->b[n[0]];

	/* second basis */
	for (m = 0; m < Dim; m++)
	    if (fabs(c[0] = INNER_PRODUCT(H[0], bxyz[m])) < 0.9)
		break;
	assert(m < Dim);
	
	H[1][0] = bxyz[m][0] - c[0] * H[0][0]; 
	H[1][1] = bxyz[m][1] - c[0] * H[0][1]; 
	H[1][2] = bxyz[m][2] - c[0] * H[0][2]; 
	norm = sqrt(INNER_PRODUCT(H[1], H[1]));
	assert(norm > 1e-10);
	H[1][0] /= norm;
	H[1][1] /= norm;
	H[1][2] /= norm;

	/* third basis */
	for (m++; m < Dim; m++)
	    if (fabs(c[0] = INNER_PRODUCT(H[0], bxyz[m])) < 0.9)
		break;
	assert(m < Dim);
	c[1] = INNER_PRODUCT(H[1], bxyz[m]);
	H[2][0] = bxyz[m][0] - c[0] * H[0][0] - c[1] * H[1][0]; 
	H[2][1] = bxyz[m][1] - c[0] * H[0][1] - c[1] * H[1][1];  
	H[2][2] = bxyz[m][2] - c[0] * H[0][2] - c[1] * H[1][2];  
	norm = sqrt(INNER_PRODUCT(H[2], H[2]));
	assert(norm > 1e-10);
	H[2][0] /= norm;
	H[2][1] /= norm;
	H[2][2] /= norm;

	*bdry = &n[0];
	*bases = &H[0][0];
	*bb = &b[0];
	return 1;
    } 
    
    if (k == 2) {
	/* fisrt basis */
	memcpy(H[0], _mb->normal + n[0]*Dim, Dim * sizeof(FLOAT));
	b[0] = _mb->b[n[0]];

	/* second basis */
	memcpy(v, _mb->normal + n[1]*Dim, Dim * sizeof(FLOAT));
	c[0] = INNER_PRODUCT(H[0], v);
	assert(fabs(c[0]) < 0.9);
		     
	H[1][0] = v[0] - c[0] * H[0][0]; 
	H[1][1] = v[1] - c[0] * H[0][1]; 
	H[1][2] = v[2] - c[0] * H[0][2]; 
	b[1] = _mb->b[n[1]] - c[0] * b[0];
	norm = sqrt(INNER_PRODUCT(H[1], H[1]));
	assert(norm > 1e-10);
	H[1][0] /= norm;
	H[1][1] /= norm;
	H[1][2] /= norm;
	b[1] /= norm;

	/* third basis */
	for (m = 0; m < Dim; m++)
	    if (fabs(c[0] = INNER_PRODUCT(H[0], bxyz[m])) < 0.9
		&& fabs(c[1] = INNER_PRODUCT(H[1], bxyz[m])) < 0.9)
		break;
	assert(m < Dim);
	H[2][0] = bxyz[m][0] - c[0] * H[0][0] - c[1] * H[1][0]; 
	H[2][1] = bxyz[m][1] - c[0] * H[0][1] - c[1] * H[1][1];  
	H[2][2] = bxyz[m][2] - c[0] * H[0][2] - c[1] * H[1][2];  
	norm = sqrt(INNER_PRODUCT(H[2], H[2]));
	assert(norm > 1e-10);
	H[2][0] /= norm;
	H[2][1] /= norm;
	H[2][2] /= norm;

	*bdry = &n[0];
	*bases = &H[0][0];
	*bb = &b[0];
	return 2;
    }

    /* never get here */
    phgError(1, "Wrong num of bdry constrains!!!\n");
    return 0;
}

static void 
trans_left_(FLOAT *A, int ncol, int lda, const FLOAT *Trans) 
{					
    int i, j, k;							
    FLOAT tmp[Dim][ncol];						
    bzero(tmp, sizeof(tmp));					
    for (i = 0; i < Dim; i++)					
	for (j = 0; j < ncol; j++)					
	    for (k = 0; k < Dim; k++)				
		tmp[i][j] += Trans[i*Dim + k] * *(A + lda*k + j);	
    for (i = 0; i < Dim; i++)					
	for (j = 0; j < ncol; j++)					
	    *(A + lda*i + j) = tmp[i][j];				
}

static void 
trans_T_left_(FLOAT *A, int ncol, int lda, const FLOAT *Trans)  
{
    int i, j, k;							
    FLOAT tmp[Dim][ncol];						
    bzero(tmp, sizeof(tmp));					
    for (i = 0; i < Dim; i++)					
	for (j = 0; j < ncol; j++)					
	    for (k = 0; k < Dim; k++)				
		tmp[i][j] += Trans[i + k*Dim] * *(A + lda*k + j);	
    for (i = 0; i < Dim; i++)					
	for (j = 0; j < ncol; j++)					
	    *(A + lda*i + j) = tmp[i][j];				
}

static void 
trans_right_(FLOAT *A, int nrow, int lda, const FLOAT *Trans) 
{					
    int i, j, k;							
    FLOAT tmp[nrow][Dim];						
    bzero(tmp, sizeof(tmp));					
    for (i = 0; i < nrow; i++)					
	for (j = 0; j < Dim; j++)					
	    for (k = 0; k < Dim; k++)				
		tmp[i][j] += *(A + lda*i + k) * Trans[k + j*Dim];	
    for (i = 0; i < nrow; i++)					
	for (j = 0; j < Dim; j++)					
	    *(A + lda*i + j) = tmp[i][j];				
}

/*
 * Set move direction constrain:
 * 1. For fixed bdry, set all bdry vertex to be constrained,
 *    and do no bases rotation.
 * 2. For constrain bdry, rotate the bases and set the direction
 *    perpendicular to bdry to be constarined.
 * 
 *  */
static void 
setMoveConstrain(MovingMesh *mmesh, DOF *logical_node_new)
{
    GRID *g = _m->g;
    INT i, n, k, n_vert_constr, size_constr;
    const BYTE *bdry;
    const FLOAT *trans, *bb;
    FLOAT *v;

    n_vert_constr = 0;
    _m->vert_constr_index = phgAlloc(g->nvert * sizeof(*_m->vert_constr_index));
    _m->vert_constr = phgAlloc((size_constr = 3) * sizeof(*_m->vert_constr));
    for (i = 0; i < g->nvert; i++) {
	_m->vert_constr_index[i] = -1;
	if (g->types_vert[i] & BDRY_MASK) {
	    //printf("bdry vert:%4d\n", i);
	    n = getBdryVertBases(mmesh, i, &bdry, &trans, &bb);
	    if (n >= 3) {
		g->types_vert[i] |= MM_CONSTR;
		continue;
	    } else if (n == 2) {
		g->types_vert[i] |= (MM_CONSTR0 | MM_CONSTR1);
	    } else if (n == 1) {
		g->types_vert[i] |= MM_CONSTR0;
	    } else {
		phgError(-1, "vert bdry type unknonw!\n");
	    }

	    if (n_vert_constr >= size_constr) {				
		_m->vert_constr = phgRealloc(_m->vert_constr,
					     (size_constr *= 2) * sizeof(*_m->vert_constr));	
	    }
	    _m->vert_constr_index[i] = n_vert_constr;
	    _m->vert_constr[n_vert_constr].index = i;
	    _m->vert_constr[n_vert_constr].n_constrain = n;
	    memcpy(_m->vert_constr[n_vert_constr].bdry,  bdry, Dim * sizeof(*bdry));
	    memcpy(_m->vert_constr[n_vert_constr].Trans, trans, Dim * Dim * sizeof(FLOAT));
	    memcpy(_m->vert_constr[n_vert_constr].bb   , bb,    Dim * sizeof(FLOAT));
	    //printf("  # vert_constr: %4d [vert:%4d]\n", n_vert_constr, i);
	    n_vert_constr++;

	    /* rotate dof values w.r.t. rotated bases */
	    v = DofVertexData(logical_node_new, i);
	    trans_left_(v, 1, 1, trans);

	    if (0) {
		/* VTK draw debug */
		FLOAT coord[Dim], end[Dim], len = 0.1;
		printf("   vert:%d, constrain:%d\n", i, n);
		memcpy(coord, g->verts + i, Dim*sizeof(FLOAT));
		for (k = 0; k < Dim; k++) {
		    end[0] = coord[0] + len * trans[k*Dim+0];
		    end[1] = coord[1] + len * trans[k*Dim+1];
		    end[2] = coord[2] + len * trans[k*Dim+2];
		    if (k < n)
			vtk.SetColor("red");
		    vtk.DrawLine(coord, end);
		    vtk.SetColor("white");
		}
		vtk.Pause();
	    }
	} /* end of bdry vert*/
    }	  /* end of vert */
 
    return;
}

static void 
freeMoveConstrain(MovingMesh *mmesh, DOF *logical_node_new)
{
    GRID *g = _m->g;
    INT i; 
    FLOAT *v;

    /* rotate dof value back w.r.t. xyz bases */
    for (i = 0; i < g->nvert; i++) 
	if (MM_ROTATED(g->types_vert[i])) {
	    VERT_CONSTRAIN *vert_constr = _m->vert_constr 
		+ _m->vert_constr_index[i];
	    FLOAT *trans = vert_constr->Trans;
	    
	    assert(_m->vert_constr_index[i] >= 0);
	    //printf("  # get vert_constr: %d\n", _m->vert_constr_index[i]);
	    assert(i == vert_constr->index);
	    v = DofVertexData(logical_node_new, i);
	    trans_T_left_(v, 1, 1, trans);
	}

    return;
}


static void
getMoveDirection(MovingMesh *mmesh)
{
    GRID *g = _m->g;
    SIMPLEX *e;
    DOF *logical_move = _m->logical_move;
    DOF *logical_node = _m->logical_node;
    DOF *move = _m->move;
    DOF *logical_node_new, *grad_logical_node;
    INT i, j, k, l0, l1;
    SOLVER *solver;
    DOF *mass_lumping;
    FLOAT *v_move, *v_mass;
    FLOAT max_vol = -1e10, min_vol = 1e10;

    /* Get new monitor */
    _m->get_monitor(mmesh); DOF_SCALE(_m->monitor);

    logical_node_new = phgDofCopy(logical_node, NULL, NULL, "new logical node coordinates");
    grad_logical_node = phgDofGradient(logical_node, NULL, NULL, "grad logical node coordinates");
    /* if (logical_node_new->DB_mask == UNDEFINED) */
    /* 	phgError(1, "g->types_vert and DB_mask is VAILD in moving mesh method!!!"); */
    
    /* Create move solver */
    phgOptionsPush();
#if 0
    phgOptionsSetOptions("-default_solver mm_ml "
			 "-mm_ml_solver gmres "
			 "-mm_ml_pc mm_ml "
			 "-mm_ml_sweep0 6 "
			 "-solver_maxit 100 "
			 "-solver_rtol 1e-12 ");
#else
    phgOptionsSetOptions("-default_solver hypre "
			 "-hypre_solver gmres "
			 "-hypre_pc boomeramg "
			 "-solver_maxit 100 "
			 "-solver_rtol 1e-12 ");
#endif
    phgOptionsSetOptions(_mp->move_opts);

    setMoveConstrain(mmesh, logical_node_new);
    solver = phgSolverCreate(SOLVER_DEFAULT, logical_node_new, NULL);
    solver->mat->mv_data = phgAlloc(sizeof(*solver->mat->mv_data));
    solver->mat->mv_data[0] = (void *) mmesh;
    phgOptionsPop();

    /* build mat */
    ForAllElements(g, e) {
	int order = 1, q;
	int N = DofGetNBas(logical_node, e);	/* number of bases in the element */
	FLOAT A[N][Dim][N][Dim], rhs[N][Dim];
	INT I[N][Dim];
	QUAD *quad;
	const FLOAT *w, *lambda;
	FLOAT vol, mon;
	
	vol = phgGeomGetVolume(g, e);
	max_vol = MAX(max_vol, vol); 	
	min_vol = MIN(min_vol, vol); 

	for (i = 0; i < N; i++)
	    for (k = 0; k < Dim; k++)
		I[i][k] = phgMapE2L(solver->rhs->map, 0, e, i * Dim + k);

	bzero(A, sizeof(A));
	bzero(rhs, sizeof(rhs));

	quad = phgQuadGetQuad3D(order);
	mon = *DofElementData(_m->monitor, e->index); 
	//printf("mon:%e\n", mon);

	lambda = quad->points;
	w = quad->weights;
	assert(quad->npoints == 1);
	for (q = 0; q < quad->npoints; q++) {
	    for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
		    const FLOAT *ggi = 
			phgQuadGetBasisGradient(e, _m->move, i, quad) + q*Dim; /* grad phi_i */
		    const FLOAT *ggj = 
			phgQuadGetBasisGradient(e, _m->move, j, quad) + q*Dim; /* grad phi_j */
		    FLOAT a = vol*(*w) * mon * INNER_PRODUCT(ggj, ggi);
			
		    for (k = 0; k < Dim; k++) {
			A[i][k][j][k] += a;
		    }
		}
	    }
	    w++; lambda += Dim + 1;
	} /* end quad point */

	for (i = 0; i < N; i++) {
	    INT vert0 = e->verts[i];
	    if (_mp->fix_bdry || MM_FIXED(g->types_vert[vert0])) {
		/* no move */
		bzero(A[i], sizeof(A[i]));

		for (k = 0; k < Dim; k++) {
		    A[i][k][i][k] = 1.;
		    rhs[i][k] = logical_node->data[vert0 * Dim + k];
		}
	    } else if (g->types_vert[vert0] & MM_CONSTR) {
		VERT_CONSTRAIN *vert_constr = _m->vert_constr + _m->vert_constr_index[vert0];
		FLOAT *trans = vert_constr->Trans;
		FLOAT *bb = vert_constr->bb;
		assert(_m->vert_constr_index[vert0] >= 0);
		assert(vert0 == vert_constr->index);

		trans_left_ (&A[i][0][0][0], Dim*N, Dim*N, trans);
		trans_right_(&A[0][0][i][0], Dim*N, Dim*N, trans);

		if (vert_constr->n_constrain == 2) {
		    bzero(A[i][0], sizeof(A[i][0]));
		    bzero(A[i][1], sizeof(A[i][1]));
		    A[i][0][i][0] = 1.;
		    A[i][1][i][1] = 1.;
		    rhs[i][0] = bb[0];
		    rhs[i][1] = bb[1];
		} else if (vert_constr->n_constrain == 1) {
		    bzero(A[i][0], sizeof(A[i][0]));
		    A[i][0][i][0] = 1.;
		    rhs[i][0] = bb[0];
		} else {
		    abort();
		}
	    }
	}

	for (i = 0; i < N; i++)
	    phgSolverAddMatrixEntries(solver, Dim, &I[i][0], N * Dim, I[0], &A[i][0][0][0]);
	phgSolverAddRHSEntries(solver, N * Dim, &I[0][0], &rhs[0][0]);
    }
    
    if (_mp->verb > 1) {
	phgPrintf(MESSAGE_HEADER2"get logical move direction: build mat ");
	elapsed_time(g, TRUE, 0.);
    }

    //phgPrintf("    vol:[%e, %e]\n", max_vol, min_vol);
    if (0) {
	VEC *x = phgMapCreateVec(solver->rhs->map, 1);
	phgMapDofToLocalData(solver->rhs->map, 1, &_m->logical_node, x->data);

	phgPrintf("### Dumping solver mat vec...\n");
	phgMatDumpMATLAB(solver->mat, "A", "A_.m");
	phgVecDumpMATLAB(solver->rhs, "b", "b_.m");
	phgVecDumpMATLAB(x, "x", "x_.m");
    }

    phgSolverSolve(solver, FALSE, logical_node_new, NULL); DOF_SCALE(logical_node_new);
    if (_mp->verb > 1) {
	phgPrintf(MESSAGE_HEADER2"get logical move direction: solve: nits = %d", solver->nits);
	elapsed_time(g, TRUE, 0.);
    }
    phgSolverDestroy(&solver);
    freeMoveConstrain(mmesh, logical_node_new);

    /* get logical move direction */
    for (i = 0; i < g->nvert; i++) {
	if (!(g->types_vert[i] == UNREFERENCED))
	    for (k = 0; k < Dim; k++)
		logical_move->data[i*Dim+k] = 
		    logical_node->data[i*Dim+k] 
		    - logical_node_new->data[i*Dim+k];
    } /*  */
    DOF_SCALE(logical_move);

    if (_mp->viz_move) {
	GRID *g_new = phgDupGrid(g, TRUE);
	FLOAT *v = logical_node_new->data;

	for (i = 0; i < g_new->nvert; i++) {
	    if (!(g_new->types_vert[i] == UNREFERENCED)) {
		for (k = 0; k < Dim; k++)
		    g_new->verts[i][k] = v[i*Dim + k];
	    }
	} /* Point loop: local */
	//phgGeomInit_(g_new, TRUE); //DOF_SCALE(g_new->geom);

	sprintf(viz_name, "Moving_mesh.logical_new_%03d.vtk", viz_step);
	if (_mp->verb > 1)
	    phgPrintf(MESSAGE_HEADER2"Output moving details to %s\n", viz_name);
	phgExportVTK(g_new, viz_name, NULL);

	sprintf(viz_name, "Moving_mesh.logical_new_%03d.plt", viz_step);
	if (_mp->verb > 1)
	    phgPrintf(MESSAGE_HEADER2"Output moving details to %s\n", viz_name);
	//phgExportTecplot(g_new, viz_name, NULL);

	/* if (_mp->verb > 1) */
	/*     phgPrintf(MESSAGE_HEADER2"Output moving details to %s\n",  */
	/* 	      "moving_mesh.logical_new.ensight"); */
	/* phgExportEnsight(g_new, "moving_mesh.logical_new", (FLOAT) viz_step, NULL); */

	phgFreeGrid(&g_new);
    }

    /* get phyical move direction from logical move direction */
    mass_lumping = phgDofNew(g, DOF_P1, 1, "mass lumping on veticies", DofNoAction);
    phgDofSetDataByValue(mass_lumping, 0.);
    phgDofSetDataByValue(move, 0.);

    DofP1ReduceBegin(move);
    ForAllElements(g, e) {
	FLOAT tmp[NVert][Dim];
	const int v0 = e->verts[0];
	const int v1 = e->verts[1];
	const int v2 = e->verts[2];
	const int v3 = e->verts[3];
	COORD *_x0 = g->verts + v0; FLOAT *xi0 = logical_node->data + v0 * Dim;
	COORD *_x1 = g->verts + v1; FLOAT *xi1 = logical_node->data + v1 * Dim;
	COORD *_x2 = g->verts + v2; FLOAT *xi2 = logical_node->data + v2 * Dim;
	COORD *_x3 = g->verts + v3; FLOAT *xi3 = logical_node->data + v3 * Dim;
	FLOAT jacobi[3][3];
#define x0 (*_x0)
#define x1 (*_x1)
#define x2 (*_x2)
#define x3 (*_x3)
	FLOAT volume = ((xi1[0] - xi0[0])*(xi2[1] - xi0[1])*(xi3[2] - xi0[2]) +
			(xi1[1] - xi0[1])*(xi2[2] - xi0[2])*(xi3[0] - xi0[0]) +
			(xi1[2] - xi0[2])*(xi2[0] - xi0[0])*(xi3[1] - xi0[1]) -
			(xi1[0] - xi0[0])*(xi2[2] - xi0[2])*(xi3[1] - xi0[1]) -
			(xi1[1] - xi0[1])*(xi2[0] - xi0[0])*(xi3[2] - xi0[2]) -
			(xi1[2] - xi0[2])*(xi2[1] - xi0[1])*(xi3[0] - xi0[0]));

	jacobi[0][0] = (( x1[0] -  x0[0])*(xi2[1] - xi0[1])*(xi3[2] - xi0[2]) +
			(xi1[1] - xi0[1])*(xi2[2] - xi0[2])*( x3[0] -  x0[0]) +
			(xi1[2] - xi0[2])*( x2[0] -  x0[0])*(xi3[1] - xi0[1]) -
			( x1[0] -  x0[0])*(xi2[2] - xi0[2])*(xi3[1] - xi0[1]) -
			(xi1[1] - xi0[1])*( x2[0] -  x0[0])*(xi3[2] - xi0[2]) -
			(xi1[2] - xi0[2])*(xi2[1] - xi0[1])*( x3[0] -  x0[0]));
	jacobi[1][0] = (( x1[1] -  x0[1])*(xi2[1] - xi0[1])*(xi3[2] - xi0[2]) +
			(xi1[1] - xi0[1])*(xi2[2] - xi0[2])*( x3[1] -  x0[1]) +
			(xi1[2] - xi0[2])*( x2[1] -  x0[1])*(xi3[1] - xi0[1]) -
			( x1[1] -  x0[1])*(xi2[2] - xi0[2])*(xi3[1] - xi0[1]) -
			(xi1[1] - xi0[1])*( x2[1] -  x0[1])*(xi3[2] - xi0[2]) -
			(xi1[2] - xi0[2])*(xi2[1] - xi0[1])*( x3[1] -  x0[1]));
	jacobi[2][0] = (( x1[2] -  x0[2])*(xi2[1] - xi0[1])*(xi3[2] - xi0[2]) +
			(xi1[1] - xi0[1])*(xi2[2] - xi0[2])*( x3[2] -  x0[2]) +
			(xi1[2] - xi0[2])*( x2[2] -  x0[2])*(xi3[1] - xi0[1]) -
			( x1[2] -  x0[2])*(xi2[2] - xi0[2])*(xi3[1] - xi0[1]) -
			(xi1[1] - xi0[1])*( x2[2] -  x0[2])*(xi3[2] - xi0[2]) -
			(xi1[2] - xi0[2])*(xi2[1] - xi0[1])*( x3[2] -  x0[2]));
	jacobi[0][1] = ((xi1[0] - xi0[0])*( x2[0] -  x0[0])*(xi3[2] - xi0[2]) +
			( x1[0] -  x0[0])*(xi2[2] - xi0[2])*(xi3[0] - xi0[0]) +
			(xi1[2] - xi0[2])*(xi2[0] - xi0[0])*( x3[0] -  x0[0]) -
			(xi1[0] - xi0[0])*(xi2[2] - xi0[2])*( x3[0] -  x0[0]) -
			( x1[0] -  x0[0])*(xi2[0] - xi0[0])*(xi3[2] - xi0[2]) -
			(xi1[2] - xi0[2])*( x2[0] -  x0[0])*(xi3[0] - xi0[0]));
	jacobi[1][1] = ((xi1[0] - xi0[0])*( x2[1] -  x0[1])*(xi3[2] - xi0[2]) +
			( x1[1] -  x0[1])*(xi2[2] - xi0[2])*(xi3[0] - xi0[0]) +
			(xi1[2] - xi0[2])*(xi2[0] - xi0[0])*( x3[1] -  x0[1]) -
			(xi1[0] - xi0[0])*(xi2[2] - xi0[2])*( x3[1] -  x0[1]) -
			( x1[1] -  x0[1])*(xi2[0] - xi0[0])*(xi3[2] - xi0[2]) -
			(xi1[2] - xi0[2])*( x2[1] -  x0[1])*(xi3[0] - xi0[0]));
	jacobi[2][1] = ((xi1[0] - xi0[0])*( x2[2] -  x0[2])*(xi3[2] - xi0[2]) +
			( x1[2] -  x0[2])*(xi2[2] - xi0[2])*(xi3[0] - xi0[0]) +
			(xi1[2] - xi0[2])*(xi2[0] - xi0[0])*( x3[2] -  x0[2]) -
			(xi1[0] - xi0[0])*(xi2[2] - xi0[2])*( x3[2] -  x0[2]) -
			( x1[2] -  x0[2])*(xi2[0] - xi0[0])*(xi3[2] - xi0[2]) -
			(xi1[2] - xi0[2])*( x2[2] -  x0[2])*(xi3[0] - xi0[0]));
	jacobi[0][2] = ((xi1[0] - xi0[0])*(xi2[1] - xi0[1])*( x3[0] -  x0[0]) +
			(xi1[1] - xi0[1])*( x2[0] -  x0[0])*(xi3[0] - xi0[0]) +
			( x1[0] -  x0[0])*(xi2[0] - xi0[0])*(xi3[1] - xi0[1]) -
			(xi1[0] - xi0[0])*( x2[0] -  x0[0])*(xi3[1] - xi0[1]) -
			(xi1[1] - xi0[1])*(xi2[0] - xi0[0])*( x3[0] -  x0[0]) -
			( x1[0] -  x0[0])*(xi2[1] - xi0[1])*(xi3[0] - xi0[0]));
	jacobi[1][2] = ((xi1[0] - xi0[0])*(xi2[1] - xi0[1])*( x3[1] -  x0[1]) +
			(xi1[1] - xi0[1])*( x2[1] -  x0[1])*(xi3[0] - xi0[0]) +
			( x1[1] -  x0[1])*(xi2[0] - xi0[0])*(xi3[1] - xi0[1]) -
			(xi1[0] - xi0[0])*( x2[1] -  x0[1])*(xi3[1] - xi0[1]) -
			(xi1[1] - xi0[1])*(xi2[0] - xi0[0])*( x3[1] -  x0[1]) -
			( x1[1] -  x0[1])*(xi2[1] - xi0[1])*(xi3[0] - xi0[0]));
	jacobi[2][2] = ((xi1[0] - xi0[0])*(xi2[1] - xi0[1])*( x3[2] -  x0[2]) +
			(xi1[1] - xi0[1])*( x2[2] -  x0[2])*(xi3[0] - xi0[0]) +
			( x1[2] -  x0[2])*(xi2[0] - xi0[0])*(xi3[1] - xi0[1]) -
			(xi1[0] - xi0[0])*( x2[2] -  x0[2])*(xi3[1] - xi0[1]) -
			(xi1[1] - xi0[1])*(xi2[0] - xi0[0])*( x3[2] -  x0[2]) -
			( x1[2] -  x0[2])*(xi2[1] - xi0[1])*(xi3[0] - xi0[0]));

	//SHOW_M(&jacobi[0][0], 3, 3);
#undef x0
#undef x1
#undef x2
#undef x3

	bzero(tmp, sizeof(tmp));
	for (j = 0; j < NVert; j++) {
	    k = e->verts[j];
	    for (l0 = 0; l0 < Dim; l0++) {
		for (l1 = 0; l1 < Dim; l1++) {
		    //move->data[k * Dim + l0] += 
		    tmp[j][l0] +=
			jacobi[l0][l1] * logical_move->data[k *Dim + l1] *
			volume / fabs(volume); 
		    /* printf("### move: %3d, add: %e\n", k * Dim + l0,  */
		    /* 	   sign_vol * jacobi[l0][l1] * logical_move->data[k *Dim + l1]); */
		}
	    }
	}
	for (j = 0; j < NVert; j++) 
	    for (k = 0; k < Dim; k++)
		DofP1ReduceAdd(move, j*Dim + k,
			       tmp[j][k]);
	DofP1ReduceAddLocal(move);
    } /* P1->P1, different vertices, need reduce */
    DofP1ReduceEnd(move);

    DofP1ReduceBegin(mass_lumping);
    ForAllElements(g, e) {
	for (j = 0; j < NVert; j++)
	    DofP1ReduceAdd(mass_lumping, j,
			   6.*phgGeomGetVolume(g, e));
	DofP1ReduceAddLocal(mass_lumping);
    } /* P0->P1, need reduce */
    DofP1ReduceEnd(mass_lumping);

    DOF_SCALE(move);
    DOF_SCALE(mass_lumping);

    /* weighted average move direction */
    v_move = move->data;
    v_mass = mass_lumping->data;
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED)
	    continue;
	for (k = 0; k < Dim; k++) 
	    v_move[i*Dim+k] /= v_mass[i];
    } /* P1->P1, same point, no reduce */
    DOF_SCALE(move);
    
    /* constrain move direction */
    v_move = move->data;
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] == UNREFERENCED)
	    continue;

	if (MM_FIXED(g->types_vert[i])) {
	    for (k = 0; k < Dim; k++) 
		v_move[i*Dim+k] = 0;
	} else if (MM_ROTATED(g->types_vert[i])) {
	    VERT_CONSTRAIN *vert_constr = _m->vert_constr 
		+ _m->vert_constr_index[i];
	    BYTE *bdry = vert_constr->bdry;
	    
	    assert(_m->vert_constr_index[i] >= 0);
	    if (vert_constr->n_constrain == 1) {
		FLOAT *n = _mb->normal + bdry[0]*Dim, 
		    *v = v_move + i*Dim;
		FLOAT ip = INNER_PRODUCT(n, v);
		
		for (k = 0; k < Dim; k++)
		    v[k] -= ip * n[k];
	    } else if (vert_constr->n_constrain == 2) {
		FLOAT *n0, *n1, t[3], norm, ip, 
		    *v = v_move + i*Dim;

		n0 = _mb->normal + bdry[0]*Dim; 
		n1 = _mb->normal + bdry[1]*Dim; 

		t[0] = n0[1]*n1[2] - n0[2]*n1[1];
		t[1] = n0[2]*n1[0] - n0[0]*n1[2];
		t[2] = n0[0]*n1[1] - n0[1]*n1[0];
		norm = sqrt(INNER_PRODUCT(t, t));
		assert(norm > 1e-10);
		t[0] /= norm;
		t[1] /= norm;
		t[2] /= norm;
		
		ip = INNER_PRODUCT(v, t);
		for (k = 0; k < Dim; k++)
		    v[k] = ip * t[k];
	    } else {
		abort();
	    }
	}
    }

    if (_mp->verb > 1) {
	phgPrintf(MESSAGE_HEADER2"get physical move direction ");
	elapsed_time(g, TRUE, 0.);
    }
    DOF_SCALE(move);

    if (_mp->viz_move) {
	sprintf(viz_name, "Moving_mesh.move_%03d.vtk", viz_step);
	if (_mp->verb > 1)
	    phgPrintf(MESSAGE_HEADER2"Output moving details to %s\n", viz_name);
	phgExportVTK(g, viz_name,
		     logical_node, logical_node_new, logical_move, move, mass_lumping, NULL);
	

	sprintf(viz_name, "Moving_mesh.move_%03d.plt", viz_step);
	/* if (_mp->verb > 1) */
	/*     phgPrintf(MESSAGE_HEADER2"Output moving details to %s\n", viz_name); */
	/* phgExportTecplot(g, viz_name, */
	/* 		 logical_node, logical_node_new, logical_move, move, mass_lumping, NULL); */

	if (_mp->verb > 1)
	    phgPrintf(MESSAGE_HEADER2"Output moving details to %s\n", 
		      "Moving_mesh.move.ensight");
	phgExportEnsightT(g, "Moving_mesh.move", (FLOAT) viz_step, viz_step,
			 logical_node, logical_node_new, logical_move, move,
			 mass_lumping, _m->monitor, NULL);

	viz_step++;
	//exit(0);
    }

    {
	static int iii = 0;
	char name[1000];
	sprintf(name, "move_dir_%d", iii++);
	phgDofDumpFile(move, name);
    }
    phgDofFree(&mass_lumping);
    phgDofFree(&logical_node_new);
    phgDofFree(&grad_logical_node);
    return;
}

/* Acoording to page 381, "A moving mesh FEM algorithm for singular... ",
 * tau_star as the leaset positive root for the following equation should be evaluated,
 *     |     1           1           1           1      |
 * det |                                                | = 0.
 *     | x_0-t*dx_0  x_1-t*dx_1  x_2-t*dx_2  x_3-t*dx_3 |
 * 
 * */
static void 
getMoveStepLength(MovingMesh *mmesh)
{
    GRID *g = _m->g;
    SIMPLEX *e;
    DOF *move = _m->move;
    INT i, j, k;
    FLOAT a[4], move_step_length = 1.;

    Unused(i);

    ForAllElements(g, e) {
	FLOAT x[Dim][Dim], mx[Dim][Dim];
	for (j = 0; j < NVert - 1 ; j++) {
	    for (k = 0; k < Dim; k++) {
		x[j][k] = g->verts[e->verts[j+1]][k] - g->verts[e->verts[0]][k];
		mx[j][k] = move->data[e->verts[j+1] * Dim + k] - 
		    move->data[e->verts[0] * Dim + k];
	    }
	}

	//SHOW_M(&x[0][0], Dim, Dim);
	//SHOW_M(&mx[0][0], Dim, Dim);
#define DETERMINANT3(r0, r1, r2)					\
	(r0[0]*r1[1]*r2[2] + r0[1]*r1[2]*r2[0] + r0[2]*r1[0]*r2[1]	\
	 - r0[0]*r1[2]*r2[1] - r0[2]*r1[1]*r2[0] - r0[1]*r1[0]*r2[2])

	a[3] =  DETERMINANT3(mx[0], mx[1], mx[2]);
	a[2] = (DETERMINANT3( x[0], mx[1], mx[2]) +
		DETERMINANT3(mx[0],  x[1], mx[2]) +
		DETERMINANT3(mx[0], mx[1],  x[2]));
	a[1] = (DETERMINANT3(mx[0],  x[1],  x[2]) +
		DETERMINANT3( x[0], mx[1],  x[2]) +
		DETERMINANT3( x[0],  x[1], mx[2]));
	a[0] =  DETERMINANT3( x[0],  x[1],  x[2]);

#undef DETERMINANT3
	//SHOW_V(&a[0], NVert); 

	if (fabs(a[3]) < 1e-13) {
	    /* 2nd order equation */
	    
	    FLOAT a_ = a[2], b_ = a[1], c_ = a[0];
	    if (fabs(a_)/(fabs(b_) + fabs(c_)) < 1.0e-04) {
		if (fabs(b_) < 1.0e-04*fabs(c_))
		    ; //do nothing
		else if (c_/b_ > 0)
		    ; //do nothing
		else
		    move_step_length = MIN(move_step_length, -c_/b_);
	    }
	    else if (b_*b_ - 4*a_*c_ < 0)
		; //do nothing
	    else {
		if (a_ < 0) {
		    a_ = -a_;
		    b_ = -b_;
		    c_ = -c_;
		}
		FLOAT d = (-b_ - sqrt(b_*b_ - 4*a_*c_))/(2*a_);
		if (d < 0) {
		    d = (-b_ + sqrt(b_*b_ - 4*a_*c_))/(2*a_);
		    if (d > 0)
			move_step_length = MIN(move_step_length, d);
		}
		else
		    move_step_length = MIN(move_step_length, d);
	    }
	} else {
	    /* 3rd order equation */

	    a[0] /= a[3]; a[1] /= a[3]; a[2] /= a[3]; a[3] = 1.0;
	    FLOAT Q = (3*a[1] - a[2]*a[2])/9; //assert(fabs(Q) > 1e-30);
	    FLOAT R = (9*a[1]*a[2] - 27*a[0] - 2*a[2]*a[2]*a[2])/54;
	    FLOAT D = R*R + Q*Q*Q;
	    if (D > 0) { /* Only one real root */
		FLOAT S = cbrt(R + sqrt(D)), T = cbrt(R - sqrt(D));
		FLOAT root = -a[2]/3 + (S + T);
		if (root > 0 && move_step_length > root) {
		    move_step_length = root;
		}
	    }
	    else { /* Three real roots */
		FLOAT theta = fabs(Q) > 1e-30 ? acos(R/sqrt(-Q*Q*Q)) : M_PI/2;
		FLOAT root[3] = {2*sqrt(-Q)*cos(theta/3) - a[2]/3,
				 2*sqrt(-Q)*cos((theta + 2*M_PI)/3) - a[2]/3,
				 2*sqrt(-Q)*cos((theta - 2*M_PI)/3) - a[2]/3};
		for (j = 0; j < 3; j ++) {
		    if (root[j] > 0 && move_step_length > root[j]) {
			move_step_length = root[j];
		    }
		}
	    }
	}
    }

    move_step_length *= 0.5;
#if USE_MPI
    {
	FLOAT msl0 = move_step_length;
	MPI_Allreduce(&msl0, &move_step_length,
		      1, MPI_DOUBLE, MPI_MIN, g->comm);
    }
#endif /* USE_MPI */

    _m->move_step_length = move_step_length;
    if (_mp->verb > 1) {
	phgPrintf(MESSAGE_HEADER2"move step length = %e ", move_step_length);
	elapsed_time(g, TRUE, 0.);
    }

    return;
}

/* Smooth monitor: elem -> vert -> elem */
void 
phgMovingMeshSmoothMonitor(MovingMesh *mmesh, int s)
{
    GRID *g = _m->g;
    SIMPLEX *e;
    INT i, j, k;
    DOF *mass_lumping, *monitor1; 
    DOF *monitor = _m->monitor;

    mass_lumping = phgDofNew(g, DOF_P1, 1, "mass lumping on vertices", DofNoAction);
    monitor1 = phgDofNew(g, DOF_P1, 1, "monitor everage on vertices", DofNoAction);

    bzero(mass_lumping->data, g->nvert * sizeof(FLOAT));
    bzero(monitor1->data, g->nvert * sizeof(FLOAT));
    
    DofP1ReduceBegin(mass_lumping);
    ForAllElements(g, e) {
	for (j = 0; j < NVert; j++)
	    DofP1ReduceAdd(mass_lumping, j,
			   6.*phgGeomGetVolume(g, e));
	DofP1ReduceAddLocal(mass_lumping);
    } /* P0->P1, need reduce */
    DofP1ReduceEnd(mass_lumping);
    DOF_SCALE(mass_lumping);
	
    /* smooth s times */
    for (i = 0; i < s; i++) {
	bzero(monitor1->data, g->nvert * sizeof(FLOAT));
	DofP1ReduceBegin(monitor1);
	ForAllElements(g, e) {
	    for (k = 0; k < NVert; k++) 
		DofP1ReduceAdd(monitor1, k,
			       monitor->data[e->index] *
			       6.*phgGeomGetVolume(g, e));
	    DofP1ReduceAddLocal(monitor1);
	} /* P0->P1, need reduce */
	DofP1ReduceEnd(monitor1);
	//DOF_SCALE(monitor1);

	for (j = 0; j < g->nvert; j++) { 
	    if (!(g->types_vert[j] == UNREFERENCED))
		monitor1->data[j] /= 4*mass_lumping->data[j];
	} /* P1->P1, same point, no reduce */
	//DOF_SCALE(monitor1);
    
	bzero(monitor->data, g->nelem * sizeof(FLOAT));
	ForAllElements(g, e) {
	    for (k = 0; k < NVert; k++)
		monitor->data[e->index] += monitor1->data[e->verts[k]];
	} /* P1->P0, element wise, no reduce */
	DOF_SCALE(monitor);
    }
    
    phgDofFree(&mass_lumping);
    phgDofFree(&monitor1);
    return;
}

/* Update mesh vertices coordinates and other geom information.  */
static void 
updateMesh(MovingMesh *mmesh)
{
    GRID *g = _m->g;
    INT i, k;
    DOF *move = _m->move;
    FLOAT msl = _m->move_step_length;

    for (i = 0; i < g->nvert; i++) {
	if (!(g->types_vert[i] == UNREFERENCED)) {
	    for (k = 0; k < Dim; k++) 
		g->verts[i][k] +=
		    msl * move->data[i * Dim + k];
	}
    }

    phgGeomInit_(g, TRUE); //DOF_SCALE(g->geom);
}


/* MoveDirection() is not needed.
 * This function is mostly called when  we want update dof value
 * after mesh moving, and use the func on all quadrature points
 * in a element,
 * so we use a DOF_P1(Dim) to keep MoveDirection values,
 * and its value on a given point in the element can be easily
 * get from phgQuadGetDofValues().
 * */


void updateDof(MovingMesh *mmesh)
{
    GRID *g = _m->g;
    SIMPLEX *e;
    DOF **dofs = _m->dofs, *grad_dofs[10];
    //int ndof = _m->ndof;
    FLOAT msl = _m->move_step_length;
    MAP *map = _m->map;		/* This map has a DB_mask of UNDEFINED */
    MAT *mat;
    SOLVER *solver;
    QUAD *quad;
    const FLOAT *w, *lambda, *u, *gu, *mv;
    INT i, j, k, l, q;

    if (_m->ndof == 0)
	return;

    /* Create dof update solver */
    mat = phgMapCreateMat(map, map);
    phgOptionsPush();
    phgOptionsSetOptions("-default_solver hypre "
			 "-hypre_solver gmres "
			 "-hypre_pc boomeramg "
			 "-solver_maxit 100 "
			 "-solver_rtol 1e-12");
    phgOptionsSetOptions(_mp->dof_opts);
    solver = phgMat2Solver(SOLVER_DEFAULT, mat);
    phgOptionsPop();

    for (l = 1; l > 0; l--) {
	ForAllElements(g, e) {
	    int N = DofGetNBas(dofs[0], e);	/* number of basises in the element */
	    int order = 3; //DofTypeOrder(dofs[0], e) * 2;
	    FLOAT A[N][N], vol;
	    INT I[N];

	    assert(N == NVert);
	    bzero(A, sizeof(A));

	    vol = phgGeomGetVolume(dofs[0]->g, e);
	    quad = phgQuadGetQuad3D(order);
	    
	    lambda = quad->points;
	    w = quad->weights;
	    for (q = 0; q < quad->npoints; q++) {
		/* Mat */
		for (i = 0; i < N; i++) {
		    for (j = 0; j < N; j++) {
			const FLOAT *gi = phgQuadGetBasisValues(e, dofs[0], i, quad) + q; /* phi_i */
			const FLOAT *gj = phgQuadGetBasisValues(e, dofs[0], j, quad) + q; /* phi_j */
			A[i][j] += vol*(*w) * (*gi) * (*gj);
		    }
		}

		w++; lambda += Dim + 1;
	    }

	    for (i = 0; i < N; i++)
		I[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (i = 0; i < N; i++)
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
	} /* end build mass mat */
	if (_mp->verb > 1) {
	    phgPrintf(MESSAGE_HEADER2"updata dof: build mat ");
	    elapsed_time(g, TRUE, 0.);
	}

	/* Note: save DB_mask of dofs, and currently set it to UNDEFINED
	 *       to remove the effect of DIRICHLET bdry nodes. */
	for (k = 0; k < _m->ndof; k++) {
	    BTYPE DB_mask0[100];
	    memcpy(DB_mask0, dofs[k]->DB_mask, dofs[k]->dim * sizeof(*DB_mask0));

	    assert(dofs[k]->type == DOF_P1);
	    phgSolverResetRHS(solver);
	    grad_dofs[k] = phgDofGradient(dofs[k], NULL, NULL, NULL);	    
	    
	    /* build rhs */
	    ForAllElements(g, e) {
		int N = DofGetNBas(dofs[0], e);	/* number of basises in the element */
		int order = 3; //DofTypeOrder(dofs[0], e) * 2;
		FLOAT rhs[N], vol;
		INT I[N];

		assert(N == NVert);
		bzero(rhs, sizeof(rhs));

		vol = phgGeomGetVolume(dofs[0]->g, e);
		quad = phgQuadGetQuad3D(order);
	    
		u = phgQuadGetDofValues(e, dofs[k], quad);       /* u[k] */
		gu = phgQuadGetDofValues(e, grad_dofs[k], quad); /* grad u[k] */
		mv = phgQuadGetDofValues(e, _m->move, quad);     /* move dir*/

		lambda = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    for (i = 0; i < N; i++) {
			const FLOAT *gi = phgQuadGetBasisValues(e, dofs[k], i, quad) + q; /* phi_i */
			rhs[i] += vol*(*w) * (*gi) * (*u + msl * INNER_PRODUCT(gu, mv)
						      );
		    }

		    u++; gu += Dim; mv += Dim;
		    w++; lambda += Dim + 1;
		}

		for (i = 0; i < N; i++)
		    I[i] = phgSolverMapE2L(solver, 0, e, i);
		phgSolverAddRHSEntries(solver, N, I, rhs);

	    } /* end build rhs */

	    phgSolverSolve(solver, FALSE, dofs[k], NULL);
	    phgDofFree(&grad_dofs[0]);
	    memcpy(dofs[k]->DB_mask, DB_mask0, dofs[k]->dim * sizeof(*DB_mask0));
	    DOF_SCALE(dofs[k]);

	    if (_mp->verb > 1) {
		phgPrintf(MESSAGE_HEADER2"update dofs:%10s, nits = %d ", dofs[k]->name, solver->nits);
		elapsed_time(g, TRUE, 0.);
	    }
	} /* end solve updated dofs  */

	phgSolverDestroy(&solver);
	phgMatDestroy(&mat);
    }

    return;
}


void
phgMovingMeshDestroy(MovingMesh **_m_ptr)
{
    MovingMesh *_m = *_m_ptr;

    if (_m == NULL)
	return;

    _m_ptr = NULL;


    phgDofFree(&_m->logical_node);
    phgDofFree(&_m->logical_move);
    phgDofFree(&_m->move);
    phgDofFree(&_m->monitor);
    phgDofFree(&_m->phi);
    phgMapDestroy(&_m->map);
    //phgVerbosity = _m->verb;

    phgFree(_m);
    return;
}

#define FORMAT "%24.14E"
void
phgDofDumpFile(DOF *dof, char *file_name)
{
    char s[1024], *fmt;
    int i, j, k, index = 0;
    GRID *g;
    FLOAT *data, *points = NULL, x, y, z;
    SIMPLEX *e = NULL;
    COORD *p0 = NULL, *p1 = NULL, *p2 = NULL, *p3 = NULL;
    DOF_TYPE *type, *t;
    int *M = NULL, M_size = 0;
    VEF_MAP *vef = NULL;
    FILE *fp;
    char file_rank[1000];


    sprintf(file_rank, "%s.%d.dat", file_name, phgRank);
    fp = fopen(file_rank, "w");

    if (dof == NULL)
	return;

    if (dof->type == DOF_ANALYTIC) {
	if (dof->userfunc_lambda != NULL)
	    fprintf(fp, "Analytic DOF \"%s\", dim = %d, lambda function: %p\n",
		    dof->name, dof->dim, dof->userfunc_lambda);
	else
	    fprintf(fp, "Analytic DOF \"%s\", dim = %d, xyz function: %p\n",
		    dof->name, dof->dim, dof->userfunc);
	return;
    }
    else if (dof->type == DOF_CONSTANT) {
	if (dof->data == NULL) {
	    fprintf(fp, "Constant DOF \"%s\", dim = %d, no data\n",
		    dof->name, dof->dim);
	    return;
	}
	fprintf(fp, "Constant DOF \"%s\", dim = %d, values:\n",
		dof->name, dof->dim);
	for (i = 0, j = 0; i < dof->dim / 4; i++, j += 4) {
	    fprintf(fp, "    "FORMAT" "FORMAT" "FORMAT" "FORMAT"\n",
		    (double)dof->data[j], (double)dof->data[j + 1],
		    (double)dof->data[j + 2], (double)dof->data[j + 3]);
	}
	switch (dof->dim % 4) {
	    case 1:
		fprintf(fp, "    "FORMAT"\n", (double)dof->data[j]);
		break;
	    case 2:
		fprintf(fp, "    "FORMAT" "FORMAT"\n", (double)dof->data[j],
			(double)dof->data[j + 1]);
		break;
	    case 3:
		fprintf(fp, "    "FORMAT" "FORMAT" "FORMAT"\n",
			(double)dof->data[j], (double)dof->data[j + 1],
			(double)dof->data[j + 2]);
		break;
	}
	return;
    }

    if (!DofIsHP(dof))
	fprintf(fp, "DOF \"%s\": dim=%dx%d, "
		    "np_vert=%d, np_edge=%d, np_face=%d, np_elem=%d\n",
		    dof->name, dof->dim, dof->type->dim,
		    dof->type->np_vert, dof->type->np_edge,
		    dof->type->np_face, dof->type->np_elem);
    else
	fprintf(fp, "DOF \"%s\": dim=%dx%d, hierachical.\n",
		    dof->name, dof->dim, DofTypeDim(dof));

    if ((g = dof->g) == NULL || (data = dof->data) == NULL) {
	fprintf(fp, "no data.\n");
	return;
    }

    type = dof->type;
    if (type == NULL)
	type = dof->hp->info->types[dof->hp->max_order];

    /*if (type->points != NULL)*/
	vef = phgDofSetupVEFMap(g, dof, EDGE_FLAG | FACE_FLAG);

    if ((k = type->np_vert) > 0) {
	fmt = (k <= 1) ? "%5d  " : ((k < 10) ? "%5d[%d]  " : "%5d[%02d]  ");
	fprintf(fp, "Vertices:\n");
	for (i = 0; i < g->nvert; i++) {
	    if (g->types_vert[i] == UNREFERENCED) {
		data += dof->count_vert;
		continue;
	    }
	    for (j = 0; j < type->np_vert; j++) {
		sprintf(s, fmt, GlobalVertex(g, i), j);
		if (type->points != NULL)
		    sprintf(s + strlen(s), "%s(%6.3lf,%6.3lf,%6.3lf) = ",
			    dof->name, (double)g->verts[i][0],
			    (double)g->verts[i][1], (double)g->verts[i][2]);
		for (k = 0; k < dof->dim; k++)
		    sprintf(s + strlen(s), ""FORMAT" ", (double)*(data++));
		fprintf(fp, "%s\n", s);
	    }
	}
    }

    if ((k = type->np_edge) > 0) {
	if (M_size < k) {
	    phgFree(M);
	    M = phgAlloc((M_size = k) * sizeof(*M));
	}
	fmt = (k <= 1) ? "%5d  " : ((k < 10) ? "%5d[%d]  " : "%5d[%02d]  ");
	fprintf(fp, "Edges:\n");
	for (i = 0; i < g->nedge; i++) {
	    if (!DofIsHP(dof))
		t = type;
	    else
	 	t = dof->hp->info->types[dof->hp->edge_order[i]];
	    if (g->types_edge[i] == UNREFERENCED) {
		if (!DofIsHP(dof))
		    data += dof->count_edge;
		else
		    data += dof->dim *
			(dof->hp->edge_index[i + 1] - dof->hp->edge_index[i]);
		continue;
	    }
	    if (t->points != NULL) {
		assert(vef != NULL && vef->Emap[i] != NULL);
		e = vef->Emap[i];
		index = vef->Eind[i];
		p0 = g->verts + e->verts[GetEdgeVertex(index, 0)];
		p1 = g->verts + e->verts[GetEdgeVertex(index, 1)];
		points = t->points + t->np_vert;
	    }
	    phgDofMapEdgeData(t, e, index, M);
	    for (j = 0; j < t->np_edge; j++) {
		sprintf(s, fmt, GlobalEdge(g, i), j);
		if (t->points != NULL) {
		    x = (*p0)[0] * points[0] + (*p1)[0] * points[1];
		    y = (*p0)[1] * points[0] + (*p1)[1] * points[1];
		    z = (*p0)[2] * points[0] + (*p1)[2] * points[1];
		    sprintf(s + strlen(s), "%s(%6.3lf,%6.3lf,%6.3lf) = ",
			    dof->name, (double)x, (double)y, (double)z);
		    points += 2;
		}
		for (k = 0; k < dof->dim; k++)
		    sprintf(s + strlen(s), ""FORMAT" ",
			    (double)*(data + dof->dim * M[j] + k));
		fprintf(fp, "%s\n", s);
	    }
	    data += dof->dim * t->np_edge;
	}
    }

    if ((k = type->np_face) > 0) {
	if (M_size < k) {
	    phgFree(M);
	    M = phgAlloc((M_size = k) * sizeof(*M));
	}
	fmt = (k <= 1) ? "%5d  " : ((k < 10) ? "%5d[%d]  " : "%5d[%02d]  ");
	fprintf(fp, "Faces:\n");
	for (i = 0; i < g->nface; i++) {
	    if (!DofIsHP(dof))
		t = type;
	    else
	 	t = dof->hp->info->types[dof->hp->face_order[i]];
	    if (g->types_face[i] == UNREFERENCED) {
		if (!DofIsHP(dof))
		    data += dof->count_face;
		else
		    data += dof->dim *
			(dof->hp->face_index[i + 1] - dof->hp->face_index[i]);
		continue;
	    }
	    if (t->points != NULL) {
		assert(vef != NULL && vef->Fmap[i] != NULL);
		e = vef->Fmap[i];
		index = vef->Find[i];
		p0 = g->verts + e->verts[GetFaceVertex(index, 0)];
		p1 = g->verts + e->verts[GetFaceVertex(index, 1)];
		p2 = g->verts + e->verts[GetFaceVertex(index, 2)];
		points = t->points + t->np_vert + t->np_edge * 2;
	    }
	    phgDofMapFaceData(t, e, index, M);
	    for (j = 0; j < t->np_face; j++) {
		sprintf(s, fmt, GlobalFace(g, i), j);
		if (t->points != NULL) {
		    x = (*p0)[0] * points[0] + (*p1)[0] * points[1] +
			(*p2)[0] * points[2];
		    y = (*p0)[1] * points[0] + (*p1)[1] * points[1] +
			(*p2)[1] * points[2];
		    z = (*p0)[2] * points[0] + (*p1)[2] * points[1] +
			(*p2)[2] * points[2];
		    sprintf(s + strlen(s), "%s(%6.3lf,%6.3lf,%6.3lf) = ",
			    dof->name, (double)x, (double)y, (double)z);
		    points += 3;
		}
		for (k = 0; k < dof->dim; k++)
		    sprintf(s + strlen(s), ""FORMAT" ",
			    (double)*(data + dof->dim * M[j] + k));
		fprintf(fp, "%s\n", s);
	    }
	    data += dof->dim * t->np_face;
	}
    }

    if ((k = type->np_elem) > 0) {
	fmt = (k <= 1) ? "%5d  " : ((k < 10) ? "%5d[%d]  " : "%5d[%02d]  ");
	fprintf(fp, "Elements:\n");
	for (i = 0; i < g->nelem; i++) {
	    if (!DofIsHP(dof))
		t = type;
	    else
	 	t = dof->hp->info->types[dof->hp->elem_order[i]];
	    if (g->types_elem[i] == UNREFERENCED) {
		if (!DofIsHP(dof))
		    data += dof->count_elem;
		else
		    data += dof->dim *
			(dof->hp->elem_index[i + 1] - dof->hp->elem_index[i]);
		continue;
	    }
	    if (t->points != NULL) {
		e = g->elems[i];
		assert(e != NULL && e->index == i);
		p0 = g->verts + e->verts[0];
		p1 = g->verts + e->verts[1];
		p2 = g->verts + e->verts[2];
		p3 = g->verts + e->verts[3];
		points = t->points + t->np_vert + t->np_edge*2 + t->np_face*3;
	    }
	    for (j = 0; j < t->np_elem; j++) {
		sprintf(s, fmt, GlobalElement(g, i), j);
		if (t->points != NULL) {
		    x = (*p0)[0] * points[0] + (*p1)[0] * points[1] +
			(*p2)[0] * points[2] + (*p3)[0] * points[3];
		    y = (*p0)[1] * points[0] + (*p1)[1] * points[1] +
			(*p2)[1] * points[2] + (*p3)[1] * points[3];
		    z = (*p0)[2] * points[0] + (*p1)[2] * points[1] +
			(*p2)[2] * points[2] + (*p3)[2] * points[3];
		    sprintf(s + strlen(s), "%s(%6.3lf,%6.3lf,%6.3lf) = ",
			    dof->name, (double)x, (double)y, (double)z);
		    points += 4;
		}
		for (k = 0; k < dof->dim; k++)
		    sprintf(s + strlen(s), ""FORMAT" ", (double)*(data++));
		fprintf(fp, "%s\n", s);
	    }
	}
    }

    phgDofFreeVEFMap(&vef);
    phgFree(M);
    fclose(fp);
}
#undef FORMAT
