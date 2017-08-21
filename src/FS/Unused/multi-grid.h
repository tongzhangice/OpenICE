#ifndef MULTI_GRID_H

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <strings.h>	/* bzero() */
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
//#include "vtk-draw.h"

/* ---------------------------------------------------------------
 * Note:
 * 1. eliminate_bdry_colume:
 *    result in symetric mat, but can not pass LINE SMOOTHING
 *    COMPATABILITY check.
 * 
 *
 *  */

//#define MG_VTK_DEUBG
//#define ELIMINATE_BDRY_COLUME   
//#define DUMP_MAT_VEC TRUE
//#define TRANS_MAT_SET
//#define TEST_LINE_SMOOTHING 1
#define DIRECTION_CONNECTED 1
#define CHECK_LINE_NEIGH 3

#ifndef DUMP_MAT_VEC
#  define DUMP_MAT_VEC 0//FALSE
#endif


/* #ifdef TRANS_MAT_SET */
/* # define phgMatHandleEntries phgMatSetEntries */
/* #else */
/* # define phgMatHandleEntries phgMatAddEntries */
/* #endif	/\* TRANS_MAT_SET *\/ */

/* *************************************************************
 * MultiGrid Data structure:
 *   see alberta.h for suggestion.
 * *********************************************************** */
typedef GRID SUB_GRID;

/* Parameters */
typedef struct MG_PARAMS_ {
    int maxit;		        /* max iteration step */
    FLOAT rtol;			/* relative tolerance */
    FLOAT rtol_coarst;          /* relative tolerance on coarset level */

    DOF_TYPE *Ctype;		/* coarse Dof type */
    char *Ctype_name;		/* coarse DOF type name */

    int max_level;		/* current num of levels */
    int min_level;		/* coarsest level.
				 * mg create [crst_level, mg_levels]
				 * */
    int n_cycle;		/* 1=V-cycle, 2=W-cycle */
    BOOLEAN reuse_mat;		/* reuse finest mat */
    BOOLEAN solve_coarst;	/* use coarse grid solver */
    char *coarst_opts;		/* coarse grid solver options */

    int n_pre_smooth;		/* # of pre smoothing */
    int n_post_smooth;		/* # of post smoothing */
    int n_coarst_smooth;	/* # of smoothing on coarset level */

    FLOAT smooth_damp;		/* smooth damping parameter */
    FLOAT correct_damp;		/* correct damping parameter */
    FLOAT line_ar_tol;		/* line smoother AR tol */

    BOOLEAN use_GSline;		/* use GS line smoother */
    BOOLEAN use_upwind;		/* use upwind for convection */
    BOOLEAN use_low_order;	/* use lower order FEM on coarse grid */
    BOOLEAN export;		/* export level mesh */
    BOOLEAN timing;		/* MG timing */
} MG_PARAMS;

/* block dof smoother */
typedef struct MG_BLOCK_DOFS_ {
    int nds[10];

    MAP *map_cnnt;
    MAT *mat_cnnt;
    INT nline;
    INT *line_dofs;
    INT *line_dof0;
    INT *line_ndof;

    INT *lngh_dofs;
    INT *lngh_dof0;
    INT *lngh_ndof;
    int neigh_size;

    INT *line_piv0;
    INT *line_pivs;
    INT *line_mat0;
    FLOAT *line_matv;

    INT *line_poc0;
    INT *line_pocs;

    INT ndof_pt;
    INT *dof_pts;

    BOOLEAN factorized;
} MG_BLOCK_DOFS ;

typedef struct MULTI_GRID_ MULTI_GRID;
typedef struct MG_LEVEL_ MG_LEVEL;

/* auxiliary transfore operator */
typedef struct MG_TRANS_OP_ {
    DOF *dof_P;			/* Dof Pn */
    MAP *map_P;			/* Map Pn */
    VEC *vec_P;			/* Vec Pn */
    MAT *mat_CtoP;		/* MAT Pn to Coarse grid Dof */
    MAT *mat_PtoF; 		/* MAT Fine grid Dof to Pn */
} MG_TRANS_OP;

struct MG_LEVEL_ {
    MULTI_GRID    *mg;		/* Multigrid */

    GRID           *grid;	/* Trees */
    SUB_GRID       *sub_grid;	/* Grid copies */

    int           level;	/* Multigrid level */
    int           ndof;	/* # of dofs */
    DOF           **dofs;	/* Dofs */
    SOLVER        *solver;	/* Solver as smoother */
    MAP           *map;		/* Maps */
    MAT           *mat;		/* Stiff mat on level l*/
    MAT           *P;		/* Prolongation mat of level l to l-1
				 *   size (n_{l} x n{l-1})*/
    MAT           *R;		/* Restriction mat of level l to l-1
				 *   size (n_{l} x n{l-1})*/
    VEC           *x;		/* Vec solution */
    VEC           *f;		/* Vec defect */
    VEC           *r;		/* Vec defect */
    BTYPE         *types_vec;	/* Mask of proc boundary */
    BOOLEAN        redunt_level;  /* redundant level */

    /* --- Multi solver --- */
    int           nslvt;		/* # of solvers types*/
    DOF           *u;	        /* solver active dof */
    DOF           *solver_dof[10];  /* solver dofs */
    DOF_TYPE      *solver_type[10]; /* solver dof types */
    SOLVER        *solver_[10];	
    MAT           *mat_[10];	
    MAP           *map_[10];   	/* Maps */
    MAT           *P_[10];		
    MAT           *R_[10];		
    VEC           *x_[10];		
    VEC           *f_[10];		
    VEC           *r_[10];		
    BTYPE         *types_vec_[10];	
    BOOLEAN        redunt_level_[10]; /* redundant level */
    BOOLEAN        dummy;	/* dummy level */

    /* line soomther */
    MG_TRANS_OP   *transfer;	/* transfer intermediate Dof */
    MG_BLOCK_DOFS *block_dofs;	/* block dofs */
    MG_BLOCK_DOFS *block_dofs_[10];	/* block dofs of each solver*/
};

struct MULTI_GRID_ {
    GRID           *grid;	/* Current active grid */
    GRID           *grids[100];	/* Trees */
    int            ngrid;	/* # of trees */

    MG_LEVEL       **ml;	/* Coarse grid */

    int             level;	/* current grid level */
    int             start_level;/* start grid level, for new proc */
    //int max_level;		/* Max grid level */
    //int min_level;		/* Min grid level */
    int             active_solver_type;      /* Active solver type*/
    int             active_solver;           /* Active solver */
    MG_PARAMS       *mg_params;	/* Parameter list */
};

typedef void (*MG_SMOOTHER)(MAT *A, VEC *x, VEC *b, int nsmooth, void *ctx);
typedef BOOLEAN (*PRLG_RSTR_CALLBACK) CB_ARGS(e);
enum {DOWN_CYCLE, UP_CYCLE, COARSEST, FINEST};

extern MG_PARAMS *mg_params;

/* multigrid main server routine */
MG_PARAMS *phgMultiGridParametersCreate();
MULTI_GRID *phgMultiGridCreate(GRID *g, MG_PARAMS *mg_params);
void phgMultiGridInit(MULTI_GRID *mg);
void phgDestroyMultiGrid(MULTI_GRID **mg_ptr);
//void phgMultiGridBuildLevel2(MULTI_GRID *mg, int level, DOF *dof, ...);


void 
phgMultiGridBuildLevel_(MULTI_GRID *mg, int level, int ndof, DOF **dof_s, 
			int nslvt, int *solver2dof, BOOLEAN copy_dof_data, 
			BOOLEAN dummy_level);
#define phgMultiGridBuildLevel(mg, level, ndof, dofs, nslvt, solver2dof, copy_dof_data) \
    phgMultiGridBuildLevel_(mg, level, ndof, dofs, nslvt, solver2dof, copy_dof_data, TRUE)
#define phgMultiGridBuildDummyLevel(mg, level, ndof, nslvt)			\
    phgMultiGridBuildLevel_(mg, level, ndof, NULL, nslvt, NULL, FALSE, FALSE)

GRID *phgMultiGridRedistBegin(MULTI_GRID *mg, int level);
void phgMultiGridRedistEnd(MULTI_GRID *mg, int level);
void phgMultiGridDofRedist(MULTI_GRID *mg, int level, DOF *u0, DOF *u1, BOOLEAN forward);
void phgMultiGridUpdate(MULTI_GRID *mg);
void phgMultiGridBuildProlongationMat(MULTI_GRID *mg, int level);
void phgMultiGridBuildRestrictionMat(MULTI_GRID *mg, int level);
void phgMultiGridSmootherRegister(void);
MG_SMOOTHER phgMultiGridGetSmoother(MULTI_GRID *mg, int level, int type);
void phgMultiGridSolve(MULTI_GRID *mg, SOLVER *solver, DOF *u);
void phgMultiGridGetDofs(MG_LEVEL *ml, DOF **dof_ptr, ...);
void recursive_MG_cycle(MULTI_GRID *mg, int level, BOOLEAN init_zero);
void mg_pc_proc(SOLVER *pc, VEC *b0, VEC **x0);
GRID *phgDupGrid2(GRID *g, INT **map_ptr, BOOLEAN flatten);
void set_active_solver(MULTI_GRID *mg, int islvt, int islv);



/* smoother */
void mg_Jacobi(MAT *A, VEC *x, VEC *b, int nsmooth, void *pctx);
void mg_Jacobi2(MAT *A, VEC *x, VEC *b, int nsmooth, void *pctx);
void mg_GaussSidel(MAT *A, VEC *x, VEC *b, int nsmooth, void *pctx);
void mg_GaussSidel2(MAT *A, VEC *x, VEC *b, int nsmooth, void *pctx);
void mg_GaussSidel_vec(MAT *A, VEC *x, VEC *b, int nsmooth, void *pctx);
void mg_GaussSidel_line(MAT *A, VEC *x, VEC *b, int nsmooth, void *pctx);


#if 1
/* line smoother */
void phgMultiGridInitLineBlock(MG_LEVEL *ml, int islv, DOF *u);
void phgMultiGridFactorizeLineBlock(MG_LEVEL *ml);
void mg_GaussSidel_line(MAT *A, VEC *x, VEC *b, int nsmooth, void *pctx);
void phgMultiGridFreeFactLineBlock(MG_LEVEL *ml);
#else
#  define phgMultiGridInitLineBlock(ml, islv, u)
#  define phgMultiGridFactorizeLineBlock(ml)
#  define mg_GaussSidel_line NULL
#  define phgMultiGridFreeFactLineBlock(ml)
#endif

#define DOF_TYPE_LAGRANGE(type)			\
    (!strncmp(type->name, "P", 1)		\
     && type->name[1] >= '0'			\
     && type->name[1] <= '9')

#define DOF_TYPE_SUBELEM(type)			\
    (!strncmp(type->name, "SUB", 3)		\
     && type->name[3] >= '0'			\
     && type->name[3] <= '9')

#define MG_PREFIX "MultiGrid: "

#define PACK_COL(ps, i) (ps[i])
#define PACK_DAT(ps, i) (ps[i] - i)
#define PACK_COL_OFFP(ps, i, nlocal) (ps[i+nlocal])
#define PACK_DAT_OFFP(ps, i, nlocal) (ps[i+nlocal] - i - nlocal)

#ifndef LIE_IN
#  define LIE_IN(x) ((x) >= (- 1E-13) && (x) <= (1. + 1E-13)) 
#  define LambdaInElement(lambda) (LIE_IN(lambda[0]) && LIE_IN(lambda[1]) && \
				   LIE_IN(lambda[2]) && LIE_IN(lambda[3]))		
#endif

#define REALLOC_VEC(vecv, vec_size, vec_alloc, add_alloc) {	\
	INT vec_alloc##0 = vec_alloc;				\
	Unused(vec_alloc##0);					\
	while ((vec_size) + (add_alloc) > vec_alloc) {		\
	    vec_alloc *= 2;					\
	    vecv = phgRealloc_(vecv, vec_alloc * sizeof(*vecv),	\
			       vec_alloc##0 * sizeof(*vecv));	\
	}							\
    }

#define MAT_EPS 1E-14
#if 1
#define REMOVE_ZERO(v, N) {			\
	int _i_;				\
	for (_i_ = 0; _i_ < N; _i_ ++) {	\
	    if (fabs(v[_i_]) < MAT_EPS)		\
		v[_i_] = 0.;			\
	    if (fabs(1. - v[_i_]) < MAT_EPS)	\
		v[_i_] = 1.;			\
	}					\
    }
#else
#define REMOVE_ZERO(v, N) 
#endif


#define ZeroVec(x) phgVecAXPBY(0., NULL, 0., &(x));
//#define MG_DEBUG(...) printf(__VA_ARGS__);
#define MG_DEBUG(...) phgInfo(1, __VA_ARGS__);
#define MG_DEBUGn(n, ...) phgInfo(n, __VA_ARGS__);
#define DEBUG_TIME_COST(...) {			\
	phgPrintf("   ");			\
	phgPrintf(__VA_ARGS__);		\
	elapsed_time(g, TRUE, 0.);	\
    }

#define DEBUG_PAUSE			{				\
	if (debug) {							\
	    if (g->rank == 0) {						\
		int i = 0;						\
		MG_DEBUG("Debug pause control     in func:%-15s, line:%d\n", \
			__FUNCTION__, __LINE__);			\
		while (0 == i)						\
		    sleep(50000);					\
		MPI_Barrier(g->comm);					\
	    } else {							\
		MG_DEBUG("Debug pause waiting ... in func:%-15s, line:%d\n", \
			__FUNCTION__, __LINE__);			\
		MPI_Barrier(g->comm);					\
	    }								\
	}								\
    }

#define MEM_DUP(dest, src, count)				\
    MG_DEBUGn(3, "# mem copy: (%s) ===> (%s)\n", #src, #dest);	\
    if(src == NULL) {						\
	dest = NULL;						\
	MG_DEBUGn(3, "#   NULL copyed\n");			\
    } else {							\
	dest = phgCalloc((count), sizeof(src[0]));		\
	memcpy(dest, src, (count)*sizeof(src[0]));		\
	MG_DEBUGn(3, "#   %d copyed\n", (count));		\
    }
/* TODO: deal with analytic dof */
#define DOF_DATA_DUP(dest, src)						\
    MEM_DUP((dest)->data, (src)->data, DofGetDataCount((src)));		\
    (dest)->data_vert = (dest)->data;					\
    (dest)->data_edge = (dest)->data_vert + DofGetVertexDataCount(src);	\
	  (dest)->data_face = (dest)->data_edge + DofGetEdgeDataCount(src); \
	  (dest)->data_elem = (dest)->data_face + DofGetFaceDataCount(src);



#define INNER_PRODUCT(p, q)			\
    (*(p    ) * *(q    ) +			\
     *(p + 1) * *(q + 1) +			\
     *(p + 2) * *(q + 2))

#define MG_DISABLED(desp)					\
    phgError(1, "MultiGrid does NOT enable %s now\n", desp);


/* Timing */
#define TIMING_BEGIN if (_mgp->timing) {        \
        MPI_Barrier(g->comm);                   \
        time0 = phgGetTime(NULL);               \
    }
#define TIMING_END(desp) if (_mgp->timing) {            \
        MPI_Barrier(g->comm);                           \
        if (g->rank == 0)                               \
        phgInfo(0, "   Slv:%d, L:%d, T:%0.8lfs, %s\n",  \
                islv, level,                            \
                phgGetTime(NULL) - time0,               \
                desp);                                  \
    }

#define MULTI_GRID_H
#endif
