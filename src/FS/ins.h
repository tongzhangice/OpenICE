#ifndef INS_H
#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <strings.h>	/* bzero() */
#include <math.h>
#include <stdarg.h>

#if USE_MG 
#  include "multi-grid.h"
#endif /* USE_MG */
#include "layers.h"

/*
 * ================================================================================
 * 
 *      Control macro
 * 
 * ================================================================================
 * */

/*
 * Note:
 * 1. face quad for curved element is NOT implemented.
 *
 *  */

/* Control Macros */
#define STEADY_STATE 1		/* Steaty or time-dependent */
#define TIME_DEP_NON 0          /* Time-dependent case, using non-linear solver  */

/* Grid */
#define USE_ISOP 0		/* Iso-parametric elements */

/* Stokes pin node */
#define PIN_AT_ROOT 0			/* Pin node options */

/* Linear solver */
#define USE_QP_ONLY 1			/* Use Qp in PCD precondtioner */
#define REUSE_MAT 0                     /* Reuse matrix, if matrix is fixed for each time step */
#define MAT_HANDLE_BDRY_EQNS FALSE
#ifndef DUMP_MAT_VEC			/* Dump mat vec for debug */
#   define DUMP_MAT_VEC 0
#endif


/* Sliding B.C. */
#define USE_NODAL_LOADS 1
#define USE_SLIDING_BC 0
#define CASE_DRAINAGE_BASIN 0
#define ZERO_TRACT_ZONE 0
#define SLIP_BDRY BC_BOTTOM
//#define SLIP_BDRY BC_BOTTOM2


/* Temp solver */
#define USE_TEMP_SDPG 1
#define USE_TEMP_TIME 1
#define USE_TEMP_CONV 1
#define USE_TEMP_DIFF 1
#define USE_TEMP_HEAT 1
#define USE_TEMP_GEOFLUX 1


/* Mulitgrid solver */
#define USE_MG 0
#define MG_SOLVER_TYPE_U  0
#define MG_SOLVER_TYPE_P  1

#define MG_SOLVER_F  0
#define MG_SOLVER_Ap 1
#define MG_SOLVER_Qp 2


/*
 * ================================================================================
 * 
 *                 Test problem
 * 
 * ================================================================================
 * */
#define DRIVEN_CAVITY 100
#define CYLINDER 101
#define KOVASZNAY 109
#define ICE_BENCH_A 102
#define ICE_BENCH_B 103
#define ICE_BENCH_C 104
#define ICE_BENCH_D 105
#define ICE_BENCH_E 106
#define ICE_BENCH_F 107
#define ESIMINT_A 112
#define ESIMINT_B 113
#define ESIMINT_C 114
#define ESIMINT_D 115
#define ESIMINT_E 116
#define ESIMINT_F 117
#define ESIMINT_G 118
#define ESIMINT_H 119
#define HEINO_A 122
#define HEINO_B 123
#define HEINO_C 124
#define HEINO_D 125
#define HEINO_E 126
#define HEINO_F 127
#define ICE_EXACT      130
#define ICE_GREEN_LAND 200
#define TEST 201
#define LAS 202
#define LAT_PRES_BC 203
#define ICE_SHELF_BUTS_2D 204
#define AMERY_ICE_SHELF 205

#define TEST_CASE AMERY_ICE_SHELF
#define NS_PROBLEM "AMERY_ICE_SHELF"


/* Bench test */
#define ICE_BENCH_TEST 0
#define ESIMINT_TEST 0
#define HEINO_TEST 0

#if TEST_CASE == ICE_BENCH_A \
    || TEST_CASE == ICE_BENCH_B \
    || TEST_CASE == ICE_BENCH_C \
    || TEST_CASE == ICE_BENCH_D \
    || TEST_CASE == ICE_BENCH_E \
    || TEST_CASE == ICE_BENCH_F
#  undef ICE_BENCH_TEST
#  define ICE_BENCH_TEST 1
#elif TEST_CASE == ESIMINT_A \
    || TEST_CASE == ESIMINT_B \
    || TEST_CASE == ESIMINT_C \
    || TEST_CASE == ESIMINT_D \
    || TEST_CASE == ESIMINT_F \
    || TEST_CASE == ESIMINT_G \
    || TEST_CASE == ESIMINT_H
#  undef ESIMINT_TEST
#  define ESIMINT_TEST 1
#elif TEST_CASE == HEINO_A \
    || TEST_CASE == HEINO_B \
    || TEST_CASE == HEINO_C \
    || TEST_CASE == HEINO_D \
    || TEST_CASE == HEINO_E \
    || TEST_CASE == HEINO_F
#  undef HEINO_TEST
#  define HEINO_TEST 1
#endif




/* Scaling */
#if ICE_BENCH_TEST
#  define EQU_SCALING 1e-8
#  define LEN_SCALING 1e3
#  define PRES_SCALING 1e5
#  define EQU_T_SCALING 1e0
#elif ESIMINT_TEST
#  define EQU_SCALING 1e-8
#  define LEN_SCALING 1e3
#  define PRES_SCALING 1e5
#  define EQU_T_SCALING 1e-5
#elif HEINO_TEST
#  define EQU_SCALING 1e-8
#  define LEN_SCALING 1e3
#  define PRES_SCALING 1e5
#  define EQU_T_SCALING 1
#elif TEST_CASE == ICE_GREEN_LAND
#  define EQU_SCALING 1
#  define LEN_SCALING 1
#  define PRES_SCALING 1
#  define EQU_T_SCALING 1e12
#elif TEST_CASE == LAS
#  define EQU_SCALING 1e-5
#  define LEN_SCALING 1
#  define PRES_SCALING 1e5
#  define EQU_T_SCALING 1e12
#elif TEST_CASE == TEST
#  define EQU_SCALING 1e-5
#  define LEN_SCALING 1
#  define PRES_SCALING 1e5
#  define EQU_T_SCALING 1e12
#elif TEST_CASE == LAT_PRES_BC
#  define EQU_SCALING 1e-5
#  define LEN_SCALING 1
#  define PRES_SCALING 1e5
#  define EQU_T_SCALING 1e12
#elif TEST_CASE == ICE_SHELF_BUTS_2D
#  define EQU_SCALING 1e-5
#  define LEN_SCALING 1
#  define PRES_SCALING 1e5
#  define EQU_T_SCALING 1e12
#elif TEST_CASE == AMERY_ICE_SHELF
#  define EQU_SCALING 1e-5
#  define LEN_SCALING 1
#  define PRES_SCALING 1e5
#  define EQU_T_SCALING 1e12
#elif TEST_CASE == ICE_EXACT
#  define EQU_SCALING 1e-8
#  define LEN_SCALING 1e3
#  define PRES_SCALING 1e5
#  define EQU_T_SCALING 1e0
#else
#  define EQU_SCALING 1
#  define LEN_SCALING 1
#  define PRES_SCALING 1
#endif

#  define LEN_SCALING2 (LEN_SCALING * LEN_SCALING)
#  define LEN_SCALING3 (LEN_SCALING * LEN_SCALING * LEN_SCALING)

/* Const parameters */
#include "parameters.h"
#include "netcdf-utils.h"


/*
 * ================================================================================
 * 
 *          discretization options
 * 
 * ================================================================================
 * */


/* Bdry type */
#define BC_TOP       (BDRY_USER1)
#define BC_BOTTOM    (BDRY_USER3)
#define BC_BOTTOM_GRD    (BDRY_USER8)
#define BC_ISHELF    (BDRY_USER7)
#define BC_LATERL    (BDRY_USER5)
#define BC_LATERL_GRD    (BDRY_USER9)
#define SETFLOW      (BC_BOTTOM)
#define BC_DIVIDE    (BDRY_USER2)
#define BC_TERMNS    (BDRY_USER4)
#define BC_BOTTOM_GRD_ADD (BDRY_USER6)

#define INFLOW  BC_BOTTOM
#define OUTFLOW BC_TOP

//#define PINNEDNODE   (BDRY_USER6)

typedef enum { PICARD, NEWTON } LTYPE;	/* linearization type */

/*
 * ================================================================================
 * 
 *          Sturcture
 * 
 * ================================================================================
 * */

/* Parameters */
typedef struct NSParams_ {
    DOF_TYPE *utype;             /* DOF type for velocity */
    DOF_TYPE *ptype;		  /* DOF type for pressure */
    DOF_TYPE *T_type;           /* DOF type for Temperature */
    char *utype_name;		  /* DOF type name for velocity */
    char *ptype_name;		  /* DOF type name for pressure */
    char *T_type_name;	          /* DOF type name for coordinates */
    
    BOOLEAN layered_mesh;	  /* Use layered mesh */
    BOOLEAN struct_mesh;	  /* Use struct mesh */
    BOOLEAN solve_temp;	          /* Solve temperature */
    BOOLEAN solve_height;	  /* Solve height */
    BOOLEAN update_bdry_type;	  /* Update boundary btype by 2D mesh */
    BOOLEAN resume;	          /* Resume step from previous step */
    BOOLEAN record;	          /* Save step data to resume */
    BOOLEAN enclosed_flow;        /* Enclosed flow */
    BOOLEAN pin_node;		  /* Pin a node */
    BOOLEAN curved_bdry;	  /* Use curved boundary */
    BOOLEAN curved_refine;	  /* Use curved refinement */
    BOOLEAN start_const_vis;	  /* Viscosity start at constant */
    BOOLEAN compensate_equ;	  /* Compensate euqations */
    BOOLEAN add_ice_shelf;	  /* add ice shelf part */
    BOOLEAN lateral_sia_pres_bc;	  /* add lateral bc with SIA pressure */
    BOOLEAN another_run_with_updated_mask;	  /* if do another_run_with_updated_mask */
    int Nx, Ny, Nz;		  /* Structured mesh size */

    INT periodicity;		  /* periodicity */
    FLOAT Re;			  /* Reynolds num */
    FLOAT nu;			  /* Viscosity */
    FLOAT time_start;		  /* Time start */
    FLOAT time_end;		  /* Time end */
    FLOAT dt0;			  /* Inintial dt */
    FLOAT eps_diagP;		  /* Diagnal entries for pressure matrix */
    FLOAT fp_scale;		  /* Dirichlet scale for Fp */

    INT slip_condition;
    FLOAT slip_beta2;
    FLOAT slip_alpha2; 
    FLOAT slip_index;

    /* Discrete scheme */
    INT height_scheme;		  /* Height solver scheme */
    FLOAT Theta;		  /* 0.5: Crank-Nicolson
				   * 1. : Backward Euler */
    BOOLEAN use_PCD;		  /* Use Pressure Convetion-Diffusion preconditioner
				   *   default: True */
    int init_temp_type;		  /* Init temperature field
				   * 0: diff, 1: interp, 2: read data */

#if USE_MG 
    BOOLEAN use_mg_F;		  /* Use Multigrid solver */
    BOOLEAN use_mg_Ap;		  /* Use Multigrid solver */
    BOOLEAN use_mg_Qp;		  /* Use Multigrid solver */
#endif /* USE_MG */
    BOOLEAN use_moc;		  /* Use method of characterics */
    BOOLEAN use_Fu;		  /* Use Fu & Qu in the preconditioner, default: True */
    BOOLEAN use_Zp;		  /* Use Zp in the preconditioner, default: True */
    BOOLEAN implicit_convect;     /* Implicit scheme for convetion term,
				   *   default True */
    BOOLEAN use_symetric;	  /* Mat and PC is symetric, for symetric check */

    /* Faster code */
    BOOLEAN extern_force;	  /* Extern force, default ture! */

    /* Solver options */
    BOOLEAN non_linear;		  /* Nonlinear iteration */
    BOOLEAN noniter_temp;	  /* nonlinear iterative of velocity and temperature */
    FLOAT non_tol;		  /* Nonlinear iteration tolerance */
    FLOAT non_sub_tol;		  /* Nonlinear iter: sub linear problem tolerance */

    INT pre_refines;		  /* Grid pre-refine level */
    INT max_tstep;		  /* Max time step */
    INT max_nonstep;		  /* Max nonlinear interation step,
				   * -1 means P2P1 is skiped */
    INT min_nonstep;		  /* Min nonlinear interation step */
    INT max_nonstep0;		  /* Max nonlinear interation step for first step,
				   * negtive means using max_nonstep instead. */
    INT newton_start;		  /* Newton start step */
    INT newton_start0;		  /* Newton start step for first step,
				   * negtive means using max_nonstep instead. */
    FLOAT u_tol0;
    FLOAT p_tol0;
    FLOAT u_tol;
    FLOAT p_tol;
    FLOAT s_tol;

    INT step_span;		  /* Step span to output geo file */
    INT step_span_resume;		  /* Step span to output geo file */
    INT mem_max;		  /* Max memory per process */
    INT n_bdry_layer;		  /* # of boundary layers */
    /* INT moc_quad_order;		  /\* MOC quad order *\/ */
    /* INT moc_quad_nelem;		  /\* MOC quad nelem *\/ */
    BOOLEAN compute_error;	  /* Compute error */

    char *fn;			  /* Mesh file */
    char *resume_mesh;		  /* Resume mesh file */
    char *resume_data;		  /* Resume data file */
    char *Stokes_opts;		  /* Solver Stokes options */
    char *F_opts;		  /* Solver F options*/
    char *Fu_opts;		  /* Solver Fu options*/
    char *Ap_opts;		  /* Solver Ap options*/
    char *Qp_opts;		  /* Solver Qp options*/
    char *Fp_opts;		  /* Solver Fp options*/
    char *Gu_opts;		  /* Grad u options*/
    char *T_opts;		  /* Solver temperature opts */

    /* 2D file */
    char *tria_file;
    char *vert_file;
    char *layer_file;
    char *nodeZ_file;
    char *dual_file;
    char *vx_txt_file;
    char *vy_txt_file;
    char *bot_txt_file;
    char *mask_txt_file;
    char *sur_txt_file;
    char *thk_txt_file;
    char *x_txt_file;
    char *y_txt_file;
    char *sur_grad_x_txt_file;
    char *sur_grad_y_txt_file;

    /* Netcdf file */
    char *nc_file;
} NSParams;

/* PCD Preconditioner */
typedef struct NSPC_ {
    DOF *pbc;			  /* pressure bdry condition
				   * for preconditioner */
    MAP *Pbcmap;		  /* pressure dof map for preconditioner */
    DOF *u1;			  /* velocity component u */
    MAP *u1map;			  /* velocity component u map */

    MAT *matFu;			  /* matrix of convection-diffusion */
    MAT *matQp;			  /* matrix of pressure mass */
    MAT *matAp;			  /* matrix of pressure diffusion */
    MAT *matFp;			  /* matrix of pressure convection-diffusion */
    VEC *rhsScale;		  /* rhs scale */

    /* DOF and map to specify boundary condition for pressure
     * convection-diffusion problem. */
    DOF *dof_inflow, *dof_outflow, *dof_nobdry;
    MAP *map_inflow, *map_outflow, *map_nobdry;

    SOLVER *solver_F;		  /* PC solver of velocity convection-diffusion */
    SOLVER *solver_Fu;		  /* PC solver of velocity convection-diffusion, seperated */
    SOLVER *solver_Ap;		  /* PC solver of pressure diffusion */
    SOLVER *solver_Qp;		  /* PC solver of pressure mass */
    SOLVER *solver_Fp;		  /* PC solver of pressure convection-diffusion */
    SOLVER *pc_F;		  /* mg shell PC solver of velocity convection-diffusion */
    SOLVER *pc_Ap;		  /* mg shell PC solver of velocity convection-diffusion */
    SOLVER *pc_Qp;		  /* mg shell PC solver of velocity convection-diffusion */
    SOLVER *pc_Zp;		  /* mg shell PC solver of velocity convection-diffusion */
} NSPC;				  

/* Surface bases */
typedef struct SURF_BAS_ {
    DOF_TYPE *type;
    DOF *dof;
    BOOLEAN *rotated;
} SURF_BAS;

/*
typedef struct GEO_INFO_{
    DOF *ice_sur;
    DOF *ice_bed;
    DOF *ice_thk;
    DOF *sur_grad_x;
    DOF *sur_grad_y;
} GEO_INFO;
*/
/* Main solver */
typedef struct NSSolver_ {
    GRID *g;			  /* Grid */
    DOF **u;			  /* velocity ptr */
    DOF **p;			  /* pressure ptr */
    DOF **T;			  /* conformation tensor */
#if STEADY_STATE || TIME_DEP_NON 
    DOF *du;			  /* delta u in non-linear iteration */
    DOF *dp;			  /* delta p in non-linear iteration */
    DOF *dT;			  /* delta C in non-linear iteration */
#endif /* STEADY_STATE || TIME_DEP_NON */
    DOF **gradu;		  /* gradient velocity ptr */
    DOF **lapu;			  /* laplace velocity ptr */
    DOF **gradp;		  /* gradient pressure ptr */
    DOF *Gradu;		  /* gradient temperature ptr */

    DOF *nodal_force;
    DOF *nodal_force_value;
    DOF *water_force;
    DOF *contact_force;
    DOF *contact_force_value;

    DOF *stress;
    DOF *stress1;
    DOF *water_pressure;
    DOF *water_pressure1;
    DOF *stress_nn;
    DOF *stress_nn1;
    DOF *mask_bot;
    DOF *avg_gu;
    DOF *viscosity;
    DOF *strain_rate;
    DOF *eu_d;
    DOF *eu_n;
        
    DOF *f;			  /* source term momentum  */
    DOF *u_queue[3];		  /* velocity at different time */
    DOF *p_queue[3];		  /* pressure at different time */
    DOF *T_queue[3];		  
    DOF *gradu_queue[3];	  /* gradient velocity at different time */
    DOF *u_shape;		  /* velocity shape DOF */
    DOF *p_shape;		  /* pressure shape DOF */
    DOF *T_shape;		  /* proj Gradient vel shape DOF */
    DOF *gn[3];			  /* Outflow bdry condition for velocity & pressure*/
    DOF *wind;			  /* predicted velocity */
    DOF *dH;			  /* coord change */
    DOF *beta;			  /* slip coef */
    DOF *coord;			  /* coord(P1) */

    INT pinned_node_id;	          /* Vert index of pinned node at rank 0
				   * -1 for no pinned. */
    MAP *Vmap;			  /* velocity dof map */
    MAP *Pmap;			  /* pressure dof map */
    MAP *T_map;

    /*          | F  Bt|
     *  matNS = |      |
     *          | B  C |
     *  */
    MAT *matNS;			  /* matrix of the coupled problem  */
    MAT *matF;			  /* matNS[0][0] */
    MAT *matBt;			  /* matNS[0][1] */
    MAT *matB;			  /* matNS[1][0] */
    MAT *matC;			  /* matNS[1][1] */

    MAT *matNS0;
    MAT *matF0;			  /* matNS[0][0] */
    MAT *matBt0;			  /* matNS[0][1] */
    MAT *matB0;			  /* matNS[1][0] */
    MAT *matC0;			  /* matNS[1][1] */
    VEC *vec_rhs0;

    MAT *matT;			  /* Mat T, bottom free */
    VEC *rhsT;
    
    SOLVER *solver_u0;		  /* solver of the coupled problem */
    SOLVER *solver_u;		  /* solver of the coupled problem */
    SOLVER *solver_T;		  /* solver of conformation tensor */
    SOLVER *pc;			  /* preconditioner */
    NSPC *pcd;			  /* PCD preconditioner  */
#if USE_MG 
    MULTI_GRID *mg;		  /* Multi grid solver */
#endif /* USE_MG */

    /* Temp solver, Constrains */
    BOOLEAN *T_mask;		/* Temperature constrains mask */
    BOOLEAN *T_actc;		/* Temperature active constrain */
    FLOAT *T_cntr;		/* Temperature constrains value */

    /* Depth */
    DOF *depth_P1;		/* Depth P1 */
    DOF *depth_P2; 		/* Depth P2 */

    DOF *surf_elev_P1;
    DOF *bot_elev_P1;
    DOF *grad_surf_elev;

    /* Variables */
    FLOAT non_res;		  /* nonlinear residual */
    FLOAT *time;		  /* time */
    FLOAT time_queue[3];	  /* time queue */
    FLOAT *dt;			  /* time step ptr
				   * dt_{n} = t_{n+1} - t_{n} */
    FLOAT dt_queue[2];		  /* time step at different time */
    INT tstep;			  /* current time step */
    int viscosity_type;	          /* Viscosity types: const, T-indep, T-dep, ... */
    int set_dirichlet_bc; /* the surface bc for the inversion */
    LTYPE ltype;	          /* Picard or Newton */
    SURF_BAS *surf_bas;

    /* Layred mesh */
    LAYERED_MESH *gL;		  /* Layered mesh */
    /* MG_BLOCK_DOFS *bk;		  /\* Block dofs for line smoothing *\/ */
    /* DOF *coord; */

    NSParams *ns_params;	  /* Parameter list */
} NSSolver;




/*
 * ================================================================================
 * 
 *                 Global parameters
 * 
 * ================================================================================
 * */
extern NSParams *ns_params;

enum { VIS_CONST  = 0,
       VIS_STRAIN = 1,
       VIS_TEMP	  = 2};

extern FLOAT _Length_;
extern FLOAT _alpha_;

extern FLOAT nu_max;
extern FLOAT nu_min;

extern FLOAT eps_height;

typedef struct TIME_LOCK_ {
    BOOLEAN locked;
    FLOAT time;
} TIME_LOCK;

extern TIME_LOCK *time_lock; 


/*
 * ================================================================================
 * 
 *                 subroutines
 * 
 * ================================================================================
 * */
# define FUNC_T_DECLARE(func_xyz)					\
    void func_xyz##_t(FLOAT x, FLOAT y, FLOAT z, FLOAT t, FLOAT *values);    

/* ins-solver.c */
NSSolver *phgNSCreate(GRID *g, NSParams *ns_params);
void phgNSFinalize(NSSolver **ns);
void phgNSTimeAdvance(NSSolver *ns, FLOAT time, int tstep);
INT phgNSPinNode(NSSolver *ns);

void phgNSInitSolverU(NSSolver *ns);
void phgNSReInitSolverU(NSSolver *ns);
void phgNSBuildSolverUMat(NSSolver *ns, INT, INT, FLOAT);
//void phgNSBuildSolverUMat(NSSolver *ns);
//void phgNSBuildSolverURHS(NSSolver *ns, GEO_INFO *geo);
void phgNSBuildSolverURHS(NSSolver *ns, INT, INT, FLOAT);
//void phgNSBuildSolverURHS(NSSolver *ns);
void phgNSSolverUAssemble(NSSolver *ns);
void phgNSDestroySolverU(NSSolver *ns);

void getPecletNum(GRID *g, DOF *u, FLOAT nu, int order);
void phgDofMaxMin(DOF *u, FLOAT *umax, FLOAT *umin);
void estimate_error(NSSolver *ns, DOF *error);
void phgResumeLogUpdate(GRID *g, FLOAT *time, int *tstep, char *mesh_file, char *data_file);
void phgResumeStage(GRID *g, FLOAT *time, int *tstep, char *mesh_file, char *data_file);

/* save load */
void save_dof_data(GRID *g, DOF *dof, const char *file);
void load_dof_data(GRID *g, DOF *dof, const char *data_file, const char *mesh_file);
void load_dof_data2(GRID *g, DOF *dof, const char *data_file, const char *mesh_file);
void save_dof_data3(GRID *g, DOF *dof, const char *file);
void load_dof_data3(GRID *g, DOF *dof, const char *data_file, const char *mesh_file);
void ns_dof_copy(NSSolver *ns, DOF *u, DOF *p);

/* ins-bc.c */
int my_bc_map(int bctype);
void setFlowParameter(FLOAT Re_in, FLOAT nu_in, FLOAT time_in);
void setFuncTime(FLOAT time_in);
void adjust_time(FLOAT delta_time);
void restore_time(void);
void set_boundary_mask(NSSolver *ns);
void func_init_params(double Lx0, double alpha0);

void func_init(FLOAT x, FLOAT y, FLOAT z, FLOAT *u);
void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *u);
void func_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *p);
void func_gradp(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradp);
void func_gradu(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradu);
void func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *f);
void func_g1(FLOAT x, FLOAT y, FLOAT z, FLOAT *g);
void func_g2(FLOAT x, FLOAT y, FLOAT z, FLOAT *g);
void func_g3(FLOAT x, FLOAT y, FLOAT z, FLOAT *g);
void func_T(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
void func_beta(FLOAT x, FLOAT y, FLOAT z, FLOAT *beta);
void func_q(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
void func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord);
void func_ice_shelf_mask(FLOAT x, FLOAT y, FLOAT z, FLOAT *ice_shelf_mask);
void func_ice_shelf_pres(FLOAT x, FLOAT y, FLOAT z, FLOAT *ice_shelf_pres);
void func_ice_sur(FLOAT x, FLOAT y, FLOAT z, FLOAT *ice_sur);
void func_ice_bed(FLOAT x, FLOAT y, FLOAT z, FLOAT *ice_bed);
void func_ice_thk(FLOAT x, FLOAT y, FLOAT z, FLOAT *ice_thk);
void func_sur_grad_x(FLOAT x, FLOAT y, FLOAT z, FLOAT *sur_grad_x);
void func_sur_grad_y(FLOAT x, FLOAT y, FLOAT z, FLOAT *sur_grad_y);
void func_xyz_(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord);
void func_s(FLOAT x, FLOAT y, FLOAT z, FLOAT *s);
void func_b(FLOAT x, FLOAT y, FLOAT z, FLOAT *b);
void func_a(FLOAT x, FLOAT y, FLOAT *a);

FUNC_T_DECLARE(func_u);
FUNC_T_DECLARE(func_p);
FUNC_T_DECLARE(func_gradu);
FUNC_T_DECLARE(func_gradp);
FUNC_T_DECLARE(func_f);
FUNC_T_DECLARE(func_T);


/* ins-sovler.c */
NSParams *phgParametersCreate(void);
NSSolver *phgNSCreate(GRID *g, NSParams *ns_params0);
void phgNSFinalize(NSSolver **ns_ptr);
void phgNSTimeAdvance(NSSolver *ns, FLOAT time, int tstep);

void phgNSInitSolverU(NSSolver *ns);
//void phgNSBuildSolverURHS(NSSolver *ns, GEO_INFO *geo);
//void phgNSBuildSolverURHS(NSSolver *ns);
//void phgNSBuildSolverUMat(NSSolver *ns);
void phgNSDestroySolverU(NSSolver *ns);
void iceSetBoundaryTypes(NSSolver *ns);

/* temp-solver.c */
void phgNSInitSolverT(NSSolver *ns);
void phgNSBuildSolverT(NSSolver *ns);
void phgNSDestroySolverT(NSSolver *ns);
void phgNSBuildSolverTMat(NSSolver *ns, BOOLEAN init_T);
void phgNSBuildSolverTRHS(NSSolver *ns, BOOLEAN init_T);
void phgNSSolverTBuildConstrain(NSSolver *ns);
void find_melt_region(NSSolver *ns);
void phgNSSolverTSolve(NSSolver *ns, BOOLEAN init_T);
void phgNSTempInit(NSSolver *ns);
void proj_gradu(NSSolver *ns, DOF *gradu);


/* ins-pcd.c */
void phgNSInitPc(NSSolver *ns);
void phgNSBuildPc(NSSolver *ns);
void pc_proc(SOLVER *pc_solver, VEC *b0, VEC **x0);
void phgNSDestroyPc(NSSolver *ns);
void estimate_error(NSSolver *ns, DOF *error);

void phgDofSetName(DOF *dof, const char *name);
double elapsed_time(GRID *g, BOOLEAN flag, double mflops);
void phgDofRemoveConstMod(GRID *g, DOF *u);

/* ins-utils.c */
void NsSolver_Options();
FLOAT dofNormL2(DOF *dof);
void dof_norm_L2(DOF *dof);
void checkBdry(GRID *g);
int my_bc_map(int bctype);
void dof_range(DOF *u);
void output_bottom_dofs(NSSolver *ns, int tstep);

/* ins-mg.c */
void build_mg_levels(NSSolver *ns, int i_slv);

/* upwind.c */
void cr1_upwind(SIMPLEX *e, DOF *wind, FLOAT nu, 
		FLOAT Samarskij_alpha, int order, FLOAT *values);
void fv_upwind(SIMPLEX *e, DOF *wind, FLOAT nu, 
		FLOAT Samarskij_alpha, int order, FLOAT *values);

/* ice-grid.c */
//GEO_INFO *ice_grid(GRID *g);
void ice_grid(GRID *g);
FLOAT get_effective_viscosity(const FLOAT *gu, FLOAT T, FLOAT p,
			      BOOLEAN initialiszed);
void ice_monitor(NSSolver *ns, int nonstep);
BOOLEAN iceParter(GRID *g, int nprocs, DOF *weights, FLOAT power);
//GEO_INFO *iceInit(GRID *g, LAYERED_MESH **gL);
void iceInit(GRID *g, LAYERED_MESH **gL);

/* slip-bdry.c */
SURF_BAS *get_surface_bases(GRID *g, DOF_TYPE *u_type);
void trans_left(FLOAT *A, int ncol, int lda, const FLOAT *Trans); 
void trans_leftT(FLOAT *A, int ncol, int lda, const FLOAT *Trans); 
void trans_rightT(FLOAT *A, int ncol, int lda, const FLOAT *Trans); 
void rotate_dof_bases(DOF *u, SURF_BAS *surf_bas, BOOLEAN forward);
void dof_set_normal_data(DOF *u_h, SURF_BAS *surf_bas);
extern FLOAT trans_eye[Dim*Dim];


/* moving-mesh.c */
void get_surf_dH(NSSolver *ns);
void get_moved_coord(NSSolver *ns, int tstep);
void move_dof(GRID *g, DOF *dz, DOF *u);
void move_mesh(NSSolver *ns);
void get_layer_height(FLOAT *H, int nv, const FLOAT *ratio, FLOAT h0, FLOAT h1);
void get_height_depth(NSSolver *ns);

/* fv-solver.c */
void fv_solver_init(const char *node_file,
		    const char *trig_file, 
		    const DOF_USER_FUNC func_f);
void fv_update(const double *H, 
	       const double *U, 
	       double *dH, 
	       double *U_vert);

/* update_surf.c */
void struct_mesh_init(GRID *g);
void struct_mesh_reinit(GRID *g);
void struct_mesh_update(NSSolver *ns, int tstep, double t);




/*
 * ================================================================================
 * 
 *                 Utils
 * 
 * ================================================================================
 * */

#define OUTPUT_DIR "./OUTPUT/"
#define dH_OUTPUT_DIR "./dH_OUTPUT/"
#define MASK_OUTPUT_DIR "./MASK_OUTPUT/"
#define X_DIR 0
#define Y_DIR 1
#define Z_DIR 2
#define DDim (Dim*Dim)



/* Quad macros */
#define NbasFace(u) (3 * (u->type->np_vert + u->type->np_edge)	\
		     + u->type->np_face)
#define SQUARE(x) ((x)*(x))
#define INNER_PRODUCT(p, q)			\
    (*(p    ) * *(q    ) +			\
     *(p + 1) * *(q + 1) +			\
     *(p + 2) * *(q + 2))

#define BasisOrder(u, e, i) (!DofIsHP(u) ? (u)->type->order :		\
			      (u)->hp->info->types[(u)->hp->elem_order[e->index]]->order)
#define Bzero(v) bzero(v, sizeof(v));
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define DATA_FILE_SURFIX     {			\
  sprintf(data_u, "%s.u", data_file);		\
  sprintf(data_p, "%s.p", data_file);		\
  sprintf(data_T, "%s.T", data_file);		\
  }

#define GET_DOF_TYPE(dof_type, dof_name) {	\
    char s[128];				\
    phgOptionsPush();				\
    sprintf(s, "-dof_type %s", dof_name);	\
    phgOptionsSetOptions(s);			\
    dof_type = DOF_DEFAULT;			\
    phgOptionsPop();				\
  }

#define phgDofBY(b, y) {			\
    	INT i_, n_ = DofGetDataCount(*y);	\
	FLOAT *p_ = DofData(*y);		\
	for (i_ = 0; i_ < n_; i_++)		\
	    *(p_++) *= b;			\
    }



/* Periodic sync */
#define PERIODIC_SYNC(dof)					\
    if (g->period != NULL) {					\
	MAP *map_ = phgMapCreate(dof, NULL);			\
	VEC *vec_ = phgMapCreateVec(map_, 1);			\
	phgMapDofToLocalData(map_, 1, &dof, vec_->data);	\
	phgMapLocalDataToDof(map_, 1, &dof, vec_->data);	\
	phgVecDestroy(&vec_);					\
	phgMapDestroy(&map_);					\
    }


/* TODO: merge isop */
#define phgGeomGetCurvedJacobianAtLambda(g, e, p, det)	\
    *det = phgGeomGetVolume(g, e) * 6.;			\

#define phgQuadGetBasisCurvedGradient(e, u, i, quad, q) \
    phgQuadGetBasisGradient(e, u, i, quad) + q*Dim






/*
 * ================================================================================
 * 
 *                 Check macros
 * 
 * ================================================================================
 * */

#if 0
#define SHOW_M(matN, mat_m, mat_n) {				\
	printf("\n### rank: %d\n", g->rank);			\
	int i, j;						\
	printf(" --- "#matN":(%3d * %3d)\n", mat_m, mat_n);	\
	for (i = 0; i < mat_m; i++) {				\
	    for (j = 0; j < mat_n; j++){			\
		printf("%14.8f, ", *(matN + i * (mat_n) + j));	\
	    }							\
	    printf("\n");					\
	}							\
    }

#define SHOW_M3(matN, mat_k, mat_m, mat_n) {				\
	printf("\n### rank: %d\n", g->rank);				\
	int i, j, k;							\
	printf("--- --- %15s :(%3d * %3d)\n", #matN, mat_m, mat_n);	\
	for (k = 0; k < mat_k; k++) {					\
	    printf("  comp: %d\n", mat_k);				\
	    for (i = 0; i < mat_m; i++) {				\
		printf("    ");						\
		for (j = 0; j < mat_n; j++){				\
		    printf("%10f, ", *(matN + k * mat_m * mat_n +	\
				       i * mat_n + j));			\
		    if (mat_n > 10 && j%5 == 4)				\
			printf("\n    ");				\
		}							\
		printf("\n");						\
	    }								\
	}								\
    }

#define SHOW_V(vec, vec_n) { int _i;		\
	printf("\n### rank: %d\n", g->rank);	\
	printf(" --- "#vec":(%3d)\n", vec_n);	\
	for (_i = 0; _i < vec_n; _i++) {	\
	    printf("%10f, ", *(vec + _i));	\
	}					\
	printf("\n");				\
    }			
#endif

#define PRINT_ELAPSED_TIME(str, ...)	\
    phgPrintf(str);				\
    elapsed_time(__VA_ARGS__);
    

#if 1
#define DOF_SCALE(u, desp) {}
#elif 0
# define DOF_SCALE(u, desp)					\
    dof_range(u);
#elif 0
# define DOF_SCALE(u, desp)					\
    phgPrintf("    %s: [%16.8e, %16.8e]\n",			\
	      (u)->name,					\
	      phgDofMinValVec(u),				\
	      phgDofMaxValVec(u)				\
	      );						
#elif 0
# define DOF_SCALE(u, description) {				\
	char trimed[100];					\
	strncpy(trimed, __FUNCTION__, 8);			\
	trimed[8]='\0';						\
	phgPrintf("   ------------------------------\n"		\
		  "   %-10s  /* %s */\n"			\
		  "     func: %-10s, line: %03d\n"		\
		  ,#u":", description,				\
		  trimed, __LINE__);				\
	phgPrintf("     [%16.8e, %16.8e] (max,min); \n"		\
		  "     [%16.8e, %16.8e] ( L1, L2); \n",	\
		  phgDofMaxValVec(u), phgDofMinValVec(u),	\
		  phgDofNormL1(u), phgDofNormL2(u)		\
		  );						\
    }
#else
# define DOF_SCALE(u, description) {				\
	char trimed[100];					\
	DOF *_tmp_grad = phgDofGradient(u, NULL, NULL, NULL);	\
	strncpy(trimed, __FUNCTION__, 8);			\
	trimed[8]='\0';						\
	phgPrintf("   ------------------------------\n"		\
		  "   %-10s  /* %s */\n"			\
		  "     func: %-10s, line: %03d\n"		\
		  ,#u":", description,				\
		  trimed, __LINE__);				\
	phgPrintf("     [%16.8e, %16.8e] (max,min); \n"		\
		  "     [%16.8e, %16.8e] ( L1, L2); \n",	\
		  phgDofMaxValVec(u), phgDofMinValVec(u),	\
		  phgDofNormL1(u), phgDofNormL2(u)		\
		  );						\
	phgPrintf("   %-10s \n"					\
		  ,"grad("#u"):");				\
	phgPrintf("     [%16.8e, %16.8e] (max,min);\n"		\
		  "     [%16.8e, %16.8e] ( L1, L2);\n\n",	\
		  phgDofMaxValVec(_tmp_grad),			\
		  phgDofMinValVec(_tmp_grad),			\
		  phgDofNormL1(_tmp_grad),			\
		  phgDofNormL2(_tmp_grad)			\
		  );						\
	phgDofFree(&_tmp_grad);					\
    }
#endif

#if 1
#define sayHello(desp)							\
    {									\
	char hostname_[256];						\
	gethostname(hostname_, sizeof(hostname_));			\
	phgInfo(1, "Host %-15s, PID %10d, "				\
		"get line:%5d, mem: %0.4lfMB, /* %-20s */\n",		\
		hostname_, getpid(), __LINE__,				\
		phgMemoryUsage(g, NULL) / (1024.0 * 1024.0), desp);	\
	fflush(stdout);							\
	MPI_Barrier(g->comm);						\
    }
#else
#define sayHello(desp)
#endif


#define ERROR_MSG(desp) \
    phgError(1, "Error: file: %s, func: %s, line: %s\n"	\
	     "   %s not avaliable", __FILE__,		\
	     __FUNCTION__, __LINE__, desp);


#define INS_H
#endif

double** read_txt_data(char *file_name);

void get_mask_bot(NSSolver *ns);
void get_sur_normals_z(GRID *g, DOF *sur_normal_z);
void get_sur_normals_y(GRID *g, DOF *sur_normal_z);
void get_sur_normals_x(GRID *g, DOF *sur_normal_z);
void get_avg_n(GRID *g, DOF *sur_normal);
INT check_u_convergence(NSSolver *ns, DOF *u, DOF *p, DOF *u_last, DOF *p_last, FLOAT, FLOAT);
INT check_u_convergence0(NSSolver *ns, DOF *u, DOF *p, DOF *u_last, DOF *p_last, FLOAT, FLOAT);
void func_smb_top(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
void func_smb_bot(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
void func_smb_slf(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
void func_bed_z(FLOAT x, FLOAT y, FLOAT z, FLOAT *bed_z);
void update_grounding_line(NSSolver *ns, int tstep);
void get_water_pressure(NSSolver *);
void get_nodal_force(NSSolver *ns);
void get_nodal_force_value(NSSolver *ns);
void get_water_force(NSSolver *ns);
void get_contact_force(NSSolver *ns);
void get_contact_force_value(NSSolver *ns);
void get_contact(NSSolver *ns);
double* read_txt_data_1D(char *file_name);
void interp_txt_data_1D(double *data, FLOAT x, FLOAT y, FLOAT z, FLOAT *a);
INT if_update_shelf_mask(NSSolver *ns);
void get_stress(NSSolver *ns, DOF *gradu, DOF *pressure);
void get_stress1(NSSolver *ns);
void get_normal_stress_value(NSSolver *ns);
void get_water_pressure_value(NSSolver *ns);
void get_water_pressure_value1(NSSolver *ns);
INT check_surf_convergence(NSSolver *, DOF *, FLOAT);
FLOAT get_ice_volume(GRID *g);
void get_avg_gu(NSSolver *ns);
DOF * compare_two_dofs(DOF *a, DOF *b);
DOF * get_dof_component(GRID *g, DOF *a, DOF_TYPE *, INT dim_all, INT dim_asked);
void get_smooth_surface_values(NSSolver *ns, DOF *dof_P1, int up_or_lower);
void load_dH_from_file(NSSolver *ns, DOF *dof_P1, int up_or_lower);
void save_free_surface_velo(NSSolver *ns, int which_dim, int up_or_lower);
void save_free_surface_elev(NSSolver *ns, int up_or_lower);
void modify_mask_bot(NSSolver *ns);
void get_surf_bot_elev(NSSolver *ns);
void get_viscosity(NSSolver *ns);
void update_viscosity_inversion(NSSolver *ns);
INT check_visc_convergence(NSSolver *ns, DOF *visc_old, FLOAT tol);
