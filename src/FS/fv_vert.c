/*
 * 
 * Vertex-based Finite volume solver for conservation law.
 *
 * By Wei Leng, Lsec, 2012.
 *
 * Input:
 *   1) nodes coordinates
 *   2) triangle to nodes
 *
 * Output:
 *   error
 *
 * TODO:
 *   1) no boundary condition is implemented
 *
 * 
 *  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "mpi.h"

static int dbg = 0;


/* ------------------------------------------------------------
 * 
 *  Finite volume solver Interface
 *
 * ------------------------------------------------------------ */
typedef void (*DOF_USER_FUNC)(double x, double y, double z, double *values);

void fv_solver_init(const char *node_file,
		    const char *trig_file, 
		    const DOF_USER_FUNC func_f);
void fv_update(const double *U, 
	       const double *H, 
	       double *dH, 
	       double *U_vert);




typedef struct FV_ELEM_ {
    int v[3];
    //double x[3], y[3];
    //double xe[3], ye[3];
    double xc, yc;
    double len[3];
    double n[3][2];		/* normal */
    double a[3][2];		/* area */
    double quad_pt[3][2][2];	/* quad point
				 * [nedge][v0|1][X|Y] */
} FV_ELEM; 


typedef struct GRID_2D_ {
    int nv;
    int ne;
    FV_ELEM *tri;		/* [ne] */
    double *X;			/* [nv] */
    double *Y;
    double *node_area;		/* [nv] */
    char *bdry_mark;		/* [nv] */
    DOF_USER_FUNC func_f;	/* source function */
} GRID_2D; 


static GRID_2D g2D;

static int edge2vert[3][2] = {
    {0, 1}, {1, 2}, {2, 0}
};

/* static double *H; */
/* static double *dH; */
/* static double *H_old; */


/* boundary condition */
static double
func_h(double x, double y, double t)
{
    return 0.;
}

static double
get_tri_area(double x0, double x1, double x2, 
	     double y0, double y1, double y2, 
	     double *xc) 
{
    double a0 = x1 - x0;
    double a1 = x2 - x1;
    double a2 = x0 - x2;

    double b0 = y1 - y0;
    double b1 = y2 - y1;
    double b2 = y0 - y2;

    double l0 = sqrt(a0*a0 + b0*b0);
    double l1 = sqrt(a1*a1 + b1*b1);
    double l2 = sqrt(a2*a2 + b2*b2);

    double s = (l0+l1+l2)/2; 
    double area = sqrt(s*(s-l0)*(s-l1)*(s-l2));

    if (xc != NULL) {
	xc[0]= (x0 + x1 + x2) / 3.;
	xc[1]= (y0 + y1 + y2) / 3.;
    }

    return area;
}


static void
get_solution(double *solu, double time, int on_bdry)
{
    FV_ELEM *fv_tri = g2D.tri;
    int nv = g2D.nv;
    int ne = g2D.ne;
    double *X = g2D.X, *Y = g2D.Y;
    int i, j, k; 


#if 0
#   warning u mean int
    if (on_bdry) {
	for (i = 0; i < nv; i++)
	    if (g2D.bdry_mark[i])
		solu[i] = 0;	/* clear bdry */
    } else {
	bzero(solu, nv * sizeof(double));
    }


    for (i = 0; i < ne; i++) {
	FV_ELEM *tri = &fv_tri[i];

	int V[3] = {tri.v[0], tri.v[1], tri.v[2]};
	double x[3] = {X[V[0]], X[V[1]], X[V[2]]};
	double y[3] = {Y[V[0]], Y[V[1]], Y[V[2]]};

	if (on_bdry && 
	    !(g2D.bdry_mark[V[0]] ||
	      g2D.bdry_mark[V[1]] ||
	      g2D.bdry_mark[V[2]]))
	    continue;

	for (k = 0; k < 3; k++) {
	    int v0 = edge2vert[k][0];
	    int v1 = edge2vert[k][1];
	    int V0 = tri.v[v0];
	    int V1 = tri.v[v1];

	    double xe = (x[v0] + x[v1]) / 2.;
	    double ye = (y[v0] + y[v1]) / 2.;

	    double xc = tri.xc;
	    double yc = tri.yc;

	    double a0 = tri.a[k][0];
	    double a1 = tri.a[k][1];

	    // x[v0], y[v0];
	    // x[v1], y[v1];

	    double xx, yy;
	    if (!on_bdry || g2D.bdry_mark[V0]) {
		xx = (xc + xe + x[v0]) / 3;
		yy = (yc + ye + y[v0]) / 3;
		solu[V0] += a0 * func_h(xx, yy, time);
	    }

	    if (!on_bdry || g2D.bdry_mark[V1]) {
		xx = (xc + xe + x[v1]) / 3;
		yy = (yc + ye + y[v1]) / 3;
		solu[V1] += a1 * func_h(xx, yy, time);
	    }
	}
    }

    for (i = 0; i < nv; i++) {
	if (!on_bdry || g2D.bdry_mark[i])
	    solu[i] /= g2D.node_area[i];
    }
    
#else
#   warning u mean v0
    for (i = 0; i < nv; i++) {
	if (!on_bdry || g2D.bdry_mark[i])
	    solu[i] = func_h(X[i], Y[i], time);
    }
#endif
}




void
fv_solver_init(const char *node_file,
	       const char *trig_file, 
	       DOF_USER_FUNC func_f)
{
    FILE *fp;
    int i, j, k, id;
    int nv, ne;
    double *X, *Y;
    int nbdry = 0;
    int rank; 
    g2D.func_f = func_f;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* nodes */
    fp = fopen(node_file, "r");
    fscanf(fp, "%d", &nv);
    X = calloc(nv, sizeof(*X));
    Y = calloc(nv, sizeof(*Y));
    g2D.nv = nv;
    g2D.X = X;
    g2D.Y = Y;
    g2D.bdry_mark = calloc(nv, sizeof(*g2D.bdry_mark));
    for (i = 0; i < nv; i++) {
	int mark;
	fscanf(fp, "%d %lf %lf %d", &id,
	       X + i, Y + i, &mark);
	g2D.bdry_mark[i] = (char) mark;
	if (g2D.bdry_mark[i])
	    nbdry ++;
    }
    fclose(fp);


    /* triangles */
    fp = fopen(trig_file, "r");
    fscanf(fp, "%d", &ne);
    FV_ELEM *fv_tri = (FV_ELEM *) calloc(ne, sizeof(*fv_tri));
    g2D.ne = ne;
    g2D.tri = fv_tri;
    for (i = 0; i < ne; i++) {
	int v0, v1, v2;
	fscanf(fp, "%d %d %d %d", &id,
	       &v0, &v1, &v2);
	fv_tri[i].v[0] = v0 - 1;
	fv_tri[i].v[1] = v1 - 1;
	fv_tri[i].v[2] = v2 - 1;
    }
    fclose(fp);

    if (rank == 0)
	printf("   Triangle mesh nv:%d ne:%d, nb:%d\n", nv, ne, nbdry);



    /*
     *
     * Bulid FV solver
     *
     * */
    double *node_area = calloc(nv, sizeof(double));
    g2D.node_area = node_area;
    double min_h = 1e20, max_h = -1e20;

    for (i = 0; i < ne; i++) {
	FV_ELEM *tri = &fv_tri[i];
	
	int V[3] = {tri->v[0], tri->v[1], tri->v[2]};
	assert(V[0] < nv);
	assert(V[1] < nv);
	assert(V[2] < nv);

	/* printf("elem %5d %5d %5d %5d\n", i, */
	/*        V[0], V[1], V[2]); */

	/* center */
	double x[3] = {X[V[0]], X[V[1]], X[V[2]]};
	double y[3] = {Y[V[0]], Y[V[1]], Y[V[2]]};

	double xc = (x[0] + x[1] + x[2]) / 3.;
	double yc = (y[0] + y[1] + y[2]) / 3.;

	tri->xc = xc;
	tri->yc = yc;
	
	{
	    double a0 = x[1] - x[0];
	    double a1 = x[2] - x[1];
	    double a2 = x[0] - x[2];

	    double b0 = y[1] - y[0];
	    double b1 = y[2] - y[1];
	    double b2 = y[0] - y[2];

	    double l0 = sqrt(a0*a0 + b0*b0);
	    double l1 = sqrt(a1*a1 + b1*b1);
	    double l2 = sqrt(a2*a2 + b2*b2);
	    double l[3] = {l0, l1, l2};

	    for (j = 0; j < 3; j++) {
		if (max_h < l[j])
		    max_h = l[j];
		if (min_h > l[j])
		    min_h = l[j];
	    }
	}

	/* for edges */
	for (k = 0; k < 3; k++) {
	    int v0 = edge2vert[k][0];
	    int v1 = edge2vert[k][1];
	    int V0 = tri->v[v0];
	    int V1 = tri->v[v1];

	    double xe = (x[v0] + x[v1]) / 2.;
	    double ye = (y[v0] + y[v1]) / 2.;

	    double xx = xe - xc;
	    double yy = ye - yc;

	    double len = sqrt(xx*xx + yy*yy);
	    double n[2] = {yy, -xx};
	    n[0] /= len; n[1] /= len;

	    double t[2] = {x[v1] - x[v0], y[v1] - y[v0]};
	    if (n[0]*t[0] + n[1]*t[1] < 0) {
		n[0] *= -1; n[1] *= -1;
	    }
	    
	    double a0 = get_tri_area(xc, x[v0], xe,
				     yc, y[v0], ye,
				     tri->quad_pt[k][0]);
	    double a1 = get_tri_area(xc, x[v1], xe,
				     yc, y[v1], ye, 
				     tri->quad_pt[k][1]);

	    tri->len[k] = len;
	    tri->n[k][0] = n[0];
	    tri->n[k][1] = n[1];
	    tri->a[k][0] = a0;
	    tri->a[k][1] = a1;

	    node_area[V0] += a0;
	    node_area[V1] += a1;
	}
    }

    if (rank == 0)
	printf("   Triangle mesh h: [%e %e]\n", min_h, max_h);


    return;
}



/*
 * FV solver update
 * Input: U [ne][3][2], 
 *        H [nv]
 * Output: dH [nv]
 * */
void
fv_update(const double *H, 
	  const double *U, 
	  double *dH, 
	  double *U_vert)
{
    FV_ELEM *fv_tri = g2D.tri;
    int nv = g2D.nv;
    int ne = g2D.ne;
    double *X = g2D.X, *Y = g2D.Y;
    double *node_area = g2D.node_area;
    int i, j, k; 
    DOF_USER_FUNC func_f = g2D.func_f;

    double *H_old = calloc(nv, sizeof(double));
    memcpy(H_old, H, nv * sizeof(double));
    bzero(dH, nv * sizeof(double));

    if (U_vert != NULL) {
	bzero(U_vert, nv * sizeof(double));
    }


    for (i = 0; i < ne; i++) {
	FV_ELEM *tri = &fv_tri[i];
	
	if (dbg) printf("\n\nelem: %d\n", i);
	int V[3] = {tri->v[0], tri->v[1], tri->v[2]};

	/* Compute flux */
	for (k = 0; k < 3; k++) {
	    if (dbg) printf("edge: %d\n", k);

	    int v0 = edge2vert[k][0];
	    int v1 = edge2vert[k][1];
	    int V0 = tri->v[v0];
	    int V1 = tri->v[v1];

	    double *a = tri->a[k];
	    double *n = tri->n[k];
	    double len = tri->len[k];
	    double h = 0;
	    double *qpt = &tri->quad_pt[k][0][0];

	    double vel_u = U[i*6 + 2*k    ];
	    double vel_v = U[i*6 + 2*k + 1];

#if 0
#   warning Flux central
	    /* Central */
	    h = .5 * (H_old[V0] + H_old[V1]);
#else
#   warning Flux upwind
	    /*
	     * 
	     * Hpwind: v0 -> v1
	     * 
	     * */
	    if (n[0] * vel_u + n[1] * vel_v > 0) {
		h = H_old[V0];
	    } else {
		h = H_old[V1];
	    }
#endif

	    double flux = len * (n[0] * vel_u + n[1] * vel_v) * h;

	    dH[V0] -= flux / node_area[V0];
	    dH[V1] += flux / node_area[V1];

	    /* source term */
	    if (func_f != NULL) {
		double f0, f1;
		func_f(qpt[0], qpt[1], 0, &f0);
		func_f(qpt[2], qpt[3], 0, &f1);

		dH[V0] += f0 * a[0] / node_area[V0];
		dH[V1] += f1 * a[1] / node_area[V1];
	    }
	    

	    /* velocity at vert */
	    if (U_vert != NULL) {
		U_vert[2*V0   ] += vel_u *  a[0] / node_area[V0];
		U_vert[2*V0 +1] += vel_v *  a[0] / node_area[V0];
		U_vert[2*V1   ] += vel_u *  a[1] / node_area[V1];
		U_vert[2*V1 +1] += vel_v *  a[1] / node_area[V1];
	    }

	    if (dbg) printf("%lf %lf %lf %lf %lf\n",
			    len, n[0], n[1], dH[V0], dH[V1]);
	} /* end flux */
    }     /* end elem */

	
    /* Set boundary */
    get_solution(dH, 0., 1);

#if 0
    fp = fopen("H_.m", "w");
    fprintf(fp, "H=[\n");
    for (i = 0; i < nv; i++)
	fprintf(fp, "%lf\n", H[i]);
    fprintf(fp, "];\n");
    fclose(fp);
#endif

    free(H_old);
    return;
}
