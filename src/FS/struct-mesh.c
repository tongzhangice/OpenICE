#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "phg.h"
#include "ins.h"

/*
 * Vertex index on the surface
 *
 * 2D (i, j)  in [0..Nx] X [0..Ny] 
 * 3D index (Nz+1) * ((i * (Ny+1) + j)) + Nz
 *
 *
 *  */

static int Nx = 0, Ny = 0, Nz = 0;

static double *x_dat = NULL,
    *u_dat = NULL, *u1_dat = NULL,
    *s_dat = NULL, *ds_dat = NULL, *ds2_dat = NULL,
    *dh_dat = NULL, *dh2_dat = NULL; 
static double dx = 0., dy = 0.;
static int *map = NULL;


/*
 * Point Data X, U, ... [Nx+1][Nx+1], 
 *   X, S, is not periodical, (X[0] != X[N]) 
 *   U, dS, is periodical, (U[0] == U[N])
 *
 * Cell Data dH, ... [Nx][Nx],
 *   X, S, is not periodical, (dH[0] != dH[N-1]) 
 *
 * Note:
 *   Keep dS periodical.
 */
#define X(i, j, k) (get_data(x_dat, 3, i, j) + k)
#define U(i, j, k) (get_data(u_dat, 3, i, j) + k)
#define U1(i, j, k) (get_data(u1_dat, 3, i, j) + k)
#define S(i, j) get_data(s_dat, 1, i, j)
#define DS(i, j) get_data(ds_dat, 1, i, j)
#define DS2(i, j) get_data(ds2_dat, 1, i, j)

#define DH(i, j) get_data2(dh_dat, 1, i, j)
#define DH2(i, j) get_data2(dh2_dat, 1, i, j)

static double *
get_data(double *dat, int dim, int i, int j)
/* length: N+1 */
{
    while (i >= Nx+1) 
	i -= Nx+1;
    while (i < 0)
	i += Nx+1;

    while (j >= Ny+1) 
	j -= Ny+1;
    while (j < 0)
	j += Ny+1;

    assert(i >= 0 && i < Nx+1);
    assert(j >= 0 && j < Ny+1);
    return dat + (i*(Ny+1) + j) * dim;
}

static double *
get_data2(double *dat, int dim, int i, int j)
/* length: N */
{
    while (i >= Nx) 
	i -= Nx;
    while (i < 0)
	i += Nx;

    while (j >= Ny) 
	j -= Ny;
    while (j < 0)
	j += Ny;

    assert(i >= 0 && i < Nx);
    assert(j >= 0 && j < Ny);
    return dat + (i*(Ny) + j) * dim;
}




void
struct_mesh_init(GRID *g)
/*
 * Note: call this function before redistribute !!!
 *  */
{
    int rank = phgRank;
    int i, j, k, p;
    
    Nx = ns_params->Nx;
    Ny = ns_params->Ny;
    Nz = ns_params->Nz;

    printf("init surf[%d]: %d, %d, %d\n", rank, Nx, Ny, Nz);


    /* Alloc data */
    x_dat = calloc(3 * (Nx+1) * (Ny+1), sizeof(*x_dat)) ;
    u_dat = calloc(3 * (Nx+1) * (Ny+1), sizeof(*u_dat)) ;
    u1_dat = calloc(3 * (Nx+1) * (Ny+1), sizeof(*u_dat)) ;
    s_dat = calloc((Nx+1) * (Ny+1), sizeof(*s_dat)) ;
    ds_dat = calloc((Nx+1) * (Ny+1), sizeof(*ds_dat)) ;
    ds2_dat = calloc((Nx+1) * (Ny+1), sizeof(*ds_dat)) ;
    dh_dat = calloc((Nx) * (Ny), sizeof(*dh_dat));
    dh2_dat = calloc((Nx) * (Ny), sizeof(*dh_dat));

    /* Init data */
    if (rank == 0) {
	assert(g->nvert == (Nx+1)*(Ny+1)*(Nz+1));

	for (i = 0; i < Nx+1; i++) 
	    for (j = 0; j < Ny+1; j++) {
		int idx = (Nz+1) * (i * (Ny+1) + j) + Nz;
		assert(g->types_vert[idx] & BC_TOP);

		/* init X and S */
		for (k = 0; k < 3; k++) 
		    *X(i, j, k) = g->verts[idx][k];
		*S(i, j) = *X(i, j, 2);
	    }
    }

    MPI_Bcast(x_dat, (Nx+1)*(Ny+1) * 3,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(s_dat, (Nx+1)*(Ny+1) ,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);

    dx = *X(1, 0, 0) - *X(0, 0, 0);
    dy = *X(0, 1, 1) - *X(0, 0, 1);
    printf("(%d, %d, %d) dx: %e, dy: %e\n", Nx, Ny, Nz, dx, dy);

    return;
}

void
struct_mesh_reinit(GRID *g)
/*
 * Note: call this function before redistribute !!!
 *  */
{
    int rank = phgRank;
    int i, j, k, l, p;

    printf("re init surf[%d]: %d, %d, %d\n", rank, Nx, Ny, Nz);

    /* Re Init data */
    for (i = 0; i < Nx+1; i++) 
	for (j = 0; j < Ny+1; j++) 
	    *S(i, j) = -1e30;

    for (p = 0; p < g->nvert; p++) {
	if (g->types_vert[p] & BC_TOP) {
	    /* idx => (i, j, l) */
	    int idx = GlobalVertex(g, p);
	    l = idx % (Nz+1);
	    idx /= (Nz+1);
	    j = idx % (Ny+1);
	    idx /= (Ny+1);
	    i = idx;
	    assert((Nz+1) * (i * (Ny+1) + j) + l == GlobalVertex(g, p));
	    idx = GlobalVertex(g, p); 

	    assert(i >= 0 && i < Nx+1);
	    assert(j >= 0 && j < Ny+1);
	    assert(l >= 0 && l < Nz+1);
	    assert(l == Nz);

	    *S(i, j) = g->verts[p][2];
	}
    }

    MPI_Allreduce(s_dat, ds_dat, (Nx+1)*(Ny+1) ,
		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    memcpy(s_dat, ds_dat, (Nx+1)*(Ny+1) * sizeof(*s_dat));

    for (i = 0; i < Nx+1; i++) 
	for (j = 0; j < Ny+1; j++) 
	    *U(i, j, 2) = *S(i, j);

    return;
}


/* ------------------------------------------------------------
 *
 * Get surface change: use kinematic relation
 *
 * ------------------------------------------------------------ */
static void
get_dS1(NSSolver *ns, int tstep)
{
    int rank = phgRank;
    GRID *g = ns->u[1]->g;
    SIMPLEX *e;
    double dt = ns->dt[0];
    int i, j, k, p, l;

    phgPrintf("Update surface: kinematic.\n");

    /* 1. Compute u, v, w on nodes.  */
    for (i = 0; i < Nx+1; i++) 
	for (j = 0; j < Ny+1; j++) {
	    for (k = 0; k < 3; k++)
		*U(i, j, k) = -1e30;
	}

    /* get U */
    for (p = 0; p < g->nvert; p++) {
	if (g->types_vert[p] & BC_TOP) {
	    /* idx => (i, j, l) */
	    int idx = GlobalVertex(g, p);
	    l = idx % (Nz+1);
	    idx /= (Nz+1);
	    j = idx % (Ny+1);
	    idx /= (Ny+1);
	    i = idx;
	    assert((Nz+1) * (i * (Ny+1) + j) + l == GlobalVertex(g, p));
	    idx = GlobalVertex(g, p); 

	    assert(i >= 0 && i < Nx+1);
	    assert(j >= 0 && j < Ny+1);
	    assert(l >= 0 && l < Nz+1);
	    assert(l == Nz);

	    for (k = 0; k < 3; k++) 
		*U(i, j, k) = ns->u[1]->data[p*3 + k];
	}
    }

    MPI_Allreduce(u_dat, u1_dat, (Nx+1)*(Ny+1)*3, /* u->u1 */
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    memcpy(u_dat, u1_dat, (Nx+1)*(Ny+1)*3 * sizeof(*u_dat)); /* u1->u */
    

    /* compute surface change */
    if (rank == 0) {
	double ds_max = -1e10, ds_min = 1e10;
	double eu_max = -1e10, eu_min = 1e10;
	char filename[1000];
	FILE *fp;

	sprintf(filename, "./output/s_%05d.dat", tstep);
	fp = fopen(filename, "w");

#if 0
	/* Smoothing 3 times. */
	int ismooth;
	for (ismooth = 0; ismooth < 3; ismooth++) {
	    for (i = 0; i < Nx+1; i++)
		for (j = 0; j < Nx+1; j++) {
		    for (k = 0; k < 3; k++) {
			/* Five point stencil */
			*U1(i, j, k) = .2 * (*U(i+1, j, k) + *U(i-1, j, k) + *U(i, j+1, k) + *U(i, j-1, k)
					     + *U(i,j, k) );
		    }
		}
	    memcpy(u_dat, u1_dat, 3*(Nx+1)*(Ny+1)*sizeof(*u_dat));
	}
#endif


    
    
	/* 2. Update height */
        /* * s_{n+1} = s_{n} + dt * (u * ds/dx + v * ds/dy - w + a)  */
	for (i = 0; i < Nx+1; i++) {
	    for (j = 0; j < Ny+1; j++) {
		double x = *X(i, j, 0);
		double y = *X(i, j, 1);
		double z = *S(i,j);
		double u0 = *U(i,j, 0), v0 = *U(i,j, 1), w0 = *U(i,j, 2);
		double dsx, dsy, exact_s, da;
#if 0
		/* central */
		//dsx = (*S(i+1, j) - *S(i-1, j)) / (*X(i+1, j, 0) - *X(i-1, j, 0));
		//dsy = (*S(i, j+1) - *S(i, j-1)) / (*X(i, j+1, 1) - *X(i, j-1, 1));
		dsx = (*S(i+1, j) - *S(i-1, j)) / (2*dx);
		dsy = (*S(i, j+1) - *S(i, j-1)) / (2*dy);
#elif 0
		/* upwind */
		if (u0 > 0) 
		    //dsx = (*S(i, j) - *S(i-1, j)) / (*X(i, j, 0) - *X(i-1, j, 0));
		    dsx = (*S(i, j) - *S(i-1, j)) / dx;
		else
		    //dsx = (*S(i+1, j) - *S(i, j)) / (*X(i+1, j, 0) - *X(i, j, 0));
		    dsx = (*S(i+1, j) - *S(i, j)) / dx;
		
		if (v0 > 0) 
		    //dsy = (*S(i, j) - *S(i, j-1)) / (*X(i, j, 1) - *X(i, j-1, 1));
		    dsy = (*S(i, j) - *S(i, j-1)) / dy;
		else
		    //dsy = (*S(i, j+1) - *S(i, j)) / (*X(i, j+1, 1) - *X(i, j, 1));
		    dsy = (*S(i, j+1) - *S(i, j)) / dy;
#else
		if (i == Nx) 
		    dsx = .5 * (*S(Nx, j) - *S(Nx-1, j)) / dx
			+ .5 * (*S(1, j) - *S(0, j)) / dx;
		else if (i == 0)
		    dsx = .5 * (*S(Nx, j) - *S(Nx-1, j)) / dx
			+ .5 * (*S(1, j) - *S(0, j)) / dx;
		else
		    dsx = (*S(i+1, j) - *S(i-1, j)) / (2*dx);
		
		if (j == Ny) 
		    dsy = .5 * (*S(i, Ny) - *S(i, Ny-1)) / dy
			+ .5 * (*S(i, 1) - *S(i, 0)) / dy;
		else if (j == 0)
		    dsy = .5 * (*S(i, Ny) - *S(i, Ny-1)) / dy
			+ .5 * (*S(i, 1) - *S(i, 0)) / dy;
		else
		    dsy = (*S(i, j+1) - *S(i, j-1)) / (2*dy);
#endif

		func_a(x*1e3, y*1e3, &da);
		double ds0= u0 * dsx + v0 * dsy - w0 + da;
		double eu = u0 * dsx + v0 * dsy - w0;

		*DS(i, j) = ds0 - da;
		
		if (ds0 > ds_max)
		    ds_max = ds0;
		if (ds0 < ds_min)
		    ds_min = ds0;

		if (eu > eu_max)
		    eu_max = eu;
		if (eu < eu_min)
		    eu_min = eu;

		func_s(x*1e3, y*1e3, 0, &exact_s); 
		exact_s /= 1e3;			   /* km */

		fprintf(fp, "%e %e %e " /* 1-3 */
			"%e %e "	/* 4-5 */
			"%e %e %e %e %e "	/* 6-8 9-10 */
			"%e\n",	/* 11-12 */
			x, y, z, 
			exact_s, z - exact_s,
			u0, v0, w0, dsx, dsy,
			u0 * dsx + v0 * dsy - w0);
	    }
	    fprintf(fp, "\n");
	}


	printf("dh: [%e %e]\n", ds_min/1e3*dt, ds_max/1e3*dt);
	printf("eu: [%e %e]\n", eu_min/1e3*dt, eu_max/1e3*dt);
	fclose(fp);

#if 0
	/* 	/\* Smoothing 3 times. *\/ */
	for (ismooth = 0; ismooth < 3; ismooth++) {
	    for (i = 0; i < Nx+1; i++)
		for (j = 0; j < Nx+1; j++) {
		    /* Four point stencil */
		    *DS2(i, j) = .2 * (*DS(i+1, j) + *DS(i-1, j) + *DS(i, j+1) + *DS(i, j-1)
				       + *DS(i,j));
		}
	    memcpy(ds_dat, ds2_dat, (Nx+1)*(Ny+1)*sizeof(*s_dat));
	}
#else
#endif


    }	/* end root */
}


/* ------------------------------------------------------------
 *
 * Get surface change: use non-conservative scheme
 *
 * ------------------------------------------------------------*/
static void
get_dS2(NSSolver *ns, int tstep)
{
    int rank = phgRank;
    GRID *g = ns->u[1]->g;
    SIMPLEX *e;
    int i, j, k, p, idx, l, s;
    double dt = ns->dt[0];

    phgPrintf("Update surface: non conservative.\n");
    
    /* 1. Compute u, v, w on nodes.  */
    for (i = 0; i < Nx+1; i++) 
	for (j = 0; j < Ny+1; j++) {
	    for (k = 0; k < 3; k++)
		*U(i, j, k) = -1e30;
	}

    /* get U */
    for (p = 0; p < g->nvert; p++) {
	if (g->types_vert[p] & BC_TOP) {
	    /* idx => (i, j, l) */
	    int idx = GlobalVertex(g, p);
	    l = idx % (Nz+1);
	    idx /= (Nz+1);
	    j = idx % (Ny+1);
	    idx /= (Ny+1);
	    i = idx;
	    assert((Nz+1) * (i * (Ny+1) + j) + l == GlobalVertex(g, p));
	    idx = GlobalVertex(g, p); 

	    assert(i >= 0 && i < Nx+1);
	    assert(j >= 0 && j < Ny+1);
	    assert(l >= 0 && l < Nz+1);
	    assert(l == Nz);

	    for (k = 0; k < 3; k++) 
		*U(i, j, k) = ns->u[1]->data[p*3 + k];
	}
    }

    MPI_Allreduce(u_dat, u1_dat, (Nx+1)*(Ny+1)*3, /* u->u1 */
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    memcpy(u_dat, u1_dat, (Nx+1)*(Ny+1)*3 * sizeof(*u_dat)); /* u1->u */


    /* Compute surface change */
    if (phgRank == 0) {
	double ds_max = -1e10, ds_min = 1e10;
	char filename[1000];
	FILE *fp;

	sprintf(filename, "./output/u_%05d.dat", tstep);
	fp = fopen(filename, "w");


	/* for cell */
	for (i = 0; i < Nx; i++) {
	    for (j = 0; j < Ny; j++) {
		double u[3] = {0, 0, 0};
		for (k = 0; k < 3; k++)
		    u[k] = (*U(i, j, k) + *U(i, j+1, k)
			    + *U(i+1, j, k) + *U(i+1, j+1, k)) / 4.;
		double dsx = .5 * (*S(i+1, j) - *S(i, j)) / dx
		    + .5 * (*S(i+1, j+1) - *S(i, j+1)) / dx;
		double dsy = .5 * (*S(i, j+1) - *S(i, j)) / dy
		    + .5 * (*S(i+1, j+1) - *S(i+1, j)) / dy;

		*DH(i, j) = u[0]*dsx + u[1]*dsy - u[2];
		fprintf(fp, "%e %e "
			"%e %e %e"
			" %e %e\n",
			*X(i, j, 0), 
			*X(i, j, 1), 
			u[0], u[1], u[2], dsx, dsy);
	    }
	    fprintf(fp, "\n");
	}
	fclose(fp);


	/* Smoothing 3 times. */
	int ismooth;
	for (ismooth = 0; ismooth < 3; ismooth++) {
	    for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++) {
		    /* Four point stencil */
		    *DH2(i, j) = .2 * (*DH(i+1, j) + *DH(i-1, j) + *DH(i, j+1) + *DH(i, j-1)
				       + *DH(i,j));
		}
	    memcpy(dh_dat, dh2_dat, (Nx)*(Ny)*sizeof(*dh_dat));
	}


	sprintf(filename, "./output/s_%05d.dat", tstep);
	fp = fopen(filename, "w");

        for (i = 0; i < Nx+1; i++) {
            for (j = 0; j < Ny+1; j++) {
		double x = *X(i, j, 0);
		double y = *X(i, j, 1);
		double z = *S(i,j);
		double exact_s, ds0;

                *DS(i, j) = (*DH(i-1, j-1) + *DH(i, j) + *DH(i-1, j) + *DH(i, j-1)) / 4.;

		ds0 = *DS(i, j); /* m/a */

		if (ds0 > ds_max)
		    ds_max = ds0;
		if (ds0 < ds_min)
		    ds_min = ds0;

		func_s(x*1e3, y*1e3, 0, &exact_s); 
		exact_s /= 1e3;			   /* km */
		fprintf(fp, "%e %e %e " /* 1-3 */
			"%e %e "	/* 4-5 */
			"%e\n",	/* 6-7 */
			x, y, z, 
			exact_s, z - exact_s,
			ds0);
	    }
	    fprintf(fp, "\n");
        }

	printf("dh: [%e %e]\n", ds_min/1e3*dt, ds_max/1e3*dt);
	fclose(fp);
    }

}


/* ------------------------------------------------------------
 *
 * Get surface change: use conservative scheme
 *
 * ------------------------------------------------------------*/
static void
get_dS3(NSSolver *ns, int tstep)
{
    int rank = phgRank;
    GRID *g = ns->u[1]->g;
    SIMPLEX *e;
    int i, j, k, idx, l, s;
    double dt = ns->dt[0];

    phgPrintf("Update surface: conservative.\n");
    
    /* Compute flux */
    memset(dh_dat, 0, Nx*Ny*sizeof(*dh_dat));
    memset(dh2_dat, 0, Nx*Ny*sizeof(*dh_dat));

    ForAllElements(g, e) {
	for (s = 0; s < NFace; s++) {
	    if (g->types_face[e->faces[s]] & BC_TOP) {
		int ie = 100000, je = 100000;		
		for (l = 0; l < 3; l++) {
		    idx = GlobalVertex(g, e->verts[GetFaceVertex(s, l)]);
		    k = idx % (Nz+1);
		    idx /= (Nz+1);
		    j = idx % (Ny+1);
		    idx /= (Ny+1);
		    i = idx;

		    /* top face on cube(ie, je, Nz-1) */
		    assert(k == Nz);
		    ie = (ie < i) ? ie : i;
		    je = (je < j) ? je : j;
		}
		i = j = k = -1;

		QUAD *quad;
		FLOAT flux, area, lambda[Dim + 1], vu[Dim];
		int v0, v1, v2, q;
		const FLOAT *w, *p, *normal;

		flux = 0;
		quad = phgQuadGetQuad2D(4);
		area = phgGeomGetFaceArea(g, e, s);
		normal = phgGeomGetFaceOutNormal(g, e, s);

		v0 = GetFaceVertex(s, 0);
		v1 = GetFaceVertex(s, 1);
		v2 = GetFaceVertex(s, 2);
		lambda[s] = 0.;

		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    lambda[v0] = *(p++);
		    lambda[v1] = *(p++);
		    lambda[v2] = *(p++);
		    phgDofEval(ns->u[1], e, lambda, vu);
		    
		    flux += area*(*w) * INNER_PRODUCT(vu, normal);
		    w++;
		} 

		*DH(ie, je) += flux / dx / dy; /* Unit: m/a */
	    } /* end top face */
	}     /* end face */
    }	      /* end elem */
    

    MPI_Allreduce(dh_dat, dh2_dat, Nx*Ny, MPI_DOUBLE, MPI_SUM, g->comm);
    memcpy(dh_dat, dh2_dat, (Nx)*(Ny) * sizeof(*dh_dat));


    /* Compute surface change */
    if (phgRank == 0) {
	double ds_max = -1e10, ds_min = 1e10;
	char filename[1000];
	FILE *fp;

	sprintf(filename, "./output/s_%05d.dat", tstep);
	fp = fopen(filename, "w");

        for (i = 0; i < Nx+1; i++) {
            for (j = 0; j < Ny+1; j++) {
		double x = *X(i, j, 0);
		double y = *X(i, j, 1);
		double z = *S(i,j);
		double exact_s = 0, ds0 = 0;

                *DS(i, j) = (*DH(i-1, j-1) + *DH(i, j) + *DH(i-1, j) + *DH(i, j-1)) / 4.;
		ds0 = *DS(i, j); /* m/a */

		if (ds0 > ds_max)
		    ds_max = ds0;
		if (ds0 < ds_min)
		    ds_min = ds0;

		func_s(x*1e3, y*1e3, 0, &exact_s); 
		exact_s /= 1e3;			   /* km */
		fprintf(fp, "%e %e %e " /* 1-3 */
			"%e %e "	/* 4-5 */
			"%e\n",	/* 6-7 */
			x, y, z, 
			exact_s, z - exact_s,
			ds0);
	    }
	    fprintf(fp, "\n");
        }

	printf("dh: [%e %e]\n", ds_min/1e3*dt, ds_max/1e3*dt);
	fclose(fp);
    }

}






 
void
struct_mesh_update(NSSolver *ns, int tstep, double t)
{
    GRID *g = ns->g;
    DOF *dof_u = ns->u[1];
    int rank = g->rank;
    int i, j, k, l, p, ismooth;
    double dt = ns->dt[0];

    phgPrintf("time: %e\n", t);

    
    if (ns_params->height_scheme == 0) /* kinematic */
	get_dS1(ns, tstep);
    else if (ns_params->height_scheme == 1) /* non-conservative */
	get_dS2(ns, tstep);
    else if (ns_params->height_scheme == 2) /* conservative */
	get_dS3(ns, tstep);
    else
	phgError(0, "Unknown scheme.\n");

    if (rank == 0) {
	/* Add to surf */
	double e_max = -1e10, e_min = 1e10;
	double t0 = t;
	FILE *fp;
	char filename[1000];

	/* 
	 * Compute accumulation height increase using a intergral scheme.
	 * Note: da = s(t + dt) - s(t) only for u \cdot n == 0.
	 * */
	double da[Nx+1][Ny+1];
#if 1
	phgPrintf("Compute accumulaiton by interal.\n");
	for (i = 0; i < Nx+1; i++) {
	    for (j = 0; j < Ny+1; j++) {
		double x = *X(i, j, 0);
		double y = *X(i, j, 1);
		func_s(x*1e3, y*1e3, 0, &da[i][j]); /* m */
	    }
	}
	setFuncTime(t0 + dt);
	for (i = 0; i < Nx+1; i++) {
	    for (j = 0; j < Ny+1; j++) {
		double x = *X(i, j, 0);
		double y = *X(i, j, 1);
		double s;
		func_s(x*1e3, y*1e3, 0, &s); /* m */
		da[i][j] = (s - da[i][j]) / dt;
	    }
	}
	setFuncTime(t0);
#else
	phgPrintf("Compute accumulaiton at t^n.\n");
	for (i = 0; i < Nx+1; i++) {
	    for (j = 0; j < Ny+1; j++) {
		double x = *X(i, j, 0);
		double y = *X(i, j, 1);
		double a;
		func_a(x*1e3, y*1e3, &a); /* m */
		da[i][j] = a;
	    }
	}
#endif



	sprintf(filename, "./output/e_%05d.dat", tstep);
	fp = fopen(filename, "w");
	setFuncTime(t0 + dt);
	for (i = 0; i < Nx+1; i++) {
	    for (j = 0; j < Ny+1; j++) {
		double x = *X(i, j, 0);
		double y = *X(i, j, 1);
		double z = *X(i, j, 2);
		double exact_s;

		*S(i, j) += dt * (*DS(i, j) + da[i][j]) / LEN_SCALING;

		func_s(x*1e3, y*1e3, 0, &exact_s); 
		exact_s /= 1e3;	/* km */

		double e = exact_s - *S(i, j);
		if (e > e_max)
		    e_max = e;
		if (e < e_min)
		    e_min = e;

		fprintf(fp, "%e %e %e " /* 1-3 */
			"%e %e %e\n",	/* 4-5 */
			x, y, z,
			exact_s, *S(i, j), e 
			);
	    }
	    fprintf(fp, "\n");
	}
	setFuncTime(t0);
	printf("e : [%e %e]\n", e_min, e_max);
	fclose(fp);
    }
    MPI_Bcast(s_dat, (Nx+1)*(Ny+1), PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);


    /* Move mesh */
    for (p = 0; p < g->nvert; p++) {
	/* idx => (i, j, l) */
	int idx = GlobalVertex(g, p);
	assert(idx >=0 && idx < g->nvert_global);
	l = idx % (Nz+1);
	idx /= (Nz+1);
	j = idx % (Ny+1);
	idx /= (Ny+1);
	i = idx;
	assert((Nz+1) * (i * (Ny+1) + j) + l == GlobalVertex(g, p));
	idx = GlobalVertex(g, p);

	assert(i >= 0 && i < Nx+1);
	assert(j >= 0 && j < Ny+1);
	assert(l >= 0 && l < Nz+1);

	if (l == 0)
	    assert(g->types_vert[p] & BC_BOTTOM);
	if (l == Nz)
	    assert(g->types_vert[p] & BC_TOP);
        
        int idx0 = p - l;
	assert(idx0 >=0 && idx0 < g->nvert);
        double z0 = g->verts[idx0][2];
        double z1 = *S(i, j);    
        g->verts[p][2] = z0 + l * (z1 - z0) / Nz;
    }
    
    phgGeomInit_(g, TRUE);
    return;
}


