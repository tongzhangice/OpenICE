/*
 * 3D interpolation of Greenland data: Greenland_5km_dev1.1.nc.
 *
 *  */
#include "phg.h"
#include <unistd.h>
#include <string.h>
#include <math.h>


#define BC_SURFACE BDRY_USER1
#define BC_BEDROCK BDRY_USER2
#define BC_LATARAL BDRY_USER3

#define OUTPUT_DIR "./output/"
#define Bzero(x) bzero(x, sizeof(x))

static int
my_bc_map(int bctype)
{
    switch (bctype) {
    case 0:
	return BC_LATARAL;
    case 1:
	return BC_BEDROCK;
    case 2:
	return BC_SURFACE;
    default:
	return DIRICHLET;
    }
}

static void
func_surf_vel(FLOAT x, FLOAT y, FLOAT z, FLOAT *u)
{
    nc_data_set_active("surfvelx");
    interp_nc_data(x, y, z, u);
    nc_data_set_active("surfvely");
    interp_nc_data(x, y, z, u+1);
    u[2] = 0;
}

static void
func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord)
{
    FLOAT Ha, Hb, Hz;
    FLOAT rz = z;		/* ratio z in [0, 1] */

    nc_data_set_active("thk");
    interp_nc_data(x, y, z, &Hz);

    nc_data_set_active("topg");
    interp_nc_data(x, y, z, &Hb);

    coord[0] = x;
    coord[1] = y; 
    coord[2] = Hb + rz * Hz;
}



#define GET_HERE printf("Get here: %s, L %d, %s\n",	\
			__FILE__, __LINE__, __FUNCTION__);
#define NbasFace(u) (3 * (u->type->np_vert + u->type->np_edge)	\
		     + u->type->np_face)


/***************************/
/* Check all kinds of bdry */
/***************************/

void
checkBdry(GRID *g)
{
    SIMPLEX *e;
    int s;
    double a[8];

    Bzero(a);
    ForAllElements(g, e) {
        a[7] += phgGeomGetVolume(g, e);
        for (s = 0; s < NFace; s++) {
            FLOAT area = phgGeomGetFaceArea(g, e, s);
            if (e->bound_type[s] & BDRY_MASK) {
                if (e->bound_type[s] & BC_LATARAL)
                    a[0] += area;
                if (e->bound_type[s] & BC_SURFACE)
                    a[1] += area;
                if (e->bound_type[s] & BC_BEDROCK)
                    a[2] += area;
            }
            else {
                a[6] += area;
            }
        }
    }

#if USE_MPI
    {
        double b[8];
        MPI_Reduce(a, b, 8, MPI_DOUBLE, MPI_SUM, 0, g->comm);
        memcpy(a, b, sizeof(b));
    }
#endif

    phgPrintf("Boundary types check:\n");
    phgPrintf("\tLateral: %20.10e\n\tSurface: %20.10e\n\tBedrock: %20.10e\n"
              "\tsphere : %20.10e\n\tCyl    : %20.10e\n"
              "\tOther  : %20.10e\n\tAll    : %20.10e\n\tVolume : %20.10e \n",
              a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]);

    return;
}



int
main(int argc, char *argv[])
{
    GRID *g;
    SIMPLEX *e;
    DOF *u_h, *coord, *thk, *lat, *lon;
    SOLVER *solver;
    char *fn = "./cube.D6.mesh", vtk_file[1000];

    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&fn);

    phgInit(&argc, &argv);
    g = phgNewGrid(-1);
    phgImportSetBdryMapFunc(my_bc_map);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    checkBdry(g);

    read_nc_data("./greenland/Greenland_5km_dev1.1.nc", 5000);
    u_h = phgDofNew(g, DOF_DEFAULT, 3, "u_h", func_surf_vel);
    coord = phgDofNew(g, DOF_P1, 3, "coord", func_ice_slab);

    nc_data_set_active("thk");
    thk = phgDofNew(g, DOF_P2, 1, "thk", interp_nc_data);

    nc_data_set_active("lat");
    lat = phgDofNew(g, DOF_P2, 1, "lat", interp_nc_data);

    nc_data_set_active("lon");
    lon = phgDofNew(g, DOF_P2, 1, "lon", interp_nc_data);


    phgPrintf("%d DOF, %d elements, %d submeshes, load imbalance: %lg\n",
	      DofGetDataCountGlobal(u_h), g->nleaf_global, g->nprocs,
	      (double)g->lif);

#if 0    
    /* ----------------------------------------
     * 
     *   Map to physical region
     * 
     * ---------------------------------------- */
    phgExportVTK(g, OUTPUT_DIR "greenland_xyz.vtk",
		 coord, u_h, thk, NULL);
    ForAllElements(g, e) {
	FLOAT *v0, *v1;
	int i;
	for (i = 0; i < NVert; i++) {
	    v1 = DofVertexData(coord, e->verts[i]);
	    v0 = g->verts[e->verts[i]];
	    memcpy(v0, v1, Dim*sizeof(*v0));
	}
    }
    GET_HERE;
#endif


    /* ----------------------------------------
     * 
     *   Initial value
     * 
     * ---------------------------------------- */
    nc_data_set_active("surfvelx");
    //nc_data_set_active("lon");
    //nc_data_set_active("thk");
    DOF_USER_FUNC func_proj = interp_nc_data;
    //DOF_USER_FUNC func_proj = func_const;
    DOF *dof_proj = phgDofNew(g, DOF_P2, 1, "DOF proj", DofNoAction);
    solver = phgSolverCreate(SOLVER_DEFAULT, dof_proj, NULL);
    phgPrintf("  DOF: %d, unknowns: %d, Dirichlet bdry: %d\n",
	      DofGetDataCountGlobal(dof_proj), solver->rhs->map->nglobal,
	      solver->rhs->map->bdry_nglobal);

    FLOAT Area = 0;
    ForAllElements(g, e) {
	int N = dof_proj->type->nbas;
	int i, j, s, q, order = 6;
	FLOAT A[N][N], E[N][N], rhs[N];
	INT I[N];
	QUAD *quad;
	FLOAT vol;
	const FLOAT *w, *p;

	bzero(A, sizeof(A));
	bzero(E, sizeof(E));
	bzero(rhs, sizeof(rhs)); 

	for (i = 0; i < N; i++) 
	    I[i] = phgMapE2L(solver->mat->cmap, 0, e, i);

	/* Face L2 projection */
	for (s = 0; s < NFace; s++) {
	    if (e->bound_type[s] & BC_SURFACE) {
		int ii, jj;
		int v0, v1, v2;
		int nbas_face = NbasFace(dof_proj);
		SHORT bases[nbas_face];
		FLOAT lambda[Dim + 1], area, x, y, z, beta;
		FLOAT vf;
		order = DofTypeOrder(dof_proj, e) * 3 - 1;

		phgDofGetBasesOnFace(dof_proj, e, s, bases);
		v0 = GetFaceVertex(s, 0);
		v1 = GetFaceVertex(s, 1);
		v2 = GetFaceVertex(s, 2);
		lambda[s] = 0.;

		/* printf("Elem %d face: %d\n", e->index, s); */
		/* for (i = 0; i < nbas_face; i++) { */
		/*     ii = bases[i]; */
		/*     printf("%d ", ii); */
		/* } */
		/* printf("\n"); */

		area = phgGeomGetFaceArea(g, e, s);
		quad = phgQuadGetQuad2D(order);
		Area += area;

		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    lambda[v0] = *(p++);
		    lambda[v1] = *(p++);
		    lambda[v2] = *(p++);
			
		    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		    func_proj(x, y, z, &vf);

		    for (i = 0; i < nbas_face; i++) {
			ii = bases[i];
			FLOAT gi_u = 
			    *dof_proj->type->BasFuncs(dof_proj, e, ii, ii + 1, lambda);
			for (j = 0; j < nbas_face; j++) { 
			    jj = bases[j];
			    FLOAT gj_u = 
				*dof_proj->type->BasFuncs(dof_proj, e, jj, jj + 1, lambda);
			    FLOAT mass_face = area*(*w) * (gj_u)*(gi_u);

			    A[ii][jj] += mass_face;
			} /* end of bas_j */
			rhs[ii] += area*(*w) * vf * (gi_u);
		    }     /* end of bas_i */
		    w++;
		}		/* end of quad point */
	    }			/* end of face outflow */
	}			/* end of all outflow face in element */
	    
	for (i = 0; i < N; i++)
	    E[i][i] = 1.;

	for (i = 0; i < N; i++) {
	    int gind, bind;
	    GTYPE gtype;
	    BTYPE btype = phgDofGetElementBasisInfo(dof_proj, e, i,
						    &gtype, &gind, &bind);
	    if (btype & BC_SURFACE) {
		phgMatAddEntries(solver->mat, 1, I + i, N, I, A[i]); 
		//assert(fabs(A[i][i]) > 1e-8);
	    } else {	/* interior node */
		phgMatAddEntries(solver->mat, 1, I + i, N, I, E[i]); 
	    }
	}

	phgVecAddEntries(solver->rhs, 0, N, I, rhs);
    }

    phgPrintf("Area: %e\n", Area);

    phgMatAssemble(solver->mat);
    phgVecAssemble(solver->rhs);

    //phgMatDumpMATLAB(solver->mat, "A", "A_.m");
    //phgVecDumpMATLAB(solver->rhs, "b", "b_.m");

    phgSolverUpdateRHS(solver);
    phgPrintf(0, "\nBuild solver\n");
    phgSolverSolve(solver, TRUE, dof_proj, NULL);
    phgPrintf("  nits = %d, ", solver->nits);










    /* export distord mesh */
    phgExportMedit(g, OUTPUT_DIR "greenland.mesh"); 
    phgExportVTK(g, OUTPUT_DIR "greenland.vtk",
    		 coord, u_h, thk, dof_proj, lat, lon, NULL);
    phgExportEnsight(g, OUTPUT_DIR "greenland",
    		     coord, u_h, thk, dof_proj, lat, lon, NULL);
    phgDofFree(&u_h);
    phgDofFree(&coord);

    phgPrintf("Exit normally.\n");
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
