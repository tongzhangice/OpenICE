
#include "ins.h"


void
get_mask_bot(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    //DOF *mask_bot = ns->mask_bot;

    //ns->mask_bot = phgDofNew(g, DOF_P1, 1, "mask bot", DofNoAction);

    INT i, j, s;

    ForAllElements(g, e)
    {
        for (s = 0; s < NFace; s++)
        {
            if (e->bound_type[s] & BC_ISHELF)
            {
                for (i = 0; i < 3; i++)
                {
                    INT v = GetFaceVertex(s, i);
                    INT local_idx = e->verts[v];

                    FLOAT *mask = DofVertexData(ns->mask_bot, local_idx);
                    mask[0] = 1;
                }
            }
            if (e->bound_type[s] == (BC_BOTTOM|BC_BOTTOM_GRD))//((e->bound_type[s] & BC_BOTTOM_GRD) && (!(e->bound_type[s] & BC_ISHELF)))
            {
                for (i = 0; i < 3; i++)
                {
                    INT v = GetFaceVertex(s, i);
                    INT local_idx = e->verts[v];

                    FLOAT *mask = DofVertexData(ns->mask_bot, local_idx);
                    mask[0] = -1;
                }
            }
        }
    }

    phgExportVTK(g, "mask.vtk", ns->mask_bot, NULL);
    //phgDofFree(&mask_bot);
}

void
modify_mask_bot(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;

    INT i, j, s;

    ForAllElements(g, e)
    {
        for (s = 0; s < NFace; s++)
        {
            if (e->bound_type[s] & BC_ISHELF)
            {
                for (i = 0; i < 3; i++)
                {
                    INT v = GetFaceVertex(s, i);
                    INT local_idx = e->verts[v];

                    if (g->types_vert[local_idx] & BC_BOTTOM_GRD)
                    {
                        e->bound_type[s] = (BC_BOTTOM_GRD|BC_BOTTOM);
                    }
                }
            }
        }
    }

    phgUpdateBoundaryTypes(g);
}

void update_grounding_line(NSSolver *ns, int tstep)
{

    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    SIMPLEX *e;

    INT s, i;

    INT if_change = 0, if_change0;

    //if (phgRank == 0)
    //    system("rm ../../mesh_gen/polyfile/gl_points.txt");

    phgPrintf("reset bound type start!\n");
    ForAllElements(g, e)
    {
        for (s = 0; s < NFace; s++)
        {
            if ((e->bound_type[s] & BC_BOTTOM))
            {
                
                QUAD *quad = phgQuadGetQuad2D(1);
                FLOAT *p = quad->points;
                FLOAT lambda[Dim + 1], x,y,z;
                int v0, v1, v2;
                FLOAT bed_z0, bed_z1, bed_z2;
                v0 = GetFaceVertex(s, 0);
                v1 = GetFaceVertex(s, 1);
                v2 = GetFaceVertex(s, 2);

                INT local_idx0 = e->verts[v0];
                INT local_idx1 = e->verts[v1];
                INT local_idx2 = e->verts[v2];

                FLOAT x0 = g->verts[local_idx0][0];
                FLOAT x1 = g->verts[local_idx1][0];
                FLOAT x2 = g->verts[local_idx2][0];

                FLOAT y0 = g->verts[local_idx0][1];
                FLOAT y1 = g->verts[local_idx1][1];
                FLOAT y2 = g->verts[local_idx2][1];

                FLOAT z0 = g->verts[local_idx0][2];
                FLOAT z1 = g->verts[local_idx1][2];
                FLOAT z2 = g->verts[local_idx2][2];

                lambda[s] = 0.;
                lambda[v0] = *(p++);
                lambda[v1] = *(p++);
                lambda[v2] = *(p++);
                phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);

                func_bed_z(x0, y0, z0, &bed_z0);
                func_bed_z(x1, y1, z1, &bed_z1);
                func_bed_z(x2, y2, z2, &bed_z2);

                if ((z0 > bed_z0 + 1e-3) || (z1 > bed_z1 + 1e-3) || (z2 > bed_z2 + 1e-3))
                {
                    e->bound_type[s] = (BC_ISHELF | BC_BOTTOM);
                    if_change = 1;
                }
                //else if (e->bound_type[s] & BC_ISHELF 
                //        && !(e->bound_type[s] & BC_BOTTOM_GRD))
                else
                {
                        // when originally floating and now becomes attched to bed
                        e->bound_type[s] = (BC_BOTTOM | BC_BOTTOM_GRD);
                }
                //else
                //    ; // do nothing             
            }
        }
    }

    phgUpdateBoundaryTypes(g);

    MPI_Allreduce(&if_change, &if_change0, 1, PHG_MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if (if_change0 == 1)
    {
        phgPrintf("\n---------------------------------------\n");
        phgPrintf("some grounded elements are now afloat !");
        phgPrintf("\n---------------------------------------\n");
    }
    else
    {
        phgPrintf("\n======================\n");
        phgPrintf("no grounded element becomes afloat !");
        phgPrintf("\n======================\n");
    }

    phgPrintf("update bound type done!\n");

	return;
}


#if 1
void get_water_pressure(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    INT s, i, j, k, l, q;


    phgPrintf("\n Compute water pressure loads !!\n");

    DOF *wp = phgDofNew(g, DOF_P2, Dim, "water pressure", DofNoAction);
    MAP *map_wp = phgMapCreate(wp, NULL);
    VEC *rhs_wp = phgMapCreateVec(map_wp, Dim);
    phgVecDisassemble(rhs_wp);

    INT M = wp->type->nbas;
    BTYPE btype[M*Dim];
    ForAllElements(g, e)
    {
        FLOAT rhs_wp_e[M][Dim];
        INT local_map_idx[M][Dim];

        memset(rhs_wp_e, 0, sizeof rhs_wp_e);

        for (s = 0; s < NFace; s++)
        {

            if (e->bound_type[s] & BC_BOTTOM)
            {

                FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
                INT M_face = 3*(wp->type->np_vert + wp->type->np_edge);
                SHORT bas_idx_e[M_face];
                INT quad_order = 5;
                phgDofGetBasesOnFace(wp, e, s, bas_idx_e);
                
                const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0;

                QUAD *quad = phgQuadGetQuad2D(quad_order);
                FLOAT *p = quad->points;
                FLOAT *w = quad->weights;
                FLOAT area = phgGeomGetFaceArea(g, e, s);

                for (q = 0; q < quad->npoints; q++)
                {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);

                    phgGeomLambda2XYZ(g, e, lambda, &lambda_x, &lambda_y, &lambda_z);
                    //const FLOAT *bas = wp->type->BasFuncs(wp, e, 0, -1, lambda);

                    for (i = 0; i < M_face; i++)
                    {
                        INT i_e = bas_idx_e[i];
                    const FLOAT bas = *wp->type->BasFuncs(wp, e, i_e, i_e+1, lambda);

                        for (k = 0; k < Dim; k++)
                        {
                            if (lambda_z < 0)
                            {
                                rhs_wp_e[i_e][k] += area*w[q]*RHO_WAT*GRAVITY*(-lambda_z)*normal[k]*bas/1.0e5;
                            }
                            else
                                rhs_wp_e[i_e][k] += 0;
                        }
                    }
                }

            }
        }

        for (i = 0; i < M; i++)
        for (k = 0; k < Dim; k++)
        {
            local_map_idx[i][k] = phgMapE2L(map_wp, 0, e, i*Dim+k);
        }
        phgVecAddEntries(rhs_wp, 0, M*Dim, local_map_idx[0], &rhs_wp_e[0][0]);
    }

    phgVecAssemble(rhs_wp);
    phgMapLocalDataToDof(map_wp, 1, &wp, rhs_wp->data);
    phgExportVTK(g, "wp.vtk", wp, NULL);


}
#endif

#if 0
void get_water_pressure(GRID *g, DOF *wp)
{

    SIMPLEX *e;
    INT s, i, j, k, l, q;



    DOF *wp1 = phgDofNew(g, DOF_P2, Dim, "wp", DofInterpolation);
    MAP *Vmap = phgMapCreate(wp, NULL);
    VEC *rhs_wp = phgMapCreateVec(Vmap, Dim);
    phgVecDisassemble(rhs_wp);


    INT M = wp->type->nbas;
    BTYPE btype[M*Dim];

    ForAllElements(g, e)
    {
        FLOAT rhs_wp_e[M][Dim];
        FLOAT lhs_wp_e[M][Dim][M][Dim];

        INT local_map_idx[M][Dim];

        memset(rhs_wp_e, 0, sizeof rhs_wp_e);
        memset(lhs_wp_e, 0, sizeof lhs_wp_e);

        for (s = 0; s < NFace; s++)
        {
            if (e->bound_type[s] & BC_BOTTOM_GRD )
            {
                FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
                INT quad_order = 5;
                QUAD *quad = phgQuadGetQuad2D(quad_order);
                FLOAT *p = quad->points;
                FLOAT *w = quad->weights;
                FLOAT area = phgGeomGetFaceArea(g, e, s);

                const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                INT M_face = 3*(wp->type->np_vert + wp->type->np_edge);
                SHORT bas_idx_e[M_face];

                phgDofGetBasesOnFace(wp, e, s, bas_idx_e);
                
                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0;

                for (q = 0; q < quad->npoints; q++)
                {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);

                    phgGeomLambda2XYZ(g, e, lambda, &lambda_x, &lambda_y, &lambda_z);
                    const FLOAT *bas = wp->type->BasFuncs(wp, e, 0, -1, lambda);

                    for (i = 0; i < M_face; i++)
                    {
                        INT i_e = bas_idx_e[i];

                        for (k = 0; k < Dim; k++)
                        {
                            //if (lambda_z < 0)
                            {
                                rhs_wp_e[i_e][k] += area*w[q]*1000*9.8*(-lambda_z)*normal[k]*bas[i_e];
                            }
                            //else
                            //    rhs_wp_e[i_e][k] += 0;
                        }
                    }
                }

            }
        }


        for (i = 0; i < M; i++)
        for (k = 0; k < Dim; k++)
        {
            local_map_idx[i][k] = phgMapE2L(Vmap, 0, e, i*Dim+k);
        }

            phgVecAddEntries(rhs_wp, 0, M*Dim, local_map_idx[0], &rhs_wp_e[0][0]);


    }

    phgVecAssemble(rhs_wp);


    phgMapLocalDataToDof(Vmap, 1, &wp, rhs_wp->data);

    phgExportVTK(g, "wp.vtk", wp, NULL);
}
#endif

#if 1
void get_water_force(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    INT s, i, j, k, l, q;


    phgPrintf("\n Compute water pressure loads !!\n");

    DOF *wp = phgDofNew(g, DOF_P2, 1, "water pressure", DofNoAction);
    MAP *map_wp = phgMapCreate(wp, NULL);
    VEC *rhs_wp = phgMapCreateVec(map_wp, 1);
    phgVecDisassemble(rhs_wp);

    INT M = wp->type->nbas;
    BTYPE btype[M];
    SOLVER *solver = phgSolverCreate(SOLVER_DEFAULT, wp, NULL);

    ForAllElements(g, e)
    {
        FLOAT A[M][M], rhs_wp_e[M], buffer[M];
        INT local_map_idx[M];

        memset(rhs_wp_e, 0, sizeof rhs_wp_e);
        memset(A, 0, sizeof A);

        for (s = 0; s < NFace; s++)
        {

            if (e->bound_type[s] & BC_BOTTOM)
            {

                FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
                INT M_face = 3*(wp->type->np_vert + wp->type->np_edge);
                SHORT bas_idx_e[M_face];
                INT quad_order = 5;
                phgDofGetBasesOnFace(wp, e, s, bas_idx_e);
                
                const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0;

                QUAD *quad = phgQuadGetQuad2D(quad_order);
                FLOAT *p = quad->points;
                FLOAT *w = quad->weights;
                FLOAT area = phgGeomGetFaceArea(g, e, s);

                for (q = 0; q < quad->npoints; q++)
                {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);

                    phgGeomLambda2XYZ(g, e, lambda, &lambda_x, &lambda_y, &lambda_z);
                    //const FLOAT *bas = wp->type->BasFuncs(wp, e, 0, -1, lambda);

                    for (i = 0; i < M_face; i++)
                    {
                        INT i_e = bas_idx_e[i];
                    const FLOAT bas_i = *wp->type->BasFuncs(wp, e, i_e, i_e+1, lambda);
                    for (j = 0; j < M_face; j++)
                    {
                        INT j_e = bas_idx_e[j];
                    const FLOAT bas_j = *wp->type->BasFuncs(wp, e, j_e, j_e+1, lambda);
                        A[i_e][j_e] += area*w[q]*bas_i*bas_j;
                    }

                        rhs_wp_e[i_e] += area*w[q]*RHO_WAT*GRAVITY*(-lambda_z)*bas_i/1.0e5;
                        //rhs_wp_e[i_e] += lambda_z*area*w[q]*bas;
                        //phgPrintf("%e\n", w[q]*bas);
                    }
                }

            }
        }

        for (i = 0; i < M; i++)
        {
            btype[i] = phgDofGetElementBoundaryType(wp, e, i);
            if (!(btype[i] & BC_BOTTOM)) 
            {
                A[i][i] = 1.;
                rhs_wp_e[i] = 0.;
            } 
        }

        for (i = 0; i < M; i++)
        {
            local_map_idx[i] = phgMapE2L(solver->rhs->map, 0, e, i);
        }

        for (i = 0; i < M; i++) 
        {
                /*
            if (0 && phgDofDirichletBC_(wp, e, i, NULL, buffer, rhs_wp_e+i, DOF_PROJ_NONE)) 
            {
                phgPrintf("no bc !\n");
                phgSolverAddMatrixEntries(solver, 1, local_map_idx + i, M, local_map_idx, buffer); 
            }
            else 
            */
            {
                phgSolverAddMatrixEntries(solver, 1, local_map_idx+i, M, local_map_idx, A[i]); 
            }
        }

        phgSolverAddRHSEntries(solver, M, local_map_idx, &rhs_wp_e[0]);
        //phgVecAddEntries(rhs_wp, 0, M, local_map_idx, &rhs_wp_e[0]);
    }

    //solver->rhs = vec_nf0;

    phgSolverSolve(solver, TRUE, wp, NULL);
    phgPrintf("Done solving water pressure!!\n");
    DOF *wp_P1 = phgDofCopy(wp, NULL, DOF_P1, NULL);
    phgExportVTK(g, "wp_P1.vtk", wp_P1, NULL);
    phgDofFree(&wp_P1);

	phgDofCopy(wp, &ns->water_force, NULL, "water_force");
    phgDofFree(&wp);
    phgVecDestroy(&rhs_wp);
    phgMapDestroy(&map_wp);


}
#endif

#if 0
void get_water_force(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    INT s, i, j, k, l, q;


    phgPrintf("\n Compute water pressure loads !!\n");

    DOF *wp = phgDofNew(g, DOF_P2, 1, "water pressure", DofNoAction);
    MAP *map_wp = phgMapCreate(wp, NULL);
    VEC *rhs_wp = phgMapCreateVec(map_wp, Dim);
    phgVecDisassemble(rhs_wp);

    INT M = wp->type->nbas;
    BTYPE btype[M*Dim];
    ForAllElements(g, e)
    {
        FLOAT rhs_wp_e[M];
        INT local_map_idx[M];

        memset(rhs_wp_e, 0, sizeof rhs_wp_e);

        for (s = 0; s < NFace; s++)
        {

            if (e->bound_type[s] & BC_BOTTOM)
            {

                FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
                INT M_face = 3*(wp->type->np_vert + wp->type->np_edge);
                SHORT bas_idx_e[M_face];
                INT quad_order = 5;
                phgDofGetBasesOnFace(wp, e, s, bas_idx_e);
                
                const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0;

                QUAD *quad = phgQuadGetQuad2D(quad_order);
                FLOAT *p = quad->points;
                FLOAT *w = quad->weights;
                FLOAT area = phgGeomGetFaceArea(g, e, s);

                for (q = 0; q < quad->npoints; q++)
                {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);

                    phgGeomLambda2XYZ(g, e, lambda, &lambda_x, &lambda_y, &lambda_z);
                    //const FLOAT *bas = wp->type->BasFuncs(wp, e, 0, -1, lambda);

                    for (i = 0; i < M_face; i++)
                    {
                        INT i_e = bas_idx_e[i];
                    const FLOAT bas = *wp->type->BasFuncs(wp, e, i_e, i_e+1, lambda);

                        rhs_wp_e[i_e] += area*w[q]*RHO_WAT*GRAVITY*(-lambda_z)*bas/1.0e5;
                        //rhs_wp_e[i_e] += lambda_z*area*w[q]*bas;
                        //phgPrintf("%e\n", w[q]*bas);
                    }
                }

            }
        }

        for (i = 0; i < M; i++)
        {
            local_map_idx[i] = phgMapE2L(map_wp, 0, e, i);
        }
        phgVecAddEntries(rhs_wp, 0, M, local_map_idx, &rhs_wp_e[0]);
    }

    phgVecAssemble(rhs_wp);
    phgMapLocalDataToDof(map_wp, 1, &wp, rhs_wp->data);
    phgExportVTK(g, "wp.vtk", wp, NULL);

	phgDofCopy(wp, &ns->water_force, DOF_P1, "water_force");
    phgDofFree(&wp);
    phgVecDestroy(&rhs_wp);
    phgMapDestroy(&map_wp);


}
#endif

#if 0
void get_water_pressure(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    INT s, i, j, k, l, q;


    phgPrintf("\n Compute water pressure loads !!\n");

    DOF *wp = phgDofNew(g, DOF_P2, 1, "water pressure", DofNoAction);
    MAP *map_wp = phgMapCreate(wp, NULL);
    VEC *rhs_wp = phgMapCreateVec(map_wp, Dim);
    phgVecDisassemble(rhs_wp);

    INT M = wp->type->nbas;
    BTYPE btype[M*Dim];
    ForAllElements(g, e)
    {
        FLOAT rhs_wp_e[M];
        INT local_map_idx[M];

        memset(rhs_wp_e, 0, sizeof rhs_wp_e);

        for (s = 0; s < NFace; s++)
        {

            if (e->bound_type[s] & BC_BOTTOM)
            {

                FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
                INT M_face = 3*(wp->type->np_vert + wp->type->np_edge);
                SHORT bas_idx_e[M_face];
                INT quad_order = 5;
                phgDofGetBasesOnFace(wp, e, s, bas_idx_e);
                
                const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0;

                QUAD *quad = phgQuadGetQuad2D(quad_order);
                FLOAT *p = quad->points;
                FLOAT *w = quad->weights;
                FLOAT area = phgGeomGetFaceArea(g, e, s);

                for (q = 0; q < quad->npoints; q++)
                {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);

                    phgGeomLambda2XYZ(g, e, lambda, &lambda_x, &lambda_y, &lambda_z);
                    //const FLOAT *bas = wp->type->BasFuncs(wp, e, 0, -1, lambda);

                    for (i = 0; i < M_face; i++)
                    {
                        INT i_e = bas_idx_e[i];
                    const FLOAT bas = *wp->type->BasFuncs(wp, e, i_e, i_e+1, lambda);

                        rhs_wp_e[i_e] += area*w[q]*RHO_WAT*GRAVITY*(-lambda_z)*bas/1.0e5;
                        //rhs_wp_e[i_e] += lambda_z*area*w[q]*bas;
                        //phgPrintf("%e\n", w[q]*bas);
                    }
                }

            }
        }

        for (i = 0; i < M; i++)
        {
            local_map_idx[i] = phgMapE2L(map_wp, 0, e, i);
        }
        phgVecAddEntries(rhs_wp, 0, M, local_map_idx, &rhs_wp_e[0]);
    }

    phgVecAssemble(rhs_wp);
    phgMapLocalDataToDof(map_wp, 1, &wp, rhs_wp->data);
    phgPrintf("Save water presure!!\n");
    phgExportVTK(g, "wp.vtk", wp, NULL);


}
#endif 


/*
void get_nodal_force(NSSolver *ns, INT IF_DB)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    INT s, i, k, q;

    MAP *Vmap, *Pmap;
    MAT *pmat[9];

    DOF *res_u = phgDofNew(ns->g, DOF_P2, Dim, "res_u", DofNoAction);
    DOF *res_p = phgDofNew(ns->g, DOF_P1, 1, "res_p", DofNoAction);


    phgNSDestroySolverU(ns);

#if 1
    ns->Vmap = phgMapCreate(ns->du, NULL);
    ns->Pmap = phgMapCreate(ns->dp, NULL);

    Vmap = ns->Vmap;
    Pmap = ns->Pmap;

    ns->matF = phgMapCreateMat(Vmap, Vmap);
    ns->matBt = phgMapCreateMat(Vmap, Pmap);
    ns->matB = phgMapCreateMat(Pmap, Vmap);
    ns->matC = phgMapCreateMat(Pmap, Pmap);

    ns->matF->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;
    if (ns->matC != NULL)
	ns->matC->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;

    ns->matB->handle_bdry_eqns = FALSE;
    ns->matBt->handle_bdry_eqns = FALSE;

    pmat[0] = ns->matF;
    pmat[1] = ns->matBt;
    pmat[2] = ns->matB;
    pmat[3] = ns->matC;

    ns->matNS = phgMatCreateBlockMatrix(g, 2, 2, pmat, NULL, NULL);
    ns->matNS->mv_data = phgAlloc(sizeof(*ns->matNS->mv_data));
    ns->matNS->mv_data[0] = (void *) ns;
    ns->matNS->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;

    phgOptionsPush();
    //phgOptionsSetOptions(_nsp->Stokes_opts);
    ns->solver_u = phgSolverCreate(SOLVER_DEFAULT, ns->du, ns->dp, NULL);
    phgOptionsPop();

    phgMatDestroy(&ns->solver_u->mat);
    ns->solver_u->mat = ns->matNS;
    ns->solver_u->rhs->mat = ns->solver_u->mat;
    //phgDofSetDataByValue(ns->du, 0.);
    //phgDofSetDataByValue(ns->dp, 0.);
        ns->ltype = PICARD;
#endif
            //phgNSDestroySolverU(ns);
    //phgNSInitSolverU(ns);
    phgNSBuildSolverURHS(ns, 0, 0);
    //phgNSBuildSolverURHS(ns);
    phgVecDumpMATLAB(ns->solver_u->rhs, "rhs1", "rhs1_.m");
    phgNSBuildSolverUMat(ns, 0, 0);
    //phgNSBuildSolverUMat(ns);
            phgMatDumpMATLAB(ns->matB, "mat1", "mat1_.m");

    MAT *all_mat[4];

    all_mat[0] = ns->matF;
    all_mat[1] = ns->matBt;
    all_mat[2] = ns->matB;
    all_mat[3] = ns->matC;

    MAT *u_mat = phgMatCreateBlockMatrix(g, 2, 2, all_mat, NULL, NULL);

    MAT *lhs = ns->matNS;
    VEC *rhs = ns->solver_u->rhs;


    VEC *res = phgVecCopy(rhs, NULL);

    DOF *sols[2] = {ns->du, ns->dp};
    //DOF *sols[2] = {velocity, pressure};

    DOF *res_u_p[2] = {res_u, res_p};

    VEC *vec_sol = phgMapCreateVec(rhs->map, 1);
        
    phgMapDofToLocalData(rhs->map, 2, &sols[0], vec_sol->data);


    //VEC *vec_sol1 = phgMapCreateVec(rhs->map, 1);

#if 0
    VEC *vec_sol1 = phgMatVec(1, 1, u_mat, rhs, 0, NULL);
    phgVecDumpMATLAB(vec_sol, "vec_sol", "vec_sol.m");
    phgVecDumpMATLAB(vec_sol1, "vec_sol1", "vec_sol1.m");

#endif
    phgMatVec(0, 1, u_mat, vec_sol, -1, &res);
    
    phgMapLocalDataToDof(rhs->map, 2, &res_u_p[0], res->data);
    phgExportVTK(ns->g, "res.vtk", res_u, res_p, NULL); 


}
*/

void get_nodal_force(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    INT s, i, j, k, l, q;
    //MAP *Vmap, *Pmap;
    //MAT *pmat[9];

    DOF *res_u = phgDofNew(ns->g, DOF_P2, Dim, "res_u", DofNoAction);
    DOF *res_p = phgDofNew(ns->g, DOF_P1, 1, "res_p", DofNoAction);

    MAT *lhs = ns->matNS0;
    VEC *rhs = ns->vec_rhs0;


    VEC *res = phgVecCopy(rhs, NULL);

    DOF *sols[2] = {ns->du, ns->dp};
    //DOF *sols[2] = {velocity, pressure};

    DOF *res_u_p[2] = {res_u, res_p};

    VEC *vec_sol = phgMapCreateVec(rhs->map, 1);
        
    phgMapDofToLocalData(rhs->map, 2, &sols[0], vec_sol->data);

    phgPrintf("\n Compute nodal loads !!\n");
    phgMatVec(0, 1, lhs, vec_sol, -1, &res);
    
    phgMapLocalDataToDof(rhs->map, 2, &res_u_p[0], res->data);
    phgExportVTK(ns->g, "res.vtk", res_u, res_p, NULL); 


    DOF *nf = phgDofCopy(res_u, NULL, NULL, NULL);
    DOF *nf0 = get_dof_component(g, nf, DOF_P2, 3, 0);

    phgPrintf("\n Compute nodal 'pressure'!!\n\n");

    DOF *wp = phgDofNew(g, DOF_P2, 1, "nodal pressure", DofNoAction);
    MAP *map_wp = phgMapCreate(wp, NULL);
    VEC *rhs_wp = phgMapCreateVec(map_wp, 1);
    phgVecDisassemble(rhs_wp);

    VEC *vec_nf0 = phgMapCreateVec(map_wp, 1);
        
    phgMapDofToLocalData(map_wp, 1, &nf0, vec_nf0->data);

    INT M = wp->type->nbas;
    BTYPE btype[M];
    SOLVER *solver = phgSolverCreate(SOLVER_DEFAULT, wp, NULL);

    ForAllElements(g, e)
    {
        FLOAT A[M][M], rhs_wp_e[M], buffer[M];
        INT local_map_idx[M];

        memset(rhs_wp_e, 0, sizeof rhs_wp_e);
        memset(A, 0, sizeof A);

        for (s = 0; s < NFace; s++)
        {

            if (e->bound_type[s] & BC_BOTTOM)
            {

                FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
                INT M_face = 3*(wp->type->np_vert + wp->type->np_edge);
                SHORT bas_idx_e[M_face];
                INT quad_order = 5;
                phgDofGetBasesOnFace(wp, e, s, bas_idx_e);
                
                const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0;

                QUAD *quad = phgQuadGetQuad2D(quad_order);
                FLOAT *p = quad->points;
                FLOAT *w = quad->weights;
                FLOAT area = phgGeomGetFaceArea(g, e, s);

                for (q = 0; q < quad->npoints; q++)
                {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);

                    phgGeomLambda2XYZ(g, e, lambda, &lambda_x, &lambda_y, &lambda_z);
                    //const FLOAT *bas = wp->type->BasFuncs(wp, e, 0, -1, lambda);

                    for (i = 0; i < M_face; i++)
                    {
                        INT i_e = bas_idx_e[i];
                    const FLOAT bas_i = *wp->type->BasFuncs(wp, e, i_e, i_e+1, lambda);
                    for (j = 0; j < M_face; j++)
                    {
                        INT j_e = bas_idx_e[j];
                    const FLOAT bas_j = *wp->type->BasFuncs(wp, e, j_e, j_e+1, lambda);
                        A[i_e][j_e] += area*w[q]*bas_i*bas_j;
                    }

                        //rhs_wp_e[i_e] += area*w[q]*RHO_WAT*GRAVITY*(-lambda_z)*bas_i/1.0e5;
                        //rhs_wp_e[i_e] += lambda_z*area*w[q]*bas;
                        //phgPrintf("%e\n", w[q]*bas);
                    }
                }

            }
        }

        for (i = 0; i < M; i++)
        {
            btype[i] = phgDofGetElementBoundaryType(wp, e, i);
            if (!(btype[i] & BC_BOTTOM)) 
            {
                A[i][i] = 1.;
                rhs_wp_e[i] = 0.;
            } 
        }

        for (i = 0; i < M; i++)
        {
            local_map_idx[i] = phgMapE2L(solver->rhs->map, 0, e, i);
        }

        for (i = 0; i < M; i++) 
        {
                /*
            if (0 && phgDofDirichletBC_(wp, e, i, NULL, buffer, rhs_wp_e+i, DOF_PROJ_NONE)) 
            {
                phgPrintf("no bc !\n");
                phgSolverAddMatrixEntries(solver, 1, local_map_idx + i, M, local_map_idx, buffer); 
            }
            else 
            */
            {
                phgSolverAddMatrixEntries(solver, 1, local_map_idx+i, M, local_map_idx, A[i]); 
            }
        }

        //phgSolverAddRHSEntries(solver, M, local_map_idx, &rhs_wp_e[0]);
        //phgVecAddEntries(rhs_wp, 0, M, local_map_idx, &rhs_wp_e[0]);
    }

    phgVecAXPBY(1, vec_nf0, 0, &solver->rhs);
    //solver->rhs = vec_nf0;

    phgSolverSolve(solver, TRUE, wp, NULL);
    //phgPrintf("\n Done solving nodal 'pressure'!!\n\n");
    DOF *wp_P1 = phgDofCopy(wp, NULL, DOF_P1, NULL);
    phgExportVTK(g, "nodal_P1.vtk", wp_P1, NULL);
    phgDofFree(&wp_P1);

	phgDofCopy(wp, &ns->nodal_force, NULL, "nodal_force");
	//phgDofCopy(res_u, &ns->nodal_force, DOF_P1, "nodal_force");
    phgDofFree(&wp);
    phgVecDestroy(&rhs_wp);
    phgMapDestroy(&map_wp);

	phgDofFree(&nf);
	phgDofFree(&nf0);

	phgDofFree(&res_u);
	phgDofFree(&res_p);
    phgVecDestroy(&vec_sol);
    phgVecDestroy(&res);
    phgVecDestroy(&rhs);
    phgMatDestroy(&lhs);


}


void get_nodal_force_value(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    INT s, i, j, k, q;

    MAP *Vmap, *Pmap;
    MAT *pmat[9];

    DOF *res_u = phgDofNew(ns->g, DOF_P2, Dim, "res_u", DofNoAction);
    DOF *res_p = phgDofNew(ns->g, DOF_P1, 1, "res_p", DofNoAction);

    MAT *lhs = ns->matNS0;
    VEC *rhs = ns->vec_rhs0;


    VEC *res = phgVecCopy(rhs, NULL);

    DOF *sols[2] = {ns->du, ns->dp};
    //DOF *sols[2] = {velocity, pressure};

    DOF *res_u_p[2] = {res_u, res_p};

    VEC *vec_sol = phgMapCreateVec(rhs->map, 1);
        
    phgMapDofToLocalData(rhs->map, 2, &sols[0], vec_sol->data);

    phgMatVec(0, 1, lhs, vec_sol, -1, &res);
    
    phgMapLocalDataToDof(rhs->map, 2, &res_u_p[0], res->data);
    //phgExportVTK(ns->g, "res.vtk", res_u, res_p, NULL); 
    //


    DOF *res_u1 = phgDofCopy(res_u, NULL, DOF_P1, NULL);
    MAP *map_res_u = phgMapCreate(res_u1, NULL);
    VEC *vec_res_u = phgMapCreateVec(map_res_u, 1);

    phgMapDofToLocalData(map_res_u, 1, &res_u1, vec_res_u->data);

    SOLVER *solver;

    DOF *nodal_force = phgDofNew(g, DOF_P1, Dim, "nodal_force", DofNoAction); 

    MAT *mat = phgMapCreateMat(map_res_u, map_res_u);

    //phgOptionsPush();
    //phgOptionsSetOptions("-solver gmres "
	//		 "-solver_maxit 100 "
	//		 "-solver_rtol 1e-10");
    solver = phgSolverCreate(SOLVER_DEFAULT, nodal_force, NULL);
    //phgOptionsPop();

    //printf("start assembly the matrix !!\n");
    ForAllElements(g, e)
    {

        int M = nodal_force->type->nbas;
        FLOAT F[M][Dim][M][Dim];
        INT Iu[M][Dim];
        const FLOAT *w, *p;
        QUAD *quad;
        BTYPE btype[M][Dim];
        Bzero(F); 
#if 0
        //int M = nodal_force->type->nbas;
        //FLOAT F[M][Dim][M][Dim], rhs[M][Dim];
        //INT Iu[M][Dim];
        //const FLOAT *w, *p;
        FLOAT vol;

        vol = phgGeomGetVolume(g, e);
        QUAD *quad1 = phgQuadGetQuad3D(5);

        Bzero(F); 
        //Bzero(rhs);

        p = quad1->points;
        w = quad1->weights;
        for (q = 0; q < quad1->npoints; q++)
        {
            for (i = 0; i < M; i++)
            {
                const FLOAT *bas_i = phgQuadGetBasisValues(e, nodal_force, i, quad1);

                for (j = 0; j < M; j++)
                {
                    const FLOAT *bas_j = phgQuadGetBasisValues(e, nodal_force, j, quad1);

                    for (k = 0; k < Dim; k++)
                    {
                          F[j][k][i][k] += vol*w[q]*bas_i[q]*bas_j[q]*EQU_SCALING; 
                    }
                }

                //for (k = 0; k < Dim; k++)
                //    rhs[i][k] += vol*w[q]*bas_i[q]*EQU_SCALING;
            }
        }
#endif

        for (s = 0; s < NFace; s++)
        {

	    int v0, v1, v2, ii, jj;
	    int nbas_face = NbasFace(nodal_force);
	    SHORT bases[nbas_face];
	    FLOAT x, y, z, lambda[Dim + 1], area, gi, gj;
	    const FLOAT *normal;

        if (e->bound_type[s] & BC_BOTTOM_GRD)
        {

	    phgDofGetBasesOnFace(nodal_force, e, s, bases);
	    v0 = GetFaceVertex(s, 0);
	    v1 = GetFaceVertex(s, 1);
	    v2 = GetFaceVertex(s, 2);
	    lambda[s] = 0.;

	    area = phgGeomGetFaceArea(g, e, s);
	    const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);
	    quad = phgQuadGetQuad2D(5);

	    /* Projection to 2D plan */
	    //area *= normal[Z_DIR];

	    p = quad->points;
	    w = quad->weights;

	    for (q = 0; q < quad->npoints; q++) 
        {

		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);

        const FLOAT *bas_v = nodal_force->type->BasFuncs(nodal_force, e, 0, -1, lambda);

        for (ii = 0; ii < nbas_face; ii++)
        {
            i = bases[ii];
            gi = bas_v[i];

        for (jj = 0; jj < nbas_face; jj++) 
        { 
            j = bases[jj];
            gj = bas_v[j];

            for (k = 0; k < Dim; k++)
            {
            F[j][k][i][k] += area*(*w) * gj*gi * EQU_SCALING;
            }
        }
        }
        }
        }
        }


#if 1
	for (i = 0; i < M; i++) 
    {
        for (k = 0; k < Dim; k++)
        {
        btype[i][k] = phgDofGetElementBoundaryType(nodal_force, e, i*Dim+k);
	    if (!(btype[i][k] & BC_BOTTOM_GRD)) 
        {
            //bzero(A[i], N*sizeof(A[0][0]));
            F[i][k][i][k] = 1.;
	    } 
        }
    }
#endif



	for (i = 0; i < M; i++)
	    for (k = 0; k < Dim; k++)
		Iu[i][k] = phgMapE2L(mat->cmap, 0, e, i * Dim + k);

    for (i = 0; i < M; i++)
    for (k = 0; k < Dim; k++)
        phgMatAddEntries(mat, 1, Iu[i] + k, M*Dim, Iu[0], &(F[i][k][0][0]));

    }
    //printf("done assembly matrix !!\n");

    phgMatAssemble(mat);
    solver->mat = mat;
    //solver = phgMat2Solver(SOLVER_SUPERLU, mat);

	solver->rhs->assembled = TRUE;
	phgVecAXPBY(1., vec_res_u, 0, &solver->rhs);
    //solver->rhs = phgVecCopy(vec_res_u, NULL);
    phgVecDump(solver->rhs, "rhs.txt");
    //printf("done copy right hand side vector!!\n");

    
    //VEC *result = phgMapCreateVec(map_res_u, 1);
    //phgSolverVecSolve(solver, FALSE, result);

    phgSolverSolve(solver, TRUE, nodal_force, NULL);

    phgExportVTK(g, "nodal_force.vtk", nodal_force, NULL);

    phgDofCopy(nodal_force, &ns->nodal_force_value, NULL, "nodal_force_value");

    phgDofFree(&nodal_force);

}

void get_contact_force_value(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    INT s, i,j,  k, l, q;

    DOF *res_u = ns->nodal_force_value;
    DOF *wp = ns->water_pressure;

    DOF *compare = phgDofNew(g, DOF_P1, 1, "compare", DofNoAction);
    ForAllElements(g, e)
    {

        for (s = 0; s < NFace; s++)
        {

            if (e->bound_type[s] & BC_BOTTOM)
            {
                const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                for (i = 0; i < 3; i++)
                {
                    INT v = GetFaceVertex(s, i);
                    INT local_idx = e->verts[v];

                    FLOAT *res_u_data = DofVertexData(res_u, local_idx);
                    FLOAT *wp_data = DofVertexData(wp, local_idx);

                    FLOAT *compare_data = DofVertexData(compare, local_idx);

                    FLOAT res_u_value = res_u_data[0]*normal[2]+res_u_data[1]*normal[0]+res_u_data[2]*normal[1];
                    //FLOAT wp_value = INNER_PRODUCT(wp_data, normal);
                    FLOAT wp_value = *wp_data; 

                    *compare_data = -fabs(res_u_value) + fabs(wp_value);

                }
            }
        }
    }

    phgExportVTK(g, "compare.vtk", compare, NULL);
	phgDofCopy(compare, &ns->contact_force_value, NULL, "contact_force");

    phgDofFree(&compare);
    //phgDofFree(&res_u);
    //phgDofFree(&wp);
}

void get_contact_force(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    INT s, i,j,  k, l, q;

    DOF *res_u = ns->nodal_force;
    DOF *wp = ns->water_force;

    DOF *compare = phgDofNew(g, DOF_P1, 1, "compare", DofNoAction);
    ForAllElements(g, e)
    {

        for (s = 0; s < NFace; s++)
        {

            if (e->bound_type[s] & BC_BOTTOM)
            {
                const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                for (i = 0; i < 3; i++)
                {
                    INT v = GetFaceVertex(s, i);
                    INT local_idx = e->verts[v];

                    FLOAT *res_u_data = DofVertexData(res_u, local_idx);
                    FLOAT *wp_data = DofVertexData(wp, local_idx);

                    FLOAT *compare_data = DofVertexData(compare, local_idx);

                    //FLOAT res_u_value = res_u_data[0]*normal[2]+res_u_data[1]*normal[0]+res_u_data[2]*normal[1];
                    //printf("diff %e\n", res_u_data[1]*normal[0]+res_u_data[2]*normal[1]);
                    FLOAT res_u_value = res_u_data[0];
                    FLOAT wp_value = wp_data[0];

                    *compare_data = -fabs(res_u_value) + fabs(wp_value);

                }
            }
        }
    }

    phgExportVTK(g, "net_nodal_loads.vtk", compare, NULL);
	phgDofCopy(compare, &ns->contact_force, NULL, "contact_force");
    phgDofFree(&compare);
    //phgDofFree(&res_u);
    //phgDofFree(&wp);
}

void get_contact(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    INT s, i,j,  k, l, q;

    MAP *Vmap, *Pmap;
    MAT *pmat[9];

    DOF *res_u = phgDofNew(ns->g, DOF_P2, Dim, "res_u", DofNoAction);
    DOF *res_p = phgDofNew(ns->g, DOF_P1, 1, "res_p", DofNoAction);


    phgNSDestroySolverU(ns);

#if 1
    ns->Vmap = phgMapCreate(ns->du, NULL);
    ns->Pmap = phgMapCreate(ns->dp, NULL);

    Vmap = ns->Vmap;
    Pmap = ns->Pmap;

    ns->matF = phgMapCreateMat(Vmap, Vmap);
    ns->matBt = phgMapCreateMat(Vmap, Pmap);
    ns->matB = phgMapCreateMat(Pmap, Vmap);
    ns->matC = phgMapCreateMat(Pmap, Pmap);

    ns->matF->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;
    if (ns->matC != NULL)
	ns->matC->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;

    ns->matB->handle_bdry_eqns = FALSE;
    ns->matBt->handle_bdry_eqns = FALSE;

    pmat[0] = ns->matF;
    pmat[1] = ns->matBt;
    pmat[2] = ns->matB;
    pmat[3] = ns->matC;

    ns->matNS = phgMatCreateBlockMatrix(g, 2, 2, pmat, NULL, NULL);
    ns->matNS->mv_data = phgAlloc(sizeof(*ns->matNS->mv_data));
    ns->matNS->mv_data[0] = (void *) ns;
    ns->matNS->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;

    phgOptionsPush();
    //phgOptionsSetOptions(_nsp->Stokes_opts);
    ns->solver_u = phgSolverCreate(SOLVER_DEFAULT, ns->du, ns->dp, NULL);
    phgOptionsPop();

    phgMatDestroy(&ns->solver_u->mat);
    ns->solver_u->mat = ns->matNS;
    ns->solver_u->rhs->mat = ns->solver_u->mat;
    //phgDofSetDataByValue(ns->du, 0.);
    //phgDofSetDataByValue(ns->dp, 0.);
        ns->ltype = PICARD;
#endif
            //phgNSDestroySolverU(ns);
    //phgNSInitSolverU(ns);
    phgNSBuildSolverURHS(ns, 0, 0, 0);
    //phgNSBuildSolverURHS(ns);
    phgVecDumpMATLAB(ns->solver_u->rhs, "rhs1", "rhs1_.m");
    phgNSBuildSolverUMat(ns, 0, 0, 0);
    //phgNSBuildSolverUMat(ns);
            phgMatDumpMATLAB(ns->matB, "mat1", "mat1_.m");

    MAT *all_mat[4];

    all_mat[0] = ns->matF;
    all_mat[1] = ns->matBt;
    all_mat[2] = ns->matB;
    all_mat[3] = ns->matC;

    MAT *u_mat = phgMatCreateBlockMatrix(g, 2, 2, all_mat, NULL, NULL);

    MAT *lhs = ns->matNS;
    VEC *rhs = ns->solver_u->rhs;


    VEC *res = phgVecCopy(rhs, NULL);

    DOF *sols[2] = {ns->du, ns->dp};
    //DOF *sols[2] = {velocity, pressure};

    DOF *res_u_p[2] = {res_u, res_p};

    VEC *vec_sol = phgMapCreateVec(rhs->map, 1);
        
    phgMapDofToLocalData(rhs->map, 2, &sols[0], vec_sol->data);


    //VEC *vec_sol1 = phgMapCreateVec(rhs->map, 1);

#if 0
    VEC *vec_sol1 = phgMatVec(1, 1, u_mat, rhs, 0, NULL);
    phgVecDumpMATLAB(vec_sol, "vec_sol", "vec_sol.m");
    phgVecDumpMATLAB(vec_sol1, "vec_sol1", "vec_sol1.m");

#endif
    phgMatVec(0, 1, u_mat, vec_sol, -1, &res);
    
    phgMapLocalDataToDof(rhs->map, 2, &res_u_p[0], res->data);
    phgExportVTK(ns->g, "res.vtk", res_u, res_p, NULL); 


    phgPrintf("\n Compute water pressure loads !!\n");

    DOF *wp = phgDofNew(g, DOF_P2, Dim, "water pressure", DofNoAction);
    MAP *map_wp = phgMapCreate(wp, NULL);
    VEC *rhs_wp = phgMapCreateVec(map_wp, Dim);
    phgVecDisassemble(rhs_wp);

    INT M = wp->type->nbas;
    BTYPE btype[M*Dim];
    ForAllElements(g, e)
    {
        FLOAT rhs_wp_e[M][Dim];
        INT local_map_idx[M][Dim];

        memset(rhs_wp_e, 0, sizeof rhs_wp_e);

        for (s = 0; s < NFace; s++)
        {

            if (e->bound_type[s] & BC_BOTTOM_GRD)
            {

                FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
                INT M_face = 3*(wp->type->np_vert + wp->type->np_edge);
                SHORT bas_idx_e[M_face];
                INT quad_order = 5;
                phgDofGetBasesOnFace(wp, e, s, bas_idx_e);
                
                const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0;

                QUAD *quad = phgQuadGetQuad2D(quad_order);
                FLOAT *p = quad->points;
                FLOAT *w = quad->weights;
                FLOAT area = phgGeomGetFaceArea(g, e, s);

                for (q = 0; q < quad->npoints; q++)
                {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);

                    phgGeomLambda2XYZ(g, e, lambda, &lambda_x, &lambda_y, &lambda_z);
                    const FLOAT *bas = wp->type->BasFuncs(wp, e, 0, -1, lambda);

                    for (i = 0; i < M_face; i++)
                    {
                        INT i_e = bas_idx_e[i];

                        for (k = 0; k < Dim; k++)
                        {
                            //if (lambda_z < 0)
                            {
                                rhs_wp_e[i_e][k] += area*w[q]*RHO_WAT*GRAVITY*(-lambda_z)*normal[k]*bas[i_e]/1.0e5;
                                //if (bas[i_e] < 0)
                                    //printf("!!!!!!!!!!!!!!!!!!! wrong bas !!!!!!!! %e\n", bas[i_e]);
                            }
                            //else
                            //    rhs_wp_e[i_e][k] += 0;
                        }
                    }
                }

            }
        }

        for (i = 0; i < M; i++)
        for (k = 0; k < Dim; k++)
        {
            local_map_idx[i][k] = phgMapE2L(map_wp, 0, e, i*Dim+k);
        }
        phgVecAddEntries(rhs_wp, 0, M*Dim, local_map_idx[0], &rhs_wp_e[0][0]);
    }

    phgVecAssemble(rhs_wp);
    phgMapLocalDataToDof(map_wp, 1, &wp, rhs_wp->data);
    phgExportVTK(g, "wp.vtk", wp, NULL);



    DOF *compare = phgDofNew(g, DOF_P1, 1, "compare", DofNoAction);
    ForAllElements(g, e)
    {

        for (s = 0; s < NFace; s++)
        {

            if (e->bound_type[s] & BC_BOTTOM)
            {
                const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                for (i = 0; i < 3; i++)
                {
                    INT v = GetFaceVertex(s, i);
                    INT local_idx = e->verts[v];

                    FLOAT *res_u_data = DofVertexData(res_u, local_idx);
                    //printf("res_u %e %e %e \n", res_u_data[0], res_u_data[1], res_u_data[2]);
                    //printf("normal %e %e %e \n", normal[0], normal[1], normal[2]);
                    FLOAT *wp_data = DofVertexData(wp, local_idx);
                    //printf("%e %e %e \n", wp_data[0], wp_data[1], wp_data[2]);

                    FLOAT *compare_data = DofVertexData(compare, local_idx);

                    FLOAT res_u_value = res_u_data[0]*normal[2]+res_u_data[1]*normal[0]+res_u_data[2]*normal[1];
                    //FLOAT res_u_value = res_u_data[0];
                    FLOAT wp_value = INNER_PRODUCT(wp_data, normal);
                    //FLOAT wp_value = wp_data[2]*normal[2];

                    //*compare_data = res_u_value - wp_value;
                    *compare_data = fabs(res_u_value) - fabs(wp_value);
                    //printf("!!!!!!!!! %e\n", *compare_data);
                    //if (*compare_data < 0)
                    //printf("%e %e\n", res_u_value, wp_value);

                }
            }
        }
    }

    phgExportVTK(g, "compare.vtk", compare, NULL);
}

void
get_stress(NSSolver *ns, DOF *gradu, DOF *pressure)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, k, q;
    DOF *stress, *stress0;
    SOLVER *solver_stress;
    
    stress = phgDofNew(g, DOF_P1, DDim, "stress", DofNoAction);
    stress0 = phgDofNew(g, DOF_P1, 1, "stress", DofNoAction);

    /* solver_Gu */
    phgOptionsPush();
    phgOptionsSetOptions("-solver gmres "
			 "-solver_maxit 100 "
			 "-solver_rtol 1e-10");
#if 0
    PetscOptionsInsertFile(MPI_COMM_WORLD, "../options/petsc_cg.opts", PETSC_TRUE); 
#endif
    //phgOptionsSetOptions(Gu_opts);
    solver_stress = phgSolverCreate(SOLVER_DEFAULT, stress0, NULL);
    solver_stress->verb = 0;
    phgVerbosity = 0;
    phgOptionsPop();

    VEC *vec[DDim];
    for (k = 0; k < DDim; k++) {
	vec[k] = phgMapCreateVec(solver_stress->rhs->map, 1) ;
	phgVecDisassemble(vec[k]);
    }

    /* Build linear system */
    ForAllElements(g, e) {
	int M = stress0->type->nbas;	/* num of bases of Velocity */
	int order = DofTypeOrder(stress0, e) * 2;
	FLOAT A[M][M], rhs[DDim][M];
	INT I[M];
	QUAD *quad;
	FLOAT vol, det;
	const FLOAT *w, *p, *gu, *pres;
    FLOAT stressD[DDim];

	Bzero(A); Bzero(rhs);
	quad = phgQuadGetQuad3D(order);
	gu = phgQuadGetDofValues(e, gradu, quad); 
	pres  = phgQuadGetDofValues(e, pressure, quad); 

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);

	    for (i = 0; i < M; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, stress0, i, quad) + q;    
		for (j = 0; j < M; j++) {
		    const FLOAT *gj = phgQuadGetBasisValues(e, stress0, j, quad) + q;       
		    FLOAT qmass = vol*(*w) * (*gj) * (*gi);
			A[i][j] += qmass;
		}
		    
        FLOAT strain_rate[DDim] = {gu[0], 0.5*(gu[1]+gu[3]), 0.5*(gu[2]+gu[6]),
                            0.5*(gu[1]+gu[3]), gu[1], 0.5*(gu[5]+gu[7]),
                            0.5*(gu[2]+gu[6]), 0.5*(gu[5]+gu[7]), gu[2]};

        FLOAT visc = get_effective_viscosity(gu, 273.15, 0, VIS_STRAIN);


        for (k = 0; k < DDim; k++)
            stressD[k] = visc*strain_rate[k];

        FLOAT sigma[DDim] = {stressD[0] - (*pres)*PRES_SCALING, stressD[1], stressD[2],
                    stressD[3], stressD[4] - (*pres)*PRES_SCALING, stressD[5],
                    stressD[6], stressD[7], stressD[8] - (*pres)*PRES_SCALING};

		for (k = 0; k < DDim; k++)
        {
		    rhs[k][i] += vol*(*w) * sigma[k] * (*gi); 
            //printf("sigma: %e\n", sigma[k]);
        }
	    }
	    gu += DDim;
	    w++; p += Dim + 1;
	}

	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    I[i] = phgMapE2L(solver_stress->mat->cmap, 0, e, i);

	/* Global res */
	for (i = 0; i < M; i++)
	    phgMatAddEntries(solver_stress->mat, 1, I + i, M, I,
			     &(A[i][0])); 

	for (k = 0; k < DDim; k++)
	    phgVecAddEntries(vec[k], 0, M, I, &rhs[k][0]);
    }				/* end element */
    

    for (k = 0; k < DDim; k++)
	phgVecAssemble(vec[k]);
    solver_stress->rhs_updated = TRUE;

    INT n = DofGetDataCount(stress0);
    for (k = 0; k < DDim; k++) {
	phgVecCopy(vec[k], &solver_stress->rhs);
	phgDofSetDataByValue(stress0, 0.);
	phgSolverSolve(solver_stress, FALSE, stress0, NULL);
	phgPrintf("      solver_stress: nits = %d, resid = %0.4lg ",
		  solver_stress->nits, solver_stress->residual);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	FLOAT *vg = stress->data, *vg0 = stress0->data;
	for (i = 0; i < n; i++, vg0++, vg += DDim)
	    vg[k] = *vg0;
    }

    if (DUMP_MAT_VEC) {
	phgPrintf("Dumping MatGu, rhsGu\n");
	phgMatDumpMATLAB(solver_stress->mat, "A_gu", "mat_gu_.m");
	phgVecDumpMATLAB(solver_stress->rhs, "b_gu", "rhs_gu_.m");
    }

    phgInfo(2, "   Destroy solver Gu\n");
    for (k = 0; k < DDim; k++)
	phgVecDestroy(&vec[k]);
    phgSolverDestroy(&solver_stress);

    phgDofFree(&stress0);

	phgDofCopy(stress, &ns->stress, NULL, "stress");
	phgDofFree(&stress);
}

void
get_water_pressure_value(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, k, l, s, q;
    int i_face, j_face;

    DOF *wp = phgDofNew(g, DOF_P1, 1, "water pressure", DofNoAction);

#if 1
    phgOptionsPush();
    phgOptionsSetOptions("-solver gmres "
			 "-gmres_pc_type jacobi "
			 "-solver_maxit 10000 "
			 "-solver_rtol 1e-12");
    SOLVER *solver = phgSolverCreate(SOLVER_DEFAULT, wp, NULL);
    phgVerbosity = 0;
    phgOptionsPop();
#endif


    ForAllElements(g, e)
    {
        int M = wp->type->nbas;	
        FLOAT A[M][M], rhs[M], buffer[M];
        INT I[M];
        BOOLEAN btype[M];
        INT local_map_idx[M][Dim];
        QUAD *quad;
        FLOAT *p, *w;

        memset(A, 0, sizeof A); 
        memset(rhs, 0, sizeof rhs); 
        memset(btype, 0, sizeof btype); 

        for (s = 0; s < NFace; s++)
        {
            int M_face = 3*(wp->type->np_vert + wp->type->np_edge);
            SHORT bases[M_face];
            FLOAT x, y, z, lambda[Dim + 1], area, vu[Dim];
            FLOAT water_pressure;
            //gi, gj, ggi[Dim], ggj[Dim], qa;
            const FLOAT *normal, *Gi, *Gj;

            if (e->bound_type[s] & BC_BOTTOM)
            {
                phgDofGetBasesOnFace(wp, e, s, bases);
                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0.;

                area = phgGeomGetFaceArea(g, e, s);
                normal = phgGeomGetFaceOutNormal(g, e, s);
                quad = phgQuadGetQuad2D(5);

                p = quad->points;
                w = quad->weights;


                for (q = 0; q < quad->npoints; q++)
                {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);

                    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);

            const FLOAT *bas = wp->type->BasFuncs(wp, e, 0, -1, lambda);
                    for (i = 0; i < M_face; i++)
                    {
                        i_face = bases[i];
                        //const FLOAT *gi = sur_normal_x->type->BasFuncs(sur_normal_x, e, i_face, i_face + 1, lambda);

                        FLOAT a1 = bas[i_face];

                        for (j = 0; j < M_face; j++)
                        {
                            j_face = bases[j];
                            //const FLOAT *gj = sur_normal_x->type->BasFuncs(sur_normal_x, e, j_face, j_face + 1, lambda);

                            //FLOAT qmass = area*(*w) * (*gj) * (*gi);
                            A[i_face][j_face] += area*(*w)*bas[i_face]*bas[j_face];

                        }
                        if (a1 != bas[i_face])
                            printf("%f %f\n", a1, bas[i_face]);

                        if (z < 0)
                            water_pressure = -RHO_WAT*GRAVITY*z;
                        else
                            water_pressure = 0;

                        rhs[i_face] += area*(*w) * water_pressure * bas[i_face]; 
                    }

                    w++;
                }
            }
        }

        for (i = 0; i < M; i++) 
        {
            btype[i] = phgDofGetElementBoundaryType(wp, e, i);
            if (!(btype[i] & BC_BOTTOM))
            {
                A[i][i] = 1.;
                rhs[i] = 0.;
            } 
        }

        for (i = 0; i < M; i++)
            I[i] = phgMapE2L(solver->mat->cmap, 0, e, i);

        for (i = 0; i < M; i++) 
            phgSolverAddMatrixEntries(solver, 1, I+i, M, I, A[i]); 

        phgSolverAddRHSEntries(solver, M, I, rhs);
    }

    phgSolverSolve(solver, TRUE, wp, NULL);
    phgExportVTK(g, "wp.vtk", wp, NULL);

	phgDofCopy(wp, &ns->water_pressure, NULL, "water_pressure");
	phgDofFree(&wp);
}

void
get_normal_stress_value(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, k, l, s, q;
    int i_face, j_face;

    DOF *sigma_nn = phgDofNew(g, DOF_P1, 1, "stress normal", DofNoAction);

#if 1
    phgOptionsPush();
    phgOptionsSetOptions("-solver gmres "
			 "-gmres_pc_type jacobi "
			 "-solver_maxit 10000 "
			 "-solver_rtol 1e-12");
    SOLVER *solver = phgSolverCreate(SOLVER_DEFAULT, sigma_nn, NULL);
    phgVerbosity = 0;
    phgOptionsPop();
#endif


    ForAllElements(g, e)
    {
        int M = sigma_nn->type->nbas;	
        FLOAT A[M][M], rhs[M], buffer[M];
        INT I[M];
        BOOLEAN btype[M];
        INT local_map_idx[M][Dim];
        QUAD *quad;
        FLOAT *p, *w;

        memset(A, 0, sizeof A); 
        memset(rhs, 0, sizeof rhs); 
        memset(btype, 0, sizeof btype); 

        for (s = 0; s < NFace; s++)
        {
            int M_face = 3*(sigma_nn->type->np_vert + sigma_nn->type->np_edge);
            SHORT bases[M_face];
            FLOAT x, y, z, lambda[Dim + 1], area, vu[Dim];
            FLOAT water_pressure, stress_nn;
            //gi, gj, ggi[Dim], ggj[Dim], qa;
            const FLOAT *normal, *Gi, *Gj;

            if (e->bound_type[s] & BC_BOTTOM)
            {
                phgDofGetBasesOnFace(sigma_nn, e, s, bases);
                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0.;

                area = phgGeomGetFaceArea(g, e, s);
                normal = phgGeomGetFaceOutNormal(g, e, s);
                quad = phgQuadGetQuad2D(5);

                p = quad->points;
                w = quad->weights;


                for (q = 0; q < quad->npoints; q++)
                {
                    FLOAT vs[DDim];
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);

                    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
                    phgDofEval(ns->stress, e, lambda, vs);

            const FLOAT *bas = sigma_nn->type->BasFuncs(sigma_nn, e, 0, -1, lambda);
                    for (i = 0; i < M_face; i++)
                    {
                        i_face = bases[i];

                        FLOAT a1 = bas[i_face];

                        for (j = 0; j < M_face; j++)
                        {
                            j_face = bases[j];

                            A[i_face][j_face] += area*(*w)*bas[i_face]*bas[j_face];

                        }
                        if (a1 != bas[i_face])
                            printf("%f %f\n", a1, bas[i_face]);

                        stress_nn = (vs[0]*normal[0]+vs[1]*normal[1]+vs[2]*normal[2])*normal[0] + (vs[3]*normal[0]+vs[4]*normal[1]+vs[5]*normal[2])*normal[1] + (vs[6]*normal[0]+vs[7]*normal[1]+vs[8]*normal[2])*normal[2];
                        if (z < 0)
                            water_pressure = -RHO_WAT*GRAVITY*z;
                        else
                            water_pressure = 0;

                        rhs[i_face] += area*(*w) * stress_nn * bas[i_face]; 
                    }

                    w++;
                }
            }
        }

        for (i = 0; i < M; i++) 
        {
            btype[i] = phgDofGetElementBoundaryType(sigma_nn, e, i);
            if (!(btype[i] & BC_BOTTOM))
            {
                A[i][i] = 1.;
                rhs[i] = 0.;
            } 
        }

        for (i = 0; i < M; i++)
            I[i] = phgMapE2L(solver->mat->cmap, 0, e, i);

        for (i = 0; i < M; i++) 
            phgSolverAddMatrixEntries(solver, 1, I+i, M, I, A[i]); 

        phgSolverAddRHSEntries(solver, M, I, rhs);
    }

    phgSolverSolve(solver, TRUE, sigma_nn, NULL);
    phgExportVTK(g, "sigma_nn.vtk", sigma_nn, NULL);

	phgDofCopy(sigma_nn, &ns->stress_nn, NULL, "sigma_nn");
	phgDofFree(&sigma_nn);
}

INT if_update_shelf_mask(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    //DOF **p;

    //p = ns->p;

    INT s, i;

    INT IF_CHANGE_MASK = 0, IF_CHANGE_MASK0;

    INT contact_method = 1;

    // after u is solved, we check the mask
    if (contact_method == 0)
    {
    //get_stress(ns, ns->gradu[1], ns->p[1]);
    //get_water_pressure_value(ns);
    //get_normal_stress_value(ns);

    get_avg_gu(ns);
    get_stress1(ns);
    //get_smooth_surface_values(ns, ns->stress_nn1, 0);
    get_water_pressure_value1(ns);

    //DOF *cw = compare_two_dofs(ns->water_pressure, ns->water_pressure1);
    //DOF *csn = compare_two_dofs(ns->stress_nn, ns->stress_nn1);
    //phgExportVTK(g, "cmp.vtk", cw, csn, NULL);


    //phgDofAXPY(1.0, ns->water_pressure, &(ns->stress_nn));
    phgDofAXPY(1.0, ns->water_pressure1, &(ns->stress_nn1));

    ns->stress_nn = ns->stress_nn1;
    // now stress_nn has the stress_nn + p_n data we need
    // stress_nn < 0, p_n > 0
    // if floating, stress_nn + p_n >= 0

    //phgExportVTK(g, "diff.vtk", ns->stress_nn, NULL);
    phgExportVTK(g, "diff1.vtk", ns->stress_nn1, ns->water_pressure1, NULL);

    }
    
    if (contact_method == 1)
    {
        get_water_force(ns);
        get_contact_force(ns);
        ns->stress_nn = ns->contact_force;
        get_smooth_surface_values(ns, ns->stress_nn, 0);
    }

    if (contact_method == 2)
    {
        get_water_pressure_value(ns);
        get_contact_force_value(ns);
        phgDofAXPY(1.0, ns->water_pressure, &(ns->contact_force_value));
        ns->stress_nn = ns->contact_force_value;
    }


    ForAllElements(g, e)
    {
        for (s = 0; s < NFace; s++)
        {
            if (e->bound_type[s] & BC_BOTTOM_GRD) 
            {
                
                QUAD *quad = phgQuadGetQuad2D(1);
                FLOAT *p = quad->points;
                FLOAT lambda[Dim + 1], x,y,z;
                int v0, v1, v2, local_idx[3];
                FLOAT bed_z;
                FLOAT *net_sigma_nn[3], sum_sigma_nn;
                v0 = GetFaceVertex(s, 0);
                v1 = GetFaceVertex(s, 1);
                v2 = GetFaceVertex(s, 2);
                local_idx[0] = e->verts[v0];
                local_idx[1] = e->verts[v1];
                local_idx[2] = e->verts[v2];

                net_sigma_nn[0] = DofVertexData(ns->stress_nn, local_idx[0]);
                net_sigma_nn[1] = DofVertexData(ns->stress_nn, local_idx[1]);
                net_sigma_nn[2] = DofVertexData(ns->stress_nn, local_idx[2]);

                sum_sigma_nn = *net_sigma_nn[0] + *net_sigma_nn[1] + *net_sigma_nn[2];

                FLOAT x0 = g->verts[local_idx[0]][0];
                FLOAT x1 = g->verts[local_idx[1]][0];
                FLOAT x2 = g->verts[local_idx[2]][0];

                FLOAT x_lim = 500e3;

                if (x0 > x_lim && x1 > x_lim && x2 > x_lim)
                {
                    // we only take care of the region of x > 480 km
                if (*net_sigma_nn[0] >= 0 &&
                    *net_sigma_nn[1] >= 0 &&
                    *net_sigma_nn[2] >= 0)
                //if ((*net_sigma_nn[0] + *net_sigma_nn[1] + *net_sigma_nn[2])/3.0 >= 0)
                {
                    e->bound_type[s] = (BC_BOTTOM | BC_ISHELF);
                    // if any of the 3 nodal values is <= 0, we re-mark this face as afloat
                    IF_CHANGE_MASK = 1;
                }
                else
                {
                    e->bound_type[s] = (BC_BOTTOM | BC_BOTTOM_GRD);
                }
                }


            }
        }
    }

    phgUpdateBoundaryTypes(g);


    MPI_Allreduce(&IF_CHANGE_MASK, &IF_CHANGE_MASK0, 1, PHG_MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if (IF_CHANGE_MASK0 == 1)
    {
        phgPrintf("\n----------------------------------------\n");
        phgPrintf("Ice shelf masks were found and updated! \n");
        phgPrintf("We need to recalculate the velocity field !\n");
        phgPrintf("------------------------------------------\n");
    }
    else
    {
        phgPrintf("\n----------------------------------------\n");
        phgPrintf("The bottom mask was remained unchanged!\n");
        phgPrintf("We do not need to recalculate the velocity field!\n");
        phgPrintf("-----------------------------------------\n");
    }

    //phgDofFree(&p[0]);
    return IF_CHANGE_MASK0;

}

void get_stress1(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;

    int i, s, k;
    INT local_idx, local_idx1[Dim];
    FLOAT *gu, *pres, *data, visc, stressD[DDim], sigma[DDim];

    DOF *gradu = ns->avg_gu;
    DOF *pressure = ns->p[1];

    DOF *stress = phgDofNew(g, DOF_P1, DDim, "stress", DofNoAction);
    DOF *stress_nn = phgDofNew(g, DOF_P1, 1, "stress", DofNoAction);

    DOF *avg_n = phgDofNew(g, DOF_P1, Dim, "avg_n", DofNoAction);
    get_avg_n(g, avg_n);

    ForAllElements(g, e)
    {
        for (i = 0; i < Dim+1; i++)
        {

            local_idx = e->verts[i];
            gu = DofVertexData(gradu, local_idx);
            pres = DofVertexData(pressure, local_idx);

            data = DofVertexData(stress, local_idx);

        FLOAT strain_rate[DDim] = {gu[0], 0.5*(gu[1]+gu[3]), 0.5*(gu[2]+gu[6]),
                            0.5*(gu[1]+gu[3]), gu[1], 0.5*(gu[5]+gu[7]),
                            0.5*(gu[2]+gu[6]), 0.5*(gu[5]+gu[7]), gu[2]};

        visc = get_effective_viscosity(gu, 273.15, 0, VIS_STRAIN);


        for (k = 0; k < DDim; k++)
            stressD[k] = visc*strain_rate[k];

        sigma[0] = stressD[0] - (*pres)*PRES_SCALING;
        sigma[1] = stressD[1];
        sigma[2] = stressD[2];
        sigma[3] = stressD[3];
        sigma[4] = stressD[4] - (*pres)*PRES_SCALING;
        sigma[5] = stressD[5];
        sigma[6] = stressD[6];
        sigma[7] = stressD[7];
        sigma[8] = stressD[8] - (*pres)*PRES_SCALING;

        for (k = 0; k < DDim; k++)
            data[k] = sigma[k];
        }


    }

    FLOAT *vs[Dim], vsn[Dim];

    ForAllElements(g, e)
    {
        for (s = 0; s < NFace; s++)
        {
            if (e->bound_type[s] & BC_BOTTOM)
            {
                //const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);
                
                FLOAT *normal[3];

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);

                local_idx1[0] = e->verts[v0];
                local_idx1[1] = e->verts[v1];
                local_idx1[2] = e->verts[v2];


                for (k = 0; k < Dim; k++)
                {

                    normal[k] = DofVertexData(avg_n, local_idx1[k]);

                vs[k] = DofVertexData(stress, local_idx1[k]);


                //vsn[k] = (vs[k][0]*normal[0]+vs[k][1]*normal[1]+vs[k][2]*normal[2])*normal[0] + 
                //    (vs[k][3]*normal[0]+vs[k][4]*normal[1]+vs[k][5]*normal[2])*normal[1] + 
                 //   (vs[k][6]*normal[0]+vs[k][7]*normal[1]+vs[k][8]*normal[2])*normal[2];
                vsn[k] = (vs[k][0]*normal[k][0]+vs[k][1]*normal[k][1]+vs[k][2]*normal[k][2])*normal[k][0] + 
                    (vs[k][3]*normal[k][0]+vs[k][4]*normal[k][1]+vs[k][5]*normal[k][2])*normal[k][1] + 
                    (vs[k][6]*normal[k][0]+vs[k][7]*normal[k][1]+vs[k][8]*normal[k][2])*normal[k][2];
                DofVertexData(stress_nn, local_idx1[k])[0] = vsn[k];

                }
            }
        }
    }

    phgDofCopy(stress, &ns->stress1, NULL, NULL);
    phgDofCopy(stress_nn, &ns->stress_nn1, NULL, NULL);
    phgExportVTK(g, "stress_nn.vtk", stress_nn, NULL);

    phgDofFree(&stress);
    phgDofFree(&stress_nn);
    phgDofFree(&avg_n);
}

void get_water_pressure_value1(NSSolver *ns)
{

    GRID *g = ns->g;
    SIMPLEX *e;
    INT local_idx;
    FLOAT z, wp;

    int i, s, k;

    DOF *water = phgDofNew(g, DOF_P1, 1, "water", DofNoAction);

    ForAllElements(g, e)
    {
        for (i = 0; i < Dim+1; i++)
        {
            local_idx = e->verts[i];
            z = g->verts[local_idx][2];

            if (z < 0)
                wp = RHO_WAT*GRAVITY*fabs(z);
            else
                wp = 0;

            DofVertexData(water, local_idx)[0]
                = wp;
        }
    }

    phgDofCopy(water, &ns->water_pressure1, NULL, NULL);
    phgDofFree(&water);
}
