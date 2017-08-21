#include "ins.h"


FLOAT trans_eye[DDim] = {
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
};

void func_eye(FLOAT x, FLOAT y, FLOAT z, FLOAT *e) {
    e[0] = 1.;
    e[1] = 0.;
    e[2] = 0.;
    e[3] = 0.;
    e[4] = 1.;
    e[5] = 0.;
    e[6] = 0.;
    e[7] = 0.;
    e[8] = 1.;
}



/* Get local bases for boundary vertex i,
 *  return num of constrained bases, 
 *  and noramlized direction of rotated bases.
 *  */
SURF_BAS *
get_surface_bases(GRID *g, DOF_TYPE *u_type)
{
    SIMPLEX *e;
    DOF *surf_dof = NULL, *norm_lat = NULL, *norm_bot = NULL;
    BOOLEAN *rotated = NULL;
    INT nrot = 0;
    surf_dof = phgDofNew(g, u_type, DDim, "Surf bases", DofNoAction);
    
    //norm_lat = phgDofNew(g, u_type, Dim, "Norm lateral", DofNoAction);
    //norm_bot = phgDofNew(g, u_type, Dim, "Norm Bottom", DofNoAction);
    //DOF *coord = phgDofNew(g, u_type, Dim, "coord", func_xyz_);
	DOF *avg_n = phgDofNew(g, DOF_P2, 3, "avg n", DofNoAction);
        get_avg_n(g, avg_n);
    
    rotated = phgCalloc(DofGetDataCount(surf_dof) / (DDim), sizeof(*rotated));
    SURF_BAS *surf_bas;

    surf_bas = phgCalloc(1, sizeof(*surf_bas));
    surf_bas->type = u_type;
    //surf_bas->dof = phgDofNew(g, u_type, DDim, "Surf bases", DofNoAction);//surf_dof;
    //surf_bas->dof = surf_dof;
    surf_bas->rotated = rotated;
    phgDofSetDataByValue(surf_dof, 0.);
    //phgDofSetDataByValue(surf_bas->dof, 0.);
    
    //phgDofSetDataByValue(norm_lat, -99.);
    //phgDofSetDataByValue(norm_bot, -99.);
    
        ForAllElements(g, e) {
	int s, ii, i, m, dof_i;
	int N = surf_dof->type->nbas;
	//int N = surf_bas->dof->type->nbas;
	FLOAT *avg_n_v;
	FLOAT normal[Dim];
	int  v[3];
	FLOAT norm;
	FLOAT norm_value_lat[N][Dim];
	FLOAT norm_value_bot[N][Dim];
	FLOAT bas_value[N][DDim],  H[Dim][Dim], c[Dim], 
	    bxyz[Dim][Dim] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

	phgDofGetElementData(surf_dof, e, &bas_value[0][0]);
	//phgDofGetElementData(surf_bas->dof, e, &bas_value[0][0]);
	
	for (s = 0; s < NFace; s++) { 
	    int nbas_face = NbasFace(surf_dof);
	    //int nbas_face = NbasFace(surf_bas->dof);
	    INT id; 
	    SHORT ibas_face[nbas_face];

        //FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

        //if (!((e->bound_type[s] & BC_BOTTOM) && !(e->bound_type[s] & BC_ISHELF)))
        if (!(e->bound_type[s] & BC_BOTTOM))
            continue;

	    phgDofGetBasesOnFace(surf_dof, e, s, ibas_face);
	    //phgDofGetBasesOnFace(surf_bas->dof, e, s, ibas_face);
	    for (ii = 0; ii < nbas_face; ii++) {
		i = ibas_face[ii];
		dof_i = phgDofMapE2D(avg_n, e, i*Dim);
		avg_n_v = DofData(avg_n);
		normal[0] = avg_n_v[dof_i + 0];
		normal[1] = avg_n_v[dof_i + 1];
		normal[2] = avg_n_v[dof_i + 2];
        /* i means the base function number in the element */
		/* Use Gram–Schmidt process to get orthogonal bases, 
		 * one constrains */

		id = phgDofMapE2D(surf_dof, e, i * (DDim)) / (DDim);
		//id = phgDofMapE2D(surf_bas->dof, e, i * (DDim)) / (DDim);
		rotated[id] = TRUE;
		nrot++;
		
		BTYPE elem_btype = phgDofGetElementBoundaryType(surf_dof, e, i*DDim);
		//BTYPE elem_btype = phgDofGetElementBoundaryType(surf_bas->dof, e, i*DDim);


		/* fisrt basis */
		memcpy(H[0], normal, Dim * sizeof(FLOAT));

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
		H[1][0] = 0;
		H[1][1] = 1;
		H[1][2] = 0;

		/* third basis */
		for (m++; m < Dim; m++)
		    if (fabs(c[0] = INNER_PRODUCT(H[0], bxyz[m])) < 0.9)
			break;
		assert(m < Dim);
		//c[1] = INNER_PRODUCT(H[1], bxyz[m]);
		c[1] = H[1][0]*(bxyz[m][0] - c[0] * H[0][0])+H[1][1]*(bxyz[m][1] - c[0] * H[0][1]) + H[1][2]*(bxyz[m][2] - c[0] * H[0][2]);
		H[2][0] = bxyz[m][0] - c[0] * H[0][0] - c[1] * H[1][0]; 
		H[2][1] = bxyz[m][1] - c[0] * H[0][1] - c[1] * H[1][1];  
		H[2][2] = bxyz[m][2] - c[0] * H[0][2] - c[1] * H[1][2];  
		H[2][0] = H[0][1]*H[1][2] - H[0][2]*H[1][1];
		H[2][1] = H[0][2]*H[1][0] - H[0][0]*H[1][2];
		H[2][2] = H[0][0]*H[1][1] - H[0][1]*H[1][0];
		norm = sqrt(INNER_PRODUCT(H[2], H[2]));
		assert(norm > 1e-10);
		H[2][0] /= norm;
		H[2][1] /= norm;
		H[2][2] /= norm;
        	
        if ((elem_btype & BC_BOTTOM_GRD) && (elem_btype & BC_DIVIDE))
        {
			
            H[1][0] = 1;
            H[1][1] = 0;
            H[1][2] = 0;

            H[2][0] = (H[0][1]*H[1][2]-H[0][2]*H[1][1]);
            H[2][1] = -(H[0][0]*H[1][2]-H[0][2]*H[1][0]);
            H[2][2] = (H[0][0]*H[1][1]-H[0][1]*H[1][0]);
            
        }


#if 0
#  warning check use only: xyz coord ----------------------------
		memcpy(bas_value[i], bxyz[0], DDim*sizeof(FLOAT));
#else
		memcpy(bas_value[i], H[0], DDim*sizeof(FLOAT));
        /* bas_value is contains all bas values of all nodes in the element. 
         * For the nodes at the boundaries, bas value is real number, otherwise
         * bas value is just 0 */
#endif
	    } /* end bas */
	}     /* end face */
	phgDofSetElementData(surf_dof, e, &bas_value[0][0]);
	//phgDofSetElementData(surf_bas->dof, e, &bas_value[0][0]);
    }	      /* end elem */

        phgDofCopy(surf_dof, &surf_bas->dof, NULL, "surf_dof");

    //INT nrot0 = nrot;
    //MPI_Allreduce(&nrot0, &nrot, 1, PHG_MPI_INT, MPI_SUM, g->comm);
    //phgPrintf("   Slip bdry Dofs: %d, ", nrot);

    //phgDofDump(surf_dof);
    PERIODIC_SYNC(surf_dof);
    //PERIODIC_SYNC(surf_bas->dof);
    phgDofFree(&avg_n);
    phgDofFree(&surf_dof);
    return surf_bas;
}






/*
 * A_{3, ncol} = Trans_{3,3} * A_{3, ncol}
 *  */
void 
trans_left(FLOAT *A, int ncol, int lda, const FLOAT *Trans) 
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


/*
 * A_{3, ncol} = Trans_{3,3}^T * A_{3, ncol}
 *  */
void 
trans_leftT(FLOAT *A, int ncol, int lda, const FLOAT *Trans)  
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


/*
 * A_{nrow, 3} = A_{nrow, 3} * Trans_{3,3}^T
 *  */
void 
trans_rightT(FLOAT *A, int nrow, int lda, const FLOAT *Trans) 
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



void 
rotate_dof_bases(DOF *u, SURF_BAS *surf_bas, BOOLEAN forward)
{
    INT i, N = DofGetDataCount(u) / Dim; 
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    FLOAT *v = DofData(u);
    const FLOAT *trans = DofData(surf_dof);

    for (i = 0; i < N; i++) { 
	if (rotated[i]) {
	    if (forward)
		trans_left(v, 1, 1, trans);
	    else
		trans_leftT(v, 1, 1, trans);
	}
	trans += DDim;
	v += Dim;
    }

    return;
}


void dof_set_normal_data(DOF *u_h, SURF_BAS *surf_bas)
{
    GRID *g = u_h->g;
    SIMPLEX *e;
    INT s;

    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);

    ForAllElements(g, e) {

    if (!(e->bound_type[s] & BC_BOTTOM_GRD))
        continue;

	int i, N = surf_dof->type->nbas;
	FLOAT vu[N][Dim], u0[Dim];

	phgDofGetElementData(u_h, e, vu[0]);
	for (i = 0; i < N; i++) {
	    INT id = phgDofMapE2D(surf_dof, e, i * (DDim)) / (DDim);
	    if (!rotated[id])
		continue;	

	    const FLOAT *trans = Trans + id*(Dim*Dim);
	    FLOAT *coord = phgDofGetElementCoordinates(u_h, e, i*Dim);
	    func_u(coord[0], coord[1], coord[2], u0);
	    trans_left  (vu[i], 1, 1, trans);
	    vu[i][0] = INNER_PRODUCT(u0, trans);
	    trans_leftT (vu[i], 1, 1, trans);
	}
	phgDofSetElementData(u_h, e, vu[0]);
    }	      /* end elem */

    return;
}
