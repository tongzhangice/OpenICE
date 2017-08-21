/*
 *  Subroutines of Building Mat & RHS of solver
 *
 *
 *  */
#include "ins.h"
#include "mat_op3.h"
#define _nsp (ns->ns_params)
#define _pcd (ns->pcd)
#define SIGN_FRICTION -1
#define SIGN_STAB -1.


static double m_idx = 1.0/3;
//static double C_mtp = 31651.7551654794; //1.0e7/pow(SEC_PER_YEAR, 1.0/3);
static double C_mtp = 31651.7551654794; //1.0e7/pow(SEC_PER_YEAR, 1.0/3);
//static double C_mtp = 1.0e4; //1.0e7/pow(SEC_PER_YEAR, 1.0/3);

FLOAT *get_bas_dot_normal(const FLOAT *bas, const FLOAT *normal, INT i_e, INT j_e)
{
	static FLOAT basdotn[Dim][Dim];

	basdotn[0][0] = normal[0]*normal[0]*bas[i_e]*bas[j_e];
	basdotn[0][1] = normal[1]*normal[0]*bas[i_e]*bas[j_e];
	basdotn[0][2] = normal[2]*normal[0]*bas[i_e]*bas[j_e];
	basdotn[1][0] = normal[0]*normal[1]*bas[i_e]*bas[j_e];
	basdotn[1][1] = normal[1]*normal[1]*bas[i_e]*bas[j_e];
	basdotn[1][2] = normal[2]*normal[1]*bas[i_e]*bas[j_e];
	basdotn[2][0] = normal[0]*normal[2]*bas[i_e]*bas[j_e];
	basdotn[2][1] = normal[1]*normal[2]*bas[i_e]*bas[j_e];
	basdotn[2][2] = normal[2]*normal[2]*bas[i_e]*bas[j_e];

	return basdotn[0];
}

/****************************************************************
 * Build RHS which is the residual of the nonlinear system.
 ***************************************************************/
void 
//phgNSBuildSolverURHS(NSSolver *ns, GEO_INFO *geo)
phgNSBuildSolverURHS(NSSolver *ns, INT IF_DB, INT nonstep, FLOAT Time)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    SOLVER *solver_u = ns->solver_u;
    int i, k, l, q, s;
    FLOAT *dt = ns->dt;
    BOOLEAN tstep_minus = (ns->u[-1] != NULL);
    VEC *vec_rhs = phgMapCreateVec(solver_u->rhs->map, 1);
#if USE_NODAL_LOADS
    ns->vec_rhs0 = phgMapCreateVec(solver_u->rhs->map, 1);
#endif
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    int viscosity_type = ns->viscosity_type;

//    DOF *ice_shelf_pres = phgDofNew(g, DOF_DG1, 1, "ice_shelf_pres", func_ice_shelf_pres);
#if USE_SLIDING_BC || CASE_DRAINAGE_BASIN
    SURF_BAS *surf_bas = ns->surf_bas;
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);
    assert(surf_dof->type == ns->du->type);
#endif
 
#if STEADY_STATE
    assert(fabs(Theta - 1) < 1e-12);
    Thet1 = 0; Unused(Thet1);
    Unused(dt);
#else
    Thet1 = 1 - Theta;
    Unused(dt);
#endif /* STEADY_STATE */

    phgPrintf("   DB_mask: [");
    for (k = 0; k < Dim; k++)
	phgPrintf("%d ", ns->du->DB_masks[k]);
    phgPrintf("]   ");

    nu_max = -1e100;
    nu_min = +1e100;

    phgVecDisassemble(vec_rhs);
#if USE_NODAL_LOADS
    phgVecDisassemble(ns->vec_rhs0);
#endif

	DOF *avg_n = phgDofNew(g, DOF_P2, 3, "avg n", DofNoAction);

    get_avg_n(g, avg_n);

    ForAllElements(g, e) {
	int M = ns->du->type->nbas;	/* num of bases of Velocity */
	int N = ns->dp->type->nbas;	/* num of bases of Pressure */
	int order = DofTypeOrder(ns->du, e) * 3 - 1; /* Note:
							*   quad order is really high here,
							*   highest order term (u \nabla u, phi)  */
	FLOAT bufu[M], bufp[N], rhsu[M][Dim], rhsp[N];
    FLOAT rhsu0[M][Dim], rhsp0[N];
	INT Iu[M][Dim], Ip[N];
	QUAD *quad;
	FLOAT vol, area, det;
	const FLOAT *w, *p, *normal,
	    **vu, *vu_queue[3],*ice_shelf_pres_value,
	    *vf[2], *gu[2], *vp[2], *vw, *vTe;
	FLOAT *vf_cache[2];
	FLOAT vp0 = 0.;

	vu = vu_queue + 1;

	quad = phgQuadGetQuad3D(order);
	vu[0] = phgQuadGetDofValues(e, ns->u[0], quad);	                /* u^{n} */
	vp[0] = phgQuadGetDofValues(e, ns->p[0], quad);	                /* p^{n} */
	gu[0] = phgQuadGetDofValues(e, ns->gradu[0], quad);             /* grad u^{n} */
	//vw = phgQuadGetDofValues(e, ns->wind, quad);                  /* wind */
	if (ns_params->noniter_temp)			                /* nonlinear temp iter */
	    vTe = phgQuadGetDofValues(e, ns->T[1], quad);	        /* T^{n+1} */
	else
	    vTe = phgQuadGetDofValues(e, ns->T[0], quad);	        /* T^{n} */

	if (tstep_minus) { 
	    vu[-1] = phgQuadGetDofValues(e, ns->u[-1], quad);            /* u^{n-1} */
	} else {
	    vu[-1] = vu[0];
	}

#if STEADY_STATE || TIME_DEP_NON
	vu[1] = phgQuadGetDofValues(e, ns->u[1], quad);          /* u^{n+1} */
	gu[1] = phgQuadGetDofValues(e, ns->gradu[1], quad);      /* grad u^{n} */
	vp[1] = phgQuadGetDofValues(e, ns->p[1], quad);          /* p^{n+1} */
#else
	TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */

	Unused(l);
	Bzero(vf); Bzero(vf_cache); 

	if (_nsp->extern_force) {
	    /* cache f values,
	     * only time_{n} */
	    for (l = 1; l < 2; l++) {
		const FLOAT *cache;
		size_t cache_size;
		setFuncTime(ns->time[l]); /* set static time in ins-test.c */

		/* cache f */
		cache_size = Dim * quad->npoints * sizeof(FLOAT);
		cache = phgQuadGetFuncValues(g, e, Dim, func_f, quad);
		vf[l] = vf_cache[l] = phgAlloc(cache_size);
		memcpy(vf_cache[l], cache, cache_size);

		phgQuadGetFuncValues(NULL, NULL, 0, NULL, NULL); /* clear cache */
	    }
	}

	/* Global Matrix */
	Bzero(rhsu); Bzero(rhsp);
    
	/* average pressure */
	vp0 = 0.;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++, w++)
	    vp0 += (*w) * vp[1][q];

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);

	    /* rhs u */
	    for (i = 0; i < M; i++) {
		/* interior node or Neumann */
		const FLOAT *gi_u = phgQuadGetBasisValues(e, ns->du, i, quad) + q;       /* phi_i */
		const FLOAT *ggi_u = phgQuadGetBasisCurvedGradient(e, ns->du, i, quad, q);    /* grad phi_i */

		for (k = 0; k < Dim; k++) {
		    nu = get_effective_viscosity(gu[1], *vTe, 0, viscosity_type);
		    FLOAT eu[DDim];

		    MAT3_SYM(gu[1], eu);
		    rhsu[i][k] += vol*(*w) * (- nu * INNER_PRODUCT(eu+k*Dim, ggi_u) 
					      + (*vp[1]) * *(ggi_u+k) * LEN_SCALING * PRES_SCALING
					      ) * EQU_SCALING;     /* left */

		    if (k == Z_DIR) { 
			const FLOAT rho = RHO_ICE;
			const FLOAT grav = GRAVITY;
			const FLOAT a = SEC_PER_YEAR;
			const FLOAT f = rho*grav * LEN_SCALING2 * EQU_SCALING; 

			Unused(a);
			rhsu[i][k] += vol*(*w) * (-f * (*gi_u) 
						  ); /* right */
		    }

		    if (_nsp->compensate_equ)
			rhsu[i][k] += vol*(*w) * EQU_SCALING * (- vf[1][k] * (*gi_u) * LEN_SCALING2);

		}
	    }

	    /* rhs p */
	    for (i = 0; i < N; i++) 
	    {
		const FLOAT *gi_p = phgQuadGetBasisValues(e, ns->dp, i, quad) + q;       /* psi_i */
		FLOAT divu1 = gu[1][0] + gu[1][4] + gu[1][8];
		//FLOAT divu0 = gu[0][0] + gu[0][4] + gu[0][8];
		rhsp[i] += vol*(*w) * (divu1 * (*gi_p));

		if (ns->du->type == DOF_P1
		    && ns->dp->type == DOF_P1) {
		    rhsp[i] -= SIGN_STAB * vol*(*w) *
			( ((*vp[1]) - vp0) * ((*gi_p) -  1./4.) ) ;
		}
	    }
	    
	    if (tstep_minus) 
		vu[-1] += Dim;

	    vu[1] += Dim;
	    gu[1] += Dim*Dim;
	    vp[1]++;
	    vu[0] += Dim;
	    gu[0] += Dim * Dim;
	    vp[0]++; 
	    //vw += Dim;
	    if (_nsp->extern_force) {
		vf[0] += Dim; vf[1] += Dim;
	    }
	    vTe++;
	    w++; p += Dim + 1;
	}

	if (_nsp->extern_force) {
	    phgFree(vf_cache[0]);
	    phgFree(vf_cache[1]);
	}

	normal = NULL; Unused(normal);
	area = 0; Unused(area);


	if (!_nsp->enclosed_flow) {
	    /* slip boundary */
	    for (s = 0; s < NFace; s++) {
		if (e->bound_type[s] == (BC_BOTTOM|BC_BOTTOM_GRD))//((e->bound_type[s] & BC_BOTTOM_GRD) && !(e->bound_type[s] & BC_ISHELF)) 
        {
		    int v0, v1, v2;
		    int nbas_face = NbasFace(ns->du);
		    SHORT bases[nbas_face];
		    FLOAT lambda[Dim + 1], x,y,z, beta;
		    order = DofTypeOrder(ns->du, e) * 3 - 1;

		    phgDofGetBasesOnFace(ns->du, e, s, bases);
		    v0 = GetFaceVertex(s, 0);
		    v1 = GetFaceVertex(s, 1);
		    v2 = GetFaceVertex(s, 2);
		    lambda[s] = 0.;
		    
		    area = phgGeomGetFaceArea(g, e, s);
		    normal = phgGeomGetFaceOutNormal(g, e, s);
		    quad = phgQuadGetQuad2D(order);
		
		    p = quad->points;
		    w = quad->weights;
		    for (q = 0; q < quad->npoints; q++) {
			FLOAT vu[Dim];
            FLOAT ub_mag, surf_elev, bot_elev, thickness, effective_pressure;
            FLOAT alpha2, beta2, mindex;

			lambda[v0] = *(p++);
			lambda[v1] = *(p++);
			lambda[v2] = *(p++);
			
			phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
			func_beta(x, y, z, &beta);
			phgDofEval(ns->u[1], e, lambda, vu);
            phgDofEval(ns->surf_elev_P1, e, lambda, &surf_elev);
            phgDofEval(ns->bot_elev_P1, e, lambda, &bot_elev);

            thickness = surf_elev - bot_elev;

            effective_pressure = RHO_ICE*GRAVITY*(thickness - MAX(0, -RHO_WAT/RHO_ICE*bot_elev));
            
            ub_mag = pow(INNER_PRODUCT(vu,vu)+1e-10,0.5);

            alpha2 = _nsp->slip_alpha2;
            beta2 = _nsp->slip_beta2;
            mindex = _nsp->slip_index;



            FLOAT C_mtp1;

           //FLOAT  C_mtp1 = C_mtp*(1 - 0.75*exp(-pow(x-522296, 2)/(2.0*225e8) - y*y/(2.0*1e8)));
           if (Time < 10000)
               C_mtp1 = C_mtp;
           else
               C_mtp1 = (_nsp->slip_beta2)*(1 - 0.75*exp(-pow(x-528000, 2)/(2.0*225e8) - y*y/(2.0*1e8)));


			const FLOAT *gi_u = 
			    ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);

			for (i = 0; i < nbas_face; i++) {
			    int ii = bases[i];
			    for (k = 0; k < Dim; k++) {
				if (_nsp->slip_condition == 0)
				    rhsu[ii][k] += -area*(*w) * beta * vu[k] * (gi_u[ii]) * EQU_SCALING * LEN_SCALING;
                else if (_nsp->slip_condition == 1)
				    rhsu[ii][k] += -area*(*w) * C_mtp1 * pow(INNER_PRODUCT(vu, vu)+1e-10, (m_idx-1)/2.0) * vu[k] * gi_u[ii] * EQU_SCALING * LEN_SCALING;
                else if (_nsp->slip_condition == 2)
				    rhsu[ii][k] += -area*(*w) * beta2 * pow(ub_mag, 1.0/mindex) * alpha2 * effective_pressure/pow(pow(beta2, mindex) * ub_mag + pow(alpha2 * effective_pressure, mindex),1.0/mindex) / ub_mag * vu[k] * gi_u[ii] * EQU_SCALING * LEN_SCALING;

                else
                    phgPrintf("Wrong sliding law ! slip_condition should be 0, 1 or 2");
			    }
			}     /* end of bas_i */
			w++;
		    }		/* end of quad point */
		}		/* end of face outflow */
	    }			/* end of all outflow face in element */

	    if (_nsp->compensate_equ) {
		/* Compensate surface boundary condition */
#if 0
		for (s = 0; s < NFace; s++) {
		    if (e->bound_type[s] & BC_TOP) { /* Note: only top surf is Neumann */
			int v0, v1, v2;
			int nbas_face = NbasFace(ns->u[1]);
			SHORT bases[nbas_face];
			FLOAT lambda[Dim + 1], x,y,z;
			order = DofTypeOrder(ns->u[1], e) * 3 - 1;

			phgDofGetBasesOnFace(ns->u[1], e, s, bases);
			v0 = GetFaceVertex(s, 0);
			v1 = GetFaceVertex(s, 1);
			v2 = GetFaceVertex(s, 2);
			lambda[s] = 0.;
		    
			area = phgGeomGetFaceArea(g, e, s);
			normal = phgGeomGetFaceOutNormal(g, e, s);
			quad = phgQuadGetQuad2D(order);
		
			p = quad->points;
			w = quad->weights;
			for (q = 0; q < quad->npoints; q++) {
			    FLOAT gn[Dim];
			    lambda[v0] = *(p++);
			    lambda[v1] = *(p++);
			    lambda[v2] = *(p++);
			
			    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
			    //if (e->bound_type[s] & BC_TOP)
			    func_g1(x, y, z, gn); /* Surf */

			    for (i = 0; i < nbas_face; i++) {
				int ii = bases[i];
				FLOAT gi_u = 
				    *ns->u[1]->type->BasFuncs(ns->u[1], e, ii, ii + 1, lambda);

				for (k = 0; k < Dim; k++) {
				    rhsu[ii][k] += area*(*w) * gn[k] * (gi_u) * LEN_SCALING
					* EQU_SCALING;
				}
			    }     /* end of bas_i */
			    w++;
			} /* end of quad point */
		    }	  /* end of face outflow */
		}	  /* end of all outflow face in element */
#endif
	    }		  /* end of compensate equations */

	                       /* end out flow boundary */

        if (_nsp->lateral_sia_pres_bc){
        for (s = 0; s < NFace; s++) {

        if ((e->bound_type[s] & BC_LATERL) && 
           !(e->bound_type[s] & BC_ISHELF))
       { 
        int v0, v1, v2;
        int nbas_face = NbasFace(ns->u[1]);
        SHORT bases[nbas_face];
        FLOAT lambda[Dim + 1], x,y,z;
        order = DofTypeOrder(ns->u[1], e) * 3 - 1;
        /* DofTypeOrder get the u type order, i.e. 2 */

        phgDofGetBasesOnFace(ns->u[1], e, s, bases);
        v0 = GetFaceVertex(s, 0);
        v1 = GetFaceVertex(s, 1);
        v2 = GetFaceVertex(s, 2);
        /* v0, v1, v2 are the vertex no in the element */
        lambda[s] = 0.;
        
        area = phgGeomGetFaceArea(g, e, s);
        normal = phgGeomGetFaceOutNormal(g, e, s);
        quad = phgQuadGetQuad2D(order);
        
        p = quad->points;
        w = quad->weights;

        for (q = 0; q < quad->npoints; q++) {
            FLOAT surf_elev;
            FLOAT grad_surf_elev[Dim];

            lambda[v0] = *(p++);
            lambda[v1] = *(p++);
            lambda[v2] = *(p++);
        
            phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);

            phgDofEval(ns->surf_elev_P1, e, lambda, &surf_elev);
            //phgDofEval(ns->grad_surf_elev, e, lambda, &grad_surf_elev[0]);

            FLOAT dsdx, dsdy;

            dsdx = 0;
            dsdy = 0;
            //func_sur_grad_x(x, y, z, &dsdx);
            //func_sur_grad_y(x, y, z, &dsdy);

            for (i = 0; i < nbas_face; i++) {
            int ii = bases[i];
            FLOAT gi_u = 
                *ns->u[1]->type->BasFuncs(ns->u[1], e, ii, ii + 1, lambda);

            for (k = 0; k < Dim; k++) {
            if (k < Dim-1) {

                rhsu[ii][k] += -area* (*w) * RHO_ICE * GRAVITY * (surf_elev - z)
                               * normal[k] * (gi_u) * LEN_SCALING * EQU_SCALING;
            }

            else {
              rhsu[ii][k] += -area*(*w) * RHO_ICE*GRAVITY* (surf_elev - z)
                           *(normal[k-2]*(dsdx)+ normal[k-1]*(dsdy))
                           * (gi_u) * LEN_SCALING*EQU_SCALING;

            }
            }
            }     /* end of bas_i */
            w++;
        } /* end of quad point */
        }
		}
        }

        if (_nsp->add_ice_shelf) {
		for (s = 0; s < NFace; s++) {
		    if ((e->bound_type[s] & BC_ISHELF)) /* && !(e->bound_type[s] & BC_LATERL)) */
            {
		
            FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
            INT quad_order = 5;
            FLOAT Ns;
            FLOAT vu[Dim], avg_n_v[Dim];
            QUAD *quad = phgQuadGetQuad2D(quad_order);
            FLOAT *p = quad->points;
            FLOAT *w = quad->weights;
            FLOAT area = phgGeomGetFaceArea(g, e, s);
            //FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);
            FLOAT normal[Dim];

            FLOAT nx, ny, nz;

            INT M_face = 3*(ns->du->type->np_vert + ns->du->type->np_edge);

            SHORT bas_idx_e[M_face];

            phgDofGetBasesOnFace(ns->du, e, s, bas_idx_e);
            
            INT v0 = GetFaceVertex(s, 0);
            INT v1 = GetFaceVertex(s, 1);
            INT v2 = GetFaceVertex(s, 2);
            lambda[s] = 0;

            if (fabs(normal[2]) < 1.0e-8)
            {
                Ns = 1.0e50;
            }
            else
                Ns = sqrt(1+SQUARE(normal[0]/normal[2])+SQUARE(normal[1]/normal[2]));

            FLOAT dt1 = ns->dt[0];

            for (q = 0; q < quad->npoints; q++)
            {
                lambda[v0] = *(p++);
                lambda[v1] = *(p++);
                lambda[v2] = *(p++);

                phgGeomLambda2XYZ(g, e, lambda, &lambda_x, &lambda_y, &lambda_z);
                phgDofEval(ns->u[1], e, lambda, vu);
#if 1
                phgDofEval(avg_n, e, lambda, avg_n_v);
                normal[0] = avg_n_v[0];
                normal[1] = avg_n_v[1];
                normal[2] = avg_n_v[2];
            if (fabs(normal[2]) < 1.0e-8)
            {
                Ns = 1.0e50;
            }
            else
                Ns = sqrt(1+SQUARE(normal[0]/normal[2])+SQUARE(normal[1]/normal[2]));

#endif
                
                const FLOAT *bas = ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);

                for (i = 0; i < M_face; i++)
                {
                    //if (phgDofGetElementBoundaryType(ns->du, e, i) & BC_BOTTOM_GRD)
                    //    continue;

                    INT i_e = bas_idx_e[i];
                    const FLOAT *bas_i = phgQuadGetBasisValues(e, ns->du, i_e, quad);
                    for (k = 0; k < Dim; k++)
                    {
                        if (lambda_z < 0)
                        {
                            rhsu[i_e][k] += area*w[q]*RHO_WAT*GRAVITY*
                                (lambda_z-(INNER_PRODUCT(vu, normal)-0)*Ns*dt1)*bas[i_e]
                                *normal[k]*EQU_SCALING;
                        }
                        else
                            rhsu[i_e][k] += 0;
                    }
                }
            } /* end of quad point */
		    }	  /* end of face outflow */
		}	  /* end of all outflow face in element */

		for (s = 0; s < NFace; s++) {
		    if (e->bound_type[s] & BC_TERMNS) {
		    //if (e->bound_type[s] == (BC_LATERL | BC_ISHELF)){
		
            FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
            INT quad_order = 5;
            FLOAT Ns;
            FLOAT vu[Dim];
            QUAD *quad = phgQuadGetQuad2D(quad_order);
            FLOAT *p = quad->points;
            FLOAT *w = quad->weights;
            FLOAT area = phgGeomGetFaceArea(g, e, s);
            const FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

            FLOAT nx1, ny1, nz1;

            INT M_face = 3*(ns->du->type->np_vert + ns->du->type->np_edge);
            SHORT bas_idx_e[M_face];

            phgDofGetBasesOnFace(ns->du, e, s, bas_idx_e);
            
            INT v0 = GetFaceVertex(s, 0);
            INT v1 = GetFaceVertex(s, 1);
            INT v2 = GetFaceVertex(s, 2);
            lambda[s] = 0;


            FLOAT dt1 = ns->dt[0];

            for (q = 0; q < quad->npoints; q++)
            {
                lambda[v0] = *(p++);
                lambda[v1] = *(p++);
                lambda[v2] = *(p++);

                phgGeomLambda2XYZ(g, e, lambda, &lambda_x, &lambda_y, &lambda_z);
                phgDofEval(ns->u[1], e, lambda, vu);
                
                const FLOAT *bas = ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);

                for (i = 0; i < M_face; i++)
                {
                    INT i_e = bas_idx_e[i];
                    const FLOAT *bas_i = phgQuadGetBasisValues(e, ns->du, i_e, quad);
                    for (k = 0; k < Dim; k++)
                    {
                        if (lambda_z < 0)
                        {
                            rhsu[i_e][k] += 2      *        area*w[q]*RHO_WAT*GRAVITY*
                                (lambda_z)*bas[i_e]
                                *normal[k]*EQU_SCALING;
                            /* I add a multiplier factor 2 for testing the melange impact */
                        }
                        else
                            rhsu[i_e][k] += 0;
                    }
                }
            } /* end of quad point */
		    }	  /* end of face outflow */
		}
	    }		  /* end of compensate equations */
	    
	}           


#if USE_SLIDING_BC || CASE_DRAINAGE_BASIN
	/* Rotate bases */
	for (i = 0; i < M; i++) {
	    INT id = phgDofMapE2D(surf_dof, e, i * (Dim*Dim)) / (Dim*Dim);
	    if (!rotated[id])
		continue;	
	    const FLOAT *trans = Trans + id*(Dim*Dim);

	    trans_left(&rhsu[i][0], 1, 1, trans);
	}
#endif

#if 1 && USE_NODAL_LOADS
    memcpy(rhsu0, rhsu, sizeof(rhsu));
    memcpy(rhsp0, rhsp, sizeof(rhsp));
#endif
    // save the rhsu0 and rhsp0 for contact problem

	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    for (k = 0; k < Dim; k++)
		Iu[i][k] = phgMapE2L(solver_u->rhs->map, 0, e, i * Dim + k);
	for (i = 0; i < N; i++)
	    Ip[i] = phgMapE2L(solver_u->rhs->map, 1, e, i);

    	/* set velocity dirichlet bdry */
	FLOAT tmp[Dim];
	for (i = 0; i < M; i++)
	    for (k = 0; k < Dim; k++) {
        if (IF_DB)
		if (phgDofDirichletBC_(ns->du, e, i*Dim+k, NULL, bufu, tmp,
				       DOF_PROJ_NONE)) {
            FLOAT *coord_e=phgDofGetElementCoordinates(ns->du, e, i*Dim+k);
            BTYPE btype_e=phgDofGetElementBoundaryType(ns->du, e, i*Dim+k);
            FLOAT *data_e=DofElementData(ns->u[1], e->index);

            
		    rhsu[i][k] = 0.;
		}
	    }


#if STEADY_STATE || TIME_DEP_NON
	/* set pressure dirichlet bdry for pinned point */
	for (i = 0; i < N; i++)
        if (IF_DB)
	    if (phgDofDirichletBC(ns->dp, e, i, NULL, bufp, &rhsp[i],
				  DOF_PROJ_NONE)) {
		if (!_nsp->enclosed_flow)
		    phgError(1, "No dirichlet bc for Unenclosed flow!\n");
		if (_nsp->pin_node) {
# if PIN_AT_ROOT 
		    if (g->rank != 0)
		    	phgError(1, "Pinned node only on rank 0!\n");
		    if (g, e->verts[i] != ns->pinned_node_id)
			phgError(1, "Build rhs: pinned node e:%d, bas:%d, [%d] and [%d] "
				 "doesn't coincide when build RHS!\n", 
				 e->index, i, e->verts[i], ns->pinned_node_id);
# else
		    if (GlobalVertex(g, e->verts[i]) != ns->pinned_node_id)
			phgError(1, "Build rhs: pinned node e:%d, bas:%d, [%d] and [%d] "
				 "doesn't coincide when build RHS!\n", 
				 e->index, i, e->verts[i], ns->pinned_node_id);
# endif /* PIN_AT_ROOT */
		}
	    } 
#else
	TIME_DEP_LINEAR; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */


#if 0 && USE_NODAL_LOADS
    memcpy(rhsu0, rhsu, sizeof(rhsu));
    memcpy(rhsp0, rhsp, sizeof(rhsp));
#endif

	/* Global res */
	phgVecAddEntries(vec_rhs, 0, M * Dim, Iu[0], &rhsu[0][0]);
	phgVecAddEntries(vec_rhs, 0, N, Ip, rhsp);

#if USE_NODAL_LOADS
	phgVecAddEntries(ns->vec_rhs0, 0, M * Dim, Iu[0], &rhsu0[0][0]);
	phgVecAddEntries(ns->vec_rhs0, 0, N, Ip, rhsp0);
#endif
    }				/* end element */

	phgDofFree(&avg_n);
    
    phgVecAssemble(vec_rhs);

#if USE_NODAL_LOADS
    phgVecAssemble(ns->vec_rhs0);
#endif

    phgVecAssemble(solver_u->rhs);
    phgVecAXPBY(1., vec_rhs, 0, &solver_u->rhs);
    solver_u->rhs_updated = FALSE;

    if (DUMP_MAT_VEC) 
    {
	phgPrintf("Dumping rhs\n");
	//phgVecDumpMATLAB(ns->solver_u->rhs, "b", "b_.m");
    }

    phgVecDestroy(&vec_rhs);

    {
	FLOAT a[2] = {nu_max, -nu_min}, b[2]; 
	MPI_Allreduce(&a, &b, 2, PHG_MPI_FLOAT, MPI_MAX, g->comm);
	nu_max = b[0];
	nu_min = -b[1];
	phgPrintf("  vis: [%8.4e, %8.4e] ", nu_min, nu_max);
    }

    return;
}


#define eu_xx eu[0]
#define eu_xy eu[1]
#define eu_xz eu[2]
#define eu_yx eu[3]
#define eu_yy eu[4]
#define eu_yz eu[5]
#define eu_zx eu[6]
#define eu_zy eu[7]
#define eu_zz eu[8]

static FLOAT * 
get_gbas_product(const FLOAT *gi, const FLOAT *gj,
		 const FLOAT *gu, LTYPE ltype) 
{
    static FLOAT prod[Dim][Dim];
    FLOAT Gi[Dim], Gj[Dim];
    FLOAT eu[DDim];
    FLOAT eps, eps2;
    int k;

    /* Picard term */
    prod[0][0] = gi[0] * gj[0] + .5 * (gi[1] * gj[1] + gi[2] * gj[2]);
    prod[0][1] = .5 * gi[0] * gj[1];
    prod[0][2] = .5 * gi[0] * gj[2];

    prod[1][0] = .5 * gi[1] * gj[0]; 
    prod[1][1] = gi[1] * gj[1] + .5 * (gi[0] * gj[0] + gi[2] * gj[2]);
    prod[1][2] = .5 * gi[1] * gj[2];

    prod[2][0] = .5 * gi[2] * gj[0]; 
    prod[2][1] = .5 * gi[2] * gj[1]; 
    prod[2][2] = gi[2] * gj[2] + .5 * (gi[0] * gj[0] + gi[1] * gj[1]);

    if (ltype == PICARD) {
	return prod[0];
    }

    /* Newton term */
    MAT3_SYM(gu, eu);
    for (k = 0; k < DDim; k++)
	eu[k] /= LEN_SCALING;
    eps = sqrt(.5) * MAT3_NORM2(eu);
    
    if (eps < MIN_EFFECTIVE_STRAIN) 
	eps = MIN_EFFECTIVE_STRAIN;

    eps2 = - (1./3.) / (eps*eps);

    Gi[0] = eu_xx * gi[0] + eu_xy * gi[1] + eu_xz * gi[2];
    Gi[1] = eu_yx * gi[0] + eu_yy * gi[1] + eu_yz * gi[2];
    Gi[2] = eu_zx * gi[0] + eu_zy * gi[1] + eu_zz * gi[2];

    Gj[0] = eu_xx * gj[0] + eu_xy * gj[1] + eu_xz * gj[2];
    Gj[1] = eu_yx * gj[0] + eu_yy * gj[1] + eu_yz * gj[2];
    Gj[2] = eu_zx * gj[0] + eu_zy * gj[1] + eu_zz * gj[2];
    
    prod[0][0] += Gi[0] * Gj[0] * eps2;
    prod[0][1] += Gi[1] * Gj[0] * eps2;
    prod[0][2] += Gi[2] * Gj[0] * eps2;

    prod[1][0] += Gi[0] * Gj[1] * eps2;
    prod[1][1] += Gi[1] * Gj[1] * eps2;
    prod[1][2] += Gi[2] * Gj[1] * eps2;

    prod[2][0] += Gi[0] * Gj[2] * eps2;
    prod[2][1] += Gi[1] * Gj[2] * eps2;
    prod[2][2] += Gi[2] * Gj[2] * eps2;

    return prod[0];
}

#undef eu_xx
#undef eu_xy
#undef eu_xz
#undef eu_yx
#undef eu_yy
#undef eu_yz
#undef eu_zx
#undef eu_zy
#undef eu_zz



/***************************************************************************
 * Build matrices *F, *Fu, *B, *Bt, *Fp, *Ap, and *Qp used * by the .
 **************************************************************************/
void
phgNSBuildSolverUMat(NSSolver *ns, INT IF_DB, INT nonstep, FLOAT Time)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, k, l, q, s;
    FLOAT *dt = ns->dt;
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    BOOLEAN use_Fu = _nsp->use_Fu;
    int viscosity_type = ns->viscosity_type;
    LTYPE ltype = ns->ltype;

#if USE_SLIDING_BC
    SURF_BAS *surf_bas = ns->surf_bas;
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);
    assert(surf_dof->type == ns->du->type);
#endif


    l = 0; Unused(l); 
    Unused(nu);
#if STEADY_STATE
    assert(fabs(Theta - 1) < 1e-12);
    Thet1 = 0; Unused(Thet1);
    Unused(dt);
#else
    Thet1 = 1 - Theta;
    Unused(dt);
#endif /* STEADY_STATE */

    if (ltype == PICARD)
	phgPrintf("   LinearType: Picard");
    else 
	phgPrintf("   LinearType: Newton");

	DOF *avg_n = phgDofNew(ns->g, DOF_P2, 3, "avg_n", DofNoAction);

    get_avg_n(ns->g, avg_n);

    ForAllElements(g, e) {


	int M = ns->du->type->nbas;	/* num of bases of Velocity */
	int N = ns->dp->type->nbas;	/* num of bases of Pressure */
	int order = DofTypeOrder(ns->du, e) * 3 - 1; /* Note:
							*   quad order is really high here,
							*   highest order term (u \nabla u, phi)  */
	FLOAT F[M][Dim][M][Dim], Fu[M][M],
	    B[N][M][Dim], Bt[M][Dim][N], C[N][N],
	    bufu[M], bufp[N],
	    rhsu[M][Dim];
	INT Iu[M][Dim], Iu1[M], Ju[Dim][M], Ip[N];
	QUAD *quad;
	FLOAT vol, det;
	const FLOAT *w, *p, *vw, *gu, *vTe;

	Unused(Iu1);
	for (i = 0; i < M; i++)
	    for (k = 0; k < Dim; k++)
		Ju[k][i] = Iu[i][k] = phgMapE2L(ns->matF->cmap, 0, e, i * Dim + k);
	if (use_Fu)
	    for (i = 0; i < M; i++)
		Iu1[i] = phgMapE2L(_pcd->matFu->cmap, 0, e, i);
	for (i = 0; i < N; i++)
	    Ip[i] = phgMapE2L(ns->matC->cmap, 0, e, i);


	quad = phgQuadGetQuad3D(order);
	//vw = phgQuadGetDofValues(e, ns->wind, quad);      /* value wind */
	gu = phgQuadGetDofValues(e, ns->gradu[1], quad);  /* grad u */
	if (ns_params->noniter_temp)
	    vTe = phgQuadGetDofValues(e, ns->T[1], quad);	  /* T^{n+1} */
	else
	    vTe = phgQuadGetDofValues(e, ns->T[0], quad);	  /* T^{n} */

	Bzero(F); Bzero(Fu); 
	Bzero(B); Bzero(Bt); Bzero(C);

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);
	    /* Test func vel type */
	    for (i = 0; i < M; i++) {
		const FLOAT *gi_u = phgQuadGetBasisValues(e, ns->du, i, quad) + q;       /* phi_i */
		const FLOAT *ggi_u = phgQuadGetBasisCurvedGradient(e, ns->du, i, quad, q);    /* grad phi_i */

		/* Mat F */
		for (j = 0; j < M; j++) {
		    const FLOAT *gj_u = phgQuadGetBasisValues(e, ns->du, j, quad) + q;       /* phi_j */
		    const FLOAT *ggj_u = phgQuadGetBasisCurvedGradient(e, ns->du, j, quad, q);    /* grad phi_i */
		    FLOAT mass = (*gj_u) * (*gi_u);
		    FLOAT diffu = INNER_PRODUCT(ggj_u, ggi_u);

		    Unused(mass);
		    nu = get_effective_viscosity(gu, *vTe, 0, viscosity_type);

		    const FLOAT *tp = get_gbas_product(ggi_u, ggj_u, gu, ltype);

		    for (k = 0; k < Dim; k++) 
			for (l = 0; l < Dim; l++) 
			    F[j][l][i][k] += vol*(*w) * EQU_SCALING * nu * tp[k+l*Dim];

		    if (use_Fu) 
			Fu[i][j] += vol*(*w) * EQU_SCALING * (nu * diffu);

		}

		/* Mat B & Bt */
		for (j = 0; j < N; j++) {
		    const FLOAT *gj_p = phgQuadGetBasisValues(e, ns->dp, j, quad) + q;       /* psi_j */
		    for (k = 0; k < Dim; k++) {
			FLOAT b = vol*(*w) * (*gj_p) * ggi_u[k];
			B[j][i][k]  -= b; /* divergence */
			Bt[i][k][j] -= EQU_SCALING * LEN_SCALING * PRES_SCALING * b; /* momentum */
		    }
		}
	    } /* end phi_i */

	    /* Test func pressure type
	     * Some solver doesn't allow zero diagnal entries in Matrix, in that case,
	     *   we perturb diagnal by a small number.
	     * */
	    for (i = 0; i < N; i++) {
		const FLOAT *gi_p = phgQuadGetBasisValues(e, ns->dp, i, quad) + q;       /* psi_i */
		C[i][i] += vol*(*w) * (_nsp->eps_diagP * (*gi_p) * (*gi_p)
				       );
	    }

	    /*
	     * P1/P1 pressure local projection stabilization.
	     * See Elman P243.
	     *
	     *  */
	    if (ns->du->type == DOF_P1
		&& ns->dp->type == DOF_P1) {
		for (i = 0; i < N; i++) {
		    const FLOAT *gi_p = phgQuadGetBasisValues(e, ns->dp, i, quad) + q;       /* psi_i */
		    for (j = 0; j < N; j++) {
			const FLOAT *gj_p = phgQuadGetBasisValues(e, ns->dp, j, quad) + q;       /* psi_j */

			C[i][j] += SIGN_STAB * vol*(*w) * ((*gj_p) * (*gi_p)
							   - (1./16.)
							   );
		    }
		}
	    }


	    /* Next quad point */
	    //vw += Dim;
	    gu += Dim*Dim;
	    vTe++;
	    w++; p += Dim+1;
	}

	if (!_nsp->enclosed_flow) {
	    /* slip boundary */
	    for (s = 0; s < NFace; s++) {
		if (e->bound_type[s] == (BC_BOTTOM|BC_BOTTOM_GRD))//((e->bound_type[s] & BC_BOTTOM_GRD) && !(e->bound_type[s] & BC_ISHELF)) 
        {
		    int v0, v1, v2;
		    int nbas_face = NbasFace(ns->du);
		    SHORT bases[nbas_face];
		    FLOAT lambda[Dim + 1], area, x, y, z, beta;
		    order = DofTypeOrder(ns->du, e) * 3 - 1;

		    phgDofGetBasesOnFace(ns->du, e, s, bases);
		    v0 = GetFaceVertex(s, 0);
		    v1 = GetFaceVertex(s, 1);
		    v2 = GetFaceVertex(s, 2);
		    lambda[s] = 0.;

		    area = phgGeomGetFaceArea(g, e, s);
		    quad = phgQuadGetQuad2D(order);

		    p = quad->points;
		    w = quad->weights;
		    for (q = 0; q < quad->npoints; q++) {
			FLOAT vu[Dim];
			lambda[v0] = *(p++);
			lambda[v1] = *(p++);
			lambda[v2] = *(p++);
            FLOAT ub_mag, surf_elev, bot_elev, thickness, effective_pressure;
            FLOAT alpha2, beta2, mindex;

            //FLOAT depth;
			
			phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
			func_beta(x, y, z, &beta);
			phgDofEval(ns->u[1], e, lambda, vu);


            phgDofEval(ns->surf_elev_P1, e, lambda, &surf_elev);
            phgDofEval(ns->bot_elev_P1, e, lambda, &bot_elev);
            //phgDofEval(ns->depth_P1, e, lambda, &depth);

            thickness = surf_elev - bot_elev;

            //phgPrintf("depth, thickness %f %f\n", depth, thickness);

            effective_pressure = RHO_ICE*GRAVITY*(thickness - MAX(0, -RHO_WAT/RHO_ICE*bot_elev));
            
            ub_mag = pow(INNER_PRODUCT(vu,vu)+1e-10,0.5);

            alpha2 = _nsp->slip_alpha2;
            beta2 = _nsp->slip_beta2;
            mindex = _nsp->slip_index;
			
            
           //FLOAT  C_mtp1 = C_mtp*(1 - 0.75*exp(-pow(x-522296, 2)/(2.0*225e8) - y*y/(2.0*1e8)));

            FLOAT C_mtp1;

            if ( Time < 10000)
                C_mtp1 = C_mtp;
            else
                C_mtp1 = C_mtp*(1 - 0.75*exp(-pow(x-528000, 2)/(2.0*225e8) - y*y/(2.0*1e8)));

			const FLOAT *gi_u = 
			    ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);
			for (i = 0; i < nbas_face; i++) {
				FLOAT mass_face;
			    int ii = bases[i];
			    for (j = 0; j < nbas_face; j++) { 
				int jj = bases[j];
                if (_nsp->slip_condition == 0)
				    mass_face = area*(*w) * beta * (gi_u[jj])*(gi_u[ii])* EQU_SCALING * LEN_SCALING;
                else if (_nsp->slip_condition == 1)
				    mass_face = area*(*w) * C_mtp1 * 
                        pow(INNER_PRODUCT(vu, vu)+1e-10, (m_idx-1)/2.0) * 
                        gi_u[jj] * gi_u[ii] * EQU_SCALING * LEN_SCALING;
                else if (_nsp->slip_condition == 2)
				    mass_face = area*(*w) * beta2 * pow(ub_mag, 1.0/mindex) * alpha2 * 
                        effective_pressure/pow(pow(beta2, mindex) * ub_mag + 
                                pow(alpha2 * effective_pressure, mindex),1.0/mindex)/ub_mag * 
                        gi_u[jj] * gi_u[ii] * EQU_SCALING * LEN_SCALING;
                else
                    phgPrintf("Wrong sliding law ! The slip_condition should be 0, 1 or 2");

				for (k = 0; k < Dim; k++) 
				{
				    F[ii][k][jj][k] -= SIGN_FRICTION * mass_face;
				}
			    } /* end of bas_j */
			}     /* end of bas_i */
			w++;
		    }		/* end of quad point */
		}		/* end of face outflow */
	    }			/* end of all outflow face in element */


    for (s = 0; s < NFace; s++)
        {
            if ((e->bound_type[s] & BC_ISHELF)) /* && !(e->bound_type[s] & BC_LATERL))*/
            {
                FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
                INT quad_order = 5;
                FLOAT Ns;
                QUAD *quad = phgQuadGetQuad2D(quad_order);
                FLOAT *w = quad->weights;
                FLOAT *p = quad->points;

                FLOAT avg_n_v[Dim];

                FLOAT nx, ny, nz;

                FLOAT area = phgGeomGetFaceArea(g, e, s);
                //FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);
                FLOAT normal[Dim];

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0;

                if (fabs(normal[2]) < 1.0e-8)
                {
                    Ns = 1.0e50;
                }
                else
                    Ns = sqrt(1+SQUARE(normal[0]/normal[2])+SQUARE(normal[1]/normal[2]));

                FLOAT dt1 = ns->dt[0];

                INT M_face = 3*(ns->du->type->np_vert + ns->du->type->np_edge);
                SHORT bas_idx_e[M_face];

                phgDofGetBasesOnFace(ns->du, e, s, bas_idx_e);

                for (q = 0; q < quad->npoints; q++)
                {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);
                    
                    phgGeomLambda2XYZ(g, e, lambda, &lambda_x, &lambda_y, &lambda_z);

                    phgDofEval(avg_n, e, lambda, avg_n_v);
                    normal[0] = avg_n_v[0];
                    normal[1] = avg_n_v[1];
                    normal[2] = avg_n_v[2];
                if (fabs(normal[2]) < 1.0e-8)
                {
                    Ns = 1.0e50;
                }
                else
                    Ns = sqrt(1+SQUARE(normal[0]/normal[2])+SQUARE(normal[1]/normal[2]));
                    
                    const FLOAT *bas = ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);
                    
                    for (i = 0; i < M_face; i++)
                    {
                        //if (phgDofGetElementBoundaryType(ns->du, e, i) & BC_BOTTOM_GRD)
                        //    continue;

                        INT i_e = bas_idx_e[i];
                        const FLOAT *bas_i = phgQuadGetBasisValues(e, ns->du, i_e, quad);
                        for (j = 0; j < M_face; j++)
                        {
                            //if (phgDofGetElementBoundaryType(ns->du, e, j) & BC_BOTTOM_GRD)
                            //    continue;

                            INT j_e = bas_idx_e[j];
                            const FLOAT *bas_j = phgQuadGetBasisValues(e, ns->du, j_e, quad);
                            
                            const FLOAT *bas_dot_n = get_bas_dot_normal(bas, normal, i_e, j_e);
                            for (k = 0; k < Dim; k++)
                            {
                                    //F[j_e][k][i_e][k] += area*w[q]*RHO_WAT*GRAVITY*bas[i_e]*bas[j_e]*normal[k]*normal[k]*Ns*dt1*EQU_SCALING;

#if 1
                                for (l = 0; l < Dim; l++)
                                {
                                    F[j_e][l][i_e][k] += area*w[q]*RHO_WAT*GRAVITY*
                                        bas_dot_n[k+l*Dim]*Ns*dt1*EQU_SCALING;
                                }
#endif

                            }
                        }
                    }
                }
            }
        }
	}		  /* end out flow boundary */

#if USE_SLIDING_BC
	/* Rotate bases */
	for (i = 0; i < M; i++) {
	    INT id = phgDofMapE2D(surf_dof, e, i * (Dim*Dim)) / (Dim*Dim);
	    if (!rotated[id])
		continue;	
	    const FLOAT *trans = Trans + id*(Dim*Dim);
		
	    //SHOW_M(trans, Dim, Dim);
	    trans_left(&F[i][0][0][0], Dim*M, Dim*M, trans);
	    trans_rightT(&F[0][0][i][0], Dim*M, Dim*M, trans);

	    trans_left(&Bt[i][0][0], N, N, trans);
	    trans_rightT(&B[0][i][0], N, Dim*M, trans);
	}
#endif


	/* Global Matrix */
	/* Mat u-p Block (1, *) */
	for (i = 0; i < M; i++) {
	    /* du = 0 at Dirichlet boundary */
	    for (k = 0; k < Dim; k++) {
        if (IF_DB)
        {
		if (phgDofDirichletBC_(ns->du, e, i*Dim+k, NULL, bufu, &rhsu[i][0],
				       DOF_PROJ_NONE)) {
		    //assert(k == 0);  /* Check for slip boundary */
		    phgMatAddEntries(ns->matF, 1, Iu[i] + k, M, Ju[k], bufu);
		    if (use_Fu && k == 0) 
			phgMatAddEntries(_pcd->matFu, 1, Iu1 + i, M, Iu1, bufu);
		}
		else {
		    phgMatAddEntries(ns->matF, 1, Iu[i] + k, M*Dim, Iu[0],
				     &(F[i][k][0][0]));
		    phgMatAddEntries(ns->matBt, 1, &Iu[i][k], N, Ip,
				     &Bt[i][k][0]);
		    if (use_Fu && k == 0) 
			phgMatAddEntries(_pcd->matFu, 1, Iu1 + i, M, Iu1, &(Fu[i][0]));
		    //SHOW_V(F[i][k][0], M*Dim);
		}
        }

#if 1 && USE_NODAL_LOADS
		if (0 && phgDofDirichletBC_(ns->du, e, i*Dim+k, NULL, bufu, &rhsu[i][0],
				       DOF_PROJ_NONE)) {
		    //assert(k == 0);  
		    phgMatAddEntries(ns->matF0, 1, Iu[i] + k, M, Ju[k], bufu);
        }
        else
        {
		    phgMatAddEntries(ns->matF0, 1, Iu[i] + k, M*Dim, Iu[0],
				     &(F[i][k][0][0]));
		    phgMatAddEntries(ns->matBt0, 1, &Iu[i][k], N, Ip,
				     &Bt[i][k][0]);
        }
#endif
	    }
	}

	/* Mat u-p (1, *) */
	for (i = 0; i < N; i++) {
        if (IF_DB)
        {
	    if (phgDofDirichletBC(ns->dp, e, i, NULL, bufp, NULL, DOF_PROJ_NONE)) {
		phgMatAddEntries(ns->matC, 1, Ip + i, N, Ip, bufp);
	    } else {
		phgMatAddEntries(ns->matB, 1, Ip + i, M * Dim, Iu[0],
				      &B[i][0][0]);
		phgMatAddEntries(ns->matC, 1, Ip + i, N, Ip,
				      &C[i][0]);
	    }
        }

#if 1 && USE_NODAL_LOADS
	    if (0 && phgDofDirichletBC(ns->dp, e, i, NULL, bufp, NULL, DOF_PROJ_NONE)) {
		phgMatAddEntries(ns->matC0, 1, Ip + i, N, Ip, bufp);
        }
        else
        {
		phgMatAddEntries(ns->matB0, 1, Ip + i, M * Dim, Iu[0],
				      &B[i][0][0]);
		phgMatAddEntries(ns->matC0, 1, Ip + i, N, Ip,
				      &C[i][0]);
        }
#endif

	} 
	/* SHOW_iV(Ip_2D, N); */
	/* SHOW_iV(Ip_3D, N); */
    }				/* end element */
	
    /*
    phgMatDestroy(&ns->matF0);
    phgMatDestroy(&ns->matB0);
    phgMatDestroy(&ns->matBt0);
    phgMatDestroy(&ns->matC0);
    */

	phgDofFree(&avg_n);

    /* mat check */
#define MAT_CHECK_DUP(mat)    {					\
	MAT_ROW *row = mat->rows;				\
	for (i = 0; i < mat->rmap->nlocal; i++, row++) {	\
	    int k_, ncol = row->ncols;				\
	    INT cols[ncol];					\
	    for (j = 0; j < ncol; j++) {			\
		cols[j] = row->cols[j];				\
		for (k_ = 0; k_ < j; k_++)			\
		    assert(cols[k_] != cols[j]);		\
	    }							\
	}							\
    }
    /* MAT_CHECK_DUP(ns->matF); */
    /* MAT_CHECK_DUP(ns->matB); */
    /* MAT_CHECK_DUP(ns->matBt); */
    /* MAT_CHECK_DUP(ns->matC); */

    if (DUMP_MAT_VEC) {
	phgPrintf("dumping F,B,Bt,C\n");
	//phgMatDumpMATLAB(ns->matF, "F", "F_.m");
	//phgMatDumpMATLAB(ns->matB, "B", "B_.m");
	//phgMatDumpMATLAB(ns->matBt, "Bt", "Bt_.m");
	//phgMatDumpMATLAB(ns->matC, "C", "C_.m");
    }


    /* Exit on checking matrix. */
    if (0) { //&& viscosity_type) {
    	phgFinalize();
    	exit(1);
    }
    return;
}



