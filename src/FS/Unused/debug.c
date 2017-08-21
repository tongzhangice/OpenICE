/****************************************************************
 * Build RHS which is the residual of the nonlinear system.
 ***************************************************************/
void 
phgNSBuildSolverURHS(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    SOLVER *solver_u = ns->solver_u;
    int i, k, l, q, s;
    FLOAT *dt = ns->dt;
    BOOLEAN tstep_minus = (ns->u[-1] != NULL);
    VEC *vec_rhs = phgMapCreateVec(solver_u->rhs->map, 1);
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    int viscosity_type = ns->viscosity_type;

    SURF_BAS *surf_bas = ns->surf_bas;
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);


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
	phgPrintf("%d ", ns->u[1]->DB_masks[k]);
    phgPrintf("]   ");

    nu_max = -1e10;
    nu_min = +1e10;

    phgVecDisassemble(vec_rhs);
    ForAllElements(g, e) {
	int M = ns->u[1]->type->nbas;	/* num of bases of Velocity */
	int N = ns->p[1]->type->nbas;	/* num of bases of Pressure */
	int order = DofTypeOrder(ns->u[1], e) * 3 - 1; /* Note:
							*   quad order is really high here,
							*   highest order term (u \nabla u, phi)  */
	FLOAT bufu[M], bufp[N], rhsu[M][Dim], rhsp[N];
	INT Iu[M][Dim], Ip[N];
	QUAD *quad;
	FLOAT vol, area, det;
	const FLOAT *w, *p, *normal,
	    **vu, *vu_queue[3],
	    *vf[2], *gu[2], *vp[2], *vw, *vT;
	FLOAT *vf_cache[2];

	vu = vu_queue + 1;

	quad = phgQuadGetQuad3D(order);
	vu[0] = phgQuadGetDofValues(e, ns->u[0], quad);	           /* u^{n} */
	vp[0] = phgQuadGetDofValues(e, ns->p[0], quad);	           /* p^{n} */
	gu[0] = phgQuadGetDofValues(e, ns->gradu[0], quad);        /* grad u^{n} */
	vw = phgQuadGetDofValues(e, ns->wind, quad);               /* wind */
	vT = phgQuadGetDofValues(e, ns->T[1], quad);	           /* T^{n} */


	if (tstep_minus) { 
	    vu[-1] = phgQuadGetDofValues(e, ns->u[-1], quad);      /* u^{n-1} */
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
	Unused(vf); Unused(vf_cache); 

	if (!_nsp->no_extern_source) {
	    /* cache f values */
	    for (l = 0; l < 2; l++) {
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
    
	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);

	    /* rhs u */
	    for (i = 0; i < M; i++) {
		/* interior node or Neumann */
		const FLOAT *gi_u = phgQuadGetBasisValues(e, ns->u[1], i, quad) + q;       /* phi_i */
		const FLOAT *ggi_u = phgQuadGetBasisCurvedGradient(e, ns->u[1], i, quad, q);    /* grad phi_i */

		for (k = 0; k < Dim; k++) {
#if ICE_BENCH_TEST
		    nu = get_effective_viscosity(gu[1], 0, 0, viscosity_type);
		    FLOAT eu[DDim];

		    MAT3_SYM(gu[1], eu);
		    rhsu[i][k] += vol*(*w) * EQU_SCALING * (- nu * INNER_PRODUCT(eu+k*Dim, ggi_u) 
							    + (*vp[1]) * *(ggi_u+k) * LEN_SCALING * PRES_SCALING
							    );     /* left */

		    if (k == Z_DIR) { 
			const FLOAT rho = RHO_ICE;
			const FLOAT grav = GRAVITY;
			const FLOAT a = SEC_PER_YEAR;
			const FLOAT f = rho*grav * EQU_SCALING * LEN_SCALING2; 

			Unused(a);
			rhsu[i][k] += vol*(*w) * (-f * (*gi_u) 
						  ); /* right */
		    }


#elif ESIMINT_TEST ||				\
    HEINO_TEST ||				\
    TEST_CASE == ICE_GREEN_LAND
		    nu = get_effective_viscosity(gu[1], *vT, 0, viscosity_type);
		    FLOAT eu[DDim];

		    MAT3_SYM(gu[1], eu);
		    rhsu[i][k] += vol*(*w) * EQU_SCALING * (- nu * INNER_PRODUCT(eu+k*Dim, ggi_u) 
							    + (*vp[1]) * *(ggi_u+k) * LEN_SCALING * PRES_SCALING
							    );     /* left */

		    if (k == Z_DIR) { 
			const FLOAT rho = RHO_ICE;
			const FLOAT grav = GRAVITY;
			const FLOAT a = SEC_PER_YEAR;
			const FLOAT f = rho*grav * EQU_SCALING * LEN_SCALING2; 

			Unused(a);
			rhsu[i][k] += vol*(*w) * (-f * (*gi_u) 
						  ); /* right */
		    }

#elif STEADY_STATE
		    rhsu[i][k] += vol*(*w) * (- nu * INNER_PRODUCT(gu[1]+k*Dim, ggi_u)
					      + (*vp[1]) * *(ggi_u+k)
					      );     /* left */
		    if (!_nsp->no_extern_source)
			rhsu[i][k] += vol*(*w) * (*(vf[1]+k) * (*gi_u)
						  ); /* right */
#elif TIME_DEP_NON
		    rhsu[i][k] -= vol*(*w) * ((vu[1][k] - vu[0][k]) * (*gi_u) / dt[0]
					      + Theta * (nu * INNER_PRODUCT(gu[1]+k*Dim, ggi_u)
							 )
					      - (*vp[1]) * *(ggi_u+k)
					      + Thet1 * (nu * INNER_PRODUCT(gu[0]+k*Dim, ggi_u)
							 )
					      );     /* left */
		    if (!_nsp->no_extern_source)
			rhsu[i][k] += vol*(*w) * (Theta * *(vf[1]+k) * (*gi_u)
						  + Thet1 * *(vf[0]+k) * (*gi_u)
						  ); /* right */
#else
		    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE */
		}
	    }

	    /* rhs p */
	    for (i = 0; i < N; i++) {
		const FLOAT *gi_p = phgQuadGetBasisValues(e, ns->p[1], i, quad) + q;       /* psi_i */
		FLOAT divu1 = gu[1][0] + gu[1][4] + gu[1][8];
		//FLOAT divu0 = gu[0][0] + gu[0][4] + gu[0][8];
		rhsp[i] += vol*(*w) * (divu1 * (*gi_p)
				       );
	    }
	    
	    if (tstep_minus) 
		vu[-1] += Dim;

#if STEADY_STATE || TIME_DEP_NON
	    vu[1] += Dim;
	    gu[1] += Dim*Dim;
	    vp[1]++;
#else
	    TIME_DEP_LINEAR; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */
	    vu[0] += Dim;
	    gu[0] += Dim * Dim;
	    vp[0]++; 
	    vw += Dim;
	    if (!_nsp->no_extern_source) {
		vf[0] += Dim; vf[1] += Dim;
	    }
	    vT++;
	    w++; p += Dim + 1;
	}

	if (!_nsp->no_extern_source) {
	    phgFree(vf_cache[0]);
	    phgFree(vf_cache[1]);
	}

	normal = NULL; Unused(normal);
	area = 0; Unused(area);

	if (!_nsp->enclosed_flow) {
	    /* slip boundary */
	    for (s = 0; s < NFace; s++) {
		if (e->bound_type[s] & INFLOW) {
		    int v0, v1, v2;
		    int nbas_face = NbasFace(ns->u[1]);
		    SHORT bases[nbas_face];
		    FLOAT lambda[Dim + 1], x,y,z, beta;
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
			FLOAT vu[Dim];
			lambda[v0] = *(p++);
			lambda[v1] = *(p++);
			lambda[v2] = *(p++);
			
			phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
			func_beta(x, y, z, &beta);
			phgDofEval(ns->u[1], e, lambda, vu);

			for (i = 0; i < nbas_face; i++) {
			    int ii = bases[i];
			    FLOAT gi_u = 
				*ns->u[1]->type->BasFuncs(ns->u[1], e, ii, ii + 1, lambda);

			    for (k = 0; k < Dim; k++) {
#if STEADY_STATE
				rhus[ii][k] += 0.;
#elif TIME_DEP_NON
#  if USE_SLIDING_BC
				abort();
				rhsu[ii][k] += SIGN_FRICTION * area*(*w) * beta * vu[k] * (gi_u)
				    * EQU_SCALING * LEN_SCALING;
#  else
				Unused(gi_u);
#  endif
#else
				TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE */
			    }
			}     /* end of bas_i */
			w++;
		    }		/* end of quad point */
		}		/* end of face outflow */
	    }			/* end of all outflow face in element */
	}                       /* end out flow boundary */


#if USE_SLIDING_BC
	/* Rotate bases */
	for (i = 0; i < M; i++) {
	    INT id = phgDofMapE2D(surf_dof, e, i * (Dim*Dim)) / (Dim*Dim);
	    if (!rotated[id])
		continue;	
	    const FLOAT *trans = Trans + id*(Dim*Dim);

	    trans_left(&rhsu[i][0], 1, 1, trans);
	}
#else
	Unused(Trans);
	Unused(rotated);
#endif


	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    for (k = 0; k < Dim; k++)
		Iu[i][k] = phgMapE2L(solver_u->rhs->map, 0, e, i * Dim + k);
	for (i = 0; i < N; i++)
	    Ip[i] = phgMapE2L(solver_u->rhs->map, 1, e, i);

	/* set velocity dirichlet bdry */
	FLOAT tmp[Dim];
	for (i = 0; i < M; i++)
	    for (k = 0; k < Dim; k++)
		if (phgDofDirichletBC_(ns->u[1], e, i*Dim+k, NULL, bufu, tmp,
				       DOF_PROJ_NONE)) {
		    rhsu[i][k] = 0.;
		}

#if STEADY_STATE || TIME_DEP_NON
	/* set pressure dirichlet bdry for pinned point */
	for (i = 0; i < N; i++)
	    if (phgDofDirichletBC(ns->p[1], e, i, NULL, bufp, &rhsp[i],
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


	/* Global res */
	phgVecAddEntries(vec_rhs, 0, M * Dim, Iu[0], &rhsu[0][0]);
	phgVecAddEntries(vec_rhs, 0, N, Ip, rhsp);
    }				/* end element */
    
    phgVecAssemble(vec_rhs);
    phgVecAssemble(solver_u->rhs);
    phgVecAXPBY(1., vec_rhs, 0, &solver_u->rhs);
    solver_u->rhs_updated = FALSE;

    if (DUMP_MAT_VEC) {
	phgPrintf("Dumping rhs\n");
	phgVecDumpMATLAB(ns->solver_u->rhs, "b", "b_.m");
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


void
phgNSBuildSolverUMat(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, k, l, q, s;
    FLOAT *dt = ns->dt;
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    BOOLEAN use_Fu = _nsp->use_Fu;
    int viscosity_type = ns->viscosity_type;

    SURF_BAS *surf_bas = ns->surf_bas;
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);


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

    ForAllElements(g, e) {
	int M = ns->u[1]->type->nbas;	/* num of bases of Velocity */
	int N = ns->p[1]->type->nbas;	/* num of bases of Pressure */
	int order = DofTypeOrder(ns->u[1], e) * 3 - 1; /* Note:
							*   quad order is really high here,
							*   highest order term (u \nabla u, phi)  */
	FLOAT F[M][Dim][M][Dim], Fu[M][M],
	    B[N][M][Dim], Bt[M][Dim][N], C[N][N],
	    bufu[M], bufp[N],
	    rhsu[M][Dim];
	INT Iu[M][Dim], Iu1[M], Ju[Dim][M], Ip[N];
	QUAD *quad;
	FLOAT vol, det;
	const FLOAT *w, *p, *vw, *gu, *vT;

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
	vw = phgQuadGetDofValues(e, ns->wind, quad);      /* value wind */
	gu = phgQuadGetDofValues(e, ns->gradu[1], quad);  /* grad u */
	vT = phgQuadGetDofValues(e, ns->T[1], quad);	  /* T^{n} */

	Bzero(F); Bzero(Fu); 
	Bzero(B); Bzero(Bt); Bzero(C);

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);
	    /* Test func vel type */
	    for (i = 0; i < M; i++) {
		const FLOAT *gi_u = phgQuadGetBasisValues(e, ns->u[1], i, quad) + q;       /* phi_i */
		const FLOAT *ggi_u = phgQuadGetBasisCurvedGradient(e, ns->u[1], i, quad, q);    /* grad phi_i */

		/* Mat F */
		for (j = 0; j < M; j++) {
		    const FLOAT *gj_u = phgQuadGetBasisValues(e, ns->u[1], j, quad) + q;       /* phi_j */
		    const FLOAT *ggj_u = phgQuadGetBasisCurvedGradient(e, ns->u[1], j, quad, q);    /* grad phi_i */
		    FLOAT mass = (*gj_u) * (*gi_u);
		    FLOAT diffu = INNER_PRODUCT(ggj_u, ggi_u);

#if ICE_BENCH_TEST
		    Unused(mass);
		    nu = get_effective_viscosity(gu, 0, 0, viscosity_type);

		    const FLOAT *tp = get_gbas_product(ggi_u, ggj_u);

		    for (k = 0; k < Dim; k++) 
			for (l = 0; l < Dim; l++) 
			    F[j][l][i][k] += vol*(*w) * EQU_SCALING * nu * tp[k+l*Dim];

		    if (use_Fu) 
			Fu[i][j] += vol*(*w) * EQU_SCALING * (nu * diffu);

#elif ESIMINT_TEST ||				\
    HEINO_TEST ||				\
    TEST_CASE == ICE_GREEN_LAND
		    Unused(mass);
		    nu = get_effective_viscosity(gu, *vT, 0, viscosity_type);

		    const FLOAT *tp = get_gbas_product(ggi_u, ggj_u);

		    for (k = 0; k < Dim; k++) 
			for (l = 0; l < Dim; l++) 
			    F[j][l][i][k] += vol*(*w) * EQU_SCALING * nu * tp[k+l*Dim];

		    if (use_Fu) 
			Fu[i][j] += vol*(*w) * EQU_SCALING * (nu * diffu);

#elif STEADY_STATE
		    Unused(mass);
		    Unused(mass);
		    for (k = 0; k < Dim; k++) {
			F[i][k][j][k] += vol*(*w) * (nu * diffu);
		    }
		    if (use_Fu) {
			Fu[i][j] += vol*(*w) * (nu * diffu);
		    }
#elif TIME_DEP_NON
		    mass /= dt[0];
		    diffu *= Theta * nu;
		    for (k = 0; k < Dim; k++) 
			F[i][k][j][k] += vol*(*w) * (mass + diffu);
		    if (use_Fu) 
			Fu[i][j] += vol*(*w) * (mass + diffu);
#else
		    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE */
		}

		/* Mat B & Bt */
		for (j = 0; j < N; j++) {
		    const FLOAT *gj_p = phgQuadGetBasisValues(e, ns->p[1], j, quad) + q;       /* psi_j */
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
		const FLOAT *gi_p = phgQuadGetBasisValues(e, ns->p[1], i, quad) + q;       /* psi_i */
		C[i][i] += vol*(*w) * (_nsp->eps_diagP * (*gi_p) * (*gi_p)
				       );
	    }

	    /* Next quad point */
	    vw += Dim;
	    gu += Dim*Dim;
	    vT++;
	    w++; p += Dim+1;
	}

	if (!_nsp->enclosed_flow) {
	    /* slip boundary */
	    for (s = 0; s < NFace; s++) {
		if (e->bound_type[s] & SLIP_BDRY) {
		    int v0, v1, v2;
		    int nbas_face = NbasFace(ns->u[1]);
		    SHORT bases[nbas_face];
		    FLOAT lambda[Dim + 1], area, x, y, z, beta;
		    order = DofTypeOrder(ns->u[1], e) * 3 - 1;

		    phgDofGetBasesOnFace(ns->u[1], e, s, bases);
		    v0 = GetFaceVertex(s, 0);
		    v1 = GetFaceVertex(s, 1);
		    v2 = GetFaceVertex(s, 2);
		    lambda[s] = 0.;

		    area = phgGeomGetFaceArea(g, e, s);
		    quad = phgQuadGetQuad2D(order);

		    p = quad->points;
		    w = quad->weights;
		    for (q = 0; q < quad->npoints; q++) {
			lambda[v0] = *(p++);
			lambda[v1] = *(p++);
			lambda[v2] = *(p++);
			
			phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
			func_beta(x, y, z, &beta);

			for (i = 0; i < nbas_face; i++) {
			    int ii = bases[i];
			    FLOAT gi_u = 
				*ns->u[1]->type->BasFuncs(ns->u[1], e, ii, ii + 1, lambda);
			    for (j = 0; j < nbas_face; j++) { 
				int jj = bases[j];
				FLOAT gj_u = 
				    *ns->u[1]->type->BasFuncs(ns->u[1], e, jj, jj + 1, lambda);
				FLOAT mass_face = area*(*w) * beta * (gj_u)*(gi_u)
				    * EQU_SCALING * LEN_SCALING;

				for (k = 0; k < Dim; k++) {
#if STEADY_STATE
				    F[ii][k][jj][k] -= SIGN_FRICTION * mass_face;
#elif TIME_DEP_NON
#  if USE_SLIDING_BC
				    F[ii][k][jj][k] -= SIGN_FRICTION * mass_face;
#  else
				    Unused(mass_face);
#  endif
#else
				    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE */
				}
			    } /* end of bas_j */
			}     /* end of bas_i */
			w++;
		    }		/* end of quad point */
		}		/* end of face outflow */
	    }			/* end of all outflow face in element */
	}                       /* end out flow boundary */

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
#else
	Unused(Trans);
	Unused(rotated);
#endif


	/* Global Matrix */
	/* Mat u-p Block (1, *) */
	for (i = 0; i < M; i++) {
	    /* du = 0 at Dirichlet boundary */
	    for (k = 0; k < Dim; k++) {
		if (phgDofDirichletBC_(ns->u[1], e, i*Dim+k, NULL, bufu, &rhsu[i][0],
				  DOF_PROJ_NONE)) {
		    //assert(k == 0);
		    phgMatAddEntries(ns->matF, 1, Iu[i] + k, M, Ju[k], bufu);
		    if (use_Fu && k == 0) 
			phgMatAddEntries(_pcd->matFu, 1, Iu1 + i, M, Iu1, bufu);
		} else {
		    phgMatAddEntries(ns->matF, 1, Iu[i] + k, M*Dim, Iu[0],
				     &(F[i][k][0][0]));
		    phgMatAddEntries(ns->matBt, 1, &Iu[i][k], N, Ip,
				     &Bt[i][k][0]);
		    if (use_Fu && k == 0) 
			phgMatAddEntries(_pcd->matFu, 1, Iu1 + i, M, Iu1, &(Fu[i][0]));
		    //SHOW_V(F[i][k][0], M*Dim);
		}
	    }
	}

	/* Mat u-p (1, *) */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(ns->p[1], e, i, NULL, bufp, NULL, DOF_PROJ_NONE)) {
		phgMatAddEntries(ns->matC, 1, Ip + i, N, Ip, bufp);
	    } else {
		phgMatAddEntries(ns->matB, 1, Ip + i, M * Dim, Iu[0],
				      &B[i][0][0]);
		phgMatAddEntries(ns->matC, 1, Ip + i, N, Ip,
				      &C[i][0]);
	    }
	} 
	/* SHOW_iV(Ip_2D, N); */
	/* SHOW_iV(Ip_3D, N); */
    }				/* end element */


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
    MAT_CHECK_DUP(ns->matF);
    MAT_CHECK_DUP(ns->matB);
    MAT_CHECK_DUP(ns->matBt);
    MAT_CHECK_DUP(ns->matC);

    if (DUMP_MAT_VEC) {
	phgPrintf("dumping F,B,Bt,C\n");
	phgMatDumpMATLAB(ns->matF, "F", "F_.m");
	phgMatDumpMATLAB(ns->matB, "B", "B_.m");
	phgMatDumpMATLAB(ns->matBt, "Bt", "Bt_.m");
	phgMatDumpMATLAB(ns->matC, "C", "C_.m");
    }

    if (0) { //&& viscosity_type) {
    	phgFinalize();
    	exit(1);
    }
    return;
}

/*
 *  Subroutines of Building Mat & RHS of Themarl solver
 *
 *
 *  */
#include "ins.h"
#include "mat_op3.h"
#define _nsp (ns->ns_params)
#define _pcd (ns->pcd)
#define DDim (Dim*Dim)


void phgNSInitSolverT(NSSolver *ns)
{
    MAP *T_map;

    /* dof copy */
    ns->T_shape = phgDofCopy(ns->T[1], NULL, NULL, "T shape");

    /* dof map */
    T_map = ns->T_map = phgMapCreate(ns->T_shape, NULL);

    /* matrices */
    ns->matT = phgMapCreateMat(T_map, T_map);
    ns->matT->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;

    /* solver_T */
    phgOptionsPush();
    phgOptionsSetOptions("-solver petsc "
			 "-petsc_ksp_type cg "
			 "-petsc_pc_type bjacobi "
			 "-solver_maxit 1000 "
			 "-solver_rtol 1e-10");

    phgOptionsSetOptions(_nsp->T_opts);
    ns->solver_T = phgMat2Solver(SOLVER_DEFAULT, ns->matT);
    phgOptionsPop();
    phgVecDisassemble(ns->solver_T->rhs);
    phgDofSetDataByValue(ns->T[1], 0.);
    return;
}


void phgNSDestroySolverT(NSSolver *ns)
{
    phgInfo(2, "   Destroy solver T\n");
    phgMatDestroy(&ns->matT);
    if (ns->T_cntr != NULL ) {
	phgFree(ns->T_mask);
	phgVecDestroy(&ns->T_cntr);
    }
    phgSolverDestroy(&ns->solver_T);
    ns->solver_T = NULL;
    phgMapDestroy(&ns->T_map);
    phgDofFree(&ns->T_shape);
}



void 
phgNSBuildSolverTMat(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    SOLVER *solver_T = ns->solver_T;
    DOF **T = ns->T;
    int i, j, q;
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    FLOAT *dt = ns->dt;
    int viscosity_type = ns->viscosity_type;

    Unused(dt);
    ForAllElements(g, e) {
	int M = T[1]->type->nbas;	/* num of bases of Velocity */
	int order = DofTypeOrder(T[1], e) * 2;
	FLOAT A[M][M], bufT[M], rhsT[M];
	INT I[M];
	QUAD *quad;
	FLOAT vol, det;
	const FLOAT *w, *p, *gu, *vu, *vT, *gT;

	Bzero(A); Bzero(bufT); Bzero(rhsT); 
	quad = phgQuadGetQuad3D(order);
	vu = phgQuadGetDofValues(e, ns->u[1], quad); 
	gu = phgQuadGetDofValues(e, ns->gradu[1], quad); 
	vT = phgQuadGetDofValues(e, ns->T[1], quad); 
	gT = phgQuadGetDofValues(e, ns->gradT[1], quad); 

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);

	    nu = get_effective_viscosity(gu, *vT, 0, viscosity_type);

	    for (i = 0; i < M; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, T[1], i, quad) + q;    
		const FLOAT *ggi = phgQuadGetBasisCurvedGradient(e, T[1], i, quad, q); 
		for (j = 0; j < M; j++) {
		    const FLOAT *gj = phgQuadGetBasisValues(e, T[1], j, quad) + q;       
		    const FLOAT *ggj = phgQuadGetBasisCurvedGradient(e, T[1], j, quad, q); 

		    FLOAT qmass = (*gj) * (*gi);
		    FLOAT conv = INNER_PRODUCT(vu, ggj) * (*gi);
		    FLOAT diff_T = INNER_PRODUCT(ggi, ggj);
		    
		    const FLOAT rho = RHO_ICE;
		    const FLOAT c = HEAT_CAPACITY;
		    const FLOAT K = THEM_CONDUCT;

		    A[i][j] += vol*(*w) * (//rho * c * qmass / dt +
					   + rho * c * conv
					   + K * diff_T
					   );
		}
	    }

	    vu += Dim;
	    gu += DDim;
	    vT++;
	    gT += Dim;
	    w++; p += Dim + 1;
	}

	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    I[i] = phgMapE2L(ns->matT->cmap, 0, e, i);

	/* Global res */
	for (i = 0; i < M; i++) {
	    if (phgDofDirichletBC_(ns->T[1], e, i, NULL, bufT, &rhsT[i],
				   DOF_PROJ_NONE)) {
		phgMatAddEntries(ns->matT, 1, I + i, M, I, bufT);
	    } else {
		phgMatAddEntries(ns->matT, 1, I + i, M, I, &(A[i][0])); 
	    }
	}
    }				/* end element */
    solver_T->rhs_updated = FALSE;

    if (DUMP_MAT_VEC) {
	phgPrintf("Dumping MatT\n");
	phgMatDumpMATLAB(ns->solver_T->mat, "A_gu", "mat_gu_.m");
    }
  
    return;
}


void 
phgNSBuildSolverTRHS(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    SOLVER *solver_T = ns->solver_T;
    DOF **T = ns->T;
    int i, j, q;
    FLOAT *dt = ns->dt;
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    VEC *vec_rhs = phgMapCreateVec(solver_T->rhs->map, 1);
    int viscosity_type = ns->viscosity_type;
    
    Unused(dt);
    Unused(j);
    phgVecDisassemble(vec_rhs);
    ForAllElements(g, e) {
	int M = T[1]->type->nbas;	/* num of bases of Velocity */
	int order = DofTypeOrder(T[1], e) * 2;
	FLOAT rhs[M], buf[M];
	INT I[M];
	QUAD *quad;
	FLOAT vol, det;
	const FLOAT *w, *p, *vu, *gu, *vT, *gT;

	Bzero(rhs);
	quad = phgQuadGetQuad3D(order);
	vu = phgQuadGetDofValues(e, ns->u[1], quad); 
	gu = phgQuadGetDofValues(e, ns->gradu[1], quad); 
	vT = phgQuadGetDofValues(e, ns->T[1], quad); 
	gT = phgQuadGetDofValues(e, ns->gradT[1], quad); 

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);

	    nu = get_effective_viscosity(gu, *vT, 0, viscosity_type);

	    for (i = 0; i < M; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, T[1], i, quad) + q;    
		const FLOAT *ggi = phgQuadGetBasisCurvedGradient(e, T[1], i, quad, q); 

		FLOAT conv = INNER_PRODUCT(vu, gT) * (*gi);
		FLOAT diff_T = INNER_PRODUCT(gT, ggi);
		FLOAT diff_u = MAT3_NORM2_2(gu) * (*gi);
		    
		const FLOAT rho = RHO_ICE;
		const FLOAT c = HEAT_CAPACITY;
		const FLOAT K = THEM_CONDUCT;
		const FLOAT q = HEAT_SOURCE;

		rhs[i] -= vol*(*w) * (//rho * c * (T[1] - T[0]) * (*gi)/ dt 
				      + rho * c * conv / LEN_SCALING
				      + K  * diff_T / LEN_SCALING2
				      - nu * diff_u / LEN_SCALING2
				      + q);
	    }

	    vu += Dim;
	    gu += DDim;
	    vT++;
	    gT += Dim;
	    w++; p += Dim + 1;
	}

	/* Geo thermal flux */
	for (s = 0; s < NFace; s++) {
	    if (e->bound_type[s] & BC_BOTTOM) {
		int v0, v1, v2;
		int nbas_face = NbasFace(ns->T[1]);
		SHORT bases[nbas_face];
		FLOAT lambda[Dim + 1];
		order = DofTypeOrder(ns->T[1], e) * 3 - 1;

		phgDofGetBasesOnFace(ns->T[1], e, s, bases);
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
		    lambda[v0] = *(p++);
		    lambda[v1] = *(p++);
		    lambda[v2] = *(p++);
			
		    for (i = 0; i < nbas_face; i++) {
			int ii = bases[i];
			FLOAT gi_T = 
			    *ns->T[1]->type->BasFuncs(ns->T[1], e, ii, ii + 1, lambda);

			rhsu[ii] -= area*(*w) * (GEOTHE_FLUX / THEM_CONDUCT) * (gi_u)
			    * EQU_SCALING * LEN_SCALING;
		    } /* end of bas_i */
		    w++;
		} /* end of quad point */
	    }	  /* end of face outflow */
	}	  /* end of faces */

	FLOAT tmp;
	for (i = 0; i < M; i++) 
	    if (phgDofDirichletBC_(ns->T[1], e, i, NULL, buf, &tmp,
				   DOF_PROJ_NONE)) 
		rhs[i] = 0.;

	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    I[i] = phgMapE2L(ns->matT->cmap, 0, e, i);

	phgVecAddEntries(vec_rhs, 0, M, I, &rhs[0]);
    }				/* end element */
    
    phgVecAssemble(vec_rhs);
    phgVecAssemble(solver_T->rhs);
    phgVecAXPBY(1., vec_rhs, 0, &solver_T->rhs);
    solver_T->rhs_updated = FALSE;

    if (DUMP_MAT_VEC) {
	phgPrintf("Dumping rhsT\n");
	phgVecDumpMATLAB(ns->solver_T->rhs, "b_gu", "rhs_gu_.m");
    }

    phgVecDestroy(&vec_rhs);
    return;
}


/*
 * See EISMINT Benchmark.
 * At the ice-bedrock interface, mixed boundary conditions are emplyerd:
 *   if T < T', dT/dz = -G/k;
 *   else,      T = T',
 *   where T' = T_0 - beta H.
 * Here we build this constrain for vector, rather than DOF,
 *   since we may use this in the solvering interation.
 *  
 * */
void 
phgNSBuildSolverTConstrain(NSSolver *ns)
{
    GRID *g;
    MAP *T_map = ns->T_map;
    BOOLEAN *T_mask;
    FLOAT *T_cntr;
    LAYERED_MESH *gL = ns->gL;
    int ii, i, j;

    phgCalloc(T_mask, T_map->nlocal);
    phgCalloc(T_cntr, T_map->nlocal);
        
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	assert(nv > 0);

	FLOAT dh = g->verts[iL[nv-1]][Z_DIR]
	    - g->verts[iL[0]][Z_DIR];
	T_mask[i] = TRUE;
	T_cntr[i] = T_0 - beta * dh;
    }

    ns->T_mask = T_mask;
    ns->T_cntr = T_cntr;
    return;
}

static void 
mg_GaussSidel2_forward(MAT *A, VEC *x, VEC *b, int nsmooth, void *ctx);

/* Solve T:
 * Contraine Vec value after each GS smooth.
 * */
void 
phgNSSolverTSolve(NSSolver *ns)
{
    SOLVER *solver = ns->solver_T;
    int max_it = ns->solver->max_it;
    FLOAT res0, res, rtol = ns->solver->rtol;
    int ndof = 1;
    DOF *dofs[] = {ns->T[1], NULL};
    int i, j;

    phgPrintf("Solver T: constrained solve\n");

    /* copy from phgSolverSolve */
    VEC *x = phgMapCreateVec(solver->rhs->map, 1);
    /* copy initial solution to x */
    phgMapDofToLocalData(x->map, ndof, dofs, x->data);

    /* copy from phgSolverVecSolve */
    assert(solver->rhs != NULL && solver->mat != NULL);

    if (solver->mat->cmap->nlocal != x->map->nlocal) {
	phgError(1, "(%s:%d) inconsistent RHS!\n", __FILE__, __LINE__);
    }

    if (solver->mat->handle_bdry_eqns == TRUE &&
	   solver->rhs->mat != solver->mat)
	phgError(1, "(%s:%d) solver->rhs->mat != solver->mat.\n",
			__FILE__, __LINE__);

    phgSolverAssemble(solver);

    if (!solver->rhs->assembled)
	phgVecAssemble(solver->rhs);

    if (!solver->rhs_updated)
	phgSolverUpdateRHS(solver);

    MAT *A = solver->mat;
    VEC *b = solver->rhs;
    VEC *r = phgVecCopy(b, NULL);
    phgMatVec(MAT_OP_N, -1., A, x, 1., &r);
    res0 = phgVecNorm2(r, 0, NULL);

    if (res0 < 1e-10)
	res0 = 1.;

    for (k = 0; k < max_it; k++) {
	mg_GaussSidel2_forward(A, x, b, 1, &ns);
	res = solver->residual;
	res /= res0;
	phgPrintf("   nit: %3d\t, res: %e\n", 
		  k, res);
	if (res < rtol)
	    break;
    }
    solver->residual = res;
    solver->max_it = k;

    /* copy solution from x */
    phgMapLocalDataToDof(x->map, ndof, dofs, x->data);

    phgVecDestroy(&x);
    phgVecDestroy(&r);
    return;
}




/*
 * Gauss-Sidel smoother2:
 *
 *   1. interior dof (!REMOTE) is smoothed, and then offp data is updated; 
 *   2. proc boundary dof (REMOTE) is smoothed, and then offp data is updated.
 *
 * Note: need types_vec.
 *  */
static void 
mg_GaussSidel2_forward(MAT *A, VEC *x, VEC *b, int nsmooth, void *ctx)
{
    NSSolver *ns = (NSSolver *) ctx;
    GRID *g = ng->g;
    SOLVER *solver = ns->solver_T;
    INT i, j, k, n, *pc, *pc_offp, *pc0, nlocal;
    FLOAT *pd, *pd0, *pd_offp, *vx, *vb;
    size_t *ps;
    FLOAT sum, omega = 1.;
    FLOAT res;
    BOOLEAN *T_mask = ns->T_mask;
    FLOAT *T_cntr = ns->T_cntr;
#if USE_MPI
    FLOAT *offp_data = NULL;
#endif	/* USE_MPI */


    MagicCheck(VEC, x);
    MagicCheck(VEC, b);
    assert(x == NULL || x->nvec == 1);
    assert(b == NULL || b->nvec == 1);
    if (x != NULL && !x->assembled)
	phgVecAssemble(x);
    if (b != NULL && !b->assembled)
	phgVecAssemble(b);
    
    assert(A->type != PHG_DESTROYED);
    if (!A->assembled)
	phgMatAssemble(A);
    assert(A->type != PHG_MATRIX_FREE);
    phgMatPack(A);

#if USE_MPI
    if (A->cmap->nprocs > 1) {
	offp_data = phgAlloc(A->cinfo->rsize * sizeof(*offp_data));
    }
#endif	/* USE_MPI */

    if (A->cmap->nlocal != A->rmap->nlocal ||
	A->cmap->nlocal != x->map->nlocal  ||
	A->cmap->nlocal != b->map->nlocal)
	phgError(1, "%s:%d: inconsistent matrix-vector.", __FILE__,
		 __LINE__);
    nlocal = A->cmap->nlocal;

#if USE_MPI
    if (A->cmap->nprocs > 1) {
	phgMapScatterBegin(A->cinfo, x->nvec, x->data, offp_data);
	phgMapScatterEnd(A->cinfo, x->nvec, x->data, offp_data);
    }
#endif	/* USE_MPI */

    /* iteration */
    for (k = 0; k < nsmooth; k++) {
	res = 0.;

	/* multiply with local data */
	vx = x->data;
	vb = b->data;

	pc0 = A->packed_cols;
	pd0 = A->packed_data;
	ps = A->packed_ind;

	/* First, interior dof (!REMOTE) is smoothed, and then offp data is updated */
	for (i = 0; i < nlocal; i++) {
	    INT jcol;
	    FLOAT aa = 0., dx;

	    pc = pc0 + PACK_COL(ps, i);
#if 0
	    if (ml->types_vec[i] & REMOTE) {
		continue;
	    }
#else
	    if (A->cmap->nprocs > 1 &&
		pc0[PACK_COL_OFFP(ps, i, nlocal)] != 0)
		continue;
#endif

	    sum = vb[i];
	    /* local data */
	    if ((n = *(pc++)) != 0) {
		pd = pd0 + PACK_DAT(ps, i);
		for (j = 0; j < n; j++) {
		    jcol = *(pc++);
		    if (jcol != i) {
			sum -= *(pd++) * vx[jcol];
		    } else {
			aa = *(pd++);
			assert(fabs(aa) > 1e-14);
		    }
		}
	    }
	    
	    /* remote data */
	    if (A->cmap->nprocs > 1) {
		pc_offp = pc0 + PACK_COL_OFFP(ps, i, nlocal);
		if ((n = *(pc_offp++)) != 0) {
		    pd_offp = pd0 + PACK_DAT_OFFP(ps, i, nlocal);
		    for (j = 0; j < n; j++) {
			jcol = *(pc_offp++);
			sum -= *(pd_offp++) * offp_data[jcol];
		    }
		}
	    }
	    
	    dx = sum / aa - vx[i];
	    vx[i] += omega * dx;
	    if (T_mask[i] && vx[i] > T_cntr[i]) {/* constrained */
		vx[i] = T_cntr[i];
	    } else {
		res += (dx*aa) * (dx*aa);        /* use last residual */
	    }
	}

#if USE_MPI
	if (A->cmap->nprocs > 1) {
	    phgMapScatterBegin(A->cinfo, x->nvec, x->data, offp_data);
	    phgMapScatterEnd(A->cinfo, x->nvec, x->data, offp_data);
	}
#endif	/* USE_MPI */

	/* Second, proc boundary dof (REMOTE) is smoothed, and then offp data is updated */
	for (i = 0; i < nlocal; i++) {
	    INT jcol;
	    FLOAT aa = 0., dx;

	    pc = pc0 + PACK_COL(ps, i);
#if 0
	    if (!(ml->types_vec[i] & REMOTE)) {
		continue;
	    }
#else
	    if (A->cmap->nprocs > 1 &&
		pc0[PACK_COL_OFFP(ps, i, nlocal)] == 0)
		continue;
#endif

	    sum = vb[i];
	    /* local data */
	    if ((n = *(pc++)) != 0) {
		pd = pd0 + PACK_DAT(ps, i);
		for (j = 0; j < n; j++) {
		    jcol = *(pc++);
		    if (jcol != i) {
			sum -= *(pd++) * vx[jcol];
		    } else {
			aa = *(pd++);
			assert(fabs(aa) > 1e-14);
		    }
		}
	    }
	    
	    /* remote data */
	    if (A->cmap->nprocs > 1) {
		pc_offp = pc0 + PACK_COL_OFFP(ps, i, nlocal);
		if ((n = *(pc_offp++)) != 0) {
		    pd_offp = pd0 + PACK_DAT_OFFP(ps, i, nlocal);
		    for (j = 0; j < n; j++) {
			jcol = *(pc_offp++);
			sum -= *(pd_offp++) * offp_data[jcol];
		    }
		}
	    }
	    
	    dx = sum / aa - vx[i];
	    vx[i] += omega * dx;
	    if (T_mask[i] && vx[i] > T_cntr[i]) { /* constrained */
		vx[i] = T_cntr[i];
	    } else {
		res += (dx*aa) * (dx*aa);         /* use last residual */
	    }
	}

#if USE_MPI
	if (A->cmap->nprocs > 1) {
	    phgMapScatterBegin(A->cinfo, x->nvec, x->data, offp_data);
	    phgMapScatterEnd(A->cinfo, x->nvec, x->data, offp_data);
	}
#endif	/* USE_MPI */
    }


#if USE_MPI
    phgFree(offp_data);
#endif	/* USE_MPI */

    solver->residual = sqrt(res);
    return;
}
