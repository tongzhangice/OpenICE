    /* Fp Ap^-1 * q */
    phgMatVec(MAT_OP_N, 1.0, _pcd->matFp, xp, 0., &xp2);
    memcpy(xp->data, xp2->data, sizeof(*xp->data) * Np);

    /* Ap^-1 * q */
    t = phgGetTime(NULL);
    solver_Ap->rhs->data = xp->data;
    solver_Ap->rhs->assembled = TRUE;
    bzero(xp2->data, sizeof(*xp2->data) * Np);
    phgSolverVecSolve(solver_Ap, FALSE, xp2); 
    memcpy(xp->data, xp2->data, sizeof(*xp->data) * Np);
    if (verb > 0)
	phgPrintf("\t    Ap: nits = %3d, residual = %0.4le [%0.4lgMB %0.4lfs]\n",
		  solver_Ap->nits, (double)solver_Ap->residual,
		  phgMemoryUsage(g, NULL) / (1024.0 * 1024.0),
		  phgGetTime(NULL) - t);
    phgVecAXPBY(0., NULL, -1.0, &xp);
