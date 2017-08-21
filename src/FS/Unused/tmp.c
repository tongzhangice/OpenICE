DOF *
phgDofCopy_(DOF *src, DOF **dest_ptr, DOF_TYPE *newtype,
	    const char *name, const char *srcfile, int srcline)
/* returns a new DOF whose type is 'newtype', and whose values are evaluated
 * using 'src' (a copy of src). */
{
    GRID *g = src->g;
    SIMPLEX *e;
    FLOAT w = 1.0, *basvalues = NULL, *a, *d, *buffer = NULL;
    int i, k, dim, dim_new, count = 0, nvalues, nd;
    INT n;
    DOF *wgts = NULL;
    DOF *dest = (dest_ptr == NULL ? NULL : *dest_ptr);
    BYTE *flags0, *flags;
    char *auto_name;
    DOF_TYPE *oldtype;
    BOOLEAN basflag = FALSE;

    MagicCheck(DOF, dest)

    if (dest != NULL && newtype != NULL && dest->type != newtype) {
	phgDofFree(&dest);
	dest = NULL;
    }

    dim = DofDim(src);
    /* the name of dest */
    if ((auto_name = (void *)name) == NULL) {
#if 0
	auto_name = phgAlloc(strlen(src->name) + 8 + 1);
	sprintf(auto_name, "copy of %s", src->name);
#else
	auto_name = strdup(src->name);
#endif
    }

    if (dest == NULL) {
	if (newtype == NULL)
	    newtype = src->type;
	dim_new = (newtype == NULL ? DofTypeDim(src) : newtype->dim);
	assert(dim % dim_new == 0);
	dest = phgDofNew_(g, newtype, newtype == NULL ? src->hp : NULL,
			  dim / dim_new, auto_name, DofNoAction,
			  srcfile, srcline);
	if (dest_ptr != NULL)
	    *dest_ptr = dest;
    }
    else {
	assert(dim == DofDim(dest));
	phgFree(dest->name); dest->name = strdup(auto_name);
	dest->srcfile = srcfile;
	dest->srcline = srcline;
	phgFree(dest->cache_func); dest->cache_func = NULL;
	newtype = dest->type;
	if (!SpecialDofType(newtype))
	    memset(dest->data, 0, DofGetDataCount(dest) * sizeof(*dest->data));
    }
    if (auto_name != name)
	phgFree(auto_name);

    phgDofClearCache(NULL, dest, NULL, NULL, FALSE);

    dest->DB_mask = src->DB_mask;
    if (src->DB_masks != NULL) {
	if (dest->DB_masks == NULL)
	    dest->DB_masks = phgAlloc(dest->dim * sizeof(*dest->DB_masks));
	memcpy(dest->DB_masks, src->DB_masks, dest->dim * sizeof(*dest->DB_masks));
    }
    dest->invariant = src->invariant;

    if (SpecialDofType(newtype)) {
	assert(newtype == src->type);
	dest->userfunc = src->userfunc;
	dest->userfunc_lambda = src->userfunc_lambda;
	if (newtype == DOF_CONSTANT)
	    memcpy(dest->data, src->data, dim * sizeof(*dest->data));
	return dest;
    }

    if ((newtype != NULL && newtype == src->type) ||
	(newtype == NULL && src->hp == dest->hp)) {
	/* simply duplicate the data */
	size_t size = DofGetDataCount(dest);
	if (size > 0)
	    memcpy(dest->data, src->data, sizeof(FLOAT) * size);
	dest->userfunc = src->userfunc;
	dest->userfunc_lambda = src->userfunc_lambda;
	return dest;
    }

    if (src->type == DOF_ANALYTIC) {
	if (src->userfunc_lambda != NULL)
	    phgDofSetDataByLambdaFunction(dest, src->userfunc_lambda);
	else
	    phgDofSetDataByFunction(dest, src->userfunc);
	return dest;
    }

    if (newtype != NULL && newtype->BasFuncs == NULL)
	phgError(1, "phgDofCopy: basis funcs for DOF type \"%s\" undefined.\n",
		 newtype->name);

    dest->userfunc = src->userfunc;
    dest->userfunc_lambda = src->userfunc_lambda;

    oldtype = src->type;
    if (oldtype == NULL)
	oldtype = src->hp->info->types[src->hp->info->min_order];

    if (!SpecialDofType(oldtype) && newtype != NULL && newtype->points != NULL
	&& !DofIsHP(src)) {
	basflag = TRUE;
	count = oldtype->nbas * oldtype->dim;
	basvalues = phgAlloc(newtype->nbas * count * sizeof(*basvalues));

	if (oldtype->invariant == TRUE)
	    get_bas_funcs(src, dest, src->g->roots, basvalues);
    }

    if (newtype == NULL)
	newtype = dest->hp->info->types[dest->hp->max_order];

    flags0 = phgCalloc((newtype->np_vert > 0 ? g->nvert : 0)
		       + (newtype->np_edge > 0 ? g->nedge : 0)
		       + (newtype->np_face > 0 ? g->nface : 0),
		       sizeof(*flags0));

    if (!SpecialDofType(oldtype) && oldtype->continuity < 0) {
	static DOF_TYPE DOF_WGTS = {DofCache,
	    "Weights", NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	    phgDofInitFuncPoint, NULL, NULL, NULL, FE_None,
	    FALSE, FALSE, -1, 0, 0, -1, 1, 0, 0, 0, 0
	};
	DOF_WGTS.np_vert = (newtype->np_vert > 0) ? 1 : 0;
	DOF_WGTS.np_edge = (newtype->np_edge > 0) ? 1 : 0;
	DOF_WGTS.np_face = (newtype->np_face > 0) ? 1 : 0;
	DOF_WGTS.nbas = DOF_WGTS.np_vert * NVert + DOF_WGTS.np_edge * NEdge +
	    DOF_WGTS.np_face * NFace;
	if (DOF_WGTS.nbas > 0) {
	    /* Other cases will be implemented later when really needed */
	    wgts = phgDofNew(g, &DOF_WGTS, 1, "weights", DofNoAction);
	    phgDofSetDataByValue(wgts, 0.0);
	    phgDofSetDataByValue(dest, 0.0);
	    /* allocate buffer for storing weighted vertex/edge/face data */
	    if ((n = DofGetDataCount(dest) - DofGetElementDataCount(dest)) > 0)
		buffer = phgCalloc(n, sizeof(*buffer));
	}
    }

    nvalues = dest->dim;
    cache_dof = src;
    bas_count = count;
    ForAllElements(g, e) {
	if (wgts != NULL) {
#if 0
	    /* use element volume as weight */
	    w = phgGeomGetVolume(g, e);
#else
	    /* use 1 as weight */
	    w = 1.0;
#endif
	}

	if (basflag && oldtype->invariant == FALSE)
	    get_bas_funcs(src, dest, e, basvalues);
	bas = basvalues;
	flags = flags0;

	if (DofIsHP(dest))
	    newtype = dest->hp->info->types[dest->hp->elem_order[e->index]];

	if (newtype->np_vert > 0) {
	    nd = nvalues * newtype->np_vert;
	    for (k = 0; k < NVert; k++) {
		if (flags[n = e->verts[k]] && wgts == NULL) {
		    /* note: count==0 and bas==NULL for variable order DOF */
		    bas += count * newtype->np_vert;
		    continue;
		}
		flags[n] = TRUE;
		a = DofVertexData(dest, n);
		newtype->InitFunc(dest, e, VERTEX, k, NULL, func, NULL, a,
					NULL);
		if (wgts != NULL) {
		    d = buffer + (a - DofData(dest));
		    for (i = 0; i < nd; i++)
			*(d++) += *(a++) * w;
		    *DofVertexData(wgts, n) += w;
		}
	    }
	    flags += g->nvert;
	}

	if (newtype->np_edge > 0) {
	    nd = nvalues * newtype->np_edge;
	    for (k = 0; k < NEdge; k++) {
		if (flags[n = e->edges[k]] && wgts == NULL) {
		    /* note: count==0 and bas==NULL for variable order DOF */
		    bas += count * newtype->np_edge;
		    continue;
		}
		flags[n] = TRUE;
		a = DofEdgeData(dest, n);
		newtype->InitFunc(dest, e, EDGE, k, NULL, func, NULL, a,
					NULL);
		if (wgts != NULL) {
		    d = buffer + (a - DofData(dest));
		    if (DofIsHP(dest))
			nd = dest->dim * (dest->hp->edge_index[n + 1] -
					  dest->hp->edge_index[n]);
		    for (i = 0; i < nd; i++)
			*(d++) += *(a++) * w;
		    *DofEdgeData(wgts, n) += w;
		}
	    }
	    flags += g->nedge;
	}

	if (newtype->np_face > 0) {
	    nd = nvalues * newtype->np_face;
	    for (k = 0; k < NFace; k++) {
		if (flags[n = e->faces[k]] && wgts == NULL) {
		    /* note: count==0 and bas==NULL for variable order DOF */
		    bas += count * newtype->np_face;
		    continue;
		}
		flags[n] = TRUE;
		a = DofFaceData(dest, n);
		newtype->InitFunc(dest, e, FACE, k, NULL, func, NULL, a,
					NULL);
		if (wgts != NULL) {
		    d = buffer + (a - DofData(dest));
		    if (DofIsHP(dest))
			nd = dest->dim * (dest->hp->face_index[n + 1] -
						 dest->hp->face_index[n]);
		    for (i = 0; i < nd; i++)
			*(d++) += *(a++) * w;
		    *DofFaceData(wgts, n) += w;
		}
	    }
	}

	if (newtype->np_elem > 0) {
	    a = DofElementData(dest, e->index);
	    newtype->InitFunc(dest, e, ELEMENT, 0, NULL, func, NULL, a,
					NULL);
	}
    }

    phgFree(basvalues);
    phgFree(flags0);

    if (wgts == NULL)
	return dest;

    if ((n = DofGetDataCount(dest) - DofGetElementDataCount(dest)) > 0) {
	memcpy(DofData(dest), buffer, n * sizeof(*buffer));
	phgFree(buffer);
    }

    if (DofIsHP(dest))
	newtype = dest->hp->info->types[dest->hp->max_order];

    a = DofData(dest);
    d = DofData(wgts);

#if USE_MPI
    /* FIXME: directly interacting with the vector assembly code in solver.c
       is more efficient */
    if (g->nprocs > 1) {
	SOLVER *solver = phgSolverCreate(SOLVER_BUILTIN, dest, NULL);
	INT K = 0, I;
	int j;

	if (newtype->np_vert > 0) {
	    nd = dest->count_vert;
	    for (n = 0; n < g->nvert; n++, d++) {
		if (g->types_vert[n] == UNREFERENCED) {
		    K += nd;
		    a += nd;
		    continue;
		}
		for (j = 0; j < nd; j++, K++, a++) {
		    I = phgSolverMapD2L(solver, 0, K);
		    assert(I >= 0 && I < solver->mat->rmap->localsize);
		    phgSolverAddRHSEntry(solver, I, *a);
		    phgSolverAddMatrixEntry(solver, I, I, *d);
		}
	    }
	}

	if (newtype->np_edge > 0) {
	    nd = nvalues * newtype->np_edge;
	    for (n = 0; n < g->nedge; n++, d++) {
		if (DofIsHP(dest))
		    nd = dest->dim * (dest->hp->edge_index[n + 1] -
					     dest->hp->edge_index[n]);
		if (g->types_edge[n] == UNREFERENCED) {
		    K += nd;
		    a += nd;
		    continue;
		}
		for (j = 0; j < nd; j++, K++, a++) {
		    I = phgSolverMapD2L(solver, 0, K);
		    assert(I >= 0 && I < solver->mat->rmap->localsize);
		    phgSolverAddRHSEntry(solver, I, *a);
		    phgSolverAddMatrixEntry(solver, I, I, *d);
		}
	    }
	}

	if (newtype->np_face > 0) {
	    nd = nvalues * newtype->np_face;
	    for (n = 0; n < g->nface; n++, d++) {
		if (DofIsHP(dest))
		    nd = dest->dim * (dest->hp->face_index[n + 1] -
					     dest->hp->face_index[n]);
		if (g->types_face[n] == UNREFERENCED) {
		    K += nd;
		    a += nd;
		    continue;
		}
		for (j = 0; j < nd; j++, K++, a++) {
		    I = phgSolverMapD2L(solver, 0, K);
		    assert(I >= 0 && I < solver->mat->rmap->localsize);
		    phgSolverAddRHSEntry(solver, I, *a);
		    phgSolverAddMatrixEntry(solver, I, I, *d);
		}
	    }
	}

	if (newtype->np_elem > 0) {
	    nd = nvalues * newtype->np_elem;
	    for (n = 0; n < g->nelem; n++) {
		if (DofIsHP(dest))
		    nd = dest->dim * (dest->hp->elem_index[n + 1] -
					     dest->hp->elem_index[n]);
		if (g->types_elem[n] == UNREFERENCED) {
		    K += nd;
		    a += nd;
		    continue;
		}
		for (j = 0; j < nd; j++, K++, a++) {
		    I = phgSolverMapD2L(solver, 0, K);
		    assert(I >= 0 && I < solver->mat->rmap->localsize);
		    phgSolverAddRHSEntry(solver, I, *a);
		    phgSolverAddMatrixEntry(solver, I, I, 1.0);
		}
	    }
	}

	phgSolverSolve(solver, TRUE, dest, NULL);
	phgSolverDestroy(&solver);

	phgDofFree(&wgts);

	return dest;
    }
#endif

    if (newtype->np_vert > 0) {
	k = nvalues * newtype->np_vert;
	for (n = 0; n < g->nvert; n++) {
	    if ((w = *(d++)) == 0.0) {
		a += k;
	    }
	    else if (k == 1) {
		*(a++) /= w;
	    }
	    else {
		w = 1.0 / w;
		for (i = 0; i < k; i++)
		    *(a++) *= w;
	    }
	}
    }

    if (newtype->np_edge > 0) {
	k = nvalues * newtype->np_edge;
	for (n = 0; n < g->nedge; n++) {
	    if (DofIsHP(dest))
		k = dest->dim * (dest->hp->edge_index[n + 1] -
					dest->hp->edge_index[n]);
	    if ((w = *(d++)) == 0.0) {
		a += k;
	    }
	    else if (k == 1) {
		*(a++) /= w;
	    }
	    else {
		w = 1.0 / w;
		for (i = 0; i < k; i++)
		    *(a++) *= w;
	    }
	}
    }

    if (newtype->np_face > 0) {
	k = nvalues * newtype->np_face;
	for (n = 0; n < g->nface; n++) {
	    if (DofIsHP(dest))
		k = dest->dim * (dest->hp->face_index[n + 1] -
					dest->hp->face_index[n]);
	    if ((w = *(d++)) == 0.0) {
		a += k;
	    }
	    else if (k == 1) {
		*(a++) /= w;
	    }
	    else {
		w = 1.0 / w;
		for (i = 0; i < k; i++)
		    *(a++) *= w;
	    }
	}
    }

    phgDofFree(&wgts);

    return dest;
}
