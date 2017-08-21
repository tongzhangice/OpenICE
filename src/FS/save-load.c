/*
 * Save & load DOF data, only for nodal DOF type.
 *
 * Usage:
 *
 * 1. Save: 
 * 1.1 phgExportMedit2: save mesh, vert coord,
 * 1.2 save_map       : save elem to vert edge face map
 * 1.3 save_dof_data  : save dof data
 *
 * 2. Load:
 * 2.0 Use exteranl program to generate text mesh file with binary
 *     coord file. 
 * 2.1 save_element_id: save element index, before redistribute
 * 2.2 load_dof_data2 : load dof date, after redistribute
 * 2.3 free_element_id
 *
 * Note:
 * Load is essentially a sequtial process, the time is mostly spend
 * on distribute mesh.
 * Using load_dof_data2, no dof data distribute is not needed, still,
 * further improvement could be made to read in only locate map and 
 * dof data. 
 *
 *  */

#include "phg.h"
#include "io.h"
#include <string.h>
#include <math.h>

#define DofVertData(dof, i) DofVertexData(dof, i)
#define GlobalVert(g,no) GlobalVertex(g, no)
//#define GlobalElem(g,no) GlobalElement(g, no)
#define NElem 0

#define GetElementVertGMap(e, M)			\
    GetElementVertices(e, 0, M[0], M[1], M[2], M[3]);
#define GetElementFaceGMap(e, M)			\
    GetElementVertices(e, 0, M[0], M[1], M[2], M[3]);

static DOF *elem_id = NULL;

static void 
GetElementEdgeGMap(SIMPLEX *e, int *M)
{
    int i, v[2], N[NVert]; 
    
    GetElementVertGMap(e, N);
    bzero(M, NEdge * sizeof(*M));
    for (i = 0; i < NEdge; i++) {
	v[0] = GetEdgeVertex(i, 0);
	v[1] = GetEdgeVertex(i, 1);
	v[0] = N[v[0]];
	v[1] = N[v[1]];
	M[i] = GetEdgeNo(v[0], v[1]);
    }
    //SHOW_V(M, NFace);
    return;
}


void 
save_map(GRID *g, const char *file) {
    SIMPLEX *e;
    int i, j, n, nEmap = NVert + NEdge + NFace + NElem;
    FILE *fp = NULL;
    MPI_Status status;
    //FLOAT *fbuffer;
    INT *ibuffer;
# define TRUNK_SIZE 65536

    /* save E2Gmap */
    if (g->rank == 0) {
	char e2g_file[1000];
	sprintf(e2g_file, "%s.map", file);
	if ((fp = fopen(e2g_file, "w+t")) == NULL) {
	    phgError(1, "cannot open output file \"%s\".\n", e2g_file);
	    return;
	}

	fwrite(&g->nelem_global, sizeof(INT), 1, fp);
	ForAllElements(g, e) {
	    int M[6];
#define WRITE_MAP(GType, gtype)					\
	    GetElement##GType##GMap(e, M);			\
	    for (i = 0; i < N##GType; i++) {			\
		INT ind = Global##GType(g, e->gtype##s[M[i]]);	\
		fwrite(&ind, sizeof(INT), 1, fp);		\
	    }

	    WRITE_MAP(Vert, vert);
	    WRITE_MAP(Edge, edge);
	    WRITE_MAP(Face, face);
	    //WRITE_MAP(Elem, elem);
#undef WRITE_MAP
	}
	/* receive elements from other processes */
	ibuffer = phgAlloc(TRUNK_SIZE * (nEmap) * sizeof(*ibuffer));
	for (j = 1; j < g->nprocs; j++) {
	    while (TRUE) {
		MPI_Probe(j, 222, g->comm, &status);
		MPI_Get_count(&status, PHG_MPI_INT, &n);
		assert(n <= ((nEmap) * TRUNK_SIZE) && n % (nEmap) == 0);
		MPI_Recv(ibuffer, n, PHG_MPI_INT, j, 222, g->comm, &status);
		/* process received vertices */
		for (i = 0; i < n; i += nEmap) {
		    fwrite(ibuffer + i, sizeof(INT), nEmap, fp);
		}
		if (n < (nEmap) * TRUNK_SIZE)
		    break;
	    }
	}
	phgFree(ibuffer);
	fclose(fp);
    }
    else {
	/* send elements to root process */
	ibuffer = phgAlloc(TRUNK_SIZE * (nEmap) * sizeof(*ibuffer));
	n = 0;
	ForAllElements(g, e) {
	    int k = 0, M[6];
#define WRITE_MAP(GType, gtype)					\
	    GetElement##GType##GMap(e, M);			\
	    for (i = 0; i < N##GType; i++) {			\
		ibuffer[(nEmap) * n + k++] =			\
		    Global##GType(g, e->gtype##s[M[i]]);	\
	    }

	    WRITE_MAP(Vert, vert);
	    WRITE_MAP(Edge, edge);
	    WRITE_MAP(Face, face);
	    //WRITE_MAP(Elem, elem);
#undef WRITE_MAP
	    if (++n >= TRUNK_SIZE) {
		MPI_Send(ibuffer, (nEmap) * n, PHG_MPI_INT, 0, 222, g->comm);
		n = 0;
	    }
	}
	/* send the last block (may be size 0), which also marks EOD */
	MPI_Send(ibuffer, (nEmap) * n, PHG_MPI_INT, 0, 222, g->comm);
	phgFree(ibuffer);
    }

    return;
}

void 
save_dof_data(GRID *g, DOF *dof, const char *file)
{
    FLOAT *buffer = NULL;
    size_t buffer_size = 0;
    INT ndata;

    if (phgIOOpen(g, file) == FALSE)
	phgError(1, "Export Dof data file %s failed!", file);
    buffer_size += DofGetDataCount(dof);
    buffer = phgCalloc(buffer_size,  sizeof(FLOAT));

#define WRITE_DATA(dof, GType, gtype)					\
    if (dof->type->np_##gtype > 0){					\
	ndata = g->n##gtype * dof->count_##gtype;			\
	memcpy(buffer, dof->data_##gtype, ndata * sizeof(FLOAT));	\
	phgIOWrite(g, buffer, sizeof(FLOAT) * dof->count_##gtype,	\
		   g->n##gtype, g->n##gtype##_global,			\
		   g->types_##gtype, g->L2Gmap_##gtype, TRUE);		\
    }

    WRITE_DATA(dof, Vertex, vert);
    WRITE_DATA(dof, Edge, edge);
    WRITE_DATA(dof, Face, face);
    //WRITE_DATA(dof, Elem, elem);
#undef WRITE_DATA
    phgFree(buffer);
    phgIOClose();

    return;
}


void 
load_dof_data(GRID *g, DOF *dof, const char *data_file, const char *mesh_file)
/* Load dof data BEFORE mesh redistribute. The loading process is sequential.
 *
 * E2Gmap is read element wize only on root proc.
 *  */
{
    SIMPLEX *e;
    FLOAT *buffer = NULL, *buf_vert, *buf_edge,	*buf_face, *buf_elem;
    int dim = dof->dim;
    size_t buffer_size = 0;
    int i, j, k, n, count, nEmap = NVert + NEdge + NFace + NElem;
    INT ndata, nelem, *E2Gmap = NULL;
    FILE *fp;
    char e2g_file[1000];

    if (phgRank == 0) {
	if ((fp = fopen(data_file, "r")) == NULL)
	    phgError(1, "read Dof data %s failed!\n", data_file);
	buffer_size += DofGetDataCount(dof);
	buffer = phgCalloc(buffer_size, sizeof(FLOAT));
	buf_vert = buffer;
	buf_edge = buf_vert + g->nvert * dof->count_vert;
	buf_face = buf_edge + g->nedge * dof->count_edge;
	buf_elem = buf_face + g->nface * dof->count_face;

#define READ_DATA(dof, GType, gtype)				\
	if (dof->type->np_##gtype > 0){				\
	    int n##gtype = g->n##gtype;			\
	    ndata = n##gtype * dof->count_##gtype;		\
	    n = fread(buf_##gtype, ndata, sizeof(FLOAT), fp);	\
	}

	READ_DATA(dof, Vertex, vert);
	READ_DATA(dof, Edge, edge);
	READ_DATA(dof, Face, face);
	//READ_DATA(dof, Elem, elem);
#undef READ_DATA
	fclose(fp);

	/* read E2Gmap */
	sprintf(e2g_file, "%s.map", mesh_file);
	if ((fp = fopen(e2g_file, "r+t")) == NULL) {
	    phgError(1, "cannot open output file \"%s\".\n", e2g_file);
	    return;
	}

	E2Gmap = phgCalloc(nEmap, sizeof(*E2Gmap));
	n = fread(&nelem, sizeof(INT), 1, fp);
	assert(g->nelem == nelem);

	/* Note:
	 * element order of E2Gmap is the same as vector g->roots, but NOT
	 * the same as g->elements, as well as ForAllElements.
	 */
	e = g->roots;
	for (j = 0; j < g->nelem; j++, e++) {
	    FLOAT *val;
	    int M[6], shift = 0;
	    
	    n = fread(E2Gmap, sizeof(INT), nEmap, fp);

	    /* Vertex data */
	    GetElementVertGMap(e, M);
	    for (i = 0; i < NVert; i++) {
		val = DofVertexData(dof, GlobalVertex(g, e->verts[M[i]]));
		for (k = 0; k < dim; k++)
		    val[k] = buf_vert[E2Gmap[i]*dim + k];
		//printf("E:%d, vert[%d]:%d\n", e->index, i, e->verts[i]);
	    }

	    /* Edge data */
	    shift += NVert;
	    count = dof->count_edge;
	    GetElementEdgeGMap(e, M);
	    for (i = 0; i < NEdge; i++) {
		val = DofEdgeData(dof, GlobalEdge(g, e->edges[M[i]]));
		for (k = 0; k < count; k++)
		    val[k] = buf_edge[E2Gmap[i+shift]*count + k];
	    }
		
	    /* Face data */
	    shift += NEdge;
	    count = dof->count_face;
	    GetElementFaceGMap(e, M);
	    for (i = 0; i < NFace; i++) {
		val = DofFaceData(dof, GlobalFace(g, e->faces[M[i]]));
		for (k = 0; k < count; k++) 
		    val[k] = buf_face[E2Gmap[i+shift]*count + k];
	    }
	}

	fclose(fp);
	phgFree(buffer);
    }
    return;
}

void 
load_dof_data2(GRID *g, DOF *dof, const char *data_file, const char *mesh_file)
/* Load dof data AFTER mesh redistribute. The loading process is parallel.
 *
 * E2Gmap is read as a whole on each proc.
 * TODO: read only the local part of E2Gmap and local part of dof data.
 *  */
{
    SIMPLEX *e;
    FLOAT *buffer = NULL, *buf_vert, *buf_edge,	*buf_face, *buf_elem;
    int dim = dof->dim;
    size_t buffer_size = 0;
    int i, k, n, count, nEmap = NVert + NEdge + NFace + NElem;
    INT ndata, nelem_global, *E2Gmap = NULL;
    FILE *fp;
    char e2g_file[1000];

    if ((fp = fopen(data_file, "r")) == NULL)
	phgError(1, "read Dof data %s failed!\n", data_file);
    buffer_size += DofGetDataCountGlobal(dof);
    buffer = phgCalloc(buffer_size, sizeof(FLOAT));
    buf_vert = buffer;
    buf_edge = buf_vert + g->nvert_global * dof->count_vert;
    buf_face = buf_edge + g->nedge_global * dof->count_edge;
    buf_elem = buf_face + g->nface_global * dof->count_face;
    phgPrintf("* Load dof: load datafile: %0.4lfMB\n", 
	      phgMemoryUsage(g, NULL) / (1024.0 * 1024.0));

#define READ_DATA(dof, GType, gtype)				\
    if (dof->type->np_##gtype > 0){				\
	int n##gtype_global = g->n##gtype##_global;		\
	ndata = n##gtype_global * dof->count_##gtype;		\
	n = fread(buf_##gtype, ndata, sizeof(FLOAT), fp);	\
    }

    READ_DATA(dof, Vertex, vert);
    READ_DATA(dof, Edge, edge);
    READ_DATA(dof, Face, face);
    //READ_DATA(dof, Elem, elem);
#undef READ_DATA
    fclose(fp);

    /* read E2Gmap */
    sprintf(e2g_file, "%s.map", mesh_file);
    if ((fp = fopen(e2g_file, "r+t")) == NULL) {
	phgError(1, "cannot open output file \"%s\".\n", e2g_file);
	return;
    }
    E2Gmap = phgCalloc(nEmap * g->nelem_global, sizeof(*E2Gmap));
    phgPrintf("* Load dof: load E2Gmap: %0.4lfMB\n", 
	      phgMemoryUsage(g, NULL) / (1024.0 * 1024.0));
    n = fread(&nelem_global, sizeof(INT), 1, fp);
    n = fread(E2Gmap, sizeof(INT), nEmap * g->nelem_global, fp);
    assert(g->nelem_global == nelem_global);
    ForAllElements(g, e) {
	FLOAT *val;
	INT *e2g_map, M[6], shift = 0, e_id;

	e_id = (INT)(*DofElementData(elem_id, e->index));
	e2g_map = E2Gmap + e_id * nEmap;

	/* Vertex data */
	GetElementVertGMap(e, M);
	count = dof->count_vert;
	for (i = 0; i < NVert; i++) {
	    val = DofVertexData(dof, e->verts[M[i]]);
	    for (k = 0; k < dim; k++)
		val[k] = buf_vert[e2g_map[i]*count + k];
	    //printf("E:%d, vert[%d]:%d\n", e->index, i, e->verts[i]);
	}

	/* Edge data */
	shift += NVert;
	count = dof->count_edge;
	GetElementEdgeGMap(e, M);
	for (i = 0; i < NEdge; i++) {
	    val = DofEdgeData(dof, e->edges[M[i]]);
	    for (k = 0; k < count; k++)
		val[k] = buf_edge[e2g_map[i+shift]*count + k];
	}
		
	/* Face data */
	shift += NEdge;
	count = dof->count_face;
	GetElementFaceGMap(e, M);
	for (i = 0; i < NFace; i++) {
	    val = DofFaceData(dof, e->faces[M[i]]);
	    for (k = 0; k < count; k++) 
		val[k] = buf_face[e2g_map[i+shift]*count + k];
	}
    }

    fclose(fp);
    phgFree(buffer);
    phgFree(E2Gmap);
    return;
}

void 
save_dof_data3(GRID *g, DOF *dof, const char *file)
/*
 * Each proc save it's own Dof data
 *  */
{
    FILE *fp = NULL;
    char fname[100];
    INT n = DofGetDataCount(dof);

    sprintf(fname, "%s.p%03d", file, g->rank);
    phgInfo(1, "* save data to %s\n", fname);
    if ((fp = fopen(fname, "w")) == NULL) {
	phgError(1, "read Dof data %s failed!\n", fname);
    } else {
	fwrite(&n, sizeof(INT), 1, fp);
	fwrite(dof->data, sizeof(FLOAT), n, fp);
	fclose(fp);
    }

    return;
}

void 
load_dof_data3(GRID *g, DOF *dof, const char *file, const char *mesh_file)
/*
 * Each proc read it's own Dof data
 *  */
{
    FILE *fp = NULL;
    char fname[100];
    INT n0 = DofGetDataCount(dof), n;
    int nr;

    Unused(mesh_file);
    sprintf(fname, "%s.p%03d", file, g->rank);
    phgInfo(1, "* read data from %s\n", fname);
    if ((fp = fopen(fname, "r")) == NULL) {
	phgError(1, "read Dof data %s failed!\n", fname);
    } else {
	nr = fread(&n, sizeof(INT), 1, fp);
	if (n0 != n) 
	    phgError(1, "read [%d] but saved [%d] !\n", n0, n);
	nr = fread(dof->data, sizeof(FLOAT), n, fp);
	fclose(fp);
    }

    return;
}


void
save_element_id(GRID *g)
{
    SIMPLEX *e;
    INT i;
    FLOAT *v;

    elem_id = phgDofNew(g, DOF_P0, 1, "element_id", DofInterpolation);
    if (phgRank == 0) {
	e = g->roots;
	for (i = 0; i < g->nelem; i++, e++) {
	    v = DofElementData(elem_id, e->index);
	    *v = i;
	}
    }
    return;
}

void 
free_element_id(GRID *g)
{
    phgDofFree(&elem_id);
    elem_id = NULL;
}


void
phgResumeStage(GRID *g, FLOAT *time, int *tstep, char *mesh_file, char *data_file)
/* Load resume info (time, tstep, mesh file name, data file name)
 * from file "resume.log".
 *
 * Note: resume should be execute by ALL ranks.
 * */
{

    FILE *fp;
    FLOAT time0;
    int tstep0, n;
    char mesh0[100];
    char data0[100];

    if ((fp = fopen("resume.log", "r")) == NULL)
	phgError(1, "read resume.log failed!\n");
    n = fscanf(fp, "%lf", &time0);
    n = fscanf(fp, "%d", &tstep0);
    n = fscanf(fp, "%s", mesh0);
    n = fscanf(fp, "%s", data0);
    n = fclose(fp);
    if (time != NULL)
	*time = time0;
    if (tstep != NULL)
	*tstep = tstep0;
    if (mesh_file != NULL)
	strcpy(mesh_file, mesh0);
    if (data_file != NULL)
	strcpy(data_file, data0);

    return;
}


/*******************/
/* Resume routines */
/*******************/
void
phgResumeLogUpdate(GRID *g, FLOAT *time, int *tstep, char *mesh_file, char *data_file)
/* Save resume info (time, tstep, mesh file name, data file name)
 * to file "resume.log".
 *
 * Resume info is stored as static variable, component is not updated
 * if input points is NULL.
 * */
{
    static FLOAT time0;
    static int tstep0;
    static char mesh0[100];
    static char data0[100];
    FILE *fp;

    if (g->rank != 0)
	return;

    if (time != NULL)
	time0 = *time;
    if (tstep != NULL)
	tstep0 = *tstep;
    if (mesh_file != NULL)
	strcpy(mesh0, mesh_file);
    if (data_file != NULL)
	strcpy(data0, data_file);

    if ((fp = fopen("resume.log", "w")) == NULL)
	phgError(1, "open resume log file failed!\n");
    fprintf(fp, "%24.12E\n%d\n%s\n%s\n", time0, tstep0, mesh0, data0);
    fclose(fp);

    return;
}

