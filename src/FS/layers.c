/*
 * Ice sheet - layers module.
 * Layers info is maintained for computaion involving surface-bottom. 
 *
 *
 * Prerequsite:
 * 1. 2D vert to 3D vert maping,
 *    L2Gmap_vert, 
 *    Z direction vert connection.
 * 2. 2D triangle to vert connection.
 * 3. 3D partation is based on 2D partation.
 *
 * Layers info:
 * 1. 2D vert to vert chain in Z dir.
 * 2. 2D triangle to elem chain in Z dir.
 * 3. Bottom Dof to Dof chain in Z difr. 
 *
 *
 *  */
#include "ins.h"
#include "layers.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "vtk-draw.h"
 
#ifdef VTK_DEBUG
#else
int vtk_verb = 3;
char vtk_tmp_str[1000];
#endif



#define VTK_VERB(VERB) {				\
	verb = (VERB < vtk_verb) ? vtk_verb : VERB;	\
	phgInfo(0, "layers verb: %d\n", verb);		\
    }


/******************/
/* GEO DEFINATION */
/******************/
#define phgElemInfoCount 6

typedef enum {
    QUADRILATERAL = 0,
    TRIANGLE	= 1,
    TETRA	= 2,
    PYRAMID	= 3,
    PRISM	= 4,
    BRICK	= 5
} ELEM_TYPE;

typedef short EDGE_INFO[2];	/* list of vertices */
typedef short FACE_INFO[5];	/* 0: # of vertices, 1..: list of vertices */
typedef void (*ADD_TETRA_FUNC)(int v0, int v1, int v2, int v3, int bound_type[]);
typedef void (*FUNC2TET)(int verts[], ADD_TETRA_FUNC func, int bound_type[]);

typedef struct {
    const char  *name;		/* name of the element type */
    EDGE_INFO   *edge_info;	/* list of edges */
    FACE_INFO   *face_info;	/* list of faces */
    short       nvert;		/* number of vertices */
    short	nedge;		/* number of edges */
    short       nface;		/* number of faces */
    FUNC2TET    func2tet;	/* function to covert to tet */
} ELEM_INFO;


/*****************/
/* Reading Utils */
/*****************/
#define READ_NUMBER						\
    if (!get_token(fp, token)) strcpy(token, "End ALL");	\
    if (isalpha((int)(token[0]))) {				\
	phgWarning("fewer entries (%d) than expected.\n", i);	\
	break;							\
    }

#define GET_TOKEN {					\
	if (!get_token(fp, token)) {			\
	    phgPrintf("error on line: %d\n", __LINE__);	\
	    goto error;					\
	}						\
    }

#define UNUSED_LINE { 				\
	fgets(token, 1024, fp);			\
    }
//	fprintf(stderr, "Unuesd: %s", token);	
#define NEXT_LINE UNUSED_LINE

#define UNUSED_LINE_CHECK(str) {				\
	fgets(token, 1024, fp);					\
	if (!strcmp(token, str))				\
	    phgError(-1, "Unmatched phase when read mesh file,"	\
		     " line:%d!\n", __LINE__);			\
    }
//	fprintf(stderr, "Unuesd: %s", token);	


static ELEM_INFO phgElemInfo[] = {
    /* name	list of edges	list of faces	nvert	nedge	nface  func2tet*/
    {"quad",	NULL,	NULL,	4,	4,	4, NULL},
    {"tria",	NULL,	NULL,	3,	3,	3, NULL},
};

static int Gambit_Elemtype[phgElemInfoCount] = {
    2,				/* Quadrilateral */
    3,				/* Triangle */
    6, 				/* Tetrahedra */
    7, 				/* Pyramid */
    5, 				/* Wedge */
    4				/* Brick */
};

static int Gambit_Vertmap[phgElemInfoCount][100] = {
    {0, 1, 2, 3},		/* Quadrilateral */
    {0, 1, 2},			/* Triangle */
    {0, 1, 2, 3},               /* Tetrahedra */
    {0, 1, 3, 2, 4}, 		/* Pyramid */
    {0, 1, 2, 3, 4, 5},		/* Wedge */
    {0, 1, 3, 2, 4, 5, 7, 6}	/* Brick */
};

static int Gambit_Facemap[phgElemInfoCount][100] = {
    {0, 1, 2, 3},		/* Quadrilateral */
    {0, 1, 2},			/* Triangle */
    {3, 2, 0, 1},               /* Tetrahedra */
    {4, 0, 1, 2, 3}, 		/* Pyramid */
    {0, 1, 2, 3, 4, 5},		/* Wedge */
    {0, 1, 2, 3, 4, 5, 6, 7}	/* Brick */
};


static BOOLEAN
get_token(FILE *fp, char *token)
{
    int c;
    char *p;

    while (TRUE) {
	if (fscanf(fp, "%s", token) != 1)
	    return FALSE;
	if (token[0] != '#')
	    break;
	/* skip to newline */
	do
	    if ((c = fgetc(fp)) == EOF)
		return FALSE;
	while (c != '\n');
    }
    if ((p = strchr(token, '#')) != NULL)
	*p = '\0';
    return TRUE;
}

#if 0
/*
 * For gambit format.
 *  */
LAYERED_MESH *
import_layered_mesh(char *file_tira,
		    char *file_layer,
		    char *file_nodeZ,
		    int nproc)
{
    LAYERED_MESH *gL;
    FILE *fp;
    char token[1024]; 
    INT i, j, n, k, Nly;
    INT *vtypes = NULL;			/* vertex types */
    INT nlistt = 0;			/* number of triangles */
    INT nlistE[phgElemInfoCount], size_listE[phgElemInfoCount],
	nlistF[phgElemInfoCount], size_listF[phgElemInfoCount];
    void *listE[phgElemInfoCount];
    INT (*listF[phgElemInfoCount])[3];
    INT (*listEtype)[2] = NULL;	/* list of all element type,
				 * and idx of that type */
    INT *nbsetss;
    INT numnp, nelem, ngrps, nbsets;
    FLOAT *R;
    int verb;

    PHG_CALLOC(gL, 1);

    Unused(vtypes);
    Unused(Gambit_Vertmap);
    bzero(nlistE, sizeof(nlistE));
    bzero(size_listE, sizeof(size_listE));
    bzero(listE, sizeof(listE));
    bzero(nlistF, sizeof(nlistF));
    bzero(size_listF, sizeof(size_listF));
    bzero(listF, sizeof(listF));

    /***********************/
    /* Control Information */
    /***********************/
    fp = fopen(file_tira, "r");
    assert(fp != NULL);

    UNUSED_LINE_CHECK("        CONTROL INFO 2.3.16");
    UNUSED_LINE_CHECK("** GAMBIT NEUTRAL FILE");
    UNUSED_LINE; /* gambit session id */
    if (!get_token(fp, token) || strcmp(token, "PROGRAM:") ||
	!get_token(fp, token) || strcmp(token, "Gambit") ||
	!get_token(fp, token) || strcmp(token, "VERSION:") ||
	!get_token(fp, token) || strcmp(token, "2.3.16")) {
    error:
	phgError(1, "invalid or unsupported input file, abort.\n");
    }
    NEXT_LINE; 
    UNUSED_LINE; /* time when gambit file is created */
    UNUSED_LINE_CHECK("     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL");
    GET_TOKEN; numnp = atoi(token);
    GET_TOKEN; nelem = atoi(token);
    GET_TOKEN; ngrps = atoi(token);
    GET_TOKEN; nbsets = atoi(token);
    gL->nvert = numnp;
    gL->ntria = nelem;
    if (!get_token(fp, token) || strcmp(token, "2") ||
	!get_token(fp, token) || strcmp(token, "2"))
	goto error;     /* 2D */
    NEXT_LINE; 
    UNUSED_LINE_CHECK("ENDOFSECTION");

    /*********************/
    /* Nodal Coordinates */
    /*********************/
    UNUSED_LINE_CHECK("   NODAL COORDINATES 2.3.16");
    n = numnp;
    phgInfo(2, "number of vertices: %d\n", n);
    PHG_CALLOC(gL->verts, n);
    for (i = 0; i < n; i++) {
	READ_NUMBER;
	GET_TOKEN; gL->verts[i][0] = atof(token);
	GET_TOKEN; gL->verts[i][1] = atof(token);
	gL->verts[i][2] = -0.1;	/* base */
	gL->verts[i][3] = 1;	/* surf */
    }
    NEXT_LINE;
    UNUSED_LINE_CHECK("ENDOFSECTION"); 

    /*****************************/
    /* Element/Cell Connectivity */
    /*****************************/
    UNUSED_LINE_CHECK("      ELEMENTS/CELLS 2.3.16");
    n = nelem;
    phgInfo(2, "number of elements: %d\n", n);
    PHG_CALLOC(listEtype, n);
    PHG_CALLOC(gL->trias, n);
    for (i = 0; i < n; i++) {
	int type, ndp, elem_type, nvert_type, nface_type;
	/* element index */
	READ_NUMBER;	
	/* element type */
	if (!get_token(fp, token)) 
	    goto error;
	type = atoi(token);
	/* element ndp */
	if (!get_token(fp, token)) 
	    goto error;
	ndp = atoi(token);

	for (j = 0; j < phgElemInfoCount; j++) 
	    if (type == Gambit_Elemtype[j]) 
		break;
	assert(j < phgElemInfoCount);
	elem_type = j;
	assert(j == TRIANGLE);		/* only triangle */
	
	nvert_type = phgElemInfo[elem_type].nvert;
	nface_type = phgElemInfo[elem_type].nface;
	assert(ndp == nvert_type);
	if (elem_type == TRIANGLE) {
	    GET_TOKEN; gL->trias[nlistt].verts[0] = atoi(token) -1;
	    GET_TOKEN; gL->trias[nlistt].verts[1] = atoi(token) -1;
	    GET_TOKEN; gL->trias[nlistt].verts[2] = atoi(token) -1;
	    listEtype[i][0] = elem_type;
	    listEtype[i][1] = nlistt;
	    nlistt++;
	} else 
	    abort();		/* only triagles */
    }
    NEXT_LINE;
    UNUSED_LINE_CHECK("ENDOFSECTION"); 



    /*****************************/
    /* Element Group Information */
    /*****************************/
    for (j = 0; j < ngrps; j++) {
	int nflags;
	UNUSED_LINE_CHECK("       ELEMENT GROUP 2.3.16");
	/* group index */
	if (!get_token(fp, token) || strcmp(token, "GROUP:") ||
	    !get_token(fp, token))
	    goto error;
	assert(j == atoi(token) - 1); 
	/* num of group record */
	if (!get_token(fp, token) || strcmp(token, "ELEMENTS:") ||
	    !get_token(fp, token))
	    goto error;
	n = atoi(token); 	
	/* materials and flags */
	if (!get_token(fp, token) || strcmp(token, "MATERIAL:") ||
	    !get_token(fp, token) || 
	    !get_token(fp, token) || strcmp(token, "NFLAGS:") ||
	    !get_token(fp, token))
	    goto error;
	nflags = atoi(token);
	NEXT_LINE;
	
	/* group name */
	UNUSED_LINE; 		
	/* solver dependent flags */
	for (i = 0; i < nflags; i++)
	    UNUSED_LINE; 	
	
	/* element in group */
	for (i = 0; i < n; i++) {
	    INT elem, elem_type, elem_idx;
   
	    GET_TOKEN;
	    elem = atoi(token) - 1;
	    elem_type = listEtype[elem][0];
	    elem_idx = listEtype[elem][1];
	    
	    if (elem_type == TRIANGLE) {
		gL->trias[elem_idx].region_mark = j;
	    } else
		abort();
	}
	NEXT_LINE;
	UNUSED_LINE_CHECK("ENDOFSECTION"); 
    }


    /****************************/
    /* Boundary Conditions Sets */
    /****************************/
    nbsetss = phgCalloc(nbsets, sizeof(*nbsetss));
    for (j = 0; j < nbsets; j++) {
	UNUSED_LINE_CHECK(" BOUNDARY CONDITIONS 2.3.16");
	/* bdry cond name */
	GET_TOKEN;
	phgInfo(2, "boundary condition %s:", token);
	/* bdry cond type */
	GET_TOKEN;
	/* Note: elem:face type is allowed,
	 *       vert type is NOT allowed! */
	assert(1 == atoi(token));
	/* num of bdry cond record */
	GET_TOKEN;
	n = atoi(token);
	phgInfo(2, "%d\n", n);
	/* bdry cond: value and code, unused */
	UNUSED_LINE;
	
	nbsetss[j] = n;
	for (i = 0; i < n; i++) {
	    int elem, face, elem_type;
	    INT (*listf)[3], nlistf;
	    /* element */
	    READ_NUMBER;
	    elem = atoi(token) - 1; 
	    /* element type */
	    GET_TOKEN;
	    /* bdry face */
	    GET_TOKEN;
	    face = atoi(token) - 1;

	    elem_type = listEtype[elem][0];
	    listf = listF[elem_type];
	    nlistf = nlistF[elem_type];
	    assert(elem_type < phgElemInfoCount);

	    if (nlistf >= size_listF[elem_type]) {
		int old_size = size_listF[elem_type];
		size_listF[elem_type] += 1024;
		listF[elem_type] = listf = 
		    phgRealloc_(listf,
				size_listF[elem_type] * sizeof(*listf),
				old_size * sizeof(*listf));
	    }

	    listf[nlistf][0] = listEtype[elem][1]; 
	    listf[nlistf][1] = Gambit_Facemap[elem_type][face];
	    listf[nlistf][2] = j;
	    nlistF[elem_type]++;
	    
	    {
		int e = listf[nlistf][0];
		int f = listf[nlistf][1];
		int b = listf[nlistf][2];
		gL->trias[e].bound_types[f] = b;
	    }
	}
	NEXT_LINE;
	UNUSED_LINE_CHECK("ENDOFSECTION"); 
    }
    fclose(fp);
    

    /**************/
    /* Layer Info */
    /**************/
    fp = fopen(file_layer, "r");
    assert(fp != NULL);
    fscanf(fp, "%d", &Nly);
    phgInfo(0, "Layer: %d\n", Nly);
    R = calloc(Nly+1, sizeof(*R));
    for (i = 0; i < Nly+1; i++) {
	fscanf(fp, "%lf", R + i);
	phgInfo(0, "   L[%2d]: %f\n", i, R[i]);
    }
    fclose(fp);
    gL->max_nlayer = Nly;
    phgFree(R);


    /******************/
    /* Vert List Info */
    /******************/
    fp = fopen(file_nodeZ, "r");
    PHG_CALLOC(gL->vert_global_lists, numnp);
    for (i = 0; i < numnp; i++) {
	int nv, id;
	fscanf(fp, "%d", &id);
        /* read 2D node id */
	assert(i == id-1);
	fscanf(fp, "%d", &nv);
        /* read # of layers */
	PHG_CALLOC(gL->vert_global_lists[i], nv + 1);
	/* fisrt as num of verts in list */
	gL->vert_global_lists[i][0] = nv; 

	for (j = 0; j < nv; j++) {
	    fscanf(fp, "%d", &id);
	    gL->vert_global_lists[i][j+1] = id -1;
	}
    }
    fclose(fp);


    /**************/
    /* Debug plot */
    /**************/
    VTK_VERB(3);
    vtkSetColor(verb, "red");
    for (i = 0; i < gL->nvert; i++) {
#if 0
	sprintf(vtk_tmp_str, "     -%d", i);
#else
	sprintf(vtk_tmp_str, "     -%d",
		gL->vert_global_lists[i][1]);
#endif
	vtkDrawPoint(verb, gL->verts[i], vtk_tmp_string);
    }
    vtkSetColor(verb, "white");
    TRIA *t = gL->trias;
    for (i = 0; i < gL->ntria; i++, t++) {
	FLOAT *v0, *v1, *v2;
	v0 = gL->verts[t->verts[0]];
	v1 = gL->verts[t->verts[1]];
	v2 = gL->verts[t->verts[2]];
	vtkDrawLine(verb, v0, v1);
	vtkDrawLine(verb, v0, v2);
	vtkDrawLine(verb, v1, v2);
    }
    vtkPause(verb);

    for (j = 0; j < phgElemInfoCount; j++) {
	if (listE[j] != NULL) {
	    phgFree(listE[j]);
	    listE[j] = NULL;
	    nlistE[j] = 0;
	}
	if (listF[j] != NULL) {
	    phgFree(listF[j]);
	    listF[j] = NULL;
	    nlistF[j] = 0;
	}
    }
    phgFree(nbsetss);
    phgFree(listEtype);
    return gL;
}
#else

FLOAT 
get_face_area(FLOAT *p0, FLOAT *p1, FLOAT *p2)
{
    FLOAT x1 = p0[0], y1 = p0[1];
    FLOAT x2 = p1[0], y2 = p1[1];
    FLOAT x3 = p2[0], y3 = p2[1];

    return .5 * fabs(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
}

void 
get_face_normals(FLOAT *p0, FLOAT *p1, FLOAT *p2, FLOAT *normals)
{
    FLOAT x1 = p0[0], y1 = p0[1];
    FLOAT x2 = p1[0], y2 = p1[1];
    FLOAT x3 = p2[0], y3 = p2[1];

    FLOAT det = .5 * (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
    char sgn = 1.;

    if (det > 0)
	sgn = 1.;
    else
	sgn = -1;

    FLOAT *n1 = normals + 0;	/* 1,2 */
    FLOAT *n2 = normals + 2;	/* 2,3 */
    FLOAT *n3 = normals + 4;	/* 3,1 */
    
    FLOAT x, y, d;

#define SET_NORM(X, Y, N)			\
    x = (X); y = (Y);				\
    d = sqrt(x*x + y*y);			\
    x /= d;  y /= d;				\
    N[0] = sgn*y; N[1] = -sgn*x;		

    SET_NORM(x2 - x1, y2 - y1, n1);
    SET_NORM(x3 - x2, y3 - y2, n2);
    SET_NORM(x1 - x3, y1 - y3, n3);
}


/* import 2D mesh, and layer
 * For simple triangle format.
 * */
LAYERED_MESH *
import_layered_mesh(char *file_tria,
		    char *file_layer,
		    char *file_nodeZ,
		    char *file_dual,
		    int nproc)
{
    LAYERED_MESH *gL;
    FILE *fp;
    char token[1024]; 
    INT i, j, n, k, Nly;
    INT numnp, nelem, ngrps, nbsets;
    FLOAT *R;
    int verb;

    PHG_CALLOC(gL, 1);
 
    /***********************/
    /* Control Information */
    /***********************/
    fp = NULL;
    if (0) {
    error:
	phgError(1, "invalid or unsupported input file, abort.\n");
    }


    /*********************/
    /* Nodal Coordinates */
    /*********************/
    fp = fopen(ns_params->vert_file, "r");
    assert(fp != NULL);
    GET_TOKEN; numnp = atoi(token);
    GET_TOKEN; GET_TOKEN; GET_TOKEN; 
    n = numnp;
        /* all 2D nodes */
    phgInfo(0, "number of vertices: %d\n", n);
    gL->verts = calloc(n, sizeof(*gL->verts));
    for (i = 0; i < n; i++) {
	READ_NUMBER;
	GET_TOKEN; gL->verts[i][0] = atof(token);
        /* x coord */
	GET_TOKEN; gL->verts[i][1] = atof(token);
        /* y coord */
	gL->verts[i][2] = -0.1;	/* base */
	gL->verts[i][3] = 1;	/* surf */
    //printf("%f, %f\n", gL->verts[i][0], gL->verts[i][1]);
	READ_NUMBER;
    }
    gL->nvert = n;
        /* # of all  2D nodes */
    UNUSED_LINE_CHECK("# end\n");
    fclose(fp);

    /*****************************/
    /* Element/Cell Connectivity */
    /*****************************/
    fp = fopen(file_tria, "r");
    assert(fp != NULL);
    GET_TOKEN; nelem = atoi(token);
    n = nelem;
        /* # of 2D triangles */
    PHG_CALLOC(gL->trias, n);
    phgInfo(0, "number of triangles: %d\n", n);
    GET_TOKEN; GET_TOKEN;
    for (i = 0; i < n; i++) {
	READ_NUMBER;
	/* c convention */
	GET_TOKEN; gL->trias[i].verts[0] = atoi(token) -1;
	GET_TOKEN; gL->trias[i].verts[1] = atoi(token) -1;
	GET_TOKEN; gL->trias[i].verts[2] = atoi(token) -1;
	gL->trias[i].area = 
	    get_face_area(gL->verts[gL->trias[i].verts[0]],
			  gL->verts[gL->trias[i].verts[1]],
			  gL->verts[gL->trias[i].verts[2]]);
        /* gL->verts[# of nodes][0;1] is x;y coord, gL->trias[i].verts[0] is the # of the first node. So we give the address of each node coord data to get_face_area and get the face area of the triangle */
	get_face_normals(gL->verts[gL->trias[i].verts[0]],
			 gL->verts[gL->trias[i].verts[1]],
			 gL->verts[gL->trias[i].verts[2]], 
			 gL->trias[i].normals);
        /* same reason of above */
	//READ_NUMBER;
    //printf("%d %d %d\n", gL->trias[i].verts[0], gL->trias[i].verts[1], gL->trias[i].verts[2]);
    }
    gL->ntria = n;
        /* # of 2D triangles */
    UNUSED_LINE_CHECK("# end\n");;
    fclose(fp);


#if 0
    /* 
     * Dual mesh
     *  edge length
     * */
    fp = fopen(file_dual, "r");
    phgInfo(0, "dual mesh edge\n");
    assert(fp != NULL);
    {
	int nelem0;
	GET_TOKEN; nelem = atoi(token);
	assert(nelem0 == nelem);
    }
    n = nelem;
    for (i = 0; i < n; i++) {
	GET_TOKEN; gL->trias[i].dual_edge[0] = atof(token);
	GET_TOKEN; gL->trias[i].dual_edge[1] = atof(token);
	GET_TOKEN; gL->trias[i].dual_edge[2] = atof(token);
	//READ_NUMBER;
    }
    UNUSED_LINE_CHECK("# end\n");;
    fclose(fp);
#else
    assert(file_dual == NULL);
#endif



    /**************/
    /* Layer Info */
    /**************/
    fp = fopen(file_layer, "r");
    assert(fp != NULL);
    fscanf(fp, "%d", &Nly);
    phgInfo(0, "Layer: %d\n", Nly);
    gL->max_nlayer = Nly;
    gL->layer_ratio = phgCalloc(Nly + 1,
				sizeof(*gL->layer_ratio));
    for (i = 0; i < Nly+1; i++) {
	fscanf(fp, "%lf", &gL->layer_ratio[i]);
    //printf("%lf\n", gL->layer_ratio[i]);
	phgInfo(0, "  layers[%3d]: %e\n", i, gL->layer_ratio[i]);
    }
        /* layer_ratio keeps the ratios of all layers, e.g. [0, 0.5, 1] for the case of 2 layers */
    fclose(fp);


    /******************/
    /* Vert List Info */
    /******************/
    fp = fopen(file_nodeZ, "r");
    PHG_CALLOC(gL->vert_global_lists, numnp);
        /* gL->vert_global_lists=calloc(numnp, sizeof **vert_global_lists) */
    for (i = 0; i < numnp; i++) {
	int nv, id;
	fscanf(fp, "%d", &id);
	assert(i == id-1);
	fscanf(fp, "%d", &nv);
	PHG_CALLOC(gL->vert_global_lists[i], nv + 1);
        /* gL->vert_global_lists[i]=calloc(nv+1, sizeof *vert_global_lists)*/	
	
	/* fisrt as num of verts in list */
	gL->vert_global_lists[i][0] = nv; 

	for (j = 0; j < nv; j++) {
	    fscanf(fp, "%d", &id);
	    gL->vert_global_lists[i][j+1] = id -1;
	}
    }
        /* the vert_global_lists is like this:
         * it is a 2D matrix, the first dim represents 2D node id and the second dim contans nv values, the first one is the layer number and the 2-nv values represent the id num of 3D nodes */ 
    fclose(fp);

    return gL;
}
#endif


void
build_layered_mesh(GRID *g, LAYERED_MESH *gL)
{
    SIMPLEX *e;
    TRIA *t;
    int i, j, k;
    int verb = 0;

    gL->g = g;

    /* ------------------------------------------------------------
     *
     *  Step 1:
     *    Build vertex relation: 
     *
     * ------------------------------------------------------------ */

    /* Step 1.1:
     *   3D mesh vertex connection */
    VEFMAP_3D *vef_3d;
    PHG_CALLOC(vef_3d, g->nvert); 
        /* vef_3d is a matrix of a size of total local verts, with all values 0 */
        /* PHG_CALLOC(p, n) p = phgCalloc(n, sizeof(*p));
         * So vef_3d is a matrix with a length of g->nvert */
    gL->vef_3d = vef_3d;
    for (i = 0; i < g->nvert; i++) {
	vef_3d[i].size = 20;
	PHG_CALLOC(vef_3d[i].v, vef_3d[i].size);
        /* vef_3d[i].v is a matrix with a length of 20 */
    }

    ForAllElements(g, e) {
	VEFMAP_3D *vf;
	int v[NVert];
	for (i = 0; i < NVert; i++)
	    v[i] = e->verts[i];
            /* NVert = 4, e->verts get the local ids of 4 verts. v is thus a matrix containing 4 values of 3D local ids */
	for (i = 0; i < NVert; i++) {
	    vf = vef_3d + v[i];
            /* get the # i VEFMAP_3D */
	    for (j = 0; j < NVert; j++) {
		for (k = 0; k < vf->nv; k++) 
		    if (v[j] == vf->v[k])
			break;
		if (k == vf->nv) {	/* new point connection
					 * v[i] ~ v[j]
					 * */
		    if (vf->nv >= vf->size) {
			int old_size = vf->size;
			vf->size += 20;
			PHG_REALLOC(vf->v,
				    vf->size, old_size)
		    }
		    vf->v[vf->nv++] = v[j];
		}
	    }
	}
    }
        /* vf(vef_3d) is a VEFMAP_3D kind of thing like this:
         * every local node has a value of vf
         * each vf has an attribute of nv, which is NVert, and an attribute of v[nv] which contains 4 values of the local index of the nodes in this element
         */

    /* debug */
    /* for (i = 0; i < g->nvert; i++) { */
    /* 	SHOW_iV(vef_3d[i].v, vef_3d[i].nv); */
    /* } */


    /*
     * Step 1.2:
     *   Build Vert L2S
     * */
    int *vert_list0, *found;
    /* map bottom verts (2D) to global (3D) */
    PHG_CALLOC(vert_list0, gL->nvert); /* only first vert,
					* the one on bottom */
    for (i = 0; i < gL->nvert; i++) {
	vert_list0[i] = gL->vert_global_lists[i][1];
        /* vert_global_lists[i][1] is the id of the 3D nodes at the bottom layer */ 
    }


    /* Map local verts (3D) to (2D) 
     * Note:
     *   assume gL->vert_global_list[.][1] is increasing !!! */
    gL->nvert_bot = 0;
    int nvb = 0;
    for (i = 0; i < g->nvert; i++) 
	if (g->types_vert[i] & BC_BOTTOM) 
	    nvb++;
            /* num of the verts at the bottom */
            /*NEED TO ADD BC_ISHELF !!!!!!!!!!!!!!!*/
    gL->nvert_bot = nvb;
    PHG_CALLOC(gL->vert_bot_Lidx, gL->nvert_bot);
        /* vert_bot_Lidx is a matrix with a length of nvb */
    PHG_CALLOC(gL->vert_bot_Gidx, gL->nvert_bot);
        /* vert_bot_Gidx is a matrix with a length of nvb */
    nvb = 0;
    for (i = 0; i < g->nvert; i++) 
	if (g->types_vert[i] & BC_BOTTOM) 
	    gL->vert_bot_Lidx[nvb++] = i;
            /* get the local index of the verts at the bottom */

    nvb = 0;
    PHG_CALLOC(gL->vert_L2S, g->nvert);
        /* vert_L2S is a matrix with a length of g->nvert */
    PHG_CALLOC(gL->vert_local_lists, g->nvert);
        /* vert_local_lists is a matrix with a length of g->nvert */
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] & BC_BOTTOM) {
	    //printf(" BC vert[%d]\n", i);
	    INT iS, iG = GlobalVertex(g, i);
            /* iG is the global index and i is the local index */
	    found = bsearch(&iG, vert_list0,
			    gL->nvert, sizeof(INT), phgCompINT);
            /* iG is what we are looking for
             * We go to vert_list0 to look for iG
             * there are gL->nvert values in vert_list0
             * every value has a length of sizeof(INI)
             * what is phgCompINI??
             * Anyway, if we find iG in vert_list0, we get the pointer found pointing to iG in vert_list0
             */
	    assert(found != NULL);
	    gL->vert_L2S[i] =
		iS = found - vert_list0;
            /* iS is the # of Shunxu counting at the bottom
             * vert_L2S is a matrix like this:
             * it has a len of g->nvert
             * if comes to bottom, its value is set to the Shunxu counting of bottom nodes, otherwise is zero
             */ 

	    INT nly, *vG, *vL;
	    gL->vert_bot_Gidx[nvb] = iS;
            /* vert_bot_Gidx keeps the # of Shunxu counting at the bottom */
	    vG = gL->vert_global_lists[iS];
            /* vG means verts of global
             * vG is a vector with a length of nvb
             * each value of vG contains nly+1 values of ids of 3D nodes
             * e.g. vG[0]={3,1,2,3}
             */
	    nly = vG[0];
            /* so nly is the # of layers */
	    //SHOW_iV(vG, nly+1);
	    PHG_CALLOC(vL, nly+1);
            /* vL is a matrix with a len of nly+1 */
	    gL->vert_local_lists[i] = vL;
            /* vert_local_lists is a 2D matrix
             * it looks like this:
             * if comes to bottom, vert_local_lists contains 4 values of data, now they are all 0, otherwise the values are set 0
             */
	    vL[0] = nly;		/* first as num */
	    vL[1] = i;			/* second as self */
	    for (j = 2; j < nly+1; j++) {
		/* In the neighbor of prev vert will
		 * be the current vert. */
		VEFMAP_3D *vf = vef_3d + vL[j-1];
		for (k = 0; k < vf->nv; k++) {
		    if (vG[j] == GlobalVertex(g, vf->v[k])) 
			break;
		}
		assert(k < vf->nv);
		vL[j] = vf->v[k];
	    }
	    //SHOW_iV(vL, nly+1);
	    //printf("    add: %x\n", vL);
	    nvb++;
	} else {
	    //printf(" Ly vert[%d]\n", i);
	    gL->vert_L2S[i] = -1;
	    gL->vert_local_lists[i] = NULL;
	}
	/* phgInfo(0, "vert local list %d: %x\n",  */
	/* 	i, gL->vert_local_lists[i] */
	/* 	); */
    }
    

    /*
     * Create gather communicator
     *
     *  */
    if (1) {
	static int pass = 0;
	int rank, *vert_bot_cnts, *vert_bot_dsps, *vert_bot_idxs, ntotal;

	phgInfo(0, "[%d] Build vert 2D gather comm: pass %d\n", 
		phgRank, pass++);
        PHG_CALLOC(vert_bot_cnts, g->nprocs);
        PHG_CALLOC(vert_bot_dsps, g->nprocs);
	MPI_Allgather(&gL->nvert_bot, 1, MPI_INT, 
		      vert_bot_cnts, 1, MPI_INT, g->comm);
	ntotal = 0;
	for (rank = 0; rank < g->nprocs; rank++) {
	    vert_bot_dsps[rank] = ntotal;
	    ntotal += vert_bot_cnts[rank];
	}
	PHG_CALLOC(vert_bot_idxs, ntotal);
	MPI_Allgatherv(gL->vert_bot_Gidx, gL->nvert_bot, PHG_MPI_INT,
		       vert_bot_idxs, vert_bot_cnts, vert_bot_dsps, 
		       MPI_INT, g->comm);

	
	gL->vert_bot_total = ntotal;
	gL->vert_bot_cnts = vert_bot_cnts;
	gL->vert_bot_dsps = vert_bot_dsps;
	gL->vert_bot_idxs = vert_bot_idxs;

	SHOW_iV_(0, gL->vert_bot_cnts, g->nprocs);
	SHOW_iV_(0, gL->vert_bot_dsps, g->nprocs);
	SHOW_iV_(0, gL->vert_bot_idxs,
		 gL->vert_bot_total);
    }


    /* ------------------------------------------------------------
     *
     *  Step 1:
     *    Build triangle relation: VEFmap_2d, mapE3to2
     *
     * ------------------------------------------------------------ */


    /* Build VEFmap_2d */
    VEFMAP_2D *vef_2d;
    PHG_CALLOC(vef_2d, gL->nvert);
    gL->vef_2d = vef_2d;
    t = gL->trias;
    for (j = 0; j < gL->ntria; j++, t++) {
	for (k = 0; k < 3; k++) {
	    i = t->verts[k];
	    assert(i < gL->nvert);
	    int ne = vef_2d[i].ne;
	    assert(ne < 12);
	    vef_2d[i].idx[ne] = j;
	    vef_2d[i].ne++;
	    //SHOW_iV(vef_2d[i].idx, vef_2d[i].ne);
	}
    }

    /* Build mapE3to2 */
    int nface_bot = 0;
    BOTTOM_FACE *face_bot, *fb;

    ForAllElements(g, e) 
	for (k = 0; k < NFace; k++) 
	    if (e->bound_type[k] & BC_BOTTOM)
		nface_bot++;
    PHG_CALLOC(face_bot, nface_bot);
    phgInfo(0, "Bottom face: %d\n", nface_bot);
    gL->nface_bot = nface_bot;
    gL->face_bot = face_bot;

    VTK_VERB(3);
    fb = face_bot;
    ForAllElements(g, e) {
	for (k = 0; k < NFace; k++) {
	    int v[3], iS[3], iL[3], ii;
	    if (!(e->bound_type[k] & BC_BOTTOM))
		continue;

	    /* Bottom face */
	    GetFaceVertices(e, k, v[0], v[1], v[2]);
	    for (i = 0; i < 3; i++) {
		iL[i] = e->verts[v[i]];
		iS[i] = gL->vert_L2S[iL[i]];
		assert(iS[i] >= 0 && iS[i] < gL->nvert);
	    }

	    /* find tria contains v[012] */
	    int n0, n1, n2, *t0, *t1, *t2,
		i0, i1, i2;
	    //printf("v: [%d %d %d]\n", v[0], v[1], v[2]);
	    n0 = vef_2d[iS[0]].ne;
	    n1 = vef_2d[iS[1]].ne;
	    n2 = vef_2d[iS[2]].ne;
	    t0 = vef_2d[iS[0]].idx;
	    t1 = vef_2d[iS[1]].idx;
	    t2 = vef_2d[iS[2]].idx;
	    //SHOW_iV(t0, n0);
	    //SHOW_iV(t1, n1);
	    //SHOW_iV(t2, n2);
#if 1
	    /* faster */
	    t0 = vef_2d[iS[0]].idx;
	    for (i0 = 0; i0 < n0; i0++, t0++) {
		t1 = vef_2d[iS[1]].idx;
		for (i1 = 0; i1 < n1; i1++, t1++) {
		    if (*t0 == *t1) {
			t2 = vef_2d[iS[2]].idx;
			for (i2 = 0; i2 < n2; i2++, t2++) {
			    if (*t0 == *t2) {
				//printf("   found: %d\n", *t0);
				break;
			    }
			}
			if (i2 < n2)	/* found */
			    break;	
		    }
		}		
		if (i1 < n1)		/* found */
		    break;
	    }
#else
	    t0 = vef_2d[iS[0]].idx;
	    for (i0 = 0; i0 < n0; i0++, t0++) {
		t1 = vef_2d[iS[1]].idx;
		for (i1 = 0; i1 < n1; i1++, t1++) {
		    t2 = vef_2d[iS[2]].idx;
		    for (i2 = 0; i2 < n2; i2++, t2++) {
			if (*t0 == *t2 && *t0 == *t1) {
			    printf("   found: %d\n", *t0);
			    break;
			}
		    }
		    if (i2 < n2)	/* found */
			break;
		}		
		if (i1 < n1)		/* found */
		    break;
	    }
#endif
        //printf("i0: %d  n0: %d\n", i0, n0);
	    assert (i0 < n0);		/* found */
	    fb->e = e;
	    fb->face = k;
	    for (j = 0; j < 3; j++)
		fb->vert[j] = iL[j];
	    fb->tria = *t0;

	    /* debug */
	    {
		FLOAT c0[3], c1[3];
		Bzero(c0);
		Bzero(c1);
		for (i = 0; i < 3; i++)
		    for (j = 0; j < 3; j++) {
			c0[j] += g->verts[iL[i]][j] / 3.;
			c1[j] += gL->verts[iS[i]][j] / 3.; 
		    }
		vtkSetColor(verb, "yellow");
		vtkDrawLine(verb, c0, c1);
	    }

	    fb++;
	}
    }



    /* ------------------------------------------------------------
     *
     *  Step 2:
     *   Build tetra relation:    VEFmap_3d, mapE2to3
     *
     * ------------------------------------------------------------ */
    phgInfo(0, "Tria to tetra relation\n");
    VTK_VERB(3);
    //assert(nface_bot == gL->ntria); /* Fixme */
    fb = face_bot;
    e = NULL;
    for (i = 0; i < nface_bot; i++, fb++) {
	SIMPLEX *e1 = fb->e, *e0 = NULL;
	int f1 = fb->face, f0 = -1;
	int *vL[3], n[3], iv[3], ne = 0;
	
	for (k = 0; k < 3; k++) {
	    assert(g->types_vert[fb->vert[k]] & BC_BOTTOM);
	    //printf(" BC vert[%d]\n", fb->vert[k]);
	    vL[k] = gL->vert_local_lists[fb->vert[k]];
	    //printf("    add: %x\n", vL[k]);
	    assert(vL[k] != NULL);
	    n[k] = vL[k][0];
	    iv[k] = 1;		/* start */
	    //SHOW_iV(vL[k], n[k]+1);
	    if (n[k] > ne)
		ne = n[k];
	}
	
	/* tets in layers */
	ne *= 3;
	PHG_CALLOC(fb->elems, ne);
	PHG_CALLOC(fb->faces, 3*ne);
	fb->ne = 0;

	/* Algorithm desp:
	 * Given tet(e0) and face(f1), find the next face(f1)
	 *   and next tet(e1) in layer.
	 * 1. find the vert not on the face, it must belongs
	 *    to one of three vertical line which contains,
	 *    say, v(j).
	 * 2. find the face opsite to v(j), the next tet
	 *    is the neigh of current tet on face v(j)
	 *
	 * */

	while (TRUE) {
	    int v[4], ff;

	    e0 = e1;
	    f0 = f1;
	    fb->elems[fb->ne++] = e0;

	    vtkSetColor(verb, "yellow");
	    vtkDrawElement(verb, e0);

	    for (k = 0; k < 3; k++) 
		v[k] = GetFaceVertex(f0, k);
	    v[3] = f0;
	    for (k = 0; k < 3; k++)
		if (iv[k]+1 < n[k]+1
		    && vL[k][iv[k]+1] == e0->verts[f0])
		    break;

	    assert(k < 3);	/* advance on line[k] */
	    for (j = 0; j < 3; j++)
		if (e0->verts[v[j]] == vL[k][iv[k]])
		    break;
	    assert(j < 3);	/* elem v[j] on line vL[k],
				 * and the next face */
	    iv[k]++;
	    ff = v[j];	

	    vtkSetColor(verb, "red");
	    vtkDrawFace(verb, e0, ff);
	    //vtkPause(0);

	    /* reach top */
	    if (e0->neighbours[ff] == NULL) {
		break;
	    } else {
		assert(e0->bound_type[ff] & INTERIOR);
		e1 = e0->neighbours[ff];
		for (k = 0; k < NFace; k++)
		    if (e1->faces[k] == e0->faces[ff])
			break;
		assert(k < NFace);
		f1 = k;
	    }
	} /* end layers */
	vtkPause(verb);


	/* Flux face */
	if (1) {
	    /* Normals */
	    TRIA *t = gL->trias + fb->tria;
	    FLOAT *normals = t->normals;
	    int ii;

	    for (ii = 0; ii < fb->ne; ii++) { /* elems */
		SIMPLEX *e = fb->elems[ii];
		
		for (j = 0; j < 3; j++) { /* tri faces */
		    FLOAT *n = normals + j * 2;
		    for (k = 0; k < NFace; k++) {
			const FLOAT *nn
			    = phgGeomGetFaceOutNormal(g, e, k);
			FLOAT nt[3] = {n[0], n[1], 0.};
			FLOAT r = INNER_PRODUCT(nn, nt);
			if (1. - r < 1e-7)
			    break;
		    }
		    if (k < NFace) {
			fb->faces[ii*3 + j] = k;

			/* check */
			int kk = 0;
			for (kk = 0; kk < j-1; kk++)
			    assert(k != fb->faces[ii*3 + kk]);
			    
		    } else {
			fb->faces[ii*3 + j] = -99;
		    }
		} /* end tria face */
	    }	  /* end elements */
	}	  /* end flux */
    }

    phgFree(vert_list0);



#if 1
    /* ------------------------------------------------------------
     * 
     *    Dual mesh
     *
     *  Build fb->edge2to3[ne][3][2]
     *    triagle edges [3] to tetra edges [ne][2]
     *    if tri_edge -> only one tet edge
     *    then edge2to3[ii][k][1] = -1 
     *
     * ------------------------------------------------------------
     * */
    fb = gL->face_bot;
    //e = NULL;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	TRIA *t = gL->trias + fb->tria;
	FLOAT area_tri = t->area;
	int ii;

	int e2v_tri[3][2] = {
	    {0, 1}, {1, 2}, {2, 0}
	};
	INT vert_tri[3] = {t->verts[0],
			   t->verts[1],
			   t->verts[2]};
    //printf("vert_tri: %d, %d, %d\n", vert_tri[0], vert_tri[1], vert_tri[2]);
    // local # of three verts
	FLOAT *x_tri[3] = {gL->verts[vert_tri[0]], 
			   gL->verts[vert_tri[1]], 
			   gL->verts[vert_tri[2]]};
	FLOAT mid_tri[3][2];
	FLOAT mid_tet[NEdge][3];

	PHG_CALLOC(fb->edge2to3, 6 * fb->ne);

	for (k = 0; k < 3; k++) { /* edge */
	    mid_tri[k][0] = .5*(  x_tri[ e2v_tri[k][0] ][0]
				+ x_tri[ e2v_tri[k][1] ][0])*1000.0;
	    mid_tri[k][1] = .5*(  x_tri[ e2v_tri[k][0] ][1]
				+ x_tri[ e2v_tri[k][1] ][1])*1000.0;
	}

	for (ii = 0; ii < fb->ne; ii++) { /* tets */
	    SIMPLEX *e = fb->elems[ii];
	    int edge3to2[NEdge];
	    int edge2to3[3][2];
	    
	    /* tetra edge */
	    for (j = 0; j < NEdge; j++) {
		int v0 = GetEdgeVertex(j, 0);
		int v1 = GetEdgeVertex(j, 1);
		FLOAT *x0 = g->verts[e->verts[v0]];
		FLOAT *x1 = g->verts[e->verts[v1]];
		
		mid_tet[j][0] = .5 * (x0[0] + x1[0]);
		mid_tet[j][1] = .5 * (x0[1] + x1[1]);
		mid_tet[j][2] = .5 * (x0[2] + x1[2]);

		for (k = 0; k < 3; k++) {
		    FLOAT d[2] = {mid_tet[j][0] - mid_tri[k][0], 
				  mid_tet[j][1] - mid_tri[k][1]};
            //printf("mid_tet(j,0)=%f, mid_tri(k,0)=%f\n", mid_tet[j][0], mid_tri[k][0]);
            //printf("mid_tet(j,1)=%f, mid_tri(k,1)=%f\n", mid_tet[j][1], mid_tri[k][1]);
		    if ((d[0]*d[0] + d[1]*d[1]) < 1e-9)
			break;
		}
		if (k < 3) {
		    edge3to2[j] = k; /* found */
		} else {
		    edge3to2[j] = -1;
		}
            //printf("j=%d, k=%d, edge3to2[j]=%d\n",j,k,edge3to2[j]);
	    }
	    
	    {
	    	/*
	    	 * check: edge3to2
	    	 *     tet edge map to tri edge
	    	 *     e.g. 2 k1, 2 k2, only 1 k3
	    	 * */
	    	int nk[3] = {0, 0, 0}, k0;
	    	for (j = 0; j < NEdge; j++)
                nk[edge3to2[j]]++;
	    	for (k = 0; k < 3; k++)
	    	    if (nk[k] == 1)
	    		break;
	    	assert(k < 3);
	    	k0 = k;
	    	for (k = 0; k < 3; k++)
	    	    if (k != k0)
	    		assert(nk[k] == 2);
	    }
	 
	    for (k = 0; k < 3; k++) {
		int v0 = -1, v1 = -1;
		for (j = 0; j < NEdge; j++)
        {
		    if (edge3to2[j] == k) {
			if (v0 == -1)
			    v0 = j;
			else 
			    v1 = j;
            } }
		assert(v0 != -1);
		if (v1 == -1) {
		    edge2to3[k][0] = v0;
		    edge2to3[k][1] = -1; /* note 2nd one */
		} else {
		    edge2to3[k][0] = v0;
		    edge2to3[k][1] = v1;
		}
	    }

	    for (k = 0; k < 3; k++) {
		fb->edge2to3[ii * 6 + k*2  ] = edge2to3[k][0];
		fb->edge2to3[ii * 6 + k*2+1] = edge2to3[k][1];
	    }
	}
    }		  /* end bot face */


    /* Volume of prism colume */
    phgInfo(0, "* Initial voulme\n");
    gL->volumes = phgCalloc(gL->ntria, sizeof(FLOAT));
    t = gL->trias;
    for (j = 0; j < gL->ntria; j++, t++) {
	int I[3] = {t->verts[0],
		    t->verts[1],
		    t->verts[2]};
	FLOAT *coord[3];

	for (i = 0; i < 3; i++) {
	    coord[i] = gL->verts[I[i]];
	}
	FLOAT x, y, z, X[3];
	x = (coord[0][0]+ coord[1][0] + coord[2][0]) / 3.;
	y = (coord[0][1]+ coord[1][1] + coord[2][1]) / 3.;
	z = 1.;			/* Note: height assume to be 1 */

	func_ice_slab(x, y, z, X);
	gL->volumes[j] = t->area * X[2];
	phgInfo(3, " X[%e %e %e], height %e\n", x, y, z, X[2]);
	phgInfo(3, "vol[%3d] %e\n", j, gL->volumes[j]);
    }    
#else
#  warning Daul mesh disabled!
#endif


#if 1
    /* Height of points */
    gL->height =
	phgCalloc(gL->nvert, sizeof(*gL->height));
    gL->height_local =
	phgCalloc(gL->nvert_bot, sizeof(*gL->height_local));
    gL->height_gather =
	phgCalloc(gL->vert_bot_total, sizeof(*gL->height_gather));
    gL->bottom=
	phgCalloc(gL->nvert, sizeof(*gL->bottom));
    gL->bottom_local =
	phgCalloc(gL->nvert_bot, sizeof(*gL->bottom_local));
    gL->bottom_gather =
	phgCalloc(gL->vert_bot_total, sizeof(*gL->bottom_gather));

    build_layered_mesh_height(g, gL);
    //check_height(g, gL);
#else
    /* Use ice-grid. */
#endif

    return;
}


void 
build_layered_mesh_height(GRID *g, LAYERED_MESH *gL)
{
    SIMPLEX *e;
    TRIA *t;
    int i, j, k;
    int verb = 0;
    
    INT nvb = 0;
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] & BC_BOTTOM) {
	    assert(gL->vert_local_lists[i] != NULL);

	    int nv = gL->vert_local_lists[i][0];
	    int *iL = &gL->vert_local_lists[i][1];
	    gL->height_local[nvb] = 
		g->verts[iL[nv-1]][Z_DIR]; /* coord z */
	    gL->bottom_local[nvb] = 
		g->verts[iL[0]][Z_DIR]; /* coord z */
	    
	    phgInfo(3, "  Asgn Bidx: %5d, Gidx: %5d height: %e\n",
		    gL->vert_bot_Gidx[nvb],
		    GlobalVertex(g, iL[nv-1]),
		    g->verts[iL[nv-1]][Z_DIR]);
	    
	    nvb++;
	}
    }
    MPI_Allgatherv(gL->height_local, gL->nvert_bot, PHG_MPI_FLOAT,
		   gL->height_gather, gL->vert_bot_cnts, gL->vert_bot_dsps, 
		   PHG_MPI_FLOAT, g->comm);
    MPI_Allgatherv(gL->bottom_local, gL->nvert_bot, PHG_MPI_FLOAT,
		   gL->bottom_gather, gL->vert_bot_cnts, gL->vert_bot_dsps, 
		   PHG_MPI_FLOAT, g->comm);
    for (i = 0; i < gL->vert_bot_total; i++) {
	INT idx = gL->vert_bot_idxs[i];
	gL->height[idx] = gL->height_gather[i];
	gL->bottom[idx] = gL->bottom_gather[i];
    }    

    /* debug */
    if (0 && g->rank == 0) {
	FILE *fp = fopen("./output/H0_.m", "w");
	fprintf(fp, "H0 = [\n");
	for (i = 0; i < gL->nvert; i++) {
	    fprintf(fp, "%e\n", gL->height[i]);
	}
	fprintf(fp, "];\n");
	fclose(fp);
    }


    phgPrintf("   Build layers using height ratio\n");
    const FLOAT *ratio = gL->layer_ratio;
    int ii;
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	int iv = gL->vert_L2S[i];
	assert(nv > 0);

	FLOAT h0, h1;
	h0 = gL->bottom[iv];
	h1 = gL->height[iv];

	FLOAT H[nv];
	get_layer_height(H, nv, ratio, h0, h1);
		
	assert(gL->max_nlayer + 1 == nv);
	for (j = 0; j < nv; j++) {
	    g->verts[iL[j]][Z_DIR] = H[j];
	}
    }

    phgGeomInit_(g, TRUE);
    phgExportVTK(g, OUTPUT_DIR "maped.vtk", NULL);
}

//#include <metis.h>
#include <parmetis.h>

#if PARMETIS_MAJOR_VERSION == 4
# define ID_TYPE idx_t
#else
# define ID_TYPE idxtype
#endif

static ID_TYPE *idx_array;

void
part_layered_mesh(GRID *g, LAYERED_MESH *gL)
{
    /* Fix me: part on the fly */
    assert(g->nleaf == g->nleaf_global);
    if (phgRank > 0)
	return;

    int ne = gL->ntria, nn = gL->nvert, etype = 1; /* triangle */
    int numflag = 0, edgecut = 0;		   /* c-style */
    ID_TYPE *epart = calloc(ne, sizeof(ID_TYPE));
    ID_TYPE *npart = calloc(nn, sizeof(ID_TYPE));
    int i, j;
    static int *part = NULL; 

#if PARMETIS_MAJOR_VERSION == 4
    ID_TYPE nprocs = phgNProcs;
    ID_TYPE *eptr = calloc(ne+1, sizeof(ID_TYPE));
    for (i = 0; i < ne+1; i++) {
	eptr[i] = 3 * i;
    }
#else
    int nprocs = phgNProcs;
#endif

    if (nprocs <= 1) {
        part = calloc(ne, sizeof(*part));
	return;
    }

    ID_TYPE *elmnts = calloc(3*ne, sizeof(ID_TYPE));
    idx_array = elmnts;
    /* bottom faces */
    for (i = 0; i < ne; i++) 
	for (j = 0; j < 3; j++) 
	    *(idx_array++) = gL->trias[i].verts[j];


#if 0
    /* Metis adaptive */
    phgPrintf("Use metis partitaion.\n");
    ParMETIS_V3_AdaptiveRepart(vtxdist, xadj, adjncy, vwgt, vsize, adjwgt,
			       &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, &ipc2redist, 
			       options, &edgecut, part, &comm);
#elif 1
    phgPrintf("Use metis dual partitaion.\n");
    /* Metis dual */
#  if PARMETIS_MAJOR_VERSION == 3
    METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &nprocs,
    		       &edgecut, epart, npart);
#  elif PARMETIS_MAJOR_VERSION == 4
    phgPrintf("PARMETIS version 4 !\n");
    idx_t ncommon = 2, objval;
    METIS_PartMeshDual(&ne, &nn, eptr, elmnts,
		       NULL, NULL, &ncommon, &nprocs,
		       NULL, NULL, &objval, epart, npart);
#  endif
#else
    /* Read in part info */
    phgPrintf("Use user partitaion.\n");
    FILE *fp = fopen("tri-part.dat", "r");
    int ne0;
    fread(&ne0, sizeof(int), 1, fp);
    assert(ne0 == ne);
    fread(epart, sizeof(*epart), ne, fp);
    fclose(fp);
#endif

    /* statistics */
    int *nep, max_ne, min_ne;
    part = calloc(ne, sizeof(*part));
    nep = calloc(nprocs, sizeof(*nep));
    //printf("\n --part: \n");
    for (i = 0; i < ne; i++) {
	part[i] = epart[i];
	//printf("   elem[%5d]: %2d\n", i, part[i]);
	nep[part[i]]++;
    }
    
    max_ne = -1e6;
    min_ne = 1e6;
    for (i = 0; i < nprocs; i++) {
	if (max_ne < nep[i])
	    max_ne = nep[i];
	if (min_ne > nep[i])
	    min_ne = nep[i];
    }
    phgPrintf("Ne per rank: [%d %d]\n", min_ne, max_ne);

#if 0
    /* output */
    FILE *fp = fopen("tri-part.dat", "w");
    fwrite(&ne, sizeof(ne), 1, fp);
    fwrite(part, sizeof(*part), ne, fp);
    fclose(fp);
#endif


    /* Bring to 3D */
    assert(gL->nface_bot == gL->ntria);
    BOTTOM_FACE *fb;
    fb = gL->face_bot;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	for (j = 0; j < fb->ne; j++) {
	    SIMPLEX *e = fb->elems[j];
	    e->mark = part[fb->tria];
	}
    }

    
    free(part);
    free(npart);
    free(epart);
    free(elmnts);
    return;
}


void 
destory_layerd_mesh(LAYERED_MESH **gL_ptr)
{
    LAYERED_MESH *gL = *gL_ptr;
    GRID *g = gL->g;
    int i, **p;

    phgInfo(0, "GRID nvert: %d\n", g->nvert);

    phgFree(gL->verts);
    phgFree(gL->trias);

    for (i = 0; i < gL->nvert; i++)
	phgFree(gL->vert_global_lists[i]);
    phgFree(gL->vert_global_lists);

    if ((p = gL->vert_local_lists) != NULL) {
	for (i = 0; i < g->nvert; i++, p++)
	    if (g->types_vert[i] & BC_BOTTOM)
		phgFree(*p);
	phgFree(gL->vert_local_lists);
    }

    phgFree(gL->vert_bot_cnts);
    phgFree(gL->vert_bot_dsps);
    phgFree(gL->vert_bot_Lidx);
    phgFree(gL->vert_bot_Gidx);
    phgFree(gL->vert_bot_idxs);
    phgFree(gL->vert_L2S);

    for (i = 0; i < gL->nface_bot; i++)
	phgFree(gL->face_bot[i].elems);
    phgFree(gL->face_bot);

    for (i = 0; i < g->nvert; i++)
	phgFree(gL->vef_3d[i].v);
    phgFree(gL->vef_3d);
    phgFree(gL->vef_2d);
    phgFree(gL);

    *gL_ptr = NULL;
    return;
}
