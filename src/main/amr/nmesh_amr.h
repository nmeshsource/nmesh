/* nmesh_amr.h */
/* (c) Wolfgang Tichy 2/2019 */

/* declarations from other parts that we need here already */
#include "../main/skeleton.h"
#include "../nMPI/nMPI_defs.h"


/* Main parts of a mesh:

tMesh:    ---------------------mesh-----------------------------------------
tPat:     |-----patch0-----|-------patch1-------|--patch2----|...
tNode:     node node ...     node node node ...  node node ...
                   ^
  The nodes shown here are the leaf nodes. They are kept in 
  linked lists (tNlist *lnodes).

  Each node has a tDat struct that can be empty if the data is on another
  proc. The tDat struct contain lists of arrays, one for each variable.

  Also, each node is part of an oct-tree. Here we only show 2 instead of 8:
level
l=0                  ___________rnode_____________
l=1          ______node______             ______node______
l=2     __node__            node     __node__            node
l=3  node      _node_                       node
l=4        node 
  We have one tree per patch.
  The ends of the tree are called leaf nodes.
*/


/* the data within a node, this should be only on one proc */
typedef struct tDAT {
  struct tNODE *node;     /* pointer to node dat is in */
  int nv;                 /* number of vars */
  int nvenabled;          /* number of enabled vars */
  struct tARRAY **v;      /* list of data pointers to vars, if v[vi]=NULL,
                             the var vi and its surfaces are not enabled */
  struct tSURFACE **s[6]; /* list of surfaces needed for data exchange,
                             e.g. s[0]=surfs in -X dir, s[3]=surfs in +Y dir,
                             if s[6][vi]=NULL var vi does not need exchange */
  struct tCOM *com[6];    /* com for each face */
} tDat;

/* surface data needed for node to neighbor node communication */
typedef struct tSURFACE {
  struct tDAT *dat;       /* pointer to dat the surface is in */
  int face;               /* face surface is on */
  int vi;                 /* var index */
  struct tARRAY *mysurf;  /* array that contains values at my surface points */
  int allocd_mysurf;      /* 1 if we need to free mysurf */
  int nnbsurf;            /* number of nbsurf */
  struct tARRAY **nbsurf; /* list of values from neighb. surfaces */
  int *allocd_nbsurf;     /* allocd_nbsurf[i]=1 if we need to free nbsurf[i] */
} tSurface;
/* NOTE:
   mysurf comes from this proc,
   nbsurf[i] can just point if nb[i] is local, otherwise we need to alloc */


/* a node */
typedef struct tNODE {
  struct tPAT *pat;       /* pointer to patch that contains node */
  struct tNODE *parent;   /* pointer to parent node */
  struct tNODE *child[8]; /* list of pointers to childeren nodes */
  struct tNODE *nb[6];    /* neighbs in +/-X,Y,Z dir: nb[+-dir], e.g.:
                             nb[4]= neigh in -Z dir, nb[1]= neigh in +X dir */
  int nfnb[6];            /* number of face neighbor nodes */
  struct tNODE **fnb[6];  /* list of neighbor nodes on face,
                             kept up to date by update_node_fnb */
  double bbox[6];         /* bounding box (in X,Y,Z) of this node */
  int patface[6];         /* whether node is at patch face 0,1,2,3,4,5 */
  int n[3];               /* number of points in X,Y,Z-directions */
  int np;                 /* np = n[0] * n[1] * n[2]; */
  int l;                  /* refinement level of this node */
  int leaf;               /* is 1 if this is a leaf node */
  int ijk;                /* node index (0-7), i.e. child number wrt. parent */
  long nid;               /* node ID, updated by update_mesh_myln_node_nid */
  //int lid;                /* local node ID */
  struct tARRAY *Dt[3];   /* transp. differentiation matrix in 3 dirs for [-1,1]
                             domain. This just points to an array in patch. */
  struct tARRAY *At[3];   /* transp. analysis matrix in 3 dirs */
  struct tARRAY *St[3];   /* transp. synthesis matrix in 3 dirs */
  struct tARRAY *Xb[3];   /* points we use e.g. Gauss-Lobatto: Xb\in[-1,1].
                             X=0.5*((a+b)+(b-a)*Xb), a=bbox[0], b=bbox[1] */
  struct tARRAY *Winteg[3];   /* integr. weights in 3 dirs */
  tDat *dat;              /* pointer to data (NULL if not on this proc) */
  int datrank;            /* rank of proc that rightfully has data */
  nMPI_Comm comm;         /* MPI_comm for this node, could contain only ranks
                             where dat is and where all neighb. have dat */
} tNode;

/* a linked list of nodes */
typedef struct tNLIST {
  tNode *node;
  struct tNLIST *next;
  struct tNLIST *prev;
} tNlist;


/* the nodes fill a patch */
typedef struct tPAT {
  double bbox[6];       /* bounding box (in X,Y,Z) of this patch */
  int p;                /* index of this patch */
  struct tMESH *mesh;   /* pointer to mesh that contains patch */
  /* funcs to compute X,Y,Z from x,y,z and vice versa: */
  int (*XYZ_Of_xyz)(struct tPAT *pat, double x, double y, double z, double *X, double *Y, double *Z);  /* func to compute X,Y,Z from x,y,z */
  int (*xyz_Of_XYZ)(struct tPAT *pat, double X, double Y, double Z, double *x, double *y, double *z);  /* func to compute x,y,z from X,Y,Z */
  tNode *rnode;         /* root node in this patch */
  int nmax;             /* max n[0],n[1],n[2] a node in this patch can have */
  struct tARRAY *(*Dt)[3]; /* list of transposed differentiation matrices
                              we store Dt[1...nmax][dir], where dir=0,1,2 */
  struct tARRAY *(*At)[3]; /* list of transposed analysis matrices */
  struct tARRAY *(*St)[3]; /* list of transposed synthesis matrices */
  struct tARRAY *(*Xb)[3]; /* list of points */
  struct tARRAY *(*Winteg)[3]; /* list of integration weights */
  //tNlist *lns;   /* start of linked list of leaf nodes in this patch */
} tPat;
/* Note: each patch should have Bfaces as in sgrid. But instead of pointlists
   we can use bounding rectangles in both adjacent bfaces, because we will
   only allow patches that are touching cubed spheres. So when a node needs
   data from the other side of a patch boundary, it can:
   1. figure out its bounding rectangle on the other side
   2. ask all nodes on the other side within the rectangle on the other side
      for data
   For this we need a node list of all leaf nodes on all faces. */


/* several patches and thus a list of leaf nodes make up the 
   computational mesh */
typedef struct tMESH {
  double dt;        /* time step */
  double time;      /* current time */
  int iteration;    /* current iteration number */
  tTodo *skel[NFUNCBINS]; // list of tTodo's from skeleton.c
  int nvdb;         /* number of variables */
  struct tVAR *vdb; /* variable data base */
  int vdb_iStart;   /* index we start at when searching for a var */
  int npdb;         /* number of mesh parameters */
  struct tPAR *pdb; /* parameter data base */
  int pdb_iStart;   /* index we start at when searching for a par */
  int npats;        /* number of patches */
  tPat **pat;       /* list of pointers to patches */
  tNlist *lns;      /* start of linked list of all leaf nodes */
  long nln;         /* total number of leaf nodes */
  int nmyln;        /* number of leaves in myln array */
  tNlist **myln;    /* elements of lns owned by this proc */
} tMesh;
/* NOTE: the list lns needs to be distributed among MPI jobs:
use space filling curve as in
http://www.speedup.ch/workshops/w42_2013/carsten.pdf
*/


/***********************************************************************/
/* other useful objects */
/***********************************************************************/
/* arrays */
typedef struct tARRAY {
  int n[3];     /* dims in all 3 dirs */
  int N;        /* N = n[0] * n[1] * n[2]; */
  double *d;    /* pointer to double data (could add one more for GPU data) */
  int d_nofree; /* a_nofree=1 if free_array should not free a */
  int ns;       /* number of segments */
  int si;       /* segment index */
//  void *par;    /* pointer to some extra pars */
} tArray;


/* variable lists */
typedef struct tVARLIST {
  double time;
  int n;
  int *index;
  tMesh *mesh;  /* pointer to mesh to which vars belong */
  void *vlPars; /* A pointer that is usually NULL, but can point to some
                   object that contains special extra pars or info. This
                   pointer is not touched by the funcs in variables.c (such
                   as vlduplicate, vlcopy, vlfree, ...). So the user has 
                   to manage it: e.g. free it, before calling vlfree. */
} tVarList;


/**************************************************************************/
/* loops */
/**************************************************************************/
#include "nmesh_amr_loops.h"

/**************************************************************************/
/* functions */
/**************************************************************************/

/* mesh.c */
tMesh *make_empty_mesh(int pr);


/* storage.c */
tArray *alloc_array_with_segs(int n[3], int ns);
tArray *alloc_array(int n[3]);
tArray *get_array_seg(tArray *array, int si);
void point_array_a_to_data(tArray *array, void *data);
void free_array(tArray *array);
tMesh *alloc_mesh(int npats);
void realloc_patlist_in_mesh(tMesh *mesh, int npats);
void free_mesh(tMesh *mesh);
tPat *alloc_patch(tMesh *mesh, int p, int nD);
void free_patch(tPat *pat);
tNode *alloc_node(void);
void free_node(tNode *node);
tNode *make_root_node(tPat *pat, int n[3], int datrank);
tNlist *make8_child_nodes(tNode *parent, int n[3]);
tDat *alloc_dat(tNode *node);
void free_dat(tDat *dat);
tNlist *alloc_nodelist(tNode *node);
tNlist *addnode_to_nodelist_after(tNlist *elem, tNode *node);
tNlist *addnode_to_nodelist_before(tNlist *elem, tNode *node);
tNlist *copy_of_nodelist(tNlist *elem);
int count_elements_nodelist(tNlist *list);
tNlist *insertnodelist_into_nodelist_after(tNlist *elem, tNlist *list);
tNlist *insertnodelist_into_nodelist_before(tNlist *elem, tNlist *list);
tNlist *replace1_in_nodelist(tNlist *elem, tNlist *list);
tNlist *first_replace1_in_nodelist(tNlist *elem, tNlist *list);
tNlist *remove1_in_nodelist(tNlist *elem, int return_next);
tNlist *first_nodelist(tNlist *list);
tNlist *last_nodelist(tNlist *list);
void free_nodelist(tNlist *elem);
long update_mesh_myln_node_nid(tMesh *mesh);
long get_node_nid(tNode *node);
tNlist *append_nodelist_to_mesh_lns_myln(tMesh *mesh, tNlist *list);
tNlist *replace1_in_mesh_lns_myln(tNlist *elem, tNlist *nlist);
tNlist *make8children_in_mesh_lns_myln(tNlist *elem, int n[3]);
void destroy8siblings_in_mesh_lns_myln(tNlist *sib);
void realloc_datvariables(tDat *dat, int nv_new);
void realloc_meshvariables(tMesh *mesh, int nvdb_new);
void enablevarcomp_innode(tNode *node, int i);
void disablevarcomp_innode(tNode *node, int i);
void enablevar_innode(tNode *node, int i);
void disablevar_innode(tNode *node, int i);
void enablevar_inpatch(tPat *pat, int i);
void disablevar_inpatch(tPat *pat, int i);
void enablevar(tMesh *mesh, int i);
void disablevar(tMesh *mesh, int i);
void enablevarlist(tVarList *vl);
void disablevarlist(tVarList *vl);

/* array.c */
void mm_array_indir(tArray *Ata, tArray *Ba, int dir, tArray *ABa);
void mm_array0(tArray *Aa, tArray *Ba, tArray *ABa);
void mm_array1(tArray *Aa, tArray *Ba, tArray *ABa);
void mm_array2(tArray *Aa, tArray *Ba, tArray *ABa);
void set_const_array(tArray *A, double c);
void copy_array_plane(tArray *A, int dir, int pA, tArray *P, int pP);

/* print.c */
void printmesh(tMesh *g);
void printpatch(tPat *box);
void printnode(tNode *n);
void printnode_and_neighbors(tNode *n);
void printnodelist_and_neighbors(tNlist *nl);
void printnodelist(tNlist *nl);
void printvar_innode(tNode *node, int vi);
void printarray(tArray *A);
void printarray_matrix0(tArray *A);
void printarray_matrix1(tArray *A);
void printarray_matrix2(tArray *A);
//void printbface(tBface *bface);
//void printbfaces(tPat *pat);

/* surface.c */
void init_all_myln_surfaces(tMesh *mesh);
void set_all_myln_mysurf(tMesh *mesh);
void request_all_myln_surfaces_exchange(tMesh *mesh);
void get_all_surfaces(tNode *node);
void get_all_myln_surfaces(tMesh *mesh);
