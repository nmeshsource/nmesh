
/* declarations from main that we need here already */
#include "../main/skeleton.h"


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
  int nv;               /* number of vars */
  int nvenabled;        /* number of enabled vars */
  struct tARRAY **v;    /* list of data pointers to 3d vars */
                        // if v[i]=NULL var v[i] and its ghosts are not enabled
  struct tARRAY **g[6]; /* list of data pointers to 2d ghost zones */
        // g[0]=ghosts in -X dir, g[3]=ghosts in +Y dir
//??? //  struct tNODE *node;   /* pointer to node where dat is in */
} tDat;


/* a node */
typedef struct tNODE {
  struct tPAT *pat;       /* pointer to patch that contains node */
  struct tNODE *parent;   /* pointer to parent node */
  struct tNODE *child[8]; /* list of pointers to childeren nodes */
  struct tNODE *nb[6][5]; /* neighbs in +/-X,Y,Z dir: nb[+-dir][neib.-index] */
             // nb[4][0]= 1st neighb in +Z dir, nb[1][4]= 3rd neighb in +X dir
             // nb[dir][k] is 0 terminated, i.e. 0 for one k in {0,1,2,3,4}
             // if e.g. nb[3][1]=0 there is no 2nd neighb. in +Y dir.
  double bbox[6];         /* bounding box (in X,Y,Z) of this node */
  int patface[6];         /* whether node is at patch face 0,1,2,3,4,5 */
  int n[3];               /* number of points in X,Y,Z-directions */
  int np;                 /* np = n[0] * n[1] * n[2]; */
  int l;                  /* refinement level of this node */
  int leaf;               /* is 1 if this is a leaf node */
  int ijk;                /* node index (0-7), i.e. child number wrt. parent */
  struct tARRAY *D[3];    /* differentiation matrix in all 3 dirs for [-1,1]
                             domain. This just points to an array in patch. */
  tDat *dat;              /* pointer to data (NULL if not on this proc) */
  int datrank;            /* rank of proc that rightfully has data */
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
  struct tARRAY ***D;   /* list of differentiation matrices */
  int nD;               /* number of diff matrices stored */
  tNlist *lns;          /* start of linked list of leaf nodes in this patch */
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
  int npats;        /* number of patches */
  int nvdb;         /* number of variables */
  int npdb;         /* number of mesh parameters */
  tTodo *skel[NFUNCBINS]; // list of tTodo's from skeleton.c
  struct tVAR *vdb; /* variable data base */
  struct tPAR *pdb; /* parameter data base */
  int pdb_iStart;   /* index we start at when searching for a par */
  int vdb_iStart;   /* index we start at when searching for a var */
  tPat **pat;       /* list of pointers to patches */
  tNlist *lns;      /* start of linked list of all leaf nodes */
} tMesh;
/* NOTE: the list lnodes needs to be distributed among MPI jobs:
use space filling curve as in
http://www.speedup.ch/workshops/w42_2013/carsten.pdf
*/


/***********************************************************************/
/* other useful objects */
/***********************************************************************/
/* arrays */
typedef struct tARRAY {
  int n[3];    /* dims in all 3 dirs */
  int N;       /* N = n[0] * n[1] * n[2]; */
  double *a;   /* pointer to double data (could add one more for GPU data) */
//  void *Owner; /* pointer to patch or node this array belongs to */
//  int tOwner;  /* type of owner, e.g. NODE OR PAT */
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
tArray *alloc_array(int n[3]); //, void *Owner, int tOwner);
void free_array(tArray *array);
tMesh *alloc_mesh(int npats);
void realloc_mesh_patches(tMesh *mesh, int npats);
void free_mesh(tMesh *mesh);
tPat *alloc_patch(tMesh *mesh, int p, int nD);
void free_patch(tPat *pat);
tNode *alloc_node();
void free_node(tNode *node);
tDat *alloc_dat(int nv);
void free_dat(tDat *dat);
tNlist *alloc_nodelist(tNode *node);
tNlist *addnode_to_nodelist(tNlist *elem, tNode *node);
tNlist *addnodelist_to_nodelist(tNlist *elem, tNlist *list);
tNlist *replace1_in_nodelist(tNlist *elem, tNlist *list);
tNlist *remove1_in_nodelist(tNlist *elem);
void free_nodelist(tNlist *elem);
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

/* print.c */
void printmesh(tMesh *g);
void printpatch(tPat *box);
void printnode(tNode *n);
//void printbface(tBface *bface);
//void printbfaces(tPat *pat);
void printvar_innode(tNode *node, char *name);
void printVarList(tVarList *vl);


