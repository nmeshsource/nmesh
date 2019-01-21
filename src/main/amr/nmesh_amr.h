
/* declarations from main that we need here already */
#include "../main/skeleton.h"


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
  struct tNODE *rnode;    /* pointer to root node */
  struct tNODE *parent;   /* pointer to parent node */
  struct tNODE *child[8]; /* list of pointers to childeren nodes */
  struct tNODE *nb[6][4]; /* neighbors in +/-X,Y,Z dirs. */
             // nb[0][0]= 1st neighb. in +X dir, nb[1][4]= 3rd neighb. in +X dir
             // if e.g. nb[3][1]=0 there in no 2nd neighb. in +Y dir.
  double bbox[6];         /* bounding box (in X,Y,Z) of this node */
  int n[3];               /* number of points in X,Y,Z-directions */
  int np;                 /* np = n[0] * n[1] * n[2]; */
  struct tARRAY *D[3];    /* differentiation matrix in all 3 dirs for [-1,1]
                             domain. This just points to an array in patch. */
  int l;                  /* refinement level of this node */
  int leaf;               /* is 1 if this is a leaf node */
  int i;                  /* node index (0-7) */
  tDat *dat;              /* pointer to data (NULL if not on this proc) */
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
  struct tARRAY **D;    /* list of differentiation matrices */
  int nD;               /* number of diff matrices stored */
  tNlist *ln;           /* linked list of leaf nodes in this patch */
} tPat;



/* several patches and thus a list of leaf nodes make up the 
   computational mesh */
typedef struct tMESH {
  double dt;        /* time step */
  int iteration;    /* current iteration number */
  int npats;        /* number of patches */
  int nvdb;         /* number of variables */
  int npdb;         /* number of mesh parameters */
  tTodo *skel[NFUNCBINS]; // list of tTodo's from skeleton.c
  struct tVAR *vdb; /* variable data base */
  struct tPAR *pdb; /* parameter data base */
  tPat **pat;       /* list of pointers to patches */
  tNlist *ln;       /* linked list of all leaf nodes */
} tMesh;
/* NOTE: the list lnodes needs to be distributed amoung MPI job:
use space filling curve as in
http://www.speedup.ch/workshops/w42_2013/carsten.pdf
*/


/***********************************************************************/
/* other useful objects */
/***********************************************************************/
/* owner types */
enum
{
  NODE,
  PAT,
  MESH,
  NOWNERS
};

/* arrays */
typedef struct tARRAY {
  int n[3];    /* dims in all 3 dirs */
  int N;       /* N = n[0] * n[1] * n[2]; */
  double *a;   /* pointer to double data (could add one more for GPU data) */
  void *Owner; /* pointer to patch or node this array belongs to */
  int tOwner;  /* type of owner, e.g. NODE OR PAT */
} tArray;


/* variable lists */
typedef struct tVARLIST {
  double time;
  int n;
  int *index;
  void *Owner;  /* pointer to patch or node this array belongs to */
  int tOwner;   /* type of owner, e.g. NODE OR PAT */
  void *vlPars; /* A pointer that is usually NULL, but can point to some
                   object that contains special extra pars or info. This
                   pointer is not touched by the funcs in variables.c (such
                   as vlduplicate, vlcopy, vlfree, ...). So the user has 
                   to manage it: e.g. free it, before calling vlfree. */
} tVarList;


/**************************************************************************/
/* functions */
/***********************************************************************/

/* mesh.c */
tMesh *make_empty_mesh(int pr);
tMesh *alloc_mesh(void);


/* storage.c */
tArray *alloc_array(int n[3], void *Owner, int tOwner);
void free_array(tArray *array);
tMesh *alloc_mesh(void);
void free_mesh(tMesh *mesh);
tPat *alloc_patch(tMesh *mesh, int p, int nD) ;
void free_patch(tPat *pat) ;
tNode *alloc_node();
void free_node(tNode *node) ;
tDat *alloc_dat(int nv);
void free_dat(tDat *dat);
void realloc_dat_varlists(tDat *dat, int nv_new);


/* print.c */
void printmesh(tMesh *g);
void printpatch(tPat *box);
void printnode(tNode *n);
//void printbface(tBface *bface);
//void printbfaces(tPat *pat);
void printvar_innode(tNode *node, char *name);
void printVarList(tVarList *vl);


