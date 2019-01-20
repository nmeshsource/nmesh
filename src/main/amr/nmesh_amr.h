
/* declarations from main that we need here already */
#include "../main/skeleton.h"


/* the data within a node, this should be only on one proc */
typedef struct tDAT {
  double **v;           /* list of data pointers to 3d variables */
  double **g[6];        /* list of data pointers to 2d ghost zones */
        // g[0]=ghosts in -X dir, g[3]=ghosts in +Y dir
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
  int l;                  /* refinement level of this node */
  int leaf;               /* is 1 if this is a leaf node */
  int i;                  /* node index (0-7) */
  tDat *dat;              /* pointer to data (NULL if not on this proc) */
} tNode;


/* the nodes fill a patch */
typedef struct tPAT {
  struct tMESH *mesh;   /* pointer to mesh that contains patch */
  int p;                /* index of this patch */
  /* funcs to compute X,Y,Z from x,y,z and vice versa: */
  int (*XYZ_Of_xyz)(struct tPAT *pat, double x, double y, double z, double *X, double *Y, double *Z);  /* func to compute X,Y,Z from x,y,z */
  int (*xyz_Of_XYZ)(struct tPAT *pat, double X, double Y, double Z, double *x, double *y, double *z);  /* func to compute x,y,z from X,Y,Z */
  double bbox[6];   /* bounding box (in X,Y,Z) of this patch */
  tNode *rnode;     /* root node in this patch */
  int np;           /* np = n[0] * n[1] * n[2]; */
  double *D[3];     /* differentiation matrix in all 3 dirs for [-1,1] domain */
} tPat;


/* a linked list of nodes */
typedef struct tNLIST {
  tNode *node;
  struct tNLIST *next;
  struct tNLIST *prev;
} tNlist;


/* several patches and thus a list of leaf nodes make up the 
   computational mesh */
typedef struct tMESH {
  double dt;        /* time step */
  int iteration;    /* current iteration number */
  int npatches;     /* number of patches */
  int nvdb;         /* number of variables */
  int npdb;         /* number of mesh parameters */
  tTodo *skel[NFUNCBINS]; // list of tTodo's from skeleton.c
  struct tVAR *vdb; /* variable data base */
  struct tPAR *pdb; /* parameter data base */
  tNlist *lnodes;   /* linked list of leaf nodes */
  tPat **pat;       /* list of pointers to patches */
} tMesh;

