
/* info that describes one leaf node (or element) */
//#define LOCSMAX 128
//typedef struct tEL {
//  int p;                  /* patch number */
//  int l;                  /* refinement level of this node */
//  char loc[LOCSMAX];      /* node location string, giving loc. in patch */
//} tEl;

//in src/main/amr/storage.c we may also want:
//@@ -1643,6 +1643,12 @@ tDat *alloc_dat(tNode *node)
//   dat = calloc(1, sizeof(tDat));
//   if(!dat) errorexit("out of memory for dat");
//
//+  /* set patch number and node location string in dat->info->nodeloc */
//+  dat->info->p = node->pat->p;
//+  dat->info->l = node->l;
//+  node_location_str(node, dat->info->loc, LOCSMAX);
//+
//+  /* set rest of dat: */
//   dat->node = node;
//   dat->nv = nv;
//   if(nv==0) return dat;



//NOTE: we need wolfGIT/c/linux_list.h


/* location of an element */
#define LOCSMAX 128
typedef struct tELOC {
  int p;                  /* patch number */
  int l;                  /* refinement level of this node */
  char loc[LOCSMAX];      /* node location string, giving loc. in patch */
} tEloc;


/* a leaf node or element */
typedef struct tELM {
  struct list_head list;  /* all elms form a linked list */
  tEloc eloc[1];          /* elm location */

  double dt;              /* time step in node */
  double time;            /* current time in node */
  int n[3];               /* number of points in X,Y,Z-directions */
  int rflag;              /* flag for refining node */
  long nid;               /* node ID, updated by update_mesh_myln_node_nid */
  int pt_typ[3];          /* e.g. pt_typ[1]=P_LGL => LGL in dir1 of node */
  int datrank;            /* rank of proc that rightfully has data */
  tDat *dat;              /* pointer to data (NULL if not on this proc) */

} tElm;



/* several patches and thus a list of leaf nodes make up the
   computational mesh */
typedef struct tMESH {
  double dt;         /* time step */
  double time;       /* current time */
  int iteration;     /* current iteration number */
  tTodo *skel[NFUNCBINS]; // list of tTodo's from skeleton.c
  int nvdb;          /* number of variables */
  struct tVAR *vdb;  /* variable data base */
  int vdb_iStart;    /* index we start at when searching for a var */
  int npdb;          /* number of mesh parameters */
  struct tPAR *pdb;  /* parameter data base */
  int pdb_iStart;    /* index we start at when searching for a par */
  tEvoSys evosys[1]; /* contains lists of VarLists and RHS for evolve */
  int npats;         /* number of patches */
  tPat **pat;        /* list of pointers to patches */
  tNlist *lns;       /* start of linked list of all leaf nodes */
  long nln;          /* total number of leaf nodes */
  tMylnodes myln[1]; /* elements of lns owned by this proc */
  tMUTEX mutex[1];   /* mutex for mesh */

  /* newamr stuff */
  struct list_head myelm_head; /* list head for elms on this proc */
  long  nmyelm;      /* number of elms on this proc */
  tElm **myelm;      /* list of pointers to elms on this proc */
                     /* myelm and myelm_head list are copies of each other */
  long  nnbelm;      /* number of nb elms on other procs */
  tElm **nbelm;      /* list of pointers to nb elms on other procs */

} tMesh;
/* NOTE: the list lns needs to be distributed among MPI jobs:
use space filling curve as in
http://www.speedup.ch/workshops/w42_2013/carsten.pdf
*/
