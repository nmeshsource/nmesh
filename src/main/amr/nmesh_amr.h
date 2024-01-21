/* nmesh_amr.h */
/* (c) Wolfgang Tichy 2/2019 */

#include <time.h>    /* for struct timespec */

/* declarations from other parts that we need here already */
#include "thread_defs.h"
#include "refine.h"
#include "../main/linux_list.h"
#include "../main/khash.h"
#include "../main/skeleton.h"
#include "../main/variables.h"
#include "../nMPI/nMPI_defs.h"
#include "../evolve/evosys.h"
#include "../coordinates/CI.h"


/* Main parts of a mesh:

tMesh:    ---------------------mesh-----------------------------------------
tPat:     |-----patch0-----|-------patch1-------|--patch2----|...
tElm:      elm elm ...      elm elm elm ...      elm elm  ...
                ^
  The elements (elm) shown here are our computational subdomains in each
  patch.

  Each elm has a tDat struct that can be empty if the data is on another
  proc. The tDat struct contains lists of arrays, one for each variable.

  The elms are the leaf nodes of an oct-tree. Here we show the structure
  in one patch starting from a root-node. Each part of the tree is a node.
  The leaf nodes are marked as elm. We only show 2 children instead of 8:
level
l=0                  _________root-node______________________
l=1          ______node_________                     ______node______
l=2     __node__        |------elm-------|       _node__     |-------elm-----|
l=3 |-elm---|   _node_                   |---elm---|---elm---|
l=4         |-elm-|-elm-|
  We have one tree per patch.
  The ends of the tree are called leaf nodes (marked as elm here).
  The elms cover our entire computational domain. Since the elms are nodes
  themselves we sometimes call them nodes as well.
*/

/* The leaf nodes are the nodes we do all our calculations in and thus the
   the most important nodes. In fact it may be good to keep only the leaf
   nodes to save memory.
   We also call a leaf nodes an element or elm for short. */

/* abbreviation for unsigned long */
typedef unsigned long ulong;

/* new element IDs (eid)s are set to this value */
#define  EID_INVALID  (ULONG_MAX)

/* location of an element (or elm) */
#define NPBYTES 21
#define NLOCS   ((8*NPBYTES)/3)
typedef struct tELOC {
  int p;              /* patch number */
  int l;              /* refinement level of this node */
  char loc[NLOCS];    /* elm location string, giving loc. in patch */
  unsigned long eid;  /* elm ID, updated by update_mesh_myelms_elm_eid_dt */
} tEloc;

/* location of an element (or elm) in packed form
   NPBYTES can be 13,21,29,37 for optimal sizeof(tEploc) */
typedef struct tEPLOC {
  unsigned short int p; // patch number (in 2 bytes)
  unsigned char l;      // refinement level of this node (in 1 byte)
  unsigned char ploc[NPBYTES];  // packed loc (3 bits per level)
  unsigned long eid;    // elm ID, updated by update_mesh_myelms_elm_eid_dt
} tEploc;

/* a leaf node or element called elm */
/* Beginning of tElm */
#define ELMHEADER \
  tEploc eploc[1];        /* elm location in packed form */ \
  double dt;              /* time step in node */ \
  double time;            /* current time in node */ \
  double bbox[6];         /* bounding box (in X,Y,Z) of this node */ \
  int n[3];               /* number of points in X,Y,Z-directions */ \
  int np;                 /* np = n[0] * n[1] * n[2]; */ \
  int rflag;              /* flag for refining node */ \
  int pt_typ[3];          /* e.g. pt_typ[1]=P_LGL => LGL in dir1 of node */ \
  int datrank;            /* rank of proc that rightfully has data */
typedef struct tELM {
  ELMHEADER
  /* stuff below this line is not copied when elm is sent to another rank */
  struct tDAT *dat;       /* pointer to data (NULL if not on this proc) */
  //nMPI_Comm comm;         // MPI_comm for node, could contain only ranks
                            // where dat is and where all neighb. have dat
  //struct tMESH *mesh;     // pointer to mesh that contains elm
  struct tPAT *pat;       // remove one day, and replace by mesh
  /* items to do with neighbors: */
  int nfnb[6];            /* number of face neighbor nodes (fnb below) */
  struct tELM **fnb[6];   /* list of neighbor nodes on face, contains info
                             condensed out of nfaces */
  struct list_head list;  /* all elms form a linked list */
  //ulong oid;              /* old elm ID */
  //char loc[NLOCS];        /* unpacked elm location */
} tElm;

/* data type for only header part of tElm */
typedef struct tELM0 {
  ELMHEADER
} tElm0;

/* the old node type tNode is now a leaf node of type tElm */
typedef tElm tNode;


/* in case we need it, we can also make more linked element lists */
/* a linked list of generic pointers */
typedef struct tGLIST {
  void *entry; /* could be entry=elm */
  struct list_head list;
} tGlist;
#define glist_entry(ptr)       list_entry(ptr, tGlist, list)->entry
#define glist_first_entry(ptr) list_first_entry(ptr, tGlist, list)->entry
#define glist_last_entry(ptr)  list_last_entry(ptr, tGlist, list)->entry


/* extra info about node state that has nothing to do with neighbor info
   or connectivity */
typedef struct tNODEINFO {
  int evo_troubled;       /* is 1 if node was troubled during RK substep */
  int trbl_score;         /* trouble score in node (e.g. due to shocks) */
  tRef trbl_ref[1];       /* refinement we use if we act on trouble score */
  int use_fv;             /* whether we currently use fin. vol. in node */
  int nlim;               /* number of consectutive evo substeps in which
                             limiter was active */
  int load_timer_running;        //1 if node load timer is running, else 0
  struct timespec load_start[1]; //time when node load timer was started
  double load_TimeIn_s;          //time spent on this node in seconds
  int desrank;            /* rank that should have this node/elm */
  int nnbinfo[6];         /* num. of eplocs in amr_elm_nbinfo, if nnbinfo<0
                             amr_elm_nbinfo needs to be updated */
} tNodeInfo;

/* the data within a node, this should be only on one proc */
typedef struct tDAT {
  struct tELM *node;//FIXME call it elm /* pointer to node dat is in */
  int nv;                 /* number of vars */
  int nvenabled;          /* number of enabled vars */
  struct tARRAY **v;      /* list of data pointers to vars, if v[vi]=NULL,
                             the var vi and its surfaces are not enabled */
  struct tSURFACE **s[6]; /* list of surfaces needed for data exchange,
                             e.g. s[0]=surfs in -X dir, s[3]=surfs in +Y dir,
                             if s[6][vi]=NULL var vi does not need exchange */
  int surfs_set;          /* whether all nb surfaces are set (via MPI) */
  int coords_set;         /* whether coord vars are enabled and set */
  struct tCOM *com[6];    /* com for each face */
  struct tINDIC **ic;     /* indicators such as min/max for each var,
                             e.g. ic[vi] has indicators for var vi  */
  struct tCOM *icom;      /* com for indc */
  struct tCOM *gcom;      /* com for ghosts */
  tNodeInfo info[1];      /* extra info about node (sent during loadbal.) */
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
  struct tARRAY *ajsurf;  /* vals from adjacent surf. interpol. to my points */
  int allocd_ajsurf;      /* 1 if we need to free ajsurf */
} tSurface;
/* NOTE:
   mysurf comes from this proc,
   nbsurf[i] can just point if nb[i] is local, otherwise we need to alloc */

/* indicator data for node to neighbor node communication */
typedef struct tINDIC {
  struct tDAT *dat;          /* pointer to dat the indicator is in */
  int vi;                    /* var index */
  int nvals;                 /* number of indic. vals, e.g. 2 for min,max */
  struct tARRAY *myindc;     /* array that contains values for my var */
  struct tARRAY **nbindc[6]; /* list of values from neighb. indicators, e.g.
                                nbindc[f][i] is indic of nb i on face f */
  int *allocd_nbindc[6]; // allocd_nbindc[f][i]=1 if need to free nbindc[f][i]
} tIndic;
/* NOTE:
   myindc comes from this proc,
   nbindc[i] can just point if nb[i] is local, otherwise we need to alloc */


/* the nodes fill a patch */
typedef struct tPAT {
  double bbox[6];       /* bounding box (in X,Y,Z) of this patch */
  double bbdiag;        /* length (in X,Y,Z-coords) of 3D-diagonal in bbox */
  int p;                /* index of this patch */
  struct tMESH *mesh;   /* pointer to mesh that contains patch */
  /* funcs to compute x,y,z from X,Y,Z and vice versa: */
  int (*xyz_of_XYZ)(struct tPAT *pat, tNode *node, int ind,
                    const double X[3], double x[3]); /* func to compute x,y,z from X,Y,Z */
  int (*XYZ_of_xyz)(struct tPAT *pat, tNode *node, int ind,
                    double X[3], const double x[3]); /* func to compute X,Y,Z from x,y,z */
  int (*dXYZ_dxyz)(struct tPAT *pat, tNode *node, int ind,
                   const double X[3], double x[3], double dXYZdxyz[3][3]);
  tCoordInfo CI[1];     /* info about coords, access e.g. as: pat->CI->xc[1] */
  int periodic[3];      /* if e.g. periodic[0]=1, patch is periodic in dir0 */
  struct tBFACE *bfaces[6]; /* 1st bface of this patch on each face */
  //tElm0 rnode[1];     //FIXME: remove!  /* root node in this patch */
} tPat;
/* Note: each patch has Bfaces as in sgrid. But instead of pointlists we use
   bounding rectangles in both adjacent bfaces. These rectangles (brct) are
   not exact but simply contain some of the region where two faces meet. In
   case the two faces are two rectangles brct is the excat intersection of
   the face rectangles. But this does not happen with e.g. the 38 cubed
   spheres from sgrid, because of the extended ranges in the A-coord lead to
   overlapping faces that are not simple rectangles. So when we connect
   nodes and set node->fnb we use bfaces only to see if both nodes are in a
   patch that share bfaces on the two node faces. */


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

  /* newamr stuff */
  struct list_head myelm_head; /* list head for elms on this proc */
  ulong nmyelm;      /* number of elms on this proc */
  tElm **myelm;      /* list of pointers to elms on this proc */
                     /* myelm and myelm_head list are copies of each other */
  ulong nnbelm;      /* number of nb elms on other procs */
  tElm **nbelm;      /* list of pointers to nb elms on other procs */
  ulong *eidlim;     /* (last eid on rank rk) = eidlim[rk]-1 */
} tMesh;
/* NOTE: the list myelm is distributed among MPI jobs:
use space filling curve as in
http://www.speedup.ch/workshops/w42_2013/carsten.pdf
*/


/* Currently unused struct. Each entry is one element of the leaf node
   list. The leaves are sorted into categories that should be based on the
   time it takes to process each in it. This can help with load
   balancing. */
typedef struct tMYLNODES {
  int nncats;      /* number of leaf node categories */
  int *ncat;       /* ncat[c] is number of leaves in category c */
  int nm;          /* max of all ncat[c] */
  tElm ***ln;    /* ln[c][i] is leaf i of category c on this proc */
} tMylnodes;


/***********************************************************************/
/* Bfaces */
/***********************************************************************/
/* Note: To include info about which patches touch, tPat also contains
         info about BCs for patches:
  struct tBFACE *bfaces[6];
where tBface is a part of a patch face that touches at most one
other patch. We use the same BC on all of tBface. */
typedef struct tBFACE {
  tPat *pat;      // patch in which our patchface is
  int f;          // face, runs from 0 to 5 (for each pat)
   // The normal vector is n^i_{a}=dX^a/dx^i, e.g. X^1=const face has n^i_{1}
  double brct[4]; // bound. rectangle of bface in the 2 coords perp. to face f
  int brct_isset; // whether bounding rectangle brct is set
  int op;         // ind. of other pat that touches or overlaps, -1 if none
  struct tBFACE *obface; // pointer to other bface that touches
  int ioC0_0;     // ind of vars in this node that contain coords in other pat
  int face2;      // 1 if we set normal derivs of field and not field itself
  int boundary;   // type of boundary: NOTBOUND=0, OUTERBOUND, INNERBOUND
  struct tBFACE *next; // next bface in this patch
  struct tBFACE *prev; // previous bface in this patch
} tBface;
/* NOTE: ioC0_0 is set to -1 when a bface is allocated with
   add_empty_bface. So -1 means we do not know it yet. */

/* possible values of bface->boundary */
enum
{
  NOTBOUND=0,
  OUTERBOUND=1, // bface is outer mesh boundary (e.g. infinity)
  INNERBOUND    // bface is inner boundary (e.g. horizon)
};


/***********************************************************************/
/* a single point */
/***********************************************************************/
/* all the info about one grid point */
typedef struct tPOINT{
  tNode *node;   // node we are on
  int ijk;       // index of point in node
} tPoint;


/***********************************************************************/
/* other useful objects */
/***********************************************************************/
/* arrays */
typedef struct tARRAY {
  int n[3];     /* dims in all 3 dirs */
  int N;        /* N = n[0] * n[1] * n[2]; */
  int Ne;       /* extra space in array beyond the N entries (usually 0) */
  size_t size;  /* size of allocated space in bytes */
  union {       /* anon. union with host data (add one more for GPU data) */
    double *d;  /* pointer to double data */
    long *l;    /* pointer to long data using same mem as double data */
    ulong *ul;  /* ptr to unsigned long data using same mem as double data */
    tEploc *eploc; /* ptr to tEploc data using same mem as double data */
    char *c;    /* pointer to char data using same mem as double data */
    int *i; };  /* pointer to int data using same mem as double data */
  int d_nofree; /* d_nofree=1 if free_array should not free d */
  int ns;       /* number of segments (usually 1) */
  int si;       /* segment index (usually 0) */
  int info;     /* space to save extra info, e.g. MPI req. numbers or tags */
  int *range[2]; // d[i+n[0]*j]!=0 only if range[0][j] <= i < range[1][j]
  void *par;    /* pointer to some extra pars (usually NULL) */
} tArray;
/* NOTE: the anon. union is used to be able to store either double or int in
   tArray */


/* variable lists */
typedef struct tVARLIST {
  int n;
  int *index;
  tMesh *mesh;  /* pointer to mesh to which vars belong */
  double time;
  void *vlPars; /* A pointer that is usually NULL, but can point to some
                   object that contains special extra pars or info. This
                   pointer is not touched by the funcs in variables.c (such
                   as vlduplicate, vlcopy, vlfree, ...). So the user has
                   to manage it: e.g. free it, before calling vlfree. */
} tVarList;


/**************************************************************************/
/* define some khash table types */
/**************************************************************************/

/* hash sets with only 32 or 64 bit ints as keys and no values */
KHASH_SET_INIT_INT(u32)
KHASH_SET_INIT_INT64(u64)
/* hash table with 32 bit ints as keys and generic pointer values */
KHASH_MAP_INIT_INT(u32_gptr, void *)


/**************************************************************************/
/* for load balancing */
/**************************************************************************/
//...

/**************************************************************************/
/* for results from timing.c */
/**************************************************************************/
typedef struct tTIMING {
  double mm1_speed;
  double mm_speed;
  double myops;  /* myops = speed * myT */
  double ops0;   /* ops0 = \sum_{r=0}^{myrank-1} myops(rank=r) */
  double allops; /* allops = \sum_{r=0}^{size-1} myops(rank=r) */
} tTiming;


/**************************************************************************/
/* useful defines */
/**************************************************************************/
#include "nmesh_amr_defs.h"

/**************************************************************************/
/* loops */
/**************************************************************************/
#include "nmesh_amr_loops.h"

/**************************************************************************/
/* amr parameters */
/**************************************************************************/

/* structure that holds global amr vars and pars */
typedef struct {
  int dir_active[3];   /* Par("amr_dir_active0"), ... */
  int sibl1to7_weight; /* Par("amr_sibl1to7_weight") */
  int MPIexchange;     /* Par("amr_MPIexchange") */
  int nghosts;         /* Par("amr_nghosts") */
  int elm_nbinfo0;     /* Ind("amr_elm_nbinfo0") */
} tAMR;


/**************************************************************************/
/* functions */
/**************************************************************************/

/* mesh.c */
tMesh *make_empty_mesh(int pr);
int amr_use_fv_if_P_UNIFORM(tMesh *mesh);
tPat *add_patch_without_rnode(tMesh *mesh, double bbox[6]);
tPat *add_patch(tMesh *mesh, double bbox[6],
                int *pt_typ_root, int nroot[3], int datrank);
int amr_setup_mesh(tMesh *mesh);
int amr_set_all_bfaces(tMesh *mesh);
int amr_set_bfaces_and_rnode_nbinfo_fnb(tMesh *mesh, int pr);

/* storage.c */
tArray *alloc_empty_array_with_segs(int n[3], int Ne, int ns);
tArray *alloc_array_with_segs(int n[3], int Ne, int ns);
tArray *alloc_array1d_with_segs(int N, int Ne, int ns);
tArray *alloc_array(int n[3]);
tArray *alloc_array1d(int N);
tArray *alloc_array2d(int n0, int n1);
tArray *alloc_empty_array2d(int n0, int n1);
tArray *get_array_seg(tArray *array, int si);
void point_array_d_to_data(tArray *array, void *data, int nofree);
int redimension_array(tArray *array, int n[3]);
int redim_array(tArray *array, int n0, int n1, int n2);
void alloc_2darray_irange_of_j(tArray *array);
void free_array(tArray *array);
void free_3_arrays(tArray *array[3]);
void *memcpy_from_array(const tArray *ar, size_t bytestride, size_t pos,
                        void *dest, size_t n);
void *memcpy_to_array(tArray *ar, size_t bytestride, size_t pos,
                      const void *src, size_t n);
void *memcpy_to_array_redim(tArray *ar, size_t bytestride, size_t pos,
                            const void *src, size_t n);
tElm *alloc_elm(tMesh *mesh);
tElm *alloc_elm_init_pat(tMesh *mesh, int p);
tElm *alloc_elm_of_elmheader(tMesh *mesh, tElm0 *elmheader);
tElm *alloc_elm_of_eploc(tMesh *mesh, tEploc *eploc);
void free_elm(tElm *elm);
tMesh *alloc_mesh(int npats);
void realloc_patlist_in_mesh(tMesh *mesh, int npats);
void free_mesh_patches_and_nodes(tMesh *mesh);
void free_mesh_contents_exceptMeshFuns(tMesh *mesh);
void free_all_mesh_contents(tMesh *mesh);
tPat *alloc_patch(tMesh *mesh, int p);
void free_patch(tPat *pat);
tNode *update_node_n_pt_typ_return_node_old(tNode *node, int *n, int *pt_typ);
void update_node_n_pt_typ_free_node_old(tNode *node, tNode *node_old);
void update_node_n_pt_typ_restore_from_node_old(tNode *node, tNode *node_old);
void update_node_n_pt_typ(tNode *node, int *n, int *pt_typ);
tDat *alloc_dat(tNode *node);
void free_dat(tDat *dat);
tElm *replace_parent_by_8children(tElm *parent, int n[3], int pt_typ[3]);
ulong update_mesh_myelms_elm_eid_dt(tMesh *mesh);
int calc_node_lid(tNode *node);
ulong calc_elm_lid(tNode *node);
ulong calc_local_elm_id(tNode *elm);
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
#define update_mesh_myln_node_nid(mesh) update_mesh_myelms_elm_eid_dt(mesh)

/* array.c */
void mm_array_indir(tArray *Ata, tArray *Ba, int dir, tArray *ABa);
void mm_array0_norestrict(tArray *Ata, tArray *Ba, tArray *ABa);
void mm_array0(tArray *Aa, tArray *Ba, tArray *ABa);
void mm_array1(tArray *Aa, tArray *Ba, tArray *ABa);
void mm_array2(tArray *Aa, tArray *Ba, tArray *ABa);
void set_const_array(tArray *A, double c);
void copy_array_data(tArray *Src, tArray *Dest);
void copy_array_plane(tArray *A, int dir, int pA, tArray *P, int pP);
void copy_array_planes(int np, tArray *A, int dir, int pA,
                       tArray *P, int pP);
double Lp_norm_array(tArray *A, double p);
void array_diff(tArray *D, tArray *A, tArray *B);
void array_reldiff(tArray *D, tArray *A, tArray *B);
double Lp_norm_array_diff(tArray *A, tArray *B, double p);
double Lp_norm_array_reldiff(tArray *A, tArray *B, double p);
double max_array(tArray *A, int *ind);
double min_array(tArray *A, int *ind);

/* print.c */
void printmesh(tMesh *g);
void printpatch(tPat *box);
void printCI(tPat *pat);
void printeloc(const tEloc *eloc);
void printeloc_s(const tEloc *eloc, const char *s);
void printeploc(const tEploc *eploc);
void printeploc_s(const tEploc *eploc, const char *s);
void printelm0(const tElm0 *e, const char *s);
void printelm(const tElm *e);
void printelmarray(long nelms, tElm **elm);
void printelmlist(struct list_head *elm_head);
void printelmglist(struct list_head *elm_head);
void printmyelms(tMesh *mesh);
void printnbelms(tMesh *mesh);
void print_amr_elm_nbinfo(const tElm *elm, int face);
void pr_nodename(tNode *node);
void printnodeinfo(const tElm *elm);
void printelm_nodeinfo(const tElm *elm);
void printvar_innode(tNode *node, int vi);
void printvar_ajsurfdiff(tNode *node, int vi);
void printvar_indc(tNode *node, int vi);
void print_matrices_innode(tNode *n);
void printarray(tArray *A);
void printarray_int(tArray *A);
void printarray_matrix0(tArray *A);
void printarray_matrix1(tArray *A);
void printarray_matrix2(tArray *A);
void printarray_eploc(tArray *A, int details);
void printthisbface(tBface *bface, const char *s);
void printbface(tBface *bface);
void printbfaces_on_f(tPat *pat, int f);
void printbfaces(tPat *pat);
void printallbfaces(tMesh *mesh);
void pr3v(const char *s, const double x[3]);
void prlarray(const char *s, long n, const long *ar);
void prbbox(double *bb, int dim);
void printcorners(tPat *pat);
void printfacecorners(tPat *pat, int  f);
void printref(tRef *ref);
void print_n_pt_typ(tNode *node);
void printpoint(tPoint *pt);
void printvar_atpoint(tPoint *pt, int vi);
void printvarlist_atpoint(tPoint *pt, tVarList *vl, const char *infostr);
void print_u32(khash_t(u32) *nbranks);
void print_ef(khash_t(u32_gptr) *ef);

/* surface.c */
int init_all_surfaces(tNode *node);
int set_all_mysurf(tNode *node);
void request_all_surfaces_exchange(tNode *node);
void get_all_surfaces(tNode *node);
void free_dat_reqs_after_Waitall_com_send(tNode *node);
void init_all_myln_surfaces(tMesh *mesh);
void free_all_myln_surfaces(tMesh *mesh);
void set_all_myln_mysurf(tMesh *mesh);
void request_all_myln_surfaces_exchange(tMesh *mesh);
void get_all_surfaces(tNode *node);
void get_all_myln_surfaces(tMesh *mesh);
int causeMPIprogress_all_myln_surfaces(tMesh *mesh, int sendrecv);
void free_all_myln_nbsurf_only(tMesh *mesh);
void init_all_vl_surfaces(tMesh *mesh, tVarList *vl);
void set_all_vl_mysurf(tNode *node, tVarList *vl);
void request_all_vl_surfaces(tNode *node, tVarList *vl);
void get_all_vl_surfaces(tNode *node, tVarList *vl);
void free_all_vl_surfaces(tNode *node, tVarList *vl);

/* load.c */
void simple_load_balance(tMesh *mesh);
void load_balance(tMesh *mesh, int strategy);
int load_balance_if_needed(tMesh *mesh);

/* amr.c */
int amr_init_global_pars(tMesh *mesh);
void MPIexchange_init_all_myln(tMesh *mesh);
void MPIexchange_set_all_myln_localdata(tMesh *mesh);
void MPIexchange_request_all_myln_data(tMesh *mesh);
void MPIexchange_get_all_myln_data(tMesh *mesh);
void MPIexchange_free_all_myln(tMesh *mesh);
int Ind_n_norm(int i, int j, int k, int n[3], int norm);
tArray *node_Xb(tNode *node, int dir);
void Xb3_pt_typ_n(int pt_typ[3], int n[3], tArray *Xb[3]);
void Xb3_n(tNode *node, int n[3], tArray *Xb[3]);
void node_Xb3(tNode *node, tArray *Xb[3]);
tArray *node_Wq(tNode *node, int dir);
void Wq3_n(tNode *node, int n[3], tArray *Wq[3]);
void node_Wq3(tNode *node, tArray *Wq[3]);
tArray *node_WL(tNode *node, int dir);
void WL3_n(tNode *node, int n[3], tArray *WL[3]);
void node_WL3(tNode *node, tArray *WL[3]);
tArray *node_Dt(tNode *node, int dir);
void Dt3_n(tNode *node, int n[3], tArray *Dt[3]);
void node_Dt3(tNode *node, tArray *Dt[3]);
tArray *node_Dpt(tNode *node, int dir);
tArray *node_Dmt(tNode *node, int dir);
tArray *node_At(tNode *node, int dir);
void At3_pt_typ_n(int pt_typ[3], int n[3], tArray *At[3]);
void At3_n(tNode *node, int n[3], tArray *At[3]);
void node_At3(tNode *node, tArray *At[3]);
tArray *node_St(tNode *node, int dir);
void St3_pt_typ_n(int pt_typ[3], int n[3], tArray *St[3]);
void St3_n(tNode *node, int n[3], tArray *St[3]);
void node_St3(tNode *node, tArray *St[3]);
double node_basis(tNode *node, int dir, int i, double x, int np);

/* bfaces.c */
tBface *first_bface_containing_point(tPat *pat, int f, double C[2]);
tBface *first_obface_of_bface_containing_point(tPat *pat, int f, double C[2]);
tBface *nbbface_of_bface_containing_point(tNode *nb,
                                          tPat *pat, int f, double C[2]);
int common_facepoints(tNode *node, int f, tNode *nb, int nb_f);

/* indicators.c */
void init_all_myln_myindc_for_vl(tMesh *mesh, tVarList  *vl, int nvals);
void free_indc(tIndic *ic);
void free_all_myln_indc_for_vl(tMesh *mesh, tVarList  *vl);
void request_all_myln_indc_exchange_for_vl(tMesh *mesh, tVarList  *vl);
void get_all_myln_indc_for_vl(tMesh *mesh, tVarList  *vl);

/* refine.c */
void hp_refine_set_n_pt_typ(tNode *pnode, tRef *ref, int *n, int *pt_typ);
int resolve_shocks_using_nlim(tMesh *mesh);
void hrefine_elms_if_rflag(tMesh *mesh, tRef *ref);
void prefine_elms_if_rflag(tMesh *mesh, tRef *ref);
void remove_elms_if_rflag(tMesh *mesh, tRef *ref);
void refine_set_rflag_forall_nodes(tMesh *mesh, int rflag);
void refine_set_use_fv_if_rflag(tMesh *mesh, int use_fv);
void refine_set_use_fv_if_pt_typ(tMesh *mesh, int pt_typ[3], int use_fv);
int refine_synchronize_ref_method(tRef *ref);
void hrefine_nodes_if_nb_finer_by_dl(tMesh *mesh, int dl, tRef *ref);
void hrefine_nodes_if_nb_finer(tMesh *mesh, tRef *ref);
void prefine_nodes_if_nb_uniform_in_any_dir(tMesh *mesh, tRef *ref);
void prefine_pat(tMesh *mesh, int p, int n[3]);
#define hrefine_nodes_if_rflag(mesh, ref) hrefine_elms_if_rflag(mesh, ref)
#define prefine_nodes_if_rflag(mesh, ref) prefine_elms_if_rflag(mesh, ref)
//#define remove_nodes_if_rflag(mesh, ref) remove_elms_if_rflag(mesh, ref)


/* connections.c */
void amr_set_elm_pat(tMesh *mesh, tElm *elm);
void amr_set_elm0_bbox(tMesh* mesh, tElm0 *elm0);
void amr_set_elm_bbox(tElm *elm);
void eloc_from_eploc(tEloc eloc[1], const tEploc eploc[1]);
void eloc_to_eploc(const tEloc eloc[1], tEploc eploc[1]);
void amr_set_sibling_elm0(const tElm *elm, int sib_ijk, tElm *sib);
char *elm_location_str(tElm *elm, char *s, int slen);
char *elmname(tElm *elm, char *s, int slen);
void eloc_from_elmname(tEloc *eloc, char *name);
tElm *elm_eid_from_elmname(tMesh *mesh, char *name, ulong *eid);
tElm *elm_from_elmname(tMesh *mesh,  char *name);
int elmname_is(tElm *elm, const char *nname);
int elm_get_ijk(tElm *elm);
tElm *elm_from_eid(tMesh *mesh, ulong eid, ulong *elmindex, int *elmrank);
tElm *elm_eploc_from_eid(tMesh *mesh, ulong eid, tEploc *eploc);
int elm_is_on_patface(tElm *elm, int f);
tElm0 *amr_alloc_get_elm0array_of_rank(tMesh *mesh, int rk, ulong *nelm0s);
ulong amr_nelms_on_rank(tMesh *mesh, int rank);
ulong amr_1st_eid_on_rank(tMesh *mesh, int rank);
void amr_elmindex_and_elmrank_of_eid(tMesh* mesh, ulong eid,
                                     ulong *elmindex, int *elmrank);
int amr_update_elm_nbinfo_if_nnbinfo_negative(tMesh *mesh);
int amr_update_elm_nbinfo_if_nnbinfo_negative_ef(tMesh *mesh,
                                                 khash_t(u32) *nbranks,
                                                 khash_t(u32_gptr) *ef);
int amr_elm_nbinfo_set_nnbinfo_mesh(tMesh *mesh, int positive);
int amr_elm_nbinfo_to_elm_fnb(tMesh *mesh);
int amr_get_nbelm_elmheaders(tMesh *mesh);

/* for compatibility with old connect.c */
#define node_location_str(node, s, slen) elm_location_str(node, s, slen)
#define nodename(node, s, slen)          elmname(node, s, slen)
#define node_location(node)              elm_location(node)
#define node_from_nodename(mesh, name)   elm_from_elmname(mesh, name)
#define nodename_is(node, nname)         elmname_is(node, nname)
