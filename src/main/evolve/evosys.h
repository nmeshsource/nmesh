/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */



/* Function bins for tasks during one evolution sub step */
enum
{
  PRESURF, //presurf; /* set vars needed before surf exchange, e.g. prims */
  SETSRC,  //setsrc;  /* set some source terms, is called before volrhs */
  VOLRHS,  //volrhs;  /* set vol. terms of RHS of evo eqns (after setsrc) */
  SURFRHS, //surfrhs; /* add RHS terms from surf. fluxes (after volrhs) */


  PRELIM,  //prelim; /* set vars that are needed early, e.g. gmunu */
  LIMDATA, //limdata; /* produce data such as min,max on each node */
                 /* NOTE: ListEntry(evosys->limdata,i)(NULL, vl)
                          must return number of data vals we need */
  LIMITER,  //limiter; /* apply limiter on node using data from limdata */

  NEVOFUNCBINS /* number of function bins in this enum */
};



/* lists of variable lists are defined in main/main/nmesh_main.h
   Here we just use an incomplete struct */
struct pVLLIST; /* pVLLIST is defined in main/main/nmesh_main.h */
struct tNODE;   /* nodes are defined in main/amr/nmesh_amr.h */

/* to evolve e.g. gravity + matter we need at least 2 evolution subsystems
   that are coupled via source terms in their RHSs */
#define NEVOTEMP 6
#define pVLL struct pVLLIST
#define pFL struct FuncPointerLIST
#define pSL struct StringLIST
typedef struct tEVOSYS {
  pVLL *u;            /* list of VarLists with evo vars of entire system */
  pVLL *rhs;          /* RHS in du/dt = func(u, t), rhs is AuxVar */
  pVLL *w;            /* temp work list, needs to be an EvoVar just like u */
  pVLL *u_p;          /* u at previous time */
  pVLL *s[NEVOTEMP];  /* temp. storage for say RK stages */
  pFL *f[NEVOFUNCBINS];      /* one list of func. pointers in each evo bin */
  pSL *f_name[NEVOFUNCBINS]; /* one list of func. names in each evo bin */
} tEvoSys;
#undef pVLL
//#undef NEVOTEMP
