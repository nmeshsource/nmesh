/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */


/* lists of variable lists are defined in main/main/nmesh_main.h
   Here we just use an incomplete struct */
struct pVLLIST; /* pVLLIST is defined in main/main/nmesh_main.h */
struct tNODE;   /* nodes are defined in main/amr/nmesh_amr.h */

/* to evolve e.g. gravity + matter we need at least 2 evolution subsystems
   that are coupled via source terms in their RHSs */
#define NEVOTEMP 6
#define pVLL struct pVLLIST
#define pFL struct FuncPointerLIST
typedef struct tEVOSYS {
  pVLL *u;            /* list of VarLists with evo vars of entire system */
  pVLL *rhs;          /* RHS in du/dt = func(u, t), rhs is AuxVar */
  pVLL *w;            /* temp work list, needs to be an EvoVar just like u */
  pVLL *u_p;          /* u at previous time */
  pVLL *s[NEVOTEMP];  /* temp. storage for say RK stages */
  /* func. pointers that are called in this order: */
  pFL *prelim;        /* set vars that are need early, e.g. gmunu */
  pFL *limdata;       /* produce data such as min,max on each node */
                      /* NOTE: ListEntry(evosys->limdata,i)(NULL, vl)
                               must return number of data vals we need */
  pFL *limiter;       /* apply limiter on node using data from limdata */
  pFL *presurf;       /* set vars needed before surf exchange, e.g. prims */
  pFL *setsrc;        /* set some source terms, is called before setrhs */
  pFL *setrhs;        /* set RHS of eve eqns, called after setsrc */
} tEvoSys;
#undef pVLL
//#undef NEVOTEMP
