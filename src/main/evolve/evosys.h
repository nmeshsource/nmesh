/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */


/* lists of variable lists are defined in main/main/nmesh_main.h
   Here we just use an incomplete struct */
struct pVLLIST; /* pVLLIST is defined in main/main/nmesh_main.h */
struct tNODE;   /* nodes are defined in main/amr/nmesh_amr.h */

/* to evolve e.g. gravity + matter we need at least 2 evolution subsystems
   that are coupled via source terms in their RHSs */
#define NUTEMP 6
#define pVLL struct pVLLIST
typedef struct tEVOSYS {
  pVLL *u;            /* list of VarLists with evo vars of entire system */
  pVLL *w;            /* temp work list, needs to be an EvoVar just like u */
  pVLL *rhs;          /* RHS in du/dt = func(u, t) */
  pVLL *u_p;          /* u at previous time */
  pVLL *s[NUTEMP];    /* temp. storage for say RK stages */
  void (*setrhs)(struct tNODE *node, pVLL *rhs, pVLL *u);
} tEvoSys;
#undef pVLL
#undef NUTEMP
