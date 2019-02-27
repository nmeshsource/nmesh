/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */

/* to evolve e.g. gravity + matter we need at least 2 evolution subsystems
   that are coupled via source terms in their RHSs */
#define NUTEMP 6
#define pVLL struct tVARLIST *
typedef struct tEVOSYS {
  pVLL *u;            /* list of VarLists with evo vars of entire system */
  pVLL *w;            /* work list, needs to be an EvoVar just like u */
  pVLL *rhs;          /* RHS in du/dt = func(u, t) */
  pVLL *u_p;          /* u at previous time */
  pVLL *uk[NUTEMP];   /* temp. storage for say RK stages */
  void (*setrhs)(pVLL *rhs, pVLL *u, double time);
} tEvoSys;
#undef pVLL
#undef NUTEMP
