/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */

/* to evolve e.g. gravity + matter we need at least 2 evolution subsystems
   that are coupled via source terms in their RHSs */
#define NUTEMP 6
typedef struct tEVOSYS {
  pVLList *u;            /* list of VarLists with evo vars of entire system */
  pVLList *w;            /* work list, needs to be an EvoVar just like u */
  pVLList *rhs;          /* RHS in du/dt = func(u, t) */
  pVLList *u_p;          /* u at previous time */
  pVLList *uk[NUTEMP];   /* temp. storage for say RK stages */
  void (*setrhs)(tNode *node, pVLList *rhs, pVLList *u, double time);
} tEvoSys;
