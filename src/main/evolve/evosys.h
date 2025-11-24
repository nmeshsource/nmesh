/* evolve.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coupled evolution systems */



/* Function bins for tasks during one evolution sub step */
enum
{
  PRESURF0, /* set some source terms that are needed early, e.g. gmunu and
               its time- and space-derivs */
  PRESURF,  /* use PRESURF0 vars to set other vars needed early */
  PRESURF2, /* set vars needed before surf exchange, e.g. prims */

  SETSRC0,  /* set some source terms needed before SETSRC */
  SETSRC,   /* set some source terms, is called before VOLRHS */

  XVOLRHS,  /* set some extra vars needed for LDG, sets vol. terms */
  XSURFRHS, /* add some extra vars needed for LDG on surf. */

  VOLRHS,   /* set vol. terms of RHS of evo eqns (after setsrc) */
  SURFRHS,  /* add RHS terms from surf. fluxes (after volrhs) */

  PRELIM0,  /* set vars that are needed very early, e.g. gmunu */
  PRELIM,   /* use PRELIM0 vars to set other vars needed early */
  LIMDATA,  /* produce data such as min,max on each node */
            /* NOTE: ListEntry(evosys->f[LIMDATA],i)(NULL, evv)
                     must return number of data vals we need */
  LIMITER,  /* apply limiter on node using data from limdata */

  TROUBLE,  /* calculate trouble score */

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
#define pFL struct EvoFuncPtrLIST
#define pSL struct constStringLIST
typedef struct tEVOSYS {
  pVLL *u;            /* list of VarLists with evo vars of entire system */
  pVLL *rhs;          /* RHS in du/dt = func(u, t), rhs is AuxVar */
  pVLL *w;            /* temp work list, needs to be an EvoVar just like u */
  pVLL *u_p;          /* u at previous time */
  pVLL *s[NEVOTEMP];  /* temp. storage for say RK stages */
  pFL *f[NEVOFUNCBINS];      /* one list of func. pointers in each evo bin */
  pSL *f_name[NEVOFUNCBINS]; /* one list of func. names in each evo bin */
  pVLL *x; /* extra vars in RHS for LDG, x is AUXVAR with surfacezones */
} tEvoSys;
#undef pVLL
//#undef NEVOTEMP



/* a struct that contains all we need to call RHS funcs */
typedef struct tEVOVARS {
  struct tVARLIST *vlu;   /* varlist of all evolved vars for one system */
  struct tVARLIST *vlr;   /* varlistof RHS vars for one system */
  struct tVARLIST *vlx;   /* vl with extra vars (for LDG) for one sys. */
  struct tVARLIST *vlu_p; /* vlu at previous time */
  void *pars;             /* extra pars we may want to pass to a RHS */
} tEvoVars;

/* macros to access tEvoVars */
#define EvoVars_vlu(evv)    (evv)->vlu
#define EvoVars_vlr(evv)    (evv)->vlr
#define EvoVars_vlx(evv)    (evv)->vlx
#define EvoVars_vlu_p(evv)  (evv)->vlu_p
#define EvoVars_pars(evv)   (evv)->pars
#define EvoVarsSet_vlu(evv, vl)  (evv)->vlu = (vl)
#define EvoVarsSet_vlr(evv, vl)  (evv)->vlr = (vl)
