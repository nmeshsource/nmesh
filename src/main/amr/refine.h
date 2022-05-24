/* refine.h */
/* Wolfgang Tichy, 5/2022 */


/**************************************************************************/
/* for mesh refinement */
/**************************************************************************/

/* refinement types: i.e. vals for var type in tRef */
enum
{
  H_REFINE,      /* do h-refinement */
  P_REFINE       /* do p-refinement */
};

/* refinement methods: i.e. vals for var method in tRef */
enum
{
  REF_METH_DONOTHING=0,  /* do nothing */
  PARENT_n,              /* use same n as parent */
  PARENT_nO2,            /* use parent->n/2 */
  PARENT_nO2_P1,         /* use parent->n/2 + 1 */
  PARENT_nO2_P1IFnG3,    /* use parent->n/2 + 1*if(parent->n>3) */
  PARENT_nO2_P1MOD,      /* use parent->n/2 + 1 or parent->n - 1 if n<=3 */
  GIVEN_n,               /* use a n given by the user */
  PARENT_n_P_LGL,        /* use same n as parent, and LGL gridpoints */
  PARENT_n_P_UNIFORM,    /* use same n as parent, and UNIFORM gridpoints */
  PARENT_nO2_P_LGL,      /* use parent->n/2, and LGL gridpoints */
  PARENT_nO2_P_UNIFORM,  /* use parent->n/2, and UNIFORM gridpoints */
  PARENT_2n_P_LGL,       /* use parent->n*2, and LGL gridpoints */
  PARENT_2n_P_UNIFORM,   /* use parent->n*2, and UNIFORM gridpoints */
  GIVEN_n_P_LGL,         /* use given n, and LGL gridpoints */
  GIVEN_n_P_UNIFORM,     /* use given n, and UNIFORM gridpoints */
  NREF_METHODS,          /* number of valid refinement methods */
  REF_METH_INVALID //invalid refinement method, add nothing beyond this!!!
};

/* struct that holds refinement method */
typedef struct tREF {
  int type;       /* H_REFINE, P_REFINE */
  int method;     /* NOREFINE, PARENT_n, ... */
  int n[3];       /* n to use if method=GIVEN_n */
} tRef;
