/* skeleton.h */
/* Wolfgang Tichy, 1/2019 */

/* for skeleton.c */
enum
{
  FIRST,            /* in main: just after libraries have been initialized */

  //timebins called once in inidata_mesh:
  POST_PARAMETERS,
  INITMESH,         /* happens only if there is no checkpoint restart */
  POST_INITMESH,    /* happens only if there is no checkpoint restart */
  PRE_COORDINATES,
  COORDINATES,
  POST_COORDINATES,
  PRE_INITIALDATA,
  INITIALDATA,
  POST_INITIALDATA,
  POST_CHECKPOINT_LOAD, /* after checkpoint_load_stage(mesh, "", 1); */

  //timebins called repeatedly in evolve_mesh while mesh->time < timemax:
  PRE_EVOLVE,
  EVOLVE,
  FILTER,
  POST_EVOLVE,
  ANALYZE, 
  OUTPUT,
  POST_OUTPUT,
  LOADBALANCING,    /* before AMR, since refinements invalidate timings */
  AMR,

  //timebins called before nmesh finishes
  FINALIZE,
  FINALIZEMESH,
  POST_FINALIZEMESH,

  NFUNCBINS
};


struct tMESH; /* incompl. struct dec. so we can use struct tMESH* already */

typedef struct tTODO {
  struct tTODO *next;
  int (*f)(struct tMESH *);
  char *name;
} tTodo;
