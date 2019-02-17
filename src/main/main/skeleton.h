/* skeleton.h */
/* Wolfgang Tichy, 1/2019 */

/* for skeleton.c */
enum
{
  FIRST,
  POST_PARAMETERS,
  INITMESH,
  PRE_COORDINATES,
  COORDINATES,
  PRE_INITIALDATA,
  INITIALDATA,
  POST_INITIALDATA,
  PRE_EVOLVE,
  EVOLVE,
  POST_EVOLVE,
  ANALYZE, 
  OUTPUT,
  POST_OUTPUT,
  AMR,
  LOADBALANCING,
  NFUNCBINS
};

struct tMESH; /* incompl. struct dec. so we can use struct tMESH* already */

typedef struct tTODO {
  struct tTODO *next;
  int (*f)(struct tMESH *);
  char *name;
} tTodo;
