/* thread_defs.h */
/* Wolfgang Tichy, 5/2019 */


/* NOTE: In C99 these two have the same effect:
   #pragma omp parallel for
   _Pragma ( "omp parallel for" )
*/
/* to use OpenMP taks we have these macros */
#ifdef USEOMP

#include <omp.h>
#define FORNODES_Pragmas
#define tMUTEX omp_lock_t
#define DECL_MESH_MUTEX(node, mutex) \
  tMesh *mesh = (node)->pat->mesh; \
  tMUTEX *mutex = mesh->mutex;
#define MUTEX_INIT(x)    omp_init_lock(x)
#define MUTEX_LOCK(x)    omp_set_lock(x)
#define MUTEX_UNLOCK(x)  omp_unset_lock(x)
#define MUTEX_DESTROY(x) omp_destroy_lock(x)
#define T_CRITICAL       _Pragma ( MSTR_OFVAL(omp critical) )
#define TASK_YIELD       _Pragma ( "omp taskyield" )
#define MAX_NTHREADS     omp_get_max_threads()

#else

#define tMUTEX int
#define DECL_MESH_MUTEX(node, mutex)
#define MUTEX_INIT(x)
#define MUTEX_LOCK(x)
#define MUTEX_UNLOCK(x)
#define MUTEX_DESTROY(x)
#define T_CRITICAL(x)
#define TASK_YIELD
#define MAX_NTHREADS 1

#endif

/* To parallelize with OpenMP we need _Pragma ( "omp parallel for" ) 
   in many places. But for different applications we want to switch 
   them on or off depending on where they are.
   FORNODES_Pragma used for omp loops all nodes
*/
/* define FORNODES_Pragma macros that allow us to include
   certain pragmas only if certain things like FORNODES_Pragmas are defined */
#ifdef FORNODES_Pragmas
#define FORNODES_Pragma(x)  _Pragma ( #x )
#else
#define FORNODES_Pragma(x)
#endif

#if defined(LEVEL6_Pragmas) || defined(TOPLEVEL_Pragmas)
#define SGRID_LEVEL6orTOP_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL6orTOP_Pragma(x)
#endif
