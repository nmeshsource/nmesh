/* thread_defs.h */
/* Wolfgang Tichy, 5/2019 */


/* NOTE: In C99 these two have the same effect:
   #pragma omp parallel for
   _Pragma ( "omp parallel for" )
*/
/* to use OpenMP taks we have these macros */
#ifdef OMPTASKS
#include <omp.h>
#define tMUTEX omp_lock_t
#define DECL_MESH_MUTEX(node, mutex) \
  tMesh *mesh = (node)->pat->mesh; \
  tMUTEX *mutex = mesh->mutex;
#define MUTEX_INIT(x)    omp_init_lock(x)
#define MUTEX_LOCK(x)    omp_set_lock(x)
#define MUTEX_UNLOCK(x)  omp_unset_lock(x)
#define MUTEX_DESTROY(x) omp_destroy_lock(x)
#define T_CRITICAL(x )   _Pragma ( MSTR(omp critical x) )
#define TASK_YIELD       _Pragma ( "omp taskyield" )
#else
#define tMUTEX int
#define DECL_MESH_MUTEX(node, mutex)
#define MUTEX_INIT(x)
#define MUTEX_LOCK(x)
#define MUTEX_UNLOCK(x)
#define MUTEX_DESTROY(x)
#define T_CRITICAL(x)
#define TASK_YIELD
#endif

/* To parallelize with OpenMP we need _Pragma ( "omp parallel for" ) 
   in many places. But for different applications we want to switch 
   them on or off depending on where they are.
   SGRID_LEVEL2_Pragma used for omp loops over a plane in a pat (2d)
   SGRID_LEVEL3_Pragma used for omp loops over all points in a pat (3d)
   SGRID_LEVEL4_Pragma used for omp loops over all pates
   SGRID_LEVEL6_Pragma used for 6d omp loops (e.g. loop over pat while interpolating onto each point)
   more can be defined easily.
   */
/* define SGRID_LEVEL2_Pragma macros that allow us to include
   certain pragmas only if certain things like LEVEL2_Pragmas are defined */
#ifdef LEVEL2_Pragmas
#define SGRID_LEVEL2_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL2_Pragma(x)
#endif

#ifdef LEVEL3_Pragmas
#define SGRID_LEVEL3_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL3_Pragma(x)
#endif

#ifdef LEVEL4_Pragmas
#define SGRID_LEVEL4_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL4_Pragma(x)
#endif

#ifdef LEVEL6_Pragmas
#define SGRID_LEVEL6_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL6_Pragma(x)
#endif

#ifdef TOPLEVEL_Pragmas
#define SGRID_TOPLEVEL_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_TOPLEVEL_Pragma(x)
#endif

#if defined(LEVEL6_Pragmas) || defined(TOPLEVEL_Pragmas)
#define SGRID_LEVEL6orTOP_Pragma(x)  _Pragma ( #x )
#else
#define SGRID_LEVEL6orTOP_Pragma(x)
#endif
