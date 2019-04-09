/* nmesh_limiter.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "limiter.h"


int nmesh_limiter(tMesh *mesh)
{
  printf("Adding limiter\n");

  /* functions */
  //AddFun(INITIALDATA, limiter_startup);
  //AddFun(ANALYZE, limiter_analyze);

  /* variables */
  //AddAuxVar("limiter_u",      "",    "test function");
   
  /* parameters */
  AddPar("limiter_name", "MRS", "name of limiter");
	     	   	   	 
  return 0;
}
