/* nmesh_dg.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "dg.h"


int nmesh_dg(tMesh *mesh)
{
  printf("Adding dg\n");

  /* functions */
  //AddFun(INITIALDATA, dg_startup);
  //AddFun(ANALYZE, dg_analyze);

  /* variables */
  //AddAuxVar("dg_u",      "",    "test function");
   
  /* parameters */
  AddPar("dg_numerical_flux", "LLF", "numerical flux [LLF]");
	     	   	   	 
  return 0;
}
