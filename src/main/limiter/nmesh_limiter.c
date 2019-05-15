/* nmesh_limiter.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "limiter.h"


int nmesh_limiter(tMesh *mesh)
{
  printf("Adding limiter\n");

  /* functions */
  AddFun(POST_PARAMETERS, limiter_init_global_par_indices);

  /* variables */
  //AddAuxVar("limiter_u",      "",    "test function");
   
  /* parameters */
  AddPar("limiter_name", "MRS", "name of limiter [MRS,minmodB]");
  AddPar("limiter_alpha", "5", "for MRS and also minmodB: [0,inf], "
         "0 is most aggressive");
  AddPar("limiter_beta", "0.5", "for minmod: [0.5,1], 0.5 is TVD");
	     	   	   	 
  return 0;
}
