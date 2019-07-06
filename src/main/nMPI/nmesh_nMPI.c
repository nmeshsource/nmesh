/* nmesh_nMPI.h */
/* (c) Wolfgang Tichy 1/2019 */

#include "nmesh.h"
#include "nMPI.h"


int nmesh_nMPI(tMesh *mesh)
{
  printf("Adding nMPI\n");

  /* functions */
  AddFun(FIRST, nMPIvars_init);
  AddFun(POST_PARAMETERS, nMPI_print_compile_info);
  AddFun(FINALIZEMESH, nMPIvars_finalize);

  /* variables */
  //AddVar("nMPI_temp1", "", "temporary variable");

  /* parameters */
  //AddPar("nMPI_lowlatency", "no", 
  //       "send many small rather than few large messages");

  return 0;
}
