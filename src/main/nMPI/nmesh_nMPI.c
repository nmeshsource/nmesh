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
  AddPar("nMPI_communicator_bits", "9", "the number of MPI communicators is "
         "2^nMPI_communicator_bits");

  return 0;
}
