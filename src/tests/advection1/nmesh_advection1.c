/* nmesh_advection1.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "advection1.h"


int nmesh_advection1(tMesh *mesh)
{
  printf("Adding advection1\n");

  /* functions */
  AddFun(INITIALDATA, advection1_init);
  AddFun(ANALYZE, advection1_analyze);

  /* variables */
  AddEvoVar("advection1_u", "",     "field we advect");
  AddAuxVar("advection1_f", "I",    "f^i = n^i u");
  AddAuxVar("advection1_f", "Ij",   "d_j f^i");
  AddAuxVar("advection1_F", "f",
            "normal comp. of numerical flux F = (f^{i*} - f^i) n^i");
  AddAuxVar("advection1_u_err", "", "error in u");

  /* parameters */
  AddPar("advection1_direction", "1 0 0", "propagation direction n^i");

  return 0;
}
