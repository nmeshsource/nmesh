/* nmesh_burgers1.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "burgers1.h"


int nmesh_burgers1(tMesh *mesh)
{
  if(!Getv(Par("physics"), "burgers1")) return 0;

  printf("Adding burgers1\n");

  /* functions */
  AddFun(INITIALDATA, burgers1_init);
  AddFun(ANALYZE, burgers1_analyze);

  /* variables */
  AddEvoVar("burgers1_u", "",     "field we advect");
  AddAuxVar("burgers1_f", "I",    "f^i = n^i u");
  AddAuxVar("burgers1_f", "Ij",   "d_j f^i");
  AddAuxVar("burgers1_u_err", "", "error in u");

  /* parameters */
  AddPar("burgers1_direction", "1 0 0", "propagation direction n^i");

  return 0;
}
