/* nmesh_GhostExchange.h */
/* (c) Wolfgang Tichy 12/2019 */

#include "nmesh.h"
#include "GhostExchange.h"


int nmesh_GhostExchange(tMesh *mesh)
{
  if(!Getv(Par("physics"), "GhostExchange")) return 0;

  printf("Adding GhostExchange\n");

  /* functions */
  AddFun(INITIALDATA, GhostExchange_test);

  /* variables */
  AddEvoVar("GhostExchange_u", "", "var that needs ghost exchange");
  AddEvoVar("GhostExchange_v", "", "var that needs ghost exchange");

  /* parameters */
  //AddPar("GhostExchange_lowlatency", "no", "hmmmm");

  return 0;
}
