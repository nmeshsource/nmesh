/* nmesh_SurfExchange.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "SurfExchange.h"


int nmesh_SurfExchange(tMesh *mesh)
{
  printf("Adding SurfExchange\n");

  /* functions */
  AddFun(INITIALDATA, SurfExchange_test);

  /* variables */
  AddEvoVar("SurfExchange_u", "", "that needs surface exchange");

  /* parameters */
  //AddPar("SurfExchange_lowlatency", "no", 
  //       "send many small rather than few large messages");

  return 0;
}
