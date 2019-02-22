/* nmesh_misc.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "misc.h"


int nmesh_misc(tMesh *mesh)
{
  printf("Adding misc\n");

  /* functions */
  AddFun(INITIALDATA, misc_test);

  /* variables */
  AddEvoVar("misc_u", "", "some var");
  AddEvoVar("misc_v", "", "another var");
  //AddVarDim("misc_us", "", "that needs surface exchange",7,0,6);

  /* parameters */
  //AddPar("misc_lowlatency", "no", 
  //       "send many small rather than few large messages");

  return 0;
}
