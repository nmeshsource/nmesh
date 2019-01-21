/* nmesh_main.c */
/* Wolfgang Tichy, April 2005 */

#include "nmesh.h"
#include "main.h"



int nmesh_main(tMesh *mesh)
{
  printf("Adding main\n");

  /* functions */

  /* variables */

  /* parameters */
  AddPar("physics", "", "what problem to solve");
  AddPar("dt"  , "0.125", "time step dt");
  
  AddPar("iterations", "0", "number of mesh iterations");
  AddPar("finaltime", "0", "iterate until mesh reaches this time");
  AddPar("iterate_parameters", "no", "whether to iterate certain parameters");

  AddPar("errorexit", "exit", "how we exit in case of error [exit,abort]");
  return 0;
}
