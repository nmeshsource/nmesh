/* nmesh_advection1.h */
/* (c) Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "advection1.h"


int nmesh_advection1(tMesh *mesh)
{
  if(!Getv(Par("physics"), "advection1")) return 0;

  printf("Adding advection1\n");

  /* functions */
  AddFun(POST_PARAMETERS, advection1_init_global_pars);
  AddFun(INITIALDATA, advection1_init);
  AddFun(ANALYZE, advection1_analyze);
  //FIXME: test hack:
  //AddFun(POST_EVOLVE, resolve_shocks_using_nlim);

  /* variables */
  AddEvoVar("advection1_u", "",     "field we advect");
  AddAuxVar("advection1_f", "I",    "f^i = n^i u");
  //AddAuxVar("advection1_f", "Ij",   "d_j f^i");
  AddAuxVar("advection1_divf", "",  "d_i f^i");
  AddAuxVar("advection1_u_err", "", "error in u");

  /* parameters */
  AddPar("advection1_profile", "sin", "initial profile [sin,square]");
  AddPar("advection1_direction", "1 0 0", "propagation direction n^i");
  AddPar("advection1_numflux", "upwind", "numerical flux [LLF,upwind]");
  AddPar("advection1_limiter", "none", "limiter [none,MRS,minmodB]");
  AddPar("advection1_outerBC_influxes", "no", " [no,yes]");

  AddPar("advection1_refine", "no", "use dynamic refinement [no,yes]");
  if(Getb(Par("advection1_refine")))
  {
    AddFun(PRE_EVOLVE, advection1_refine);
  }

  return 0;
}
