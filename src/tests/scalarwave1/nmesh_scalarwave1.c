/* nmesh_scalarwave1.h */
/* (c) Wolfgang Tichy 9/2019 */

#include "nmesh.h"
#include "scalarwave1.h"


int nmesh_scalarwave1(tMesh *mesh)
{
  if(!Getv(Par("physics"), "scalarwave1")) return 0;

  printf("Adding scalarwave1\n");

  /* functions */
  AddFun(POST_PARAMETERS, scalarwave1_init_global_pars);
  AddFun(INITIALDATA, scalarwave1_init);
  AddFun(ANALYZE, scalarwave1_analyze);

  /* variables */
  AddEvoVar("scalarwave1_phi",   "", "scalar field phi");
  AddEvoVar("scalarwave1_pi",    "", "pi  = d_t phi");
  AddEvoVar("scalarwave1_c",    "i", "c_i = d_i phi");
  AddAuxVar("scalarwave1_f_pi", "I", "flux for pi:  f_{pi}^i = -g^{ij} c_j");
  AddAuxVar("scalarwave1_f_c", "iJ", "flux for c_i: f_{c_i}^j= -g_i^j pi");
  AddAuxVar("scalarwave1_divf_pi", "", "d_i f_{pi}^i");
  AddAuxVar("scalarwave1_divf_c", "i", "d_j f_{c_i}^j");

AddAuxVar("scalarwave1_phi_r",   "", "scalar field phi");

  /* parameters */
  AddPar("scalarwave1_profile", "sin", "initial profile [sin,square]");
  AddPar("scalarwave1_k", "0.8 -0.5 0.33166247903554", "wave vector k^i");
  AddPar("scalarwave1_numflux", "upwind", "numerical flux [upwind,LLF]");
  AddPar("scalarwave1_limiter", "none", "limiter [none,MRS,minmodB]");

  return 0;
}
