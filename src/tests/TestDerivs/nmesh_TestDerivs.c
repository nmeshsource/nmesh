/* nmesh_TestDerivs.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "TestDerivs.h"


int nmesh_TestDerivs(tMesh *mesh)
{
  if (!Getv(Par("physics"), "TestDerivs")) return 0;
  printf("Adding TestDerivs\n");

  /* functions */
  AddFun(INITIALDATA, TestDerivs_startup);
  AddFun(ANALYZE, TestDerivs_analyze);

  /* variables */
  AddAuxVar("TestDerivs_u",      "",   "test function");
  AddAuxVar("TestDerivs_Err_du", "i",  "Error in 1st deriv of test function");
  AddAuxVar("TestDerivs_Err_ddu","ij", "Error in 2nd deriv of test function");
   
  /* parameters */
  AddPar("TestDerivs_A",     "1",   "amplitude of wave");
  AddPar("TestDerivs_sigmax","1",   "sigmax of Gaussian wavepacket");
  AddPar("TestDerivs_sigmay","1.3", "sigmay of Gaussian wavepacket");
  AddPar("TestDerivs_sigmaz","0.7", "sigmaz of Gaussian wavepacket");
  AddPar("TestDerivs_x0",    "1",   "x-location of Gaussian wavepacket");
  AddPar("TestDerivs_y0",    "0.3", "y-location of Gaussian wavepacket");
  AddPar("TestDerivs_z0",    "0.1", "z-location of Gaussian wavepacket");
	     	   	   	 
  return 0;
}
