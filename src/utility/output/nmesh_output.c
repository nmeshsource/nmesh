/* nmesh_output.c */
/* Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "output.h"



int nmesh_output(tMesh *mesh)
{
  printf("Adding output\n");

  /* functions */
  AddFun(OUTPUT, mesh_output);

  /* variables */
  //AddAuxVar("output_temp1", "", "temporary variable(for vol. integrals)");

  /* parameters */
  AddPar("0doutiter", "-1", "when to output based on iterations");
  AddPar("1doutiter", "-1", "when to output based on iterations");
  AddPar("2doutiter", "-1", "when to output based on iterations");
  AddPar("3doutiter", "-1", "when to output based on iterations");

  AddPar("0douttime", "-1", "when to output based on time");
  AddPar("1douttime", "-1", "when to output based on time");
  AddPar("2douttime", "-1", "when to output based on time");
  AddPar("3douttime", "-1", "when to output based on time");

  AddPar("0doutput_normtype", "integral", 
	 "how we compute norms such as rms [integral,L2norm]");

  AddPar("0doutput", "", "variables to output norms for");
  AddPar("1doutput", "", "variables to output along axes");
  AddPar("2doutput", "", "variables to output on coordinate planes");
  AddPar("3doutput", "", "variables to output in 3d");

  AddPar("2dformat", "gnuplot",
	 "format for 2d output [gnuplot,vtk,text,binary,float,double]"); 
  AddPar("3dformat", "vtk binary float", "format for 3d output "
	 "[vtk,text,binary,float,double]"); 

  AddPar("outputX0", "0", "origin for output in X");
  AddPar("outputY0", "0", "origin for output in Y");
  AddPar("outputZ0", "0", "origin for output in Z");
  AddPar("outputX0Y0Z0coord", "X", "meaning of X in outputX0 [X,Xb]");

  return 0;
}
