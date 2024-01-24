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
  AddPar("cooutiter", "-1", "when to output based on iterations");

  AddPar("0douttime", "-1", "when to output based on time");
  AddPar("1douttime", "-1", "when to output based on time");
  AddPar("2douttime", "-1", "when to output based on time");
  AddPar("3douttime", "-1", "when to output based on time");
  AddPar("coouttime", "-1", "when to output based on time");

  AddPar("0doutput_normtype", "integral",
	 "how we compute norms such as rms [integral,L2norm]");
  AddPar("0doutput_per_patch", "128", "output per patch [yes,no,#], here"
         "the # means yes, but only if patchnumber does not exceed #");
  AddPar("0doutput_add_xyz", "yes", "add location x,y,z of min/max to"
         "datafiles after min/max values [yes,no]");

  AddPar("0doutput", "", "variables to output norms for");
  AddPar("1doutput", "", "variables to output along axes");
  AddPar("2doutput", "", "variables to output on coordinate planes");
  AddPar("3doutput", "", "variables to output in 3d");
  AddPar("cooutput", "", "variables coeffs to output");

  AddPar("0doutputall", "yes", "whether to output all components");
  AddPar("1doutputall", "yes", "whether to output all components");
  AddPar("2doutputall", "yes", "whether to output all components");
  AddPar("3doutputall", "yes", "whether to output all components");
  AddPar("cooutputall", "yes", "whether to output all components");

  AddPar("1dformat", "gnuplot pernode",
	 "format for 1d output [gnuplot;pernode]");
  AddPar("2dformat", "gnuplot",
	 "format for 2d output [xdmf,gnuplot,vtk,text,binary,float,double]");
  AddPar("3dformat", "vtk binary float", "format for 3d output "
	 "[xdmf,vtk,text,binary,float,double]");
  AddPar("coformat", "vtk text arrange_as_1d", "format for coef output "
	 "[vtk,text,binary,float,double]");

  /* pars to limit 2d or 3d output to certain regions */
  AddPar("2doutputregion", "all", "do 2d output only inside certain regions "
         "[all,sphere0,sphere1,sphere2]");
  AddPar("3doutputregion", "all", "do 3d output only inside certain regions "
         "[all,sphere0,sphere1,sphere2]");
  /* sphere0,sphere1,sphere2 are centered on pt0,pt1,pt2 */
  AddPar("output_sphere0_radius", "40", "radius of sphere0  3doutputregion");
  AddPar("output_sphere1_radius", "10", "radius of sphere1 in 3doutputregion");
  AddPar("output_sphere2_radius", "10", "radius of sphere2 in 3doutputregion");

  /* pars for 1d output */
  AddPar("outputX0", "0", "origin for output in X");
  AddPar("outputY0", "0", "origin for output in Y");
  AddPar("outputZ0", "0", "origin for output in Z");
  AddPar("outputX0Y0Z0coord", "X", "meaning of X in outputX0 [X,Xb]");

  AddPar("output_xcoord", "x", "var used as Cartesian x-coordinate");

  return 0;
}
