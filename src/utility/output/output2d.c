/* output2d.c */
/* Wolfgang Tichy, Feb. 2019 */

#include "nmesh.h"
#include "output.h"



/* 2d output */
void output2d_vl(tVarList *vl, int It, double T)
{
  tMesh *mesh = vl->mesh;
  int fmt = Par("2dformat");
  int gnuplot = Getv(fmt, "gnuplot");
  int xdmf    = Getv(fmt, "xdmf");
/*
  int vtk     = Getv(fmt, "vtk");
  int text    = Getv(fmt, "text");
  int binary  = Getv(fmt, "binary");
  int flt     = Getv(fmt, "float");
  int dbl     = Getv(fmt, "double");
*/

  TIMER_START;

  if(gnuplot)
  {
    int vli;
    for(vli=0; vli<vl->n; vli++)
    {
      int vi = vl->index[vli];
      char *vname = VarName(vi);
      gnuplot_output2d_meshvar(mesh, vname, It, T);
    }
  }
  if(xdmf)
  {
    output2d_xdmf(vl, It, T);
  }

  TIMER_STOP;
}
