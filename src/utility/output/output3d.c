/* output3d.c */
/* Wolfgang Tichy, Feb. 2019 */

#include "nmesh.h"
#include "output.h"




/* 3d output */
void output3d_vl(tVarList *vl, int It, double T)
{
  tMesh *mesh = vl->mesh;
  int fmt = Par("3dformat");
  int xdmf = Getv(fmt, "xdmf");
  int vtk  = Getv(fmt, "vtk");
/*
  int gnuplot = Getv(fmt, "gnuplot");
  int text    = Getv(fmt, "text");
  int binary  = Getv(fmt, "binary");
  int flt     = Getv(fmt, "float");
  int dbl     = Getv(fmt, "double");
*/

  TIMER_START;

  if(vtk)
  {
    int vli;
    for(vli=0; vli<vl->n; vli++)
    {
      int vi = vl->index[vli];
      char *vname = VarName(vi);
      vtk_output3d_meshvar(mesh, vname, It, T);
    }
  }
  if(xdmf)
  {
    output3d_xdmf(vl, It, T);
  }

  TIMER_STOP;
}
