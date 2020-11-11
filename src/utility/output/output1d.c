/* output1d.c */
/* Wolfgang Tichy, Feb. 2019 */

#include "nmesh.h"
#include "output.h"



/* 1d output */
void output1d_vl(tVarList *vl, int It, double T)
{
  if(vl)
  {
    tMesh *mesh = vl->mesh;
    /*
    int fmt = Par("1dformat");
    int gnuplot = Getv(fmt, "gnuplot");
    int xdmf    = Getv(fmt, "xdmf");
    int vtk     = Getv(fmt, "vtk");
    int text    = Getv(fmt, "text");
    int binary  = Getv(fmt, "binary");
    int flt     = Getv(fmt, "float");
    int dbl     = Getv(fmt, "double");
    */
    int gnuplot=1, xdmf=0;

    TIMER_START;

    if(gnuplot)
    {
      int vli;
      for(vli=0; vli<vl->n; vli++)
      {
        int vi = vl->index[vli];
        char *vname = VarName(vi);
        gnuplot_output1d_meshvar(mesh, vname, It, T);
      }
    }
    if(xdmf)
    {
      //....
    }

    TIMER_STOP;
  }
}
