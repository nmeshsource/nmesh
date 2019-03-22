/* output2d.c */
/* Wolfgang Tichy, Feb. 2019 */

#include "nmesh.h"
#include "output.h"



/* 2d output */
void output2d_meshvar(tMesh *mesh, char *name, int It, double T)
{
  tNode *node;
  int gnuplot = Getv(Par("2dformat"), "gnuplot");
/*
  int vtk     = Getv(Par("2dformat"), "vtk");
  int text    = Getv(Par("2dformat"), "text");
  int binary  = Getv(Par("2dformat"), "binary");
  int flt     = Getv(Par("2dformat"), "float");
  int dbl     = Getv(Par("2dformat"), "double");
*/
  int vi = Ind(name);
  FILE *fXY, *fXZ, *fYZ;
  char XYfil[1000];
  char XZfil[1000];
  char YZfil[1000];

  TIMER_START;

  /* loop over all nodes */
  forlnodes(mesh, node)
  {
    if(node->dat)
    if(node->dat->v[vi])
    {
      int p = node->pat->p;
      char ns[100];
      int ijk[3];
      
      //TODO: use diffferent Xb0 for diff patches
      double X0[] = { Getd(Par("outputX0")),
                      Getd(Par("outputY0")),
                      Getd(Par("outputZ0")) };
      
      /* find indices of nearest, if all are negative, node does not have
         outputX0, outputY0, outputZ0 */
      if(Getv(Par("outputX0Y0Z0coord"), "Xb"))
        nearest_ijk_of_XbYbZb(node, ijk, X0);
      else
        nearest_ijk_of_XYZ(node, ijk, X0);

      /* find string that idetifies node */
      node_location_str(node, ns,100);

      /* write files */
      /* XY-plane:  Z = Z0 */
      if(ijk[2]>=0)
      {
        if(gnuplot)
        {
          snprintf(XYfil, 999, "%s/%s.%02dXY%s",
                   Gets(Par("outdir")),name, p, ns);
          fXY = fopen(XYfil, "a");
          if(!fXY) errorexits("failed opening %s", XYfil);
          write_plane_ascii(node, fXY, 2, ijk, VarA(node, vi), It,T);
          fclose(fXY);
        }
      }

      /* XZ-plane:  Y = Y0 */
      if(ijk[1]>=0)
      {
        if(gnuplot)
        {
          snprintf(XZfil, 999, "%s/%s.%02dXZ%s",
                   Gets(Par("outdir")),name, p, ns);
          fXZ = fopen(XZfil, "a");
          if(!fXZ) errorexits("failed opening %s", XZfil);
          write_plane_ascii(node, fXZ, 1, ijk, VarA(node, vi), It,T);
          fclose(fXZ);
        }
      }

      /* YZ-plane:  X = X0 */
      if(ijk[0]>=0)
      {
        if(gnuplot)
        {
          snprintf(YZfil, 999, "%s/%s.%02dYZ%s",
                   Gets(Par("outdir")),name, p, ns);
          fYZ = fopen(YZfil, "a");
          if(!fYZ) errorexits("failed opening %s", YZfil);
          write_plane_ascii(node, fYZ, 0, ijk, VarA(node, vi), It,T);
          fclose(fYZ);
        }
      }
    }

    /* sysnchronize, so that we write only one node at a time */
    nMPI_barrier();
  } endforlnodes;

  TIMER_STOP;
}


/* output on patch planes */
void outputPatchPlanes_meshvar(tMesh *mesh, char *name, int It, double T)
{
  tNode *node;
  int gnuplot = Getv(Par("2dformat"), "gnuplot");
/*
  int vtk     = Getv(Par("2dformat"), "vtk");
  int text    = Getv(Par("2dformat"), "text");
  int binary  = Getv(Par("2dformat"), "binary");
  int flt     = Getv(Par("2dformat"), "float");
  int dbl     = Getv(Par("2dformat"), "double");
*/
  int vi = Ind(name);
  FILE *fpl;
  char plfil[1000];

  TIMER_START;

  /* loop over all nodes */
  forlnodes(mesh, node)
  {
    if(node->dat)
    if(node->dat->v[vi])
    {
      char ns[100];
      int ijk[3];
      int f;

      /* find string that idetifies node */
      nodename(node, ns,100);

      /* write files */
      /* pl-plane:  Z = Z0 */
      for(f=0; f<6; f++)
      {
        int norm = f/2;

        ijk[norm] = (node->n[norm]-1)*(f%2);

        if(gnuplot)
        {
          snprintf(plfil, 999, "%s/%s.%sf%d",
                   Gets(Par("outdir")),name, ns, f);
          fpl = fopen(plfil, "a");
          if(!fpl) errorexits("failed opening %s", plfil);
          write_plane_ascii(node, fpl, norm, ijk, VarA(node, vi), It,T);
          fclose(fpl);
        }
      }
    }

    /* sysnchronize, so that we write only one node at a time */
    nMPI_barrier();
  } endforlnodes;

  TIMER_STOP;
}
