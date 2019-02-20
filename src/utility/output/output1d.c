/* output1d.c */
/* Wolfgang Tichy, Feb. 2019 */

#include "nmesh.h"
#include "output.h"



/* 1d output */
void output1d_meshvar(tMesh *mesh, char *name, int It, double T)
{
  tNode *node;
  int gnuplot = 1; /* we always use the same format in 1d */
  int vi = Ind(name);
  FILE *fX, *fY, *fZ;
  char Xfil[1000];
  char Yfil[1000];
  char Zfil[1000];

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
      /* X-axis:  Y = Y0, Z = Z0 */
      if(ijk[1]>=0 && ijk[2]>=0)
      {
        if(gnuplot)
        {
          snprintf(Xfil, 999, "%s/%s_%02d.X%s",
                   Gets(Par("outdir")),name, p, ns);
          fX = fopen(Xfil, "a");
          if(!fX) errorexits("failed opening %s", Xfil);
          write_line_ascii(node, fX, 0, ijk, vi, It,T);
          fclose(fX);
        }
      }

      /* Y-axis:  X = X0, Z = Z0 */
      if(ijk[0]>=0 && ijk[2]>=0)
      {
        if(gnuplot)
        {
          snprintf(Yfil, 999, "%s/%s_%02d.Y%s",
                   Gets(Par("outdir")),name, p, ns);
          fY = fopen(Yfil, "a");
          if(!fY) errorexits("failed opening %s", Yfil);
          write_line_ascii(node, fY, 1, ijk, vi, It,T);
          fclose(fY);
        }
      }

      /* Z-axis:  X = X0, Y = Y0 */
      if(ijk[0]>=0 && ijk[1]>=0)
      {
        if(gnuplot)
        {
          snprintf(Zfil, 999, "%s/%s_%02d.Z%s",
                   Gets(Par("outdir")),name, p, ns);
          fZ = fopen(Zfil, "a");
          if(!fZ) errorexits("failed opening %s", Zfil);
          write_line_ascii(node, fZ, 2, ijk, vi, It,T);
          fclose(fZ);
        }
      }
    }

    /* sysnchronize, so that we write only one node at a time */
    nMPI_barrier();
  } endforlnodes;
}
