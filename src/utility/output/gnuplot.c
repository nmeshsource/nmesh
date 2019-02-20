/* gnuplot.c */
/* Wolfgang Tichy, Feb. 2019 */

#include "nmesh.h"
#include "output.h"

/* the functions in this file write data in a text format that is understood
   by gnuplot and tgraph */


/* write a single plane */
void write_plane_ascii(tNode *node, FILE *fp, int normal, int plane[], int iv,
                       int Iter, double Time)
{
  tArray *p1;
  tArray *p2;
  double Xb[3], X[3];
  double *pv = GetVarDpointer(node, iv);
  tArray *va = GetVarArray(node, iv);
  int i,j,k;
  int imin, jmin, kmin;
  int imax, jmax, kmax;

  if(pv==NULL) return;

  fprintf(fp, "# \"time = %.15g\"", Time);
  
  imin = jmin = kmin=0;
  imax = va->n[0] - 1;
  jmax = va->n[1] - 1;
  kmax = va->n[2] - 1;

  if(normal==0)
  {
    p1 = node->Xb[1];
    p2 = node->Xb[2];
    if(plane[0]>imax) imin = imax;
    else              imin = imax = plane[0];
    XYZ_of_ijk(node, plane[0],0,0, X);
    fprintf(fp, ", i=%d, X=%.15g\n", plane[0], X[0]);
  }
  else if(normal==1)
  {
    p1 = node->Xb[0];
    p2 = node->Xb[2];
    if(plane[1]>jmax) jmin = jmax;
    else              jmin = jmax = plane[1];
    XYZ_of_ijk(node, 0,plane[1],0, X);
    fprintf(fp, ", j=%d, Y=%.15g\n", plane[1], X[1]);
  }
  else
  {
    p1 = node->Xb[0];
    p2 = node->Xb[1];
    if(plane[2]>kmax) kmin = kmax;
    else              kmin = kmax = plane[2];
    XYZ_of_ijk(node, 0,0,plane[2], X);
    fprintf(fp, ", k=%d, Z=%.15g\n", plane[2], X[2]);
  }

  /* go over plane, with normal */
  for(k=kmin; k<=kmax; k++)
  {
    for(j=jmin; j<=jmax; j++)
    {
      for(i=imin; i<=imax; i++)
      {
        int dir1 = Dir1_norm(normal);
        int dir2 = Dir2_norm(normal);
        int ind1 = Ind1_norm(i,j,k, normal);
        int ind2 = Ind2_norm(i,j,k, normal);
        int indv = Ind_n(i,j,k, va->n);

        Xb[dir1] = p1->d[ind1];
        Xb[dir2] = p2->d[ind2];
        XYZ_of_XbYbZb(node, Xb, X);
        //fprintf(fp, "%d %d %d: %d %d  %d %d: ", i,j,k, ind1,ind2,dir1,dir2);
        //fprintf(fp, "%g %g  %g %g: \t\t", Xb[dir1],Xb[dir2], X[dir1],X[dir2]);
        //fprintf(fp, "%d %d %d: %d %d: %d: ", i,j,k, ind1,ind2, indv);
        fprintf(fp, "%.15g %.15g %.15g\n", X[dir1], X[dir2], pv[indv]);
      }
      if(normal==2) fprintf(fp, "\n");
    }
    if(normal==1 || normal==0) fprintf(fp, "\n");
  }
  fprintf(fp,"\n");
}


/* write along a line in direc. dir */
void write_line_ascii(tNode *node, FILE *fp, int dir, int axis[], int iv,
                      int Iter, double Time)
{
  tArray *p1;
  double Xb[3], X[3];
  double *pv = GetVarDpointer(node, iv);
  tArray *va = GetVarArray(node, iv);
  int i,j,k;
  int imin, jmin, kmin;
  int imax, jmax, kmax;

  if(pv==NULL) return;

  fprintf(fp, "# \"time = %.15g\"", Time);
  
  imin = jmin = kmin=0;
  imax = va->n[0] - 1;
  jmax = va->n[1] - 1;
  kmax = va->n[2] - 1;
  
  switch(dir)
  {
  case 0:
    p1 = node->Xb[0];
    if(axis[1]>jmax) jmin = jmax;
    else             jmin = jmax = axis[1];
    if(axis[2]>kmax) kmin = kmax;
    else             kmin = kmax = axis[2];
    XYZ_of_ijk(node, 0,jmin,kmin, X);
    fprintf(fp, ", j=%d, k=%d, Y=%.15g, Z=%.15g\n", jmin, kmin, X[1], X[2]);
    break;
  case 1:
    p1 = node->Xb[1];
    if(axis[0]>imax) imin = imax;
    else             imin = imax = axis[0];
    if(axis[2]>kmax) kmin = kmax;
    else             kmin = kmax = axis[2];
    XYZ_of_ijk(node, imin,0,kmin, X);
    fprintf(fp, ", i=%d, k=%d, X=%.15g, Z=%.15g\n", imin, kmin, X[0], X[2]);
    break;
  case 2:
    p1 = node->Xb[2];
    if(axis[0]>imax) imin = imax;
    else             imin = imax = axis[0];
    if(axis[1]>jmax) jmin = jmax;
    else             jmin = jmax = axis[1];
    XYZ_of_ijk(node, imin,jmin,0, X);
    fprintf(fp, ", i=%d, j=%d, X=%.15g, Y=%.15g\n", imin, jmin, X[0], X[1]);
    break;
  default:
    errorexit("dir has to be 0,1,2");
  }

  /* go over line */
  for(k=kmin; k<=kmax; k++)
  {
    for(j=jmin; j<=jmax; j++)
    {
      for(i=imin; i<=imax; i++)
      {
        int ind1 = Ind_dir(i,j,k, dir);
        int indv = Ind_n(i,j,k, va->n);

        Xb[dir] = p1->d[ind1];
        XYZ_of_XbYbZb(node, Xb, X);
        //fprintf(fp, "%d %d %d: %d %d: ", i,j,k, ind1,dir);
        //fprintf(fp, "%g %g: ", Xb[dir], X[dir]);
        //fprintf(fp, "%d %d %d: %d: %d: ", i,j,k, ind1, indv);
        fprintf(fp, "%.15g %.15g \n", X[dir], pv[indv]);
      }
    }
  }
  fprintf(fp,"\n");
}
