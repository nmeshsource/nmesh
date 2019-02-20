/* gnuplot2d.c */
/* Wolfgang Tichy, Feb. 2019 */

#include "nmesh.h"
#include "output.h"



void write_plane_ascii(tNode *node, FILE *fp, int normal, int plane[], int iv,
                       int Iter, double Time)
{
  double *p1;
  double *p2;
  double Xb[3], X[3];
  double *pv = GetVarDpointer(node, iv);
  tArray *va = GetVarArray(node, iv);
  int i,j,k;
  int imin, jmin, kmin;
  int imax, jmax, kmax;
  int index;

  if(pv==NULL) return;

  fprintf(fp, "# \"time = %.16g\"", Time);
  
  imin = jmin = kmin=0;
  imax = node->n[0] - 1;
  jmax = node->n[1] - 1;
  kmax = node->n[2] - 1;
  
  if(normal==0)
  {
    p1 = node->Xb[1]->d;
    p2 = node->Xb[2]->d;
    imin = imax = plane[0];
    XbYbZb_of_ijk(node, plane[0],0,0, Xb);
    XYZ_of_XbYbZb(node, Xb, X);
    fprintf(fp, ", i=%d, X=%.16g\n", plane[0], X[0]);
  }
  else if(normal==1)
  {
    p1 = node->Xb[0]->d;
    p2 = node->Xb[2]->d;
    jmin = jmax = plane[1];
    XbYbZb_of_ijk(node, 0,plane[1],0, Xb);
    XYZ_of_XbYbZb(node, Xb, X);
    fprintf(fp, ", j=%d, Y=%.16g\n", plane[1], X[1]);
  }
  else
  {
    p1 = node->Xb[0]->d;
    p2 = node->Xb[1]->d;
    kmin = kmax = plane[2];
    XbYbZb_of_ijk(node, 0,0,plane[2], Xb);
    XYZ_of_XbYbZb(node, Xb, X);
    fprintf(fp, ", k=%d, Z=%.16g\n", plane[2], X[2]);
  }

  /* go over plane, with normal */
  for(k=kmin; k<=kmax; k++)
  {
    for(j=jmin; j<=jmax; j++)
    {
      for(i=imin; i<=imax; i++)
      {
        index = Ind_n(i,j,k, va->n);
        fprintf(fp, "%.16g %.16g %.16g\n", p1[index], p2[index], pv[index]);
      }
      if(normal==3) fprintf(fp, "\n");
    }
    if(normal==2 || normal==1) fprintf(fp, "\n");
  }
  fprintf(fp,"\n");
}


/* 2d gnuplot */
void gnuplot_out2d_meshvar(tMesh *mesh, char *name, int It, double T)
{
  tNode *node;
  int vi = Ind(name);
  FILE *fXY, *fXZ, *fYZ;
  char XYfil[1000];
  char XZfil[1000];
  char YZfil[1000];

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
        snprintf(XYfil, 999, "%s/%s_%2d.XY%s",
                 Gets(Par("outdir")),name, p, ns);
        fXY = fopen(XYfil, "a");
        if(!fXY) errorexits("failed opening %s", XYfil);
        write_plane_ascii(node, fXY, 2, ijk, vi, It,T);
        fclose(fXY);
      }

      /* XZ-plane:  Y = Y0 */
      if(ijk[1]>=0)
      {
        snprintf(XZfil, 999, "%s/%s_%2d.XZ%s",
                 Gets(Par("outdir")),name, p, ns);
        fXZ = fopen(XZfil, "a");
        if(!fXZ) errorexits("failed opening %s", XZfil);
        write_plane_ascii(node, fXZ, 1, ijk, vi, It,T);
        fclose(fXZ);
      }

      /* YZ-plane:  X = X0 */
      if(ijk[0]>=0)
      {
        snprintf(YZfil, 999, "%s/%s_%2d.YZ%s",
                 Gets(Par("outdir")),name, p, ns);
        fYZ = fopen(YZfil, "a");
        if(!fYZ) errorexits("failed opening %s", YZfil);
        write_plane_ascii(node, fYZ, 0, ijk, vi, It,T);
        fclose(fYZ);
      }
    }

    /* sysnchronize, so that we write only one node at a time */
    nMPI_barrier();
  } endforlnodes;
}
