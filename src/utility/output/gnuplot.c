/* gnuplot.c */
/* Wolfgang Tichy, Feb. 2019 */

#include "nmesh.h"
#include "output.h"

/* the functions in this file write data in a text format that is understood
   by gnuplot and tgraph */


/* write a single plane */
void write_plane_ascii(tNode *node, FILE *fp, int normal, int plane[],
                       tArray *va, int Iter, double Time)
{
  tArray *p1;
  tArray *p2;
  double Xb[3], X[3];
  double *pv = va->d;
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
        int i1 = i1_norm(i,j,k, normal);
        int i2 = i2_norm(i,j,k, normal);
        int indv = Ind_n(i,j,k, va->n);

        Xb[dir1] = p1->d[i1];
        Xb[dir2] = p2->d[i2];
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
void write_line_ascii(tNode *node, FILE *fp, int dir, int axis[],
                      tArray *va, int Iter, double Time)
{
  tArray *p1;
  double Xb[3], X[3];
  double *pv = va->d;
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
    p1=NULL;
    errorexit("dir has to be 0,1,2");
  }

  /* go over line */
  for(k=kmin; k<=kmax; k++)
  {
    for(j=jmin; j<=jmax; j++)
    {
      for(i=imin; i<=imax; i++)
      {
        int i0 = i0_norm(i,j,k, dir);
        int indv = Ind_n(i,j,k, va->n);

        Xb[dir] = p1->d[i0];
        XYZ_of_XbYbZb(node, Xb, X);
        //fprintf(fp, "%d %d %d: %d %d: ", i,j,k, i0,dir);
        //fprintf(fp, "%g %g: ", Xb[dir], X[dir]);
        //fprintf(fp, "%d %d %d: %d: %d: ", i,j,k, i0, indv);
        fprintf(fp, "%.15g %.15g \n", X[dir], pv[indv]);
      }
    }
  }
  fprintf(fp,"\n");
}


/* 1d output in gnuplot format for one var */
void gnuplot_output1d_meshvar(tMesh *mesh, char *name, int It, double T)
{
  tNode *node;
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
        snprintf(Xfil, 999, "%s/%s.%02dX%s",
                 Gets(Par("outdir")),name, p, ns);
        fX = fopen(Xfil, "a");
        if(!fX) errorexits("failed opening %s", Xfil);
        write_line_ascii(node, fX, 0, ijk, VarA(node, vi), It,T);
        fclose(fX);
      }

      /* Y-axis:  X = X0, Z = Z0 */
      if(ijk[0]>=0 && ijk[2]>=0)
      {
        snprintf(Yfil, 999, "%s/%s.%02dY%s",
                 Gets(Par("outdir")),name, p, ns);
        fY = fopen(Yfil, "a");
        if(!fY) errorexits("failed opening %s", Yfil);
        write_line_ascii(node, fY, 1, ijk, VarA(node, vi), It,T);
        fclose(fY);
      }

      /* Z-axis:  X = X0, Y = Y0 */
      if(ijk[0]>=0 && ijk[1]>=0)
      {
        snprintf(Zfil, 999, "%s/%s.%02dZ%s",
                 Gets(Par("outdir")),name, p, ns);
        fZ = fopen(Zfil, "a");
        if(!fZ) errorexits("failed opening %s", Zfil);
        write_line_ascii(node, fZ, 2, ijk, VarA(node, vi), It,T);
        fclose(fZ);
      }
    }

    /* sysnchronize, so that we write only one node at a time */
    nMPI_barrier();
  } endforlnodes;
}

/* 2d output in gnuplot format for one var */
void gnuplot_output2d_meshvar(tMesh *mesh, char *name, int It, double T)
{
  tNode *node;
  int outd = Par("outdir");
  char *outdir = Gets(outd);
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

      //TODO: use different Xb0 for diff patches
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
        snprintf(XYfil, 999, "%s/%s.%02dXY%s",
                 outdir,name, p, ns);
        fXY = fopen(XYfil, "a");
        if(!fXY) errorexits("failed opening %s", XYfil);
        write_plane_ascii(node, fXY, 2, ijk, VarA(node, vi), It,T);
        fclose(fXY);
      }

      /* XZ-plane:  Y = Y0 */
      if(ijk[1]>=0)
      {
        snprintf(XZfil, 999, "%s/%s.%02dXZ%s",
                 outdir,name, p, ns);
        fXZ = fopen(XZfil, "a");
        if(!fXZ) errorexits("failed opening %s", XZfil);
        write_plane_ascii(node, fXZ, 1, ijk, VarA(node, vi), It,T);
        fclose(fXZ);
      }

      /* YZ-plane:  X = X0 */
      if(ijk[0]>=0)
      {
        snprintf(YZfil, 999, "%s/%s.%02dYZ%s",
                 outdir,name, p, ns);
        fYZ = fopen(YZfil, "a");
        if(!fYZ) errorexits("failed opening %s", YZfil);
        write_plane_ascii(node, fYZ, 0, ijk, VarA(node, vi), It,T);
        fclose(fYZ);
      }
    }

    /* sysnchronize, so that we write only one node at a time */
    nMPI_barrier();
  } endforlnodes;
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
}
