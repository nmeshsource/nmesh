/* VTK_out.c */
/* Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "output.h"

/* for mkdir */
#include <sys/stat.h>
#include <sys/types.h>



/* open file for vtk writing */
FILE *fopen_vtk(char *varname, char *outdir, char *suffix,
                int p, char *nstr, int series)
{
  char filename[1000];
  FILE *fp;

  /* make subdirectory if it does not exist */
  snprintf(filename,999, "%s/%s.%02d%s%s_vtk", outdir,
           varname, p, suffix, nstr);
  fp = fopen(filename, "r");
  if(!fp)
    mkdir(filename, 0777);
  else
    fclose(fp);

  /* open file */
  snprintf(filename,999, "%s/%s.%02d%s%s_vtk/%s.%02d%s%s_%08d.vtk",
	   outdir, varname, p, suffix, nstr,
	   varname, p, suffix, nstr, series);
  fp = fopen(filename, "wb");
  if(!fp)
    errorexits("failed opening %s", filename);

  //PRF;printf(": %s", filename); fflush(stdout);

  /* return non-null file pointer */
  return fp;
}


/* write data part of a VTK file, i.e without the ASCII header */
void write_raw_vtk_data(FILE *fp, double *buffer, int n,
                        int stride, int offset, tOutpars *par)
{
  int i;

  if(par->text) /* ascii vtk data */
  {
    if(par->dbl)
      for(i = 0; i < n; i++)
	fprintf(fp, "%.16g\n", buffer[stride*i+offset]);
    else
      for(i = 0; i < n; i++)
	fprintf(fp, "%.7g\n", (float) buffer[stride*i+offset]);
  }
  else /* assume binary vtk data */
  {
    if(par->dbl)
    {
      double xdouble;
      for(i = 0; i < n; i++)
      {
        xdouble = buffer[stride*i+offset];
        fwrite_big(&xdouble, sizeof(double), 1, fp);
      }
    }
    else
    {
      float xfloat;
      for(i = 0; i < n; i++)
      {
        xfloat = buffer[stride*i+offset];
        fwrite_big(&xfloat, sizeof(float), 1, fp);
      }
    }
  } /* end else */
}


/* 3d output of one array */
void write3d_vtk(tNode *node, FILE *fp, tArray *va, int Iter, double Time,
                 int series, tOutpars *par)
{
  tArray *X[3];
  tArray *Xb[] = { node->Xb[0], node->Xb[1], node->Xb[2] };
  double *pX, *pY, *pZ;
  double *pV = va->d;
  int *n = va->n;
  int n0, n1, n2;

  /* return if var has no memory */
  if(pV==NULL) return;

  /* make room for X,Y,Z */
  X[0] = alloc_array(Xb[0]->n);
  X[1] = alloc_array(Xb[1]->n);
  X[2] = alloc_array(Xb[2]->n);

  /* get arrays with X,Y,Z */
  array_XYZ_of_XbYbZb(node, Xb, X);
  pX = X[0]->d;
  pY = X[1]->d;
  pZ = X[2]->d;

  if(par->arrange_as_1d) /* pretend that all points are along X-dir */
  {
    double X0,Y0,Z0, dX,dY,dZ;

    n0 = n[0]*n[1]*n[2];
    n1 = n2 = 1;
    X0 = Y0 = Z0 = 0.0;
    dX = dY = dZ = 1.0;

    /* write header */
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "variable %s, patch %d, node %s, time %.15g, "
            "Note: uniform grid spacings below are fake\n",
            par->name, par->p, par->nodeloc, Time);
    fprintf(fp, par->binary ? "BINARY" : "ASCII\n");
    fprintf(fp, "\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", n0, n1, n2);
    fprintf(fp, "ORIGIN  %16.9e %16.9e %16.9e\n", X0, Y0, Z0);
    fprintf(fp, "SPACING %16.9e %16.9e %16.9e\n", dX, dY, dZ);
    fprintf(fp, "\n");
    fprintf(fp, "POINT_DATA %d\n", n0*n1*n2);
    fprintf(fp, "SCALARS scalars %s\n", par->dbl ? "double" : "float");
    fprintf(fp, "LOOKUP_TABLE default\n");
    write_raw_vtk_data(fp, pV,n0*n1*n2, 1,0, par);
  }
  else
  {
    n0 = fmin(n[0], X[0]->n[0]);
    n1 = fmin(n[1], X[1]->n[0]);
    n2 = fmin(n[2], X[2]->n[0]);

    /* write header */
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "variable %s, patch %d, node %s, time %.15g\n",
            par->name, par->p, par->nodeloc, Time);
    fprintf(fp, par->binary ? "BINARY" : "ASCII\n");
    fprintf(fp, "\n");
    fprintf(fp, "DATASET RECTILINEAR_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", n0, n1, n2);
    fprintf(fp, "X_COORDINATES %d %s\n", n0, par->dbl ? "double" : "float");
    write_raw_vtk_data(fp, pX,n0, 1,0, par);
    fprintf(fp, "\n\n");
    fprintf(fp, "Y_COORDINATES %d %s\n", n1, par->dbl ? "double" : "float");
    write_raw_vtk_data(fp, pY,n1, 1,0, par);
    fprintf(fp, "\n\n");
    fprintf(fp, "Z_COORDINATES %d %s\n", n2, par->dbl ? "double" : "float");
    write_raw_vtk_data(fp, pZ,n2, 1,0, par);
    fprintf(fp, "\n\n");
    fprintf(fp, "POINT_DATA %d\n", n0*n1*n2);
    fprintf(fp, "SCALARS scalars %s\n", par->dbl ? "double" : "float");
    fprintf(fp, "LOOKUP_TABLE default\n");
    write_raw_vtk_data(fp, pV,n0*n1*n2, 1,0, par);
  }

  free_array(X[2]);
  free_array(X[1]);
  free_array(X[0]);
}
