/* VTK_out.c */
/* Wolfgang Tichy 2/2019 */

#include "nmesh.h"
#include "output.h"

/* for mkdir */
#include <sys/stat.h>
#include <sys/types.h>


/* write data part of a VTK file, i.e without the ASCII header */
void write_raw_vtk_data(FILE *fp, double *buffer, int n, int stride, int offset,
		        int dbl, int flt, int text, int binary)
{
  int i;

  if(!text && !binary)
    errorexit("write_raw_vtk_data: pick text or binary format");
  if(text && binary)
    errorexit("write_raw_vtk_data: pick either text or binary format");
  if(dbl && flt)
    errorexit("write_raw_vtk_data: pick either double or float format");

  if(text)
  {
    if(dbl) 
      for (i = 0; i < n; i++)
	fprintf(fp, "%.16g\n", buffer[stride*i+offset]);
    else
      for (i = 0; i < n; i++)
	fprintf(fp, "%.7g\n", (float) buffer[stride*i+offset]);
  }

  if(binary)
  {
    if(dbl)
    {
      double xdouble;
      for(i = 0; i < n; i++)
      {
        xdouble = buffer[stride*i+offset];
        fwrite_big(&xdouble, sizeof(double), 1, fp);
      }
    }
    if(flt)
    {
      float xfloat;
      for (i = 0; i < n; i++)
      {
        xfloat = buffer[stride*i+offset];
        fwrite_big(&xfloat, sizeof(float), 1, fp);
      }
    }
  }
}

/* open file for vtk writing */
FILE *fopen_vtk(char *varname, char *suffix, char *outdir,
                int p, char *nstr, int series)
{
  char filename[1000];
  FILE *fp;

  /* make subdirectory if it does not exist */
  snprintf(filename,999, "%s/%s.%s%d_vtk", outdir, varname, suffix, p);
  fp = fopen(filename, "r");
  if(!fp)
    mkdir(filename, 0777);  
  else
    fclose(fp);

  /* open file */
  snprintf(filename,999, "%s/%s.%s%d_vtk/%s.%s%d_%08d.vtk",
	   outdir, varname, suffix, p,
	   varname, suffix, p, series);
  fp = fopen(filename, "wb");
  if(!fp) 
    errorexits("failed opening %s", filename);
  
  /* return non-null file pointer */
  return fp;
}



/* 3d output of one array */
void write3d_vtk(tNode *node, FILE *fp, tArray *va, int Iter, double Time)
{
  tGrid *grid = box->grid; 
  char *outdir      = Gets(Par("outdir"));
  int text          = Getv(Par("3dformat"), "text");
  int binary        = Getv(Par("3dformat"), "binary");
  int vtk           = Getv(Par("3dformat"), "vtk");
  int arrange_as_1d = Getv(Par("3dformat"), "arrange_as_1d");
  int addpoints     = Getv(Par("3dformat"), "addpoints");
  int flt           = Getv(Par("3dformat"), "float");
  int dbl           = Getv(Par("3dformat"), "double");
  FILE *fp;
  int nseries;
  tArray *X[3];
  tArray *Xb[] = { node->Xb[0], node->Xb[1], node->Xb[2] };
  double *pV = va->d;

  /* return if var has no memory */
  if(pV==NULL) return;  

  /* make room for X,Y,Z */
  X[0] = alloc_array(Xb[0]->n);
  X[1] = alloc_array(Xb[1]->n);
  X[2] = alloc_array(Xb[2]->n);

  /* get arrays with X,Y,Z */
  array_XYZ_of_XbYbZb(node, Xb, X);

  /* parameter defaults */
  if (!flt && !dbl) flt = 1;
  if (!text && !binary) binary = 1;

  /* VTK */
  /* output one file per time step in separate subdirectories */
  if(vtk)
  {
    /* open file (returns non-null file pointer) */
    fp = fopen_vtk(name, "XYZ", box->b, nseries-1);

    /* FIXME: only fakepoints works so far!!! */
    fakepoints=1;

    if(fakepoints) /* put data on a fake grid with uniform grid spacings dX,dY,dZ */
    {
      double X0 = box->bbox[0];
      double Y0 = box->bbox[2];
      double Z0 = box->bbox[4];
      double X1 = box->bbox[1];
      double Y1 = box->bbox[3];
      double Z1 = box->bbox[5];
      double dX = fabs(X1-X0)/(n1);
      double dY = fabs(Y1-Y0)/(n2);
      double dZ = fabs(Z1-Z0)/(n3);

      if(arrange_as_1d) /* pretend that all points are along X-dir */
      {
        n1 = n1*n2*n3;
        n2 = n3 = 1;
        X0 = Y0 = Z0 = 0.0;
        dX = dY = dZ = 1.0;
      }

      /* write header */
      fprintf(fp, "# vtk DataFile Version 2.0\n");
      fprintf(fp, "variable %s, box %d, time %.16g, "
              "Note: uniform grid spacings below are fake\n", 
              name, box->b, grid->time);
      fprintf(fp, binary ? "BINARY" : "ASCII\n");
      fprintf(fp, "\n");
      fprintf(fp, "DATASET STRUCTURED_POINTS\n");
      fprintf(fp, "DIMENSIONS %d %d %d\n", n1, n2, n3);
      fprintf(fp, "ORIGIN  %16.9e %16.9e %16.9e\n", X0, Y0, Z0);
      fprintf(fp, "SPACING %16.9e %16.9e %16.9e\n", dX, dY, dZ);
      fprintf(fp, "\n");
      fprintf(fp, "POINT_DATA %d\n", n1*n2*n3);
      fprintf(fp, "SCALARS scalars %s\n", dbl ? "double" : "float");
      fprintf(fp, "LOOKUP_TABLE default\n");
      /* write data,
         has to be in the file right after the \n of the last header line */
      write_raw_vtk_data(fp, pV, n1*n2*n3,1,0, dbl, flt, text, binary);
    }
    else
    {
      if(addpoints) /* add grid in X,Y,Z coords. */
      {
        errorexit("implement 3d vtk output with real grid point coordinates");
      }
      else /* use RECTILINEAR_GRID */
      {
        /* write header */
        fprintf(fp, "# vtk DataFile Version 2.0\n");
        fprintf(fp, "variable %s, box %d, time %.16g\n", 
                name, box->b, grid->time);
        fprintf(fp, binary ? "BINARY" : "ASCII\n");
        fprintf(fp, "\n");
        fprintf(fp, "DATASET RECTILINEAR_GRID\n");
        fprintf(fp, "DIMENSIONS %d %d %d\n", n1, n2, n3);
        fprintf(fp, "X_COORDINATES %d %s\n", n1, dbl ? "double" : "float");
        write_raw_vtk_data(fp, pX, n1,1,0, 1, 0, text, binary);
        fprintf(fp, "\n");
        fprintf(fp, "Y_COORDINATES %d %s\n", n2, dbl ? "double" : "float");
        write_raw_vtk_data(fp, pY, n2,n1,0, 1, 0, text, binary);
        fprintf(fp, "\n");
        fprintf(fp, "Z_COORDINATES %d %s\n", n3, dbl ? "double" : "float");
        write_raw_vtk_data(fp, pZ, n3,n1*n2,0, 1, 0, text, binary);
        fprintf(fp, "\n");
        fprintf(fp, "POINT_DATA %d\n", n1*n2*n3);
        fprintf(fp, "SCALARS scalars %s\n", dbl ? "double" : "float");
        fprintf(fp, "LOOKUP_TABLE default\n");
        /* write data,
           has to be in the file right after the \n of the last header line */
        write_raw_vtk_data(fp, pV, n1*n2*n3,1,0, dbl, flt, text, binary);
        errorexit("this is not in proper vtk format");
      }
    }
    fclose(fp);
  }
}
