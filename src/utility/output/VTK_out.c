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
void write3d_vtk(tNode *node, FILE *fp, tArray *va, int Iter,
                 double Time, int series, tOutpars *par)
{
  tArray *X[3];
  tArray *Xb[3];
  double *pX, *pY, *pZ;
  double *pV = va->d;
  int *n = va->n;
  int n0, n1, n2;

  /* return if var has no memory */
  if(pV==NULL) return;

  if(node)
  {
    Xb[0] = node->Xb[0];
    Xb[1] = node->Xb[1];
    Xb[2] = node->Xb[2];

    /* make room for X,Y,Z */
    X[0] = alloc_array(Xb[0]->n);
    X[1] = alloc_array(Xb[1]->n);
    X[2] = alloc_array(Xb[2]->n);

    /* get arrays with X,Y,Z */
    array_XYZ_of_XbYbZb(node, Xb, X);
  }
  else
  {
    int d, m;
    int N[] = { 1,1,1 };
    /* put index ranges in X */
    for(d=0; d<3; d++)
    {
      N[0] = n[d];
      X[d] = alloc_array(N);
      for(m=0; m<n[d]; m++) X[d]->d[m] = m;
    }
  }

  /* pointers to X data */
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
    fprintf(fp, "variable %s, patch %d, node %s, time %.15g\n",
            par->name, par->p, par->nodeloc, Time);
    fprintf(fp, par->text ? "ASCII\n" : "BINARY");
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
    fprintf(fp, par->text ? "ASCII\n" : "BINARY");
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

/* 3d vtk output */
void vtk_output3d_meshvar(tMesh *mesh, char *name, int It, double T)
{
  tNode *node;
  int vi = Ind(name);
  FILE *fp;
  int nseries;
  int vtk      = Getv(Par("3dformat"), "vtk");
  char *outdir = Gets(Par("outdir"));
  tOutpars par[1];

  /* pars we may need for vtk or others */
  par->name          = name;
  par->text          = Getv(Par("3dformat"), "text");
  par->arrange_as_1d = Getv(Par("3dformat"), "arrange_as_1d");
  par->flt           = Getv(Par("3dformat"), "float");
  par->dbl           = Getv(Par("3dformat"), "double");

  /* a number that counts the output */
  nseries = TimeForMeshOutput_di_dt(mesh,Geti(Par("3doutiter")),
                                    Getd(Par("3douttime")));
  /* loop over all nodes */
  forlnodes(mesh, node)
  {
    if(node->dat)
    if(node->dat->v[vi])
    {
      int p = node->pat->p;
      char ns[100];

      /* find string that identifies node */
      node_location_str(node, ns,100);

      /* set some more pars */
      par->nodeloc = ns;
      par->p       = p;

      /* write files */
      if(vtk || 1) /* can do only VTK right now */
      {
        /* VTK output: one file per time step in separate subdirectories */
        fp = fopen_vtk(name, outdir, "XYZ", p, ns, nseries-1);
        write3d_vtk(node, fp, VarA(node, vi), It,T, nseries-1, par);
        fclose(fp);
      }
    }
    /* sysnchronize, so that we write only one node at a time */
    nMPI_barrier();
  } endforlnodes;
}

/* quick array output in vtk format */
void write_array(tNode *node, tArray *va, char *name, int as_1d,
                 int fake_it, double fake_t)
{
  tMesh *mesh;
  FILE *fp;
  int nseries;
  char *outdir;
  tOutpars par[1];
  int p;
  char ns[100];

  if(node)
  {
    mesh = node->pat->mesh;
    outdir = Gets(Par("outdir"));
    p = node->pat->p;

    /* find string that identifies node */
    node_location_str(node, ns,100);
  }
  else
  {
    outdir = ".";
    p = -1;
    ns[0] = 0;
  }

  /* pars we may need for vtk or others */
  par->name          = name;
  par->text          = 1;
  par->arrange_as_1d = as_1d;
  par->flt           = 0;
  par->dbl           = 1;
  par->nodeloc       = ns;
  par->p             = p;

  /* a number that counts the output */
  nseries = fake_it + 1;

  /* write files */
  fp = fopen_vtk(name, outdir, "XYZ", p, ns, nseries-1);
  write3d_vtk(node, fp, va, fake_it,fake_t, nseries-1, par);
  fclose(fp);
}

/* quick var output in vtk format */
void write_var(tNode *node, char *name, int as_1d,
               int fake_it, double fake_t)
{
  tMesh *mesh = node->pat->mesh;
  tArray *va = VarA(node, Ind(name));

  if(va) write_array(node, va, name, as_1d, fake_it,fake_t);
}

/* quick varlist output in vtk format */
void write_vl(tNode *node, tVarList *vl, int as_1d,
              int fake_it, double fake_t)
{
  tMesh *mesh = node->pat->mesh;
  int vli;

  for(vli=0; vli<vl->n; vli++)
  {
    int vi = vl->index[vli];
    char *vname = VarName(vi);
    write_var(node, vname, as_1d, fake_it, fake_t);
  }
}


/* 2d output of one array */
void write_plane_vtk(tNode *node, FILE *fp, int normal, int plane[],
                     tArray *va, int Iter, double Time,
                     int series, tOutpars *par)
{
  tArray *X[3];
  tArray *Xb[] = { node->Xb[0], node->Xb[1], node->Xb[2] };
  double *pX, *pY, *pZ, *buf;
  int *n = va->n;
  int n0, n1, n2, off0, off1, off2, i,j,k, cnt;

  /* return if var has no memory */
  if(!va || va->d==NULL) return;

  /* make room for X,Y,Z and buf for var-plane data */
  X[0] = alloc_array(Xb[0]->n);
  X[1] = alloc_array(Xb[1]->n);
  X[2] = alloc_array(Xb[2]->n);
  buf = dmalloc(va->N);

  /* get arrays with X,Y,Z */
  array_XYZ_of_XbYbZb(node, Xb, X);
  pX = X[0]->d;
  pY = X[1]->d;
  pZ = X[2]->d;

  /* offsets */
  off0 = off1 = off2 = 0;

  /* number of points */
  n0 = fmin(n[0], X[0]->n[0]);
  n1 = fmin(n[1], X[1]->n[0]);
  n2 = fmin(n[2], X[2]->n[0]);

  /* write var plane data in buf */
  cnt = 0;
  if(normal==0)
  {
    n0 = 1;
    off0 = plane[0];
    for(k=0; k<n2; k++)
      for(j=0; j<n1; j++)
        buf[cnt++] = va->d[Ind_n(off0,j,k, n)];
  }
  else if(normal==1)
  {
    n1 = 1;
    off1 = plane[1];
    for(k=0; k<n2; k++)
      for(i=0; i<n0; i++)
        buf[cnt++] = va->d[Ind_n(i,off1,k, n)];
  }
  else
  {
    n2 = 1;
    off2 = plane[2];
    for(j=0; j<n1; j++)
      for(i=0; i<n0; i++)
        buf[cnt++] = va->d[Ind_n(i,j,off2, n)];
  }

  /* write header and data */
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "variable %s, patch %d, node %s, time %.15g\n",
          par->name, par->p, par->nodeloc, Time);
  fprintf(fp, par->text ? "ASCII\n" : "BINARY");
  fprintf(fp, "\n");
  fprintf(fp, "DATASET RECTILINEAR_GRID\n");
  fprintf(fp, "DIMENSIONS %d %d %d\n", n0, n1, n2);
  fprintf(fp, "X_COORDINATES %d %s\n", n0, par->dbl ? "double" : "float");
  write_raw_vtk_data(fp, pX,n0, 1,off0, par);
  fprintf(fp, "\n\n");
  fprintf(fp, "Y_COORDINATES %d %s\n", n1, par->dbl ? "double" : "float");
  write_raw_vtk_data(fp, pY,n1, 1,off1, par);
  fprintf(fp, "\n\n");
  fprintf(fp, "Z_COORDINATES %d %s\n", n2, par->dbl ? "double" : "float");
  write_raw_vtk_data(fp, pZ,n2, 1,off2, par);
  fprintf(fp, "\n\n");
  fprintf(fp, "POINT_DATA %d\n", n0*n1*n2);
  fprintf(fp, "SCALARS scalars %s\n", par->dbl ? "double" : "float");
  fprintf(fp, "LOOKUP_TABLE default\n");
  write_raw_vtk_data(fp, buf,n0*n1*n2, 1,0, par);

  free(buf);
  free_array(X[2]);
  free_array(X[1]);
  free_array(X[0]);
}

/* 2d vtk output */
void vtk_output2d_meshvar(tMesh *mesh, char *name, int It, double T)
{
  tNode *node;
  int vi = Ind(name);
  FILE *fp;
  int nseries;
  //int vtk      = Getv(Par("2dformat"), "vtk");
  char *outdir = Gets(Par("outdir"));
  tOutpars par[1];

  /* pars we may need for vtk or others */
  par->name          = name;
  par->text          = Getv(Par("2dformat"), "text");
  par->arrange_as_1d = 0; // Getv(Par("2dformat"), "arrange_as_1d");
  par->flt           = Getv(Par("2dformat"), "float");
  par->dbl           = Getv(Par("2dformat"), "double");

  /* a number that counts the output */
  nseries = TimeForMeshOutput_di_dt(mesh,Geti(Par("2doutiter")),
                                    Getd(Par("2douttime")));
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

      /* set some more pars */
      par->nodeloc = ns;
      par->p       = p;

      /* write files */
      /* XY-plane:  Z = Z0 */
      if(ijk[2]>=0)
      {
        /* VTK output: one file per time step in separate subdirectories */
        fp = fopen_vtk(name, outdir, "XY", p, ns, nseries-1);
        write_plane_vtk(node, fp, 2,ijk, VarA(node, vi), It,T, nseries-1, par);
        fclose(fp);
      }
      /* XZ-plane:  Y = Y0 */
      if(ijk[1]>=0)
      {
        /* VTK output: one file per time step in separate subdirectories */
        fp = fopen_vtk(name, outdir, "XZ", p, ns, nseries-1);
        write_plane_vtk(node, fp, 1,ijk, VarA(node, vi), It,T, nseries-1, par);
        fclose(fp);
      }
      /* YZ-plane:  X = X0 */
      if(ijk[0]>=0)
      {
        /* VTK output: one file per time step in separate subdirectories */
        fp = fopen_vtk(name, outdir, "YZ", p, ns, nseries-1);
        write_plane_vtk(node, fp, 0,ijk, VarA(node, vi), It,T, nseries-1, par);
        fclose(fp);
      }
    }
    /* sysnchronize, so that we write only one node at a time */
    nMPI_barrier();
  } endforlnodes;
}
