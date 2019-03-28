/* XDMF_out.c */
/* Wolfgang Tichy 3/2019 */

#include "nmesh.h"
#include "output.h"


/* VisIt and others (maybe paraview) can read XDMF format.
   Our output in XDMF format consists of a .xml file and a .bin file. The
   actual data is in the .bin file, while the .xml file only contains a
   description of the data in XML format.

   Some description about XDMF is on
   http://www.xdmf.org
   http://www.xdmf.org/index.php/XDMF_Model_and_Format
   http://www.paraview.org/Wiki/ParaView/Data_formats

   What we have here is quite similar to what is in bamps 4.0. But I added
   "Time Value" and "Spatial" to B_spatial. Now VisIt shows the correct time
   and paraview can show more than t=0. */



/* XML format strings to make .xmf files using fprintf,
   based on bamps and https://www.paraview.org/Wiki/ParaView/Data_formats */
char *B_head =
  "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
  "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.1\">\n"
  "  <Domain>\n";
char *E_head =
  "  </Domain>\n"
  "</Xdmf>\n";

char *B_temporal =
  "    <Grid CollectionType=\"Temporal\" GridType=\"Collection\" Name=\"TCollection\">\n"
  "      <Geometry Type=\"None\"/>\n"
  "      <Topology Dimensions=\"0\" Type=\"NoTopology\"/>\n";
char *E_temporal =
  "    </Grid>\n";

char *B_spatial =
  "      <Grid CollectionType=\"Spatial\" GridType=\"Collection\" Name=\"SCollection\">\n"
  "        <Time Value=\"%.9f\"/>\n"
  "        <Geometry Type=\"None\"/>\n"
  "        <Topology Dimensions=\"0\" Type=\"NoTopology\"/>\n";
char *E_spatial =
  "      </Grid>\n";

char *B_E_grid =
  "        <Grid Name=\"%s\">\n"
  "          <Time Value=\"%.9f\"/>\n"
  "          <Geometry Type=\"XYZ\">\n"
  "            <DataItem DataType=\"Float\" Dimensions=\"%d %d\" Format=\"%s\" Seek=\"%d\" Precision=\"4\">\n"
  "              %s\n"
  "            </DataItem>\n"
  "          </Geometry>\n"
  "          <Topology Dimensions=\"%d %d %d\" Type=\"3DSMesh\"/>\n"
  "          <Attribute Center=\"Node\" Name=\"%s\" Type=\"Scalar\">\n"
  "            <DataItem DataType=\"Float\" Dimensions=\"%d %d %d\" Format=\"%s\" Seek=\"%d\" Precision=\"4\">\n"
  "              %s\n"
  "            </DataItem>\n"
  "          </Attribute>\n"
  "        </Grid>\n";



/* open file xmf file with XML description and position file pointer */
FILE *fopen_xdmf_xmf(char *varname, char *outdir, char *suffix, double time)
{
  FILE *fp;
  char fname[1000];
  long offset;

  /* name of XML file */
  snprintf(fname, 1000, "%s/%s.%s.xmf", outdir, varname, suffix);

  /* open file such that we can append and seek backwards */
  fp = fopen(fname, "r+");

  /* create file if it could not be opened */
  if(!fp)
  {
    fp = fopen(fname, "w+");
    if(!fp) errorexits("Cannot open %s for writing", fname);

    /* write fixed part of XML header into new file */
    fprintf(fp, "%s", B_head);
    fprintf(fp, "%s", B_temporal);
    fprintf(fp, "%s", E_temporal);
    fprintf(fp, "%s", E_head);
  }

  /* we want to append more data, which requires us to remove
     the last E_temporal and E_head */
  offset = strlen(E_temporal) + strlen(E_head);
  fseek(fp, -offset, SEEK_END);

  /* start new collection of nodes */
  fprintf(fp, B_spatial, time);

  return fp;
}

/* write ends and close .xmf file */
void fclose_xdmf_xmf(FILE *fp)
{
  /* close collections and header */
  fprintf(fp, "%s", E_spatial);
  fprintf(fp, "%s", E_temporal);
  fprintf(fp, "%s", E_head);
  fclose(fp);
}

/* open file to add more nodes still with the same Time Value */
FILE *fopen_add_spatial_xdmf_xmf(char *varname, char *outdir, char *suffix)
{
  FILE *fp;
  char fname[1000];
  long offset;

  /* name of XML file */
  snprintf(fname, 1000, "%s/%s.%s.xmf", outdir, varname, suffix);

  /* open file such that we can append and seek backwards */
  fp = fopen(fname, "r+");
  if(!fp) errorexit("cannot add if file was never created with fopen_xdmf_xmf");

  /* remove E_spatial, E_temporal, E_head */
  offset = strlen(E_spatial) + strlen(E_temporal) + strlen(E_head);
  fseek(fp, -offset, SEEK_END);

  return fp;
}

/* open a .bin file with raw binary data */
FILE *fopen_bin(char *varname, char *outdir, char *suffix)
{
  FILE *fp;
  char fname[1000];
  snprintf(fname, 1000, "%s/%s.%s.bin", outdir, varname, suffix);

  fp = fopen(fname, "a");
  if(!fp) errorexits("Cannot open %s for writing", fname);

  return fp;
}



/* write XML grid description into .xmf file with file pointer fp */
void write_xdmf_xmf(FILE *fp, int voffset, int xyzoffset,
                    char *vname, char *suffix,
		    char *nodename, double time,
		    int n[3], int bin, int dbl)
{
  char fname[1000];
  char fname_xyz[1000];
  int np = n[0] * n[1] * n[2];
  char *format = (bin) ? "Binary" : "XML";

  /* filenames for field and also xyz data*/
  snprintf(fname, 1000, "%s.%s.bin", vname, suffix);
  snprintf(fname_xyz, 1000, "xyz.%s.bin", suffix);

  /* <Grid> info to .xmf file */
  fprintf(fp, B_E_grid,
	  nodename,
	  time,
	  np, 3,
	  format,
	  xyzoffset,
	  fname_xyz,
	  n[2], n[1], n[0],
	  vname,
	  n[2], n[1], n[0],
	  format,
	  voffset,
          fname);
}


/* 2d output */
void output2d_xdmf(tVarList *vl, int It, double T)
{
  tMesh *mesh = vl->mesh;
  int outd = Par("outdir");
  char *outdir = Gets(outd);

  write_plane_xdmf(vl, 0, outdir, T);
  write_plane_xdmf(vl, 1, outdir, T);
  write_plane_xdmf(vl, 2, outdir, T);
}

/* output varlist in XDMF format in one plane */
void write_plane_xdmf(tVarList *vl, int norm, char *outdir, double Time)
{
  tMesh *mesh = vl->mesh;
  int bin = 1; /* we can only do binary output right now */
  int dbl = 0; /* we output float not double */
  int voffset, xyzoffset;
  FILE *fpxmf, *fpbin;
  char ndname[100];
  char *suffix[] = { "yz", "xz", "xy" };
  double X0[] = { Getd(Par("outputX0")),
                  Getd(Par("outputY0")),
                  Getd(Par("outputZ0")) };
  int vli;

  /* loop over varlist */
  for(vli=0; vli<vl->n; vli++)
  {
    int vi = vl->index[vli];
    char *vname = VarName(vi);
    int rk;
    int myid;

    /* MPI motivated loop to assign work */
    for(rk=0; rk<nMPI_size(); rk++)
    {
      /* do work when it is my turn */
      if(rk == nMPI_rank())
      {
        if(Rank0) /* open xmf to start a new spatial series */
          fpxmf = fopen_xdmf_xmf(vname, outdir, suffix[norm], Time);
        else /* just add to the same spatial series */
          fpxmf = fopen_add_spatial_xdmf_xmf(vname, outdir, suffix[norm]);

        /* open binary file */
        fpbin = fopen_bin(vname, outdir, suffix[norm]);

        /* loop over all leaf nodes */
        formylnodes_noomp(mesh, myid)
        {
          tNode *node = MyNode(mesh, myid);

          /* do something only if this proc has dat */
          if(node->dat)
          if(node->dat->v[vi])
          {
            int n[] = { node->n[0], node->n[1], node->n[2] };
            int plane[3];
            intList *plist;
            int normal;

            /* find indices of nearest, if all are negative, node does not have
               outputX0, outputY0, outputZ0 */
            nearest_lowernode_ijk_of_XYZ(node, plane, X0);
            normal = approxXYZnormal_of_xyznormal(node, norm);

            if(normal>=0)
            if(plane[normal]>=0)
            {
              /* node name and n for plane */
              nodename(node, ndname,99);
              n[normal] = 1;

              /* list of points in plane */
              plist = pointindexList_plane(node, normal, plane);

              /* write binary data in var */
              voffset = ftell(fpbin);
              write_buffer_idx(Vard(node,vi), plist, dbl, fpbin);

              /* write point's x,y,z coordinates */
              if(vli==0) /* write xyz only for first var in list */
              {
                int ix = Ind("x");
                double *px = Vard(node, ix);
                double *py = Vard(node, ix+1);
                double *pz = Vard(node, ix+2);
                FILE *fpxyz = fopen_bin("xyz", outdir, suffix[norm]);

                write_3buffers_idx(px,py,pz, plist, dbl, fpxyz);
                fclose(fpxyz);
              }

              /* we wrote 3 things (x,y,z) for each var */
              xyzoffset = 3 * voffset;

              /* write grid information into xmf file */
              write_xdmf_xmf(fpxmf, voffset, xyzoffset, vname, suffix[norm],
                             ndname, Time, n, bin, dbl);

              free_intList(plist);
            }
          }
        }
        /* close files on this proc now */
        fclose(fpbin);
        fclose_xdmf_xmf(fpxmf);
      }
      /* wait until everyone is here */
      nMPI_barrier();
    } /* end rk-loop */
  }
}


/* output varlist in XDMF form for the entire volume */
void output3d_xdmf(tVarList *vl, int It, double Time)
{
  tMesh *mesh = vl->mesh;
  int outd = Par("outdir");
  char *outdir = Gets(outd);
  int bin = 1; /* we can only do binary output right now */
  int dbl = 0; /* we output float not double */
  int voffset, xyzoffset;
  FILE *fpxmf, *fpbin;
  char ndname[100];
  char *suffix = "xyz";
  int vli;

  /* loop over varlist */
  for(vli=0; vli<vl->n; vli++)
  {
    int vi = vl->index[vli];
    char *vname = VarName(vi);
    int rk;
    int myid;

    /* MPI motivated loop to assign work */
    for(rk=0; rk<nMPI_size(); rk++)
    {
      /* do work when it is my turn */
      if(rk == nMPI_rank())
      {
        if(Rank0) /* open xmf to start a new spatial series */
          fpxmf = fopen_xdmf_xmf(vname, outdir, suffix, Time);
        else /* just add to the same spatial series */
          fpxmf = fopen_add_spatial_xdmf_xmf(vname, outdir, suffix);

        /* open binary file */
        fpbin = fopen_bin(vname, outdir, suffix);

        /* loop over all leaf nodes */
        formylnodes_noomp(mesh, myid)
        {
          tNode *node = MyNode(mesh, myid);

          /* do something only if this proc has dat */
          if(node->dat)
          if(node->dat->v[vi])
          {
            int np = node->np;
            int n[] = { node->n[0], node->n[1], node->n[2] };

            {
              /* node name and n for plane */
              nodename(node, ndname,99);

              /* write binary data in var */
              voffset = ftell(fpbin);
              write_buffer(Vard(node,vi), np, dbl, fpbin);

              /* write point's x,y,z coordinates */
              if(vli==0) /* write xyz only for first var in list */
              {
                int ix = Ind("x");
                double *px = Vard(node, ix);
                double *py = Vard(node, ix+1);
                double *pz = Vard(node, ix+2);
                FILE *fpxyz = fopen_bin("xyz", outdir, suffix);

                write_3buffers(px,py,pz, np, dbl, fpxyz);
                fclose(fpxyz);
              }

              /* we wrote 3 things (x,y,z) for each var */
              xyzoffset = 3 * voffset;

              /* write grid information into xmf file */
              write_xdmf_xmf(fpxmf, voffset, xyzoffset, vname, suffix,
                             ndname, Time, n, bin, dbl);
            }
          }
        }
        /* close files on this proc now */
        fclose(fpbin);
        fclose_xdmf_xmf(fpxmf);
      }
      /* wait until everyone is here */
      nMPI_barrier();
    } /* end rk-loop */
  }
}



/******************************************************************/
/* General functions to do output for an entire buffer */
/******************************************************************/

/* write a double buffer as double or float */
size_t write_buffer(const double *buf, int buflen, int dbl, FILE *fp)
{
  int i;
  size_t cnt;

  if(dbl)
    cnt = fwrite(buf, sizeof(double), buflen, fp);
  else
  {
    cnt = 0;
    for(i=0; i<buflen; i++)
    {
      float fval = buf[i];
      cnt += fwrite(&fval, sizeof(float), 1, fp);
    }
  }
  return cnt;
}

/* write 3 double buffers as double or float */
size_t write_3buffers(const double *b1, const double *b2, const double *b3,
                      int buflen, int dbl, FILE *fp)
{
  int i;
  size_t cnt;

  cnt = 0;
  for(i=0; i<buflen; i++)
  {
    if(dbl)
    {
      double dval;
      dval = b1[i];
      cnt += fwrite(&dval, sizeof(double), 1, fp);
      dval = b2[i];
      cnt += fwrite(&dval, sizeof(double), 1, fp);
      dval = b3[i];
      cnt += fwrite(&dval, sizeof(double), 1, fp);
    }
    else
    {
      float fval;
      fval = b1[i];
      cnt += fwrite(&fval, sizeof(float), 1, fp);
      fval = b2[i];
      cnt += fwrite(&fval, sizeof(float), 1, fp);
      fval = b3[i];
      cnt += fwrite(&fval, sizeof(float), 1, fp);
    }
  }
  return cnt;
}


/******************************************************************/
/* General functions to do output for a point index list */
/******************************************************************/

/* write a buffer at all indices in idx */
size_t write_buffer_idx(const double *buf, intList *idx, int dbl, FILE *fp)
{
  int i;
  size_t cnt;

  cnt = 0;
  forList(idx, i)
  {
    int ind = ListEntry(idx, i);
    double dval = buf[ind];
    float  fval = dval;

    if(dbl)
      cnt += fwrite(&dval, sizeof(double), 1, fp);
    else
      cnt += fwrite(&fval, sizeof(float), 1, fp);
  }
  return cnt;
}

/* write 3 buffers at all indices in idx */
size_t write_3buffers_idx(const double *b1, const double *b2, const double *b3,
                          intList *idx, int dbl, FILE *fp)
{
  int i;
  size_t cnt;

  cnt = 0;
  forList(idx, i)
  {
    int ind = (ListEntry(idx, i));
    if(dbl)
    {
      double dval;
      dval = b1[ind];
      cnt += fwrite(&dval, sizeof(double), 1, fp);
      dval = b2[ind];
      cnt += fwrite(&dval, sizeof(double), 1, fp);
      dval = b3[ind];
      cnt += fwrite(&dval, sizeof(double), 1, fp);
    }
    else
    {
      float fval;
      fval = b1[ind];
      cnt += fwrite(&fval, sizeof(float), 1, fp);
      fval = b2[ind];
      cnt += fwrite(&fval, sizeof(float), 1, fp);
      fval = b3[ind];
      cnt += fwrite(&fval, sizeof(float), 1, fp);
    }
  }
  return cnt;
}


/******************************************************************/
/* not needed ??? */
/******************************************************************/

/* write a buffer at all indices in idx */
size_t fwrite_buffer_idx(const void *ptr, size_t size,
                         intList *idx, FILE *fp)
{
  const char *buf = ptr;
  int i, cnt;

  cnt = 0;
  forList(idx, i)
  {
    int ind = (ListEntry(idx, i))*size;
    cnt += fwrite(buf+ind, size, 1, fp);
  }
  return cnt;
}

/* write 3 buffers at all indices in idx */
size_t fwrite_3buffers_idx(const void *p1, const void *p2, const void *p3,
                           size_t size, intList *idx, FILE *fp)
{
  const char *buf1 = p1;
  const char *buf2 = p2;
  const char *buf3 = p3;
  int i;
  size_t ind, cnt;

  cnt = 0;
  forList(idx, i)
  {
    ind = (ListEntry(idx, i))*size;
    cnt += fwrite(buf1+ind, size, 1, fp);
    cnt += fwrite(buf2+ind, size, 1, fp);
    cnt += fwrite(buf3+ind, size, 1, fp);
  }
  return cnt;
}
