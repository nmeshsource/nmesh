/* XDMF_out.c */
/* Wolfgang Tichy 3/2019 */

#include "nmesh.h"
#include "output.h"

#include <unistd.h>      /* for truncate */

/* VisIt and ParaView can read XDMF format.
   Our output in XDMF format consists of a .xml file and a .bin file. The
   actual data is in the .bin file, while the .xml file only contains a
   description of the data in XML format.

   Some description about XDMF is on
   http://www.xdmf.org
   http://www.xdmf.org/index.php/XDMF_Model_and_Format
   http://www.paraview.org/Wiki/ParaView/Data_formats */

/* include _xmf and _smf format strings: */
#include "XDMF_out.h"

/* pointers that can point to the _xmf or _smf strings in XDMF_out.h */
static const char *B_head;
static const char *E_head;
static const char *B_temporal;
static const char *E_temporal;
static const char *B_spatial;
static const char *E_spatial;
static const char *B_E_grid;
/* set the pointers to the format strings above */
void xdmf_set_format_strings(int mode)
{
  switch(mode)
  {
  case 1: /* use WT's simplified text format */
    B_head     = B_head_smf;
    E_head     = E_head_smf;
    B_temporal = B_temporal_smf;
    E_temporal = E_temporal_smf;
    B_spatial  = B_spatial_smf;
    E_spatial  = E_spatial_smf;
    B_E_grid   = B_E_grid_smf;
    break;
  case 0:
  default: /* use regular xmf format */
    B_head     = B_head_xmf;
    E_head     = E_head_xmf;
    B_temporal = B_temporal_xmf;
    E_temporal = E_temporal_xmf;
    B_spatial  = B_spatial_xmf;
    E_spatial  = E_spatial_xmf;
    B_E_grid   = B_E_grid_xmf;
  }
}


/* write filename of xml file into string fname */
void filename_xdmf_xmf(char *varname, const char *outdir, const char *suffix,
                       int nfname, char *fname)
{
  const char *ext;

  /* set file extension */
  if(B_E_grid==B_E_grid_xmf)      ext = "xmf";
  else if(B_E_grid==B_E_grid_smf) ext = "txt";
  else errorexit("B_E_grid is not set correctly");

  /* name of XML file */
  snprintf(fname, nfname, "%s/%s.%s.%s", outdir, varname, suffix, ext);
}

/* open file xmf file with XML description and position file pointer */
FILE *fopen_xdmf_xmf(char *varname, const char *outdir, const char *suffix,
                     double time, char *IObuf, size_t IObufsiz)
{
  FILE *fp;
  char fname[1000];
  long offset;

  /* put name of XML file into fname */
  filename_xdmf_xmf(varname, outdir, suffix, 1000, fname);

  /* open file such that we can append and seek backwards */
  fp = fopen(fname, "r+"); //with "a" fseek below would fail
  if(fp) /* file actually exists */
  {
    /* attach IO buffer */
    if(IObufsiz) setvbuf(fp, IObuf, _IOFBF, IObufsiz);

    /* we want to append more data, which requires us to remove
       the last E_temporal and E_head */
    offset = strlen(E_temporal) + strlen(E_head);
    fseek(fp, -offset, SEEK_END);
  }
  else /* if file does not exist yet */
  {
    fp = fopen(fname, "w");
    if(!fp) errorexits("Cannot open %s for writing", fname);

    /* attach IO buffer */
    if(IObufsiz) setvbuf(fp, IObuf, _IOFBF, IObufsiz);

    /* write fixed part of XML header into new file */
    fprintf(fp, "%s", B_head);
    fprintf(fp, "%s", B_temporal);
  }

  /* start new collection of nodes */
  fprintf(fp, B_spatial, time);

  return fp;
}

/* write ends and close .xmf file */
void fclose_xdmf_xmf(FILE *fp, int syncmode, int E_markers)
{
  /* close collections and header */
  if(E_markers)
  {
    fprintf(fp, "%s", E_spatial);
    fprintf(fp, "%s", E_temporal);
    fprintf(fp, "%s", E_head);
  }
  fclose_sync_mode(fp, syncmode);
}

/* Remove an empty B_spatial,E_spatial block:
   In case we have written B_spatial and E_spatial with nothing in between
   it is better to remove it. */
void rm_empty_spatial_xdmf_xmf(char *varname,
                               const char *outdir, const char *suffix,
                               double time, int syncmode)
{
  FILE *fp;
  char fname[1000];
  long nbytes, rmbytes;
  char str[10000];

  /* put name of XML file into fname */
  filename_xdmf_xmf(varname, outdir, suffix, 1000, fname);

  /* current number of bytes in file */
  nbytes = nbytes_infile_name(fname);
  //PRF;printf(": %s has %ld bytes\n", fname, nbytes);
  /* make str with stuff we want to remove */
  rmbytes = 0;
  rmbytes += sprintf(str+rmbytes, B_spatial, time);
  rmbytes += sprintf(str+rmbytes, "%s", E_spatial);
  rmbytes += sprintf(str+rmbytes, "%s", E_temporal);
  rmbytes += sprintf(str+rmbytes, "%s", E_head);

  /* remove last rmbytes */
  nbytes -= rmbytes;
  truncate(fname, nbytes);
  //PRF;printf(": remove %ld bytes, i.e. truncate to %ld\n", rmbytes, nbytes);

  /* open file and append the now missing end markers */
  fp = fopen(fname, "a");
  if(!fp) errorexits("cannot add to %s if file was never created with "
                     "fopen_xdmf_xmf", fname);
  fprintf(fp, "%s", E_temporal);
  fprintf(fp, "%s", E_head);
  fclose_sync_mode(fp, syncmode);
}

/* open file to add more nodes still with the same Time Value */
FILE *fopen_add_spatial_xdmf_xmf(char *varname,
                                 const char *outdir, const char *suffix,
                                 char *IObuf, size_t IObufsiz)
{
  FILE *fp;
  char fname[1000];
  const char *ext;

  /* set file extension */
  if(B_E_grid==B_E_grid_xmf)      ext = "xmf";
  else if(B_E_grid==B_E_grid_smf) ext = "txt";
  else errorexit("B_E_grid is not set correctly");

  /* name of XML file */
  snprintf(fname, 1000, "%s/%s.%s.%s", outdir, varname, suffix, ext);

  /* open file such that we can append and seek backwards */
  fp = fopen(fname, "a");
  if(!fp) errorexits("cannot add to %s if file was never created with "
                     "fopen_xdmf_xmf", fname);

  /* attach IO buffer */
  if(IObufsiz) setvbuf(fp, IObuf, _IOFBF, IObufsiz);

  return fp;
}

/* open a .bin file with raw binary data */
FILE *fopen_bin(const char *varname, const char *outdir, const char *suffix,
                char *IObuf, size_t IObufsiz)
{
  FILE *fp;
  char fname[1000];
  snprintf(fname, 1000, "%s/%s.%s.bin", outdir, varname, suffix);

  fp = fopen(fname, "a");
  if(!fp) errorexits("Cannot open %s for writing", fname);

  /* attach IO buffer */
  if(IObufsiz) setvbuf(fp, IObuf, _IOFBF, IObufsiz);

  return fp;
}



/* write XML grid description into .xmf file with file pointer fp */
void write_xdmf_xmf(FILE *fp, long voffset, long xyzoffset,
                    char *vname, const char *suffix,
		    char *nodename, double time,
		    int n[3], int bin, int dbl)
{
  char fname[1000];
  char fname_xyz[1000];
  int np = n[0] * n[1] * n[2];
  const char *format = (bin) ? "Binary" : "XML";

  /* filenames for field and also xyz data*/
  snprintf(fname, 1000, "%s.%s.bin", vname, suffix);
  snprintf(fname_xyz, 1000, "xyz.%s.bin", suffix);

  /* write <Grid> info to .xmf file in normal xdmf format */
  if(B_E_grid==B_E_grid_xmf)
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
  /* write <Grid> info to .xmf file in simplified format */
  else if(B_E_grid==B_E_grid_smf)
    fprintf(fp, B_E_grid, nodename, n[0], n[1], n[2], xyzoffset, voffset);
  else
    errorexit("B_E_grid is not set correctly");
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
void write_plane_xdmf(tVarList *vl, int norm, const char *outdir,
                      double Time)
{
  tMesh *mesh = vl->mesh;
  int bin = 1; /* we can only do binary output right now */
  int dbl = 0; /* we output float not double */
  int vbytes = ((dbl) ? sizeof(double) : sizeof(float));
  long voffset, xyzoffset;
  long xyzoffset0=0, xyztotal, xyzcnt, vtotal, vcnt, Vtotal;
  FILE *fpxmf, *fpbin, *fpxyz=NULL;
  char ndname[100];
  int ix = Ind( Gets(Par("output_xcoord")) );
  const char *suffix[] = { "yz", "xz", "xy" };
  double X0[] = { Getd(Par("outputX0")),
                  Getd(Par("outputY0")),
                  Getd(Par("outputZ0")) };
  int vli;
  int bufsize  = Geti(Par("fwrite_bufsize"));
  int syncmode = get_file_sync_mode(mesh); /* uses Par("file_sync") */
  char *bufxmf = cmalloc(bufsize); /* larger buffers for write */
  char *bufbin = cmalloc(bufsize);
  char *bufxyz = cmalloc(bufsize);
  tOutpars outpars[1] = {0};

  /* write region info into outpars */
  outpars->outputregion = Par("2doutputregion");
  output_set_regions_in_outpars(mesh, outpars);

  /* loop over varlist */
  for(vli=0; vli<vl->n; vli++)
  {
    int vi = vl->index[vli];
    char *vname = VarName(vi);
    int rk;
    int size = nMPI_size();

    /* num of bytes written for var so far */
    vtotal = 0;

    /* MPI motivated loop to assign work */
    for(rk=0; rk<size; rk++)
    {
      /* do work when it is my turn */
      if(rk == nMPI_rank())
      {
        if(Rank0) /* open xmf to start a new spatial series */
          fpxmf = fopen_xdmf_xmf(vname, outdir, suffix[norm], Time,
                                 bufxmf,bufsize);
        else /* just add to the same spatial series */
          fpxmf = fopen_add_spatial_xdmf_xmf(vname, outdir, suffix[norm],
                                             bufxmf,bufsize);

        /* open binary files */
        fpbin = fopen_bin(vname, outdir, suffix[norm], bufbin,bufsize);
        xyztotal = 0; /* num of bytes written for xyz */
        if(vli==0) /* write xyz only for first var in list */
        {
          fpxyz = fopen_bin("xyz", outdir, suffix[norm], bufxyz,bufsize);
          xyzoffset0 = ftell(fpxyz);
        }

        /* loop over all leaf nodes */
        formylnodes_noomp(mesh)
        {
          tNode *node = MyLnode;

          /* do something only if this proc has dat */
          if(node->dat)
          {
            int n[] = { node->n[0], node->n[1], node->n[2] };
            int plane[3];
            intList *plist;
            int normal;

            /* Check if we want to output this node. Skip node if it
               is outside the output region */
            if( !output_keep_elm(node, outpars) ) continue;

            /* find indices of nearest, if all are negative, node does not have
               outputX0, outputY0, outputZ0 */
            nearest_ijk_of_XYZplus(node, plane, X0);
            normal = approxXYZnormal_of_xyznormal(node, norm);

            if(normal>=0)
            if(plane[normal]>=0)
            {
              /* node name and n for plane */
              nodename(node, ndname,99);
              n[normal] = 1;

              /* list of points in plane */
              plist = pointindexList_plane(node, normal, plane);

              /* write point's x,y,z coordinates */
              if(vli==0) /* write xyz only for first var in list */
              {
                double *px = Vard(node, ix);
                double *py = Vard(node, ix+1);
                double *pz = Vard(node, ix+2);
                xyzoffset = ftell(fpxyz);
                xyzcnt = write_3buffers_idx(px,py,pz, plist, dbl, fpxyz);
              }
              else
              {
                xyzoffset = xyzoffset0 + xyztotal;
                /* we wrote 3 things (x,y,z) for vli=0 */
                xyzcnt = 3*(plist->n); //number of items written
              }
              xyzcnt *= vbytes;   //number of bytes written
              xyztotal += xyzcnt; //byte total outputted for xyz so far

              if(node->dat->v[vi])
              {
                /* write binary data in var */
                voffset = ftell(fpbin);
                vcnt = write_buffer_idx(Vard(node,vi), plist, dbl, fpbin);
                vcnt *= vbytes; //number of bytes written
                vtotal += vcnt; //byte total outputted for var so far

                /* write grid information into xmf file */
                write_xdmf_xmf(fpxmf, voffset, xyzoffset, vname, suffix[norm],
                               ndname, Time, n, bin, dbl);
              }
              intList_free(plist);
            }
          }
        }
        /* close files on this proc now */
        if(vli==0) fclose_sync_mode(fpxyz, syncmode); /* done only for first var in list */
        fclose_sync_mode(fpbin, syncmode);
        fclose_xdmf_xmf(fpxmf, syncmode, rk==size-1); /* last rank puts end markers */
        fs_sync(mesh); /* make sure every MPI proc flushes buffers to disk */
      }
      /* wait until everyone is here */
      MCK( nMPI_barrier() );
    } /* end rk-loop */

    /* get number of bytes written for binary data of vname */
    Vtotal = vtotal;
    MCK( nMPI_Allreduce(&vtotal, &Vtotal, 1, nMPI_LONG, nMPI_SUM) );
    //PRF;printf(": Vtotal=%ld\n", Vtotal);
    /* if we have written no binary data for vname at all, the last rank
       removes the entire xmf entry from the file */
    if( (Vtotal==0) && (nMPI_rank()==size-1) )
      rm_empty_spatial_xdmf_xmf(vname, outdir, suffix[norm], Time, syncmode);
  }
  free(bufxyz);
  free(bufbin);
  free(bufxmf);
}


/* output varlist in XDMF form for the entire volume */
void output3d_xdmf(tVarList *vl, int It, double Time)
{
  tMesh *mesh = vl->mesh;
  int outd = Par("outdir");
  char *outdir = Gets(outd);
  int bin = 1; /* we can only do binary output right now */
  int dbl = 0; /* we output float not double */
  int vbytes = ((dbl) ? sizeof(double) : sizeof(float));
  long voffset, xyzoffset;
  long xyzoffset0=0, xyztotal, xyzcnt, vtotal, vcnt, Vtotal;
  FILE *fpxmf, *fpbin, *fpxyz=NULL;
  char ndname[100];
  int ix = Ind( Gets(Par("output_xcoord")) );
  const char *suffix = "xyz";
  int vli;
  int bufsize  = Geti(Par("fwrite_bufsize"));
  int syncmode = get_file_sync_mode(mesh); /* uses Par("file_sync") */
  char *bufxmf = cmalloc(bufsize); /* larger buffers for write */
  char *bufbin = cmalloc(bufsize);
  char *bufxyz = cmalloc(bufsize);
  tOutpars outpars[1] = {0};

  /* write region info into outpars */
  outpars->outputregion = Par("3doutputregion");
  output_set_regions_in_outpars(mesh, outpars);

  /* loop over varlist */
  for(vli=0; vli<vl->n; vli++)
  {
    int vi = vl->index[vli];
    char *vname = VarName(vi);
    int rk;
    int size = nMPI_size();

    /* num of bytes written for var so far */
    vtotal = 0;

    /* MPI motivated loop to assign work */
    for(rk=0; rk<size; rk++)
    {
      /* do work when it is my turn */
      if(rk == nMPI_rank())
      {
        if(Rank0) /* open xmf to start a new spatial series */
          fpxmf = fopen_xdmf_xmf(vname, outdir, suffix, Time, bufxmf,bufsize);
        else /* just add to the same spatial series */
          fpxmf = fopen_add_spatial_xdmf_xmf(vname, outdir, suffix,
                                             bufxmf,bufsize);
        /* open binary files */
        fpbin = fopen_bin(vname, outdir, suffix, bufbin,bufsize);
        xyztotal = 0; /* num of bytes written for xyz */
        if(vli==0) /* write xyz only for first var in list */
        {
          fpxyz = fopen_bin("xyz",outdir,suffix, bufxyz,bufsize);
          xyzoffset0 = ftell(fpxyz);
        }

        /* loop over all leaf nodes */
        formylnodes_noomp(mesh)
        {
          tNode *node = MyLnode;

          /* do something only if this proc has dat */
          if(node->dat)
          {
            int np = node->np;
            int n[] = { node->n[0], node->n[1], node->n[2] };

            /* Check if we want to output this node. Skip node if it
               is outside the output region */
            if( !output_keep_elm(node, outpars) ) continue;

            /* write point's x,y,z coordinates */
            if(vli==0) /* write xyz only for first var in list */
            {
              double *px = Vard(node, ix);
              double *py = Vard(node, ix+1);
              double *pz = Vard(node, ix+2);
              xyzoffset = ftell(fpxyz);
              xyzcnt = write_3buffers(px,py,pz, np, dbl, fpxyz);
            }
            else
            {
              xyzoffset = xyzoffset0 + xyztotal;
              /* we wrote 3 things (x,y,z) for vli=0 */
              xyzcnt = 3*np;      //number of items written
            }
            xyzcnt *= vbytes;   //number of bytes written
            xyztotal += xyzcnt; //byte total outputted for xyz so far

            if(node->dat->v[vi])
            {
              /* get node name */
              nodename(node, ndname,99);

              /* write binary data in var */
              voffset = ftell(fpbin);
              vcnt = write_buffer(Vard(node,vi), np, dbl, fpbin);
              vcnt *= vbytes; //number of bytes written
              vtotal += vcnt; //byte total outputted for var so far

              /* write grid information into xmf file */
              write_xdmf_xmf(fpxmf, voffset, xyzoffset, vname, suffix,
                             ndname, Time, n, bin, dbl);
            }
          }
        }
        /* close files on this proc now */
        if(vli==0) fclose_sync_mode(fpxyz, syncmode); /* done only for first var in list */
        fclose_sync_mode(fpbin, syncmode);
        fclose_xdmf_xmf(fpxmf, syncmode, rk==size-1); /* last rank puts end markers */
        fs_sync(mesh); /* make sure every MPI proc flushes buffers to disk */
      }
      /* wait until everyone is here */
      MCK( nMPI_barrier() );
    } /* end rk-loop */

    /* get number of bytes written for binary data of vname */
    Vtotal = vtotal;
    MCK( nMPI_Allreduce(&vtotal, &Vtotal, 1, nMPI_LONG, nMPI_SUM) );
    //PRF;printf(": Vtotal=%ld\n", Vtotal);
    /* if we have written no binary data for vname at all, the last rank
       removes the entire xmf entry from the file */
    if( (Vtotal==0) && (nMPI_rank()==size-1) )
      rm_empty_spatial_xdmf_xmf(vname, outdir, suffix, Time, syncmode);
  }
  free(bufxyz);
  free(bufbin);
  free(bufxmf);
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
  return write_3buffers_Fwrite(b1,b2,b3,buflen, dbl, fp, fwrite);
}

/* write 3 buffers as double or float
   using Fwrite = fwrite, fwrite_big or fwrite_little */
size_t write_3buffers_Fwrite(const double *b1, const double *b2,
                             const double *b3, int buflen, int dbl, FILE *fp,
                             size_t (*Fwrite)(const void *ptr, size_t size,
                                              size_t nmemb, FILE *fp))
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
      cnt += Fwrite(&dval, sizeof(double), 1, fp);
      dval = b2[i];
      cnt += Fwrite(&dval, sizeof(double), 1, fp);
      dval = b3[i];
      cnt += Fwrite(&dval, sizeof(double), 1, fp);
    }
    else
    {
      float fval;
      fval = b1[i];
      cnt += Fwrite(&fval, sizeof(float), 1, fp);
      fval = b2[i];
      cnt += Fwrite(&fval, sizeof(float), 1, fp);
      fval = b3[i];
      cnt += Fwrite(&fval, sizeof(float), 1, fp);
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
