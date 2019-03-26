/* XDMF_out.c */
/* Wolfgang Tichy 3/2019 */

#include "nmesh.h"
#include "output.h"



/* XML format strings to make .xmf files using printf,
   based on bamps and Paraview */
char *B_head = 
  "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
  "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.1\">\n"
  "  <Domain>\n";
char *E_head =
  "  </Domain>\n"
  "</Xdmf>\n";

char *B_temporal = 
  "    <Grid CollectionType=\"Temporal\" GridType=\"Collection\" Name=\"Collection\">\n"
  "      <Geometry Type=\"None\"/>\n"
  "      <Topology Dimensions=\"0\" Type=\"NoTopology\"/>\n";
char *E_temporal = 
  "    </Grid>\n";

char *B_spatial =
  "      <Grid CollectionType=\"None\" GridType=\"Collection\" Name=\"Collection\">\n"
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
FILE *fopen_xdmf_xmf(char *varname, char *outdir, char *suffix)
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

  /* we want to append more data, which requires us overwrite to remove
     the last E_temporal and E_head */
  offset = strlen(E_temporal) + strlen(E_head);
  fseek(fp, -offset, SEEK_END);

  /* start new collection of nodes */
  fprintf(fp, "%s", B_spatial);

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


/* output in XDMF in one plane */
void write_plane_xdmf(tNode *node, int normal, int plane[], tArray *va,
                      char *vname, char *outdir, char *suffix, double Time,
                      int write_xyz)
{
  tMesh *mesh = node->pat->mesh;
  intList *plist = pointindexList_plane(node, normal, plane);
  int bin = 1; /* we can only do binary output right now */
  int dbl = 0; /* we output float not double */
  int voffset, xyzoffset;
  FILE *fpxmf, *fpbin;
  char ndname[100];
  int n[] = { node->n[0], node->n[1], node->n[2] };

  /* open xmf and bin files */
  fpxmf = fopen_xdmf_xmf(vname, outdir, suffix);
  fpbin = fopen_bin(vname, outdir, suffix);

  /* write data in array va */
  voffset = ftell(fpbin);
  write_buffer_idx(Arrd(va), plist, dbl, fpbin);  

  /* node name and n for plane */
  nodename(node, ndname,99);
  n[normal] = 1;

  /* write point's x,y,z coordinates */
  if(write_xyz)
  {
    int ix = Ind("x");
    double *px = Vard(node, ix);
    double *py = Vard(node, ix+1);
    double *pz = Vard(node, ix+2);
    FILE *fpxyz = fopen_bin("xyz", outdir, suffix);

    write_3buffers_idx(px,py,pz, plist, dbl, fpxyz);
    fclose(fpxyz);
  }

  /* we write 3 things (x,y,z) for each var */
  xyzoffset = 3 * voffset;

  /* write grid information into xmf file */
  write_xdmf_xmf(fpxmf, voffset, xyzoffset, vname, suffix, ndname, Time,
                 n, bin, dbl);

  /* close files and free list */
  fclose(fpbin);
  fclose_xdmf_xmf(fpxmf);
  free_intList(plist);
}


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
