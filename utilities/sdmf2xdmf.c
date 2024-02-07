/* Parse output of diff and check for differences in floating point numbers */
/* (c) Wolfgang Tichy 2002 */

#define STRLEN 65536

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../src/utility/output/XDMF_out.h"


int main(int argc, char *argv[])
{
  char *textname;
  char xmfname[STRLEN];
  FILE *ftxt;
  FILE *fxmf;
  int slen;
  char line[STRLEN];
  char *pstr;
  char str[STRLEN];
  char str1[STRLEN];
  char str2[STRLEN];
  char vname[1000];
  char suffix[1000];
  char fname[STRLEN];
  char fname_xyz[STRLEN];
  char TopologyType[STRLEN];
  char AttributeCenter[STRLEN];
  int first_time;

  /* read args */
  if(argc!=2)
  {
    printf("Usage: sdmf2xdmf file.txt\n");
    return -1;
  }
  textname = argv[1];

  /* get vname, suffix */
  strcpy(vname, textname);
  pstr = strstr(vname, ".");
  pstr[0] = 0;
  pstr++;
  strcpy(suffix, pstr);
  pstr = strstr(suffix, ".");
  pstr[0] = 0;
  pstr++;
  if(strcmp(pstr, "txt"))
  {
    printf("file extension of %s is not txt\n", textname);
    return -1;
  }

  /* open textname */
  ftxt = fopen(textname, "r");
  if(!ftxt)
  {
    printf("could not open %s\n", textname);
    return -1;
  }
  /* read header in textname */
  fgets(line,STRLEN, ftxt);
  if(strncmp(line, B_head_smf, 4))
  {
    printf("sdmf header missing in %s\n", textname);
    return -1;
  }
  pstr = strstr(line, "TopologyType:");
  sscanf(pstr,"%s %s", str1, TopologyType);
  pstr = strstr(line, "AttributeCenter:");
  sscanf(pstr,"%s %s", str1, AttributeCenter);
  fgets(line,STRLEN, ftxt);

  /* construct xmfname */
  strcpy(xmfname, textname);
  slen = strlen(xmfname);
  xmfname[slen-3] = 0;
  strcat(xmfname, "xmf");

  /* open xmfname */
  fxmf = fopen(xmfname, "w");
  if(!ftxt)
  {
    printf("could not open %s\n", xmfname);
    return -1;
  }
  /* write header */
  fprintf(fxmf, "%s", B_head_xmf);
  fprintf(fxmf, "%s", B_temporal_xmf);


  /* filenames for field and also xyz data */
  snprintf(fname, STRLEN, "%s.%s.bin", vname, suffix);
  snprintf(fname_xyz, STRLEN, "xyz.%s.bin", suffix);

  /* now go over textname */
  first_time = 1;
  while(fgets(line,STRLEN, ftxt))
  {
    char nodename[STRLEN];
    double time;
    int np;
    int n[3];
    int bin = 1;
    const char *format = (bin) ? "Binary" : "XML";
    long xyzoffset, voffset;

    if(line[0]=='\n')
    {
      if(first_time)
        first_time = 0;
      else /* close previous time */
        fprintf(fxmf, E_spatial_xmf);

      /* new time series starts, so get its time */
      fgets(line,STRLEN, ftxt);
      sscanf(line,"%s %s %s %lf", str, str1, str2, &time);
      fprintf(fxmf, B_spatial_xmf, time);
      fgets(line,STRLEN, ftxt);
    }

    /* read line with node info */
    sscanf(line,"%s  %d %d %d  %ld %ld",
           nodename, &n[0],&n[1],&n[2], &xyzoffset, &voffset);
    np = n[0]*n[1]*n[2];

    /* print node info, as in write_xdmf_xmf of src/utility/output/XDMF_out.c */
    fprintf(fxmf, B_E_grid_xmf,
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

  /* close previous time */
  fprintf(fxmf, E_spatial_xmf);
  fprintf(fxmf, E_temporal_xmf);

  /* write end marker */
  fprintf(fxmf, "%s", E_head_xmf);

  fclose(fxmf);
  fclose(ftxt);

  return 0;
}
