/* convert nmesh's sdmf output into xdmf output */
/* (c) Wolfgang Tichy 2/2024 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../src/utility/output/XDMF_out.h"

#define STRLEN 65536


int main(int argc, char *argv[])
{
  char *textname;
  char xmfname[STRLEN];
  FILE *ftxt;
  FILE *fxmf;
  int slen;
  char line[STRLEN];
  char *pstr, *p1, *p2;
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
    printf("Usage: sdmf2xdmf file.xyz.txt\n");
    return -1;
  }
  textname = argv[1];

  /* find location of two last dots in file name */
  strcpy(vname, textname);
  p1 = p2 = NULL;
  pstr = vname;
  while(pstr = strstr(pstr, "."))
  {
    p1 = p2;   /* pointer to 2nd to last . */
    p2 = pstr; /* pointer to last . */
    pstr++;    /* go to char after . */
  }
  if(!p2)
  {
    printf("file %s has no extension after suffix\n", textname);
    return -1;
  }
  if(!p1)
  {
    printf("file %s has no suffix\n", textname);
    return -1;
  }
  /* set vname, suffix */
  //printf("p1=%s\n", p1);
  //printf("p2=%s\n", p2);
  p1[0] = 0;
  p2[0] = 0;
  strcpy(suffix, p1+1);
  //printf("vname=%s\n", vname);
  //printf("suffix=%s\n", suffix);
  if(strcmp(p2+1, "txt"))
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

    /* skip empty lines */
    if(line[0]=='\n') continue;

    /* look for comments */
    if(line[0]=='#')
    {
      /* skip comments without "time */
      sscanf(line,"%s %s %s %lf", str, str1, str2, &time);
      if(strcmp(str1, "\"time")) continue;

      /* get next data line, but break if there is nothing */
      if(fgets(line,STRLEN, ftxt)==NULL) break;

      /* is this the 1st time we found the string "time ? */
      if(first_time)
        first_time = 0;
      else /* end previous time */
        fprintf(fxmf, E_spatial_xmf);

      /* start new time */
      fprintf(fxmf, B_spatial_xmf, time);
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
