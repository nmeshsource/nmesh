/* output0d.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "output.h"



/* 0d output */
void output0d_vl(tVarList *vl, int It, double T)
{
  tMesh *mesh = vl->mesh;
  int pi;

  TIMER_START;

  /* 0d output for mesh */
  output0d_mesh_vl(vl,NULL, It, T);

  /* 0d output for each patch separately */
  forpatches(mesh, pi)
    output0d_mesh_vl(vl, mesh->pat[pi], It, T);

  TIMER_STOP;
}


/* 0d output for one patch or entire mesh if pat==NULL */
void output0d_mesh_vl(tVarList *vl, tPat *pat, int It, double T)
{
  tMesh *mesh = vl->mesh;
  char filename[1000];
  double max, min, maxAbs, rms, mean, VolInt;
  double Vol;
  int vli;

  for(vli=0; vli<vl->n; vli++)
  {
    int vi = vl->index[vli];
    char *name = VarName(vi);

    /* get volume only once */
    if(vli==0) Vol = MeshVolumeIntegral(mesh,pat, vi, 0.,0);

    VolInt = MeshVolumeIntegral(mesh,pat, vi, 1.,0);
    mean   = VolInt/Vol;
    rms    = sqrt(MeshVolumeIntegral(mesh,pat, vi, 2.,0) / Vol);

    min    = MeshMin(mesh,pat, vi);
    max    = MeshMax(mesh,pat, vi);
    maxAbs = max2(fabs(min), fabs(max));

    if(Rank0)
    {
      /* output max, min, maxAbs, rms, mean, VolInt */
      output0d_filename(mesh, filename,999, name, "VolInt", pat);
      output0d_value(filename, T, VolInt);

      output0d_filename(mesh, filename,999, name, "mean", pat);
      output0d_value(filename, T, mean);

      output0d_filename(mesh, filename,999, name, "rms", pat);
      output0d_value(filename, T, rms);

      output0d_filename(mesh, filename,999, name, "min", pat);
      output0d_value(filename, T, min);

      output0d_filename(mesh, filename,999, name, "max", pat);
      output0d_value(filename, T, max);

      output0d_filename(mesh, filename,999, name, "maxAbs", pat);
      output0d_value(filename, T, maxAbs);
    }
  }
}

/* filename for 0d output */
void output0d_filename(tMesh *mesh, char *filename, int len,
                       char *name, char *type, tPat *pat)
{
  char *outdir = Gets(Par("outdir"));

  if(pat)
  {
    char fmt[100];
    snprintf(fmt,99, "%%s/%%s_%%s.%%0%ddt", (int) log10(mesh->npats)+1);
    snprintf(filename,len, fmt, outdir, name, type, pat->p);
  }
  else
  {
    snprintf(filename,len, "%s/%s_%s.t", outdir, name, type);
  }
}

/* output one value */
void output0d_value(char *filename, double time, double val)
{
  FILE *fp;
  
  /* open file */
  fp = fopen(filename, "a");
  if(!fp) errorexits("failed opening %s", filename);

  /* write value */
  fprintf(fp, "%.15g %.15g\n", time, val);

  /* close file */
  fclose(fp);
}
