/* output0d.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "output.h"

extern char outpt[];

/* make output visible */
extern tOutput output[1];


/* 0d output */
void output0d_vl(tVarList *vl, int It, double T)
{
  if(vl)
  {
    tMesh *mesh = vl->mesh;
    int per_patch = Par("0doutput_per_patch");
    int per_patch_limit, patch_output;
    int pi;

    TIMER_START;

    /* 0d output for mesh */
    output0d_mesh_vl(vl,NULL, It, T);

    /* 0d output for each patch separately */
    per_patch_limit = Geti(per_patch); /* try to read number from par */
    if(per_patch_limit==0)
      patch_output = Getb(per_patch); /* if limit is zero, use yes/no */
    else
      patch_output = (mesh->npats <= per_patch_limit); /* use limit */
    if(patch_output)
      forpatches(mesh, pi)
        output0d_mesh_vl(vl, mesh->pat[pi], It, T);

    TIMER_STOP;
  }
}


/* 0d output for one patch or entire mesh if pat==NULL */
void output0d_mesh_vl(tVarList *vl, tPat *pat, int It, double T)
{
  tMesh *mesh = vl->mesh;
  char filename[1000];
  double max, min, maxAbs, rms, mean, VolInt;
  int add_xyz = Getb(Par("0doutput_add_xyz"));
  int interp_np = Geti(Par("output_interpolation_order"));
  int p, ijk, ipt;
  char nodeloc[105]; /* node location string inside patch */
  double X[3], xmin[3], xmax[3], *xmaxAbs;
  double Vol;
  int Npt; /* number of points at which we output */
  double val_pt[Noutptmax];
  double x_pt[Noutptmax][3], Xb_pt[Noutptmax][3];
  int have_pt[Noutptmax];
  tNode *node_pt[Noutptmax];
  int vli;

  /* set special points x_pt[], ... at which we output: */
  if(pat) /* we output only in one patch */
  {
    /* turn point output off, as a point is not in every patch */
    Npt = 0;
  }
  else /* we output on the entire mesh */
  {
    /* set number of points x_pt[] on which we output */
    Npt = output->Noutpt;
    if(Npt>Noutptmax) errorexit("output->Noutpt > (Noutptmax in tOutput)");
  }
  /* loop over special points */
  for(ipt=0; ipt<Npt; ipt++)
  {
    int d;
    for(d=0; d<3; d++) x_pt[ipt][d] = output->xpt[ipt][d];
    /* FIXME: use interp_var_xyz instead: */
    node_pt[ipt] = node_XYZ_of_xyz_mesh(mesh, X, x_pt[ipt]);
    /* NOTE: node_XYZ_of_xyz_mesh finds the 1st node that has x_pt.
             But what if several nodes have this point??? */
    //PRF;pr3v(": x_pt[ipt]", x_pt[ipt]);pr3v("X", X);
    //printf("node_pt[%d]=%p\n", ipt, node_pt[ipt]);
    if(node_pt[ipt])
      XbYbZb_of_XYZ(node_pt[ipt], Xb_pt[ipt], X);  /* set Xb in node */
  }

  /* loop over varlist vl */
  for(vli=0; vli<vl->n; vli++)
  {
    int vi = vl->index[vli];
    char *name = VarName(vi);

    /* get volume only once, and use var X to the power 0.0 to compute it */
    if(vli==0) Vol = MeshVolumeIntegral(mesh,pat, Ind("X"),0.0 ,0);

    VolInt = MeshVolumeIntegral(mesh,pat, vi, 1.,0);
    mean   = VolInt/Vol;
    rms    = sqrt( fabs(MeshVolumeIntegral(mesh,pat, vi, 2.,0) / Vol) );

    /* min, max and their positions xmin, xmax */
    min = MeshExtremumLoc(mesh,pat, vi, 0, &p, nodeloc, &ijk, X, xmin);
    max = MeshExtremumLoc(mesh,pat, vi, 1, &p, nodeloc, &ijk, X, xmax);

    /* if max<min there was no node with any data, then we skip the output */
    if(finit(min) && finit(max))
    {
      if(max<min)
      {
        if(!finit(rms))
        {
          PRF;printf(": skipping %s: "
                     "No elm with valid data found (max<min). "
                     "Probably all NaN (rms=%g).\n", name, rms);
          //printf("min=%g max=%g\n", min, max);
        }
        continue;
      }
    }

    /* maxAbs and its pos. */
    if(fabs(max)>fabs(min))
    {
      maxAbs  = fabs(max);
      xmaxAbs = &(xmax[0]);
    }
    else
    {
      maxAbs  = fabs(min);
      xmaxAbs = &(xmin[0]);
    }

    /* get value of var at some points */
    for(ipt=0; ipt<Npt; ipt++)
      have_pt[ipt] = interpolate_var_ok(node_pt[ipt], vi, Xb_pt[ipt],
                                 interp_np, INTERP_LAGRANGE, &(val_pt[ipt]));
                                  /* ^-- FIXME: use interp_var_xyz instead */
                                  /* NOTE: use INTERP_WENO if P_UNIFORM & fv */
    /* output is done by rank0 */
    if(Rank0)
    {
      /* output max, min, maxAbs, rms, mean, VolInt */
      output0d_filename(mesh, filename,999, name, "VolInt", pat);
      output0d_value(filename, T, VolInt, 0, NULL);

      output0d_filename(mesh, filename,999, name, "mean", pat);
      output0d_value(filename, T, mean, 0, NULL);

      output0d_filename(mesh, filename,999, name, "rms", pat);
      output0d_value(filename, T, rms, 0, NULL);

      output0d_filename(mesh, filename,999, name, "min", pat);
      output0d_value(filename, T, min, add_xyz, xmin);

      output0d_filename(mesh, filename,999, name, "max", pat);
      output0d_value(filename, T, max, add_xyz, xmax);

      output0d_filename(mesh, filename,999, name, "maxAbs", pat);
      output0d_value(filename, T, maxAbs, add_xyz, xmaxAbs);

      /* output value of var at some points */
      for(ipt=0; ipt<Npt; ipt++)
      {
        if(have_pt[ipt])
        {
          char typestr[99];
          //PRF;printf(": val_pt[ipt]=%g\n", val_pt[ipt]);
          snprintf(typestr,99, "%s%d", outpt, ipt);
          output0d_filename(mesh, filename,999, name, typestr, pat);
          output0d_value(filename, T, val_pt[ipt], 1, x_pt[ipt]);
        }
      }
    }
  } /* end loop over vli */
}

/* filename for 0d output */
void output0d_filename(tMesh *mesh, char *filename, int len,
                       char *name, const char *type, tPat *pat)
{
  char *outdir = Gets(Par("outdir"));

  if(pat) /* we output only in one patch */
  {
    char fmt[100];
    snprintf(fmt,99, "%%s/%%s_%%s.%%0%ddt", (int) log10(mesh->npats)+1);
    snprintf(filename,len, fmt, outdir, name, type, pat->p);
  }
  else /* we output on the entire mesh */
  {
    snprintf(filename,len, "%s/%s_%s.t", outdir, name, type);
  }
}

/* output one value */
void output0d_value(char *filename, double time, double val,
                    int coords, double x[3])
{
  FILE *fp;
  
  /* open file */
  fp = fopen(filename, "a");
  if(!fp) errorexits("failed opening %s", filename);

  /* write value */
  fprintf(fp, "%.15g %.15g", time, val);
  if(coords) fprintf(fp, " %.15g %.15g %.15g", x[0],x[1],x[2]);
  fprintf(fp, "\n");

  /* close file */
  fclose(fp);
}
