/* checkpoint_load.c */
/* Wolfgang Tichy, 8/2019 */


#include "nmesh.h"
#include "checkpoint.h"

/******************************************************************/
/* some functions to load nmesh data for checkpoints  */
/******************************************************************/

/******************************************************************/
/* functions to load patches */
/******************************************************************/
/* load patch info */
int checkpoint_load_patches(tMesh *mesh, char *fname)
{
  FILE *fp;
  char buf[1000];
  tPat *pat = NULL;
  char par[1000], val[1000];
  int clabel=0, nmax=0, useF=0;
  int n[] = { 0,0,0 };
  double bbox[] = { 0.,0.,0.,0.,0.,0. };

  /* open source file */
  fp = fopen(fname, "rb");
  if(!fp) errorexits("failed opening %s", fname);

  /* read file line by line */
  while(fgets(buf,999, fp))
  {
    /* CI-> signifies that we know all to create a new patch */
    if(strcmp(buf, " CI->\n")==0)
    {
      pat = add_patch(mesh, bbox, n, nmax);
      useF = 0;
    }

    /* read various info pieces */
    if(get_par_from_str(buf, par, val, 999))
    {
      if(strcmp(par, "time")==0)      mesh->time = atof(val);
      if(strcmp(par, "iteration")==0) mesh->iteration = atoi(val);
      if(strcmp(par, "dt")==0)        mesh->dt = atof(val);

      if(strcmp(par, "bbox[0]")==0) bbox[0] = atof(val);
      if(strcmp(par, "bbox[1]")==0) bbox[1] = atof(val);
      if(strcmp(par, "bbox[2]")==0) bbox[2] = atof(val);
      if(strcmp(par, "bbox[3]")==0) bbox[3] = atof(val);
      if(strcmp(par, "bbox[4]")==0) bbox[4] = atof(val);
      if(strcmp(par, "bbox[5]")==0) bbox[5] = atof(val);
      if(strcmp(par, "nmax")==0) nmax = atoi(val);
      if(strcmp(par, "rnode->n[0]")==0) n[0] = atoi(val);
      if(strcmp(par, "rnode->n[1]")==0) n[1] = atoi(val);
      if(strcmp(par, "rnode->n[2]")==0) n[2] = atoi(val);
      if(strcmp(par, "s[0]")==0)  pat->CI->s[0] = atof(val);
      if(strcmp(par, "s[1]")==0)  pat->CI->s[1] = atof(val);
      if(strcmp(par, "s[2]")==0)  pat->CI->s[2] = atof(val);
      if(strcmp(par, "s[3]")==0)  pat->CI->s[3] = atof(val);
      if(strcmp(par, "s[4]")==0)  pat->CI->s[4] = atof(val);
      if(strcmp(par, "s[5]")==0)  pat->CI->s[5] = atof(val);
      if(strcmp(par, "xc[0]")==0) pat->CI->xc[0] = atof(val);
      if(strcmp(par, "xc[1]")==0) pat->CI->xc[1] = atof(val);
      if(strcmp(par, "xc[2]")==0) pat->CI->xc[2] = atof(val);
      if(strcmp(par, "dom")==0)   pat->CI->dom = atoi(val);
      if(strcmp(par, "type")==0)  pat->CI->type = atoi(val);
      if(strcmp(par, "use_FSurf")==0)  useF = atoi(val);
      if(strcmp(par, "coordinates_get_label(pat)")==0)
      {
        /* this signifies end of patch info, so now set final pat info */
        clabel = atoi(val);
        switch(clabel)
        {
        case CubedSphere:
          set_1_CubedSphere_pat(pat, 0, useF);
          break;
        case stretchedCubedSphere:
          set_1_CubedSphere_pat(pat, 1, useF);
          break;
        }
      }
    } /* end of if get_par_from_str */
  }

  fclose(fp);
  return 0;
}


/******************************************************************/
/* functions to load nodes */
/******************************************************************/
/* load node info */
int checkpoint_load_nodes(tMesh *mesh, char *fname)
{
  FILE *fp;
  char buf[1000];
  tPat *pat;
  tNode *parent;
  int n[3];

  /* open source file */
  fp = fopen(fname, "rb");
  if(!fp) errorexits("failed opening %s", fname);

  /* read file line by line */
  while(fgets(buf,999, fp))
  {
    /* all node names contain an '_' */
    if(strstr(buf, "_"))
    {
      //pat = ...;
      //parent = ...;

      /* read n for child0 */
      fscanf(fp, "%d", &(n[0]));
      fscanf(fp, "%d", &(n[1]));
      fscanf(fp, "%d", &(n[2]));

      /* make 8 child nodes */
      //...
    }
  }

  fclose(fp);
  return 0;
}


/******************************************************************/
/* functions to load variables */
/******************************************************************/

/* load all EvoVars */
int checkpoint_load_Vars(tMesh *mesh, char *fname)
{
  tVarList *vl;
  FILE *fp;

  /* open file */
  fp = fopen(fname, "rb");
  if(!fp) errorexits("failed opening %s", fname);

  /* get varlist from file */
  vl = checkpoint_make_vl(fp, mesh);

  /* now read data for vars */
  checkpoint_read_vl(fp, vl, 0);

  vlfree(vl);
  fclose(fp);
  return 0;
}

/* make varlist that we find in checkpoint file */
tVarList *checkpoint_make_vl(FILE *fp, tMesh *mesh)
{
  tVarList *vl = vlalloc(mesh);
  char buf[1000];
  int i, nvars;

  /* read first comment line */
  fgets(buf,999, fp);

  /* read number of vars */
  fscanf(fp, "%d", &nvars);

  /* read names and make varlist */
  for(i=0; i<nvars; i++)
  {
    fscanf(fp, "%s", buf);
    vlpushone(vl, Ind(buf));
  }

  return vl;
}

/* output varlist on each node */
void checkpoint_read_vl(FILE *fp, tVarList *vl, int read_big)
{
  //tMesh *mesh = vl->mesh;
  char buf[1000];
  int rk;

  /* MPI motivated loop to assign work */
  for(rk=0; rk<nMPI_size(); rk++)
  {
    /* do work when it is my turn */
    if(rk == nMPI_rank())
    {
      /* find string node */
      while(fgets(buf,999, fp))
      {
        if(strcmp(buf, "node\n"))
        {
          tNode *node;
          char name[256];
          int np, vli, vi;

          /* read var info */
          //fgets(name,255, fp);
          fscanf(fp, "%s", name);
          fscanf(fp, "%d", &np);
          fscanf(fp, "%d", &vli);
          vi = Vind(vl, vli);

          /* find node from its name */
          //...

          /* check if this var needs to be read on this node */
          if(node->dat)
          {
            double *v;

            /* switch on var on this node and get pointer to its data */
            enablevarcomp_innode(node, vi);
            v = Vard(node, vi);

            /* read var as raw binary */
            if(read_big) fread_big(v, sizeof(double), np, fp);
            else         fread_little(v, sizeof(double), np, fp);
          }
          else /* just over-read data */
          {
            long offset = sizeof(double) * np;
            fseek(fp, offset, SEEK_CUR);
          }
        } /* end if(strcmp(buf, "node\n")) */
      }
    }
    /* wait until everyone is here */
    nMPI_barrier();
  } /* end rk-loop */
}
