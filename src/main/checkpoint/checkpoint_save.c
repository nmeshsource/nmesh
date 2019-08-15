/* checkpoint.c */
/* Wolfgang Tichy, 8/2019 */


#include "nmesh.h"
#include "checkpoint.h"

/******************************************************************/
/* some functions to save nmesh data for checkpoints  */
/******************************************************************/

int checkpoint_save(tMesh *mesh)
{
  char *outdir = Gets(Par("outdir"));

  int cl = strlen(outdir) + 20;
  char *cdir = cmalloc(cl);

  int pl = cl + 40;
  char *pars  = cmalloc(pl);
  char *pats  = cmalloc(pl);
  char *nodes = cmalloc(pl);
  char *vars  = cmalloc(pl);

  /* output filenames */
  snprintf(cdir,cl, "%s/%s", outdir, "checkpoint");
  snprintf(pars,pl,  "%s/%s", cdir, "save_pars.txt");
  snprintf(pats,pl,  "%s/%s", cdir, "patches.txt");
  snprintf(nodes,pl, "%s/%s", cdir, "nodes.txt");
  snprintf(vars,pl,  "%s/%s", cdir, "variables.bin");

  /* save checkpoint in various files */
  checkpoint_save_patches(mesh, pats);
  checkpoint_save_nodes(mesh, nodes);
  checkpoint_save_EvoVars(mesh, vars);

  /* free strings */
  free(vars);
  free(nodes);
  free(pats);
  free(pars);
  free(cdir);
}


/******************************************************************/
/* functions to save patches */
/******************************************************************/
/* save patch info */
int checkpoint_save_patches(tMesh *mesh, char *fname)
{
  FILE *fp;
  int p;

  /* open destination file */
  fp = fopen(fname, "wb");
  if(!fp) errorexits("failed opening %s", fname);

  /* header with some mesh info */
  fprintf(fp, "mesh->\n");
  fprintf(fp, " time = %.19g\n", mesh->time);
  fprintf(fp, " iteration = %d\n", mesh->iteration);
  fprintf(fp, " dt = %.19g\n", mesh->dt);
  fprintf(fp, "\n");

  /* write data of each patch */
  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];

    checkpoint_write_pat(fp, pat)
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}

/* write non-pointer part of tPat */
void checkpoint_write_pat(FILE *fp, tPat *pat)
{
  int d, f;

  fprintf(fp, "patch%d: pat->\n", pat->p);

  for(f=0; f<6; f++)
    fprintf(fp, " bbox[%d] = %.19g\n", f, pat->bbox[f]);

  fprintf(fp,   " nmax = %d\n", pat->nmax);

  for(d=0; d<3; d++)
    fprintf(fp, " periodic[%d] = %d\n", d, pat->periodic[d]);

  for(d=0; d<3; d++)
    fprintf(fp, " rnode->n[%d] = %d\n", d, rnode->n[d]);

  checkpoint_write_CI(fp, pat->CI);
}


/* write non-pointer part of tCoordInfo */
void checkpoint_write_CI(FILE *fp, tCoordInfo *CI)
{
  int d, f;

  fprintf(fp, " CI->\n");

  for(f=0; f<6; f++)
    fprintf(fp,   "  iSurf[%d] = %d\n", f, CI->iSurf[f]);

  for(f=0; f<6; f++)
    for(d=0; d<3; d++)
      fprintf(fp, "  idSurfdX[%d][%d] = %d\n", f,d, CI->idSurfdX[f][d]);

  for(f=0; f<6; f++)
    fprintf(fp,   "  s[%d] = %.19g\n", f, CI->s[f]);

  for(d=0; d<3; d++)
    fprintf(fp,   "  xc[%d] = %.19g\n", d, CI->xc[d]);

  fprintf(fp,     "  dom = %d\n", CI->dom);
  fprintf(fp,     "  type = %d\n", CI->type);
}


/******************************************************************/
/* functions to save nodes */
/******************************************************************/
/* save node info */
int checkpoint_save_nodes(tMesh *mesh, char *fname)
{
  FILE *fp;
  tNlist *rnlist, *el;
  int p;

  /* make list of all root nodes, and save first elem. in rnlist */
  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];
    tNode *rnode = pat->rnode;

    if(p==0) rnlist = el = alloc_nodelist(rnode);
    else     el = addnode_to_nodelist_after(el, rnode);
  }

  /* open destination file */
  fp = fopen(fname, "wb");
  if(!fp) errorexits("failed opening %s", fname);

  /* write all nodes */
  checkpoint_write_nodetrees(fp, rnlist);
  fprintf(fp, "\n");

  fclose(fp);
  free_nodelist(rnlist);
  return 0;
}

/* write nodelist and all their children */
void checkpoint_write_nodetrees(FILE *fp, tNlist *rnlist)
{
  tNlist *nlist  = rnlist;
  tNlist *cnlist = NULL;
  int d;

  fprintf(fp, "nodetrees in mesh:\n");

  /* write all nodes in rnlist anf their children */
  while(nlist)
  {
    checkpoint_write_nodes_with_child0(fp, nlist);

    /* make list of all children of nlist, and then update nlist */
    free_nodelist(cnlist);
    cnlist = childnodelist_of_nodelist(nlist);
    nlist = cnlist;
  }
}

/* write non-pointer part of tPat */
void checkpoint_write_1_nodetree(FILE *fp, tNode *rnode)
{
  tNlist *nlist, *cnlist;
  int d;

  fprintf(fp,   "nodetree in patch%d: rnode->\n");
  fprintf(fp,   " pat->p = %.19g\n", rnode->pat->p);

  /* put root node into a nodelist */
  nlist = alloc_nodelist(rnode);

  /* write all nodes in nlist */
  while(nlist)
  {
    checkpoint_write_nodes_with_child0(fp, nlist);

    /* make list of all children of nlist, and then update nlist */
    cnlist = childnodelist_of_nodelist(nlist);
    free_nodelist(nlist);
    nlist = cnlist;
  }
}

/* write info about all nodes in a list */
void checkpoint_write_nodes_with_child0(FILE *fp, tNlist *nlist)
{
  tNlist *elem;

  fornodelist(nlist, elem)
  {
    tNode *node = elem->node;
    tNode *child0 = node->child[0];

    if(child0)
    {
      checkpoint_write_node(fp, node);
      fprintf(fp, "\n");
    }
  }
}

/* write info about one node */
void checkpoint_write_node(FILE *fp, tNode *node)
{
  tNode *child0 = node->child[0];
  char name[256];
  int d;

  nodename(node, name,255);

  /* info in node struct */
  fprintf(fp, "nodename = %s\n", name);

  // Not needed:
  /*
  fprintf(fp,   " dt = %.19g\n", node->dt);
  fprintf(fp,   " time = %.19g\n", node->time);

  for(d=0; d<3; d++)
    fprintf(fp, " n[%d] = %d\n", d, node->n[d]);

  fprintf(fp,   " ijk = %d\n", node->ijk);
  fprintf(fp,   " l = %d\n", node->l);
  fprintf(fp,   " leaf = %d\n", node->leaf);
  */

  if(child0)
  {
    /* add info about child0, so that we can easily re-create children */
    for(d=0; d<3; d++)
      fprintf(fp, " child[0]->n[%d] = %d\n", d, child0->n[d]);
  }
}


/******************************************************************/
/* functions to save variables */
/******************************************************************/

/* save all EvoVars */
int checkpoint_save_EvoVars(tMesh *mesh, char *fname)
{
  tVarList *vl = vlalloc(mesh);
  FILE *fp;

  /* loop over all vars and put anything that not AUXVAR into vl */
  for(vi=0; vi<nvdb; vi++)
    if(MeshVarType(mesh, vi) != AUXVAR) vlpushone(vl, vi)

  /* open destination file */
  fp = fopen(fname, "wb");
  if(!fp) errorexits("failed opening %s", fname);

  /* write var list vl in little endian format */
  checkpoint_write_vl(fp, vl, 0);

  fclose(fp);
  vlfree(vl);
  return 0;
}

/* output varlist on each node */
void checkpoint_write_vl(FILE *fp, tVarList *vl, int write_big)
{
  tMesh *mesh = vl->mesh;
  char name[256];
  int vli, rk;

  /* write all var names */
  if(Rank0)
  {
    fprintf(fp, "variable list in this file:\n");
    fprintf(fp, "vl->n = %d\n", vl->n);
    forvl(vl, vli)
    {
      int vi = Vind(vl, vli);
      fprintf(fp, "%s\n", VarName(vi));
    }
    fprintf(fp, "\n");
    fprintf(fp, "node data:\n");
  }

  /* MPI motivated loop to assign work */
  for(rk=0; rk<nMPI_size(); rk++)
  {
    /* do work when it is my turn */
    if(rk == nMPI_rank())
    {
      /* loop over all leaf nodes */
      formylnodes_noomp(mesh)
      {
        tNode *node = MyLnode;

        /* node name and number of points np on this node */
        nodename(node, name,255);
        fprintf(fp, "node\n");
        fprintf(fp, "%s\n", name);
        fprintf(fp, "%d\n", node->np);

        /* write data for all in vl on this node */
        forvl(vl, vli)
        {
          /* do something only if this proc has dat */
          if(node->dat)
          {
            int vi = Vind(vl, vli);
            tArray *va = node->dat->v[vi]
            if(va)
            {
              /* print only index in varlist */
              fprintf(fp, "%d\n", vli);

              /* write var array in raw binary */
              if(write_big)
                fwrite_big(Arrd_(va), sizeof(double), node->np, fp);
              else
                fwrite_little(Arrd_(va), sizeof(double), node->np, fp);
              fprintf(fp, "\n");
            }
          }
        }
        fprintf(fp, "\n");
      } /* end node-loop */
    }
    /* wait until everyone is here */
    nMPI_barrier();
  } /* end rk-loop */
}
