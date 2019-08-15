/* checkpoint.c */
/* Wolfgang Tichy, 8/2019 */


#include "nmesh.h"
#include "checkpoint.h"

/******************************************************************/
/* some functions to save nmesh data for checkpoints  */
/******************************************************************/


/******************************************************************/
/* functions to save patches */
/******************************************************************/
/* save patch info */
int checkpoint_save_patches(tMesh *mesh, char *fname)
{
  int p;

  /* open destination file */
  out = fopen(fname, "wb");
  if(!out) errorexits("failed opening %s", fname);

  /* write data of each patch */
  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];

    checkpoint_write_pat(out, pat)
    fprintf(fp, "\n");
  }

  fclose(out);
  return 0;
}

/* write non-pointer part of tPat */
void checkpoint_write_pat(FILE *fp, tPat *pat)
{
  int d, f;

  fprintf(fp, "patch%d: pat->\n", pat->p);

  for(f=0; f<6; f++)
    fprintf(fp, " bbox[%d] = %.19g\n", f, pat->bbox[f]);

  fprintf(fp,   " bbdiag = %.19g\n", pat->bbdiag);
  fprintf(fp,   " nmax = %d\n", pat->nmax);

  for(d=0; d<3; d++)
    fprintf(fp, " periodic[%d] = %d\n", d, pat->periodic[d]);

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
  int p;

  /* open destination file */
  out = fopen(fname, "wb");
  if(!out) errorexits("failed opening %s", fname);

  /* write data of each patch */
  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];

    checkpoint_write_nodetree(out, pat->rnode)
    fprintf(fp, "\n");
  }

  fclose(out);
  return 0;
}

/* write non-pointer part of tPat */
void checkpoint_write_nodetree(FILE *fp, tNode *rnode)
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
    checkpoint_write_nodelist(fp, nlist);

    /* make list of all children of nlist, and then update nlist */
    cnlist = childnodelist_of_nodelist(nlist);
    free_nodelist(nlist);
    nlist = cnlist;
  }
}

/* write info about all nodes in a list */
void checkpoint_write_nodelist(FILE *fp, tNlist *nlist)
{
  tNlist *elem;

  fornodelist(nlist, elem)
  {
    checkpoint_write_node(fp, elem->node);
    fprintf(fp, "\n");
  }
}

/* write info about one node */
void checkpoint_write_node(FILE *fp, tNode *node)
{
  char loc[100];
  char ploc[100];
  int d;

  node_location_str(node, loc,99);
  node_location_str(node->parent, ploc,99);

  /* info in node struct */
  fprintf(fp, "node_location_str = %s\n", loc);
  fprintf(fp, "parent_location_str = %s\n", ploc);

  fprintf(fp,   " dt = %.19g\n", node->dt);
  fprintf(fp,   " time = %.19g\n", node->time);

  for(d=0; d<3; d++)
    fprintf(fp, " n[%d] = %d\n", d, node->n[d]);

  fprintf(fp,   " ijk = %d\n", node->ijk);
  fprintf(fp,   " l = %d\n", node->l);
  fprintf(fp,   " leaf = %d\n", node->leaf);

  /* add info about child0, so that we can easily create children */
  if(node->child[0])
    for(d=0; d<3; d++)
      fprintf(fp, " child[0]->n[%d] = %d\n", d, node->child[0]->n[d]);
}
