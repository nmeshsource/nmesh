/* storage.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"




/**************************************************************************/
/* basic memory management */
/**************************************************************************/

/**********************************************************************/
/* storage for arrays */
/**********************************************************************/
/* allocate an array */
tArray *alloc_array(int n[3]) //, void *Owner, int tOwner)
{
  int i;
  tArray *array = calloc(1, sizeof(tArray));
  if(!array) errorexit("out of memory");

  array->N = n[0] * n[1] * n[2];
  for(i=0; i<3; i++)  array->n[i] = n[i];

  array->a = calloc(array->N, sizeof(array->a[0]));
  if(!array->a) errorexit("out of memory for array->a");

  //array->Owner  = Owner;
  //array->tOwner = tOwner;

  return array;
}

/* free an array */
void free_array(tArray *array)
{
  free(array->a);
  free(array);
}


/**************************************************************************/
/* mesh, patch, and node storage */
/**************************************************************************/

/* allocate mesh */
tMesh *alloc_mesh(void)
{
  tMesh *mesh;

  mesh = calloc(1, sizeof(tMesh));
  if (!mesh) errorexit("out of memory");

  return mesh;
}

/* free mesh */
void free_mesh(tMesh *mesh)
{
  int i;

  if (!mesh) return;
  for(i = 0; i < mesh->npats; i++)
    free_patch(mesh->pat[i]);
  free(mesh);
}


/* allocate patch */
tPat *alloc_patch(tMesh *mesh, int p, int nD) 
{
  tPat *pat;

  pat = calloc(1, sizeof(tPat));
  if (!pat) errorexit("out of memory");

  pat->mesh = mesh;
  pat->p = p;
  pat->nD = nD;
    
  /* allocate storage for data pointers, they default to NULL */

  /* get mem. for diff. matrices */
  pat->D = calloc(nD, sizeof(pat->D[0]));
  if( !(pat->D) )
    errorexit("out of memory for diff. matrices");

  return pat;
} 

/* free pat, currently leaves mesh untouched */
void free_patch(tPat *pat) 
{
  int i;

  if (!pat) return;

  //for (i = 0; i < pat->mesh->nvariables; i++)
  //  disablevarcomp_inpat(pat, i);
  //free(pat->v);

  //free_all_bfaces(pat);

  free(pat);
}

tNode *alloc_node()
{
  tNode *node = calloc(1, sizeof(tNode));
  if(!node) errorexit("out of memory");

  return node;
}

/* free one node */
void free_node(tNode *node) 
{
  int i;

  if (!node) return;

  //for (i = 0; i < pat->mesh->nvariables; i++)
  //  disablevarcomp_inpat(pat, i);
  //free(pat->v);

  //free_all_bfaces(pat);

  free(node);
}


/**********************************************************************/
/* storage for dat lists in the nodes */
/**********************************************************************/
/* allocate for nv variables that can be enabled or disabled */
tDat *alloc_dat(int nv)
{
  tDat *dat;
  int j;

  dat = calloc(1, sizeof(tDat));
  if(!dat) errorexit("out of memory for dat");

  dat->nv = nv;

  dat->v = calloc(nv, sizeof(tArray *));
  if(!dat->v) errorexit("out of memory for dat->v");

  for(j=0; j<6; j++)
  {
    dat->g[j] = calloc(nv, sizeof(tArray *));
    if(!dat->g[j]) errorexit("out of memory for dat->g[j]");
  }
  return dat;
}
/* free dat and all arrays within it */
void free_dat(tDat *dat)
{
  int i,j;

  if(!dat) return;

  for(i=0; i<dat->nv; i++)
  {
    free_array(dat->v[i]);
    for(j=0; j<6; j++) free_array(dat->g[j][i]);
  }

  free(dat->v);
  for(j=0; j<6; j++) free(dat->g[j]);
  free(dat);
}

/* change dat->nv  to  dat->nv=nv_new */
void realloc_datvariables(tDat *dat, int nv_new)
{
  int i,j;

  dat->v = realloc(dat->v, nv_new*sizeof(tArray *));
  if(!dat->v) errorexit("out of memory for dat->v");

  for(j=0; j<6; j++)
  {
    dat->g[j] = realloc(dat->g, nv_new*sizeof(tArray *));
    if(!dat->g[j]) errorexit("out of memory for dat->g");
  }
  if(nv_new<dat->nv) errorexit("implement var removal");

  /* set newly added var pointers to NULL */
  for(i=dat->nv; i<nv_new; i++)
  {
    dat->v[i] = NULL;
    for(j=0; j<6; j++) dat->g[j][i] = NULL; 
  }
  dat->nv = nv_new;
}

/**********************************************************************/
/* storage for variable data base vdb in mesh */
/**********************************************************************/
void realloc_nodevariables(tNode *node, int nvdb_new)
{
  tDat *dat = node->dat;
  if(dat) realloc_datvariables(dat, nvdb_new);
}

void realloc_meshvariables(tMesh *mesh, int nvdb_new)
{
  tNode *node;
  if(1) printf("realloc_meshvariables from %d to %d\n", 
                 mesh->nvdb, nvdb_new);

  fornodelist(mesh->lnodes, node)
    realloc_nodevariables(node, nvdb_new);

  mesh->nvdb = nvdb_new;
}

/**********************************************************************/
/* storage for variables in the nodes */
/**********************************************************************/

/* enable one component of variable i on one node */
void enablevarcomp_innode(tNode *node, int i)
{
  tDat *dat = node->dat;

  /* do nothing if no data is stored on this node */
  if(dat==NULL) return;

  if(i>=dat->nv) errorexiti("var comp %i does not exist", i);
  if(!dat->v[i])
  {
    dat->v[i] = alloc_array(node->n); //, node, NODE);
    dat->nvenabled++;
  }
}

/* disable one component of variable */
void disablevarcomp_innode(tNode *node, int i) 
{
  tDat *dat = node->dat;

  /* do nothing if no data is stored on this node */
  if(dat==NULL) return;

  if(i>=dat->nv) errorexiti("var comp %i does not exist", i);
  if(!dat->v[i])
  {
    free_array(dat->v[i]);
    dat->nvenabled--;
  }
}

/* enable all components of a variable on one node */
void enablevar_innode(tNode *node, int i)
{
  tMesh *mesh = node->pat->mesh;
  int j, n = MeshVarNComponents(mesh, i);
  for(j=0; j<n; j++) enablevarcomp_innode(node, i+j);
}

/* disable all components of a variable on one node */
void disablevar_innode(tNode *node, int i)
{
  tMesh *mesh = node->pat->mesh;
  int j, n = MeshVarNComponents(mesh, i);
  for(j=0; j<n; j++) disablevarcomp_innode(node, i+j);
}



/* enable all components of a variable on one pat */
void enablevar_inpatch(tPat *pat, int i)
{
  //fornodes(pat->ln)
  //  enablevar_innode(node, i);
}

/* disable all components of a variable on one pat */
void disablevar_inpatch(tPat *pat, int i)
{
  //fornodes(pat->ln)
  //  disablevar_innode(node, i);
}


/* enable all components of a variable on one mesh */
void enablevar(tMesh *mesh, int i)
{
  //forallpats(mesh)
  //  enablevar_inpatch(pat, i)
}

/* disable all components of a variable on one mesh */
void disablevar(tMesh *mesh, int i)
{
  //forallpats(mesh)
  //  disablevar_inpatch(pat, i)
}


/* enable variable list */
void enablevarlist_innode(tNode *node, tVarList *vl)
{
  int i;
  if(vl) for(i=0; i<vl->n; i++) enablevarcomp_innode(node, vl->index[i]);
}

/* disable variable list */
void disablevarlist_innode(tNode *node, tVarList *vl)
{
  int i;
  if(vl) for(i=0; i<vl->n; i++) disablevarcomp_innode(node, vl->index[i]);
}

// may need enablevarlist rather for patch or mesh ???




/* Bfaces stuff */
// /* add a bface to a pat, f denotes the pat face (0 to 5),
//    return index in bface list */
// int add_empty_bface(tPat *pat, int f)
// {
//   int fi = pat->nbfaces; /* add bface in this pos. in bface list */
//   void *ret;
// 
//   /* increase size of bface list */
//   ret = realloc( pat->bface, (sizeof( *(pat->bface) ))*(fi+1) );
//   if(ret==NULL)  errorexit("add_bface: not enough memory for pat->bface");
//   pat->bface = ret;
// 
//   /* mem for new bface */
//   ret = calloc( 1, sizeof( *(pat->bface[fi]) ) );
//   if(ret==NULL) errorexit("add_bface: not enough memory for pat->bface[n]");
//   /* add new bface */
//   pat->bface[fi] = ret;
//   pat->nbfaces = fi+1;
// 
//   /* set some bface info */
//   pat->bface[fi]->mesh = pat->mesh;
//   pat->bface[fi]->b    = pat->b;
//   pat->bface[fi]->f    = f;
//   pat->bface[fi]->fi   = fi;
//   pat->bface[fi]->ob   = -1; /* other pat is not known yet */
//   pat->bface[fi]->ofi  = -1; /* fi in other pat is not known yet */
//   pat->bface[fi]->oXi  = -1; /* var indices of other coords not known yet */
//   pat->bface[fi]->oYi  = -1; /* var indices of other coords not known yet */
//   pat->bface[fi]->oZi  = -1; /* var indices of other coords not known yet */
//   return fi;
// }
// 
// /* add a point ijk on face f to a bface with index fi, returns fi */
// /* if called with fi<0, it first calls add_empty_bface and returns the new fi */
// int add_point_to_bface_inpat(tPat *pat, int fi, int ijk, int f)
// {
//   tBface *bface;
//   /* make new bface when needed */
//   if(fi<0 || fi>=pat->nbfaces) fi = add_empty_bface(pat, f);
//   bface = pat->bface[fi];
//   /* make PointList id needed */
//   if(bface->fpts==NULL)
//     bface->fpts = AllocatePointList(pat->mesh);
//   /* add point ijk */
//   AddToPointList(bface->fpts, pat->b, ijk);
//   /* check if we have more than one face on this bface */
//   if(bface->f!=f) bface->f = -1;
//   return fi;
// }
// 
// /* duplicate bface without pointlist fpts */
// tBface *duplicate_bface_without_fpts_for_mesh(tBface *bface0, tMesh *mesh)
// {
//   tBface *bface;
//   void *ret;
//   if(bface0==NULL) return NULL;
// 
//   /* mem for new bface */
//   ret = calloc( 1, sizeof( *(bface) ) );
//   if(ret==NULL)
//     errorexit("duplicate_bface_without_fpts_for_mesh: not enough memory");
//   bface = ret;
// 
//   /* make a shallow copy of the struct */
//   *bface = *bface0;
//   /* now set mesh pointer */
//   bface->mesh = mesh;
//   /* remove pointer to bface0->fpts */
//   bface->fpts = NULL;
//   
//   return bface;
// }
// 
// /* duplicate bface with pointlist fpts */
// tBface *duplicate_bface_for_mesh(tBface *bface0, tMesh *mesh)
// {
//   tBface *bface = duplicate_bface_without_fpts_for_mesh(bface0, mesh);
//   /* copy fpts */
//   bface->fpts = DuplicatePointList_for_mesh(bface0->fpts, mesh);
//   return bface;
// }
// 
// /* free a bface */
// void free_bface(tBface *bface)
// {
//   if(bface!=NULL)
//   {
//     /* free the point lists of the faces */ 
//     FreePointList(bface->fpts);
//     free(bface);
//   }
// }
// 
// /* free all bfaces in pat */
// void free_all_bfaces(tPat *pat)
// {
//   int i;
//   for(i=0; i<pat->nbfaces; i++)  free_bface(pat->bface[i]);
//   free(pat->bface);
//   pat->bface = NULL;
//   pat->nbfaces = 0;
// }
// 
// /* remove a bface with index fi from a pat, return number of bfaces removed */
// int remove_bface(tPat *pat, int fi)
// {
//   void *ret;
//   int nbfaces = pat->nbfaces;
//   int i;
// 
//   /* return 0 if bface does not exist */
//   if(fi<0 || fi>=nbfaces) return 0;
// 
//   /* free the bface */
//   free_bface(pat->bface[fi]);
// 
//   /* shift bfaces behind fi one position to the front */
//   for(i=fi; i<nbfaces-1; i++)
//   {
//     pat->bface[i] = pat->bface[i+1];
//     pat->bface[i]->fi = i;
//   }
// 
//   /* reduce size of bface list */
//   nbfaces--;
//   ret = realloc( pat->bface, (sizeof( *(pat->bface) ))*(nbfaces) );
//   if(ret==NULL && nbfaces!=0)
//     errorexit("remove_bface: not enough memory for pat->bface");
//   pat->bface = ret;
//   pat->nbfaces = nbfaces;
// 
//   return 1;
// }
// 
// /* look for empty bfaces and remove them */
// int remove_bfaces_with_NULL_fpts(tPat *pat)
// {
//   int n, fi;
//   n=0;
//   for(fi=0; fi<pat->nbfaces; fi++)
//     if(pat->bface[fi]->fpts == NULL)
//     {
//       int r = remove_bface(pat, fi); /* this decreases pat->nbfaces */
//       fi=fi-r;  /* go back by one in fi so that for-loop covers fi again */
//       n++;      /* count number of bfaces removed */
//     }
//   return n;
// }


