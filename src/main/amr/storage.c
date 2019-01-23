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

/* allocate one node*/
tNode *alloc_node()
{
  tNode *node = calloc(1, sizeof(*node));
  if(!node) errorexit("out of memory");

  return node;
}

/* free one node */
void free_node(tNode *node) 
{
  if(!node) return;
  free_dat(node->dat);
  free(node);
}

/* make root node */
Node *make_root_node(tPat *pat, int n[3], int withdat)
{
  tNode *node = alloc_node();
  int i, nvdb;

  /* fill in info */
  node->pat = pat;
  /* node->nb is left uninitialized here !!! */
  /* we assume that it has no neighbors */

  /* take bounding boxes from pat */
  for(i=0; i<6; i++) node->bbox[i] = pat->bbox[i];

  for(i=0; i<3; i++) node->n[i] = n[i];
  node->np = n[0] * n[1] * n[2];
  node->l = 0;
  node->leaf = 1;    /* make this a leaf node */
  nvdb = pat->mesh->nvdb;

  /* get node->D from patch */
  for(i=0; i<3; i++) node->D[i] = node->pat->D[n[i]][i];

  /* if parent has dat the child will have it too */
  if(withdat)
    node->dat = alloc_dat(nvdb);

  return node;
}

/* make a child node */
Node *make_child_node(tNode *parent, int n[3], int ijk)
{
  tNode *node = alloc_node();
  double mid[3];
  int i,j,k, nvdb;
  int ns[] = {2,2,2};

  /* register this child with the parent */
  parent->child[ijk] = node;
  parent->leaf = 0;  /* parent is now no longer a leaf node */

  /* node coords from node index ijk */
  k = kOfInd_n(ijk, n);
  j = jOfInd_n_k(ijk, n,k);
  i = iOfInd_n_jk(ijk,n,j,k);

  /* mid point in parent node */
  for(i=0; i<3; i++)
    mid[i] = 0.5*(parent->bbox[2*i] + parent->bbox[2*i+1]);

  /* set new bounding boxes */
  for(i=0; i<6; i++) node->bbox[i] = parent->bbox[i];

  if(i%2) node->bbox[0] = mid[0];
  else    node->bbox[1] = mid[0];

  if(j%2) node->bbox[2] = mid[1];
  else    node->bbox[3] = mid[1];

  if(k%2) node->bbox[4] = mid[2];
  else    node->bbox[5] = mid[2];

  /* fill in info */
  node->pat    = parent->pat;
  node->parent = parent;
  /* node->nb is left uninitialized here !!! */

  for(i=0; i<3; i++) node->n[i] = n[i];
  node->np = n[0] * n[1] * n[2];

  node->l = parent->l + 1;
  node->leaf = 1;    /* make this a leaf node */
  node->ijk = ijk;
  nvdb = node->pat->mesh->nvdb;

  /* get node->D from patch */
  for(i=0; i<3; i++) node->D[i] = node->pat->D[n[i]][i];

  /* if parent has dat the child will have it too */
  if(parent->dat)
  {
    node->dat = alloc_dat(nvdb);
    /* enable same vars in this dat as in parent->dat */
    for(i=0; i<nvdb; i++)
      if(parent->dat->v[i])  enablevarcomp_innode(node, i);

    /* fill node->dat with interpolation data from parent */
    // still TODO
  }

  return node;
}

/* make 8 childern and return them in a short list */
tNlist *make8_child_nodes(tNode *parent, int n[3])
{
  tNlist *nlist = alloc_nodelist();
  tNode *node;
  int ijk;

  for(ijk=0; ijk<7; ijk++)
  {
    node = make_child_node(parent, n, ijk);
    addto_nodelist(nlist, node);
  }
  /* fill in neighbor info, as fas as these 8 are concerned */
  // still TODO
  return nlist;
}

/* replace current entry in leaf node list with its 8 childern */
void insert8_childnodes_asleaves(tNlist *lnodes, int n[3])
{
  tNode *parent = lnodes->node;
  tNlist *children = make8_child_nodes(parent, n);
  replace1_in_nodelist(lnodes, children);
}


/* allocate patch */
tPat *alloc_patch(tMesh *mesh, int p, int nD) 
{
  tPat *pat;

  pat = calloc(1, sizeof(*pat));
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



/* allocate mesh */
tMesh *alloc_mesh(int npatches)
{
  tMesh *mesh;

  mesh = calloc(1, sizeof(*mesh));
  if(!mesh) errorexit("out of memory for mesh");

  /* alloc list of pointers to patches */
  mesh->pat = calloc(npatches, sizeof(mesh->pat[0]));
  if(mesh->pat) errorexit("out of memory for mesh->pat");

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

/**********************************************************************/
/* storage for lists of nodes */
/**********************************************************************/
tNlist *alloc_nodelist(void)
{
  tNlist *nlist;
  nlist = calloc(1, sizeof(*nlist));
}

void addto_nodelist(tNlist *nlist, tNode *node)
{
  if(nlist->node)
  {
    tNlist *tmp = alloc_nodelist();
    tmp->next = nlist->next;
    tmp->prev = nlist;
    nlist->next = tmp;
  }
  else
    nlist->node = node;
}

void remove1_in_nodelist(tNlist *nlist)
{
  tNode *left = nlist->prev;
  tNode *right= nlist->next;
  if(right) right->prev = left;
  if(left)  left->next = right;
  free(nlist);
}

free_nodelist(tNlist *nlist)
{
  tNode *start=nlist->right;
  tNode *tmp;

  for(tmp=start; tmp->next; tmp=tmp->next)
    remove1_in_nodelist(tmp);

}

/**********************************************************************/
/* storage for dat lists in the nodes */
/**********************************************************************/
/* allocate room for nv variables that can be enabled or disabled */
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
  int nvdb_old = mesh->nvdb;
  tNode *node;
  if(1) printf("realloc_meshvariables from %d to %d\n", 
                 mesh->nvdb, nvdb_new);

  /* realloc list on mesh struct */
  mesh->vdb = realloc(mesh->vdb, sizeof(tVar)*(nvdb_new));

  /* set newly added stuff to 0 */
  if(nvdb_new > mesh->nvdb)
  {
    tVar *newv0 = &(mesh->vdb[nvdb_old]);
    memset(newv0, 0, sizeof(tVar)*(nvdb_new-nvdb_old));
  }
  mesh->nvdb = nvdb_new;

  /* now make sure dat in nodes is also reallocated */
  fornodelist(mesh->lnodes, node)
    realloc_nodevariables(node, nvdb_new);
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

/* enable onr component of a variable on one pat */
void enablevarcomp_inpatch(tPat *pat, int i)
{
  tNode *node;

  fornodelist(pat->lnodes, node)
    enablevarcomp_innode(node, i);
}

/* disable one component of a variable on one pat */
void disablevarcomp_inpatch(tPat *pat, int i)
{
  tNode *node;

  fornodelist(pat->lnodes, node)
    disablevarcomp_innode(node, i);
}

/* enable all components of a variable on one pat */
void enablevar_inpatch(tPat *pat, int i)
{
  tNode *node;

  fornodelist(pat->lnodes, node)
    enablevar_innode(node, i);
}

/* disable all components of a variable on one pat */
void disablevar_inpatch(tPat *pat, int i)
{
  tNode *node;

  fornodelist(pat->lnodes, node)
    disablevar_innode(node, i);
}

/* enable all components of a variable on one mesh */
void enablevar(tMesh *mesh, int i)
{
  int pi;

  forpatches(mesh, pi)
    enablevar_inpatch(mesh->pat[pi], i);
}

/* disable all components of a variable on one mesh */
void disablevar(tMesh *mesh, int i)
{
  int pi;

  forpatches(mesh, pi)
    disablevar_inpatch(mesh->pat[pi], i);
}


/* enable variable list in a node */
void enablevarlist_innode(tNode *node, tVarList *vl)
{
  int i;
  if(vl) for(i=0; i<vl->n; i++) enablevarcomp_innode(node, vl->index[i]);
}

/* disable variable list in a node */
void disablevarlist_innode(tNode *node, tVarList *vl)
{
  int i;
  if(vl) for(i=0; i<vl->n; i++) disablevarcomp_innode(node, vl->index[i]);
}

/* enable variable list */
void enablevarlist(tVarList *vl)
{
  tMesh *mesh = vl->mesh;
  int i, pi;
  if(vl)
  {
    for(i=0; i<vl->n; i++)
      forpatches(mesh, pi)
        enablevarcomp_inpatch(mesh->pat[pi], vl->index[i]);
  }
}

/* disable variable list */
void disablevarlist(tVarList *vl)
{
  tMesh *mesh = vl->mesh;
  int i, pi;
  if(vl)
  {
    for(i=0; i<vl->n; i++)
      forpatches(mesh, pi)
        disablevarcomp_inpatch(mesh->pat[pi], vl->index[i]);
  }
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


