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
/* node storage */
/**************************************************************************/

/* allocate one node*/
tNode *alloc_node(void)
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
tNode *make_root_node(tPat *pat, int n[3], int datrank)
{
  tNode *node = alloc_node();
  int i, nvdb;

  /* fill in info */
  node->pat = pat;
  /* node->nb is left uninitialized here !!! */
  /* we assume that it has no neighbors */

  /* take bounding boxes from pat, and mark it as touching all patch faces */
  for(i=0; i<6; i++)
  {
    node->bbox[i] = pat->bbox[i];
    node->patface[i] = 1;
  }

  for(i=0; i<3; i++) node->n[i] = n[i];
  node->np = n[0] * n[1] * n[2];
  node->l = 0;
  node->leaf = 1;    /* make this a leaf node */
  nvdb = pat->mesh->nvdb;

  /* get node->D from patch */
  for(i=0; i<3; i++) node->D[i] = node->pat->D[n[i]][i];

  /* see where dat needs to be allocated */
  node->datrank = datrank;
  if(nMPI_rank()==datrank)
    node->dat = alloc_dat(nvdb);

  return node;
}

/* make a child node */
tNode *make_child_node(tNode *parent, int n[3], int ijk)
{
  tNode *node = alloc_node();
  double mid[3];
  int i,j,k, nvdb;
  int ns[] = {2,2,2};

  /* register this child with the parent */
  parent->child[ijk] = node;
  parent->leaf = 0;  /* parent is now no longer a leaf node */

  /* node coords from node index ijk */
  k = kOfInd_n(ijk, ns);
  j = jOfInd_n_k(ijk, ns,k);
  i = iOfInd_n_jk(ijk, ns,j,k);

  /* mid point in parent node */
  for(i=0; i<3; i++)
    mid[i] = 0.5*(parent->bbox[2*i] + parent->bbox[2*i+1]);

  /* set new bounding boxes */
  /* at first take bbox from pat, patface from parent */
  for(i=0; i<6; i++)
  {
    node->bbox[i]    = parent->bbox[i];
    node->patface[i] = parent->patface[i];
  }
  if(i%2) { node->bbox[0] = mid[0]; node->patface[0] = 0; }
  else    { node->bbox[1] = mid[0]; node->patface[1] = 0; }
  if(j%2) { node->bbox[2] = mid[1]; node->patface[2] = 0; }
  else    { node->bbox[3] = mid[1]; node->patface[3] = 0; }
  if(k%2) { node->bbox[4] = mid[2]; node->patface[4] = 0; }
  else    { node->bbox[5] = mid[2]; node->patface[5] = 0; }

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
    node->datrank = parent->datrank;
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
  tNlist *nlist, el;
  tNlist *elem = NULL;
  tNode *node, *onode;
  tNode *narray[8];
  int ijk;

  for(ijk=0; ijk<7; ijk++)
  {
    node = make_child_node(parent, n, ijk);
    elem = addnode_to_nodelist_after(elem, node);
    if(ijk==0) nlist = elem; // save list head
    narray[ijk] = node; /* save nodes also in an array */
  }
  /* fill in neighbor info, as fas as these 8 are concerned */
  connect8_siblings(narray);

  connect8_with_neighbors(narray);

  return nlist;
}

/* remove leaves */
tNode *remove8_leaf_nodes(tNode *leaf0)
{
  tNode *parent = leaf0->parent;
  tNode *node;

  if(!leaf0->leaf) errorexit("argument leaf0 is not a leaf node");

  /* parent is now a leaf node */
  parent->leaf = 1;

  /* update lns lists in patch here ??? NO! make separate func */
  //...

  /* update neighbor info */
  //...

  //for node \in {parent's children}
  //  free_node(node);

  return parent;
}


/* replace current entry in leaf node list with its 8 childern */
void insert8_childnodes_asleaves(tNlist *elem, int n[3])
{
  tNode *parent = elem->node;
  tNlist *children = make8_child_nodes(parent, n);
  replace1_in_nodelist(elem, children);
}

/**************************************************************************/
/* patch and mesh storage */
/**************************************************************************/

/* allocate patch */
tPat *alloc_patch(tMesh *mesh, int p, int nD)
{
  int i;
  tPat *pat;

  pat = calloc(1, sizeof(*pat));
  if(!pat) errorexit("out of memory");

  pat->mesh = mesh;
  pat->p = p;
  pat->nD = nD;

  /* get mem. for diff. matrices */
  pat->D = calloc(nD, sizeof(pat->D[0][0]));
  if(!(pat->D) )
    errorexit("out of memory for diff. matrices");

  /* Bfaces */

  return pat;
}

/* free pat, currently leaves mesh untouched */
void free_patch(tPat *pat)
{
  int i;

  if (!pat) return;

PRF;printf(" isn't working yet!!!\n");
  //for (i = 0; i < pat->mesh->nvariables; i++)
  //  disablevarcomp_inpat(pat, i);
  //free(pat->v);

  //free_all_bfaces(pat);

  //free_nodesinlist(tNlist *elem)

  free(pat);
}



/* allocate mesh */
tMesh *alloc_mesh(int npats)
{
  tMesh *mesh;

  mesh = calloc(1, sizeof(*mesh));
  if(!mesh) errorexit("out of memory for mesh");

  realloc_patlist_in_mesh(mesh, npats);

  return mesh;
}

/* make room for more patches */
void realloc_patlist_in_mesh(tMesh *mesh, int npats)
{
  int opats = mesh->npats;

//printf("npats=%d sss=%d\n", npats, sizeof(mesh->pat[0]));

  /* alloc list of pointers to patches */
  if(npats > opats)
  {
    mesh->pat = realloc(mesh->pat, npats*sizeof(mesh->pat[0]));
    if(!mesh->pat) errorexit("out of memory for mesh->pat");

    /* zero newly allocated part */
    memset(&(mesh->pat[opats]), 0, (npats-opats)*sizeof(mesh->pat[0]));
  }
  if(npats < opats)
  {
    int p;
    for(p=npats; p<opats; p++)
      free_patch(mesh->pat[p]);
    mesh->pat = realloc(mesh->pat, npats*sizeof(mesh->pat[0]));
    if(!mesh->pat) errorexit("cannot shrink mesh->pat");
  }
  mesh->npats = npats;
}

/* free mesh */
void free_mesh(tMesh *mesh)
{
  int i;

  if(!mesh) return;

  /* free patches */
  for(i = 0; i < mesh->npats; i++)
    free_patch(mesh->pat[i]);

  /* free vdb and pdb in mesh */
  free(mesh->vdb);
  free(mesh->pdb);

  /* node list in mesh */
  free_nodelist(mesh->lns);

  free(mesh);
}

/**********************************************************************/
/* storage for lists of nodes */
/**********************************************************************/
/* allocate an empty node list */
tNlist *alloc_nodelist(tNode *node)
{
  tNlist *nlist;
  nlist = calloc(1, sizeof(*nlist));
  if(!nlist) errorexit("out of memory for nlist");
  nlist->node = node;
  return nlist;
}

/* add one node to nodelist after elem, and return new nodelist element
   that now contains the node */
tNlist *addnode_to_nodelist_after(tNlist *elem, tNode *node)
{
  tNlist *after = alloc_nodelist(node);
  return insertnodelist_into_nodelist_after(elem, after);
}
/* add one node to nodelist before elem, and return new nodelist element
   that now contains the node */
tNlist *addnode_to_nodelist_before(tNlist *elem, tNode *node)
{
  tNlist *before = alloc_nodelist(node);
  return insertnodelist_into_nodelist_before(elem, before);
}

/* insert nodelist "list" into another nodelist after elem,
   and return the end of "list" */
tNlist *insertnodelist_into_nodelist_after(tNlist *elem, tNlist *list)
{
  tNlist *elem2;
  tNlist *lend;
  tNlist *lbeg;

  /* find end and beginning of tNlist *list */
  for(lend=list; lend->next; lend=lend->next) ;
  for(lbeg=list; lbeg->prev; lbeg=lbeg->prev) ;

  if(!elem) return lend;

  elem2 = elem->next;
  lend->next = elem2;
  lbeg->prev = elem;
  elem->next = lbeg;
  if(elem2) elem2->prev = lend;
  return lend;
}

/* insert nodelist "list" into another nodelist before elem,
   and return the first of "list" */
tNlist *insertnodelist_into_nodelist_before(tNlist *elem, tNlist *list)
{
  tNlist *elem2;
  tNlist *lend;
  tNlist *lbeg;

  /* find end and beginning of tNlist *list */
  for(lend=list; lend->next; lend=lend->next) ;
  for(lbeg=list; lbeg->prev; lbeg=lbeg->prev) ;

  if(!elem) return lbeg;

  elem2 = elem->prev;
  lend->next = elem;
  lbeg->prev = elem2;
  elem->prev = lend;
  if(elem2) elem2->next = lbeg;
  return lbeg;
}


tNlist *replace1_in_nodelist(tNlist *elem, tNlist *list)
{
  tNlist *left;
  tNlist *right;
  tNlist *lend;
  tNlist *lbeg;

  /* find end and beginning of tNlist *list */
  for(lend=list; lend->next; lend=lend->next) ;
  for(lbeg=list; lbeg->prev; lbeg=lbeg->prev) ;

  if(!elem) return lend;

  left = elem->prev;
  right= elem->next;

  lend->next = right;
  lbeg->prev = left;
  if(right) right->prev = lend;
  if(left)  left->next = lbeg;
  return lbeg;
}

/* remove 1 element from nodelist, and return element after elem */
tNlist *remove1_in_nodelist(tNlist *elem)
{
  tNlist *left;
  tNlist *right;
  if(!elem) return 0;

  left = elem->prev;
  right= elem->next;
  if(right) right->prev = left;
  if(left)  left->next = right;
  free(elem);
  return right;
}

/* remove all from nodelist and free it */
void free_nodelist(tNlist *elem)
{
  tNlist *tmp;

  if(!elem) return;

  /* remove all after elem */
  while(tmp=elem->next)
    remove1_in_nodelist(tmp);

  /* remove all before elem */
  while(tmp=elem->prev)
    remove1_in_nodelist(tmp);

  /* remove elem */
  remove1_in_nodelist(elem);
}

void free_nodesinlist(tNlist *elem)
{
  tNlist *tmp;

  /* free nodes in elem and all after in */
  for(tmp=elem; tmp; tmp=tmp->next)
    free_node(tmp->node);

  /* free nodes in all before elem */
  for(tmp=elem->prev; tmp; tmp=tmp->prev)
    free_node(tmp->node);
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
  tNlist *elem;
  tNode *node;
  if(0) printf("realloc_meshvariables from %d to %d\n",
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
  forlnodes(mesh, node) {
    realloc_nodevariables(node, nvdb_new);
  } endforlnodes;
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

  forlnodes(pat, node) {
    enablevarcomp_innode(node, i);
  } endforlnodes;
}

/* disable one component of a variable on one pat */
void disablevarcomp_inpatch(tPat *pat, int i)
{
  tNode *node;

  forlnodes(pat, node) {
    disablevarcomp_innode(node, i);
  } endforlnodes;
}

/* enable all components of a variable on one pat */
void enablevar_inpatch(tPat *pat, int i)
{
  tNode *node;

  forlnodes(pat, node) {
    enablevar_innode(node, i);
  } endforlnodes;
}

/* disable all components of a variable on one pat */
void disablevar_inpatch(tPat *pat, int i)
{
  tNode *node;

  forlnodes(pat, node) {
    disablevar_innode(node, i);
  } endforlnodes;
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


