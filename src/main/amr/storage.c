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

/* replace pointer to data with something else */
void point_array_a_to_data(tArray *array, void *data)
{
  free(array->a);
  array->a = data;
}

/* free an array */
void free_array(tArray *array)
{
  if(!array) return;
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

/* free one node only, leaves children hanging */
void free_this_node_only(tNode *node)
{
  tNode *parent = node->parent;
  int ijk, face;

  if(!node) return;

  /* free surface neigbhor list */
  for(face=0; face<6; face++) free(node->fnb[face]);

  /* remove parent's pointer to it */
  ijk = node->ijk;
  if(parent) parent->child[ijk] = NULL;

  /* should we also remove pointer of neighbors to it??? */

  /* free variable data */
  free_dat(node->dat);
  free(node);
}

/* free one node and all its children */
void free_node(tNode *node)
{
  tNode *chld;
  int ijk;

  //PRF;printf(": nid=%d l=%d node=%p\n", get_node_nid(node), node->l, node);
  //printnode(node);
  if(!node) return;

  for(ijk=0; ijk<8; ijk++)
  {
    chld = node->child[ijk];
    if(chld) free_node(chld);
  }
  //PRF;printf(": nid=%d l=%d node=%p\n", get_node_nid(node), node->l, node);
  free_this_node_only(node);
}

/* set node pointers for node->Dt ... . Point them to arrays from patch */
void point_nodearrays_to_patarrays(tPat *pat, tNode *node)
{
  int *n = node->n;
  int dir;
  /* get node->Dt ... from patch */
  for(dir=0; dir<3; dir++)
  {
    node->Dt[dir] = node->pat->Dt[n[dir]][dir];
    node->At[dir] = node->pat->At[n[dir]][dir];
    node->St[dir] = node->pat->St[n[dir]][dir];
    node->Xb[dir] = node->pat->Xb[n[dir]][dir];
    node->Winteg[dir] = node->pat->Winteg[n[dir]][dir];
  }
}


/* make root node */
tNode *make_root_node(tPat *pat, int n[3], int datrank)
{
  tNode *node = alloc_node();
  int i;

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

  /* get node->Dt ... from patch */
  point_nodearrays_to_patarrays(pat, node);

  /* see where dat needs to be allocated */
  node->datrank = datrank;
  if(nMPI_rank()==datrank)
    node->dat = alloc_dat(node);

  /* initialize surface neigbhor list in dat */
  update_node_and_neighbors_fnb(node);

  return node;
}

/* make a child node */
tNode *make_child_node(tNode *parent, int n[3], int ijk)
{
  tNode *node = alloc_node();
  double mid[3];
  int i,j,k, d, nvdb;
  int ns[] = {2,2,2};

  /* register this child with the parent */
  parent->child[ijk] = node;
  parent->leaf = 0;  /* parent is now no longer a leaf node */

  /* node coords from node index ijk */
  k = kOfInd_n(ijk, ns);
  j = jOfInd_n_k(ijk, ns,k);
  i = iOfInd_n_jk(ijk, ns,j,k);

  /* mid point in parent node */
  for(d=0; d<3; d++)
    mid[d] = 0.5*(parent->bbox[2*d] + parent->bbox[2*d+1]);

  /* set new bounding boxes */
  /* at first take bbox from pat, patface from parent */
  for(d=0; d<6; d++)
  {
    node->bbox[d]    = parent->bbox[d];
    node->patface[d] = parent->patface[d];
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

  for(d=0; d<3; d++) node->n[d] = n[d];
  node->np = n[0] * n[1] * n[2];

  node->l = parent->l + 1;
  node->leaf = 1;    /* make this a leaf node */
  node->ijk = ijk;
  nvdb = node->pat->mesh->nvdb;

  /* get node->Dt ... from patch */
  point_nodearrays_to_patarrays(node->pat, node);

  /* default is same proc as parent */
  node->datrank = parent->datrank;

  /* if parent has dat the child will have it too */
  if(parent->dat)
  {
    node->dat = alloc_dat(node);
    /* enable same vars in this dat as in parent->dat */
    for(d=0; d<nvdb; d++)
      if(parent->dat->v[d])  enablevarcomp_innode(node, d);

    /* fill node->dat with interpolation data from parent */
    // still TODO
  }

  return node;
}

/* make 8 childern and return them in a short list */
tNlist *make8_child_nodes(tNode *parent, int n[3])
{
  tNlist *nlist;
  tNlist *elem = NULL;
  tNode *node;
  tNode *narray[8];
  int ijk;

  for(ijk=0; ijk<8; ijk++)
  {
    node = make_child_node(parent, n, ijk);
    elem = addnode_to_nodelist_after(elem, node);
    if(ijk==0) nlist = elem; // save list head
    narray[ijk] = node; /* save nodes also in an array */
  }
  /* fill in neighbor info, as far as these 8 are concerned */
  connect8_with_neighbors(narray, 1);

  /* update fnb of all in narray and their neighbors */
  for(ijk=0; ijk<8; ijk++)
    update_node_and_neighbors_fnb(narray[ijk]);

  /* free all data on parent */
  free_dat(parent->dat);
  parent->dat = NULL;

  return nlist;
}

/* remove children */
tNode *destroy_children(tNode *parent)
{
  tNode *narray[8];
  tNode *child0 = parent->child[0];
  int ijk;

  for(ijk=0; ijk<8; ijk++)
    narray[ijk] = parent->child[ijk]; /* save children in an array */

  /* set parents datrank to the same as child0 */
  parent->datrank = child0->datrank;
  if(child0->dat)
  {
    int nvdb = parent->pat->mesh->nvdb;
    int d;

    if(!parent->dat) parent->dat = alloc_dat(parent);
    /* enable same vars in this dat as in parent->dat */
    for(d=0; d<nvdb; d++)
      if(child0->dat->v[d])  enablevarcomp_innode(parent, d);

    /* fill node->dat with interpolation data from parent */
    // still TODO
  }

  /* update neighbor info */
  /* set neighbor info to NULL, as far as these 8 are concerned */
  connect8_with_neighbors(narray, 0);

  /* free child nodes */
  for(ijk=0; ijk<8; ijk++)
  {
    if(narray[ijk]->child[0])
      errorexit("cannot destroy child that itself has child[0]");
    free_node(narray[ijk]);
    parent->child[ijk] = NULL;
  }

  /* update fnb on parent and its neighbors */
  update_node_and_neighbors_fnb(parent);

  /* parent is now a leaf node */
  parent->leaf = 1;

  return parent;
}

/**************************************************************************/
/* patch and mesh storage */
/**************************************************************************/

/* allocate patch */
tPat *alloc_patch(tMesh *mesh, int p, int nmax)
{
  int n[3];
  int i, d;
  tPat *pat;

  pat = calloc(1, sizeof(*pat));
  if(!pat) errorexit("out of memory");

  pat->mesh = mesh;
  pat->p = p;
  pat->nmax = nmax;

  /* get mem. for diff. matrices */
  pat->Dt = calloc(nmax+1, sizeof(pat->Dt[0]));
  if(!(pat->Dt) )
    errorexit("out of memory for diff. matrices");
  pat->At = calloc(nmax+1, sizeof(pat->At[0]));
  if(!(pat->At) )
    errorexit("out of memory for ana. matrices");
  pat->St = calloc(nmax+1, sizeof(pat->St[0]));
  if(!(pat->St) )
    errorexit("out of memory for syn. matrices");
  pat->Xb = calloc(nmax+1, sizeof(pat->Xb[0]));
  if(!(pat->Xb) )
    errorexit("out of memory for points");
  pat->Winteg = calloc(nmax+1, sizeof(pat->Winteg[0]));
  if(!(pat->Winteg) )
    errorexit("out of memory for integr. weights");
  for(d=1; d<=nmax; d++)
  {
    n[0] = n[1] = d;
    n[2] = 1;
    for(i=0; i<3; i++)
    {
      pat->Dt[d][i] = alloc_array(n);
      pat->At[d][i] = alloc_array(n);
      pat->St[d][i] = alloc_array(n);
    }
    n[0] = d;
    n[1] = n[2] = 1;
    for(i=0; i<3; i++)
    {
      pat->Xb[d][i] = alloc_array(n);
      pat->Winteg[d][i] = alloc_array(n);
    }
  }

  /* Bfaces */

  return pat;
}

/* free pat, currently leaves mesh untouched */
void free_patch(tPat *pat)
{
  tMesh *mesh = pat->mesh;
  tNlist *elem;
  int d, i;

  if (!pat) return;

  PRFs(":\n");

  /* free diff matrices and such */
  for(d=1; d<=pat->nmax; d++)
    for(i=0; i<3; i++)
    {
      free_array(pat->Dt[d][i]);
      free_array(pat->At[d][i]);
      free_array(pat->St[d][i]);
      free_array(pat->Xb[d][i]);
      free_array(pat->Winteg[d][i]);
    }
  free(pat->Dt);
  free(pat->At);
  free(pat->St);
  free(pat->Xb);
  free(pat->Winteg);

  //free_all_bfaces(pat);
  PRF;printf(": implement free_all_bfaces!!!\n");

  /* remove all nodes from this patch from mesh->lns */
  for(elem=mesh->lns; elem; )
  {
    if(elem->node->pat == pat) elem = remove1_in_nodelist(elem, 1);
    else                       elem = elem->next;
  }

  mesh->lns = first_nodelist(elem);
  update_mesh_myln_node_nid(mesh);

  /* free root node and all its children ... */
  free_node(pat->rnode);

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

  PRFs(":\n");

  /* free patches */
  for(i = 0; i < mesh->npats; i++)
    free_patch(mesh->pat[i]);

  /* now patch list array */
  free(mesh->pat);

  /* node list in mesh */
  free_nodelist(mesh->lns);
  mesh->lns = NULL;
  free(mesh->myln);
  mesh->myln = NULL;
  mesh->nmyln = 0;

  /* free vdb and pdb in mesh */
  free_mesh_vdb_contents(mesh);
  free(mesh->vdb);
  free_mesh_pdb_contents(mesh);
  free(mesh->pdb);

  /* free skeleton in mesh */
  remove_all_MeshFuns(mesh);

  free(mesh);
}

/**********************************************************************/
/* storage for lists of nodes */
/**********************************************************************/
/* allocate a node list with one node */
/* NOTE: we can also add to a nodelist that is NULL, so alloc_nodelist
   is not always needed */
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

/* make a copy */
tNlist *copy_of_nodelist(tNlist *elem)
{
  tNlist *dest = NULL;
  tNlist *src  = first_nodelist(elem);
  tNlist *el;

//printnodelist(src);
  for(el=src; el; el=el->next)
    dest = addnode_to_nodelist_after(dest, el->node);
//printnodelist(dest);

  return dest;
}

/* count num. of elem. in list */
int count_elements_nodelist(tNlist *list)
{
  tNlist *beg = first_nodelist(list);
  tNlist *el;
  int count=0;

  /* count elem. in list */
  for(el=beg; el; el=el->next) count++;

  return count;
}

/* replace element elem in nodelist by node */
tNlist *replacenode_in_nodelist(tNlist *elem, tNode *node)
{
  tNlist *repl = alloc_nodelist(node);
  return replace1_in_nodelist(elem, repl);
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

/* replace 1 element in a nodelist by a list and then free the element */
tNlist *replace1_in_nodelist(tNlist *elem, tNlist *list)
{
  tNlist *left;
  tNlist *right;
  tNlist *lend;
  tNlist *lbeg;

  /* find end and beginning of tNlist *list */
  for(lend=list; lend->next; lend=lend->next) ;
  for(lbeg=list; lbeg->prev; lbeg=lbeg->prev) ;

  if(!elem) return lbeg;

  left = elem->prev;
  right= elem->next;

  lend->next = right;
  lbeg->prev = left;
  if(right) right->prev = lend;
  if(left)  left->next = lbeg;

  free(elem);
  return lbeg;
}

/* replace 1 element in a nodelist by a list, then free the element,
   return the very first element of the new list */
tNlist *first_replace1_in_nodelist(tNlist *elem, tNlist *list)
{
  tNlist *newlist = replace1_in_nodelist(elem, list);
  return first_nodelist(newlist);
}

/* remove 1 element from nodelist, and
   return element after elem if return_next=1 */
tNlist *remove1_in_nodelist(tNlist *elem, int return_next)
{
  tNlist *left;
  tNlist *right;
  if(!elem) return 0;

  left = elem->prev;
  right= elem->next;
  if(right) right->prev = left;
  if(left)  left->next = right;

  free(elem);
  if(return_next) return right;
  else		  return left;
}

/* return 1st element in a nodelist */
tNlist *first_nodelist(tNlist *list)
{
  tNlist *lbeg;

  if(!list) return NULL;

  /* find beginning of tNlist *list */
  for(lbeg=list; lbeg->prev; lbeg=lbeg->prev) ;
  return lbeg;
}
/* return last element in a nodelist */
tNlist *last_nodelist(tNlist *list)
{
  tNlist *lend;

  if(!list) return NULL;

  /* find end of tNlist *list */
  for(lend=list; lend->next; lend=lend->next) ;
  return lend;
}


/* remove all from nodelist and free it */
void free_nodelist(tNlist *elem)
{
  tNlist *tmp;

  if(!elem) return;

  /* remove all after elem */
  for(tmp=elem->next; tmp; )
    tmp = remove1_in_nodelist(tmp, 1);

  /* remove all before elem */
  for(tmp=elem->prev; tmp; )
    tmp = remove1_in_nodelist(tmp, 0);

  /* remove elem */
  remove1_in_nodelist(elem, 0);
}

/* free all nodes in a list */
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
/* functions to update the nodelist and node array in mesh */
/**********************************************************************/
/* update array of leaf nodes on this proc */
int update_mesh_myln_node_nid(tMesh *mesh)
{
  tNlist *elem;
  int allocd = mesh->nmyln;
  int ainc = 256;
  int nid = 0;
  int nmyln = 0;

  /* go over leaves if  mesh->lns is not NULL */
  if(mesh->lns) fornodelist(mesh->lns, elem)
  {
    tNode *node = elem->node;
    tNode *parent = node->parent;

    if(node->dat)
    {
      if(nmyln >= allocd)
      {
        mesh->myln = realloc(mesh->myln, sizeof(mesh->myln[0])*(allocd+ainc));
        allocd += ainc;
      }
      mesh->myln[nmyln++] = elem;
    }
    /* set nid and invalidate parent's nid */
    node->nid = nid++;
    if(parent) parent->nid = -nid;
  }
  else /* mesh->lns is NULL, so free myln array */
  {
    free(mesh->myln);
    mesh->myln = NULL;
  }
  mesh->nmyln = nmyln;
  return nid;
}

/* return nid or -1 */
int get_node_nid(tNode *node)
{
  return node ? node->nid : -1;
}

/* append a node list to mesh->lns and also update mesh->myln */
tNlist *append_nodelist_to_mesh_lns_myln(tMesh *mesh, tNlist *list)
{
  tNlist *lnl;
  if(mesh->lns)
  {
    lnl = last_nodelist(mesh->lns); /* last elem. in mesh->lns */
    lnl = insertnodelist_into_nodelist_after(lnl, list);
  }
  else
    mesh->lns = first_nodelist(list);

  update_mesh_myln_node_nid(mesh);
  return lnl;
}

/* replace elem in mesh->lns by nlist, return first of nlist */
tNlist *replace1_in_mesh_lns_myln(tNlist *elem, tNlist *nlist)
{
  tNlist *nlist_beg;
  tMesh *mesh = NULL;

  if(elem) mesh = elem->node->pat->mesh;
  else     errorexit("elem is NULL!!!");

  nlist_beg = replace1_in_nodelist(elem, nlist);
  mesh->lns = first_nodelist(nlist_beg);
  update_mesh_myln_node_nid(mesh);
  return nlist_beg;
}

/* replace current entry in leaf node list with its 8 new childern,
   return element with 0th child */
tNlist *make8children_in_mesh_lns_myln(tNlist *elem, int n[3])
{
  tNode *parent = elem->node;
  tNlist *children = make8_child_nodes(parent, n);
  return replace1_in_mesh_lns_myln(elem, children);
}

/* replace siblings at element sib of mesh->lns by parent,
   node with sibling 0 is returned so we can destroy it later */
tNode *remove8siblings_in_mesh_lns_myln(tNlist *sib)
{
  tNode *parent, *node0;
  tNlist *elem, *elem0;
  tMesh *mesh = NULL;
  int ijk;

  if(sib==NULL) errorexit("sib is NULL!!!");
  mesh = sib->node->pat->mesh;
  parent = sib->node->parent;
  if(parent==NULL) errorexit("parent is NULL!!!");

  /* find sibling 0 */
  elem0=sib;
  for(ijk=sib->node->ijk; ijk>0; ijk--)
    elem0=elem0->prev;
  node0 = elem0->node;
  if(node0->parent != parent || node0->ijk != 0)
    errorexit("elem0 has wrong parent!");

  /* set elem to sibling 1 and remove the 7 after sibling 0 */
  elem=elem0->next;
  for(ijk=1; ijk<8; ijk++)
  {
    if(elem->node->parent != parent) errorexit("elem has wrong parent!");
    elem = remove1_in_nodelist(elem, 1);
  }

  /* replace sibling 0 by parent in list */
  elem0 = replacenode_in_nodelist(elem0, parent);

  /* reset mesh lists */
  mesh->lns = first_nodelist(elem0);
  update_mesh_myln_node_nid(mesh);
  return parent;
}

/* replace siblings at element sib of mesh->lns by their parent, and then
   destroy the 8 siblings */
void destroy8siblings_in_mesh_lns_myln(tNlist *sib)
{
  tNode *parent = remove8siblings_in_mesh_lns_myln(sib);
  destroy_children(parent);
}

/**********************************************************************/
/* allocate and fill surfaces for vars that need it */
/**********************************************************************/
/* empty surface that we need to fill in */
tSurface *alloc_empty_surface(int nnb)
{
  tSurface *s = calloc(1, sizeof(*s));
  s->nbsurf = calloc(nnb, sizeof(s->nbsurf[0]));
  s->recv_req = calloc(nnb, sizeof(s->recv_req[0]));
  s->send_req = calloc(nnb, sizeof(s->send_req[0]));
  return s;
}

/* free all we need to in a surface */
void free_surface(tSurface *s)
{
  tNode *node;
  tDat *dat;
  int f, i;

  if(!s) return;
  dat = s->dat;
  f = s->face;
  node = dat->node;

  /* free content of lists */
  if(s->allocd_mysurf) free_array(s->mysurf);
  for(i=0; i<node->nfnb[f]; i++)
    if(!node->fnb[f][i]->dat) free_array(s->nbsurf[i]);

  /* free lists */
  free(s->nbsurf);
  free(s->recv_req);
  free(s->send_req);
}


/* initialize a surface for var vi at face with nnb neighbors */
tSurface *init_surface(tNode *node, int vi, int face)
{
  int dir = face/2;
  int zones;
  tDat *dat;
  int i;
  tSurface *s;
  int n[3];
  int alloc_mysurf;

  /* do nothing if no data on this node */
  if(!node->dat) return NULL;
  dat = node->dat;

  /* do nothing if var is not enabled */
  if(!dat->v[vi]) return NULL;

  /* do nothing if ghost zone width is 0 for this var */
  zones = MeshVarSurfacezones(node->pat->mesh, vi);
  if(zones==0) return NULL;

  /* prep. */
  s = alloc_empty_surface(node->nfnb[face]);
  s->dat = dat;
  s->face = face;
  s->vi = vi;

  /* set n */
  for(i=0; i<3; i++) n[i] = node->n[i];
  alloc_mysurf = 1;
  if(n[dir] == zones) alloc_mysurf = 0;
  else                n[dir] = zones;

  /* allocate my surface array */
  if(alloc_mysurf) s->mysurf = alloc_array(n);
  else             s->mysurf = dat->v[vi];
  s->allocd_mysurf = alloc_mysurf;

  return s;
}

void init_all_surfaces(tNode *node)
{
  int face, ni, j;

  for(face=0; face<6; face++)
  {

  /*
    for(j=0; j<6; j++)
    {
      dat->s[j] = calloc(dat->nv, sizeof(tSurface *));
      if(!dat->s[j]) errorexit("out of memory for dat->s[j]");
    }
  */

  //s = alloc_empty_surface(nnb);
  }
}


/**********************************************************************/
/* storage for dat lists in the nodes */
/**********************************************************************/
/* allocate room for nv variables that can be enabled or disabled */
tDat *alloc_dat(tNode *node)
{
  int nv = node->pat->mesh->nvdb;
  tDat *dat;
  int j;

  dat = calloc(1, sizeof(tDat));
  if(!dat) errorexit("out of memory for dat");

  dat->node = node;
  dat->nv = nv;
  if(nv==0) return dat;

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
  int face, i,j;

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

  if(nv_new<dat->nv) errorexit("implement var removal");

  if(nv_new)
  {
    dat->v = realloc(dat->v, nv_new*sizeof(tArray *));
    if(!dat->v) errorexit("out of memory for dat->v");
  }
  else
  {
    free(dat->v);
    return;
  }

  for(j=0; j<6; j++)
  {
    dat->g[j] = realloc(dat->g, nv_new*sizeof(tArray *));
    if(!dat->g[j]) errorexit("out of memory for dat->g");
  }

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
  int lni;

  if(0) printf("realloc_meshvariables from %d to %d\n",
               mesh->nvdb, nvdb_new);

  /* realloc list on mesh struct */
  if(nvdb_new)
    mesh->vdb = realloc(mesh->vdb, sizeof(tVar)*(nvdb_new));
  else
    free(mesh->vdb);

  /* set newly added stuff to 0 */
  if(nvdb_new > nvdb_old)
  {
    tVar *newv0 = &(mesh->vdb[nvdb_old]);
    memset(newv0, 0, sizeof(tVar)*(nvdb_new-nvdb_old));
  }
  mesh->nvdb = nvdb_new;

  /* now make sure dat in nodes is also reallocated */
  formylnodes(mesh, lni)
  {
    tNode *node = mesh->myln[lni]->node;
    realloc_nodevariables(node, nvdb_new);
  }
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
  tMesh *mesh = pat->mesh;
  tNode *node;

  forlnodes(mesh, node)
  {
    if(node->pat == pat) enablevarcomp_innode(node, i);
  } endforlnodes;
}

/* disable one component of a variable on one pat */
void disablevarcomp_inpatch(tPat *pat, int i)
{
  tMesh *mesh = pat->mesh;
  tNode *node;

  forlnodes(mesh, node)
  {
    if(node->pat == pat) disablevarcomp_innode(node, i);
  } endforlnodes;
}

/* enable all components of a variable on one pat */
void enablevar_inpatch(tPat *pat, int i)
{
  tMesh *mesh = pat->mesh;
  tNode *node;

  forlnodes(mesh, node)
  {
    if(node->pat == pat) enablevar_innode(node, i);
  } endforlnodes;
}

/* disable all components of a variable on one pat */
void disablevar_inpatch(tPat *pat, int i)
{
  tMesh *mesh = pat->mesh;
  tNode *node;

  forlnodes(mesh, node)
  {
    if(node->pat == pat) disablevar_innode(node, i);
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


