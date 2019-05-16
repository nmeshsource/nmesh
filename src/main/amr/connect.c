/* connect.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"




/* if connect=1:
    we set neighbor connection info of 8 siblings with themselves and with
    nodes from different parents. We assume that all priorly created nodes
    already have complete neighbor info.
   else:
    we set the pointers to neighbors to NULL */
void connect8_with_neighbors(tNode *narray[8], int connect)
{
  int ind;
  int ns[] = {2,2,2};
  tNode *parent = narray[0]->parent;

  /* the root node has no neigbors in its patch */
  if(!parent) return;

  /* fill in neighbor info for the 8 in narray */
  for(ind=0; ind<8; ind++)
  {
    int fs,fo;           /* sibling and other face index */
    tNode *node = narray[ind];
    tNode *parentnb;     /* parent's neighbor */
    tNode *othernb;      /* other neighbor */
    int ijk = node->ijk;
    /* set node's i,j,k */
    int k = kOfInd_n(ijk, ns);
    int j = jOfInd_n_k(ijk, ns,k);
    int i = iOfInd_n_jk(ijk, ns,j,k);
    /* set index of neighbor on same level */
    int inb = i^1;   /* i index of neighbor in same level */
    int jnb = j^1;   /* j index of neighbor in same level */
    int knb = k^1;   /* k index of neighbor in same level */

    /* X-dir: nb[0] and nb[1] */
    fo = i;   /* face index in direction of other neighbor */
    fs = inb; /* face index in direction of sibling */
    /* set sibling neighbor */
    node->nb[fs] = parent->child[Ind_n(inb,j,k, ns)];
    /* set other neighbor at same level if it exists */
    parentnb = parent->nb[fo];
    if(parentnb)
    {
      othernb = parentnb->child[Ind_n(inb,j,k, ns)];
      node->nb[fo] = connect ? othernb : NULL;
      if(othernb) othernb->nb[fs] = connect ? node : NULL;
    }
    else
      node->nb[fo] = NULL;

    /* Y-dir: nb[2] and nb[3] */
    fo = 2 + j;   /* face index in direction of other neighbor */
    fs = 2 + jnb; /* face index in direction of sibling */
    /* set sibling neighbor */
    node->nb[fs] = parent->child[Ind_n(i,jnb,k, ns)];
    /* set other neighbor at same level if it exists */
    parentnb = parent->nb[fo];
    if(parentnb)
    {
      othernb = parentnb->child[Ind_n(i,jnb,k, ns)];
      node->nb[fo] = connect ? othernb : NULL;
      if(othernb) othernb->nb[fs] = connect ? node : NULL;
    }
    else
      node->nb[fo] = NULL;

    /* Z-dir: nb[4] and nb[5] */
    fo = 4 + k;   /* face index in direction of other neighbor */
    fs = 4 + knb; /* face index in direction of sibling */
    /* set sibling neighbor */
    node->nb[fs] = parent->child[Ind_n(i,j,knb, ns)];
    /* set other neighbor at same level if it exists */
    parentnb = parent->nb[fo];
    if(parentnb)
    {
      othernb = parentnb->child[Ind_n(i,j,knb, ns)];
      node->nb[fo] = connect ? othernb : NULL;
      if(othernb) othernb->nb[fs] = connect ? node : NULL;
    }
    else
      node->nb[fo] = NULL;
  }
}


/* enter neighbor info as far as the 8 children of one parent are concerned */
/* this operates on a node array indexed by ijk */
void connect8_siblings(tNode *narray[8])
{
  int ijk;

  errorexit("connect8_siblings is not needed because it does the same as "
            "connect8_with_neighbors");
  errorexit("this function needs to be tested");

  /* fill in neighbor info, as far as these 8 are concerned */
  for(ijk=0; ijk<8; ijk++)
  {
    tNode *node = narray[ijk];
    switch(node->ijk)
    {
    case 0: // i,j,k=0,0,0
        node->nb[1] = narray[1]; // neig. in +X has i,j,k=1,0,0: ijk=1
        node->nb[3] = narray[2]; 
        node->nb[5] = narray[4];
        break;
    case 1: // i,j,k=1,0,0
        node->nb[0] = narray[0]; // neig. in -X has i,j,k=0,0,0: ijk=0
        node->nb[3] = narray[3];
        node->nb[5] = narray[5];
        break;
    case 2: // i,j,k=0,1,0
        node->nb[1] = narray[3];
        node->nb[2] = narray[0];
        node->nb[5] = narray[6];
        break;
    case 3: // i,j,k=1,1,0
        node->nb[0] = narray[2]; // neig. in -X has i,j,k=0,1,0: ijk=2
        node->nb[2] = narray[1]; // neig. in -Y has i,j,k=1,0,0: ijk=1
        node->nb[5] = narray[7]; // neig. in +Z has i,j,k=1,1,1: ijk=7
        break;
    case 4: // i,j,k=0,0,1
        node->nb[1] = narray[5];
        node->nb[3] = narray[6];
        node->nb[4] = narray[0];
        break;
    case 5: // i,j,k=1,0,1
        node->nb[0] = narray[4];
        node->nb[3] = narray[7];
        node->nb[4] = narray[1];
        break;
    case 6:
        node->nb[1] = narray[7];
        node->nb[2] = narray[4];
        node->nb[4] = narray[2];
        break;
    case 7:
        node->nb[0] = narray[6];
        node->nb[2] = narray[5];
        node->nb[4] = narray[3];
        break;
    }
  }
}

/* Old version:
   construct a unique number that describes node location in patch: 
   e.g. 1347 in octal is node on level 3 that has
   ijk=3 on l3, ijk=4 on l2, ijk=7 on l1, the leading 1 is only there to
   distinguish case like: 1000 (ijk=0 on l1-3) from 100 (ijk=0 on l1-2). */
unsigned long node_location0(tNode *node)
{
  tNode *anc;
  long loc;
  
  for(loc=1, anc=node; anc->parent; anc = anc->parent)
  {
    loc = loc<<3;
    loc += anc->ijk;
  }
  return loc;
}
void node_location0_str(tNode *node, char *s, int slen)
{
  int n, i;
  unsigned long loc = node_location0(node);
  unsigned long lloc;

  /* n is how mant times n we can shift loc right (by 3) so that is is zero */
  for(n=0, lloc=loc;  lloc;  n++, lloc=lloc>>3) ;

  /* shift loc and write into string s */
  for(i=2; i<=n && i<=slen; i++)
  {
    lloc = loc>>(3*(n-i));
    s[i-2] = (lloc & 7) + '0';
  }
  s[i-2]=0;
}

/* construct a unique string that describes node location in patch: 
   e.g. 743 in octal is node on level 3 that has
    ijk=7 on l1,  ijk=4 on l2,  ijk=3 on l3   */
char *node_location_str(tNode *node, char *s, int slen)
{
  tNode *anc;
  int l;

  if(!node)
  {
    if(slen>0) s[0] = 0;
    return s;
  }

  l = node->l;

  if(slen<=l) errorexit("slen is not big enough");
  s[l--] = 0;
  for(anc=node; anc->parent; anc = anc->parent)
    s[l--] = anc->ijk + '0';

  return s;
}

/* convert string from node_location_str into a long int */
long node_location(tNode *node)
{
  tNode *anc;
  int il;
  long loc = 0;

  if(!node) return 0;

  for(il=0, anc=node; anc->parent; anc = anc->parent)
  {
    loc |= (anc->ijk)<<(il*3);
    il++;
  }
  /* add a leading 1 for the root node */
  loc |= (1)<<(il*3);

  /* in case overflow occurred in bit shifts */
  if(loc<0) loc=0;

  return loc;
}

/* use node_location_str to make a unique node name that also contains the
   patch number */
char *nodename(tNode *node, char *s, int slen)
{
  char loc[100];
  if(node)
  {
    node_location_str(node, loc,99);
    snprintf(s,slen, "%d_%s", node->pat->p, loc);
  }
  else
    snprintf(s,slen, "-");
  return s;
}

/* is node on a face? */
int node_is_at_face(tNode *node, int face)
{
  int ns[] = {2,2,2};
  int ijk = node->ijk;
  /* set node's i,j,k */
  int k = kOfInd_n(ijk, ns);
  int j = jOfInd_n_k(ijk, ns,k);
  int i = iOfInd_n_jk(ijk, ns,j,k);

  switch(face)
  {
  case 0:  return i^1;
  case 1:  return i;
  case 2:  return j^1;
  case 3:  return j;
  case 4:  return k^1;
  case 5:  return k;
  default: errorexit("face must be 0,1,2,3,4,5");
  }
  return -1;
}

/* count children of a node, should be 0 or 8 */
int count_children(tNode *node)
{
  int i, nc;
  for(nc=0, i=0; i<8; i++) if(node->child[i]) nc++;
  return nc;
}

/* return all leaf descendants along face, in ndescends we return how many
   times we have descended, the list it returns has to be freed
   by caller */
/* DO NOT USE THIS. Instead use: leafdescendants_along_face */
tNlist *ldescendants_along_face(tNlist *nl, int face, int *ndescends)
{
  tNlist *elem;
  tNlist *nl2 = copy_of_nodelist(nl);
  tNlist *children;
  tNode *node, *child;

  *ndescends = 0; /* number of replcements made in nl2 */

  /* loop over nl2 and make replacements with children in face */
  fornodelist(nl2, elem)
  {
  nl2_loop_start:
    children = NULL;
    node = elem->node;

    /* make list of children */    
    if(node->child[0])
    {
      int i;

      children = NULL;
      for(i=0; i<8; i++)
      {
        child = node->child[i];
        if(node_is_at_face(child, face))
          children = addnode_to_nodelist_after(children, child);
      }
      /* insert children into nl2, replacing parent in elem */
      elem = replace1_in_nodelist(elem, children, 0);
      nl2 = elem;
      (*ndescends)++;
      if(elem) goto nl2_loop_start;
      else     errorexit("elem should not be NULL after inserting children");
    }
  }
  /* results */
  nl2 = first_nodelist(nl2);

  return nl2;
}

/* return all leaf descendants of node along face. Call it like this:
   leafdesc = leafdescendants_along_face(node, f, NULL);
   the list leafdesc it returns has to be freed by caller */
tNlist *leafdescendants_along_face(tNode *node, int face, tNlist *leafdesc)
{
  tNode *child;
  int i;
  tNode *locker = get_node_nc_lock(node);

  /* check if node has children but only if it is finished (unlocked) */
  if( (locker==NULL) && (node->child[0]) )
  {
    for(i=0; i<8; i++)
    {
      child = node->child[i];
      if(node_is_at_face(child, face))
        leafdesc = leafdescendants_along_face(child, face, leafdesc);
    }
  }
  else /* no children, so add the node to leaf descendants list */
    leafdesc = addnode_to_nodelist_after(leafdesc, node);

  return leafdesc;
}

/*************************************************************************/
/* functions to set nfaces and fnb */
/*************************************************************************/

/* find nfaces within this patch, this adds these nfaces to the node and
   its neighbors */
int add_nfaces_within_patch(tNode *node, int face)
{
  tNlist *nblist, *elem;
  tNode *anc, *nb;
  int nc, nnb;
  int nbface;

  /* no neighb. if on patch face */
  if(node->patface[face]) return 0;

  /* find neighbor at same level or on lower level
     note: root node has no patch neighbors */
  for(nb=NULL, anc=node; anc->parent; anc=anc->parent)
  {
    nb = anc->nb[face];
    if(nb) break;
  }
  /* if no neighbor at all is found return just 0 */
  if(!nb) return 0;

  /* face where neighbors are */
  nbface = face^1;

  /* so now we have a neighbor, but is it childless? */
  nc = count_children(nb);
  if(nc==0) /* neighbor has 0 children */
  {
    add_nface(node, face, nb, nbface);
    return 1;  /* there is only one neighbor */
  }
  if(nc!=8) errorexiti("nb has %d children, not 8!!!", nc);

  /* ok so this neighbor has 8 children, who also may have children */
  nblist = leafdescendants_along_face(nb, nbface, NULL);
  nblist = first_nodelist(nblist);

  /* add all in nblist as nfaces */
  nnb = 0;
  fornodelist(nblist, elem)
  {
    nb = elem->node;
    add_nface(node, face, nb, nbface);
    nnb++; /* count neighbors */
  }

  /* free node lists */
  free_nodelist(nblist);

  return nnb;
}

/* find nfaces outside this patch (using bfaces), this adds these nfaces to
   the node and its neighbors */
int add_nfaces_outside_patch(tNode *node, int face)
{
  tPat *pat = node->pat;
  tBface *bface;
  tNlist *nbl, *nblist1, *elem;
  tNode *nb;
  int nc, nb_f;
  int nnb = 0;   /* number of nfaces added */

  /* no outside neighb. if not on patch face */
  if(!node->patface[face]) return 0;

  //PRF;printf(":\n");

  /* loop over all bfaces on face and find nb */
  forbfacesonface(pat, face, bface)
  {
    tBface *obface = bface->obface;
    int touch;

    /* do nothing if no other patch face */
    if(!obface) continue;

    /* root node in other patch */
    nb = obface->pat->rnode;
    nb_f = obface->f;

    /* so now we have a neighbor, but is it childless? */
    nc = count_children(nb);
    if(nc==0) /* neighbor has 0 children */
    {
      nblist1 = alloc_nodelist(nb);
    }
    else
    {
      if(nc!=8) errorexiti("nb has %d children, not 8!!!", nc);

      /* find nblist1 with all leaves on face nb_f */
      nblist1 = leafdescendants_along_face(nb, nb_f, NULL);
    }

    /* beginning of nblist1 */
    nbl = first_nodelist(nblist1);

//char s[100];
//if(strcmp(nodename(node,s,99),"0_1")==0 && face==1)
//{
//PRFs(":\n");printnodelist(nbl);
//}

    /* go over nbl and remove all who do not have common face points
       with the node */
    nblist1 = NULL;
    fornodelist(nbl, elem)
    {
    nbl_loop_start:

      /* get neigh. and check if node and nb have common points */
      nb = elem->node;
      touch = common_facepoints(node,face, nb,nb_f);
      if(touch)
      {
        nblist1 = elem; /* save elem that touches our node */
        continue;
      }

//if(strcmp(nodename(node,s,99),"0_1")==0 && face==1)
//{
//PRFs(": remove:\n");printnodelist(elem);
//}

      /* remove nb=elem->node from nbl */
      elem = remove1_in_nodelist(elem, 1); /* now elem has the next one */
      if(elem) goto nbl_loop_start;
      else     break;
    }

//if(strcmp(nodename(node,s,99),"0_1")==0 && face==1)
//{
//PRFs(": final:\n");printnodelist(nblist1);
//}

    /* rewind nblist1 so that the fornodelist loop blow works */
    nblist1 = first_nodelist(nblist1);

    /* add all in nblist1 as nfaces */
    fornodelist(nblist1, elem)
    {
      nb = elem->node;
      add_nface(node, face, nb, nb_f);
      nnb++; /* count neighbors */
    }

    /* free node lists */
    free_nodelist(nblist1);
  }

//  char s[100];
//  PRF;printf(": %s f%d\n", nodename(node,s,99), face);
//  printnodelist(nblist1); fflush(stdout);

  return nnb;
}

/* all all nfaces on one node face */
int add_nfaces_at_nodeface(tNode *node, int face)
{
  int num;

  num  = add_nfaces_within_patch(node, face);
  num += add_nfaces_outside_patch(node, face);
  return num;
}

/* update all nfaces of a node. This also updates the neighbor faces */
void update_node_nfaces(tNode *node)
{
  tNode *parent = node->parent;
  int face;

  /* remove all nfaces of node and its parent */
  remove_all_nfaces(node);
  remove_all_nfaces(parent);

  /* go over all faces and add nefaces */
  for(face=0; face<6; face++)
    add_nfaces_at_nodeface(node, face);
}


/* set surface neigbhor list in node on face:
   save neighbors on all 6 faces in list node->fnb[face] */
void set_node_fnb_from_nfaces(tNode *node, int face)
{
  tNface *nface;
  int nfnb, ni;

  /* count number of nb */
  nfnb = 0;
  for(nface=node->nfaces[face]; nface; nface=nface->next) nfnb++;

  /* first free and then allocate room for neighbors */
  free(node->fnb[face]);
  node->fnb[face] = NULL;
  if(nfnb)
    node->fnb[face] = calloc(nfnb, sizeof(node->fnb[face][0]));

  /* now set fnb on face */
  //PRF;printf(" face=%d\n", face);
  //printnfaces_on_f(node, face);
  node->nfnb[face] = nfnb;
  ni = 0;
  for(nface=node->nfaces[face]; nface; nface=nface->next)
  {
    //printnface(nface);
    node->fnb[face][ni] = nface->onface->node;
    ni++;
  }
}

/* initialize the 6 surface neigbhor lists in node and in the faces of
   other nodes that touch it */
//void update_node_and_neighbors_fnb__new(tNode *node)
void update_node_and_neighbors_nfaces_fnb(tNode *node)
{
  DECL_MESH_MUTEX(node, mutex);
  int face;

  /* set nfaces on node and neighbors */
  update_node_nfaces(node);

//char s[100];
//if(strcmp(nodename(node,s,99),"0_1")==0)
//{
//PRF;printnfaces(node);
//}

  /* update fnb on this node */
  for(face=0; face<6; face++)
    set_node_fnb_from_nfaces(node, face);

  /* now update fnb on all its neighbors */
  MUTEX_LOCK(mutex); /* use mutex, as fnb is used in locking nodes */
  for(face=0; face<6; face++)
  {
    tNface *nface;

    for(nface=node->nfaces[face]; nface; nface=nface->next)
    {
      tNface *onface = nface->onface;
      tNode *nb = onface->node;
      int nb_f  = onface->f;

      set_node_fnb_from_nfaces(nb, nb_f);
    }
  }
  MUTEX_UNLOCK(mutex);
}

// Not sure if this function is useful ...
/* initialize the 6 surface neigbhor lists in node and in the faces of
   other nodes that touch it, but set only node's fnb */
void update_node_nfaces_fnb(tNode *node)
{
  int face;

  /* set nfaces on node and neighbors */
  update_node_nfaces(node);

  /* update fnb on this node */
  for(face=0; face<6; face++)
    set_node_fnb_from_nfaces(node, face);
}

/* init surface neighbor list for all root nodes in the mesh */
void update_all_rnode_nfaces_fnb(tMesh *mesh)
{
  int p;
  forpatches(mesh, p)
    update_node_and_neighbors_nfaces_fnb(mesh->pat[p]->rnode);
}

/*************************************************************************/
/* functions to set fnb without nfaces */
/*************************************************************************/

/* find leaf node neighbors within this patch, this allocates the nodelist
   containing them, which has to be freed by caller */
tNlist *make_patch_neighbor_list(tNode *node, int face)
{
  tNlist *nblist;
  tNode *anc, *nb;
  int nc;
  int nbface;

  /* no neighb. if on patch face */
  if(node->patface[face]) return NULL;

  /* find neighbor at same level or on lower level
     note: root node has no patch neighbors */
  for(nb=NULL, anc=node; anc->parent; anc=anc->parent)
  {
    nb = anc->nb[face];
    if(nb) break;
  }
  /* if no neighbor at all is found return just 0 */
  if(!nb) return NULL;

  /* so now we have a neighbor, but is it childless? */
  nc = count_children(nb);
  if(nc==0) /* neighbor has 0 children */
  {
    nblist = alloc_nodelist(nb);
    return nblist;  /* there is only one neighbor */
  }
  if(nc!=8) errorexiti("nb has %d children, not 8!!!", nc);

  /* ok so this neighbor has 8 children, who also may have children */
  nbface = face^1; /* face where neighbors are */
  nblist = leafdescendants_along_face(nb, nbface, NULL);

  return nblist;
}

/* find leaf node neighbors outside this patch (using bfaces) */
/* this allocates the nodelist containing them, which has to be freed by
   caller */
tNlist *make_outside_neighbor_list(tNode *node, int face)
{
  tPat *pat = node->pat;
  tBface *bface;
  tNlist *nblist;
  tNlist *nbl, *nblist1, *elem;
  tNode *nb;
  int nc, nb_f;

  /* no outside neighb. if not on patch face */
  if(!node->patface[face]) return NULL;

  //PRF;printf(":\n");

  /* loop over all bfaces on face and find nb */
  nblist = NULL;
  forbfacesonface(pat, face, bface)
  {
    tBface *obface = bface->obface;
    int touch;

    /* do nothing if no other patch face */
    if(!obface) continue;

    /* root node in other patch */
    nb = obface->pat->rnode;
    nb_f = obface->f;

    /* so now we have a neighbor, but is it childless? */
    nc = count_children(nb);
    if(nc==0) /* neighbor has 0 children */
    {
      nblist1 = alloc_nodelist(nb);
    }
    else
    {
      if(nc!=8) errorexiti("nb has %d children, not 8!!!", nc);

      /* find nblist1 with all leaves on face nb_f */
      nblist1 = leafdescendants_along_face(nb, nb_f, NULL);
    }

    /* beginning of nblist1 */
    nbl = first_nodelist(nblist1);

    /* go over nbl and remove all whose do not have common face points
       with the node */
    nblist1 = NULL;
    fornodelist(nbl, elem)
    {
    nbl_loop_start:

      /* get neigh. and check if node and nb have common points */
      nb = elem->node;
      touch = common_facepoints(node,face, nb,nb_f);
      if(touch)
      {
        nblist1 = elem; /* save elem that touches our node */
        continue;
      }

      /* remove nb=elem->node from nbl */
      elem = remove1_in_nodelist(elem, 1); /* now elem has the next one */
      if(elem) goto nbl_loop_start;
      else     break;
    }

    /* add nblist1 to nblist */
    nblist = insertnodelist_into_nodelist_after(nblist, nblist1);
  }

  return nblist;
}


/* find all leaf node neighbors of node in mesh, this allocates the
   nodelist containing them, which has to be freed by caller */
tNlist *make_mesh_neighbor_list(tNode *node, int face)
{
  tNlist *nblist, *nblist1;

  /* first find all neighbors inside patch */
  nblist1 = make_patch_neighbor_list(node, face);

  /* find leaf node neighbors outside this patch (using bfaces) */
  nblist = make_outside_neighbor_list(node, face);

  /* combine result with nblist1, i.e add nblist to nblist1 */
  nblist = insertnodelist_into_nodelist_after(nblist1, nblist);

  nblist = first_nodelist(nblist);

  return nblist;
}


/* initialize surface neigbhor list in node:
   save neighbors on all 6 faces in list node->fnb[face] */
void update_node_fnb_only(tNode *node)
{
  tNlist *nblist, *elem;
  int face, nfnb, ni;

  for(face=0; face<6; face++)
  {
    /* find neighbors */
    nblist = make_mesh_neighbor_list(node, face);
    nfnb = count_elements_nodelist(nblist);
    node->nfnb[face] = nfnb;

    /* first free and then allocate room for neighbors */
    free(node->fnb[face]);
    node->fnb[face] = NULL;
    if(nfnb)
      node->fnb[face] = calloc(nfnb, sizeof(node->fnb[face][0]));

    /* add neighbors to node */
    ni = 0;
    fornodelist(nblist, elem)
    {
      /* save this neighbor */
      node->fnb[face][ni] = elem->node;
      ni++; /* inc node counter */
    }
    free_nodelist(nblist);
  }
}

/* remove nb on face f with index ni from node->fnb */
void remove_node_fnb_only(tNode *node, int f, int ni)
{
  int nfnb = node->nfnb[f];
  int i;

  for(i=ni+1; i<nfnb; i++) node->fnb[f][i-1] = node->fnb[f][i];
  node->nfnb[f] = nfnb - 1;
}

/* remove nb on face f with index ni from node->fnb */
void remove_unpaired_node_fnbs_only(tNode *node)
{
  int f, ni;

  for(f=0; f<6; f++)
    for(ni=0; ni<node->nfnb[f]; ni++)
    {
      tNode *nb = node->fnb[f][ni];
      int found, nb_f, nb_ni;

      /* remove nb with index ni if it does not also have node as neighb. */
      found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
      if(!found)
      {
        remove_node_fnb_only(node, f, ni);
        ni--;

if(node->nid==72)
{
Yo(3.0);
printnode(node);
printnode(nb);
printf("remove:"); printnd(nb);
}
      }
    }
}


/* same as update_node_fnb, but do it for neighbors as well */
//void update_node_and_neighbors_fnb__old(tNode *node)
void update_node_and_neighbors_fnb_only(tNode *node)
{
  int face;

  /* update this node */
  update_node_fnb_only(node);

  /* now update all its neighbors */
  for(face=0; face<6; face++)
  {
    int ni;
    for(ni=0; ni<node->nfnb[face]; ni++)
    {
      tNode *nb = node->fnb[face][ni];
      update_node_fnb_only(nb);
    }
  }
  // no unpaired nodes should exist!!!
  // /* remove all that is not paired on this node */
  // remove_unpaired_node_fnbs(node);

  // /* now  remove all that is not paired on neighbors */
  // for(face=0; face<6; face++)
  // {
  //   int ni;
  //   for(ni=0; ni<node->nfnb[face]; ni++)
  //   {
  //     tNode *nb = node->fnb[face][ni];
  //     remove_unpaired_node_fnbs(nb);
  //   }
  // }
}

/* init surface neighbor list for all root nodes in the mesh */
void update_all_rnode_fnb_only(tMesh *mesh)
{
  int p;
  forpatches(mesh, p)
    update_node_and_neighbors_fnb_only(mesh->pat[p]->rnode);
}

/* find node facenb in the node->fnb lists, returns 1  if it is there
   returns face and nb index in vars: face and ni */
int locate_facenb_in_fnbs(tNode *node, tNode *facenb, int *face, int *ni)
{
  int found = 0;
  int f;
  for(f=0; f<6; f++)
  {
    int nfnb = node->nfnb[f];
    int i;
    for(i=0; i<nfnb; i++)
    {
      if(node->fnb[f][i] == facenb)
      {
        found = 1;
        *face = f;
        *ni   = i;
        break;
      }
    }
    if(found) break;
  }
  return found;
}


/*******************************************************************/
/* nc_lock helper functions */
/*******************************************************************/

/* return first locker on node and neighbors */
tNode *node_and_fnbs_nc_locked(tNode *node)
{
  int face;

  /* is this node locked */
  if(node->nc_lock) return node->nc_lock;

  /* are any neighbors locked */
  for(face=0; face<6; face++)
  {
    int ni;
    for(ni=0; ni<node->nfnb[face]; ni++)
    {
      tNode *nb = node->fnb[face][ni];
      if(nb->nc_lock) return nb->nc_lock;
    }
  }

  return NULL;
}

/* set nc_lock on node and its fnb neighbors to nc_lock = locker */
void set_nc_lock_on_node_and_fnbs(tNode *node, tNode *locker)
{
  int face;

  /* mark node as locked by locker */
  node->nc_lock = locker;
  //char ns[100], ls[100];
  //PRF;printf(": %s locked by %s\n",
  //           nodename(node,ns,99), nodename(locker,ls,99));
  //fflush(stdout);

  /* mark neighbors as locked by locker */
  for(face=0; face<6; face++)
  {
    int ni;
    for(ni=0; ni<node->nfnb[face]; ni++)
    {
      tNode *nb = node->fnb[face][ni];
      nb->nc_lock = locker;
      //PRF;printf(":   %s locked by %s\n",
      //           nodename(nb,ns,99), nodename(locker,ls,99));
      //fflush(stdout);
    }
  }
}

/* return first locker on parent or neighbors of a nodearray */
tNode *parent_and_fnbs_nc_locked(tNode *narray[8])
{
  tNode *parent = narray[0]->parent;
  int ns[] = {2,2,2};
  int ind;

  /* check parent */
  if(parent)
    if(parent->nc_lock) return parent->nc_lock;

  /* check fnbs */
  for(ind=0; ind<8; ind++)
  {
    int ni;
    int fo;   /* other face index */
    tNode *node = narray[ind];
    int ijk = node->ijk;
    /* set node's i,j,k */
    int k = kOfInd_n(ijk, ns);
    int j = jOfInd_n_k(ijk, ns,k);
    int i = iOfInd_n_jk(ijk, ns,j,k);

    /* return locker of 1st nodearray neighbor that is locked */
    /* X-dir: */
    fo = i;   /* face index in direction of other neighbor */
    for(ni=0; ni<node->nfnb[fo]; ni++)
    {
      tNode *nb = node->fnb[fo][ni];
      if(nb->nc_lock) return nb->nc_lock;
    }

    /* Y-dir: */
    fo = 2 + j;   /* face index in direction of other neighbor */
    for(ni=0; ni<node->nfnb[fo]; ni++)
    {
      tNode *nb = node->fnb[fo][ni];
      if(nb->nc_lock) return nb->nc_lock;
    }

    /* Z-dir: */
    fo = 4 + k;   /* face index in direction of other neighbor */
    for(ni=0; ni<node->nfnb[fo]; ni++)
    {
      tNode *nb = node->fnb[fo][ni];
      if(nb->nc_lock) return nb->nc_lock;
    }
  }

  return NULL;
}

/* set nc_lock on parent and neighbors of narray to nc_lock = locker */
void set_nc_lock_on_parent_and_fnbs_of_nodearray(tNode *narray[8],
                                                 tNode *locker)
{
  tNode *parent = narray[0]->parent;
  int ns[] = {2,2,2};
  int ind;

  if(parent) parent->nc_lock = locker;

  for(ind=0; ind<8; ind++)
  {
    int ni;
    int fo;   /* other face index */
    tNode *node = narray[ind];
    int ijk = node->ijk;
    /* set node's i,j,k */
    int k = kOfInd_n(ijk, ns);
    int j = jOfInd_n_k(ijk, ns,k);
    int i = iOfInd_n_jk(ijk, ns,j,k);

    /* mark nodearray neighbors as locked by locker */
    /* X-dir: */
    fo = i;   /* face index in direction of other neighbor */
    for(ni=0; ni<node->nfnb[fo]; ni++)
    {
      tNode *nb = node->fnb[fo][ni];
      nb->nc_lock = locker;
    }

    /* Y-dir: */
    fo = 2 + j;   /* face index in direction of other neighbor */
    for(ni=0; ni<node->nfnb[fo]; ni++)
    {
      tNode *nb = node->fnb[fo][ni];
      nb->nc_lock = locker;
    }

    /* Z-dir: */
    fo = 4 + k;   /* face index in direction of other neighbor */
    for(ni=0; ni<node->nfnb[fo]; ni++)
    {
      tNode *nb = node->fnb[fo][ni];
      nb->nc_lock = locker;
    }
  }
}

/*******************************************************************/
/* lock or unlock */
/*******************************************************************/

/* wait until exclusive lock is acquired by node */
void node_and_fnbs_lock(tNode *node)
{
  DECL_MESH_MUTEX(node, mutex);

  MUTEX_LOCK(mutex);
  {
    /* if node is locked by another node or not locked at all */
    if(node->nc_lock != node)
    {
      /* wait until node and nbs are no longer locked */
      while(node_and_fnbs_nc_locked(node))
      {
        /* if node cannot get the nc_lock release mutex and yield */
        MUTEX_UNLOCK(mutex);
        TASK_YIELD;

        /* get mutex back */
        MUTEX_LOCK(mutex);
      }
    }
    /* lock this node and its nbs */
    set_nc_lock_on_node_and_fnbs(node, node);
  }
  MUTEX_UNLOCK(mutex);
}

/* unlock node and its neighbors if node locked it */
void node_and_fnbs_unlock(tNode *node)
{
  DECL_MESH_MUTEX(node, mutex);

  MUTEX_LOCK(mutex);
  {
    /* unlock this node and its nbs */
    if(node->nc_lock == node)
      set_nc_lock_on_node_and_fnbs(node, NULL);
  }
  MUTEX_UNLOCK(mutex);
}

/* wait until exclusive lock is acquired on fnbs of narray */
void parent_and_fnbs_lock(tNode *narray[8], tNode *locker)
{
  DECL_MESH_MUTEX(locker, mutex);

  MUTEX_LOCK(mutex);
  {
    if(parent_and_fnbs_nc_locked(narray) != locker)
    {
      /* wait until nbs are no longer locked */
      while(parent_and_fnbs_nc_locked(narray))
      {
        /* if locker cannot get the nc_lock release mutex and yield */
        MUTEX_UNLOCK(mutex);
        TASK_YIELD;

        /* get mutex back */
        MUTEX_LOCK(mutex);
      }
    }
    /* lock nbs */
    set_nc_lock_on_parent_and_fnbs_of_nodearray(narray, locker);
  }
  MUTEX_UNLOCK(mutex);
}

/* unlock fnbs if locker locked it */
void parent_and_fnbs_unlock(tNode *narray[8], tNode *locker)
{
  DECL_MESH_MUTEX(locker, mutex);

  MUTEX_LOCK(mutex);
  {
    /* unlock fnbs */
    if(parent_and_fnbs_nc_locked(narray) == locker)
      set_nc_lock_on_parent_and_fnbs_of_nodearray(narray, NULL);
  }
  MUTEX_UNLOCK(mutex);
}


/*******************************************************************/
/* check lock */
/*******************************************************************/

/* return locker on node, protected by mutex */
tNode *get_node_nc_lock(tNode *node)
{
  tNode *locker;
  DECL_MESH_MUTEX(node, mutex);

  MUTEX_LOCK(mutex);
  {
    locker = node->nc_lock;
  }
  MUTEX_UNLOCK(mutex);

  return locker;
}



/*******************************************************************/
/* Everything below this line is untested and may not work
   Probably it should be removed !!!!! */
/*******************************************************************/

/* Go up in ancestors until we reach one that has either a different
   i, j, or k. We look for the difference in direction dir=0,1,2 and return
   the ancestor's up this point in anclist. */
tNlist *ancestors_alongBoundary(tNode *node, int dir)
{
  tNode *anc; /* ancestor */
  int ijk = node->ijk;
  int ns[] = {2,2,2};
  int k = kOfInd_n(ijk, ns);       /* set node's i,j,k */
  int j = jOfInd_n_k(ijk, ns,k);
  int i = iOfInd_n_jk(ijk, ns,j,k);
  tNlist *anclist = alloc_nodelist(node);

  /* loop up ancestors */
  for(anc=node->parent; anc; anc=anc->parent)
  {
    /* get ancestor's i,j,k */
    int aijk = anc->ijk;
    int ak = kOfInd_n(aijk, ns);
    int aj = jOfInd_n_k(aijk, ns,ak);
    int ai = iOfInd_n_jk(aijk, ns,aj,ak);

    anclist = addnode_to_nodelist_before(anclist, anc);

    if(dir==0)
    {
      if(ai!=i) break;
    }
    else if(dir==1)
    {
      if(aj!=j) break;
    }
    else if(dir==2)
    {
      if(ak!=k) break;
    }
    else
      errorexit("dir needs to be 0,1,2");
  }
  /* now also add parent of last real anc */
  if(anc)
    addnode_to_nodelist_before(anclist, anc->parent);
  else
    addnode_to_nodelist_before(anclist, NULL);

  /* return the list of ancestors */
  return anclist;
}
/* go back down again to find neighbor */
void descend_alongBoundary(tNlist *anclist, int dir, tNode *nb)
{
  int i;
  int ns[] = {2,2,2};
  tNode *anc0 = anclist->node;
  tNode *anc1;
  tNode *a0nb;
  tNode *a1nb;

  if(!anc0) return;

  a0nb = anc0;
  fornodes(anclist->next, anc1)
  {
    /* get descendant's i,j,k */
    int ijk1 = anc1->ijk;
    int k1 = kOfInd_n(ijk1, ns);
    int j1 = jOfInd_n_k(ijk1, ns,k1);
    //int i1 = iOfInd_n_jk(ijk1, ns,j1,k1);

    // X dir only!!!
    i=0; // need to set i to in anc1 ???
    a1nb = a0nb->child[Ind_n(i,j1,k1,  ns)];
    if(a1nb->leaf)
      break; // nb = a1nb;

    a0nb = a1nb;
  } endfornodes;
  nb = a1nb;
}
