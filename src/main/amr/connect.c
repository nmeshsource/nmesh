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
  int l = node->l;

  if(slen<=l) errorexit("slen is not big enough");
  s[l--] = 0;
  for(anc=node; anc->parent; anc = anc->parent)
    s[l--] = anc->ijk + '0';

  return s;
}

/* use node_location_str to make a unique node name that also contains the
   patch number */
char *nodename(tNode *node, char *s, int slen)
{
  char loc[100];
  node_location_str(node, loc,99);
  snprintf(s,slen, "p%d_%s", node->pat->p, loc);
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
      elem = replace1_in_nodelist(elem, children);
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

  /* check if node has children */
  if(node->child[0])
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
  double brct[4]; // bound. rect. of node

  /* no outside neighb. if not on patch face */
  if(!node->patface[face]) return NULL;

  //PRF;printf(":\n");

  /* set bound. rect. of node */
  brct_nodeface(node, face/2, brct);

  /* loop over all bfaces on face and find nb */
  nblist = NULL;
  forbfacesonface(pat, face, bface)
  {
    tBface *obface = bface->obface;
    double nbrct[4]; // bound. rect. of nb
    int problem;     // is set to 1 if there is a problem finding nbrct
    double irct[4];
    int isec;

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
/*
if(node->nid==28)
{
Yo(1);
printnodelist(nbl);
}
*/
    /* go over nbl and remove all whose face does not intersect
       with the node bounding rectangle brct */
    nblist1 = NULL;
    fornodelist(nbl, elem)
    {
    nbl_loop_start:
/*
if(node->nid==28)
{
Yo(2.1);
printnodelist(nblist1);
//printnode(node);
//printnode(nb);
}
*/
      /* get neigh. bound. rect. in its own X coords */
      nb = elem->node;
      brct_nodeface(nb, nb_f/2, nbrct);
/*
if(node->nid==28)
{
Yo(2.12);
prbbox(brct,2);
prbbox(nbrct,2);
printf("\n");
}
*/
      /* transform nbrct from nb coords to node coords */
      problem = brctpat2_of_brctpat1(nb->pat, nb_f, nbrct,
                                     node->pat, face, nbrct);

      /* does brct intersect nbrct? */
      isec = intersection_brct1_brct2(brct, nbrct, irct);
      if(isec && !problem)
      {
        nblist1 = elem; /* save elem that touches our node */
        continue;
      }
/*
if(node->nid==28)
{
Yo(2.2);
printnodelist(nblist1);
printnodelist(elem);
prbbox(brct,2);
prbbox(nbrct,2);
prbbox(irct,2);
printf("isec=%d problem=%d\n", isec, problem);
}
*/
      /* remove nb=elem->node from nbl */
      elem = remove1_in_nodelist(elem, 1); /* now elem has the next one */
      if(elem) goto nbl_loop_start;
      else     break;
    }
/*
if(node->nid==28)
{
Yo(2.9);
printnodelist(nblist1);
}
*/
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
void update_node_fnb(tNode *node)
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

/* same as update_node_fnb, but do it for neighbors as well */
void update_node_and_neighbors_fnb(tNode *node)
{
  int face;

  /* update this node */
  update_node_fnb(node);

  /* now update all its neighbors */
  for(face=0; face<6; face++)
  {
    int ni;
    for(ni=0; ni<node->nfnb[face]; ni++)
    {
      tNode *nb = node->fnb[face][ni];
      update_node_fnb(nb);
    }
  }
}

/* init surface neighbor list for all root nodes in th e mesh */
void update_all_rnode_fnb(tMesh *mesh)
{
  int p;
  forpatches(mesh, p)
    update_node_and_neighbors_fnb(mesh->pat[p]->rnode);
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
