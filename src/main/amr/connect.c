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

  errorexit("connect8_siblings is not neded it does the same as "
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

/* return all descendants along face */
tNlist *all_descendants_along_face(tNlist *nl, int face)
{
  tNlist *elem;
  tNlist *nl2 = copy_of_nodelist(nl);

  fornodelist(nl2, elem)
  {
    tNlist *children = NULL;
    tNode *node = elem->node;
    tNode *child;
    int i;

    /* make list of children */    
    if(node->child[0])
    {
      children = NULL;
      for(i=0; i<8; i++)
      {
        child = node->child[i];
        if(node_is_at_face(child, face))
          children = addnode_to_nodelist_after(children, child);
      }
      /* insert children into nl2, replacing parent in elem */
      elem = first_replace1_in_nodelist(elem, children);
    }
  }
  return first_nodelist(elem);
}


/* find leaf node neighbors within this patch, this allocates the nodelist
   containing them */
tNlist *find_patch_neighbors(tNode *node, int face)
{
  tNlist *nbl;
  tNlist *nblist;
  tNode *anc, *nb;
  int nc;
  int nbface;

  /* empty list */
  nblist = NULL;

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
  if(!nb) return NULL;

  /* so now we have a neighbor, but it it childless? */
  nc = count_children(nb);
  if(nc==0) /* neighbor has 0 children */
  {
    nblist = alloc_nodelist(nb);
    return nblist;  /* there is only one neighbor */
  }
  if(nc!=8) errorexiti("nb has %d children, not 8!!!", nc);

  /* ok so this neighbor has 8 children, who also may have children */
  nbl = alloc_nodelist(nb);
  nbface = face^1; /* face where neighbors are */
  nblist = all_descendants_along_face(nbl, nbface);
  free_nodelist(nbl);

  return nblist;
}

/* find leaf node neighbors outside this patch (using bfaces) */
//... TODO

/* find all neighbors of a leaf node */



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
