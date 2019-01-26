/* connect.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"




// We need a func that just connects one newly added, whether if has
// siblings or not!!!



/* Go up in ancestors until we reach one that has either a different
   i, j, or k. We look for the difference in direction dir=0,1,2 and return
   the ancestor node on &ancestor. */
void ancestor_with_diff_i_j_k(tNode *node, int dir, tNode **ancestor)
{
  int ijk = node->ijk;
  int ns[] = {2,2,2};
  int k = kOfInd_n(ijk, ns);       /* set node's i,j,k */
  int j = jOfInd_n_k(ijk, ns,k);
  int i = iOfInd_n_jk(ijk, ns,j,k);

  /* loop up ancestors */
  for(anc=node->parent; anc; anc=anc->parent)
  {
    /* get ancestor's i,j,k */
    aijk = anc->ijk;
    ak = kOfInd_n(aijk, ns);
    aj = jOfInd_n_k(aijk, ns,ak);
    ai = iOfInd_n_jk(aijk, ns,aj,ak);

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
      errorexit("dir needs to be 0,1,2")
  }
  /* return anc */
  &ancestor = anc;
}


tNlist *ancestors_alongBoundary(tNode *node, int dir)
{
  tNode *anc; /* ancestor */
  int ijk = node->ijk;
  int ns[] = {2,2,2};
  int k = kOfInd_n(ijk, ns);       /* set node's i,j,k */
  int j = jOfInd_n_k(ijk, ns,k);
  int i = iOfInd_n_jk(ijk, ns,j,k);
  tNlist *anclist = alloc_nodelist(tNode *node);

  /* loop up ancestors */
  for(anc=node->parent; anc; anc=anc->parent)
  {
    anclist = addnode_to_nodelist_before(anclist, anc);
    /* get ancestor's i,j,k */
    aijk = anc->ijk;
    ak = kOfInd_n(aijk, ns);
    aj = jOfInd_n_k(aijk, ns,ak);
    ai = iOfInd_n_jk(aijk, ns,aj,ak);

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
      errorexit("dir needs to be 0,1,2")
  }
  /* now also add parent of last real anc */
  if(anc)
    addnode_to_nodelist_before(anclist, anc->parent);
  else
    addnode_to_nodelist_before(anclist, NULL);

  /* return the list of ancestors */
  return anclist;
}


void descend_alongBoundary(tNlist *anclist, int dir, tNode **descendant)
{
  int count;
  int ns[] = {2,2,2};
  tNode *anc0 = anclist->node;

  if(!anc0) return;

  /* set node's i,j,k */
  k = kOfInd_n(ijk, ns);
  j = jOfInd_n_k(ijk, ns,k);
  i = iOfInd_n_jk(ijk, ns,j,k);
  /* set index of neighbor on same level */
  inb = i^1;   /* i index of hypothetical neighbor in same level */
  jnb = j^1;   /* j index of hypothetical neighbor in same level */
  knb = k^1;   /* k index of hypothetical neighbor in same level */

//
  count = 0;
  a0nb = anc0;
  fornodelist(anclist->next, anc1)
  {
    /* get descendant's i,j,k */
    int ijk1 = anc1->ijk;
    int k1 = kOfInd_n(ijk1, ns);
    int j1 = jOfInd_n_k(ijk1, ns,k1);
    int i1 = iOfInd_n_jk(ijk1, ns,j1,k1);

    // X dir only!!!
    a1nb = a0nb->child[Ind_n(i,j1,k1,  ns)];
    if(a1nb->leaf)
      break; // nb[fc] = a1nb;

    a0nb = a1nb;
  }
  nb[fc] = a1nb;


*******

  /* Y-dir: nb[2] and nb[3] */
  fs = 2 + jnb; /* face index in direction of sibling */
  fc = 2 + j;   /* face index in direction of cousin */
  /* set sibling neighbor */
  nb[fs] = parent->child[Ind_n(i,jnb,k, ns)];
  /* get cousin neighbor from parent's neighbor */
  if(grandp)
  {
      parentnb = grandp->child[Ind_n(pi,j,pk,  ns)];
    if(parentnb->leaf)
      nb[fc] = parentnb;
    else
      nb[fc] = parentnb->child[Ind_n(i,jnb,k,  ns)];
  }
  else
    nb[fc] = NULL;

*******


}







/* Return neighbors on the faces of a node in the node array nb.
   nb[f] is NULL if the neighbor does not exist in this node tree, i.e.
   in the patch the node is in. */
void neighbors_6faces(tNode *node, tNode *nb[6])
{
  int ijk = node->ijk;
  int ns[] = {2,2,2};
  int i,j,k;           /* node's i,j,k */
  int inb,jnb,knb;     /* neighbor i,j,k */
  int fs,fc;           /* sibling and cousin face index */
  tNode *parent = node->parent;
  int pijk;        /* parent's ijk */
  int pi,pj,pk;    /* parent's i,j,k */
  tNode *parentnb; /* parent's neighbor */
  tNode *grandp;   /* grandparent */
  int gijk;        /* grandparent's ijk */
  int gi,gj,gk;    /* grandparent's i,j,k */
  tNode *grandpnb; /* grandparent's neighbor */
  tNode *ggrandp = NULL;  /* great grandparent */
  tNode *gopnb;      /* grandparent or parent neighbor */

  /* init all neighbor pointers to NULL */
  for(i=0; i<6; i++) nb[i] = NULL;

  /* the root node has no neigbors in its patch */
  if(!parent) return;

  /* set parent's i,j,k */
  pijk = parent->ijk;
  pk = kOfInd_n(pijk, ns);
  pj = jOfInd_n_k(pijk, ns,pk);
  pi = iOfInd_n_jk(pijk, ns,pj,pk);

  /* grandparent info */
  grandp = parent->parent;
  if(grandp)
  {
    gijk = grandp->ijk;
    gk = kOfInd_n(gijk, ns);
    gj = jOfInd_n_k(gijk, ns,gk);
    gi = iOfInd_n_jk(gijk, ns,gj,gk);
    ggrandp = grandp->parent;
  }
  
  /* set node's i,j,k */
  k = kOfInd_n(ijk, ns);
  j = jOfInd_n_k(ijk, ns,k);
  i = iOfInd_n_jk(ijk, ns,j,k);

  /* set index of neighbor on same level */
  inb = i^1;   /* i index of hypothetical neighbor in same level */
  jnb = j^1;   /* j index of hypothetical neighbor in same level */
  knb = k^1;   /* k index of hypothetical neighbor in same level */




  /* Y-dir: nb[2] and nb[3] */
  dir = 1; 
  fs = 2*dir + jnb; /* face index in direction of sibling */
  fc = 2*dir + j;   /* face index in direction of cousin */
  /* set sibling neighbor */
  nb[fs] = parent->child[Ind_n(i,jnb,k, ns)];

  /* get other (e.g. cousin) neighbor from ancestor's neighbor */
  ancestor_with_diff_i_j_k(node, dir, &ancestor);

  /* now back down from ancestor's neighbor to same or closer level */

  if(ancestor)
  {
      parentnb = grandp->child[Ind_n(pi,j,pk,  ns)];
    if(parentnb->leaf)
      nb[fc] = parentnb;
    else
      nb[fc] = parentnb->child[Ind_n(i,jnb,k,  ns)];
  }
  else
    nb[fc] = NULL;



  if(grandp)
  {
      parentnb = grandp->child[Ind_n(pi,j,pk,  ns)];
    if(parentnb->leaf)
      nb[fc] = parentnb;
    else
      nb[fc] = parentnb->child[Ind_n(i,jnb,k,  ns)];
  }
  else
    nb[fc] = NULL;
















  /* X-dir: nb[0] and nb[1] */
  fs = inb; /* face index in direction of sibling */
  fc = i;   /* face index in direction of cousin */
  /* set sibling neighbor */
  nb[fs] = parent->child[Ind_n(inb,j,k, ns)];
  /* get cousin neighbor from parent's neighbor */
  if(grandp)
  {
    if(pi!=i) /* simple case: grandparent has child that is parentnb */
      gopnb = parentnb = grandp->child[Ind_n(i,pj,pk,  ns)];
    else /* the parent neighbor is not child of grandp */
    {
      if(ggrandp) /* look for great grandparent */
      {
        gopnb = grandpnb = ggrandp->child[Ind_n(inb,gj,gk,  ns)];
        if(grandpnb)
        {
          parentnb = grandpnb->child[Ind_n(inb,pj,pk,  ns)];
          if(parentnb) gopnb = parentnb;
        }
      }
      else /* give up and assume that there is no neighbor at all */
        gopnb = NULL;// => nb[fc] = NULL;
    }
    /* gopnb is now one of these: parentnb, grandpnb, NULL */
    if(gopnb)
    {
      if(gopnb->leaf)
        nb[fc] = gopnb;
      else
        nb[fc] = gopnb->child[Ind_n(inb,j,k,  ns)];
    }
    else
      nb[fc] = NULL;
  }
  else
    nb[fc] = NULL;
  
  /* Y-dir: nb[2] and nb[3] */
  fs = 2 + jnb; /* face index in direction of sibling */
  fc = 2 + j;   /* face index in direction of cousin */
  /* set sibling neighbor */
  nb[fs] = parent->child[Ind_n(i,jnb,k, ns)];
  /* get cousin neighbor from parent's neighbor */
  if(grandp)
  {
      parentnb = grandp->child[Ind_n(pi,j,pk,  ns)];
    if(parentnb->leaf)
      nb[fc] = parentnb;
    else
      nb[fc] = parentnb->child[Ind_n(i,jnb,k,  ns)];
  }
  else
    nb[fc] = NULL;
  
  /* Z-dir: nb[4] and nb[5] */
  fs = 4 + jnb; /* face index in direction of sibling */
  fc = 4 + j;   /* face index in direction of cousin */
  /* set sibling neighbor */
  nb[fs] = parent->child[Ind_n(i,j,knb, ns)];
  /* get cousin neighbor from parent's neighbor */
  if(grandp)
  {
      parentnb = grandp->child[Ind_n(pi,pj,k,  ns)];
    if(parentnb->leaf)
      nb[fc] = parentnb;
    else
      nb[fc] = parentnb->child[Ind_n(i,j,knb,  ns)];
  }
  else
    nb[fc] = NULL;
}


/* enter neighbor info as far as the 8 children of one parent are concerned */
/* this operates on a node array indexed by ijk */
void connect8_siblings(tNode *narray[])
{
  int ijk;

  /* fill in neighbor info, as far as these 8 are concerned */
  for(ijk=0; ijk<7; ijk++)
  {
    tNode *node = narray[ijk];
    switch(node->ijk)
    {
    case 0: // i,j,k=0,0,0
        node->nb[1][0] = narray[1]; // neig. in +X has i,j,k=1,0,0: ijk=1
        node->nb[3][0] = narray[2]; 
        node->nb[5][0] = narray[4];
        /* the 3 above are the only neigbors in these directions */
        node->nb[1][1] = node->nb[3][1] = node->nb[5][1] = NULL;
        break;
    case 1: // i,j,k=1,0,0
        node->nb[0][0] = narray[0]; // neig. in -X has i,j,k=0,0,0: ijk=0
        node->nb[3][0] = narray[3];
        node->nb[5][0] = narray[5];
        /* the 3 above are the only neigbors in these directions */
        node->nb[0][1] = node->nb[3][1] = node->nb[5][1] = NULL;
        break;
    case 2: // i,j,k=0,1,0
        node->nb[1][0] = narray[3];
        node->nb[2][0] = narray[0];
        node->nb[5][0] = narray[6];
        /* the 3 above are the only neigbors in these directions */
        node->nb[1][1] = node->nb[2][1] = node->nb[5][1] = NULL;
        break;
    case 3: // i,j,k=1,1,0
        node->nb[0][0] = narray[2]; // neig. in -X has i,j,k=0,1,0: ijk=2
        node->nb[2][0] = narray[1]; // neig. in -Y has i,j,k=1,0,0: ijk=1
        node->nb[5][0] = narray[7]; // neig. in +Z has i,j,k=1,1,1: ijk=7
        /* the 3 above are the only neigbors in these directions */
        node->nb[0][1] = node->nb[2][1] = node->nb[5][1] = NULL;
        break;
    case 4: // i,j,k=0,0,1
        node->nb[1][0] = narray[5];
        node->nb[3][0] = narray[6];
        node->nb[4][0] = narray[0];
        /* the 3 above are the only neigbors in these directions */
        node->nb[1][1] = node->nb[3][1] = node->nb[4][1] = NULL;
        break;
    case 5: // i,j,k=1,0,1
        node->nb[0][0] = narray[4];
        node->nb[3][0] = narray[7];
        node->nb[4][0] = narray[1];
        /* the 3 above are the only neigbors in these directions */
        node->nb[0][1] = node->nb[3][1] = node->nb[4][1] = NULL;
        break;
    case 6:
        node->nb[1][0] = narray[7];
        node->nb[2][0] = narray[4];
        node->nb[4][0] = narray[2];
        /* the 3 above are the only neigbors in these directions */
        node->nb[1][1] = node->nb[2][1] = node->nb[4][1] = NULL;
        break;
    case 7:
        node->nb[0][0] = narray[6];
        node->nb[2][0] = narray[5];
        node->nb[4][0] = narray[3];
        /* the 3 above are the only neigbors in these directions */
        node->nb[0][1] = node->nb[2][1] = node->nb[4][1] = NULL;
        break;
    }
  }
}

/* set neighbor connection info of 8 siblings with nodes from different
   parents. We assume that all priorly created nodes already have
   complete neighbor info. */
void connect8_____(tNode *narray[])
{
  tNode *parent = narray[0]->parent;
  tNode *oparent = narray[0]->parent;
  int ns[] = {2,2,2};

  /* level 0 and 1 nodes have only their siblings as neighbors */
  if(narray[0]->l <= 1) return;

  /* -X dir => i=0 */
/*
  i=0;
  // i,j,k = 0,0,0
  // loop over j,k
  ijk = Ind_n(i,j,k, ns)
  node = narray[ijk];

  parentnb = parent->nb[0]; // neig. in -X
  node->nb[0] = parentnb;
  if(parentnb)
  {
    ijk2 = Ind_n(i^1,j,k, ns);
    nbchild = parentnb->child[ijk2]
    if(nbchild)
    {
      node->nb[0] = nbchild;
      nbchild->nb[1] = node;
    }
  }
*/
}
