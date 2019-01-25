/* connect.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"




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
