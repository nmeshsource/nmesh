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
