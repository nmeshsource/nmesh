/* SurfExchange.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "SurfExchange.h"

#define PR 1



/* print some compile info */
int SurfExchange_test(tMesh *mesh)
{
  int ui = Ind("SurfExchange_u");
  int li, ijk;

  PRF;printf(": Hmmm.\n");
  enablevar(mesh, ui);

  formylnodes(mesh, li)
  {
    tNode *node = GetMyNode(mesh, li);
    tArray *ua = GetVarArray(node, ui);

    /* set particular patter in u */
    forarray(ua, ijk)
    {
      int k = kOfInd_n(ijk, ua->n);
      int j = jOfInd_n_k(ijk, ua->n, k);
      int i = iOfInd_n_jk(ijk, ua->n, j,k);
      ua->d[ijk] = node->nid + 0.0007;
      if(i==0)          ua->d[ijk] += 0.1; 
      if(i==ua->n[0]-1) ua->d[ijk] += 0.2;
      if(j==0)          ua->d[ijk] += 0.03; 
      if(j==ua->n[1]-1) ua->d[ijk] += 0.04;
      if(k==0)          ua->d[ijk] += 0.005; 
      if(k==ua->n[2]-1) ua->d[ijk] += 0.006;
    }
  }

  /* print var */
  formylnodes(mesh, li)
  {
    tNode *node = GetMyNode(mesh, li);
    printnode(node);
    printvar_innode(node, ui);
  }

  /* now exchange surfaces */
  // ...

  /* print var again */
  formylnodes(mesh, li)
  {
    tNode *node = GetMyNode(mesh, li);
    printnode(node);
    printvar_innode(node, ui);
  }

  return 0;
}
