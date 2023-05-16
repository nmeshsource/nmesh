/* GhostExchange.c */
/* Wolfgang Tichy, 12/2019 */

#include "nmesh.h"
#include "GhostExchange.h"

#define PR 1



/* exchange some ghosts for testing */
int GhostExchange_test(tMesh *mesh)
{
  //tNode *nd;
  int ui = Ind("GhostExchange_u");
  int vi = Ind("GhostExchange_v");
  int l;

  PRF;printf(": Hmmm.\n");
  enablevar(mesh, ui);
  enablevar(mesh, vi);

  formylnodes(mesh)
  {
    int ijk;
    tNode *node = MyLnode;
    tArray *va = VarA(node, vi);

    /* set particular pattern in v */
    forarray(va, ijk)
    {
      int k = kOfInd_n(ijk, va->n);
      int j = jOfInd_n_k(ijk, va->n, k);
      int i = iOfInd_n_jk(ijk, va->n, j,k);
      va->d[ijk] = Node_eid(node) + 0.0001*(ijk%9+1);
      if(i==0)          va->d[ijk] += 0.1;
      if(i==va->n[0]-1) va->d[ijk] += 0.2;
      if(j==0)          va->d[ijk] += 0.03;
      if(j==va->n[1]-1) va->d[ijk] += 0.04;
      if(k==0)          va->d[ijk] += 0.005;
      if(k==va->n[2]-1) va->d[ijk] += 0.006;
    }
  }

  /* print var */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printnode(node);
    printvar_innode(node, vi);
  }

//  /* print var in one node again */
//  nd = Lnode_myid(mesh, 0); /* my first node */
//  printnode(nd);
//  printvar_innode(nd, vi);

  /* exchange ghosts */
  prdivider('G');
  PRF;printf(": exchange ghosts\n");
  MPIexchange_init_all_myln(mesh);
  MPIexchange_set_all_myln_localdata(mesh);

  /* do exchange 4 times similar to RK4 */
  for(l=0; l<4; l++)
  {
    MPIexchange_request_all_myln_data(mesh);

    /* here useful work can be done */

    MPIexchange_get_all_myln_data(mesh);
  }

  /* print var again */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printnode(node);
    printvar_innode(node, vi);
  }

//  /* print var in one node yet again */
//  nd = Lnode_myid(mesh, 0); /* my first node */
//  printnode(nd);
//  printvar_innode(nd, vi);

  /* after we have printed them, we no longer need the ghosts */
  MPIexchange_free_all_myln(mesh);
  return 0;
}
