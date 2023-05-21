/* SurfExchange.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "SurfExchange.h"

#define PR 1



/* exchange some surfaces for testing */
int SurfExchange_test(tMesh *mesh)
{
  tNode *nd;
  int ui = Ind("SurfExchange_u");
  int l;

  PRF;printf(": Hmmm.\n");
  enablevar(mesh, ui);

  formylnodes(mesh)
  {
    int ijk;
    tNode *node = MyLnode;
    tArray *ua = VarA(node, ui);

    /* set particular pattern in u */
    forarray(ua, ijk)
    {
      int k = kOfInd_n(ijk, ua->n);
      int j = jOfInd_n_k(ijk, ua->n, k);
      int i = iOfInd_n_jk(ijk, ua->n, j,k);
      ua->d[ijk] = Node_eid(node) + 0.0007;
      if(i==0)          ua->d[ijk] += 0.1; 
      if(i==ua->n[0]-1) ua->d[ijk] += 0.2;
      if(j==0)          ua->d[ijk] += 0.03; 
      if(j==ua->n[1]-1) ua->d[ijk] += 0.04;
      if(k==0)          ua->d[ijk] += 0.005; 
      if(k==ua->n[2]-1) ua->d[ijk] += 0.006;
    }
  }

  /* print var */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    printelm(node);
    printvar_innode(node, ui);
  }

  /* print var in one node again */
  nd = Lnode_myid(mesh, 0); /* my first node */
  printelm(nd);
  printvar_innode(nd, ui);

  /* exchange surfaces */
  prdivider('S');
  PRF;printf(": exchange surfaces\n");
  init_all_myln_surfaces(mesh);
  /* instead of:
       set_all_myln_mysurf(mesh);
       request_all_myln_surfaces_exchange(mesh);
       get_all_myln_surfaces(mesh);
     we do the exchange nodewise */

  /* do exchange 4 times similar to RK4 */
  for(l=0; l<4; l++)
  {
    formylnodes(mesh)
    {
      tNode *node = MyLnode;

      set_all_mysurf(node);
      request_all_surfaces_exchange(node);
    }

    /* Here we can do work. MPI is now busy sending buffers */

    /* now get the surfaces and wait for buffers if necessary */
    formylnodes(mesh)
    {
      tNode *node = MyLnode;

      get_all_surfaces(node);
    }
  }

  /* get_all_myln_surfaces sets ajsurf via interpolation */

  /* print var in one node yet again with surfaces */
  nd = Lnode_myid(mesh, 0); /* my first node */
  printelm(nd);
  printvar_innode(nd, ui);

  /* after we have printed them, we no longer need the surfaces */
  free_all_myln_surfaces(mesh);
  return 0;
}
