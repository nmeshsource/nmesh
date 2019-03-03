/* coordinates.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "coordinates.h"

#define PR 1
#define Vd(node,i) GetVarDpointer(node,i)



/* try to enable all coord vars, return 1 if we have dat */
int coordinates_coordvars_enabled(tNode *node)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  int iX   = Ind("X");
  int idXd = Ind("dXdx");
  int ix   = Ind("x");

  if(!node->dat) return 0;

  /* give all these memory */
  enablevar_innode(node, iX);
  enablevar_innode(node, iX+1);
  enablevar_innode(node, iX+2);
  enablevar_innode(node, idXd);
  enablevar_innode(node, idXd+3);
  enablevar_innode(node, idXd+6);
  enablevar_innode(node, ix);
  enablevar_innode(node, ix+1);
  enablevar_innode(node, ix+2);
  return 1;
}

/* (re)initialize coordinates in a node */
int coordinates_init_node(tNode *node)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  int *n = node->n;
  int i,j,k, d,e;
  int vars_on = coordinates_coordvars_enabled(node);
  int iX   = Ind("X");
  int idXd = Ind("dXdx");
  int ix   = Ind("x");
  double *pX[] = { Vd(node,iX), Vd(node,iX+1), Vd(node,iX+2) };
  double *px[] = { Vd(node,ix), Vd(node,ix+1), Vd(node,ix+2) };
  double *pdXd[3][3] = { {Vd(node,idXd),   Vd(node,idXd+1), Vd(node,idXd+2)},
                         {Vd(node,idXd+3), Vd(node,idXd+4), Vd(node,idXd+5)},
                         {Vd(node,idXd+6), Vd(node,idXd+7), Vd(node,idXd+8)} };
  PRF;printf(":\n");
  if(!vars_on) return 0;

  /* set coords */
  forijk(i,j,k, n)
  {
    double Xb[] = { node->Xb[0]->d[i], node->Xb[1]->d[j], node->Xb[2]->d[k] };
    double X[3], x[3], dXd[3][3];
    int ijk = Ind_n(i,j,k, n);

    /* get X from Xb */
    XYZ_of_XbYbZb(node, Xb, X);
    for(d=0; d<3; d++) px[d][ijk] = pX[d][ijk] = X[d];

    /* now set x, dXdx */
    if(pat->dXYZ_dxyz)
    {
      pat->dXYZ_dxyz(node, -1, X, x, dXd);
      for(d=0; d<3; d++)
      {
        px[d][ijk] = x[d];
        for(e=0; e<3; e++) pdXd[d][e][ijk] = dXd[d][e];
      }
    }
  }

  return 0;
}

/* initialize coordinates in each patch */
int coordinates_init(tMesh *mesh)
{
  int myid;

  PRF;printf(":\n");

  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    coordinates_init_node(node);
  }

  return 0;
}
