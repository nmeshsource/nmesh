/* evolve.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "evolve.h"

#define PR 0

/* Determine and set trouble score in a node,
   i.e. set node->dat->info->trouble .
   Note evosys->f[TROUBLE] returns a score of 1, 0, or -1:
    1 trouble  =>  switch to more robust method (e.g. fv)
    0 node ok with current method  =>  do nothing
   -1 all ok   =>  switch to more accurate but delicate method (e.g. dg)
   examples: if shock in dg return 1
             if shock in fv return 0
             if no shock in dg return 0
             if no shock in fv return -1
   Here we also accumulate the 1 and -1 scores and return them */
int evolve_set_trouble_score(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  tEvoSys *evosys = mesh->evosys;
  pVLList *u_p = evosys->u_p;
  pVLList *u   = evosys->u;
  pVLList *r   = evosys->rhs;
  int max_trouble = 1073741824;
  int troubled = -1; /* default is to assume all is very well */
  int i;

  if(node->dat == NULL) errorexit("node->dat is NULL");

  if(PR) PRFs(":\n");

  /* check all evo systems for trouble and accumulate result in troubled */
  forList(u, i)
  {
    tVarList *vlu_p = ListEntry(u_p,i);
    tVarList *vlu   = ListEntry(u,i);
    tVarList *vlr   = ListEntry(r,i);
    int tr = 0;

    /* run TROUBLE func */
    if(troubled<=0) /* need to check only if there no trouble yet */
      if(ListEntry(evosys->f[TROUBLE],i))
        tr = ListEntry(evosys->f[TROUBLE],i)(node, vlr, vlu, vlu_p);

    /* set troubled flag to 1, 0, or -1 (default above) */
    if(tr>0) /* new trouble found, switch to fv */
    {
      troubled = 1;
    }
    else if(tr==0) /* ok, keep as is */
    {
      if(troubled<0) troubled = 0;
    }
  }

  /* set node->dat->info->trouble */
  if(troubled>0) /* i.e. there is trouble now */
  {
    /* if there was no trouble, we set the trouble score to 1 */
    if(node->dat->info->trouble<=0) node->dat->info->trouble  = 1;
    /* if there was trouble before, we continuously increase the score */
    else                            node->dat->info->trouble += 1;
  }
  else if(troubled==0) /* is ok, keep node as is */
  {
    node->dat->info->trouble = 0;
  }
  else /* all is good, can switch back to dg */
  {
    /* if there was trouble, we set the trouble score to -1 */
    if(node->dat->info->trouble>=0) node->dat->info->trouble  = -1;
    /* if there was no trouble, we continuously lower the score */
    else                            node->dat->info->trouble -= 1;
  }

  /* make sure trouble score does not become too large or small */
  if(node->dat->info->trouble < -max_trouble)
    node->dat->info->trouble = -max_trouble;
  if(node->dat->info->trouble > max_trouble)
    node->dat->info->trouble = max_trouble;

  return node->dat->info->trouble;
}

/* determine and set trouble score in each node,
   i.e. set node->dat->info->trouble */
int evolve_set_trouble_score_mesh(tMesh *mesh)
{
  int Max_trb, max_trb=INT_MIN;
  int Min_trb, min_trb=INT_MAX;
  if(PR) PRFs(":\n");

  /* loop over all nodes, check for trouble, and node-info trouble score */
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    evolve_set_trouble_score(node);
  }

  /* now find rank-local max of node->dat->info->trouble */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    int trb = node->dat->info->trouble;
    if(trb>max_trb) max_trb = trb; /* max trouble */
    if(trb<min_trb) min_trb = trb; /* min trouble */
  }

  /* Max over all ranks */
  Max_trb = max_trb;
  nMPI_Allreduce(&max_trb, &Max_trb, 1, nMPI_INT, nMPI_MAX);

  /* Min over all ranks */
  Min_trb = min_trb;
  nMPI_Allreduce(&min_trb, &Min_trb, 1, nMPI_INT, nMPI_MIN);

  /* if there is bad trouble somewhere */
  if(Max_trb>0)
    return Max_trb; /* returns max. trouble value of all nodes */
  else
    return Min_trb; /* no trouble, return smallest to signal need for dg */
}
