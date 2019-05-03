/* refine.c */
/* Wolfgang Tichy, 5/2019 */


#include "nmesh.h"
#include "amr.h"

#define PR 0



/* refine nodes with nids in array, we assume long *nid is sorted in ascending
   order. We do not update nid's in here */
void refine_nodes_without_nid_update__old(tMesh *mesh, long nnodes, long *nid)
{
  tNlist *elem = mesh->lns;
  long i;

  for(i=0; i<nnodes; i++)
  {
    tNode *parent;
    tNlist *children;
    tNlist *lastchild;
    int *n;

    /* forward to node with nid[i] */
    for(; elem->node->nid != nid[i]; elem = elem->next) ;

    /* make children */
    parent = elem->node;
    n = elem->node->n; /* pick n */
    children = make8_child_nodes(parent, n);

    /* update mesh->lns if needed and add children to list */
    if(elem == mesh->lns) mesh->lns = first_nodelist(children);
    lastchild = replace1_in_nodelist(elem, children, 1);

    /* set elem to last child we added */
    elem = lastchild;
  }
}


/* refine nodes with nids in array, we assume long *nid is sorted in ascending
   order. We do not update nid's in here */
void refine_nodes_without_nid_update(tMesh *mesh, long nnodes, long *nid)
{
  tNlist **replace  = calloc(nnodes, sizeof(replace[0]));
  tNlist **children = calloc(nnodes, sizeof(children[0]));
  long i;

  if(!replace || !children) errorexit("no memory for replace, children");

  FORNODES_Pragma(omp parallel)
  {
    tNlist *elem = mesh->lns;

    FORNODES_Pragma(omp for)
    for(i=0; i<nnodes; i++)
    {
      tNode *parent;
      int *n;

      /* forward to node with nid[i] */
      for(; elem->node->nid != nid[i]; elem = elem->next) ;

      /* make children */
      parent = elem->node;
      n = elem->node->n; /* pick n */
      children[i] = make8_child_nodes(parent, n);

      /* save elem that has to be replaced by children[i] later */
      replace[i] = elem;
    }
  }

  /* update mesh->lns using info in replace[i], children[i] */
  for(i=0; i<nnodes; i++)
  {
    tNlist *elem = replace[i];

    /* update mesh->lns if needed and add children to list */
    if(elem == mesh->lns) mesh->lns = first_nodelist(children[i]);
    replace1_in_nodelist(elem, children[i], 1); // can return last child
  }
  free(children);
  free(replace);
}
