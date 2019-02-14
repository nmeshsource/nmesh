/* load.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "amr.h"


/* functions to move nodes between procs */

/**********************************************************************/
/*  */
/**********************************************************************/

void simple_load_balance(tMesh *mesh)
{
  long nnodes = mesh->nln;
  long nid;
  tNlist *elem;
  tNode *node;
  int size = nMPI_size();
  double nperproc = nnodes/size;
  int desrank;
  tCom *com = alloc_com(sizeof(double), 1);

  fornodelist(mesh->lns, elem)
  {
    node = elem->node;
    nid = node->nid;
    desrank = (nid/nperproc);
    if(node->datrank != desrank)
      move_node_to_rank(com, node, desrank);
  }
  nMPI_Waitall_in_com(com);
  /* now unpack the buffers */
  //...
  /* free com with all its buffers */
  free_com(com);
}

void move_node_to_rank(tCom *com, tNode *node, int desrank)
{
  int rn;
  int slen=1, rlen=1;
  
  /* alloc buffers */
  double *sbuf = calloc(slen, sizeof(double));
  double *rbuf = calloc(rlen, sizeof(double));

  /* fill send buffer to be sent to desrank */
  //..

  /* put buffers in com */
  rn = com->n_rq;
  realloc_com_reqs(com, rn + 1);
  put_buffers_in_com(com, rn, sbuf,slen, rbuf,rlen);

  /* now call MPI */
  nMPI_Isend_Irecv_double_com(com, com->n_rq - 1,
                              desrank, node->nid, node->nid);
}
