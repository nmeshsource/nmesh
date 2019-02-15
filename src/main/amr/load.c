/* load.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 1

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
  int nperproc = nnodes/size +  nnodes%size;
  int desrank;
  tCom *scom = alloc_com(sizeof(double), 1);
  tCom *rcom = alloc_com(sizeof(double), 1);

  PRF;printf(": nperproc=%d\n", nperproc);

  fornodelist(mesh->lns, elem)
  {
    node = elem->node;
    nid = node->nid;
    desrank = (nid/nperproc);
    if(node->datrank != desrank)
    {
      move_node_to_rank(node, desrank, scom, rcom, 1);
    }
  }
  nMPI_Waitall_com_send(scom);
  nMPI_Waitall_com_recv(rcom);

  set_com_counters(rcom, 0);
  fornodelist(mesh->lns, elem)
  {
    node = elem->node;
    nid = node->nid;
    desrank = (nid/nperproc);
    if(node->datrank != desrank)
    {
      move_node_to_rank(node, desrank, scom, rcom, 0);
    }
  }

  /* free com with all their buffers */
  free_com(scom);
  free_com(rcom);
  update_mesh_myln_node_nid(mesh);
}

void move_node_to_rank(tNode *node, int desrank,
                       tCom *scom, tCom *rcom, int setbufs)
{
  int slen=1, rlen=1;
  double *sbuf, *rbuf;
  int rank = nMPI_rank();
  int other, rq;

  if(setbufs) /* setup buffers and fill them */
  {
    if(PR) printf("nid%ld datrank%d rank%d desrank%d\n",
                  node->nid, node->datrank, rank, desrank);
    if(rank == node->datrank)
    {
      /* alloc buffer */
      sbuf = calloc(slen, sizeof(double));
      /* fill send buffer to be sent to desrank */
      //..
      //test:
      sbuf[0] = node->nid + 0.1234;
      other = desrank;
      /* put buffers in com */
      rq = append_buffers_to_com(scom, sbuf,slen, NULL,0);
      print_com(scom);

      /* send */
      nMPI_Isend_double_com(scom, rq, other, node->nid);
    }
    if(rank == desrank)
    {
      /* alloc buffer */
      rbuf = calloc(rlen, sizeof(double));
      other = node->datrank;
      /* put buffers in com */
      rq = append_buffers_to_com(rcom, NULL,0, rbuf,rlen);
      print_com(rcom);
      /* receive */
      nMPI_Irecv_double_com(rcom, rq, other, node->nid);
    }
  }
  else /* retrieve data from buffers */
  {
    double *rbuf;
    if(rank == desrank)
    {
      /* now unpack the buffers */
      //...
      //test:
      rbuf = get_next_com_recv_buf(rcom);
      if(node->dat) errorexit("destination node should not have dat yet");
      node->dat = alloc_dat(node);
      PRF;printf(": nid%ld rank%d recv %g\n", node->nid, rank, rbuf[0]);
    }
    else
    {
      free_dat(node->dat);
      node->dat = NULL;
    }
    node->datrank = desrank;
  }
}
