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

  /* fill MPI send and recv buffers */
  fornodelist(mesh->lns, elem)
  {
    node = elem->node;
    nid = node->nid;
    desrank = (nid/nperproc);
    if(node->datrank != desrank)
      move_node_to_rank(node, desrank, scom, rcom, 1);
  }
  nMPI_Waitall_com_send(scom);
  free_com(scom);  /* free scom with all its buffers */
  nMPI_Waitall_com_recv(rcom);

  /* get var data out of recv buffer */
  set_com_counters(rcom, 0,0);
  fornodelist(mesh->lns, elem)
  {
    node = elem->node;
    nid = node->nid;
    desrank = (nid/nperproc);
    if(node->datrank != desrank)
      move_node_to_rank(node, desrank, scom, rcom, 0);
  }

  free_com(rcom);
  update_mesh_myln_node_nid(mesh);
}

/* return: number of variables and number of doubles inside dat */
int nvars_ndoubles_in_dat(tDat *dat, int *ndoubles)
{
  tMesh *mesh;
  int nvars, vi;

  if(!dat)
  {
    *ndoubles = 0;
    return 0;
  }
  mesh = dat->node->pat->mesh;

  /* find amount of data */
  for(nvars=0, *ndoubles=0, vi=0; vi<dat->nv; vi++)
    if(dat->v[vi])
    {
      *ndoubles += (dat->v[vi]->N) * (MeshVarType(mesh, vi)!=1);
      if(PR) { PRF;printf(": vi=%d ndoubles=%d\n", vi, *ndoubles); }
      nvars++;
    }
  return nvars;
}

/* pack all (non-aux) vars of dat into a buffer which contains:
  |nvars||varind1|npoints1|<--data1-->||varind2|npoints2|<--data2-->||...
  here nvars is the number of variables that we send,
  the buffer has to be freed by caller later */
double *buffer_with_all_needed_dat_vars(tDat *dat, int *buflen)
{
  tMesh *mesh = dat->node->pat->mesh;
  int len;
  double *buf;
  int vi, datlen, nvars, bi, N;

  /* find amount of data */
  nvars = nvars_ndoubles_in_dat(dat, &datlen);
  if(PR) { PRF;printf(": nvars=%d datlen=%d\n", nvars, datlen); }

  /* alloc buffer */
  len = 1 + nvars*2 + datlen;
  buf = calloc(len, sizeof(double));

  /* fill buffer */
  buf[0] = nvars;
  for(bi=1, vi=0; vi<dat->nv; vi++)
  {
    /* add to buffer if eneabled and not auxiliary var */
    if(dat->v[vi])
    {
      N = (dat->v[vi]->N) * (MeshVarType(mesh, vi)!=1);
      buf[bi++] = vi;
      buf[bi++] = N;
      memcpy(buf+bi, dat->v[vi]->d, N * sizeof(double));
      bi += N;
      if(PR) { PRF;printf(": vi=%d bi=%d\n", vi, bi); }
    }
  }
  /* return pointer to this buffer, and its length */
  *buflen = bi;
  return buf;
}
/* put buffer back into dat and enable all vars needed */
int write_buffer_into_dat_vars(tDat *dat, double *buf)
{
  tNode *node = dat->node;
  int i, vi, nvars, bi, N;

  /* write buffer into vars */
  nvars = buf[0];
  for(bi=1, i=0; i<nvars; i++)
  {
    vi = buf[bi++];
    N  = buf[bi++];
    enablevarcomp_innode(node, vi);
    memcpy(dat->v[vi]->d, buf+bi, N * sizeof(double));
    bi += N;
    if(PR) { PRF;printf(": vi=%d bi=%d\n", vi, bi); }
  }
  return bi;
}

/* move data on one node between 2 ranks: the buffer we send/recv contains
  |nvars||varind1|npoints1|<--data1-->||varind2|npoints2|<--data2-->||...
 */
void move_node_to_rank(tNode *node, int desrank,
                       tCom *scom, tCom *rcom, int setbufs)
{
  tDat *dat = node->dat;
  int slen, rlen;
  double *sbuf, *rbuf;
  int rank = nMPI_rank();
  int other, rq;

  if(setbufs) /* setup buffers and fill them */
  {
    if(PR) { PRF;printf(": nid%ld datrank%d rank%d desrank%d\n",
                        node->nid, node->datrank, rank, desrank); }
    if(rank == node->datrank)
    {
      /* alloc and fill buffer */
      sbuf = buffer_with_all_needed_dat_vars(dat, &slen);
      if(PR)
      {
        PRF;printf(": sbuf =");
        for(int i=0; i<min2(3, slen); i++) printf(" %g", sbuf[i]);
        printf("\n");
        if(node->dat) printf("nv=%d\n", node->dat->nv);
      }
      other = desrank;
      /* put buffers in com */
      rq = append_buffers_to_com(scom, sbuf,slen, NULL,0);
      //print_com(scom);

      /* send */
      nMPI_Isend_double_com(scom, rq, other, node->nid, WORLD);
    }
    if(rank == desrank)
    {
      tMesh *mesh = node->pat->mesh;

      /* alloc buffer */
      rlen = 1 + (mesh->nvdb) * (2 + node->np); /* size to hold all vars */
      rbuf = calloc(rlen, sizeof(double));
      other = node->datrank;
      /* put buffers in com */
      rq = append_buffers_to_com(rcom, NULL,0, rbuf,rlen);
      //print_com(rcom);
      /* receive */
      nMPI_Irecv_double_com(rcom, rq, other, node->nid, WORLD);

      /* allocate space already and init some stuff */
      if(node->dat) errorexit("destination node should not have dat yet");
      node->dat = alloc_dat(node);

      if(PR) { PRF;printf(": calling coordinates_init_node\n"); }
      coordinates_init_node(node);
    }
  }
  else /* retrieve data from buffers */
  {
    double *rbuf;
    if(rank == desrank)
    {
      /* now unpack the buffers */
      rbuf = get_com_recv_i_buf(rcom);
      write_buffer_into_dat_vars(node->dat, rbuf);
      /* we can free the buffer here already */
      free_com_recv_i_buf(rcom);
      inc_com_recv_i(rcom);
    }
    else
    {
      free_dat(node->dat);
      node->dat = NULL;
    }
    node->datrank = desrank;
  }
}

/* move all nodes in list to rank */
void move_nodelist_to_rank(tNlist *list, int desrank)
{
  tNlist *elem, *list0;
  tNode *node = list->node;
  tMesh *mesh = node->pat->mesh;
  tCom *scom = alloc_com(sizeof(double), 1);
  tCom *rcom = alloc_com(sizeof(double), 1);

  if(PR) { PRF;printf(": desrank=%d\n", desrank); }

  /* find element 0 in list */
  list0 = first_nodelist(list);

  /* fill MPI send and recv buffers */
  fornodelist(list0, elem)
  {
    node = elem->node;
    if(node->datrank != desrank)
      move_node_to_rank(node, desrank, scom, rcom, 1);
  }
  nMPI_Waitall_com_send(scom);
  free_com(scom);  /* free scom with all its buffers */
  nMPI_Waitall_com_recv(rcom);

  /* get var data out of recv buffer */
  set_com_counters(rcom, 0,0);
  fornodelist(list0, elem)
  {
    node = elem->node;
    if(node->datrank != desrank)
      move_node_to_rank(node, desrank, scom, rcom, 0);
  }

  free_com(rcom);
  update_mesh_myln_node_nid(mesh);
}
