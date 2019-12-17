/* main/amr/indicators.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 0


/* functions to exchange ghost zone data on a uniform mesh */


/***************************************************************************/
/* neighbor indices */
/***************************************************************************/

/* indices of neighbors, next-nearest and next to next-nearest neighbors */
enum
{
  mcc=0, pcc, cmc, cpc, ccm, ccp, // nearest neighbors on faces 0-6
  mmc, mpc, pmc, ppc,             // next-nearest nbs on face 0,1 in Y-dir
  mcm, pcm, cmm, cpm,             // next-nearest nbs on face 4
  mcp, pcp, cmp, cpp,             // next-nearest nbs on face 5
  mmm, pmm, mpm, ppm,             // next to next-nearest nbs on face 4
  mmp, pmp, mpp, ppp              // next to next-nearest nbs on face 5
};

/***************************************************************************/
/* find neighbors */
/***************************************************************************/

/* find nearest + next-nearest + next-next-nearest neighbors */
void find_nb26(tNode *node, tNode *nb26[26])
{
  /* nearest neighbors */
  nb26[mcc] = node->fnb[mcc][0];
  nb26[pcc] = node->fnb[pcc][0];
  nb26[cmc] = node->fnb[cmc][0];
  nb26[cpc] = node->fnb[cpc][0];
  nb26[ccm] = node->fnb[ccm][0];
  nb26[ccp] = node->fnb[ccp][0];

  /* next-nearest nbs on face 0,1 */
  if(nb26[mcc])
  {
    nb26[mmc] = nb26[mcc]->fnb[cmc][0];
    nb26[mpc] = nb26[mcc]->fnb[cpc][0];
  }
  else
  {
    nb26[mmc] = nb26[mpc] = NULL;
  }
  if(nb26[pcc])
  {
    nb26[pmc] = nb26[pcc]->fnb[cmc][0];
    nb26[ppc] = nb26[pcc]->fnb[cpc][0];
  }
  else
  {
    nb26[pmc] = nb26[ppc] = NULL;
  }

  /* next-nearest nbs on face 4 */
  if(nb26[ccm])
  {
    nb26[mcm] = nb26[ccm]->fnb[mcc][0];
    nb26[pcm] = nb26[ccm]->fnb[pcc][0];
    nb26[cmm] = nb26[ccm]->fnb[cmc][0];
    nb26[cpm] = nb26[ccm]->fnb[cpc][0];
  }
  else
  {
    nb26[mcm] = nb26[pcm] = nb26[cmm] = nb26[cpm] = NULL;
  }

  /* next-nearest nbs on face 5 */
  if(nb26[ccp])
  {
    nb26[mcp] = nb26[ccp]->fnb[mcc][0];
    nb26[pcp] = nb26[ccp]->fnb[pcc][0];
    nb26[cmp] = nb26[ccp]->fnb[cmc][0];
    nb26[cpp] = nb26[ccp]->fnb[cpc][0];
  }
  else
  {
   nb26[mcp] = nb26[pcp] = nb26[cmp] = nb26[cpp] = NULL;
  }

  /* next to next-nearest nbs on face 4 */
  if(nb26[cmm])
  {
    nb26[mmm] = nb26[cmm]->fnb[mcc][0];
    nb26[pmm] = nb26[cmm]->fnb[pcc][0];
  }
  else
  {
    nb26[mmm] = nb26[pmm] = NULL;
  }
  if(nb26[cpm])
  {
    nb26[mpm] = nb26[cpm]->fnb[mcc][0];
    nb26[ppm] = nb26[cpm]->fnb[pcc][0];
  }
  else
  {
    nb26[mpm] = nb26[ppm] = NULL;
  }

  /* next to next-nearest nbs on face 5 */
  if(nb26[cmp])
  {
    nb26[mmp] = nb26[cmp]->fnb[mcc][0];
    nb26[pmp] = nb26[cmp]->fnb[pcc][0];
  }
  else
  {
    nb26[mmp] = nb26[pmp] = NULL;
  }
  if(nb26[cpp])
  {
    nb26[mpp] = nb26[cpp]->fnb[mcc][0];
    nb26[ppp] = nb26[cpp]->fnb[pcc][0];
  }
  else
  {
    nb26[mpp] = nb26[ppp] = NULL;
  }
}

/* calculate neighbor index nb_ni of thes node as seen from neighbor ni */
int get_nb_ni(int ni)
{
  switch(ni)
  {
  case mcc: return pcc;
  case pcc: return mcc;
  case cmc: return cpc;
  case cpc: return cmc;
  case ccm: return ccp;
  case ccp: return ccm;
  case mmc: return ppc;
  case mpc: return pmc;
  case pmc: return mpc;
  case ppc: return mmc;
  case mcm: return pcp;
  case pcm: return mcp;
  case cmm: return cpp;
  case cpm: return cmp;
  case mcp: return pcm;
  case pcp: return mcm;
  case cmp: return cpm;
  case cpp: return cmm;
  case mmm: return ppp;
  case pmm: return mpp;
  case mpm: return pmp;
  case ppm: return mmp;
  case mmp: return ppm;
  case pmp: return mpm;
  case mpp: return pmm;
  case ppp: return mmm;
  default: errorexit("ni has unexpected value!");
  }
}

/* select start value sta, start value of ghosts gsta and thickness nc
  of layer to copy, and also start value nb_sta on nb */
void set_ghoststart_start_nc(int ni, int n[3], int nghosts, int gsta[3],
                             int sta[3], int nc[3], int nb_sta[3])
{
  int n2gho = nghosts*2;

  /* select start values for copy */
  gsta[0] = 0;
  gsta[1] = 0;
  gsta[2] = 0;
  sta[0] = nghosts;
  sta[1] = nghosts;
  sta[2] = nghosts;
  nb_sta[0] = n[0] - n2gho;
  nb_sta[1] = n[1] - n2gho;
  nb_sta[2] = n[2] - n2gho;
  switch(ni)
  {
  case pcc:
  case pmc:
  case pcm:
  case pmm:
    gsta[0] = n[0] - nghosts;
    sta[0] = n[0] - n2gho;
    nb_sta[0] = nghosts;
    break;

  case cpc:
  case mpc:
  case cpm:
  case mpm:
    gsta[0] = n[0] - nghosts;
    sta[1] = n[1] - n2gho;
    nb_sta[1] = nghosts;
    break;

  case ccp:
  case mcp:
  case cmp:
  case mmp:
    gsta[0] = n[0] - nghosts;
    sta[2] = n[2] - n2gho;
    nb_sta[2] = nghosts;
    break;

  case ppc:
  case ppm:
    gsta[0] = n[0] - nghosts;
    gsta[1] = n[1] - nghosts;
    sta[0] = n[0] - n2gho;
    sta[1] = n[1] - n2gho;
    nb_sta[0] = nb_sta[1] = nghosts;
    break;

  case pcp:
  case pmp:
    gsta[0] = n[0] - nghosts;
    gsta[2] = n[2] - nghosts;
    sta[0] = n[0] - n2gho;
    sta[2] = n[2] - n2gho;
    nb_sta[0] = nb_sta[2] = nghosts;
    break;

  case cpp:
  case mpp:
    gsta[1] = n[1] - nghosts;
    gsta[2] = n[2] - nghosts;
    sta[1] = n[1] - n2gho;
    sta[2] = n[2] - n2gho;
    nb_sta[1] = nb_sta[2] = nghosts;
    break;

  case ppp:
    gsta[0] = n[0] - nghosts;
    gsta[1] = n[1] - nghosts;
    gsta[2] = n[2] - nghosts;
    sta[0] = n[0] - n2gho;
    sta[1] = n[1] - n2gho;
    sta[2] = n[2] - n2gho;
    nb_sta[0] = nb_sta[1] = nb_sta[2] = nghosts;
    break;
  }

  /* select thinkness of layer to copy */
  nc[0] = n[0];
  nc[1] = n[1];
  nc[2] = n[2];
  switch(ni)
  {
  case mcc:
  case pcc:
    nc[0] = nghosts;
    break;

  case cmc:
  case cpc:
    nc[1] = nghosts;
    break;

  case ccm:
  case ccp:
    nc[2] = nghosts;
    break;

  case mmc:
  case mpc:
  case pmc:
  case ppc:
    nc[0] = nc[1] = nghosts;
    break;

  case mcm:
  case mcp:
  case pcm:
  case pcp:
    nc[0] = nc[2] = nghosts;
    break;

  case cmm:
  case cmp:
  case cpm:
  case cpp:
    nc[1] = nc[2] = nghosts;
    break;

  default:
    nc[0] = nc[1] = nc[2] = nghosts;
  }
}


/***************************************************************************/
/* put ghost data into buffers and start MPI send/recv to get data */
/***************************************************************************/

/* setup ghost exchange */
void request_ghostdata_for_vl(tNode *node, tVarList  *vl)
{
  int nghosts = 2;
  tNode *nb26[26];
  tDat *dat = node->dat;
  int *rqs;
  int vi0 = Vind(vl,0); /* 1st var */
  int ni;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  /* alloc memory to store request numbers for all 26 nbs */
  if(VarA_(node, vi0)->par) free(VarA_(node, vi0)->par);
  rqs = imalloc(26);
  VarA_(node, vi0)->par = (void *) rqs;

  /* find 26 neighbors */
  find_nb26(node, nb26);

  /* loop over 26 neighbors */
  for(ni=0; ni<26; ni++)
  {
    tNode *nb = nb26[ni];
    int *n = node->n;
    int *nb_n = nb->n;
    int vli;
    int gsta[3], sta[3], nc[3], nb_sta[3];

    /* goto next neighbor if nb is NULL */
    if(!nb) continue;

    /* select start value and thickness of layer to copy */
    set_ghoststart_start_nc(ni, n, nghosts, gsta, sta, nc, nb_sta);

    /* is nb local? */
    if(nb->dat)
    {
      /* nb is local so just copy */
      for(vli=0; vli<vl->n; vli++)
      {
        int vi = Vind(vl, vli);
        double *src = Vard_(nb, vi);
        double *des = Vard_(node, vi);
        int i0,j0,k0;

        /* copy from nb into ghosts on node */
        for(k0=0; k0<nc[2]; k0++)
        for(j0=0; j0<nc[1]; j0++)
        for(i0=0; i0<nc[0]; i0++)
        {
          int k = k0 + gsta[2];
          int j = j0 + gsta[1];
          int i = i0 + gsta[0];
          int ijk = Ind_n(i,j,k, n);
          int nb_k = k0 + nb_sta[2];
          int nb_j = j0 + nb_sta[1];
          int nb_i = i0 + nb_sta[0];
          int nb_ijk = Ind_n(nb_i,nb_j,nb_k, nb_n);
          des[ijk] = src[nb_ijk];
        }
      }
    }
    else
    {
      /* nb is on other process so use MPI to exchange data */
      int nvars = vl->n;
      int buflen = nvars * nc[0]*nc[1]*nc[2];
      int rq, nb_rank, s_tag, r_tag;
      nMPI_Comm s_comm, r_comm;
      tCom *com = dat->gcom;
      double *sbuf, *rbuf; /* buffers for MPI */
      int si, nb_ni, lid, nb_lid;

      /* find my index in nb */
      nb_ni = get_nb_ni(ni);

      /* set up tags */
      nb_rank = nb->datrank;
      lid = calc_node_lid(node);
      nb_lid = calc_node_lid(nb);
      s_tag = nb_lid*26 + nb_ni;
      r_tag = lid*26 + ni;
      if(r_tag<0) r_tag = -r_tag;
      if(s_tag<0) s_tag = -s_tag;
      r_comm = nb->comm;
      s_comm = node->comm;

      /* alloc buffers for send/recv */
      rbuf = dmalloc(buflen);
      sbuf = dmalloc(buflen);

      /* save buffers in com */
      rq = append_buffers_to_com(com, sbuf,buflen, rbuf,buflen);

      /* save MPI request rq in par of 1st var */
      rqs[ni] = rq;

      /* copy var data into send buffer */
      for(si=0, vli=0; vli<vl->n; vli++)
      {
        int vi = Vind(vl, vli);
        double *v = Vard_(node, vi);
        int i0,j0,k0;

        /* copy from var on node into sbuf */
        for(k0=0; k0<nc[2]; k0++)
        for(j0=0; j0<nc[1]; j0++)
        for(i0=0; i0<nc[0]; i0++)
        {
          int k = k0 + sta[2];
          int j = j0 + sta[1];
          int i = i0 + sta[0];
          int ijk = Ind_n(i,j,k, n);

          sbuf[si++] = v[ijk];
        }
      }

      /* use MPI to recv/send buffers */
      nMPI_Isend_Irecv_double_com(com, rq, nb_rank, s_tag,r_tag,
                                  s_comm,r_comm);
    }
  } /* end loop over neighbors */
}


/* request ghost exchanges on all my nodes in the mesh
   Note: We need to call this! If we call request_ghostdata_for_vl(n1)
   for only node n1, MPI deadlocks because the other nodes are not sending
   to n1 or receiving from n1 */
void request_all_myln_ghostdata_for_vl(tMesh *mesh, tVarList  *vl)
{

  TIMER_START;

  /* If we want threads in this loop, we need MPI_Init_thread with
     MPI_THREAD_MULTIPLE, instead of just MPI_Init in main. */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    request_ghostdata_for_vl(node, vl);
  }
  TIMER_STOP;
}


/**********************************************************************/
/* get the ghost data out of the MPI buffers */
/**********************************************************************/

/* get all ghost data from neighbors */
void get_ghostdata_for_vl(tNode *node, tVarList  *vl)
{
  int nghosts = 2;
  tNode *nb26[26];
  tDat *dat = node->dat;
  tCom *com = dat->gcom;
  int ni, nvars, rq;
  int vi0 = Vind(vl,0); /* 1st var */
  int *rqs;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  /* do nothing if com is empty */
  if(com->n_rq == 0) return;

  /* get MPI rq numbers */
  rqs = (int *) VarA_(node, vi0)->par;

  /* find 26 neighbors */
  find_nb26(node, nb26);

  /* loop over 26 neighbors */
  for(ni=0; ni<26; ni++)
  {
    tNode *nb = nb26[ni];
    int *n = node->n;
    double *rbuf; /* buffer for recv */
    int si, vli;
    int gsta[3], sta[3], nc[3], nb_sta[3];

    /* goto next neighbor if nb is NULL */
    if(!nb) continue;

    /* if nb is local we have already exchanged info  */
    if(nb->dat) continue;

    /* nb is on other process, we have used MPI to exchange data */

    /* number of vars that exchanged ghost */
    nvars = vl->n;

    /* do nothing if there are no vars that exchanged ghost */
    if(!nvars) return;

    /* select start value and thickness of layer to copy */
    set_ghoststart_start_nc(ni, n, nghosts, gsta, sta, nc, nb_sta);

    /* get MPI request number, and buffer */
    rq = rqs[ni];
    rbuf = get_com_recv_buf(com, rq);

    /* wait for MPI buffer */
    nMPI_Wait_com_recv(com, rq);

    /* copy data from recv buffer into ghost zones */
    for(si=0, vli=0; vli<vl->n; vli++)
    {
      int vi = Vind(vl, vli);
      double *v = Vard_(node, vi);
      int i0,j0,k0;

      /* copy from rbuf into var */
      for(k0=0; k0<nc[2]; k0++)
      for(j0=0; j0<nc[1]; j0++)
      for(i0=0; i0<nc[0]; i0++)
      {
        int k = k0 + gsta[2];
        int j = j0 + gsta[1];
        int i = i0 + gsta[0];
        int ijk = Ind_n(i,j,k, n);

        v[ijk] = rbuf[si++];
      }
    }
  } /* end loop over neighbors */

  /* free rqs in par of 1st var */
  if(rqs)
  {
    free(rqs);
    VarA_(node, vi0)->par = NULL;
  }
}


/* free req and send arrays after all has been sent */
void free_dat_gcom_reqs_after_Waitall_com_send(tNode *node)
{
  tDat *dat = node->dat;

  if(!dat) return;

  /* to be sure, wait again for all recvs */
  nMPI_Waitall_com_recv(dat->gcom);

  /* wait until all has been sent, then free all buffers for this face */
  nMPI_Waitall_com_send(dat->gcom);
  realloc_com_reqs(dat->gcom, 0); /* free req and send arrays */
}

/* get ghostdata for all nodes out of buffers and free the buffers */
void get_all_myln_ghostdata_for_vl(tMesh *mesh, tVarList  *vl)
{
  TIMER_START;

  /* If we want threads in this loop, we need MPI_Init_thread with
     MPI_THREAD_MULTIPLE, instead of just MPI_Init in main. */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    get_ghostdata_for_vl(node, vl);
  }

  /* postpone Waitall until we have finished all nodefaces. This could have
     been already called in get_all_ghost_for_vl to free mem earlier. */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    free_dat_gcom_reqs_after_Waitall_com_send(node);
  }

  TIMER_STOP;
}
