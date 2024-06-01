/* main/amr/indicators.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 0


/* get global pars for amr */
extern tAMR amr[1];


/* functions to exchange ghost zone data on a uniform mesh */


/***************************************************************************/
/* neighbor indices */
/***************************************************************************/

/* indices of nearest neighbors on faces 0-5 */
enum
{
  f_mcc=0, f_pcc, f_cmc, f_cpc, f_ccm, f_ccp // nearest neighbors on faces 0-5
};

/* indices of node (ccc), nearest neighbors, next-nearest and
   next to next-nearest neighbors */
enum
{
  mmm, cmm, pmm,
  mcm, ccm, pcm,
  mpm, cpm, ppm,

  mmc, cmc, pmc,
  mcc, ccc, pcc,
  mpc, cpc, ppc,

  mmp, cmp, pmp,
  mcp, ccp, pcp,
  mpp, cpp, ppp
};
/* Here we have:
  mcc, pcc, cmc, cpc, ccm, ccp,   // nearest neighbors on faces 0-5
  mmc, mpc, pmc, ppc,             // next-nearest nbs on face 0,1 in Y-dir
  mcm, pcm, cmm, cpm,             // next-nearest nbs on face 4
  mcp, pcp, cmp, cpp,             // next-nearest nbs on face 5
  mmm, pmm, mpm, ppm,             // next to next-nearest nbs on face 4
  mmp, pmp, mpp, ppp              // next to next-nearest nbs on face 5 */


/***************************************************************************/
/* find neighbors */
/***************************************************************************/

/* find nearest + next-nearest + next-next-nearest neighbors */
void find_nb27(tNode *node, tNode *nb27[27])
{
  /* self */
  nb27[ccc] = node;

  /* nearest neighbors */
  nb27[mcc] = node->fnb[f_mcc] ? node->fnb[f_mcc][0] : NULL;
  nb27[pcc] = node->fnb[f_pcc] ? node->fnb[f_pcc][0] : NULL;
  nb27[cmc] = node->fnb[f_cmc] ? node->fnb[f_cmc][0] : NULL;
  nb27[cpc] = node->fnb[f_cpc] ? node->fnb[f_cpc][0] : NULL;
  nb27[ccm] = node->fnb[f_ccm] ? node->fnb[f_ccm][0] : NULL;
  nb27[ccp] = node->fnb[f_ccp] ? node->fnb[f_ccp][0] : NULL;

  /* next-nearest nbs on face 0,1 */
  if(nb27[mcc])
  {
    nb27[mmc] = nb27[mcc]->fnb[f_cmc] ? nb27[mcc]->fnb[f_cmc][0] : NULL;
    nb27[mpc] = nb27[mcc]->fnb[f_cpc] ? nb27[mcc]->fnb[f_cpc][0] : NULL;
  }
  else
  {
    nb27[mmc] = nb27[mpc] = NULL;
  }
  if(nb27[pcc])
  {
    nb27[pmc] = nb27[pcc]->fnb[f_cmc] ? nb27[pcc]->fnb[f_cmc][0] : NULL;
    nb27[ppc] = nb27[pcc]->fnb[f_cpc] ? nb27[pcc]->fnb[f_cpc][0] : NULL;
  }
  else
  {
    nb27[pmc] = nb27[ppc] = NULL;
  }

  /* next-nearest nbs on face 4 */
  if(nb27[ccm])
  {
    nb27[mcm] = nb27[ccm]->fnb[f_mcc] ? nb27[ccm]->fnb[f_mcc][0] : NULL;
    nb27[pcm] = nb27[ccm]->fnb[f_pcc] ? nb27[ccm]->fnb[f_pcc][0] : NULL;
    nb27[cmm] = nb27[ccm]->fnb[f_cmc] ? nb27[ccm]->fnb[f_cmc][0] : NULL;
    nb27[cpm] = nb27[ccm]->fnb[f_cpc] ? nb27[ccm]->fnb[f_cpc][0] : NULL;
  }
  else
  {
    nb27[mcm] = nb27[pcm] = nb27[cmm] = nb27[cpm] = NULL;
  }

  /* next-nearest nbs on face 5 */
  if(nb27[ccp])
  {
    nb27[mcp] = nb27[ccp]->fnb[f_mcc] ? nb27[ccp]->fnb[f_mcc][0] : NULL;
    nb27[pcp] = nb27[ccp]->fnb[f_pcc] ? nb27[ccp]->fnb[f_pcc][0] : NULL;
    nb27[cmp] = nb27[ccp]->fnb[f_cmc] ? nb27[ccp]->fnb[f_cmc][0] : NULL;
    nb27[cpp] = nb27[ccp]->fnb[f_cpc] ? nb27[ccp]->fnb[f_cpc][0] : NULL;
  }
  else
  {
   nb27[mcp] = nb27[pcp] = nb27[cmp] = nb27[cpp] = NULL;
  }

  /* next to next-nearest nbs on face 4 */
  if(nb27[cmm])
  {
    nb27[mmm] = nb27[cmm]->fnb[f_mcc] ? nb27[cmm]->fnb[f_mcc][0] : NULL;
    nb27[pmm] = nb27[cmm]->fnb[f_pcc] ? nb27[cmm]->fnb[f_pcc][0] : NULL;
  }
  else
  {
    nb27[mmm] = nb27[pmm] = NULL;
  }
  if(nb27[cpm])
  {
    nb27[mpm] = nb27[cpm]->fnb[f_mcc] ? nb27[cpm]->fnb[f_mcc][0] : NULL;
    nb27[ppm] = nb27[cpm]->fnb[f_pcc] ? nb27[cpm]->fnb[f_pcc][0] : NULL;
  }
  else
  {
    nb27[mpm] = nb27[ppm] = NULL;
  }

  /* next to next-nearest nbs on face 5 */
  if(nb27[cmp])
  {
    nb27[mmp] = nb27[cmp]->fnb[f_mcc] ? nb27[cmp]->fnb[f_mcc][0] : NULL;
    nb27[pmp] = nb27[cmp]->fnb[f_pcc] ? nb27[cmp]->fnb[f_pcc][0] : NULL;
  }
  else
  {
    nb27[mmp] = nb27[pmp] = NULL;
  }
  if(nb27[cpp])
  {
    nb27[mpp] = nb27[cpp]->fnb[f_mcc] ? nb27[cpp]->fnb[f_mcc][0] : NULL;
    nb27[ppp] = nb27[cpp]->fnb[f_pcc] ? nb27[cpp]->fnb[f_pcc][0] : NULL;
  }
  else
  {
    nb27[mpp] = nb27[ppp] = NULL;
  }

  if(PR)
  {
    int ni;
    char s[100];
    int nn[] = { 3,3,3 }; /* pcc - mcc = 3 */
    int i,j,k;


    PRF;printf(": %s\n", nodename(node, s,99));
    for(ni=0; ni<27; ni++)
    {
      k = kOfInd_n(ni, nn);
      j = jOfInd_n_k(ni, nn,k);
      i = iOfInd_n_jk(ni, nn,j,k);
      printf("ni=%d:=%+d%+d%+d: %s\n", ni, i-1,j-1,k-1,
             nodename(nb27[ni], s,99));
    }
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
                             int sta[3], int nc[3])
{
  int nn[] = { 3,3,3 }; /* pcc - mcc = 3 */
  int n2gho = nghosts*2;
  int i,j,k;

  /* select start values for copy */
  k = kOfInd_n(ni, nn);
  j = jOfInd_n_k(ni, nn,k);
  i = iOfInd_n_jk(ni, nn,j,k);

  gsta[0] = nghosts;
  gsta[1] = nghosts;
  gsta[2] = nghosts;
  sta[0] = nghosts;
  sta[1] = nghosts;
  sta[2] = nghosts;

  /* cases m??, ?m?, ??,m */
  if(i==0)
  {
    gsta[0] = 0;
    sta[0] = nghosts;
  }
  if(j==0)
  {
    gsta[1] = 0;
    sta[1] = nghosts;
  }
  if(k==0)
  {
    gsta[2] = 0;
    sta[2] = nghosts;
  }

  /* cases p??, ?p?, ??,p */
  if(i==2)
  {
    gsta[0] = n[0] - nghosts;
    sta[0] = n[0] - n2gho;
  }
  if(j==2)
  {
    gsta[1] = n[1] - nghosts;
    sta[1] = n[1] - n2gho;
  }
  if(k==2)
  {
    gsta[2] = n[2] - nghosts;
    sta[2] = n[2] - n2gho;
  }

  /* select thinkness of layer to copy */
  nc[0] = n[0] - n2gho;
  nc[1] = n[1] - n2gho;
  nc[2] = n[2] - n2gho;

  if(i!=1) nc[0] = nghosts;
  if(j!=1) nc[1] = nghosts;
  if(k!=1) nc[2] = nghosts;

  if(PR)
  {
    PRF;printf(": ni=%d=%+d%+d%+d: %d %d %d, %d %d %d, %d %d %d\n",
               ni,i-1,j-1,k-1,
               gsta[0],gsta[1],gsta[2], sta[0],sta[1],sta[2], nc[0],nc[1],nc[2]);
  }
}


/***************************************************************************/
/* put ghost data into buffers and start MPI send/recv to get data */
/***************************************************************************/

/* setup ghost exchange */
void request_ghostdata_for_vl(tNode *node, tVarList  *vl)
{
  tMesh *mesh = node->pat->mesh;
  int nghosts;
  tNode *nb27[27];
  tDat *dat = node->dat;
  int *rqs;
  int vi0 = Vind(vl,0); /* 1st var */
  int ni;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  nghosts = Geti(amr->nghosts);

  /* alloc memory to store request numbers for all 27 nbs */
  if(VarA_(node, vi0)->par) free(VarA_(node, vi0)->par);
  rqs = imalloc(27);
  VarA_(node, vi0)->par = (void *) rqs;

  /* find 27 neighbors */
  find_nb27(node, nb27);

  /* loop over 26 neighbors */
  for(ni=0; ni<27; ni++)
  {
    tNode *nb = nb27[ni];
    int *n = node->n;
    int *nb_n = nb->n;
    int nb_ni, vli;
    int gsta[3], sta[3], nc[3], nb_gsta[3], nb_sta[3],  nb_nc[3];

    /* do nothing for self */
    if(ni==ccc) continue;

    /* goto next neighbor if nb is NULL */
    if(!nb) continue;

    /* select start value and thickness of layer to copy */
    set_ghoststart_start_nc(ni, n, nghosts, gsta, sta, nc);

    /* is nb local? */
    if(nb->dat)
    {
      /* nb is local so just copy */
      nb_ni = get_nb_ni(ni);
      set_ghoststart_start_nc(nb_ni, nb_n, nghosts, nb_gsta, nb_sta, nb_nc);
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
      int buflen = (vl->n) * nc[0]*nc[1]*nc[2];
      long s_ltag, r_ltag;
      int rq, nb_rank, s_tag, r_tag, ci;
      nMPI_Comm s_comm, r_comm;
      tCom *com = dat->gcom;
      double *sbuf, *rbuf; /* buffers for MPI */
      int si, lid, nb_lid;

      /* find my index in nb */
      nb_ni = get_nb_ni(ni);

      /* set up tags */
      nb_rank = nb->datrank;
      lid = calc_node_lid(node);
      nb_lid = calc_node_lid(nb);
      s_ltag = lid*27 + ni;
      r_ltag = nb_lid*27 + nb_ni;
      nMPI_long_tag_to_commi_tag(s_ltag, &ci, &s_tag);
      s_comm = nMPIvars_get_comm(ci);
      nMPI_long_tag_to_commi_tag(r_ltag, &ci, &r_tag);
      r_comm = nMPIvars_get_comm(ci);

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
      MCK(
      nMPI_Isend_Irecv_double_com(com, rq, nb_rank, s_tag,r_tag,
                                  s_comm,r_comm)
      );
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

/* call request_all_myln_ghostdata_for_vl for all EvoVars */
void request_all_myln_ghostdata(tMesh *mesh)
{
  int vi;
  tVarList *vl = vlalloc(mesh);

  /* make varlist with all non-Aux. Vars */
  for(vi=0; vi<mesh->nvdb; vi++)
    if(MeshVarType(mesh, vi) != AUXVAR) vlpushone(vl, vi);

  request_all_myln_ghostdata_for_vl(mesh, vl);
  vlfree(vl);
}


/**********************************************************************/
/* funcs to free stuff like MPI com */
/**********************************************************************/

/* free req and send arrays after all has been sent */
void free_dat_gcom_reqs_after_Waitall_com_send(tNode *node)
{
  tDat *dat = node->dat;

  if(!dat) return;

  /* to be sure, wait again for all recvs */
  MCK( nMPI_Waitall_com_recv(dat->gcom) );

  /* wait until all has been sent, then free all buffers for this face */
  MCK( nMPI_Waitall_com_send(dat->gcom) );
  realloc_com_reqs(dat->gcom, 0); /* free req and send arrays */
}


/**********************************************************************/
/* get the ghost data out of the MPI buffers */
/**********************************************************************/

/* get all ghost data from neighbors */
void get_ghostdata_for_vl(tNode *node, tVarList  *vl)
{
  tMesh *mesh = node->pat->mesh;
  int nghosts;
  tNode *nb27[27];
  tDat *dat = node->dat;
  tCom *com = dat->gcom;
  int ni, nvars, rq;
  int vi0 = Vind(vl,0); /* 1st var */
  int *rqs;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  /* number of vars that exchanged ghost */
  nvars = vl->n;

  /* get MPI rq numbers */
  rqs = (int *) VarA_(node, vi0)->par;

  /* do something only if we have requests and vars */
  if(com->n_rq != 0 && nvars != 0)
  {
    nghosts = Geti(amr->nghosts);

    /* find 27 neighbors */
    find_nb27(node, nb27);

    /* loop over 26 neighbors */
    for(ni=0; ni<27; ni++)
    {
      tNode *nb = nb27[ni];
      int *n = node->n;
      double *rbuf; /* buffer for recv */
      int si, vli;
      int gsta[3], sta[3], nc[3];

      /* do nothing for self */
      if(ni==ccc) continue;

      /* goto next neighbor if nb is NULL */
      if(!nb) continue;

      /* if nb is local we have already exchanged info  */
      if(nb->dat) continue;

      /* nb is on other process, we have used MPI to exchange data */


      /* select start value and thickness of layer to copy */
      set_ghoststart_start_nc(ni, n, nghosts, gsta, sta, nc);

      /* get MPI request number, and buffer */
      rq = rqs[ni];
      rbuf = get_com_recv_buf(com, rq);

      /* wait for MPI buffer */
      MCK( nMPI_Wait_com_recv(com, rq) );

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
  }

  /* free rqs in par of 1st var */
  if(rqs)
  {
    free(rqs);
    VarA_(node, vi0)->par = NULL;
  }
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

/* call get_all_myln_ghostdata_for_vl for all EvoVars */
void get_all_myln_ghostdata(tMesh *mesh)
{
  int vi;
  tVarList *vl = vlalloc(mesh);

  /* make varlist with all non-Aux. Vars */
  for(vi=0; vi<mesh->nvdb; vi++)
    if(MeshVarType(mesh, vi) != AUXVAR) vlpushone(vl, vi);

  get_all_myln_ghostdata_for_vl(mesh, vl);
  vlfree(vl);
}
