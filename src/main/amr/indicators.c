/* main/amr/indicators.c */
/* Wolfgang Tichy, 4/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 0

/* functions to exchange indicator data */

/**********************************************************************/
/* allocate and fill indcs for vars that need it */
/**********************************************************************/
/* empty indc that we need to fill in */
tIndic *alloc_empty_indc(tNode *node)
{
  int f;
  tIndic *ic = calloc(1, sizeof(*ic));

  for(f=0; f<6; f++)
  {
    int nnb = node->nfnb[f];
    if(nnb)
    {
      ic->nbindc[f] = calloc(nnb, sizeof(ic->nbindc[0][0]));
      ic->allocd_nbindc[f] = calloc(nnb, sizeof(ic->allocd_nbindc[0][0]));
    }
  }

  return ic;
}

/* free all we need to in a indc */
void free_indc(tIndic *ic)
{
  tNode *node;
  int f, ni;

  if(!ic) return;
  node = ic->dat->node;

  /* free content of lists */

  /* free myindc */
  free_array(ic->myindc);

  /* free nbindc[f][i] only if it is allocd  */
  for(f=0; f<6; f++)
  {
    int nnb = node->nfnb[f];

    for(ni=0; ni<nnb; ni++)
      if(ic->allocd_nbindc[f][ni]) free_array(ic->nbindc[f][ni]);

    /* free lists */
    free(ic->nbindc[f]);
    free(ic->allocd_nbindc[f]);
  }

  /* free indc */
  free(ic);
}

/* free indc on node */
void free_all_indc_for_vl(tNode *node, tVarList  *vl)
{
  tDat *dat = node->dat;
  int vi;

  if(!dat) return;

  /* free all indc */
  for(vi=0; vi<dat->nv; vi++)
  {
    free_indc(dat->ic[vi]);
    dat->ic[vi] = NULL;
  }
}

/* free all indc on all nodes in the mesh */
void free_all_myln_indc_for_vl(tMesh *mesh, tVarList  *vl)
{
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    free_all_indc_for_vl(node, vl);
  }
}


/* initialize a indc for var vi at face with nnb neighbors */
tIndic *init_myindc_NULL(tNode *node, int vi, int nvals)
{
  tDat *dat = node->dat;
  tIndic *ic;

  /* do nothing if no data on this node */
  if(!dat) return NULL;

  /* do nothing if var is not enabled */
  if(!dat->v[vi]) return NULL;

  /* prep. */
  ic = alloc_empty_indc(node);
  ic->dat = dat;
  ic->vi = vi;
  ic->nvals = nvals;

  /* do not allocate my indc array yet */
  ic->myindc = NULL;

  return ic;
}

/* init indicators of a node for a varlist, allowing for nvals values
   used as indicators, e.g. 3 if 3 indic. vals: min, max, av */
int init_myindc_for_vl(tNode *node, tVarList  *vl, int nvals)
{
  tDat *dat = node->dat;
  int vli, vi;
  tArray *myindc0; /* segmented array containing myindc for all vi */

  if(!dat) return 0;

  for(vli=0; vli<vl->n; vli++)
  {
    vi = vl->index[vli];
    free_indc(dat->ic[vi]);
    dat->ic[vi] = init_myindc_NULL(node, vi, nvals);

    /* allocate segmented array for indc of first needed var */
    if(vli==0) myindc0 = alloc_array1d_with_segs(nvals, 0, vl->n);

    /* set my indc array with mem in myindc0*/
    dat->ic[vi]->myindc = get_array_seg(myindc0, vli);
  }
  return vli;
}

/* init all indc on all nodes in the mesh for a varlist */
void init_all_myln_myindc_for_vl(tMesh *mesh, tVarList  *vl, int nvals)
{
  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    init_myindc_for_vl(node, vl, nvals);
  }
}


/***************************************************************************/
/* put myindc into buffers and start MPI send/recv to get data into nbindc */
/***************************************************************************/

/* get all indcs from neighbors */
void request_indc_exchange_for_vl(tNode *node, tVarList  *vl)
{
  tNode *nb;
  tDat *dat = node->dat;
  int f, ni;
  int vli, vi;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  /* loop over all neighbors in fnb */
  for(f=0; f<6; f++)
  {
    for(ni=0; ni<node->nfnb[f]; ni++)
    {
      nb = node->fnb[f][ni];

      /* is nb local? */
      if(nb->dat)
      {
        /* nb is local so just point ic->nbindc[f][ni] to its data */
        for(vli=0; vli<vl->n; vli++)
        {
          tIndic *my_ic;
          tIndic *nb_ic;

          vi = vl->index[vli];
          my_ic = dat->ic[vi];
          nb_ic = nb->dat->ic[vi];

          if(!my_ic) continue;
          if(nb_ic) my_ic->nbindc[f][ni] = nb_ic->myindc;
        }
      }
      else
      {
        /* nb is on other process so use MPI to exchange data */
        long s_ltag, r_ltag;
        int rq, nb_rank, s_tag, r_tag, ci;
        nMPI_Comm s_comm, r_comm;
        tCom *com = dat->icom;
        int nvals;
        double *sbuf, *rbuf; /* buffers for MPI */
        tArray *nbindc_f_ni; /* segmented array with nbindc[f][ni] for all vi */
        tIndic *ic;
        int nvars = vl->n;
        int si;
        int found, nb_f, nb_ni, lid, nb_lid;

        /* we get nb_f, nb_ni only so that we can make unique tags */
        found = locate_facenb_in_fnbs(nb, node, &nb_f, &nb_ni);
        if(!found) errorexit("couldn't find nb face!!!");

        /* use MPI to recv nb->dat->ic[vi]->myindc in 
           dat->ic[vi]->nbindc[f][ni], and also send dat->ic[vi]->myindc to 
           nb->dat->ic[vi]->nbindc[nb_nbi] */
        nb_rank = nb->datrank;
        //r_tag = (node->nid);
        //s_tag = (nb->nid);
        //r_comm = nb->comm;
        //s_comm = node->comm;
        lid = calc_node_lid(node);
        nb_lid = calc_node_lid(nb);
        s_ltag = (nb_lid*64 + nb_ni)*6 + nb_f;
        r_ltag = (lid*64 + ni)*6 + f;
        nMPI_long_tag_to_commi_tag(s_ltag, &ci, &s_tag);
        s_comm = nMPIvars_get_comm(ci);
        nMPI_long_tag_to_commi_tag(r_ltag, &ci, &r_tag);
        r_comm = nMPIvars_get_comm(ci);

        /* alloc one segmented array for nbindc of all nvars variables needed */
        vi = vl->index[0]; /* first vi */
        ic = dat->ic[vi];  /* indc of 1st var that needs it */
        if(!ic) errorexit("1st var that needs indc, has no indc mem!");
        nvals = ic->nvals;
        if(!ic->nbindc[f][ni])
          nbindc_f_ni = alloc_array1d_with_segs(nvals, 0, nvars);
        else
          nbindc_f_ni = ic->nbindc[f][ni];

        /* point send and recv buffers to data in arrays */
        sbuf = ic->myindc->d;  /* use segmented array as sbuf */
        rbuf = nbindc_f_ni->d; /* use segmented array as rbuf */

        /* save buffers in com, and set flag to not free them, because
           they get freed when we free the indicator arrays */
        rq = append_buffers_to_com(com, sbuf,nvars*nvals, rbuf,nvars*nvals);
        set_free_buf_in_com(com, 0);

        /* point nbindc to correct place */
        for(si=0, vli=0; vli<vl->n; vli++)
        {
          vi = vl->index[vli];
          ic = dat->ic[vi];
          /* do nothing if there is no ic */
          if(ic && dat->v[vi])
          {
            /* point indc data to nbindc_f_ni=rbuf to later recv
               neighbor data */
            if(!ic->nbindc[f][ni])
            {
              ic->nbindc[f][ni] = get_array_seg(nbindc_f_ni, si);
              si++;                         // inc segment index
              ic->allocd_nbindc[f][ni] = 1; // flag that we allocd
            }

            /* save MPI request number in the array */
            ic->nbindc[f][ni]->info = rq;
          }
        }
        /* now call MPI */
        nMPI_Isend_Irecv_double_com(com, rq, nb_rank, s_tag,r_tag,
                                    s_comm,r_comm);
      }
    }
  } /* end loop over neighbors */
}


/* request indc exchanges on all my nodes in the mesh
   Note: We need to call this! If we call request_all_indc_exchange(n1)
   for only node n1, MPI deadlocks because the other nodes are not sending
   to n1 or receiving from n1 */
void request_all_myln_indc_exchange_for_vl(tMesh *mesh, tVarList  *vl)
{

  TIMER_START;

  /* If we want threads in this loop, we need MPI_Init_thread with
     MPI_THREAD_MULTIPLE, instead of just MPI_Init in main. */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    request_indc_exchange_for_vl(node, vl);
  }
  TIMER_STOP;
}


/**********************************************************************/
/* get the nbindc data out of the MPI buffers */
/**********************************************************************/

/* get all indc from neighbor with index ni at face */
void get_indc_for_vl(tNode *node, tVarList  *vl, int face, int ni)
{
  tNode *nb = node->fnb[face][ni];
  tDat *dat = node->dat;
  tCom *com = dat->icom;
  int vi, nvars, rq;

  /* do nothing if this node is on other proc */
  if(!dat) return;

  /* do nothing if com is empty */
  if(com->n_rq == 0) return;

  /* if nb is local we have already exchanged info  */
  if(nb->dat) return;

  /* nb is on other process, we have used MPI to exchange data */

  /* number of vars that exchanged indc */
  nvars = vl->n;

  /* do nothing if there are no vars that exchanged indc */
  if(!nvars) return;

  /* get MPI request number */
  vi = vl->index[0]; // first var in vl
  rq = dat->ic[vi]->nbindc[face][ni]->info;

  /* wait for MPI buffer */
  nMPI_Wait_com_recv(com, rq);
}


/* free req and send arrays after all has been sent */
void free_dat_icom_reqs_after_Waitall_com_send(tNode *node)
{
  tDat *dat = node->dat;

  if(!dat) return;

  /* to be sure, wait again for all recvs */
  nMPI_Waitall_com_recv(dat->icom);

  /* wait until all has been sent, then free all buffers for this face */
  nMPI_Waitall_com_send(dat->icom);
  realloc_com_reqs(dat->icom, 0); /* free req and send arrays */
}

/* get nbindc from all faces and variables for this node out of buffers */
void get_all_indc_for_vl(tNode *node, tVarList  *vl)
{
  int face, ni;
  tDat *dat = node->dat;

  if(!dat) return;

  /* get nbindc for each face and neighbor */
  for(face=0; face<6; face++)
    for(ni=0; ni<node->nfnb[face]; ni++)
      get_indc_for_vl(node, vl, face, ni);
}

/* get nbindc for all nodes out of buffers and free the buffers */
void get_all_myln_indc_for_vl(tMesh *mesh, tVarList  *vl)
{
  TIMER_START;

  /* If we want threads in this loop, we need MPI_Init_thread with
     MPI_THREAD_MULTIPLE, instead of just MPI_Init in main. */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    get_all_indc_for_vl(node, vl);
  }

  /* postpone Waitall until we have finished all nodefaces. This could have
     been already called in get_all_indc_for_vl to free mem earlier. */
  formylnodes_noomp(mesh)
  {
    tNode *node = MyLnode;
    free_dat_icom_reqs_after_Waitall_com_send(node);
  }

  TIMER_STOP;
}
