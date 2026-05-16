/* amr.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "amr.h"

#define PR 0


/* use gridpoints from basis/gridpoints.c */
extern tGridPoints gridpoints[1];

/*************************************************************************/
/* global pars */
/*************************************************************************/

/* global pars for amr */
tAMR amr[1];

/* func to init frequently used pars */
int amr_init_global_pars(tMesh *mesh)
{
  int amr_N0 = Par("amr_N0");
  int amr_N1 = Par("amr_N1");
  int amr_N2 = Par("amr_N2");
  int d;

  PRFs(":\n");

  /* set amr globals */
  amr->elm_nbinfo0 = Ind("amr_elm_nbinfo0");

  /* set amr globals */
  amr->dir_active[0] = Par("amr_dir_active0");
  amr->dir_active[1] = Par("amr_dir_active1");
  amr->dir_active[2] = Par("amr_dir_active2");
  amr->sibl1to7_weight = Par("amr_sibl1to7_weight");
  amr->MPIexchange = Par("amr_MPIexchange");
  amr->nghosts     = Par("amr_nghosts");
  amr->nbsearch_n  = Par("amr_nbsearch_n");
  amr->Lagrange_interp_order = Par("amr_Lagrange_interp_order");
  amr->WENO_interp_order     = Par("amr_WENO_interp_order");

  /* print global vars and pars */
  printf(" amr->elm_nbinfo0 = var_%04d : VarName(amr->elm_nbinfo0) = %s\n",
         amr->elm_nbinfo0, VarName(amr->elm_nbinfo0));

  printf(" amr->dir_active =\n");
  for(d=0; d<3; d++)
    printf("      par_%04d : Getb(amr->dir_active[%d]) = %d\n",
           amr->dir_active[d], d, Getb(amr->dir_active[d]));
  printf(" amr->sibl1to7_weight = par_%04d"
         " : Getd(amr->sibl1to7_weight) = %g\n",
         amr->sibl1to7_weight, Getd(amr->sibl1to7_weight));
  printf(" amr->MPIexchange = par_%04d : Geti(amr->MPIexchange) = %d\n",
         amr->MPIexchange, Geti(amr->MPIexchange));
  printf(" amr->nghosts = par_%04d :     Geti(amr->nghosts) = %d\n",
         amr->nghosts, Geti(amr->nghosts));
  printf(" amr->nbsearch_n = par_%04d :  Geti(amr->nbsearch_n) = %d\n",
         amr->nbsearch_n, Geti(amr->nbsearch_n));
  printf(" amr->Lagrange_interp_order=par_%04d : Geti(amr->Lagrange_interp_order) = %d\n",
         amr->Lagrange_interp_order, Geti(amr->Lagrange_interp_order));
  printf(" amr->WENO_interp_order = par_%04d :   Geti(amr->WENO_interp_order) = %d\n",
         amr->WENO_interp_order, Geti(amr->WENO_interp_order));

  /* set amr_n0,... from amr_N0, ... */
  if(Geti(amr_N0)>0)
  {
    int nghosts = Geti(amr->nghosts);
    int lmax = ceil( log(nMPI_size())/log(8.) );
    int n0i = Geti(amr_N0)/pow(2.,lmax) + 0.5;
    int n1i = Geti(amr_N1)/pow(2.,lmax) + 0.5;
    int n2i = Geti(amr_N2)/pow(2.,lmax) + 0.5;

    /* set some amr pars */
    Seti(Par("amr_luni"), lmax);
    Seti(Par("amr_n0"), n0i + 2*nghosts);
    Seti(Par("amr_n1"), n1i + 2*nghosts);
    Seti(Par("amr_n2"), n2i + 2*nghosts);
    PRF;printf(": Setting:\n");
    printf("  amr_luni = %d\n", Geti(Par("amr_luni")));
    printf("  amr_n0 = %d\n", Geti(Par("amr_n0")));
    printf("  amr_n1 = %d\n", Geti(Par("amr_n1")));
    printf("  amr_n2 = %d\n", Geti(Par("amr_n2")));
  }
  return 0;
}

/* print some info about things in amr */
int amr_print_thread_info(tMesh *mesh)
{
  PRFs(":\n");
  system2("echo", "$OMP_NUM_THREADS");
  printf("MAX_NTHREADS = %d\n", MAX_NTHREADS);
#ifdef USEOMP
  printf("OpenMP pragmas are active.\n");
  printf("  OMP_VERSION = _OPENMP = %d\n", OMP_VERSION);
#ifdef USE_OMP_TASKLOOP
  printf("  USE_OMP_TASKLOOP is defined.\n");
#endif
#else
  printf("OpenMP pragmas are off.\n");
#endif
  return 0;
}

/* set patch group flags in tPat */
void amr_set_patgroup_npg_pg0(tMesh *mesh, int npg, int pg0)
{
  int off;
  for(off=0; off<npg; off++)
  {
    int p = pg0 + off;
    tPat *pat = mesh->pat[p];
    pat->npg = npg;
    pat->pg0 = pg0;
  }
}

/* zero patch group info for all patches on mesh */
int amr_zero_all_patgroup_npg_pg0(tMesh *mesh)
{
  int pi;
  forpatches(mesh, pi)
  {
    tPat *pat = mesh->pat[pi];
    pat->npg = pat->pg0 = 0;
  }
  return 0;
}


/*************************************************************************/
/* funcs for MPI exchange */
/*************************************************************************/

/* init exchange. Could be replaced by:
   formyelms(mesh) MPIexchange_init(MyElm); */
void MPIexchange_init_all_myln(tMesh *mesh)
{
  if(PR) PRFs(":\n");

  switch(Geti(amr->MPIexchange))
  {
  case 1:
    init_all_myln_surfaces(mesh);
    break;
  case 2:
    break;
  default:
    errorexit("unknown value in amr_MPIexchange");
  }
}
/* init exchange on just one elm */
void MPIexchange_init(tElm *elm)
{
  tMesh *mesh = Elm_mesh(elm);
  if(PR) PRFs(":\n");
  switch(Geti(amr->MPIexchange))
  {
  case 1:
    init_all_surfaces(elm);
    break;
  case 2:
    break;
  default:
    errorexit("unknown value in amr_MPIexchange");
  }
}

/* set some local data. Could be replaced by:
   formyelms(mesh) MPIexchange_set_localdata(MyElm); */
void MPIexchange_set_all_myln_localdata(tMesh *mesh)
{
  if(PR) PRFs(":\n");

  switch(Geti(amr->MPIexchange))
  {
  case 1:
    set_all_myln_mysurf(mesh);
    break;
  case 2:
    break;
  default:
    errorexit("unknown value in amr_MPIexchange");
  }
}
/* set some local data on just one elm */
void MPIexchange_set_localdata(tElm *elm)
{
  tMesh *mesh = Elm_mesh(elm);
  if(PR) PRFs(":\n");
  switch(Geti(amr->MPIexchange))
  {
  case 1:
    set_all_mysurf(elm);
    break;
  case 2:
    break;
  default:
    errorexit("unknown value in amr_MPIexchange");
  }
}

/* request exchange */
void MPIexchange_request_all_myln_data(tMesh *mesh)
{
  if(PR) PRFs(":\n");

  switch(Geti(amr->MPIexchange))
  {
  case 1:
    request_all_myln_surfaces_exchange(mesh);
    break;
  case 2:
    request_all_myln_ghostdata(mesh);
    break;
  default:
    errorexit("unknown value in amr_MPIexchange");
  }
}

/* get buffers */
void MPIexchange_get_all_myln_data(tMesh *mesh)
{
  if(PR) PRFs(":\n");

  switch(Geti(amr->MPIexchange))
  {
  case 1:
    get_all_myln_surfaces(mesh);
    break;
  case 2:
    get_all_myln_ghostdata(mesh);
    break;
  default:
    errorexit("unknown value in amr_MPIexchange");
  }
}

/* clean up after MPIexchange */
void MPIexchange_free_all_myln(tMesh *mesh)
{
  if(PR) PRFs(":\n");

  switch(Geti(amr->MPIexchange))
  {
  case 1:
    free_all_myln_surfaces(mesh);
    break;
  case 2:
    break;
  default:
    errorexit("unknown value in amr_MPIexchange");
  }
}


/*************************************************************************/
/* various untilities */
/*************************************************************************/

/* index, but reduced by dimension along norm */
int Ind_n_norm(int i, int j, int k, int n[3], int norm)
{
  int N0,N1, II,JJ,KK;

  N0 = n[0];
  N1 = n[1];
  II = i;
  JJ = j;
  KK = k;

  switch(norm)
  {
  case 0:
    N0 = 1;
    II = 0;
    break;
  case 1:
    N1 = 1;
    JJ = 0;
    break;
  case 2:
    KK = 0;
    break;
  default:
    errorexit("norm needs to be 0,1,2");
  }
  return II + N0*(JJ + N1*KK);
}


/* get Xb in dir for a node out of gridpoints */
tArray *node_Xb(tNode *node, int dir)
{
  return gridpoints->Xb[node->pt_typ[dir]][node->n[dir]];
}
/* get Xb[3] for pt_typ and n out of gridpoints */
void Xb3_pt_typ_n(int pt_typ[3], int n[3], tArray *Xb[3])
{
  int d;
  for(d=0; d<3; d++) Xb[d] = gridpoints->Xb[pt_typ[d]][n[d]];
}
/* get Xb[3] for n (assuming point type of node) out of gridpoints */
void Xb3_n(tNode *node, int n[3], tArray *Xb[3])
{
  int d;
  for(d=0; d<3; d++) Xb[d] = gridpoints->Xb[node->pt_typ[d]][n[d]];
}
/* get Xb[3] for a node out of gridpoints */
void node_Xb3(tNode *node, tArray *Xb[3])
{
  Xb3_n(node, node->n, Xb);
}

/* get Wq in dir for a node out of gridpoints */
tArray *node_Wq(tNode *node, int dir)
{
  return gridpoints->Wq[node->pt_typ[dir]][node->n[dir]];
}
/* get Wq[3] for n (assuming point type of node) out of gridpoints */
void Wq3_n(tNode *node, int n[3], tArray *Wq[3])
{
  int d;
  for(d=0; d<3; d++) Wq[d] = gridpoints->Wq[node->pt_typ[d]][n[d]];
}
/* get Wq[3] for a node out of gridpoints */
void node_Wq3(tNode *node, tArray *Wq[3])
{
  Wq3_n(node, node->n, Wq);
}

/* get WL in dir for a node out of gridpoints */
tArray *node_WL(tNode *node, int dir)
{
  return gridpoints->WL[node->pt_typ[dir]][node->n[dir]];
}
/* get WL[3] for n (assuming point type of node) out of gridpoints */
void WL3_n(tNode *node, int n[3], tArray *WL[3])
{
  int d;
  for(d=0; d<3; d++) WL[d] = gridpoints->WL[node->pt_typ[d]][n[d]];
}
/* get WL[3] for a node out of gridpoints */
void node_WL3(tNode *node, tArray *WL[3])
{
  WL3_n(node, node->n, WL);
}

/* get Dt in dir for a node out of gridpoints */
tArray *node_Dt(tNode *node, int dir)
{
  return gridpoints->Dt[node->pt_typ[dir]][node->n[dir]];
}
/* get Dt[3] for n (assuming point type of node) out of gridpoints */
void Dt3_n(tNode *node, int n[3], tArray *Dt[3])
{
  int d;
  for(d=0; d<3; d++) Dt[d] = gridpoints->Dt[node->pt_typ[d]][n[d]];
}
/* get Dt[3] for a node out of gridpoints */
void node_Dt3(tNode *node, tArray *Dt[3])
{
  Dt3_n(node, node->n, Dt);
}

/* get Dpt in dir for a node out of gridpoints */
tArray *node_Dpt(tNode *node, int dir)
{
  return gridpoints->Dpt[node->pt_typ[dir]][node->n[dir]];
}

/* get Dmt in dir for a node out of gridpoints */
tArray *node_Dmt(tNode *node, int dir)
{
  return gridpoints->Dmt[node->pt_typ[dir]][node->n[dir]];
}

/* get At in dir for a node out of gridpoints */
tArray *node_At(tNode *node, int dir)
{
  return gridpoints->At[node->pt_typ[dir]][node->n[dir]];
}
/* get At[3] for pt_typ and n out of gridpoints */
void At3_pt_typ_n(int pt_typ[3], int n[3], tArray *At[3])
{
  int d;
  for(d=0; d<3; d++) At[d] = gridpoints->At[pt_typ[d]][n[d]];
}
/* get At[3] for n (assuming point type of node) out of gridpoints */
void At3_n(tNode *node, int n[3], tArray *At[3])
{
  int d;
  for(d=0; d<3; d++) At[d] = gridpoints->At[node->pt_typ[d]][n[d]];
}
/* get At[3] for a node out of gridpoints */
void node_At3(tNode *node, tArray *At[3])
{
  At3_n(node, node->n, At);
}

/* get St in dir for a node out of gridpoints */
tArray *node_St(tNode *node, int dir)
{
  return gridpoints->St[node->pt_typ[dir]][node->n[dir]];
}
/* get St[3] for pt_typ and n out of gridpoints */
void St3_pt_typ_n(int pt_typ[3], int n[3], tArray *St[3])
{
  int d;
  for(d=0; d<3; d++) St[d] = gridpoints->St[pt_typ[d]][n[d]];
}
/* get St[3] for n (assuming point type of node) out of gridpoints */
void St3_n(tNode *node, int n[3], tArray *St[3])
{
  int d;
  for(d=0; d<3; d++) St[d] = gridpoints->St[node->pt_typ[d]][n[d]];
}
/* get St[3] for a node out of gridpoints */
void node_St3(tNode *node, tArray *St[3])
{
  St3_n(node, node->n, St);
}

/* basis in of node in dir */
double node_basis(tNode *node, int dir,
                  int i, double x, int np)
{
  int typ = node->pt_typ[dir];
  return gridpoints->basis[typ](i, x, np);
}

/* check if all pt_typ[3] are are of type typ */
int pt_typ_all(int pt_typ[3], int typ)
{
  int d;
  int is_typ=1;
  for(d=0; d<3; d++) if(pt_typ[d]!=typ) { is_typ=0; break; }
  return is_typ;
}

/* check if pt_typ[3] has some of type typ */
int pt_typ_has(int pt_typ[3], int typ)
{
  int d;
  int has_typ=0;
  for(d=0; d<3; d++) if(pt_typ[d]==typ) { has_typ=1; break; }
  return has_typ;
}
