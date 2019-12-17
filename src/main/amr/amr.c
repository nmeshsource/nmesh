/* amr.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "amr.h"


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

  /* set amr globals */
  amr->MPIexchange = Par("amr_MPIexchange");
  amr->nghosts     = Par("amr_nghosts");

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

/*************************************************************************/
/* funcs for MPI exchange */
/*************************************************************************/

/* init exchange */
void MPIexchange_init_all_myln(tMesh *mesh)
{
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

/* set some local data */
void MPIexchange_set_all_myln_localdata(tMesh *mesh)
{
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

/* request exchange */
void MPIexchange_request_all_myln_data(tMesh *mesh)
{
  switch(Geti(amr->MPIexchange))
  {
  case 1:
    request_all_myln_surfaces_exchange(mesh);
    break;
  case 2:
    break;
  default:
    errorexit("unknown value in amr_MPIexchange");
  }
}

/* get buffers */
void MPIexchange_get_all_myln_data(tMesh *mesh)
{
  switch(Geti(amr->MPIexchange))
  {
  case 1:
    get_all_myln_surfaces(mesh);
    break;
  case 2:
    break;
  default:
    errorexit("unknown value in amr_MPIexchange");
  }
}

/* clean up after MPIexchange */
void MPIexchange_free_all_myln(tMesh *mesh)
{
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
  int N0,N1, I,J,K;

  N0 = n[0];
  N1 = n[1];
  I = i;
  J = j;
  K = k;

  switch(norm)
  {
  case 0:
    N0 = 1;
    I = 0;
    break;
  case 1:
    N1 = 1;
    J = 0;
    break;
  case 2:
    K = 0;
    break;
  default:
    errorexit("norm needs to be 0,1,2");
  }
  return I + N0*(J + N1*K);
}
