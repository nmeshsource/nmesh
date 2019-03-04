/* amr.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "amr.h"



/*************************************************************************/
/* various untilities */
/*************************************************************************/

/* index, but reduced by dimension along norm */
int Ind_n_norm(int i, int j, int k, int n[3], int norm)
{
  int N0,N1, I,J,K;
  switch(norm)
  {
  case 0:
    N0 = 1;
    N1 = n[1];
    I = 0;
    J = j;
    K = k;
    break;
  case 1:
    N0 = n[0];
    N1 = 1;
    I = i;
    J = 0;
    K = k;
    break;
  case 2:
    N0 = n[0];
    N1 = n[1];
    I = i;
    J = j;
    K = 0;
    break;
  default:
    errorexit("norm needs to be 0,1,2");
  }
  return I + N0*(J + N1*K);
}

