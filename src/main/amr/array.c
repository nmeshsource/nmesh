/* array.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"


/* NOTE about nmesh variables and arrays:
   var arrays a = v_{ijk} are indexed like this:
   ind = i + n[0]*(j + n[1]*k)
   So rows in X-dir can be accessed sequentially as a[i]. I.e. v is 
   stored like J major (~column major), where J = j + n[1]*k.
   Thus we can simply store the diff matrices and such in row-major form.
   The matrix multiplication simplifies */


/* Multiply two matricies A and B:  AB = A B ,   AB_ij = A_il B_lj
   Aa contains an (Aa->n[0]) x (Aa->n[1]) matrix stored in row-major
   form. We should have Aa->n[2] = 1.
   Ba contains a (Ba->n[0]) x (Ba->n[1] * Ba->n[2]) matrix stored in
   column major form. ABa will contain AB stored in column major form. */
void mm_RMarray_times_CMarray0(tArray *Aa, tArray *Ba, tArray *ABa)
{
  double *restrict A  =  Aa->a;
  double *restrict B  =  Ba->a;
  double *restrict AB = ABa->a;
  int an0 = Aa->n[0];
  int an1 = Aa->n[1];
  int bn0 = Ba->n[0];
  int bn1 = Ba->n[1] * Ba->n[2];
  int i,l,j;

  if(an1 != bn0) errorexit("Aa->n[1] != Ba->n[0]");

  for(i=0; i<an0; i++)
    for(j=0; j<bn1; j++)
    {
      double sum=0.0;
      for(l=0; l<an1; l++)
        sum += A[an1*i + l] * B[l + an1*j];
        // A is row major and B is col major
      AB[i + an0*j] = sum; // AB is col major
    }
}

/* Multiply two matricies A and B:  AB = A B ,   AB_ij = A_il B_lj
   Aa contains an (Aa->n[0]) x (Aa->n[1]) matrix stored in row-major
   form. We should have Aa->n[2] = 1.
   Ba will first be copied into a form Ta where j is the least major index.
   Next we call mm_RMarray_times_CMarray0 with Ta to get AT.
   Then we construct AB from AT */
void mm_RMarray_times_CMarray1(tArray *Aa, tArray *Ba, tArray *ABa)
{
  int n0 = Ba->n[0];
  int n1 = Ba->n[1];
  int n2 = Ba->n[2];
  int nt0 = n1;
  int nt1 = n0;
  int nT[] = { nt0, nt1, n2 };
  tArray *Ta = alloc_array(nT);
  int nat0 = Aa->n[0];
  int nat1 = nt1;
  int nAT[] = { nat0, nat1, n2 };
  tArray *ATa = alloc_array(nAT);
  int nB[] = { n0,n1,n2 };
  int nAB[] = { ABa->n[0],ABa->n[1],ABa->n[2] };
  double *restrict B = Ba->a;
  double *restrict T = Ta->a;
  double *restrict AB = ABa->a;
  double *restrict AT = ATa->a;
  int i,j,k;

  /* copy Ba into Ta */
  for(k=0; k<n2; k++)
    for(i=0; i<n0; i++)
      for(j=0; j<n1; j++)
        T[Ind_n(j,i,k, nT)] = B[Ind_n(i,j,k, nB)];

  /* now multiply */
  mm_RMarray_times_CMarray0(Aa, Ta, ATa);

  /* copy ATa into ABa */
  for(k=0; k<nAB[2]; k++)
    for(j=0; j<nAB[1]; j++)
      for(i=0; i<nAB[0]; i++)
        AB[Ind_n(i,j,k, nAB)] = AT[Ind_n(j,i,k, nAT)];

  free_array(ATa);
  free_array(Ta);
}

/* Multiply two matricies A and B:  AB = A B ,   AB_ij = A_il B_lj
   Aa contains an (Aa->n[0]) x (Aa->n[1]) matrix stored in row-major
   form. We should have Aa->n[2] = 1.
   Ba will first be copied into a form Ta where k is the least major index.
   Next we call mm_RMarray_times_CMarray0 with Ta to get AT.
   Then we construct AB from AT */
void mm_RMarray_times_CMarray2(tArray *Aa, tArray *Ba, tArray *ABa)
{
  int n0 = Ba->n[0];
  int n1 = Ba->n[1];
  int n2 = Ba->n[2];
  int nt0 = n2;
  int nt1 = n1;
  int nT[] = { nt0, nt1, n0 };
  tArray *Ta = alloc_array(nT);
  int nat0 = Aa->n[0];
  int nat1 = nt1;
  int nAT[] = { nat0, nat1, n0 };
  tArray *ATa = alloc_array(nAT);
  int nB[] = { n0,n1,n2 };
  int nAB[] = { ABa->n[0],ABa->n[1],ABa->n[2] };
  double *restrict B = Ba->a;
  double *restrict T = Ta->a;
  double *restrict AB = ABa->a;
  double *restrict AT = ATa->a;
  int i,j,k;

  /* copy Ba into Ta */
  for(i=0; i<n0; i++)
    for(j=0; j<n1; j++)
      for(k=0; k<n2; k++)
        T[Ind_n(k,j,i, nT)] = B[Ind_n(i,j,k, nB)];

  /* now multiply */
  mm_RMarray_times_CMarray0(Aa, Ta, ATa);

  /* copy ATa into ABa */
  for(k=0; k<nAB[2]; k++)
    for(j=0; j<nAB[1]; j++)
      for(i=0; i<nAB[0]; i++)
        AB[Ind_n(i,j,k, nAB)] = AT[Ind_n(k,j,i, nAT)];

  free_array(ATa);
  free_array(Ta);
}
