/* array.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"


/* NOTE about nmesh variables and arrays:
   var arrays a[ind] = v_{ijk} are indexed like this:
   ind = i + n[0]*(j + n[1]*k)
   So rows in X-dir can be accessed sequentially as a[i + n[0]*J]. I.e. v is
   stored in J major (~column major) form, where J = j + n[1]*k.
   So we will store all other matrices also in column major form. But to speed
   up matrix products A B, we will store always the transpose of A, since
   computing (A^t)^t B is faster than A B.
   Thus we simply store the transposes of diff matrices and such to speed up
   matrix multiplication!!! */

/* multiply matrix in Ata with 3d array Ba in a direction,
   store result in ABa */
void mm_array_indir(tArray *Ata, tArray *Ba, int dir, tArray *ABa)
{
  switch(dir)
  {
  case 0:
    mm_array0(Ata, Ba, ABa);  break;
  case 1:
    mm_array1(Ata, Ba, ABa);  break;
  case 2:
    mm_array2(Ata, Ba, ABa);  break;
  default:
    errorexit("dir must be 0,1,2");
  }
}

/* Multiply two matricies A and B:
   AB = A B ,   AB_ij = A_il B_lj = At_li B_lj
   Ata contains an (Ata->n[0]) x (Ata->n[1]) matrix that is the transpose
   of A, stored in column major form. We should have Ata->n[2] = 1.
   Ba contains a (Ba->n[0]) x (Ba->n[1] * Ba->n[2]) matrix stored in
   column major form. ABa will contain AB stored in column major form. */
void mm_array0(tArray *Ata, tArray *Ba, tArray *ABa)
{
  double *restrict At = Ata->d;
  double *restrict B  =  Ba->d;
  double *restrict AB = ABa->d;
  int atn0 = Ata->n[0];
  int atn1 = Ata->n[1] * Ata->n[2];
  int bn0 = Ba->n[0];
  int bn1 = Ba->n[1] * Ba->n[2];
  int i,l,j;

  if(atn0 != bn0) errorexit("Ata->n[0] != Ba->n[0]");

  for(j=0; j<bn1; j++)
    for(i=0; i<atn1; i++)
    {
      double sum=0.0;
      for(l=0; l<atn0; l++)
        sum += At[l + atn0*i] * B[l + atn0*j];
        // At is col major transpose of A, and B is col major
      AB[i + atn1*j] = sum; // AB is col major
    }
}

/* Multiply two matricies A and B:  AB = A B ,   AB_ij = A_il B_lj
   Ata contains an (Ata->n[0]) x (Ata->n[1]) matrix stored in row-major
   form. We should have Ata->n[2] = 1.
   Ba will first be copied into a form Ca where j is the least major index.
   Next we call mm_array0 with Ca to get AC.
   Then we construct AB from AC */
void mm_array1(tArray *Ata, tArray *Ba, tArray *ABa)
{
  int n0 = Ba->n[0];
  int n1 = Ba->n[1];
  int n2 = Ba->n[2];
  int nt0 = n1;
  int nt1 = n0;
  int nC[] = { nt0, nt1, n2 };
  tArray *Ca = alloc_array(nC);
  int nat1 = Ata->n[1];
  int nAC[] = { nat1, nt1, n2 };
  tArray *ACa = alloc_array(nAC);
  int nB[] = { n0,n1,n2 };
  int nAB[] = { nt1, nat1, n2 };
  double *B = Ba->d;
  double *AB = ABa->d;
  double *restrict C = Ca->d;
  double *restrict AC = ACa->d;
  int i,j,k;

  /* copy Ba into Ca */
  for(k=0; k<n2; k++)
    for(i=0; i<n0; i++)
      for(j=0; j<n1; j++)
        C[Ind_n(j,i,k, nC)] = B[Ind_n(i,j,k, nB)];

  /* now multiply */
  mm_array0(Ata, Ca, ACa);

  /* copy ACa into ABa */
  for(k=0; k<nAC[2]; k++)
    for(i=0; i<nAC[0]; i++)
      for(j=0; j<nAC[1]; j++)
        AB[Ind_n(j,i,k, nAB)] = AC[Ind_n(i,j,k, nAC)];

  free_array(ACa);
  free_array(Ca);
}

/* Multiply two matricies A and B:  AB = A B ,   AB_ij = A_il B_lj
   Ata contains an (Ata->n[0]) x (Ata->n[1]) matrix stored in row-major
   form. We should have Ata->n[2] = 1.
   Ba will first be copied into a form Ca where k is the least major index.
   Next we call mm_array0 with Ca to get AC.
   Then we construct AB from AC */
void mm_array2(tArray *Ata, tArray *Ba, tArray *ABa)
{
  int n0 = Ba->n[0];
  int n1 = Ba->n[1];
  int n2 = Ba->n[2];
  int nt0 = n2;
  int nt1 = n1;
  int nC[] = { nt0, nt1, n0 };
  tArray *Ca = alloc_array(nC);
  int nat1 = Ata->n[1];
  int nAC[] = { nat1, nt1, n0 };
  tArray *ACa = alloc_array(nAC);
  int nB[] = { n0,n1,n2 };
  int nAB[] = { n0, nt1, nat1 };
  double *B = Ba->d;
  double *AB = ABa->d;
  double *restrict C = Ca->d;
  double *restrict AC = ACa->d;
  int i,j,k;

  /* copy Ba into Ca */
  for(i=0; i<n0; i++)
    for(j=0; j<n1; j++)
      for(k=0; k<n2; k++)
        C[Ind_n(k,j,i, nC)] = B[Ind_n(i,j,k, nB)];

  /* now multiply */
  mm_array0(Ata, Ca, ACa);

  /* copy ACa into ABa */
  for(i=0; i<nAC[0]; i++)
    for(j=0; j<nAC[1]; j++)
      for(k=0; k<nAC[2]; k++)
        AB[Ind_n(k,j,i, nAB)] = AC[Ind_n(i,j,k, nAC)];

  free_array(ACa);
  free_array(Ca);
}

/* set entire array to value c */
void set_const_array(tArray *A, double c)
{
  int i;
  for(i=0; i<A->N; i++)  A->d[i] = c;
}

/* take plane pA with normal dir from A and copy it into P to plane pP */
void copy_array_plane(tArray *A, int dir, int pA, tArray *P, int pP)
{
  int i,j,k;
  switch(dir)
  {
  case 0:
    forplane0(i,j,k, A->n, pA)
      P->d[Ind_n(pP,j,k, P->n)] = A->d[Ind_n(i,j,k, A->n)];
    break;
  case 1:
    forplane1(i,j,k, A->n, pA)
      P->d[Ind_n(i,pP,k, P->n)] = A->d[Ind_n(i,j,k, A->n)];
    break;
  case 2:
    forplane2(i,j,k, A->n, pA)
      P->d[Ind_n(i,j,pP, P->n)] = A->d[Ind_n(i,j,k, A->n)];
    break;
  default:
    errorexit("dir must be 0,1,2");
  }
}

/* norm of array */
double Lp_norm_array(tArray *A, double p)
{
  int i;
  double sum = 0.;

  for(i=0; i<A->N; i++) sum += pow(fabs(A->d[i]), p);

  return pow(sum, 1./p);
}

/* D = A - B */
void array_diff(tArray *D, tArray *A, tArray *B)
{
  int i;
  for(i=0; i<D->N; i++) D->d[i] = A->d[i] - B->d[i];
}

/* norm of A-B */
double Lp_norm_array_diff(tArray *A, tArray *B, double p)
{
  tArray *D = alloc_array(A->n);
  double norm;

  array_diff(D, A, B);
  norm = Lp_norm_array(D, p);
  free_array(D);
  return norm;
}
