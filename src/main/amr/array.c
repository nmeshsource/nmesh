/* array.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"


#define PR 0

/* Temporary arrays in mm_array1/2.
   We can use arrays on the heap (1) or stack (0). */
//#define MM_TEMP_HEAP_ARRAYS 1
#define MM_TEMP_HEAP_ARRAYS 0


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

/************************************************************************/
/* macros and functions for matrices in raw C-arrays */
/************************************************************************/

/* Multiply two matricies A and B:
   AB = A B ,   AB_ij = A_il B_lj = At_li B_lj
   At contains a atn0 x atn1 matrix that is the transpose of A, stored in
   column major form. Ba contains a bn0 x bn1 matrix stored in column major
   form. ABa will contain AB stored in column major form. */
#define MM_SET_AB \
  for(int j=0; j<bn1; j++) \
    for(int i=0; i<atn1; i++) \
    { \
      double sum=0.0; \
      for(int l=0; l<atn0; l++) \
        sum += At[l + atn0*i] * B[l + atn0*j]; \
        /* At is col major transpose of A, and B is col major */ \
      AB[i + atn1*j] = sum; /* AB is col major */ \
    }
/* same as MM_SET_AB, but restrict sum over l to lrange,
   because At_li = 0 outside this range */
#define MM_SET_AB_IN_lrange \
  for(int j=0; j<bn1; j++) \
    for(int i=0; i<atn1; i++) \
    { \
      double sum=0.0; \
      int l0=lrange[0][i], l1=lrange[1][i]; \
      for(int l=l0; l<l1; l++) \
        sum += At[l + atn0*i] * B[l + atn0*j]; \
        /* At is col major transpose of A, and B is col major */ \
      AB[i + atn1*j] = sum; /* AB is col major */ \
    }

/* Set AB = A B ,   AB_ij = A_il B_lj = At_li B_lj
   Assume that A, B and AB point to different memory ararys */
void mm_Carray(int atn0, int atn1, int bn1,
               double *restrict At, double *restrict B,
               double *restrict AB)
{
  MM_SET_AB
}
/* Set AB = A B ,   AB_ij = A_il B_lj = At_li B_lj */
void mm_Carray_norestrict(int atn0, int atn1, int bn1,
                          double *At, double *B,
                          double *AB)
{
  MM_SET_AB
}

/* Set AB = A B ,   AB_ij = A_il B_lj = At_li B_lj
   but restrict sum over l to lrange, because At_li = 0 outside this range */
void mm_Carray_lrange(int atn0, int atn1, int bn1, int *lrange[2],
                      double *restrict At, double *restrict B,
                      double *restrict AB)
{
  MM_SET_AB_IN_lrange
}
/* Set AB = A B ,   AB_ij = A_il B_lj = At_li B_lj
   but restrict sum over l to lrange, because At_li = 0 outside this range */
void mm_Carray_lrange_norestrict(int atn0, int atn1, int bn1, int *lrange[2],
                                 double *At, double *B,
                                 double *AB)
{
  MM_SET_AB_IN_lrange
}


/************************************************************************/
/* functions for for arrays of type tArray* */
/************************************************************************/

/* multiply matrix in Ata with 3d array Ba in a direction,
   store result in ABa */
void mm_array_indir(tArray *Ata, tArray *Ba, int dir, tArray *ABa)
{
  if(PR)
  {
    PRF;printf(": dir=%d\n", dir);
    printf("Ata");printarray(Ata);
    printf("Ba");printarray(Ba);
  }

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
#define GET_At_A_AB_atn0_atn1_bn0_bn1 \
  double *At = Ata->d; \
  double *B  =  Ba->d; \
  double *AB = ABa->d; \
  int atn0 = Ata->n[0]; \
  int atn1 = Ata->n[1] * Ata->n[2]; \
  int bn0 = Ba->n[0]; \
  int bn1 = Ba->n[1] * Ba->n[2]; \
  if(atn0 != bn0) \
  { \
    printf("cannot multiply a %dx%d with %dx%d matrix:\n", atn1,atn0, bn0,bn1); \
    printf("Ata");printarray(Ata); \
    printf("Ba");printarray(Ba); \
    errorexit("Ata->n[0] != Ba->n[0]"); \
  }
/* Calculate AB = A B */
void mm_array0(tArray *Ata, tArray *Ba, tArray *ABa)
{
  GET_At_A_AB_atn0_atn1_bn0_bn1

  /* set AB */
  if(Ata->range[0] && Ata->range[1])
    mm_Carray_lrange(atn0,atn1, bn1, Ata->range, At, B, AB);
  else
    mm_Carray(atn0,atn1, bn1, At, B, AB);
}
/* same as mm_array0 but without restrict */
void mm_array0_norestrict(tArray *Ata, tArray *Ba, tArray *ABa)
{
  GET_At_A_AB_atn0_atn1_bn0_bn1

  /* set AB */
  if(Ata->range[0] && Ata->range[1])
    mm_Carray_lrange_norestrict(atn0,atn1, bn1, Ata->range, At, B, AB);
  else
    mm_Carray_norestrict(atn0,atn1, bn1, At, B, AB);
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
  int nat1 = Ata->n[1];
  int nAC[] = { nat1, nt1, n2 };
#if MM_TEMP_HEAP_ARRAYS == 1
  tArray *Ca = alloc_array(nC);
  tArray *ACa = alloc_array(nAC);
#else
  DECL_STACK_ARRAY(Ca, nC);
  DECL_STACK_ARRAY(ACa, nAC);
#endif
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

  if(PR)
  {
    PRFs(":\n");
    printf("Ba");printarray(Ba);
    printf("Ca");printarray(Ca);
  }

  /* now multiply */
  mm_array0(Ata, Ca, ACa);

  /* copy ACa into ABa */
  for(k=0; k<nAC[2]; k++)
    for(i=0; i<nAC[0]; i++)
      for(j=0; j<nAC[1]; j++)
        AB[Ind_n(j,i,k, nAB)] = AC[Ind_n(i,j,k, nAC)];

#if MM_TEMP_HEAP_ARRAYS == 1
  free_array(ACa);
  free_array(Ca);
#endif
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
  int nat1 = Ata->n[1];
  int nAC[] = { nat1, nt1, n0 };
#if MM_TEMP_HEAP_ARRAYS == 1
  tArray *Ca = alloc_array(nC);
  tArray *ACa = alloc_array(nAC);
#else
  DECL_STACK_ARRAY(Ca, nC);
  DECL_STACK_ARRAY(ACa, nAC);
#endif
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

  if(PR)
  {
    PRFs(":\n");
    printf("Ba");printarray(Ba);
    printf("Ca");printarray(Ca);
  }

  /* now multiply */
  mm_array0(Ata, Ca, ACa);

  /* copy ACa into ABa */
  for(i=0; i<nAC[0]; i++)
    for(j=0; j<nAC[1]; j++)
      for(k=0; k<nAC[2]; k++)
        AB[Ind_n(k,j,i, nAB)] = AC[Ind_n(i,j,k, nAC)];

#if MM_TEMP_HEAP_ARRAYS == 1
  free_array(ACa);
  free_array(Ca);
#endif
}

/* set entire array to value c */
void set_const_array(tArray *A, double c)
{
  int i;
  for(i=0; i<A->N; i++)  A->d[i] = c;
}

/* copy values from Src into Dest */
void copy_array_data(tArray *Src, tArray *Dest)
{
  double *dest = Dest->d;
  double *src  = Src->d;
  if(dest != src)
  {
    int sN = Src->N;
    int dN = Dest->N;
    int N = (sN < dN) ? sN : dN; // min
    int n = N * sizeof(Src->d[0]);
    memcpy(dest, src, n);
  }
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

/* copy np planes starting at pA from A into P starting at plane pP */
void copy_array_planes(int np, tArray *A, int dir, int pA,
                       tArray *P, int pP)
{
  int p;
  for(p=0; p<np; p++)  copy_array_plane(A, dir, pA+p, P, pP+p);
}

/* sum of array */
double sum_array(tArray *A, double p)
{
  int i;
  double sum = 0.;

  for(i=0; i<A->N; i++) sum += A->d[i];

  return sum;
}

/* sum of abs(array)^p */
double abs_p_sum_array(tArray *A, double p)
{
  int i;
  double sum = 0.;

  for(i=0; i<A->N; i++) sum += pow(fabs(A->d[i]), p);

  return sum;
}

/* norm of array */
double Lp_norm_array(tArray *A, double p)
{
  double sum = abs_p_sum_array(A, p);
  return pow(sum, 1./p);
}

/* D = A - B */
void array_diff(tArray *D, tArray *A, tArray *B)
{
  int i;
  for(i=0; i<D->N; i++) D->d[i] = A->d[i] - B->d[i];
}

/* D = (A - B)/B */
void array_reldiff(tArray *D, tArray *A, tArray *B)
{
  int i;
  for(i=0; i<D->N; i++)
  {
    double b = B->d[i];
    if(b!=0.) D->d[i] = (A->d[i] - b)/b;
    else      D->d[i] = (A->d[i] - b);
  }
}

/* norm of A-B */
double Lp_norm_array_diff(tArray *A, tArray *B, double p)
{
  tArray *D;
  double norm;

  if(!A) errorexit("A is NULL");

  D = alloc_array(A->n);
  array_diff(D, A, B);
  norm = Lp_norm_array(D, p);
  free_array(D);
  return norm;
}

/* norm of (A-B)/B */
double Lp_norm_array_reldiff(tArray *A, tArray *B, double p)
{
  tArray *D;
  double norm;

  if(!A) errorexit("A is NULL");

  D = alloc_array(A->n);
  array_reldiff(D, A, B);
  norm = Lp_norm_array(D, p);
  free_array(D);
  return norm;
}

/* max in an array */
double max_array(tArray *A, int *ind)
{
  return max_in_1d_array(A->d, A->N, ind);
}

/* min in an array */
double min_array(tArray *A, int *ind)
{
  return min_in_1d_array(A->d, A->N, ind);
}


/******************************************************************/
/* functions for array IO */
/******************************************************************/
/* write an array into a file */
int array_write(tMesh *mesh, tArray *array, char *fname)
{
  int *n = Arrn(array);
  int Nt = array->N + array->Ne;
  int nmemb = Nt * (array->ns);
  int size1 = sizeof(array->d[0]);
  FILE *fp;
  int IObufsz = ((mesh) ? Geti(Par("fwrite_bufsize")) : 0);
  char *IObuf; /* larger buffer for write */

  /* open destination file */
  fp = fopen_buf(fname, "wb", &IObuf,IObufsz);
  if(!fp) errorexits("failed opening %s", fname);

  /* header with some array info */
  fprintf(fp, "%d\n", n[0]);
  fprintf(fp, "%d\n", n[1]);
  fprintf(fp, "%d\n", n[2]);
  fprintf(fp, "%d\n", array->Ne);
  fprintf(fp, "%d\n", array->ns);
  /* write double data */
  fwrite_little(array->d, size1, nmemb, fp);

  fclose_buf(fp, &IObuf);
  return 0;
}

/* allocate and read an array form file*/
tArray *array_alloc_read(tMesh *mesh, char *fname, int fatal)
{
  tArray *array;
  int n[3];
  int Ne, ns, Nt, nmemb;
  int size1 = sizeof(array->d[0]);
  char str[1000];
  FILE *fp;
  int IObufsz = ((mesh) ? Geti(Par("fread_bufsize")) : 0);
  char *IObuf; /* larger buffer for write */

  /* open destination file */
  fp = fopen_buf(fname, "rb", &IObuf,IObufsz);
  if(!fp)
  {
    if(fatal) errorexits("failed opening %s", fname);
    else      return NULL;
  }

  /* read header with n, Ne, ns */
  fgets(str,999, fp);
  n[0] = atoi(str);
  fgets(str,999, fp);
  n[1] = atoi(str);
  fgets(str,999, fp);
  n[2] = atoi(str);
  fgets(str,999, fp);
  Ne = atoi(str);
  fgets(str,999, fp);
  ns = atoi(str);

  /* alloc array */
  array = alloc_array_with_segs(n, Ne, ns);

  /* read double data */
  Nt = array->N + array->Ne;
  nmemb = Nt * (array->ns);
  fread_little(array->d, size1, nmemb, fp);

  fclose_buf(fp, &IObuf);
  return 0;
}
