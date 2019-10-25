/* tensors.c */
/* Wolfgang Tichy, 1/2019 & Bernd Bruegmann, 11/02 */

#include "nmesh.h"
#include <ctype.h>

#define ilistSTRLEN 8


/* decode tensor index string
   return list of indices into variable list
   return sign under reflections for symmetry boundaries
   note that we treat 3d indices i,j,k,... and 4d indices a,b,c,...
   and 2d indices q,r,s,...
   ilist[0], ilist[1], ... need to be freed by caller
   we could make a real parser to avoid the "if" for each case below */
void tensorindexlist(const char *tensind, int *nilist, char **ilist, int *sym)
{
  /* name of coordinates, could be made variable */
  char *coord[]  = {"x", "y", "z"};
  char *coord4[] = {"t", "x", "y", "z"};
  char *coord2[] = {"u", "v"}; /* for 2d tensors with 2 coords e.g. Y,Z */
  char *face[] = {"0","1","2","3","4","5"}; /* for tensors on node faces */
  int i, j, k, l, f;
  int n = 0;
  char *tensorindices = strdup(tensind);
  char *t;

  /* convert local copy to lower case since we ignore co/contra-variance */
  for(t = tensorindices; *t; t++)
    *t = tolower(*t);

  /* initialize symmetries */
  for(i = 0; i < 3*NINDEXLIST; i++)
    sym[i] = 1;

  /* now we treat each case */

  /* scalar */
  if(strcmp(tensorindices, "") == 0) {
    ilist[n] = calloc(ilistSTRLEN, sizeof(char));
    sprintf(ilist[n++], "%s", "");
  }

  /* 3d indices, i,j,k,... */
  if(strcmp(tensorindices, "i") == 0) {
    for(i = 0; i < 3; i++) {
      sym[3*n+i] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s", coord[i]);
    }
  }

  if(strcmp(tensorindices, "ij") == 0) {
    for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++) {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s", coord[i], coord[j]);
    }
  }

  if(strcmp(tensorindices, "ij+ji") == 0 ||
     strcmp(tensorindices, "(ij)" ) == 0) {
    for(i = 0; i < 3; i++)
    for(j = i; j < 3; j++) {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s", coord[i], coord[j]);
    }
  }

  if(strcmp(tensorindices, "ij-ji") == 0 ||
     strcmp(tensorindices, "[ij]") == 0) {
    for(i = 0; i < 3; i++)
    for(j = i+1; j < 3; j++) {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s", coord[i], coord[j]);
    }
  }

  if(strcmp(tensorindices, "ijk") == 0) {
    for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
    for(k = 0; k < 3; k++) {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }

  if(strcmp(tensorindices, "ijk+ikj") == 0 ||
     strcmp(tensorindices, "i(jk)") == 0 ) {
    for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
    for(k = j; k < 3; k++) {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }

  if(strcmp(tensorindices, "(ij)k") == 0 ||
     strcmp(tensorindices, "ijk+jik") == 0) {
    for(i = 0; i < 3; i++)
    for(j = i; j < 3; j++)
    for(k = 0; k < 3; k++) {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }

  if(strcmp(tensorindices, "ijk-jik") == 0 ||
      strcmp(tensorindices, "[ij]k") == 0) {
    for(i = 0; i < 3; i++)
    for(j = i+1; j < 3; j++)
    for(k = 0; k < 3; k++) {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }

  if(strcmp(tensorindices, "ijkl") == 0) {
    for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
    for(k = 0; k < 3; k++)
    for(l = 0; l < 3; l++)
    {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      sym[3*n+l] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s%s", coord[i], coord[j], coord[k], coord[l]);
    }
  }

  if(strcmp(tensorindices, "ijkl+ijlk+jikl+jilk") == 0 ||
     strcmp(tensorindices, "(ij)(kl)") == 0) {
    for(i = 0; i < 3; i++)
    for(j = i; j < 3; j++)
    for(k = 0; k < 3; k++)
    for(l = k; l < 3; l++)
    {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      sym[3*n+l] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s%s", coord[i], coord[j], coord[k], coord[l]);
    }
  }

  /* 4d indices, a,b,c,... */
  if(strcmp(tensorindices, "a") == 0) {
    for(i = 0; i <= 3; i++) {
      if(i > 0) sym[3*n+i-1] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s", coord4[i]);
    }
  }

  if(strcmp(tensorindices, "ab") == 0) {
    for(i = 0; i <= 3; i++)
    for(j = 0; j <= 3; j++) {
      if(i > 0) sym[3*n+i-1] *= -1;
      if(j > 0) sym[3*n+j-1] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s", coord4[i], coord4[j]);
    }
  }

  if(strcmp(tensorindices, "ai") == 0) {
    for(i = 0; i <= 3; i++)
    for(j = 0; j <  3; j++) {
      if(i > 0) sym[3*n+i-1] *= -1;
      sym[3*n+j] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s", coord4[i], coord[j]);
    }
  }

  if(strcmp(tensorindices, "ab+ba") == 0 ||
     strcmp(tensorindices, "(ab)") == 0) {
    for(i = 0; i <= 3; i++)
    for(j = i; j <= 3; j++) {
      if(i > 0) sym[3*n+i-1] *= -1;
      if(j > 0) sym[3*n+j-1] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s", coord4[i], coord4[j]);
    }
  }

  if(strcmp(tensorindices, "abc+acb") == 0 ||
     strcmp(tensorindices, "a(bc)") == 0) {
    for(i = 0; i <= 3; i++)
    for(j = 0; j <= 3; j++)
    for(k = j; k <= 3; k++) {
      if(i > 0) sym[3*n+i-1] *= -1;
      if(j > 0) sym[3*n+j-1] *= -1;
      if(k > 0) sym[3*n+k-1] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord4[i], coord4[j], coord4[k]);
    }
  }

  if(strcmp(tensorindices, "aij+aji") == 0 ||
     strcmp(tensorindices, "a(ij)") == 0) {
    for(i = 0; i <= 3; i++)
    for(j = 0; j <  3; j++)
    for(k = j; k <  3; k++) {
      if(i > 0) sym[3*n+i-1] *= -1;
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord4[i], coord[j], coord[k]);
    }
  }

  if(strcmp(tensorindices, "abc+bac") == 0 ||
     strcmp(tensorindices, "(ab)c")   == 0) {
    for(i = 0; i <= 3; i++)
    for(j = i; j <= 3; j++)
    for(k = 0; k <= 3; k++) {
      if(i > 0) sym[3*n+i-1] *= -1;
      if(j > 0) sym[3*n+j-1] *= -1;
      if(k > 0) sym[3*n+k-1] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord4[i], coord4[j], coord4[k]);
    }
  }

  if(strcmp(tensorindices, "abi+bai")   == 0 ||
     strcmp(tensorindices, "(ab)i")   == 0) {
    for(i = 0; i <= 3; i++)
    for(j = i; j <= 3; j++)
    for(k = 0; k <  3; k++) {
      if(i > 0) sym[3*n+i-1] *= -1;
      if(j > 0) sym[3*n+j-1] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord4[i], coord4[j], coord[k]);
    }
  }

  if(strcmp(tensorindices, "iab+iba")   == 0 ||
     strcmp(tensorindices, "i(ab)")   == 0) {
    for(i = 0; i <  3; i++)
    for(j = 0; j <= 3; j++)
    for(k = j; k <= 3; k++) {
      sym[3*n+i] *= -1;
      if(j > 0) sym[3*n+j-1] *= -1;
      if(k > 0) sym[3*n+k-1] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord[i], coord4[j], coord4[k]);
    }
  }

  if(strcmp(tensorindices, "abij+abji+baij+baji") == 0 ||
     strcmp(tensorindices, "(ab)(ij)") == 0) {
    for(i = 0; i <= 3; i++)
    for(j = i; j <= 3; j++)
    for(k = 0; k < 3; k++)
    for(l = k; l < 3; l++)
    {
      if(i > 0) sym[3*n+i-1] *= -1;
      if(j > 0) sym[3*n+j-1] *= -1;
      sym[3*n+k] *= -1;
      sym[3*n+l] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s%s", coord4[i], coord4[j], coord[k], coord[l]);
    }
  }

  if(strcmp(tensorindices, "ijab+ijba+jiab+jiba") == 0 ||
     strcmp(tensorindices, "(ij)(ab)") == 0) {
    for(i = 0; i < 3; i++)
    for(j = i; j < 3; j++)
    for(k = 0; k <= 3; k++)
    for(l = k; l <= 3; l++)
    {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      if(k > 0) sym[3*n+k-1] *= -1;
      if(l > 0) sym[3*n+l-1] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s%s", coord[i], coord[j], coord4[k], coord4[l]);
    }
  }

  /* 2d indices q,r,s,... \in [1,2] */
  if(strcmp(tensorindices, "q") == 0) {
    for(i = 0; i < 2; i++) {
      sym[3*n+i] *= -1;   /* FIXME: these syms probably need to be changed */
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s", coord2[i]);
    }
  }

  if(strcmp(tensorindices, "qr") == 0) {
    for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++) {
      sym[3*n+i] *= -1;   /* FIXME: these syms probably need to be changed */
      sym[3*n+j] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s", coord2[i], coord2[j]);
    }
  }

  if(strcmp(tensorindices, "qr+rq") == 0 ||
     strcmp(tensorindices, "(qr)" ) == 0) {
    for(i = 0; i < 2; i++)
    for(j = i; j < 2; j++) {
      sym[3*n+i] *= -1;   /* FIXME: these syms probably need to be changed */
      sym[3*n+j] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s", coord2[i], coord2[j]);
    }
  }

  if(strcmp(tensorindices, "qrs+rqs") == 0 ||
     strcmp(tensorindices, "(qr)s" ) == 0) {
    for(i = 0; i < 2; i++)
    for(j = i; j < 2; j++)
    for(k = 0; k < 2; k++) {
      sym[3*n+i] *= -1;   /* FIXME: these syms probably need to be changed */
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord2[i], coord2[j], coord2[k]);
    }
  }

  if(strcmp(tensorindices, "qrst+qrts+rqst+rqts") == 0 ||
     strcmp(tensorindices, "(qr)(st)") == 0) {
    for(i = 0; i < 2; i++)
    for(j = i; j < 2; j++)
    for(k = 0; k < 2; k++)
    for(l = k; l < 2; l++)
    {
      sym[3*n+i] *= -1;   /* FIXME: these syms probably need to be changed */
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      sym[3*n+l] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s%s",
              coord2[i], coord2[j], coord2[k], coord2[l]);
    }
  }

  if(strcmp(tensorindices, "iq") == 0) {
    for(i = 0; i < 3; i++)
    for(j = 0; j < 2; j++) {
      sym[3*n+i] *= -1;  /* FIXME: these syms probably need to be changed */
      sym[3*n+j] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s", coord[i], coord2[j]);
    }
  }

  if(strcmp(tensorindices, "iqr") == 0) {
    for(i = 0; i < 3; i++)
    for(j = 0; j < 2; j++)
    for(k = 0; k < 2; k++) {
      sym[3*n+i] *= -1;  /* FIXME: these syms probably need to be changed */
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord[i], coord2[j], coord2[k]);
    }
  }

  if(strcmp(tensorindices, "iqr+irq") == 0 ||
     strcmp(tensorindices, "i(qr)") == 0 ) {
    for(i = 0; i < 3; i++)
    for(j = 0; j < 2; j++)
    for(k = j; k < 2; k++) {
      sym[3*n+i] *= -1;  /* FIXME: these syms probably need to be changed */
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", coord[i], coord2[j], coord2[k]);
    }
  }

  /* node surfaces (denoted by @) plus 3d indices i,j,... */
  if(strcmp(tensorindices, "@") == 0) {
    for(f=0; f<6; f++)
    {
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s", face[f]);
    }
  }

  if(strcmp(tensorindices, "@i") == 0) {
    for(f=0; f<6; f++)
    for(i = 0; i < 3; i++) {
      sym[3*n+i] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s", face[f], coord[i]);
    }
  }

  if(strcmp(tensorindices, "@ij") == 0) {
    for(f=0; f<6; f++)
    for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++) {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", face[f], coord[i], coord[j]);
    }
  }

  if(strcmp(tensorindices, "@ij+@ji") == 0 ||
     strcmp(tensorindices, "@(ij)" ) == 0) {
    for(f=0; f<6; f++)
    for(i = 0; i < 3; i++)
    for(j = i; j < 3; j++) {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", face[f], coord[i], coord[j]);
    }
  }

  if(strcmp(tensorindices, "@ij-@ji") == 0 ||
     strcmp(tensorindices, "@[ij]") == 0) {
    for(f=0; f<6; f++)
    for(i = 0; i < 3; i++)
    for(j = i+1; j < 3; j++) {
      sym[3*n+i] *= -1;
      sym[3*n+j] *= -1;
      ilist[n] = calloc(ilistSTRLEN, sizeof(char));
      sprintf(ilist[n++], "%s%s%s", face[f], coord[i], coord[j]);
    }
  }


  /* error */
  if(n == 0) {
    printf("Unknown tensor index string %s.\n", tensorindices);
    printf("Implemented (besides the empty string) are:\n");
    printf("i, ij, (ij), [ij], ijk, i(jk), (ij)k [ij]k ijkl (ij)(kl)\n");
    printf("ij+ji, ij-ji, ijk+ikj, ijk+jik\n");
    printf("a, ab, ab+ba, abc+acb\n");
    printf("(ab), a(bc), (ab)c, (ab)i, (ab)(ij)\n");
    printf("q, qr, qr+rq, qrs+rqs, qrst+qrts+rqst+rqts\n");
    printf("q, qr, (qr), (qr)s, (qr)(st)\n");
    printf("iq, iqr, i(qr)\n");
    printf("@i, @ij, @(ij), @[ij]\n");
    printf("Other combinations can be added to main/tensors.c .\n");
    errorexits("Error in tensor index string %s", tensorindices);
  }

  *nilist = n;
  free(tensorindices);
}
