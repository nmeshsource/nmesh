/* nan_checker.c */
/* Wolfgang Tichy, 3/2019 */


#include "nmesh.h"


/* check if an array is finite */
int array_finite(tArray *a, char *name, int ijk[3])
{
  int ind;
  double *d = Arrd(a);
  int *n = Arrn(a);

  if(!d) return -1;

  //printf("%s  %p", name, a->d);
  ////printarray(a);
  //printf("\n");

  forarray(a, ind)
  {
    if(!isfinite(d[ind]))
    {
      int k = kOfInd_n(ind,n);
      int j = jOfInd_n_k(ind,n,k);
      int i = iOfInd_n_jk(ind,n,j,k);
      PRF;printf(":  %s->n[] = {%d,%d,%d}\n", name, n[0],n[1],n[2]);
      printf("d=%g at ind=%d i,j,k = %d %d %d\n", d[ind], ind, i,j,k);
      printf("%s", name);
      printarray(a);
      //errorexiti("not finite at %d", ind);
      ijk[0] = i;
      ijk[1] = j;
      ijk[2] = k;
      return ind;
    }
  }
  return -1;
}

/* check if a var is finite */
int array_in_var_finite(tNode *node, tArray *a, char *name, int ijk[3])
{
  int ind;

  /* first check array of var */
  ind = array_finite(a, name, ijk);
  if(ind>=0)
  {
    double Xb[3];
    
    if(node->np == a->N)
    {
      int ai = ijk[0];
      int aj = ijk[1];
      int ak = ijk[2];
      Xb[0] = node->Xb[0]->d[ai];
      Xb[1] = node->Xb[1]->d[aj];
      Xb[2] = node->Xb[2]->d[ak];
      printf("problem at Xb = %g %g %g\n", Xb[0], Xb[1], Xb[2]);
    }
    if(a->n[0] == 1)
    {
      int aj = ijk[1];
      int ak = ijk[2];
      Xb[1] = node->Xb[1]->d[aj];
      Xb[2] = node->Xb[2]->d[ak];
      printf("problem in plane0 at Xb = # %g %g\n", Xb[1], Xb[2]);
    }
    if(a->n[1] == 1)
    {
      int ai = ijk[0];
      int ak = ijk[2];
      Xb[0] = node->Xb[0]->d[ai];
      Xb[2] = node->Xb[2]->d[ak];
      printf("problem in plane1 at Xb = %g # %g\n", Xb[0], Xb[2]);
    }
    if(a->n[2] == 1)
    {
      int ai = ijk[0];
      int aj = ijk[1];
      Xb[0] = node->Xb[0]->d[ai];
      Xb[1] = node->Xb[1]->d[aj];
      printf("problem in plane2 at Xb = %g %g #\n", Xb[0], Xb[1]);
    }
  }

  return ind;
}


/* check if a var is finite */
int var_finite(tNode *node, int vi)
{
  tMesh *mesh = node->pat->mesh;
  char *vname = VarName(vi);
  char str[1000];
  char nname[100];
  int ind, ijk[3], f;
  tArray *a = VarA(node, vi);
  if(!a) return -1;

  /* set nname */
  nodename(node,nname,99);

  /* first check main array of var */
  snprintf(str,999, "%s: %s", nname, vname);
  ind = array_in_var_finite(node, a, str, ijk);
  if(ind>=0)
  {
    printf("node: %s,  var%d %s\n",
           nname, vi, vname);
    errorexiti("not finite at ind=%d", ind);
  }

  /* check surfaces */
  for(f=0; f<6; f++)
  {
    tDat *dat = node->dat;
    tSurface *sf = dat->s[f][vi];
    tArray *msa  = sf ? sf->mysurf : NULL;
    tArray *asa  = sf ? sf->ajsurf : NULL;
    tArray **nsa = sf ? sf->nbsurf : NULL;
    if(msa)
    {
      snprintf(str,999, "%s: %s:f%d:%s", nname, vname, f, "mysurf");
      ind = array_in_var_finite(node, msa, str, ijk);
      if(ind>=0)
      {
        printf("node: %s,  var%d %s, f%d mysurf\n",
               nname, vi, vname, f);
        errorexiti("not finite at ind=%d", ind);
      }
    }
    if(nsa)
    {
      int ni;
      for(ni=0; ni<sf->nnbsurf; ni++)
      {
        snprintf(str,999, "%s: %s:f%d:%s[%d]", nname, vname, f, "nbsurf", ni);
        ind = array_in_var_finite(node, nsa[ni], str, ijk);
        if(ind>=0)
        {
          printf("node: %s,  var%d %s, f%d nbsurf[%d]\n",
                 nname, vi, vname, f, ni);
          errorexiti("not finite at ind=%d", ind);
        }
      }
    }
    if(asa)
    {
      snprintf(str,999, "%s: %s:f%d:%s", nname, vname, f, "ajsurf");
      ind = array_in_var_finite(node, asa, str, ijk);
      if(ind>=0)
      {
        printf("node: %s,  var%d %s, f%d ajsurf\n",
               nname, vi, vname, f);
        errorexiti("not finite at ind=%d", ind);
      }
    }
  }

  return ind;
}
