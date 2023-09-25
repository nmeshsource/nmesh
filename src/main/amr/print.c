/* print.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"

/* mesh made in main */
extern tMesh *main_mesh;


/* print information about various objects */

void printmesh(tMesh *mesh)
{
  int size = nMPI_size();
  int p;

  if(mesh==main_mesh) printf("mesh=main_mesh: ");
  else                printf("mesh=%p: ", (void *) mesh);

  printf("npats=%d npdb=%d nvdb=%d ", mesh->npats, mesh->npdb, mesh->nvdb);
  printf("nln=%lu ", mesh->eidlim[size-1]);
  printf("dt=%g\n", mesh->dt);

  forpatches(mesh, p)
    printpatch(mesh->pat[p]);
  printf("leaf nodes on this rank: %lu \n", mesh->nmyelm);
  printmyelms(mesh);
  printf("leaf nodes on nearby ranks: %lu\n", mesh->nnbelm);
  printnbelms(mesh);
}

void printpatch(tPat *pat)
{
  printf("p%d: [%g,%g]x[%g,%g]x[%g,%g]\n",
         pat->p, pat->bbox[0], pat->bbox[1], pat->bbox[2], pat->bbox[3],
         pat->bbox[4],pat->bbox[5]);
}

void printCI(tPat *pat)
{
  tCoordInfo *CI = pat->CI;
  int i, j, p;

  printf("p%d: CI->dom=%d CI->type=%d\n",
         pat->p, CI->dom, CI->type);

  /* s, xc */
  printf("CI->s[] = [ ");
  for(i=0;i<6;i++) printf("%g ", CI->s[i]);
  printf("]\n");
  printf("CI->xc[] = [ ");
  for(i=0;i<3;i++) printf("%g ", CI->xc[i]);
  printf("]\n");

  /* iSurf */
  for(i=p=0; i<6; i++) if(CI->iSurf[i]) p++;
  if(p)
  {
    //printf("CI->iFS[] = [ ");
    //for(i=0;i<6;i++) printf("%d ", CI->iFS[i]);
    //printf("]\n");
    printf("CI->iSurf[] = [ ");
    for(i=0;i<6;i++) printf("%d ", CI->iSurf[i]);
    printf("]\n");
    printf("idSurfdX[][] = [");
    for(i=0;i<6;i++)
    {
      for(j=1;j<3;j++) printf(" %d", CI->idSurfdX[i][j]);
      printf(";");
    }
    printf("]\n");
  }

  /* FSurf */
  printf("CI->FSurf[] = [ ");
  for(i=0;i<6;i++) printf("%p ", (void *) (size_t) CI->FSurf[i]);
  printf("]\n");
  for(i=p=0; i<6; i++) if(CI->dFSurfdC[i]) p++;
  if(p)
  {
    printf("dFSurfdC[] = [");
    for(i=0;i<6;i++) printf(" %p", (void *) (size_t) CI->dFSurfdC[i]);
    printf("]\n");
  }
}

void printeloc(const tEloc *eloc)
{
  char ch;
  int i;
  //int len = strnlen(eloc->loc, NLOCS);
  printf("%d_", eloc->p);
  for(i=0; i<eloc->l; i++)
  {
    ch = eloc->loc[i];
    if(ch) printf("%c", ch);
    else   printf(".");
  }
  //printf("|");
  //if(eloc->l > len) printf(">");
  //if(eloc->l < len) printf("<");
  #if 0
  for(i=eloc->l; i<len; i++)
  {
    ch = eloc->loc[i];
    if(ch) printf("%c", ch);
    else   printf(".");
  }
  #endif
  //printf(" eid%lu", eloc->eid);
}

void printeloc_s(const tEloc *eloc, const char *s)
{
  printeloc(eloc);
  printf("%s", s);
}

void printeploc(const tEploc *eploc)
{
  tEloc eloc[1];
  eloc_from_eploc(eloc, eploc);
  printeloc(eloc);
}
void printeploc_s(const tEploc *eploc, const char *s)
{
  tEloc eloc[1];
  eloc_from_eploc(eloc, eploc);
  printeloc_s(eloc, s);
}

void printelm0(const tElm0 *e, const char *s)
{
  printeploc(e->eploc);
  printf(": eid%lu [%g,%g]x[%g,%g]x[%g,%g]",
         e->eploc->eid,
         e->bbox[0],e->bbox[1], e->bbox[2],e->bbox[3], e->bbox[4],e->bbox[5]);
  printf(" n=%dx%dx%d=%d", e->n[0],e->n[1],e->n[2], e->np);
  //printf(" rflag=%d", e->rflag);
  printf(" datrank=%d", e->datrank);
  printf("%s", s);
}

void printelm(const tElm *e)
{
  tDat *dat = e->dat;
  int nnbinfo_neg = 0;
  int nfnb_tot    = 0;
  int i;
  union { const tElm *elm; tElm0 *elm0; } e2e0;
  e2e0.elm = e;
  printelm0(e2e0.elm0, "");
  printf(" dat:%s\n", dat ? "yes" : "no");
  printf(" fnb =");
  for(i=0; i<6; i++)
  {
    int j;
    printf(" {");
    if(dat && dat->info->nnbinfo[i]<0)
    {
      printf("? ");
      nnbinfo_neg=1;
    }
    //else
    {
      for(j=0; j<e->nfnb[i]; j++)
      {
        if(e->fnb[i][j]) printeploc(e->fnb[i][j]->eploc);
        else             printf("nil");
        if(j<e->nfnb[i]-1) printf(" ");
        nfnb_tot++;
      }
    }
    printf("}");
  }
  printf("\n");
  if(nnbinfo_neg || nfnb_tot==0)
  {
    printf(" nbi =");
    for(i=0; i<6; i++)
    {
      printf(" {");
      print_amr_elm_nbinfo(e, i);
      printf("}");
    }
    printf("\n");
  }
}

void printelmarray(long nelms, tElm **elm)
{
  long i;
  printf("%ld elms: ", nelms);
  for(i=0; i<nelms; i++)
    printeploc_s(elm[i]->eploc, " ");
  printf("\n");
}

void printelmlist(struct list_head *elm_head)
{
  struct list_head *pos;
  printf("%ld elms: ", list_count_nodes(elm_head));
  list_for_each(pos, elm_head)
  {
    tElm *elm = list_entry(pos, tElm, list);
    printeploc_s(elm->eploc, " ");
  }
  printf("\n");
}

void printelmglist(struct list_head *elm_head)
{
  struct list_head *pos;
  printf("%ld elms: ", list_count_nodes(elm_head));
  list_for_each(pos, elm_head)
  {
    tElm *elm = glist_entry(pos);
    printeploc_s(elm->eploc, " ");
  }
  printf("\n");
}

void printelmsinglist(struct list_head *elm_head)
{
  struct list_head *pos;
  printf("%ld elms:\n", list_count_nodes(elm_head));
  list_for_each(pos, elm_head)
  {
    tElm *elm = glist_entry(pos);
    printelm(elm);
  }
}

void printmyelms(tMesh *mesh)
{
  formyelms_noomp(mesh)
  {
    tElm *elm = MyElm;
    printelm(elm);
  }
}

void printnbelms(tMesh *mesh)
{
  int ei;
  for(ei=0; ei < mesh->nnbelm; ei++)
  {
    tElm *elm = mesh->nbelm[ei];
    printelm(elm);
  }
}

void printeidlim(tMesh *mesh)
{
  int size = nMPI_size();
  ulong *eidlim = mesh->eidlim;
  int rk;
  for(rk=0; rk < size; rk++)
  {
    printf("eidlim[%d]=%lu\n", rk, eidlim[rk]);
  }
}

void print_amr_elm_nbinfo(const tElm *elm, int face)
{
  tMesh *mesh = elm->pat->mesh;
  int i_nbinfo = Ind("amr_elm_nbinfo0") + face;
  tArray *nbinfo = VarA(elm, i_nbinfo);
  tDat *dat = elm->dat;
  if(dat && dat->info->nnbinfo[face]<0) printf("- ");
  if(!nbinfo) printf("nil");
  printarray_eploc(nbinfo, 0);
}

/* print nodename */
void pr_nodename(tNode *node)
{
  char str[100];
  printf("%s", nodename(node, str,99));
}

/* print info inside dat of node */
void printnodeinfo(const tElm *elm)
{
  tDat *dat = elm->dat;
  if(dat)
  {
    tNodeInfo *info = dat->info;
    printf("info->  ");
    printf("evo_troubled=%d  trbl_score=%d  use_fv=%d\n",
           info->evo_troubled, info->trbl_score, info->use_fv);
    printf("info->  trbl_ref: ");
    printref(info->trbl_ref);
  }
  else
  {
    printf("dat=NULL\n");
  }
}

void printelm_nodeinfo(const tElm *elm)
{
  union { const tElm *elm; tElm0 *elm0; } e2e0;
  e2e0.elm = elm;
  printelm0(e2e0.elm0, "\n");
  printnodeinfo(elm);
}

/* print a variable in a node */
void printvar_innode(tNode *node, int vi)
{
  int f;
  tMesh *mesh = node->pat->mesh;
  char *name = VarName(vi);
  int zones =  MeshVarSurfacezones(mesh, vi);
  int type =  MeshVarType(mesh, vi);
  tDat *dat = node->dat;
  tArray *va = dat ? dat->v[vi] : NULL;
  char s[100];

  printf("%s %s Ind=%d type=%d zones=%d Array",
         nodename(node,s,99), name, vi, type, zones);
  printarray(va);

  if(va)
    for(f=0; f<6; f++)
    {
      tSurface *sf = dat->s[f][vi];
      tArray *msa  = sf ? sf->mysurf : NULL;
      tArray *asa  = sf ? sf->ajsurf : NULL;
      tArray **nsa = sf ? sf->nbsurf : NULL;
      if(msa)
      {
        printf("f%d allocd=%d mysurf", f, sf->allocd_mysurf);
        printarray(msa);
      }
      if(asa)
      {
        printf("f%d allocd=%d ajsurf", f, sf->allocd_ajsurf);
        printarray(asa);
      }
      if(nsa)
      {
        int ni;
        for(ni=0; ni<sf->nnbsurf; ni++)
        {
          printf("f%d fnb[f][%d]=%s allocd=%d nbsurf[%d]",
                 f, ni, nodename(node->fnb[f][ni],s,99),
                 sf->allocd_nbsurf[ni], ni);
          printarray(nsa[ni]);
        }
      }
    }
}

void printvar_ajsurfdiff(tNode *node, int vi)
{
  int f;
  tMesh *mesh = node->pat->mesh;
  char *name = VarName(vi);
  int zones =  MeshVarSurfacezones(mesh, vi);
  int type =  MeshVarType(mesh, vi);
  tDat *dat = node->dat;
  tArray *va = dat ? dat->v[vi] : NULL;
  char s[100];

  printf("%s Ind=%d type=%d zones=%d\n", name, vi, type, zones);

  if(va)
    for(f=0; f<6; f++)
    {
      tSurface *sf = dat->s[f][vi];
      tArray *msa  = sf ? sf->mysurf : NULL;
      tArray *asa  = sf ? sf->ajsurf : NULL;
      tArray **nsa = sf ? sf->nbsurf : NULL;
      if(msa)
      {
        printf("f%d allocd=%d mysurf\n", f, sf->allocd_mysurf);
      }
      if(nsa)
      {
        int ni;
        for(ni=0; ni<sf->nnbsurf; ni++)
        {
          printf("f%d fnb[f][%d]=%s allocd=%d nbsurf[%d]\n",
                 f, ni, nodename(node->fnb[f][ni],s,99),
                 sf->allocd_nbsurf[ni], ni);
        }
      }
      if(asa)
      {
        tArray *diff = alloc_array(asa->n);
        printf("f%d allocd=%d ajsurf, diff", f, sf->allocd_ajsurf);
        array_reldiff(diff, asa, msa);
        printarray(diff);
        free_array(diff);
      }
    }
}

void printvar_indc(tNode *node, int vi)
{
  int f;
  tMesh *mesh = node->pat->mesh;
  char *name = VarName(vi);
  int zones =  MeshVarSurfacezones(mesh, vi);
  int type =  MeshVarType(mesh, vi);
  tDat *dat = node->dat;
  tArray *va = dat ? dat->v[vi] : NULL;
  char s[100];

  printf("%s Ind=%d type=%d zones=%d\n", name, vi, type, zones);

  if(va)
  {
    tIndic *ic = dat->ic[vi];
    tArray *mia  = ic ? ic->myindc : NULL;
    if(mia)
    {
      printf("myindc");
      printarray(mia);
    }
    for(f=0; f<6; f++)
    {
      tArray **nia = ic ? ic->nbindc[f] : NULL;
      if(nia)
      {
        int ni;
        for(ni=0; ni<node->nfnb[f]; ni++)
        {
          printf("f%d fnb[f][%d]=%s allocd=%d nbindc[%d][%d]",
                 f, ni, nodename(node->fnb[f][ni],s,99),
                 ic->allocd_nbindc[f][ni], f, ni);
          printarray(nia[ni]);
        }
      }
    }
  }
}

/* print important arrays in a node */
void print_matrices_innode(tNode *n)
{
  int d;
  char s[100];

  printf("eid%lu:  %s  l%d  rflag=%d  datrank=%d  dat: %s\n",
          Elm_eid(n), nodename(n, s,99), Elm_l(n), n->rflag,
          n->datrank, n->dat ? "yes" : "no");
  printf(" ijk%d  [%g,%g]x[%g,%g]x[%g,%g]  np=%dx%dx%d=%d ",
          elm_get_ijk(n),
          n->bbox[0], n->bbox[1], n->bbox[2],
          n->bbox[3], n->bbox[4],n->bbox[5], n->n[0], n->n[1], n->n[2], n->np);

  printf(" arrays:\n");
  for(d=0; d<3; d++)
  {
    printf("Xb[%d]", d);
    printarray(node_Xb(n,d));
    printf("Wq[%d]", d);
    printarray(node_Wq(n,d));
    printf("WL[%d]", d);
    printarray(node_WL(n,d));

    printf("D=(Dt)^t: Dt[%d]", d);
    printarray(node_Dt(n,d));
    //printf("\n");
    printf("A=(At)^t: At[%d]", d);
    printarray(node_At(n,d));
    //printf("\n");
    printf("S=(St)^t: St[%d]", d);
    printarray(node_St(n,d));
    //printf("\n");
  }
  /*
  DECL_STACK_ARRAY(fa, node_Xb(n,0)->n);
  DECL_STACK_ARRAY(ca, node_Xb(n,0)->n);
  DECL_STACK_ARRAY(ga, node_Xb(n,0)->n);
  fa->d[0] = 1.;
  fa->d[1] = 0.25;
  fa->d[2] = 0.;
  fa->d[3] = 0.25;
  fa->d[4] = 1.;

  printf("fa");printarray(fa);
  mm_array0(node_At(n,0), fa, ca);
  printf("ca");printarray(ca);
  mm_array0(node_St(n,0), ca, ga);
  printf("ga");printarray(ga);
  */
}


/* print an array */
void printarray_sel(tArray *A, int dbl)
{
  int i,j,k;
  if(!A)
  {
    printf(" <NULL array>\n\n");
    return;
  }
  printf("->n[] = {%d,%d,%d}", A->n[0],A->n[1],A->n[2]);
  if(A->si) printf("  si=%d", A->si);
  printf("\n");
  if(A->range[0])
  {
    printf(" range =");
    for(j=0; j<A->n[1]; j++)
      printf(" %d,%d", A->range[0][j], A->range[1][j]);
    printf("\n");
  }
  for(k=0; k<A->n[2]; k++)
  {
    for(j=0; j<A->n[1]; j++)
    {
      for(i=0; i<A->n[0]; i++)
      {
        if(dbl) printf(" %g", A->d[Ind_n(i,j,k, A->n)]);
        else    printf(" %d", A->i[Ind_n(i,j,k, A->n)]);
      }
      printf("\n");
    }
    printf("\n");
  }
}
void printarray(tArray *A)
{
  printarray_sel(A, 1);
}
void printarray_int(tArray *A)
{
  printarray_sel(A, 0);
}

/* print a matrix in an array */
void printarray_matrix0(tArray *A)
{
  int i,J;
  int nJ = A->n[1] * A->n[2];
  int dJ = A->n[0];
  printf("->n[] = {%d,%d,%d} => %dx%d\n", A->n[0],A->n[1],A->n[2], A->n[0],nJ);
  for(i=0; i<A->n[0]; i++)
  {
    for(J=0; J<nJ; J++)
      printf(" %g", A->d[i + dJ*J]);
    printf("\n");
  }
}

/* print a matrix in an array */
void printarray_matrix1(tArray *A)
{
  int j,J;
  int nJ = A->n[0] * A->n[2];
  int dJ = A->n[1];
  printf("->n[] = {%d,%d,%d} => %dx%d\n", A->n[0],A->n[1],A->n[2], A->n[1],nJ);
  for(j=0; j<A->n[1]; j++)
  {
    for(J=0; J<nJ; J++)
      printf(" %g", A->d[j + dJ*J]);
    printf("\n");
  }
}

/* print an array */
void printarray_matrix2(tArray *A)
{
  int k,J;
  int nJ = A->n[1] * A->n[0];
  int dJ = A->n[2];
  printf("->n[] = {%d,%d,%d} => %dx%d\n", A->n[0],A->n[1],A->n[2], A->n[2],nJ);
  for(k=0; k<A->n[2]; k++)
  {
    for(J=0; J<nJ; J++)
      printf(" %g", A->d[k + dJ*J]);
    printf("\n");
  }
}

/* print an array */
void printarray_eploc(tArray *A, int details)
{
  int neplocs;
  int i;
  if(!A)
  {
    if(details) printf(" <NULL array>\n");
    return;
  }
  neplocs = array_Neplocs(A); //num. of nbs we already have

  if(details)
  {
    printf("->n[] = {%d,%d,%d}", A->n[0],A->n[1],A->n[2]);
    if(A->si) printf("  si=%d", A->si);
    printf(" => neplocs=%d", neplocs);
    printf("\n ");
  }

  for(i=0; i<neplocs; i++)
  {
    printeploc(&(A->eploc[i]));
    if(i<neplocs-1) printf(" ");
  }
  if(details) printf("\n");
}


/* print one bface */
void printthisbface(tBface *bface, const char *s)
{
  if(!bface) return;
  printf(" %s p%d f%d op%d", s, bface->pat->p,bface->f, bface->op);
  if(bface->brct_isset)
    printf(" [%g,%g]x[%g,%g]",
           bface->brct[0], bface->brct[1], bface->brct[2], bface->brct[3]);
  else
    printf(" [ , ]x[ , ]");
  printf(" ioC0_0=%d flags=%1d%1d\n",
         bface->ioC0_0,
         bface->face2, bface->boundary);
}

/* print one bface and its pair */
void printbface(tBface *bface)
{
  tBface *obface;

  if(!bface) return;

  obface = bface->obface;

  if(!obface)
  {
    printthisbface(bface, " ");
    return;
  }

  if(obface->obface == bface)
  {
    printthisbface(bface, "/");
    printthisbface(obface, "\\");
  }
  else
  {
    printf("  WARNING: bfaces are not properly linked!!\n");
    printf("  bface=%p bface->obface=%p:\n",
            (void *) bface, (void *) bface->obface);
    printthisbface(bface, " ");
    printf("  obface=%p obface->obface%p:\n",
           (void *) obface, (void *) obface->obface);
    printthisbface(obface, " ");
  }
}

/* print bfaces on face f */
void printbfaces_on_f(tPat *pat, int f)
{
  tBface *bf;

  printf("bfaces on f%d of pat->p=%d\n", f, pat->p);
  forbfacesonface(pat, f, bf)
    printbface(bf);
}

/* print all patch bfaces */
void printbfaces(tPat *pat)
{
  tBface *bf;

  printf("bfaces on pat->p=%d\n", pat->p);
  forbfaces(pat, bf)
    printbface(bf);
}

/* print all bfaces in mesh */
void printallbfaces(tMesh *mesh)
{
  int p;
  forpatches(mesh, p)
    printbfaces(mesh->pat[p]);
}

/* print 3-vec */
void pr3v(const char *s, const double x[3])
{
  printf("%s=%g %g %g  ", s, x[0],x[1],x[2]);
}

/* print an array of long ints */
void prlarray(const char *s, long n, const long *ar)
{
  long i;
  printf("%s[%ld] =", s, n);
  for(i=0; i<n; i++) printf(" %ld", ar[i]);
  printf("\n");
}

void prbbox(double *bb, int dim)
{
  int d;
  printf("[%g,%g]", bb[0], bb[1]);
  for(d=1; d<dim; d++)
    printf("x[%g,%g]", bb[2*d], bb[2*d+1]);
  printf(" ");
}

/* print corners of patch */
void printcorners(tPat *pat)
{
  double *bb = pat->bbox;
  int n[] = { 2,2,2 };
  int i,j,k;

  printf(" p%d\n", pat->p);
  forijk(i,j,k, n)
  {
    double X[] = { bb[i], bb[2+j], bb[4+k] };
    double x[3];

    set_xyz(pat, 0,-1, X, x);
    printf("   ijk%d%d%d", i,j,k); pr3v("x",x); pr3v("X",X);
    printf("\n");
  }
}

void printfacecorners(tPat *pat, int  f)
{
  double *bb = pat->bbox;
  int n[] = { 2,2,2 };
  int dir = f/2;
  int pl  = f%2;
  int i,j,k;

  printf("# p%d  f%d\n", pat->p, f);
  forplaneN(dir, i,j,k, n, pl)
  {
    double X[] = { bb[i], bb[2+j], bb[4+k] };
    double x[3];

    set_xyz(pat, 0,-1, X, x);
    printf("%g %g %g    %g %g %g\n", x[0],x[1],x[2], X[0],X[1],X[2]);
  }
}

/* print tRef contents */
void printref(tRef *ref)
{
  int d;
  printf("ref->  ");
  printf("type=%d  method=%d", ref->type, ref->method);
  printf("  n=");
  for(d=0; d<3; d++) printf("%d ", ref->n[d]);
  printf("\n");
}

/* print point numbers and types */
void print_n_pt_typ(tNode *node)
{
  int d;
  pr_nodename(node);
  printf("  n=");
  for(d=0; d<3; d++) printf("%d ", node->n[d]);
  printf(" pt_typ=");
  for(d=0; d<3; d++) printf("%d ", node->pt_typ[d]);
  printf("\n");
}

/* print info about one point */
void printpoint(tPoint *pt)
{
  tNode *node = pt->node;
  int ijk = pt->ijk;
  int *n = node->n;
  int k = kOfInd_n(ijk, n);
  int j = jOfInd_n_k(ijk, n,k);
  int i = iOfInd_n_jk(ijk, n,j,k);
  double X[3], x[3];

  XYZ_of_ind(node, ijk, X);
  set_xyz(NULL,node, ijk, X, x);

  pr_nodename(node);
  printf(" ijk=%d", ijk);
  //printf(" X=%g,%g,%g x=%g,%g,%g\n", X[0],X[1],X[2], x[0],x[1],x[2]);
  printf("=%d/%d,%d/%d,%d/%d x=%g,%g,%g\n", i,n[0], j,n[1], k,n[2],
         x[0],x[1],x[2]);
}

/* print var at one point */
void printvar_atpoint(tPoint *pt, int vi)
{
  tNode *node = pt->node;
  int ijk = pt->ijk;
  int np;
  tMesh *mesh;
  char *name;
  tDat *dat;
  tArray *va;

  if(!node)
  {
    printf(" node=NULL ijk=%d\n", ijk);
    return;
  }

  np = node->np;
  mesh = node->pat->mesh;
  name = VarName(vi);
  dat = node->dat;
  va = dat ? dat->v[vi] : NULL;

  pr_nodename(node);
  printf(" ijk=%d/%d: %s", ijk,np, name);
  if(ijk<np && ijk>=0)
  {
    if(va) printf("=%.18g\n", Arrd(va)[ijk]);
    else   printf(": NULL\n");
  }
  else
  {
    printf(": ijk outside range\n");
  }
}

/* print varlist at one point */
void printvarlist_atpoint(tPoint *pt, tVarList *vl)
{
  if(pt->node==NULL) return;
  if(vl)
  {
    int i, vi;

    PRF;printf(": nvars=%d\n", VLn(vl));
    forvl(vl, i)
    {
      vi = Vind(vl, i);
      printf(" %d vi=%d ", i, vi);
      printvar_atpoint(pt, vi);
    }
  }
  else
  {
    PRF;printf(": NULL\n");
  }
}

/* print all in nbranks */
void print_u32(khash_t(u32) *nbranks)
{
  khiter_t ki;
  printf("nbranks:");
  forkhiter(nbranks, ki)
    printf(" %u", kh_key(nbranks, ki));
  printf("\n");
}

/* print all in ef */
void print_ef(khash_t(u32_gptr) *ef)
{
  khiter_t ki;
  forkhiter(ef, ki)
  {
    struct list_head *fhead = kh_val(ef, ki);
    int f;
    printf("rank%u %p\n", kh_key(ef, ki), kh_val(ef, ki));
    for(f=0; f<6; f++)
    {
      printf("f%d ", f);
      printelmglist(&(fhead[f]));
    }
  }
}
