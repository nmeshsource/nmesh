/* variables.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"




/* add a variable to data base that is located in mesh */
void AddMeshVar(tMesh *mesh, char *name,
                char *tensorindices, char *description)
{
  tVar *newv;
  char fullname[100];
  int nilist;
  char *ilist[NINDEXLIST];
  int sym[3*NINDEXLIST];
  char *symsigns[] = {"-", "0", "+"};
  char **ss = symsigns+1;
  int i, j;

  if(0) {PRF;printf(": name=%s tensorindices=%s\n", name, tensorindices);}

  /* get list with tensor indices, this allocates mem in ilist[j] */
  tensorindexlist(tensorindices, &nilist, ilist, sym);

  /* for each tensor index */
  for(j = 0; j < nilist; j++)
  {
    /* construct name of variable */
    snprintf(fullname, 100, "%s%s", name, ilist[j]);
    free(ilist[j]); /* free string allocated in tensorindexlist */

    /* make sure that this variable does not exist yet */
    for(i = 0; i < mesh->nvdb; i++)
      if(!strcmp(mesh->vdb[i].name, fullname))
        break; /* this var is there already */

    if(i < mesh->nvdb) /* we found a var that exists already */
      errorexits("variable \"%s\" exists already\n", fullname);

    /* print name of variable */
    printf("  variable  %s\n", fullname);
    if(0) printf("%13s%s%s\n", ss[sym[3*j]], ss[sym[3*j+1]], ss[sym[3*j+2]]);

    /* variable does not exist, so add a new entry to data base */
    mesh->vdb = realloc(mesh->vdb, sizeof(tVar)*(mesh->nvdb+1));

    /* initialize and fill in structure */
    newv = &(mesh->vdb[mesh->nvdb]);
    memset(newv, 0, sizeof(tVar));
    newv->name          = strdup(fullname);
    newv->tensorindices = strdup(tensorindices);
    newv->description   = strdup(description);
    newv->index         = mesh->nvdb;
    newv->ncomponents   = nilist;
    newv->component     = j;
    newv->sym[0]        = sym[3*j];
    newv->sym[1]        = sym[3*j+1];
    newv->sym[2]        = sym[3*j+2];

    mesh->nvdb++;
  }

  /* ensure that allocation also happens in the dat structs of the nodes */
  realloc_meshvariables(mesh, mesh->nvdb);
}

/* free strings in vdb */
void free_mesh_vdb_contents(tMesh *mesh)
{
  tVar *var;
  int i;

  for(i = 0; i < mesh->nvdb; i++)
  {
    var = &(mesh->vdb[i]);
    free(var->name);
    free(var->tensorindices);
    free(var->description);
  }
}

/* add auxiliary variable to data base */
void AddAuxMeshVar(tMesh *mesh, char *name,
                   char *tensorindices, char *description)
{
  int nvdb  = mesh->nvdb;
  AddMeshVar(mesh, name, tensorindices, description);
  MeshVarSetType(mesh, nvdb, 1);
}

/* add variable with surfaces */
void AddEvoMeshVar(tMesh *mesh, char *name,
                   char *tensorindices, char *description)
{
  int nvdb  = mesh->nvdb;
  AddMeshVar(mesh, name, tensorindices, description);
  MeshVarSetSurfInfo(mesh, nvdb, 1);
}

/* add auxiliary variable to data base */
void AddMeshVarDim(tMesh *mesh, char *name,
                   char *tensorindices, char *description,
                   int n_special0, int n_special1, int n_special2)
{
  int nvdb  = mesh->nvdb;
  AddMeshVar(mesh, name, tensorindices, description);
  MeshVarSetSpecial(mesh, nvdb, n_special0, n_special1, n_special2);
}


/* return index of variable or -1 if it was not found */
int MeshVarIndLax(tMesh *mesh, char *name)
{
  tVar *vdb = mesh->vdb;
  int nvdb  = mesh->nvdb;
  int i, iS = mesh->vdb_iStart;

  if( (iS < 0) || (iS >= nvdb) ) iS = 0; /* make sure first i is in range */

  for(i = iS; i < nvdb; i++)
    if(!strcmp(vdb[i].name, name))
    {
      if(0) printf("index(%s) = %d\n", name, vdb[i].index);
      return vdb[i].index;
    }

  for(i = 0; i < iS; i++)
    if(!strcmp(vdb[i].name, name))
    {
      if(0) printf("index(%s) = %d\n", name, vdb[i].index);
      return vdb[i].index;
    }

  return -1;
}

/* return index of variable */
int MeshVarInd(tMesh *mesh, char *name)
{
  int i = MeshVarIndLax(mesh, name);
  if(i<0) errorexits("variable \"%s\" does not exist\n", name);
  return i;
}

/* return index of variable given pointer */
int VarIndFromPtr(tNode *node, double *p)
{
  tMesh *mesh = node->pat->mesh;
  tVar *vdb = mesh->vdb;
  int nvdb  = mesh->nvdb;
  tDat *dat = node->dat;
  int i;

  if(dat)
    for(i = 0; i < nvdb; i++)
      if(dat->v[i]->d == p)
        return vdb[i].index;
  return -1;
}

/* set the global var vdb_iStart to the index of variable "name" */
int Set_vdb_iStart_AtVar(tMesh *mesh, char *name)
{
  int i = MeshVarIndLax(mesh, name);
  if(i<0) errorexits("Ind: variable \"%s\" does not exist\n", name);
  mesh->vdb_iStart = i;
  return i;
}


/* return name given index */
char *MeshVarName(tMesh *mesh, int i)
{
  tVar *vdb = mesh->vdb;
  int nvdb  = mesh->nvdb;
  if (i < 0 || i >= nvdb)
    errorexit("VarName: index out of range");

  return vdb[i].name;
}

/* return number of components */
int MeshVarNComponents(tMesh *mesh, int i)
{
  tVar *vdb = mesh->vdb;
  int nvdb  = mesh->nvdb;
  if (i < 0 || i >= nvdb)
    errorexit("VarNComponents: index out of range");
  if (vdb[i].component != 0)
  {
    /* errorexit("VarNComponents: you have to use index of zeroth component"); */
    i = MeshVarIndComponent0(mesh, i);
  }
  return vdb[i].ncomponents;
}

/* return component */
int MeshVarComponent(tMesh *mesh, int i)
{
  tVar *vdb = mesh->vdb;
  return vdb[i].component;
}

/* return index of component 0 */
int MeshVarIndComponent0(tMesh *mesh, int i)
{
  tVar *vdb = mesh->vdb;
  return i - vdb[i].component;
}

/* return name of component 0 for a given name */
char *MeshVarNameComponent0(tMesh *mesh, char *name)
{
  return MeshVarName(mesh,
                     MeshVarIndComponent0(mesh, MeshVarInd(mesh, name)));
}

/* return string with tensor indices */
char *MeshVarTensorIndices(tMesh *mesh, int i)
{
  tVar *vdb = mesh->vdb;
  return vdb[i].tensorindices;
}

/* set information on how variable behaves at Boundary*/
void MeshVarNameSetBoundaryInfo(tMesh *mesh, char *name,
			        double farlimit, double falloff)
{
  tVar *vdb = mesh->vdb;
  int i = MeshVarInd(mesh, name);
  vdb[i].farlimit = farlimit;
  vdb[i].falloff = falloff;
}

/* set information on how variable behaves at Boundary*/
void MeshVarSetType(tMesh *mesh, int i, int type)
{
  tVar *vdb = mesh->vdb;
  int j, i0 = MeshVarIndComponent0(mesh, i);
  int n = MeshVarNComponents(mesh, i0);

  for(j = 0; j < n; j++)
  {
    vdb[j+i0].type = 1;
    if(0) printf("  setting %s type\n", vdb[j+i0].name);
  }
}

/* set information how many surface/ghost zones this var has */
void MeshVarSetSurfInfo(tMesh *mesh, int i, int surfacezones)
{
  tVar *vdb = mesh->vdb;
  int j, i0 = MeshVarIndComponent0(mesh, i);
  int n = MeshVarNComponents(mesh, i0);

  for(j = 0; j < n; j++)
  {
    vdb[i0+j].surfacezones = surfacezones;
  }
}
void MeshVarNameSetSurfInfo(tMesh *mesh, char *name, int surfacezones)
{
  int i = MeshVarInd(mesh, name);
  MeshVarSetSurfInfo(mesh, i, surfacezones);
}

void MeshVarSetSpecial(tMesh *mesh, int i,  int ns0, int ns1, int ns2)
{
  tVar *vdb = mesh->vdb;
  int j,d, i0 = MeshVarIndComponent0(mesh, i);
  int n = MeshVarNComponents(mesh, i0);
  int n_special[3] = { ns0, ns1, ns2 };

  for(j = 0; j < n; j++)
    for(d = 0; d < 3; d++)
    {
      vdb[i0+j].n_special[d] = n_special[d];
    }
}

/* return various pieces of information, e.g. boundary information */
double MeshVarFallOff(tMesh *mesh, int i)
{ return mesh->vdb[i].falloff; }

double MeshVarFarLimit(tMesh *mesh, int i)
{ return mesh->vdb[i].farlimit; }

int MeshVarSymmetry(tMesh *mesh, int i, int dir)
{ return mesh->vdb[i].sym[dir]; }

int MeshVarSurfacezones(tMesh *mesh, int i)
{ return mesh->vdb[i].surfacezones; }

int MeshVarType(tMesh *mesh, int i)
{ return mesh->vdb[i].type; }

int *MeshVar_n_special(tMesh *mesh, int i)
{ return mesh->vdb[i].n_special; }


/************************************************************************/
/* utility functions for variable lists */

/* print variable list */
void prvarlist(tVarList *v)
{
  tMesh *mesh = v->mesh;
  int i, j;

  printf("VarList=%p  time=%g  n=%d\n", v, v->time, v->n);
  for(i=0; i<v->n; i++)
  {
    j = v->index[i];
    printf(" %d  VarIndex=%d  %s\n", i, j, VarName(j));
  }
}

/* allocate an empty variable list */
tVarList *vlalloc(tMesh *mesh)
{
  tVarList *u;

  u = calloc(1, sizeof(tVarList));

  u->mesh = mesh;
  if(mesh) u->time = mesh->time;
  u->vlPars = NULL; /* set special pointer to NULL */
  return u;
}

/* free a variable list */
void vlfree(tVarList *u)
{
  if(u)
  {
    if(u->index) free(u->index);
    free(u);
  }
}

/* add a variable (one component) to a variable list */
void vlpushone(tVarList *v, int vi)
{
  v->n += 1;
  v->index = realloc(v->index, sizeof(int) * v->n);
  v->index[v->n-1] = vi;
}

/* add a variable with all its components to a variable list */
void vlpush(tVarList *v, int vi)
{
  tMesh *mesh = v->mesh;
  int i, n = MeshVarNComponents(mesh, vi);

  if(MeshVarIndComponent0(mesh, vi)!=vi)
    errorexit("vi needs to be index of component 0. "
              "Consider using vlpushone.");
  v->n += n;
  v->index = realloc(v->index, sizeof(int) * v->n);
  for(i = 0; i < n; i++)
    v->index[v->n-n+i] = vi + i;
}

/* add a variable list to a variable list */
void vlpushvl(tVarList *v, tVarList *u)
{
  int i;

  if(!v || !u) return;
  v->n += u->n;
  v->index = realloc(v->index, sizeof(int) * v->n);
  for(i = 0; i < u->n; i++)
    v->index[v->n - u->n + i] = u->index[i];
}

/* drop a variable (one component) from a variable list */
void vldropone(tVarList *v, int vi)
{
  int i;

  for(i = 0; i < v->n; i++)
    if(v->index[i] == vi)
    {
      v->n -= 1;
      for (; i < v->n; i++)
	v->index[i] = v->index[i+1];
      break;
    }
}

/* drop a variable with all its components from a variable list */
void vldrop(tVarList *v, int vi)
{
  tMesh *mesh = v->mesh;
  int i, n = MeshVarNComponents(mesh, vi);
  for(i = 0; i < n; i++) vldropone(v, vi+i);
}

/* drop last n variables from a variable list */
void vldropn(tVarList *v, int n)
{
  if(n <= 0) return;
  if(n >= v->n) v->n  = 0;
  else          v->n -= n;
}

/* duplicate variable list */
tVarList *vlduplicate(tVarList *v)
{
  int i;
  tVarList *u = vlalloc(v->mesh);

  u->time = v->time;
  for (i = 0; i < v->n; i++) vlpushone(u, v->index[i]);

  return u;
}

/* enable all variables in a variable list */
void vlenable(tVarList *v)
{
  enablevarlist(v);
}

void vlenablemesh(tMesh *mesh, tVarList *v)
{
  v->mesh = mesh;
  enablevarlist(v);
}

/* disable all variables in a variable list */
void vldisable(tVarList *v)
{
  disablevarlist(v);
}


/* create, enable, return pointer for a 1 variable VarList */
tVarList *VLPtrEnable1(tMesh *mesh, char *varname)
{
  tVarList *vl = vlalloc(mesh);
  int i = MeshVarInd(mesh, varname);

  enablevar(mesh, i);
  vlpush(vl, i);
  return vl;
}

/* disable variables in a VarList and free VarList */
void VLDisableFree(tVarList *vl)
{
  disablevarlist(vl);
  vlfree(vl);
}

/* add variables based on an existing variable list and a postfix
   note that we add each component as a scalar first but then fix it because
   we want gxx_p, gxy_p, ...  and not gxx_pxx, gxx_pxy ...
   We copy all properties if type<0, surfacezones<0, but if one is
   non-negative we set it this this value. */
tVarList *AddDuplicate(tVarList *vl, char *postfix, int type, int surfacezones)
{
  tMesh *mesh = vl->mesh;
  char name[1000];
  int i, j;
  tVarList *newvl;
  tVar *var, *newvar;

  if(mesh==NULL) errorexit("vl->mesh is NULL");

  /* new variable list with same number of indices */
  newvl = vlduplicate(vl);

  /* for all scalar variables in list */
  for(i = 0; i < vl->n; i++)
  {
    /* construct new name */
    var = &(mesh->vdb[vl->index[i]]);
    snprintf(name, 1000, "%s%s", var->name, postfix);

    /* if variable already exists, don't add it again */
    /* note that we nevertheless return a corresponding variable list */
    if ((j = MeshVarIndLax(mesh, name)) >= 0)
    {
      newvl->index[i] = j;
      continue;
    }

    /* add scalar variable with new name to variable database */
    AddMeshVar(mesh, name, "", var->description);

    /* get index of new variable and overwrite index in duplicate */
    newvl->index[i] = MeshVarInd(mesh, name);

    /* get pointer to old variable again since AddMeshVar reallocates vdb */
    var = &(mesh->vdb[vl->index[i]]);

    /* set structure in variable data base */
    newvar = &(mesh->vdb[newvl->index[i]]);
    free(newvar->tensorindices);
    newvar->tensorindices = strdup(var->tensorindices);
    newvar->component     = var->component;
    newvar->ncomponents   = var->ncomponents;
    newvar->farlimit      = var->farlimit;
    newvar->falloff       = var->falloff;
    newvar->type          = var->type;
    newvar->surfacezones  = var->surfacezones;
    for (j = 0; j < 3; j++)
    {
      newvar->sym[j]       = var->sym[j];
      newvar->n_special[j] = var->n_special[j];
    }
    /* set particular properties */
    if(type>=0)         newvar->type         = type;
    if(surfacezones>=0) newvar->surfacezones = surfacezones;
  }
  if(0) printf("mesh->nvdb is now %d\n", mesh->nvdb);

  return newvl;
}

/* add duplicate and enable variables */
tVarList *AddDuplicateEnable(tVarList *vl, char *postfix,
                             int type, int surfacezones)
{
  tVarList *newvl;

  newvl = AddDuplicate(vl, postfix, type, surfacezones);
  enablevarlist(newvl);
  return newvl;
}

/********************************************************************/
/* functions to deal with the var contents within the VarLists,
   some are node based, some work on the entire mesh */
/********************************************************************/

/* set: u = c */
void vlsetconstant(tVarList *u, const double c)
{
  tMesh *mesh = u->mesh;
  int myid;

  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    double *pu;
    int i, n, ui;

    for(n=0; n<u->n; n++)
    {
      ui = u->index[n];
      pu = GetVarDpointer(node, ui);
      forvari(node, ui, i)
        pu[i] = c;
    }
  }
}

/* copy: v = u on one node */
void vlcopy_node(tNode *node, tVarList *v, tVarList *u)
{
  double *pu, *pv;
  int i, n, vi,ui;

  for(n=0; n<v->n; n++)
  {
    ui = u->index[n];
    vi = v->index[n];
    pu = GetVarDpointer(node, ui);
    pv = GetVarDpointer(node, vi);
    forvari(node, vi, i)
      pv[i] = pu[i];
  }
}

/* copy: v = u */
void vlcopy(tVarList *v, tVarList *u)
{
  tMesh *mesh = v->mesh;
  int myid;

  /* copy time */
  v->time = u->time;

  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    vlcopy_node(node, v, u);
  }
}

void vlcopymesh(tMesh *mesh, tVarList *v, tVarList *u)
{
  if (!mesh || !v || !u) return;
  v->mesh = u->mesh = mesh;
  vlcopy(v, u);
}

/* wrapper for single variable: v = u (iv/u is index of v/u) */
void varcopy(tMesh *mesh, int iv, int iu)
{
  tVarList *v = vlalloc(mesh);
  tVarList *u = vlalloc(mesh);
  vlpushone(v, iv);
  vlpushone(u, iu);
  vlcopy(v, u);
  vlfree(u);
  vlfree(v);
}


/* swap v and u */
void vlswap(tVarList *v, tVarList *u)
{
  tMesh *mesh = v->mesh;
  int myid;
  double temp;

  /* swap time */
  temp = v->time;
  v->time = u->time;
  u->time = temp;

  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    double *pu, *pv;
    int i, n;

    for(n=0; n<v->n; n++)
    {
      int ui = u->index[n];
      int vi = v->index[n];
      pu = GetVarDpointer(node, ui);
      pv = GetVarDpointer(node, vi);
      forvari(node, vi, i)
      {
        temp  = pv[i];
        pv[i] = pu[i];
        pu[i] = temp;
      }
    }
  }
}

/* wrapper for single variable: swap u and v (iv/u is index of v/u) */
void varswap(tMesh *mesh, int iv, int iu)
{
  tVarList *v = vlalloc(mesh);
  tVarList *u = vlalloc(mesh);
  vlpushone(v, iv);
  vlpushone(u, iu);
  vlswap(v, u);
  vlfree(u);
  vlfree(v);
}

/* average: r=(a+b)/2 */
void vlaverage(tVarList *r, tVarList *a, tVarList *b)
{
  tMesh *mesh = r->mesh;
  int myid;
  double c = 0.5;

  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    double *pr, *pa, *pb;
    int i, n;

    for(n=0; n<r->n; n++)
    {
      int ri = r->index[n];
      int ai = a->index[n];
      int bi = b->index[n];
      pr = GetVarDpointer(node, ri);
      pa = GetVarDpointer(node, ai);
      pb = GetVarDpointer(node, bi);

      forvari(node, ri, i)
        pr[i] = c * (pa[i] + pb[i]);
    }
  }
  /* average times as well */
  r->time = c * (a->time + b->time);
}

/* subtract two var lists: r = a - b
   can be called as vlsubtract(r,a,b); or vlsubtract(a,a,b); */
void vlsubtract(tVarList *r, tVarList *a, tVarList *b)
{
  tMesh *mesh = r->mesh;
  int myid;

  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    double *pr, *pa, *pb;
    int i, n;

    for(n=0; n<r->n; n++)
    {
      int ri = r->index[n];
      int ai = a->index[n];
      int bi = b->index[n];
      pr = GetVarDpointer(node, ri);
      pa = GetVarDpointer(node, ai);
      pb = GetVarDpointer(node, bi);

      forvari(node, ri, i)
        pr[i] = pa[i] - pb[i];
    }
  }
  /* subtract times as well */
  r->time = a->time - b->time;
}

/* linear combination of two var lists: r = ca*a + cb*b
   catch special cases like ca=0 or cb=0
   if coefficient is zero memory is not accessed */
void vladd_node(tNode *node,
                tVarList *r, double ca, tVarList *a, double cb, tVarList *b)
{
  double *pr, *pa, *pb;
  int i, n;

  for(n=0; n<r->n; n++)
  {
    int ri = r->index[n];
    int ai = a->index[n];
    int bi = b->index[n];
    pr = GetVarDpointer(node, ri);
    if(ca!=0)  pa = GetVarDpointer(node, ai);
    if(cb!=0)  pb = GetVarDpointer(node, bi);

    if (ca == 0 && cb == 0)
    {
      forvari(node, ri, i) pr[i] = 0;
    }
    else if(ca == 0)
    {
      forvari(node, ri, i) pr[i] = cb * pb[i];
    }
    else if(cb == 0)
    {
      forvari(node, ri, i) pr[i] = ca * pa[i];
    }
    else
    {
      forvari(node, ri, i) pr[i] = ca * pa[i] + cb * pb[i];
    }
  }
}

/* add for all my leaf nodes */
void vladd(tVarList *r, double ca, tVarList *a, double cb, tVarList *b)
{
  tMesh *mesh = r->mesh;
  int myid;

  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    vladd_node(node, r, ca,a, cb,b);
  }
  /* add times as well */
  if (ca == 0 && cb == 0) r->time = 0.0;
  else if (ca == 0)	  r->time = cb * b->time;
  else if (cb == 0)	  r->time = ca * a->time;
  else			  r->time = ca * a->time + cb * b->time;
}


/* wrapper for single variable: r = ca*a + cb*b (ia/b/r is index of a/b/r) */
void varadd(tMesh *mesh, int ir, double ca, int ia, double cb, int ib)
{
  tVarList *a = vlalloc(mesh);
  tVarList *b = vlalloc(mesh);
  tVarList *r = vlalloc(mesh);
  vlpushone(a, ia);
  vlpushone(b, ib);
  vlpushone(r, ir);
  vladd(r, ca,a, cb,b);
  vlfree(a);
  vlfree(b);
  vlfree(r);
}

/* add second var list to first on one node: r += ca*a
   with special treatment for ca = 1 and ca = -1 */
void vladdto_node(tNode *node, tVarList *r, const double ca, tVarList *a)
{
  int i, n;

  if (ca == 0) return;

  /* loop over vars in varlists */
  for(n=0; n<r->n; n++)
  {
    int ri = r->index[n];
    int ai = a->index[n];
    double *pr = GetVarDpointer(node, ri);
    double *pa = GetVarDpointer(node, ai);

    if (ca == 1)
    {
      forvari(node, ri, i) pr[i] += pa[i];
    }
    else if(ca == -1)
    {
      forvari(node, ri, i) pr[i] -= pa[i];
    }
    else
    {
      forvari(node, ri, i) pr[i] += ca * pa[i];
    }
  }
}

/* add second var list to first: r += ca*a
   special treatment for ca = 1 and ca = -1 */
void vladdto(tVarList *r, const double ca, tVarList *a)
{
  tMesh *mesh = r->mesh;
  int myid;

  if (ca == 0) return;

  formylnodes(mesh, myid)
  {
    tNode *node = GetMyNode(mesh, myid);
    vladdto_node(node, r, ca, a);
  }

  /* add times as well */
  r->time += ca * a->time;
}
