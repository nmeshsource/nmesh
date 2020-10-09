/* nmesh_dg.h */
/* Wolfgang Tichy, April 2019 */


/* structure that contains info needed for DG flux */
typedef struct tDGINFO {
  tNode *node;
  int face;
  int i;
  int j;
  int k;
  tVarList *vlu;  /* varlist with cons. vars */
  double *ui;     /* array with cons vars inside node at i,j,k */
  double *ua;     /* array with cons vars on adjacent side of i,j,k */
  double *fi;     /* array of inside fluxes */
  double *fa;     /* array of adjacent fluxes */
  double *lami;   /* array of inside eigenvals */
  double *lama;   /* array of adjacent eigenvals */
  double *fnum;   /* array of numerical fluxes f^* = f^i* n_i */
  double Ffac;    /* factor (usually 1) by which we multiply surface terms */
  tVarList *vls;  /* varlist with sources needed */
  int info;       /* user defined extra info, 0 means default */
} tDGinfo;


/* numflux1d.c */
void numflux1d_scalarGodunov(tDGinfo *d);
void numflux1d_upwind(tDGinfo *d);
void numflux1d_LLF(tDGinfo *d);
void numflux1d_HLL(tDGinfo *d);

/* dg.c */
tDGinfo *alloc_DGinfo(tVarList *vlu, tVarList *vls);
void free_DGinfo(tDGinfo *dgi);
int dg_add_surface_fluxes(tMesh *mesh, tVarList *vlr, tVarList *vlu,
                          tVarList *vls,
                          void (*u_f_lam)(tDGinfo *d),
                          void (*numflux)(tDGinfo *d));

/* rec1d.c */
double rec1d_p_0(int n, const double *u, int im);
double rec1d_m_0(int n, const double *u, int im);
