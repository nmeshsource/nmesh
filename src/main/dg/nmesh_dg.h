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
  double *ui;     /* arrays with cons vars inside node at i,j,k */
  double *ua;     /* arrays with cons vars on adjacent side of i,j,k */
  double *fi;     /* inside fluxes */
  double *fa;     /* adjacent fluxes */
  double *lami;   /* inside eigenvals */
  double *lama;   /* adjacent eigenvals */
  double *fnum;   /* numerical flux f^* */
} tDGinfo;


/* numflux1d.c */
void numflux1d_scalarGodunov(tDGinfo *d);
void numflux1d_upwind(tDGinfo *d);
void numflux1d_LLF(tDGinfo *d);
void numflux1d_HLL(tDGinfo *d);

/* dg.c */
int dg_add_surface_fluxes(tMesh *mesh, tVarList *vlr, tVarList *vlu,
                          void (*u_f_lam)(tDGinfo *d),
                          void (*numflux)(tDGinfo *d));
