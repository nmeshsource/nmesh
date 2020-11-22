/* nmesh_dg.h */
/* Wolfgang Tichy, April 2019 */



/* reconstruction modes that we can have in DGglobals->fv_rec_mode */
enum
{
  FV_REC_0=0,              /* reconstruction mode not set */
  FV_REC_1,                /* 1st order accurate reconstruction */
  FV_REC_WENO3if2away_1,   /* WENO3 if 2 away from face, else 1st order */
  FV_REC_WENO3if1away_1,   /* WENO3 if 1 away from face, else 1st order */
  FV_REC_WENO3_2           /* WENO3 if 1 away from face, else 2nd? order*/
};

/* structure that holds global dg pars */
typedef struct {
  /* frequently used pars */
  int fv_rec_mode;  /* reconstruction mode based on par fv_rec */
} tDGglobals;



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
void printDGinfo(tDGinfo *d);

/* rec1d.c */
double rec1d_p_1(int n, const double *u, int im, double u_scale);
double rec1d_m_1(int n, const double *u, int im, double u_scale);
double rec1d_p_WENO3_if2away(int n, const double *u, int im, double u_scale);
double rec1d_m_WENO3_if2away(int n, const double *u, int im, double u_scale);
double rec1d_p_WENO3_if1away(int n, const double *u, int im, double u_scale);
double rec1d_m_WENO3_if1away(int n, const double *u, int im, double u_scale);
double rec1d_p_WENO3_2(int n, const double *u, int im, double u_scale);
double rec1d_m_WENO3_2(int n, const double *u, int im, double u_scale);
