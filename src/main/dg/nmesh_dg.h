/* nmesh_dg.h */
/* Wolfgang Tichy, April 2019 */



/* reconstruction modes that we can have in DGglobals->fv_rec_mode */
enum
{
  FV_REC_0=0,             /* reconstruction mode not set */
  FV_REC_1,               /* 1st order accurate reconstruction */
  FV_REC_WENOm3_2,        /* non-standard WENO3 with uniform opt. weights */
  FV_REC_WENOm5_2,        /* non-standard WENO5 with uniform opt. weights */
  FV_REC_WENOmZ_2,        /* non-standard WENOZ with uniform opt. weights */
  FV_REC_WENO3if2away_1,  /* WENO3 if 2 away from face, else 1st order */
  FV_REC_WENO3if1away_1,  /* WENO3 if 1 away from face, else 1st order */
  FV_REC_WENO3_2,         /* WENO3 if 1 away from face, else 2nd? order */
  FV_REC_WENO3_2g         /* WENO3 if 1 away, else 2nd? order w. ghost */
};

/* structure that holds global dg pars */
typedef struct {
  /* frequently used par values */
  double outerBC_flux_fac[3]; /* values from par dg_outerBC_flux_fac */
  int fv_rec_mode;            /* reconstruction mode based on par fv_rec */
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
  double *si;     /* array with source vars inside node at i,j,k */
  double *sa;     /* array with source vars on adjacent side of i,j,k */
  int info;       /* extra info, 0=DGINFO_NULL means default */
} tDGinfo;

/* meaning of bits in info field of DGINFO: */
enum
{
  DGINFO_NULL  = 0,  /* default */
  DGINFO_MIDPT = 1,  /* i,j,k is interpreted as a midpoint */
};


/* structure that contains info needed to reconstruct cons vars u */
typedef struct tFVINFO {
  int nq;        /* num of q-vars, (q can be u) */
  double **qc;   /* qc[l][i0] = val. of var q_l at grid point i0 */
  int npts;      /* number of grid points, i.e. n[dir] */
  int q_scale;   /* scale of q */
  double (*rec1d_p)(int n, const double *q, int im, double q_scale);
  double (*rec1d_m)(int n, const double *q, int im, double q_scale);
  int im;        /* midpoint where we reconstruct vars */
  double *qm_p;  /* qm_p[l] is reconstructed q in positive direction (p) */
  double *qm_m;  /* qm_m[l] is reconstructed q in negative direction (m) */
} tFVinfo;




/* numflux1d.c */
void numflux1d_scalarGodunov(tDGinfo *d);
void numflux1d_upwind(tDGinfo *d);
void numflux1d_LLF(tDGinfo *d);
void numflux1d_HLL(tDGinfo *d);

/* dg.c */
tDGinfo *alloc_DGinfo(tVarList *vlu, tVarList *vls);
void free_DGinfo(tDGinfo *dgi);
int dg_add_surface_fluxes_sign(tMesh *mesh, double sign, tVarList *vldf,
                               tVarList *vlu, tVarList *vls,
                               void (*u_f_lam)(tDGinfo *d),
                               void (*numflux)(tDGinfo *d));
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
double rec1d_p_WENO3_2g(int n, const double *u, int im, double u_scale);
double rec1d_m_WENO3_2g(int n, const double *u, int im, double u_scale);
void rec1d_uface_to_uin_1(tNode *node, tVarList *vlu, int forward);
void rec1d_uface_to_uin_1_mesh(tMesh *mesh, tVarList *vlu, int forward);
double rec1d_p_WENOm3_2(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOm3_2(int n, const double *u, int im, double u_scale);
double rec1d_p_WENOm5_2(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOm5_2(int n, const double *u, int im, double u_scale);
double rec1d_p_WENOmZ_2(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOmZ_2(int n, const double *u, int im, double u_scale);

/* fv.c */
void fv_rec1d_q_midpt(tFVinfo *fv);
void fv_divf(tNode *node, tVarList *vldivf, tVarList *vlq, tVarList *vlu,
             void (*rec1d_u_f_lam_midpt)(tFVinfo *f, tDGinfo *d),
             void (*numflux)(tDGinfo *d));
