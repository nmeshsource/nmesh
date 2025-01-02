/* nmesh_dg.h */
/* Wolfgang Tichy, April 2019 */



/* reconstruction modes that we can have in DGglobals->fv_rec_mode,
   as well as extrapolation modes in DGglobals->fv_divf_extrap_mode */
enum
{
  FV_REC_0=0,             /* reconstruction mode not set */
  FV_REC_1,               /* 1st order accurate reconstruction */
  FV_REC_WENOm3_2,        /* non-standard WENO3 with uniform opt. weights */
  FV_REC_WENOm5_2,        /* non-standard WENO5 with uniform opt. weights */
  FV_REC_WENOmT_2,        /* similar to WENOmZ_2 */
  FV_REC_WENOmZ_2,        /* non-standard WENOZ with uniform opt. weights */
  FV_REC_WENOm3_1,        /* non-standard WENO3 with uniform opt. weights */
  FV_REC_WENOm5_1,        /* non-standard WENO5 with uniform opt. weights */
  FV_REC_WENOmZ_1,        /* non-standard WENOZ with uniform opt. weights */
  FV_REC_WENOm3,          /* WENOm3_2 or WENOm3_1 */
  FV_REC_WENOm5,          /* WENOm5_2 or WENOm5_1 */
  FV_REC_WENOmZ,          /* WENOmZ_2 or WENOmZ_1 */
  FV_REC_WENO3if2away_1,  /* WENO3 if 2 away from face, else 1st order */
  FV_REC_WENO3if1away_1,  /* WENO3 if 1 away from face, else 1st order */
  FV_REC_WENO3_2,         /* WENO3 if 1 away from face, else 2nd? order */
  FV_REC_WENO3_2g,        /* WENO3 if 1 away, else 2nd? order w. ghost */
  FV_NO_EXTRAP,           /* do not iterpolate div(flux) */
  FV_DIVF_EXTRAP1,        /* extrapolate all of div(flux) = d_i f^i */
  FV_DNFN_EXTRAP1,        /* in d_i f^i extrap. only term along face normal */
  FV_2DINTERP_LINEAR,     /* use linear interp. on fv node surface */
  FV_2DINTERP_PARAB       /* use parabolic interp. on fv node surface */
};

/* structure that holds global dg pars */
typedef struct {
  /* frequently used par values */
  double outerBC_flux_fac[3]; /* values from par dg_outerBC_flux_fac */
  int fv_rec_mode;            /* reconstruction mode based on par fv_rec */
  int fv_rec_mode_WENOm;      /* 1 if WENOm3, WENOm5, or WENOmZ*/
  double fv_rec_WENOm_s1;     /* value of par fv_rec_WENOm_s1 */
  double fv_rec_WENOm_s2;     /* value of par fv_rec_WENOm_s2 */
  int fv_rec_WENOm_opt;       /* value of par fv_rec_WENOm_opt */
  int fv_flux_is_fnum_minus_fi; //whether we subtract inner from num. flux
  int fv_divf_extrap_mode;    /* div(f) extrap. mode from fv_divf_extrap */
  int fv_divf_adds_surface_fluxes; //if 1 fv_divf adds surf. fluxes
  int fv_surface_interp_mode; /* surf. interp. mode from fv_surface_interp */
  double fv_WENOm3_optw[2];   /* optimal weights for WENOm3 */
  double fv_divf_extrap_s1;   /* value of par fv_divf_extrap_s1 */
  double fv_divf_extrap_s2;   /* value of par fv_divf_extrap_s2 */
  int fv_divf_extrap_opt;     /* value of par fv_divf_extrap_opt */
  int fv_divf_use_only_right_flux; //Boolean val of par fv_divf_use_only_right_flux
  int fv2dg_interp_use_extrap1;    //Boolean val of par fv2dg_interp_use_extrap1
  tVarList *fv2dg_interp_use_extrap1_vl; //vl for fv2dg_interp_use_extrap1
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
  DGINFO_NULL      = 0,  /* default */
  DGINFO_MIDPTNORM = 1,  /* use midpoint normal */
  DGINFO_INNONLY   = 2,  /* set lami, ui, fi only */
  DGINFO_ADJONLY   = 4,  /* set lama, ua, fa only */
};


/* structure that contains info needed to reconstruct cons vars u */
typedef struct tFVINFO {
  int nq;        /* num of q-vars, (q can be u) */
  double **qc;   /* qc[l][i0] = val. of var q_l at grid point i0 */
  int npts;      /* number of grid points, i.e. n[dir] */
  double q_scale; /* scale of q */
  double (*rec1d_p)(int n, const double *q, int im, double q_scale);
  double (*rec1d_m)(int n, const double *q, int im, double q_scale);
  int im;        /* midpoint where we reconstruct vars */
  double *qm_p;  /* qm_p[l] is reconstructed q in positive direction (p) */
  double *qm_m;  /* qm_m[l] is reconstructed q in negative direction (m) */
  unsigned *stat; //rec1d_u_f_lam_midpt can write a status for each qm-var
} tFVinfo;


/* possible status bits written into stat of tFVinfo by the func
   rec1d_u_f_lam_midpt passed to fv_divf */
enum
{
  FV_REC_OK               = 0,  //must be zero (because fv_divf uses memset for init)
  FV_REC_NO_LEFT_EXTRAP1  = 2,  //do not extrap rhs on left
  FV_REC_NO_RIGHT_EXTRAP1 = 4,  //do not extrap rhs on right
  FV_REC_SHARP_ON_LEFT    = 8,  //qc sharp on left
  FV_REC_SHARP_ON_RIGHT   = 16  //qc sharp on right
};


/* numflux1d.c */
void numflux1d_scalarGodunov(tDGinfo *d);
void numflux1d_upwind(tDGinfo *d);
void numflux1d_LLF(tDGinfo *d);
void numflux1d_HLL(tDGinfo *d);

/* dg.c */
int dg_add_surface_fluxes_sign_fvflag(tNode *node, double sign,
                          tVarList *vldf, tVarList *vlu, tVarList *vls,
                          void (*u_f_lam)(tDGinfo *d),
                          void (*numflux)(tDGinfo *d),
                          int use_fv);
int dg_add_surface_fluxes_sign(tNode *node, double sign, tVarList *vldf,
                               tVarList *vlu, tVarList *vls,
                               void (*u_f_lam)(tDGinfo *d),
                               void (*numflux)(tDGinfo *d));
int dg_add_surface_fluxes(tNode *node, tVarList *vlr, tVarList *vlu,
                          tVarList *vls,
                          void (*u_f_lam)(tDGinfo *d),
                          void (*numflux)(tDGinfo *d));
tDGinfo *alloc_DGinfo(tVarList *vlu, tVarList *vls);
void free_DGinfo(tDGinfo *dgi);
void copy_nonallocd_DGinfo(tDGinfo *dsrc, tDGinfo *ddest);
void printDGinfo(tDGinfo *d);
double node_normal_from_DGinfo(tDGinfo *d, double nrm[3]);

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
void rec1d_uface_to_uin_1_Carray(int n, double *u, int forward,
                                 double u_scale, double s1, double s2,
                                 int opt);
void rec1d_uface_to_uin_1(tNode *node, tVarList *vlu, int forward);
void rec1d_uface_to_uin_1_mesh(tMesh *mesh, tVarList *vlu, int forward);
void rec1d_fv_uface_to_uin_1_if_rflag(tMesh *mesh, tVarList *vlu,
                                      int forward);
void rec1d_uface_to_uin_1_var(tNode *node, int vi, int forward);
double rec1d_p_WENOm3_2(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOm3_2(int n, const double *u, int im, double u_scale);
double rec1d_p_WENOm3_1(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOm3_1(int n, const double *u, int im, double u_scale);
double rec1d_u_in1_weightfac(int n, const double *u, int i0, double u_scale,
                             double s1, double s2, int opt);
double rec1d_p_WENOm3(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOm3(int n, const double *u, int im, double u_scale);
double rec1d_p_WENOm5_2(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOm5_2(int n, const double *u, int im, double u_scale);
double rec1d_p_WENOm5_1(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOm5_1(int n, const double *u, int im, double u_scale);
double rec1d_p_WENOm5(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOm5(int n, const double *u, int im, double u_scale);
double rec1d_p_WENOmT_2(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOmT_2(int n, const double *u, int im, double u_scale);
double rec1d_p_WENOmZ_2(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOmZ_2(int n, const double *u, int im, double u_scale);
double rec1d_p_WENOmZ_1(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOmZ_1(int n, const double *u, int im, double u_scale);
double rec1d_p_WENOmZ(int n, const double *u, int im, double u_scale);
double rec1d_m_WENOmZ(int n, const double *u, int im, double u_scale);

/* fv.c */
void fv_rec1d_q_midpt(tFVinfo *fv);
void fv_stat_WENOm_1or2(tFVinfo *fv);
void fv_rec1d_q_midpt_WENOm_1or2(tFVinfo *fv);
void fv_divf(tNode *node, tVarList *vldivf, tVarList *vlq,
             tVarList *vlu, tVarList *vls,
             void (*rec1d_u_f_lam_midpt)(tFVinfo *f, tDGinfo *d),
             void (*u_f_lam)(tDGinfo *d),
             void (*numflux)(tDGinfo *d));
void printFVinfo(tFVinfo *fv);

/* dissipation.c */
void dissipation_add_KO4(tNode *node, tVarList *vlr, tVarList *vlu,
                         double dissfac);
void dissipation_add_KO4_mesh(tMesh *mesh, tVarList *vlr, tVarList *vlu,
                              double dissfac);
void dissipation_add_KO_order(tNode *node, tVarList *vlr, tVarList *vlu,
                              double dissfac, int order);
void dissipation_add_taperedKO_order_cf(tNode *node, tVarList *vlr,
                                        tVarList *vlu, double dissfac,
                                        int order, double *cf);
void dissipation_add_taperedKO_order(tNode *node, tVarList *vlr, tVarList *vlu,
                                     double dissfac, int order);
void dissipation_add_taperedKO_order_min(tNode *node,
                                         tVarList *vlr, tVarList *vlu,
                                         double dissfac, int order,
                                         int min_order);
void dissipation_add_WTmodKO_order(tNode *node, tVarList *vlr, tVarList *vlu,
                                   double dissfac, int order);
