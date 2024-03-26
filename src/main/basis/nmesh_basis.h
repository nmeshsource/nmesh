/* nmesh_basis.h */
/* (c) Wolfgang Tichy 1/2019 */
/* global header file */


/* types of interpolation we can do in 1d */
enum
{
  INTERP_NOT_SET=0, /* interpolation mode is not set */
  INTERP_LAGRANGE,  /* Lagrange interpolation */
  INTERP_WENO,      /* WENO interpolation */
};


/* structure that holds global coordinates vars and pars */
typedef struct {
  int expfilter_JacobianPower; /* Par("basis_expfilter_JacobianPower") */
  int filter_fv;               /* Par("basis_filter_fv") */
} tbasis;


/* types of grid points we can have in 1d */
enum
{
  P_LGL,     /* Legendre Gauss-Lobatto points */
  P_UNIFORM, /* uniform grid spacing */
  P_NTYPES   /* number of point types */
};


/* structure to hold points, weights and matrices for each point type */
typedef struct tGRIDPOINTS {
  int nmax;                     /* max n we can have */
  struct tARRAY **Xb[P_NTYPES]; /* list of points (often LGL in [-1,1]) */
  struct tARRAY **Wq[P_NTYPES]; /* list of quadrature weights for Xb */
  struct tARRAY **WL[P_NTYPES]; /* list of Lagrange interp. weights */
  struct tARRAY **Dt[P_NTYPES]; /* list of transp. differentiation matrices
                                 we store Dt[typ][1...nmax], typ=P_LGL,... */
  struct tARRAY **Dpt[P_NTYPES]; /* transp. forward diff matrices */
  struct tARRAY **Dmt[P_NTYPES]; /* transp. backward diff matrices */
  struct tARRAY **At[P_NTYPES]; /* list of transposed analysis matrices */
  struct tARRAY **St[P_NTYPES]; /* list of transposed synthesis matrices */
  double (*basis[P_NTYPES])(int l, double Xb, int np);//basis related to At,St
} tGridPoints;


/* structure to hold options for deriv taking */
typedef struct tDERIVOPT {
  int lop;  /* lopsided stencil?: 0=centered, 1=forward, -1=backward fd */
} tDerivOpt;



/* global functions */

/* basis.c */
void basis_array_deriv1(tNode *node, int dir, tArray *var, tArray *dvar,
                        tDerivOpt *opt);
void basis_array_derivs(tNode *node, tArray *var, tArray *dvar[3],
                        tDerivOpt *opt);
int basis_var_deriv1(tNode *node, int dir, int vi, int dvi, tDerivOpt *opt);
int basis_var_analysis3(tNode *node, int ui, int ci);
int basis_var_synthesis3(tNode *node, int ui, int ci);
void basis_array_analysis3(tNode *node, tArray *u, tArray *c);
void basis_array_synthesis3(tNode *node, tArray *u, tArray *c);
void basis_array_analysis3_At(tArray *At[3], tArray *u, tArray *c);
void basis_array_synthesis3_St(tArray *St[3], tArray *u, tArray *c);
double basis_array_interpolate(tNode *node, tArray *coef, double Xb[3]);
tArray *array_GLquadrature1(tNode *node, int dir, tArray *var, tArray *Ivar);
double array_GLquadrature3(tNode *node, tArray *var);
double array_nodeaverage(tNode *node, tArray *var);
double var_GLquadrature3(tNode *node, int ui);
double array_GLquadratureXYZ3(tNode *node, tArray *var);
double var_nodeaverage(tNode *node, int ui);
double var_GLquadratureXYZ3(tNode *node, int ui);

/* Legendre.c */
double basis_LegendreP(int l, double x, int np);
double basis_normLegendreP(int i, double x, int np);
void LGL_x_wquad(int npoints, double *x, double *w);
void LGL_AT_ST_matrices(int n, const double *x, const double *w,
                        double *AT, double *ST);
void LG_x_wquad(int npoints, double *x, double *w);
void LG_set_Xb_Wq(tArray *Xb, tArray *Wq);
void LG_2Sphere_get_zi_theta_phi(tArray *Zb, int nphi, int i, int j,
                                 double *z_i, double *theta_i, double *phi_j);
void LG_2SphereIntegral(tArray *auijk, tArray *Wq, tArray *aUk);
void Gauss_wquad_from_symm_x(int npoints, const double *x, double *w);
void uniform_x_wGaussquad(int npoints, double *x, double *w);
void uniform_x_wTrapez(int npoints, double *x, double *w);
void Legendre_AT_ST_matrices(int n, const double *x, const double *w,
                             double *AT, double *ST);

/* Lagrange.c */
void Lagrange_winterp(int n, const double *x, double *w_interp);
void Lagrange_DT(int n, const double *x, const double *w_interp, double *DT);
double Lagrange_of_x(int k, double x, int np,
                     const double *x_p, const double *w_interp);
double Lagrange_prod1(int l, double x, int np, const double *x_p);
double Lagrange_prod2(int l, int m, double x, int np, const double *x_p);
void fd_lopderiv_DT_uniform(int n, const double *x, int ssz, int lop,
                            double *DT, int *range[2]);

/* interpolate.c */
double basis_pw_const(int k, double x, int np,
                      const double *x_p, const double *w_interp);
double basis_pw_linear(int k, double x, int np,
                       const double *x_p, const double *w_interp);
double basis_pw_parab(int k, double x, int np,
                      const double *x_p, const double *w_interp);
double basis_array_interp(tNode *node, tArray *var, double Xb[3],
                          double Basis(int k, double x, int np,
                                       const double *x_p,
                                       const double *w_interp));
double basis_array_interp2d(tNode *node, tArray *var, int dir, int p,
                            double Cb[2],
                            double Basis(int k, double x, int np,
                                         const double *x_p,
                                         const double *w_interp));
void array_1d1d1d_coords_to_3d_coords(tArray *X1d[3], tArray *Xp[3]);
void fill_3arrays_with_nodepoints(tNode *node, tArray *Xp[3]);
void fill_2arrays_with_nodepoints(tNode *node, int dir, tArray *Cp[2]);
void basis_interp_topoints(tNode *node, tArray *var,
                           tArray *Xp[3], tArray *interp,
                           double Basis(int k, double x, int np,
                                        const double *x_p,
                                        const double *w_interp));
void basis_interp_toIpoints(tNode *node, tArray *var,
                            tArray *Xp[3], tArray *Ip, tArray *interp,
                            double Basis(int k, double x, int np,
                                         const double *x_p,
                                         const double *w_interp));
void basis_interp2d_topoints(tNode *node, tArray *var, int dir, int p,
                             tArray *Cp[2], tArray *interp,
                             double Basis(int k, double x, int np,
                                          const double *x_p,
                                          const double *w_interp));
void basis_interp2d_toIpoints(tNode *node, tArray *var, int dir,int p,
                              tArray *Cp[2], tArray *Ip,
                              tArray *interp,
                              double Basis(int k, double x, int np,
                                           const double *x_p,
                                           const double *w_interp));
void insert_array_inplane(tArray *var, int dir, int p, tArray *interp2d);
double interp_to_Xb0(tElm *elm, tArray *var, double Xb0[3], int np[3],
                     int scheme, double vscal);
void interpolate_topoints(tElm *elm, tArray *var, tArray *Xp[3],
                          int np[3], int scheme, double vscal,
                          tArray *interp);
void interpolate_toIpoints(tElm *elm, tArray *var, tArray *Xp[3], tArray *Ip,
                           int np[3], int scheme, double vscal,
                           tArray *interp);
void interp_topoints(tElm *elm, tArray *var, tArray *Xp[3],
                     int npts, int scheme, double vscal,
                     tArray *interp);
void interp_toIpoints(tElm *elm, tArray *var, tArray *Xp[3], tArray *Ip,
                      int npts, int scheme, double vscal, tArray *interp);
void interp_to_pt_typ(tNode *node, int iu, int pt_typ[3],
                      int npts, int scheme, double vscal, tArray *interp);
int interpolate_var_ok(tNode *node, int vi, double Xb[3],
                       int npts, int scheme, double *vinterp);
double basis_var_interp_xyz(tMesh *mesh, int ivar, double xyz[3],
                            double Basis(int k, double x, int np,
                                         const double *x_p,
                                         const double *w_interp));
int interp_VL_xp(tMesh *mesh, tVarList *vl, tArray *xp[3],
                 int npts, int scheme, double vscal, tArray *Value);
double interp_var_xyz(tMesh *mesh, int ivar, const double xyz[3],
                      int np, int scheme, double vscal);
double interp_var_x_y_z(tMesh *mesh, int ivar, double x,double y,double z,
                        int np, int scheme, double vscal);

/* SphericalHarmonics.c */
double *alloc_Plm_Tab(int lmax);
void set_YlmTabs(int lmax, double th, double ph, double *ReYtab, double *ImYtab);
void Ylm_from_Tabs(int lmax, double *ReYtab, double *ImYtab, int l, int m,
                  double *ReYlm, double *ImYlm);
void Re_Im_Ylm(int l, int m, double theta, double phi,
               double *ReYlm, double *ImYlm);
void SphHarm_dphi_forRealFunc(double *c, double *cdphi, int lmax);
void SphHarm_sin_theta_dtheta_forRealFunc(double *c, double *csdth, int lmax);

/* get_coords.c */
void nearest_ijk_of_XYZ(tNode *node, int ijk[3], const double X0[3]);

/* filter.c */
int expfilter_var(tNode *node, int ui, double alp[3], double s[3]);
int has_expfalloff_var(tNode *node, int ui, double alp[3], double s[3]);
void expfilter_vl(tVarList *vl, double af, double sf);
double linear_fit_result(double beta[4], int i, int j, int k);
void unfiltered_range_of_expfilter(int n[3], double alp[3], double s[3],
                                   double f_unfilt, int n_unfilt[3]);
double fit_n_coefflogs(tArray *ca, int n_fit[3], double beta[4]);
double fit_unfiltered_coefflogs(tArray *ca, double alp[3], double s[3],
                                double f_unfilt, int n_unfilt[3],
                                double beta[4]);
int topcoeff_has_expfalloff_array(tNode *node, tArray *ua,
                                  double alp[3], double s[3],
                                  double f_unfilt, double fac);
int topcoeff_has_expfalloff_var(tNode *node, int ui,
                                double alp[3], double s[3],
                                double f_unfilt, double fac);

/* SpinWeightedSphericalHarmonics.c */
double Re_sYlm(int l, int m, int s, double theta, double phi);
double Im_sYlm(int l, int m, int s, double theta, double phi);
void ReIm_sYlm(int l, int m, int s, double theta, double phi,
               double *Re_sYlm, double *Im_sYlm);

/* sphericalDF.c */
void sphericalDF_theta_phi(int i, int j, const int *n,
                           double *theta_i, double *phi_j);
void sphericalDF_2dIntegral(tArray *auijk, tArray *aUk);
void sphericalDF_copy_to_doubleCoveredPoints(tArray *AsDF);
