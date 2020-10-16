/* nmesh_basis.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for global functions */


/* basis.c */
void basis_array_derivs(tNode *node, tArray *var, tArray *dvar[3]);
int basis_var_derivs(tNode *node, int vi, int dvi[3]);
int basis_var_analysis3(tNode *node, int ui, int ci);
int basis_var_synthesis3(tNode *node, int ui, int ci);
void basis_array_analysis3(tNode *node, tArray *u, tArray *c);
void basis_array_synthesis3(tNode *node, tArray *u, tArray *c);
double basis_array_interpolate(tNode *node, tArray *coef, double Xb[3]);
double basis_var_interpolate(tNode *node, int vi, double Xb[3]);
//tNode *basis_var_interpolate_mesh(tMesh *mesh, int ui, const double x[3],
//                                  double *val);
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
void Gauss_wquad_from_symm_x(int npoints, const double *x, double *w);
void uniform_x_wGaussquad(int npoints, double *x, double *w);
void Legendre_AT_ST_matrices(int n, const double *x, const double *w,
                             double *AT, double *ST);

/* Lagrange.c */
void Lagrange_winterp(int n, const double *x, double *w_interp);
void Lagrange_DT(int n, const double *x, const double *w_interp, double *DT);
double Lagrange_of_x(int k, double x, int np,
                     const double *x_p, const double *w_interp);
double Lagrange_array_interpolate(tNode *node, tArray *var, double Xb[3]);
double Lagrange_array_interpolate2d(tNode *node, tArray *var, int dir, int p,
                                    double Cb[2]);
void fill_3arrays_with_nodepoints(tNode *node, tArray *Xp[3]);
void fill_2arrays_with_nodepoints(tNode *node, int dir, tArray *Cp[2]);
void Lagrange_interpolate_topoints(tNode *node, tArray *var,
                                   tArray *Xp[3], tArray *interp);
void Lagrange_interpolate_toIpoints(tNode *node, tArray *var,
                                    tArray *Xp[3], tArray *Ip, tArray *interp);
void Lagrange_interpolate2d_topoints(tNode *node, tArray *var, int dir, int p,
                                     tArray *Cp[2], tArray *interp);
void Lagrange_interpolate2d_toIpoints(tNode *node, tArray *var, int dir,int p,
                                      tArray *Cp[2], tArray *Ip,
                                      tArray *interp);
void insert_array_inplane(tArray *var, int dir, int p, tArray *interp2d);

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
