/* nmesh_basis.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for global functions */


/* basis.c */

/* Legendre.c */
double basis_LegendreP(int l, double x, int np);
double basis_normLegendreP(int i, double x, int np);
void LGL_x_wquad(int npoints, double *x, double *w);
void LGL_AT_ST_matrices(int n, double *x, double *w, double *AT, double *ST);
void basis_array_analysis3(tNode *node, tArray *u, tArray *c);
void basis_array_synthesis3(tNode *node, tArray *u, tArray *c);
double basis_array_interpolate(tNode *node, tArray *coef, double Xb[3]);;

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

/* get_coords.c */
void nearest_ijk_of_XYZ(tNode *node, int ijk[3], const double X0[3]);
