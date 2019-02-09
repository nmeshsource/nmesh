/* nmesh_basis.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for global functions */


/* basis.c */
void LGL_x_winteg(int npoints, double *x, double *w);
void Lagrange_winterp(int n, const double *x, double *w_interp);
void LGL_DT(int n, const double *x, const double *w_interp, double *DT);
void LGL_AT_ST_matrices(int n, double *x, double *w, double *AT, double *ST);
