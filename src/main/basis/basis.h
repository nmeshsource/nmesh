/* basis.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for basis local functions */


/* basis.c */
int basis_init_globals(tMesh *mesh);
void Lagrange_InterpMatT(tArray *Xb, tArray *WL, tArray *Yb, tArray *Mt);
void Inverse_InterpMatT(tArray *Pt, tArray *wt, tArray *rt, tArray *Rt);
int IndexRange_Xb0_get(tNode *node, int dir, double Xb0, int n,
                       int CenterOnXb0, int *i0, int *ni);

/* gridpoints.c */
int gridpoints_init(tMesh *mesh);
int gridpoints_alloc(tMesh *mesh);
int gridpoints_free(tMesh *mesh);

/* Fourier.c */
void Fourier_coeffs(int N, const double u[], double c[]);
void set_TrafoArray(tArray *At,
                    void (*get_coeffs)(int N, const double *u, double *c));

/* Lagrange.c */
double Lagrange_interp_barycentric2_ds(double x, int n, const double *x_p,
                                       const double *w_interp,
                                       const double *f, int ds, double fscal);
void Lagrange_InterpMatrixT(int nx, const double *x, const double *w_interp,
                            int ny, const double *y, double *MT);

/* WENOinterp.c */
double interpolate_WENO_n_ds(double x, int n, const double *x_p,
                             const double *w_interp,
                             const double *f, int ds, double fscal);

/* Chebyshev.c */
void ChebyshevExtrema_DT(int np, double *DT);
void ChebyshevExtrema_AT_ST(int np, double *AT, double *ST);
void ChebyshevExtrema_x(int np, double *x);
double Chebyshev_basisfunc(int n, double X, int np);
