/* basis.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for basis local functions */


/* basis.c */
int basis_init_globals(tMesh *mesh);

/* gridpoints.c */
int init_gridpoints(tMesh *mesh);
int free_gridpoints(tMesh *mesh);

/* Fourier.c */
void Fourier_coeffs(int N, const double u[], double c[]);
void set_TrafoArray(tArray *At,
                    void (*get_coeffs)(int N, const double *u, double *c));

/* Lagrange.c */
double Lagrange_interp_barycentric2_ds(double x, int n, const double *x_p,
                                       const double *w_interp,
                                       const double *f, int ds, double fscal);

/* WENOinterp.c */
double interpolate_WENO_n_ds(double x, int n, const double *x_p,
                             const double *w_interp,
                             const double *f, int ds, double fscal);
