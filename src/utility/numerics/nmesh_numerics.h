/* nmesh_numerics.h */
/* Wolfgang Tichy, June 2019 */


/* newton1d_brak.c */
int newton1d_brak(double *x0,
                  void (*fdf)(double x, void *par, double *f, double *df),
                  double x1, double x2, void *par, int maxits, double xacc,
                  int pr);
int rtbisect(double *x0, double (*func)(double,void *par),
             double x1, double x2, void *par, int maxits, double xacc, int pr);

/* rtbrent_brak.c */
int rtbrent_brak(double *x0, double (*func)(double,void *par),
                 double x1, double x2, void *par, int maxits, double xacc,
                 int pr);

/* newton1d_fd.c */
int newton1d_fd_region(double *x0, double (*func)(double x, void *par),
                       double x1, double x2, void *par,
                       int maxits, double xacc, int pr);
int find_2roots_region(double x0[2],
                       double (*func)(double x, void *par),
                       double x1, double x2, void *par,
                       int maxits, double xacc, int pr);

/* matrix_inv.c */
int gaussjordan_inv(int n, double a[]);
int Minv_from_M_gaussjordan(int n, const double M[], double Minv[]);
