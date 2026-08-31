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
int rtbrent_brak_fdf(double *x0,
                     void (*fdf)(double x, void *par, double *f, double *df),
                     double x1, double x2, void *par, int maxits, double xacc,
                     int pr);

/* widen_brak.c */
int widen_brak(double (*func)(double,void *par),
               double *x1, double *x2, void *par, int ntries, int pr);

/* newton1d_fd.c */
int newton1d_fd_region(double *x0, double (*func)(double x, void *par),
                       double x1, double x2, void *par,
                       int maxits, double xacc, int pr);
int find_2roots_region(double x0[2],
                       double (*func)(double x, void *par),
                       double x1, double x2, void *par,
                       int maxits, double xacc, int pr);

/* matrix.c */
int gaussjordan_inv(int n, double a[]);
int M_to_Minv_gaussjordan(int n, const double M[], double Minv[]);
void M_to_Mtranspose(int n, int m, const double M[], double MT[]);

/* binarysearch.c */
const void *binarysearch(const void *key, const void *base0,
                         size_t *base0offset, size_t *num, size_t size,
                         int (*compar)(const void *, const void *, void *),
                         void *arg);
void binarysearchmore(const void *key, const void *base0,
                      size_t nmemb, size_t size, const void *result,
                      int (*compar)(const void *, const void *, void *),
                      void *arg,
                      unsigned *more, unsigned *atBoundary);
const void *bisectionsearch(const void *key, const void *base0,
                            size_t *base0offset, size_t *num, size_t size,
                            int (*compar)(const void *, const void *, void *),
                            void *arg);
