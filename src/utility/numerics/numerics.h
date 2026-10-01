/* numerics.h */
/* Wolfgang Tichy, June 2019 */

/* from main/main/utilities.c */
int finit(double x);


/* widen_brak.c */
int widen_brak_f(double (*func)(double,void *par),
                 double *x1, double *x2, void *par, double Fac, int ntries,
                 int pr);
