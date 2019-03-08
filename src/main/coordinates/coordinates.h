/* coordinates.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coordinates local functions */

#define LDOUBLE double


/* coordinates.c */


/* coordtrans_CubedSphere.c */
double CubedSphere_sigma(tPat *pat, tNode *node, int si, int ind,
                         double A, double B);
double CubedSphere_dsigma_dA(tPat *pat, tNode *node, int si, int ind,
                             double A, double B);
double CubedSphere_dsigma_dB(tPat *pat, tNode *node, int si, int ind,
                             double A, double B);

/* setup_CubedSpheres.c */
int convert_1pat_to_cube(tMesh *mesh, int b0, double *xc, double dout);
int convert_6pats_to_CubedSphere(tMesh *mesh, int p0, int type, int stretch,
                                  double *xc, double *Din, double *Dout);
void set_AB_min_max_from_Din(int dom, double *Din,
                             double *Amin, double *Amax,
                             double *Bmin, double *Bmax);
