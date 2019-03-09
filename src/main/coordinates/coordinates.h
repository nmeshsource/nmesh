/* coordinates.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coordinates local functions */

#define LDOUBLE double


/* coordinates.c */


/* coordtrans_CubedSphere.c */
double CubedSphere_sigma(tPat *pat, tNode *node, int si, int ind,
                         double A, double B);
void CubedSphere_dsigma_dAB(tPat *pat, tNode *node, int si, int ind,
                            double A, double B, double dSig[2]);
int ThetaPhi_of_AB_CubSph(tPat *pat, double A, double B,
                          double *Theta, double *Phi);
int ThetaPhi_dThetaPhidAB_of_AB_CubSph(tPat *pat, double A, double B,
                                       double *Theta,    double *Phi,
                                       double *dThetadA, double *dThetadB,
                                       double *dPhidA,   double *dPhidB);

/* setup_CubedSpheres.c */
int add_1cube_pat(tMesh *mesh, double *xc, double dout);
int add_6CubedSphere_pats(tMesh *mesh, int type, int stretch, int r_is_const,
                          double *xc, double *Din, double *Dout);
void set_AB_min_max_from_Din(int dom, double *Din,
                             double *Amin, double *Amax,
                             double *Bmin, double *Bmax);
