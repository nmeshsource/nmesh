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
int xyz_of_lamAB_CubSph(tPat *pat, tNode *node, int ind,
                        const double lamAB[3], double xyz[3]);
int lamAB_of_xyz_CubSph(tPat *pat, tNode *node, int ind,
                        double lamAB[3], const double xyz[3]);
int dlamAB_dxyz_CubSph(tPat *pat, tNode *node, int ind, const double lamAB[3],
                       double xyz[3], double dlamABdxyz[3][3]);
int xyz_of_rhoAB_CubSph(tPat *pat, tNode *node, int ind,
                        const double rhoAB[3], double xyz[3]);
int rhoAB_of_xyz_CubSph(tPat *pat, tNode *node, int ind,
                        double rhoAB[3], const double xyz[3]);
int drhoAB_dxyz_CubSph(tPat *pat, tNode *node, int ind, const double rhoAB[3],
                       double xyz[3], double drhoABdxyz[3][3]);
int xyz_of_rh2AB_CubSph(tPat *pat, tNode *node, int ind,
                        const double rh2AB[3], double xyz[3]);
int rh2AB_of_xyz_CubSph(tPat *pat, tNode *node, int ind,
                        double rh2AB[3], const double xyz[3]);
int drh2AB_dxyz_CubSph(tPat *pat, tNode *node, int ind, const double rh2AB[3],
                       double xyz[3], double drh2ABdxyz[3][3]);

/* setup_CubedSpheres.c */
int add_6CubedSphere_pats(tMesh *mesh, int type, int stretch, int r_is_const,
                          double *xc, double *Din, double *Dout);
int add_N_CubedSphere_pats(tMesh *mesh, int N,
                           int type, int stretch, int r_is_const,
                           double *xc, double *Din, double *Dout);
void set_AB_min_max_from_Din(int dom, double *Din,
                             double *Amin, double *Amax,
                             double *Bmin, double *Bmax);
int add_1_CubedSphere_pat(tMesh *mesh, int dom, int type,
                          int stretch, int r_is_const, double *xc,
                          double Din, double Dout, double ABrct[4]);

/* setup_Boxes.c */
int add_1box_pat(tMesh *mesh, double xc[3], double dout[3]);
int add_1cube_pat(tMesh *mesh, double *xc, double dout);
