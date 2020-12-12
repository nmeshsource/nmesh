/* coordinates.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for coordinates local functions */

#define LDOUBLE double


/* structure that holds global coordinates vars and pars */
typedef struct {
  int itmp1;        /* Ind("coordinates_tmp1") */
  int idet_dXbdx;   /* Ind("det_dXbdx") */
} tcoordinates;


/* coordinates.c */
void coordinates_set_sqrtgdiag_var(tNode *node, int idXdx, int igxx,
                                   int isqrtgdiagx);
void coordinates_set_sqrtgdiag_array(tNode *node, tArray *AdXdx[3][3],
                                     tArray *Ag[6], tArray *Asqrtgdiag[3]);

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
