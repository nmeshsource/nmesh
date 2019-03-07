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
