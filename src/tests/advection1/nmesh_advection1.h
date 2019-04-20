/* nmesh_advection1.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global functions */

void advection1_fluxes_pt(tNode *node, int face, int i, int j, int k,
                          tVarList *vlu,
                          double *ui, double *ua,
                          double *fi,  double *fa,
                          double *lami, double *lama);
