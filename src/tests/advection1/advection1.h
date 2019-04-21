/* advection1.h */
/* (c) Wolfgang Tichy 1/2019 */
/* header file for advection1 local functions */



/* advection1.c */
void advection1_fluxes_pt(tNode *node, int face, int i, int j, int k,
                          tVarList *vlu,
                          double *ui, double *ua,
                          double *fi,  double *fa,
                          double *lami, double *lama);
int advection1_init(tMesh *mesh);
int advection1_analyze(tMesh *mesh);
