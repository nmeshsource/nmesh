/* nmesh_dg.h */
/* Wolfgang Tichy, April 2019 */


/* numflux1d.c */
void numflux1d_LLF(tMesh *mesh, int nf, double *fnum,
                   double *uL, double *uR, double *fL, double *fR,
                   double *lamL, double *lamR);

/* dg.c */
int dg_add_surface_fluxes(tMesh *mesh, tVarList *vlr, tVarList *vlu,
                          void (*u_f_lam)(tNode *node, int face,
                                          int i, int j, int k,
                                          tVarList *vlu,
                                          double *ui, double *ua,
                                          double *fi,  double *fa,
                                          double *lami, double *lama),
                          void (*numflux)(tMesh *mesh, int nf, double *fnum,
                                          double *uL, double *uR,
                                          double *fL, double *fR,
                                          double *lamL, double *lamR));
