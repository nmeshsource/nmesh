/* nmesh_dg.h */
/* Wolfgang Tichy, April 2019 */


/* numflux1d.c */
void numflux1d_LLF(tMesh *mesh, int nf, double *fnum,
                   double *uL, double *uR, double *fL, double *fR,
                   double *lamL, double *lamR);

/* dg.c */
int dg_add_surface_fluxes(tMesh *mesh, tVarList *vlr, tVarList *vlu);
