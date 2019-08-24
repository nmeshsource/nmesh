/* nmesh_checkpoint.h */
/* Wolfgang Tichy, 8/2019 */


/* checkpoint.c */
int checkpoint_exists(tMesh *mesh, const char *outdir_suffix,
                      const char *Dir_suffix);
int checkpoint_load_stage(tMesh *mesh, const char *outdir_suffix, int stage);
int checkpoint_save_if_needed(tMesh *mesh, int always);

/* checkpoint_load.c */
int checkpoint_load_pars(tMesh *mesh, char *fname);
