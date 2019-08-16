/* nmesh_checkpoint.h */
/* Wolfgang Tichy, 8/2019 */


/* checkpoint.c */
int checkpoint_exists(tMesh *mesh, const char *outdir_suffix,
                      const char *Dir_suffix);
int checkpoint_load(tMesh *mesh, const char *outdir_suffix);
int checkpoint_save_if_needed(tMesh *mesh, int always);
