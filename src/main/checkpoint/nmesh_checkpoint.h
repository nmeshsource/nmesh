/* nmesh_checkpoint.h */
/* Wolfgang Tichy, 8/2019 */



/* parts that an existing checkpoint can have (or not) */
enum
{
  CHECKPOINT_PATS=2,
  CHECKPOINT_NBINFO=8,
  CHECKPOINT_VARS=16,
};


/* checkpoint.c */
int checkpoint_exists(tMesh *mesh, const char *outdir_suffix,
                      const char *Dir_suffix);
int checkpoint_load_stage(tMesh *mesh, const char *outdir_suffix, int stage);
int checkpoint_save_if_needed(tMesh *mesh, int always);
