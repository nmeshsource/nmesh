/* nmesh_checkpoint.h */
/* Wolfgang Tichy, 8/2019 */


/* checkpoint.c */
int checkpoint_load(tMesh *mesh);
int checkpoint_save_if_needed(tMesh *mesh, int always);
