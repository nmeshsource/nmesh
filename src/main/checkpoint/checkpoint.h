/* checkpoint.h */
/* Wolfgang Tichy, 8/2019 */


/* checkpoint_save.c */
int checkpoint_save_pars(tMesh *mesh, char *fname);
int checkpoint_save_patches(tMesh *mesh, char *fname);
void checkpoint_write_pat(FILE *fp, tPat *pat);
void checkpoint_write_CI(FILE *fp, tCoordInfo *CI);
int checkpoint_save_elms(tMesh *mesh, char *fname);
int checkpoint_write_elms(tMesh *mesh, FILE *fp);
int checkpoint_save_EvoVars(tMesh *mesh, char *fname);
void checkpoint_write_vl(FILE *fp, tVarList *vl, int write_big);

/* checkpoint_load.c */
int checkpoint_load_patches(tMesh *mesh, char *fname);
int checkpoint_load_elms(tMesh *mesh, char *fname);
int checkpoint_load_Vars(tMesh *mesh, char *fname);
tVarList *checkpoint_make_vl(FILE *fp, tMesh *mesh);
void checkpoint_read_vl(tNode *node, char *buffer, long nbuffer,
                        tVarList *vl);
char *checkpoint_make_nodebuffer(FILE *fp, tVarList *vl, int read_big,
                                 long *nbuffer, char *name);
