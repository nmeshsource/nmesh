/* checkpoint.h */
/* Wolfgang Tichy, 8/2019 */


/* checkpoint.c */
void checkpoint_set_CI_Fcoef_filename(tPat *pat, int f, const char *pats,
                                      char *Fcoef_filename);
void checkpoint_set_nbinfo_fnb_nbelm_loadbal(tMesh *mesh, int reset_nbinfo);
void checkpoint_keep_previous(char *outdir, int pl, char *dirp, int nprev);

/* checkpoint_save.c */
int checkpoint_save_pars(tMesh *mesh, char *fname);
int checkpoint_save_patches(tMesh *mesh, char *fname);
void checkpoint_write_pat(FILE *fp, tPat *pat, char *fname);
void checkpoint_write_CI(FILE *fp, tPat *pat, char *fname);
void checkpoint_write_CI_Fcoef(tPat *pat, char *fname);
int checkpoint_save_elms(tMesh *mesh, char *fname);
int checkpoint_write_elms(tMesh *mesh, FILE *fp);
int checkpoint_save_EvoVars(tMesh *mesh, char *fname);
int checkpoint_save_nbinfoVars(tMesh *mesh, char *fname);
int checkpoint_save_VL(tMesh *mesh, char *fname, tVarList *vl,
                       int write_native);
void checkpoint_write_vl(FILE *fp, tVarList *vl, int write_native);
int checkpoint_save_CRCs(tMesh *mesh, char *fname);

/* checkpoint_load.c */
int checkpoint_load_patches(tMesh *mesh, char *fname);
void checkpoint_load_CI_Fcoef(tPat *pat, char *fname);
int checkpoint_load_elms(tMesh *mesh, char *fname);
int checkpoint_load_Vars(tMesh *mesh, char *fname, int read_native);
tVarList *checkpoint_make_vl(FILE *fp, tMesh *mesh);
void checkpoint_read_vl(tNode *node, char *buffer, long nbuffer,
                        tVarList *vl);
char *checkpoint_make_nodebuffer(FILE *fp, tVarList *vl, int read_native,
                                 long *nbuffer, char *name);
int checkpoint_load_CRCs(tMesh *mesh, char *fname);
