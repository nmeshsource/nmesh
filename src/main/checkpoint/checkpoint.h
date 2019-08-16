/* checkpoint.h */
/* Wolfgang Tichy, 8/2019 */


/* checkpoint_save.c */
int checkpoint_save_pars(tMesh *mesh, char *fname);
int checkpoint_save_patches(tMesh *mesh, char *fname);
void checkpoint_write_pat(FILE *fp, tPat *pat);
void checkpoint_write_CI(FILE *fp, tCoordInfo *CI);
int checkpoint_save_nodes(tMesh *mesh, char *fname);
void checkpoint_write_nodetrees(FILE *fp, tNlist *rnlist);
void checkpoint_write_nodes_with_child0(FILE *fp, tNlist *nlist);
void checkpoint_write_node(FILE *fp, tNode *node);
int checkpoint_save_EvoVars(tMesh *mesh, char *fname);
void checkpoint_write_vl(FILE *fp, tVarList *vl, int write_big);
