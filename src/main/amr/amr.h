/* MemoryMan.h */
/* Wolfgang Tichy, 1/2019 */


/* world comm from main */
extern nMPI_Comm main_comm;


/* storage.c */
void realloc_dat_reqs(tDat *dat, int n_rq_new, int f);
int realloc_myln_nncats(tMylnodes *myln, int nncats);
int addto_myln_ln_c(tMylnodes *myln, int c, tNlist *elem);

/* mesh.c */
int setup_mesh(tMesh *mesh);
int setup_l2_mesh(tMesh *mesh);
int setup_3patchl2_mesh(tMesh *mesh);
int setup_test_mesh(tMesh *mesh);

/* connect.c */
void connect8_with_neighbors(tNode *narray[8], int connect);
void connect8_siblings(tNode *narray[8]);
tNlist *all_descendants_along_face(tNlist *nl, int face, int *ndescends);
tNlist *make_patch_neighbor_list(tNode *node, int face);
tNlist *make_mesh_neighbor_list(tNode *node, int face);
void update_node_fnb(tNode *node);
void update_node_and_neighbors_fnb(tNode *node);
int locate_facenb_in_fnbs(tNode *node, tNode *facenb, int *face, int *ni);

/* surface.c */
void free_surface(tSurface *s);
tSurface *init_surface(tNode *node, int vi, int face);
void set_mysurf(tSurface *s);
void free_dat_reqs_after_Waitall_com_send(tNode *node);
void set_ajsurf_forall_vars(tNode *node, int f);
void free_nbsurf_only_forall_vars(tNode *node, int f);

/* load.c */
void simple_load_balance(tMesh *mesh);
void move_node_to_rank(tNode *node, int desrank,
                       tCom *scom, tCom *rcom, int setbufs);

/* bfaces.c */
int Coordinates_set_bfaces_WT(tMesh *mesh);
int set_bfaces_on_patface(tPat *pat, int f);
