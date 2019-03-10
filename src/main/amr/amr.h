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
int setup_CubedSphere_mesh(tMesh *mesh);
int setup_l2_mesh(tMesh *mesh);
int setup_3patchl2_mesh(tMesh *mesh);
int setup_test_mesh(tMesh *mesh);

/* connect.c */
void connect8_with_neighbors(tNode *narray[8], int connect);
void connect8_siblings(tNode *narray[8]);
tNlist *ldescendants_along_face(tNlist *nl, int face, int *ndescends);
tNlist *make_patch_neighbor_list(tNode *node, int face);
tNlist *make_mesh_neighbor_list(tNode *node, int face);
void update_node_fnb(tNode *node);
void update_node_and_neighbors_fnb(tNode *node);
void update_all_rnode_fnb(tMesh *mesh);
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
void remove_all_bfaces(tPat *pat);
int amr_set_all_bfaces(tMesh *mesh);
int set_bfaces_on_patface(tPat *pat, int f);
void expand_bfaces_to_patch_edges(tMesh *mesh);
tBface *ith_bface_on_f_with_obface(tPat *pat, int f, int i);
int nbfaces_on_f_with_obface(tPat *pat, int f);
void mark_all_bfaces_without_op_as_outerbound(tMesh *mesh);
int set_consistent_flags_in_all_bfaces(tMesh *mesh);
int zero_face2_flag_in_all_bfaces(tMesh *mesh);
int toggle_face2_flag_in_faces4_5_of_cubes(tMesh *mesh);
int toggle_face2_flag_of_CubSph_doms_0_4_and_1_5(tMesh *mesh);
