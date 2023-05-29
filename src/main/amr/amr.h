/* MemoryMan.h */
/* Wolfgang Tichy, 1/2019 */



/* world comm from main */
extern nMPI_Comm main_comm;



/* storage.c */
tNode *destroy_children(tNode *parent);
tNlist *remove8siblings_in_mesh_lns(tNlist *sib);
int realloc_myln_nncats(tMylnodes *myln, int nncats);
int addto_myln_ln_c(tMylnodes *myln, int c, tNlist *elem);

/* mesh.c */
int amr_set_use_fv_flag(tMesh *mesh);
void remove_all_patches(tMesh *mesh);
int setup_box_mesh(tMesh *mesh);
int setup_CubedSphere_mesh(tMesh *mesh);
int setup_Shell_mesh(tMesh *mesh);
int setup_l2_mesh(tMesh *mesh);
int setup_3patchl2_mesh(tMesh *mesh);
int setup_test_mesh(tMesh *mesh);

/* test_mesh.c */
int test_mesh(tMesh *mesh);

/* setup_Boxes.c */
int add_1box_pat(tMesh *mesh, double xc[3], double dout[3]);

/* connect.c */
void connect8_with_neighbors(tNode *narray[8], int connect);
void connect8_siblings(tNode *narray[8]);
tNlist *ldescendants_along_face(tNlist *nl, int face, int *ndescends);
void update_node_and_neighbors_nfaces_fnb(tNode *node);
void update_all_rnode_nfaces_fnb(tMesh *mesh);
tNlist *make_patch_neighbor_list(tNode *node, int face);
tNlist *make_mesh_neighbor_list(tNode *node, int face);
void update_node_fnb_only(tNode *node);
int locate_facenb_in_fnbs(tNode *node, tNode *facenb, int *face, int *ni);
void node_and_fnbs_lock(tNode *node);
void node_and_fnbs_unlock(tNode *node);
void parent_and_fnbs_lock(tNode *narray[8], tNode *locker);
void parent_and_fnbs_unlock(tNode *narray[8], tNode *locker);
tNode *get_node_nc_lock(tNode *node);

/* surface.c */
void free_surface(tSurface *s);
tSurface *init_surface(tNode *node, int vi, int face);
void set_mysurf(tSurface *s);
void free_dat_reqs_after_Waitall_com_send(tNode *node);
void set_ajsurf_forall_vars(tNode *node, int f);
void free_nbsurf_only_forall_vars(tNode *node, int f);
void copy_ajsurf_from_nbsurf0(tNode *node, int f, int nb_f,
                              int intrch, int rev1, int rev2);

/* load.c */
void simple_load_balance(tMesh *mesh);
void move_node_to_rank(tNode *node, int desrank,
                       tCom *scom, tCom *rcom, int setbufs);

/* bfaces.c */
void remove_all_bfaces(tPat *pat);
int amr_set_all_bfaces(tMesh *mesh);
void find_external_faces_of_pat(tPat *pat, double Lmin,
                                int *extface, int inclOuterBound);
int set_bfaces_on_patface(tPat *pat, double Lmin, int f);
void expand_bfaces_to_patch_edges(tMesh *mesh);
tBface *ith_bface_on_f_with_obface(tPat *pat, int f, int i);
int nbfaces_on_f_with_obface(tPat *pat, int f);
void mark_all_bfaces_without_op_as_OUTERBOUND(tMesh *mesh);
int set_consistent_flags_in_all_bfaces(tMesh *mesh);
int zero_face2_flag_in_all_bfaces(tMesh *mesh);
int toggle_face2_flag_in_faces4_5_of_cubes(tMesh *mesh);
int toggle_face2_flag_of_CubSph_doms_0_4_and_1_5(tMesh *mesh);
int facepoint_in_bfacepair(tBface *bface, tNode *node, int ijk, double C[2],
                           int ofaces[6], double oX[3]);
int add_nface(tNode *node, int f, tNode *nb, int nb_f);
int remove_nface(tNface *nface);
void remove_all_nfaces(tNode *node);

/* refine.c */
void hrefine_mesh_to_level(tMesh *mesh, int l);
void hrefine_mesh_to_level_loadbalance(tMesh *mesh, int l);
void hcoarsen_mesh_to_level(tMesh *mesh, int l);
void hrefine_pat(tMesh *mesh, int p);
void hcoarsen_pat(tMesh *mesh, int p);
void hrefine_sphere_loadbalance(tMesh *mesh, double radius, double xc[3],
                                int levels);

/* ghostzones.c */
void request_all_myln_ghostdata(tMesh *mesh);
void get_all_myln_ghostdata(tMesh *mesh);

/* amr.c */
int amr_print_thread_info(tMesh *mesh);

/* timing.c  */
int timing_mm_speed(tMesh *mesh);
