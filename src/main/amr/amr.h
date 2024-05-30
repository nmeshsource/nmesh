/* MemoryMan.h */
/* Wolfgang Tichy, 1/2019 */



/* world comm from main */
extern nMPI_Comm main_comm;



/* storage.c */
int array_Neplocs(tArray *ar);
int redim_array_Neplocs(tArray *array, int Neplocs);
ulong alloc_and_set_mesh_myelm(tMesh *mesh);
void make_and_add_root_elm(tPat *pat, int n[3], int pt_typ[3], int datrank);
tElm *replace_8localchildren_by_parent(tElm *child0, int n[3], int pt_typ[3],
                                       struct list_head *ch_head);

/* mesh.c */
int amr_set_use_fv_flag(tMesh *mesh);
void remove_all_patches(tMesh *mesh);
int amr_set_bfaces_and_rnode_nbinfo_fnb(tMesh *mesh, int pr);
int setup_box_mesh(tMesh *mesh, double boxes_rc[3]);
int setup_CubedSphere_mesh(tMesh *mesh, double BoxMesh_rc[3],
                           double *CubedSphere_r);
int setup_Shell_mesh(tMesh *mesh, double CubedSphere_r);
int setup_elm_mesh1(tMesh *mesh);
int setup_l2_mesh(tMesh *mesh);
int setup_3patchl2_mesh(tMesh *mesh);
int setup_test_mesh(tMesh *mesh);

/* test_mesh.c */
int test_mesh(tMesh *mesh);

/* setup_Boxes.c */
int add_1box_pat(tMesh *mesh, double xc[3], double dout[3]);


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
void load_exchange_dat_after_moving_elms(tMesh *mesh);
void load_balance_elms(tMesh *mesh);

/* bfaces.c */
void remove_all_bfaces(tPat *pat);
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
void remove_all_nfaces(tNode *node);

/* refine.c */
void hp_refine_elms_if_rflag(tMesh *mesh, tRef *ref);
void set_children_nbinfo_ef(tElm *child0, tElm *parent,
                            khash_t(u32_gptr) *ef);
void set_parent_nbinfo_ef(tElm *parent, struct list_head *ch_head,
                          khash_t(u32_gptr) *ef);
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
void printTiming(void);
double timing_load(void);
int timing_print_load(tMesh *mesh);
int timing_mm_speed(tMesh *mesh);
double timing_get_mm_speed(tMesh *mesh);
double timing_get_elm_load_TimeIn_s(tElm *elm);
int timing_set_myops_ops0_allops(tMesh *mesh);
int timing_set_maxops(tMesh *mesh);

/* connections.c */
int connections_loc_on_patchface(int l, const char loc[NLOCS],
                                 int patface[6]);
void test_eploc(void);
int eploc1_eploc2_agree_upto_l_max(const tEploc *eploc1,
                                   const tEploc *eploc2, int l_max);
void amr_init_elm0_from_eploc(tMesh* mesh, tEploc *eploc, tElm0 *elm0);
int amr_get_8elms_at_elm_start(tElm *elm_start, void *ptr_elm);
int amr_get_8elms_at_myid(tMesh *mesh, ulong myid, void *ptr_elm);
int amr_elms_are_siblings(int n, void *ptr_elm);
int locate_facenb_in_fnbs(tNode *node, tNode *facenb, int *face, int *ni);
int amr_set_child_eloc(tEloc *parentloc, int ijk, tEloc *eloc);
int amr_set_child_eploc(tEploc *parenteploc, int ijk, tEploc *eploc);
int amr_set_parent_eploc(tEploc *eploc, tEploc *parenteploc);
tElm **amr_elmarray_bsearch_eloc(ulong narr, tElm **arr, tEloc *eloc);
void amr_set_intersibling_nbinfo_nnbinfo(tElm *sib0);
void amr_elm_nbinfo_add_nbeploc(tElm *elm, int face,
                                int nnb, const tEploc nbeploc[nnb]);
void amr_elm_nbinfo_redim_according_to_nnbinfo(tElm *elm);
void amr_elm_nbinfo_update_eid_locally_using_fnb_mesh(tMesh *mesh);
int amr_make_fnb_list(tElm *elm, int elmface, long narr, tElm **arr,
                      struct list_head *fnb_head);
void amr_erase_all_elm_fnb(tMesh *mesh);
int amr_get_elm0_for_eids(tMesh *mesh, ulong neids, ulong *eidarr,
                          tElm0 *elm0);
int amr_invalidate_nbinfo_of_all_nbs(tElm *elm, int Keep_nbs_fnb);
void amr_remove_mesh_nbelm(tMesh *mesh, int Keep_nbs_fnb);

int amr_khset_add_nb_ranks(tMesh *mesh, khash_t(u32) *nbranks);
void amr_khmap_add_elm_face_for_rank(khash_t(u32_gptr) *ef, int rank,
                                     tElm *elm, int face);
void amr_khmap_free_all_lists(khash_t(u32_gptr) *ef);
void amr_khmap_add_negelm_forallfaces(khash_t(u32_gptr) *ef, tElm *elm);
void amr_khmap_add_negchildren_forallparentfaces(khash_t(u32_gptr) *ef,
                                                 tElm *child0, tElm *parent);
void amr_khmap_add_negparent_forallchildrenfaces(khash_t(u32_gptr) *ef,
                                                 tElm *parent,
                                                 struct list_head *ch_head);
void amr_remove_mesh_nbelm_ef(tMesh *mesh, int Keep_nbs_fnb,
                              khash_t(u32_gptr) *ef);
