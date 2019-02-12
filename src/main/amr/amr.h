/* MemoryMan.h */
/* Wolfgang Tichy, 1/2019 */


/* storage.c */

/* mesh.c */
int setup_test_mesh(tMesh *mesh);

/* connect.c */
void connect8_with_neighbors(tNode *narray[8], int connect);
void connect8_siblings(tNode *narray[8]);
tNlist *all_descendants_along_face(tNlist *nl, int face, int *ndescends);
tNlist *make_patch_neighbor_list(tNode *node, int face);
tNlist *make_mesh_neighbor_list(tNode *node, int face);
void update_node_fnb(tNode *node);
void update_node_and_neighbors_fnb(tNode *node);
