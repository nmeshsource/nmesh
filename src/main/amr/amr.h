/* MemoryMan.h */
/* Wolfgang Tichy, 1/2019 */


/* storage.c */

/* mesh.c */
int setup_test_mesh(tMesh *mesh);

/* connect.c */
void connect8_with_neighbors(tNode *narray[8], int connect);
void connect8_siblings(tNode *narray[8]);
tNlist *find_patch_neighbors(tNode *node, int face);
