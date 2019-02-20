/* nmesh_SurfExchange.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global functions */



/* coordinates.c */
int coordinates_init(tMesh *mesh);
int coordinates_init_node(tNode *node);

/* get_coords.c */
void nearest_ijk_of_XbYbZb(tNode *node, int ijk[3], const double Xb0[3]);
void XbYbZb_of_ijk(tNode *node, int i, int j, int k, double Xb[3]);
void XbYbZb_of_ind(tNode *node, int ind, double Xb[3]);
void XYZ_of_XbYbZb(tNode *node, const double Xb[3], double X[3]);
void XbYbZb_of_XYZ(tNode *node, double Xb[3], const double X[3]);
