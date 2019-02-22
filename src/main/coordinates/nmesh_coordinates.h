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
void X_of_Xb_indir(tNode *node, int dir, double Xb, double *X);
void array_XYZ_of_XbYbZb(tNode *node, tArray *aXb[3], tArray *aX[3]);
void XYZ_of_ijk(tNode *node, int i, int j, int k, double X[3]);
void XbYbZb_of_XYZ(tNode *node, double Xb[3], const double X[3]);
void Xb_of_X_indir(tNode *node, int dir, double *Xb, const double X);
void array_XbYbZb_of_XYZ(tNode *node, tArray *aXb[3], tArray *aX[3]);
