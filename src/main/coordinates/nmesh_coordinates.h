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
void dXYZ_dXbYbZb(tNode *node, double dXdXb[3]);
void dXbYbZb_dXYZ(tNode *node, double dXbdX[3]);
void X_of_Xb_indir(tNode *node, int dir, double Xb, double *X);
void array_XYZ_of_XbYbZb(tNode *node, tArray *aXb[3], tArray *aX[3]);
void array_Xplane_of_Xb(tNode *node, int dir, tArray *aCb[2], tArray *aC[2]);
void XYZ_of_ijk(tNode *node, int i, int j, int k, double X[3]);
void XbYbZb_of_XYZ(tNode *node, double Xb[3], const double X[3]);
void Xb_of_X_indir(tNode *node, int dir, double *Xb, const double X);
void array_XbYbZb_of_XYZ(tNode *node, tArray *aXb[3], tArray *aX[3]);
void array_Xbplane_of_X(tNode *node, int dir, tArray *aCb[2], tArray *aC[2]);
int XYZ_is_in_node(tNode *node, double X[3]);
void array_find_XYZ_in_node(tNode *node, tArray *aXP[3], tArray *aI);
void array_find_Xplane_in_node(tNode *node,int dir, tArray *aCP[2], tArray *aI);
int p_XYZ_of_xyz(tPat *pat, double X[3], const double x[3]);
int p_XYZ_of_xyz_inpatlist(tMesh *mesh, intList *pl,
                           double X[3], const double x[3]);
int XYZ_on_face(tPat *pat, int *face, const double X[3]);

/* derivs.c */
int cart_partials(tNode *node, int ui, int dui);

/* utilities.c */
double patch_normal_at_XYZ(tPat *pat, int f, const double X[3], double n[3]);
double node_normal_at_ijk(tNode *node, int f, int ijk, double n[3]);
double det_3Dmatrix(const double M[3][3]);
double inv3Dmat_from_3Dmat(const double M[3][3], double invM[3][3]);
double smallest_pat_size(tMesh *mesh);
