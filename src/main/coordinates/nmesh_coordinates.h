/* nmesh_SurfExchange.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global functions */



/* For coordtrans_CubedSphere.c :
   Type of cubed sphere or rather sphered cube coord transform */
enum
{
  CoordInfoNotSet,  /* if box->CI->type is not set, box->CI->type=0 */
  PyramidFrustum,   /* both inner & outer surfaces are flat */  
  innerCubedSphere, /* inner surface is curved, but outer surface is flat */
  outerCubedSphere, /* outer surface is curved, but inner surface is flat */
  CubedShell        /* both inner & outer surfaces are curved */
};


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
int set_xyz(tPat *pat, tNode *node, int ind, const double X[3], double x[3]);
int set_XYZ(tPat *pat, tNode *node, int ind, double X[3], const double x[3]);
int set_xyz_dXYZdxyz(tPat *pat, tNode *node, int ind,
                     const double X[3], double x[3], double dXYZdxyz[3][3]);
void brct_nodeface(tNode *node, int norm, double brct[4]);
void expand_brct_to_include_X(double brct[4], int norm,
                              const double X[3], int expand);
int intersection_brct1_brct2(const double brct1[4], const double brct2[4],
                             double brct[4]);
void C_from_X_on_face(const double X[3], int face, double C[2]);
void X_from_C_on_face(tPat *pat, int face, const double C[2], double X[3]);
int brctpat2_of_brctpat1(tPat *pat1, int f1, const double brct1[4],
                         tPat *pat2, int f2, double brct2[4]);
int C_in_brct(const double brct[4], double C[2]);
int fnb_containing_point(tNode *node, int f,
                         tPat *o_pat, int o_f, double C[2]);
void mark_points_in_fnb_f_ni(tNode *node, int f, int ni, tArray *aC[2],
                             tArray *aoC[2], tArray *aI);
void mark_points_in_nb_f(tNode *node, int f, tArray *aC[2],
                         tNode *nb, int nb_f, tArray *aoC[2], tArray *aI);
void array_nbXface_of_Xface(tNode *node, int f,
                            tNode *nb, int nb_f, tArray *nbC[2]);
void array_find_nbXface_of_Xface(tNode *node, int f, tNode *nb, int nb_f,
                                 tArray *nbC[2], tArray *nbI);

/* derivs.c */
int cart_partials(tNode *node, int ui, int dui);

/* utilities.c */
double patch_normal_at_XYZ(tPat *pat, int f, const double X[3], double n[3]);
double node_normal_at_ijk(tNode *node, int f, int ijk, double n[3]);
double det_3Dmatrix(double M[3][3]);
double inv3Dmat_from_3Dmat(double M[3][3], double invM[3][3]);
double smallest_pat_size(tMesh *mesh);

/* ComplexFunctions.c */
double BaseAngle(double p, double peri, double p0);
double Arg(double x, double y);
double Arg_plus(double x, double y);
double dArgdx(double x, double y);
double dArgdy(double x, double y);
double ddArgdxdx(double x, double y);
double ddArgdxdy(double x, double y);
double ddArgdydy(double x, double y);

/* setup_CubedSpheres.c */
int arrange_12CubSph_into_empty_cube(tMesh *mesh, double *xc,
                                     double din, double dmid, double dout);
int arrange_1pat12CubSph_into_full_cube(tMesh *mesh, double *xc,
                                        double din, double dmid, double dout);
int two_full_cubes_touching_at_x0(tMesh *mesh, double dc,
                                  double din1, double dmid1,
                                  double din2, double dmid2);
int sphere_around_two_full_cubes_touching_at_x0(tMesh *mesh,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0);
int two_spheres_around_two_full_cubes(tMesh *mesh,
        double dc, double din1, double dmid1, double din2, double dmid2,
        double r0, double r1);
int two_spheres_around_box_at_x0(tMesh *mesh, double dc[3],
                                 double r0, double r1);
int two_spheres_around_empty_box_at_x0(tMesh *mesh, double dc[3],
                                       double r0, double r1);
int two_wegdes_touching_1_wedge(tMesh *mesh, double dc, double r0, double r1);
int two_diff_wegdes_touching_1_wedge(tMesh *mesh, double dc,
                                     double r0, double r1);

/* coordtrans_CubedSphere.c */
int xyz_of_lamAB_CubSph(tPat *pat, tNode *node, int ind,
                        const double lamAB[3], double xyz[3]);
int lamAB_of_xyz_CubSph(tPat *pat, tNode *node, int ind,
                        double lamAB[3], const double xyz[3]);
int dlamAB_dxyz_CubSph(tPat *pat, tNode *node, int ind, const double lamAB[3],
                       double xyz[3], double dlamABdxyz[3][3]);
int xyz_of_rhoAB_CubSph(tPat *pat, tNode *node, int ind,
                        const double rhoAB[3], double xyz[3]);
int rhoAB_of_xyz_CubSph(tPat *pat, tNode *node, int ind,
                        double rhoAB[3], const double xyz[3]);
int drhoAB_dxyz_CubSph(tPat *pat, tNode *node, int ind, const double rhoAB[3],
                       double xyz[3], double drhoABdxyz[3][3]);
