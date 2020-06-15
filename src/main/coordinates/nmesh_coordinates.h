/* nmesh_SurfExchange.h */
/* (c) Wolfgang Tichy 2/2019 */
/* header file for global functions */



/* labels for the different coord trafos that can be in pat->xyz_of_XYZ */
enum
{
  Cartesian=0,          /* Cartesian coords */
  CubedSphere,          /* cubed sphere coords xyz_of_lamAB_CubSph */
  stretchedCubedSphere  /* stretched cubed sphere coords xyz_of_rhoAB_CubSph */
};

/* Extras for coordtrans_CubedSphere.c : Type of cubed sphere.
   NOTE: Cartesian means box->CI->type=0, so we start with PyramidFrustum=1 */
enum
{
  /* Cartesian=0, */
  PyramidFrustum=1, /* both inner & outer surfaces are flat */
  innerCubedSphere, /* inner surface is curved, but outer surface is flat */
  outerCubedSphere, /* outer surface is curved, but inner surface is flat */
  CubedShell        /* both inner & outer surfaces are curved */
};


/* coordinates.c */
int coordinates_init(tMesh *mesh);
int coordinates_init_node(tNode *node);
int coordinates_get_label(tPat *pat);

/* get_coords.c */
void nearest_ijk_of_XbYbZb(tNode *node, int ijk[3], const double Xb0[3]);
void nearest_ijk_of_XYZ(tNode *node, int ijk[3], const double X0[3]);
void nearest_lowernode_ijk_of_XYZ(tNode *node, int ijk[3], const double X0[3]);
double nearest_ijk_of_xyz(tNode *node, int ijk[3], const double x0[3]);
double nearest_ijk_of_xyz_inplaneN(tNode *node, int N, int pl,
                                   int ijk[3], const double x0[3]);
double nearest_corner_of_xyz_inplaneN(tNode *node, int N, int pl,
                                      int ijk[3], const double x0[3]);
int approxXYZnormal_of_xyznormal(tNode *node, int cartN);
double magnitude_xyz(const double x[3]);
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
int p_XYZ_of_xyz_mesh(tMesh *mesh, double X[3], const double x[3]);
tNode *node_XYZ_of_xyz_mesh(tMesh *mesh, double X[3], const double x[3]);
long l_XYZ_of_xyz(tNode *node, int ind, double X[3], const double x[3]);
int XYZ_on_face(tPat *pat, int *face, const double X[3]);
int ind_on_nodeface(tNode *node, int ind, int *face);
int ind_on_outerbound(tNode *node, int ind);
int set_xyz(tPat *pat, tNode *node, int ind, const double X[3], double x[3]);
int set_XYZ(tPat *pat, tNode *node, int ind, double X[3], const double x[3]);
int set_xyz_dXYZdxyz(tPat *pat, tNode *node, int ind,
                     const double X[3], double x[3], double dXYZdxyz[3][3]);
void brct_nodeface(tNode *node, int norm, double brct[4]);
void resize_brct(double brct[4], double eps);
void expand_brct_to_include_X(double brct[4], int norm,
                              const double X[3], int expand);
int intersection_brct1_brct2(const double brct1[4], const double brct2[4],
                             double brct[4]);
int touch_or_intersect_bb1_bb2(const double bb1[6], const double bb2[6],
                               double bb[6]);
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
void array_find_nbXface_of_Xface(tNode *node, int f, tNode *nb, int nb_f,
                                 tArray *nbC[2], tArray *nbI);

/* derivs.c */
int cart_partials(tNode *node, int ui, int dui[3]);
void cart_partials_dU_di(tNode *node, int U, int dUx);
void cart_partials_dUi_dj(tNode *node, int Ux, int dUxx);
void cart_partials_dSij_dk(tNode *node, int Sxx, int dSxxx);
void cart_partials_dUij_dk(tNode *node, int Uxx, int dUxxx);
void cart_partials_diUj(tNode *node, int Ux, int dUxx);
void cart_partials_diSjk(tNode *node, int Sxx, int dSxxx);
void cart_partials_diUjk(tNode *node, int Uxx, int dUxxx);
void cart_partials_dUa_di(tNode *node, int Ut, int dUtx);
void cart_partials_dSab_di(tNode *node, int Stt, int dSttx);
void cart_partials_dUab_di(tNode *node, int Utt, int dUttx);
#define cart_partials_diU(node, U, dU) cart_partials_dU_di((node), (U), (dU))
void cart_partials_diUa(tNode *node, int Ut, int dUxt);
void cart_partials_diSab(tNode *node, int Stt, int dSxtt);
void cart_partials_diUab(tNode *node, int Utt, int dUxtt);
void cart_partials_diTjab(tNode *node, int Txtt, int dTxxtt);
void cart_partials_dU_di_dU_dij(tNode *node, int U, int dUx, int ddUxx);
void cart_partials_dUi_dj_dUi_djk(tNode *node, int Ux, int dUxx, int ddUxxx);
void cart_partials_dSij_dk_dSij_dkl(tNode *node, int Sxx,
                                    int dSxxx, int ddSxxxx);
int cart_di_Ui(tNode *node, int Ux, int divUi);

/* utilities.c */
double patch_normal_at_XYZ(tPat *pat, int f, const double X[3], double n[3]);
double node_normal_at_ijk(tNode *node, int f, int ijk, double n[3]);
double det_3Dmatrix(double M[3][3]);
double inv3Dmat_from_3Dmat(double M[3][3], double invM[3][3]);
double det_3Dsymmmat(double M11, double M12, double M13,
                     double M22, double M23, double M33);
double inv3Dmat_from_3Dsymmmat(double M11, double M12, double M13,
                               double M22, double M23, double M33,
                               double *i11, double *i12, double *i13,
                               double *i22, double *i23, double *i33);
double det_4Dsymmmat(double s11,double s12,double s13,double s14,
                     double s22, double s23, double s24,
                     double s33, double s34, double s44);
double inv4Dmat_from_4Dsymmmat(double s11,double s12,double s13,double s14,
                               double s22, double s23, double s24,
                               double s33, double s34, double s44,
                               double *i11,double *i12,double *i13,double *i14,
                               double *i22, double *i23, double *i24,
                               double *i33, double *i34, double *i44);
double invmat_from_symmmat(int n, const double *S, double *I);
int index_symmmat(int n, int a, int b);
double matel_symmmat(int n, const double *M, int a, int b);
void symmmat_times_vec(int n, const double *M, const double *v, double *Mv);
double vec_times_vec(int n, const double *u, const double *v);
double mag2_vector_metric(int n, const double *g, const double *v);
double BM_symmmat(int n, const double *B, const double *M);
void BMB_symmmat(int n, const double *B, const double *M, double *BMB);
void BBBMa_symm_bc(int n, const double *B, const double *M, double *BBBMa,
                   int a);
void BBBBMab_symm_ab_cd(int n, const double *B, const double *M,
                        double *BBBBMab, int a, int b);
double invertNxNx1symm_array(tArray *a);
double symmmat3D_times_vec(double gxx, double gxy, double gxz,
                           double gyy, double gyz, double gzz,
                           double Vx, double Vy, double Vz,
                           double *gVx, double *gVy, double *gVz);
double mag2_vector_3Dmetric(double gxx, double gxy, double gxz,
                            double gyy, double gyz, double gzz,
                            double Vx, double Vy, double Vz);
double smallest_pat_size(tMesh *mesh);
void array_dXdx(tNode *node, int ind, tArray *aJ);
void array_dxdX(tNode *node, int ind, tArray *aJ);
void array_2dxdX(tNode *node, int ind, int norm, tArray *aJ);
void array_2dxdXb(tNode *node, int ind, int norm, tArray *aJ);
void array_3Dmat(double M[3][3], tArray *aM);
double det_2_2_array(tArray *aM);
void var_to_var_times_JtoPower(tNode *node, int ui, int Jpower);
double find_hmin(tNode *node, int *ijk1, int *ijk2);

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
int sphere_around_empty_box_at_xc(tMesh *mesh, double xc[3], double dc[3],
                                  double r0);
int sphere_around_full_box_at_xc(tMesh *mesh, double xc[3], double dc[3],
                                 double r0);
int two_spheres_around_box_at_xc(tMesh *mesh, double xc[3], double dc[3],
                                 double r0, double r1);
int two_spheres_around_empty_box_at_xc(tMesh *mesh, double xc[3],
                                       double dc[3], double r0, double r1);
int CubedSphere_shell_at_xc(tMesh *mesh, double xc[3],
                            double rin, double rout);
int two_wegdes_touching_1_wedge(tMesh *mesh, double dc, double r0, double r1);
int two_diff_wegdes_touching_1_wedge(tMesh *mesh, double dc,
                                     double r0, double r1);
int add_1_CubedSphere_pat(tMesh *mesh, int dom, int type,
                          int stretch, int SigFunc, double *xc,
                          double Din, double Dout, double ABrct[4]);
int set_1_CubedSphere_pat(tPat *pat, int stretch, int SigFunc);

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

/* setup_Boxes.c */
int add_Nbox_pats_indir(tMesh *mesh, double xc[3], double dout[3],
                        int N, int dir);
int arrange_box_pats_inBox(tMesh *mesh, double xc[3], double dout[3], int N[3]);

/* pointlists.c */
intList *pointindexList_line(tNode *node, int dir, int axis[]);
intList *pointindexList_plane(tNode *node, int normal, int plane[]);
intList *pointindexList_node(tNode *node);

/* integrals.c */
double NodeVolumeIntegral(tNode *node, int vind, double power, int mode);
