/* dg.h */
/* Wolfgang Tichy, April 2019 */



/* for WENO.c: */
/* Structure for weights needed for WENO3 at one midpoint.
   All weights are ordered from left to right along the increasing
   point index:
   final rec. = optw[0] lin.rec0 + optw[1] lin.rec1
   where lin.rec0 = lw[0][0]*u_0 + lw[0][1]*u_1
         lin.rec1 = lw[1][0]*u_1 + lw[1][1]*u_2  */
typedef struct WENO3weight {
  double lw[2][2]; //normalized linear weights
  double optw[2];  //non-normalized optimal WENO3-weights
} tWENO3weight;

/* object to hold WENO weights for many points */
typedef struct WENOweights {
  int n;                  /* number of points for which we have weights */
  tWENO3weight **p_WENO3; /* array with weights for WENO3 in pos dir (p) */
  tWENO3weight **m_WENO3; /* array with weights for WENO3 in neg dir (m) */
} tWENOweights;

/* for WENO5 and WENOZ */
typedef struct WENO5weight {
  double qw[3][3]; //normalized quadratic weights
  double optw[3];  //non-normalized optimal WENO5-weights
} tWENO5weight;


/* funcs in WENO.c */
int WENOweights_init_globals(tMesh *mesh);
int WENOweights_free_globals(tMesh *mesh);
tWENO3weight *WENOweights_global_p_WENO3_at_(int i);
tWENO3weight *WENOweights_global_m_WENO3_at_(int i);
tWENO3weight *WENOweights_global_p_WENO3_at_last_minus_(int l);
tWENO3weight *WENOweights_global_m_WENO3_at_last_minus_(int l);
int pr_weight_ratios(tMesh *mesh);


/* funcs in dg.c */
int dg_set_DGglobals(tMesh *mesh);
int dg_free_DGglobals(tMesh *mesh);
int dg_print_DGglobals(tMesh *mesh);


/* funcs in tests.c */
int fv_tests(tMesh *mesh);


/* funcs in rec1d.c */
void rec1d_LR_uface_to_uin_1_Carray(int n, double *u, int right, int forward,
                                    double u_scale, double s1, double s2,
                                    int opt);

/* penalty.c */
double dg_scale_penalty_bamps(tDGinfo *dgi);
