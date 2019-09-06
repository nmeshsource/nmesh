/* filter.c */
/* Wolfgang Tichy, 8/2019 */

#include "nmesh.h"
#include "basis.h"


#define PR 0


/***********************************************************************/
/* functions to filter */
/***********************************************************************/

/* Apply exponential filter as in https://arxiv.org/abs/1804.02003
   In 1D: c_i -> c_i * e^{-alp (i/(n0-1))^s} */
void expfilter_coeff_array(tArray *ca, double alp[3], double s[3])
{
  int *n = ca->n;
  double N0 = n[0] - 1;
  double N1 = n[1] - 1;
  double N2 = n[2] - 1;
  double *c = Arrd_(ca);
  double alp0 = alp[0];
  double alp1 = alp[1];
  double alp2 = alp[2];
  double s0 = s[0];
  double s1 = s[1];
  double s2 = s[2];
  double f0, f1, f2;
  int i,j,k, ijk;

  /* make sure we don't divide by zero */
  if(N0==0.) N0 = 1.;
  if(N1==0.) N1 = 1.;
  if(N2==0.) N2 = 1.;

  /* loop over array and let c -> c * f0*f1*f2 */
  for(k = 0; k < n[2]; k++)
  {
    f2 = exp( -alp2 * pow(k/N2, s2) );
    for(j = 0; j < n[1]; j++)
    {
      f1 = exp( -alp1 * pow(j/N1, s1) );
      for(i = 0; i < n[0]; i++)
      {
        f0 = exp( -alp0 * pow(i/N0, s0) );
        ijk = Ind_n(i,j,k, n);
        c[ijk] = c[ijk] * f0 * f1 * f2;
      }
    }
  }
}

/* get coeffs ca of array ua, filter and then reset ua */
void expfilter_array(tNode *node, tArray *ua, double alp[3], double s[3])
{
  DECL_STACK_ARRAY(ca, ua->n);

  basis_array_analysis3(node, ua, ca);
  expfilter_coeff_array(ca, alp, s);
  basis_array_synthesis3(node, ua, ca);
}

/* filter var with index ui */
int expfilter_var(tNode *node, int ui, double alp[3], double s[3])
{
  tArray *ua;
  tDat *dat = node->dat;
  if(dat)
  {
    var_to_var_times_JtoPower(node, ui, 1);
    ua = dat->v[ui];
    expfilter_array(node, ua, alp, s);
    var_to_var_times_JtoPower(node, ui, -1);
    return 1;
  }
  else
    return 0;
}

/* filter varlist on entire mesh */
void expfilter_mesh_vl(tVarList *vl, double alp[3], double s[3])
{
  tMesh *mesh = vl->mesh;

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int vli;

    forvl(vl, vli)
    {
      expfilter_var(node, vli, alp, s);
    }
  }
}

/* filter varlist on entire mesh with same filter pars in all 3 directions */
void expfilter_vl(tVarList *vl, double af, double sf)
{
  double alp[] = { af, af, af }; /* use same filter pars in all 3 directions */
  double s[]   = { sf, sf, sf };
  expfilter_mesh_vl(vl, alp, s);
}

/***********************************************************************/
/* functions to check coeff falloff */
/***********************************************************************/

/* is array ca falling off? */
int has_expfalloff_coeff_array(tArray *ca, double alp[3], double s[3])
{
  int *n = ca->n;
  double N0 = n[0] - 1;
  double N1 = n[1] - 1;
  double N2 = n[2] - 1;
  double *c = Arrd_(ca);
  double absc;
  double alp0 = alp[0];
  double alp1 = alp[1];
  double alp2 = alp[2];
  double s0 = s[0];
  double s1 = s[1];
  double s2 = s[2];
  double f0, f1, f2;
  double norm;
  int i,j,k, ijk;

  /* make sure we don't divide by zero */
  if(N0==0.) N0 = 1.;
  if(N1==0.) N1 = 1.;
  if(N2==0.) N2 = 1.;

  /* L1 norm of coeffs */
  norm = Lp_norm_array(ca, 1.);

  /* loop over array and check if coeffs are falling off exponentially */
  for(k = 0; k < n[2]; k++)
  {
    f2 = exp( -alp2 * pow(k/N2, s2) );
    for(j = 0; j < n[1]; j++)
    {
      f1 = exp( -alp1 * pow(j/N1, s1) );
      for(i = 0; i < n[0]; i++)
      {
        f0 = exp( -alp0 * pow(i/N0, s0) );
        ijk = Ind_n(i,j,k, n);
        absc = fabs(c[ijk]);
        //printf("%d %d %d: absc=%g\n", i,j,k, absc);

        /* if absc is too big we do not have proper exp. falloff */
        if( (absc > norm * f0) || (absc > norm * f1) || (absc > norm * f2) )
        {
          if(PR)
          {
            printf("%d %d %d: ", i,j,k);
            printf("absc=%g norm=%g f0=%g f1=%g f2=%g\n", absc,norm, f0,f1,f2);
          }
          //abort();
          return 0;
        }
      } /* end for i */
    }
  }
  //if(absc!=0.) abort();
  return 1;
}

/* get coeffs ca of array ua, and check for falloff */
int has_expfalloff_array(tNode *node, tArray *ua, double alp[3], double s[3])
{
  DECL_STACK_ARRAY(ca, ua->n);
  basis_array_analysis3(node, ua, ca);
  return has_expfalloff_coeff_array(ca, alp, s);
}

/* check var for exp falloff */
int has_expfalloff_var(tNode *node, int ui, double alp[3], double s[3])
{
  tArray *ua;
  tDat *dat = node->dat;
  if(dat)
  {
    ua = dat->v[ui];
    return has_expfalloff_array(node, ua, alp, s);
  }
  return 1; /* return 1 if we do not have var ui on this proc */
}



/***********************************************************************/
/* functions for linear fit */
/***********************************************************************/

/*
1D least squares fit

data: (x_i, y_i),  i=1, ... , n
x_i is indep. and y_i dependent variable

fit parameters \beta = \beta_j,  j=1, ... , m

residual:  r_i = y_i - f(x_i, \beta),

least squares minimizes: S = r_i r_i
\partial S / \partial\beta_j
 = -2 r_i \partial f(x_i, \beta) / \partial\beta_j = 0

Assume f(x_i, \beta) = \beta_j \phi_j(x_i) = \beta_j X_{ij}
where X_{ij} := \phi_j(x_i)

(usually m=2: \phi_1(x_i)=1, \phi_2(x_i)=x_i)

min. cond. gives
[y_i - f(x_i, \beta)] X_{ij} = 0
y_i X_{ij} - \beta_j X_{ij} X_{ij} = 0

X^T \vec{y} - X^T X \beta = 0

==> \beta = (X^T X)^{-1} X^T \vec{y}
---------------------------------------------

here X is the n*m matrix:

     ( X_{11} X_{12} ... )
 X = ( X_{21} X_{22} ... )
     ( ...               )
*/

/*
3D least squares fit

data: (x_i, y_i, z_i, v_i),  i=1, ... , n
\vec{x}_i = (x_i, y_i, z_i) are indep. and v_i dependent variables

fit parameters: \beta_0, \beta_1, \beta_2, \beta_3,   j=0,1,2,3

residual:  r_i = v_i - f(\vec{x}_i, \beta}),

least squares minimizes: S = r_i r_i

\partial S / \partial\beta_j
 = -2 r_i \partial f(\vec{x}_i, \beta) / \partial\beta_j = 0

Assume f(\vec{x}_i, \beta) = \beta_0 x_i + \beta_1 y_i + \beta_2 z_i + \beta_3
                           = \beta_j X_{ij}
here X is the n*4 matrix:

     ( x_1  y_1  z_1  1 )
 X = ( x_2  y_2  z_2  1 )
     ( ...              )

min. cond. gives
[v_i - f(\vec{x}_i, \beta)] X_{ij} = 0
v_i X_{ij} - \beta_j X_{ij} X_{ij} = 0

X^T \vec{v} - X^T X \vec{\beta} = 0

==> \beta = (X^T X)^{-1} X^T \vec{v}

*/

/* 3d linear fit to array c_ijk: */
void linear_fit_array(tArray *c, double beta[4])
{
  tArray *X   = alloc_array2d(c->N, 4);
  tArray *XTX = alloc_array2d(4, 4);
  tArray *XTc = alloc_array2d(4, 1);
  tArray *b = alloc_empty_array2d(4, 1);
  int n[3] = { c->n[0], c->n[1], c->n[2] }; /* save c->n */


  /* make c into 1d array, i.e. a column vector */
  redim_array(c, c->N,1,1);

  /* set N*4 matrix X */
  //... FIXME!!!

  /* set 4*4 matrix XTX = X^T X */
  mm_array0_norestrict(X, X, XTX);

  /* get inverse of XTX */
  invert4x4x1symm_array(XTX); /* XTX now contains inverse of X^T X */

  /* get X^T \vec{c} */
  mm_array0(X, c, XTc);

  /* b = (X^T X)^{-1} XTc = ((X^T X)^T)^{-1} XTc = ((X^T X)^{-1})^T XTc */
  point_array_a_to_data(b, beta, 1);
  mm_array0(XTX, XTc, b);

  /* restore dims of c */
  redim_array(c, n[0],n[1],n[2]);

  /* free all arrays */
  free_array(b);
  free_array(XTX);
  free_array(X);
}
