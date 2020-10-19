/* setup_Boxes.c */
/* Wolfgang Tichy, April 2019 */
/* create various arrangements of Cartesian boxes */

#include "nmesh.h"
#include "coordinates.h"



/* add 1 pat as a box centered at xc[i], returns the index of the pat */
int add_1box_pat(tMesh *mesh, double xc[3], double dout[3])
{
  int amr_n0 = Geti(Par("amr_n0"));
  int amr_n1 = Geti(Par("amr_n1"));
  int amr_n2 = Geti(Par("amr_n2"));
  int n[] = { amr_n0, amr_n1, amr_n2 };
  double bbox[6];
  tPat *pat;
  int d;

  /* set min/max in each direction */
  for(d=0; d<3; d++)
  {
    bbox[2*d]   = xc[d] - dout[d];
    bbox[2*d+1] = xc[d] + dout[d];
  }

  /* make new patch */
  pat = add_patch(mesh, bbox, n, 0);

  return pat->p; /* return pat index */
}

/* add 1 pat as a cube centered at xc[i], returns the index of the pat */
int add_1cube_pat(tMesh *mesh, double *xc, double rout)
{
  double dout[] = { rout, rout, rout };
  return add_1box_pat(mesh, xc, dout);
}

/* add N box patches in a line in dir */
int add_Nbox_pats_indir(tMesh *mesh, double xc[3], double dout[3],
                        int N, int dir)
{
  double x[] = { xc[0], xc[1], xc[2] };
  double mid = xc[dir];
  double L = 2.*(dout[dir]);
  double c = 0.5 * (!(N%2));
  int s = N/2;
  int i, ret=-1;

  for(i=0; i<N; i++)
  {
    x[dir] = mid + L*(i - s + c);
    ret = add_1box_pat(mesh, x, dout);
  }
  return ret;
}

/* arrange N[0]*N[1]*N[2] box patches into a bigger box, centered on xc
   with side lengths dout[0], dout[1], dout[2] */
int arrange_box_pats_inBox(tMesh *mesh, double xc[3], double dout[3], int N[3])
{
  double x[] = { xc[0], xc[1], xc[2] };
  double d[] = { dout[0], dout[1], dout[2] };
  //{ dout[0], dout[1]/N[1], dout[2]/N[2] };
  int j, k, ret=-1;

  for(k = 0; k < N[2]; k++)
  {
    double L2 = 2.*(dout[2]);
    double c2 = 0.5 * (!(N[2]%2));
    int s2 = N[2]/2;

    x[2] = xc[2] + L2*(k - s2 + c2);

    for(j = 0; j < N[1]; j++)
    {
      double L1 = 2.*(dout[1]);
      double c1 = 0.5 * (!(N[1]%2));
      int s1 = N[1]/2;

      x[1] = xc[1] + L1*(j - s1 + c1);

      ret = add_Nbox_pats_indir(mesh, x, d, N[0], 0);
    }
  }
  return ret;
}
