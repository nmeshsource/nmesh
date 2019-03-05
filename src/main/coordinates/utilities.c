/* utilities.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "coordinates.h"



/* find normal vector (n[0],n[1],n[2]) of patch face f at X,Y,Z */
double patch_normal_at_XYZ(tPat *pat, int f, const double X[3], double n[3])
{
  int dir = f/2;         /* get direction */
  int sig = 2*(f%2) - 1; /* get sign for outward direction */
  double smag;
  int j;

  if(f>5 || f<0) errorexit("f must be 0,1,2,3,4,5");

  if(pat->dXYZ_dxyz) /* not Cartesian */
  {
    double x[3], dXYZdxyz[3][3];

    /* set dXYZdxyz */
    pat->dXYZ_dxyz(pat, NULL,-1, X, x, dXYZdxyz);
    for(j=0; j<3; j++) n[j] = dXYZdxyz[dir][j];
  }
  else  /* use Cartesian normal */
  {
    n[0] = n[1] = n[2] = 0.0;
    n[dir] = 1.0;
  }

  /* normalize and set sign from sig */
  smag = sig * sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if(smag == 0.0) smag = sig;
  n[0] /= smag;
  n[1] /= smag;
  n[2] /= smag;
  return smag;
}

/* find normal vector (n[0],n[1],n[2]) of pat face f at point ijk */
double node_normal_at_ijk(tNode *node, int f, int ijk, double n[3])
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  tDat *dat = node->dat;
  int dir = f/2;     /* get direction */
  int sig = 2*(f%2) - 1; /* get sign for outward direction */
  int idXd = Ind("dXdx");
  double smag;

  if(f>5 || f<0) errorexit("f must be 0,1,2,3,4,5");

  if(!dat) return 0.;

  /* do we need to init. coords? */
  if(!(dat->coords_set)) coordinates_init_node(node);

  /* get normal from derivs */
  {
    double *pdXd[3][3]
            = { {Vard(node,idXd),   Vard(node,idXd+1), Vard(node,idXd+2)},
                {Vard(node,idXd+3), Vard(node,idXd+4), Vard(node,idXd+5)},
                {Vard(node,idXd+6), Vard(node,idXd+7), Vard(node,idXd+8)} };
    n[0] = pdXd[dir][0][ijk];
    n[1] = pdXd[dir][1][ijk];
    n[2] = pdXd[dir][2][ijk];
  }

  /* normalize and set sign from sig */
  smag = sig * sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if(smag == 0.0) smag = sig;
  n[0] /= smag;
  n[1] /= smag;
  n[2] /= smag;
  return smag;
}

/* get transpose of a 3x3 matrix. */
void transp3Dmat_from_3Dmat(const double M[3][3], double transpM[3][3])
{
  transpM[0][0] = M[0][0];
  transpM[0][1] = M[1][0];
  transpM[0][2] = M[2][0];
  transpM[1][0] = M[0][1];
  transpM[1][1] = M[1][1];
  transpM[1][2] = M[2][1];
  transpM[2][0] = M[0][2];
  transpM[2][1] = M[1][2];
  transpM[2][2] = M[2][2];
}

/* determinant of 3x3 matrix */
double det_3Dmatrix(const double M[3][3])
{
  /* M = {{m11,m12,m13},{m21,m22,m23},{m31,m32,m33}}
     Det[M] = m11*m22*m33 - m13*m22*m31 + m12*m23*m31 +
              m13*m21*m32 - m11*m23*m32 - m12*m21*m33   */
  LDOUBLE DetM =M[0][0]*M[1][1]*M[2][2] -
                M[0][2]*M[1][1]*M[2][0] +
                M[0][1]*M[1][2]*M[2][0] +
                M[0][2]*M[1][0]*M[2][1] -
                M[0][0]*M[1][2]*M[2][1] -
                M[0][1]*M[1][0]*M[2][2];
  return DetM;
}

/* 3d matrix inverse and return det. */
double inv3Dmat_from_3Dmat(const double M[3][3], double invM[3][3])
{
  /* M = {{m11,m12,m13},{m21,m22,m23},{m31,m32,m33}}
     Det[M] = m11*m22*m33 - m13*m22*m31 + m12*m23*m31 +
              m13*m21*m32 - m11*m23*m32 - m12*m21*m33
     Inverse[M]*Det[M] =
       {{-(m23*m32) + m22*m33, m13*m32 - m12*m33, -(m13*m22) + m12*m23},
        {m23*m31 - m21*m33, -(m13*m31) + m11*m33, m13*m21 - m11*m23},
        {-(m22*m31) + m21*m32, m12*m31 - m11*m32, -(m12*m21) + m11*m22}}  */
  LDOUBLE DetM = det_3Dmatrix(M);
  LDOUBLE sum;
  if(DetM==0.0) DetM=dequaleps*dequaleps*dequaleps;
  sum        = (-(M[1][2]*M[2][1]) + M[1][1]*M[2][2])/DetM;
  invM[0][0] = sum;
  sum        = (  M[0][2]*M[2][1]  - M[0][1]*M[2][2])/DetM;
  invM[0][1] = sum;
  sum        = (-(M[0][2]*M[1][1]) + M[0][1]*M[1][2])/DetM;
  invM[0][2] = sum;
  sum        = (  M[1][2]*M[2][0]  - M[1][0]*M[2][2])/DetM;
  invM[1][0] = sum;
  sum        = (-(M[0][2]*M[2][0]) + M[0][0]*M[2][2])/DetM;
  invM[1][1] = sum;
  sum        = (  M[0][2]*M[1][0]  - M[0][0]*M[1][2])/DetM;
  invM[1][2] = sum;
  sum        = (-(M[1][1]*M[2][0]) + M[1][0]*M[2][1])/DetM;
  invM[2][0] = sum;
  sum        = (  M[0][1]*M[2][0]  - M[0][0]*M[2][1])/DetM;
  invM[2][1] = sum;
  sum        = (-(M[0][1]*M[1][0]) + M[0][0]*M[1][1])/DetM;
  invM[2][2] = sum;

  return DetM;
}



/* find pat size L of pat */
double find_pat_size(tPat *pat)
{
  int i,j,n;
  double dx,dy,dz;
  double sum, L;
  double *bb = pat->bbox;
  double X[8][3], x[8][3];

  /* max and min X */
  X[0][0] = bb[0];  X[0][1] = bb[2];  X[0][2] = bb[4];
  X[1][0] = bb[1];  X[1][1] = bb[2];  X[1][2] = bb[4];
  X[2][0] = bb[0];  X[2][1] = bb[3];  X[2][2] = bb[4];
  X[3][0] = bb[1];  X[3][1] = bb[3];  X[3][2] = bb[4];
  X[4][0] = bb[0];  X[4][1] = bb[2];  X[4][2] = bb[5];
  X[5][0] = bb[1];  X[5][1] = bb[2];  X[5][2] = bb[5];
  X[6][0] = bb[0];  X[6][1] = bb[3];  X[6][2] = bb[5];
  X[7][0] = bb[1];  X[7][1] = bb[3];  X[7][2] = bb[5];

  /* x at max and min X */
  for(i=0; i<8; i++)
    pat->xyz_of_XYZ(pat, NULL,-1, X[i], x[i]);

  sum=0.0;
  n=0;
  for(i=1; i<8; i++)
  for(j=0; j<i; j++)
  {
    dx = x[i][0] - x[j][0];
    dy = x[i][1] - x[j][1];
    dz = x[i][2] - x[j][2];

    sum += dx*dx + dy*dy + dz*dz;
    n++;
  }

  L = sqrt(sum/n);
  return L;
}

/* find pat size L of smallest pat */
double smallest_pat_size(tMesh *mesh)
{
  int i;
  double L=1e300;
  forpatches(mesh, i)
  {
    tPat *pat = mesh->pat[i];
    double Li = find_pat_size(pat);
    if(Li<L) L=Li;
  }
  return L;
}

/* determine which patch faces touch or overlap other patches,
   or touch the outside of our pat, returns extface */
/* extface[f]=1 means face f is external, i.e. needs BC */
/* Faces with e.g. periodic coordinates are marked as extface[f]=0. */
/* If inclOuterBound=1 we mark faces that are not in contact with any other
   pat as external (because they probably need an outer BC) */
void find_external_faces_of_pat(tPat *pat, int *extface, int inclOuterBound)
{
  tMesh *mesh = pat->mesh;
  int n[] = { 16,16,16 };     /* we use 16 points */
  double X0[3], LX[3], dX[3]; /* grid of points */
  int nd, dir, dd;
  int i,j,k, pl, f;
  intList *opl; /* list that contains other patches */
  int op;
  double oX[3];
  double L;
  double EPS = 1e-5;

  /* mark faces in periodic dirs with not external, i.e. extface[f]=0  */
  for(f=0; f<6; f++)
  {
    if(pat->periodic[f/2]) extface[f] = 0;
    else                   extface[f] = 1;
  }

  if(pat->dXYZ_dxyz==NULL) /* Cartesian pat */
  {
    if(inclOuterBound) return; /* do nothing else for a Cartesian pat */
  }

  /* find pat size L of smallest pat */
  L = smallest_pat_size(mesh);

  opl = alloc_intList(); /* list that will contain other patches */

  /* make opl that contains all pats except this one */
  forpatches(mesh, i) if(i!=pat->p) unionpush_intList(opl, i);

  /* make a grid of points, that excludes endpoints */
  for(dir=0; dir<3; dir++)
  {
    X0[dir] = pat->bbox[2*dir];
    LX[dir] = pat->bbox[2*dir+1] - X0[dir];
    dX[dir] = LX[dir]/(n[dir]);
    X0[dir] += dX[dir] * 0.5;
  }

  /* go over directions */
  for(dir=0; dir<3; dir++)
  {
    nd = n[dir];

    /* go over faces, but not edges */
    for(pl=0; pl<nd; pl+=nd-1)
    {
      /* pick face index */
      if(pl==0) f=2*dir;
      else      f=2*dir+1;

      /* do nothing if not an external face */
      if(!extface[f]) continue;

      /* look for points in other pates just outside this pat */
      forinnerplaneN(dir, i,j,k, n, pl)
      {
        double X[3], x[3];
        double Ndir[3];
        double ox[3], dx[3], d0, ret;

        /* point grid */
        X[0] = X0[0] + dX[0] * i;
        X[1] = X0[1] + dX[1] * j;
        X[2] = X0[2] + dX[2] * k;

        /* pick one of X,Y,Z on boundary */
        X[dir] = pat->bbox[f];

        /* get x from X */
        if(pat->xyz_of_XYZ)
          pat->xyz_of_XYZ(pat, NULL,-1, X, x);
        else
          { x[0]=X[0]; x[1]=X[1]; x[2]=X[2]; }

        /* use normal vector to find point ox slightly outside pat */
        patch_normal_at_XYZ(pat, f, X, Ndir);
        for(d0=0., dd=0; dd<3; dd++) d0 += x[dd]*x[dd];
        d0 = sqrt(d0);
        for(dd=0; dd<3; dd++)
        {
          dx[dd] = Ndir[dd]*(L+d0)*EPS;
          ox[dd] = x[dd] + dx[dd];
        }

        /* find point in this pat */
        ret = p_XYZ_of_xyz(pat, oX, ox);
        if(ret>=0)
          for(dd=0; dd<3; dd++)
          {
            if(!(pat->periodic[dd]))
              if(oX[dd] < pat->bbox[2*dd] || pat->bbox[2*dd+1] < oX[dd])
                ret=-1;
          }
        if(ret>=0 && d0<1e60) /* point is in this pat and x,y,z is not inf  */
        {
          extface[f]=0; /* mark face as not external */
          //if(extface[f]==0) errorexit("1");
          goto endplaneloop; /* break; does not work for nested loop */
        }

        /* find point in other pates */
        op = p_XYZ_of_xyz_inpatlist(mesh, opl, oX, x);

        /* if we find one point in another pat this face is external */
        if(op>=0) break; /* leave plane loop if face is external*/
      } endplaneloop:

      /* check about including outer boundaries */
      if(inclOuterBound)
      {
        /* if op=p the other pat is the pat itself,
           so it's not an external face */
        if(op==pat->p) extface[f]=0;
        /* NOTE currently opl does not contain b, so ob=b cannot happen!!! */
        //if(extface[f]==0) errorexit("2");
      }
      else
      {
        /* if op<0, we found no other pat face and f is not external */
        if(op<0) extface[f]=0;
        //if(extface[f]==0) errorexit("3");
      }
    }
  } /* end loop over directions */

  free_intList(opl);
}
