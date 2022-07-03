/* tests.c */
/* Wolfgang Tichy, 7/2022 */

#include "nmesh.h"
#include "coordinates.h"


/* global vars */
extern tcoordinates coordinates[1];




/* Compute \partial_{\bar{i}} (J \sqrt{g^{\bar{i}\bar{i}}} n^{\bar{i}}_i)
   for each i. Save in var divb_J_sqrtgdiag_nx.
   This should converge to zero, if sqrtgdiag is obtained by transforming
   the flat metric. */
int coordinates_set_divb_J_sqrtgdiag_n(tNode *node)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  //tDat *dat = node->dat;

  //int midpoint_data = coordinates->midpoint_data;
  //int sqrtdet2g_o_det3gamma = coordinates->sqrtdet2g_o_det3gamma;
  //int idXdx = coordinates->idXdx;
  int itmp1 = coordinates->itmp1;
  int itmp2 = itmp1 + 1;
  int idivb_J_sqrtgdiag_nx = Ind("divb_J_sqrtgdiag_nx");
  int idet_dXbdx = coordinates->idet_dXbdx;
  //int isqrtdet2g_o_det3gamma0 = coordinates->isqrtdet2g_o_det3gamma0;
  int isqrtgdiagx = coordinates->isqrtgdiagx;
  //int iXm_det_dXbdx = coordinates->iXm_det_dXbdx;
  //int iYm_det_dXbdx = coordinates->iYm_det_dXbdx;
  //int iZm_det_dXbdx = coordinates->iZm_det_dXbdx;
  //int iXm_sqrtgdiagx = coordinates->iXm_sqrtgdiagx;
  //int iYm_sqrtgdiagy = coordinates->iYm_sqrtgdiagy;
  //int iZm_sqrtgdiagz = coordinates->iZm_sqrtgdiagz;

  tVarList *vltmp2 = vlalloc(mesh);
  tVarList *vldivb_J_sqrtgdiag_n = vlalloc(mesh);
  int ib, i, ijk;

  vlpush(vltmp2, itmp2);

  /* loop over i of n^{\bar{i}_i */
  for(i=0; i<3; i++)
  {
    vlpushone(vldivb_J_sqrtgdiag_n, idivb_J_sqrtgdiag_nx+i);
    vlsetconstant_node(node, vldivb_J_sqrtgdiag_n, 0.);

    /* loop over ib=\bar{i} of n^{\bar{i}_i */
    for(ib=0; ib<3; ib++)
    {
      double *tmp1 = Vard(node, itmp1);
      double *sqrtgdiag = Vard(node, isqrtgdiagx+ib);
      double *ooJ = Vard(node, idet_dXbdx);

      /* set tmp1 = \sqrt{g^{\bar{i}\bar{i}}} J * n^{\bar{i}}_i
         for each \bar{i}=ib */
      forpoints(node,ijk)
      {
        int f = 2*ib + 1; /* results in right pointing normal */
        double norm[3];
        node_normal_at_ijk(node, f, ijk, norm);
        tmp1[ijk] = (sqrtgdiag[ijk] / ooJ[ijk]) * norm[i];
      }

      /* tmp2 = \partial_{\bar{i}} tmp1 */
      basis_var_deriv1(node, ib, itmp1, itmp2, NULL);
      /* divb_J_sqrtgdiag_n += tmp2 */
      vladdto_node(node, vldivb_J_sqrtgdiag_n, 1., vltmp2);
    }
    vldropn(vldivb_J_sqrtgdiag_n, 1);
  }

  vlfree(vldivb_J_sqrtgdiag_n);
  vlfree(vltmp2);
  return 0;
}

/* Compute (1/J) D_{\bar{i}} (J \sqrt{g^{\bar{i}\bar{i}}} n^{\bar{i}}_i)
   for each i. Here Db = D_{\bar{i}} is defined as the difference between
   the values at the two midpoints aroound each grid point.
   Save in var ooJ_Db_J_sqrtgdiag_nx.
   This should converge to zero, if sqrtgdiag is obtained by transforming
   the flat metric. */
int coordinates_set_ooJ_Db_J_sqrtgdiag_n(tNode *node)
{
  tPat *pat = node->pat;
  tMesh *mesh = pat->mesh;
  //tDat *dat = node->dat;
  int *n = node->n;

  //int midpoint_data = coordinates->midpoint_data;
  //int sqrtdet2g_o_det3gamma = coordinates->sqrtdet2g_o_det3gamma;
  //int idXdx = coordinates->idXdx;
  //int itmp1 = coordinates->itmp1;
  int idet_dXbdx = coordinates->idet_dXbdx;
  //int isqrtdet2g_o_det3gamma0 = coordinates->isqrtdet2g_o_det3gamma0;
  int isqrtgdiagx = coordinates->isqrtgdiagx;
  int iXm_det_dXbdx = coordinates->iXm_det_dXbdx;
  int iYm_det_dXbdx = coordinates->iYm_det_dXbdx;
  int iZm_det_dXbdx = coordinates->iZm_det_dXbdx;
  int im_det_dXbdx[] = { iXm_det_dXbdx, iYm_det_dXbdx, iZm_det_dXbdx };
  int iXm_sqrtgdiagx = coordinates->iXm_sqrtgdiagx;
  int iYm_sqrtgdiagy = coordinates->iYm_sqrtgdiagy;
  int iZm_sqrtgdiagz = coordinates->iZm_sqrtgdiagz;
  int im_sqrtgdiag[] = { iXm_sqrtgdiagx, iYm_sqrtgdiagy, iZm_sqrtgdiagz };
  int iooJ_Db_J_sqrtgdiag_nx = Ind("ooJ_Db_J_sqrtgdiag_nx");

  tVarList *vlooJ_Db_J_sqrtgdiag_n = vlalloc(mesh);
  int maxn = max3(n[0],n[1],n[2]);
  double *Xbm = dmalloc(maxn);
  double *dXb = dmalloc(maxn);
  int ii, dir;

  /* loop over ii=i of n^{\bar{i}_i */
  for(ii=0; ii<3; ii++)
  {
    /* zero var w. index iooJ_Db_J_sqrtgdiag_nx+ii */
    vlpushone(vlooJ_Db_J_sqrtgdiag_n, iooJ_Db_J_sqrtgdiag_nx+ii);
    vlsetconstant_node(node, vlooJ_Db_J_sqrtgdiag_n, 0.);

    /* loop over dir=\bar{i} of n^{\bar{i}_i */
    for(dir=0; dir<3; dir++)
    {
      double *sqrtgdiagm = Vard(node, im_sqrtgdiag[dir]);
      double *ooJm = Vard(node, im_det_dXbdx[dir]);
      double *sqrtgdiag = Vard(node, isqrtgdiagx+dir);
      double *ooJ = Vard(node, idet_dXbdx);
      double *ooJ_Db_J_sqrtgdiag_n
               = Vard(node, Vind(vlooJ_Db_J_sqrtgdiag_n, 0));
      int i,j,k;

      /* get midpoints */
      set_nm_nodemidpoints_Xb_dir(node, n[dir]-1,0, dir, Xbm);
      set_nm_nodemidpoint_distsXb_dir(node, dir, Xbm, dXb);

      /* set ooJ_Db_J_sqrtgdiag_nx =
           (1/J) D_{\bar{i}} (J \sqrt{g^{\bar{i}\bar{i}}} n^{\bar{i}}_i)
         for each \bar{i}=dir */
      //forinnerplaneN(dir, i,j,k, n, 0)
      forplaneN(dir, i,j,k, n, 0)
      {
        int i1 = i1_norm(i,j,k, dir); /* 1st and 2nd index in plane */
        int i2 = i2_norm(i,j,k, dir);
        int i0;                       /* index orthogonal to plane */

        /* i0 runs orth. to plane */
        for(i0=0; i0<n[dir]; i0++)
        {
          int ic,jc,kc, ccc;          /* index of gridpoints */
          int im0,   cccR;            /* index of right midpoints */
          int im0m1, cccL;            /* index of left midpoints */
          int im,jm,km;
          double i0g0, i0lN, wc, tmp;
          double normR[3], normL[3];
          double norm[3];
          double Jgd_R, Jgd_L;

//if(i0<=0 || i0>=n[dir]-1) continue;

          /* set 1d index of left and right midpoint and some flags if we
             are at endpoints */
          if(i0>0) { i0g0=1; im0m1 = i0-1; }
          else     { i0g0=0; im0m1 = i0; /* safe value */ }
          if(i0<n[dir]-1) { i0lN=1; im0 = i0; }
          else            { i0lN=0; im0 = i0-1; /* safe value */ }

          /* gridpoint index and weight */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);
          wc = dXb[i0];

          if(i0g0 && i0lN) /* in middle */
          {
            /* set right midpoint index */
            ijk_inplaneN(dir, im,jm,km, i1,i2,im0);
            cccR = Ind_n(im,jm,km, n);

            /* set left midpoint index */
            ijk_inplaneN(dir, im,jm,km, i1,i2,im0m1);
            cccL = Ind_n(im,jm,km, n);

            node_normal_at_midpt_ijk(node, 2*dir+1, cccR, normR);
            node_normal_at_midpt_ijk(node, 2*dir, cccL, normL);
            //// DO I have correct cccL in norm of fv_divf???
            //This does the same:
            //node_normal_at_midpt_right_of_ijk(node, 2*dir+1, ccc, normR);
            //node_normal_at_midpt_left_of_ijk(node, 2*dir, ccc, normL);

            Jgd_R = sqrtgdiagm[cccR] / ooJm[cccR];
            Jgd_L = sqrtgdiagm[cccL] / ooJm[cccL];
          }
          else if(i0g0==0) /* left end */
          {
            /* set right midpoint index */
            ijk_inplaneN(dir, im,jm,km, i1,i2,im0);
            cccR = Ind_n(im,jm,km, n);

            /* left midpoint is i=0 gridpoint */
            cccL = ccc;

            node_normal_at_midpt_ijk(node, 2*dir+1, cccR, normR);
            node_normal_at_ijk(node, 2*dir, ccc, normL);
            //This does the same:
            //node_normal_at_midpt_right_of_ijk(node, 2*dir+1, ccc, normR);
            //node_normal_at_midpt_left_of_ijk(node, 2*dir, ccc, normL);

            Jgd_R = sqrtgdiagm[cccR] / ooJm[cccR];
            Jgd_L = sqrtgdiag[ccc] / ooJ[ccc];
          }
          else /* right end */
          {
            /* right midpoint is i=n[dir]-1 gridpoint */
            cccR = ccc;

            /* set left midpoint index */
            ijk_inplaneN(dir, im,jm,km, i1,i2,im0m1);
            cccL = Ind_n(im,jm,km, n);

            node_normal_at_ijk(node, 2*dir+1, cccR, normR);
            node_normal_at_midpt_ijk(node, 2*dir, cccL, normL);
            //This does the same:
            //node_normal_at_midpt_right_of_ijk(node, 2*dir+1, ccc, normR);
            //node_normal_at_midpt_left_of_ijk(node, 2*dir, ccc, normL);

            Jgd_R = sqrtgdiag[ccc] / ooJ[ccc];
            Jgd_L = sqrtgdiagm[cccL] / ooJm[cccL];
          }

          wc = dXb[i0];
          tmp =  Jgd_R * normR[ii] + Jgd_L * normL[ii];
          tmp *= ooJ[ccc]/wc;

//JUNK:
          node_normal_at_ijk(node, 2*dir+1, cccR, norm);

          tmp = ooJm[cccR] - ooJ[cccR];
          tmp = sqrtgdiagm[cccR] - sqrtgdiag[cccR];

          tmp = normR[ii] -  norm[ii];

          tmp =  (sqrtgdiagm[cccR] / ooJm[cccR]) * normR[ii]
                -(sqrtgdiag[cccR] / ooJ[cccR]) * norm[ii];

          node_normal_at_ijk(node, 2*dir+1, cccL, norm);
          tmp =  (sqrtgdiagm[cccL] / ooJm[cccL]) * normL[ii]
                +(sqrtgdiag[cccL] / ooJ[cccL]) * norm[ii];

//          node_normal_at_ijk(node, 2*dir+1, ccc, norm);
//          tmp =  (sqrtgdiagm[cccL] / ooJm[cccL]) * normL[ii]
//                +(sqrtgdiag[ccc] / ooJ[ccc]) * norm[ii];

          tmp =  (sqrtgdiagm[cccR] / ooJm[cccR]) * normR[ii]
                +(sqrtgdiagm[cccL] / ooJm[cccL]) * normL[ii];

          tmp =  (sqrtgdiagm[cccR] / ooJm[cccR]) * normR[ii]
                +(sqrtgdiagm[cccL] / ooJm[cccL]) * normL[ii];
          tmp *= ooJ[ccc]/wc;
//End JUNK



          tmp =  Jgd_R * normR[ii] + Jgd_L * normL[ii];
          /* ^--this term should be extrapolated to the boundary if
             one of the summands was not constructed on a real midpoint */
          tmp *= ooJ[ccc]/wc;

          ooJ_Db_J_sqrtgdiag_n[ccc] += tmp;
        }
      } /* end forplaneN */
    }
    vldropn(vlooJ_Db_J_sqrtgdiag_n, 1);
  }

  free(dXb);
  free(Xbm);
  vlfree(vlooJ_Db_J_sqrtgdiag_n);
  return 0;
}


/*************************************************************************/
/* run all the tests */
/*************************************************************************/

/* run the tests */
int coordinates_tests(tMesh *mesh)
{
  enablevar(mesh, Ind("divb_J_sqrtgdiag_nx"));
  enablevar(mesh, Ind("ooJ_Db_J_sqrtgdiag_nx"));

  formylnodes(mesh)
  {
    tNode *node = MyLnode;

    coordinates_set_divb_J_sqrtgdiag_n(node);
    coordinates_set_ooJ_Db_J_sqrtgdiag_n(node);
  }
  return 0;
}
