/* fv.c */
/* Wolfgang Tichy, Dec 2020 */


#include "nmesh.h"
#include "dg.h"

/* get glabal amr vars */
extern tAMR amr[1];

/* use DGglobals */
extern tDGglobals DGglobals[1];

/***********************************************************************/
/* funcs needed for finite volume method in nmesh */
/***********************************************************************/

/* Function that reconstructs q-vars at midpoint im,
   directly from arrays qc that contain q-vars at gridpoints.
   Note: fv->qc[l][i], qm_p[l], qm_m[l] are allocated in fv_divf.
   In: nq, qc, npts, im, q_scale, rec1d_p,rec1d_m. Out: qm_p, qm_m */
void fv_rec1d_q_midpt(tFVinfo *fv)
{
  int nq      = fv->nq;
  double **qc = fv->qc;   // qc[0..nvars-1][0..npts-1]
  int npts    = fv->npts;
  int im      = fv->im;   // im = 0..npts-2, im is midpt to right of grdpt im
  int q_scale = fv->q_scale;
  int l;

  /* interpolate fields qc towards the midpoint at im */
  for(l=0; l<nq; l++)
  {
    /* reconstruct from both sides of midpoint at im */
    fv->qm_p[l] = fv->rec1d_p(npts, qc[l], im, q_scale);
    fv->qm_m[l] = fv->rec1d_m(npts, qc[l], im, q_scale);
  }
}

/* compute d_i f^i with finite vol. methods on one node.
   In:
     vlq are vars we reconstruct (we can use cons. vars u here)
     vlu are cons. vars u
     vls are extra source vars we may need (can be NULL)
     rec1d_u_f_lam_midpt = rec. cons u, fluxes f & eiegnvals at a midpoint
     numflux             = numerical flux we want
   Out:
     vldivf = div(f(u)) on all inner gridpoints and a piece of div(f(u)) on
              face points */
void fv_divf(tNode *node, tVarList *vldivf, tVarList *vlq,
             tVarList *vlu, tVarList *vls,
             void (*rec1d_u_f_lam_midpt)(tFVinfo *f, tDGinfo *d),
             void (*u_f_lam)(tDGinfo *d),
             void (*numflux)(tDGinfo *d))
{
  tMesh *mesh = vlq->mesh;
  int norms_and_sqrtgdiag_on_midpoints = 0;
  int nqvars = vlq->n;
  int nfvars = vldivf->n;
  int sqrtgdiagx = Ind("sqrtgdiagx");
  int iXm_sqrtgdiagx, iYm_sqrtgdiagx, iZm_sqrtgdiagx;
  /* func ptrs for reconstruction */
  double (*rec1d_p)(int n, const double *u, int im, double u_scale);
  double (*rec1d_m)(int n, const double *u, int im, double u_scale);
  double q_scale = 1.; /* typical order of magnitude of fields */
  int nghosts;         /* number of ghost points on each end */
  int add_surface_fluxes; /* whether we set all of divf on faces */

  if(norms_and_sqrtgdiag_on_midpoints)
  {
    iXm_sqrtgdiagx = Ind("Xm_sqrtgdiagx");
    iYm_sqrtgdiagx = Ind("Ym_sqrtgdiagx");
    iZm_sqrtgdiagx = Ind("Zm_sqrtgdiagx");
  }
  else
  {
    iXm_sqrtgdiagx = iYm_sqrtgdiagx = iZm_sqrtgdiagx = sqrtgdiagx;
  }

  /* set func ptrs for rec. */
  switch(DGglobals->fv_rec_mode)
  {
  case FV_REC_1:
    /* reconstruct from both sides of midpoint at i0m */
    rec1d_p = rec1d_p_1;
    rec1d_m = rec1d_m_1;
    nghosts = 0;
    break;
  /* use WENO3 from both sides of midpoint at i0m */
  case FV_REC_WENOm3_2:
    rec1d_p = rec1d_p_WENOm3_2;
    rec1d_m = rec1d_m_WENOm3_2;
    nghosts = 0;
    break;
  case FV_REC_WENOm5_2:
    rec1d_p = rec1d_p_WENOm5_2;
    rec1d_m = rec1d_m_WENOm5_2;
    nghosts = 0;
    break;
  case FV_REC_WENOmZ_2:
    rec1d_p = rec1d_p_WENOmZ_2;
    rec1d_m = rec1d_m_WENOmZ_2;
    nghosts = 0;
    break;
  /* WENO3 experiments that didn't show good convergence: */
  case FV_REC_WENO3if2away_1:
    rec1d_p = rec1d_p_WENO3_if2away;
    rec1d_m = rec1d_m_WENO3_if2away;
    nghosts = 0;
    break;
  case FV_REC_WENO3if1away_1:
    rec1d_p = rec1d_p_WENO3_if1away;
    rec1d_m = rec1d_m_WENO3_if1away;
    nghosts = 0;
    break;
  case FV_REC_WENO3_2:
    rec1d_p = rec1d_p_WENO3_2;
    rec1d_m = rec1d_m_WENO3_2;
    nghosts = 0;
    break;
  case FV_REC_WENO3_2g:
    rec1d_p = rec1d_p_WENO3_2g;
    rec1d_m = rec1d_m_WENO3_2g;
    nghosts = 1;
    break;
  default:
    errorexit("unknown DGglobals->fv_rec_mode");
  }

  /* do we add in the surface fluxes here already? */
  add_surface_fluxes = DGglobals->fv_divf_adds_surface_fluxes;

  /* set var list for div of fluxes to zero */
  vlsetconstant_node(node, vldivf, 0.);

  /* RHS */
  {
    tDGinfo *d = alloc_DGinfo(vlu, vls);
    double *m_sqrtgdiag[3][3] =
      { { Vard(node, iXm_sqrtgdiagx), Vard(node, iXm_sqrtgdiagx+1),
                                           Vard(node, iXm_sqrtgdiagx+2) },
        { Vard(node, iYm_sqrtgdiagx), Vard(node, iYm_sqrtgdiagx+1),
                                           Vard(node, iYm_sqrtgdiagx+2) },
        { Vard(node, iZm_sqrtgdiagx), Vard(node, iZm_sqrtgdiagx+1),
                                           Vard(node, iZm_sqrtgdiagx+2) } };
    int *n = node->n;
    int maxn = max3(n[0],n[1],n[2]);
    double *Xbm = dmalloc(maxn);
    double *dXb = dmalloc(maxn);
    double *qc[nqvars];          //pointers to data of the q-fields
    double *fnumR[nfvars];       //pointers to data of the fluxes
    int npg = maxn + 2*nghosts;  //number of points in qcg[l]
    int npe = maxn + 2;          //number of points in fnumRe[l]
    double (*qcg)[npg] = dtensor(nqvars*npg);     //array for the q-fields
    double (*fnumRe)[npe] = calloc(nfvars, sizeof *fnumRe); //for fluxes
    double (*di0fi0)[maxn] = dtensor(nfvars*maxn); //array for d_i flux^i
    double *qm_p = dmalloc(nqvars); // array for rec u at one point
    double *qm_m = dmalloc(nqvars);
    int l; /* field index */
    int dir;

    /* set qc to part of qcg without ghosts */
    for(l=0; l<nqvars; l++) qc[l] = &(qcg[l][nghosts]);
    /* NOTE: now qc[l][-1] = qcg[l][0] i.e. ghost on left */
    for(l=0; l<nfvars; l++) fnumR[l] = &(fnumRe[l][1]);

    /* write node into d because numflux needs this */
    d->node = node;
    if(norms_and_sqrtgdiag_on_midpoints)
      d->info = DGINFO_MIDPT;
    else
      d->info = DGINFO_NULL;

    /* get nbsurf and ajsurf already */
    if(nghosts || add_surface_fluxes) get_all_surfaces(node);

    /* add fluxes in each direction to RHS */
    for(dir=0; dir<3; dir++)
    {
      int dir_active = Getb(amr->dir_active[dir]);
      double *sqrtgdiag = m_sqrtgdiag[dir][dir];
      int i,j,k;

      /* get midpoints */
      set_nm_nodemidpoints_Xb_dir(node, n[dir]-1,0, dir, Xbm);
      set_nm_nodemidpoint_distsXb_dir(node, dir, Xbm, dXb);

      /* loop over plane */
      forplaneN(dir, i,j,k, n, 0)
      {
        int i1 = i1_norm(i,j,k, dir); /* 1st and 2nd index in plane */
        int i2 = i2_norm(i,j,k, dir);
        int i0;                       /* index orthogonal to plane */
        int ic,jc,kc, ccc;

        /* fill field arrays qc, i0 runs orth. to plane */
        for(i0=0; i0<n[dir]; i0++)
        {
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          forvl(vlq, l)
          {
            double *q = Vard(node, Vind(vlq, l));
            qc[l][i0] = q[ccc];
          }
        }
        /* fill in ghost points in qc with value from adjacent node */
        if(nghosts)
        {
          int JK = Ind_n_norm(i,j,k, n, dir);
          forvl(vlq, l)
          {
            int vi = Vind(vlq, l);
            double *qaj;
            /* put adj. val on left face in left ghost: qc[l][-1] */
            qaj = Varaj(node, vi, dir*2);
            if(qaj) qc[l][-1] = qaj[JK];
            else    qc[l][-1] = 1e30; //large value that WENO should ignore
            /* get adj. val on right face in right ghost: qc[l][n] */
            qaj = Varaj(node, vi, dir*2+1);
            if(qaj) qc[l][n[dir]] = qaj[JK];
            else    qc[l][n[dir]] = 1e30;
          }
        }

        /* loop over points in dir */
        for(i0=0; i0<n[dir]; i0++)
        {
          int im0, im0m1, im,jm,km, cccm1;
          double i0g0, i0lN, gd_ow_m, gd_ow_m1, wm;

          /* if we have only 1 point do nothing, as there are no midpoints */
//          if(n[dir]<=1) break;

          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          /* set 1d index of right and left midpoint and some flags if we
             are at endpoints */
          if(i0>0) { i0g0=1; im0m1 = i0-1; }
          else     { i0g0=0; im0m1 = i0; /* safe value */ }
          if(i0<n[dir]-1) { i0lN=1; im0 = i0; }
          else            { i0lN=0; im0 = i0-1; /* safe value */ }

          /* set left midpoint index */
          ijk_inplaneN(dir, im,jm,km, i1,i2,im0m1);
          cccm1 = Ind_n(im,jm,km, n);

          /* set index and face of the right midpoint */
          d->i = ic;
          d->j = jc;
          d->k = kc;
          d->face = dir*2 + 1;

          /* if i0 has a midpoint to its right */
          if(i0<n[dir]-1)
          {
            tFVinfo fv[1];

            /* set fv  */
            fv->nq = nqvars;
            fv->qc = qc;
            fv->npts = n[dir];
            fv->im = im0;
            fv->q_scale = q_scale;
            fv->rec1d_p = rec1d_p;
            fv->rec1d_m = rec1d_m;
            fv->qm_p = qm_p;
            fv->qm_m = qm_m;

            /* reconstruct q,u and then set fluxes and eigenvalues in d */
            rec1d_u_f_lam_midpt(fv, d);

            /* compute numerical flux directly after rec1d_u_f_lam_midpt,
               if not set already in rec1d_u_f_lam_midpt */
            if(numflux) numflux(d);

            /* save d->fnum in fnumR for each field and point */
            forvl(vldivf, l)
              fnumR[l][im0] = d->fnum[l];
            /* we could now also compute fnumL with normL=-normR, but I think
               this results in fnumL = -fnumR */
            //printDGinfo(d);
          }

          /* factors in flux terms on RHS at right and left midpoint */
          wm  = dXb[i0];
          if(norms_and_sqrtgdiag_on_midpoints)
          {
            gd_ow_m  = sqrtgdiag[ccc]/wm;
            gd_ow_m1 = sqrtgdiag[cccm1]/wm;
          }
          else /* get sqrtgdiag on grid points */
          {
            gd_ow_m  = sqrtgdiag[ccc]/wm;
            gd_ow_m1 = sqrtgdiag[ccc]/wm;
          }

          /* include flux terms on facepoints */
          if(add_surface_fluxes && dir_active)
          {
            int d_face_sav = d->face; /* save parts of d we may alter */
            int d_info_sav = d->info;

            if(i0 == 0)
            {
              /* set non-zero weight on left boundary */
              gd_ow_m1 = sqrtgdiag[ccc]/wm;

              /* compute numerical fluxes on the left side of node */
              d->face = dir*2;        /* normal points to the left */
              d->info = DGINFO_NULL;  /* facepoint is grid point */
              u_f_lam(d);
              numflux(d);

              /* save fluxes of left face in fnumR, use fnumR = -fnumL */
              forvl(vldivf, l) fnumR[l][-1] = -( d->fnum[l] );
            }
            if(i0 == n[dir]-1)
            {
              /* set non-zero weight on right boundary */
              gd_ow_m = sqrtgdiag[ccc]/wm;

              /* compute numerical fluxes on the right side of node */
              d->face = dir*2 + 1;   /* normal points to the right */
              d->info = DGINFO_NULL; /* facepoint is grid point */
              u_f_lam(d);
              numflux(d);

              /* save fluxes of right face in fnumR */
              forvl(vldivf, l) fnumR[l][i0] = d->fnum[l];
            }

            /* restore altered parts of d */
            d->face = d_face_sav;
            d->info = d_info_sav;
          }
          else
          {
            /* Set factors in flux on faces to zero.
               Note: i0g0=0 on left face, i0lN=0 on right face */
            gd_ow_m  = i0lN * gd_ow_m;
            gd_ow_m1 = i0g0 * gd_ow_m1;
          }

          //printf("i0=%d im0=%d im0m1=%d: wm=%g gd_ow_m=%g gd_ow_m1=%g\n",
          //i0, im0, im0m1, wm, gd_ow_m, gd_ow_m1);

          /* get piece of div(flux) in direction dir with FV method */
          forvl(vldivf, l)
          {
            double *fnum = fnumR[l];
            double *df = di0fi0[l];
            df[i0] = fnum[i0]*gd_ow_m - fnum[i0-1]*gd_ow_m1;
          }
        } /* end i0 loop */

        /* final loop over points in dir */
        for(i0=0; i0<n[dir]; i0++)
        {
          /* set points and their index */
          ijk_inplaneN(dir, ic,jc,kc, i1,i2, i0);
          ccc = Ind_n(ic,jc,kc, n);

          /* add to divf */
          forvl(vldivf, l)
          {
            int idivf = Vind(vldivf,l);
            double *divf = Vard_(node, idivf);
            double *df = di0fi0[l];

            // interpol. df at ends
            // ...
            divf[ccc] += df[i0];
          }
        }
      } /* end plane loop */
    } /* end dir-loop*/

    /* release mem */
    free(qm_m);
    free(qm_p);
    free(di0fi0);
    free(fnumRe);
    free(qcg);
    free(dXb);
    free(Xbm);
    free_DGinfo(d);
  }
}
