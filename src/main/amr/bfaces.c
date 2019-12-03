/* bfaces.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "amr.h"

#define EPS1 1e-5
#define EPS2 1e-11
#define NPOINTS 16



/*************************************************************************/
/* funcs to add and remove bfaces */
/*************************************************************************/

/* get mem for 1 bface */
tBface *alloc_bface(tPat *pat, int f)
{
  tBface *bface = calloc(1, sizeof(bface[0]));
  if(!bface) errorexit("not enough memory");
  bface->pat = pat;
  bface->f   = f;
  return bface;
}

/* add a bface to a pat, f denotes the pat face (0 to 5),
   return the new bface */
tBface *add_empty_bface(tPat *pat, int f)
{
  //tMesh *mesh = pat->mesh;
  tBface *bface = alloc_bface(pat, f);
  tBface *bface0 = pat->bfaces[f];
  tBface *bf;

  /* add bface to end of list in pat */
  if(bface0)
  {
    for(bf=bface0; bf->next; bf=bf->next) ;
    bf->next = bface;
    bface->prev = bf;
  }
  else
    pat->bfaces[f] = bface;

  /* set some bface info */
  bface->op = -1; /* other patch not known yet */
  //bface->ioC0_0 = Ind("oC0_0"); /* var indices of other coords */
  bface->ioC0_0 = -1; /* do not use oC0_0 !!! */

  return bface;
}

/* remove a bface from a pat, return number of bfaces removed */
int remove_bface(tBface *bface)
{
  tPat *pat;
  tBface *obface, *bn, *bp;
  int f;

  if(!bface) return 0;

  pat = bface->pat;
  obface = bface->obface;
  f = bface->f;

  /* remove bface from list */
  bp = bface->prev;
  bn = bface->next;
  if(bn) bn->prev = bp;
  if(bp) bp->next = bn;
  else   pat->bfaces[f] = bn;

  /* now free bface */
  free(bface);

  /* remove pointer to this bface in other bface */
  if(obface) obface->obface = NULL;
  return 1 + remove_bface(obface);
}

/* free all bfaces in pat */
void remove_all_bfaces(tPat *pat)
{
  int f;

  for(f=0; f<6; f++)
  {
    tBface *bface0 = pat->bfaces[f];
    tBface *bf, *bft;

    for(bf=bface0; bf;)
    {
      bft = bf;
      bf  = bf->next;
      remove_bface(bft);
    }
  }
}

/* remove all bfaces without bounding rectangle */
void remove_bfaces_without_brct(tPat *pat)
{
  int f;

  for(f=0; f<6; f++)
  {
    tBface *bface0 = pat->bfaces[f];
    tBface *bf, *bft;

    for(bf=bface0; bf;)
    {
      bft = bf;
      bf  = bf->next;
      if(!bft->brct_isset) remove_bface(bft);
    }
  }
}

/* make bounding rectangle large enough to fit the point X[3] */
void expand_bface_to_include_X(tBface *bface, const double X[3])
{
  int f, dir;

  if(!bface) return;

  f = bface->f;
  dir = f/2;

  if(bface->brct_isset)
  {
    /* expand rectangle */
    expand_brct_to_include_X(bface->brct, dir, X, 1);
  }
  else
  {
    expand_brct_to_include_X(bface->brct, dir, X, 0);
    bface->brct_isset = 1;
  }
}

/* expand bface to cover face edges */
void expand_bface_to_pat_bbox(tBface *bface)
{
  tPat *pat;
  tBface *obface;
  int f,dir, d[2], dd, j;
  double A[2],B[2], L[2], C, X[3], x[3], oX[3];
  double rct[4];

  if(!bface) return;

  pat = bface->pat;
  obface = bface->obface;
  f = bface->f;
  dir = f/2;
  d[0] = Dir1_norm(dir);
  d[1] = Dir2_norm(dir);

  /* copy rectangle of bface */
  for(j=0; j<4; j++)  rct[j] = bface->brct[j];

  /* save adjusted rectangle of bface in rct */
  for(dd=0; dd<2; dd++)
  {
    A[dd] = pat->bbox[2*d[dd]];
    B[dd] = pat->bbox[2*d[dd]+1];
    L[dd] = B[dd] - A[dd];

    if(bface->brct[2*dd] - A[dd] < L[dd]/NPOINTS)   rct[2*dd] = A[dd];
    if(B[dd] - bface->brct[2*dd+1] < L[dd]/NPOINTS) rct[2*dd+1] = B[dd];
  }

  /* expand other bface to same point */
  if(obface)
  {
    tPat *opat = obface->pat;
    int ret;

    C = pat->bbox[f];
    for(j=0; j<2; j++)
    {
      switch(dir)
      {
      case 0:
        X[0] = C;
        X[1] = rct[j];
        X[2] = rct[2+j];
        break;
      case 1:
        X[0] = rct[j];
        X[1] = C;
        X[2] = rct[2+j];
        break;
      case 2:
        X[0] = rct[j];
        X[1] = rct[2+j];
        X[2] = C;
        break;
      }
      set_xyz(pat, 0,-1, X, x);
      /* we use p_XYZ_of_xyz not set_XYZ(opat, 0,-1, oX, x); because
         p_XYZ_of_xyz also rounds oX te be inside opat) */
      ret = p_XYZ_of_xyz(opat, oX, x);
      if(ret>=0)
      {
        expand_bface_to_include_X(bface, X);
        expand_bface_to_include_X(obface, oX);
      }
      if(0)
      {
        prdivider(0);
        printf("pat:\n");
        printpatch(pat);
        printCI(pat);
        printf("bface:\n");
        printbface(bface);
        pr3v("X",X); pr3v("x",x); pr3v("oX",oX);
        printf("\n*** looked for oX(x(X)) on f=%d -> ret=%d:\n", f, ret);
        printf("opat:\n");
        printpatch(opat);
        printCI(opat);
        printf("obface:\n");
        printbface(obface);
        errorexit("x is not in opat!!!");
      }
      //pr3v("X",X);pr3v("x",x);pr3v("oX",oX);
    }
  }
  else /* if there is no obface we can just expand the bface */
  {
    /* put expanded rectangle into bface */
    for(j=0; j<4; j++)  bface->brct[j] = rct[j];
  }
}

/* find first bface with a particular op */
tBface *first_bface_with_op(tPat *pat, int op)
{
  tBface *bf;

  forbfaces(pat, bf)
    if(bf->op == op) return bf;

  return NULL;
}

/* find first bface with a particular op and f */
tBface *first_bface_with_op_f(tPat *pat, int op, int f)
{
  tBface *bf;

  for(bf=pat->bfaces[f]; bf; bf=bf->next)
    if(bf->op == op) return bf;

  return NULL;
}

/* find ith bface on face f with an obface, the first one has i=0 */
tBface *ith_bface_on_f_with_obface(tPat *pat, int f, int i)
{
  tBface *bf;
  int k = 0;

  for(bf=pat->bfaces[f]; bf; bf=bf->next)
    if(bf->obface)
    {
      if(k==i) return bf;
      k++;
    }

  return NULL;
}

/* find number of bfaces on face f with an obface */
int nbfaces_on_f_with_obface(tPat *pat, int f)
{
  tBface *bf;
  int k = 0;

  for(bf=pat->bfaces[f]; bf; bf=bf->next)
    if(bf->obface) k++;

  return k;
}

/* find first bface on face f that contains a point C */
tBface *first_bface_containing_point(tPat *pat, int f, double C[2])
{
  tBface *bf;

  forbfacesonface(pat, f, bf)
    if(C_in_brct(bf->brct , C)) return bf;

  return NULL;
}

/* find first obface on other side of face f that contains a point C */
tBface *first_obface_of_bface_containing_point(tPat *pat, int f, double C[2])
{
  tBface *bf, *obf;

  forbfacesonface(pat, f, bf)
    if(C_in_brct(bf->brct , C))
    {
      obf = bf->obface;
      if(obf)
        if(obf->pat) return obf;
    }

  return NULL;
}

/* find bface of the neighbor node nb on the other side of face f,
   that contains a point described by pat,f,C */
tBface *brct_nbbface_of_bface_containing_point(tNode *nb,
                                               tPat *pat, int f, double C[2])
{
  tBface *bf, *obf;

  forbfacesonface(pat, f, bf)
    if(C_in_brct(bf->brct, C))
    {
      obf = bf->obface;
      if(obf)
        if(obf->pat == nb->pat) return obf;
    }

  return NULL;
}

/* find bface of the neighbor node nb on the other side of face f,
   that contains a point described by pat,f,C */
tBface *nbbface_of_bface_containing_point(tNode *nb,
                                          tPat *pat, int f, double C[2])
{
  tBface *bf, *obf;

  forbfacesonface(pat, f, bf)
  {
    obf = bf->obface;
    if(obf)
      if(obf->pat == nb->pat)
      {
        int ofaces[6];
        double oX[3];

        if(facepoint_in_bfacepair(bf, NULL,-1, C, ofaces, oX))
          return obf;
      }
  }

  return NULL;
}


/*************************************************************************/
/* funcs to set bface info */
/*************************************************************************/

/* set bfaces for each box on the mesh with old algorithm.
   This one is not general and fails in many cases. */
int amr_set_all_bfaces(tMesh *mesh)
{
  double Lmin;
  int inclOuterBound = 1;
  int pr=1, p;

  TIMER_START;
  PRFs(":\n");

  /* find pat size L of smallest pat */
  Lmin = smallest_pat_size(mesh);

  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];
    int extface[6];  /* extface[f]=1  means face f is external, i.e. needs BC */
    int f;

    remove_all_bfaces(pat);

    /* find all external faces of pat */
    find_external_faces_of_pat(pat, Lmin, extface, inclOuterBound);
    if(pr)
    {
      printf("external faces on pat%d: ", p);
      for(f=0; f<6; f++) if(extface[f]) printf(" %d", f);
      printf("\n");
    }
    /* set the bfaces for each external face we found */
    for(f=0; f<6; f++)
    {
      if(extface[f])
      {
        set_bfaces_on_patface(pat, Lmin, f);
      }
    }
    if(0) printbfaces(pat);
  }

  if(0)
  {
    printf("------before expand_bfaces_to_patch_edges:\n");
    printallbfaces(mesh);
    printf("------\n");
  }

  /* now expand bfaces to cover the edges as well */
  expand_bfaces_to_patch_edges(mesh);

  /* if there is no other patch on the other side of a bface, mark it as
     outer boundary */
  mark_all_bfaces_without_op_as_OUTERBOUND(mesh);

//  /* set ofi and bit fields for all bfaces */
//  if(pr)
//    printf("Coordinates_set_bfaces_oldWT:  setting ofi and bit fields for all bfaces\n");
//  set_ofi_in_all_bfaces(mesh);
//  set_bits_in_all_bfaces(mesh);
//
//  /* fix some issues */
//  if(remove_bfacepoints_that_cause_inconsistent_touch_bits(mesh))
//    set_bits_in_all_bfaces(mesh);


  /* this sets the setnormalderiv again */
  set_consistent_flags_in_all_bfaces(mesh);

//  /* set some var indices if we need to interpolate */
//  set_oXi_oYi_oZi_in_all_bfaces(mesh);

  TIMER_STOP;
  return 0;
}

/* determine which patch faces touch or overlap other patches,
   or touch the outside of our pat, returns extface */
/* extface[f]=1 means face f is external, i.e. needs BC */
/* Faces with e.g. periodic coordinates are marked as extface[f]=0. */
/* If inclOuterBound=1 we mark faces that are not in contact with any other
   pat as external (because they probably need an outer BC) */
void find_external_faces_of_pat(tPat *pat, double L,
                                int *extface, int inclOuterBound)
{
  tMesh *mesh = pat->mesh;
  int n[] = { NPOINTS,NPOINTS,NPOINTS };     /* we use NPOINTS points */
  double X0[3], LX[3], dX[3]; /* grid of points */
  int nd, dir, dd;
  int i,j,k, pl, f;
  intList *opl; /* list that contains other patches */
  int op;
  double oX[3];

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

  /* list that will contain other patches */
  opl = alloc_intList();

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

    /* go over faces, but not edges, because our grid doesn't contain edges */
    for(pl=0; pl<nd; pl+=nd-1)
    {
      op = -2; /* -2 means op is not set */

      /* pick face index */
      if(pl==0) f=2*dir;
      else      f=2*dir+1;

      /* do nothing if not an external face */
      if(!extface[f]) continue;

      /* look for points in other patches just outside this pat */
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
        set_xyz(pat, NULL,-1, X, x);

        /* use normal vector to find point ox slightly outside pat */
        patch_normal_at_XYZ(pat, f, X, Ndir);
        for(d0=0., dd=0; dd<3; dd++) d0 += x[dd]*x[dd];
        d0 = sqrt(d0);
        for(dd=0; dd<3; dd++)
        {
          dx[dd] = Ndir[dd]*(L*EPS1 + d0*EPS2);
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

        /* find point in other patches */
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


/* find and set all bfaces on an external pat face f */
/* The idea is to loop over the pat faces and then move out using the Cartesian
   normal vectors. Then we check if we are in an other pat. */
/* This only creates a bface if there is some contact with another pat,
   or if we seem to be at an other boundary.
   It returns the number of new bfaces made for face f in this pat. So if it
   returns 0 nothing was done at all. */
int set_bfaces_on_patface(tPat *pat, double L, int f)
{
  tMesh *mesh = pat->mesh;
  int p = pat->p;
  int dir = f/2;
  tBface *bface, *obface;
  tPat *opat;
  int n[] = { NPOINTS,NPOINTS,NPOINTS };     /* we use NPOINTS points */
  double X0[3], LX[3], dX[3]; /* grid of points */
  int dd;
  int i,j,k, plane, li, ret0, ret, of;
  int face[6];
  intList *opl = alloc_intList(); /* list that contains other patches*/
  int op, nbfaces;
  double oX[3];

  /* make opl that contains all patches except p,
     and add one bface for each of the other patches */
  forpatches(mesh, op) if(op!=p)
  {
    opat = mesh->pat[op];

    /* if pat and other pat are Cartesian, we can know already if
       if they touch or not */
    if( (pat->dXYZ_dxyz==NULL) && (opat->dXYZ_dxyz==NULL) )
    {
      double bb[6];
      int touch = touch_or_intersect_bb1_bb2(pat->bbox, opat->bbox, bb);
      /* skip patch op if they do not touch */
      if(!touch) continue;
    }

    /* we want one bface for each other pat */
    bface = first_bface_with_op_f(pat, op, f);
    /* add empty bface for each other pat where we don't have one already */
    if(!bface)
    {
      bface = add_empty_bface(pat, f);
      bface->op = op;
    }
    unionpush_intList(opl, op);
  }
  /* add one more empty bface for outer boundary points that
     are not in contact with any other pat */
  bface = add_empty_bface(pat, f);
  bface->op = -1;

  /* make a grid of points, that excludes endpoints */
  for(dd=0; dd<3; dd++)
  {
    X0[dd] = pat->bbox[2*dd];
    LX[dd] = pat->bbox[2*dd+1] - X0[dd];
    dX[dd] = LX[dd]/(n[dd]);
    X0[dd] += dX[dd] * 0.5;
  }

  /* look for points in other patches just outside this pat */
  plane = (n[dir] - 1) * (f%2);
  forplaneN(dir, i,j,k, n, plane)
  {
    double X[3], x[3];
    double N[3];
    double ox[3], dx[3], d0;

    /* point grid, that never includes edges */
    X[0] = X0[0] + dX[0] * i;
    X[1] = X0[1] + dX[1] * j;
    X[2] = X0[2] + dX[2] * k;

    /* pick one of X,Y,Z on boundary */
    X[dir] = pat->bbox[f];

    /* get x,y,z of point from which we move out, sometimes this is
       one point in from the edge */
    set_xyz(pat, NULL,-1, X, x);

    /* get outward vector N[0],N[1],N[2] */
    patch_normal_at_XYZ(pat, f, X, N);

    /* use vector N to find point ox,oy,oz slightly outside pat */
    for(d0=0., dd=0; dd<3; dd++) d0 += x[dd]*x[dd];
    d0 = sqrt(d0);
    for(dd=0; dd<3; dd++)
    {
      dx[dd] = N[dd]*(L*EPS1 + d0*EPS2);
      ox[dd] = x[dd] + dx[dd];
    }

    /* mark other pat as non-existent by default */
    op = -1;
    opat = NULL;

    /* find point in other patches */
    for(li=0; li<opl->n; li++)
    {
      int bi = opl->e[li];
      opat = mesh->pat[bi];

      /* check if point is in opat */
      ret = p_XYZ_of_xyz(opat, oX, ox);
      //printf("p=%d f=%d: opat->p=%d", p,f,opat->p);
      //prbbox(opat->bbox,3);pr3v("X",X);
      //pr3v("ox",ox);pr3v("oX",oX);
      //printf(": ret=%d\n", ret);
      if(ret>=0)
        for(dd=0; dd<3; dd++)
        {
          if(!(pat->periodic[dd]))
            if(oX[dd] < opat->bbox[2*dd] || opat->bbox[2*dd+1] < oX[dd])
              ret=-1;
        }
      if(ret>=0)
      {
        op = opat->p;
        break; /* we found point in pat op */
      }
    }

    /* the bface that should have this point has other patch op and face f */
    bface = first_bface_with_op_f(pat, op, f);

    /* now try to find this point also in the other bface */
    if(op>=0)
    {
      /* if we get here the point is also in the other patch opat */

      /* get ox very close to face of opat and recalc oX */
      for(dd=0; dd<3; dd++)  ox[dd] = x[dd] + dx[dd]*1e-8;
      ret0 = p_XYZ_of_xyz(opat, oX, ox);
      //printf("ret0=%d oX[2]=%.15g\n", ret0, oX[2]);

      /* find point in other bface and see what face it is on */
      ret = XYZ_on_face(opat, face, oX);
      /* if point is only on one face of the other patch, we add it
         to the other bface */
      if(ret==1)
      {
        for(of=0; face[of]==0; of++) ;
        obface = first_bface_with_op_f(opat, p, of);
        if(!obface)
        {
          obface = add_empty_bface(opat, of);
          obface->op = p;
        }

        expand_bface_to_include_X(obface, oX);

        /* link bface and obface */
        if(!bface->obface)
          bface->obface = obface;
        if(bface->obface != obface)
        {
          printf("bface=%p  bface->obface=%p\n", bface, bface->obface);
          printbface(bface);
          printf("obface=%p  obface->obface=%p\n", obface, obface->obface);
          printbface(obface);
          errorexit("what happened???");
        }
        if(!obface->obface)
          obface->obface = bface;
        if(obface->obface != bface)
        {
          printf("bface=%p  bface->obface=%p\n", bface, bface->obface);
          printbface(bface);
          printf("obface=%p  obface->obface=%p\n", obface, obface->obface);
          printbface(obface);
          errorexit("what happened???");
        }
      }
      else if(ret>1)
      {
        continue; /* do not use this point as it lies on two faces of opat */
      }
      if(ret==0)
      {
        int fi;
        printf("dir=%d plane=%d i,j,k=%d,%d,%d ret0=%d\n",
               dir, plane, i,j,k, ret0);
        printpatch(pat);
        printCI(pat);
        printbface(bface);
        pr3v("X", X);
        pr3v("x", x);
        printf("\n");
        printpatch(opat);
        printCI(opat);
        pr3v("ox", ox);
        pr3v("oX", oX);
        printf("\n");
        printf("face=");
        for(fi=0; fi<6; fi++) printf("%d ", face[fi]);
        printf("-> ret=%d\n", ret);
        errorexiti("oX was supposed to be on 1 face, not %d faces!!!", ret);
      }
    }

    /* add point to the bface with the correct op and f */
    expand_bface_to_include_X(bface, X);
  }

  /* remove empty bfaces */
  remove_bfaces_without_brct(pat);

  /* count num of bfaces */
  nbfaces = 0;
  for(bface=pat->bfaces[f]; bface; bface=bface->next) nbfaces++;

  free_intList(opl);

  return nbfaces;
}


/* find point pat,f,C in opat, computes ofaces and oX,
   returns -1 if point is not in opat, otherwise it returns
   the number of opat faces the point is on */
int facepoint_in_opat(tPat *pat, int f, double C[2],
                      tPat *opat, int ofaces[6], double oX[3])
{
  int dir = f/2;
  int od0 = Dir1_norm(dir);
  int od1 = Dir2_norm(dir);
  int op, j, nofaces;
  double X[3], x[3];

  /* pick one of X,Y,Z on pat boundary, and the other two are from C */
  X[dir] = pat->bbox[f];
  X[od0] = C[0];
  X[od1] = C[1];

  /* find X in Cart coords */
  set_xyz(pat, NULL,-1, X, x);

  /* check if point is in opat */
  op = p_XYZ_of_xyz(opat, oX, x);

  /* now try to find this point also on a face of opat */
  if(op>=0)
    nofaces = XYZ_on_face(opat, ofaces, oX);
  else
  {
    for(j=0; j<6; j++) ofaces[j] = 0;
    nofaces = -1;
  }

  return nofaces;
}

/* check if the point bface,node,ijk,C is in a bface pair,
   returns 0 if no, or number of other faces if yes,
   it also writes the point in opat coords into oX */
int facepoint_in_bfacepair(tBface *bface, tNode *node, int ijk, double C[2],
                           int ofaces[6], double oX[3])
{
  tMesh *mesh;
  tPat *pat, *opat;
  tBface *obface;
  int f, op, of;

  if(!bface) return 0;

  pat = bface->pat;
  mesh = pat->mesh;
  f = bface->f;
  obface = bface->obface;

  if(!obface) return 0;

  opat = obface->pat;
  of = obface->f;

  /* if we do not have a grid point, transform from C to x to oX */
  if(node==NULL || ijk<0)
    facepoint_in_opat(pat,f,C, opat,ofaces,oX);
  else
  {
    int ix = Ind("x");
    double x[] = { Vard(node,ix)[ijk],
                   Vard(node,ix+1)[ijk], Vard(node,ix+2)[ijk] };

    /* check if point is in opat */
    op = p_XYZ_of_xyz(opat, oX, x);

    /* now try to find this point also on a face of opat */
    if(op>=0)
      XYZ_on_face(opat, ofaces, oX);
    else
      ofaces[of] = 0;
  }
  return ofaces[of];
}


/* expand all bfaces to cover patch edges */
void expand_bfaces_to_patch_edges(tMesh *mesh)
{
  int p;

  /* now expand bfaces to cover the edges as well */
  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];
    tBface *bface;

    forbfaces(pat, bface)
      expand_bface_to_pat_bbox(bface);
  }
}

/* mark all bfaces that have no other patch as outer boundary */
void mark_all_bfaces_without_op_as_OUTERBOUND(tMesh *mesh)
{
  int p;

  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];
    tBface *bface;

    forbfaces(pat, bface)
    {
      /* if op=-1 there is no other box, mark as outer boundary */
      if(bface->op == -1)
        bface->boundary = OUTERBOUND;
    }
  }
}

/* make sure bit fields in all bfaces are consitent.
   Right now we just set bface->face2, which has the same meaning as sgrid's
   setnormalderiv */
int set_consistent_flags_in_all_bfaces(tMesh *mesh)
{
  int p;
  int bface_options = Par("bface_options");
  int forder1 = Getv(bface_options,"face2_order1");
  int forder2 = Getv(bface_options,"face2_order2");
  int forder3 = Getv(bface_options,"face2_order3");
  int forder4 = Getv(bface_options,"face2_order4");

  /* face2_order4 implies face2_order3 */
  forder3 = forder3 || forder4;

  /* set face2=0 in all bfaces if we want a particluar order */
  if(forder1 || forder2 || forder3)
    zero_face2_flag_in_all_bfaces(mesh);

  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];
    tBface *bface;

    forbfaces(pat, bface)
    {
      tBface *obface = bface->obface;
      //tPat *opat;
      //int op  = bface->op;
      //int f = bface->f;
      //int of;

       /* do nothing if there no other bface */
      if(!obface) continue;

      /* other pat on corresponding bface */
      //opat = obface->pat;
      //of = obface->f;

      /* unlike in sgrid, the two bfaces are always touching */
      if(obface->obface == bface) /* we have 2 paired bfaces */
      {
        /* set consistent face2 flag */
        if(forder2)
        {
        /* note forder2 makes templates_GMRES_with_BlockJacobi_precon fail
             with 6 or more cubed spheres */
          if(obface->face2 == 0) bface->face2 = 1;
          else                   bface->face2 = 0;
        }
        else /* use forder1 */
        {
          if(bface->face2 == 0) obface->face2 = 1;
          else                  obface->face2 = 0;
        }
      }
      else /* obface doesn't refer to this bface */
      {
        if(obface->face2 == 0) bface->face2 = 1;
        else                   bface->face2 = 0;
      }
    } /* end forbfaces */
  }
  /* Note: forder3 it supposed to fix failures in
     templates_GMRES_with_BlockJacobi_precon. It fails if there is a Neumann
     BC on all faces of a pat. This can happen on pat0 with forder2.
     Or it can happen on pat5&6 with forder1 if the BC on face1 of pat5
     or 6 is also a Neumann BC. The corresponding block is then a singular
     matrix, as can be seen with GridIterators_verbose = yes */
  if(forder3) toggle_face2_flag_in_faces4_5_of_cubes(mesh);

  /* Note: forder4 it supposed to fix more failures in
     templates_GMRES_with_BlockJacobi_precon. It fails for an odd number of
     points in patches 11,12 and 24,25, i.e. CubSph domains 4,5. Presumably
     this happens because even with forder3 dom 4 and 5 have Neuman BCs
     everywhere except for face1. */
  if(forder4) toggle_face2_flag_of_CubSph_doms_0_4_and_1_5(mesh);

  return 0;
}

/* set face2=0 in all bfaces */
int zero_face2_flag_in_all_bfaces(tMesh *mesh)
{
  int p;

  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];
    tBface *bface;

    forbfaces(pat, bface) bface->face2 = 0;
  }
  return 0;
}

/* Toggle face2 flag of face 4 & 5 in all Cartesian cubes.
   We want to avoid that BCs are the same on all sides of the cube. */
int toggle_face2_flag_in_faces4_5_of_cubes(tMesh *mesh)
{
  int p;

  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];
    tBface *bface;

    /* look at cubes only */
    if(pat->CI->type == 0)
    {
      /* loop over bfaces */
      forbfaces(pat, bface)
      {
        int f = bface->f;
        //int op  = bface->op;
        //tPat *opat;
        tBface *obface = bface->obface;

        /* do something only on face 4 and 5 */
        if((f == 4) || (f == 5))
        {
          if(bface->face2) bface->face2 = 0;
          else             bface->face2 = 1;
        }
        else /* otherwise go to next bface */
          continue;

        /* do nothing more if there is no other bface */
        if(!obface) continue;

        /* other pat from corresponding bface */
        //opat = obface->pat;

        /* the 2 are always touching */
        if(obface->obface == bface) /* we have 2 paired bfaces */
        {
          /* set consistent face2 flag */
          if(bface->face2 == 0) obface->face2 = 1;
          else                  obface->face2 = 0;
        }
      } /* end forbfaces */
    } /* end cube case */
  }
  return 0;
}

/* Toggle face2 flag between pat->CI->dom=0 and pat->CI->dom=4
   patches, as well between pat->CI->dom=1 and pat->CI->dom=5 patches in
   all types of cubed spheres.
   We want to avoid too many Neuman BCs in dom 4 and 5. */
int toggle_face2_flag_of_CubSph_doms_0_4_and_1_5(tMesh *mesh)
{
  int p;

  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];
    tBface *bface;
    int domain = pat->CI->dom;

    /* look at cubed spheres only */
    if( (0 < pat->CI->type) && (pat->CI->type <= CubedShell) )
    {
      /* loop over bfaces */
      forbfaces(pat, bface)
      {
        //int f = bface->f;
        //int op  = bface->op;
        tPat *opat = NULL;
        tBface *obface = bface->obface;
        int odomain = -1;

        /* do nothing more if there is no other bface */
        if(!obface) continue;

        /* other pat from corresponding bface */
        opat = obface->pat;
        odomain = opat->CI->dom;

        /* do something only if domain,odomain is 0,4 or 1,5 */
        if((domain == 0 && odomain == 4) || (domain == 1 && odomain == 5))
        {
          if(bface->face2) bface->face2 = 0;
          else             bface->face2 = 1;
        }
        else /* otherwise go to next bface */
          continue;

        /* the 2 are always touching */
        if(obface->obface == bface) /* we have 2 paired bfaces */
        {
          /* set consistent face2 flag */
          if(bface->face2 == 0) obface->face2 = 1;
          else                  obface->face2 = 0;
        }
      } /* end forbfaces */
    } /* end cubed sph. case */
  }
  return 0;
}


/*************************************************************************/
/* funcs to add and remove nfaces */
/*************************************************************************/

/* get mem for 1 nface */
tNface *alloc_nface(tNode *node, int f)
{
  tNface *nface = calloc(1, sizeof(nface[0]));
  if(!nface) errorexit("not enough memory");
  nface->node = node;
  nface->f    = f;
  return nface;
}

/* add a nface to a node, f denotes the node face (0 to 5),
   return the new nface */
tNface *add_empty_nface(tNode *node, int f)
{
  //tMesh *mesh = node->mesh;
  tNface *nface = alloc_nface(node, f);
  tNface *nface0 = node->nfaces[f];
  tNface *nf;

  /* add nface to end of list in node */
  if(nface0)
  {
    for(nf=nface0; nf->next; nf=nf->next) ;
    nf->next = nface;
    nface->prev = nf;
  }
  else
    node->nfaces[f] = nface;

  return nface;
}

/* add a pair of node + neighbor that touch via face f and nb_f */
int add_nface(tNode *node, int f, tNode *nb, int nb_f)
{
  tNface *nface  = add_empty_nface(node, f);
  tNface *onface = add_empty_nface(nb, nb_f);

  /* set info in new nface in node */
  nface->node   = node;
  nface->f      = f;
  nface->onface = onface;

  /* set info in new nface in nb */
  onface->node   = nb;
  onface->f      = nb_f;
  onface->onface = nface;

//char s[100];
//printf("add: node %s f=%d,  ", nodename(node,s,99), f);
//printf("nb %s nb_f=%d\n", nodename(nb,s,99), nb_f);
//printnface(nface);
  return 0;
}


/* remove a nface from a node, return number of nfaces removed */
int remove_nface(tNface *nface)
{
  tNode *node;
  tNface *onface, *nxt, *prv;
  int f;

  if(!nface) return 0;

  node = nface->node;
  onface = nface->onface;
  f = nface->f;

  /* remove nface from list */
  prv = nface->prev;
  nxt = nface->next;
  if(nxt) nxt->prev = prv;
  if(prv) prv->next = nxt;
  else    node->nfaces[f] = nxt;

  /* now free nface */
  free(nface);

//char s[100];
//printf("rem: node %s\n", nodename(node,s,99));

  /* remove pointer to this nface in other nface */
  if(onface) onface->onface = NULL;
  return 1 + remove_nface(onface);
}

/* free all nfaces in node */
void remove_all_nfaces(tNode *node)
{
  int f;

  if(!node) return;

//PRFs(" Begin\n");printnfaces(node);

  for(f=0; f<6; f++)
  {
    tNface *nface0 = node->nfaces[f];
    tNface *nf, *nft;

    for(nf=nface0; nf;)
    {
      nft = nf;
      nf  = nf->next;
      remove_nface(nft);
    }
  }
//PRFs(" End\n");printnfaces(node);
}



/* find out if any node points on face f are on face nb_f of node nb */
int find_nodefacepoints_in_nbface(tNode *node, int f, tNode *nb, int nb_f)
{
  int dir = f/2;
  int n[] = { 3,3,3 };        /* we use 3 points */
  double X0[3], LX[3], dX[3]; /* grid of points */
  int dd;
  int i,j,k, plane, ret0, ret;

  /* make a grid of points, that excludes endpoints */
  for(dd=0; dd<3; dd++)
  {
    X0[dd] = node->bbox[2*dd];
    LX[dd] = node->bbox[2*dd+1] - X0[dd];
    dX[dd] = LX[dd]/(n[dd]);
    X0[dd] += dX[dd] * 0.5;
  }

  /* loop over points */
  plane = (n[dir] - 1) * (f%2);
  forplaneN(dir, i,j,k, n, plane)
  {
    double X[3], x[3], oX[3];
    int nbface[6];

    /* point grid, that never includes edges */
    X[0] = X0[0] + dX[0] * i;
    X[1] = X0[1] + dX[1] * j;
    X[2] = X0[2] + dX[2] * k;

    /* pick one of X,Y,Z on boundary */
    X[dir] = node->bbox[f];

    /* get x,y,z of X,Y,Z and then oX,oY,oZ in nb */
    set_xyz(NULL, node,-1, X, x);
    ret0 = l_XYZ_of_xyz(nb,-1, oX, x);

    /* try another point, if this one is not in nb */
    if(ret0<0) continue;

    /* check if this point is on nb face */
    ret = XYZ_on_face(nb->pat, nbface, oX);

    /* if this point is only in face nb_f of nb we are done */
    if(ret==1 && nbface[nb_f]) return 1;

    /* if this point is in several faces try another point */
    if(ret>1) continue;

    if(ret==0)
    {
      errorexiti("oX was supposed to be on 1 face, not %d faces!!!", ret);
    }
  }

  return 0;
}

/* check if node and nb has common points on faces f and nb_f
   since res in find_nodefacepoints_in_nbface is low, try it both ways */
int common_facepoints(tNode *node, int f, tNode *nb, int nb_f)
{
  int f1, f2;

  f1 = find_nodefacepoints_in_nbface(node,f, nb,nb_f);
  if(f1) return 1;

  f2 = find_nodefacepoints_in_nbface(nb,nb_f, node,f);

  return f2 || f1;
}
