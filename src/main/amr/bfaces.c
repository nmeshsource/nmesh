/* bfaces.c */
/* Wolfgang Tichy, 3/2019 */

#include "nmesh.h"
#include "amr.h"

#define EPS 1e-5
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
  tBface *bface = alloc_bface(pat, f);
  tBface *bface0 = pat->bface0;
  tBface *bf;

  /* add bface to end of list in pat */
  if(bface0)
  {
    for(bf=bface0; bf->next; bf=bf->next) ;
    bf->next = bface;
    bface->prev = bf;
  }
  else
    pat->bface0 = bface;

  /* set some bface info */
  bface->op = -1; /* other patch not known yet */
  bface->ioX[0] = -1; /* var indices of other coords not known yet */
  bface->ioX[1] = -1; /* var indices of other coords not known yet */
  bface->ioX[2] = -1; /* var indices of other coords not known yet */
  return bface;
}

/* remove a bface from a pat, return number of bfaces removed */
int remove_bface(tBface *bface)
{
  tPat *pat;
  tBface *obface, *bn, *bp;

  if(!bface) return 0;

  pat = bface->pat;
  obface = bface->obface;

  /* remove bface from list */
  bp = bface->prev;
  bn = bface->next;
  if(bn) bn->prev = bp;
  if(bp) bp->next = bn;
  else   pat->bface0 = bn;

  /* now free bface */
  free(bface);

  /* remove pointer to this bface in other bface */
  if(obface) obface->obface = NULL;
  return 1 + remove_bface(obface);
}

/* free all bfaces in pat */
void remove_all_bfaces(tPat *pat)
{
  tBface *bface0 = pat->bface0;
  tBface *bf, *bft;

  for(bf=bface0; bf;)
  {
    bft = bf;
    bf  = bf->next;
    remove_bface(bft);
  }
}

/* remove all bfaces without bounding rectangle */
void remove_bfaces_without_brct(tPat *pat)
{
  tBface *bface0 = pat->bface0;
  tBface *bf, *bft;

  for(bf=bface0; bf;)
  {
    bft = bf;
    bf  = bf->next;
    if(!bft->brct_isset) remove_bface(bft);
  }
}

/* make bounding rectangle large enough to fit the point X[3] */
void expand_bface_to_include_X(tBface *bface, const double X[3])
{
  int f, dir, d;
  double C[2]; /* point coords in face */

  if(!bface) return;

  f = bface->f;
  dir = f/2;

  switch(dir)
  {
  case 0:
    C[0] = X[1];  C[1] = X[2];
    break;
  case 1:
    C[0] = X[0];  C[1] = X[2];
    break;
  case 2:
    C[0] = X[0];  C[1] = X[1];
    break;
  }
  if(bface->brct_isset)
  {
    /* expand rectangle */
    for(d=0; d<2; d++)
    {
      if(C[d] < bface->brct[2*d])   bface->brct[2*d]   = C[d];
      if(C[d] > bface->brct[2*d+1]) bface->brct[2*d+1] = C[d];
    }
  }
  else
  {
    bface->brct[1] = bface->brct[0] = C[0];
    bface->brct[3] = bface->brct[2] = C[1];
    bface->brct_isset = 1;
   } 
}

/* expand bface to cover face edges */
void expand_bface_to_pat_bbox(tBface *bface)
{
  tPat *pat;
  tBface *obface;
  int f,dir, d[2], dd;
  double A[2],B[2], L[2], C, X[3], x[3], oX[3];

  if(!bface) return;

  pat = bface->pat;
  obface = bface->obface;
  f = bface->f;
  dir = f/2;
  d[0] = Dir1_norm(dir);
  d[1] = Dir2_norm(dir);

  /* adjust rectangle of bface */
  for(dd=0; dd<2; dd++)
  {
    A[dd] = pat->bbox[2*d[dd]];
    B[dd] = pat->bbox[2*d[dd]+1];
    L[dd] = B[dd] - A[dd];

    if(bface->brct[2*dd] - A[dd] < L[dd]/NPOINTS)
      bface->brct[2*dd] = A[dd];

    if(B[dd] - bface->brct[2*dd+1] < L[dd]/NPOINTS)
      bface->brct[2*dd+1] = B[dd];
  }

  /* expand other bface to same point */
  if(obface)
  {
    tPat *opat = obface->pat;
    int j;

    C = pat->bbox[f];
    for(j=0; j<2; j++)
    {
      switch(dir)
      {
      case 0:
        X[0] = C;
        X[1] = bface->brct[j];
        X[2] = bface->brct[2+j];
        break;
      case 1:
        X[0] = bface->brct[j];
        X[1] = C;
        X[2] = bface->brct[2+j];
        break;
      case 2:
        X[0] = bface->brct[j];
        X[1] = bface->brct[2+j];
        X[2] = C;
        break;
      }
      set_xyz(pat, 0,-1, X, x);
      set_XYZ(opat, 0,-1, oX, x);
      //pr3v("X",X);pr3v("x",x);pr3v("oX",oX);
      expand_bface_to_include_X(obface, oX);
    }
  }
}

/* find first bface with a particular op */
tBface *first_bface_with_op(tPat *pat, int op)
{
  tBface *bf;

  for(bf=pat->bface0; bf; bf=bf->next)
    if(bf->op == op) return bf;

  return NULL;
}

/* find first bface with a particular op and f */
tBface *first_bface_with_op_f(tPat *pat, int op, int f)
{
  tBface *bf;

  for(bf=pat->bface0; bf; bf=bf->next)
    if(bf->op == op && bf->f == f) return bf;

  return NULL;
}

/*************************************************************************/
/* funcs to set bface info */
/*************************************************************************/

/* set bfaces for each box on the mesh with old algorithm.
   This one is not general and fails in many cases. */
int amr_set_all_bfaces(tMesh *mesh)
{
  int inclOuterBound = 1;
  int pr=1, p;

  PRFs(":\n");

  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];
    int extface[6];  /* extface[f]=1  means face f is external, i.e. needs BC */
    int f;

    remove_all_bfaces(pat);

    /* find all external faces of pat */
    find_external_faces_of_pat(pat, extface, inclOuterBound);
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
        set_bfaces_on_patface(pat, f);
      }
    }
    if(0) printbfaces(pat);
  }

  /* now expand bfaces to cover the edges as well */
  forpatches(mesh, p)
  {
    tPat *pat = mesh->pat[p];
    tBface *bface;

    forbfaces(pat, bface)
      expand_bface_to_pat_bbox(bface);
  }

//expand_bface_to_pat_bbox(mesh->pat[1]->bface0->next);

//  /* set ofi and bit fields for all bfaces */
//  if(pr)
//    printf("Coordinates_set_bfaces_oldWT:  setting ofi and bit fields for all bfaces\n");
//  set_ofi_in_all_bfaces(mesh);
//  set_bits_in_all_bfaces(mesh);
//
//  /* fix some issues */
//  if(remove_bfacepoints_that_cause_inconsistent_touch_bits(mesh))
//    set_bits_in_all_bfaces(mesh);
//
//  /* set outer boundary flag */
//  mark_all_bfaces_without_ob_as_outerbound(mesh);
//
//  /* this sets the setnormalderiv again */
//  set_consistent_flags_in_all_bfaces(mesh);
//
//  /* set some var indices if we need to interpolate */
//  set_oXi_oYi_oZi_in_all_bfaces(mesh);

  return 0;
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
  int n[] = { NPOINTS,NPOINTS,NPOINTS };     /* we use NPOINTS points */
  double X0[3], LX[3], dX[3]; /* grid of points */
  int nd, dir, dd;
  int i,j,k, pl, f;
  intList *opl; /* list that contains other patches */
  int op;
  double oX[3];
  double L;

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
        set_xyz(pat, NULL,-1, X, x);

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


/* find and set all bfaces on an external pat face f */
/* The idea is to loop over the pat faces and then move out using the Cartesian
   normal vectors. Then we check if we are in an other pat. */
/* This only creates a bface if there is some contact with another pat,
   or if we seemto be at an other boundary.
   It returns the number of new bfaces made for face f in this pat. So if it
   returns 0 nothing was done at all. */
int set_bfaces_on_patface(tPat *pat, int f)
{
  tMesh *mesh = pat->mesh;
  int p = pat->p;
  int dir = f/2;
  tBface *bface, *obface;
  tPat *opat;
  int n[] = { NPOINTS,NPOINTS,NPOINTS };     /* we use NPOINTS points */
  double X0[3], LX[3], dX[3]; /* grid of points */
  int dd;
  int i,j,k, plane, li, ret, of;
  int face[6];
  intList *opl = alloc_intList(); /* list that contains other pates*/
  int op, nbfaces;
  double oX[3];
  double L;

  /* find pat size L of smallest pat */
  L = smallest_pat_size(mesh);

  /* make opl that contains all patches except p,
     and add one bface for each of the other pates */
  forpatches(mesh, op) if(op!=p)
  {
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

  /* look for points in other pates just outside this pat */
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
      dx[dd] = N[dd]*(L+d0)*EPS;
      ox[dd] = x[dd] + dx[dd];
    }

    /* mark other pat as non-existent by default */
    op = -1;

    /* find point in other patches */
    for(li=0; li<opl->n; li++)
    {
      int bi = opl->e[li];
      opat = mesh->pat[bi];

      /* check if point is in opat */
      ret = p_XYZ_of_xyz(opat, oX, ox);
      //printf("p=%d f=%d: opat->p=%d", p,f,opat->p);
      //prbbox(opat->bbox,3);pr3v("X",X);
      //pr3v("oX",oX);pr3v("ox",ox);
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

    /* add point to the bface with the correct op and f */
    bface = first_bface_with_op_f(pat, op, f);
    expand_bface_to_include_X(bface, X);

    /* now try to put this also in the other bface */
    if(op>=0)
    {
      /* get ox very close to face of opat and recalc oX */
      for(dd=0; dd<3; dd++)  ox[dd] = x[dd] + dx[dd]*1e-8;
      ret = p_XYZ_of_xyz(opat, oX, ox);
      ret = XYZ_on_face(opat, face, oX);
      if(ret!=1) errorexit("oX was supposed to be on one face...");
      for(of=0; face[of]==0; of++) ;
      obface = first_bface_with_op_f(opat, p, of);
      if(!obface)
      {
        obface = add_empty_bface(opat, of);
        obface->op = p;
      }
      expand_bface_to_include_X(obface, oX);

      /* link bface and obface */
      if(!bface->obface)          bface->obface = obface;
      if(bface->obface != obface) errorexit("what happened???");
      if(!obface->obface)         obface->obface = bface;
      if(obface->obface != bface) errorexit("what happened???");
    }
  }

  /* remove empty bfaces */
  remove_bfaces_without_brct(pat);

  /* count num of bfaces */
  nbfaces = 0;
  for(bface=pat->bface0; bface; bface=bface->next) nbfaces++;

  free_intList(opl);

  return nbfaces;
}
