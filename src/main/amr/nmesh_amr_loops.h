/* nmesh_amr_loops.h */
/* Wolfgang Tichy, 1/2019 */


/* decide if we use "omp taskloop" of OpenMP 4.5 from Nov. 2015 */
#if OMP_VERSION >= 201511
#define USE_OMP_TASKLOOP
#endif
// for now we do not use "omp taskloop" because it messes up ham & mom
#undef USE_OMP_TASKLOOP

/****************************************************************************/
/* Loops are performed by macros so that the user has to know very little
   about the implementation details. */
/****************************************************************************/

/****************************************************************************/
/* loops that should be used in most modules                                */
/****************************************************************************/

/* loop over my leaf nodes on this proc without OpenMP */
#define formyelms_noomp(mesh) \
  for(int ei_=0; ei_ < mesh->nmyelm; ei_++)

/* for easy access */
#define MyID   ei_
#define MyElm  mesh->myelm[ei_]
#define MyElm0 mesh->myelm[0]
#define Elm_MyID(mesh, myid)  mesh->myelm[myid]

#define formyelms_s_n(mesh, s, n) \
  for(int ei_=s; (ei_ < mesh->nmyelm) && (ei_ < s+n); ei_++)

/* do we use omp taskloop? */
#ifdef USE_OMP_TASKLOOP

/* we use OpenMP to parallelize the 2nd loop in formyelms_noomp */
#define formyelms(mesh) \
  NODELEVEL_Pragma(omp parallel) \
  NODELEVEL_Pragma(omp master) \
  NODELEVEL_Pragma(omp taskloop) \
  for(int ei_=0; ei_ < mesh->nmyelm; ei_++)

/* this one will have no "parallel" on its own */
/*
#define formyelms_ompfor(mesh) \
  NODELEVEL_Pragma(omp master) \
  NODELEVEL_Pragma(omp taskloop) \
  for(int ei_=0; ei_ < mesh->nmyelm; ei_++)
*/
/* THIS DIDN'T WORK in dg_add_surface_fluxes!!! It caused a race for e.g.
   dgi->node = node; ... It seems tasks get processed by arbitrary threads so
   that allocating dgi on a per thread basis is not good enough... */
/* to start tasks formyelms_ompfor has to be inside a:
   #pragma omp parallel {  } */

#else

/* we use OpenMP to parallelize the 2nd loop in formyelms_noomp */
#define formyelms(mesh) \
  NODELEVEL_Pragma(omp parallel for) \
  for(int ei_=0; ei_ < mesh->nmyelm; ei_++)
#endif


/* this one will have no "parallel" on its own */
#define formyelms_ompfor(mesh) \
  NODELEVEL_Pragma(omp for) \
  for(int ei_=0; ei_ < mesh->nmyelm; ei_++)
/* to start tasks formyelms_ompfor has to be inside a:
   #pragma omp parallel {  } */


/* for compatibility */
#define formylnodes_noomp(mesh)  formyelms_noomp(mesh)
#define formylnodes(mesh)        formyelms(mesh)
#define formylnodes_ompfor(mesh) formyelms_ompfor(mesh)
#define MyLnode                  MyElm
#define MyLnode0                 MyElm0
#define Lnode_myid(mesh, myid)   Elm_MyID(mesh, myid)


/****************************************************************************/
/* Old node macros. Do not use!!! */
/****************************************************************************/

/* we use OpenMP to parallelize the 2nd loop in formylnodes_noomp */
#define old_taskloop_formylnodes(mesh) \
  for(int cat_=0; cat_ < mesh->myln->nncats; cat_++) \
  NODELEVEL_Pragma(omp parallel) \
  NODELEVEL_Pragma(omp master) \
  NODELEVEL_Pragma(omp taskloop) \
  for(int li_=0; li_ < mesh->myln->ncat[cat_]; li_++)

/* we use OpenMP to parallelize the 2nd loop in formylnodes_noomp */
#define old_formylnodes(mesh) \
  for(int li_, cat_=0; cat_ < mesh->myln->nncats; cat_++) \
  NODELEVEL_Pragma(omp parallel for) \
  for(li_=0; li_ < mesh->myln->ncat[cat_]; li_++)

/* get leaf node from mesh, cat_ and li_ */
#define old_MyLnode mesh->myln->ln[cat_][li_]->node

/* get 1st leaf node on this proc from mesh, using cat_=0 and li_=0 */
#define old_MyLnode0 (mesh->myln->ln ? mesh->myln->ln[0][0]->node : 0)

/* get node from myid */
#define old_Lnode_myid(mesh, myid) \
  mesh->myln->ln[(myid) / mesh->myln->nm][(myid) % mesh->myln->nm]->node

/****************************************************************************/
/* Macros that should be used in most modules to access certain structs */
/****************************************************************************/

/* get parts of tElm struct */
#define Elm_p(elm)      (elm)->eploc->p
#define Elm_l(elm)      (elm)->eploc->l
#define Elm_eid(elm)    (elm)->eploc->eid
#define Elm_np(elm)     (elm)->np
#define Elm_mesh(elm)   (elm)->pat->mesh

/* tests for patch-face, boundary, OUTERBOUND, INNERBOUND */
#define Elm_patface(elm, face) elm_is_on_patface(elm, face)
#define Elm_on_BOUND(elm, face) \
  ( Elm_patface(elm, face) && \
    (elm)->pat->bfaces[face] && \
    (elm)->pat->bfaces[face]->boundary )
#define Elm_on_OUTERBOUND(elm, face) \
  ( Elm_patface(elm, face) && \
    (elm)->pat->bfaces[face] && \
    ((elm)->pat->bfaces[face]->boundary==OUTERBOUND) )
#define Elm_on_INNERBOUND(elm, face) \
  ( Elm_patface(elm, face) && \
    (elm)->pat->bfaces[face] && \
    ((elm)->pat->bfaces[face]->boundary==INNERBOUND) )


/* get parts of tNode struct */
#define Node_l(node)      Elm_l(node)
#define Node_eid(node)    Elm_eid(node)
#define Node_np(node)     Elm_np(elm)
#define Node_pat(node)    (node)->pat
#define Node_mesh(node)   Elm_mesh(elm)

/* tests for patch-face, boundary, OUTERBOUND, INNERBOUND */
#define Node_patface(node, face)       Elm_patface(node, face)
#define Node_on_BOUND(node, face)      Elm_on_BOUND(node, face)
#define Node_on_OUTERBOUND(node, face) Elm_on_OUTERBOUND(node, face)
#define Node_on_INNERBOUND(node, face) Elm_on_INNERBOUND(node, face)


/* marcos to start and stop node load timers */
#define LOAD_START loadtimer_start(node)
#define LOAD_STOP  loadtimer_stop(node)

/* loop over all points in a node */
#define forpoints(node,ijk)  for(ijk=0; ijk < node->np; ijk++)

/* get pointer to variable array */
#define VarA(node, varindex) \
  ((node->dat) ? node->dat->v[(varindex)] : 0)

/* loop over array */
#define forarray(array,k) \
  for(k=0; k<array->N; k++)

/* access Array */
#define Arrd(Arr) ((Arr) ? Arr->d : 0)
#define Arrn(Arr) ((Arr) ? Arr->n : 0)
#define ArrN(Arr) ((Arr) ? Arr->N : 0)
#define ArrNe(Arr) ((Arr) ? Arr->Ne : 0)
#define Arri(Arr) ((Arr) ? Arr->i : 0)

/* get double pointer to data in a variable */
#define Vard(node, varindex) \
  ((node->dat) ? node->dat->v[(varindex)]->d : 0)

/* get n of a variable */
#define Varn(node, varindex) \
  ((node->dat) ? node->dat->v[(varindex)]->n : 0)

/* get array pointer to surface data in a variable */
#define VarAaj(node, varindex, face)  ( (node->dat->s[(face)][(varindex)]) ? \
  node->dat->s[(face)][(varindex)]->ajsurf : 0 )

/* get double pointer to surface data in a variable */
#define Varaj(node, varindex, face) ( (VarAaj((node), (varindex), (face))) ? \
  node->dat->s[(face)][(varindex)]->ajsurf->d : 0 )

/* loop over one variable (it is in an array) */
//#define forvari(node,varindex, k) if(node->dat) forarray(node->dat->v[(varindex)], k)
#define forvari(node,varindex, k) \
  forarray(node->dat->v[(varindex)], k)

/* loop over var indices of a var list */
#define forvl(vl, vli) \
  for(vli=0; vli<vl->n; vli++)

/* get global var index from entry vli in VarList vl */
#define Vind(vl, vli) vl->index[vli]

/* get number of vars in VarList */
#define VLn(vl) vl->n

/* get integration weights for direc. dir on node */
#define Wquad(node, dir) node_Wq(node,dir)->d

/* grid points in Xb-coords in direc. dir on node */
#define Xbpts(node, dir) node_Xb(node,dir)->d

/***************************************************************************/
/* macros that should be used only in very particular advanced cases       */
/***************************************************************************/
/* get double pointer to data in a variable without checking if node has dat */
#define Vard_(node, varindex) node->dat->v[(varindex)]->d

/* get array pointer without check */
#define VarA_(node, varindex) node->dat->v[(varindex)]

/* access Array, without check */
#define Arrd_(Arr) (Arr->d)

/* get double pointer to surface data in a variable, without check */
#define Varaj_(node,varindex,face)  node->dat->s[(face)][(varindex)]->ajsurf->d

/* natural loop over my leaf nodes on this proc */
#define formylnodes_cat_i(mesh, cat, i) \
  for(cat=0; cat < mesh->myln->nncats; cat++) \
  for(i=0; i < mesh->myln->ncat[cat]; i++)

/* get node number i in cat.c out of nodelist on this proc */
#define GetMyNode_cat_i(mesh, c, i) mesh->myln->ln[c][i]->node

/* loop over all patches */
#define forpatches(mesh,patindex) \
  for(patindex=0; patindex < mesh->npats; patindex++)

/* loop over all bfaces in a patch */
#define forbfaces(pat, bface) \
  for(int face_=0; face_<6; face_++) \
  for(bface=pat->bfaces[face_]; bface; bface=bface->next)

/* loop over all bfaces in a patch on face f */
#define forbfacesonface(pat, f, bface) \
  for(bface=pat->bfaces[f]; bface; bface=bface->next)

/* loop over a node list nlist (type tNlist) */
#define fornodelist(nlist, elem) \
  for(elem=(nlist); elem; elem=elem->next)

/* loop over nodes in a node list */
#define fornodes(nlist, listnode) { \
  tNlist *elem_; \
  for(elem_ = nlist,            listnode = elem_ ?  elem_->node : 0; \
      elem_; elem_=elem_->next, listnode = elem_ ?  elem_->node : 0)
#define endfornodes }
/* loop over nodes in the node list in mesh */
#define forlnodes(meshORpat, listnode) { \
  tNlist *elem_; \
  for(elem_ = meshORpat->lns,   listnode = elem_ ?  elem_->node : 0; \
      elem_; elem_=elem_->next, listnode = elem_ ?  elem_->node : 0)
#define endforlnodes }


/****************************************************************************/
/* useful loops over point indices */
/****************************************************************************/

/* Indices */
#define Ind_n(i,j,k,n) ( (i) + (n[0])*( (j) + (n[1])*(k) ) )
#define i0_norm(i,j,k,norm) ( (i)*(norm==0) + (j)*(norm==1) + (k)*(norm==2) )
/* ijk = i + n0*j + n0*n1*k, thus:
   ijk/(n0*n1) = k
   (ijk - n0*n1*k)/n0 = j
   (ijk - n0*n1*k - n0*j ) = i   */
#define kOfInd_n(ijk,n)        ((ijk)/((n[0])*(n[1])))
#define jOfInd_n_k(ijk,n,k)    (((ijk) - (n[0])*(n[1])*(k))/(n[0]))
#define iOfInd_n_jk(ijk,n,j,k) ((ijk) - (n[0])*(n[1])*(k) - (n[0])*(j))
/* first and second index of i,j,k when we are in plane normal to norm */
#define i1_norm(i,j,k,norm) ( (i)*(norm==2 || norm==1) + (j)*(norm==0) )
#define i2_norm(i,j,k,norm) ( (k)*(norm==0 || norm==1) + (j)*(norm==2) )
/* direc. 1 and 2 that are normal to norm=0,1,2 */
#define Dir1_norm(norm) (norm==0)
#define Dir2_norm(norm) ( (2)*(norm==0 || norm==1) + (norm==2) )

#define ijk_inplaneN(N, i,j,k, i1,i2,pl) \
  switch((N)) { \
    case 0:  i=(pl); j=(i1); k=(i2); break; \
    case 1:  i=(i1); j=(pl); k=(i2); break; \
    case 2:  i=(i1); j=(i2); k=(pl); break; \
    default: errorexit("N has to be 0,1,2");  }

/* loop over i,j,k */
#define forijk(i,j,k, n) \
  for (k = 0; k < n[2]; k++) \
  for (j = 0; j < n[1]; j++) \
  for (i = 0; i < n[0]; i++)

/* loop over i,j,k but only include face points */
#define forfacepoints(i,j,k, n) \
  for (k = 0; k < n[2]; k++) \
  for (j = 0; j < n[1]; j++) \
  for (i = 0; i < n[0]; i++) \
  if( (i==0 || i==n[0]-1) || (j==0 || j==n[1]-1) || (k==0 || k==n[2]-1) )

/* loop over planes e.g. i=p plane */
#define forplane0(i,j,k, n, p) \
  for(i=(p), k = 0; k < (n[2]); k++) \
    for(     j = 0; j < (n[1]); j++)

#define forplane1(i,j,k, n, p) \
  for(j=(p), k = 0; k < (n[2]); k++) \
    for(     i = 0; i < (n[0]); i++)

#define forplane2(i,j,k, n, p) \
  for(k=(p), j = 0; j < (n[1]); j++) \
    for(     i = 0; i < (n[0]); i++)

/* same as forplane0/1/2, but we can specify the plane number N */
#define forplaneN(N, i,j,k, n, p) \
  for(k=(p)*((N)==2); ( k<(n[2]) ) && ( ( k==(p) ) || (N)!=2 ); k++) \
  for(j=(p)*((N)==1); ( j<(n[1]) ) && ( ( j==(p) ) || (N)!=1 ); j++) \
  for(i=(p)*((N)==0); ( i<(n[0]) ) && ( ( i==(p) ) || (N)!=0 ); i++)

/* like forplaneN but loop only over 4 corners */
#define forcornersN(N, i,j,k, n, p) \
  for(k=(p)*((N)==2); (k<n[2]) && (k==(p) || (N)!=2); k+=(n[2])-(n[2]>1)) \
  for(j=(p)*((N)==1); (j<n[1]) && (j==(p) || (N)!=1); j+=(n[1])-(n[1]>1)) \
  for(i=(p)*((N)==0); (i<n[0]) && (i==(p) || (N)!=0); i+=(n[0])-(n[0]>1))

/* same as forplaneN, but omit edges */
#define forinnerplaneN(N, i,j,k, n, p) \
  for(k=1+((p)-1)*((N)==2); ( k<(n[2])-((N)!=2) ) && ( ( k==(p) ) || (N)!=2 ); k++) \
  for(j=1+((p)-1)*((N)==1); ( j<(n[1])-((N)!=1) ) && ( ( j==(p) ) || (N)!=1 ); j++) \
  for(i=1+((p)-1)*((N)==0); ( i<(n[0])-((N)!=0) ) && ( ( i==(p) ) || (N)!=0 ); i++)

/* loop over planes smoothly without any jumping in i,j or k */
#define forplane0_nojump(i,j,k, n, p) \
  for(i=(p), k = 0; k < (n[2]); k++) \
    for(j =((n[1])-1)*(k%2); j < (n[1]) && j >= 0; j=j+1-2*(k%2))

#define forplane1_nojump(i,j,k, n, p) \
  for(j=(p), k = 0; k < (n[2]); k++) \
    for(i =((n[0])-1)*(k%2); i < (n[0]) && i >= 0; i=i+1-2*(k%2))

#define forplane2_nojump(i,j,k, n, p) \
  for(k=(p), j = 0; j < (n[1]); j++) \
    for(i =((n[0])-1)*(j%2); i < (n[0]) && i >= 0; i=i+1-2*(j%2))


/****************************************************************************/
/* defines useful for finite differences on Cartesian grids */
/****************************************************************************/
typedef struct tIDX {
  int i;  /* index of point */
  int f ; /* index of face where neighbor is (0-5), -1 means in my node */
} tIdx;

/* declaration of neighbor indices:
     ccc is the index of the center point 'i,j,k'
     mcc is the index of the i-1,j,k point
     pcc is the index of the i+1,j,k point
     Mcc is the index of the i-2,j,k point
     Pcc is the index of the i+2,j,k point
     m3cc is the index of the i-3,j,k point
     p3cc is the index of the i+3,j,k point */
#define tIdx_ptinds_faces   tIdx mcc[1], pcc[1], cmc[1], cpc[1], ccm[1], ccp[1]
#define tIdx_ptinds_faces2  tIdx Mcc[1], Pcc[1], cMc[1], cPc[1], ccM[1], ccP[1]
#define tIdx_ptinds_faces3 \
  tIdx m3cc[1], p3cc[1], cm3c[1], cp3c[1], ccm3[1], ccp3[1]
#define tIdx_ptinds_corners \
  tIdx mmm[1], mmp[1], mpm[1], mpp[1], pmm[1], pmp[1], ppm[1], ppp[1]
#define tIdx_ptinds_edges \
  tIdx mmc[1], mcm[1], cmm[1], ppc[1], pcp[1], cpp[1], \
       mpc[1], mcp[1], cmp[1], pmc[1], pcm[1], cpm

#define tIdx_ptinds_7   tIdx ccc[1];      tIdx_ptinds_faces
#define tIdx_ptinds2_13 tIdx_ptinds_7;    tIdx_ptinds_faces2
#define tIdx_ptinds3_19 tIdx_ptinds2_13;  tIdx_ptinds_faces3
#define tIdx_ptinds_19  tIdx_ptinds_7;    tIdx_ptinds_edges
#define tIdx_ptinds2_25 tIdx_ptinds_19;   tIdx_ptinds_faces2
#define tIdx_ptinds_27  tIdx_ptinds_19;   tIdx_ptinds_corners

/* 7 point stencil where some points are in the neighbors,
   here na[4] contains adjacent n from neighbor 4 */
#define set_ptinds_check_7(i,j,k, n, na)				\
  ccc->f = -1;								\
  ccc->i = Ind_n(i,j,k, n);						\
  if(i>0) { mcc->f =-1;  mcc->i = ccc->i-1; }				\
  else    { mcc->f = 0;  mcc->i = Ind_n(na[0][0]-1,j,k, na[0]); }	\
  if(j>0) { cmc->f =-1;  cmc->i = ccc->i-n[0]; }			\
  else    { cmc->f = 2;  cmc->i = Ind_n(i,na[2][1]-1,k, na[2]); }	\
  if(k>0) { ccm->f =-1;  ccm->i = ccc->i-n[0]*n[1]; }			\
  else    { ccm->f = 4;  ccm->i = Ind_n(i,j,na[4][2]-1, na[4]); }	\
  if(i+1<n[0]) { pcc->f =-1;  pcc->i = ccc->i+1; }			\
  else         { pcc->f = 1;  pcc->i = Ind_n(0,j,k, na[1]); }		\
  if(j+1<n[1]) { cpc->f =-1;  cpc->i = ccc->i+n[0]; }			\
  else         { cpc->f = 3;  cpc->i = Ind_n(i,0,k, na[3]); }		\
  if(k+1<n[2]) { ccp->f =-1;  ccp->i = ccc->i+n[0]*n[1]; }		\
  else         { ccp->f = 5;  ccp->i = Ind_n(i,j,0, na[5]); }

/* 7+6 point stencil where some points are in the neighbors */
#define set_ptinds2_check_6(i,j,k, n, na)				\
  if(i>1) { Mcc->f =-1;  Mcc->i = mcc->i-1; }				\
  else    { Mcc->f = 0;  Mcc->i = Ind_n(na[0][0]-2,j,k, na[0]); }	\
  if(j>1) { cMc->f =-1;  cMc->i = cmc->i-n[0]; }			\
  else    { cMc->f = 2;  cMc->i = Ind_n(i,na[2][1]-2,k, na[2]); }	\
  if(k>1) { ccM->f =-1;  ccM->i = ccm->i-n[0]*n[1]; }			\
  else    { ccM->f = 4;  ccM->i = Ind_n(i,j,na[4][2]-2, na[4]); }	\
  if(i+2<n[0]) { Pcc->f =-1;  Pcc->i = pcc->i+1; }			\
  else         { Pcc->f = 1;  Pcc->i = Ind_n(1,j,k, na[1]); }		\
  if(j+2<n[1]) { cPc->f =-1;  cPc->i = cpc->i+n[0]; }			\
  else         { cPc->f = 3;  cPc->i = Ind_n(i,1,k, na[3]); }		\
  if(k+2<n[2]) { ccP->f =-1;  ccP->i = ccp->i+n[0]*n[1]; }		\
  else         { ccP->f = 5;  ccP->i = Ind_n(i,j,1, na[5]); }

/* 7+6+6 point stencil where some points are in the neighbors */
#define set_ptinds3_check_6(i,j,k, n, na)				\
  if(i>2) { m3cc->f =-1;  m3cc->i = Mcc->i-1; }				\
  else    { m3cc->f = 0;  m3cc->i = Ind_n(na[0][0]-3,j,k, na[0]); }	\
  if(j>2) { cm3c->f =-1;  cm3c->i = cMc->i-n[0]; }			\
  else    { cm3c->f = 2;  cm3c->i = Ind_n(i,na[2][1]-3,k, na[2]); }	\
  if(k>2) { ccm3->f =-1;  ccm3->i = ccM->i-n[0]*n[1]; }			\
  else    { ccm3->f = 4;  ccm3->i = Ind_n(i,j,na[4][2]-3, na[4]); }	\
  if(i+3<n[0]) { p3cc->f =-1;  p3cc->i = Pcc->i+1; }			\
  else         { p3cc->f = 1;  p3cc->i = Ind_n(2,j,k, na[1]); }		\
  if(j+3<n[1]) { cp3c->f =-1;  cp3c->i = cPc->i+n[0]; }			\
  else         { cp3c->f = 3;  cp3c->i = Ind_n(i,2,k, na[3]); }		\
  if(k+3<n[2]) { ccp3->f =-1;  ccp3->i = ccP->i+n[0]*n[1]; }		\
  else         { ccp3->f = 5;  ccp3->i = Ind_n(i,j,2, na[5]); }
