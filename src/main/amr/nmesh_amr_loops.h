/* nmesh_amr_loops.h */
/* Wolfgang Tichy, 1/2019 */

/****************************************************************************/
/* Loops are performed by macros so that the user has to know very little
   about the implementation details. */
/****************************************************************************/

/****************************************************************************/
/* loops that should be used in most modules                                */
/****************************************************************************/
/* loop over my leaf nodes on this proc without OpenMP */
/* Note: for node in ln[c][i]: myid = myln->nm*c + i */
#define formylnodes_noomp(mesh, myid) \
  for(int li_, cat_=0; cat_ < mesh->myln->nncats; cat_++) \
  for(myid=mesh->myln->nm*cat_, li_=0; \
      li_ < mesh->myln->ncat[cat_]; myid++, li_++)

/* we could use OpenMP to parallelize the 2nd loop in formylnodes_noomp,
   but for now formylnodes and formylnodes_noomp do the same: */
#define formylnodes(mesh, myid) formylnodes_noomp(mesh, myid)

/* get node from myid */
#define MyNode(mesh, myid) \
  mesh->myln->ln[myid / mesh->myln->nm][myid % mesh->myln->nm]->node

/* loop over all points in a node */
#define forpoints(node,ijk)  for(ijk=0; ijk < node->np; ijk++)

/* get pointer to variable array */
#define VarA(node, varindex) \
  ((node->dat) ? node->dat->v[(varindex)] : 0)

/* loop over array */
#define forarray(array,k) \
  for(k=0; k<array->N; k++)

/* get double pointer to data in a variable */
#define Vard(node, varindex) \
  ((node->dat) ? node->dat->v[(varindex)]->d : 0)

/* get double pointer to surface data in a variable */
#define Varaj(node, varindex, face)  ( (node->dat->s[(face)][(varindex)]) ? \
  node->dat->s[(face)][(varindex)]->ajsurf->d : 0 )

/* loop of one variable (it is in an array*/
//#define forvari(node,varindex, k) if(node->dat) forarray(node->dat->v[(varindex)], k)
#define forvari(node,varindex, k) \
  forarray(node->dat->v[(varindex)], k)

/* get integration weights for direc. dir on node */
#define Wquad(node, dir) node->Wq[dir]->d

/****************************************************************************/
/* loops that should be used only in very particular advanced cases         */
/****************************************************************************/
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

/* loop over i,j,k */
#define forijk(i,j,k, n) \
  for (k = 0; k < n[2]; k++) \
  for (j = 0; j < n[1]; j++) \
  for (i = 0; i < n[0]; i++)

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
