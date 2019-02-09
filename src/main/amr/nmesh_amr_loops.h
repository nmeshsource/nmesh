/* nmesh_amr_loops.h */
/* Wolfgang Tichy, 1/2019 */

/* Loops are performed by macros so that the user has to know very little
   about the implementation details. */

/****************************************************************************/
/* loop over all patches */
#define forpatches(mesh,patindex) \
  for(patindex=0; patindex < mesh->npats; patindex++)

///* loop over all points in a node */
//#define forpoints(node,ijk) \
//  for(ijk = 0; ijk < node->np; ijk++)

#define forarray(array,k) \
  for(k=0; k<array->N; k++)

#define GetVarDpointer(node, varindex) \
  ((node->dat) ? node->dat->v[(varindex)]->a : 0)

//#define forvari(node,varindex, k) \
//  if(node->dat) forarray(node->dat->v[(varindex)], k)
#define forvari(node,varindex, k) \
  forarray(node->dat->v[(varindex)], k)

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
  
/* to loop over leaf node array on this proc */
#define formylnodes(mesh, lni) \
  for(lni=0; lni < mesh->nmyln; lni++)


/* loop over planes e.g. i=p plane */
#define forplane1(i,j,k, n, p) \
  for(i=(p), k = 0; k < (n[2]); k++) \
    for(     j = 0; j < (n[1]); j++)

#define forplane2(i,j,k, n, p) \
  for(j=(p), k = 0; k < (n[2]); k++) \
    for(     i = 0; i < (n[0]); i++)

#define forplane3(i,j,k, n, p) \
  for(k=(p), j = 0; j < (n[1]); j++) \
    for(     i = 0; i < (n[0]); i++)

/* same as forplane1/2/3, but we can specify the plane number N */
#define forplaneN(N, i,j,k, n, p) \
  for(k=(p)*((N)==3); ( k<(n[2]) ) && ( ( k==(p) ) || (N)!=3 ); k++) \
  for(j=(p)*((N)==2); ( j<(n[1]) ) && ( ( j==(p) ) || (N)!=2 ); j++) \
  for(i=(p)*((N)==1); ( i<(n[0]) ) && ( ( i==(p) ) || (N)!=1 ); i++)

/* same as forplaneN, but omit edges */
#define forinnerplaneN(N, i,j,k, n, p) \
  for(k=1+((p)-1)*((N)==3); ( k<(n[2])-((N)!=3) ) && ( ( k==(p) ) || (N)!=3 ); k++) \
  for(j=1+((p)-1)*((N)==2); ( j<(n[1])-((N)!=2) ) && ( ( j==(p) ) || (N)!=2 ); j++) \
  for(i=1+((p)-1)*((N)==1); ( i<(n[0])-((N)!=1) ) && ( ( i==(p) ) || (N)!=1 ); i++)

/* loop over planes smoothly without any jumping in i,j or k */
#define forplane1_nojump(i,j,k, n, p) \
  for(i=(p), k = 0; k < (n[2]); k++) \
    for(j =((n[1])-1)*(k%2); j < (n[1]) && j >= 0; j=j+1-2*(k%2))

#define forplane2_nojump(i,j,k, n, p) \
  for(j=(p), k = 0; k < (n[2]); k++) \
    for(i =((n[0])-1)*(k%2); i < (n[0]) && i >= 0; i=i+1-2*(k%2))

#define forplane3_nojump(i,j,k, n, p) \
  for(k=(p), j = 0; j < (n[1]); j++) \
    for(i =((n[0])-1)*(j%2); i < (n[0]) && i >= 0; i=i+1-2*(j%2))

/* loop over bfaces in a patch */
//#define forbfaces(pat,fi)  for(fi = 0; fi<pat->nbfaces; fi++)
