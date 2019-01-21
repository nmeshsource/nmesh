/* storage.c */
/* Wolfgang Tichy, 1/2019 */

#include "nmesh.h"
#include "amr.h"




/**************************************************************************/
/* basic memory management */

/* allocate mesh */
tMesh *alloc_mesh(void)
{
  tMesh *mesh;

  mesh = calloc(1, sizeof(tMesh));
  if (!mesh) errorexit("out of memory");

  return mesh;
}


/* allocate patch */
tPat *alloc_patch(tMesh *mesh, int p, int nD) 
{
  tPat *pat;

  pat = calloc(1, sizeof(tPat));
  if (!pat) errorexit("out of memory");

  pat->mesh = mesh;
  pat->p = p;
  pat->nD = nD;
    
  /* allocate storage for data pointers, they default to NULL */

  /* get mem. for diff. matrices */
  pat->D = calloc(nD, sizeof(pat->D[0]));
  if( !(pat->D) )
    errorexit("out of memory for diff. matrices");

  return pat;
} 


/* free pat, currently leaves mesh untouched */
void free_pat(tPat *pat) 
{
  int i;

  if (!pat) return;

  //for (i = 0; i < pat->mesh->nvariables; i++)
  //  disablevarcomp_inpat(pat, i);
  //free(pat->v);

  //free_all_bfaces(pat);

  free(pat);
}


/* free mesh */
void free_mesh(tMesh *mesh)
{
  int i;

  if (!mesh) return;
  for (i = 0; i < mesh->npatches; i++)
    free_pat(mesh->pat[i]);
  free_mesh_only(mesh);
}




 
// /* add a bface to a pat, f denotes the pat face (0 to 5),
//    return index in bface list */
// int add_empty_bface(tPat *pat, int f)
// {
//   int fi = pat->nbfaces; /* add bface in this pos. in bface list */
//   void *ret;
// 
//   /* increase size of bface list */
//   ret = realloc( pat->bface, (sizeof( *(pat->bface) ))*(fi+1) );
//   if(ret==NULL)  errorexit("add_bface: not enough memory for pat->bface");
//   pat->bface = ret;
// 
//   /* mem for new bface */
//   ret = calloc( 1, sizeof( *(pat->bface[fi]) ) );
//   if(ret==NULL) errorexit("add_bface: not enough memory for pat->bface[n]");
//   /* add new bface */
//   pat->bface[fi] = ret;
//   pat->nbfaces = fi+1;
// 
//   /* set some bface info */
//   pat->bface[fi]->mesh = pat->mesh;
//   pat->bface[fi]->b    = pat->b;
//   pat->bface[fi]->f    = f;
//   pat->bface[fi]->fi   = fi;
//   pat->bface[fi]->ob   = -1; /* other pat is not known yet */
//   pat->bface[fi]->ofi  = -1; /* fi in other pat is not known yet */
//   pat->bface[fi]->oXi  = -1; /* var indices of other coords not known yet */
//   pat->bface[fi]->oYi  = -1; /* var indices of other coords not known yet */
//   pat->bface[fi]->oZi  = -1; /* var indices of other coords not known yet */
//   return fi;
// }
// 
// /* add a point ijk on face f to a bface with index fi, returns fi */
// /* if called with fi<0, it first calls add_empty_bface and returns the new fi */
// int add_point_to_bface_inpat(tPat *pat, int fi, int ijk, int f)
// {
//   tBface *bface;
//   /* make new bface when needed */
//   if(fi<0 || fi>=pat->nbfaces) fi = add_empty_bface(pat, f);
//   bface = pat->bface[fi];
//   /* make PointList id needed */
//   if(bface->fpts==NULL)
//     bface->fpts = AllocatePointList(pat->mesh);
//   /* add point ijk */
//   AddToPointList(bface->fpts, pat->b, ijk);
//   /* check if we have more than one face on this bface */
//   if(bface->f!=f) bface->f = -1;
//   return fi;
// }
// 
// /* duplicate bface without pointlist fpts */
// tBface *duplicate_bface_without_fpts_for_mesh(tBface *bface0, tMesh *mesh)
// {
//   tBface *bface;
//   void *ret;
//   if(bface0==NULL) return NULL;
// 
//   /* mem for new bface */
//   ret = calloc( 1, sizeof( *(bface) ) );
//   if(ret==NULL)
//     errorexit("duplicate_bface_without_fpts_for_mesh: not enough memory");
//   bface = ret;
// 
//   /* make a shallow copy of the struct */
//   *bface = *bface0;
//   /* now set mesh pointer */
//   bface->mesh = mesh;
//   /* remove pointer to bface0->fpts */
//   bface->fpts = NULL;
//   
//   return bface;
// }
// 
// /* duplicate bface with pointlist fpts */
// tBface *duplicate_bface_for_mesh(tBface *bface0, tMesh *mesh)
// {
//   tBface *bface = duplicate_bface_without_fpts_for_mesh(bface0, mesh);
//   /* copy fpts */
//   bface->fpts = DuplicatePointList_for_mesh(bface0->fpts, mesh);
//   return bface;
// }
// 
// /* free a bface */
// void free_bface(tBface *bface)
// {
//   if(bface!=NULL)
//   {
//     /* free the point lists of the faces */ 
//     FreePointList(bface->fpts);
//     free(bface);
//   }
// }
// 
// /* free all bfaces in pat */
// void free_all_bfaces(tPat *pat)
// {
//   int i;
//   for(i=0; i<pat->nbfaces; i++)  free_bface(pat->bface[i]);
//   free(pat->bface);
//   pat->bface = NULL;
//   pat->nbfaces = 0;
// }
// 
// /* remove a bface with index fi from a pat, return number of bfaces removed */
// int remove_bface(tPat *pat, int fi)
// {
//   void *ret;
//   int nbfaces = pat->nbfaces;
//   int i;
// 
//   /* return 0 if bface does not exist */
//   if(fi<0 || fi>=nbfaces) return 0;
// 
//   /* free the bface */
//   free_bface(pat->bface[fi]);
// 
//   /* shift bfaces behind fi one position to the front */
//   for(i=fi; i<nbfaces-1; i++)
//   {
//     pat->bface[i] = pat->bface[i+1];
//     pat->bface[i]->fi = i;
//   }
// 
//   /* reduce size of bface list */
//   nbfaces--;
//   ret = realloc( pat->bface, (sizeof( *(pat->bface) ))*(nbfaces) );
//   if(ret==NULL && nbfaces!=0)
//     errorexit("remove_bface: not enough memory for pat->bface");
//   pat->bface = ret;
//   pat->nbfaces = nbfaces;
// 
//   return 1;
// }
// 
// /* look for empty bfaces and remove them */
// int remove_bfaces_with_NULL_fpts(tPat *pat)
// {
//   int n, fi;
//   n=0;
//   for(fi=0; fi<pat->nbfaces; fi++)
//     if(pat->bface[fi]->fpts == NULL)
//     {
//       int r = remove_bface(pat, fi); /* this decreases pat->nbfaces */
//       fi=fi-r;  /* go back by one in fi so that for-loop covers fi again */
//       n++;      /* count number of bfaces removed */
//     }
//   return n;
// }


/**********************************************************************/
/* storage for arrays in the nodes */
/**********************************************************************/
tDat *alloc_dat(int listlen)
{
  tDat *dat;

  dat = calloc(1, sizeof(tDat));
  if (!dat) errorexit("out of memory for dat");
  dat->v = calloc(listlen, sizeof(*tArray));
  if (!dat) errorexit("out of memory for dat->v");
  dat->g = calloc(listlen, sizeof(*tArray));
  if (!dat) errorexit("out of memory for dat->g");

  return dat;
}


/**********************************************************************/
/* storage for variables in the nodes */
/**********************************************************************/

/* enable one component of variable on one node */
void enablevarcomp_innode(tNode *node, int i)
{
  if (!node->v[i]) {
    node->v[i] = calloc( node->nnodes, sizeof(double) );
    if (!node->v[i]) {
      printf("enable storage for variable %d = %s: out of memory\n",
	     i, VarName(i));
      errorexit("");
    }
    nenabled++;
    if (storage_verbose) {
      printf("enabling variable %d = %s in b=%d, nenabled=%d => ",
	     i, VarName(i), node->b, nenabled);
      printf("%d*%d = %d bytes\n", 
	     node->nnodes, (int) sizeof(double), 
	     node->nnodes * ((int) (sizeof(double))) );
    } 
  }
}


/* disable one component of variable */
void disablevarcomp_innode(tNode *node, int i) 
{
  if (node->v[i]) {
    free(node->v[i]);
    node->v[i] = 0;
    nenabled--;
    if (storage_verbose) {
      printf("disabling variable %d = %s in b=%d, nenabled=%d => ",
             i, VarName(i), node->b, nenabled);
      printf("%d*%d = %d bytes\n", 
	     node->nnodes, (int) sizeof(double), 
	     node->nnodes * ((int) (sizeof(double))) );
    } 
  }
}


/* enable all components of a variable on one node */
void enablevar_innode(tNode *node, int i)
{
  int j, n = VarNComponents(i);

  for (j = 0; j < n; j++)
    enablevarcomp_innode(node, i+j);
}


/* disable all components of a variable on one node */
void disablevar_innode(tNode *node, int i)
{
  int j, n = VarNComponents(i);

  for (j = 0; j < n; j++)
    disablevarcomp_innode(node, i+j);
}






/* enable variable list */
void enablevarlist(tVarList *v) 
{
  int i;

  if (v) 
    for (i = 0; i < v->n; i++)
      enablevarcomp(v->grid, v->index[i]);
}




/* disable variable list */
void disablevarlist(tVarList *v) 
{
  int i;

  if (v)
    for (i = 0; i < v->n; i++)
      disablevarcomp(v->grid, v->index[i]);
}




/* enable all variables found on given grid */
void enablesamevars(tGrid *grid, tGrid *newgrid)
{
  int i, b;

  if (grid->nvariables != newgrid->nvariables)
  {
    printf("nvariables: old %d, new %d\n", 
	   grid->nvariables, newgrid->nvariables);
    errorexit("enablesamevars: need same number of variables");
  }
  if (grid->nboxes != newgrid->nboxes)
  {
    printf("nboxes: old %d, new %d\n", 
	   grid->nboxes, newgrid->nboxes);
    errorexit("enablesamevars: need same number of boxes");
  }

  for (i = 0; i < newgrid->nvariables; i++)
    for (b = 0; b < newgrid->nboxes; b++)
    {
      tNode *newbox = newgrid->box[b];
      if(grid->box[b]->v[i])
        enablevarcomp_inbox(newbox, i);
    }
}



/* create room for more variables in one box */
void realloc_boxvariables(tNode *box, int nvariables)
{
  tGrid *grid = box->grid;
  int i;
  int box_nvariables = grid->nvariables;

  if (PR) printf(" realloc_boxvariables in box%d from %d to %d\n", 
		 box->b, grid->nvariables, nvariables);

  if(box->v == NULL) box_nvariables=0;
  box->v = (double **) realloc(box->v, sizeof(double *) * nvariables);
  if(nvariables > box_nvariables) 
    for(i=box_nvariables; i<nvariables; i++)  box->v[i] = 0;
}

/* create room for more variables */
void realloc_gridvariables(tGrid *grid, int nvariables)
{
  int b;

  if (PR) printf("realloc_gridvariables from %d to %d\n", 
		 grid->nvariables, nvariables);

  for (b = 0; b < grid->nboxes; b++)
  {
    tNode *box = grid->box[b];
    realloc_boxvariables(box, nvariables);
  }
  grid->nvariables = nvariables;
}

