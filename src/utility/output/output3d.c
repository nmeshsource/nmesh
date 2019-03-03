/* output3d.c */
/* Wolfgang Tichy, Feb. 2019 */

#include "nmesh.h"
#include "output.h"





/* 3d output */
void output3d_meshvar(tMesh *mesh, char *name, int It, double T)
{
  tNode *node;
  int vi = Ind(name);
  FILE *fp;
  int nseries;
  int vtk      = Getv(Par("3dformat"), "vtk");
  char *outdir = Gets(Par("outdir"));
  tOutpars par[1];

  /* pars we may need for vtk or others */
  par->name          = name;
  par->text          = Getv(Par("3dformat"), "text");
  par->binary        = Getv(Par("3dformat"), "binary");
  par->arrange_as_1d = Getv(Par("3dformat"), "arrange_as_1d");
  par->flt           = Getv(Par("3dformat"), "float");
  par->dbl           = Getv(Par("3dformat"), "double");

  /* a number that counts the output */
  nseries = TimeForMeshOutput_di_dt(mesh,Geti(Par("3doutiter")),
                                    Getd(Par("3douttime")));
  /* loop over all nodes */
  forlnodes(mesh, node)
  {
    if(node->dat)
    if(node->dat->v[vi])
    {
      int p = node->pat->p;
      char ns[100];
      
      /* find string that identifies node */
      node_location_str(node, ns,100);

      /* set some more pars */
      par->nodeloc = ns;
      par->p       = p;

      /* write files */
      if(vtk || 1) /* can do only VTK right now */
      {
        /* VTK output: one file per time step in separate subdirectories */
        fp = fopen_vtk(name, outdir, "XYZ", p, ns, nseries-1);
        write3d_vtk(node, fp, VarA(node, vi), It,T, nseries-1, par);
        fclose(fp);
      }
    }
    /* sysnchronize, so that we write only one node at a time */
    nMPI_barrier();
  } endforlnodes;
}
