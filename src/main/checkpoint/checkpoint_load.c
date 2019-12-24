/* checkpoint_load.c */
/* Wolfgang Tichy, 8/2019 */


#include "nmesh.h"
#include "checkpoint.h"

/******************************************************************/
/* some functions to load nmesh data for checkpoints  */
/******************************************************************/

/******************************************************************/
/* functions to load patches */
/******************************************************************/
/* load patch info */
int checkpoint_load_patches(tMesh *mesh, char *fname)
{
  int rk;

  /* MPI motivated loop to assign work */
  for(rk=0; rk<nMPI_size(); rk++)
  {
    /* do work when it is my turn */
    if(rk == nMPI_rank())
    {
      FILE *fp;
      char buf[1000];
      tPat *pat = NULL;
      char par[1000], val[1000];
      int clabel=0, nmax=0, useF=0;
      int n[] = { 0,0,0 };
      double bbox[] = { 0.,0.,0.,0.,0.,0. };

      /* open source file */
      fp = fopen(fname, "rb");
      if(!fp) errorexits("failed opening %s", fname);

      /* read file line by line */
      while(fgets(buf,999, fp))
      {
        /* CI-> signifies that we know all we need to create a new patch */
        if(strcmp(buf, " CI->\n")==0)
        {
          //prbbox(bbox, 3);
          //printf("n = %d %d %d\n", n[0],n[1],n[2]);
          /* We set datrank=-1 to save memory. No dat is allocated
             anywhere! */
          pat = add_patch(mesh, bbox, n, nmax, -1);
          useF = 0;
        }
/*
char ttt[] = "  domu=read 99 \n";
printf("ttt = |%s|\n", ttt);
get_par_from_str(ttt, par,val,999);
printf("|%s| = |%s|\n", par,val);
exit(8);
*/
        /* read various info pieces */
        get_par_from_str(buf, par, val, 999);
        if(val[0])
        {
          //printf("%s = %s\n", par, val);
          if(strcmp(par, "time")==0)      mesh->time = atof(val);
          if(strcmp(par, "iteration")==0) mesh->iteration = atoi(val);
          if(strcmp(par, "dt")==0)        mesh->dt = atof(val);

          if(strcmp(par, "bbox[0]")==0) bbox[0] = atof(val);
          if(strcmp(par, "bbox[1]")==0) bbox[1] = atof(val);
          if(strcmp(par, "bbox[2]")==0) bbox[2] = atof(val);
          if(strcmp(par, "bbox[3]")==0) bbox[3] = atof(val);
          if(strcmp(par, "bbox[4]")==0) bbox[4] = atof(val);
          if(strcmp(par, "bbox[5]")==0) bbox[5] = atof(val);
          if(strcmp(par, "nmax")==0) nmax = atoi(val);
          if(strcmp(par, "rnode->n[0]")==0) n[0] = atoi(val);
          if(strcmp(par, "rnode->n[1]")==0) n[1] = atoi(val);
          if(strcmp(par, "rnode->n[2]")==0) n[2] = atoi(val);
          if(strcmp(par, "s[0]")==0)  pat->CI->s[0] = atof(val);
          if(strcmp(par, "s[1]")==0)  pat->CI->s[1] = atof(val);
          if(strcmp(par, "s[2]")==0)  pat->CI->s[2] = atof(val);
          if(strcmp(par, "s[3]")==0)  pat->CI->s[3] = atof(val);
          if(strcmp(par, "s[4]")==0)  pat->CI->s[4] = atof(val);
          if(strcmp(par, "s[5]")==0)  pat->CI->s[5] = atof(val);
          if(strcmp(par, "xc[0]")==0) pat->CI->xc[0] = atof(val);
          if(strcmp(par, "xc[1]")==0) pat->CI->xc[1] = atof(val);
          if(strcmp(par, "xc[2]")==0) pat->CI->xc[2] = atof(val);
          if(strcmp(par, "dom")==0)   pat->CI->dom = atoi(val);
          if(strcmp(par, "type")==0)  pat->CI->type = atoi(val);
          if(strcmp(par, "use_FSurf")==0)  useF = atoi(val);
          if(strcmp(par, "coordinates_get_label(pat)")==0)
          {
            /* this signifies end of patch info, so now set final pat info */
            clabel = atoi(val);
            switch(clabel)
            {
            case CubedSphere:
              set_1_CubedSphere_pat(pat, 0, useF);
              break;
            case stretchedCubedSphere:
              set_1_CubedSphere_pat(pat, 1, useF);
              break;
            }
          }
        } /* end of if get_par_from_str */
      }

      fclose(fp);
    }
    /* wait until everyone is here */
    nMPI_barrier();
  } /* end rk-loop */

  /* do not load balance root nodes here! datrank=-1 saves memory only
     until simple_load_balance(mesh) is called for the 1st time! */

  /* setup all bfaces and root node connections */
  amr_set_bfaces_and_rnode_nfaces_fnb(mesh, 0);
  printallbfaces(mesh);

  return 0;
}


/******************************************************************/
/* functions to load nodes */
/******************************************************************/
/* load node info */
int checkpoint_load_nodes(tMesh *mesh, char *fname)
{
  int rk;

  /* MPI motivated loop to assign work */
  for(rk=0; rk<nMPI_size(); rk++)
  {
    /* do work when it is my turn */
    if(rk == nMPI_rank())
    {
      FILE *fp;
      char buf[1000];
      tNode *parent;
      int n[3];
      int p = -1; /* patch number read from file */
      int p_prev; /* previous patch number read from file */
      int lp = -1; /* ref. level of parent */
      int lp_prev; /* ref. level of previous parent */
      tNlist *elem = mesh->lns; /* first element in leaf node list */

      /* open source file */
      fp = fopen(fname, "rb");
      if(!fp) errorexits("failed opening %s", fname);

      /* read file line by line */
      while(fgets(buf,999, fp))
      {
        /* all node names contain an '_' */
        if(strstr(buf, "_"))
        {
          tNlist *children, *lastchild;

          /* strip trailing '\n' from buf */
          buf[strlen(buf)-1] = 0;

          /* find parent node from its name in buf */
          parent = node_from_nodename(mesh, buf);
          p_prev = p;
          p = parent->pat->p;
          lp_prev = lp;
          lp = parent->l;
          //printf("buf=%s\n", buf);
          //printnode(parent);

          /* read n for child0 */
          fscanf(fp, "%d", &(n[0]));
          fscanf(fp, "%d", &(n[1]));
          fscanf(fp, "%d", &(n[2]));
          //printf("n = %d %d %d\n", n[0],n[1],n[2]);
          //fflush(stdout);

          /* update starting point in leaf node list of the 'for'-loop below */
          if(p<p_prev || !elem->next || lp!=lp_prev)
          {
            elem = mesh->lns;
            /* We could load balance here, BUT ONLY if all MPI-proc. go
               through the file in parallel! In any case, it's better to
               do load balance only once after all nodes are read. */
          }
          /* find element in nodelist with parent: */
          /* we could start search loop at elem=mesh->lns,
             but this would be slower */
          for( ; (elem) && (elem->node != parent); elem = elem->next) ;
          if(0 && !elem) /* try again from beginning if parent is not found */
          {
            for(elem = mesh->lns; (elem) && (elem->node != parent);
                elem = elem->next) ;
          }
          if(!elem)
          {
            char ns[100];
            nodename(parent, ns,99);
            errorexits("cannot find parent %s in mesh->lns", ns);
          }

          /* make 8 child nodes */
          children = make8_child_nodes(parent, n);

          /* update mesh->lns if needed and add children to leaf node list */
          if(elem == mesh->lns) mesh->lns = first_nodelist(children);
          lastchild = replace1_in_nodelist(elem, children, 1);

          /* set elem to last child we added */
          elem = lastchild;
        }
      }
      fclose(fp);
      PRF;printf(": mesh->iteration=%d mesh->time=%g\n",
                 mesh->iteration, mesh->time);
      fflush(stdout);
    }
    /* wait until everyone is here */
    nMPI_barrier();
  } /* end rk-loop */

  /* make sure all nodes have new current nids */
  update_mesh_myln_node_nid(mesh);

  /* load balance all leaf nodes */
  simple_load_balance(mesh);
  //printmesh(mesh);
  PRF;printf(": total number of leaf nodes mesh->nln=%ld\n", mesh->nln);
  fflush(stdout);

  return 0;
}


/******************************************************************/
/* functions to load variables */
/******************************************************************/
/* load all EvoVars */
int checkpoint_load_Vars(tMesh *mesh, char *fname)
{
  int rk;

  /* MPI motivated loop to assign work */
  for(rk=0; rk<nMPI_size(); rk++)
  {
    /* do work when it is my turn */
    if(rk == nMPI_rank())
    {
      tVarList *vl;
      FILE *fp;

      /* open file */
      fp = fopen(fname, "rb");
      if(fp)
      {
        /* get varlist from file */
        vl = checkpoint_make_vl(fp, mesh);

        /* now read data for vars in little endian format */
        checkpoint_read_vl(fp, vl, 0);
        PRF;printf(": finished reading varlist.\n");
        prvarlist(vl);
        fflush(stdout);

        vlfree(vl);
        fclose(fp);
      }
      else
      {
        errorexits("failed opening %s", fname);
        //PRF;printf(": failed opening %s\n", fname);
        fflush(stdout);
      }
    }
    /* wait until everyone is here */
    nMPI_barrier();
  } /* end rk-loop */

  return 0;
}

/* make varlist that we find in checkpoint file */
tVarList *checkpoint_make_vl(FILE *fp, tMesh *mesh)
{
  tVarList *vl = vlalloc(mesh);
  char buf[1000];
  int i, nvars;

  /* read first comment line */
  fgets(buf,999, fp);

  /* read number of vars */
  fscanf(fp, "%d", &nvars);

  /* read names and make varlist */
  for(i=0; i<nvars; i++)
  {
    fscanf(fp, "%s", buf);
    vlpushone(vl, Ind(buf));
  }

  return vl;
}

/* output varlist on each node */
void checkpoint_read_vl(FILE *fp, tVarList *vl, int read_big)
{
  tMesh *mesh = vl->mesh;
  char buf[1000];

  /* find string node */
  while(fgets(buf,999, fp))
  {
    tNode *node;
    char name[256];
    int np, vli, vi, found_node;

    if(strcmp(buf, "{\n")==0)
    {
      /* read node info */
      fscanf(fp, "%s", name);
      fscanf(fp, "%d", &np);
      fgets(buf,999, fp); /* use fgets to also read the '\n' after np */

      /* find node from its name */
      node = node_from_nodename(mesh, name);
      found_node = 1;
    }
    else
    {
      found_node = 0;
    }
    while(found_node)
    {
      /* check for end / read var info */
      fgets(buf,999, fp); /* use fgets to read '}' or vli plus '\n' */
      if(strcmp(buf, "}\n")==0) break;
      vli = atoi(buf);
      vi  = Vind(vl, vli);
      //printf("name=%s np=%d vli=%d\n", name, np, vli);
      //printf("%s nid=%ld vi=%d\n", nodename(node,buf,99), node->nid, vi);

      /* check if this var needs to be read on this node */
      if(node->dat)
      {
        double *v;

        /* switch on var on this node and get pointer to its data */
        enablevarcomp_innode(node, vi);
        v = Vard(node, vi);

        /* read var as raw binary */
        if(read_big) fread_big(v, sizeof(double), np, fp);
        else         fread_little(v, sizeof(double), np, fp);
        //if(node->nid==28) printf("v[]=%g\n", v[0]);
      }
      else /* just over-read data */
      {
        long offset = sizeof(double) * np;
        fseek(fp, offset, SEEK_CUR);
      }
      fgets(buf,999, fp); /* use fgets to also read the '\n' after var. */
    } /* end while(found_node) */
  }
}
