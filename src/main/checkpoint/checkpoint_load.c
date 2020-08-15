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
  char *buffer=NULL;
  long nbuffer;
  char buf[1000];
  long off, len;
  tPat *pat = NULL;
  char par[1000], val[1000];
  int clabel=0, nmax=0, useF=0;
  int n[] = { 0,0,0 };
  double bbox[] = { 0.,0.,0.,0.,0.,0. };

/*
long len, off=0;
char ttt[] = "  domu=read 99 newl=wwww 6    6  = qqqq3=";
char *str=cmalloc(999);
printf("ttt = |%s|:%d\n", ttt, sizeof(ttt)-1);
off=0;
while((off = str_from_buf(ttt, sizeof(ttt)-1, off, '=', str, 999, &len))>=0L)
{
printf("|%s|: %ld\n", str, off);
}
printf("off=%ld: %d\n", off, ttt[off-1]);
exit(8);
*/

  /* get file into mem buffer on proc0 */
  if(Rank0)
  {
    FILE *fp;

    /* open file */
    fp = fopen(fname, "rb");
    if(!fp) errorexits("failed opening %s", fname);

    /* find number of bytes in file and write them into buffer */
    nbuffer = nbytes_infile(fp);
    buffer = cmalloc(nbuffer);
    fread(buffer, sizeof(char), nbuffer, fp);

    fclose(fp);
  }

  /* broadcast buffer to all MPI ranks */
  nMPI_Bcast(&nbuffer,1, nMPI_LONG, 0);
  if(!Rank0) buffer = cmalloc(nbuffer);
  nMPI_Bcast(buffer,nbuffer, nMPI_CHAR, 0);

  /* read buffer line by line */
  off = 0;
  while((off = str_from_buf(buffer,nbuffer, off, '\n', buf,999, &len))>=0)
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

    /* read various info pieces */
    if(get_par_from_str(buf, par, "=", val, 999) && val[0])
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

  /* do not load balance root nodes here! datrank=-1 saves memory only
     until simple_load_balance(mesh) is called for the 1st time! */

  free(buffer);

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
  FILE *fp=NULL;
  char *buffer=NULL;
  long nbuffer;
  long first = 524288;
  long extra = 1024;
  long nbuffer_allocd = first + extra;
  int file_end;

  /* alloc buffer */
  buffer = cmalloc(nbuffer_allocd);

  /* open file on rank0 */
  if(Rank0)
  {
    fp = fopen(fname, "rb");
    if(!fp) errorexits("failed opening %s", fname);
  }

  /* get part of file into mem buffer on proc0 */
  file_end = 0;
  do
  {
    char buf[1000];
    char str[1000];
    tNode *parent;
    int n[3];
    int p = -1; /* patch number read from file */
    int p_prev; /* previous patch number read from file */
    int lp = -1; /* ref. level of parent */
    int lp_prev; /* ref. level of previous parent */
    tNlist *elem = mesh->lns; /* first element in leaf node list */

    /* read a chunk from file into buffer */
    if(Rank0)
    {
      char *tailbuf;
      long pos;
      int lastone, newl, c;

      /* read a certain number of parent nodes and their child0->n */
      nbuffer = fread(buffer, sizeof(char), first, fp);
      if(nbuffer<first) file_end = 1;

      tailbuf = buffer + nbuffer;
      pos = 0;
      lastone = newl = 0;
      while((c=fgetc(fp)) != EOF)
      {
        tailbuf[pos++] = c;
        if(c=='_') lastone = 1; /* from now on we read 4 more \n */
        if(lastone && c=='\n') newl++;
        if(newl>=4) break;
      }
      nbuffer += pos;
      if(c == EOF) file_end = 1;
    }

    /* broadcast buffer to all MPI ranks */
    nMPI_Bcast(&file_end,1, nMPI_INT, 0);
    nMPI_Bcast(&nbuffer,1, nMPI_LONG, 0);
    if(!Rank0) buffer = cmalloc(nbuffer);
    nMPI_Bcast(buffer,nbuffer, nMPI_CHAR, 0);

//PRF;printf(": nbuffer=%ld\n", nbuffer);

    {
      long off, len;

      /* read buffer line by line */
      off = 0;
      while((off = str_from_buf(buffer,nbuffer, off, '\n', buf,999, &len))>=0)
      {
        /* all node names contain an '_' */
        if(strstr(buf, "_"))
        {
          int d;
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
          for(d=0; d<3; d++)
          {
            off = str_from_buf(buffer,nbuffer, off, '\n', str,999, &len);
            n[d] = atoi(str);
          }
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
    }

  } while(!file_end);

  /* close file */
  if(Rank0) fclose(fp);

  PRF;printf(": mesh->iteration=%d mesh->time=%g\n",
             mesh->iteration, mesh->time);
  fflush(stdout);

  /* make sure all nodes have new current nids */
  update_mesh_myln_node_nid(mesh);

  /* load balance all leaf nodes */
  simple_load_balance(mesh);
  //printmesh(mesh);
  PRF;printf(": total number of leaf nodes mesh->nln=%ld\n", mesh->nln);
  fflush(stdout);

  free(buffer);
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
