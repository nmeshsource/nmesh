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
  int useF=0;
  int n[] = { 0,0,0 };
  int pt_typ[] = { 0,0,0 };
  double bbox[] = { 0.,0.,0.,0.,0.,0. };

/*
char ttt[] = "  domu=read 99 newl=wwww 6    6  = qqqq3=xyz";
int sl = 20;
char *str=cmalloc(sl);
printf("ttt = |%s|:%ld\n", ttt, sizeof(ttt)-1);
off=0;
while((off = str_from_buf(ttt,sizeof(ttt)-1, off, '=', str,sl, &len))>=0L)
{
printf("|%s|: %ld %ld\n", str, off,len);
}
printf("off=%ld\n", off);
exit(8);
*/

  /* get file into mem buffer on proc0 */
  if(Rank0)
  {
    FILE *fp;
    int IObufsz = Geti(Par("fread_bufsize"));
    char *IObuf; /* larger buffer for file read */

    /* open file */
    fp = fopen_buf(fname, "rb", &IObuf,IObufsz);
    if(!fp) errorexits("failed opening %s", fname);

    /* find number of bytes in file and write them into buffer */
    nbuffer = nbytes_infile(fp);
    buffer = cmalloc(nbuffer);
    fread(buffer, sizeof(char), nbuffer, fp);

    fclose_buf(fp, &IObuf);
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
      pat = add_patch(mesh, bbox, pt_typ, n, -1);
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
      //if(strcmp(par, "nmax")==0) nmax = atoi(val);
      if(strcmp(par, "rnode->n[0]")==0) n[0] = atoi(val);
      if(strcmp(par, "rnode->n[1]")==0) n[1] = atoi(val);
      if(strcmp(par, "rnode->n[2]")==0) n[2] = atoi(val);
      if(strcmp(par, "rnode->pt_typ[0]")==0) pt_typ[0] = atoi(val);
      if(strcmp(par, "rnode->pt_typ[1]")==0) pt_typ[1] = atoi(val);
      if(strcmp(par, "rnode->pt_typ[2]")==0) pt_typ[2] = atoi(val);
      if(strcmp(par, "s[0]")==0)  pat->CI->s[0] = atof(val);
      if(strcmp(par, "s[1]")==0)  pat->CI->s[1] = atof(val);
      if(strcmp(par, "s[2]")==0)  pat->CI->s[2] = atof(val);
      if(strcmp(par, "s[3]")==0)  pat->CI->s[3] = atof(val);
      if(strcmp(par, "s[4]")==0)  pat->CI->s[4] = atof(val);
      if(strcmp(par, "s[5]")==0)  pat->CI->s[5] = atof(val);
      if(strcmp(par, "xc[0]")==0) pat->CI->xc[0] = atof(val);
      if(strcmp(par, "xc[1]")==0) pat->CI->xc[1] = atof(val);
      if(strcmp(par, "xc[2]")==0) pat->CI->xc[2] = atof(val);
      if(strcmp(par, "co[0]")==0) pat->CI->co[0] = atof(val);
      if(strcmp(par, "co[1]")==0) pat->CI->co[1] = atof(val);
      if(strcmp(par, "co[2]")==0) pat->CI->co[2] = atof(val);
      if(strcmp(par, "co[3]")==0) pat->CI->co[3] = atof(val);
      if(strcmp(par, "dom")==0)   pat->CI->dom = atoi(val);
      if(strcmp(par, "type")==0)  pat->CI->type = atoi(val);
      if(strcmp(par, "use_FSurf")==0)  useF = atoi(val);
      if(strcmp(par, "label")==0 ||
         strcmp(par, "coordinates_get_label(pat)")==0)
      {
        /* this signifies end of patch info, so now set final pat info */
        int clabel = atoi(val);
        pat->CI->label = clabel;
        switch(clabel)
        {
        case CubedSphere:
          set_1_CubedSphere_pat(pat, 0, useF);
          break;
        case CubedSphere_Stretch1:
          set_1_CubedSphere_pat(pat, 1, useF);
          break;
        case CubedSphere_Stretch2:
          set_1_CubedSphere_pat(pat, 2, useF);
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
  int IObufsz = Geti(Par("fread_bufsize"));
  char *IObuf;

  /* alloc buffer */
  buffer = cmalloc(nbuffer_allocd);

  /* open file on rank0 */
  if(Rank0)
  {
    fp = fopen_buf(fname, "rb", &IObuf,IObufsz);
    if(!fp) errorexits("failed opening %s", fname);
  }

  /* get part of file into mem buffer on proc0 */
  file_end = 0;
  do
  {
    char buf[1000];
    char str[1000];
    tNode *parent;
    int cn[8][3];      /* number of points in each child */
    int cp_typ[8][3];  /* point type we set for each child */
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
      int lastone, c;

      /* read a certain number of parent nodes and their child0->n: */
      /* start with reading a certain number of bytes */
      nbuffer = fread(buffer, sizeof(char), first, fp);
      if(nbuffer<first) file_end = 1;

      /* read more bytes until info for last parent in buffer is complete */
      tailbuf = buffer + nbuffer;
      pos = 0;
      lastone = 0;
      while((c=fgetc(fp)) != EOF)
      {
        tailbuf[pos++] = c;
        if(c=='_') lastone = 1; /* from now on we read until we get \n\n */
        if(lastone && c=='\n' && tailbuf[pos-2]=='\n') break;
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
    /* now use the info in buffer to create nodes */
    {
      long off, len;

      /* read buffer line by line */
      off = 0;
      while((off = str_from_buf(buffer,nbuffer, off, '\n', buf,999, &len))>=0)
      {
        /* all node names contain an '_' */
        if(strstr(buf, "_"))
        {
          int chld;
          int d;
          tNlist *children, *lastchild, *child;

          /* strip trailing '\n' from buf */
          buf[strlen(buf)-1] = 0;

          /* find parent node from its name in buf */
          parent = node_from_nodename(mesh, buf);
          p_prev = p;
          p = parent->pat->p;
          lp_prev = lp;
          lp = Elm_l(parent);
          //printf("buf=%s\n", buf);
          //printnode(parent);

          /* init children n and pt_type to same as parent */
          for(chld=0; chld<8; chld++)
            for(d=0; d<3; d++)
            {
              cn[chld][d]     = parent->n[d];
              cp_typ[chld][d] = parent->pt_typ[d];
            }

          /* read n and point type for child0, and potentially other
             children */
          for(chld=0; chld<8; chld++)
          {
            int ret;
            for(d=0; d<3; d++)
            {
              ret = str_from_buf(buffer,nbuffer, off, '\n', str,999, &len);
              if(str[0]=='\n') /* if we get \n, read nothing */
              {
                cn[chld][d]     = cn[0][d];     /* take info from child0 */
                cp_typ[chld][d] = cp_typ[0][d];
                continue;
              }
              off = ret; /* advance buffer offset off */
              ret = sscanf(str, "%d%d", &(cn[chld][d]), &(cp_typ[chld][d]));
              if(ret<2) cp_typ[chld][d] = cp_typ[0][d];
            }
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
          //FIXME:
          //children = make8_child_nodes(parent, cp_typ[0], cn[0]);
          errorexit("call new elm ref function");



          /* set n and pt_typ for each child if needed */
          chld=0;
          fornodelist(children, child)
          {
            tNode *cnode = child->node;
            int update;
            for(update=0, d=0; d<3; d++)
              if(cnode->n[d]      != cn[chld][d] ||
                 cnode->pt_typ[d] != cp_typ[chld][d]) { update=1; break; }
            if(update) update_node_n_pt_typ(cnode, cn[chld], cp_typ[chld]);
            chld++;
          }

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
  if(Rank0) fclose_buf(fp, &IObuf);

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
  tVarList *vl;
  int nvars=0;
  int *ibuffer;
  FILE *fp=NULL;
  char *buffer;
  long nbuffer;
  int i;
  int IObufsz = Geti(Par("fread_bufsize"));
  char *IObuf; /* larger buffer for file read */

  /* open file on rank0 */
  if(Rank0)
  {
    fp = fopen_buf(fname, "rb", &IObuf,IObufsz);
    if(!fp) errorexits("failed opening %s", fname);

    /* get varlist from file */
    vl = checkpoint_make_vl(fp, mesh);
    nvars = vl->n;
  }
  else
  {
    vl = vlalloc(mesh);
  }

  /* broadcast vl to all MPI ranks */
  nMPI_Bcast(&nvars,1, nMPI_INT, 0);
  ibuffer = imalloc(nvars);
  if(Rank0)
    for(i=0; i<nvars; i++) ibuffer[i] = vl->index[i];
  nMPI_Bcast(ibuffer,nvars, nMPI_INT, 0);
  if(!Rank0)
    for(i=0; i<nvars; i++) vlpushone(vl, ibuffer[i]);
  free(ibuffer);

  /* loop over all nodes */
  buffer = NULL;
  do
  {
    int datrank;

    /* make sure to free everything in buffer */
    free(buffer);
    buffer = NULL;

    /* get data (in little endian format) for 1 node into buffer on rank0 */
    if(Rank0)
      buffer = checkpoint_make_nodebuffer(fp, vl, 0, &nbuffer, &datrank);
/*
if(Rank0)
{
FILE *out=fopen("out", "w");
fwrite(buffer,1,nbuffer, out);
fclose(out);
}
nMPI_barrier();
exit(9);
*/
    /* broadcast nbuffer and stop if it's empty */
    nMPI_Bcast(&nbuffer,1, nMPI_LONG, 0);
    nMPI_Bcast(&datrank,1, nMPI_INT, 0);
    //PRF;printf(": nbuffer=%ld datrank=%d\n", nbuffer, datrank);
    if(!nbuffer) break; /* break do-loop if no more node-data */

    /* if the buffer needs to go to another rank, send it there */
    if(datrank != 0)
    {
      /* send buffer from 0 to rank with dat */
      if(Rank0)
        nMPI_Send(buffer, nbuffer, nMPI_CHAR, datrank, 11);

      if(nMPI_rank()==datrank)
      {
        buffer = cmalloc(nbuffer);
        nMPI_Recv(buffer, nbuffer, nMPI_CHAR, 0, 11);
      }
    }

    /* read var data from buffer for this one node */
    if(nMPI_rank()==datrank) checkpoint_read_vl(buffer, nbuffer, vl);

    /* wait for all to get here */
    nMPI_barrier();
  } while(nbuffer);
  free(buffer);

  PRF;printf(": finished reading varlist.\n");
  prvarlist(vl);
  fflush(stdout);

  vlfree(vl);
  if(Rank0) fclose_buf(fp, &IObuf);
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

/* read varlist on each node */
void checkpoint_read_vl(char *buffer, long nbuffer, tVarList *vl)
{
  tMesh *mesh = vl->mesh;
  char buf[1000];
  long off, len;

  /* find string node */
  off = 0;
  while((off = str_from_buf(buffer,nbuffer, off, '\n', buf,999, &len))>=0)
  {
    tNode *node;
    char name[256];
    int np, vli, vi, found_node;

    if(strcmp(buf, "{\n")==0)
    {
      /* read node info */
      off = str_from_buf(buffer,nbuffer, off, '\n', buf,999, &len);
      sscanf(buf, "%s", name);
      off = str_from_buf(buffer,nbuffer, off, '\n', buf,999, &len);
      sscanf(buf, "%d", &np);

      /* find node from its name */
      node = node_from_nodename(mesh, name);
      found_node = 1;
    }
    else
    {
      found_node = 0;
    }
    while(found_node && (off>=0))
    {
      /* check for end ("}\n") or read var info */
      off = str_from_buf(buffer,nbuffer, off, '\n', buf,999, &len);
      if(strcmp(buf, "}\n")==0) break;
      vli = atoi(buf);
      vi  = Vind(vl, vli);
      len = sizeof(double) * np;
      //PRF;printf(": name=%s np=%d vli=%d node->datrank=%d\n",
      //           name, np, vli, node->datrank);
      //printf("%s nid=%ld vi=%d\n", nodename(node,buf,99), Node_eid(node), vi);

      /* check if this var needs to be read on this node */
      if(node->dat)
      {
        double *v;

        /* switch on var on this node and get pointer to its data */
        enablevarcomp_innode(node, vi);
        v = Vard(node, vi);

        /* read var as raw binary */
        memcpy(v, buffer+off, len);
        //if(Node_eid(node)==28) printf("v[]=%g\n", v[0]);
      }
      off += len; /* go len bytes further in buffer */

      /* also read the '\n' after var. data */
      off = str_from_buf(buffer,nbuffer, off, '\n', buf,999, &len);
    } /* end while(found_node && (off>=0)) */
  }
}

/* Add buf to buffer. To be used e.g. like this:
     buffer = append_buf(buffer,nbuffer, buf,strlen(buf));  */
char *append_buf(char *buffer, long *nbuffer, const char *buf, long nbuf)
{
  long len = *nbuffer + nbuf;
  buffer = realloc(buffer, len * sizeof(char));
  memcpy(buffer + *nbuffer, buf, nbuf);
  *nbuffer = len;
  return buffer;
}

/* read var info for one node into buffer */
char *checkpoint_make_nodebuffer(FILE *fp, tVarList *vl, int read_big,
                                 long *nbuffer, int *datrank)
{
  tMesh *mesh = vl->mesh;
  char *buffer;
  char buf[1000];
  char *s;
  double *v = NULL;

  /* read line by line into the buffer */
  buffer = NULL;
  *nbuffer = 0;
  while((s = fgets(buf,999, fp)))
  {
    char name[256];
    int np, found_node;
    tNode *node;

    if(strcmp(buf, "{\n")==0)
    {
      buffer = append_buf(buffer,nbuffer, buf,strlen(buf)); /* app "{\n" */

      fgets(buf,999, fp);
      buffer = append_buf(buffer,nbuffer, buf,strlen(buf)); /* app "nodename\n" */
      sscanf(buf, "%s", name);   /* find node name */
      node = node_from_nodename(mesh, name);
      *datrank = node->datrank;

      s = fgets(buf,999, fp); /* s=NULL at EOF */
      buffer = append_buf(buffer,nbuffer, buf,strlen(buf)); /* app "np\n" */
      np = atoi(buf);

      found_node = 1;
      free(v);
      v = dmalloc(np);
    }
    else
    {
      found_node = 0;
    }
    while(found_node && s)
    {
      //PRF;printf(": %s %d found_node=%d\n", name, np, found_node);
      /* check for end / read var info */
      fgets(buf,999, fp); /* use fgets to read "}\n" or vli plus '\n' */
      buffer = append_buf(buffer,nbuffer, buf,strlen(buf)); /* app buf */
      if(strcmp(buf, "}\n")==0) break;

      /* read var as raw binary */
      if(read_big) fread_big(v, sizeof(double), np, fp);
      else         fread_little(v, sizeof(double), np, fp);
      buffer = append_buf(buffer,nbuffer, (char *) v,np*sizeof(*v)); /* app v */

      /* use fgets to also read the '\n' after var. */
      s = fgets(buf,999, fp); /* s=NULL at EOF */
      buffer = append_buf(buffer,nbuffer, buf,strlen(buf)); /* app buf */
    } /* end while(found_node && s) */
    free(v);
    v = NULL;
    if(found_node) break; /* stop after we found a node */
  }
  return buffer;
}
