/* checkpoint.c */
/* Wolfgang Tichy, 8/2019 */


#include "nmesh.h"
#include "checkpoint.h"


/* checkpoint filenames */
char chckpt_dir[] = "checkpoint";
char save_pars_file[] = "save_pars.par";
char patches_file[]   = "patches.txt";
char nodes_file[]     = "nodes.txt";
char variables_file[] = "variables.bin";

/* make and allocate complete checkpoint pathnames, return allocated chars */
int checkpoint_create_pathnames(tMesh *mesh, char **Dir, const char *suffix,
                                char **Pars, char **Pats,
                                char **Nodes, char **Vars)
{
  char *outdir = Gets(Par("outdir"));

  int nl = strlen(outdir) + 40;
  char *dir = cmalloc(nl);

  int pl = nl + 40;
  char *pars  = cmalloc(pl);
  char *pats  = cmalloc(pl);
  char *nodes = cmalloc(pl);
  char *vars  = cmalloc(pl);

  /* filenames */
  snprintf(dir,nl, "%s/%s%s", outdir, chckpt_dir, suffix);
  snprintf(pars,pl,  "%s/%s", dir, save_pars_file);
  snprintf(pats,pl,  "%s/%s", dir, patches_file);
  snprintf(nodes,pl, "%s/%s", dir, nodes_file);
  snprintf(vars,pl,  "%s/%s", dir, variables_file);

  /* save filenames */
  *Dir   = dir;
  *Pars  = pars;
  *Pats  = pats;
  *Nodes = nodes;
  *Vars  = vars;
  return pl;
}


/******************************************************************/
/* some functions to load nmesh data from checkpoints  */
/******************************************************************/

/* read a checkpoint */
int checkpoint_load(tMesh *mesh)
{
  char *dir;
  char *pars;
  char *pats;
  char *nodes;
  char *vars;

  checkpoint_create_pathnames(mesh, &dir,"", &pars, &pats, &nodes, &vars);

  /* load checkpoint from the various files */
  checkpoint_load_pars(mesh, pars);
  //checkpoint_load_patches(mesh, pats);
  //checkpoint_load_nodes(mesh, nodes);
  //checkpoint_load_EvoVars(mesh, vars);

  /* free strings */
  free(vars);
  free(nodes);
  free(pats);
  free(pars);
  free(dir);
  return 0;
}

/******************************************************************/
/* some functions to save nmesh data for checkpoints  */
/******************************************************************/

/* write a checkpoint */
int checkpoint_save(tMesh *mesh)
{
  char *outdir = Gets(Par("outdir"));
  char *dirn;
  char *pars;
  char *pats;
  char *nodes;
  char *vars;
  int pl = checkpoint_create_pathnames(mesh, &dirn,"_new",
                                       &pars, &pats, &nodes, &vars);
  char *dir = cmalloc(pl);
  char *dirp = cmalloc(pl);

  /* output filenames */
  snprintf(dir,pl, "%s/%s", outdir, chckpt_dir);
  snprintf(dirp,pl, "%s_previous", dir);

  /* make new dir */
  if(Rank0) system2("mkdir", dirn);

  /* save checkpoint in various files */
  checkpoint_save_pars(mesh, pars);
  checkpoint_save_patches(mesh, pats);
  checkpoint_save_nodes(mesh, nodes);
  checkpoint_save_EvoVars(mesh, vars);

  /* wait until all get here */
  nMPI_barrier();

  /* rename checkpoint and remove old one */
  if(Rank0)
  {
    system3("mv", dir, dirp);
    system3("mv", dirn, dir);
    system2("rm -rf", dirp);
  }

  /* free strings */
  free(dirp);
  free(dir);
  free(vars);
  free(nodes);
  free(pats);
  free(pars);
  free(dirn);
  return 0;
}

/* save a checkpoint if the time is right for it */
int checkpoint_save_if_needed(tMesh *mesh, int always)
{
  static double last_checkpoint_time = 0.;
  double hours      = Getd(Par("checkpoint_hours"));
  double hours_quit = Getd(Par("checkpoint_hours_quit"));
  double time       = getTimeIn_s()/3600.;
  double time_since_checkpoint;
  int do_checkpoint = 0;

  /* test if it is time */
  time_since_checkpoint = time - last_checkpoint_time;
  if(Rank0)
  {
    /* test based on walltime */
    if((hours      >= 0. && hours      <= time_since_checkpoint) ||
       (hours_quit >= 0. && hours_quit <= time) ||
       (hours      >= 0. && always))
    {
      do_checkpoint = 1; /* yes, we want to save a checkpoint */
      last_checkpoint_time = time;
    }
  }

  /* broadcast do_checkpoint from rank0 to all others */
  nMPI_Bcast(&do_checkpoint, 1, nMPI_INT, 0);

  /* now do it if needed */
  if(do_checkpoint)
  {
    double ntime = getTimeIn_s()/3600.;

    nMPI_Bcast(&ntime, 1, nMPI_INT, 0);
    printf("checkpoint_save by walltime, after %g hours\n",
           time_since_checkpoint);

    checkpoint_save(mesh);

    printf("It took %g minutes to checkpoint %g hours into the run.\n",
           60.*(ntime - time), time);

    /* check if we want to kill nmesh */
    if(0. < hours_quit && hours_quit <= ntime)
    {
      /*
      if(Rank0 && Getb(Par("checkpoint_resub")))
      {
        printf("resub simulation\n");
        system2("", Gets(Par("checkpoint_resub_command")));
      }
      */

      printf("Now kill nmesh before the queuing system does!\n");
      //finalize_mesh(mesh);
      //nMPI_Comm_free(&(main_comm));
      nMPI_Finalize();
    }
    return 1; /* signal that checkpoint was saved */
  }

  return 0;
}
