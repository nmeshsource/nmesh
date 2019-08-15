/* checkpoint.c */
/* Wolfgang Tichy, 8/2019 */


#include "nmesh.h"
#include "checkpoint.h"


/* checkpoint filenames */
char chckpt_dir[] = "checkpoint";
char save_pars_file[] = "save_pars.txt";
char patches_file[]   = "patches.txt";
char nodes_file[]     = "nodes.txt";
char variables_file[] = "variables.bin";


/******************************************************************/
/* some functions to save nmesh data for checkpoints  */
/******************************************************************/

/* write a checkpoint */
int checkpoint_save(tMesh *mesh)
{
  char *outdir = Gets(Par("outdir"));

  int cl = strlen(outdir) + 20;
  char *cdir = cmalloc(cl);

  int pl = cl + 40;
  char *pars  = cmalloc(pl);
  char *pats  = cmalloc(pl);
  char *nodes = cmalloc(pl);
  char *vars  = cmalloc(pl);

  /* output filenames */
  snprintf(cdir,cl, "%s/%s", outdir, chckpt_dir);
  snprintf(pars,pl,  "%s/%s", cdir, save_pars_file);
  snprintf(pats,pl,  "%s/%s", cdir, patches_file);
  snprintf(nodes,pl, "%s/%s", cdir, nodes_file);
  snprintf(vars,pl,  "%s/%s", cdir, variables_file);

  /* save checkpoint in various files */
  if(Rank0) checkpoint_save_patches(mesh, pats);
  if(Rank0) checkpoint_save_nodes(mesh, nodes);
  checkpoint_save_EvoVars(mesh, vars);

  /* free strings */
  free(vars);
  free(nodes);
  free(pats);
  free(pars);
  free(cdir);
  return 0;
}

/* save a checkpoint if the time is right for it */
int checkpoint_save_if_needed(tMesh *mesh)
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
    if((hours      > 0. && hours      <= time_since_checkpoint) ||
       (hours_quit > 0. && hours_quit <= time))
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
  }

  return 0;
}
