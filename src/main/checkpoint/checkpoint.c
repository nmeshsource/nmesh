/* checkpoint.c */
/* Wolfgang Tichy, 8/2019 */


#include "nmesh.h"
#include "checkpoint.h"


/* checkpoint filenames */
char chckpt_dir[] = "checkpoint";
char save_pars_file[] = "save_pars.par";
char patches_file[]   = "patches.txt";
char elms_file[]      = "elms.txt";
char nbinfo_file[]    = "nbinfo.bin";
char variables_file[] = "variables.bin";

/* make and allocate complete checkpoint pathnames, return allocated chars */
int checkpoint_create_pathnames(tMesh *mesh, const char *outdir_suffix,
                                char **Dir, const char *Dir_suffix,
                                char **Pars, char **Pats,
                                char **Elms, char **Nbinfo, char **Vars)
{
  char *outdir = Gets(Par("outdir"));

  int nl = strlen(outdir) + 40;
  char *dir = cmalloc(nl);

  int pl = nl + 40;
  char *pars  = cmalloc(pl);
  char *pats  = cmalloc(pl);
  char *elms  = cmalloc(pl);
  char *nbinfo= cmalloc(pl);
  char *vars  = cmalloc(pl);

  /* filenames */
  snprintf(dir,nl, "%s%s/%s%s", outdir,outdir_suffix, chckpt_dir,Dir_suffix);
  snprintf(pars,pl,  "%s/%s", dir, save_pars_file);
  snprintf(pats,pl,  "%s/%s", dir, patches_file);
  snprintf(elms,pl,  "%s/%s", dir, elms_file);
  snprintf(nbinfo,pl,"%s/%s", dir, nbinfo_file);
  snprintf(vars,pl,  "%s/%s", dir, variables_file);

  /* save filenames */
  *Dir   = dir;
  *Pars  = pars;
  *Pats  = pats;
  *Elms  = elms;
  *Nbinfo= nbinfo;
  *Vars  = vars;
  return pl;
}

/* return 0 if there is no checkpoint,
   return ret otherwise:
   ret has bit1=2^1 set if checkpoint has patches_file
   ret has bit3=2^4 set if checkpoint has variables_file */
int checkpoint_exists(tMesh *mesh, const char *outdir_suffix,
                      const char *Dir_suffix)
{
  char *dir;
  char *pars;
  char *pats;
  char *elms;
  char *nbinfo;
  char *vars;
  int ret=0;

  checkpoint_create_pathnames(mesh, outdir_suffix, &dir, Dir_suffix,
                              &pars, &pats, &elms, &nbinfo, &vars);
  /* check if files exist */
  if(Rank0)
  {
    FILE *fp;

    fp = fopen(pats, "r");
    if(fp)
    {
      ret |= CHECKPOINT_PATS;
      fclose(fp);
    }

    fp = fopen(nbinfo, "r");
    if(fp)
    {
      ret |= CHECKPOINT_NBINFO;
      fclose(fp);
    }

    fp = fopen(vars, "r");
    if(fp)
    {
      ret |= CHECKPOINT_VARS;
      fclose(fp);
    }
  }
  /* broadcast ret from rank0 to all others */
  nMPI_Bcast(&ret, 1, nMPI_INT, 0);

  /* free strings */
  free(vars);
  free(nbinfo);
  free(elms);
  free(pats);
  free(pars);
  free(dir);
  return ret;
}

/******************************************************************/
/* some functions to load nmesh data from checkpoints  */
/******************************************************************/

/* read a checkpoint
   stage=0 loads only the mesh, stage>0 loads only the vars */
int checkpoint_load_stage(tMesh *mesh, const char *outdir_suffix,
                          int stage)
{
  int checkpoint = Par("checkpoint");
  int chkpt_on;
  char *dir;
  char *pars;
  char *pats;
  char *elms;
  char *nbinfo;
  char *vars;
  int ret=0;
  double time, ntime;

  /* is checkpointing on? */
  chkpt_on = (Getb(checkpoint) || Getv(checkpoint, "load_mesh"));
  if(!chkpt_on) return 0;

  time = getTimeIn_s()/60.;

  checkpoint_create_pathnames(mesh, outdir_suffix, &dir,"",
                              &pars, &pats, &elms, &nbinfo, &vars);
  prdivider(1);
  PRF;printf(": loading stage%d\n", stage);
  fflush(stdout);

  /* load checkpoint from the various files */
  if(stage==0)
  {
    int chkpt_exists = checkpoint_exists(mesh, "", "");

    nmesh_load_parameters(mesh, pars, 0, 1);
    PRF;printf(": finished loading pars.\n");
    fflush(stdout);
    checkpoint_load_patches(mesh, pats);
    PRF;printf(": finished loading patches.\n");
    fflush(stdout);
    checkpoint_load_elms(mesh, elms);
    PRF;printf(": finished loading elms.\n");
    fflush(stdout);
    if(chkpt_exists & 8)
    {
      checkpoint_load_Vars(mesh, nbinfo, 1); /* load with native endianess */
      PRF;printf(": finished loading nbinfo.\n");
      checkpoint_set_nbinfo_fnb_nbelm_loadbal(mesh, 0);
    }
    else
    {
      checkpoint_set_nbinfo_fnb_nbelm_loadbal(mesh, 1);
    }
    fflush(stdout);
  }
  else
  {
    checkpoint_load_Vars(mesh, vars, 0); /* load as little endian */
    PRF;printf(": finished loading variables.\n");
    fflush(stdout);
  }
  ntime = getTimeIn_s()/60.;
  PRF;printf(": finished stage%d in %g minutes.\n", stage, ntime-time);
  prdivider(1);
  fflush(stdout);

  /* free strings */
  free(vars);
  free(nbinfo);
  free(elms);
  free(pats);
  free(pars);
  free(dir);
  return ret;
}

/* setup nb related stuff and do load bal */
void checkpoint_set_nbinfo_fnb_nbelm_loadbal(tMesh *mesh, int reset_nbinfo)
{
  //PRFs(": 1.\n");
  //printmesh(mesh);

  /* either keep the nbinfo we loaded or recreate it*/
  if(reset_nbinfo)
    amr_update_elm_nbinfo_if_nnbinfo_negative(mesh);

  /* set fnb stuff */
  amr_elm_nbinfo_to_elm_fnb(mesh);
  amr_elm_nbinfo_set_nnbinfo_mesh(mesh, 1); //make nnbinfo positive */
  amr_get_nbelm_elmheaders(mesh);

  /* load balance all leaf nodes */
  simple_load_balance(mesh);
  amr_elm_nbinfo_set_nnbinfo_mesh(mesh, 1); //make nnbinfo positive */

  //PRFs(": 2.\n");
  //printmesh(mesh);
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
  char *elms;
  char *nbinfo;
  char *vars;
  int pl = checkpoint_create_pathnames(mesh, "", &dirn,"_new",
                                       &pars, &pats, &elms, &nbinfo, &vars);
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
  checkpoint_save_elms(mesh, elms);
  checkpoint_save_nbinfoVars(mesh, nbinfo);
  checkpoint_save_EvoVars(mesh, vars);

  /* wait until all get here */
  nMPI_barrier();

  /* rename checkpoint and remove old one */
  if(Rank0)
  {
    system2("rm -rf", dirp);
    system3("mv", dir, dirp);
    system3("mv", dirn, dir);
    checkpoint_keep_previous(outdir, pl, dirp,
                             Geti(Par("checkpoint_n_previous")));
  }

  /* free strings */
  free(dirp);
  free(dir);
  free(nbinfo);
  free(vars);
  free(elms);
  free(pats);
  free(pars);
  free(dirn);
  return 0;
}

/* keep nprev previous checkpoints */
void checkpoint_keep_previous(char *outdir, int pl, char *dirp, int nprev)
{
  char *dir_i   = cmalloc(pl+12);
  char *dir_ip1 = cmalloc(pl+12);
  int i;

  /* if we keep 1 or more, remove oldest in case it's there due to errors */
  if(nprev>0)
  {
    snprintf(dir_ip1,pl, "%s/%s-0LDE5T.tmp", outdir, chckpt_dir);
    system2("rm -rf", dir_ip1);
  }

  /* increase suffix part of all previous checkpoints, e.g.:
     checkpoint-1 -> checkpoint-2 */
  for(i=nprev; i>0; i--)
  {
    snprintf(dir_i,  pl, "%s/%s-%d", outdir, chckpt_dir, i);
    if(i==nprev) snprintf(dir_ip1,pl, "%s/%s-0LDE5T.tmp", outdir, chckpt_dir);
    else         snprintf(dir_ip1,pl, "%s/%s-%d", outdir, chckpt_dir, i+1);
    system3("mv", dir_i, dir_ip1);
  }

  if(nprev>0)
  {
    /* rename checkpoint_previous to checkpoint-1 */
    snprintf(dir_i,  pl, "%s/%s-%d", outdir, chckpt_dir, 1);
    system3("mv", dirp, dir_i);
    /* if we keep 1 or more, remove oldest one */
    snprintf(dir_ip1,pl, "%s/%s-0LDE5T.tmp", outdir, chckpt_dir);
    system2("rm -rf", dir_ip1);
  }
  else
  {
    /* if nprev=0 we keep nothing and thus delete checkpoint_previous */
    system2("rm -rf", dirp);
  }

  free(dir_ip1);
  free(dir_i);
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

  /* is checkpointing on? */
  if(!Getb(Par("checkpoint"))) return 0;

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

  /* broadcast do_checkpoint, last_checkpoint_time from rank0 to all others */
  nMPI_Bcast(&do_checkpoint, 1, nMPI_INT, 0);
  nMPI_Bcast(&last_checkpoint_time, 1, nMPI_DOUBLE, 0);

  /* now do it if needed */
  if(do_checkpoint)
  {
    double ntime;

    prdivider(1);
    printf("checkpoint_save after %g hours\n", time_since_checkpoint);
    fflush(stdout);

    checkpoint_save(mesh);
    ntime = getTimeIn_s()/3600.;

    printf("It took %g minutes to checkpoint %g hours into the run.\n",
           60.*(ntime - time), time);
    prdivider(1);
    fflush(stdout);

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

      printf("Now exit nmesh before the queuing system kills it!\n");
      RunFun(FINALIZE); /* hook for funcs to run for a clean return */
      prTimeIn_s("WallTime just before finalize_all_and_exit: ");
      finalize_all_and_exit(mesh, 0);
    }
    return 1; /* signal that checkpoint was saved */
  }

  return 0;
}
