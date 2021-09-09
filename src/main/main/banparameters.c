/* banparameters.c */
/* Wolfgang Tichy, 9/2021 */


#include "nmesh.h"


/* structure for one banned parameter */
typedef struct tBANNEDPAR {
  const char *name;
  const char *ban_reason;
} tBannedPar;

/***************************************************************************/
/* compile lists with entries of type (tBannedPar *) */
/***************************************************************************/
typedef tBannedPar *pBanned;  /* list_templates.h only works with numbers */
#define TYP pBanned           /* the pointer pBanned is a number */
#include "list_templates.c"
#undef TYP


/* global list with banned pars */
pBannedList *bannedpars = NULL;


/**************************************************************************/
/* functions to deal with parameters that should not be used anymore,
   because they have been removed, or are considered bad somehow */
/**************************************************************************/

void prFreeBannedParList(void)
{
  int k;
  if(bannedpars)
    forList(bannedpars, k)
      printf("%s %s\n", ListEntry(bannedpars,k)->name,
                        ListEntry(bannedpars,k)->ban_reason);
}

/* compare name in name_ob with name in bp,
   return 1 if names are equal, return 0 otherwise */
int prop_samenames(const void *name_ob, tBannedPar *bp)
{
  const char *name = (const char *) name_ob;
  printf("strcmp(%s, %s);\n", bp->name, name);
  return ( strcmp(name, bp->name)==0 );
}

/* put an outdated par in a blacklist */
void BanPar(const char *name, const char *ban_reason)
{
  tBannedPar *bp;
  int i;

  if(!bannedpars) bannedpars = alloc_pBannedList();

  if(0) {PRF;printf("%s %s\n", name, ban_reason);}

  /* look for name in bannedpars */
  i = index_prop_pBannedList(bannedpars,0, prop_samenames, name);
  if(i>=0)
    errorexit("cannot ban the same par twice");

  /* add par name to bannedpars */
  bp = calloc(1, 1);
  bp->name = name;
  bp->ban_reason = ban_reason;
  push_pBannedList(bannedpars, bp);
}


/* look for par name in bannedpars */
int ExitIfParBanned(const char *name)
{
  //prFreeBannedParList();
  if(bannedpars)
  {
    /* look for name in bannedpars */
    int i = index_prop_pBannedList(bannedpars,0, prop_samenames, name);
    PRF;printf(": %s -> i=%d\n", name, i);
    if(i>=0)
    {
      printf("Parameter %s is banned:\n", name);
      printf("%s\n", bannedpars->e[i]->ban_reason);
      errorexits("Parameter %s is banned", name);
    }
    /* i=-1 means par name was not found */
  }
  return 0;
}

/* free all in bannedpars */
int FreeBannedParList(tMesh *mesh)
{
  freeall_pBannedList(bannedpars, free, NULL);
  bannedpars = NULL;
  return 0;
}
