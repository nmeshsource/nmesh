/* output.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "output.h"


/* output par prefixes */
char *outpre[] = {"0d", "1d", "2d", "3d", "co"};
int Nout = sizeof(outpre)/sizeof(outpre[0]); /* Nout=5 */

/* various types of output */
char *outpt[] = {"pt"}; /* prefix of the Noutpt special 0doutput points */
char *out0d[] = {"max", "min", "maxAbs", "VolInt", "rms", "meanAbs"};
char *out1d[] = {"X", "Y", "Z"};
char *out2d[] = {"XY", "XZ", "YZ"};
char *out3d[] = {"XYZ"};
char *outco[] = {"co"};

#define LEN0d sizeof(out0d)/sizeof(out0d[0])
#define LEN1d sizeof(out1d)/sizeof(out1d[0])
#define LEN2d sizeof(out2d)/sizeof(out2d[0])
#define LEN3d sizeof(out3d)/sizeof(out3d[0])
#define LENco sizeof(outco)/sizeof(outco[0])

/* number of output types (Noutpt is in output->Noutpt) */
int Nout0d = LEN0d;
int Nout1d = LEN1d;
int Nout2d = LEN2d;
int Nout3d = LEN3d;
int Noutco = LENco;
//int NoutAll = LEN0d_x + LEN0d + LEN1d + LEN2d + LEN3d + LENco;




/* is it time to output all nodes? */
int TimeForMeshOutput_di_dt(tMesh *mesh, int di, double dt) 
{
  double Time = fabs(mesh->time);
  double dT = fabs(mesh->dt);
  long Iter = mesh->iteration;

  /* Note: this function could check if all nodes are aligned in time and then
     call TimeForNodeOutput_di_dt to test if all nodes want to output */

  /* time for output based on number of iterations */
  if(di > 0 && Iter % di == 0)
    return Iter/di + 1;

  /* time for output based on time, assumes Time >= 0 */
  if(dt > 0)
  {
    int i = (0.5 + Time/dt); /* (Time + dequaleps)/dt; */
    /* if(dequal(Time-i*dt, 0)) */
    if( Time-i*dt > -0.5*dT && Time-i*dt <= 0.5*dT )
      return i + 1;
  }
    
  /* not time for output */
  return 0;
}



/* registered function that calls write_mesh to do output */
int mesh_output(tMesh *mesh)
{
  return write_mesh(mesh, mesh->iteration, mesh->time);
}

/* function that controls all writing, can be called directly, and
   generates files with Iter & Time as labels, no matter
   what iteration & time the mesh has */
int write_mesh(tMesh *mesh, int Iter, double Time)
{
  int di[Nout];
  double dt[Nout];
  char *ou[Nout];
  int all[Nout];
  char str[1000];
  int start;
  int d, vi, vi0;
  tVarList *vl[Nout];    /* varlists for 0d,1d,2d,3d,... */

  /* varlists of 0d,1d,2d,3d output */  
  for(d = 0; d < Nout; d++) vl[d] = vlalloc(mesh);

  /* par values in strings */
  for(d = 0; d < Nout; d++)
  {
    snprintf(str,999, "%soutiter", outpre[d]);
    di[d] = Geti(Par(str));
    snprintf(str,999, "%souttime", outpre[d]);
    dt[d] = Getd(Par(str));
    snprintf(str,999, "%soutput", outpre[d]);
    ou[d] = Gets(Par(str));
    snprintf(str,999, "%soutputall", outpre[d]);
    all[d] = Getb(Par(str));
  }

  /* d=0: 0d output, d=1: 1d output, d=2: 2d output, ... */
  for(d=0; d<Nout; d++)
  {
    if(TimeForMeshOutput_di_dt(mesh, di[d], dt[d]))
    {
      //printf("2dout ... |%s|\n", ou[d]);
      start=0;
      while(sscanf(ou[d]+start, "%s", str)==1)
      {
        start += strlen(str);
        if(ou[d][start]==' ') start++;
        //printf("2dout |%s|\n", str);

        /* check if str has an index that exists */
        vi = MeshVarIndLax(mesh, str);
        if(vi<0) continue;

        /* ?doutputall */
        if(all[d])
        {
          tVarList *vltmp = vlalloc(mesh);
          intList *iltmp  = vl2intList(vltmp);
          intList *il     = vl2intList(vl[d]);

          vi0 = MeshVarIndComponent0(mesh, vi);
          vlpush(vltmp, vi0);
          /* use unionpushlist_intList to add all in vltmp only once */
          unionpushlist_intList(il, iltmp);
          vlfree(vltmp);
        }
        else
        {
          /* use unionpush_intList to add var-comp only once */
          intList *il = vl2intList(vl[d]);
          unionpush_intList(il, vi);
        }
      } //end: while loop
    }
  } // end: for d loop

  /* output the varlists we just created */
  output0d_vl(vl[0], Iter, Time);
  output1d_vl(vl[1], Iter, Time);
  output2d_vl(vl[2], Iter, Time);
  output3d_vl(vl[3], Iter, Time);
  outputco_vl(vl[4], Iter, Time);

  /* free varlists */
  for(d = 0; d < Nout; d++) vlfree(vl[d]);

  return 0;
}

/*******************************************************************/
/* more funcs to test if it is time to output */
/*******************************************************************/

/* is it time to output a variable with index vindex? */
int TimeForMeshOutput_vindex(tMesh *mesh, int vindex)
{
  //int Noutput = NextAll + Nout; // really just need 5 from 0d-,1d-,2d-,... output
  int Noutput = Nout;
  int di[Noutput];
  double dt[Noutput];
  char output[Noutput][128];
  char *name, s[128];
  int d, n;

  //errorexit("this function is not tested yet!");

  name = VarName(vindex);

  /* read all pars for all dims */
  for(d = 0; d < Nout; d++)
  {
    sprintf(s, "%soutiter", outpre[d]);
    di[d] = Geti(Par(s));
    sprintf(s, "%souttime", outpre[d]);
    dt[d] = Getd(Par(s));
    sprintf(output[d], "%soutput", outpre[d]);
  }

  /* check if "name" is contained in any output par */
  for(n = 0; n < Noutput; n++)
  {
    if(Getv(Par(output[n]), name))
      if(TimeForMeshOutput_di_dt(mesh, di[n], dt[n])) return 1;
  }
  return 0;
}

/* is it time to output any variable in a VarList ? */
int TimeForMeshOutput_vl(tMesh *mesh, tVarList *vl)
{
  int i;

  for(i = 0; i < vl->n; i++)
    if(TimeForMeshOutput_vindex(mesh, vl->index[i])) return 1;
  return 0;
}



/*******************************************************************/
/* Everything below this line is untested and may not work
   It may not be needed!!! */
/*******************************************************************/

/* Is it time to output? Yes if Iter is a multiple of di, or 
   if Time has exceeded a multiple of dt */
int TimeForNodeOutput_di_dt(tNode *node, int di, double dt) 
{
  tMesh *mesh = node->pat->mesh;  // NOTE: later each node might have its
  double Time = fabs(mesh->time); // own time
  double dT = fabs(mesh->dt);
  long Iter = mesh->iteration;

  errorexit("this function is not tested yet!");

  /* time for output based on number of iterations */
  if(di > 0 && Iter % di == 0)
    return Iter/di + 1;

  /* time for output based on time, assumes Time >= 0 */
  if(dt > 0)
  {
    int i = (0.5 + Time/dt); /* (Time + dequaleps)/dt; */
    /* if(dequal(Time-i*dt, 0)) */
    if( Time-i*dt > -0.5*dT && Time-i*dt <= 0.5*dT )
      return i + 1;
  }

  /* not time for output */
  return 0;
}

/* is it time to output a variable with index vindex? */
int TimeForNodeOutput_vindex(tNode *node, int vindex)
{
  tMesh *mesh = node->pat->mesh;
  //int Noutput = NextAll + Nout; // really just need Nout from 0d-,1d-,2d-,... output
  int Noutput = Nout;
  int di[Noutput];
  double dt[Noutput];
  char output[Noutput][64];
  char *name, s[64];
  int d, n;

  errorexit("this function is not tested yet!");

  name = VarName(vindex);

  /* read all pars for all dims */
  for(d = 0; d < Nout; d++)
  {
    sprintf(s, "%soutiter", outpre[d]);
    di[d] = Geti(Par(s));
    sprintf(s, "%souttime", outpre[d]);
    dt[d] = Getd(Par(s));
    sprintf(output[d], "%soutput", outpre[d]);
  }

  /* check if "name" is contained in any output par */
  for(n = 0; n < Noutput; n++)
  {
    if(Getv(Par(output[n]), name))
      if(TimeForNodeOutput_di_dt(node, di[n], dt[n])) return 1;
  }
  return 0;
}

/* is it time to output any variable in a VarList ? */
int TimeForNodeOutput_vl(tNode *node, tVarList *vl)
{
  int i;

  errorexit("this function is not tested yet!");

  for(i = 0; i < vl->n; i++)
    if(TimeForNodeOutput_vindex(node, vl->index[i])) return 1;
  return 0;
}

/* is it time for output in any output in any format? */
int TimeForNodeOutput_any(tNode *node)
{
  tMesh *mesh = node->pat->mesh;
  int di;
  double dt;
  char str[1000];
  int d, ret;

  errorexit("this function is not tested yet!");

  ret = 0;
  for(d = 0; d < Nout; d++)
  {
    snprintf(str,999, "%soutiter", outpre[d]);
    di = Geti(Par(str));
    snprintf(str,999, "%souttime", outpre[d]);
    dt = Getd(Par(str));
    ret = ret | TimeForNodeOutput_di_dt(node, di, dt);
  }
  return ret;
}
