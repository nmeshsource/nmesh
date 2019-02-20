/* output.c */
/* Wolfgang Tichy, 2/2019 */

#include "nmesh.h"
#include "output.h"


/* file extensions to indice various types of output*/
char *extn0d[] = {"max", "min", "maxAbs", "VolInt", "rms", "meanAbs"};
char *extn1d[] = {"X", "Y", "Z"};
char *extn2d[] = {"XY", "XZ", "YZ"};
char *extn3d[] = {"XYZ"};

#define LEN0d sizeof(extn0d)/sizeof(extn0d[0])
#define LEN1d sizeof(extn1d)/sizeof(extn1d[0])
#define LEN2d sizeof(extn2d)/sizeof(extn2d[0])
#define LEN3d sizeof(extn1d)/sizeof(extn1d[0])

/* number of extensions */
int Nextn0d = LEN0d;
int Nextn1d = LEN1d;
int Nextn2d = LEN2d;
int Nextn3d = LEN3d;
int NextAll = LEN0d + LEN1d + LEN2d + LEN3d;




/* is it time to output all nodes? */
int TimeForMeshOutput_di_dt(tMesh *mesh, int di, double dt) 
{
  double Time = fabs(mesh->time);
  double dT = fabs(mesh->dt);
  long Iter = mesh->iteration;

  /* Note: this function could if all nodes are aligned in time and then
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
  int di[4];
  double dt[4];
  char *ou[4];
  char s[32];
  char str[1000];
  int start;
  int d, i, vi, vi0;
Yo(1);
  /* par values in strings */
  for(d = 0; d <= 3; d++)
  {
    snprintf(s,999, "%ddoutiter", d);
    di[d] = Geti(Par(s));
    snprintf(s,999, "%ddouttime", d);
    dt[d] = Getd(Par(s));
    snprintf(s,999, "%ddoutput", d);
    ou[d] = Gets(Par(s));
  }

  /* 0d output */
  d = 0;
  if(TimeForMeshOutput_di_dt(mesh, di[d], dt[d]))
  {
    errorexit("0d output not implemented");
  }

  /* 1d output */
  d = 1;
  if(TimeForMeshOutput_di_dt(mesh, di[d], dt[d]))
  {
    errorexit("1d output not implemented");
  }

  /* 2d output */
  d = 2;
  if(TimeForMeshOutput_di_dt(mesh, di[d], dt[d]))
  {
    printf("2dout ... |%s|\n", ou[d]);
    start=0;
    while(sscanf(ou[d]+start, "%s", str)==1)
    {
      start += strlen(str);
      if(ou[d][start]==' ') start++;
      printf("2dout |%s|\n", str);

      /* check if str has an index that exists */
      vi = MeshVarIndLax(mesh, str);
      if(vi<0) continue;

      /* we do 2doutputall */
      vi0 = MeshVarIndComponent0(mesh, vi);
      for(i=0; i<MeshVarNComponents(mesh, vi0); i++)
        gnuplot_out2d_meshvar(mesh, VarName(vi0+i), Iter, Time);
    }
  }

  /* 3d output */
  d = 3;
  if(TimeForMeshOutput_di_dt(mesh, di[d], dt[d]))
  {
    errorexit("3d output not implemented");
  }
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
  //int Noutput = NextAll + 4; // really just need 4 from 0d-,1d-,2d-,3d-output 
  int Noutput = 4;
  int di[Noutput];
  double dt[Noutput];
  char output[Noutput][64];
  char *name, s[64];
  int d, n;

  errorexit("this function is not tested yet!");

  name = VarName(vindex);

  /* read all pars for all dims */
  for(d = 0; d <= 3; d++)
  {
    sprintf(s, "%ddoutiter", d);
    di[d] = Geti(Par(s));
    sprintf(s, "%ddouttime", d);
    dt[d] = Getd(Par(s));
    sprintf(output[d], "%ddoutput", d);
  }

  /* check if we "name" is contained in any output par */
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

  errorexit("this function is not tested yet!");

  return 
  TimeForNodeOutput_di_dt(node,Geti(Par("0doutiter")),Getd(Par("0douttime"))) ||
  TimeForNodeOutput_di_dt(node,Geti(Par("1doutiter")),Getd(Par("1douttime"))) ||
  TimeForNodeOutput_di_dt(node,Geti(Par("2doutiter")),Getd(Par("2douttime"))) ||
  TimeForNodeOutput_di_dt(node,Geti(Par("3doutiter")),Getd(Par("3douttime")));
}
