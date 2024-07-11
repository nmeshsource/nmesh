########################################################################
# Initilize code generator:
########################################################################
# get file SymPyToC.py directly into here (like when we include in C)
codegen = '../../Math/SymPyToC/SymPyToC.py'
exec(compile(open(codegen).read(), codegen, 'exec'))

########################################################################
# indices we need, and their min and max value in brackets
########################################################################
i = sympy.Idx('i', (1, 3))
j = sympy.Idx('j', (1, 3))
k = sympy.Idx('k', (1, 3))
l = sympy.Idx('l', (1, 3))
m = sympy.Idx('m', (1, 3))
n = sympy.Idx('n', (1, 3))

########################################################################
# Program Text and Eqations (in one long string tuple)
########################################################################
tocompute = (
r''':Text =
/* Some defines to make it compile without nmesh: */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define tVarList int
#define tElm int
#define MyElm 0
#define formyelms(mesh)      for(int elm; elm<1; elm++)
#define NPTS 3
#define forpoints(elm, ijk)  for(ijk=0; ijk<NPTS; ijk++)
#define Vard(elm, vi) Mem+NPTS*(vi)
#define Vind(vl, vli)  vl[vli]
double Mem[1000]; // instead of the data arrays of nmesh

/* To compile with nmesh remove the above #includes and #defines
   INSTEAD uncomment this: */
//#include "nmesh.h"

#define POWER(x,y) pow(x,y)
#define POW2(x)    ((x)*(x))
#define POW3(x)    ((x)*(x)*(x))
#define SQRT(x)    sqrt(x)
#define CBRT(x)    cbrt(x)


/**************************************************/
/* function to compute v_i = M_ij u_j on one elm */
/**************************************************/
void MatrixMul(tElm *elm, tVarList *vlv, int iM, tVarList *vlu)
{
  double e = 0;// Getd(Par("FOCCZ4_e"));
  int ijk; /* grid point index */
''',
''':Decl: u[i]; :
          Access={VAR}{COMP}[ijk] : Format =
  double *{VAR}{COMP} = Vard( elm, Vind(vlu,{LI}) );
''',
''':Decl: v[i]; :
          Access={VAR}{COMP}[ijk] : Format =
  double *{VAR}{COMP} = Vard( elm, Vind(vlv,{LI}) );
''',
''':Decl: M[i,j]; :
          Access={VAR}{COMP}[ijk] : Format =
  double *{VAR}{COMP} = Vard( elm, iM+{LI} ); // name={VAR} ID={VARID} compind={CI} listind={LI}
''',
''':Decl: AUTOVARS :
          Access={VAR}{COMP} : DeclFunc=make_DeclList : Format =
  double {VAR}{COMP}; // name={VAR} ID={VARID} compind={CI} listind={LI}
''',

r''':Text =

  /* To compile with nmesh remove these printfs */
  printf("v1 -Mem=%zu\n", v1-Mem);
  printf("M11-Mem=%zu\n", M11-Mem);
  printf("M12-Mem=%zu\n", M12-Mem);
  printf("u1 -Mem=%zu\n", u1-Mem);
  printf("\n");
  printf("M12[0]=%g\n", M12[0]);

  forpoints(elm, ijk)
  {
''',

    # actual matrix multiplication
    'v[i] = M[i,j] * u[j]',

r''':Text =
  }
}



/* To compile with nmesh remove the entire main below: */

//test main so that it compiles and runs without nmesh
int main()
{
  int *elm = NULL;
  tVarList vlu[] = { 0,1,2 };
  tVarList vlv[] = { 3,4,5 };
  int iM = 6;

  Mem[NPTS*0 + 0] = 2.; //vlu[0] at point0
  Mem[NPTS*1 + 0] = 1.; //vlu[1] at point0
  Mem[NPTS*6 + 0] = 4.; //M11 at point0
  Mem[NPTS*7 + 0] = 5.; //M12 at point0
  
  printf("vlu[1] at point0 = Mem[NPTS*1 + 0] = %g\n", Mem[NPTS*1 + 0]);
  MatrixMul(elm, vlv, iM, vlu);
  printf("\n");
  printf("vlv[0] at point0 = Mem[NPTS*3 + 0] = %g\n", Mem[NPTS*3 + 0]);
}
''',
)



if __name__ == '__main__':
  ########################################################################
  # Create pool of processes. Can use multiprocessing or multiprocess
  ########################################################################
  import multiprocessing as mp
  #mp.set_start_method('fork')
  #mp.set_start_method('spawn')
  #mp.set_start_method('forkserver')
  #import multiprocess as mp
  SymPyToC_pool = mp.Pool()

  ########################################################################
  # get all Eqn components and the undeclared vars that only appear on LHS
  ########################################################################
  Eqs = read_LHS_RHS_strlist(tocompute)
  allEqs, AUTOVARS = assemble_all_EqnComponents(Eqs)

  ########################################################################
  # declare all symmetries for all variables we use
  ########################################################################
  symmetries = {
      M[i,j] :   {  '+': ( [j,i], ) }, # make M be symmetric
  }

  ########################################################################
  # apply symmetries and remove duplicates. This step takes the longest.
  ########################################################################
  allEqs = apply_symmetries_to_all_EqnComponents(symmetries, Eqs, allEqs,
                                                 AUTOVARS)

  ########################################################################
  # run more sympy simplification operations on RHSs (this is optional)
  ########################################################################
  #allEqs = simplify_all_EqnComponents(sympy.simplify, allEqs)
  #allEqs = simplify_all_EqnComponents(sympy.expand, allEqs)
  allEqs = simplify_all_EqnComponents(sympy.N, allEqs)

  ########################################################################
  # try to find and remove unused components (this is optional)
  ########################################################################
  allEqs = remove_unused_Components(allEqs)

  ########################################################################
  # now write all into a .c file
  ########################################################################
  cfilename = __file__.replace('.py', '.c')
  write_Eqs(cfilename, allEqs, AUTOVARS)

  ########################################################################
  # close pool of processes
  ########################################################################
  if SymPyToC_pool != None:
    SymPyToC_pool.close()
