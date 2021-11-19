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
/* stuff to make it compile without nmesh */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define tVarList int
#define tNode int
#define MyLnode 0
#define formylnodes(mesh)      for(int node; node<1; node++)
#define forpoints(node, ijk)   for(ijk=0; ijk<10; ijk++)
double Mem[100]; // instead of the data arrays of nmesh


//#include "nmesh.h"

#define POW2(x)      (x)*(x)
#define POWER(x, y)  pow(x,y)



void exmpl_func(tVarList *vlu)
{
  double e = 0;// Getd(Par("FOCCZ4_e"));
  int ijk; /* grid point index */
''',
''':Decl: W[i,j,k,l]; T[i,k,l]; S[i,j]; A[i,j]; SS[i,j,k]; AA[i,j,k];
          R[i,j]; :
          Access={VAR}{COMP}[ijk] : Format =
  double *{VAR}{COMP} = Mem; // Vard( node, Vind(vlu,{LI}) );
''',
''':Decl: cA; cS; : Access={VAR}{COMP}: Format =
  double {VAR}{COMP} = 1.; // name={VAR} ID={VARID} compind={CI} listind={LI}
''',
''':Decl: AUTOVARS :
          Access={VAR}{COMP} : DeclFunc=make_DeclList : Format =
  double {VAR}{COMP}; // name={VAR} ID={VARID} compind={CI} listind={LI}
''',

r''':Text =

  formylnodes(mesh)
  {
    tNode *node = MyLnode;
    int ijk;

    forpoints(node, ijk)
    {

''',

    'R2a[i,j] = cS*W[i,j,k,l]*S[k,l] + cA*W[i,j,k,l]*A[k,l]',
    'R2b[i,j] = cS*W[k,j,i,l]*S[k,l] + cA*W[k,j,i,l]*A[k,l]',
    'R[i,j]   = R2a[i,j] + R2b[i,j]',

    'R1a[i] = cS*W[i,j,k,l]*SS[j,k,l] + cA*W[i,j,k,l]*AA[j,k,l]',

    'R2a[i,j] = T[i,k,l]*SS[j,k,l] + T[i,k,l]*AA[j,k,l]',
    'R2b[i,j] = T[k,l,i]*SS[j,k,l] + T[k,l,i]*AA[j,k,l]',

    'R[i,j] = R2a[i,j] + R2b[i,j]',

r''':Text =
    }
  }
}


/* add main so that it compiles without nmesh  */
int main()
{
  int i = 2;
  Mem[0] = 42.;
  exmpl_func(&i);
  printf("Mem[0]=%g\n", Mem[0]);
}
''',
)

########################################################################
# get all Eqn components and the undeclared vars that only appear on LHS
########################################################################
Eqs = read_LHS_RHS_strlist(tocompute)
allEqs, AUTOVARS = assemble_all_EqnComponents(Eqs)

########################################################################
# declare all symmetries for all variables we use
########################################################################
symmetries = {
    S[i,j] :   {  '+': ( [j,i], ) },
    A[i,j] :   {  '-': ( [j,i], ) },
    T[k,l,m] : {   '+': ( [l,k,m], [k,m,l], ) },
    SS[k,l,m] : {  '+': ( [l,k,m], [k,m,l], ) },
    AA[k,l,m] : {  '-': ( [l,k,m], [k,m,l], ) },
##    T[i,j,k,l,m] : {  '+': ( [j,i,k,l,m], [i,j,l,k,m] ),
##                      '-': ( [i,j,k,m,l], ) }
}


########################################################################
# apply symmetries and remove duplicates. This step takes the longest.
########################################################################
allEqs = apply_symmetries_to_all_EqnComponents(symmetries, Eqs, allEqs,
                                               AUTOVARS)

########################################################################
# run more sympy simplification operations on RHSs (this is optional)
########################################################################
allEqs = simplify_all_EqnComponents(sympy.simplify, allEqs)
#allEqs = simplify_all_EqnComponents(sympy.expand, allEqs)
allEqs = simplify_all_EqnComponents(sympy.N, allEqs)

########################################################################
# now write all into a .c file
########################################################################
cfilename = __file__.replace('.py', '.c')
write_Eqs(cfilename, allEqs, AUTOVARS)
