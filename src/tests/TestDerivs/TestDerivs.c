/* TestDerivs.c */
/* Wolfgang Tichy 2005 */


#include "nmesh.h"
#include "TestDerivs.h"

#define Power pow

/* initialize TestDerivs */
int TestDerivs_startup(tMesh *mesh)
{
  double A         = Getd(Par("TestDerivs_A"));
  double sigmax    = Getd(Par("TestDerivs_sigmax"));
  double sigmay    = Getd(Par("TestDerivs_sigmay"));
  double sigmaz    = Getd(Par("TestDerivs_sigmaz"));
  double x0        = Getd(Par("TestDerivs_x0"));
  double y0        = Getd(Par("TestDerivs_y0"));
  double z0        = Getd(Par("TestDerivs_z0"));
  
  printf("Initializing TestDerivs:\n");

  /* enable all vars */
  coordinates_init(mesh);
  enablevar(mesh, Ind("TestDerivs_u"));
  enablevar(mesh, Ind("TestDerivs_Err_dux"));
  enablevar(mesh, Ind("TestDerivs_Err_dduxx"));

  /* set initial data in nodes */
  formylnodes(mesh)
  {  
    tNode *node = MyLnode;
    int i;
    //double *pX = Vard(node,Ind("X"));
    //double *pY = Vard(node,Ind("Y"));
    //double *pZ = Vard(node,Ind("Z"));
    double *px = Vard(node,Ind("x"));
    double *py = Vard(node,Ind("y"));
    double *pz = Vard(node,Ind("z"));
    double *u = Vard(node,Ind("TestDerivs_u"));

    forpoints(node,i)
    {
      double x = px[i];
      double y = py[i];
      double z = pz[i];

      /* Gaussian */
      u[i] = A*exp(-0.5*( 
               (x-x0)*(x-x0)/(sigmax*sigmax) + (y-y0)*(y-y0)/(sigmay*sigmay)
             + (z-z0)*(z-z0)/(sigmaz*sigmaz)   ) );
    }
  }

  return 0;
}


/* compute absolute error in ANALYSIS */
int TestDerivs_analyze(tMesh *mesh)
{
  double A         = Getd(Par("TestDerivs_A"));
  double sigmax    = Getd(Par("TestDerivs_sigmax"));
  double sigmay    = Getd(Par("TestDerivs_sigmay"));
  double sigmaz    = Getd(Par("TestDerivs_sigmaz"));
  double x0        = Getd(Par("TestDerivs_x0"));
  double y0        = Getd(Par("TestDerivs_y0"));
  double z0        = Getd(Par("TestDerivs_z0"));
  
  printf("TestDerivs: computing absolute error in derivatives\n");
  
  /* set initial data in nodees */
  formylnodes(mesh)
  {  
    tNode *node = MyLnode;
    int i;
    //double *pX = Vard(node,Ind("X"));
    //double *pY = Vard(node,Ind("Y"));
    //double *pZ = Vard(node,Ind("Z"));
    double *px = Vard(node,Ind("x"));
    double *py = Vard(node,Ind("y"));
    double *pz = Vard(node,Ind("z"));
    double *u  = Vard(node,Ind("TestDerivs_u"));
    double *ux = Vard(node,Ind("TestDerivs_Err_dux"));
    double *uy = Vard(node,Ind("TestDerivs_Err_duy"));
    double *uz = Vard(node,Ind("TestDerivs_Err_duz"));
    double *uxx= Vard(node,Ind("TestDerivs_Err_dduxx"));
    double *uxy= Vard(node,Ind("TestDerivs_Err_dduxx")+1);
    double *uxz= Vard(node,Ind("TestDerivs_Err_dduxx")+2);
    double *uyy= Vard(node,Ind("TestDerivs_Err_dduxx")+3);
    double *uyz= Vard(node,Ind("TestDerivs_Err_dduxx")+4);
    double *uzz= Vard(node,Ind("TestDerivs_Err_dduxx")+5);

    /* compute the derivs */
    cart_partials(node, Ind("TestDerivs_u"), Ind("TestDerivs_Err_dux"));
    cart_partials(node, Ind("TestDerivs_Err_duz"), Ind("TestDerivs_Err_dduyy"));
    cart_partials(node, Ind("TestDerivs_Err_duy"), Ind("TestDerivs_Err_dduxz"));
    cart_partials(node, Ind("TestDerivs_Err_dux"), Ind("TestDerivs_Err_dduxx"));

    /* subtract true values */
    forpoints(node,i)
    {
      double x = px[i];
      double y = py[i];
      double z = pz[i];

      if(fabs(x)+fabs(y)+fabs(z)>1e299) continue; /* give up if x,y,z is inf */
      /*
      u[i] = A*exp(-0.5*( 
               (x-x0)*(x-x0)/(sigmax*sigmax) + (y-y0)*(y-y0)/(sigmay*sigmay)
             + (z-z0)*(z-z0)/(sigmaz*sigmaz)   ) );

      $Assumptions = sigmax>0 && sigmay>0 && sigmaz>0
      u=A*Exp[-( (x-x0)*(x-x0)/(sigmax*sigmax) + (y-y0)*(y-y0)/(sigmay*sigmay)
                +(z-z0)*(z-z0)/(sigmaz*sigmaz))/2]

      CForm[Simplify[D[u,x]/u]]
       = A*(-x + x0)/Power(sigmax,2)
      CForm[Simplify[D[u,y]/u]]
       = A*(-y + y0)/Power(sigmay,2)
      CForm[Simplify[D[u,z]/u]]
       = A*(-z + z0)/Power(sigmaz,2)
       
      CForm[Simplify[D[u,x,x]/u]]
       = A*(-Power(sigmax,2) + Power(x - x0,2))/Power(sigmax,4)
      CForm[Simplify[D[u,x,y]/u]]
       = A*((x - x0)*(y - y0))/(Power(sigmax,2)*Power(sigmay,2))
      CForm[Simplify[D[u,x,z]/u]]
       = A*((x - x0)*(z - z0))/(Power(sigmax,2)*Power(sigmaz,2))
      CForm[Simplify[D[u,y,y]/u]]
       = A*(-Power(sigmay,2) + Power(y - y0,2))/Power(sigmay,4)
      CForm[Simplify[D[u,y,z]/u]]
       = A*((y - y0)*(z - z0))/(Power(sigmay,2)*Power(sigmaz,2))
      CForm[Simplify[D[u,z,z]/u]]
       = A*(-Power(sigmaz,2) + Power(z - z0,2))/Power(sigmaz,4)
      */
      //printf("z0=%f ",z0);
      ux[i] -= A*(-x + x0)/Power(sigmax,2) * u[i];
      uy[i] -= A*(-y + y0)/Power(sigmay,2) * u[i];
      uz[i] -= A*(-z + z0)/Power(sigmaz,2) * u[i];
      uxx[i]-= A*(-Power(sigmax,2) + Power(x - x0,2))/Power(sigmax,4) * u[i];
      uxy[i]-= A*((x - x0)*(y - y0))/(Power(sigmax,2)*Power(sigmay,2)) * u[i];
      uxz[i]-= A*((x - x0)*(z - z0))/(Power(sigmax,2)*Power(sigmaz,2)) * u[i];
      uyy[i]-= A*(-Power(sigmay,2) + Power(y - y0,2))/Power(sigmay,4) * u[i];
      uyz[i]-= A*((y - y0)*(z - z0))/(Power(sigmay,2)*Power(sigmaz,2)) * u[i];
      uzz[i]-= A*(-Power(sigmaz,2) + Power(z - z0,2))/Power(sigmaz,4) * u[i];
    }
  }
  return 0;
}
