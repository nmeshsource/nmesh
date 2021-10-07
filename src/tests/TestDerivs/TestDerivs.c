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
    double *uyx= Vard(node,Ind("TestDerivs_Err_dduxx")+3);
    double *uyy= Vard(node,Ind("TestDerivs_Err_dduxx")+4);
    double *uyz= Vard(node,Ind("TestDerivs_Err_dduxx")+5);
    double *uzx= Vard(node,Ind("TestDerivs_Err_dduxx")+6);
    double *uzy= Vard(node,Ind("TestDerivs_Err_dduxx")+7);
    double *uzz= Vard(node,Ind("TestDerivs_Err_dduxx")+8);

    /* compute the derivs */
    if(mesh->iteration % 2 == 0)
    {
      /* use ddu_ij = d_j (d_i u)  when mesh->iteration is even */
      cart_partials_dTensor_di(node, Ind("TestDerivs_u"),
                                     Ind("TestDerivs_Err_dux"), NULL);
      cart_partials_dTensor_di(node, Ind("TestDerivs_Err_dux"),
                                     Ind("TestDerivs_Err_dduxx"), NULL);
    }
    else
    {
      /* use ddu_ij = d_i (d_j u)  when mesh->iteration is odd */
      cart_partials_diTensor(node, Ind("TestDerivs_u"),
                                   Ind("TestDerivs_Err_dux"), NULL);
      cart_partials_diTensor(node, Ind("TestDerivs_Err_dux"),
                                   Ind("TestDerivs_Err_dduxx"), NULL);
    }

    /* subtract true values */
    forpoints(node,i)
    {
      double x = px[i];
      double y = py[i];
      double z = pz[i];
      double U, Ux,Uy,Uz, Uxx,Uxy,Uxz, Uyx,Uyy,Uyz, Uzx, Uzy,Uzz;

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
      /* analytic derivs */
      U = u[i];
      Ux  = A*(-x + x0)/Power(sigmax,2) * U;
      Uy  = A*(-y + y0)/Power(sigmay,2) * U;
      Uz  = A*(-z + z0)/Power(sigmaz,2) * U;
      Uxx = A*(-Power(sigmax,2) + Power(x - x0,2))/Power(sigmax,4) * U;
      Uxy = A*((x - x0)*(y - y0))/(Power(sigmax,2)*Power(sigmay,2)) * U;
      Uxz = A*((x - x0)*(z - z0))/(Power(sigmax,2)*Power(sigmaz,2)) * U;
      Uyx = Uxy;
      Uyy = A*(-Power(sigmay,2) + Power(y - y0,2))/Power(sigmay,4) * U;
      Uyz = A*((y - y0)*(z - z0))/(Power(sigmay,2)*Power(sigmaz,2)) * U;
      Uzx = Uxz;
      Uzy = Uyz;
      Uzz = A*(-Power(sigmaz,2) + Power(z - z0,2))/Power(sigmaz,4) * U;

      /* subtract analytic derivs from numerical ones */
      ux[i] -= Ux;
      uy[i] -= Uy;
      uz[i] -= Uz;
      uxx[i]-= Uxx;
      uxy[i]-= Uxy;
      uxz[i]-= Uxz;
      uyx[i]-= Uyx;
      uyy[i]-= Uyy;
      uyz[i]-= Uyz;
      uzx[i]-= Uzx;
      uzy[i]-= Uzy;
      uzz[i]-= Uzz;
    }
  }
  return 0;
}
