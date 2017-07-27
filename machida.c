#include "copyright.h"
/*============================================================================*/
/*! \file larson.c
 *  \brief Problem generator for Machida Outflow Generation Model(Masahiro N Machida (ApJ 796(2014))
 *  config  = --with-problem=machida --with-eos=isothermal --with-gas=hydro --with-gravity=fft --enable-fft
/*============================================================================*/


#include <math.h>
#include <stdio.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

 Real omega0 = 3.469e-2; //Rotation
static Real vphi(const Real x1,const Real x2) {
  return sqrt(x1*x1+x2*x2)*omega0;
}


void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,no=0;
  Real x1,x2,x3,rad,drat,cs,diag,temp,G;
  #ifdef MHD
  Real B0;
  #endif /* MHD */
  /* Ambient Density,Pressure,velocity components and magnetic field */
  
  Real Rc,omega0,d0,dmedium,Rmax,u0,v0,w0;  
  
  /* Spherical Density,Pressure,Velocity components and magnetic field */
  
  //Real dr,pr,ur,vr,wr,br;//
  
   
  rad  = par_getd("problem","radius"); //Cloud Radius in terms of Rc
  
  /* Set paramters in ambient medium */
  
 Rc = 0.2;  //Critical Ebert Radius
 omega0 = 3.469e-2; //Rotation
 d0    = 3.92*1.5;     //Central Density
 dmedium = 0.04; //Density of ISM
 Rmax = Rc*rad ;   //Cloud Radius
 #ifdef MHD
 B0  = 0.123 ;     // Magnetic Field Strength in Z direction
 #endif
 G   = 3.941e-3 ;
 #ifndef ISOTHERMAL
  cs=sqrt(Gamma*p0/d0);
 #else
    cs=Iso_csound;
 #endif
  /* Initialize the grid */

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        diag = sqrt(x1*x1 + x2*x2 + x3*x3);
        temp = sqrt(x1*x1+x2*x2);
        
        u0= -1.0*vphi(x1,x2)*(x2/temp);
        v0= -1.0*vphi(x1,x2)*(x1/temp);
        if(diag<=Rc) {

         printf("%f, %f,%f \n",x1,x2,x3); 
         pGrid->U[k][j][i].d= d0;
         
         pGrid->U[k][j][i].M1= d0*u0;
         pGrid->U[k][j][i].M2= d0*v0;
         pGrid->U[k][j][i].M3= d0*w0;
         
         
         #ifdef MHD
          
          pGrid->B3i[k][j][i] = B0;
          
          pGrid->U[k][j][i].B3c = B0;
         #endif
         
        #ifdef ADIABATIC
          pGrid->U[k][j][i].E = p0/Gamma_1
        #ifdef MHD
           pGrid->U[k][j][i].E+ = 0.5*(B0*B0);
        #endif
            pGrid->U[k][j][i].E+ = 0.5*d0*(u0*u0+v0*v0);
        #endif 
        }
        
      else if(diag>Rc && diag<=Rmax) {
         pGrid->U[k][j][i].d=  d0*Rc*Rc/pow(diag,2);
         pGrid->U[k][j][i].M1= d0*Rc*Rc/pow(diag,2)*u0;
         pGrid->U[k][j][i].M2= d0*Rc*Rc/pow(diag,2)*v0;
         pGrid->U[k][j][i].M3= d0*Rc*Rc/pow(diag,2)*w0;
        
        
          #ifdef MHD
          
          pGrid->B3i[k][j][i] = B0;
          
          pGrid->U[k][j][i].B3c = B0;
         #endif
         
        #ifdef ADIABATIC
          pGrid->U[k][j][i].E = p0/Gamma_1
        #ifdef MHD
           pGrid->U[k][j][i].E+ = 0.5*(B0*B0);
        #endif
            pGrid->U[k][j][i].E+ = 0.5*(cs**2)/(2*PI*G*pow(diag,2))*(u0*u0+v0*v0);
        #endif 
         }

       else {
          pGrid->U[k][j][i].d= dmedium;
         
         pGrid->U[k][j][i].M1= 0.0;
         pGrid->U[k][j][i].M2= 0.0;
         pGrid->U[k][j][i].M3= 0.0;
  
         #ifdef MHD
          
          pGrid->B3i[k][j][i] = B0;
          
          pGrid->U[k][j][i].B3c = B0;
         #endif
         
        #ifdef ADIABATIC
          pGrid->U[k][j][i].E = p0/Gamma_1;
          
        #ifdef MHD
           pGrid->U[k][j][i].E += 0.5*(B0*B0);
        #endif
        #endif
        }
         }
        }
       }
  #ifdef SELF_GRAVITY
/* Set gravity constant*/
  four_pi_G = (4*PI)*3.941e-3 ;
  grav_mean_rho = d0;

  #endif /* SELF_GRAVITY */ 
  
  return;
  
  }
  
  /*=============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 
 *----------------------------------------------------------------------------*/
 
 void problem_write_restart(MeshS *pM, FILE *fp)
{ static Real d0    = 3.92*1.5;
  four_pi_G = (4*PI)*3.941e-3 ;
  grav_mean_rho = d0;

  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{ static Real d0    = 3.92*1.5;
  four_pi_G = (4*PI)*3.941e-3 ;
  grav_mean_rho = d0;

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}    
