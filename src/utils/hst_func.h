#ifndef HST_FUNC_H
#define HST_FUNC_H


#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "../athena.h"
#include "../ath_array.c"
#include "code_units.h"

// need defined and finished
//#include "../eos/eos.h"
//#include "../field/field.h"
//#include "../hydro/hydro.h"




  return cold_gas_mass;

double cold_gas(MeshBlock *pmb, int out){
  double cold_gad_mass = 0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmg->ks, ke=pmb->ke;
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
            
        double rho = pmb->phydro->u(IDN,k,j,i);
        double prs = pmb->phydro->w(IPR,k,j,i);

        double temp = (prs / rho) * KELVIN * mu ;

        if (temp <= T_cold){
          cold_gas_mass += rho*pmb->pcoord->GetCellVolume(k,j,i);
        }   

      }   
    }   
  }   
  return cold_gas_mass;
}

double rho_sum(MeshBlock *pmb, int iout){

  double rho_sum    = 0;

  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        
        double rho = pmb->phydro->u(IDN,k,j,i);

        rho_sum += rho;

      }
    }
  }

  return rho_sum;
}

double rho_sq_sum(MeshBlock *pmb, int iout){

  double rho_sq_sum = 0;

  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        
        double rho = pmb->phydro->u(IDN,k,j,i);
        rho_sq_sum += rho*rho;

      }
    }
  }

  return rho_sq_sum;
}

double c_s_sum(MeshBlock *pmb, int iout){

  double c_s_sum = 0;
  double gamma = pmb->peos->GetGamma();

  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        
        double rho = pmb->phydro->u(IDN,k,j,i);
        double prs = pmb->phydro->w(IPR,k,j,i);

        c_s_sum += sqrt(gamma*prs/rho);

      }
    }
  }

  return c_s_sum;
}

double Pth_sum(MeshBlock *pmb, int iout){

  double Pth_sum = 0;
  

  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;


  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {

        Pth_sum += pmb->phydro->w(IPR,k,j,i);
        
      }
    }
  }
    
  

  return Pth_sum;

}

double PB_sum(MeshBlock *pmb, int iout){

  double PB_sum = 0;
  
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;  

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        PB_sum += pmb->pfield->b.x1f(k,j,i)*pmb->pfield->b.x1f(k,j,i);
        PB_sum += pmb->pfield->b.x2f(k,j,i)*pmb->pfield->b.x2f(k,j,i);
        PB_sum += pmb->pfield->b.x3f(k,j,i)*pmb->pfield->b.x3f(k,j,i);
      }
    }
  }

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      PB_sum += pmb->pfield->b.x1f(k,j,ie+1)*pmb->pfield->b.x1f(k,j,ie+1);
    }
  }

  for (int k=ks; k<=ke; ++k) {
    for (int i=is; i<=ie; ++i) {
      PB_sum += pmb->pfield->b.x2f(k,je+1,i)*pmb->pfield->b.x2f(k,je+1,i);
    }
  }

  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      PB_sum += pmb->pfield->b.x3f(ke+1,j,i)*pmb->pfield->b.x3f(ke+1,j,i);
    }
  }

  return PB_sum/2.0;

}

double Bx_sum(MeshBlock *pmb, int iout){

  double Bx_sum = 0;
  
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;  

  for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          Bx_sum += pmb->pfield->b.x1f(k,j,i);
        }
      }
    }

  return Bx_sum;

}

double By_sum(MeshBlock *pmb, int iout){

  double By_sum = 0;
  
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          By_sum += pmb->pfield->b.x2f(k,j,i);
        }
      }
    }

  return By_sum;

}

double Bz_sum(MeshBlock *pmb, int iout){

  double Bz_sum = 0;
  
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          Bz_sum += pmb->pfield->b.x3f(k,j,i);
        }
      }
    }

  return Bz_sum;

}


#endif
