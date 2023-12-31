#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "defs.h"
#include <stdbool.h>
#include <mpi.h>
//#include "./utils/hst_func.hpp"
#include "./utils/code_units.h"
#include "./utils/tcm3.h"

//algorithm and memory and omp.h headers

Real tfloor;
Real tnotcool;
Real tcut_hst;
Real r_drop;
Real Lambda_fac;
Real Lambda_fac_time;
Real total_cooling = 0.0;

typedef struct T{
    int value;
} T;

//T *make_unique_T(int value) {
  //  T *ptr = (T *)malloc(sizeof(T));
   // if (ptr) {
     //   ptr->value = value;
   // }
    //return ptr;
//}


//T *obj = make_unique_T(T);
  //  if (obj) {
    //    printf("Value: %d\n", obj->value);
      //  free(obj);
   // } else {
     //   printf("Memory allocation failed\n");
    //}

    //return 0;
//}

//typedef struct cooler{
  //  int value;
//} cooler;

//cooler *make_unique_cooler(int value) {
  //  cooling *ptr = (cooler *)malloc(sizeof(cooler));
    //if (ptr) {
      //  ptr->value = value;
//    }
  //  return ptr;
//}


//cooler *obj = make_unique_cooler(cooler);
  //  if (obj) {
    //    printf("Value: %d\n", obj->value);
      //  free(obj);
  //  } else {
    //    printf("Memory allocation failed\n");
    //}

    //return 0;
//}

Real L1 = 0.0;
Real L2 = 0.0;
Real L3 = 0.0;

int nx1 = 1;
int nx2 = 1;
int nx3 = 1;

int cooling_flag = 1;
int shift_flag  = 1;

Real x3min = 0.0;
Real x3max = 0.0;
Real amb_rho = 1.0;
Real front_thick = 2.5;
Real v_shear = 100 * (1.023*1e-3);
Real v_shift0 = 0.0;
Real front_velocity = 0.0;
Real knx_KH = 1.0;
Real kny_KH = 1.0;
Real amp_KH = 0.01;
Real B_x = 0.0;
Real B_y = 0.0;
Real B_z = 1.0;
Real front_posn_old = 0.0;
_Bool cooling_flag_print_count = false;
_Bool shift_flag_print_count = false;
_Bool DEBUG_FLAG_MIX = false;
Real time_last = -1.0;
Real shift_start = 0.0;
Real debug_hst = 0.0;
//static Real T_floor;
//static Real T_cold;
//static Real T_ceil;
//static Real T_hot;
//static Real T_cut_mul;
// new PARAMETERS
//static Real Kbolts = 1.38e-16
Real velx,vely,velz;
Real rho;
Real g;
//static Real mu;
//static Real muH;
//static Real mue;
//static Real Kboltz;


void read_input(DomainS *pDomain){
  L1 = par_getd("domain1", "x1max")-par_getd("domain1", "x1min");
  L2 = par_getd("domain1", "x2max")-par_getd("domain1", "x2min");
  L3 = par_getd("domain1", "x3max")-par_getd("domain1", "x3min");

  x3min = par_getd("domain1", "x3min");
  x3max = par_getd("domain1", "x3max");

  nx1 = par_getd("domain1", "Nx1");
  nx2 = par_getd("domain1", "Nx2");  //domain2?
  nx3 = par_getd("domain1", "Nx3");  //domain3?

  cooling_flag = par_geti("problem","cooling_flag");
  shift_flag = par_geti("problem", "shift_flag");

  shift_start = par_getd("problem", "shift_start");

  //amb_rho, front_thick, v_shift, Thot, Tcold, Tceil, Tfloor, T_cut_mul, Xsol, Zsol, Lambda_fac,
  v_shift0 = par_getd("problem", "v_shift");
  amb_rho = par_getd("problem", "amb_rho");
  v_shear = par_getd("problem", "v_shear"); 
  DEBUG_FLAG_MIX = par_geti("problem", "DEBUG_FLAG");
  
  //bfield = par_getd("problem", "b0");
#ifdef MHD
  B_x = par_getd("problem", "B_x");
  B_y = par_getd("problem", "B_y");  
  B_z = par_getd("problem", "B_z");
#endif
  Xsol = par_getd("problem", "Xsol");
  Zsol = par_getd("problem", "Zsol");
  X = Xsol * 0.7381;
  Z = Zsol * 0.0134;
  Y = 1-X-Z;

  mu = 1.0/(2.*X+3.*(1.-X-Z)/4.+Z/2.);
  mue = 2.0/(1.0+X);
  muH = 1.0/X;

  T_floor = par_getd("problem", "T_floor");
  T_ceil = par_getd("problem", "T_ceil");
  T_hot = par_getd("problem", "T_hot");
  T_cold = par_getd("problem", "T_cold");
  T_cut_mul = par_getd("problem", "T_cut_mul");

  g = par_getd("problem", "gamma");
  printf("input file read\n");

  return;
}

void townsend_cooling(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0, j=0, k=0;
  int is,ie,js,je,ks,ke,iprob;
  Real x1,x2,x3;
  Real time, dt;
  printf("townsend cooling initiated \n");
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  //iprob = par_geti("problem","iprob");
  //if (iprob ==2) {
  Real a=0.05;
  Real sigma=0.2;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        debug_hst += 1.0;
////////////// UNSURE IF NEEDED HERE
        pGrid->U[k][j][i].d = rho;
        pGrid->U[k][j][i].M1 = v_shear*tanh(x2/a);
	pGrid->U[k][j][i].M2 = amb_rho*sin(2.0*PI*x1)*exp(-(x2*x2)/(sigma*sigma));
        pGrid->U[k][j][i].M3 = 0.0;

	Real temp = (g-1.0)*rho*mu/(Kboltz*mue*muH);
	  // Real temp = pres/rho *Kboltz*mu;

	if (isnan(temp)){
	  printf("temp \n",temp);
	  printf("rho \n", rho);
	    
	    //Real E_kin = (Kbolts * temp * rho) / ((gamma-1.0)*mp*mh);
	  Real E_kin = 0.5 * (SQR(pGrid->U[k][j][i].M1)+SQR(pGrid->U[k][j][i].M2)+SQR(pGrid->U[k][j][i].M3))* pGrid->U[k][j][i].d;
	  Real E_mag = 0.0;
	    
#ifdef MHD
	  E_mag = 0.5 *B_x*B_x + B_y*B_y +B_z*B_z;
#endif
	  pGrid->U[k][j][i].E = E_kin + E_mag + (T_floor/(Kboltz*mu))*pGrid->U[k][j][i].d/(g-1.0);
        }
        if (temp > T_floor) {
	  if (temp<T_cut){
	    Real rho = rho;
	    Real temp_cgs = temp;
	    Real rho_cgs = rho * unit_density; // NEED CODE_UNITS
            Real dt_cgs = dt *unit_time;
	    Real cLfac = Lambda_fac;
//////////////// DEFINE TOWNSEND!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	    Real temp_new = fmax(Cooling_townsend(temp_cgs,rho_cgs,dt_cgs,cLfac), T_floor);
	    if (isnan(temp_new)) {
	      printf("temp \n", temp_cgs);
		
	      Real ccool = ((T_floor-temp)/(Kboltz*mu))*pGrid->U[k][j][i].d/(g-1.0);
	      pGrid->U[k][j][i].E += ccool;
	    }
	    else{
	      Real ccool = ((temp_new-temp)/(Kboltz*mu))*pGrid->U[k][j][i].d/(g-1.0);
   	      pGrid->U[k][j][i].E += ccool;
              printf("temp < T_floor\n");
	      if (DEBUG_FLAG_MIX) {
		total_cooling -=1.0;
	      }
	      else{
		total_cooling -= ccool;
	      }
	    }
          }
        }
        else{
	  Real ccool = ((T_floor-temp)/(Kboltz*mu))* pGrid->U[k][j][i].d/(g-1.0);
	  pGrid->U[k][j][i].E += ccool;
	}
      }
    }
  }
}


void frame_shift(DomainS *pDomain) {
  GridS *pGrid = pDomain->Grid;
  Real time, dt, g;
  printf("frame shift initiated\n");
  
  if (time_last==time){
    return;
  }
  else{
    time_last = time;
  }

  Real local_cold_mass = 0.0;
  Real global_cold_mass = 0.0;
  
  Real global_v3_sum = 0.0;
  Real local_v3_sum = 0.0;
  
  double front_posn_new = 0.0;
  Real v_shift_t_new = 0.0;
  
  Real sim_dt = dt;
  
  g = par_getd("problem", "gamma");
  Real c_s = sqrt(g * T_floor/(mu*Kboltz));
  Real c_s_cap = 100.0;
  int i=0, j=0, k=0;
  int is,ie,js,je,ks,ke,iprob;
  Real x1,x2,x3;
 // Real time, dt;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        Real dvol= pGrid->dx1*pGrid->dx2*pGrid->dx3;
        Real rho = pGrid->U[k][j][i].d;
	Real prs = amb_rho * T_floor * Kboltz/muH*mue;
	Real temp = (prs/rho)*Kboltz*mu;
	  
	if (temp <=T_cold){
	  local_cold_mass += rho * dvol;
	}
	local_v3_sum += velz;
      }
    }
  }
  MPI_Allreduce(&local_cold_mass, &global_cold_mass, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_v3_sum, &global_v3_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  // MPI_Allreduce ??????????????????????????/
  Real v3_avg = global_v3_sum/(nx1*nx2*nx3);
  front_posn_new = x3min + global_cold_mass/(L1*L2*amb_rho*T_hot/T_floor);
  v_shift_t_new = -2.0 * (front_posn_new-front_posn_old)/sim_dt;
    
  if (fabs(v_shift_t_new) > (c_s_cap*c_s)) {
    v_shift_t_new /= fabs(v_shift_t_new);
    v_shift_t_new *= c_s*c_s_cap;
  }
  Real front_velocity = -1.0*v_shift_t_new;

//////////////////NGHOST?????
  int il = pGrid->is-nghost;
  int iu = pGrid->ie+nghost;
  int jl = pGrid->js-nghost;
  int ju = pGrid->je+nghost;
  int kl = pGrid->ks-nghost;
  int ku = pGrid->ke+nghost;
    
  if (time <= shift_start){
    front_posn_old = front_posn_new;
  }
  else {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
	  pGrid->U[k][j][i].M3 += pGrid->U[k][j][i].d*v_shift_t_new;
	  Real dE_kin = pow(velz+v_shift_t_new,2)- pow(velz,2);
	  pGrid->U[k][j][i].E += 0.5 * pGrid->U[k][j][i].d * dE_kin;
	}
      } 
    }
  }
}

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  
  read_input(pDomain);
  printf("Line 343\n");
  townsend_cooling(pDomain);
  printf("townsend cooling finished.\n");
  frame_shift(pDomain);
  printf("frame shift done.\n");
 
  Real cloud_chi = T_hot/T_floor;
  Real rho_cold = amb_rho * cloud_chi;

  Real A_KH = amp_KH * v_shear; /////////////!!!!!!!!!!!!!!!!!!!!!111
  Real k_x = 2*PI*knx_KH/L1;
  Real k_y = 2*PI*kny_KH/L2;
  
  Real x1,x2,x3;
  Real x1f,x2f,x3f;

  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        pGrid->U[k][j][i].d = rho_cold + (0.5*(amb_rho - rho_cold)) * (1 + tanh(x3/front_thick));
	
	pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d * (v_shear/2) * (1 + tanh(x3/front_thick));
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d * A_KH * exp(-1.0*x3*x3/front_thick/front_thick)*sin(k_x*x1) * sin(k_y*x2);
	pGrid->U[k][j][i].M3 += pGrid->U[k][j][i].d * v_shift0;

#ifndef BAROTROPIC
	pGrid->U[k][j][i].E = pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1;
	pGrid->U[k][j][i].E += pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2;
	pGrid->U[k][j][i].E += pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3;
	pGrid->U[k][j][i].E /= 2.0*pGrid->U[k][j][i].d;
	pGrid->U[k][j][i].E += (T_hot/(Kboltz*mu))*amb_rho/(g-1);
#endif
#if (NSCALARS > 0)
        pGrid->U[k][j][i].s[0] = 0.0;
        if (x2 > 0) pGrid->U[k][j][i].s[0] = 1.0;
#endif
      }
    }
  }
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
	pGrid->B1i[k][j][i] = B_x;
      }
    }
  }
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][j][i] = B_y;
      }
    }
  }
  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B3i[k][j][i] = B_z;
      }
    }
  }
#ifndef BAROTROPIC
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].E += 0.5*(SQR(B_x)+SQR(B_y)+SQR(B_z));
      }
    }
  }
#endif   
#endif
}

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

//static Real color(const GridS *pG, const int i, const int j, const int k)
//{
 // return pG->U[k][j][i].s[0]/pG->U[k][j][i].d;
////}
void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
#if (NSCALARS > 0)
  if(strcmp(expr,"color")==0) return color;
#endif
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(DomainS *pD)
{
  GridS *pGrid = pD->Grid;
  Real g;
  g = par_getd("problem","gamma");
  int i=0, j=0, k=0;
  int is,ie,js,je,ks,ke;
 //Real x1,x2,x3;
  Real dt;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	Real lum_cell = 0.0;
	Real prs = amb_rho * T_floor * Kboltz/muH*mue;
	Real temp = (prs/pGrid->U[k][j][i].d) * Kboltz * mu;
	if (temp > T_floor) {
	  if (temp < T_cut) {
	    Real rho = pGrid->U[k][j][i].d;
	    Real temp_cgs = temp;
	    Real rho_cgs = rho * unit_density;
	    Real dt_cgs = dt*unit_time;
	    Real cLfac = Lambda_fac;
	 
	    Real temp_new = fmax(Cooling_townsend(temp_cgs,rho_cgs,dt_cgs,cLfac),T_floor);
	    Real ccool_2 = ((temp_new-temp)/(Kboltz*mu))*pGrid->U[k][j][i].d/(g-1.0);
	    lum_cell -= ccool_2;
	  }
	}
	else{
	  Real ccool_2 = ((T_floor-temp)/(Kboltz*mu))*pGrid->U[k][j][i].d/(g-1.0);
	}
	return;
      }
    }
  }
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}
