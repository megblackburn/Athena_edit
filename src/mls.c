#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "defs.h"

//#include "./utils/hst_func.hpp"
#include "./utils/code_units.h"
#include "./utils/t_c_m.h"

//algorithm and memory and omp.h headers

static. Real tfloor;
stati/c Real tnotcool;
static Real tcut_hst;
static Real r_drop;
static Real Lambda_fac;
static Real Lambda_fac_time;
static Real total_cooling = 0.0;

typedef struct T{
    int value;
} T;

T *make_unique_T(int value) {
    T *ptr = (T *)malloc(sizeof(T));
    if (ptr) {
        ptr->value = value;
    }
    return ptr;
}


T *obj = make_unique_T(T);
    if (obj) {
        printf("Value: %d\n", obj->value);
        free(obj);
    } else {
        printf("Memory allocation failed\n");
    }

    return 0;
}

typedef struct cooler{
    int value;
} cooler;

cooler *make_unique_cooler(int value) {
    cooling *ptr = (cooler *)malloc(sizeof(cooler));
    if (ptr) {
        ptr->value = value;
    }
    return ptr;
}


cooler *obj = make_unique_cooler(cooler);
    if (obj) {
        printf("Value: %d\n", obj->value);
        free(obj);
    } else {
        printf("Memory allocation failed\n");
    }

    return 0;
}

static Real L1 = 0.0;
static Real L2 = 0.0;
static Real L3 = 0.0;

static int nx1 = 1;
static int nx2 = 1;
static int nx3 = 1;


static int cooling_flag = 1;
static int shift_flag  = 1;

static Real x3min = 0.0;
static Real x3max = 0.0;
static Real amb_rho = 1.0;
static Real front_thick = 2.5;
static Real v_shear = 100 * (1.023*1e-3);
static Real v_shift0 = 0.0;
static Real front_velocity = 0.0;
static Real knx_KH = 1.0;
static Real kny_KH = 1.0;
static Real amp_KH = 0.01;
static Real B_x = 0.0;
static Real B_y = 0.0;
static Real B_z = 1.0;
static Real front_posn_old = 0.0;
static bool cooling_flag_print_count = false;
static bool shift_flag_print_count = false;
static Real time_last = -1.0;
static Real shift_start = 0.0;
static Real debug_hst = 0.0;
static Real T_floor;
static Real T_cold;
static Real T_ceil;
static Real T_hot;
static Real T_cut_mul;
// new PARAMETERS
//static Real Kbolts = 1.38e-16
static Real velx,vely,velz;
static Real rho;
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

  mu = 1.0/(2.*X+3.*(1.-X-Z)/4.+Z.2);
  mue = 2.0/(1.0+X);
  muH = 1.0/X;

  T_floor = par_getd("problem", "T_floor");
  T_ceil = par_getd("problem", "T_ceil");
  T_hot = par_getd("problem", "T_hot");
  T_cold = par_getd("problem", "T_cold");
  T_cut_mul = par_getd("problem", "T_cut_mul");

  gamma = par_getd("problem", "gamma");
  printf("input file read\n");

  return;
}

void townsend_cooling(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0, j=0, k=0;
  int is,ie,js,je,ks,ke,iprob;
  Real time, dt

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  iprob = par_geti("problem","iprob");
  if (iprob ==1) {
    Real a=0.05;
    Real sigma=0.2;
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
	  cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          debug_hst += 1.0;
////////////// UNSURE IF NEEDED HERE
          pGrid->U[k][j][i].d = rho;
          pGrid->U[k][j][i].M1 = vflow*tanh(x2/a);
	  pGrid->U[k][j][i].M2 = amb_rho*sin(2.0*PI*x1)*exp(-(x2*x2)/(sigma*sigma));
          pGrid->U[k][j][i].M3 = 0.0;

	  Real temp = (gamma-1.0)*rho*mu/(Kboltz*mue*muH)
	  // Real temp = pres/rho *Kboltz*mu;

	  if (isnan(temp)){
	    printf("temp \n",temp);
	    printf("rho \n", rho);
	    
	    //Real E_kin = (Kbolts * temp * rho) / ((gamma-1.0)*mp*mh);
	    Real E_kin = 0.5 * (SQR(pGrid->U[k][j][i].M1)+SQR(pGrid->U[k][j][i].M2)+SQR(pGrid->U[k][j][i].M3))* pGrid->U[k][j][i].d;
	    Real E_mag = 0.0;
	    
#ifdef MHD
	    E_mag = 0.5 *B_x*B_x + B_y*B_y +B_z*B_z
#endif
	    pGrid->U[k][j][i].E = E_kin + E_mag + (T_floor/(Kboltz*mu))*pGrid->U[k][j][i].d/(gamma-1.0);
	  }
          if (temp > T_floor) {
	    if (temp<T_cut){
	      Real rho = rho;
	      Real temp_cgs = temp;
	      Real rho_cgs = rho * unit_density; // NEED CODE_UNITS
	      Real dt_cgs = dt *unit_time;
	      Real cLfac = Lambda_fac;
//////////////// DEFINE TOWNSEND!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	      Real temp_new = MAX(cooling_ctr(townsend(temp_cgs,rho_cgs,dt_cgs,cLfac), T_floor));
	      if (isnan(temp_new)) {
		printf("temp \n", temp_cgs);
		
		Real ccool = ((T_floor-temp)/(Kboltz*mu))*pGrid->U[k][j][i].d/(gamma-1.0);
		pGrid->U[k][j][i].E += ccool;
	      else {
		Real ccool = ((temp_new-temp)/(Kboltz*mu))*pGrid->U[k][j][i].d/(gamma-1.0);
		pGrid->U[k][j][i].E += ccool;
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
	    Real ccool = ((T_floor-temp)/(Kboltz*mu))* pGrid->U[k][j][i].d/(gamma-1.0);
	    pGrid->U[k][j][i].E += ccool;
	  }
        }
      }
    }
  }
}

void frame_shift(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  Real time, dt, g;
  
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
  
  Real front_posn_new = 0.0;
  Real v_shift_t_new = 0.0;
  
  Real sim_dt = pGrid->dt;
  
  g = par_getd("problem", "gamma")
  Real c_s = sqrt(g * T_floor/(mu*Kboltz));
  Real c_s_cap = 100.0;
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) 
	cc_pos(pG,i,j,k,&x1,&x2,&x3);
        Real dvol= pG->dx1*pG->dx2*pG->dx3;
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
  MPI_Allreduce(&local_cold_mass, &global_cold_mass, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_v3_sum, &global_v3_sum, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  // MPI_Allreduce ??????????????????????????/
  Real v3_avg = global_v3_sum/Real(nx1*nx2*nx3);
  front_posn_new = x3min + global_cold_mass/(L1*L2*amb_rho*T_hot/T_floor);
  v_shift_t_new = -2.0 * (front_posn_new-front_posn_old)/sim_dt;
    
  if (abs(v_shift_t_new) > (c_s_cap*c_s)){
    v_shift_t_new /= abs(v_shift_t_new);
    v_shift_t_new *= c_s*c_s_cap;
  }
  front_velocity = -1.0*v_shift_t_new;

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
	  U[k][j][i].M3 += U[k][j][i].d*v_shift_t_new;
	  Real dE_kin = pow(velz+v_shift_t_new,2)- pow(velz,2);
	  U[k][j][i].E += 0.5 * U[k][j][i].d * dE_kin;
	}
      } 
    }
  }
}

void problem(DomainS *pDomain)
{
  Grids *pGrid = pDomain->Grid;
  
  read_input(pDomain);
  
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

#ifndef NONBAROTROPIC
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
#ifndef NONBAROTROPIC
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].E += 0.5*(SQR(B_x)+SQR(B_y)+SQR(B_z));
      }
    }
  }
  
}

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

static Real color(const Grids *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].s[0]/pG->U[k][j][i].d;
}

void Userwork_in_loop(Mesh *pM)
{
  Real gamma;
  gamma = par_getd("problem","gamma");
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	Real lum_cell = 0.0;
	Real prs = amb_rho * T_floor * Kboltz/muH*mue;
	Real temp = (prs/pM->U[k][j][i].d) * Kboltz * mu;
	if (temp > T_floor) {
	  if (temp < T_cut) {
	    Real rho = pM->U[k][j][i].d;
	    Real temp_cgs = temp;
	    Real rho_cgs = rho * unit_density;
	    Real dt_cgs = dt*unit_time;
	    Real cLfac = Lambda_fac;
	 
	    Real temp_new = MAX(cooling_ctr(townsend(temp_cgs,rho_cgs,dt_cgs,cLfac),T_floor));
	    Real ccool_2 = ((temp_new-temp)/(Kboltz*mu))*pM->U[k][j][i].d/(gamma-1.0);
	    lum_cell -= ccool_2;
	  }
	}
	else{
	  Real ccool_2 = ((T_floor-temp)/(Kboltz*mu))*pM->U[k][j][i].d/(gamma-1.0);
	}
      }
      return lum_cell/dt;
    }
  }
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}
