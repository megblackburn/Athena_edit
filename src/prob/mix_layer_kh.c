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

// parameters

Real tfloor;
Real tnotcool;
Real tcut_hst;
Real r_drop;
Real Lambda_fac;
Real Lambda_fac_time;
Real total_cooling = 0.0;

Real L1 = 0.0;
Real L2 = 0.0;
Real L3 = 0.0;

int nx1 = 1;
int nx2 = 1;
int nx3 = 1;

// ==1 -> cooling/shift included
int cooling_flag = 1;
int shift_flag = 1;

Real x3min = 0.0;
Real x3max = 0.0;

Real amb_rho = 1.0;
Real front_thick = 2.5;
Real v_shift0 = 0.0;
Real knx_KH = 1.0;
Real kny_KH = 1.0;
Real amp_KH = 0.01;

// Magnetic fields
Real B_x = 0.0;
Real B_y = 0.0;
Real B_z = 1.0;

Real front_posn_old = 0.0;

// debugging
_Bool cooling_flag_print_count = false;
_Bool shift_flag_print_count = false;
_Bool DEBUG_FLAG_MIX = false;

Real time_last = -1.0;
Real shift_start = 0.0;
Real debug_hst = 0.0;

// velocities
Real velx;
//Real vely;
//Real velz;

// initial densities in cold and hot layers
//Real rho_cold;
//Real rho_hot;
Real rho;

// C_p/C_v
Real g; //gamma

// scaling factors
Real a;
Real sigma;

void read_input(Domain *pDomain){
  L1 = par_getd("domain1", "x1max")-par_getd("domain1", "x1min");
  L2 = par_getd("domain1", "x2max")-par_getd("domain1", "x2min");
  L3 = par_getd("domain1", "x3max")-par_getd("domain1", "x3min");

  x3min = par_getd("domain1", "x3min");
  x3max = par_getd("domain1", "x3max");

  nx1 = par_getd("domain1", "Nx1");
  nx2 = par_getd("domain1", "Nx2"); 
  nx3 = par_getd("domain1", "Nx3");

  cooling_flag = par_geti("problem","cooling_flag");
  shift_flag = par_geti("problem", "shift_flag");
  shift_start = par_getd("problem", "shift_start");
  v_shift0 = par_getd("problem", "v_shift");
  amb_rho = par_getd("problem", "amb_rho");
 // v_shear = par_getd("problem", "v_shear"); 
  DEBUG_FLAG_MIX = par_geti("problem", "DEBUG_FLAG");

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

  rho_hot = par_getd("problem", "rho_hot");
  rho_cold = par_getd("problem", "rho_cold");

  velx = par_getd("problem", "velx");
  //vely = par_getd("problem", "vely");
  //velz = par_getd("problem", "velz");
  
  g = par_getd("problem", "gamma");

  // scaling factors
  a = par_getd("problem", "a");
  sigma = par_getd("problem", "sigma");
  printf("input file read! \n");

  return;
}

void townsend_cooling(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;

  int i=0, j=0, k=0;
  int is, ie, js, je, ks, ke;
  Real x1,x2,x3;
  Real time, dt;
  Real temp;
  Real ccool;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (k=ks; k<=ke; k++) {
        cc_pos(pGrid, i, j, k, &x1, &x2, &x3);
        debug_hst +=1.0;
	rho = rho_hot * rho_scale;
        pGrid->U[k][j][i].d = rho_hot * rho_scale;
        temp = T_hot;
        pGrid->U[k][j][i].M1 = velx * tan(x2/a); // v_x
        pGrid->U[k][j][i].M2 = velx * amb_rho*sin(2.0*PI*x1)*exp(-(x2*x2)/(sigma*sigma)); //v_y
        pGrid->U[k][j][i].M3 = velx;
          
        if (fabs(x2) < (L1/2)) {   
	  rho = rho_cold * rho_scale;
          pGrid->U[k][j][i].d = rho_cold *rho_scale;
          temp = T_cold;
          pGrid->U[k][j][i].M1 = -rho_cold * velx * tan(x2/a); // v_x
          pGrid->U[k][j][i].M2 = rho_cold * amb_rho*sin(2.0*PI*x1)*exp(-(x2*x2)/(sigma*sigma)); //v_y
        }
          
        if (isnan(temp)) { 
          printf("temp\n", temp);
          printf("rho\n", rho);

          Real E_kin = 0.5 * (SQR(pGrid->U[k][j][i].M1)+SQR(pGrid->U[k][j][i].M2)+SQR(pGrid->U[k][j][i].M3)*pGrid->U[k][j][i].d);
          Real E_mag = 0.0;
#ifdef MHD
          E_mag = 0.5 * SQR(B_x)+SQR(B_y)+SQR(B_z);
#endif
          pGrid->U[k][j][i].E = E_kin + E_mag + (T_floor/(Kboltz*mu)) * (pGrid->U[k][j][i].d/(g-1.0));
        }

        if (temp > T_floor) { 
          if (temp < T_cut) { 
            Real rho_cgs = rho;
            Real temp_cgs = rho_cgs * unit_density;
            Real dt_cgs = dt * unit_time;
            Real cLfac = Lambda_fac;
            Real temp_new = fmax(Cooling_townsend(temp_cgs, rho_cgs, dt_cgs, cLfac), T_floor);

            if (isnan(temp_new)) {
              printf("temp \n", temp_cgs);
              ccool = ((T_floor-temp)/(Kboltz*mu))*pGrid->U[k][j][i].d/(g-1.0);
              pGrid->U[k][j][i].E += ccool;
            }
            else{
              ccool = ((temp_new - temp)/(Kboltz*mu))*pGrid->U[k][j][i].d/(g-1.0);
              pGrid->U[k][j][i].E += ccool;
              printf("temp < T_floor\n");
              if (DEBUG_FLAG_MIX) { 
                total_cooling -= 1.0;
              }
              else{
                total_cooling -=ccool;
              }
            }  
          }
        }
        else{
          ccool = ((T_floor-temp)/(Kboltz*mu))*pGrid->U[k][j][i].d/(g-1.0);
          pGrid->U[k][j][i].E += ccool;
        }
      }
    }
  }
}

void frame_shift(DomainS *pDomain) { 
  GridS *pGrid = pDomain->Grid;

  printf("frame shift initiated\n");
  Real time, dt, g;
  if (time_last==time) {
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

  Real c_s = sqrt(g * T_floor/(mu*Kboltz));
  Real c_s_cap = 100.0;
  int i=0, j=0, k=0;
  int is,ie,js,je,ks,ke,iprob;
  Real x1,x2,x3;

  g = par_getd("problem", "gamma");

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

        Real dvol = pGrid->dx1*pGrid->dx2*pGrid->dx3;
        Real rho = pGrid->U[k][j][i].d;
        Real prs = amb_rho * T_floor * Kboltz/(muH*mue); //rho or amb_rho??????????
        Real temp = (prs/rho) * Kboltz * mu;

        if (temp <=T_cold){
          local_cold_mass += rho * dvol;
	      }
	      local_v3_sum += velx;
      }
    }
  }

  MPI_Allreduce(&local_cold_mass, &global_cold_mass, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_v3_sum, &global_v3_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  Real v3_avg = global_v3_sum/(nx1*nx2*nx3);
  front_posn_new = x3min + global_cold_mass/(L1*L2*amb_rho*T_hot/T_floor);
  v_shift_t_new = -2.0 * (front_posn_new-front_posn_old)/sim_dt;

  if (fabs(v_shift_t_new) > (c_s_cap*c_s)) {
    v_shift_t_new /= fabs(v_shift_t_new);
    v_shift_t_new *= c_s*c_s_cap;
  }
  Real front_velocity = -1.0*v_shift_t_new;

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
	  Real dE_kin = pow(velx+v_shift_t_new,2)- pow(velx,2);
	  pGrid->U[k][j][i].E += 0.5 * pGrid->U[k][j][i].d * dE_kin;
	}
      } 
    }
  }
}

void problem(DomainsS *pDomain)
{
  GridS *pGrid = pDomain->Grid;

// gets input parameters from the input file
  read_input(pDomain);

// parameter definitions
  Real cloud_chi = T_hot/T_floor;
  Real rho_cold = amb_rho * cloud_chi;

  Real A_KH = amp_KH * velx;
  Real k_x = 2*PI*knx_KH/L1;
  Real k_y = 2*PI*kny_KH/L2;
  
  Real x1,x2,x3;
  Real x1f,x2f,x3f;

  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
// start of problem
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	if (fabs(x2) < (L1/2)){
	  pGrid->U[k][j][i].d = rho_cold + (0.5*(amb_rho - rho_cold)) * (1 + tanh(x3/front_thick));
	  pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d * (velx/2) * (1 + tanh(x3/front_thick));
	  pGrid->U[k][j][i].M2 = 0.0;
	  pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d * A_KH * exp(-1.0*x3*x3/front_thick/front_thick)*sin(k_x*x1) * sin(k_y*x2);
	  pGrid->U[k][j][i].M3 += pGrid->U[k][j][i].d * v_shift0;
	}
	if (fabs(x2) > (L1/2)) {
	  pGrid->U[k][j][i].d = rho_hot + (0.5*(amb_rho - rho_hot)) * (1 + tanh(x3/front_thick));
	  pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d * (velx/2) * (1 + tanh(x3/front_thick));
	  pGrid->U[k][j][i].M2 = 0.0;
	  pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d * A_KH * exp(-1.0*x3*x3/front_thick/front_thick)*sin(k_x*x1) * sin(k_y*x2);
	  pGrid->U[k][j][i].M3 += pGrid->U[k][j][i].d * v_shift0;
	}
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
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}


            

