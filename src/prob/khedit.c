#include "copyright.h"
/*============================================================================*/
/*! \file kh.c
 *  \brief Problem generator for KH instability. 
 *
 * PURPOSE: Problem generator for KH instability.  Sets up two versions:
 * - iprob=1: slip surface with random perturbations
 * - iprob=2: tanh profile at interface, with single-mode perturbation	      */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()
 *============================================================================*/

static double ran2(long int *idum);
#ifdef MHD
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

/*--------------- Read Input -----------------------------------*/
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
  
  rho_cold = par_getd("problem", "rho_cold");
  rho_hot = par_getd("problem", "rho_hot");

  g = par_getd("problem", "gamma");
  printf("input file read\n");

  return;
}

/*------------------ Townsend Cooling ---------------------*/
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


/*--------------- Frame Shift ----------------------------*/
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
                                 
/*-------------------- problem: -------------------------------*/

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,iprob;
  Real amp,drat,vflow,b0,a,sigma,x1,x2,x3;
  long int iseed = -1;
  static int frst=1;  /* flag so new history variables enrolled only once */

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Read problem parameters */

  iprob = par_geti("problem","iprob");
  vflow = par_getd("problem","vflow");
  drat = par_getd("problem","drat");
  amp = par_getd("problem","amp");
#ifdef MHD
  b0  = par_getd("problem","b0");
#endif

/* iprob=1.  Two uniform streams moving at +/- vflow, random perturbations */

  if (iprob == 1) {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          pGrid->U[k][j][i].d = rho_hot;
          pGrid->U[k][j][i].M1 = -vflow;// + amp*(ran2(&iseed) - 0.5);
          pGrid->U[k][j][i].M2 = amp*(ran2(&iseed) - 0.5);
          pGrid->U[k][j][i].M3 = 0.0;
          if (fabs(x2) < (L2/2.0)) {
  	    pGrid->U[k][j][i].d = rho_cold;
            pGrid->U[k][j][i].M1 = vflow + amp*(ran2(&iseed) - 0.5);
           // pGrid->U[k][j][i].M2 = drat*amp*(ran2(&iseed) - 0.5);
          }
/* Pressure scaled to give a sound speed of 1 with gamma=1.4 */
#ifndef BAROTROPIC
          pGrid->U[k][j][i].E = 2.5/Gamma_1
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* BAROTROPIC */
#ifdef MHD
          pGrid->B1i[k][j][i] = b0;
          pGrid->U[k][j][i].B1c = b0;
#ifndef BAROTROPIC
          pGrid->U[k][j][i].E += 0.5*b0*b0;
#endif /* BAROTROPIC */
#endif /* MHD */
        }
#ifdef MHD
      pGrid->B1i[k][j][ie+1] = b0;
#endif
      }
    }
  }

/* iprob=2.  Test suggested by E. Zweibel, based on Ryu & Jones.
 * Two uniform density flows with single mode perturbation
 */

  if (iprob == 2) {
    a = 0.05;
    sigma = 0.2;
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          pGrid->U[k][j][i].d = 1.0;
          pGrid->U[k][j][i].M1 = vflow*tanh(x2/a);
          pGrid->U[k][j][i].M2 = amp*sin(2.0*PI*x1)*exp(-(x2*x2)/(sigma*sigma));
          pGrid->U[k][j][i].M3 = 0.0;
#ifndef BAROTROPIC
          pGrid->U[k][j][i].E = 1.0/Gamma_1
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* BAROTROPIC */
#ifdef MHD
          pGrid->B1i[k][j][i] = b0;
          pGrid->U[k][j][i].B1c = b0;
#ifndef BAROTROPIC
          pGrid->U[k][j][i].E += 0.5*b0*b0;
#endif /* BAROTROPIC */
#endif /* MHD */
/* Use passive scalar to keep track of the fluids, since densities are same */
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = 0.0;
          if (x2 > 0) pGrid->U[k][j][i].s[0] = 1.0;
#endif
        }
#ifdef MHD
      pGrid->B1i[k][j][ie+1] = b0;
#endif
      }
    }
  }

/* With viscosity and/or resistivity, read diffusion coeffs */

#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  Q_Hall  = par_getd_def("problem","Q_H",0.0);
  Q_AD    = par_getd_def("problem","Q_AD",0.0);
#endif
#ifdef VISCOSITY
  nu_iso = par_getd_def("problem","nu_iso",0.0);
  nu_aniso = par_getd_def("problem","nu_aniso",0.0);
#endif

/* enroll new history variables, only once  */

  if (frst == 1) {
#ifdef MHD
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
#endif /* MHD */
    frst = 0;
  }

}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  Q_Hall  = par_getd_def("problem","Q_H",0.0);
  Q_AD    = par_getd_def("problem","Q_AD",0.0);
#endif
#ifdef VISCOSITY
  nu_iso = par_getd_def("problem","nu_iso",0.0);
  nu_aniso = par_getd_def("problem","nu_aniso",0.0);
#endif
  return;
}

#if (NSCALARS > 0)
/*! \fn static Real color(const GridS *pG, const int i, const int j,const int k)
 *  \brief Returns first passively advected scalar s[0] */
static Real color(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].s[0]/pG->U[k][j][i].d;
}
#endif

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

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*------------------------------------------------------------------------------
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum)
 *  \brief  Extracted from the Numerical Recipes in C (version 2) code. Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. 

 */
double ran2(long int *idum){
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX

/*------------------------------------------------------------------------------
 * MHD history variables
 * hst_Bx, etc.: Net flux, and Maxwell stress, added as history variables
 */

#ifdef MHD
/*! \fn static Real hst_Bx(const GridS *pG, const int i,const int j,const int k)
 *  \brief x-component of B-field */
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B1c;
}

/*! \fn static Real hst_By(const GridS *pG, const int i,const int j,const int k)
 *  \brief y-component of B-field */
static Real hst_By(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B2c;
}

/*! \fn static Real hst_Bz(const GridS *pG, const int i,const int j,const int k)
 *  \brief z-component of B-field */
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B3c;
}
#endif
