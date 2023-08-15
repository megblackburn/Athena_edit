#ifndef CCOLING_CLASS
#define CCOLING_CLASS
#include "../athena.h"
#include "../ath_array.c"
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>


//typedef int bool;
//#define true 1
//#define false 0

typedef double Real;

bool DEBUG_FLAG = false;

typedef struct Cooling{
   Real tfloor;
   Real (*townsend)(Real temp, Real rho, Real const dt, Real Lambda_fac);
   Real (*Lambda)(Real temp);
   Real (*tcool)(Real temp, Real rho);
   Real const mh;
   Real const kb;
   int nbins;
   Real const_factor;
   Real *cool_t;
   Real *cool_tef;
   Real *cool_coef;
   Real *cool_index;
} Cooling;

Cooling* Cooling_new();
void Cooling_delete(Cooling* cooling);
Real Cooling_townsend(Real temp, Real rho, Real dt, Real Lambda_fac);
Real Cooling_Lambda(Real temp);
Real Cooling_tcool(Real temp, Real rho);


//Cooling* Cooling_new() {
  // Cooling* Cooling = ((Cooling*)malloc(sizeof(Cooling));
//}
double cooling_ctr(double cooling()){
   if (DEBUG_FLAG) {
      printf("Start of Cooling::Cooling!\n");
   }

    Real g  = 5.0/3.0; 
    Real X = 0.7; auto Z = 0.02; 
    Real mu   = 1.0/(2.0*X + 0.75*(1.0-X-Z) + Z/2.0);
    Real mu_e = 2.0/(1.0+X);
    Real mu_h = 1.0/X;
    cooling->const_factor = (1.0e-23)*(g-1.0)*mu/(cooling->kb*mu_e*mu_h*cooling->mh);
  
    cooling->nbins = 40;
    cooling->cool_t = (Real*)malloc(cooling->nbins * sizeof(Real));
    cooling->cool_tef = (Real*)malloc(cooling->nbins * sizeof(Real));
    cooling->cool_coef = (Real*)malloc(cooling->nbins * sizeof(Real));
    cooling->cool_index = (Real*)malloc(cooling->nbins * sizeof(Real));
  
    cooling->cool_t[0]  = 10000.0;
    cooling->cool_t[1]  = 12600.0;
    cooling->cool_t[2]  = 15890.0;
    cooling->cool_t[3]  = 20020.0;
    cooling->cool_t[4]  = 25240.0;
    cooling->cool_t[5]  = 31810.0;
    cooling->cool_t[6]  = 40090.0;
    cooling->cool_t[7]  = 50530.0;
    cooling->cool_t[8]  = 63680.0;
    cooling->cool_t[9]  = 80260.0;
    cooling->cool_t[10] = 101200.0;
    cooling->cool_t[11] = 127500.0;
    cooling->cool_t(12) = 160700.0;
    cooling->cool_t(13) = 202600.0;
    cooling->cool_t(14) = 255300.0;
    cooling->cool_t(15) = 321800.0;
    cooling->cool_t(16) = 405600.0;
    cooling->cool_t(17) = 511100.0;
    cooling->cool_t(18) = 644200.0;
    cooling->cool_t(19) = 812000.0;
    cooling->cool_t(20) = 1000000.0;
    cooling->cool_t(21) = 1259000.0;
    cooling->cool_t(22) = 1585000.0;
    cooling->cool_t(23) = 1995000.0;
    cooling->cool_t(24) = 2512000.0;
    cooling->cool_t(25) = 3162000.0;
    cooling->cool_t(26) = 3981000.0;
    cooling->cool_t(27) = 5012000.0;
    cooling->cool_t(28) = 6310000.0;
    cooling->cool_t(29) = 7943000.0;
    cooling->cool_t(30) = 10000000.0;
    cooling->cool_t(31) = 12590000.0;
    cooling->cool_t(32) = 15850000.0;
    cooling->cool_t(33) = 19950000.0;
    cooling->cool_t(34) = 25120000.0;
    cooling->cool_t(35) = 31620000.0;
    cooling->cool_t(36) = 39810000.0;
    cooling->cool_t(37) = 50120000.0;
    cooling->cool_t(38) = 63100000.0;
    cooling->cool_t(39) = 79430000.0; 

    cooling->cool_coef(0) = 1.6408984689285624;
    cooling->cool_coef(1) = 5.789646575948292;
    cooling->cool_coef(2) = 18.797203755396648;
    cooling->cool_coef(3) = 16.7384754689852;
    cooling->cool_coef(4) = 11.274384717759935;
    cooling->cool_coef(5) = 9.95038422958871;
    cooling->cool_coef(6) = 11.302144847043829;
    cooling->cool_coef(7) = 15.819149070534786;
    cooling->cool_coef(8) = 25.224636283348048;
    cooling->cool_coef(9) = 38.02107089248533;
    cooling->cool_coef(10) = 43.98219098299675;
    cooling->cool_coef(11) = 41.277704007796586;
    cooling->cool_coef(12) = 41.95311185975414;
    cooling->cool_coef(13) = 45.260670345801;
    cooling->cool_coef(14) = 47.275626188961176;
    cooling->cool_coef(15) = 32.21420131907784;
    cooling->cool_coef(16) = 24.350976818250636;
    cooling->cool_coef(17) = 23.383616480583676;
    cooling->cool_coef(18) = 18.333394532081098;
    cooling->cool_coef(19) = 14.89691888284402;
    cooling->cool_coef(20) = 14.392505898454834;
    cooling->cool_coef(21) = 13.027915287005817;
    cooling->cool_coef(22) = 11.671262753284271;
    cooling->cool_coef(23) = 9.070904785425046;
    cooling->cool_coef(24) = 6.489695397654223;
    cooling->cool_coef(25) = 4.766239129792971;
    cooling->cool_coef(26) = 3.7811870710765074;
    cooling->cool_coef(27) = 3.313622783657129;
    cooling->cool_coef(28) = 3.0600313080475674;
    cooling->cool_coef(29) = 2.9993768457216112;
    cooling->cool_coef(30) = 2.9491332141250552;
    cooling->cool_coef(31) = 2.744653611808266;
    cooling->cool_coef(32) = 2.3449511265716;
    cooling->cool_coef(33) = 2.0169621177549892;
    cooling->cool_coef(34) = 1.8907205849384978;
    cooling->cool_coef(35) = 1.91584885606706;
    cooling->cool_coef(36) = 2.056870288868004;
    cooling->cool_coef(37) = 2.233680315878366;
    cooling->cool_coef(38) = 2.4097186710383474;
    cooling->cool_coef(39) = 2.5502102007949023;

    cooling->cool_index(0) = 5.455488390256632;
    cooling->cool_index(1) = 5.076170519863754;
    cooling->cool_index(2) = -0.5020655826640291;
    cooling->cool_index(3) = -1.7055659800651979;
    cooling->cool_index(4) = -0.5399688186820728;
    cooling->cool_index(5) = 0.550609170202909;
    cooling->cool_index(6) = 1.4527662908446985;
    cooling->cool_index(7) = 2.0172644735605223;
    cooling->cool_index(8) = 1.773197476674277;
    cooling->cool_index(9) = 0.6282445620956022;
    cooling->cool_index(10) = -0.2747076405016009;
    cooling->cool_index(11) = 0.07013182420220869;
    cooling->cool_index(12) = 0.32752568568776125;
    cooling->cool_index(13) = 0.1883881016798681;
    cooling->cool_index(14) = -1.6570303755459093;
    cooling->cool_index(15) = -1.209120245966656;
    cooling->cool_index(16) = -0.17533183860418153;
    cooling->cool_index(17) = -1.0512755674245657;
    cooling->cool_index(18) = -0.896664392554265;
    cooling->cool_index(19) = -0.16540667885641686;
    cooling->cool_index(20) = -0.43250361812273735;
    cooling->cool_index(21) = -0.4775539072045259;
    cooling->cool_index(22) = -1.0956186284443203;
    cooling->cool_index(23) = -1.453147878451421;
    cooling->cool_index(24) = -1.3412596915753237;
    cooling->cool_index(25) = -1.0051719479026813;
    cooling->cool_index(26) = -0.573142729390977;
    cooling->cool_index(27) = -0.3457087236213044;
    cooling->cool_index(28) = -0.08698732111048613;
    cooling->cool_index(29) = -0.07335511773234596;
    cooling->cool_index(30) = -0.3119882060952377;
    cooling->cool_index(31) = -0.6835132944311395;
    cooling->cool_index(32) = -0.6549261784681947;
    cooling->cool_index(33) = -0.2804886559029823;
    cooling->cool_index(34) = 0.05737205818565948;
    cooling->cool_index(35) = 0.30836313806582183;
    cooling->cool_index(36) = 0.3580735000106496;
    cooling->cool_index(37) = 0.3293929876114671;
    cooling->cool_index(38) = 0.24620665148692336;
    cooling->cool_index(39) = 0.10953503955831644;


  cool_tef[nbins-1] = 0.0;
  for (int i=nbins-2; i>=0; i--) {
    double t_n    = cool_t[nbins-1];
    double coef_n = cool_coef[nbins-1];
    double t_i   = cool_t[i];
    double tef_i = cool_tef[i];
    double coef  = cool_coef[i];
    double slope = cool_index[i];
    double sm1  = slope - 1.0;
    double step = (coef_n/coef)*(t_i/t_n)*(pow(t_i/cool_t[i+1],sm1)-1)/sm1;
    cool_tef[i] = cool_tef[i+1] - step;
  }
  tfloor = cool_t[0];
  if (DEBUG_FLAG){
    printf("End of cooling:Cooling!\n")
  }

}
void cooling_dtr(struct cooling* obj) {
  free(cool_t);
  free(cool_tef);
  free(cool_coef);
  free(cool_index);
}

double cooling_ctr(Real townsend(Real temp, Real rho, Real const dt, Real Lambda_fac)) {
  if (DEBUG_FLAG){
    printf("Start of Cooling:Townsend!\n");
  if (DEBUG_FLAG){
    printf("Get Reference values from the last bin\n");
  }
  double t_n    = cool_t[nbins-1];
  double coef_n = cool_coef[nbins-1] * Lambda_fac;
  if (DEBUG_FLAG){
    printf("Get the index of the right temperature bin\n");
  }
  int idx = 0;
  while ((idx < nbins-2) && (cool_t[idx+1] < temp)) { idx += 1; }
  if (DEBUG_FLAG){
    printf("Look up the corresponding slope and coefficient\n");
  }
  double t_i   = cool_t[idx];
  double tef_i = cool_tef[idx];
  double coef  = cool_coef[idx] * Lambda_fac;
  double slope = cool_index[idx];

  if (DEBUG_FLAG){
    printf("Compute the temporal evolution function Y(T)\n");
  }
  double sm1 = slope - 1.0;
  double tef = tef_;
  if (DEBUG_FLAG){
    printf("compute the adjusted TEF for new timestep\n");
  }
  double tef_adj = tef + rho*coef_n*const_factor*dt/t_n;
  
  if (DEBUG_FLAG){
    print("TEF is a strictly decreasing function and new_tef > tef\n");
  }
  while ((idx > 0) && (tef_adj > cool_tef[idx])) {
    idx -= 1;
    t_i   = cool_t[idx];
    tef_i = cool_tef[idx];
    coef  = cool_coef[idx] * Lambda_fac;
    slope = cool_index[idx];
  }
  if (DEBUG_FLAG){
    print("compute the inverse temporal evolution function Y^{-1}(Y)\n");
  }
  double oms  = 1.0-slope;
  double tnew = t_i*pow(1-oms*(coef/coef_n)*(t_n/t_i)*(tef_adj-tef_i),1/oms);
  if (DEBUG_FLAG){
    printf("End of Cooling:townsend!/n");
  }

  return tnew;
}



double cooling_ctr(double Lambda(double temp)) {
  if (DEBUG_FLAG){
    printf("Start of Cooling:Lambda\n");
  }
  if((temp< cool_t[0]) || (temp > cool_t[nbins-1]))
    return NAN;
  int idx = 0;
  while ((idx < nbins-2) && (cool_t[idx+1] < temp)) { idx += 1; }
  
  double t_i = cool_t[idx];
  double coef = cool_coef[idx];
  double slope = cool_index[idx];
  
  if (DEBUG_FLAG){
    print("End of Cooling:Lambda\n");
  }
  return coef * pow(temp / t_i, slope);
}


double cooling_ctr(double tcool(double temp, double rho)){
  if (DEBUG_FLAG){
    printf("start of Cooling:tcool\n");
  }
  if((temp < cool_t[0]) || (temp > cool_t[nbins - 1]))
      return NAN;
  double L = Lambda(temp);
  
  if (DEBUG_FLAG){
    printf("end of Cooling:tcool\n");
  }
  return temp / (const_factor * rho * L);
}
#endif


