#ifndef CODE_UNITS_H
#define CODE_UNITS_H
#include <stdio.h>
//#include <math.h>
#include "code_units.h"
#include <string.h>
#include <stdlib.h>
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#include "../defs.h"
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <float.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>

//#include "../athena.h"
//#include "../ath_array.c"


static double CONST_pc, CONST_yr, CONST_amu, CONST_kB;



static Real CONST_pc  = 3.086e18;
static Real CONST_yr  = 3.154e7;
static Real CONST_amu = 1.66053886e-24;
static Real CONST_kB  = 1.3806505e-16;

static Real unit_length = 3.086e18*1e3; // 1 kpc
static Real unit_time   = 3.154e7*1e6; // 1 Myr
static Real unit_density = 1.66053886e-24;   // 1 mp/cm-3

Real unit_velocity = (3.086e18*1e3)/(3.154e7*1e6); // in kpc/Myr
//float unit_q = (1.66053886e-24 * pow(((3.086e18*1e3)/(3.154e7*1e6)),3.0))/(3.086e18*1e3);

Real Kboltz = ((3.086e18*1e3)/(3.154e7*1e6))*((3.086e18*1e3)/(3.154e7*1e6))*1.66053886e-24/1.3806505e-16;

// in terms of solar abundances
static Real Zsol = 1.0;
static Real Xsol = 1.0;

double X = 1.0 * 0.7381;
double Z = 1.0 * 0.0134;
double Y = 1 - 0.7381 - 0.0134;

double mu  = 1.0/(2.*0.7381+ 3.*(1.-0.7381-0.0134)/4.+ 0.0134/2.);
double mue = 2.0/(1.0+0.7381);
double muH = 1.0/0.7381;

// Relevant temperatures
double T_floor   = 1.e4;   // in
double T_ceil    = 1.e8;   // in K
double T_hot     = 1.e7;
double T_hot_req = 1.e7;
double T_cold    = 1.e4;
double T_cut_mul = 0.6;
double T_cut     = 0.6*1.e7;


#endif
