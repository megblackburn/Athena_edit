#ifndef EOS_EOS_H_
#define EOS_EOS_H_

#include <limits.h>
#include "../athena.h"
#include "../athena_arrays.h"
#include "../coordinates/coordinates.h"
#include "../utils/interp_table.h"

struct Hydro;
struct ParameterInput;
struct FaceField;

struct EquationOfState {
    struct Hydro *pmy_block_;
    Real iso_sound_speed_;
    Real gamma_;
    Real density_floor_;
    Real pressure_floor_;
    Real energy_floor_;
    Real scalar_floor_;
    Real sigma_max_;
    Real beta_min_;
    Real gamma_max_;
    Real rho_min_;
    Real rho_pow_;
    Real pgas_min_;
    Real pgas_pow_;
    Real rho_unit_;
    Real inv_rho_unit_;
    Real egas_unit_;
    Real inv_egas_unit_;
    Real vsqr_unit_;
    Real inv_vsqr_unit_;
    struct AthenaArray g_;
    struct AthenaArray g_inv_;
    struct AthenaArray normal_dd_;
    struct AthenaArray normal_ee_;
    struct AthenaArray normal_mm_;
    struct AthenaArray normal_bb_;
    struct AthenaArray normal_tt_;
};

void EquationOfState_InitEosConstants(struct EquationOfState *eos, struct ParameterInput *pin);

void EquationOfState_ConservedToPrimitive(
    struct EquationOfState *eos, struct AthenaArray *cons, const struct AthenaArray *prim_old,
    const struct FaceField *b, struct AthenaArray *prim, struct AthenaArray *bcc,
    struct Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);

void EquationOfState_PrimitiveToConserved(
    struct EquationOfState *eos, struct AthenaArray *prim, struct AthenaArray *bc,
    struct AthenaArray *cons, 
    struct Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);

void EquationOfState_ConservedToPrimitiveCellAverage(
    struct EquationOfState *eos, struct AthenaArray *cons, struct AthenaArray *prim_old,
    struct FaceField *b, struct AthenaArray *prim, struct AthenaArray *bcc, 
    struct Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);

void EquationOfState_PassiveScalarConservedToPrimitive(
    struct EquationOfState *eos, struct AthenaArray *cons, struct AthenaArray *prim_old,
    struct AthenaArray *prim, struct AthenaArray *bcc,
    struct FaceField *b,
    struct Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);

void EquationOfState_PassiveScalarPrimitiveToConserved(
    struct EquationOfState *eos, struct AthenaArray *r, struct AthenaArray *u, 
    struct AthenaArray *s,
    struct Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);

void EquationOfState_PassiveScalarConservedToPrimitiveCellAverage(
    struct EquationOfState *eos, struct AthenaArray *s, struct AthenaArray *r_old, 
    struct AthenaArray *r
    struct Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);

void EquationOfState_ApplyPrimitiveFloors(
    struct EquationOfState *eos, struct AthenaArray *s, int n, int k, int j, int i);

void EquationOfState_ApplyPassiveScalarFloors(
    struct EquationOfState *eos, struct AthenaArray *s, int n, int k, int j, int i);

void EquationOfState_ApplyPassiveScalarPrimitiveConservedFloors(
    struct EquationOfState *eos, struct AthenaArray *s, struct AthenaArray *w, 
    struct AthenaArray *r, int n, int k, int j, int i);


// Other function declarations...

Real EquationOfState_PresFromRhoEg(struct EquationOfState *eos, Real rho, Real egas);
Real EquationOfState_EgasFromRhoP(struct EquationOfState *eos, Real rho, Real pres);
Real EquationOfState_AsqFromRhoP(struct EquationOfState *eos, Real rho, Real pres);
Real EquationOfState_GetIsoSoundSpeed(struct EquationOfState *eos);
Real EquationOfState_GetDensityFloor(struct EquationOfState *eos);
Real EquationOfState_GetPressureFloor(struct EquationOfState *eos);
struct EosTable* ptable;
#if GENERAL_EOS
Real EquationOfState_GetGamma(struct EquationOfState *eos);
#else
Real EquationOfState_GetGamma(struct EquationOfState *eos);
#endif

#endif // EOS_EOS_H_
