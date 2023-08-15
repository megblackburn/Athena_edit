#ifndef EOS_EOS_H
#define EOS_EOS_H


#include <limits.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "../athena.h"
//#include "../ath_arrays.c"

struct Hydro{};
struct ParameterInput{};
struct FaceField{};

typedef struct EquationOfState{
	struct Hydro* pmb_block_;
	struct ParameterInput* pin 
} EquationOfState;

void ConservedToPrimitive(
	EquationOfState* eos, float* cons, const float* prim_old, const struct FaceField* b,
	float* prim, float* bcc,
	struct Coordinates* pco, int il, int iu, int jl, int ju, int kl, int ku);


struct Hydro;
struct ParameterInput;
struct FaceField;
struct EquationOfState {
  friend struct Hydro;
 public:
  EquationOfState(struct MeshBlock *pmb, struct ParameterInput *pin);
  void ConservedToPrimitive(
      struct AthenaArray_Real *cons, const struct AthenaArray_Real *prim_old, const struct FaceField *b,
      struct AthenaArray_Real *prim, struct AthenaArray_Real *bcc,
      struct Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);
  void PrimitiveToConserved(const struct AthenaArray_Real *prim, const struct AthenaArray_Real *bc,
                            struct AthenaArray_Real *cons, struct Coordinates *pco,
                            int il, int iu, int jl, int ju, int kl, int ku);
};

