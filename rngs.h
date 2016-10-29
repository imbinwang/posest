/////////////////////////////////////////////////////////////////////////////////
// 
//  Reentrant random number generators (RNGs)
//
//  Copyright (C) 2002 - 2013  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _RNGS_H
#define _RNGS_H

#ifdef _MSC_VER
typedef unsigned __int32 uint32_t;
#else // assume a Unix-like system
#include <stdint.h>
#endif

#define RNG_USE_KISS32      1 // use KISS32
#define RNG_USE_KISS32a     0 // use KISS32a
#define RNG_USE_JKISS       0 // use JKISS
#define RNG_USE_LFSR113     0 // use LFSR113

#if RNG_USE_KISS32 + RNG_USE_KISS32a + RNG_USE_JKISS + RNG_USE_LFSR113 != 1
#error Exactly one of the RNG_USE_XXX macros should be defined as 1!
#endif


/* repeatable random sequences or not? */
#define RNG_REPEATABLE      0
#define RNG_NONREPEATABLE   1


#if RNG_USE_KISS32==1
#define rng_state rng_kiss32_state
#define rng_next  rng_kiss32
#elif RNG_USE_KISS32a==1
#define rng_state rng_kiss32a_state
#define rng_next  rng_kiss32a
#elif RNG_USE_JKISS==1
#define rng_state rng_jkiss_state
#define rng_next  rng_jkiss
#elif RNG_USE_LFSR113==1
#define rng_state rng_lfsr113_state
#define rng_next  rng_lfsr113
#else
#error if/else error in rngs.h!
#endif

// random integer between 0 and n-1 inclusive
#define rng_rint(state, n)    (rng_next(state)%(n))
#define rng_rint2(state, n)    ((int)((double)(n)*rng_next(state)*2.3283064365387E-10)) // 1/(2^32); uses floating point arithmetic

// random double in [0, 1); better variant as rng_stduniform()
#define rng_rdbl01(state)     (rng_next(state)*2.3283064365387E-10) // 1/(2^32)

#ifdef __cplusplus
extern "C" {
#endif

struct rng_lfsr113_state;
struct rng_kiss32_state;
struct rng_kiss32a_state;
struct rng_jkiss_state;

/* user-callable functions */
void rng_init(int howtoseed, struct rng_state *state);
uint32_t rng_next(struct rng_state *state);

/* misc functions */
double rng_stduniform(struct rng_state *state);
double rng_stdnormal (struct rng_state *state);
#define rng_uniform(state, low, high) ( (low) + rng_stduniform(state)*((high) - (low)) )
void rng_shuffle(struct rng_state *state, int arr[], int n);


/* actual RNGs */

/* Ecuyer's generator */
struct rng_lfsr113_state{
  uint32_t z1, z2, z3, z4;
  int seeded;
};
void rng_lfsr113_seed(struct rng_lfsr113_state *state);
uint32_t  rng_lfsr113(struct rng_lfsr113_state *state);

/* Marsaglia's KISS generator */
struct rng_kiss32_state{
  uint32_t mwc_up, mwc_lo, shr3, cong;
  int seeded;
};
void rng_kiss32_seed(struct rng_kiss32_state *state);
uint32_t  rng_kiss32(struct rng_kiss32_state *state);

/* faster Marsaglia's KISS variant */
struct rng_kiss32a_state{
  uint32_t x, y, z, w, c;
  int seeded;
};
void rng_kiss32a_seed(struct rng_kiss32a_state *state);
uint32_t  rng_kiss32a(struct rng_kiss32a_state *state);

/* Jones JKISS() Marsaglia variant */
struct rng_jkiss_state{
  unsigned int x, y, z, c;
  int seeded;
};

void    rng_jkiss_seed(struct rng_jkiss_state *state);
unsigned int rng_jkiss(struct rng_jkiss_state *state);


#ifdef __cplusplus
}
#endif

#endif /* _RNGS_H */
