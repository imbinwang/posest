/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////
//
#ifndef _COMPILER_H
#define _COMPILER_H

/* C compiler specific constants & macros */

/* lapack f77 routines name mangling */
#define F77_FUNC(func)    func ## _

/* note: intel's icc defines both __ICC & __INTEL_COMPILER.
 * Also, some compilers other than gcc define __GNUC__,
 * therefore gcc should be checked last
 */
#ifdef _MSC_VER
#define inline __inline // MSVC
#elif !defined(__ICC) && !defined(__INTEL_COMPILER) && !defined(__GNUC__)
#define inline // other than MSVC, ICC, GCC: define empty
#endif

#ifdef _MSC_VER
#define restrict __restrict // MSVC
#elif !defined(__ICC) && !defined(__INTEL_COMPILER) && !defined(__GNUC__)
#define restrict // other than MSVC, ICC, GCC: define empty
#elif !defined(__INTEL_COMPILER) && defined(__GNUC__)
#define restrict __restrict__ // GCC: enable as GCC extension
//#else // ICC, leave intact
#endif

#ifdef _MSC_VER
#include <float.h>
#define POSEST_FINITE _finite // MSVC
#elif defined(__ICC) || defined(__INTEL_COMPILER) || defined(__GNUC__)
#define POSEST_FINITE finite // ICC, GCC
#else
#define POSEST_FINITE finite // other than MSVC, ICC, GCC, let's hope this will work
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923  /* pi/2 */
#endif

#ifdef _MSC_VER // MSVC: use the fsincos assembly instruction
// SINCOS(double x, double *snptr_, double *csptr_)
#define SINCOS(x, snptr_, csptr_)                             \
do {                                                          \
  double *snptr=snptr_, *csptr=csptr_;                        \
__asm {                                                       \
  __asm fld QWORD PTR [x]                                     \
  __asm fsincos                                               \
  __asm mov ebx,[csptr]       /* get the pointer into ebx */  \
  __asm fstp QWORD PTR [ebx]  /* store through pointer    */  \
  __asm mov ebx,[snptr]                                       \
  __asm fstp QWORD PTR [ebx]                                  \
  }                                                           \
} while (0)
#elif defined(__ICC) || defined(__INTEL_COMPILER) || (defined(__GNUC__) && !defined(__MINGW32__)) // ICC, GCC: use sincos
#define SINCOS(x, sn, cs) sincos(x, sn, cs)
#else
#define SINCOS(x, sn, cs) do { *(sn)=sin(x); *(cs)=cos(x); } while (0) // other than MSVC, ICC, GCC: use sin & cos
#endif

#ifdef HAVE_CBRT
#define CBRT(x) cbrt(x)
#else
#define CBRT(x) pow((double)(x), 0.333333)
#endif /* HAVE_CBRT */

#if defined(__ICC) || defined(__INTEL_COMPILER) || defined(__GNUC__)  // ICC, GCC
#define _flushall() fflush(NULL)
#endif

#endif /* _COMPILER_H */
