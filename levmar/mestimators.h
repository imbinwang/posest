////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Prototypes and definitions for the common M-estimators
//  It can be integrated into Lev-Mar least-squares as reweighted least-squares problem.
//  http://research.microsoft.com/en-us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html
//  Copyright (C) 2016  Bin Wang (binwangsdu at gmail dot com)
//  School of Computer Science and Technology, Shandong University,
//  Jinan, P. R. China.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _MESTIMATORS_H_
#define _MESTIMATORS_H_

/* specifies whether double precision routines will be compiled or not */
#define ME_DBL_PREC
/* specifies whether single precision routines will be compiled or not */
#define ME_SNGL_PREC

//#define ME_L2 0
//#define ME_L1 1
//#define ME_L1L2 2
#define ME_FAIR 3
//#define ME_HUBER 4
//#define ME_CAUCHY 5
#define ME_GEMANMCCLURE 6
//#define ME_WELSCH 7
#define ME_TUKEY 8

#ifdef __cplusplus
extern "C" {
#endif

#ifdef ME_DBL_PREC
	extern double dl2_cost(double x);
	extern double dl2_influence(double x);
	extern double dl2_weight(double x);

	extern double dl1_cost(double x);
	extern double dl1_influence(double x);
	extern double dl1_weight(double x);

	extern double dl1l2_cost(double x);
	extern double dl1l2_influence(double x);
	extern double dl1l2_weight(double x);

	extern double dfair_cost(double x, double c);
	extern double dfair_influence(double x, double c);
	extern double dfair_weight(double x, double c);

	extern double dhuber_cost(double x, double c);
	extern double dhuber_influence(double x, double c);
	extern double dhuber_weight(double x, double c);

	extern double dcauchy_cost(double x, double c);
	extern double dcauchy_influence(double x, double c);
	extern double dcauchy_weight(double x, double c);

	extern double dgemanmcclure_cost(double x, double c);
	extern double dgemanmcclure_influence(double x, double c);
	extern double dgemanmcclure_weight(double x, double c);

	extern double dwelsch_cost(double x, double c);
	extern double dwelsch_influence(double x, double c);
	extern double dwelsch_weight(double x, double c);

	extern double dtukey_cost(double x, double c);
	extern double dtukey_influence(double x, double c);
	extern double dtukey_weight(double x, double c);
#endif

#ifdef ME_SNGL_PREC
	extern float sl2_cost(float x);
	extern float sl2_influence(float x);
	extern float sl2_weight(float x);

	extern float sl1_cost(float x);
	extern float sl1_influence(float x);
	extern float sl1_weight(float x);

	extern float sl1l2_cost(float x);
	extern float sl1l2_influence(float x);
	extern float sl1l2_weight(float x);

	extern float sfair_cost(float x, float c);
	extern float sfair_influence(float x, float c);
	extern float sfair_weight(float x, float c);

	extern float shuber_cost(float x, float c);
	extern float shuber_influence(float x, float c);
	extern float shuber_weight(float x, float c);

	extern float scauchy_cost(float x, float c);
	extern float scauchy_influence(float x, float c);
	extern float scauchy_weight(float x, float c);

	extern float sgemanmcclure_cost(float x, float c);
	extern float sgemanmcclure_influence(float x, float c);
	extern float sgemanmcclure_weight(float x, float c);

	extern float swelsch_cost(float x, float c);
	extern float swelsch_influence(float x, float c);
	extern float swelsch_weight(float x, float c);

	extern float stukey_cost(float x, float c);
	extern float stukey_influence(float x, float c);
	extern float stukey_weight(float x, float c);
#endif

#ifdef __cplusplus
}
#endif


#endif //_MESTIMATORS_H_
