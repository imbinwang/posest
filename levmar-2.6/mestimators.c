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

/********************************************************************************
* The common M-estimators. The same core code is used with
* appropriate #defines to derive single and double precision versions, see
* also mestimators_core.c
********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "mestimators.h"
#include "compiler.h"

#define SIGN(x) (((x)>=0.0)? 1 : -1)

#if !defined(ME_DBL_PREC) && !defined(ME_SNGL_PREC)
#error At least one of ME_DBL_PREC, ME_SNGL_PREC should be defined!
#endif


#ifdef ME_SNGL_PREC
/* L2 */
float sl2_cost(float x)
{
	return 0.5*x*x;
}
float sl2_influence(float x)
{
	return x;
}
float sl2_weight(float x)
{
	return 1;
}

/* L1 */
float sl1_cost(float x)
{
	return abs(x);
}
float sl1_influence(float x)
{
	return SIGN(x);
}
float sl1_weight(float x)
{
	return 1 / SIGN(x);
}

/* L1L2 */
float sl1l2_cost(float x)
{
	return 2 * (sqrtf((1 + 0.5*x*x)) - 1);
}
float sl1l2_influence(float x)
{
	return x / sqrtf((1 + 0.5*x*x));
}
float sl1l2_weight(float x)
{
	return 1 / sqrtf((1 + 0.5*x*x));
}

/* Fair */
float sfair_cost(float x, float c)
{
	return c * c * (abs(x) / c - logf((1 + abs(x) / c)));
}
float sfair_influence(float x, float c)
{
	return x / (1 + abs(x) / c);
}
float sfair_weight(float x, float c)
{
	return 1 / (1 + abs(x) / c);
}

/* Huber */
float shuber_cost(float x, float c)
{
	if (abs(x) <= c)
	{
		return 0.5*x*x;
	}
	else
	{
		return c*(abs(x) - 0.5*c);
	}
}
float shuber_influence(float x, float c)
{
	if (abs(x) <= c)
	{
		return x;
	}
	else
	{
		return c*SIGN(x);
	}
}
float shuber_weight(float x, float c)
{
	if (abs(x) <= c)
	{
		return 1;
	}
	else
	{
		return c / abs(x);
	}
}

/* CAUCHY */
float scauchy_cost(float x, float c)
{
	return 0.5*c*c*logf((1 + (x / c)*(x / c)));
}
float scauchy_influence(float x, float c)
{
	return x / (1 + (x / c)*(x / c));
}
float scauchy_weight(float x, float c)
{
	return 1 / (1 + (x / c)*(x / c));
}

/* GEMANMCCLURE */
float sgemanmcclure_cost(float x, float c)
{
	return 0.5*x*x / (1 + x*x);
}
float sgemanmcclure_influence(float x, float c)
{
	return x / ((1 + x*x)*(1 + x*x));
}
float sgemanmcclure_weight(float x, float c)
{
	return 1 / ((1 + x*x)*(1 + x*x));
}

/* WELSCH */
float swelsch_cost(float x, float c)
{
	return 0.5*c*c*(1 - expf((-(x / c)*(x / c))));
}
float swelsch_influence(float x, float c)
{
	return x * expf((-(x / c)*(x / c)));
}
float swelsch_weight(float x, float c)
{
	return expf((-(x / c)*(x / c)));
}

/* TUKEY */
float stukey_cost(float x, float c)
{
	if (abs(x) <= c)
	{
		return (c*c / 6)*(1 - (1 - (x / c)*(x / c))*(1 - (x / c)*(x / c))*(1 - (x / c)*(x / c)));
	}
	else
	{
		return c*c / 6;
	}
}
float stukey_influence(float x, float c)
{
	if (abs(x) <= c)
	{
		return x*(1 - (x / c)*(x / c))*(1 - (x / c)*(x / c));
	}
	else
	{
		return 0;
	}
}
float stukey_weight(float x, float c)
{
	if (abs(x) <= c)
	{
		return (1 - (x / c)*(x / c))*(1 - (x / c)*(x / c));
	}
	else
	{
		return 0;
	}
}
#endif /* ME_SNGL_PREC */

#ifdef ME_DBL_PREC
/* L2 */
double dl2_cost(double x)
{
	return 0.5*x*x;
}
double dl2_influence(double x)
{
	return x;
}
double dl2_weight(double x)
{
	return 1;
}

/* L1 */
double dl1_cost(double x)
{
	return abs(x);
}
double dl1_influence(double x)
{
	return SIGN(x);
}
double dl1_weight(double x)
{
	return 1 / SIGN(x);
}

/* L1L2 */
double dl1l2_cost(double x)
{
	return 2 * (sqrt((1 + 0.5*x*x)) - 1);
}
double dl1l2_influence(double x)
{
	return x / sqrt((1 + 0.5*x*x));
}
double dl1l2_weight(double x)
{
	return 1 / sqrt((1 + 0.5*x*x));
}

/* Fair */
double dfair_cost(double x, double c)
{
	return c * c * (abs(x) / c - log((1 + abs(x) / c)));
}
double dfair_influence(double x, double c)
{
	return x / (1 + abs(x) / c);
}
double dfair_weight(double x, double c)
{
	return 1 / (1 + abs(x) / c);
}

/* Huber */
double dhuber_cost(double x, double c)
{
	if (abs(x) <= c)
	{
		return 0.5*x*x;
	}
	else
	{
		return c*(abs(x) - 0.5*c);
	}
}
double dhuber_influence(double x, double c)
{
	if (abs(x) <= c)
	{
		return x;
	}
	else
	{
		return c*SIGN(x);
	}
}
double dhuber_weight(double x, double c)
{
	if (abs(x) <= c)
	{
		return 1;
	}
	else
	{
		return c / abs(x);
	}
}

/* CAUCHY */
double dcauchy_cost(double x, double c)
{
	return 0.5*c*c*log((1 + (x / c)*(x / c)));
}
double dcauchy_influence(double x, double c)
{
	return x / (1 + (x / c)*(x / c));
}
double dcauchy_weight(double x, double c)
{
	return 1 / (1 + (x / c)*(x / c));
}

/* GEMANMCCLURE */
double dgemanmcclure_cost(double x, double c)
{
	return 0.5*x*x / (1 + x*x);
}
double dgemanmcclure_influence(double x, double c)
{
	return x / ((1 + x*x)*(1 + x*x));
}
double dgemanmcclure_weight(double x, double c)
{
	return 1 / ((1 + x*x)*(1 + x*x));
}

/* WELSCH */
double dwelsch_cost(double x, double c)
{
	return 0.5*c*c*(1 - exp((-(x / c)*(x / c))));
}
double dwelsch_influence(double x, double c)
{
	return x * exp((-(x / c)*(x / c)));
}
double dwelsch_weight(double x, double c)
{
	return exp((-(x / c)*(x / c)));
}

/* TUKEY */
double dtukey_cost(double x, double c)
{
	if (abs(x) <= c)
	{
		return (c*c / 6)*(1 - (1 - (x / c)*(x / c))*(1 - (x / c)*(x / c))*(1 - (x / c)*(x / c)));
	}
	else
	{
		return c*c / 6;
	}
}
double dtukey_influence(double x, double c)
{
	if (abs(x) <= c)
	{
		return x*(1 - (x / c)*(x / c))*(1 - (x / c)*(x / c));
	}
	else
	{
		return 0;
	}
}
double dtukey_weight(double x, double c)
{
	if (abs(x) <= c)
	{
		return (1 - (x / c)*(x / c))*(1 - (x / c)*(x / c));
	}
	else
	{
		return 0;
	}
}
#endif /* ME_DBL_PREC */
