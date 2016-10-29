/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

/* functions computing roots of polynomials. Return the number of real roots */

#include <stdio.h>
#include <math.h>

#include "polysolve.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* see http://mathworld.wolfram.com/QuadraticFormula.html */
int solve_deg2(double a, double b, double c, double *x0, double *x1)
{
double delta = b * b - 4 * a * c;
double inv_2a, sqrt_delta;

  if (delta < 0) return 0;

  if (a == 0){
    /* solve first order system */
    *x1 = 0;
    if (b != 0){
      *x0 = -c / b;
      return 1;
    }

    *x0 = 0;
    return 0;
  }

  inv_2a = 0.5 / a;

  if (delta == 0)
  {
    *x0 = -b * inv_2a;
    *x1 = *x0;
    return 1;
  }

  sqrt_delta = sqrt(delta);
  *x0 = (-b + sqrt_delta) * inv_2a;
  *x1 = (-b - sqrt_delta) * inv_2a;
  return 2;
}


/* see http://mathworld.wolfram.com/CubicEquation.html */
int solve_deg3(double a, double b, double c, double d, 
               double *x0, double *x1, double *x2)
{
double inv_a, b_a, b_a2, c_a, d_a;
double Q, R, Q3, D, b_a_3;
double AD, BD;

  if (a == 0)
  {
    /* solve second order system */
    if (b == 0)
    {
      /* solve first order system */
      if (c == 0) 
        return 0;

      *x0 = -d / c;
      return 1;
    }

    *x2 = 0;
    return solve_deg2(b, c, d, x0, x1);
  }

  /* calculate the normalized form x^3 + a2 * x^2 + a1 * x + a0 = 0 */
  inv_a = 1. / a;
  b_a = inv_a * b; b_a2 = b_a * b_a;
  c_a = inv_a * c;
  d_a = inv_a * d;

  /* solve the cubic equation */
  Q = (3 * c_a - b_a2) / 9;
  R = (9 * b_a * c_a - 27 * d_a - 2 * b_a * b_a2) / 54;
  Q3 = Q * Q * Q;
  D = Q3 + R * R;
  b_a_3 = (1. / 3.) * b_a;

  if (Q == 0) {
    if(R == 0) {
      *x0 = *x1 = *x2 = - b_a_3;
      return 3;
    }
    else
    {
      *x0 = pow(2 * R, 1 / 3.0) - b_a_3;
      return 1;
    }
  }

  if (D <= 0)
  {
    /* three real roots */
    double theta = acos(R / sqrt(-Q3));
    double sqrt_Q = sqrt(-Q);
    *x0 = 2 * sqrt_Q * cos(theta             / 3.0) - b_a_3;
    *x1 = 2 * sqrt_Q * cos((theta + 2 * M_PI)/ 3.0) - b_a_3;
    *x2 = 2 * sqrt_Q * cos((theta + 4 * M_PI)/ 3.0) - b_a_3;

    return 3;
  }

  /* D > 0, only one real root */
  AD = pow(fabs(R) + sqrt(D), 1.0 / 3.0) * (R > 0 ? 1 : (R < 0 ? -1 : 0));
  BD = (AD == 0) ? 0 : -Q / AD;

  /* calculate the sole real root */
  *x0 = AD + BD - b_a_3;

  return 1;
}

/* see http://mathworld.wolfram.com/QuarticEquation.html */
int solve_deg4(double a, double b, double c, double d, double e,
               double *x0, double *x1, double *x2, double *x3)
{
double inv_a, b2, bc, b3, b_4;
double r0, r1, r2;
int n, nb_real_roots;
double R2, R_2, R, inv_R;
double D2, E2;

  if (a == 0) 
  {
    *x3 = 0;
    return solve_deg3(b, c, d, e, x0, x1, x2);
  }

  /* normalize coefficients */
  inv_a = 1. / a;
  b *= inv_a; c *= inv_a; d *= inv_a; e *= inv_a;
  b2 = b * b; bc = b * c; b3 = b2 * b;

  /* solve resultant cubic */
  n = solve_deg3(1, -c, d * b - 4 * e, 4 * c * e - d * d - b2 * e, &r0, &r1, &r2);
  if (n == 0) return 0;

  /* calculate R^2 */
  R2 = 0.25 * b2 - c + r0;
  if (R2 < 0)
    return 0;

  R = sqrt(R2);
  inv_R = 1. / R;

  nb_real_roots = 0;

  /* calculate D^2 and E^2 */
  if (R < 10E-12)
  {
    double temp = r0 * r0 - 4 * e;
    if (temp < 0)
      D2 = E2 = -1;
    else
    {
      double sqrt_temp = sqrt(temp);
      D2 = 0.75 * b2 - 2 * c + 2 * sqrt_temp;
      E2 = D2 - 4 * sqrt_temp;
    }
  }
  else
  {
    double u = 0.75 * b2 - 2 * c - R2, v = 0.25 * inv_R * (4 * bc - 8 * d - b3);
    D2 = u + v;
    E2 = u - v;
  }

  b_4 = 0.25 * b; R_2 = 0.5 * R;
  if (D2 >= 0) {
    double D = sqrt(D2);
    double D_2 = 0.5 * D;
    nb_real_roots = 2;
    *x0 = R_2 + D_2 - b_4;
    *x1 = *x0 - D;
  }

  /* calculate E^2 */
  if (E2 >= 0) {
    double E = sqrt(E2);
    double E_2 = 0.5 * E;
    if (nb_real_roots == 0)
    {
      *x0 = - R_2 + E_2 - b_4;
      *x1 = *x0 - E;
      nb_real_roots = 2;
    }
    else
    {
      *x2 = - R_2 + E_2 - b_4;
      *x3 = *x2 - E;
      nb_real_roots = 4;
    }
  }

  return nb_real_roots;
}
