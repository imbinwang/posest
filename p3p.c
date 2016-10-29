/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011-12  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

//#include <sam.h>

#include "compiler.h"

#include "polysolve.h"
#include "util.h"
#include "p3p.h"

/* set intrinsics provided in K */
void p3p_set_calib(struct p3p_calib_params *p, double K[9])
{
  p->fx = K[0];
  p->fy = K[4];
  p->cx = K[2];
  p->cy = K[5];
  p->s = K[1];

  p->inv_fx = 1. / p->fx;
  p->inv_fy = 1. / p->fy;
  p->cy_fy = p->cy / p->fy;
  p->s_fxfy = p->s / (p->fx*p->fy);
  p->scy_cxfy_fxfy = (p->s*p->cy - p->cx*p->fy)/(p->fx*p->fy);
}


/* Given the squared 3D distances between three points and cosines of 3 angles at the apex,
 * calculates the lengths of the line segments connecting projection center and the
 * three 3D points (A, B, C) by applying the direct algorithm of Grunert
 * [Grunert-1841], as reviewed in [Haralick-IJCV1994]
 */
static int lengths_Grunert(double lengths[4][3], double sqdistances[3], double cosines[3])
{
double aa, bb, cc; // squared lengths
double ca, cb, cg; // angles
double caca, cbcb, cgcg;
double q1, q2, q3, q4;
double A, B, C, D, E, real_roots[4];
int n, nb_solutions=0;
register int i;

  aa=sqdistances[0];
  bb=sqdistances[1];
  cc=sqdistances[2];

  ca=cosines[0];
  cb=cosines[1];
  cg=cosines[2];

  /* roots of fourth order polynomial, eq. (9) */
  caca=ca*ca; cbcb=cb*cb; cgcg=cg*cg;

  q1=(aa - cc)/bb;
  q2=(aa + cc)/bb;
  q3=(bb - cc)/bb;
  q4=(bb - aa)/bb;

  A=(q1-1.)*(q1-1.) - 4.*cc*caca/bb;
  B=4.*(q1*(1.-q1)*cb - (1.-q2)*ca*cg + 2.*cc*caca*cb/bb);
  C=2.*(q1*q1 - 1. + 2.*q1*q1*cbcb + 2.*q3*caca - 4.*q2*ca*cb*cg + 2.*q4*cgcg);
  D=4.*(-q1*(1.+q1)*cb + 2.*aa*cgcg*cb/bb - (1.-q2)*ca*cg);
  E=(1.+q1)*(1.+q1) - 4.*aa*cgcg/bb;

  n=solve_deg4(A, B, C, D, E, real_roots, real_roots+1, real_roots+2, real_roots+3);

  if(n==0) return 0;

  nb_solutions=0;
  /* examine each solution v */
  for(i=0; i<n; ++i) {
    double u, v=real_roots[i];
    double s1, s1sq;


    /* substitution for u in eq. (8),
     *              for s1 in eq. (5) and
     *              for s2, s3 in eq. (4)
     */
    u=((q1-1.)*v*v - 2.*q1*cb*v + 1. + q1) / (2.*(cg-v*ca));

    /* three alternative expressions for s1, attempt to find a finite one */
    s1sq=aa/(u*u+v*v-2.*u*v*ca); // s1_1
    if(!POSEST_FINITE(s1sq)){
      s1sq=bb/(1.+v*v-2.*v*cb); // s1_2
      if(!POSEST_FINITE(s1sq)){
        s1sq=cc/(1.+u*u-2.*u*cg); // s1_3
        if(!POSEST_FINITE(s1sq)) // all three are infinite, no solution
          continue;
      }
    }
    s1=sqrt(s1sq);

    lengths[nb_solutions][0]=s1;
    lengths[nb_solutions][1]=u*s1; // s2
    lengths[nb_solutions][2]=v*s1; // s3

    ++nb_solutions;
  }

  return nb_solutions;
}

/* Given the squared 3D distances between three points and the cosines of 3 angles at the apex,
 * calculates the lengths of the line segments connecting projection center (P) and the three
 * 3D points (A, B, C). Returned distances are for |PA|, |PB|, |PC| respectively.
 * Only the solution to the main branch.
 * Reference : X.S. Gao, X.-R. Hou, J. Tang, H.-F. Chang; "Complete Solution Classification for the Perspective-Three-Point Problem"
 * IEEE Trans. on PAMI, vol. 25, No. 8, August 2003
 * \param lengths3D Lengths of line segments up to four solutions.
 * \param dist3D Distance between 3D points in pairs |BC|, |AC|, |AB|.
 * \param cosines Cosine of the angles /_BPC, /_APC, /_APB.
 * \returns Number of solutions.
 * WARNING: NOT ALL THE DEGENERATE CASES ARE IMPLEMENTED
 */

static int solve_for_lengths(double lengths[4][3], double sqdistances[3], double cosines[3])
{
  register int i;
  double p, q, r, inv_d22, a, b;
  double a2, b2, p2, q2, r2, pr, pqr, ab;
  double a_2, A, a_4, B, C, D, E, temp, b0;
  double real_roots[4];
  int n, nb_solutions;
  double r3, pr2, r3q, inv_b0;

  p = cosines[0] * 2;
  q = cosines[1] * 2;
  r = cosines[2] * 2;

  inv_d22 = 1. / sqdistances[2];
  a = inv_d22 * sqdistances[0];
  b = inv_d22 * sqdistances[1];

  a2 = a * a; b2 = b * b; p2 = p * p; q2 = q * q; r2 = r * r;
  pr = p * r; pqr = q * pr;

  /* check reality condition (the four points should not be coplanar) */
  if (p2 + q2 + r2 - pqr - 1 == 0)
    return 0;

  ab = a * b; a_2 = 2*a;

  A = -2 * b + b2 + a2 + 1 + ab*(2 - r2) - a_2;

  /* check reality condition */
  if (A == 0) return 0;

  a_4 = 4*a;

  B = q*(-2*(ab + a2 + 1 - b) + r2*ab + a_4) + pr*(b - b2 + ab);
  C = q2 + b2*(r2 + p2 - 2) - b*(p2 + pqr) - ab*(r2 + pqr) + (a2 - a_2)*(2 + q2) + 2;
  D = pr*(ab-b2+b) + q*((p2-2)*b + 2 * (ab - a2) + a_4 - 2);
  E = 1 + 2*(b - a - ab) + b2 - b*p2 + a2;

  temp = (p2*(a-1+b) + r2*(a-1-b) + pqr - a*pqr);
  b0 = b * temp * temp;
  /* check reality condition */
  if (b0 == 0)
    return 0;

  n = solve_deg4(A, B, C, D, E,  real_roots, real_roots+1, real_roots+2, real_roots+3);

  if (n == 0)
    return 0;

  nb_solutions = 0;
  r3 = r2*r; pr2 = p*r2; r3q = r3 * q;
  inv_b0 = 1. / b0;

  /* for each solution of x */
  for(i = 0; i < n; i++) {
    double x = real_roots[i];
    double x2, b1, y, v;
    double X, Y, Z;

    /* check reality condition */
    if (x <= 0)
      continue;

    x2 = x*x;

    b1 = 
      ((1-a-b)*x2 + (q*a-q)*x + 1 - a + b) * 
      (((r3*(a2 + ab*(2 - r2) - a_2 + b2 - 2*b + 1)) * x + 

      (r3q*(2*(b-a2) + a_4 + ab*(r2 - 2) - 2) + pr2*(1 + a2 + 2*(ab-a-b) + r2*(b - b2) + b2))) * x2 +

      (r3*(q2*(1-2*a+a2) + r2*(b2-ab) - a_4 + 2*(a2 - b2) + 2) + r*p2*(b2 + 2*(ab - b - a) + 1 + a2) + pr2*q*(a_4 + 2*(b - ab - a2) - 2 - r2*b)) * x + 

      2*r3q*(a_2 - b - a2 + ab - 1) + pr2*(q2 - a_4 + 2*(a2 - b2) + r2*b + q2*(a2 - a_2) + 2) + 
      p2*(p*(2*(ab - a - b) + a2 + b2 + 1) + 2*q*r*(b + a_2 - a2 - ab - 1)));

    /* check reality condition */
    if (b1 <= 0)
      continue;

    y = inv_b0 * b1;
    v = x2 + y*y - x*y*r;

    if (v <= 0)
      continue;

    Z = sqrt(sqdistances[2] / v);
    X = x * Z;
    Y = y * Z;

    lengths[nb_solutions][0] = X;
    lengths[nb_solutions][1] = Y;
    lengths[nb_solutions][2] = Z;

    nb_solutions++;
  }

  return nb_solutions;
}


/* disambiguate P3P using a 4th point and 2D reprojection errors.
 *
 * Returns 0 on failure, nonzero otherwise
 */
int p3p_solve4_2Derr(struct p3p_calib_params *cp,
           double m[4][2], double M[4][3],
           double R[3][3], double t[3])
{
  register int i, j;
  double Rs[4][3][3], ts[4][3];
  int n, ns;
  double mu3, mv3, X3, Y3, Z3;
  double min_reproj;

  n = p3p_solve3(cp, m, M, Rs, ts); // solve using first 3 points

  if (n == 0) return 0;

  mu3=m[3][0]; mv3=m[3][1];
  X3=M[3][0]; Y3=M[3][1]; Z3=M[3][2];

  ns = 0;
  min_reproj = DBL_MAX;
  for(i = 0; i < n; i++)
  {
    double X3p, Y3p, Z3p, mu3p, mv3p, reproj;

#if 0
    printf("P3P solution #%d\n", i);
    printf("%g %g %g\n", Rs[i][0][0], Rs[i][0][1], Rs[i][0][2]);
    printf("%g %g %g\n", Rs[i][1][0], Rs[i][1][1], Rs[i][1][2]);
    printf("%g %g %g\n", Rs[i][2][0], Rs[i][2][1], Rs[i][2][2]);
    printf("%g %g %g\n", ts[i][0], ts[i][1], ts[i][2]);
#endif

    X3p = Rs[i][0][0] * X3 + Rs[i][0][1] * Y3 + Rs[i][0][2] * Z3 + ts[i][0];
    Y3p = Rs[i][1][0] * X3 + Rs[i][1][1] * Y3 + Rs[i][1][2] * Z3 + ts[i][1];
    Z3p = Rs[i][2][0] * X3 + Rs[i][2][1] * Y3 + Rs[i][2][2] * Z3 + ts[i][2];
    mv3p = Y3p / Z3p;
    mu3p = cp->cx + cp->fx * X3p / Z3p + cp->s * mv3p;
    mv3p = cp->cy + cp->fy * mv3p;
    reproj = (mu3p - mu3) * (mu3p - mu3) + (mv3p - mv3) * (mv3p - mv3);
    if (min_reproj > reproj)
    {
      ns = i;
      min_reproj = reproj;
    }
  }

  for(i = 0; i < 3; i++)
  {
    for(j = 0; j < 3; j++)
      R[i][j] = Rs[ns][i][j];
    t[i] = ts[ns][i];
  }

  return 1+ns;
}

/* disambiguate P3P using a 4th point and 3D errors 
 *
 * Returns 0 on failure, nonzero otherwise
 */
int p3p_solve4_3Derr(struct p3p_calib_params *cp,
           double m[4][2], double M[4][3], double plnorm[3],
           double R[3][3], double t[3])
{
  register int i, j;
  double Rs[4][3][3], ts[4][3];
  int n, ns;
  double mu3, mv3, X0, Y0, Z0, X3, Y3, Z3, rp[3];
  double min_dist;

  /* Note: plnorm could be computed as cross(M1-M0, M2-M0);
   * it is passed as an argument to save some calculations
   */

  n = p3p_solve3(cp, m, M, Rs, ts); // solve using first 3 points

  if (n == 0) return 0;

  mu3=m[3][0]; mv3=m[3][1];
  X0=M[0][0]; Y0=M[0][1]; Z0=M[0][2];
  X3=M[3][0]; Y3=M[3][1]; Z3=M[3][2];

  /* normalized image plane point corresponding to (mu3, mv3) */
  rp[0]=cp->inv_fx * mu3 - cp->s_fxfy * mv3 + cp->scy_cxfy_fxfy;
  rp[1]=cp->inv_fy * mv3 - cp->cy_fy;
  rp[2]=1.0;

  ns = 0;
  min_dist = DBL_MAX;
  for(i = 0; i < n; i++)
  {
    double nn[3], X0p, Y0p, Z0p, X3p, Y3p, Z3p, X3i, Y3i, Z3i, u, dist;

#if 0
    printf("P3P solution #%d\n", i);
    printf("%g %g %g\n", Rs[i][0][0], Rs[i][0][1], Rs[i][0][2]);
    printf("%g %g %g\n", Rs[i][1][0], Rs[i][1][1], Rs[i][1][2]);
    printf("%g %g %g\n", Rs[i][2][0], Rs[i][2][1], Rs[i][2][2]);
    printf("%g %g %g\n", ts[i][0], ts[i][1], ts[i][2]);
#endif

    /* transform the equation of the plane defined by M0, M1, M2 with Rs[i], ts[i] */
    /* nn=R*n */
    nn[0] = Rs[i][0][0]*plnorm[0] + Rs[i][0][1]*plnorm[1] + Rs[i][0][2]*plnorm[2];
    nn[1] = Rs[i][1][0]*plnorm[0] + Rs[i][1][1]*plnorm[1] + Rs[i][1][2]*plnorm[2];
    nn[2] = Rs[i][2][0]*plnorm[0] + Rs[i][2][1]*plnorm[1] + Rs[i][2][2]*plnorm[2];

    /* transform M[0] according to Rs[i], ts[i] */
    X0p = Rs[i][0][0]*X0 + Rs[i][0][1]*Y0 + Rs[i][0][2]*Z0 + ts[i][0];
    Y0p = Rs[i][1][0]*X0 + Rs[i][1][1]*Y0 + Rs[i][1][2]*Z0 + ts[i][1];
    Z0p = Rs[i][2][0]*X0 + Rs[i][2][1]*Y0 + Rs[i][2][2]*Z0 + ts[i][2];

    /* compute the point of intersection p of the ray through rp
     * with the plane having normal nn and going through [X0p, Y0p, Z0p]
     *
     * p=u*rp with u=(nn'*X0p)/(nn'*rp)
     */
    u=(nn[0]*X0p + nn[1]*Y0p + nn[2]*Z0p)/(nn[0]*rp[0] + nn[1]*rp[1] + nn[2]*rp[2]);
    X3i=u*rp[0]; Y3i=u*rp[1]; Z3i=u*rp[2];

    /* transform M[3] according to Rs[i], ts[i] */
    X3p = Rs[i][0][0]*X3 + Rs[i][0][1]*Y3 + Rs[i][0][2]*Z3 + ts[i][0];
    Y3p = Rs[i][1][0]*X3 + Rs[i][1][1]*Y3 + Rs[i][1][2]*Z3 + ts[i][1];
    Z3p = Rs[i][2][0]*X3 + Rs[i][2][1]*Y3 + Rs[i][2][2]*Z3 + ts[i][2];

    dist=(X3p-X3i)*(X3p-X3i) + (Y3p-Y3i)*(Y3p-Y3i) + (Z3p-Z3i)*(Z3p-Z3i);
    if (min_dist > dist)
    {
      ns = i;
      min_dist = dist;
    }
  }

  for(i = 0; i < 3; i++)
  {
    for(j = 0; j < 3; j++)
      R[i][j] = Rs[ns][i][j];
    t[i] = ts[ns][i];
  }

  return 1+ns;
}

int p3p_solve3(struct p3p_calib_params *cp,
           double m[3][2], double M[3][3],
           double R[4][3][3], double t[4][3])
{
  register int i;
  double mk0, mk1, mk2;
  double norm;
  double sqdistances[3];
  double cosines[3]; 
  double lengths[4][3];
  int n, nb_solutions;
  double mu0, mv0, mu1, mv1, mu2, mv2;
  double X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2;

  mu0=m[0][0]; mv0=m[0][1];
  mu1=m[1][0]; mv1=m[1][1];
  mu2=m[2][0]; mv2=m[2][1];
  X0=M[0][0]; Y0=M[0][1]; Z0=M[0][2];
  X1=M[1][0]; Y1=M[1][1]; Z1=M[1][2];
  X2=M[2][0]; Y2=M[2][1]; Z2=M[2][2];

  mu0 = cp->inv_fx * mu0 - cp->s_fxfy * mv0 + cp->scy_cxfy_fxfy;
  mv0 = cp->inv_fy * mv0 - cp->cy_fy;
  norm = sqrt(mu0 * mu0 + mv0 * mv0 + 1);
  mk0 = 1. / norm; mu0 *= mk0; mv0 *= mk0; 

  mu1 = cp->inv_fx * mu1 - cp->s_fxfy * mv1 + cp->scy_cxfy_fxfy;
  mv1 = cp->inv_fy * mv1 - cp->cy_fy;
  norm = sqrt(mu1 * mu1 + mv1 * mv1 + 1);
  mk1 = 1. / norm; mu1 *= mk1; mv1 *= mk1;

  mu2 = cp->inv_fx * mu2 - cp->s_fxfy * mv2 + cp->scy_cxfy_fxfy;
  mv2 = cp->inv_fy * mv2 - cp->cy_fy;
  norm = sqrt(mu2 * mu2 + mv2 * mv2 + 1);
  mk2 = 1. / norm; mu2 *= mk2; mv2 *= mk2; 

  sqdistances[0] = ( (X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2) + (Z1 - Z2) * (Z1 - Z2) );
  sqdistances[1] = ( (X0 - X2) * (X0 - X2) + (Y0 - Y2) * (Y0 - Y2) + (Z0 - Z2) * (Z0 - Z2) );
  sqdistances[2] = ( (X0 - X1) * (X0 - X1) + (Y0 - Y1) * (Y0 - Y1) + (Z0 - Z1) * (Z0 - Z1) );

  /* calculate angles */
  cosines[0] = mu1 * mu2 + mv1 * mv2 + mk1 * mk2;
  cosines[1] = mu0 * mu2 + mv0 * mv2 + mk0 * mk2;
  cosines[2] = mu0 * mu1 + mv0 * mv1 + mk0 * mk1;

  //n = solve_for_lengths(lengths, sqdistances, cosines);
  n = lengths_Grunert(lengths, sqdistances, cosines);

  nb_solutions = 0;
  for(i = 0; i < n; i++) 
  {
    double M_orig[3][3];
    register double len;

    len=lengths[i][0];
    M_orig[0][0] = len * mu0;
    M_orig[0][1] = len * mv0;
    M_orig[0][2] = len * mk0;

    len=lengths[i][1];
    M_orig[1][0] = len * mu1;
    M_orig[1][1] = len * mv1;
    M_orig[1][2] = len * mk1;

    len=lengths[i][2];
    M_orig[2][0] = len * mu2;
    M_orig[2][1] = len * mv2;
    M_orig[2][2] = len * mk2;

    if(posest_align3Pts(M_orig, M, R[nb_solutions], t[nb_solutions])) continue;
    //if(posest_alignNPts(M, M_orig, 3, (double *)R[nb_solutions], t[nb_solutions])) continue;
    //if(sam_absorq(M, M_orig, NULL, 3, (double *)R[nb_solutions], t[nb_solutions], NULL)) continue;

    nb_solutions++;
  }

  return nb_solutions;
}


/******************************** Kneip's P3P start ********************************/

/* Converted to C from P3p.cpp by Manolis Lourakis, May 2012 */


/* see http://www.asl.ethz.ch/people/kneipl/personal/p3p_code_final.zip */

/*
 * Copyright (c) 2011, Laurent Kneip, ETH Zurich
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ETH Zurich nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * P3p.cpp
 *
 *  Created on: Nov 2, 2010
 *      Author: Laurent Kneip
 * Description: Compute the absolute pose of a camera using three 3D-to-2D correspondences
 *   Reference: A Novel Parametrization of the P3P-Problem for a Direct Computation of
 *              Absolute Camera Position and Orientation
 *
 *       Input: imgPoints: 3x2 matrix with NON-UNITARY 2D image points (each column is a vector)
 *              worldPoints: 3x3 matrix with corresponding 3D world points (each column is a point)
 *              Rs, ts: matrices that will contain the solutions
 *      Output: int: 0 if correct execution
 *                  -1 if world points aligned
 */


#define _CROSSPROD(v, x, y){ (v)[0]=(x)[1]*(y)[2] - (x)[2]*(y)[1]; (v)[1]=(x)[2]*(y)[0] - (x)[0]*(y)[2]; (v)[2]=(x)[0]*(y)[1] - (x)[1]*(y)[0]; }


#if 0
/* C=A*B, for 3x3 matrices */
static void mat3x3Mult(double *C, double *A, double *B)
{
register int i, i3;

  for(i=0; i<3; ++i){
    i3=i*3;
    C[i3+0]=A[i3+0]*B[0*3+0] + A[i3+1]*B[1*3+0] + A[i3+2]*B[2*3+0];
    C[i3+1]=A[i3+0]*B[0*3+1] + A[i3+1]*B[1*3+1] + A[i3+2]*B[2*3+1];
    C[i3+2]=A[i3+0]*B[0*3+2] + A[i3+1]*B[1*3+2] + A[i3+2]*B[2*3+2];
  }
}
#endif

int p3p_Kneip(struct p3p_calib_params *cp, double imgPoints[3][2], double worldPoints[3][3], double Rs[4][3][3], double ts[4][3])
{
double *P1, *P2, *P3;
double *f1, *f2, *f3;
double temp1[3], temp2[3], Tf3[3];
double T[9], N[9], factors[5];
double *const e1=T, *const e2=T+3, *const e3=T+6; // T=[e1; e2; e3]
double *const n1=N, *const n2=N+3, *const n3=N+6; // N=[n1; n2; n3]
double norm1, roots[4];
double d_12, f_1, f_2, p_1, p_2, cos_beta, b;
double f_1_pw2, f_2_pw2, p_1_pw2, p_1_pw3, p_1_pw4, p_2_pw2, p_2_pw3, p_2_pw4, d_12_pw2, b_pw2;
register int i, j;
int nposes;
double nfeatureVectors[3][3];
double mu, mv, norm;

 // Normalization of image points

  mu=imgPoints[0][0]; mv=imgPoints[0][1];
  mu = cp->inv_fx * mu - cp->s_fxfy * mv + cp->scy_cxfy_fxfy;
  mv = cp->inv_fy * mv - cp->cy_fy;
  norm = sqrt(mu * mu + mv * mv + 1);
  nfeatureVectors[0][2] = 1. / norm; nfeatureVectors[0][0] = mu*nfeatureVectors[0][2]; nfeatureVectors[0][1] = mv*nfeatureVectors[0][2];

  mu=imgPoints[1][0]; mv=imgPoints[1][1];
  mu = cp->inv_fx * mu - cp->s_fxfy * mv + cp->scy_cxfy_fxfy;
  mv = cp->inv_fy * mv - cp->cy_fy;
  norm = sqrt(mu * mu + mv * mv + 1);
  nfeatureVectors[1][2] = 1. / norm; nfeatureVectors[1][0] = mu*nfeatureVectors[1][2]; nfeatureVectors[1][1] = mv*nfeatureVectors[1][2];

  mu=imgPoints[2][0]; mv=imgPoints[2][1];
  mu = cp->inv_fx * mu - cp->s_fxfy * mv + cp->scy_cxfy_fxfy;
  mv = cp->inv_fy * mv - cp->cy_fy;
  norm = sqrt(mu * mu + mv * mv + 1);
  nfeatureVectors[2][2] = 1. / norm; nfeatureVectors[2][0] = mu*nfeatureVectors[2][2]; nfeatureVectors[2][1] = mv*nfeatureVectors[2][2];

	// Extraction of world points

  P1=worldPoints[0];
  P2=worldPoints[1];
  P3=worldPoints[2];

#if 0
	// Verification that world points are not collinear

  temp1[0]=P2[0] - P1[0]; temp1[1]=P2[1] - P1[1]; temp1[2]=P2[2] - P1[2];
  temp2[0]=P3[0] - P1[0]; temp2[1]=P3[1] - P1[1]; temp2[2]=P3[2] - P1[2];
  _CROSSPROD(e2, temp1, temp2);
  if(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2] < 1E-15) return 0;
#endif

	// Extraction of feature vectors

  f1=nfeatureVectors[0];
  f2=nfeatureVectors[1];
  f3=nfeatureVectors[2];

	// Creation of intermediate camera frame

  e1[0]=f1[0]; e1[1]=f1[1]; e1[2]=f1[2];
  _CROSSPROD(e3, f1, f2);
  norm1=1.0/sqrt(e3[0]*e3[0] + e3[1]*e3[1] + e3[2]*e3[2]);
  e3[0]*=norm1; e3[1]*=norm1; e3[2]*=norm1;
  _CROSSPROD(e2, e3, e1);

#if 0
  T[0]=e1[0]; T[1]=e1[1]; T[2]=e1[2];
  T[3]=e2[0]; T[4]=e2[1]; T[5]=e2[2];
  T[6]=e3[0]; T[7]=e3[1]; T[8]=e3[2];
#endif

  /* temp1=T*f3 */
  temp1[0]=T[0]*f3[0] + T[1]*f3[1] + T[2]*f3[2];
  temp1[1]=T[3]*f3[0] + T[4]*f3[1] + T[5]*f3[2];
  temp1[2]=T[6]*f3[0] + T[7]*f3[1] + T[8]*f3[2];

	// Reinforce that temp1[2] > 0 for having theta in [0;pi]
	if(temp1[2]>0.0){
    f1=nfeatureVectors[1];
    f2=nfeatureVectors[0];
    f3=nfeatureVectors[2];

    e1[0]=f1[0]; e1[1]=f1[1]; e1[2]=f1[2];
    _CROSSPROD(e3, f1, f2);
    norm1=1.0/sqrt(e3[0]*e3[0] + e3[1]*e3[1] + e3[2]*e3[2]);
    e3[0]*=norm1; e3[1]*=norm1; e3[2]*=norm1;
    _CROSSPROD(e2, e3, e1);

#if 0
    T[0]=e1[0]; T[1]=e1[1]; T[2]=e1[2];
    T[3]=e2[0]; T[4]=e2[1]; T[5]=e2[2];
    T[6]=e3[0]; T[7]=e3[1]; T[8]=e3[2];
#endif

    /* temp1=T*f3 */
    temp1[0]=T[0]*f3[0] + T[1]*f3[1] + T[2]*f3[2];
    temp1[1]=T[3]*f3[0] + T[4]*f3[1] + T[5]*f3[2];
    temp1[2]=T[6]*f3[0] + T[7]*f3[1] + T[8]*f3[2];

    P1=worldPoints[1];
    P2=worldPoints[0];
    P3=worldPoints[2];
	}
  Tf3[0]=temp1[0]; Tf3[1]=temp1[1]; Tf3[2]=temp1[2];

	// Creation of intermediate world frame

  n1[0]=P2[0]-P1[0]; n1[1]=P2[1]-P1[1]; n1[2]=P2[2]-P1[2];
  norm1=1.0/sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
  n1[0]*=norm1; n1[1]*=norm1; n1[2]*=norm1;

  temp1[0]=P3[0]-P1[0]; temp1[1]=P3[1]-P1[1]; temp1[2]=P3[2]-P1[2];
  _CROSSPROD(n3, n1, temp1);
  norm1=1.0/sqrt(n3[0]*n3[0] + n3[1]*n3[1] + n3[2]*n3[2]);
  n3[0]*=norm1; n3[1]*=norm1; n3[2]*=norm1;

  _CROSSPROD(n2, n3, n1);

#if 0
  N[0]=n1[0]; N[1]=n1[1]; N[2]=n1[2];
  N[3]=n2[0]; N[4]=n2[1]; N[5]=n2[2];
  N[6]=n3[0]; N[7]=n3[1]; N[8]=n3[2];
#endif

	// Extraction of known parameters

  //temp1[0]=P3[0]-P1[0]; temp1[1]=P3[1]-P1[1]; temp1[2]=P3[2]-P1[2];
  /* N*(P3-P1) */
  temp2[0]=N[0]*temp1[0] + N[1]*temp1[1] + N[2]*temp1[2];
  temp2[1]=N[3]*temp1[0] + N[4]*temp1[1] + N[5]*temp1[2];
  temp2[2]=N[6]*temp1[0] + N[7]*temp1[1] + N[8]*temp1[2];
  P3[0]=temp2[0]; P3[1]=temp2[1]; P3[2]=temp2[2];
  
  temp1[0]=P2[0]-P1[0]; temp1[1]=P2[1]-P1[1]; temp1[2]=P2[2]-P1[2];
	d_12=sqrt(temp1[0]*temp1[0] + temp1[1]*temp1[1] + temp1[2]*temp1[2]); 
	f_1=Tf3[0]/Tf3[2];
	f_2=Tf3[1]/Tf3[2];
	p_1=P3[0];
	p_2=P3[1];

	cos_beta=f1[0]*f2[0] + f1[1]*f2[1] + f1[2]*f2[2];
	b=1./(1.-cos_beta*cos_beta) - 1.;

  b=(cos_beta<0.0)? -sqrt(b) : sqrt(b);

	// Definition of temporary variables for avoiding multiple computation

	f_1_pw2=f_1*f_1;
	f_2_pw2=f_2*f_2;
	p_1_pw2=p_1*p_1;
	p_1_pw3=p_1_pw2*p_1;
	p_1_pw4=p_1_pw3*p_1;
	p_2_pw2=p_2*p_2;
	p_2_pw3=p_2_pw2*p_2;
	p_2_pw4=p_2_pw3*p_2;
	d_12_pw2=d_12*d_12;
	b_pw2=b*b;

	// Computation of factors of 4th degree polynomial

	factors[0] = -f_2_pw2*p_2_pw4
				 -p_2_pw4*f_1_pw2
				 -p_2_pw4;

	factors[1] = 2*p_2_pw3*d_12*b
				 +2*f_2_pw2*p_2_pw3*d_12*b
				 -2*f_2*p_2_pw3*f_1*d_12;

	factors[2] = -f_2_pw2*p_2_pw2*p_1_pw2
			     -f_2_pw2*p_2_pw2*d_12_pw2*b_pw2
				 -f_2_pw2*p_2_pw2*d_12_pw2
				 +f_2_pw2*p_2_pw4
				 +p_2_pw4*f_1_pw2
				 +2*p_1*p_2_pw2*d_12
				 +2*f_1*f_2*p_1*p_2_pw2*d_12*b
				 -p_2_pw2*p_1_pw2*f_1_pw2
				 +2*p_1*p_2_pw2*f_2_pw2*d_12
				 -p_2_pw2*d_12_pw2*b_pw2
				 -2*p_1_pw2*p_2_pw2;

	factors[3] = 2*p_1_pw2*p_2*d_12*b
				 +2*f_2*p_2_pw3*f_1*d_12
				 -2*f_2_pw2*p_2_pw3*d_12*b
				 -2*p_1*p_2*d_12_pw2*b;

	factors[4] = -2*f_2*p_2_pw2*f_1*p_1*d_12*b
				 +f_2_pw2*p_2_pw2*d_12_pw2
				 +2*p_1_pw3*d_12
				 -p_1_pw2*d_12_pw2
				 +f_2_pw2*p_2_pw2*p_1_pw2
				 -p_1_pw4
				 -2*f_2_pw2*p_2_pw2*p_1*d_12
				 +p_2_pw2*f_1_pw2*p_1_pw2
				 +f_2_pw2*p_2_pw2*d_12_pw2*b_pw2;

	// Computation of roots

  nposes=solve_deg4(factors[0], factors[1], factors[2], factors[3], factors[4], roots, roots+1, roots+2, roots+3);
  //printf("%d quartic roots: %g %g %g %g\n", n, roots[0], roots[1], roots[2], roots[3]);

	// Backsubstitution of each solution

	for(i=0; i<nposes; ++i){
		double cot_alpha = (-f_1*p_1/f_2-roots[i]*p_2+d_12*b)/(-f_1*roots[i]*p_2/f_2+p_1-d_12);

		double cos_theta = roots[i];
		double sin_theta = sqrt(1.-roots[i]*roots[i]);
		double sin_alpha = sqrt(1./(cot_alpha*cot_alpha+1.));
		double cos_alpha = sqrt(1.-sin_alpha*sin_alpha);

    double *R, *t, RN[9];

		if(cot_alpha<0.0)
			cos_alpha = -cos_alpha;

    /* NOTE: Kneip's original code computes Rc,tc aligning the camera to the world frame.
     * Below I compute Rw,tw that align the world to the camera frame, ie Rw=Rc', tw=-Rc'*tc
     */

    R=(double *)Rs[i];
		R[0]=-cos_alpha; R[1]=-sin_alpha*cos_theta; R[2]=-sin_alpha*sin_theta;
		R[3]=sin_alpha;  R[4]=-cos_alpha*cos_theta; R[5]=-cos_alpha*sin_theta;
		R[6]=0.0;        R[7]=-sin_theta;           R[8]=cos_theta;

    /* R = (N'*R'*T)' = T'*R*N */
    /* RN=R*N */
    for(j=0; j<3; ++j){
      int j3=j*3;
      RN[j3+0]=R[j3+0]*N[0*3+0] + R[j3+1]*N[1*3+0] + R[j3+2]*N[2*3+0];
      RN[j3+1]=R[j3+0]*N[0*3+1] + R[j3+1]*N[1*3+1] + R[j3+2]*N[2*3+1];
      RN[j3+2]=R[j3+0]*N[0*3+2] + R[j3+1]*N[1*3+2] + R[j3+2]*N[2*3+2];
    }

    /* R = T'*RN */
    for(j=0; j<3; ++j){
      int j3=j*3;
      R[j3+0]=T[0*3+j]*RN[0*3+0] + T[1*3+j]*RN[1*3+0] + T[2*3+j]*RN[2*3+0];
      R[j3+1]=T[0*3+j]*RN[0*3+1] + T[1*3+j]*RN[1*3+1] + T[2*3+j]*RN[2*3+1];
      R[j3+2]=T[0*3+j]*RN[0*3+2] + T[1*3+j]*RN[1*3+2] + T[2*3+j]*RN[2*3+2];
    }

    /* temp1 = C */
    temp1[0]=d_12*cos_alpha*(sin_alpha*b+cos_alpha);
    temp1[1]=cos_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha);
    temp1[2]=sin_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha);

    /* C = -R*(P1 + N'*C) */    
    temp2[0]=P1[0] + N[0]*temp1[0] + N[3]*temp1[1] + N[6]*temp1[2];
    temp2[1]=P1[1] + N[1]*temp1[0] + N[4]*temp1[1] + N[7]*temp1[2];
    temp2[2]=P1[2] + N[2]*temp1[0] + N[5]*temp1[1] + N[8]*temp1[2];

    t=ts[i];
    t[0]=-(R[0]*temp2[0] + R[1]*temp2[1] + R[2]*temp2[2]);
    t[1]=-(R[3]*temp2[0] + R[4]*temp2[1] + R[5]*temp2[2]);
    t[2]=-(R[6]*temp2[0] + R[7]*temp2[1] + R[8]*temp2[2]);
	}

	return nposes;
}

#if 0
// test main for Kneip's P3P
int main()
{
int i, nsol;
double R[4][3][3], t[4][3];
struct p3p_calib_params cal;

/*
  double K[9]={443.59549, 0, 344.89962, 0, 444.751606, 207.652054, 0, 0, 1};
  double m[5][2]={
    {345.4272, 328.7803},
    {368.4228, 325.2764},
    {407.2719, 323.1395},
    {397.0703, 323.8156},
    {386.1493, 325.5412},
  };

  double M[5][3]={
    {0.3790, 0.8374, 66.2935},
    {0.4744, -0.0844, 65.4697},
    {0.6153, -1.8516, 64.3976},
    {0.4823, -1.3676, 64.7460},
    {0.4998, -1.0749, 65.0080},
  };
  */

  double K[9]={17647, 0, 1296, 0, 17647, 972, 0, 0, 1};
  double m[4][2]={
    {1147.56, 268.815},
    {1174.91, 876.597},
    {1226.25, 473.397},
    {1141.05, 1035.68},
  };

  double M[4][3]={
    {-0.03296, -0.060915, 3.5076},
    {0.002355, 0.026614, 3.6041},
    {0.012907, -0.053868, 3.5422},
    {-0.003344, 0.060906, 3.6277}
  };

  p3p_set_calib(&cal, K);

  nsol=p3p_Kneip(&cal, m, M, R, t); // only first 3 rows of n, M used here
  printf("P3P: %d\n", nsol);
  for(i=0; i<nsol; ++i){
    printf("\nSol %d\n", i);
    printf("%g %g %g\n", R[i][0][0], R[i][0][1], R[i][0][2]);
    printf("%g %g %g\n", R[i][1][0], R[i][1][1], R[i][1][2]);
    printf("%g %g %g\n", R[i][2][0], R[i][2][1], R[i][2][2]);

    printf("\n%g %g %g\n", t[i][0], t[i][1], t[i][2]);
  }

  exit(0);
}
#endif
/********************************** Kneip's P3P end **********************************/

#if 0 
// test main for Grunert's P3P/P4P
main()
{
int i, nsol;
double R[3][3], t[3];
double Rs[4][3][3], ts[4][3];
struct p3p_calib_params cal;

  double K[9]={17647, 0, 1296, 0, 17647, 972, 0, 0, 1};
  double m[4][2]={
    {1147.56, 268.815},
    {1174.91, 876.597},
    {1226.25, 473.397},
    {1141.05, 1035.68},
  };

  double M[4][3]={
    {-0.03296, -0.060915, 3.5076},
    {0.002355, 0.026614, 3.6041},
    {0.012907, -0.053868, 3.5422},
    {-0.003344, 0.060906, 3.6277}
  };

/*
  double K[9]={ 949.358704, 0.0, 683.450378, 0.0, 962.721877, 496.401459, 0.0, 0.0, 1.0};
  double m[4][2]={
    {828.815, 602.923},
    {567.755, 646.787},
    {794.484, 583.056},
    {580.417, 587.523},
  };

  double M[4][3]={
    {1350, 1350, 0},
    {270, 1350, 0}, 
    {1350, 1170, 0},
    {540, 990, 0}
  };
*/

  p3p_set_calib(&cal, K);

  nsol=p3p_solve3(&cal, m, M, Rs, ts);
  printf("P3P: %d\n", nsol);
  for(i=0; i<nsol; ++i){
    printf("\nSol %d\n", i);
    printf("%g %g %g\n", Rs[i][0][0], Rs[i][0][1], Rs[i][0][2]);
    printf("%g %g %g\n", Rs[i][1][0], Rs[i][1][1], Rs[i][1][2]);
    printf("%g %g %g\n", Rs[i][2][0], Rs[i][2][1], Rs[i][2][2]);

    printf("\n%g %g %g\n", ts[i][0], ts[i][1], ts[i][2]);
  }

  nsol=p3p_solve4_2Derr(&cal, m, M, R, t);
  printf("\nP4P: %d\n", nsol);
  if(nsol){
    printf("%g %g %g\n", R[0][0], R[0][1], R[0][2]);
    printf("%g %g %g\n", R[1][0], R[1][1], R[1][2]);
    printf("%g %g %g\n", R[2][0], R[2][1], R[2][2]);

    printf("\n%g %g %g\n", t[0], t[1], t[2]);
  }
 
  exit(0);
}
#endif
