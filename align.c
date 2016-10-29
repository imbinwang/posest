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
#include <math.h>
#include <string.h>

#include "compiler.h"
#include "sam.h"

static int jacobi_4x4(double *A, double *D, double *U)
{
  int iter;
  register int i, j, k;
  double B[4], Z[4];
  static double Id[16] = {1., 0., 0., 0., 
                          0., 1., 0., 0., 
                          0., 0., 1., 0., 
                          0., 0., 0., 1.};

  memcpy(U, Id, 16 * sizeof(double));

  B[0] = A[0]; B[1] = A[5]; B[2] = A[10]; B[3] = A[15];
  memcpy(D, B, 4 * sizeof(double));
  memset(Z, 0, 4 * sizeof(double));

  for(iter = 0; iter < 50; iter++)
  {
    double sum, thresh;

    sum = fabs(A[1]) + fabs(A[2]) + fabs(A[3]) + fabs(A[6]) + fabs(A[7]) + fabs(A[11]);

    if (sum == 0.0)
      return 1;

    thresh =  (iter < 3) ? 0.2 * sum / 16. : 0.0;
    for(i = 0; i < 3; i++)
    {
      double *pAij = A + 5 * i + 1;
      for(j = i + 1 ; j < 4; j++)
      {
        double Aij = *pAij;
        double eps_machine = 100.0 * fabs(Aij);

        if ( iter > 3 && fabs(D[i]) + eps_machine == fabs(D[i]) && fabs(D[j]) + eps_machine == fabs(D[j]) )
          *pAij = 0.0;
        else if (fabs(Aij) > thresh) 
        {
          double h = D[j] - D[i], t;
          double c, s, tau;

          if (fabs(h) + eps_machine == fabs(h))
            t = Aij / h;
          else 
          {
            double theta = 0.5 * h / Aij;
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0) t = -t;
          }

          h = t * Aij;
          Z[i] -= h;
          Z[j] += h;
          D[i] -= h;
          D[j] += h;
          *pAij = 0.0;

          c = 1.0 / sqrt(1 + t * t);
          s = t * c;
          tau = s / (1.0 + c);
          for(k = 0; k <= i - 1; k++) 
          { 
            double g = A[k * 4 + i], h = A[k * 4 + j];
            A[k * 4 + i] = g - s * (h + g * tau);
            A[k * 4 + j] = h + s * (g - h * tau);
          }
          for(k = i + 1; k <= j - 1; k++)
          {
            double g = A[i * 4 + k], h = A[k * 4 + j];
            A[i * 4 + k] = g - s * (h + g * tau);
            A[k * 4 + j] = h + s * (g - h * tau);
          }
          for(k = j + 1; k < 4; k++)
          {
            double g = A[i * 4 + k], h = A[j * 4 + k];
            A[i * 4 + k] = g - s * (h + g * tau); 
            A[j * 4 + k] = h + s * (g - h * tau); 
          }
          for(k = 0; k < 4; k++)
          {
            double g = U[k * 4 + i], h = U[k * 4 + j];
            U[k * 4 + i] = g - s * (h + g * tau); 
            U[k * 4 + j] = h + s * (g - h * tau); 
          }
        }
        pAij++;
      }
    }

    for(i = 0; i < 4; i++) B[i] += Z[i];
    memcpy(D, B, 4 * sizeof(double));
    memset(Z, 0, 4 * sizeof(double));
  } 

  return 0;
}

/* align two 3-point sets.
 * Returns 0 on success
 */
int posest_align3Pts(double M_end[3][3], 
                     double XYZ[3][3],
                     double R[3][3], double T[3])
{
  register int i, j;
  double X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2;
  double C_start[3], C_end[3];
  double s[3 * 3];
  double Qs[16], evs[4], U[16];
  int i_ev;
  double ev_max, q[4];
  double q02, q12, q22, q32, q0_1, q0_2, q0_3, q1_2, q1_3, q2_3;

  X0=XYZ[0][0]; Y0=XYZ[0][1]; Z0=XYZ[0][2];
  X1=XYZ[1][0]; Y1=XYZ[1][1]; Z1=XYZ[1][2];
  X2=XYZ[2][0]; Y2=XYZ[2][1]; Z2=XYZ[2][2];

  /* centroids: */
  C_end[0] = (M_end[0][0] + M_end[1][0] + M_end[2][0]) / 3.0;
  C_end[1] = (M_end[0][1] + M_end[1][1] + M_end[2][1]) / 3.0;
  C_end[2] = (M_end[0][2] + M_end[1][2] + M_end[2][2]) / 3.0;

  C_start[0] = (X0 + X1 + X2) / 3.0;
  C_start[1] = (Y0 + Y1 + Y2) / 3.0;
  C_start[2] = (Z0 + Z1 + Z2) / 3.0;

  /* covariance matrix s: */
  for(j = 0; j < 3; j++)
  {
    s[0 * 3 + j] = (X0 * M_end[0][j] + X1 * M_end[1][j] + X2 * M_end[2][j]) / 3.0 - C_end[j] * C_start[0];
    s[1 * 3 + j] = (Y0 * M_end[0][j] + Y1 * M_end[1][j] + Y2 * M_end[2][j]) / 3.0 - C_end[j] * C_start[1];
    s[2 * 3 + j] = (Z0 * M_end[0][j] + Z1 * M_end[1][j] + Z2 * M_end[2][j]) / 3.0 - C_end[j] * C_start[2];
  }

  Qs[0 * 4 + 0] = s[0 * 3 + 0] + s[1 * 3 + 1] + s[2 * 3 + 2];
  Qs[1 * 4 + 1] = s[0 * 3 + 0] - s[1 * 3 + 1] - s[2 * 3 + 2];
  Qs[2 * 4 + 2] = s[1 * 3 + 1] - s[2 * 3 + 2] - s[0 * 3 + 0];
  Qs[3 * 4 + 3] = s[2 * 3 + 2] - s[0 * 3 + 0] - s[1 * 3 + 1];

  Qs[1 * 4 + 0] = Qs[0 * 4 + 1] = s[1 * 3 + 2] - s[2 * 3 + 1]; 
  Qs[2 * 4 + 0] = Qs[0 * 4 + 2] = s[2 * 3 + 0] - s[0 * 3 + 2];
  Qs[3 * 4 + 0] = Qs[0 * 4 + 3] = s[0 * 3 + 1] - s[1 * 3 + 0];
  Qs[2 * 4 + 1] = Qs[1 * 4 + 2] = s[1 * 3 + 0] + s[0 * 3 + 1];
  Qs[3 * 4 + 1] = Qs[1 * 4 + 3] = s[2 * 3 + 0] + s[0 * 3 + 2];
  Qs[3 * 4 + 2] = Qs[2 * 4 + 3] = s[2 * 3 + 1] + s[1 * 3 + 2];

  jacobi_4x4(Qs, evs, U);

  /* looking for the largest eigen value: */
  i_ev = 0;
  ev_max = evs[i_ev];
  for(i = 1; i < 4; i++)
    if (evs[i] > ev_max)
      ev_max = evs[i_ev = i];

  /* quaternion: */
/***
  for(i = 0; i < 4; i++)
    q[i] = U[i * 4 + i_ev];
***/
  q[0] = U[0 * 4 + i_ev];
  q[1] = U[1 * 4 + i_ev];
  q[2] = U[2 * 4 + i_ev];
  q[3] = U[3 * 4 + i_ev];

  q02 = q[0] * q[0]; q12 = q[1] * q[1]; q22 = q[2] * q[2]; q32 = q[3] * q[3];
  q0_1 = q[0] * q[1]; q0_2 = q[0] * q[2]; q0_3 = q[0] * q[3];
  q1_2 = q[1] * q[2]; q1_3 = q[1] * q[3];
  q2_3 = q[2] * q[3];

  R[0][0] = q02 + q12 - q22 - q32;
  R[0][1] = 2. * (q1_2 - q0_3);
  R[0][2] = 2. * (q1_3 + q0_2);

  R[1][0] = 2. * (q1_2 + q0_3);
  R[1][1] = q02 + q22 - q12 - q32;
  R[1][2] = 2. * (q2_3 - q0_1);

  R[2][0] = 2. * (q1_3 - q0_2);
  R[2][1] = 2. * (q2_3 + q0_1);
  R[2][2] = q02 + q32 - q12 - q22;

/****
  for(i = 0; i < 3; i++)
    T[i] = C_end[i] - (R[i][0] * C_start[0] + R[i][1] * C_start[1] + R[i][2] * C_start[2]);
****/

  T[0] = C_end[0] - (R[0][0] * C_start[0] + R[0][1] * C_start[1] + R[0][2] * C_start[2]);
  T[1] = C_end[1] - (R[1][0] * C_start[0] + R[1][1] * C_start[1] + R[1][2] * C_start[2]);
  T[2] = C_end[2] - (R[2][0] * C_start[0] + R[2][1] * C_start[1] + R[2][2] * C_start[2]);

  return 0;
}


/* LAPACK eigenvalues/eigenvectors */
extern int F77_FUNC(dsyev)(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);

/* Absolute orientation with the unit quaternion method
 *
 * see Horn, "Closed-form Solution of Absolute Orientation Using Unit Quaternions",
 * JOSAA (4):4, 1987, pp.629
 *
 * Given two corresponding point sets pts0, pts1, computes R, t such that
 * pts1 = R*pts0 + t
 *
 * Returns 0 on success, nonzero otherwise
 */
int posest_alignNPts(double (*pts0)[3], double (*pts1)[3], int npts, double R[9], double t[3])
{
register int i;
double pc[3], qc[3], p[3], q[3], M[4][4], *dptr, one_over_npts;
double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
double eigvals[4], *eigvecs=(double *)M;
register double *x;
double x0x0, x1x1, x2x2, x3x3, x0x1, x0x2, x0x3, x1x2, x1x3, x2x3;
int nvars=4, info;

  if(npts<3){
    fprintf(stderr, "At least 3 points are necessary for estimating absolute orientation in posest_alignNPts()! [%d]\n", npts); fflush(stderr);
    //memset(R, 0, 9*sizeof(double));
    //memset(t, 0, 3*sizeof(double));

    return 1;
  }

  pc[0]=pc[1]=pc[2]=
  qc[0]=qc[1]=qc[2]=0.0;
  one_over_npts=1.0/(double)(npts);
  Sxx=Syy=Szz=Syz=Szy=Szx=Sxz=Sxy=Syx=0.0;
	
  /* compute centroids */
  for(i=0; i<npts; i++){
    dptr=pts0[i];
    pc[0]+=*dptr++; pc[1]+=*dptr++; pc[2]+=*dptr;

    dptr=pts1[i];
    qc[0]+=*dptr++; qc[1]+=*dptr++; qc[2]+=*dptr;
  }

  pc[0]*=one_over_npts; pc[1]*=one_over_npts; pc[2]*=one_over_npts;
  qc[0]*=one_over_npts; qc[1]*=one_over_npts; qc[2]*=one_over_npts;

  /* compute the S,, from p', q' */
  for(i=0; i<npts; i++){
    dptr=pts0[i];
    p[0]=*dptr++ - pc[0];
    p[1]=*dptr++ - pc[1];
    p[2]=*dptr   - pc[2];

    dptr=pts1[i];
    q[0]=*dptr++ - qc[0];
    q[1]=*dptr++ - qc[1];
    q[2]=*dptr   - qc[2];

    Sxx+=p[0]*q[0]; Sxy+=p[0]*q[1]; Sxz+=p[0]*q[2];
    Syx+=p[1]*q[0]; Syy+=p[1]*q[1]; Syz+=p[1]*q[2];
    Szx+=p[2]*q[0]; Szy+=p[2]*q[1]; Szz+=p[2]*q[2];
  }

  /* M is symmetric, hence only its upper triangle is computed. Note that in column-major, this corresponds to the lower triangle */
  M[0][0]=Sxx + Syy + Szz; M[0][1]=Syz - Szy;       M[0][2]=Szx - Sxz;        M[0][3]=Sxy - Syx;
                           M[1][1]=Sxx - Syy - Szz; M[1][2]=Sxy + Syx;        M[1][3]=Szx + Sxz;
                                                    M[2][2]=-Sxx + Syy - Szz; M[2][3]=Syz + Szy;
                                                                              M[3][3]=-Sxx - Syy + Szz;

  {
    double work[136];
    int lwork=136; // optimal size determined by querying with lwork==-1

    F77_FUNC(dsyev)("V", "L", &nvars, (double *)M, &nvars, eigvals, work, &lwork, &info);
  }

  if(info<0){
    fprintf(stderr, "LAPACK error: illegal value for argument %d of dspev/dsyev in posest_alignNPts()\n", -info);
    exit(1);
  }
  else if(info>0){
    fprintf(stderr, "LAPACK error: dspev/dsyev failed to converge in posest_alignNPts();\n%d %s", info,
        "off-diagonal elements of an intermediate tridiagonal form did not converge to zero\n");
    //memset(R, 0, 9*sizeof(double));
    //memset(t, 0, 3*sizeof(double));

    return 1;
  }

  /* eigenvalues are returned in ascending order, therefore the largest is the last one */
  x=eigvecs + 4*3;

  /* x points to the computed quaternion; by computation it has a unit length */

  /* compute the rotation matrix corresponding to the estimated quaternion */
  x0x0=x[0]*x[0];
  x1x1=x[1]*x[1];
  x2x2=x[2]*x[2];
  x3x3=x[3]*x[3];

  x0x1=x[0]*x[1];
  x0x2=x[0]*x[2];
  x0x3=x[0]*x[3];
  x1x2=x[1]*x[2];
  x1x3=x[1]*x[3];
  x2x3=x[2]*x[3];

  R[0]=x0x0+x1x1-x2x2-x3x3;
  R[1]=2*(x1x2-x0x3);
  R[2]=2*(x1x3+x0x2);

  R[3]=2*(x1x2+x0x3);
  R[4]=x0x0+x2x2-x1x1-x3x3;
  R[5]=2*(x2x3-x0x1);

  R[6]=2*(x1x3-x0x2);
  R[7]=2*(x2x3+x0x1);
  R[8]=x0x0+x3x3-x1x1-x2x2;

  if(sam_det3x3(R)<0){ // left-handed R, negate 3rd row
    R[6]=-R[6];
    R[7]=-R[7];
    R[8]=-R[8];
  }

  /* t = qc - R *pc */
  t[0]=qc[0] - (R[0]*pc[0] + R[1]*pc[1] + R[2]*pc[2]);
  t[1]=qc[1] - (R[3]*pc[0] + R[4]*pc[1] + R[5]*pc[2]);
  t[2]=qc[2] - (R[6]*pc[0] + R[7]*pc[1] + R[8]*pc[2]);

  return 0;
}
