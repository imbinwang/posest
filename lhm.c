/*////////////////////////////////////////////////////////////////////////////////
// 
//  Pose estimation using the algorithm of Lu, Hager and Mjolsness (LHM)
//
//  Copyright (C) 2002-12  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//////////////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "compiler.h"
#include "sam.h"

#include "svd3.h"

#include "posest.h" 

#define USE_LQS_FIT    1 // use LQS (LMedS) if 1
#define USE_RANSAC_FIT 0 // use RANSAC if 1
#define USE_PROSAC_FIT 0 // use PROSAC if 1

#if USE_LQS_FIT + USE_RANSAC_FIT + USE_PROSAC_FIT != 1
#error Exactly one of the USE_XXX_FIT macros should be defined as 1!
#endif

#define LHM_RANSAC_OUTL_THRESH  1E-03 // relevant to RANSAC/PROSAC only. This is in normalized coordinates!

#define LHM_TOLERANCE     1E-5 /* convergence tolerance */
#define LHM_EPSILON       1E-8 /* lower bound of objective function */
#define LHM_MAXITER       35   /* maximum number of iterations */

#define _SQR(x) ((x)*(x))


static int lhm_estPose(double (*proj2D)[2], double (*pts3D)[3], double (*wkpts3D)[3], int *ptsidx, int npts, double R[9], double t[3], int verbose);
static int lhm_calcRt(double (*pts0)[3], double (*pts1)[3], int *indx, int npts, double (*V)[9], double Tfact[9], double R[9], double t[3]);
static void lhm_tFromR(double (*pts3D)[3], int *ptsidx, int npts, double (*V)[9], double R[9], double Tfact[9], double t[3]);
static double lhm_objSpaceErr(double (*pts3D)[3], int *ptsidx, int npts, double (*V)[9]);
static double lhm_imgSpaceErr(double (*pts2D)[2], double (*pts3D)[3], int *ptsidx, int npts, double R[9], double t[3]);

/* estimate 3D pose based on the 3D - 2D correspondences specified by pts3D, proj2D.
 * 2D points are assumed to be expressed in normalized coordinates, i.e. the effect of
 * the camera intrinsics has been canceled out.
 * wkpts3D is an auxiliary work array allocated by the caller
 * R should contain an initial estimate for the rotation
 *
 * ptsidx is an index that defines the points that are to be used
 * in the estimation. If ptsidx==NULL all points are employed
 *
 * Upon return, R and t contain the rotation & translation of the estimated
 * rigid motion.
 *
 * The algorithm used for estimating R, t is that by Lu et al.
 * Lu, Hager & Mjolsness "Fast and Globally Convergent Pose Estimation
 * from Video Images", PAMI (22):6, 2000, pp.613
 *
 * returns 0 on success, nonzero otherwise
 *
 */
static int lhm_estPose(double (*proj2D)[2], double (*pts3D)[3], double (*wkpts3D)[3], int *ptsidx, int npts, double R[9], double t[3], int verbose)
{
double (*V)[9]; /* holds the projection matrices of Eq.(18). Note that V does not use ptsidx for indexing! */
double Tfact[9], tmpm[9], tmpv[3];

register int i, j, k;
int niter;
register double *Vp, *p2, *p3, *wp3;
double mag;
double cur_err=DBL_MAX, old_err;

  if((V=(double (*)[9])malloc(npts*sizeof(double[9])))==NULL){
    fprintf(stderr, "memory allocation request failed in lhm_estPose()\n");
    exit(1);
  }

  /* clear tmpm */
  memset(tmpm, 0, 9*sizeof(double));

  /* compute the Vi and tmpm as \sum_k Vk */
  if(ptsidx==NULL){ /* no index, use all points */
    for(i=0; i<npts; ++i){
      p2=proj2D[i];
      mag=_SQR(p2[0]) + _SQR(p2[1]) + 1.0;
      mag=1.0/mag;

      Vp=V[i];
      tmpm[0]+=Vp[0]=p2[0]*p2[0]*mag;
      tmpm[1]+=Vp[1]=p2[0]*p2[1]*mag;
      tmpm[2]+=Vp[2]=p2[0]* 1.0 *mag;

      tmpm[3]+=Vp[3]=p2[1]*p2[0]*mag;
      tmpm[4]+=Vp[4]=p2[1]*p2[1]*mag;
      tmpm[5]+=Vp[5]=p2[1]* 1.0 *mag;

      tmpm[6]+=Vp[6]=p2[0]*mag;
      tmpm[7]+=Vp[7]=p2[1]*mag;
      tmpm[8]+=Vp[8]= 1.0 *mag;
    }
  }
  else{
    for(j=0; j<npts; ++j){
      i=ptsidx[j];

      p2=proj2D[i];
      mag=_SQR(p2[0]) + _SQR(p2[1]) + 1.0;
      mag=1.0/mag;

      Vp=V[j];
      tmpm[0]+=Vp[0]=p2[0]*p2[0]*mag;
      tmpm[1]+=Vp[1]=p2[0]*p2[1]*mag;
      tmpm[2]+=Vp[2]=p2[0]* 1.0 *mag;

      tmpm[3]+=Vp[3]=p2[1]*p2[0]*mag;
      tmpm[4]+=Vp[4]=p2[1]*p2[1]*mag;
      tmpm[5]+=Vp[5]=p2[1]* 1.0 *mag;

      tmpm[6]+=Vp[6]=p2[0]*mag;
      tmpm[7]+=Vp[7]=p2[1]*mag;
      tmpm[8]+=Vp[8]= 1.0 *mag;
    }
  }

  /* compute the first coefficient of Eq.(20):
   * 1/n (I - 1/n \sum_k Vk)^-1
   */

  /* tmpm is multiplied by n so that its inverse (computed below) is
   * implicitly multiplied by 1/n:  tmpm = n*(I - 1/n \sum_k Vk)
   */
  tmpm[0]=npts - tmpm[0]; tmpm[1]=     - tmpm[1]; tmpm[2]=     - tmpm[2];
  tmpm[3]=     - tmpm[3]; tmpm[4]=npts - tmpm[4]; tmpm[5]=     - tmpm[5];
  tmpm[6]=     - tmpm[6]; tmpm[7]=     - tmpm[7]; tmpm[8]=npts - tmpm[8]; 

  /* Tfact = 1/n*(I - 1/n \sum_k Vk)^-1 */
  sam_inv3x3(tmpm, Tfact);

  /* compute the initial estimate of translation */
  lhm_tFromR(pts3D, ptsidx, npts, V, R, Tfact, t);

  niter=0;
  do{
    /* transform 3D points using current R, t:  wkpts3D = R*pts3D + t */
    if(ptsidx==NULL){
      for(i=0; i<npts; ++i){
        p3=pts3D[i]; wp3=wkpts3D[i];
        wp3[0]=R[0]*p3[0] + R[1]*p3[1] + R[2]*p3[2] + t[0];
        wp3[1]=R[3]*p3[0] + R[4]*p3[1] + R[5]*p3[2] + t[1];
        wp3[2]=R[6]*p3[0] + R[7]*p3[1] + R[8]*p3[2] + t[2];
      }
    }
    else{
      for(j=0; j<npts; ++j){
        i=ptsidx[j];

        p3=pts3D[i]; wp3=wkpts3D[i];
        wp3[0]=R[0]*p3[0] + R[1]*p3[1] + R[2]*p3[2] + t[0];
        wp3[1]=R[3]*p3[0] + R[4]*p3[1] + R[5]*p3[2] + t[1];
        wp3[2]=R[6]*p3[0] + R[7]*p3[1] + R[8]*p3[2] + t[2];
      }
    }

    old_err=cur_err;
    cur_err=lhm_objSpaceErr(wkpts3D, ptsidx, npts, V);

    /* wkpts3D[k] = V[k] * (R*pts3D[k]+t) */
    if(ptsidx==NULL){
      for(k=0; k<npts; ++k){
        /* tmpv = V[k] * wkpts3D[k] */
        Vp=V[k]; wp3=wkpts3D[k];
        tmpv[0]=Vp[0]*wp3[0] + Vp[1]*wp3[1] + Vp[2]*wp3[2];
        tmpv[1]=Vp[3]*wp3[0] + Vp[4]*wp3[1] + Vp[5]*wp3[2];
        tmpv[2]=Vp[6]*wp3[0] + Vp[7]*wp3[1] + Vp[8]*wp3[2];

        wp3[0]=tmpv[0]; wp3[1]=tmpv[1]; wp3[2]=tmpv[2];
      }
    }
    else{
      for(j=0; j<npts; ++j){
        k=ptsidx[j];
        /* tmpv = V[j] * wkpts3D[k] */
        Vp=V[j]; wp3=wkpts3D[k];
        tmpv[0]=Vp[0]*wp3[0] + Vp[1]*wp3[1] + Vp[2]*wp3[2]; 
        tmpv[1]=Vp[3]*wp3[0] + Vp[4]*wp3[1] + Vp[5]*wp3[2]; 
        tmpv[2]=Vp[6]*wp3[0] + Vp[7]*wp3[1] + Vp[8]*wp3[2]; 

        wp3[0]=tmpv[0]; wp3[1]=tmpv[1]; wp3[2]=tmpv[2];
      }
    }

    /* update R,t */
    lhm_calcRt(pts3D, wkpts3D, ptsidx, npts, V, Tfact, R, t);
    niter++;
  } while ((fabs((old_err-cur_err)/old_err) > LHM_TOLERANCE) && (cur_err > LHM_EPSILON) && niter<LHM_MAXITER);

  if(verbose){
    old_err=sqrt(lhm_imgSpaceErr(proj2D, pts3D, ptsidx, npts, R, t)/(double)npts);
    fprintf(stderr, "\n%d iterations,  obj space error %.8g -- img space error %.8g\n", niter, sqrt(cur_err/(double)npts), old_err);
  }

  free(V);

#if 0
  if(npts==3){
    double K[9]={
        3252.33, 0, 650.247, 0, 3256.29, 414.877, 0, 0, 1
    };
    double x, y, Z, xx, yy;


    printf("R, t\n");
    printf("%g %g %g\n", R[0], R[1], R[2]);
    printf("%g %g %g\n", R[3], R[4], R[5]);
    printf("%g %g %g\n", R[6], R[7], R[8]);
    printf("   %g %g %g\n\n", t[0], t[1], t[2]);
    for(j=0; j<npts; ++j){
      i=ptsidx[j];

      printf("Point %d: %g %g %g\n", i, pts3D[i][0], pts3D[i][1], pts3D[i][2]);
      /* transform point */
      x=R[0]*pts3D[i][0] + R[1]*pts3D[i][1] + R[2]*pts3D[i][2] + t[0];
      y=R[3]*pts3D[i][0] + R[4]*pts3D[i][1] + R[5]*pts3D[i][2] + t[1];
      Z=R[6]*pts3D[i][0] + R[7]*pts3D[i][1] + R[8]*pts3D[i][2] + t[2];
      x/=Z; y/=Z;

      x=K[0]*x+K[1]*y+K[2];
      y=       K[4]*y+K[5];

      xx=K[0]*proj2D[i][0]+K[1]*proj2D[i][1]+K[2];
      yy=                  K[4]*proj2D[i][1]+K[5];

      printf("  %g %g -- %g %g\n", xx, yy, x, y);
    }
    printf("\n");
  }
#endif

  return POSEST_OK;
}

/* Estimate the rigid transformation between two point sets using SVD
 * If the resulting translation has a T.z<0, R,t are refined to ensure T.z>0
 *
 * See
 * Arun, Huang & Blostein "Least-Squares Fitting of Two 3-D Point Sets",
 * PAMI 9:(5), 1987, p. 698-700
 * and
 * Lu, Hager & Mjolsness "Fast and Globally Convergent Pose Estimation
 * from Video Images", PAMI (22):6, 2000, pp.613
 *
 * An indx!=NULL specifies the indices of corresponding points that should be
 * considered; all points are used otherwise (no mismatches are tolerated)
 *
 * returns 0 on success, nonzero otherwise
 */
static int lhm_calcRt(double (*pts0)[3], double (*pts1)[3], int *indx, int npts, double (*V)[9], double Tfact[9],
                      double R[9], double t[3])
{
register int i;
double pc[3], qc[3], p[3], q[3], M[3][3], *dptr, one_over_npts, dpsum, dqsum;
double S[3], U[3][3], Vt[3][3];
extern int F77_FUNC(dgesvd)(char *jobu, char *jobvt, int *m, int *n,
                            double *a, int *lda, double *s, double *u, int *ldu,
                            double *vt, int *ldvt, double *work, int *lwork, int *info);

  if(npts<3){
    fprintf(stderr, "At least 3 points are necessary for estimating R, t in lhm_calcRt()! [%d]\n", npts); fflush(stderr);
    memset(R, 0, 9*sizeof(double));
    memset(t, 0, 3*sizeof(double));

    return 1;
  }

  pc[0]=pc[1]=pc[2]=
  qc[0]=qc[1]=qc[2]=0.0;
  one_over_npts=1.0/(double)(npts);
  dpsum=dqsum=0.0;
  M[0][0]=M[0][1]=M[0][2]=
  M[1][0]=M[1][1]=M[1][2]=
  M[2][0]=M[2][1]=M[2][2]=0.0;
	
  /* Note: LHM erroneously define M in Eq.(13) as \sum_i q'_i*p't_i; the former is
   * the transpose of M, M is actually \sum_i p'_i*q't_i. See also Arun et al.
   */
  if(!indx){
    /* compute centroids */
    for(i=0; i<npts; i++){
      dptr=pts0[i];
      pc[0]+=*dptr++; pc[1]+=*dptr++; pc[2]+=*dptr;

      dptr=pts1[i];
      qc[0]+=*dptr++; qc[1]+=*dptr++; qc[2]+=*dptr;
    }

    pc[0]*=one_over_npts; pc[1]*=one_over_npts; pc[2]*=one_over_npts;
    qc[0]*=one_over_npts; qc[1]*=one_over_npts; qc[2]*=one_over_npts;

    /* compute M from p', q' */
    for(i=0; i<npts; i++){
      dptr=pts0[i];
      p[0]=*dptr++ - pc[0];
      p[1]=*dptr++ - pc[1];
      p[2]=*dptr   - pc[2];
      dpsum+=p[0]*p[0] + p[1]*p[1] + p[2]*p[2];

      dptr=pts1[i];
      q[0]=*dptr++ - qc[0];
      q[1]=*dptr++ - qc[1];
      q[2]=*dptr   - qc[2];
      dqsum+=q[0]*q[0] + q[1]*q[1] + q[2]*q[2];

      /* compute \sum_i q'_i*p't_i which is M in column major */
      M[0][0]+=q[0]*p[0]; M[0][1]+=q[0]*p[1]; M[0][2]+=q[0]*p[2];
      M[1][0]+=q[1]*p[0]; M[1][1]+=q[1]*p[1]; M[1][2]+=q[1]*p[2];
      M[2][0]+=q[2]*p[0]; M[2][1]+=q[2]*p[1]; M[2][2]+=q[2]*p[2];
    }
  }
  else{
    /* compute centroids */
    for(i=0; i<npts; i++){
      dptr=pts0[indx[i]];
      pc[0]+=*dptr++; pc[1]+=*dptr++; pc[2]+=*dptr;

      dptr=pts1[indx[i]];
      qc[0]+=*dptr++; qc[1]+=*dptr++; qc[2]+=*dptr;
    }

    pc[0]*=one_over_npts; pc[1]*=one_over_npts; pc[2]*=one_over_npts;
    qc[0]*=one_over_npts; qc[1]*=one_over_npts; qc[2]*=one_over_npts;

    /* compute M from p', q' */
    for(i=0; i<npts; i++){
      dptr=pts0[indx[i]];
      p[0]=*dptr++ - pc[0];
      p[1]=*dptr++ - pc[1];
      p[2]=*dptr   - pc[2];
      dpsum+=p[0]*p[0] + p[1]*p[1] + p[2]*p[2];

      dptr=pts1[indx[i]];
      q[0]=*dptr++ - qc[0];
      q[1]=*dptr++ - qc[1];
      q[2]=*dptr   - qc[2];
      dqsum+=q[0]*q[0] + q[1]*q[1] + q[2]*q[2];

      /* compute \sum_i q'_i*p't_i which is M in column major */
      M[0][0]+=q[0]*p[0]; M[0][1]+=q[0]*p[1]; M[0][2]+=q[0]*p[2];
      M[1][0]+=q[1]*p[0]; M[1][1]+=q[1]*p[1]; M[1][2]+=q[1]*p[2];
      M[2][0]+=q[2]*p[0]; M[2][1]+=q[2]*p[1]; M[2][2]+=q[2]*p[2];
    }
  }

  //printf("lhm_calcRt() scale: %g\n", sqrt(dqsum/dpsum)); // see also Eq.(39) in Lu et al.

  /* note that M has been computed above in column-major order */

#if 1 // use custom 3x3 SVD
    svd3(U[0], S, Vt[0], M[0]);
    /* svd3 actually returns V; for compatibility with LAPACK which returns V^T,
     * the computed Vt is transposed in place to yield the true Vt
     */
    { double tmp;
    tmp=Vt[0][1]; Vt[0][1]=Vt[1][0]; Vt[1][0]=tmp;
    tmp=Vt[0][2]; Vt[0][2]=Vt[2][0]; Vt[2][0]=tmp;
    tmp=Vt[1][2]; Vt[1][2]=Vt[2][1]; Vt[2][1]=tmp;
    }
#else // use generic SVD
  {
  double work[32]; /* SVD stuff */
  int info, three=3, lwork=32; /* probably too much here...*/

  F77_FUNC(dgesvd)("A", "A", &three, &three, M[0], &three, S, U[0], &three, Vt[0], &three, work, &lwork, &info);
  if(info<0){
    fprintf(stderr, "LAPACK error: illegal value for argument %d of dgesvd in lhm_calcRt()\n", -info);
    exit(1);
  }
  else if(info>0){
    fprintf(stderr, "LAPACK error: dbdsqr failed to converge in lhm_calcRt();\n%d %s", info,
        "superdiagonals of an intermediate bidiagonal form did not converge to zero\n");
    memset(R, 0, 9*sizeof(double));
    memset(t, 0, 3*sizeof(double));

    return 1;
  }
  }
#endif

  /* compute R as V * U^t. Remember that Vt & U are in column major */
  R[0]=Vt[0][0]*U[0][0] + Vt[0][1]*U[1][0] + Vt[0][2]*U[2][0];
  R[1]=Vt[0][0]*U[0][1] + Vt[0][1]*U[1][1] + Vt[0][2]*U[2][1];
  R[2]=Vt[0][0]*U[0][2] + Vt[0][1]*U[1][2] + Vt[0][2]*U[2][2];

  R[3]=Vt[1][0]*U[0][0] + Vt[1][1]*U[1][0] + Vt[1][2]*U[2][0];
  R[4]=Vt[1][0]*U[0][1] + Vt[1][1]*U[1][1] + Vt[1][2]*U[2][1];
  R[5]=Vt[1][0]*U[0][2] + Vt[1][1]*U[1][2] + Vt[1][2]*U[2][2];

  R[6]=Vt[2][0]*U[0][0] + Vt[2][1]*U[1][0] + Vt[2][2]*U[2][0];
  R[7]=Vt[2][0]*U[0][1] + Vt[2][1]*U[1][1] + Vt[2][2]*U[2][1];
  R[8]=Vt[2][0]*U[0][2] + Vt[2][1]*U[1][2] + Vt[2][2]*U[2][2];

  /* currently no collinearity check for the p, q. They are collinear
   * when two of the singular values are equal
   */

  if(sam_det3x3(R)<0){ // R is a reflection, see p. 699 in Arun et al.
    if(S[2]>S[0]*1E-3){ // case b): smallest singular value (S[2]) is nonzero
      memset(R, 0, 9*sizeof(double));
      memset(t, 0, 3*sizeof(double));

      return 2;
    }

    /* case a): change the sign of V's 3rd colun and recompute R */
    // R = [V(:, 1:2) -V(:, 3)] * U^t
    Vt[0][2]=-Vt[0][2]; Vt[1][2]=-Vt[1][2]; Vt[2][2]=-Vt[2][2];
    R[0]=Vt[0][0]*U[0][0] + Vt[0][1]*U[1][0] + Vt[0][2]*U[2][0];
    R[1]=Vt[0][0]*U[0][1] + Vt[0][1]*U[1][1] + Vt[0][2]*U[2][1];
    R[2]=Vt[0][0]*U[0][2] + Vt[0][1]*U[1][2] + Vt[0][2]*U[2][2];

    R[3]=Vt[1][0]*U[0][0] + Vt[1][1]*U[1][0] + Vt[1][2]*U[2][0];
    R[4]=Vt[1][0]*U[0][1] + Vt[1][1]*U[1][1] + Vt[1][2]*U[2][1];
    R[5]=Vt[1][0]*U[0][2] + Vt[1][1]*U[1][2] + Vt[1][2]*U[2][2];

    R[6]=Vt[2][0]*U[0][0] + Vt[2][1]*U[1][0] + Vt[2][2]*U[2][0];
    R[7]=Vt[2][0]*U[0][1] + Vt[2][1]*U[1][1] + Vt[2][2]*U[2][1];
    R[8]=Vt[2][0]*U[0][2] + Vt[2][1]*U[1][2] + Vt[2][2]*U[2][2];
  }

compute_trans:
  lhm_tFromR(pts0, indx, npts, V, R, Tfact, t);

  if(t[2]<0.0){ // t should be inverted so that object is in front of the camera
    // R = -[V(:, 1:2) -V(:, 3)] * U^t
    Vt[0][2]=-Vt[0][2]; Vt[1][2]=-Vt[1][2]; Vt[2][2]=-Vt[2][2];
    R[0]=-(Vt[0][0]*U[0][0] + Vt[0][1]*U[1][0] + Vt[0][2]*U[2][0]);
    R[1]=-(Vt[0][0]*U[0][1] + Vt[0][1]*U[1][1] + Vt[0][2]*U[2][1]);
    R[2]=-(Vt[0][0]*U[0][2] + Vt[0][1]*U[1][2] + Vt[0][2]*U[2][2]);

    R[3]=-(Vt[1][0]*U[0][0] + Vt[1][1]*U[1][0] + Vt[1][2]*U[2][0]);
    R[4]=-(Vt[1][0]*U[0][1] + Vt[1][1]*U[1][1] + Vt[1][2]*U[2][1]);
    R[5]=-(Vt[1][0]*U[0][2] + Vt[1][1]*U[1][2] + Vt[1][2]*U[2][2]);

    R[6]=-(Vt[2][0]*U[0][0] + Vt[2][1]*U[1][0] + Vt[2][2]*U[2][0]);
    R[7]=-(Vt[2][0]*U[0][1] + Vt[2][1]*U[1][1] + Vt[2][2]*U[2][1]);
    R[8]=-(Vt[2][0]*U[0][2] + Vt[2][1]*U[1][2] + Vt[2][2]*U[2][2]);

    /* recompute t */
    goto compute_trans; // new t[2]>0 so goto will not return to this branch
  }

  return 0;
}

/* compute t given R -- see Eq.(20)
 * The first part of this equation is given by Tfact
 */
static void lhm_tFromR(double (*pts3D)[3], int *ptsidx, int npts, double (*V)[9], double R[9], double Tfact[9], double t[3])
{
double sumv[3], tmpv1[3], tmpv2[3], tmpm[9];
register int k, i;
register double *p3, *Vp;

  sumv[0]=sumv[1]=sumv[2]=0.0;

  /* sumv = \sum_k (Vk - I)*R*pk */
  if(ptsidx==NULL){
    for(k=0; k<npts; ++k){
      /* tmpv1 = R*pk */
      p3=pts3D[k];
      tmpv1[0]=R[0]*p3[0] + R[1]*p3[1] + R[2]*p3[2];
      tmpv1[1]=R[3]*p3[0] + R[4]*p3[1] + R[5]*p3[2];
      tmpv1[2]=R[6]*p3[0] + R[7]*p3[1] + R[8]*p3[2];

      /* tmpm = (Vk - I) */
      Vp=V[k];
      tmpm[0]=Vp[0] - 1.0; tmpm[1]=Vp[1]      ; tmpm[2]=Vp[2]      ;
      tmpm[3]=Vp[3]      ; tmpm[4]=Vp[4] - 1.0; tmpm[5]=Vp[5]      ;
      tmpm[6]=Vp[6]      ; tmpm[7]=Vp[7]      ; tmpm[8]=Vp[8] - 1.0; 

      /* tmpv2 = tmpm*tmpv1 = (Vk - I)*R*pk */
      tmpv2[0]=tmpm[0]*tmpv1[0] + tmpm[1]*tmpv1[1] + tmpm[2]*tmpv1[2];
      tmpv2[1]=tmpm[3]*tmpv1[0] + tmpm[4]*tmpv1[1] + tmpm[5]*tmpv1[2];
      tmpv2[2]=tmpm[6]*tmpv1[0] + tmpm[7]*tmpv1[1] + tmpm[8]*tmpv1[2];

      sumv[0]+=tmpv2[0]; sumv[1]+=tmpv2[1]; sumv[2]+=tmpv2[2];
   }
  }
   else{
    for(i=0; i<npts; i++){
      k=ptsidx[i];

      /* tmpv1 = R*pk */
      p3=pts3D[k];
      tmpv1[0]=R[0]*p3[0] + R[1]*p3[1] + R[2]*p3[2];
      tmpv1[1]=R[3]*p3[0] + R[4]*p3[1] + R[5]*p3[2];
      tmpv1[2]=R[6]*p3[0] + R[7]*p3[1] + R[8]*p3[2];

      /* tmpm = (Vi - I) */
      Vp=V[i];
      tmpm[0]=Vp[0] - 1.0; tmpm[1]=Vp[1]      ; tmpm[2]=Vp[2]      ;
      tmpm[3]=Vp[3]      ; tmpm[4]=Vp[4] - 1.0; tmpm[5]=Vp[5]      ;
      tmpm[6]=Vp[6]      ; tmpm[7]=Vp[7]      ; tmpm[8]=Vp[8] - 1.0; 

      /* tmpv2 = tmpm*tmpv1 = (Vi - I)*R*pk */
      tmpv2[0]=tmpm[0]*tmpv1[0] + tmpm[1]*tmpv1[1] + tmpm[2]*tmpv1[2];
      tmpv2[1]=tmpm[3]*tmpv1[0] + tmpm[4]*tmpv1[1] + tmpm[5]*tmpv1[2];
      tmpv2[2]=tmpm[6]*tmpv1[0] + tmpm[7]*tmpv1[1] + tmpm[8]*tmpv1[2];

      sumv[0]+=tmpv2[0]; sumv[1]+=tmpv2[1]; sumv[2]+=tmpv2[2];
    }
   }

   /* t=Tfact*sumv */
   t[0]=Tfact[0]*sumv[0] + Tfact[1]*sumv[1] + Tfact[2]*sumv[2];
   t[1]=Tfact[3]*sumv[0] + Tfact[4]*sumv[1] + Tfact[5]*sumv[2];
   t[2]=Tfact[6]*sumv[0] + Tfact[7]*sumv[1] + Tfact[8]*sumv[2];
}

/* calculate the object-space collinearity error of Eqs.(17) - (19)
 * The 3D points in pts3D are assumed to have been rigidly transformed
 * with R, t
 */
static double lhm_objSpaceErr(double (*pts3D)[3], int *ptsidx, int npts, double (*V)[9])
{
double tmpv[3], tmpm[9];
double sum;
register int i, j;
register double *p3, *Vp;

  if(ptsidx==NULL){
    for(i=0, sum=0.0; i<npts; ++i){
      /* tmpm = (I - Vi) */
      Vp=V[i];
      tmpm[0]=1.0 - Vp[0]; tmpm[1]=    - Vp[1]; tmpm[2]=    - Vp[2];
      tmpm[3]=    - Vp[3]; tmpm[4]=1.0 - Vp[4]; tmpm[5]=    - Vp[5];
      tmpm[6]=    - Vp[6]; tmpm[7]=    - Vp[7]; tmpm[8]=1.0 - Vp[8];

      /* tmpv = (I - Vi)*Pi */
      p3=pts3D[i];
      tmpv[0]=tmpm[0]*p3[0] + tmpm[1]*p3[1] + tmpm[2]*p3[2];
      tmpv[1]=tmpm[3]*p3[0] + tmpm[4]*p3[1] + tmpm[5]*p3[2];
      tmpv[2]=tmpm[6]*p3[0] + tmpm[7]*p3[1] + tmpm[8]*p3[2];

      sum+=_SQR(tmpv[0]) + _SQR(tmpv[1]) + _SQR(tmpv[2]);
    }
  }
  else{
    for(j=0, sum=0.0; j<npts; j++){
      i=ptsidx[j];

      /* tmpm = (I - Vj) */
      Vp=V[j];
      tmpm[0]=1.0 - Vp[0]; tmpm[1]=    - Vp[1]; tmpm[2]=    - Vp[2];
      tmpm[3]=    - Vp[3]; tmpm[4]=1.0 - Vp[4]; tmpm[5]=    - Vp[5];
      tmpm[6]=    - Vp[6]; tmpm[7]=    - Vp[7]; tmpm[8]=1.0 - Vp[8];

      /* tmpv = (I - Vj)*Pi */
      p3=pts3D[i];
      tmpv[0]=tmpm[0]*p3[0] + tmpm[1]*p3[1] + tmpm[2]*p3[2];
      tmpv[1]=tmpm[3]*p3[0] + tmpm[4]*p3[1] + tmpm[5]*p3[2];
      tmpv[2]=tmpm[6]*p3[0] + tmpm[7]*p3[1] + tmpm[8]*p3[2];

      sum+=_SQR(tmpv[0]) + _SQR(tmpv[1]) + _SQR(tmpv[2]);
    }
  }

  return sum;
}

/* calculate the image-space error of Eq.(8)
 * The 3D points in pts3D should be first
 * rigidly transformed with R, t. Recall that
 * pts2D are normalized
 */
static double lhm_imgSpaceErr(double (*pts2D)[2], double (*pts3D)[3], int *ptsidx, int npts, double R[9], double t[3])
{
register int i, j;
register double sum, x, y, Z;
register double *p3, *p2;

  sum=0.0;
  if(ptsidx==NULL){
    for(i=0; i<npts; ++i){
      /* transform 3D point */
      p3=pts3D[i];
      x=R[0]*p3[0] + R[1]*p3[1] + R[2]*p3[2] + t[0];
      y=R[3]*p3[0] + R[4]*p3[1] + R[5]*p3[2] + t[1];
      Z=R[6]*p3[0] + R[7]*p3[1] + R[8]*p3[2] + t[2];
      /* project 3D point */
      x/=Z;
      y/=Z;

      p2=pts2D[i];
      x-=p2[0];
      y-=p2[1];

      sum+=_SQR(x) + _SQR(y);
    }
  }
  else{
    for(j=0; j<npts; ++j){
      i=ptsidx[j];

      /* transform 3D point */
      p3=pts3D[i];
      x=R[0]*p3[0] + R[1]*p3[1] + R[2]*p3[2] + t[0];
      y=R[3]*p3[0] + R[4]*p3[1] + R[5]*p3[2] + t[1];
      Z=R[6]*p3[0] + R[7]*p3[1] + R[8]*p3[2] + t[2];
      /* project 3D point */
      x/=Z;
      y/=Z;

      p2=pts2D[i];
      x-=p2[0];
      y-=p2[1];

      /* this converts the error to pixels:
      x=K[0]*x + K[1]*y;
      y=         K[4]*y;
      */

      sum+=_SQR(x) + _SQR(y);
    }
  }

  return sum;
}


/*********************************** Robust pose estimation with the LHM algorithm ***************************************/
#include "lqs.h"
#include "ransac.h"
#include "prosac.h"


/* variables used by various estimation routines */
struct LHMPoseData {
  double (*pts2D)[2], (*pts3D)[3], (*wrk)[3];
  double R0[9]; // initial estimate for rotation, used only if supplied externally
  int *inliersidx, numInliers, howtoinit;
};


/* estimate pose. p2=K[R t]*P3, with p2, P3 specified by ptsidx */
static int estLHMPose(double *Rt, int npts, int *ptsidx, void *adata)
{
register int i, j;
int n, verb=0;
struct LHMPoseData *dat=(struct LHMPoseData *)adata;
double (*pts2D)[2]=dat->pts2D, (*pts3D)[3]=dat->pts3D, (*wrk)[3]=dat->wrk;
register double *p2;
double dummy;

  switch(dat->howtoinit){
    case LHM_INITR_WEAKP:
    /* initialize rotation using the weak perspective approximation:
     * consider the set of image points as coplanar 3D points 
     */
    if(ptsidx==NULL){ /* no index, use all points */
      for(i=0; i<npts; ++i){
        p2=pts2D[i];
        wrk[i][0]=p2[0];
        wrk[i][1]=p2[1];
        wrk[i][2]=1.0;
      }
    }
    else{
      for(j=0; j<npts; ++j){
        i=ptsidx[j];

        p2=pts2D[i];
        wrk[i][0]=p2[0];
        wrk[i][1]=p2[1];
        wrk[i][2]=1.0;
      }
    }
    n=sam_absorq(pts3D, wrk, ptsidx, npts, Rt, NULL, &dummy); // t need not be estimated here, use of "dummy" avoids scale warnings
    if(n) return 0; // could not initialize
    break;

    case LHM_INITR_IDENT:
    default:
    /* trivial initialization: eye(3) */
    Rt[0]=1.0; Rt[1]=0.0; Rt[2]=0.0;
    Rt[3]=0.0; Rt[4]=1.0; Rt[5]=0.0;
    Rt[6]=0.0; Rt[7]=0.0; Rt[8]=1.0;
    break;

    case LHM_INITR_SUPPLIED:
    /* use supplied R */
    memcpy(Rt, dat->R0, 9*sizeof(double));
    break;
  }

  n=lhm_estPose(pts2D, pts3D, wrk, ptsidx, npts, Rt, Rt+9, verb);

  return n==POSEST_OK;
}

/* compute the geometric residuals corresponding to Rt as the (squared) reprojection error */
static void poseRT_LHM_ResidualsGeom(double *Rt, int numres, void *adata, double *resid)
{
register int i;
double *R, *t, X, Y, Z, dx, dy;
struct LHMPoseData *dat=(struct LHMPoseData *)adata;
double (*pts2D)[2]=dat->pts2D, (*pts3D)[3]=dat->pts3D;
register double *pt; // pts3D[i] or pts2D[i]

  R=Rt; t=R+9;
  for(i=0; i<numres; ++i){
    /* transform point */
    pt=pts3D[i];
    X=(R[0]*pt[0] + R[1]*pt[1]) + (R[2]*pt[2] + t[0]);
    Y=(R[3]*pt[0] + R[4]*pt[1]) + (R[5]*pt[2] + t[1]);
    Z=(R[6]*pt[0] + R[7]*pt[1]) + (R[8]*pt[2] + t[2]);

    /* note that K is not needed below since pts2D have been normalized! */
    pt=pts2D[i];
    dx=pt[0] - X/Z;
    dy=pt[1] - Y/Z;

    resid[i]=_SQR(dx) + _SQR(dy);
  }
}

/* Robust pose estimation from "nmatches" matched 2D-3D point correspondences, possibly
 * including outliers, using the LHM algorithm. "pts2D", "pts3D" contain the matched 2D-3D
 * point coordinates, "inlPcent" is the expected percentage of inliers (>=0.5),
 * "p" contains the estimated parameters upon return, idxOutliers" points to sufficiently
 * large memory which upon return is set to the indices of the detected outlying points (pass
 * NULL if don't care), "nbOutliers" contains the number of outliers, "verbose" specifies the
 * verbosity level
 *
 * Returns 1 in case of error, 0 if successfull
 *
 */
int posestRT_LHM(double (*pts2D)[2], double (*pts3D)[3], int nmatches, double inlPcent, double K[9],
                 int howtoinit, double p[NUM_RTPARAMS], int *idxOutliers, int *nbOutliers, int verbose)
{
register int i;
const int min_nmatches=3;
int ret;
struct LHMPoseData dat;
double Rt[9+3]; // layout is R (3x3), t (1x3)
double (*wrk)[3], K1[9], *p2;

  if(nmatches<min_nmatches) return 1;  // too few matches

  /* invert K in K1 */
  K1[0]=1.0/K[0]; K1[1]=-K[1]/K[0]/K[4]; K1[2]=(K[1]*K[5]-K[2]*K[4])/K[0]/K[4];
  K1[3]=0.0;      K1[4]=1.0/K[4];        K1[5]=-K[5]/K[4];
  //K1[6]=0.0;      K1[7]=0.0;            K1[8]=1.0; // last row of K assumed [0 0 1]

  /* use K^-1 to convert the supplied pixel coordinates to normalized ones.
   * Note that the input 2D points are overwritten (and restored before return)
   */
  for(i=0; i<nmatches; ++i){
    p2=pts2D[i];
    p2[0]=K1[0]*p2[0] + K1[1]*p2[1] + K1[2];
    p2[1]=              K1[4]*p2[1] + K1[5]; // K1[3]==0
  }

  /* working memory */
  if((wrk=(double (*)[3])malloc(nmatches*sizeof(double[3])))==NULL){
    fprintf(stderr, "memory allocation request failed in posestRT_LHM()\n");
    exit(1);
  }

  dat.pts2D=pts2D; dat.pts3D=pts3D;
  dat.wrk=wrk;
  dat.inliersidx=NULL;
  dat.howtoinit=howtoinit;
  if(howtoinit==LHM_INITR_SUPPLIED) sam_vec2rotmat(p, dat.R0); // store supplied initial R

  if(inlPcent!=1.0){ // robust version
    register int j;
    double gate=2.0, premResid=-1.0, outlierThresh=LHM_RANSAC_OUTL_THRESH; // outlierThresh: distance threshold for outliers (RANSAC/PROSAC only)
    int *outliersMap, **sets=NULL, nbSets=0;
    const int isSqr=1, maxNbSol=1, nparams=9+3;
    int  (*estimator)(double *Rt, int npts, int *ptsidx, void *adata);
    void (*residuals)(double *Rt, int numres, void *adata, double *resid);

    estimator=estLHMPose;
    residuals=poseRT_LHM_ResidualsGeom;

    if(!(outliersMap=(int *)malloc(nmatches*sizeof(int)))){
      fprintf(stderr, "Error: not enough memory for 'outliersMap' in posestRT_LHM()\n");
      exit(1);
    }
    verbose=verbose>1;

#if USE_LQS_FIT==1
    j=lqsfit(nmatches, min_nmatches, sets, nbSets, residuals, estimator,
              isSqr, verbose, maxNbSol, gate, premResid, nparams, inlPcent, (void *)&dat,
              Rt, NULL, outliersMap, nbOutliers, &outlierThresh);
#elif (USE_RANSAC_FIT==1) || (USE_PROSAC_FIT==1)
    gate=premResid=0; /* -Wall */

#if USE_RANSAC_FIT==1
    j=ransacfit(nmatches, min_nmatches, sets, nbSets, residuals, estimator,
              isSqr, verbose, maxNbSol, outlierThresh, 0, nparams, inlPcent, (void *)&dat,
              Rt, NULL, outliersMap, nbOutliers);
#else
    j=prosacfit(nmatches, min_nmatches, residuals, estimator,
              isSqr, verbose, maxNbSol, outlierThresh, 0, nparams, inlPcent, (void *)&dat,
              Rt, NULL, outliersMap, nbOutliers);
#endif
#endif /* USE_LQS_FIT */

    if(verbose){
      fprintf(stderr, "Outlier threshold: %g\n", outlierThresh);
      fprintf(stderr, "posestRT_LHM(): robust fit returned %d, %d outliers [out of %d]\n", j, *nbOutliers, nmatches);
    }

    if(sets) lqs_freesets(sets);

    dat.numInliers=nmatches - *nbOutliers;
    if(j!=0){
      dat.inliersidx=(int *)malloc(dat.numInliers*sizeof(int));
      if(!dat.inliersidx){
        fprintf(stderr, "Error: not enough memory for 'dat.inliersidx' in posestRT_LHM()\n");
        exit(1);
      }

      for(i=j=0; i<nmatches; ++i)
        if(!outliersMap[i]) dat.inliersidx[j++]=i;

      /* LS estimation on inliers */
      estimator(Rt, dat.numInliers, dat.inliersidx, (void *)&dat);

      /* expose outliers */
      if(idxOutliers!=NULL)
        for(i=j=0; i<nmatches; ++i)
          if(outliersMap[i]) idxOutliers[j++]=i;

#if 0
      if(verbose){
        fputs("Outliers: ", stderr);
        for(i=j=0; i<nmatches; ++i)
          if(outliersMap[i]) fprintf(stderr, "%d ", i);
        fputc('\n', stderr);
    }
#endif

      ret=POSEST_OK;

#if 0
      /* include the following code fragment to print the matching 3D-2D point pairs found to be inlying */
      for(i=0; i<dat.numInliers; ++i){
        printf("%.6lf %.6lf %.6lf   %.4lf %.4lf\n", pts3D[dat.inliersidx[i]][0], pts3D[dat.inliersidx[i]][1], pts3D[dat.inliersidx[i]][2], 
                                                    pts2D[dat.inliersidx[i]][0], pts2D[dat.inliersidx[i]][1]);
      }
#endif

    }
    else{ /* robust fit failed */
      memset(Rt, 0, nparams*sizeof(double));
      *nbOutliers=nmatches;
      dat.numInliers=0;
      ret=POSEST_ERR;
    }

    if(dat.inliersidx) free(dat.inliersidx);
    free(outliersMap);
  }
  else{ // use all matches
    *nbOutliers=0;
    //ret=lhm_estPose(pts2D, pts3D, wrk, NULL, nmatches, Rt, Rt+9, verbose);
    ret=(estLHMPose(Rt, nmatches, NULL, (void *)&dat)==1)? POSEST_OK : POSEST_ERR; // as above plus R initialization
  }

#if 0
  if(sam_det3x3(Rt)<0.0) // R is a reflection, negate Rt
    for(i=0; i<12; ++i) Rt[i]=-Rt[i];
#endif

  /* convert to rotation vector and copy */
  sam_rotmat2vec(Rt, p); // R
  p[3]=Rt[9]; p[4]=Rt[10]; p[5]=Rt[11]; // t

  /* undo the normalization with K^-1 and restore input points */
  for(i=0; i<nmatches; ++i){
    p2=pts2D[i];
    p2[0]=K[0]*p2[0] + K[1]*p2[1] + K[2];
    p2[1]=             K[4]*p2[1] + K[5]; // K[3]==0
  }

  free(wrk);

  return ret;
}

/*********************************** Robust pose estimation with the LHM algorithm ***************************************/

#if 0

main()
{
clock_t start_time, end_time;

#if 0
#define N 25
double pts2D[N][2]={
  -695.433721377876, 1105.9792763002,
  206.524561488808, 666.245576686066,
  206.458530505325, 917.498753384862,
  228.891354266207, 489.747375125415,
  -798.53305274582, 1663.78347211428,
  -536.235256193386, 923.334023406567,
  -167.273389130415, 201.303060501723,
  105.985487441951, 349.326654443714,
  22.326819944374, 482.104332456162,
  126.15302427832, 902.198507139287,
  -44.6299914800215, 392.744863236346,
  -26.9224900580949, 405.723119735403,
  -30194.2108405043, 10681.3215715422,
  -785.692461363758, 824.649666201522,
  133.155015689277, 227.772052062493,
  39.3948648858497, 102.023493766942,
  95.1366349699194, 477.304114250089,
  142.125057435593, 737.190443614491,
  -403.188742015687, 474.272676732881,
  -162.183920738888, 222.638796949274,
  -115.049703946065, 92.3402791937369,
  -288.603106211099, 580.862649144927,
  192.009174920585, 883.094665293152,
  -450.283861419754, 117.802601486909,
  -426.03173907901, 784.028945507316,
};

double pts3D[N][3]={
  {68.1335, 66.5823, 13.4718},
  {2.2493, 26.2199, 11.6515},
  {6.9318, 85.2930, 18.0331},
  {3.2419, 73.3926, 53.6517},
  {27.6030, 36.8458, 1.2886},
  {88.9206, 86.6021, 25.4247},
  {56.9481, 15.9265, 59.4364},
  {33.1100, 65.8613, 86.3634},
  {56.7623, 98.0481, 79.1832},
  {15.2594, 83.3027, 19.1863},
  {63.8987, 66.9000, 77.2088},
  {37.9818, 44.1585, 48.3060},
  {60.8106, 17.5996, 0.2026},
  {79.0224, 51.3609, 21.3229},
  {10.3450, 15.7337, 40.7515},
  {40.7757, 5.2693, 94.1815},
  {14.9972, 38.4374, 31.1059},
  {16.8534, 89.6648, 32.2724},
  {73.3996, 41.0904, 39.9794},
  {50.5522, 16.9306, 52.4745},
  {64.1203, 1.6197, 83.6852},
  {80.3462, 69.7785, 46.1888},
  {8.2613, 82.0717, 19.3020},
  {44.5355, 1.2958, 30.8742},
  {87.5351, 83.5259, 33.3095},
};

double K[9]={
  394.2054307, 0, 256,
  0, 394.2054307, 192,
  0, 0, 1};
#endif

#if 1
// outliers included
#define N 195
double pts2D[N][2]={
 {3.2042383e+02, 6.2678540e+02},
 {3.1401236e+02, 6.3092230e+02},
 {3.1171075e+02, 6.4624512e+02},
 {3.3232681e+02, 6.4502295e+02},
 {3.3020486e+02, 6.9652472e+02},
 {3.2563989e+02, 6.7685614e+02},
 {2.9381891e+02, 6.4272961e+02},
 {3.5625665e+02, 6.8932257e+02},
 {3.4925406e+02, 6.6554376e+02},
 {3.2724851e+02, 6.1521417e+02},
 {3.3409604e+02, 6.7707788e+02},
 {3.2529886e+02, 6.5170502e+02},
 {3.4318426e+02, 6.6778339e+02},
 {3.3701782e+02, 6.5472070e+02},
 {3.0310031e+02, 6.3624560e+02},
 {3.0657693e+02, 6.7567895e+02},
 {3.1699576e+02, 6.9266046e+02},
 {3.3685211e+02, 6.0570935e+02},
 {3.4216580e+02, 6.1189099e+02},
 {3.2025241e+02, 6.1559100e+02},
 {3.1345343e+02, 6.1299756e+02},
 {3.3054767e+02, 6.5632629e+02},
 {3.1545798e+02, 6.6174536e+02},
 {3.0293677e+02, 6.4467120e+02},
 {3.3008374e+02, 6.0739557e+02},
 {3.3409604e+02, 6.7707788e+02},
 {3.0862061e+02, 7.0099652e+02},
 {2.9381891e+02, 6.4272961e+02},
 {3.1975537e+02, 6.0522882e+02},
 {3.0891788e+02, 6.3529474e+02},
 {3.0889908e+02, 5.5874518e+02},
 {2.9868109e+02, 6.3358142e+02},
 {3.3271948e+02, 6.1943835e+02},
 {3.4736084e+02, 6.5829608e+02},
 {3.0703326e+02, 7.1017310e+02},
 {3.4853845e+02, 5.6752997e+02},
 {3.1336847e+02, 6.9782343e+02},
 {3.2966022e+02, 6.8083484e+02},
 {3.3964591e+02, 6.8725427e+02},
 {3.1806491e+02, 6.9879712e+02},
 {3.5404248e+02, 5.6730499e+02},
 {3.2905530e+02, 7.1225915e+02},
 {3.2782172e+02, 5.6464624e+02},
 {3.3046835e+02, 6.6444000e+02},
 {3.1655829e+02, 5.6227228e+02},
 {3.3055633e+02, 6.5135004e+02},
 {3.4966153e+02, 6.7499005e+02},
 {3.1919693e+02, 6.3709265e+02},
 {2.9120688e+02, 5.6582013e+02},
 {3.0922400e+02, 5.4196887e+02},
 {3.1744049e+02, 5.7249072e+02},
 {3.4854700e+02, 6.3645520e+02},
 {3.0710266e+02, 6.8067859e+02},
 {3.1069778e+02, 6.7161621e+02},
 {2.8079245e+02, 6.9544336e+02},
 {2.8811664e+02, 6.4389496e+02},
 {3.5671747e+02, 5.4821320e+02},
 {2.9905328e+02, 6.6020221e+02},
 {3.5177985e+02, 6.4100006e+02},
 {3.6234598e+02, 7.0864594e+02},
 {3.3275732e+02, 5.7158179e+02},
 {3.0478000e+02, 5.7235840e+02},
 {3.3818259e+02, 6.7807660e+02},
 {2.9742987e+02, 5.4604706e+02},
 {3.4702957e+02, 5.6159509e+02},
 {3.0819955e+02, 5.4896008e+02},
 {3.3902472e+02, 6.7255835e+02},
 {2.8545526e+02, 6.9836090e+02},
 {3.7468390e+02, 5.5854596e+02},
 {3.2247501e+02, 6.6835083e+02},
 {2.8772791e+02, 6.3964679e+02},
 {3.4348346e+02, 5.5340656e+02},
 {3.0552896e+02, 6.3926263e+02},
 {3.4112259e+02, 6.8338068e+02},
 {3.3632153e+02, 6.6132684e+02},
 {3.3392188e+02, 6.4179590e+02},
 {3.0062262e+02, 6.7347705e+02},
 {2.9929916e+02, 5.4101776e+02},
 {3.5742258e+02, 5.5263238e+02},
 {3.0200619e+02, 6.9704858e+02},
 {3.1131580e+02, 7.0798523e+02},
 {3.4671033e+02, 6.9059589e+02},
 {3.2782172e+02, 5.6464624e+02},
 {3.0982031e+02, 5.6740216e+02},
 {3.4239987e+02, 6.3997351e+02},
 {2.7749432e+02, 6.1947479e+02},
 {3.3860113e+02, 6.2327942e+02},
 {2.9868109e+02, 6.3358142e+02},
 {2.9442288e+02, 6.4995001e+02},
 {3.0560556e+02, 6.1146405e+02},
 {3.3447528e+02, 6.5814612e+02},
 {3.2195261e+02, 5.4570483e+02},
 {3.0819955e+02, 5.4896008e+02},
 {3.1345496e+02, 6.8214825e+02},
 {2.8255286e+02, 6.3410712e+02},
 {3.0981473e+02, 6.6740155e+02},
 {3.0560556e+02, 6.1146405e+02},
 {3.4352942e+02, 6.7876251e+02},
 {2.8811664e+02, 6.4389496e+02},
 {3.7260306e+02, 5.6670593e+02},
 {3.4662683e+02, 5.8306384e+02},
 {3.1964047e+02, 6.2182214e+02},
 {2.7961423e+02, 6.1221368e+02},
 {3.2150476e+02, 6.8882629e+02},
 {2.9516232e+02, 6.3807636e+02},
 {3.0215289e+02, 6.7969733e+02},
 {3.3672607e+02, 6.1477490e+02},
 {3.1395694e+02, 6.3588501e+02},
 {2.7663184e+02, 5.8517346e+02},
 {3.4043384e+02, 7.0203040e+02},
 {3.0724634e+02, 6.5447119e+02},
 {3.6636533e+02, 5.3574853e+02},
 {2.8093637e+02, 6.8819153e+02},
 {3.5571454e+02, 6.5088477e+02},
 {3.0294681e+02, 6.0510394e+02},
 {3.2905530e+02, 7.1225915e+02},
 {3.5385602e+02, 5.6320612e+02},
 {3.5365463e+02, 5.3879767e+02},
 {3.3816861e+02, 6.6522333e+02},
 {3.6127588e+02, 6.8524799e+02},
 {3.5350275e+02, 6.4497766e+02},
 {3.2286850e+02, 6.8203656e+02},
 {2.9281427e+02, 6.7437628e+02},
 {2.8741901e+02, 5.7095758e+02},
 {3.0237738e+02, 5.5141748e+02},
 {3.2660388e+02, 5.4013934e+02},
 {2.8369141e+02, 6.1055249e+02},
 {3.4456964e+02, 6.0082080e+02},
 {3.4141470e+02, 5.6941278e+02},
 {3.4991058e+02, 6.0041425e+02},
 {3.3745078e+02, 6.4732587e+02},
 {3.1921762e+02, 5.8370081e+02},
 {2.8622006e+02, 5.3373578e+02},
 {3.1529044e+02, 6.5607520e+02},
 {3.2601947e+02, 6.4342548e+02},
 {3.7881073e+02, 6.5801526e+02},
 {3.0784564e+02, 6.8975128e+02},
 {3.2171594e+02, 5.5308154e+02},
 {3.2196332e+02, 6.3178729e+02},
 {3.4141470e+02, 5.6941278e+02},
 {3.2735788e+02, 5.9574438e+02},
 {3.2130029e+02, 7.2458795e+02},
 {3.3271948e+02, 6.1943835e+02},
 {2.8910123e+02, 5.8542786e+02},
 {3.5825842e+02, 5.7821906e+02},
 {3.0660849e+02, 7.0571149e+02},
 {3.7294101e+02, 6.4086340e+02},
 {2.7827399e+02, 6.2680115e+02},
 {3.5432459e+02, 6.7079297e+02},
 {3.0887598e+02, 6.1638953e+02},
 {3.5667783e+02, 6.6229138e+02},
 {3.0552896e+02, 6.3926263e+02},
 {3.6397946e+02, 6.7287109e+02},
 {3.4699707e+02, 7.0789081e+02},
 {2.9258255e+02, 7.0148499e+02},
 {3.2945776e+02, 6.7224768e+02},
 {2.9521484e+02, 6.5433508e+02},
 {3.6158402e+02, 6.5581457e+02},
 {3.5108124e+02, 6.8881219e+02},
 {3.0147253e+02, 6.4020117e+02},
 {3.0101294e+02, 7.0917554e+02},
 {2.9905328e+02, 6.6020221e+02},
 {3.6754645e+02, 6.9929956e+02},
 {2.7977234e+02, 6.4974133e+02},
 {3.3507929e+02, 5.6304645e+02},
 {3.6181000e+02, 5.2883441e+02},
 {3.6094900e+02, 7.1197308e+02},
 {3.1744049e+02, 5.7249072e+02},
 {2.9460745e+02, 6.6891455e+02},
 {3.3685211e+02, 6.0570935e+02},
 {3.6854752e+02, 6.4123303e+02},
 {2.9682944e+02, 6.6570093e+02},
 {3.0487048e+02, 6.2119745e+02},
 {3.5523142e+02, 6.0085962e+02},
 {2.9322839e+02, 6.9692511e+02},
 {3.3818259e+02, 6.7807660e+02},
 {3.1493994e+02, 7.1734668e+02},
 {3.5617236e+02, 7.1616809e+02},
 {3.6284427e+02, 5.6077014e+02},
 {3.5667783e+02, 6.6229138e+02},
 {3.2601947e+02, 6.4342548e+02},
 {3.6284427e+02, 5.6077014e+02},
 {3.2221655e+02, 5.7764081e+02},
 {3.3020126e+02, 6.8505212e+02},
 {3.5373328e+02, 6.5507617e+02},
 {3.3793070e+02, 6.0976270e+02},
 {3.1655829e+02, 5.6227228e+02},
 {3.0966904e+02, 6.5029047e+02},
 {2.8409195e+02, 6.8052234e+02},
 {3.5089386e+02, 5.8800153e+02},
 {3.1960120e+02, 7.1313361e+02},
 {3.2660388e+02, 5.4013934e+02},
 {3.6430533e+02, 6.4127118e+02},
 {3.1311206e+02, 5.8081616e+02},
 {3.0311478e+02, 5.6341443e+02},
};
double pts3D[N][3]={
 {1.1957120e+00, -4.0135284e+01, 6.3842944e+02},
 {3.3878620e+00, -4.0098549e+01, 6.4183081e+02},
 {5.5019790e+00, -3.1898710e+01, 6.4988635e+02},
 {-2.2910940e+00, -2.8695370e+01, 6.4369269e+02},
 {1.4763700e+00, -3.5323320e+00, 6.6727795e+02},
 {2.0386210e+00, -1.4241404e+01, 6.5966498e+02},
 {1.1362709e+01, -3.8331661e+01, 6.5271008e+02},
 {-8.2470770e+00, -1.5710780e+00, 6.5696802e+02},
 {-7.1297110e+00, -1.4866354e+01, 6.4795648e+02},
 {-2.2594330e+00, -4.4907253e+01, 6.3099609e+02},
 {-9.2456400e-01, -1.2416084e+01, 6.5736536e+02},
 {9.7604000e-01, -2.6798077e+01, 6.4865979e+02},
 {-4.0575720e+00, -1.3533005e+01, 6.5108826e+02},
 {-3.1646470e+00, -2.2976480e+01, 6.4649091e+02},
 {7.5149790e+00, -3.9426750e+01, 6.4741046e+02},
 {8.3003150e+00, -1.9550343e+01, 6.6413275e+02},
 {5.9321090e+00, -8.4653140e+00, 6.6924963e+02},
 {-6.0604580e+00, -4.6894360e+01, 6.2486786e+02},
 {-7.6545100e+00, -4.2842793e+01, 6.2598999e+02},
 {2.8097100e-01, -4.5809998e+01, 6.3358704e+02},
 {2.6634710e+00, -4.8491188e+01, 6.3427844e+02},
 {-1.0932860e+00, -2.3556313e+01, 6.4898712e+02},
 {5.3827720e+00, -2.3373871e+01, 6.5552576e+02},
 {8.6336950e+00, -3.4949787e+01, 6.5131067e+02},
 {-3.6308010e+00, -4.7933022e+01, 6.2736688e+02},
 {-1.0672540e+00, -1.2380124e+01, 6.5730420e+02},
 {9.5325000e+00, -5.6595310e+00, 6.7569312e+02},
 {1.1383117e+01, -3.8392998e+01, 6.5269958e+02},
 {-9.4150400e-01, -5.1905003e+01, 6.2912494e+02},
 {5.0987200e+00, -3.8941105e+01, 6.4516846e+02},
 {1.9170237e+01, -5.3881077e+01, 6.1388287e+02},
 {9.2712120e+00, -4.1717407e+01, 6.4748785e+02},
 {-4.3291960e+00, -4.1676567e+01, 6.3189563e+02},
 {-7.5228070e+00, -1.9586622e+01, 6.4494653e+02},
 {1.0671677e+01, -1.6058060e+00, 6.7976758e+02},
 {7.4804500e-01, -4.7057396e+01, 6.0617200e+02},
 {7.3001540e+00, -7.3790740e+00, 6.7184998e+02},
 {6.8364100e-01, -1.1380725e+01, 6.6042151e+02},
 {-3.0679930e+00, -6.7632110e+00, 6.6029211e+02},
 {6.0250100e+00, -5.2038570e+00, 6.7168787e+02},
 {-9.7998400e-01, -4.5882355e+01, 6.0484637e+02},
 {3.0097210e+00, 4.3225030e+00, 6.7468842e+02},
 {1.0281110e+01, -5.1048782e+01, 6.1015759e+02},
 {-6.7795200e-01, -1.9781914e+01, 6.5271545e+02},
 {1.4269083e+01, -5.3250347e+01, 6.1313074e+02},
 {-9.7375100e-01, -2.5889183e+01, 6.4664062e+02},
 {-6.9334240e+00, -1.0328885e+01, 6.5212207e+02},
 {1.0099120e+00, -3.6434826e+01, 6.4342108e+02},
 {2.2558277e+01, -6.0800293e+01, 6.2031982e+02},
 {2.9096621e+01, -4.8981224e+01, 6.0709302e+02},
 {8.0551770e+00, -5.5851124e+01, 6.1675989e+02},
 {-8.2943220e+00, -2.9385485e+01, 6.3474744e+02},
 {8.8630780e+00, -1.6804825e+01, 6.6651141e+02},
 {8.0687850e+00, -1.9481894e+01, 6.6130585e+02},
 {1.1421287e+01, -3.8222885e+01, 5.9463593e+02},
 {1.2848510e+01, -3.8929886e+01, 6.5485400e+02},
 {8.7680860e+00, -4.0396423e+01, 5.9648712e+02},
 {1.0433314e+01, -2.8556152e+01, 6.5931555e+02},
 {-1.0111902e+01, -2.7132530e+01, 6.3594586e+02},
 {2.2855070e+00, 5.4050310e+00, 6.3514258e+02},
 {4.0148520e+00, -5.1808365e+01, 6.1176001e+02},
 {1.3464607e+01, -5.8495922e+01, 6.1853131e+02},
 {-2.3341520e+00, -1.0855775e+01, 6.5674518e+02},
 {3.1185076e+01, -5.3111397e+01, 6.1164569e+02},
 {4.8998020e+00, -4.5420170e+01, 6.0392499e+02},
 {2.5299435e+01, -5.1447960e+01, 6.1027289e+02},
 {-3.1303050e+00, -1.3294635e+01, 6.5394397e+02},
 {-2.3143650e+00, -5.3713470e+01, 6.1453412e+02},
 {3.7666027e+01, -4.2763802e+01, 6.0561438e+02},
 {4.6379997e+01, -4.1466412e+01, 6.1695435e+02},
 {1.3252425e+01, -1.4388830e+01, 6.0156372e+02},
 {1.0118769e+01, -4.4295284e+01, 6.0262604e+02},
 {3.5615414e+01, -5.0588364e+01, 6.0861035e+02},
 {-3.5745570e+00, -7.4441700e+00, 6.5814667e+02},
 {3.1989712e+01, -2.8460140e+00, 6.7153015e+02},
 {-9.9710200e+00, -4.6604023e+01, 6.2059448e+02},
 {2.7775215e+01, -5.4107014e+01, 6.3119751e+02},
 {3.3281689e+01, -5.1447258e+01, 6.0950415e+02},
 {3.6114918e+01, -3.2863544e+01, 6.3979736e+02},
 {1.1321325e+01, -9.9219590e+00, 6.7541418e+02},
 {7.6703030e+00, -3.1059110e+00, 6.7776050e+02},
 {5.0912000e-02, -6.8960010e+00, 6.2760290e+02},
 {1.0024904e+01, -5.0828506e+01, 6.1047052e+02},
 {1.3072802e+01, -5.6382530e+01, 6.1691980e+02},
 {-6.1620400e+00, -2.9161976e+01, 6.3802771e+02},
 {-1.1276275e+01, -3.1510538e+01, 6.2154303e+02},
 {-6.6372910e+00, -3.8658726e+01, 6.3210767e+02},
 {9.3478290e+00, -4.1501156e+01, 6.4744202e+02},
 {1.1316812e+01, -3.5610653e+01, 6.5550568e+02},
 {3.8034115e+01, 8.8882760e+00, 6.4447766e+02},
 {-5.2486460e+00, -5.4970989e+01, 6.1677600e+02},
 {2.2966030e+01, -4.6906868e+01, 6.0508545e+02},
 {2.5456758e+01, -5.1552525e+01, 6.1017731e+02},
 {7.8151080e+00, -1.3251060e+01, 6.6570575e+02},
 {2.9001411e+01, 2.1328740e+00, 6.7953357e+02},
 {6.9999430e+00, -2.2967831e+01, 6.5968970e+02},
 {2.2393669e+01, 4.1042830e+00, 6.2230395e+02},
 {-4.4789460e+00, -9.8277090e+00, 6.5551910e+02},
 {3.0729513e+01, -5.0094494e+01, 6.0770050e+02},
 {1.3673960e+00, -5.4183334e+01, 6.2957483e+02},
 {1.7523211e+01, -5.4392067e+01, 6.1439587e+02},
 {1.4277203e+01, -6.2339240e+00, 6.0881006e+02},
 {1.6102396e+01, 2.9714685e+01, 6.3834515e+02},
 {3.2134820e+00, -1.0399455e+01, 6.6622601e+02},
 {1.0645333e+01, -4.0495251e+01, 6.5042822e+02},
 {1.0332906e+01, -1.3351426e+01, 6.7001123e+02},
 {-5.4555070e+00, -4.2678547e+01, 6.2881134e+02},
 {6.8105000e+00, -2.6656120e+00, 6.2233728e+02},
 {2.2926653e+01, -2.4212904e+01, 6.0096979e+02},
 {-6.5287300e-01, 2.6047660e+00, 6.6798163e+02},
 {6.9835700e+00, -2.9988213e+01, 6.5422644e+02},
 {9.0851740e+00, -1.9289347e+01, 6.0444971e+02},
 {1.0846106e+01, -3.7523087e+01, 5.9382123e+02},
 {-9.7982370e+00, -2.0308395e+01, 6.3961176e+02},
 {3.9277100e+00, -5.6242554e+01, 6.3340601e+02},
 {4.6235107e+01, 8.0706940e+00, 6.6273676e+02},
 {4.5181782e+01, -4.0432648e+01, 6.1675165e+02},
 {1.9415985e+01, -8.0555570e+00, 6.0908124e+02},
 {3.2495220e+01, -4.4451950e+00, 6.6964008e+02},
 {-1.0313950e+01, -2.3967320e+00, 6.5367535e+02},
 {-1.1291746e+01, -2.3232410e+01, 6.3652106e+02},
 {2.7293680e+00, -1.2806765e+01, 6.6236322e+02},
 {3.9711050e+00, -3.0956135e+01, 6.4998456e+02},
 {1.9821436e+01, -6.2529953e+01, 6.2405237e+02},
 {2.6377678e+01, -5.3777893e+01, 6.1305920e+02},
 {1.6825657e+01, -6.5533530e+00, 6.8354218e+02},
 {2.3650757e+01, -4.7432388e+01, 6.4372516e+02},
 {2.2417715e+01, -2.3335573e+01, 6.0073737e+02},
 {1.9240930e+00, -4.9489803e+01, 6.0890869e+02},
 {-1.2231112e+01, -4.7714554e+01, 6.1879749e+02},
 {3.5954334e+01, -4.0486240e+01, 6.3172400e+02},
 {6.5514160e+00, -4.3820122e+01, 6.0206018e+02},
 {3.2140144e+01, -2.9822987e+01, 6.4643866e+02},
 {1.9048473e+01, -6.2032310e+01, 6.2376312e+02},
 {2.9387342e+01, 7.8378310e+00, 6.3329614e+02},
 {2.5002831e+01, -5.2604542e+01, 6.3747174e+02},
 {3.3548958e+01, -1.8076929e+01, 6.5529144e+02},
 {1.7383673e+01, -4.9904396e+01, 6.0797437e+02},
 {8.1414300e+00, -5.3559734e+01, 6.1412244e+02},
 {-2.4197580e+00, 6.2627680e+00, 6.4226184e+02},
 {1.0841356e+01, -1.2991453e+01, 6.6957245e+02},
 {-1.4121530e+00, -1.0254523e+01, 6.5844190e+02},
 {3.8409225e+01, -5.5149193e+01, 6.1449878e+02},
 {3.6448238e+01, -8.0971170e+00, 6.3010394e+02},
 {4.5190517e+01, 1.1090660e+01, 6.6708948e+02},
 {1.0295213e+01, -4.2058170e+00, 6.7786359e+02},
 {3.9334061e+01, -4.3772991e+01, 6.2571173e+02},
 {1.7851017e+01, -5.8081432e+01, 6.4056799e+02},
 {-8.3638380e+00, -1.0868964e+01, 6.4893860e+02},
 {4.4318700e+00, -4.8087921e+01, 6.3693207e+02},
 {3.7786442e+01, -5.2254932e+01, 6.1027515e+02},
 {1.9024252e+01, 1.5508677e+01, 6.2696808e+02},
 {2.5287127e+01, -6.2094269e+01, 6.2286047e+02},
 {4.2373741e+01, -1.6968431e+01, 6.3009082e+02},
 {3.6320637e+01, -5.6797451e+01, 6.1680688e+02},
 {-8.2500000e-01, -1.5181301e+01, 6.5518310e+02},
 {1.0868499e+01, -3.2948746e+01, 6.5769318e+02},
 {-8.5824790e+00, 1.6066044e+01, 6.6912390e+02},
 {-6.5053190e+00, -2.9361480e+00, 6.5814136e+02},
 {2.9709064e+01, 4.4977720e+00, 6.8299768e+02},
 {1.3880058e+01, -2.5750660e+00, 6.8142120e+02},
 {-1.1112975e+01, -2.9244797e+01, 6.2164417e+02},
 {2.1478964e+01, -3.9749640e+00, 6.8460168e+02},
 {1.7272957e+01, -4.1991856e+01, 5.9900256e+02},
 {1.9197670e+01, -4.9596207e+01, 6.0889838e+02},
 {2.4219570e+01, -1.2366125e+01, 6.7444104e+02},
 {-8.9452970e+00, -4.5092960e+01, 6.0573499e+02},
 {7.6647590e+00, -5.6033855e+01, 6.1729669e+02},
 {1.2149939e+01, -2.5788202e+01, 6.6449359e+02},
 {1.2440146e+01, -4.4262848e+01, 6.0197162e+02},
 {2.6771618e+01, 3.3123295e+01, 6.4756818e+02},
 {2.4033689e+01, -3.0445013e+01, 6.5757544e+02},
 {-8.4404750e+00, 5.9622330e+00, 6.6181116e+02},
 {2.2821976e+01, 1.8610540e+00, 6.2030457e+02},
 {1.4668765e+01, -1.1635214e+01, 6.7754181e+02},
 {9.1276120e+00, 9.6380100e-01, 6.2330194e+02},
 {2.1705940e+01, -2.1195372e+01, 6.0172687e+02},
 {2.3463015e+01, -4.3640938e+01, 6.0088483e+02},
 {1.2552464e+01, -4.5547527e+01, 6.0325647e+02},
 {-9.5817620e+00, 1.7096622e+01, 6.6067895e+02},
 {3.6740570e+00, -1.1430763e+01, 6.1850842e+02},
 {3.1907507e+01, -3.0536860e+01, 6.4655194e+02},
 {-1.1165501e+01, 2.0843527e+01, 6.6509375e+02},
 {4.8502000e-01, -9.5868400e+00, 6.6172052e+02},
 {-9.3589000e+00, -1.9583582e+01, 6.4185590e+02},
 {-6.7344600e+00, -4.5171566e+01, 6.2610541e+02},
 {1.8780005e+01, 1.9959452e+01, 6.3005810e+02},
 {2.6827921e+01, -2.9166050e+00, 6.7853418e+02},
 {-4.1270010e+00, 2.3161541e+01, 6.5979340e+02},
 {4.2250782e+01, 2.2553558e+01, 6.6037756e+02},
 {7.6629320e+00, -5.2160168e+01, 6.1213672e+02},
 {3.6158085e+01, -4.9726440e+01, 6.0656281e+02},
 {2.1011604e+01, -1.5177703e+01, 6.0566943e+02},
 {-1.5076400e+00, -3.4270485e+01, 6.0528302e+02},
 {1.7737875e+01, -3.7233231e+01, 5.9281580e+02},
};
#endif

double K[9]={
  3252.33, 0, 650.247, 0, 3256.29, 414.877, 0, 0, 1
};
double rt[6];
int noutl;


  start_time=clock();
  posestRT_LHM(pts2D, pts3D, N, 0.4, K, LHM_INITR_WEAKP, rt, NULL, &noutl, 1);
  end_time=clock();

  fprintf(stdout, "\nEstimated motion ([rv t]) [%d outliers, %.2lf%%]\n", noutl, (double)(100.0*noutl)/N);
  printf("%g %g %g %g %g %g\n", rt[0], rt[1], rt[2], rt[3], rt[4], rt[5]);

  fprintf(stdout, "\nElapsed time: %.2lf seconds, %.2lf msecs\n", ((double) (end_time - start_time)) / CLOCKS_PER_SEC,
                                ((double) (end_time - start_time)) / (CLOCKS_PER_SEC/1000.0));

  /*
  {
  double R[9], t[3];
  posestLHM(pts2D, pts3D, N, K, R, t);
  printf("%g %g %g\n", R[0], R[1], R[2]);
  printf("%g %g %g\n", R[3], R[4], R[5]);
  printf("%g %g %g\n", R[6], R[7], R[8]);
  printf("\n%g %g %g\n", t[0], t[1], t[2]);
  }
  */
}
#endif
