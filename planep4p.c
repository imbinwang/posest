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

#include "compiler.h"

#include "sam.h"

#include "polysolve.h"
#include "util.h"
#include "p3p.h"
#include "planep4p.h"

#include "svd3.h"


/******************************* F & B coplanar P4P start *******************************/

/* P4P for coplanar points. see Appendix B in Fischler and Bolles */

#define _CROSSPROD(v, x, y){ (v)[0]=(x)[1]*(y)[2] - (x)[2]*(y)[1]; (v)[1]=(x)[2]*(y)[0] - (x)[0]*(y)[2]; (v)[2]=(x)[0]*(y)[1] - (x)[1]*(y)[0]; }

/* compute the rotation and translation aligning
 * the plane of three arbitrary points with the XY
 * plane. This is equivalent to transforming the
 * plane of the points to Z=0.
 * see http://www.mathworks.ch/matlabcentral/newsreader/view_thread/306882
 */
static void alignXY(double p1[3], double p2[3], double p3[3], double R[9], double t[3])
{
double n[3], mag1, th;
double tmp1[3], tmp2[3];
double ax[3];
static const double nz[3]={0.0, 0.0, 1.0};

  /* find the normal to the plane of p1,p2,p3 */
  tmp1[0]=p2[0]-p1[0]; tmp1[1]=p2[1]-p1[1]; tmp1[2]=p2[2]-p1[2];
  tmp2[0]=p3[0]-p1[0]; tmp2[1]=p3[1]-p1[1]; tmp2[2]=p3[2]-p1[2];
  _CROSSPROD(n, tmp1, tmp2);
  mag1=1.0/sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0]*=mag1; n[1]*=mag1; n[2]*=mag1;

  /* find the dihedral angle between the plane of p1,p2,p3 and Z=0:
   * n'*[0 0 1]=n(3)
   */
  th=acos(n[2]);

  /* calculate R for axis n x [0 0 1] and angle th */
  _CROSSPROD(ax, n, nz);
  mag1=ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2];
  if(mag1>0.0)
    sam_axangle2rotmat(ax, th, R);
  else{
    /* zero axis, i.e. plane is already aligned with XY */
    memset(R, 0, 9*sizeof(double));
    R[0]=R[4]=R[8]=1.0;
    t[0]=t[1]=t[2]=0.0;
  }

  /* calculate t as (0, 0, -[R*p1]_z). To bring the plane's
   * origin at p1, t should be -R*p1
   */
  t[0]=0.0; //-(R[0]*p1[0] + R[1]*p1[1] + R[2]*p1[2]);
  t[1]=0.0; //-(R[3]*p1[0] + R[4]*p1[1] + R[5]*p1[2]);
  t[2]=-(R[6]*p1[0] + R[7]*p1[1] + R[8]*p1[2]);
}

#if 0
/* project Q in the plane defined by the three points in pts */
static void projectPtOnPlane(double pts[3*3], double Q[3], double proj[3])
{
double *p1=pts, *p2=pts+3, *p3=pts+6;
double tmp1[3], tmp2[3], n[3], mag;

  /* find the plane normal */
  tmp1[0]=p2[0]-p1[0]; tmp1[1]=p2[1]-p1[1]; tmp1[2]=p2[2]-p1[2];
  tmp2[0]=p3[0]-p1[0]; tmp2[1]=p3[1]-p1[1]; tmp2[2]=p3[2]-p1[2];
  _CROSSPROD(n, tmp1, tmp2);
  mag=1.0/sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0]*=mag; n[1]*=mag; n[2]*=mag;

  /* tmp1=p1-Q */
  tmp1[0]=p1[0]-Q[0]; tmp1[1]=p1[1]-Q[1]; tmp1[2]=p1[2]-Q[2];
  /* (p1-Q)'*n */
  mag=tmp1[0]*n[0] + tmp1[1]*n[1] + tmp1[2]*n[2];
  
  /* proj=Q + (p1-Q)'*n)*n */
  proj[0]=Q[0] + mag*n[0];
  proj[1]=Q[1] + mag*n[1];
  proj[2]=Q[2] + mag*n[2];

  printf("PROJ: %g %g %g\n", proj[0], proj[1], proj[2]);
}
#endif

/* calculate in H the collineation mapping XY-aligned planar object points
 * xM to normalized image plane points K^-1*m. K are the camera intrinsics
 * contained in cp and xM contains the points of M in a coordinate system
 * such that the object plane is parallel to the XY plane, i.e. has equation Z=0.
 * See Fischler & Bolles, p. 395
 *
 * returns 0 if successful, 1 otherwise
 */
static int W2Icollineation(struct p3p_calib_params *cp, double m[4][2], double M[4][3], double xM[4][2], double H[9], double *R, double *t)
{
register int i;
double Ra[9], ta[3]; // transform M to xM: xM=R*M+t (Z coordinates are zero)
double Q[4*3], P[4*3]; // 4x3
double Q1[9], P1[9], r[3], v[3], *pq4;
double w1, w2;
double wQ1t[9];
double Hn[9], K1[9];

  /* Note: m contains image points in pixel coordinates (not normalized).
   * H is computed from xM, m and then left multiplied by K^-1 to
   * correspond to normalized image coordinates, as Fischler & Bolles
   * assume. This has been observed to produce more accurate results
   * compared to the alternative of computing H from normalized image
   * points directly.
   */

  if(!R) R=Ra;
  if(!t) t=ta;

  /* transform object points to Z=0 */
  alignXY(M[0], M[1], M[2], R, t);

  /* compute H */
  /* Q=R*M+t */
  for(i=0; i<4; ++i){
    double *q=M[i];

    Q[i*3+0]=xM[i][0]=R[0]*q[0] + R[1]*q[1] + R[2]*q[2] + t[0];
    Q[i*3+1]=xM[i][1]=R[3]*q[0] + R[4]*q[1] + R[5]*q[2] + t[1];
                    //R[6]*q[0] + R[7]*q[1] + R[8]*q[2] + t[2]; // this should be zero
    Q[i*3+2]=1.0; // homogeneous cordinates

    P[i*3+0]=m[i][0];
    P[i*3+1]=m[i][1];
    P[i*3+2]=1.0; // homogeneous cordinates
  }

  /* notice that only the first 3 rows are used in the inversions below! */
  if(sam_inv3x3(P, P1)) return 1;
  if(sam_inv3x3(Q, Q1)) return 1;

  /* v=inv(P)'*P4, r=inv(Q)'*Q4 */
  pq4=P+3*3; /* P4 */
  v[0]=P1[0]*pq4[0] + P1[3]*pq4[1] + P1[6]*pq4[2];
  v[1]=P1[1]*pq4[0] + P1[4]*pq4[1] + P1[7]*pq4[2];
  v[2]=P1[2]*pq4[0] + P1[5]*pq4[1] + P1[8]*pq4[2];
  
  pq4=Q+3*3; /* Q4 */
  r[0]=Q1[0]*pq4[0] + Q1[3]*pq4[1] + Q1[6]*pq4[2];
  r[1]=Q1[1]*pq4[0] + Q1[4]*pq4[1] + Q1[7]*pq4[2];
  r[2]=Q1[2]*pq4[0] + Q1[5]*pq4[1] + Q1[8]*pq4[2];

  w1=v[0]*r[2]/(r[0]*v[2]);
  w2=v[1]*r[2]/(r[1]*v[2]);

  /* compute H as P'*diag([w1, w2, 1])*inv(Q)' */
  /* wQ1t=W*inv(Q)' */
  wQ1t[0]=w1*Q1[0]; wQ1t[1]=w1*Q1[3]; wQ1t[2]=w1*Q1[6];
  wQ1t[3]=w2*Q1[1]; wQ1t[4]=w2*Q1[4]; wQ1t[5]=w2*Q1[7];
  wQ1t[6]=   Q1[2]; wQ1t[7]=   Q1[5]; wQ1t[8]=   Q1[8];

  /* H=P'*wQ1t */
  for(i=0; i<3; ++i){
    int i3=i*3;

    H[i3+0]=P[0*3+i]*wQ1t[0*3+0] + P[1*3+i]*wQ1t[1*3+0] + P[2*3+i]*wQ1t[2*3+0];
    H[i3+1]=P[0*3+i]*wQ1t[0*3+1] + P[1*3+i]*wQ1t[1*3+1] + P[2*3+i]*wQ1t[2*3+1];
    H[i3+2]=P[0*3+i]*wQ1t[0*3+2] + P[1*3+i]*wQ1t[1*3+2] + P[2*3+i]*wQ1t[2*3+2];
  }


  K1[0]=cp->inv_fx;  K1[1]=-cp->s_fxfy;  K1[2]=cp->scy_cxfy_fxfy;
  K1[3]=0.0;         K1[4]=cp->inv_fy;   K1[5]=-cp->cy_fy;
  K1[6]=0.0;         K1[7]=0.0;          K1[8]=1.0;

  /* Hn=K1*H */
  for(i=0; i<3; ++i){
    int i3=i*3;

    Hn[i3+0]=K1[i3+0]*H[0*3+0] + K1[i3+1]*H[1*3+0] + K1[i3+2]*H[2*3+0];
    Hn[i3+1]=K1[i3+0]*H[0*3+1] + K1[i3+1]*H[1*3+1] + K1[i3+2]*H[2*3+1];
    Hn[i3+2]=K1[i3+0]*H[0*3+2] + K1[i3+1]*H[1*3+2] + K1[i3+2]*H[2*3+2];
  }

  /* H=Hn */
  H[0]=Hn[0]; H[1]=Hn[1]; H[2]=Hn[2];
  H[3]=Hn[3]; H[4]=Hn[4]; H[5]=Hn[5];
  H[6]=Hn[6]; H[7]=Hn[7]; H[8]=Hn[8];

  return 0;
}


/* Analytic solution for the P4P problem for coplanar points.
 * m are the image points in pixel coordinates and M their
 * corresponding 3D points.
 * See Appendix B in Fischler and Bolles.
 *
 * returns whether a solution was found, i.e. 0 on failure and
 * 1 on success
 */
int coplanarP4P_FB(struct p3p_calib_params *cp, double m[4][2], double M[4][3], double R[3][3], double t[3])
{
double xM[4][2]; // points on plane Z=0
double T[9], T1[9]; // Z=0 plane to image homography and its inverse
double a1, a2, a3;
double di, th, d0, pan, dcp;
double b1, b2, b3;
double c1, c2, c3;
double xsgn, ysgn, xcp, ycp, zcp;

/* auxiliaries */
double bc, sinth, dcpsinth, maga12, len;
double M_orig[3][3];
double norm, mu, mv, mk;

  /* compute the collineation from world to image plane and its inverse */
  if(W2Icollineation(cp, m, M, xM, T, NULL, NULL)) return 0;
  if(sam_inv3x3(T, T1)) return 0;

  /* (b): VLI=[a1, a2, a3]=T^-t*[0 0 1]' */
  a1=T1[6]; a2=T1[7]; a3=T1[8];
//printf("A: %g %g %g\n", a1, a2, a3);
  maga12=a1*a1 + a2*a2;
  if(maga12<1E-17){ // parallel image and object planes
    fprintf(stderr, "Parallel image & object planes in coplanarP4P_FB, special case not yet implemented\n");
    th=0.0;

    return 0;
//CHECKME: a1==a2==0?
  }

  /* (c) DI */
  di=fabs(a3)/sqrt(maga12);

  /* (d) theta */
  th=atan2(1.0, di); // foc==1.0

  /* (e) VLO=[b1, b2, b3]=T^t*[0 0 1]' */
  b1=T[6]; b2=T[7]; b3=T[8];
//printf("B: %g %g %g\n", b1, b2, b3);

  /* NOTE: following calculation in the paper has a typo and uses T^-t instead of T^-1 */
  /* (f) PPO=[c1, c2, c3]=T^-1*[0 0 1]' */
  c1=T1[2]; c2=T1[5]; c3=T1[8];
//printf("C: %g %g %g\n", c1, c2, c3);

  /* (g) DO */
  bc=b1*c1 + b2*c2 + b3*c3;
  d0=fabs(bc/(c3*sqrt(b1*b1 + b2*b2)));

  /* (h) $ (pan) */
  pan=atan2(-b2, b1);

//printf("di th d0 pan %g %g %g %g\n", di, th, d0, pan);

  /* (i) xsgn, ysgn */
  xsgn=(bc/(b1*c3)<0)? 1 : -1;
  ysgn=(bc/(b2*c3)<0)? 1 : -1;
//printf("xsgn, ysgn %g %g\n", xsgn, ysgn);

  /* (j) DCP, XCP, YCP, ZCP */
  /* CP: center of projection */
  sinth=sin(th);
  dcp=d0*sinth; // distance of CP to the point of the plane pierced by the optical axis
  dcpsinth=dcp*sinth;
  xcp=xsgn*fabs(dcpsinth*cos(pan)) + c1/c3; // CP x
  ycp=ysgn*fabs(dcpsinth*sin(pan)) + c2/c3; // CP y
  zcp=dcp*cos(th); // CP z

//printf("DCP, XCP, YCP, ZCP %g %g %g %g\n", dcp, xcp, ycp, zcp);

  /* Z coordinates of xM[i] are all zero */
  len=sqrt((xM[0][0]-xcp)*(xM[0][0]-xcp) + (xM[0][1]-ycp)*(xM[0][1]-ycp) + zcp*zcp);
  mu=m[0][0]; mv=m[0][1];
  mu = cp->inv_fx * mu - cp->s_fxfy * mv + cp->scy_cxfy_fxfy;
  mv = cp->inv_fy * mv - cp->cy_fy;
  norm=sqrt(mu*mu + mv*mv + 1.0);
  mk=1./norm; mu*=mk; mv*=mk;
  M_orig[0][0]=len*mu;
  M_orig[0][1]=len*mv;
  M_orig[0][2]=len*mk;

  len=sqrt((xM[1][0]-xcp)*(xM[1][0]-xcp) + (xM[1][1]-ycp)*(xM[1][1]-ycp) + zcp*zcp);
  mu=m[1][0]; mv=m[1][1];
  mu = cp->inv_fx * mu - cp->s_fxfy * mv + cp->scy_cxfy_fxfy;
  mv = cp->inv_fy * mv - cp->cy_fy;
  norm=sqrt(mu*mu + mv*mv + 1.0);
  mk=1./norm; mu*=mk; mv*=mk;
  M_orig[1][0]=len*mu;
  M_orig[1][1]=len*mv;
  M_orig[1][2]=len*mk;

  len=sqrt((xM[2][0]-xcp)*(xM[2][0]-xcp) + (xM[2][1]-ycp)*(xM[2][1]-ycp) + zcp*zcp);
  mu=m[2][0]; mv=m[2][1];
  mu = cp->inv_fx * mu - cp->s_fxfy * mv + cp->scy_cxfy_fxfy;
  mv = cp->inv_fy * mv - cp->cy_fy;
  norm=sqrt(mu*mu + mv*mv + 1.0);
  mk=1./norm; mu*=mk; mv*=mk;
  M_orig[2][0]=len*mu;
  M_orig[2][1]=len*mv;
  M_orig[2][2]=len*mk;

  /*
  len=sqrt((xM[3][0]-xcp)*(xM[3][0]-xcp) + (xM[3][1]-ycp)*(xM[3][1]-ycp) + zcp*zcp);
  mu=m[3][0]; mv=m[3][1];
  mu = cp->inv_fx * mu - cp->s_fxfy * mv + cp->scy_cxfy_fxfy;
  mv = cp->inv_fy * mv - cp->cy_fy;
  norm=sqrt(mu*mu + mv*mv + 1.0);
  mk=1./norm; mu*=mk; mv*=mk;
  M_orig[3][0]=len*mu;
  M_orig[3][1]=len*mv;
  M_orig[3][2]=len*mk;
  */

  if(posest_align3Pts(M_orig, M, R, t))
    return 0;

  return 1;
}

#if 0
// coplanar P4P test main()
main()
{
register int i;
struct p3p_calib_params cal;

double R[3][3], t[3];
double H[9];

double M[4][3]={
{0.3790, 0.8374, 66.2935},
{0.4744, -0.0844, 65.4697},
{0.6153, -1.8516, 64.3976},
{0.571487, -1.36497, 64.7534},
};

double m[4][2]={
{345.426935338125, 328.779623340388},
{368.422545223704, 325.275726411889},
{407.271884900435, 323.138646597495},
{395.931588268555, 324.275448815518},
};

double K[9]={
443.59549, 0, 344.89962,
0, 444.751606, 207.652054,
0, 0, 1,
};

  p3p_set_calib(&cal, K);
  i=coplanarP4P_FB(&cal, m, M, R, t);
  printf("P4Pcp returned %d\n", i);

  printf("%g %g %g\n", R[0][0], R[0][1], R[0][2]);
  printf("%g %g %g\n", R[1][0], R[1][1], R[1][2]);
  printf("%g %g %g\n", R[2][0], R[2][1], R[2][2]);

  printf("\n%g %g %g\n", t[0], t[1], t[2]);
}
#endif
/******************************** F & B coplanar P4P end ********************************/


/******************************* Zhang coplanar P4P start *******************************/
#define _CROSSPROD(v, x, y){ (v)[0]=(x)[1]*(y)[2] - (x)[2]*(y)[1]; (v)[1]=(x)[2]*(y)[0] - (x)[0]*(y)[2]; (v)[2]=(x)[0]*(y)[1] - (x)[1]*(y)[0]; }

/* compute the camera pose [R|t] from a known world-image homography H.
 * H is assumed to map scene to *normalized* image points (i.e., calibration
 * K implicitly assumed known, H=K^-1*Hu, Hu the homography for points in pixelxi
 * coordinates).
 *
 * see HZ2 p.234, Zhang "Flexible Camera Calibration By Viewing a Plane From Unknown Orientations". 
 *
 * [R|t]=[r1, r2, r1xr2, t], where
 * r1=h1*lam, r2=h2*lam, t=h3*lam
 * and lam=1/||h1||=1/||h2||
 *
 * Note that this results from H being [r1 r2 t] and the plane is at Z=0
 *
 * returns 0 if successful, 1 otherwise
 */
static int poseFromNormHomog(double H[9], double R[3], double t[3])
{
double r1[3], r2[3], r3[3], tmp[9], lam1, lam2;
register int i, j, k;
/* remaining variables are needed in the SVD */
const int three=3;
double svals[3], U[9], Vt[9], sum;

  /* lam */
  lam1=sqrt(H[0]*H[0] + H[3]*H[3] + H[6]*H[6]);
  lam2=sqrt(H[1]*H[1] + H[4]*H[4] + H[7]*H[7]);
  lam1=2.0/(lam1+lam2); // average
  /* [r1, r2, t] */
  r1[0]=H[0]*lam1; r1[1]=H[3]*lam1; r1[2]=H[6]*lam1;
  r2[0]=H[1]*lam1; r2[1]=H[4]*lam1; r2[2]=H[7]*lam1;
   t[0]=H[2]*lam1;  t[1]=H[5]*lam1;  t[2]=H[8]*lam1;

  /* R=[r1, r2, r1xr2] */
  _CROSSPROD(r3, r1, r2);
  R[0]=r1[0]; R[1]=r2[0]; R[2]=r3[0];
  R[3]=r1[1]; R[4]=r2[1]; R[5]=r3[1];
  R[6]=r1[2]; R[7]=r2[2]; R[8]=r3[2];

  /* Due to errors during H's estimation, R might not be a rotation matrix (i.e. not orthonormal).
   * Therefore, we need to compute the best rotation matrix R' approximating R. Best is in the 
   * sense of the minimum Frobenius norm of the difference R'-R. This can be achieved using SVD:
   * if R=U S V^t then R'=U V^t
   *
   * see "A Flexible New Technique for Camera Calibration", Z. Zhang, MSR-TR-98-71
  */

  /* copy R (row major) to tmp */
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      tmp[i+j*3]=R[i*3+j];

#if 1 // use custom 3x3 SVD
  svd3(U, svals, Vt, tmp);
  /* svd3 actually returns V; for compatibility with LAPACK which returns V^T,
   * the computed Vt is transposed in place to yield the true Vt
   */
  lam1=Vt[1]; Vt[1]=Vt[3]; Vt[3]=lam1;
  lam2=Vt[2]; Vt[2]=Vt[6]; Vt[6]=lam2;
  lam1=Vt[5]; Vt[5]=Vt[7]; Vt[7]=lam1;
#else // use generic SVD
  {
  double work[64];
  const int worksz=64; /* more than needed */
  int info, iwork[24];

/* Singular Value Decomposition (SVD) */
extern int F77_FUNC(dgesvd)(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu,
                       double *vt, int *ldvt, double *work, int *lwork, int *info);
/* lapack 3.0 routine, faster than dgesvd() */
extern int F77_FUNC(dgesdd)(char *jobz, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt,
                       double *work, int *lwork, int *iwork, int *info);

  //F77_FUNC(dgesvd)("A", "A", (int *)&three, (int *)&three, tmp, (int *)&three, svals, U, (int *)&three, Vt, (int *)&three, work, (int *)&worksz, &info);
  F77_FUNC(dgesdd)("A", (int *)&three, (int *)&three, tmp, (int *)&three, svals, U, (int *)&three, Vt, (int *)&three, work, (int *)&worksz, iwork, &info);
	if(info!=0){
		if(info<0){
			fprintf(stderr, "argument %d of dgesdd illegal in poseFromNormHomog()\n", -info);
			exit(1);
		}
		else{
			fprintf(stderr, "dgesdd (dbdsdc) failed to converge; updating process failed in poseFromNormHomog() [%d]\n", info);
			return 1;
		}
	}

  }
#endif

  /* compute the orthonormal matrix 
   *
   * note that H and thus R is defined up to scale and that the "orthonormalization" process
   * influences R's scale. Scale information for the original R is contained in the
   * singular values. Therefore, we compensate for it by scaling t using the 
   * largest singular value
   */
  lam1=1.0/svals[0];
  for(i=0; i<3; i++){
    for(j=0; j<3; j++){
      for(k=0, sum=0.0; k<3; k++)
        sum+=U[i+k*3]*Vt[j*3+k];
      R[i*3+j]=sum;
    }
    t[i]*=lam1;
  }

  /* finally, we ensure that R's determinant is equal to 1 */
  /* sum=m11 m22 m33 + m12 m23 m31 + m13 m21 m32 - m11 m23 m32 - m12 m21 m33 - m13 m22 m31 */
  sum=R[0]*R[4]*R[8] + R[1]*R[5]*R[6] + R[2]*R[3]*R[7] - R[0]*R[5]*R[7] - R[1]*R[3]*R[8] - R[2]*R[4]*R[6];
  /* alternative formula for the determinant: */
  //sum=R[0]*(R[4]*R[8] - R[7]*R[5]) - R[1]*(R[3]*R[8] - R[5]*R[6]) + R[2]*(R[3]*R[7] - R[6]*R[4]);
  lam1=1.0/sum;
    
  if(lam1==1.0) return 0; /* nothing to do */

  /* rescale R & t */
  for(i=0; i<3; i++) // multiplication of 1st row by lam1 multiplies the determinant by lam1
    R[i]*=lam1;

  for(i=0; i<3; i++)
    t[i]*=lam1;

  return 0;
}


/* Solution for the P4P problem for coplanar points.
 * m are the image points in pixel coordinates and M
 * their corresponding 3D points.
 * See Zhang's PAMI paper
 *
 * returns whether a solution was found, i.e. 0 on failure
 * and 1 on success
 */
int coplanarP4P_Zhang(struct p3p_calib_params *cp, double m[4][2], double M[4][3], double R[3][3], double t[3])
{
register int i;
double xM[4][2]; // points on plane Z=0
double H[9]; // Z=0 plane to image homography
double Ra[9], ta[3]; // transformation aligning with Z=0 plane
double Rxy[3*3], txy[3]; // pose wrt to plane at Z=0

  /* compute the collineation from world to image plane */
  if(W2Icollineation(cp, m, M, xM, H, Ra, ta)) return 0;

  if(poseFromNormHomog(H, (double *)Rxy, txy)) return 0;

  /* R=Rxy*Ra */
  /* R=K1*H */
  for(i=0; i<3; ++i){
    int i3=i*3;

    R[i][0]=Rxy[i3+0]*Ra[0*3+0] + Rxy[i3+1]*Ra[1*3+0] + Rxy[i3+2]*Ra[2*3+0];
    R[i][1]=Rxy[i3+0]*Ra[0*3+1] + Rxy[i3+1]*Ra[1*3+1] + Rxy[i3+2]*Ra[2*3+1];
    R[i][2]=Rxy[i3+0]*Ra[0*3+2] + Rxy[i3+1]*Ra[1*3+2] + Rxy[i3+2]*Ra[2*3+2];
  }

  /* t=txy + Rxy*ta */
  t[0]=txy[0] + Rxy[0]*ta[0] + Rxy[1]*ta[1] + Rxy[2]*ta[2];
  t[1]=txy[1] + Rxy[3]*ta[0] + Rxy[4]*ta[1] + Rxy[5]*ta[2];
  t[2]=txy[2] + Rxy[6]*ta[0] + Rxy[7]*ta[1] + Rxy[8]*ta[2];

  return 1;
}
/******************************** Zhang coplanar P4P end ********************************/

