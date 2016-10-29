/*////////////////////////////////////////////////////////////////////////////////
// 
//  Routines related to motion and structure estimation
//
//  Copyright (C) 2002-13  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//////////////////////////////////////////////////////////////////////////////// */

#ifndef _SAM_H
#define _SAM_H

#ifdef __cplusplus
extern "C" {
#endif

/* sam.c */
extern int sam_inv3x3(double *a, double *a1);
extern double sam_det3x3(double a[9]);
extern void sam_vec2rotmat(double r[3], double R[9]);
extern void sam_axangle2rotmat(double n[3], double phi, double R[9]);
extern void sam_rotmat2vec(double R[9], double r[3]);
extern void sam_quat2rotmat(double q[4], double R[9]);
extern void sam_rotmat2quat(double R[9], double q[4]);
extern void sam_vec2quat(double r[3], double q[4]);
extern void sam_quat2vec(double q[4], double r[3]);
extern void sam_rvecnorm(double r[3]);
extern void sam_composeMotions(double R1[9], double t1[3], double R2[9], double t2[3], double R[9], double t[3]);
extern void sam_composeMotionsr(double r1[3], double t1[3], double r2[3], double t2[3], double r[3], double t[3]);
extern void sam_composeMotionsq(double q1[4], double t1[3], double q2[4], double t2[3], double q[4], double t[3]);
extern void sam_composeMotionsvec(double r1[3], double t1[3], double r2[3], double t2[3], double r[3], double t[3]);
extern void sam_relMotion(double R1[9], double t1[3], double R2[9], double t2[3], double R[9], double t[3]);
extern void sam_inverseMotion(double R[9], double t[3], double Ri[9], double ti[3]);
extern void sam_inverseMotionq(double qrt[7], double qrti[7]);
extern void sam_inverseMotionr(double rt[6], double rti[6]);
extern void sam_xformPoints(double (*pts3D)[3], int npts, double rta[6], double rtb[6], double (*xpts3D)[3]);
extern void sam_PfromKrt(double K[9], double r[3], double t[3], double P[12]);
extern void sam_PfromKRt(double K[9], double R[9], double t[3], double P[12]);
extern void sam_FfromPs(double P0[12], double P1[12], double F[9]);
extern int sam_KRtfromP(double P[12], double K[9], double R[9], double t[3]);
extern void sam_KfromSensor(double foc, double senssz[2], int imgsz[2], double ppt[2], double sk, double K[9]);

/* absor.c */
extern int sam_absorq  (double (*pts0)[3], double (*pts1)[3], int *indx, int npts, double R[9], double t[3], double *sc);
extern int sam_absorsvd(double (*pts0)[3], double (*pts1)[3], int *indx, int npts, double R[9], double t[3], double *sc);
extern int sam_absorqRob(double (*pts0)[3], double (*pts1)[3], int nmatches, int samplesz, double inlPcent, double outlThresh,
                         double R[9], double t[3], double *scale, int *idxOutliers, int *nbOutliers, int verbose);
extern double sam_absordist(double (*pts0)[3], double (*pts1)[3], int *indx, int npts, double R[9], double t[3], double sc);

extern void sam_rtErr(double r[3], double t[3], double rg[3], double tg[3], double *terr, double *aerr);
extern void sam_rtErrRel(double r[3], double t[3], double rg[3], double tg[3], double *terr, double *rerr);

#ifdef __cplusplus
}
#endif

#endif /* _SAM_H */
