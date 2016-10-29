/*
/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////
*/

#ifndef _POSEST_H
#define _POSEST_H

/* define the following if you want to build a DLL with MSVC */
/**
#define DLL_BUILD 
**/

#ifdef __cplusplus
extern "C" {
#endif

#define POSEST_VERSION    "1.1 (May 2015)"

/* non-linear refinement cost functions */
#define POSEST_REPR_ERR_NO_NLN_REFINE     0 /* reprojection error without non-linear refinement */
#define POSEST_REPR_ERR_NLN_REFINE        1 /* reprojection error with non-linear refinement */
#define POSEST_REPR_ERR_NLN_MLSL_REFINE   2 /* reprojection error with non-linear refinement & multistart scheme */
#define POSEST_OBJSPC_ERR_LHM             3 /* object space error with LHM */

#define NUM_RTPARAMS      6 /* #params involved in rotation + translation (pose) */
#define NUM_RTFPARAMS     7 /* #params involved in rotation + translation+ focal length */
#define NUM_PPARAMS       12 /* #params involved in projection matrix */

/* use as: extern POSEST_API_MOD int POSEST_CALL_CONV func(...) */
#if defined(DLL_BUILD) && defined(_MSC_VER) /* build DLLs with MSVC only! */
#define POSEST_API_MOD    __declspec(dllexport)
#define POSEST_CALL_CONV  __cdecl
#else /* define empty */
#define POSEST_API_MOD 
#define POSEST_CALL_CONV
#endif /* DLL_BUILD && _MSC_VER */

#define POSEST_ERR     -1
#define POSEST_OK       0

/* rotation initializations for LHM */
#define LHM_INITR_WEAKP       0 // weak perspective approximation
#define LHM_INITR_IDENT       1 // I3
#define LHM_INITR_SUPPLIED    2 // supplied R

/* posest.c */
extern int posest(double (*pts2D)[2], double (*pts3D)[3], int nmatches, double inlPcent, double K[9],
                  double *pp, int npp, int NLrefine, int *idxOutliers, int *nbOutliers, int verbose);

extern int posestBinoc(double (*pts2DL)[2], double (*pts3DL)[3], int nmatchesL, double PextL[NUM_PPARAMS],
                       double (*pts2DR)[2], double (*pts3DR)[3], int nmatchesR, double PextR[NUM_PPARAMS],
                       double inlPcent, double inscale, double rtLs[NUM_RTPARAMS+1], int estScale,
                       int NLrefine, int *idxOutliersLR, int *nbOutliersLR, int verbose);

/* RMS & RMedS errors for a P */
extern void posest_RMS_RMedS(double (*inpts2D)[2], double (*inpts3D)[3], int nmatches, double P[NUM_PPARAMS],
                             double *rms, double *rmeds);

/* first, second and third quartiles for a P */
extern void posest_quartiles(double (*inpts2D)[2], double (*inpts3D)[3], int nmatches, double P[NUM_PPARAMS],
                             double *Q1, double *Q2, double *Q3);

/* specified percentile for a P */
extern double posest_percentile(double (*pts2D)[2], double (*pts3D)[3], int nmatches, double P[NUM_PPARAMS], double pct);

extern void posest_PfromKRt(double P[NUM_PPARAMS], double K[9], double rt[NUM_RTPARAMS]);

/* RMS & RMedS errors for the binocular case */
extern void posestBinoc_RMS_RMedS(double (*pts2DL)[2], double (*pts3DL)[3], int nmatchesL, double PextL[NUM_PPARAMS],
                                  double (*pts2DR)[2], double (*pts3DR)[3], int nmatchesR, double PextR[NUM_PPARAMS],
                                  double rtLs[NUM_RTPARAMS+1],
                                  double *rms, double *rmeds);

/* first, second and third quartiles for the binocular case */
extern void posestBinoc_quartiles(double (*pts2DL)[2], double (*pts3DL)[3], int nmatchesL, double PextL[NUM_PPARAMS],
                                  double (*pts2DR)[2], double (*pts3DR)[3], int nmatchesR, double PextR[NUM_PPARAMS],
                                  double rtLs[NUM_RTPARAMS+1],
                                  double *Q1, double *Q2, double *Q3);

/* specified percentile for the binocular case */
extern double posestBinoc_percentile(double (*pts2DL)[2], double (*pts3DL)[3], int nmatchesL, double PextL[NUM_PPARAMS],
                                     double (*pts2DR)[2], double (*pts3DR)[3], int nmatchesR, double PextR[NUM_PPARAMS],
                                     double rtLs[NUM_RTPARAMS+1],
                                     double pct);

/* lhm.c */
extern int posestRT_LHM(double (*pts2D)[2], double (*pts3D)[3], int nmatches, double inlPcent, double K[9],
                        int howtoinit, double p[NUM_RTPARAMS], int *idxOutliers, int *nbOutliers, int verbose);

#ifdef __cplusplus
}
#endif

#endif /* _POSEST_H */
