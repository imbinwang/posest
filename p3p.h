/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _P3P_H
#define _P3P_H

/* P3P pose computation.
 * Computes a 3D pose from known correspondences and camera calibration
 */

#ifdef __cplusplus
extern "C" {
#endif

/* maximum number of solutions to the P3P problem */
#define MAX_NUM_P3P_SOL   4

struct p3p_calib_params{
  double fx, fy, cx, cy, s;
  double inv_fx, inv_fy, cy_fy, s_fxfy, scy_cxfy_fxfy; /* hold elements of K^-1 to avoid recalculations */
};

extern void p3p_set_calib(struct p3p_calib_params *p, double K[9]);

extern int p3p_solve3(struct p3p_calib_params *cp,
                      double m[3][2], double M[3][3],
                      double R[4][3][3], double t[4][3]); 

extern int p3p_solve4_2Derr(
                      struct p3p_calib_params *cp,
                      double m[4][2], double M[4][3],
                      double R[3][3], double t[3]);

extern int p3p_solve4_3Derr(
                      struct p3p_calib_params *cp,
                      double m[4][2], double M[4][3], double norm[3],
                      double R[3][3], double t[3]);

extern int p3p_Kneip(struct p3p_calib_params *cp,
                     double m[3][2], double M[3][3],
                     double R[4][3][3], double t[4][3]);

#ifdef __cplusplus
}
#endif

#endif /* _P3P_H */
