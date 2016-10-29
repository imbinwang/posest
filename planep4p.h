/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _PLANEP4P_H
#define _PLANEP4P_H

/* PLANEP4P pose computation.
 * Computes a 3D pose from known coplanar correspondences and camera calibration
 */

#ifdef __cplusplus
extern "C" {
#endif

extern int coplanarP4P_FB(struct p3p_calib_params *cp,
                          double m[4][2], double M[4][3],
                          double R[3][3], double t[3]);

extern int coplanarP4P_Zhang(struct p3p_calib_params *cp,
                             double m[4][2], double M[4][3],
                             double R[3][3], double t[3]);
#ifdef __cplusplus
}
#endif

#endif /* _PLANEP4P_H */
