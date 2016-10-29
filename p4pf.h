/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef P4PF_H
#define P4PF_H

/* P4PF pose computation.
 * Computes a 3D pose from known correspondences and camera calibration
 */

#ifdef __cplusplus
extern "C" {
#endif

/* maximum number of solutions to the P4PF problem */
#define MAX_NUM_P4PF_SOL    8

extern int p4pf_solve(double m[4][2], double M[4][3],
                      double R[MAX_NUM_P4PF_SOL][3][3],
                      double t[MAX_NUM_P4PF_SOL][3],
                      double foc[MAX_NUM_P4PF_SOL]);

#ifdef __cplusplus
}
#endif

#endif /* P4PF_H */
