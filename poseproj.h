
/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _POSEPROJ_H
#define _POSEPROJ_H

extern void calc_poseProjRT(double K[9], double rt[6], double M[3], double m[2]);
extern void calc_poseProjRTJac(double K[9], double rt[6], double M[3], double m_grad0[6], double m_grad1[6]);

extern void calc_poseProjRTF(double rtf[7], double u0v0[2], double M[3], double m[2]);
extern void calc_poseProjRTFJac(double rtf[7], double u0v0[2], double M[3], double m_grad0[7], double m_grad1[7]);

extern void calc_poseProjRTBinoc(double KR[9], double rtL[6], double stereo_rt[6], double M[3], double m[2]);
extern void calc_poseProjRTBinocJac(double KR[9], double rtL[6], double stereo_rt[6], double M[3], double mR_grad0[6], double mR_grad1[6]);

extern void calc_poseProjRTScaleJac(double K[9], double rt[6], double scl[1], double M[3], double m_grad0[7], double m_grad1[7]);
extern void calc_poseProjRTBinocScaleJac(double KR[9], double rtL[6], double scl[1], double stereo_rt[6], double M[3],
                                         double mR_grad0[7], double mR_grad1[7]);

extern void calc_poseProjPJac(double p[12],double M[3],double m_grad[2][12]);
extern void calc_posePRTJac(double K[9],double rt[6],double p_grad[12][6]);

#endif /* _POSEPROJ_H */
