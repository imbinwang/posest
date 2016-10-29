/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _UTIL_H
#define _UTIL_H

/* buckets.c */
extern int posest_genRandomSetsNoBuckets(int sizeSet, int nbData, int nbSets, int **sets);
extern int posest_genRandomSetsWithBuckets(double (*pts)[2], int sizeSet, int nbData, int nbSets, int **sets);

/* align.c */
extern int posest_align3Pts(double M_end[3][3], double XYZ[3][3], double R[3][3], double T[3]);
extern int posest_alignNPts(double (*pts0)[3], double (*pts1)[3], int npts, double R[9], double t[3]);

#endif /* _UTIL_H */
