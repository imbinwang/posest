/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2012  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _PROSAC_H
#define _PROSAC_H

extern int prosac_numtries(int p, double i, double P);
extern double prosac_getthresh(double sigma, int dof);

extern int prosacfit(int nbData, int sizeSet,
                  void (*residual)(double *x, int nb, void *adata, double *resid),
                  int (*estimator)(double *x, int nb, int *indexes, void *adata),
                  int isResidualSqr, int verbose, int maxNbSol, double consensusThresh,
                  int prematureConsensus, int dim, double percentageOfGoodData, void *adata,
                  double *estimate, int *bestSet, int *outliersMap, int *nbOutliers);

extern void prosac_deal(int n, int k, int *a);

#endif /* _PROSAC_H */
