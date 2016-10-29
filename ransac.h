/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2011  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _RANSAC_H
#define _RANSAC_H

struct rng_state *state;

extern int ransac_numtries(int p, double i, double P);
extern int **ransac_allocsets(int sizeSet, int nbSets);
extern void ransac_freesets(int **sets);
extern void ransac_initrandom(struct rng_state *state);
extern void ransac_setrepeatablerandom(int flag);
extern double ransac_getthresh(double sigma, int dof);

extern int ransacfit(int nbData, int sizeSet, int **sets, int nbSets,
                  void (*residual)(double *x, int nb, void *adata, double *resid),
                  int (*estimator)(double *x, int nb, int *indexes, void *adata),
                  int isResidualSqr, int verbose, int maxNbSol, double consensusThresh,
                  int prematureConsensus, int dim, double percentageOfGoodData, void *adata,
                  double *estimate, int *bestSet, int *outliersMap, int *nbOutliers);

extern int mlesacfit(int nbData, int sizeSet, int **sets, int nbSets,
                  void (*residual)(double *x, int nb, void *adata, double *resid),
                  int (*estimator)(double *x, int nb, int *indexes, void *adata),
                  int isResidualSqr, int verbose, int maxNbSol, double consensusThresh,
                  int prematureConsensus, int dim, double percentageOfGoodData, void *adata,
                  double *estimate, int *bestSet, int *outliersMap, int *nbOutliers);

extern int presacfit(int nbData, int sizeSet, int **sets, int nbSets,
                  void (*residualb)(double *x, int nx, int obsno, void *adata, double *resid),
                  void (*residuald)(double *x, int nb, void *adata, double *resid),
                  int (*estimator)(double *x, int nb, int *indexes, void *adata),
                  int isResidualSqr, int verbose, int maxNbSol, double consensusThresh,
                  int dim, double percentageOfGoodData, int blocksz, void *adata,
                  double *estimate, int *outliersMap, int *nbOutliers);

#endif /* _RANSAC_H */
