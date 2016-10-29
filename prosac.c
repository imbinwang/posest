/////////////////////////////////////////////////////////////////////////////////
// 
//  Non-linear calibrated camera pose estimation from 3D - 2D correspondences
//  Copyright (C) 2012  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////


/* 
 * PROgressive SAmple Consensus (PROSAC)
 * see Chum & Matas CVPR 2005
 * 
 * and
 * 
 * RANdom SAmple Consensus (RANSAC)
 * see Fischler & Bolles, CACM v.24 n.6, 1981 and HZ pp. 101
 *
 * Adapted from F. Devernay's prosac.c by Manolis Lourakis, May 2012
 *
 */


/*
 prosac.c
 
 version 1.3 (Sep. 22, 2011)
 
 Author: Frederic Devernay <Frederic.Devernay@inria.fr>

 Description: a sample implementation of the PROSAC sampling algorithm, derived from RANSAC.
 
 Reference:
 O. Chum and J. Matas.
 Matching with PROSAC - progressive sample consensus.
 Proc. of Conference on Computer Vision and Pattern Recognition (CVPR), volume 1, pages 220-226,
 Los Alamitos, California, USA, June 2005.
 ftp://cmp.felk.cvut.cz/pub/cmp/articles/matas/chum-prosac-cvpr05.pdf
 
 Note:
 In the above article, the test is step 2 of Algorithm 1 seem to be reversed, and the test in step 1
 is not consistent with the text. They were corrected in this implementation.

 History:
 version 1.0 (Apr. 14, 2009): initial version
 version 1.1 (Apr. 22, 2009): fix support computation, "The hypotheses are veriﬁed against all data"
 version 1.2 (Mar. 16, 2011): Add comments about beta, psi, and set eta0 to its original value (0.05 rather than 0.01)
 version 1.3 (Sep. 22, 2011): Check that k_n_star is never nore than T_N
 version 1.4 (Sep. 24, 2011): Don't stop until we have found at least the expected number of inliers (improvement over original PROSAC).
 version 1.5 (Oct. 10, 2011): Also stop if t > T_N (maximum number of iterations given the apriori proportion of outliers).
 version 1.6 (Oct. 10, 2011): Rewrite niter_RANSAC() and also use it to update k_n_star.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <float.h>

#include "prosac.h"
#include "compiler.h"

#define PROBABILITY_CLOSE_TO_ONE 0.991 // P_GOOD_SAMPLE, probability that at least one of the random samples picked up by RANSAC is free of outliers


/* NOTE: the PROSAC article sets T_N (the number of iterations before PROSAC becomes RANSAC) to 200000,
 but that means :
 - only 535 correspondences out of 1000 will be used after 2808 iterations (60% outliers)
 -      395                                                 588            (50%)
 -      170                                                 163            (40%)
 (the # of iterations is the # of RANSAC draws for a 0.99 confidence
 of finding the right model given the percentage of outliers)
 
 QUESTION: Is it more reasonable to set it to the maximum number of iterations we plan to
 do given the percentage of outlier?
 
 MY ANSWER: If you know that you won't draw more than XX samples (say 2808, because you only
 accept 60% outliers), then it's more reasonable to set N to that value, in order to give
 all correspondences a chance to be drawn (even if that chance is very small for the last ones).
 Anyway, PROSAC should find many inliers in the first rounds and stop right away.
 
 T_N=2808 gives:
 - only 961 correspondences out of 1000 will be used after 2808 iterations (60% outliers)
 -      595                                                 588            (50%)
 -      177                                                 163            (40%)
 
 */

// beta is the probability that a match is declared inlier by mistake, i.e. the ratio of the "inlier"
// surface by the total surface. The inlier surface is a disc with radius 1.96s for
// homography/displacement computation, or a band with width 1.96*s*2 for epipolar geometry (s is the
// detection noise), and the total surface is the surface of the image.
// YOU MUST ADJUST THIS VALUE, DEPENDING ON YOUR PROBLEM!
#define BETA 0.01

// eta0 is the maximum probability that a solution with more than In_star inliers in Un_star exists and was not found
// after k samples (typically set to 5%, see Sec. 2.2 of [Chum-Matas-05]).
#define ETA0 0.05

/// Computation of the Maximum number of iterations for Ransac
/// with the formula from [HZ] Section: "How many samples" p.119
static inline
int niter_RANSAC(double p, // probability that at least one of the random samples picked up by RANSAC is free of outliers
                 double epsilon, // proportion of outliers
                 int s, // sample size
                 int Nmax) // upper bound on the number of iterations (-1 means INT_MAX)
{
    // compute safely N = ceil(log(1. - p) / log(1. - exp(log(1.-epsilon) * s)))
    double logarg, logval, N;
    if (Nmax == -1) {
        Nmax = INT_MAX;
    }
    assert(Nmax >= 1);
    if (epsilon <= 0.) {
        return 1;
    }
    // logarg = -(1-epsilon)^s
    logarg = -exp(s*log(1.-epsilon)); // C++/boost version: logarg = -std::pow(1.-epsilon, s);
    // logval = log1p(logarg)) = log(1-(1-epsilon)^s)
    logval = log(1.+logarg); // C++/boost version: logval = boost::math::log1p(logarg)
    N = log(1.-p) / logval;
    if (logval  < 0. && N < Nmax) {
        return (int)ceil(N);
    }
    return Nmax;
}

// * Non-randomness: eq. (7) states that i-m (where i is the cardinal of the set of inliers for a wrong
// model) follows the binomial distribution B(n,beta). http://en.wikipedia.org/wiki/Binomial_distribution
// For n big enough, B(n,beta) ~ N(mu,sigma^2) (central limit theorem),
// with mu = n*beta and sigma = sqrt(n*beta*(1-beta)).
// psi, the probability that In_star out of n_star data points are by chance inliers to an arbitrary
// incorrect model, is set to 0.05 (5%, as in the original paper), and you must change the Chi2 value if
// you chose a different value for psi.
static inline
int Imin(int m, int n, double beta) {
    const double mu = n*beta;
    const double sigma = sqrt(n*beta*(1-beta));
    // Imin(n) (equation (8) can then be obtained with the Chi-squared test with P=2*psi=0.10 (Chi2=2.706)
    return (int)ceil(m + mu + sigma*sqrt(2.706));
}

/********************************* RANSAC */

/*
 * number of tries, i.e. m, according to
 * 1 - (1 - (i)^p)^m > P
 */
int prosac_numtries(
  int p,    /* I: size of a subsample */
  double i, /* I: expected fraction of inliers */
  double P  /* I: probability to have a good subsample */
)
{
double E;
register int j;

    if(i<0.0 || i>1.0 || p<=0 || P>=1.0 || P<=0.0){
      fprintf(stderr, "Error: invalid arguments to `prosac_numtries'!\n");
      exit(1);
    }

    /* lines below compute ((int) ceil(log(1.0 - P) / log(1.0 - pow((double)i, (double)p)))); */

    for(j=1, E=i; j<p; ++j)
      E*=i;
    
    if(E>=.9999999999) return 0; // E close to 1.0, avoid division by -infty
    return (int) ceil(log(1.0 - P)/log(1.0 - E));
}


/*
 * Chi-square inverse cumulative distribution function
 * values for P=99% and various DOFs.
 * Values were computed with the following Maple script:
 * with(stats);
 * 
 * for i from 1 to 16 do
 *  statevalf[icdf, chisquare[i]](0.99);
 * od;
 *
 * The matlab fnction chi2inv(p, dof) can also be used.
 */

#define PROSAC_MAXDOF   16

static double chi2inv[PROSAC_MAXDOF+1]={0.0, // note this is just a placeholder
  6.634896601, 9.210340372, 11.34486673, 13.27670414,
  15.08627247, 16.81189383, 18.47530691, 20.09023503,
  21.66599433, 23.20925116, 24.72497031, 26.21696731, 
  27.68824961, 29.14123774, 30.57791417, 31.99992691

};

/* following values correspond to P=95%:
static double chi2inv[PROSAC_MAXDOF+1]={0.0,
  3.841458821, 5.991464547, 7.814727903, 9.487729037,
  11.07049769, 12.59158724, 14.06714045, 15.50731306,
  16.91897760, 18.30703805, 19.67513757, 21.02606982
  22.36203249, 23.68479130, 24.99579014, 26.29622760
};
*/

/* determine the error threshold for a constraint that is
 * to be used for separating outliers from inliers.
 * 
 * sigma is the standard deviation of noise in measurements,
 * dof is the model codimension (i.e., # of squared terms)
 */
double prosac_getthresh(double sigma, int dof)
{
  if(dof>0 && dof<=PROSAC_MAXDOF)
    return sqrt(sigma*sigma*chi2inv[dof]);

  fprintf(stderr, "only DOFs between 1 and %d are accepted by prosac_getthresh()!\n", PROSAC_MAXDOF);
  return 0.0;
}

/********************************* RANSAC */


/* 
 * PROSAC technique
 *
 * returns 1 on success, 0 on failure
 */

int prosacfit(
  int nbData,       /* I: the number of Data, which is assumed ordered in descending order according to a quality function */
  int sizeSet,      /* I: size of each randomly selected subset of data */
  void (*residual)( /* I: function which calculates the residuals */
    double *x,      /* I: parameter vector */
    int nb,         /* I: number of residuals to be computed */
    void *adata,    /* I: additional data */
    double *resid), /* O: computed residuals */
  int (*estimator)( /* I: function which estimate the parameters returning the number of solutions */
    double *x,      /* O: estimated parameter vector */
    int nb,         /* I: number of data to be used */
    int *indexes,   /* I: indexes of data to be used */
    void *adata),   /* I: additional data */
  int isResidualSqr,/* I: set 1 if `residual' computes the squared residuals */
  int verbose,      /* I: verbose mode */
  int maxNbSol,     /* I: maximum number of solutions given by "estimator" */
  double consensusThresh, /* I: threshold for outlier detection. squared internally if 'isResidualSqr' is set */
  int prematureConsensus, /* I: if > 0, PROSAC stops when a subset gives a consensus set larger than `prematureConsensus' */
  int dim,          /* I: dimension of the parameter vector */
  double percentageOfGoodData,  /* I: the percentage (between epsilon and 1.0) of good data should be at least this. often 0.5 */
  void *adata,      /* I: pointer to additional data, passed uninterpreted to residual() and estimator() */
  double *estimate, /* O: the corresponding estimate of parameters */
  int *bestSet,     /* O: best subset retained, can be set to NULL */
  int *outliersMap, /* O: contains 1 in indexes of detected outliers, 0 in indexes of inliners */
  int *nbOutliers   /* O: the number of detected outliers */
)
{
double *sols;
double *resGood, *resNew;
double errGood=DBL_MAX, errNew;
int consNew, nbSol;
int *set, best=-1, premexit=0;
register int i, j, k;
int ret;

const double maxoutlratio=1.0-percentageOfGoodData; // MAX_OUTLIERS_PROPORTION, max allowed outliers proportion in the input data: used to compute T_N (can be as high as 0.95)
// Note on MAX_OUTLIERS_PROPORTION: in this implementation, PROSAC won't stop before having reached the
// corresponding inliers rate on the complete data set.

    const int N = nbData;
    const int m = sizeSet;
    const int T_N = niter_RANSAC(PROBABILITY_CLOSE_TO_ONE, maxoutlratio, sizeSet, -1);
    const double beta = BETA;
    int n_star; // termination length (see sec. 2.2 Stopping criterion)
    int I_n_star; // number of inliers found within the first n_star data points
    int I_N_best; // best number of inliers found so far (store the model that goes with it)
    const int I_N_min = (int)(percentageOfGoodData*N); // the minimum number of total inliers
    int t; // iteration number
    int n; // we draw samples from the set U_n of the top n data points
    double T_n; // average number of samples {M_i}_{i=1}^{T_N} that contain samples from U_n only
    int T_n_prime; // integer version of T_n, see eq. (4)
    int k_n_star; // number of samples to draw to reach the maximality constraint

    if(nbData<sizeSet){
	    fprintf(stderr,	"Error: the number of data received by PROSAC cannot be less than the size of random subsets!\n");
      return 0;
    }

    /* memory allocation */
    if(!(sols=(double *)malloc(maxNbSol*dim*sizeof(double))) || !(resGood=(double *)malloc(2*nbData*sizeof(double)))){
      fprintf(stderr, "Error: Not enough memory in `prosacfit'!\n");
      exit(1);
    }

    resNew=resGood + nbData;

    if(!(set=(int *)malloc(sizeSet*sizeof(int)))){
      fprintf(stderr, "Error: Not enough memory in `prosacfit'!\n");
      exit(1);
    }

    /* validate proportion of inliers */
    if(percentageOfGoodData<=0.0 || percentageOfGoodData>=1.0){
	    fprintf(stderr,	"\nBad argument in prosacfit(): percentageOfGoodData must be between 0.0 and 1.0\n"); 
	    exit(1);
    }
    /* validate consensus threshold */
    if(consensusThresh<=0.0){
      fprintf(stderr,	"\nBad argument in prosacfit(): consensusThresh must be > 0.0\n");
      exit(1);
    }

    /* make sure that we do not have a premature stop if provided value is <=0 */
    if(prematureConsensus<=0) prematureConsensus=INT_MAX;

    if(isResidualSqr) consensusThresh*=consensusThresh;

    if(verbose){
      printf("PROSAC\n");
      printf("number of correspondences (N): %d\n", N);
      printf("sample size (m): %d\n", m);
    }

    n_star = N;
    I_n_star = 0;
    I_N_best = 0;
    t = 0;
    n = m;
    T_n = T_N;
    for(i=0; i<m; i++) {
        T_n *= (double)(n-i)/(N-i);
    }
    T_n_prime = 1;
    k_n_star = T_N;
    // Note: the condition (I_N_best < I_N_min) was not in the original paper, but it is reasonable:
    // we shouldn't stop if we haven't found the expected number of inliers
    while(((I_N_best < I_N_min) || t <= k_n_star) && t < T_N) {
        int I_N; // total number of inliers for that sample
        
        // Choice of the hypothesis generation set
        t = t + 1;
        if(verbose>1) printf("Iteration t=%d, ", t);
        // from the paper, eq. (5) (not Algorithm1):
        // "The growth function is then deﬁned as
        //  g(t) = min {n : T′n ≥ t}"
        // Thus n should be incremented if t > T'n, not if t = T'n as written in the algorithm 1
        if ((t > T_n_prime) && (n < n_star)) {
            double T_nplus1 = (T_n * (n+1)) / (n+1-m);
            n = n+1;
            T_n_prime = T_n_prime + (int)(ceil(T_nplus1 - T_n));
            if(verbose>1)
              printf("g(t)=n=%d, n_star=%d, T_n-1>=%d, T_n>=%d, T_n'=%d...",
                    n, n_star, (int)ceil(T_n), (int)ceil(T_nplus1), T_n_prime);
            T_n = T_nplus1;
        }
        else {
            if(verbose>1)
              printf("g(t)=n=%d, n_star=%d, T_n>=%d, T_n'=%d: ",
                    n, n_star, (int)ceil(T_n), T_n_prime);
        }
        // Draw semi-random sample (note that the test condition from Algorithm1 in the paper is reversed):
        if (t > T_n_prime) {
            // during the finishing stage (n== n_star && t > T_n_prime), draw a standard RANSAC sample
            // The sample contains m points selected from U_n at random
            if(verbose>1) printf("Draw %d points from U_%d... ", m, n);
            prosac_deal(n, m, set); // draw m poins from U_n
        }
        else {
            // The sample contains m-1 points selected from U_{n−1} at random and u_n
            if(verbose>1) printf("Draw %d points from U_%d and point u_%d... ", m-1, n-1, n);
            prosac_deal(n-1, m-1, set); // draw m-1 poins from U_{n-1}
            set[m-1]=n; // add un
        }
        //****** SET set according to draw! ********************
        // INSERT Compute model parameters p_t from the sample M_t
        //printf("Model parameter estimation... ");
	      /* estimate the parameters for this sample */
	      nbSol=(*estimator)(sols, sizeSet, set, adata);

        // INSERT (OPTIONAL): Test for degenerate model configuration (DEGENSAC)
        //                    (i.e. discard the sample if more than 1 model is consistent with the sample)
        // ftp://cmp.felk.cvut.cz/pub/cmp/articles/matas/chum-degen-cvpr05.pdf

        // Find support of the model with parameters p_t
        // From first paragraph of section 2: "The hypotheses are veriﬁed against all data"

	      /* now test the solutions on the whole set of data */
	      for(k=0; k<nbSol; ++k){
	        (*residual)(sols + k*dim, nbData, adata, resNew);
	        /* use the residuals to find size of consensus set */
          if(isResidualSqr){
	          for(j=consNew=0, errNew=0.0; j<nbData; ++j){
              if(resNew[j]<=consensusThresh){
                ++consNew;
                errNew+=resNew[j];
              }
            }
          }
          else{ /* use the absolute values */
	          for(j=consNew=0, errNew=0.0; j<nbData; ++j){
		          resNew[j]=(resNew[j]>=0.0)? resNew[j] : -resNew[j]; /*fabs(resNew[j]);*/
              if(resNew[j]<=consensusThresh){
                ++consNew;
                errNew+=resNew[j];
              }
            }
          }
          I_N = consNew;

          if(verbose) printf("iteration %d, solution %d, consensus set size %d\n", t, k, consNew);

          if(verbose>1) printf("found %d inliers, sol %d!\n", I_N, k);

          /* check if we have a larger consensus set or a tie with decreased error */
          if (I_N > I_N_best || (I_N==I_N_best && errNew<errGood)) {
            int n_best; // best value found so far in terms of inliers ratio
            int I_n_best; // number of inliers for n_best

            // INSERT (OPTIONAL): Do local optimization, and recompute the support (LO-RANSAC)
            // http://cmp.felk.cvut.cz/~matas/papers/chum-dagm03.pdf
            // for the fundamental matrix, the normalized 8-points algorithm performs very well:
            // http://axiom.anu.edu.au/~hartley/Papers/fundamental/ICCV-final/fundamental.pdf
            // ...
            // I_N = findSupport(/* model, sample, */ N, isInlier);
            
            I_N_best = I_N;
            errGood=errNew;
            best=t;
		        memcpy(resGood, resNew, nbData*sizeof(double));
		        memcpy(estimate, sols+k*dim, dim*sizeof(double));
            if(bestSet)
              memcpy(bestSet, set, nbData*sizeof(int));

            // INSERT: Store the best model
		        memcpy(estimate, sols+k*dim, dim*sizeof(double));

            if(I_N_best>=prematureConsensus){ /* size of consensus set >= prematureConsensus */
              premexit=1;
              goto finish;
            }

            // Select new termination length n_star if possible, according to Sec. 2.2.
            // Note: the original paper seems to do it each time a new sample is drawn,
            // but this really makes sense only if the new sample is better than the previous ones.
            n_best = N;
            I_n_best = I_N;
            if (1) { // change to if(0) to disable n_star optimization (i.e. draw the same # of samples as RANSAC)
                int n_test; // test value for the termination length
                int I_n_test; // number of inliers for that test value
                double epsilon_n_best = (double)I_n_best/n_best;

                for(n_test = N, I_n_test = I_N; n_test > m; n_test--) { 
                    // Loop invariants:
                    // - I_n_test is the number of inliers for the n_test first correspondences
                    // - n_best is the value between n_test+1 and N that maximizes the ratio I_n_best/n_best
                    assert(n_test >= I_n_test);

                    // * Non-randomness : In >= Imin(n*) (eq. (9))
                    // * Maximality: the number of samples that were drawn so far must be enough
                    // so that the probability of having missed a set of inliers is below eta=0.01.
                    // This is the classical RANSAC termination criterion (HZ 4.7.1.2, eq. (4.18)),
                    // except that it takes into account only the n first samples (not the total number of samples).
                    // kn_star = log(eta0)/log(1-(In_star/n_star)^m) (eq. (12))
                    // We have to minimize kn_star, e.g. maximize I_n_star/n_star
                    //printf("n_best=%d, I_n_best=%d, n_test=%d, I_n_test=%d\n",
                    //        n_best,    I_n_best,    n_test,    I_n_test);
                    // a straightforward implementation would use the following test:
                    //if (I_n_test > epsilon_n_best*n_test)
                    // However, since In is binomial, and in the case of evenly distributed inliers,
                    // a better test would be to reduce n_star only if there's a significant improvement in
                    // epsilon. Thus we use a Chi-squared test (P=0.10), together with the normal approximation
                    // to the binomial (mu = epsilon_n_star*n_test, sigma=sqrt(n_test*epsilon_n_star*(1-epsilon_n_star)).
                    // There is a significant difference between the two tests (e.g. with the findSupport
                    // functions provided above).
                    // We do the cheap test first, and the expensive test only if the cheap one passes.
                    if (( I_n_test * n_best > I_n_best * n_test ) &&
                        ( I_n_test > epsilon_n_best*n_test + sqrt(n_test*epsilon_n_best*(1.-epsilon_n_best)*2.706) )) {
                        if (I_n_test < Imin(m,n_test,beta)) {
                            // equation 9 not satisfied: no need to test for smaller n_test values anyway
                            break; // jump out of the for(n_test) loop
                        }
                        n_best = n_test;
                        I_n_best = I_n_test;
                        epsilon_n_best = (double)I_n_best/n_best;
                    }

                    // prepare for next loop iteration
                    //I_n_test -= isInlier[n_test-1];
                    I_n_test -= (resNew[n_test-1]<=consensusThresh);
                } // for(n_test ...
            } // n_star optimization

            // is the best one we found even better than n_star?
            if ( I_n_best * n_star > I_n_star * n_best ) {
                assert(n_best >= I_n_best);
                // update all values
                n_star = n_best;
                I_n_star = I_n_best;
                k_n_star = niter_RANSAC(1.-ETA0, 1.-I_n_star/(double)n_star, m, T_N);
                if(verbose>1) printf("new values: n_star=%d, k_n_star=%d, I_n_star=%d, I_min=%d\n", n_star, k_n_star, I_n_star, Imin(m,n_best,beta));
            }
        } // if (I_N > I_N_best)
      } // for(k=0 ...
    } // while(t <= k_n_star ...
finish:
    
    if(best==-1){		/* failure */
      fprintf(stderr, "PROSAC failed. Among all trials, no solution was found!\n");
      ret=0;
      goto cleanup;
    }

    if(verbose>1){
      printf("PROSAC finished, reason:\n");
      if(premexit)
        printf("premature exit t=%d\n", t);
      else if (t > T_N)
          printf("t=%d > T_N=%d (k_n_star=%d)\n", t, T_N, k_n_star);
      else if (t > k_n_star)
          printf("t=%d > k_n_star=%d (T_N=%d)\n", t, k_n_star, T_N);
    }
    
    if(verbose){
      printf("PROSAC result:\n Consensus set size = %d\n", I_N_best) ;
      printf(" Subset %d : ", best) ;
      if(bestSet)
        for(i=0; i<sizeSet; ++i)
          printf("%d ", bestSet[i]);
      printf("\n Estimate: ");
      for(i=0; i<dim; ++i)
      printf("%g ", estimate[i]);
      printf("\n");
      fflush(stdout);
    }

    /* Now detect the outliers */
    for(i=j=0; i<nbData; ++i){
	    if(resGood[i]>consensusThresh){
        outliersMap[i]=1;
	      ++j;
	    }
      else
        outliersMap[i]=0;
    }
    *nbOutliers=j;
    ret=1;

    if(verbose){
      for(i=0; i<nbData; ++i)
        if(outliersMap[i])
          printf(" Data %d is an outlier!\n", i);
    }

cleanup:
    /* free memory */
    free(sols);
    free(resGood);
    free(set);

    return ret;
}
