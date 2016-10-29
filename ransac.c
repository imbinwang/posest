/////////////////////////////////////////////////////////////////////////////////
// 
//  RANdom SAmple Consensus (RANSAC)
//  M.A. Fischler and R.C. Bolles, Communications of the ACM, v.24 n.6, 1981
//  and HZ pp. 101
//
//  Copyright (C) 2002 - 2013  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "compiler.h"
#include "rngs.h"
#include "ransac.h"

#undef GENERATE_RANDOM_SETS_ON_FLY // define this to generate random sets on the fly when none is given
#define SELECT_UNIQUE_RANDOM_SETS // undef this to accept non-unique sets; irrelevant if GENERATE_RANDOM_SETS_ON_FLY defined
#undef ATTEMPT_EARLY_TERMINATION // undef this to disable attempts to end RANSAC/MLESAC early


#define PROBABILITY_CLOSE_TO_ONE 0.991       
/* seems to work in most cases. It is only used if the number of sets is not given by the user. */

/* If the number of combinations is less than MINIMUM_NUMBER_OF_CONFIGURATIONS
 * times the number of sets needed to reach the probability of
 * PROBABILITY_CLOSE_TO_ONE, then generate all possibilities */
#define MINIMUM_NUMBER_OF_CONFIGURATIONS 6

/* If the fraction of inliers at a particular iteration is less than MINIMUM_INLIERS_ADAPT,
 * then no attempt is made to adaptively re-estimate RANSAC's maximum number of iterations
 */
#define MINIMUM_INLIERS_ADAPT 0.07 // 7%

/* Each RANSAC invocation should consider at least this fraction of samples
 * (wrt the supplied or theoretical ones)
 */
#define MINIMUM_ITERATIONS_FRAC 0.10 // 10%
/*************************************************************************************************/
/*
  set the flag specifying if random has to be initialized with a fixed seed. 
  If flag is set to nonzero (default), several calls to the random routine
  will yield the same result, otherwise not.
  This function should be called prior to calling any other ransacxxx() one!
 */
static int _repeatablerandom_=1;

void ransac_setrepeatablerandom(int flag)
{
  _repeatablerandom_=flag;
}

/*
 * initialize the random sequence
 */
void ransac_initrandom(struct rng_state *state)
{
  /* rng_init() will do nothing if seeded already */
  rng_init(_repeatablerandom_? RNG_REPEATABLE : RNG_NONREPEATABLE, state);
}
/*************************************************************************************************/


/* sort routines for ints, avoiding the use of qsort */

static void shellsort_int(int *v, int n)
{
register int gap, i, j, jg;
register int tmp;

  for(gap=n>>1; gap>0; gap>>=1) // gap=n/2, gap/=2
    for(i=gap; i<n; ++i)
      for(j=i-gap; j>=0 && v[j]>v[jg=(j+gap)]; j-=gap){
        tmp=v[j];
        v[j]=v[jg];
        v[jg]=tmp;
      }
}

inline static void insertionsort_int(int *v, int n)
{
register int i, j, tmp, *v1=v+1;

  for(i=0; i<n; ++i){
    for(j=i-1, tmp=v[i]; j>=0; --j){
      if(v[j]<=tmp) break;
      v1[j]=v[j];
    }
    v1[j]=tmp;
  }
}

/* initialize arr with "inside-out" Fisher-Yates shuffle */
static inline void ransac_shuffle(struct rng_state *rstate, int *arr, int n)
{
register int i, j;

  for(i=0; i<n; ++i){
    j=(int)rng_rint(state, i+1); // random int in {0 ... i}
    arr[i]=arr[j];
    arr[j]=i;
  }
}

#if 0
/*
 * generate one subset with indices from 0 to n-1, making sure that all its elements are distinct
 */
static inline void ransac_gensubset(struct rng_state *rstate, int n, int sizeSet, int subset[])
{
register int i, j, r;

  for(i=0; i<sizeSet; ++i){
resample:
    r=subset[i]=(int)rng_rint(rstate, n); // int in [0, n-1]
    for(j=0; j<i; ++j) 
      if(r==subset[j]) goto resample; /* element already in subset, try another one */
  }
}
#endif

/* Select sizeSet items in the range 0 to n-1 (inclusive), ensuring non repeated items
 *
 * work should point to working memory containing all n elements from {0 ... n-1} in some order.
 * The function performs a Fisher-Yates shuffle for sizeSet elements at the beginning of work
 *
 * This is more efficient than the version above
 */
static inline void ransac_gensubset(struct rng_state *rstate, int n, int sizeSet, int subset[], int work[])
{
register int i, r;

  for(i=0; i<sizeSet; ++i){
    r=i + (int)rng_rint(state, n-i); // select index from {i ... n-1} ...
    subset[i]=work[r]; work[r]=work[i]; work[i]=subset[i]; // ...and swap it with work[i]
  }
}

/*
 * generate unique random subsets. For non-unique sets, see ransac_genrandomsets() below
 *
 * returns the number of random sets really generated
 */
static int ransac_genuniquerandomsets(struct rng_state *rstate, int sizeSet, int nbData, int nbSets, int **sets)
{
register int i, j, k;
int tries, *wrk;
void (*howtosort)(int *v, int n);

    ransac_initrandom(rstate);

    if(!(wrk=(int *)malloc(nbData*sizeof(int)))){
	    fprintf(stderr, "Error: Not enough memory in `ransac_genuniquerandomsets'!\n");
	    exit(1);
    }
    ransac_shuffle(rstate, wrk, nbData); // init wrk

    howtosort=(sizeSet<=20)? insertionsort_int : shellsort_int; /* insertion sort for sort lists, shell sort otherwise */

    for(i=0, tries=(3*nbSets)>>1; i<nbSets; ++i){ // tries=3*nbSets/2
resample:
      ransac_gensubset(rstate, nbData, sizeSet, sets[i], wrk);
      (*howtosort)(sets[i], sizeSet);

      /* check if set i is already contained in sets */
      for(j=0; j<i; ++j){
        for(k=0; k<sizeSet; ++k)
          if(sets[i][k]!=sets[j][k]) goto unequalsets; /* sets contain at least one different element */

        /* sets i, j are identical */
        if(--tries==0){ /* premature termination, too many tries failed */
	        /* fprintf(stderr,"Warning: premature termination, too many tries failed in ransac_genuniquerandomsets!\n");*/
          goto terminate;
	      }
        goto resample; /* discard set i */

unequalsets:
        continue; /* dummy statement to prevent MSVC from complaining */
      }
    }

terminate:

    free(wrk);
    return i;
}

/*
 * generate (non-unique) random subsets
 *
 * returns the number of random sets generated
 */
static int ransac_genrandomsets(struct rng_state *rstate, int sizeSet, int nbData, int nbSets, int **sets)
{
register int i;
int *wrk;

    ransac_initrandom(rstate);

    if(!(wrk=(int *)malloc(nbData*sizeof(int)))){
	    fprintf(stderr, "Error: Not enough memory in `ransac_genrandomsets'!\n");
	    exit(1);
    }
    ransac_shuffle(rstate, wrk, nbData); // init wrk

    for(i=0; i<nbSets; ++i){
      ransac_gensubset(rstate, nbData, sizeSet, sets[i], wrk);

      /* don't bother to check if set i is already contained in sets */
    }

    free(wrk);
    return i;
}

/*
 * generate combinations recursively
 */
static void ransac_gencombrec(int sizeSet, int ns[], int is[], int depth, int *nbComb, int **sets)
{
    if(depth==sizeSet-1){
      register int i;

      is[depth]=(depth==0)? 0 : is[depth-1]+1;
      for( ; is[depth]<ns[depth]; ++(is[depth])){
        for(i=0; i<sizeSet; ++i)
		      sets[*nbComb][i]=is[i];
	        ++(*nbComb);
	      }
    }
    else{
      is[depth]=(depth==0)? 0 : is[depth-1]+1;
      for( ; is[depth]<ns[depth]; ++(is[depth]))
	    ransac_gencombrec(sizeSet, ns, is, depth+1, nbComb, sets);
    }
}
    
/*
 * generate all combinations,
 * returns the number of combinations
 */
static int ransac_gencomb(
    int maxInt, /* maximum integer */
    int sizeSet, /* size of a combination */
    int **sets	/* Output: combinations generated */
)
{
int *is, *ns, nbComb;
register int i;

    if(!(is=(int *)malloc(2*sizeSet*sizeof(int)))){
	    fprintf(stderr, "Error: Not enough memory in `ransac_gencomb'!\n");
	    exit(1);
    }
    ns=is+sizeSet;

    for(i=0; i<sizeSet; ++i) 
      ns[i]=maxInt-sizeSet + i+1;

    nbComb=0;
    ransac_gencombrec(sizeSet, ns, is, 0, &nbComb, sets);

    free(is);

    return nbComb;
}

/*
 * number of combinations
 */
static int ransac_ncomb(int n, int m)
{
double nb;
int nbInt;
register int i;

  if(m<0 || m>n){
    fprintf(stderr, "\nError in argument m of `ransac_ncomb'!\n");
    return 0;
  }

  nb=1.0;
  for(i=0; i<m; ++i, --n)
    nb*=n;
  while(m>1)
    nb/=m--;
  nbInt=(nb>INT_MAX)? (INT_MAX) : ((int)nb);

  return nbInt;
}

/*
 * number of tries, i.e. m, according to
 * 1 - (1 - (i)^p)^m > P
 */
int ransac_numtries(
  int p,    /* I: size of a subsample */
  double i, /* I: expected fraction of inliers */
  double P  /* I: probability to have a good subsample */
)
{
double E;
register int j;

    if(i<0.0 || i>1.0 || p<=0 || P>=1.0 || P<=0.0){
      fprintf(stderr, "Error: invalid arguments to `ransac_numtries'!\n");
      exit(1);
    }

    /* lines below compute ((int) ceil(log(1.0 - P) / log(1.0 - pow((double)i, (double)p)))); */

    for(j=1, E=i; j<p; ++j)
      E*=i;
    
    if(E>=.9999999999) return 0; // E close to 1.0, avoid division by -infty
    if(E<=1E-08) return INT_MAX; // E close to 0.0, avoid division by zero
    return (int) ceil(log(1.0 - P)/log(1.0 - E));
}

/*
 * allocate necessary subsets for RANSAC
 */
int **ransac_allocsets(
    int sizeSet, /* I: size of a subset */
		int nbSets	/* I: number of subsets */
)
{
int **sets;
register int i;

    if(!(sets=(int **)malloc(nbSets * sizeof(int *)))){
	    fprintf(stderr, "Error: Not enough memory in `ransac_allocsets'!\n");
	    exit(1);
    }

    if(!(sets[0]=(int *)malloc(sizeSet*nbSets*sizeof(int)))){
      fprintf(stderr, "Error: Not enough memory in `ransac_allocsets'!\n");
      exit(1);
    }

    for(i=1; i<nbSets; ++i)
      sets[i]=sets[i-1] + sizeSet;

    return sets;
}

/*
 * free allocated sets
 */
void ransac_freesets(int **sets)
{
    if(sets==NULL) return;

    free(sets[0]);
    free(sets);
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

#define RANSAC_MAXDOF   16

static double chi2inv[RANSAC_MAXDOF+1]={0.0, // note this is just a placeholder
  6.634896601, 9.210340372, 11.34486673, 13.27670414,
  15.08627247, 16.81189383, 18.47530691, 20.09023503,
  21.66599433, 23.20925116, 24.72497031, 26.21696731, 
  27.68824961, 29.14123774, 30.57791417, 31.99992691

};

/* following values correspond to P=95%:
static double chi2inv[RANSAC_MAXDOF+1]={0.0,
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
double ransac_getthresh(double sigma, int dof)
{
  if(dof>0 && dof<=RANSAC_MAXDOF)
    return sqrt(sigma*sigma*chi2inv[dof]);

  fprintf(stderr, "only DOFs between 1 and %d are accepted by ransac_getthresh()!\n", RANSAC_MAXDOF);
  return 0.0;
}

/* 
 * RANSAC technique
 *
 * returns 1 on success, 0 on failure
 */

int ransacfit(
  int nbData,       /* I: the number of Data, which is ordered from 0 */
  int sizeSet,      /* I: size of each randomly selected subset of data */
  int **sets,       /* I: randomly selected subsets, can be set to NULL forcing a default routine to be used for generating subsets */
  int nbSets,       /* I: the number of subsets, set to 0 if you have no idea */
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
  int prematureConsensus, /* I: if > 0, RANSAC stops when a subset gives a consensus set larger than `prematureConsensus' */
  int dim,          /* I: dimension of the parameter vector */
  double percentageOfGoodData,  /* I: the percentage (between epsilon and 1.0) of good data. often 0.5 */
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
int consGood=0, consNew;
int nbSol;
int freeSets=0, *setbuf=NULL;
int best=-1;
register int i, j, k;
int ret;
struct rng_state rstate={0};
#ifdef ATTEMPT_EARLY_TERMINATION
int minNbSamples;
#endif

    if(nbData<sizeSet){
	    fprintf(stderr,	"Error: the number of data received by RANSAC cannot be less than the size of random subsets!\n");
      return 0;
    }

    /* validate proportion of inliers */
    if(percentageOfGoodData<=0.0 || percentageOfGoodData>=1.0){
	    fprintf(stderr,	"\nBad argument in ransacfit(): percentageOfGoodData must be between 0.0 and 1.0\n"); 
	    exit(1);
    }
    /* validate consensus threshold */
    if(consensusThresh<=0.0){
      fprintf(stderr,	"\nBad argument in ransacfit(): consensusThresh must be > 0.0\n");
      exit(1);
    }

    /* make sure that we won't have a premature stop if provided value is <=0 */
    if(prematureConsensus<=0) prematureConsensus=INT_MAX;

    /* memory allocation */
    if(!(sols=(double *)malloc(maxNbSol*dim*sizeof(double))) || !(resGood=(double *)malloc(2*nbData*sizeof(double)))){
      fprintf(stderr, "Error: Not enough memory in `ransacfit'!\n");
      exit(1);
    }

    resNew=resGood + nbData;

    /* if sets == NULL, use the default routine to generate random subsets (unless GENERATE_RANDOM_SETS_ON_FLY is defined) */
    if(!sets){
      int allComb;

      allComb=ransac_ncomb(nbData, sizeSet);
      if(nbSets<=0)
        nbSets=ransac_numtries(sizeSet, percentageOfGoodData, PROBABILITY_CLOSE_TO_ONE);
      if(allComb<MINIMUM_NUMBER_OF_CONFIGURATIONS*nbSets) { /* generate all combinations */
        sets=ransac_allocsets(sizeSet, allComb);
        freeSets=1;
        nbSets=ransac_gencomb(nbData, sizeSet, sets);
      }
      else{
#ifndef GENERATE_RANDOM_SETS_ON_FLY // randomly generate a few subsets
        sets=ransac_allocsets(sizeSet, nbSets);
        freeSets=1;
# ifdef SELECT_UNIQUE_RANDOM_SETS
        nbSets=ransac_genuniquerandomsets(&rstate, sizeSet, nbData, nbSets, sets);
# else
        nbSets=ransac_genrandomsets(&rstate, sizeSet, nbData, nbSets, sets);
# endif /* SELECT_UNIQUE_RANDOM_SETS */
#else // generate sets on the fly
        if(!(setbuf=(int *)malloc((nbData+sizeSet)*sizeof(int)))){ // working memory for ransac_gensubset() + set storage
          fprintf(stderr, "Error: Not enough memory in ransacfit()!\n");
          exit(1);
        }
        ransac_shuffle(&rstate, setbuf, nbData); // init (part of) setbuf
#endif /* GENERATE_RANDOM_SETS_ON_FLY */
      }
    }

    if(isResidualSqr) consensusThresh*=consensusThresh;

#ifdef ATTEMPT_EARLY_TERMINATION
    minNbSamples=(int)(MINIMUM_ITERATIONS_FRAC*nbSets);
#endif
    for(i=0; i<nbSets; ++i){
	    /* estimate the parameters for the ith subset */
      if(sets)
	      nbSol=(*estimator)(sols, sizeSet, sets[i], adata);
      else{
        ransac_gensubset(&rstate, nbData, sizeSet, setbuf+nbData, setbuf);
        nbSol=(*estimator)(sols, sizeSet, setbuf+nbData, adata);
      }
	    /* now test the solutions on the whole set of data */
	    for(k=0; k<nbSol; ++k){
	      (*residual)(sols + k*dim, nbData, adata, resNew);
	      /* use the residuals to score the consensus set */

        /* MSAC cost function uses the redescending M-estimator suggested by Torr & Zisserman in their MLESAC paper:
         * errNew=sum_j MSAC_RHO(resNew[j], consensusThresh), with MSAC_RHO(e, T) defined as ( ((e)<(T))? (e) : (T) ) 
         */
        if(isResidualSqr){
	        for(j=consNew=0, errNew=0.0; j<nbData; ++j){
            if(resNew[j]<=consensusThresh){
              ++consNew;
              errNew+=resNew[j];
            }
            else
              errNew+=consensusThresh;
          }
        }
        else{ /* use the absolute values */
	        for(j=consNew=0, errNew=0.0; j<nbData; ++j){
		        resNew[j]=(resNew[j]>=0.0)? resNew[j] : -resNew[j]; /*fabs(resNew[j]);*/
            if(resNew[j]<=consensusThresh){
              ++consNew;
              errNew+=resNew[j];
            }
            else
              errNew+=consensusThresh;
          }
        }

        if(verbose) printf("iteration %d of %d, consensus set size %d, data fit error %g\n", i, nbSets, consNew, errNew);
        /* next line is the classical RANSAC scoring test: is this a larger consensus set or a tie with decreased error? */
	      //if(consNew>consGood || (consNew==consGood && errNew<errGood))
	      if(errNew<errGood){ // MSAC test: score inliers according to their fitness to the model
          /* keep the result */
          errGood=errNew;
          consGood=consNew;
          best=i;
          memcpy(resGood, resNew, nbData*sizeof(double));
          memcpy(estimate, sols+k*dim, dim*sizeof(double));

#ifdef ATTEMPT_EARLY_TERMINATION
          if(i<minNbSamples) continue;

          /* adaptively estimate the number of samples */
          {
          int adaptNbSamples;
          double r;
          r=consGood/(double)(nbData);
          adaptNbSamples=
            (r>=MINIMUM_INLIERS_ADAPT)? /* check that enough inliers have been found */
            ransac_numtries(sizeSet, r, PROBABILITY_CLOSE_TO_ONE) :
            nbSets;

          /* check for termination: */
	        if(i+1>=adaptNbSamples ||                  /* number of samples >= than that adaptively estimated */
            consGood>=prematureConsensus ||          /* size of consensus set >= prematureConsensus */
            consGood>=percentageOfGoodData*nbData){  /* size of consensus set >= expected number of inliers */
              i++;
              goto breakout; /* premature termination */
          }
          }
#endif /* ATTEMPT_EARLY_TERMINATION */
        }
	    }
    }

#ifdef ATTEMPT_EARLY_TERMINATION
breakout:
#endif

    if(best==-1){		/* failure */
      fprintf(stderr, "RANSAC failed. Among all trials, no solution was found!\n");
      ret=0;
      goto cleanup;
    }


    if(verbose){
      printf("RANSAC result:\n%d samples drawn, consensus set size = %d\n", i, consGood);
      if(sets){
        printf(" Subset %d : ", best);
        for(i=0; i<sizeSet; ++i)
          printf("%d ", sets[best][i]);
      }
      printf("\n Estimate: ");
      for(i=0; i<dim; ++i)
      printf("%g ", estimate[i]);
      printf("\n");
      fflush(stdout);
    }

    if(bestSet)
      *bestSet=best;

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
    if(freeSets)
	    ransac_freesets(sets);
    if(setbuf) free(setbuf);

    return ret;
}


/* 
 * MLESAC technique
 * see Torr & Zisserman: MLESAC: A New Robust Estimator with Application to Estimating Image Geometry. CVIU, 2000
 *
 * returns 1 on success, 0 on failure
 */

int mlesacfit(
  int nbData,       /* I: the number of Data, which is ordered from 0 */
  int sizeSet,      /* I: size of each randomly selected subset of data */
  int **sets,       /* I: randomly selected subsets, can be set to NULL forcing a default routine to be used for generating subsets */
  int nbSets,       /* I: the number of subsets, set to 0 if you have no idea */
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
  int prematureConsensus, /* I: if > 0, MLESAC stops when a subset gives a consensus set larger than `prematureConsensus' */
  int dim,          /* I: dimension of the parameter vector */
  double percentageOfGoodData,  /* I: the percentage (between epsilon and 1.0) of good data. often 0.5 */
  void *adata,      /* I: pointer to additional data, passed uninterpreted to residual() and estimator() */
  double *estimate, /* O: the corresponding estimate of parameters */
  int *bestSet,     /* O: best subset retained, can be set to NULL */
  int *outliersMap, /* O: contains 1 in indexes of detected outliers, 0 in indexes of inliners */
  int *nbOutliers   /* O: the number of detected outliers */
)
{
double *sols;
double *resGood, *resNew;
int consGood=0, consNew;
int nbSol;
int freeSets=0, *setbuf=NULL;
int best=-1;
register int i, j, k;
int ret;
struct rng_state rstate={0};

double dMinPenalty=DBL_MAX;
double dMix, dResidInlierProb, dResidOutlierProb, meanMix, dCurPenalty;
const double GlobMaxResid=2000.0; // u pg. 142, CHECKME
const double dSigma=1.0;
const double dblSigma_sq=2.0*dSigma*dSigma, dSigma_sqrt2pi=dSigma*sqrt(2.0*M_PI);
register int ii;
#ifdef ATTEMPT_EARLY_TERMINATION
int minNbSamples;
#endif

    if(nbData<sizeSet){
	    fprintf(stderr,	"Error: the number of data received by MLESAC cannot be less than the size of random subsets!\n");
      return 0;
    }

    /* validate proportion of inliers */
    if(percentageOfGoodData<=0.0 || percentageOfGoodData>=1.0){
	    fprintf(stderr,	"\nBad argument in mlesacfit(): percentageOfGoodData must be between 0.0 and 1.0\n"); 
	    exit(1);
    }
    /* validate consensus threshold */
    if(consensusThresh<=0.0){
      fprintf(stderr,	"\nBad argument in mlesacfit(): consensusThresh must be > 0.0\n");
      exit(1);
    }

    /* make sure that we won't have a premature stop if provided value is <=0 */
    if(prematureConsensus<=0) prematureConsensus=INT_MAX;

    /* memory allocation */
    if(!(sols=(double *)malloc(maxNbSol*dim*sizeof(double))) || !(resGood=(double *)malloc(2*nbData*sizeof(double)))){
      fprintf(stderr, "Error: Not enough memory in `mlesacfit'!\n");
      exit(1);
    }

    resNew=resGood + nbData;

    /* if sets == NULL, use the default routine to generate random subsets (unless GENERATE_RANDOM_SETS_ON_FLY is defined) */
    if(!sets){
      int allComb;

      allComb=ransac_ncomb(nbData, sizeSet);
      if(nbSets<=0)
        nbSets=ransac_numtries(sizeSet, percentageOfGoodData, PROBABILITY_CLOSE_TO_ONE);
      if(allComb<MINIMUM_NUMBER_OF_CONFIGURATIONS*nbSets) { /* generate all combinations */
        sets=ransac_allocsets(sizeSet, allComb);
        freeSets=1;
        nbSets=ransac_gencomb(nbData, sizeSet, sets);
      }
      else{
#ifndef GENERATE_RANDOM_SETS_ON_FLY // randomly generate a few subsets
        sets=ransac_allocsets(sizeSet, nbSets);
        freeSets=1;
# ifdef SELECT_UNIQUE_RANDOM_SETS
        nbSets=ransac_genuniquerandomsets(&rstate, sizeSet, nbData, nbSets, sets);
# else
        nbSets=ransac_genrandomsets(&rstate, sizeSet, nbData, nbSets, sets);
# endif /* SELECT_UNIQUE_RANDOM_SETS */
#else // generate sets on the fly
        if(!(setbuf=(int *)malloc((nbData+sizeSet)*sizeof(int)))){ // working memory for ransac_gensubset() + set storage
          fprintf(stderr, "Error: Not enough memory in mlesacfit()!\n");
          exit(1);
        }
        ransac_shuffle(&rstate, setbuf, nbData); // init (part of) setbuf
#endif /* GENERATE_RANDOM_SETS_ON_FLY */
      }
    }

    if(isResidualSqr) consensusThresh*=consensusThresh;

#ifdef ATTEMPT_EARLY_TERMINATION
    minNbSamples=(int)(MINIMUM_ITERATIONS_FRAC*nbSets);
#endif
    for(i=0; i<nbSets; ++i){
	    /* estimate the parameters for the ith subset */
      if(sets)
	      nbSol=(*estimator)(sols, sizeSet, sets[i], adata);
      else{
        ransac_gensubset(&rstate, nbData, sizeSet, setbuf+nbData, setbuf);
        nbSol=(*estimator)(sols, sizeSet, setbuf+nbData, adata);
      }
	    /* now test the solutions on the whole set of data */
	    for(k=0; k<nbSol; ++k){
	      (*residual)(sols + k*dim, nbData, adata, resNew);
	      /* use the residuals to score the consensus set */

        /* find the mixing parameter (gamma) using EM; see eqs.(18) - (21) & MLESAC.m in the GML toolbox */
        if(isResidualSqr)
          for(ii=0, dMix=0.5; ii<5; ++ii){
            for(j=0, meanMix=0.0; j<nbData; ++j){
              dResidInlierProb=dMix*exp(-resNew[j]/dblSigma_sq) / dSigma_sqrt2pi;
              dResidOutlierProb=(1.0 - dMix)/GlobMaxResid;
              meanMix+= dResidInlierProb/(dResidInlierProb + dResidOutlierProb);
            }
            dMix=meanMix/(double)nbData;
          }
        else // as above but with squared resNew[]
          for(ii=0, dMix=0.5; ii<5; ++ii){
            for(j=0, meanMix=0.0; j<nbData; ++j){
              dResidInlierProb=dMix*exp(-(resNew[j]*resNew[j])/dblSigma_sq) / dSigma_sqrt2pi;
              dResidOutlierProb=(1.0 - dMix)/GlobMaxResid;
              meanMix+= dResidInlierProb/(dResidInlierProb + dResidOutlierProb);
            }
            dMix=meanMix/(double)nbData;
          }

        /* find log-likehood of the model; see eq. (10) */
        if(isResidualSqr)
          for(j=consNew=0, dCurPenalty=0.0; j<nbData; ++j){
            dResidInlierProb=dMix*exp(-resNew[j]/dblSigma_sq) / dSigma_sqrt2pi;
            dResidOutlierProb=(1.0 - dMix)/GlobMaxResid;
            dCurPenalty-=log(dResidInlierProb + dResidOutlierProb);

            consNew+=(dResidInlierProb>dResidOutlierProb);
          }
        else // as above but with squared resNew[]
          for(j=consNew=0, dCurPenalty=0.0; j<nbData; ++j){
            dResidInlierProb=dMix*exp(-(resNew[j]*resNew[j])/dblSigma_sq) / dSigma_sqrt2pi;
            dResidOutlierProb=(1.0 - dMix)/GlobMaxResid;
            dCurPenalty-=log(dResidInlierProb + dResidOutlierProb);

            consNew+=(dResidInlierProb>dResidOutlierProb);
          }

        if(verbose) printf("iteration %d of %d, consensus set size %d, mix %g, model log-likehood %g\n", i, nbSets, consNew, dMix, dCurPenalty);
	      if(dCurPenalty<dMinPenalty){
          /* keep the result */
          dMinPenalty=dCurPenalty;
          consGood=consNew;
          best=i;
          memcpy(resGood, resNew, nbData*sizeof(double));
          memcpy(estimate, sols+k*dim, dim*sizeof(double));

#ifdef ATTEMPT_EARLY_TERMINATION
          if(i<minNbSamples) continue;

          /* adaptively estimate the number of samples */
          {
          int adaptNbSamples;
          double r;
          r=consGood/(double)(nbData);
          adaptNbSamples=
            (r>=MINIMUM_INLIERS_ADAPT)? /* check that enough inliers have been found */
            ransac_numtries(sizeSet, r, PROBABILITY_CLOSE_TO_ONE) :
            nbSets;

          /* check for termination: */
	        if(i+1>=adaptNbSamples ||                  /* number of samples >= than that adaptively estimated */
            consGood>=prematureConsensus ||          /* size of consensus set >= prematureConsensus */
            consGood>=percentageOfGoodData*nbData){  /* size of consensus set >= expected number of inliers */
              i++;
              goto breakout; /* premature termination */
          }
          }
#endif /* ATTEMPT_EARLY_TERMINATION */
        }
	    }
    }

#ifdef ATTEMPT_EARLY_TERMINATION
breakout:
#endif

    if(best==-1){		/* failure */
      fprintf(stderr, "MLESAC failed. Among all trials, no solution was found!\n");
      ret=0;
      goto cleanup;
    }


    if(verbose){
      printf("MLESAC result:\n%d samples drawn, consensus set size = %d\n", i, consGood);
      if(sets){
        printf(" Subset %d : ", best);
        for(i=0; i<sizeSet; ++i)
          printf("%d ", sets[best][i]);
      }
      printf("\n Estimate: ");
      for(i=0; i<dim; ++i)
      printf("%g ", estimate[i]);
      printf("\n");
      fflush(stdout);
    }

    if(bestSet)
      *bestSet=best;

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
    if(freeSets)
	    ransac_freesets(sets);
    if(setbuf) free(setbuf);

    return ret;
}


/* 
 * Preemptive RANSAC technique
 * see D. Nistér, Preemptive RANSAC for Live Structure and Motion Estimation. ICCV'03
 *
 * returns 1 on success, 0 on failure
 */

/* k-th smallest out of n doubles */
static double quickSelect(/*struct rng_state *rstate,*/ double a[], int n, int k) 
{
register int i, j;
int l, r, s;
register double pivot, temp;

  //if(k<1 || k>=n) return -1; // k out of range 

  l=0; r=n-1;
  while(1){
    //if(l==r) return a[l];
#if 1 // a[r] as pivot
    pivot=a[r];
#else // random pivot
    i=(int)rng_uniform(rstate, l, r); // rstate should be passed as arg to quickSelect()...
    pivot=a[i]; a[i]=a[r]; a[r]=pivot;
#endif
    for(i=j=l; j<r; ++j)
      if(a[j]<=pivot){
        temp=a[i]; a[i]=a[j]; a[j]=temp; 
        ++i;
      }

    temp=a[r]; a[r]=a[i]; a[i]=temp; 

    s=i-l+1;

    if(k<s) r=i-1;
    else if(k>s) {l=i+1; k-=s;}
    else return a[i];
  }
}

/* block size; lower implies heavier preemption */
#define BLCK_SZ   100

int presacfit(
  int nbData,       /* I: the number of Data, which is ordered from 0 */
  int sizeSet,      /* I: size of each randomly selected subset of data */
  int **sets,       /* I: randomly selected subsets, can be set to NULL forcing a default routine to be used for generating subsets */
  int nbSets,       /* I: the number of subsets, set to 0 if you have no idea */
  void (*residualb)(/* I: function which calculates the residuals for several hypotheses and one observation [breadth first order, differs from ransacfit()!] */
    double *x,      /* I: array of parameter vectors */
    int nx,         /* I: number of parameter vectors */
    int datano,     /* I: observation for which residuals are to be computed */
    void *adata,    /* I: additional data */
    double *resid), /* O: computed residuals */
  void (*residuald)(/* I: function which calculates the residuals a la ransacfit() (depth first order), NULL uses residualb in a less efficient manner */
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
  int dim,          /* I: dimension of the parameter vector */
  double percentageOfGoodData,  /* I: the percentage (between epsilon and 1.0) of good data. often 0.5 */
  int blocksz,      /* I: block size for preemption. a non-positive value indicates the default; lower implies heavier preemption */
  void *adata,      /* I: pointer to additional data, passed uninterpreted to residual() and estimator() */
  double *estimate, /* O: the corresponding estimate of parameters */
  int *outliersMap, /* O: contains 1 in indexes of detected outliers, 0 in indexes of inliners */
  int *nbOutliers   /* O: the number of detected outliers */
)
{
double *sols, *scores;
register double *solsptr, err;
double *resNew;
int *dataperm, nsols, fi, fi_1;
int nbSol;
int freeSets=0, *setbuf=NULL;
register int i, j, k;
int ret;
struct rng_state rstate={0};

    if(nbData<sizeSet){
	    fprintf(stderr,	"Error: the number of data received by RANSAC cannot be less than the size of random subsets!\n");
      return 0;
    }

    /* validate proportion of inliers */
    if(percentageOfGoodData<=0.0 || percentageOfGoodData>=1.0){
	    fprintf(stderr,	"\nBad argument in presacfit(): percentageOfGoodData must be between 0.0 and 1.0\n"); 
	    exit(1);
    }
    /* validate consensus threshold */
    if(consensusThresh<=0.0){
      fprintf(stderr,	"\nBad argument in presacfit(): consensusThresh must be > 0.0\n");
      exit(1);
    }

    /* if sets == NULL, use the default routine to generate random subsets (unless GENERATE_RANDOM_SETS_ON_FLY is defined) */
    if(!sets){
      int allComb;

      allComb=ransac_ncomb(nbData, sizeSet);
      if(nbSets<=0)
        nbSets=ransac_numtries(sizeSet, percentageOfGoodData, PROBABILITY_CLOSE_TO_ONE);
      if(allComb<MINIMUM_NUMBER_OF_CONFIGURATIONS*nbSets) { /* generate all combinations */
        sets=ransac_allocsets(sizeSet, allComb);
        freeSets=1;
        nbSets=ransac_gencomb(nbData, sizeSet, sets);
      }
      else{
#ifndef GENERATE_RANDOM_SETS_ON_FLY // randomly generate a few subsets
        sets=ransac_allocsets(sizeSet, nbSets);
        freeSets=1;
# ifdef SELECT_UNIQUE_RANDOM_SETS
        nbSets=ransac_genuniquerandomsets(&rstate, sizeSet, nbData, nbSets, sets);
# else
        nbSets=ransac_genrandomsets(&rstate, sizeSet, nbData, nbSets, sets);
# endif /* SELECT_UNIQUE_RANDOM_SETS */
#else // generate sets on the fly
        if(!(setbuf=(int *)malloc((nbData+sizeSet)*sizeof(int)))){ // working memory for ransac_gensubset() + set storage
          fprintf(stderr, "Error: Not enough memory for 'setbuf' in presacfit()!\n");
          exit(1);
        }
        ransac_shuffle(&rstate, setbuf, nbData); // init (part of) setbuf
#endif /* GENERATE_RANDOM_SETS_ON_FLY */
      }
    }

    if(isResidualSqr) consensusThresh*=consensusThresh;

    if(!(dataperm=(int *)malloc(nbData*sizeof(int)))){
      fprintf(stderr, "Error: Not enough memory #1 in presacfit()!\n");
      exit(1);
    }

    /* permute observations using the Fisher-Yates (a.k.a. Knuth) shuffle */
    for(i=0; i<nbData; ++i) dataperm[i]=i;
    for(i=nbData; i-->1;  ){  // nbData-1 downto 1
      j=rng_rint(&rstate, i+1); // random int in [0, i]
      k=dataperm[j];
      dataperm[j]=dataperm[i];
      dataperm[i]=k;
    }

    /* allocate memory for hypotheses; at most maxNbSol solutions per set */
    if(!(sols=(double *)malloc(nbSets*maxNbSol*dim*sizeof(double)))){
      fprintf(stderr, "Error: Not enough memory #2 in `presacfit'!\n");
      exit(1);
    }

    /* generate hypotheses by estimating the parameters for each subset.
     * Note that the observations in each subset are not permuted
     */
    if(sets){
      for(i=nsols=0, solsptr=sols; i<nbSets; ++i){
	      nbSol=(*estimator)(solsptr, sizeSet, sets[i], adata);
        nsols+=nbSol;
        solsptr+=nbSol*dim;
      }
    }
    else{
      int *set;

      set=setbuf+nbData;
      for(i=nsols=0, solsptr=sols; i<nbSets; ++i){
        ransac_gensubset(&rstate, nbData, sizeSet, set, setbuf);

        nbSol=(*estimator)(solsptr, sizeSet, set, adata);
        nsols+=nbSol;
        solsptr+=nbSol*dim;
      }
    }

    /* allocate memory for hypotheses scores */
    if(!(scores=(double *)malloc(nsols*sizeof(double)))){
      fprintf(stderr, "Error: Not enough memory #3 in `presacfit'!\n");
      exit(1);
    }
    memset(scores, 0, nsols*sizeof(double)); // init L_{-1}(h) to zero
    k=(nbData>=nsols)? nbData : nsols; // max(nbData, nsols)
    if(!(resNew=(double *)malloc(k*sizeof(double)))){
      fprintf(stderr, "Error: Not enough memory #4 in `presacfit'!\n");
      exit(1);
    }

    if(blocksz<=0) blocksz=BLCK_SZ;
    for(i=0, fi=nsols; i<nbData; ++i){
      double kth;

      fi_1=fi; // f(i-1)
      fi=nsols/(1<<(i/blocksz)); // f(i)

      if(fi<fi_1){ // note that 'if' condition fails for i==0!
        /* reorder hypotheses */
        memcpy(resNew, scores, fi_1*sizeof(double));
        kth=quickSelect(/*&rstate,*/ resNew, fi_1, fi); //  select fi-th smallest score
        for(j=k=0; j<fi_1; ++j)
          if(scores[j]<=kth){
            register int l;
            register double *solsrc;

            scores[k]=scores[j];
            for(l=dim, solsptr=sols+k*dim, solsrc=sols+j*dim; l-->0;  )
              *solsptr++=*solsrc++;
            k++;
          }
        if(verbose) printf("iteration %d of %d, %d hypotheses, cutoff score %g\n", i, nbData, fi, kth);
      }

      if(fi==1) break;

      /* compute L_i(h) from L_{i-1}(h) */ 
      (*residualb)(sols, fi, dataperm[i], adata, resNew);

      /* MSAC cost function uses the redescending M-estimator suggested by Torr & Zisserman in their MLESAC paper:
       * scores[j]=scores[j] + MSAC_RHO(resNew[j], consensusThresh), with MSAC_RHO(e, T) defined as ( ((e)<(T))? (e) : (T) ) 
       */
      if(isResidualSqr){
        for(j=0; j<fi; ++j){
          err=resNew[j];
          scores[j]+=(err<=consensusThresh)? err : consensusThresh;
        }
      }
      else{ /* use the absolute values */
        for(j=0; j<fi; ++j){
          err=(resNew[j]>=0.0)? resNew[j] : -resNew[j]; /*fabs(resNew[j]);*/
          scores[j]+=(err<=consensusThresh)? err : consensusThresh;
        }
      }
      if(verbose) printf("iteration %d of %d, retaining %d from %d hypotheses\n", i, nbData, fi, nsols);
    }

    if(i==0){		/* failure */
      fprintf(stderr, "Preemptive RANSAC failed. Among all trials, no solution was found!\n");
      ret=0;
      goto cleanup;
    }

    /* find best hypothesis */
    for(j=0, err=DBL_MAX; j<fi; ++j)
      if(scores[j]<err){
        err=scores[j];
        k=j;
      }
    memcpy(estimate, sols+k*dim, dim*sizeof(double));

    if(verbose){
      printf("Preemptive RANSAC result:\nretained %d hypotheses, evaluated %d out of %d observations\n", fi, i, nbData);
      printf("Estimate: ");
      for(i=0; i<dim; ++i)
        printf("%g ", estimate[i]);
      printf("\n");
      fflush(stdout);
    }

#if 0
    /* print best remaining hypotheses */
    printf("\n");
    for(i=0; i<fi; ++i){
      printf("sol %d (%g): ", i, scores[i]);
      for(j=0; j<dim; ++j)
        printf("%g ", sols[i*dim+j]);
      printf("\n");
    }
#endif

    /* compute the residuals for the best scoring hypothesis */
    if(residuald)
      (*residuald)(estimate, nbData, adata, resNew);
    else{
      for(i=0; i<nbData; ++i)
        (*residualb)(estimate, 1, i, adata, resNew+i); // resisuals computed inefficiently observation after observation
      }

    /* Now detect the outliers */
    for(i=j=0; i<nbData; ++i){
	    if(resNew[i]>consensusThresh || resNew[i]<-consensusThresh){ // absolute value not necessary for squared residuals...
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
    free(dataperm);
    free(sols);
    free(resNew);
    free(scores);
    if(freeSets)
	    ransac_freesets(sets);
    if(setbuf) free(setbuf);

    return ret;
}
