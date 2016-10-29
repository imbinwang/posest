/////////////////////////////////////////////////////////////////////////////////
// 
//  Least quantile of squares (LMedS for a quantile of 50%)
//  P.J. Rousseeuw, Least Median of Squares Regression, Journal of the American
//  Statistics Association, Vol. 79, No. 388, pp. 871-880, Dec. 1984.
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
#include "lqs.h"

#undef GENERATE_RANDOM_SETS_ON_FLY // define this to generate random sets on the fly when none is given
#define SELECT_UNIQUE_RANDOM_SETS // undef this to accept non-unique sets; irrelevant if GENERATE_RANDOM_SETS_ON_FLY defined


#define ROUND(x) ((int)((x) + 0.5))

#define PROBABILITY_CLOSE_TO_ONE 0.991       
/* seems to work in most cases. Only used when the number of sets is not supplied by the user */

/* If the number of combinations is less than MINIMUM_NUMBER_OF_CONFIGURATIONS
   times the number of sets needed to reach the probability of
   PROBABILITY_CLOSE_TO_ONE, then generate all possibilities */
#define MINIMUM_NUMBER_OF_CONFIGURATIONS 6


/*************************************************************************************************/
/*
  set the flag specifying if random has to be initialized with a fixed seed. 
  If flag is set to nonzero (default), several calls to the random routine
  will yield the same result, otherwise not.
  This function should be called prior to calling any other lqsxxx() one!
 */
static int _repeatablerandom_=1;

void lqs_setrepeatablerandom(int flag)
{
  _repeatablerandom_=flag;
}

/*
 * initialize the random sequence
 */
void lqs_initrandom(struct rng_state *state)
{
  /* rng_init() will do nothing if seeded already */
  rng_init(_repeatablerandom_? RNG_REPEATABLE : RNG_NONREPEATABLE, state);
}
/*************************************************************************************************/

/* this array of precomputed values corresponds to the inverse
   cumulative function for a normal distribution. For more information
   consult the literature (Robust Regression for Outlier Detection,
   rouseeuw-leroy) or the maple routine stats/statevalf[icdf,normald] 
   The values are computed with a step of 2.5% starting from 50% */
static double NICDF[21] = {
        1.4e16,      15.94723940, 7.957896558, 5.287692054, 
			  3.947153876, 3.138344200, 2.595242369, 2.203797543, 
			  1.906939402, 1.672911853, 1.482602218, 1.323775627, 
			  1.188182950, 1.069988721, 0.9648473415, 0.8693011162, 
			  0.7803041458, 0.6946704675, 0.6079568319, 0.5102134568,
			  0.3236002672};

#if 1
 /*
  * The following function has been adapted from
  * http://ndevilla.free.fr/median/median/src/quickselect.c
  *
  *  This Quickselect routine is based on the algorithm described in
  *  "Numerical recipes in C", Second Edition,
  *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
  *  This code by Nicolas Devillard - 1998. Public domain.
  */
 
#define _SWAP(a,b) { _tmp=(a);(a)=(b);(b)=_tmp; }
double lqs_kth_smallest(double *arr, int nelems, int select)
{
register double _tmp;
register int low, high, middle, ll, hh;
 
     low = 0;
     high = nelems - 1;
 
     for (;;) {
         if (high <= low) /* One element only */
             return arr[select];
 
         if (high == low + 1) {  /* Two elements only */
             if (arr[low] > arr[high])
                 _SWAP(arr[low], arr[high]);
             return arr[select];
         }
 
         /* Find median of low, middle and high items; swap into position low */
         middle = (low + high) >> 1; // (low + high) / 2
         if (arr[middle] > arr[high])    _SWAP(arr[middle], arr[high]);
         if (arr[low]    > arr[high])    _SWAP(arr[low],    arr[high]);
         if (arr[middle] > arr[low])     _SWAP(arr[middle], arr[low]);
 
         /* Swap low item (now in position middle) into position (low+1) */
         _SWAP(arr[middle], arr[low+1]);
 
         /* Nibble from each end towards middle, swapping items when stuck */
         ll = low + 1;
         hh = high;
         for (;;) {
             do ++ll; while (arr[low] > arr[ll]);
             do --hh; while (arr[hh]  > arr[low]);
 
             if (hh < ll)
                 break;
 
             _SWAP(arr[ll], arr[hh]);
         }
 
         /* Swap middle item (in position low) back into correct position */
         _SWAP(arr[low], arr[hh]);
 
         /* Re-set active partition */
         if (hh <= select)
             low = ll;
         if (hh >= select)
             high = hh - 1;
     }
}
#undef _SWAP

#else

/*
 * Algorithm from N. Wirth's book, implementation by N. Devillard.
 * This code in public domain.
 * See http://ndevilla.free.fr/median/median/median.html
 */


/*---------------------------------------------------------------------------
   Function :   lqs_kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median. 

                Reference:

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/

double lqs_kth_smallest(double *a, int n, int k)
{
register int i, j, l, m;
register double x, tmp;

    l=0; m=n-1;
    while(l<m){
        x=a[k];
        i=l;
        j=m;
        do{
            while(a[i]<x) ++i;
            while(x<a[j]) --j;
            if(i<=j){
                /* swap a[i], a[j] */
                tmp=a[i];
                a[i]=a[j];
                a[j]=tmp;
                ++i; --j;
            }
        } while(i<=j);
        if(j<k) l=i;
        if(k<i) m=j;
    }
    return a[k];
}

#endif /* 0 */

/******
#define median(a,n) lqs_kth_smallest(a, n, (((n)&1)? ((n)/2) :  (((n)/2)-1)) )
******/

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
static inline void lqs_shuffle(struct rng_state *rstate, int *arr, int n)
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
static inline void lqs_gensubset(struct rng_state *rstate, int n, int sizeSet, int subset[])
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
static inline void lqs_gensubset(struct rng_state *rstate, int n, int sizeSet, int subset[], int work[])
{
register int i, r;

  for(i=0; i<sizeSet; ++i){
    r=i + (int)rng_rint(state, n-i); // select index from {i ... n-1} ...
    subset[i]=work[r]; work[r]=work[i]; work[i]=subset[i]; // ...and swap it with work[i]
  }
}

/*
 * generate unique random subsets. For non-unique sets, see lqs_genrandomsets() below
 *
 * returns the number of random sets really generated
 */
static int lqs_genuniquerandomsets(struct rng_state *rstate, int sizeSet, int nbData, int nbSets, int **sets)
{
register int i, j, k;
int tries, *wrk;
void (*howtosort)(int *v, int n);

    lqs_initrandom(rstate);

    if(!(wrk=(int *)malloc(nbData*sizeof(int)))){
	    fprintf(stderr, "Error: Not enough memory in `lqs_genuniquerandomsets'!\n");
	    exit(1);
    }
    lqs_shuffle(rstate, wrk, nbData); // init wrk

    howtosort=(sizeSet<=20)? insertionsort_int : shellsort_int; /* insertion sort for sort lists, shell sort otherwise */

    for(i=0, tries=(3*nbSets)>>1; i<nbSets; ++i){ // tries=3*nbSets/2
resample:
      lqs_gensubset(rstate, nbData, sizeSet, sets[i], wrk);
      (*howtosort)(sets[i], sizeSet);

      /* check if set i is already contained in sets */
      for(j=0; j<i; ++j){
        for(k=0; k<sizeSet; ++k)
          if(sets[i][k]!=sets[j][k]) goto unequalsets; /* sets contain at least one different element */

        /* sets i, j are identical */
        if(--tries==0){ /* premature termination, too many tries failed */
	        /* fprintf(stderr,"Warning: premature termination, too many tries failed in lqs_genuniquerandomsets!\n");*/
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
static int lqs_genrandomsets(struct rng_state *rstate, int sizeSet, int nbData, int nbSets, int **sets)
{
register int i;
int *wrk;

    lqs_initrandom(rstate);

    if(!(wrk=(int *)malloc(nbData*sizeof(int)))){
	    fprintf(stderr, "Error: Not enough memory in `lqs_genrandomsets'!\n");
	    exit(1);
    }
    lqs_shuffle(rstate, wrk, nbData); // init wrk

    for(i=0; i<nbSets; ++i){
      lqs_gensubset(rstate, nbData, sizeSet, sets[i], wrk);

      /* don't bother to check if set i is already contained in sets */
    }

    free(wrk);
    return i;
}

/*
 * generate combinations recursively
 */
static void lqs_gencombrec(int sizeSet, int ns[], int is[], int depth, int *nbComb, int **sets)
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
	    lqs_gencombrec(sizeSet, ns, is, depth+1, nbComb, sets);
    }
}
    
/*
 * generate all combinations,
 * returns the number of combinations
 */
static int lqs_gencomb(
    int maxInt, /* maximum integer */
    int sizeSet, /* size of a combination */
    int **sets	/* Output: combinations generated */
)
{
int *is, *ns, nbComb;
register int i;

    if(!(is=(int *)malloc(2*sizeSet*sizeof(int)))){
	    fprintf(stderr, "Error: Not enough memory in `lqs_gencomb'!\n");
	    exit(1);
    }
    ns=is+sizeSet;

    for(i=0; i<sizeSet; ++i) 
      ns[i]=maxInt-sizeSet + i+1;

    nbComb=0;
    lqs_gencombrec(sizeSet, ns, is, 0, &nbComb, sets);

    free(is);

    return nbComb;
}

/*
 * number of combinations
 */
static int lqs_ncomb(int n, int m)
{
double nb;
int nbInt;
register int i;

  if(m<0 || m>n){
    fprintf(stderr, "\nError in argument m of `lqs_ncomb'!\n");
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
 * 1 - (1 - (1-e)^p)^m > P
 */
int lqs_numtries(
  int p,    /* I: size of a subsample */
  double i, /* I: expected fraction of inliers */
  double P  /* I: probability to have a good subsample */
)
{
double E;
register int j;

    if(i<0.0 || i>1.0 || p<=0 || P>=1.0 || P<=0.0){
      fprintf (stderr, "Error: invalid arguments to `lqs_numtries'!\n");
      exit(1);
    }

    for(j=1, E=i; j<p; ++j)
      E*=i;
    
    if(E>=.9999999999) return 0; // E close to 1.0, avoid division by -infty
    if(E<=1E-08) return INT_MAX; // E close to 0.0, avoid division by zero
    return (int) ceil(log(1.0 - P)/log(1.0 - E));
}

/*
 * allocate necessary subsets for LQS
 */
int **lqs_allocsets(
    int sizeSet, /* I: size of a subset */
		int nbSets	/* I: number of subsets */
)
{
int **sets;
register int i;

    if(!(sets=(int **)malloc(nbSets * sizeof(int *)))){
	    fprintf(stderr, "Error: Not enough memory in `lqs_allocsets'!\n");
	    exit(1);
    }

    if(!(sets[0]=(int *)malloc(sizeSet*nbSets*sizeof(int)))){
      fprintf(stderr, "Error: Not enough memory in `lqs_allocsets'!\n");
      exit(1);
    }

    for(i=1; i<nbSets; ++i)
      sets[i]=sets[i-1] + sizeSet;

    return sets;
}

/*
 * free allocated sets
 */
void lqs_freesets(int **sets)
{
    if(sets==NULL) return;

    free(sets[0]);
    free(sets);
}


/* 
 * Least quantile squares fitting technique
 *
 * returns 1 on success, 0 on failure
 */
int lqsfit(
  int nbData,       /* I: the number of Data, which is ordered from 0 */
  int sizeSet,      /* I: size of each randomly selected subset of data */
  int **sets,       /* I: randomly selected subsets, can be set to NULL forcing a default routine to be used for generating subsets */
  int nbSets,       /* I: the number of subsets, set to 0 if you have no idea */
  void (*residual)( /* I: function which calculates the residuals */
    double *x,        /* I: parameter vector */
    int nb,           /* I: number of residuals to be computed */
    void *adata,      /* I: additional data */
    double *resid),   /* O: computed residuals */
  int (*estimator)( /* I: function which estimate the parameters returning the number of solutions */
    double *x,        /* O: estimated parameter vector */
    int nb,           /* I: number of data to be used */
    int *indexes,     /* I: indexes of data to be used */
    void *adata),     /* I: additional data */
  int isResidualSqr,/* I: set 1 if `residual' computes the squared residuals */
  int verbose,      /* I: verbose mode */
  int maxNbSol,     /* I: maximum number of solutions given by "estimator" */
  double gate,      /* I: gate value for outlier detection, a value in [2, 4], often 2.5 */
  double prematureResidual, /* I: if > 0.0, LQS stops when a subset gives a quantile residual less than `prematureResidual'.
                                   * The specified value must be in accordance with `residual' */
  int dim,          /* I: dimension of the parameter vector */
  double percentageOfGoodData,  /* I: the percentage (between epsilon and 1.0) of good data. often 0.5 */
  void *adata,      /* I: pointer to additional data, passed uninterpreted to residual() and estimator() */
  double *estimate, /* O: the corresponding estimate of parameters */
  int *bestSet,     /* O: best subset retained, can be set to NULL */
  int *outliersMap, /* O: contains 1 in indexes of detected outliers, 0 in indexes of inliers */
  int *nbOutliers,  /* O: the number of detected outliers */
  double *outlierThresh /* O: residual threshold for determining outliers */
)
{
double *sols;
double *resGood, *resNew, *resTmp;
double quantMin=DBL_MAX, quantNew, sigma;
int quantPos, nbSol;
int freeSets=0, *setbuf=NULL;
int best=-1;
register int i, j, k;
int ret;
struct rng_state rstate={0};

    if(nbData<sizeSet){
      fprintf(stderr, "Error: the number of data received by LQS cannot be less than the size of random subsets!\n");
      return 0;
    }

    /* memory allocation */
    if(!(sols=(double *)malloc(maxNbSol*dim*sizeof(double))) || !(resGood=(double *)malloc(3*nbData*sizeof(double)))){
      fprintf(stderr, "Error: Not enough memory in `lqsfit'!\n");
      exit(1);
    }

    resNew=resGood + nbData;
    resTmp=resNew + nbData;

    /* position of the quantile. If set to zero (meaningless) set it to 0.5 */
    if(percentageOfGoodData<=0.0 || percentageOfGoodData>=1.0){
	    fprintf(stderr,	"\nBad argument in lqsfit: percentageOfGoodData must be between 0.0 and 1.0\n");
      percentageOfGoodData=0.5;
    }
    quantPos=(int)(nbData*percentageOfGoodData); 
    if(quantPos<sizeSet - 1) quantPos=sizeSet - 1;

    /* if sets == NULL, use the default routine to generate random subsets (unless GENERATE_RANDOM_SETS_ON_FLY is defined) */
    if(!sets){
	    int allComb;

	    allComb=lqs_ncomb(nbData, sizeSet);
	    if(nbSets<=0)
	      nbSets=lqs_numtries(sizeSet, percentageOfGoodData, PROBABILITY_CLOSE_TO_ONE);
	    if(allComb<MINIMUM_NUMBER_OF_CONFIGURATIONS*nbSets){ /* generate all combinations */
	      sets=lqs_allocsets(sizeSet, allComb);
	      freeSets=1;
	      nbSets=lqs_gencomb(nbData, sizeSet, sets);
	    }
	    else{
#ifndef GENERATE_RANDOM_SETS_ON_FLY // randomly generate a few subsets
	      sets=lqs_allocsets(sizeSet, nbSets);
	      freeSets=1;
# ifdef SELECT_UNIQUE_RANDOM_SETS
        nbSets=lqs_genuniquerandomsets(&rstate, sizeSet, nbData, nbSets, sets);
# else
        nbSets=lqs_genrandomsets(&rstate, sizeSet, nbData, nbSets, sets);
# endif /* SELECT_UNIQUE_RANDOM_SETS */
#else // generate sets on the fly
        if(!(setbuf=(int *)malloc((nbData+sizeSet)*sizeof(int)))){ // working memory for lqs_gensubset() + set storage
          fprintf(stderr, "Error: Not enough memory in lqsfit()!\n");
          exit(1);
        }
        lqs_shuffle(&rstate, setbuf, nbData); // init (part of) setbuf
#endif /* GENERATE_RANDOM_SETS_ON_FLY */
	    }
    }

    for(i=0; i<nbSets; ++i){
	    /* estimate the parameters for the ith subset */
      if(sets)
	      nbSol=(*estimator)(sols, sizeSet, sets[i], adata);
      else{
        lqs_gensubset(&rstate, nbData, sizeSet, setbuf+nbData, setbuf);
        nbSol=(*estimator)(sols, sizeSet, setbuf+nbData, adata);
      }
	    /* now test the solutions on the whole set of data */
	    for(k=0; k<nbSol; ++k){
	      (*residual)(sols + k*dim, nbData, adata, resNew);
        if(isResidualSqr){
          memcpy(resTmp, resNew, nbData*sizeof(double));
        }
        else{
	        /* we use the absolute residuals */
	        for(j=0; j<nbData; ++j)
		        resTmp[j]=resNew[j]=(resNew[j]>=0.0)? resNew[j] : -resNew[j]; /*fabs(resNew[j]);*/
        }

        /* find the desired quantile of the set of residuals */
        quantNew=lqs_kth_smallest(resTmp, nbData, quantPos);
        /* save the result */
        if(quantNew<quantMin){
          quantMin=quantNew;
          best=i;
          memcpy(resGood, resNew, nbData*sizeof(double));
          memcpy(estimate, (void *)(sols+k*dim), dim*sizeof(double));

          if(quantMin<=prematureResidual)
            i=nbSets;	/* premature termination */
        }
	    }
    }

    if(best==-1){	/* failure, usually does not happen */
      fprintf(stderr, "LQS failed. Among all trials, no solution was found!\n");
      ret=0;
      goto cleanup;
    }

    if(verbose){
      printf("LQS result:\n Min quantile residual = %g\n", quantMin);
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

    if(gate<2.0 || gate>4.0){
	    fprintf(stderr,	"\nBad argument in lqsfit: gate must be between 2.0 and 4.0\n");
      gate=2.5;
    }

    /* compute the robust standard deviation estimate  */
    sigma=NICDF[(int)(ROUND(percentageOfGoodData*20.0))]*(1.0 + 5.0 / (double)(nbData - sizeSet));
    if(isResidualSqr)
	    sigma*=sigma;
    sigma*=gate*quantMin;
    if(isResidualSqr)
	    sigma*=gate;

    *outlierThresh=(isResidualSqr)? sqrt(sigma) : sigma;

    /* Now detect the outliers (use sqrt(KSQR) sigma) */
    for(i=j=0; i<nbData; ++i){
	    if(resGood[i]>sigma){ /* an outlier */
        outliersMap[i]=1;
        ++j;
      }
      else
        outliersMap[i]=0;
    }
    *nbOutliers=j;
    ret=1; /* success */

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
	    lqs_freesets(sets);
    if(setbuf) free(setbuf);

    return ret;
}
