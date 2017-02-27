/* Multilevel single-linkage global optimization, adapted from nlopt-2.3, M. Lourakis Dec. 2012 */


/* Copyright (c) 2007-2012 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

/* Multi-Level Single-Linkage (MLSL) algorithm for
   global optimization by random local optimizations (a multistart algorithm
   with "clustering" to avoid repeated detection of the same local minimum), 
   modified to optionally use a Sobol' low-discrepancy sequence (LDS) instead 
   of pseudorandom numbers.  See:

   A. H. G. Rinnooy Kan and G. T. Timmer, "Stochastic global optimization
   methods," Mathematical Programming, vol. 39, p. 27-78 (1987).
       [ actually 2 papers -- part I: clustering methods (p. 27), then 
                              part II: multilevel methods (p. 57) ]

   and also:

   Sergei Kucherenko and Yury Sytsko, "Application of deterministic
   low-discrepancy sequences in global optimization," Computational
   Optimization and Applications, vol. 30, p. 297-318 (2005).

   Compared to the above papers, I made a couple other modifications
   to the algorithm besides the use of a LDS.

   1) first, the original algorithm suggests discarding points
      within a *fixed* distance of the boundaries or of previous
      local minima; I changed this to a distance that decreases with,
      iteration, proportional to the same threshold radius used
      for clustering.  (In the case of the boundaries, I make
      the proportionality constant rather small as well.)

   2) Kan and Timmer suggest using fancy "spiral-search" algorithms
      to quickly locate nearest-neighbor points, reducing the
      overall time for N sample points from O(N^2) to O(N log N)
      However, recent papers seem to show that such schemes (kd trees,
      etcetera) become inefficient for high-dimensional spaces (d > 20),
      and finding a better approach seems to be an open question.  Therefore,
      since I am mainly interested in MLSL for high-dimensional problems
      (in low dimensions we have DIRECT etc.), I punted on this question
      for now.  In practice, O(N^2) (which does not count the #points
      evaluated in local searches) does not seem too bad if the objective
      function is expensive.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "../levmar/levmar.h"

#include "mt19937ar.h"
#include "sobol.h"
#include "redblack.h"
#include "mlsl.h"

/* define the following macro to use a fixed random seed */
//#define MLSL_REPEATABLE_RANDOM

/* data structure for each random/quasirandom point in the population */
typedef struct
{
  double f;			/* function value at x */
  int minimized;		/* if we have already minimized starting from x */
  double closest_pt_d;		/* distance^2 to closest pt with smaller f */
  double closest_lm_d;		/* distance^2 to closest lm with smaller f */
  double p[1];			/* array of length n (K&R struct hack) */
}
pt;

/* all of the data used by the various mlsl routines...it's
   not clear in hindsight that we need to put all of this in a data
   structure since most of the work occurs in a single routine,
   but it doesn't hurt us */
typedef struct
{
  int m; /* # dimensions */
  int n; /* # measurements */
  const double *lb, *ub;
  /* stopping criteria */
  double minf_max;
  int nevals, maxeval;
  int nsamples, minsamples;

  rb_tree pts;			/* tree of points (k == pt), sorted by f */
  rb_tree lms;			/* tree of local minimizers, sorted by function value
				   (k = array of length d+1, [0] = f, [1..d] = x) */

  nlopt_sobol s;		/* sobol data for LDS point generation, or NULL
				   to use pseudo-random numbers */

  double R_prefactor, dlm, dbound, gamma;	/* parameters of MLSL */
  int N;			/* number of pts to add per iteration */
}
mlsl_data;

/* comparison routines to sort the red-black trees by function value */
static int pt_compare (rb_key p1_, rb_key p2_)
{
  pt *p1 = (pt *) p1_;
  pt *p2 = (pt *) p2_;
  if (p1->f < p2->f)
    return -1;
  if (p1->f > p2->f)
    return +1;
  return 0;
}

static int lm_compare (double *k1, double *k2)
{
  if (*k1 < *k2)
    return -1;
  if (*k1 > *k2)
    return +1;
  return 0;
}

/* Euclidean distance |x1 - x2|^2 */
static double distance2 (int n, const double *x1, const double *x2)
{
  int i;
  double d = 0.;

  for (i = 0; i < n; ++i){
      double dx = x1[i] - x2[i];
      d += dx * dx;
  }
  return d;
}

/* find the closest pt to p with a smaller function value;
   this function is called when p is first added to our tree */
static void find_closest_pt (int n, rb_tree * pts, pt * ppt)
{
  rb_node *node = rb_tree_find_lt (pts, (rb_key) ppt);
  double closest_d = HUGE_VAL;

  while (node) {
      double d = distance2 (n, ppt->p, ((pt *) node->k)->p);
      if (d < closest_d)
	      closest_d = d;
      node = rb_tree_pred (node);
  }
  ppt->closest_pt_d = closest_d;
}

/* find the closest local minimizer (lm) to p with a smaller function value;
   this function is called when p is first added to our tree */
static void find_closest_lm (int n, rb_tree * lms, pt * ppt)
{
  rb_node *node = rb_tree_find_lt (lms, &ppt->f);
  double closest_d = HUGE_VAL;

  while (node){
      double d = distance2 (n, ppt->p, node->k + 1);
      if (d < closest_d)
	      closest_d = d;
      node = rb_tree_pred (node);
  }
  ppt->closest_lm_d = closest_d;
}

/* given newpt, which presumably has just been added to our
   tree, update all pts with a greater function value in case
   newpt is closer to them than their previous closest_pt ...
   we can ignore already-minimized points since we never do
   local minimization from the same point twice */
static void pts_update_newpt (int n, rb_tree * pts, pt * newpt)
{
  rb_node *node = rb_tree_find_gt (pts, (rb_key) newpt);

  while (node){
      pt *ppt = (pt *) node->k;
      if (!ppt->minimized){
        double d = distance2 (n, newpt->p, ppt->p);
        if (d < ppt->closest_pt_d)
          ppt->closest_pt_d = d;
      }
      node = rb_tree_succ (node);
  }
}

/* given newlm, which presumably has just been added to our
   tree, update all pts with a greater function value in case
   newlm is closer to them than their previous closest_lm ...
   we can ignore already-minimized points since we never do
   local minimization from the same point twice */
static void pts_update_newlm (int n, rb_tree * pts, double *newlm)
{
  pt tmp_pt;
  rb_node *node;

  tmp_pt.f = newlm[0];
  node = rb_tree_find_gt (pts, (rb_key) & tmp_pt);
  while (node) {
      pt *ppt = (pt *) node->k;
      if (!ppt->minimized){
        double d = distance2 (n, newlm + 1, ppt->p);
        if (d < ppt->closest_lm_d)
          ppt->closest_lm_d = d;
      }
      node = rb_tree_succ (node);
  }
}

static int is_potential_minimizer (mlsl_data * mlsl, pt * ppt,
			double dpt_min, double dlm_min, double dbound_min)
{
  int i, m = mlsl->m;
  const double *lb = mlsl->lb;
  const double *ub = mlsl->ub;
  const double *p = ppt->p;

  if (ppt->minimized)
    return 0;

  if (ppt->closest_pt_d <= dpt_min * dpt_min)
    return 0;

  if (ppt->closest_lm_d <= dlm_min * dlm_min)
    return 0;

  for (i = 0; i < m; ++i)
    if ((p[i] - lb[i] <= dbound_min || ub[i] - p[i] <= dbound_min)
	&& ub[i] - lb[i] > dbound_min)
      return 0;

  return 1;
}

#define K2PI (6.2831853071795864769252867665590057683943388)

/* compute Gamma(1 + n/2)^{1/n} ... this is a constant factor
   used in MLSL as part of the minimum-distance cutoff*/
static double gam (int n)
{
  /* use Stirling's approximation:
     Gamma(1 + z) ~ sqrt(2*pi*z) * z^z / exp(z)
     ... this is more than accurate enough for our purposes
     (error in gam < 15% for d=1, < 4% for d=2, < 2% for d=3, ...)
     and avoids overflow problems because we can combine it with
     the ^{1/n} exponent */
  double z = n / 2;
  return sqrt (pow (K2PI * z, 1.0 / n) * z) * exp (-0.5);
}

static pt * alloc_pt (int n)
{
  pt *p = (pt *) malloc (sizeof (pt) + (n - 1) * sizeof (double));
  if (p){
      p->minimized = 0;
      p->closest_pt_d = HUGE_VAL;
      p->closest_lm_d = HUGE_VAL;
    }
  return p;
}

static void
get_minf (mlsl_data * d, double *minf, double *p)
{
  rb_node *node = rb_tree_min (&d->pts);
  if (node) {
      *minf = ((pt *) node->k)->f;
      memcpy (p, ((pt *) node->k)->p, sizeof (double) * d->m);
  }
  node = rb_tree_min (&d->lms);

  if (node && node->k[0] < *minf) {
      *minf = node->k[0];
      memcpy (p, node->k + 1, sizeof (double) * d->m);
  }
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define MLSL_SIGMA 2.		/* MLSL sigma parameter, using value from the papers */
#define MLSL_GAMMA 0.3		/* MLSL gamma parameter (FIXME: what should it be?) */


/* returns the squared L2 norm of |x-y|; if x==NULL it is assumed 0.
 * Note: faster version exists as LEVMAR_L2NRMXMY
 */
static double L2sqnrm(const double *x, const double *y, int n)
{
register int i;
double sumsq;

  sumsq=0.0;
  if(x)
    for(i=0; i<n; ++i){
      double e;

      e=x[i]-y[i];
      sumsq+=e*e;
    }
  else
    for(i=0; i<n; ++i){
      sumsq+=y[i]*y[i];
    }

  return sumsq;
}

static int mlsl_initialized=0; // has random gen for MLSL been initialized?

/* L-M embedded in a MLSL global optimization scheme.
 * Arguments up to adata are as in dlevmar_bc_der()
 *
 * minf: error at the computed global minimum
 * Nsamples: #samples/iteration (0=default)
 * lds: use random or low-discrepancy seq. (lds)
 */
int mlsl_dlevmar_der(
    void (*func)(double *p, double *hx, int m, int n, void *adata),
    void (*jacf)(double *p, double *j, int m, int n, void *adata),
    double *p, /* in: initial guess, out: minimizer */
    const double *x, int m, int n, 
    const double *lb, const double *ub,	/* bounds */
    const double *scl, int itmax, double *opts,
    double *info, double *work, double *covar, void *adata,
    double *minf, int Nsamples,
    int lds, int verbose)
{
int ret = NLOPT_SUCCESS;
mlsl_data d;
register int i, j;
pt *ppt;
double *hx=NULL, locinfo[LM_INFO_SZ];

  /* make sure that the random genenerator is initialized only the first time mlsl_dlevmar_der() is called */
  if(!mlsl_initialized){
    mlsl_initialized=1;

#ifdef MLSL_REPEATABLE_RANDOM
    nlopt_init_genrand(5489UL); /* a default initial seed */
#else
    nlopt_init_genrand((unsigned long)time(NULL));
#endif
  }

  if(Nsamples < 1)
    d.N = 4;			/* FIXME: what is good number of samples per iteration? */
  else
    d.N = Nsamples;

  if(!info) info=locinfo; // ensure non-zero info below

  d.m = m;
  d.n = n;
  hx = (double *) malloc (sizeof (double) * n);
  if (!hx) {
    ret = NLOPT_OUT_OF_MEMORY;
    goto done;
  }
  d.lb = lb;
  d.ub = ub;
  d.minf_max = -1.0; // constraint not enforced, CHECKME
  d.maxeval = 50000; //30000; // at most that many func evals (in combination with minsamples), CHECKME
  d.nevals = 0;
  d.minsamples = 250; //200; // at least that many samples considered, CHECKME
  d.nsamples = 0;
  rb_tree_init (&d.pts, pt_compare);
  rb_tree_init (&d.lms, lm_compare);
  d.s = lds ? nlopt_sobol_create ((unsigned) m) : NULL;
  d.gamma = MLSL_GAMMA;

  d.R_prefactor = sqrt (2. / K2PI) * pow (gam (m) * MLSL_SIGMA, 1.0 / m);
  for (i = 0; i < m; ++i)
    d.R_prefactor *= pow (ub[i] - lb[i], 1.0 / m);

  /* MLSL also suggests setting some minimum distance from points
     to previous local minimiza and to the boundaries; I don't know
     how to choose this as a fixed number, so I set it proportional
     to R; see also the comments at top.  dlm and dbound are the
     proportionality constants. */
  d.dlm = 1.0;			/* min distance/R to local minima (FIXME: good value?) */
  d.dbound = 1e-6;		/* min distance/R to ub/lb boundaries (good value?) */

  ppt = alloc_pt (m);
  if (!ppt) {
      ret = NLOPT_OUT_OF_MEMORY;
      goto done;
    }

  /* FIXME: how many sobol points to skip, if any? */
  nlopt_sobol_skip (d.s, (unsigned) (10 * m + d.N), ppt->p);

  memcpy (ppt->p, p, m * sizeof (double));
  (*func)(p, hx, m, n, adata);
  ppt->f = L2sqnrm(x, hx, n);
  d.nevals++;
  if (!rb_tree_insert (&d.pts, (rb_key) ppt)) {
      free (ppt);
      ret = NLOPT_OUT_OF_MEMORY;
  }
  if (d.nevals>=d.maxeval)
    ret = NLOPT_MAXEVAL_REACHED;
  else if (ppt->f < d.minf_max)
    ret = NLOPT_MINF_MAX_REACHED;

  while (ret == NLOPT_SUCCESS) {
      rb_node *node;
      double R;

      get_minf (&d, minf, p);

      /* sampling phase: add random/quasi-random points */
      for (i = 0; i < d.N && ret == NLOPT_SUCCESS; ++i){
	      ppt = alloc_pt (m);
        if (!ppt) {
            ret = NLOPT_OUT_OF_MEMORY;
            goto done;
        }
        if (d.s)
          nlopt_sobol_next (d.s, ppt->p, lb, ub);
        else {			/* use random points instead of LDS */
            for (j = 0; j < m; ++j)
              ppt->p[j] = nlopt_urand (lb[j], ub[j]);
        }
        d.nsamples++;
        (*func)(ppt->p, hx, m, n, adata);
        ppt->f = L2sqnrm(x, hx, n);
        d.nevals++;
        if (!rb_tree_insert (&d.pts, (rb_key) ppt)) {
            free (ppt);
            ret = NLOPT_OUT_OF_MEMORY;
        }
        if (d.nsamples>=d.minsamples && d.nevals>=d.maxeval)
          ret = NLOPT_MAXEVAL_REACHED;
        else if (ppt->f < d.minf_max)
          ret = NLOPT_MINF_MAX_REACHED;
        else {
            find_closest_pt (m, &d.pts, ppt);
            find_closest_lm (m, &d.lms, ppt);
            pts_update_newpt (m, &d.pts, ppt);
        }
	    }

      /* distance threshhold parameter R in MLSL */
      R = d.R_prefactor * pow (log ((double) d.pts.N) / d.pts.N, 1.0 / m);

      /* local search phase: do local opt. for promising points */
      node = rb_tree_min (&d.pts);
      for (i = (int) (ceil (d.gamma * d.pts.N) + 0.5); node && i > 0 && ret == NLOPT_SUCCESS; --i) {
	      ppt = (pt *) node->k;
        if (is_potential_minimizer (&d, ppt, R, d.dlm * R, d.dbound * R)) {
            int lret;
            double *lm;

            if (d.nsamples>=d.minsamples && d.nevals>=d.maxeval) {
              ret = NLOPT_MAXEVAL_REACHED;
              break;
            }
            lm = (double *) malloc (sizeof (double) * (m + 1));
            if (!lm) {
              ret = NLOPT_OUT_OF_MEMORY;
              goto done;
            }
            memcpy(lm + 1, ppt->p, sizeof (double) * m);
            // unconstrained
            //lret=dlevmar_der(func, jacf, lm + 1, (double *)x, m, n, itmax, opts, info, work, NULL, adata);

            // with box constraints
            lret=dlevmar_bc_der(func, jacf, lm + 1, (double *)x, m, n, (double *)lb, (double *)ub, (double *)scl, itmax, opts, info, work, NULL, adata);
            //printf("local search: "); for(j=1; j<=m; ++j) printf("%g ", lm[j]); putchar('\n');
            d.nevals+=(int)info[7]; // add function evaluations by levmar
#if 0
            (*func)(lm + 1, hx, m, n, adata);
            lm[0]=L2sqnrm(x, hx, n);
#else
            // avoid function call
            lm[0]=info[1];
#endif
            //printf("\t%g\n", lm[0]); // sum of squared residuals for the minimizer in lm[1:m]
            ppt->minimized = 1;
            if (lret == LM_ERROR) {
              free (lm);
              ret = NLOPT_FAILURE;
              goto done;
            }
            if (!rb_tree_insert (&d.lms, lm)) {
              free (lm);
              ret = NLOPT_OUT_OF_MEMORY;
            }
            else if (*lm < d.minf_max)
              ret = NLOPT_MINF_MAX_REACHED;
            else if (d.nsamples>=d.minsamples && d.nevals>=d.maxeval)
              ret = NLOPT_MAXEVAL_REACHED;
            else
              pts_update_newlm (m, &d.pts, lm);
        }

	  /* TODO: additional stopping criteria based
	     e.g. on improvement in function values, etc? */

	      node = rb_tree_succ (node);
	    }
    }

  get_minf (&d, minf, p);

  /* compute the covariance from the Jacobian at p */
  if(covar){
    double *jac;
    /* declarations from levmar/misc.h (this is ugly!) */
    extern void dlevmar_trans_mat_mat_mult(double *a, double *b, int n, int m);
    extern int dlevmar_covar(double *JtJ, double *C, double sumsq, int m, int n);

    if(!(jac=(double *)malloc(sizeof(double)*n*m))){
      ret=NLOPT_OUT_OF_MEMORY;
      goto done;
    }
    (*jacf)(p, jac, m, n, adata);
    dlevmar_trans_mat_mat_mult(jac, covar, n, m); /* covar = J^T J */
    dlevmar_covar(covar, covar, *minf, m, n);
    free(jac);
  }

done:
  if(hx) free(hx);

  nlopt_sobol_destroy (d.s);
  rb_tree_destroy_with_keys (&d.lms);
  rb_tree_destroy_with_keys (&d.pts);

  if(verbose) printf("MLSL return reason %d, evals %d, samples %d, minf %g\n", ret, d.nevals, d.nsamples, *minf);

  return ret;
}
