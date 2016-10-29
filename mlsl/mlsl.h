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

#ifndef MLSL_H
#define MLSL_H

#define NLOPT_FAILURE  -1 /* generic failure code */
#define NLOPT_INVALID_ARGS  -2
#define NLOPT_OUT_OF_MEMORY  -3
#define NLOPT_ROUNDOFF_LIMITED  -4
#define NLOPT_FORCED_STOP  -5
#define NLOPT_SUCCESS  1 /* generic success code */
#define NLOPT_STOPVAL_REACHED  2
#define NLOPT_FTOL_REACHED  3
#define NLOPT_XTOL_REACHED  4
#define NLOPT_MAXEVAL_REACHED  5
#define NLOPT_MAXTIME_REACHED  6
#define NLOPT_MINF_MAX_REACHED NLOPT_STOPVAL_REACHED


#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

int mlsl_dlevmar_der(
      void (*func)(double *p, double *hx, int m, int n, void *adata),
      void (*jacf)(double *p, double *j, int m, int n, void *adata),
      double *p, const double *x, int m, int n,
      const double *lb, const double *ub, const double *scl,
      int itmax, double *opts, double *info, double *work, double *covar, void *adata,
      double *minf, int Nsamples,
      int lds, int verbose);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif

