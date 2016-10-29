#ifndef _SOBOL_H
#define _SOBOL_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sobol' low-discrepancy-sequence generation */

typedef struct nlopt_soboldata_s *nlopt_sobol;
extern nlopt_sobol nlopt_sobol_create(unsigned sdim);
extern void nlopt_sobol_destroy(nlopt_sobol s);
extern void nlopt_sobol_next01(nlopt_sobol s, double *x);
extern void nlopt_sobol_next(nlopt_sobol s, double *x,
			    const double *lb, const double *ub);
extern void nlopt_sobol_skip(nlopt_sobol s, unsigned n, double *x);

#ifdef __cplusplus
}
#endif

#endif /* _SOBOL_H */
