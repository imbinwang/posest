#ifndef _MT19937AR_H
#define _MT19937AR_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/* pseudorandom number generation by Mersenne twister algorithm */
extern void nlopt_init_genrand(unsigned long s);
extern double nlopt_urand(double a, double b);
extern int nlopt_iurand(int n);
extern double nlopt_nrand(double mean, double stddev);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
