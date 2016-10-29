/////////////////////////////////////////////////////////////////////////////////
// 
//  Reentrant random number generators (RNGs)
//
//  Copyright (C) 2002 - 2013  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#if defined(_WIN32) && defined(_MSC_VER) // Visual Studio
#include <process.h>
#define GETPID  _getpid
#else // !defined(_WIN32) || !defined(_MSC_VER), e.g. GCC, Mingw, ...
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#define GETPID  getpid
#endif /* _WIN32 */

#include "rngs.h"


/*
   32-bits Random number generator U[0, 2^32): lfsr113
   Author: Pierre L'Ecuyer,
   Source: http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps
   ---------------------------------------------------------
*/

/**** VERY IMPORTANT **** :
  The initial seeds z1, z2, z3, z4  MUST be larger than
  1, 7, 15, and 127 respectively.
****/

#define _Zx_SEED 987654321

void rng_lfsr113_seed(struct rng_lfsr113_state *state)
{
//printf("seeding lfsr113...\n");
  state->z1=_Zx_SEED;
  state->z2=_Zx_SEED;
  state->z3=_Zx_SEED;
  state->z4=_Zx_SEED;
}

uint32_t rng_lfsr113(struct rng_lfsr113_state *state)
{
uint32_t b;
uint32_t z1, z2, z3, z4;

  z1=state->z1; z2=state->z2;
  z3=state->z3; z4=state->z4;

  b =((z1 << 6) ^ z1) >> 13;
  z1=((z1 & 4294967294UL) << 18) ^ b;
  b =((z2 << 2) ^ z2) >> 27;
  z2=((z2 & 4294967288UL) << 2) ^ b;
  b =((z3 << 13) ^ z3) >> 21;
  z3=((z3 & 4294967280UL) << 7) ^ b;
  b =((z4 << 3) ^ z4) >> 12;
  z4=((z4 & 4294967168UL) << 13) ^ b;

   /* copy back */
  state->z1=z1; state->z2=z2;
  state->z3=z3; state->z4=z4;

  return (z1 ^ z2 ^ z3 ^ z4); // * 2.3283064365386963e-10; // do not convert to U[0, 1)
}


/* Marsaglia's KISS generator (with multiply-with-carry) */
void rng_kiss32_seed(struct rng_kiss32_state *state)
{
//printf("seeding kiss32...\n");
  state->mwc_up=362436069U;
  state->mwc_lo=521288629U;
  state->shr3  =362436000U;
  state->cong  =380116160U;
}

uint32_t rng_kiss32(struct rng_kiss32_state *state)
{
uint32_t mwc_up, mwc_lo, shr3, cong;
uint32_t mwc;

  mwc_up=state->mwc_up;
  mwc_lo=state->mwc_lo;
  shr3=state->shr3;
  cong=state->cong;

  /* multiply-with-carry generator */
  mwc_up=36969*(mwc_up&65535)+(mwc_up>>16);
  mwc_lo=18000*(mwc_lo&65535)+(mwc_lo>>16);
  mwc=(mwc_up<<16)+mwc_lo;
  //mwc=(mwc_up<<16) + (mwc_up>>16) + mwc_low;

  /* 3-shift shift-register generator */
  shr3^=(shr3<<13);
  shr3^=(shr3>>17);
  shr3^=(shr3<<5);

  /* congruential generator */
  cong=69069U*cong+12345U;

  /* copy back */
  state->mwc_up=mwc_up;
  state->mwc_lo=mwc_lo;
  state->shr3=shr3;
  state->cong=cong;

  return (mwc^cong)+shr3;
}


/*  Marsaglia's add-with-carry generator (no multiply); called JKISS32() by Jones */
void rng_kiss32a_seed(struct rng_kiss32a_state *state)
{
//printf("seeding kiss32a...\n");
  state->x=123456789;
  state->y=234567891;
  state->z=345678912;
  state->w=456789123;
  state->c=0;
}

uint32_t rng_kiss32a(struct rng_kiss32a_state *state)
{
int t;
uint32_t x, y, z, w, c;

  x=state->x;
  y=state->y;
  z=state->z;
  w=state->w;
  c=state->c;

  y^=(y<<5); y^=(y>>7); y^=(y<<22); 
  t=z+w+c; z=w; c=t<0; w=t&2147483647;
  x += 1411392427;

  /* copy back */
  state->x=x;
  state->y=y;
  state->z=z;
  state->w=w;
  state->c=c;

  return x + y + w;
}

/* Jones JKISS generator */
void rng_jkiss_seed(struct rng_jkiss_state *state)
{
//printf("seeding jkiss...\n");
  state->x=123456789;
  state->y=987654321;
  state->z=43219876;
  state->c=6543217;
}

unsigned int rng_jkiss(struct rng_jkiss_state *state)
{
unsigned int x, y, z, c;
unsigned long long t;

  x=state->x;
  y=state->y;
  z=state->z;
  c=state->c;

  x=314527869 * x + 1234567;
  y ^= y << 5; y ^= y >> 7; y ^= y << 22;
  t=4294584393ULL * z + c; c=(unsigned int)(t >> 32); z=(unsigned int)t;

  /* copy back */
  state->x=x;
  state->y=y;
  state->z=z;
  state->c=c;

  return x + y + z;
}


static unsigned long int rng_seedgen(void);
static unsigned long int rng_rdtscrand(void);

/* initialize a RNG, also specifying whether a repeatable random sequence is desired or not */

#define __RSEED  rng_seedgen
//#define __RSEED  rng_rdtscrand
void rng_init(int howtoseed, struct rng_state *state)
{
register int i;
int n;

  if(state->seeded) return;
  state->seeded=1;

  /* seed */
#if RNG_USE_KISS32==1
  rng_kiss32_seed(state);
#elif RNG_USE_KISS32a==1
  rng_kiss32a_seed(state);
#elif RNG_USE_JKISS==1
    rng_jkiss_seed(state);
#elif RNG_USE_LFSR113==1
  rng_lfsr113_seed(state);
#else
#error Error: do not know how to initialize RNG in rng_init()!
#endif

  n=300;
  if(howtoseed==RNG_NONREPEATABLE){
    n+=(__RSEED()&511); // random int in [300, 812)
    //n+=(__RSEED()&1023); // random int in [300, 1324)
    //n+=(__RSEED()&2047); // random int in [300, 2348)
  }
 
  /* warm up */
  for(i=n; i-->0;  )
    rng_next(state);
}
#undef __RSEED

/*** misc functions ***/


#ifdef _MSC_VER
/* gettimeofday function for windows.
 * code adapted from:
 * http://www.openasthra.com/c-tidbits/gettimeofday-function-for-windows/
 */

#include <time.h>
#include <windows.h>

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone{
  int tz_minuteswest; /* minutes W of Greenwich */
  int tz_dsttime;     /* type of dst correction */
};

static int gettimeofday(struct timeval *tv, struct timezone *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres=0;
  static int tzflag;
  int aux;

  if(NULL!=tv){
    GetSystemTimeAsFileTime(&ft);

    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;

    /* converting file time to unix epoch */
    tmpres/=10;  /* convert into microseconds */
    tmpres-=DELTA_EPOCH_IN_MICROSECS;
    tv->tv_sec=(long)(tmpres / 1000000UL);
    tv->tv_usec=(long)(tmpres % 1000000UL);
  }

  if(NULL!=tz){
    if(!tzflag){
      _tzset();
      tzflag++;
    }
    _get_timezone(&aux);
    tz->tz_minuteswest=aux / 60;
    _get_daylight(&aux);
    tz->tz_dsttime=(aux>0);
  }

  return 0;
}
#endif /* _MSC_VER */


/* generate a seed using a hash function of time and PID
 * see Katzgrabber "Random Numbers in Scientific Computing: An Introduction"
 * http://arxiv.org/pdf/1005.4117v1.pdf
 */
static unsigned long int rng_seedgen(void)
{
unsigned long long int s;
unsigned long int seed, pid;

#if 0
  time_t seconds;

  s=(long long int)time(&seconds); /* get seconds since the Epoch */
#else
  struct timeval tv;

  gettimeofday(&tv, NULL);
  s=tv.tv_sec*1000000 + tv.tv_usec; /* get microseconds since the Epoch */
#endif
  pid=GETPID(); /* get process ID */

  seed=abs((int)(((s*181)*((pid-83)*359))%104729));

  return seed;
}

/* return random integer from /dev/urandom
 *
 * Returns 0 on error (i.e., /dev/urandom cannot be read)
 */
static unsigned long int rng_devrand(void)
{
FILE *fp;
unsigned long int r;

  if((fp=fopen("/dev/urandom", "rb"))==NULL){
    fprintf(stderr, "opening /dev/urandom failed!\n");
    return 0;
  };

  if(fread(&r, sizeof(unsigned long int), 1, fp)!=1)
    r=0;

  fclose(fp);

  return r;
}

/* random integer via the processor time stamp counter.
 *
 * Notice that the high and low bits are XORed to yield a 32 bit int
 */
unsigned long int rng_rdtscrand()
{
# ifdef _MSC_VER
unsigned long long i64 = __rdtsc();

  return (unsigned long int)((i64>>32) ^ (i64&0x00000000FFFFFFFFULL));

# else // everything else

unsigned int lo, hi;

  __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
  return (unsigned long int)(hi ^ lo);

# endif /* _MSC_VER */
}

/* standard uniform and Gaussian distributions
 * see Jones "Good Practice in (Pseudo) Random Number Generation for Bioinformatics Applications"
 * http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf
 */

/* generate a random double in [0, 1) */
double rng_stduniform(struct rng_state *state)
{
double x;
unsigned int a, b;

  a=rng_next(state) >> 6; /* Upper 26 bits */  
  b=rng_next(state) >> 5; /* Upper 27 bits */
  x=(a * 134217728.0 + b) / 9007199254740992.0;

  return x;
}

/* generate Gaussian deviate with mean 0 and stdev 1 (aka standard/unit normal distribution) */
double rng_stdnormal(struct rng_state *state)
{
double x, y, r;

  /* Box-Mueller transform */
  do{
    x=2.0*rng_stduniform(state) - 1.0;
    y=2.0*rng_stduniform(state) - 1.0;
    r=x*x + y*y;
  }while (r==0.0 || r>=1.0);

  r=sqrt((-2.0*log(r))/r);

  return x*r;
}

/* shuffle an array of n elements using the Fisher–Yates (a.k.a. Knuth) shuffle */
void rng_shuffle(struct rng_state *state, int arr[], int n)
{
register int i, j;
register int tmp;
const int repeatable=1;

  /* rng_init() will do nothing if seeded already */
  rng_init(repeatable? RNG_REPEATABLE : RNG_NONREPEATABLE, state);

  for(i=n; i-->1;  ){  // n−1 downto 1 
    j=rng_rint(state, i+1); // random int in [0, i]
    tmp=arr[j];
    arr[j]=arr[i];
    arr[i]=tmp;
  }
}


#if 0
/* This is a test main program */
int main (void)
{
struct rng_state state={0};

  rng_init(RNG_REPEATABLE, &state);
  rng_next(&state);

  return 0;
}
#endif
