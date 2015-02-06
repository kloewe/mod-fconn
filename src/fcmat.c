/*----------------------------------------------------------------------
  File    : fcmat.c
  Contents: data type for functional connectivity matrix
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
#ifndef _WIN32                  /* if Linux/Unix system */
#define _POSIX_C_SOURCE 200809L /* needed for clock_gettime() */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <time.h>
#include <assert.h>
#include <math.h>

#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __AVX__
#include <immintrin.h>
#endif

#include "cpuinfo.h"
#include "binarize.h"
#include "pcc.h"
#include "clamp.h"
#include "tetracc.h"
#include "fcmat.h"
// #include "stats.h"
// #include "mex.h"

// #define FORCENAIVE
#define PAIRSPLIT
#define RECTGRID
#define SAFETHREAD

/*----------------------------------------------------------------------
  Data Type Definition / Recursion Handling
----------------------------------------------------------------------*/
#ifdef REAL                     /* if REAL is defined, */
#  undef  _FCM_PASS             /* ensure _FCM_PASS is undefined */
#  define _FCM_PASS 0           /* define macro for single pass */
#  ifndef SUFFIX                /* function names get no suffix */
#  define SUFFIX                /* (only single set of functions) */
#  endif
#elif !defined _FCM_PASS        /* if in first pass of two */
#  undef  _FCM_PASS             /* ensure _FCM_PASS is undefined */
#  define _FCM_PASS 1           /* define macro for first pass */
#  define REAL      float       /* first pass: single precision */
#  define SUFFIX    _flt        /* function name suffix is '_flt' */
#else                           /* if in second pass of two */
#  undef  _FCM_PASS             /* ensure _FCM_PASS is undefined */
#  define _FCM_PASS 2           /* define macro for second pass */
#  define REAL      double      /* second pass: double precision */
#  define SUFFIX    _dbl        /* function name suffix is '_dbl' */
#endif

/*--------------------------------------------------------------------*/
#define float  1                /* to check the definition of REAL */
#define double 2

#if   REAL==float               /* if single precision data */
#undef  REAL_IS_DOUBLE
#define REAL_IS_DOUBLE  0       /* clear indicator for double */
#elif REAL==double              /* if double precision data */
#undef  REAL_IS_DOUBLE
#define REAL_IS_DOUBLE  1       /* set   indicator for double */
#else
#error "REAL must be either 'float' or 'double'"
#endif

#undef float                    /* delete definitions */
#undef double                   /* used for type checking */
/*--------------------------------------------------------------------*/
#define int         1           /* to check definitions */
#define long        2           /* for certain types */
#define ptrdiff_t   3

#if   DIM==int
#ifndef DIM_FMT
#define DIM_FMT     "d"         /* printf format code for int */
#endif
#ifndef strtodim
#define strtodim(s,p)   (int)strtol(s,p,0)
#endif

#elif DIM==long
#ifndef DIM_FMT
#define DIM_FMT     "ld"        /* printf format code for long */
#endif
#ifndef strtodim
#define strtodim(s,p)   strtol(s,p,0)
#endif

#elif DIM==ptrdiff_t
#ifndef DIM_FMT
#  ifdef _MSC_VER
#  define DIM_FMT   "Id"        /* printf format code for ptrdiff_t */
#  else
#  define DIM_FMT   "td"        /* printf format code for ptrdiff_t */
#  endif                        /* MSC still does not support C99 */
#endif
#ifndef strtodim
#define strtodim(s,p)   (ptrdiff_t)strtoll(s,p,0)
#endif

#else
#error "DIM must be either 'int', 'long', or 'ptrdiff_t'"
#endif

#undef int                      /* remove preprocessor definitions */
#undef long                     /* needed for the type checking */
#undef ptrdiff_t
/*--------------------------------------------------------------------*/

#ifndef SFXNAME                 /* macros to generate function names */
#define SFXNAME(n)      SFXNAME_1(n,SUFFIX)
#define SFXNAME_1(n,s)  SFXNAME_2(n,s)
#define SFXNAME_2(n,s)  n##s    /* the two step recursion is needed */
#endif                          /* to ensure proper expansion */

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define R2Z_MAX     18.3684002848385504    /* atanh(1-epsilon) */
#define TILE_MIN    16          /* minimum tile size for comp. */

#ifdef FORCENAIVE               /* if to enforce the naive version */
#  define INIT_PCC SFXNAME(init_naive)
#  define PAIR_PCC SFXNAME(pair_naive)
#elif defined __AVX__           /* if AVX instructions available */
#  define INIT_PCC SFXNAME(init_avx)
#  define PAIR_PCC SFXNAME(pair_avx)
#elif defined __SSE2__          /* if SSE2 instructions available */
#  define INIT_PCC SFXNAME(init_sse2)
#  define PAIR_PCC SFXNAME(pair_sse2)
#else                           /* if neither extension available */
#  define INIT_PCC SFXNAME(init_naive)
#  define PAIR_PCC SFXNAME(pair_naive)
#endif

#if defined __POPCNT__ && defined __SSE4_1__
#  define PCAND_TCC pcand_m128i /* use m128i implementation */
#  define BPI       128         /* 128 bits per integer */
#else
#  define PCAND_TCC pcand_lut16 /* use lut16 implemenation */
#  define BPI       32          /* 32 bits per integer */
#endif

#define INDEX(i,j,N)    ((size_t)(i)*((size_t)(N)+(size_t)(N) \
                        -(size_t)(i)-3)/2-1+(size_t)(j))

#ifndef GET_THREAD              /* if not yet defined */
#if USERTHD                     /* if to respect user flag */
#define GET_THREAD(var) ((var) & PCC_THREAD)
#else                           /* copy the threading flag */
#define GET_THREAD(var) PCC_THREAD
#endif                          /* otherwise force threading flag */
#endif

#ifndef THREAD_OK               /* if not yet defined */
#ifdef _WIN32                   /* if Microsoft Windows system */
#define THREAD_OK       0       /* return value is DWORD */
#define WORKERDEF(n,p)  DWORD WINAPI SFXNAME(n) (LPVOID p)
#else                           /* if Linux/Unix system */
#define THREAD_OK       NULL    /* return value is void* */
#define WORKERDEF(n,p)  void*        SFXNAME(n) (void* p)
#endif                          /* definition of a worker function */
#endif

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- thread worker data --- */
  DIM ra, rb;                   /* row    index range */
  DIM ca, cb;                   /* column index range */
  DIM cm;                       /* reference column for mirroring */
  SFXNAME(FCMAT)    *fcm;       /* underlying f.c. matrix object */
  SFXNAME(FCMGETFN) *get;       /* element computation function */
} SFXNAME(WORK);                /* (thread worker data) */

#ifndef WORKERTYPE              /* if not yet defined */
#define WORKERTYPE              /* define worker function type */
#ifdef _WIN32                   /* if Microsoft Windows system */
typedef DWORD WINAPI WORKER (LPVOID);
#else                           /* if Linux/Unix system */
typedef void*        WORKER (void*);
#endif                          /* worker for parallel execution */
#endif

/*----------------------------------------------------------------------
  Cache Filling Functions
----------------------------------------------------------------------*/

#define clamp       clamp1      /* choose clamping implementation */

/*--------------------------------------------------------------------*/

static inline REAL SFXNAME(r2z) (REAL r)
{                               /* --- Fisher's r to z transform */
  if (r <= (REAL)-1) return (REAL)-R2Z_MAX;
  if (r >= (REAL)+1) return (REAL)+R2Z_MAX;
  return (REAL)atanh(r);        /* compute tangens hyperbolicus */
}  /* r2z() */

/*--------------------------------------------------------------------*/

static REAL SFXNAME(pcc_pure) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- compute Pearson corr. coeff. */
  REAL r;                       /* correlation coefficient */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col > row) && (col < fcm->V));
  r = PAIR_PCC((REAL*)fcm->data +(size_t)row *(size_t)fcm->X,
               (REAL*)fcm->data +(size_t)col *(size_t)fcm->X,
               (int)fcm->T);    /* compute Pearson correlation coeff. */
  return SFXNAME(clamp)(r, (REAL)-1, (REAL)+1); /* clamp to [-1,1] */
}  /* pcc_pure() */

/*--------------------------------------------------------------------*/

static REAL SFXNAME(pcc_r2z) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get Pearson corr. coeff. */
  REAL r;                       /* correlation coefficient */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col > row) && (col < fcm->V));
  r = PAIR_PCC((REAL*)fcm->data +(size_t)row *(size_t)fcm->X,
               (REAL*)fcm->data +(size_t)col *(size_t)fcm->X,
               (int)fcm->T);    /* compute Pearson correlation coeff. */
  return SFXNAME(r2z)(r);       /* apply Fisher's r to z transform */
}  /* pcc_r2z() */

/*--------------------------------------------------------------------*/

static REAL SFXNAME(tcc_pure) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- compute tetrachoric corr. cf. */
  int n;                        /* number of 11 configurations */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col > row) && (col < fcm->V));
  n = PCAND_TCC((uint32_t*)fcm->data +(size_t)row *(size_t)fcm->X,
                (uint32_t*)fcm->data +(size_t)col *(size_t)fcm->X,
                (int)fcm->X);   /* count 11 configurations and */
  return fcm->cmap[n];          /* compute tetrachoric corr. coeff. */
}  /* tcc_pure() */

/*--------------------------------------------------------------------*/

static REAL SFXNAME(tcc_r2z) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- compute tetrachoric corr. cf. */
  int n;                        /* number of 11 configurations */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col > row) && (col < fcm->V));
  n = PCAND_TCC((uint32_t*)fcm->data +(size_t)row *(size_t)fcm->X,
                (uint32_t*)fcm->data +(size_t)col *(size_t)fcm->X,
                (int)fcm->X);   /* compute tetrachoric corr. coeff. */
  return SFXNAME(r2z)(fcm->cmap[n]);
}  /* tcc_r2z() */              /* apply Fisher's r to z transform */

/*--------------------------------------------------------------------*/
#ifdef PAIRSPLIT                /* --- split rectangle into 2 parts */

static void SFXNAME(rec_rct) (SFXNAME(WORK) *w,
                              DIM ra, DIM rb, DIM ca, DIM cb)
{                               /* --- process a rectangle */
  DIM i, j;                     /* loop variables */

  assert(w                      /* check the function arguments */
  &&    (ra >= 0) && (ra <  w->fcm->V)
  &&    (rb > ra) && (rb <= w->fcm->V)
  &&    (ca >= 0) && (ca <  w->fcm->V)
  &&    (cb > ca) && (cb <= w->fcm->V));
  i = rb-ra; j = cb-ca;         /* compute row and column range */
  if ((i > TILE_MIN)            /* if at least one direction */
  ||  (j > TILE_MIN)) {         /* is larger than the minimum size */
    if (i > j) {                /* if there are more rows */
      i = (ra+rb)/2;            /* split the tile horizontally */
      SFXNAME(rec_rct)(w, ra, i, ca, cb);
      SFXNAME(rec_rct)(w, i, rb, ca, cb); }
    else {                      /* if there are more columns */
      j = (ca+cb)/2;            /* split the tile vertically */
      SFXNAME(rec_rct)(w, ra, rb, ca, j);
      SFXNAME(rec_rct)(w, ra, rb, j, cb);
    } }                         /* process the parts recursively */
  else {                        /* if no larger than minimum size */
    REAL *cache = w->fcm->cache;/* get the cache array */
    DIM  r = w->fcm->ra;        /* and the coordinates */
    DIM  c = w->fcm->ca;        /* of the reference element */
    for (i = ra; i < rb; i++) { /* traverse the rows */
      for (j = ca; j < cb; j++) /* traverse the columns */
        cache[(size_t)(i-r) *(size_t)w->fcm->tile +(size_t)(j-c)]
          = w->get(w->fcm, i, j);
    }                           /* compute correlation coefficient */
  }                             /* and store it in the cache */
}  /* rec_rct() */

/*--------------------------------------------------------------------*/
#else                           /* --- split rectangle into 4 parts */

static void SFXNAME(rec_rct) (SFXNAME(WORK) *w,
                              DIM ra, DIM rb, DIM ca, DIM cb)
{                               /* --- process a rectangle */
  assert(w                      /* check the function arguments */
  &&    (ra >= 0) && (ra <  w->fcm->V)
  &&    (rb > ra) && (rb <= w->fcm->V)
  &&    (ca >= 0) && (ca <  w->fcm->V)
  &&    (cb > ca) && (cb <= w->fcm->V));
  if (cb-ca > TILE_MIN) {       /* if larger than minimum size */
    DIM i = (ra+rb)/2;          /* halven the tile size and */
    DIM j = (ca+cb)/2;          /* process parts recursively */
    SFXNAME(rec_rct)(w, ra, i, ca, j);
    SFXNAME(rec_rct)(w, ra, i, j, cb);
    SFXNAME(rec_rct)(w, i, rb, j, cb);
    SFXNAME(rec_rct)(w, i, rb, ca, j); }
  else {                        /* if no larger than minimum size */
    REAL *cache = w->fcm->cache;/* get the cache array */
    DIM  i, r = w->fcm->ra;     /* and the coordinates */
    DIM  j, c = w->fcm->ca;     /* of the reference element */
    for (i = ra; i < rb; i++) { /* traverse the rows */
      for (j = ca; j < cb; j++) /* traverse the columns */
        cache[(size_t)(i-r) *(size_t)w->fcm->tile +(size_t)(j-c)]
          = w->get(w->fcm, i, j);
    }                           /* compute correlation coefficient */
  }                             /* and store it in the cache */
}  /* rec_rct() */

#endif
/*--------------------------------------------------------------------*/

static void SFXNAME(rec_trg) (SFXNAME(WORK) *w, DIM a, DIM b)
{                               /* --- recursively process a triangle */
  assert(w                      /* check the function arguments */
  &&    (a >= 0) && (a <= w->fcm->V)
  &&    (b >= a) && (b <= w->fcm->V));
  if (b-a > TILE_MIN) {         /* if larger than minimal size */
    DIM i = (a+b)/2;            /* halven the area size and */
    SFXNAME(rec_trg)(w, a, i);  /* split into three parts */
    SFXNAME(rec_rct)(w, a, i, i, b);
    SFXNAME(rec_trg)(w, i, b); }
  else {                        /* if no larger than min. tile size */
    REAL *cache = w->fcm->cache;/* get the cache array */
    DIM  i, r = w->fcm->ra;     /* and the coordinates */
    DIM  j, c = w->fcm->ca;     /* of the reference element */
    for (i = a; i < b; i++) {   /* traverse the rows */
      for (j = i+1; j < b; j++) /* traverse the columns */
        cache[(size_t)(i-r) *(size_t)w->fcm->tile +(size_t)(j-c)]
          = w->get(w->fcm, i, j);
    }                           /* compute the correlation coeff. */
  }                             /* and store it in the cache */
}  /* rec_trg() */

/*--------------------------------------------------------------------*/
#ifdef PAIRSPLIT                /* --- split rectangle into 2 parts */

static WORKERDEF(work_trg, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  DIM a;                        /* start index for second strip */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    if (w->ca > w->ra)          /* process rectangle part of strip */
      SFXNAME(rec_rct)(w, w->ra, w->ca, w->ca, w->cb);
    SFXNAME(rec_trg)(w, w->ca, w->cb); /* process triangle part */
    if (w->ca > w->cm/2) break; /* if second strip done, abort */
    a     = w->cm -w->cb;       /* get start of opposite stripe */
    if (w->cb > a)       break; /* if no opposite strip, abort */
    w->cb = w->cm -w->ca;       /* get start and end index */
    w->ca = a;                  /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* work_trg() */

/*--------------------------------------------------------------------*/

static WORKERDEF(work_rct, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  assert(p);                    /* check the function argument */
  SFXNAME(rec_rct)(w, w->ra, w->rb, w->ca, w->cb);
  return THREAD_OK;             /* return a dummy result */
}  /* work_rct() */

/*--------------------------------------------------------------------*/
#else                           /* --- split rectangle into 4 parts */

static WORKERDEF(work_trg, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  DIM a, b, k;                  /* loop variables */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    SFXNAME(rec_trg)(w, w->ca, w->cb);
    k = w->cb -w->ca;           /* compute leading triangle */
    for (a = w->ra; a < w->ca; a += k) {
      b = (a+k < w->ca) ? a+k : w->ca;
      SFXNAME(rec_rct)(w, a, b, w->ca, w->cb);
    }                           /* split the strip into squares */
    if (w->ca > w->cm/2) break; /* if second strip done, abort */
    a     = w->cm -w->cb;       /* get start of opposite stripe */
    if (w->cb > a)       break; /* if no opposite strip, abort */
    w->cb = w->cm -w->ca;       /* get start and end index */
    w->ca = a;                  /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* work_trg() */

/*--------------------------------------------------------------------*/

static WORKERDEF(work_rct, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  DIM a, b, k;                  /* loop variables */

  assert(p);                    /* check the function argument */
  k = w->rb -w->ra;             /* get the size of the rectangles */
  for (a = w->ca; a < w->cb; a += k) {
    b = (a+k < w->cb) ? a+k : w->cb;
    SFXNAME(rec_rct)(w, w->ra, w->rb, a, b);
  }                             /* split the strip into squares */
  return THREAD_OK;             /* return a dummy result */
}  /* work_rct() */

#endif
/*--------------------------------------------------------------------*/

static int SFXNAME(fcm_fill) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- fill the cache */
  #ifdef RECTGRID               /* if to split rectangle into grid */
  DIM    i, j, n;               /* loop variables for threads */
  DIM    x, y;                  /* number of grid rows/columns */
  DIM    dx, dy;                /* number of voxels per grid row/col. */
  #else                         /* if to split rectangle into stripes */
  DIM    i, n;                  /* loop variables for threads */
  #endif
  DIM    k;                     /* number of rows/columns */
  DIM    m;                     /* reference column for mirroring */
  THREAD *threads;              /* thread handles */
  WORKER *worker;               /* worker for parallel execution */
  SFXNAME(WORK) *w;             /* data   for worker thread */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  DWORD  thid;                  /* dummy for storing the thread id */
  #endif                        /* (not really needed here) */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row > col) {              /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  fcm->ra = row -(row % fcm->tile);
  fcm->rb = fcm->ra +fcm->tile; /* compute the new row    range */
  if (fcm->rb > fcm->V) fcm->rb = fcm->V;
  fcm->ca = col -(col % fcm->tile);
  fcm->cb = fcm->ca +fcm->tile; /* compute the new column range */
  if (fcm->cb > fcm->V) fcm->cb = fcm->V;

  w = fcm->work;                /* and the data for the workers */
  if (fcm->ra >= fcm->ca) {     /* if to process a triangle */
    worker = SFXNAME(work_trg); /* get the worker function */
    m =   fcm->cb +fcm->ca;     /* reference index for mirroring */
    k = ((fcm->cb -fcm->ca)/2 +(DIM)fcm->nthd-1) /(DIM)fcm->nthd;
    if (k <= 0) k = 1;          /* compute the number of voxels */
    for (n = 0; n < fcm->nthd; n++) {
      w[n].ra = fcm->ra;        /* traverse the threads */
      w[n].rb = fcm->rb;        /* note the row range and */
      w[n].ca = fcm->ca +n*k;   /* compute the column start index */
      if (w[n].ca >= m/2)       /* if beyond half of the columns, */
        break;                  /* all columns are provided for */
      w[n].cb = w[n].ca +k;     /* compute and store end index */
      if (w[n].cb >= m/2) w[n].cb = fcm->cb -n*k;
      w[n].cm = m;              /* ensure no duplicate computations */
    } }                         /* and note reference for mirroring */
  #ifdef RECTGRID               /* if to split rectangle into grid */
  else {                        /* if to process a rectangle */
    worker = SFXNAME(work_rct); /* get the worker function */
    x = fcm->gc;                /* get the number of grid  columns */
    k = fcm->cb -fcm->ca;       /* and the number of voxel columns */
    if (x > k) { x = k; dx = 1;}/* limit grid to voxel columns */
    else dx = (k +x-1) /x;      /* compute columns per grid column */
    y = fcm->gr;                /* get the number of grid  rows */
    k = fcm->rb -fcm->ra;       /* and the number of voxel rows */
    if (y > k) { y = k; dy = 1;}/* limit grid to voxel rows */
    else dy = (k +y-1) /y;      /* compute rows    per grid row */
    n  = 0;                     /* initialize the thread index */
    for (i = 0; i < x; i++) {   /* traverse the grid rows */
      for (j = 0; j < y; j++) { /* traverse the grid columns */
        w[n].ra = fcm->ra +j*dy;/* compute the row range */
        w[n].rb = w[n].ra +dy;  /* for the next grid cell */
        if (w[n].rb > fcm->rb) w[n].rb = fcm->rb;
        w[n].ca = fcm->ca +i*dx;/* compute the column range */
        w[n].cb = w[n].ca +dx;  /* for the next grid cell */
        if (w[n].cb > fcm->cb) w[n].cb = fcm->cb;
        n++;                    /* advance the thread index */
      }
    }                           /* initialize the grid cells */
  }                             /* to be processed by the threads */
  #else                         /* if to split rectangle into stripes */
  else {                        /* if to process a rectangle */
    worker = SFXNAME(work_rct); /* get the worker function */
    k = ((fcm->rb-fcm->ra) +fcm->nthd-1) /fcm->nthd;
    if (k <= 0) k = 1;          /* compute the number of voxels */
    for (n = 0; n < fcm->nthd; n++) {
      w[n].ca = fcm->ca;        /* traverse the threads */
      w[n].cb = fcm->cb;        /* note the column range and */
      w[n].ra = fcm->ra +n*k;   /* compute the row start index */
      if (w[n].ra >= fcm->rb)   /* if beyond the last row, */
        break;                  /* all rows are provided for */
      w[n].rb = w[n].ra +k;     /* compute and store end index */
      if (w[n].rb >  fcm->rb) w[n].rb = fcm->rb;
    }                           /* initialize the row/column ranges */
  }                             /* to be processed by the threads */
  #endif
  fcm->err = 0;                 /* clear the error status */
  if (n <= 1) {                 /* if there is only one thread, */
    worker(w); return 0; }      /* execute the worker directly */

  threads = fcm->thrds;         /* get the array of thread handles */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  for (i = 0; i < n; i++) {     /* traverse the threads */
    threads[i] = CreateThread(NULL, 0, worker, w+i, 0, &thid);
    if (!threads[i]) { fcm->err = -1; break; }
  }                             /* create a thread for each strip */
  WaitForMultipleObjects(i, threads, TRUE, INFINITE);
  while (--i >= 0)              /* wait for threads to finish, */
    CloseHandle(threads[i]);    /* then close all thread handles */
  #else                         /* if Linux/Unix system */
  for (i = 0; i < n; i++)       /* traverse the threads */
    if (pthread_create(threads+i, NULL, worker, w+i) != 0) {
      fcm->err = -1; break; }   /* create a thread for each strip */
  while (--i >= 0)              /* wait for threads to finish */
    pthread_join(threads[i], NULL);
  #endif                        /* (join threads with this one) */
  return fcm->err;              /* return the error status */
}  /* fcm_fill() */

/*----------------------------------------------------------------------
  Retrieval Functions
----------------------------------------------------------------------*/

static REAL SFXNAME(fcm_pccotf) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get Pearson cc. on the fly */
  REAL r;                       /* correlation coefficient */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return (REAL)+1;            /* always return +1.0 */
  if (row >  col) {             /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  r = PAIR_PCC((REAL*)fcm->data +(size_t)row*(size_t)fcm->X,
               (REAL*)fcm->data +(size_t)col*(size_t)fcm->X,
               (int)fcm->T);    /* compute Pearson correlation coeff. */
  return SFXNAME(clamp)(r, (REAL)-1, (REAL)+1); /* clamp to [-1,1] */
}  /* fcm_pccotf() */

/*--------------------------------------------------------------------*/

static REAL SFXNAME(fcm_pccr2z) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get Pearson cc. on the fly */
  REAL r;                       /* correlation coefficient */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return (REAL)+R2Z_MAX;      /* return atanh(1-epsilon) */
  if (row > col) {              /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  r = PAIR_PCC((REAL*)fcm->data +(size_t)row*(size_t)fcm->X,
               (REAL*)fcm->data +(size_t)col*(size_t)fcm->X,
               (int)fcm->T);    /* compute Pearson correlation coeff. */
  return SFXNAME(r2z)(r);       /* apply Fisher's r to z transform */
}  /* fcm_pccr2z() */

/*--------------------------------------------------------------------*/

static REAL SFXNAME(fcm_tccotf) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get tetrachoric cc. on the fly */
  int n;                        /* number of 11 configurations */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return (REAL)+1.0;          /* always return +1.0 */
  if (row > col) {              /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  n = PCAND_TCC((uint32_t*)fcm->data +(size_t)row*(size_t)fcm->X,
                (uint32_t*)fcm->data +(size_t)col*(size_t)fcm->X,
                (int)fcm->X);   /* count number of 11 configs. and */
  return fcm->cmap[n];          /* compute tetrachoric corr. coeff. */
}  /* fcm_tccotf() */

/*--------------------------------------------------------------------*/

static REAL SFXNAME(fcm_tccr2z) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get tetrachoric cc. on the fly */
  int n;                        /* number of 11 configurations */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return (REAL)+R2Z_MAX;      /* return atanh(1-epsilon) */
  if (row > col) {              /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  n = PCAND_TCC((uint32_t*)fcm->data +(size_t)row*(size_t)fcm->X,
                (uint32_t*)fcm->data +(size_t)col*(size_t)fcm->X,
                (int)fcm->X);   /* compute tetrachoric corr. coeff. */
  return SFXNAME(r2z)(fcm->cmap[n]);
}  /* fcm_tccr2z() */           /* apply Fisher's r to z transform */

/*--------------------------------------------------------------------*/

static REAL SFXNAME(fcm_cache) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get cached correlation coeff. */
  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return fcm->diag;           /* return a fixed value */
  if (row > col) {              /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  #ifdef SAFETHREAD             /* if to guarantee fill success */
  if (((row < fcm->ra) || (row >= fcm->rb)
  ||   (col < fcm->ca) || (col >= fcm->cb))
  &&  (SFXNAME(fcm_fill)(fcm, row, col) != 0)) {
    DIM n = fcm->nthd;          /* if outside cache, fill cache; */
    fcm->nthd = 1;              /* on failure retry with one thread */
    SFXNAME(fcm_fill)(fcm, row, col);
    fcm->nthd = n;              /* (a single thread cannot fail) */
  }                             /* and restore the number of threads */
  #else                         /* if to allow a thread failure */
  if (((row < fcm->ra) || (row >= fcm->rb)
  ||   (col < fcm->ca) || (col >= fcm->cb))
  &&  (SFXNAME(fcm_fill)(fcm, row, col) != 0))
    return -INFINITY;           /* if outside cache, fill cache */
  #endif
  row -= fcm->ra;               /* compute row and column in cache */
  col -= fcm->ca;               /* and retrieve correlation coeff. */
  return fcm->cache[(size_t)row *(size_t)fcm->tile +(size_t)col];
}  /* fcm_cache() */

/*--------------------------------------------------------------------*/

static REAL SFXNAME(fcm_full) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get corr.c. from full rep. */
  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return fcm->diag;           /* return a fixed value */
  return fcm->cache[(row > col) /* retrieve the correlation coeff. */
                    ? INDEX(col, row, fcm->V)
                    : INDEX(row, col, fcm->V)];
}  /* fcm_full() */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/

SFXNAME(FCMAT)* SFXNAME(fcm_create) (REAL *data, DIM V, DIM T,
                                     int mode, ...)
{                               /* --- create a func. connect. matrix */
  SFXNAME(FCMAT)    *fcm;       /* func. con. matrix to be created */
  SFXNAME(FCMGETFN) *get;       /* element computation function */
  SFXNAME(WORK)     *w;         /* to initialize the worker data */
  int     i;                    /* loop variable for threads */
  va_list args;                 /* list of variable arguments */
  size_t  k, z;                 /* loop variable, cache size */
  #ifdef RECTGRID               /* if to split rectangle into grid */
  DIM     g;                    /* loop variable for grid size */
  #endif

  assert(data && (V > 1) && (T > 1)); /* check the function arguments */
  fcm = (SFXNAME(FCMAT)*)malloc(sizeof(SFXNAME(FCMAT)));
  if (!fcm) return NULL;        /* allocate the base structure */
  fcm->V     = V;               /* note the number of voxels */
  fcm->T     = T;               /* and  the number of scans */
  fcm->X     = T;               /* default: data blocks like scans */
  fcm->tile  = 0;               /* default: compute on the fly */
  fcm->mode  = mode;            /* note the processing mode */
  fcm->nthd  = proccnt();       /* default: use all processors */
  fcm->data  = NULL;            /* clear normalized/binarized data, */
  fcm->cmap  = NULL;            /* cosine map and cache */
  fcm->cache = NULL;            /* for easier cleanup */
  fcm->thrds = NULL;
  fcm->work  = NULL;
  fcm->diag  = (REAL)((mode & FCM_R2Z) ? R2Z_MAX : 1.0);
  fcm->ra    = 0; fcm->rb = V;
  fcm->ca    = 0; fcm->cb = V;  /* init. the coordinate ranges */
  fcm->row   = fcm->col = 0;    /* and the coordinates and value */
  fcm->value = fcm->diag;       /* of the current element */
  fcm->err   = 0;               /* clear the error status */

  va_start(args, mode);         /* start variable arguments */
  if (mode & FCM_CACHE)         /* if to use a cached version */
    fcm->tile = va_arg(args, DIM); /* get the tile/cache size */
  assert((fcm->tile == 0) || (fcm->tile <= V));
  if (mode & FCM_THREAD)        /* if to use a threaded version */
    fcm->nthd = va_arg(args, int); /* get the number of threads */
  va_end(args);                 /* end variable arguments */
  if (fcm->nthd < 1) fcm->nthd = 1;

  mode &= FCM_CORR;             /* get the correlation type */
  if      (mode == FCM_PCC) {   /* if Pearson's correlation coeff. */
    #ifdef __AVX__              /* use AVX if possible */
    #if REAL_IS_DOUBLE          /* if to use double precision values */
    fcm->X = (((int)T+3) & ~3); /* process blocks with 4 numbers */
    #else                       /* if to use single precision values */
    fcm->X = (((int)T+7) & ~7); /* process blocks with 8 numbers */
    #endif
    #else                       /* otherwise fall back to SSE2 */
    #if REAL_IS_DOUBLE          /* if to use double precision values */
    fcm->X = (((int)T+1) & ~1); /* process blocks with 2 numbers */
    #else                       /* if to use single precision values */
    fcm->X = (((int)T+3) & ~3); /* process blocks with 4 numbers */
    #endif                      /* get the data block size */
    #endif                      /* allocate memory for norm.ed data */
    if (posix_memalign(&fcm->data, 32,
            (size_t)V *(size_t)fcm->X *sizeof(REAL))) {
      SFXNAME(fcm_delete)(fcm); return NULL; }
    // fcm->data = malloc((size_t)V*(size_t)fcm->X *sizeof(REAL));
    // if (!fcm->data) { SFXNAME(fcm_delete)(fcm); return NULL; }
    INIT_PCC(data, (int)V, (int)T, fcm->data, (int)fcm->X); }
  else if (mode == FCM_TCC) {   /* if tetrachoric correlation coeff. */
    #if defined __POPCNT__ \
    &&  defined __SSE4_1__      /* use 128 bit and popcnt if possible */
    fcm->X = 4 *(((int)T+127) >> 7);
    #else                       /* otherwise fall back to LUT16 */
    fcm->X = ((int)T+31) >> 5;  /* get block size of binarized data */
    #endif                      /* allocate memory for binarized data */
    fcm->data = SFXNAME(binarize)(data, (int)V, (int)T, BIN_MEDIAN,BPI);
    if (!fcm->data) { SFXNAME(fcm_delete)(fcm); return NULL; };
    fcm->cmap = SFXNAME(make_cmap)((int)T); /* create cosine map */
    if (!fcm->cmap) { SFXNAME(fcm_delete)(fcm); return NULL; };
    init_popcnt(); }            /* initialize bit count table */
  else {                        /* if unknown correlation variant */
    fprintf(stderr, "fcm_create: unknown correlation variant\n");
    free(fcm); return NULL;     /* print an error message */
  }                             /* and abort the function */

  if      (fcm->tile <= 0) {    /* if computation on the fly */
    if (!(fcm->mode & FCM_R2Z)) /* if pure correlation coefficients */
      fcm->get = (mode == FCM_PCC)
               ? SFXNAME(fcm_pccotf) : SFXNAME(fcm_tccotf);
    else {                      /* if to apply Fisher's r to z trans. */
      fcm->get = (mode == FCM_PCC)
               ? SFXNAME(fcm_pccr2z) : SFXNAME(fcm_tccr2z);
    } }                         /* get the element retrieval function */
  else if (fcm->tile < V) {     /* if to cache smaller areas */
    z = (size_t)fcm->tile *(size_t)fcm->tile;
    fcm->cache = (REAL*)  malloc(z *sizeof(REAL));
    fcm->thrds = (THREAD*)malloc((size_t)fcm->nthd *sizeof(THREAD));
    fcm->work  = w = malloc((size_t)fcm->nthd *sizeof(SFXNAME(WORK)));
    if (!fcm->cache || !fcm->thrds || !fcm->work)  {
      SFXNAME(fcm_delete)(fcm); return NULL; }
    if (!(fcm->mode & FCM_R2Z)) /* if to compute pure corr. coeffs. */
      get = (mode == FCM_PCC) ? SFXNAME(pcc_pure) : SFXNAME(tcc_pure);
    else                        /* if to apply Fisher's r to z trans. */
      get = (mode == FCM_PCC) ? SFXNAME(pcc_r2z)  : SFXNAME(tcc_r2z);
    for (i = 0; i < fcm->nthd; i++) {
      w[i].fcm = fcm;           /* store the f.c. matrix object */
      w[i].get = get;           /* and the function that */
    }                           /* computes a matrix element */
    #ifdef RECTGRID             /* if to split rectangle into grid */
    for (g = (DIM)floor(sqrt((double)fcm->nthd)); g > 1; g--)
      if (g *(fcm->nthd/g) == fcm->nthd)
        break;                  /* compute factors as close as */
    fcm->gc = g;                /* possible to the square root */
    fcm->gr = fcm->nthd /g;     /* for the rectangle grid */
    #endif
    fcm->get  = SFXNAME(fcm_cache);
    fcm->ra = fcm->rb = -1;     /* get element retrieval function and */
    fcm->ca = fcm->cb = -1; }   /* invalidate row and column range */
  else {                        /* if to cache the whole matrix */
    z = (size_t)V *(size_t)(V-1)/2; /* cache for upper triangle */
    fcm->cache = (REAL*)malloc(z *sizeof(REAL));
    if (!fcm->cache) { SFXNAME(fcm_delete)(fcm); return NULL; }
    if (mode == FCM_PCC)        /* if Pearson correlation cofficient */
      #ifdef FORCENAIVE
      SFXNAME(pccx)    (data, fcm->cache, (int)V, (int)T, PCC_NAIVE);
      #else
      SFXNAME(pccx)    (data, fcm->cache, (int)V, (int)T,
                        PCC_AUTO|PCC_THREAD, fcm->nthd);
      #endif
    else                        /* if tetrachoric correlation coeff. */
      #ifdef FORCENAIVE
      SFXNAME(tetraccx)(data, fcm->cache, (int)V, (int)T, TCC_LUT16);
      #else
      SFXNAME(tetraccx)(data, fcm->cache, (int)V, (int)T,
                        TCC_AUTO|TCC_THREAD, fcm->nthd);
      #endif
    fcm->get = SFXNAME(fcm_full);
    if (fcm->mode & FCM_R2Z)    /* if Fisher's r-to-z transform, */
      for (k = 0; k < z; k++)   /* transform the matrix elements */
        fcm->cache[k] = SFXNAME(r2z)(fcm->cache[k]);
  }                             /* set the element retrieval function */

  return fcm;                   /* return created FC matrix */
} /* fcm_create() */

/*--------------------------------------------------------------------*/

void SFXNAME(fcm_delete) (SFXNAME(FCMAT) *fcm)
{                               /* --- delete a func. connect. matrix */
  assert(fcm);                  /* check the function argument */
  if (fcm->work)  free(fcm->work);
  if (fcm->thrds) free(fcm->thrds);
  if (fcm->cache) free(fcm->cache);
  if (fcm->cmap)  free(fcm->cmap);
  free(fcm->data);              /* delete cache, cosine map, data, */
  free(fcm);                    /* and the base structure */
} /* fcm_delete() */

/*--------------------------------------------------------------------*/

int SFXNAME(fcm_first) (SFXNAME(FCMAT) *fcm)
{                               /* --- get first matrix element */
  assert(fcm);                  /* check the function argument */
  fcm->err = 0;                 /* clear the error status */
  fcm->row = 0; fcm->col = 1;   /* start with first off-diag. element */
  fcm->value = fcm->get(fcm, 0, 1);
  return (fcm->err) ? -1 : +1;   /* store the value and return status */
}  /* fcm_first() */

/*--------------------------------------------------------------------*/

int SFXNAME(fcm_next) (SFXNAME(FCMAT) *fcm)
{                               /* --- get next matrix element */
  assert(fcm);                  /* check the function argument */
  if (fcm->col >= fcm->V)       /* if traversal is already finished, */
    return 0;                   /* abort the function with failure */
  fcm->col += 1;                /* go to the next column */
  if (fcm->col >= fcm->cb) {    /* if at the end of a cache row, */
    fcm->col  = fcm->ca;        /* return to the start of a row */
    fcm->row += 1;              /* and go to the next row */
    if (fcm->col <= fcm->row)   /* if in the lower triangle, */
      fcm->col = fcm->row+1;    /* go to the upper triangle */
    if (fcm->col >= fcm->cb) {  /* if at the end of a strip */
      fcm->row = 0;             /* start a new strip */
      fcm->col = fcm->cb;       /* (row 0, first column after strip) */
      if (fcm->col >= fcm->V) return 0;
    }                           /* if whole matrix is traversed, */
  }                             /* abort the function with failure */
  fcm->value = fcm->get(fcm, fcm->row, fcm->col);
  return (fcm->err) ? -1 : +1;  /* store the value and return status */
}  /* fcm_next() */

/*--------------------------------------------------------------------*/

void SFXNAME(fcm_show) (SFXNAME(FCMAT) *fcm)
{                               /* --- show funct. connect. matrix */
  assert(fcm);                  /* check the function argument */
  printf("     ");              /* start the matrix output */
  for (DIM col = 0; col < fcm->V; col++)
    printf(" %-8"DIM_FMT, col); /* print the header line */
  printf("\n");                 /* (column indices) */
  for (DIM row = 0; row < fcm->V; row++) {
    printf("%-5"DIM_FMT, row);  /* print the row index */
    for (DIM col = 0; col < fcm->V; col++)
      printf(" %-8.3g", SFXNAME(fcm_get)(fcm, row, col));
    printf("\n");               /* traverse the matrix rows */
  }                             /* print the matrix elements */
}  /* fcm_show() */

/*----------------------------------------------------------------------
  Recursion Handling
----------------------------------------------------------------------*/
#if _FCM_PASS == 1              /* if in first of two passes */
#undef REAL
#undef SUFFIX
#include "fcmat.c"              /* process file recursively */
#endif

/*----------------------------------------------------------------------
  Main Function
----------------------------------------------------------------------*/
#ifdef FCMAT_MAIN
#undef REAL                     /* delete definitions */
#undef SUFFIX                   /* of precision selectors */
#if 1
#define REAL    float           /* test case: single precision */
#define SUFFIX  _flt            /* function name suffix is '_flt' */
#else
#define REAL    double          /* test case: double precision */
#define SUFFIX  _dbl            /* function name suffix is '_dbl' */
#endif

/*--------------------------------------------------------------------*/
#if 1

static double timer (void)
{                               /* --- get current time */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  return 0.001 *(double)timeGetTime();
  #else                         /* if Linux/Unix system */
  struct timespec tp;           /* POSIX time specification */
  clock_gettime(CLOCK_MONOTONIC, &tp);
  return (double)tp.tv_sec +1e-9 *(double)tp.tv_nsec;
  #endif                        /* return time in seconds */
}  /* timer() */

#endif
/*--------------------------------------------------------------------*/

int main (int argc, char* argv[])
{                               /* --- main function for testing */
  int  t;                       /* indicator for another element */
  DIM  V, T;                    /* number of voxels and time points */
  DIM  C;                       /* tile size for caching */
  long S;                       /* seed value for random numbers */
  REAL *data;                   /* data array */
  REAL *corr;                   /* correlation coefficients */
  DIM  r, c;                    /* loop variables (row and column) */
  REAL a, b;                    /* to compare correlation coeffs. */
  int  diff;                    /* indicator for a difference */
  double t0;                    /* timer for measurements */
  SFXNAME(FCMAT) *fcm;          /* functional connectivity matrix */

  if ((argc < 3) || (argc > 5)){/* check the program arguments */
    fprintf(stderr, "usage: %s V T [C [seed]]\n", argv[0]); return 0; }
  V = strtodim(argv[1], NULL);  /* get the number of voxels */
  T = strtodim(argv[2], NULL);  /* and the number of scans */
  if (argc <= 3) C = V;         /* get the tile size for caching */
  else C = strtodim(argv[3], NULL);
  if (argc <= 4) S = time(NULL);/* get the seed value */
  else S = strtol(argv[4], NULL, 0);
  srand((unsigned)S);           /* seed the random number generator */

  data = malloc((size_t)V *(size_t)T       *sizeof(REAL));
  if (!data) {                  /* allocate and generate data */
    fprintf(stderr, "%s: not enough memory\n", argv[0]); return -1; }
  for (size_t i = 0; i < (size_t)(T*V); i++)
    data[i] = (REAL)(rand()/((double)RAND_MAX+1));
  
  /* --- test correctness --- */
  corr = malloc((size_t)V *(size_t)(V-1)/2 *sizeof(REAL));
  if (!corr) {                  /* compute correlation coefficients */
    fprintf(stderr, "%s: not enough memory\n", argv[0]); return -1; }
  SFXNAME(pccx)(data, corr, (int)V, (int)T, PCC_AUTO); //PCC_NAIVE);
  
  printf("test #1 ... ");
  fcm = SFXNAME(fcm_create)(data, V, T, FCM_PCC|FCM_CACHE, C);
  if (!fcm) {                   /* create functional connect. matrix */
    fprintf(stderr, "%s: not enough memory\n", argv[0]); return -1; }
  diff = 0;                     /* initialize the difference flag */
  for (t = SFXNAME(fcm_first)(fcm); t > 0; t = SFXNAME(fcm_next)(fcm)) {
    r = SFXNAME(fcm_row)(fcm);  /* traverse the matrix elements */
    c = SFXNAME(fcm_col)(fcm);  /* and get their row and column */
    a = SFXNAME(fcm_value)(fcm);
    b = corr[INDEX(r,c,V)];     /* get the correlation coefficients */
    if (a == b) continue;       /* and compare them */
    printf("%6"DIM_FMT" %6"DIM_FMT": % 18.16f % 18.16f\n", r,c,a,b);
    diff = 1;                   /* print any difference and */
  }                             /* set the difference flag */
  if (t < 0) {                  /* check for a computation error */
    fprintf(stderr, "%s: thread error\n", argv[0]); return -1; }
  printf("%s\n", (diff) ? "[FAILED]" : "[PASSED]");
  SFXNAME(fcm_delete)(fcm);     /* delete the func. connect. matrix */
  
  #if 0
  printf("test #2 ... ");
  fcm = SFXNAME(fcm_create)(data, V, T, FCM_PCC|FCM_CACHE, C);
  if (!fcm) {                   /* create functional connect. matrix */
    fprintf(stderr, "%s: not enough memory\n", argv[0]); return -1; }
  diff = 0;                     /* initialize the difference flag */
  for (DIM i = 0; i < V; i++) { /* traverse rows and cols */
    for (DIM j = i+1; j < V; j++) {
      a = SFXNAME(fcm_get)(fcm,i,j);
      b = corr[INDEX(i,j,V)];   /* get the correlation coefficients */
      if (a == b) continue;     /* and compare them */
      printf("%6"DIM_FMT" %6"DIM_FMT": % 18.16f % 18.16f\n", r,c,a,b);
      diff = 1;                 /* print any difference and */
    }
  }
  printf("%s\n", (diff) ? "[FAILED]" : "[PASSED]");
  SFXNAME(fcm_delete)(fcm);     /* delete the func. connect. matrix */
  #endif
          
  /* --- test perfomance --- */
  int n = 1;
  
  t0 = timer();
  fcm = SFXNAME(fcm_create)(data, V, T, FCM_PCC|FCM_CACHE, C);
  if (!fcm) {                   /* create functional connect. matrix */
    fprintf(stderr, "%s: not enough memory\n", argv[0]); return -1; }

  for (t = SFXNAME(fcm_first)(fcm); t > 0; t = SFXNAME(fcm_next)(fcm)) {
    r = SFXNAME(fcm_row)(fcm);  /* traverse the matrix elements */
    c = SFXNAME(fcm_col)(fcm);  /* and get their row and column */
    a = SFXNAME(fcm_value)(fcm);
  }
  if (t < 0) {                  /* check for a computation error */
    fprintf(stderr, "%s: thread error\n", argv[0]); return -1; }
  SFXNAME(fcm_delete)(fcm);     /* delete the func. connect. matrix */
  printf("time (fcm_next): %8.2f\n", (timer()-t0)/(double)n);

  #if 0
  t0 = timer();
  fcm = SFXNAME(fcm_create)(data, V, T, FCM_PCC|FCM_CACHE, C);
  if (!fcm) {                   /* create functional connect. matrix */
    fprintf(stderr, "%s: not enough memory\n", argv[0]); return -1; }
  for (DIM i = 0; i < SFXNAME(fcm_dim)(fcm); i++) {    /* traverse rows */
    for (DIM j = i+1; j < SFXNAME(fcm_dim)(fcm); j++) {/* traverse cols */
      a = SFXNAME(fcm_get)(fcm,i,j);
    }
  }
  SFXNAME(fcm_delete)(fcm);     /* delete the func. connect. matrix */
  printf("time (fcm_get):  %8.2f\n", (timer()-t0)/(double)n);
  #endif  
  
  free(corr);                   /* delete correlation coefficients */
  free(data);                   /* delete the data */
  return 0;                     /* return 'ok' */
}  /* main() */

#endif  /* #ifdef FCMAT_MAIN */

#undef FCMAT_MAIN
