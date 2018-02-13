/*----------------------------------------------------------------------------
  File    : edgestats.c
  Contents: edge-level statistics
  Authors : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdarg.h>
#include <time.h>
#include <pthread.h>
#include <assert.h>
#include "cpuinfo.h"
#include "stats.h"
#include "fcmat.h"
#include "matrix.h"
#include "dot.h"
#include "edgestats.h"

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE)
#  include "mex.h"
#endif

#ifndef NDEBUG
#  line __LINE__ "edgestats.c"
#endif

/*----------------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------------*/
#define THREAD    pthread_t                 // use the POSIX thread type
#define THREAD_OK NULL                      // return value is void*

/*----------------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------------*/
typedef struct {                            // --- thread worker data ---
  FCMAT      **fcm;                         // FC matrices
  int        n;                             // number of FC matrices
  Func       *f;                            // function to be used
  const int  *an;                           // ptr to an array of n-vals
  const REAL *arg;                          // additional argument
  MATRIX     *mos;                          // result: matrix of statistics
  int        s, e;                          // index of start and end series
  int        err;                           // error indicator
} WORK;

typedef void* WORKER (void*);               // worker

/*----------------------------------------------------------------------------
  Auxiliary Functions
----------------------------------------------------------------------------*/
REAL corrv_w (const REAL* a, const int *n, const REAL *arg, REAL *buf)
{
  assert(a && n && arg && buf);

  REAL sqr = 0;                             // pre-normalize a
  REAL m = mean(a, *n);                     // (it's assumed that arg has been
  for (int k = 0; k < *n; k++) {            // pre-normalized by the caller)
    buf[k] = a[k] - m;
    sqr += buf[k]*buf[k]; }
  sqr = sqrt(sqr);
  sqr = (sqr > 0) ? 1/sqr : 0;
  for (int k = 0; k < *n; k++)
    buf[k] *= sqr;

  return dot(buf, arg, *n);                 // return corr. coefficient
}

/*----------------------------------------------------------------------------
  Functions
----------------------------------------------------------------------------*/

/* fcm_uni_single
 * --------------
 * compute statistics across functinal connectomes (single-threaded version)
 *
 * parameters
 * fcm   set of FC matrices
 * n     number of FC matrices
 * f     function to be used (see typedefs in edgestats.h)
 * an    ptr to an array of n-vals, will be passed to f as 2nd arg
 * arg   additional arg. that will be passed to f.u.f2 if relevant
 * mos   result: matrix of statistics
 *
 * returns
 * 0 on success
 */
static int fcm_uni_single(FCMAT **fcm, int n, Func f, const int *an,
                          const REAL *arg, MATRIX *mos)
{
  assert(fcm && (n > 0) && an && mos);

  // allocate memory for temp. array 1
  void *mem = malloc((size_t)n *sizeof(REAL) +31);
  if (!mem) {
    DBGMSG("ERROR: malloc failed");
    return -1; }                            // return 'failure'
  REAL *a = (REAL*)(((uintptr_t)mem +31) & ~(uintptr_t)31);

  if      (f.kind == 1) {                   // --- comp. stats using f1
    #if 1                                   // variant 1: use iterator
    int t = 0;                              // init status variable t
    for (int k = 0; k < n; k++)             // initiate iterations
      t += fcm_first(fcm[k]);
    if (t < 0) {                            // check status
      DBGMSG("ERROR: failure when calling fcm_first");
      return -1;                            // return 'failure'
    }
    while (t == 0) {                        // while iterations not completed
      for (int k = 0; k < n; k++)           // collect vals from all matrices
        a[k] = fcm_value(fcm[k]);           // in a and compute the statistic
      mat_set(mos, fcm_row(fcm[0]), fcm_col(fcm[0]), (*f.u.f1)(a, an));
      for (int k = 0; k < n; k++)           // for the current pair, then
        t += fcm_next(fcm[k]);              // move on to the next elements
    }
    free(mem);
    if (t < n) {                            // check status
      DBGMSG("ERROR: failure when calling fcm_next");
      return -1;                            // return 'failure'
    } }
    #else                                   // variant 2: use fcm_get
    int N = fcm_dim(*fcm);                  // get number of nodes
    for (int i = 0; i < N; i++) {           // traverse pairs of nodes
      for (int j = i+1; j < N; j++) {       // (upper triangular matrix)
        for (int k = 0; k < n; k++)         // collect vals from all matrices
          a[k] = fcm_get(fcm[k], i, j);     // in a and compute the statistic
        mat_set(mos, i, j, (*f.u.f1)(a, an));
      }                                     // for the current pair (i,j)
    }
    free(mem); }
    #endif

  else if (f.kind == 2) {                   // --- comp. stats using f2
    assert(arg);
    int N = fcm_dim(*fcm);
    for (int i = 0; i < N; i++) {
      for (int j = i+1; j < N; j++) {
        for (int k = 0; k < n; k++)         // collect vals from all matrices
          a[k] = fcm_get(fcm[k], i, j);
        mat_set(mos, i, j, (*f.u.f2)(a, an, arg, mem));
      }
    }
    free(mem); }

  else
    return -1;                              // return 'failure'

  return 0;                                 // return 'ok'
}  // fcm_uni_single()

/*--------------------------------------------------------------------------*/

/* fcm_uni_wrk
 * -----------
 * worker function for fcm_uni_multi
 *
 * parameters
 * p  pointer to WORK
 *
 * returns
 * THREAD_OK
 */
static void* fcm_uni_wrk(void* p)
{
  assert(p);

  WORK *w = p;

  // allocate aligned memory for a temp. array (for FC values)
  void *mem = malloc((size_t)(w->n) *sizeof(REAL) +31);
  if (!mem) {
    DBGMSG("ERROR: malloc failed");
    w->err = -1;                            // set error indicator
    return THREAD_OK; }                     // return a dummy result
  REAL *a = (REAL*)(((uintptr_t)mem +31) & ~(uintptr_t)31);

  int N = fcm_dim(*(w->fcm));               // get number of nodes

  if      ((*(w->f)).kind == 1) {           // --- comp. stats using f1
    while (1) {                             // process two strip parts
      int i;
      for (i = w->s; i < w->e; i++) {       // traverse pairs of nodes
        for (int j = i+1; j < N; j++) {     // (upper triangular matrix)
          for (int k = 0; k < w->n; k++)    // traverse matrices in order to
            a[k] = fcm_get(w->fcm[k], i, j);// populate tmp array
          mat_set(w->mos, i, j, ((*(w->f)).u.f1)(a, w->an));
        }                                   // comp. stat. for current pair
      }
      if (w->s > N/2) break;                // if second strip done, abort
      i    = N -w->e;                       // get start of opposite stripe
      if (w->e > i)   break;                // if no opposite strip, abort
      w->e = N -w->s;                       // get start and end index
      w->s = i;                             // of the opposite strip
    } }

  else if ((*(w->f)).kind == 2) {           // --- comp. stats using f2
    while (1) {                             // process two strip parts
      int i;
      for (i = w->s; i < w->e; i++) {       // traverse pairs of nodes
        for (int j = i+1; j < N; j++) {     // (upper triangular matrix)
          for (int k = 0; k < w->n; k++)    // traverse matrices in order to
            a[k] = fcm_get(w->fcm[k], i, j);// populate tmp array
          mat_set(w->mos, i, j, ((*(w->f)).u.f2)(a, w->an, w->arg, mem));
        }                                   // comp. stat. for current pair
      }
      if (w->s > N/2) break;                // if second strip done, abort
      i    = N -w->e;                       // get start of opposite stripe
      if (w->e > i)   break;                // if no opposite strip, abort
      w->e = N -w->s;                       // get start and end index
      w->s = i;                             // of the opposite strip
    } }

  else
    w->err = -1;                            // set error indicator

  free(mem);

  return THREAD_OK;                         // return a dummy result
}  // fcm_uni_wrk()

/*--------------------------------------------------------------------------*/

/* fcm_uni_multi
 * -------------
 * compute statistics across functinal connectomes (multi-threaded version)
 *
 * parameters
 * fcm   set of FC matrices
 * n     number of FC matrices
 * f     function to be used (see typedefs in edgestats.h)
 * an    ptr to an array of n-vals, will be passed to f.u.f1/f2 as 2nd arg
 * arg   additional arg. that will be passed to f.u.f2 if relevant
 * mos   result: matrix of statistics
 * nthd  number of threads
 *
 * returns
 * 0 on success
 */
static int fcm_uni_multi(FCMAT **fcm, int n, Func f, const int *an,
                         const REAL *arg, MATRIX *mos, int nthd)
{
  assert(fcm && (n > 0) && an && mos && (nthd > 0));

  int N = fcm_dim(*fcm);                    // number of nodes

  // thread handles, data and worker
  int k = (N/2 +nthd-1) /nthd;              // compute the number of series
  if (k <= 0) k = 1;                        // to be processed per thread

  THREAD *threads = malloc((size_t)nthd *sizeof(THREAD));
  if (!threads) {
    DBGMSG("ERROR: malloc failed");
    return -1; }                            // return 'failure'

  WORK *w = malloc((size_t)nthd *sizeof(WORK));
  if (!w) {
    DBGMSG("ERROR: malloc failed");
    free(threads);
    return -1; }                            // return 'failure'

  WORKER *worker = fcm_uni_wrk;

  // execute the threads
  int i;
  for (i = 0; i < nthd; i++) {              // traverse the threads
    w[i].fcm = fcm;
    w[i].n   = n;
    w[i].f   = &f;
    w[i].an  = an;
    w[i].arg = arg;
    w[i].mos = mos;
    w[i].s   = i*k;                         // compute and store start index
    if (w[i].s >= N/2) break;               // if beyond half, already done
    w[i].e   = w[i].s +k;                   // compute and store end index
    if (w[i].e >= N/2) w[i].e = N -w[i].s;
    w[i].err = 0;                           // init error indicator
    if (pthread_create(threads+i, NULL, worker, w+i)) {
      DBGMSG("ERROR: could not create thread");
      w[i].err = -1;                        // create a thread for each strip
      break; }                              // to comp. the strips in parallel
  }
  int r = 0;                                // error status
  while (--i >= 0) {                        // wait for threads to finish
    pthread_join(threads[i], NULL);         // join threads with this one
    r |= w[i].err;                          // join the error indicators
  }

  free(threads);
  free(w);

  return r;                                 // return error status
}  // fcm_uni_multi()

/*--------------------------------------------------------------------------*/

/* fcm_uni
 * -------
 * compute statistics across functinal connectomes
 *
 * mandatory parameters
 * fcm   set of FC matrices
 * n     number of FC matrices
 * f     function to be used (see above typedefs)
 * an    ptr to an array of n-vals, will be passed to f.u.f1/f2 as 2nd arg
 * arg   additional arg. that will be passed to f.u.f2 if relevant
 * mos   result: matrix of statistics
 * mode  contains bit flags
 *          0           use defaults (no optional parameters)
 *          FCM_THREAD  set number of threads (optional parameter #1)
 *
 * optional parameters
 * #1    number of threads
 *         -1          auto-determine
 *         0          use single-threaded version
 *         1-p        use multi-threaded version with p threads
 * returns
 * 0 on success
 */
int fcm_uni(FCMAT **fcm, int n, Func f, const int *an,
            const REAL *arg, MATRIX *mos, int mode, ...)
{
  assert(fcm && (n > 0) && an && mos);

  int N = fcm[0]->V;                        // number of nodes
  int P = -1;                               // number of threads
  int C = fcm[0]->tile;                     // cache/tile size
  DBGMSG("N: %d  P: %d  C: %d\n", N, P, C);
  assert(N > 0);

  // get optional input
  if (mode & FCM_THREAD) {
   va_list args;
   va_start(args, mode);
   P = va_arg(args, int);
   va_end(args);
  }
  DBGMSG("P: %d\n", P);
  assert((P == -1) || (P >= 0));

  // auto-determine number of threads
  if (P == -1) {
    if (C == 0 || C == N) {
      int nprocs = proccnt();
      P = (nprocs > 1) ? nprocs : 0; }
    else
      P = 0;
  }
  DBGMSG("P: %d\n", P);
  assert(P >= 0);

  // compute statistics
  if      (P >  0)
    return fcm_uni_multi (fcm, n, f, an, arg, mos, P);
  else if (P == 0)
    return fcm_uni_single(fcm, n, f, an, arg, mos);
  else
    return -1;
}  // fcm_uni()
