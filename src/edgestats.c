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

#if defined(MATLAB_MEX_FILE)
#  include "mex.h"
#endif

#ifndef NDEBUG
#  line __LINE__ "edgestats.c"
#endif

#include "real-is-double.inc"
#if REAL_IS_DOUBLE
#  define dot ddot
#else
#  define dot sdot
#endif

/*----------------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------------*/
#define THREAD    pthread_t               // use the POSIX thread type
#define THREAD_OK NULL                    // return value is void*

/*----------------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------------*/
typedef struct {                          // --- thread worker data ---
  FCMAT    **fcm;                         // sample: set of FC matrices
  int      n;                             // sample size
  MATRIX   *mos;                          // result: matrix of statistics
  void     *var;                          // additional variable
  STATFUNC *func;                         // function pointer
  int      s, e;                          // index of start and end series
  int      err;                           // error indicator
} WORK;

typedef void* WORKER (void*);

typedef void (*fn_ptr)(void);             // function pointer

/*----------------------------------------------------------------------------
  Functions
----------------------------------------------------------------------------*/

/* fcm_uni_single
 * --------------
 * comp. descriptive statistics for univariate data across func. connectomes
 *   (single-threaded version)
 *
 * parameters
 * fcm   data: set of matrices
 * n     number of matrices
 * f     function pointer
 * mos   result: matrix of statistics
 *
 * returns
 * 0 on success
 */
static int fcm_uni_single(FCMAT **fcm, int n, fn_ptr f, MATRIX *mos)
{
  assert(fcm && f && mos && (n > 0));

  STATFUNC *func = (STATFUNC*)f;          // statistical function

  // allocate aligned memory for a temp. array (for FC values)
  void *mem = malloc((size_t)n *sizeof(REAL) +31);
  if (!mem) {
    DBGMSG("ERROR: malloc failed");
    return -1; }                          // return 'failure'
  REAL *a = (REAL*)(((uintptr_t)mem +31) & ~(uintptr_t)31);

  // compute statistics
  #if 1                                   // variant 1: use iterator
  int t = 0;                              // init status variable t
  for (int k = 0; k < n; k++)             // initiate iterations
    t += fcm_first(fcm[k]);
  if (t < 0) {                            // check status
    DBGMSG("ERROR: failure when calling fcm_first");
    return -1;                            // return 'failure'
  }
  while (t == 0) {                        // while iterations not completed
    for (int k = 0; k < n; k++)           // collect values from all matrices
      a[k] = fcm_value(fcm[k]);           // in a and compute the statistic
    mat_set(mos, fcm_row(fcm[0]), fcm_col(fcm[0]), (*func)(a, n));
    for (int k = 0; k < n; k++)           // for the current pair, then
      t += fcm_next(fcm[k]);              // move on to the next elements
  }
  free(mem);
  if (t < n) {                            // check status
    DBGMSG("ERROR: failure when calling fcm_next");
    return -1;                            // return 'failure'
  }
  #else                                   // variant 2: use fcm_get
  int N = fcm_dim(*fcm);                  // determine number of nodes
  for (int i = 0; i < N; i++) {           // traverse pairs of nodes 
    for (int j = i+1; j < N; j++) {       // (upper triangular matrix)
      for (int k = 0; k < n; k++)         // collect values from all matrices
        a[k] = fcm_get(fcm[k], i, j);     // in a and compute the statistic
      mat_set(mos, i, j, (*func)(a, n));  // for the current pair (i,j)
    }
  }
  free(mem);
  #endif

  return 0;                               // return 'ok'
}  // fcm_uni_single()

/*--------------------------------------------------------------------------*/

/* fcm_uni_wrk
 * -----------
 * worker function for fcm_uni_multi
 *
 * parameters
 * p     pointer to the data
 *
 * returns
 * THREAD_OK
 */
static void* fcm_uni_wrk(void* p)
{
  assert(p);

  WORK *w = p;

  // allocate aligned memory for a temp. array (for FC values)
  void *mem = malloc((size_t)w->n *sizeof(REAL) +31);
  if (!mem) {
    DBGMSG("ERROR: malloc failed");
    w->err = -1;                          // set error indicator
    return THREAD_OK; }                   // return a dummy result
  REAL *a = (REAL*)(((uintptr_t)mem +31) & ~(uintptr_t)31);

  // compute statistics
  int N = fcm_dim(*(w->fcm));             // determine number of nodes
  while (1) {                             // process two strip parts
    int i;
    for (i = w->s; i < w->e; i++) {       // traverse pairs of nodes
      for (int j = i+1; j < N; j++) {     // (upper triangular matrix)
        for (int k = 0; k < w->n; k++)    // traverse matrices in order to
          a[k] = fcm_get(w->fcm[k], i, j);// populate tmp array
        mat_set(w->mos, i, j, (*(w->func))(a, w->n));
      }                                   // comp. statistic for current pair
    }
    if (w->s > N/2) break;                // if second strip done, abort
    i    = N -w->e;                       // get start of opposite stripe
    if (w->e > i)   break;                // if no opposite strip, abort
    w->e = N -w->s;                       // get start and end index
    w->s = i;                             // of the opposite strip
  }

  free(mem);

  return THREAD_OK;                       // return a dummy result
}  // fcm_uni_wrk()

/*--------------------------------------------------------------------------*/

/* fcm_uni_multi
 * -------------
 * comp. descriptive statistics for univariate data across func. connectomes
 *   (multi-threaded version)
 *
 * parameters
 * fcm   data: set of matrices
 * n     number of matrices
 * f     function pointer
 * mos   result: matrix of statistics
 * nthd  number of threads
 *
 * returns
 * 0 on success
 */
static int fcm_uni_multi(FCMAT **fcm, int n, fn_ptr f, MATRIX *mos, int nthd)
{
  assert(fcm && f && mos && (nthd > 0));

  int N = fcm_dim(*fcm);                  // number of nodes
  STATFUNC *func = (STATFUNC*)f;          // statistical function

  // thread handles, data and worker
  int k = (N/2 +nthd-1) /nthd;            // compute the number of series
  if (k <= 0) k = 1;                      // to be processed per thread
  THREAD *threads = malloc((size_t)nthd *sizeof(THREAD));
  if (!threads) {
    DBGMSG("ERROR: malloc failed");
    return -1; }                          // return 'failure'
  WORK *w = malloc((size_t)nthd *sizeof(WORK));
  if (!w) {
    DBGMSG("ERROR: malloc failed");
    free(threads);
    return -1; }                          // return 'failure'
  WORKER *worker = fcm_uni_wrk;

  // execute the threads
  int i;
  for (i = 0; i < nthd; i++) {            // traverse the threads
    w[i].fcm  = fcm;
    w[i].n    = n;
    w[i].mos  = mos;
    w[i].var  = NULL;
    w[i].func = (*func);
    w[i].s    = i*k;                      // compute and store start index
    if (w[i].s >= N/2) break;             // if beyond half, already done
    w[i].e    = w[i].s +k;                // compute and store end index
    if (w[i].e >= N/2) w[i].e = N -w[i].s;
    w[i].err  = 0;                        // init error indicator
    if (pthread_create(threads+i, NULL, worker, w+i)) {
      DBGMSG("ERROR: could not create thread");
      w[i].err = -1;                      // create a thread for each strip
      break; }                            // to comp. the strips in parallel
  }
  int r = 0;                              // error status
  while (--i >= 0) {                      // wait for threads to finish
    pthread_join(threads[i], NULL);       // join threads with this one
    r |= w[i].err;                        // join the error indicators
  }

  free(threads);
  free(w);

  return r;                               // return error status
}  // fcm_uni_multi()

/*--------------------------------------------------------------------------*/

/* fcm_uni
 * -------
 * comp. descriptive statistics for univariate data across func. connectomes
 *
 * mandatory parameters
 * fcm   data: set of matrices
 * n     number of matrices
 * func  function pointer
 * mos   result: matrix of statistics
 * mode  contains bit flags
 *         0           use defaults (no optional parameters)
 *         FCM_THREAD  set number of threads (optional parameter #1)
 *
 * optional parameters
 * #1    number of threads
 *         -1          auto-determine
 *          0          use single-threaded version
 *          1-p        use multi-threaded version with p threads
 *
 * returns
 * 0 on success
 */
int fcm_uni(FCMAT **fcm, int n, STATFUNC *func, MATRIX *mos, int mode, ...)
{
  assert(fcm && func && mos && (n > 0));

  int N = fcm[0]->V;                      // number of nodes
  int P = -1;                             // number of threads
  int C = fcm[0]->tile;                   // cache/tile size
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
    return fcm_uni_multi (fcm, n, (fn_ptr)func, mos, P);
  else if (P == 0)
    return fcm_uni_single(fcm, n, (fn_ptr)func, mos);
  else
    return -1;
}  // fcm_uni()

/*--------------------------------------------------------------------------*/

/* fcm_corr_single
 * ---------------
 * compute correlation coefficients across functional connectomes
 *   (single-threaded version)
 *
 * parameters
 * fcm   data: set of matrices
 * n     number of matrices
 * v     additional variable
 * mos   result: matrix of correlation coefficients
 *
 * returns
 * 0 on success
 */
static int fcm_corr_single(FCMAT **fcm, int n, REAL *v, MATRIX *mos)
{
  assert(fcm && v && mos && (n > 0));

  // allocate (aligned) memory for the pre-normalized values (add. variable)
  void *mem1 = malloc((size_t)n *sizeof(REAL) +31);
  if (!mem1) {
    DBGMSG("ERROR: malloc failed");
    return -1; }                          // return 'failure'
  REAL *a = (REAL*)(((uintptr_t)mem1 +31) & ~(uintptr_t)31);

  // allocate (aligned) memory for the pre-normalized values (FC values)
  void *mem2 = malloc((size_t)n *sizeof(REAL) +31);
  if (!mem2) {
    DBGMSG("ERROR: malloc failed");
    free(mem1);
    return -1; }                          // return 'failure'
  REAL *b = (REAL*)(((uintptr_t)mem2 +31) & ~(uintptr_t)31);

  // pre-normalize the add. variable
  REAL sqr = 0;
  REAL ma = mean(v, n);
  for (int k = 0; k < n; k++) {
    a[k] = v[k] - ma;
    sqr += a[k]*a[k]; }
  sqr = sqrt(sqr);
  sqr = (sqr > 0) ? 1/sqr : 0;
  for (int k = 0; k < n; k++)
    a[k] *= sqr;

  // compute correlation coefficients for all pairs (i,j)
  int N = fcm_dim(*fcm);
  for (int i = 0; i < N; i++) {
    for (int j = i+1; j < N; j++) {

      // pre-normalize the FC values
      sqr = 0;
      for (int k = 0; k < n; k++)
        b[k] = fcm_get(fcm[k], i, j);
      REAL mb = mean(b, n);
      for (int k = 0; k < n; k++) {
        b[k] -= mb;
        sqr += b[k]*b[k]; }
      sqr = sqrt(sqr);
      sqr = (sqr > 0) ? 1/sqr : 0;
      for (int k = 0; k < n; k++)
        b[k] *= sqr;

      mat_set(mos, i, j, dot(a, b, n));
    }
  }

  free(mem1);
  free(mem2);

  return 0;                               // return 'ok'
}  // fcm_corr_single()

/*--------------------------------------------------------------------------*/

/* fcm_corr_wrk
 * ------------
 * worker function for fcm_corr_multi
 *
 * parameters
 * p     pointer to the data
 *
 * returns
 * THREAD_OK
 */
static void* fcm_corr_wrk(void* p)
{
  assert(p);

  WORK *w = p;

  // allocate (aligned) memory for the pre-normalized values (add. variable)
  void *mem1 = malloc((size_t)(w->n) *sizeof(REAL) +31);
  if (!mem1) {
    DBGMSG("ERROR: malloc failed");
    w->err = -1;                          // set error indicator
    return THREAD_OK; }                   // return a dummy result
  REAL *a = (REAL*)(((uintptr_t)mem1 +31) & ~(uintptr_t)31);

  // allocate (aligned) memory for the pre-normalized values (FC values)
  void *mem2 = malloc((size_t)(w->n) *sizeof(REAL) +31);
  if (!mem2) {
    DBGMSG("ERROR: malloc failed");
    w->err = -1;                          // set error indicator
    free(mem1);
    return THREAD_OK; }                   // return a dummy result
  REAL *b = (REAL*)(((uintptr_t)mem2 +31) & ~(uintptr_t)31);

  // pre-normalize the add. variable
  REAL sqr = 0;
  REAL ma = mean(w->var, w->n);
  for (int k = 0; k < w->n; k++) {
    a[k] = *((REAL *)(w->var)+k) - ma;
    sqr += a[k]*a[k]; }
  sqr = sqrt(sqr);
  sqr = (sqr > 0) ? 1/sqr : 0;
  for (int k = 0; k < w->n; k++)
    a[k] *= sqr;

  // compute correlation coefficients
  int N = fcm_dim(*(w->fcm));             // determine number of nodes
  while (1) {                             // process two strip parts
    int i;
    for (i = w->s; i < w->e; i++) {       // traverse row indices
      for (int j = i+1; j < N; j++) {     // traverse column indices

        // pre-normalize the FC values
        sqr = 0;
        for (int k = 0; k < w->n; k++)
          b[k] = fcm_get(w->fcm[k], i, j);
        REAL mb = mean(b, w->n);
        for (int k = 0; k < w->n; k++) {
          b[k] -= mb;
          sqr += b[k]*b[k]; }
        sqr = sqrt(sqr);
        sqr = (sqr > 0) ? 1/sqr : 0;
        for (int k = 0; k < w->n; k++)
          b[k] *= sqr;

        mat_set(w->mos, i, j, dot(a, b, w->n));
      }
    }
    if (w->s > N/2) break;                // if second strip done, abort
    i    = N -w->e;                       // get start of opposite stripe
    if (w->e > i)   break;                // if no opposite strip, abort
    w->e = N -w->s;                       // get start and end index
    w->s = i;                             // of the opposite strip
  }

  free(mem1);
  free(mem2);

  return THREAD_OK;                       // return a dummy result
}  // fcm_corr_wrk()

/*--------------------------------------------------------------------------*/

/* fcm_corr_multi
 * --------------
 * compute correlation coefficients across functional connectomes
 *   (multi-threaded version)
 *
 * parameters
 * fcm   data: set of matrices
 * n     number of matrices
 * v     additional variable
 * mos   result: matrix of correlation coefficients
 * nthd  number of threads
 *
 * returns
 * 0 on success
 */
static int fcm_corr_multi(FCMAT **fcm, int n, REAL *v, MATRIX *mos, int nthd)
{
  assert(fcm && v && mos && (nthd > 0));

  int N = fcm_dim(*fcm);                  // number of nodes

  // thread handles, data and worker
  int k = (N/2 +nthd-1) /nthd;            // compute the number of series
  if (k <= 0) k = 1;                      // to be processed per thread
  THREAD *threads = malloc((size_t)nthd *sizeof(THREAD));
  if (!threads) {
    DBGMSG("ERROR: malloc failed");
    return -1; }                          // return 'failure'
  WORK *w = malloc((size_t)nthd *sizeof(WORK));
  if (!w) {
    DBGMSG("ERROR: malloc failed");
    free(threads);
    return -1; }                          // return 'failure'
  WORKER *worker = fcm_corr_wrk;

  // execute the threads
  int i;
  for (i = 0; i < nthd; i++) {            // traverse the threads
    w[i].fcm  = fcm;
    w[i].n    = n;
    w[i].mos  = mos;
    w[i].var  = v;
    w[i].func = NULL;
    w[i].s    = i*k;                      // compute and store start index
    if (w[i].s >= N/2) break;             // if beyond half, already done
    w[i].e    = w[i].s +k;                // compute and store end index
    if (w[i].e >= N/2) w[i].e = N -w[i].s;
    w[i].err  = 0;                        // init error indicator
    if (pthread_create(threads+i, NULL, worker, w+i)) {
      DBGMSG("ERROR: could not create thread");
      w[i].err = -1;                      // create a thread for each strip
      break; }                            // to comp. the strips in parallel
  }
  int r = 0;                              // error status
  while (--i >= 0) {                      // wait for threads to finish
    pthread_join(threads[i], NULL);       // join threads with this one
    r |= w[i].err;                        // join the error indicators
  }

  free(threads);
  free(w);

  return r;                               // return error status
}  // fcm_corr_multi()

/*--------------------------------------------------------------------------*/

/* fcm_corr
 * --------
 * compute correlation coefficients across functional connectomes
 *
 * mandatory parameters
 * fcm   data: set of matrices
 * n     number of matrices
 * v     additional variable
 * mos   result: matrix of correlation coefficients
 * mode  contains bit flags
 *         0           use defaults (no optional parameters)
 *         FCM_THREAD  set number of threads (optional parameter #1)
 *
 * optional parameters
 * #1    number of threads
 *         -1          auto-determine
 *          0          use single-threaded version
 *          1-p        use multi-threaded version with p threads
 *
 * returns
 * 0 on success
 */
int fcm_corr(FCMAT **fcm, int n, REAL *v, MATRIX *mos, int mode, ...)
{
  assert(fcm && v && mos && (n > 0));

  int N = fcm[0]->V;                      // number of nodes
  int P = -1;                             // number of threads
  int C = fcm[0]->tile;                   // cache/tile size
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

  // compute correlation coefficients
  if      (P >  0)
    return fcm_corr_multi (fcm, n, v, mos, P);
  else if (P == 0)
    return fcm_corr_single(fcm, n, v, mos);
  else
    return -1;
}  // fcm_corr()

/*--------------------------------------------------------------------------*/

/* fcm_tstat2_single
 * -----------------
 * compute t statistics across functional connectomes
 *   (single-threaded version)
 *
 * parameters
 * fcm   data: set of matrices (representing two independent samples)
 * n     number of matrices
 * g     binary vector of length n indicating sample membership:
 *       0 -> sample #1;  1 -> sample #2
 * mos   result: matrix of t statistics
 *
 * returns
 * 0 on success
 */
static int fcm_tstat2_single(FCMAT **fcm, int n, int *g, MATRIX *mos)
{
  assert(fcm && g && mos && (n > 0));

  // determine sample sizes
  int n1 = 0, n2 = 0;
  for (int k = 0; k < n; k++) {
    if (g[k] == 0) n1++;
    else           n2++;
  }
  assert((n1 + n2) == n);

  // split fcm into fcm1 and fcm2 based on sample membership
  FCMAT **fcm1 = malloc((size_t)n *sizeof(FCMAT*));
  if (!fcm1) return -1;
  FCMAT **fcm2 = fcm1 + n1;
  for (int k = 0; k < n; k++) {
    if (g[k] == 0)  *(fcm1++) = fcm[k];
    else            *(fcm2++) = fcm[k];
  }
  fcm1 -= n1; fcm2 -= n2;                 // reset pointers

  // allocate aligned memory for two temp. arrays (for FC values)
  void *mem1 = malloc((size_t)n1 *sizeof(REAL) +31);
  if (!mem1) {
    DBGMSG("ERROR: malloc failed");
    free(fcm1);
    return -1; }                          // return 'failure'
  REAL *tmp1 = (REAL*)(((uintptr_t)mem1 +31) & ~(uintptr_t)31);

  void *mem2 = malloc((size_t)n2 *sizeof(REAL) +31);
  if (!mem2) {
    DBGMSG("ERROR: malloc failed");
    free(mem1);
    free(fcm1);
    return -1; }                          // return 'failure'
  REAL *tmp2 = (REAL*)(((uintptr_t)mem2 +31) & ~(uintptr_t)31);

  // compute t statistics
  int N = fcm_dim(*fcm);                  // determine number of nodes
  for (int i = 0; i < N; i++) {           // traverse pairs of nodes
    for (int j = i+1; j < N; j++) {       // (upper triangular matrix)
      for (int k = 0; k < n1; k++)        // traverse samples in order to
        tmp1[k] = fcm_get(fcm1[k], i, j);
      for (int k = 0; k < n2; k++)        // populate temp. arrays
        tmp2[k] = fcm_get(fcm2[k], i, j);
      mat_set(mos, i, j, tstat2(tmp1, tmp2, n1, n2));
    }                                     // compute t statistic
  }

  free(fcm1);
  free(mem1);
  free(mem2);

  return 0;                               // return 'ok'
}  // fcm_tstat2_single()

/*--------------------------------------------------------------------------*/

/* fcm_tstat2_wrk
 * --------------
 * worker function for fcm_tstat2_multi
 *
 * parameters
 * p     pointer to the data
 *
 * returns
 * THREAD_OK
 */
static void* fcm_tstat2_wrk(void* p)
{
  assert(p);

  WORK *w = p;

  int *g = (int *)(w->var);               // sample membership

  // determine sample sizes
  int n1 = 0, n2 = 0;
  for (int k = 0; k < w->n; k++) {
    if (g[k] == 0) n1++;
    else           n2++;
  }
  assert((n1 + n2) == w->n);

  // split fcm into fcm1 and fcm2 based on sample membership
  FCMAT **fcm1 = malloc((size_t)(w->n) *sizeof(FCMAT*));
  if (!fcm1) {
    DBGMSG("ERROR: malloc failed");
    w->err = -1;                          // set error indicator
    return THREAD_OK; }                   // return a dummy result
  FCMAT **fcm2 = fcm1 + n1;
  for (int k = 0; k < w->n; k++) {
    if (g[k] == 0)  *(fcm1++) = w->fcm[k];
    else            *(fcm2++) = w->fcm[k];
  }
  fcm1 -= n1; fcm2 -= n2;                 // reset pointers

  // allocate aligned memory for two temp. arrays (for FC values)
  void *mem1 = malloc((size_t)n1 *sizeof(REAL) +31);
  if (!mem1) {
    DBGMSG("ERROR: malloc failed");
    w->err = -1;                          // set error indicator
    free(fcm1);
    return THREAD_OK; }                   // return a dummy result
  REAL *tmp1 = (REAL*)(((uintptr_t)mem1 +31) & ~(uintptr_t)31);
  void *mem2 = malloc((size_t)n2 *sizeof(REAL) +31);
  if (!mem2) {
    DBGMSG("ERROR: malloc failed");
    w->err = -1;                          // set error indicator
    free(mem1);
    free(fcm1);
    return THREAD_OK; }                   // return a dummy result
  REAL *tmp2 = (REAL*)(((uintptr_t)mem2 +31) & ~(uintptr_t)31);

  // compute t statistics
  int N = fcm_dim(*(w->fcm));             // determine number of nodes
  while (1) {                             // process two strip parts
    int i;
    for (i = w->s; i < w->e; i++) {       // traverse row indices
      for (int j = i+1; j < N; j++) {     // traverse column indices
        for (int k = 0; k < n1; k++)      // traverse samples in order to
          tmp1[k] = fcm_get(fcm1[k], i, j);
        for (int k = 0; k < n2; k++)      // populate temp. arrays
          tmp2[k] = fcm_get(fcm2[k], i, j);
        mat_set(w->mos, i, j, tstat2(tmp1, tmp2, n1, n2));
      }                                   // compute t statistic
    }
    if (w->s > N/2) break;                // if second strip done, abort
    i    = N -w->e;                       // get start of opposite stripe
    if (w->e > i)   break;                // if no opposite strip, abort
    w->e = N -w->s;                       // get start and end index
    w->s = i;                             // of the opposite strip
  }

  free(fcm1);
  free(mem1);
  free(mem2);

  return THREAD_OK;                       // return a dummy result
}  // fcm_tstat2_wrk()

/*--------------------------------------------------------------------------*/

/* fcm_tstat2_multi
 * ---------------
 * compute t statistics across functional connectomes
 *   (multi-threaded version)
 *
 * parameters
 * fcm   data: set of matrices (representing two independent samples)
 * n     number of matrices
 * g     binary vector of length n indicating sample membership:
 *       0 -> sample #1;  1 -> sample #2
 * mos   result: matrix of t statistics
 * nthd  number of threads
 *
 * returns
 * 0 on success
 */
static int fcm_tstat2_multi(FCMAT **fcm, int n, int *g, MATRIX *mos, int nthd)
{
  assert(fcm && g && mos && (nthd > 0));

  int    N = fcm_dim(*fcm);               // number of nodes

  // thread handles, data and worker
  int k = (N/2 +nthd-1) /nthd;            // compute the number of series
  if (k <= 0) k = 1;                      // to be processed per thread
  THREAD *threads = malloc((size_t)nthd *sizeof(THREAD));
  if (!threads) {
    DBGMSG("ERROR: malloc failed");
    return -1; }                          // return 'failure'
  WORK *w = malloc((size_t)nthd *sizeof(WORK));
  if (!w) {
    DBGMSG("ERROR: malloc failed");
    free(threads);
    return -1; }                          // return 'failure'
  WORKER *worker = fcm_tstat2_wrk;

  // execute the threads
  int i;
  for (i = 0; i < nthd; i++) {            // traverse the threads
    w[i].fcm  = fcm;
    w[i].n    = n;
    w[i].mos  = mos;
    w[i].var  = g;
    w[i].func = NULL;
    w[i].s    = i*k;                      // compute and store start index
    if (w[i].s >= N/2) break;             // if beyond half, already done
    w[i].e    = w[i].s +k;                // compute and store end index
    if (w[i].e >= N/2) w[i].e = N -w[i].s;
    w[i].err  = 0;                        // init error indicator
    if (pthread_create(threads+i, NULL, worker, w+i)) {
      DBGMSG("ERROR: could not create thread");
      w[i].err = -1;                      // create a thread for each strip
      break; }                            // to comp. the strips in parallel
  }
  int r = 0;                              // error status
  while (--i >= 0) {                      // wait for threads to finish
    pthread_join(threads[i], NULL);       // join threads with this one
    r |= w[i].err;                        // join the error indicators
  }

  free(threads);
  free(w);

  return r;                               // return error status
}  // fcm_tstat2_multi()

/*--------------------------------------------------------------------------*/

/* fcm_tstat2
 * ---------
 * compute t statistics across functional connectomes
 *
 * mandatory parameters
 * fcm   data: set of matrices (representing two independent samples)
 * n     number of matrices
 * g     binary vector of length n indicating sample membership:
 *       0 -> sample #1;  1 -> sample #2
 * mos   result: matrix of t statistics
 * mode  contains bit flags
 *         0           use defaults (no optional parameters)
 *         FCM_THREAD  set number of threads (optional parameter #1)
 *
 * optional parameters
 * #1    number of threads
 *         -1          auto-determine
 *          0          use single-threaded version
 *          1-p        use multi-threaded version with p threads
 *
 * returns
 * 0 on success
 */
int fcm_tstat2(FCMAT **fcm, int n, int *g, MATRIX *mos, int mode, ...)
{
  assert(fcm && g && mos && (n > 0));

  int N = fcm[0]->V;                      // number of nodes
  int P = -1;                             // number of threads
  int C = fcm[0]->tile;                   // cache/tile size
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

  // compute t statistics
  if      (P >  0)
    return fcm_tstat2_multi (fcm, n, g, mos, P);
  else if (P == 0)
    return fcm_tstat2_single(fcm, n, g, mos);
  else
    return -1;
}  // fcm_tstat2()
