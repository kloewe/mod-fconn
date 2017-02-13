/*----------------------------------------------------------------------------
  File    : nodedeg.c
  Contents: compute the degree of each node in a func. connectivity graph
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdarg.h>
#include <pthread.h>
#include <assert.h>
#include "cpuinfo.h"
#include "fcmat.h"
#include "nodedeg.h"

#ifndef NDEBUG
#line __LINE__ "nodedeg.c"
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
  FCMAT  *fcm;                            // FC matrix
  DIM    s, e;                            // index of start and end series
  REAL   thr;                             // FC threshold
  DIM    *res;                            // result: node degrees
} WORK;

typedef void* WORKER (void*);

/*----------------------------------------------------------------------------
  Functions
----------------------------------------------------------------------------*/

/* fcm_nodedeg_single
 * ---------------------
 * compute node degrees based on a FC matrix and a FC threshold
 *   (single-threaded version)
 *
 * parameters
 * fcm  FC matrix
 * thr  FC threshold
 * res  result: node degrees
 *
 * returns
 * 0 on success
 */
int fcm_nodedeg_single(FCMAT *fcm, REAL thr, DIM *res)
{
  assert(fcm && res);

  int N = fcm_dim(fcm);                   // get number of nodes

  for (int i = 0; i < N; i++)             // traverse nodes
    res[i] = 0;                           // initialize values

  #if 1
  int t = fcm_first(fcm);
  if (t < 0) {                            // check status
    DBGMSG("ERROR: failure when calling fcm_first");
    return -1;                            // return 'failure'
  }
  for ( ; t == 0; t = fcm_next(fcm)) {    // traverse matrix elements
    if (fcm_value(fcm) > thr) {           // increment the appropriate
      res[fcm_row(fcm)]++;                // degree values, if the obtained
      res[fcm_col(fcm)]++;                // FC estimate exceeds the threshold
    }
  }
  if (t < 1) {                            // check status
    DBGMSG("ERROR: failure when calling fcm_next");
    return -1;                            // return 'failure'
  }
  #else
  for (int i = 0; i < N; i++) {           // traverse row indices
    for (int j = i+1; j < N; j++) {       // traverse column indices
      if (fcm_get(fcm,i,j) > thr) {       // increment the appropriate
        res[i]++;                         // degree values, if the obtained
        res[j]++;                         // FC estimate exceeds the threshold
      }
    }
  }
  #endif

  return 0;                               // return 'ok'
}  // fcm_nodedeg_single()

/*--------------------------------------------------------------------------*/

/* fcm_nodedeg_wrk
 * ---------------
 * worker function for fmc_nodedeg_multi
 *
 * parameters
 * p     pointer to the data
 *
 * returns
 * THREAD_OK
 */
static void* fcm_nodedeg_wrk(void* p)
{
  assert(p);

  WORK *w = p;

  int n = fcm_dim(w->fcm);                // get number of nodes

  for (int i = 0; i < n; i++)             // traverse nodes
    w->res[i] = 0;                        // initialize values

  while (1) {                             // process two strip parts
    int i;
    for (i = w->s; i < w->e; i++) {       // traverse row indices
      for (int j = i+1; j < n; j++) {     // traverse column indices
        if (fcm_get(w->fcm,i,j) > w->thr) {
          w->res[i]++;                    // increment degree values,
          w->res[j]++;                    // if the obtained FC estimate
        }                                 // exceeds the threshold
      }
    }
    if (w->s > n/2) break;                // if second strip done, abort
    i    = n -w->e;                       // get start of opposite stripe
    if (w->e > i)   break;                // if no opposite strip, abort
    w->e = n -w->s;                       // get start and end index
    w->s = i;                             // of the opposite strip
  }

  return THREAD_OK;                       // return a dummy result
}  // fcm_nodedeg_wrk()

/*--------------------------------------------------------------------------*/

/* fcm_nodedeg_multi
 * -----------------
 * compute node degrees based on a FC matrix and a FC threshold
 *   (multi-threaded version)
 *
 * parameters
 * fcm   FC matrix
 * thr   FC threshold
 * res   result: node degrees
 * nthd  number of threads
 *
 * returns
 * 0 on success
 */
int fcm_nodedeg_multi(FCMAT *fcm, REAL thr, DIM *res, int nthd)
{
  assert(fcm && res && (nthd > 0));

  int n = fcm_dim(fcm);                   // number of nodes

  // thread, handles, data and worker
  int k = (n/2 +nthd-1) /nthd;            // compute the number of series
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
  WORKER *worker = fcm_nodedeg_wrk;

  // execute the threads
  int i;                                  // loop variable
  int r = 0;                              // error status
  for (i = 0; i < nthd; i++) {            // traverse the threads
    w[i].fcm = fcm;                       // FC matrix
    w[i].thr = thr;                       // FC threshold
    w[i].res = malloc((size_t)n*sizeof(DIM)); // partial result
    w[i].s   = (DIM)i*k;                  // compute and store start index
    if (w[i].s >= n/2) break;             // if beyond half, already done
    w[i].e    = w[i].s +k;                // compute and store end index
    if (w[i].e >= n/2) w[i].e = n -w[i].s;
    if (pthread_create(threads+i, NULL, worker, w+i)) {
      DBGMSG("ERROR: could not create thread");
      r = -1;                             // create a thread for each strip
      break; }                            // to compute the strips in parallel
  }
  for (int j = 0; j < n; j++)             // result: traverse nodes
    res[j] = 0;                           // initialize values
  while (--i >= 0) {                      // wait for threads to finish
    pthread_join(threads[i], NULL);       // (join threads with this one)
    for (int j = 0; j < n; j++)           // combine partial results
      res[j] += w[i].res[j];              // from the individual threads
    free(w[i].res);                       // deallocate partial results
  }

  free(threads);
  free(w);

  return r;                               // return error status
}  // fcm_nodedeg_multi()

/*--------------------------------------------------------------------------*/

/* fcm_nodedeg
 * -----------
 * compute node degrees based on a FC matrix and a FC threshold
 *
 * mandatory parameters
 * fcm   FC matrix
 * thr   FC threshold
 * res   result: node degrees
 * mode  contains bit flags
 *         0           use defaults (no optional parameters)
 *         FCM_THREAD  set number of threads (optional parameter #1)
 *
 * optional parameters
 * #1    number of threads
 *        -1    auto-determine
 *         0    use single-threaded version
 *         1-p  use multi-threaded version with p threads
 *
 * returns
 * 0 on success
 */
int fcm_nodedeg(FCMAT *fcm, REAL thr, DIM *res, int mode, ...)
{
  assert(fcm && res);

  int N = fcm->V;                         // number of nodes
  int P = -1;                             // number of threads
  int C = fcm->tile;                      // cache/tile size
  DBGMSG("N: %d  P: %d  C: %d\n", N, P, C);

  // get optional input
  if (mode & FCM_THREAD) {
    va_list args;
    va_start(args, mode);
    P = va_arg(args, int);
    va_end(args);
  }
  DBGMSG("N: %d  P: %d  C: %d\n", N, P, C);

  // auto-determine # threads
  if (P == -1) {
    if (C == 0 || C == N) {
      int nprocs = proccnt();
      P = (nprocs > 1) ? nprocs : 0; }
    else
      P = 0;
  }
  DBGMSG("N: %d  P: %d  C: %d\n", N, P, C);

  // compute degrees
  if      (P >  0)
    return fcm_nodedeg_multi (fcm, thr, res, P);
  else if (P == 0)
    return fcm_nodedeg_single(fcm, thr, res);
  else
    return -1;
}  // fcm_nodedeg()
