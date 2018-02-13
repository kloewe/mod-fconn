/*----------------------------------------------------------------------------
  File    : edgestats.h
  Contents: edge-level statistics
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
#ifndef EDGESTATS_H
#define EDGESTATS_H

#include "fcmat.h"
#include "matrix.h"

/*----------------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------------*/
typedef REAL Func1 (const REAL* a, const int *n);
typedef REAL Func2 (const REAL* a, const int *n, const REAL *arg, REAL *buf);

typedef struct {
  int kind;
  union {
    Func1 *f1;
    Func2 *f2;
  } u;
} Func;

/*----------------------------------------------------------------------------
  Auxiliary Functions
----------------------------------------------------------------------------*/
extern REAL corrv_w (const REAL* a, const int *n, const REAL *arg, REAL *buf);

/*----------------------------------------------------------------------------
  Functions
----------------------------------------------------------------------------*/

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
extern int fcm_uni (FCMAT **fcm, int n, Func f, const int *an,
                    const REAL *arg, MATRIX *mos, int mode, ...);

#endif  /* #ifndef EDGESTATS_H */
