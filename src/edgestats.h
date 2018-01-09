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
typedef REAL STATFUNC (const REAL* array, int len);

/*----------------------------------------------------------------------------
  Functions
----------------------------------------------------------------------------*/

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
 *       0          -> use defaults (no optional parameters)
 *       FCM_THREAD -> set number of threads (optional parameter 'nthd')
 *
 * optional parameters
 * #1    number of threads
 *         -1          auto-determine
 *          0          use single-threaded version
 *          1-p        use multi-threaded version with p threads
 * returns
 * 0 on success
 */
extern int fcm_uni (FCMAT **fcm, int n, STATFUNC *func, MATRIX *mos,
                    int mode, ...);

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
 *       0          -> use defaults (no optional parameters)
 *       FCM_THREAD -> set number of threads (optional parameter 'nthd')
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
extern int fcm_corr (FCMAT **fcm, int n, REAL *v, MATRIX *mos, int mode, ...);

/* fcm_tstat2
 * ----------
 * compute t statistics across functional connectomes
 *
 * mandatory parameters
 * fcm   data: set of matrices (representing two independent samples)
 * n     number of matrices
 * g     binary vector of length n indicating sample membership:
 *       0 -> sample #1;  1 -> sample #2
 * mos   result: matrix of t statistics
 * mode  contains bit flags
 *       0          -> use defaults (no optional parameters)
 *       FCM_THREAD -> set number of threads (optional parameter 'nthd')
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
extern int fcm_tstat2(FCMAT **fcm, int n, int *g, MATRIX *mos, int mode, ...);

#endif  /* #ifndef EDGESTATS_H */
