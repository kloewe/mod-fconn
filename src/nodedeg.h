/*----------------------------------------------------------------------------
  File    : nodedeg.h
  Contents: compute the degree of each node in a func. connectivity graph
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
#ifndef NODEDEG_H
#define NODEDEG_H

#include "fcmat.h"

/*----------------------------------------------------------------------------
  Functions
----------------------------------------------------------------------------*/

/* fcm_nodedeg
 * -----------------
 * compute node degrees based on a FC matrix and a FC threshold
 *
 * mandatory parameters
 * fcm   FC matrix
 * thr   FC threshold
 * res   result: node degrees
 * mode  contains bit flags
 *       0          -> use defaults (no optional parameters)
 *       FCM_THREAD -> set number of threads (optional parameter 'nthd')
 *
 * optional parameters
 * nthd  number of threads
 * 
 * returns
 * 0 on success
 */
extern int fcm_nodedeg(FCMAT *fcm, REAL thr, DIM *res, int mode, ...);


#endif  /* #ifndef NODEDEG_H */
