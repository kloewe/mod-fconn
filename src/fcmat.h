/*----------------------------------------------------------------------------
  File    : fcmat.h
  Contents: data type for functional connectivity matrix
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
#ifndef __FCMAT__
#define __FCMAT__

/*----------------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------------*/
#ifndef REAL
#define REAL float
#endif
#ifndef DIM
#define DIM size_t
#endif

#ifndef FCM_CORR
#define FCM_CORR  0x0f /* mask for correlation estimation method */
#define FCM_PCC      1 /* Pearson correlation coefficient  */
#define FCM_TCC      2 /* tetrachoric correlation coefficient */

#define FCM_R2Z   0x10 /* Fisher r-to-z transform */
#endif

/*----------------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------------*/
typedef struct {
  DIM        V;        /* number of voxels */
  DIM        T;        /* number of scans */
  REAL       *data;    /* voxel time series data */
  int        flags;
  int        X;        /* if pcc: size of padded data arrays
                        * if tcc: binarized array size */
  void       *buf[2];  /* for precomp. subexpressions etc. */
  REAL       last[3];  /* cache for last element (i,j,val) */
} FCMAT;

/*----------------------------------------------------------------------------
  Functions
----------------------------------------------------------------------------*/
extern FCMAT* fcm_create (REAL *data, DIM V, DIM T, int flags);
extern void   fcm_delete (FCMAT *fcm);
extern REAL   fcm_get    (FCMAT *fcm, DIM i, DIM j);
extern void   fcm_print  (FCMAT *fcm);

/*----------------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------------*/
#define fcm_dim(m)        ((m)->V)

#endif
