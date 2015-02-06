/*----------------------------------------------------------------------
  File    : fcmat.h
  Contents: data type for functional connectivity matrix
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
#ifndef __FCMAT__

#ifdef _MSC_VER
#define uint32_t   unsigned __int32
#define uint64_t   unsigned __int64
#define inline     __inline
#else                           /* MSC still does not support C99 */
#include <stdint.h>
#endif

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <pthread.h>
#endif

/*----------------------------------------------------------------------
  Data Type Definition / Recursion Handling
----------------------------------------------------------------------*/
#ifndef DIM                     /* if matrix dimension is not defined */
#define DIM ptrdiff_t           /* use ptrdiff_t as the default */
#endif

#ifdef REAL                     /* if REAL is defined */
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

#ifndef SFXNAME                 /* macros to generate function names */
#define SFXNAME(n)      SFXNAME_1(n,SUFFIX)
#define SFXNAME_1(n,s)  SFXNAME_2(n,s)
#define SFXNAME_2(n,s)  n##s    /* the two step recursion is needed */
#endif                          /* to ensure proper expansion */

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#ifndef FCM_CORR                /* if not yet defined */
#define FCM_CORR    0x000f      /* mask for correlation method */
#define FCM_PCC     0x0001      /* Pearson correlation coefficient  */
#define FCM_TCC     0x0002      /* tetrachoric correlation coeff. */

#define FCM_R2Z     0x0010      /* Fisher r-to-z transform */

#define FCM_CACHE   0x0100      /* use a cache (needs tile size) */
#define FCM_THREAD  0x0200      /* threads (needs thread count) */
#endif

#ifndef THREAD                  /* if not yet defined */
#ifdef _WIN32                   /* if Microsoft Windows system */
#define THREAD      HANDLE      /* threads are identified by handles */
#else                           /* if Linux/Unix system */
#define THREAD      pthread_t   /* use the POSIX thread type */
#endif                          /* definition of a worker function */
#endif

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
struct SFXNAME(fcmat);          /* --- element retrieval function */
typedef REAL SFXNAME(FCMGETFN) (struct SFXNAME(fcmat) *fcmat,
                                DIM row, DIM col);

typedef struct SFXNAME(fcmat) { /* --- a func. connectivity matrix */
  DIM    V;                     /* number of voxels */
  DIM    T;                     /* number of scans */
  DIM    X;                     /* size of padded/binarized data */
  DIM    tile;                  /* size of tiles/blocks for caching */
  int    mode;                  /* processing mode (e.g. FCM_PCC) */
  DIM    nthd;                  /* number of threads for computation */
  void   *data;                 /* normalized/binarized data */
  void   *mem;                  /* allocated memory block */
  REAL   *cmap;                 /* map from n_11 to cosine values */
  REAL   *cache;                /* cached rectangle/triangle */
  REAL   diag;                  /* value of diagonal element */
  DIM    ra, rb;                /* row    range of cached area */
  DIM    ca, cb;                /* column range of cached area */
  DIM    gr, gc;                /* rows/columns of grid */
  DIM    row, col;              /* current row and column */
  REAL   value;                 /* value of current matrix element */
  int    err;                   /* error status */
  THREAD *thrds;                /* thread handles for parallelization */
  void   *work;                 /* data for worker threads */
  SFXNAME(FCMGETFN) *get;       /* element retrieval function */

#ifdef FCMAT_MAIN
  double  sum;                  /* total thread execution time */
  double  beg;                  /* loss due to thread start times */
  double  end;                  /* loss due to thread end   times */
  int     cnt;                  /* number of times cache filled */
#endif  

} SFXNAME(FCMAT);               /* (functional connectivity matrix) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern SFXNAME(FCMAT)*
            SFXNAME(fcm_create) (REAL *data, DIM V, DIM T,
                                 int mode, ...);
extern void SFXNAME(fcm_delete) (SFXNAME(FCMAT) *fcm);
extern DIM  SFXNAME(fcm_dim)    (SFXNAME(FCMAT) *fcm);

extern int  SFXNAME(fcm_first)  (SFXNAME(FCMAT) *fcm);
extern int  SFXNAME(fcm_next)   (SFXNAME(FCMAT) *fcm);
extern DIM  SFXNAME(fcm_row)    (SFXNAME(FCMAT) *fcm);
extern DIM  SFXNAME(fcm_col)    (SFXNAME(FCMAT) *fcm);
extern REAL SFXNAME(fcm_value)  (SFXNAME(FCMAT) *fcm);
extern REAL SFXNAME(fcm_error)  (SFXNAME(FCMAT) *fcm);

extern REAL SFXNAME(fcm_get)    (SFXNAME(FCMAT) *fcm, DIM row, DIM col);

extern void SFXNAME(fcm_show)   (SFXNAME(FCMAT) *fcm);

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#if   _FCM_PASS <= 0
#define fcm_dim(m)          ((m)->V)
#define fcm_row(m)          ((m)->row)
#define fcm_col(m)          ((m)->col)
#define fcm_value(m)        ((m)->value)
#define fcm_error(m)        ((m)->err)
#define fcm_get(m,r,c)      ((m)->get(m,r,c))
#elif _FCM_PASS <= 1
#define fcm_dim_flt(m)      ((m)->V)
#define fcm_row_flt(m)      ((m)->row)
#define fcm_col_flt(m)      ((m)->col)
#define fcm_value_flt(m)    ((m)->value)
#define fcm_error_flt(m)    ((m)->err)
#define fcm_get_flt(m,r,c)  ((m)->get(m,r,c))
#elif _FCM_PASS <= 2
#define fcm_dim_dbl(m)      ((m)->V)
#define fcm_row_dbl(m)      ((m)->row)
#define fcm_col_dbl(m)      ((m)->col)
#define fcm_value_dbl(m)    ((m)->value)
#define fcm_error_dbl(m)    ((m)->err)
#define fcm_get_dbl(m,r,c)  ((m)->get(m,r,c))
#endif

/*----------------------------------------------------------------------
  Recursion Handling
----------------------------------------------------------------------*/
#if   _FCM_PASS == 1            /* if in first of two passes */
#undef REAL
#undef SUFFIX
#include "fcmat.h"              /* process file recursively */
#elif _FCM_PASS == 2
#undef REAL
#endif

#undef SUFFIX
#undef SFXNAME
#undef SFXNAME_1
#undef SFXNAME_2
#undef REAL_IS_DOUBLE

#undef  _FCM_PASS
#define __FCMAT__
#endif
