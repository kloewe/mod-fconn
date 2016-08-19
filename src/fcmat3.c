/*----------------------------------------------------------------------------
  File    : fcmat3.c
  Contents: data type for functional connectivity matrix (half-stored)
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
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

#include "stats.h"
#include "cpuinfo.h"
#include "binarize.h"
#include "pcc.h"
#include "tetracc.h"
#include "fcmat3.h"
// #include "mex.h"

/*----------------------------------------------------------------------------
  Data Type Definition / Recursion Handling
----------------------------------------------------------------------------*/
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

/*--------------------------------------------------------------------------*/
#define float  1                /* to check the definition of REAL */
#define double 2

#if   REAL == float             /* if single precision data */
#undef  REAL_IS_DOUBLE
#define REAL_IS_DOUBLE  0       /* clear indicator for double */
#elif REAL == double            /* if double precision data */
#undef  REAL_IS_DOUBLE
#define REAL_IS_DOUBLE  1       /* set   indicator for double */
#else
#error "REAL must be either 'float' or 'double'"
#endif

#undef float                    /* delete definitions */
#undef double                   /* used for type checking */
/*--------------------------------------------------------------------------*/
#define int         1           /* to check definitions */
#define long        2           /* for certain types */
#define ptrdiff_t   3

#if   DIM == int
#ifndef DIM_FMT
#define DIM_FMT     "d"         /* printf format code for int */
#endif
#ifndef strtodim
#define strtodim(s,p)   (int)strtol(s,p,0)
#endif

#elif DIM == long
#ifndef DIM_FMT
#define DIM_FMT     "ld"        /* printf format code for long */
#endif
#ifndef strtodim
#define strtodim(s,p)   strtol(s,p,0)
#endif

#elif DIM == ptrdiff_t
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
/*--------------------------------------------------------------------------*/

#ifndef SFXNAME                 /* macros to generate function names */
#define SFXNAME(n)      SFXNAME_1(n,SUFFIX)
#define SFXNAME_1(n,s)  SFXNAME_2(n,s)
#define SFXNAME_2(n,s)  n##s    /* the two step recursion is needed */
#endif                          /* to ensure proper expansion */

/*----------------------------------------------------------------------------
  Function Prototypes
----------------------------------------------------------------------------*/
extern REAL SFXNAME(fcm_full) (SFXNAME(FCMAT) *fcm, DIM row, DIM col);

/*----------------------------------------------------------------------------
  Functions
----------------------------------------------------------------------------*/

SFXNAME(FCMAT)* SFXNAME(fcm_create) (REAL *data, DIM V, DIM T, int mode, ...)
{                               /* --- create a func. connect. matrix */
  SFXNAME(FCMAT) *fcm;          /* func. con. matrix to be created */
  va_list args;                 /* list of variable arguments */
  size_t  k, z;                 /* loop variable, cache size */

  assert(data && (V > 1) && (T > 1)); /* check the function arguments */
  fcm = (SFXNAME(FCMAT)*)malloc(sizeof(SFXNAME(FCMAT)));
  if (!fcm) return NULL;        /* allocate the base structure */
  fcm->V       = V;             /* note the number of voxels */
  fcm->T       = T;             /* and  the number of scans */
  fcm->X       = T;             /* default: data blocks like scans */
  fcm->mode    = mode;          /* note the processing mode */
  fcm->nthd    = proccnt();     /* default: use all processors */
  fcm->mem     = NULL;          /* clear memory block, */
  fcm->cmap    = NULL;          /* cosine map and cache */
  fcm->cache   = NULL;          /* for easier cleanup */
  fcm->diag    = (REAL)((mode & FCM_R2Z) ? R2Z_MAX : 1.0);
  fcm->value   = fcm->diag;     /* set the value and coordinates */
  fcm->row     = fcm->col = 0;  /* of the current element */
  fcm->ra      = 0; fcm->rb = V;
  fcm->ca      = 0; fcm->cb = V;/* init. the coordinate ranges */
  fcm->err     = 0;             /* clear the error status */

  va_start(args, mode);         /* start variable arguments */
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
    fcm->mem = malloc((size_t)V*(size_t)fcm->X *sizeof(REAL) +31);
    if (!fcm->mem) { SFXNAME(fcm_delete)(fcm); return NULL; }
    fcm->data = (REAL*)(((uintptr_t)fcm->mem +31) & ~(uintptr_t)31);
    INIT_PCC(data, (int)V, (int)T, fcm->data, (int)fcm->X);
  }
  else if (mode == FCM_TCC) {   /* if tetrachoric correlation coeff. */
    #if defined __POPCNT__ \
    &&  defined __SSE4_1__      /* use 128 bit and popcnt if possible */
    fcm->X = 4 *(((int)T+127) >> 7);
    #else                       /* otherwise fall back to LUT16 */
    fcm->X = ((int)T+31) >> 5;  /* get block size of binarized data */
    #endif
    fcm->mem  =                 /* allocate memory for binarized data */
    fcm->data = SFXNAME(binarize)(data, (int)V, (int)T, BIN_MEDIAN, BPI);
    if (!fcm->data) { SFXNAME(fcm_delete)(fcm); return NULL; };
    fcm->cmap = SFXNAME(make_cmap)((int)T); /* create cosine map */
    if (!fcm->cmap) { SFXNAME(fcm_delete)(fcm); return NULL; };
    init_popcnt(); }            /* initialize bit count table */
  else {                        /* if unknown correlation variant */
    fprintf(stderr, "fcm_create: unknown correlation variant\n");
    free(fcm); return NULL;     /* print an error message */
  }                             /* and abort the function */

  z = (size_t)V *(size_t)(V-1)/2; /* cache for upper triangle */
  fcm->cache = (REAL*)malloc(z *sizeof(REAL));
  if (!fcm->cache) { SFXNAME(fcm_delete)(fcm); return NULL; }
  if (mode == FCM_PCC)          /* if Pearson correlation coefficient */
    SFXNAME(pccx)    (data, fcm->cache, (int)V, (int)T,
                  // PCC_AVX|PCC_COBL|PCC_THREAD, 0, fcm->nthd);
                      PCC_AUTO|PCC_THREAD, fcm->nthd);
  else                          /* if tetrachoric correlation coeff. */
    SFXNAME(tetraccx)(data, fcm->cache, (int)V, (int)T,
                      TCC_AUTO|TCC_THREAD, fcm->nthd);

  fcm->cget = SFXNAME(fcm_full);/* retrieval function for fcm_next() */
  fcm->get = SFXNAME(fcm_full); /* retrieval function for fcm_get() */

  if (fcm->mode & FCM_R2Z)      /* if Fisher's r-to-z transform, */
    for (k = 0; k < z; k++)     /* transform the matrix elements */
      fcm->cache[k] = SFXNAME(fisher_r2z)(fcm->cache[k]);
                                /* set the element retrieval function */

  return fcm;                   /* return created FC matrix */
}  /* fcm_create() */

/*--------------------------------------------------------------------------*/

void SFXNAME(fcm_delete) (SFXNAME(FCMAT) *fcm)
{                               /* --- delete a func. connect. matrix */
  assert(fcm);                  /* check the function argument */
  if (fcm->cache)  free(fcm->cache);
  if (fcm->cmap)   free(fcm->cmap);
  if (fcm->mem)    free(fcm->mem);/* delete cache, cosine map, */
  free(fcm);                    /* data, and the base structure */
}  /* fcm_delete() */

/*--------------------------------------------------------------------------*/

int SFXNAME(fcm_first) (SFXNAME(FCMAT) *fcm)
{                               /* --- get first matrix element */
  assert(fcm);                  /* check the function argument */
  fcm->err = 0;                 /* clear the error status */
  fcm->row = 0; fcm->col = 1;   /* start with first off-diag. element */
  fcm->value = fcm->cget(fcm, 0, 1);
  return (fcm->err) ? -1 : +1;  /* store the value and return status */
}  /* fcm_first() */

/*--------------------------------------------------------------------------*/

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
  fcm->value = fcm->cget(fcm, fcm->row, fcm->col);
  return (fcm->err) ? -1 : +1;  /* store the value and return status */
}  /* fcm_next() */

/*--------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
  Recursion Handling
----------------------------------------------------------------------------*/
#if _FCM_PASS == 1              /* if in first of two passes */
#undef REAL
#undef SUFFIX
#include "fcmat3.c"              /* process file recursively */
#endif
