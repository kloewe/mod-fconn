/*----------------------------------------------------------------------------
  File    : fcmat3.h
  Contents: data type for functional connectivity matrix (half-stored)
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
#ifndef FCMAT3_H

#include "stats.h"
#include "pcc.h"
#include "tetracc.h"
#include "fcmat.h"

/*----------------------------------------------------------------------------
  Data Type Definition / Recursion Handling
----------------------------------------------------------------------------*/
#ifndef DIM                     /* if matrix dimension is not defined */
#define DIM int                 /* use int as the default */
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

/*--------------------------------------------------------------------------*/

#ifndef SFXNAME                 /* macros to generate function names */
#define SFXNAME(n)      SFXNAME_1(n,SUFFIX)
#define SFXNAME_1(n,s)  SFXNAME_2(n,s)
#define SFXNAME_2(n,s)  n##s    /* the two step recursion is needed */
#endif                          /* to ensure proper expansion */

/*----------------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------------*/
#if   defined __AVX__           /* if AVX instructions available */
#  define INIT_PCC SFXNAME(init_avx)
#elif defined __SSE2__          /* if SSE2 instructions available */
#  define INIT_PCC SFXNAME(init_sse2)
#else                           /* if neither extension available */
#  define INIT_PCC SFXNAME(init_naive)
#endif                          /* fall back to naive computations */

#if defined __POPCNT__ && defined __SSE4_1__
// #  define PCAND_TCC pcand_m128i /* use m128i implementation */
#  define BPI       128         /* 128 bits per integer */
#else
// #  define PCAND_TCC pcand_lut16 /* use lut16 implemenation */
#  define BPI       32          /* 32 bits per integer */
#endif

#define INDEX(i,j,N)    ((size_t)(i)*((size_t)(N)+(size_t)(N) \
                         -(size_t)(i)-3)/2-1+(size_t)(j))

/*----------------------------------------------------------------------------
  Inline Retrieval Functions
----------------------------------------------------------------------------*/

inline REAL SFXNAME(fcm_full) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get corr.c. from full rep. */
  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return fcm->diag;           /* return a fixed value */
  return fcm->cache[(row > col) /* retrieve the correlation coeff. */
                    ? INDEX(col, row, fcm->V)
                    : INDEX(row, col, fcm->V)];
}  /* fcm_full() */

/*----------------------------------------------------------------------------
  Recursion Handling
----------------------------------------------------------------------------*/
#if   _FCM_PASS == 1            /* if in first of two passes */
#undef REAL
#undef SUFFIX
#include "fcmat3.h"             /* process file recursively */
#elif _FCM_PASS == 2
#undef REAL
#endif

#undef SUFFIX
#undef SFXNAME
#undef SFXNAME_1
#undef SFXNAME_2

#undef  _FCM_PASS

#define FCMAT3_H
#endif  /* #ifndef FCMAT3_H */
