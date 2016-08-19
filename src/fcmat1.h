/*----------------------------------------------------------------------------
  File    : fcmat1.h
  Contents: data type for functional connectivity matrix (on-demand)
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
#ifndef FCMAT1_H

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
#  define PAIR_PCC SFXNAME(pair_avx)
#elif defined __SSE2__          /* if SSE2 instructions available */
#  define INIT_PCC SFXNAME(init_sse2)
#  define PAIR_PCC SFXNAME(pair_sse2)
#else                           /* if neither extension available */
#  define INIT_PCC SFXNAME(init_naive)
#  define PAIR_PCC SFXNAME(pair_naive)
#endif                          /* fall back to naive computations */

#if defined __POPCNT__ && defined __SSE4_1__
#  define PCAND_TCC pcand_m128i /* use m128i implementation */
#  define BPI       128         /* 128 bits per integer */
#else
#  define PCAND_TCC pcand_lut16 /* use lut16 implemenation */
#  define BPI       32          /* 32 bits per integer */
#endif

/*----------------------------------------------------------------------------
  Inline Retrieval Functions
----------------------------------------------------------------------------*/

inline REAL SFXNAME(fcm_pccotf) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get Pearson cc. on the fly */
  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return (REAL)+1;            /* always return +1.0 */
  if (row >  col) {             /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  return PAIR_PCC((REAL*)fcm->data +(size_t)row*(size_t)fcm->X,
                  (REAL*)fcm->data +(size_t)col*(size_t)fcm->X,
                  (int)fcm->T); /* compute Pearson correlation coeff. */
}  /* fcm_pccotf() */

/*--------------------------------------------------------------------------*/

inline REAL SFXNAME(fcm_pccr2z) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get Pearson cc. on the fly */
  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return (REAL)+R2Z_MAX;      /* return atanh(1-epsilon) */
  if (row > col) {              /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  REAL r = PAIR_PCC((REAL*)fcm->data +(size_t)row*(size_t)fcm->X,
                    (REAL*)fcm->data +(size_t)col*(size_t)fcm->X,
                    (int)fcm->T);/* compute Pearson correlation coeff. */
  return fisher_r2z(r);         /* apply Fisher's r to z transform */
}  /* fcm_pccr2z() */

/*--------------------------------------------------------------------------*/

inline REAL SFXNAME(fcm_tccotf) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get tetrachoric cc. on the fly */
  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return (REAL)+1.0;          /* always return +1.0 */
  if (row > col) {              /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  int n = PCAND_TCC((uint32_t*)fcm->data +(size_t)row*(size_t)fcm->X,
                    (uint32_t*)fcm->data +(size_t)col*(size_t)fcm->X,
                    (int)fcm->X);/* count number of 11 configs. and */
  return fcm->cmap[n];          /* compute tetrachoric corr. coeff. */
}  /* fcm_tccotf() */

/*--------------------------------------------------------------------------*/

inline REAL SFXNAME(fcm_tccr2z) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get tetrachoric cc. on the fly */
  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return (REAL)+R2Z_MAX;      /* return atanh(1-epsilon) */
  if (row > col) {              /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  int n = PCAND_TCC((uint32_t*)fcm->data +(size_t)row*(size_t)fcm->X,
                    (uint32_t*)fcm->data +(size_t)col*(size_t)fcm->X,
                    (int)fcm->X);/* compute tetrachoric corr. coeff. */
  return fisher_r2z(fcm->cmap[n]);
}  /* fcm_tccr2z() */           /* apply Fisher's r to z transform */


/*----------------------------------------------------------------------------
  Recursion Handling
----------------------------------------------------------------------------*/
#if   _FCM_PASS == 1            /* if in first of two passes */
#undef REAL
#undef SUFFIX
#include "fcmat1.h"             /* process file recursively */
#elif _FCM_PASS == 2
#undef REAL
#endif

#undef SUFFIX
#undef SFXNAME
#undef SFXNAME_1
#undef SFXNAME_2

#undef  _FCM_PASS

#define FCMAT1_H
#endif  /* #ifndef FCMAT1_H */
