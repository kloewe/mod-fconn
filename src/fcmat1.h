/*----------------------------------------------------------------------------
  File    : fcmat1.h
  Contents: data type for functional connectivity matrix (on-demand)
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
#ifndef FCMAT1_H

#include "stats.h"
#include "pcc.h"
#include "binarize.h"
#include "tetracc.h"
#include "fcmat.h"

/*----------------------------------------------------------------------------
  Data Type Definition / Recursion Handling
----------------------------------------------------------------------------*/
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

#ifndef SFXNAME                 /* macros to generate function names */
#define SFXNAME(n)      SFXNAME_1(n,SUFFIX)
#define SFXNAME_1(n,s)  SFXNAME_2(n,s)
#define SFXNAME_2(n,s)  n##s    /* the two step recursion is needed */
#endif                          /* to ensure proper expansion */

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
  Inline Iteration Functions
----------------------------------------------------------------------------*/

inline int SFXNAME(fcm_first) (SFXNAME(FCMAT) *fcm)
{                               /* --- get first matrix element */
  assert(fcm);                  /* check the function argument */
  fcm->err = 0;                 /* clear the error status */
  fcm->row = 0; fcm->col = 1;   /* start with first off-diag. element */
  fcm->value = fcm->cget(fcm, 0, 1);
  return (fcm->err) ? -1 : +1;  /* store the value and return status */
}  /* fcm_first() */

/*--------------------------------------------------------------------------*/

inline int SFXNAME(fcm_next) (SFXNAME(FCMAT) *fcm)
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

/*----------------------------------------------------------------------------
  Functions
----------------------------------------------------------------------------*/

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
