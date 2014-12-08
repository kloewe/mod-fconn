/*----------------------------------------------------------------------------
  File    : fcmat.c
  Contents: data type for functional connectivity matrix
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>

#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __AVX__
#include <immintrin.h>
#endif
#ifdef __SSE4_1__
#include <smmintrin.h>
#endif
#ifdef __POPCNT__
#include <popcntintrin.h>
#endif

#include "pcc.h"
#include "tetracc.h"
#include "binarize.h"
#include "stats.h"
#include "fcmat.h"

/*----------------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------------*/
#define rssd buf[0]               /* roots of sums of squared deviats. */
#define diff buf[1]               /* (differences to) mean values */
#define cmap buf[0]               /* map for cosine function */
#define bits buf[1]               /* binarized data */

#ifdef __AVX__                    /* if AVX instructions available */
#define INIT_PCC init_avx
#define PAIR_PCC pair_avx
#else
#ifdef __SSE2__                   /* if SSE2 instructions available */
#define INIT_PCC init_sse2
#define PAIR_PCC pair_sse2
#else
#error "SSE2 not supported!"      /* else it's time to buy a new computer */
#endif
#endif

#if defined __POPCNT__ && defined __SSE4_1__
  #define PCAND_TCC pcand_m128i   /* use m128i implementation */
  #define BPI 128                 /* 128 bits per integer */
#else
  #define PCAND_TCC pcand_lut16   /* use lut16 implemenation */
  #define BPI 32                  /* 32 bits per integer */
#endif

/*----------------------------------------------------------------------------
  Functions
----------------------------------------------------------------------------*/

FCMAT* fcm_create (REAL *data, DIM V, DIM T, int flags)
{
  assert((V > 1) && (T > 1));
  
  FCMAT *fcm;                     /* FC matrix to be created */
  
  fcm = (FCMAT*)malloc(sizeof(FCMAT));
  if (!fcm) return NULL;
  fcm->V = V;
  fcm->T = T;
  fcm->data = data;
  fcm->flags = flags;
  
  if        ((flags & FCM_CORR) == FCM_PCC) {
    #ifdef __AVX__                /* if available use AVX */
    fcm->X = (((int)T+7) & ~7);
    #else                         /* else fallback to SSE2 */
    fcm->X = (((int)T+3) & ~3);
    #endif
    fcm->rssd = malloc(           /* malloc and compute rssd and diff */
             ((size_t)V*(size_t)fcm->X +(size_t)V) *sizeof(REAL) +31);
    fcm->diff = (REAL*)(((ptrdiff_t)((REAL*)(fcm->rssd) +V) +31) & ~31);
    INIT_PCC(data, (int)V, (int)T, fcm->diff, fcm->X, fcm->rssd);
    
  } else if ((flags & FCM_CORR) == FCM_TCC) {
    #if defined __POPCNT__ \
            && defined __SSE4_1__ /* if possible use 128 bit version */
    fcm->X = 4 *(((int)T+127) >> 7);
    #else                         /* else fallback to LUT16 version */
    fcm->X = ((int)T+31) >> 5;
    #endif
    fcm->bits = (uint32_t *)binarize(data, (int)V, (int)T, BIN_MEDIAN, BPI);
    fcm->cmap = make_cmap((int)T);/* make map from n_11 to r (avoid cos */
    if (!fcm->cmap) { free(fcm->bits); return NULL; }; /* computations) */
    init_popcnt();                /* initialize bit count table */
    
  } else {
    printf("Error: Unexpected flag.\n");
    return NULL;
  }
  return fcm;                     /* return created FC matrix */
} /* fcm_create() */

/*--------------------------------------------------------------------------*/

void fcm_delete (FCMAT *fcm)
{
  if ((fcm->flags & FCM_CORR) == FCM_TCC) {
    free(fcm->buf[0]);            /* tetrachoric corr. coef. */
    free(fcm->buf[1]);
  } else                          /* Pearson corr. coef. */
    free(fcm->buf[0]);
  free(fcm);
} /* fcm_delete() */

/*--------------------------------------------------------------------------*/

REAL fcm_get(FCMAT *fcm, DIM i, DIM j)
{
  REAL r;
  if (i==j)
    r = 1.0f;
  else {
    if ((fcm->flags & FCM_CORR) == FCM_TCC) {
      r = ((REAL*)fcm->cmap)[     /* tetrachoric corr. coef. */
              PCAND_TCC((uint32_t*)(fcm->bits)+i*(size_t)(fcm->X),
                        (uint32_t*)(fcm->bits)+j*(size_t)(fcm->X), fcm->X)];
    } else {                      /* Pearson corr. coef. */
      r = PAIR_PCC((REAL*)(fcm->diff)+i*(size_t)(fcm->X),
                      (REAL*)(fcm->diff)+j*(size_t)(fcm->X),
                      (int)(fcm->T))
                        /(((REAL*)fcm->rssd)[i] * ((REAL*)(fcm->rssd))[j]);
    }
  }
  
  if (fcm->flags & FCM_R2Z) {     /* apply r-to-z transform */
    return fisher_r2z(r);
  } else {                        /* or return raw correlation */
    return r;
  }
} /* fcm_get() */

/*--------------------------------------------------------------------------*/

void fcm_print(FCMAT *fcm)
{
  printf("     ");
  for (DIM j = 0; j < fcm->V; j++)
    printf("%-8zu ", j);
  printf("\n");
  for (DIM i = 0; i < fcm->V; i++) {
    printf("%-5zu", i);
    for (DIM j = 0; j < fcm->V; j++) {
      printf("%-8.3g ", fcm_get(fcm, i, j));
    }
    printf("\n");
  }
} /* fcm_print() */

/*--------------------------------------------------------------------------*/
#ifdef FCMAT_MAIN
int main (int argc, char* argv[])
{                                 /* --- main function for testing */
  
  int      V, T;                  /* number of voxels and time points */
  REAL     *data;                 /* data and result arrays */
  if (argc != 3) {                /* check the function arguments */
    fprintf(stderr, "usage: %s V T\n", argv[0]); return 0; }
  V = atoi(argv[1]);              /* get the number of voxels */
  T = atoi(argv[2]);              /* and the number of time points */
  
#if 1
  srand((unsigned)time(NULL));    /* seed the random number generator */
#else
  srand(1);                       /* seed the random number generator */
#endif
  
  data = malloc((size_t)V *(size_t)T *sizeof(REAL));
  if (!data) {                    /* malloc data array and generate */
    fprintf(stderr, "%s: not enough memory\n", argv[0]); return -1; }
  for (DIM i = 0; i < (DIM)(T*V); i++) /* random numbers */
    data[i] = (REAL)(rand()/((double)RAND_MAX+1));
  
  FCMAT *fcm1, *fcm2;
  fcm1 = fcm_create(data, (DIM)V, (DIM)T, FCM_PCC);
  fcm2 = fcm_create(data, (DIM)V, (DIM)T, FCM_TCC);
  fcm_print(fcm1);
  fcm_print(fcm2);
  fcm_delete(fcm1);
  fcm_delete(fcm2);
  
  free(data);
}  /* main() */

#endif  /* #ifdef FCMAT_MAIN */

#undef FCMAT_MAIN
