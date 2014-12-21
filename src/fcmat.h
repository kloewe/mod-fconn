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
