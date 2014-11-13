/*----------------------------------------------------------------------
  File    : matrix.h
  Contents: data type for a square matrix
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
#ifndef __MATRIX__
#define __MATRIX__

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#ifndef REAL
#define REAL float              /* type of matrix elements */
#endif
#ifndef DIM
#define DIM  size_t             /* type of matrix dimension(s) */
#endif

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- a square matrix --- */
  DIM    n;                     /* matrix size (n-by-n matrix) */
  REAL   *elems;                /* matrix elements (upper triangle?) */
  size_t *idxs;                 /* linear index array if sparse */
  size_t cnt;                   /* current number of indices */
  size_t size;                  /* current size of the index array */
} MATRIX;                       /* (square matrix) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern MATRIX* mat_create (DIM n, size_t size); /* n = nrows = ncols */
extern void    mat_delete (MATRIX *mat);
extern DIM     mat_dim    (const MATRIX *mat);
extern REAL    mat_get    (const MATRIX *mat, DIM i, DIM j);
extern REAL    mat_set    (MATRIX *mat, DIM i, DIM j, REAL val);
extern void    mat_print  (const MATRIX *mat);

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define mat_dim(m)        ((m)->n)

#endif
