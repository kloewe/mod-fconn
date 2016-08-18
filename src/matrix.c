/*----------------------------------------------------------------------
  File    : matrix.c
  Contents: data type for a square matrix
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include "matrix.h"

#ifndef INFINITY
#define INFINITY    (DBL_MAX+DBL_MAX)
#endif                          /* MSC still does not support C99 */
#ifndef NAN
#define NAN         (INFINITY-INFINITY)
#endif                          /* MSC still does not support C99 */

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#ifndef LOWER                   /* whether the linear index is based */
#define LOWER       0           /* on a lower triangular matrix */
#endif

#define BLKSIZE    256          /* block size for arrays */

#if LOWER                       /* if lower triangular matrix */
#define INDEX(i,j,N)    ((size_t)(i)*(size_t)((i)-1)/2+(size_t)(j))
#else                           /* if upper triangular matrix */
#define INDEX(i,j,N)    ((size_t)(i)*((size_t)(N)+(size_t)(N) \
                        -(size_t)(i)-3)/2-1+(size_t)(j))
#endif                          /* index computation for result */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/

MATRIX* mat_create (DIM n, size_t size)
{                               /* --- create a square matrix */
  MATRIX *mat;                  /* created matrix */

  mat = (MATRIX*)malloc(sizeof(MATRIX));
  if (!mat) return NULL;        /* create the matrix body */
  mat->size = size;             /* note matrix size/type and */
  mat->n    = n;                /* matrix dimension (n-by-n matrix) */
  if (size <= 0) {              /* if full matrix */
    size = (size_t)n *(size_t)(n-1)/2;
    mat->idxs = NULL; }         /* compute full matrix size */
  else {                        /* if sparse matrix */
    mat->idxs = (size_t*)malloc(size *sizeof(size_t));
    if (!mat->idxs) { free(mat); return NULL; }
  }                             /* create the index array */
  mat->elems = (REAL*)malloc(size *sizeof(REAL));
  if (!mat->elems) { mat_delete(mat); return NULL; }
  mat->cnt = 0;                 /* create the matrix elements */
  return mat;                   /* return the created matrix */
}  /* mat_create() */

/*--------------------------------------------------------------------*/

void mat_delete (MATRIX *mat)
{                               /* --- create a square matrix */
  if (mat->idxs) free(mat->idxs);
  free(mat->elems);             /* delete index array and elements */
  free(mat);                    /* delete the base structure */
}  /* mat_delete() */

/*--------------------------------------------------------------------*/

REAL mat_get (const MATRIX *mat, DIM i, DIM j)
{                               /* --- get a matrix element */
  size_t k;                     /* linear index */
  size_t l, r, m;               /* index array indices */

  if (j == i) return (REAL)NAN; /* diag. elems are not defined */ 
  #if LOWER
  if (j > i)                    /* compute the linear index */
    k = INDEX(j, i, mat->n);    /* A[j][i] = A[i][j]; */
  #else                         /* index computation depends on */
  if (j < i)                    /* appropriate combinations of i and j */
    k = INDEX(j, i, mat->n);    /* wrt. upper vs. lower triangle */
  #endif
  else
    k = INDEX(i, j, mat->n);
  if (mat->size <= 0)           /* if full matrix */
    return mat->elems[k];       /* get the element directly */
  l = 0; r = mat->cnt;          /* init. bounds for binary search */
  while (l < r) {               /* while search range is not empty */
    m = (l+r)/2;                        /* compare the given index */
    if      (k > mat->idxs[m]) l = m+1; /* to the middle element */
    else if (k < mat->idxs[m]) r = m;   /* adapt the search range */
    else return mat->elems[m];          /* according to the result */
  }                             /* if match found, return element */
  return (REAL)NAN;             /* return 'element not found' */
}  /* mat_get() */

/*--------------------------------------------------------------------*/

REAL mat_set (MATRIX *mat, DIM i, DIM j, REAL val)
{                               /* --- set a matrix element */
  size_t k;                     /* linear index */
  size_t n;                     /* new array size */
  void   *p;                    /* buffer for reallocation */

  #if LOWER
  if (j > i)                    /* compute the linear index */
    k = INDEX(j, i, mat->n);    /* A[j][i] = A[i][j]; */
  #else                         /* index computation depends on */
  if (j < i)                    /* appropriate combinations of i and j */
    k = INDEX(j, i, mat->n);    /* wrt. upper vs. lower triangle */
  #endif
  else
    k = INDEX(i, j, mat->n);
  if (mat->size <= 0)           /* if full matrix */
    mat->elems[k] = val;        /* set the element directly */
  else {                        /* if sparse matrix */
    assert((mat->cnt <= 0) || (k > mat->idxs[mat->cnt-1]));
    n = mat->size;              /* get the index array size */
    if (mat->cnt >= n) {        /* if the index array is full */
      n += (n > BLKSIZE) ? n >> 1 : BLKSIZE;
      p = realloc(mat->idxs, n *sizeof(size_t));
      if (!p) return (REAL)NAN; /* enlarge the index array */
      mat->idxs = (size_t*)p;   /* and set the new array */
      p = realloc(mat->elems, n *sizeof(REAL));
      if (!p) return (REAL)NAN; /* enlarge the element array */
      mat->elems = (REAL*)p;    /* and set the new array */
      mat->size  = n;           /* note the new array size */
    }
    mat->idxs [mat->cnt]   = k; /* store the linear index */
    mat->elems[mat->cnt++] = val;
  }                             /* store the matrix element */
  return val;                   /* return the value that was set */
}  /* mat_set() */

/*--------------------------------------------------------------------*/

void mat_print (const MATRIX *mat)
{                               /* --- print a square matrix */
  size_t i, j;                  /* matrix element indices */
  size_t k, p;                  /* linear index, position */
  REAL   x;                     /* matrix element to print */

  if (mat->size <= 0) {         /* if full matrix */
    // for (i = k = 0; i < mat->n; i++) {
    //   #if LOWER
    //   for (j = 0; j < i; j++) { /* traverse rows, then columns */
    //     if (j > 0) fputc(' ', stdout);
    //   #else
    //   for (j = i; j < mat->n; j++) {
    //     if (j > i) fputc(' ', stdout);
    //   #endif
    //   printf("%g", mat->elems[k++]);
    //   }                         /* print the matrix elements */
    // }
    printf("\n     ");
    for (j = 0; j < mat->n; j++)
      printf("%-8zu ", j);
    printf("\n");
    for (i = 0; i < mat->n; i++) {
      printf("%-5zu", i);
      for (j = 0; j < mat->n; j++) {
        printf("%-8.3g ", mat_get(mat, i, j));
      }
      printf("\n");
    }
  } else {                    /* if sparse matrix */
    for (i = k = p = 0; i < mat->n; i++) {
      #if LOWER
      for (j = 0; j < i; j++) { /* traverse rows, then columns */
        if (j > 0) fputc(' ', stdout);
      #else
      for (j = i; j < mat->n; j++) {
        if (j > i) fputc(' ', stdout);
      #endif
        x = ((p < mat->cnt) && (mat->idxs[p] == k))
            ? mat->elems[k] : (REAL)NAN;
        printf("%g", x); k++;   /* get the matrix elements */
      }                         /* or NAN if it does not exist */
    }                           /* and then print the value */
  }
}  /* mat_print() */
