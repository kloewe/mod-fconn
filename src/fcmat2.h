/*----------------------------------------------------------------------------
  File    : fcmat2.h
  Contents: data type for functional connectivity matrix (cache-based)
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------------*/
#ifndef FCMAT2_H

#include "fcmat1.h"

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

/*----------------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------------*/
#define TILE_MIN    16          /* minimum tile size for comp. */
#define RECTANGLE   1           /* shape id for rectangular area */
#define TRIANGLE    2           /* shape id for triangular  area */
// #define PAIRSPLIT
// #define RECTGRID
// #define SAFETHREAD

/*--------------------------------------------------------------------------*/
#ifndef THREAD_OK               /* if not yet defined */
#ifdef _WIN32                   /* if Microsoft Windows system */
#define THREAD_OK       0       /* return value is DWORD */
#define WORKERDEF(n,p)  DWORD WINAPI SFXNAME(n) (LPVOID p)
#else                           /* if Linux/Unix system */
#define THREAD_OK       NULL    /* return value is void* */
#define WORKERDEF(n,p)  void*        SFXNAME(n) (void* p)
#endif                          /* definition of a worker function */
#endif

/*----------------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------------*/
typedef struct {                /* --- thread worker data --- */
  int    work;                  /* flag for assigned work */
  DIM    ra, rb;                /* row    index range */
  DIM    ca, cb;                /* column index range */
  DIM    cm;                    /* reference column for mirroring */
  SFXNAME(FCMAT)    *fcm;       /* underlying f.c. matrix object */
  SFXNAME(FCMGETFN) *get;       /* element computation function */
  #ifdef FCM_BENCH              /* if to do some benchmarking */
  double beg;                   /* start time of thread */
  double end;                   /* end   time of thread */
  #endif                        /* (for time loss computation) */
} SFXNAME(WORK);                /* (thread worker data) */

#ifndef WORKERTYPE              /* if not yet defined */
#define WORKERTYPE              /* define worker function type */
#ifdef _WIN32                   /* if Microsoft Windows system */
typedef DWORD WINAPI WORKER (LPVOID);
#else                           /* if Linux/Unix system */
typedef void*        WORKER (void*);
#endif                          /* worker for parallel execution */
#endif

/*----------------------------------------------------------------------------
  Timer Function
----------------------------------------------------------------------------*/
#if (defined FCM_BENCH || defined FCMAT_MAIN) && !defined TIMER
#define TIMER

static double timer (void)
{                               /* --- get current time */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  return 0.001 *(double)timeGetTime();
  #else                         /* if Linux/Unix system */
  struct timespec tp;           /* POSIX time specification */
  clock_gettime(CLOCK_MONOTONIC, &tp);
  return (double)tp.tv_sec +1e-9 *(double)tp.tv_nsec;
  #endif                        /* return time in seconds */
}  /* timer() */

#endif

/*----------------------------------------------------------------------------
  Cache Filling Functions
----------------------------------------------------------------------------*/

inline REAL SFXNAME(pcc_pure) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- compute Pearson corr. coeff. */
  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col > row) && (col < fcm->V));
  return SFXNAME(clamp)(dot((REAL*)fcm->data +(size_t)row *(size_t)fcm->X,
                  (REAL*)fcm->data +(size_t)col *(size_t)fcm->X,
                  (int)fcm->T), (REAL)-1, (REAL)+1); /* compute Pearson cc. */
}  /* pcc_pure() */

/*--------------------------------------------------------------------------*/

inline REAL SFXNAME(pcc_r2z) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get Pearson corr. coeff. */
  REAL r;                       /* correlation coefficient */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col > row) && (col < fcm->V));
  r = SFXNAME(clamp)(dot((REAL*)fcm->data +(size_t)row *(size_t)fcm->X,
               (REAL*)fcm->data +(size_t)col *(size_t)fcm->X,
               (int)fcm->T), (REAL)-1, (REAL)+1); /* compute Pearson cc. */
  return fr2z(r);         /* apply Fisher's r to z transform */
}  /* pcc_r2z() */

/*--------------------------------------------------------------------------*/

inline REAL SFXNAME(tcc_pure) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- compute tetrachoric corr. cf. */
  int n;                        /* number of 11 configurations */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col > row) && (col < fcm->V));
  n = PCAND_TCC((uint32_t*)fcm->data +(size_t)row *(size_t)fcm->X,
                (uint32_t*)fcm->data +(size_t)col *(size_t)fcm->X,
                (int)fcm->X);   /* count 11 configurations and */
  return fcm->cmap[n];          /* compute tetrachoric corr. coeff. */
}  /* tcc_pure() */

/*--------------------------------------------------------------------------*/

inline REAL SFXNAME(tcc_r2z) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- compute tetrachoric corr. cf. */
  int n;                        /* number of 11 configurations */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col > row) && (col < fcm->V));
  n = PCAND_TCC((uint32_t*)fcm->data +(size_t)row *(size_t)fcm->X,
                (uint32_t*)fcm->data +(size_t)col *(size_t)fcm->X,
                (int)fcm->X);   /* compute tetrachoric corr. coeff. */
  return fr2z(fcm->cmap[n]);
}  /* tcc_r2z() */              /* apply Fisher's r to z transform */

/*--------------------------------------------------------------------------*/
#ifdef PAIRSPLIT                /* --- split rectangle into 2 parts */

inline void SFXNAME(rec_rct) (SFXNAME(WORK) *w,
                              DIM ra, DIM rb, DIM ca, DIM cb)
{                               /* --- process a rectangle */
  DIM i, j;                     /* loop variables */

  assert(w                      /* check the function arguments */
  &&    (ra >= 0) && (ra <  w->fcm->V)
  &&    (rb > ra) && (rb <= w->fcm->V)
  &&    (ca >= 0) && (ca <  w->fcm->V)
  &&    (cb > ca) && (cb <= w->fcm->V));
  i = rb-ra; j = cb-ca;         /* compute row and column range */
  if ((i > TILE_MIN)            /* if at least one direction */
  ||  (j > TILE_MIN)) {         /* is larger than the minimum size */
    if (i > j) {                /* if there are more rows */
      i = (ra+rb)/2;            /* split the tile horizontally */
      SFXNAME(rec_rct)(w, ra, i, ca, cb);
      SFXNAME(rec_rct)(w, i, rb, ca, cb); }
    else {                      /* if there are more columns */
      j = (ca+cb)/2;            /* split the tile vertically */
      SFXNAME(rec_rct)(w, ra, rb, ca, j);
      SFXNAME(rec_rct)(w, ra, rb, j, cb);
    } }                         /* process the parts recursively */
  else {                        /* if no larger than minimum size */
    REAL *cache = w->fcm->cache;/* get the cache array */
    DIM  r = w->fcm->ra;        /* and the coordinates */
    DIM  c = w->fcm->ca;        /* of the reference element */
    for (i = ra; i < rb; i++) { /* traverse the rows */
      for (j = ca; j < cb; j++) /* traverse the columns */
        cache[(size_t)(i-r) *(size_t)w->fcm->tile +(size_t)(j-c)]
          = w->get(w->fcm, i, j);
    }                           /* compute correlation coefficient */
  }                             /* and store it in the cache */
}  /* rec_rct() */

/*--------------------------------------------------------------------------*/
#else                           /* --- split rectangle into 4 parts */

inline void SFXNAME(rec_rct) (SFXNAME(WORK) *w,
                              DIM ra, DIM rb, DIM ca, DIM cb)
{                               /* --- process a rectangle */
  assert(w                      /* check the function arguments */
  &&    (ra >= 0) && (ra <  w->fcm->V)
  &&    (rb > ra) && (rb <= w->fcm->V)
  &&    (ca >= 0) && (ca <  w->fcm->V)
  &&    (cb > ca) && (cb <= w->fcm->V));
  if (cb-ca > TILE_MIN) {       /* if larger than minimum size */
    DIM i = (ra+rb)/2;          /* halven the tile size and */
    DIM j = (ca+cb)/2;          /* process parts recursively */
    SFXNAME(rec_rct)(w, ra, i, ca, j);
    SFXNAME(rec_rct)(w, ra, i, j, cb);
    SFXNAME(rec_rct)(w, i, rb, j, cb);
    SFXNAME(rec_rct)(w, i, rb, ca, j); }
  else {                        /* if no larger than minimum size */
    REAL *cache = w->fcm->cache;/* get the cache array */
    DIM  i, r = w->fcm->ra;     /* and the coordinates */
    DIM  j, c = w->fcm->ca;     /* of the reference element */
    for (i = ra; i < rb; i++) { /* traverse the rows */
      for (j = ca; j < cb; j++) /* traverse the columns */
        cache[(size_t)(i-r) *(size_t)w->fcm->tile +(size_t)(j-c)]
          = w->get(w->fcm, i, j);
    }                           /* compute correlation coefficient */
  }                             /* and store it in the cache */
}  /* rec_rct() */

#endif
/*--------------------------------------------------------------------------*/

inline void SFXNAME(rec_trg) (SFXNAME(WORK) *w, DIM a, DIM b)
{                               /* --- recursively process a triangle */
  assert(w                      /* check the function arguments */
  &&    (a >= 0) && (a <= w->fcm->V)
  &&    (b >= a) && (b <= w->fcm->V));
  if (b-a > TILE_MIN) {         /* if larger than minimal size */
    DIM i = (a+b)/2;            /* halven the area size and */
    SFXNAME(rec_trg)(w, a, i);  /* split into three parts */
    SFXNAME(rec_rct)(w, a, i, i, b);
    SFXNAME(rec_trg)(w, i, b); }
  else {                        /* if no larger than min. tile size */
    REAL *cache = w->fcm->cache;/* get the cache array */
    DIM  i, r = w->fcm->ra;     /* and the coordinates */
    DIM  j, c = w->fcm->ca;     /* of the reference element */
    for (i = a; i < b; i++) {   /* traverse the rows */
      for (j = i+1; j < b; j++) /* traverse the columns */
        cache[(size_t)(i-r) *(size_t)w->fcm->tile +(size_t)(j-c)]
          = w->get(w->fcm, i, j);
    }                           /* compute the correlation coeff. */
  }                             /* and store it in the cache */
}  /* rec_trg() */

/*--------------------------------------------------------------------------*/
#ifdef PAIRSPLIT                /* --- split rectangle into 2 parts */

inline WORKERDEF(fill, p)
{                               /* --- fill (part of) the cache */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  DIM a;                        /* start index for second strip */

  assert(p);                    /* check the function argument */
  #ifndef _WIN32                /* not yet available for Windows */
  while (1) {                   /* work loop for blocked threads */
    if (!w->fcm->join) {        /* if to block and signal threads */
      pthread_mutex_lock(&w->fcm->mutex);
      w->fcm->idle++;           /* signal that thread is idle */
      pthread_cond_signal(&w->fcm->cond_idle);
      while (!w->work)          /* wait until work is assigned */
        pthread_cond_wait(&w->fcm->cond_work, &w->fcm->mutex);
      pthread_mutex_unlock(&w->fcm->mutex);
      if (w->work < 0) break;   /* if processing was stopped, */
    }                           /* terminate the thread */
  #endif
    #ifdef FCM_BENCH            /* if to do some benchmarking */
    w->beg = timer();           /* note the start time of the thread */
    #endif
    if (w->work == RECTANGLE)   /* if to process a rectangle */
      SFXNAME(rec_rct)(w, w->ra, w->rb, w->ca, w->cb);
    else {                      /* if to process a triangle */
      while (1) {               /* process two strip parts */
        if (w->ca > w->ra)      /* process rectangle part of strip */
          SFXNAME(rec_rct)(w, w->ra, w->ca, w->ca, w->cb);
        SFXNAME(rec_trg)(w, w->ca, w->cb); /* process triangle part */
        if (w->ca > w->cm/2) break;  /* if second strip done, abort */
        a     = w->cm -w->cb;   /* get start of opposite strip */
        if (w->cb > a)       break;  /* if no opposite strip, abort */
        w->cb = w->cm -w->ca;   /* get start and end index */
        w->ca = a;              /* of the opposite strip */
      }
    }
    #ifdef FCM_BENCH            /* if to do some benchmarking */
    w->end = timer();           /* note the end time of the thread */
    #endif
  #ifndef _WIN32                /* not yet available for Windows */
    if (w->fcm->join) break;    /* if to join the threads */
    w->work = 0;                /* after the work is done, abort, */
  }                             /* otherwise mark work as done */
  #endif
  return THREAD_OK;             /* return a dummy result */
}  /* fill() */

/*--------------------------------------------------------------------------*/
#else                           /* --- split rectangle into 4 parts */

inline WORKERDEF(fill, p)
{                               /* --- fill (part of) the cache */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  DIM a, b, k;                  /* loop variables */

  assert(p);                    /* check the function argument */
  #ifndef _WIN32                /* not yet available for Windows */
  while (1) {                   /* work loop for blocked threads */
    if (!w->fcm->join) {        /* if to block and signal threads */
      pthread_mutex_lock(&w->fcm->mutex);
      w->fcm->idle++;           /* signal that thread is idle */
      pthread_cond_signal(&w->fcm->cond_idle);
      while (!w->work)          /* wait until work is assigned */
        pthread_cond_wait(&w->fcm->cond_work, &w->fcm->mutex);
      pthread_mutex_unlock(&w->fcm->mutex);
      if (w->work < 0) break;   /* if processing was stopped, */
    }                           /* terminate the thread */
  #endif
    #ifdef FCM_BENCH            /* if to do some benchmarking */
    w->beg = timer();           /* note the start time of the thread */
    #endif
    if (w->work == RECTANGLE) { /* if to process a rectangle */
      k = w->rb -w->ra;         /* get the size of the rectangles */
      for (a = w->ca; a < w->cb; a += k) {
        b = (a+k < w->cb) ? a+k : w->cb;
        SFXNAME(rec_rct)(w, w->ra, w->rb, a, b);
      } }                       /* split the strip into squares */
    else {                      /* if to process a triangle */
      while (1) {               /* process two strip parts */
        SFXNAME(rec_trg)(w, w->ca, w->cb);
        k = w->cb -w->ca;       /* compute leading triangle */
        for (a = w->ra; a < w->ca; a += k) {
          b = (a+k < w->ca) ? a+k : w->ca;
          SFXNAME(rec_rct)(w, a, b, w->ca, w->cb);
        }                       /* split the strip into squares */
        if (w->ca > w->cm/2) break; /* if second strip done, abort */
        a     = w->cm -w->cb;   /* get start of opposite strip */
        if (w->cb > a)       break; /* if no opposite strip, abort */
        w->cb = w->cm -w->ca;   /* get start and end index */
        w->ca = a;              /* of the opposite strip */
      }
    }
    #ifdef FCM_BENCH            /* if to do benchmarking */
    w->end = timer();           /* note the end time of the thread */
    #endif
  #ifndef _WIN32                /* not yet available for Windows */
    if (w->fcm->join) break;    /* if to join the threads */
    w->work = 0;                /* after the work is done, abort, */
  }                             /* otherwise mark work as done */
  #endif
  return THREAD_OK;             /* return a dummy result */
}  /* fill() */

#endif
/*--------------------------------------------------------------------------*/

inline int SFXNAME(fcm_fill) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- fill the cache */
  #ifdef RECTGRID               /* if to split rectangle into grid */
  DIM    i, j, n;               /* loop variables for threads */
  DIM    x, y;                  /* number of grid rows/columns */
  DIM    dx, dy;                /* number of voxels per grid row/col. */
  #else                         /* if to split rectangle into strips */
  DIM    i, n;                  /* loop variables for threads */
  #endif
  int    shape;                 /* shape of area to cache */
  DIM    k;                     /* number of rows/columns */
  DIM    m;                     /* reference column for mirroring */
  WORKER *worker;               /* worker for parallel execution */
  SFXNAME(WORK) *w;             /* data   for worker thread */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  DWORD  thid;                  /* dummy for storing the thread id */
  #endif                        /* (not really needed here) */

  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row > col) {              /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  fcm->ra = row -(row % fcm->tile);
  fcm->rb = fcm->ra +fcm->tile; /* compute the new row    range */
  if (fcm->rb > fcm->V) fcm->rb = fcm->V;
  fcm->ca = col -(col % fcm->tile);
  fcm->cb = fcm->ca +fcm->tile; /* compute the new column range */
  if (fcm->cb > fcm->V) fcm->cb = fcm->V;
  #ifdef FCM_BENCH              /* if to do some benchmarking */
  fcm->cnt += 1;                /* count the number of times */
  #endif                        /* that the cache was filled */
  w = fcm->work;                /* get the data for the workers */
  if (fcm->ra >= fcm->ca) {     /* if to process a triangle */
    shape = TRIANGLE;           /* note shape of area to cache */
    m =   fcm->cb +fcm->ca;     /* reference index for mirroring */
    k = ((fcm->cb -fcm->ca)/2 +(DIM)fcm->nthd-1) /(DIM)fcm->nthd;
    if (k <= 0) k = 1;          /* compute the number of voxels */
    for (n = 0; n < fcm->nthd; n++) {
      w[n].ra = fcm->ra;        /* traverse the threads */
      w[n].rb = fcm->rb;        /* note the row range and */
      w[n].ca = fcm->ca +n*k;   /* compute the column start index */
      if (w[n].ca >= m/2)       /* if beyond half of the columns, */
        break;                  /* all columns are provided for */
      w[n].cb = w[n].ca +k;     /* compute and store end index */
      if (w[n].cb >= m/2) w[n].cb = fcm->cb -n*k;
      w[n].cm = m;              /* ensure no duplicate computations */
    } }                         /* and note reference for mirroring */
  #ifdef RECTGRID               /* if to split rectangle into grid */
  else {                        /* if to process a rectangle */
    shape = RECTANGLE;          /* note shape of area to cache */
    x = fcm->gc;                /* get the number of grid  columns */
    k = fcm->cb -fcm->ca;       /* and the number of voxel columns */
    if (x > k) { x = k; dx = 1;}/* limit grid to voxel columns */
    else dx = (k +x-1) /x;      /* compute columns per grid column */
    y = fcm->gr;                /* get the number of grid  rows */
    k = fcm->rb -fcm->ra;       /* and the number of voxel rows */
    if (y > k) { y = k; dy = 1;}/* limit grid to voxel rows */
    else dy = (k +y-1) /y;      /* compute rows    per grid row */
    n  = 0;                     /* initialize the thread index */
    for (i = 0; i < x; i++) {   /* traverse the grid rows */
      for (j = 0; j < y; j++) { /* traverse the grid columns */
        w[n].ra = fcm->ra +j*dy;/* compute the row range */
        w[n].rb = w[n].ra +dy;  /* for the next grid cell */
        if (w[n].rb > fcm->rb) w[n].rb = fcm->rb;
        w[n].ca = fcm->ca +i*dx;/* compute the column range */
        w[n].cb = w[n].ca +dx;  /* for the next grid cell */
        if (w[n].cb > fcm->cb) w[n].cb = fcm->cb;
        n++;                    /* advance the thread index */
      }
    }                           /* initialize the grid cells */
  }                             /* to be processed by the threads */
  #else                         /* if to split rectangle into strips */
  else {                        /* if to process a rectangle */
    shape = RECTANGLE;          /* note shape of area to cache */
    k = ((fcm->rb-fcm->ra) +fcm->nthd-1) /fcm->nthd;
    if (k <= 0) k = 1;          /* compute the number of voxels */
    for (n = 0; n < fcm->nthd; n++) {
      w[n].ca = fcm->ca;        /* traverse the threads */
      w[n].cb = fcm->cb;        /* note the column range and */
      w[n].ra = fcm->ra +n*k;   /* compute the row start index */
      if (w[n].ra >= fcm->rb)   /* if beyond the last row, */
        break;                  /* all rows are provided for */
      w[n].rb = w[n].ra +k;     /* compute and store end index */
      if (w[n].rb >  fcm->rb) w[n].rb = fcm->rb;
    }                           /* initialize the row/column ranges */
  }                             /* to be processed by the threads */
  #endif
  fcm->err = 0;                 /* clear the error status */
  worker   = SFXNAME(fill);     /* get the worker function */
  if (n <= 1) {                 /* if there is only one thread, */
    w[0].work = shape;          /* note shape of area to cache and */
    worker(w); return 0;        /* execute the worker directly */
  }
  #ifdef _WIN32                 /* if Microsoft Windows system */
  for (i = 0; i < n; i++) {     /* traverse the threads */
    w[i].work = shape;          /* set the area shape identifier */
    fcm->threads[i] = CreateThread(NULL, 0, worker, w+i, 0, &thid);
    if (!fcm->threads[i]) { fcm->err = -1; break; }
  }                             /* create a thread for each strip */
  WaitForMultipleObjects(n, fcm->threads, TRUE, INFINITE);
  for (i = 0; i < n; i++)          /* wait for threads to finish, */
    CloseHandle(fcm->threads[i]);  /* then close all thread handles */
  #else                         /* if Linux/Unix system */
  if (fcm->join) {              /* if to create new threads */
    for (i = 0; i < n; i++) {   /* traverse the threads */
      w[i].work = shape;        /* set the area shape identifier */
      if (pthread_create(fcm->threads+i, NULL, worker, w+i) != 0) {
        fcm->err = -1; break; } /* create a thread for each part */
    }                           /* of the cache to fill */
    for (i = 0; i < n; i++) {   /* wait for threads to finish */
      pthread_join(fcm->threads[i], NULL);
    } }                         /* (join threads with this one) */
  else {                        /* if to block and signal threads */
    pthread_mutex_lock(&fcm->mutex);
    for (i = 0; i < n; i++)     /* traverse the treads and */
      w[i].work = shape;        /* assign work to each thread */
    fcm->idle = 0;              /* signal that work was assigned */
    pthread_cond_broadcast(&fcm->cond_work);
    while (fcm->idle < n)       /* wait for threads to finish */
      pthread_cond_wait(&fcm->cond_idle, &fcm->mutex);
    pthread_mutex_unlock(&fcm->mutex);
  }                             /* now all threads are blocked */
  #endif
  #ifdef FCM_BENCH              /* if to do some benchmarking */
  double min = +INFINITY;       /* initialize minimal start time */
  double max = -INFINITY;       /* and maximal end time of a thread */
  for (i = 0; i < n; i++) {     /* traverse the threads */
    if (w[i].beg < min) min = w[i].beg;
    if (w[i].end > max) max = w[i].end;
  }                             /* find min. start/max. end time */
  fcm->sum += max -min;         /* find time span for threads */
  for (i = 0; i < n; i++) {     /* traverse the threads again */
    fcm->beg += w[i].beg -min;  /* compute loss at start */
    fcm->end += max -w[i].end;  /* compute loss at end */
  }                             /* of each thread */
  #endif
  return fcm->err;              /* return the error status */
}  /* fcm_fill() */


/*----------------------------------------------------------------------------
  Retrieval Functions
----------------------------------------------------------------------------*/

inline REAL SFXNAME(fcm_cache) (SFXNAME(FCMAT) *fcm, DIM row, DIM col)
{                               /* --- get cached correlation coeff. */
  assert(fcm                    /* check the function arguments */
  &&    (row >= 0) && (row < fcm->V) && (col >= 0) && (col < fcm->V));
  if (row == col)               /* if diagonal element, */
    return fcm->diag;           /* return a fixed value */
  if (row > col) {              /* ensure col >= row (upper triangle) */
    DIM t = row; row = col; col = t; }
  #ifdef SAFETHREAD             /* if to guarantee fill success */
  if (((row < fcm->ra) || (row >= fcm->rb)
  ||   (col < fcm->ca) || (col >= fcm->cb))
  &&  (SFXNAME(fcm_fill)(fcm, row, col) != 0)) {
    DIM n = fcm->nthd;          /* if outside cache, fill cache; */
    fcm->nthd = 1;              /* on failure retry with one thread */
    fprintf(stderr, "safethread\n");
    SFXNAME(fcm_fill)(fcm, row, col);
    fcm->nthd = n;              /* (a single thread cannot fail) */
  }                             /* and restore the number of threads */
  #else                         /* if to allow a thread failure */
  if (((row < fcm->ra) || (row >= fcm->rb)
  ||   (col < fcm->ca) || (col >= fcm->cb))
  &&  (SFXNAME(fcm_fill)(fcm, row, col) != 0))
    return -INFINITY;           /* if outside cache, fill cache */
  #endif
  row -= fcm->ra;               /* compute row and column in cache */
  col -= fcm->ca;               /* and retrieve correlation coeff. */
  return fcm->cache[(size_t)row *(size_t)fcm->tile +(size_t)col];
}  /* fcm_cache() */

/*----------------------------------------------------------------------------
  Recursion Handling
----------------------------------------------------------------------------*/
#if   _FCM_PASS == 1            /* if in first of two passes */
#undef REAL
#undef SUFFIX
#include "fcmat2.h"             /* process file recursively */
#elif _FCM_PASS == 2
#undef REAL
#endif

#undef SUFFIX
#undef SFXNAME
#undef SFXNAME_1
#undef SFXNAME_2

#undef  _FCM_PASS

#define FCMAT2_H
#endif  /* #ifndef FCMAT2_H */
