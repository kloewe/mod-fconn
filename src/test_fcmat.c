/*----------------------------------------------------------------------
  File    : test_fcmat.c
  Contents: tests
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
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

#include "cpuinfo.h"
#include "binarize.h"
#include "pcc.h"
#include "tetracc.h"
#include "fcmat.h"

/*----------------------------------------------------------------------
  Preprocessor definitions
----------------------------------------------------------------------*/
#define PRGNAME     "test_fcmat"
#define DESCRIPTION "test functional connectivity matrix implementation"
#define VERSION     "v20160819         " \
                    "(c) 2014-2016   Kristian Loewe/Christian Borgelt"

/*--------------------------------------------------------------------*/
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

/*--------------------------------------------------------------------*/
#define INDEX(i,j,N)    ((size_t)(i)*((size_t)(N)+(size_t)(N) \
                        -(size_t)(i)-3)/2-1+(size_t)(j))

/*--------------------------------------------------------------------*/
#define E_NONE         0        /* no error */
#define E_NOMEM      (-1)       /* not enough memory */
#define E_FOPEN      (-2)       /* cannot open file */
#define E_FREAD      (-3)       /* read error on file */
#define E_FWRITE     (-4)       /* write error on file */
#define E_STDIN      (-5)       /* double assign. of standard input */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* wrong number of arguments */
#define E_THREAD     (-9)       /* thread computation error */

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
static const char *prgname;     /* program name for error messages */

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
static const char *errmsgs[] = {
  /* E_NONE      0 */  "no error",
  /* E_NOMEM    -1 */  "not enough memory",
  /* E_FOPEN    -2 */  "cannot open file %s",
  /* E_FREAD    -3 */  "read error on file %s",
  /* E_FWRITE   -4 */  "write error on file %s",
  /* E_STDIN    -5 */  "double assignment of standard input",
  /* E_OPTION   -6 */  "unknown option -%c",
  /* E_OPTARG   -7 */  "missing option argument",
  /* E_ARGCNT   -8 */  "wrong number of arguments",
  /* E_THREAD   -9 */  "thread computation error",
  /*           -10 */  "unknown error"
};                              /* list of error messages */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
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

/*--------------------------------------------------------------------*/

static int error (int code, ...)
{                               /* --- print an error message */
  int        k;                 /* maximal error code */
  va_list    args;              /* list of variable arguments */
  const char *msg;              /* error message */

  assert(prgname);              /* check the program name */
  va_start(args, code);         /* start variable arguments */
  if      (code > 0) {          /* if an error message is given, */
    msg = va_arg(args, const char*);        /* print it directly */
    if (msg) fprintf(stderr, "\n%s: %s\n", prgname, msg); }
  else if (code < 0) {          /* if code and arguments are given */
    k = 1-(int)(sizeof(errmsgs)/sizeof(*errmsgs));
    if (code < k) code = k;     /* check and adapt the error code */
    msg = errmsgs[-code];       /* get the error message format */
    if (!msg) msg = errmsgs[-k];/* check and adapt the message */
    fprintf(stderr, "\n%s: ", prgname);
    vfprintf(stderr, msg, args);/* print the error message and */
    fputc('\n', stderr);        /* terminate the output line */
  }
  va_end(args);                 /* end variable arguments */
  exit(abs(code));              /* abort the program */
}  /* error() */

/*----------------------------------------------------------------------
  Main Function
----------------------------------------------------------------------*/
int main (int argc, char* argv[])
{                               /* --- main function for testing */
  int     k = 0;                /* counter */
  char    *s;                   /* to traverse the options */
  char    **optarg = NULL;      /* option argument */
  DIM     V     = 0;            /* number of voxels */
  DIM     T     = 0;            /* number of time points */
  DIM     C     = 8192;         /* tile size for caching */
  int     P     = proccnt();    /* number of threads */
  int     mode                  /* computation mode */
            = FCM_PCC|FCM_THREAD|FCM_CACHE;
  int     valid = 0;            /* flag for result validation */
  long    S     = time(NULL);   /* seed value for random numbers */
  size_t  E     = 0;            /* number of edges of the graph */
  REAL    *data;                /* data array */
  REAL    *corr;                /* correlation coefficients */
  int     t;                    /* indicator for another element */
  DIM     r, c;                 /* loop variables (row and column) */
  REAL    a, b;                 /* to compare correlation coeffs. */
  int     diff;                 /* indicator for a difference */
  int     n = 1;                /* number of repetitions */
  double  t0;                   /* timer for measurements */
  FCMAT   *fcm;                 /* functional connectivity matrix */

  prgname = argv[0];            /* get program name for error msgs. */

  /* --- print usage message --- */
  if (argc > 1) {               /* if arguments are given */
    fprintf(stderr, "%s - %s\n", PRGNAME, DESCRIPTION);
    fprintf(stderr, VERSION); } /* print a startup message */
  else {                        /* if no argument is given */
    printf("usage: %s [options] V T\n", argv[0]);
    printf("%s\n", DESCRIPTION);
    printf("%s\n", VERSION);
    printf("-v       validate result (compare to pcc)         "
           "(default: performance)\n");
    printf("-s#      seed value for random number generator   "
           "(default: time)\n");
    printf("-t#      number of parallel threads               "
           "(default: %d)\n", P);
    printf("-c#      size of cache (number of rows/columns)   "
           "(default: 0)\n");
    printf("-j       join and re-create threads               "
           "(default: blk+sig)\n");
    printf("V        number of voxels\n");
    printf("T        number of time points\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */

  /* --- evaluate arguments --- */
  for (int i = 1; i < argc; i++) {  /* traverse the arguments */
    s = argv[i];                /* get an option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse the options */
        switch (*s++) {         /* evaluate the options */
          case 'v': valid  = 1;                     break;
          case 's': S      =      strtol(s, &s, 0); break;
          case 't': P      = (int)strtol(s, &s, 0); break;
          case 'c': C      =      strtodim(s, &s);  break;
          case 'j': mode  |= FCM_JOIN;              break;
          #if 0
          case 'x': optarg = &arg;                  break;
          #endif
          default : error(E_OPTION, *--s);          break;
        }                       /* set the option variables */
        if (optarg && *s) { *optarg = s; optarg = NULL; break; }
      } }                       /* get an option argument */
    else {                      /* -- if argument is no option */
      switch (k++) {            /* evaluate non-options */
        case  0: V = strtodim(s, NULL); break;
        case  1: T = strtodim(s, NULL); break;
        default: error(E_ARGCNT);       break;
      }                         /* note number of voxels */
    }                           /* and  number of time points */
  }
  if (optarg) error(E_OPTARG);  /* check option arguments */
  if (k != 2) error(E_ARGCNT);  /* and number of arguments */
  if (P <  0) error(-10);       /* and the number of threads */
  if (C <  0) error(-10);       /* get the tile size for caching */
  if (S <  0) error(-10);       /* get the seed value and */
  srand((unsigned)S);           /* seed the random number generator */
  E = (size_t)V*(size_t)(V-1)/2;/* compute the number of edges */
  fprintf(stderr, "\n");        /* terminate the startup message */

  /* --- generate test data --- */
  data = malloc((size_t)V *(size_t)T *sizeof(REAL));
  if (!data) error(E_NOMEM);    /* allocate and generate data */
  for (size_t i = 0; i < (size_t)(T*V); i++)
    data[i] = (REAL)(rand()/((double)RAND_MAX+1));

  /* --- validate computation results --- */
  if (valid) {                  /* if to validate the results */
    fprintf(stderr, "computing reference result using pcc ... ");
    corr = malloc((size_t)V *(size_t)(V-1)/2 *sizeof(REAL));
    if (!corr) error(E_NOMEM);  /* allocate memory for corr. coeffs. */
    pccx(data, corr, (int)V, (int)T, PCC_AUTO);
    fprintf(stderr, "done.\n"); /* compute correlation coefficients */

    fprintf(stderr, "test (fcm_get) ... ");
    fcm = fcm_create(data, V, T, mode, P, C);
    if (!fcm) error(E_NOMEM);   /* create functional connect. matrix */
    diff = 0;                   /* initialize the difference counter */
    for (DIM i = 0; i < V; i++) { /* traverse rows and cols */
      for (DIM j = i+1; j < V; j++) {
        a = fcm_get(fcm,i,j);
        b = corr[INDEX(i,j,V)]; /* get correlation coefficients */
        if (a == b) continue;   /* and compare them */
        if (!diff) fprintf(stderr, "\n");
        fprintf(stderr, "%6"DIM_FMT" %6"DIM_FMT, i, j);
        fprintf(stderr, ": % 18.16f % 18.16f\n", a, b);
        diff += 1;              /* print any difference and */
      }                         /* count the number of differences */
    }
    fcm_delete(fcm);   /* delete the func. connect. matrix */
    if (diff) fprintf(stderr, "failed [%d].\n", diff);
    else      fprintf(stderr, "passed.\n");

    fprintf(stderr, "test (fcm_next) ... ");
    fcm = fcm_create(data, V, T, mode, P, C);
    if (!fcm) error(E_NOMEM); /* create functional connect. matrix */
    diff = 0;                 /* initialize the difference flag */
    for (t = fcm_first(fcm); t > 0; t = fcm_next(fcm)) {
      r = fcm_row(fcm); /* traverse the matrix elements */
      c = fcm_col(fcm); /* and get their row and column */
      a = fcm_value(fcm);
      b = corr[INDEX(r,c,V)]; /* get correlation coefficients */
      if (a == b) continue;   /* and compare them */
      if (!diff) fprintf(stderr, "\n");
      fprintf(stderr, "%6"DIM_FMT" %6"DIM_FMT, r, c);
      fprintf(stderr, ": % 18.16f % 18.16f\n", a, b);
      diff += 1;              /* print any difference and */
    }                         /* count the number of differences */
    if (t < 0) error(E_THREAD);  /* check for a computation error */
    fcm_delete(fcm); /* delete the func. connect. matrix */
    if (diff) fprintf(stderr, "failed [%d].\n", diff);
    else      fprintf(stderr, "passed.\n");

    free(corr);                 /* delete correlation coefficients */
  }

  /* --- test performance --- */
  else {                        /* if to test performance */
    fprintf(stderr, "perf (fcm_get) ... ");
    t0 = timer();               /* start the timer */
    fcm = fcm_create(data, V, T, mode, P, C);
    if (!fcm) error(E_NOMEM);   /* create functional connect. matrix */
    for (DIM i = 0; i < fcm_dim(fcm); i++)
      for (DIM j = i+1; j < fcm_dim(fcm); j++)
        a = fcm_get(fcm,i,j);
    fcm_delete(fcm);   /* delete the func. connect. matrix */
    t0 = (timer()-t0)/(double)n;
    fprintf(stderr, "done.\n");
    fprintf(stderr, "time:   %8.2fs\n", t0);
    fprintf(stderr, "Mccf/s: %8.2f\n", (double)E/t0/1e6);

    fprintf(stderr, "perf (fcm_next) ... ");
    t0 = timer();               /* start the timer */
    fcm = fcm_create(data, V, T, mode, P, C);
    if (!fcm) error(E_NOMEM); /* create functional connect. matrix */
    for (t = fcm_first(fcm); t > 0; t = fcm_next(fcm)) {
      r = fcm_row(fcm);  /* traverse the matrix elements */
      c = fcm_col(fcm);  /* and get their row and column */
      a = fcm_value(fcm);
    }                         /* get the matrix elements */
    if (t < 0) error(E_THREAD);  /* check for a computation error */
    fcm_delete(fcm); /* delete the func. connect. matrix */
    t0 = (timer()-t0)/(double)n;
    fprintf(stderr, "done.\n");
    fprintf(stderr, "time:   %8.2fs\n", t0);
    fprintf(stderr, "Mccf/s: %8.2f\n", (double)E/t0/1e6);
  }

  free(data);                   /* delete the generated data */
  return 0;                     /* return 'ok' */
}  /* main() */
