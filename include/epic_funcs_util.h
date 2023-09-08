/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * All material in this file is covered by the                     *
 * GNU General Public License.                                     *
 *                                                                 *
 * Copyright (C) 1998-2009 Timothy E. Dowling, except as noted.    *
 *                                                                 *
 * The fcmp() function is derived from the fcmp() function that is *
 * Copyright (c) 1998-2000 Theodore C. Belding                     *
 * University of Michigan Center for the Study of Complex Systems. *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License     *
 * as published by the Free Software Foundation; either version 2  *
 * of the License, or (at your option) any later version.          *
 * A copy of this License is in the file:                          *
 *   $EPIC_PATH/License.txt                                        *
 *                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
 *                                                                 *
 * You should have received a copy of the GNU General Public       *
 * License along with this program; if not, write to the Free      *
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     *
 * Boston, MA 02110-1301, USA.                                     *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef EPIC_FUNCS_UTIL_H
#define EPIC_FUNCS_UTIL_H

/* * * * * * * * * epic_funcs_util.h * * * * * * * * * * * * * * * * * * 
 *                                                                     *
 * Header file for utility functions that are essentially independent  *
 * of the EPIC model; the source code is in                            *
 *   $EPIC_PATH/src/shared/epic_funcs_util.c.                          *
 *                                                                     *
 * NOTE: One dependency with the EPIC model is the floating-point      *
 *       precision, EPIC_PRECISION, which is an environment variable.  *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h> 
#include <stdlib.h>   
#include <errno.h>    
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#include "epic_datatypes.h"
#define FLOAT EPIC_FLOAT

/*
 * Parameters.
 */
#define N_STR 256

/*
 * Logical:
 */
#define TRUE    1
#define FALSE   0
#define SUCCESS TRUE
#define FAILURE FALSE

#define SETUP_UTIL   0
#define RUN_UTIL     1
#define CLEANUP_UTIL 2

/*
 * Mathematical:
 */
#if !defined(M_E)
#  define M_E         2.7182818284590452354
#  define M_LOG2E     1.4426950408889634074
#  define M_LOG10E    0.43429448190325182765
#  define M_LN2       0.69314718055994530942
#  define M_LN10      2.30258509299404568402
#  define M_PI        3.14159265358979323846
#  define M_PI_2      1.57079632679489661923
#  define M_PI_4      0.78539816339744830962
#  define M_1_PI      0.31830988618379067154
#  define M_2_PI      0.63661977236758134308
#  define M_2_SQRTPI  1.12837916709551257390
#  define M_SQRT2     1.41421356237309504880
#  define M_SQRT1_2   0.70710678118654752440
#endif

#define DEG (M_PI/180.)

#undef MIN
#define MIN(x,y) ({ \
         const FLOAT _x = (FLOAT)(x); \
         const FLOAT _y = (FLOAT)(y); \
         _x < _y ? _x : _y; })

#undef MAX
#define MAX(x,y) ({ \
         const FLOAT _x = (FLOAT)(x); \
         const FLOAT _y = (FLOAT)(y); \
         _x > _y ? _x : _y; })

#undef LIMIT_RANGE
#define LIMIT_RANGE(min,x,max) ({ \
         const FLOAT _min = (FLOAT)(min); \
         const FLOAT _x   = (FLOAT)(x);   \
         const FLOAT _max = (FLOAT)(max); \
         _x < _min ? _min : ( _x > _max ? _max : _x ); })

#undef IMIN
#define IMIN(i,j) ({ \
         const int _i = (int)(i); \
         const int _j = (int)(j); \
         _i < _j ? _i : _j; })

#undef IMAX
#define IMAX(i,j) ({ \
         const int _i = (int)(i); \
         const int _j = (int)(j); \
         _i > _j ? _i : _j; })

#undef NINT
#define NINT(x) ({ \
         const FLOAT _x = (FLOAT)(x); \
         _x > 0. ? (int)(_x+.5) : (int)(_x-.5); })

#undef SIGN
#define SIGN(x) ({ \
         const FLOAT _x = (FLOAT)(x); \
         _x == 0. ? 0. : (_x > 0. ? 1. : -1.); })

#undef NR_SIGN
#define NR_SIGN(a,b) ((b) > 0. ? fabs(a) : -fabs(a))

#undef SQR
#define SQR(x) ({ \
          const FLOAT _x = (FLOAT)(x); \
          _x*_x; })

/* 
 * DS stands for Dennis and Schnabel (1996).
 */
#define DS_EPS      machine_epsilon()
#define DS_SQRT_EPS sqrt(DS_EPS)
#define DS_TOL_X    pow(DS_EPS,2./3.)
#define DS_TOL_F    pow(DS_EPS,1./3.)
#define DS_TOL_MIN  pow(DS_EPS,2./3.)
#define DS_ALPHA    1.e-4
#define DS_MAX_STEP 100.

#define DS_LINE_STEP   0
#define DS_HOOK_STEP   1    /* Not yet implemented. */
#define DS_DOGLEG_STEP 2
/*
 * Values for status (retcode in DS96):
 */
#define DS_X_ACCEPTED        0
#define DS_X_NO_PROGRESS     1
#define DS_REDUCE_DELTA      2
#define DS_INCREASE_DELTA    3
#define DS_MAX_IT_EXCEEDED   4
#define DS_MAX_TAKEN_5       5
#define DS_INITIAL           6
#define DS_SINGULAR_JACOBIAN 7

/*
 * Shift macros:
 */
#define QT(i,j)       qt[j+(i)*n]
#define R(i,j)         r[j+(i)*n]
#define A_EIG(r,c) a_eig[(c)+klen*(r)]
#define AA_EIG(r,c) aa_eig[(c)+klen*(r)]

/* * * * * * * * * * * * *
 *                       *
 * Function prototypes.  *
 *                       *
 * * * * * * * * * * * * */

/* 
 * Memory allocation. 
 */
int *ivector(int   nl,
             int   nh,
             char *calling_func);

void free_ivector(int  *m,
                  int   nl,
                  int   nh,
                  char *calling_func);

FLOAT *fvector(int   nl,
               int   nh,
               char *calling_func);

void free_fvector(FLOAT *m,
                  int    nl,
                  int    nh,
                  char  *calling_func);

double *dvector(int   nl,
                int   nh,
                char *calling_func);

void free_dvector(double *m,
                  int     nl,
                  int     nh,
                  char   *calling_func);

float_triplet *ftriplet(int   nl,
                        int   nh,
                        char *calling_func);

void free_ftriplet(float_triplet *m,
                   int            nl,
                   int            nh,
                   char          *calling_func);

/*
 *  The following functions are adapted from 
 *  Numerical Recipes in C:
 */
void spline(int            n,
            float_triplet *table,
            FLOAT          yp0, 
            FLOAT          ypn);

FLOAT splint(register FLOAT          xx, 
             register float_triplet *table,
             register FLOAT          dx);

FLOAT linint(register FLOAT          xx,
             register float_triplet *table,
             register FLOAT          dx);

void spline_pchip(int            n,
                  float_triplet *table);

FLOAT splint_pchip(FLOAT          xx,
                   float_triplet *table,
                   FLOAT          h);

void periodic_spline_pchip(int            n,
                           float_triplet *table);

FLOAT pchst(FLOAT arg1,
            FLOAT arg2);

FLOAT lagrange_interp(register FLOAT *f,
                      register FLOAT *x,
                      register int order);

FLOAT gamma_nr(FLOAT  xx);

FLOAT sech2(FLOAT);

void exp_integral_setup(float_triplet *exp3table,
                        float_triplet *exp4table,
                        int            numdatapoints);

FLOAT normed_legendre(int   l,
                      int   m,
                      FLOAT x);

FLOAT machine_epsilon(void);

int find_root(FLOAT  x1,
              FLOAT  x2,
              FLOAT  xacc,
              FLOAT *x_root,
              FLOAT  (*func)(FLOAT));

int broyden_root(int    n,
                 FLOAT *x,
                 void  (*vecfunc)(int,FLOAT *,FLOAT *),
                 FLOAT tol_f,
                 int   max_it);

int global_step(int     n,
                FLOAT *x_old,
                FLOAT  f_old,
                FLOAT *g,
                FLOAT *r,
                FLOAT *sn,
                FLOAT  max_step,
                FLOAT *delta,
                int     step_type,
                int    *status,
                FLOAT *x,
                FLOAT *f,
                FLOAT *fvec,
                void   (*vecfunc)(int,FLOAT *,FLOAT *));

int line_search(int    n,
                FLOAT *x_old,
                FLOAT  f_old,
                FLOAT *g,
                FLOAT *sn,
                FLOAT  max_step,
                int   *status,
                FLOAT *x,
                FLOAT *f,
                FLOAT *fvec,
                void  (*vecfunc)(int,FLOAT *,FLOAT *));

int dogleg_driver(int    n,
                  FLOAT *x_old,
                  FLOAT  f_old,
                  FLOAT *g,
                  FLOAT *r,
                  FLOAT *sn,
                  FLOAT  max_step,
                  FLOAT *delta,
                  int   *status,
                  FLOAT *x,
                  FLOAT *f,
                  FLOAT *fvec,
                  void  (*vecfunc)(int,FLOAT *,FLOAT *));

int dogleg_step(int    n,
                FLOAT *g,
                FLOAT *r,
                FLOAT *sn,
                FLOAT  newt_length,
                FLOAT  max_step,
                FLOAT *delta,
                int   *first_dog,
                FLOAT *s_hat,
                FLOAT *nu_hat,
                FLOAT *s);

int trust_region(int    n,
                 FLOAT *x_old,
                 FLOAT  f_old,
                 FLOAT *g,
                 FLOAT *s,
                 int    newt_taken,
                 FLOAT  max_step,
                 int    step_type,
                 FLOAT *r,
                 FLOAT *delta,
                 int   *status,
                 FLOAT *x_prev,
                 FLOAT *f_prev,
                 FLOAT *x,
                 FLOAT *f,
                 FLOAT *fvec,
                 void  (*vecfunc)(int,FLOAT *,FLOAT *));

int qr_decompose(int    n,
                 FLOAT *r,
                 FLOAT *c,
                 FLOAT *d);

void qr_update(int    n,
               FLOAT *r,
               FLOAT *qt,
               FLOAT *u,
               FLOAT *v);

void qr_rotate(int    n,
               FLOAT *r,
               FLOAT *qt,
               int    i,
               FLOAT  a,
               FLOAT  b);

void lu_decompose(int   n,
                 FLOAT *a,
                 int   *index,
                 FLOAT *d);

void lu_backsub(int    n,
                FLOAT *a,
                int   *index,
                FLOAT *b);

void lu_improve(int    n,
                FLOAT *a,
                FLOAT *alu,
                int   *index,
                FLOAT *b,
                FLOAT *x);

int find_place_in_table(int            n,
                        float_triplet *table,
                        FLOAT          x,
                        FLOAT         *dx);

int hunt_place_in_table(int            n,
                        float_triplet *table,
                        FLOAT          x,
                        FLOAT         *dx,
                        int            il);

/*
 * Choices for pivot_type in tridiag():
 */
#define WITHOUT_PIVOTING 0
#define WITH_PIVOTING    1

void tridiag(int    n,
             FLOAT *a,
             FLOAT *b,
             FLOAT *c,
             FLOAT *r,
             FLOAT *u,
             int    pivot_type);

void band_decomp(int     n,
                 int     m1,
                 int     m2,
                 FLOAT  *a,
                 FLOAT  *al,
                 int    *index,
                 FLOAT  *d);

void band_back_sub(int     n,
                   int     m1,
                   int     m2,
                   FLOAT  *a,
                   FLOAT  *al,
                   int    *index,
                   FLOAT  *b);

void band_multiply(int     n,
                   int     m1,
                   int     m2,
                   FLOAT  *a,
                   FLOAT  *x,
                   FLOAT  *b);

void band_improve(int     n,
                  int     m1,
                  int     m2,
                  FLOAT  *aorig,
                  FLOAT  *a,
                  FLOAT  *al,
                  int    *index,
                  FLOAT  *b,
                  FLOAT  *x);

FLOAT poly_interp(int     n,
                  FLOAT  *xa,
                  FLOAT  *ya,
                  FLOAT   x,
                  FLOAT  *dy);

void poly_coeff(int    n,
                FLOAT *x,
                FLOAT *y,
                FLOAT *coeff);

FLOAT nth_trapezoidal(int   n,
                      FLOAT (*func)(FLOAT),
                      FLOAT a,
                      FLOAT b);

FLOAT romberg_integral(FLOAT (*func)(FLOAT),
                       FLOAT a,
                       FLOAT b,
                       FLOAT tol);

void crank_nicolson(int    n,
                    FLOAT  dt,
                    FLOAT *z,
                    FLOAT *A,
                    FLOAT *mu,
                    FLOAT *rho,
                    FLOAT *ANS);

void compact_differentiation(int    action,
                             int    n,
                             FLOAT *z,
                             FLOAT *f,
                             FLOAT *dfdz);

void compact_integration(int    n,
                         FLOAT *z,
                         FLOAT *dfdz,
                         FLOAT *f);

void hqr(FLOAT *a,
         int    n,
         FLOAT *wr,
         FLOAT *wi);

void quicksort(FLOAT *mag, 
               int    left, 
               int    right);

void swap(FLOAT *mag,
          int    i,
          int    j);

void four1(FLOAT          *,
           unsigned long, 
           int              );

void realft(FLOAT          *,
            unsigned long,
            int             );

/*
 * NOTE: For LINUX, with -D_BSD_SOURCE cabs() is defined, such that
 * a type-mismatch error occurs if we name the function c_abs() as cabs().
 */
complex c_num(FLOAT x, FLOAT y);
complex c_mult(complex z1,complex z2);
complex c_add(complex z1,complex z2);
complex c_sub(complex z1,complex z2);
complex c_exp(complex z);
FLOAT   c_abs(complex z);
FLOAT   c_real(complex z);
FLOAT   c_imag(complex z);

int fcmp(double x1,
         double x2);

void least_squares(FLOAT *x,
                   FLOAT *y,
                   int    n,
                   FLOAT *a);

void savitzky_golay(FLOAT *c,
                    int    np,
                    int    nl,
                    int    nr,
                    int    ld,
                    int    m);

FLOAT random_number(long *idum);

FLOAT lat_centric_to_graphic(FLOAT lat,
                             FLOAT rerp);

FLOAT lat_graphic_to_centric(FLOAT lat,
                             FLOAT rerp);

FLOAT surface_area_oblate(FLOAT a,
                          FLOAT c);

void util_error(char *calling_function,
                char *Message);

/* * * * * * * * * *  end of epic.h  * * * * * * * * * * * * * * * * * * * * */ 
#endif

