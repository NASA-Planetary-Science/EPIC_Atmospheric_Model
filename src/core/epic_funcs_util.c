/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * All material in this file is covered by the                     *
 * GNU General Public License.                                     *
 *                                                                 *
 * Copyright (C) 1998-2023 Timothy E. Dowling, except as noted.    *
 *                                                                 *
 * The fcmp() function is derived from the fcmp() function that is *
 * Copyright (c) 1998-2000 Theodore C. Belding                     *
 * University of Michigan Center for the Study of Complex Systems. *
 *                                                                 *
 * The monotonic spline, spline_pchip(), is derived from the PCHIP *
 * (Piecewise Cubic Hermite Interpolating Polynomial) code written *
 * by Fred Fritsch, which was tranlated into C by John Burkardt,   *
 * and adapted for use here by Tim Dowling.                        *
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

/* * * * * * * * *  epic_funcs_util.c  * * * * * * * * * * * * * * * * * * * * 
 *                                                                           *
 *       These functions for math and utility operations should  not         *
 *       reference EPIC variables by name or index, such that they can be    *
 *       used outside the model.                                             *
 *                                                                           *
 *       NOTE: One dependency with the EPIC model is the environment         *
 *             variable EPIC_PRECISION, which specifies the floating-point   *
 *             precision.                                                    *
 *                                                                           *
 *       This file includes the following:                                   *
 *                                                                           *
 *           ivector(),free_ivector()                                        *
 *           fvector(),free_fvector()                                        *
 *           dvector(),free_dvector()                                        *
 *           ftriplet(),free_ftriplet()                                      *
 *           spline(),splint()                                               *
 *           linint()                                                        *
 *           spline_pchip(), splint_pchip()                                  *
 *           pchst()                                                         *
 *           lagrange_interp()                                               *
 *           gamma_nr()                                                      *
 *           sech2()                                                         *
 *           exp_integral(), exp_integral_setup()                            *
 *           machine_epsilon()                                               *
 *           find_root()                                                     *
 *           broyden_root()                                                  *
 *           global_step()                                                   *
 *           line_search()                                                   *
 *           dogleg_driver()                                                 *
 *           dogleg_step()                                                   *
 *           trust_region()                                                  *
 *           qr_decompose()                                                  *
 *           qr_update()                                                     *
 *           qr_rotate()                                                     *
 *           lu_decompose()                                                  *
 *           lu_backsub()                                                    *
 *           lu_improve()                                                    *
 *           find_place_in_table()                                           *
 *           hunt_place_in_table()                                           *
 *           tridiag()                                                       *
 *           band_decomp()                                                   *
 *           band_back_sub()                                                 *
 *           band_multiply()                                                 *
 *           band_improve()                                                  *
 *           poly_interp()                                                   *
 *           poly_coeff()                                                    *
 *           nth_trapezoidal()                                               *
 *           romberg_integral()                                              *
 *           compact_integration()                                           *
 *           compact_differentiation()                                       *
 *           crank_nicolson()                                                *
 *           hqr()                                                           *
 *           quicksort()                                                     *
 *           swap()                                                          *
 *           four1(),realft()                                                *
 *           c_num(),c_mult(),c_add(),c_sub()                                *
 *           c_exp(),c_abs(),c_real(),c_imag                                 *
 *           fcmp()                                                          *
 *           least_squares()                                                 *
 *           savitzky_golay()                                                *
 *           random_number()                                                 *
 *           lat_centric_to_graphic(), lat_graphic_to_centric()              *
 *           surface_area_oblate()                                           *
 *           util_error()                                                    *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic_funcs_util.h>

/*
 * Global declarations.
 */
char
  Message[256];

/*
 * The following macro is useful to insert into code while trying to corral a problem.
 */
#undef  DEBUG_MILESTONE
#define DEBUG_MILESTONE(comment) fprintf(stderr,"%s, %2d: "#comment"\n", \
                                         dbmsname,++idbms);fflush(stderr);

/*======================= ivector() ============================================*/
      
/*
 *  Allocates memory for a 1D int array 
 *  with range [nl..nh].
 */

#undef  DEBUG

int *ivector(int  nl, 
             int  nh,
             char *calling_func)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  int         
    *m;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="ivector";

  if (nh < nl) {
    sprintf(Message,"called by %s, range (%d,%d)",calling_func,nl,nh);
    util_error(dbmsname,Message);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (int *)calloc(len_safe,sizeof(int));
  if (!m) {
    sprintf(Message,"called by %s, nl=%d,nh=%d,len_safe=%d",calling_func,nl,nh,len_safe);
    util_error(dbmsname,Message);
  }
  m -= nl_safe;

#if defined(DEBUG)
  fprintf(stderr,"ivector() called by %s \n",calling_func);
#endif

  return m;
}

/*======================= end of ivector() ====================================*/

/*======================= free_ivector() ======================================*/

#undef  DEBUG

/*
 *  Frees memory allocated by ivector().
 */

void free_ivector(int  *m, 
                  int  nl, 
                  int  nh,
                  char *calling_func)
{
  int  
    nl_safe;

  if (m) {
    nl_safe = (nl < 0) ? nl : 0;
    m += nl_safe;
    free(m);
  }

#if defined(DEBUG)
  if (m) {
    fprintf(stderr,"free_ivector() called by %s \n",calling_func);
  }
  else {
    fprintf(stderr,"free_ivector() called by %s for NULL pointer\n",calling_func);
  }
#endif

  return;
}

/*======================= end of free_ivector() ================================*/

/*======================= fvector() ============================================*/
      
/*
 *  Allocates memory for a 1D FLOAT array 
 *  with range [nl..nh].
 */

#undef  DEBUG 

FLOAT *fvector(int  nl, 
               int  nh,
               char *calling_func)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  FLOAT         
    *m;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="fvector";

#if defined(DEBUG)
  fprintf(stderr,"fvector() called by %s \n",calling_func);
  fflush(stderr);
#endif

  if (nh < nl) {
    sprintf(Message,"range (%d,%d)",nl,nh);
    util_error(dbmsname,Message);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (FLOAT *)calloc(len_safe,sizeof(FLOAT));

  if (!m) {
    sprintf(Message,"called by %s, nl=%d,nh=%d,len_safe=%d",calling_func,nl,nh,len_safe);
    util_error(dbmsname,Message);
  }
  m -= nl_safe;

  return m;
}

/*======================= end of fvector() ====================================*/

/*======================= free_fvector() ======================================*/

#undef  DEBUG

/*
 *  Frees memory allocated by fvector().
 */

void free_fvector(FLOAT *m, 
                  int     nl, 
                  int     nh,
                  char   *calling_func)
{
  int  
    nl_safe;

  if (m) {
    nl_safe = (nl < 0) ? nl : 0;
    m += nl_safe;
    free(m);
  }

#if defined(DEBUG)
  if (m) {
    fprintf(stderr,"free_fvector() called by %s \n",calling_func);
  }
  else {
    fprintf(stderr,"free_fvector() called by %s for NULL pointer\n",calling_func);
  }
#endif

  return;
}

/*======================= end of free_fvector() ===============================*/

/*======================= dvector() ============================================*/
      
/*
 *  Allocates memory for a 1D double array 
 *  with range [nl..nh].
 */

#undef  DEBUG 

double *dvector(int  nl, 
                int  nh,
                char *calling_func)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  double         
    *m;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="dvector";

#if defined(DEBUG)
  fprintf(stderr,"dvector() called by %s \n",calling_func);
  fflush(stderr);
#endif

  if (nh < nl) {
    sprintf(Message,"called by %s, range (%d,%d)",calling_func,nl,nh);
    util_error(dbmsname,Message);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (double *)calloc(len_safe,sizeof(double));

  if (!m) {
    sprintf(Message,"called by %s, nl=%d,nh=%d,len_safe=%d",calling_func,nl,nh,len_safe);
    util_error(dbmsname,Message);
  }
  m -= nl_safe;

  return m;
}

/*======================= end of dvector() ====================================*/

/*======================= free_dvector() ======================================*/

#undef  DEBUG

/*
 *  Frees memory allocated by dvector().
 */

void free_dvector(double *m, 
                  int     nl, 
                  int     nh,
                  char   *calling_func)
{
  int  
    nl_safe;

  if (m) {
    nl_safe = (nl < 0) ? nl : 0;
    m += nl_safe;
    free(m);
  }

#if defined(DEBUG)
  if (m) {
    fprintf(stderr,"free_dvector() called by %s \n",calling_func);
  }
  else {
    fprintf(stderr,"free_dvector() called by %s for NULL pointer\n",calling_func);
  }
#endif

  return;
}

/*======================= end of free_dvector() ===============================*/

/*======================= ftriplet() ==========================================*/
      
/*
 *  Allocates memory for a 1D float_triplet array 
 *  with range [nl..nh].
 */

/*
 * Define DEBUG here to enable this state.
 */
#undef  DEBUG 

float_triplet *ftriplet(int  nl, 
                        int  nh,
                        char *calling_func)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  float_triplet         
    *m;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="ftriplet";

#if defined(DEBUG)
  fprintf(stderr,"ftriplet() called by %s \n",calling_func);
  fflush(stderr);
#endif

  if (nh < nl) {
    sprintf(Message,"called by %s, range (%d,%d)",calling_func,nl,nh);
    util_error(dbmsname,Message);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (float_triplet *)calloc(len_safe,sizeof(float_triplet));

  if (!m) {
    sprintf(Message,"called by %s",calling_func);
    util_error(dbmsname,Message);
  }
  m -= nl_safe;

  return m;
}

/*======================= end of ftriplet() ===================================*/

/*======================= free_ftriplet() =====================================*/

#undef  DEBUG

/*
 *  Frees memory allocated by ftriplet().
 */

void free_ftriplet(float_triplet *m, 
                  int             nl, 
                  int             nh,
                  char           *calling_func)
{
  int  
    nl_safe;

  if (m) {
    nl_safe = (nl < 0) ? nl : 0;
    m += nl_safe;
    free(m);
  }

#if defined(DEBUG)
  if (m) {
    fprintf(stderr,"free_ftriplet() called by %s \n",calling_func);
  }
  else {
    fprintf(stderr,"free_ftriplet() called by %s for NULL pointer\n",calling_func);
  }
#endif

  return;
}

/*======================= end of free_ftriplet() ==============================*/

/*======================= spline() ============================================*/

/* 
 * Cubic spline routine. 
 * Adapted from Numerical Recipes in C, p. 115-116. 
 * The tridiagonal system is solved with pivoting
 * to avoid numerical instability.
 * Assumes zero-offset arrays.
 *
 * NOTE: We stripe the data into one array with float_triplet to get
 *       a cache-aware memory layout.
 */

void spline(int            n,
            float_triplet *table,
            FLOAT          y1_bot, 
            FLOAT          y1_top)
{
  int     
    j,jm1,jp1;
  FLOAT 
     dx_a,dx_b,dx_c,  
    *a,
    *b,
    *c,
    *r,
    *y2;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="spline";

  /*
   * Check validity of n.
   */
  if (n < 3) {
    sprintf(Message,"n=%d < 3",n);
    util_error(dbmsname,Message);
  }
  
  /*
   * Allocate memory.
   */
  a  = fvector(0,n-1,dbmsname);
  b  = fvector(0,n-1,dbmsname);
  c  = fvector(0,n-1,dbmsname);
  r  = fvector(0,n-1,dbmsname);
  y2 = fvector(0,n-1,dbmsname);

  for (j = 0; j < n; j++) {
    y2[j] = (table+j)->z;
  }

  for (j = 0; j < n; j++) {
    jm1 = j-1;
    jp1 = j+1;
    if (jm1 >= 0) {
      dx_a = (table+j)->x - (table+jm1)->x; 
      if (dx_a <= 0.) {
        sprintf(Message,"x[%d]=%g x[%d]=%g",jm1,(table+jm1)->x,j,(table+j)->x);
        util_error(dbmsname,Message);
      }
      if (jp1 < n) {
        dx_b = (table+jp1)->x - (table+jm1)->x;
        if (dx_b <= 0.) {
          sprintf(Message,"x[%d]=%g x[%d]=%g",jm1,(table+jm1)->x,jp1,(table+jp1)->x);
          util_error(dbmsname,Message);
        }
      }
    }
    if (jp1 < n) {
      dx_c = (table+jp1)->x - (table+j)->x;
      if (dx_c <= 0.) {
        sprintf(Message,"x[%d]=%g x[%d]= %g",j,(table+j)->x,jp1,(table+jp1)->x);
        util_error(dbmsname,Message);
      }
    }

    if (j == 0) {
      if (y1_bot > 0.99e+30) {
        y2[j] = 0.;
      }
      else {
        b[j] = dx_c/3.;
        c[j] = dx_c/6.;
        r[j] = ((table+jp1)->y - (table+j)->y)/dx_c-y1_bot;
      }
    }
    else if (j == n-1) {
      if (y1_top > 0.99e+30) {
        y2[j] = 0.;
      }
      else {
        a[j] = dx_a/6.;
        b[j] = dx_a/3.;
        r[j] = y1_top-((table+j)->y - (table+jm1)->y)/dx_a;
      }
    }
    else {
      a[j] = dx_a/6.;
      b[j] = dx_b/3.;
      c[j] = dx_c/6.;
      r[j] = ((table+jp1)->y - (table+j)->y)/dx_c - ((table+j)->y - (table+jm1)->y)/dx_a;
    }
  }

  if (y1_bot > 0.99e+30) {
    if (y1_top > 0.99e+30) {
      /* y2 = 0 on both ends. */
      tridiag(n-2,a+1,b+1,c+1,r+1,y2+1,WITH_PIVOTING);
    }
    else {
      /* y2 = 0 at start. */
      tridiag(n-1,a+1,b+1,c+1,r+1,y2+1,WITH_PIVOTING);
    }
  }
  else {
    if (y1_top > 0.99e+30) {
      /* y2 = 0 at end. */
      tridiag(n-1,a,b,c,r,y2,WITH_PIVOTING);
    }
    else {
      /* y2 needed at both ends */
      tridiag(n,a,b,c,r,y2,WITH_PIVOTING);
    }
  }

  for (j = 0; j < n; j++) {
    (table+j)->z = y2[j];
  }

  /* 
   * Free allocated memory.
   */
  free_fvector(y2,0,n-1,dbmsname);
  free_fvector(r, 0,n-1,dbmsname);
  free_fvector(c, 0,n-1,dbmsname);
  free_fvector(b, 0,n-1,dbmsname);
  free_fvector(a, 0,n-1,dbmsname);
  
  return;
}

/*======================= end of spline() ===================================*/

/*======================= splint() ==========================================*/

/*  
 *  Evaluates cubic-spline interpolations.
 *  This version assumes you have already found the correct position
 *  in the tables, unlike the Numerical Recipes version. 
 *  The function find_place_in_table() may be used to find the position.
 *
 *  NOTE: We stripe the data into one array with float_triplet to get
 *        a cache-aware memory layout.
 */

FLOAT splint(register FLOAT          xx, 
             register float_triplet *table,
             register FLOAT          dx)
{
  register FLOAT  
    a, 
    b,
    ans;

  a = ( (table+1)->x - xx       )/dx;
  b = (     xx       - table->x )/dx;

  ans = a*table->y + b*(table+1)->y +( (a*a*a-a)*(table  )->z
                                      +(b*b*b-b)*(table+1)->z )*dx*dx/6;
    
  return ans;
}

/*======================= end of splint() ===================================*/

/*======================= linint() ==========================================*/

/*
 * Evaluates linear interpolation using same arguments as splint(),
 * to make it easy to switch between the two.
 */

FLOAT linint(register FLOAT          xx,
             register float_triplet *table,
             register FLOAT          dx)
{
  FLOAT
    ans;

  ans = table->y + ((table+1)->y - table->y)*(xx - table->x)/dx;
  
  return ans;
}

/*======================= end of linint() ===================================*/

/*=========================== spline_pchip() ================================*/

/*
 *  Purpose:
 *
 *    Sets derivatives for a piecewise cubic Hermite interpolant.
 *
 *  Discussion:
 *
 *    This routine computes what would normally be called a Hermite 
 *    interpolant.  However, the user is only required to supply function
 *    values, not derivative values as well.  This routine computes
 *    "suitable" derivative values, so that the resulting Hermite interpolant
 *    has desirable shape and monotonicity properties.
 *
 *    The interpolant will have an extremum at each point where
 *    monotonicity switches direction.
 *
 *    The resulting piecewise cubic Hermite function may be evaluated
 *    by splint_pchip().
 *
 *    This routine was originally called "PCHIM".
 *
 *    The acronym PCHIP means Piecewise Cubic Hermite Interpolating Polynomial.
 *
 *  Modified:
 *
 *    25 August 2007
 *
 *  Author:
 *
 *    Fred Fritsch,
 *    Mathematics and Statistics Division,
 *    Lawrence Livermore National Laboratory.
 *
 *    C translation by John Burkardt.
 *
 *    Adapted for EPIC by Tim Dowling.
 *
 *  References:
 *
 *    Fred Fritsch, Ralph Carlson,
 *    Monotone Piecewise Cubic Interpolation,
 *    SIAM Journal on Numerical Analysis,
 *    Volume 17, Number 2, April 1980, pages 238-246.
 *
 *    Fred Fritsch, Judy Butland,
 *    A Method for Constructing Local Monotone Piecewise 
 *    Cubic Interpolants,
 *    SIAM Journal on Scientific and Statistical Computing,
 *    Volume 5, Number 2, 1984, pages 300-304.
 *
 *  Parameters:
 *
 *    Input, int n, the number of data points; n must be at least 2.
 *
 *    Input, 
 *      float_triplet table[n]:
 *      table[n].x are the strictly increasing independent variable values.
 *
 *      table[n].y are the dependent variable values to be interpolated. 
 *      This routine is designed for monotonic data, but it will work for
 *      any data. It will force extrema at points where monotonicity switches
 *      direction.
 *
 *    Output, 
 *      table[n].z are the derivative values at the data points.  If the
 *      data are monotonic, these values will determine a monotone cubic
 *      Hermite function.
 */

void spline_pchip(int            n,
                  float_triplet *table)
{
  int
    i,nless1;
  register FLOAT
    del1,del2,
    dmax,dmin,
    drat1,drat2,
    h1,h2,hsum,hsumt3,
    w1,w2;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="spline_pchip";

  /*
   * Check the arguments.
   */
  if ( n < 2 ) {
    sprintf(Message,"number of data points, n=%d, is less than 2",n);
    util_error(dbmsname,Message);
  }

  for ( i = 1; i < n; i++ ) {
    if ( table[i].x <= table[i-1].x ) {
      sprintf(Message,"x array is not strictly increasing, table[%d].x=%g <= table[%d].x=%g\n",i,table[i].x,i-1,table[i-1].x);
      util_error(dbmsname,Message);
    }
  }

  nless1 = n - 1;
  h1     =   table[1].x - table[0].x;
  del1   = ( table[1].y - table[0].y ) / h1;

  /*
   * Special case n = 2, use linear interpolation.
   */
  if ( n == 2 ) {
    table[0  ].z = del1;
    table[n-1].z = del1;
    return;
  }

  /*
   *  Normal case, 3 <= n.
   */
  h2   =   table[2].x - table[1].x;
  del2 = ( table[2].y - table[1].y ) / h2;

  /*
   * Set table[0].z via non-centered three point formula, adjusted to be shape preserving.
   */
  hsum       = h1 + h2;
  w1         = ( h1 + hsum ) / hsum;
  w2         = -h1 / hsum;
  table[0].z = w1 * del1 + w2 * del2;

  if ( pchst( table[0].z, del1 ) <= 0.0 ) {
    table[0].z = 0.0;
  }
  else if ( pchst( del1, del2 ) < 0.0 ) {
    /*
     *  Need to do this check only if monotonicity switches.
     */
    dmax = 3.0 * del1;

    if ( fabs( dmax ) < fabs( table[0].z ) ) {
      table[0].z = dmax;
    }
  }

  /*
   *  Loop through interior points.
   */
  h2   = h1;
  del2 = del1;
  for ( i = 2; i <= nless1; i++ ) {
    h1   = h2;
    h2   = table[i].x - table[i-1].x;
    hsum = h1 + h2;
    del1 = del2;
    del2 = ( table[i].y - table[i-1].y ) / h2;

    /*
     *  Set table[i-1].z = 0 unless data are strictly monotonic.
     */
    table[i-1].z = 0.0;
    if ( pchst( del1, del2 ) > 0.0 ) {
     /*
      *  Use Brodlie modification of Butland formula.
      */
      hsumt3       = 3.0 * hsum;
      w1           = ( hsum + h1 ) / hsumt3;
      w2           = ( hsum + h2 ) / hsumt3;
      dmax         = MAX ( fabs ( del1 ), fabs ( del2 ) );
      dmin         = MIN ( fabs ( del1 ), fabs ( del2 ) );
      drat1        = del1 / dmax;
      drat2        = del2 / dmax;
      table[i-1].z = dmin / ( w1 * drat1 + w2 * drat2 );
    }
  }

  /*
   *  Set table[n-1].z via non-centered three point formula, adjusted to be
   *  shape preserving.
   */
  w1           = -h2 / hsum;
  w2           = ( h2 + hsum ) / hsum;
  table[n-1].z = w1 * del1 + w2 * del2;

  if ( pchst( table[n-1].z, del2 ) <= 0.0 ) {
    table[n-1].z = 0.0;
  }
  else if ( pchst( del1, del2 ) < 0.0 ) {
    /*
     *  Need to do this check only if monotonicity switches.
     */
    dmax = 3.0 * del2;

    if ( fabs( dmax ) < fabs( table[n-1].z ) ) {
      table[n-1].z = dmax;
    }
  }

  return;
}

/*========================= end of spline_pchip() ==============================*/

/*=========================== periodic_spline_pchip() ==========================*/

/*
 *  Purpose:
 *
 *    Sets derivatives for a piecewise cubic Hermite interpolant, assuming periodic input.
 *    As for spline_pchip(), assume a zero-offset input vector of length n, x[0] to x[n-1],
 *    with periodic boundary conditions x[0] = x[n-1].
 *
 *  NOTE: The endpoints span the full periodic interval in order for the full range to
 *        behave as an interpolation.
 *
 *  For Discussion, Authorship, References and Parameters, see the leading comments
 *  in spline_pchip() above.
 */

void periodic_spline_pchip(int            n,
                           float_triplet *table)
{
  int
    i,nless1;
  register FLOAT
    del1,del2,
    dmax,dmin,
    drat1,drat2,
    h1,h2,hsum,hsumt3,
    w1,w2;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="periodic_spline_pchip";

  /*
   * Check the arguments.
   */
  if (n < 2) {
    sprintf(Message,"number of data points, n=%d, is less than 2",n);
    util_error(dbmsname,Message);
  }

  /*
   * Check for correct periodicity of table[].y
   */
  if (table[0].y != table[n-1].y) {
    sprintf(Message,"table[0].y=%g != table[%d].y=%g; these need to be periodic",table[0].y,n-1,table[n-1].y);
    util_error(dbmsname,Message);
  }

  /*
   * Check that table[].x is strictly increasing.
   */
  for (i = 1; i < n; i++) {
    if (table[i].x <= table[i-1].x) {
      sprintf(Message,"x array is not strictly increasing, table[%d].x=%g <= table[%d].x=%g\n",i,table[i].x,i-1,table[i-1].x);
      util_error(dbmsname,Message);
    }
  }

  nless1 = n-1;

  /*
   * Special case n = 2, use linear interpolation.
   */
  if (n == 2) {
    h1           =  table[1].x-table[0].x;
    del1         = (table[1].y-table[0].y)/h1;
    table[0  ].z = del1;
    table[n-1].z = del1;
    return;
  }

  /*
   * Set table[0].z assuming periodicity.
   */
  h1   = table[n-1].x-table[n-2].x;
  h2   = table[1  ].x-table[0  ].x;
  hsum = h1+h2;
  del1 = (table[0].y-table[n-2].y)/h1;
  del2 = (table[1].y-table[0  ].y)/h2;

  /*
   *  Set table[0].z = 0. unless data are strictly monotonic.
   */
  table[0].z = 0.;
  if (pchst(del1,del2) > 0.) {
   /*
    *  Use Brodlie modification of Butland formula.
    */
    hsumt3     = 3.*hsum;
    w1         = (hsum+h1)/hsumt3;
    w2         = (hsum+h2)/hsumt3;
    dmax       = MAX(fabs(del1),fabs(del2));
    dmin       = MIN(fabs(del1),fabs(del2));
    drat1      = del1/dmax;
    drat2      = del2/dmax;
    table[0].z = dmin/(w1*drat1+w2*drat2);
  }

  /*
   *  Loop through interior points.
   */
  for (i = 2; i <= nless1; i++) {
    h1   = h2;
    h2   = table[i].x-table[i-1].x;
    hsum = h1+h2;
    del1 = del2;
    del2 = (table[i].y-table[i-1].y)/h2;

    /*
     *  Set table[i-1].z = 0 unless data are strictly monotonic.
     */
    table[i-1].z = 0.;
    if (pchst(del1,del2) > 0.) {
     /*
      *  Use Brodlie modification of Butland formula.
      */
      hsumt3       = 3.*hsum;
      w1           = (hsum+h1)/hsumt3;
      w2           = (hsum+h2)/hsumt3;
      dmax         = MAX(fabs(del1),fabs(del2));
      dmin         = MIN(fabs(del1),fabs(del2));
      drat1        = del1/dmax;
      drat2        = del2/dmax;
      table[i-1].z = dmin/(w1*drat1+w2*drat2);
    }
  }

  /*
   *  Set table[n-1].z assuming periodicity.
   */
  table[nless1].z = table[0].z;

  return;
}

/*========================= end of periodic_spline_pchip() =====================*/

/*========================= splint_pchip() =====================================*/

/*
 *  Purpose:
 *
 *    Evaluates a cubic polynomial given in Hermite form.
 *
 *  Discussion:
 *
 *    This routine evaluates a cubic polynomial given in Hermite form, and
 *    is designed to work with spline_pchip().
 *
 *    The cubic polynomial is determined by function values
 *    table[0].y, table[1].y and derivatives table[0].z, table[1].z
 *    on the interval [table[0].x,table[1].x].
 *
 *    This routine was originally called "CHFEV".
 *
 *  Modified:
 *
 *    25 August 2007
 *
 *  Author:
 *
 *    Fred Fritsch,
 *    Mathematics and Statistics Division,
 *    Lawrence Livermore National Laboratory.
 *
 *    C translation by John Burkardt.
 *
 *    Modified by Tim Dowling.
 *
 *  References:
 *
 *    Fred Fritsch, Ralph Carlson, 
 *    Monotone Piecewise Cubic Interpolation,
 *    SIAM Journal on Numerical Analysis,
 *    Volume 17, Number 2, April 1980, pages 238-246.
 *
 *    David Kahaner, Cleve Moler, Steven Nash,
 *    Numerical Methods and Software,
 *    Prentice Hall, 1989,
 *    ISBN: 0-13-627258-4,
 *    LC: TA345.K34.
 *
 *  Parameters:
 *
 *    Input, FLOAT xx, the position of the interpolation
 *
 *    Input, float_triplet *table:
 *       The use of the data type float_triplet facilitates
 *       cache-aware memory allocation. 
 *
 *       NOTE: This function does not search the table,
 *             but assumes xx is between the given inputs as
 *             returned by find_place_in_table() or hunt_place_in_table().
 *
 *       table[0].x and table[1].x are the endpoints of the interval of
 *       definition of the cubic; they must be distinct.
 *
 *       table[0].y and table[1].y are the corresponding values of the function.
 *
 *       table[0].z and table[1].z are the corresponding derivative values.
 *
 *    Input, h = dx, the x interval, as returned by find_place_in_table().
 *    This is included here to make it easy to switch with splint().
 *
 *    Returns the value of the cubic function at the point xx.
 */

FLOAT splint_pchip(FLOAT          xx,
                   float_triplet *table,
                   FLOAT          h)
{
  register FLOAT 
    c2,c3,
    del1,del2,delta,
    x,f;

  /*
   *  Compute cubic coefficients expanded about table[0].x.
   */
  delta =  (table[1].y-table[0].y)/h;
  del1  =  (table[0].z-delta     )/h;
  del2  =  (table[1].z-delta     )/h;
  c3    =  (del1+del2            )/h;
  c2    = -(del1+del1+del2);

  /*
   *  Evaluation.
   */
  x = xx-table[0].x;
  f = table[0].y+x*(table[0].z+x*(c2+x*c3));

  return f;
}

/*========================= end of splint_pchip() ==============================*/

/*========================= pchst() ============================================*/

/*
 *  Purpose:
 *
 *    Sign-testing routine.
 *
 *  Discussion:
 *
 *    This routine essentially computes the sign of arg1 * arg2.
 *
 *    The object is to do this without multiplying arg1 * arg2, to avoid
 *    possible over/underflow problems.
 *
 *  Modified:
 *
 *    25 August 2007
 *
 *  Author:
 *
 *    Fred Fritsch,
 *    Mathematics and Statistics Division,
 *    Lawrence Livermore National Laboratory.
 *
 *    C translation by John Burkardt.
 *
 *    Adapted by Tim Dowling.
 *
 *  Reference:
 *
 *    Fred Fritsch, Ralph Carlson, 
 *    Monotone Piecewise Cubic Interpolation,
 *    SIAM Journal on Numerical Analysis,
 *    Volume 17, Number 2, April 1980, pages 238-246.
 *
 *  Parameters:
 *
 *    Input, FLOAT arg1, arg2, two values to check.
 *
 *    Output, FLOAT pchst,
 *    -1.0, if arg1 and arg2 are of opposite sign.
 *     0.0, if either argument is zero.
 *    +1.0, if arg1 and arg2 are of the same sign.
 */

FLOAT pchst(FLOAT arg1,
            FLOAT arg2)
{
  FLOAT
    value;

  if (arg1 == 0.) {
    value = 0.;
  }
  else if (arg1 < 0.) {
    if (arg2 < 0.) {
      value = 1.;
    }
    else if (arg2 == 0.) {
      value = 0.;
    }
    else {
      value = -1.;
    }
  }
  else {
    if (arg2 < 0.) {
      value = -1.;
    }
    else if (arg2 == 0.) {
      value = 0.;
    }
    else {
      value = 1.;
    }
  }

  return value;
}

/*======================= end of pchst() ====================================*/

/*======================= lagrange_interp() =================================*/
/*
 * Returns the lagrangian interpolation of f(x).
 *
 * Input: "order" is the order of polynomial interpolation 
 *        (0=copy,1=linear,2=quadradic,3=cubic,etc.).
 *
 *        f[order+2] where f[0] is empty,
 *                         f[1] = f(x1),
 *                         f[2] = f(x2),
 *                         f[3] = f(x3), etc.
 *
 *        x[order+2] where x[0] is the point of interpolation, i.e., "x" in f(x),
 *                          x[1] = x1,
 *                          x[2] = x2,
 *                          x[3] = x3, etc.
 */

FLOAT lagrange_interp(register FLOAT *f,
                      register FLOAT *x,
                      register int order)
{
  register int
    jj, il;
  FLOAT 
    lp[ order+2 ];

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lagrange_interp";

  if ( order < 0 ) {
    sprintf(Message,"order=%d is not in valid (must be >= 0)",order);
    util_error(dbmsname,Message);
  }

  for (jj=1; jj<=order+1; jj++) {
    lp[jj] = 1.;
    for (il=1; il<=order+1; il++) {
      if (il != jj)
      lp[jj] *= (x[0]-x[il])/(x[jj]-x[il]);
    }
  }

  f[0] = 0.0;
  for (jj=1; jj<=order+1; jj++) {
    f[0] += f[jj]*lp[jj];
  }

  return f[0];
}

/*==================== end of lagrange_interp() =============================*/

/*=============================== gamma_nr() ================================*/

/*
 * CJP 06/2003 *A*  
 * Based on "Numerical Recipes in C", 2nd Ed., p 213-214, Cambridge.
 */

FLOAT gamma_nr(FLOAT xx_input)
{
  register int
    j;
  register double 
    x,xx,
    y,tmp,ser,loggamma;
  static double 
    cof[6]={76.18009172947146,     -86.505320032941677,
            24.01409824083091,      -1.231739572450155,
             0.1208650973866179e-2  -0.5395239384953e-5};
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="gamma_nr";
    
  if (xx_input <= 0.) {
    fprintf(stderr,"%s: Negative or zero argument: xx=%e\n",dbmsname,xx_input);
  }  
  /*
   * Use double precision for internal calculations.
   */
  xx   = (double)xx_input;
  y    = x = xx;
  tmp  = x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser  = 1.000000000190015;
  for  (j = 0; j <= 5; j++) {
    ser += cof[j]/++y;
  }
  loggamma = -tmp+log(2.5066282746310005*ser/x);

  return (FLOAT)exp(loggamma);
}

/*============================end of gamma_nr() =============================*/

/*======================= sech2() ===========================================*/

FLOAT sech2(FLOAT xx)
     /*
      *  Evaluates the square of the hyperbolic secant
      */
{
  FLOAT 
    a;

  a = 1./cosh(xx);
  return a*a;
}

/*======================= end of sech2() =====================================*/

/*===================== exp_integral_setup() ================================*/

/*
 * Populate the exponential integral tables with x values ranging
 * from 10^-8 to 10^9. Similarly, populate y values in exp3table
 * with corresponding third-order exponential integral values, and
 * populate exp4table with corresponding fourth-order exponential
 * integral values.
 */

void exp_integral_setup(float_triplet *exp3table,
                        float_triplet *exp4table,
                        int num_datapoints)
{
  int 
    q;
  EPIC_FLOAT
    xtable[501]={0.0,1.0000000133514320E-10, 1.0966806129791248E-10, 1.2027083508263903E-10, 1.3189869138090339E-10, 1.4465073578345911E-10,
 1.5863565546887220E-10, 1.7397264555716406E-10, 1.9079242502387863E-10, 2.0923835083332934E-10, 2.2946763978691522E-10,
 2.5165270850046318E-10, 2.7598264293138146E-10, 3.0266480998057774E-10, 3.3192662490501456E-10, 3.6401748960476869E-10,
 3.9921091830484256E-10, 4.3780686874916542E-10, 4.8013429877582553E-10, 5.2655397006352521E-10, 5.7746152294592025E-10,
 6.3329084850085110E-10, 6.9451778665517647E-10, 7.6166418182459659E-10, 8.3530233065515997E-10, 9.1605985977508146E-10,
 1.0046250751305239E-09, 1.1017528284983645E-09, 1.2082709511769238E-09, 1.3250873096897659E-09, 1.4531975436391641E-09,
 1.5936935516598661E-09, 1.7477727977999503E-09, 1.9167485176481519E-09, 2.1020609111956811E-09, 2.3052894178306228E-09,
 2.5281661780861389E-09, 2.7725907968785406E-09, 3.0406465340640516E-09, 3.3346180603082896E-09, 3.6570109296038233E-09,
 4.0105729344023886E-09, 4.3983175253740625E-09, 4.8235494954027595E-09, 5.2898931467257979E-09, 5.8013231812893488E-09,
 6.3621985776019489E-09, 6.9772997428224800E-09, 7.6518692567342347E-09, 8.3916565548707874E-09, 9.2029669316332006E-09,
 1.0092715281057908E-08, 1.1068485033274243E-08, 1.2138592788973871E-08, 1.3312159202779461E-08, 1.4599186719659832E-08,
 1.6010644826948210E-08, 1.7558563548576651E-08, 1.9256135978388971E-08, 2.1117830726435721E-08, 2.3159515236644231E-08,
 2.5398591026914819E-08, 2.7854142004310340E-08, 3.0547097119447167E-08, 3.3500408746409748E-08, 3.6739248308542308E-08,
 4.0291220817458961E-08, 4.4186600153813794E-08, 4.8458587095156477E-08, 5.3143592290075464E-08, 5.8281546590450547E-08,
 6.3916241386813387E-08, 7.0095701847534227E-08, 7.6872596242996580E-08, 8.4304684843477915E-08, 9.2455312216746503E-08,
 1.0139394712128353E-07, 1.1119677459670605E-07, 1.2194734529784548E-07, 1.3373728760683390E-07, 1.4666708859261249E-07,
 1.6084695047407276E-07, 1.7639772988656679E-07, 1.9345196795726538E-07, 2.1215501996881855E-07, 2.3266629423957396E-07,
 2.5516061077947784E-07, 2.7982969130163863E-07, 3.0688379328910225E-07, 3.3655350204418053E-07, 3.6909169599418045E-07,
 4.0477570200406743E-07, 4.4390965906604659E-07, 4.8682711051206733E-07, 5.3389384684298286E-07, 5.8551102340413489E-07,
 6.4211857947967827E-07, 7.0419898794699363E-07, 7.7228136744992618E-07, 8.4694599213940009E-07, 9.2882923741845013E-07,
 1.0186290038448420E-06, 1.1171106654198511E-06, 1.2251135929611383E-06, 1.3435583081592535E-06, 1.4734543292925594E-06,
 1.6159087754706149E-06, 1.7721358027409086E-06, 1.9434669524840887E-06, 2.1313625002981784E-06, 2.3374239020997728E-06,
 2.5634074435216748E-06, 2.8112392089425298E-06, 3.0830314977311616E-06, 3.3811008276232187E-06, 3.7079876786751165E-06,
 4.0664781460751599E-06, 4.4596276863615106E-06, 4.8907871594386008E-06, 5.3636313883514040E-06, 5.8821904802356204E-06,
 6.4508841753962389E-06, 7.0745595172753768E-06, 7.7585321643753403E-06, 8.5086316962429967E-06, 9.3312512996638049E-06,
 1.0233402258546995E-05, 1.1222773711925170E-05, 1.2307798189393386E-05, 1.3497723482553353E-05, 1.4802691465032227E-05,
 1.6233824532867583E-05, 1.7803320402000974E-05, 1.9524556070849389E-05, 2.1412201834040869E-05, 2.3482346319065706E-05,
 2.5752633611546517E-05, 2.8242413637861015E-05, 3.0972907086846724E-05, 3.3967386276235338E-05, 3.7251373505361269E-05,
 4.0852858584729810E-05, 4.4802537396472919E-05, 4.9134073518973902E-05, 5.3884385145515931E-05, 5.9093959742398745E-05,
 6.4807199128388601E-05, 7.1072797916657727E-05, 7.7944158544715630E-05, 8.5479846429682393E-05, 9.3744089128247821E-05,
 1.0280732375571144E-04, 1.1274679732982590E-04, 1.2364722515623901E-04, 1.3560151286703782E-04, 1.4871154826640747E-04,
 1.6308906973240602E-04, 1.7885661857633812E-04, 1.9614858347680736E-04, 2.1511234589027454E-04, 2.3590953620059347E-04,
 2.5871741131384024E-04, 2.8373036543982775E-04, 3.1116158693688902E-04, 3.4124487534139663E-04, 3.7423663406879336E-04,
 4.1041805577013955E-04, 4.5009751897025143E-04, 4.9361321641423490E-04, 5.4133603752414873E-04, 5.9367272953329364E-04,
 6.5106936424085362E-04, 7.1401513993446476E-04, 7.8304655088491095E-04, 8.5875195995008626E-04, 9.4177661326101810E-04,
 1.0328281397307640E-03, 1.1326825822590880E-03, 1.2421910120376470E-03, 1.3622867823301259E-03, 1.4939934835522084E-03,
 1.6384336674534199E-03, 1.7968384147580847E-03, 1.9705578278115688E-03, 2.1610725376617892E-03, 2.3700063236523081E-03,
 2.5991399530851727E-03, 2.8504263589107107E-03, 3.1260072848053769E-03, 3.4282315395057742E-03, 3.7596750159825123E-03,
 4.1231626460798733E-03, 4.5217924777431065E-03, 4.9589620800464555E-03, 5.4383975010751946E-03, 5.9641850254728689E-03,
 6.5408060023273845E-03, 7.1731750402377710E-03, 7.8666818951030563E-03, 8.6272374076474222E-03, 9.4613238822129995E-03,
 1.0376050336204224E-02, 1.1379213091081926E-02, 1.2479362221330467E-02, 1.3685874427751224E-02, 1.5009032956189414E-02,
 1.6460115242851604E-02, 1.8051489033224326E-02, 1.9796717793827894E-02, 2.1710676315241241E-02, 2.3809677491696869E-02,
 2.6111611357802020E-02, 2.8636097567414004E-02, 3.1404652614264568E-02, 3.4440873219573204E-02, 3.7770637449685812E-02,
 4.1422325277886245E-02, 4.5427060470256728E-02, 4.9818975857206307E-02, 5.4635504251611081E-02, 5.9917697493093176E-02,
 6.5710576337692991E-02, 7.2063514175081694E-02, 7.9030657843785754E-02, 8.6671389131079385E-02, 9.5050830890959809E-02,
 1.0424040209391441E-01, 1.1431842653922947E-01, 1.2537080041797166E-01, 1.3749172441635343E-01, 1.5078450659929446E-01,
 1.6536244291724325E-01, 1.8134978284092879E-01, 1.9888278835425116E-01, 2.1811089533124428E-01, 2.3919798719565935E-01,
 2.6232379171871056E-01, 2.8768541286007188E-01, 3.1549901070817171E-01, 3.4600164383812820E-01, 3.7945328978993892E-01,
 4.1613906088772379E-01, 4.5637163428569016E-01, 5.0049391695240575E-01, 5.4888196830732794E-01, 6.0194820541959659E-01,
 6.6014491808736819E-01, 7.2396812375709874E-01, 7.9396179513875598E-01, 8.7072249654938250E-01, 9.5490446850114963E-01,
 1.0472252038703913E+00, 1.1484715631740312E+00, 1.2595064810745713E+00, 1.3812763212740260E+00, 1.5148189424833491E+00,
 1.6612725442146634E+00, 1.8218853677899436E+00, 1.9980263352490628E+00, 2.1911967168337654E+00, 2.4030429264910280E+00,
 2.6353705544533632E+00, 2.8901597564972823E+00, 3.1695821310445811E+00, 3.4760192279518076E+00, 3.8120828467406520E+00,
 4.1806372972732557E+00, 4.5848238126030676E+00, 5.0280873220746400E+00, 5.5142058128628451E+00, 6.0473225302030693E+00,
 6.6319812907591071E+00, 7.2731652101086235E+00, 7.9763391744248295E+00, 8.7474964183460351E+00, 9.5932096060213148E+00,
 1.0520686850702461E+01, 1.1537833150343213E+01, 1.2653317762829364E+01, 1.3876648095086191E+01, 1.5218250735829306E+01,
 1.6689560322609790E+01, 1.8303117000578911E+01, 2.0072673303624612E+01, 2.2013311368839549E+01, 2.4141571483354774E+01,
 2.6475593059157713E+01, 2.9035269237440552E+01, 3.1842416440188416E+01, 3.4920960314117600E+01, 3.8297139651784413E+01,
 4.1999730027909592E+01, 4.6060289056994883E+01, 5.0513425362594113E+01, 5.5397093550694514E+01, 6.0752917701299218E+01,
 6.6626546135371413E+01, 7.3068040480857874E+01, 8.0132302353849823E+01, 8.7879541291530984E+01, 9.6375787925165838E+01,
 1.0569345676695771E+02, 1.1591196340747925E+02, 1.2711840138411567E+02, 1.3940828448955776E+02, 1.5288636084712476E+02,
 1.6766750569137136E+02, 1.8387770046326364E+02, 2.0165510656485884E+02, 2.2115124292523461E+02, 2.4253227741418871E+02,
 2.6598044311059613E+02, 2.9169559149643311E+02, 3.1989689581453536E+02, 3.5082471910799188E+02, 3.8474266286271524E+02,
 4.2193981371400304E+02, 4.6273320736602017E+02, 5.0747053072450581E+02, 5.5653308527327761E+02, 6.1033903695176207E+02,
 6.6934698023260739E+02, 7.3405984677647257E+02, 8.0502919197788958E+02, 8.8285989593695876E+02, 9.6821531892377038E+02,
 1.0618229552762391E+03, 1.1644806339201693E+03, 1.2770633183592686E+03, 1.4005305640923107E+03, 1.5359346970178228E+03,
 1.6844297825318645E+03, 1.8472814617634883E+03, 2.0258777387833895E+03, 2.2217408107263032E+03, 2.4365400416567641E+03,
 2.6721061907558815E+03, 2.9304470160978017E+03, 3.2137643870086317E+03, 3.5244730508583489E+03, 3.8652212142375070E+03,
 4.2389131139343335E+03, 4.6487337700875087E+03, 5.0981761324883810E+03, 5.5910708514041062E+03, 6.1316189266616602E+03,
 6.7244275132649846E+03, 7.3745491887203516E+03, 8.0875250167504792E+03, 8.8694317744343189E+03, 9.7269337452946711E+03,
 1.0667339519772848E+04, 1.1698664287206975E+04, 1.2829698150236314E+04, 1.4070081043882592E+04, 1.5430384897853410E+04,
 1.6922203742346006E+04, 1.8558252525334428E+04, 2.0352475483570488E+04, 2.2320164990958459E+04, 2.4478091897261977E+04,
 2.6844648468034145E+04, 2.9440005144066621E+04, 3.2286282456435365E+04, 3.5407739562396717E+04, 3.8830981009045281E+04,
 4.2585183487007402E+04, 4.6702344506816473E+04, 5.1217555117470103E+04, 5.6169298991578718E+04, 6.1599780426246827E+04,
 6.7555284055276163E+04, 7.4086569338557456E+04, 8.1249303190935170E+04, 8.9104534437888156E+04, 9.7719214141871387E+04,
 1.0716676623411621E+05, 1.1752771331545799E+05, 1.2889036295996004E+05, 1.4135156037078530E+05, 1.5501751380328296E+05,
 1.7000469979047761E+05, 1.8644085588630088E+05, 2.0446606938782654E+05, 2.2423397131582952E+05, 2.4591304583009542E+05,
 2.6968806623979105E+05, 2.9576166984816740E+05, 3.2435608505419997E+05, 3.5571502543144312E+05, 3.9010576692756987E+05,
 4.2782142588879075E+05, 4.6918345732504985E+05, 5.1454439470899681E+05, 5.6429085466034501E+05, 6.1884683212493081E+05,
 6.7867731413364271E+05, 7.4429224294174288E+05, 8.1625086232684145E+05, 8.9516648408958630E+05, 9.8171171538247541E+05,
 1.0766241914199635E+06, 1.1807128624303939E+06, 1.2948648884342024E+06, 1.4200532006134253E+06, 1.5573447937189180E+06,
 1.7079098201924718E+06, 1.8730315635140955E+06, 2.0541173757784809E+06, 2.2527106727228863E+06, 2.4705040884415871E+06,
 2.7093539019058025E+06, 2.9712958582483879E+06, 3.2585625196597050E+06, 3.5736022937785028E+06, 3.9191003017590186E+06,
 4.2980012638752516E+06, 4.7135345977193331E+06, 5.1692419429083746E+06, 5.6690073468967983E+06, 6.2170903691708911E+06,
 6.8181623859767225E+06, 7.4773464050108921E+06, 8.2002607294193190E+06, 8.9930668432579134E+06, 9.8625219265475534E+06,
 1.0816036447516667E+07, 1.1861737322895968E+07, 1.3008537184588330E+07, 1.4266210343081746E+07, 1.5645476095050527E+07,
 1.7158090085185066E+07, 1.8816944500938669E+07, 2.0636177954161089E+07, 2.2631295986155499E+07, 2.4819303223235812E+07,
 2.7218848309162449E+07, 2.9850382849732768E+07, 3.2736335724229291E+07, 3.5901304249404423E+07, 3.9372263825310849E+07,
 4.3178797849818990E+07, 4.7353349861406229E+07, 5.1931500059262104E+07, 5.6952268557522751E+07, 6.2458447958305635E+07,
 6.8496968078106850E+07, 7.5119295936159685E+07, 8.2381874413911983E+07, 9.0346603324361861E+07, 9.9081366991464972E+07,
 1.0866061283623478E+08, 1.1916598590089288E+08, 1.3068702471919313E+08, 1.4332192446390998E+08, 1.5717837387587246E+08,
 1.7237447310780799E+08, 1.8903974030587351E+08, 2.0731621550809187E+08, 2.2735967126834959E+08, 2.4934094032424560E+08,
 2.7344737162467134E+08, 2.9988442712698895E+08, 3.2887743297353560E+08, 3.6067349997291082E+08, 3.9554362975454479E+08,
 4.3378502454755139E+08, 4.7572362027037841E+08, 5.2171686452110040E+08, 5.7215676314543593E+08, 6.2747322134882343E+08,
 6.8813770782919073E+08, 7.5466727316027880E+08, 8.2762895667466843E+08, 9.0764461940689456E+08, 9.9539624428838527E+08,
 1.0916317487684593E+09, 1.1971713594029174E+09, 1.3129146027417550E+09, 1.4398479721000414E+09, 1.5790533355568018E+09,
 1.7317171568442633E+09, 1.8991406077181740E+09, 2.0827506579982133E+09, 2.2841122377999897E+09, 2.5049415756190591E+09,
 2.7471208259488106E+09, 3.0127141111052170E+09, 3.3039851139848003E+09, 3.6234163717009315E+09, 3.9737304345405865E+09,
 4.3579130705813999E+09, 4.7792387137452745E+09, 5.2412983721848993E+09, 5.7480302348697958E+09, 6.3037532372354250E+09,
 6.9132038719793253E+09, 7.5815765587469902E+09, 8.3145679167834330E+09, 9.1184253178905945E+09, 1.0000000133514320E+10},
    e3table[501]={0.5,4.9999999989999999E-01, 4.9999999989033195E-01, 4.9999999987972915E-01, 4.9999999986810129E-01, 4.9999999985534926E-01,
 4.9999999984136434E-01, 4.9999999982602733E-01, 4.9999999980920756E-01, 4.9999999979076165E-01, 4.9999999977053233E-01,
 4.9999999974834730E-01, 4.9999999972401737E-01, 4.9999999969733516E-01, 4.9999999966807340E-01, 4.9999999963598252E-01,
 4.9999999960078906E-01, 4.9999999956219315E-01, 4.9999999951986568E-01, 4.9999999947344603E-01, 4.9999999942253848E-01,
 4.9999999936670914E-01, 4.9999999930548222E-01, 4.9999999923833582E-01, 4.9999999916469767E-01, 4.9999999908394016E-01,
 4.9999999899537495E-01, 4.9999999889824714E-01, 4.9999999879172907E-01, 4.9999999867491268E-01, 4.9999999854680244E-01,
 4.9999999840630643E-01, 4.9999999825222724E-01, 4.9999999808325152E-01, 4.9999999789793914E-01, 4.9999999769471065E-01,
 4.9999999747183388E-01, 4.9999999722740923E-01, 4.9999999695935360E-01, 4.9999999666538203E-01, 4.9999999634298919E-01,
 4.9999999598942724E-01, 4.9999999560168268E-01, 4.9999999517645072E-01, 4.9999999471010714E-01, 4.9999999419867713E-01,
 4.9999999363780179E-01, 4.9999999302270076E-01, 4.9999999234813131E-01, 4.9999999160834413E-01, 4.9999999079703389E-01,
 4.9999998990728572E-01, 4.9999998893151615E-01, 4.9999998786140859E-01, 4.9999998668784246E-01, 4.9999998540081531E-01,
 4.9999998398935763E-01, 4.9999998244143934E-01, 4.9999998074386748E-01, 4.9999997888217346E-01, 4.9999997684048969E-01,
 4.9999997460141493E-01, 4.9999997214586511E-01, 4.9999996945291136E-01, 4.9999996649960143E-01, 4.9999996326076385E-01,
 4.9999995970879374E-01, 4.9999995581341727E-01, 4.9999995154143378E-01, 4.9999994685643268E-01, 4.9999994171848328E-01,
 4.9999993608379434E-01, 4.9999992990434089E-01, 4.9999992312745489E-01, 4.9999991569537633E-01, 4.9999990754476092E-01,
 4.9999989860614041E-01, 4.9999988880333007E-01, 4.9999987805277996E-01, 4.9999986626286214E-01, 4.9999985333309060E-01,
 4.9999983915326379E-01, 4.9999982360252643E-01, 4.9999980654833859E-01, 4.9999978784534660E-01, 4.9999976733414414E-01,
 4.9999974483991344E-01, 4.9999972017093558E-01, 4.9999969311695636E-01, 4.9999966344739427E-01, 4.9999963090937577E-01,
 4.9999959522557946E-01, 4.9999955609187308E-01, 4.9999951317472130E-01, 4.9999946610834312E-01, 4.9999941449159468E-01,
 4.9999935788455024E-01, 4.9999929580475333E-01, 4.9999922772310468E-01, 4.9999915305935350E-01, 4.9999907117715198E-01,
 4.9999898137863280E-01, 4.9999888289846173E-01, 4.9999877489731503E-01, 4.9999865645472769E-01, 4.9999852656124888E-01,
 4.9999838410984004E-01, 4.9999822788644138E-01, 4.9999805655962637E-01, 4.9999786866925666E-01, 4.9999766261404022E-01,
 4.9999743663788676E-01, 4.9999718881494554E-01, 4.9999691703319571E-01, 4.9999661897645215E-01, 4.9999629210463192E-01,
 4.9999593363211364E-01, 4.9999554050400619E-01, 4.9999510937012454E-01, 4.9999463655645082E-01, 4.9999411803383897E-01,
 4.9999354938369500E-01, 4.9999292576034321E-01, 4.9999224184975682E-01, 4.9999149182430291E-01, 4.9999066929311647E-01,
 4.9998976724768301E-01, 4.9998877800216635E-01, 4.9998769312797592E-01, 4.9998650338201728E-01, 4.9998519862801788E-01,
 4.9998376775025966E-01, 4.9998219855898840E-01, 4.9998047768669729E-01, 4.9997859047440579E-01, 4.9997652084697242E-01,
 4.9997425117638433E-01, 4.9997176213187028E-01, 4.9996903251556718E-01, 4.9996603908235543E-01, 4.9996275634234022E-01,
 4.9995915634431576E-01, 4.9995520843838664E-01, 4.9995087901574975E-01, 4.9994613122344905E-01, 4.9994092465170820E-01,
 4.9993521499121885E-01, 4.9992895365751505E-01, 4.9992208737929228E-01, 4.9991455774723381E-01, 4.9990630071958497E-01,
 4.9989724608036135E-01, 4.9988731684569299E-01, 4.9987642861338683E-01, 4.9986448885033236E-01, 4.9985139611187490E-01,
 4.9983703918674072E-01, 4.9982129616050491E-01, 4.9980403338995233E-01, 4.9978510437998025E-01, 4.9976434855393737E-01,
 4.9974158990746659E-01, 4.9971663553502849E-01, 4.9968927401731683E-01, 4.9965927365673068E-01, 4.9962638054693936E-01,
 4.9959031646135682E-01, 4.9955077654403174E-01, 4.9950742678504217E-01, 4.9945990126097428E-01, 4.9940779911943206E-01,
 4.9935068128479515E-01, 4.9928806686058080E-01, 4.9921942920180457E-01, 4.9914419162863460E-01, 4.9906172275043675E-01,
 4.9897133136697752E-01, 4.9887226091112652E-01, 4.9876368339486093E-01, 4.9864469281775997E-01, 4.9851429799447927E-01,
 4.9837141475496111E-01, 4.9821485746837502E-01, 4.9804332983904864E-01, 4.9785541491998236E-01, 4.9764956428699941E-01,
 4.9742408631425006E-01, 4.9717713348974091E-01, 4.9690668870791571E-01, 4.9661055047519614E-01, 4.9628631696396647E-01,
 4.9593136885090883E-01, 4.9554285087712852E-01, 4.9511765207036140E-01, 4.9465238457405647E-01, 4.9414336103459933E-01,
 4.9358657050680843E-01, 4.9297765284953848E-01, 4.9231187159829665E-01, 4.9158408532082432E-01, 4.9078871748529457E-01,
 4.8991972489988383E-01, 4.8897056481786622E-01, 4.8793416084498453E-01, 4.8680286783672733E-01, 4.8556843603341610E-01,
 4.8422197475194123E-01, 4.8275391603575646E-01, 4.8115397876097321E-01, 4.7941113380710743E-01, 4.7751357102794256E-01,
 4.7544866890218651E-01, 4.7320296790637290E-01, 4.7076214883471612E-01, 4.6811101749297523E-01, 4.6523349741588005E-01,
 4.6211263249974782E-01, 4.5873060170193569E-01, 4.5506874823425147E-01, 4.5110762596352888E-01, 4.4682706602379946E-01,
 4.4220626693149129E-01, 4.3722391176680830E-01, 4.3185831622564913E-01, 4.2608761153810498E-01, 4.1988996636801590E-01,
 4.1324385182378703E-01, 4.0612835358956278E-01, 3.9852353488559467E-01, 3.9041085344073284E-01, 3.8177363485401367E-01,
 3.7259760357756505E-01, 3.6287147120625168E-01, 3.5258757974751093E-01, 3.4174259500414561E-01, 3.3033824208308599E-01,
 3.1838207130441309E-01, 3.0588823841659518E-01, 2.9287827805286404E-01, 2.7938184386807208E-01, 2.6543738292665964E-01,
 2.5109270590999549E-01, 2.3640540892363354E-01, 2.2144309758931485E-01, 2.0628336031730438E-01, 1.9101343592939429E-01,
 1.7572952203145092E-01, 1.6053567569440941E-01, 1.4554226810700291E-01, 1.3086397083742859E-01, 1.1661727388854144E-01,
 1.0291756512262613E-01, 8.9875836485256003E-02, 7.7595123434194960E-02, 6.6166827677873560E-02, 5.5667115836853728E-02,
 4.6153622919623584E-02, 3.7662713182032861E-02, 3.0207555176368986E-02, 2.3777246044091743E-02, 1.8337167517958171E-02,
 1.3830671092128384E-02, 1.0182075730820113E-02, 7.3008277891896574E-03, 5.0865354148132971E-03, 3.4344697042746376E-03,
 2.2410455992717877E-03, 1.4087775683716691E-03, 8.5026075467549238E-04, 4.9085599031933500E-04, 2.6993893983773612E-04,
 1.4077682307248664E-04, 6.9279727591009345E-05, 3.1999094669105979E-05, 1.3789367958869503E-05, 5.5080007414986386E-06,
 2.0247856906796911E-06, 6.7965714342261872E-07, 2.0653095948282847E-07, 5.6280733067289823E-08, 1.3611603529708075E-08,
 2.8886273651071479E-09, 5.3123143879146993E-10, 8.3509606735025826E-11, 1.1054064464640461E-11, 1.2119341833395253E-12,
 1.0808141032977697E-13, 7.6863288236566162E-15, 4.2650871116602505E-16, 1.8030425841131320E-17, 5.6568659765699998E-19,
 1.2798460366328279E-20, 2.0233103042674566E-22, 2.1591124288850023E-24, 1.4973629146008960E-26, 6.4737287251839906E-29,
 1.6670381546465146E-31, 2.4320400406683830E-34, 1.9028203615618248E-37, 7.5177614556341401E-41, 1.4040286684565986E-44,
 1.1529750863415696E-48, 3.8453779450546871E-53, 4.7743784475376273E-58, 2.0057622026503778E-63, 2.5676862973101026E-69,
 8.9293710268395008E-76, 7.4371884136663089E-83, 1.2921359169093143E-90, 4.0245755717631016E-99,1.9031941050362146E-108,
1.1388287842191206E-118,7.0608255450046438E-130,3.6432747192463110E-142,1.2302193534643889E-155,2.0885727119947121E-170,
1.3352019789713973E-186,2.3409396842212579E-204,7.9502920574365268E-224,3.5721562286206707E-245,1.3977419734875481E-268,
3.0109907506062429E-294,2.1738888417014848E-322, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0},
    e4table[501]={3.3333333333333333E-01,3.3333333328333331E-01, 3.3333333327849929E-01, 3.3333333327319792E-01, 3.3333333326738396E-01, 3.3333333326100795E-01,
 3.3333333325401548E-01, 3.3333333324634701E-01, 3.3333333323793712E-01, 3.3333333322871411E-01, 3.3333333321859948E-01,
 3.3333333320750697E-01, 3.3333333319534197E-01, 3.3333333318200092E-01, 3.3333333316737002E-01, 3.3333333315132457E-01,
 3.3333333313372787E-01, 3.3333333311442986E-01, 3.3333333309326618E-01, 3.3333333307005630E-01, 3.3333333304460255E-01,
 3.3333333301668788E-01, 3.3333333298607443E-01, 3.3333333295250123E-01, 3.3333333291568212E-01, 3.3333333287530337E-01,
 3.3333333283102079E-01, 3.3333333278245691E-01, 3.3333333272919785E-01, 3.3333333267078968E-01, 3.3333333260673453E-01,
 3.3333333253648656E-01, 3.3333333245944691E-01, 3.3333333237495905E-01, 3.3333333228230289E-01, 3.3333333218068861E-01,
 3.3333333206925020E-01, 3.3333333194703790E-01, 3.3333333181301006E-01, 3.3333333166602430E-01, 3.3333333150482786E-01,
 3.3333333132804682E-01, 3.3333333113417457E-01, 3.3333333092155859E-01, 3.3333333068838672E-01, 3.3333333043267172E-01,
 3.3333333015223404E-01, 3.3333332984468345E-01, 3.3333332950739875E-01, 3.3333332913750507E-01, 3.3333332873184990E-01,
 3.3333332828697571E-01, 3.3333332779909086E-01, 3.3333332726403697E-01, 3.3333332667725385E-01, 3.3333332603374005E-01,
 3.3333332532801102E-01, 3.3333332455405168E-01, 3.3333332370526547E-01, 3.3333332277441818E-01, 3.3333332175357600E-01,
 3.3333332063403814E-01, 3.3333331940626271E-01, 3.3333331805978522E-01, 3.3333331658312948E-01, 3.3333331496370983E-01,
 3.3333331318772375E-01, 3.3333331124003424E-01, 3.3333330910404091E-01, 3.3333330676153855E-01, 3.3333330419256174E-01,
 3.3333330137521466E-01, 3.3333329828548486E-01, 3.3333329489703811E-01, 3.3333329118099447E-01, 3.3333328710568150E-01,
 3.3333328263636491E-01, 3.3333327773495219E-01, 3.3333327235966809E-01, 3.3333326646469846E-01, 3.3333325999979979E-01,
 3.3333325290987104E-01, 3.3333324513448392E-01, 3.3333323660736802E-01, 3.3333322725584580E-01, 3.3333321700021329E-01,
 3.3333320575306047E-01, 3.3333319341852680E-01, 3.3333317989148376E-01, 3.3333316505663890E-01, 3.3333314878755343E-01,
 3.3333313094556427E-01, 3.3333311137860233E-01, 3.3333308991989657E-01, 3.3333306638655241E-01, 3.3333304057799301E-01,
 3.3333301227424972E-01, 3.3333298123408733E-01, 3.3333294719294781E-01, 3.3333290986069591E-01, 3.3333286891914599E-01,
 3.3333282401935022E-01, 3.3333277477862455E-01, 3.3333272077728732E-01, 3.3333266155508179E-01, 3.3333259660725417E-01,
 3.3333252538025115E-01, 3.3333244726700217E-01, 3.3333236160174562E-01, 3.3333226765435453E-01, 3.3333216462411397E-01,
 3.3333205163289703E-01, 3.3333192771768028E-01, 3.3333179182233691E-01, 3.3333164278863531E-01, 3.3333147934636848E-01,
 3.3333130010252821E-01, 3.3333110352943407E-01, 3.3333088795171317E-01, 3.3333065153202307E-01, 3.3333039225539285E-01,
 3.3333010791205198E-01, 3.3332979607859858E-01, 3.3332945409734754E-01, 3.3332907905368231E-01, 3.3332866775121789E-01,
 3.3332821668456303E-01, 3.3332772200944971E-01, 3.3332717950997570E-01, 3.3332658456268116E-01, 3.3332593209715394E-01,
 3.3332521655282665E-01, 3.3332443183159993E-01, 3.3332357124588702E-01, 3.3332262746163782E-01, 3.3332159243585835E-01,
 3.3332045734809296E-01, 3.3331921252528735E-01, 3.3331784735939274E-01, 3.3331635021701139E-01, 3.3331470834031435E-01,
 3.3331290773838990E-01, 3.3331093306809984E-01, 3.3330876750343147E-01, 3.3330639259223499E-01, 3.3330378809913208E-01,
 3.3330093183326115E-01, 3.3329779945939958E-01, 3.3329436429086118E-01, 3.3329059706241471E-01, 3.3328646568130038E-01,
 3.3328193495423791E-01, 3.3327696628811709E-01, 3.3327151736184257E-01, 3.3326554176656198E-01, 3.3325898861124292E-01,
 3.3325180209027605E-01, 3.3324392100946376E-01, 3.3323527826641108E-01, 3.3322580028095578E-01, 3.3321540637086305E-01,
 3.3320400806756095E-01, 3.3319150836619893E-01, 3.3317780090377680E-01, 3.3316276905850556E-01, 3.3314628496292592E-01,
 3.3312820842261187E-01, 3.3310838573153334E-01, 3.3308664837432950E-01, 3.3306281160484735E-01, 3.3303667288933342E-01,
 3.3300801020161019E-01, 3.3297658015642989E-01, 3.3294211596596252E-01, 3.3290432520304319E-01, 3.3286288735336311E-01,
 3.3281745113723987E-01, 3.3276763157993994E-01, 3.3271300680773824E-01, 3.3265311454499541E-01, 3.3258744828549625E-01,
 3.3251545310913466E-01, 3.3243652111274352E-01, 3.3234998642145908E-01, 3.3225511974449118E-01, 3.3215112243654416E-01,
 3.3203712002342323E-01, 3.3191215514759531E-01, 3.3177517988666355E-01, 3.3162504739493193E-01, 3.3146050281550671E-01,
 3.3128017340779092E-01, 3.3108255783285057E-01, 3.3086601453708164E-01, 3.3062874917300961E-01, 3.3036880099506832E-01,
 3.3008402816803128E-01, 3.2977209192663282E-01, 3.2943043952711226E-01, 3.2905628593525682E-01, 3.2864659420143683E-01,
 3.2819805448155892E-01, 3.2770706167439256E-01, 3.2716969166097204E-01, 3.2658167615152817E-01, 3.2593837617044985E-01,
 3.2523475424125153E-01, 3.2446534537236649E-01, 3.2362422699223736E-01, 3.2270498803996533E-01, 3.2170069748729319E-01,
 3.2060387265065043E-01, 3.1940644775020832E-01, 3.1809974328829449E-01, 3.1667443695406650E-01, 3.1512053691698017E-01,
 3.1342735855006870E-01, 3.1158350582706107E-01, 3.0957685886564834E-01, 3.0739456934362080E-01, 3.0502306579408162E-01,
 3.0244807108877880E-01, 2.9965463474085124E-01, 2.9662718299368623E-01, 2.9334959000189109E-01, 2.8980527374016318E-01,
 2.8597732057863656E-01, 2.8184864271491256E-01, 2.7740217282376473E-01, 2.7262110033732145E-01, 2.6748915365525916E-01,
 2.6199093225086456E-01, 2.5611229202094954E-01, 2.4984078625179609E-01, 2.4316616316149697E-01, 2.3608091904672340E-01,
 2.2858090352821242E-01, 2.2066597018001219E-01, 2.1234066189011228E-01, 2.0361491561152567E-01, 1.9450476575376355E-01,
 1.8503301943512301E-01, 1.7522987036157631E-01, 1.6513341153789624E-01, 1.5479000081337055E-01, 1.4425442805132027E-01,
 1.3358982929495919E-01, 1.2286729265963482E-01, 1.1216510391938465E-01, 1.0156758804109277E-01, 9.1163517349781439E-02,
 8.1044078440569670E-02, 7.1300418765016846E-02, 6.2020829604116284E-02, 5.3287663498961323E-02, 4.5174128276847310E-02,
 3.7741142385459128E-02, 3.1034471510695102E-02, 2.5082387570383009E-02, 1.9894090866834677E-02, 1.5459108033764636E-02,
 1.1747818401212730E-02, 8.7131697148918486E-03, 6.2935267560691038E-03, 4.4164659425391783E-03, 3.0032042187271880E-03,
 1.9732540528024818E-03, 1.2488502587001519E-03, 7.5871461878397356E-04, 4.4081499768156777E-04, 2.4392524575220451E-04,
 1.2797359725202879E-04, 6.3343077300475448E-05, 2.9419655025585636E-05, 1.2745365392302772E-05, 5.1169335599020947E-06,
 1.8901741251535226E-06, 6.3740959843300738E-07, 1.9454465428425228E-07, 5.3235307748363305E-08, 1.2925791306150350E-08,
 2.7532843110198589E-09, 5.0811401790392261E-10, 8.0138377804403672E-11, 1.0640562274610211E-11, 1.1699708776515703E-12,
 1.0462097780832971E-13, 7.4589950649813616E-15, 4.1486642985185821E-16, 1.7576585184704076E-17, 5.5256748401722395E-19,
 1.2525196093227465E-20, 1.9835664295728714E-22, 2.1201207776495630E-24, 1.4725180052948681E-26, 6.3750956542489416E-29,
 1.6437286717753684E-31, 2.4008472015154221E-34, 1.8804448508520476E-37, 7.4367468049092510E-41, 1.3901683296400136E-44,
 1.1425524663149736E-48, 3.8135576771636335E-53, 4.7382249921302625E-58, 1.9918674183610088E-63, 2.5514182629099523E-69,
 8.8776430462140520E-76, 7.3978039880032954E-83, 1.2858821359350879E-90, 4.0067768641364156E-99,1.8955044384018749E-108,
1.1346257060496344E-118,7.0370251809213380E-130,3.6320602408953402E-142,1.2267617615690636E-155,2.0832135805381402E-170,
1.3320744648091175E-186,2.3359346409749428E-204,7.9347779296212160E-224,3.5657946327740084E-245,1.3954704267022407E-268,
3.0065256255961158E-294,2.1738888417014848E-322, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0, 0.0, 0.0};

  for (q=0;q<501;q++) {
    exp3table[q].x = xtable[q];
    exp4table[q].x = xtable[q];
    exp3table[q].y = e3table[q];
    exp4table[q].y = e4table[q];
  }

  spline_pchip(num_datapoints,exp3table);
  spline_pchip(num_datapoints,exp4table);

  return;
}

/*================== end of exp_integral_setup() ============================*/

/*======================= machine_epsilon() ==================================*/

/*
 * Calculate machine's smallest meaningful floating-point number.
 * From Dennis and Schnabel (1996), Algorithm A1.3.1.
 */

FLOAT machine_epsilon(void)
{
  static int
    initialized = FALSE;
  static FLOAT
    eps;

  if (!initialized) {
    eps = 1.;
    while (1.+eps != 1.) {
      eps *= .5;
    }
    eps *= 2.;

    initialized = TRUE;
  }

  return eps;
}

/*======================= end of machine_epsilon() ===========================*/

/*======================= find_root() ========================================*/

/*
 * Adapted from zridder(), Numerical Recipes in C, 2nd ed., p.358.
 * Returns 0 if root is found,
 *        -1 if fabs(func(x1)) <  fabs(func(x2)) and zero is not bracketed,
 *         1 if fabs(func(x2)) <= fabs(func(x1)) and zero is not bracketed.
 */

#undef  MAX_IT
#define MAX_IT 100

#undef  UNLIKELY_VAL
#define UNLIKELY_VAL -1.11111e+30

int find_root(FLOAT  x1,
              FLOAT  x2,
              FLOAT  xacc,
              FLOAT *x_root,
              FLOAT  (*func)(FLOAT))
{
  register int
    iter,
    compare;
  register FLOAT
    fh,fl,fm,fnew,
    s,xh,xl,xm,xnew;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="find_root";

  fl = (*func)(x1);
  fh = (*func)(x2);
  if ((fl > 0. && fh < 0.) || (fl < 0. && fh > 0.)) {
    xl      = x1;
    xh      = x2;
    /* Set *x_root to an unlikely value: */
    *x_root = UNLIKELY_VAL;

    for (iter = 0; iter < MAX_IT; iter++) {
      xm = 0.5*(xl+xh);
      fm = func(xm);
      s  = sqrt(fm*fm-fl*fh);
      if (s == 0.) {
        return 0;
      }
      xnew = xm+(xm-xl)*((fl > fh ? 1. : -1.)*fm/s);

      if (fabs(xnew-*x_root) <= xacc) {
        return 0;
      }
      *x_root = xnew;

      fnew    = func(*x_root);
      if (fnew == 0.) {
        return 0;
      }

      if ((fnew > 0. ? fabs(fm) : -fabs(fm)) != fm) {
        xl = xm;
        fl = fm;
        xh = *x_root;
        fh = fnew;
      }
      else if ((fnew > 0. ? fabs(fl) : -fabs(fl)) != fl) {
        xh = *x_root;
        fh = fnew;
      }
      else if ((fnew > 0. ? fabs(fh) : -fabs(fh)) != fh) {
        xl = *x_root;
        fl = fnew;
      }
      else {
        sprintf(Message,"should never get here");
        util_error(dbmsname,Message);
      }
      if (fabs(xh-xl) <= xacc) {
        return 0;
      }
    }
    sprintf(Message,"exceeded MAX_IT = %d, current root calc = %e",MAX_IT,*x_root);
    util_error(dbmsname,Message);
 }
  else {
    if (fl == 0.) {
      *x_root = x1;
      return 0;
    }
    if (fh == 0.) {
      *x_root = x2;
      return 0;
    }

    /*
     * Return unsuccessfully, but assign closest value to x_root.
     */
    compare = fcmp(fabs(fl),fabs(fh));
    if (compare < 0) {
      *x_root = x1;
      return -1;
    }
    else {
      *x_root = x2;
      return 1;
    }
  }

  /* Should never get here. */
  sprintf(Message,"should never reach here");
  util_error(dbmsname,Message);
  return 2;
}

/*======================= end of find_root() =================================*/

/*======================= broyden_root() =====================================*/

/*
 * Finds a vector root, x[], of the vector function vecfunc().
 * This globally convergent algorithm is adapted from 
 * Numerical Recipes in C, 2nd ed., p. 389-392, which is based on
 * Dennis and Schnabel (1996).
 *
 * Call with an initial guess for x[].
 * Assumes zero-based indexing.
 *
 * The function global_step() performs the globally-convergent step. 
 * The argument step_type specifies the type of step taken. Currently, the
 * valid choices are DS_LINE_STEP and DS_DOGLEG_STEP, the latter being more
 * sophisticated and reliable (the "DS" stands for Dennis and Schnabel 1996).
 *
 * Returns 0 on normal execution and an error code if
 * the routine has failed, has converged to a local minimum, or can make
 * no further progress, in which case one should retry with a
 * different initial guess.
 */

int broyden_root(int    n,
                 FLOAT *x,
                 void  (*vecfunc)(int,FLOAT *,FLOAT *),
                 FLOAT  tol_f,
                 int    max_it)
{
  int
    k,j,i,
    restart,singular,skip,
    num_bytes,
    old_max_taken,
    it              = 0,
    max_taken       = FALSE,
    count_max_taken = 0,
    status          = DS_X_ACCEPTED;
  static int
    nold = 0;
  FLOAT
    sum,denom,
    f,f_old,
    max_step,
    test,h,
    tmp,
    delta = -1;
  static FLOAT
    *c,*d,
    *x_old,
    *fvec,*fvec2,*fvec_old,
    *g,*sn,*s,*t,*w,
    *qt,*r;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="broyden_root";

  if (n > nold) {
    if (nold != 0) {
      /* 
       * Free previously allocated memory: 
       */
      free_fvector(r,       0,nold*nold-1,dbmsname);
      free_fvector(qt,      0,nold*nold-1,dbmsname);
      free_fvector(w,       0,nold-1,dbmsname);
      free_fvector(t,       0,nold-1,dbmsname);
      free_fvector(s,       0,nold-1,dbmsname);
      free_fvector(sn,      0,nold-1,dbmsname);
      free_fvector(g,       0,nold-1,dbmsname);
      free_fvector(fvec_old,0,nold-1,dbmsname);
      free_fvector(fvec2,   0,nold-1,dbmsname);
      free_fvector(fvec,    0,nold-1,dbmsname);
      free_fvector(x_old,   0,nold-1,dbmsname);
      free_fvector(d,       0,nold-1,dbmsname);
      free_fvector(c,       0,nold-1,dbmsname);
    }
    /*
     * Allocate memory: 
     */
    c        = fvector(0,n-1,dbmsname);
    d        = fvector(0,n-1,dbmsname);
    x_old    = fvector(0,n-1,dbmsname);
    fvec     = fvector(0,n-1,dbmsname);
    fvec2    = fvector(0,n-1,dbmsname);
    fvec_old = fvector(0,n-1,dbmsname);
    g        = fvector(0,n-1,dbmsname);
    sn       = fvector(0,n-1,dbmsname);
    s        = fvector(0,n-1,dbmsname);
    t        = fvector(0,n-1,dbmsname);
    w        = fvector(0,n-1,dbmsname);
    qt       = fvector(0,n*n-1,dbmsname);
    r        = fvector(0,n*n-1,dbmsname);

    nold = n;
  }
  else {
    /* Clear working memory: */
    num_bytes = n*sizeof(FLOAT);
    memset(c,       0,num_bytes);
    memset(d,       0,num_bytes);
    memset(x_old,   0,num_bytes);
    memset(fvec,    0,num_bytes);
    memset(fvec2,   0,num_bytes);
    memset(fvec_old,0,num_bytes);
    memset(g,       0,num_bytes);
    memset(sn,      0,num_bytes);
    memset(s,       0,num_bytes);
    memset(t,       0,num_bytes);
    memset(w,       0,num_bytes);
    num_bytes *= n;
    memset(qt,      0,num_bytes);
    memset(r,       0,num_bytes);
  }

  /*
   * Calculate fvec[].
   */
  (*vecfunc)(n,x,fvec);

  f = 0.;
  for (i = 0; i < n; i++) {
    f += fvec[i]*fvec[i];
  }
  f *= 0.5;

  /*
   * Test if initial guess is a root.
   * NOTE: We do not compare to the more stringent 0.01*tol_f used in 
   *       NR and DS96 at iteration zero, so that our solutions are 
   *       consistent.
   */
  test = 0.;
  for (i = 0; i < n; i++) {
    test = MAX(test,fabs(fvec[i]));
  }
  if (test < tol_f) {
    /* initial x[] is a root. */
    return status;
  }

  /*
   * Calculate max_step for globally convergent step.
   */
  sum = 0.;
  for (i = 0; i < n; i++) {
    sum += x[i]*x[i];
  }
  max_step = DS_MAX_STEP*MAX(sqrt(sum),(FLOAT)n);

  /*
   * Main iteration loop.
   */
  restart = TRUE;
  for (it = 0; it < max_it; it++) {
    if (restart == TRUE) {
      /*
       * Compute forward-difference approximation to Jacobian.
       */
      for (j = 0; j < n; j++) {
        tmp  = x[j];
        h    = DS_SQRT_EPS*fabs(tmp);
        if (h == 0.) {
          h = DS_SQRT_EPS;
        }
        x[j] = tmp+h;

        (*vecfunc)(n,x,fvec2);

        h    = x[j]-tmp;
        x[j] = tmp;
        for (i = 0; i < n; i++) {
          R(i,j) = (fvec2[i]-fvec[i])/h;
        }
      }

      /*
       * Calculate QR decomposition.
       */
      singular = qr_decompose(n,r,c,d);
      if (singular) {
        status = DS_SINGULAR_JACOBIAN;
        return status;
      }


      /* Compute transpose, QT. */
      memset(qt,0,n*n*sizeof(FLOAT));
      for (i = 0; i < n; i++) {
        QT(i,i) = 1.;
      }
      for (k = 0; k < n; k++) {
        if (c[k]) {
          for (j = 0; j < n; j++) {
            sum = 0.;
            for (i = k; i < n; i++) {
              sum += R(i,k)*QT(i,j);
            }
            sum /= c[k];
            for (i = k; i < n; i++) {
              QT(i,j) -= sum*R(i,k);
            }
          }
        }
      }
      /* Form R explicitly. */
      for (i = 0; i < n; i++) {
        R(i,i) = d[i];
        for (j = 0; j < i; j++) {
          R(i,j) = 0.;
        }
      }
    }
    else {
      for (i = 0; i < n; i++) {
        s[i] = x[i]-x_old[i];
      }
      for (i = 0; i < n; i++) {
        sum = 0.;
        for (j = i; j < n; j++) {
          sum += R(i,j)*s[j];
        }
        t[i] = sum;
      }
      skip = TRUE;
      for (i = 0; i < n; i++) {
        sum = 0.;
        for (j = 0; j < n; j++) {
          sum += QT(j,i)*t[j];
        }
        w[i] = fvec[i]-fvec_old[i]-sum;
        if (fabs(w[i]) >= DS_EPS*(fabs(fvec[i])+fabs(fvec_old[i]))) {
          skip = FALSE;
        }
        else {
          w[i] = 0.;
        }
      }
      if (skip == FALSE) {
        for (i = 0; i < n; i++) {
          sum = 0.;
          for (j = 0; j < n; j++) {
            sum += QT(i,j)*w[j];
          }
          t[i] = sum;
        }

        denom = 0.;
        for (i = 0; i < n; i++) {
          denom += s[i]*s[i];
        }
        /* Store s/(s.s) in s. */
        for (i = 0; i < n; i++) {
          s[i] /= denom;
        }

        /* 
         * Update r and qt. 
         */
        qr_update(n,r,qt,t,s);

        for (i = 0; i < n; i++) {
          if (R(i,i) == 0.) {
            sprintf(Message,"R(%d,%d) singular",i,i);
            util_error(dbmsname,Message);
          }
          d[i] = R(i,i);
        }
      }
    }

    for (i = 0; i < n; i++) {
      sum = 0.;
      for (j = 0; j < n; j++) {
        sum += QT(i,j)*fvec[j];
      }
      g[i] = sum;
    }
    for (i = n-1; i >= 0; i--) {
      sum = 0.;
      for (j = 0; j <= i; j++) {
        sum += R(j,i)*g[j];
      }
      g[i] = sum;
    }

    /*
     * Store old x,fvec,f.
     */
    for (i = 0; i < n; i++) {
      /*
       * Screen for NaN: 
       */
      if (!isfinite(x[i]) || !isfinite(fvec[i])) {
        sprintf(Message,"**x[%d]=%g fvec[%d]=%g",i,x[i],i,fvec[i]);
        util_error(dbmsname,Message);
      }
      x_old[   i] = x[   i];
      fvec_old[i] = fvec[i];
    }
    f_old = f;

    /*
     * Compute right-hand side of linear equations, sn[].
     */
    for (i = 0; i < n; i++) {
      sum = 0.;
      for (j = 0; j < n; j++) {
        sum += QT(i,j)*fvec[j];
      }
      sn[i] = -sum;
    }

    /*
     * Solve R.x = sn.
     * See Numerical Recipes in C, 2nd ed., p. 100.
     */
    sn[n-1] /= d[n-1];
    for (i = n-2; i >= 0; i--) {
      sum = 0.;
      for (j = i+1; j < n; j++) {
        sum += R(i,j)*sn[j];
      }
      sn[i] = (sn[i]-sum)/d[i];
    }

    /*
     * Calculate new x,f,fvec[].
     */
    old_max_taken = max_taken;
    max_taken = global_step(n,x_old,f_old,g,r,sn,max_step,&delta,
                            DS_DOGLEG_STEP,&status,x,&f,fvec,vecfunc);

    /*
     * Screen for NaN: 
     */
    for (i = 0; i < n; i++) {
      if (!isfinite(x[i]) || !isfinite(fvec[i])) {
        sprintf(Message,"after global_step(): x[%d]=%g fvec[%d]=%g",i,x[i],i,fvec[i]);
        util_error(dbmsname,Message);
      }
    }

    /*
     * Screen for repeated maximum steps.
     */
    if (max_taken == TRUE) {
      if (old_max_taken == TRUE) {
        count_max_taken++;
      }
      else {
        count_max_taken = 1;
      }
    }
    else {
      count_max_taken = 0;
    }

    /*
     * Test for convergence.
     */
    test = 0.;
    for (i = 0; i < n; i++) {
      test = MAX(test,fabs(fvec[i]));
    }
    if (test < tol_f) {
      status = DS_X_ACCEPTED;
      return status;
    }

    if (count_max_taken >= 5) {
      status = DS_MAX_TAKEN_5;
      return status;
    }

    if (status == DS_X_NO_PROGRESS) {
      if (restart == TRUE) {
        return status;
      }
      else {
        test  = 0.;
        denom = MAX(f,0.5*(FLOAT)n);
        for (i = 0; i < n; i++) {
          tmp  = fabs(g[i])*MAX(fabs(x[i]),1.)/denom;
          test = MAX(test,tmp);
        }
        if (test < DS_TOL_MIN) {
          return status;
        }
        else {
          /*
           * Try reinitializing the Jacobian.
           */
          restart = TRUE;
        }
      }
    }
    else {
      restart = FALSE;
      test    = 0.;
      for (i = 0; i < n; i++) {
        tmp  = (fabs(x[i]-x_old[i]))/MAX(fabs(x[i]),1.);
        test = MAX(test,tmp);
        if (test < DS_TOL_X) {
          /* 
           * Convergence. 
           */
          return status;
        }
      }
    }
  }

  status = DS_MAX_IT_EXCEEDED;
  return status;
}

/*======================= end of broyden_root() ==============================*/

/*======================= global_step() ======================================*/


/*
 * Take a globally convergent step towards a vector root.
 * Returns max_taken.
 *
 * The function line_search() backtracks along the quasi-Newton direction.
 *
 * The function dogleg_driver() uses the model trust-region approach, where
 * delta is the radius of the trust region. A step is taken 
 * in the steepest descent direction for small delta, in the quasi_Newton 
 * direction for large delta, and on the connecting line segment for  
 * intermediate delta.
 * See Dennis and Schnabel (1996, DS96).
 *
 * The value of step_type can be DS_LINE_STEP, DS_HOOK_STEP, or DS_DOGLEG_STEP.
 *
 * NOTE: Unlike in DS96, here R is R of QR, not the transpose of R. 
 */

int global_step(int    n,
                FLOAT *x_old,
                FLOAT  f_old,
                FLOAT *g,
                FLOAT *r,
                FLOAT *sn,
                FLOAT  max_step,
                FLOAT *delta,
                int    step_type,
                int   *status,
                FLOAT *x,
                FLOAT *f,
                FLOAT *fvec,
                void  (*vecfunc)(int,FLOAT *,FLOAT *))
{
  int
    max_taken = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="global_step";

  if (step_type == DS_LINE_STEP) {
    max_taken = line_search(n,x_old,f_old,g,sn,max_step,
                            status,x,f,fvec,vecfunc);
  }
  else if (step_type == DS_HOOK_STEP) {
    sprintf(Message,"DS_HOOK_STEP not yet implemented");
    util_error(dbmsname,Message);
  }
  else if (step_type == DS_DOGLEG_STEP) {
    max_taken = dogleg_driver(n,x_old,f_old,g,r,sn,max_step,delta,
                              status,x,f,fvec,vecfunc);
  }
  else {
    sprintf(Message,"unrecognized step_type=%d",step_type);
    util_error(dbmsname,Message);
  }

  return max_taken;
}

/*======================= end of global_step() ===============================*/

/*======================= line_search() ======================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p. 385-386, which is
 * an implementation of Algorithm A6.3.1 of Dennis and Schnabel (1996).
 * Returns max_taken.
 * Assumes zero-based indexing.
 */

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
                void  (*vecfunc)(int,FLOAT *,FLOAT *))
{
  int
    i,
    max_taken = FALSE;
  FLOAT
    a,b,
    lambda,lambda_prev,lambda_min,
    disc,f_prev,
    rhs1,rhs2,
    initial_slope,
    newt_length,
    rel_step_length,
    tmp,
    tmp_lambda; 
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="line_search";

  *status = DS_X_ACCEPTED;

  /* Screen for bad input values. */
  if (n < 1) {
    sprintf(Message,"n=%d < 1",n);
    util_error(dbmsname,Message);
  }

  tmp = 0.;
  for (i = 0; i < n; i++) {
    tmp += sn[i]*sn[i];
  }
  newt_length = sqrt(tmp);

  if (newt_length > max_step) {
    tmp = max_step/newt_length;
    for (i = 0; i < n; i++) {
      sn[i] *= tmp;
    }
  }

  initial_slope = 0.;
  for (i = 0; i < n; i++) {
    initial_slope += g[i]*sn[i];
  }

  rel_step_length = 0.;
  for (i = 0; i < n; i++) {
    rel_step_length = MAX( rel_step_length,
                           fabs(sn[i])/MAX(fabs(x_old[i]),1.) );
  }

  lambda_min = DS_TOL_X/rel_step_length;
  lambda     = 1.;

  /*
   * Iteration loop.
   */
  while (TRUE) {
    for (i = 0; i < n; i++) {
      x[i] = x_old[i]+lambda*sn[i];
    }

    /*
     * Calculate fvec[].
     */
    (*vecfunc)(n,x,fvec);

    *f = 0.;
    for (i = 0; i < n; i++) {
      *f += fvec[i]*fvec[i];
    }
    *f *= .5;

    if (lambda < lambda_min) {
      /* 
       * Convergence on dx. 
       * For zero finding, calling program should verify the convergence.
       */
      for (i = 0; i < n; i++) {
        x[i] = x_old[i];
      }
      *status = DS_X_NO_PROGRESS;
      return max_taken;
    }
    else if (*f <= f_old+DS_ALPHA*lambda*initial_slope) {
      /* Sufficient function decrease. */
      if (lambda == 1. && newt_length > 0.99*max_step) {
        max_taken = TRUE;
      }
      return max_taken;
    }
    else {
      if (lambda == 1.) {
        /* First backtrack uses a quadratic model. */
        tmp_lambda = -initial_slope/(2.*(*f-f_old-initial_slope));
      }
      else {
        /* Subsequent backtracks use a cubic model. */
        rhs1 =     *f-initial_slope*lambda     -f_old;
        rhs2 = f_prev-initial_slope*lambda_prev-f_old;
        a    = (rhs1/(lambda*lambda)-rhs2/(lambda_prev*lambda_prev))
               /(lambda-lambda_prev);
        b    = ( -lambda_prev*rhs1/(lambda*lambda)
                 +lambda*rhs2/(lambda_prev*lambda_prev) )/(lambda-lambda_prev);
        if (a == 0.) {
          tmp_lambda = -initial_slope/(2.*b);
        }
        else {
          disc = b*b-3.*a*initial_slope;
          if (disc < 0.) {
            sprintf(Message,"roundoff problem");
            util_error(dbmsname,Message);
          }
          else {
            tmp_lambda = (-b+sqrt(disc))/(3.*a);
          }
        }
        if (tmp_lambda > .5*lambda) {
          tmp_lambda = .5*lambda;
        }
      }
    }
    lambda_prev = lambda;
    f_prev      = *f;
    lambda      = MAX(tmp_lambda,.1*lambda);
  }

  /* Never get here. */
}

/*======================= end of line_search() ===============================*/

/*======================= dogleg_driver() ====================================*/

/*
 * Adapted from Dennis and Schnabel (1996), Appendix A, Algorithm A6.4.3.
 * Returns max_taken;
 * Assumes zero-based indexing.
 */
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
                  void   (*vecfunc)(int,FLOAT *,FLOAT *))
{
  int
    i,
    max_taken,
    newt_taken,
    first_dog = TRUE;
  FLOAT
    newt_length,
    f_prev,
    tmp,
   *s,
   *s_hat,
   *nu_hat,
   *x_prev;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="dogleg_driver";

  /*
   * Allocate memory.
   */
  s      = fvector(0,n-1,dbmsname);
  s_hat  = fvector(0,n-1,dbmsname);
  nu_hat = fvector(0,n-1,dbmsname);
  x_prev = fvector(0,n-1,dbmsname);

  *status = DS_INITIAL;

  tmp = 0.;
  for (i = 0; i < n; i++) {
    tmp += sn[i]*sn[i];
  }
  newt_length = sqrt(tmp);

  while (*status >= DS_REDUCE_DELTA) {
    /*
     * Find new step.
     */
    newt_taken = dogleg_step(n,g,r,sn,newt_length,max_step,
                             delta,&first_dog,s_hat,nu_hat,s);
    /*
     * Check new point and update trust region.
     */
    max_taken = trust_region(n,x_old,f_old,g,s,newt_taken,max_step,
                             DS_DOGLEG_STEP,r,delta,status,x_prev,&f_prev,
                             x,f,fvec,vecfunc);
  }

  /*
   * Free allocated memory.
   */
  free_fvector(x_prev,0,n-1,dbmsname);
  free_fvector(nu_hat,0,n-1,dbmsname);
  free_fvector(s_hat, 0,n-1,dbmsname);
  free_fvector(s,     0,n-1,dbmsname);

  return max_taken;
}

/*======================= end of dogleg_driver() =============================*/

/*======================= dogleg_step() ======================================*/

/*
 * Adapted from Dennis and Schnabel (1996), Appendix A, Algorithm A6.4.4.
 * Returns newt_taken.
 * Assumes zero-based indexing.
 */
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
                FLOAT *s)
{
  int
    i,j,
    newt_taken;
  static int
    eta_warned = FALSE;
  FLOAT
    alpha,beta,al_be,
    lambda,
    tmp,tmp_nu,tmp_cauchy;
  static FLOAT
    eta,
    cauchy_length;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="dogleg_step";

  if (newt_length <= *delta) {
    /*
     * s is Newton step.
     */
    newt_taken = TRUE;
    for (i = 0; i < n; i++) {
      s[i] = sn[i];
    }
    *delta = newt_length;
  }
  else {
    /*
     * Newton step is too long, find s on double-dogleg curve.
     */
    newt_taken = FALSE;
    if (*first_dog == TRUE) {
      /*
       * Calculate double-dogleg curve.
       */
      *first_dog = FALSE;

      alpha      = 0.;
      for (i = 0; i < n; i++) {
        alpha += g[i]*g[i];
      }

      beta = 0.;
      for (i = 0; i < n; i++) {
        tmp = 0.;
        for (j = i; j < n; j++) {
          /*
           * NOTE: Unlike DS96, here R is R of QR, not the transpose of R.
           */
          tmp += R(i,j)*g[j];
        }
        beta += tmp*tmp;
      }
      al_be         = alpha/beta;
      cauchy_length = al_be*sqrt(alpha);

      tmp = 0.;
      for (i = 0; i < n; i++) {
        s_hat[i] = -g[i]*al_be;
        tmp     +=  g[i]*sn[i];
      }
      eta = 0.2+0.8*alpha*al_be/fabs(tmp);

      /*
       * Check range of eta, which should be [0,1]:
       */
      if (eta > 1.) {
        /* 
         * If eta > 1.0 print a one-time warning. 
         * This behavior is known to be caused by optimization (-O) 
         * for gcc (at least versions 2.7.2.3.f.1 and 2.96).
         */
        if (!eta_warned) {
          fprintf(stderr,"eta=%g > 1, setting eta=1. This warning will not be repeated.\n"
                          "Check for optimization error (-O vs -g).",eta);
          eta_warned = TRUE;
        }
        eta = 1.;
      }

      for (i = 0; i < n; i++) {
        nu_hat[i] = eta*sn[i]-s_hat[i];
      }

      if (*delta == -1) {
        /*
         * First iteration, and no initial trust region was
         * provided by the user.
         */
        *delta = MIN(cauchy_length,max_step);
      }
    }

    if (eta*newt_length <= *delta) {
      /*
       * Take partial step in Newton direction.
       */
      for (i = 0; i < n; i++) {
        s[i] = ((*delta)/newt_length)*sn[i];
      }
    }
    else if ((cauchy_length) >= (*delta)) {
      /*
       * Take step in steepest descent direction.
       */
      for (i = 0; i < n; i++) {
        s[i] = ((*delta)/(cauchy_length))*s_hat[i];
        /*
         * Screen for NaN.
         */
        if (!isfinite(s[i])) {
          sprintf(Message,"s[%d]=%g,*delta=%g,cauchy_length=%g,s_hat=%g",i,s[i],*delta,cauchy_length,s_hat[i]);
          util_error(dbmsname,Message);
        }
      }
    }
    else {
      /*
       * Take convex-combination step.
       */
      tmp    = 0.;
      tmp_nu = 0.;
      for (i = 0; i < n; i++) {
        tmp    += nu_hat[i]*s_hat[i];
        tmp_nu += nu_hat[i]*nu_hat[i];
      }
      tmp_cauchy = cauchy_length*cauchy_length-(*delta)*(*delta);
      lambda     = (-tmp+sqrt(tmp*tmp-tmp_nu*tmp_cauchy))/tmp_nu;
      for (i = 0; i < n; i++) {
        s[i] = s_hat[i]+lambda*nu_hat[i];
        /*
         * Screen for NaN.
         */
        if (!isfinite(s[i])) {
          sprintf(Message,"s[%d]=%g,lambda=%g,nu_hat=%g,s_hat=%g",i,s[i],lambda,nu_hat[i],s_hat[i]);
          util_error(dbmsname,Message);
        }
      }
    }
  }


  return newt_taken;
}


/*======================= end of dogleg_step() ===============================*/

/*======================= trust_region() =====================================*/

/*
 * Adapted from Dennis and Schnabel (1996), Appendix A, Algorithm A6.4.5.
 * Returns max_taken.
 * Assumes zero-based indexing.
 */

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
                 void  (*vecfunc)(int,FLOAT *,FLOAT *))
{
  int
    i,j,
    max_taken = FALSE;
  FLOAT
    initial_slope,
    step_length,
    rel_step_length,
    delta_tmp,
    tmp;
  FLOAT
     df,
     df_tol,
     df_pred;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="trust_region";

  /*
   * Screen for NaN.
   */
  for (i = 0; i < n; i++) {
    if (!isfinite(s[i])) {
      sprintf(Message,"s[%d]=%g",i,s[i]);
      util_error(dbmsname,Message);
    }
  }

  tmp = 0.;
  for (i = 0; i < n; i++) {
    tmp += s[i]*s[i];
  }
  step_length = sqrt(tmp);

  /* Take step. */
  for (i = 0; i < n; i++) {
    x[i] = x_old[i]+s[i];
  }

  /* Calculate fvec[]. */
  (*vecfunc)(n,x,fvec);

  /* Compute f. */
  *f = 0.;
  for (i = 0; i < n; i++) {
    *f += fvec[i]*fvec[i];
  }
  *f *= .5;

  df = *f-f_old;

  initial_slope = 0.;
  for (i = 0; i < n; i++) {
    initial_slope += g[i]*s[i];
  }

  if (*status != DS_INCREASE_DELTA) {
    *f_prev = 0.;
  }

  df_tol = DS_ALPHA*initial_slope;

  if (*status == DS_INCREASE_DELTA && (*f >= *f_prev || df > df_tol)) {
    /*
     * Retreat.
     */
    *status = DS_X_ACCEPTED;
    for (i = 0; i < n; i++) {
      x[i] = x_prev[i];
    }
    *f      = *f_prev;
    *delta *= .5;
  }
  else if (df >= df_tol) {
    /*
     * The value of f is too large.
     */
    rel_step_length = 0.;
    for (i = 0; i < n; i++) {
      rel_step_length = MAX( rel_step_length,
                             fabs(s[i])/MAX(fabs(x[i]),1.) );
    }

    if (rel_step_length < DS_TOL_X) {
      /*
       * The step is too small.
       */
      *status = DS_X_NO_PROGRESS;
      for (i = 0; i < n; i++) {
        x[i] = x_old[i];
      }
    }
    else {
      /*
       * Reduce delta.
       */
      *status   = DS_REDUCE_DELTA;
      delta_tmp = (-initial_slope*step_length)/(2.*(df-initial_slope));
      if (delta_tmp < (*delta)*0.1) {
        *delta *= 0.1;
      }
      else if (delta_tmp > (*delta)*0.5) {
        *delta *= 0.5;
      }
      else {
        *delta = delta_tmp;
      }
    }
  }
  else {
    /*
     * The value of f is sufficiently small.
     */
    df_pred = initial_slope;

    if (step_type == DS_HOOK_STEP) {
      sprintf(Message,"DS_HOOK_STEP not yet implemented");
      util_error(dbmsname,Message);
    }
    else if (step_type == DS_DOGLEG_STEP) {
      for (i = 0; i < n; i++) {
        tmp = 0.;
        for (j = i; j < n; j++) {
          /*
           * NOTE: Unlike DS96, here R is R of QR, not the transpose of R.
           */
          tmp += R(i,j)*s[j];
        }
        df_pred += tmp*tmp*.5;
      }
    }
    else {
      sprintf(Message,"unrecognized step_type=%d",step_type);
      util_error(dbmsname,Message);
    }

    if ( ((*status) != DS_REDUCE_DELTA && fabs(df_pred-df) <= 0.1*fabs(df)) ||
         (df <= initial_slope && newt_taken == FALSE && (*delta) <= 0.99*max_step) ) {
      /*
       * Double delta.
       */
      *status = DS_INCREASE_DELTA;
      for (i = 0; i < n; i++) {
        x_prev[i] = x[i];
      }
      *f_prev = *f;
      *delta  = MIN((*delta)*2.,max_step);
    }
    else {
      /*
       * Accept x, choose delta for next iteration.
       */
      *status = DS_X_ACCEPTED;
      if (step_length > 0.99*max_step) {
        max_taken = TRUE;
      }
      if (df >= 0.1*df_pred) {
        /*
         * Decrease delta.
         */
        *delta *= 0.5;
      }
      else if (df <= 0.75*df_pred) {
        /*
         * Increase delta.
         */
        *delta = MIN((*delta)*2.,max_step);
      }
      else {
        /*
         * Leave delta unchanged.
         */
        ;
      }
    }
  }

  return max_taken;
}

/*======================= end of trust_region() ==============================*/

/*======================= qr_decompose() =====================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p.99.
 * Returns FALSE if normal and TRUE if singular.
 * Assumes zero-based indexing.
 */

int qr_decompose(int    n,
                 FLOAT *r,
                 FLOAT *c,
                 FLOAT *d)
{
  int
    i,j,k,
    singular;
  FLOAT
    sigma,sum,tau,scale;

  singular = FALSE;

  for (k = 0; k < n-1; k++) {
    /* 
     * Put scale=0. inside k loop as in Dennis & Schnabel (1996), 
     * Algorithm A3.2.1, p.305, rather than outside the loop as in 
     * Numerical Recipes' qrdcmp(), p. 99.
     */
    scale = 0.;
    for (i = k; i < n; i++) {
      scale = MAX(scale,fabs(R(i,k)));
    }
    if (scale == 0.) {
      /*
       * Singular case.
       */
      c[k]     = 0.;
      d[k]     = 0.;
      singular = TRUE;
    }
    else {
      for (i = k; i < n; i++) {
        R(i,k) /= scale;
      }
      sum = 0.;
      for (i = k; i < n; i++) {
        sum += R(i,k)*R(i,k);
      }
      sigma = NR_SIGN(sqrt(sum),R(k,k)); 
      R(k,k) += sigma;
      c[k]    = sigma*R(k,k);
      d[k]    = -scale*sigma;
      for (j = k+1; j < n; j++) {
        sum = 0.;
        for (i = k; i < n; i++) {
          sum += R(i,k)*R(i,j);
        }
        tau = sum/c[k];
        for (i = k; i < n; i++) {
          R(i,j) -= tau*R(i,k);
        }
      }
    }
  }
  d[n-1] = R(n-1,n-1);

  if (d[n-1] == 0.) {
    singular = TRUE;
  }

  return singular;
}

/*======================= end of qr_decompose() ==============================*/

/*======================= qr_update() ========================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p. 101.
 * Assumes zero-based indexing.
 */

void qr_update(int    n,
               FLOAT *r,
               FLOAT *qt,
               FLOAT *u,
               FLOAT *v)
{
  int
    i,j,k;
  FLOAT
    tmp;

  /* Find largest k such that u[k] != 0. */
  for (k = n-1; k > 0; k--) {
    if (u[k]) {
      break;
    }
  }

  for (i = k-1; i >= 0; i--) {
    qr_rotate(n,r,qt,i,u[i],-u[i+1]);
    if (u[i] == 0.) {
      u[i] = fabs(u[i+1]);
    }
    else if (fabs(u[i]) > fabs(u[i+1])) {
      tmp  = u[i+1]/u[i];
      u[i] = fabs(u[i])*sqrt(1.+tmp*tmp);
    }
    else {
      tmp  = u[i]/u[i+1];
      u[i] = fabs(u[i+1])*sqrt(1.+tmp*tmp);
    }
  }

  for (j = 0; j < n; j++) {
    R(0,j) += u[0]*v[j];
  }
  for (i = 0; i < k; i++) {
    qr_rotate(n,r,qt,i,R(i,i),-R(i+1,i));
  }

  return;
}

/*======================= end of qr_update() =================================*/

/*======================= qr_rotate() ========================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p. 101.
 * Assumes zero-based indexing.
 */

void qr_rotate(int    n,
               FLOAT *r,
               FLOAT *qt,
               int    i,
               FLOAT  a,
               FLOAT  b)
{
  int
    j;
  FLOAT
    c,factor,
    s,w,y;

  if (a == 0.) {
    c = 0.;
    s = b > 0. ? 1. : -1.;
  }
  else if (fabs(a) > fabs(b)) {
    factor = b/a;
    c      = NR_SIGN(1./sqrt(1.+factor*factor),a);
    s      = factor*c;
  }
  else {
    factor = a/b;
    s      = NR_SIGN(1./sqrt(1.+factor*factor),b);
    c      = factor*s;
  }

  for (j = 0; j < n; j++) {
    y        = R(i,  j);
    w        = R(i+1,j);
    R(i,  j) = c*y-s*w;
    R(i+1,j) = s*y+c*w;
  }

  for (j = 0; j < n; j++) {
    y         = QT(i,  j);
    w         = QT(i+1,j);
    QT(  i,j) = c*y-s*w;
    QT(i+1,j) = s*y+c*w;
  }

  return;
}

/*======================= end of qr_rotate() =================================*/

/*======================= lu_decompose() =====================================*/

/*
 * Adapted from Numerical Recipes in C, pp. 46-47.
 * Assumes zero-based indexing, with (i,j) referring to (row,column).
 */

#undef  A
#define A(i,j) a[(j)+n*(i)]

void lu_decompose(int   n,
                 FLOAT *a,
                 int   *index,
                 FLOAT *d)
{
  int
    i,imax,j,k;
  static int
    n_max = 0;
  FLOAT
    tiny = 1.e-20,
    big,dum,sum,temp;
  static FLOAT
   *vv;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lu_decompose";

  /*
   * Allocate memory.
   */
  if (n_max == 0) {
    n_max = n;
    vv    = fvector(0,n_max-1,dbmsname);
  }
  else if (n > n_max) {
    free_fvector(vv,0,n_max-1,dbmsname);
    n_max = n;
    vv    = fvector(0,n_max-1,dbmsname);
  }

  *d = 1.;
  for (i = 0; i < n; i++) {
    big = 0.;
    for (j = 0; j < n; j++) {
      temp = fabs(A(i,j));
      if (temp > big) {
        big = temp;
      }
    }
    if (big == 0.) {
      sprintf(Message,"matrix singular");
      util_error(dbmsname,Message);
    }
    vv[i] = 1./big;
  }

  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = A(i,j);
      for (k = 0; k < i; k++) {
        sum -= A(i,k)*A(k,j);
      }
      A(i,j) = sum;
    }
    big = 0.;
    for (i = j; i < n; i++) {
      sum = A(i,j);
      for (k = 0; k < j; k++) {
        sum -= A(i,k)*A(k,j);
      }
      A(i,j) = sum;
 
      dum = vv[i]*fabs(sum);
      if (dum >= big) {
        big  = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum       = A(imax,k);
        A(imax,k) = A(j,   k);
        A(j,   k) = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    index[j] = imax;

    if (A(j,j) == 0.) {
      A(j,j) = tiny;
      fprintf(stderr,"Warning, %s(), matrix is singular\n",dbmsname);
    }

    if (j != n-1) {
      dum = 1./A(j,j);
      for (i = j+1; i < n; i++) {
        A(i,j) *= dum;
      }
    }
  }
 
  return;
}

/*======================= end of lu_decompose() ==============================*/

/*======================= lu_backsub() =======================================*/

/*
 * Adapted from Numerical Recipes in C, p. 47.
 * Assumes zero-based indexing, with (i,j) referring to (row,column).
 */

#undef  A
#define A(i,j) a[(j)+n*(i)]

void lu_backsub(int    n,
                FLOAT *a,
                int   *index,
                FLOAT *b)
{
  int
    ii = -1,
    i,ip,j; 
  FLOAT
    sum;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lu_backsub";

  for (i = 0; i < n; i++) {
    ip    = index[i];
    sum   = b[ip];
    b[ip] = b[i];
    if (ii != -1) {
      for (j = ii; j <= i-1; j++) {
        sum -= A(i,j)*b[j];
      }
    }
    else if (sum != 0.) {
      ii = i;
    }
    b[i] = sum;
  } 

  for (i = n-1; i >= 0; i--) {
    sum = b[i];
    for (j = i+1; j < n; j++) {
      sum -= A(i,j)*b[j];
    }
    b[i] = sum/A(i,i);
  }

  return;
}

/*======================= end of lu_backsub() ================================*/

/*======================= lu_improve() =======================================*/

/*
 * Adapted from Numerical Recipes in C, p. 56.
 * Assumes zero-based indexing.
 *
 * NOTE: May not have much effect if lu_decompose(), lu_backsub()
 *       are in double precision.
 */

void lu_improve(int    n,
                FLOAT *a,
                FLOAT *alu,
                int   *index,
                FLOAT *b,
                FLOAT *x)
{
  int
    j,i;
  static int
    n_max = 0;
  static FLOAT
   *r;
  double
    sdp;  /* NOTE: sdp must be double precision. */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lu_improve";

  /*
   * Allocate memory.
   */
  if (n_max == 0) {
    n_max = n;
    r     = fvector(0,n_max-1,dbmsname);
  }
  else if (n > n_max) {
    free_fvector(r,0,n_max-1,dbmsname);
    n_max = n;
    r     = fvector(0,n_max-1,dbmsname);
  }

  for (i = 0 ; i < n; i++) {
    sdp = -b[i];
    for (j = 0; j < n; j++) {
      sdp += A(i,j)*x[j];
    }
    r[i] = sdp;
  }
  lu_backsub(n,alu,index,r);
  for (i = 0; i < n; i++) {
    x[i] -= r[i];
  }

  return;
}

/*======================= end of lu_improve() ================================*/

/*======================= find_place_in_table() ==============================*/

/*
 * Find place in table using bisection.
 * Adapted from Numerical Recipes in C, p. 117.
 *
 * NOTE: Finds place relative to the x component of the float_triplet table.
 */
int find_place_in_table(int            n,
                        float_triplet *table,
                        FLOAT          x,
                        FLOAT         *dx)
{
  register int
    il,im,iu,
    ascend,
    cmplo,cmphi;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="find_place_in_table";

  il     = -1;
  iu     =  n;
  ascend = (table[n-1].x >= table[0].x);
  cmplo  = fcmp(x,table[0  ].x);
  cmphi  = fcmp(x,table[n-1].x);

  /*
   * Check x against table endpoints.
   */
  if (ascend) {
    if (cmplo <= 0) {
      il  = 0;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
    if (cmphi >= 0) {
      il  = n-2;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
  }
  else {
    if (cmplo >= 0) {
      il  = 0;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
    if (cmphi <= 0) {
      il  = n-2;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
  }

  /*
   * Use bisection to search table.
   */
  while (iu-il > 1) {
    im = (iu+il) >> 1;
    if (fcmp(x,table[im].x) >= 0 == ascend) {
      il = im;
    }
    else {
      iu = im;
    }
  }

  *dx = table[il+1].x-table[il].x;
  return il;
}

/*======================= end of find_place_in_table() =======================*/

/*======================= hunt_place_in_table() ==============================*/

/*
 * Find place in table using bisection.
 * Adapted from Numerical Recipes in C, p. 118-119.
 *
 * NOTE: Hunts place relative to the x component of the float_triplet table.
 *       Requires an initial guess (il) to location within table, which should be
 *       populated by the calling function as the value returned previously.
 */
int hunt_place_in_table(int            n,
                        float_triplet *table,
                        FLOAT          x,
                        FLOAT         *dx,
                        int            il)
{
  register int
    im,iu,inc,
    ascend,
    cmplo,cmphi;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="hunt_place_in_table";

  ascend = (table[n-1].x >= table[0].x);
  cmplo  = fcmp(x,table[0  ].x);
  cmphi  = fcmp(x,table[n-1].x);

  /*
   * Check x against table endpoints.
   */
  if (ascend) {
    if (cmplo <= 0) {
      il  = 0;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
    if (cmphi >= 0) {
      il  = n-2;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
  }
  else {
    if (cmplo >= 0) {
      il  = 0;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
    if (cmphi <= 0) {
      il  = n-2;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
  }

  /*
   * Begin the hunt based on initial guess...
   */

  /* jhi = iu    *jlo = il     jm = im  */

  if (il <= -1 || il > n-1) {
    /* 
     * Our guess wasn't useful here since it's off the end
     * of the table. Skip straight to bisection.
     */
    il = -1;
    iu = n;
  }
  else {
    inc = 1;
    if (fcmp(x,table[il].x) >= 0 == ascend){
      /* We're hunting up... */
      if (il == n-1) {
        il = n-2;
        return il;
      }
      iu = il + 1;
      while (fcmp(x,table[iu].x) >= 0 == ascend) {
        /* Not done hunting... */
        il = iu;
        inc += inc;
        /* so double the increment */
        iu = il + inc;
        /* and set the upper bound inc above the lower bound */
        if (iu > n-1){
          /* Oops, fell off the end of the table... */
          iu = n;
          break;
        }
      }
    }
    else {
      /* We're hunting down... */
      if (il == 0){
        il = 0;
        return il;
      }
      /* decrement lower bound by one and set the upper bound to that */
      iu = il--;
      while (fcmp(x,table[il].x) <= 0 == ascend) {
        /* Not done hunting... */
        iu = il;
        inc <<= 1;
        /* bitwise shift and assign to double the increment */
        if (inc >= iu){
          il = -1;
          break;
        }
        else {
          il = iu-inc;
        }
      }
    }
  }

  /*
   * Hunt is done, so use bisection to search table.
   */
  while (iu-il > 1) {
    im = (iu+il) >> 1;
    if (fcmp(x,table[im].x) >= 0 == ascend) {
      il = im;
    }
    else {
      iu = im;
    }
  }

  *dx = table[il+1].x-table[il].x;
  return il;
}

/*======================= end of hunt_place_in_table() ======================*/

/*======================= tridiag() =========================================*/

/*
 * Solve a tridiagnonal matrix system.
 * Adapted from Numerical Recipes in C, 2nd ed., pp. 51-54.
 * If pivot_type = WITH_PIVOTING, use band_decomp() and band_back_sub().
 * Assumes zero-based indexing.
 */

#undef  AA
#define AA(i,j) aa[(m1+m2+1)*(i)+(j)]
#undef  AAORIG
#define AAORIG(i,j) aaorig[(m1+m2+1)*(i)+(j)]
#undef  AAL
#define AAL(i,j) aal[m1*(i)+(j)]

void tridiag(int    n,
             FLOAT *a,
             FLOAT *b,
             FLOAT *c,
             FLOAT *r,
             FLOAT *u,
             int    pivot_type)
{
  int
    j,
    m1,m2,mm,
    *index;
  FLOAT
    bet,
    *gam,
    *aa,
    *aaorig,
    *aal,
     d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="tridiag";

  /*
   * Check validity of n.
   */
  if (n <= 0) {
    sprintf(Message,"n = %d",n);
    util_error(dbmsname,Message);
  }

  if (pivot_type == WITHOUT_PIVOTING) {
    /* Allocate memory. */
    gam = fvector(0,n-1,dbmsname);

    if (b[0] == 0.0) {
      sprintf(Message,"b[0] = 0\n"
                     "Rewrite equations as a set of order n-1, with u[1] trivially eliminated");
      util_error(dbmsname,Message);
    }
    bet  = b[0];
    u[0] = r[0]/bet;
    for (j = 1; j < n; j++) {
      gam[j] = c[j-1]/bet;
      bet    = b[j]-a[j]*gam[j];
      if (bet == 0.) {
        /* 
         * Encountered a zero pivot. 
         * Try again using pivot_type = WITH_PIVOTING.
         */
        fprintf(stderr,"Warning: tridiag(): retrying with pivoting.\n");
        /* Free allocated memory. */
        free_fvector(gam,0,n-1,dbmsname);
        tridiag(n,a,b,c,r,u,WITH_PIVOTING);
        return;
      }
      u[j] =(r[j]-a[j]*u[j-1])/bet;
    }

    /* Backsubstitution: */
    for (j = n-2; j >= 0; j--) {
      u[j] -= gam[j+1]*u[j+1];
    }
    /* Free allocated memory. */
    free_fvector(gam,0,n-1,dbmsname);
    return;
  }
  else if (pivot_type == WITH_PIVOTING) {
    /*
     * Use band_decomp() and band_back_sub().
     */
    m1 = 1;
    m2 = 1;
    mm = m1+m2+1;
    /*
     * Allocate memory.
     */
    aa     = fvector(0,n*mm-1,dbmsname);
    aaorig = fvector(0,n*mm-1,dbmsname);
    aal    = fvector(0,n*m1-1,dbmsname);
    index  = ivector(0,n-1,dbmsname);

    /*
     * Load matrix AA and keep copy AAORIG.
     */
    for (j = 0; j < n; j++) {
      AA(j,m1+1) = AAORIG(j,m1+1) = c[j];
      AA(j,m1  ) = AAORIG(j,m1  ) = b[j];
      AA(j,m1-1) = AAORIG(j,m1-1) = a[j];
    }
    
    band_decomp(n,m1,m2,aa,aal,index,&d);

    /* 
     * Since tridiag() does not overwrite the input rhs vector, r,
     * with the answer, u, but band_back_sub() does, copy r into u
     * before calling band_back_sub().
     */
    for (j = 0; j < n; j++) {
      u[j] = r[j];
    }

    band_back_sub(n,m1,m2,aa,aal,index,u);

    /*
     *  Reduce roundoff errors with call to band_improve().
     */
    band_improve(n,m1,m2,aaorig,aa,aal,index,r,u);

    /*
     * Free allocated memory.
     */
    free_fvector(aa,    0,n*mm-1,dbmsname);
    free_fvector(aaorig,0,n*mm-1,dbmsname);
    free_fvector(aal,   0,n*m1-1,       dbmsname);
    free_ivector(index, 0,n-1,          dbmsname);

    return;
  }
  else {
    sprintf(Message,"unrecognized pivot_type=%d",pivot_type);
    util_error(dbmsname,Message);
  }
}

/*======================= end of tridiag() ===================================*/

/*======================= band_decomp() ======================================*/

/*
 * Decompose a banded matrix.
 * Adapted from Numerical Recipes in C, 2nd ed., p. 53.
 * The input matrix must be stored in compact form, as described on p. 52.
 * Assumes zero-based indexing.
 */

#undef  SWAP
#define SWAP(a,b) {tmp=(a);(a)=(b);(b)=tmp;}
#undef  TINY
#define TINY 1.e-20

#undef  A
#define A(i,j) a[(m1+m2+1)*(i)+(j)]
#undef  AL
#define AL(i,j) al[m1*(i)+(j)]

void band_decomp(int     n,
                 int     m1,
                 int     m2,
                 FLOAT  *a,
                 FLOAT  *al,
                 int    *index,
                 FLOAT  *d)
{
  int
    i,j,k,l,
    mm;
  FLOAT
    tmp;
  static int
    warned = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="band_decomp";

  mm = m1+m2+1;
  l  = m1;
  for (i = 0; i < m1; i++) {
    for (j = m1-i; j < mm; j++) {
      A(i,j-l) = A(i,j);
    }
    l--;
    for (j = mm-l-1; j < mm; j++) {
      A(i,j) = 0.;
    }
  }
  *d = 1.;
  l  = m1;
  for (k = 0; k < n; k++) {
    tmp = A(k,0);
    i   = k;
    if (l < n) {
      l++;
    }
    for (j = k+1; j < l; j++) {
      if (fabs(A(j,0)) > fabs(tmp)) {
        tmp = A(j,0);
        i   = j;
      }
    }
    index[k] = i;
    if (tmp == 0.) {
      /*
       * Matrix is algorithmically singular.
       */
      A(k,0) = TINY;
      if (!warned) {
        fprintf(stderr,"**warning: %s, matrix is algorithmically singular (future warnings suppressed)\n",dbmsname);
        warned = TRUE;
      }
    }
    if (i != k) {
      *d = -(*d);
      for (j = 0; j < mm; j++) {
        SWAP(A(k,j),A(i,j))
      }
    }
    for (i = k+1; i < l; i++) {
      tmp          = A(i,0)/A(k,0);
      AL(k,i-k-1) = tmp;
      for (j = 1; j < mm; j++) {
        A(i,j-1) = A(i,j)-tmp*A(k,j);
      }
      A(i,mm-1) = 0.;
    }
  }

  return;
}

/*======================= end of band_decomp() ===============================*/

/*======================= band_back_sub() ====================================*/

/*
 * Back substitute for a banded matrix using the output from band_decomp().
 * Adapted from Numerical Recipes in C, 2nd ed., p. 54.
 * Assumes zero-based indexing.
 */

#undef  SWAP
#define SWAP(a,b) {tmp=(a);(a)=(b);(b)=tmp;}

#undef  A
#define A(i,j) a[(m1+m2+1)*(i)+(j)]
#undef  AL
#define AL(i,j) al[m1*(i)+(j)]

void band_back_sub(int     n,
                   int     m1,
                   int     m2,
                   FLOAT  *a,
                   FLOAT  *al,
                   int    *index,
                   FLOAT  *b)
{
  int
    i,k,l,
    mm;
  FLOAT
    tmp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="band_back_sub";

  mm = m1+m2+1;
  l  = m1;
  for (k = 0; k < n; k++) {
    i = index[k];
    if (i != k) {
      SWAP(b[k],b[i]);
    }
    if (l < n) {
      l++;
    }
    for (i = k+1; i < l; i++) {
      b[i] -= AL(k,i-k-1)*b[k];
    }
  }
  l = 1;
  for (i = n-1; i >= 0; i--) {
    tmp = b[i];
    for (k = 1; k < l; k++) {
      tmp -= A(i,k)*b[k+i];
    }
    b[i] = tmp/A(i,0);
    if (l < mm) {
      l++;
    }
  }

  return;
}

/*======================= end of band_back_sub() =============================*/

/*======================= band_multiply() ====================================*/
/*
 *  Compute matrix muliplication b = A.x, where A is in the compact-storage
 *  form of a band-diagonal matrix. 
 *  Based on Numerical Recipes in C, banmul(), p. 52.
 *  Assumes zero-based indexing.
 */

#undef  A
#define A(i,j) a[(m1+m2+1)*(i)+(j)]

void band_multiply(int     n,
                   int     m1,
                   int     m2,
                   FLOAT  *a,
                   FLOAT  *x,
                   FLOAT  *b)
{
  int
    i,j,k,tmploop;

  for (i = 0; i < n; i++) {
    k       = i-m1;
    tmploop = IMIN(m1+m2+1,n-k);
    b[i]    = 0.;
    for (j = IMAX(0,-k); j < tmploop; j++) {
      b[i] += A(i,j)*x[j+k];
    }
  }

  return;
}
                   
/*======================= end of band_multiply() =============================*/

/*======================= band_improve() =====================================*/
/*
 * Based on Numerical Recipes in C, Secion 2.5, Iterative Improvement of
 * a Solution to Linear Equations.
 * This is for the band-diagnonal matrix case, and is analogous to lu_improve().
 * AORIG is the original matrix in compact form, whereas A and AL are the
 * matrices returned from band_decomp().
 * NOTE: The functionality of band_multiply() is echoed here because of the
 *       requirement of double precision.
 * Assumes zero-based indexing.
 */

#undef  AORIG
#define AORIG(i,j) aorig[(m1+m2+1)*(i)+(j)]

void band_improve(int     n,
                  int     m1,
                  int     m2,
                  FLOAT  *aorig,
                  FLOAT  *a,
                  FLOAT  *al,
                  int    *index,
                  FLOAT  *b,
                  FLOAT  *x)
{
  int
    k,j,i,tmploop;
  static int
    n_max = 0;
  static FLOAT
   *r;
  double
    sdp;  /* NOTE: sdp must be double precision. */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="band_improve";

  /*
   * Allocate memory.
   */
  if (n_max == 0) {
    n_max = n;
    r     = fvector(0,n_max-1,dbmsname);
  }
  else if (n > n_max) {
    free_fvector(r,0,n_max-1,dbmsname);
    n_max = n;
    r     = fvector(0,n_max-1,dbmsname);
  }

  /*
   * The band-diagonal indexing is as in band_multiply().
   */
  for (i = 0; i < n; i++) {
    k       = i-m1;
    tmploop = IMIN(m1+m2+1,n-k);
    sdp     = -b[i];
    for (j = IMAX(0,-k); j < tmploop; j++) {
      sdp += AORIG(i,j)*x[j+k];
    }
    r[i] = sdp;
  }

  band_back_sub(n,m1,m2,a,al,index,r);

  for (i = 0; i < n; i++) {
    x[i] -= r[i];
  }

  return;
}

/*======================= end of band_improve() ==============================*/

/*======================= poly_interp() ======================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd Ed., pp. 109-110.
 * Assumes zero-offset arrays.
 */

#undef  N_MAX
#define N_MAX 10

FLOAT poly_interp(int     n,
                  FLOAT  *xa,
                  FLOAT  *ya,
                  FLOAT   x,
                  FLOAT  *dy)
{
  int
    i,m,ns;
  static int
    initialized=0;
  FLOAT
    den,dif,dift,
    ho,hp,w,
    y;
  static FLOAT
    *c,
    *d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="poly_interp";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    c = fvector(0,N_MAX-1,dbmsname);
    d = fvector(0,N_MAX-1,dbmsname);
    initialized = 1;
  }

  dif = FLOAT_MAX, ns = INT_MAX;
  for (i = 0; i < n; i++) {
    dift = fabs(x-xa[i]);
    if (dift < dif) {
      ns  = i;
      dif = dift;
    }
  }
  y = ya[ns--];

  if (dif == 0.) {
    *dy = 0.;
  }
  else {
    if (n > N_MAX) {
      /*
       * Some systems have trouble with too many calls to calloc(),
       * so we allocate memory once for this function, assuming
       * n <= N_MAX.  Screen for this limit.
       */
      sprintf(Message,"n=%d > N_MAX=%d, increase N_MAX",n,N_MAX);
      util_error(dbmsname,Message);
    }

    for (i = 0; i < n; i++) {
      c[i] = ya[i];
      d[i] = ya[i];
    }

    for (m = 0; m <= n-2; m++) {
      for (i = 0; i <= n-2-m; i++) {
        ho  = xa[i    ]-x;
        hp  = xa[i+m+1]-x;
        w   = c[i+1]-d[i];
        den = ho-hp;
        if (den == 0.) {
          sprintf(Message,"xa[%d]=%g, xa[%d]=%g",i,xa[i],i+m+1,xa[i+m+1]);
          util_error(dbmsname,Message);
        }
        den  = w/den;
        d[i] = hp*den;
        c[i] = ho*den;
      }
      *dy  = (2*ns < n-m-3) ? c[ns+1] : d[ns--];
      y   += *dy;
    }
  }

  return y;
}

/*======================= end of poly_interp() ===============================*/

/*======================= poly_coeff() =======================================*/

/*
 * Compute polynomial coefficients.
 * Assumes zero-offset arrays [0..n].
 * Adapted from polcoe() in Numerical Recipes in C, 2nd ed, p. 121.
 * Due to G.B. Rybicki.
 */
void poly_coeff(int    n,
                FLOAT *x,
                FLOAT *y,
                FLOAT *coeff)
{
  register int
    k,j,i;
  static int
    nmax = 0;
  register FLOAT
    phi,ff,b;
  static FLOAT
   *s = NULL;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="poly_coeff";

  if (n > nmax) {
    /* Allocate working array. */
    if (s) free_fvector(s,0,nmax,dbmsname);
    s    = fvector(0,n,dbmsname);
    nmax = n;
  }
  
  /* Zero arrays s and coeff */
  memset(s,    0,(n+1)*sizeof(FLOAT));
  memset(coeff,0,(n+1)*sizeof(FLOAT));

  s[n] = -x[0];
  for (i = 1; i <= n; i++) {
    for (j = n-i; j <= n-1; j++) {
      s[j] -= x[i]*s[j+1];
    }
    s[n] -= x[i];
  }

  for (j = 0; j <= n; j++) {
    phi = n+1;
    for (k = n; k >= 1; k--) {
      phi = (double)k*s[k]+x[j]*phi;
    }
    ff = y[j]/phi;
    b  = 1.;
    for (k = n; k >= 0; k--) {
      coeff[k] += b*ff;
      b         = s[k]+x[j]*b;
    }
  }

  return;
}

/*======================= end of poly_coeff() ================================*/

/*======================= nth_trapezoidal() ==================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd Ed., p. 137.
 * Usage note: Call with successively increasing n, starting at n = 1.
 */

FLOAT nth_trapezoidal(int   n,
                      FLOAT (*func)(FLOAT),
                      FLOAT a,
                      FLOAT b)
{
  register int
    it,j;
  register FLOAT
    x,tnm,sum,dx;
  static FLOAT
    s;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="nth_trapezoidal";

  if (n < 1) {
    sprintf(Message,"n=%d",n);
    util_error(dbmsname,Message);
  }
  else if (n == 1) {
    s = 0.5*(b-a)*(func(a)+func(b));
    return s;
  }
  else {
    it  = (1 << (n-2));
    tnm = (FLOAT)it;
    dx  = (b-a)/tnm;
    x   = a+0.5*dx;
    sum = 0.;
    for (j = 0; j < it; j++) {
      sum += func(x);
      x   += dx;
    }
    s = 0.5*(s+(b-a)*sum/tnm);
    return s;
  }

  /* Never get here. */
  sprintf(Message,"should never get here");
  util_error(dbmsname,Message);
  return FLOAT_MAX;
}

/*======================= end of nth_trapezoidal() ===========================*/

/*======================= romberg_integral() =================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd Ed., p. 140.
 * The argument tol is the fractional accuracy (tolerance) desired.
 */
#undef  MAX_IT
#define MAX_IT  30
#undef  NUM_PTS
#define NUM_PTS 5

FLOAT romberg_integral(FLOAT (*func)(FLOAT),
                       FLOAT a,
                       FLOAT b,
                       FLOAT tol)
{
  int
    it,
    n;
  FLOAT
    ss,
    dss,
    s[MAX_IT+1],
    h[MAX_IT+1];

  if (a == b) {
    return 0.;
  }

  h[0]   = 1.;
  for (it = 0; it < MAX_IT; it++) {
    n = it+1;
    s[it] = nth_trapezoidal(n,func,a,b);
    if (n >= NUM_PTS) {
      ss = poly_interp(NUM_PTS,h+(n-NUM_PTS),s+(n-NUM_PTS),0.,&dss);
      if (fabs(ss) < tol && fabs(dss) < tol) {
        return ss;
      }
      else if (fabs(dss) <= tol*fabs(ss)) {
        return ss;
      }
    }
    s[n] = s[it];
    h[n] = 0.25*h[it];
  }

  fprintf(stderr,"Warning: romberg_integral(): reached it=%d, a,b,tol=%g %g %g ans=%g\n",
                  MAX_IT,a,b,tol,ss);
  return ss;
}

/*======================= end of romberg_integral() ==========================*/

/*======================= compact_integration() ==============================*/

/*
 * 3rd-order accurate integration on a non-staggered grid
 * (endpoints are 2nd-order accurate).
 *
 * Input:
 * z[0 to n+1]: the independent-variable 
 * dfdz[0 to n+1]: the derivative
 * f[0]: the value of the integral at the starting position
 *
 * Output:
 * f[0 to n+1]: the integral of df/dz
 *
 * This scheme is based on the one described and illustrated in Fig. 8(d) of
 *    Leslie LM,  Purser RJ, 1992, A comparative study of the performance of various
 *       vertical discretization schemes, Meteorol. Atmos. Phys. 50, 61-73.
 */
void compact_integration(int    n,
                         FLOAT *z,
                         FLOAT *dfdz,
                         FLOAT *f)
{
  int
    k;
  FLOAT
    dz30,dz31,dz32,dz20,dz21,dz10,
    df0,df1,df2,df3,
    df;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="compact_integration";

  /*  Centered step from start. */
  k      = 0;
  df     = .5*(dfdz[k]+dfdz[k+1])*(z[k+1]-z[k]);
  f[k+1] = f[k]+df;

  for (k = 1; k < n; k++) {
    dz10   = z[k  ]-z[k-1];
    dz21   = z[k+1]-z[k  ];
    dz32   = z[k+2]-z[k+1];
    dz20   = dz21+dz10;
    dz30   = dz32+dz20;
    dz31   = dz32+dz21;
    df0    = dfdz[k-1];
    df1    = dfdz[k  ];
    df2    = dfdz[k+1];
    df3    = dfdz[k+2];
    df     = -(df0*(dz21*dz21*dz21*(dz31+dz32))/(dz10*dz20*dz30)
              +df1*(dz21*(dz21*(dz21+2.*(dz10-dz31))-6.*dz10*dz31))/(dz10*dz31)
              +df2*(dz21*(dz21*(dz21+2.*(dz32-dz20))-6.*dz20*dz32))/(dz20*dz32)
              +df3*(dz21*dz21*dz21*(dz10+dz20))/(dz30*dz31*dz32)               )/12.;
    f[k+1] = f[k]+df;
  }

  /*  Centered step to end. */
  k      = n;
  df     = .5*(dfdz[k]+dfdz[k+1])*(z[k+1]-z[k]);
  f[k+1] = f[k]+df;

  return;
}

/*======================= end of compact_integration() =======================*/

/*======================= compact_differentiation() ==========================*/

#undef  ACOMPACT
#define ACOMPACT(i,j) acompact[(m1+m2+1)*(i)+(j)]
#undef  AORIG
#define AORIG(i,j) aorig[(m1+m2+1)*(i)+(j)]
#undef  AL
#define AL(i,j) al[m1*(i)+(j)]

/*
 * 2nd-order accurate compact differentiation on a non-staggered grid.
 * See Leslie and Purser (1992) as cited above.
 *
 * Input:
 * action: SETUP_UTIL, RUN_UTIL, and CLEANUP_UTIL
 * z[0 to n+1]: the independent-variable 
 * f[0 to n+1]: the dependent variable
 * dfdz[0], the derivative df/dz at the starting position.
 *          Set to 1.e+20 to calculate this internally assuming a parabolic fit.
 *          
 * NOTE: Call with SETUP_UTIL and CLEANUP_UTIL only once for a given z[] vector.
 *
 * Output:
 * dfdz[0 to n+1], df/dz
 *
 * This is the exact inverse of the 3rd-order accurate integration scheme
 * used in compact_integration(); it involves a 4-band matrix inversion.
 */
void compact_differentiation(int    action,
                             int    n,
                             FLOAT *z,
                             FLOAT *f,
                             FLOAT *dfdz)
{
  const int
    m1 = 2,
    m2 = 1;
  FLOAT
    denom;
  static FLOAT
    d,
    dz10,dz20,dz21,dz30,dz31,dz32,
    term1,term2,
   *acompact,
   *aorig,
   *al,
   *b,
   *borig;
  int
    i;
  static int
   *index;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="compact_differentiation";

  if (action == SETUP_UTIL) {
    /* Allocate memory. */
    acompact = fvector(0,n*(m1+m2+1)-1,dbmsname);
    aorig    = fvector(0,n*(m1+m2+1)-1,dbmsname);
    al       = fvector(0,n*m1-1,dbmsname);
    b        = fvector(0,n-1,dbmsname);
    borig    = fvector(0,n-1,dbmsname);
    index    = ivector(0,n-1,dbmsname);

    for (i = 0; i < n-1; i++) {
      dz10 = z[i+1]-z[i  ];
      dz21 = z[i+2]-z[i+1];
      dz32 = z[i+3]-z[i+2];
      dz20 = dz21+dz10;
      dz30 = dz32+dz20;
      dz31 = dz32+dz21;

      ACOMPACT(i,0) = AORIG(i,0) = -(dz21*dz21*dz21*(dz31+dz32))/(dz10*dz20*dz30)/12.;
      ACOMPACT(i,1) = AORIG(i,1) = -(dz21*(dz21*(dz21+2.*(dz10-dz31))-6.*dz10*dz31))/(dz10*dz31)/12.;
      ACOMPACT(i,2) = AORIG(i,2) = -(dz21*(dz21*(dz21+2.*(dz32-dz20))-6.*dz20*dz32))/(dz20*dz32)/12.;
      ACOMPACT(i,3) = AORIG(i,3) = -(dz21*dz21*dz21*(dz10+dz20))/(dz30*dz31*dz32)/12.;
    }
    /* 
     * The values ACOMPACT(0,0), ACOMPACT(0,1), and ACOMPACT(1,0) are needed because
     * they are used below for b[0] and b[1].
     */

    i = n-1;
    ACOMPACT(i,0) = AORIG(i,0) = 0.;
    ACOMPACT(i,1) = AORIG(i,1) = .5*(z[i+2]-z[i+1]);
    ACOMPACT(i,2) = AORIG(i,2) = .5*(z[i+2]-z[i+1]);

    band_decomp(n,m1,m2,acompact,al,index,&d);

    /*
     * Set up terms for internal calculation of derivative at starting point.
     */
    denom = (z[1]-z[0])*(z[2]-z[0])*(z[2]-z[1]);
    term1 =  SQR(z[2]-z[0])/denom;
    term2 = -SQR(z[1]-z[0])/denom;
  }
  else if (action == RUN_UTIL) {
    if (dfdz[0] > .99e+20) {
      /*
       * Calculate the derivative at the starting point assuming a parabolic fit.
       */
      dfdz[0] = term1*(f[1]-f[0])+term2*(f[2]-f[0]);
    }

    /*
     * dfdz[1] comes from a simple centered step.
     */
    dfdz[1] = (f[1]-f[0])/(.5*(z[1]-z[0]))-dfdz[0];

    /*
     * Load the right-hand side vector.
     */
    b[0] = borig[0] = f[2]-f[1]-AORIG(0,0)*dfdz[0]-AORIG(0,1)*dfdz[1];
    b[1] = borig[1] = f[3]-f[2]-AORIG(1,0)*dfdz[1];

    for (i = 2; i < n; i++) {
      b[i] = borig[i] = f[i+2]-f[i+1];
    }

    band_back_sub(n,m1,m2,acompact,al,index,b);
    band_improve(n,m1,m2,aorig,acompact,al,index,borig,b);

    for (i = 0; i < n; i++) {
      dfdz[i+2] = b[i];
    }
  }
  else if (action == CLEANUP_UTIL) {
    /* Free allocated memory */
    free(index);
    free(borig);
    free(b);
    free(al);
    free(aorig);
    free(acompact);
  }
  else {
    sprintf(Message,"action=%d not recognized",action);
    util_error(dbmsname,Message);
  }

  return;
}

/*======================= end of compact_differentiation() ===================*/

/*======================= crank_nicolson() ===================================*/

/*
 * Apply diffusion to the 1D input vector, A[k], for one timestep using the Crank-Nicolson scheme.
 *
 * The equation integrated in time is:
 *
 *      dA     1   d   /     dA  \
 *      --- = --- --- (  mu ----  )
 *      dt    rho  dz  \     dz  /
 *
 * The algorithm does not assume a regularly spaced grid, but
 * it does assume that the grid position, z[k], does not change during the timestep.
 * See Blottner (1980, Computer and Fluids 8, 421-434) for a description of how the variable grid
 * is handled.
 *
 * For input A[k], k = 0 and k = n+1 are boundary points.  Set to 1.e+20
 * to signal a no-flux boundary, otherwise set to the boundary value.
 *
 * One may specify the case rho[k] = 1.0 for all k by inputting NULL for rho.
 *
 * The input vector mu[k] may vary spatially, but is assumed to not change during the timestep.
 * The position of mu[k] is centered between A[k] and A[k+1].
 *
 * Input vectors are assumed to be continguous 1D vectors with memory
 * allocations in the range [0,n+1]. The grid positions z are expected to be monotonically
 * increasing, and an error is generated if this is not the case.
 * The result is returned in ANS[k] for k = 1 to n, such that 
 * ANS[0] and ANS[n+1] are not altered. Memory for ANS should not overlap A.
 *
 * NOTE: Assumes no MPI domain decomposition in active dimension.
 */

void crank_nicolson(int    n,
                    FLOAT  dt,
                    FLOAT *z,
                    FLOAT *A,
                    FLOAT *mu,
                    FLOAT *rho,
                    FLOAT *ANS)
{
  int
    k,kstart,kend;
  FLOAT
    coeff1,coeff3,factor;
  static int
    nold = -1;
  static FLOAT
   *a,
   *b,
   *c,
   *r;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="crank_nicolson";

  /*
   * Screen for invalid n.
   */
  if (n < 1) {
    sprintf(Message,"called with n=%d < 1",n);
    util_error(dbmsname,Message);
  }

  /* 
   * Screen for non-monotonically increasing z.
   */
  if (A[0] > .99e+20) {
    /* z[1]-z[0] is not used for no-flux b.c. */
    kstart = 2;
  }
  else {
    kstart = 1;
  }
  if (A[n+1] > .99e+20) {
    /* z[n+1]-z[n] is not used for no-flux b.c. */
    kend = n;
  }
  else {
    kend = n+1;
  }
  for (k = kstart; k <= kend; k++) {
    if (z[k] <= z[k-1]) {
      sprintf(Message,"z[%d]=%g <= z[%d]=%g",k,z[k],k-1,z[k-1]);
      util_error(dbmsname,Message);
    }
  }

  /* Allocate memory for local vectors. */
  if (nold == -1) {
    a    = fvector(0,n-1,dbmsname);
    b    = fvector(0,n-1,dbmsname);
    c    = fvector(0,n-1,dbmsname);
    r    = fvector(0,n-1,dbmsname);
    nold = n;
  }
  else if (nold < n) {
    free_fvector(a,0,nold-1,dbmsname);
    free_fvector(b,0,nold-1,dbmsname);
    free_fvector(c,0,nold-1,dbmsname);
    free_fvector(r,0,nold-1,dbmsname);
    a    = fvector(0,n-1,dbmsname);
    b    = fvector(0,n-1,dbmsname);
    c    = fvector(0,n-1,dbmsname);
    r    = fvector(0,n-1,dbmsname);
    nold = n;
  }
  
  if (A[0] > .99e+20) {
    /* No-flux boundary condition. */
    coeff3 = 0.;
  }
  else {
    coeff3 = mu[0]/(z[1]-z[0]);
  }
  for (k = 1; k <= n; k++) {
    coeff1 = coeff3;
    if (k == 1) {
      coeff3 = mu[k]/(z[k+1]-z[k]);
      if (rho) {
        factor = dt/(rho[k]*(z[k+1]-z[k-1]));
      }
      else {
        factor = dt/(z[k+1]-z[k-1]);
      }
      b[k-1] = 1.+(coeff1+coeff3)*factor;
      c[k-1] = -coeff3*factor;
      r[k-1] = 2.*A[k-1]*coeff1*factor
                 +A[k  ]*(1.-(coeff1+coeff3)*factor)
                 +A[k+1]*coeff3*factor;
    }
    else if (k == n) {
      if (A[n+1] > .99e+20) {
        /* No-flux boundary condition. */
        coeff3 = 0.;
      }
      else {
        coeff3 = mu[k]/(z[k+1]-z[k]);
      }
      if (rho) {
        factor = dt/(rho[k]*(z[k+1]-z[k-1]));
      }
      else {
        factor = dt/(z[k+1]-z[k-1]);
      }
      a[k-1] = -coeff1*factor;
      b[k-1] = 1.+(coeff1+coeff3)*factor;
      r[k-1] =    A[k-1]*coeff1*factor
                 +A[k  ]*(1.-(coeff1+coeff3)*factor)
              +2.*A[k+1]*coeff3*factor;
    }
    else {
      coeff3 = mu[k]/(z[k+1]-z[k]);
      if (rho) {
        factor = dt/(rho[k]*(z[k+1]-z[k-1]));
      }
      else {
        factor = dt/(z[k+1]-z[k-1]);
      }
      a[k-1] = -coeff1*factor;
      b[k-1] = 1.+(coeff1+coeff3)*factor;
      c[k-1] = -coeff3*factor;
      r[k-1] = A[k-1]*coeff1*factor
              +A[k  ]*(1.-(coeff1+coeff3)*factor)
              +A[k+1]*coeff3*factor;
    }
  }

  tridiag(n,a,b,c,r,ANS+1,WITH_PIVOTING);

  return;
}

/*======================= end of crank_nicolson() ============================*/

/*======================= hqr() ==============================================*/

/* 
 * From Numerical Recipes in C, 2nd Ed., p. 491-492.
 * Finds eigenvalues of a real Hessenberg matrix. 
 * Assumes zero-based indexing.
 */

void hqr(FLOAT *a_eig, 
         int    n, 
         FLOAT *wr, 
         FLOAT *wi)
{
  int 
    nn,m,l,k,j,
    klen,its,i,mmin;
  FLOAT 
    z,y,x,w,v,u,
    t,s,r,q,p,
    anorm;

  /* Need klen for A_EIG() macro. */
  klen  = n;

  anorm = fabs(A_EIG(0,0));
  for (i = 1; i < n; i++) {
    for (j= i-1; j < n; j++) {
      anorm += fabs(A_EIG(i,j));
    }
  }

  nn = n-1;
  t  = 0.;
  while (nn >= 0) {
    its = 0;
    do {
      for (l = nn; l >= 1; l--) {
	s = fabs(A_EIG(l-1,l-1))+fabs(A_EIG(l,l));
	if (s == 0.0) {
          s = anorm;
        }
	if ((fabs(A_EIG(l,l-1))+s) == s) {
          break;
        }
      }
      x = A_EIG(nn,nn);
      if (l == nn) {
	wr[nn  ] = x+t;
	wi[nn--] = 0.0;
      } 
      else {
	y = A_EIG(nn-1,nn-1);
	w = A_EIG(nn,nn-1)*A_EIG(nn-1,nn);
	if (l == nn-1) {
	  p = 0.5*(y-x);
	  q = p*p+w;
	  z = sqrt(fabs(q));
	  x += t;
	  if (q >= 0.) {
	    z        = p+NR_SIGN(z,p);
	    wr[nn-1] = wr[nn]=x+z;
	    if (z) {
              wr[nn] = x-w/z;
            }
	    wi[nn-1] = wi[nn] = 0.0;
	  } 
          else {
	    wr[nn-1] =   wr[nn] = x+p;
	    wi[nn-1] = -(wi[nn] = z);
	  }
	  nn -= 2;
	} 
        else {
	  if (its == 30){
	    fprintf(stderr, "Warning: hqr(): too many iterations\n");
	    return;
	  }
	  if (its == 10 || its == 20) {
	    t += x;
	    for (i = 0; i < nn; i++) {
              A_EIG(i,i) -= x;
            }
	    s = fabs(A_EIG(nn,nn-1))+fabs(A_EIG(nn-1,nn-2));
	    y = x = 0.75*s;
	    w = -0.4375*s*s;
	  }
	  its++;
	  for (m = nn-2; m >= 0; m--) {
	    z = A_EIG(m,m);
	    r = x-z;
	    s = y-z;
	    p = (r*s-w)/A_EIG(m+1,m)+A_EIG(m,m+1);
	    q = A_EIG(m+1,m+1)-z-r-s;
	    r = A_EIG(m+2,m+1);
	    s = fabs(p)+fabs(q)+fabs(r);
	    p /= s;
	    q /= s;
	    r /= s;
	    if (m == l) {
              break;
            }
	    u = fabs(A_EIG(m,m-1))*(fabs(q)+fabs(r));
	    v = fabs(p)*(fabs(A_EIG(m-1,m-1))+fabs(z)+fabs(A_EIG(m+1,m+1)));
	    if ((u+v) == v) {
              break;
            }
	  }
	  for (i = m+2; i <= nn; i++) {
	    A_EIG(i,i-2) = 0.;
	    if (i != (m+2)) {
              A_EIG(i,i-3) = 0.;
            }
	  }
	  for (k = m; k <= nn-1; k++) {
	    if (k != m) {
	      p = A_EIG(k,  k-1);
	      q = A_EIG(k+1,k-1);
	      r = 0.0;
	      if (k != nn-1) {
                r = A_EIG(k+2,k-1);
              }
	      if ((x = fabs(p)+fabs(q)+fabs(r)) != 0.) {
		p /= x;
		q /= x;
		r /= x;
	      }
	    }
	    if ((s = NR_SIGN(sqrt(p*p+q*q+r*r),p)) != 0.) {
	      if (k == m) {
		if (l != m) {
		  A_EIG(k,k-1) = -A_EIG(k,k-1);
                }
	      } 
              else {
		A_EIG(k,k-1) = -s*x;
              }
	      p += s;
	      x  = p/s;
	      y  = q/s;
	      z  = r/s;
	      q /= p;
	      r /= p;
	      for (j = k; j <= nn; j++) {
		p = A_EIG(k,j)+q*A_EIG(k+1,j);
		if (k != nn-1) {
		  p            += r*A_EIG(k+2,j);
		  A_EIG(k+2,j) -= p*z;
		}
		A_EIG(k+1,j) -= p*y;
		A_EIG(k,  j) -= p*x;
	      }
	      mmin = (nn < k+3) ? nn : k+3;
	      for (i = l; i <= mmin; i++) {
		p = x*A_EIG(i,k)+y*A_EIG(i,k+1);
		if (k != nn-1) {
		  p            += z*A_EIG(i,k+2);
		  A_EIG(i,k+2) -= p*r;
		}
		A_EIG(i,k+1) -= p*q;
		A_EIG(i,k  ) -= p;
	      }
	    }
	  }
	}
      }
    } while (l < nn-1);
  }
}

/*======================= end of hqr() =======================================*/

/*======================= quicksort() ========================================*/
/* 
 * From K & R 2nd ed., p. 87: 
 */
void quicksort(FLOAT *mag, 
               int    left, 
               int    right)
{
  int 
    i,last;
  void 
    swap(FLOAT *mag,int i,int j);

  if (left >= right) {
    return;
  }

  swap(mag,left,(left+right)/2);
  last = left;
  for (i = left+1; i <= right; i++) {
    if (mag[i] > mag[left]) {
      swap(mag, ++last, i);
    }
  }
  swap(mag,left,last);

  quicksort(mag,left,  last-1);
  quicksort(mag,last+1,right);

  return;
}

/*====================== end of quicksort() =================================*/

/*====================== swap() =============================================*/

void swap(FLOAT *mag,
          int    i,
          int    j)
{
  FLOAT 
    temp;

  temp   = mag[i];
  mag[i] = mag[j];
  mag[j] = temp;

  return;
}

/*======================= end of swap() ======================================*/

/*======================= four1() ============================================*/

      /*
       *  FFT routine.
       *  Numerical Recipes in C, 2nd ed, p. 507.
       *  Assumes data length is a power of two.
       */

#undef  SWAP
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(FLOAT         data[], 
           unsigned long nn, 
           int           isign)
{
  unsigned long 
    n,mmax,m,j,istep,i;
  FLOAT 
    tempr,tempi;
  FLOAT 
    wtemp,wr,wpr,wpi,wi,theta;

  n = nn << 1;
  j = 1;
  for (i = 1; i < n; i+= 2) {
    if (j > i) {
       SWAP(data[j  ],data[i  ]);
       SWAP(data[j+1],data[i+1]);
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    theta = isign*(2*M_PI/mmax);
    wtemp = sin(0.5*theta);
    wpr   = -2.0*wtemp*wtemp;
    wpi   = sin(theta);
    wr    = 1.0;
    wi    = 0.0;
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
        j = i + mmax;
        tempr      = wr*data[j]-wi*data[j+1];
        tempi      = wr*data[j+1]+wi*data[j];
        data[j  ]  = data[i  ]-tempr;
        data[j+1]  = data[i+1] - tempi;
        data[i  ] += tempr;
        data[i+1] += tempi;
      }
      wr = (wtemp = wr)*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
  }
}
      
/*======================= end of four1() ====================================*/

/*======================= realft() ==========================================*/

      /*
       *  Real FFT routine.
       *  Numerical Recipes in C, 2nd ed, p. 513.
       *  Assumes data length is a power of 2; assumes unit-based array.
       *  Result of inverse must be multiplied by 2/n.
       */
void realft(FLOAT data[], unsigned long n, int isign)
{
  void 
    four1(FLOAT data[], unsigned long nn, int isign);
  unsigned long 
    i, i1, i2, i3, i4, np3;
  FLOAT 
    c1=0.5,c2,h1r,h1i,h2r,h2i;
  FLOAT 
    wr,wi,wpr,wpi,wtemp,theta;

  theta = M_PI/(FLOAT) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data, n>>1, 1);
  } 
  else {
    c2 = 0.5;
    theta = -theta;
  }
  wtemp = sin(0.5*theta);
  wpr   = -2.0*wtemp*wtemp;
  wpi   = sin(theta);
  wr    = 1.0+wpr;
  wi    = wpi;
  np3   = n+3;
  for (i = 2; i <= (n>>2); i++) {
    i4       =  1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r      =  c1*(data[i1]+data[i3]);
    h1i      =  c1*(data[i2]-data[i4]);
    h2r      = -c2*(data[i2]+data[i4]);
    h2i      =  c2*(data[i1]-data[i3]);
    data[i1] =  h1r+wr*h2r-wi*h2i;
    data[i2] =  h1i+wr*h2i+wi*h2r;
    data[i3] =  h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr       = (wtemp=wr)*wpr-wi*wpi+wr;
    wi       =  wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r - data[2];
  }
  else {
    data[1] = c1*((h1r=data[1])+data[2]);
    data[2] = c1*(h1r-data[2]);
    four1(data, n>>1, -1);
  }
}

/*======================= end of realft() ====================================*/

/*======================= c_num() ============================================*/

complex c_num(FLOAT x, FLOAT y)
{
  complex
    ans;

  ans.x = x;
  ans.y = y;

  return ans;
}

/*======================= end of c_num() =====================================*/

/*======================= c_mult() ===========================================*/

complex c_mult(complex z1,complex z2)
{
  complex
    ans;

  ans.x = (z1.x)*(z2.x)-(z1.y)*(z2.y);
  ans.y = (z1.x)*(z2.y)+(z1.y)*(z2.x);

  return ans;
}
/*======================= end of c_mult() ====================================*/

/*======================= c_add() ============================================*/

complex c_add(complex z1,complex z2)
{
  complex
    ans;

  ans.x = (z1.x)+(z2.x);
  ans.y = (z1.y)+(z2.y);

  return ans;
}

/*======================= end of c_add() =====================================*/

/*======================= c_sub() ============================================*/

complex c_sub(complex z1,complex z2)
{
  complex
    ans;

  ans.x = (z1.x)-(z2.x);
  ans.y = (z1.y)-(z2.y);

  return ans;
}

/*======================= end of c_sub() =====================================*/

/*======================= c_exp() ============================================*/

complex c_exp(complex z)
{
  complex 
    ans;

  ans.x = exp(z.x)*cos(z.y);
  ans.y = exp(z.x)*sin(z.y);

  return ans;
}

/*======================= end of c_exp() =====================================*/

/*======================= c_abs() ============================================*/

/*
 * NOTE: For LINUX with -D_BSD_SOURCE, cabs() is defined, such that a 
 * type-mismatch error occurs if we call this function cabs().
 */
 
FLOAT c_abs(complex z)
{

  return sqrt((z.x)*(z.x)+(z.y)*(z.y));
}

/*======================= end of c_abs() =====================================*/

/*======================= c_real() ===========================================*/

FLOAT c_real(complex z) 
{
  return (z.x);
}

/*======================= end of c_real() ====================================*/

/*======================= c_imag() ===========================================*/

FLOAT c_imag(complex z)
{
  return (z.y);
}

/*======================= end of c_imag() ====================================*/

/*======================= fcmp() =============================================*/

/*
 * Derived from fcmp(), version 1.2.2, 
 * Copyright (c) 1998-2000 Theodore C. Belding
 * University of Michigan Center for the Study of Complex Systems
 * <mailto:Ted.Belding@umich.edu>
 * <http://fcmp.sourceforge.net>
 *
 * The major modification we have made is to remove the "epsilon" argument
 * and set epsilon inside the fcmp() function.
 *
 * Description:
 *   It is generally not wise to compare two floating-point values for
 *   exact equality, for example using the C == operator.  The function
 *   fcmp() implements Knuth's suggestions for safer floating-point
 *   comparison operators, from:
 *   Knuth, D. E. (1998). The Art of Computer Programming.
 *   Volume 2: Seminumerical Algorithms. 3rd ed. Addison-Wesley.
 *   Section 4.2.2, p. 233. ISBN 0-201-89684-2.
 *
 * Input parameters:
 *   x1, x2: numbers to be compared
 *
 * This routine may be used for both single and double precision.
 *
 * Returns:
 *   -1 if x1 < x2
 *    0 if x1 == x2
 *    1 if x1 > x2		
 */

int fcmp(double x1, double x2) {
  int 
    exponent;
  double
    delta,
    difference;
#if EPIC_PRECISION == DOUBLE_PRECISION
  const double
    epsilon = DBL_EPSILON;
#else
  const double
    epsilon = FLT_EPSILON;
#endif
  
  /* 
   * Get exponent(max(fabs(x1),fabs(x2))) and store it in exponent. 
   *
   * If neither x1 nor x2 is 0,
   * this is equivalent to max(exponent(x1),exponent(x2)).
   *
   * If either x1 or x2 is 0, its exponent returned by frexp would be 0,
   * which is much larger than the exponents of numbers close to 0 in
   * magnitude. But the exponent of 0 should be less than any number
   * whose magnitude is greater than 0.
   *
   * So we only want to set exponent to 0 if both x1 and x2 are 0. 
   * Hence, the following works for all x1 and x2. 
   */
  frexp(fabs(x1) > fabs(x2) ? x1 : x2,&exponent);

  /* 
   * Do the comparison.
   *
   * delta = epsilon*pow(2,exponent)
   *
   * Form a neighborhood around x2 of size delta in either direction.
   * If x1 is within this delta neighborhood of x2, x1 == x2.
   * Otherwise x1 > x2 or x1 < x2, depending on which side of
   * the neighborhood x1 is on.
   */
  delta      = ldexp(epsilon,exponent); 
  difference = x1-x2;

  if (difference > delta) {
    /* x1 > x2 */
    return 1;
  }
  else if (difference < -delta) {
    /* x1 < x2 */
    return -1;
  }
  else  {
    /* -delta <= difference <= delta */
    return 0;  /* x1 == x2 */
  }
}

/*======================= end of fcmp() ======================================*/

/*======================= least_squares() ====================================*/

/*
 * Code for fitting data to a straight line using the least squares technique.
 *   Csaba J. Palotai *A*  4/14/2004
 * Adapted from "Numerical recipes in C," p665.
 *
 * Given a set of data points x[0..n-1],y[0..n-1] the code fits them to a 
 * straight line y= a+bx by minimizing X^2. Returned are a,b.
 *
 * Commented out is code that can be used to calculate the uncertainties 
 * siga and sigb, and the chi-square, chi2.
 */

void least_squares(FLOAT *x,
                   FLOAT *y,
                   int    n,
                   FLOAT *a)

{
  
  register int
    i;
  register FLOAT 
    siga,sigb,chi2,q,t,
    sxoss,ss,sigdat,
    sx  = 0.0,
    sy  = 0.0,
    st2 = 0.0;

  a[1] = 0.0;

  for (i = 0; i < n; i++) {
    sx += x[i];
    sy += y[i];
  }
  
  ss    = n;
  sxoss = sx/ss;
  
  for (i = 0; i < n; i++) {
    t     = x[i]-sxoss;
    st2  += t*t;
    a[1] += t*y[i];
  }
  
  a[1] /= st2;
  a[0]  = (sy-sx*a[1])/ss;

  /*  
   * Code for chi2 and uncertainties for a and b.
   *
  siga = sqrt((1.+sx*sx/(ss*st2))/ss);
  sigb = sqrt(1./st2);
  chi2 = 0.0;
  q    = 1.0;
  for (i = 0; i < n; i++) {
    t      = (y[i]-(*a)-(*b)*x[i]);
    chi2  += t*t;
    sigdat = sqrt((chi2)/(n-2));
    siga  *= sigdat;
    sigb  *= sigdat;
  }
   *
   *
   */

  return;
}  

/*======================= end of least_squares() =============================*/

/*======================= savitzky_golay() ===================================*/

/*
 * Calculate Savitzky-Golay weighting coefficients to compute a smooth value
 * or smooth derivative of noisy data.
 *
 * Based on Numerical Recipes in C, p. 652.
 * Note the wrap-around ordering of the coefficients returned in c[].
 *
 * Assumes zero-based arrays.
 */

#undef  A
#define A(i,j) a[j+(m+1)*i]
 
void savitzky_golay(FLOAT *c,
                    int    np,
                    int    nl,
                    int    nr,
                    int    ld,
                    int    m)
{
  int
    imj,ipj,j,k,kk,mm;
  int
   *index;
  FLOAT
    d,fac,sum;
  FLOAT
    *a,
    *b;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="savitzky_golay";

  /*
   * Screen for inconsistent arguments.
   */
  if (np < nl+nr+1) {
    sprintf(Message,"np=%d < nl+nr+1=%d",np,nl+nr+1);
    util_error(dbmsname,Message);
  }
  else if (nl < 0) {
    sprintf(Message,"nl=%d < 0",nl);
    util_error(dbmsname,Message);
  }
  else if (nr < 0) {
    sprintf(Message,"nr=%d < 0",nr);
    util_error(dbmsname,Message);
  }
  else if (ld > m) {
    sprintf(Message,"ld=%d > m=%d",ld,m);
    util_error(dbmsname,Message);
  }
  else if (nl+nr < m) {
    sprintf(Message,"nl+nr=%d < m=%d",nl+nr,m);
    util_error(dbmsname,Message);
  }

  /* Allocate memory */
  index = ivector(0,m,dbmsname);
  a     = fvector(0,(m+1)*(m+1)-1,dbmsname);
  b     = fvector(0,m,dbmsname);

  for (ipj = 0; ipj <= (m << 1); ipj++) {
    sum = (ipj ? 0. : 1.);
    for (k = 1; k <= nr; k++) {
      sum += pow((double)k,(double)ipj);
    }
    for (k = 1; k <= nl; k++) {
      sum += pow((double)-k,(double)ipj);
    }
    mm = IMIN(ipj,2*m-ipj);
    for (imj = -mm; imj <= mm; imj+=2) {
      A((ipj+imj)/2,(ipj-imj)/2) = sum;
    }
  }

  lu_decompose(m+1,a,index,&d);

  for (j = 0; j < m+1; j++) {
    b[j] = 0.;
  }
  b[ld] = 1.0;

  lu_backsub(m+1,a,index,b);

  for (kk = 0; kk < np; kk++) {
    c[kk] = 0.;
  }

  for (k = -nl; k <= nr; k++) {
    sum = b[0];
    fac = 1.;
    for (mm = 0; mm < m; mm++) {
      sum += b[mm+1]*(fac *= k);
    }
    kk = ((np-k)%np);
    c[kk] = sum;
  }

  /* Free allocated memory. */
  free_ivector(index,0,m,dbmsname);
  free_fvector(a,0,(m+1)*(m+1)-1,dbmsname);
  free_fvector(b,0,m,dbmsname);

  return;
}

/*======================= end of savitzky_golay() ============================*/

/*======================= random_number() ====================================*/

/*
 * Based on ran1() in Numerical Recipes in C, p. 280.
 * 
 * Returns a uniform deviate between 0.0 and 1.0 (exclusive of the endpoint values).
 *
 * Call with idum a negative number to initialize; thereafter, do not alter
 * idum between successive deviates in a sequence.
 */

#undef  IA
#define IA   16807
#undef  IM
#define IM   2147483647
#undef  AM
#define AM   (1./IM)
#undef  IQ
#define IQ   127773
#undef  IR
#define IR   2836
#undef  NTAB
#define NTAB 32
#undef  NDIV
#define NDIV (1+(IM-1)/NTAB)
#undef  EPS
#define EPS   1.2e-7
#undef  RNMX
#define RNMX (1.-EPS)

FLOAT random_number(long *idum)
{
  int
    j;
  long
    k;
  static long
    iy=0,
    iv[NTAB];
  FLOAT
    temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) {
      *idum = 1;
    }
    else {
      *idum = -(*idum);
    }
    for (j = NTAB+7; j >= 0; j--) {
      k     = (*idum)/IQ;
      *idum = IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j     = iy/NDIV;
  iy    = iv[j];
  iv[j] = *idum;
  if ((temp = AM*iy) > RNMX) {
    return RNMX;
  }
  else {
    return temp;
  }
}

/*======================= end of random_number() =============================*/

#undef DEG
#define DEG (M_PI/180.)

/*======================= lat_centric_to_graphic() ===========================*/

/*
 * lat: [Deg]
 * rerp: (equatorial radius)/(polar radius)
 */
FLOAT lat_centric_to_graphic(FLOAT lat,
                             FLOAT rerp)
{
  return (fabs(lat) == 90.) ? lat : atan(rerp*rerp*tan(lat*DEG))/DEG;
}

/*======================= end of lat_centric_to_graphic() ====================*/

/*======================= lat_graphic_to_centric() ===========================*/

/*
 * lat: [Deg]
 * rerp: (equatorial radius)/(polar radius)
 */
FLOAT lat_graphic_to_centric(FLOAT lat,
                             FLOAT rerp)
{
  return (fabs(lat) == 90.) ? lat : atan(tan(lat*DEG)/(rerp*rerp))/DEG;
}

/*======================= end of lat_graphic_to_centric() ====================*/

/*======================= surface_area_oblate() ==============================*/

/*
 * Surface area of an oblate spheroid, given equatorial and polar radii a and c,
 * respectively.
 * See http://mathworld.wolfram.com/OblateSpheroid.html.
 */

FLOAT surface_area_oblate(FLOAT a,
                          FLOAT c)
{
  register FLOAT
    e;

  e = sqrt(1.-(c/a)*(c/a));

  if (fcmp(e,0.) == 0) {
    return 4.*M_PI*a*a;
  }
  else if (fcmp(e,1.) == 0) {
    return 2.*M_PI*a*a;
  }
  else {
    return M_PI*(2.*a*a+c*c*(log((1.+e)/(1.-e))/e));
  }
}

/*======================= end of surface_area_oblate() =======================*/

/*======================= util_error() =======================================*/
/*
 * Prints calling function name, and Message to stderr, then aborts.
 */
void util_error(char *calling_function,
                char *Message)
{
  fprintf(stderr,"\n** Error: %s(): %s\n",calling_function,Message);
  fflush(stderr);

#if defined(EPIC_MPI)
  MPI_Abort(MPI_COMM_WORLD,1);
#endif

  exit(1);
}

/*======================= end of util_error() ================================*/

/* * * * * * * * * * *  end of epic_funcs_util.c  * * * * * * * * * * * * * * */

