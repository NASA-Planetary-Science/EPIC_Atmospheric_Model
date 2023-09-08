/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2008 Kunio M. Sayanagi                            *
 *                    Timothy E. Dowling                           *
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
 * Software Foundation, Inc., 59 Temple Place - Suite 330,         *
 * Boston, MA  02111-1307, USA.                                    *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * *  savitzky-golay_filter.c  * * * * * * * * * * * * 
 *                                                                 *
 *       Kunio M. Sayanagi                                         *
 *       Timothy E. Dowling                                        *
 *								   *
 * This program filters applies the Savitzky-Golay filter to any   *
 * source data with the following format:                          *
 *   1. the first line is the total data points in file (integer)  *
 *   2. the subsequent lines have independent and dependent vars   *
 *      separated by spaces, one ind-dep pair per line             *
 *                                                                 *
 * This filter can be used for both smoothing and numerical        *
 * differentiation -- see the Numerical Recipes for details.       *
 *                                                                 *
 * Note that this implementation dynamically adjusts the size of   *
 * "moving window" near the edges of the data.                     *
 *								   *
 * To complile:                                                    *
 *  gcc -lm -o savitzky-golay_filter.x savitzky-golay_filter.c     *
 * ****************************************************************/


#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define FLOAT double

static int iminarg1,iminarg2;

#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define N_STR 256		      

void sg_driver(FLOAT *data_in,    // input:  source data
 	       FLOAT *data_out,   // output: smoothed data
	       int    ndat,       // input:  total size of the input/output data
	       int    nl_def,     // input:  size of "moving window" to the left
	       int    nr_def,     // input:  size of "moving window" to the right
	       int    ld,         // input:  order of derivative  
	       int    moment);    // input:  largest moment of the data to conserve 

void savitzky_golay(FLOAT *c,
                    int    np,
                    int    nl,
                    int    nr,
                    int    ld,
                    int    m);

void lu_decompose(int   n,
                 FLOAT *a,
                 int   *index,
                 FLOAT *d);

void lu_backsub(int    n,
                FLOAT *a,
                int   *index,
                FLOAT *b);

int *ivector(int  nl, 
             int  nh );

void free_ivector(int  *m, 
                  int  nl, 
                  int  nh);

FLOAT *fvector(int  nl, 
               int  nh);


void free_fvector(FLOAT *m, 
                  int     nl, 
                  int     nh);

FLOAT input_float(char  prompt[N_STR], 
		  FLOAT def);

int input_int(char prompt[N_STR], 
	      int  def);

void input_string(char prompt[N_STR], 
		  char def[N_STR], 
		  char *ans); 


/*======================= main() =================================================*/
int main()
{
  // Things needed for reading source data...
  int 
    ndat;
  FLOAT 
    *Xdata, 
    *Ydata_source,
    *Ydata_smooth;
  FILE
    *data_source,
    *data_smooth;
  int
    ikount, jkount;
  char
    filename[N_STR];

  // Parameters for Savitzky-Golay
  int   
    nl_def,     // size of "moving window" to the left
    nr_def,     // size of "moving window" to the right
    ld,         // order of derivative -- 0 = simple smoothing, 1 = smoothed first derivative, etc.
    moment;     // largest moment of the data to conserve 

  // set the input parameters for SG filter...
  nl_def = 6;
  nr_def = 6;
  ld = 0;
  moment = 4;

  nl_def = input_int("Enter window size to the left,  nl", nl_def);
  nr_def = input_int("Enter window size to the right, nr", nr_def);
  ld     = input_int("Enter the order of derivative to output (0 = simple smoothing, max=5)", ld);
  moment = input_int("Enter the largest moment of data to preserve", moment);

  // Read in the u_vs_lat data to be filtered...
  strcpy(filename, "source.dat");
  input_string("Enter the source file name", filename, filename);
  data_source = fopen(filename,"r");
  fscanf(data_source, "%d", &ndat);

  Xdata         = fvector(0, ndat-1);
  Ydata_source  = fvector(0, ndat-1);
  Ydata_smooth  = fvector(0, ndat-1);  

  for (ikount = 0; ikount < ndat; ikount++) {
    fscanf(data_source, "%lf %lf", (Xdata+ikount), (Ydata_source+ikount));
  }
  fclose(data_source);    

  // Call sg_driver to smooth the data, to be stored in Ydata_smooth...
  sg_driver( Ydata_source,
	     Ydata_smooth,
	     ndat, nl_def, nr_def, ld, moment);

  // Write the data to file...
  strcpy(filename, "smoothed.dat");
  input_string("Enter the output file name", filename, filename);
  
  data_smooth = fopen(filename,"w");
  fprintf(data_smooth, "%d \n", ndat);
  for (ikount = 0; ikount < ndat; ikount++)
  {
    fprintf(data_smooth, "%6.5f  %6.5f \n", Xdata[ikount], Ydata_smooth[ikount]);
  }

  fclose(data_smooth);

  // Free memory and wrap things up...
  // printf("Finished calculating stuff ... free memory and exit \n");
  free_fvector(Xdata,        0, ndat-1);
  free_fvector(Ydata_source, 0, ndat-1);
  free_fvector(Ydata_smooth, 0, ndat-1);

  return 0;
}
/*======================= end of main() =================================================*/

/*======================= sg_driver() =================================================*/
#define SG(i,j) sg[j+i*(np_def)]
void sg_driver(FLOAT *data_in,    // input:  source data
	       FLOAT *data_out,   // output: smoothed data
	       int    ndat,       // input:  total size of the input/output data
	       int    nl_def,     // input:  size of "moving window" to the left
	       int    nr_def,     // input:  size of "moving window" to the right
	       int    ld,         // input:  order of derivative -- 0 = simple smoothing, 1 = smoothed first derivative, etc.
	       int    moment)     // input:  largest moment of the data to conserve 
{
// ====================================================================================
// This implements the Savitzky-Golay filter using savitzky_golay() function.
// See the Numerical Recipes for the details -- the variable names are kept the same.
//
// This function can either do a simple smoothing (with ld=0) or calculate smoothed
// derivatives with higher values of ld.
//
// Note, in this implementation, the size of the moving window are adjusted 
// near the edge of the input data.
//
// Kunio 11-05-2008
// ====================================================================================

  // Things needed for S-G filter...
  // NOTE: these are adjusted near the edges of the data array.
  int
    np_def, 
    nl, // "left" window
    nr, // "right" window
    nn, // pointer for the SGcoeffs
    np; // num coefficients
//    ld, // the order of the derivative to calculate using SG filter
//    moment;  // highest moment to conserve in the smoothing process:
        //   m = 0 -- area under the curve
	//   m = 1 -- mean position of the curve
	//   m = 2 -- width of peaks
	//   m = ...
  
  FLOAT
    checksum,      // total of the SG coeffs should equal approx 1.0
    *sg; // SG coefficients array

  int
    ikount, jkount;

  np_def = nl_def +nr_def +1;
  sg     = fvector(0,ndat*np_def-1);

  /*
   * Calculate the Savitsky-Golay coefficients.
   *    g(i) = sum(c(i, j)*f(i), j=-nl..nr)
   *    where: g(i)    -- smoothed function
   *           f(i)    -- original input function
   *           c(i, *) -- S-G coefficients to calculate the i-th element
   */
  for (ikount = 0; ikount < ndat; ikount++) 
  {
    nl = (ikount < nl_def)             ? ikount           : nl_def;
    nr = ( (ndat -ikount-1) < nr_def ) ? (ndat -ikount-1) : nr_def;
    np = nl +nr +1;
    checksum = 0.0;

    //printf("Calculating sg coeffs, ikount = %d \n", ikount);
    savitzky_golay(sg+ikount*np_def,np,nl,nr,ld,moment);
    
    // Below lines are for debugging
/*
    printf("Printing the SG coeffs -- ikount = %d, nl = %d, nr = %d, np = %d \n", ikount, nl, nr, np);
    for (jkount = 0; jkount <= nl; jkount++)
    {
      printf("%f, ", SG(ikount, jkount) );
      checksum += SG(ikount, jkount);
    }
    for (jkount = 1; jkount <= nr; jkount++)
    {
      printf("%f, ", SG(ikount, np-jkount) );
      checksum += SG(ikount, np-jkount);
    }
    
    printf("checksum = %f \n", checksum);
*/    
  }

  // Now, use the coeffs to calculate the smoothed profile  
  for (ikount = 0; ikount < ndat; ikount++)
  { 
    data_out[ikount] = 0; // zero the array element...

    nl = ( ikount < nl_def)            ? ikount           : nl_def;
    nr = ( (ndat -ikount-1) < nr_def ) ? (ndat -ikount-1) : nr_def;
    np = nl +nr+1;

    // First, sum the elements to the left...
    for (jkount = 0; jkount <= nl; jkount++)
    {
      data_out[ikount] += SG(ikount, jkount)*data_in[ikount -jkount];
    }

    // Second, sum the elements to the right...
    for (jkount = 1; jkount <= nr; jkount++)
    {
      data_out[ikount] += SG(ikount, np-jkount)*data_in[ikount +jkount];
    }
  }

  // Free memory and return...
  free_fvector(sg,           0, ndat*np_def-1);

  return;
}	

/*======================= end of sg_driver() ============================*/

/*======================= savitzky_golay() ===================================*/

/*
 * Calculate Savitzky-Golay weighting coefficients to compute a smooth value
 * or smooth derivative of noisy data.
 *
 * Based on Numerical Recipes in C, p. 652.
 * Note the wrap-around ordering of the coefficients returned in c[].
 *
 * Assumes zero-based arrays.
 *
 * Kunio 11-05-2008 -- took this from EPIC4 package -- 
 *   Took out debugging milestones in the process, but in retrospect I should have kept them...
 */

#undef  A
#define A(i,j) a[j+(m+1)*i]
 
void savitzky_golay(FLOAT *c, // 
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
   * Screen for inconsistent arguments.
   */
  if (np < nl+nr+1) {
    fprintf(stderr,"**error: np=%d < nl+nr+1=%d \n",np,nl+nr+1);
    exit(1);
  }
  else if (nl < 0) {
    fprintf(stderr,"**error: nl=%d < 0 \n",nl);
    exit(1);
  }
  else if (nr < 0) {
    fprintf(stderr,"**error: nr=%d < 0 \n",nr);
    exit(1);
  }
  else if (ld > m) {
    fprintf(stderr,"**error: ld=%d > m=%d \n",ld,m);
    exit(1);
  }
  else if (nl+nr < m) {
    fprintf(stderr,"**error: nl+nr=%d < m=%d \n",nl+nr,m);
    exit(1);
  }

  /* Allocate memory */
  index = ivector(0,m);
  a     = fvector(0,(m+1)*(m+1)-1);
  b     = fvector(0,m);

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
  free_ivector(index,0,m);
  free_fvector(a,0,(m+1)*(m+1)-1);
  free_fvector(b,0,m);

  return;
}

/*======================= end of savitzky_golay() ============================*/

/*======================= ivector() ============================================*/
      
/*
 *  Allocates memory for a 1D int array 
 *  with range [nl..nh].
 */

#undef  DEBUG

int *ivector(int  nl, 
             int  nh )
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  int         
    *m;

  if (nh < nl) {
    fprintf(stderr,"**error: ivector() range (%d,%d)\n",nl,nh);
    exit(1);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (int *)calloc(len_safe,sizeof(int));
  if (!m) {
    fprintf(stderr,"**error ivector() \n");
    exit(1);
  }
  m -= nl_safe;

#if defined(DEBUG)
  fprintf(stderr,"ivector() error \n");
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
                  int  nh)
{
  int  
    nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m += nl_safe;
  free(m);

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
               int  nh)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  FLOAT         
    *m;

  if (nh < nl) {
    fprintf(stderr,"**error:fvector, range (%d,%d)\n",nl,nh);
    exit(1);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (FLOAT *)calloc(len_safe,sizeof(FLOAT));

  if (!m) {
    fprintf(stderr,"**error: fvector(), nl=%d,nh=%d,len_safe=%d\n",nl,nh,len_safe);
    exit(1);
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
                  int     nh)
{
  int  
    nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m += nl_safe;
  free(m);

  return;
}

/*======================= end of free_fvector() ===============================*/

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
   * Allocate memory.
   */
  if (n_max == 0) {
    n_max = n;
    vv    = fvector(0,n_max-1);
  }
  else if (n > n_max) {
    free_fvector(vv,0,n_max-1);
    n_max = n;
    vv    = fvector(0,n_max-1);
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
      fprintf(stderr,"**error:lu_decompose(), matrix singular\n");
      exit(1);
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
      fprintf(stderr,"Warning, lu_decompose(), matrix is singular\n");
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

/*======================= input_float() =======================================*/

/* 
 * Read in float, or set to default if input is a return ('\n').
 * C.Santori, T.Dowling 
 */

FLOAT input_float(char prompt[N_STR], FLOAT def) 
{
  char  
    c,
    buffer[N_STR];
  int   
    len;
  FLOAT 
    ans;

  fprintf(stderr,"%s[%g]: ",prompt,def);
  for (len = 0; (c = getchar()) != '\n' && len < N_STR; len++) {
    buffer[len]=c;
  }
  buffer[len] = '\0';
  if (len == 0) {
    ans = def;
  }
  else {
    sscanf(buffer,"%f",&ans);
  }

  return ans;
}

/*====================== end input_float() =====================================*/


/*====================== input_int() ========================================*/

/* 
 * Read in int, or set to default if input is a return ('\n').
 * C.Santori, T.Dowling 
 */

int input_int(char prompt[N_STR], int def) 
{
  char  
    c,
    buffer[N_STR];
  int 
    ans,
    len;

  fprintf(stderr,"%s[%d]: ",prompt,def);
  for (len = 0; (c = getchar()) != '\n' && len < N_STR; len++) {
    buffer[len]=c;
  }
  buffer[len] = '\0';
  if (len == 0) {
    ans = def;
  }
  else {
    sscanf(buffer,"%d",&ans);
  }

  return ans;
}

/*====================== end input_int() ====================================*/

/*====================== input_string() =====================================*/

/* 
 * Read in a string, or set to default if input is a return ('\n').
 *
 * C.Santori, T.Dowling 
 */

void input_string(char prompt[N_STR], char def[N_STR], char *ans) 
{
  char 
    c,
    buffer[N_STR];
  int  
    len;

  fprintf(stderr,"%s[%s]: ",prompt,def);
  for (len = 0; (c = getchar()) != '\n' && len < N_STR; len++) {
    buffer[len]=c;
  }
  buffer[len] = '\0';
  if (len == 0) {
    strcpy(ans,def);
  }
  else {
    strcpy(ans,buffer);
    strcpy(def,buffer);
  }

}

/*====================== end input_string() ==================================*/
