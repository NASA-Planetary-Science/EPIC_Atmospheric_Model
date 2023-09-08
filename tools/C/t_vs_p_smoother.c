/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2019 Timothy E. Dowling                           *
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

/*
 * t_vs_p_smoother.c
 *
 * This tool smooths an EPIC t_vs_p ascii file by applying
 * diffusion via crank_nicolson(). It is useful for getting rid
 * of small superadiabatic regions that are plaguing an otherwise 
 * useful profile.
 *
 * Example compilation:
 *   %gcc -lm -O -g -o t_vs_p_smoother t_vs_p_smoother.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define TRUE    1
#define FALSE   0

#undef  N_STR
#define N_STR 256

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

/*
 * Function prototypes.
 */
void input_string(char   *prompt, 
                  char   *def, 
                  char   *ans);

int *ivector(int   nl,
             int   nh,
             char *calling_func);

void free_ivector(int  *m,
                  int   nl,
                  int   nh,
                  char *calling_func);

double *dvector(int   nl,
                int   nh,
                char *calling_func);

void free_dvector(double *m,
                  int     nl,
                  int     nh,
                  char   *calling_func);

/*
 * Choices for pivot_type in tridiag():
 */
#define WITHOUT_PIVOTING 0
#define WITH_PIVOTING    1

/*
 * The following macro is useful to insert into code while trying to corral a problem.
 */
#undef  DEBUG_MILESTONE
#define DEBUG_MILESTONE(comment) fprintf(stderr,"%s, %2d: "#comment"\n", \
                                         dbmsname,++idbms);fflush(stderr);

void tridiag(int     n,
             double *a,
             double *b,
             double *c,
             double *r,
             double *u,
             int     pivot_type);

void band_decomp(int      n,
                 int      m1,
                 int      m2,
                 double  *a,
                 double  *al,
                 int     *index,
                 double  *d);

void band_back_sub(int      n,
                   int      m1,
                   int      m2,
                   double  *a,
                   double  *al,
                   int     *index,
                   double  *b);

void band_improve(int      n,
                  int      m1,
                  int      m2,
                  double  *aorig,
                  double  *a,
                  double  *al,
                  int     *index,
                  double  *b,
                  double  *x);

void crank_nicolson(int     n,
                    double  dt,
                    double *z,
                    double *A,
                    double *mu,
                    double *rho,
                    double *ANS);

void util_error(char *calling_function,
                char *Message);

/*
 * Global declarations.
 */
char
  Message[256];


/*======================= main() ==============================================*/

int main(int   argc,
         char *argv[])
{
  int
    ntp,nn;
  char
    input_file[N_STR],
    output_file[N_STR],
    header[N_STR];
  double
   *pdat,
   *tdat,
   *dtdat,
   *zeedat,
   *sm_tdat,
   *mudat;
  FILE
   *infile,
   *outfile;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="main";

  /*
   * Prompt for input file, and open it.
   */
  sprintf(header,"");
  input_string("Input t_vs_p file:",header,input_file);
  infile = fopen(input_file,"r");
  if (!infile) {
    sprintf(Message,"Unable to open input file %s",input_file);
    util_error(dbmsname,Message);
  }
  
  sprintf(output_file,"smoothed_%s",input_file);
  outfile = fopen(output_file,"w");

  /* Read and write 6-line header */
  for (nn = 0; nn < 6; nn++) {
    fgets(header,N_STR,infile);
    if (nn == 4) {
      if (strcmp(header," .\n") == 0) {
        fprintf(outfile," Smoothed via t_vs_p_smoother (epic/tools/C).\n");
      }
      else {
        fprintf(outfile," Smoothed via t_vs_p_smoother (epic/tools/C). %s",header);
      }
    }
    else {
      fprintf(outfile,"%s",header);
    }
  }
  /* Read and write number of data points */
  fscanf(infile,"%d",&ntp);
  fprintf(outfile,"%d\n",ntp);

  /* Allocate memory */
  pdat    = dvector(0,ntp-1,dbmsname);
  tdat    = dvector(0,ntp-1,dbmsname);
  dtdat   = dvector(0,ntp-1,dbmsname);
  zeedat  = dvector(0,ntp-1,dbmsname);
  sm_tdat = dvector(0,ntp-1,dbmsname);
  mudat   = dvector(0,ntp-1,dbmsname);

  /* 
   * Store in order of increasing altitude.
   */
  for (nn = ntp-1; nn >= 0; nn--) { 
    fscanf(infile,"%lf %lf %lf",pdat+nn,tdat+nn,dtdat+nn);

    zeedat[nn] = -log(pdat[nn]);
    /*
     * NOTE: Smoothing effect is proportional to the size of mudat.
     */
    mudat[nn] = 0.01;
  }

  fclose(infile);

  /*
   * Smooth t_vs_p profile
   */
  sm_tdat[0]     = tdat[0];
  sm_tdat[ntp-1] = tdat[ntp-1];
  crank_nicolson(ntp-2,1.,zeedat,tdat,mudat,NULL,sm_tdat);

  /*
   * Output smoothed t_vs_p profile
   */

  for (nn = ntp-1; nn >= 0; nn--) { 
    fprintf(outfile," %10.4e    %7.2f     %.2f\n",pdat[nn],sm_tdat[nn],dtdat[nn]);
  }

  fclose(outfile);

  /*
   * Free allocated memory.
   */
  free_dvector(pdat,   0,ntp-1,dbmsname);
  free_dvector(tdat,   0,ntp-1,dbmsname);
  free_dvector(dtdat,  0,ntp-1,dbmsname);
  free_dvector(zeedat, 0,ntp-1,dbmsname);
  free_dvector(sm_tdat,0,ntp-1,dbmsname);
  free_dvector(mudat,  0,ntp-1,dbmsname);

  return 0;
}

/*======================= end of main() =======================================*/

/*----------------------------------------------------------------------------*
 * The following i/o functions are adapted from epic/src/core/epic_funcs_io.c *
 *----------------------------------------------------------------------------*/

/*====================== input_string() ======================================*/

/* 
 * Read in a string, or set to default if input is a return ('\n').
 *
 * C.Santori, T.Dowling 
 */

void input_string(char *prompt, 
                  char *def, 
                  char *ans) 
{
  char 
    c,
    buffer[N_STR];
  int  
    len;

  fprintf(stdout,"%s[%s]: ",prompt,def);
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

/*-------------------------------------------------------------------------------*
 * The following math functions are adapted from epic/src/core/epic_funcs_util.c *
 *-------------------------------------------------------------------------------*/

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

void crank_nicolson(int     n,
                    double  dt,
                    double *z,
                    double *A,
                    double *mu,
                    double *rho,
                    double *ANS)
{
  int
    k,kstart,kend;
  double
    coeff1,coeff3,factor;
  static int
    nold = -1;
  static double
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
    a    = dvector(0,n-1,dbmsname);
    b    = dvector(0,n-1,dbmsname);
    c    = dvector(0,n-1,dbmsname);
    r    = dvector(0,n-1,dbmsname);
    nold = n;
  }
  else if (nold < n) {
    free_dvector(a,0,nold-1,dbmsname);
    free_dvector(b,0,nold-1,dbmsname);
    free_dvector(c,0,nold-1,dbmsname);
    free_dvector(r,0,nold-1,dbmsname);
    a    = dvector(0,n-1,dbmsname);
    b    = dvector(0,n-1,dbmsname);
    c    = dvector(0,n-1,dbmsname);
    r    = dvector(0,n-1,dbmsname);
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

void tridiag(int     n,
             double *a,
             double *b,
             double *c,
             double *r,
             double *u,
             int     pivot_type)
{
  int
    j,
    m1,m2,mm,
    *index;
  double
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
    gam = dvector(0,n-1,dbmsname);

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
        free_dvector(gam,0,n-1,dbmsname);
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
    free_dvector(gam,0,n-1,dbmsname);
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
    aa     = dvector(0,n*mm-1,dbmsname);
    aaorig = dvector(0,n*mm-1,dbmsname);
    aal    = dvector(0,n*m1-1,dbmsname);
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
    free_dvector(aa,    0,n*mm-1,dbmsname);
    free_dvector(aaorig,0,n*mm-1,dbmsname);
    free_dvector(aal,   0,n*m1-1,       dbmsname);
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

void band_decomp(int      n,
                 int      m1,
                 int      m2,
                 double  *a,
                 double  *al,
                 int     *index,
                 double  *d)
{
  int
    i,j,k,l,
    mm;
  double
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

void band_back_sub(int      n,
                   int      m1,
                   int      m2,
                   double  *a,
                   double  *al,
                   int     *index,
                   double  *b)
{
  int
    i,k,l,
    mm;
  double
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
                  double *aorig,
                  double *a,
                  double *al,
                  int    *index,
                  double *b,
                  double *x)
{
  int
    k,j,i,tmploop;
  static int
    n_max = 0;
  static double
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
    r     = dvector(0,n_max-1,dbmsname);
  }
  else if (n > n_max) {
    free_dvector(r,0,n_max-1,dbmsname);
    n_max = n;
    r     = dvector(0,n_max-1,dbmsname);
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

/*======================= ivector() ============================================*/
      
/*
 *  Allocates memory for a 1D int array 
 *  with range [nl..nh].
 */

#undef  DEBUG

int *ivector(int   nl, 
             int   nh,
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
                  int   nl, 
                  int   nh,
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

/*======================= dvector() ============================================*/
      
/*
 *  Allocates memory for a 1D double array 
 *  with range [nl..nh].
 */

#undef  DEBUG 

double *dvector(int   nl, 
                int   nh,
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

/*======================= util_error() =======================================*/
/*
 * Prints calling function name, and Message to stderr, then aborts.
 */
void util_error(char *calling_function,
                char *Message)
{
  fprintf(stderr,"\n** Error: %s(): %s\n",calling_function,Message);
  fflush(stderr);

  exit(1);
}

/*======================= end of util_error() ================================*/




