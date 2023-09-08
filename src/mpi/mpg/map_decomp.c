/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998 Joseph Matarese                              *
 *                                                                 *
 * Modified as indicated below by Tim Dowling                      *
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "mpg.h"
#include "qa.h"

#define TRUE 1

/* Function prototypes: */
void balance (int  numproc, 
              int  dim, 
              int *dims, 
              int *length);

/* Replaces factor(); T. Dowling, 13 Feb 2019 */
int prime_factors(int  num,
                  int *factors);

int  comp (const void *i, 
           const void *j);

/*=================== MPG_Cart_decomp() ===========================*/

void  MPG_Cart_decomp(MPI_Comm  comm,
                      int       dim,
                      int      *length,
                      int      *dims, 
		      int      *periods,
                      int      *pad,
                      int      *start,
                      int      *end,
		      int      *stride,
                      MPI_Comm *comm_cart) 
{
  int  
    i,numproc,
    *coords;

  /* Allocate memory for coords: */
  if ((coords = (int *)malloc(dim*sizeof(int))) == NULL) {
    fprintf (stderr,"Cannot allocate space\n");
    MPI_Abort(MPI_COMM_WORLD,0);
  }

  MPI_Comm_size(comm,&numproc);

  balance(numproc,dim,dims,length);

  for (i = 0; i < dim; i++) {
    debug ('v',"i = %d; ",i);
    debug ('v',"dims[i] = %d\n",dims[i]);
  }

  MPI_Dims_create(numproc,dim,dims);
  MPI_Cart_create(comm,dim,dims,periods,TRUE,comm_cart);
  MPI_Cart_get(*comm_cart,dim,dims,periods,coords);

  for (i = 0; i < dim; i++) {
    MPE_Decomp1d(length[i],dims[i],coords[i],&start[i],&end[i]);

    stride[i] = end[i]-start[i]+1+2*pad[i];
    start[i]--;  
    end[  i]--;

    debug ('v',"i = %d; ",               i);
    debug ('v',"dims[i] = %d; ",    dims[i]);
    debug ('v',"start[i] = %d; ",  start[i]);
    debug ('v',"end[i] = %d; ",      end[i]);
    debug ('v',"stride[i] = %d\n",stride[i]);
  }

  (void)free ((char *)coords);

  return;
}

/*=================== end of MPG_Cart_decomp() ====================*/

/*=================== balance() ===================================*/

void balance (int  numproc, 
              int  dim, 
              int *dims, 
              int *length) 
{
  int  
    i,j,nf,
    max_dim,
    max_length,
    numleft,
   *dims_copy,
   *length_copy,
   *factors;

  /* Allocate memory for dims_copy,length_copy: */
  if ((dims_copy   = (int *)malloc(dim*sizeof(int))) == NULL ||
      (length_copy = (int *)malloc(dim*sizeof(int))) == NULL) {
    fprintf (stderr, "mpg:map_decomp:cannot allocate memory\n");
    MPI_Abort (MPI_COMM_WORLD, 0);
  }

/*---------------------------------------------------------------------------*
 * if dims[i] == 0 for any i, just return
 * if dims[i] <  0, "fill it in" by rotating bisections (greatest dimlen)
 *---------------------------------------------------------------------------*/

  for (i = 0; i < dim; i++) {
    if (dims[i] == 0) {
      /* Free allocated memory: */
      (void)free((char *)dims_copy);
      (void)free((char *)length_copy);
      return;
    }
  }

  numleft = numproc;
  for (i = 0; i < dim; i++) {
    dims_copy[  i] = dims[  i];
    length_copy[i] = length[i];
    if (dims[i] < 0) {
      dims[i] = 1;
    }
    numleft /= dims[i];
    if (numleft < 1) {
      fprintf(stderr,"mpg:map_decomp:cannot subdivide grid further\n");
      MPI_Abort(MPI_COMM_WORLD,0);
    }
  }

  /* Determine number of prime factors in numleft */
  nf = prime_factors(numleft,NULL);

  /* Allocate memory for factors[] */
  factors = (int *)malloc(nf*sizeof(int));

  if (!factors) {
    fprintf(stderr,"Cannot allocate space for factors[]\n");
    MPI_Abort(MPI_COMM_WORLD,0);
  }

  /* Store factors */
  prime_factors(numleft,factors);

  for (j = 0; j < nf; j++) {
    max_dim    = -1;  
    max_length =  0;
    for (i = 0; i < dim; i++) {
      if (dims_copy[i] < 0 && max_length <= length_copy[i]) {
        max_dim    = i;
        max_length = length_copy[i];
      }
    }
    if (max_dim != -1) {
      length_copy[max_dim] /= factors[j]; 
      if (length_copy[max_dim] < 1) {
        fprintf (stderr,"mpg:map_decomp:cannot subdivide grid further\n");
        MPI_Abort(MPI_COMM_WORLD,0);
      }
      dims[max_dim] *= factors[j];
    }
  }

  /* Free allocated memory: */
  (void)free((char *)dims_copy);
  (void)free((char *)length_copy);
  (void)free((char *)factors);

  return;
}

/*=================== end of balance() ============================*/

/*=================== prime_factors() =============================*/

/*
 * Adapted from an algorithm by Vishwas Garg.
 *   https://www.geeksforgeeks.org/print-all-prime-factors-of-a-given-number/
 *
 * Call initially as
 *   nf = prime_factors(num,NULL)
 * to count the number of factors, nf.
 * Allocate memory in the calling function for factors[], and call as
 *   prime_factors(num,factors)
 * to store the factors.
 *
 * NOTE: This replaces the original program factor(), which did not return
 *       correct results for many input values, including 5, 7, 10, 14, and 20.
 *
 * T. Dowling, 13 Feb 2019.
 */

int prime_factors(int  num,
                  int *factors)
{
  int
    i,n,count;

  /*
   * Special case of num = 1.
   */
  if (num == 1) {
    if (factors) {
      factors[0] = 1;
    }
    
    return 1;
  }

  n     = num;
  count = 0;
  while (n%2 == 0) {
    if (factors) {
      factors[count] = 2;
    }
    count++;
    n /= 2;
  }

  for (i = 3; i <= (int)sqrt((double)n); i+=2) {
    while (n%i == 0) {
      if (factors) {
        factors[count] = i;
      }
      count++;
      n /= i;
    }
  }

  if (n > 2) {
    if (factors) {
      factors[count] = n;
    }
    count++;
  }

  if (factors) {
    qsort(factors,count,sizeof(int),comp);
  }

  return count;
}

/*=================== end of prime_factors() ======================*/

/*=================== comp() ======================================*/

int comp(const void *i,
         const void *j) 
{
  return (*(int *)j - *(int *)i);
}

/*=================== end of comp () ==============================*/

