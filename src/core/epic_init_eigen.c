/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2023 Timothy E. Dowling                      *
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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Charles Santori, Feb 1997; originally named low_mode().         *
 * Terrestrial boundary conditions added Jan 1999, T. Dowling.     *
 * Estimate of eigenvectors added and name changed to              *
 * vertical_modes(), Nov 2006, T. Dowling.                         *
 *                                                                 *
 * This function finds the eigenvalues, c, and eigenvectors, psi,  *
 * of the equation:                                                *
 *                                                                 *
 *     z  d   / -z    1    d(psi) \         1                      *
 *    e * -- ( e  * -----  ------  )  = - ----- (psi) ,            *
 *        dz  \     (NH)^2   dz   /        c^2                     *
 *                                                                 *
 *    where z = -ln(p/p0), and writes the results to a file.       *
 *                                                                 *
 * See Achterberg and Ingersoll (1989, J.Atmos.Sci 46, 2448-2462). *
 *                                                                 *
 *              Gas-giant planets   Terrestrial planets            *
 *    Top b.c.:   d(psi)/dz = 0        d(psi)/dz = 0               *
 * Bottom b.c.:     psi     = 0        d(psi)/dz = 0               *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*===================== vertical_modes() ==========================*/

#undef  TOL
#define TOL      1.e-4
#undef  MAX_VECS
#define MAX_VECS 10

void vertical_modes(planetspec *planet,
                    int         J,
                    int         I)
{
  int 
    K,kk,klen,k,
    kstart,
    row,n,iter,
    nvecs,
    istart,
    max_starts = 5,
    max_iters  = 10*grid.nk,
    converged[MAX_VECS];
  int
    *index;
  static long
    seed=-42355;
  EPIC_FLOAT 
    press,temperature,
    rgas,brunt2,hscale,nh2,
    clow,mult,cn,
    bb,change,max_y,d,
    c_ans[MAX_VECS],
    psi_ans[MAX_VECS][KHI];
  EPIC_FLOAT
    *z,*dz,*b,
    *wr,*wi,*a_eig,*aa_eig,
    *bee,*oldbee;
  FILE
    *outfile;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="vertical_modes";

  rgas = planet->rgas;
  klen = KHI;

  /* 
   * Allocate memory: 
   */
  z      = fvector(1,2*KHI+1,dbmsname);
  dz     = fvector(2,2*KHI+1,dbmsname);
  b      = fvector(2,2*KHI+1,dbmsname);
  wr     = fvector(0,klen-1, dbmsname);
  wi     = fvector(0,klen-1, dbmsname);
  bee    = fvector(0,klen-1, dbmsname);
  oldbee = fvector(0,klen-1, dbmsname);
  index  = ivector(0,klen-1, dbmsname);
  a_eig  = fvector(0,klen*klen-1,dbmsname);
  aa_eig = fvector(0,klen*klen-1,dbmsname);

  /* Calculate z = -log(p/p0): */
  for (kk = 1; kk <= 2*KHI+1; kk++) {
    press = get_p(planet,P2_INDEX,kk,J,I);
    z[kk] = -log(press/planet->p0);
  }

  /* Calculate dz, b = e^-z/((NH)^2*dz): */
  for (kk = 2; kk <= 2*KHI+1; kk++) {
    K = kk/2;
    if (kk <= 2*KHI) {
     dz[kk] = z[kk-1]-z[kk+1];
    }
    else {
      dz[kk] = dz[kk-1];
    }
    if (kk%2 == 0) {
      /* layer value */
      temperature = T2(K,J,I);
    }
    else {
      /* interface value */
      temperature = T3(K,J,I);
    }

    brunt2 = get_brunt2(planet,kk,J,I);
    /*
     * Limit N^2 to be positive-definite for this calculation.
     */
    brunt2 = MAX(brunt2,1.e-16);

    hscale = rgas*temperature/grid.g[kk][2*J+1];
    nh2    = brunt2*hscale*hscale;
    b[kk]  = exp(-z[kk])/(nh2*dz[kk]);
  }

  /*
   * The eigenvector psi is aligned inside the layers. 
   * Apply the top boundary condition, d(psi)/dz = 0.
   */
  row = 0;
  kk  = 2*row+2;
  mult             =  exp(z[kk])/dz[kk];
  A_EIG(row,row  ) = -mult*b[kk+1];
  A_EIG(row,row+1) =  mult*b[kk+1];

  /* 
   * Apply bottom boundary condition: 
   */
  row = klen-1;
  kk  = 2*row+2;
  mult = exp(z[kk])/dz[kk];
  if (strcmp(planet->type,"gas-giant") == 0 || grid.coord_type == COORD_ISENTROPIC) {
    /* 
     * Assume psi = 0. in planet interior. 
     */
    A_EIG(row,row-1) =  mult*(b[kk-1]        );
    A_EIG(row,row  ) = -mult*(b[kk-1]+b[kk+1]);
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    /* 
     * d(psi)/dz = 0. at surface.
     */
    A_EIG(row,row-1) =  mult*b[kk-1];
    A_EIG(row,row  ) = -mult*b[kk-1];
  }
  else {
    sprintf(Message,"Unknown planet->type = %s",planet->type);
    epic_error(dbmsname,Message);
  }

  /*
   * Calculate interior points:
   */
  for (row = 1; row <= klen-2; row++) {
    kk = 2*row+2;
    mult             =  exp(z[kk])/dz[kk];
    A_EIG(row,row-1) =  mult*(b[kk-1]        );
    A_EIG(row,row  ) = -mult*(b[kk-1]+b[kk+1]);
    A_EIG(row,row+1) =  mult*(        b[kk+1]);
  }

  /* 
   * Find eigenvalues: 
   */
  hqr(a_eig,klen,wr,wi);

  /* 
   * Sort in increasing order: 
   */
  quicksort(wr,0,klen-1);

  /*
   * Loop over modes, store real eigenvalues, c.
   */
  nvecs = 0;
  /*
   * Skip to first baroclinic mode.
   */
  kstart = 1;
  for (k = kstart; k < klen; k++) {
    /* Want only real c. */
    if (wr[k] < 0. && fabs(wi[k]/wr[k]) <= .001) {
      c_ans[nvecs++] = 1./sqrt(-wr[k]);
    }      
    if (nvecs == MAX_VECS || nvecs >= klen/2) break;
  }

  /*
   * Loop over first nvecs real eigenvalues.
   * Use inverse iteration to calculate the corresponding eigenvectors.
   *
   * NOTE: We need a better way to find the eigenvectors, because the
   *       inverse iteration method converges for only about half of the
   *       eigenvalues determined above.  This is outside the scope of
   *       Numerical Recipes in C because our matrix, while only tridiagonal,
   *       is not symmetric.
   */
  memcpy(aa_eig,a_eig,klen*klen*sizeof(EPIC_FLOAT));
  for (n = 0; n < nvecs; n++) {
    converged[n] = FALSE;

    /*
     * Try a handful (max_starts) of different random starting vectors
     * to increase the chance of getting a converged solution.
     */
    for (k = 0; k < klen; k++) {
      AA_EIG(k,k) =  A_EIG(k,k)-(-1./(c_ans[n]*c_ans[n]));
    }
    lu_decompose(klen,aa_eig,index,&d);

    for (istart = 0; istart < max_starts; istart++) {
      /*
       * Randomly initialize the righthand-side vector, bee.
       */
      for (k = 0; k < klen; k++) {
        bee[k]    =  random_number(&seed);
        oldbee[k] =  1.;
      }

      for (iter = 0; iter < max_iters; iter++) {
        /*
         * Normalize the righthand-side vector.
         */
        bb = 0.;
        for (k = 0; k < klen; k++) {
          bb += bee[k]*bee[k];
        }
        bb = sqrt(bb);
        for (k = 0; k < klen; k++) {
          bee[k] /= bb;
        }

        /*
         * Check for convergence.
         */
        change = 0.;
        for (k = 0; k < klen; k++) {
          change    += SQR(bee[k]-oldbee[k]);
          oldbee[k]  = bee[k];
        }
        change = sqrt(change);
        if (change <= TOL) {
          converged[n] = TRUE;
          break;
        }

        if (change == 2.) {
          /* Break out of a non-converging limit cycle. */
          break;
        }
        
        /*
         * Solve for y, assign y to bee, and loop to the next iteration.
         */
        lu_backsub(klen,aa_eig,index,bee);
      }
      if (converged[n]) {
        break;
      }
    }

    for (k = 0; k < klen; k++) {
      psi_ans[n][k] = bee[k];
    }
  }

  /* 
   * Write results to vertical_modes.dat file.
   */
  outfile = fopen("vertical_modes.dat","w");
  fprintf(outfile,"  EPIC model output from vertical_modes(). \n");
  fprintf(outfile,"  planet->name = %s.\n",planet->name);
  fprintf(outfile,"  Eigenvectors are normalized and headed with their eigenvalues, c [m/s].\n");
  fprintf(outfile,"  The first %d baroclinic modes are listed.\n",nvecs);
  fprintf(outfile,"  A column of zeros implies that the eigenvector solution did not converge.\n\n");
  fprintf(outfile,"   p[hPa]  ");
  for (n = 0; n < nvecs; n++) {
    fprintf(outfile,"  %7.2f",c_ans[n]);
  }
  fprintf(outfile,"\n");
  for (k = 0; k < klen; k++) {
    fprintf(outfile,"  %9.3e",P2(k+1,J,I)*.01);
    for (n = 0; n < nvecs; n++) {
      if (converged[n]) {
        fprintf(outfile,"  %7.4f",psi_ans[n][k]);
      }
      else {
        fprintf(outfile,"  %7.4f",0.);
      }
    }
    fprintf(outfile,"\n");
  }

  fclose(outfile);

  /* 
   * Free allocated memory: 
   */
  free_fvector(aa_eig,0,klen*klen-1,dbmsname);
  free_fvector(a_eig, 0,klen*klen-1,dbmsname);
  free_ivector(index, 0,klen-1,     dbmsname);
  free_fvector(oldbee,0,klen-1,     dbmsname);
  free_fvector(bee,   0,klen-1,     dbmsname);
  free_fvector(wi,    0,klen-1,     dbmsname);
  free_fvector(wr,    0,klen-1,     dbmsname);
  free_fvector(b,     2,2*KHI+1,    dbmsname);
  free_fvector(dz,    2,2*KHI+1,    dbmsname);
  free_fvector(z,     1,2*KHI+1,    dbmsname);
}

/*===================== end of vertical_modes() =====================*/

/* * * * * * * * * * * * end of epic_init_eigen.c  * * * * * * * * * */
