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

/* * * * * * * * *  epic_funcs_thermo.c  * * * * * * * * * * * * * * * *
 *                                                                     *
 *  Timothy E. Dowling                                                 *
 *                                                                     *
 *  Thermodynamic functions that do not reference EPIC model variables *
 *  by name or index (functions may reference planetspec).             *
 *                                                                     *
 *  This file includes the following:                                  *
 *                                                                     *
 *      thermo_setup()                                                 *
 *      return_temp()                                                  *
 *      alt_return_temp()                                              *
 *      return_density()                                               *
 *      return_theta()                                                 *
 *      return_press()                                                 *
 *      return_enthalpy()                                              *
 *      return_fpe()                                                   *
 *      return_cp()                                                    *
 *      enthalpy()                                                     *
 *      blackbody_fraction()                                           *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*====================== thermo_setup() ======================================*/

/*
 * This function initializes the thermodynamics functions.
 * Adapted from Peter Gierasch's Fortran subroutines setup(),
 * numbers(), h2properties(), hydrogen(), theta2t_table() and trgrid().
 *
 * NOTE: cpr_out is the low-temperature-limit of cp/rgas.
 */

void thermo_setup(planetspec *planet,
                  EPIC_FLOAT *cpr_out)
{
  int
    i,ii,j,m,n,
    jmax = 50;
  int
    jn[2];
  EPIC_FLOAT
    xh2,xhe,x3,cpr,
    t0,p0,
    c1,c2,
    p,theta,temperature,
    ho,hp,ff,so,sp,
    pottempo,pottempp,
    den,term,y,thetaln;
  EPIC_FLOAT
    temp[MDIM_THERMO],
    tho[MDIM_THERMO],
    thp[MDIM_THERMO],
    tvector[MDIM_THERMO],
    thvector[MDIM_THERMO],
    aa[MDIM_THERMO],
    a[8],
    z[3][2],
    ndegeneracy[]={3.,1.};

  xh2      = planet->x_h2;
  xhe      = planet->x_he;
  x3       = planet->x_3;
  if (planet->x_h2 > 0) {
    *cpr_out = cpr = (CPRH2*xh2+CPRHE*xhe+CPR3*x3)/(xh2+xhe+x3);
  }
  else {
    *cpr_out = cpr = planet->cpr;
  }

  /* 
   * Calculate hydrogen (H_2) properties.
   * The subscript "o" refers to ortho-hydrogen, "p" to para-hydrogen.
   *
   * In Gierasch's Fortran code, this segment is contained in the subroutine
   * h2properties().
   */

  /* 
   * Reference temperature and pressure used to define the mean 
   * potential temperature.
   * See Dowling et al (1998), Appendix A.
   */
  t0 = 1.;
  if (planet->p0 < 1.e+5) {
    p0 = planet->p0;
  }
  else {
    p0 = 1.e+5; /* 1 bar */
  }

  /*
   * See Dowling et al (1998), eq. (A.12).
   */
  c1    = log(K_B*t0/p0*pow(2.*M_PI*(M_PROTON/H_PLANCK)*(K_B*t0/H_PLANCK),1.5));
  c2    = log(9.);
  p     = p0;
  theta = 87.567;

  /*
   * The array a[0:7] has entries
   *  a[0]:      equilibrium para fraction
   *  a[1],a[2]: ortho, para rotational internal energy per particle over K_B
   *  a[3],a[4]: ortho, para rotational cp per particle, units K_B
   *  a[5],a[6]: ortho, para rotational -Helmholtz free energy per particle over K_B*T
   *  a[7]:      equilibrium H2 (converting) cp per particle, units K_B
   */

  for (i = 0; i < MDIM_THERMO; i++) {
    temperature = 500.*(EPIC_FLOAT)(i+1)/(EPIC_FLOAT)MDIM_THERMO;
    if (temperature < 10.) {
      /*
       * Real hydrogen is not an ideal gas at this low-T limit.
       * These values are placeholders to avoid blow-ups during testing.
       */
      ho       = CPRH2*temperature;
      hp       = CPRH2*temperature;
      pottempo = temperature;
      pottempp = temperature;
      a[0]     = 1.;
      a[1]     = 175.1340;
      a[2]     = 0.;
      a[3]     = 0.;
      a[4]     = 0.;
      a[5]     = 0.;
      a[6]     = 0.; 
      a[7]     = 0.;
      ff       = 0.;
    }
    else {
      /*
       * In Gierasch's Fortran, this segment is contained in the subroutine
       * hydrogen().
       */
      y = theta/temperature;
      y = MIN(y,30.);
      for (n = 0; n < 2; n++) {
        for (m = 0; m < 3; m++) {
          z[m][n] = 0.;
        }
      }
      for (j = 1; j <= jmax; j++) {
        jn[0] = 2*j-1;
        jn[1] = jn[0]-1;
        for (n = 0; n < 2; n++) {
          term = ndegeneracy[n]*(2*jn[n]+1)*exp(-jn[n]*(jn[n]+1)*y);
          for (m = 0; m < 3; m++) {
            z[m][n] += term;
            if (m < 2) {
              term *= jn[n]*(jn[n]+1);
            }
          }
        }
        if (j > 1 && term < 1.e-20) break;
      }
      den  = z[0][0]+z[0][1];

      a[0] = z[0][1]/den;
      for (n = 0; n < 2; n++) {
        a[n+1] = theta*z[1][n]/z[0][n]; 
        a[n+3] = y*y*(z[0][n]*z[2][n]-z[1][n]*z[1][n])/(z[0][n]*z[0][n]);
        a[n+5] = log(z[0][n]);
      }
      a[7] = (1.-a[0])*a[3]+a[0]*a[4]
            +(a[2]-a[1])*y/temperature*
                     (z[1][1]*z[0][0]-z[0][1]*z[1][0])/(den*den);
      /*
       * End of segment contained in Gierasch's Fortran subroutine hydrogen().
       */
      ho = a[1]+CPRH2*temperature-2.*theta;
      hp = a[2]+CPRH2*temperature;
      /*
       * entropies normalized at p0, T-->0, per particle divided by K_B
       */
      so = -log(p/p0)+2.5*log(temperature)
                     +1.5*M_LN2+c1+(ho+2.*theta)/temperature+a[5];
      sp = -log(p/p0)+2.5*log(temperature)
                     +1.5*M_LN2+c1+hp/temperature+a[6];
      /*
       * potential temperatues, equal T as T-->0
       */
       pottempo = exp(0.4*(so-c2-1.5*M_LN2-c1-2.5));
       pottempp = exp(0.4*(sp   -1.5*M_LN2-c1-2.5));
      /*
       * curly F, equals -free energy difference, normalized at T=0
       */
      ff = -(so-c2-1.5*M_LN2-c1-2.5-ho/temperature)
            +(sp-1.5*M_LN2-c1-2.5-hp/temperature);
      ff *= temperature;
    }

    /*
     * Save T, ortho and para enthalpies (offset so h(T=0)=0), ortho and
     * para entropies, ortho and para potential temperatures, and curly F.
     * Units are per particle, divided by Boltzmann constant, K_B.  Potential
     * temperatures and curly F are degrees K and degrees K per particle
     * over K_B.
     */
    temp[i]= temperature;
    tho[i] = pottempo;
    thp[i] = pottempp;

    thermo.array[0][i]       = ho;
    thermo.array[1][i]       = hp;
    thermo.array[2][i]       = ff;
    thermo.array[3][i]       = a[0];
    thermo.array[4][i]       = a[1]-a[2];
    thermo.t_grid[i]         = temperature;
    thermo.theta_array[0][i] = pottempo;
    thermo.theta_array[1][i] = pottempp;
  }
  /*
   * End of segment contained in Gierasch's Fortran subroutine h2properties().
   */

  /*
   * In Gierasch's Fortran, this segment is contained in the subroutine
   * theta2t_table().
   */

  for (m = 0; m < MDIM_THERMO; m++) {
    thermo.theta_grid[m] = (THLO_THERMO*(EPIC_FLOAT)(MDIM_THERMO-(m+1))
                           +THHI_THERMO*(EPIC_FLOAT)(m))/(EPIC_FLOAT)(MDIM_THERMO-1);
  }
  for (n = 0; n < NDIM_THERMO; n++) {
    thermo.fpdat[n] = ((EPIC_FLOAT)(n))/(EPIC_FLOAT)(NDIM_THERMO-1);
    for (i = 0; i < MDIM_THERMO; i++) {
      thetaln = ( planet->x_h2*CPRH2*((1.-thermo.fpdat[n])*log(tho[i])
                                     +(   thermo.fpdat[n])*log(thp[i]))
                +(planet->x_he*CPRHE
                 +planet->x_3*CPR3 )*log(temp[i]) )/cpr;
      thvector[i] = exp(thetaln);
      tvector[i]  = temp[i];
    }
    /*
     * In Gierasch's Fortran, this segment is contained in the subroutine
     * trgrid().
     */
    for (j = 0; j < MDIM_THERMO; j++) {
      aa[j] = tvector[j];
    }
    tvector[MDIM_THERMO-1] = aa[MDIM_THERMO-1];
    for (i = 0; i < MDIM_THERMO; i++) {
      for (j = 1; j < MDIM_THERMO; j++) {
        if (thvector[j] >= thermo.theta_grid[i]) {
          break;
        }
      }
      tvector[i] = aa[j-1]+(aa[j]-aa[j-1])*
                         (thermo.theta_grid[i]-thvector[j-1])/
                         (thvector[j]         -thvector[j-1]);
    }
    /*
     * End of segment contained in Gierasch's Fortran subroutine trgrid().
     */

    for (m = 0; m < MDIM_THERMO; m++) {
      thermo.t[n][m] = tvector[m];
    }
  }
  /*
   * End of segment contained in Gierasch's Fortran subroutine theta2t_table().
   */

  return;
}

/*====================== end of thermo_setup() ===============================*/

/*====================== return_temp() =======================================*/

/*
 * MAX_IT stands for "maximum number of iterations."
 */
#undef  MAX_IT
#define MAX_IT 10

/*
 * Adapted from Gierasch's Fortran subroutine get_temperature().
 * See Dowling et al (1998), Appendix A, eq. (A.15).
 * Table resolution has been improved since the Dowling et al paper.
 *
 * NOTE:  The calling program must initialize this function with a call
 *        to thermo_setup().
 */

EPIC_FLOAT return_temp(planetspec *planet,
                       EPIC_FLOAT  fp,
                       EPIC_FLOAT  p,
                       EPIC_FLOAT  theta)
{
  int
    m,n,it,
    error_flag;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    p0,kappa;
  EPIC_FLOAT
    temperature,
    theta1,
    t1,t2,ttol,
    em,en,
    fract_fp,fract_theta;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="return_temp";

#if EPIC_CHECK == 1
  /* Sanity checks: */
  if (p <= 0.) {
    sprintf(Message,"p=%e <= 0",p);
    epic_error(dbmsname,Message);
  }
  if (theta <= 0.) {
    sprintf(Message,"theta=%e <= 0",theta);
    epic_error(dbmsname,Message);
  }
#endif

  if (!initialized) {
    kappa = planet->kappa;
    if (planet->p0 < 1.e+5) {
      p0 = planet->p0;
    }
    else {
      p0 = 1.e+5; /* 1 bar */
    }
    initialized = TRUE;
  }

  if (planet->x_h2 == 0.) {
    /*
     * No hydrogen, so use the standard relationship between
     * pressure, temperature, and theta.
     */
    temperature = theta*pow(p/p0,kappa);
  }
  else {

#if EPIC_CHECK == 1
    /* Sanity checks: */
    if (fp < 0. || fp > 1.) {
      sprintf(Message,"fp=%e",fp);
      epic_error(dbmsname,Message);
    }
#endif

    theta1 = theta*pow(p/p0,kappa);
  
    if (!isfinite(theta1)) {
      sprintf(Message,"theta=%g K; pressure=%g hPa, fp=%g",theta,p/100.,fp);
      epic_error(dbmsname,Message);
    }
    else if (theta1 <= THLO_THERMO) {
      temperature = theta1;
    }
    else if (theta1 >= THHI_THERMO) {
      temperature = exp( (planet->cpr*log(theta1)
                         -CPRH2*planet->x_h2*((1.-fp)*CCOLN_THERMO+fp*CCPLN_THERMO))/
                         (CPRHE*planet->x_he+CPR3*(1.-planet->x_he)) );
    }
    else {
      /* 0. < en < NDIM_THERMO-1 */
      en = (EPIC_FLOAT)(NDIM_THERMO-1)*
                (fp                         -thermo.fpdat[0])/
                (thermo.fpdat[NDIM_THERMO-1]-thermo.fpdat[0]);
      n  = (int)en;
      if (n > NDIM_THERMO-2) {
        /* 0 < n < nmax-1 */
        n        = NDIM_THERMO-2;
        fract_fp = 1.;
      }
      else if (n < 0) {
        n        = 0;
        fract_fp = 0.;
      }
      else {
        /* 0. < fract_fp < 1. */
        fract_fp = fmod(en,1.);
      }

      em = (EPIC_FLOAT)(MDIM_THERMO-1)*
                (theta1                          -thermo.theta_grid[0])/
                (thermo.theta_grid[MDIM_THERMO-1]-thermo.theta_grid[0]);
      m  = (int)em;
      if (m > MDIM_THERMO-2) {
        m           = MDIM_THERMO-2;
        fract_theta = 1.;
      }
      else if (m < 0) {
        m           = 0;
        fract_theta = 0.;
      }
      else {
        fract_theta = fmod(em,1.);
      }
      /*
       * NOTE: Using bilinear interpolation.  It may be possible to improve
       *       the accuracy with a more sophisticated two-variable interpolation scheme.
       */
      temperature = thermo.t[n  ][m  ]*(1.-fract_theta)*(1.-fract_fp)
                   +thermo.t[n+1][m  ]*(1.-fract_theta)*(   fract_fp)
                   +thermo.t[n  ][m+1]*(   fract_theta)*(1.-fract_fp)
                   +thermo.t[n+1][m+1]*(   fract_theta)*(   fract_fp);
    }
  }

  if (strcmp(grid.eos,"ideal") == 0) {
    return temperature;
  }
  else if (strcmp(grid.eos,"virial") == 0) {
      /*
       * Iterate to get temperature that satisfies
       * theta-return_theta(temperature) = 0.
       */
    THMTH_fp     = fp;
    THMTH_p      = p;
    THMTH_theta  = theta;
    THMTH_planet = planet;

    /* Initial guess: */
    ttol  = pow(machine_epsilon(),2./3.);
    t1    = temperature*0.9;
    t2    = temperature*1.1;

    for (it = 0; it < MAX_IT; it++) {
      error_flag = find_root(t1,t2,ttol,&temperature,th_minus_th_t);
      if (error_flag == 0) {
        /* Convergence */
        return temperature;
      }
      /* Try a wider interval. */
      t1 *= .5;
      t2 *= 2.;
    }

    sprintf(Message,"exceeded MAX_IT = %d",MAX_IT);
    epic_error(dbmsname,Message);
  }
  else {
    sprintf(Message,"unrecognized grid.eos: %s",grid.eos);
    epic_error(dbmsname,Message);
  }

  /* Should never get here.*/
  sprintf(Message,"should never get here");
  epic_error(dbmsname,Message);
  return 0.;
}

/*======================= end of return_temp() ==============================*/

/*======================= alt_return_temp() =================================*/

/*
 * MAX_IT stands for "maximum number of iterations."
 */
#undef  MAX_IT
#define MAX_IT 10

EPIC_FLOAT alt_return_temp(planetspec *planet,
                           EPIC_FLOAT  fp,
                           EPIC_FLOAT  p,
                           EPIC_FLOAT  mu,
                           EPIC_FLOAT  density)
{
  int
    it,
    error_flag;
  EPIC_FLOAT
    temperature,
    t1,t2,ttol;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="alt_return_temp";

  temperature = (p*mu)/(density*R_GAS);

  if (strcmp(grid.eos,"ideal") == 0) {
    return temperature;
  }
  else if (strcmp(grid.eos,"virial") == 0) {
    /*
     * Iterate to get temperature that satisfies
     * density-return_density(temperature) = 0.
     */
    RHOMRHO_fp      = fp;
    RHOMRHO_p       = p;
    RHOMRHO_mu      = mu;
    RHOMRHO_density = density;
    RHOMRHO_planet  = planet;

    ttol        = pow(machine_epsilon(),2./3.);
    t1          = temperature*0.9;
    t2          = temperature*1.1;

    for (it = 0; it < MAX_IT; it++) {
      error_flag = find_root(t1,t2,ttol,&temperature,rho_minus_rho);
      if (error_flag == 0) {
        /* Convergence */
        return temperature;
      }
      /* Try a wider interval. */
      t1 *= .5;
      t2 *= 2.;
    }

    sprintf(Message,"exceeded MAX_IT = %d",MAX_IT);
    epic_error(dbmsname,Message);
  }

  /* Should never get here.*/
  sprintf(Message,"should never get here");
  epic_error(dbmsname,Message);
  return 0.;
}

/*======================= end of alt_return_temp() ==========================*/

/*======================= rho_minus_rho() ===================================*/

/*
 * For use with find_root().
 */

EPIC_FLOAT rho_minus_rho(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    ans;

  ans = RHOMRHO_density-return_density(RHOMRHO_planet,
                                       RHOMRHO_fp,
                                       RHOMRHO_p,
                                       temperature,
                                       RHOMRHO_mu,
                                       PASSING_T);
  return ans;
}

/*======================= end of rho_minus_rho_p() ==========================*/


/*======================= return_density() ==================================*/

EPIC_FLOAT return_density(planetspec *planet, 
                          EPIC_FLOAT  fp,
                          EPIC_FLOAT  p,
                          EPIC_FLOAT  theta,
                          EPIC_FLOAT  mu,
                          int         temp_type)
{
  EPIC_FLOAT 
    temperature,
    density,
    b,z_comp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="return_density";

  if (temp_type == PASSING_THETA) {
    temperature = return_temp(planet,fp,p,theta);
  }
  else if (temp_type == PASSING_T) {
    temperature = theta;
  }
  else {
    sprintf(Message,"unknown temp_type = %d",temp_type);
    epic_error(dbmsname,Message);
  }

  density = p*mu/(R_GAS*temperature);

  if (strcmp(grid.eos,"virial") == 0) {
    /* 
     * Make non-ideal equation of state correction:
     */
    b        = sum_xx(planet,b_vir,temperature);
    z_comp   = 1.+b*p;
    density /= z_comp;
  }

  return density;
}

/*======================= end of return_density() ===========================*/

/*======================= p_from_t_rho_mu() =================================*/

EPIC_FLOAT p_from_t_rho_mu(planetspec *planet,
                           EPIC_FLOAT  temperature,
                           EPIC_FLOAT  rho,
                           EPIC_FLOAT  mu)
{
  EPIC_FLOAT
    p;

  p = rho*(R_GAS/mu)*temperature;

  if (strcmp(grid.eos,"virial") == 0) {
    /* 
     * Make non-ideal equation of state correction:
     */
    p /= 1.-p*sum_xx(planet,b_vir,temperature);
  }

  return p;
}

/*======================= end of p_from_t_rho_mu() ==========================*/

/*======================= return_theta() ====================================*/

/*
 * Adapted from Peter Gierasch's Fortran subroutine get_theta().
 *
 * NOTE:  The calling program must initialize this function with a call
 *        to thermo_setup().
 */

EPIC_FLOAT return_theta(planetspec *planet,
                        EPIC_FLOAT  fp,
                        EPIC_FLOAT  p,
                        EPIC_FLOAT  temperature,
                        EPIC_FLOAT *theta_ortho,
                        EPIC_FLOAT *theta_para)
{
  int
    j,m;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    p0,kappa;
  EPIC_FLOAT
    b,b1,tmp,
    theta,thetaln,
    cc,tt,pp,
    em,fract;
  EPIC_FLOAT
    thermo_vector[2];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="return_theta";

  if (!initialized) {
    kappa = planet->kappa;
    if (planet->p0 < 1.e+5) {
      p0 = planet->p0;
    }
    else {
      p0 = 1.e+5; /* 1 bar */
    }
    initialized = TRUE;
  }

  if (planet->x_h2 == 0.) {
    /*
     * No hydrogen, so use standard definition of theta.
     */
    theta = temperature*pow(p0/p,kappa);
  }
  else {
    /*
     * Use mean theta as defined in Dowling et al (1998), to handle ortho/para hydrogen.
     */
    if (!isfinite(temperature)) {
      sprintf(Message,"temperature=%g K; pressure=%g hPa, fp=%g",temperature,p/100.,fp);
      epic_error(dbmsname,Message);
    }
    else if (temperature <= 20.) {
      theta        = temperature;
      *theta_ortho = temperature;
      *theta_para  = temperature;
    }
    else if (temperature > 500.) {
      cc           = planet->x_h2*2.5*((1.-fp)*CCOLN_THERMO+fp*CCPLN_THERMO);
      theta        = exp(cc/planet->cpr)*
                      pow(temperature,(( 3.5*planet->x_h2  /* 3.5 since high T */
                                        +2.5*planet->x_he
                                        +3.5*planet->x_3 )/planet->cpr));
      tt           = pow(temperature,3.5/2.5);
      *theta_ortho = 0.12175*tt;
      *theta_para  = 0.18892*tt;
    }
    else {
      /* 0 < em < MDIM_THERMO-1 */
      em = (EPIC_FLOAT)(MDIM_THERMO-1)*
                         (temperature                 -thermo.t_grid[0])/
                         (thermo.t_grid[MDIM_THERMO-1]-thermo.t_grid[0]);
      m  = (int)em;
      /*  0 < m < MDIM_THERMO-2 */
      if (m == MDIM_THERMO-1) {
        m     -= 1;
        fract  = 1.;
      }
      else {
        fract = fmod(em,1.);
      }
      for (j = 0; j < 2; j++) {
        thermo_vector[j] = (1.-fract)*thermo.theta_array[j][m  ]
                          +(   fract)*thermo.theta_array[j][m+1];
      }

      thetaln = (planet->x_h2)*( (1.-fp)*log(thermo_vector[0])
                                +(   fp)*log(thermo_vector[1]) )
               +(planet->x_he*2.5+planet->x_3*3.5)*log(temperature)/planet->cpr;

      theta        = exp(thetaln);
      *theta_ortho = thermo_vector[0];
      *theta_para  = thermo_vector[1];
    }
    pp            = pow(p0/p,kappa);
    theta        *= pp;
    *theta_ortho *= pp;
    *theta_para  *= pp;
  }

  if (strcmp(grid.eos,"virial") == 0) {
    /* 
     * Make non-ideal equation of state corrections.
     *
     * NOTE: Need to check validity of these formulas.
     *       The formula in Dymon and Smith (1980) on p.x is confusing,
     *       and the implementation below may be in error.
     */
    kappa         = planet->kappa;
    b             = sum_xx(planet,b_vir, temperature);
    b1            = sum_xx(planet,b1_vir,temperature);
    tmp           = exp(-p*(b+b1)*kappa);
    theta        *= tmp;

    kappa         = planet->kappa*R_GAS/(2.016*planet->rgas);
    b             = b_vir( "H_2","H_2",temperature);
    b1            = b1_vir("H_2","H_2",temperature);
    tmp           = exp(-p*(b+b1)*kappa);
    *theta_ortho *= tmp;
    *theta_para  *= tmp;
  }

  return theta;
}

/*======================= end of return_theta() =============================*/

/*======================= return_press() ====================================*/

/*
 * MAX_IT stands for "maximum number of iterations."
 */
#undef  MAX_IT
#define MAX_IT 10

/*
 * NOTE:  The calling program must initialize this function with a call
 *        to thermo_setup().
 */

EPIC_FLOAT return_press(planetspec *planet,
                        EPIC_FLOAT  fp,
                        EPIC_FLOAT  temperature,
                        EPIC_FLOAT  theta)
{
  int
    it,
    error_flag;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    p0;
  EPIC_FLOAT
    press,p1,p2,
    ptol,p_root;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="return_press";

#if EPIC_CHECK == 1
  /* Sanity checks: */
  if (temperature <= 0.) {
    sprintf(Message,"temperature = %e <= 0",temperature);
    epic_error(dbmsname,Message);
  }
  if (theta <= 0.) {
    sprintf(Message,"theta = %e <= 0",theta);
    epic_error(dbmsname,Message);
  }
#endif

  if (!initialized) {
    if (planet->p0 < 1.e+5) {
      p0 = planet->p0;
    }
    else {
      p0 = 1.e+5; /* 1 bar */
    }
    initialized = TRUE;
  }

  if (planet->x_h2 == 0.) {
    /*
     * No hydrogen, so use the standard relationship between 
     * pressure, temperature, and theta.
     */
    press = p0*pow(temperature/theta,planet->cpr);
    return press;
  }
  else {
    /*
     * Iterate to get pressure that satisfies theta-return_theta(pressure) = 0.
     */

#if EPIC_CHECK == 1
    /* Sanity checks: */
    if (fp <= 0.) {
      sprintf(Message,"fp = %e <= 0",fp);
      epic_error(dbmsname,Message);
    }
#endif

    THMTH_fp          = fp;
    THMTH_temperature = temperature;
    THMTH_theta       = theta;
    THMTH_planet      = planet;

    /* Initial guess: */
    press = p0*pow(temperature/theta,planet->cpr);
    ptol  = pow(machine_epsilon(),2./3.);
    p1    = press*0.9;
    p2    = press*1.1;

    for (it = 0; it < MAX_IT; it++) {
      error_flag = find_root(p1,p2,ptol,&p_root,th_minus_th_p);
      if (error_flag == 0) {
        /* 
         * Convergence.
         */
        return p_root;
      }
      /* Try a wider interval. */
      p1 *= .5;
      p2 *= 2.;
    }

    sprintf(Message,"exceeded MAX_IT = %d",MAX_IT);
    epic_error(dbmsname,Message);
  }

  /* Should never get here.*/
  sprintf(Message,"should never get here");
  epic_error(dbmsname,Message);
  return 0.;
}

/*======================= end of return_press() =============================*/

/*======================= th_minus_th_p() ===================================*/

/*
 * For use with find_root().
 */

EPIC_FLOAT th_minus_th_p(EPIC_FLOAT p)
{
  EPIC_FLOAT
    theta_ortho,theta_para,
    ans;

  ans = THMTH_theta-return_theta(THMTH_planet,
                                 THMTH_fp,
                                 p,
                                 THMTH_temperature,
                                 &theta_ortho,&theta_para);
  return ans;
}

/*======================= end of th_minus_th_p() ============================*/

/*======================= th_minus_th_t() ===================================*/

/*
 * For use with find_root().
 */

EPIC_FLOAT th_minus_th_t(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    theta_ortho,theta_para,
    ans;

  ans = THMTH_theta-return_theta(THMTH_planet,
                                 THMTH_fp,
                                 THMTH_p,
                                 temperature,
                                 &theta_ortho,&theta_para);
  return ans;
}

/*======================= end of th_minus_th_t() ============================*/

/*======================= return_enthalpy() =================================*/

EPIC_FLOAT return_enthalpy(planetspec *planet,
                           EPIC_FLOAT  fp,
                           EPIC_FLOAT  pressure,
                           EPIC_FLOAT  temperature,
                           EPIC_FLOAT *fgibb,
                           EPIC_FLOAT *fpe,
                           EPIC_FLOAT *uoup)
{
  /*
   * Adapted from Peter Gierasch's Fortran subroutine get_enthalpy().
   *
   * NOTE:  The calling program must initialize this function with a call
   *        to thermo_setup().
   */
  int
    j,m;
  EPIC_FLOAT
    b,b1,em,
    rgas,
    ho,hp,enthalpy,
    fract;
  EPIC_FLOAT
    thermo_vector[5];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="return_enthalpy";

  if (planet->x_h2 == 0.) {
    /*
     * No hydrogen, so assume cp is constant and set enthalpy = cp*T.
     */
    enthalpy = planet->cpr*temperature;
    *fgibb   = 0.;
    *fpe     = 0.;
    *uoup    = 0.;
  }
  else {
    if (!isfinite(temperature)) {
      sprintf(Message,"temperature=%g K; pressure=%g hPa, fp=%g",temperature,pressure/100.,fp);
      epic_error(dbmsname,Message);
    }
    else if (temperature <= 20.) {
      enthalpy = planet->cpr*temperature;
      *fgibb   = 0.;
      *fpe     = 1.;
      *uoup    = 175.1340;
    }
    else if (temperature > 500.) {
      ho       = 1545.3790+3.5*(temperature-500.);
      hp       = 1720.3776+3.5*(temperature-500.);
      /*
       * NOTE: Should replace "planet->x_3" with a loop
       *       over condensable species.
       */
      enthalpy = (planet->x_h2)*((1.-fp)*ho+fp*hp)
                +(planet->x_he*2.5+planet->x_3*3.5)*temperature;
      *fgibb   = planet->x_h2*2.5*(CCPLN_THERMO-CCOLN_THERMO)*temperature
                    -planet->x_h2*(hp-ho);
      *fpe     = 0.25;
      *uoup    = 0.;
    }
    else {
      /* 0 < em < MDIM_THERMO-1 */
      em = (EPIC_FLOAT)(MDIM_THERMO-1)*
                 (temperature                 -thermo.t_grid[0])/
                 (thermo.t_grid[MDIM_THERMO-1]-thermo.t_grid[0]);
      m = (int)em;
      /* 0 < m < MDIM_THERMO-2 */
      if (m == MDIM_THERMO-1) {
        m--;
        fract = 1.;
      }
      else {
        fract = fmod(em,1.);
      }
      for (j = 0; j < 5; j++) {
        thermo_vector[j] = (1.-fract)*thermo.array[j][m  ]
                          +(   fract)*thermo.array[j][m+1];
      }
      enthalpy = (planet->x_h2)*((1.-fp)*thermo_vector[0]
                               +(    fp)*thermo_vector[1])
                +(planet->x_he*2.5+planet->x_3*3.5)*temperature;
      *fgibb   = planet->x_h2*thermo_vector[2];
      *fpe     = thermo_vector[3];
      *uoup    = thermo_vector[4];
    }
  }

  if (strcmp(grid.eos,"virial") == 0) {
    /* 
     * Make non-ideal equation of state corrections. Use enthalpy-correction 
     * equation on p.x of Dymond and Smith (1980); see also p.xiv.
     *
     * The quantities fgibb and uoup are differences between ortho and para
     * hydrogen.  Since we are not distinguishing these in the non-ideal 
     * equation of state, we make no corrections to fgibb and uoup.
     */
    b         = sum_xx(planet,b_vir, temperature);
    b1        = sum_xx(planet,b1_vir,temperature);
    enthalpy += pressure*temperature*(b-b1);
  }

  rgas      = planet->rgas;
  enthalpy *= rgas;
  *fgibb   *= rgas;
  *uoup    *= rgas;

  return enthalpy;
}

/*======================= end of return_enthalpy() ==========================*/

/*======================= return_fpe() ======================================*/

/*
* Calculate equilibrium fraction of para hydrogen.
* Adapted from Peter Gierasch's hydrogen() Fortran subroutine.
*/

EPIC_FLOAT return_fpe(EPIC_FLOAT temperature) {
  int
    n,j;
  double
    y,term,
    z[2],jn[2],
    ndegen[2] = {3.,1.};
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="return_fpe";

  if (temperature > 800.) return 0.25;

  y = MIN(87.567/temperature,30.);

  z[0] = z[1] = 0.;
  for (j = 1; j <= 12; j++) {
    jn[0] = 2.*(double)j-1.;
    jn[1] = jn[0]-1.;
    for (n = 0; n < 2; n++) {
      term  = ndegen[n]*(2.*jn[n]+1.)*exp(-jn[n]*(jn[n]+1.)*y);
      z[n] += term;
    }
    if (j > 1 && term < 1.e-20) break;
  }

  return (EPIC_FLOAT)(z[1]/(z[0]+z[1]));
}

/*======================= end of return_fpe() ===============================*/

/*======================= return_cp() =======================================*/

/*  
 * Returns the specific heat at constant pressure. This is calculated as a 
 * derivative of the enthalpy, cp = (denthalpy/dT)_p, a partial derivative  
 * at constant p. All other state variables, such as fpara, water amount, etc,
 * are also to be held constant.  By using this method, cp will be correct
 * even if there are changes to the thermodynamics in return_enthalpy().
 *   -- A.P. Showman, 8/31/1999.
 * 
 * NOTE: thermo_setup() must have already been called at initialization.
 */

EPIC_FLOAT return_cp(planetspec *planet,
                     EPIC_FLOAT  fp,
                     EPIC_FLOAT  p,
                     EPIC_FLOAT  temp)
{
  EPIC_FLOAT
    cp,h1,h2,
    deltaT,     
    fgibb,fpe,uoup,
    epsilon = 1.e-6;

  /*
   * Handle special cases.
   */
  if (strcmp(planet->name,"Held_Suarez") == 0) {
    cp = 1004.;
  }
  else if (strcmp(planet->name,"Goldmine") == 0) {
    /*  Goldmine is a testbed case  */
    cp = 1004.;
  }
  else if (strcmp(planet->name,"Mars") == 0) {
    return isobaric_specific_heat(temp,CO_2_INDEX);
  }
  else {
    /*
     * Handle general case.
     */   
    deltaT = temp*epsilon;
    h2     = return_enthalpy(planet,fp,p,temp+deltaT,&fgibb,&fpe,&uoup);
    h1     = return_enthalpy(planet,fp,p,temp-deltaT,&fgibb,&fpe,&uoup);
    cp     = (h2-h1)/(2.*deltaT);
  }

  return cp;
}

/*======================= end of return_cp() ================================*/

/*======================= isobaric_specific_heat() ==========================*/

/*
 * Return cp [J kg-1 K-1] for the given temperature and species.
 */
EPIC_FLOAT isobaric_specific_heat(EPIC_FLOAT temp,
                                  int        index)
{
  EPIC_FLOAT
    ans;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="isobaric_specific_heat";

  switch (index) {
    case CO_2_INDEX:{
      /*
       * Data from H.W. Woolley (1954) "Thermodynamic functions for carbon dioxide
       * in the ideal gas state", J. Res. Nat. Bureau Stand. 52, 289-292.
       */
      EPIC_FLOAT
        t_d;
      const int
        ndat = 141;
      static int
        i = -2,
        initialized = FALSE;
      static float_triplet
        *CO2_table;


      if (!initialized) {
        /* Allocate memory */
        CO2_table = ftriplet(0,ndat-1,dbmsname);

        CO2_table[  0].x =   50., CO2_table[  0].y = 3.5001;
        CO2_table[  1].x =   60., CO2_table[  1].y = 3.5002;
        CO2_table[  2].x =   70., CO2_table[  2].y = 3.5006;
        CO2_table[  3].x =   80., CO2_table[  3].y = 3.5020;
        CO2_table[  4].x =   90., CO2_table[  4].y = 3.5055;

        CO2_table[  5].x =  100., CO2_table[  5].y = 3.5128;
        CO2_table[  6].x =  110., CO2_table[  6].y = 3.5249;
        CO2_table[  7].x =  120., CO2_table[  7].y = 3.5432;
        CO2_table[  8].x =  130., CO2_table[  8].y = 3.5680;
        CO2_table[  9].x =  140., CO2_table[  9].y = 3.5995;

        CO2_table[ 10].x =  150., CO2_table[ 10].y = 3.6372;
        CO2_table[ 11].x =  160., CO2_table[ 11].y = 3.6804;
        CO2_table[ 12].x =  170., CO2_table[ 12].y = 3.7282;
        CO2_table[ 13].x =  180., CO2_table[ 13].y = 3.7800;
        CO2_table[ 14].x =  190., CO2_table[ 14].y = 3.8347;

        CO2_table[ 15].x =  200., CO2_table[ 15].y = 3.8916;
        CO2_table[ 16].x =  210., CO2_table[ 16].y = 3.9502;
        CO2_table[ 17].x =  220., CO2_table[ 17].y = 4.0097;
        CO2_table[ 18].x =  230., CO2_table[ 18].y = 4.0695;
        CO2_table[ 19].x =  240., CO2_table[ 19].y = 4.1296;

        CO2_table[ 20].x =  250., CO2_table[ 20].y = 4.1892;
        CO2_table[ 21].x =  260., CO2_table[ 21].y = 4.2484;
        CO2_table[ 22].x =  270., CO2_table[ 22].y = 4.3068;
        CO2_table[ 23].x =  280., CO2_table[ 23].y = 4.3643;
        CO2_table[ 24].x =  290., CO2_table[ 24].y = 4.4208;

        CO2_table[ 25].x =  300., CO2_table[ 25].y = 4.4763;
        CO2_table[ 26].x =  310., CO2_table[ 26].y = 4.5307;
        CO2_table[ 27].x =  320., CO2_table[ 27].y = 4.5840;
        CO2_table[ 28].x =  330., CO2_table[ 28].y = 4.6361;
        CO2_table[ 29].x =  340., CO2_table[ 29].y = 4.6871;

        CO2_table[ 30].x =  350., CO2_table[ 30].y = 4.7371;
        CO2_table[ 31].x =  360., CO2_table[ 31].y = 4.7859;
        CO2_table[ 32].x =  370., CO2_table[ 32].y = 4.8335;
        CO2_table[ 33].x =  380., CO2_table[ 33].y = 4.8801;
        CO2_table[ 34].x =  390., CO2_table[ 34].y = 4.9257;

        CO2_table[ 35].x =  400., CO2_table[ 35].y = 4.9704;
        CO2_table[ 36].x =  410., CO2_table[ 36].y = 5.0140;
        CO2_table[ 37].x =  420., CO2_table[ 37].y = 5.0566;
        CO2_table[ 38].x =  430., CO2_table[ 38].y = 5.0983;
        CO2_table[ 39].x =  440., CO2_table[ 39].y = 5.1392;

        CO2_table[ 40].x =  450., CO2_table[ 40].y = 5.1792;
        CO2_table[ 41].x =  460., CO2_table[ 41].y = 5.2183;
        CO2_table[ 42].x =  470., CO2_table[ 42].y = 5.2566;
        CO2_table[ 43].x =  480., CO2_table[ 43].y = 5.2942;
        CO2_table[ 44].x =  490., CO2_table[ 44].y = 5.3310;

        CO2_table[ 45].x =  500., CO2_table[ 45].y = 5.3671;
        CO2_table[ 46].x =  510., CO2_table[ 46].y = 5.4024;
        CO2_table[ 47].x =  520., CO2_table[ 47].y = 5.4371;
        CO2_table[ 48].x =  530., CO2_table[ 48].y = 5.4711;
        CO2_table[ 49].x =  540., CO2_table[ 49].y = 5.5044;

        CO2_table[ 50].x =  550., CO2_table[ 50].y = 5.5371;
        CO2_table[ 51].x =  560., CO2_table[ 51].y = 5.5691;
        CO2_table[ 52].x =  570., CO2_table[ 52].y = 5.6006;
        CO2_table[ 53].x =  580., CO2_table[ 53].y = 5.6315;
        CO2_table[ 54].x =  590., CO2_table[ 54].y = 5.6618;

        CO2_table[ 55].x =  600., CO2_table[ 55].y = 5.6915;
        CO2_table[ 56].x =  610., CO2_table[ 56].y = 5.7207;
        CO2_table[ 57].x =  620., CO2_table[ 57].y = 5.7494;
        CO2_table[ 58].x =  630., CO2_table[ 58].y = 5.7775;
        CO2_table[ 59].x =  640., CO2_table[ 59].y = 5.8052;

        CO2_table[ 60].x =  650., CO2_table[ 60].y = 5.8324;
        CO2_table[ 61].x =  660., CO2_table[ 61].y = 5.8591;
        CO2_table[ 62].x =  670., CO2_table[ 62].y = 5.8853;
        CO2_table[ 63].x =  680., CO2_table[ 63].y = 5.9110;
        CO2_table[ 64].x =  690., CO2_table[ 64].y = 5.9363;

        CO2_table[ 65].x =  700., CO2_table[ 65].y = 5.9611;
        CO2_table[ 66].x =  710., CO2_table[ 66].y = 5.9855;
        CO2_table[ 67].x =  720., CO2_table[ 67].y = 6.0094;
        CO2_table[ 68].x =  730., CO2_table[ 68].y = 6.0329;
        CO2_table[ 69].x =  740., CO2_table[ 69].y = 6.0559;

        CO2_table[ 70].x =  750., CO2_table[ 70].y = 6.0786;
        CO2_table[ 71].x =  760., CO2_table[ 71].y = 6.1009;
        CO2_table[ 72].x =  770., CO2_table[ 72].y = 6.1228;
        CO2_table[ 73].x =  780., CO2_table[ 73].y = 6.1442;
        CO2_table[ 74].x =  790., CO2_table[ 74].y = 6.1653;

        CO2_table[ 75].x =  800., CO2_table[ 75].y = 6.1860;
        CO2_table[ 76].x =  810., CO2_table[ 76].y = 6.2064;
        CO2_table[ 77].x =  820., CO2_table[ 77].y = 6.2264;
        CO2_table[ 78].x =  830., CO2_table[ 78].y = 6.2460;
        CO2_table[ 79].x =  840., CO2_table[ 79].y = 6.2653;

        CO2_table[ 80].x =  850., CO2_table[ 80].y = 6.2843;
        CO2_table[ 81].x =  860., CO2_table[ 81].y = 6.3029;
        CO2_table[ 82].x =  870., CO2_table[ 82].y = 6.3212;
        CO2_table[ 83].x =  880., CO2_table[ 83].y = 6.3392;
        CO2_table[ 84].x =  890., CO2_table[ 84].y = 6.3569;

        CO2_table[ 85].x =  900., CO2_table[ 85].y = 6.3742;
        CO2_table[ 86].x =  910., CO2_table[ 86].y = 6.3913;
        CO2_table[ 87].x =  920., CO2_table[ 87].y = 6.4080;
        CO2_table[ 88].x =  930., CO2_table[ 88].y = 6.4244;
        CO2_table[ 89].x =  940., CO2_table[ 89].y = 6.4406;

        CO2_table[ 90].x =  950., CO2_table[ 90].y = 6.4565;
        CO2_table[ 91].x =  960., CO2_table[ 91].y = 6.4721;
        CO2_table[ 92].x =  970., CO2_table[ 92].y = 6.4874;
        CO2_table[ 93].x =  980., CO2_table[ 93].y = 6.5025;
        CO2_table[ 94].x =  990., CO2_table[ 94].y = 6.5173;

        CO2_table[ 95].x = 1000., CO2_table[ 95].y = 6.5318;
        CO2_table[ 96].x = 1050., CO2_table[ 96].y = 6.601;
        CO2_table[ 97].x = 1100., CO2_table[ 97].y = 6.664;
        CO2_table[ 98].x = 1150., CO2_table[ 98].y = 6.723;
        CO2_table[ 99].x = 1200., CO2_table[ 99].y = 6.776;

        CO2_table[100].x = 1250., CO2_table[100].y = 6.826;
        CO2_table[101].x = 1300., CO2_table[101].y = 6.872;
        CO2_table[102].x = 1350., CO2_table[102].y = 6.913;
        CO2_table[103].x = 1400., CO2_table[103].y = 6.952;
        CO2_table[104].x = 1450., CO2_table[104].y = 6.988;

        CO2_table[105].x = 1500., CO2_table[105].y = 7.021;
        CO2_table[106].x = 1600., CO2_table[106].y = 7.082;
        CO2_table[107].x = 1700., CO2_table[107].y = 7.134;
        CO2_table[108].x = 1800., CO2_table[108].y = 7.180;
        CO2_table[109].x = 1900., CO2_table[109].y = 7.222;

        CO2_table[110].x = 2000., CO2_table[110].y = 7.258;
        CO2_table[111].x = 2100., CO2_table[111].y = 7.291;
        CO2_table[112].x = 2200., CO2_table[112].y = 7.320;
        CO2_table[113].x = 2300., CO2_table[113].y = 7.347;
        CO2_table[114].x = 2400., CO2_table[114].y = 7.371;

        CO2_table[115].x = 2500., CO2_table[115].y = 7.393;
        CO2_table[116].x = 2600., CO2_table[116].y = 7.414;
        CO2_table[117].x = 2700., CO2_table[117].y = 7.433;
        CO2_table[118].x = 2800., CO2_table[118].y = 7.451;
        CO2_table[119].x = 2900., CO2_table[119].y = 7.468;

        CO2_table[120].x = 3000., CO2_table[120].y = 7.484;
        CO2_table[121].x = 3100., CO2_table[121].y = 7.499;
        CO2_table[122].x = 3200., CO2_table[122].y = 7.513;
        CO2_table[123].x = 3300., CO2_table[123].y = 7.526;
        CO2_table[124].x = 3400., CO2_table[124].y = 7.539;

        CO2_table[125].x = 3500., CO2_table[125].y = 7.551;
        CO2_table[126].x = 3600., CO2_table[126].y = 7.563;
        CO2_table[127].x = 3700., CO2_table[127].y = 7.575;
        CO2_table[128].x = 3800., CO2_table[128].y = 7.586;
        CO2_table[129].x = 3900., CO2_table[129].y = 7.597;

        CO2_table[130].x = 4000., CO2_table[130].y = 7.608;
        CO2_table[131].x = 4100., CO2_table[131].y = 7.618;
        CO2_table[132].x = 4200., CO2_table[132].y = 7.628;
        CO2_table[133].x = 4300., CO2_table[133].y = 7.638;
        CO2_table[134].x = 4400., CO2_table[134].y = 7.647;

        CO2_table[135].x = 4500., CO2_table[135].y = 7.657;
        CO2_table[136].x = 4600., CO2_table[136].y = 7.666;
        CO2_table[137].x = 4700., CO2_table[137].y = 7.676;
        CO2_table[138].x = 4800., CO2_table[138].y = 7.685;
        CO2_table[139].x = 4900., CO2_table[139].y = 7.694;

        CO2_table[140].x = 5000., CO2_table[140].y = 7.702;

        spline_pchip(ndat,CO2_table);
  
        initialized = TRUE;
      }

      if (temp < CO2_table[0].x || temp > CO2_table[ndat-1].x) {
        sprintf(Message,"temp=%g is outside data-table range [%g,%g]",temp,CO2_table[0].x,CO2_table[ndat-1].x);
        epic_error(dbmsname,Message);
      }
      else {
        i   = hunt_place_in_table(ndat,CO2_table,temp,&t_d,i);
        ans = splint_pchip(temp,CO2_table+i,t_d)*R_GAS/44.;
      }
    } break;
    default:
      sprintf(Message,"index=%d not yet implemented",index);
      epic_error(dbmsname,Message);
    break;
  }

  return ans;
}

/*======================= end of isobaric_specific_heat() ===================*/

/*======================= blackbody_fraction() ==============================*/

/*
 * For the given temperature [K], computes the fraction of blackbody emission
 * (the partial area under the Planck curve),  from 0 to the given wavelength [m].
 * Differencing this function for two different wavelengths gives the band-emission
 * fraction.
 *
 * The series used here, which converges rapidly, is from
 *
 *  Chang SL, Rhee KT, 1984, Blackbody radiation functions, Int. Comm. Heat Mass
 *      Transfer, 11, 451-455.
 *
 * Inputs: 
 *
 *    wavelength  [m]
 *    temperature [K]
 */

/* 
 * Second radiation constant, C2 = hc/k_b [m K]
 * FACTOR = 15.0/pi^4
 */
#define C2        1.438777e-2
#define FACTOR    0.15398973382026502783729
#define TOL       1.e-8
#define MAX_ITER  200

EPIC_FLOAT blackbody_fraction(EPIC_FLOAT wavelength,
                              EPIC_FLOAT temperature)
{
  int
    n;
  EPIC_FLOAT
    x,term,sum;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="blackbody_fraction";

  x = wavelength*temperature;
  if (x < 0.) {
    sprintf(Message,"wavelength*temperature=%g < 0, wavelength=%g, temperature=%g",
            x,wavelength,temperature);
    epic_error(dbmsname,Message);
  }
  else if (x == 0.) {
    return 0.;
  }
  x = C2/x;

  sum = 0.;
  for (n = 1; n <= MAX_ITER; n++) {
    term  = (6./(n*n*n)+x*(6./(n*n)+x*(3./n+x)))*exp(-(double)n*x)/n;
    sum  += term;
    if (term <= TOL) {
      return sum*FACTOR;
    }
  }

  /* Not enough terms allowed */
  sprintf(Message,"MAX_ITER=%d reached",MAX_ITER);
  epic_warning(dbmsname,Message);

  return sum*FACTOR;
}

#undef C2
#undef FACTOR
#undef TOL
#undef MAX_ITER

/*======================= end of blackbody_fraction() =======================*/

/* * * * * * * * * * * end of epic_funcs_diag.c  * * * * * * * * * * * * * * */



