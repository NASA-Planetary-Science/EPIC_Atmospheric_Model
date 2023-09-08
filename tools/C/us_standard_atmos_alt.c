/**************************************************************
 *
 *  us_standard_atmos_alt.c:  
 *  A program to produce standard atmosphere
 *  profile for Earth, including potential temperature, Theta, and
 *  the NH profile (product of buoyancy frequency
 *  and pressure scale height).
 *
 *  Example compilation:         clang -lm -o us_standard_atmos_alt us_standard_atmos_alt.c
 *  For debugging compile with:  clang -lm -g -o us_standard_atmos_alt us_standard_atmos_alt.c
 *
 *  TE Dowling, Summer 2013: Modified to include Theta and NH
 *  ME Bradley, Summer 2010: us_standard_atmos.c
 *
 ***************************************************************/
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define T0     288.15      /* Surface value in Kelvin       */
#define P0    1013.25      /* Surface pressure in hPa       */
#define Rho0     1.225     /* Surface air density in kg/m^3 */
#define g0       9.80665   /* standard Earth gravity        */
#define kappa    0.286     /* Rdry/c_p for Earth            */
#define Rdry   287.047     /* Rdry for Earth                */

int main() {
  int 
    k,
    num_levels;
  double 
    h,dh,dhfeet,
    dlntheta,dz;
  double
   *Z,
   *T,
   *p,
   *rho,
   *theta,
   *N2,
   *NH;
  FILE
   *sounding;
  
  fprintf(stdout,"Input distance between each level dh [m]:  ");
  fscanf(stdin,"%lf",&dh);
  fprintf(stdout,"%.2f meters\n",dh);
  
  /* Convert meters to feet in order to use formula */
  dhfeet = dh*3.2808399;
  
  /*
   * Number of levels needed to range 0 ft to 278,386 ft
   */
  num_levels = (int)278386/dhfeet;
  fprintf(stdout,"number of levels %i \n",num_levels);
    
  /* Allocate memory */
  Z     = (double *)calloc(num_levels+1,sizeof(double));
  T     = (double *)calloc(num_levels+1,sizeof(double));
  p     = (double *)calloc(num_levels+1,sizeof(double));
  rho   = (double *)calloc(num_levels+1,sizeof(double));
  theta = (double *)calloc(num_levels+1,sizeof(double));
  N2    = (double *)calloc(num_levels+1,sizeof(double));
  NH    = (double *)calloc(num_levels+1,sizeof(double));

  for (k = 0; k <= num_levels; k++){
    /*  Generate h in feet, where the bottom is a height of -4003ft = -1220m */
    h = (double)(k)*dhfeet;
 
    /*  Store height in km, not feet! */
    Z[k] = (h*.3048)/1000.;
    
    /* 
     * Formulas for T, p, and rho from -1220m = -4003ft up to 278,386 ft
     * derived from website:  www.atmosculator.com/The Standard Atmosphere.html
     */

    /* Up to end of "ground" level; includes below sea level values of h */
    if (h <= 36089.){
      T[k]   = T0*(1.-h/145442.);
      p[k]   = P0*pow((1.-h/145442.),5.255876);
      rho[k] = Rho0*pow((1.-h/145442.),4.255876);
    }
    /*  Isothermal level 1*/
    else if (h <= 65617.){
      T[k]   = T0*.751865;
      p[k]   = P0*.223361*exp(-(h-36089.)/20806.);
      rho[k] = Rho0*.297076*exp(-(h-36089.)/20806.);
    }
    /* Inversion level 1 */
    else if (h <= 104987.){     
      T[k]   = T0*(.682457 + h/945374.);
      p[k]   = P0*pow((.988626 + h/652600.),-34.16320);
      rho[k] = Rho0*pow((.978261 + h/659515.),-35.16320);
    }
    /* Inversion level 2 */
    else if (h <= 154199.){
      T[k]   = T0*(.482561 + h/337634.);
      p[k]   = P0*pow((.898309 + h/181373.),-12.20114);
      rho[k] = Rho0*pow((.857003 + h/190115.),-13.20114);
    }
    /* Isothermal level 2 */
    else if (h <= 167323.){
      T[k]   = T0*.939268;
      p[k]   = P0*.00109456*exp(-(h-154199.)/25992.);
      rho[k] = Rho0*.00116533*exp(-(h-154199.)/25992.);
    }
    else if (h <= 232940.){
      T[k]   = T0*(1.434843-h/337634.);
      p[k]   = P0*pow((.838263 - h/577922.),12.20114);
      rho[k] = Rho0*pow((.798990 -h/606330.),11.20114);
    }
    else if (h <= 278386.){
      T[k]   = T0*(1.237723-h/472687.);
      p[k]   = P0*pow((.917131-h/637919),17.08160);
      rho[k] = Rho0*pow((.900194-h/649922.),16.08160);
    }
    else {
      fprintf(stdout,"ERROR: Elevation is Off The Chart!");
      return 1;
    }

    theta[k] = T[k]*pow(1000./p[k],kappa);
  }

  /*
   * Calculate square of buoyancy (Brunt-Vaisaila) frequency, N2 = g(d/dz)ln(theta).
   * Use centered differences for the interior points.
   */
  for (k = 1; k < num_levels; k++) {
    dlntheta = log(theta[k+1]/theta[k-1]);
    dz       = (Z[k+1]-Z[k-1])*1000.;       /* converting from [km] to [m] */
    N2[k]    = g0*dlntheta/dz;
  }
  /* Use linear extrapolation to fill in the top and bottom points. */
  k = 0;
  N2[k] = 2.*N2[k+1]-N2[k+2];
  k = num_levels;
  N2[k] = 2.*N2[k-1]-N2[k-2];

  /* 
   * Calculate NH (product of buoyancy frequency and pressure scale height).
   */
  for (k = 0; k <= num_levels; k++) {
    NH[k] = sqrt(N2[k])*Rdry*T[k]/g0;
  }

  /* Write us_standard_atmos.dat */
  sounding = fopen("./us_standard_atmos.dat","w");

  fprintf(sounding,"U.S. Standard Atmosphere (1976), calculated using epic/tools/C/us_standard_atmos_alt.c\n");
  fprintf(sounding,"Based on formulas from http://www.atmosculator.com/The Standard Atmosphere.html\n\n");
  fprintf(sounding,"  Z[km]  Temp[K]  Press[hPa]  Dens[kg/m^3]  Theta[K]  N2[1/s^2]  NH[m/s]\n");
  for (k = 0; k <= num_levels; k++){
    fprintf(sounding," %6.3f  %7.2f   %9.3e     %4.3e   %7.2f  %9.3e %8.2f\n",
                     Z[k],T[k],p[k],rho[k],theta[k],N2[k],NH[k]);
  }

  fclose(sounding);

  /* Free allocated memory */
  free(Z),free(T),free(p),free(rho);
  free(theta),free(N2),free(NH);

  return 0;
} 
  
