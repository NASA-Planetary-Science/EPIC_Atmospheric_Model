/**************************************************************
 *
 *  us_standard_atmos.c:  A program to produce standard atmosphere
 *  data (temperature, pressure and density) for Earth, for any number of
 *  points desired. For example, to give smoother initial data input
 *  for EPIC, we may tell this program how large to make the distance
 *  in height for these quantities, then have as many initial data
 *  points as desired.  This is as opposed to using table input,
 *  with predetermined gap in the height of the points).
 *
 *  Two output files are produced, sounding.dat and t_vs_p.us_standard.
 *  Each is a format used by EPIC for sounding data.
 *
 *  Example compilation:  clang -lm -o us_standard_atmos us_standard_atmos.c
 *  For debugging compile with:  clang -lm -g -o us_standard_atmos us_standard_atmos.c
 *
 *  ME Bradley, Summer 2010
 *
 ***************************************************************/
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define T0    288.15      /* Surface value in Kelvin */
#define P0   1013.25      /* Surface pressure in hPa */
#define RHO0    1.225     /* Surface air density in kg/m^3 */
#define g0      9.80665   /* standard Earth gravity */

int main() {
  int 
    nk,
    num_levels;
  double 
    h,
    dh,
    dhfeet,
    *T, *p, *rho,*H;
  static int initialized = 0;
  
  FILE *sounding, *t_vs_p;
  
  sounding = fopen("sounding.dat","w");
  t_vs_p = fopen("t_vs_p.us_standard","w");
  
  fprintf(stdout,"Input distance between each level dh[m]:  ");
  fscanf(stdin,"%lf",&dh);
  fprintf(stdout,"%.2f meters\n",dh);
  
  /* Convert meters to feet in order to use formula */
  dhfeet = dh*3.2808399;
  fprintf(stdout,"%.2f feet\n",dhfeet);
  
  /*  Number of levels needed to range -4003 ft to 278,386 ft 
      (i.e. a total of 282389 ft) */
  num_levels = (int)282389/dhfeet;
  fprintf(stdout,"number of levels %i \n",num_levels);
 
  /*  Write headers for the files sounding.dat at t_vs_p.us_standard */
  fprintf(sounding,"U.S. Standard Atmosphere (1976),\n");
  fprintf(sounding,"calculated using epic/tools/C/us_standard_atmos.c\n");
  fprintf(sounding,"based on formulas from http://www.atmosculator.com/The Standard Atmosphere.html\n\n");
  fprintf(sounding,"# Geopot. Ht[km]  Temp[K]  Press[hPa]  Dens[kg/m^3]\n");
  fprintf(sounding,"%i\n",num_levels+1);
  
  fprintf(t_vs_p,"U.S. Standard Atmosphere (1976)\n");
  fprintf(t_vs_p,"calculated using epic/tools/C/us_standard_atmos.c\n");
  fprintf(t_vs_p,"based on formulas from http://www.atmosculator.com/The Standard Atmosphere.html\n.\n.\n");
  fprintf(t_vs_p,"#  p[hPa]     T[K]    dT[K]\n");
  fprintf(t_vs_p,"%i\n",num_levels+1);
 
  if (!initialized){
    H   = (double *)calloc(num_levels+1,sizeof(double));
    T   = (double *)calloc(num_levels+1,sizeof(double));
    p   = (double *)calloc(num_levels+1,sizeof(double));
    rho = (double *)calloc(num_levels+1,sizeof(double));

    initialized = 1;
  }

  for(nk = 0; nk <= num_levels; nk++){
    /*  Generate h in feet, where the bottom is a height of -4003ft = -1220m */
    h = (double)(nk)*dhfeet -4003;
 
    /*  Store height in km, not feet! */
    H[nk] = (h*.3048)/1000.;
    
    /* 
     * Formulas for T, p, and rho from 0 up to 278,386 ft
     * derived from website:  www.atmosculator.com/The Standard Atmosphere.html
     */

    if (h <= 36089.){
      T[nk] = T0*(1. - h/145442.);
      p[nk] = P0*pow((1. - h/145442.),5.255876);
      rho[nk] = RHO0*pow((1.-h/145442.),4.255876);
    }
    /*  Isothermal level 1*/
    else if (h <= 65617.){
      T[nk] = T0*.751865;
      p[nk] = P0*.223361*exp(-(h-36089.)/20806.);
      rho[nk] = RHO0*.297076*exp(-(h-36089.)/20806.);
    }
    /* Inversion level 1 */
    else if (h <= 104987.){     
      T[nk] = T0*(.682457 + h/945374.);
      p[nk] = P0*pow((.988626 + h/652600.),-34.16320);
      rho[nk] = RHO0*pow((.978261 + h/659515.),-35.16320);
    }
    /* Inversion level 2 */
    else if (h <= 154199.){
      T[nk] = T0*(.482561 + h/337634.);
      p[nk] = P0*pow((.898309 + h/181373.),-12.20114);
      rho[nk] = RHO0*pow((.857003 + h/190115.),-13.20114);
    }
    /* Isothermal level 2 */
    else if (h <= 167323.){
      T[nk] = T0*.939268;
      p[nk] = P0*.00109456*exp(-(h-154199.)/25992.);
      rho[nk] = RHO0*.00116533*exp(-(h-154199.)/25992.);
    }
    else if (h <= 232940.){
      T[nk] = T0*(1.434843 - h/337634.);
      p[nk] = P0*pow((.838263 - h/577922.),12.20114);
      rho[nk] = RHO0*pow((.798990 -h/606330.),11.20114);
    }
    else if (h <= 278386.){
      T[nk] = T0*(1.237723 - h/472687.);
      p[nk] = P0*pow((.917131 - h/637919),17.08160);
      rho[nk] = RHO0*pow((.900194 -h/649922.),16.08160);
    }
    else {
      fprintf(stdout,"ERROR:  Elevation is Off The Chart!");
      return 1;
    }
  }

  for (nk = 0; nk <= num_levels; nk++){
    fprintf(sounding,"    %6.3f       %7.2f   %9.3e     %4.3e\n",H[nk],T[nk],p[nk],rho[nk]);
  }
  for (nk = num_levels; nk >= 0; nk--){
    fprintf(t_vs_p," %9.3e    %5.2f   %2.1f \n",p[nk],T[nk],0.0);
  }

  fclose(sounding);
  fclose(t_vs_p);
} 
  
