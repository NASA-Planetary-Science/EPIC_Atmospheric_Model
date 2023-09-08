/*
 * Compile with: clang -lm -o extend_sounding extend_sounding.c
 */

#include <stdio.h>
#include <math.h>
      
int main() {
  double
   T[]   = {    93.5027,    92.9007},  /* K      */
   z[]   = {     3.,       506.0   },  /* m      */
   p[]   = {146645.,    143150.0   },  /* Pa     */
   rho[] = {     5.3446,     5.2432};  /* kg/m^3 */
  double
   gam,x,R,
   temp,dens,press,zee;

  gam = (T[0]-T[1])/(z[1]-z[0]);

  x = log(p[0]/p[1])/log(T[0]/T[1]);

  R = p[0]/(rho[0]*T[0]);

  while (zee < 1000.) {
    printf("Input z[m]\n");
    scanf("%lf",&zee);

    temp  = T[0]-gam*(zee-z[0]);
    press = p[0]*pow(temp/T[0],x);
    dens  = press/(R*temp);

    fprintf(stdout,"z=%.1f p=%.1f T=%.4f rho=%.7f\n",zee,press,temp,dens);
  }

  return 0;
}


