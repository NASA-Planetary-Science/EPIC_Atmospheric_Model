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

/* * * * * * * epic_funcs_astron.c * * * * * * * * * * * * * * * * * * 
 *                                                                   *
 * Astronomical and geophysical functions that do not reference      *
 * EPIC model prognostic variables.                                  * 
 *                                                                   *
 *   solar_longitude()                                               *
 *   solar_declination()                                             *
 *   equation_of_time()                                              *
 *   east_longitude_of_noon()                                        *
 *   season_string()                                                 *
 *   radius_vector()                                                 *
 *   eccentric_anomaly()                                             *
 *   true_anomaly()                                                  *
 *   polar_radius()                                                  *
 *   equatorial_radius()                                             *
 *   solar_fraction()                                                *
 *                                                                   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*
 * Variables global to this file.
 */
double
  Kepler_M,
  Kepler_e;

/*
 * Function prototypes.
 */
double kepler_zero(double E);

/*====================== solar_longitude() =========================*/

/*
 * Converts calendar date/time [C data type time_t, sec] into L_s [deg], 
 * the planetocentric longitude of the Sun.
 */

double solar_longitude(planetspec *planet,
                       time_t      date)
{
  double
    tJ2000,
    M,Mrad,E,nu,
    alpha_fms,
    l_s,
    tmp;
  static time_t
    J2000;
  struct tm
    epoch2000;
  static int
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="solar_longitude";

  if (!initialized) {
    /* Trigger timezone evaluation */
    epoch2000.tm_year  =  0;
    epoch2000.tm_mon   =  0;
    epoch2000.tm_mday  =  0;
    epoch2000.tm_hour  =  0;
    epoch2000.tm_min   =  0;
    epoch2000.tm_sec   =  0;
    epoch2000.tm_isdst = -1;
    mktime(&epoch2000);

    /* Set the standard epoch J2000.0 in terms of UTC */
    epoch2000.tm_year = 2000-1900;
    epoch2000.tm_mon  = 0;
    epoch2000.tm_mday = 1;
    epoch2000.tm_hour = 11-timezone/3600;
    epoch2000.tm_min  = 58;
    epoch2000.tm_sec  = 56;

    J2000 = mktime(&epoch2000);

    initialized = TRUE;
  }

  /*
   * Elapsed time since the J2000 epoch [day].
   */
  tJ2000 = difftime(date,J2000)/86400.;

  switch(planet->index) {
    case MARS_INDEX:
      /*
       * Allison M, 1997, Accurate analytical representations of solar time and seasons on Mars
       * with applications to the Pathfinder/Surveyor missions, Geophys. Res. Lett. 24, 1967-1970
       */
      Mrad      = (19.41+0.5240212*tJ2000)*DEG;
      alpha_fms = 270.39+0.5240384*tJ2000;
      l_s       = alpha_fms+(10.691+3.7e-7*tJ2000)*sin(   Mrad)
                                            +0.623*sin(2.*Mrad)
                                            +0.050*sin(3.*Mrad)
                                            +0.005*sin(4.*Mrad);
    break;
    default:
      /*
       * Use Keplerian elements to determine L_s.
       */

      /* Calculate the mean anomaly [deg] */
      M = (planet->mean_lon-planet->lon_perihelion)+(360./planet->orbit_period)*(tJ2000/365.25);

      /* Solve Kepler's equation for the eccentric anomaly, E [deg] */
      E = eccentric_anomaly(M,planet->e);

      /* Solve for the true anomaly, nu [deg]. */
      nu = true_anomaly(E,planet->e);

      /* Calculate L_s */
      l_s = nu-planet->vernal_equinox_anomaly;
    break;
  }

  /* Map l_s to [0.,360.]. */
  l_s = 360.*modf(1.+modf(l_s/360.,&tmp),&tmp);

  return l_s;
}


/*====================== end of solar_longitude() ==================*/

/*====================== solar_declination() =======================*/

/*
 * Returns the declination of the Sun [deg], also called the
 * sub-solar latitude.
 *
 * gurn: need to decide if IAU vs dynamic north flips this
 *
 * See Meeus (2005), (13.4), with the ecliptic latitude of the Sun,
 * beta = 0., and the solar longitude, l_s, which yields
 *  sin(dec_s) = sin(obliquity)*sin(l_s)
 */


double solar_declination(planetspec *planet,
                         double      l_s)
{
  double
    dec_s,
    tmp;

  dec_s = asin(sin(planet->obliquity*DEG)*sin(l_s*DEG))/DEG;

  return dec_s;
}

/*====================== end of solar_declination() ================*/

/*====================== equation_of_time() ========================*/

/*
 * Return the difference between the planet's apparent solar time and its
 * mean solar time as expressed in west longitude [deg] (positive to
 * the right on the sky as seen from the planet), given a
 * calender date and the solar longitude [deg], l_s in the range [0.,360.].
 *
 * If l_s is negative, calculate it internally given the
 * planet and date.  We take l_s as an argument to save computational
 * time if it has already been determined.
 */

double equation_of_time(planetspec *planet,
                        time_t      date,
                        double      l_s)
{
  double
    ans,alpha_s,
    M,l_s_mean,
    tJ2000,
    tmp;
  static time_t
    J2000;
  struct tm
    epoch2000;
  static int
    initialized = FALSE;

  if (!initialized) {
    /* Trigger timezone evaluation. */
    epoch2000.tm_year  =  0;
    epoch2000.tm_mon   =  0;
    epoch2000.tm_mday  =  0;
    epoch2000.tm_hour  =  0;
    epoch2000.tm_min   =  0;
    epoch2000.tm_sec   =  0;
    epoch2000.tm_isdst = -1;
    mktime(&epoch2000);

    /* Set the standard epoch J2000.0 in terms of UTC. */
    epoch2000.tm_year = 2000-1900;
    epoch2000.tm_mon  = 0;
    epoch2000.tm_mday = 1;
    epoch2000.tm_hour = 11-timezone/3600;
    epoch2000.tm_min  = 58;
    epoch2000.tm_sec  = 56;

    J2000 = mktime(&epoch2000);

    initialized = TRUE;
  }

  if (l_s < 0.) {
    /*
     * User has asked for l_s to be calculated here.
     */
    l_s = solar_longitude(planet,date);
  }
  else {
    /* 
     * Make sure l_s is mapped to [0.,360.].
     */
    l_s = 360.*modf(1.+modf(l_s/360.,&tmp),&tmp);
  }

  /*
   * Elapsed time since the J2000 epoch [day].
   */
  tJ2000 = difftime(date,J2000)/86400.;

  /* 
   * Calculate the mean anomaly [deg]
   */
  M = (planet->mean_lon-planet->lon_perihelion)+(360./planet->orbit_period)*(tJ2000/365.25);

  /* Map M to ecliptic longitude of mean Sun [0.,360.]. */
  l_s_mean = M-planet->vernal_equinox_anomaly;
  /* map to [0.,360.] */
  l_s_mean = 360.*modf(1.+modf(l_s_mean/360.,&tmp),&tmp);

  /*
   * The solar longitude, l_s, and the right ascention of the Sun, alpha_s,
   * are the same concept in the ecliptic and equatorial planes, respectively.
   * See Meeus (2005, Astronomical Algorithms), Chap. 13, for the
   * transformation between ecliptic and equatorial coordinates.
   * The Sun's ecliptic latitude, beta_s, is zero, hence Meeus (13.3) yields
   *   tan(alpha_s) = cos(obliquity)*tan(l_s)
   */

  /*
   * NOTE: Care must be taken to use the correct branch of the atan function
   *       to ensure a continuous result.  The Wikipedia page "Equation of time"
   *       gives a useful discussion of this issue.
   */
  if (l_s < 90.) {
    alpha_s = atan(cos(planet->obliquity*DEG)*tan(l_s*DEG))/DEG;
  }
  else if (l_s == 90.) {
    alpha_s = l_s;
  }
  else if (l_s < 270.) {
    alpha_s = atan(cos(planet->obliquity*DEG)*tan(l_s*DEG))/DEG+180.;
  }
  else if (l_s == 270.) {
    alpha_s = l_s;
  }
  else {
    alpha_s = atan(cos(planet->obliquity*DEG)*tan(l_s*DEG))/DEG+360.;
  }

  ans = l_s_mean-alpha_s;

  /* Remove aliasing. */
  while (ans < -180.) {
    ans += 360.;
  }
  while (ans > 180.) {
    ans -= 360.;
  }

  return ans;
}

/*====================== end of equation_of_time() =================*/

/*====================== east_longitude_of_solar_noon() ============*/

/*
 * Returns the eastward-increasing longitude [-180, 180] used in 
 * the EPIC model, where it is solar noon (where the Sun is on the meridian)
 * for the specified planet, as a function of calendar time.
 * The zero of east longitude is aligned with the zero of west longitude [0, 360].
 *
 * If l_s is negative, calculate it internally given the
 * planet and date.  We take l_s as an argument to save computational
 * time if it has already been determined.
 */

double east_longitude_of_solar_noon(planetspec *planet,
                                    time_t      date,
                                    double      l_s)
{
  double
    west_lon,east_lon,
    tJ2000,
    tmp;
  static time_t
    J2000;
  struct tm
    epoch2000;
  static int
    initialized = FALSE;

  if (!initialized) {
    /* Trigger timezone evaluation */
    epoch2000.tm_year  =  0;
    epoch2000.tm_mon   =  0;
    epoch2000.tm_mday  =  0;
    epoch2000.tm_hour  =  0;
    epoch2000.tm_min   =  0;
    epoch2000.tm_sec   =  0;
    epoch2000.tm_isdst = -1;
    mktime(&epoch2000);

    /* Set the standard epoch J2000.0 in terms of UTC */
    epoch2000.tm_year = 2000-1900;
    epoch2000.tm_mon  = 0;
    epoch2000.tm_mday = 1;
    epoch2000.tm_hour = 11-timezone/3600;
    epoch2000.tm_min  = 58;
    epoch2000.tm_sec  = 56;

    J2000 = mktime(&epoch2000);

    initialized = TRUE;
  }

  /*
   * Elapsed time since the J2000 epoch [s].
   */
  tJ2000 = difftime(date,J2000);

  if (l_s < 0.) {
    /*
     * User has asked for l_s to be calculated here.
     */
    l_s = solar_longitude(planet,date);
  }
  else {
    /* 
     * Make sure l_s is mapped to [0.,360.].
     */
    l_s = 360.*modf(1.+modf(l_s/360.,&tmp),&tmp);
  }

  /*
   * Determine west longitude of clock (mean) noon.
   */
  west_lon = planet->wlon_noon_J2000+tJ2000*planet->omega_synodic/DEG;

  /*
   * Use the equation of time to shift clock noon to solar noon.
   */
  west_lon += equation_of_time(planet,date,l_s);

  /* Map west longitude to [0.,360.] */
  west_lon = 360.*modf(1.+modf(west_lon/360.,&tmp),&tmp);

  /* 
   * Calculate the corresponding east longitude [-180.,180.],
   * which is the longitude system used by EPIC.
   */
  if (west_lon < 180.) {
    east_lon = -west_lon;
  }
  else {
    east_lon = 360.-west_lon;
  }

  return east_lon;
}

/*====================== end of east_longitude_of_solar_noon() =====*/

/*====================== season_string() ===========================*/

void season_string(double  l_s, 
                   char   *outstring)
{

  if (fcmp(l_s,0.) == 0 || fcmp(l_s,360.) == 0) {
    sprintf(outstring,"vernal equinox");
  }
  else if (fcmp(l_s,90.) < 0) {
    sprintf(outstring,"northern spring");
  }
  else if (fcmp(l_s,90.) == 0) {
    sprintf(outstring,"summer solstice");
  }
  else if (fcmp(l_s,180.) < 0) {
    sprintf(outstring,"northern summer");
  }
  else if (fcmp(l_s,180.) == 0) {
    sprintf(outstring,"autumnal equinox");
  }
  else if (fcmp(l_s,270.) < 0) {
    sprintf(outstring,"northern autumn");
  }
  else if (fcmp(l_s,270.) == 0) {
    sprintf(outstring,"winter solstice");
  }
  else {
    sprintf(outstring,"northern winter");
  }

  return;
}

/*====================== end of season_string() ====================*/

/*====================== radius_vector() ===========================*/

/*
 *  Input: calendar time; see comments for solar_longitude()
 * Output: distance between planet and the Sun [AU]
 *
 * See Chap. 30 Equation of Kepler, p. 195, in 
 *   Meeus J, 2005, Astronomical Algorithms, Willmann-Bell
 */

double radius_vector(planetspec *planet,
                     time_t      date)
{
  double
    tJ2000,
    M,E;
  static time_t
    J2000;
  struct tm
    epoch2000;
  static int
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="radius_vector";

  if (!initialized) {
    /* Trigger timezone evaluation. */
    epoch2000.tm_year  =  0;
    epoch2000.tm_mon   =  0;
    epoch2000.tm_mday  =  0;
    epoch2000.tm_hour  =  0;
    epoch2000.tm_min   =  0;
    epoch2000.tm_sec   =  0;
    epoch2000.tm_isdst = -1;
    mktime(&epoch2000);

    /* Set the standard epoch J2000.0 in terms of UTC. */
    epoch2000.tm_year = 2000-1900;
    epoch2000.tm_mon  = 0;
    epoch2000.tm_mday = 1;
    epoch2000.tm_hour = 11-timezone/3600;
    epoch2000.tm_min  = 58;
    epoch2000.tm_sec  = 56;

    J2000 = mktime(&epoch2000);

    initialized = TRUE;
  }

  /*
   * Elapsed time since the J2000 epoch [day].
   */
  tJ2000 = difftime(date,J2000)/86400.;

  /* 
   * Calculate the mean anomaly [deg]
   */
  M = (planet->mean_lon-planet->lon_perihelion)+(360./planet->orbit_period)*(tJ2000/365.25);

  /* Solve Kepler's equation for the eccentric anomaly, E [deg] */
  E = eccentric_anomaly(M,planet->e);

  return  planet->a*(1.-planet->e*cos(E*DEG));
}
                     
/*====================== end of _radius_vector() ===================*/

/*====================== kepler_zero() =============================*/

double kepler_zero(double E)
{
  return (Kepler_M-E)*DEG+Kepler_e*sin(E*DEG);
}

/*====================== end of kepler_zero() ======================*/

/*====================== eccentric_anomaly() =======================*/

/*
 * Solve Kepler's equation, M = E - e sin E, for eccentric anomaly, E,
 * given the mean anomaly, M, and the orbital eccentricity, e.
 */

double eccentric_anomaly(double M,
                         double e)
{
  double
    E,
    tol = 1.e-6,
    tmp;
  int
    error_flag;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="eccentric_anomaly";

  /* Map M to [0.,360.]. */
  M = 360.*modf(1.+modf(M/360.,&tmp),&tmp);

  /*
   * Handle special cases.
   */
  if (fcmp(M,  0.) == 0) return   0.;
  if (fcmp(M,180.) == 0) return 180.;
  if (fcmp(M,360.) == 0) return 360.;

  /* Set global variable to communicate with function kepler_zero(). */
  Kepler_M = M;
  Kepler_e = e;

  error_flag = find_root(M-180.,M+180.,tol,&E,kepler_zero);
  if (error_flag) {
    if (error_flag) {
      sprintf(Message,"Error solving Kepler's equation, find_root(): error_flag=%d",error_flag);
      epic_error(dbmsname,Message);
    }
  }

  /* Map E to [0.,360.]. */
  E = 360.*modf(1.+modf(E/360.,&tmp),&tmp);

  return E;
}

/*====================== end of eccentric_anomaly() ================*/

/*====================== true_anomaly() ============================*/

/*
 * Returns true anomaly, nu [deg], given eccentric anomaly, E [deg],
 * and orbital eccentricity, e.
 */

double true_anomaly(double E,
                    double e)
{
  double
    Erad_2,
    nu,
    tmp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="true_anomaly";

  /* Map E to [0.,360.]. */
  E = 360.*modf(1.+modf(E/360.,&tmp),&tmp);

  /*
   * Handle special cases.
   */
  if (fcmp(E,  0.) == 0) return   0.;
  if (fcmp(E,180.) == 0) return 180.;
  if (fcmp(E,360.) == 0) return 360.;

  if (fcmp(e,1.) < 0) {
    if (E < 180.) {
      nu = 2.*atan(sqrt((1.+e)/(1.-e))*tan(E*.5*DEG))/DEG;
    }
    else {
      nu = 2.*atan(sqrt((1.+e)/(1.-e))*tan((E-360.)*.5*DEG))/DEG+360.;
    }
  }
  else {
    sprintf(Message,"e=%g >= 1.",e);
    epic_error(dbmsname,Message);
  }

  return nu;
}

/*====================== end of true_anomaly() =====================*/

/*======================= polar_radius() ===========================*/

/*
 * Returns the polar radius, c, given the equatorial radius, a, and J2, GM, and omega.
 * Assumes a uniformly rotating fluid planet, and uses a root finder to solve the
 * cubic equation for ep = c/a.
 */

EPIC_FLOAT polar_radius(EPIC_FLOAT a,
                        EPIC_FLOAT J2,
                        EPIC_FLOAT GM,
                        EPIC_FLOAT omega)
{
  int
    error_flag;
  double
    a3w2_GM,
    ep,ep0,ep1,ep2,
    eptol;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="polar_radius";

  /* Error tolerance */
  eptol = pow(machine_epsilon(),2./3.);

  /* Nondimensional rotational parameter */
  a3w2_GM = a*omega;
  a3w2_GM = a3w2_GM*a3w2_GM*(a/GM);

  /* 
   * Initial guess.
   * See for example, Section 5-4 on the geopotential in
   * Turcotte and Schubert's 1982 text "Geodynamics".
   */
  ep0 = 1.-1.5*(double)J2-.5*(double)a3w2_GM;
  ep1 = MAX(0.,2.*ep0-1.);
  ep2 = 1.;

  PRC_J2      = J2;
  PRC_a3w2_GM = a3w2_GM;

  error_flag = find_root(ep1,ep2,eptol,&ep,polar_radius_cubic);

  if (error_flag) {
    sprintf(Message,"a,J2,GM,omega=%g %g %g %g; find_root() error_flag=%d",a,J2,GM,omega,error_flag);
    epic_error(dbmsname,Message);
  }

  return (EPIC_FLOAT)(a*ep);
}

/*======================= end of polar_radius() =============================*/

/*======================= polar_radius_cubic() ==============================*/

/* 
 * Cubic equation involving ratio of polar to equatorial radius, ep = c/a, 
 * used with polar_radius().
 * See for example (5-57) of Turcotte and Schubert's 1982 text "Geodynamics".
 */

EPIC_FLOAT polar_radius_cubic(EPIC_FLOAT ep)
{
  return ep*(ep*(ep*(1.+.5*PRC_J2+.5*PRC_a3w2_GM)-1.))+PRC_J2;
}

/*======================= end of polar_radius_cubic() =======================*/

/*======================= equatorial_radius() ===============================*/

/*
 * Returns the equatorial radius, a, given the polar radius, c, and J2, GM, and omega.
 * Assumes a uniformly rotating fluid planet, and uses a root finder to solve the
 * cubic equation for ep = c/a.
 */

EPIC_FLOAT equatorial_radius(EPIC_FLOAT c,
                             EPIC_FLOAT J2,
                             EPIC_FLOAT GM,
                             EPIC_FLOAT omega)
{
  int
    error_flag;
  double
    c3w2_GM,
    ep,ep0,ep1,ep2,
    eptol;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="equatorial_radius";

  /* Error tolerance */
  eptol = pow(machine_epsilon(),2./3.);

  /* Nondimensional rotational parameter */
  c3w2_GM = c*omega;
  c3w2_GM = c3w2_GM*c3w2_GM*(c/GM);

  /* 
   * Initial guess.
   * See for example, Section 5-4 on the geopotential in
   * Turcotte and Schubert's 1982 text "Geodynamics".
   */
  ep0 = 1.-1.5*(double)J2-.5*(double)c3w2_GM;
  ep1 = MAX(0.,2.*ep0-1.);
  ep2 = 1.;

  ERC_J2      = J2;
  ERC_c3w2_GM = c3w2_GM;

  error_flag = find_root(ep1,ep2,eptol,&ep,equatorial_radius_cubic);

  if (error_flag) {
    sprintf(Message,"c,J2,GM,omega=%g %g %g %g; find_root() error_flag=%d",c,J2,GM,omega,error_flag);
    epic_error(dbmsname,Message);
  }

  return (EPIC_FLOAT)(c/ep);
}

/*======================= end of equatorial_radius() ========================*/

/*======================= equatorial_radius_cubic() =========================*/

/* 
 * Cubic equation involving ratio of polar to equatorial radius, ep = c/a, 
 * used with equatorial_radius().
 * See for example (5-57) of Turcotte and Schubert's 1982 text "Geodynamics".
 */

EPIC_FLOAT equatorial_radius_cubic(EPIC_FLOAT ep)
{
  return ep*(ep*(ep*(1.+.5*ERC_J2)-1.))+ERC_J2+.5*ERC_c3w2_GM;
}


/*======================= end of equatorial_radius_cubic() ==================*/

/*====================== solar_fraction() ===================================*/

/*
 * Returns solar mixing ratio of the least-abundant element in
 * the given species name (divided by its stochiometric count).
 * Choices for the type argument: MASS, MOLAR. 
 * The character string min_element should be 4 bytes.
 *
 * NOTE: Here, "mixing ratio" has the total in the denominator, and hence
 *       in the language of meteorology is analogous to "specific humidity", rather
 *       than to a "mixing ratio" that has the dry component in the denominator.
 */

EPIC_FLOAT solar_fraction(char *species,
                          int   type,
                          char *min_element)
{
  int
    i,ii,ii_min,
    min_count;
  EPIC_FLOAT
    ratio,
    mu,
    min_abundance,
    abundance,tmp;
  static int
    initialized  = FALSE,
    num_elements = 0,
    *counts      = NULL;
  static EPIC_FLOAT
    total_number,
    total_mass,
    n_H_2;
  static char
    **symbols    = NULL;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="solar_fraction";

  if (!initialized) {
    /* 
     * Add up total number and total mass:
     */
    i = 1;
    abundance = pow(10.,Element[i].solar_abundance);
    /* 
     * Assume hydrogen is molecular (H_2) rather than atomic (H)
     * and that C, N, O, and S are in their reduced forms.
     */
    n_H_2 = .5*(pow(10.,Element[ 1].solar_abundance)
                    -4.*pow(10.,Element[ 6].solar_abundance)
                    -3.*pow(10.,Element[ 7].solar_abundance)
                    -2.*pow(10.,Element[ 8].solar_abundance)
                    -2.*pow(10.,Element[16].solar_abundance));

    total_number = n_H_2;

    total_mass = abundance*(Element[1].molar_mass);

    for (i = 2; i <= LAST_ATOMIC_NUMBER; i++) {
      abundance     = pow(10.,Element[i].solar_abundance);
      total_number += abundance;
      total_mass   += abundance*(Element[i].molar_mass);
    }

    initialized = TRUE;
  }

  /* Check for null string: */
  if (species == NULL || *species == '\0') {
    return 0.;
  }

  parse_species_name(species,&num_elements,&symbols,&counts);

  /*
   * Return ratio = 0. if num_elements is zero.
   */
  if (num_elements == 0) {
    ratio = 0.;
    return ratio;
  }

  /*
   * Sum molar masses of components to get total, mu, and 
   * find abundance of least-abundant element in species:
   */

  mu            =  0.;
  min_abundance =  Element[1].solar_abundance;
  ii_min        = -1;
  for (ii = 0; ii < num_elements; ii++) {
    /* Identify element */
    for (i = 1; i <= LAST_ATOMIC_NUMBER; i++) {
      if (strcmp(Element[i].symbol,symbols[ii]) == 0) {
        mu  += Element[i].molar_mass*counts[ii];
        tmp  = Element[i].solar_abundance;
        if (tmp <= min_abundance) {
          min_abundance = tmp;
          ii_min        = ii;
        }
        break;
      }
    }
  }

#if EPIC_CHECK == 1
  /* Sanity check on ii_min: */
  if (ii_min < 0) {
    epic_error(dbmsname,"ii_min < 0");
  }
#endif

  min_count = counts[ii_min];
  strcpy(min_element,symbols[ii_min]);

  abundance = pow(10.,min_abundance)/min_count;

  if (type == MOLAR) {
    if (strcmp(min_element,"H") == 0) {
      ratio = n_H_2/total_number;
    }
    else {
      ratio = abundance/total_number;
    }
  }
  else if (type == MASS) {
    if (strcmp(min_element,"H") == 0) {
      ratio = n_H_2*2.*Element[1].molar_mass;
    }
    else {
      ratio = abundance*mu/total_mass;
    }
  }
  else {
    sprintf(Message,"Unknown type %d",type);
    epic_error(dbmsname,Message);
  }

  return ratio;
}

/*====================== end of solar_fraction() ===================*/

/* * * * * * * * * * * * end of epic_funcs_astron.c * * * * * * * * */
