/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2023 Timothy E. Dowling                      *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License     *
 * as published by the Free Software Foundation; either version 2  *
 * of the License, or (at your option) any later version.          *
 * A copy of this License is in the file:                          *
 *   $EPIC_PAT/License.txt                                         *
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

/* * * * * * * * * * * * * * * epic_globals.c  * * * * * * * * * * *
 *                                                                 * 
 *  Global declarations and related parameter data go here.        *
 *  These are referenced elsewhere, such as in epic.h, with        *
 *  the "extern" modifier.                                         *
 *                                                                 *
 *  Data types are defined in epic_datatypes.h, including          *
 *  structures defined using typedef.                              *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*
 * Globally defined structures.
 */
char
  Message[N_STR];
gridspec 
  grid;
thermospec
  thermo;
variablespec
  var;

/*
 *  The planetspec structure is defined in epic_datatypes.h.
 *  It contains the following members (mks SI units, unless otherwise noted):
 *    {index,
 *     name,type,
 *     orbital_epoch,
 *     re[m],rp[m],
 *     obliquity[deg],
 *     omega_sidereal,omega_synodic,wlon_noon_J2000,
 *     cp,cpr=cp/rgas,rgas,p0,kappa=1/cpr,
 *     GM,J2,
 *     x_h2,x_he,x_3,
 *     a[AU],e,i[deg],
 *     lon_ascending_node[deg],lon_perihelion[deg],mean_lon[deg],
 *     orbit_period[yrs],
 *     vernal_equinox_anomaly [deg],
 *     kinvisc[m^2/s],dynvisc[kg/m/s],k_a[J/m/s/K],
 *     u(p,lat[deg]),
 *     cloud[MAX_NUM_SPECIES]}
 *
 *  The Keplarian elements are taken from p. 316, Table 5.8.1 of Seidelmann, 1992, 
 *  the Explanatory Supplement to the Astronomical Almanac, University Science Books.
 *  The elements for satellites are of their parent planet.
 *
 *  The dynamic-west longitude [0,360] where it is mean-sun noon at time J2000.0
 *  (JD 2451545.0 TT), wlon_noon_J2000, was determined using the JPL Horizons web interface,
 *  with the observer at the Sun, and each planet in turn as the target.  Horizons
 *  gives the subsolar longitude as seen from the Sun, to which we add the
 *  one-way light travel time multipled by the planet's synotic rotation rate,
 *  to obtain the longitude of local noon on the planet at J2000.0.  We then use the equation
 *  of time to adjust it to the mean-sun position.  The one variation is Earth
 *  (which would also apply to the Moon, and the Sun itself), for which IAU uses east longitude [-180,180]; EPIC
 *  does, too, but the parameter wlon_noon_J2000 is west longitude even for Earth, to be consistent.
 *
 *  Only data for the constants in the structure are specified here;
 *  the function pointer u() is set elsewhere.
 *
 *  The planetspec values of cpr (= 1/kappa) and cp are low-temperature-limit
 *  reference values, and should not be used otherwise.  Use return_cp()
 *  to get cp for a given thermodynamical state.
 *
 *  NORTH-SOUTH CONVENTION:
 *  Venus, Uranus and Pluto are confusing, because IAU south equals dynamic north, 
 *  and IAU east equals dynamic west.
 *  There is no simple answer because the IAU convention is more convenient when talking about 
 *  planetary orbits, but the dynamic convention is more convenient when talking about the
 *  weather on the planet.  Because EPIC is an atmospheric dynamics model, we use dynamic north and east,
 *  which equals IAM south and west when the IAU obliquity is greater than 90deg.
 *  Thus, in EPIC the Coriolis parameter is always positive in the nothern hemisphere and the Sun
 *  always rises in the east.
 *    Consequences of IAU obliquity greater than 90deg:
 *      planet->obliquity              = 180deg-obliquity(IAU)
 *      planet->vernal_equinox_anomaly = 180+vernal_equinox_anomaly(IAU)
 *                              => L_s = 180+L_s(IAU)
 *  The ascending node of the planet's orbit is defined in the IAU sense,
 *  regardless of the planet's obliquity, as is the longitude of the ascending node.
 *  The longitude of perihelion and the mean longitude are both compound angles, meaning
 *  they are a sum of angles measured in different planes, and hence are listed with IAU-convention values.
 *
 *  The gravitational parameters GM and J2 are taken from
 *  The Planetary Scientist's Companion (Lodders and Fegley, 1998). 
 *
 *  Typical kinematic viscosities, kinvisc, are from D.C. Wilcox (2000)
 *  Basic Fluid Mechanics, DCW Industries, 689. These numbers correspond
 *  to 15C for the main component of each atmosphere.
 *
 *  The dynamic viscosity for Earth is taken from Fowler et al (1996, J. Cli.).
 *  For the gas giants, we use the value for hydrogen listed in
 *  Weast et al (1987, CRC Handbook of Chemistry and Physics). For CO_2 and N_2
 *  atmospheres, one can use the gas viscosity calculator at
 *  http://lmnoeng.com/Flow/GasViscosity.htm.
 *  The parameter k_a is a crude value for the thermal conductivity of dry air [W/m/K]. 
 *  NOTE: thermal conductivity should be a temperature-dependent function. 
 *  
 *  The quantity planet->rgas should refer to the dry-air gas constant,
 *  R_GAS/mu_dry, where mu_dry is the typical molar mass of dry air for
 *  the planet.
 *
 *  The member cloud[] is of type microphysics_spec, and its elements are set by 
 *  functions in the subdirectory epic/src/clouds.
 *
 * NOTE: All parameters are subject to revision and some are not yet accurate.
 *
 * NOTE: We declare the planet pointer globally, although it is usually passed
 *       as an argument as well.
 */
planetspec
  *planet,
   venus   = {VENUS_INDEX,
              "Venus","terrestrial",
              "J2000",
              6052.e+3,6052.e+3,
              2.7,                             /* IAU: 177.3 */
              2.9924e-7,6.2289e-7,200.51,
              860.,4.50,191.,92.e+5,0.222,
              0.32486e+15,6.e-6,
              0.,0.,1.,
              0.72333199,0.00677323,3.39471,
              76.68069,131.53298,181.97973,
              0.61521,
              286.29,                          /* IAU 106.29 */
              7.84e-6,1.456e-5,.1},
   earth   = {EARTH_INDEX,
              "Earth","terrestrial",
              "J2000",
              6378.e+3,6357.e+3,
              23.45,
              7.2921e-5,7.272205e-5,359.96,
              1004.,3.48,287.1,1.e+5,0.287,
              0.39860e+15,1082.636e-6,
              0.,0.,1.,
              1.00000011,0.01671022,0.00005,
              -11.26064,102.94719,100.46435,
              1.00004,
              77.059,
              14.60e-6,1.718e-5,2.43e-2},
   mars    = {MARS_INDEX,
              "Mars","terrestrial",
              "J2000",
              3396.19e+3,3376.2e+3,
              25.19,
              7.0882e-5,7.0776e-5,43.26,
              735.,3.83,192.,610.,0.261,   /* cp = 735 at 200 K, varies with temp.; R. Haberle, email 8/10/15 */
              0.042828375e+15,1960.454e-6,
              0.,0.,1.,
              1.52366231,0.09341233,1.85061,
              49.57854,336.04084,355.45332,
              1.88089,
              109.021,
              7.84e-6,1.456e-5,.1},
   jupiter = {JUPITER_INDEX,
              "Jupiter","gas-giant",
              "J2000",
              71492.e+3,66852.e+3,
              3.12,
              1.7585e-4,1.7584e-4,69.54,
              9092.8,2.50,3637.,1.e+5,.4,
              126.6865e+15,14697.0e-6,
              0.864,0.136,0.,
              5.20336301,0.04839266,1.30530,
              100.55615,14.75385,34.40438,
              11.8623,
              302.592,
              101.00e-6,3.0e-6,0.2},
   saturn  = {SATURN_INDEX,
              "Saturn","gas-giant",
              "J2000",
              60330.e+3,54180.e+3,
              26.73,
              1.65117e-4,1.65117e-4,298.79,  /* Using System IIIw rotation period = 10h34m13s, not Voyager SKR System III = 10h39m24s */
              10976.,2.939,3735.,1.e+5,.34,  /* Using avg_mu = 2*.877+4*.118 = 2.226 */
              37.931208e+15,16291.9e-6,      /* GM Jacobson et al (2006, Astron. J.), J2 after Helled et al (2009, Icarus), with System IIIw */
              0.877,0.118,0.005,             /* Conrath et al (2000, Icarus 144, 124-134) */
              9.53707032,0.05415060,2.48446,
              113.71504,92.43194,49.94432,
              29.458,
              80.998,
              101.00e-6,3.0e-6,0.2},
   titan   = {TITAN_INDEX,
              "Titan","terrestrial",
              "J2000",
              2574.8e+3,2574.8e+3,
              0.27,                            /* Titan obliquity from http://www.gtti.it/GTTI08/papers/persi.pdf */
              4.5608e-6,4.5608e-6,87.45, 
              1044.,3.60,290.,1.467e+5,.2778,
              0.008978e+15,33.599e-6,          /* J2 from Iess et al. (2012), Table 2, SOL1a. */
              0.,0.,1.,
              9.53707032,0.05415060,2.48446,
              113.71504,92.43194,49.94432,
              29.458,
              80.998,
              1.49e-6,0.805e-5,.02},
   uranus  = {URANUS_INDEX,
              "Uranus","gas-giant",
              "J2000",
              25560.e+3,24972.e+3,
              82.14,                            /* IAU: 97.86 */
              1.0124e-4,1.0124e-4,355.62,
              9280.,2.67,3480.,1.e+5,.375,
              5.79395e+15,3516.0e-6,
              0.85,0.15,0.,
              19.19126393,0.04716771,0.76986,
              74.22988,170.96424,313.23218,
              84.01,
              356.546,                          /* IAU: 176.546 */
              101.00e-6,3.0e-6,0.1},
   neptune = {NEPTUNE_INDEX,
              "Neptune","gas-giant",
              "J2000",
              24764.e+3,24343.e+3,
              29.56,
              1.0834e-4,1.0834e-4,153.64,
              9280.,2.67,3480.,1.e+5,.375,
              6.83473e+15,3538.0e-6,
              0.81,0.19,0.,
              30.06896348,0.00858587,1.76917,
              131.72169,44.97135,304.88003,
              164.79,
              0.515,
              101.00e-6,3.0e-6,0.1},
   pluto   = {PLUTO_INDEX,
              "Pluto","terrestrial",
              "J2000",
              1152.e+3,1152.e+3,
              62.0,                             /* IAU: 118.0 NOTE: Need to update obliquity for Pluto. */
              1.1386e-5,1.1386e-5,133.68,
              1.,1.,1.,0.3,1.,                      /* NOTE: Need data for Pluto. */
              0.000884e+15,0.,                  /* NOTE: Need J2 for Pluto. */        
              0.,0.,1.,
              39.48168677,.24880766,17.14175,
              110.30347,224.06676,238.92881,
              248.6,
              355.545,                          /* IAU: 175.545 */
              14.50e-6,1.724e-5,0.02},
   hot_jupiter = {HOT_JUPITER_INDEX,            /* NOTE: These are just Jupiter placeholder values */
              "Hot_Jupiter","gas-giant",
              "J2000",
              71492.e+3,66852.e+3,
              3.12,
              1.7585e-4,1.7584e-4,69.54,
              9092.8,2.50,3637.,1.e+5,.4,
              126.6865e+15,14697.0e-6,
              0.864,0.136,0.,
              5.20336301,0.04839266,1.30530,
              100.55615,14.75385,34.40438,
              11.8623,
              302.592,
              101.00e-6,3.0e-6,0.3};
/*
 * The held_suarez case is based on
 *   Held, IM and MJ Suarez, 1994, A proposal for the intercomparison of the dynamical
 *      cores of atmospheric general circulation models, Bull. Am. Meteorol. Soc. 
 *      73, 1825-1830.
 */
planetspec 
  held_suarez = {HELD_SUAREZ_INDEX,
                 "Held_Suarez","terrestrial",
                 "J2000",
                 6371.e+3,6371.e+3,
                 23.45,
                 7.2921e-5,7.272205e-5,359.13,
                 1004.,3.5,286.86,1.e+5,0.2857,
                 0.39860e+15,1082.636e-6,
                 0.,0.,1.,
                 1.00000011,0.01671022,0.00005,
                 -11.26064,102.94719,100.46435,
                 1.0,
                 77.059,
                 14.60e-6,1.718e-5,2.43e-2};
/*
 * The goldmine testbed case is based on Held-Suarez model
 */
planetspec 
  goldmine = {GOLDMINE_INDEX,
              "Goldmine","terrestrial",
              "J2000",
              6371.e+3,6371.e+3,
              23.45,
              7.2921e-5,7.272205e-5,359.13,
              1004.,3.5,286.86,1.e+5,0.2857,
              0.39860e+15,1082.636e-6,
              0.,0.,1.,
              1.00000011,0.01671022,0.00005,
              -11.26064,102.94719,100.46435,
              1.0,
              77.059,
              14.60e-6,1.718e-5,2.43e-2};
/*
 * The venus_llr05 case is based on
 *   Lee, C., S.R. Lewis and P.L. Read, 2005, A numerical model of the atmosphere of Venus,
 *       Adv. Space Res. 36, 2142-2145.
 * The gas constant and specific heat used in this study are not listed in the 2005 paper.
 * In an email from Chris Lee, 9/30/05 <leec@atm.ox.ac.uk>, he mentioned he used 
 * rgas = 185., which implies mu = 44.94, and is what we used in Herrnstein and Dowling (2007, JGR).
 * We subsequently discovered (1/17/08) that this is about 3% higher than is consistent with
 * the T and p data in epic/data/venus/VIRA.dat (Venus International Reference Atmosphere; Herrnstein and Dowling
 * used the equatorial data in epic/data/venus/archive/Venus_book.dat, which has the same problem),
 * and causes a 3-4% bias between the input T(p) and the model's calculation of T(p,rho).
 * We have now switched to mu = 43.53 and the corresponding rgas = 191.0 for venus_llr05.
 * We have also switched from Lee's cp = 840. to 860., to match our Venus case.
 */
planetspec
  venus_llr05 = {VENUS_LLR05_INDEX,
                 "Venus_LLR05","terrestrial",
                 "J2000",
                 6052.e+3,6052.e+3,
                 2.7,                     /* IAU: 177.3 */
                 2.9924e-7,6.2289e-7,200.51,
                 860.,4.50,191.,92.e+5,0.222,
                 0.32486e+15,6.e-6,
                 0.,0.,1.,
                 0.72333199,0.00677323,3.39471,
                 76.68069,131.53298,181.97973,
                 .61521,
                 286.29,                  /* IAU: 106.29 */
                 7.84e-6,1.456e-5,.1};

const chem_element
  /*
   * Chemical element data.
   * Solar abundances from
   *   Grevesse, N. and A.J. Sauval, 2001, Solar Abundances, in Encyclopedia of 
   *     Astronomy & Astrophysics, Nature Publishing Group.
   * Abundances are in the form A_el =log N_el/N_H+12.0, such that A_H is 12.00
   * by definition and N_el is the abundance by number.
   */
  ElementGS01[LAST_ATOMIC_NUMBER+1] = {
    {  0," ", " ",             0.,  -10.00},
    {  1,"H", "hydrogen",      1.008,12.00},
    {  2,"He","helium",        4.003,10.93},
    {  3,"Li","lithium",       6.941, 1.10},
    {  4,"Be","beryllium",     9.012, 1.40},
    {  5,"B", "boron",        10.81,  2.55},
    {  6,"C", "carbon",       12.01,  8.52},
    {  7,"N", "nitrogen",     14.01,  7.92},
    {  8,"O", "oxygen",       16.00,  8.83},
    {  9,"F", "fluorine",     19.00,  4.56},
    { 10,"Ne","neon",         20.18,  8.08},
    { 11,"Na","sodium",       22.99,  6.33},
    { 12,"Mg","magnesium",    24.30,  7.58},
    { 13,"Al","aluminum",     26.98,  6.47},
    { 14,"Si","silicon",      28.09,  7.55},
    { 15,"P", "phosphorus",   30.97,  5.45},
    { 16,"S", "sulfur",       32.07,  7.33},
    { 17,"Cl","chlorine",     35.45,  5.50},
    { 18,"Ar","argon",        39.95,  6.40},
    { 19,"K", "potassium",    39.10,  5.12},
    { 20,"Ca","calcium",      40.08,  6.36},
    { 21,"Sc","scandium",     44.96,  3.17},
    { 22,"Ti","titanium",     47.88,  5.02},
    { 23,"V", "vanadium",     50.94,  4.00},
    { 24,"Cr","chromium",     52.00,  5.67},
    { 25,"Mn","manganese",    54.94,  5.39},
    { 26,"Fe","iron",         55.85,  7.50},
    { 27,"Co","cobalt",       58.93,  4.92},
    { 28,"Ni","nickel",       58.69,  6.25},
    { 29,"Cu","copper",       63.55,  4.21},
    { 30,"Zn","zinc",         65.39,  4.60},
    { 31,"Ga","gallium",      69.72,  2.88},
    { 32,"Ge","germanium",    72.61,  3.41},
    /* Next 4 abundances are meteorite values: */
    { 33,"As","arsenic",      74.92,  2.37},
    { 34,"Se","selenium",     78.96,  3.41},
    { 35,"Br","bromine",      79.90,  2.63},
    { 36,"Kr","krypton",      83.80,  3.31},
    { 37,"Rb","rubidium",     85.47,  2.60},
    { 38,"Sr","strontium",    87.62,  2.97},
    { 39,"Y", "yttrium",      88.91,  2.24},
    { 40,"Zr","zirconium",    91.22,  2.60},
    { 41,"Nb","niobium",      92.91,  1.42},
    { 42,"Mo","molybdenum",   95.94,  1.92},
    { 43,"Tc","technetium",   98.91,-10.00},
    { 44,"Ru","ruthenium",   101.1,   1.84},
    { 45,"Rh","rhodium",     102.9,   1.12},
    { 46,"Pd","palladium",   106.4,   1.69},
    { 47,"Ag","silver",      107.9,   0.94},
    { 48,"Cd","cadmium",     112.4,   1.77},
    { 49,"In","indium",      114.8,   1.66},
    { 50,"Sn","tin",         118.7,   2.0 },
    { 51,"Sb","antimony",    121.8,   1.0 },
    /* Next 4 abundances are meteorite values: */
    { 52,"Te","tellurium",   127.6,   2.24},
    { 53,"I", "iodine",      126.9,   1.51},
    { 54,"Xe","xenon",       131.3,   2.17},
    { 55,"Cs","cesium",      132.9,   1.13},
    { 56,"Ba","barium",      137.3,   2.13},
    { 57,"La","lanthanum",   138.9,   1.17},
    { 58,"Ce","cerium",      140.1,   1.58},
    { 59,"Pr","praseodymium",140.9,   0.71},
    { 60,"Nd","neodymium",   144.2,   1.50},
    { 61,"Pm","promethium",  144.9, -10.00},
    { 62,"Sm","samarium",    150.4,   1.01},
    { 63,"Eu","europium",    152.0,   0.51},
    { 64,"Gd","gadolinium",  157.2,   1.12},
    { 65,"Tb","terbium",     158.9,  -0.1 },
    { 66,"Dy","dysprosium",  162.5,   1.14},
    { 67,"Ho","holmium",     164.9,   0.26},
    { 68,"Er","erbium",      167.3,   0.93},
    { 69,"Tm","thulium",     168.9,   0.00},
    { 70,"Yb","ytterbium",   173.0,   1.08},
    { 71,"Lu","lutetium",    175.0,   0.06},
    { 72,"Hf","hafnium",     178.5,   0.88},
    /* Next abundance is the meteorite value: */
    { 73,"Ta","tantalum",    180.9,  -0.13},
    { 74,"W", "tungsten",    183.8,   1.11},
    /* Next abundance is the meteorite value: */
    { 75,"Re","rhenium",     186.2,   0.28},
    { 76,"Os","osmium",      190.2,   1.45},
    { 77,"Ir","iridium",     192.2,   1.35},
    { 78,"Pt","platinum",    195.1,   1.8 },
    { 79,"Au","gold",        197.0,   1.01},
    /* Next abundance is the meteorite value: */
    { 80,"Hg","mercury",     200.6,   1.13},
    { 81,"Tl","thallium",    204.4,   0.9 },
    { 82,"Pb","lead",        207.2,   1.95},
    /* Next abundances is the meteorite value: */
    { 83,"Bi","bismuth",     209.0,   0.71},
    { 84,"Po","polonium",    210.0, -10.00},
    { 85,"At","astatine",    210.0, -10.00},
    { 86,"Rn","radon",       222.0, -10.00},
    { 87,"Fr","francium",    223.0, -10.00},
    { 88,"Ra","radium",      226.0, -10.00},
    { 89,"Ac","actinium",    227.0, -10.00},
	/* Next abundance is the meteorite value: */
    { 90,"Th","thorium",     232.0,   0.09},
    { 91,"Pa","protactinium",231.0, -10.00},
    /* Next abundance is the meteorite value: */
    { 92,"U", "uranium",     238.0,  -0.50},
    { 93,"Np","neptunium",   237.0, -10.00},
    { 94,"Pu","plutonium",   239.1, -10.00},
    { 95,"Am","americium",   243.1, -10.00},
    { 96,"Cm","curium",      247.1, -10.00},
    { 97,"Bk","berkelium",   247.1, -10.00},
    { 98,"Cf","californium", 252.1, -10.00},
    { 99,"Es","einsteinium", 252.1, -10.00},
    {100,"Fm","fermium",     257.1, -10.00},
    {101,"Md","mendelevium", 256.1, -10.00},
    {102,"No","nobelium",    259.1, -10.00},
    {103,"Lr","lawrencium",  260.1, -10.00}
  },
  /*
   * Chemical element data.
   * Solar abundances from
   *   Anders, E. and N. Grevesse, 1989, Abundances of the elements:
   *     Meteoritic and solar, Geochim. Cosmochim. Acta 53, 197-214
   */
  ElementAG89[LAST_ATOMIC_NUMBER+1] = {
    {  0," ", " ",             0.,  -10.00},
    {  1,"H", "hydrogen",      1.008,12.00},
    {  2,"He","helium",        4.003,10.99},
    {  3,"Li","lithium",       6.941, 1.16},
    {  4,"Be","beryllium",     9.012, 1.15},
    {  5,"B", "boron",        10.81,  2.6 },
    {  6,"C", "carbon",       12.01,  8.56},
    {  7,"N", "nitrogen",     14.01,  8.05},
    {  8,"O", "oxygen",       16.00,  8.93},
    {  9,"F", "fluorine",     19.00,  4.56},
    { 10,"Ne","neon",         20.18,  8.09},
    { 11,"Na","sodium",       22.99,  6.33},
    { 12,"Mg","magnesium",    24.30,  7.58},
    { 13,"Al","aluminum",     26.98,  6.47},
    { 14,"Si","silicon",      28.09,  7.55},
    { 15,"P", "phosphorus",   30.97,  5.45},
    { 16,"S", "sulfur",       32.07,  7.21},
    { 17,"Cl","chlorine",     35.45,  5.5 },
    { 18,"Ar","argon",        39.95,  6.56},
    { 19,"K", "potassium",    39.10,  5.12},
    { 20,"Ca","calcium",      40.08,  6.36},
    { 21,"Sc","scandium",     44.96,  3.10},
    { 22,"Ti","titanium",     47.88,  4.99},
    { 23,"V", "vanadium",     50.94,  4.00},
    { 24,"Cr","chromium",     52.00,  5.67},
    { 25,"Mn","manganese",    54.94,  5.39},
    { 26,"Fe","iron",         55.85,  7.67},
    { 27,"Co","cobalt",       58.93,  4.92},
    { 28,"Ni","nickel",       58.69,  6.25},
    { 29,"Cu","copper",       63.55,  4.21},
    { 30,"Zn","zinc",         65.39,  4.60},
    { 31,"Ga","gallium",      69.72,  2.88},
    { 32,"Ge","germanium",    72.61,  3.41},
    /* Next 4 abundances are meteorite values: */
    { 33,"As","arsenic",      74.92,  2.37},
    { 34,"Se","selenium",     78.96,  3.35},
    { 35,"Br","bromine",      79.90,  2.63},
    { 36,"Kr","krypton",      83.80,  3.23},
    { 37,"Rb","rubidium",     85.47,  2.60},
    { 38,"Sr","strontium",    87.62,  2.90},
    { 39,"Y", "yttrium",      88.91,  2.24},
    { 40,"Zr","zirconium",    91.22,  2.60},
    { 41,"Nb","niobium",      92.91,  1.42},
    { 42,"Mo","molybdenum",   95.94,  1.92},
    { 43,"Tc","technetium",   98.91,-10.00},
    { 44,"Ru","ruthenium",   101.1,   1.84},
    { 45,"Rh","rhodium",     102.9,   1.12},
    { 46,"Pd","palladium",   106.4,   1.69},
    { 47,"Ag","silver",      107.9,   0.94},
    { 48,"Cd","cadmium",     112.4,   1.86},
    { 49,"In","indium",      114.8,   1.66},
    { 50,"Sn","tin",         118.7,   2.0 },
    { 51,"Sb","antimony",    121.8,   1.0 },
    /* Next 4 abundances are meteorite values: */
    { 52,"Te","tellurium",   127.6,   2.24},
    { 53,"I", "iodine",      126.9,   1.51},
    { 54,"Xe","xenon",       131.3,   2.23},
    { 55,"Cs","cesium",      132.9,   1.12},
    { 56,"Ba","barium",      137.3,   2.13},
    { 57,"La","lanthanum",   138.9,   1.22},
    { 58,"Ce","cerium",      140.1,   1.55},
    { 59,"Pr","praseodymium",140.9,   0.71},
    { 60,"Nd","neodymium",   144.2,   1.50},
    { 61,"Pm","promethium",  144.9, -10.00},
    { 62,"Sm","samarium",    150.4,   1.00},
    { 63,"Eu","europium",    152.0,   0.51},
    { 64,"Gd","gadolinium",  157.2,   1.12},
    { 65,"Tb","terbium",     158.9,  -0.1 },
    { 66,"Dy","dysprosium",  162.5,   1.1 },
    { 67,"Ho","holmium",     164.9,   0.26},
    { 68,"Er","erbium",      167.3,   0.93},
    { 69,"Tm","thulium",     168.9,   0.00},
    { 70,"Yb","ytterbium",   173.0,   1.08},
    { 71,"Lu","lutetium",    175.0,   0.76},
    { 72,"Hf","hafnium",     178.5,   0.88},
    /* Next abundance is the meteorite value: */
    { 73,"Ta","tantalum",    180.9,  -0.13},
    { 74,"W", "tungsten",    183.8,   1.11},
    /* Next abundance is the meteorite value: */
    { 75,"Re","rhenium",     186.2,   0.27},
    { 76,"Os","osmium",      190.2,   1.45},
    { 77,"Ir","iridium",     192.2,   1.35},
    { 78,"Pt","platinum",    195.1,   1.8 },
    { 79,"Au","gold",        197.0,   1.01},
    /* Next abundance is the meteorite value: */
    { 80,"Hg","mercury",     200.6,   1.09},
    { 81,"Tl","thallium",    204.4,   0.9 },
    { 82,"Pb","lead",        207.2,   1.85},
    /* Next abundance is the meteorite value: */
    { 83,"Bi","bismuth",     209.0,   0.71},
    { 84,"Po","polonium",    210.0, -10.00},
    { 85,"At","astatine",    210.0, -10.00},
    { 86,"Rn","radon",       222.0, -10.00},
    { 87,"Fr","francium",    223.0, -10.00},
    { 88,"Ra","radium",      226.0, -10.00},
    { 89,"Ac","actinium",    227.0, -10.00},
    { 90,"Th","thorium",     232.0,   0.12},
    { 91,"Pa","protactinium",231.0, -10.00},
    /* Next abundance is the meteorite value: */
    { 92,"U", "uranium",     238.0,  -0.49},
    { 93,"Np","neptunium",   237.0, -10.00},
    { 94,"Pu","plutonium",   239.1, -10.00},
    { 95,"Am","americium",   243.1, -10.00},
    { 96,"Cm","curium",      247.1, -10.00},
    { 97,"Bk","berkelium",   247.1, -10.00},
    { 98,"Cf","californium", 252.1, -10.00},
    { 99,"Es","einsteinium", 252.1, -10.00},
    {100,"Fm","fermium",     257.1, -10.00},
    {101,"Md","mendelevium", 256.1, -10.00},
    {102,"No","nobelium",    259.1, -10.00},
    {103,"Lr","lawrencium",  260.1, -10.00}
  },
   /*
    * Assign which solar-abundance data set to use.
    */
  *Element = ElementAG89;

/*
 * Declare globally the following shift integers to speed up multidimensional 
 * array referencing. Set Ishift = ILO-IPAD, etc. in make_arrays().
 */
int
  Ishift,Jshift,Kshift,
  Iadim,Jadim,Kadim,
  Nelem2d,Nelem3d,
  Shift2d,Shift3d,Shiftkj;

/*
 * Global variables used to communicate with rho_minus_rho().
 */
EPIC_FLOAT
  RHOMRHO_fp,
  RHOMRHO_p,
  RHOMRHO_mu,
  RHOMRHO_density;
planetspec
  *RHOMRHO_planet;

/*
 * Global variables used to communicate with th_minus_th_p() and
 * th_minus_th_t().
 */
EPIC_FLOAT
  THMTH_fp,
  THMTH_temperature,
  THMTH_theta,
  THMTH_p;
planetspec
  *THMTH_planet;

/*
 * Global variables used to communicate with fpe_minus_fpe().
 */
EPIC_FLOAT
  FPEMFPE_p,
  FPEMFPE_theta;
planetspec
  *FPEMFPE_planet;

/*
 * Global variables used to communicate with sgth_minus_sgth().
 */
EPIC_FLOAT
  SGTHMSGTH_theta,
  SGTHMSGTH_sigmatheta;

/*
 * Global variables used to communicate with equatorial_radius_cubic()
 */
double
  ERC_J2,
  ERC_c3w2_GM;

/*
 * Global variables used to communicate with polar_radius_cubic()
 */
double
  PRC_J2,
  PRC_a3w2_GM;

/* * * * * * * * * * * * end of epic_globals.c * * * * * * * * * * * * * * * */
