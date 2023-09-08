/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2002-2019 Csaba Palotai *A*                       *
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

/* * * * * * * * * *  epic_microphysics_funcs.c  * * * * * * * * * *
 *                         v.2.1                                   *        
 *                                                                 *
 *       Functions governing the hydrological cycle.               *
 *       This file includes the following functions:               *
 *                                                                 *
 *           set_species_thermo_data()                             *
 *           set_microphysics_params()                             *
 *           enthalpy_change()                                     *
 *           read_enthalpy_change_data()                           *
 *           enthalpy_change_H_2O(), etc.                          *
 *           sat_vapor_p_H_2O(), etc.                              *
 *           dynvisc()                                             *
 *           conductivity()                                        *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*======================= set_species_thermo_data() =============================*/

void set_species_thermo_data(void)
{
  /*
   * Set species vapor-liquid-solid triple point temperature [K]
   * and pressure [Pa].
   *
   * grub:gurn:burlap: Need values for FLOAT_MAX placeholders,
   *                   and need to check other values; need references
   *
   * Lf values:   Constant latent heat of fusion for substance
   *              Used for melting/freezing processes
   *              Unit: J/kg
   */

  /* 
   * http://www.trgn.com/database/cryogen.htm
   */
  T_triple_pt(H_2O_INDEX) =  273.16;
  p_triple_pt(H_2O_INDEX) =  607.95;
  Lf(H_2O_INDEX)          =  333.5e+3;  /* [J/kg] */
  Lv(H_2O_INDEX)          = 2498.7e+3; 
  Ls(H_2O_INDEX)          = 2832.2e+3;  /* Average of values from data table. */

  /* http://www.trgn.com/database/cryogen.htm */
  T_triple_pt(NH_3_INDEX) =  195.5;      /* 195.49 */
  p_triple_pt(NH_3_INDEX) =   66.2;
  Lf(NH_3_INDEX)          =  339.0e+3;
  Lv(NH_3_INDEX)          = 1490.2e+3;
  Ls(NH_3_INDEX)          = 1829.2e+3;

  /* 
   * Lodders and Fegley (1998), inferred from Table 1.20
   * http://www.encyclopedia.airliquide.com
   */
  T_triple_pt(H_2S_INDEX)   = 187.61;
  p_triple_pt(H_2S_INDEX)   = 23295.;
  Lf(H_2S_INDEX)            =    69.75e+3;
  Lv(H_2S_INDEX)            = FLOAT_MAX;
  Ls(H_2S_INDEX)            = FLOAT_MAX;

  /* http://www.trgn.com/database/cryogen.htm */
  T_triple_pt(CH_4_INDEX)   = 88.7;
  p_triple_pt(CH_4_INDEX)   = 10031.;
  Lf(CH_4_INDEX)            =    58.7e+3;
  Lv(CH_4_INDEX)            = FLOAT_MAX;
  Ls(CH_4_INDEX)            = FLOAT_MAX;

  /* http://encyclopedia.airliquide.com */
  T_triple_pt(C_2H_2_INDEX) = 192.55;
  p_triple_pt(C_2H_2_INDEX) = 1.282e+5;
  Lf(C_2H_2_INDEX)          = FLOAT_MAX;
  Lv(C_2H_2_INDEX)          = 801.9e+3;;
  Ls(C_2H_2_INDEX)          = FLOAT_MAX;

  /* 
   * http://www.nist.gov/data/PDFfiles/jpcrd294.pdf (Jahangiri et al, 1986, J. Phys. Chem. Ref. Data 15
   */
  T_triple_pt(C_2H_4_INDEX) = 103.986;
  p_triple_pt(C_2H_4_INDEX) = 122.5;
  Lf(C_2H_4_INDEX)          = 119.37e+3;
  Lv(C_2H_4_INDEX)          = 482.86e+3;
  Ls(C_2H_4_INDEX)          = 602.23e+3;

  /* http://encyclopedia.airliquide.com */
  T_triple_pt(C_2H_6_INDEX) = FLOAT_MAX;
  p_triple_pt(C_2H_6_INDEX) = FLOAT_MAX;
  Lf(C_2H_6_INDEX)          = 94.977e+3;
  Ls(C_2H_6_INDEX)          = 488.76e+3;

  /*
   * http://www.trgn.com/database/cryogen.htm
   * Ls, Forget et al. (1998) 
   */
  T_triple_pt(CO_2_INDEX)   = 216.6;
  p_triple_pt(CO_2_INDEX)   = 518800.;
  Lf(CO_2_INDEX)            = FLOAT_MAX;
  Lv(CO_2_INDEX)            = FLOAT_MAX;
  Ls(CO_2_INDEX)            =    590.e+3;

  /*
   * Lodders and Fegley (1998), inferred from Table 1.20
   * NOTE: Need to check this, T=317 may just be the top range
   *       of validity of the solid-vapor curve fit.
   */
  T_triple_pt(NH_4SH_INDEX) = 317.;
  p_triple_pt(NH_4SH_INDEX) = 98600.;
  Lf(NH_4SH_INDEX)          = FLOAT_MAX;
  Lv(NH_4SH_INDEX)          = FLOAT_MAX;
  Ls(NH_4SH_INDEX)          = FLOAT_MAX;

  /* http://www.e-cats.com/databook/Page%2022.htm */
  T_triple_pt(O_3_INDEX)    = 80.65;
  p_triple_pt(O_3_INDEX)    = 11400.;
  Lf(O_3_INDEX)             = FLOAT_MAX;
  Lv(O_3_INDEX)             = FLOAT_MAX;
  Ls(O_3_INDEX)             = FLOAT_MAX;

  /* http://www.trgn.com/database/cryogen.htm */
  T_triple_pt(N_2_INDEX)    = 63.2;
  p_triple_pt(N_2_INDEX)    = 12870.;
  Lf(N_2_INDEX)             = FLOAT_MAX;
  Lv(N_2_INDEX)             = FLOAT_MAX;
  Ls(N_2_INDEX)             = FLOAT_MAX;

  /*
   * http://webbook.nist.gov/cgi/cbook.cgi?ID=C7803512&Units=SI&Mask=7
   * Fluck, The Chemistry of Phosphine 
   */
  T_triple_pt(PH_3_INDEX)   = 139.41;
  p_triple_pt(PH_3_INDEX)   = FLOAT_MAX;
  Lf(PH_3_INDEX)            = 160.89e+3;
  Lv(PH_3_INDEX)            = FLOAT_MAX;
  Ls(PH_3_INDEX)            = FLOAT_MAX;
 
  return;
}

/*======================= end of set_species_thermo_data() ======================*/

/*======================= set_microphysics_params() =============================*/

  /*
   * Constants and coefficients for microphysical processes.
   * Call after set_species_thermo_data().
   *
   * NOTE: If a value is not known, flag it by setting to "1./0." (not to "NaN").
   *
   * Note: 10/2/03 CJP
   *       Some of them are just "made up" values thus no citation is available
   *       more data or sensitivity tests are needed for the proper value
   *
   * densities of rain [kg/m3]
   * densities of snow [kg/m3]  (Using 50% of bulk ice values)
   * T_00 values: Supercooled liquid threshold temperatures
   *              H_20         : -20C
   *              other species:  use triple point
   */

void set_microphysics_params(int planet_index)
{
  int
    is;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_microphysics_params";

  switch(planet_index) {
    case VENUS_INDEX:
    case VENUS_LLR05_INDEX:
      is = CO_2_INDEX;
      T_00(is) = T_triple_pt(is);
    break;

    case EARTH_INDEX:
      is = H_2O_INDEX;
      RHO_RAIN(  is) =  1000.;        /* Fowler et al. (1996) */
      RHO_SNOW(  is) =   917.0*0.5;   /* Perry's Chem.Eng.Handbook p.2-304. */
      COEFF_XI(  is) =   14900.0;     /* Hong et al. 2004 */
      COEFF_YI(  is) =     1.31;
      COEFF_XS(  is) =    11.72;      /* Dudhia 1989 */ 
      COEFF_YS(  is) =     0.41; 
      COEFF_XR(  is) =  842.0;        /* Dudhia 1989 */
      COEFF_YR(  is) =    0.8; 
      COEFF_A(   is) =   13.16;
      COEFF_B(   is) =    0.16;
      COEFF_C(   is) =    5.38e+7;
      COEFF_D(   is) =    0.75;
      COEFF_M(   is) =    7.06165e-3; /* HDC04 */ 
      COEFF_N(   is) =    2.0;
      N_0R(      is) =    8.e+6;      /* Intercept value in raindrop size distribution [m-4], FRR96 p526, D89 p3103 */
      E_R(       is) =    1.0;        /* Accretion efficiency for rain, D89 p3104 */
      D_Icrit(   is) =  500.e-6;      /* Maximum diameter of cloud ice crystal H04 p108 */
      f1r(       is) =    0.78;       /* Coefficient in the ventillation factor, D89 p3103 */
      f2r(       is) =    0.308;      /* Coefficient in the ventillation factor, D89 p3103 */
      f1s(       is) =    0.65;       /* Coefficient in the ventillation factor, D89 p3103 */
      f2s(       is) =    0.44;       /* Coefficient in the ventillation factor, D89 p3103 */
      P_EXP_LIQ( is) =    0.4;        /* Exponent in the (p0/p)^x term for rain: liquid cloud particles do not fall at this point */
      P_EXP_ICE( is) =    0.4;        /* Exponent in the (p0/p)^x term for cloud ice */
      P_EXP_SNOW(is) =    0.4;        /* Exponent in the (p0/p)^x term for snow */
      ALPHA_RAUT(is) =    0.001;      /* Rate coeff. for autoconversion of cloud water [sec-1]. Fowler et al (1996) and Kessler (1969) */
      Q_LIQ_0(   is) =    0.00025;    /* Threshold value for rain autoconversion [kg kg-1]. Fowler et al. (1996) */
      T_00(      is) =  273.16;
    break;

    case HELD_SUAREZ_INDEX:
      is = H_2O_INDEX;
      RHO_RAIN(is) = 1./0.;
      RHO_SNOW(is) = 1./0.;
      T_00(    is) = 273.16;
    break;

    case GOLDMINE_INDEX:
      is = H_2O_INDEX;
      RHO_RAIN(is) = 1./0.;
      RHO_SNOW(is) = 1./0.;
      T_00(    is) = 273.16;
    break;

    case MARS_INDEX:
      is = CO_2_INDEX;
      RHO_RAIN(is) = 1./0.;
      RHO_SNOW(is) = 1./0.;  /* ice denisty: ~1630kg/m3 at 140K : Forget et al. 1998 */
      T_00(    is) = T_triple_pt(is);
    break;

    case JUPITER_INDEX:
      is = CH_4_INDEX;
        T_00(    is) = T_triple_pt(is);

      is = C_2H_2_INDEX;
        T_00(    is) = T_triple_pt(is);

      is = C_2H_6_INDEX;
        T_00(    is) = T_triple_pt(is);

      is = NH_3_INDEX;
      RHO_RAIN(  is) =   733.0;      /* From Bureau of Standards circ.No.142. See hardcopy for chart */ 
      RHO_SNOW(  is) =   786.8*0.5;  /* Yurtseven and Salihoglu, 2002. p.420 */
      COEFF_XI(  is) =   702.46;
      COEFF_YI(  is) =     0.8695;
      COEFF_XS(  is) =    53.065;
      COEFF_YS(  is) =     0.46179; 
      COEFF_XR(  is) =  2479.0;
      COEFF_YR(  is) =     0.76217;
      COEFF_A(   is) =     16.45;
      COEFF_B(   is) =     0.16;
      COEFF_C(   is) =     5.38e+7; 
      COEFF_D(   is) =     0.75;
      COEFF_M(   is) =     7.06165e-3*RHO_SNOW(is)/(0.5*917.0); /* HDC04, Perry's Chem. Eng. Handbook p.2-304 */
      COEFF_N(   is) =     2.0;
      N_0R(      is) =     8.e+6;    /* Intercept value in raindrop size distribution [m-4], FRR96 p526, D89 p3103 */
      E_R(       is) =     1.0;      /* Accretion efficiency for rain, D89 p3104 */
      D_Icrit(   is) =   500.e-6;    /* Maximum diameter of cloud ice crystal H04 p108 */
      f1r(       is) =     0.78;     /* Coefficient in the ventillation factor, D89 p3103 */
      f2r(       is) =     0.308;    /* Coefficient in the ventillation factor, D89 p3103 */
      f1s(       is) =     0.65;     /* Coefficient in the ventillation factor, D89 p3103 */
      f2s(       is) =     0.44;     /* Coefficient in the ventillation factor, D89 p3103 */
      P_EXP_LIQ( is) =     0.324;    /* Exponent in the (p0/p)^x term for rain: liquid cloud particles do not fall at this point */
      P_EXP_ICE( is) =     0.238;    /* Exponent in the (p0/p)^x term for cloud ice */
      P_EXP_SNOW(is) =     0.318;    /* Exponent in the (p0/p)^x term for snow. */
      ALPHA_RAUT(is) =     0.001;    /* Rate coeff. for autoconversion of cloud water [sec-1]. Fowler et al (1996) and Kessler (1969) */
      Q_LIQ_0(   is) =     0.00025;  /* Threshold value for rain autoconversion [kg kg-1]. Fowler et al. (1996) */
      T_00(      is) = T_triple_pt(is);

      is = H_2S_INDEX;
      T_00(    is) = T_triple_pt(is);

      is = NH_4SH_INDEX;
      RHO_RAIN(is) = 1./0.;
      RHO_SNOW(is) = 1./0.;
      T_00(    is) = T_triple_pt(is);

      is = H_2O_INDEX;
      RHO_RAIN(  is) =  1000.;        /* Fowler et al. (1996) */
      RHO_SNOW(  is) =   917.0*0.5;   /* Perry's Chem.Eng.Handbook p.2-304. */
      COEFF_XI(  is) =   735.5; 
      COEFF_YI(  is) =     0.86273;
      COEFF_XS(  is) =    54.96; 
      COEFF_YS(  is) =     0.45254; 
      COEFF_XR(  is) =  2615.1; 
      COEFF_YR(  is) =    0.74245; 
      COEFF_A(   is) =   13.16;
      COEFF_B(   is) =    0.16;
      COEFF_C(   is) =    5.38e+7;
      COEFF_D(   is) =    0.75;
      COEFF_M(   is) =    7.06165e-3; /* HDC04 */ 
      COEFF_N(   is) =    2.0;
      N_0R(      is) =    8.e+6;      /* Intercept value in raindrop size distribution [m-4], FRR96 p526, D89 p3103 */
      E_R(       is) =    1.0;        /* Accretion efficiency for rain, D89 p3104 */
      D_Icrit(   is) =  500.e-6;      /* Maximum diameter of cloud ice crystal H04 p108 */
      f1r(       is) =    0.78;       /* Coefficient in the ventillation factor, D89 p3103 */
      f2r(       is) =    0.308;      /* Coefficient in the ventillation factor, D89 p3103 */
      f1s(       is) =    0.65;       /* Coefficient in the ventillation factor, D89 p3103 */
      f2s(       is) =    0.44;       /* Coefficient in the ventillation factor, D89 p3103 */
      P_EXP_LIQ( is) =    0.331;      /* Exponent in the (p0/p)^x term for rain: liquid cloud particles do not fall at this point */
      P_EXP_ICE( is) =    0.238;      /* Exponent in the (p0/p)^x term for cloud ice */
      P_EXP_SNOW(is) =    0.318;      /* Exponent in the (p0/p)^x term for snow */
      ALPHA_RAUT(is) =    0.001;      /* Rate coeff. for autoconversion of cloud water [sec-1]. Fowler et al (1996) and Kessler (1969) */
      Q_LIQ_0(   is) =    0.00025;    /* Threshold value for rain autoconversion [kg kg-1]. Fowler et al. (1996) */
      T_00(      is) =  273.16;
    break;

    case SATURN_INDEX:
      is = NH_3_INDEX;
      RHO_RAIN(  is) =   733.0;      /* From Bureau of Standards circ.No.142. See hardcopy for chart */ 
      RHO_SNOW(  is) =   786.8*0.5;  /* Yurtseven and Salihoglu, 2002. p.420*/
      COEFF_XI(  is) =   365.32;
      COEFF_YI(  is) =     0.89124;
      COEFF_XS(  is) =    32.18;
      COEFF_YS(  is) =     0.4962; 
      COEFF_XR(  is) =  1652.3;
      COEFF_YR(  is) =     0.77046;
      COEFF_A(   is) =     16.45;
      COEFF_B(   is) =     0.16;
      COEFF_C(   is) =     5.38e+7; 
      COEFF_D(   is) =     0.75;
      COEFF_M(   is) =     7.06165e-3*RHO_SNOW(is)/(0.5*917.0); /* HDC04, Perry's Chem. Eng. Handbook p.2-304 */
      COEFF_N(   is) =     2.0;
      N_0R(      is) =     8.e+6;    /* Intercept value in raindrop size distribution [m-4], FRR96 p526, D89 p3103 */
      E_R(       is) =     1.0;      /* Accretion efficiency for rain, D89 p3104 */
      D_Icrit(   is) =   500.e-6;    /* Maximum diameter of cloud ice crystal H04 p108 */
      f1r(       is) =     0.78;     /* Coefficient in the ventillation factor, D89 p3103 */
      f2r(       is) =     0.308;    /* Coefficient in the ventillation factor, D89 p3103 */
      f1s(       is) =     0.65;     /* Coefficient in the ventillation factor, D89 p3103 */
      f2s(       is) =     0.44;     /* Coefficient in the ventillation factor, D89 p3103 */
      P_EXP_LIQ( is) =     0.319;    /* Exponent in the (p0/p)^x term for rain: liquid cloud particles do not fall at this point */
      P_EXP_ICE( is) =     0.202;    /* Exponent in the (p0/p)^x term for cloud ice */
      P_EXP_SNOW(is) =     0.298;    /* Exponent in the (p0/p)^x term for snow. */
      ALPHA_RAUT(is) =     0.001;    /* Rate coeff. for autoconversion of cloud water [sec-1]. Fowler et al (1996) and Kessler (1969) */
      Q_LIQ_0(   is) =     0.00025;  /* Threshold value for rain autoconversion [kg kg-1]. Fowler et al. (1996) */
      T_00(      is) = T_triple_pt(is);

      is = H_2O_INDEX;
      RHO_RAIN(  is) =  1000.;        /* Fowler et al. (1996) */
      RHO_SNOW(  is) =   917.0*0.5;   /* Perry's Chem.Eng.Handbook p.2-304. */
      COEFF_XI(  is) =   395.78; 
      COEFF_YI(  is) =     0.88525;
      COEFF_XS(  is) =    33.184; 
      COEFF_YS(  is) =     0.48591; 
      COEFF_XR(  is) =  1734.8; 
      COEFF_YR(  is) =    0.74991; 
      COEFF_A(   is) =   13.16;
      COEFF_B(   is) =    0.16;
      COEFF_C(   is) =    5.38e+7;
      COEFF_D(   is) =    0.75;
      COEFF_M(   is) =    7.06165e-3; /* HDC04 */ 
      COEFF_N(   is) =    2.0;
      N_0R(      is) =    8.e+6;      /* Intercept value in raindrop size distribution [m-4], FRR96 p526, D89 p3103 */
      E_R(       is) =    1.0;        /* Accretion efficiency for rain, D89 p3104 */
      D_Icrit(   is) =  500.e-6;      /* Maximum diameter of cloud ice crystal H04 p108 */
      f1r(       is) =    0.78;       /* Coefficient in the ventillation factor, D89 p3103 */
      f2r(       is) =    0.308;      /* Coefficient in the ventillation factor, D89 p3103 */
      f1s(       is) =    0.65;       /* Coefficient in the ventillation factor, D89 p3103 */
      f2s(       is) =    0.44;       /* Coefficient in the ventillation factor, D89 p3103 */
      P_EXP_LIQ( is) =    0.33;       /* Exponent in the (p0/p)^x term for rain: liquid cloud particles do not fall at this point */
      P_EXP_ICE( is) =    0.22;       /* Exponent in the (p0/p)^x term for cloud ice */
      P_EXP_SNOW(is) =    0.313;      /* Exponent in the (p0/p)^x term for snow */
      ALPHA_RAUT(is) =    0.001;      /* Rate coeff. for autoconversion of cloud water [sec-1]. Fowler et al (1996) and Kessler (1969) */
      Q_LIQ_0(   is) =    0.00025;    /* Threshold value for rain autoconversion [kg kg-1]. Fowler et al. (1996) */
      T_00(      is) =  273.16;
    break;

    case TITAN_INDEX:
      /*
       * Placeholder for Titan's clouds.
       */
    break;

    case URANUS_INDEX:
      is = CH_4_INDEX;
      RHO_RAIN(  is) =   448.9; 
      RHO_SNOW(  is) =   505.0*0.5;
      COEFF_XI(  is) =   179.6;
      COEFF_YI(  is) =     0.87644;
      COEFF_XS(  is) =    14.145;
      COEFF_YS(  is) =     0.47194; 
      COEFF_XR(  is) =  707.94;
      COEFF_YR(  is) =     0.74922;
      COEFF_A(   is) =     16.45;
      COEFF_B(   is) =     0.16;
      COEFF_C(   is) =     5.38e+7; 
      COEFF_D(   is) =     0.75;
      COEFF_M(   is) =     7.06165e-3*RHO_SNOW(is)/(0.5*917.0);
      COEFF_N(   is) =     2.0;
      N_0R(      is) =     8.e+6;     /* Intercept value in raindrop size distribution [m-4], FRR96 p526, D89 p3103 */
      E_R(       is) =     1.0;       /* Accretion efficiency for rain, D89 p3104 */
      D_Icrit(   is) =   500.e-6;     /* Maximum diameter of cloud ice crystal H04 p108 */
      f1r(       is) =     0.78;      /* Coefficient in the ventillation factor, D89 p3103 */
      f2r(       is) =     0.308;     /* Coefficient in the ventillation factor, D89 p3103 */
      f1s(       is) =     0.65;      /* Coefficient in the ventillation factor, D89 p3103 */
      f2s(       is) =     0.44;      /* Coefficient in the ventillation factor, D89 p3103 */
      P_EXP_LIQ( is) =     0.321;     /* Exponent in the (p0/p)^x term for rain: liquid cloud particles do not fall at this point */ 
      P_EXP_ICE( is) =     0.268;     /* Exponent in the (p0/p)^x term for cloud ice */ 
      P_EXP_SNOW(is) =     0.302;     /* Exponent in the (p0/p)^x term for snow */
      ALPHA_RAUT(is) =     0.001;     /* Rate coeff. for autoconversion of cloud water [sec-1]. Fowler et al (1996) and Kessler (1969) */
      Q_LIQ_0(   is) =     0.00025;   /* Threshold value for rain autoconversion [kg kg-1]. Fowler et al. (1996) */
      T_00(      is) = T_triple_pt(is);
    break;

    case NEPTUNE_INDEX:
      is = CH_4_INDEX;
      RHO_RAIN(  is) =   448.9; 
      RHO_SNOW(  is) =   505.0*0.5;
      COEFF_XI(  is) =   203.99;
      COEFF_YI(  is) =     0.85989;
      COEFF_XS(  is) =    14.684;
      COEFF_YS(  is) =     0.44883; 
      COEFF_XR(  is) =  727.71;
      COEFF_YR(  is) =     0.72708;
      COEFF_A(   is) =     16.45;
      COEFF_B(   is) =     0.16;
      COEFF_C(   is) =     5.38e+7; 
      COEFF_D(   is) =     0.75;
      COEFF_M(   is) =     7.06165e-3*RHO_SNOW(is)/(0.5*917.0);
      COEFF_N(   is) =     2.0;
      N_0R(      is) =     8.e+6;     /* Intercept value in raindrop size distribution [m-4], FRR96 p526, D89 p3103 */
      E_R(       is) =     1.0;       /* Accretion efficiency for rain, D89 p3104 */
      D_Icrit(   is) =   500.e-6;     /* Maximum diameter of cloud ice crystal H04 p108 */
      f1r(       is) =     0.78;      /* Coefficient in the ventillation factor, D89 p3103 */
      f2r(       is) =     0.308;     /* Coefficient in the ventillation factor, D89 p3103 */
      f1s(       is) =     0.65;      /* Coefficient in the ventillation factor, D89 p3103 */
      f2s(       is) =     0.44;      /* Coefficient in the ventillation factor, D89 p3103 */
      P_EXP_LIQ( is) =     0.306;     /* Exponent in the (p0/p)^x term for rain: liquid cloud particles do not fall at this point */ 
      P_EXP_ICE( is) =     0.247;     /* Exponent in the (p0/p)^x term for cloud ice */ 
      P_EXP_SNOW(is) =     0.312;     /* Exponent in the (p0/p)^x term for snow */
      ALPHA_RAUT(is) =     0.001;     /* Rate coeff. for autoconversion of cloud water [sec-1]. Fowler et al (1996) and Kessler (1969) */
      Q_LIQ_0(   is) =     0.00025;   /* Threshold value for rain autoconversion [kg kg-1]. Fowler et al. (1996) */
      T_00(      is) = T_triple_pt(is);
    break;

    default:
      sprintf(Message,"not yet implemented for planet_index=%d",planet_index);
      epic_error(dbmsname,Message);
    break;
  }

  return;
}

/*======================= end of set_microphysics_params() ======================*/

/*======================= read_enthalpy_change_data() ===========================*/

/*
 * Allocate memory for data tables and read in values.
 */

int read_enthalpy_change_data(int             species_index,
                              float_triplet **pt_hi,
                              float_triplet **pt_hf,
                              float_triplet **pt_hg)
{
  int
    j,
    ndat;
  char
    infile[N_STR],
    header[N_STR];
  float_triplet
    *hi,
    *hf,
    *hg;
#if defined(EPIC_MPI)
  EPIC_FLOAT
    *buffer;
#endif
  FILE
    *enth_vs_t;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_enthalpy_change_data";

  if (IAMNODE == NODE0) {
    sprintf(infile,EPIC_PATH"/data/chemistry/enthalpy/enth_vs_t.%s",
                   var.species[species_index].info[0].name);
    enth_vs_t = fopen(infile,"r");
    if (!enth_vs_t) {
      sprintf(Message,"failed to open %s",infile);
      epic_error(dbmsname,Message);
    }
    /* Skip over header: */
    for (j = 0; j < 6; j++) {
      fgets(header,100,enth_vs_t);
    }
    /* Input number of data points: */
    fscanf(enth_vs_t,"%d",&ndat);
  }

#if defined(EPIC_MPI)
  MPI_Bcast(&ndat,1,MPI_INT,NODE0,para.comm);
#endif

    /*
     * Allocate memory.
     *
     * NOTE: This function should only be called once per species.
     */
    *pt_hi = hi = ftriplet(0,ndat-1,dbmsname);
    *pt_hf = hf = ftriplet(0,ndat-1,dbmsname);
    *pt_hg = hg = ftriplet(0,ndat-1,dbmsname);

    if (IAMNODE == NODE0) {
      /* Input enthalpies(T): */
      for (j = 0; j < ndat; j++) {
#       if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(enth_vs_t,"%lf %lf %lf %lf",&hi[j].x,&hi[j].y,&hf[j].y,&hg[j].y);
#       else
          fscanf(enth_vs_t,"%f %f %f %f",&hi[j].x,&hi[j].y,&hf[j].y,&hg[j].y);
#       endif
        /* Fill in temperature column for hf and hg. */
        hf[j].x = hi[j].x;
        hg[j].x = hi[j].x;
      }
      fclose(enth_vs_t);
    }

#   if defined(EPIC_MPI)
      /*
       * Pack buffer.
       */
      buffer = fvector(0,4*ndat-1,dbmsname);
      for (j = 0; j < ndat; j++) {
        buffer[j       ] = hi[j].x;
        buffer[j+1*ndat] = hi[j].y;
        buffer[j+2*ndat] = hf[j].y;
        buffer[j+3*ndat] = hg[j].y;
      }
#     if EPIC_PRECISION == DOUBLE_PRECISION
         MPI_Bcast(buffer,4*ndat,MPI_DOUBLE,NODE0,para.comm);
#     else
         MPI_Bcast(buffer,4*ndat,MPI_FLOAT,NODE0,para.comm);
#     endif
      /*
       * Unpack buffer.
       */
      for (j = 0; j < ndat; j++) {
        hi[j].x = hf[j].x = hg[j].x = buffer[j];
        hi[j].y = buffer[j+1*ndat];
        hf[j].y = buffer[j+2*ndat];
        hg[j].y = buffer[j+3*ndat];
      }
      spline_pchip(ndat,hi);
      spline_pchip(ndat,hf);
      spline_pchip(ndat,hg);

      free_fvector(buffer,0,4*ndat-1,dbmsname);
#   endif

  return ndat;
}

/*======================= end of read_enthalpy_change_data() ====================*/

/*======================= enthalpy_change() =====================================*/

/*
 * Generic enthalpy change calculation between two states of the given species.
 * C.J.Palotai *A* 06.2003.
 *
 * The enthalpy data tables are in kJ/kg, which is converted here to J/kg.
 */
EPIC_FLOAT enthalpy_change(int            species_index,
                           int            init_phase,
			   int            final_phase,
			   EPIC_FLOAT     temperature,
                           int            ndat,
                           float_triplet *hi,
                           float_triplet *hf,
                           float_triplet *hg)
{
  int
    j;
  static int
    warned_once=FALSE;
  EPIC_FLOAT
    enth_change,t_d,
    enth_init,enth_final;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change";

  /* Check whether temperature is out of range? */
  if (temperature < hg[0].x) {
    temperature = hg[0].x;
    if (!warned_once) {
      sprintf(Message,"%s, init temperature=%g < hg[0].x=%g",
                      var.species[species_index].info[0].name,temperature,hg[0].x);
      epic_warning(dbmsname,Message);
      warned_once = TRUE;
    }
  }
  else if (temperature > hg[ndat-1].x) {
    temperature = hg[ndat-1].x;
    if (!warned_once) {
      sprintf(Message,"%s, init temperature=%g > hg[ndat-1].x=%g",
                      var.species[species_index].info[0].name,temperature,hg[ndat-1].x);
      epic_warning(dbmsname,Message);
      warned_once = TRUE;
    }
  }
  
  /* ======= Interpolate: ======= */
  switch(init_phase) {
    case SOLID:
      j         = find_place_in_table(ndat,hi,temperature,&t_d);
      enth_init = splint_pchip(temperature,hi+j,t_d);
    break;
    case LIQUID:
      j         = find_place_in_table(ndat,hf,temperature,&t_d);
      enth_init = splint_pchip(temperature,hf+j,t_d);
    break;
    case VAPOR:
      j         = find_place_in_table(ndat,hg,temperature,&t_d);
      enth_init = splint_pchip(temperature,hg+j,t_d);
    break;
    default:
      sprintf(Message,"init_phase=%d not recognized",init_phase);
      epic_error(dbmsname,Message);
    break;
  }

  switch (final_phase) {
    case SOLID:
      j          = find_place_in_table(ndat,hi,temperature,&t_d);
      enth_final = splint_pchip(temperature,hi+j,t_d);
    break;
    case LIQUID:
      j          = find_place_in_table(ndat,hf,temperature,&t_d);
      enth_final = splint_pchip(temperature,hf+j,t_d);
    break;
    case VAPOR:
      j          = find_place_in_table(ndat,hg,temperature,&t_d);
      enth_final = splint_pchip(temperature,hg+j,t_d);
    break;
    default:
      sprintf(Message,"final_phase=%d not recognized",final_phase);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * Check for NaN.
   */
  if (!isfinite(enth_init)) {
    sprintf(Message,"species=%s, temperature=%g, init_phase=%d enth_init=%g",
                     var.species[species_index].info[0].name,temperature,init_phase,enth_init);
    epic_error(dbmsname,Message);
  }
  if (!isfinite(enth_final)) {
    sprintf(Message,"species=%s, final_phase=%d, temperature=%g enth_final=%g",
                     var.species[species_index].info[0].name,final_phase,temperature,enth_final);
    epic_error(dbmsname,Message);
  }
  enth_change = enth_final-enth_init;

  /* Convert from kJ/kg to J/kg: */
  enth_change *= 1.e+3;

  return enth_change;
}

/*======================= end of enthalpy_change() ==============================*/

/*======================= enthalpy_change_H_2O() ===============================*/

/*
 * Returns the enthalpy change of water vapor between two temperatures, in J/kg.
 * C.J.Palotai *A* 03.2003.
 */

EPIC_FLOAT enthalpy_change_H_2O(int         init_phase,
                                int         final_phase, 
	                        EPIC_FLOAT  temperature)

{
  static int
    species_index = H_2O_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_H_2O";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);
    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_H_2O ==========================*/

/*======================= enthalpy_change_NH_3() ===============================*/

/*
 * Returns the enthalpy change of ammonia vapor between two temperatures, in J/kg.
 * C.J.Palotai *A* 03.2003.
 */

EPIC_FLOAT enthalpy_change_NH_3(int         init_phase,
                                int         final_phase, 
	                        EPIC_FLOAT  temperature)
{
  static int
    species_index = NH_3_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_NH_3";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_NH3 ============================*/

/*======================= enthalpy_change_H_2S() ===============================*/

/*
 * Returns the enthalpy change of hydrogen sulfide vapor between two temperatures, in J/kg.
 * C.J.Palotai *A* 03.2003.
 */

EPIC_FLOAT enthalpy_change_H_2S(int         init_phase,
                                int         final_phase, 
                                EPIC_FLOAT  temperature)
{
  static int
    species_index = H_2S_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_H_2S";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_H_2S ===========================*/

/*======================= enthalpy_change_CH_4() ===============================*/

/*
 * Returns the enthalpy change of methane vapor between two temperatures, in J/kg.
 * C.J.Palotai *A* 03.2003.
 */

EPIC_FLOAT enthalpy_change_CH_4(int         init_phase,
                                int         final_phase, 
                                EPIC_FLOAT  temperature)
{
  static int
    species_index = CH_4_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_CH_4";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_CH_4 ===========================*/

/*======================= enthalpy_change_C_2H_2() ===============================*/

/*
 * Returns the enthalpy change of acetylene vapor between two temperatures, in J/kg.
 * NOTE: Need to add data.
 */

EPIC_FLOAT enthalpy_change_C_2H_2(int         init_phase,
                                  int         final_phase, 
                                  EPIC_FLOAT  temperature)
{
  static int
    species_index = C_2H_2_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_C_2H_2";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_C_2H_2 =========================*/

/*======================= enthalpy_change_C_2H_4() ===============================*/

/*
 * Returns the enthalpy change of ethylene vapor between two temperatures, in J/kg.
 * NOTE: Need to add data.
 */

EPIC_FLOAT enthalpy_change_C_2H_4(int         init_phase,
                                  int         final_phase, 
                                  EPIC_FLOAT  temperature)
{
  static int
    species_index = C_2H_4_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_C_2H_4";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_C_2H_2 =========================*/

/*======================= enthalpy_change_C_2H_6() ===============================*/

/*
 * Returns the enthalpy change of ethane vapor between two temperatures, in J/kg.
 * NOTE: Need to add data.
 */

EPIC_FLOAT enthalpy_change_C_2H_6(int         init_phase,
                                  int         final_phase, 
                                  EPIC_FLOAT  temperature)
{
  static int
    species_index = C_2H_6_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_C_2H_6";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_C_2H_6 =========================*/

/*======================= enthalpy_change_CO_2() ================================*/

/*
 * Returns the enthalpy change of carbon dioxide vapor between two temperatures, in J/kg.
 * C.J.Palotai *A* 03.2003.
 */

EPIC_FLOAT enthalpy_change_CO_2(int         init_phase,
                                int         final_phase, 
	                        EPIC_FLOAT  temperature)
{
  static int
    species_index = CO_2_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_CO_2";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_CO_2 ===========================*/

/*======================= enthalpy_change_NH_4SH() ==============================*/

/*
 * Returns the enthalpy change of ammonium hydrosulfide vapor between two temperatures, in J/kg.
 * C.J.Palotai *A* 03.2003.
 */

EPIC_FLOAT enthalpy_change_NH_4SH(int         init_phase,
                                  int         final_phase, 
                                  EPIC_FLOAT  temperature)
{
  static int
    species_index = NH_4SH_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_NH_4SH";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_NH_4SH =========================*/

/*======================= enthalpy_change_O_3() =================================*/

/*
 * Returns the enthalpy change of ozone vapor between two temperatures, in J/kg.
 * C.J.Palotai *A* 03.2003.
 */

EPIC_FLOAT enthalpy_change_O_3(int         init_phase,
                               int         final_phase, 
	                       EPIC_FLOAT  temperature)
{
  static int
    species_index = O_3_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_O_3";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_O_3 ============================*/

/*======================= enthalpy_change_N_2() =================================*/

/*
 * Returns the enthalpy change of nitrogen vapor between two temperatures, in J/kg.
 * C.J.Palotai *A* 03.2003.
 */

EPIC_FLOAT enthalpy_change_N_2(int         init_phase,
                               int         final_phase, 
                               EPIC_FLOAT  temperature)
{
  static int
    species_index = N_2_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_N_2";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_N_2 ============================*/

/*======================= enthalpy_change_PH_3() ================================*/

/*
 * Returns the enthalpy change of phosphine vapor between two temperatures, in J/kg.
 * Need data for this function.
 */

EPIC_FLOAT enthalpy_change_PH_3(int         init_phase,
                                int         final_phase, 
                                EPIC_FLOAT  temperature)
{
  static int
    species_index = PH_3_INDEX,
    ndat,
    initialized = FALSE;
  static float_triplet
    *hi,
    *hf,
    *hg;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy_change_PH_3";

  /*
   * Initialization:
   */
  if (!initialized) {
    ndat = read_enthalpy_change_data(species_index,&hi,&hf,&hg);

    initialized = TRUE;
  }
  /* End of initialization. */

  return enthalpy_change(species_index,init_phase,final_phase,temperature,ndat,hi,hf,hg);
}

/*======================= end of enthalpy_change_PH_3 ===========================*/

/*======================= sat_vapor_p_H_2O() ====================================*/

/*
 * From Table 1.20 of Lodders and Fegley, 1998, 
 * The Planetary Scientist's Companion, Oxford University Press.
 * The parameter a has been modified to yield pressure in Pa instead of bar.
 * Temperature in Kelvin.
 *
 * 06/2003 CJP *A*
 * Added case for mixed phase: Supercooled liquid data taken from:
 * Rogers&Yau: A Short Course in Cloud Physics,3rd Ed., p16. 
 */

EPIC_FLOAT sat_vapor_p_H_2O(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    t_c,omega,
    sat_vapor_p;
  static EPIC_FLOAT
    a_s =    12.610,
    b_s = -2681.18,
    a_l =    11.079,
    b_l = -2261.10;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_H_2O";

  if (temperature <= T_triple_pt(H_2O_INDEX)) {
    sat_vapor_p = pow(10.,a_s+b_s/temperature);
  }
  else if (temperature <  T_triple_pt(H_2O_INDEX) &&
           temperature >  T_00(H_2O_INDEX)) {
    omega       = (temperature-T_triple_pt(H_2O_INDEX))/
                  (T_triple_pt(H_2O_INDEX)-T_00(H_2O_INDEX));
    t_c         = temperature-273.15;
    t_c         = 17.67*t_c/(t_c+243.5);
    sat_vapor_p = 100*omega*6.112*pow(10.,t_c)+(1.-omega)*pow(10.,a_l+b_l/temperature);		 
  }	   
  else {
    sat_vapor_p = pow(10.,a_l+b_l/temperature);
  }

  return sat_vapor_p;
}

/*======================= end of sat_vapor_p_H_2O() =============================*/

/*======================= sat_vapor_p_NH_3() ====================================*/

/*
 * From Table 1.20 of Lodders and Fegley, 1998, 
 * The Planetary Scientist's Companion, Oxford University Press.
 * The parameter a has been modified to yield pressure in Pa instead of bar.
 * Temperature in Kelvin.
 */

EPIC_FLOAT sat_vapor_p_NH_3(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  static EPIC_FLOAT
    a_s =    11.900,
    b_s = -1588.,
    a_l =    10.201,
    b_l = -1248.;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_NH_3";

  if (temperature < T_00(NH_3_INDEX)) {
    sat_vapor_p = pow(10.,a_s+b_s/temperature);
  }
  else {
    sat_vapor_p = pow(10.,a_l+b_l/temperature);
  }

  return sat_vapor_p;
}

/*======================= end of sat_vapor_p_NH_3() =============================*/

/*======================= sat_vapor_p_H_2S() ====================================*/

/*
 * From Table 1.20 of Lodders and Fegley, 1998, 
 * The Planetary Scientist's Companion, Oxford University Press.
 * The parameter a has been modified to yield pressure in Pa instead of bar.
 * Temperature in Kelvin.
 */

EPIC_FLOAT sat_vapor_p_H_2S(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  static EPIC_FLOAT
    a_s =    10.610,
    b_s = -1171.2,
    a_l =     9.780,
    b_l = -1015.5;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_H_2S";

  if (temperature < T_triple_pt(H_2S_INDEX)) {
    sat_vapor_p = pow(10.,a_s+b_s/temperature);
  }
  else {
    sat_vapor_p = pow(10.,a_l+b_l/temperature);
  }

  return sat_vapor_p;
}

/*======================= end of sat_vapor_p_H_2S() =============================*/

/*======================= sat_vapor_p_CH_4() ====================================*/

/*
 * From Table 1.20 of Lodders and Fegley, 1998, 
 * The Planetary Scientist's Companion, Oxford University Press.
 * The parameter a has been modified to yield pressure in Pa instead of bar.
 * Temperature in Kelvin.
 */

EPIC_FLOAT sat_vapor_p_CH_4(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  static EPIC_FLOAT
    a_s =     9.283,
    b_s =  -475.6,
    a_l =     9.092,
    b_l =  -459.8;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_CH_4";

  if (temperature < T_triple_pt(CH_4_INDEX)) {
    sat_vapor_p = pow(10.,a_s+b_s/temperature);
  }
  else {
    sat_vapor_p = pow(10.,a_l+b_l/temperature);
  }

  return sat_vapor_p;
}

/*======================= end of sat_vapor_p_CH_4() =============================*/

/*======================= sat_vapor_p_C_2H_2() ==================================*/

/*
 * NOTE: Need this function.
 */

EPIC_FLOAT sat_vapor_p_C_2H_2(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_C_2H_2";

  /* Placeholder */
  sat_vapor_p = 0.;

  return sat_vapor_p;
}

/*======================= end of sat_vapor_p_C_2H_2() ===========================*/

/*======================= sat_vapor_p_C_2H_4() ==================================*/

/*
 * NOTE: Need this function.
 */

EPIC_FLOAT sat_vapor_p_C_2H_4(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_C_2H_4";

  /* Placeholder */
  sat_vapor_p = 0.;

  return sat_vapor_p;
}

/*======================= end of sat_vapor_p_C_2H_4() ===========================*/

/*======================= sat_vapor_p_C_2H_6() ==================================*/

/*
 * NOTE: Need this function.
 */

EPIC_FLOAT sat_vapor_p_C_2H_6(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_C_2H_6";

  /* Placeholder */
  sat_vapor_p = 0.;

  return sat_vapor_p;
}

/*======================= end of sat_vapor_p_C_2H_6() ===========================*/

/*======================= sat_vapor_p_CO_2() ====================================*/

/*
 * From Table 1.20 of Lodders and Fegley, 1998, 
 * The Planetary Scientist's Companion, Oxford University Press.
 * The parameter a has been modified to yield pressure in Pa instead of bar.
 * Temperature in Kelvin.
 */

EPIC_FLOAT sat_vapor_p_CO_2(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  static EPIC_FLOAT
    a_s =     12.025,
    b_s =  -1336.,
    a_l =     11.045,
    b_l =  -1201.;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_CO_2";

  if (temperature < T_triple_pt(CO_2_INDEX)) {
    sat_vapor_p = pow(10.,a_s+b_s/temperature);
  }
  else {
    sat_vapor_p = pow(10.,a_l+b_l/temperature);
  }

  return sat_vapor_p;
}

/*======================= end of sat_vapor_p_CO_2() =============================*/

/*======================= sat_vapor_p_NH_4SH() ==================================*/

/*
 * From Table 1.20 of Lodders and Fegley, 1998, 
 * The Planetary Scientist's Companion, Oxford University Press.
 * The parameter a has been modified to yield pressure in Pa instead of bar.
 * Temperature in Kelvin.
 */

EPIC_FLOAT sat_vapor_p_NH_4SH(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  static EPIC_FLOAT
    a_s =     12.60,
    b_s =  -2411.2,
    a_l =     12.60,   /* No liquid value listed; here a_l = a_s. */
    b_l =  -2411.2;    /* No liquid value listed; here b_l = b_s. */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_NH_4SH";

  if (temperature < T_triple_pt(NH_4SH_INDEX)) {
    sat_vapor_p = pow(10.,a_s+b_s/temperature);
  }
  else {
    sat_vapor_p = pow(10.,a_l+b_l/temperature);
  }

  return sat_vapor_p;
}

/*======================= end of sat_vapor_p_NH_4SH() ===========================*/

/*======================= sat_vapor_p_O_3() =====================================*/

/*
 * From Table 1.20 of Lodders and Fegley, 1998, 
 * The Planetary Scientist's Companion, Oxford University Press.
 * The parameter a has been modified to yield pressure in Pa instead of bar.
 * Temperature in Kelvin.
 */

EPIC_FLOAT sat_vapor_p_O_3(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  static EPIC_FLOAT
    a_s =     8.912, /* No solid value listed; here a_s = a_l. */
    b_s =  -632.4,   /* No solid value listed; here b_s = b_l. */
    a_l =     8.912,
    b_l =  -632.4;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_O_3";

  if (temperature < T_triple_pt(O_3_INDEX)) {
    sat_vapor_p = pow(10.,a_s+b_s/temperature);
  }
  else {
    sat_vapor_p = pow(10.,a_l+b_l/temperature);
  }

  return sat_vapor_p;
}

/*======================= end of sat_vapor_p_O_3() ==============================*/

/*======================= sat_vapor_p_N_2() =====================================*/

/*
 * From Table 1.20 of Lodders and Fegley, 1998, 
 * The Planetary Scientist's Companion, Oxford University Press.
 * The parameter a has been modified to yield pressure in Pa instead of bar.
 * Temperature in Kelvin.
 */

EPIC_FLOAT sat_vapor_p_N_2(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  static EPIC_FLOAT
    a_s =     9.798,
    b_s =  -360.2,
    a_l =     8.944,
    b_l =  -305.;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_N_2";

  if (temperature < T_triple_pt(N_2_INDEX)) {
    sat_vapor_p = pow(10.,a_s+b_s/temperature);
  }
  else {
    sat_vapor_p = pow(10.,a_l+b_l/temperature);
  }

  return sat_vapor_p;
}

/*======================== end of sat_vapor_p_N_2() ==============================*/

/*======================= sat_vapor_p_PH_3() =====================================*/

/*
 * Need data for this function.
 */

EPIC_FLOAT sat_vapor_p_PH_3(EPIC_FLOAT temperature)
{
  EPIC_FLOAT
    sat_vapor_p;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="sat_vapor_p_PH_3";

  /* Need to write this function */

  return sat_vapor_p;
}

/*======================== end of sat_vapor_p_N_2() ==============================*/

/*=============================== dynvisc() ======================================*/

/*
 * Code for calculating dynamic viscosity
 *           CJP *A*   4/12/2004
 *        Jupiter: C.F. Hansen (1979)
 */

EPIC_FLOAT dynvisc(      char   *globe,
                   EPIC_FLOAT    temp)
{
  EPIC_FLOAT 
    x1, x2,n1,n2,q1,q2,q3,r1,r2,nu;
  static EPIC_FLOAT 
    viscx[65],viscy[65],visca[2];
  int 
    i,
    n = 64;
  EPIC_FLOAT 
    temperature;
  static int
    initialized=0;
 /*
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="dynvisc";
    
  
  if (strcmp(globe,"Jupiter") == 0 ) {
    if (!initialized) {
      initialized = TRUE;
  
      /*
       * Note: Unit conversion between micropoise and kg/m/s : 1 micropoise = 1.e-7 kg/m/s
       */

      x1 = planet->x_h2;                    /* Mol fraction of H_2   */
      x2 = planet->x_he;                    /* Mol fraction of He    */

      for (i=0; i<=n; i++) {   
        
	temperature = 100. + i*(500.-100.)/n;
	
	n1  = 90.6*pow(temperature/300.,.6658);
        n1  = n1/(1.+4./temperature);                 /*  viscosity of H_2 in micropoise */
        n2  = 191.6*pow(temperature/300.,.7176);
        n2  = n2/(1.-11.4/temperature);               /*  viscosity of He in micropoise */

        q1  = 32.3*(1.+4./temperature)*pow(300./temperature,.1658);
        q2  = 21.5*(1.-11.4/temperature)*pow(300./temperature,.2176);
        q3  = (sqrt(q1)+sqrt(q2))/2.;
        q3  = q3*q3;

        r1  = 1.+.7967*(x2/x1)*(q3/q1);
        r2  = 1.+.5634*(x1/x2)*(q3/q2);
        nu  = n1/r1 + n2/r2 ;                  /*  viscosity of atmosphere in micropoise */
        nu  = nu*1.e-7;                         /*  viscosity of atmosphere in kg/m/s     */
      
        viscx[i] = temperature;
	viscy[i] = nu;
      }
      least_squares(viscx,viscy,n,visca);
    }  
    /* End of initialization. */
    
    nu = visca[0] + visca[1]*temp;   

  }  /* End of Jupiter case */

  else {
    nu = planet->dynvisc;
  }

 
  return(nu);
}


/*============================== end of dynvisc() ================================*/

/*=============================== conductivity() ======================================*/

/*
 * Code for calculating thermal coductivity
 *           CJP *A*   4/12/2004
 *        Jupiter: C.F. Hansen (1979)
 */

EPIC_FLOAT conductivity(      char   *globe,
                        EPIC_FLOAT    temp)
{
  EPIC_FLOAT 
    x1, x2,n1,n2,q1,q2,q3,r1,r2,k,k1,k2,v;
  static EPIC_FLOAT 
    condx[65],condy[65],conda[2];
  int 
    i,
    n = 64;
  EPIC_FLOAT 
    temperature;
  static int
    initialized=0;
 /*
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="conductivity";
    
  
  if (strcmp(globe,"Jupiter") == 0 ) {
    if (!initialized) {
      initialized = TRUE;
  
      /*
       * Note: Unit conversion between micropoise and kg/m/s : 1 micropoise = 1.e-7 kg/m/s
       */

      x1 = planet->x_h2;                    /* Mol fraction of H_2   */
      x2 = planet->x_he;                    /* Mol fraction of He    */

      for (i=0; i<=n; i++) {   
        
	temperature = 100. + i*(500.-100.)/n;
	
        k1 = .11*pow(temperature/300.,.6983)/(1.+49.4/temperature);
        v  = 3079.5/temperature;
        v  = 2.*v/(exp(v)-exp(-v));
        k1 = k1*(4.75+v*v);                                 /*  conductivity of H_2 in millical/cm/deg/sec */
        k2 = .3418*pow(temperature/300.,.7412)/(1.-13.74/temperature); /*  conductivity of He in millical/cm/deg/sec */
        q1 = 26.1*(1.+49.4/temperature)*pow(300./temperature,.1983);
        q2 = 5.96*(1.-13.74/temperature)*pow(300./temperature,.2412);
        q3 = (sqrt(q1)+sqrt(q2))/2.;
        q3 = q3*q3;
        r1 = 1.+0.7698*(x2/x1)*(q3/q1);
        r2 = 1.+1.0887*(x1/x2)*(q3/q2);
        k  = k1/r1 + k2/r2 ;              /*  conductivity of atmosphere in millical/cm/deg/sec */
        k  = .4186 * k;                   /*  conductivity of atmosphere in J/m/K/s             */
      
        condx[i] = temperature;
	condy[i] = k;
      }
      least_squares(condx,condy,n,conda);

    }  
    /* End of initialization. */
    
    k = conda[0] + conda[1]*temp;   

  }  /* End of Jupiter case */

  else {
    k = K_a;         /* K_a defined in epic_microphysics.h */
                     /* grub : might be better to incorporate it into the planet-structure */
  }

 
  return(k);
}


/*=========================== end of conductivity() =============================*/

/* * * * * * * * * * * * * * * end of epic_microphysics_funcs.c * * * * * * * * * */
