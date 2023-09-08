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

/* * * * * * * * *  epic_funcs_diag.c  * * * * * * * * * * * * * * * * * * * *
 *                                                                           *
 *       Timothy E. Dowling                                                  *
 *                                                                           *
 *       Functions that reference data or EPIC model variables directly      *
 *       by name or index.                                                   *
 *                                                                           *
 *       This file includes the following:                                   *
 *                                                                           *
 *           set_var_props()                                                 *
 *           free_var_props()                                                *
 *           make_arrays()                                                   *
 *           free_arrays()                                                   *
 *           return_sigmatheta()                                             *
 *           f_sigma()                                                       *
 *           g_sigma()                                                       *
 *           set_lonlat()                                                    *
 *           set_fmn()                                                       *
 *           set_gravity()                                                   *
 *           set_dsgth()                                                     *
 *           get_sigma()                                                     *
 *           get_p_sigma()                                                   *
 *           calc_h()                                                        *
 *           get_p()                                                         *
 *           molar_mixing_ratio()                                            *
 *           get_var()                                                       *
 *           get_var_mean2d()                                                *
 *           onto_kk()                                                       *
 *           fpe_minus_fpe()                                                 *
 *           get_kin()                                                       *
 *           get_brunt2()                                                    *
 *           get_richardson()                                                *
 *           set_p2_etc()                                                    *
 *           store_pgrad_vars()                                              *
 *           store_diag()                                                    *
 *           divergence()                                                    *
 *           vorticity()                                                     *
 *           phi_from_u()                                                    *
 *           relative_humidity()                                             *
 *           source_sink()                                                   *
 *           cfl_dt()                                                        *
 *           time_mod()                                                      *
 *           avg_molar_mass()                                                *
 *           sync_x_to_q()                                                   *
 *           molar_mass()                                                    *
 *           diffusivity()                                                   *
 *           timeplane_bookkeeping()                                         *
 *           check_nan()                                                     *
 *           u_venus, u_jupiter(), etc.                                      *
 *           u_amp()                                                         *
 *           pioneer_venus_u()                                               *
 *           galileo_u()                                                     *
 *           cassini_cirs_u()                                                *
 *           p_sigmatheta()                                                  *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*======================= set_var_props() ===================================*/

/*
 * Wind variables: 
 *   3rd-Order Adams-Bashforth timestep:  
 *     Need 3 time planes of tendencies, 
 *     the oldest 2 of which are written to epic.nc files.
 *   Leapfrog timestep:
 *     Need 2 time planes for variable, both of which are 
 *     written to epic.nc files.
 */
#define SET_WIND(iwind,wind,the_standard_name,the_long_name,the_units) \
          if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) { \
            var.wind.info[0].index = iwind; \
            var.wind.info[0].name = (char *)malloc(strlen(#wind)+1); \
            strcpy(var.wind.info[0].name,#wind); \
            var.wind.info[0].standard_name = (char *)malloc(strlen(#the_standard_name)+1); \
            strcpy(var.wind.info[0].standard_name,#the_standard_name); \
            var.wind.info[0].long_name = (char *)malloc(strlen(#the_long_name)+1); \
            strcpy(var.wind.info[0].long_name,#the_long_name); \
            var.wind.info[0].units = (char *)malloc(strlen(#the_units)+1); \
            strcpy(var.wind.info[0].units,#the_units); \
            var.wind.info_tend[0].index = iwind*1000000+1; \
            var.wind.info_tend[0].name = (char *)malloc(strlen("d"#wind"_IT_MINUS1dt")+1); \
            strcpy(var.wind.info_tend[0].name,"d"#wind"_IT_MINUS1dt"); \
            var.wind.info_tend[0].units = (char *)malloc(strlen(#the_units"/s")+1); \
            strcpy(var.wind.info_tend[0].units,#the_units"/s"); \
            var.wind.info_tend[0].standard_name = (char *)malloc(strlen("tendency_of_"#the_standard_name)+1); \
            strcpy(var.wind.info_tend[0].standard_name,"tendency_of_"#the_standard_name); \
            var.wind.info_tend[0].long_name = (char *)malloc(strlen("tendency of "#the_long_name)+1); \
            strcpy(var.wind.info_tend[0].long_name,"tendency of "#the_long_name); \
            var.wind.info_tend[1].index = iwind*1000000+2; \
            var.wind.info_tend[1].name = (char *)malloc(strlen("d"#wind"_IT_MINUS2dt")+1); \
            strcpy(var.wind.info_tend[1].name,"d"#wind"_IT_MINUS2dt"); \
            var.wind.info_tend[1].units = (char *)malloc(strlen(#the_units"/s")+1); \
            strcpy(var.wind.info_tend[1].units,#the_units"/s"); \
            var.wind.info_tend[1].standard_name = (char *)malloc(strlen("tendency_of_"#the_standard_name)+1); \
            strcpy(var.wind.info_tend[1].standard_name,"tendency_of_"#the_standard_name); \
            var.wind.info_tend[1].long_name = (char *)malloc(strlen("tendency of "#the_long_name)+1); \
            strcpy(var.wind.info_tend[1].long_name,"tendency of "#the_long_name); \
          } \
          else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) { \
            var.wind.info[0].index = iwind; \
            var.wind.info[0].name = (char *)malloc(strlen(#wind)+1); \
            strcpy(var.wind.info[0].name,#wind); \
            var.wind.info[0].standard_name = (char *)malloc(strlen(#the_standard_name)+1); \
            strcpy(var.wind.info[0].standard_name,#the_standard_name); \
            var.wind.info[0].long_name = (char *)malloc(strlen(#the_long_name)+1); \
            strcpy(var.wind.info[0].long_name,#the_long_name); \
            var.wind.info[0].units = (char *)malloc(strlen(#the_units)+1); \
            strcpy(var.wind.info[0].units,#the_units); \
            var.wind.info[1].index = iwind*1000000+1; \
            var.wind.info[1].name = (char *)malloc(strlen(#wind"_IT_MINUS1")+1); \
            strcpy(var.wind.info[1].name,#wind"_IT_MINUS1"); \
            var.wind.info[1].standard_name = (char *)malloc(strlen(#the_standard_name)+1); \
            strcpy(var.wind.info[1].standard_name,#the_standard_name); \
            var.wind.info[1].long_name = (char *)malloc(strlen(#the_long_name)+1); \
            strcpy(var.wind.info[1].long_name,#the_long_name); \
            var.wind.info[1].units = (char *)malloc(strlen(#the_units)+1); \
            strcpy(var.wind.info[1].units,#the_units); \
          }

#define SET_THERMO(ithermo,thermo,the_standard_name,the_long_name,the_units) \
          var.thermo.info[0].index = ithermo; \
          var.thermo.info[0].name  = (char *)malloc(strlen(#thermo)+1); \
          strcpy(var.thermo.info[0].name,#thermo); \
          var.thermo.info[0].standard_name = (char *)malloc(strlen(#the_standard_name)+1); \
          strcpy(var.thermo.info[0].standard_name,#the_standard_name); \
          var.thermo.info[0].long_name = (char *)malloc(strlen(#the_long_name)+1); \
          strcpy(var.thermo.info[0].long_name,#the_long_name); \
          var.thermo.info[0].units = (char *)malloc(strlen(#the_units)+1); \
          sprintf(var.thermo.info[0].units,"%s",#the_units);

/*
 * Species variables are distinguished by the fact that separate memory is needed
 * to handle each phase.
 */
#define SET_SPECIES(ispecies,the_species,the_standard_name,the_long_name) \
          var.species[ispecies].info[0].index = ispecies; \
          var.species[ispecies].info[0].name = (char *)malloc(strlen(#the_species)+1); \
          strcpy(var.species[ispecies].info[0].name,#the_species); \
          var.species[ispecies].info[0].standard_name = (char *)malloc(strlen(#the_standard_name)+1); \
          strcpy(var.species[ispecies].info[0].standard_name,#the_standard_name); \
          var.species[ispecies].info[0].long_name = (char *)malloc(strlen(#the_long_name)+1); \
          strcpy(var.species[ispecies].info[0].long_name,#the_long_name); \
          for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) { \
            ; /* Mass mixing ratios */ \
            var.species[ispecies].phase[ip].info[MASS].index = ispecies*1000000+(ip+1)*1000+MASS; \
            var.species[ispecies].phase[ip].info[MASS].name  = (char *)malloc(strlen(#the_species"_")+ \
                                                               strlen(Phase_Name[ip])+1); \
            sprintf(var.species[ispecies].phase[ip].info[MASS].name,"%s_%s",#the_species,Phase_Name[ip]); \
            var.species[ispecies].phase[ip].info[MASS].long_name = (char *)malloc(strlen(#the_long_name" ")+ \
                                                                   strlen(Long_Phase_Name[ip])+strlen(" mass mixing ratio")+1); \
            sprintf(var.species[ispecies].phase[ip].info[MASS].long_name,"%s %s mass mixing ratio",#the_long_name,Long_Phase_Name[ip]); \
            var.species[ispecies].phase[ip].info[MASS].units = (char *)malloc(strlen("kg/kg")+1); \
            strcpy(var.species[ispecies].phase[ip].info[MASS].units,"kg/kg"); \
            ; /* Mole fractions */ \
            var.species[ispecies].phase[ip].info[MOLAR].index = ispecies*1000000+(ip+1)*1000+MOLAR; \
            var.species[ispecies].phase[ip].info[MOLAR].name  = (char *)malloc(strlen(#the_species"_")+ \
                                                                strlen(Phase_Name[ip])+1); \
            sprintf(var.species[ispecies].phase[ip].info[MOLAR].name,"%s_%s",#the_species,Phase_Name[ip]); \
            var.species[ispecies].phase[ip].info[MOLAR].long_name = (char *)malloc(strlen(#the_long_name" ")+ \
                                                                    strlen(Long_Phase_Name[ip])+strlen(" mole fraction")+1); \
            sprintf(var.species[ispecies].phase[ip].info[MOLAR].long_name,"%s %s mole fraction",#the_long_name,Long_Phase_Name[ip]); \
            var.species[ispecies].phase[ip].info[MOLAR].units = (char *)malloc(strlen("kmol/kmol")+1); \
            strcpy(var.species[ispecies].phase[ip].info[MOLAR].units,"kmol/kmol"); \
          }
/*
 * Diagnostic variables are not associated with any tendency or moment data.
 */
#define SET_DIAG(idiag,diag,the_standard_name,the_long_name,the_units) \
            var.diag.info[0].index = idiag; \
            var.diag.info[0].name  = (char *)malloc(strlen(#diag)+1); \
            strcpy(var.diag.info[0].name,#diag); \
            var.diag.info[0].standard_name = (char *)malloc(strlen(#the_standard_name)+1); \
            strcpy(var.diag.info[0].standard_name,#the_standard_name); \
            var.diag.info[0].long_name = (char *)malloc(strlen(#the_long_name)+1); \
            strcpy(var.diag.info[0].long_name,#the_long_name); \
            var.diag.info[0].units = (char *)malloc(strlen(#the_units)+1) ; \
            strcpy(var.diag.info[0].units,#the_units);

void set_var_props(planetspec *planet) 
{
  int
    im,is,ip;
  const char
    *Phase_Name[MAX_NUM_PHASES] 
      = {"vapor","liquid","solid","rain","snow"},
    *Long_Phase_Name[MAX_NUM_PHASES]
      = {"vapor","cloud liquid","cloud ice","rain","snow"};
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_var_props";

  switch(planet->index) {
    case VENUS_INDEX:
      planet->u = u_venus;
    break;
    case EARTH_INDEX:
      planet->u = u_earth;
    break;
    case MARS_INDEX:
      planet->u = u_mars;
    break;
    case JUPITER_INDEX:
      planet->u = u_jupiter;
    break;
    case SATURN_INDEX:
      planet->u = u_saturn;
    break;
    case TITAN_INDEX:
      planet->u = u_titan;
    break;
    case URANUS_INDEX:
      planet->u = u_uranus;
    break;
    case NEPTUNE_INDEX:
      planet->u = u_neptune;
    break;
    case TRITON_INDEX:
      planet->u = u_triton;
    break;
    case PLUTO_INDEX:
      planet->u = u_pluto;
    break;
    case HOT_JUPITER_INDEX:
      planet->u = u_hot_jupiter;
    break;
    case HELD_SUAREZ_INDEX:
      planet->u = u_null;
    break;
    case VENUS_LLR05_INDEX:
      planet->u = u_null;
    break;
    case GOLDMINE_INDEX:
      planet->u = u_null;
    break;
    default:
      sprintf(Message,"unrecognized planet->index=%d",planet->index);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * Solar longitude, L_s.
   */
  var.l_s.info.name = (char *)malloc(strlen("L_s")+1);
  strcpy(var.l_s.info.name,"L_s");
  var.l_s.info.long_name = (char *)malloc(strlen("solar longitude (planetocentric)")+1);
  strcpy(var.l_s.info.long_name,"solar longitude (planetocentric)");
  var.l_s.info.units = (char *)malloc(strlen("degrees")+1);
  strcpy(var.l_s.info.units,"degrees");

  /* 
   * Wind variables.
   */
  SET_WIND(U_INDEX,u,eastward_wind,zonal wind,m/s);
  SET_WIND(V_INDEX,v,northward_wind,meridional wind,m/s);

  /*
   * Thermo variables.
   */
  SET_THERMO(H_INDEX,h,air_hybrid_density,total air hybrid density,kg/m^2/K);
  SET_THERMO(THETA_INDEX,theta,air_potential_temperature,potential temperature,K);
  SET_THERMO(FPARA_INDEX,fpara,mole_fraction_of_para_in_hydrogen,para fraction of molecular hydrogen,kmol/kmol);

  /* 
   * Optional species (Qs).
   */
  SET_SPECIES(H_2O_INDEX,H_2O,water,water);
  var.species[H_2O_INDEX].enthalpy_change        = enthalpy_change_H_2O;
  var.species[H_2O_INDEX].sat_vapor_p            = sat_vapor_p_H_2O;
  var.species[H_2O_INDEX].HITRAN_index           = 1;
  
  SET_SPECIES(NH_3_INDEX,NH_3,ammonia,ammonia);
  var.species[NH_3_INDEX].enthalpy_change        = enthalpy_change_NH_3;
  var.species[NH_3_INDEX].sat_vapor_p            = sat_vapor_p_NH_3;
  var.species[NH_3_INDEX].HITRAN_index           = 11;

  SET_SPECIES(H_2S_INDEX,H_2S,hydrogen_sulfide,hydrogen sulfide);
  var.species[H_2S_INDEX].enthalpy_change        = enthalpy_change_H_2S;
  var.species[H_2S_INDEX].sat_vapor_p            = sat_vapor_p_H_2S;
  var.species[H_2S_INDEX].HITRAN_index           = 31;
  
  SET_SPECIES(CH_4_INDEX,CH_4,methane,methane);
  var.species[CH_4_INDEX].enthalpy_change        = enthalpy_change_CH_4;
  var.species[CH_4_INDEX].sat_vapor_p            = sat_vapor_p_CH_4;
  var.species[CH_4_INDEX].HITRAN_index           = 6;

  SET_SPECIES(C_2H_2_INDEX,C_2H_2,acetylene,acetylene);
  var.species[C_2H_2_INDEX].enthalpy_change      = enthalpy_change_C_2H_2;
  var.species[C_2H_2_INDEX].sat_vapor_p          = sat_vapor_p_C_2H_2;
  var.species[C_2H_2_INDEX].HITRAN_index         = 26;

  SET_SPECIES(C_2H_4_INDEX,C_2H_4,ethylene,ethylene);
  var.species[C_2H_4_INDEX].enthalpy_change     = enthalpy_change_C_2H_4;
  var.species[C_2H_4_INDEX].sat_vapor_p         = sat_vapor_p_C_2H_4;
  var.species[C_2H_4_INDEX].HITRAN_index        = 38;

  SET_SPECIES(C_2H_6_INDEX,C_2H_6,ethane,ethane);
  var.species[C_2H_6_INDEX].enthalpy_change      = enthalpy_change_C_2H_6;
  var.species[C_2H_6_INDEX].sat_vapor_p          = sat_vapor_p_C_2H_6;
  var.species[C_2H_6_INDEX].HITRAN_index         = 27;
  
  SET_SPECIES(CO_2_INDEX,CO_2,carbon_dioxide,carbon dioxide);
  var.species[CO_2_INDEX].enthalpy_change        = enthalpy_change_CO_2;
  var.species[CO_2_INDEX].sat_vapor_p            = sat_vapor_p_CO_2;
  var.species[CO_2_INDEX].HITRAN_index           = 2;
 
  SET_SPECIES(NH_4SH_INDEX,NH_4SH,ammonium_hydrosulfide,ammonium hydrosulfide);
  var.species[NH_4SH_INDEX].enthalpy_change      = enthalpy_change_NH_4SH;
  var.species[NH_4SH_INDEX].sat_vapor_p          = sat_vapor_p_NH_4SH;
  var.species[NH_4SH_INDEX].HITRAN_index         = -1; /* needs a value */
  
  SET_SPECIES(O_3_INDEX,O_3,ozone,ozone);
  var.species[O_3_INDEX].enthalpy_change         = enthalpy_change_O_3;
  var.species[O_3_INDEX].sat_vapor_p             = sat_vapor_p_O_3;
  var.species[O_3_INDEX].HITRAN_index            = 3;
  
  SET_SPECIES(N_2_INDEX,N_2,nitrogen,nitrogen);
  var.species[N_2_INDEX].enthalpy_change         = enthalpy_change_N_2;
  var.species[N_2_INDEX].sat_vapor_p             = sat_vapor_p_N_2;
  var.species[N_2_INDEX].HITRAN_index            = 22;

  SET_SPECIES(PH_3_INDEX,PH_3,phosphine,phosphine);
  var.species[PH_3_INDEX].enthalpy_change        = enthalpy_change_PH_3;
  var.species[PH_3_INDEX].sat_vapor_p            = sat_vapor_p_PH_3;
  var.species[PH_3_INDEX].HITRAN_index           = 28;

  /*
   * Set turbulence-model variables.
   */
  if (strcmp(grid.turbulence_scheme,"on")               == 0 ||
      strcmp(grid.turbulence_scheme,"on_vertical_only") == 0)  {
    SET_THERMO(NU_TURB_INDEX,nu_turb,Spalart_Allmaras_turbulent_viscosity,Spalart-Allmaras turbulent viscosity,m^2/s);
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  /* 
   * Diagnostic variables.
   */
  SET_DIAG(HDRY2_INDEX,hdry2,air_dry_hybrid_density,dry-air hybrid density,kg/m^2/K);
  SET_DIAG(HDRY3_INDEX,hdry3,air_dry_hybrid_density,dry-air hybrid density,kg/m^2/K);
  SET_DIAG(H3_INDEX,h3,air_hybrid_density,total hybrid density,kg/m^2/K);
  SET_DIAG(PDRY3_INDEX,pdry3,air_dry_pressure,dry-air pressure,Pa);
  SET_DIAG(P2_INDEX,p2,air_pressure,total air pressure,Pa);
  SET_DIAG(P3_INDEX,p3,air_pressure,total air pressure,Pa);
  SET_DIAG(THETA2_INDEX,theta2,air_potential_temperature,potential temperature,K);
  SET_DIAG(T2_INDEX,t2,air_temperature,temperature,K);
  SET_DIAG(T3_INDEX,t3,air_temperature,temperature,K);
  SET_DIAG(RHO2_INDEX,rho2,air_density,total density,kg/m^3);
  SET_DIAG(RHO3_INDEX,rho3,air_density,total density,kg/m^3);
  SET_DIAG(EXNER2_INDEX,exner2,exner_function,Exner function,J/kg/K);
  SET_DIAG(EXNER3_INDEX,exner3,exner_function,Exner function,J/kg/K);
  SET_DIAG(FGIBB2_INDEX,fgibb2,gibbs_term_hydrogen,Gibbs term (ortho-para H2),J/kg);
  SET_DIAG(PHI2_INDEX,phi2,geopotential,geopotential,m^2/s^2);
  SET_DIAG(PHI3_INDEX,phi3,geopotential,geopotential,m^2/s^2);
  SET_DIAG(MONT2_INDEX,mont2,montgomery_potential,Montgomery potential,m^2/s^2);
  SET_DIAG(HEAT3_INDEX,heat3,heating_rate,heating rate per mass,W/kg);
  SET_DIAG(PV2_INDEX,pv2,potential_vorticity,potential vorticity,m^2/s K/kg);
  SET_DIAG(EDDY_PV2_INDEX,eddy_pv2,eddy_potential_vorticity,eddy potential vorticity,m^2/s K/kg);
  SET_DIAG(MOLAR_MASS3_INDEX,molar_mass3,mean_molar_mass,mean molar mass,kg/kmol);
  SET_DIAG(RI2_INDEX,ri2,local_richardson_number,local Richardson number,s^2/s^2);
  SET_DIAG(REL_VORT2_INDEX,rel_vort2,atmosphere_relative_vorticity,relative vorticity,1/s);
  SET_DIAG(EDDY_REL_VORT2_INDEX,eddy_rel_vort2,eddy_relative_vorticity,eddy relative vorticity,1/s);
  SET_DIAG(ABS_VORT2_INDEX,abs_vort2,atmosphere_absolute_vorticity,absolute vorticity,1/s);
  SET_DIAG(KIN2_INDEX,kinetic_energy2,specific_kinetic_energy_of_air,horizontal kinetic energy per mass,m^2/s^2);
  SET_DIAG(DIV_UV2_INDEX,div_uv2,divergence_of_wind,horizontal divergence,1/s);

  /*
   * For the isobaric-coordinate case, the external version of w3 is omega = D(p)/Dt, with units [Pa/s].
   * For the isentropic-coordinate case, w3 is D(theta)/Dt both internally and externally, with units [K/s],
   * which are also the units for the hybrid-isentropic case.
   */
  switch(grid.coord_type) {
    case COORD_ISENTROPIC:
    case COORD_HYBRID:
      SET_DIAG(W3_INDEX,w3,coordinate_upward_air_velocity,coordinate vertical velocity,K/s);
    break;
    case COORD_ISOBARIC:
      SET_DIAG(W3_INDEX,w3,lagrangian_tendency_of_air_pressure,omega,Pa/s);
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  SET_DIAG(Z2_INDEX,z2,geopotential_height,geopotential height,m);
  SET_DIAG(Z3_INDEX,z3,geopotential_height,geopotential height,m);
  SET_DIAG(DZDT2_INDEX,dzdt2,upward_air_velocity,vertical velocity,m/s);
  SET_DIAG(DIFFUSION_COEF_UV_INDEX,diffusion_coef_uv,wind_diffusion,wind eddy diffusion,m^2/s);
  SET_DIAG(DIFFUSION_COEF_THETA_INDEX,diffusion_coef_theta,potential_temperature_diffusion,potential temperature eddy diffusion,m^2/s);
  SET_DIAG(DIFFUSION_COEF_MASS_INDEX,diffusion_coef_mass,mass_diffusion,mass eddy diffusion,m^2/s);
  SET_DIAG(PHI_SURFACE_INDEX,phi_surface,surface_geopotential,surface geopotential,m^2/s^2);
  SET_DIAG(GRAVITY2_INDEX,gravity2,gravity,acceleration of gravity,m/s^2);
  SET_DIAG(PBOT_INDEX,pbot,pbot,pressure bottom boundary condition,Pa);

  /*
   * Set species molar_mass values.
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    var.species[is].molar_mass = molar_mass(is);
  }

  /*
   * These functions are defined in the subdirectory epic/src/clouds
   */
  set_species_thermo_data();
  set_microphysics_params(planet->index);
  
  return;
}

/*======================= end of set_var_props() ============================*/

/*======================= free_var_props() ==================================*/
/*
 * Free name-string memory allocated by set_var_props().
 *
 * NOTE: For convenience, keep the same argument list for the FREE_*() macros 
 *       as is used in the SET_*() above.
 */
#define FREE_WIND(wind) \
          if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) { \
            free(var.wind.info[0].name); \
            free(var.wind.info[0].standard_name); \
            free(var.wind.info[0].long_name); \
            free(var.wind.info[0].units); \
            free(var.wind.info_tend[0].name); \
            free(var.wind.info_tend[0].units); \
            free(var.wind.info_tend[1].name); \
            free(var.wind.info_tend[1].units); \
          } \
          else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) { \
            free(var.wind.info[0].name); \
            free(var.wind.info[0].standard_name); \
            free(var.wind.info[0].long_name); \
            free(var.wind.info[0].units); \
            free(var.wind.info[1].name); \
            free(var.wind.info[1].units); \
          }

#define FREE_THERMO(thermo) \
          free(var.thermo.info[0].name); \
          free(var.thermo.info[0].standard_name); \
          free(var.thermo.info[0].long_name); \
          free(var.thermo.info[0].units);

#define FREE_SPECIES(ispecies) \
          free(var.species[ispecies].info[0].name); \
          free(var.species[ispecies].info[0].standard_name); \
          free(var.species[ispecies].info[0].long_name); \
          for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) { \
            free(var.species[ispecies].phase[ip].info[MASS ].name); \
            free(var.species[ispecies].phase[ip].info[MASS ].units); \
            free(var.species[ispecies].phase[ip].info[MOLAR].name); \
            free(var.species[ispecies].phase[ip].info[MOLAR].units); \
          }

#define FREE_DIAG(diag) \
          free(var.diag.info[0].name); \
          free(var.diag.info[0].standard_name); \
          free(var.diag.info[0].long_name); \
          free(var.diag.info[0].units);

/*
 */
void free_var_props(planetspec *planet)
{
  int
    im,is,ip;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="free_var_props";

  /*
   * Solar longitude, L_s.
   */
  free(var.l_s.info.name);
  free(var.l_s.info.long_name);
  free(var.l_s.info.units);

  /* 
   * Wind variables.
   */
  FREE_WIND(u);
  FREE_WIND(v);
  /*
   * Thermo variables.
   */
  FREE_THERMO(h);
  FREE_THERMO(theta);
  FREE_THERMO(fpara);
  /* 
   * Optional species.
   */
  FREE_SPECIES(H_2O_INDEX);
  FREE_SPECIES(NH_3_INDEX);
  FREE_SPECIES(H_2S_INDEX);
  FREE_SPECIES(CH_4_INDEX);
  FREE_SPECIES(C_2H_2_INDEX);
  FREE_SPECIES(C_2H_4_INDEX);
  FREE_SPECIES(C_2H_6_INDEX);
  FREE_SPECIES(CO_2_INDEX);
  FREE_SPECIES(NH_4SH_INDEX);
  FREE_SPECIES(O_3_INDEX);
  FREE_SPECIES(N_2_INDEX);
  FREE_SPECIES(PH_3_INDEX);

  /*
   * Turbulence-model variables.
   */
  if (strcmp(grid.turbulence_scheme,"on")               == 0 ||
      strcmp(grid.turbulence_scheme,"on_vertical_only") == 0)  {
    FREE_THERMO(nu_turb);
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  /* 
   * Diagnostic variables.
   */
  FREE_DIAG(hdry2);
  FREE_DIAG(hdry3);
  FREE_DIAG(h3);
  FREE_DIAG(pdry3);
  FREE_DIAG(p2);
  FREE_DIAG(p3);
  FREE_DIAG(theta2);
  FREE_DIAG(t2);
  FREE_DIAG(t3);
  FREE_DIAG(rho2);
  FREE_DIAG(rho3);
  FREE_DIAG(exner2);
  FREE_DIAG(exner3);
  FREE_DIAG(fgibb2);
  FREE_DIAG(phi2);
  FREE_DIAG(phi3);
  FREE_DIAG(mont2);
  FREE_DIAG(heat3);
  FREE_DIAG(pv2);
  FREE_DIAG(eddy_pv2);
  FREE_DIAG(molar_mass3);
  FREE_DIAG(ri2);
  FREE_DIAG(rel_vort2);
  FREE_DIAG(eddy_rel_vort2);
  FREE_DIAG(abs_vort2);
  FREE_DIAG(kinetic_energy2);
  FREE_DIAG(div_uv2);
  FREE_DIAG(w3);
  FREE_DIAG(z2);
  FREE_DIAG(z3);
  FREE_DIAG(dzdt2);
  FREE_DIAG(diffusion_coef_uv);
  FREE_DIAG(diffusion_coef_theta);
  FREE_DIAG(diffusion_coef_mass);
  FREE_DIAG(phi_surface);
  FREE_DIAG(gravity2);
  FREE_DIAG(pbot);

  return;
}

/*======================= end of free_var_props() ===========================*/

/*======================= make_arrays() =====================================*/

void make_arrays(planetspec *planet)
/*
 * Allocate memory for variables.
 *
 * NOTE: Call var_read() with portion = SIZE_DATA
 *       before calling make_arrays(), to get size of model 
 *       before allocating memory.
 */
{
  int    
    is,ip,kk,
    itmp;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="make_arrays";

#if defined(EPIC_MPI)
  /* Initialize the parallel-bookkeeping structure */
  mpispec_init();

  /*
   * One should use the macros IS_NPOLE and IS_SPOLE,
   * but just in case, synchronize grid.is_npole to para.is_npole and
   * grid.is_spole to para.is_spole.
   */
  grid.is_npole = para.is_npole;
  grid.is_spole = para.is_spole;
#else
  grid.we_num_nodes = 1;

  /* Indicate whether the single-processor job has poles in its range: */
  if (strcmp(grid.geometry,"globe") == 0) {
    if (fcmp(grid.globe_latbot,-90.) == 0) {
      grid.is_spole = TRUE;
    }
    else {
      grid.is_spole = FALSE;
    }
    if (fcmp(grid.globe_lattop,90.) == 0) {
      grid.is_npole = TRUE;
    }
    else {
      grid.is_npole = FALSE;
    }

  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    grid.is_spole = FALSE;
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      grid.is_npole = FALSE;
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      grid.is_npole = TRUE;
    }
  }
  else {
    grid.is_spole = FALSE;
    grid.is_npole = FALSE;
  }
#endif

  /* 
   * Set the shift integers used in the multidimensional-array shift macros:
   */
  Ishift  = ILO-IPAD;
  Jshift  = JLO-JPAD;
  Kshift  = KLO-KPAD;
  Iadim   = IADIM;
  Jadim   = JADIM;
  Kadim   = KADIM;
  Nelem2d = NELEM2D;
  Nelem3d = NELEM3D;
  Shift2d = Ishift+Iadim*Jshift;
  Shift3d = Ishift+Iadim*Jshift+Nelem2d*Kshift;
  Shiftkj = Jshift+Jadim*Kshift;

  /* 
   * Three tendency (d/dt) time planes are used in the 3rd Order 
   * Adams-Bashforth timestep for u and v.
   *
   * Pointers to the time planes.
   *
   * NOTE: IT_ZERO should be initialized to zero here.
   */
  IT_ZERO   = 0;
  IT_MINUS1 = 1;
  IT_MINUS2 = 2;

  /*
   * Allocate memory for prognostic variables and tendencies.
   *
   * NOTE:  The on_list[] variable is a tri-state variable with values: 
   *          NOT_LISTED, LISTED_AND_OFF, LISTED_AND_ON
   *        Hence, it should not be tested like a boolean.
   */

  if (var.on_list[U_INDEX] == LISTED_AND_ON) {
    var.u.on       = TRUE;
    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      var.u.value    = fvector(0,  Nelem3d-1,dbmsname);
      var.u.tendency = fvector(0,3*Nelem3d-1,dbmsname);
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      var.u.value    = fvector(0,2*Nelem3d-1,dbmsname);
      var.u.tendency = fvector(0,  Nelem3d-1,dbmsname);
    }
    if (var.extract_on_list[U_INDEX] == LISTED_AND_ON) var.u.extract_on = TRUE;
  }

  if (var.on_list[V_INDEX] == LISTED_AND_ON) {
    var.v.on       = TRUE;
    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      var.v.value    = fvector(0,  Nelem3d-1,dbmsname);
      var.v.tendency = fvector(0,3*Nelem3d-1,dbmsname);
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      var.v.value    = fvector(0,2*Nelem3d-1,dbmsname);
      var.v.tendency = fvector(0,  Nelem3d-1,dbmsname);
    }
    if (var.extract_on_list[V_INDEX] == LISTED_AND_ON) var.v.extract_on = TRUE;
  }

  if (var.on_list[H_INDEX] == LISTED_AND_ON) {
    var.h.on    = TRUE;
    var.h.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[H_INDEX] == LISTED_AND_ON) var.h.extract_on = TRUE;
  }

  if (var.on_list[THETA_INDEX] == LISTED_AND_ON) {
    var.theta.on    = TRUE;
    var.theta.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[THETA_INDEX] == LISTED_AND_ON) var.theta.extract_on = TRUE;
  }

  if (var.on_list[FPARA_INDEX] == LISTED_AND_ON) {
    var.fpara.on    = TRUE;
    var.fpara.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[FPARA_INDEX] == LISTED_AND_ON) var.fpara.extract_on = TRUE;
  }

  /*
   * Turn on phases appropriate to choice of physics package.
   * The phases are switched on here, whether or not the species are invoked. 
   */
  turn_on_phases(planet);
  
  /*
   * NOTE: Do not use grid.nq yet, it is set below.
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.on_list[is] == LISTED_AND_ON) {
      var.species[is].on = TRUE;
      if (var.extract_on_list[is] == LISTED_AND_ON) var.species[is].extract_on = TRUE;
      for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
        /*
         * The appropriate phases for the chosen physics package should already be turned on.
         */
        if (var.species[is].phase[ip].on) {
          /* 
           * Allocate memory for mass mixing ratio, Q [density_i/density_dry_air],
           * and number fraction, X [n_i/n_total].
           */
          var.species[is].phase[ip].q = fvector(0,Nelem3d-1,dbmsname);
          var.species[is].phase[ip].x = fvector(0,Nelem3d-1,dbmsname);

          if (var.extract_on_list[is] == LISTED_AND_ON) {
            if (grid.cloud_microphysics == ACTIVE ||
                ip == VAPOR) {
              var.species[is].phase[ip].extract_on = TRUE;
            }
            else {
              var.species[is].phase[ip].extract_on = FALSE;
            }
          }
          else {
            var.species[is].phase[ip].extract_on = FALSE;
          }
        }
      }
    }
  }
  if (var.on_list[NU_TURB_INDEX] == LISTED_AND_ON) {
    var.nu_turb.on    = TRUE;
    var.nu_turb.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[NU_TURB_INDEX] == LISTED_AND_ON) var.nu_turb.extract_on = TRUE;
  }

  /*
   * Count total number of active species-phases, grid.nq.
   */
  grid.nq = 0;
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
        if (var.species[is].phase[ip].on) {
          if (grid.cloud_microphysics == ACTIVE ||
              ip == VAPOR) {
            grid.nq++;
          }
        }
      }
    }
  }

  if (grid.nq > 0) {
    /* 
     * Allocate memory for grid.is[], grid.ip[] arrays.
     */
    grid.is = ivector(0,grid.nq-1,dbmsname);
    grid.ip = ivector(0,grid.nq-1,dbmsname);
    /*
     * Assign active species and phase index arrays.
     */
    itmp = 0;
    for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
      if (var.species[is].on) {
        for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
          if (var.species[is].phase[ip].on) {
            if (grid.cloud_microphysics == ACTIVE ||
                ip == VAPOR) {
              grid.is[itmp] = is;
              grid.ip[itmp] = ip;
              itmp++;
            }
          }
        }
      }
    }
  }

  /*
   * Allocate memory for 2D arrays.
   */
  grid.rln = (EPIC_FLOAT **)calloc(2*grid.nk+2,sizeof(EPIC_FLOAT *));
  if (!grid.rln) epic_error(dbmsname,"calloc error allocating grid.rln[kk]");

  grid.rlt = (EPIC_FLOAT **)calloc(2*grid.nk+2,sizeof(EPIC_FLOAT *));
  if (!grid.rlt) epic_error(dbmsname,"calloc error allocating grid.rlt[kk]");

  grid.m =   (EPIC_FLOAT **)calloc(2*grid.nk+2,sizeof(EPIC_FLOAT *));
  if (!grid.m)   epic_error(dbmsname,"calloc error allocating grid.m[kk]");

  grid.n =   (EPIC_FLOAT **)calloc(2*grid.nk+2,sizeof(EPIC_FLOAT *));
  if (!grid.n)   epic_error(dbmsname,"calloc error allocating grid.n[kk]" );

  grid.mn =  (EPIC_FLOAT **)calloc(2*grid.nk+2,sizeof(EPIC_FLOAT *));
  if (!grid.mn)  epic_error(dbmsname,"calloc error allocating grid.mn[kk]");

  grid.beta = (EPIC_FLOAT **)calloc(2*grid.nk+2,sizeof(EPIC_FLOAT *));
  if (!grid.beta) epic_error(dbmsname,"calloc error allocating grid.beta[kk]");

  grid.g =   (EPIC_FLOAT **)calloc(2*grid.nk+2,sizeof(EPIC_FLOAT *));
  if (!grid.g)   epic_error(dbmsname,"calloc error allocating grid.g[kk]");

  for (kk = 0; kk <= 2*grid.nk+1; kk++) {
    grid.rln[ kk] = fvector(0,2*(grid.nj+1),dbmsname);
    grid.rlt[ kk] = fvector(0,2*(grid.nj+1),dbmsname);
    grid.m[   kk] = fvector(0,2*(grid.nj+1),dbmsname);
    grid.n[   kk] = fvector(0,2*(grid.nj+1),dbmsname);
    grid.mn[  kk] = fvector(0,2*(grid.nj+1),dbmsname);
    grid.beta[kk] = fvector(0,2*(grid.nj+1),dbmsname);
    grid.g[   kk] = fvector(0,2*(grid.nj+1),dbmsname);
  }

  /*
   * Allocate memory for 1D arrays.
   */
  grid.lon        = fvector(0,2*(grid.ni+1),dbmsname);

  grid.lat        = fvector(0,2*(grid.nj+1),dbmsname); 
  grid.f          = fvector(0,2*(grid.nj+1),dbmsname);
  grid.f2         = fvector(0,2*(grid.nj+1),dbmsname);

  grid.sigmatheta = dvector(0,2*(grid.nk+1)+1,dbmsname);
  grid.dsgth      = fvector(0,2*(grid.nk+1)+1,dbmsname);
  grid.dsgth_inv  = fvector(0,2*(grid.nk+1)+1,dbmsname);
  grid.p_ref      = fvector(0,2*(grid.nk+1)+1,dbmsname);
  grid.t_ref      = fvector(0,2*(grid.nk+1)+1,dbmsname);
  grid.rho_ref    = fvector(0,2*(grid.nk+1)+1,dbmsname);
  grid.theta_ref  = fvector(0,2*(grid.nk+1)+1,dbmsname);
  grid.h_min      = fvector(0,   grid.nk+1,   dbmsname);
  grid.re         = fvector(0,   grid.nk+1,   dbmsname);
  grid.rp         = fvector(0,   grid.nk+1,   dbmsname);

  if (var.ntp > 0) {
    var.pdat  = fvector(0,var.ntp-1,dbmsname);
    var.tdat  = fvector(0,var.ntp-1,dbmsname);
    var.dtdat = fvector(0,var.ntp-1,dbmsname);
  }

  /*
   * Allocate memory for diagnostic arrays and parameter arrays.
   */
  if (var.n_t_cool > 0) {
    var.t_cool_table = ftriplet(0,var.n_t_cool-1,dbmsname);
  }

  var.gravity2.on    = TRUE;
  var.gravity2.value = fvector(0,KADIM*JADIM-1,dbmsname); 

  if (strcmp(planet->type,"gas-giant") == 0) {
    var.phi_surface.on                     = FALSE;
    var.extract_on_list[PHI_SURFACE_INDEX] = NOT_LISTED;
    var.phi_surface.extract_on             = FALSE;

    var.pbot.on    = TRUE;
    var.pbot.value = fvector(0,NELEM2D-1,dbmsname);
    if (var.extract_on_list[PBOT_INDEX] == LISTED_AND_ON) var.pbot.extract_on = TRUE;
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    var.phi_surface.on    = TRUE;
    var.phi_surface.value = fvector(0,NELEM2D-1,dbmsname);
    if (var.extract_on_list[PHI_SURFACE_INDEX] == LISTED_AND_ON) var.phi_surface.extract_on = TRUE;

    var.pbot.on                     = FALSE;
    var.extract_on_list[PBOT_INDEX] = NOT_LISTED;
    var.pbot.extract_on             = FALSE;
  }
  else {
    sprintf(Message,"unrecognized planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }

  var.hdry2.on = TRUE;
  if (grid.nq > 0) {
    var.hdry2.value = fvector(0,Nelem3d-1,dbmsname);
  }
  else {
    /*
     * No difference between HDRY2 and H, so assign them the same memory.
     */
    var.hdry2.value = var.h.value;
  }
  if (var.extract_on_list[HDRY2_INDEX] == LISTED_AND_ON) var.hdry2.extract_on = TRUE;

  var.hdry3.on       = TRUE;
  var.hdry3.value    = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[HDRY3_INDEX] == LISTED_AND_ON) var.hdry3.extract_on = TRUE;

  var.h3.on = TRUE;
  if (grid.nq > 0) {
    var.h3.value = fvector(0,Nelem3d-1,dbmsname);
  }
  else {
    /*
     * No difference between HDRY3 and H3, so assign them the same memory.
     */
    var.h3.value = var.hdry3.value;
  }
  if (var.extract_on_list[H3_INDEX] == LISTED_AND_ON) var.h3.extract_on = TRUE;

  var.p2.on          = TRUE;
  var.p2.value       = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[P2_INDEX] == LISTED_AND_ON) var.p2.extract_on = TRUE;

  var.pdry3.on       = TRUE;
  var.pdry3.value    = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[PDRY3_INDEX] == LISTED_AND_ON) var.pdry3.extract_on = TRUE;
 
  var.p3.on          = TRUE;
  if (grid.nq > 0) {
    var.p3.value = fvector(0,Nelem3d-1,dbmsname);
  }
  else {
    /*
     * No difference between PDRY3 and P3, so assign them the same memory.
     */
    var.p3.value = var.pdry3.value;
  }
  if (var.extract_on_list[P3_INDEX] == LISTED_AND_ON) var.p3.extract_on = TRUE;

  var.theta2.on      = TRUE;
  var.theta2.value   = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[THETA2_INDEX] == LISTED_AND_ON) var.theta2.extract_on = TRUE;

  var.t2.on          = TRUE;
  var.t2.value       = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[T2_INDEX] == LISTED_AND_ON) var.t2.extract_on = TRUE;

  var.t3.on          = TRUE;
  var.t3.value       = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[T3_INDEX] == LISTED_AND_ON) var.t3.extract_on = TRUE;

  var.rho2.on        = TRUE;
  var.rho2.value     = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[RHO2_INDEX] == LISTED_AND_ON) var.rho2.extract_on = TRUE;

  var.rho3.on        = TRUE;
  var.rho3.value     = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[RHO3_INDEX] == LISTED_AND_ON) var.rho3.extract_on = TRUE;

  var.exner2.on      = TRUE;
  var.exner2.value   = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[EXNER2_INDEX] == LISTED_AND_ON) var.exner2.extract_on = TRUE;

  var.exner3.on      = TRUE;
  var.exner3.value   = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[EXNER3_INDEX] == LISTED_AND_ON) var.exner3.extract_on = TRUE;

  var.phi2.on        = TRUE;
  var.phi2.value     = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[PHI2_INDEX] == LISTED_AND_ON) var.phi2.extract_on = TRUE;

  var.phi3.on        = TRUE;
  var.phi3.value     = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[PHI3_INDEX] == LISTED_AND_ON) var.phi3.extract_on = TRUE;

  var.mont2.on       = TRUE;
  var.mont2.value    = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[MONT2_INDEX] == LISTED_AND_ON) var.mont2.extract_on = TRUE;

  var.heat3.on       = TRUE;
  var.heat3.value    = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[HEAT3_INDEX] == LISTED_AND_ON) var.heat3.extract_on = TRUE;

  var.pv2.on         = TRUE;
  var.pv2.value      = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[PV2_INDEX] == LISTED_AND_ON) var.pv2.extract_on = TRUE;

  var.ri2.on         = TRUE;
  var.ri2.value      = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[RI2_INDEX] == LISTED_AND_ON) var.ri2.extract_on = TRUE;

  var.div_uv2.on     = TRUE;
  var.div_uv2.value  = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[DIV_UV2_INDEX] == LISTED_AND_ON) var.div_uv2.extract_on = TRUE;

  /* 
   * The following diagnostic variables are available for output, but are
   * not assigned permanent memory.
   */
  if (var.extract_on_list[EDDY_PV2_INDEX] == LISTED_AND_ON) {
    var.eddy_pv2.extract_on = TRUE;
  }
  if (var.extract_on_list[REL_VORT2_INDEX] == LISTED_AND_ON) {
    var.rel_vort2.extract_on = TRUE;
  }
  if (var.extract_on_list[EDDY_REL_VORT2_INDEX] == LISTED_AND_ON) {
    var.eddy_rel_vort2.extract_on = TRUE;
  }
  if (var.extract_on_list[ABS_VORT2_INDEX] == LISTED_AND_ON) {
    var.abs_vort2.extract_on = TRUE;
  }
  if (var.extract_on_list[KIN2_INDEX] == LISTED_AND_ON) {
    var.kinetic_energy2.extract_on = TRUE;
  }
  if (var.extract_on_list[MOLAR_MASS3_INDEX] == LISTED_AND_ON) {
    var.molar_mass3.extract_on = TRUE;
  }

  var.w3.on          = TRUE;
  var.w3.value       = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[W3_INDEX] == LISTED_AND_ON) var.w3.extract_on = TRUE;

  var.z2.on          = TRUE;
  var.z2.value       = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[Z2_INDEX] == LISTED_AND_ON) var.z2.extract_on = TRUE;

  var.z3.on          = TRUE;
  var.z3.value       = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[Z3_INDEX] == LISTED_AND_ON) var.z3.extract_on = TRUE;

  var.dzdt2.on       = TRUE;
  var.dzdt2.value    = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[DZDT2_INDEX] == LISTED_AND_ON) var.dzdt2.extract_on = TRUE;

  var.diffusion_coef_uv.on    = TRUE;
  var.diffusion_coef_uv.value = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[DIFFUSION_COEF_UV_INDEX] == LISTED_AND_ON) var.diffusion_coef_uv.extract_on = TRUE;

  var.diffusion_coef_theta.on    = TRUE;
  var.diffusion_coef_theta.value = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[DIFFUSION_COEF_THETA_INDEX] == LISTED_AND_ON) var.diffusion_coef_theta.extract_on = TRUE;

  var.diffusion_coef_mass.on     = TRUE;
  var.diffusion_coef_mass.value  = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[DIFFUSION_COEF_MASS_INDEX] == LISTED_AND_ON) var.diffusion_coef_mass.extract_on = TRUE;

  var.fgibb2.on    = TRUE;
  var.fgibb2.value = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[FGIBB2_INDEX] == LISTED_AND_ON) var.fgibb2.extract_on = TRUE;

  /*
   * Allocate turbulence-model memory.
   */
  if (strcmp(grid.turbulence_scheme,"on")               == 0 ||
      strcmp(grid.turbulence_scheme,"on_vertical_only") == 0)  {
    make_arrays_subgrid();
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  return;
}

/*====================== end of make_arrays() ===============================*/

/*====================== free_arrays() ======================================*/

void free_arrays(planetspec *planet)
/*
 * Free memory allocated by make_arrays().
 */
{
  int    
    iq,kk;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="free_arrays";

  if (var.u.on) {
    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      free_fvector(var.u.value,   0,  Nelem3d-1,dbmsname);
      free_fvector(var.u.tendency,0,3*Nelem3d-1,dbmsname);
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      free_fvector(var.u.value,   0,2*Nelem3d-1,dbmsname);
      free_fvector(var.u.tendency,0,  Nelem3d-1,dbmsname);
    }
  }
  if (var.v.on) {
    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      free_fvector(var.v.value,   0,  Nelem3d-1,dbmsname);
      free_fvector(var.v.tendency,0,3*Nelem3d-1,dbmsname);
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      free_fvector(var.v.value,   0,2*Nelem3d-1,dbmsname);
      free_fvector(var.v.tendency,0,  Nelem3d-1,dbmsname);
    }
  }
  if (var.h.on) {
    free_fvector(var.h.value,0,Nelem3d-1,dbmsname);
  }
  if (var.theta.on) {
    free_fvector(var.theta.value,0,Nelem3d-1,dbmsname);
  }
  if (var.fpara.on) {
    free_fvector(var.fpara.value,0,Nelem3d-1,dbmsname);
  }
  for (iq = 0; iq < grid.nq; iq++) {
    free_fvector(var.species[grid.is[iq]].phase[grid.ip[iq]].q,0,Nelem3d-1,dbmsname);
    free_fvector(var.species[grid.is[iq]].phase[grid.ip[iq]].x,0,Nelem3d-1,dbmsname);
  }
  if (var.nu_turb.on) {
    free_fvector(var.nu_turb.value,0,Nelem3d-1,dbmsname);
  }

  if (grid.nq > 0) {
    free_ivector(grid.is,0,grid.nq-1,dbmsname);
    free_ivector(grid.ip,0,grid.nq-1,dbmsname);
  }

  for (kk = 0; kk <= 2*grid.nk+1; kk++) {
    free_fvector(grid.rln[ kk],0,2*(grid.nj+1),dbmsname);
    free_fvector(grid.rlt[ kk],0,2*(grid.nj+1),dbmsname);
    free_fvector(grid.m[   kk],0,2*(grid.nj+1),dbmsname);
    free_fvector(grid.n[   kk],0,2*(grid.nj+1),dbmsname);
    free_fvector(grid.mn[  kk],0,2*(grid.nj+1),dbmsname);
    free_fvector(grid.beta[kk],0,2*(grid.nj+1),dbmsname);
    free_fvector(grid.g[   kk],0,2*(grid.nj+1),dbmsname);
  }
  free(grid.rln );
  free(grid.rlt );
  free(grid.m   );
  free(grid.n   );
  free(grid.mn  );
  free(grid.beta);
  free(grid.g   );

  free_fvector(grid.lon,       0,2*(grid.ni+1),  dbmsname);
  free_fvector(grid.lat,       0,2*(grid.nj+1),  dbmsname); 
  free_fvector(grid.f,         0,2*(grid.nj+1),  dbmsname);
  free_fvector(grid.f2,        0,2*(grid.nj+1),  dbmsname);
  free_dvector(grid.sigmatheta,0,2*(grid.nk+1)+1,dbmsname);
  free_fvector(grid.dsgth,     0,2*(grid.nk+1)+1,dbmsname);
  free_fvector(grid.dsgth_inv, 0,2*(grid.nk+1)+1,dbmsname);
  free_fvector(grid.p_ref,     0,2*(grid.nk+1)+1,dbmsname);
  free_fvector(grid.t_ref,     0,2*(grid.nk+1)+1,dbmsname);
  free_fvector(grid.rho_ref,   0,2*(grid.nk+1)+1,dbmsname);
  free_fvector(grid.theta_ref, 0,2*(grid.nk+1)+1,dbmsname);
  free_fvector(grid.h_min,     0,   grid.nk+1,   dbmsname);
  free_fvector(grid.re,        0,   grid.nk+1,   dbmsname);
  free_fvector(grid.rp,        0,   grid.nk+1,   dbmsname);

  if (var.ntp > 0) {
    free_fvector(var.pdat, 0,var.ntp-1,dbmsname);
    free_fvector(var.tdat, 0,var.ntp-1,dbmsname);
    free_fvector(var.dtdat,0,var.ntp-1,dbmsname);
  }

  if (var.n_t_cool > 0) free_ftriplet(var.t_cool_table,0,var.n_t_cool-1,dbmsname);

  free_fvector(var.gravity2.value,0,KADIM*JADIM-1,dbmsname); 

  if (strcmp(planet->type,"gas-giant") == 0) {
    free_fvector(var.pbot.value,0,NELEM2D-1,dbmsname);
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    free_fvector(var.phi_surface.value,0,NELEM2D-1,dbmsname);
  }

  free_fvector(var.pdry3.value,               0,Nelem3d-1,dbmsname);
  free_fvector(var.p2.value,                  0,Nelem3d-1,dbmsname);
  free_fvector(var.theta2.value,              0,Nelem3d-1,dbmsname);
  if (grid.nq > 0) {
    free_fvector(var.hdry2.value,             0,Nelem3d-1,dbmsname);
    free_fvector(var.h3.value,                0,Nelem3d-1,dbmsname);
    free_fvector(var.p3.value,                0,Nelem3d-1,dbmsname);
  }
  free_fvector(var.t2.value,                  0,Nelem3d-1,dbmsname);
  free_fvector(var.t3.value,                  0,Nelem3d-1,dbmsname);
  free_fvector(var.rho2.value,                0,Nelem3d-1,dbmsname);
  free_fvector(var.rho3.value,                0,Nelem3d-1,dbmsname);
  free_fvector(var.exner2.value,              0,Nelem3d-1,dbmsname);
  free_fvector(var.exner3.value,              0,Nelem3d-1,dbmsname);
  free_fvector(var.phi2.value,                0,Nelem3d-1,dbmsname);
  free_fvector(var.phi3.value,                0,Nelem3d-1,dbmsname);
  free_fvector(var.mont2.value,               0,Nelem3d-1,dbmsname);
  free_fvector(var.heat3.value,               0,Nelem3d-1,dbmsname);
  free_fvector(var.pv2.value,                 0,Nelem3d-1,dbmsname);
  free_fvector(var.ri2.value,                 0,Nelem3d-1,dbmsname);
  free_fvector(var.div_uv2.value,             0,Nelem3d-1,dbmsname);
  free_fvector(var.w3.value,                  0,Nelem3d-1,dbmsname);
  free_fvector(var.z2.value,                  0,Nelem3d-1,dbmsname);
  free_fvector(var.z3.value,                  0,Nelem3d-1,dbmsname);
  free_fvector(var.dzdt2.value,               0,Nelem3d-1,dbmsname);
  free_fvector(var.diffusion_coef_uv.value,   0,Nelem3d-1,dbmsname);
  free_fvector(var.diffusion_coef_theta.value,0,Nelem3d-1,dbmsname);
  free_fvector(var.diffusion_coef_mass.value, 0,Nelem3d-1,dbmsname);
  if (var.fpara.on) {
    free_fvector(var.fgibb2.value,0,Nelem3d-1,dbmsname);
  }

  /*
   * Free turbulence-model memory.
   */
  if (strcmp(grid.turbulence_scheme,"on")               == 0 ||
      strcmp(grid.turbulence_scheme,"on_vertical_only") == 0)  {
    free_arrays_subgrid();
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  return;
}

/*======================= end free_arrays() =================================*/

/*======================= return_sigmatheta() ===============================*/

/*
 * Generalized coordinate, sigmatheta, as a function of pressure and potential
 * temperature. Use double precision to increase accuracy of diagnostic theta
 * calculations.
 */

double return_sigmatheta(register double theta,
                         register double p,
                         register double pbot)
{
  register double
    sigma,
    sigmatheta;

  if (p <= pbot) {
    sigma      = get_sigma(pbot,p);
    sigmatheta = f_sigma(sigma)+g_sigma(sigma)*theta;
  }
  else {
    sprintf(Message,"pbot=%g, ptop=%g; p=%g out of range",pbot,grid.ptop,p);
    epic_error("return_sigmatheta",Message);
  }

  return sigmatheta;
}

/*======================= end of return_sigmatheta() ========================*/

/*======================= f_sigma() =========================================*/
/*
 * Hybrid vertical coordinate definition:
 *   sigmatheta = f(sigma)+g(sigma)*theta.
 * Use double precision to increase accuracy of diagnostic theta calculations.
 */

double f_sigma(double sigma)
{
  switch(grid.coord_type) {
    case COORD_ISENTROPIC:
      return 0.;
    break;
    case COORD_HYBRID:
      if (sigma < 0.) {
        return grid.zeta0;
      }
      else if (sigma < 1.) {
        return (1.-g_sigma(sigma))*(grid.zeta0+sigma*(grid.zeta1-grid.zeta0));
      }
      else {
        return 0.;
      }
    break;
    case COORD_ISOBARIC:
      return grid.zeta0+sigma*(grid.zeta1-grid.zeta0);
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error("f_sigma",Message);
    break;
  }

  epic_error("f_sigma","not supposed to reach this point");
  return DBL_MAX;
}

/*======================= end of f_sigma() ==================================*/

/*======================= g_sigma() =========================================*/
/*
 * Part of the hybrid vertical coordinate definition,
 * sigmatheta = f(sigma)+g(sigma)*theta.
 * Use double precision to increase accuracy of diagnostic theta calculations.
 */

double g_sigma(double sigma)
{
  const double
    a = tanh(-grid.hybrid_alpha*grid.sigma_sigma),
    b = 1./(tanh(grid.hybrid_alpha*(1.-grid.sigma_sigma))-a);

  switch(grid.coord_type) {
    case COORD_ISENTROPIC:
      return 1.;
    break;
    case COORD_HYBRID:
      return LIMIT_RANGE(0.,(tanh(grid.hybrid_alpha*(sigma-grid.sigma_sigma))-a)*b,1.);
    break;
    case COORD_ISOBARIC:
      return 0.;
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error("g_sigma",Message);
    break;
  }

  epic_error("g_sigma","not supposed to reach this point");
  return DBL_MAX;
}

/*======================= end of g_sigma() ==================================*/

/*======================= set_lonlat() ======================================*/
  
/*
 *  Set longitude and latitude.
 */

void set_lonlat(void)
{
  int
    nj,ni,
    jj,ii;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_lonlat";
  
  nj = grid.nj;
  ni = grid.ni;

  /*
   *  Compute lon:
   */
  if (strcmp(grid.geometry,"globe") == 0) {
    grid.lon[0] = grid.globe_lonbot-grid.dln;
  }
  else if (strcmp(grid.geometry,"f-plane") == 0 && strcmp(grid.f_plane_map,"cartesian") == 0) {
    grid.lon[0] = -90.-grid.dln;
  }
  else {
    grid.lon[0] = -180.-grid.dln;
  }
  for (ii = 1; ii <= 2*(ni+1); ii++) {
    grid.lon[ii] = grid.lon[ii-1]+grid.dln*.5;
  }

  /*
   * Compute lat:
   */
  if (strcmp(grid.geometry,"globe") == 0) {
    if (grid.globe_latbot == -90.) {
      grid.lat[0] = -90.+grid.dlt*sqrt(.5);
    }
    else {
      (grid.lat)[0] = grid.globe_latbot;
    }
  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      grid.lat[0] = -90.-grid.dlt;
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      grid.lat[0] = 0.-grid.dlt;
    }
  }
  else {
    sprintf(Message,"unrecognized geometry %s",grid.geometry);
    epic_error(dbmsname,Message);
  }

  for (jj = 1; jj <= 2*(nj+1); jj++)  {      
    grid.lat[jj] = grid.lat[jj-1]+grid.dlt*.5;
  }

  if (strcmp(grid.geometry,"globe") == 0) {
    if (grid.globe_latbot == -90.) {
      /* poles are offset by extra dlt*sqrt(.5) */
      grid.lat[       0] = -90.;   
    }
    if (grid.globe_lattop == 90.) {  
      grid.lat[2*(nj+1)] =  90.; 
    }
  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    if (strcmp(grid.f_plane_map,"polar") == 0) {
      grid.lat[2*(nj+1)] = 90.;
    }
  }

  if (strcmp(grid.geometry,"globe") == 0) {
    /*
     * Improve mirror-image symmetry across the equator when applicable.
     */
    int
      jjn = 2*(nj+1),
      jjs = 0;

    while (grid.lat[jjs] < 0. && grid.lat[jjn] > 0.) {
      for (; jjs < jjn; jjs++) {
        if (fabs(grid.lat[jjs]+grid.lat[jjn]) < 1.e-3) {
          grid.lat[jjn--] = -grid.lat[jjs++];
          break;            
        }
      }
    }
  }

  return;
}

/*======================= end of set_lonlat() ===============================*/

/*======================= set_fmn() =========================================*/
  
/*
 * Compute the Coriolis parameter, f (and the second Coriolis parameter, f2), 
 * and the geometric map factors, m, n, etc.
 * The map factors are functions of vertical index as well as meridional index.
 * Also compute beta = df/dy.
 */
void set_fmn(planetspec *planet)
{
  int
    K,kk,jj,
    nj,ni;
  EPIC_FLOAT
    omega,re,rp,
    dlnr,dltr,lat,
    rln,rlt,
    dx,dy,tmp,
    lat0,m0,n0;
  float_triplet
   *buff_triplet;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_fmn";

  nj = grid.nj;
  ni = grid.ni;

  omega = planet->omega_sidereal;

  /*
   * Set the Coriolis parameters, which depends on latitude but not height.
   */
  if (strcmp(grid.geometry,"globe") == 0) {
    for (jj = 1; jj < 2*(nj+1); jj++) {
      lat  = DEG*(grid.lat)[jj];
      grid.f[ jj] = 2.*omega*sin(lat);
      grid.f2[jj] = 2.*omega*cos(lat);
    }
    /*
     *  The Arakawa and Lamb (1981) scheme calls for a 
     *  special (3/2)*dlt spacing for n next to the poles,
     *  as illustrated by their Fig. A2; their equation
     *  (A40) specifies mn at the poles.
     *
     *  NOTE: We find that the special spacing next to the poles is more accurately
     *  given by (1.+sqrt(.5))*dlt.  
     */
    if (grid.globe_latbot == -90.) {
      /* south pole */
      lat = DEG*(grid.lat)[0];
      grid.f[ 2*0] = 2.*omega*sin(lat);
      grid.f2[2*0] = 2.*omega*cos(lat);
    }
    else {
      jj = 0;
      lat  = DEG*(grid.lat)[jj];
      grid.f[ jj] = 2.*omega*sin(lat);
      grid.f2[jj] = 2.*omega*cos(lat); 
    }
    if (grid.globe_lattop == 90.) {
      /* north pole */
      lat = DEG*(grid.lat)[2*(nj+1)];
      grid.f[ 2*(nj+1)] = 2.*omega*sin(lat);
      grid.f2[2*(nj+1)] = 2.*omega*cos(lat);
    }
    else {
      jj = 2*(nj+1);
      lat  = DEG*(grid.lat)[jj];
      grid.f[ jj] = 2.*omega*sin(lat);
      grid.f2[jj] = 2.*omega*cos(lat);
    }
  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      lat = DEG*(grid.f_plane_lat0);
      for (jj = 0; jj <= 2*(nj+1); jj++) {
        grid.f[ jj] = 2.*omega*sin(lat);
        grid.f2[jj] = 2.*omega*cos(lat);
      }
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      lat  = DEG*(grid.f_plane_lat0);
      for (jj = 0; jj < 2*(nj+1); jj++) {
        grid.f[ jj] = 2.*omega*sin(lat);
        grid.f2[jj] = 2.*omega*cos(lat);
      }
      /* pole */
      grid.f[ 2*(nj+1)] = 2.*omega*sin(lat);
      grid.f2[2*(nj+1)] = 2.*omega*cos(lat);
    }
  }

  /* 
   * Set the map factors, which depend on latitude and height.
   */
  for (K = 0; K <= grid.nk; K++) {
    kk = 2*K+1;
    if (strcmp(grid.geometry,"globe") == 0) {
      dlnr  = grid.dln*DEG;
      dltr  = grid.dlt*DEG;
      re    = grid.re[K];
      rp    = grid.rp[K];
      for (jj = 1; jj < 2*(nj+1); jj++) {
        lat  = DEG*(grid.lat)[jj];
        rln  = re/sqrt( 1.+SQR(rp/re*tan(lat)) );
        rlt  = rln/( cos(lat)*(SQR(sin(lat))+SQR(re/rp*cos(lat))) );

        grid.rln[kk][jj] = rln;
        grid.rlt[kk][jj] = rlt;
        grid.m[  kk][jj] = 1./(rln*dlnr);
        grid.n[  kk][jj] = 1./(rlt*dltr);
        grid.mn[ kk][jj] = grid.m[kk][jj]*grid.n[kk][jj];
      }
      /*
       *  The Arakawa and Lamb (1981) scheme calls for a 
       *  special (3/2)*dlt spacing for n next to the poles,
       *  as illustrated by their Fig. A2; their equation
       *  (A40) specifies mn at the poles.
       *
       *  NOTE: We find that the special spacing next to the poles is more accurately
       *  given by (1.+sqrt(.5))*dlt.  
       */
      if (grid.globe_latbot == -90.) {
        /* south pole */
        grid.rln[kk][2*0  ]  = 0.;
        grid.rlt[kk][2*0  ]  = re*re/rp;

        /* 
         * The pole value of m is not defined, but is never needed.
         * The pole value of n should not be needed, but is included for completeness.
         */
        tmp                  = grid.n[kk][2*0+1];
        grid.n[  kk][2*0+1]  = tmp/(1.+sqrt(.5));
        grid.n[  kk][2*0  ]  = tmp/(1.+sqrt(2.)); 
        /* wedge shaped area: */
        grid.mn[ kk][2*0+1]  = grid.n[kk][2*0+1]*grid.n[kk][2*0+1]/(.5*dlnr);
        grid.mn[ kk][2*0  ]  = 2.*grid.mn[kk][2*0+1];
      }
      else {
        jj = 0;
        lat  = DEG*(grid.lat)[jj];
        rln  = re/sqrt( 1.+ pow(rp/re*tan(lat),2.) );
        rlt  = rln/( cos(lat)*(pow(sin(lat),2.)+
                     pow(re/rp*cos(lat),2.)) );

        grid.rln[kk][jj] = rln;
        grid.rlt[kk][jj] = rlt;
        grid.m[  kk][jj] = 1./(rln*dlnr);
        grid.n[  kk][jj] = 1./(rlt*dltr);
        grid.mn[ kk][jj] = grid.m[kk][jj]*grid.n[kk][jj]; 
      }
      if (grid.globe_lattop == 90.) {
        /* north pole */
        grid.rln[kk][2*(nj+1)  ]  = 0.;
        grid.rlt[kk][2*(nj+1)  ]  = re*re/rp;

        /* 
         * The pole value of m is not defined, but is never needed.
         * The pole value of n should not be needed, but is included for completeness.
         */
        tmp                       = grid.n[kk][2*(nj+1)-1];
        grid.n[  kk][2*(nj+1)-1]  = tmp/(1.+sqrt(.5));
        grid.n[  kk][2*(nj+1)  ]  = tmp/(1.+sqrt(2.));
        /* wedge shaped area: */
        grid.mn[kk][2*(nj+1)-1]  = grid.n[kk][2*(nj+1)-1]*grid.n[kk][2*(nj+1)-1]/(.5*dlnr);
        grid.mn[kk][2*(nj+1)  ]  = 2.*grid.mn[kk][2*(nj+1)-1];
      }
      else {
        jj = 2*(nj+1);
        lat  = DEG*(grid.lat)[jj];
        rln  = re/sqrt( 1.+ pow(rp/re*tan(lat),(EPIC_FLOAT)2.) );
        rlt  = rln/( cos(lat)*( pow(sin(lat),(EPIC_FLOAT)2.)+
               pow(re/rp*cos(lat),(EPIC_FLOAT)2.)));

        grid.rln[kk][jj] = rln;
        grid.rlt[kk][jj] = rlt;
        grid.m[  kk][jj] = 1./(rln*dlnr);
        grid.n[  kk][jj] = 1./(rlt*dltr);
        grid.mn[ kk][jj] = grid.m[kk][jj]*grid.n[kk][jj];
      }
    }
    else if (strcmp(grid.geometry,"f-plane") == 0) {
      if (strcmp(grid.f_plane_map,"cartesian") == 0) {
        dx  = 2.*(grid.f_plane_half_width)/ni;
        dy  = dx;

        for (jj = 0; jj <= 2*(nj+1); jj++) {
          grid.rln[kk][jj] = 1.;
          grid.rlt[kk][jj] = 1.;
          grid.m[  kk][jj] = 1./dx;
          grid.n[  kk][jj] = 1./dy;
          grid.mn[ kk][jj] = 1./(dx*dy);
        }
      }
      else if (strcmp(grid.f_plane_map,"polar") == 0) {
        dlnr = grid.dln*DEG;
        dy   = grid.f_plane_half_width*(grid.dlt/90.);
        rln  = grid.f_plane_half_width+dy;

        for (jj = 0; jj < 2*(nj+1); jj++) {
          dx = rln*dlnr;
          grid.rln[kk][jj] = rln;
          grid.rlt[kk][jj] = 1.;
          grid.m[  kk][jj] = 1./dx;
          grid.n[  kk][jj] = 1./dy;
          grid.mn[ kk][jj] = grid.m[kk][jj]*grid.n[kk][jj];
          rln -= dy/2.;
        }
        /* pole */
        grid.rln[kk][2*(nj+1)  ]  = 0.;
        grid.rlt[kk][2*(nj+1)  ]  = 1.;

        /* 
         * The pole value of m is not defined, but is never needed.
         * The pole value of n should not be needed, but is included for completeness.
         */
        tmp                       = grid.n[kk][2*(nj+1)-1];
        grid.n[  kk][2*(nj+1)-1]  = tmp/(1.+sqrt(.5));
        grid.n[  kk][2*(nj+1)  ]  = tmp/(1.+sqrt(2.));
        /* wedge shaped area: */
        grid.mn[kk][2*(nj+1)-1]  = grid.n[kk][2*(nj+1)-1]*grid.n[kk][2*(nj+1)-1]/(.5*dlnr);
        grid.mn[kk][2*(nj+1)  ]  = 2.*grid.mn[kk][2*(nj+1)-1];
      }
    }
  }

  /*
   * Fill in layer values using a smooth, monotonic spline.
   */

  /* Allocate memory */
  buff_triplet = ftriplet(0,KHI,dbmsname);
  for (K = 0; K <= KHI; K++) {
    kk                = 2*K+1;
    buff_triplet[K].x = (double)(kk);
  }
  for (jj = 0; jj <= 2*(grid.nj+1); jj++) {
    /* grid.m */
    for (K = 0; K <= KHI; K++) {
      kk                = 2*K+1;
      buff_triplet[K].y = grid.m[kk][jj];
    }
    spline_pchip(grid.nk+1,buff_triplet);
    for (K = KLO; K <= KHI; K++) {
      kk              = 2*K;
      grid.m[kk][jj] = splint_pchip((double)kk,buff_triplet+K-KLO,2.);
    }

    /* grid.n */
    for (K = 0; K <= KHI; K++) {
      kk                = 2*K+1;
      buff_triplet[K].y = grid.n[kk][jj];
    }
    spline_pchip(grid.nk+1,buff_triplet);
    for (K = KLO; K <= KHI; K++) {
      kk              = 2*K;
      grid.n[kk][jj] = splint_pchip((double)kk,buff_triplet+K-KLO,2.);
    }

    /* grid.mn */
    for (K = 0; K <= KHI; K++) {
      kk                = 2*K+1;
      buff_triplet[K].y = grid.mn[kk][jj];
    }
    spline_pchip(grid.nk+1,buff_triplet);
    for (K = KLO; K <= KHI; K++) {
      kk              = 2*K;
      grid.mn[kk][jj] = splint_pchip((double)kk,buff_triplet+K-KLO,2.);
    }

    /* grid.rln */
    for (K = 0; K <= KHI; K++) {
      kk                = 2*K+1;
      buff_triplet[K].y = grid.rln[kk][jj];
    }
    spline_pchip(grid.nk+1,buff_triplet);
    for (K = KLO; K <= KHI; K++) {
      kk               = 2*K;
      grid.rln[kk][jj] = splint_pchip((double)kk,buff_triplet+K-KLO,2.);
    }

    /* grid.rlt */
    for (K = 0; K <= KHI; K++) {
      kk                = 2*K+1;
      buff_triplet[K].y = grid.rlt[kk][jj];
    }
    spline_pchip(grid.nk+1,buff_triplet);
    for (K = KLO; K <= KHI; K++) {
      kk               = 2*K;
      grid.rlt[kk][jj] = splint_pchip((double)kk,buff_triplet+K-KLO,2.);
    }
  }

  /*
   * Compute beta = df/dy
   */
  for (kk = 1; kk <= 2*(grid.nk)+1; kk++) {
    for (jj = 2*grid.jlo+1; jj < 2*(grid.nj+1); jj++) {
      grid.beta[kk][jj] = (grid.f[jj+1]-grid.f[jj-1])*grid.n[kk][jj];
    }
    jj = 2*grid.jlo;
    if (grid.lat[jj] == -90.) {
      grid.beta[kk][jj] = 0.;
      /* 
       * Replace next-to-pole value with interpolation
       * to avoid a glitch from the grid-spacing change.
       */
      grid.beta[kk][jj+1] = grid.beta[kk][jj]+(grid.beta[kk][jj+2]-grid.beta[kk][jj])
                                             *(grid.lat[jj+1]-grid.lat[jj])/(grid.lat[jj+2]-grid.lat[jj]);
    }
    else {
      /*
       * Simple extrapolation for southern boundary.
       */
      grid.beta[kk][jj] = 2.*grid.beta[kk][jj+1]-grid.beta[kk][jj+2];
    }
    jj = 2*(grid.nj+1);
    if (grid.lat[jj] == 90.) {
      grid.beta[kk][jj] = 0.;
      /* 
       * Replace next-to-pole value with interpolation
       * to avoid a glitch from the grid-spacing change.
       */
      grid.beta[kk][jj-1] = grid.beta[kk][jj]+(grid.beta[kk][jj-2]-grid.beta[kk][jj])
                                             *(grid.lat[jj-1]-grid.lat[jj])/(grid.lat[jj-2]-grid.lat[jj]);
    }
    else {
      /*
       * Simple extrapolation for northern boundary.
       */
      grid.beta[kk][jj] = 2.*grid.beta[kk][jj-1]-grid.beta[kk][jj-2];
    }
  }

  /* Free allocated memory */
  free_ftriplet(buff_triplet,0,KHI,dbmsname);

  return;
}

/*======================= end of set_fmn() ==================================*/

/*======================= set_gravity() =====================================*/

/*
 * Calculate gravity, g [m/s^2], as a function of planetographic latitude.
 *
 * NOTE: In this version, gravity does not vary in the vertical dimension.
 *
 * See eqn (41) of Yoder C, 1995, Global Earth Physics: A handbook of physical constants,
 *   http://www.agu.org/reference/gephys/4_yoder.pdf
 */

void set_gravity(planetspec *planet)
{
  register int
    K,J,
    kk,jj;
  register double
    a,c,ge,gp,
    spin_factor,g_avg,
    sinlat2,coslat2;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_gravity";

  a           = planet->re;
  c           = planet->rp;
  spin_factor = (planet->omega_sidereal*a)*(planet->omega_sidereal*a)*(a/planet->GM);
  ge          = (1.+1.5*planet->J2-spin_factor)*planet->GM/(a*a);
  gp          = (1.-3.*(a/c)*(a/c)*planet->J2)*planet->GM/(c*c);

  if (strcmp(grid.geometry,"f-plane") == 0) {
    sinlat2  = sin(grid.f_plane_lat0*DEG);
    sinlat2 *= sinlat2;
    coslat2  = cos(grid.f_plane_lat0*DEG);
    coslat2 *= coslat2;
    for (jj = 0; jj <= 2*(grid.nj+1); jj++) {
      for (kk = 1; kk <= 2*grid.nk+1; kk++) {
        grid.g[kk][jj] = (a*ge*coslat2+c*gp*sinlat2)/sqrt(a*a*coslat2+c*c*sinlat2);
      }
    }
  }
  else {
    for (jj = 0; jj <= 2*(grid.nj+1); jj++) {
      sinlat2  = sin(grid.lat[jj]*DEG);
      sinlat2 *= sinlat2;
      coslat2  = cos(grid.lat[jj]*DEG);
      coslat2 *= coslat2;
      for (kk = 1; kk <= 2*grid.nk+1; kk++) {
        grid.g[kk][jj] = (a*ge*coslat2+c*gp*sinlat2)/sqrt(a*a*coslat2+c*c*sinlat2);
      }
    }
  }

  /*
   * Assign GRAVITY2(K,J), which is written to .nc files.
   */
  for (J = JLO; J <= JHI; J++) {
    for (K = KLO; K <= KHI; K++) {
      GRAVITY2(K,J) = grid.g[2*K][2*J+1];
    }
    GRAVITY2(KLO-1,J) = grid.g[1      ][2*J+1];
    GRAVITY2(KHI+1,J) = grid.g[2*KHI+1][2*J+1];
  }

  return;
}

/*======================= end of set_gravity() ==============================*/

/*======================= set_dsgth() =======================================*/

void set_dsgth(void)
{
  register int
    K,kk;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_dsgth";

  /*
   * Calculate differential and its reciprocal for sigmatheta.
   */

  /* odd kk correspond to interface values */
  for (K = KLO-1; K <= KHI; K++) {
    kk = 2*K;
    if (K == 0) {
      grid.dsgth[kk+1] = grid.sigmatheta[kk+1]-grid.sigmatheta[kk+2];
    }
    else if (K == KHI) {
      grid.dsgth[kk+1] = grid.sigmatheta[kk  ]-grid.sigmatheta[kk+1];
    }
    else {
      grid.dsgth[kk+1] = grid.sigmatheta[kk  ]-grid.sigmatheta[kk+2];
    }

    if (grid.dsgth[kk+1] > 0.) {
      grid.dsgth_inv[kk+1] = 1./grid.dsgth[kk+1];
    }
    else {
      sprintf(Message,"grid.dsgth[kk=%d]=%g; check input T(p)",kk+1,grid.dsgth[kk+1]);
      epic_error(dbmsname,Message);
    }
  }

  /* even kk correspond to layer values */
  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;
    grid.dsgth[kk] = grid.sigmatheta[kk-1]-grid.sigmatheta[kk+1];

    if (grid.dsgth[kk] > 0.) {
      grid.dsgth_inv[kk] = 1./grid.dsgth[kk];
    }
    else{
      sprintf(Message,"grid.dsgth[kk=%d]=%g; check input T(p)",kk,grid.dsgth[kk]);
      epic_error(dbmsname,Message);
    }
  }

  return;
}

/*======================= end of set_dsgth() ================================*/

/*======================= get_sigma() =======================================*/
/*
 * Use a function to calculate sigma so that its definition can be 
 * easily changed.
 *
 * NOTE: If this function is changed, then it is necessary to also modify
 *       its inverse function, get_p_sigma().
 *
 * NOTE: We use log p instead of p because sigma defined with the latter stays
 *       close to 0 too long with respect to height when pbot is large, such as 
 *       for Venus or deep gas-giant models.
 */

double get_sigma(double pbot,
                 double p)
{
  return log(p/pbot)/log(grid.ptop/pbot);
}

/*======================= end of get_sigma() ================================*/

/*======================= get_p_sigma() =====================================*/
/*
 * Inverse of get_sigma() function.  If one is changed, the other should
 * be matched accordingly.
 */
double get_p_sigma(double pbot,
                   double sigma)
{
  return pbot*exp(sigma*log(grid.ptop/pbot));
}

/*======================= end of get_p_sigma() ==============================*/

/*======================= calc_h() ===========================================*/

/*
 * Calculates h = -(1/g) dp/dsgth.
 *
 * The input pressure, p, and the output h, are vertical columns that use
 * the doubled kk notation for which the top of the model is kk = 2*KLO-1 and
 * the bottom is kk = 2*KHI+1.
 */

void calc_h(int         jj,
            EPIC_FLOAT *p,
            EPIC_FLOAT *h)
{
  int
    K,kk;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="calc_h";

  for (kk = 2; kk <= 2*KHI; kk++) {
    h[kk] = (p[kk+1]-p[kk-1])*grid.dsgth_inv[kk]/grid.g[kk][jj];
  }
  h[1]       = h[2]*h[2]/h[3];
  h[2*KHI+1] = h[2*KHI]*h[2*KHI]/h[2*KHI-1];

  /*
   * Maintain h_min.
   */
  K = 0;
  h[2*K+1] = MAX(h[2*K+1],grid.h_min[K]);
  for (K = KLO; K <= KHI; K++) {
    h[2*K  ] = MAX(h[2*K  ],grid.h_min[K]);
    h[2*K+1] = MAX(h[2*K+1],sqrt(grid.h_min[K]*grid.h_min[K+1]));
  }

  return;
}

/*======================= end of calc_h() ====================================*/

/*======================= get_p() ============================================*/
/*
 * Return total or partial gas pressure, depending on the input index.
 */

EPIC_FLOAT get_p(planetspec *planet,
                 int         index,
                 int         kk,
                 int         J,
                 int         I) 
{
  register int
    K;
  register EPIC_FLOAT
    pressure;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_p";

#if EPIC_CHECK == 1
  /*
   * Check validity of kk:
   */
  if (kk < 1) {
    sprintf(Message,"kk = %d < 1",kk);
    epic_error(dbmsname,Message);
  }
  else if (kk > 2*(grid.nk+1)) {
    sprintf(Message,"kk = %d > 2*(nk+1) = %d",kk,2*(grid.nk+1));
    epic_error(dbmsname,Message);
  }
#endif

  if (kk%2 == 0) {
    /* 
     * Layer value.
     *
     * NOTE: We get better results if we calculate the layer p, P2(K),
     *       based only on P3(K) and P3(K-1), rather than mixing in some 
     *       p = p(sigmatheta,theta).
     */
    K = kk/2;
    pressure = onto_kk(planet,P2_INDEX,P3(K-1,J,I),P3(K,J,I),kk,J,I);
    if (index == P2_INDEX || index == P3_INDEX) {
      return pressure;
    }
    else if (index >= FIRST_SPECIES && index <= LAST_SPECIES) { 
      return pressure*.5*(X(index,VAPOR,K,J,I)+X(index,VAPOR,K-1,J,I));
    }
    else {
      sprintf(Message,"index = %d unknown",index);
      epic_error(dbmsname,Message);
    }
  }
  else {
    /* 
     * Interface value.
     */
    K = (kk-1)/2;
    pressure = P3(K,J,I);
    if (index == P2_INDEX || index == P3_INDEX) {
      return pressure;
    }
    else if (index >= FIRST_SPECIES && index <= LAST_SPECIES) {
      return pressure*X(index,VAPOR,K,J,I); 
    }
    else {
      sprintf(Message,"index = %d unknown",index);
      epic_error(dbmsname,Message);
    }
  }

  /* Should never get here.*/
  sprintf(Message,"should never get here");
  epic_error(dbmsname,Message);
  return 0.;
}

/*======================= end of get_p() =====================================*/

/*======================= molar_mixing_ratio() ===============================*/

/*
 * The molar mass for dry air is assumed to be R_GAS/planet->rgas.
 */

EPIC_FLOAT molar_mixing_ratio(planetspec *planet,
                              int         is,
                              int         ip,
                              int         kk,
                              int         J,
                              int         I)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="molar_mixing_ratio";

  /*
   * Return 0. if species and phase are not invoked.
   */
  if (!var.species[is].phase[ip].on) {
    return 0.;
  }

  return get_var(planet,is,ip,grid.it_h,kk,J,I)*R_GAS/
         (planet->rgas*var.species[is].molar_mass);
}

/*======================= end of molar_mixing_ratio() ========================*/

/*======================= get_var() ==========================================*/
/*
 * Returns the value of a prognostic variable at the point IT,kk,J,I
 * referenced by its index or species/phase index pair.  
 *
 * Error if called for a diagnostic variable.
 *
 * For prognostic variables that are not species, species_index is taken 
 * to be the index, and phase_index is ignored.
 *
 * To illustrate the use of the kk index, note that the layer value of layer
 * K = 3 corresponds to kk = 6, the lower-altitude interface value corresponds
 * to kk = 7, and the upper-altitude interface value corresponds to kk = 5.
 *
 * NOTE: Do not rely on diagnostic arrays like P2(K,J,I) or THETA2(J,K,I), 
 *       because they are not necessarily set prior to calling get_var(). 
 *
 * NOTE: High-level subroutines should call get_var() to evaluate a variable
 *       at a vertical position that is different than its natural position, 
 *       not onto_kk(), which is just an averager and does not handle 
 *       boundary cases.
 */

EPIC_FLOAT get_var(planetspec *planet,
                   int         species_index,
                   int         phase_index,
                   int         IT,
                   int         kk,
                   int         J,
                   int         I)
{
  int
    K;
  EPIC_FLOAT
    x,y,dy,
    xa[3],ya[3];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_var";

  if (species_index > MAX_NUM_PROGS-1) {
    /*
     * Error if called for a diagnostic variable.
     */
    sprintf(Message,"index=%d does not refer to a prognostic variable",species_index);
    epic_error(dbmsname,Message);
  }
  else if (species_index >= FIRST_SPECIES && species_index <= LAST_SPECIES) {
    /*
     * Return 0. if referring to a species or phase that is not on.
     */
    if (!var.species[species_index].on) {
      return 0.;
    }
    else if (!var.species[species_index].phase[phase_index].on) {
      return 0.;
    }
  }

  if (kk < 1) {
    K = 0;
    switch(species_index) {
      case U_INDEX:
        return U(IT,K,J,I);
      break;
      case V_INDEX:
        return V(IT,K,J,I);
      break;
      case H_INDEX:
        return H(K,J,I);
      break;
      case NU_TURB_INDEX:
        return NU_TURB(K,J,I);
      break;
      case FPARA_INDEX:
        if (var.fpara.on) {
          return FPARA(K,J,I);
        }
        else {
          return return_fpe(T2(K,J,I));
        }
      break;
      default:
        if (species_index >= FIRST_SPECIES && species_index <= LAST_SPECIES) {
          return Q(species_index,phase_index,K,J,I);
        }
        else {
          sprintf(Message,"need implementation for species_index=%d",species_index);
          epic_error(dbmsname,Message);
        }
      break;
    }
  }
  else if (kk > 2*grid.nk+1) {
    K = grid.nk+1;
    switch(species_index) {
      case U_INDEX:
        return U(IT,K,J,I);
      break;
      case V_INDEX:
        return V(IT,K,J,I);
      break;
      case H_INDEX:
        return H(K,J,I);
      break;
      case NU_TURB_INDEX:
        return NU_TURB(K,J,I);
      break;
      case FPARA_INDEX:
        if (var.fpara.on) {
          return FPARA(K,J,I);
        }
        else {
          return return_fpe(T2(K,J,I));
        }
      break;
      default:
        if (species_index >= FIRST_SPECIES && species_index <= LAST_SPECIES) {
          return Q(species_index,phase_index,K,J,I);
        }
        else {
          sprintf(Message,"need implementation for species_index=%d",species_index);
          epic_error(dbmsname,Message);
        }
      break;
    }
  }
  else if (kk%2 == 0) {
    /*
     * Layer value.
     */
    K = kk/2;

   /*
    * Handle pass-through cases.
    */
    switch(species_index) {
      case U_INDEX:
        return U(IT,K,J,I);
      break;
      case V_INDEX:
        return V(IT,K,J,I);
      break;
      case H_INDEX:
        return H(K,J,I);
      break;
      case NU_TURB_INDEX:
        return NU_TURB(K,J,I);
      break;
    }

    return onto_kk(planet,species_index,get_var(planet,species_index,phase_index,IT,kk-1,J,I),
                                        get_var(planet,species_index,phase_index,IT,kk+1,J,I),kk,J,I);
  }
  else {
    /*
     * Interface value.
     */
    K = (kk-1)/2;

    /*
     * Handle pass-through cases.
     */
    switch (species_index) {      
      case THETA_INDEX:
        return THETA(K,J,I);
      break;
      case FPARA_INDEX:
        if (var.fpara.on) {
          return FPARA(K,J,I);
        }
        else {
          return return_fpe(T3(K,J,I));
        }
      break;
      default:
        if (species_index >= FIRST_SPECIES && species_index <= LAST_SPECIES) {
          return Q(species_index,phase_index,K,J,I);
        }
      break;
    }

    return onto_kk(planet,species_index,get_var(planet,species_index,phase_index,IT,kk-1,J,I),
                                        get_var(planet,species_index,phase_index,IT,kk+1,J,I),kk,J,I);
  }

  /* Should not get here. */
  sprintf(Message,"should not get to end of this function");
  epic_error(dbmsname,Message);
  return FLT_MAX;
}

/*======================= end of get_var() ===================================*/

/*======================= get_var_mean2d() ===================================*/

/*
 * Return the global mean of a 2 dimensional layer variable.
 *
 * Aaron Herrnstein, January 2006.
 */

EPIC_FLOAT get_var_mean2d(EPIC_FLOAT *a,
                          int         index,
                          int         kk)
{
  register int
    I,J,jbot,jay;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    *da_u,
    *da_v,
     sum_da_u,
     sum_da_v;
  EPIC_FLOAT
    *da,
     sum_da,
     var_mean;
#if defined(EPIC_MPI)
  EPIC_FLOAT
    mpi_tmp;
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_var_mean2d";

  if (!initialized) {
    /* Allocate memory. */
    da_u = fvector(0,JHI-JLO,  dbmsname);
    da_v = fvector(0,JHI+1-JLO,dbmsname);

    sum_da_u = 0.;
    for (J = JLO; J <= JHI; J++) { 
      da_u[J-JLO] = 1./grid.mn[kk][2*J+1];
      sum_da_u   += da_u[J-JLO];
    }
    sum_da_u *= grid.ni;
#if defined(EPIC_MPI)
    mpi_tmp = sum_da_u;
    MPI_Allreduce(&mpi_tmp,&sum_da_u,1,float_type,MPI_SUM,para.comm);
#endif

    sum_da_v = 0.;
    for (J = JFIRST; J <= JHI; J++) { 
      da_v[J-JLO] = 1./grid.mn[kk][2*J];
      sum_da_v   += da_v[J-JLO]; 
    }
    if (JLO == grid.jlo) {
      J = JLO;
      da_v[J-JLO] = 1./grid.mn[kk][2*J];
      sum_da_v   += da_v[J-JLO];
    }
    if (JHI == grid.nj) {
      J = grid.nj+1;
      da_v[J-JLO] = 1./grid.mn[kk][2*J];
      sum_da_v   += da_v[J-JLO];
    }
    sum_da_v *= grid.ni;
#if defined(EPIC_MPI)
    mpi_tmp = sum_da_v;
    MPI_Allreduce(&mpi_tmp,&sum_da_v,1,float_type,MPI_SUM,para.comm);
#endif

    initialized = TRUE;
  }

  if (index == V_INDEX || index == PV2_INDEX) {
    jbot   = JFIRST;
    da     = da_v;
    sum_da = sum_da_v;
  } 
  else {
    jbot   = JLO;
    da     = da_u;
    sum_da = sum_da_u;
  }
  
  var_mean = 0.0;  
  for (J = jbot; J <= JHI; J++) {
    jay = J-JLO;
    for (I = ILO; I <= IHI; I++) {
      var_mean += A(J,I)*da[jay];
    }
  }
  /* No need for bc_lateral() here. */

  if (index == PV2_INDEX) {
    /* Include poles or channel edges. */
    if (JLO == grid.jlo) {
      jay = grid.jlo-JLO;
      for (I = ILO; I <= IHI; I++) {
        var_mean += A(J,I)*da[jay];
      }
    }
    if (JHI == grid.nj) {
      jay = grid.nj+1-JLO;
      for (I = ILO; I <= IHI; I++) {
        var_mean += A(J,I)*da[jay];
      }
    }
  }

#if defined(EPIC_MPI)
  mpi_tmp = var_mean;
  MPI_Allreduce(&mpi_tmp,&var_mean,1,float_type,MPI_SUM,para.comm);
#endif
  var_mean /= sum_da;
  
  return var_mean;
}

/*==================== end of get_var_mean2d() ===============================*/

/*======================= onto_kk() ==========================================*/

/*
 * NOTE: High-level subroutines should call get_var() to evaluate a variable
 *       at a vertical position that is different than its natural position, 
 *       not onto_kk(), which is just an averager and does not handle 
 *       boundary cases.
 *
 * Evaluate variable onto vertical position kk, the doubled 
 * vertical index.  Having all the vertical averaging schemes together
 * in one place makes it easier to manage them.
 *
 * An illustration of the kk index convention is:
 *   kk = 3 => K = 1.5, the bottom of layer 1, 
 *   kk = 4 => K = 2.0, the middle of layer 2, etc.
 *
 * The input parameters topval and botval are the higher-altitude and
 * lower-altitude values to be averaged onto level kk, in other words the
 * lower K and higher K endpoints, respectively.
 */

EPIC_FLOAT onto_kk(planetspec *planet,
                   int         index,
                   EPIC_FLOAT  topval,
                   EPIC_FLOAT  botval,
                   int         kk,
                   int         J,
                   int         I)
{
  register EPIC_FLOAT
    kappap1,
    topwt,botwt;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="onto_kk";

  /*
   * Check validity of kk.
   */
  if (kk < 1) {
    sprintf(Message,"kk=%d < 1",kk);
    epic_error(dbmsname,Message);
  }
  else if (kk > 2*(grid.nk+1)) {
    sprintf(Message,"kk = %d > 2*nk+1=%d",kk,2*grid.nk+1);
    epic_error(dbmsname,Message);
  }

  if (kk%2 == 0) {
    /* Layer value. */

    switch(index) {
      case U_INDEX:
      case V_INDEX:
      case H_INDEX:
        sprintf(Message,"kk=%d, index=%d carried in layer, use directly",kk,index);
        epic_error(dbmsname,Message);
      break;
      case P2_INDEX:
        /*
         * NOTE: For many years we used the complicated inside-the-layer definition of pressure 
         *       of Hsu and Arakawa 1990 (5.43) (also Konor and Arakawa (3.34)), which we preserve here
         *       commented out.  However, we have switched to a simple geometric average, because
         *       it produces smoother results.
         */
        /**** Hsu and Arakawa (1990) 
        if (fcmp(botval,topval) == 0) {
          return .5*(botval+topval);
        }
        else if (topval < 0. || botval < 0.) {
          sprintf(Message,"index=P2_INDEX, kk=%d, topval=%g, botval=%g",
                           kk,topval,botval);
          epic_error(dbmsname,Message);
        }
        else {
          kappap1 = planet->kappa+1.;
          return pow( (pow(botval,kappap1)-pow(topval,kappap1))/
                          (kappap1*(botval-topval)),1./planet->kappa);
        }
        ******/
        return sqrt(botval*topval);
      break;
      case RHO2_INDEX:
      case HDRY2_INDEX:
        /*
         * Use a geometric average for variables that vary exponentially with height.
         */
        return sqrt(botval*topval);
      break;
      case THETA2_INDEX:
        /*
         * Use plain averaging as in Konor and Arakawa (1997, eqn 3.26).
         */
        return .5*(botval+topval);
      break;
      default:
        /*
         * By default, use delta-sigmatheta weighting.
         */
        topwt = grid.dsgth[kk-1];
        botwt = grid.dsgth[kk+1];
        return (topwt*topval+botwt*botval)/(topwt+botwt);
      break;
    } /* end switch */
  }
  else {
    /* Interface value. */

    switch(index) {
      case P3_INDEX:
      case THETA_INDEX:
      case FPARA_INDEX:
      case W3_INDEX:
        sprintf(Message,"kk=%d, index=%d defined on interface, use directly",kk,index);
        epic_error(dbmsname,Message);
      break;
      case U_INDEX:
        return .5*(botval+topval);
      break;
      case V_INDEX:
        return .5*(botval+topval);
      break;
      default:
        /*
         * By default, use delta-sigmatheta weighting.
         */
        if (kk == 1) {
          topwt = 0.;
          botwt = 1.;
        }
        else if (kk == 2*grid.nk+1) {
          topwt = 1.;
          botwt = 0.;
        }
        else {
          topwt = grid.dsgth[kk-1];
          botwt = grid.dsgth[kk+1];
        }
        return (topwt*topval+botwt*botval)/(topwt+botwt);
      break;
    }  /* end switch */
  }

  /* Should never get here.*/
  sprintf(Message,"should never get here");
  epic_error(dbmsname,Message);
  return 0.;
}

/*======================= end of onto_kk() ===================================*/

/*====================== fpe_minus_fpe() =====================================*/

EPIC_FLOAT fpe_minus_fpe(EPIC_FLOAT fpe)
{
  EPIC_FLOAT
    p,theta,temp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="fpe_minus_fpe";

  p     = FPEMFPE_p;
  theta = FPEMFPE_theta;
  temp  = return_temp(planet,fpe,p,theta);

  return fpe-return_fpe(temp);
}

/*====================== fpe_minus_fpe() =====================================*/

/*======================= get_kin() ==========================================*/

/*
 * Return kinetic energy per mass. 
 * Kinetic energy resides on the h-grid.
 * See comments on Hollingsworth-Kallberg non-cancellation instability below.
 */

EPIC_FLOAT get_kin(planetspec *planet,
                   EPIC_FLOAT *u2d,
                   EPIC_FLOAT *v2d,
                   int         kk,
                   int         J,
                   int         I)
{
  register int
    jj;
  static int
    j_periodic  = FALSE,
    initialized = FALSE;
  register EPIC_FLOAT
    kin,kin_c,kin_s,
    u2,u4,v2,v4,
    alpha;
  static EPIC_FLOAT
    *mn_inv;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_kin";

  if (!initialized) {
    /* Allocate memory. */
    mn_inv = fvector(2*grid.jlo,2*(grid.nj+1),dbmsname);

    for (jj = 2*grid.jlo; jj <= 2*(grid.nj+1); jj++) {
      mn_inv[jj] = 1./grid.mn[kk][jj];
    }

    if (strcmp(grid.geometry,"f-plane") == 0 &&
        strcmp(grid.f_plane_map,"cartesian") == 0) {
      j_periodic = TRUE;
    }

    initialized = TRUE;
  }

  /*
   * Check validity of J, I.
   */
  if (J < JLO || J > JHI) {
    sprintf(Message,"J=%d out of range [%d,%d]",J,JLO,JHI);
    epic_error(dbmsname,Message);
  }
  if (I < ILO || I > IHI) {
    sprintf(Message,"I=%d out of range [%d,%d]",I,ILO,IHI);
    epic_error(dbmsname,Message);
  }

  /*
   * The Hollingsworth-Kallberg instability arises when there is incomplete 
   * cancellation of terms in the discrete horizontal momentum equations when 
   * written in vector-invarient form.  See
   *
   *   Hollingsworth A, Kallberg P, Renner V, Burridge DM, 1983, An internal symmetric
   *       computational instability, Quart. J. R. Met. Soc., 109, 417--428.
   *
   * Basically, the term -u*du/dy from (pv)*uh on the left-hand side of the dvdt
   * equation does not exactly cancel the -d(u*u/2)/dy term from dK/dy on the
   * right-hand side, and the same for the (pv)*vh and dK/dx terms in the dudt
   * equation.
   *
   * One remedy is to modify the kinetic energy calculation, specifically to extend
   * its stencil to have the same reach as the Arakawa C-grid's averaging of the Coriolis terms.
   * This was done for example by
   *
   *   Suarez MJ, Takacs LL, 1995, NASA Technical Memorandum 104606, Vol. 5.
   * 
   * They define K_s to be a kinetic energy with a longer stencil, which augments
   * the original K_c, where the 'c' refers energy conservation.
   * Then, K = alpha*K_c + (1.-alpha)*K_s.
   * Suarez and Takacs analyzed their algorithm with a linearized approach and determined
   * alpha = 5/6. However, their scheme is different than the two algorithms available in EPIC, 
   * PV_SCHEME = SADOURNY_1975 and PV_SCHEME = ARAKAWA_LAMB_1981, hence we determined 
   * alpha for each of our cases. 
   *
   * We did this empirically, as follows.  We modified a copy of EPIC to turn off all physics
   * except uv_core(), set H = 1, V = 0, and reduced to a 2D model with ni = 1, thus isolating
   * the meridional terms in question. To handle map factors, we monitored
   *
   *   TE = zeta*u*dy + (1/r^2) d(r^2 K) ,
   *
   * where the r^2 weighting on K takes into account the u*u*dr/r term.
   *
   * The interior points and the boundary points are different, because K_s as defined for the interior
   * points extends past the poles if applied on the boundary rows (J = 0 and nj for the h-grid).
   * We experimented with different zonal-wind profiles, u, including an actual profile from 
   * a Uranus spin-up experiment, and idealed profiles.  A linear profile of u vs lat does not
   * vary enough, but a quadratic u yields results similar to a realistic u profile.
   * For interior points, our empirically determined results are 
   *
   *   ___PV_SCHEME___      ___alpha___
   *
   *    ARAKAWA_LAMB_1981       1/2
   *    SADOURNY_1975           2/3
   *
   * The quality of the cancellation does depend on the curvature of u, however the optimal alpha
   * value does not depend on the magnitude of u (the truncation error does scale as u^2).
   * For the pole rows, we modify K_s by just setting u=0 for the two points that would extend
   * past the poles (the other two points extend the stencil as needed, since the zeta term at the
   * poles is computed via the circulation theorem).  We find that the same values for alpha work
   * reasonably well.
   *
   */

#if PV_SCHEME == ARAKAWA_LAMB_1981
  /*
   * Arakawa and Lamb (1981)
   */
  alpha = 1./2.;
#elif PV_SCHEME == SADOURNY_1975
  alpha = 2./3.;
#else
  sprintf(Message,"unrecognized PV_SCHEME.")
  epic_error(dbmsname,Message);
#endif

  jj    = 2*J;
  u2    = U2D(J,  I  );
  u4    = U2D(J,  I+1);
  v2    = V2D(J,  I  );
  v4    = V2D(J+1,I  );
  kin_c = (                 u2*u2
                           +u4*u4
            +grid.mn[kk][jj+1]*(v2*v2*mn_inv[jj  ]
                           +v4*v4*mn_inv[jj+2]) )*.25;

  if (!j_periodic && (J == grid.jlo  || J == grid.nj)) {
    /*
     * Set u values that are past pole position to zero for K_s.
     */
    if (J == grid.nj){
      u2 = .5*(0.+U2D(J-1,I  ));
      u4 = .5*(0.+U2D(J-1,I+1));
    } 
    else {
      u2 = .5*(U2D(J+1,I  )+0.);
      u4 = .5*(U2D(J+1,I+1)+0.);
    }
  }
  else {
    u2 = .5*(U2D(J+1,I  )+U2D(J-1,I  ));
    u4 = .5*(U2D(J+1,I+1)+U2D(J-1,I+1));
  }

  v2    = .5*(V2D(J,  I+1)+V2D(J,  I-1));
  v4    = .5*(V2D(J+1,I+1)+V2D(J+1,I-1));
  kin_s = (                     u2*u2
                               +u4*u4
            +grid.mn[kk][jj+1]*(v2*v2*mn_inv[jj  ]
                               +v4*v4*mn_inv[jj+2]) )*.25;

  kin   = alpha*kin_c+(1.-alpha)*kin_s;

#if EPIC_CHECK == 1
  /*
   * Screen for nan.
   */
  if (!isfinite(kin)) {
    sprintf(Message,"JI=%d %d, kin=%g, u2=%g u4=%g v2=%g v4=%g kin_s=%g kin_c=%g; PV_SCHEME=SADOURNY_1975",
                    J,I,kin,u2,u4,v2,v4,kin_s,kin_c);
    epic_error(dbmsname,Message);
  }
#endif

  return kin;
}

/*======================= end of get_kin() ===================================*/

/*======================= get_brunt2() =======================================*/

/*
 * A.P. Showman, 8/31/99.
 * See notes dated 8/31/99.
 *
 * Option of smooth derivative of Drho_Dp added by T. Dowling, 11/03/05.
 *
 * Calculates and returns the squared Brunt-Vaisala (buoyancy) frequency at
 * position kk/2,J,I.  The formula used holds for any equation of state. It
 * incorporates a dry adiabatic lapse rate assuming no chemical reactions or
 * condensation. The environmental density structure takes into account the 
 * effects of vertical gradients of molar mass and entropy as well as 
 * compressibility. 
 *
 * The equation used is
 *
 *   N^2 = g^2{-[drho/dp]_T+(T/rho^2 cp)([drho/dT]_p)^2+Drho/Dp}      (1)
 *
 * where [drho/dx]_y is the partial derivative of rho with respect to x 
 * at const y, and Drho/Dp is the total derivative of rho along the 
 * environmental profile. This form is more accurate than 
 *
 *   N^2 = (g/theta)dtheta/dz                                         (2)
 *
 * which assumes the ideal gas law with no molar mass gradients and only
 * works for the traditional definition of theta.  EPIC uses a more 
 * general mean-theta for hydrogen (ortho and para) such that (2) does
 * not yield the correct value for N^2.
 */

#undef C
#define C(i,j) c[j+i*(np)]

EPIC_FLOAT get_brunt2(planetspec *planet,
                      int         kk,
                      int         J,
                      int         I)
{
  register int
    K,
    kayk,
    nl,nr,nn;
  static int
    np,nhw,
    initialized = FALSE;
  static EPIC_FLOAT
    *c,
    *rho,
    *press;
  EPIC_FLOAT
    brunt2,        /* squared Brunt-Vaisala frequency, 1/s^2      */
    mu,            /* molar mass                                  */
    pressure,
    temperature,
    density,
    fpara,
    deltap,
    deltaT,
    drho_dp_T,     /* partial deriv of rho w/r to p at const T       */
    drho_dT_p,     /* partial deriv of rho w/r to T at const p       */
    Drho,Dp,
    Drho_Dp,       /* total deriv of environmental rho profile wrt p */
    cp,            /* specific heat at constant pressure             */
    g;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_brunt2";

  if (!initialized) {
    rho   = fvector(0,2*KHI+1,dbmsname);
    press = fvector(0,2*KHI+1,dbmsname);

    /*
     * Set half-width, nhw, for Savitsky-Golay smoothing of Drho_Dp.
     */
    /****Currently setting np = 1, which turns off smoothing.
    np  = 2*(grid.nk/30)+1;
    ****/
    np  = 1;

    nhw = (np-1)/2;

    if (np > 3) {
      c   = fvector(0,np*np-1,dbmsname);
      /*
       * Calculate Savitsky-Golay coefficients for first derivative.
       * NOTE: the second-to-last argument in savitzky_golay() indicates
       *       the order of the derivative, which is 1 in this case.
       */
      for (nl = 0; nl < np; nl++) {
        nr = np-nl-1;
        savitzky_golay(c+nl*np,np,nl,nr,1,2);
      }
    }

    initialized = TRUE;
  }

  K  = kk/2;

  if (kk%2 == 0) {
    /* 
     * Get values in layer:
     */
    pressure    = P2(  K,J,I);  
    temperature = T2(  K,J,I);  
    density     = RHO2(K,J,I);  
  }
  else {
    /* 
     * Get values at interface:
     */
    pressure    = P3(  K,J,I);  
    temperature = T3(  K,J,I);  
    density     = RHO3(K,J,I);  
  }
  if (var.fpara.on) {
    fpara = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,kk,J,I);
  }
  else {
    fpara = return_fpe(temperature);
  }
  cp    = return_cp(planet,fpara,pressure,temperature);
  mu    = avg_molar_mass(planet,kk,J,I);

  deltap = 0.001*pressure;
  deltaT = 0.001*temperature;

  drho_dp_T = (return_density(planet,fpara,pressure+deltap,temperature,mu,PASSING_T) 
              -return_density(planet,fpara,pressure-deltap,temperature,mu,PASSING_T))/
              (2.*deltap);

  drho_dT_p = (return_density(planet,fpara,pressure,temperature+deltaT,mu,PASSING_T)
              -return_density(planet,fpara,pressure,temperature-deltaT,mu,PASSING_T))/
              (2.*deltaT);

  /*
   * Calculate Drho_Dp.
   */
  if (np <= 3) {
    /*
     * Do regular differencing without smoothing.
     */
    if (kk%2 == 0) {
      Drho_Dp = (RHO3(K,J,I)-RHO3(K-1,J,I))/
                (P3(  K,J,I)-P3(  K-1,J,I));
    }
    else {
      if (K == KHI) {
        Drho_Dp = (RHO3(K,J,I)-RHO2(K,J,I))/
                  (P3(  K,J,I)-P2(  K,J,I));
      }
      else if (K == KLO-1) {
        Drho_Dp = (RHO2(K+1,J,I)-RHO3(K,J,I))/
                  (P2(  K+1,J,I)-P3(  K,J,I));
      }
      else {
         Drho_Dp = (RHO2(K+1,J,I)-RHO2(K,J,I))/
                   (P2(  K+1,J,I)-P2(  K,J,I));
      }
    }
  }
  else {
    /*
     * Use Savitzky-Golay smoothing of the Drho/Dp derivative
     * to control computational mode that may arise in low-N^2 regions.
     */
    if (kk <= nhw) {
      nl = kk-1;
      nr = np-nl-1;
    }
    else if (kk >= 2*KHI+2-nhw) {
      nr = 2*KHI+1-kk;
      nl = np-nr-1;
    }
    else {
      nl = nr = nhw;
    }

    for (K = (kk-nl)/2; K <= (kk+nr)/2; K++) {
      kayk = 2*K;
      rho[  kayk  ] = RHO2(K,J,I);
      rho[  kayk+1] = RHO3(K,J,I);
      press[kayk  ] = P2(  K,J,I);
      press[kayk+1] = P3(  K,J,I);
    }

    Drho = C(nl,0)*rho[  kk];
    Dp   = C(nl,0)*press[kk];
    for (nn = 1; nn <= nl; nn++) {
      Drho += C(nl,nn)*rho[  kk-nn];
      Dp   += C(nl,nn)*press[kk-nn];
    }
    for (nn = 1; nn <= nr; nn++) {
      Drho += C(nl,np-nn)*rho[kk+nn];
      Dp   += C(nl,np-nn)*press[kk+nn];
    }
    Drho_Dp = Drho/Dp;
  }

  g = grid.g[kk][2*J+1];

  brunt2 = g*g*(Drho_Dp-drho_dp_T
                +temperature/(cp*density*density)*drho_dT_p*drho_dT_p);

#if EPIC_CHECK == 1
  if (!isfinite(brunt2)) {
    sprintf(Message,"brunt2=%g, temperature=%g, pressure=%g, density=%g, cp=%g, mu=%g",
                     brunt2,temperature,pressure,density,cp,mu);
    epic_error(dbmsname,Message);
  }
#endif

  return brunt2;
}

/*======================= end of get_brunt2() ================================*/

/*====================== get_richardson() ====================================*/
        
/*
 * Returns Richardson number, Ri = N^2/(du/dz)^2.
 * Calculates Ri on the h-grid or p3-grid.
 */

/*
 * Set a large but finite positive limit for the amplitude of the answer.
 */
#define MAX_RI (1.e+10)
 
EPIC_FLOAT get_richardson(planetspec *planet,
                          int         kk,
                          int         J,
                          int         I)
{
  int
    K;
  EPIC_FLOAT
    dudz2,
    dvdz2,
    dveldz2;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_richardson";

  if (kk < 1 || kk > grid.nk*2+1) {
    /*
     * kk is out of bounds.
     */
    sprintf(Message,"kk=%d out of range",kk);
    epic_error(dbmsname,Message);
  }
  else if (kk == 1) {
    /*
     * Top of model.
     */
    return get_richardson(planet,kk+1,J,I);
  }
  else if (kk == 2*grid.nk+1) {
    /*
     * Bottom of model.
     */
    return get_richardson(planet,kk-1,J,I);
  }
  else if (kk%2 == 0) {
    /*
     * Layer value.
     */
    K = kk/2;

    dudz2 = .5*(U(grid.it_uv,K-1,J,I)+U(grid.it_uv,K-1,J,I+1)-U(grid.it_uv,K+1,J,I)-U(grid.it_uv,K+1,J,I+1))
              /(Z2(K-1,J,I)-Z2(K+1,J,I));

    dvdz2 = .5*(V(grid.it_uv,K-1,J,I)+V(grid.it_uv,K-1,J+1,I)-V(grid.it_uv,K+1,J,I)-V(grid.it_uv,K+1,J+1,I))
              /(Z2(K-1,J,I)-Z2(K+1,J,I));
  }
  else {
    /*
     * Interface value.
     */
    K     = (kk-1)/2;

    dudz2 = .5*(U(grid.it_uv,K,J,I)+U(grid.it_uv,K,J,I+1)-U(grid.it_uv,K+1,J,I)-U(grid.it_uv,K+1,J,I+1))
              /(Z2(K,J,I)-Z2(K+1,J,I));

    dvdz2 = .5*(V(grid.it_uv,K,J,I)+V(grid.it_uv,K,J+1,I)-V(grid.it_uv,K+1,J,I)-V(grid.it_uv,K+1,J+1,I))
              /(Z2(K,J,I)-Z2(K+1,J,I));
  }

  dudz2  *= dudz2;
  dvdz2  *= dvdz2;
  dveldz2 = dudz2+dvdz2;

  if (fcmp(dveldz2,0.) != 0) {
    /*
     * Limit answer to signed MAX_RI value.
     */
    return LIMIT_RANGE(-MAX_RI,get_brunt2(planet,kk,J,I)/dveldz2,MAX_RI);
  }
  else {
    /*
     * Return signed MAX_RI for Ri = infinity case.
     */
    return NR_SIGN(MAX_RI,get_brunt2(planet,kk,J,I));
  }
}

#undef MAX_RI

/*====================== end of get_richardson() =============================*/

/*======================= set_p2_etc() =======================================*/
/*
 * The diagnostic variables P2, P3, PDRY3, HDRY2, HDRY3, H3, and THETA2 are updated
 * here rather than in store_diag(), because they are directly 
 * related to the prognostic variables H and THETA and need to be updated
 * more often than the other diagnostic variables.
 *
 * NOTE: Currently not using argument Buff2D.
 */

void set_p2_etc(planetspec  *planet,
                int          theta_switch,
                EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    kk,jj,iq;
  register EPIC_FLOAT
    sigma,
    sum,
    pbot;
  EPIC_FLOAT
    tmp,dp,p_diag,
    x,gsg;
  static int
    initialized   = FALSE,
    theta2_is_set = FALSE;
  static EPIC_FLOAT
    *h,*p;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_p2_etc";

  if (!initialized) {
    p = fvector(0,2*grid.nk+1,dbmsname);
    h = fvector(0,2*grid.nk+1,dbmsname);

    initialized = TRUE;
  }

  switch(grid.coord_type) {
    case COORD_ISOBARIC:
      ;
    break;
    case COORD_ISENTROPIC:
    case COORD_HYBRID:
      /*
       * Ensure nonnegative H.
       */
      restore_mass(planet,H_INDEX,NO_PHASE);
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  if (grid.coord_type != COORD_ISOBARIC) {
    /*
     * Integrate H to get P3.
     */
    for (K = KLO; K <= KHI; K++) {
      kk   = 2*K;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        jj = 2*J+1;
        for (I = ILOPAD; I <= IHIPAD; I++) {
          P3(K,J,I) = P3(K-1,J,I)+grid.g[kk][jj]*H(K,J,I)*grid.dsgth[kk];
        }
      }
    }
    /* No need to call bc_lateral() here. */

    K = KHIPAD;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        P3(K,J,I) = SQR(P3(K-1,J,I))/P3(K-2,J,I);
      }
    }
    /* No need to call bc_lateral() here. */
  }

  for (J = JLOPAD; J <=JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /*
       * Calculate P2 from P3.
       */
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K;
        /*
         * P2 depends on P3.
         */
        P2(K,J,I) = get_p(planet,P2_INDEX,kk,J,I);
      }

      /*
       * These extrapolations should be cosmetic only.
       */
      K = KLOPAD;
      P2(K,J,I) = SQR(P3(K,J,I))/P2(K+1,J,I);

      K = KHIPAD;
      P2(K,J,I) = SQR(P3(K-1,J,I))/P2(K-1,J,I);

      /*
       * Calculate H3 from H.
       *
       * NOTE: If this algorithm is changed, it should also be changed in advection().
       */
      for (K = KLO; K < KHI; K++) {
        H3(K,J,I) = sqrt(H(K,J,I)*H(K+1,J,I));
      }
      K = KLO-1;
      H3(K,J,I) = SQR(H(K+1,J,I))/H3(K+1,J,I);
      K = KHI;
      H3(K,J,I) = SQR(H(K,J,I))/H3(K-1,J,I);
      K = KHI+1;
      H3(K,J,I) = SQR(H3(K-1,J,I))/H3(K-2,J,I);
    }
  }
  /* No need to call bc_lateral() here. */

  if (grid.nq > 0) {
    /*
     * Calculate HDRY3, PDRY3 (dry hybrid density, pressure) from H3, P3 and Q_i (mass mixing ratios,
     * which are carried on the interfaces).
     *
     * NOTE: The grid.coord == COORD_ISOBARIC case with grid.nq > 0 needs to calculate PDRY3, too,
     *       since P3 is the coordinate but PDRY3 is not. 
     */
    for (K = KLO-1; K <= KHI; K++) {
      kk = 2*K+1;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          sum = 1.;
          for (iq = 0; iq < grid.nq; iq++) {
            sum += Q(grid.is[iq],grid.ip[iq],K,J,I);
          }
          HDRY3(K,J,I) = H3(K,J,I)/sum;

          PDRY3(K,J,I) = P3(K,J,I);
          for (iq = 0; iq < grid.nq; iq++) {
            PDRY3(K,J,I) -= get_p(planet,grid.is[iq],kk,J,I);
          }
        }
      }
    }

    /*
     * Calculate HDRY2 from HDRY3.
     */
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        for (K = KLO; K <= KHI; K++) {
          HDRY2(K,J,I) = sqrt(HDRY3(K,J,I)*HDRY3(K-1,J,I));
        }
        K = KLO-1;
        HDRY2(K,J,I) = SQR(HDRY3(K,J,I))/HDRY2(K+1,J,I);
        K = KHI+1;
        HDRY2(K,J,I) = SQR(HDRY3(K-1,J,I))/HDRY2(K-1,J,I);
      }
    }
  }

  if (!theta2_is_set) {
    /*
     * Calculate THETA2.
     */
    if (theta_switch == UPDATE_THETA) {
      for (J = JLOPAD; J <=JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          K = KLO-1;
          THETA2(K,J,I) = THETA(K,J,I);
          for (K = KLO; K <= KHI+1; K++) {
            THETA2(K,J,I) = .5*(THETA(K,J,I)+THETA(K-1,J,I));
          }
        }
      }
      /* No need to apply bc_lateral() here. */
    }

    if (grid.coord_type == COORD_ISENTROPIC) {
      theta2_is_set = TRUE;
    }
  }

  return;
}

/*======================= end of set_p2_etc() ================================*/

/*==================== store_pgrad_vars() ====================================*/

/*
 * Calculates and stores diagnostic variables that are typically used in the
 * pressure-gradient terms of the momentum equations.  Doing this in a
 * separate function facilitates the "poor man's implicit" timestep.
 *
 * An illustration of the notation with respect to the vertical index K:
 *
 *              ---------------------------------------------------
 *     Layer K:   T2  RHO2  EXNER2  FGIBB2  PHI2  MONT2  Z2 DZDT2(*)      
 * Interface K: --T3--RHO3--EXNER3----------PHI3---------Z3--------
 *
 * (*) For action == STEP_PROGS_AND_SYNC_DIAGS, this is a convenient place
 *     to start the calculation of the regular vertical velocity, DZDT2, in m/s.
 *     This variable is finished in calc_w(), after the hybrid vertical velocity
 *     is calculated.
 *
 * If phi3nk_status == PASSING_PHI3NK, PHI3(KHI,J,I) is assumed to be set
 * before calling this function; if CALC_PHI3NK, it is calculated here.
 *
 * NOTE: This function should not modify any prognostic variables. 
 *       In addition, it should not modify P2, P3, or THETA2.
 *
 * NOTE: This function must be called from all nodes.
 */

void store_pgrad_vars(planetspec  *planet,
                      EPIC_FLOAT **Buff2D,
                      int          action,
                      int          phi3nk_status) 
{
  register int
    K,J,I;
  int
    K0,J0,I0,
    kstart,kend,klen,
    is,ip,
    kk,kay,nkay,
    jj,order;
  register EPIC_FLOAT
    theta,
    fpara,
    pressure,
    temperature,
    enthalpy,
    mu,g_inv,
    avg,pbot,sigma,
    dt_inv,
    phi,sgth,
    z0,z1,z2,
    dtheta,mont3,mont0,phi0,
    small = 1.e-6;
  double
    nu;
  EPIC_FLOAT
    fgibb,fpe,uoup,
    fptol,
    dsgth,dz,tmp,g0_inv,
    nu_nondim,
    interval;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    *u1d,*phi1d,*mont1d,*a;
  float_triplet
    table[KHI-KLO+2];
#if defined(EPIC_MPI)
  EPIC_FLOAT
    mpi_tmp;
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="store_pgrad_vars";

  if (!initialized) {
    /* Allocate memory */
    a = fvector(0,KHI,dbmsname);

    if (strcmp(planet->type,"gas-giant") == 0) {
      /*
       * Calculate bottom geopotential based on gradient balance with
       * bottom wind profile.
       */

      /* Allocate memory */
      u1d = fvector(0,JADIM-1,dbmsname);

      for (J = JLOPAD; J <= JHIPAD; J++) {
        /*
         * Calculate zonal averages.
         */
        U1D(J) = 0.;
        for (I = ILO; I <= IHI; I++) {
          U1D(J) += U(grid.it_uv,KHI+1,J,I);
        }

#if defined(EPIC_MPI)
        mpi_tmp = U1D(J);
        MPI_Allreduce(&mpi_tmp,&U1D(J),1,float_type,MPI_SUM,para.comm_JLO);
#endif

        U1D(J) /= grid.ni;
      }

      /*
       * NOTE: Assuming the abyssal zonal wind is steady and barotropic,
       *       starting at the bottom of layer KHI.
       *
       * NOTE: Need to fix phi0 such that Z = 0 is close to p = 1000 hPa.
       */
      phi0 = 0.;

      /* Allocate memory */
      if (grid.coord_type == COORD_ISENTROPIC) {
        mont1d = fvector(0,JADIM-1,dbmsname);
        K      = KHI;
        mont0  = phi0+planet->cp*grid.t_ref[2*K+1];
        mont_from_u(planet,2*KHI+1,u1d,mont1d,grid.jtp,mont0);
      }
      else {
        /*
         * NOTE: Need to fix phi0 such that Z = 0 is close to p = 1000 hPa.
         */
        phi1d = fvector(0,JADIM-1,dbmsname);
        phi_from_u(planet,2*KHI+1,u1d,phi1d,grid.jtp,0.);
      }
    }

    initialized = TRUE;
  }

  if (action == STEP_PROGS_AND_SYNC_DIAGS) {
    /*
     * Store old value of Z2 in DZDT2, as the first
     * step to calculating the partial derivative of Z2 with respect to time.
     */
    for (K = KLO; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          DZDT2(K,J,I) = Z2(K,J,I);
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  /*
   * Calculate T3, EXNER3, RHO3.
   *
   * Global variable used to communicate with fpe_minus_fpe():
   */
  FPEMFPE_planet = planet;

  fptol = pow(machine_epsilon(),2./3.);

  if (var.fpara.on) {
    for (K = 0; K <= KHI; K++) {
      kk = 2*K+1;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          theta         = THETA(K,J,I);  
          pressure      = P3(K,J,I);

          fpara         = FPARA(K,J,I);

          T3(K,J,I)     = temperature = return_temp(planet,fpara,pressure,theta);
          EXNER3(K,J,I) = planet->cp*temperature/theta;
          mu            = avg_molar_mass(planet,kk,J,I);
          RHO3(K,J,I)   = return_density(planet,fpara,pressure,temperature,mu,PASSING_T);
        }
      }
      /* No need to apply bc_lateral() here. */
    }
  }
  else {
    for (K = 0; K <= KHI; K++) {
      kk = 2*K+1;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          theta         = THETA(K,J,I);
          pressure      = P3(K,J,I);

          /*
           * Use a root finder to solve for fp using the implicit equation:
           *   fp = fpe(T(p,theta,fp))
           * Global variables used to communicate with fpe_minus_fpe():
           */
          FPEMFPE_theta = theta;
          FPEMFPE_p     = pressure;
          find_root(0.25,1.00,fptol,&fpe,fpe_minus_fpe);
          fpara = fpe;

          T3(K,J,I)     = temperature = return_temp(planet,fpara,pressure,theta);
          EXNER3(K,J,I) = planet->cp*temperature/theta;
          mu            = avg_molar_mass(planet,kk,J,I);
          RHO3(K,J,I)   = return_density(planet,fpara,pressure,temperature,mu,PASSING_T);
        }
      }
      /* No need to apply bc_lateral() here. */
    }
  }

  /*
   * Calculate T2, EXNER2, RHO2 (and FGIBB2 if FPARA is on)
   */
  if (var.fpara.on) {
    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          theta         = THETA2(K,J,I); 
          pressure      = P2(K,J,I);

          fpara         = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,kk,J,I); 

          T2(K,J,I)     = temperature = return_temp(planet,fpara,pressure,theta);
          EXNER2(K,J,I) = planet->cp*temperature/theta;
          mu            = avg_molar_mass(planet,kk,J,I);
          RHO2(K,J,I)   = return_density(planet,fpara,pressure,temperature,mu,PASSING_T);

          return_enthalpy(planet,fpara,pressure,temperature,&fgibb,&fpe,&uoup);
          FGIBB2(K,J,I) = fgibb;
        }
      }
      /* No need to apply bc_lateral() here. */
    }
  }
  else {
    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          theta         = THETA2(K,J,I);    
          pressure      = P2(K,J,I);

          /*
           * Use a root finder to solve for fp using the implicit equation:
           *   fp = fpe(T(p,theta,fp))
           * Global variables used to communicate with fpe_minus_fpe():
           */
          FPEMFPE_theta = theta;
          FPEMFPE_p     = pressure;
          find_root(0.25,1.00,fptol,&fpe,fpe_minus_fpe);
          fpara = fpe;

          T2(K,J,I)     = temperature = return_temp(planet,fpara,pressure,theta);
          EXNER2(K,J,I) = planet->cp*temperature/theta;
          mu            = avg_molar_mass(planet,kk,J,I);
          RHO2(K,J,I)   = return_density(planet,fpara,pressure,temperature,mu,PASSING_T);
        }
      }
      /* No need to apply bc_lateral() here. */
    }
  }

  /*
   * Assume constant ratios for density.
   * Assume temperature does not change with height above top of model.
   */
  K = 0;
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      RHO2(K,J,I) = RHO3(K,J,I)*RHO3(K,J,I)/RHO2(K+1,J,I);
      T2(  K,J,I) = T3(K,J,I);
    } 
  }
  /* No need to apply bc_lateral() here. */

  /*
   * Calculate the geopotentials, PHI3 and PHI2, the Montgomery potential, MONT2,
   * and geopotential heights, Z3 and Z2.
   *
   * NOTE: There is no map factor in the hydrostatic integration. For a good
   *       reference on this subtle issue, see 
   *  Ambaum MHP, 2008, General relationships between pressure, weight and mass of a hydrostatic fluid,
   *    Proc. R. Soc. A, 464 no. 2092, 943--950, doi:10.1098/rspa.2007.0148.
   */

  g0_inv = 1./grid.g[2*KHI+1][2*(grid.nj/2)+1];

  /*
   * Integrate hydrostatic balance to get MONT2, PHI3 and PHI2.
   *
   * NOTE: It is not easy to pick an integration scheme for the hydrostatic equation
   *       that yields a numerically stable model.
   */

  /*
   * Start with bottom interface, PHI3(KHI,J,I).
   *   phi3nk_status == PASSING_PHI3NK: the values are assumed to have been calculated
   *                                    and stored before calling this function
   *
   *   phi3nk_status == CALC_PHI3NK:    the values are calculated here
   */

  if (phi3nk_status == CALC_PHI3NK) {
    K = KHI;
    if (strcmp(planet->type,"gas-giant") == 0) {
      if (grid.coord_type == COORD_ISENTROPIC) {
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            fpara       = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,2*K+1,J,I);
            PHI3(K,J,I) = MONT1D(J)-planet->cp*T3(K,J,I);
          }
        }
      }
      else {
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            PHI3(K,J,I) = PHI1D(J);
          }
        }
      }
    }
    else if (strcmp(planet->type,"terrestrial") == 0) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          PHI3(K,J,I) = PHI_SURFACE(J,I);
        }
      }
      /* No need to apply bc_lateral() here. */
    }
    else {
      sprintf(Message,"not implemented for planet->type=%s",planet->type);
      epic_error(dbmsname,Message);
    }
  }
  else if (phi3nk_status == PASSING_PHI3NK) {
    ;
  }
  else {
    sprintf(Message,"unrecognized phi3nk_status=%d",phi3nk_status);
    epic_error(dbmsname,Message);
  }

  switch(grid.coord_type) {
    case COORD_ISENTROPIC:
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          for (K = KHI; K >= KLO; K--) {
            /*
             * See Konor and Arakawa (1997), eqn. (3.3), and the accompanying discussion.
             */
            PHI3(K-1,J,I) = PHI3(K,J,I)+(EXNER2(K,J,I)-EXNER3(K-1,J,I))*THETA(K-1,J,I)
                                       +(EXNER3(K,J,I)-EXNER2(K,  J,I))*THETA(K,  J,I);

            PHI2(K,J,I) = PHI3(K-1,J,I)-(EXNER2(K,J,I)-EXNER3(K-1,J,I))*THETA(K-1,J,I);
          }
          /* 
           * Do a linear extension above the model for PHI2.
           */
          PHI2(KLO-1,J,I) = 2.*PHI3(KLO-1,J,I)-PHI2(KLO,J,I);
        }
      }
    break;

    case COORD_ISOBARIC:
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          for (K = KHI; K >= KLO; K--) {
            /*
             * Avoid a computational mode by using centered differencing.
             * See Konor and Arakawa (1997), eqn. (3.3), and the accompanying discussion.
             *
             * NOTE: Versions of hydrostatic balance that explicitly contain the density tend to not
             *       behave as well as the version dphi/dexner = -theta.
             */
            PHI3(K-1,J,I) = PHI3(K,J,I)+(EXNER3(K,J,I)-EXNER3(K-1,J,I))*THETA2(K,J,I);

            /*
             * Use the same equation as in the COORD_HYBRID case.
             */
            PHI2(K,J,I) = PHI3(K-1,J,I)-(EXNER2(K,J,I)-EXNER3(K-1,J,I))*THETA(K-1,J,I);
          }
          /* 
           * Do a linear extension above the model for PHI2.
           */
          PHI2(KLO-1,J,I) = 2.*PHI3(KLO-1,J,I)-PHI2(KLO,J,I);
        }
      }
    break;

    case COORD_HYBRID:
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          for (K = KHI; K >= KLO; K--) {
            /*
             * See Konor and Arakawa (1997), eqn. (3.3), and the accompanying discussion.
             */
            PHI3(K-1,J,I) = PHI3(K,J,I)+(EXNER2(K,J,I)-EXNER3(K-1,J,I))*THETA(K-1,J,I)
                                       +(EXNER3(K,J,I)-EXNER2(K,  J,I))*THETA(K,  J,I);

            PHI2(K,J,I) = PHI3(K-1,J,I)-(EXNER2(K,J,I)-EXNER3(K-1,J,I))*THETA(K-1,J,I);
          }
          /* 
           * Do a linear extension above the model for PHI2.
           */
          PHI2(KLO-1,J,I) = 2.*PHI3(KLO-1,J,I)-PHI2(KLO,J,I);
        }
      }
    break;

    default:
      sprintf(Message,"grid.coord_type = %d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * Calculate MONT2, and the geopotential heights, Z2, Z3, from geopotentials, PHI2, PHI3.
   */
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      for (K = KLO-1; K <= KHI; K++) {
        fpara        = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,2*K,J,I);
        MONT2(K,J,I) = PHI2(K,J,I)+planet->cp*T2(K,J,I);
      }
      for (K = KLOPAD; K <= KHIPAD; K++) {
        Z2(K,J,I) = PHI2(K,J,I)*g0_inv;
        Z3(K,J,I) = PHI3(K,J,I)*g0_inv;
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  if (action == STEP_PROGS_AND_SYNC_DIAGS) {
    /*
     * Complete the calculation of the partial derivative of Z2 with respect to time, which is
     * being stored in DZDT2.
     *
     * NOTE: The rest of the DZDT2 calculation, the advection terms, are in calc_w().
     */

    /* Negative because we need (new Z - old Z)/dt */
    dt_inv = -1./DT;
    for (K = KLO; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          DZDT2(K,J,I) -= Z2(K,J,I);
          DZDT2(K,J,I) *= dt_inv;
        }
      }
    }
    /* Apply bc_lateral() here, even though calc_w() also takes care of it. */
    bc_lateral(var.dzdt2.value,THREEDIM);
  }

  return;
}

/*==================== end of store_pgrad_vars() =============================*/

/*======================= store_diag() =======================================*/

/*
 * Calculates and stores commonly used diagnostic variables other than
 * those covered by set_p2_etc() and store_pgrad_vars().
 *
 * NOTE: This function must be called from all nodes.
 */

void store_diag(planetspec *planet) 
{
  register int
    K,J,I,
    kk,
    shift;
  EPIC_FLOAT
   *div;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="store_diag";

  /*
   * Calculate potential vorticity.
   */
  for (K = KLO; K <= KHI; K++) {
    shift = (K-Kshift)*Nelem2d;

    /* 
     * Calculate potential vorticity in the layers, PV.
     * Use the total hybrid density, H (see Schubert et al 2001, JAS 58, 3148-57).
     *
     * NOTE: For isentropic coordinates, and for the isentropic region of hybrid coordinates,
     *       this PV is Ertel's isentropic potential vorticity. For isobaric coordinates, it is
     *       proportional to absolute vorticity on pressure levels, since the thickness density, h, is constant.
     *
     * NOTE: The same type of thickness density, h, used for PV here must also be used 
     *       for UH and VH in uv_core() to get the correct Coriolis terms.
     */
    vorticity(ON_SIGMATHETA,POTENTIAL,2*K,
              var.u.value+shift+grid.it_uv*Nelem3d,
              var.v.value+shift+grid.it_uv*Nelem3d,
              var.h.value+shift+grid.it_h*Nelem3d,
              var.pv2.value+shift);
  }

  /*
   * Calculate horizontal divergence.
   */
  for (K = KLO; K <= KHI; K++) {
    shift = (K-Kshift)*Nelem2d;
    div   = var.div_uv2.value+shift;

    divergence(2*K,var.u.value+shift+grid.it_uv*Nelem3d,
                   var.v.value+shift+grid.it_uv*Nelem3d,
                   div);
  }

  /*
   * Calculate Richardson number array.
   * If not finite, set to FLOAT_MAX.
   */
  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      for (K = KLO; K <= KHI; K++) {
        kk         = 2*K;
        RI2(K,J,I) = get_richardson(planet,kk,J,I);
        if (!isfinite(RI2(K,J,I))) {
          RI2(K,J,I) = FLOAT_MAX;
        }
      }   
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.ri2.value,THREEDIM);

  /*
   * Calculate turbulence-model variables.
   */
  if (strcmp(grid.turbulence_scheme,"on")               == 0 ||
      strcmp(grid.turbulence_scheme,"on_vertical_only") == 0)  {
    set_diffusion_coef(planet);
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  return;
}

/*==================== end of store_diag() ===================================*/

/*==================== divergence() ==========================================*/

/*
 * Compute the divergence of the vector (uu,vv) using map factors. 
 * UU, VV, and DI, are assumed to be 2D arrays on the staggered C grid.
 */

void divergence(int         kk,
                EPIC_FLOAT *uu,
                EPIC_FLOAT *vv,
                EPIC_FLOAT *di)
{
  int
    J,I,
    jj;
  EPIC_FLOAT
    m_2j_inv,m_2jp2_inv,
    n_2jp1,n_2jp1_inv,
    mn_2jp1;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="divergence";

  for (J = JLO; J <= JHI; J++) {
    jj = 2*J+1;
    /*
     * Map factors needed for divergence calculation.
     * NOTE: mn != m*n at the poles, because the area is triangular. 
     */
    n_2jp1     = grid.n[kk][jj];
    n_2jp1_inv = 1./n_2jp1;
    mn_2jp1    = grid.mn[kk][jj];

    if (J == grid.jlo && IS_SPOLE) {
      m_2j_inv = 0.;
    }
    else {
      m_2j_inv = 1./grid.m[kk][jj-1];
    }

    if (J == grid.nj && IS_NPOLE) {
      m_2jp2_inv = 0.;
    }
    else {
      m_2jp2_inv = 1./grid.m[kk][jj+1];
    }

    for (I = ILO; I <= IHI; I++) {
      DI(J,I) = mn_2jp1*( (UU(J,  I+1)*n_2jp1_inv-UU(J,I)*n_2jp1_inv)
                         +(VV(J+1,I  )*m_2jp2_inv-VV(J,I)*m_2j_inv  ) );
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(di,TWODIM);

  return;
}

/*==================== end of divergence() ===================================*/

/*==================== vorticity() ===========================================*/

/*
 * Calculates vorticity for C-grid in a JI (horizontal) plane.
 * The valid types are POTENTIAL, ABSOLUTE, and RELATIVE.
 *
 * NOTE: Currently, only surface_type == ON_SIGMATHETA is implemented.
 *       We may wish to implement ON_THETA in the future.
 */

void vorticity(int         surface_type,
               int         type,
               int         kk,
               EPIC_FLOAT *uu,
               EPIC_FLOAT *vv,
               EPIC_FLOAT *hh,
               EPIC_FLOAT *pv2d)
{
  register int
    J,I;
  register EPIC_FLOAT
    f_2j,m_2j,n_2j,mn_pv,
    m_2jm1_inv,m_2jp1_inv,
    mn_2jm1_inv,mn_2jp1_inv,
    pv_pole,pvbot,pvtop;
  EPIC_FLOAT
    ze,zetabot,zetatop,
    h_pv;
#if defined(EPIC_MPI)
  EPIC_FLOAT
    mpi_tmp;
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="vorticity";

  /* 
   * Check validity of types.
   */
  if (surface_type != ON_SIGMATHETA) {
    sprintf(Message,"surface_type=%d not implemented",surface_type);
    epic_error(dbmsname,Message);
  }

  if (type == POTENTIAL) {
    if (!hh) {
      sprintf(Message,"hh=NULL");
      epic_error(dbmsname,Message);
    }
  }
  else if (type != RELATIVE && type != ABSOLUTE) {
    sprintf(Message,"unrecognized type=%d",type);
    epic_error(dbmsname,Message);
  }

  /* 
   * Calculate interior points.
   */ 
  for (J = JFIRST; J <= JHI; J++) {
    if (type == RELATIVE) {
      f_2j = 0.;
    }
    else {
      f_2j = (grid.f)[2*J];
    }
    m_2j        =    grid.m[ kk][2*J  ];
    n_2j        =    grid.n[ kk][2*J  ];
    m_2jm1_inv  = 1./grid.m[ kk][2*J-1];
    m_2jp1_inv  = 1./grid.m[ kk][2*J+1];
    mn_2jm1_inv = 1./grid.mn[kk][2*J-1];
    mn_2jp1_inv = 1./grid.mn[kk][2*J+1];
    mn_pv       = .5/(mn_2jm1_inv+mn_2jp1_inv);
    for (I = ILO; I <= IHI; I++) {
      ze        =  m_2j*( (VV(J,I)-VV(J,I-1)) +
                   n_2j*(UU(J-1,I)*m_2jm1_inv-UU(J,I)*m_2jp1_inv) );
      PV2D(J,I) = (ze+f_2j);
    }
    if (type == POTENTIAL) {
      for (I = ILO; I <= IHI; I++) {
        h_pv = ( (HH(J,  I)+HH(J,  I-1))*mn_2jp1_inv  
                +(HH(J-1,I)+HH(J-1,I-1))*mn_2jm1_inv )*mn_pv;
        PV2D(J,I) /= h_pv;
      }
    }
  }

  if (strcmp(grid.geometry,"f-plane")  == 0 &&
      strcmp(grid.f_plane_map,"polar") == 0) {
    if (JLO == 0 && grid.nj > grid.jlo) {
      /*  Apply channel boundary condition for pv: */
      J           = 1;
      m_2j        =    grid.m[kk][2*J  ];
      n_2j        =    grid.n[kk][2*J  ];
      m_2jm1_inv  = 1./grid.m[kk][2*J-1];
      m_2jp1_inv  = 1./grid.m[kk][2*J+1];
      /* 
       * Calculate average zeta in next row.
       */
      zetabot = 0.;
      for (I = ILO; I <= IHI; I++) {
        zetabot += m_2j*( (VV(J,I)-VV(J,I-1)) +
                   n_2j*(UU(J-1,I)*m_2jm1_inv-UU(J,I)*m_2jp1_inv) );
      }

#if defined(EPIC_MPI)
      mpi_tmp = zetabot;
      MPI_Allreduce(&mpi_tmp,&zetabot,1,float_type,MPI_SUM,para.comm_JLO);
#endif

      zetabot /= grid.ni;

      if (type == RELATIVE) {
        pvbot = zetabot;
      }
      else if (type == ABSOLUTE) {
        pvbot = zetabot+grid.f[2*J];
      }
      else if (type == POTENTIAL) {
        h_pv = 0.;
        J   = 0;
        for (I = ILO; I <= IHI; I++) {
          h_pv  += HH(J,I);
        }

#if defined(EPIC_MPI)
        mpi_tmp = h_pv;
        MPI_Allreduce(&mpi_tmp,&h_pv,1,float_type,MPI_SUM,para.comm_JLO);
#endif

        h_pv  /= grid.ni;

        pvbot  = (zetabot+grid.f[2*J])/h_pv;
      }
      else {
        sprintf(Message,"unrecognized type=%d",type);
        epic_error(dbmsname,Message);
      }

      J = 0;
      for (I = ILO; I <= IHI; I++) {
        PV2D(J,I) = pvbot;
      }
    }
    if (IS_NPOLE) {
      /* Calculate "north pole" pv */
      J   = grid.nj;
      ze  = 0.;
      for (I = ILO; I <= IHI; I++) {
        ze += UU(J,I);
      }

#if defined(EPIC_MPI)
      mpi_tmp = ze;
      MPI_Allreduce(&mpi_tmp,&ze,1,float_type,MPI_SUM,para.comm_JLO);
#endif

      ze *= grid.mn[kk][2*(J+1)]/(grid.m[kk][2*J+1]*(EPIC_FLOAT)grid.ni);

      if (type == RELATIVE) {
        pv_pole = ze;
      }
      else if (type == ABSOLUTE) {
        pv_pole = ze+grid.f[2*(J+1)];
      }
      else if (type == POTENTIAL) {
        h_pv = 0.;
        for (I = ILO; I <= IHI; I++) {
          h_pv += HH(J,I);
        }

#if defined(EPIC_MPI)
        mpi_tmp = h_pv;
        MPI_Allreduce(&mpi_tmp,&h_pv,1,float_type,MPI_SUM,para.comm_JLO);
#endif

        h_pv    /= grid.ni;

        pv_pole  = (ze+grid.f[2*(J+1)])/h_pv;
      }
      else {
        sprintf(Message,"unrecognized type=%d",type);
        epic_error(dbmsname,Message);
      }

      J = grid.nj+1;
      for (I = ILO; I <= IHI; I++) {
        PV2D(J,I) = pv_pole;
      }
    }
  }
  else if (strcmp(grid.geometry,"globe") == 0) {
    if (grid.globe_latbot == -90.) {
      if (IS_SPOLE) {
        /* Calculate pv at the south pole: */
        J  = 0;
        ze = 0.;
        for (I = ILO; I <= IHI; I++) {
          ze -= UU(J,I);  /*  Beware of southern circulation sign. */
        }

#if defined(EPIC_MPI)
        mpi_tmp = ze;
        MPI_Allreduce(&mpi_tmp,&ze,1,float_type,MPI_SUM,para.comm_JLO);
#endif

        ze *= grid.mn[kk][0]/((grid.m[kk][1])*(EPIC_FLOAT)grid.ni); 

        if (type == RELATIVE) {
          pv_pole = ze;
        }
        else if (type == ABSOLUTE) {
          pv_pole = ze+grid.f[0];
        }
        else if (type == POTENTIAL) {
          h_pv = 0.;
          for (I = ILO; I <= IHI; I++) {
            h_pv += HH(J,I);
          }

#if defined(EPIC_MPI)
          mpi_tmp = h_pv;
          MPI_Allreduce(&mpi_tmp,&h_pv,1,float_type,MPI_SUM,para.comm_JLO);
#endif

          h_pv    /= grid.ni;

          pv_pole  = (ze+grid.f[0])/h_pv;
        }
        else {
          sprintf(Message,"unrecognized type=%d",type);
          epic_error(dbmsname,Message);
        }
        for (I = ILO; I <= IHI; I++) {
          PV2D(J,I) = pv_pole;
        }
      }
    }
    else {
      if (JLO == 0) {
        /*  Apply channel boundary condition for pv: */
        if (grid.globe_latbot == 0.) {
          /* special case at equator */
          pvbot = 0.;
        }
        else if (grid.nj == grid.jlo) {
          pvbot = 0.;
        }
        else {
          /* Calculate average zeta in next row */
          J           = 1;
          m_2j        =    grid.m[kk][2*J  ];
          n_2j        =    grid.n[kk][2*J  ];
          m_2jm1_inv  = 1./grid.m[kk][2*J-1];
          m_2jp1_inv  = 1./grid.m[kk][2*J+1];

          zetabot = 0.;
          for (I = ILO; I <= IHI; I++) {
            zetabot +=  m_2j*( (VV(J,I)-VV(J,I-1)) +
                        n_2j*(UU(J-1,I)*m_2jm1_inv-UU(J,I)*m_2jp1_inv) );
          }

#if defined(EPIC_MPI)
          mpi_tmp = zetabot;
          MPI_Allreduce(&mpi_tmp,&zetabot,1,float_type,MPI_SUM,para.comm_JLO);
#endif

          zetabot /= grid.ni;

          if (type == RELATIVE) {
            pvbot = zetabot;
          }
          else if (type == ABSOLUTE) {
            J     = 0;
            pvbot = zetabot+grid.f[2*J];
          }
          else if (type == POTENTIAL) {
            h_pv = 0.;
            J    = 0;
            for (I = ILO; I <= IHI; I++) {
              h_pv += HH(J,I);
            }

#if defined(EPIC_MPI)
            mpi_tmp = h_pv;
            MPI_Allreduce(&mpi_tmp,&h_pv,1,float_type,MPI_SUM,para.comm_JLO);
#endif

            h_pv  /= grid.ni;

            J      = 0;
            pvbot  = (zetabot+grid.f[2*J])/h_pv;
          }
          else {
            sprintf(Message,"unrecognized type=%d",type);
            epic_error(dbmsname,Message);
          }
        }

        J = 0;
        for (I = ILO; I <= IHI; I++) {
          PV2D(J,I) = pvbot;
        }
      }
    }
    if (grid.globe_lattop == 90.) {
      if (IS_NPOLE) {
        /* Calculate pv at the north pole: */
        J   = grid.nj;
        ze  = 0.;
        for (I = ILO; I <= IHI; I++) {
          ze += UU(J,I);
        }

#if defined(EPIC_MPI)
        mpi_tmp = ze;
        MPI_Allreduce(&mpi_tmp,&ze,1,float_type,MPI_SUM,para.comm_JLO);
#endif

        ze *= grid.mn[kk][2*(J+1)]/(grid.m[kk][2*J+1]*(EPIC_FLOAT)grid.ni);

        if (type == RELATIVE) {
          pv_pole = ze;
        }
        else if (type == ABSOLUTE) {
          pv_pole = ze+grid.f[2*(J+1)];
        }
        else if (type == POTENTIAL) {
          h_pv = 0.;
          for (I = ILO; I <= IHI; I++) {
            h_pv += HH(J,I);
          }

#if defined(EPIC_MPI)
          mpi_tmp = h_pv;
          MPI_Allreduce(&mpi_tmp,&h_pv,1,float_type,MPI_SUM,para.comm_JLO);
#endif

          h_pv    /= grid.ni;

          pv_pole  = (ze+grid.f[2*(J+1)])/h_pv;
        }
        else {
          sprintf(Message,"unrecognized type=%d",type);
          epic_error(dbmsname,Message);
        }

        J = grid.nj+1;
        for (I = ILO; I <= IHI; I++) {
          PV2D(J,I) = pv_pole;
        }
      }
    }
    else {
      if (JHI == grid.nj) {
        /*  Apply channel boundary condition for pv: */
        if (grid.globe_lattop == 0.) {
          /* special case at equator */
          pvtop = 0.;
        }
        else if (grid.nj == grid.jlo) {
          pvtop = 0.;
        }
        else {
          /* Calculate average zeta in next row */
          J           = grid.nj;
          m_2j        =    grid.m[kk][2*J  ];
          n_2j        =    grid.n[kk][2*J  ];
          m_2jm1_inv  = 1./grid.m[kk][2*J-1];
          m_2jp1_inv  = 1./grid.m[kk][2*J+1];

          zetatop = 0.;
          for (I = ILO; I <= IHI; I++) {
            zetatop +=  m_2j*( (VV(J,I)-VV(J,I-1)) +
                        n_2j*(UU(J-1,I)*m_2jm1_inv-UU(J,I)*m_2jp1_inv) );
          }

#if defined(EPIC_MPI)
          mpi_tmp = zetatop;
          MPI_Allreduce(&mpi_tmp,&zetatop,1,float_type,MPI_SUM,para.comm_JLO);
#endif

          zetatop /= grid.ni;

          if (type == RELATIVE) {
            pvtop = zetatop;
          }
          else if (type == ABSOLUTE) {
            J     = grid.nj+1;
            pvtop = zetatop+grid.f[2*J];
          }
          else if (type == POTENTIAL) {
            J    = grid.nj;
            h_pv = 0.;
            for (I = ILO; I <= IHI; I++) {
              h_pv += HH(J,I);
            }

#if defined(EPIC_MPI)
            mpi_tmp = h_pv;
            MPI_Allreduce(&mpi_tmp,&h_pv,1,float_type,MPI_SUM,para.comm_JLO);
#endif

            h_pv  /= grid.ni;

            J      = grid.nj+1;
            pvtop  = (zetatop+grid.f[2*J])/h_pv;
          }
          else {
            sprintf(Message,"unrecognized type=%d",type);
            epic_error(dbmsname,Message);
          }
        }

        J = grid.nj+1;
        for (I = ILO; I <= IHI; I++) {
          PV2D(J,I) = pvtop;
        }
      }
    }
  }

  /* Need to apply bc_lateral() here. */
  bc_lateral(pv2d,TWODIM);

  return;
}

/*======================= end of vorticity() ===================================*/

/*======================= phi_from_u() =========================================*/

/*
 * Calculate geopotential, phi, assuming gradient balance with a zonally 
 * symmetric zonal wind, u, defined on an isobaric surface.
 * The input parameter phi0 is the value at J = j0.
 *
 * NOTE: K is used internally here to make the PV macros work, so do not assign
 *       to K the value corresponding to the input kk.
 */
#include <epic_pv_schemes.h>

void phi_from_u(planetspec *planet,
                int         kk,
                EPIC_FLOAT *u1d,
                EPIC_FLOAT *phi1d,
                int         j0,
                EPIC_FLOAT  phi0)
{
  register int
    K,J,I;
  EPIC_FLOAT
    kin,kin0;
  register EPIC_FLOAT
    uu;
  static int
    initialized=FALSE;
  static EPIC_FLOAT
    *u2d,
    *v2d,
    *phi2d,
    *udy,
    *pvhudy,
    *bern,
    *sendbuf;
#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="phi_from_u";

  if (!initialized) {
    /* 
     * Verify vertical coordinate is not isentropic.
     */
    if (grid.coord_type == COORD_ISENTROPIC) {
      sprintf(Message,"Calling %s with grid.vertical_coordinate=%s\n",dbmsname,grid.vertical_coordinate);
      epic_error(dbmsname,Message);
    }

    /* Allocate memory. */
    udy     = fvector(0,grid.nj+1,dbmsname);
    pvhudy  = fvector(0,grid.nj+1,dbmsname);
    bern    = fvector(0,grid.nj+1,dbmsname);
    sendbuf = fvector(0,grid.nj+1,dbmsname);
    u2d     = fvector(0,Nelem2d-1,dbmsname);
    v2d     = fvector(0,Nelem2d-1,dbmsname);
    phi2d   = fvector(0,Nelem2d-1,dbmsname);

    initialized = TRUE;
  }

  /*
   * This procedure assumes zonal symmetry.
   */
  for (J = JLOPAD; J <= JHIPAD; J++) {
    uu = U1D(J);
    for (I = ILOPAD; I <= IHIPAD; I++) {
      U2D(J,I) = uu;
    }
  }

  /*
   * Calculate absolute vorticity.
   * Store in PV(nk+1,J,I), which is otherwise not used, 
   * so that epic_pv_schemes.h macros (GA_V, DE_V, etc) work.
   *
   * NOTE: kk and K are not related here.
   */
  K = grid.nk+1;
  vorticity(ON_SIGMATHETA,ABSOLUTE,kk,u2d,v2d,NULL,var.pv2.value+(K-Kshift)*Nelem2d);

  kin0 = 0.;
  if (ILO == grid.ilo) {
    if (j0 >= JLO && j0 <= JHI) {
      kin0 = get_kin(planet,u2d,v2d,kk,j0,grid.ilo);
    }

    /*
     * Calculate udy.
     */
    memset(udy,0,(grid.nj+2)*sizeof(EPIC_FLOAT));
    for (J = JLOPAD; J <= JHI; J++) {
      udy[J] = U1D(J)/grid.n[kk][2*J+1];
    }

    /*
     * Calculate (zeta+f)*u*dy = pvhudy.
     */
    memset(pvhudy,0,(grid.nj+2)*sizeof(EPIC_FLOAT));
    I = grid.ilo;
    for (J = JFIRST; J <= JHI; J++) {
      /* 
       * Don't multiply by grid.n to leave in dy factor:
       *
       * NOTE: the GA_V, etc. coefficients involve a stencil on PV2(K,J,I),
       *       such that both K and I must be set.
       */
      pvhudy[J] = (GA_V*udy[J  ]+DE_V*udy[J  ]
                  +AL_V*udy[J-1]+BE_V*udy[J-1])*PV_COEF;
    }
  }

#if defined(EPIC_MPI)
  /* 
   * Broadcast kin0, which is the value at J=j0, I=grid.ilo.
   */
  sendbuf[0] = kin0;
  MPI_Allreduce(sendbuf,&kin0,1,float_type,MPI_SUM,para.comm);

  /*
   *  Fill in pvhudy for each node:
   */ 
  memcpy(sendbuf,pvhudy,(grid.nj+2)*sizeof(EPIC_FLOAT));
  MPI_Allreduce(sendbuf,pvhudy,grid.nj+2,float_type,MPI_SUM,para.comm);
#endif

  /*
   * Calculate bern by integrating -pvhudy.
   *
   * For convenience, this global-spanning integral is calculated 
   * redundantly on every node.
   */
  bern[j0] = phi0+kin0;
  for (J = j0+1; J <= grid.nj; J++) {
    bern[J] = bern[J-1]-pvhudy[J];
  }
  for (J = j0-1; J >= grid.jlo; J--) {
    bern[J] = bern[J+1]+pvhudy[J+1];
  }

  /*
   * Calculate phi = bern-kin.
   * Store in PHI2D, so that bc_lateral() may be applied.
   */
  for (J = JLO; J <= JHI; J++) {
    kin          = get_kin(planet,u2d,v2d,kk,J,ILO);
    PHI2D(J,ILO) = bern[J]-kin;
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(phi2d,TWODIM);

  for (J = JLOPAD; J <= JHIPAD; J++) {
    PHI1D(J) = PHI2D(J,ILO);
  }

  return;
}

/*======================= end of phi_from_u() ==================================*/

/*======================= mont_from_u() ========================================*/

/*
 * Calculate Montgomery potential, mont, assuming gradient balance with a zonally 
 * symmetric zonal wind, u, defined on an isentropic surface.
 * The input parameter mont0 is the value at J = j0.
 *
 * NOTE: K is used internally here to make the PV macros work, so do not assign
 *       to K the value corresponding to the input kk.
 */
#include <epic_pv_schemes.h>

void mont_from_u(planetspec *planet,
                 int         kk,
                 EPIC_FLOAT *u1d,
                 EPIC_FLOAT *mont1d,
                 int         j0,
                 EPIC_FLOAT  mont0)
{
  register int
    K,J,I;
  EPIC_FLOAT
    kin,kin0;
  register EPIC_FLOAT
    uu;
  static int
    initialized=FALSE;
  static EPIC_FLOAT
    *u2d,
    *v2d,
    *mont2d,
    *udy,
    *pvhudy,
    *bern,
    *sendbuf;
#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="mont_from_u";

  if (!initialized) {
    /* 
     * Verify vertical coordinate is isentropic.
     */
    if (grid.coord_type != COORD_ISENTROPIC) {
      sprintf(Message,"Calling %s with grid.vertical_coordinate=%s\n",dbmsname,grid.vertical_coordinate);
      epic_error(dbmsname,Message);
    }

    /* Allocate memory. */
    udy      = fvector(0,grid.nj+1,dbmsname);
    pvhudy   = fvector(0,grid.nj+1,dbmsname);
    bern     = fvector(0,grid.nj+1,dbmsname);
    sendbuf  = fvector(0,grid.nj+1,dbmsname);
    u2d      = fvector(0,Nelem2d-1,dbmsname);
    v2d      = fvector(0,Nelem2d-1,dbmsname);
    mont2d   = fvector(0,Nelem2d-1,dbmsname);

    initialized = TRUE;
  }

  /*
   * This procedure assumes zonal symmetry.
   */
  for (J = JLOPAD; J <= JHIPAD; J++) {
    uu = U1D(J);
    for (I = ILOPAD; I <= IHIPAD; I++) {
      U2D(J,I) = uu;
    }
  }

  /*
   * Calculate absolute vorticity.
   * Store in PV(nk+1,J,I), which is otherwise not used, 
   * so that epic_pv_schemes.h macros (GA_V, DE_V, etc) work.
   *
   * NOTE: kk and K are not related here.
   */
  K = grid.nk+1;
  vorticity(ON_SIGMATHETA,ABSOLUTE,kk,u2d,v2d,NULL,var.pv2.value+(K-Kshift)*Nelem2d);

  kin0 = 0.;
  if (ILO == grid.ilo) {
    if (j0 >= JLO && j0 <= JHI) {
      kin0 = get_kin(planet,u2d,v2d,kk,j0,grid.ilo);
    }

    /*
     * Calculate udy.
     */
    memset(udy,0,(grid.nj+2)*sizeof(EPIC_FLOAT));
    for (J = JLOPAD; J <= JHI; J++) {
      udy[J] = U1D(J)/grid.n[kk][2*J+1];
    }

    /*
     * Calculate (zeta+f)*u*dy = pvhudy.
     */
    memset(pvhudy,0,(grid.nj+2)*sizeof(EPIC_FLOAT));
    I = grid.ilo;
    for (J = JFIRST; J <= JHI; J++) {
      /* 
       * Don't multiply by grid.n to leave in dy factor:
       *
       * NOTE: the GA_V, etc. macro coefficients involve a stencil on PV2(K,J,I),
       *       such that both K and I must be set.
       */
      pvhudy[J] = (GA_V*udy[J  ]+DE_V*udy[J  ]
                  +AL_V*udy[J-1]+BE_V*udy[J-1])*PV_COEF;
    }
  }

#if defined(EPIC_MPI)
  /* 
   * Broadcast kin0, which is the value at J=j0, I=grid.ilo.
   */
  sendbuf[0] = kin0;
  MPI_Allreduce(sendbuf,&kin0,1,float_type,MPI_SUM,para.comm);

  /*
   *  Fill in pvhudy for each node:
   */ 
  memcpy(sendbuf,pvhudy,(grid.nj+2)*sizeof(EPIC_FLOAT));
  MPI_Allreduce(sendbuf,pvhudy,grid.nj+2,float_type,MPI_SUM,para.comm);
#endif

  /*
   * Calculate bern by integrating -pvhudy.
   *
   * For convenience, this global-spanning integral is calculated 
   * redundantly on every node.
   */
  bern[j0] = mont0+kin0;
  for (J = j0+1; J <= grid.nj; J++) {
    bern[J] = bern[J-1]-pvhudy[J];
  }
  for (J = j0-1; J >= grid.jlo; J--) {
    bern[J] = bern[J+1]+pvhudy[J+1];
  }

  /*
   * Calculate mont = bern-kin.
   * Store in MONT2D, so that bc_lateral() may be applied.
   */
  for (J = JLO; J <= JHI; J++) {
    kin           = get_kin(planet,u2d,v2d,kk,J,ILO);
    MONT2D(J,ILO) = bern[J]-kin;
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(mont2d,TWODIM);

  for (J = JLOPAD; J <= JHIPAD; J++) {
    MONT1D(J) = MONT2D(J,ILO);
  }

  return;
}

/*======================= end of mont_from_u() =================================*/

/*======================= relative_humidity() ==================================*/

/*
 * Calculate the relative humidity (RH) at the bottom of layer K.
 *
 * See Bohren and Albrecht (1998, Atmospheric Thermodynamics, Oxford, p.186)
 * for an interesting discussion about different definitions of relative humidity.
 *
 * We use the WMO definition because microphysical equations are usually
 * written in terms of mixing ratio rather than partial pressure. 
 */
void relative_humidity(planetspec *planet,
                       int         is,
                       EPIC_FLOAT *rh,
                       int         K)
{
  register int
    J,I;
  EPIC_FLOAT
    pp,psat;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="relative_humidity";

#if EPIC_CHECK == 1
  /*
   * Check that 'is' is valid.
   */
  if (is < FIRST_SPECIES) {
    sprintf(Message,"is=%d < FIRST_SPECIES=%d",is,FIRST_SPECIES);
    epic_error(dbmsname,Message);
  }
  if (is > LAST_SPECIES) {
    sprintf(Message,"is=%d > LAST_SPECIES=%d",is,LAST_SPECIES);
    epic_error(dbmsname,Message);
  }
#endif

  /*
   * Return 0. if the species is not on.
   */
  if (!var.species[is].on) {
    memset(rh,0,Nelem2d*sizeof(EPIC_FLOAT));
    return;
  }

#if EPIC_CHECK == 1
  /*
   * Check that the phase VAPOR is on for the species.
   */
  if (!var.species[is].phase[VAPOR].on) {
    sprintf(Message,"var.species[%d].phase[VAPOR].on is not on",is);
    epic_error(dbmsname,Message);
  }
#endif

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /*
       * Start with the "old school" definition of RH as the ratio of
       * partial pressure to saturation pressure.
       */
      pp      = get_p(planet,is,2*K+1,J,I);
      psat    = var.species[is].sat_vapor_p(T3(K,J,I));
      RH(J,I) = pp/psat;

     /*
      * Now factor in the term that converts this to the ratio of mixing ratio to
      * saturation mixing ratio.
      */
      RH(J,I) *= 1.-(psat-pp)/PDRY3(K,J,I);
    }
  }
  /* No need to apply bc_lateral() here. */

  return;
}

/*======================= end of relative_humidity() ============================*/

/*======================= source_sink() =========================================*/

/*
 * Apply source and sink terms to h-grid variables.
 * Use a forward timestep to integrate rates.
 */

void source_sink(planetspec  *planet,
                 EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I;
  register EPIC_FLOAT
    dfpdt,
    fpara,pressure,temperature,
    time_fp_inv;
  EPIC_FLOAT
    fgibb,fpe,uoup;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="source_sink";

  if (var.fpara.on) {
    /*
     * Advance FPARA, 
     */
    for (K = KLO; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          fpara       = FPARA(K,J,I);
          pressure    = P3(   K,J,I);
          temperature = T3(   K,J,I);

          return_enthalpy(planet,fpara,pressure,temperature,&fgibb,&fpe,&uoup);

          /*
           * The ortho-para conversion rate [1/s] is from Huestis (2008, Planet. Space Sci. 56, 1733--1743), 
           * Fig. 7. This is an updated fit to data, and at 1 bar is about 7 times faster than the rate
           * used in Dowling et al (1998).
           *
           * The input parameter var.fpara_rate_scaling is nominally 1.0.
           */
          time_fp_inv  = var.fpara_rate_scaling
                        *(1./5.2e7)*(pressure/1.e+5)*(planet->x_h2/0.85)*(50./temperature)
                        *(1.+7.82*exp(-173./temperature));
          dfpdt        = (fpe-fpara)*time_fp_inv;
          FPARA(K,J,I) += DT*dfpdt;
        }
      }
    }
    bc_lateral(var.fpara.value,THREEDIM);
  }

  /* 
   * Add sources and sinks from microphysical processes.
   *
   * The function cloud_microphysics() assumes that the vapor, liquid, and solid
   * phases of each active species is turned on in epic_initial.c.
   *
   */
  if (grid.cloud_microphysics == ACTIVE) {
    cloud_microphysics(planet,Buff2D);
  }

  /*
   * Advance THETA where it is a prognostic variable.
   */
  switch(grid.coord_type) {
    case COORD_ISENTROPIC:
      ;
    break;
    case COORD_ISOBARIC:
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          /*
           * Apply heating term.
           */
          for (K = KLO; K <= KHI; K++) {
            THETA(K,J,I) += DT*HEAT3(K,J,I)/EXNER3(K,J,I);
          }
          /*
           * Extrapolate THETA to the top of the model.
           */
          THETA(KLO-1,J,I) = THETA(KLO,J,I)+(THETA(KLO,J,I)-THETA(KLO+1,J,I))
                                           *grid.dsgth[2*KLO]*grid.dsgth_inv[2*(KLO+1)];
        }
        /* No need to apply bc_lateral() here. */
      }
    break;
    case COORD_HYBRID:
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          /*
           * Apply heating term.
           */
          for (K = KLO; K <= KHI; K++) {
            THETA(K,J,I) += DT*HEAT3(K,J,I)/EXNER3(K,J,I);
          }
        }
        /* No need to apply bc_lateral() here. */
      }
    break;
  }

  return;
}

/*======================= end of source_sink() ==================================*/

/*======================= calc_heating() ========================================*/

/*
 * Calculate HEAT3, in W/kg.
 *
 * NOTE: Latent heating from clouds is handled outside this function.
 */

void calc_heating(planetspec *planet)
{
  register int
    K,J,I;
  register EPIC_FLOAT
    cprh2,
    fpara,pressure,temperature,
    time_fp_inv,
    dfpdt;
  EPIC_FLOAT
    fgibb,fpe,uoup,
    theta_o,theta_p;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="calc_heating";

  /*
   * Zero HEAT array.
   */
  memset(var.heat3.value,0,Nelem3d*sizeof(EPIC_FLOAT));

  /*
   * Add radiative heating/cooling to HEAT3.
   * The scheme used is specified by grid.radiation_scheme.
   */
  radiative_heating(planet,EPIC_APPLY);

  /*
   * Add ortho-para hydrogen heating to HEAT3.
   */
  if (var.fpara.on) {
    cprh2 = CPRH2*planet->rgas;
    for (K = KLO-1; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          fpara       = FPARA(K,J,I);
          pressure    = P3(    K,J,I);
          temperature = T3(    K,J,I);

          return_enthalpy(planet,fpara,pressure,temperature,&fgibb,&fpe,&uoup);

          /*
           * The ortho-para conversion rate [1/s] is from Huestis (2008, Planet. Space Sci. 56, 1733--1743), 
           * Fig. 7. This is an updated fit to data, and at 1 bar is about 7 times faster than the rate
           * used in Dowling et al (1998).
           *
           * The input parameter var.fpara_rate_scaling is nominally 1.0.
           */
          time_fp_inv  = var.fpara_rate_scaling
                        *(1./5.2e7)*(pressure/1.e+5)*(planet->x_h2/0.85)*(50./temperature)
                        *(1.+7.82*exp(-173./temperature));
          dfpdt        = (fpe-fpara)*time_fp_inv;
          /* 
           * Call return_theta() to get theta_o,theta_p. 
           * These only depend on p,T.
           */
          return_theta(planet,fpara,pressure,temperature,&theta_o,&theta_p);

          HEAT3(K,J,I) += (planet->x_h2)*(uoup+cprh2*temperature)*log(theta_p/theta_o)*dfpdt;  
        }
      }
    }
  }

  /*
   * Apply lateral boundary conditions to HEAT3 array.
   */
  bc_lateral(var.heat3.value,THREEDIM);

  return;
}

/*======================= end of calc_heating() =================================*/

/*======================= cfl_dt() ==============================================*/

/* 
 * Estimate CFL timestep for numerical stability.
 * Use sound speed as an upper bound on gravity-wave speed.
 *
 * NOTE: For spin-up experiments, this estimate of dt_cfl will initially be
 *       based on zero winds, and may be too large to keep developing jets
 *       numerically stable.
 */

int cfl_dt(planetspec *planet)
{
  register int
    K,J,I,
    min_cfl_dt;
  register EPIC_FLOAT
    u,v,w,cs,
    dx,dy,dz;
  EPIC_FLOAT
    cflx,cfly,cflz,
    cfl_dt,
    tmp;

#if defined(EPIC_MPI)
  int
    itmp;
# if EPIC_PRECISION == DOUBLE_PRECISION
    MPI_Datatype
      float_type = MPI_DOUBLE;
# else
    MPI_Datatype
      float_type = MPI_FLOAT;
# endif
#endif

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="cfl_dt";

  /* 
   * Analyze each direction in turn, since the
   * advection schemes operate this way.
   */
  
  /*
   * Vertical direction.
   */
  cflz = FLOAT_MAX;
  for (K = KLO; K <= KHI; K++) {
    dz = grid.dsgth[2*K];
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        /*
         * Neglect vertical component of gravity-wave speed.
         */
        w = W3(K-1,J,I);
        if (w > 0.) {
          cflz = MIN(cflz,dz/w);
        }

        w = W3(K,J,I);
        if (w < 0.) {
          cflz = MIN(cflz,-dz/w);
        }
      }
    }
  }

#if defined(EPIC_MPI)
  /* determine global minimum */
  tmp = cflz;
  MPI_Allreduce(&tmp,&cflz,1,float_type,MPI_MIN,para.comm);
#endif

  /*
   * Meridional direction.
   */
  cfly = FLOAT_MAX;
  for (J = JFIRST; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      for (K = KLO; K <= KHI; K++) {
        dy = 1./grid.n[2*K][2*J+1];
        /* 
         * Use speed of sound to estimate fastest horizontal gravity-wave speed.
         */
        cs = sqrt(planet->rgas*.5*(T2(K,J,I)+T2(K,J-1,I))/(1.-planet->kappa));

        v = V(grid.it_uv,K,J+1,I)+cs;
        if (v > 0.) {
          cfly = MIN(cfly,dy/v);
        }
        v = V(grid.it_uv,K,J,I)-cs;
        if (v < 0.) {
          cfly = MIN(cfly,-dy/v);
        }
      }
    }
  }

#if defined(EPIC_MPI)
  /* determine global minimum */
  tmp = cfly;
  MPI_Allreduce(&tmp,&cfly,1,float_type,MPI_MIN,para.comm);
#endif

  /*
   * Zonal direction.
   */
  cflx = FLOAT_MAX;
  for (J = JLO; J <= JHI; J++) {
    if (fabs(grid.lat[2*J+1]) <= LAT0) {
      for (I = ILO; I <= IHI; I++) {
        for (K = KLO; K <= KHI; K++) {
          dx = 1./grid.m[2*K][2*J+1];
          /* 
           * Use speed of sound to estimate fastest gravity-wave speed.
           */
          cs = sqrt(planet->rgas*.5*(T2(K,J,I)+T2(K,J,I-1))/(1.-planet->kappa));

          u = U(grid.it_uv,K,J,I+1)+cs;
          if (u > 0.) {
            cflx = MIN(cflx,dx/u);
          }
          u = U(grid.it_uv,K,J,I)-cs;
          if (u < 0.) {
            cflx = MIN(cflx,-dx/u);
          }
        }
      }
    }
  }

#if defined(EPIC_MPI)
  /* determine global minimum */
  tmp = cflx;
  MPI_Allreduce(&tmp,&cflx,1,float_type,MPI_MIN,para.comm);
#endif

  cfl_dt = MIN(cflx,cfly);
  cfl_dt = MIN(cfl_dt,cflz);

  min_cfl_dt = (int)cfl_dt;

  return min_cfl_dt;
}

/*======================= end of cfl_dt() ======================================*/

/*======================= time_mod() ===========================================*/

/* 
 * Return integer remainder of time/step.
 * Avoids directly forming YEAR*time.years to stay below INT_MAX.
 */

int time_mod(int *time,
             int  step) {
  int
    ans;

  ans = (time[0]%step+(YEAR%step)*(time[1]%step))%step;

  return ans;
}

/*====================== end of time_mod() ===================================*/

/*====================== avg_molar_mass() ====================================*/

/*
 * Computes average molar mass at position kk/2,j,i.
 */

EPIC_FLOAT avg_molar_mass(planetspec *planet,
                          int         kk,
                          int         J,
                          int         I) 
{
  register int
    K,
    iq;
  const EPIC_FLOAT
    mu_dry = R_GAS/planet->rgas;
  register EPIC_FLOAT
    ans;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="avg_molar_mass";

  if (kk%2 == 0) {
    /* 
     * Layer value.
     */
    ans = .5*(avg_molar_mass(planet,kk-1,J,I)+
              avg_molar_mass(planet,kk+1,J,I));
    return ans;
  }
  else {
    /*
     * Interface value.
     */
    K = (kk-1)/2;
    ans = mu_dry;
    for (iq = 0; iq < grid.nq; iq++) {
      ans += X(grid.is[iq],grid.ip[iq],K,J,I)*(var.species[grid.is[iq]].molar_mass-mu_dry);
    }
    return ans;
  }
}

/*====================== end of avg_molar_mass() =============================*/

/*====================== sync_x_to_q() =======================================*/

/*
 * Synchronize all number fractions, X_i, to mass mixing ratios, Q_i.
 *
 * NOTE: Assumes the bookkeeping arrays grid.is[] and grid.ip[] have been set up.
 */

void sync_x_to_q(planetspec *planet)
{
  int
    K,J,I,iq;
  const EPIC_FLOAT
    mu_dry_inv = planet->rgas/R_GAS;
  register EPIC_FLOAT
    sum;

  for (K = KLOPAD; K <= KHIPAD; K++) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        sum = mu_dry_inv;
        for (iq = 0; iq < grid.nq; iq++) {
          sum += Q(grid.is[iq],grid.ip[iq],K,J,I)/var.species[grid.is[iq]].molar_mass;
        }
        for (iq = 0; iq < grid.nq; iq++) {
          X(grid.is[iq],grid.ip[iq],K,J,I) = Q(grid.is[iq],grid.ip[iq],K,J,I)/(var.species[grid.is[iq]].molar_mass*sum);
        }
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  return;
}

/*====================== end of sync_x_to_q() ================================*/

/*====================== molar_mass() ========================================*/

/*
 * Returns the molar mass (molecular weight) for the indicated substance.
 * Units are kg/kmol, which is the same as g/mol.
 */
EPIC_FLOAT molar_mass(int index)
{
  register EPIC_FLOAT 
    mu;
  register int
    i,ii;
  static int
    num_elements=0,
    *counts     =NULL;
  static char
    **symbols   =NULL;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="molar_mass";

  if (index == HDRY3_INDEX || index == HDRY2_INDEX) {
    mu = R_GAS/planet->rgas;
  }
  else {
    parse_species_name(var.species[index].info[0].name,&num_elements,&symbols,&counts);
    mu = 0.;
    for (i = 0; i < num_elements; i++) {
      ii = 1;
      /* Identify element */
      while(strcmp(Element[ii].symbol,symbols[i]) != 0) {
        ii++;
      };
      mu += counts[i]*(Element[ii].molar_mass);
    }
  }

  return mu;
}

/*====================== end of molar_mass() =================================*/

/*====================== mass_diffusivity() ==================================*/

/*
 * Calculate the binary-gas mass diffusivity [m^2/s] based on the
 * Fuller-Schettler-Giddings correlation, as described in 
 * Perry's Chemical Engineers' Handbook (1997), Table 5-14 and 5-16.
 */

EPIC_FLOAT mass_diffusivity(planetspec *planet,
                            int         vapor_index,
                            EPIC_FLOAT  temperature,
                            EPIC_FLOAT  pressure)
{
  static EPIC_FLOAT
    sumv_h2,sumv_he,
    sumv[LAST_INDEX+1];
  const EPIC_FLOAT
    one_third = 1./3.;
  static int
    initialized = FALSE;
  EPIC_FLOAT
    diff,tmp,sqrt_mab;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="mass_diffusivity";

  if (!initialized) {
    /*
     * Assign data. The parameter sumv [cm^3/mol] is given in Table 5-16 of Perry.
     * For molecules not listed, we sum the structural diffusion-volumes if available (for ozone,
     * we take three-halves of the O_2 value).
     *
     * NOTE: For critical applications, this calculation should be improved if possible.
     */
    sumv_h2            =  7.07;
    sumv_he            =  2.88;
    sumv[H_2O_INDEX]   = 12.7;
    sumv[NH_3_INDEX]   = 14.9;
    sumv[H_2S_INDEX]   = 1.98*2.+17.0;
    sumv[CH_4_INDEX]   = 16.5+1.98*4.;
    sumv[C_2H_2_INDEX] = 16.5*2.+1.98*2.;
    sumv[C_2H_4_INDEX] = 16.5*2.+1.98*4.;
    sumv[C_2H_6_INDEX] = 16.5*2.+1.98*6.;
    sumv[CO_2_INDEX]   = 26.9;
    sumv[NH_4SH_INDEX] = 5.69+1.98*4.+17.0+1.98;
    sumv[O_3_INDEX]    = 16.6*3./2.;
    sumv[N_2_INDEX]    = 17.9;
    sumv[PH_3_INDEX]   = 20.;   /* NOTE: 20. is a placeholder; need actual data for phosphine */
    /*
     * For planets where sumv[HDRY3_INDEX] has not
     * been measured, assume mole-fraction weighting.
     *
     * NOTE: Many of the mole fractions are rough and should be 
     *       updated when possible. 
     */
    if (strcmp(planet->type,"gas-giant") == 0) {
      sumv[HDRY2_INDEX] = planet->x_h2*sumv_h2+planet->x_he*sumv_he;
      sumv[HDRY3_INDEX] = planet->x_h2*sumv_h2+planet->x_he*sumv_he;
    }
    else if (strcmp(planet->type,"terrestrial") == 0) {
      switch(planet->index) {
        case VENUS_INDEX:
        case VENUS_LLR05_INDEX:
        case MARS_INDEX:
          sumv[HDRY2_INDEX] = sumv[CO_2_INDEX];
          sumv[HDRY3_INDEX] = sumv[CO_2_INDEX];
        break;
        case EARTH_INDEX:
        case HELD_SUAREZ_INDEX:
        case GOLDMINE_INDEX:
          /*  Goldmine is a testbed that uses the same numbers as Held_Suarez  */
          sumv[HDRY2_INDEX] = 20.1;
          sumv[HDRY3_INDEX] = 20.1;
        case TITAN_INDEX:
          sumv[HDRY2_INDEX] = sumv[N_2_INDEX];
          sumv[HDRY3_INDEX] = sumv[N_2_INDEX];
        break;
        default:
          sprintf(Message,"planet->name=%s not yet set up",planet->name);
          epic_error(dbmsname,Message);
        break;
      }
    }
    else {
      sprintf(Message,"planet->type=%s not yet set up",planet->type);
      epic_error(dbmsname,Message);
    }

    initialized = TRUE;
  }

  switch (vapor_index) {
    case HDRY3_INDEX:
    case H_2O_INDEX:
    case NH_3_INDEX:
    case H_2S_INDEX:
    case CH_4_INDEX:
    case C_2H_2_INDEX:
    case C_2H_4_INDEX:
    case C_2H_6_INDEX:
    case CO_2_INDEX:
    case NH_4SH_INDEX:
    case O_3_INDEX:
    case N_2_INDEX:
    case PH_3_INDEX:
      tmp = pow(sumv[HDRY3_INDEX],one_third)+pow(sumv[vapor_index],one_third);
    break;
    case HDRY2_INDEX:
      tmp = pow(sumv[HDRY2_INDEX],one_third)+pow(sumv[vapor_index],one_third);
    break;
    default:
      sprintf(Message,"vapor_index=%d not yet implemented",vapor_index);
      epic_error(dbmsname,Message);
    break;
  }

  sqrt_mab = sqrt(1./molar_mass(vapor_index)+planet->rgas/R_GAS);
  /* Convert pressure from Pa to atm. */
  pressure /= 1.0133e+5;
  /* Diffusivity is in cm^2/s. */
  diff = 0.001*pow(temperature,1.75)*sqrt_mab/(pressure*tmp*tmp);

  /* Convert to m^2/s. */
  return diff*1.e-4;
}

/*====================== end of mass_diffusivity() ===========================*/

/*======================= timeplane_bookkeeping() ===========================*/

/*
 * Set the time index for the current time for the prognostic variables
 * and their tendencies, depending on the chosen time-marching algorithms.
 */

void timeplane_bookkeeping(void)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="timeplane_bookkeeping";

  if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
    grid.it_uv      = 0;
    grid.it_uv_dis  = 0;
    grid.it_uv_tend = IT_ZERO;
  }
  else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
    grid.it_uv = IT_ZERO;
    if (U(IT_MINUS1,KLO,JLO,ILO) == FLOAT_MAX) {
      /* 
       * First timestep is a forward (Euler) step.
       */
      grid.it_uv_dis = IT_ZERO;
    }
    else {
      /*
       * Dissipative tendencies are lagged for numerical stability.
       */
      grid.it_uv_dis = IT_MINUS1;
    }
    grid.it_uv_tend = 0;
  }
  else {
    sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
    epic_error(dbmsname,Message);
  }
  grid.it_h = 0;

  return;
}

/*======================= end of timeplane_bookkeeping() ====================*/

/*====================== check_nan() =========================================*/

/*
 * Debugging tool to screen for occurrances of nan (not-a-number).
 * Input NULL to get silent mode.
 */

#define CHECK_NAN_WIND(U) \
  for (K = KLOPAD; K <= KHIPAD; K++) { \
    for (J = JLOPAD; J <= JHIPADPV; J++) { \
      for (I = ILOPAD; I <= IHIPAD; I++) { \
        if (!isfinite(U(grid.it_uv,K,J,I))) { \
          fprintf(stderr,"%s(%2d,%2d,%2d,%2d)=%g ", \
                          #U,grid.it_uv,K,J,I,U(grid.it_uv,K,J,I));fflush(stderr); \
          have_nan = TRUE; \
        } \
        if (!isfinite(DUDT(grid.it_uv_tend,K,J,I))) { \
          fprintf(stderr,"D%sDT(%2d,%2d,%2d,%2d)=%g ", \
                          #U,grid.it_uv_tend,K,J,I,DUDT(grid.it_uv_tend,K,J,I));fflush(stderr); \
          have_nan = TRUE; \
        } \
      } \
    } \
  } \
  if (!have_nan && message) { \
    fprintf(stderr,"%s ",#U);fflush(stderr); \
  }


#define CHECK_NAN_H(H) \
  for (K = KLOPAD; K <= KHIPAD; K++) { \
    for (J = JLOPAD; J <= JHIPADPV; J++) { \
      for (I = ILOPAD; I <= IHIPAD; I++) { \
        if (!isfinite(H(K,J,I))) { \
          fprintf(stderr,"%s(%2d,%2d,%2d)=%g ", \
                       #H,K,J,I,H(K,J,I));fflush(stderr); \
          have_nan = TRUE; \
        } \
      } \
    } \
  } \
  if (!have_nan && message) { \
    fprintf(stderr,"%s ",#H);fflush(stderr); \
  }

#define CHECK_NAN_DIAG(T2) \
  for (K = KLOPAD; K <= KHIPAD; K++) { \
    for (J = JLOPAD; J <= JHIPADPV; J++) { \
      for (I = ILOPAD; I <= IHIPAD; I++) { \
        if (!isfinite(T2(K,J,I))) { \
          fprintf(stderr,"%s(%2d,%2d,%2d)=%g ", \
                       #T2,K,J,I,T2(K,J,I));fflush(stderr); \
          have_nan = TRUE; \
        } \
      } \
    } \
  } \
  if (!have_nan && message) { \
    fprintf(stderr,"%s ",#T2);fflush(stderr); \
  }

#define CHECK_NAN_SPECIES(is) \
  for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) { \
    if (var.species[is].phase[ip].on) { \
      for (K = KLOPAD; K <= KHIPAD; K++) { \
        for (J = JLOPAD; J <= JHIPADPV; J++) { \
          for (I = ILOPAD; I <= IHIPAD; I++) { \
            if (!isfinite(Q(is,ip,K,J,I))) { \
              fprintf(stderr,"%s(%2d,%2d,%2d)=%g ", \
                             var.species[is].phase[ip].info[MASS].name,K,J,I,Q(is,ip,K,J,I));fflush(stderr); \
              have_nan = TRUE; \
            } \
          } \
        } \
      } \
      if (!have_nan && message) { \
        fprintf(stderr,"%s ",var.species[is].phase[ip].info[MASS].name);fflush(stderr); \
      } \
    } \
  }

#define CHECK_NAN_2D(PHI_SURFACE) \
  for (J = JLOPAD; J <= JHIPADPV; J++) { \
    for (I = ILOPAD; I <= IHIPAD; I++) { \
      if (!isfinite(PHI_SURFACE(J,I))) { \
        fprintf(stderr,"%s(%2d,%2d)=%g ", \
                     #PHI_SURFACE,J,I,PHI_SURFACE(J,I));fflush(stderr); \
        have_nan = TRUE; \
      } \
    } \
  } \
  if (!have_nan && message) { \
    fprintf(stderr,"%s ",#PHI_SURFACE);fflush(stderr); \
  }

void check_nan(char *message)
{
  int
    K,J,I,
    is,ip,im,
    ip_first,ip_last,
    have_nan = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="check_nan";
#ifdef EPIC_MPI
  int
    mpi_itmp;
#endif

  if (IAMNODE == NODE0 && message) {
    fprintf(stderr,"\n\nIT=%lu %s: no nan:",grid.itime,message);fflush(stderr);
  }

  /*
   * Check prognostic variables.
   */
  if (var.u.on)     {CHECK_NAN_WIND(U)};
  if (var.v.on)     {CHECK_NAN_WIND(V)};
  if (var.h.on)     {CHECK_NAN_H(H)};
  if (var.theta.on) {CHECK_NAN_H(THETA)};
  if (var.nu_turb.on) {CHECK_NAN_H(NU_TURB)};
  if (var.fpara.on) {CHECK_NAN_H(FPARA)};
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {CHECK_NAN_SPECIES(is)};
  }

  /*
   * Check 3D diagnostic variables.
   */
  if (var.hdry2.on)                {CHECK_NAN_DIAG(HDRY2)};
  if (var.hdry3.on)                {CHECK_NAN_DIAG(HDRY3)};
  if (var.h3.on)                   {CHECK_NAN_DIAG(H3)};
  if (var.pdry3.on)                {CHECK_NAN_DIAG(PDRY3)};
  if (var.p2.on)                   {CHECK_NAN_DIAG(P2)};
  if (var.p3.on)                   {CHECK_NAN_DIAG(P3)};
  if (var.theta2.on)               {CHECK_NAN_DIAG(THETA2)};
  if (var.t2.on)                   {CHECK_NAN_DIAG(T2)};
  if (var.t3.on)                   {CHECK_NAN_DIAG(T3)};
  if (var.rho2.on)                 {CHECK_NAN_DIAG(RHO2)};
  if (var.rho3.on)                 {CHECK_NAN_DIAG(RHO3)};
  if (var.exner2.on)               {CHECK_NAN_DIAG(EXNER2)};
  if (var.exner3.on)               {CHECK_NAN_DIAG(EXNER3)};
  if (var.fgibb2.on)               {CHECK_NAN_DIAG(FGIBB2)};
  if (var.phi2.on)                 {CHECK_NAN_DIAG(PHI2)};
  if (var.phi3.on)                 {CHECK_NAN_DIAG(PHI3)};
  if (var.mont2.on)                {CHECK_NAN_DIAG(MONT2)};
  if (var.heat3.on)                {CHECK_NAN_DIAG(HEAT3)};
  if (var.pv2.on)                  {CHECK_NAN_DIAG(PV2)};
  if (var.ri2.on)                  {CHECK_NAN_DIAG(RI2)};
  if (var.div_uv2.on)              {CHECK_NAN_DIAG(DIV_UV2)};
  if (var.w3.on)                   {CHECK_NAN_DIAG(W3)};
  if (var.z2.on)                   {CHECK_NAN_DIAG(Z2)};
  if (var.z3.on)                   {CHECK_NAN_DIAG(Z3)};
  if (var.dzdt2.on)                {CHECK_NAN_DIAG(DZDT2)};
  if (var.diffusion_coef_mass.on)  {CHECK_NAN_DIAG(DIFFUSION_COEF_MASS)};
  if (var.diffusion_coef_theta.on) {CHECK_NAN_DIAG(DIFFUSION_COEF_THETA)};
  if (var.diffusion_coef_uv.on)    {CHECK_NAN_DIAG(DIFFUSION_COEF_UV)};

  /*
   * Check 2D diagnostic variables.
   */
  if (var.phi_surface.on) {CHECK_NAN_2D(PHI_SURFACE)};

  if (var.pbot.on) {CHECK_NAN_2D(PBOT)};

  if (IAMNODE == NODE0 && message) {
    fprintf(stderr,"\n");fflush(stderr);
  }

#ifdef EPIC_MPI
  /*
   * Make status of have_nan global.
   */
  mpi_itmp = have_nan;
  MPI_Allreduce(&mpi_itmp,&have_nan,1,MPI_INT,MPI_SUM,para.comm);
#endif

  if (have_nan && IAMNODE == NODE0) {
    if (message) {
      sprintf(Message,"%s: aborting run",message);
    }
    else {
      sprintf(Message,"aborting run");
    }
    epic_error(dbmsname,Message);
  }

  return;
}

/*====================== end of check_nan() ==================================*/

/*======================= u_venus() ==========================================*/

/* 
 * Calculate u(p,lat) for Venus model, where lat is in degrees.
 */

EPIC_FLOAT u_venus(EPIC_FLOAT p,
                   EPIC_FLOAT lat)
{
  EPIC_FLOAT
    tmp,
    u,
    u0;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="u_venus";

  u0  = -117.374;
  tmp = cos(lat*DEG);
  u   = u_amp(planet,p)*u0*tmp*tmp;

  return u;
}

/*======================= end of u_venus() ==================================*/

/*======================= u_earth() =========================================*/

/*
 * Set u(p,lat) for Earth, where lat is in degrees.
 */

EPIC_FLOAT u_earth(EPIC_FLOAT p,
                   EPIC_FLOAT lat)
{
  int 
    j;
  static int
    ndat,
    initialized = FALSE;
  EPIC_FLOAT 
    u,
    lat_d;
  EPIC_FLOAT
    *latdat,
    *udat;
  static float_triplet
    *u_table;
  char
    header[N_STR];
  FILE
    *u_dat;
    
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="u_earth";

  if (!initialized) {
    if (IAMNODE == NODE0) {
      /* Look in local directory first. */
      u_dat = fopen("./u_vs_lat.Earth","r");
      if (!u_dat) {
        u_dat = fopen(EPIC_PATH"/data/Earth/u_vs_lat.Earth","r");
      }
      if (!u_dat) {
        sprintf(Message,"Failed to open file %s",EPIC_PATH"/data/Earth/u_vs_lat.Earth");
        epic_error(dbmsname,Message);
      }
      /* Skip 6-line header. */
      for (j = 0; j < 6; j++) {
        fgets(header,N_STR,u_dat);
      }
      fscanf(u_dat,"%d",&ndat);
    }

    /* Allocate memory. */
    latdat  = fvector( 0,ndat-1,dbmsname);
    udat    = fvector( 0,ndat-1,dbmsname);
    u_table = ftriplet(0,ndat-1,dbmsname);

    for (j = 0; j < ndat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(u_dat,"%lf %lf",latdat+j,udat+j);
#else
      fscanf(u_dat,"%f %f",latdat+j,udat+j);
#endif

    }
    fclose(u_dat);

    for (j = 0; j < ndat; j++) {
      u_table[j].x = latdat[j];
      u_table[j].y = udat[j];
    }

    /* Free allocated memory. */
    free_fvector(latdat,0,ndat-1,dbmsname);
    free_fvector(udat,  0,ndat-1,dbmsname);

    spline_pchip(ndat,u_table);

    initialized = TRUE;
  }
  /* End initialization */
  
  j = find_place_in_table(ndat,u_table,lat,&lat_d);
  u = splint_pchip(lat,u_table+j,lat_d);

  return u_amp(planet,p)*u;
}

/*======================= end of u_earth() ==================================*/

/*======================= u_mars() ==========================================*/

/*
 * Set u(p,lat) for Mars, where lat is in degrees.
 */

EPIC_FLOAT u_mars(EPIC_FLOAT p,
                  EPIC_FLOAT lat)
{
  static int
    initialized = FALSE;
  EPIC_FLOAT
    u;
 
  if (!initialized){
    fprintf(stderr,"Note: u_mars(): currently u(p,lat) = 0. \n");
    initialized = TRUE;
  }

  u = 0.;

  return u;
}

/*======================= end of u_mars() ===================================*/

/*======================= u_jupiter() =======================================*/

/*  
 *  Calculate u(p,lat) for Jupiter, where lat is in degrees.
 */

EPIC_FLOAT u_jupiter(EPIC_FLOAT p,
                     EPIC_FLOAT lat)
{
  int 
    j;
  static int
    ndat,
    initialized = FALSE;
  EPIC_FLOAT 
    u,
    lat_d;
  EPIC_FLOAT
    *latdat,
    *udat;
  static float_triplet
    *u_table;
  char
    header[N_STR];
  FILE
    *u_dat;
    
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="u_jupiter";

  if (!initialized) {
    /* Look in local directory first. */
    u_dat = fopen("./u_vs_lat.Jupiter","r");
    if (!u_dat) {
      u_dat = fopen(EPIC_PATH"/data/Jupiter/u_vs_lat.Jupiter","r");
    }
    if (!u_dat) {
      sprintf(Message,"Failed to open file %s",EPIC_PATH"/data/Jupiter/u_vs_lat.Jupiter");
      epic_error(dbmsname,Message);
    }
    /* Skip 6-line header. */
    for (j = 0; j < 6; j++) {
      fgets(header,N_STR,u_dat);
    }
    fscanf(u_dat,"%d",&ndat);

    /* Allocate memory. */
    latdat  = fvector( 0,ndat-1,dbmsname);
    udat    = fvector( 0,ndat-1,dbmsname);
    u_table = ftriplet(0,ndat-1,dbmsname);

    for (j = 0; j < ndat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(u_dat,"%lf %lf",latdat+j,udat+j);
#else
      fscanf(u_dat,"%f %f",latdat+j,udat+j);
#endif

    }
    fclose(u_dat);

    for (j = 0; j < ndat; j++) {
      u_table[j].x = latdat[j];
      u_table[j].y = udat[j];
    }

    /* Free allocated memory. */
    free_fvector(latdat,0,ndat-1,dbmsname);
    free_fvector(udat,  0,ndat-1,dbmsname);

    spline_pchip(ndat,u_table);

    initialized = TRUE;
  }
  /* End initialization */
  
  j = find_place_in_table(ndat,u_table,lat,&lat_d);
  u = splint_pchip(lat,u_table+j,lat_d);

  return u_amp(planet,p)*u;
}

/*======================= end of u_jupiter() ================================*/

/*======================= u_saturn() ========================================*/

/*  
 *  Calculate u(p,lat) for Saturn, where lat is in degrees.
 */

EPIC_FLOAT u_saturn(EPIC_FLOAT p,
                    EPIC_FLOAT lat)
{
  int 
    j;
  static int
    ndat,
    initialized = FALSE;
  EPIC_FLOAT 
    u,
    lat_d;
  EPIC_FLOAT
    *latdat,
    *udat;
  static float_triplet
    *u_table;
  char
    header[N_STR];
  FILE
    *u_dat;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="u_saturn";

  if (!initialized) {
    /* Look in local directory first. */
    u_dat = fopen("./u_vs_lat.Saturn","r");
    if (!u_dat) {
      u_dat = fopen(EPIC_PATH"/data/Saturn/u_vs_lat.Saturn","r");
    }
    if (!u_dat) {
      sprintf(Message,"Failed to open file %s",EPIC_PATH"/data/Saturn/u_vs_lat.Saturn");
      epic_error(dbmsname,Message);
    }
    /* Skip 6-line header. */
    for (j = 0; j < 6; j++) {
      fgets(header,N_STR,u_dat);
    }
    fscanf(u_dat,"%d",&ndat);

    /* Allocate memory. */
    latdat  = fvector( 0,ndat-1,dbmsname);
    udat    = fvector( 0,ndat-1,dbmsname);
    u_table = ftriplet(0,ndat-1,dbmsname);

    for (j = 0; j < ndat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(u_dat,"%lf %lf",latdat+j,udat+j);
#else
      fscanf(u_dat,"%f %f",latdat+j,udat+j);
#endif

    }
    fclose(u_dat);

    for (j = 0; j < ndat; j++) {
      u_table[j].x = latdat[j];
      u_table[j].y = udat[j];
    }

    /* Free allocated memory. */
    free_fvector(latdat,0,ndat-1,dbmsname);
    free_fvector(udat,  0,ndat-1,dbmsname);

    spline_pchip(ndat,u_table);

    initialized = TRUE;
  }
  /* End initialization */
  
  j = find_place_in_table(ndat,u_table,lat,&lat_d);
  u = splint_pchip(lat,u_table+j,lat_d);

  return u_amp(planet,p)*u;
}

/*======================= end of u_saturn() =================================*/

/*======================= u_titan() =========================================*/

/*
 * Set u(p,lat) for Titan, where lat is in degrees.
 */

EPIC_FLOAT u_titan(EPIC_FLOAT p,
                   EPIC_FLOAT lat)
{
  static int
    initialized = FALSE;
  EPIC_FLOAT
    u;
 
  if (!initialized){
    fprintf(stderr,"Note: u_titan(): currently u(p,lat) = 0. \n");
    initialized = TRUE;
  }

  u = 0.;

  return u;
}

/*======================= end of u_titan() ==================================*/

/*====================== u_uranus() =========================================*/

/*
 *  Calculate Uranus u(p,lat), where lat is in degrees.
 *
 *  Legendre polynomials (e.g. LeBeau and Dowling 1998, Icarus) fit to
 *  unbinned cloud measurment data from Sromovsky (2005, Icarus),
 *  by Michael Sussman.
 */

EPIC_FLOAT u_uranus(EPIC_FLOAT p,
                    EPIC_FLOAT lat)
{
  static int
    ndat,
    initialized = FALSE;
  int
    j;
  EPIC_FLOAT
    u,
    lat_d,
   *latdat,
   *udat;
  static float_triplet
    *u_table;
  char
    header[N_STR];
  FILE
    *u_dat;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
   static char
    dbmsname[]="u_uranus";

  if (!initialized) {
    /* Look in local directory first. */
    u_dat = fopen("./u_vs_lat.Uranus","r");
    if (!u_dat) {
      u_dat = fopen(EPIC_PATH"/data/Uranus/u_vs_lat.Uranus","r");
    }
    if (!u_dat) {
      sprintf(Message,"Failed to open file %s",EPIC_PATH"/data/Uranus/u_vs_lat.Uranus");
      epic_error(dbmsname,Message);
    }
    /* Skip 6-line header. */
    for (j = 0; j < 6; j++) {
      fgets(header,N_STR,u_dat);
    }
    fscanf(u_dat,"%d",&ndat);

    /* Allocate memory. */
    latdat  = fvector( 0,ndat-1,dbmsname);
    udat    = fvector( 0,ndat-1,dbmsname);
    u_table = ftriplet(0,ndat-1,dbmsname);

    for (j = 0; j < ndat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(u_dat,"%lf %lf",latdat+j,udat+j);
#else
      fscanf(u_dat,"%f %f",latdat+j,udat+j);
#endif

    }
    fclose(u_dat);

    for (j = 0; j < ndat; j++) {
      u_table[j].x = latdat[j];
      u_table[j].y = udat[j];
    }

    /* Free allocated memory. */
    free_fvector(latdat,0,ndat-1,dbmsname);
    free_fvector(udat,  0,ndat-1,dbmsname);

    spline_pchip(ndat,u_table);

    initialized = TRUE;
  }
  /* End initialization */

  j = find_place_in_table(ndat,u_table,lat,&lat_d);
  u = splint_pchip(lat,u_table+j,lat_d);

  return u_amp(planet,p)*u;
}

/*======================= end of u_uranus() =================================*/

/*====================== u_neptune() ========================================*/

/*
 * Calculate Neptune u(p,lat), where lat is in degrees.
 */

EPIC_FLOAT u_neptune(EPIC_FLOAT p,
                     EPIC_FLOAT lat)
{
  static int
    ndat,
    initialized = FALSE;
  int
    j;
  EPIC_FLOAT   
    u,
    lat_d,
   *latdat,
   *udat;
  static float_triplet
    *u_table;
  char
    header[N_STR];
  FILE
    *u_dat;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="u_neptune";

  if (!initialized) {
    /* Look in local directory first. */
    u_dat = fopen("./u_vs_lat.Neptune","r");
    if (!u_dat) {
      u_dat = fopen(EPIC_PATH"/data/Neptune/u_vs_lat.Neptune","r");
    }
    if (!u_dat) {
      sprintf(Message,"Failed to open file %s",EPIC_PATH"/data/Neptune/u_vs_lat.Neptune");
      epic_error(dbmsname,Message);
    }
    /* Skip 6-line header. */
    for (j = 0; j < 6; j++) {
      fgets(header,N_STR,u_dat);
    }
    fscanf(u_dat,"%d",&ndat);

    /* Allocate memory. */
    latdat  = fvector( 0,ndat-1,dbmsname);
    udat    = fvector( 0,ndat-1,dbmsname);
    u_table = ftriplet(0,ndat-1,dbmsname);

    for (j = 0; j < ndat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(u_dat,"%lf %lf",latdat+j,udat+j);
#else
      fscanf(u_dat,"%f %f",latdat+j,udat+j);
#endif

    }
    fclose(u_dat);

    for (j = 0; j < ndat; j++) {
      u_table[j].x = latdat[j];
      u_table[j].y = udat[j];
    }

    /* Free allocated memory. */
    free_fvector(latdat,0,ndat-1,dbmsname);
    free_fvector(udat,  0,ndat-1,dbmsname);

    spline_pchip(ndat,u_table);

    initialized = TRUE;
  }
  /* End initialization */

  j = find_place_in_table(ndat,u_table,lat,&lat_d);
  u = splint_pchip(lat,u_table+j,lat_d);

  return u_amp(planet,p)*u;
}

/*======================= end of u_neptune() ================================*/

/*======================= u_triton() ========================================*/

EPIC_FLOAT u_triton(EPIC_FLOAT p,
                    EPIC_FLOAT lat)
{
  EPIC_FLOAT
    u;

  u = u_amp(planet,p)*cos(lat*DEG);

  return u;
}

/*======================= end of u_triton() =================================*/

/*======================= u_pluto() =========================================*/

EPIC_FLOAT u_pluto(EPIC_FLOAT p,
                   EPIC_FLOAT lat)
{
  static int
    initialized = FALSE;
  EPIC_FLOAT
    u;
 
  if (!initialized){
    fprintf(stderr,"Note: u_pluto(): currently u = 0. \n");
    initialized = TRUE;
  }

  u = 0.;

  return u;
}

/*======================= end of u_pluto() ==================================*/

/*======================= u_hot_jupiter() ===================================*/

EPIC_FLOAT u_hot_jupiter(EPIC_FLOAT p,
                         EPIC_FLOAT lat)
{
  static int
    initialized = FALSE;
   EPIC_FLOAT
    u;

  if (!initialized) {
    fprintf(stderr,"Note: u_hot_jupiter(): currently u = 0. \n");
    initialized = TRUE;
  }

  u = 0.;

  return u;
}

/*======================= end of u_hot_jupiter() ============================*/

/*======================= u_null() ==========================================*/

EPIC_FLOAT u_null(EPIC_FLOAT p,
                  EPIC_FLOAT lat)
{
  EPIC_FLOAT
    u;
 
  u = 0.;

  return u;
}

/*======================= end of u_null() ===================================*/

/*======================= u_amp() ===========================================*/

/*
 * Returns nondimensional u(p) amplitude; units of p are Pa.
 * Useful for specifying the initial variation of zonal wind with pressure in a 
 * separable manner.
 *
 * NOTE: The actual vertical shears in a given planet are not likely to be separable,
 *       so this is just used as a convenient way to specify an initial condition 
 *       for the zonal wind.
 *
 * NOTE: The case grid.du_vert == 0. returns the pure barotropic value, 1.
 */

EPIC_FLOAT u_amp(planetspec *planet,
                 EPIC_FLOAT  p)
{
  EPIC_FLOAT
    p0,u;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="u_amp";

  if (grid.du_vert == 0.) {
    return 1.;
  }
  else {
    switch(planet->index) {
      case VENUS_INDEX:
        p0 = 87.47*100.;
        /*
         * Pioneer Venus profile.
         * Normalize to 87.47 hPa (location of max u in profile).
         */
        return grid.du_vert*(pioneer_venus_u(p)/pioneer_venus_u(p0)-1.)+1.;
      break;
      case JUPITER_INDEX:
        p0 = 680.*100.;
        /* 
         *  p < 680 mb: u_amp is set to follow the thermal-wind decay 
         *     determined by Gierasch et al (1986, Icarus 67, 456-483).
         *
         *  p > 680 mb: u_amp is the Galileo Probe Doppler wind profile,
         *              normalized at 680 hPa and scaled by grid.du_vert.
         */
        if (p <= p0) {
          return galileo_u(p);
        }
        else {
          return grid.du_vert*(galileo_u(p)-1.)+1.;
        }
      break;
      case SATURN_INDEX:
        return grid.du_vert*(cassini_cirs_u(p)-1.)+1.;
      break;
      default:
        sprintf(Message,"planet=%s not yet implemented",planet->name);
        epic_error(dbmsname,Message);
      break;
    }
  }

  /* Should never get here.*/
  sprintf(Message,"should never get here");
  epic_error(dbmsname,Message);
  return 0.;
}

/*======================= end of u_amp() ====================================*/

/*====================== pioneer_venus_u() ==================================*/

EPIC_FLOAT pioneer_venus_u(EPIC_FLOAT pressure)
{
  char   
    header[N_STR],
    infile[N_STR];
  int
    nn;
  static int
    nup,       
    initialized = FALSE;
  EPIC_FLOAT
    neg_log_p,
    u,
    p_up_d;
  EPIC_FLOAT
    *pdat,
    *udat;
  static float_triplet
    *up_table; 
  FILE
    *u_vs_p; 
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="pioneer_venus_u";

  if (!initialized) {
    /* 
     * Read in u vs p data.
     */

    /* Look in local directory first. */
    sprintf(infile,"./u_vs_p.Venus.Pioneer");
    u_vs_p = fopen(infile,"r");
    if (!u_vs_p) {
      sprintf(infile,EPIC_PATH"/data/Venus/u_vs_p.Venus.Pioneer");
      u_vs_p = fopen(infile,"r");
    }
    if (!u_vs_p) {
      sprintf(Message,"Failed to open file %s",infile);
      epic_error(dbmsname,Message);
    }
    for (nn = 0; nn < 6; nn++) {
      fgets(header,128,u_vs_p);  
    }
    /* input number of data points */
    fscanf(u_vs_p,"%d",&nup);  

    /* Allocate memory. */
    pdat     = fvector( 0,nup-1,dbmsname);
    udat     = fvector( 0,nup-1,dbmsname);
    up_table = ftriplet(0,nup-1,dbmsname);

    /* In order of increasing sigmatheta. */
    for (nn = nup-1; nn >= 0; nn--) { 

#if EPIC_PRECISION == DOUBLE_PRECISION 
      fscanf(u_vs_p,"%lf %lf",pdat+nn,udat+nn);
#else
      fscanf(u_vs_p,"%f %f",  pdat+nn,udat+nn);
#endif

      /* convert from hPa to Pa */
      pdat[nn] *= 100.;
      /*
       * Spline on -log p.
       */
      pdat[nn] = -log(pdat[nn]);
    }
    fclose(u_vs_p);

    for (nn = 0; nn < nup; nn++) {
      up_table[nn].x = pdat[nn];
      up_table[nn].y = udat[nn];
    }
    /* Free allocated memory. */
    free_fvector(pdat,0,nup-1,dbmsname);
    free_fvector(udat,0,nup-1,dbmsname);

    spline_pchip(nup,up_table);

    initialized = TRUE;
  }
  /* End initialization. */
        
  /*
   *  Interpolate to get zonal wind:
   */
  neg_log_p = -log(pressure);
  /*
   * Test whether out of table range.
   */
  if (neg_log_p < up_table[0].x || neg_log_p > up_table[nup-1].x) {
    sprintf(Message,"pressure=%g bar out of up_table range [%g,%g] \n",
            pressure*1.e-5,exp(-up_table[nup-1].x)*1.e-5,exp(-up_table[0].x)*1.e-5);
    epic_error(dbmsname,Message);
  }
  else {
    nn = find_place_in_table(nup,up_table,neg_log_p,&p_up_d);
    u  = splint_pchip(neg_log_p,up_table+nn,p_up_d);
  }

  return u;
}

/*====================== end of pioneer_venus_u() ===========================*/

/*====================== galileo_u() ========================================*/

EPIC_FLOAT galileo_u(EPIC_FLOAT pressure) 
{
  char   
    header[N_STR],
    infile[N_STR];
  int
    nn;
  static int
    initialized=FALSE,
    nup;
  EPIC_FLOAT
    neg_log_p,neg_log_p0,
    u,u0,
    p_up_d;
  EPIC_FLOAT
    *pdat,
    *udat;
  static float_triplet
    *up_table;
  FILE
    *u_vs_p;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="galileo_u";

  if (!initialized) {
    /* Look in local directory first. */
    sprintf(infile,"./u_vs_p.Jupiter.GalileoProbe");
    u_vs_p = fopen(infile,"r");
    if (!u_vs_p) {
      sprintf(infile,EPIC_PATH"/data/Jupiter/u_vs_p.Jupiter.GalileoProbe");
      u_vs_p = fopen(infile,"r");
    }
    if (!u_vs_p) {
      sprintf(Message,"Failed to open file %s",infile);
      epic_error(dbmsname,Message);
    }
    for (nn = 0; nn < 7; nn++) {
      fgets(header,100,u_vs_p); 
    }
    /* input number of data points */
    fscanf(u_vs_p,"%d",&nup); 
    for (nn = 0; nn < 4; nn++) {
      fgets(header,100,u_vs_p); 
    }

    /* Allocate memory: */
    pdat     = fvector( 0,nup-1,dbmsname);
    udat     = fvector( 0,nup-1,dbmsname);
    up_table = ftriplet(0,nup-1,dbmsname);

    /* In order of increasing sigmatheta. */
    for (nn = nup-1; nn >= 0;  nn--) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(u_vs_p, "%*f %*f %lf %lf",pdat+nn,udat+nn);
#else
      fscanf(u_vs_p, "%*f %*f %f %f",pdat+nn,udat+nn);
#endif

      /* Convert from bar to Pa. */
      pdat[nn] *= 1.e+5;
    }
    fclose(u_vs_p);

    for (nn = 0; nn < nup; nn++) {
      /* spline on neg log p */
      up_table[nn].x = -log(pdat[nn]);
      up_table[nn].y = udat[nn];
    }
    /* Free allocated memory. */
    free_fvector(pdat,0,nup-1,dbmsname);
    free_fvector(udat,0,nup-1,dbmsname);

    spline_pchip(nup,up_table);

    initialized = TRUE;
  }
  /* End of initialization. */

  /*
   *  Interpolate to get zonal wind:
   */
  neg_log_p  = -log(pressure);
  neg_log_p0 = -log(680.*100.);

  nn = find_place_in_table(nup,up_table,neg_log_p0,&p_up_d);

  /* Normalization factor */
  u0  = splint_pchip(neg_log_p0,up_table+nn,p_up_d);

  if (neg_log_p > neg_log_p0) {
    /* 
     * Thermal-wind decay determined by
     * Gierasch et al (1986, Icarus 67, 456-483).
     */
    u = 1.0 + (1.0/2.4) * log(pressure/(680.*100.));
    if (u < 0.) {
      u = 0.;
    }
  }
  else if (neg_log_p < neg_log_p0) {
    /*
     * Handle cases that are out of the table's range
     * by using the appropriate end-member value.
     */
    if (neg_log_p < up_table[0].x) {
      u = up_table[0].y/u0;
    }
    else if (neg_log_p > up_table[nup-1].x) {
      u = up_table[nup-1].y/u0;
    }
    else {
      nn = find_place_in_table(nup,up_table,neg_log_p,&p_up_d);
      u  = splint_pchip(neg_log_p,up_table+nn,p_up_d)/u0;
    }
  }
  else {
    u  = 1.0;
  }
  return u;
}

/*====================== end of galileo_u() =================================*/

/*====================== cassini_cirs_u() ===================================*/

/*
 * Cassini CIRS profile for Saturn: dashed line in Fig. 12 of 
 * Perez-Hoyos and Sanchez-Lavega (2006, Icarus 180, 161-175)
 * normalized to unity at p_ref.
 * Truncates to yield u >= 0.;
 * Assumes no shear for pressures deeper than p1.
 */
EPIC_FLOAT cassini_cirs_u(EPIC_FLOAT pressure)
{
  EPIC_FLOAT
    u;
  const EPIC_FLOAT
    u0    =  290.,
    dudh0 =   12.6,
    dudh1 =   43.,
    p0    =  100.*100.,
    p_ref =  500.*100.,
    p1    = 1000.*100.;
  static EPIC_FLOAT
    u_ref;
  static int
    initialized = FALSE;

  if (!initialized) {
    u_ref = u0+dudh1*log(p_ref/p0);

    initialized = TRUE;
  }

  if (pressure <= p0) {
    u = u0+dudh0*log(pressure/p0);
  }
  else if (pressure <= p1) {
    u = u0+dudh1*log(pressure/p0);
  }
  else {
    /*
     * Assume no shear for p > p1.
     */
    u = u0+dudh1*log(p1/p0);
  }

  return MAX(u/u_ref,0.);
}

/*====================== end of cassini_cirs_u() ============================*/

/*======================= p_sigmatheta() ====================================*/

/*
 * Returns the value of pressure consistent with 
 * theta, sigmatheta, pbot and ptop.
 */

EPIC_FLOAT p_sigmatheta(EPIC_FLOAT  theta,
                        EPIC_FLOAT  sigmatheta,
                        EPIC_FLOAT  pbot,
                        EPIC_FLOAT  hi_p,
                        EPIC_FLOAT  lo_p)
{
  int
    error_flag;
  const EPIC_FLOAT
    sgtol = pow(machine_epsilon(),2./3.);
  EPIC_FLOAT
    sg_root,
    p;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="p_sigmatheta";

  switch(grid.coord_type) {
    case COORD_HYBRID:
      SGTHMSGTH_theta      = theta;
      SGTHMSGTH_sigmatheta = sigmatheta;

      error_flag = find_root(get_sigma(pbot,hi_p),get_sigma(pbot,lo_p),sgtol,&sg_root,sgth_minus_sgth);
      if (error_flag) {
        /* 
         * The input range did not bracket the root, so before giving up, try the full range.
         */
        error_flag = find_root(0.,1.,sgtol,&sg_root,sgth_minus_sgth);
      }
      if (error_flag) {
        /*
         * Signal error by returning negative pressure.
         */
        return -1.;
      }

      p = get_p_sigma(pbot,sg_root);
    break;
    case COORD_ISOBARIC:
      p = get_p_sigma(pbot,sigmatheta);
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  return p;
}

/*======================= end of p_sigmatheta() =============================*/

/*======================= sgth_minus_sgth() =================================*/

EPIC_FLOAT sgth_minus_sgth(EPIC_FLOAT sigma)
{
  EPIC_FLOAT
    theta,sigmatheta;

  theta      = SGTHMSGTH_theta;
  sigmatheta = SGTHMSGTH_sigmatheta;

  return (sigmatheta-f_sigma(sigma)-g_sigma(sigma)*theta);
}

/*======================= end of sgth_minus_sgth_p() ========================*/

/*======================= turn_on_phases() ==================================*/

/*
 * Turn on phases appropriate to choice of cloud microphysics package.
 * The phases are switched on, whether or not the species are invoked.
 *
 * The Palotai & Dowling cloud microphysics scheme uses five phases for each species.
 *
 */
void turn_on_phases(planetspec *planet)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="turn_on_phases";

  if (grid.cloud_microphysics != OFF) {
    switch(planet->index) {
      case JUPITER_INDEX:
      case SATURN_INDEX:
        var.species[CH_4_INDEX].phase[VAPOR ].on = TRUE;
        var.species[CH_4_INDEX].phase[LIQUID].on = FALSE;
        var.species[CH_4_INDEX].phase[ICE   ].on = FALSE;
        var.species[CH_4_INDEX].phase[RAIN  ].on = FALSE;
        var.species[CH_4_INDEX].phase[SNOW  ].on = FALSE;

        var.species[NH_3_INDEX].phase[VAPOR ].on = TRUE;
        var.species[NH_3_INDEX].phase[LIQUID].on = TRUE;
        var.species[NH_3_INDEX].phase[ICE   ].on = TRUE;
        var.species[NH_3_INDEX].phase[RAIN  ].on = TRUE;
        var.species[NH_3_INDEX].phase[SNOW  ].on = TRUE;

        var.species[H_2O_INDEX].phase[VAPOR ].on = TRUE;
        var.species[H_2O_INDEX].phase[LIQUID].on = TRUE;
        var.species[H_2O_INDEX].phase[ICE   ].on = TRUE;
        var.species[H_2O_INDEX].phase[RAIN  ].on = TRUE;
        var.species[H_2O_INDEX].phase[SNOW  ].on = TRUE;
      break;
      case URANUS_INDEX:
      case NEPTUNE_INDEX:
        var.species[CH_4_INDEX].phase[VAPOR ].on = TRUE;
        var.species[CH_4_INDEX].phase[LIQUID].on = TRUE;
        var.species[CH_4_INDEX].phase[ICE   ].on = TRUE;
        var.species[CH_4_INDEX].phase[RAIN  ].on = TRUE;
        var.species[CH_4_INDEX].phase[SNOW  ].on = TRUE;

        var.species[H_2O_INDEX].phase[VAPOR ].on = TRUE;
        var.species[H_2O_INDEX].phase[LIQUID].on = TRUE;
        var.species[H_2O_INDEX].phase[ICE   ].on = TRUE;
        var.species[H_2O_INDEX].phase[RAIN  ].on = TRUE;
        var.species[H_2O_INDEX].phase[SNOW  ].on = TRUE;
      break;
      default:
        sprintf(Message,"need to implement cloud phases for planet->index=%d",planet->index);
        epic_error(dbmsname,Message);
      break;
    }
  }

  var.species[C_2H_2_INDEX].phase[VAPOR ].on = TRUE;
  var.species[C_2H_2_INDEX].phase[LIQUID].on = FALSE;
  var.species[C_2H_2_INDEX].phase[ICE   ].on = FALSE;
  var.species[C_2H_2_INDEX].phase[RAIN  ].on = FALSE;
  var.species[C_2H_2_INDEX].phase[SNOW  ].on = FALSE;

  var.species[C_2H_4_INDEX].phase[VAPOR ].on = TRUE;
  var.species[C_2H_4_INDEX].phase[LIQUID].on = FALSE;
  var.species[C_2H_4_INDEX].phase[ICE   ].on = FALSE;
  var.species[C_2H_4_INDEX].phase[RAIN  ].on = FALSE;
  var.species[C_2H_4_INDEX].phase[SNOW  ].on = FALSE;

  var.species[C_2H_6_INDEX].phase[VAPOR ].on = TRUE;
  var.species[C_2H_6_INDEX].phase[LIQUID].on = FALSE;
  var.species[C_2H_6_INDEX].phase[ICE   ].on = FALSE;
  var.species[C_2H_6_INDEX].phase[RAIN  ].on = FALSE;
  var.species[C_2H_6_INDEX].phase[SNOW  ].on = FALSE;

  var.species[PH_3_INDEX].phase[VAPOR ].on = TRUE;
  var.species[PH_3_INDEX].phase[LIQUID].on = FALSE;
  var.species[PH_3_INDEX].phase[ICE   ].on = FALSE;
  var.species[PH_3_INDEX].phase[RAIN  ].on = FALSE;
  var.species[PH_3_INDEX].phase[SNOW  ].on = FALSE;

  return;
}

/* * * * * * * * * * * end of epic_funcs_diag.c  * * * * * * * * * * * * * * */
