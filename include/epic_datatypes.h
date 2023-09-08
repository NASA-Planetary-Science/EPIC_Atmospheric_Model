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

#ifndef EPIC_DATATYPES_H
#define EPIC_DATATYPES_H
/* * * * * * * * * * * * * epic_datatypes.h  * * * * * * * * * * * * * * * 
 *                                                                       *
 *       Timothy E. Dowling                                              *
 *                                                                       *
 *       Header file containing the EPIC model data and structure        *
 *       type definitions.                                               *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "epic_microphysics.h"

#if defined(EPIC_MPI)
#  include "mpi.h"
#  include "mpg.h"
#endif

/*
 * Set floating-point precision.
 */
#define SINGLE_PRECISION 4
#define DOUBLE_PRECISION 8

/*
 * Array types
 */
#define EPIC_FLOAT_ARRAY      1
#define FLOAT_TRIPLET_ARRAY   2

#if EPIC_PRECISION == DOUBLE_PRECISION
#  define EPIC_FLOAT double
#  define FLOAT_MAX  DBL_MAX
#  define FLOAT_MIN  DBL_MIN
#elif EPIC_PRECISION == SINGLE_PRECISION
#  define EPIC_FLOAT float
#  define FLOAT_MAX  FLT_MAX
#  define FLOAT_MIN  FLT_MIN
#else
#  error Unrecognized value for EPIC_PRECISION environment variable.
#endif

#define VAR_NM_SZ 64


/*
 * System index, ordered by mass, then by special cases (benchmarks).
 */
#define NO_INDEX             0
#define SUN_INDEX            1
#define JUPITER_INDEX        2
#define SATURN_INDEX         3
#define NEPTUNE_INDEX        4
#define URANUS_INDEX         5
#define EARTH_INDEX          6
#define VENUS_INDEX          7
#define MARS_INDEX           8
#define MERCURY_INDEX        9
#define GANYMEDE_INDEX      10
#define TITAN_INDEX         11
#define CALLISTO_INDEX      12
#define IO_INDEX            13
#define MOON_INDEX          14
#define EUROPA_INDEX        15
#define TRITON_INDEX        16
#define ERIS_INDEX          17
#define PLUTO_INDEX         18
#define HOT_JUPITER_INDEX   19
#define HELD_SUAREZ_INDEX   20
#define VENUS_LLR05_INDEX   21
#define GOLDMINE_INDEX      22

#define MAX_NUM_SYSTEMS     22   /* update as needed */

/* 
 * Core prognostic variables. 
 * The first valid index value should be zero.
 */
#define FIRST_INDEX     1
#define FIRST_PROG      FIRST_INDEX

#define U_INDEX         1
#define V_INDEX         2
#define H_INDEX         3
#define THETA_INDEX     4

#define NU_TURB_INDEX   5

/* Para hydrogen fraction */
#define FPARA_INDEX     6    

/* Species indices */
#define FIRST_SPECIES   7

#define H_2O_INDEX      7
#define NH_3_INDEX      8
#define H_2S_INDEX      9
#define CH_4_INDEX     10
#define C_2H_2_INDEX   11
#define C_2H_4_INDEX   12
#define C_2H_6_INDEX   13
#define CO_2_INDEX     14
#define NH_4SH_INDEX   15
#define O_3_INDEX      16
#define N_2_INDEX      17
#define PH_3_INDEX     18

#define LAST_SPECIES   18

#define MAX_NUM_SPECIES (LAST_SPECIES-FIRST_SPECIES+1) 

#define LAST_PROG        LAST_SPECIES

#define MAX_NUM_PROGS   (LAST_PROG-FIRST_PROG+1)

/*
 * Diagnostic variables. 
 *
 * Suffix "2" is a variable carried in the layer.
 * Suffix "3" is a variable carried on the lower interface.
 */
#define H3_INDEX                   (LAST_PROG+ 1)
#define HDRY2_INDEX                (LAST_PROG+ 2)
#define HDRY3_INDEX                (LAST_PROG+ 3)
#define P2_INDEX                   (LAST_PROG+ 4)
#define P3_INDEX                   (LAST_PROG+ 5)
#define PDRY3_INDEX                (LAST_PROG+ 6)
#define THETA2_INDEX               (LAST_PROG+ 7)
#define T2_INDEX                   (LAST_PROG+ 8)
#define T3_INDEX                   (LAST_PROG+ 9)
#define RHO2_INDEX                 (LAST_PROG+10)
#define RHO3_INDEX                 (LAST_PROG+11)
#define EXNER2_INDEX               (LAST_PROG+12)
#define EXNER3_INDEX               (LAST_PROG+13)
#define FGIBB2_INDEX               (LAST_PROG+14)
#define PHI2_INDEX                 (LAST_PROG+15)
#define PHI3_INDEX                 (LAST_PROG+16)
#define MONT2_INDEX                (LAST_PROG+17)
#define HEAT3_INDEX                (LAST_PROG+18)
#define PV2_INDEX                  (LAST_PROG+19)
#define EDDY_PV2_INDEX             (LAST_PROG+20)
#define MOLAR_MASS3_INDEX          (LAST_PROG+21)
#define RI2_INDEX                  (LAST_PROG+22)
#define REL_VORT2_INDEX            (LAST_PROG+23)
#define EDDY_REL_VORT2_INDEX       (LAST_PROG+24)
#define ABS_VORT2_INDEX            (LAST_PROG+25)
#define KIN2_INDEX                 (LAST_PROG+26)
#define DIV_UV2_INDEX              (LAST_PROG+27)
#define W3_INDEX                   (LAST_PROG+28)
#define Z2_INDEX                   (LAST_PROG+29)
#define Z3_INDEX                   (LAST_PROG+30)
#define DZDT2_INDEX                (LAST_PROG+31)

/*
 * Turbulence-model variables.
 */
#define DIFFUSION_COEF_UV_INDEX    (LAST_PROG+32)
#define DIFFUSION_COEF_THETA_INDEX (LAST_PROG+33)
#define DIFFUSION_COEF_MASS_INDEX  (LAST_PROG+34)

/*
 * 3D parameters.
 */

/*
 * 2D parameters.
 */
#define PHI_SURFACE_INDEX          (LAST_PROG+35)
#define GRAVITY2_INDEX             (LAST_PROG+36)
#define PBOT_INDEX                 (LAST_PROG+37)

#define LAST_INDEX                 (LAST_PROG+37)

#define FILE_STR 256
#define GEOM_STR        16   /* geometry string length                          */
#define TOPDIM           3   /* 3 for x,y,z                                     */

/*
 * Dimension numbers.
 */
#define ONEDIM   1
#define TWODIM   2
#define THREEDIM 3
#define FOURDIM  4

#define N_STR        256
#define MAX_NU_ORDER 8

/*
 * Data structures.
 */
typedef struct {
  EPIC_FLOAT
    x,y;
} complex;

typedef struct {
  EPIC_FLOAT
    x,y;
} float_pair;

typedef struct {
  EPIC_FLOAT
    x,y,z;
} float_triplet;

typedef struct {
  EPIC_FLOAT
    e,w,n,s;
} float_quartet;

/*
 *  See epic_globals.c for details on planetspec members.
 */
typedef struct {
  int
    index;                  /* unique integer index, biggest first                */
  char
    name[32],               /* name of planet                                     */
    type[16],               /* gas-giant or terrestrial                           */
    orbital_epoch[8];       /* Epoch for Keplerian orbital elements               */
  EPIC_FLOAT  
    re,                     /* equatorial radius, m                               */
    rp,                     /* polar radius, m                                    */
    obliquity,              /* angle between rotational and orbital axes [deg]    */
    omega_sidereal,         /* angular vel. of rotation, 1/s, rel. to stars       */
    omega_synodic,          /* angular vel. of rotation, 1/s, rel. to solar day   */
    wlon_noon_J2000,        /* west longitude of mean-sun noon at J2000.0         */
    cp,                     /* specific heat at constant pressure                 */
    cpr,                    /* nondimensional reference cp, from thermo_setup()   */
    rgas,                   /* gas constant                                       */
    p0,                     /* ref. surface pressure [Pa], for potential temp.    */
    kappa,                  /* rgas/cp                                            */
    GM,                     /* gravitational constant times total mass [m^3/s^2]  */
    J2,                     /* gravitational zonal harmonic                       */
    x_h2,                   /* number fraction of molecular hydrogen              */
    x_he,                   /* number fraction of helium                          */
    x_3,                    /* number fraction of remaining constituents          */
    a,                      /* orbit semimajor axis [AU]                          */
    e,                      /* orbit eccentricity                                 */
    i,                      /* orbit inclination [deg]                            */
    lon_ascending_node,     /* longitude of ascending node [deg]                  */
    lon_perihelion,         /* longitude of perihelion [deg]                      */
    mean_lon,               /* orbit mean longitude [deg]                         */
    orbit_period,           /* orbit period, in years                             */
    vernal_equinox_anomaly, /* true anomaly of vernal equinox [deg]               */
    kinvisc,                /* typical laminar kinematic viscosity [m^2/s]        */
    dynvisc,                /* typical laminar dynamic viscosity [kg/m/s]         */
    k_a;                    /* thermal conductivity of dry air [J/m/s/K]          */
  EPIC_FLOAT
    (*u)(EPIC_FLOAT p, EPIC_FLOAT lat); /* zonal-wind profile                     */
  microphysics_spec
    cloud[MAX_NUM_SPECIES]; /* Top structure for cloud microphysics parameters    */
} planetspec;

/*
 * The gridspec structure contains model dimension and bookkeeping information.
 */
typedef struct {
  char
    geometry[GEOM_STR],
    vertical_coordinate[N_STR],
    advection_scheme[N_STR],
    uv_timestep_scheme[N_STR],
    radiation_scheme[N_STR],
    turbulence_scheme[N_STR],
    extract_append[N_STR];
  int
    coord_type;
  EPIC_FLOAT
    epic_version,
    globe_lonbot,
    globe_lontop,
    globe_latbot,
    globe_lattop;
  char
    f_plane_map[GEOM_STR];
  EPIC_FLOAT
    f_plane_lat0,
    f_plane_half_width;
  int
    dt,                            /* timestep, s                                        */
    cfl_dt,                        /* CFL timestep                                       */
    nk,                            /* number of vertical layers                          */
    nj,                            /* number of grid points in latitude                  */
    ni,                            /* number of grid points in longitude                 */
    jtp,                           /* j value associated with T(p) probe data            */
    wrap[TOPDIM],                  /* periodicity flags                                  */
    pad[TOPDIM],                   /* boundary pad widths                                */
    jlo,                           /* low-end index for j, 0 for globe geometry          */
    jfirst,                        /* used to handle staggered C-grid                    */
    jlast,                         /* used to handle staggered C-grid                    */
    ilo,                           /* low-end index for i, typically 1                   */
    we_num_nodes;                  /* number of nodes on computer running the model      */
  EPIC_FLOAT 
    dln,                           /* longitudinal grid spacing, deg                     */
    dlt,                           /* latitudinal  grid spacing, deg                     */
    sgth_bot,                      /* sigmatheta for bottom of model                     */
    sgth_top,                      /* sigmatheta for top of model                        */
    ptop,                          /* reference top pressure                             */
    pbot,                          /* reference bottom pressure                          */
    thetatop,                      /* reference top potential temperature [K]            */
    thetabot,                      /* reference bottom potential temperature [K]         */
    phi0;                          /* reference geopotential, where altitude z = 0       */
  double                           /* double to improve diag. theta calculation          */
    zeta0,                         /* zeta at sigma = 0, used in f_sigma()               */
    zeta1,                         /* zeta=theta at sigma = 1, used in f_sigma()         */
    hybrid_alpha,                  /* sigma-to-theta transition parameter                */
    sigma_sigma;                   /* sigma value below which vert. coord. is f(sigma)   */
  int
    k_sponge,                      /* number of top sponge layers                        */
    j_sponge,                      /* number of lateral sponge layers                    */
    k_sigma,                       /* top layer in pure-sigma region (g_sigma = 0)       */
    newt_cool_adjust,              /* 1 sets layer avg of Newtonian cooling to zero      */
    cloud_microphysics,            /* ACTIVE enables microphysics (phase changes, etc.)  */
    include_nontrad_accel,         /* TRUE includes non-trad Coriolis, spherical accels. */
    radiation_index,               /* index to distinguish radiation schemes             */
    zonal_average_rt,              /* TRUE zonally averages radiative transfer heating   */
    extract_species_fraction_type; /* mass mixing ratio or mole fraction                 */
  EPIC_FLOAT
    du_vert;                       /* used in u_amp() to set vert. profile of zonal wind */
  char
    eos[8];                        /* equation of state: "ideal", "virial"               */
  int    
    aux_a,                         /* for any use                                        */
    aux_b,                         /* for any use                                        */
    aux_c;                         /* for any use                                        */
  EPIC_FLOAT 
    aux_fa,                        /* for any use                                        */
    aux_fb,                        /* for any use                                        */
    aux_fc;                        /* for any use                                        */
  int
     nq,                           /* total number of active species-phases              */
    *is,                           /* array of active species indices                    */
    *ip;                           /* array of active phase indices                      */
  double                           /* double to improve diag. theta calculation          */
    *sigmatheta;                   /* hybrid sigma-theta array                           */ 
  EPIC_FLOAT
    *p_ref,                        /* typical pressure values, kk index                  */
    *t_ref,                        /* typical temperature values, kk index               */
    *rho_ref,                      /* typical density values [kg/m^3], kk index          */
    *theta_ref,                    /* typical theta values, kk index                     */
    *h_min,                        /* minimum layer thickness parameter                  */
    *re,                           /* equatorial radius of bottom interface of layer K   */
    *rp;                           /* polar radius of bottom interface of layer K        */
  EPIC_FLOAT
    **rln,                         /* longitudinal map factor, r                         */
    **rlt,                         /* latitudinal map factor, R                          */
    **m,                           /* longitudinal map factor, 1/dx = 1/(r*dln*DEG)      */
    **n,                           /* latitudinal  map factor, 1/dy = 1/(R*dlt*DEG)      */
    **mn,                          /* map factor, 1/d(area)                              */
    **g,                           /* gravity, m/s^2                                     */
    *f,                            /* Coriolis parameter, 2*Omega*sin(lat), 1/s          */
    *f2,                           /* Second Coriolis parameter, 2*Omega*cos(lat), 1/s   */
    **beta,                        /* beta = df/dy                                       */
    *lat,                          /* latitude,  deg                                     */
    *lon,                          /* longitude, deg                                     */
    *dsgth,                        /* differential of sigmatheta                         */
    *dsgth_inv;                    /* reciprical of dsgth                                */
  int   
    is_spole,                      /* south pole flag                                    */
    is_npole,                      /* north pole flag                                    */
    nu_order;                      /* Hyperviscosity order                               */
  EPIC_FLOAT 
    ab[3],                         /* Adams-Bashforth coefficients                       */
    nudiv_nondim,                  /* Divergence damping coefficient, nondimensional     */
    nu_nondim;                     /* Non-dimensional hyperviscosity coefficient         */
  double
    nu_hyper;                      /* dimensional value of hyperviscosity coefficient    */
  int   
    itback,                        /* num. of timesteps bet. backups to disk (same file) */
    itsave,                        /* num. of timesteps bet. outputs to disk (diff file) */
    itextract,                     /* num. of timesteps bet. writes to extract.nc        */
    itrun;                         /* num. of timesteps to run the model                 */
  int
    it_uv,                         /* time index for momentum variables u,v              */
    it_uv_dis,                     /* time index for u,v dissipation (allows lagging)    */ 
    it_uv_tend,                    /* time index for u,v tendency fields                 */
    it_h;                          /* time index for h-grid variables                    */
  unsigned long 
    itime;                         /* index for main time-incrementing loop              */
} gridspec;


/*
 * Structure to hold OpenMARS (formerly MACDA) .nc file information.
 * This is used when running "change -openmars".
 */
typedef struct {
  int
    nk,
    nj,
    ni,
    ntime;
  char
    project[64],
    infile_history[128];
  int
   *MY;
  EPIC_FLOAT
   *lon,
   *lat,
   *sigma,
   *time,
   *Ls,
   *ps,
   *tsurf,
   *temp,
   *u,
   *v,
   *theta;
} openmars_gridspec;

/*
 * Structure to hold EMARS (emars_v1.0_back_*.nc) file information.
 * This is used when running "change -emars".
 *
 * NOTE: Where OpenMARS and EMARS variables are the same except for different names,
 *       using EPIC-internal variable names adopted for openmars, as indicated below.
 *
 * NOTE: The EMARS surface geopotential is obtained from emars_v1.0_anal_mean_MY24_Ls090-120.nc,
 *       whereas for OpenMARS we have to constuct it via EPIC.
 */
typedef struct {
  int
    nk,
    nj,
    ni,
    ntime;
  int
   *MY;
  EPIC_FLOAT
   *lon,
   *lat,
   *ak,          /* emars uses a hybrid sigma-p coordinate               */
   *bk,          /*   "                             "                    */
   *time,        /* emars: macda_sol+mars_hour/24.                       */
   *Ls,
   *ps,
   *tsurf,       /* emars: ts                                            */
   *phi_surface, /* emars: Surface_geopotential, in emars_v1.0_anal_*.nc */
   *temp,        /* emars: t                                             */
   *u,
   *v,
   *p,
   *theta;
} emars_gridspec;


/*
 * Structure to hold Weizmann Institute gas-giant .nc file information.
 * This is used when running "change -weizmann".
 */
typedef struct {
  int
    nk,
    nj,
    ni,
    ntime;
  char
    title[64],
    infile_history[128];
  EPIC_FLOAT
   *lon,
   *lat,
   *p,
   *time,
   *temp,
   *u,
   *v,
   *rho,
   *theta;
} weizmann_gridspec;


/*
 * The thermospec structure is used in the thermodynamics routines
 * that have been adapted from Peter Gierasch's original Fortran routines.
 */

#define MDIM_THERMO   128
#define NDIM_THERMO    16

#define THLO_THERMO    20.
#define THHI_THERMO   600.
#define CCOLN_THERMO   -2.105769
#define CCPLN_THERMO   -1.666421
#define CPRH2           2.5       /* nondim low T ref. cp for H_2 */
#define CPRHE           2.5       /* nondim low T ref. cp for He */
#define CPR3            3.5       /* nondim low T ref. cp for non H_2,He component */

typedef struct {
  EPIC_FLOAT
    t_grid[MDIM_THERMO],
    theta_grid[MDIM_THERMO],
    array[5][MDIM_THERMO],
    theta_array[2][MDIM_THERMO],
    fpdat[NDIM_THERMO],
    t[NDIM_THERMO][MDIM_THERMO];
} thermospec;

typedef struct {
  int
    index,
    id,
    dim,
    dimid[FOURDIM],
    coorid[FOURDIM];
  char
    *name,
    *standard_name,
    *long_name,
    *units;
} id_information;

/*
 * Set indices for volatile-species phases.
 *
 * NOTE: FIRST_PHASE should correspond to VAPOR.
 *
 * NOTE: The string names ("rain," "snow") for additional phases need to be 
 *       added to the Phase_Names string array where it is declared in
 *       *.c files.
 *       
 */
#define NO_PHASE    -1
#define FIRST_PHASE  0

#define VAPOR        FIRST_PHASE
#define VAPOUR       VAPOR
#define GAS          VAPOR

#define LIQUID       1
#define CLOUD_LIQUID LIQUID

#define SOLID        2
#define ICE          SOLID
#define CLOUD_ICE    SOLID

#define RAIN         3
#define SNOW         4

#define LAST_PHASE   4

#define FIRST_NONPRECIP  VAPOR
#define LAST_NONPRECIP   SOLID
#define FIRST_PRECIP     RAIN
#define LAST_PRECIP      SNOW

#define MAX_NUM_PHASES (LAST_PHASE-FIRST_PHASE+1)

typedef struct {
  int
    on,
    extract_on;
    /*
     * q is mass mixing ratio [density_i/density_dry_air]
     * NOTE: Many textbooks use "q" for specific humidity [density_i/density_total]
     *       and "w" for mass mixing ratio, but w is also vertical velocity.
     *
     * x is number fraction [n_i/n_total], aka mole fraction or amount fraction
     * (and for ideal gases, the volume mixing ratio)
     * NOTE: The letter x is traditionally used for solids and liquids, and y for gases,
     *       but we use x for all three phases.
     */
  EPIC_FLOAT
    *q,
    *x;
  id_information
    info[2];
} phasespec;

typedef struct {
  int
    on,
    extract_on;
  EPIC_FLOAT
    *value,
    *tendency;
  id_information
    info[2],        /* Need one for each variable timeplane stored to disk. */
    info_tend[2];   /* Need one for each tendency timeplane stored to disk. */
} wind_variable;

typedef struct {
  int
    on,
    extract_on;
  EPIC_FLOAT
    *value;
  char
    advection_scheme[N_STR];
  id_information
    info[1];
} thermo_variable;

typedef struct {
  int
    on,
    extract_on;
  EPIC_FLOAT
    *value;
  id_information
    info[1];
} diagnostic_variable;

/*
 * Include in species_variable any information that is needed about
 * each chemical species (H_2O, NH_3, etc.).
 */
typedef struct {
  int
    on,
    extract_on;
  int
    HITRAN_index;
  phasespec
    phase[MAX_NUM_PHASES];
  EPIC_FLOAT
    molar_mass,
    triple_pt_t,
    triple_pt_p,
    Lf,
    Lv,
    Ls;    
  EPIC_FLOAT
    (*enthalpy_change)(int        init_phase,
                       int        final_phase,
                       EPIC_FLOAT temperature),
    (*sat_vapor_p)(EPIC_FLOAT temperature);
  char
    advection_scheme[N_STR];
  id_information
    info[1];
} species_variable;

typedef struct {
  EPIC_FLOAT
    value;
  id_information
    info;
} time_variable;

typedef struct {
  int
    on_list[LAST_INDEX+1],
    extract_on,
    extract_on_list[LAST_INDEX+1],
    ntp,
    n_t_cool;
  time_t
    start_time,
    model_time;
  size_t
    extract_time_index;
  wind_variable
    u,
    v;
  thermo_variable
    h,
    theta,
    fpara,
    nu_turb;
  /*
   * Since FIRST_SPECIES > 0, the memory for the species_variables
   * with indices 0 to FIRST_SPECIES-1 allocated next is not used.
   * However, we find this to be convenient, since it allows us to
   * use the species indices directly. (Note that this is different than
   * the way the cloud[] structure is treated in the microphysics code, for
   * which the macros convert the index "is" to "is-FIRST_SPECIES".)
   */
  species_variable
    species[LAST_SPECIES+1];
  /*
   * The following variables are diagnostic, meaning they are calculated
   * from the prognostic variables.
   */
  EPIC_FLOAT
    *pdat,                /* Pressure in sounding profile T(p), e.g. from t_vs_p.jupiter    */
    *tdat,                /* Temperature in sounding profile T(p), e.g. from t_vs_p.jupiter */
    *dtdat;               /* Temperature difference, used in some Newtonian cooling schemes */
  diagnostic_variable
    h3,                   /* total hybrid density in layer                                  */
    hdry2,hdry3,          /* hydrid density of dry air in layer, interface                  */
    pdry3,                /* partial pressure of dry air on interface                       */
    p2,p3,                /* pressure in layer, interface                                   */
    theta2,               /* potential temperature in layer                                 */
    t2,t3,                /* temperature in layer, interface                                */
    rho2,rho3,            /* density in layer, interface                                    */
    exner2,exner3,        /* exner=cp*T/theta in layer, interface                           */
    fgibb2,               /* ortho-para Gibbs function, F(T), in layer                      */
    phi2,                 /* geopotential in layer                                          */
    phi3,                 /* geopotential on interface                                      */
    mont2,                /* Montgomery potential in layer                                  */
    heat3,                /* heating rate on bottom interface of layer K                    */
    pv2,                  /* potential vorticity in layer                                   */
    eddy_pv2,             /* pv minus zonal average in layer                                */
    molar_mass3,          /* average molar mass [kg/kmol] on interface                      */
    ri2,                  /* Richardson number in layer                                     */
    rel_vort2,            /* relative vorticity in layer                                    */
    eddy_rel_vort2,       /* relative vorticty minus zonal average in layer                 */
    abs_vort2,            /* absolute vorticity in layer                                    */
    kinetic_energy2,      /* horizontal kinetic energy per mass in layer                    */
    div_uv2,              /* horizontal divergence of (u,v), layer                          */
    w3,                   /* hybrid vertical velocity at layer interface                    */
    z2,                   /* geopotential height of layer middle [m]                        */
    z3,                   /* geopotential height of bottom interface of layer [m]           */
    dzdt2,                /* regular vertical velocity, in m/s, in layer                    */
    diffusion_coef_uv,    /* turbulence-model diffusion coefficient for (u,v), in layer     */
    diffusion_coef_theta, /* turbulence-model diffusion coefficient for theta, on interface */
    diffusion_coef_mass,  /* turbulence-model diffusion coefficient for mass, in layer      */
    phi_surface,          /* geopotential at bottom of model                                */
    gravity2,             /* gravity as a function of layer and latitude [m/s^2]            */
    pbot;                 /* pressure at bottom of gas giant, used as a boundary condition  */
  float_triplet
    *t_cool_table;        /* Profile for time constant used in Newtonian cooling            */
  time_variable
    l_s;                  /* planetocentric solar longitude [deg]                           */ 
  id_information
    info[1];              /* used for the common time dimension                             */ 
  EPIC_FLOAT
    fpara_rate_scaling;   /* Nominal = 1.0                                                  */
} variablespec;

#if defined(EPIC_MPI)
typedef struct {
  int 
    iamnode,
    nproc,
    nelem2d,
    nelem3d,
    ndim,
    jfirst,
    jlast,
    is_npole,
    is_spole;
  int
    npad[TOPDIM],
    wrap[TOPDIM],
    nstart[TOPDIM],
    nend[TOPDIM],
    dimlen[TOPDIM],
    mylo[TOPDIM],
    myhi[TOPDIM],
    arraydims[TOPDIM],
    nprocs[TOPDIM];
  MPI_Comm 
    comm,
    comm_JLO,
    comm_ijk,
    comm_ij;
} mpispec;
mpispec para;
#endif

/*
 * Default-parameter structures:
 */
typedef struct {
  char  
    geometry[GEOM_STR],
    vertical_coordinate[N_STR],
    system_id[32],
    f_plane_map[GEOM_STR],
    eos[8],
    extract_str[N_STR],
    species_str[N_STR],
    layer_spacing_dat[N_STR],
    turbulence_scheme[N_STR];
  int 
    nk,nj,ni,dt,
    start_date_input_type,
    newt_cool_adjust,
    cloud_microphysics,
    on[LAST_SPECIES+1],
    spacing_type,
    coord_type,
    uv_timestep_scheme,
    radiation_index,
    zonal_average_rt,
    nu_order,
    k_sponge,
    j_sponge,
    extract_species_fraction_type;
  time_t
    start_time;
  struct tm
    UTC_start;
  EPIC_FLOAT
    globe_lonbot,  /* keep globe_lonbot as first floating-point variable */
    globe_lontop,globe_latbot,globe_lattop,
    f_plane_lat0,f_plane_half_width,
    lat_tp,
    ptop,pbot,p_sigma,
    thetatop,thetabot,
    nudiv_nondim,
    u_scale,du_vert,
    mole_fraction[LAST_SPECIES+1],
    mole_fraction_over_solar[LAST_SPECIES+1],
    rh_max[LAST_SPECIES+1],
    nu_nondim,
    fpara_rate_scaling;
} init_defaultspec;

/*
 * NOTE: When a value stored in epic.nc should be used as the prompt,
 *       there is no need to include it here in change_defaultspec.
 */
typedef struct {
  char
    infile[ N_STR],
    outfile[N_STR],
    extract_str[N_STR],
    openmars_infile[N_STR],
    emars_infile[N_STR],
    weizmann_infile[N_STR];
  int
    uv_timestep_scheme;
} change_defaultspec;

#define LAST_ATOMIC_NUMBER 103

typedef struct {
  int 
    atomic_number;
  char
    symbol[4],
    name[16];
  EPIC_FLOAT
    molar_mass,
    solar_abundance;
} chem_element;

typedef struct {
  char
    filename[100];
  int
    header_exists;
} write_stats_file_info;

/* * * * * * * * * *  end of epic_datatypes.h  * * * * * * * * * * * * * * */ 
#endif
