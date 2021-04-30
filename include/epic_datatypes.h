/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2009 Timothy E. Dowling                      *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License     *
 * as published by the Free Software Foundation; either version 2  *
 * of the License, or (at your option) any later version.          *
 * A copy of this License is in the file:                          *
 *   $EPIC4_PATH/License.txt                                       *
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
#define EPIC_FLOAT_ARRAY      0
#define FLOAT_TRIPLET_ARRAY   1

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

#define NO_INDEX       -1

/*
 * System index, ordered by mass.
 */
#define SUN_INDEX            0
#define JUPITER_INDEX        1
#define SATURN_INDEX         2
#define NEPTUNE_INDEX        3
#define URANUS_INDEX         4
#define EARTH_INDEX          5
#define VENUS_INDEX          6
#define MARS_INDEX           7
#define MERCURY_INDEX        8
#define GANYMEDE_INDEX       9
#define TITAN_INDEX         10
#define CALLISTO_INDEX      11
#define IO_INDEX            12
#define MOON_INDEX          13
#define EUROPA_INDEX        14
#define TRITON_INDEX        15
#define ERIS_INDEX          16
#define PLUTO_INDEX         17
#define HOT_JUPITER_INDEX   18
#define HELD_SUAREZ_INDEX   19
#define VENUS_LLR05_INDEX   20

#define MAX_NUM_SYSTEMS     21   /* update as needed */

/* 
 * Core prognostic variables. 
 * The first valid index value should be zero.
 */
#define FIRST_INDEX     0

#define U_INDEX         0
#define V_INDEX         1
#define HDRY_INDEX      2
#define THETA_INDEX     3

#define NU_TURB_INDEX   4

/* Para hydrogen fraction: */
#define FPARA_INDEX     5 

/* Species indices: */
#define FIRST_SPECIES   6
#define H_2O_INDEX      6
#define NH_3_INDEX      7
#define H_2S_INDEX      8
#define CH_4_INDEX      9
#define C_2H_2_INDEX   10
#define C_2H_6_INDEX   11
#define CO_2_INDEX     12
#define NH_4SH_INDEX   13
#define O_3_INDEX      14
#define N_2_INDEX      15
#define LAST_SPECIES   15
#define MAX_NUM_SPECIES (LAST_SPECIES-FIRST_SPECIES+1) 

#define LAST_PROG      15
#define MAX_NUM_PROGS  (LAST_PROG+1)

/*
 * Diagnostic variables. 
 *
 * A "2" at the end of the name implies the variable is located  
 * in the layer, and a  "3" implies it is carried on the lower 
 * interface of the layer.
 */
#define HDRY3_INDEX                (MAX_NUM_PROGS+ 0)
#define PDRY3_INDEX                (MAX_NUM_PROGS+ 1)
#define P2_INDEX                   (MAX_NUM_PROGS+ 2)
#define P3_INDEX                   (MAX_NUM_PROGS+ 3)
#define THETA2_INDEX               (MAX_NUM_PROGS+ 4)
#define H2_INDEX                   (MAX_NUM_PROGS+ 5)
#define H3_INDEX                   (MAX_NUM_PROGS+ 6)
#define T2_INDEX                   (MAX_NUM_PROGS+ 7)
#define T3_INDEX                   (MAX_NUM_PROGS+ 8)
#define RHO2_INDEX                 (MAX_NUM_PROGS+ 9)
#define RHO3_INDEX                 (MAX_NUM_PROGS+10)
#define EXNER3_INDEX               (MAX_NUM_PROGS+11)
#define FGIBB3_INDEX               (MAX_NUM_PROGS+12)
#define GZ3_INDEX                  (MAX_NUM_PROGS+13)
#define MONT3_INDEX                (MAX_NUM_PROGS+14)
#define HEAT3_INDEX                (MAX_NUM_PROGS+15)
#define PV3_INDEX                  (MAX_NUM_PROGS+16)
#define EDDY_PV3_INDEX             (MAX_NUM_PROGS+17)
#define RI2_INDEX                  (MAX_NUM_PROGS+18)
#define VORT3_INDEX                (MAX_NUM_PROGS+19)
#define DIV_UV3_INDEX              (MAX_NUM_PROGS+20)
#define W3_INDEX                   (MAX_NUM_PROGS+21)
#define DZDT3_INDEX                (MAX_NUM_PROGS+22)
/*
 * Turbulence-model variables.
 */
#define DIFFUSION_COEF_UV_INDEX    (MAX_NUM_PROGS+23)
#define DIFFUSION_COEF_THETA_INDEX (MAX_NUM_PROGS+24)
#define DIFFUSION_COEF_MASS_INDEX  (MAX_NUM_PROGS+25)
/*
 * 3D parameters.
 */
#define U_SPINUP_INDEX             (MAX_NUM_PROGS+26)
/*
 * 2D parameters.
 */
#define GZ_SURFACE_INDEX            (MAX_NUM_PROGS+27)
/*
 * NMT physics parameters
 */
#define DRY_ENTROPY_INDEX           (MAX_NUM_PROGS+28)
#define MOIST_ENTROPY_INDEX         (MAX_NUM_PROGS+29)
#define SAT_MOIST_ENTROPY_INDEX     (MAX_NUM_PROGS+30)
#define THE_FLUX_INDEX              (MAX_NUM_PROGS+31)
#define RT_FLUX_INDEX               (MAX_NUM_PROGS+32)
#define U_FLUX_INDEX                (MAX_NUM_PROGS+33)
#define V_FLUX_INDEX                (MAX_NUM_PROGS+34)
#define CONVTHROT_INDEX             (MAX_NUM_PROGS+35)
#define FLUXTHROT_INDEX             (MAX_NUM_PROGS+36)
#define RAIN_RATE_INDEX             (MAX_NUM_PROGS+37)

#define LAST_INDEX                  (MAX_NUM_PROGS+37)

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
typedef struct{
  EPIC_FLOAT
    x,y;
} complex;

typedef struct{
  EPIC_FLOAT
    x,y;
} float_pair;

typedef struct{
  EPIC_FLOAT
    x,y,z;
} float_triplet;

/*
 *  The planetspec structure contains the following members.
 *  See epic_globals.c for details.
 *    {index,
 *     name,type,
 *     orbital_epoch,
 *     re,rp,omega_sidereal,omega_synodic,
 *     cp,cpr=cp/rgas,rgas,kappa=1/cpr,g,
 *     x_h2,x_he,x_3,
 *     a[AU],e,i[deg],
 *     lon_ascending_node[deg],lon_perihelion[deg],mean_lon[deg],
 *     orbit_period[yrs],
 *     kinvisc[m^2/s],dynvisc[kg/m/s],k_a[J/m/s/K]
 *     u(p,lat[deg])}
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
    omega_sidereal,         /* angular vel. of rotation, 1/s, rel. to stars       */
    omega_synodic,          /* angular vel. of rotation, 1/s, rel. to solar day   */
    cp,                     /* specific heat at constant pressure                 */
    cpr,                    /* nondimensional reference cp, from thermo_setup()   */
    rgas,                   /* gas constant                                       */
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
    advection_scheme[N_STR],
    uv_timestep_scheme[N_STR],
    turbulence_scheme[N_STR],
    stability_factor[N_STR],
    extract_append[N_STR];
  EPIC_FLOAT
    data_version,
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
    dt,                   /* timestep, s                                        */
    cfl_dt,               /* CFL timestep                                       */
    nk,                   /* number of vertical layers                          */
    nj,                   /* number of grid points in latitude                  */
    ni,                   /* number of grid points in longitude                 */
    jtp,                  /* j value associated with T(p) probe data            */
    wrap[TOPDIM],         /* periodicity flags                                  */
    pad[TOPDIM],          /* boundary pad widths                                */
    jlo,                  /* lower index for j, 0 for globe geometry            */
    jfirst,               /* used to handle staggered C-grid                    */
    jlast,                /* used to handle staggered C-grid                    */
    we_num_nodes;         /* number of nodes on computer running the model      */
  EPIC_FLOAT 
    dln,                  /* longitudinal grid spacing, deg                     */
    dlt,                  /* latitudinal  grid spacing, deg                     */
    dy0,                  /* dy at lat = LAT0, m                                */
    sgth_bot,             /* sigmatheta for bottom of model                     */
    sgth_top,             /* sigmatheta for top of model                        */
    pbot,                 /* reference bottom pressure                          */
    press0,               /* reference pressure                                 */
    mont0;                /* reference Montgomery potential                     */
  double                  /* double to improve diag. theta calculation          */
    zeta0,                /* zeta at sigma = 0, used in f_sigma()               */
    zeta1,                /* zeta=theta at sigma = 1, used in f_sigma()         */
    hybrid_alpha,         /* sigma-to-theta transition parameter                */
    sigma_sigma;          /* sigma value below which vert. coord. is f(sigma)   */
  int
    k_sponge,             /* number of sponge layers                            */
    k_sigma,              /* top layer in pure-sigma region (g_sigma = 0)       */
    newt_cool_on,         /* 1 enables Newtonian cooling                        */
    newt_cool_adjust,     /* 1 sets layer avg of Newtonian cooling to zero      */
    microphysics_on,      /* 1 enables microphysics (phase changes, etc.)       */
    nmt_physics_on,       /* 1 enables nmt_physics (New Mexico Tech physics)    */
    diffusion_direction;  /* 1=horiz.+vert., 2=just horiz., 3=just vert.        */
  EPIC_FLOAT
    du_vert;              /* used in u_amp() to set vert. profile of zonal wind */
  char
    eos[8];               /* equation of state: "ideal", "virial"               */
  int    
    aux_a,                /* for any use                                        */
    aux_b,                /* for any use                                        */
    aux_c;                /* for any use                                        */
  EPIC_FLOAT 
    aux_fa,               /* for any use                                        */
    aux_fb,               /* for any use                                        */
    aux_fc;               /* for any use                                        */
  int
     nq,                  /* total number of active species-phases              */
    *is,                  /* array of active species indices                    */
    *ip;                  /* array of active phase indices                      */
  double                  /* double to improve diag. theta calculation          */
    *sigmatheta;          /* hybrid sigma-theta array                           */ 
  EPIC_FLOAT
    *p_ref,               /* typical pressure values, kk index                  */
    *theta_ref,           /* typical theta values, kk index                     */
    *h_min;               /* minimum layer thickness parameter                  */
  EPIC_FLOAT
    *rln,                 /* longitudinal map factor, r                         */
    *rlt,                 /* latitudinal map factor, R                          */
    *f,                   /* Coriolis parameter, 1/s                            */
    *m,                   /* longitudinal map factor, 1/(r*dln*DEG)             */
    *n,                   /* latitudinal  map factor, 1/(R*dlt*DEG)             */
    *mn,                  /* m*n                                                */
    *g,                   /* gravity vs latitutde, m/s^2                        */
    *lat,                 /* latitude,  deg                                     */
    *lon,                 /* longitude, deg                                     */
    *dsgth,               /* differential of sigmatheta                         */
    *dsgth_inv;           /* reciprical of dsgth                                */
  int   
    is_spole,             /* south pole flag                                    */
    is_npole,             /* north pole flag                                    */
    drag_v,               /* true if Rayleigh drag is to be applied to v        */
    drag_zonal_avg;       /* Indicates if Rayleigh drag operates on zonal avg   */
  EPIC_FLOAT 
    *t_sponge_inv,        /* Rayleigh damping time for sponge at top            */
    topo_scale,           /* artificial topography scaling                      */
    prandtl,              /* Ratio tau_rad/tau_drag                             */
    tau_drag,             /* Rayleigh drag timescale; nu0 = 1/tau_drag          */
    nudiv_nondim,         /* Divergence damping coefficient, nondimensional     */
    ab[3],                /* Adams-Bashforth coefficients                       */
    nu_nondim[MAX_NU_ORDER+1]; /*non-dimensional viscosity coefficients         */
  double
    nu[MAX_NU_ORDER+1];   /* viscosity coefficients, declared as doubles        */
  int   
    itback,               /* num. of timesteps bet. backups to disk (same file) */
    itsave,               /* num. of timesteps bet. outputs to disk (diff file) */
    itextract,            /* num. of timesteps bet. writes to extract.nc        */
    itrun;                /* num. of timesteps to run the model                 */
  int
    it_uv,                /* time index for momentum variables u,v              */
    it_uv_dis,            /* time index for u,v dissipation (allows lagging)    */ 
    it_uv_tend,           /* time index for u,v tendency fields                 */
    it_h;                 /* time index for h-grid variables                    */
  unsigned long 
    itime;                /* index for main time-incrementing loop              */
} gridspec;

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
  EPIC_FLOAT
    *q;
  id_information
    info[1];
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
  phasespec
    phase[MAX_NUM_PHASES];
  EPIC_FLOAT
    molar_mass,
    triple_pt_t,
    triple_pt_p,
    Lf,
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
  wind_variable
    u,
    v;
  thermo_variable
    hdry,
    theta,
    fpara,
    nu_turb;
  /*
   * Since FIRST_SPECIES > 0, the memory for the species_variables
   * with indices 0 to FIRST_SPECIES-1 allocated next is not used.
   * However, we find this to be convenient, since it allows us to
   * use the species indices directly.
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
    hdry3,                /* hydrid density on interface                                    */
    pdry3,                /* partial pressure of dry air on interface                       */
    p2,p3,                /* pressure in layer                                              */
    theta2,               /* potential temperature in layer                                 */
    h2,h3,                /* total hybrid density in layer, interface                       */
    t2,t3,                /* temperature in layer, interface                                */
    rho2,rho3,            /* density in layer, interface                                    */
    exner3,               /* exner=cp*T/theta interface                                     */
    fgibb3,               /* ortho-para Gibbs function, F(T), interface                     */
    gz3,                  /* geopotential on interface                                      */
    mont3,                /* Montgomery potential on interface                              */
    heat3,                /* heating rate on bottom interface of layer K                    */
    pv3,                  /* potential vorticity on interface                               */
    eddy_pv3,             /* pv3 minus zonal average                                        */
    ri2,                  /* Richardson number in layer                                     */
    vort3,                /* relative vorticity on interface                                */
    div_uv3,              /* horizontal divergence of (u,v), interface                      */
    w3,                   /* hybrid vertical velocity through layer interface               */
    dzdt3,                /* regular vertical velocity, in m/s, through layer interface     */
    u_spinup,             /* Rayleigh drag zonal-wind profile (k,j,i)                       */
    diffusion_coef_uv,    /* turbulence-model diffusion coefficient for (u,v), in layer     */
    diffusion_coef_theta, /* turbulence-model diffusion coefficient for theta, in layer     */
    diffusion_coef_mass,  /* turbulence-model diffusion coefficient for mass, in layer      */
    gz_surface,           /* surface elevation at bottom of model, multiplied by g          */
    dry_entropy,          /* nmt physics output (J/kg/K)                                    */
    moist_entropy,        /* nmt physics output (J/kg/K)                                    */
    sat_moist_entropy,    /* nmt physics output (J/kg/K)                                    */
    the_flux,             /* nmt physics: surface theta_e flux (K/kg/m^2)                   */ 
    rt_flux,              /* nmt physics: surface total water flux, evap. (g/g/kg/m^2)      */ 
    u_flux,               /* nmt physics: surface zonal momentum flux (m/s/kg/m^2)          */  
    v_flux,               /* nmt physics: surface meridional momentum flux (m/s/kg/m^2)     */  
    convthrot,            /* nmt physics: convective throttle value                         */  
    fluxthrot,            /* nmt physics: flux throttle value                               */  
    rain_rate;            /* nmt physics: rainfall rate (kg/m^2/s)                          */
  float_triplet
    *t_cool_table;        /* Profile for time constant used in Newtonian cooling            */
  time_variable
    l_s;                  /* planetocentric solar longitude [deg]                           */ 
  id_information
    info[1];              /* used for the common time dimension                             */ 
  EPIC_FLOAT
    time_fp_bar;          /* Nominal = 3.e+8*1000.*100. [sec*1bar, mks]                     */
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
    system_id[32],
    f_plane_map[GEOM_STR],
    eos[8],
    extract_str[N_STR],
    species_str[N_STR],
    layer_spacing_dat[N_STR];
  int 
    nk,nj,ni,dt,
    start_date_input_type,
    newt_cool_on,
    newt_cool_adjust,
    microphysics_on,
    nmt_physics_on,
    diffusion_direction,
    on[LAST_SPECIES+1],
    spacing_type,
    hgrid_mass_advection_scheme,
    hgrid_nonmass_advection_scheme,
    uv_timestep_scheme,
    turbulence_scheme,
    stability_factor,
    drag_v,
    drag_zonal_avg,
    k_sponge;
  int
    nmt_land,
    nmt_dorad,
    nmt_frad;
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
    prandtl,
    tau_drag,
    nudiv_nondim,
    u_scale,u_spinup_scale,du_vert,
    topo_scale,
    mole_fraction[LAST_SPECIES+1],
    mole_fraction_over_solar[LAST_SPECIES+1],
    rh_max[LAST_SPECIES+1],
    nu_nondim[MAX_NU_ORDER+1],
    time_fp;
  EPIC_FLOAT
    nmt_radcool,
    nmt_tpause,
    nmt_radbrk,
    nmt_cvc,
    nmt_cvs,
    nmt_cvp,
    nmt_cve,
    nmt_theslop,
    nmt_pstiff,
    nmt_pscale,
    nmt_pbltop,
    nmt_cdrag,
    nmt_wscale,
    nmt_sst,
    nmt_lfrac,
    nmt_eflux0,
    nmt_cld,
    nmt_cfract;
} init_defaultspec;

/*
 * NOTE: When a value stored in epic.nc should be used as the prompt,
 *       there is no need to include it here in change_defaultspec.
 */
typedef struct {
  char
    infile[ N_STR],
    outfile[N_STR],
    extract_str[N_STR];
  int
    hgrid_mass_advection_scheme,
    hgrid_nonmass_advection_scheme,
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
