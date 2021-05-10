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

/* * * * * * * * *  epic_funcs_diag.c  * * * * * * * * * * * * * * * * * * * *
 *                                                                           *
 *       Timothy E. Dowling                                                  *
 *                                                                           *
 *       Functions that reference data or EPIC model variables               *
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
 *           set_sponge()                                                    *
 *           get_sigma()                                                     *
 *           get_p_sigma()                                                   *
 *           get_h()                                                         *
 *           get_p()                                                         *
 *           alt_get_p()                                                     *
 *           state_from_exner()                                              *
 *           molar_mixing_ratio()                                            *
 *           get_var()                                                       *
 *           get_var_mean2d()                                                *
 *           onto_kk()                                                       *
 *           get_kin()                                                       *
 *           get_brunt2()                                                    *
 *           get_richardson()                                                *
 *           get_sounding()                                                  *
 *           return_cp()                                                     *
 *           set_p2etc()                                                     *
 *           store_pgrad_vars()                                              *
 *           store_diag()                                                    *
 *           divergence()                                                    *
 *           vorticity()                                                     *
 *           gz_from_u()                                                     *
 *           relative_humidity()                                             *
 *           source_sink()                                                   *
 *           cfl_dt()                                                        *
 *           time_mod()                                                      *
 *           b_vir(),b1_vir(),b2_vir()                                       *
 *           sum_xx()                                                        *
 *           avg_molar_mass()                                                *
 *           molar_mass()                                                    *
 *           diffusivity()                                                   *
 *           parse_species_name()                                            *
 *           solar_fraction()                                                *
 *           thermo_setup()                                                  *
 *           return_temp()                                                   *
 *           alt_return_temp()                                               *
 *           return_density()                                                *
 *           return_theta()                                                  *
 *           return_press()                                                  *
 *           return_enthalpy()                                               *
 *           return_fpe()                                                    *
 *           enthalpy()                                                      *
 *           timeplane_bookkeeping()                                         *
 *           check_periodic()                                                *
 *           check_nan()                                                     *
 *           u_venus, u_jupiter(), etc.                                      *
 *           u_amp()                                                         *
 *           galileo_u()                                                     *
 *           pioneer_venus_u()                                               *
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
#define SET_SPECIES(ispecies,the_species,the_standard_name,the_long_name,the_units) \
          var.species[ispecies].info[0].index = ispecies; \
          var.species[ispecies].info[0].name = (char *)malloc(strlen(#the_species)+1); \
          strcpy(var.species[ispecies].info[0].name,#the_species); \
          var.species[ispecies].info[0].standard_name = (char *)malloc(strlen(#the_standard_name)+1); \
          strcpy(var.species[ispecies].info[0].standard_name,#the_standard_name); \
          var.species[ispecies].info[0].long_name = (char *)malloc(strlen(#the_long_name)+1); \
          strcpy(var.species[ispecies].info[0].long_name,#the_long_name); \
          for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) { \
            var.species[ispecies].phase[ip].info[0].index = ispecies*1000000+(ip+1)*1000; \
            var.species[ispecies].phase[ip].info[0].name  = (char *)malloc(strlen(#the_species"_")+ \
                                                            strlen(Phase_Name[ip])+1); \
            sprintf(var.species[ispecies].phase[ip].info[0].name,"%s_%s",#the_species,Phase_Name[ip]); \
            var.species[ispecies].phase[ip].info[0].long_name = (char *)malloc(strlen(#the_species" ")+ \
                                                            strlen(Long_Phase_Name[ip])+1); \
            sprintf(var.species[ispecies].phase[ip].info[0].long_name,"%s %s",#the_species,Long_Phase_Name[ip]); \
            var.species[ispecies].phase[ip].info[0].units = (char *)malloc(strlen(#the_units)+1); \
            strcpy(var.species[ispecies].phase[ip].info[0].units,#the_units); \
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

void set_var_props(void) 
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
  SET_THERMO(HDRY_INDEX,hdry,air_dry_hybrid_density,dry-air hybrid density,kg/m^2/K);
  SET_THERMO(THETA_INDEX,theta,air_potential_temperature,potential temperature,K);
  SET_THERMO(FPARA_INDEX,fpara,mole_fraction_of_para_in_hydrogen,para fraction of molecular hydrogen,mol/mol);

  /* 
   * Optional species (Qs).
   */
  SET_SPECIES(H_2O_INDEX,H_2O,mixing_ratio_of_water,water mixing ratio,kg/kg);
  var.species[H_2O_INDEX].enthalpy_change        = enthalpy_change_H_2O;
  var.species[H_2O_INDEX].sat_vapor_p            = sat_vapor_p_H_2O;
  
  SET_SPECIES(NH_3_INDEX,NH_3,mixing_ratio_of_ammonia,ammonia mixing ratio,kg/kg);
  var.species[NH_3_INDEX].enthalpy_change        = enthalpy_change_NH_3;
  var.species[NH_3_INDEX].sat_vapor_p            = sat_vapor_p_NH_3;

  SET_SPECIES(H_2S_INDEX,H_2S,mixing_ratio_of_hydrogen_sulfide,hydrogen sulfide mixing ratio,kg/kg);
  var.species[H_2S_INDEX].enthalpy_change        = enthalpy_change_H_2S;
  var.species[H_2S_INDEX].sat_vapor_p            = sat_vapor_p_H_2S;
  
  SET_SPECIES(CH_4_INDEX,CH_4,mixing_ratio_of_methane,methane mixing ratio,kg/kg);
  var.species[CH_4_INDEX].enthalpy_change        = enthalpy_change_CH_4;
  var.species[CH_4_INDEX].sat_vapor_p            = sat_vapor_p_CH_4;

  SET_SPECIES(C_2H_2_INDEX,C_2H_2,mixing_ratio_of_acetylene,acetylyne mixing ratio,kg/kg);
  var.species[C_2H_2_INDEX].enthalpy_change      = enthalpy_change_C_2H_2;
  var.species[C_2H_2_INDEX].sat_vapor_p          = sat_vapor_p_C_2H_2;

  SET_SPECIES(C_2H_6_INDEX,C_2H_6,mixing_ratio_of_ethane,ethane mixing ratio,kg/kg);
  var.species[C_2H_6_INDEX].enthalpy_change      = enthalpy_change_C_2H_6;
  var.species[C_2H_6_INDEX].sat_vapor_p          = sat_vapor_p_C_2H_6;
  
  SET_SPECIES(CO_2_INDEX,CO_2,mixing_ratio_of_carbon_dioxide,carbon dioxide mixing ratio,kg/kg);
  var.species[CO_2_INDEX].enthalpy_change        = enthalpy_change_CO_2;
  var.species[CO_2_INDEX].sat_vapor_p            = sat_vapor_p_CO_2;
  
  SET_SPECIES(NH_4SH_INDEX,NH_4SH,mixing_ratio_of_ammonium_hydrosulfide,ammonium hydrosulfide mixing ratio,kg/kg);
  var.species[NH_4SH_INDEX].enthalpy_change      = enthalpy_change_NH_4SH;
  var.species[NH_4SH_INDEX].sat_vapor_p          = sat_vapor_p_NH_4SH;
  
  SET_SPECIES(O_3_INDEX,O_3,mixing_ratio_of_ozone,ozone mixing ratio,kg/kg);
  var.species[O_3_INDEX].enthalpy_change         = enthalpy_change_O_3;
  var.species[O_3_INDEX].sat_vapor_p             = sat_vapor_p_O_3;
  
  SET_SPECIES(N_2_INDEX,N_2,mixing_ratio_of_nitrogen,nitrogen mixing ratio,kg/kg);
  var.species[N_2_INDEX].enthalpy_change         = enthalpy_change_N_2;
  var.species[N_2_INDEX].sat_vapor_p             = sat_vapor_p_N_2;

  /*
   * Set turbulence-model variables.
   */
  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
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
  SET_DIAG(HDRY3_INDEX,hdry3,air_dry_hybrid_density,dry-air hybrid density (interface),kg/m^2/K);
  SET_DIAG(PDRY3_INDEX,pdry3,air_dry_pressure,dry-air pressure (interface),Pa);
  SET_DIAG(P2_INDEX,p2,air_pressure,total air pressure (layer),Pa);
  SET_DIAG(P3_INDEX,p3,air_pressure,total air pressure (interface),Pa);
  SET_DIAG(THETA2_INDEX,theta2,air_potential_temperature,potential temperature (layer),K);
  SET_DIAG(H2_INDEX,h2,air_hybrid_density,total hybrid density (layer),kg/m^2/K);
  SET_DIAG(H3_INDEX,h3,air_hybrid_density,total hybrid density (interface),kg/m^2/K);
  SET_DIAG(T2_INDEX,t2,air_temperature,temperature (layer),K);
  SET_DIAG(T3_INDEX,t3,air_temperature,temperature (interface),K);
  SET_DIAG(RHO2_INDEX,rho2,air_density,total density (layer),kg/m^3);
  SET_DIAG(RHO3_INDEX,rho3,air_density,total density (interface),kg/m^3);
  SET_DIAG(EXNER3_INDEX,exner3,exner_function,Exner function (interface),J/kg);
  SET_DIAG(FGIBB3_INDEX,fgibb3,gibbs_term_hydrogen,Gibbs term (ortho-para H2),J/kg);
  SET_DIAG(GZ3_INDEX,gz3,geopotential,geopotential (interface),m^2/s^2);
  SET_DIAG(MONT3_INDEX,mont3,montgomery_potential,Montgomery potential (interface),m^2/s^2);
  SET_DIAG(HEAT3_INDEX,heat3,heating_rate,heating rate per mass,W/kg);
  SET_DIAG(PV3_INDEX,pv3,ertel_potential_vorticity,potential vorticity,m^2/s K/kg);
  SET_DIAG(EDDY_PV3_INDEX,eddy_pv3,eddy_ertel_potential_vorticity,eddy potential vorticity,m^2/s K/kg);
  SET_DIAG(RI2_INDEX,ri2,local_richardson_number,local Richardson number (layer),s^2/s^2);
  SET_DIAG(VORT3_INDEX,vort3,atmosphere_relative_vorticity,relative vorticity (layer),1/s);
  SET_DIAG(DIV_UV3_INDEX,div_uv3,divergence_of_wind,horizontal divergence,1/s);
  SET_DIAG(W3_INDEX,w3,hybrid_upward_air_velocity,hybrid vertical velocity,K/s);
  SET_DIAG(DZDT3_INDEX,dzdt3,upward_air_velocity,vertical velocity,m/s);
  SET_DIAG(DIFFUSION_COEF_UV_INDEX,diffusion_coef_uv,wind_diffusion,wind eddy diffusion,m^2/s);
  SET_DIAG(DIFFUSION_COEF_THETA_INDEX,diffusion_coef_theta,potential_temperature_diffusion,potential temperature eddy diffusion,m^2/s);
  SET_DIAG(DIFFUSION_COEF_MASS_INDEX,diffusion_coef_mass,mass_diffusion,mass eddy diffusion,m^2/s);
  SET_DIAG(U_SPINUP_INDEX,u_spinup,rayleigh_drag_eastward_wind,Rayleigh-drag zonal wind profile,m/s);
  SET_DIAG(GZ_SURFACE_INDEX,gz_surface,surface_geopotential,surface geopotential,m^2/s^2);

  if (grid.nmt_physics_on == 1) {
    /*
     * nmt_physics diagnostic variables.
     */
    SET_DIAG(DRY_ENTROPY_INDEX,dry_entropy,air_dry_entropy,dry entropy,J/kg/K);
    SET_DIAG(MOIST_ENTROPY_INDEX,moist_entropy,air_moist_entropy,moist entropy,J/kg/K);
    SET_DIAG(SAT_MOIST_ENTROPY_INDEX,sat_moist_entropy,air_saturation_moist_entropy,saturation moist entropy,J/kg/K);
    SET_DIAG(THE_FLUX_INDEX,the_flux,the_flux,the flux,K/kg/m^2);
    SET_DIAG(RT_FLUX_INDEX,rt_flux,rt_flux,RT flux,g/g/kg/m^2);
    SET_DIAG(U_FLUX_INDEX,u_flux,eastward_wind_flux,zonal-wind flux,m/s/kg/m^2);
    SET_DIAG(V_FLUX_INDEX,v_flux,northward_wind_flux,meridional-wind flux,m/s/kg/m^2);
    SET_DIAG(CONVTHROT_INDEX,convthrot,convective_throttle,convective throttle,nodim);
    SET_DIAG(FLUXTHROT_INDEX,fluxthrot,flux_throttle,flux throttle,nodim);
    SET_DIAG(RAIN_RATE_INDEX,rain_rate,rainfall_flux,rainfall flux,kg/m^2/s);
  }

  /*
   * Set species molar_mass values.
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    var.species[is].molar_mass = molar_mass(is);
  }

  /*
   * These functions are defined in epic/src/shared/microphysics/epic_microphysics.c.
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
            free(var.species[ispecies].phase[ip].info[0].name); \
            free(var.species[ispecies].phase[ip].info[0].units); \
          }

#define FREE_DIAG(diag) \
          free(var.diag.info[0].name); \
          free(var.diag.info[0].standard_name); \
          free(var.diag.info[0].long_name); \
          free(var.diag.info[0].units);

/*
 */
void free_var_props(void)
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
  FREE_THERMO(hdry);
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
  FREE_SPECIES(C_2H_6_INDEX);
  FREE_SPECIES(CO_2_INDEX);
  FREE_SPECIES(NH_4SH_INDEX);
  FREE_SPECIES(O_3_INDEX);
  FREE_SPECIES(N_2_INDEX);

  /*
   * Turbulence-model variables.
   */
  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
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
  FREE_DIAG(hdry3);
  FREE_DIAG(pdry3);
  FREE_DIAG(p2);
  FREE_DIAG(p3);
  FREE_DIAG(theta2);
  FREE_DIAG(h2);
  FREE_DIAG(h3);
  FREE_DIAG(t2);
  FREE_DIAG(t3);
  FREE_DIAG(rho2);
  FREE_DIAG(rho3);
  FREE_DIAG(exner3);
  FREE_DIAG(fgibb3);
  FREE_DIAG(gz3);
  FREE_DIAG(mont3);
  FREE_DIAG(heat3);
  FREE_DIAG(pv3);
  FREE_DIAG(eddy_pv3);
  FREE_DIAG(ri2);
  FREE_DIAG(vort3);
  FREE_DIAG(div_uv3);
  FREE_DIAG(w3);
  FREE_DIAG(dzdt3);
  FREE_DIAG(diffusion_coef_uv);
  FREE_DIAG(diffusion_coef_theta);
  FREE_DIAG(diffusion_coef_mass);
  FREE_DIAG(u_spinup);
  FREE_DIAG(gz_surface);

  if (grid.nmt_physics_on == 1) {
    /*
     * nmt_physics diagnostic variables.
     */
    FREE_DIAG(gz_surface);
    FREE_DIAG(dry_entropy);
    FREE_DIAG(moist_entropy);
    FREE_DIAG(sat_moist_entropy);
    FREE_DIAG(the_flux);
    FREE_DIAG(rt_flux);
    FREE_DIAG(u_flux);
    FREE_DIAG(v_flux);
    FREE_DIAG(convthrot);
    FREE_DIAG(fluxthrot);
    FREE_DIAG(rain_rate);
  }

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
    is,ip,
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

  if (strcmp(grid.geometry,"globe") == 0) {
    if(para.nproc > grid.nj) {
      if (IAMNODE == 0) {
        fprintf(stderr,"This model must be run on <= %d processors\n",grid.nj);
      }
      exit(1);
    }
  }
#else
  grid.we_num_nodes = 1;

  /* Indicate whether the processor has poles in its range: */
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
   * Allocate memory for prognostic variables, tendencies, and moments.
   *
   * NOTE:  The on_list[] variable is a tri-state variable (-1, 0, 1), not a boolean (0, 1).
   *        Hence, it should not be tested like a boolean. The states are:
   *            1 => turned on
   *            0 => turned off
   *           -1 => not printed as an option
   */

  if (var.on_list[U_INDEX] == 1) {
    var.u.on       = TRUE;
    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      var.u.value    = fvector(0,  Nelem3d-1,dbmsname);
      var.u.tendency = fvector(0,3*Nelem3d-1,dbmsname);
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      var.u.value    = fvector(0,2*Nelem3d-1,dbmsname);
      var.u.tendency = fvector(0,  Nelem3d-1,dbmsname);
    }
    if (var.extract_on_list[U_INDEX] == 1) var.u.extract_on = TRUE;
  }

  if (var.on_list[V_INDEX] == 1) {
    var.v.on       = TRUE;
    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      var.v.value    = fvector(0,  Nelem3d-1,dbmsname);
      var.v.tendency = fvector(0,3*Nelem3d-1,dbmsname);
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      var.v.value    = fvector(0,2*Nelem3d-1,dbmsname);
      var.v.tendency = fvector(0,  Nelem3d-1,dbmsname);
    }
    if (var.extract_on_list[V_INDEX] == 1) var.v.extract_on = TRUE;
  }

  if (var.on_list[HDRY_INDEX] == 1) {
    var.hdry.on    = TRUE;
    var.hdry.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[HDRY_INDEX] == 1) var.hdry.extract_on = TRUE;
  }

  if (var.on_list[THETA_INDEX] == 1) {
    var.theta.on    = TRUE;
    var.theta.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[THETA_INDEX] == 1) var.theta.extract_on = TRUE;
  }

  if (var.on_list[FPARA_INDEX] == 1) {
    var.fpara.on    = TRUE;
    var.fpara.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[FPARA_INDEX] == 1) var.fpara.extract_on = TRUE;
  }

  /*
   * Turn on phases appropriate to choice of physics package.
   * The phases are switched on here, whether or not the species are invoked. 
   */
  var.species[C_2H_2_INDEX].phase[VAPOR ].on = TRUE;
  var.species[C_2H_2_INDEX].phase[ICE   ].on = FALSE;
  var.species[C_2H_2_INDEX].phase[LIQUID].on = FALSE;
  var.species[C_2H_2_INDEX].phase[RAIN  ].on = FALSE;
  var.species[C_2H_2_INDEX].phase[SNOW  ].on = FALSE;

  var.species[C_2H_6_INDEX].phase[VAPOR ].on = TRUE;
  var.species[C_2H_6_INDEX].phase[ICE   ].on = FALSE;
  var.species[C_2H_6_INDEX].phase[LIQUID].on = FALSE;
  var.species[C_2H_6_INDEX].phase[RAIN  ].on = FALSE;
  var.species[C_2H_6_INDEX].phase[SNOW  ].on = FALSE;

  if (grid.microphysics_on == TRUE) {
    /*
     * The Palotai & Dowling cloud microphysics scheme uses five phases for each species.
     * Currently, water and ammonia are implemented.
     */
    var.species[H_2O_INDEX].phase[VAPOR ].on = TRUE;
    var.species[H_2O_INDEX].phase[ICE   ].on = TRUE;
    var.species[H_2O_INDEX].phase[LIQUID].on = TRUE;
    var.species[H_2O_INDEX].phase[RAIN  ].on = TRUE;
    var.species[H_2O_INDEX].phase[SNOW  ].on = TRUE;

    var.species[NH_3_INDEX].phase[VAPOR ].on = TRUE;
    var.species[NH_3_INDEX].phase[ICE   ].on = TRUE;
    var.species[NH_3_INDEX].phase[LIQUID].on = TRUE;
    var.species[NH_3_INDEX].phase[RAIN  ].on = TRUE;
    var.species[NH_3_INDEX].phase[SNOW  ].on = TRUE;
  }
  else if (grid.nmt_physics_on == TRUE) {
    /*
     * The New Mexico Tech (NMT) physics scheme uses total advected water 
     * (vapor + cloud ice + cloud liquid) as a single, quasi-conserved 
     * prognostic variable. We use the VAPOR phase to hold this variable,
     * and turn the others off.
     */
    var.species[H_2O_INDEX].phase[VAPOR ].on = TRUE;
    var.species[H_2O_INDEX].phase[ICE   ].on = FALSE;
    var.species[H_2O_INDEX].phase[LIQUID].on = FALSE;
    var.species[H_2O_INDEX].phase[RAIN  ].on = FALSE;
    var.species[H_2O_INDEX].phase[SNOW  ].on = FALSE;
  }

  /*
   * NOTE: Do not use grid.nq yet, it is set below.
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.on_list[is] == 1) {
      var.species[is].on = TRUE;
      if (var.extract_on_list[is] == 1) var.species[is].extract_on = TRUE;
      for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
        /*
         * The appropriate phases for the chosen physics package should already be turned on.
         */
        if (var.species[is].phase[ip].on) {
          var.species[is].phase[ip].q = fvector(0,Nelem3d-1,dbmsname);
          if (var.extract_on_list[is] == 1) var.species[is].phase[ip].extract_on = TRUE;
        }
      }
    }
  }
  if (var.on_list[NU_TURB_INDEX] == 1) {
    var.nu_turb.on    = TRUE;
    var.nu_turb.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[NU_TURB_INDEX] == 1) var.nu_turb.extract_on = TRUE;
  }

  /*
   * Count total number of active species-phases, grid.nq.
   */
  grid.nq = 0;
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
        if (var.species[is].phase[ip].on) {
          grid.nq++;
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
             grid.is[itmp] = is;
             grid.ip[itmp] = ip;
             itmp++;
          }
        }
      }
    }
  }

  /*
   * Allocate memory for 1D arrays.
   */
  grid.lon        = fvector(0,2*(grid.ni+1),dbmsname);
  grid.lat        = fvector(0,2*(grid.nj+1),dbmsname); 
  grid.rln        = fvector(0,2*(grid.nj+1),dbmsname);
  grid.rlt        = fvector(0,2*(grid.nj+1),dbmsname);  
  grid.f          = fvector(0,2*(grid.nj+1),dbmsname);
  grid.m          = fvector(0,2*(grid.nj+1),dbmsname);
  grid.n          = fvector(0,2*(grid.nj+1),dbmsname);
  grid.mn         = fvector(0,2*(grid.nj+1),dbmsname);
  grid.g          = fvector(0,2*(grid.nj+1),dbmsname);
  grid.sigmatheta = dvector(0,2*(grid.nk+1),dbmsname);
  grid.dsgth      = fvector(0,2*(grid.nk+1),dbmsname);
  grid.dsgth_inv  = fvector(0,2*(grid.nk+1),dbmsname);
  grid.p_ref      = fvector(0,2*(grid.nk+1),dbmsname);
  grid.theta_ref  = fvector(0,2*(grid.nk+1),dbmsname);
  grid.h_min      = fvector(0,   grid.nk+1, dbmsname);

  if (grid.k_sponge > 0) {
    grid.t_sponge_inv = fvector(0,grid.k_sponge,dbmsname);
  }

  if (var.ntp > 0) {
    var.pdat  = fvector(0,var.ntp-1,dbmsname);
    var.tdat  = fvector(0,var.ntp-1,dbmsname);
    var.dtdat = fvector(0,var.ntp-1,dbmsname);
  }

  if (var.n_t_cool > 0) var.t_cool_table = ftriplet(0,var.n_t_cool-1,dbmsname);

  /*
   * Allocate memory for diagnostic arrays and parameter arrays.
   */
  var.gz_surface.on    = TRUE;
  var.gz_surface.value = fvector(0,NELEM2D-1,dbmsname);
  if (var.extract_on_list[GZ_SURFACE_INDEX] == 1) var.gz_surface.extract_on = TRUE;

  var.hdry3.on       = TRUE;
  var.hdry3.value    = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[HDRY3_INDEX] == 1) var.hdry3.extract_on = TRUE;

  var.pdry3.on       = TRUE;
  var.pdry3.value    = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[PDRY3_INDEX] == 1) var.pdry3.extract_on = TRUE;
 
  var.p2.on          = TRUE;
  var.p2.value       = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[P2_INDEX] == 1) var.p2.extract_on = TRUE;

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
  if (var.extract_on_list[P3_INDEX] == 1) var.p3.extract_on = TRUE;

  var.theta2.on      = TRUE;
  var.theta2.value   = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[THETA2_INDEX] == 1) var.theta2.extract_on = TRUE;

  var.h2.on          = TRUE;
  if (grid.nq > 0) {
    var.h2.value = fvector(0,Nelem3d-1,dbmsname);
  }
  else {
    /*
     * No difference between HDRY and H2, so assign them the same memory.
     */
    var.h2.value = var.hdry.value;
  }
  if (var.extract_on_list[H2_INDEX] == 1) var.h2.extract_on = TRUE;

  var.h3.on          = TRUE;
  if (grid.nq > 0) {
    var.h3.value = fvector(0,Nelem3d-1,dbmsname);
  }
  else {
    /*
     * No difference between HDRY3 and H3, so assign them the same memory.
     */
    var.h3.value = var.hdry3.value;
  }
  if (var.extract_on_list[H3_INDEX] == 1) var.h3.extract_on = TRUE;

  var.t2.on          = TRUE;
  var.t2.value       = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[T2_INDEX] == 1) var.t2.extract_on = TRUE;

  var.t3.on          = TRUE;
  var.t3.value       = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[T3_INDEX] == 1) var.t3.extract_on = TRUE;

  var.rho2.on        = TRUE;
  var.rho2.value     = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[RHO2_INDEX] == 1) var.rho2.extract_on = TRUE;

  var.rho3.on        = TRUE;
  var.rho3.value     = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[RHO3_INDEX] == 1) var.rho3.extract_on = TRUE;

  var.exner3.on      = TRUE;
  var.exner3.value   = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[EXNER3_INDEX] == 1) var.exner3.extract_on = TRUE;

  var.gz3.on         = TRUE;
  var.gz3.value      = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[GZ3_INDEX] == 1) var.gz3.extract_on = TRUE;

  var.mont3.on       = TRUE;
  var.mont3.value    = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[MONT3_INDEX] == 1) var.mont3.extract_on = TRUE;

  var.heat3.on       = TRUE;
  var.heat3.value    = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[HEAT3_INDEX] == 1) var.heat3.extract_on = TRUE;

  var.pv3.on         = TRUE;
  var.pv3.value      = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[PV3_INDEX] == 1) var.pv3.extract_on = TRUE;

  var.ri2.on         = TRUE;
  var.ri2.value      = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[RI2_INDEX] == 1) var.ri2.extract_on = TRUE;

  var.div_uv3.on     = TRUE;
  var.div_uv3.value  = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[DIV_UV3_INDEX] == 1) var.div_uv3.extract_on = TRUE;

  /* 
   * The following diagnostic variables are available for output, but are
   * not assigned permanent memory.
   */
  if (var.extract_on_list[VORT3_INDEX] == 1) {
    var.vort3.extract_on = TRUE;
  }
  if (var.extract_on_list[EDDY_PV3_INDEX] == 1) {
    var.eddy_pv3.extract_on = TRUE;
  }


  var.w3.on          = TRUE;
  var.w3.value       = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[W3_INDEX] == 1) var.w3.extract_on = TRUE;

  var.dzdt3.on       = TRUE;
  var.dzdt3.value    = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[DZDT3_INDEX] == 1) var.dzdt3.extract_on = TRUE;

  var.diffusion_coef_uv.on    = TRUE;
  var.diffusion_coef_uv.value = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[DIFFUSION_COEF_UV_INDEX] == 1) var.diffusion_coef_uv.extract_on = TRUE;

  var.diffusion_coef_theta.on    = TRUE;
  var.diffusion_coef_theta.value = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[DIFFUSION_COEF_THETA_INDEX] == 1) var.diffusion_coef_theta.extract_on = TRUE;

  var.diffusion_coef_mass.on     = TRUE;
  var.diffusion_coef_mass.value  = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[DIFFUSION_COEF_MASS_INDEX] == 1) var.diffusion_coef_mass.extract_on = TRUE;

  var.u_spinup.on    = TRUE;
  var.u_spinup.value = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[U_SPINUP_INDEX] == 1) var.u_spinup.extract_on = TRUE;

  var.fgibb3.on    = TRUE;
  var.fgibb3.value = fvector(0,Nelem3d-1,dbmsname);
  if (var.extract_on_list[FGIBB3_INDEX] == 1) var.fgibb3.extract_on = TRUE;

  if (grid.nmt_physics_on == 1) {
    /*
     * nmt_physics diagnostic variables
     */
    var.dry_entropy.on    = TRUE;
    var.dry_entropy.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[DRY_ENTROPY_INDEX] == 1) var.dry_entropy.extract_on = TRUE;

    var.moist_entropy.on    = TRUE;
    var.moist_entropy.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[MOIST_ENTROPY_INDEX] == 1) var.moist_entropy.extract_on = TRUE;

    var.sat_moist_entropy.on    = TRUE;
    var.sat_moist_entropy.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[SAT_MOIST_ENTROPY_INDEX] == 1) var.sat_moist_entropy.extract_on = TRUE;

    var.the_flux.on    = TRUE;
    var.the_flux.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[THE_FLUX_INDEX] == 1) var.the_flux.extract_on = TRUE;

    var.rt_flux.on    = TRUE;
    var.rt_flux.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[RT_FLUX_INDEX] == 1) var.rt_flux.extract_on = TRUE;

    var.u_flux.on    = TRUE;
    var.u_flux.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[U_FLUX_INDEX] == 1) var.u_flux.extract_on = TRUE;

    var.v_flux.on    = TRUE;
    var.v_flux.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[V_FLUX_INDEX] == 1) var.v_flux.extract_on = TRUE;

    var.convthrot.on    = TRUE;
    var.convthrot.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[CONVTHROT_INDEX] == 1) var.convthrot.extract_on = TRUE;

    var.fluxthrot.on    = TRUE;
    var.fluxthrot.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[FLUXTHROT_INDEX] == 1) var.fluxthrot.extract_on = TRUE;

    var.rain_rate.on    = TRUE;
    var.rain_rate.value = fvector(0,Nelem3d-1,dbmsname);
    if (var.extract_on_list[RAIN_RATE_INDEX] == 1) var.rain_rate.extract_on = TRUE;
  }

  /*
   * Allocate turbulence-model memory.
   */
  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
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
    iq;
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
  if (var.hdry.on) {
    free_fvector(var.hdry.value,0,Nelem3d-1,dbmsname);
  }
  if (var.theta.on) {
    free_fvector(var.theta.value,0,Nelem3d-1,dbmsname);
  }
  if (var.fpara.on) {
    free_fvector(var.fpara.value,0,Nelem3d-1,dbmsname);
  }
  for (iq = 0; iq < grid.nq; iq++) {
    free_fvector(var.species[grid.is[iq]].phase[grid.ip[iq]].q,0,Nelem3d-1,dbmsname);
  }
  if (var.nu_turb.on) {
    free_fvector(var.nu_turb.value,0,Nelem3d-1,dbmsname);
  }

  if (grid.nq > 0) {
    free_ivector(grid.is,0,grid.nq-1,dbmsname);
    free_ivector(grid.ip,0,grid.nq-1,dbmsname);
  }

  free_fvector(grid.lon,       0,2*(grid.ni+1),dbmsname);
  free_fvector(grid.lat,       0,2*(grid.nj+1),dbmsname); 
  free_fvector(grid.rln,       0,2*(grid.nj+1),dbmsname);
  free_fvector(grid.rlt,       0,2*(grid.nj+1),dbmsname);  
  free_fvector(grid.f,         0,2*(grid.nj+1),dbmsname);
  free_fvector(grid.m,         0,2*(grid.nj+1),dbmsname);
  free_fvector(grid.n,         0,2*(grid.nj+1),dbmsname);
  free_fvector(grid.mn,        0,2*(grid.nj+1),dbmsname);
  free_fvector(grid.g,         0,2*(grid.nj+1),dbmsname);
  free_dvector(grid.sigmatheta,0,2*(grid.nk+1),dbmsname);
  free_fvector(grid.dsgth,     0,2*(grid.nk+1),dbmsname);
  free_fvector(grid.dsgth_inv, 0,2*(grid.nk+1),dbmsname);
  free_fvector(grid.p_ref,     0,2*(grid.nk+1),dbmsname);
  free_fvector(grid.theta_ref, 0,2*(grid.nk+1),dbmsname);
  free_fvector(grid.h_min,     0,   grid.nk+1, dbmsname);

  if (grid.k_sponge > 0) {
    free_fvector(grid.t_sponge_inv,0,grid.k_sponge,dbmsname);
  }

  if (var.ntp > 0) {
    free_fvector(var.pdat, 0,var.ntp-1,dbmsname);
    free_fvector(var.tdat, 0,var.ntp-1,dbmsname);
    free_fvector(var.dtdat,0,var.ntp-1,dbmsname);
  }

  if (var.n_t_cool > 0) free_ftriplet(var.t_cool_table,0,var.n_t_cool-1,dbmsname);

  free_fvector(var.gz_surface.value,0,NELEM2D-1,dbmsname);

  free_fvector(var.hdry3.value,               0,Nelem3d-1,dbmsname);
  free_fvector(var.pdry3.value,               0,Nelem3d-1,dbmsname);
  free_fvector(var.p2.value,                  0,Nelem3d-1,dbmsname);
  free_fvector(var.theta2.value,              0,Nelem3d-1,dbmsname);
  if (grid.nq > 0) {
    free_fvector(var.h2.value,                0,Nelem3d-1,dbmsname);
    free_fvector(var.h3.value,                0,Nelem3d-1,dbmsname);
    free_fvector(var.p3.value,                0,Nelem3d-1,dbmsname);
  }
  free_fvector(var.t2.value,                  0,Nelem3d-1,dbmsname);
  free_fvector(var.t3.value,                  0,Nelem3d-1,dbmsname);
  free_fvector(var.rho2.value,                0,Nelem3d-1,dbmsname);
  free_fvector(var.rho3.value,                0,Nelem3d-1,dbmsname);
  free_fvector(var.exner3.value,              0,Nelem3d-1,dbmsname);
  free_fvector(var.gz3.value,                 0,Nelem3d-1,dbmsname);
  free_fvector(var.mont3.value,               0,Nelem3d-1,dbmsname);
  free_fvector(var.heat3.value,               0,Nelem3d-1,dbmsname);
  free_fvector(var.pv3.value,                 0,Nelem3d-1,dbmsname);
  free_fvector(var.ri2.value,                 0,Nelem3d-1,dbmsname);
  free_fvector(var.div_uv3.value,             0,Nelem3d-1,dbmsname);
  free_fvector(var.w3.value,                  0,Nelem3d-1,dbmsname);
  free_fvector(var.dzdt3.value,               0,Nelem3d-1,dbmsname);
  free_fvector(var.diffusion_coef_uv.value,   0,Nelem3d-1,dbmsname);
  free_fvector(var.diffusion_coef_theta.value,0,Nelem3d-1,dbmsname);
  free_fvector(var.diffusion_coef_mass.value, 0,Nelem3d-1,dbmsname);
  free_fvector(var.u_spinup.value,0,Nelem3d-1,dbmsname);
  if (var.fpara.on) {
    free_fvector(var.fgibb3.value,0,Nelem3d-1,dbmsname);
  }

  if (grid.nmt_physics_on == 1) {
    /*
     * nmt_physics diagnostic variables.
     */
    free_fvector(var.dry_entropy.value,         0,Nelem3d-1,dbmsname);
    free_fvector(var.moist_entropy.value,       0,Nelem3d-1,dbmsname);
    free_fvector(var.sat_moist_entropy.value,   0,Nelem3d-1,dbmsname);
    free_fvector(var.the_flux.value,            0,Nelem3d-1,dbmsname);
    free_fvector(var.rt_flux.value,             0,Nelem3d-1,dbmsname);
    free_fvector(var.u_flux.value,              0,Nelem3d-1,dbmsname);
    free_fvector(var.v_flux.value,              0,Nelem3d-1,dbmsname);
    free_fvector(var.convthrot.value,           0,Nelem3d-1,dbmsname);
    free_fvector(var.fluxthrot.value,           0,Nelem3d-1,dbmsname);
    free_fvector(var.rain_rate.value,           0,Nelem3d-1,dbmsname);
  }

  /*
   * Free turbulence-model memory.
   */
  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
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
 * Hybrid coordinate, sigmatheta, as a function of pressure and potential temperature.
 * Use double precision to increase accuracy of diagnostic theta calculations.
 */

double return_sigmatheta(register double theta,
                                register double p,
                                register double pbot,
                                register double ptop)
{
  register double
    sigma,
    sigmatheta;

  if (p <= ptop) {
    sigmatheta = theta;
  }
  else if (p <= pbot) {
    sigma      = get_sigma(pbot,p,ptop);
    sigmatheta = f_sigma(sigma)+g_sigma(sigma)*theta;
  }
  else {
    sprintf(Message,"pbot=%g, ptop=%g; p=%g out of range",pbot,ptop,p);
    epic_error("return_sigmatheta",Message);
  }

  return sigmatheta;
}

/*======================= end of return_sigmatheta() ========================*/

/*======================= f_sigma() =========================================*/
/*
 * Part of hybrid vertical coordinate definition,
 * sigmatheta = f(sigma)+g(sigma)*theta.
 * Use double precision to increase accuracy of diagnostic theta calculations.
 */

double f_sigma(double sigma)
{
  if (sigma < 0.) {
    return grid.zeta0;
  }
  else if (sigma <= grid.sigma_sigma) {
    /*
     * NOTE: If this formula is changed, then the assignment of grid.sigma_sigma
     *       needs a corresponding change.
     */
    return grid.zeta0+sigma*(grid.zeta1-grid.zeta0);
  }
  else if (sigma < 1.) {
    return (1.-g_sigma(sigma))*(grid.zeta0+sigma*(grid.zeta1-grid.zeta0));
  }
  else {
    return 0.;
  }
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
    coeff = 1./(1-exp(-grid.hybrid_alpha*(1.-grid.sigma_sigma)));

  if (sigma <= grid.sigma_sigma) {
    return 0.;
  }
  else if (sigma >= 1.) {
    return 1.;
  }
  else {
    return coeff*(1.-exp(-grid.hybrid_alpha*(sigma-grid.sigma_sigma)));
  }
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
      grid.lat[0] = -180.-grid.dlt;
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
 *  Compute Coriolis parameter, f, and geometric map factors, m, n, etc.
 */
void set_fmn(planetspec *planet)
{
  int
    jj,
    nj,ni;
  EPIC_FLOAT
    omega,re,rp,
    dlnr,dltr,lat,
    rln,rlt,
    dx,dy;
  EPIC_FLOAT
    lat0,m0,n0;
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

  if (strcmp(grid.geometry,"globe") == 0) {
    dlnr  = grid.dln*DEG;
    dltr  = grid.dlt*DEG;
    re    = planet->re;
    rp    = planet->rp;
    for (jj = 1; jj < 2*(nj+1); jj++) {
      lat  = DEG*(grid.lat)[jj];
      rln  = re/sqrt( 1.+pow(rp/re*tan(lat),2.) );
      rlt  = rln/( cos(lat)*( pow(sin(lat),2.)+
             pow(re/rp*cos(lat),2.)));

      (grid.rln)[jj] = rln;
      (grid.rlt)[jj] = rlt;
      (grid.f  )[jj] = 2.*omega*sin(lat);
      (grid.m  )[jj] = 1./(rln*dlnr);
      (grid.n  )[jj] = 1./(rlt*dltr);
      (grid.mn )[jj] = (grid.m)[jj]*(grid.n)[jj]; 
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
      (grid.rln)[2*0  ]  = 0.;
      (grid.rlt)[2*0  ]  = re*re/rp;
      (grid.f  )[2*0  ]  = 2.*omega*sin(lat);
      (grid.n  )[2*0+1] /= 1.+sqrt(.5);  
      /* wedge shaped area: */
      (grid.mn )[2*0+1]  = grid.n[2*0+1]*grid.n[2*0+1]/(.5*dlnr);
      (grid.mn )[2*0  ]  = 2.*(grid.mn)[2*0+1];
    }
    else {
      jj = 0;
      lat  = DEG*(grid.lat)[jj];
      rln  = re/sqrt( 1.+ pow(rp/re*tan(lat),2.) );
      rlt  = rln/( cos(lat)*(pow(sin(lat),2.)+
                   pow(re/rp*cos(lat),2.)) );

      (grid.rln)[jj] = rln;
      (grid.rlt)[jj] = rlt;
      (grid.f  )[jj] = 2.*omega*sin(lat);
      (grid.m  )[jj] = 1./(rln*dlnr);
      (grid.n  )[jj] = 1./(rlt*dltr);
      (grid.mn )[jj] = (grid.m)[jj]*(grid.n)[jj]; 
    }
    if (grid.globe_lattop == 90.) {
      /* north pole */
      lat = DEG*(grid.lat)[2*(nj+1)];
      (grid.rln)[2*(nj+1)  ]  = 0.;
      (grid.rlt)[2*(nj+1)  ]  = re*re/rp;
      (grid.f  )[2*(nj+1)  ]  = 2.*omega*sin(lat);
      (grid.n  )[2*(nj+1)-1] /= 1.+sqrt(.5);
      /* wedge shaped area: */
      (grid.mn )[2*(nj+1)-1]  = grid.n[2*(nj+1)-1]*grid.n[2*(nj+1)-1]/(.5*dlnr);
      (grid.mn )[2*(nj+1)  ]  = 2.*(grid.mn)[2*(nj+1)-1];
    }
    else {
      jj = 2*(nj+1);
      lat  = DEG*(grid.lat)[jj];
      rln  = re/sqrt( 1.+ pow(rp/re*tan(lat),(EPIC_FLOAT)2.) );
      rlt  = rln/( cos(lat)*( pow(sin(lat),(EPIC_FLOAT)2.)+
             pow(re/rp*cos(lat),(EPIC_FLOAT)2.)));

      (grid.rln)[jj] = rln;
      (grid.rlt)[jj] = rlt;
      (grid.f  )[jj] = 2.*omega*sin(lat);
      (grid.m  )[jj] = 1./(rln*dlnr);
      (grid.n  )[jj] = 1./(rlt*dltr);
      (grid.mn )[jj] = grid.m[jj]*grid.n[jj]; 
    }
  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      lat = DEG*(grid.f_plane_lat0);
      dx  = 2.*(grid.f_plane_half_width)/ni;
      dy  = dx;

      for (jj = 0; jj <= 2*(nj+1); jj++) {
        (grid.rln)[jj] = 1.;
        (grid.rlt)[jj] = 1.;
        (grid.f  )[jj] = 2.*omega*sin(lat);
        (grid.m  )[jj] = 1./dx;
        (grid.n  )[jj] = 1./dy;
        (grid.mn )[jj] = 1./(dx*dy);
      }
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      lat  = DEG*(grid.f_plane_lat0);
      dlnr = grid.dln*DEG;
      dy   = grid.f_plane_half_width*(grid.dlt/90.);
      rln  = grid.f_plane_half_width+dy;

      for (jj = 0; jj < 2*(nj+1); jj++) {
        dx = rln*dlnr;
        (grid.rln)[jj] = rln;
        (grid.rlt)[jj] = 1.;
        (grid.f  )[jj] = 2.*omega*sin(lat);
        (grid.m  )[jj] = 1./dx;
        (grid.n  )[jj] = 1./dy;
        (grid.mn )[jj] = grid.m[jj]*grid.n[jj];
        rln -= dy/2.;
      }
      /* pole */
      (grid.rln)[2*(nj+1)  ]  = 0.;
      (grid.rlt)[2*(nj+1)  ]  = 1.;
      (grid.f  )[2*(nj+1)  ]  = 2.*omega*sin(lat);
      (grid.n  )[2*(nj+1)-1] /= 1.+sqrt(.5);
      /* wedge shaped area: */
      (grid.mn )[2*(nj+1)-1]  = grid.n[2*(nj+1)-1]*grid.n[2*(nj+1)-1]/(.5*dlnr);
      (grid.mn )[2*(nj+1)  ]  = 2.*(grid.mn)[2*(nj+1)-1];
    }
  }

  /*
   * Set grid.dy0 = dy at LAT0.
   */
  if (strcmp(grid.geometry,"globe") == 0) {
    lat0 = LAT0*DEG;  
    rln  = planet->re/sqrt( 1.+ pow(planet->rp/planet->re*tan(lat0),2.) );
    rlt  = rln/( cos(lat0)*( pow(sin(lat0),2.) +
                pow(planet->re/planet->rp*cos(lat0),2.) ) );
    grid.dy0 =rlt*grid.dlt*DEG;
  }
  else if (strcmp(grid.geometry,"f-plane")  == 0) { 
    if (strcmp(grid.f_plane_map,"polar") == 0) {
      grid.dy0 = grid.f_plane_half_width/grid.nj;
    }
    else {
      grid.dy0 = grid.f_plane_half_width/grid.ni;
    }
  }
  else {
    sprintf(Message,"unrecognized grid.geometry %s",grid.geometry);
    epic_error(dbmsname,Message);
  }

  return;
}

/*======================= end of set_fmn() ==================================*/

/*======================= set_gravity() =====================================*/

/*
 * Calculate gravity, g [m/s^2], as a function of planetographic latitude
 * on the reference surface.
 *
 * See eqn (41) of Yoder C, 1995, Global Earth Physics: A handbook of physical constants,
 *   http://www.agu.org/reference/gephys/4_yoder.pdf
 */

void set_gravity(planetspec *planet)
{
  register int
    jj;
  register double
    a,b,ge,gp,
    spin_factor,
    sinlat2,coslat2;

  a = planet->re;
  b = planet->rp;

  spin_factor = (planet->omega_sidereal*a)*(planet->omega_sidereal*a)*(a/planet->GM);
  ge          = (1.+1.5*planet->J2-spin_factor)*planet->GM/(a*a);
  gp          = (1.-3.*(a/b)*(a/b)*planet->J2)*planet->GM/(b*b);

  for (jj = 0; jj <= 2*(grid.nj+1); jj++) {
    sinlat2  = sin(grid.lat[jj]*DEG);
    sinlat2 *= sinlat2;
    coslat2  = cos(grid.lat[jj]*DEG);
    coslat2 *= coslat2;

    grid.g[jj] = (a*ge*coslat2+b*gp*sinlat2)/sqrt(a*a*coslat2+b*b*sinlat2);
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
      sprintf(Message,"grid.dsgth[kk=%d]=%g\n",kk+1,grid.dsgth[kk+1]);
      epic_error(dbmsname,Message);
    }
  }

  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;
    grid.dsgth[kk] = grid.sigmatheta[kk-1]-grid.sigmatheta[kk+1];

    if (grid.dsgth[kk] > 0.) {
      grid.dsgth_inv[kk] = 1./grid.dsgth[kk];
    }
    else{
      sprintf(Message,"grid.dsgth[kk=%d]=%g",kk,grid.dsgth[kk]);
      epic_error(dbmsname,Message);
    }
  }

  return;
}

/*======================= end of set_dsgth() ================================*/

/*======================= set_sponge() ======================================*/

void set_sponge(void)
{
  int
    K;
  const EPIC_FLOAT
    /* 
     * Minimum relaxation time (strongest effect) in units of timesteps.
     *
     * NOTE: Use >= 5.0 to avoid numerical instability.
     *       Actually, even 5.0 can be too strong in practice.
     *       A value of 30.0 seems to behave better.
     */
    min_relax_time = 30.;
  EPIC_FLOAT
    t0_inv,
    tmp;
  FILE
    *outfile;

  outfile = fopen("sponge_params.dat","w");
  fprintf(outfile,"Parameters related to the EPIC model's sponge in the top layers.\n");
  fprintf(outfile,"This sponge uses Rayleigh drag in uv_drag() to relax u to its initial value and v to zero. \n");
  fprintf(outfile," Layer Stiffness [1/s] \n");

  if (grid.k_sponge > 0) {
    t0_inv = 1./(min_relax_time*(EPIC_FLOAT)grid.dt);
    for (K = 0; K <= grid.k_sponge; K++) {
      tmp                  = (EPIC_FLOAT)(grid.k_sponge-K)/(grid.k_sponge);
      grid.t_sponge_inv[K] = t0_inv*.5*(1.-cos(M_PI*tmp));
      fprintf(outfile,"  %2d     %10.3e\n",K,grid.t_sponge_inv[K]);
    }
  }
  fclose(outfile);

  return;
}

/*======================= end of set_sponge() ===============================*/

/*======================= get_sigma() =======================================*/
/*
 * Use a function to calculate sigma so that its definition can be 
 * easily changed.
 *
 * NOTE: If this function is changed, then it is necessary to also modify
 *       its inverse function, get_p_sigma().
 *
 * The value of sigma should range monotonically from 0 at the bottom of layer
 * K = grid.nk to 1 at the bottom of layer K = 0.
 *
 * NOTE: We use log p instead of p because sigma defined with the latter stays
 *       close to 0 too long with respect to height when pbot is large, such as 
 *       for Venus or deep Jupiter models.
 */

double get_sigma(double pbot,
                        double p,
                        double ptop)
{
  return log(p/pbot)/log(ptop/pbot);
}

/*======================= end of get_sigma() ================================*/

/*======================= get_p_sigma() =====================================*/
/*
 * Inverse of get_sigma() function.  If one is changed, the other should
 * be matched accordingly.
 */
double get_p_sigma(double pbot,
                          double sigma,
                          double ptop)
{
  return pbot*exp(sigma*log(ptop/pbot));
}

/*======================= end of get_p_sigma() ==============================*/

/*======================= get_h() ===========================================*/

/*
 * Calculates h = -1/g dp/dsgth.
 */

EPIC_FLOAT get_h(planetspec *planet,
                 int         kk,
                 int         J,
                 int         I,
                 int         type)
{
  register int
    K,
    iq;
  register double
    h,tmp,sum;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_h";

#if EPIC_CHECK == 1
  /* Check validity of kk. */
  if (kk < 1 || kk > 2*grid.nk+2) {
    sprintf(Message,"kk = %d",kk);
    epic_error(dbmsname,Message);
  }
#endif

  if (kk%2 == 0) {
    /* 
     * Layer value for h. 
     */
    K = kk/2;
    if (K > grid.nk) {
      sprintf(Message,"K > grid.nk");
      epic_error(dbmsname,Message);
    }
    else {
      h = ((double)P3(K,J,I)-(double)P3(K-1,J,I))/((double)grid.g[2*J+1]*(double)grid.dsgth[kk]);
    }
  }
  else {
    /* 
     * Interface value for h. 
     * Konor and Arakawa favor averaging instead of 
     *     h = g_inv*(P2(K+1,J,I)-P2(K,J,I))*grid.dsgth_inv[kk];
     */
    K = (kk-1)/2;
    if (K == 0) {
      /* KA97 (B.14) */
      h = get_h(planet,kk+1,J,I,TOTAL);
    }
    else if (K < KHI) {
      h = onto_kk(planet,H3_INDEX,get_h(planet,kk-1,J,I,TOTAL),
                                  get_h(planet,kk+1,J,I,TOTAL),kk,J,I);
    }
    else if (K == KHI) {
      /* KA97 (B.14) */
      h = get_h(planet,kk-1,J,I,TOTAL);
    }
  }

  if (type == DRY) {
    /*
     * hdry = h/(1.+sum_i (Q_i))
     */
    sum = 1.;
    for (iq = 0; iq < grid.nq; iq++) {
      sum += get_var(planet,grid.is[iq],grid.ip[iq],grid.it_h,kk,J,I);
    }
    h /= sum;
  }

  return (EPIC_FLOAT)h;
}

/*======================= end of get_h() =====================================*/

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
    K,is,ip;
  register EPIC_FLOAT
    mu_dry_inv,
    pressure,
    partial_pressure,
    mole_fraction;
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

  if (kk > 2*grid.nk+1) {
    sprintf(Message,"kk > 2*grid.nk+1");
    epic_error(dbmsname,Message);
  }

  mu_dry_inv = planet->rgas/R_GAS;

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
      mole_fraction = mu_dry_inv;
      for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
        if (var.species[is].on) {
          mole_fraction += .5*(Q(is,VAPOR,K,J,I)+Q(is,VAPOR,K-1,J,I))/var.species[is].molar_mass;
        }
      }
      mole_fraction    = .5*(Q(index,VAPOR,K,J,I)+Q(index,VAPOR,K-1,J,I))/
                         (var.species[index].molar_mass*mole_fraction);
      partial_pressure = pressure*mole_fraction;  
      return partial_pressure;
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
      mole_fraction = mu_dry_inv;
      for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
        if (var.species[is].on) {
          mole_fraction += Q(is,VAPOR,K,J,I)/var.species[is].molar_mass;
        }
      }
      mole_fraction    = Q(index,VAPOR,K,J,I)/(var.species[index].molar_mass*mole_fraction);
      partial_pressure = pressure*mole_fraction;  
      return partial_pressure;
    }
    else {
      sprintf(Message,"index = %d unknown",index);
      epic_error(dbmsname,Message);
    }
  }
}

/*======================= end of get_p() =====================================*/

/*======================= alt_get_p() ========================================*/

/*
 * Return pressure consistent with hybrid vertical coordinate level 
 * and input temperature. The input integer kk is the model's doubled-K index
 * that facilitates reference to half levels.
 *
 * NOTE: This is different than return_press(), which is purely a 
 *       thermodynamical function, and is more closely related
 *       to get_p().
 */

/*
 * MAX_IT stands for "maximum number of iterations."
 */
#undef  MAX_IT
#define MAX_IT 10

#undef  MAX_PRESSURE_HYBRID_REGION   
#define MAX_PRESSURE_HYBRID_REGION  0.9999999*get_p_sigma(P3(KHI,J,I), grid.sigma_sigma, P3(0,J,I)) 

EPIC_FLOAT alt_get_p(planetspec *planet,
                     int         kk,
                     int         J,
                     int         I,
                     EPIC_FLOAT  temperature) 
{
  int
    it,
    error_flag;
  EPIC_FLOAT
    theta,
    pressure,
    p1,p2,ptol;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="alt_get_p";

  if (kk >= 2*grid.k_sigma-1) {
    /*
     * Pressure does not depend on temperature in the bottom
     * portion of the model.
     *
     * NOTE: P2_INDEX and P3_INDEX have the same effect in get_p().
     */
    return get_p(planet,P2_INDEX,kk,J,I);
  }
  else {
    /*
     * Iterate to get pressure that satisfies
     *   temperature-return_temp(fpara,pressure,theta(sgth,pressure))
     */
    if (var.fpara.on) {
      TMT_fp = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,kk,J,I);
    }
    else {
      TMT_fp = 0.25;
    }
    TMT_planet      = planet;
    TMT_kk          = kk;
    TMT_J           = J;
    TMT_I           = I;
    TMT_temperature = temperature;

    p1         = .5*get_p(planet,P2_INDEX,kk,J,I);
    p2         = 4.*p1;

    p1         = MAX( p1, 0.0 );
    p2         = MIN( p2, MAX_PRESSURE_HYBRID_REGION );

    ptol       = pow(machine_epsilon(),2./3.);

    for (it = 0; it < MAX_IT; it++) {
      error_flag = find_root(p1,p2,ptol,&pressure,t_minus_t_p);
      if (error_flag == 0) {
        /* Convergence */
        return pressure;
      }
      /* Try a wider interval. */
      p1 *= 0.5;
      p2 *= 2.0;

      p2  = MIN( p2, MAX_PRESSURE_HYBRID_REGION );
    }

    sprintf(Message,"exceeded MAX_IT = %d,  T0=%e\n",MAX_IT,temperature);
    epic_error(dbmsname,Message);
  }
}

/*======================= end of alt_get_p() =================================*/

/*======================= t_minus_t_p() ======================================*/

/*
 * For use with find_root().
 */

EPIC_FLOAT t_minus_t_p(EPIC_FLOAT p) {
  EPIC_FLOAT
    theta,
    sigma,
    ans;

  sigma = get_sigma(P3(grid.nk,TMT_J,TMT_I),p,P3(0,TMT_J,TMT_I));
  theta = (grid.sigmatheta[TMT_kk]-f_sigma(sigma))/g_sigma(sigma);
  ans   = TMT_temperature-return_temp(TMT_planet,TMT_fp,p,theta);

  return ans;
}

/*======================= end of t_minus_t_p() ===============================*/

/*======================= state_from_exner() =================================*/

/*
 * Determine thermodynamic state given vertical coordinate and exner = cp*temperature/theta.
 * Similar to alt_get_p().
 */

/*
 * MAX_IT stands for "maximum number of iterations."
 */
#undef  MAX_IT
#define MAX_IT 10

#undef  MAX_PRESSURE_HYBRID_REGION   
#define MAX_PRESSURE_HYBRID_REGION  0.9999999*get_p_sigma(P3(KHI,J,I), grid.sigma_sigma, P3(0,J,I)) 

void state_from_exner(planetspec *planet,
                      int         kk,
                      int         J,
                      int         I,
                      EPIC_FLOAT  exner,
                      EPIC_FLOAT *pressure,
                      EPIC_FLOAT *theta)
{
  int
    it,
    error_flag;
  EPIC_FLOAT
    p,
    p1,p2,ptol;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="state_from_exner";

  if (kk >= 2*grid.k_sigma-1) {
    sprintf(Message,"not implemented in sigma-coordinate region of model");
    epic_error(dbmsname,Message);
  }
  else {
    /*
     * Iterate to get pressure that satisfies
     *   exner-planet->cp*return_temp(fpara,pressure,theta(sgth,pressure))/theta = 0.
     */
    if (var.fpara.on) {
      EME_fp = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,kk,J,I);
    }
    else {
      EME_fp = 0.25;
    }
    EME_planet = planet;
    EME_kk     = kk;
    EME_J      = J;
    EME_I      = I;
    EME_exner  = exner;

    p1         = .5*get_p(planet,P2_INDEX,kk,J,I);
    p2         = 4.*p1;

    p1         = MAX( p1, 0.0 );
    p2         = MIN( p2, MAX_PRESSURE_HYBRID_REGION );

    ptol       = pow(machine_epsilon(),2./3.);

    for (it = 0; it < MAX_IT; it++) {
      error_flag = find_root(p1,p2,ptol,&p,exner_minus_exner_p);
      if (error_flag == 0) {
        /* Convergence */
        *pressure = p;
        *theta    = EME_theta;
        return;
      }
      /* Try a wider interval. */
      p1 *= 0.5;
      p2 *= 2.0;

      p2  = MIN( p2, MAX_PRESSURE_HYBRID_REGION );
    }

    sprintf(Message,"exceeded MAX_IT = %d,  exner=%e\n",MAX_IT,exner);
    epic_error(dbmsname,Message);
  }
}

/*======================= end of state_from_exner() ==========================*/

/*======================= exner_minus_exner_p() ==============================*/

/*
 * For use with find_root().
 */

EPIC_FLOAT exner_minus_exner_p(EPIC_FLOAT p) {
  EPIC_FLOAT
    temperature,
    sigma,
    ans;

  sigma       = get_sigma(P3(grid.nk,EME_J,EME_I),p,P3(0,EME_J,EME_I));
  EME_theta   = (grid.sigmatheta[EME_kk]-f_sigma(sigma))/g_sigma(sigma);
  temperature = return_temp(EME_planet,EME_fp,p,EME_theta);
  ans         = EME_exner-EME_planet->cp*temperature/EME_theta;  

  return ans;
}

/*======================= end of exner_minus_exner_p() =======================*/

/*======================= molar_mixing_ratio() ===============================*/

/*
 * The molar mass for dry air is R_GAS/planet->rgas.
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
   * Return 0. if species or phase is not on.
   */
  if (!var.species[is].on) {
    return 0.;
  }
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

#if EPIC_CHECK == 1
  if (species_index > MAX_NUM_PROGS-1) {
    /*
     * Error if called for a diagnostic variable.
     */
    sprintf(Message,"called with index=%d",species_index);
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
#endif

  if (kk%2 == 0) {
    /*
     * Layer value.
     */
    K = kk/2;

    /*
     * Handle pass-through cases.
     */
    switch (species_index) {
      case HDRY_INDEX:
        return HDRY(K,J,I);
      break;
      default:
        return onto_kk(planet,species_index,
                       get_var(planet,species_index,phase_index,IT,kk-1,J,I),
                       get_var(planet,species_index,phase_index,IT,kk+1,J,I),kk,J,I);
      break;
    } /* end of switch */
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
      case U_INDEX:
        return U(IT,K,J,I);
      break;
      case V_INDEX:
        return V(IT,K,J,I);
      break;
      case THETA_INDEX:
      case THETA2_INDEX:
        return THETA(K,J,I);
      break;
      case FPARA_INDEX:
        if (var.fpara.on) {
          return FPARA(K,J,I);
        }
        else {
          return 0.25;
        }
      break;
      case NU_TURB_INDEX:
        return NU_TURB(K,J,I);
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
    } /* end of switch */

    if (K == 0) {
      /*
       * Handle top of top layer.
       * Use linear extrapolation with respect to sigmatheta for HDRY.
       */
      switch (species_index) {
        case HDRY_INDEX:
          return HDRY(1,J,I)+(HDRY(1,J,I)-HDRY(2,J,I))*(grid.sigmatheta[1]-grid.sigmatheta[2])*
                          grid.dsgth_inv[3];
        break;
        default:
          sprintf(Message,"kk=%d, species_index=%d not recognized",kk,species_index);
          epic_error(dbmsname,Message);
        break;
      }  /* end of switch */
    }
    else if (K == KHI) {
      /*
       * Handle bottom of bottom layer.
       * Use linear extrapolation with respect to sigmatheta for HDRY.
       */
      switch (species_index) {
        case HDRY_INDEX:
          return HDRY(K,J,I)+(HDRY(K,J,I)-HDRY(K-1,J,I))*(grid.sigmatheta[kk]-grid.sigmatheta[kk-1])*
                          grid.dsgth_inv[kk-2];
        break;
        default:
          sprintf(Message,"kk=%d; species_index=%d not recognized",kk,species_index);
          epic_error(dbmsname,Message);
        break;
      }  /* end of switch */
    }
    else {
      /*
       * Handle interior interfaces.
       */
      return onto_kk(planet,species_index,
                     get_var(planet,species_index,phase_index,IT,kk-1,J,I),
                     get_var(planet,species_index,phase_index,IT,kk+1,J,I),kk,J,I);
    }
  }
}

/*======================= end of get_var() ===================================*/

/*======================= get_var_mean2d() ===================================*/

/*
 * Return the global mean of a 2 dimensional layer variable.
 *
 * Aaron Herrnstein, January 2006.
 */

EPIC_FLOAT get_var_mean2d(EPIC_FLOAT *a,
                          int         index)
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
      da_u[J-JLO] = 1./grid.mn[2*J+1];
      sum_da_u   += da_u[J-JLO];
    }
    sum_da_u *= grid.ni;
#if defined(EPIC_MPI)
    mpi_tmp = sum_da_u;
    MPI_Allreduce(&mpi_tmp, &sum_da_u, 1, float_type, MPI_SUM, para.comm);
#endif

    sum_da_v = 0.;
    for (J = JFIRST; J <= JHI; J++) { 
      da_v[J-JLO] = 1./grid.mn[2*J];
      sum_da_v   += da_v[J-JLO]; 
    }
    if (JLO == grid.jlo) {
      J = JLO;
      da_v[J-JLO] = 1./grid.mn[2*J];
      sum_da_v   += da_v[J-JLO];
    }
    if (JHI == grid.nj) {
      J = grid.nj+1;
      da_v[J-JLO] = 1./grid.mn[2*J];
      sum_da_v   += da_v[J-JLO];
    }
    sum_da_v *= grid.ni;
#if defined(EPIC_MPI)
    mpi_tmp = sum_da_v;
    MPI_Allreduce(&mpi_tmp, &sum_da_v, 1, float_type, MPI_SUM, para.comm);
#endif

    initialized = TRUE;
  }

  if (index == V_INDEX || index == PV3_INDEX) {
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

  if (index == PV3_INDEX) {
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
  MPI_Allreduce(&mpi_tmp, &var_mean, 1, float_type, MPI_SUM, para.comm);
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
 * vertical index.  We find that having all the vertical averaging schemes
 * together in one place makes it easier to manage them.
 *
 * Inside-the-layer pressure has a complicated definition.
 * See Hsu and Arakawa 1990 (5.43) or Konor and Arakawa (3.34).
 *
 * An illustration of the kk index convention is:
 *   kk = 3 => K = 1.5, the bottom of layer 1, and
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
  register int
    K;
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

#if EPIC_CHECK == 1
  /*
   * Check validity of kk.
   *
   * NOTE: EPIC_CHECK == 1 when the environment variable
   *       EPIC_CFLAG is set to -g.
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

  if (kk > 2*grid.nk+1) {
    sprintf(Message,"kk > 2*grid.nk+1");
    epic_error(dbmsname,Message);
  }

  if (kk%2 == 0) {
    /* Layer value. */
    K = kk/2;

    switch(index) {
      case P2_INDEX:
        /*
         * Inside-the-layer pressure has a complicated definition.
         * See Hsu and Arakawa 1990 (5.43) or Konor and Arakawa (3.34).
         */
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
      break;
      case U_INDEX:
        return .5*(botval+topval);
      break;
      case V_INDEX:
        return .5*(botval+topval);
      break;
      case HDRY_INDEX:
        sprintf(Message,"kk=%d, index=%d defined in layer, use directly",kk,index);
        epic_error(dbmsname,Message);
      break;
      case THETA_INDEX:
      case THETA2_INDEX:
        /*
         * If averaging is called for, use plain averaging as in Konor and Arakawa (1997, eqn 3.26).
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
    K = (kk-1)/2;

    switch(index) {
      case U_INDEX:
      case V_INDEX:
      case P3_INDEX:
      case THETA_INDEX:
        sprintf(Message,"kk=%d, index=%d defined on interface, use directly",kk,index);
        epic_error(dbmsname,Message);
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
}

/*======================= end of onto_kk() ===================================*/

/*======================= get_kin() ==========================================*/

/*
 * Return kinetic energy per mass. 
 * Kinetic energy resides on the h-grid.
 */

EPIC_FLOAT get_kin(planetspec *planet,
                   EPIC_FLOAT *u2d,
                   EPIC_FLOAT *v2d,
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
    u2,u4,v2,v4;
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
      mn_inv[jj] = 1./(grid.mn)[jj];
    }

    if (strcmp(grid.geometry,"f-plane") == 0 &&
        strcmp(grid.f_plane_map,"cartesian")) {
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

#if PV_SCHEME == ARAKAWA_LAMB_1981
  /*
   * See Arakawa and Lamb (1981), eq. (A33).
   *
   * NOTE: The Hollingsworth-Kallberg instability is not corrected in this 
   *       implementation, see below.
   */
  jj    = 2*J;
  u2    = U2D(J,  I  );
  u4    = U2D(J,  I+1);
  v2    = V2D(J,  I  );
  v4    = V2D(J+1,I  );
  kin_c = (                 u2*u2
                           +u4*u4
            +grid.mn[jj+1]*(v2*v2*mn_inv[jj  ]
                           +v4*v4*mn_inv[jj+2]) )*.25;
  kin   = kin_c;

#if EPIC_CHECK == 1
  /*
   * Screen for nan.
   */
  if (!isfinite(kin)) {
    sprintf(Message,"JI=%d %d, kin=%g, u2=%g u4=%g v2=%g v4=%g; PV_SCHEME=ARAKAWA_LAMB_1981",
                    J,I,kin,u2,u4,v2,v4);
    epic_error(dbmsname,Message);
  }
#endif

#elif PV_SCHEME == SADOURNEY_1975
  /*
   * The Hollingsworth-Kallberg instability arises when there is incomplete 
   * cancellation of terms in the discrete horizontal momentum equations when 
   * written in vector-invarient form. The term -u*du/dy from (pv)*uh on the 
   * left-hand side of the dvdt equation does not exactly cancel the -d(u*u/2)/dy 
   * term from dK/dy on the right-hand side, and the same for the (pv)*vh and dK/dx 
   * terms in the dudt equation.
   *
   * We follow the method described by Suarez and Takacs (1995, NASA Technical 
   * Memorandum 104606, Vol. 5) for eliminating this Hollingsworth-Kallberg
   * instability when using the Sadourney (1975) averaging scheme for pv
   * in the (pv)*uh and (pv)*vh terms.
   */
  jj    = 2*J;
  u2    = U2D(J,  I  );
  u4    = U2D(J,  I+1);
  v2    = V2D(J,  I  );
  v4    = V2D(J+1,I  );
  kin_c = (                 u2*u2
                           +u4*u4
            +grid.mn[jj+1]*(v2*v2*mn_inv[jj  ]
                           +v4*v4*mn_inv[jj+2]) )*.25;

  if (!j_periodic && (J == grid.jlo  || J == grid.nj)) {
    /*
     * For simplicity, at poles or channel boundaries use kin = kin_c.
     */
    kin = kin_c;

#if EPIC_CHECK == 1
    /*
     * Screen for nan.
     */
    if (!isfinite(kin)) {
      sprintf(Message,"JI=%d %d, kin=%g, u2=%g u4=%g v2=%g v4=%g; PV_SCHEME=SADOURNEY_1975",
                      J,I,kin,u2,u4,v2,v4);
      epic_error(dbmsname,Message);
    }
#endif

  }
  else {
    u2    = .5*(U2D(J+1,I  )+U2D(J-1,I  ));
    u4    = .5*(U2D(J+1,I+1)+U2D(J-1,I+1));
    v2    = .5*(V2D(J,  I+1)+V2D(J,  I-1));
    v4    = .5*(V2D(J+1,I+1)+V2D(J+1,I-1));
    kin_s = (                 u2*u2
                             +u4*u4
              +grid.mn[jj+1]*(v2*v2*mn_inv[jj  ]
                             +v4*v4*mn_inv[jj+2]) )*.25;
    kin   = (5./6.)*kin_c+(1./6.)*kin_s;

#if EPIC_CHECK == 1
    /*
     * Screen for nan.
     */
    if (!isfinite(kin)) {
      sprintf(Message,"JI=%d %d, kin=%g, u2=%g u4=%g v2=%g v4=%g kin_s=%g kin_c=%g; PV_SCHEME=SADOURNEY_1975",
                      J,I,kin,u2,u4,v2,v4,kin_s,kin_c);
      epic_error(dbmsname,Message);
    }
#endif

  }

#else
  sprintf(Message,"unrecognized PV_SCHEME");
  epic_error(dbmsname,Message);
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
    fpara = 0.25;
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

  g = grid.g[2*J+1];

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
    dz_inv;
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
    K      = kk/2;
    dz_inv = grid.g[2*J+1]/(GZ3(K-1,J,I)-GZ3(K,J,I));

    dudz2  = .5*(U(grid.it_uv,K-1,J,I)+U(grid.it_uv,K-1,J,I+1)
                -U(grid.it_uv,K,  J,I)-U(grid.it_uv,K,  J,I+1));

    dvdz2  = .5*(V(grid.it_uv,K-1,J,I)+V(grid.it_uv,K-1,J+1,I)
                -V(grid.it_uv,K,  J,I)+V(grid.it_uv,K,  J+1,I));
  }
  else {
    /*
     * Interior interface value.
     */
    K      = (kk-1)/2;
    dz_inv = grid.g[2*J+1]/(GZ3(K-1,J,I)-GZ3(K+1,J,I));

    dudz2  = .5*(U(grid.it_uv,K-1,J,I)+U(grid.it_uv,K-1,J,I+1)
                -U(grid.it_uv,K+1,J,I)-U(grid.it_uv,K+1,J,I+1));

    dvdz2  = .5*(V(grid.it_uv,K-1,J,I)+V(grid.it_uv,K-1,J+1,I)
                -V(grid.it_uv,K+1,J,I)+V(grid.it_uv,K+1,J+1,I));
  }

  dudz2 *= dudz2;
  dvdz2 *= dvdz2;

  if (fcmp(dudz2,0.) != 0) {
    return get_brunt2(planet,kk,J,I)/(dudz2+dvdz2);
  }
  else {
    /*
     * Return a large but finite number for Ri = infinity case.
     */
    return 1.e+6;
  }
}

/*====================== end of get_richardson() =============================*/

/*====================== get_sounding() ======================================*/
/*
 * Smooth, monotonic interpolation of t_vs_p data table.
 *
 * Valid output_name values: "temperature" 
 *                           "theta"
 *
 * For example, to get a theta(p) via interpolation, use 
 *   get_sounding(planet,p,"theta",&theta);
 *
 * Interpolation uses -log(p).
 *
 * NOTE: If input_value is below the bottom of the data table, then the bottom
 *       value of the data table is used. 
 */
void get_sounding(planetspec *planet,
                  EPIC_FLOAT  pressure,
                  char       *output_name,
                  EPIC_FLOAT *pt_output_value)
{
  int
    ki;
  static int
    initialized = FALSE;
  EPIC_FLOAT
    fpara,fgibb,fpe,uoup,
    theta_ortho,theta_para,
    x,x_d,
    tmp;
  static float_triplet
    *tdat,
    *thetadat;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_sounding";

  if (!initialized) {
    /* 
     * Allocate memory: 
     */
    tdat     = ftriplet(0,var.ntp-1,dbmsname);
    thetadat = ftriplet(0,var.ntp-1,dbmsname);

    /*
     * Set temperature data table.
     */
    for (ki = 0; ki < var.ntp; ki++) {
      tdat[ki].x = -log(var.pdat[ki]);
      tdat[ki].y = var.tdat[ki];
    }

    spline_pchip(var.ntp,tdat);

    /*
     * Set theta data table.
     */
    for (ki = 0; ki < var.ntp; ki++) {
      thetadat[ki].x = -log(var.pdat[ki]);
      if (var.fpara.on) {
        fpara = return_fpe(var.tdat[ki]);
      }
      else {
        fpara = .25;
      }
      thetadat[ki].y = return_theta(planet,fpara,var.pdat[ki],var.tdat[ki],&theta_ortho,&theta_para);
    }

    spline_pchip(var.ntp,thetadat);

    /*
     * From the top down, iron out negative-slope regions in theta.
     * Adjust temperature accordingly.
     */
    if (var.ntp-2 >= 0) {
      for (ki = var.ntp-2; ki >= 0; ki--) {
        tmp = thetadat[ki+1].y;
        if (thetadat[ki].y > tmp) {
          thetadat[ki].y = tmp;
          if (var.fpara.on) {
            fpara = return_fpe(var.tdat[ki]);
          }
          else {
            fpara = .25;
          }
          /*
           * Adjust tdat.
           */
          tdat[ki].y = return_temp(planet,fpara,var.pdat[ki],thetadat[ki].y);
        }
      }
    } 
    initialized = TRUE;
  } 
  /* End of initialization. */

  /*
   * Interpolate using -log(p).
   */
  x = -log(pressure);

  /*
   * Restrict y(x) to not fall below bottom of data table.
   */
  if (x < tdat[0].x) {
    if (strcmp(output_name,"temperature") == 0) {
      *pt_output_value = tdat[0].y;
    }
    else if (strcmp(output_name,"theta") == 0) {
      *pt_output_value = thetadat[0].y;
    }
    else {
      sprintf(Message,"unrecognized output_name=%s",output_name);
      epic_error(dbmsname,Message);
    }
    sprintf(Message,"-ln(p)=%e < %e; setting %s[0]=%e \n",
                     x,tdat[0].x,output_name,*pt_output_value);
    epic_warning(dbmsname,Message);
  }
  else {
    if (strcmp(output_name,"temperature") == 0) {
      ki               = find_place_in_table(var.ntp,tdat,x,&x_d);
      *pt_output_value = splint_pchip(x,tdat+ki,x_d);
    }
    else if (strcmp(output_name,"theta") == 0) {
      ki               = find_place_in_table(var.ntp,thetadat,x,&x_d);
      *pt_output_value = splint_pchip(x,thetadat+ki,x_d);
    }
  }

  return;
}

/*======================= end of get_sounding() ==============================*/

/*======================= return_cp() ========================================*/

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
  if (strcmp(planet->name,"held_suarez") == 0) {
    cp = 1004.;
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

/*======================= end of return_cp() =================================*/

/*======================= set_p2etc() ========================================*/

/*
 * The diagnostic variables HDRY3, PDRY3, P2, P3, and THETA2 are updated
 * here rather than in store_diag(), because they are directly 
 * related to the prognostic variables HDRY and THETA and need to be updated
 * more often than the other diagnostic variables.
 */

void set_p2etc(planetspec *planet,
               int         theta_switch)
{
  register int
    K,J,I,
    kk,iq;
  int
    itmp;
  register EPIC_FLOAT
    sigma,
    sum,
    ptop,pbot,
    tmp,
    g;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_p2etc";

  if (grid.nq > 0) {
    /*
     * Calculate H2 (total hybrid density) from HDRY and Q_i (mass mixing ratios).
     */
    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          sum = 1.;
          for (iq = 0; iq < grid.nq; iq++) {
            /*
             * NOTE: Q is carried on layer interfaces, so use get_var().
             */
            sum += get_var(planet,grid.is[iq],grid.ip[iq],grid.it_h,kk,J,I);
          }
          H2(K,J,I) = HDRY(K,J,I)*sum;
        }
      }
    }
  }

  /*
   * NOTE: In the dry case, grid.nq = 0 and the memory for H2 is the same as HDRY; 
   *       this also applies to P3 and PDRY3.
   */

  for (J = JLOPAD; J <= JHIPAD; J++) {
    g = grid.g[2*J+1];
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /*
       * P3 is a top-down integral of H2; the model's top pressure is fixed.
       * whereas pbot varies. 
       */
      for (K = KLO; K <= KHI; K++) {
        kk        = 2*K;
        P3(K,J,I) = P3(K-1,J,I)+grid.dsgth[kk]*g*H2(K,J,I);
      }

      /*
       * Apply the diagnostic value of pressure to P3 in the pure-sigma portion of the model.
       */
      ptop = P3(0,  J,I);
      pbot = P3(KHI,J,I);
      for (K = grid.k_sigma-1; K < KHI; K++) {
        sigma     = (grid.sigmatheta[2*K+1]-grid.zeta0)/(grid.zeta1-grid.zeta0);
        P3(K,J,I) = get_p_sigma(pbot,sigma,ptop);
      }

      /*
       * Update H2.
       */
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K;
        H2(K,J,I) = MAX(grid.h_min[K],get_h(planet,kk,J,I,TOTAL));
      }
      if (grid.nq > 0) {
        /*
         * HDRY and H2 are separate memory, so update HDRY.
         */
        for (K = KLO; K <= KHI; K++) {
          kk = 2*K;
          HDRY(K,J,I) = get_h(planet,kk,J,I,DRY);
        }
      }
    }
  }

  if (grid.nq > 0) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      g = grid.g[2*J+1];
      for (I = ILOPAD; I <= IHIPAD; I++) {
        /*
         * PDRY3 is a straightforward integral of HDRY.
         */
        for (K = KLO; K <= KHI; K++) {
          kk = 2*K;
          PDRY3(K,J,I) = PDRY3(K-1,J,I)+grid.dsgth[kk]*g*HDRY(K,J,I);
        }
      }
    }
  }

  if (theta_switch == UPDATE_THETA) {
    /*
     * Apply the diagnostic value of potential temperature to THETA outside the pure-sigma region.
     */
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        ptop = P3(0,  J,I);
        pbot = P3(KHI,J,I);
        for (K = 0; K < grid.k_sigma-1; K++) {
          sigma        = get_sigma(pbot,P3(K,J,I),ptop);
          THETA(K,J,I) = (EPIC_FLOAT)(((double)grid.sigmatheta[2*K+1]-(double)f_sigma(sigma))/(double)g_sigma(sigma));
        }
      }
    }
  }

  /*
   * Calculate P2, H3, HDRY3.
   */
  for (J = JLOPAD; J <=JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K;
        /*
         * P2 depends on P3.
         */
        P2(K,J,I) = get_p(planet,P2_INDEX,kk,J,I);
      }

      /*
       * Top and bottom of model.
       */
      K = 0;
      P2(K,J,I) = P3(K,J,I)*P3(K,J,I)/P2(K+1,J,I);

      for (K = KLO; K <= KHI; K++) {
        kk        = 2*K+1;
        H3(K,J,I) = get_h(planet,kk,J,I,TOTAL);
      }
      H3(KLO-1,J,I) = H3(KLO,J,I);
 
      if (grid.nq > 0) {
        for (K = KLO; K <= KHI; K++) {
          kk           = 2*K+1;
          HDRY3(K,J,I) = get_h(planet,kk,J,I,DRY);
        }
        HDRY3(KLO-1,J,I) = HDRY3(KLO,J,I);
      }
    }
  }

  /*
   * Set THETA2.
   */
  if (theta_switch == UPDATE_THETA) {
    int
      kay,
      nkay = KHI-(grid.k_sigma-1)+1;
    float_triplet
      table[nkay];
    FLOAT
      sgth,h;

    for (K = KHI; K >= grid.k_sigma-1; K--) {
      kay          = KHI-K;
      table[kay].x = grid.sigmatheta[2*K+1];
    }

    for (J = JLOPAD; J <=JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        ptop = P3(0,  J,I);
        pbot = P3(KHI,J,I);
        for (K = KLO; K < grid.k_sigma-1; K++) {
          /*
           * THETA2 is known diagnostically.
           *
           * NOTE: The value of g_sigma(sigma) just above
           *       sigma_sigma can sometimes become too small for accurate division.
           *       We find that averaging THETA to get THETA2(grid.k_sigma-1) in practice 
           *       alleviates the problem (whereas, extending the spline up to cover
           *       this we find generates a numerical instability).
           */
          sigma         = get_sigma(pbot,P2(K,J,I),ptop);
          THETA2(K,J,I) = (EPIC_FLOAT)(((double)grid.sigmatheta[2*K]-(double)f_sigma(sigma))/(double)g_sigma(sigma));
        }
        /*
         * This is the averaging described in the previous comment.
         */
        K = grid.k_sigma-1;
        THETA2(K,J,I) = .5*(THETA(K,J,I)+THETA(K-1,J,I));

        /*
         * NOTE: We find that using a spline on live variables like THETA, even a smooth, monotonic
         *       one like spline_pchip(), tends to expose the model to unwanted computational modes. 
         *       Instead, we use linear interpolation, which is more localized and less prone to numerical
         *       instability.
         */
        for (K = KHI; K >= grid.k_sigma-1; K--) {
          kay          = KHI-K;
          table[kay].y = THETA(K,J,I);
        }

        for (K = grid.k_sigma; K <= KHI; K++) {
          sgth          = grid.sigmatheta[2*K];
          kay           = find_place_in_table(nkay,table,sgth,&h);
          THETA2(K,J,I) = linint(sgth,table+kay,h);
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  return;
}

/*======================= end of set_p2etc() =================================*/

/*==================== store_pgrad_vars() ====================================*/

/*
 * Calculates and stores diagnostic variables that are typically used in the
 * pressure-gradient terms of the momentum equations.  Doing this in a
 * separate function facilitates the "poor man's implicit" timestep.
 *
 * This is also a convenient place to calculate the regular vertical velocity,
 * DZDT3, in m/s (the hybrid vertical velocity is calculated in calc_w()).
 *
 * An illustration of the notation with respect to the vertical index K:
 *
 *              --------------------------------------------------
 *     Layer K:   T2   RHO2          
 * Interface K: --T3---RHO3---EXNER3---GZ3--MONT3--FGIBB3--DZDT3--
 *
 * NOTE: This function should not modify any prognostic variables. 
 *       In addition, it should not modify P2, P3, or THETA2.
 *
 * NOTE: This function must be called from all nodes.
 */

void store_pgrad_vars(planetspec *planet) 
{
  register int
    K,J,I,
    J0,I0,
    is,ip,
    kk,kay,nkay,
    jj;
  register EPIC_FLOAT
    theta,
    fpara,
    pressure,
    temperature,
    enthalpy,
    mu,
    avg,pbot,ptop,sigma,
    dtg_inv,gdx_inv3,gdy_inv2,gdy_inv4,
    dz30,dz31,dz32,dz20,dz21,dz10,
    df0,df1,df2,df3,
    gz,dgz,
    sgth;
  EPIC_FLOAT
    fgibb,fpe,uoup,
    dsgth;
  float_triplet
    table[KHI-KLO+2];
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    *int_exner_m,
    *exner_m,
    *sendbuf,
    *u1d,*p1d,*rho1d,*gz1d;

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
    dbmsname[]="store_pgrad_vars";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    int_exner_m = fvector(grid.jlo,grid.nj,dbmsname);
    exner_m     = fvector(grid.jlo,grid.nj,dbmsname);
    sendbuf     = fvector(grid.jlo,grid.nj,dbmsname);
    u1d         = fvector(0,JADIM-1,dbmsname);
    p1d         = fvector(0,JADIM-1,dbmsname);
    rho1d       = fvector(0,JADIM-1,dbmsname);
    gz1d        = fvector(0,JADIM-1,dbmsname);

    initialized = TRUE;
  }

  /*
   * Store original value of GZ3 in DZDT3, as the first
   * step to calculating the partial derivative of z with respect to time.
   */
  if (grid.itime > 1) {
    for (K = KLO; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          DZDT3(K,J,I) = GZ3(K,J,I);
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }
  else {
    /*
     * Zero out DZDT3(K,J,I).
     */
    memset(var.dzdt3.value,0,Nelem3d*sizeof(EPIC_FLOAT));
  }

  /*
   * Calculate T3, EXNER3, RHO3.
   *
   * NOTE: We use EXNER = cp*T/THETA, which is more general than cp*(p/p0)^kappa,
   *       where cp = planet->cp is a reference value.
   */
  if (var.fpara.on) {
    for (K = 0; K <= KHI; K++) {
      kk = 2*K+1;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          theta         = THETA(K,J,I);  
          fpara         = FPARA(K,J,I);
          pressure      = P3(K,J,I);
          T3(K,J,I)     = temperature = return_temp(planet,fpara,pressure,theta);
          EXNER3(K,J,I) = planet->cp*T3(K,J,I)/theta;
          mu            = avg_molar_mass(planet,kk,J,I);
          RHO3(  K,J,I) = return_density(planet,fpara,pressure,temperature,mu,PASSING_T);
          return_enthalpy(planet,fpara,pressure,temperature,&fgibb,&fpe,&uoup);
          FGIBB3(K,J,I) = fgibb;
        }
      }
      /* No need to apply bc_lateral() here. */
    }
  }
  else {
    fpara = 0.25;
    for (K = 0; K <= KHI; K++) {
      kk = 2*K+1;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          theta         = THETA(K,J,I);    
          pressure      = P3(K,J,I);
          T3(K,J,I)     = temperature = return_temp(planet,fpara,pressure,theta);
          EXNER3(K,J,I) = planet->cp*T3(K,J,I)/theta;
          mu            = avg_molar_mass(planet,kk,J,I);
          RHO3(  K,J,I) = return_density(planet,fpara,pressure,temperature,mu,PASSING_T);
        }
      }
      /* No need to apply bc_lateral() here. */
    }
  }

  /*
   * Calculate T2, which depends on FPARA, P2, THETA2.
   * Calculate RHO2, which depends on FPARA, P2, T2.
   */
  if (var.fpara.on) {
    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          fpara         = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,kk,J,I);
          pressure      = P2(K,J,I);
          theta         = THETA2(K,J,I);
          T2(    K,J,I) = temperature = return_temp(planet,fpara,pressure,theta);
          mu            = avg_molar_mass(planet,kk,J,I);
          RHO2(  K,J,I) = return_density(planet,fpara,pressure,temperature,mu,PASSING_T);
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }
  else {
    fpara = 0.25;
    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          pressure      = P2(K,J,I);
          theta         = THETA2(K,J,I);
          T2(    K,J,I) = temperature = return_temp(planet,fpara,pressure,theta);
          mu            = avg_molar_mass(planet,kk,J,I);
          RHO2(  K,J,I) = return_density(planet,fpara,pressure,temperature,mu,PASSING_T);
        }
      }
      /* No need to apply bc_lateral() here. */
    }
  }

  /*
   * Assume constant ratios for density.
   * Assume temperature does not change with height above top of model.
   */
  K  = 0;
  kk = 2*K;
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      RHO2(  K,J,I) = RHO3(K,J,I)*RHO3(K,J,I)/RHO2(K+1,J,I);
      THETA2(K,J,I) = theta = grid.sigmatheta[kk];
      temperature   = T3(K,J,I);
      T2(    K,J,I) = temperature;
    } 
  }
  /* No need to apply bc_lateral() here. */

  if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * For gas giants, update GZ_SURFACE, which is determined from
     * gradient balance with U at the bottom interface of the model,
     * based on zonal averages.
     */
    K = grid.nk;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      /*
       * Calculate zonal averages.
       *
       * NOTE: Not ready for cut in I direction.
       */
      U1D(  J) = 0.;
      P1D(  J) = 0.;
      RHO1D(J) = 0.;
      for (I = ILO; I <= IHI; I++) {
        U1D(  J) += U(grid.it_uv,K,J,I);
        P1D(  J) += P3(  K,J,I);
        RHO1D(J) += RHO3(K,J,I);
      }
      U1D(  J) /= grid.ni;
      P1D(  J) /= grid.ni;
      RHO1D(J) /= grid.ni;
    }
    gz_from_u(planet,u1d,p1d,rho1d,gz1d,grid.jtp,0.);
    for (J = JLOPAD; J <= JHIPAD; J++) {
      gz = GZ1D(J);
      for (I = ILOPAD; I <= IHIPAD; I++) {
        GZ_SURFACE(J,I) = gz;
      }
    }
  }

  /*
   * Calculate the geopotential, GZ3, assuming hydrostatic balance,
   * and the corresponding Montgomery potential MONT3.
   *
   * Start with the surface.
   */

  K = KHI;
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /* 
       * The artificial parameter grid.topo_scale can be used to
       * suppress the topography.
       */
      GZ3(K,J,I)   = grid.topo_scale*GZ_SURFACE(J,I);

      fpara        = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,2*K+1,J,I);
      MONT3(K,J,I) = GZ3(K,J,I)+return_enthalpy(planet,fpara,P3(K,J,I),T3(K,J,I),&fgibb,&fpe,&uoup);
    }
  }
  /* No need to apply bc_lateral() here. */

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /* Short centered step up from bottom. */
      K = KHI;
      GZ3(K-1,J,I)   = GZ3(K,J,I)+((P3(K,J,I)+P3(K-1,J,I))/(RHO3(K,J,I)+RHO3(K-1,J,I)))*log(P3(K,J,I)/P3(K-1,J,I));

      fpara          = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,2*(K-1)+1,J,I);
      MONT3(K-1,J,I) = GZ3(K-1,J,I)+return_enthalpy(planet,fpara,P3(K-1,J,I),T3(K-1,J,I),&fgibb,&fpe,&uoup);

      for (K = KHI-1; K >= KLO+1; K--) {
        /* 
         * This accurate integration scheme is based on the one described and illustrated in Fig. 8(d) of 
         *   Leslie LM,  Purser RJ, 1992,  A comparative study of the performance of various
         *   vertical discretization schemes, Meteorol. Atmos. Phys. 50, 61-73.
         */
        dz10 = log(P3(K+1,J,I)/P3(K,  J,I));
        dz21 = log(P3(K,  J,I)/P3(K-1,J,I));
        dz32 = log(P3(K-1,J,I)/P3(K-2,J,I));
        dz20 = dz21+dz10;
        dz30 = dz32+dz20;
        dz31 = dz32+dz21;

        df0  = P3(K+1,J,I)/RHO3(K+1,J,I);
        df1  = P3(K  ,J,I)/RHO3(K,  J,I);
        df2  = P3(K-1,J,I)/RHO3(K-1,J,I);
        df3  = P3(K-2,J,I)/RHO3(K-2,J,I);

        dgz = -(df0*(dz21*dz21*dz21*(dz31+dz32))/(dz10*dz20*dz30)
               +df1*(dz21*(dz21*(dz21+2.*(dz10-dz31))-6.*dz10*dz31))/(dz10*dz31)
               +df2*(dz21*(dz21*(dz21+2.*(dz32-dz20))-6.*dz20*dz32))/(dz20*dz32)
               +df3*(dz21*dz21*dz21*(dz10+dz20))/(dz30*dz31*dz32)               )/12.;

        GZ3(K-1,J,I) = GZ3(K,J,I)+dgz;

        fpara          = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,2*(K-1)+1,J,I);
        MONT3(K-1,J,I) = GZ3(K-1,J,I)+return_enthalpy(planet,fpara,P3(K-1,J,I),T3(K-1,J,I),&fgibb,&fpe,&uoup);
      }

      /* Short centered step to top. */
      K = KLO;
      GZ3(K-1,J,I)   = GZ3(K,J,I)+((P3(K,J,I)+P3(K-1,J,I))/(RHO3(K,J,I)+RHO3(K-1,J,I)))*log(P3(K,J,I)/P3(K-1,J,I));

      fpara          = get_var(planet,FPARA_INDEX,NO_PHASE,grid.it_h,2*(K-1)+1,J,I);
      MONT3(K-1,J,I) = GZ3(K-1,J,I)+return_enthalpy(planet,fpara,P3(K-1,J,I),T3(K-1,J,I),&fgibb,&fpe,&uoup);
    }
  }

  /*
   * Screen for non-monotonic GZ.
   */
  for (K = KLO; K <= KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        if (GZ3(K,J,I) >= GZ3(K-1,J,I)) {
          sprintf(Message,"GZ3(%d,%d,%d)=%g >= GZ3(%d,%d,%d)=%g",
                           K,J,I,GZ3(K,J,I),K-1,J,I,GZ3(K-1,J,I));
          epic_error(dbmsname,Message);
        }
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  /*
   * Finish the calculation of the partial derivative of z with respect to time.
   */
  if (grid.itime > 1) {
    for (K = KLO; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        dtg_inv = -1./(DT*grid.g[2*J+1]);
        for (I = ILO; I <= IHI; I++) {
          DZDT3(K,J,I) -= GZ3(K,J,I);
          DZDT3(K,J,I) *= dtg_inv;
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  /*
   * Add to DZDT3 the hybrid vertical velocity term. 
   */
  for (K = KLO; K <= KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      dtg_inv = -1./(DT*grid.g[2*J+1]);
      for (I = ILO; I <= IHI; I++) {
        DZDT3(K,J,I) += W3(K,J,I)*H3(K,J,I)/RHO3(K,J,I);
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  /*
   * Finish the calculation of DZDT3 by adding the partial derivatives
   * with respect to latitude and longitude.
   */
  for (K = KLO; K <= KHI; K++) {
    kk = 2*K+1;
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J;
      gdx_inv3 = .5*grid.m[jj+1]/grid.g[jj+1];
      gdy_inv2 = .5*grid.n[jj  ]/grid.g[jj  ];
      gdy_inv4 = .5*grid.n[jj+2]/grid.g[jj+2];
      for (I = ILO; I <= IHI; I++) {
        DZDT3(K,J,I) += get_var(planet,V_INDEX,NO_PHASE,grid.it_uv,kk,J,  I)*(GZ3(K,J,  I)-GZ3(K,J-1,I))*gdy_inv2
                       +get_var(planet,V_INDEX,NO_PHASE,grid.it_uv,kk,J+1,I)*(GZ3(K,J+1,I)-GZ3(K,J,  I))*gdy_inv4
                       +get_var(planet,U_INDEX,NO_PHASE,grid.it_uv,kk,J,  I)*(GZ3(K,J,  I)-GZ3(K,J,I-1))*gdx_inv3
                       +get_var(planet,U_INDEX,NO_PHASE,grid.it_uv,kk,J,I+1)*(GZ3(K,J,I+1)-GZ3(K,J,  I))*gdx_inv3;
      }
    }
  }
  /* Need to call bc_lateral() here. */
  bc_lateral(var.dzdt3.value,THREEDIM);

  return;
}

/*==================== end of store_pgrad_vars() =============================*/

/*======================= store_diag() =======================================*/

/*
 * Calculates and stores commonly used diagnostic variables other than
 * those covered by set_p2etc() and store_pgrad_vars().
 *
 * NOTE: This function must be called from all nodes.
 */

void store_diag(planetspec *planet) 
{
  register int
    K,J,I,
    kk;
  static EPIC_FLOAT
    *buff2d;
  static int
    initialized = FALSE;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="store_diag";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    buff2d = fvector(0,Nelem2d-1,dbmsname);

    initialized = TRUE;
  }

  /*
   * Calculate potential vorticity.
   */
  for (K = KLO-1; K <= KHI; K++) {
    /* 
     * Calculate potential vorticity on the layer interfaces, PV3.
     * Use the total hybrid density, H3 (see Schubert et al 2001, JAS 58, 3148-57).
     *
     * Before calculating vorticity, apply zonal_filter() to a copy of V to
     * control numerical instability for dv/dx.
     *
     * NOTE: The same type of hybrid density used for PV3 here must also be used 
     * for UH and VH in uv_core() to get the correct Coriolis terms.
     */
    memcpy(buff2d,var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,Nelem2d*sizeof(EPIC_FLOAT));
    zonal_filter(V_INDEX,buff2d,NULL,TWODIM);
    vorticity(planet,ON_SIGMATHETA,POTENTIAL,
              var.u.value  +(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
              buff2d,
              var.h3.value +(K-Kshift)*Nelem2d+grid.it_h*Nelem3d,
              var.pv3.value+(K-Kshift)*Nelem2d);
  }

  /*
   * Calculate horizontal divergence.
   * Apply zonal_filter() to a copy of U to control numerical instability of du/dx.
   */
  for (K = KLO-1; K <= KHI; K++) {
    memcpy(buff2d,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,Nelem2d*sizeof(EPIC_FLOAT));
    zonal_filter(U_INDEX,buff2d,NULL,TWODIM);
    divergence(buff2d,
               var.v.value      +(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
               var.div_uv3.value+(K-Kshift)*Nelem2d);
  }

  /*
   * Calculate Richardson number array.
   */
  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      for (K = KLO; K <= KHI; K++) {
        kk         = 2*K;
        RI2(K,J,I) = get_richardson(planet,kk,  J,I);
      }   
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.ri2.value,THREEDIM);

  /*
   * Calculate turbulence-model variables.
   */
  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
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

void divergence(EPIC_FLOAT *uu,
                EPIC_FLOAT *vv,
                EPIC_FLOAT *di)
{
  int
    J,I,
    jj;
  EPIC_FLOAT
    m_2jp1,m_2j_inv,m_2jp2_inv,
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
    m_2jp1     = (grid.m)[jj];
    n_2jp1     = (grid.n)[jj];
    n_2jp1_inv = 1./n_2jp1;
    mn_2jp1    = (grid.mn)[jj];

    if (J == grid.jlo && IS_SPOLE) {
      m_2j_inv = 0.;
    }
    else {
      m_2j_inv = 1./(grid.m)[jj-1];
    }

    if (J == grid.nj && IS_NPOLE) {
      m_2jp2_inv = 0.;
    }
    else {
      m_2jp2_inv = 1./(grid.m)[jj+1];
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

void vorticity(planetspec *planet,
               int         surface_type,
               int         type,
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
    ze,zetabot,zetatop,
    h_pv,pv_pole,pvbot,pvtop;
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
    m_2j        =    (grid.m )[2*J  ];
    n_2j        =    (grid.n )[2*J  ];
    m_2jm1_inv  = 1./(grid.m )[2*J-1];
    m_2jp1_inv  = 1./(grid.m )[2*J+1];
    mn_2jm1_inv = 1./(grid.mn)[2*J-1];
    mn_2jp1_inv = 1./(grid.mn)[2*J+1];
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
      m_2j        =    (grid.m )[2*J  ];
      n_2j        =    (grid.n )[2*J  ];
      m_2jm1_inv  = 1./(grid.m )[2*J-1];
      m_2jp1_inv  = 1./(grid.m )[2*J+1];
      /* Calculate average zeta in next row */
      /* 
       * NOTE: This global sum is fine so long as the zonal direction is
       *       not decomposed.
       */
      zetabot = 0.;
      for (I = ILO; I <= IHI; I++) {
        zetabot += m_2j*( (VV(J,I)-VV(J,I-1)) +
                   n_2j*(UU(J-1,I)*m_2jm1_inv-UU(J,I)*m_2jp1_inv) );
      }
      zetabot /= (IHI-ILO+1);

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
        h_pv  /= (IHI-ILO+1);
        pvbot  = (zetabot+grid.f[2*J])/h_pv;
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
        ze  += UU(J,I);
      }
      ze *= (grid.mn)[2*(J+1)]/((grid.m)[2*J+1]*(EPIC_FLOAT)(grid.ni));

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
        h_pv    /= (EPIC_FLOAT)(grid.ni);
        pv_pole  = (ze+grid.f[2*(J+1)])/h_pv;
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
          ze  -= UU(J,I);  /*  Beware of southern circulation sign. */
        }
        ze *= (grid.mn)[0]/((grid.m[1])*(EPIC_FLOAT)(grid.ni)); 

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
          h_pv    /= (EPIC_FLOAT)(grid.ni);
          pv_pole  = (ze+grid.f[0])/h_pv;
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
          m_2j        =    (grid.m )[2*J  ];
          n_2j        =    (grid.n )[2*J  ];
          m_2jm1_inv  = 1./(grid.m )[2*J-1];
          m_2jp1_inv  = 1./(grid.m )[2*J+1];
          /* 
           * NOTE: This global sum is fine so long as the zonal direction is
           *       not decomposed.
           */
          zetabot = 0.;
          for (I = ILO; I <= IHI; I++) {
            zetabot +=  m_2j*( (VV(J,I)-VV(J,I-1)) +
                        n_2j*(UU(J-1,I)*m_2jm1_inv-UU(J,I)*m_2jp1_inv) );
          }
          zetabot /= (IHI-ILO+1);

          if (type == RELATIVE) {
            pvbot = zetabot;
          }
          else if (type == ABSOLUTE) {
            pvbot = zetabot+grid.f[2*J];
          }
          else if (type == POTENTIAL) {
            h_pv = 0.;
            J    = 0;
            for (I = ILO; I <= IHI; I++) {
              h_pv += HH(J,I);
            }
            h_pv  /= (IHI-ILO+1);
            pvbot  = (zetabot+grid.f[2*J])/h_pv;
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
        ze *= (grid.mn)[2*(J+1)]/((grid.m)[2*J+1]*(EPIC_FLOAT)(grid.ni));

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
          h_pv    /= (EPIC_FLOAT)(grid.ni);
          pv_pole  = (ze+grid.f[2*(J+1)])/h_pv;
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
          m_2j        =    (grid.m )[2*J  ];
          n_2j        =    (grid.n )[2*J  ];
          m_2jm1_inv  = 1./(grid.m )[2*J-1];
          m_2jp1_inv  = 1./(grid.m )[2*J+1];
          /* 
           * NOTE: This global sum is fine so long as the zonal direction is
           *       not decomposed.
           */
          zetatop = 0.;
          for (I = ILO; I <= IHI; I++) {
            zetatop +=  m_2j*( (VV(J,I)-VV(J,I-1)) +
                        n_2j*(UU(J-1,I)*m_2jm1_inv-UU(J,I)*m_2jp1_inv) );
          }
          zetatop /= (IHI-ILO+1);

          if (type == RELATIVE) {
            pvtop = zetatop;
          }
          else if (type == ABSOLUTE) {
            pvtop = zetatop+grid.f[2*J];
          }
          else if (type == POTENTIAL) {
            h_pv = 0.;
            for (I = ILO; I <= IHI; I++) {
              h_pv += HH(J,I);
            }
            h_pv  /= (IHI-ILO+1);
            pvtop  = (zetatop+grid.f[2*J])/h_pv;
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

/*======================= gz_from_u() ==========================================*/

/*
 * Calculate geopotential, gz, assuming gradient balance with a zonally 
 * symmetric zonal wind, u, defined on a hybrid-coordinate surface.
 * The input parameter gz0 is the value at J = j0.
 * Used to calculate gz_surface for gas giants.
 *
 */
#include <epic_pv_schemes.h>

void gz_from_u(planetspec *planet,
               EPIC_FLOAT *u1d,
               EPIC_FLOAT *p1d,
               EPIC_FLOAT *rho1d,
               EPIC_FLOAT *gz1d,
               int         j0,
               EPIC_FLOAT  gz0)
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
    *gz2d,
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
    dbmsname[]="gz_from_u";

  if (!initialized) {
    /* Allocate memory. */
    udy     = fvector(0,grid.nj+1,dbmsname);
    pvhudy  = fvector(0,grid.nj+1,dbmsname);
    bern    = fvector(0,grid.nj+1,dbmsname);
    sendbuf = fvector(0,grid.nj+1,dbmsname);
    u2d     = fvector(0,Nelem2d-1,dbmsname);
    v2d     = fvector(0,Nelem2d-1,dbmsname);
    gz2d    = fvector(0,Nelem2d-1,dbmsname);

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

  I = ILO;

  /*
   * Calculate absolute vorticity.
   * Store in PV3(nk+1,J,I) so that epic_pv_schemes.h macros (GA_V, DE_V, etc) work.
   */
  K = grid.nk+1;
  vorticity(planet,ON_SIGMATHETA,ABSOLUTE,
            u2d,v2d,NULL,var.pv3.value+(K-Kshift)*Nelem2d);

  /*
   * Calculate udy. 
   */
  for (J = JLOPAD; J <= JHI; J++) {
    udy[J-Jshift] = U2D(J,I)/grid.n[2*J+1];
  }

  /*
   * Calculate local (zeta+f)*u*dy = pvhudy.
   */
  memset(pvhudy,0,(grid.nj+2)*sizeof(EPIC_FLOAT));
  for (J = JFIRST; J <= JHI; J++) {
    /* 
     * Don't multiply by grid.n[2*J] to leave in dy factor:
     */
    pvhudy[J] = (GA_V*udy[J  -Jshift]+DE_V*udy[J  -Jshift]
                +AL_V*udy[J-1-Jshift]+BE_V*udy[J-1-Jshift])*PV_COEF;
    /*
     * Add in pressure-gradient term.
     * P1D and RHO1D are on the P3-grid.
     */
    pvhudy[J] += (P1D(J)-P1D(J-1))*2./(RHO1D(J)+RHO1D(J-1));
  }

#if defined(EPIC_MPI)
  /*
   *  Fill in global-spanning pvuh for each node:
   */ 
  memcpy(sendbuf,pvhudy,(grid.nj+2)*sizeof(EPIC_FLOAT));
  MPI_Allreduce(sendbuf,pvhudy,grid.nj+2,float_type,MPI_SUM,para.comm);
#endif

  /* 
   * Broadcast reference kin. 
   * This is the value at J=j0, I=ILO.
   */
  kin0 = 0.;
  if (j0 >= JLO && j0 <= JHI) {
    kin0 = get_kin(planet,u2d,v2d,j0,ILO);
  }

#if defined(EPIC_MPI)
  sendbuf[0] = kin0;
  MPI_Allreduce(sendbuf,&kin0,1,float_type,MPI_SUM,para.comm);
#endif

  /*
   * Calculate bern by integrating -pvhudy.
   *
   * For convenience, this global-spanning integral is calculated 
   * redundantly on every node.
   */
  bern[j0] = gz0+kin0;
  for (J = j0+1; J <= grid.nj; J++) {
    bern[J] = bern[J-1]-pvhudy[J];
  }
  for (J = j0-1; J >= grid.jlo; J--) {
    bern[J] = bern[J+1]+pvhudy[J+1];
  }

  /*
   * Calculate gz = bern-kin.
   * Store in GZ2D, so that bc_lateral() may be applied.
   */
  for (J = JLO; J <= JHI; J++) {
    kin         = get_kin(planet,u2d,v2d,J,ILO);
    GZ2D(J,ILO) = bern[J]-kin;
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(gz2d,TWODIM);

  for (J = JLOPAD; J <= JHIPAD; J++) {
    GZ1D(J) = GZ2D(J,ILO);
  }

  return;
}
/*======================= end of gz_from_u() ===================================*/

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
   * Check that is is valid.
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

          time_fp_inv  = pressure/var.time_fp_bar;
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
  /* Question: cloud_microphysics(): do we need to separate sync'ing diags from advancing progs? */
  if (grid.microphysics_on) {
    cloud_microphysics(planet,Buff2D);
  }

  /* Question: nmt_apply_sources(): do we need to separate sync'ing diags from advancing progs? */
  if (grid.nmt_physics_on) {
   /* 
    * Apply nmt_physics source terms.
    */
    nmt_apply_sources(planet);
  }

  /*
   * Advance THETA where it is a prognostic variable.
   */
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /*
       * Apply heating term.
       */
      for (K = grid.k_sigma-1; K <= KHI; K++) {
        THETA(K,J,I) += DT*HEAT3(K,J,I)/EXNER3(K,J,I);
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  return;
}

/*======================= end of source_sink() ==================================*/

/*======================= calc_heating() ========================================*/

/*
 * Calculate HEAT3, in W/kg.
 *
 * NOTE: Cloud heating, which contributes to HEAT3, is currently added in elsewhere.
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
   * Zero HEAT3 array.
   */
  memset(var.heat3.value,0,Nelem3d*sizeof(EPIC_FLOAT));

  /*
   * Add Newtonian cooling to HEAT3.
   */
  if (grid.newt_cool_on) {
    newtonian_cooling(planet);
  }

  /*
   * Add solar heating to HEAT3.
   */
  solar_insolation(planet);

  /*
   * Add perturbation heating to HEAT3.
   */
  perturbation_heating(planet);

  /*
   * Add ortho-para hydrogen heating to HEAT3.
   */
  if (var.fpara.on) {
    cprh2 = CPRH2*planet->rgas;
    for (K = KLO; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          fpara       = FPARA(K,J,I);
          pressure    = P3(   K,J,I);
          temperature = T3(   K,J,I);

          return_enthalpy(planet,fpara,pressure,temperature,&fgibb,&fpe,&uoup);

          time_fp_inv  = pressure/var.time_fp_bar;
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
 * NOTE: For spinup experiments, this estimate of dt_cfl will initially be
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
    dy = 1./grid.n[2*J+1];
    for (I = ILO; I <= IHI; I++) {
      for (K = KLO; K <= KHI; K++) {
        /* 
         * Use speed of sound to estimate fastest horizontal gravity-wave speed.
         */
        cs = sqrt(planet->rgas*.5*(T3(K,J,I)+T3(K,J-1,I))/(1.-planet->kappa));

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
      dx = 1./grid.m[2*J+1];
      for (I = ILO; I <= IHI; I++) {
        for (K = KLO; K <= KHI; K++) {
          /* 
           * Use speed of sound to estimate fastest gravity-wave speed.
           */
          cs = sqrt(planet->rgas*.5*(T3(K,J,I)+T3(K,J,I-1))/(1.-planet->kappa));

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

inline int time_mod(int *time,
                    int  step) {
  int
    ans;

  ans = (time[0]%step+(YEAR%step)*(time[1]%step))%step;

  return ans;
}

/*====================== end of time_mod() ===================================*/

/*====================== b_vir() =============================================*/

/*
 * Returns 2nd virial coefficient B_{ab} as a function of temperature.
 * For example, b_vir("H_2","H_2",temperature) returns B for pure molecular
 * hydrogen, H_2, whereas b_vir("H_2","He",temperature) returns the cross-term  
 * B for H_2+He. Data are taken from "The virial coefficients of pure gases
 * and mixtures," by J.H. Dymond and E.B. Smith (1980, Oxford), and converted
 * for a pressure expansion in mks units.
 */

#define MAX_CHEM_PAIRS 8

EPIC_FLOAT b_vir(char       *chem_a,
                 char       *chem_b,
                 EPIC_FLOAT  temperature) 
{
  char
    header[N_STR],
    infile[N_STR],
    chem_ab[16];
  /* Memory for up to MAX_CHEM_PAIRS different chemical pairs: */
  static char
    list[MAX_CHEM_PAIRS][16];
  int
    i,j;
  static int
    ndat[MAX_CHEM_PAIRS],
    count=0;
  EPIC_FLOAT 
    b,t_d,
    *buffer;
  static float_triplet
    *b_table[MAX_CHEM_PAIRS];
  static EPIC_FLOAT
    first_dbdt[MAX_CHEM_PAIRS],
    last_dbdt[ MAX_CHEM_PAIRS];
  FILE
    *input;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="b_vir";

  /* 
   * Determine list index for chemical pair: 
   */
  sprintf(chem_ab,"%s%s",chem_a,chem_b);
  i = 0;
  while (strcmp(chem_ab,list[i]) != 0 && i < count && i < MAX_CHEM_PAIRS) {
    i++;
  }

  if (i >= MAX_CHEM_PAIRS) {
    sprintf(Message,"b_vir() exceeded MAX_CHEM_PAIRS = %d",MAX_CHEM_PAIRS);
    epic_error(dbmsname,Message);
  }
  else if (i == count) {
    /* 
     * Chemical pair not on list. 
     * Add to list and input data. 
     */
    if (IAMNODE == 0) {
      sprintf(list[count],"%s",chem_ab);
      if (strcmp(chem_a,chem_b) == 0) {
        sprintf(infile,EPIC4_PATH"/data/chemistry/virial/b_vs_t.%s",chem_a);
        input = fopen(infile,"r");
        if (!input) {
          fprintf(stderr,"Warning: b_vir() cannot find %s \n",infile);
          fprintf(stderr,"         Defaulting to ideal equation of state for %s. \n",chem_a);
        }
      }
      else {
        sprintf(infile,EPIC4_PATH"/data/chemistry/virial/b_vs_t.%sx%s",chem_a,chem_b);
        input = fopen(infile,"r");
        if (!input) {
          /* Try the names reversed: */
          sprintf(infile,EPIC4_PATH"/data/chemistry/virial/b_vs_t.%sx%s",chem_b,chem_a);
          input = fopen(infile,"r");
          if (!input) {
            fprintf(stderr,"Warning: b_vir() cannot find %s \n",infile);
            fprintf(stderr,"         Defaulting to ideal equation of state for %s.\n",chem_ab);
          }
        }
      }
      if (input) {
        /* Skip over header: */
        for (j = 0; j < 6; j++) {
          fgets(header,100,input);  
        }
        /* Input number of data points: */
        fscanf(input,"%d",ndat+count); 
        /* 
         * Allocate memory: 
         */
        b_table[count] = ftriplet(0,ndat[count]-1,dbmsname);

        /* Input B(T): */
        for (j = 0; j < ndat[count]; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(input,"%lf %lf %*lf",&b_table[count][j].x,&b_table[count][j].y);
#else
          fscanf(input,"%f %f %*f",&b_table[count][j].x,&b_table[count][j].y);
#endif

          /* Convert for pressure expansion in mks units: */
          b_table[count][j].y /= 1.e+3*R_GAS*b_table[count][j].x;
        }
        fclose(input);
      }
      else {
        /*
         * When a data file is not available, default to ideal equation of state
         * by setting bdat to zero:
         */
        ndat[count] = 3;
        /*
         * Allocate memory.
         */
        b_table[count] = ftriplet(0,ndat[count]-1,dbmsname);
        for (j = 0; j < ndat[count]; j++) {
          b_table[count][j].x = (EPIC_FLOAT)(j+1)*100.;
          b_table[count][j].y = 0.;
        }
      }
      count++;
    }

#if defined(EPIC_MPI)
    MPI_Bcast(&count,                   1,MPI_INT,   NODE0,para.comm);
    MPI_Bcast(list[count-1],            8,MPI_CHAR,  NODE0,para.comm);
    MPI_Bcast(ndat+(count-1),           1,MPI_INT,   NODE0,para.comm);
    /* pack buffer */
    buffer = fvector(0,2*ndat[count-1]-1,dbmsname);
    for (j = 0; j < ndat[count-1]; j++) {
      buffer[j              ] = b_table[count-1][j].x;
      buffer[j+ndat[count-1]] = b_table[count-1][j].y;
    }
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Bcast(buffer,2*ndat[count-1],MPI_DOUBLE,NODE0,para.comm);
#  else
     MPI_Bcast(buffer,2*ndat[count-1],MPI_FLOAT,NODE0,para.comm);
#  endif
    /* unpack buffer */
    for (j = 0; j < ndat[count-1]; j++) {
      b_table[count-1][j].x = buffer[j              ];
      b_table[count-1][j].y = buffer[j+ndat[count-1]];
    }
    free_fvector(buffer,0,2*ndat[count-1]-1,dbmsname);
#endif

    /* Set endpoint slopes: */
    first_dbdt[count-1] = 
          (b_table[count-1][1].y-b_table[count-1][0].y)/
          (b_table[count-1][1].x-b_table[count-1][0].x);
    last_dbdt[count-1]  = 
          (b_table[count-1][ndat[count-1]-1].y-b_table[count-1][ndat[count-1]-2].y)/
          (b_table[count-1][ndat[count-1]-1].x-b_table[count-1][ndat[count-1]-2].x);
    /* Prepare for monotonic spline interpolation: */
    spline_pchip(ndat[count-1],b_table[count-1]);
  }

  /* 
   * Main function evaluation:
   */

  /* Use cubic-spline interpolation: */
  if (temperature <= b_table[i][0].x) {
    /* At or before start of table. */
    b = b_table[i][0].y+first_dbdt[i]*(temperature-b_table[i][0].x);
  }
  else if (temperature >= b_table[i][ndat[i]-1].x) {
    /* At or after end of table. */
    b = b_table[i][ndat[i]-1].y+last_dbdt[i]*(temperature-b_table[i][ndat[i]-1].x);
  }
  else {
    j = find_place_in_table(ndat[i],b_table[i],temperature,&t_d);
    b = splint_pchip(temperature,b_table[i]+j,t_d);
  }

  return b;
}

/*====================== end of b_vir() ======================================*/

/*====================== b1_vir() ============================================*/

/*
 * Returns B1 = T dB/dT.
 */
EPIC_FLOAT b1_vir(char       *chem_a,
                  char       *chem_b,
                  EPIC_FLOAT  temperature)
{
  EPIC_FLOAT
    b1,tt;
  static EPIC_FLOAT
    dt=1.;

  tt = temperature/dt;
  b1 = (b_vir(chem_a,chem_b,temperature+.5*dt)-
        b_vir(chem_a,chem_b,temperature-.5*dt))*tt;

  return b1;
}

/*====================== end of b1_vir() =====================================*/

/*====================== b2_vir() ============================================*/

/*
 * Returns B2 = T^2 (d/dT)^2 B.
 */
EPIC_FLOAT b2_vir(char       *chem_a,
                  char       *chem_b,
                  EPIC_FLOAT  temperature)
{
  EPIC_FLOAT
    b2,tt;
  static EPIC_FLOAT
    dt=1.;

  tt = temperature/dt;
  b2 = (    b_vir(chem_a,chem_b,temperature+dt)
        -2.*b_vir(chem_a,chem_b,temperature   )
           +b_vir(chem_a,chem_b,temperature-dt))*tt*tt;

  return b2;
}

/*====================== end of b2_vir() =====================================*/

/*====================== sum_xx() ============================================*/
/*
 * Returns sum of 2nd virial coefficient, or related function,
 * with quadratic mole-fraction weighting appropriate to specified planet.
 */
EPIC_FLOAT sum_xx(planetspec *planet,
                  EPIC_FLOAT (*b_func)(char *,
                               char *,
                               EPIC_FLOAT),
                  EPIC_FLOAT   temperature)
{
  static int
    initialized = FALSE;
  EPIC_FLOAT
    b_sum,x_sum;
  static EPIC_FLOAT
    x_H_2,x_He;

  if (strcmp(grid.eos,"ideal") == 0) {
    /*
     * B = 0 for ideal case:
     */
    b_sum = 0.;
  }
  else if (strcmp(grid.eos,"virial") == 0) {
    if (strcmp(planet->type,"gas-giant") == 0) {
      if (!initialized) {
        /* NOTE: currently including only H_2+He: */
        x_sum  = planet->x_h2+planet->x_he;
        x_H_2  = planet->x_h2/x_sum;
        x_He   = planet->x_he/x_sum;

        initialized = TRUE;
      }
      b_sum =    (*b_func)("H_2","H_2",temperature)*x_H_2*x_H_2
             +2.*(*b_func)("H_2","He", temperature)*x_H_2*x_He
                +(*b_func)("He", "He", temperature)*x_He *x_He;
    }
    else {
      switch(planet->index) {
        case VENUS_INDEX:
        case VENUS_LLR05_INDEX:
        case MARS_INDEX:
          b_sum = (*b_func)("CO_2","CO_2",temperature);
        break;
        case EARTH_INDEX:
        case TITAN_INDEX:
          b_sum = (*b_func)("N_2","N_2",temperature);
        break;
        default:
          /* Default to ideal case. */
          if (!initialized) {
            if (IAMNODE == 0) {
              fprintf(stderr,"Warning: sum_xx(): equation of state = %s not defined for %s, "
                             "using ideal equation of state.\n",grid.eos,planet->name);
            }
            initialized = TRUE;
          }
          b_sum = 0.;
        break;
      }
    }
  }
  else {
    fprintf(stderr,"Unrecognized equation of state = %s in sum_xx()\n",grid.eos);
    exit(1);
  }

  return b_sum;
}

/*====================== end of sum_xx() =====================================*/

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
  EPIC_FLOAT
    num,denom,
    mu_avg,mu_dry;
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
    return .5*(avg_molar_mass(planet,kk-1,J,I)+
               avg_molar_mass(planet,kk+1,J,I));
  }
  else {
    /*
     * Interface value.
     */
    K = (kk-1)/2;
    num   = 1.;
    denom = planet->rgas/R_GAS;
    for (iq = 0; iq < grid.nq; iq++) {
      num   += Q(grid.is[iq],grid.ip[iq],K,J,I);
      denom += Q(grid.is[iq],grid.ip[iq],K,J,I)/var.species[grid.is[iq]].molar_mass;
    }

    return num/denom;
  }
}

/*====================== end of avg_molar_mass() =============================*/

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

  if (index == HDRY_INDEX || index == HDRY3_INDEX) {
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
    sumv[LAST_SPECIES+1];
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
     * For molecules not listed, we sum the structural diffusion-volumes (for ozone,
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
    sumv[C_2H_6_INDEX] = 16.5*2.+1.98*6.;
    sumv[CO_2_INDEX]   = 26.9;
    sumv[NH_4SH_INDEX] = 5.69+1.98*4.+17.0+1.98;
    sumv[O_3_INDEX]    = 16.6*3./2.;
    sumv[N_2_INDEX]    = 17.9;
    /*
     * For planets where sumv[HDRY_INDEX] has not
     * been measured, assume mole-fraction weighting.
     *
     * NOTE: Many of the mole fractions are rough and should be 
     *       updated when possible. 
     */
    if (strcmp(planet->type,"gas-giant") == 0) {
      sumv[HDRY_INDEX] = planet->x_h2*sumv_h2+planet->x_he*sumv_he;
    }
    else if (strcmp(planet->type,"terrestrial") == 0) {
      switch(planet->index) {
        case VENUS_INDEX:
        case VENUS_LLR05_INDEX:
        case MARS_INDEX:
          sumv[HDRY_INDEX] = sumv[CO_2_INDEX];
        break;
        case EARTH_INDEX:
        case HELD_SUAREZ_INDEX:
          sumv[HDRY_INDEX] = 20.1;
        break;
        case TITAN_INDEX:
          sumv[HDRY_INDEX] = sumv[N_2_INDEX];
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
    case HDRY_INDEX:
    case H_2O_INDEX:
    case NH_3_INDEX:
    case H_2S_INDEX:
    case CH_4_INDEX:
    case C_2H_2_INDEX:
    case C_2H_6_INDEX:
    case CO_2_INDEX:
    case NH_4SH_INDEX:
    case O_3_INDEX:
    case N_2_INDEX:
      tmp = pow(sumv[HDRY_INDEX],one_third)+pow(sumv[vapor_index],one_third);
    break;
    case HDRY3_INDEX:
      tmp = pow(sumv[HDRY_INDEX],one_third)+pow(sumv[HDRY_INDEX],one_third);
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

/*====================== parse_species_name() ================================*/

/*
 * Takes a species name string (e.g., "NH_4SH") and returns  
 * the number of distinct elements, their symbols, and how many atoms 
 * of each are present.  Reallocates the necessary memory to hold this
 * information.
 *
 * NOTE: num_elements and **symbols shoud be declared static in the calling 
 *       function, with their input values equal to the last call, 
 *       in order properly reallocate memory.
 *
 */
void parse_species_name(char   *name,
                        int    *num_elements,
                        char ***symbols,
                        int   **counts)
{
  int
    num_caps,
    i,ii;
  char
    *ptr,
    subscript[4],
    format[4];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="parse_species_name";

  /*
   * Free previous memory:
   */
  for (i = 0; i < *num_elements; i++) {
    free((*symbols)[i]);
  }

  /*
   * Count number of capital letters to determine working array sizes:
   */
  num_caps = 0;
  ptr      = name;
  while (*ptr != '\0') {
    if (isupper(*ptr)) num_caps++;
    ptr++;
  }

  /*
   * Reallocate memory:
   */
  *counts  = (int   *)realloc(*counts, num_caps*sizeof(int   ));
  *symbols = (char **)realloc(*symbols,num_caps*sizeof(char *));
  for (i = 0; i < num_caps; i++) {
    (*symbols)[i] = (char *)calloc(4,sizeof(char));
  }

  /*
   * Determine symbols and counts:
   */
  i   = 0;
  ptr = name;
  while (*ptr != '\0') {
    if (isupper(*ptr)) {
      if (islower(*(ptr+1))) {
        /* Element symbol has two letters: */
        strncpy((*symbols)[i],ptr,2);
        ptr+=2;
      }
      else {
        /* Element symbol has one letter: */
        strncpy((*symbols)[i],ptr,1);
        ptr+=1;
      }
      if (*ptr == '_') {
        subscript[0] = '\0';
        /* Determine subscript's number of places: */
        ii = 0;
        while(isdigit(*(++ptr))) {
          ii++;
        }
        if (ii == 0) {
          sprintf(Message," \"_\" not followed by digit: %s",name);
          epic_error(dbmsname,Message);
        }
        else {
          sprintf(format,"%%%dd",ii);
          sscanf(ptr-ii,format,*counts+i);
        }
      }
      else {
        (*counts)[i] = 1;
      }
      i++;
    }
    else {
      ptr++;
    }
  }
  /*
   * Trim arrays to refer to distinct elements:
   */
  *num_elements = num_caps;
  for (i = 0; i < num_caps; i++) {
    for (ii = i+1; ii < num_caps; ii++) {
      if ((*counts)[ii] > 0 && strcmp((*symbols)[i],(*symbols)[ii]) == 0) {
        (*counts)[i ] += (*counts)[ii];
        (*counts)[ii]  = 0;
        (*num_elements)--;
      }
    }
  }
  /* Remove zero entries by shifting: */
  for (i = 0; i < num_caps; i++) {
    if ((*counts)[i] == 0) {
      for (ii = i; ii < num_caps-1; ii++) {
        (*counts)[ii] = (*counts)[ii+1];
        strcpy((*symbols)[ii],(*symbols)[ii+1]);
      }
      (*counts)[num_caps-1] = 0;
    }
  }
  /* Trim allocated memory */
  for (i = (*num_elements); i < num_caps; i++) {
    free((*symbols)[i]);
  }
  *counts  = (int   *)realloc(*counts, (*num_elements)*sizeof(int   ));
  *symbols = (char **)realloc(*symbols,(*num_elements)*sizeof(char *));

  return;
}

/*====================== end of parse_species_name() =========================*/

/*====================== solar_fraction() ====================================*/

/*
 * Returns solar mixing ratio of the least-abundant element in
 * the given species name (divided by its stochiometric count).
 * Choices for the type argument: MASS, MOLAR. 
 * The character string min_element should be 4 bytes.
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
  min_count     =  1;
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

/*====================== end of solar_fraction() =============================*/
  
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
  p0 = 1.e+5;

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

  if (planet->x_h2 == 0.) {
    /*
     * No hydrogen, so use the standard relationship between
     * pressure, temperature, and theta.
     */
    temperature = theta*pow(p/grid.press0,planet->kappa);
  }
  else {

#if EPIC_CHECK == 1
    /* Sanity checks: */
    if (fp < 0. || fp > 1.) {
      sprintf(Message,"fp=%e",fp);
      epic_error(dbmsname,Message);
    }
#endif

    theta1 = theta*pow(p/grid.press0,planet->kappa);
  
    if (theta1 <= THLO_THERMO) {
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

  sprintf(Message,"should never get here");
  epic_error(dbmsname,Message);
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
  EPIC_FLOAT
    b,b1,kappa,tmp,
    theta,thetaln,
    cc,tt,pp,
    em,fract;
  EPIC_FLOAT
    thermo_vector[2];

  if (planet->x_h2 == 0.) {
    /*
     * No hydrogen, so use standard definition of theta.
     */
    theta = temperature*pow(grid.press0/p,planet->kappa);
  }
  else {
    /*
     * Use mean theta as defined in Dowling et al (1998), to handle ortho/para hydrogen.
     */
    if (temperature <= 20.) {
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
    pp            = pow(grid.press0/p,planet->kappa);
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

  if (planet->x_h2 == 0.) {
    /*
     * No hydrogen, so use the standard relationship between 
     * pressure, temperature, and theta.
     */
    press = grid.press0*pow(temperature/theta,planet->cpr);
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
    press = grid.press0*pow(temperature/theta,planet->cpr);
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
    if (temperature <= 20.) {
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

/*======================= enthalpy() ========================================*/

/* NOTE: This function needs to be tested before it is put into use. */

/*
 * Returns absolute enthalpy, defined as
 *   ...eqns, units [J/kg]
 *  See the book "The Properties of geses and Liquids", Fifth Edition.Bruce E. Poling, John M. Prausnitz, John P. O'Connell, 2000.
 *  Appendix, Section C, page A.35-46
 *
 * Nino Ghurtskaia, Feb 2007.
 */

 /* The absolute enthalpy of substance H(T) is defined in terms of its formation enthalpy and its heat content as follows:
    
	     H_T = DHf(298K) + [HT-H298] */

/* A common six-term equation for isobaric heat content:

 HT - H298 = A1*T + B1*T*T + C1*(1./T) + D*sqrt(T) + E*T*T*T + F  	

      Database from FREED software for:
			      
      CH4    Carbon Tetrahydride     - Methane
      C2H6   Carbon Trihydride dimer  - Ethane   
      CO2    Carbon Dioxide 
      H2S    Hydrogen Sulfide     
      H2SO4  Hydrogen Sulfate        - Sulfuric Acid
      NH3    Nitrogen Trihydride     - Ammonia
      N2     Nitrogen dimer         
      O3     Oxygen Trimer           - Ozone

*/
 
  
EPIC_FLOAT enthalpy(EPIC_FLOAT T,
                    int        INDEX)
{
  
  double
    ans=0;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="enthalpy";

  switch (INDEX) {
    
    case H_2O_INDEX:  
      
	
	if ((50. <= T) && (T <= 1000.)) {
	  const double a0=4.395, a1=-4.186e-3, a2=1.405e-5, a3=-1.564e-8, a4=0.632e-11, C= -5.7717e+05; // C is an integration constant
	               ans = ((a0 + (a1/2 + (a2/3 + (a3/4 + (a4*T)/5)*T)*T)*T)*T)*R_GAS/18 + C;		       
	  }	 
	  	
	//else if ((3000. < T) && (T <= 5000.)) {
	  //const double A2=13.92501, B2=6.35227e-05, C2=8672045, D2=.0, E2=.0, F2=-15036.0;
            //           ans = T*(T*((T*E2) + B2) + A2) + C2*(1./T) + D2*sqrt(T) + F2 - 57795.0;	
	//}
	
	else {
	  
	    if (T < 50.) sprintf(Message,"Error, temperature should be more than 50");
	    else
	    {
	      sprintf(Message,"T=%f exceeds limit=%f",T,1000.);
	     
	  }
	 epic_error(dbmsname,Message);
	}			
    
    break;

    case NH_3_INDEX:
      		 	    
	  	
	if ((50. <= T ) && (T <= 1000.)) {
	  const double a0=4.238, a1=-4.215e+3, a2=2.041e+5, a3=-2.126e+8, a4=0.761e+11, C=-5.9836e+05;	  	  
	               ans = ((a0 + (a1/2 + (a2/3 + (a3/4 + (a4*T)/5)*T)*T)*T)*T)*R_GAS/17 + C;
		       
	  }
	  
		
	else {
	  sprintf(Message,"T=%f exceeds limit=%f",T,1000.);
	  epic_error(dbmsname,Message);
	}	  
   
    break;
	
    case H_2S_INDEX:
	    	    
	
	if ((50 <= T) && (T <= 1000.)) {
	  const double a0=4.266, a1=-3.438e+3, a2=1.319e+5, a3=-1.331e+8, a4=0.488e+11, C=-2.9685e+05;	  	  
	               ans = ((a0 + (a1/2 + (a2/3 + (a3/4 + (a4*T)/5)*T)*T)*T)*T)*R_GAS/34 + C;
	              
	  }
	
	         
	else {
	  sprintf(Message,"T=%f exceeds limit=%f",T,1000.);
	  epic_error(dbmsname,Message);
	}
   
    break;
	 
	 	 
   // case H_2SO_4_INDEX:
	
	    		
   //   if ((50. <= T) && (T <= 1000.)) {
	//  const double A1= 52.137 , B1=-0.00142266, C1=-178200, D1=-1150.52, E1=1.51184e-07, F1=5041.0;
	//	       ans = T*(T*(T*E1 + B1) + A1) + C1*(1./T) + D1*sqrt(T) + F1 - 175700.0;	
	//
	//  }
	     
	 	  
	//else {
	//  sprintf(Message,"T=%f exceeds limit=%f",T,1000.);
	//  epic_error(dbmsname,Message);
	//}
	
// break;
	 	 
    case CH_4_INDEX:
			
				   
	if ((50 <= T) && (T <= 1000.)) {
	  const double a0=4.568, a1=-8.975e+3, a2=3.631e+5, a3=-3.407e+8, a4=1.091e+11, C=-6.3945e+05;	  	  
	               ans = ((a0 + (a1/2 + (a2/3 + (a3/4 + (a4*T)/5)*T)*T)*T)*T)*R_GAS/16+ C; 
		      
	
	  }
					  
	//else if((1500. < T) && (T <= 4500.)) {
	 // const double A2=32.64010799, B2=-0.000215464, C2=4381338.687, D2=-726.1264935,  F2= -4593.0;	  
		// ans  =  T*(B2*T +A2) + C2*(1./T) + D2*sqrt(T) + F2 - 17880.0 ; // E2=0.
	
	  //}
  	
	else {
	  sprintf(Message,"T=%f exceeds limit=%f",T,1000.);
	  epic_error(dbmsname,Message);
	}
	
    break;

   // case C_2H_2_INDEX:
	
			  	  
	// if ((50. <= T) && (T <= 1000.)) {
	  // need this function 
		      
	  // }
	 
	
       // else {
	  // sprintf(Message,"T=%f exceeds limit=%f",T,1000.);
	  // epic_error(dbmsname,Message);
	// }
	
    // break;

      
   // case C_2H_6_INDEX:
	
			  	  
	// if ((50. <= T) && (T <= 1000.)) {
	  // const double a0=4.178, a1=-4.427e+3, a2=5.66e+5, a3=-6.651e+8, a4=2.487e+11; C=-6.6294e+05;
	 	      // ans = ((a0 + (a1/2 + (a2/3 + (a3/4 + (a4*T)/5)*T)*T)*T)*T)*R_GAS/30+ C; 
		      
	  // }
	 
	
       // else {
	  // sprintf(Message,"T=%f exceeds limit=%f",T,1000.);
	  // epic_error(dbmsname,Message);
	// }
	
    // break;
	
    case CO_2_INDEX:
    
	  	  
        if ((50 <= T) && (T <= 1000.)) {
	  const double a0=3.259, a1=1.356e+3, a2=1.502e+5, a3=-2.374e+8, a4=1.056e+11, C=-2.2115e+05;	  	  
	               ans = ((a0 + (a1/2 + (a2/3 + (a3/4 + (a4*T)/5)*T)*T)*T)*T)*R_GAS/44 + C;
		      
	
	  }
	  
   	else {
	  sprintf(Message,"T=%f exceeds limit=%f",T,1000.);
	  epic_error(dbmsname,Message);
	}
	
	 break;
	 
     // case NH_4SH_INDEX:                       
	  
	//  double a0=, a1=, a2=, a3=, a4=, a5=;
	  	  
	  //if ((298.15 <= T) && (T <= 1500)) {
	  
	  	//	ans = ((a0 + (a1/2 + (a2/3 + (a3/4 + (a4*T)/5)*T)*T)*T)*T)*R_GAS/35 + C;	
	
	  //}
	  
						  
	  //else ((1500 < T) && (T <= 4500)) {
	       //  double a0=, a1=, a2=, a3=, a4=, a5=;
	      //   ans = (a0 + (a1/2 + (a2/3 + (a3/4 + (a4*T)/5)*T)*T)*T)*T + C;	
	
	  //} 
	  //else {
	  //sprintf(Message,"T=%f exceeds limit=%f",T,45000.);
	  //epic_error(dbmsname,Message);
	//} 
	
    //break;
    	
    case N_2_INDEX:
      
        	  
	  if ((50. <= T) && (T <= 1000.)) {
	    const double a0=3.539, a1=-0.261e+3, a2=0.007e+5, a3=-0.157e+8, a4=-0.099e+11, C=-3.0410e+05;	  	  
	               ans = ((a0 + (a1/2 + (a2/3 + (a3/4 + (a4*T)/5)*T)*T)*T)*T)*R_GAS/28 + C;
		      
	
	  }
	  
	   				  		  
	//  else if((3000. < T) && (T <= 5000.)) {
	//    const double A2=8.511, B2=0.0000565, F2=-3883.0;   
	//	   ans  =  T*(B2*T +A2) + F2;	//  C2=0., D2=0., E2=0. 
	
	//  }
	  
	  else {
	  sprintf(Message,"T=%f exceeds limit=%f",T,1000.);
	  epic_error(dbmsname,Message);
	} 
	
    break;
    
 case O_3_INDEX:
      
        	  
	  if ((50. <= T) && (T <= 1000.)) {
	    const double a0=4.106, a1=-3.809e+3, a2=3.131e+5, a3=-4.3e+8, a4=1.813e+11, C=-2.0976e+05;	  	  
	               ans = ((a0 + (a1/2 + (a2/3 + (a3/4 + (a4*T)/5)*T)*T)*T)*T)*R_GAS/48 + C;
		      
	
	  }
	 	  
	  else {
	  sprintf(Message,"T=%f exceeds limit=%f",T,1000.);
	  epic_error(dbmsname,Message);
	} 
	
    break;

	
    default:
      sprintf(Message,"index=%d not yet defined",INDEX);
      epic_error(dbmsname, Message);
    break;
  }
  
  return ans;
}


/*======================= end of enthalpy() =================================*/

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

/*====================== check_periodic() ====================================*/

/*
 * Debugging tool to check whether periodicity, bc_lateral(), is being applied
 * properly.
 *
 * NOTE: Not ready for MPI decomposition in I direction.
 */

#define CHECK_PERIODIC_WIND(U) \
  periodic = TRUE; \
  for (K = KLOPAD; K <= KHIPAD; K++) { \
    for (J = JLOPAD; J <= JHIPADPV; J++) { \
      if (U(grid.it_uv,K,J,ILO-1) != U(grid.it_uv,K,J,IHI)) { \
        fprintf(stderr,"%s(%2d,%2d,%2d,%2d)=%e != %s(%2d,%2d,%2d,%2d)=%e\n", \
                       #U,grid.it_uv,K,J,ILO-1,U(grid.it_uv,K,J,ILO-1), \
                       #U,grid.it_uv,K,J,IHI,U(grid.it_uv,K,J,IHI));fflush(stderr); \
        periodic = FALSE; \
      } \
      if (U(grid.it_uv,K,J,IHI+1) != U(grid.it_uv,K,J,ILO)) { \
        fprintf(stderr,"%s(%2d,%2d,%2d,%2d)=%e != %s(%2d,%2d,%2d,%2d)=%e\n", \
                       #U,grid.it_uv,K,J,IHI+1,U(grid.it_uv,K,J,IHI+1), \
                       #U,grid.it_uv,K,J,ILO,U(grid.it_uv,K,J,ILO));fflush(stderr); \
        periodic = FALSE; \
      } \
    } \
  } \
  if (periodic) { \
    fprintf(stderr,"%s ",#U);fflush(stderr); \
  } \
  else { \
    sprintf(Message,"\n itime=%d",grid.itime); \
    epic_error(message,Message); \
  }

#define CHECK_PERIODIC_H(HDRY) \
  periodic = TRUE; \
  for (K = KLOPAD; K <= KHIPAD; K++) { \
    for (J = JLOPAD; J <= JHIPADPV; J++) { \
      if (HDRY(K,J,ILO-1) != HDRY(K,J,IHI)) { \
        fprintf(stderr,"%s(%2d,%2d,%2d)=%e != %s(%2d,%2d,%2d)=%e ", \
                       #HDRY,K,J,ILO-1,HDRY(K,J,ILO-1),#HDRY,K,J,IHI,HDRY(K,J,IHI));fflush(stderr); \
        periodic = FALSE; \
      } \
      if (HDRY(K,J,IHI+1) != HDRY(K,J,ILO)) { \
        fprintf(stderr,"%s(%2d,%2d,%2d)=%e != %s(%2d,%2d,%2d)=%e ", \
                       #HDRY,K,J,IHI+1,HDRY(K,J,IHI+1),#HDRY,K,J,ILO,HDRY(K,J,ILO));fflush(stderr); \
        periodic = FALSE; \
      } \
    } \
  } \
  if (periodic) { \
    fprintf(stderr,"%s ",#HDRY);fflush(stderr); \
  } \
  else { \
    sprintf(Message,"\n itime=%d",grid.itime); \
    epic_error(message,Message); \
  }

#define CHECK_PERIODIC_DIAG(T2) \
  periodic = TRUE; \
  for (K = KLOPAD; K <= KHIPAD; K++) { \
    for (J = JLOPAD; J <= JHIPADPV; J++) { \
      if (T2(K,J,ILO-1) != T2(K,J,IHI)) { \
        fprintf(stderr,"%s(%2d,%2d,%2d)=%e != %s(%2d,%2d,%2d)=%e ", \
                       #T2,K,J,ILO-1,T2(K,J,ILO-1),#T2,K,J,IHI,T2(K,J,IHI));fflush(stderr); \
        periodic = FALSE; \
      } \
      if (T2(K,J,IHI+1) != T2(K,J,ILO)) { \
        fprintf(stderr,"%s(%2d,%2d,%2d)=%e != %s(%2d,%2d,%2d)=%e ", \
                       #T2,K,J,IHI+1,T2(K,J,IHI+1),#T2,K,J,ILO,T2(K,J,ILO));fflush(stderr); \
        periodic = FALSE; \
      } \
    } \
  } \
  if (periodic) { \
    fprintf(stderr,"%s ",#T2);fflush(stderr); \
  } \
  else { \
    sprintf(Message,"\n itime=%d",grid.itime); \
    epic_error(message,Message); \
  }

#define CHECK_PERIODIC_SPECIES(is) \
  for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) { \
    if (var.species[is].phase[ip].on) { \
      periodic = TRUE; \
      for (K = KLOPAD; K <= KHIPAD; K++) { \
        for (J = JLOPAD; J <= JHIPADPV; J++) { \
          if (Q(is,ip,K,J,ILO-1) != Q(is,ip,K,J,IHI)) { \
            fprintf(stderr,"%s(%2d,%2d,%2d)=%e != %s(%2d,%2d,%2d)=%e ", \
                           var.species[is].phase[ip].info[0].name,K,J,ILO-1,Q(is,ip,K,J,ILO-1), \
                           var.species[is].phase[ip].info[0].name,K,J,IHI,  Q(is,ip,K,J,IHI  ));fflush(stderr); \
            periodic = FALSE; \
          } \
          if (Q(is,ip,K,J,IHI+1) != Q(is,ip,K,J,ILO)) { \
            fprintf(stderr,"%s(%2d,%2d,%2d)=%e != %s(%2d,%2d,%2d)=%e ", \
                           var.species[is].phase[ip].info[0].name,K,J,IHI+1,Q(is,ip,K,J,IHI+1), \
                           var.species[is].phase[ip].info[0].name,K,J,ILO,  Q(is,ip,K,J,ILO  ));fflush(stderr); \
            periodic = FALSE; \
          } \
        } \
      } \
      if (periodic) { \
        fprintf(stderr,"%s ",var.species[is].phase[ip].info[0].name);fflush(stderr); \
      } \
      else { \
        sprintf(Message,"\n itime=%d",grid.itime); \
        epic_error(message,Message); \
      } \
    } \
  } 

void check_periodic(char *message)
{
  int
    K,J,
    is,ip,im,
    ip_first,ip_last,
    periodic;

  if (IAMNODE == NODE0) {
    fprintf(stderr,"IT=%lu %s:periodic:",grid.itime,message);fflush(stderr);
  }

  /*
   * Check prognostic variables.
   */
  if (var.u.on)     {CHECK_PERIODIC_WIND(U)};
  if (var.v.on)     {CHECK_PERIODIC_WIND(V)};
  if (var.hdry.on)  {CHECK_PERIODIC_H(HDRY)};
  if (var.theta.on) {CHECK_PERIODIC_H(THETA)};
  if (var.nu_turb.on) {CHECK_PERIODIC_H(NU_TURB)};
  if (var.fpara.on) {CHECK_PERIODIC_H(FPARA)};
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {CHECK_PERIODIC_SPECIES(is)};
  }

  /*
   * Check 3D diagnostic variables.
   */
  if (var.hdry3.on)    {CHECK_PERIODIC_DIAG(HDRY3)};
  if (var.pdry3.on)    {CHECK_PERIODIC_DIAG(PDRY3)};
  if (var.p2.on)       {CHECK_PERIODIC_DIAG(P2)};
  if (var.p3.on)       {CHECK_PERIODIC_DIAG(P3)};
  if (var.theta2.on)   {CHECK_PERIODIC_DIAG(THETA2)};
  if (var.h2.on)       {CHECK_PERIODIC_DIAG(H2)};
  if (var.h3.on)       {CHECK_PERIODIC_DIAG(H3)};
  if (var.t2.on)       {CHECK_PERIODIC_DIAG(T2)};
  if (var.t3.on)       {CHECK_PERIODIC_DIAG(T3)};
  if (var.rho2.on)     {CHECK_PERIODIC_DIAG(RHO2)};
  if (var.rho3.on)     {CHECK_PERIODIC_DIAG(RHO3)};
  if (var.exner3.on)   {CHECK_PERIODIC_DIAG(EXNER3)};
  if (var.fgibb3.on)   {CHECK_PERIODIC_DIAG(FGIBB3)};
  if (var.gz3.on)      {CHECK_PERIODIC_DIAG(GZ3)};
  if (var.mont3.on)    {CHECK_PERIODIC_DIAG(MONT3)};
  if (var.heat3.on)    {CHECK_PERIODIC_DIAG(HEAT3)};
  if (var.pv3.on)      {CHECK_PERIODIC_DIAG(PV3)};
  if (var.ri2.on)      {CHECK_PERIODIC_DIAG(RI2)};
  if (var.div_uv3.on)  {CHECK_PERIODIC_DIAG(DIV_UV3)};
  if (var.w3.on)       {CHECK_PERIODIC_DIAG(W3)};
  if (var.dzdt3.on)    {CHECK_PERIODIC_DIAG(DZDT3)};

  if (IAMNODE == NODE0) {
    fprintf(stderr,"\n");fflush(stderr);
  }

  return;
}

/*====================== end of check_periodic() =============================*/

/*====================== check_nan() =========================================*/

/*
 * Debugging tool to screen for occurrances of nan (not-a-number).
 */

#define CHECK_NAN_WIND(U) \
  no_nan = TRUE; \
  for (K = KLOPAD; K <= KHIPAD; K++) { \
    for (J = JLOPAD; J <= JHIPADPV; J++) { \
      for (I = ILOPAD; I <= IHIPAD; I++) { \
        if (!isfinite(U(grid.it_uv,K,J,I))) { \
          fprintf(stderr,"%s(%2d,%2d,%2d,%2d)=%g ", \
                          #U,grid.it_uv,K,J,I,U(grid.it_uv,K,J,I));fflush(stderr); \
          no_nan = FALSE; \
        } \
      } \
    } \
  } \
  if (no_nan) { \
    fprintf(stderr,"%s ",#U);fflush(stderr); \
  } \
  else { \
    sprintf(Message,"\n itime=%d",grid.itime); \
    epic_error(message,Message); \
  }


#define CHECK_NAN_H(HDRY) \
  no_nan = TRUE; \
  for (K = KLOPAD; K <= KHIPAD; K++) { \
    for (J = JLOPAD; J <= JHIPADPV; J++) { \
      for (I = ILOPAD; I <= IHIPAD; I++) { \
        if (!isfinite(HDRY(K,J,I))) { \
          fprintf(stderr,"%s(%2d,%2d,%2d)=%g ", \
                       #HDRY,K,J,I,HDRY(K,J,I));fflush(stderr); \
          no_nan = FALSE; \
        } \
      } \
    } \
  } \
  if (no_nan) { \
    fprintf(stderr,"%s ",#HDRY);fflush(stderr); \
  } \
  else { \
    sprintf(Message,"\n itime=%d",grid.itime); \
    epic_error(message,Message); \
  }

#define CHECK_NAN_DIAG(T2) \
  no_nan = TRUE; \
  for (K = KLOPAD; K <= KHIPAD; K++) { \
    for (J = JLOPAD; J <= JHIPADPV; J++) { \
      for (I = ILOPAD; I <= IHIPAD; I++) { \
        if (!isfinite(T2(K,J,I))) { \
          fprintf(stderr,"%s(%2d,%2d,%2d)=%g ", \
                       #T2,K,J,I,T2(K,J,I));fflush(stderr); \
          no_nan = FALSE; \
        } \
      } \
    } \
  } \
  if (no_nan) { \
    fprintf(stderr,"%s ",#T2);fflush(stderr); \
  } \
  else { \
    sprintf(Message,"\n itime=%d",grid.itime); \
    epic_error(message,Message); \
  }

#define CHECK_NAN_SPECIES(is) \
  for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) { \
    if (var.species[is].phase[ip].on) { \
      no_nan = TRUE; \
      for (K = KLOPAD; K <= KHIPAD; K++) { \
        for (J = JLOPAD; J <= JHIPADPV; J++) { \
          for (I = ILOPAD; I <= IHIPAD; I++) { \
            if (!isfinite(Q(is,ip,K,J,I))) { \
              fprintf(stderr,"%s(%2d,%2d,%2d)=%g ", \
                             var.species[is].phase[ip].info[0].name,K,J,I,Q(is,ip,K,J,I));fflush(stderr); \
              no_nan = FALSE; \
            } \
          } \
        } \
      } \
      if (no_nan) { \
        fprintf(stderr,"%s ",var.species[is].phase[ip].info[0].name);fflush(stderr); \
      } \
      else { \
        sprintf(Message,"\n itime=%d",grid.itime); \
        epic_error(message,Message); \
      } \
    } \
  } 

void check_nan(char *message)
{
  int
    K,J,I,
    is,ip,im,
    ip_first,ip_last,
    no_nan;

  if (IAMNODE == NODE0) {
    fprintf(stderr,"IT=%lu %s: no nan:",grid.itime,message);fflush(stderr);
  }

  /*
   * Check prognostic variables.
   */
  if (var.u.on)     {CHECK_NAN_WIND(U)};
  if (var.v.on)     {CHECK_NAN_WIND(V)};
  if (var.hdry.on)  {CHECK_NAN_H(HDRY)};
  if (var.theta.on) {CHECK_NAN_H(THETA)};
  if (var.nu_turb.on) {CHECK_NAN_H(NU_TURB)};
  if (var.fpara.on) {CHECK_NAN_H(FPARA)};
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {CHECK_NAN_SPECIES(is)};
  }

  /*
   * Check 3D diagnostic variables.
   */
  if (var.hdry3.on)      {CHECK_NAN_DIAG(HDRY3)};
  if (var.pdry3.on)      {CHECK_NAN_DIAG(PDRY3)};
  if (var.p2.on)         {CHECK_NAN_DIAG(P2)};
  if (var.p3.on)         {CHECK_NAN_DIAG(P3)};
  if (var.theta2.on)     {CHECK_NAN_DIAG(THETA2)};
  if (var.h2.on)         {CHECK_NAN_DIAG(H2)};
  if (var.h3.on)         {CHECK_NAN_DIAG(H3)};
  if (var.t2.on)         {CHECK_NAN_DIAG(T2)};
  if (var.t3.on)         {CHECK_NAN_DIAG(T3)};
  if (var.rho2.on)       {CHECK_NAN_DIAG(RHO2)};
  if (var.rho3.on)       {CHECK_NAN_DIAG(RHO3)};
  if (var.exner3.on)     {CHECK_NAN_DIAG(EXNER3)};
  if (var.fgibb3.on)     {CHECK_NAN_DIAG(FGIBB3)};
  if (var.gz3.on)        {CHECK_NAN_DIAG(GZ3)};
  if (var.mont3.on)      {CHECK_NAN_DIAG(MONT3)};
  if (var.heat3.on)      {CHECK_NAN_DIAG(HEAT3)};
  if (var.pv3.on)        {CHECK_NAN_DIAG(PV3)};
  if (var.ri2.on)        {CHECK_NAN_DIAG(RI2)};
  if (var.div_uv3.on)    {CHECK_NAN_DIAG(DIV_UV3)};
  if (var.w3.on)         {CHECK_NAN_DIAG(W3)};
  if (var.dzdt3.on)      {CHECK_NAN_DIAG(DZDT3)};

  if (IAMNODE == NODE0) {
    fprintf(stderr,"\n");fflush(stderr);
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
      u_dat = fopen("./u_vs_lat.earth","r");
      if (!u_dat) {
        u_dat = fopen(EPIC4_PATH"/data/earth/u_vs_lat.earth","r");
      }
      if (!u_dat) {
        sprintf(Message,"Failed to open file %s",EPIC4_PATH"/data/earth/u_vs_lat.earth");
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
    u_dat = fopen("./u_vs_lat.jupiter","r");
    if (!u_dat) {
      u_dat = fopen(EPIC4_PATH"/data/jupiter/u_vs_lat.jupiter","r");
    }
    if (!u_dat) {
      sprintf(Message,"Failed to open file %s",EPIC4_PATH"/data/jupiter/u_vs_lat.jupiter");
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
    u_dat = fopen("./u_vs_lat.saturn","r");
    if (!u_dat) {
      u_dat = fopen(EPIC4_PATH"/data/saturn/u_vs_lat.saturn","r");
    }
    if (!u_dat) {
      sprintf(Message,"Failed to open file %s",EPIC4_PATH"/data/saturn/u_vs_lat.saturn");
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
    u_dat = fopen("./u_vs_lat.uranus","r");
    if (!u_dat) {
      u_dat = fopen(EPIC4_PATH"/data/uranus/u_vs_lat.uranus","r");
    }
    if (!u_dat) {
      sprintf(Message,"Failed to open file %s",EPIC4_PATH"/data/uranus/u_vs_lat.uranus");
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
    u_dat = fopen("./u_vs_lat.neptune","r");
    if (!u_dat) {
      u_dat = fopen(EPIC4_PATH"/data/neptune/u_vs_lat.neptune","r");
    }
    if (!u_dat) {
      sprintf(Message,"Failed to open file %s",EPIC4_PATH"/data/neptune/u_vs_lat.neptune");
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
 * Useful for specifying the variation of zonal wind with pressure.
 */

EPIC_FLOAT u_amp(planetspec *planet,
                 EPIC_FLOAT  p)
{
  EPIC_FLOAT
    p0;
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
      case VENUS_INDEX:
        p0 = 87.47*100.;
        /*
         * Pioneer Venus profile.
         * Normalize to 87.47 hPa (location of max u in profile).
         */
        return grid.du_vert*(pioneer_venus_u(p)/pioneer_venus_u(p0)-1.)+1.;
      break;
      default:
        sprintf(Message,"planet=%s not yet implemented",planet->name);
        epic_error(dbmsname,Message);
      break;
    }
  }
}

/*======================= end of u_amp() ====================================*/

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
    sprintf(infile,EPIC4_PATH"/data/jupiter/u_vs_p.jupiter.GalileoProbe");
    u_vs_p = fopen(infile,"r");
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
      fscanf(u_vs_p, "%*lf %*lf %lf %lf",pdat+nn,udat+nn);
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
      u = up_table[0].y;
    }
    else if (neg_log_p > up_table[nup-1].x) {
      u = up_table[nup-1].y;
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
    sprintf(infile,EPIC4_PATH"/data/venus/u_vs_p.venus.Pioneer");
    u_vs_p = fopen(infile,"r");
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

/*======================= p_sigmatheta() ====================================*/

/*
 * Returns the value of pressure consistent with 
 * theta, sigmatheta, pbot and ptop.
 */

EPIC_FLOAT p_sigmatheta(EPIC_FLOAT  theta,
                        EPIC_FLOAT  sigmatheta,
                        EPIC_FLOAT  pbot,
                        EPIC_FLOAT  ptop)
{
  int
    error_flag;
  EPIC_FLOAT
    sgtol,sg1,sg2,sg_root,
    p;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="p_sigmatheta";

  if (theta < sigmatheta) {
    sprintf(Message,"theta=%g < sigmatheta=%g",theta,sigmatheta);
    epic_error(dbmsname,Message);
  }

  SGTHMSGTH_theta      = theta;
  SGTHMSGTH_sigmatheta = sigmatheta;

  /*
   * Specify tolerance and search range.
   */
  sgtol = pow(machine_epsilon(),2./3.);

  sg1   = 0.;
  sg2   = 1.;

  error_flag = find_root(sg1,sg2,sgtol,&sg_root,sgth_minus_sgth);
  if (error_flag != 0){
    sprintf(Message,"no root found between sg1,sg2=%11.3g %11.3g",sg1,sg2);
    epic_error(dbmsname,Message);
  }

  p = get_p_sigma(pbot,sg_root,ptop);

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

/* * * * * * * * * * * end of epic_funcs_diag.c  * * * * * * * * * * * * * * */




















