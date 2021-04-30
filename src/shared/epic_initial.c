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

/* * * * * * * * * * epic_initial.c  * * * * * * * * * * * * * * * * * * * * 
 *                                                                         *
 *  Timothy E. Dowling                                                     *
 *                                                                         *
 *  Creates zonally-symmetric initial epic.nc input file for the           *
 *  EPIC atmospheric model.                                                *
 *                                                                         *
 *  The horizontal differencing scheme uses the staggered C grid.          *
 *                                                                         *
 *  For the full-globe geometry, the numbering is:                         *
 *                                                                         *
 *                            -180 deg             +180 deg                *
 *                               |                   |                     *
 *       j = nj+1   -- pv---0----pv---0----pv---0---pv -- north pole       *
 *                     |    :    |    :    |    :    |                     *
 *                     |    :    |    :    |    :    |                     *
 *       j = nj+1/2 -- u...(h)...u...(h)...u...(h)...u                     *
 *                     |    :    |    :    |    :    |                     *
 *       j = nj     -- pv---v----pv---v----pv---v---pv                     *
 *                     |    :    |    :    |    :    |                     *
 *                     u...(h)...u11.(h)11.u...(h)...u -- equator          *
 *                     |    :    |    :    |    :    |                     *
 *       j = 1      -- pv---v---pv11--v11--pv---v---pv                     *
 *                     |    :    |    :    |    :    |                     *
 *       j = 1/2    -- u...(h)...u...(h)...u...(h)...u                     *
 *                     |    :    |    :    |    :    |                     *
 *                     |    :    |    :    |    :    |                     *
 *       j = 0      -- pv---0----pv---0----pv---0---pv -- south pole       *
 *                     |         |         |         |                     *
 *                   i = 0     i = 1     i = ni    i = ni+1                *
 *                                                                         *
 *   The EPIC model's vertical coordinate is based in part on the work of  *
 *   Konor and Arakawa (1997), who generalize Hsu and Arakawa (1990)       *
 *   with a hybrid sigma-theta vertical coordinate. This combines the      *
 *   benefits of isentropic coordinates aloft with the terrain-following   *
 *   sigma coordinate near the surface.                                    *
 *                                                                         *
 *   In the EPIC model, the hybrid coordinate is called "sigmatheta" and   *
 *   sometimes abbreviated "sgth."                                         *
 *                                                                         *
 *   We carry u and v on the interfaces, and h in the layer, which is      *
 *   indicated in the figure above by writing "(h)" instead of "h".        *
 *   We also carry theta and p on the interfaces (Charney-Phillips grid).  *
 *                                                                         *
 *     k = 1/2    .......u,v,p,theta..... sigmatheta[1]      = sgth_top    *
 *                                                                         *
 *     k = 1             h                sigmatheta[2]                    *
 *                                                                         *
 *     k = 3/2    .......u,v,p,theta..... sigmatheta[3]                    *
 *                                                                         *
 *     k = nk            h                sigmatheta[2*nk]                 * 
 *                                                                         *
 *     k = nk+1/2 .......u,v,p,theta..... sigmatheta[2*nk+1] = sgth_bot    *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

#define G_SIGMA_DIAG 0.1

/*
 *  Data from the books "Venus," "Mars," "Saturn," "Uranus," 
 *  University of Arizona Press; Lindal et al (1983, Icarus) for Titan. 
 *
 *  NOTE: planet->cp and planet->kappa are reassigned according to the value of cpr 
 *        returned from the initialization routine thermo_setup(), which is the
 *        low-temperature-limit value of cp/rgas, in order to maintain self-consistency.
 *        Use the function return_cp() rather than planet->cp unless you actually want
 *        the low-temperature reference value.  
 */

init_defaultspec 
  defaults;

/*
 * Function prototypes:
 */
void read_defaults(init_defaultspec *def);

void write_defaults(init_defaultspec *def);

/*
 * The growth factor 1.+ALPHA_BOUNDARY is used
 * for geometrically decreasing boundary-layer spacing 
 * in the top sponge and bottom plaentary boundary layer (PBL).
 */
#undef  ALPHA_BOUNDARY
#define ALPHA_BOUNDARY 0.50
   
/*======================= main() ============================================*/

int main(int   argc,
         char *argv[])
{
  char   
    header[N_STR],       /*  character string buffer                     */
    system_id[32],       /*  Name of planet or system to study           */
    infile[FILE_STR],    /*  used to input data                          */
    outfile[FILE_STR],   /*  used to output data                         */
    sflag[80],           /*  string to hold command-line flags           */
    out_str[128];        /*  string to hold user inquiries               */ 
  register int    
    spacing_type,        /*  choice of layer spacing                     */
    ki,ii,index,ip,      /*  utility indices                             */
    K,J,I,kk,            /*  counters                                    */
    floor_tp,            /*  bottom index of shortened t_vs_p data       */
    ceiling_tp,          /*  top index of shortened t_vs_p data          */
    is,                  /*  species index                               */
    count,
    nk_sponge,nk_reg,nk_pbl;
  unsigned int
    time_index = 0;
  int
    itmp,
    num_advection_schemes   = 3,
    num_uv_timestep_schemes = 2,
    num_turbulence_schemes  = 2,
    num_stability_factors   = 2;
  /*
   * NOTE: Store the 1st-order "Upwind" advection scheme last on the list, and 
   *       do not offer it as an option, since it is highly dissipative.
   */
  const char
    *advection_scheme[]
      = {"Predictor-corrector (Hsu, Konor, and Arakawa)",
         "Monotonized, 4th-order centered (Akima)",
         "Upwind"},
    *uv_timestep_scheme[]
      = {"3rd-order Adams-Bashforth",
         "Leapfrog (Asselin filtered)"},
    *turbulence_scheme[]
      = {"off",
         "Spalart-Allmaras DES"},
    *stability_factor[]
      = {"off",
         "Collins et al. (2004)"},
    *Phase_Name[MAX_NUM_PHASES] 
      = {"vapor","liquid","solid","rain","snow"};
  EPIC_FLOAT
    fgibb,fpe,uoup,
    theta,theta_ortho,theta_para,
    temperature,
    ptop,pbot;
  double
    max_nu[MAX_NU_ORDER+1];
  register EPIC_FLOAT
    tmp,tmp2,            /*  temporary storage                          */
    fpara,pressure,
    sigma,sg1,sg2,slope,theta_knee,sigma_knee,
    sgth,log_sgth,neglogp,negp,
    dlnp;
  EPIC_FLOAT  
    *neglogpdat,         /*  -log(pdat)                                 */
    *thetadat,           /*  theta corresponding to t_vs_p data         */
     p_d,sgth_d;
  float_triplet
    *t_cool_table;
  FILE
    *t_cool_vs_p; 
  float_triplet
    *buff_triplet;
  struct tm
    date_start;
  time_t
    current_time;
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
    dbmsname[]="epic_initial";

#if defined(EPIC_MPI)
  MPI_Init(&argc,&argv);
  para.comm = MPI_COMM_WORLD;
  MPI_Errhandler_set(para.comm,MPI_ERRORS_RETURN);
  MPI_Comm_rank(para.comm,&para.iamnode);
  MPI_Comm_size(para.comm,&para.nproc);
  para.ndim = NINT(log((EPIC_FLOAT)para.nproc)/log(2.));
#endif
  
  /* Model version number: */
  grid.data_version = 4.30;

  declare_copyright();

  /* 
   * Interpret command-line arguments: 
   */
  if (argc > 1) {
    for (count = 1; count < argc; count++) {
      sscanf(argv[count],"%s",sflag);
      if (strcmp(sflag,"-help") == 0 ||
          strcmp(sflag,"-h")    == 0) {
        /* Print help, exit: */
        system("more "EPIC4_PATH"/help/epic_initial.help");
        exit(1);
      }
      else {
        fprintf(stderr,"Unrecognized epic_initial command-line flag: %s \n",sflag);
        exit(1);
      }
    }
  }

  /* 
   * Use 1 bar for the reference pressure.  
   * NOTE: This value is assumed in the thermodynamics subroutines.
   */
  grid.press0 = 1.e+5;  /* 1 bar in Pa */

  /*
   *  Read default parameter settings:
   */
  read_defaults(&defaults);

  if (IAMNODE == NODE0) {
    /*  Print welcome statement: */
    fprintf(stdout,"\n");    
    fprintf(stdout,"        ______/      ___       __   ___/        __     \n");
    fprintf(stdout,"       /            /     /       /          /   __/ \n");  
    fprintf(stdout,"      ______/      ______/       /          /       \n");
    fprintf(stdout,"     /            /             /          (       \n");    
    fprintf(stdout," _________/   ___/       _________/    \\_______/\n");
    fprintf(stdout," \n\n             WELCOME TO THE EPIC MODEL   \n\n");
    fprintf(stdout,"             Version: %4.2f\n",grid.data_version);
#if EPIC_PRECISION == DOUBLE_PRECISION
    fprintf(stdout,"      Floating-point: double precision\n");
#else
    fprintf(stdout,"      Floating-point: single precision\n");
#endif
    fprintf(stdout,"          Geometries: globe, f-plane \n");
    fprintf(stdout,"         Atmospheres: Venus, Earth, Jupiter, Saturn, Titan, Uranus, Neptune \n");
    fprintf(stdout,"          Benchmarks: Held_Suarez, Venus_LLR05,\n");
    fprintf(stdout,"  Under construction: Mars, Triton, Pluto, Brown Dwarfs, Extrasolar Planets \n\n");
  }

  /*
   * Choose atmosphere or benchmark to initialize:
   */
  ii = FALSE;
  while (ii == FALSE) {
    char
      *ptr;

    input_string("Choose atmosphere or benchmark to initialize \n",defaults.system_id,system_id);

    /* 
     * Convert system_id to lower case.
     */
    ptr = system_id;
    while (*ptr) {
      *ptr = (char)tolower((int)*ptr);
      ptr++;
    }

    /* 
     * Set up system.
     *
     * NOTE: When adding a new system, be sure to include the following in the list of
     *       files requiring modification: epic.h, epic_globals.c,
     *       and epic_funcs_diag.c. The easiest way is to start with a similar
     *       existing system and edit it appropriately. 
     */
    ii = TRUE;
    if      (strcmp(system_id,"venus")   == 0) {
      planet                  = &venus;
    }
    else if (strcmp(system_id,"earth") == 0) {
      planet                  = &earth;
    }
    else if (strcmp(system_id,"mars") == 0) {
      planet                  = &mars;
    }
    else if (strcmp(system_id,"jupiter") == 0) {
      planet                  = &jupiter;
    }
    else if (strcmp(system_id,"saturn")  == 0) {
      planet                  = &saturn;
    }
    else if (strcmp(system_id,"titan") == 0) {
      planet                  = &titan;
    }
    else if (strcmp(system_id,"uranus")  == 0) {
      planet                  = &uranus;
    }
    else if (strcmp(system_id,"neptune") == 0) {
      planet                  = &neptune;
    }
    else if (strcmp(system_id,"pluto") == 0) {
      planet                  = &pluto;
    }
    else if (strcmp(system_id,"hot_jupiter") == 0) {
      planet                  = &hot_jupiter;
    }
    else if (strcmp(system_id,"held_suarez") == 0) {
      planet                  = &held_suarez;
    }
    else if (strcmp(system_id,"venus_llr05") == 0) {
      planet                  = &venus_llr05;
    }
    else {
      /* unimplemented system */
      ii = FALSE;
      strcpy(defaults.system_id,"Jupiter");
      if (IAMNODE == 0) {
        fprintf(stderr,"\n\"%s\" is not defined.\n\n", system_id);
      }
    }
  }

  /*
   *  Initialize time.
   */
  strftime(header,N_STR,"%Y_%m_%d_%H:%M:%S (UTC)",gmtime(&defaults.start_time));
  if (defaults.start_date_input_type == -1) {
    defaults.start_date_input_type = 0;
    sprintf(Message,"Input starting date and time for model's time dimension:\n"
                    " 0 => default, %s",header);
  }
  else {
    sprintf(Message,"Input starting date and time for model's time dimension:\n"
                    " 0 => previous, %s",header);
  }
  fprintf(stdout,"%s\n",Message);
  current_time = time(NULL);
  strftime(header,N_STR,"%Y_%m_%d_%H:%M:%S (UTC)",gmtime(&current_time));
  sprintf(Message," 1 => current date and time, %s",header);
  fprintf(stdout,"%s\n",Message);
  defaults.start_date_input_type = 
    input_int(" 2 => prompt for values (UTC)\n",defaults.start_date_input_type);
  if (defaults.start_date_input_type == 0) {
    var.start_time = defaults.start_time;
  }
  else if (defaults.start_date_input_type == 1) {
    var.start_time = time(NULL);
  }
  else if (defaults.start_date_input_type == 2) {
    /* Trigger timezone evaluation. */
    mktime(&date_start);

    defaults.UTC_start = *gmtime(&defaults.start_time);
    date_start.tm_year = input_int("  Input year (YYYY):",defaults.UTC_start.tm_year+1900)-1900;
    date_start.tm_mon  = input_int("   Input month (MM):",defaults.UTC_start.tm_mon+1)-1;
    date_start.tm_mday = input_int("     Input day (DD):",defaults.UTC_start.tm_mday);
    date_start.tm_hour = input_int("  Input hour (0-23):",defaults.UTC_start.tm_hour);
    date_start.tm_min  = input_int("Input minute (0-59):",defaults.UTC_start.tm_min);
    date_start.tm_sec  = input_int("Input second (0-59):",defaults.UTC_start.tm_sec);
    /*
     * Shift UTC input into local time (because "struct tm" does not keep track of timezone and assumes time is local).
     */
    date_start.tm_hour -= timezone/3600;  /* timezone is defined in time.h */

    /*
     * Make and store calendar time for starting point.
     */
    var.start_time = mktime(&date_start);
  }
  else {
    sprintf(Message,"unrecognized defaults.start_date_input_type=%d",defaults.start_date_input_type);
    epic_error(dbmsname,Message);
  }
  var.model_time = defaults.start_time = var.start_time;

  /*
   * Set solar longitude, L_s [deg], which is a function of time.
   */
  L_s = solar_longitude(planet,var.model_time);

  /*
   * NOTE: Currently, we have the New Mexico Tech (NMT) physics option for Earth
   *       and the Palotai & Dowling cloud microphysics option otherwise. 
   *       These choices will evolve as these physics components are further developed.
   */
  if (strcmp(planet->name,"earth") == 0) {
    defaults.microphysics_on = grid.microphysics_on = FALSE;
    defaults.nmt_physics_on  = grid.nmt_physics_on =
                            input_int("NMT physics on or off? [1 or 0]\n",defaults.nmt_physics_on);
  }
  else {
    defaults.nmt_physics_on  = grid.nmt_physics_on = FALSE;
    defaults.microphysics_on = grid.microphysics_on =
                           input_int("Cloud microphysics on (latent heat, phase changes) or off (passive)? [1 or 0]\n",defaults.microphysics_on);
  }

 /*
  * Prompt for initial winds.
  *
  * u_scale - initial wind scaling factor
  * du_vert - a flag: 1 = use vert. wind profile; 0 = du is constant w/ height
  */
  switch(planet->index) {
    case EARTH_INDEX:
      defaults.u_scale = 0.;
      grid.du_vert     = defaults.du_vert = 0.;
    break;
    case HELD_SUAREZ_INDEX:
      defaults.u_scale = 0.;
      grid.du_vert     = defaults.du_vert = 0.;
    break;
    case VENUS_LLR05_INDEX:
      defaults.u_scale = 0.;
      grid.du_vert     = defaults.du_vert = 0.;
    break;
    case TITAN_INDEX:
      defaults.u_scale = 0.;
      grid.du_vert     = defaults.du_vert = 0.;
    break;
    default:
      defaults.u_scale = input_float("Input initial-wind scaling factor [zero=0., full strength=1.]\n",defaults.u_scale);
      if (defaults.u_scale != 0.) {
        /*
         * Prompt for vertical variation of zonal wind.
         */
        if (strcmp(planet->name,"jupiter") == 0 ||
            strcmp(planet->name,"venus")   == 0) {
          fprintf(stdout,"Vertical wind profile: 1.0 => use full probe profile \n");
          fprintf(stdout,"                       0.0 => constant with height \n");
          grid.du_vert = defaults.du_vert = input_float("",defaults.du_vert);
        }
        else {
          grid.du_vert = defaults.du_vert = 0.;
        }
      }
      else {
        grid.du_vert = defaults.du_vert = 0.;
      }
    break;
  }
    
  /*
   * Inquire about Newtonian cooling and internal Rayleigh drag:
   *
   *     grid.newt_cool_on - flag for newtonian cooling
   * grid.newt_cool_adjust - 1: sets layer avg. of newt. cooling to 0
   *          grid.prandtl - relaxation ratio (tau_rad/tau_drag); 0 = no drag
   *         grid.tau_drag - ray. drag timescl (s) (nu0 = 1/tau_drag) (1e20 => inf.)
   *           grid.drag_v - 1: rayleigh drag is applied to v
   *   grid.drag_zonal_avg - 1: rayleigh drag operates on zonal avg.
   */
  if (strcmp(system_id,"held_suarez") == 0) {
    defaults.newt_cool_on = grid.newt_cool_on = TRUE;
  }
  else if (strcmp(system_id,"venus_llr05") == 0) {
    defaults.newt_cool_on = grid.newt_cool_on = TRUE;
  }
  else if (strcmp(system_id,"earth") == 0) {
    defaults.newt_cool_on = grid.newt_cool_on = FALSE;
  }
  else {
    defaults.newt_cool_on = grid.newt_cool_on 
      = input_int("Newtonian cooling on or off? [1 or 0]\n",defaults.newt_cool_on);
  }

  if (grid.newt_cool_on == TRUE) {
    if (strcmp(system_id,"held_suarez") == 0) {
      defaults.newt_cool_adjust = grid.newt_cool_adjust = FALSE;
      defaults.prandtl          = grid.prandtl          = 0.;
    }
    else if (strcmp(system_id,"venus_llr05") == 0) {
      defaults.newt_cool_adjust = grid.newt_cool_adjust = FALSE;
      defaults.prandtl          = grid.prandtl          = 0.;
    }
    else {
      defaults.newt_cool_adjust = grid.newt_cool_adjust =
        input_int("Adjust layer average of Newtonian cooling to zero? [1=yes, 0=no]\n",
                   defaults.newt_cool_adjust);
      defaults.prandtl = grid.prandtl = 
        input_float("Input relaxation ratio [tau_rad/tau_drag, 0. => no drag]\n",defaults.prandtl);
    }

    if (defaults.prandtl == 0.) {
      /* For consistency. */
      defaults.tau_drag   = 1.e+20;
      defaults.drag_v     = grid.drag_v = FALSE;
      grid.drag_zonal_avg = defaults.drag_zonal_avg;
    }
    else {
      defaults.drag_v = grid.drag_v = input_int("Apply Rayleigh drag to v? [1=yes, 0=no]\n",defaults.drag_v);
      defaults.drag_zonal_avg = grid.drag_zonal_avg = 
        input_int("Apply Rayleigh drag to zonal average velocity components? [1=yes, 0=no]\n",defaults.drag_zonal_avg);
    }

    /*
     * Set up t_cool_table.
     * Input t_cool_vs_p, looking first in the local directory.
     */
    sprintf(infile,"./t_cool_vs_p.%s",planet->name);
    t_cool_vs_p = fopen(infile,"r");
    if (!t_cool_vs_p) {
      sprintf(infile,EPIC4_PATH"/data/%s/t_cool_vs_p.%s",planet->name,planet->name);
      t_cool_vs_p = fopen(infile,"r");
    }
    if (t_cool_vs_p) {
      EPIC_FLOAT
        t_spinup;
      FILE
        *relax_times;

      for (ki = 0; ki < 6; ki++) {
        fgets(header,N_STR,t_cool_vs_p);
      }
      /* input number of data points */
      fscanf(t_cool_vs_p, "%d", &var.n_t_cool); 

      /* 
       * Allocate memory for t_cool_table.
       * We'll calculate the values and store them in t_cool_table,
       * and then copy them into var.t_cool_table after calling make_arrays() below.
       */
      t_cool_table = ftriplet(0,var.n_t_cool-1,dbmsname);

      relax_times = fopen("./relax_times.dat","w");
      fprintf(relax_times," %s\n",planet->name);
      fprintf(relax_times,"   p[hPa]        trad[s]       tdrag[s] \n");

      /* stored in order of decreasing p (increasing -log p) */
      for (ki = var.n_t_cool-1; ki >= 0;  ki--) {  

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(t_cool_vs_p,"%lf %*lf %lf",&t_cool_table[ki].x,&t_cool_table[ki].y);
#else
        fscanf(t_cool_vs_p,"%f %*f %f",&t_cool_table[ki].x,&t_cool_table[ki].y);
#endif

        /* convert from hPa to Pa */
        t_cool_table[ki].x *= 100.;

        /* spline on -log p  */
        t_cool_table[ki].x = -log(t_cool_table[ki].x);  

        /* Output table of relaxation times vs pressure: */
        if (grid.prandtl == 0.) {
          t_spinup = 1.e+20;
        }
        else {
          t_spinup = t_cool_table[ki].y/(grid.prandtl);
        }
        fprintf(relax_times," %e  %e  %e \n",
                exp(-t_cool_table[ki].x)/100.,t_cool_table[ki].y,t_spinup);
      }
      fclose(t_cool_vs_p);
      fclose(relax_times);
    }
    else {
      /*
       * Unable to read t_cool_vs_p data file.
       */
      var.n_t_cool = 0;
      if (strcmp(planet->name,"held_suarez") == 0) {
        /*
         * The Held-Suarez reference case does not use t_cool_table.
         */
        ;
      }
      else {
        sprintf(Message,"cannot find %s",infile);
        epic_error(dbmsname,Message);
      }
    }
  }
  else {
    if (strcmp(system_id,"earth") == 0) {
      /*
       * Turn Rayleigh drag off for Earth case.
       */
      defaults.tau_drag = grid.tau_drag = 1.e+20;
      defaults.prandtl  = grid.prandtl  = 0.;
    }
    else {
      defaults.tau_drag = grid.tau_drag = 
        input_float("Rayleigh drag timescale [s, 1.e+20 => infinity]\n",defaults.tau_drag);
    }
    if (grid.tau_drag < 1.e+20) {
      defaults.drag_v = grid.drag_v = input_int("Apply Rayleigh drag to v? [1=yes, 0=no]\n",defaults.drag_v);
      defaults.drag_zonal_avg = grid.drag_zonal_avg = 
        input_int("Apply Rayleigh drag to zonal average velocity components? [1=yes, 0=no]\n",defaults.drag_zonal_avg);
    }
    else {
      defaults.drag_v     = grid.drag_v = FALSE;
      grid.drag_zonal_avg = defaults.drag_zonal_avg;
    }
  }

  /*
   * Inquire about strength of spinup wind.
   *
   * u_spinup_scale -  scaling factor for u_spinup. If zero, u_spinup = 0
   *                   if 1, u_spinup = u, else, u_spinup = u*u_spinup_scale.
   */
  switch(planet->index) {
    case HELD_SUAREZ_INDEX:
      /* defaults.u_spinup_scale not used in this case. */
      defaults.u_spinup_scale = 0.;
    break;
    case VENUS_LLR05_INDEX:
      defaults.u_spinup_scale = 0.;
    break;
    case TITAN_INDEX:
      defaults.u_spinup_scale = 0.;
    break;
    default:
      if (defaults.prandtl != 0. || defaults.tau_drag < 1.e+20) {
        defaults.u_spinup_scale = input_float("Input spinup zonal-wind scaling factor [1.=initial u, 0.=zero]\n",defaults.u_spinup_scale);
      }
      else {
        /*
         * Set u_spinup to be the same as the initial U.
         */
        defaults.u_spinup_scale = 1.;
      }
    break;
  }

  /* 
   * Prompt for timestep scheme for momentum variables (u,v).
   *
   * grid.uv_timestep_scheme - name of timestep scheme
   */
  fprintf(stdout,"\n");
  sprintf(header,"Input timestep scheme for momentum variables (u,v): 0 => %s \n",
                 uv_timestep_scheme[0]);
  for (ii = 1; ii < num_uv_timestep_schemes; ii++) {
    sprintf(Message,"                                                   %2d => %s \n",
            ii,uv_timestep_scheme[ii]);
    strcat(header,Message);
  }
  defaults.uv_timestep_scheme = input_int(header,defaults.uv_timestep_scheme);
  strcpy(grid.uv_timestep_scheme,uv_timestep_scheme[defaults.uv_timestep_scheme]);

  /* 
   * Prompt for turbulence scheme.
   *
   * grid.turbulence_scheme - name of turb. scheme
   */
  fprintf(stdout,"\n");
  sprintf(header,"Input turbulence scheme: 0 => %s \n",
                 turbulence_scheme[0]);
  for (ii = 1; ii < num_turbulence_schemes; ii++) {
    sprintf(Message,"                        %2d => %s \n",
            ii,turbulence_scheme[ii]);
    strcat(header,Message);
  }
  defaults.turbulence_scheme = input_int(header,defaults.turbulence_scheme);
  strcpy(grid.turbulence_scheme,turbulence_scheme[defaults.turbulence_scheme]);


  /*
   * Inquire whether diffusion should be applied to both horizontal and vertical directions, 
   * just horizontal, or just vertical.
   *
   * grid.diffusion_direction - 1=horiz.+vert.;2=horiz.; 3=vert. 
   *
   * Inquire about Richardson-number-based (Ri-based) stability factor.
   */
  if (strcmp(grid.turbulence_scheme,"off") == 0) {
    defaults.diffusion_direction = grid.diffusion_direction = 0;
    strcpy(grid.stability_factor,"off");
  }
  else {
    if (defaults.diffusion_direction == 0) {
      defaults.diffusion_direction = 1;
    }
    defaults.diffusion_direction = grid.diffusion_direction =
                             input_int("Diffusion applied to both horizontal and vertical directions, just horizontal, or just vertical? [1,2,3]\n",
                                        defaults.diffusion_direction);
    fprintf(stdout,"\n");
    sprintf(header,"Input Ri-based stability factor scheme: 0 => %s \n",
                    stability_factor[0]);
    for (ii = 1; ii < num_stability_factors; ii++) {
      sprintf(Message,"                                       %2d => %s \n",
                      ii,stability_factor[ii]);
      strcat(header,Message);
    }
    defaults.stability_factor = input_int(header,defaults.stability_factor);
    strcpy(grid.stability_factor,stability_factor[defaults.stability_factor]);
  }

  /*
   * Grab the vertical u profile, if it exists.
   * Malloc the string space, and name the indices, and units for each 
   * prognostic and diagnostic into var.(var_name).info[0].index,name,units
   * Also set enthalpy change and saturation vapor pressure for species.
   * Also, set molar mass and triple-point values for species.
   */
  set_var_props();

  /*
   * Prompt for advection scheme for h-grid mass variables (HDRY, Qs).
   *
   * NOTE: This must come before make_arrays().
   *
   * NOTE: The 1st-order "Upwind" advection scheme at the end of the list is not
   *       offered as an option, since it is highly dissipative.
   */
  fprintf(stdout,"\n");
  sprintf(header,"Input advection scheme for mass: 0 => %s \n",advection_scheme[0]);
  for (ii = 1; ii < num_advection_schemes-1; ii++) {
    sprintf(Message,"                                %2d => %s \n",
            ii,advection_scheme[ii]);
    strcat(header,Message);
  }
  defaults.hgrid_mass_advection_scheme = input_int(header,defaults.hgrid_mass_advection_scheme);
  strcpy(var.hdry.advection_scheme,advection_scheme[defaults.hgrid_mass_advection_scheme]);
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    strcpy(var.species[is].advection_scheme,advection_scheme[defaults.hgrid_mass_advection_scheme]);
  }

  /*
   * Prompt for advection scheme for h-grid non-mass variables (THETA, FPARA, NU_TURB).
   */
  fprintf(stdout,"\n");
  sprintf(header,"Input advection scheme for theta: 0 => %s \n",advection_scheme[0]);
  for (ii = 1; ii < num_advection_schemes-1; ii++) {
    sprintf(Message,"                                 %2d => %s \n",
            ii,advection_scheme[ii]);
    strcat(header,Message);
  }
  defaults.hgrid_nonmass_advection_scheme = input_int(header,defaults.hgrid_nonmass_advection_scheme);
  strcpy(var.theta.advection_scheme,  advection_scheme[defaults.hgrid_nonmass_advection_scheme]);
  strcpy(var.fpara.advection_scheme,  advection_scheme[defaults.hgrid_nonmass_advection_scheme]);
  strcpy(var.nu_turb.advection_scheme,advection_scheme[defaults.hgrid_nonmass_advection_scheme]);

  /*
   * Set equation of state to ideal.
   * NOTE: We have code available in the model for handling the virial equation of state,
   *       but we have not run across any case in the solar system that benefits signficantly from it,
   *       so we hardwire the model to use the ideal eos.
   */
  sprintf(grid.eos,"ideal");
  sprintf(defaults.eos,"ideal");

  /*
   * Inquire about vertical spacing between layers:
   */
  if (IAMNODE == 0) {
    fprintf(stdout,"\n");
    sprintf(header,"Layer intervals even in: %1d => log pressure, with increased resolution in boundary layers\n"
                   "                         %1d => log pressure\n"
                   "                         %1d => pressure \n"
                   "                         %1d => from file \n",
                   SPACING_LOGP_W_BOUNDARIES,SPACING_LOGP,SPACING_P,SPACING_FROM_FILE);
  }
  defaults.spacing_type = spacing_type = input_int(header,defaults.spacing_type);

  /*
   * Inquire about the vertical range and number of vertical layers.
   *
   * Grab the file, if user wants to read one.
   * read_spacing_file() is a simple template for reading ascii data files.
   * defaults.ptop - p at k=1/2 (Pa) (used below in determining zeta)
   * grid.pbot = p at k=nk+1/2 (Pa)
   * grid.nk = user's number of vert. layers
   */
  if (defaults.spacing_type == SPACING_FROM_FILE) {
    input_string("File containing layer spacing data (type 'initial -h' to see an example)\n",
                 defaults.layer_spacing_dat,defaults.layer_spacing_dat);
    read_spacing_file(&defaults,SIZE_DATA);
  }
  else {
    /* 
     * Determine ptop (external units are hPa, internal are Pa):
     */
    defaults.ptop = 100.*input_float(
                           "Target pressure at k = 1/2 (model's top) [hPa]\n",
                           defaults.ptop/100.);

    if (strcmp(planet->type,"gas-giant") == 0) {
      /*
       * Determine pbot for gas giant.
       */
      defaults.pbot = grid.pbot = 100.*
                                  input_float("Target pressure at k = nk+1/2 (model's bottom) [hPa]\n",
                                              defaults.pbot/100.);
    }

    if (defaults.nk == -1) {
      /* 
       * Smart toggle case.
       * Prompt nk ~ 2*(model range in scale heights).
       */
      if (strcmp(planet->type,"gas-giant") == 0) {
        defaults.nk = 2*NINT(log(defaults.pbot/defaults.ptop));
      }
      else {
        defaults.nk = 20;
      }
    }
    defaults.nk = grid.nk = input_int("\nInput the number of vertical layers, nk\n",defaults.nk);
  }

  /*
   * Compute hyperviscosity coefficients.
   *
   * max_nu - an array of maximum numerically stable hyperviscosity values
   */
  set_max_nu(max_nu);

  /*
   * Choose geometry.
   *
   *           grid.geometry - name of geometry
   *             grid.wrap[] - periodicity flags for each dim. (globe: 0,0,1)
   *              grid.pad[] - boundary pad widths for each dim. (globe: 1,1,1)
   *                grid.jlo - lower index for j (globe: jlo = 0)  
   *       grid.globe_latbot - lowest latitude value
   *       grid.globe_lattop - highest latitude value
   *       grid.globe_lonbot - lowest longitude value
   *       grid.globe_lontop - highest longitude value
   *       grid.f_plane_lat0 - central latitude of f_plane
   * grid.f_plane_half_width - half of the f_plane width
   */
  if (strcmp(system_id,"held_suarez") == 0) {
    sprintf(grid.geometry,    "globe");
    sprintf(defaults.geometry,"globe");
  }
  else if (strcmp(system_id,"venus_llr05") == 0) {
    sprintf(grid.geometry,    "globe");
    sprintf(defaults.geometry,"globe");
  }
  else {
    input_string("Choose geometry \n",defaults.geometry,grid.geometry);
  }

  if (strcmp(grid.geometry,"f-plane") == 0) {
    grid.wrap[2] = 0;
    grid.wrap[0] = 1;
    grid.pad[2]  = 1;
    grid.pad[0]  = 1;
    grid.f_plane_lat0 = defaults.f_plane_lat0
                      = input_float("Latitude of f-plane [deg] \n",
                                     defaults.f_plane_lat0);
    input_string("Choose mapping: cartesian or polar \n",
                 defaults.f_plane_map,grid.f_plane_map);
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      grid.wrap[1] = 1;
      grid.pad[1]  = 1;
      grid.jlo     = 1;
      grid.f_plane_half_width = defaults.f_plane_half_width
                              = input_float("Half-width [km] \n",
                                             defaults.f_plane_half_width);
      /* convert km to m */
      grid.f_plane_half_width *= 1000.;
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      grid.wrap[1] = 0;
      grid.pad[1]  = 1;
      grid.jlo     = 0;
      grid.f_plane_half_width = defaults.f_plane_half_width
                              = input_float("Radius [km] \n",
                                            defaults.f_plane_half_width);
      /* convert km to m */
      grid.f_plane_half_width *= 1000.;
    }
    else {
      fprintf(stderr,"Unrecognized f-plane mapping. \n");
      exit(1);
    }
  }
  else if (strcmp(grid.geometry,"globe") == 0) {
    grid.wrap[2] = 0;
    grid.wrap[1] = 0;
    grid.wrap[0] = 1;
    grid.pad[2]  = 1;
    grid.pad[1]  = 1;
    grid.pad[0]  = 1;
    grid.jlo     = 0;

    grid.globe_latbot=defaults.globe_latbot
      = input_float("Lowest latitude [deg] \n",defaults.globe_latbot);
    if (grid.globe_latbot < -90.) {
      fprintf(stderr,"latbot must be >= -90.  Setting latbot to -90. \n");
      grid.globe_latbot = defaults.globe_latbot = -90.;
    }
    grid.globe_lattop=defaults.globe_lattop
      = input_float("Highest latitude [deg] \n",defaults.globe_lattop);
    /* Sanity check: */
    if (grid.globe_lattop < grid.globe_latbot) {
      sprintf(Message,"lattop=%f < latbot=%f",grid.globe_lattop,grid.globe_latbot);
      epic_error(dbmsname,Message);
    }
    if (grid.globe_lattop > 90.) {
      fprintf(stderr,"lattop must be <= 90. Setting lattop to 90. \n");
      grid.globe_lattop = defaults.globe_lattop = 90.;
    }

    grid.globe_lonbot=defaults.globe_lonbot
      = input_float("Lowest longitude [deg] \n",defaults.globe_lonbot);
    grid.globe_lontop=defaults.globe_lontop
      = input_float("Highest longitude [deg] \n",defaults.globe_lontop);
    /* Sanity check: */
    if (grid.globe_lontop < grid.globe_lonbot) {
      sprintf(Message,"lontop=%f < lonbot=%f",grid.globe_lontop,grid.globe_lonbot);
      epic_error(dbmsname,Message);
    }
  }
  else {
    fprintf(stderr,"initial: unrecognized geometry: %s \n",grid.geometry);
    exit(1);
  }

  /*
   * Set up thermodynamics subroutines.
   *
   * This fills in the 'thermo' struct with planet-relevant values.
   * It returns planet->cpr, which is the planet's reference cp.
   * The (->) operator allows us to reference the member of a pointed-to struct.
   * Here, we are passing the address of planet->cpr, so that its value can be set
   * by the subroutine. 
   * 
   *   planet->cpr - nondimensional ref. cp
   *    planet->cp - spec. heat at const. p
   *  planet->rgas - gas constant
   * planet->kappa - rgas/cp
   */
  thermo_setup(planet,&planet->cpr);
  /* Assign thermodynamics function's reference cpr to planet->cp */
  planet->cp    = planet->cpr*planet->rgas;
  planet->kappa = 1./planet->cpr;

  /* 
   * Turn on core prognostic variables, u,v,hdry,theta.
   *
   * var.(prog_name).on - a flag to tell us it is on
   * var.on_list[INDEX] - an array of flags that does the same thing
   */
  var.u.on     = var.on_list[U_INDEX    ] = TRUE;
  var.v.on     = var.on_list[V_INDEX    ] = TRUE;
  var.hdry.on  = var.on_list[HDRY_INDEX ] = TRUE;
  var.theta.on = var.on_list[THETA_INDEX] = TRUE;

  /*
   * Turn on turbulence-model variables.
   */
  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
    var.nu_turb.on = var.on_list[NU_TURB_INDEX] = TRUE;
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  /*
   * Turn on diagnostic variables that are used
   * in the model.  Variables not turned on here can
   * still be written to extract.nc.
   */
  var.hdry3.on      = var.on_list[HDRY3_INDEX     ] = TRUE;
  var.p2.on         = var.on_list[P2_INDEX        ] = TRUE;
  var.p3.on         = var.on_list[P3_INDEX        ] = TRUE;
  var.theta2.on     = var.on_list[THETA2_INDEX    ] = TRUE;
  var.h2.on         = var.on_list[H2_INDEX        ] = TRUE;
  var.h3.on         = var.on_list[H3_INDEX        ] = TRUE;
  var.t2.on         = var.on_list[T2_INDEX        ] = TRUE;
  var.t3.on         = var.on_list[T3_INDEX        ] = TRUE;
  var.rho2.on       = var.on_list[RHO2_INDEX      ] = TRUE;
  var.rho3.on       = var.on_list[RHO3_INDEX      ] = TRUE;
  var.exner3.on     = var.on_list[EXNER3_INDEX    ] = TRUE;
  var.gz3.on        = var.on_list[GZ3_INDEX       ] = TRUE;
  var.mont3.on      = var.on_list[MONT3_INDEX     ] = TRUE;
  var.heat3.on      = var.on_list[HEAT3_INDEX     ] = TRUE;
  var.pv3.on        = var.on_list[PV3_INDEX       ] = TRUE;
  var.w3.on         = var.on_list[W3_INDEX        ] = TRUE;
  var.dzdt3.on      = var.on_list[DZDT3_INDEX     ] = TRUE;
  var.u_spinup.on   = var.on_list[U_SPINUP_INDEX  ] = TRUE;
  if (strcmp(planet->type,"terrestrial") == 0) {
    var.gz_surface.on = var.on_list[GZ_SURFACE_INDEX] = TRUE;
  }

  if (grid.nmt_physics_on) {
    var.dry_entropy.on       = var.on_list[DRY_ENTROPY_INDEX      ] = TRUE;
    var.moist_entropy.on     = var.on_list[MOIST_ENTROPY_INDEX    ] = TRUE;
    var.sat_moist_entropy.on = var.on_list[SAT_MOIST_ENTROPY_INDEX] = TRUE;
    var.the_flux.on          = var.on_list[THE_FLUX_INDEX         ] = TRUE;
    var.rt_flux.on           = var.on_list[RT_FLUX_INDEX          ] = TRUE;
    var.u_flux.on            = var.on_list[U_FLUX_INDEX           ] = TRUE;
    var.v_flux.on            = var.on_list[V_FLUX_INDEX           ] = TRUE;
    var.convthrot.on         = var.on_list[CONVTHROT_INDEX        ] = TRUE;
    var.fluxthrot.on         = var.on_list[FLUXTHROT_INDEX        ] = TRUE;
    var.rain_rate.on         = var.on_list[RAIN_RATE_INDEX        ] = TRUE;
  }

  fprintf(stdout,"\n Core prognostic variables: \n");

  fprintf(stdout,"  %2d  ><  %-s \n",U_INDEX,    var.u.info[    0].name);
  fprintf(stdout,"  %2d  ><  %-s \n",V_INDEX,    var.v.info[    0].name);
  fprintf(stdout,"  %2d  ><  %-s \n",HDRY_INDEX, var.hdry.info[ 0].name);
  fprintf(stdout,"  %2d  ><  %-s \n",THETA_INDEX,var.theta.info[0].name); 

  /*
   * Inquire about optional prognostic variables.
   *
   * prompt_species_on()  - a function to inquire about the species, and turn 
   *                        them on, as above.
   */
  if (strcmp(system_id,"held_suarez") == 0) {
    ;
  }
  else {
    count = prompt_species_on(planet,defaults.species_str);
  }

  /*
   * Turn on phases appropriate to choice of physics package.
   * The phases are switched on here, whether or not the species are invoked. 
   */
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

    var.species[CH_4_INDEX].phase[VAPOR ].on = TRUE;
    var.species[CH_4_INDEX].phase[ICE   ].on = TRUE;
    var.species[CH_4_INDEX].phase[LIQUID].on = TRUE;
    var.species[CH_4_INDEX].phase[RAIN  ].on = TRUE;
    var.species[CH_4_INDEX].phase[SNOW  ].on = TRUE;

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

    var.species[H_2S_INDEX].phase[VAPOR ].on = TRUE;
    var.species[H_2S_INDEX].phase[ICE   ].on = TRUE;
    var.species[H_2S_INDEX].phase[LIQUID].on = TRUE;
    var.species[H_2S_INDEX].phase[RAIN  ].on = TRUE;
    var.species[H_2S_INDEX].phase[SNOW  ].on = TRUE;

    var.species[NH_4SH_INDEX].phase[VAPOR ].on = TRUE;
    var.species[NH_4SH_INDEX].phase[ICE   ].on = TRUE;
    var.species[NH_4SH_INDEX].phase[LIQUID].on = TRUE;
    var.species[NH_4SH_INDEX].phase[RAIN  ].on = TRUE;
    var.species[NH_4SH_INDEX].phase[SNOW  ].on = TRUE;
  }
  else if (grid.nmt_physics_on == TRUE) {
    /*
     * The New Mexico Tech physics scheme uses total advected water (vapor + cloud ice + cloud liquid)
     * as a single, quasi-conserved prognostic variable. We use the VAPOR phase to hold this variable,
     * and turn the others off.
     */
    var.species[H_2O_INDEX].phase[VAPOR ].on = TRUE;
    var.species[H_2O_INDEX].phase[ICE   ].on = FALSE;
    var.species[H_2O_INDEX].phase[LIQUID].on = FALSE;
    var.species[H_2O_INDEX].phase[RAIN  ].on = FALSE;
    var.species[H_2O_INDEX].phase[SNOW  ].on = FALSE;
    /*
     * Placeholders.
     */
    var.species[CO_2_INDEX].phase[VAPOR ].on = FALSE;
    var.species[CO_2_INDEX].phase[ICE   ].on = FALSE;
    var.species[CO_2_INDEX].phase[LIQUID].on = FALSE;
    var.species[CO_2_INDEX].phase[RAIN  ].on = FALSE;
    var.species[CO_2_INDEX].phase[SNOW  ].on = FALSE;

    var.species[O_3_INDEX ].phase[VAPOR ].on = FALSE;
    var.species[O_3_INDEX ].phase[ICE   ].on = FALSE;
    var.species[O_3_INDEX ].phase[LIQUID].on = FALSE;
    var.species[O_3_INDEX ].phase[RAIN  ].on = FALSE;
    var.species[O_3_INDEX ].phase[SNOW  ].on = FALSE;
  }

  /* 
   * Source-sink parameters.
   *
   * var.fpara.on - flag to indicate fpara prog. is on (turned on by user in prompt_species_on, above)
   * var.time_fp_bar - fpara conversion time at 1 bar
   */
  if (var.fpara.on) {
    defaults.time_fp = 
      input_float("Ortho-para conversion time at 1 bar [s, nominal = 3.e+8]?\n",defaults.time_fp);
    var.time_fp_bar = defaults.time_fp*1000.*100.;
  }

  /*
   * Inquire about number of latitude grid points.
   *
   * NOTE: nj equals the number of interior v or q points, counting from 1 to nj,
   *       whereas nj+1 equals the number of u or h points, counting from 0 to nj
   *       for channel models and globes. 
   *       One consequence is that for a 1D model, nj = 0.
   *
   * grid.jlo    - low index for j (globe: jlo = 0)
   * grid.nj     - number of lat. grid points
   * grid.dlt    - lat. grid spacing (deg)
   */
  if (strcmp(grid.geometry,"globe") == 0) {
    sprintf(Message,"\nLatitude indexing: nj   = no. v or q pts. interior to the boundaries, from 1 to nj, whereas"
                    "\n                   nj+1 = no. u or h points (they have no boundary pts), from 0 to nj."
                    "\nInput the number of latitude gridpoints (nj, 0 => one u,h point)\n");
    defaults.nj = grid.nj = input_int(Message,defaults.nj);
    if (grid.globe_latbot == -90. &&
        grid.globe_lattop ==  90.) {
      /* 
       * Poles are offset by extra sqrt(.5)*dlt.
       * NOTE: This is more accurate than the .5*dlt prescribed by Arakawa and Lamb.
       */
      grid.dlt = 180./((grid.nj+1)+sqrt(.5)+sqrt(.5));
    }
    else if (grid.globe_latbot == -90.) {
      grid.dlt = (grid.globe_lattop+90.)/((grid.nj+1)+sqrt(.5));
    }  
    else if (grid.globe_lattop ==  90.) {
      grid.dlt = (90.-grid.globe_latbot)/((grid.nj+1)+sqrt(.5));
    }
    else {
      grid.dlt = (grid.globe_lattop-grid.globe_latbot)/(grid.nj+1);
    }
  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      defaults.nj = grid.nj =
        input_int("\nInput the number of gridpoints on a side\n",defaults.nj);
      grid.dlt = 360./(grid.nj);
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      /* 
       * "latitude" runs from 0 at the edge (j = 1) 
       *  to 90 in the center (j = nj+1) */
      defaults.nj = grid.nj =
        input_int("\nInput the number of radial gridpoints\n",defaults.nj);
      grid.dlt   = 90./((grid.nj+1-1)+sqrt(.5));
    }
  }
  else {
    sprintf(Message,"unrecognized grid.geometry=%s",grid.geometry);
    epic_error(dbmsname,Message);
  }
  /*
   * To handle staggered C-grid.
   */
  grid.jfirst = IMAX(grid.jlo,1);
  grid.jlast  = grid.nj-1;

  /*
   * Inquire about number of longitude gridpoints.
   *
   * frexp()  - a standard c math.h function.
   * grid.ni  - number of long. grid points (must be power of 2)
   * grid.dln - long. grid spacing (deg)
   */
  if (strcmp(grid.geometry,"globe") == 0) {
    if (IAMNODE == 0) {
      fprintf(stdout,"\nInput the number of longitude gridpoints (ni), \n");
    }
    defaults.ni = grid.ni = 
      input_int("which must be an integer power of two\n",defaults.ni);
    /*
     * Verify that ni is a power of 2.
     */
    if (frexp((double)grid.ni,&itmp) != 0.5) {
      sprintf(Message,"ni=%d is not a power of 2",grid.ni);
      epic_error(dbmsname,Message);
    }

    grid.dln = (grid.globe_lontop-grid.globe_lonbot)/(grid.ni);
  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      defaults.ni = grid.ni = grid.nj;
      grid.dln    = grid.dlt;
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      if (IAMNODE == 0) {
        fprintf(stdout,"\nInput the number of longitude gridpoints, \n");
      }
      defaults.ni = grid.ni = 
        input_int("which must be an integer power of two\n",defaults.ni);
      /*
       * Verify that ni is a power of 2.
       */
      if (frexp((double)grid.ni,&itmp) != 0.5) {
        sprintf(Message,"ni=%d is not a power of 2",grid.ni);
        epic_error(dbmsname,Message);
      }
      grid.dln = 360./(grid.ni);
    }
  }

  /*
   * Input dt.
   */
  defaults.dt = grid.dt = 
    input_int("Input timestep, dt [s, int to prevent roundoff]\n",defaults.dt);

  if(grid.nmt_physics_on == 1) {
    /*
     * Set dt for nmt_physics (ks).
     */
    nmt.dt = grid.dt/1000.;
  }

  /*
   * Setup sponge in top of model.
   * grid.k_sponge - number of sponge layers
   *
   * NOTE: this needs to be done after grid.dt is set.
   */
  if (strcmp(planet->name,"held_suarez") == 0) {
    defaults.k_sponge = grid.k_sponge = 0;
  }
  else if (defaults.k_sponge == -1) {
    /* Smart toggle case. */
    if (grid.nk >= 5) {
      defaults.k_sponge = NINT((EPIC_FLOAT)grid.nk*.2);
      defaults.k_sponge = IMAX(IMIN(defaults.k_sponge,4),1);
    }
    else {
      /* Default to no sponge for nk <= 4 */
      defaults.k_sponge = 0;
    }
    defaults.k_sponge = grid.k_sponge = input_int("Input k_sponge (0 = no effect)\n",defaults.k_sponge);
  }
  else {
    defaults.k_sponge = grid.k_sponge = input_int("Input k_sponge (0 = no effect)\n",defaults.k_sponge);
  }

  var.ntp = read_t_vs_p(planet,SIZE_DATA);

  /* * * * * * * * * * * * * * * * * * *
   *                                   *
   *  Allocate memory for variables.   *
   *                                   *
   * * * * * * * * * * * * * * * * * * */
  make_arrays(planet);

  /*
   * NOTE: these must come after make_arrays().
   */
  read_t_vs_p(planet,POST_SIZE_DATA);
  set_sponge();

  if (var.n_t_cool > 0) {
    /*
     * Load var.t_cool_table.
     */
    for (ki = 0; ki < var.n_t_cool; ki++) {
      var.t_cool_table[ki].x = t_cool_table[ki].x;
      var.t_cool_table[ki].y = t_cool_table[ki].y;
    }
    free_ftriplet(t_cool_table,0,var.n_t_cool-1,dbmsname);
  }

  if (spacing_type == SPACING_FROM_FILE) {
    read_spacing_file(&defaults,ALL_DATA);
  }

  if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
    /*
     * Flag DUDT and DVDT for timeplanes IT_MINUS2 and IT_MINUS1 with FLOAT_MAX
     * to indicate this is the initial timestep.
     */
    for (K = 1; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          DUDT(IT_MINUS2,K,J,I) = FLOAT_MAX;
          DUDT(IT_MINUS1,K,J,I) = FLOAT_MAX;
          DVDT(IT_MINUS2,K,J,I) = FLOAT_MAX;
          DVDT(IT_MINUS1,K,J,I) = FLOAT_MAX;
        }
      }
    }
  }
  else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
    /*
     * Flag U and V for timeplane IT_MINUS1 with FLOAT_MAX
     * to indicate this is the initial timestep.
     */
    for (K = 1; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          U(IT_MINUS1,K,J,I) = FLOAT_MAX;
          V(IT_MINUS1,K,J,I) = FLOAT_MAX;
        }
      }
    }
  }
  else {
    sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
    epic_error(dbmsname,Message);
  }

  /*
   * NOTE: timeplane_bookkeeping() should come after setting the initial-step
   *       FLOAT_MAX flags in U or DUDT.
   */
  timeplane_bookkeeping();

  /*
   * Set lon, lat, etc.
   *
   * set_lonlat():
   * this propogates grid.lon, grid.lat (the dim arrays) with proper values 
   * according to ni, nj, and latbot, lattop, dlt, dln, etc.
   *
   * set_fmn():
   * compute Coriolis parameter and geometric map factors according to geometry.
   * grid.rln - long. map factor, r
   * grid.rlt - lat. map factor, R
   * grid.f   - Coriolis parameter
   * grid.m   - long. map factor 1/(r*dln*DEG) = 1./dx
   * grid.n   - lat. map factor 1/(R*dlt*DEG)  = 1./dy
   * grid.mn  - 1/(grid-box area)
   * grid.dy0 - dy at lat. = LAT0
   *
   * set_gravity():
   * Compute gravity as a function of latitude on the reference surface.
   */
  set_lonlat();
  set_fmn(planet);
  set_gravity(planet);

  /* Call set_dsgth() below after sigmatheta is defined. */

  /* 
   * Allocate memory: 
   */
  thetadat     = fvector( 0,var.ntp-1,dbmsname);
  neglogpdat   = fvector( 0,var.ntp-1,dbmsname);
  buff_triplet = ftriplet(0,var.ntp-1,dbmsname);

  /*
   * Interpolate on -log p:
   */
  for (ki = 0; ki < var.ntp; ki++) {
    neglogpdat[ki] = -log(var.pdat[ki]);
  }

  for (ki = 0; ki < var.ntp; ki++) {
    if (var.fpara.on) {
      fpara = return_fpe(var.tdat[ki]);
    }
    else {
      fpara = .25;
    }
    thetadat[ki] = return_theta(planet,fpara,var.pdat[ki],var.tdat[ki],
                                &theta_ortho,&theta_para);
  }

  /************************************************
   *                                              *
   * Set the pressure at the bottom of the model. *
   *                                              *
   ************************************************/

  if (strcmp(planet->type,"terrestrial") == 0) {
    /*
     * Inquire about artificial parameter grid.topo_scale.
     */
    if (strcmp(system_id,"earth") == 0) {
      defaults.topo_scale = grid.topo_scale = 1.;
    }
    else if (strcmp(system_id,"held_suarez") == 0) {
      defaults.topo_scale = grid.topo_scale = 1.;
    }
    else {
      defaults.topo_scale = grid.topo_scale = 
        input_float("Input artificial topo_scale factor (1.0 = no effect)\n",defaults.topo_scale);
    }
  }
  else {
    defaults.topo_scale = grid.topo_scale = 1.;
  }

  /*
   * Set the J index (latitude) to apply T(p) sounding data. This position is also used
   * to set the zero for gz_surface for gas giants.
   */
  if (strcmp(planet->name,"venus") == 0) {
    defaults.lat_tp = 4.;
  }
  else if (strcmp(planet->name,"jupiter") == 0) {
    defaults.lat_tp = input_float("Latitude to apply T(p) sounding data [deg]\n",defaults.lat_tp);
  }
  else {
    defaults.lat_tp = 0.;
  }
  /* Limit to valid range */
  defaults.lat_tp = LIMIT_RANGE(grid.globe_latbot,defaults.lat_tp,grid.globe_lattop);

  grid.jtp = 0;
  for (J = JLO; J <= JHI; J++) {
    if (fabs(grid.lat[2*J+1]-defaults.lat_tp) <= .5*grid.dlt) {
      grid.jtp = J;
      break;
    }
  }

  if (strcmp(planet->type,"terrestrial") == 0) {
    /*
     * For terrestrial planets, set the surface geopotential, and then 
     * find corresponding surface pressure.
     */
    init_gz_surface(planet,&defaults);

    /*
     * Set P for K=nk.
     * Reset defaults.pbot, grid.pbot to maximum surface pressure.
     */
    K = grid.nk;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        P3(K,J,I) = p_gz(planet,J,grid.topo_scale*GZ_SURFACE(J,I));
      }
    }
    /* No need to apply bc_lateral() here. */

    defaults.pbot = -FLOAT_MAX;
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        defaults.pbot = MAX(P3(K,J,I),defaults.pbot);
      }
    }
    grid.pbot = defaults.pbot;
    /* No need to apply bc_lateral() here. */
  }
  else if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * For gas giants, let the bottom of the model be a constant
     * pressure surface, and then find the corresponding surface geopotential
     * via gradient balance with the bottom zonal wind or zero.
     *
     * NOTE: We may want to use a different bottom boundary condition
     *       for models confined to the stratosphere, perhaps a constant-theta
     *       bottom condition to preserve the advantages of isentropic coordinates.
     */
    K = grid.nk;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        P3(K,J,I) = grid.pbot;
      }
    }
    /* No need to apply bc_lateral() here. */

    init_gz_surface(planet,&defaults);
  }
  else {
    sprintf(Message,"unrecognized planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }

  /*
   * Prompt for the target pressure for the top of the pure-sigma portion of the model.
   */
  if (defaults.p_sigma < 0.) {
    switch(planet->index) {
      case VENUS_INDEX:
      case VENUS_LLR05_INDEX:
        defaults.p_sigma = 500.*100.;
      break;
      case EARTH_INDEX:
      case HELD_SUAREZ_INDEX:
        defaults.p_sigma = 900.*100.;
      break;
      case MARS_INDEX:
        defaults.p_sigma = 0.2*100.;
      break;
      case JUPITER_INDEX:
        defaults.p_sigma = 680.*100.;
      break;
      case SATURN_INDEX:
        defaults.p_sigma = 500.*100.;
      break;
      case TITAN_INDEX:
        defaults.p_sigma = 100.*100.;
      break;
      case URANUS_INDEX:
        defaults.p_sigma = 1000.*100.;
      break;
      case NEPTUNE_INDEX:
        defaults.p_sigma = 1000.*100.;
      break;
      case HOT_JUPITER_INDEX:
        defaults.p_sigma = 680.*100.;  /*NOTE: Jupiter value used as a placeholder. */
      break;
      default:
        defaults.p_sigma = sqrt(defaults.pbot*defaults.ptop);
      break;
    }
  }
  defaults.p_sigma = 100.*input_float("Target pressure for top of pure-sigma vertical coordinate [hPa]\n",
                                      defaults.p_sigma/100.);

  /*
   * Determine t_vs_p data indices corresponding to pbot, ptop.
   */
  floor_tp = 0;
  tmp      = FLOAT_MAX;
  for (ki = 0; ki < var.ntp; ki++) {
    tmp2 = fabs(var.pdat[ki]-defaults.pbot);
    /* 
     * Need var.pdat[floor_tp] <= pbot.
     */
    if (tmp2 <= tmp && var.pdat[ki] <= defaults.pbot) {
      floor_tp = ki;
      tmp      = tmp2;
    }
  }
  ceiling_tp = var.ntp-1;
  tmp        = FLOAT_MAX;
  for (ki = var.ntp-1; ki >= 0; ki--) {
    tmp2 = fabs(var.pdat[ki]-defaults.ptop);
    if (tmp2 < tmp && var.pdat[ki] <= defaults.ptop) {
      ceiling_tp = ki;
      tmp        = tmp2;
    }
  }

  /*
   * Set grid.zeta0, grid.zeta1, grid.hybrid_alpha
   *
   * NOTE: These values may need to be adjusted to yield a hybrid coordinate, zeta,
   *       that is monotonic and reasonable.
   *
   * These parameters control the behavior of the function
   * f_sigma() used to define the hybrid coordinate, zeta.
   *
   * The value of the coordinate at the bottom of the model atmosphere, sigma = 0., is grid.zeta0.
   * This is a free parameter that should be chosen such that zeta is monotonic and smooth.
   */
  switch(planet->index) {
    case VENUS_INDEX:
      /*
       * Grace Lee 3/22/05
       *
       * Recommended p_sigma = 500 hPa (for pbot > 500 hPa).
       */
      grid.zeta0        = 420.-53.333*log10(defaults.pbot/100.);
      grid.zeta1        = 420.-53.333*log10(defaults.ptop/100.);
      grid.hybrid_alpha = 50.;
    break;
    case VENUS_LLR05_INDEX:
      /*
       * Aaron Herrnstein 6/26/06
       */
      grid.zeta0        = 165.;
      grid.zeta1        = 480.;
      grid.hybrid_alpha = 20.;
    break;
    case EARTH_INDEX:
    case HELD_SUAREZ_INDEX:
      /*
       * Recommended p_sigma = 900 hPa.
       */
      grid.zeta0        = 260.;
      grid.zeta1        = 700.;
      grid.hybrid_alpha = 50.;
    break;
    case MARS_INDEX:
      /*
       * TD 5/4/09
       * Recommended p_sigma = 0.2 hPa.
       */
      grid.zeta0        = 700.-300.*log10(defaults.pbot/100.);
      grid.zeta1        = 700.-300.*log10(defaults.ptop/100.);
      grid.hybrid_alpha = 50.;
    break;
    case JUPITER_INDEX:
      /*
       * CJP 2/9/05
       * The tangent line is written in the form of
       *
       *     zeta = 412.5 - 87.5*log10(p)
       *
       * Recommended p_sigma = 680 hPa (for pbot > 680 hPa).
       *
       * grid.zeta0 and grid.zeta1 are being determined from this equation
       * substituting p_bot and p_top respectively
       */
      grid.zeta0        = 412.5-87.5*log10(defaults.pbot/100.);
      grid.zeta1        = 412.5-87.5*log10(defaults.ptop/100.);
      grid.hybrid_alpha = 50.;
      /*
       * NOTE: A check that zeta does not cross the theta curve might be useful.
       */
    break;
    case SATURN_INDEX:
      /*
       * Csaba Palotai 4/23/07
       * The tangent line is written in the form of
       *
       *     zeta = 230.0 - 35.0*log10(p)
       *
       * Recommended p_sigma = 500 hPa (for pbot > 500 hPa).
       *
       * grid.zeta0 and grid.zeta1 are being determined from this equation
       * substituting p_bot and p_top respectively
       */
      grid.zeta0        = 230.0-35.0*log10(defaults.pbot/100.);
      grid.zeta1        = 230.0-35.0*log10(defaults.ptop/100.);
      grid.hybrid_alpha = 50.;
      /*
       * NOTE: A check that zeta does not cross the theta curve might be useful.
       */
    break;
    case TITAN_INDEX:
       /*
        * Recommended p_sigma = 100 hPa.
        */
      grid.zeta0        = 30.;
      grid.zeta1        = grid.zeta0+(133.-grid.zeta0)*log(defaults.ptop/defaults.pbot)/log(10000./defaults.pbot);
      grid.hybrid_alpha = 50.;
    break;
    case URANUS_INDEX:
       /*
        * Recommended p_sigma = 100 hPa.
        */
      grid.zeta0        = 280.5-83.915*log10(defaults.pbot/100.);
      grid.zeta1        = 280.5-83.915*log10(defaults.ptop/100.);
      grid.hybrid_alpha = 50.;
    break;
    case NEPTUNE_INDEX:
       /*
        * Recommended p_sigma = 1000 hPa.
        */
      grid.zeta0        = 138.7-23.2*log10(defaults.pbot/100.);
      grid.zeta1        = 138.7-23.2*log10(defaults.ptop/100.);
      grid.hybrid_alpha = 50.;
    break;
    case HOT_JUPITER_INDEX:     /* NOTE: Jupiter values used as placeholders */
      /*
       * CJP 2/9/05
       * The tangent line is written in the form of
       *
       *     zeta = 412.5 - 87.5*log10(p)
       *
       * Recommended p_sigma = 680 hPa (for pbot > 680 hPa).
       *
       * grid.zeta0 and grid.zeta1 are being determined from this equation
       * substituting p_bot and p_top respectively
       */
      grid.zeta0        = 412.5-87.5*log10(defaults.pbot/100.);
      grid.zeta1        = 412.5-87.5*log10(defaults.ptop/100.);
      grid.hybrid_alpha = 50.;
      /*
       * NOTE: A check that zeta does not cross the theta curve might be useful.
       */
    break;
    default:
      tmp = FLOAT_MAX;
      for (ki = floor_tp; ki <= ceiling_tp; ki++) {
        tmp = MIN(tmp,thetadat[ki]);
      }
      grid.zeta0 = tmp*0.7;

      /*
       * Find the highest tangent point from grid.zeta0 to the theta profile,
       * call it theta_knee, sigma_knee.
       */
      tmp = FLOAT_MAX;
      for (ki = floor_tp; ki <= ceiling_tp; ki++) {
        if (var.pdat[ki] < defaults.ptop) {
          break;
        }
        tmp2 = fabs(var.pdat[ki]-defaults.ptop);
        if (tmp2 < tmp) {
          tmp        = tmp2;
          theta_knee = thetadat[ki];
          sigma_knee = get_sigma(defaults.pbot,var.pdat[ki],defaults.ptop);
        }
      }
      for (ki = ceiling_tp; ki > floor_tp; ki--) {
        if (var.pdat[ki] < defaults.ptop) {
          continue;
        }
        sigma = get_sigma(defaults.pbot,var.pdat[ki],defaults.ptop);
        slope = (thetadat[ki]-grid.zeta0)/(sigma-0.);
        tmp   = get_sigma(defaults.pbot,var.pdat[ki+1],defaults.ptop);
        tmp2  = grid.zeta0+slope*tmp;
        if (tmp2 <= thetadat[ki+1]) {
          if (thetadat[ki] < theta_knee) {
            theta_knee = thetadat[ki];
            sigma_knee = sigma;
          }
        }
        else {
          break;
        }
      }
      if (sigma_knee > 0.) {
        grid.zeta1  = (theta_knee-grid.zeta0*(1.-sigma_knee))/sigma_knee;
      }
      else {
        sprintf(Message,"sigma_knee=%g",sigma_knee);
        epic_error(dbmsname,Message);
      }
      grid.hybrid_alpha = 50.;
    break;
  }
  fprintf(stdout,"hybrid-coordinate parameters: grid.zeta0=%g, grid.zeta1=%g, grid.hybrid_alpha=%g\n",
                  grid.zeta0,grid.zeta1,grid.hybrid_alpha);
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   *                                                         *
   * Set grid.p_ref[kk], grid.theta_ref[kk], and             *
   * the vertical coordinate grid.sigmatheta[kk].            *
   *                                                         *
   * grid.p_ref[kk] - typical pressure values                *
   * grid.theta_ref[kk] - typical theta values               *
   * grid.sigmatheta[kk] - hybrid sigmatheta array           *
   * grid.sgth_bot - sgth at bottom of model                 *
   * grid.sgth_top - sgth at top of model                    *
   *                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  /*
   * Top and bottom values:
   */
  grid.p_ref[      1] = defaults.ptop;
  grid.p_ref[2*KHI+1] = defaults.pbot;

  for (ki = 0; ki < var.ntp; ki++) {
    buff_triplet[ki].x = -log(var.pdat[ki]);
    buff_triplet[ki].y = thetadat[ki];
  }
  spline_pchip(var.ntp,buff_triplet);

  neglogp           = -log(defaults.ptop);
  ki                = find_place_in_table(var.ntp,buff_triplet,neglogp,&p_d);
  grid.theta_ref[1] = splint_pchip(neglogp,buff_triplet+ki,p_d);

  neglogp                 = -log(defaults.pbot);
  ki                      = find_place_in_table(var.ntp,buff_triplet,neglogp,&p_d);
  grid.theta_ref[2*KHI+1] = splint_pchip(neglogp,buff_triplet+ki,p_d);

  grid.sgth_top = grid.theta_ref[1];
  grid.sgth_bot = return_sigmatheta(grid.theta_ref[2*KHI+1],defaults.pbot,defaults.pbot,defaults.ptop);

  grid.sigmatheta[0] = grid.sigmatheta[1] = grid.sgth_top;

  if (strcmp(planet->type,"gas-giant") == 0) {
    grid.p_ref[    2*(KHI+1)]  = grid.p_ref[    2*KHI+1];
    grid.theta_ref[2*(KHI+1)]  = grid.theta_ref[2*KHI+1];
    grid.sigmatheta[2*KHI+1  ] = grid.sgth_bot;
    grid.sigmatheta[2*(KHI+1)] = grid.sgth_bot;
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    grid.sigmatheta[2*KHI+1] = grid.sgth_bot;
  }
  else {
    sprintf(Message,"unrecognized planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }

  fpara      = .25;
  grid.mont0 = return_enthalpy(planet,fpara,grid.p_ref[2*KHI+1],var.tdat[floor_tp],&fgibb,&fpe,&uoup);

  /* 
   * Set grid.p_ref[kk], grid.theta_ref[kk], grid.sigmatheta[kk] for interior of model.
   */
  if (spacing_type == SPACING_P) {
    for (kk = 2; kk < 2*KHI+1; kk++) {
      grid.p_ref[kk] = defaults.ptop+(EPIC_FLOAT)(kk-1)/(EPIC_FLOAT)(2*KHI)*(defaults.pbot-defaults.ptop);
    }
  }
  else if (spacing_type == SPACING_LOGP_W_BOUNDARIES && ALPHA_BOUNDARY != 0.) {
    if (strcmp(planet->type,"terrestrial") == 0) {
      nk_sponge = grid.k_sponge;
      /*
       * Specify number of layers in planetary boundary layer with proportional spacing, nk_pbl.
       */
      nk_pbl    = grid.nk/5;
      nk_reg    = grid.nk-nk_pbl-nk_sponge;
      /*
       * Calculate regular spacing, dlnp.
       */
      dlnp = log(defaults.pbot/defaults.ptop)/
             ((1.-pow(1.+ALPHA_BOUNDARY,-nk_sponge))/ALPHA_BOUNDARY
             +(EPIC_FLOAT)nk_reg
             +(1.-pow(1.+ALPHA_BOUNDARY,-nk_pbl))/ALPHA_BOUNDARY);

      grid.p_ref[2*KLO+1] = grid.p_ref[1]*exp(dlnp*pow(1.+ALPHA_BOUNDARY,-nk_sponge));
      for (K = KLO+1; K <= grid.k_sponge; K++) {
        /*
         * Geometrically decreasing spacing in the top sponge.
         */
        grid.p_ref[2*K+1] = grid.p_ref[2*(K-1)+1]*pow(grid.p_ref[2*(K-1)+1]/grid.p_ref[2*(K-2)+1],1.+ALPHA_BOUNDARY);
      }
      for (K = grid.k_sponge+1; K <= grid.nk-nk_pbl; K++) {
        /* 
         * Regular spacing.
         */
        grid.p_ref[2*K+1] = grid.p_ref[2*(K-1)+1]*exp(dlnp);
      }
      for (K = grid.nk-nk_pbl+1; K < grid.nk; K++) {
        /*
         * Geometrically decreasing spacing in the planetary boundary layer (PBL).
         */
        grid.p_ref[2*K+1] = grid.p_ref[2*(K-1)+1]*pow(grid.p_ref[2*(K-1)+1]/grid.p_ref[2*(K-2)+1],1./(1.+ALPHA_BOUNDARY));
      }
      /*
       * Fill in layer values.
       */
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K;
        grid.p_ref[kk] = onto_kk(planet,P2_INDEX,grid.p_ref[kk-1],grid.p_ref[kk+1],kk,JLO,ILO);
      }
    }
    else if (strcmp(planet->type,"gas-giant") == 0) {
      nk_sponge = grid.k_sponge;
      nk_reg    = grid.nk-nk_sponge;
      /*
       * Calculate regular spacing, dlnp.
       */
      dlnp      = log(defaults.pbot/defaults.ptop)/
                  ((1.-pow(1.+ALPHA_BOUNDARY,-nk_sponge))/ALPHA_BOUNDARY
                  +(EPIC_FLOAT)nk_reg);

      grid.p_ref[2*KLO+1] = grid.p_ref[1]*exp(dlnp*pow(1.+ALPHA_BOUNDARY,-nk_sponge));
      for (K = KLO+1; K <= grid.k_sponge; K++) {
        /*
         * Geometrically decreasing spacing in the top sponge.
         */
        grid.p_ref[2*K+1] = grid.p_ref[2*(K-1)+1]*pow(grid.p_ref[2*(K-1)+1]/grid.p_ref[2*(K-2)+1],1.+ALPHA_BOUNDARY);
      }
      for (K = grid.k_sponge+1; K <= grid.nk; K++) {
        /* 
         * Regular spacing.
         */
        grid.p_ref[2*K+1] = grid.p_ref[2*(K-1)+1]*exp(dlnp);
      }
      /*
       * Fill in layer values.
       */
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K;
        grid.p_ref[kk] = onto_kk(planet,P2_INDEX,grid.p_ref[kk-1],grid.p_ref[kk+1],kk,JLO,ILO);
      }
    }
    else {
      sprintf(Message,"unrecognized planet->type=%s",planet->type);
      epic_error(dbmsname,Message);
    }
  }
  else if (spacing_type == SPACING_LOGP || (spacing_type == SPACING_LOGP_W_BOUNDARIES && ALPHA_BOUNDARY == 0.)) {
    for (K = KLO; K < KHI; K++) {
      neglogp        = -log(defaults.ptop)
                       +(EPIC_FLOAT)(K)/(EPIC_FLOAT)(KHI)*(-log(defaults.pbot)+log(defaults.ptop));
      grid.p_ref[2*K+1] = exp(-neglogp);
    }
    /*
     * Fill in layer values.
     */
    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;
      grid.p_ref[kk] = onto_kk(planet,P2_INDEX,grid.p_ref[kk-1],grid.p_ref[kk+1],kk,JLO,ILO);
    }
  }
  else if (spacing_type == SPACING_FROM_FILE) {
    /*
     * Have grid.p_ref[K] from a file.
     */
    ;
  }
  else {
    sprintf(Message,"unrecognized spacing_type=%d",spacing_type);
    epic_error(dbmsname,Message);
  }

  /*
   * Set grid.k_sigma.
   * The model's vertical coordinate is a pure sigma coordinate (with
   * units of Kelvins) for K = KHI up to the top of layer K = grid.k_sigma.
   */
  grid.k_sigma = -1;
  tmp2       = FLOAT_MAX;
  for (K = KLO; K <= KHI; K++) {
    tmp = fabs(grid.p_ref[2*(K-1)+1]-defaults.p_sigma);
    if (tmp < tmp2) {
      grid.k_sigma = K;
      tmp2       = tmp;
    }
  }
  if (grid.k_sigma == -1) {
    sprintf(Message,"Error calculating grid.k_sigma with p_sigma=%g hPa",
                    defaults.p_sigma/100.);
    epic_error(dbmsname,Message);
  }
  else if (grid.k_sigma < 2) {
    sprintf(Message,"Need grid.k_sigma=%d >= 2",grid.k_sigma);
    epic_error(dbmsname,Message);
  }

  /*
   * Set corresponding grid.sigma_sigma, which is the value of sigma
   * at the top of layer K = grid.k_sigma.
   *
   * NOTE: grid.sigmatheta[kk] is not yet set at this point.
   */
  grid.sigma_sigma = get_sigma(grid.p_ref[2*KHI+1],grid.p_ref[2*(grid.k_sigma-1)+1],grid.p_ref[1]);

  /*
   * Set grid.theta_ref[K].
   */
  if (spacing_type == SPACING_P) {
    for (ki = 0; ki < var.ntp; ki++) {
      buff_triplet[ki].x = -var.pdat[ki];
      buff_triplet[ki].y = thetadat[ki];
    }
    spline_pchip(var.ntp,buff_triplet);
    for (kk = 2; kk < 2*KHI+1; kk++) {
      negp               = -grid.p_ref[kk];
      ki                 = find_place_in_table(var.ntp,buff_triplet,negp,&p_d);
      grid.theta_ref[kk] = splint_pchip(negp,buff_triplet+ki,p_d);
    }
  }
  else {
    for (ki = 0; ki < var.ntp; ki++) {
      buff_triplet[ki].x = neglogpdat[ki];
      buff_triplet[ki].y = thetadat[ki];
    }
    spline_pchip(var.ntp,buff_triplet);
    for (kk = 2; kk < 2*KHI+1; kk++) {
      negp               = -log(grid.p_ref[kk]);
      ki                 = find_place_in_table(var.ntp,buff_triplet,negp,&p_d);
      grid.theta_ref[kk] = splint_pchip(negp,buff_triplet+ki,p_d);
    }
  }

  /*
   * Set grid.sigmatheta[kk] for interior points.
   */
  ptop = grid.p_ref[      1];
  pbot = grid.p_ref[2*KHI+1];
  for (kk = 2; kk < 2*KHI+1; kk++) {
    pressure            = grid.p_ref[kk];
    theta               = grid.theta_ref[kk];
    grid.sigmatheta[kk] = return_sigmatheta(theta,pressure,pbot,ptop);
  }

  /*
   * Now, reset grid.sigma_sigma to reduce roundoff error.
   */
  grid.sigma_sigma = (grid.sigmatheta[2*grid.k_sigma-1]-grid.zeta0)/(grid.zeta1-grid.zeta0);

  /*
   * Set grid.dsgth and grid.dsgth_inv.
   */
  set_dsgth();

  /*
   * Set grid.h_min[K], the minimum thickness (a.k.a. hybrid density) for
   * each layer. Konor and Arakawa use delta p_min = 1 hPa to find h_min, 
   * but this will not work for planets that have small atmospheric pressures
   * (for example, the surface pressure of Mars is 6 hPa). 
   * Instead, we set h_min equal to a fraction of the reference value of h in each layer.
   */  
  tmp = 0.02/grid.g[2*grid.jtp+1];
  for (K = KLO; K <= KHI; K++) {
    kk            = 2*K;
    grid.h_min[K] = tmp*(grid.p_ref[2*K+1]-grid.p_ref[2*(K-1)+1])*grid.dsgth_inv[kk];
  }
  grid.h_min[0] = grid.h_min[KLO];

  if (var.fpara.on) {
   /*
    * Preliminary initialization of para-hydrogen fraction, fpara.
    * Set to fpe(T), for representative T for layer.
    *
    * NOTE: The function init_fpara_as_fpe() cannot be used here because
    *       P2 and THETA2 are not known yet.
    */
    for (K = KLO; K <= KHI; K++) {
      pressure = grid.p_ref[2*K+1];
      get_sounding(planet,pressure,"temperature",&temperature);
      fpe = return_fpe(temperature);
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FPARA(K,J,I) = fpe;
        }
      }
    }
    /* Assume fpara doesn't change with height above top of model. */
    K = 0;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        FPARA(K,J,I) = FPARA(K+1,J,I);
      }
    }
    bc_lateral(var.fpara.value,THREEDIM);
  }

  /*
   * Initialize u,v,hdry,theta, and mixing ratios of optional species.
   *
   * Start with wind-free case to establish an approximate molar_mass vs pressure function
   * that includes the contributions from condensables, such that setup_mu_p() will work.
   */
  if (grid.nmt_physics_on) {
    /* set the number of cumulus cells to the number of EPIC layers */
    nmt.nz = grid.nk;

    sprintf(out_str,"Do you want this model over land? [1 = land, 0 = water]\n");
    defaults.nmt_land = nmt.land = input_int(out_str,defaults.nmt_land);

    sprintf(out_str,"Radiation routine on or off? [1 = on, 0 = off]\n");
    defaults.nmt_dorad = nmt.dorad = input_int(out_str,defaults.nmt_dorad);

    sprintf(out_str,"Fixed radiative cooling? [1 = yes, 0 = no]\n");
    defaults.nmt_frad = nmt.frad = input_int(out_str,defaults.nmt_frad);

    sprintf(out_str,"Radiative cooling rate? (K/day)\n");
    defaults.nmt_radcool = nmt.radcool = input_float(out_str,defaults.nmt_radcool);

    sprintf(out_str,"Height of the tropopause? (km)\n");
    defaults.nmt_tpause = nmt.tpause = input_float(out_str,defaults.nmt_tpause);

    sprintf(out_str,"Fraction of tropopause where fixed radiative cooling begins to decrease?\n");
    defaults.nmt_radbrk = nmt.radbrk = input_float(out_str,defaults.nmt_radbrk);

    sprintf(out_str,"Convective mixing rate parameter? (ks^-1)\n");
    defaults.nmt_cvc = nmt.cvc = input_float(out_str,defaults.nmt_cvc);

    sprintf(out_str,"Stratiform rain rate parameter? (ks^-1)\n");
    defaults.nmt_cvs = nmt.cvs = input_float(out_str,defaults.nmt_cvs);

    sprintf(out_str,"Convective rain rate parameter? (ks^-1)\n");
    defaults.nmt_cvp = nmt.cvp = input_float(out_str,defaults.nmt_cvp);

    sprintf(out_str,"Rain evaporation rate parameter? (ks^-1)\n");
    defaults.nmt_cve = nmt.cve = input_float(out_str,defaults.nmt_cve);

    sprintf(out_str,"Range of theta_e going from suppressed to full deep convection? (K)\n");
    defaults.nmt_theslop = nmt.theslop = input_float(out_str,defaults.nmt_theslop);

    sprintf(out_str,"Stiffness parameter in convective precipitation generation?\n");
    defaults.nmt_pstiff = nmt.pstiff = input_float(out_str,defaults.nmt_pstiff);

    sprintf(out_str,"Shape parameter for convective precipitation production? (km)\n");
    defaults.nmt_pscale = nmt.pscale = input_float(out_str,defaults.nmt_pscale);

    sprintf(out_str,"Top of the planetary boundary layer? (km)\n");
    defaults.nmt_pbltop = nmt.pbltop = input_float(out_str,defaults.nmt_pbltop);

    sprintf(out_str,"Surface drag / thermodynamics coefficient? \n");
    defaults.nmt_cdrag = nmt.cdrag = input_float(out_str,defaults.nmt_cdrag);

    sprintf(out_str,"Gustiness parameter for surface fluxes? (m/s)\n");
    defaults.nmt_wscale = nmt.wscale = input_float(out_str,defaults.nmt_wscale);

    sprintf(out_str,"Sea surface temperature? (K)\n");
    defaults.nmt_sst = nmt.sst = input_float(out_str,defaults.nmt_sst);

    sprintf(out_str,"Latent heat fraction out of total over land?\n");
    defaults.nmt_lfrac = nmt.lfrac = input_float(out_str,defaults.nmt_lfrac);

    sprintf(out_str,"Theta_e flux over land?\n");
    defaults.nmt_eflux0 = nmt.eflux0 = input_float(out_str,defaults.nmt_eflux0);

    sprintf(out_str,"Cloud radiative absorption parameter? (m^3/kg/km)\n");
    defaults.nmt_cld = nmt.cld = input_float(out_str,defaults.nmt_cld);

    sprintf(out_str,"Cloud fractional area coverage?\n");
    defaults.nmt_cfract = nmt.cfract = input_float(out_str,defaults.nmt_cfract);
  }

  /*
   * Use reference profiles to estimate mu(p), in case it is needed.
   */
  init_with_ref(planet);
  init_species(planet,&defaults,USE_PROMPTS);
  setup_mu_p(planet);

  if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * NOTE: If gradient-balance in the meridional plane is not desired (init_with_u())
     *       then add an appropriate initialization scheme.
     */
    init_with_u(planet,&defaults);
  }
  else if (strcmp(planet->name,"held_suarez") == 0) {
    init_with_u(planet,&defaults);
  }
  else {
    if (fcmp(defaults.u_scale,   0.) == 0 &&
        fcmp(defaults.topo_scale,0.) != 0   ) {
      init_with_hydrostatic(planet);
    }
    else {
      /*
       * NOTE: If gradient-balance in the meridional plane is not desired (init_with_u())
       *       then add an appropriate initialization scheme.
       */
      init_with_u(planet,&defaults);
    }
  }

  init_species(planet,&defaults,USE_DEFAULTS);

  /*
   * Update fpara as fpe:
   */
  if (var.fpara.on) {
    init_fpara_as_fpe(planet);
  }

  /* 
   * Initialize hyperviscosity coefficients.
   */
  grid.nudiv_nondim = defaults.nudiv_nondim;
  for (ii = 0; ii <= MAX_NU_ORDER; ii++) {
    grid.nu_nondim[ii] = defaults.nu_nondim[ii];
  }
  init_viscosity(planet);
  defaults.nudiv_nondim = grid.nudiv_nondim;
  for (ii = 0; ii <= MAX_NU_ORDER; ii++) {
    defaults.nu_nondim[ii] = grid.nu_nondim[ii];
  }

  /*
   * Initialize turbulence-model variables.
   */
  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
    /*
     * NOTE: init_subgrid() must be called here before synching diagnostic variables.
     */
    init_subgrid(planet);
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  /*
   * Calculate most commonly used diagnostic variables, in case they are needed below.
   */
  set_p2etc(planet,UPDATE_THETA);
  store_pgrad_vars(planet);
  store_diag(planet);

  if (grid.nmt_physics_on == 1) {
   /*
    * Initialize total cloud water for nmt_physics
    */
    fprintf(stdout,"Initialize NMT Physics...");

    nmt_init_water();

    fprintf(stdout,"done.\n");
  }

  /*
   * set_u_spinup() assumes U, defaults.u_spinup_scale, and grid.k_sponge are set.
   */
  set_u_spinup(planet,&defaults);

  /*
   * Prompt for which variables to write to extract.nc.
   */
  prompt_extract_on(defaults.extract_str);

  /*
   * Write defaults file:
   */
  write_defaults(&defaults);

  /* 
   * Print out zonal-wind information: 
   */
  print_zonal_info(planet);

  /* 
   * Print out vertical information: 
   */
  print_vertical_column(planet,JLO,ILO,"vertical.dat");

  /*
   * Call vertical_modes() to write vertical eigenvalues and eigenvectors to a file.
   */
  if (KHI >= 4) {
    vertical_modes(planet,JLO,ILO);
  }

  /* 
   * Output epic.nc file (netCDF = Network Common Data Form):
   */
  sprintf(outfile,"epic.nc");
  time_index = 0;
  var_write(planet,outfile,ALL_DATA,time_index,0);

  /* 
   * Free allocated memory: 
   */
  free_arrays(planet); 
  free_ftriplet(buff_triplet,0,var.ntp-1,         dbmsname);
  free_fvector(neglogpdat,   0,var.ntp-1,         dbmsname);
  free_fvector(thetadat,     0,var.ntp-1,         dbmsname);
  free_var_props();
 
  return 0;
}

/*======================= end of main() =====================================*/

/*======================= read_defaults() ===================================*/

/*
 * NOTE: To use the READ* macros in epic_io_macros.h, dummy
 *       io and fd variables are declared.
 */

void read_defaults(init_defaultspec *def) 
{
  int
    index,
    nc_id,nc_err;
  EPIC_FLOAT
    solar;
  char
    min_element[4];
  static char
    **gattname=NULL,
    **varname =NULL;
  static int
    ngatts   =0,
    num_vars =0;
  nc_type
    the_nc_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_defaults";

  nc_err = lookup_netcdf("init_defaults.nc",
                         &nc_id,&ngatts,&gattname,&num_vars,&varname);

  if (nc_err == NC_NOERR) {
    READI(&def->start_date_input_type,def_start_date_input_type,1);
    READTIME(&def->start_time,def_start_time);  

    READC(def->geometry,def_geometry,GEOM_STR);
    READC(def->system_id,def_system_id,32);
    READC(def->f_plane_map,def_f_plane_map,GEOM_STR);
    READC(def->eos,def_eos,8);
    READC(def->extract_str,def_extract_str,N_STR);
    READC(def->species_str,def_species_str,N_STR);
    READC(def->layer_spacing_dat,def_layer_spacing_dat,N_STR);

    READI(&def->nk,def_nk,1);
    READI(&def->nj,def_nj,1);
    READI(&def->ni,def_ni,1);
    READI(&def->dt,def_dt,1);
    READI(&def->newt_cool_on,def_newt_cool_on,1);
    READI(&def->newt_cool_adjust,def_newt_cool_adjust,1);
    READI(&def->microphysics_on,def_microphysics_on,1);
    READI(&def->nmt_physics_on,def_nmt_physics_on,1);
    READI(&def->diffusion_direction,def_diffusion_direction,1);

    READI(def->on,def_on,LAST_SPECIES+1);

    READF(&def->du_vert,def_du_vert,1);
    READI(&def->hgrid_mass_advection_scheme,def_hgrid_mass_advection_scheme,1);
    READI(&def->hgrid_nonmass_advection_scheme,def_hgrid_nonmass_advection_scheme,1);
    READI(&def->uv_timestep_scheme,def_uv_timestep_scheme,1);
    READI(&def->turbulence_scheme,def_turbulence_scheme,1);
    READI(&def->stability_factor,def_stability_factor,1);
    READI(&def->spacing_type,def_spacing_type,1);
    READI(&def->drag_v,def_drag_v,1);
    READI(&def->drag_zonal_avg,def_drag_zonal_avg,1);
    READI(&def->k_sponge,def_k_sponge,1);

    READF(&def->globe_lonbot,def_globe_lonbot,1);
    READF(&def->globe_lontop,def_globe_lontop,1);
    READF(&def->globe_latbot,def_globe_latbot,1);
    READF(&def->globe_lattop,def_globe_lattop,1);
    READF(&def->lat_tp,def_lat_tp,1);
    READF(&def->f_plane_half_width,def_f_plane_half_width,1);
    READF(&def->ptop,def_ptop,1);
    READF(&def->pbot,def_pbot,1);
    READF(&def->p_sigma,def_p_sigma,1);
    READF(&def->prandtl,def_prandtl,1);
    READF(&def->tau_drag,def_tau_drag,1);
    READF(&def->nudiv_nondim,def_nudiv_nondim,1);
    READF(&def->u_scale,def_u_scale,1);
    READF(&def->u_spinup_scale,def_u_spinup_scale,1);
    READF(&def->topo_scale,def_topo_scale,1);
    READF(def->nu_nondim,def_nu_nondim,MAX_NU_ORDER+1);
    READF(&def->time_fp,def_time_fp,1);

    READF(def->mole_fraction,def_mole_fraction,LAST_SPECIES+1);
    READF(def->mole_fraction_over_solar,def_mole_fraction_over_solar,LAST_SPECIES+1);
    READF(def->rh_max,def_rh_max,LAST_SPECIES+1);

    /*
     * nmt_physics
     */
    READI(&def->nmt_land,def_nmt_land,1);
    READI(&def->nmt_dorad,def_nmt_dorad,1);
    READI(&def->nmt_frad,def_nmt_frad,1);
    READF(&def->nmt_radcool,def_nmt_radcool,1);
    READF(&def->nmt_tpause,def_nmt_tpause,1);
    READF(&def->nmt_radbrk,def_nmt_radbrk,1);
    READF(&def->nmt_cvc,def_nmt_cvc,1);
    READF(&def->nmt_cvs,def_nmt_cvs,1);
    READF(&def->nmt_cvp,def_nmt_cvp,1);
    READF(&def->nmt_cve,def_nmt_cve,1);
    READF(&def->nmt_theslop,def_nmt_theslop,1);
    READF(&def->nmt_pstiff,def_nmt_pstiff,1);
    READF(&def->nmt_pscale,def_nmt_pscale,1);
    READF(&def->nmt_pbltop,def_nmt_pbltop,1);
    READF(&def->nmt_cdrag,def_nmt_cdrag,1);
    READF(&def->nmt_wscale,def_nmt_wscale,1);
    READF(&def->nmt_sst,def_nmt_sst,1);
    READF(&def->nmt_lfrac,def_nmt_lfrac,1);
    READF(&def->nmt_eflux0,def_nmt_eflux0,1);
    READF(&def->nmt_cld,def_nmt_cld,1);
    READF(&def->nmt_cfract,def_nmt_cfract,1);
  }
  else {
    /*
     * If the file is not readable, use standard defaults.
     */
    strcpy(def->geometry,"globe");
    strcpy(def->system_id,  "Jupiter");
    strcpy(def->f_plane_map,"polar");
    strcpy(def->eos,        "ideal");
    strcpy(def->extract_str,"none");
    strcpy(def->species_str,"none");
    strcpy(def->layer_spacing_dat,"layer_spacing.dat");

    def->start_date_input_type   = -1; /* toggles a smart prompt */
    def->start_time              = (time_t)(-INT_MAX)+1566847;  /* round up to January 1, 1902 */

    def->nk                      =  -1;  /* toggles a smart prompt */
    def->nj                      =  64;
    def->ni                      =  128;
    def->dt                      =  120;
    def->newt_cool_on            =  TRUE;
    def->newt_cool_adjust        =  FALSE;
    def->microphysics_on         =  FALSE;
    def->nmt_physics_on          =  FALSE;
    def->diffusion_direction     =  HORIZONTAL_AND_VERTICAL;

    def->on[    U_INDEX]    =  TRUE;
    def->on[    V_INDEX]    =  TRUE;
    def->on[ HDRY_INDEX]    =  TRUE;
    def->on[THETA_INDEX]    =  TRUE;
    for (index = THETA_INDEX+1; index <= LAST_SPECIES; index++) {
      def->on[index]   = FALSE;
    }

    def->du_vert                        = .0;
    def->hgrid_mass_advection_scheme    =  0;
    def->hgrid_nonmass_advection_scheme =  0;
    def->uv_timestep_scheme             =  0;
    def->turbulence_scheme              =  1;
    def->stability_factor               =  1;
    def->spacing_type                   =  SPACING_LOGP_W_BOUNDARIES;
    def->drag_v                         =  FALSE;
    def->drag_zonal_avg                 =  TRUE;
    def->k_sponge                       = -1;     /* toggles a smart prompt */

    def->globe_lonbot       = -180.;
    def->globe_lontop       =  180.;
    def->globe_latbot       = -90.;
    def->globe_lattop       =  90.;
    def->lat_tp             =  0.;
    def->f_plane_half_width =  90.;
    def->ptop               =  .1*100.;
    def->pbot               =  10000.*100.;
    def->p_sigma            = -1.;        /* toggles a smart prompt */
    def->prandtl            =  0.;
    def->tau_drag           =  1.e+20;
    def->nudiv_nondim       =  0.;
    def->u_scale            =  1.;
    def->u_spinup_scale     =  1.;
    def->topo_scale         =  1.;
    def->nu_nondim[2]       =  0.;
    def->nu_nondim[4]       =  0.;
    def->nu_nondim[6]       =  0.;
    def->nu_nondim[8]       =  0.5;
    def->time_fp            =  3.e+8;
    for (index = FIRST_SPECIES; index <= LAST_SPECIES; index++) {
      solar = solar_fraction(var.species[index].info[0].name,MOLAR,min_element);

      def->mole_fraction[index]            = solar;
      def->mole_fraction_over_solar[index] = 1.;
      def->rh_max[index]                   = 1.;
    }

    /*
     * nmt_physics
     */
    def->nmt_land           =  0;
    def->nmt_dorad          =  1;
    def->nmt_frad           =  0;
    def->nmt_radcool        =  2.;
    def->nmt_tpause         =  15.;
    def->nmt_radbrk         =  0.8;
    def->nmt_cvc            =  0.01;
    def->nmt_cvs            =  0.1;
    def->nmt_cvp            =  0.0004;
    def->nmt_cve            =  100.0;
    def->nmt_theslop        =  4.0;
    def->nmt_pstiff         =  4.0;
    def->nmt_pscale         =  6.0;
    def->nmt_pbltop         =  1.5;
    def->nmt_cdrag          =  0.001;
    def->nmt_wscale         =  3.0;
    def->nmt_sst            =  300.;
    def->nmt_lfrac          =  1.;
    def->nmt_eflux0         =  0.;
    def->nmt_cld            =  75000.;
    def->nmt_cfract         =  0.;
  }

  return;
}

/*======================= end of read_defaults() ============================*/

/*======================= write_defaults() ==================================*/

/*
 * NOTE: To use the WRITE* macros in epic_io_macros.h, dummy
 *       io and fd variables are declared.
 */

void write_defaults(init_defaultspec *def)
{
  int
    nc_id,nc_err;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="write_defaults";

  nc_err = nc_create("init_defaults.nc",NC_CLOBBER,&nc_id);

  WRITEI(&def->start_date_input_type,def_start_date_input_type,1);
  WRITETIME(&def->start_time,def_start_time);  

  WRITEC(def->geometry,def_geometry,GEOM_STR);
  WRITEC(def->system_id,def_system_id,32);
  WRITEC(def->f_plane_map,def_f_plane_map,GEOM_STR);
  WRITEC(def->eos,def_eos,8);
  WRITEC(def->extract_str,def_extract_str,N_STR);
  WRITEC(def->species_str,def_species_str,N_STR);
  WRITEC(def->layer_spacing_dat,def_layer_spacing_dat,N_STR);

  WRITEI(&def->nk,def_nk,1);
  WRITEI(&def->nj,def_nj,1);
  WRITEI(&def->ni,def_ni,1);
  WRITEI(&def->dt,def_dt,1);
  WRITEI(&def->newt_cool_on,def_newt_cool_on,1);
  WRITEI(&def->newt_cool_adjust,def_newt_cool_adjust,1);
  WRITEI(&def->microphysics_on,def_microphysics_on,1);
  WRITEI(&def->nmt_physics_on,def_nmt_physics_on,1);
  WRITEI(&def->diffusion_direction,def_diffusion_direction,1);
  WRITEI(def->on,def_on,LAST_SPECIES+1);
  WRITEF(&def->du_vert,def_du_vert,1);
  WRITEI(&def->hgrid_mass_advection_scheme,   def_hgrid_mass_advection_scheme,   1);
  WRITEI(&def->hgrid_nonmass_advection_scheme,def_hgrid_nonmass_advection_scheme,1);
  WRITEI(&def->uv_timestep_scheme,def_uv_timestep_scheme,1);
  WRITEI(&def->turbulence_scheme,def_turbulence_scheme,1);
  WRITEI(&def->stability_factor,def_stability_factor,1);
  WRITEI(&def->spacing_type,def_spacing_type,1);
  WRITEI(&def->drag_v,def_drag_v,1);
  WRITEI(&def->drag_zonal_avg,def_drag_zonal_avg,1);
  WRITEI(&def->k_sponge,def_k_sponge,1);

  WRITEF(&def->globe_lonbot,def_globe_lonbot,1);
  WRITEF(&def->globe_lontop,def_globe_lontop,1);
  WRITEF(&def->globe_latbot,def_globe_latbot,1);
  WRITEF(&def->globe_lattop,def_globe_lattop,1);
  WRITEF(&def->lat_tp,def_lat_tp,1);
  WRITEF(&def->f_plane_half_width,def_f_plane_half_width,1);
  WRITEF(&def->ptop,def_ptop,1);
  WRITEF(&def->pbot,def_pbot,1);
  WRITEF(&def->p_sigma,def_p_sigma,1);
  WRITEF(&def->prandtl,def_prandtl,1);
  WRITEF(&def->tau_drag,def_tau_drag,1);
  WRITEF(&def->nudiv_nondim,def_nudiv_nondim,1);
  WRITEF(&def->u_scale,def_u_scale,1);
  WRITEF(&def->u_spinup_scale,def_u_spinup_scale,1);
  WRITEF(&def->topo_scale,def_topo_scale,1);
  WRITEF(def->nu_nondim,def_nu_nondim,MAX_NU_ORDER+1);
  WRITEF(&def->time_fp,def_time_fp,1);
  WRITEF(def->mole_fraction,def_mole_fraction,LAST_SPECIES+1);
  WRITEF(def->mole_fraction_over_solar,def_mole_fraction_over_solar,LAST_SPECIES+1);
  WRITEF(def->rh_max,def_rh_max,LAST_SPECIES+1);

  /*
   * nmt_physics
   */
  WRITEI(&def->nmt_land,def_nmt_land,1);
  WRITEI(&def->nmt_dorad,def_nmt_dorad,1);
  WRITEI(&def->nmt_frad,def_nmt_frad,1);
  WRITEF(&def->nmt_radcool,def_nmt_radcool,1);
  WRITEF(&def->nmt_tpause,def_nmt_tpause,1);
  WRITEF(&def->nmt_radbrk,def_nmt_radbrk,1);
  WRITEF(&def->nmt_cvc,def_nmt_cvc,1);
  WRITEF(&def->nmt_cvs,def_nmt_cvs,1);
  WRITEF(&def->nmt_cvp,def_nmt_cvp,1);
  WRITEF(&def->nmt_cve,def_nmt_cve,1);
  WRITEF(&def->nmt_theslop,def_nmt_theslop,1);
  WRITEF(&def->nmt_pstiff,def_nmt_pstiff,1);
  WRITEF(&def->nmt_pscale,def_nmt_pscale,1);
  WRITEF(&def->nmt_pbltop,def_nmt_pbltop,1);
  WRITEF(&def->nmt_cdrag,def_nmt_cdrag,1);
  WRITEF(&def->nmt_wscale,def_nmt_wscale,1);
  WRITEF(&def->nmt_sst,def_nmt_sst,1);
  WRITEF(&def->nmt_lfrac,def_nmt_lfrac,1);
  WRITEF(&def->nmt_eflux0,def_nmt_eflux0,1);
  WRITEF(&def->nmt_cld,def_nmt_cld,1);
  WRITEF(&def->nmt_cfract,def_nmt_cfract,1);

  nc_close(nc_id);

  return;
}

/*======================= end of write_defaults() ===========================*/

/* * * * * * * * * * * *  end of epic_initial.c  * * * * * * * * * * * * * * */







