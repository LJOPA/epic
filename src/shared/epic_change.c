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

/* * * * * * * * * * epic_change.c * * * * * * * * * * * * * * * * * 
 *                                                                 *
 *       Makes changes to parameters in epic.nc                    *
 *                                                                 *
 *       Use "-spots spots.dat" to add vortices.                   *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*
 * Function prototypes:
 */
void add_spots(planetspec *planet,
               char       *spots_file);

void add_noise(planetspec *planet,
               char       *noise_file);

void read_defaults(change_defaultspec  *def);

void write_defaults(change_defaultspec *def);


/*======================= main() =====================================*/

int main(int   argc,
         char *argv[])
{
/*
 *  NOTE: structures planet, grid, and var are declared globally in epic.h.
 */
  char   
    spots_file[FILE_STR], /*  added-spot locations and sizes         */
    noise_file[FILE_STR], /*  spots for velocity noise perturbation  */
    sflag[80],            /*  string to hold command-line flags      */
    infile[ FILE_STR],
    outfile[FILE_STR],
    buffer[16];
  int    
    time_index,
    nk,nj,ni,
    K,J,I,
    is,itmp,
    count,index,ii,
    spots      = FALSE,
    noise      = FALSE,
    stretch_ni = FALSE;
  EPIC_FLOAT  
    dx0,
    dt;
  double
    max_nu[MAX_NU_ORDER+1];
  change_defaultspec
    defaults;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="epic_change";

  declare_copyright();

  /* 
   * Interpret command-line arguments: 
   */
  /* Start with defaults: */
  sprintf(spots_file,"none");
  sprintf(noise_file,"none");
  if (argc > 1) {
    /* Read flags: */
    for (count = 1; count < argc; count++) {
      sscanf(argv[count],"%s",sflag);
      if (strcmp(sflag,"-spots") == 0) {
        sscanf(argv[++count],"%s",spots_file);
        spots = TRUE;
      }
      else if (strcmp(sflag,"-noise") == 0) {
        sscanf(argv[++count],"%s",noise_file);
        noise = TRUE;
      }
      else if (strcmp(sflag,"-stretch_ni") == 0) {
        sscanf(argv[++count],"%d",&stretch_ni);
        /*
         * Verify that stretch_ni is a power of 2.
         */
        if (frexp((double)stretch_ni,&itmp) != 0.5) {
          sprintf(Message,"-stretch_ni %d is not a power of 2",stretch_ni);
          epic_error(dbmsname,Message);
        }
      }
      else if (strcmp(sflag,"-help") == 0 ||
               strcmp(sflag,"-h")    == 0) {
        /* Print help, exit: */
        system("more "EPIC4_PATH"/help/epic_change.help");
        exit(1);
      }
      else {
        sprintf(Message,"Unrecognized change command-line flag: %s \n",sflag);
        epic_error(dbmsname,Message);
      }
    }
  }

  /* Allocate memory */
  if((planet=( planetspec *)malloc(sizeof(planetspec))) == 0) {
    sprintf(Message,"allocating space for planetspec \n");
    epic_error(dbmsname,Message);
  }

#if defined(EPIC_MPI)
  MPI_Init(&argc,&argv);
  para.comm = MPI_COMM_WORLD;
  MPI_Errhandler_set(para.comm,MPI_ERRORS_RETURN);
  MPI_Comm_rank(para.comm,&para.iamnode);
  MPI_Comm_size(para.comm,&para.nproc);
  para.ndim = NINT(log((EPIC_FLOAT)para.nproc)/log(2.));
#endif

  /*
   *  Read in default parameter settings:
   */
  read_defaults(&defaults);

  /* Time-plane index. */
  time_index = 0;

  /*
   * Determine model size from epic.nc, allocate memory for arrays,
   * and read in the rest of the data:
   */
  input_string("Input file [netCDF format]\n",defaults.infile,infile);
  /* NOTE: time_index is not used for SIZE_DATA */
  var_read(planet,infile,SIZE_DATA,time_index);

  set_var_props();
  make_arrays(planet);

  var_read(planet,infile,POST_SIZE_DATA,time_index);

  /* timeplane_bookkeeping() must come after reading in variables. */
  timeplane_bookkeeping();

  /* 
   * Set lon, lat, etc. 
   */
  set_lonlat();
  set_fmn(planet);
  set_gravity(planet);
  set_dsgth();
  set_sponge();

  /*
   * Set up thermodynamics subroutines.
   *
   * NOTE: The return value cpr is a low-temperature reference value, and
   *       should not be used otherwise.  Use return_cp() for a given
   *       thermodynamical state.
   */
  thermo_setup(planet,&planet->cpr);

  /* 
   * Store diagnostic variables. 
   */
  fprintf(stdout,"\nCalculating and storing diagnostic variables...");

  /*
   * Reconstitute the prognostic variable HDRY from input P3 and any Qs.
   */
  for (K = KLO; K <= KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        HDRY(K,J,I) = get_h(planet,2*K,J,I,DRY);
      }
    }
  }
  bc_lateral(var.hdry.value,THREEDIM);

  set_p2etc(planet,UPDATE_THETA);
  store_pgrad_vars(planet);
  store_diag(planet);

  fprintf(stdout,"done\n");

  nk = grid.nk;
  nj = grid.nj;
  ni = grid.ni;

  /*
   * Compute hyperviscosity coefficients.
   */
  set_max_nu(max_nu);

  /*
   * Print out a select listing of model parameters.
   */
  print_model_description(planet);

  /* 
   *  Print out vertical information:
   */
  print_vertical_column(planet,JLO,ILO,"vertical.dat");

  /*
   *  Change parameters as instructed.
   */
  grid.dt = input_int("\nInput timestep\n", grid.dt);

  /*
   * Inquire about Newtonian cooling and internal Rayleigh drag:
   */
  grid.newt_cool_on = input_int("Newtonian cooling on or off? [1 or 0]\n",grid.newt_cool_on);

  if (grid.newt_cool_on == TRUE) {
    grid.newt_cool_adjust = input_int("Adjust layer average of Newtonian cooling to zero? [1=yes, 0=no]\n",
                                       grid.newt_cool_adjust);
    grid.prandtl = input_float("Input relaxation ratio [tau_rad/tau_drag, 0. => no drag]\n",grid.prandtl);
  }
  else {
    if (grid.tau_drag == 0. || grid.tau_drag == 1.e+20) {
      grid.tau_drag = 
        input_float("Rayleigh drag timescale [s, 1.e+20 => infinity]\n",1.e+20);
    }
    else {
      grid.tau_drag = 
        input_float("Rayleigh drag timescale [s, 1.e+20 => infinity]\n",grid.tau_drag);
    }
    if (grid.tau_drag > 0. && grid.tau_drag < 1.e20) {
      grid.drag_v         = input_int("Apply Rayleigh drag to v? [1=yes, 0=no]\n",grid.drag_v);
      grid.drag_zonal_avg = input_int("Apply Rayleigh drag to zonal average velocity components? [1=yes, 0=no]\n",grid.drag_zonal_avg);
    }
  }

  /*
   * Inquire about microphysics (phase changes, latent heating, etc.).
   */
  count = 0;
  if (var.fpara.on) {
    count++;
  }
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      count++;
    }
  }

  /*
   * Set sponge:
   */
  grid.k_sponge = input_int("Input k_sponge (0 = no effect):\n",grid.k_sponge);

  /*
   * Recompute hyperviscosity coefficients (since dt may have changed).
   */
  set_max_nu(max_nu);

  /* grid.nu[2] is not used. */
  grid.nu[2] = 0.;
  sprintf(Message,"Divergence damping (fraction of max) \n");
  grid.nudiv_nondim = input_float(Message,grid.nudiv_nondim);
  for (ii = 4; ii <= MAX_NU_ORDER; ii+=2) {
    if (max_nu[ii] > 0.) {
      sprintf(Message,"Input nu[%d] (fraction of max) \n",ii);
      grid.nu_nondim[ii] = input_float(Message,grid.nu_nondim[ii]);
      grid.nu[ii]        = max_nu[ii]*grid.nu_nondim[ii];
    }
  }
  
  /*
   * Add vortices if requested:
   */
  if (spots) {
    add_spots(planet,spots_file);
  }

  /*
   * Add noise to (u,v) if requested:
   */
  if (noise) {
    add_noise(planet,noise_file);
  }

  /*
   * Prompt for which variables to write to extract.nc.
   */
  defaults.extract_str[0] = '\0';
  for (index = FIRST_INDEX; index <= LAST_INDEX; index++) {
    if (var.extract_on_list[index] == 1) {
      sprintf(buffer," %d",index);
      strcat(defaults.extract_str,buffer);
    }
  }
  prompt_extract_on(defaults.extract_str);

  /*
   * Write epic.nc file.
   */
  input_string("Name of output file\n",defaults.outfile,outfile);

  var_write(planet,outfile,ALL_DATA,time_index,stretch_ni);

  /* Write defaults file: */
  write_defaults(&defaults);

  free_arrays(planet);
  free(planet);
  free_var_props();

  return 0;
}

/*======================= end of main() =====================================*/

/*======================= add_spots() =======================================*/

/*
 * Add vortices via a perturbation streamfunction.
 *
 * NOTE: Not MPI ready. 
 */
#undef  PERT
#define PERT(k,j,i) pert[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  UG
#define UG(k,j,i) ug[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  VG
#define VG(k,j,i) vg[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  ZETAG
#define ZETAG(j,i) zetag[i+(j)*Iadim-Shift2d]

#undef  KING
#define KING(j,i) king[i+(j)*Iadim-Shift2d]

/*
 * Implemented styles for the vortex perturbation streamfunction:
 */
#define LEBEAU_DOWLING_1998            -1    /* Special case, see below. */
#define GAUSSIAN_ELLIPSOID              0
#define POLYNOMIAL_GAUSSIAN_ELLIPSOID   1
#define POLYNOMIAL_GAUSSIAN_ORDER       2.0

/*
 * Choose streamfunction style:
 */
#define SPOT_MONT  GAUSSIAN_ELLIPSOID

/*
 * Specify whether or not to apply the McIntyre and Roulstone (2002)
 * gradient-balance correction term (TRUE or FALSE).
 */
#define MCINTYRE_CORRECTION TRUE

void add_spots(planetspec *planet,
               char       *spots_file)
{
  register int
    K,kk,
    J,jj,
    I,i,ispot;
  int
    nspots=0;
  EPIC_FLOAT
    rr,xspot,yspot,zspot,fspot,
    exnerspot,alphaspot,tspot,
    uspot,vspot,
    sigma,theta,dtheta,
    t3kji,mont3,dmont,mu,fpara,pressure,
    dgz,dp,
    theta_ortho,theta_para,
    *ampspot,
    *lonspot,*latspot,*pspot,
    *aspot,*bspot,*cspot_up,*cspot_down,
    *pert,
    *ug,
    *vg,
    *zetag,
    *king,
     factor,
    lon_width,
    lon_half_width;
  char
    *char_pt,
    buffer[FILE_STR];
  FILE
    *spots;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="add_spots";

  if (strcmp(spots_file,"none") == 0) {
    /* Return if there is nothing to do: */
    return;
  }

  if ( SPOT_MONT != LEBEAU_DOWLING_1998 ) {
    /* Read vortex description file: */
    lon_width      = grid.globe_lontop - grid.globe_lonbot;
    lon_half_width = 0.5 * lon_width;
    nspots = number_spots_in_file( spots_file );
  }
  else if (SPOT_MONT == LEBEAU_DOWLING_1998) {
    nspots = 1;
  }
  else {
    sprintf(Message,"Unrecognized SPOT_MONT=%d",SPOT_MONT);
    epic_error(dbmsname,Message);
  }

  if (nspots == 1) {
    fprintf(stdout,"Generating spot...\n"); fflush(stdout);
  } else {
    fprintf(stdout,"Generating spots...\n"); fflush(stdout);
  }

  /* 
   * Allocate memory:
   */
  lonspot    = fvector(0,nspots-1, dbmsname);
  latspot    = fvector(0,nspots-1, dbmsname);
  pspot      = fvector(0,nspots-1, dbmsname);
  aspot      = fvector(0,nspots-1, dbmsname);
  bspot      = fvector(0,nspots-1, dbmsname);
  cspot_up   = fvector(0,nspots-1, dbmsname);
  cspot_down = fvector(0,nspots-1, dbmsname);
  ampspot    = fvector(0,nspots-1, dbmsname);

  pert       = fvector(0,Nelem3d-1,dbmsname);
  ug         = fvector(0,Nelem3d-1,dbmsname);
  vg         = fvector(0,Nelem3d-1,dbmsname);

  zetag      = fvector(0,Nelem2d-1,dbmsname);
  king       = fvector(0,Nelem2d-1,dbmsname);

  if ( SPOT_MONT != LEBEAU_DOWLING_1998 ) {
    /* 
     * Read in vortex information: 
     */
    read_spots_file( spots_file, 
                     ampspot,
                     lonspot,
                     latspot,
                     pspot,
                     aspot,
                     bspot,
                     cspot_up,
                     cspot_down,
                     ADJUST_SPOT_AMPLITUDE );
  }

  /* Clear PERT memory. */
  memset(pert,0,Nelem3d*sizeof(EPIC_FLOAT));

  /* 
   * Calculate streamfunction for vortices, and store
   * in PERT memory.
   */
  if ( SPOT_MONT != LEBEAU_DOWLING_1998 ) {
    for (K = KLO; K < KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          pressure = P3(K,J,I);
          for (ispot = 0; ispot < nspots; ispot++) {
            /* Account for periodicity in x-direction: */
            xspot  = (grid.lon[2*I+1]-lonspot[ispot]); 
            if (xspot > lon_half_width) { 
              xspot -= lon_width; 
            } else if (xspot < -lon_half_width) { 
              xspot += lon_width; 
            } 
            xspot /= aspot[ispot]; 

            yspot  = (grid.lat[2*J+1]-latspot[ispot])/bspot[ispot];

	    rr     = xspot*xspot+yspot*yspot;
            if (pressure <= pspot[ispot]){
              zspot = -log(pressure/pspot[ispot])/cspot_up[ispot];
	    }
            else{
              zspot =  log(pressure/pspot[ispot])/cspot_down[ispot];
            }
            rr += zspot*zspot; 
      
#if (SPOT_MONT == GAUSSIAN_ELLIPSOID)
            PERT(K,J,I) += ampspot[ispot]*exp(-rr);
#elif (SPOT_MONT == POLYNOMIAL_GAUSSIAN_ELLIPSOID)
            rr = sqrt(rr);
            if ( rr <= pow(2.0,1./POLYNOMIAL_GAUSSIAN_ORDER) ) {
              rr = pow(rr,POLYNOMIAL_GAUSSIAN_ORDER);
              PERT(K,J,I) += ampspot[ispot]*(1.+(1.-rr)*exp(-rr+2.))/(1.+exp(2.));
            }
#endif
          }
        }
      }
    }
  }
  else if (SPOT_MONT == LEBEAU_DOWLING_1998) {
    /*
     * NOTE: This is the specialized initial vortex that was 
     *       used in the LeBeau & Dowling (1998, Icarus) study of Neptune's GDS,
     *       specifically Run 90; it is not designed for general use.
     *       The nominal case should be run with nk = 10, nj = 128, ni = 256,
     *       latbot = -90., lattop = 0., lonbot = -90., lontop = 90.
     * 
     */
    int
      j0;
    EPIC_FLOAT
      AR     = 1.5,
      ampK[] = {0.,0.,0.,0.,0.,.2,.6,1.,1.,.5,0.},
      mm0,aa0,aa2,a0,f0,
      lin,lout,rel,
      rln,rlt,ratio;

    /*
     * Nominal: NINT((EPIC_FLOAT)(43*grid.nj)/64.);
     *
     *    j0 = NINIT((EPIC_FLOAT)(41*grid.nj)/64.); gives good results.
     */
    j0 = NINT((EPIC_FLOAT)(41*grid.nj)/64.);

    lonspot[0] = grid.lon[2*grid.ni/2+1];
    latspot[0] = grid.lat[2*j0+1];

    f0  = grid.f[2*j0+1];
    rln = planet->re/sqrt(1.+pow(planet->rp/planet->re*tan(latspot[0]*DEG),2.));
    rlt = rln/(cos(latspot[0]*DEG)*(pow(sin(latspot[0]*DEG),2.)+
                                    pow(planet->re/planet->rp*cos(latspot[0]*DEG),2.)));

    /*
     * NOTE: The original LeBeau and Dowling (1998) Run 90 source code writes lout (Rlimit) 
     *       as 50/grid.n[j0] instead of 50./grid.n[2*j0+1] and lin (Rinner) as
     *       14./grid.n[j0] instead of 14./grid.n[2*j0+1]. 
     *       Here we use the original source code's resultant values.  
     */      
    lout  = 1.51714e+7;
    a0    = lout/3.;
    lin   = 4.24798e+6;

    ratio = lin/a0;
    mm0   = 35.e+3;
    aa2   = -ratio*sech2(ratio)*tanh(ratio);
    aa0   = sech2(ratio)-aa2;

    for (K = 5; K <= 9; K++) {
      for (J = JLO; J <= JHI; J++) {
        yspot  = rlt*(grid.lat[2*J+1]-latspot[0])*DEG;
        factor = ampK[K]*mm0*grid.f[2*J+1]/f0;
        for (I = ILO; I <= IHI; I++) {
          xspot  = rln*(grid.lon[2*I+1]-lonspot[0])*DEG;
          rel    = sqrt(xspot*xspot+AR*AR*yspot*yspot);
          if (rel <= lin) {
            PERT(K,J,I) = aa0+aa2*rel*rel/(lin*lin);
          }
          else if (rel <= lout) {
            PERT(K,J,I) = sech2(rel/a0);
          }
          else {
            PERT(K,J,I) = 0.;
          }
          PERT(K,J,I) *= factor;
        }
      }
    }
  }
  else {
    sprintf(Message,"Unrecognized SPOT_MONT=%d",SPOT_MONT);
    epic_error(dbmsname,Message);
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(pert,THREEDIM);

  /*
   * Estimate the geostrophic wind (UG,VG) of the perturbed system.
   * Remove the EXNER Grad ( THETA ) component from the hybrid layers, in preparation for adding in the new one.
   *
   * NOTE: Need a similar term for the sigma region if the bottom pressure is not constant (otherwise, it is
   *       like a pressure-coordinate model, with no second term in the pressure-gradient force).
   */
  for (J = JLO+1; J <= JHI-1; J++) {
    jj = 2*J+1;
    if (fabs(grid.lat[jj]) > 10.) {
      for (I = ILO; I <= IHI; I++) {
        for (K = KLO; K < KHI; K++) {
          UG(K,J,I) = U(grid.it_uv,K,J,I)-(PERT(K,J+1,I-1)+PERT(K,J+1,I)-PERT(K,J-1,I-1)-PERT(K,J-1,I))*.25*grid.n[jj]/grid.f[jj];
        }
        for (K = KLO; K < grid.k_sigma-1; K++) {
          UG(K,J,I) -= (EXNER3(K,J,I)+EXNER3(K,J,I-1))
                      *(THETA(K,J+1,I-1)+THETA(K,J+1,I)-THETA(K,J-1,I-1)-THETA(K,J-1,I))*.125*grid.n[jj]/grid.f[jj];
        }
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  for (J = JFIRST; J <= JHI; J++) {
    jj = 2*J;
    if (fabs(grid.lat[jj]) > 10.) {
      for (I = ILO; I <= IHI; I++) {
        for (K = KLO; K < KHI; K++) {
          VG(K,J,I) = V(grid.it_uv,K,J,I)+(PERT(K,J,I+1)+PERT(K,J-1,I+1)-PERT(K,J,I-1)-PERT(K,J-1,I-1))*.25*grid.m[jj]/grid.f[jj];
        }
        for (K = KLO; K < grid.k_sigma-1; K++) {
          VG(K,J,I) += (EXNER3(K,J,I)+EXNER3(K,J-1,I))
                      *(THETA(K,J,I+1)+THETA(K,J-1,I+1)-THETA(K,J,I-1)-THETA(K,J-1,I-1))*.125*grid.m[jj]/grid.f[jj];
        }
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  /* 
   * Modify P3 and THETA.
   */
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /*
       * Hybrid region.
       */
      for (K = KLO; K <= grid.k_sigma-2; K++) {
        dtheta         = THETA(K-1,J,I)-THETA(K+1,J,I);
        dmont          = PERT( K-1,J,I)-PERT( K+1,J,I);
        /* 
         * NOTE: This is approximate, because theta changes with the introduction of the spot.
         */
        exnerspot      = dmont/dtheta;
        EXNER3(K,J,I) += exnerspot;
        state_from_exner(planet,2*K+1,J,I,EXNER3(K,J,I),&P3(K,J,I),&THETA(K,J,I));
      }
  
      /*
       * Sigma region.
       *
       * NOTE: Currently assuming bottom pressure is constant, such that the sigma region is actually
       *       a pressure-coordinate region.
       */
      for (K = grid.k_sigma-2; K < KHI; K++) {
        GZ3(K,J,I) += PERT(K,J,I);
      }
      for (K = grid.k_sigma-1; K < KHI; K++) {
        fpara        = (var.fpara.on) ? FPARA(K,J,I) : 0.25;
        mu           = avg_molar_mass(planet,2*K+1,J,I);
        RHO3(K,J,I)  = P3(K,J,I)*log(P3(K+1,J,I)/P3(K-1,J,I)) / (GZ3(K-1,J,I)-GZ3(K+1,J,I));
        T3(  K,J,I)  = alt_return_temp(planet,fpara,P3(K,J,I),mu,RHO3(K,J,I));
        THETA(K,J,I) = return_theta(planet,fpara,P3(K,J,I),T3(K,J,I),&theta_ortho,&theta_para);
      }
    }
  }

  /*
   * Update variables.
   */
  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        THETA2(K,J,I) = get_var(planet,THETA2_INDEX,NO_PHASE,grid.it_h,kk,J,I);
        H2(  K,J,I)   = get_h(planet,kk,J,I,TOTAL);
        HDRY(K,J,I)   = get_h(planet,kk,J,I,DRY);
      }
    }
    bc_lateral( var.theta2.value  +(K-Kshift)*Nelem2d, TWODIM );
    bc_lateral( var.h2.value      +(K-Kshift)*Nelem2d, TWODIM );
    bc_lateral( var.hdry.value    +(K-Kshift)*Nelem2d, TWODIM );
  }

  set_p2etc(planet,UPDATE_THETA);
  store_pgrad_vars(planet);
  store_diag(planet);

  /*
   * Add the new EXNER Grad( THETA ) component to hybrid layers.
   */
  for (J = JLO+1; J <= JHI-1; J++) {
    jj = 2*J+1;
    if (fabs(grid.lat[jj]) > 10.) {
      for (I = ILO; I <= IHI; I++) {
        for (K = KLO; K < grid.k_sigma-1; K++) {
          UG(K,J,I) += (EXNER3(K,J,I)+EXNER3(K,J,I-1))
                      *(THETA(K,J+1,I-1)+THETA(K,J+1,I)-THETA(K,J-1,I-1)-THETA(K,J-1,I))*.125*grid.n[jj]/grid.f[jj];
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(ug,THREEDIM);

  for (J = JFIRST; J <= JHI; J++) {
    jj = 2*J;
    if (fabs(grid.lat[jj]) > 10.) {
      for (I = ILO; I <= IHI; I++) {
        for (K = KLO; K < grid.k_sigma-1; K++) {
          VG(K,J,I) -= (EXNER3(K,J,I)+EXNER3(K,J-1,I))
                      *(THETA(K,J,I+1)+THETA(K,J-1,I+1)-THETA(K,J,I-1)-THETA(K,J-1,I-1))*.125*grid.m[jj]/grid.f[jj];
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(vg,THREEDIM);
  
#if MCINTYRE_CORRECTION
  /*
   * The gradient-balance correction term is based on the paper:
   *     McIntyre and Roulstone, 2002, Large-Scale Atmosphere-Ocean Dynamics. II Geometric Methods
   *         and Models. Cambridge Univ. Press, Ch. 8.
   */
  for (K = KLO; K < KHI; K++) {
    /*
     * Calculate geostrophic relative vorticity and kinetic energy per mass.
     */
    vorticity(planet,ON_SIGMATHETA,RELATIVE,ug+(K-Kshift)*Nelem2d,vg+(K-Kshift)*Nelem2d,NULL,zetag);
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        KING(J,I) = get_kin(planet,ug+(K-Kshift)*Nelem2d,vg+(K-Kshift)*Nelem2d,J,I);
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(king,TWODIM);

    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      if (fabs(grid.lat[jj]) > 10.) {
        for (I = ILO; I <= IHI; I++) {
          UG(K,J,I) += (-.5*(ZETAG(J,I)+ZETAG(J+1,I))*UG(K,J,I)
                        -.25*grid.n[jj]*(KING(J+1,I)+KING(J+1,I-1)-KING(J-1,I)-KING(J-1,I-1)))/grid.f[jj];
        }
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      jj = 2*J;
      if (fabs(grid.lat[jj]) > 10.) {
        for (I = ILO; I <= IHI; I++) {
          VG(K,J,I) += (-.5*(ZETAG(J,I)+ZETAG(J,I+1))*VG(K,J,I)
                        +.25*grid.m[jj]*(KING(J,I+1)+KING(J-1,I+1)-KING(J,I-1)-KING(J-1,I-1)))/grid.f[jj];
        }
      }
    }
  }
#endif

  for (K = KLO; K < KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      if (fabs(grid.lat[jj]) > 10.) {
        for (I = ILO; I <= IHI; I++) {
          U(grid.it_uv,K,J,I) = UG(K,J,I);
        }
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      jj = 2*J;
      if (fabs(grid.lat[jj]) > 10.) {
        for (I = ILO; I <= IHI; I++) {
          V(grid.it_uv,K,J,I) = VG(K,J,I);
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);
  bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);
   
  fprintf(stdout,"done.\n"); fflush(stdout);

  /* Free allocated memory: */
  free_fvector(zetag,     0,Nelem2d-1,dbmsname);
  free_fvector(king,      0,Nelem2d-1,dbmsname);

  free_fvector(ug,        0,Nelem3d-1,dbmsname);
  free_fvector(vg,        0,Nelem3d-1,dbmsname);
  free_fvector(pert,      0,Nelem3d-1,dbmsname);

  free_fvector(ampspot,   0,nspots-1, dbmsname);
  free_fvector(cspot_down,0,nspots-1, dbmsname);
  free_fvector(cspot_up,  0,nspots-1, dbmsname);
  free_fvector(bspot,     0,nspots-1, dbmsname);
  free_fvector(aspot,     0,nspots-1, dbmsname);
  free_fvector(pspot,     0,nspots-1, dbmsname);
  free_fvector(latspot,   0,nspots-1, dbmsname);
  free_fvector(lonspot,   0,nspots-1, dbmsname);

  return;
}

/*======================= end of add_spots() ================================*/

/*======================= add_noise() =======================================*/

/*
* Add noise to the wind field.
* Kunio Sayanagi, 11-07-07
* This function perturbs the wind velocity field.
* It takes the same input file as add_spots, and is based on add_spots
*
* NOTE: Not MPI ready.
*/
#undef  PERT_U
#define PERT_U(k,j,i) pert_u[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  PERT_V
#define PERT_V(k,j,i) pert_v[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

void add_noise(planetspec *planet,
               char       *noise_file)
{
 register int
   K,J,I,jj,ispot,i;
 int
   nspots=0;
 EPIC_FLOAT
   rr,xspot,yspot,zspot,fspot,
   pressure,
   *ampspot,
   *lonspot,*latspot,*pspot,
   *aspot,*bspot,*cspot_up,*cspot_down,
   *pert_u,
   *pert_v,
    factor,
   lon_width,
   lon_half_width;
 char
   *char_pt,
   buffer[FILE_STR];
 FILE
   *spots;
 /*
  * The following are part of DEBUG_MILESTONE(.) statements:
  */
 int
   idbms=0;
 static char
   dbmsname[]="add_spots";

 int
  kount1, kount2, kount3;

 if (strcmp(noise_file,"none") == 0) {
   /* Return if there is nothing to do: */
   return;
 }

 /* Read vortex description file: */
 lon_width      = grid.globe_lontop - grid.globe_lonbot;
 lon_half_width = 0.5 * lon_width;
 nspots = number_spots_in_file(noise_file);


 fprintf(stdout,"add_noise: Perturbing the wind field ...\n"); fflush(stdout);

 /*
  * Allocate memory:
  */
 lonspot    = fvector(0,nspots-1, dbmsname);
 latspot    = fvector(0,nspots-1, dbmsname);
 pspot      = fvector(0,nspots-1, dbmsname);
 aspot      = fvector(0,nspots-1, dbmsname);
 bspot      = fvector(0,nspots-1, dbmsname);
 cspot_up   = fvector(0,nspots-1, dbmsname);
 cspot_down = fvector(0,nspots-1, dbmsname);
 ampspot    = fvector(0,nspots-1, dbmsname);
 pert_u     = fvector(0,Nelem3d-1,dbmsname);
 pert_v     = fvector(0,Nelem3d-1,dbmsname);

 /*
  * Read in vortex information:
  */
 read_spots_file(noise_file,
                 ampspot,
                 lonspot,
                 latspot,
                 pspot,
                 aspot,
                 bspot,
                 cspot_up,
                 cspot_down,
                 DONT_ADJUST_SPOT_AMPLITUDE);

 /* Clear PERT memory. */
 memset(pert_u,0,Nelem3d*sizeof(EPIC_FLOAT));
 memset(pert_v,0,Nelem3d*sizeof(EPIC_FLOAT));

 /*
  * Calculate the perturbation fields, and store
  * in PERT memory.
  */
 for (K = KLO; K < KHI; K++) {
   for (J = JLO; J <= JHI; J++) {
     for (I = ILO; I <= IHI; I++) {
       pressure = P3(K,J,I);
       for (ispot = 0; ispot < nspots; ispot++) {
         /* Account for periodicity in x-direction: */
         xspot  = (grid.lon[2*I+1]-lonspot[ispot]);
         if (xspot > lon_half_width) {
           xspot -= lon_width;
         } else if (xspot < -lon_half_width) {
           xspot += lon_width;
         }
         xspot /= aspot[ispot];

         yspot  = (grid.lat[2*J+1]-latspot[ispot])/bspot[ispot];

         rr     = xspot*xspot+yspot*yspot;
         if (pressure <= pspot[ispot]){
           zspot = -log(pressure/pspot[ispot])/cspot_up[ispot];
         }
         else{
           zspot =  log(pressure/pspot[ispot])/cspot_down[ispot];
         }
         rr += zspot*zspot;

         /*
          * Kind of kludgie, but when ispot = even, perturb U and when odd, perturb V
          */
         if ((ispot%2) == 0) {
           PERT_U(K,J,I) += ampspot[ispot]*exp(-rr);
         }
         else {
           PERT_V(K,J,I) += ampspot[ispot]*exp(-rr);
         }
       }
     }
   }
 }

 /* Need to apply bc_lateral() here. */
 bc_lateral(pert_u,THREEDIM);
 bc_lateral(pert_v,THREEDIM);

 /*
  * Modify U and V.
  *
  * NOTE: Not MPI ready.
  */
 printf("add_noise: Modifing U and V ... \n");
 for (K = KLO; K < KHI; K++) {
   for (J = JLO; J <= JHI; J++) {
     jj = 2*J+1;
     for (I = ILO; I <= IHI; I++) {
       U(grid.it_uv,K,J,I) += PERT_U(K,J,I);
     }
   }
   /* Need to apply bc_lateral() here. */
   bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);

   for (J = JFIRST; J <= JHI; J++) {
     jj = 2*J;
     for (I = ILO; I <= IHI; I++) {
       V(grid.it_uv,K,J,I) += PERT_V(K,J,I);
     }
   }

   /* Need to apply bc_lateral() here. */
   bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);

 } /* (end loop over K) */


 fprintf(stdout,"done.\n"); fflush(stdout);

 /* Free allocated memory: */
 free_fvector(pert_u,    0,Nelem3d-1,dbmsname);
 free_fvector(pert_v,    0,Nelem3d-1,dbmsname);
 free_fvector(ampspot,   0,nspots-1, dbmsname);
 free_fvector(cspot_down,0,nspots-1, dbmsname);
 free_fvector(cspot_up,  0,nspots-1, dbmsname);
 free_fvector(bspot,     0,nspots-1, dbmsname);
 free_fvector(aspot,     0,nspots-1, dbmsname);
 free_fvector(pspot,     0,nspots-1, dbmsname);
 free_fvector(latspot,   0,nspots-1, dbmsname);
 free_fvector(lonspot,   0,nspots-1, dbmsname);

 return;
}

/*======================= end of add_noise() ================================*/

/*======================= read_defaults() ===================================*/

void read_defaults(change_defaultspec *def) 
{
  int
    nc_id,
    nc_err,
    index;
  char
    min_element[4];
  static char
    **gattname=NULL,
    **varname =NULL;
  static int
    ngatts    =0,
    num_progs =0;
  nc_type
    the_nc_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_defaults";

  nc_err = lookup_netcdf("change_defaults.nc",
                         &nc_id,&ngatts,&gattname,&num_progs,&varname);
  if (nc_err == NC_NOERR) {
    READI(&def->hgrid_mass_advection_scheme,def_hgrid_mass_advection_scheme,1);
    READI(&def->hgrid_nonmass_advection_scheme,def_hgrid_nonmass_advection_scheme,1);
    READI(&def->uv_timestep_scheme,def_uv_timestep_scheme,1);
    READC(def->infile,def_infile,N_STR);
    READC(def->outfile,def_outfile,N_STR);
    READC(def->extract_str,def_extract_str,N_STR);
  }
  else {
    /*
     * If the file is not readable, use standard defaults.
     */
    def->hgrid_mass_advection_scheme    = 0;
    def->hgrid_nonmass_advection_scheme = 0;
    def->uv_timestep_scheme             = 0;

    strcpy(def->infile,"epic.nc");
    strcpy(def->outfile,"epic.nc");
    strcpy(def->extract_str,"");
  }

  return;
}

/*======================= end of read_defaults() ============================*/

/*======================= write_defaults() ==================================*/

void write_defaults(change_defaultspec *def)
{
  int
    nc_id,
    nc_err;
  nc_type
    the_nc_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="write_defaults";

  nc_err = nc_create("change_defaults.nc",NC_CLOBBER,&nc_id);

  WRITEI(&def->hgrid_mass_advection_scheme,def_hgrid_mass_advection_scheme,1);
  WRITEI(&def->hgrid_nonmass_advection_scheme,def_hgrid_nonmass_advection_scheme,1);
  WRITEI(&def->uv_timestep_scheme,def_uv_timestep_scheme,1);
  WRITEC(def->infile,def_infile,N_STR);
  WRITEC(def->outfile,def_outfile,N_STR);
  WRITEC(def->extract_str,def_extract_str,N_STR);

  nc_close(nc_id);

  return;
}

/*======================= end of write_defaults() ===========================*/


/* * * * * * * * * * * * end of epic_change.c * * * * * * * * * * * * * * * * */

