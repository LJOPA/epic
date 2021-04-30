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

/* * * * * * * * * epic_funcs_init.c * * * * * * * * * * * * * * * * * * * * * 
 *                                                                           *
 *       Include here functions used to initialize or change the model       *
 *       that are not called while the model is running.  These functions    *
 *       do not need to be MPI ready.                                        *
 *                                                                           *
 *       This file contains the following functions:                         *
 *                                                                           *
 *           init_gz_surface()                                               *
 *           init_with_u()                                                   *
 *           init_with_hydrostatic()                                         *
 *           init_with_ref()                                                 *
 *           init_fpara_as_fpe()                                             *
 *           fpe_minus_fpe()                                                 *
 *           init_species()                                                  *
 *           init_species_via_deep_vmr()                                     *
 *           init_species_via_data()                                         *
 *           setup_mu_p(), mu_p()                                            *
 *           t_yp(),fp_yp()                                                  *
 *           p_gz()                                                          *
 *           rho_gz()                                                        *
 *           gz_p()                                                          *
 *           set_u_spinup()                                                  *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*====================== init_gz_surface() ===================================*/

/*
 * Set surface geopotential, gz.
 * The level gz = 0 corresponds to an appropriate reference radius.
 *
 * For terrestrial planets, use spherical harmonic gravity and surface data.
 *
 * For gas giants, use gradient balance with the zonal wind on a 
 * constant-pressure surface.
 *
 * The elevation can be suppressed with the artificial parameter
 * grid.topo_scale, which is applied when GZ_SURFACE(J,I) is used,
 * but not to the copy of GZ_SURFACE(J,I) stored in epic.nc.
 */

void init_gz_surface(planetspec       *planet,
                     init_defaultspec *def)
{
  int
    K,J,I,
    jlocal,ilocal,
    l,m,
    max_l_g = -1, /*g for gravity coeffs; r for radius coeffs */
    max_l_r = -1,
    max_l   = -1, 
    i;
  int
    num_files,nc_err,nc_id,
    nc_grid_ni,
    nc_grid_nj,
    file_match,
    node,num_nodes,
    jdim,idim,
    start[TWODIM],
    end[TWODIM];
  size_t
    index[1];
  EPIC_FLOAT
    nc_grid_globe_lonbot,
    nc_grid_globe_lontop,
    nc_grid_globe_latbot,
    nc_grid_globe_lattop;
  char
    nc_grid_geometry[GEOM_STR],
    nc_planet_name[N_STR],  /* NOTE: If [16], on Darwin next variable sometimes overwritten */
    gz_surface_nc[N_STR];
  EPIC_FLOAT
    lonr,latr,sin_lat,cos_lat,
    sum,r,gz,phi0,
    rerp,
    r0 = 0.,
    two_l,lplusm,lminusm;
  EPIC_FLOAT
   ***p_lm,
    **cos_mlon,
    **sin_mlon,
    **c_lm_g,  /*g for gravity coeffs; r for radius coeffs */
    **c_lm_r,
    **s_lm_g,
    **s_lm_r,
     *r0_r;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    *u1d,*p1d,*rho1d,*gz1d;
  char
    data_file_name[N_STR],
    header[N_STR];
  nc_type
    nc_float_type;
  FILE
    *infile_r,
    *infile_g;
  struct dirent
    **namelist;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_gz_surface";

#if EPIC_PRECISION == DOUBLE_PRECISION
   nc_float_type  = NC_DOUBLE;
#else
   nc_float_type  = NC_FLOAT;
#endif

  if (!initialized) {
    if (strcmp(planet->type,"gas-giant") == 0) {
      /*
       * Allocate memory.
       */
      u1d    = fvector(0,JADIM-1,dbmsname);
      p1d    = fvector(0,JADIM-1,dbmsname);
      rho1d  = fvector(0,JADIM-1,dbmsname);
      gz1d   = fvector(0,JADIM-1,dbmsname);
    }
    else if (strcmp(planet->type,"terrestrial") == 0) {
      ;
    }
    else {
      /*
       * Screen for unimplemented planet type.
       */
      sprintf(Message,"unrecognized planet->type=%s",planet->type);
      epic_error(dbmsname,Message);
    }

    initialized = TRUE;
  }

  /* * * * * * * * * * * * * * *
   *                           *
   *  Planet type: gas-giant   *
   *                           *
   * * * * * * * * * * * * * * */

  if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * Set gz_surface for gas giants using gradient balance with the bottom zonal wind.
     * Place the zero at J = grid.jtp.
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

    return;
  }

  /* * * * * * * * * * * * * * * *
   *                             *
   *  Planet type: terrestrial   *
   *                             *
   * * * * * * * * * * * * * * * */

  /*
   * Handle special cases first.
   */
  if (strcmp(planet->name,"held_suarez") == 0) {
    /*
     * Held-Suarez test case has no topography.
     * Set GZ_SURFACE(J,I) = 0. and return.
     */
    memset(var.gz_surface.value,0,sizeof(EPIC_FLOAT)*Nelem2d);
    return;
  }

  /* 
   * file_match flags whether an appropriate gz_surface.nc file 
   * already exists. 
   */
  file_match = FALSE;

  if (IAMNODE == NODE0) {
    /*
     * Search to see if appropriate gz_surface.nc data file
     * already exists.
     */

    /*
     * On sun4, rs6000, and sp2 platforms, we find that calling the function alphasort 
     * in the last argument of scandir() generates errors, even though it is
     * described in their man pages.  Since we don't need the directory
     * entries to be alphabetized, we fall back to using the NULL argument.
     */
    num_files = scandir(".",&namelist,is_gz_surface_file,NULL);

    if (num_files < 0) {
      /*
       * Problem encountered reading the files.
       */
      perror("init_gz_surface():scandir");
      exit(1);
    }
    else if (num_files > 0) {
      /*
       * There are a total of num_files gz_surface netCDF data files.
       * Search them to see if one matches what is needed.
       */
      for (l = 0; l < num_files; l++) {
        /*
         * Open file.
         */
        sprintf(gz_surface_nc,"./%s",namelist[l]->d_name);
        nc_err = nc_open(gz_surface_nc,NC_NOWRITE,&nc_id);
        if (nc_err != NC_NOERR) {
          sprintf(Message,"%s, %s",nc_strerror(nc_err),gz_surface_nc);
          epic_error(dbmsname,Message);
        }

        /*
         * Check to see if this file matches what is needed.
         */

        /* name */
        nc_err = nc_get_att_text(nc_id,NC_GLOBAL,"planet_name",nc_planet_name);
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_gz_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (strcmp(nc_planet_name,planet->name) != 0) {
          nc_close(nc_id); 
          continue;
        }

        /* geometry */
        nc_err = nc_get_att_text(nc_id,NC_GLOBAL,"grid_geometry",nc_grid_geometry);
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_gz_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (strcmp(nc_grid_geometry,grid.geometry) != 0) {
          nc_close(nc_id); 
          continue;
        }

        /* ni */
        nc_err = nc_get_att_int(nc_id,NC_GLOBAL,"grid_ni",&nc_grid_ni);
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_gz_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_ni != grid.ni) {
          nc_close(nc_id); 
          continue;
        }

        /* nj */
        nc_err = nc_get_att_int(nc_id,NC_GLOBAL,"grid_nj",&nc_grid_nj);
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_gz_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_nj != grid.nj) {
          nc_close(nc_id); 
          continue;
        }

        /* latbot */
#if EPIC_PRECISION == DOUBLE_PRECISION
        nc_err = nc_get_att_double(nc_id,NC_GLOBAL,"grid_globe_latbot",&nc_grid_globe_latbot);
#else
        nc_err = nc_get_att_float(nc_id,NC_GLOBAL,"grid_globe_latbot",&nc_grid_globe_latbot);
#endif
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_gz_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_globe_latbot != grid.globe_latbot) {
          nc_close(nc_id); 
          continue;
        }

        /* lattop */
#if EPIC_PRECISION == DOUBLE_PRECISION
        nc_err = nc_get_att_double(nc_id,NC_GLOBAL,"grid_globe_lattop",&nc_grid_globe_lattop);
#else
        nc_err = nc_get_att_float(nc_id,NC_GLOBAL,"grid_globe_lattop",&nc_grid_globe_lattop);
#endif
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_gz_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_globe_lattop != grid.globe_lattop) {
          nc_close(nc_id); 
          continue;
        }

        /* lonbot */
#if EPIC_PRECISION == DOUBLE_PRECISION
        nc_err = nc_get_att_double(nc_id,NC_GLOBAL,"grid_globe_lonbot",&nc_grid_globe_lonbot);
#else
        nc_err = nc_get_att_float(nc_id,NC_GLOBAL,"grid_globe_lonbot",&nc_grid_globe_lonbot);
#endif
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_gz_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_globe_lonbot != grid.globe_lonbot) {
          nc_close(nc_id); 
          continue;
        }

        /* lontop */
#if EPIC_PRECISION == DOUBLE_PRECISION
        nc_err = nc_get_att_double(nc_id,NC_GLOBAL,"grid_globe_lontop",&nc_grid_globe_lontop);
#else
        nc_err = nc_get_att_float(nc_id,NC_GLOBAL,"grid_globe_lontop",&nc_grid_globe_lontop);
#endif
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_gz_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_globe_lontop != grid.globe_lontop) {
          nc_close(nc_id); 
          continue;
        }

        /* 
         * Everything matches. 
         */
        file_match = TRUE;
        nc_close(nc_id);
        break;
      }
    }
  }
  
  if (file_match == TRUE) {
    /*
     * Input GZ_SURFACE(J,I) from file.
     */
    if (IAMNODE == NODE0) {
      fprintf(stderr,"\nReading GZ_SURFACE(J,I) from %s \n\n",
                     gz_surface_nc);
      /* Open file. */
      nc_err = nc_open(gz_surface_nc,NC_NOWRITE,&nc_id);
      if (nc_err != NC_NOERR) {
        fprintf(stderr,"Cannot find input file %s \n",gz_surface_nc);
        exit(1);
      }
    }
    num_nodes = setup_read_array();
    /*
     * Load start and end vectors.
     */
    start[0]  = 1;
    end[  0]  = grid.ni;

    for (node = 0; node < num_nodes; node++) {
      get_jlohi(node,num_nodes,start+1,end+1);
      read_array(node,TWODIM,start,end,var.gz_surface.info[0].name,
                 var.gz_surface.info[0].index,var.gz_surface.value,EPIC_FLOAT_ARRAY,nc_id);
    }
    if (IAMNODE == NODE0) {
      /* Close file. */
      nc_close(nc_id);
    }
    /* Need bc_lateral() here. */
    bc_lateral(var.gz_surface.value,TWODIM);

    return;
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   *                                                                       *
   * The required GZ_SURFACE does not exist, so we now construct it.       *
   *                                                                       *
   * We prefer to construct the topography from spherical harmonic data,   *
   * rather than gridded data, because the former provides a good          *
   * way to interpolate smoothly onto our user-defined model grid.         *
   *                                                                       *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  /*
   * Spherical harmonic data are truncated to maximum degree max_l_g for gravity data
   * and max_l_r for radius data.
   * Read in max_l_g and/or max_l_r.
   */
  if (IAMNODE == NODE0) {
    if (strncmp(planet->name,"venus",5) == 0) {
      /*
       * Venus topography data:
       *   Rappaport, N.J., A.S. Konopliv, A.B. Kucinskas, P.G. Ford, 1999, 
       *     An improved 360 degree and order model of Venus topography, Icarus 139, 19-31.
       */
      strcpy(data_file_name,EPIC4_PATH"/data/venus/TOPOCOEF-gtdr3.2");
      infile_r = fopen(data_file_name,"r");
      if (!infile_r) {
        sprintf(Message,"could not open %s",data_file_name);
        epic_error(dbmsname,Message);
      }
      fprintf(stderr,"Reading spherical harmonic (l,m) data from %s \n",data_file_name);
      fscanf(infile_r,"%d",&max_l_r);
    }
    else if (strncmp(planet->name,"earth",5) == 0) {
      /*
       * Earth topography data:
       *   Wieczorek MA, 2007, Gravity and topography of the terrestrial planets, 
       *     Treatise on Geophysics 10, 165-205.
       */
      strcpy(data_file_name,EPIC4_PATH"/data/earth/srtmp720.msl");
      infile_r = fopen(data_file_name,"r");
      if (!infile_r) {
        sprintf(Message,"could not open %s",data_file_name);
        epic_error(dbmsname,Message);
      }
      fprintf(stderr,"Reading spherical harmonic (l,m) data from %s \n",data_file_name);
      /* Skip header */
      for (i = 0; i < 7; i++) {
        fgets(header,N_STR,infile_r);
      }
      /* Input max_l_r */
      fscanf(infile_r,"%d",&max_l_r);
      /* Skip column heading */
      for (i = 0; i < 2; i++) {
        fgets(header,N_STR,infile_r);
      }
    }
    else if (strncmp(planet->name,"mars",4) == 0) {
      /*
       * Radius data from the Mars Orbiter Laser Altimeter (MOLA) on the Mars Global Surveyor (MGS) spacecraft.
       */  
      strcpy(data_file_name,EPIC4_PATH"/data/mars/gtm090aa_sha.txt");
      infile_r = fopen(data_file_name,"r");
      fprintf(stderr,"Reading spherical harmonic (l,m) radius data from %s \n",data_file_name);
      fscanf(infile_r," %*f, %*f, %*f, %d, %*d, %*d, %*f, %*f",&max_l_r);

      /*
       * Gravity data from:
       *   Zuber, M., Mars Reconnaissance Orbiter Derived Gravity Data, NASA Planetary 
       *     Data System, MRO-M-RSS-5-SDP-V1.0, 2008.
       */  
      strcpy(data_file_name,EPIC4_PATH"/data/mars/jgmro_095a_sha.tab");
      infile_g = fopen(data_file_name,"r");
      fprintf(stderr,"Reading spherical harmonic (l,m) gravity data from %s \n",data_file_name);

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(infile_g," %lf, %*lf, %*lf, %d, %*d, %*d, %*lf, %*lf",&r0,&max_l_g);
#else
      fscanf(infile_g," %f, %*f, %*f, %d, %*d, %*d, %*f, %*f",&r0,&max_l_g);
#endif

      /* Convert reference radius from km to m. */
      r0 *= 1.e+3;  
    }
    else {
      sprintf(Message,"failed to open spherical-harmonic file(s) for %s",planet->name);
      epic_error(dbmsname,Message);
    }
  }
  
  /*
   * Allocate memory for (l,m) spherical-harmonic coefficients.
   */
  max_l = IMAX(max_l_g,max_l_r);

  p_lm = (EPIC_FLOAT ***)calloc(JHI-JLO+1,sizeof(EPIC_FLOAT **));
  if (!p_lm) {
    sprintf(Message,"calloc error allocating p_lm");
    epic_error(dbmsname,Message);
  }
  for (J = JLO; J <= JHI; J++) {
    p_lm[J-JLO] = (EPIC_FLOAT **)calloc(max_l+1,sizeof(EPIC_FLOAT *));
    if (!p_lm[J-JLO]) {
      sprintf(Message,"calloc error allocating p_lm[%d-%d]",J,JLO);
      epic_error(dbmsname,Message);
    }
    for (l = 0; l <= max_l; l++) {
      p_lm[J-JLO][l] = fvector(0,l,dbmsname);
    }
  }
  cos_mlon = (EPIC_FLOAT **)calloc(IHI-ILO+1,sizeof(EPIC_FLOAT *));
  if (!cos_mlon) {
    sprintf(Message,"calloc error allocating cos_mlon");
    epic_error(dbmsname,Message);
  }
  sin_mlon = (EPIC_FLOAT **)calloc(IHI-ILO+1,sizeof(EPIC_FLOAT *));
  if (!sin_mlon) {
    sprintf(Message,"calloc error allocating sin_mlon");
    epic_error(dbmsname,Message);
  }
  for (I = ILO; I <= IHI; I++) {
    ilocal = I-ILO;
    cos_mlon[ilocal] = fvector(0,max_l,dbmsname);
    sin_mlon[ilocal] = fvector(0,max_l,dbmsname);
  }

  if (max_l_r >= 0) {
    c_lm_r = (EPIC_FLOAT **)calloc(max_l_r+1,sizeof(EPIC_FLOAT *));
    if (!c_lm_r) {
      sprintf(Message,"calloc error allocating c_lm_r");
      epic_error(dbmsname,Message);
    }

    s_lm_r = (EPIC_FLOAT **)calloc(max_l_r+1,sizeof(EPIC_FLOAT *));
    if (!s_lm_r) {
      sprintf(Message,"calloc error allocating s_lm_r");
      epic_error(dbmsname,Message);
    }

    for (l = 0; l <= max_l_r; l++) {
      c_lm_r[l] = fvector(0,l,dbmsname);
      s_lm_r[l] = fvector(0,l,dbmsname);
    }
  }

  if (max_l_g >= 0) {
    c_lm_g = (EPIC_FLOAT **)calloc(max_l_g+1,sizeof(EPIC_FLOAT *));
    if (!c_lm_g) {
      sprintf(Message,"calloc error allocating c_lm_g");
      epic_error(dbmsname,Message);
    }

    s_lm_g = (EPIC_FLOAT **)calloc(max_l_g+1,sizeof(EPIC_FLOAT *));
    if (!s_lm_g) {
      sprintf(Message,"calloc error allocating s_lm_g");
      epic_error(dbmsname,Message);
    }

    for (l = 0; l <= max_l_g; l++) {
      c_lm_g[l] = fvector(0,l,dbmsname);
      s_lm_g[l] = fvector(0,l,dbmsname);
    }
    r0_r = fvector(0,max_l_g,dbmsname);
  }

  /*
   * Input spherical-harmonic coefficients.
   */
  if (IAMNODE == NODE0) {
    if (strncmp(planet->name,"venus",5) == 0) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(infile_r," %lf",&r0);
#else
      fscanf(infile_r," %f",&r0);
#endif

      /* Convert reference radius from km to m. */
      r0 *= 1.e+3;
      for (l = 0; l <= max_l_r; l++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile_r," %lf",c_lm_r[l]);
#else
          fscanf(infile_r," %f",c_lm_r[l]);
#endif

        for (m = 1; m <= l; m++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile_r," %lf",c_lm_r[l]+m);
          fscanf(infile_r," %lf",s_lm_r[l]+m);
#else
          fscanf(infile_r," %f",c_lm_r[l]+m);
          fscanf(infile_r," %f",s_lm_r[l]+m);
#endif

        }
      }
      fclose(infile_r);
    }
    else if (strncmp(planet->name,"earth",5) == 0) {
      /* Data gives elevation rather than radius; signal this with r0 = 0. */
      r0 = 0.;
      for (l = 0; l <= max_l_r; l++) {
        for (m = 0; m <= l; m++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile_r,"%*d %*d %lf %lf",c_lm_r[l]+m,s_lm_r[l]+m);
#else
          fscanf(infile_r,"%*d %*d %f %f",c_lm_r[l]+m,s_lm_r[l]+m);
#endif

        }
      }
      fclose(infile_r);
    }
    else if (strncmp(planet->name,"mars",4) == 0) {
      for (l = 0; l <= max_l_r; l++) {
        for (m = 0; m <= l; m++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile_r," %*d, %*d, %lf, %lf, %*lf, %*lf",c_lm_r[l]+m,s_lm_r[l]+m);
#else
          fscanf(infile_r," %*d, %*d, %f, %f, %*f, %*f",c_lm_r[l]+m,s_lm_r[l]+m);
#endif

        }
      }
      fclose(infile_r);


      l = 0, m = 0;
      c_lm_g[l][m] = 1.;
      s_lm_g[l][m] = 0.;
      for (l = 1; l <= max_l_g; l++) {
        for (m = 0; m <= l; m++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile_g," %*d, %*d, %lf, %lf, %*lf, %*lf",c_lm_g[l]+m,s_lm_g[l]+m);
#else
          fscanf(infile_g," %*d, %*d, %f, %f, %*f, %*f",c_lm_g[l]+m,s_lm_g[l]+m);
#endif

        }
      }
      fclose(infile_g);
    }
    else {
      sprintf(Message,"failed to read spherical-harmonic topography-coefficient file for %s",planet->name);
      epic_error(dbmsname,Message);
    }
  }

  /*
   * Calculate GZ_SURFACE(J,I).
   */
  if (IAMNODE == NODE0) {
    fprintf(stderr,"Calculating topography: \n");
  }

  /*
   * Precalculate cos, sin factors.
   */
  if (IAMNODE == NODE0) {
    fprintf(stderr,"  Calculating and storing cos and sin factors.\n");
  }
  for (I = ILO; I <= IHI; I++) {
    ilocal = I-ILO;
    lonr   = grid.lon[2*I+1]*DEG;
    for (m = 0; m <= max_l; m++) {
      cos_mlon[ilocal][m] = cos((EPIC_FLOAT)m*lonr);
      sin_mlon[ilocal][m] = sin((EPIC_FLOAT)m*lonr);
    }
  }

  /*
   * Precalculate Legendre factors.
   *
   * We use a recursive formula by 
   *   Colombo C, 1981, Numerical methods for harmonic analysis on the sphere. Rep 310, 
   *     Dept. of Geodetic Science and Surveying, The Ohio State University, Columbus
   * which is cited and reproduced as eqn.(11) in Holmes and Featherstone (2002), 
   * and also described on the MITGCM website: http://mitgcm.org/~mlosch/geoidcookbook/node11.html.
   *
   * NOTE: We no longer use normed_legendre() from epic_funcs_util.c, which is based on the Numerical Recipes in C
   *       algorithm plgndr(), because it is slower and less accurate.
   */
  if (IAMNODE == NODE0) {
    fprintf(stderr,"  Calculating and storing normalized Legendre factors:   0%%");
  }

  rerp = planet->re/planet->rp;

  for (J = JLO; J <= JHI; J++) {
    jlocal  = J-JLO;
    latr    = lat_graphic_to_centric(grid.lat[2*J+1],rerp)*DEG;
    sin_lat = sin(latr);
    cos_lat = cos(latr);
    /*
     * Start with sectorial (l == m) values. These serve as seeds for the non-sectorial (l != m) recursion.
     * We divide out the cos_lat factor to avoid underflow for high degree and order, following 
     *   Holmes SA, Featherstone WE, 2002,  A unified approach to the Clenshaw summation and the recursive
     *     computation of very high degree and order normalised associated Legendre functions, 
     *     J. Geodesy 76, 279-299   
     */
    p_lm[jlocal][0][0] = 1.;
    p_lm[jlocal][1][1] = sqrt(3.);
    for (m = 2; m <= max_l; m++) {
      /*
       * The cos_lat factor is deliberately left off, and brought in below, following Holmes and Featherstone (2002).
       */
      p_lm[jlocal][m][m] = sqrt(1.+1./(2.*(double)m))*p_lm[jlocal][m-1][m-1];
    }
    /*
     * Calculate non-sectorial values.
     */
    for (l = 1; l <= max_l; l++) {
      two_l = 2.*(double)l;
      m     = l-1;
      p_lm[jlocal][l][m] = sqrt(two_l+1.)*p_lm[jlocal][l-1][m]*sin_lat;
      for (m = 0; m < l-1; m++) {
        lplusm  = (double)(l+m);
        lminusm = (double)(l-m);
        p_lm[jlocal][l][m] = p_lm[jlocal][l-1][m]*sqrt((two_l-1.)*(two_l+1.)/(lminusm*lplusm))*sin_lat
                            -p_lm[jlocal][l-2][m]*sqrt((two_l+1.)*(lplusm-1.)*(lminusm-1.)/(lminusm*lplusm*(two_l-3.)));
      }
    }
    if (IAMNODE == NODE0) {
      fprintf(stderr,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)(jlocal+1)/(JHI-JLO+1)));
    }
  }
  if (IAMNODE == NODE0) {
    fprintf(stderr,"\n");
  }

  if (IAMNODE == NODE0) {
    fprintf(stderr,"  Summing to get GZ_SURFACE(J,I):   0%%");
  }
  if (strncmp(planet->name,"venus",5) == 0) {
    /*
     * Coefficients yield radius divided by average radius, which we convert to elevation, z, and then to gz.
     */
    for (J = JLO; J <= JHI; J++) {
      jlocal = J-JLO;
      latr    = lat_graphic_to_centric(grid.lat[2*J+1],rerp)*DEG;
      cos_lat = cos(latr);
      for (I = ILO; I <= IHI; I++) {
        ilocal = I-ILO;
        /*
         * The cos_lat factor is folded back in here, following Holmes and Featherstone (2002).
         */
        sum = 0.;
        for (m = max_l_r; m >= 0; m--) {
          sum *= cos_lat;
          for (l = m; l <= max_l_r; l++) {
            sum += p_lm[jlocal][l][m]*(c_lm_r[l][m]*cos_mlon[ilocal][m]
                                      +s_lm_r[l][m]*sin_mlon[ilocal][m]);
          }
        }
        GZ_SURFACE(J,I) = grid.g[2*J+1]*r0*(sum-1.);
      }
      if (IAMNODE == NODE0) {
        fprintf(stderr,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)(jlocal+1)/(JHI-JLO+1)));
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(var.gz_surface.value,TWODIM);
    if (IAMNODE == NODE0) {
      fprintf(stderr,"\n\n");
    }
  }
  else if (strncmp(planet->name,"earth",5) == 0) {
    /*
     * Coefficients yield elevation, z, which we convert to gz.
     */
    for (J = JLO; J <= JHI; J++) {
      jlocal = J-JLO;
      latr    = lat_graphic_to_centric(grid.lat[2*J+1],rerp)*DEG;
      cos_lat = cos(latr);
      for (I = ILO; I <= IHI; I++) {
        ilocal = I-ILO;
        /*
         * The cos_lat factor is folded back in here, following Holmes and Featherstone (2002).
         */
        sum = 0.;
        for (m = max_l_r; m >= 0; m--) {
          sum *= cos_lat;
          for (l = m; l <= max_l_r; l++) {
            sum += p_lm[jlocal][l][m]*(c_lm_r[l][m]*cos_mlon[ilocal][m]
                                      +s_lm_r[l][m]*sin_mlon[ilocal][m]);
          }
        }
        /*
         * NOTE: Earth topography data set currently includes ocean-floor topography,
         *       which will need to be set to sea level for atmospheric applications.
         */
        GZ_SURFACE(J,I) = grid.g[2*J+1]*sum;
      }
      if (IAMNODE == NODE0) {
        fprintf(stderr,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)(jlocal+1)/(JHI-JLO+1)));
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(var.gz_surface.value,TWODIM);
    if (IAMNODE == NODE0) {
      fprintf(stderr,"\n\n");
    }
  }
  else if (strncmp(planet->name,"mars",4) == 0) {
    /*
     * For each horizontal position on the grid, we first sum to get the local surface radius,
     * then use this in the sum to get the local gravity.  Finally, we subtract the reference value to get the
     * surface gz in m^2/s^2 (or equivalently, J/kg).
     *
     * First, calculate the reference geopotential at the equator, phi0.
     */
    phi0 = 0.;
    for (l = 0; l <= max_l_g; l++) {
      phi0 += normed_legendre(l,0,0.)*c_lm_g[l][0];
    }
    phi0 = phi0*planet->GM/r0+.5*SQR(planet->omega_sidereal*r0);

    for (J = JLO; J <= JHI; J++) {
      jlocal = J-JLO;
      latr    = lat_graphic_to_centric(grid.lat[2*J+1],rerp)*DEG;
      cos_lat = cos(latr);
      for (I = ILO; I <= IHI; I++) {
        ilocal = I-ILO;
        /*
         * The cos_lat factor is folded back in here, following Holmes and Featherstone (2002).
         */
        r = 0.;
        for (m = max_l_r; m >= 0; m--) {
          r *= cos_lat;
          for (l = m; l <= max_l_r; l++) {
            r += p_lm[jlocal][l][m]*(c_lm_r[l][m]*cos_mlon[ilocal][m]
                                    +s_lm_r[l][m]*sin_mlon[ilocal][m]);
          }
        }

        /*
         * Calculate (r0/r)^l as an array.
         */
        r0_r[0] = 1.;
        for (l = 1; l <= max_l_g; l++) {
          r0_r[l] *= (r0/r);
        }

        sum = 0.;
        for (m = max_l_g; m >= 0; m--) {
          sum *= cos_lat;
          for (l = m; l <= max_l_g; l++) {
            sum += r0_r[l]*p_lm[jlocal][l][m]*(c_lm_g[l][m]*cos_mlon[ilocal][m]
                                              +s_lm_g[l][m]*sin_mlon[ilocal][m]);
          }
        }

        /*
         * Need to add in the rotational potential, .5*SQR(planet->omega_sidereal*r*cos(lat)).
         */
        GZ_SURFACE(J,I) = -(planet->GM*sum/r+.5*SQR(planet->omega_sidereal*r*cos_lat)-phi0);
      }
      if (IAMNODE == NODE0) {
        fprintf(stderr,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)(jlocal+1)/(JHI-JLO+1)));
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(var.gz_surface.value,TWODIM);
    if (IAMNODE == NODE0) {
      fprintf(stderr,"\n\n");
    }
  }

  /*
   * Free allocated memory.
   */
  if (max_l_r >= 0) {
    for (l = max_l_r; l >= 0; l--) {
      free_fvector(s_lm_r[l],0,l,dbmsname);
      free_fvector(c_lm_r[l],0,l,dbmsname);
    }
    free(s_lm_r);
    free(c_lm_r);
  }

  if (max_l_g >= 0) {
    for (l = max_l_g; l >= 0; l--) {
      free_fvector(s_lm_g[l],0,l,dbmsname);
      free_fvector(c_lm_g[l],0,l,dbmsname);
    }
    free(s_lm_g);
    free(c_lm_g);
    free_fvector(r0_r,0,max_l_g,dbmsname);
  }

  for (I = ILO; I <= IHI; I++) {
    ilocal = I-ILO;
    free_fvector(sin_mlon[ilocal],0,max_l,dbmsname);
    free_fvector(cos_mlon[ilocal],0,max_l,dbmsname);
  }
  free(sin_mlon);
  free(cos_mlon);

  for (J = JLO; J <= JHI; J++) {
    jlocal = J-JLO;
    for (l = 0; l <= max_l; l++) {
      free_fvector(p_lm[jlocal][l],0,max_l,dbmsname);
    }
    free(p_lm[jlocal]);
  }
  free(p_lm);

  /*
   * Write GZ_SURFACE(J,I) data to a gz_surface.nc file.
   */

  if (IAMNODE == NODE0) {
    /*
     * Enter define mode:
     */
    sprintf(gz_surface_nc,"./%s_%02d_gz_surface.nc",
            planet->name,num_files);
    nc_err = nc_create(gz_surface_nc,NC_CLOBBER,&nc_id);
    if (nc_err != NC_NOERR) {
      fprintf(stderr,"Cannot write %s \n",gz_surface_nc);
      exit(1);
    }
    WRITEC(planet->name,planet_name,32);
    WRITEI(&grid.nj,grid_nj,1);
    WRITEI(&grid.ni,grid_ni,1);
    WRITEC(grid.geometry,grid_geometry,GEOM_STR);
    WRITEF(&grid.globe_lonbot,grid_globe_lonbot,1);
    WRITEF(&grid.globe_lontop,grid_globe_lontop,1);
    WRITEF(&grid.globe_latbot,grid_globe_latbot,1);
    WRITEF(&grid.globe_lattop,grid_globe_lattop,1);

    /*
     * GZ_SURFACE(J,I) resides on the p-grid.
     */
    jdim = grid.nj-grid.jlo+1;
    idim = grid.ni;

    /* 
     * lon (I direction): 
     */
    nc_def_dim(nc_id,"lon_p",idim,&var.gz_surface.info[0].dimid[NETCDF_I_INDEX]);
    nc_def_var(nc_id,"lon_p",nc_float_type,ONEDIM,
               &var.gz_surface.info[0].dimid[NETCDF_I_INDEX],
               &var.gz_surface.info[0].coorid[NETCDF_I_INDEX]);
    nc_put_att_text(nc_id,var.gz_surface.info[0].coorid[NETCDF_I_INDEX],"units",
                    strlen("degrees_east")+1,"degrees_east");
    /* 
     * lat (J direction):
     */
    nc_def_dim(nc_id,"lat_p",jdim,&var.gz_surface.info[0].dimid[NETCDF_J_INDEX]);
    nc_def_var(nc_id,"lat_p",nc_float_type,ONEDIM,
               &var.gz_surface.info[0].dimid[NETCDF_J_INDEX],
               &var.gz_surface.info[0].coorid[NETCDF_J_INDEX]);
    nc_put_att_text(nc_id,var.gz_surface.info[0].coorid[NETCDF_J_INDEX],"units",
                    strlen("degrees_north")+1,"degrees_north");
    nc_put_att_text(nc_id,var.gz_surface.info[0].coorid[NETCDF_J_INDEX],"mapping",
                    strlen("planetographic")+1,"planetographic");

    /* gz_surface */
    nc_def_var(nc_id,"gz_surface",nc_float_type,TWODIM,
               &var.gz_surface.info[0].dimid[NETCDF_J_INDEX],
               &var.gz_surface.info[0].id);
    nc_put_att_text(nc_id,var.gz_surface.info[0].id,"units",
                    strlen("J/kg")+1,"J/kg");

    /*
     * Leave define mode:
     */
    nc_enddef(nc_id);

    /*
     * Assign values to coordinates.
     * Longitude:
     */
    for (I = 1; I <= grid.ni; I++) {
      index[0] = I-1;
#if EPIC_PRECISION == DOUBLE_PRECISION
      nc_put_var1_double(nc_id,var.gz_surface.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I+1]));
#else
      nc_put_var1_float(nc_id,var.gz_surface.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I+1]));
#endif
    }
    /*
     * Latitude:
     */
    for (J = grid.jlo; J <= grid.nj; J++) {
      index[0] = J-grid.jlo;
#if EPIC_PRECISION == DOUBLE_PRECISION
      nc_put_var1_double(nc_id,var.gz_surface.info[0].coorid[NETCDF_J_INDEX],
                         index,&(grid.lat[2*J+1]));
#else
      nc_put_var1_float(nc_id,var.gz_surface.info[0].coorid[NETCDF_J_INDEX],
                        index,&(grid.lat[2*J+1]));
#endif
    }
  }

  /*
   * Setup write_array().
   */
  num_nodes = setup_write_array();
  /*
   * Load start and end vectors.
   */
  start[0] = 1;
  end[  0] = grid.ni;

  for (node = 0; node < num_nodes; node++) {
    get_jlohi(node,num_nodes,start+1,end+1);
    write_array(node,TWODIM,start,end,0,var.gz_surface.info[0].name,
                var.gz_surface.info[0].index,var.gz_surface.value,EPIC_FLOAT_ARRAY,nc_id);
  }

  if (IAMNODE == NODE0) {
    /* Close file. */
    nc_close(nc_id);
  }

  return;
}

/*====================== end of init_gz_surface() ============================*/

/*====================== is_gz_surface_file() ================================*/

/*
 * Returns TRUE if entry->d_name contains the string "gz_surface.nc",
 * otherwise returns FALSE.
 */

int is_gz_surface_file(const struct dirent *entry) 
{
  if (strstr(entry->d_name,"gz_surface.nc")) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}

/*====================== end of is_gz_surface_file() =========================*/

/*====================== init_with_u() =======================================*/

/* 
 * Initialize the prognostic variables given meridional-plane
 * zonal wind profile and T(p) sounding.
 *
 * NOTE: Assumes bottom P is already assigned.
 * NOTE: Not MPI ready.
 */

#include <epic_pv_schemes.h>

#define GZ_FINE(iitp,J)    gz_fine[iitp][J-JLO]
#define FP_FINE(iitp,J)    fp_fine[iitp][J-JLO]
#define T_FINE(iitp,J)     t_fine[iitp][J-JLO]
#define THETA_FINE(iitp,J) theta_fine[iitp][J-JLO]
#define RHO_FINE(iitp,J)   rho_fine[iitp][J-JLO]
#define U_FINE(iitp,J)     u_fine[iitp][J-JLO]

#undef  MAX_ITER
#define MAX_ITER 100

void init_with_u(planetspec       *planet,
                 init_defaultspec *def)
{
  int
    K,J,I,kk,
    iitp,iitp0,ii,
    count,
    nbuff,
    theta_monotonic;
  EPIC_FLOAT
    pressure,temperature,density,
    fpara,mu,theta,
    fgibb,fpe,uoup,
    theta_ortho,theta_para,
    ptop,pbot,sigma,
    lat,
    x,dx,dy,dp,dgz,
    neglogp,h,tmp;
  EPIC_FLOAT
    *mudat,
    *rhodat,
    *thetadat,
    *fpdat,
    *udy,
    *pvhudy,
    *bern,
    **gz_fine,
    **theta_fine,
    **u_fine;
  float_triplet
    *buff_triplet;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_with_u";

  /*
   * Allocate memory:
   */
   mudat      = fvector(0,var.ntp-1,    dbmsname);
   rhodat     = fvector(0,var.ntp-1,    dbmsname);
   thetadat   = fvector(0,var.ntp-1,    dbmsname);
   fpdat      = fvector(0,var.ntp-1,    dbmsname);
   udy        = fvector(0,grid.nj+1,dbmsname);
   pvhudy     = fvector(0,grid.nj+1,dbmsname);
   bern       = fvector(0,grid.nj+1,dbmsname);

   buff_triplet = ftriplet(0,var.ntp,dbmsname);

   /*
    * Allocate 2D arrays for gz_fine, etc.
    */
   gz_fine    = (EPIC_FLOAT **)calloc(var.ntp,sizeof(EPIC_FLOAT *));
   if (!gz_fine) {
     sprintf(Message,"unable to allocate gz_fine");
     epic_error(dbmsname,Message);
   }
   theta_fine = (EPIC_FLOAT **)calloc(var.ntp,sizeof(EPIC_FLOAT *));
   if (!theta_fine) {
     sprintf(Message,"unable to allocate theta_fine");
     epic_error(dbmsname,Message);
   }
   u_fine = (EPIC_FLOAT **)calloc(var.ntp,sizeof(EPIC_FLOAT *));
   for (iitp = 0; iitp < var.ntp; iitp++) {
     gz_fine[iitp]    = fvector(0,JHI-JLO,dbmsname);
     theta_fine[iitp] = fvector(0,JHI-JLO,dbmsname);
     u_fine[iitp]     = fvector(0,JHI-JLO,dbmsname);
   }

  /*
   * Set fpdat.
   *
   * NOTE: For planets without hydrogen, fp (fraction of para hydrogen) has no effect. 
   */
  if (var.fpara.on) {
    /* Use equilibrium value, fpe(T).*/
    for (iitp = 0; iitp < var.ntp; iitp++) {
      fpdat[iitp] = return_fpe(var.tdat[iitp]);
    }
  }
  else {
    /* Use high-temperature value. */
    for (iitp = 0; iitp < var.ntp; iitp++) {
      fpdat[iitp] = 0.25;
    }
  }

  /*
   * Set thetadat.
   */
  for (iitp = 0; iitp < var.ntp; iitp++) {
    thetadat[iitp] = return_theta(planet,fpdat[iitp],var.pdat[iitp],var.tdat[iitp],&theta_ortho,&theta_para);
  }
  /*
   * Iron out thetadat from top down, adjust var.tdat accordingly.
   */
  for (iitp = var.ntp-2; iitp >= 0; iitp--) {
    thetadat[iitp] = MIN(thetadat[iitp],thetadat[iitp+1]*.999999);
    var.tdat[iitp] = return_temp(planet,fpdat[iitp],var.pdat[iitp],thetadat[iitp]);
  }

  /*
   * Set mudat.  The function mu_p() contains an estimate of the effects of 
   * optional condensables.
   */
  for (iitp = 0; iitp < var.ntp; iitp++) {
    mudat[iitp] = mu_p(var.pdat[iitp]);
  } 

  /*
   * Set rhodat.
   */
  for (iitp = 0; iitp < var.ntp; iitp++) {
    rhodat[iitp] = return_density(planet,fpdat[iitp],var.pdat[iitp],var.tdat[iitp],mudat[iitp],PASSING_T);
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   *                                                           *
   * 1. Specify T(p) at one latitude and determine GZ          *
   *    via hydrostatic balance.                               *
   *                                                           *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  /*
   * Determine fine (high vertical resolution) GZ at J = grid.jtp, I = ILO.
   * Find which iitp corresponds best to pbot, and call it iitp0.
   */
  tmp = FLOAT_MAX;
  for (iitp = 0; iitp < var.ntp; iitp++) {
    dp = fabs(var.pdat[iitp]-def->pbot);
    if (dp < tmp) {
      iitp0 = iitp;
      tmp   = dp;
    }
  }

  iitp = iitp0;
  J    = grid.jtp;
  I    = ILO;
  if (strcmp(planet->type,"terrestrial") == 0) {
    /*
     * The artificial parameter grid.topo_scale can be used
     * to suppress the topography.
     */
    GZ_FINE(iitp,J) = grid.topo_scale*GZ_SURFACE(J,I);
  }
  else if (strcmp(planet->type,"gas-giant") == 0) {
    GZ_FINE(iitp,J) = 0.;
  }
  else {
    sprintf(Message,"unrecognized planet type %s",planet->type);
    epic_error(dbmsname,Message);
  }

  /*
   * Vertically integrate the hydrostatic balance equation
   * down and up from iitp0. 
   *
   * NOTE: We tried the dgz/dlnp = p/rho form, but it sometimes generates problems like nonmonotonic GZ3,
   *       so we are sticking with this version. 
   */
  for (iitp = iitp0-1; iitp >= 0; iitp--) {
    GZ_FINE(iitp,J) = GZ_FINE(iitp+1,J)+.5*(1./rhodat[iitp]+1./rhodat[iitp+1])*(var.pdat[iitp+1]-var.pdat[iitp]);
  }
  for (iitp = iitp0+1; iitp < var.ntp; iitp++) {
    GZ_FINE(iitp,J) = GZ_FINE(iitp-1,J)-.5*(1./rhodat[iitp]+1./rhodat[iitp-1])*(var.pdat[iitp]-var.pdat[iitp-1]);
  }

  /* * * * * * * * * * * * * * * * * * * *
   *                                     *
   * 2. Specify target U_FINE.           *
   *                                     *
   * * * * * * * * * * * * * * * * * * * */

  if (def->u_scale == 0.) {
    for (iitp = 0; iitp < var.ntp; iitp++) {
      for (J = JLO; J <= JHI; J++) {
        U_FINE(iitp,J) = 0.;
      }
    }
  }
  else {
    for (iitp = 0; iitp < var.ntp; iitp++) {
      for (J = JLO; J <= JHI; J++) {
        U_FINE(iitp,J) = def->u_scale*planet->u(var.pdat[iitp],grid.lat[2*J+1]);
      }
    }
  }

  /*
   * Temporarily use the KLO plane of the U and PV arrays, in order to
   * employ the same algorithms for (zeta+f)*u*dy as in the DVDT code.
   */
  K = KLO;

  /*
   * Calculate the top level of GZ_FINE.
   */
  iitp = var.ntp-1;

  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      U(grid.it_uv,K,J,I) = U_FINE(iitp,J);
      V(grid.it_uv,K,J,I) = 0.;
    }
  }
  bc_lateral(var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,TWODIM);
  bc_lateral(var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,TWODIM);
  /* 
   * Calculate absolute vorticity.
   * Store in PV3 so that epic_pv_schemes.h macros work.
   * Use the same scheme as for dvdt in timestep().
   */
  vorticity(planet,ON_SIGMATHETA,ABSOLUTE,
            var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
            var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
            NULL,
            var.pv3.value+(K-Kshift)*Nelem2d);
  /* 
   * Form udy.
   */
  I = ILO;
  for (J = JLO; J <= JHI; J++) {
    dy     = 1./grid.n[2*J+1];
    udy[J] = U(grid.it_uv,K,J,I)*dy;
  }
  /*
   * Form (zeta+f)*u*dy = pvhudy.
   */
  for (J = JFIRST; J <= JHI; J++) {
    pvhudy[J] = (GA_V*udy[J  ]+DE_V*udy[J  ]
                +AL_V*udy[J-1]+BE_V*udy[J-1])*PV_COEF;
  }
  /*
   * Calculate bern by integrating -pvhudy.
   * Calculate GZ_FINE = bern-kin.
   */
  bern[grid.jtp] = GZ_FINE(iitp,grid.jtp)+get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                                         var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,grid.jtp,I);
  for (J = grid.jtp+1; J <= JHI; J++) {
    bern[J]         = bern[J-1]-pvhudy[J];
    GZ_FINE(iitp,J) = bern[J]-get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                             var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,J,I);
  }
  for (J = grid.jtp-1; J >= JLO; J--) {
    bern[J]         = bern[J+1]+pvhudy[J+1];
    GZ_FINE(iitp,J) = bern[J]-get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                             var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,J,I);
  }

  /*
   * Set THETA_FINE for iitp = var.ntp-1 to FLOAT_MAX to facilitate the non-monotonic logic.
   */
  for (J = JLO; J <= JHI; J++) {
    THETA_FINE(iitp,J) = FLOAT_MAX;
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
   *                                                                         *
   * 3. Proceed through the rest of the pressure levels.                     *
   * If a non-monotonic theta develops, relax the local zonal wind towards   *
   * no vertical shear and retry until theta is monotonic.                   *
   *                                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  for (iitp = var.ntp-2; iitp >= 0; iitp--) {
    for (count = 1; count <= MAX_ITER; count++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          U(grid.it_uv,K,J,I) = U_FINE(iitp,J);
        }
      }
      bc_lateral(var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,TWODIM);

      /*
       * Calculate absolute vorticity.
       * Store in PV3 so that epic_pv_schemes.h macros work.
       * Use the same scheme as for dvdt in timestep().
       */
      vorticity(planet,ON_SIGMATHETA,ABSOLUTE,
                var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                NULL,
                var.pv3.value+(K-Kshift)*Nelem2d);
      /*
       * Form udy. 
       */
      I = ILO;
      for (J = JLO; J <= JHI; J++) {
        dy     = 1./grid.n[2*J+1];
        udy[J] = U(grid.it_uv,K,J,I)*dy;
      }
      /*
       * Form (zeta+f)*u*dy = pvhudy.
       */
      for (J = JFIRST; J <= JHI; J++) {
        pvhudy[J] = (GA_V*udy[J  ]+DE_V*udy[J  ]
                    +AL_V*udy[J-1]+BE_V*udy[J-1])*PV_COEF;
      }
      /*
       * Calculate bern by integrating -pvhudy.
       * Calculate GZ_FINE = bern-kin.
       */
      bern[grid.jtp] = GZ_FINE(iitp,grid.jtp)+get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                                             var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,grid.jtp,I);
      for (J = grid.jtp+1; J <= JHI; J++) {
        bern[J]         = bern[J-1]-pvhudy[J];
        GZ_FINE(iitp,J) = bern[J]-get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                                 var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,J,I);
      }
      for (J = grid.jtp-1; J >= JLO; J--) {
        bern[J]         = bern[J+1]+pvhudy[J+1];
        GZ_FINE(iitp,J) = bern[J]-get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                                 var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,J,I);
      }

      /*
       * Obtain density from GZ via hydrostatic balance.
       * Get THETA_FINE(y,p) = theta(p,rho).
       */
      dp                 = var.pdat[iitp+1]-var.pdat[iitp];
      mu                 = .5*(mudat[iitp]+mudat[iitp+1]);
      fpara              = .5*(fpdat[iitp]+fpdat[iitp+1]);
      /*
       * NOTE: Using onto_kk() to get in-the-layer style averaging for pressure tends to cause
       *       a bias since in most of the code we interpolate with respect to log p.
       */
      pressure = sqrt(var.pdat[iitp]*var.pdat[iitp+1]);

      for (J = JLO; J <= JHI; J++) {
        dgz                = GZ_FINE(iitp+1,J)-GZ_FINE(iitp,J);
        density            = -dp/dgz;
        temperature        = alt_return_temp(planet,fpara,pressure,mu,density);
        THETA_FINE(iitp,J) = return_theta(planet,fpara,pressure,temperature,&theta_ortho,&theta_para);
      }

      /*
       * Check whether THETA_FINE is monotonic or not.
       */
      theta_monotonic = TRUE;
      for (J = JLO; J <= JHI; J++) {
        if (THETA_FINE(iitp,J) > THETA_FINE(iitp+1,J)) {
          /*
           * THETA_FINE is not monotonic.
           */
          theta_monotonic = FALSE;
          break;
        }
      }

      if (theta_monotonic) {
        break;
      }
      else {
        if (count == MAX_ITER) {
          break;
        }
        else if (count == MAX_ITER-1) {
          /*
           * Eliminate vertical shear for last iteration.
           */
          for (J = JLO; J <= JHI; J++) {
            U_FINE(iitp,J) = U_FINE(iitp+1,J);
          }
        }
        else {
          /*
           * Lessen vertical shear for next iteration.
           */
          for (J = JLO; J <= JHI; J++) {
            U_FINE(iitp,J) = 0.9*U_FINE(iitp,J)+0.1*U_FINE(iitp+1,J);
          }
        }
      }
    }
  }

  /*
   * Shift THETA_FINE by a half spacing so that all the pressure-coordinate variables are on the same level.
   * Recall THETA_FINE is currently shifted relative to the var.pdat levels because it is a function of density, 
   * and density is obtained by a vertical derivative of the geopotential.  The centering of this difference 
   * lands density, and hence THETA_FINE, between the var.pdat levels.
   */
  for (J = JLO; J <= JHI; J++) {
    for (ii = 0; ii <= var.ntp-2; ii++) {
      pressure           = sqrt(var.pdat[ii]*var.pdat[ii+1]);
      buff_triplet[ii].x = -log(pressure);
      buff_triplet[ii].y = THETA_FINE(ii,J);
    }
 
    spline_pchip(var.ntp-1,buff_triplet);

    /* 
     * Assume constant theta for bottom half-space shift: THETA_FINE(0,J) -> THETA_FINE(0,J).
     * Use monotonic spline to shift the rest.
     */
    for (iitp = 1; iitp <= var.ntp-2; iitp++) {
      neglogp            = -log(var.pdat[iitp]);
      ii                 = find_place_in_table(var.ntp-1,buff_triplet,neglogp,&h);
      THETA_FINE(iitp,J) = splint_pchip(neglogp,buff_triplet+ii,h);
    }
  }

  /*
   * We didn't shift the last THETA_FINE, so trim that level off and proceed.
   */
  nbuff = var.ntp-1;

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   *                                                           *
   * 4. The fine-resolution, pressure-coordinate variables are *
   *    now ready to be interpolated onto the hybrid coord.    *
   *                                                           *
   *    Interpolate P3, THETA onto sigmatheta surfaces.        *
   *    Fill in P2, etc                                        *
   *                                                           *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  fprintf(stderr,"Interpolating variables onto sigmatheta surfaces.\n");

  for (J = JLO; J <= JHI; J++) {
    lat = grid.lat[2*J+1];
    /*
     * Loop over I to handle terrain-following coordinate.
     *
     * Find P3 using sigmatheta as the independent variable.
     */
    for (I = ILO; I <= IHI; I++) {
      /*
       * Determine ptop that is consistent with the top sigmatheta value.
       * Use THETA_FINE as the independent variable, which should be monotonic
       * for altitudes above p_sigma.
       */
      count = 0;
      for (iitp = 0; iitp < nbuff; iitp++) {
        if (var.pdat[iitp] <= def->p_sigma) {
          buff_triplet[count].x = THETA_FINE(iitp,J);
          if (def->spacing_type == SPACING_P) {
            buff_triplet[count].y = var.pdat[iitp];
          }
          else {
            buff_triplet[count].y = -log(var.pdat[iitp]);
          }
          count++;
        }
      }
      spline_pchip(count,buff_triplet);

      K    = 0;
      x    = grid.sigmatheta[2*K+1];
      /*
       * Check for x out of range.
       */
      if (x <= buff_triplet[count-1].x) {
        iitp = find_place_in_table(count,buff_triplet,x,&dx);
      }
      else {
        sprintf(Message,"\n** Please rerun with ptop >= %.2g mbar, or extend \n** the top of t_vs_p.%s.",
                         var.pdat[nbuff-2]/100.,planet->name);
        epic_error(dbmsname,Message);
      }

      if (def->spacing_type == SPACING_P) {
        P3(K,J,I) = splint_pchip(x,buff_triplet+iitp,dx);
      }
      else {
        P3(K,J,I) = exp(-splint_pchip(x,buff_triplet+iitp,dx));
      }

      ptop = P3(0,  J,I);
      pbot = P3(KHI,J,I);

      /*
       * Start with p=pbot table entry.
       */
      buff_triplet[0].x = grid.zeta0;
      if (def->spacing_type == SPACING_P) {
        buff_triplet[0].y = pbot;
      }
      else {
        buff_triplet[0].y = -log(pbot);
      }
      count = 1;
      for (iitp = 0; iitp < nbuff; iitp++) {
        if (var.pdat[iitp] < pbot) {
          buff_triplet[count].x = return_sigmatheta(THETA_FINE(iitp,J),var.pdat[iitp],pbot,ptop);
          /*
           * Check that x is monotonic.
           */
          if (buff_triplet[count].x <= buff_triplet[count-1].x) {
            sprintf(Message,"J=%2d, iitp=%2d, sgth=%g %g not monotonically increasing (maybe try a different p_sigma)",
                             J,iitp,buff_triplet[count-1].x,buff_triplet[count].x);
            epic_error(dbmsname,Message);
          }
          if (def->spacing_type == SPACING_P) {
            buff_triplet[count].y = var.pdat[iitp];
          }
          else {
            buff_triplet[count].y = -log(var.pdat[iitp]);
          }
          count++;
        }
      }
      spline_pchip(count,buff_triplet);

      for (K = KLO-1; K < grid.k_sigma-1; K++) {
        x           = grid.sigmatheta[2*K+1];
        iitp        = find_place_in_table(count,buff_triplet,x,&dx);

        if (def->spacing_type == SPACING_P) {
          P3(K,J,I) = splint_pchip(x,buff_triplet+iitp,dx);
        }
        else {
          P3(K,J,I) = exp(-splint_pchip(x,buff_triplet+iitp,dx));
        }
      }

      /*
       * Fill in P3 in pure-sigma layers, where it is known diagnostically given
       * pbot, ptop and sigmatheta.
       */
      for (K = grid.k_sigma-1; K < KHI; K++) {
        sigma     = (grid.sigmatheta[2*K+1]-grid.zeta0)/(grid.zeta1-grid.zeta0);
        P3(K,J,I) = get_p_sigma(pbot,sigma,ptop);
      }

      /*
       * Set THETA to its diagnostic value.
       */
      for (K = KLO; K < grid.k_sigma-1; K++) {
        sigma        = get_sigma(pbot,P3(K,J,I),ptop);
        THETA(K,J,I) = (grid.sigmatheta[2*K+1]-f_sigma(sigma))/g_sigma(sigma);
      }

      /*
       * Interpolate to get THETA in the pure-sigma region, using P3 as the independent variable.
       */
      if (def->spacing_type == SPACING_P) {
        for (iitp = 0; iitp < nbuff; iitp++) {
          ii = nbuff-1-iitp;
          buff_triplet[iitp].x = var.pdat[ii];
          buff_triplet[iitp].y = THETA_FINE(ii,J);
        }
        spline_pchip(nbuff,buff_triplet);

        for (K = grid.k_sigma-1; K <= KHI; K++) {
          x            = P3(K,J,I);
          iitp         = find_place_in_table(nbuff,buff_triplet,x,&dx);
          THETA(K,J,I) = splint_pchip(x,buff_triplet+iitp,dx);
        }
      }
      else {
        for (iitp = 0; iitp < nbuff; iitp++) {
          /*
           * Need independent variable to be increasing for spline_pchip().
           */
          ii = nbuff-1-iitp;
          buff_triplet[iitp].x = log(var.pdat[ii]);
          buff_triplet[iitp].y = THETA_FINE(ii,J);
        }

        spline_pchip(nbuff,buff_triplet);

        for (K = grid.k_sigma-1; K <= KHI; K++) {
          x            = log(P3(K,J,I));
          iitp         = find_place_in_table(nbuff,buff_triplet,x,&dx);
          THETA(K,J,I) = splint_pchip(x,buff_triplet+iitp,dx);
        }
      } 
    } /* loop over I */
  } /* loop over J */

  bc_lateral(var.p3.value,   THREEDIM);
  bc_lateral(var.theta.value,THREEDIM);

  /* * * * * * * * * *
   *                 *
   * 5. Set HDRY, U  *
   *                 *
   * * * * * * * * * */

  /*
   * Set HDRY.
   */
  for (K = 1; K <= grid.nk; K++) {
    kk = 2*K;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        HDRY(K,J,I) = get_h(planet,kk,J,I,DRY);
      }
    }
    /* There is no need to apply bc_lateral(). */
  }
  restore_mass(planet,HDRY_INDEX,NO_PHASE);

  /*
   * Set P2, P3, H3, THETA2, etc.
   */
  set_p2etc(planet,UPDATE_THETA);

  /*
   * Interpolate to get U, using pressure as the independent variable.
   *
   * NOTE: This routine assumes zonal symmetry for the zonal wind,
   *       and is not appropriate for terrestrial initialization except
   *       in the case of zero-wind initial-condition experiments.
   */
  for (J = JLO; J <= JHI; J++) {
    if (def->spacing_type == SPACING_P) {
      for (iitp = 0; iitp < nbuff; iitp++) {
        /*
         * Need independent variable to be increasing for spline_pchip().
         */
        ii = nbuff-1-iitp;
        buff_triplet[iitp].x = var.pdat[ii];
        buff_triplet[iitp].y = U_FINE(ii,J);
      }
      spline_pchip(nbuff,buff_triplet);

      for (K = KLO-1; K <= KHI; K++) {
        x                     = P3(K,J,ILO);
        iitp                  = find_place_in_table(nbuff,buff_triplet,x,&dx);
        U(grid.it_uv,K,J,ILO) = splint_pchip(x,buff_triplet+iitp,dx);
        for (I = ILO+1; I <= IHI; I++) {
          U(grid.it_uv,K,J,I) = U(grid.it_uv,K,J,ILO);
        }
      }
    }
    else {
      for (iitp = 0; iitp < nbuff; iitp++) {
        ii = nbuff-1-iitp;
        buff_triplet[iitp].x = log(var.pdat[ii]);
        buff_triplet[iitp].y = U_FINE(ii,J);
      }

      spline_pchip(nbuff,buff_triplet);

      for (K = KLO-1; K <= KHI; K++) {
        x                     = log(P3(K,J,ILO));
        iitp                  = find_place_in_table(nbuff,buff_triplet,x,&dx);
        U(grid.it_uv,K,J,ILO) = splint_pchip(x,buff_triplet+iitp,dx);
        for (I = ILO+1; I <= IHI; I++) {
          U(grid.it_uv,K,J,I) = U(grid.it_uv,K,J,ILO);
        }
      }
    }
  }

  bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);

  /*
   * Free allocated memory:
   */
  free_ftriplet(buff_triplet,0,var.ntp,dbmsname);
  for (iitp = 0; iitp < var.ntp; iitp++) {
    free_fvector(theta_fine[iitp],0,JHI-JLO,dbmsname);
    free_fvector(gz_fine[iitp],   0,JHI-JLO,dbmsname);
    free_fvector(u_fine[iitp],    0,JHI-JLO,dbmsname);
  }
  free(theta_fine);
  free(gz_fine);
  free(u_fine);
  free_fvector(bern,      0,grid.nj+1,dbmsname);
  free_fvector(pvhudy,    0,grid.nj+1,dbmsname);
  free_fvector(udy,       0,grid.nj+1,dbmsname);
  free_fvector(fpdat,     0,var.ntp-1,    dbmsname);
  free_fvector(thetadat,  0,var.ntp-1,    dbmsname);
  free_fvector(rhodat,    0,var.ntp-1,    dbmsname);
  free_fvector(mudat,     0,var.ntp-1,    dbmsname);

  return;
}

/*====================== end of init_with_u() ================================*/

/*====================== init_with_hydrostatic() =============================*/

/*
 * Initialize the prognostic variables assuming U = 0, V = 0,  
 * p = p(z) and z = z(p). Assumes the surface pressure has already been set.
 *
 * In the sigma region, we know P3 on the interfaces, and hence we can use 
 * the mapping z = z(p), the inverse of the p = p(z) map used to generate
 * the surface pressure, to calculate GZ3.  Then, if possible, we use the precise 
 * algorithm the model employs for its hydrostatic integration and invert it to
 * generate the corresponding THETA. 
 *
 * In the hybrid region the problem is overconstrained, so our approach
 * is approximate.  
 */

void init_with_hydrostatic(planetspec *planet) 
{
  int
    K,J,I,
    kk,
    error_flag,
    count;
  EPIC_FLOAT
    pbot,ptop,sigma,
    plo,phi,ptol,
    fpara,mu,press,
    theta_ortho,theta_para,
    fgibb,fpe,uoup,
    p_guess,gz_guess;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_with_hydrostatic";

  /* * * * * * * * *
   * Sigma region. *
   * * * * * * * * */

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      K = 0;
      ptop = P3(K,J,I) = grid.p_ref[2*K+1];
      K = grid.nk;
      pbot = P3(K,J,I);
      /*
       * Set P3 in sigma portion of model using diagnostic value.
       */
      for (K = grid.nk-1; K >= grid.k_sigma-1; K--) {
        sigma     = (grid.sigmatheta[2*K+1]-grid.zeta0)/(grid.zeta1-grid.zeta0);
        P3(K,J,I) = get_p_sigma(pbot,sigma,ptop);
      }

      /*
       * Pick P3 at K = grid.k_sigma-2, the value at the interface just above the seam.
       */
      K = grid.k_sigma-2;
      P3(K,J,I) = SQR(P3(K+1,J,I))/P3(K+2,J,I);
    }
  }
  /* No need to apply bc_lateral() here. */

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      for (K = grid.nk; K >= grid.k_sigma-2; K--) {
        GZ3(K,J,I) = gz_p(planet,J,P3(K,J,I));
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  /*
   * Determine RHO3. 
   * If possible, the algorithm should match the discrete form used for the hydrostatic integration.
   * Given P3 and RHO3, determine first T3 and then THETA.
   */

  /*
   * Surface value.
   */
  K = grid.nk;
  if (strcmp(planet->type,"gas-giant") == 0) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        /*
         * NOTE: Assuming fpara and molar mass are accurately known at bottom of model at this point,
         *       which may not be true.
         */
        fpara       = (var.fpara.on) ? FPARA(K,J,I) : 0.25;
        mu          = avg_molar_mass(planet,2*K+1,J,I);
        RHO3(K,J,I) = return_density(planet,fpara,grid.p_ref[2*K+1],grid.theta_ref[2*K+1],mu,PASSING_THETA);
      }
    }
  }
  else {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        RHO3(K,J,I) = rho_gz(planet,J,GZ3(K,J,I));
      }
    }
  }

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      for (K = KHI-1; K >= grid.k_sigma-1; K--) {
        RHO3(K,J,I) = rho_gz(planet,J,GZ3(K,J,I));
      }
      for (K = KHI; K >= grid.k_sigma-1; K--) {
        fpara        = (var.fpara.on) ? FPARA(K,J,I) : 0.25;
        mu           = avg_molar_mass(planet,2*K+1,J,I);
        T3(K,J,I)    = alt_return_temp(planet,fpara,P3(K,J,I),mu,RHO3(K,J,I));
        THETA(K,J,I) = return_theta(planet,fpara,P3(K,J,I),T3(K,J,I),&theta_ortho,&theta_para);
      }
    }
  }

  /* * * * * * * * * *
   * Hybrid region.  *
   * * * * * * * * * */

  /*
   * Unlike the sigma region, the goal of achieving p = p(z) exactly on the grid appears
   * to be overcontrained in the hybrid region, so we settle for an approximate approach.
   */

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /*
       * P3 at K = grid.k_sigma-2 is set above.
       * Determine corresponding THETA.
       */
      K = grid.k_sigma-2;
      sigma        = get_sigma(P3(KHI,J,I),P3(K,J,I),P3(0,J,I));
      THETA(K,J,I) = (EPIC_FLOAT)(((double)grid.sigmatheta[2*K+1]-(double)f_sigma(sigma))/(double)g_sigma(sigma));

      for (K = grid.k_sigma-3; K >= KLO; K--) {
        /*
         * Set P3 and determine corresponding THETA.
         */
        P3(K,J,I)    = grid.p_ref[2*K+1];
        sigma        = get_sigma(P3(KHI,J,I),P3(K,J,I),P3(0,J,I));
        THETA(K,J,I) = (EPIC_FLOAT)(((double)grid.sigmatheta[2*K+1]-(double)f_sigma(sigma))/(double)g_sigma(sigma));
      }

      /*
       * Top of model.
       */
      K = 0;
      THETA(K,J,I) = grid.sigmatheta[2*K+1];
    }
  }

  /*
   * Set HDRY. Only need P3 for get_h() here.
   */
  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        HDRY(K,J,I) = get_h(planet,kk,J,I,DRY);
      }
    }
  }

  /*
   * Start with U = 0, V = 0, which requires no action.
   */

  return;
}

/*====================== end of init_with_hydrostatic() ======================*/

/*====================== init_with_ref() =====================================*/
    
/* 
 * Initializes P, H using grid.p_ref[kk] and THETA using grid.theta_ref[kk],
 * but make use of the diagnostic values where they are known.
 * Assumes the surface pressure has already been set.
 */

void init_with_ref(planetspec *planet)
{
  int
    K,J,I,
    kk;
  EPIC_FLOAT
    pbot,ptop,sigma;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_with_ref";

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      K = 0;
      ptop = P3(K,J,I) = grid.p_ref[2*K+1];
      pbot = P3(grid.nk,J,I);
      /*
       * Set P3 in sigma portion of model using diagnostic value.
       */
      for (K = grid.nk-1; K >= grid.k_sigma-1; K--) {
        sigma     = (grid.sigmatheta[2*K+1]-grid.zeta0)/(grid.zeta1-grid.zeta0);
        P3(K,J,I) = get_p_sigma(pbot,sigma,ptop);
      }
      /*
       * Set P3 in hybrid portion of model using the reference value.
       * Take care that P3 decreases with height in the transition zone.
       */
      for (K = grid.k_sigma-2; K > 0; K--) {
        P3(K,J,I) = MIN(grid.p_ref[2*K+1],.99999*P3(K+1,J,I));
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  /*
   * Set HDRY. Only need P3 for get_h() here.
   */
  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        HDRY(K,J,I) = get_h(planet,kk,J,I,DRY);
      }
    }
  }

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /*
       * Set THETA in the sigma portion of the model using the reference value.
       */
      for (K = grid.k_sigma-1; K <= KHI; K++) {
        THETA(K,J,I) = grid.theta_ref[2*K+1];
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  /*
   * Set P2, P3, H3, THETA2, etc.
   */
  set_p2etc(planet,UPDATE_THETA);

  /*
   * Start with U = 0, V = 0, which requires no action.
   */

  return;
}

/*====================== end of init_with_ref() ==============================*/

/*====================== init_fpara_as_fpe() =================================*/

/*
 * Iterate to get fpara = fpe(temp(fpara,pressure,theta)):
 */

#undef  MAX_IT
#define MAX_IT 10

void init_fpara_as_fpe(planetspec *planet)
{
  int
    K,J,I,
    it,error_flag;
  EPIC_FLOAT
    fpara,fptol,fp1,fp2,
    fpe;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_fpara_as_fpe";

  /* Check that fpara is turned on. */
  if (!var.fpara.on) {
    fprintf(stderr,"Warning: init_fpara_as_fpe(): var.fpara.on is off \n");
    return;
  }

  FPEMFPE_planet = planet;

  fptol = pow(machine_epsilon(),2./3.);

  for (K = KLO-1; K <= KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        FPEMFPE_p     = P3(   K,J,I);
        FPEMFPE_theta = THETA(K,J,I);

        /* Initial guess. */
        fpara = FPARA(K,J,I);
        fp1   = FPARA(K,J,I)*0.9;
        fp2   = FPARA(K,J,I)*1.1;

        for (it = 0; it < MAX_IT; it++) {
          error_flag = find_root(fp1,fp2,fptol,&fpe,fpe_minus_fpe);
          if (error_flag == 0) {
            /* Convergence. */
            break;
          }
          /* Try a wider interval. */
          fp1 *= .5;
          fp2 *= 2.;
        }
        if (it == MAX_IT) {
          sprintf(Message,"exceeded MAX_IT = %d",MAX_IT);
          epic_error(dbmsname,Message);
        }

        FPARA(K,J,I) = fpe;
      }
    }
  }

  return;
}

/*====================== end of init_fpara_as_fpe() ==========================*/

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

/*====================== init_species() ======================================*/

/*
 * Top function for initializing optional species (condensables) in the atmosphere.
 *
 * NOTE: Need to calculate required diagnostic variables like T3 here,
 *       since this function is called before they are set.
 *
 */

void init_species(planetspec       *planet,
                  init_defaultspec *def,
                  int               prompt_mode)
{
  int
    K,J,I,
    is,iq;
  EPIC_FLOAT
    factor,
    mu_dry,
    pressure,theta,fpara,temperature;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_species";

  /*
   * Make sure pressure-related diagnostic variables are set.
   */
  set_p2etc(planet,DONT_UPDATE_THETA);

  /*
   * Calculate T3.
   */
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      for (K = KLO; K <= KHI; K++) {
        pressure = P3(K,J,I);
        theta    = THETA(K,J,I);
        if (var.fpara.on) {
          fpara = FPARA(K,J,I);
        }
        else {
          fpara = 0.25;
        }
        T3(K,J,I) = temperature = return_temp(planet,fpara,pressure,theta);
      }
    }
  }

  /*
   * Use Q initially to store volume mixing ratios.
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      switch (is) {
        case C_2H_2_INDEX:
          init_vmr_via_data(planet,C_2H_2_INDEX);
        break;
        case C_2H_6_INDEX:
          init_vmr_via_data(planet,C_2H_6_INDEX);
        break;
        default:
          init_vmr_via_deep_value(planet,is,&def->mole_fraction[is],
                                  &def->mole_fraction_over_solar[is],&def->rh_max[is],prompt_mode);
        break;
      }
    }
  }

  /*
   * Convert from volume mixing ratio (vmr, aka mole fraction) to mass mixing ratio.
   */
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      for (K = KHI; K >= KLO; K--) {
        mu_dry = R_GAS/planet->rgas;
        factor  = 1.;
        for (iq = 0; iq < grid.nq; iq++) {
          factor -= Q(grid.is[iq],grid.ip[iq],K,J,I);
        }
        factor = 1./(mu_dry*factor);
        for (iq = 0; iq < grid.nq; iq++) {
          Q(grid.is[iq],grid.ip[iq],K,J,I) *= var.species[grid.is[iq]].molar_mass*factor;
        }
      }

      /*
       * Calculate HDRY, HDRY3 to accomodate the added species while keeping H2, H3 constant.
       * This maintains consistency with the use of total pressure.
       */
      for (K = KLO; K <= KHI; K++) {
        /*
         * Calculate HDRY.
         */
        factor = 1.;
        for (iq = 0; iq < grid.nq; iq++) {
          factor += .5*(Q(grid.is[iq],grid.ip[iq],K,J,I)+Q(grid.is[iq],grid.ip[iq],K-1,J,I));
        }
        HDRY(K,J,I) = H2(K,J,I)/factor;

        /*
         * Calculate HDRY3.
         */
        factor = 1.;
        for (iq = 0; iq < grid.nq; iq++) {
          factor += Q(grid.is[iq],grid.ip[iq],K,J,I);
        }
        HDRY3(K,J,I) = H3(K,J,I)/factor;
      }
    }
  }

  /*
   * Update pressure-related diagnostic variables.
   */
  set_p2etc(planet,UPDATE_THETA);

  return;
}

/*====================== end of init_species() ===============================*/

/*====================== init_vmr_via_deep_value() ===========================*/

/*
 * Assign volume mixing ratio (vmr, aka mole fraction) starting with a deep value,
 * which is assumed to be bounded by the saturation curve and a cold trap effect.
 * Stores in Q (converted to mass mixing ratio elsewhere).
 */

void init_vmr_via_deep_value(planetspec *planet,
                             int         is,
                             EPIC_FLOAT *mole_fraction,
                             EPIC_FLOAT *mole_fraction_over_solar,
                             EPIC_FLOAT *rh_max,
                             int         prompt_mode)
{
  register int
    K,J,I;
  register EPIC_FLOAT
    solar,
    sat_vapor_p,
    x_sat,
    mole_frac0;
  char
    min_element[4],
    prompt[64];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_vmr_via_deep_value";

  /*
   * Establish deep-atmosphere mole fractions and store in *mole_fraction.
   */
  if (prompt_mode == USE_PROMPTS) {
    solar = solar_fraction(var.species[is].info[0].name,MOLAR,min_element);
    fprintf(stdout,"Solar mole fraction of %s is %e \n",var.species[is].info[0].name,solar);
    sprintf(prompt,"Input %s mole fraction relative to solar [1.=solar]:\n",var.species[is].info[0].name);
    *mole_fraction_over_solar = input_float(prompt,*mole_fraction_over_solar);
    *mole_fraction            = *mole_fraction_over_solar*solar;
    *rh_max                   = input_float("Input maximum initial relative humidity [1.=sat]\n",*rh_max);
  }

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /* 
       * Start at lowest layer with mole fractions set to specified value. 
       */
      mole_frac0 = *mole_fraction;

      for (K = KHI; K >= KLO; K--) {
        /* 
         * Calculate saturation mole fraction.
         * Assume the mole fraction is bounded by the saturation curve,
         * and store it in the mass variable, Q.
         */
        sat_vapor_p       = var.species[is].sat_vapor_p(T3(K,J,I));
        x_sat             = sat_vapor_p/P3(K,J,I);
        Q(is,VAPOR,K,J,I) = MIN(mole_frac0,x_sat*(*rh_max));

        /*
         * Assume that as we travel upwards, mole fraction stays 
         * reduced once it has been limited by saturation (cold trap effect).
         */
        mole_frac0 = Q(is,VAPOR,K,J,I);
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  return;
}

/*====================== end of init_vmr_via_deep_value() ====================*/

/*====================== init_vmr_via_data() =================================*/

/*
 * Assign volume mixing ratio (vmr, aka mole fraction) using data from a file.
 * Stores in Q (converted to mass mixing ratio elsewhere).
 */

void init_vmr_via_data(planetspec *planet,
                       int         is)
{
  int
    K,J,I;
  EPIC_FLOAT
    p,lat,
    temp,C_2H_2,C_2H_6;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_vmr_via_data";

  switch (is) {
    case C_2H_2_INDEX:
      for (K = KLO; K <= KHI; K++) {
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            HCSdata(planet,P3(K,J,I),grid.lat[2*J+1],&temp,&C_2H_2,&C_2H_6);
            Q(is,VAPOR,K,J,I) = C_2H_2;
          }
        }
      }
      /* Need to apply bc_lateral() here */
      bc_lateral(var.species[is].phase[VAPOR].q,THREEDIM);
    return;
    case C_2H_6_INDEX:
      for (K = KLO; K <= KHI; K++) {
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            HCSdata(planet,P3(K,J,I),grid.lat[2*J+1],&temp,&C_2H_2,&C_2H_6);
            Q(is,VAPOR,K,J,I) = C_2H_6;
          }
        }
      }
      /* Need to apply bc_lateral() here */
      bc_lateral(var.species[is].phase[VAPOR].q,THREEDIM);
    return;
    default:
      sprintf(Message,"unrecognized species index %d",is);
      epic_error(dbmsname,Message);
    return;
  }
}

/*====================== end of init_vmr_via_data() ==========================*/

/*====================== HCSdata() ===========================================*/

/*
 * Specific function to return values of temperature, and volume mixing ratios
 * of C_2H_2 and C_2H_6 from meridional-plane Cassini HCS data, given the
 * pressure [Pa] and latitude [deg planetographic].
 * Uses bilinear interpolation.
 * Returns the nearest table value if outside the range of table.
 */

#define DATA_TEMPS(k,j)  data_temps[ j+(k)*nlat_temps ]
#define DATA_C_2H_2(k,j) data_C_2H_2[j+(k)*nlat_C_2H_2]
#define DATA_C_2H_6(k,j) data_C_2H_6[j+(k)*nlat_C_2H_6]

void HCSdata(planetspec *planet,
             EPIC_FLOAT  p,
             EPIC_FLOAT  lat,
             EPIC_FLOAT  *temp,
             EPIC_FLOAT  *C_2H_2,
             EPIC_FLOAT  *C_2H_6)
{
  int
    j,k;
  EPIC_FLOAT
    frac_lat,frac_p;
  static EPIC_FLOAT
    *p_temps,*p_C_2H_2,*p_C_2H_6,
    *lat_temps,*lat_C_2H_2,*lat_C_2H_6,
    *data_temps,*data_C_2H_2,*data_C_2H_6;
  static int
    np_temps,np_C_2H_2,np_C_2H_6,
    nlat_temps,nlat_C_2H_2,nlat_C_2H_6,
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="HCSdata()";

  if (!initialized) {
    if (strcmp(planet->name,"jupiter") == 0) {
      /* Determine size of data files. */
      read_meridional_plane(planet,EPIC4_PATH"/data/jupiter/temps.jupiter.HCS",SIZE_DATA,
                            &np_temps,&nlat_temps,NULL,NULL,NULL);
      read_meridional_plane(planet,EPIC4_PATH"/data/jupiter/C_2H_2.jupiter.HCS",SIZE_DATA,
                            &np_C_2H_2,&nlat_C_2H_2,NULL,NULL,NULL);
      read_meridional_plane(planet,EPIC4_PATH"/data/jupiter/C_2H_6.jupiter.HCS",SIZE_DATA,
                            &np_C_2H_6,&nlat_C_2H_6,NULL,NULL,NULL);

      /* Allocate memory */
      p_temps     = fvector(0,np_temps   -1,dbmsname);
      p_C_2H_2    = fvector(0,np_C_2H_2  -1,dbmsname);
      p_C_2H_6    = fvector(0,np_C_2H_6  -1,dbmsname);
      lat_temps   = fvector(0,nlat_temps -1,dbmsname);
      lat_C_2H_2  = fvector(0,nlat_C_2H_2-1,dbmsname);
      lat_C_2H_6  = fvector(0,nlat_C_2H_6-1,dbmsname);
      data_temps  = fvector(0,np_temps*nlat_temps  -1,dbmsname);
      data_C_2H_2 = fvector(0,np_C_2H_2*nlat_C_2H_2-1,dbmsname);
      data_C_2H_6 = fvector(0,np_C_2H_6*nlat_C_2H_6-1,dbmsname);

      /* Read in data */
      read_meridional_plane(planet,EPIC4_PATH"/data/jupiter/temps.jupiter.HCS",POST_SIZE_DATA,
                            &np_temps,&nlat_temps,p_temps,lat_temps,data_temps);
      read_meridional_plane(planet,EPIC4_PATH"/data/jupiter/C_2H_2.jupiter.HCS",POST_SIZE_DATA,
                            &np_C_2H_2,&nlat_C_2H_2,p_C_2H_2,lat_C_2H_2,data_C_2H_2);
      read_meridional_plane(planet,EPIC4_PATH"/data/jupiter/C_2H_6.jupiter.HCS",POST_SIZE_DATA,
                            &np_C_2H_6,&nlat_C_2H_6,p_C_2H_6,lat_C_2H_6,data_C_2H_6);
    }
    else {
      sprintf(Message,"data not available for %s",planet->name);
      epic_error(dbmsname,Message);
    }

    initialized = TRUE;
  }

  /*
   * Use bilinear interpolation to get values.
   * First, temperature.
   */
  if (lat < lat_temps[0]) {
    j        = 0;
    frac_lat = 0.;
  }
  else if (lat > lat_temps[nlat_temps-1]) {
    j        = nlat_temps-2;
    frac_lat = 1.;
  }
  else {
    for (j = 1; j < nlat_temps; j++) {
      if (lat <= lat_temps[j]) {
        j--;
        frac_lat = (lat-lat_temps[j])/(lat_temps[j+1]-lat_temps[j]);
        break;
      }
    }  
  }

  if (p > p_temps[0]) {
    k      = 0;
    frac_p = 0.;
  }
  else if (p < p_temps[np_temps-1]) {
    k      = np_temps-2;
    frac_p = 1.;
  }
  else {
    for (k = 1; k < np_temps; k++) {
      if (p >= p_temps[k]) {
        k--;
        frac_p = (p-p_temps[k])/(p_temps[k+1]-p_temps[k]);
        break;
      }
    }
  }  

  *temp = DATA_TEMPS(k,  j  )*(1.-frac_p)*(1.-frac_lat)
         +DATA_TEMPS(k+1,j  )*(   frac_p)*(1.-frac_lat)
         +DATA_TEMPS(k,  j+1)*(1.-frac_p)*(   frac_lat)
         +DATA_TEMPS(k+1,j+1)*(   frac_p)*(   frac_lat);

  /*
   * C_2H_2
   */
  if (lat < lat_C_2H_2[0]) {
    j        = 0;
    frac_lat = 0.;
  }
  else if (lat > lat_C_2H_2[nlat_C_2H_2-1]) {
    j        = nlat_C_2H_2-2;
    frac_lat = 1.;
  }
  else {
    for (j = 1; j < nlat_C_2H_2; j++) {
      if (lat <= lat_C_2H_2[j]) {
        j--;
        frac_lat = (lat-lat_C_2H_2[j])/(lat_C_2H_2[j+1]-lat_C_2H_2[j]);
        break;
      }
    }  
  }

  if (p > p_C_2H_2[0]) {
    k      = 0;
    frac_p = 0.;
  }
  else if (p < p_C_2H_2[np_C_2H_2-1]) {
    k      = np_C_2H_2-2;
    frac_p = 1.;
  }
  else {
    for (k = 1; k < np_C_2H_2; k++) {
      if (p >= p_C_2H_2[k]) {
        k--;
        frac_p = (p-p_C_2H_2[k])/(p_C_2H_2[k+1]-p_C_2H_2[k]);
        break;
      }
    }
  }  

  *C_2H_2 = DATA_C_2H_2(k,  j  )*(1.-frac_p)*(1.-frac_lat)
           +DATA_C_2H_2(k+1,j  )*(   frac_p)*(1.-frac_lat)
           +DATA_C_2H_2(k,  j+1)*(1.-frac_p)*(   frac_lat)
           +DATA_C_2H_2(k+1,j+1)*(   frac_p)*(   frac_lat);

  /*
   * C_2H_6
   */
  if (lat < lat_C_2H_6[0]) {
    j        = 0;
    frac_lat = 0.;
  }
  else if (lat > lat_C_2H_6[nlat_C_2H_6-1]) {
    j        = nlat_C_2H_6-2;
    frac_lat = 1.;
  }
  else {
    for (j = 1; j < nlat_C_2H_6; j++) {
      if (lat <= lat_C_2H_6[j]) {
        j--;
        frac_lat = (lat-lat_C_2H_6[j])/(lat_C_2H_6[j+1]-lat_C_2H_6[j]);
        break;
      }
    }  
  }

  if (p > p_C_2H_6[0]) {
    k      = 0;
    frac_p = 0.;
  }
  else if (p < p_C_2H_6[np_C_2H_6-1]) {
    k      = np_C_2H_6-2;
    frac_p = 1.;
  }
  else {
    for (k = 1; k < np_C_2H_6; k++) {
      if (p >= p_C_2H_6[k]) {
        k--;
        frac_p = (p-p_C_2H_6[k])/(p_C_2H_6[k+1]-p_C_2H_6[k]);
        break;
      }
    } 
  } 

  *C_2H_6 = DATA_C_2H_6(k,  j  )*(1.-frac_p)*(1.-frac_lat)
           +DATA_C_2H_6(k+1,j  )*(   frac_p)*(1.-frac_lat)
           +DATA_C_2H_6(k,  j+1)*(1.-frac_p)*(   frac_lat)
           +DATA_C_2H_6(k+1,j+1)*(   frac_p)*(   frac_lat);

  return;
}

/*====================== end of HCSdata() ====================================*/

/*====================== setup_mu_p() ========================================*/

/*
 * Global declarations for mu_p().
 */
static float_triplet
  *mu_p_table;
static int
   n_mu_p;

/*
 * Assumes a preliminary call to init_with_ref() and init_species() has been made. 
 * Used to establish an approximate function relating molar mass to pressure
 * to improve init_with_u()'s gradient balance for models with condensables.
 */
void setup_mu_p(planetspec *planet)
{
  int
    K,kk,ii;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="setup_mu_p";

  n_mu_p = 2*grid.nk+1;

  /*
   * Allocate memory.
   */
  mu_p_table = ftriplet(0,n_mu_p-1,dbmsname);

  K  = 0;
  kk = 2*K+1;
  ii = 0;
  mu_p_table[ii].x = P3(K,JLO,ILO);
  mu_p_table[ii].y = avg_molar_mass(planet,kk,JLO,ILO);
  for (K = 1; K <= KHI; K++) {
    kk = 2*K;
    ii++;
    mu_p_table[ii].x = P2(K,JLO,ILO);
    mu_p_table[ii].y = avg_molar_mass(planet,kk,JLO,ILO);

    kk++;
    ii++;
    mu_p_table[ii].x = P3(K,JLO,ILO);
    mu_p_table[ii].y = avg_molar_mass(planet,kk,JLO,ILO);
  }
  spline_pchip(n_mu_p,mu_p_table);

  return;
}

/*====================== end of setup_mu_p() =================================*/

/*====================== mu_p() ==============================================*/

/*
 * Relate molar mass to pressure to help improve init_with_u()'s 
 * gradient balance for models with condensables.
 *
 * Requires initial call to setup_mu_p().
 */

EPIC_FLOAT mu_p(EPIC_FLOAT pressure)
{
  int
    i;
  EPIC_FLOAT
    mu,p,dp;

  p  = LIMIT_RANGE(mu_p_table[0].x,pressure,mu_p_table[n_mu_p-1].x);
  i  = find_place_in_table(n_mu_p,mu_p_table,p,&dp);
  mu = splint_pchip(p,mu_p_table+i,dp);

  return mu;
}

/*====================== end of mu_p() =======================================*/

/*====================== t_yp() ==============================================*/

/*
 * Interpolates T(y,p) data table using a bicubic spline.
 */

#define T_DATA(k,j)    t_data[   j+(k)*nlat]
#define T_DATA_Y2(k,j) t_data_y2[j+(k)*nlat]

EPIC_FLOAT t_yp(EPIC_FLOAT        p, 
                EPIC_FLOAT        latr, 
                int               mode,
                init_defaultspec *def)
{
  char
    string[N_STR];
  int
    k,j;
  static int
    initialized=FALSE,
    nlat,
    npress;
  EPIC_FLOAT 
    z,temp,lat_d,
    log_p,p_d,
    a4 = 2.;
  static EPIC_FLOAT
    *lat,
    *press,
    *t_data,
    *t_data_y2,
    *t_k,
    *t_k_y2;
  static float_triplet
    *buff_triplet;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="t_yp";

  if (!initialized) {
    char
      header[80];
    FILE
      *infile;

    sprintf(string,EPIC4_PATH"/data/%s/iris/temperature.dat",planet->name);
    infile = fopen(string,"r");
    if (!infile) {
      fprintf(stderr,"Cannot find input file %s \n",string);
      exit(1);
    }
    fscanf(infile,"%d %d",&nlat,&npress);
    fscanf(infile,"%s",header);
    /* Allocate memory: */
    lat          = fvector( 0,nlat-1,             dbmsname);
    press        = fvector( 0,npress-1,           dbmsname);
    t_k          = fvector( 0,npress-1,           dbmsname);
    t_k_y2       = fvector( 0,npress-1,           dbmsname);
    t_data       = fvector( 0,nlat*npress-1,      dbmsname);
    t_data_y2    = fvector( 0,nlat*npress-1,      dbmsname);
    buff_triplet = ftriplet(0,IMAX(nlat,npress)-1,dbmsname);

    for (k = 0; k < npress; k++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(infile,"%lf",press+k);
#else
      fscanf(infile,"%f",press+k);
#endif

      /* convert press to mks */
      press[k] *= 100.;
      /* spline on log p */
      press[k] = log(press[k]);
    }
    fscanf(infile,"%s",header);
    for (j = 0; j < nlat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(infile,"%lf",lat+j);
#else
      fscanf(infile,"%f",lat+j);
#endif

      /* convert lat to radians */
      lat[j] *= DEG;
      for (k = 0; k < npress; k++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile,"%lf",&(T_DATA(k,j)));
#else
        fscanf(infile,"%f",&(T_DATA(k,j)));
#endif

      }
    }
    fclose(infile);
    /* Compute spline information for each row: */

    for (j = 0; j < nlat; j++) {
      buff_triplet[j].x = lat[j];
    }

    for (k = 0; k < npress; k++) {
      /*
       * Pack float_triplet buffer.
       */
      for (j = 0; j < nlat; j++) {
        buff_triplet[j].y = T_DATA(k,j);
      }
      spline_pchip(nlat,buff_triplet);
      /* 
       * Unpack y2. 
       */
      for (j = 0; j < nlat; j++) {
        T_DATA_Y2(k,j) = buff_triplet[j].z;
      }
    }
    if (IAMNODE == 0) {
      /* Echo range limits: */
      fprintf(stderr,"Range limits: lat[%2.1f,%2.1f]  p[%2.1f,%2.1f] \n",
                     lat[0]/DEG,lat[nlat-1]/DEG,
                     exp(press[npress-1])/100.,exp(press[0])/100.);
    }

    if (def) {
      /* Set default ranges */
      /* latbot, lattop in degrees */
      def->globe_latbot = lat[0     ]/DEG;
      def->globe_lattop = lat[nlat-1]/DEG;
      /* ptop, pbot in mbar */
      def->ptop         = exp(press[0       ])/100.;
      def->pbot         = exp(press[npress-1])/100.;
    }

    initialized = TRUE;
  }
  /* End initialization. */

  // if (p == 0) return;

  if (mode == 1) {
    /* Synthetic data: */
    z    = -log(p/grid.press0);
    temp = a4*z*cos(4.*latr);
  }
  else if (mode == 0) {
    /* Data */
    if (latr < lat[0]) {
      sprintf(Message,"lat = %f below data range = %f",latr/DEG,lat[0]/DEG);
      epic_error(dbmsname,Message);
    }
    else if (latr > lat[nlat-1]) {
      sprintf(Message,"lat = %f above data range = %f",latr/DEG,lat[nlat-1]/DEG);
      epic_error(dbmsname,Message);
    }
    for (j = 0; j < nlat; j++) {
      buff_triplet[j].x = lat[j];
    }
    j = find_place_in_table(nlat,buff_triplet,latr,&lat_d);
    /* Make column at latr */
    for (k = 0; k < npress; k++) {
      /*
       * Pack float_triplet buffer segment.
       */
      buff_triplet[j  ].y = T_DATA(   k,j  );
      buff_triplet[j+1].y = T_DATA(   k,j+1);
      buff_triplet[j  ].z = T_DATA_Y2(k,j  );
      buff_triplet[j+1].z = T_DATA_Y2(k,j+1);

      t_k[k] = splint_pchip(latr,buff_triplet+j,lat_d);
    }
    /*
     * Pack float_triplet buffer.
     */
    for (k = 0; k < npress; k++) {
      buff_triplet[k].x = press[k];
      buff_triplet[k].y = t_k[k];
    }
    spline_pchip(npress,buff_triplet);

    log_p = log(p);
    if (log_p >= press[0] && log_p <= press[npress-1]) {
      k    = find_place_in_table(npress,buff_triplet,log_p,&p_d);
      temp = splint_pchip(log_p,buff_triplet+k,p_d);
    }
    else if (log_p > press[0]) {
      sprintf(Message,"p(%5.1f) = %e out of bounds (increase floor_tp)",
                      latr/DEG,p/100.);
      epic_error(dbmsname,Message);
    }
    else {
      sprintf(Message,"p(%5.1f) = %e out of bounds (decrease ceiling_tp)",
                      latr/DEG,p/100.);
      epic_error(dbmsname,Message);
    }
  }
  else {
    sprintf(Message,"unrecognized mode");
    epic_error(dbmsname,Message);
  }

  return temp;
}

/*====================== end of t_yp() ========================================*/

/*====================== fp_yp() ==============================================*/

/*
 * Interpolates fp data table using a bicubic spline.
 * Uses data-boundary value when pressure out of range.
 */

#define FP_DATA(k,j)    fp_data[   j+(k)*nlat]
#define FP_DATA_Y2(k,j) fp_data_y2[j+(k)*nlat]

EPIC_FLOAT fp_yp(EPIC_FLOAT p, 
                 EPIC_FLOAT latr, 
                 int        mode)
{
  char
    string[N_STR];
  int
    k,j;
  static int
    initialized=FALSE,
    nlat,
    npress;
  EPIC_FLOAT 
    z,fp,lat_d,
    log_p,p_d,
    a1 =  0.5,
    a2 =  0.1,
    a3 = -0.05;
  static EPIC_FLOAT
    *lat,
    *press,
    *fp_data,
    *fp_data_y2,
    *fp_k,
    *fp_k_y2;
  static float_triplet
    *buff_triplet;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="fp_yp";

  if (!initialized) {
    char
      header[80];
    FILE
      *infile;

    sprintf(string,EPIC4_PATH"/data/%s/iris/fpara.dat",planet->name);
    infile = fopen(string,"r");
    fscanf(infile,"%d %d",&nlat,&npress);
    fscanf(infile,"%s",header);
    /* 
     * Allocate memory 
     */
    lat          = fvector( 0,nlat-1,             dbmsname);
    press        = fvector( 0,npress-1,           dbmsname);
    fp_k         = fvector( 0,npress-1,           dbmsname);
    fp_k_y2      = fvector( 0,npress-1,           dbmsname);
    fp_data      = fvector( 0,nlat*npress-1,      dbmsname);
    fp_data_y2   = fvector( 0,nlat*npress-1,      dbmsname);
    buff_triplet = ftriplet(0,IMAX(nlat,npress)-1,dbmsname);

    for (k = 0; k < npress; k++) {

#if EPIC_PRECSION == DOUBLE_PRECISION
      fscanf(infile,"%lf",press+k);
#else
      fscanf(infile,"%lf",press+k);
#endif

      /* convert press to mks */
      press[k] *= 100.;
      /* spline on log p */
      press[k] = log(press[k]);
    }
    fscanf(infile,"%s",header);
    for (j = 0; j < nlat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(infile,"%lf",lat+j);
#else
      fscanf(infile,"%f",lat+j);
#endif

      /* convert lat to radians */
      lat[j] *= DEG;
      for (k = 0; k < npress; k++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile,"%lf",&(FP_DATA(k,j)));
#else
        fscanf(infile,"%f",&(FP_DATA(k,j)));
#endif

      }
    }
    fclose(infile);
    /* Compute spline information for each row: */
    /*
     * Pack float_triplet buffer. 
     */
    for (j = 0; j < nlat; j++) {
      buff_triplet[j].x = lat[j];
    }
    for (k = 0; k < npress; k++) {
      for (j = 0; j < nlat; j++) {
        buff_triplet[j].y = FP_DATA(k,j);
      }
      spline_pchip(nlat,buff_triplet);
      for (j = 0; j < nlat; j++) {
        FP_DATA_Y2(k,j) = buff_triplet[j].z;
      }
    }
    initialized = TRUE;
  }
  /* End initialization. */

  if (mode == 1) {
    /* Synthetic data */
    z  = -log(p/grid.press0);
    fp = a1+a2*z+a3*(1.-cos(4.*latr));
  }
  else if (mode == 0) {
    /* Data */
    if (latr < lat[0]) {
      sprintf(Message,"lat = %f below data range = %f",latr/DEG,lat[0]/DEG);
      epic_error(dbmsname,Message);
    }
    else if (latr > lat[nlat-1]) {
      sprintf(Message,"lat = %f above data range = %f",latr/DEG,lat[nlat-1]/DEG);
      epic_error(dbmsname,Message);
    }

    for (j = 0; j < nlat; j++) {
      buff_triplet[j].x = lat[j];
    }
    j = find_place_in_table(nlat,buff_triplet,latr,&lat_d);
    /* Make column at latr */
    for (k = 0; k < npress; k++) {
      /*
       * Pack float_triplet buffer segment.
       */
      buff_triplet[j  ].y = FP_DATA(   k,j  );
      buff_triplet[j+1].y = FP_DATA(   k,j+1);
      buff_triplet[j  ].z = FP_DATA_Y2(k,j  );
      buff_triplet[j+1].z = FP_DATA_Y2(k,j+1);

      fp_k[k] = splint_pchip(latr,buff_triplet+j,lat_d);
    }
    /*
     * Pack float_triplet buffer.
     */
    for (k = 0; k < npress; k++) {
      buff_triplet[k].x = press[k];
      buff_triplet[k].y = fp_k[k];
    }
    spline_pchip(npress,buff_triplet);

    log_p = log(p);
    if (log_p < press[0]) {
      log_p = press[0];
    }
    else if (log_p > press[npress-1]) {
      log_p = press[npress-1];
    }
    k  = find_place_in_table(npress,buff_triplet,log_p,&p_d);
    fp = splint_pchip(log_p,buff_triplet+k,p_d);
  }
  else {
    sprintf(Message,"unrecognized mode");
    epic_error(dbmsname,Message);
  }

  return MIN(fp,1.);

}

/*====================== end of fp_yp() ======================================*/

/*====================== p_gz() ==============================================*/
/*
 * Returns pressure at the bottom of the model given gz,
 * based on p vs z data file for specified planet.
 * For gas giants, just returns grid.pbot;
 */
EPIC_FLOAT p_gz(planetspec *planet,
                int         J,
                EPIC_FLOAT  gz)
{
  int
    i,
    ntp;
  static int
    ndat,
    initialized=FALSE;
  char
    header[N_STR],
    data_file[FILE_STR];
  EPIC_FLOAT
    p,
    temperature,
    tmp,
    gz_d,
    rttop;
  static EPIC_FLOAT
    pbot0,
    *pdat,
    *tdat,
    *zdat,
    *neglogpdat;
  static float_triplet
    *buff_triplet;
  FILE
    *infile;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="p_gz";

  if (!initialized) {
    if (strcmp(planet->type,"gas-giant") == 0) {
      pbot0 = grid.pbot;
    }
    else if (strcmp(planet->name,"earth") == 0) {
      pbot0 = 1013.25*100.;
    }
    else if (strcmp(planet->name,"held_suarez") == 0) {
      pbot0 = 1000.*100.;
    }
    else if (strcmp(planet->name,"titan") == 0) {
      /*
       * Set the reference surface pressure to be 1.5 bar.
       */
      pbot0 = 1500.*100;
    }
    else if (strncmp(planet->name,"venus",5) == 0) {
      /*
       * This case includes all systems with a name that starts with "venus".
       */
      sprintf(data_file,EPIC4_PATH"/data/venus/VIRA.dat");
      infile = fopen(data_file,"r");
      if (!infile) {
        sprintf(Message,"unable to read %s",data_file);
        epic_error(dbmsname,Message);
      }
      /* Skip header. */
      for (i = 0; i < 5; i++) {
        fgets(header,N_STR,infile);
      }
      /* Input number of data points. */
      fscanf(infile,"%d",&ndat);
      /*
       * Allocate memory.
       */
      pdat         = fvector( 0,ndat-1,dbmsname);
      zdat         = fvector( 0,ndat-1,dbmsname);
      neglogpdat   = fvector( 0,ndat-1,dbmsname);
      buff_triplet = ftriplet(0,ndat-1,dbmsname);
      /*
       * Input data.
       */
      for (i = 0; i < ndat; i++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile," %lf %lf %lf %lf",zdat+i,&temperature,pdat+i,&tmp);
#else
        fscanf(infile," %f %f %f %f",zdat+i,&temperature,pdat+i,&tmp);
#endif

        /* Convert km to m. */
        zdat[i] *= 1000.;
        /* Convert mbar to Pa. */
        pdat[i] *= 100.;
        /* Interpolate on -log(p). */
        neglogpdat[i] = -log(pdat[i]);
        /* Discard the rest of the line. */
        fgets(header,N_STR,infile);
      }
      fclose(infile);
    }
    else if (strncmp(planet->name,"mars",4) == 0) {
      /*
       * This case includes all systems with a name that starts with "mars".
       */
      sprintf(data_file,EPIC4_PATH"/data/mars/LMD.Hellas.dat");
      infile = fopen(data_file,"r");
      if (!infile) {
        sprintf(Message,"unable to read %s",data_file);
        epic_error(dbmsname,Message);
      }
      /* Skip header. */
      for (i = 0; i < 5; i++) {
        fgets(header,N_STR,infile);
      }
      /* Input number of data points. */
      fscanf(infile,"%d",&ndat);
      /*
       * Allocate memory.
       */
      pdat         = fvector( 0,ndat-1,dbmsname);
      zdat         = fvector( 0,ndat-1,dbmsname);
      neglogpdat   = fvector( 0,ndat-1,dbmsname);
      buff_triplet = ftriplet(0,ndat-1,dbmsname);
      /*
       * Input data.
       */
      for (i = 0; i < ndat; i++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile," %lf %lf %lf %lf",zdat+i,&temperature,pdat+i,&tmp);
#else
        fscanf(infile," %f %f %f %f",zdat+i,&temperature,pdat+i,&tmp);
#endif

        /* Convert km to m. */
        zdat[i] *= 1000.;
        /* Convert hPa to Pa. */
        pdat[i] *= 100.;
        /* Interpolate on -log(p). */
        neglogpdat[i] = -log(pdat[i]);
        /* Discard the rest of the line. */
        fgets(header,N_STR,infile);
      }
      fclose(infile);

    }
    else {
      sprintf(Message,"%s not yet implemented",planet->name);
      epic_error(dbmsname,Message);
    }

    initialized = TRUE;
  }
  /* End of initialization. */

  if (strcmp(planet->type,"gas-giant") == 0) {
    p = pbot0;
  }
  else if (strcmp(planet->name,"earth") == 0) {
    /*
     * Currently, GZ_SURFACE(J,I) = 0. everywhere.
     */
    if (gz == 0.) {
      p = pbot0;
    }
    else {
      sprintf(Message,"gz=%f; case gz != 0. not yet implemented for Earth\n",gz);
      epic_error(dbmsname,Message);
    }
  }
  else if (strcmp(planet->name,"held_suarez") == 0) {
    /*
     * GZ_SURFACE(J,I) = 0. everywhere.
     */
    if (gz == 0.) {
      p = pbot0;
    }
    else {
      sprintf(Message,"gz=%f; case gz != 0. not implemented for Held_Suarez\n",gz);
      epic_error(dbmsname,Message);
    }
  }
  else if (strcmp(planet->name,"titan") == 0) {
    if (gz == 0.) {
      p = pbot0;
    }
    else {
      sprintf(Message,"gz=%f; case gz != 0. not implemented for Titan\n",gz);
      epic_error(dbmsname,Message);
    }
  }
  else {
    if (gz < grid.g[2*J+1]*zdat[0]) {
      sprintf(Message,"gz=%g < g*zdat[0]=%g",gz,grid.g[2*J+1]*zdat[0]);
      epic_error(dbmsname,Message);
    }
    else if (gz >= grid.g[2*J+1]*zdat[ndat-1]) {
      /* Assume a constant scale height above data table. */
      rttop = grid.g[2*J+1]*(zdat[ndat-1]-zdat[ndat-2])/(neglogpdat[ndat-1]-neglogpdat[ndat-2]);
      p     = pdat[ndat-1]*exp(-(gz-grid.g[2*J+1]*zdat[ndat-1])/rttop);
    }
    else {
      /*
       * Use smooth, monotonic interpolation.
       */
      for (i = 0; i < ndat; i++) {
        buff_triplet[i].x = grid.g[2*J+1]*zdat[i];
        buff_triplet[i].y = neglogpdat[i];
      }
      spline_pchip(ndat,buff_triplet);
      i   = find_place_in_table(ndat,buff_triplet,gz,&gz_d);
      tmp = splint_pchip(gz,buff_triplet+i,gz_d);
      p   = exp(-tmp);
    }
  }

  /*
   * Screen for NaN.
   */
  if (!isfinite(p) || !isfinite(gz)) {
    sprintf(Message,"gz=%g, p=%g",gz,p);
    epic_error(dbmsname,Message);
  }

  return p;
}

/*====================== end of p_gz() ======================================*/

/*====================== rho_gz() ===========================================*/
/*
 * Returns density given gz, 
 * based on rho vs z data for specified planet.
 * For gas giants, returns value consistent with p_ref, theta_ref at K = grid.nk.
 */
EPIC_FLOAT rho_gz(planetspec *planet,
                  int         J,
                  EPIC_FLOAT  gz)
{
  int
    i,
    ntp;
  static int
    ndat,
    initialized=FALSE;
  char
    header[N_STR],
    data_file[FILE_STR];
  EPIC_FLOAT
    rho,
    fpara,mu,tmp,
    gz_d;
  static EPIC_FLOAT
    rttop,
    rho0,
    *rhodat,
    *neglogrhodat,
    *zdat;
  static float_triplet
    *buff_triplet;
  FILE
    *infile;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="rho_gz";

  if (!initialized) {
    if (strcmp(planet->type,"gas-giant") == 0) {
      sprintf(Message,"%s case not yet implemented",planet->type);
      epic_error(dbmsname,Message);
    }
    else {
      if (strncmp(planet->name,"venus",5) == 0) {
        /*
         * This case includes all systems with a name that starts with "venus".
         */
        sprintf(data_file,EPIC4_PATH"/data/venus/VIRA.dat");
        infile = fopen(data_file,"r");
        if (!infile) {
          sprintf(Message,"unable to read %s",data_file);
          epic_error(dbmsname,Message);
        }
        /* Skip header. */
        for (i = 0; i < 5; i++) {
          fgets(header,N_STR,infile);
        }
        /* Input number of data points. */
        fscanf(infile,"%d",&ndat);
        /*
         * Allocate memory.
         */
        rhodat       = fvector( 0,ndat-1,dbmsname);
        neglogrhodat = fvector( 0,ndat-1,dbmsname);
        zdat         = fvector( 0,ndat-1,dbmsname);
        buff_triplet = ftriplet(0,ndat-1,dbmsname);
        /*
         * Input data.
         */
        for (i = 0; i < ndat; i++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile," %lf %lf %lf %lf",zdat+i,&tmp,&tmp,rhodat+i);
#else
          fscanf(infile," %f %f %f %f",zdat+i,&tmp,&tmp,rhodat+i);
#endif

          /* Convert km to m. */
          zdat[i] *= 1000.;

          /* Interpolate on log density. */
          neglogrhodat[i] = -log(rhodat[i]);

          /* Discard the rest of the line. */
          fgets(header,N_STR,infile);
        }
        fclose(infile);
      }
      else if (strncmp(planet->name,"mars",4) == 0) {
        /*
         * This case includes all systems with a name that starts with "mars".
         */
        sprintf(data_file,EPIC4_PATH"/data/mars/LMD.Hellas.dat");
        infile = fopen(data_file,"r");
        if (!infile) {
          sprintf(Message,"unable to read %s",data_file);
          epic_error(dbmsname,Message);
        }
        /* Skip header. */
        for (i = 0; i < 5; i++) {
          fgets(header,N_STR,infile);
        }
        /* Input number of data points. */
        fscanf(infile,"%d",&ndat);
        /*
         * Allocate memory.
         */
        rhodat       = fvector( 0,ndat-1,dbmsname);
        neglogrhodat = fvector( 0,ndat-1,dbmsname);
        zdat         = fvector( 0,ndat-1,dbmsname);
        buff_triplet = ftriplet(0,ndat-1,dbmsname);
        /*
         * Input data.
         */
        for (i = 0; i < ndat; i++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile," %lf %lf %lf %lf",zdat+i,&tmp,&tmp,rhodat+i);
#else
          fscanf(infile," %f %f %f %f",zdat+i,&tmp,&tmp,rhodat+i);
#endif

          /* Convert km to m. */
          zdat[i] *= 1000.;

          /* Interpolate on log density. */
          neglogrhodat[i] = -log(rhodat[i]);

          /* Discard the rest of the line. */
          fgets(header,N_STR,infile);
        }
        fclose(infile);
      }
      else {
        sprintf(Message,"%s not yet implemented",planet->name);
        epic_error(dbmsname,Message);
      }
    }

    initialized = TRUE;
  }
  /* End of initialization. */

  if (strcmp(planet->type,"gas-giant") == 0) {
    rho = rho0;
  }
  else {
    if (gz < grid.g[2*J+1]*zdat[0]) {
      sprintf(Message,"gz=%g < g*zdat[0]=%g",gz,grid.g[2*J+1]*zdat[0]);
      epic_error(dbmsname,Message);
    }
    else if (gz >= grid.g[2*J+1]*zdat[ndat-1]) {
      /* Assume a constant scale height above data table. */
      rttop = grid.g[2*J+1]*(zdat[ndat-1]-zdat[ndat-2])/log(rhodat[ndat-2]/rhodat[ndat-1]);
      rho   = rhodat[ndat-1]*exp(-(gz-grid.g[2*J+1]*zdat[ndat-1])/rttop);
    }
    else {
      /*
       * Use smooth, monotonic interpolation.
       */
      for (i = 0; i < ndat; i++) {
        buff_triplet[i].x = grid.g[2*J+1]*zdat[i];
        buff_triplet[i].y = neglogrhodat[i];
      }
      spline_pchip(ndat,buff_triplet);
      i   = find_place_in_table(ndat,buff_triplet,gz,&gz_d);
      tmp = splint_pchip(gz,buff_triplet+i,gz_d);
      rho = exp(-tmp);
    }
  }

  /*
   * Screen for NaN.
   */
  if (!isfinite(rho) || !isfinite(gz)) {
    sprintf(Message,"gz=%g, rho=%g",gz,rho);
    epic_error(dbmsname,Message);
  }

  return rho;
}

/*====================== end of rho_gz() ====================================*/

/*====================== gz_p() =============================================*/
/*
 * Returns gz given pressure, based on p vs z data file for specified planet.
 * Used when setting up a hydrostatic initial condition.
 * For gas giants, just returns 0.
 */
EPIC_FLOAT gz_p(planetspec *planet,
                int         J,
                EPIC_FLOAT  p)
{
  int
    i,
    ntp;
  static int
    ndat,
    initialized=FALSE;
  char
    header[N_STR],
    data_file[FILE_STR];
  EPIC_FLOAT
    gz,neglogp,
    tmp,
    p_d,
    rttop;
  static EPIC_FLOAT
    gz0,
    *pdat,
    *tdat,
    *zdat,
    *neglogpdat;
  static float_triplet
    *buff_triplet;
  FILE
    *infile;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="gz_p";

  if (!initialized) {
    if (strcmp(planet->type,"gas-giant") == 0) {
      gz0 = 0.;
    }
    else if (strncmp(planet->name,"venus",5) == 0) {
      /*
       * This case includes all systems with a name that starts with "venus".
       */
      sprintf(data_file,EPIC4_PATH"/data/venus/VIRA.dat");
      infile = fopen(data_file,"r");
      if (!infile) {
        sprintf(Message,"unable to read %s",data_file);
        epic_error(dbmsname,Message);
      }
      /* Skip header. */
      for (i = 0; i < 5; i++) {
        fgets(header,N_STR,infile);
      }
      /* Input number of data points. */
      fscanf(infile,"%d",&ndat);
      /*
       * Allocate memory.
       */
      pdat         = fvector( 0,ndat-1,dbmsname);
      zdat         = fvector( 0,ndat-1,dbmsname);
      neglogpdat   = fvector( 0,ndat-1,dbmsname);
      buff_triplet = ftriplet(0,ndat-1,dbmsname);
      /*
       * Input data.
       */
      for (i = 0; i < ndat; i++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile," %lf %lf %lf %lf",zdat+i,&tmp,pdat+i,&tmp);
#else
        fscanf(infile," %f %f %f %f",zdat+i,&tmp,pdat+i,&tmp);
#endif

        /* Convert km to m. */
        zdat[i] *= 1000.;
        /* Convert p from mbar to Pa. */
        pdat[i] *= 100.;
        /* Interpolate on -log(p). */
        neglogpdat[i] = -log(pdat[i]);
        /* Discard the rest of the line. */
        fgets(header,N_STR,infile);
      }
      fclose(infile);
    }
    else if (strncmp(planet->name,"mars",4) == 0) {
      /*
       * This case includes all systems with a name that starts with "mars".
       */
      sprintf(data_file,EPIC4_PATH"/data/mars/LMD.Hellas.dat");
      infile = fopen(data_file,"r");
      if (!infile) {
        sprintf(Message,"unable to read %s",data_file);
        epic_error(dbmsname,Message);
      }
      /* Skip header. */
      for (i = 0; i < 5; i++) {
        fgets(header,N_STR,infile);
      }
      /* Input number of data points. */
      fscanf(infile,"%d",&ndat);
      /*
       * Allocate memory.
       */
      pdat         = fvector( 0,ndat-1,dbmsname);
      zdat         = fvector( 0,ndat-1,dbmsname);
      neglogpdat   = fvector( 0,ndat-1,dbmsname);
      buff_triplet = ftriplet(0,ndat-1,dbmsname);
      /*
       * Input data.
       */
      for (i = 0; i < ndat; i++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile," %lf %lf %lf %lf",zdat+i,&tmp,pdat+i,&tmp);
#else
        fscanf(infile," %f %f %f %f",zdat+i,&tmp,pdat+i,&tmp);
#endif

        /* Convert km to m. */
        zdat[i] *= 1000.;
        /* Convert p from mbar to Pa. */
        pdat[i] *= 100.;
        /* Interpolate on -log(p). */
        neglogpdat[i] = -log(pdat[i]);
        /* Discard the rest of the line. */
        fgets(header,N_STR,infile);
      }
      fclose(infile);
    }
    else {
      sprintf(Message,"not yet implemented for planet->name=%s",planet->name);
      epic_error(dbmsname,Message);
    }

    initialized = TRUE;
  }
  /* End of initialization. */

  if (strcmp(planet->type,"gas-giant") == 0) {
    gz = gz0;
  }
  else {
    if (p > pdat[0]) {
      sprintf(Message,"p=%g > pdat[0]=%g",p,pdat[0]);
      epic_error(dbmsname,Message);
    }
    else if (p <= pdat[ndat-1]) {
      /* Assume a constant scale height above data table. */
      rttop = grid.g[2*J+1]*(zdat[ndat-1]-zdat[ndat-2])/(neglogpdat[ndat-1]-neglogpdat[ndat-2]);
      gz    = grid.g[2*J+1]*zdat[ndat-1]+rttop*log(pdat[ndat-1]/p);
    }
    else {
      /*
       * Use smooth, monotonic interpolation.
       *
       * NOTE: This could be made more efficient, since multiple input p could be used
       *       with the same spline.  However, this is currently intended for initialization
       *       purposes only, so we have left it as is.
       */
      for (i = 0; i < ndat; i++) {
        buff_triplet[i].x = neglogpdat[i];
        buff_triplet[i].y = grid.g[2*J+1]*zdat[i];
      }
      spline_pchip(ndat,buff_triplet);
      neglogp = -log(p);
      i       = find_place_in_table(ndat,buff_triplet,neglogp,&p_d);
      gz      = splint_pchip(neglogp,buff_triplet+i,p_d);
    }
  }

  /*
   * Screen for NaN.
   */
  if (!isfinite(gz) || !isfinite(p)) {
    sprintf(Message,"gz=%g, p=%g",gz,p);
    epic_error(dbmsname,Message);
  }

  return gz;
}

/*====================== end of gz_p() ======================================*/

/*====================== set_u_spinup() =====================================*/

void set_u_spinup(planetspec       *planet,
                  init_defaultspec *def)
{
  register int
    K,J,I;
  register EPIC_FLOAT
    lat;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_u_spinup";

  if (def->u_spinup_scale == 0.) {
    memset(var.u_spinup.value,0,Nelem3d*sizeof(EPIC_FLOAT));
  }
  else if (def->u_spinup_scale == 1.) {
    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          U_SPINUP(K,J,I) = U(grid.it_uv,K,J,I);
        }
      }
    }
  }
  else {
    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          U_SPINUP(K,J,I) = def->u_spinup_scale*U(grid.it_uv,K,J,I);
        }
      }
    }
  }

  return;
}

/*====================== end of set_u_spinup() ==============================*/

/* * * * * * * * * * * * end of epic_funcs_init.c * * * * * * * * * * * * * * */









