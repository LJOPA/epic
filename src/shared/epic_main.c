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

/* * * * * * * * * * * * * * * * epic_main.c * * * * * * * * * * * * * * * * * 
 *                                                                           *
 *     Explicit Planetary Isentropic-Coordinate (EPIC) Model                 *
 *                                                                           *
 *     Model-development contributions from:                                 *
 *         E.  Charrette                                                     *
 *         E.  Colon                                                         *
 *         T.  Dowling, PI                                                   *
 *         K.  Emanuel                                                       *
 *         A.  Fischer                                                       *
 *         N.  Ghurtskaia                                                    *
 *         P.  Gierasch                                                      *
 *         T.  Greathouse                                                    *
 *         J.  Harrington                                                    *
 *         M.  Herman                                                        *
 *         A.  Herrnstein                                                    *
 *         R.  Lebeau                                                        *
 *         G.  Lee                                                           *
 *         J.  Matarese                                                      *
 *         R.  Morales-Juberias                                              *
 *         Cs. Palotai                                                       *
 *         D.  Raymond                                                       *
 *         C.  Santori                                                       *
 *         K.  Sayanagi                                                      *
 *         A.  Showman                                                       *
 *         M.  Sussman                                                       *
 *         E.  Thompson                                                      *
 *                                                                           *
 *     Flags of EPIC-model developers and users:                             *
 *         Basque Country                                                    *
 *         Belgium                                                           *
 *         Georgia                                                           *
 *         Hungary                                                           *
 *         Japan                                                             *
 *         Spain                                                             *
 *         United States                                                     * 
 *                                                                           *
 *     Important references related to EPIC model development:               *
 *                                                                           * 
 *     Arakawa and Lamb, 1981, A potential enstrophy and energy conserving   *
 *       scheme for the shallow water equations,  Monthly Weather Review     *
 *       109: 18-36.                                                         *
 *     Dennis JE & Schnabel RB, 1996, Numerical Methods for                  *
 *       Unconstrained Optimization and Nonlinear equations, Classics        *
 *       in Applied Mathematics, Vol 16, SIAM.                               *
 *     Dowling TE et al, 1998, The explicit planetary isentropic             *
 *        coordinate (EPIC) atmospheric model, Icarus 132, 221-238           *
 *     Dowling TE et al, 2006, The EPIC atmospheric model with an            *
 *        isentropic/terrain-following hybrid vertical coordinate,           *
 *        Icarus 182, 259-273                                                *
 *     Durran, 1991, The third-order Adams-Bashforth method: an attractive   *
 *       alternative to leapfrog time differencing, Monthly Weather Review   *
 *       119: 702-720                                                        *                   
 *     Hsu and Arakawa, 1990, Numerical modeling of the atmosphere with an   *
 *       isentropic vertical coordinate, Monthly Weather Review              *
 *       118: 1933-1959                                                      *
 *     Konor and Arakawa, 1997, Design of an atmospheric model based on a    *
 *       a generalized vertical coordinate, Monthly Weather Review           *
 *       125: 1649-1673                                                      *
 *     Press WH, Teukolsky SA, Vetterling WT, and Flannery BP,               *
 *       1992, Numerical Recipes in C, 2nd Ed., Cambridge.                   *
 *     Tannehill JC, Anderson DA, and Pletcher RH, 1997, Computational       *
 *       Fluid Mechanics and Heat Transfer, 2nd Ed., Taylor & Francis        *
 *     Wang, Y, 1995, An inverse balance equation in sigma coordinates for   *
 *       model initialization. Mon. Wea. Rev. 123, 482-488.                  *
 *                                                                           *
 *     The model has nk active vertical layers.  The vertical coordinate is  *
 *     a hybrid combination of potential temperature, theta, and a terrain   *
 *     following pressure coordinate, sigma.                                 * 
 *     The top and bottom layers are denoted by K = 1 and K = grid.nk,       *
 *     respectively. 
 *                                                                           *
 *     SI (mks) units are used unless otherwise noted.                       *
 *     Latitude is planetographic.                                           *
 *                                                                           *
 *     The input file epic.nc is generated by epic_initial.c.                *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*======================= main() ============================================*/

int main(int   argc,
         char *argv[]) 
{ 
  /*
   * NOTE: structures planet, grid, and var are declared globally in epic.h.
   */
  char  
    infile[80],            /*  input file                                */
    outfile[80],           /*  name of output data file, labeled by time */
    savdir[80],            /*  directory in which to save data           */
    sflag[80],             /*  string to hold command-line flag          */
    *ptr;
  int
    K,J,I,
    count,
    stop_request = FALSE,
    nc_err,nc_id;
  static int
    initialized = FALSE;
  unsigned int
    time_index         = 0;
  size_t
    extract_time_index = 0;
  static EPIC_FLOAT
    *Buff2D[NUM_WORKING_BUFFERS];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  static char
    dbmsname[]="epic_main";
  int
    idbms=0;

#ifdef EPIC_MPI
  MPI_Init(&argc,&argv);
  para.comm = MPI_COMM_WORLD;
  MPI_Comm_set_errhandler(para.comm,MPI_ERRORS_RETURN);
  MPI_Comm_rank(para.comm,&para.iamnode);
  MPI_Comm_size(para.comm,&para.nproc);
#endif

  /* 
   * Interpret command-line arguments: 
   */
  /* Start with defaults: */
  grid.itback            = INT_MAX;
  grid.itsave            = INT_MAX;
  grid.itextract         = INT_MAX;
  grid.itrun             = INT_MAX;
  grid.extract_append[0] = '\0';
  if (argc == 1) {
    /* Print help, exit: */
    /* Declare copyright: */
    declare_copyright();
    system("more "EPIC4_PATH"/help/epic_main.help");
    exit(1);
  }
  else {  
    /* Read flags: */
    for (count = 1; count < argc; count++) {
      sscanf(argv[count],"%s",sflag);
      if (strcmp(sflag,"-help") == 0 ||
          strcmp(sflag,"-h")    == 0) {
        /* Print help, exit: */
        /* Declare copyright: */
        declare_copyright();
        system("more "EPIC4_PATH"/help/epic_main.help");
        exit(1);
      }
      else if (strcmp(sflag,"-itback") == 0) {
        sscanf(argv[++count],"%d",&grid.itback);
      }
      else if (strcmp(sflag,"-itsave") == 0) {
        sscanf(argv[++count],"%d",&grid.itsave);
      }
      else if (strcmp(sflag,"-itextract") == 0) {
        sscanf(argv[++count],"%d",&grid.itextract);
      }
      else if (strcmp(sflag,"-append") == 0) {
        sscanf(argv[++count],"%s",grid.extract_append);
        /*
         * Open this extract file.
         */
        nc_err = nc_open(grid.extract_append,NC_NOWRITE,&nc_id);
        if (nc_err != NC_NOERR) {
          sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),grid.extract_append);
          epic_error(dbmsname,Message);
        }
        /*
         * Get the id for the time dimension, which is the unlimited dimension.
         */
        nc_err = nc_inq_unlimdim(nc_id,&var.info[0].coorid[NETCDF_T_INDEX]);
        if (nc_err != NC_NOERR) {
          sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),grid.extract_append);
          epic_error(dbmsname,Message);
        }
        /*
         * Read the number of time planes already stored.
         */
        nc_err = nc_inq_dimlen(nc_id,var.info[0].coorid[NETCDF_T_INDEX],&extract_time_index);
        if (nc_err != NC_NOERR) {
          sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),grid.extract_append);
          epic_error(dbmsname,Message);
        }
      }
      else if (strcmp(sflag,"-itrun") == 0) {
        sscanf(argv[++count],"%d",&grid.itrun);
      }
      else {
        if (count == argc-1) {
          /* The last command-line argument is the input filename: */
          sscanf(argv[argc-1],"%s",infile);
        }
        else {
          fprintf(stderr,"Unrecognized epic command-line flag: %s \n",sflag);
          exit(1);
        }
      }
    }
    if (grid.itextract < INT_MAX && grid.itback == INT_MAX) {
      grid.itback = grid.itextract;
    }
    if (grid.itrun < INT_MAX && grid.itsave == INT_MAX) {
      /*
       * If grid.itrun < 1, set grid.itsave to 1 to avoid dividing by zero.
       */
      grid.itsave = IMAX(1,grid.itrun);
    }
  }

  /* Allocate memory */
  if((planet=( planetspec *)malloc(sizeof( planetspec))) == 0) {
    sprintf(Message,"allocating space for planet");
    epic_error(dbmsname,Message);
  }

  /* 
   * Read model size, allocate memory for arrays,
   * read in rest of data:
   */

  /* NOTE: time_index not used for SIZE_DATA */
  var_read(planet,infile,SIZE_DATA,time_index);
  set_var_props();
  make_arrays(planet);

  for (I = 0; I < NUM_WORKING_BUFFERS; I++) {
    Buff2D[I] = fvector(0,Nelem2d-1,dbmsname);
  }

  /*
   * Time index 0 corresponds to IT_ZERO.
   */
  time_index = 0;
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
  /*
   * Set U_SPINUP.
   *
   * NOTE: We do not update U_SPINUP during a run, because
   *       that tends to cause positive feedbacks.
   */
  set_sponge();

  /*
   *  Setup thermodynamic functions.
   *
   *  NOTE: the return value cpr is the low-temperature limit of cp/rgas.
   */
  thermo_setup(planet,&planet->cpr);

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

  /*
   * Synchronize all the diagnostic variables with the prognostic variables
   * (i.e. prime the pump).
   */
  timestep(planet,SYNC_DIAGS_ONLY,Buff2D);

  if (grid.nmt_physics_on == 1) {
   /*
    * Initialize diagnostics for nmt_physics.
    */
    nmt_init_diag();
  }

  /* 
   * cd to infile directory so that saved files go there.
   */
  strcpy(savdir,infile);
  ptr = strrchr(savdir,'/');
  if (ptr != NULL) {
    *ptr = '\0';
    if(chdir(savdir) == -1){
      fprintf(stderr,"Warning: epic_main: Couldn't cd to save directory. \n");
      fprintf(stderr,"%s",savdir);
    }
  }

  if (IAMNODE == NODE0) {
    if (var.extract_on) {
      /*
       * Write EXTRACT_HEADER_DATA to extract file.
       */
      if (grid.extract_append[0] == '\0') {
        var_write(planet,"extract.nc",EXTRACT_HEADER_DATA,extract_time_index,0);
      }
    }
  }

  grid.itime = 0;
  while (grid.itime <= grid.itrun) {
    if (var.extract_on) {
      /* 
       * Write data to extract file periodically. 
       */
      if ((grid.itime)%(grid.itextract) == 0) {
        /*
         * The counter var.extract_time_index is incremented
         * after each use here.
         */
        if (grid.extract_append[0] != '\0') {
          if (grid.itime > 0) {
            var_write(planet,grid.extract_append,EXTRACT_DATA,extract_time_index++,0);
          }
        }
        else {
          var_write(planet,"extract.nc",EXTRACT_DATA,extract_time_index++,0);
        }
      }
    }

    /* Save data when requested. */
    if ( (grid.itime)%(grid.itsave) == 0 && grid.itime > 0) {
      strftime(outfile,80,"epic%Y_%m_%d_%H:%M:%S.nc",gmtime(&var.model_time));
      /*
       * Time index 0 corresponds to IT_ZERO.
       */
      time_index = 0;
      var_write(planet,outfile,ALL_DATA,time_index,0);

      /*
       * Check if early stop is requested.
       */
      if (IAMNODE == NODE0) {
        if (fopen("epic_stop","r")) {
          stop_request = TRUE;
        }
      }
#if defined(EPIC_MPI)
      MPI_Bcast(&stop_request,1,MPI_INT,NODE0,para.comm);
#endif
    } 

    /* Write out backup frame periodically; always have a good one. */
    if (  (grid.itime)%(grid.itback) == 0 && grid.itime > 0) {
      sprintf(outfile,"epic_back1.nc");
      if (IAMNODE == 0) { 
        sprintf(Message,"%s",outfile);
        if (0 != rename(Message,"epic_back2.nc") ) {
          if (errno != ENOENT){
            perror(outfile);
          } 
          else {
            errno = 0;
          }
        }
      }
      /*
       * Time index 0 corresponds to IT_ZERO.
       */
      time_index = 0;
      var_write(planet,outfile,ALL_DATA,time_index,0);

      /*
       * Check if early stop is requested.
       */
      if (IAMNODE == NODE0) {
        if (fopen("epic_stop","r")) {
          stop_request = TRUE;
        }
      }
#if defined(EPIC_MPI)
      MPI_Bcast(&stop_request,1,MPI_INT,NODE0,para.comm);
#endif
    }

    if (stop_request) {
      if (IAMNODE == NODE0) {
        fprintf(stdout,"Heeding stop request.\n");
      }
      break;
    }

    if (grid.itime < grid.itrun) {
      /* 
       * Take a step: 
       */
      timestep(planet,STEP_PROGS_AND_SYNC_DIAGS,Buff2D);
    }
    grid.itime++;
  }

  for (I = 0; I < NUM_WORKING_BUFFERS; I++) {
    free_fvector(Buff2D[I],0,Nelem2d-1,dbmsname);
  }

  free_arrays(planet); 
  free(planet);
  free_var_props();

#if defined(EPIC_MPI)
  MPI_Barrier(para.comm);
  MPI_Finalize();
#endif

  return 0;
}

/*======================= end of main()  ====================================*/

/* * * * * * * * * * * * end of epic_main.c  * * * * * * * * * * * * * * * * */
