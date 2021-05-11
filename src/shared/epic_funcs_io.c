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

/* * * * * * * * *  epic_funcs_io.c  * * * * * * * * * * * * * * * * 
 *                                                                 *
 *       Functions for EPIC model input and output.                *
 *       This file includes the following:                         *
 *                                                                 *
 *           setup_read_array()                                    *
 *           setup_write_array()                                   *
 *           read_array()                                          *
 *           write_array()                                         *
 *           var_read()                                            *
 *           var_write()                                           *
 *           number_spots_in_file()                                *
 *           read_spots_file()                                     *
 *           lookup_netcdf()                                       *
 *           define_netcdf()                                       *
 *           get_jlohi()                                           *
 *           prompt_extract_on()                                   *
 *           prompt_species_on()                                   *
 *           bcast_char(),bcast_int(),                             *
 *           bcast_float(),bcast_double()                          *
 *           read_spacing_file()                                   *
 *           read_t_vs_p()                                         *
 *           read_meridional_plane()                               *
 *           input_float()                                         *
 *           input_int(),input_string()                            *
 *           print_model_description()                             *
 *           print_zonal_info()                                    *
 *           print_vertical_column()                               *
 *           node0_barrier_print()                                 *
 *           scdswap()                                             *
 *           epic_error()                                          *
 *           epic_warning()                                        *
 *           declare_copyright()                                   *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*======================= setup_read_array() ================================*/
/*
 * A function that uses read_array() must first call this setup function.
 *
 * Returns the number of nodes on the machine running the model.
 */
int setup_read_array(void)
{
  int
    node,
    num_nodes,
   *jlo,
   *jhi;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="setup_read_array";

  num_nodes = grid.we_num_nodes;

  /* 
   * Allocate memory.
   */
  jlo = ivector(0,num_nodes-1,dbmsname);
  jhi = ivector(0,num_nodes-1,dbmsname);

  /*
   * Exchange jlohi information.
   */
  for (node = 0; node < num_nodes; node++) {
    if (node == IAMNODE) {
      jlo[node] = JLO;
      jhi[node] = JHI;
    }

#if defined(EPIC_MPI)
    MPI_Bcast(jlo+node,1,MPI_INT,node,para.comm);
    MPI_Bcast(jhi+node,1,MPI_INT,node,para.comm);
#endif

  }

  /*
   * Setup get_jlohi().
   */
  get_jlohi(SETUP_GET_JLOHI,num_nodes,jlo,jhi);

  /*
   * Free allocated memory.
   */
  free_ivector(jlo,0,num_nodes-1,dbmsname);
  free_ivector(jhi,0,num_nodes-1,dbmsname);

  return num_nodes;
}

/*======================= end of setup_read_array() =========================*/

/*======================= setup_write_array() ===============================*/
/*
 * A function that calls write_array() must first call this setup function.
 *
 * Returns the number of nodes on the machine running the model.
 */

int setup_write_array(void)
{
  int
    node,
    num_nodes;
  static int
    initialized=FALSE,
   *jlo,
   *jhi;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="setup_write_array";

  if (!initialized) {

    /* 
     * Allocate memory.
     */
    jlo = ivector(0,grid.we_num_nodes-1,dbmsname);
    jhi = ivector(0,grid.we_num_nodes-1,dbmsname);

    initialized = TRUE;
  }

  num_nodes = grid.we_num_nodes;

  /*
   * Exchange jlohi information.
   */
  for (node = 0; node < num_nodes; node++) {
    if (node == IAMNODE) {
      jlo[node] = JLO;
      jhi[node] = JHI;
    }

#if defined(EPIC_MPI)
    MPI_Bcast(jlo+node,1,MPI_INT,node,para.comm);
    MPI_Bcast(jhi+node,1,MPI_INT,node,para.comm);
#endif

  }

  /*
   * Setup get_jlohi().
   */
  get_jlohi(SETUP_GET_JLOHI,num_nodes,jlo,jhi);
 
  return num_nodes;
}

/*======================= end of setup_write_array() ========================*/

/*======================= read_array() ======================================*/

/*
 *   Basic steps:
 *
 *   1) NODE0 reads subarray data from the file into its buff_subarray. 
 *
 *   2) NODE0 sends data from its buff_subarray to the target node's buff_subarray.
 *
 *   3) The target node transfers data from its buff_subarray to the appropriate place.
 *
 * A function that uses read_array() must first call setup_read_array().
 * Call from all nodes. 
 *
 * For array types that have striped data (array_type != EPIC_FLOAT_ARRAY), pass
 * the name of the first float array that will be read in, and the naming
 * convention will be reconstructed. 
 *
 * NOTE: Does not apply boundary conditions.
 */

void read_array(int   node,
                int   dim,
                int  *start,
                int  *end,
                char *name,
                int   index,
                void *array,
                int   array_type,
                int   nc_id)
{
  register int
    n,
    K,J,I,
    nk,nj,ni,
    nlo,nhi,n_len,nkji_len,
    klo,khi,k_len,kji_len,
    jlo,jhi,j_len,ji_len,
    ilo,ihi,i_len,
    shift_buff,offset,
    i_bytes,
    al=FOURDIM-dim;
  int
    nc_varid,
    nc_err;
  char
    the_name[VAR_NM_SZ];
  size_t
    nc_start[FOURDIM], 
    nc_count[FOURDIM];
  EPIC_FLOAT
    *buff_subarray,
    *epic_float_array;

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
    dbmsname[]="read_array";

  if (IAMNODE != node  &&
      IAMNODE != NODE0) {
    /*
     * Return if not a participating node.
     */
    return;
  }

  /*
   * Cast pointers based on array_type.
   */
  nlo  = 0;
  if (array_type == EPIC_FLOAT_ARRAY) {
    epic_float_array = (EPIC_FLOAT *)array;
    nhi              = 1-1;
  }
  else {
    sprintf(Message,"unrecognized array_type=%d",array_type);
    epic_error(dbmsname,Message);
  }

  khi = klo = 0;
  jhi = jlo = 0;
  ihi = ilo = 0;

  nk  = grid.nk;
  nj  = grid.nj;
  ni  = grid.ni;

  if (dim >= TWODIM) {
    ilo = start[0];
    ihi = end[  0]; 

    jlo = start[1];
    jhi = end[  1];
   
    if (dim >= THREEDIM) {
      klo = start[2];
      khi = end[  2];
    }
  }
  else {
    sprintf(Message,"dim = %d not recognized",dim);
    epic_error(dbmsname,Message);
  }

  n_len      = nhi-nlo+1;
  k_len      = khi-klo+1;
  j_len      = jhi-jlo+1;
  i_len      = ihi-ilo+1;
  i_bytes    = i_len*sizeof(EPIC_FLOAT);
  ji_len     = j_len*i_len;
  kji_len    = k_len*ji_len;
  nkji_len   = n_len*kji_len;
  shift_buff = ilo+(jlo)*i_len+(klo)*ji_len+(nlo)*kji_len;

  /*
   * Allocate memory for buff_subarray:
   */
  buff_subarray = fvector(0,nkji_len-1,dbmsname);

 /*
  * 1) NODE0 reads subarray data and stores it in buff_subarray.
  */
  if (IAMNODE == NODE0) {
    nc_count[NETCDF_T_INDEX] = 1;
    nc_count[NETCDF_K_INDEX] = 1; 
    nc_count[NETCDF_J_INDEX] = 1;  
    nc_count[NETCDF_I_INDEX] = ihi-ilo+1; 

    nc_start[NETCDF_T_INDEX] = start[3];
    nc_start[NETCDF_I_INDEX] = 0;

    for (n = nlo; n <= nhi; n++) {
      /*
       * Get the netCDF variable ID.
       * First, reconstruct its name.
       */
      if (array_type == EPIC_FLOAT_ARRAY) {
        strcpy(the_name,name);
      }
      else {
        sprintf(Message,"unrecognized array_type=%d",array_type);
        epic_error(dbmsname,Message);
      }
      nc_err = nc_inq_varid(nc_id,the_name,&nc_varid);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }

      for (K = klo; K <= khi; K++) {
        nc_start[NETCDF_K_INDEX] = K-klo;
        for (J = jlo; J <= jhi; J++) {
          nc_start[NETCDF_J_INDEX] = J-grid.jlo;
 
#if EPIC_PRECISION == DOUBLE_PRECISION 
          nc_err = nc_get_vara_double(nc_id,nc_varid,nc_start+al,nc_count+al,
                                      &BUFF_SUBARRAY(n,K,J,ILO));
#else
          nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start+al,nc_count+al,
                                     &BUFF_SUBARRAY(n,K,J,ILO));
#endif

          if (nc_err != NC_NOERR) {
            sprintf(Message,"%s, K,J=%d,%d; %s",name,K,J,nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
        }
      }
    }
  }

 /*
  * 2) NODE0 sends data from its buff_subarray to the target node's buff_subarray.
  */

  if (node != NODE0) {
    /*
     * Send data from NODE0 to node.
     */

#if defined(EPIC_MPI)
    if (IAMNODE == NODE0) {
      MPI_Send(buff_subarray,nkji_len,float_type,node,index,para.comm);
    }
    else if (IAMNODE == node) {
      int
        count;
      MPI_Status
        status;

      MPI_Recv(buff_subarray,nkji_len,float_type,NODE0,index,para.comm,&status);

      /* 
       * Verify number of items received.
       */
      MPI_Get_count(&status,float_type,&count);
      if (count != nkji_len) {
        sprintf(Message,"count=%d != nkji_len=%d",count,nkji_len);
        epic_error(dbmsname,Message);
      }
    }
#endif

  }
  /*
   * Free allocated memory and return if not the target node.
   */
  if (IAMNODE != node) {
    free_fvector(buff_subarray,0,nkji_len-1,dbmsname);
    return;
  }

 /*
  * 3) The target node transfers data from its buff_subarray to the appropriate place.
  */
  for (K = klo; K <= khi; K++) {
    for (J = jlo; J <= jhi; J++) {
      if (array_type == EPIC_FLOAT_ARRAY) {
        /*
         * Use a fast string copy, since the data are contiguous.
         */
        offset = ilo+(J)*Iadim+(K)*Nelem2d-Shift3d;
        memcpy(epic_float_array+offset,&BUFF_SUBARRAY(0,K,J,ILO),i_bytes);
      }
      else {
        sprintf(Message,"unrecognized array_type=%d",array_type);
        epic_error(dbmsname,Message);
      }
    }
  }

  /*
   * Free allocated memory.
   */
  free_fvector(buff_subarray,0,nkji_len-1,dbmsname);

  return;
}

/*======================= end of read_array() ===============================*/

/*======================= write_array() =====================================*/

/*
 *   Basic steps:
 *
 *   1) Source node transfers data from the appropriate place to its buff_subarray.
 *
 *   2) Source node sends data from its buff_subarray to NODE0's buff_subarray.
 *
 *   3) NODE0 writes source node's data from its buff_subarray to the file.
 *
 * A function that uses write_array() must first call setup_write_array().
 * Call from all nodes. 
 *
 * For array types that have striped data (array_type != EPIC_FLOAT_ARRAY), pass
 * the name of the first float array that will be read in, and the naming
 * convention will be deduced from this.
 *
 */

void write_array(int    node,
                 int    dim,
                 int   *start,
                 int   *end,
                 int    stretch_ni,
                 char  *name,
                 int    index,
                 void  *array,
                 int    array_type,
                 int    nc_id)
{
  register int
    n,K,J,I,
    nk,nj,ni,
    nlo,nhi,n_len,nkji_len,
    klo,khi,k_len,kji_len,
    jlo,jhi,j_len,ji_len,
    ilo,ihi,i_len,
    shift_buff,offset,
    i_bytes,
    al=FOURDIM-dim;
  int
    nc_varid,
    nc_err;
  char
    the_name[VAR_NM_SZ];
  size_t
    nc_start[FOURDIM], 
    nc_count[FOURDIM];
  EPIC_FLOAT
    *buff_subarray,
    *buffer,
    *epic_float_array,
     tmp;

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
    dbmsname[]="write_array";

  /*
   * Return if not a participating node.
   */
  if (IAMNODE != node && 
      IAMNODE != NODE0) {
    return;
  }

  /*
   * Cast pointers based on array_type.
   */
  nlo  = 0;
  if (array_type == EPIC_FLOAT_ARRAY) {
    epic_float_array = (EPIC_FLOAT *)array;
    nhi              = 1-1;
  }
  else {
    sprintf(Message,"unrecognized array_type=%d",array_type);
    epic_error(dbmsname,Message);
  }

  khi = klo = 0;
  jhi = jlo = 0;
  ihi = ilo = 0;

  nk  = grid.nk;
  nj  = grid.nj;
  ni  = grid.ni;

  if (dim >= TWODIM) {
    ilo = start[0];
    ihi = end[  0]; 
    jlo = start[1];
    jhi = end[  1];
    if (dim >= THREEDIM) {
      klo = start[2];
      khi = end[  2];
    }
  }
  else {
    sprintf(Message,"dim = %d not recognized",dim);
    epic_error(dbmsname,Message);
  }

  n_len      = nhi-nlo+1;
  k_len      = khi-klo+1;
  j_len      = jhi-jlo+1;
  i_len      = ihi-ilo+1;
  i_bytes    = i_len*sizeof(EPIC_FLOAT);
  ji_len     = j_len*i_len;
  kji_len    = k_len*ji_len;
  nkji_len   = n_len*kji_len;
  shift_buff = ilo+(jlo)*i_len+(klo)*ji_len+(nlo)*kji_len;

  /*
   * Allocate memory.
   */
  buff_subarray = fvector(0,nkji_len-1,dbmsname);
  if (stretch_ni) {
    buffer = fvector(0,stretch_ni-ILO,dbmsname);
  }

  /*
   * 1) node obtains subarray data and stores it in buff_subarray.
   */
  if (IAMNODE == node) {
    for (K = klo; K <= khi; K++) {
      for (J = jlo; J <= jhi; J++) {
        if (array_type == EPIC_FLOAT_ARRAY) {
          /*
           * Use a fast string copy, since the data are contiguous.
           */
          offset = ilo+(J)*Iadim+(K)*Nelem2d-Shift3d;
          memcpy(&BUFF_SUBARRAY(0,K,J,ILO),epic_float_array+offset,i_bytes);
        }
        else {
          sprintf(Message,"unrecognized array_type=%d",array_type);
          epic_error(dbmsname,Message);
        }
      }
    }
  }

  /*
   * 2) node sends data from its buff_subarray to NODE0's buff_subarray.
   */

#if defined(EPIC_MPI)
  if (node != NODE0) {
    if (IAMNODE == node) {
      MPI_Send(buff_subarray,nkji_len,float_type,NODE0,index,para.comm);
    }
    else if (IAMNODE == NODE0) {
      int
        count;
      MPI_Status
        status;

      MPI_Recv(buff_subarray,nkji_len,float_type,node,index,para.comm,&status);

      /* 
       * Verify number of items received.
       */
      MPI_Get_count(&status,float_type,&count);
      if (count != nkji_len) {
        sprintf(Message,"count=%d != nkji_len=%d",count,nkji_len);
        epic_error(dbmsname,Message);
      }
    }
  }
#endif

  /*
   * 3) NODE0 writes data from its buff_subarray.
   */
  if (IAMNODE == NODE0) {
    nc_count[NETCDF_T_INDEX] = 1;
    nc_count[NETCDF_K_INDEX] = 1; 
    nc_count[NETCDF_J_INDEX] = 1;
    if (stretch_ni) {
      nc_count[NETCDF_I_INDEX] = stretch_ni-ilo+1;
    }
    else {
      nc_count[NETCDF_I_INDEX] = ihi-ilo+1;
    }

    nc_start[NETCDF_T_INDEX] = start[3];
    nc_start[NETCDF_I_INDEX] = 0;

    for (n = nlo; n <= nhi; n++) {
      /*
       * Get the netCDF variable ID.
       * First, reconstruct its name.
       */
      if (array_type == EPIC_FLOAT_ARRAY) {
        strcpy(the_name,name);
      }
      else {
        sprintf(Message,"unrecognized array_type=%d",array_type);
        epic_error(dbmsname,Message);
      }
      nc_err = nc_inq_varid(nc_id,the_name,&nc_varid);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s, %s",the_name,nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }

      for (K = klo; K <= khi; K++) {
        nc_start[NETCDF_K_INDEX] = K-klo;
        for (J = jlo; J <= jhi; J++) {
          nc_start[NETCDF_J_INDEX] = J-grid.jlo;

#if EPIC_PRECISION == DOUBLE_PRECISION
          if (stretch_ni) {
            /*
             * Assume zonal symmetry and copy the ILO value into an I buffer.
             */
            tmp = BUFF_SUBARRAY(n,K,J,ILO);
            for (I = ILO; I <= stretch_ni; I++) {
              buffer[I-ILO] = tmp;
            }
            nc_err = nc_put_vara_double(nc_id,nc_varid,nc_start+al,nc_count+al,buffer);
          }
          else {
            nc_err = nc_put_vara_double(nc_id,nc_varid,nc_start+al,nc_count+al,&BUFF_SUBARRAY(n,K,J,ILO));
          }
          if (nc_err != NC_NOERR) {
            sprintf(Message,"nc_put_vara_double(),%s",nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
#else
          if (stretch_ni) {
            /*
             * Assume zonal symmetry and copy the ILO value into an I buffer.
             */
            tmp = BUFF_SUBARRAY(n,K,J,ILO);
            for (I = ILO; I <= stretch_ni; I++) {
              buffer[I-ILO] = tmp;
            }
            nc_err = nc_put_vara_float(nc_id,nc_varid,nc_start+al,nc_count+al,buffer);

          }
          else {
            nc_err = nc_put_vara_float(nc_id,nc_varid,nc_start+al,nc_count+al,
                                       &BUFF_SUBARRAY(n,K,J,ILO));
          }
          if (nc_err != NC_NOERR) {
            sprintf(Message,"array=%s, nc_put_vara_float(), %s",
                            name,nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
#endif
        }
      }
    }
  }

  /*
   * Free allocated memory for buff_subarray:
   */
  free_fvector(buff_subarray,0,nkji_len-1,dbmsname);
  if (stretch_ni) {
    free_fvector(buffer,0,stretch_ni-ILO,dbmsname);
  }

  return;
}

/*======================= end of write_array() ==============================*/

/*======================= var_read() ========================================*/

/*
 * Read in the variables and applies boundary conditions.
 *
 * NOTE: When adding new variables, remember to apply boundary conditions
 *       here.
 */

void var_read(planetspec   *planet,
              char         *infile,
              int           portion,
              unsigned int  time_index)
{
  int 
    K,J,I,
    it,itlo,ithi,
    nk,
    node,
    is,iq,
    i;
  int
    start[FOURDIM],
    end[FOURDIM],
    num_nodes;
  static char
    **gattname=NULL,
    **varname =NULL;
  static int
    ngatts    =0,
    num_progs =0;
  int
    nc_err,nc_id;
  nc_type
    the_nc_type;     /* NOTE: Used in i/o macros. */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="var_read";

  /*
   * NOTE: Call lookup_netcdf() from all nodes so that
   *       its calls to MPI_Bcast() will work.
   */
  nc_err = lookup_netcdf(infile,&nc_id,&ngatts,&gattname,&num_progs,&varname);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),infile);
    epic_error(dbmsname,Message);
  }

  /*
   * Open NODE0's data connection.
   */
  if (IAMNODE == NODE0) {
    if (portion == SIZE_DATA) {
      fprintf(stdout,"Reading SIZE_DATA from %s \n",infile);
    }
    else if (portion == POST_SIZE_DATA) {
      fprintf(stdout,"Reading POST_SIZE_DATA from %s \n",infile);
    }
    else if (portion == HEADER_DATA) {
      fprintf(stdout,"Reading HEADER_DATA from %s \n",infile);
    }
    else if (portion == EXTRACT_HEADER_DATA) {
      fprintf(stdout,"Reading EXTRACT_HEADER_DATA from %s \n",infile);
    }
    else if (portion == POST_HEADER_DATA) {
      fprintf(stdout,"Reading POST_HEADER_DATA from %s \n",infile);
    }
    else if (portion == EXTRACT_DATA) {
      fprintf(stdout,"Reading EXTRACT_DATA from %s \n",infile);
    }
    else if (portion == ALL_DATA) {
      fprintf(stdout,"Reading ALL_DATA from %s \n",infile);
    }
    else {
      sprintf(Message,"unrecognized portion = %d",portion);
      epic_error(dbmsname,Message);
    }
    fflush(stderr);
  }

  /*
   * Read in size of model and other data needed by make_arrays().
   */
  if (portion == SIZE_DATA           || 
      portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA ||
      portion == ALL_DATA              ) {
    READF(&grid.data_version,grid_data_version,1);

    READTIME(&var.start_time,var_start_time);

    READI(&planet->index,planet_index,1);
    READC(planet->name,planet_name,32);
    READC(planet->type,planet_type,16);
    READC(planet->orbital_epoch,planet_orbital_epoch,8);
    READF(&planet->re,planet_re,1);
    READF(&planet->rp,planet_rp,1);
    READF(&planet->omega_sidereal,planet_omega_sidereal,1);
    READF(&planet->omega_synodic,planet_omega_synodic,1);
    READF(&planet->cp,planet_cp,1);
    READF(&planet->rgas,planet_rgas,1);

    planet->cpr   = planet->cp/planet->rgas;
    planet->kappa = 1./planet->cpr;

    READF(&planet->GM,planet_GM,1);
    READF(&planet->J2,planet_J2,1);
    READF(&planet->x_he,planet_x_he,1);
    READF(&planet->x_h2,planet_x_h2,1);
    READF(&planet->x_3,planet_x_3,1);
    READF(&planet->a,planet_a,1);
    READF(&planet->e,planet_e,1);
    READF(&planet->i,planet_i,1);
    READF(&planet->lon_ascending_node,planet_lon_ascending_node,1);
    READF(&planet->lon_perihelion,planet_lon_perihelion,1);
    READF(&planet->mean_lon,planet_mean_lon,1);
    READF(&planet->orbit_period,planet_orbit_period,1);
    READF(&planet->vernal_equinox_anomaly,planet_vernal_equinox_anomaly,1);
    READF(&planet->kinvisc,planet_kinvisc,1);
    READF(&planet->dynvisc,planet_dynvisc,1);
    READF(&planet->k_a,planet_k_a,1);

    READI(&grid.nk,grid_nk,1);
    READI(&grid.nj,grid_nj,1);
    READI(&grid.ni,grid_ni,1);

    READI(&grid.jtp,grid_jtp,1);

    READC(grid.geometry,grid_geometry,GEOM_STR);
    READC(grid.uv_timestep_scheme,grid_uv_timestep_scheme,N_STR);
    READC(grid.turbulence_scheme,grid_turbulence_scheme,N_STR);
    READC(grid.stability_factor,grid_stability_factor,N_STR);

    READC(var.hdry.advection_scheme,var_hdry_advection_scheme,N_STR);
    READC(var.theta.advection_scheme,var_theta_advection_scheme,N_STR);
    READC(var.fpara.advection_scheme,var_fpara_advection_scheme,N_STR);
    READC(var.nu_turb.advection_scheme,var_nu_turb_advection_scheme,N_STR);
    READC(var.species[FIRST_SPECIES].advection_scheme,var_species_advection_scheme,N_STR);
    for (is = FIRST_SPECIES+1; is <= LAST_SPECIES; is++) {
      sprintf(var.species[is].advection_scheme,"%s",var.species[FIRST_SPECIES].advection_scheme);
    }

    READF(&grid.globe_lonbot,grid_globe_lonbot,1);
    READF(&grid.globe_lontop,grid_globe_lontop,1);
    READF(&grid.globe_latbot,grid_globe_latbot,1);
    READF(&grid.globe_lattop,grid_globe_lattop,1);
    READC(grid.f_plane_map,grid_f_plane_map,GEOM_STR);
    READF(&grid.f_plane_lat0,grid_f_plane_lat0,1);
    READF(&grid.f_plane_half_width,grid_f_plane_half_width,1);

    READI(grid.wrap,grid_wrap,TOPDIM);
    READI(grid.pad,grid_pad,TOPDIM);
    READI(&grid.jlo,grid_jlo,1);
    READI(&grid.jfirst,grid_jfirst,1);
    READI(&grid.jlast,grid_jlast,1);
    READI(&grid.k_sponge,grid_k_sponge,1);
    READI(&grid.k_sigma,grid_k_sigma,1);

    READF(&grid.du_vert,grid_du_vert,1);

    READI(&grid.newt_cool_on,grid_newt_cool_on,1);
    READI(&grid.microphysics_on,grid_microphysics_on,1);
    READI(&grid.nmt_physics_on,grid_nmt_physics_on,1);
    READI(var.on_list,var_on_list,LAST_INDEX-FIRST_INDEX+1);
    READI(&var.extract_on,var_extract_on,1);
    READI(var.extract_on_list,var_extract_on_list,LAST_INDEX-FIRST_INDEX+1);

    READI(&var.ntp,var_ntp,1);
    READI(&var.n_t_cool,var_n_t_cool,1);

    if (grid.nmt_physics_on == 1) {
      READI(&nmt.nz,nmt_nz,1);
      READI(&nmt.land,nmt_land,1);
      READI(&nmt.dorad,nmt_dorad,1);
      READI(&nmt.frad,nmt_frad,1);
      READF(&nmt.radcool,nmt_radcool,1);
      READF(&nmt.tpause,nmt_tpause,1);
      READF(&nmt.radbrk,nmt_radbrk,1);
      READF(&nmt.cvc,nmt_cvc,1);
      READF(&nmt.cvs,nmt_cvs,1);
      READF(&nmt.cvp,nmt_cvp,1);
      READF(&nmt.cve,nmt_cve,1);
      READF(&nmt.theslop,nmt_theslop,1);
      READF(&nmt.pstiff,nmt_pstiff,1);
      READF(&nmt.pscale,nmt_pscale,1);
      READF(&nmt.pbltop,nmt_pbltop,1);
      READF(&nmt.cdrag,nmt_cdrag,1);
      READF(&nmt.wscale,nmt_wscale,1);
      READF(&nmt.sst,nmt_sst,1);
      READF(&nmt.lfrac,nmt_lfrac,1);
      READF(&nmt.eflux0,nmt_eflux0,1);
      READF(&nmt.cld,nmt_cld,1);
      READF(&nmt.cfract,nmt_cfract,1);
      READF(&nmt.dt,nmt_dt,1);
    }
  }

  if (portion == SIZE_DATA) {
    if (IAMNODE == NODE0) {
      nc_close(nc_id);
    }
    return;
  }

  /********************************************** 
   *                                            *
   * Seam between SIZE_DATA and POST_SIZE_DATA. *
   *                                            * 
   **********************************************/

  if (portion == POST_SIZE_DATA      ||
      portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA ||
      portion == ALL_DATA              ) {
    EPIC_FLOAT
      *buffer;

    nk = grid.nk;

    READF(&grid.dlt,grid_dlt,1);
    READF(&grid.dln,grid_dln,1);
    READI(&grid.dt,grid_dt,1);

    READF(&var.time_fp_bar,var_time_fp_bar,1);

    if (var.ntp > 0) {
      /*
       * Read the reference temperature sounding profile.
       */
      READF(var.pdat,var_pdat,var.ntp);
      READF(var.tdat,var_tdat,var.ntp);
      READF(var.dtdat,var_dtdat,var.ntp);
      /*
       * Convert input pressure from hPa to Pa.
       */
      for (i = 0; i < var.ntp; i++) {
        var.pdat[i] *= 100.;
      }
    }

    if (var.n_t_cool > 0) {
      /*
       * Read Newtonian-cooling time constant profile.
       * Since float_triplet memory is striped, use a buffer.
       */
      buffer = fvector(0,var.n_t_cool-1,dbmsname);
      READF(buffer,var_t_cool_p,var.n_t_cool);
      for (i = 0; i < var.n_t_cool; i++) {
        var.t_cool_table[i].x = -log(buffer[i]*100.);
      }
      READF(buffer,var_t_cool_tau,var.n_t_cool);
      for (i = 0; i < var.n_t_cool; i++) {
        var.t_cool_table[i].y = buffer[i];
      }
      free_fvector(buffer,0,var.n_t_cool-1,dbmsname);
    }

    READD(grid.sigmatheta,grid_sigmatheta,2*(nk+1)+1);

    READF(grid.p_ref,grid_p_ref,2*(nk+1)+1);
    READF(grid.theta_ref,grid_theta_ref,2*(nk+1)+1);
    READF(grid.h_min,grid_h_min,nk+1);

    READF(&grid.sgth_bot,grid_sgth_bot,1);
    READF(&grid.sgth_top,grid_sgth_top,1);
    /*
     * Double precision to improve calculation of diagnostic theta.
     */
    READD(&grid.zeta0,grid_zeta0,1);
    READD(&grid.zeta1,grid_zeta1,1);
    READD(&grid.hybrid_alpha,grid_hybrid_alpha,1);
    READD(&grid.sigma_sigma,grid_sigma_sigma,1);

    READI(&grid.newt_cool_adjust,grid_newt_cool_adjust,1);
    READI(&grid.diffusion_direction,grid_diffusion_direction,1);
    READC(grid.eos,grid_eos,8);
    READF(&grid.prandtl,grid_prandtl,1);
    READF(&grid.tau_drag,grid_tau_drag,1);
    READI(&grid.drag_v,grid_drag_v,1);
    READI(&grid.drag_zonal_avg,grid_drag_zonal_avg,1);

    READF(&grid.pbot,grid_pbot,1);
    READF(&grid.press0,grid_press0,1);
    READF(&grid.mont0,grid_mont0,1);
    READF(&grid.topo_scale,grid_topo_scale,1);

    READF(&grid.nudiv_nondim,grid_nudiv_nondim,1);
    READF(grid.nu_nondim,grid_nu_nondim,MAX_NU_ORDER+1);
    READD(grid.nu,grid_nu,MAX_NU_ORDER+1);
  }

  if (portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA   ) {
    if (IAMNODE == NODE0) {
      nc_close(nc_id);
    }
    return;
  }

  /********************************************************************* 
   *                                                                   *
   * Seam between HEADER_DATA/EXTRAC_HEADER_DATA and POST_HEADER_DATA. *
   *                                                                   * 
   *********************************************************************/

  /*
   * Read non-array data that is time dependent.
   */
  READTIME(&var.model_time,var_model_time);

  READI(&grid.cfl_dt,grid_cfl_dt,1);

  READI(&grid.aux_a,grid_aux_a,1);
  READI(&grid.aux_b,grid_aux_b,1);
  READI(&grid.aux_c,grid_aux_c,1);
  READF(&grid.aux_fa,grid_aux_fa,1);
  READF(&grid.aux_fb,grid_aux_fb,1);
  READF(&grid.aux_fc,grid_aux_fc,1);

  if (IAMNODE == 0) {
    fprintf(stdout,"  0%%");
    fflush(stdout);
  }

  /*
   * Setup read_array().
   */
  num_nodes = setup_read_array();

  /*
   * Load start and end vectors.
   */
  start[0]  = 1;
  end[  0]  = grid.ni;
  /*
   * Read 0 to grid.nk, such that the top and bottom interfaces are read in
   * for arrays defined on the interfaces (which is the majority).
   */
  start[2]  = 0;
  end[  2]  = grid.nk;
  start[3]  = time_index;
  end[  3]  = time_index;

  for (node = 0; node < num_nodes; node++) {
    if (IAMNODE == NODE0) {
      fprintf(stdout,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)node/num_nodes));
      fflush(stdout);
    }

    get_jlohi(node,num_nodes,start+1,end+1);

    /*
     * Read parameter arrays.
     */
    if (strcmp(planet->type,"terrestrial") == 0) {
      if (var.gz_surface.on) {
        /*
         * Read surface geopotential.
         */
        read_array(node,TWODIM,start,end,var.gz_surface.info[0].name,
                   var.gz_surface.info[0].index,var.gz_surface.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.u_spinup.on) {
      /*
       * Read spinup zonal-wind profile.  
       */
      read_array(node,FOURDIM,start,end,var.u_spinup.info[0].name,
                 var.u_spinup.info[0].index,var.u_spinup.value,EPIC_FLOAT_ARRAY,nc_id);
    }

    if (portion == EXTRACT_DATA) {
      sprintf(Message,"portion == EXTRACT_DATA not yet implemented");
      epic_error(dbmsname,Message);
    }

    /* 
     * Read variables, and their associated tendency and moment data as appropriate.
     */
    if (var.u.on) {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        if (portion != EXTRACT_DATA) {
          read_array(node,FOURDIM,start,end,var.u.info[0].name,
                     var.u.info[0].index,var.u.value,EPIC_FLOAT_ARRAY,nc_id);
        }
        else {
          sprintf(Message,"unrecognized portion=%d",portion);
          epic_error(dbmsname,Message);
        }
        if (portion != VAR_DATA) {
          read_array(node,THREEDIM,start,end,var.u.info_tend[0].name,
                     var.u.info_tend[0].index,var.u.tendency+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          read_array(node,THREEDIM,start,end,var.u.info_tend[1].name,
                     var.u.info_tend[1].index,var.u.tendency+IT_MINUS2*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        if (portion != EXTRACT_DATA) {
          read_array(node,FOURDIM,start,end,var.u.info[0].name,
                     var.u.info[0].index,var.u.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          read_array(node,FOURDIM,start,end,var.u.info[1].name,
                     var.u.info[1].index,var.u.value+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
        else {
          sprintf(Message,"unrecognized portion=%d",portion);
          epic_error(dbmsname,Message);
        }
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }

    if (var.v.on) {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        if (portion != EXTRACT_DATA) {
          read_array(node,FOURDIM,start,end,var.v.info[0].name,
                     var.v.info[0].index,var.v.value,EPIC_FLOAT_ARRAY,nc_id);
        }
        else {
          sprintf(Message,"unrecognized portion=%d",portion);
          epic_error(dbmsname,Message);
        }
        if (portion != VAR_DATA) {
          read_array(node,THREEDIM,start,end,var.v.info_tend[0].name,
                     var.v.info_tend[0].index,var.v.tendency+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          read_array(node,THREEDIM,start,end,var.v.info_tend[1].name,
                     var.v.info_tend[1].index,var.v.tendency+IT_MINUS2*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        if (portion != EXTRACT_DATA) {
          read_array(node,FOURDIM,start,end,var.v.info[0].name,
                     var.v.info[0].index,var.v.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          read_array(node,FOURDIM,start,end,var.v.info[1].name,
                     var.v.info[1].index,var.v.value+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
        else {
          sprintf(Message,"unrecognized portion=%d",portion);
          epic_error(dbmsname,Message);
        }
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }

    if (var.p3.on) {
      if ((portion == EXTRACT_DATA && var.p3.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        read_array(node,FOURDIM,start,end,var.p3.info[0].name,
                   var.p3.info[0].index,var.p3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.theta.on) {
      if ((portion == EXTRACT_DATA && var.theta.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        read_array(node,FOURDIM,start,end,var.theta.info[0].name,
                   var.theta.info[0].index,var.theta.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.fpara.on) {
      if ((portion == EXTRACT_DATA && var.fpara.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        read_array(node,FOURDIM,start,end,var.fpara.info[0].name,
                   var.fpara.info[0].index,var.fpara.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    for (iq = 0; iq < grid.nq; iq++) {
      if ((portion == EXTRACT_DATA && var.species[grid.is[iq]].phase[grid.ip[iq]].extract_on) ||
          (portion != EXTRACT_DATA                                        )   ) {
        read_array(node,FOURDIM,start,end,var.species[grid.is[iq]].phase[grid.ip[iq]].info[0].name,
                   var.species[grid.is[iq]].phase[grid.ip[iq]].info[0].index,
                   var.species[grid.is[iq]].phase[grid.ip[iq]].q,
                   EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.nu_turb.on) {
      if ((portion == EXTRACT_DATA && var.nu_turb.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        read_array(node,FOURDIM,start,end,var.nu_turb.info[0].name,
                   var.nu_turb.info[0].index,var.nu_turb.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }
  }

  /*
   * Close NODE0's data connection.
   */
  if (IAMNODE == NODE0) {
    nc_close(nc_id);
    fprintf(stdout,"\b\b\b\b%3d%%\n",100);
    fflush(stdout);
  }

  /*
   * Set solar longitude, L_s [deg], which is a function of time.
   */
  L_s = solar_longitude(planet,var.model_time);

  /*
   * Apply lateral boundary conditions.
   *
   * NOTE: bc_lateral() must be called from all nodes.
   */

  if (var.u_spinup.on) {
    /*
     * Assume u_spinup above model top does not change with height.
     */
    K = 0;
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        U_SPINUP(K,J,I) = U_SPINUP(K+1,J,I);
      }
    }
    bc_lateral(var.u_spinup.value,THREEDIM);
  }

  if (var.u.on) {
    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      bc_lateral(var.u.value+0*Nelem3d,THREEDIM);
      if (portion != VAR_DATA) {
        bc_lateral(var.u.tendency+IT_MINUS2*Nelem3d,THREEDIM);
        bc_lateral(var.u.tendency+IT_MINUS1*Nelem3d,THREEDIM);
      }
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      bc_lateral(var.u.value+IT_MINUS1*Nelem3d,THREEDIM);
      bc_lateral(var.u.value+IT_ZERO*Nelem3d,  THREEDIM);
    }
    else {
      sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
      epic_error(dbmsname,Message);
    }
  }

  if (var.v.on) {
    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      bc_lateral(var.v.value+0*Nelem3d,THREEDIM);
      if (portion != VAR_DATA) {
        bc_lateral(var.v.tendency+IT_MINUS2*Nelem3d,THREEDIM);
        bc_lateral(var.v.tendency+IT_MINUS1*Nelem3d,THREEDIM);
      }
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      bc_lateral(var.v.value+IT_MINUS1*Nelem3d,THREEDIM);
      bc_lateral(var.v.value+IT_ZERO*Nelem3d,  THREEDIM);
    }
    else {
      sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
      epic_error(dbmsname,Message);
    }
  }

  if (var.p3.on) {
    bc_lateral(var.p3.value,THREEDIM);
  }

  if (var.theta.on) {
    bc_lateral(var.theta.value,THREEDIM);
  }

  if (var.fpara.on) {
    bc_lateral(var.fpara.value,THREEDIM);
  }

  for (iq = 0; iq < grid.nq; iq++) {
    bc_lateral(var.species[grid.is[iq]].phase[grid.ip[iq]].q,THREEDIM);
  }

  if (var.nu_turb.on) {
    bc_lateral(var.nu_turb.value,THREEDIM);
  }

  if (strcmp(planet->type,"terrestrial") == 0) {
    if (var.gz_surface.on) {
      bc_lateral(var.gz_surface.value,TWODIM);
    }
  }

  /*
   * Make sure v = 0 at poles and channel walls.
   */
  if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
    itlo = 0;
    ithi = 0;
  }
  else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
    itlo = 0;
    ithi = 1;
  }
  else {
    sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
    epic_error(dbmsname,Message);
  }
  if (!grid.wrap[1]) {
    J = JLO;
    if (J == grid.jlo) {
      for (it = itlo; it <= ithi; it++) {
        for (K = KLOPAD; K <= KHIPAD; K++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            if (V(it,K,J,I) != 0.) {
              V(it,K,J,I) = 0.;
            }
          }
        }
      }
    }
    J = JHI+1;
    if (J == grid.nj+1) {
      for (it = itlo; it <= ithi; it++) {
        for (K = KLOPAD; K <= KHIPAD; K++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            if (V(it,K,J,I) != 0.) {
              V(it,K,J,I) = 0.;
            }
          }
        }
      }
    }
  }

  return;
}

/*======================= end of var_read() =================================*/

/*======================= var_write() =======================================*/

/*
 * Write the variables u, v, p, theta, etc.
 *
 * If ni == 1 and stretch_ni > 1, stretch the model from 2D to 3D on output.
 */

void var_write(planetspec   *planet,
               char         *outfile,
               int           portion,
               unsigned int  time_index,
               int           stretch_ni)
{
  int 
    K,J,I,kk,
    nk,
    node,
    num_nodes,
    iq,
    i;
  int
    start[FOURDIM],
    end[FOURDIM],
    nc_err,
    nc_id;
  static int
    initialized = FALSE;
  size_t
    t_index[1];
  EPIC_FLOAT
    the_time[1],
    the_L_s[1],
    avg;
  static EPIC_FLOAT
    *buff3d;
  nc_type
    the_nc_type;     /* NOTE: Used in i/o macros. */
  char 
    history[N_STR];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="var_write";

  if (!initialized) {
    /* Allocate memory. */
    buff3d = fvector(0,Nelem3d-1,dbmsname);

    initialized = TRUE;
  }

  /*
   * Open NODE0's data connection.
   */
  if (IAMNODE == NODE0) {
    if (portion == SIZE_DATA           ||
        portion == HEADER_DATA         ||
        portion == EXTRACT_HEADER_DATA ||
        portion == ALL_DATA              ) {
      /*
       * Create netCDF file and define variables and attributes.
       */
      handle_file_compression(outfile);
      nc_err = nc_create(outfile,NC_CLOBBER,&nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),outfile);
        epic_error(dbmsname,Message);
      }

      if (stretch_ni) {
        /*
         * Modify grid.dln to be consistent with stretch_ni,
         * before call to define_netcdf().
         */
        if (strcmp(grid.geometry,"globe") == 0) {
          grid.dln = (grid.globe_lontop-grid.globe_lonbot)/stretch_ni;
        }
        else if (strcmp(grid.geometry,"f-plane") == 0) {
          if (strcmp(grid.f_plane_map,"cartesian") == 0) {
            sprintf(Message,"-stretch_ni for case %s %s not yet implemented",grid.geometry,grid.f_plane_map);
            epic_error(dbmsname,Message);
          }
          else if (strcmp(grid.f_plane_map,"polar") == 0) {
            grid.dln = 360./(stretch_ni);
          }
        }
      }

      define_netcdf(planet,outfile,portion,stretch_ni,nc_id);
    }
    else {
      /*
       * Open existing netCDF file.
       */
      handle_file_compression(outfile);
      nc_err = nc_open(outfile,NC_WRITE,&nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),outfile);
        epic_error(dbmsname,Message);
      }

      /*
       * Put into define mode.
       */
      nc_err = nc_redef(nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"nc_redef(), %s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
    }

    if (portion == SIZE_DATA) {
      fprintf(stdout,"Writing SIZE_DATA to %s \n",outfile);
    }
    else if (portion == POST_SIZE_DATA) {
      fprintf(stdout,"Writing POST_SIZE_DATA to %s \n",outfile);
    }
    else if (portion == HEADER_DATA) {
      fprintf(stdout,"Writing HEADER_DATA to %s \n",outfile);
    }
    else if (portion == EXTRACT_HEADER_DATA) {
      fprintf(stdout,"Writing EXTRACT_HEADER_DATA to %s \n",outfile);
    }
    else if (portion == POST_HEADER_DATA) {
      fprintf(stdout,"Writing POST_HEADER_DATA to %s \n",outfile);
    }
    else if (portion == VAR_DATA) {
      fprintf(stdout,"Writing VAR_DATA to %s for timestep %lu\n",outfile,grid.itime);
    }
    else if (portion == EXTRACT_DATA) {
      fprintf(stdout,"Writing EXTRACT_DATA to %s for timestep %lu\n",outfile,grid.itime);
    }
    else if (portion == ALL_DATA) {
      fprintf(stdout,"Writing ALL_DATA to %s for timestep %lu\n",outfile,grid.itime);
    }
    else {
      sprintf(Message,"unrecognized portion = %d",portion);
      epic_error(dbmsname,Message);
    }
    fflush(stderr);
  }

  if (portion == SIZE_DATA           || 
      portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA ||
      portion == ALL_DATA              ) {

    /* 
     * Write global attributes.
     * We are following the Climate and Forecast (CF) Metadata Convention.
     */
    sprintf(Message,"CF-1.0");
    nc_put_att_text(nc_id,NC_GLOBAL,"Conventions",strlen(Message)+1,Message);
    if (portion == EXTRACT_HEADER_DATA) {
      sprintf(Message,"Explicit Planetary Isentropic-Coordinate (EPIC) Atmospheric Model: extract data file");
      nc_put_att_text(nc_id,NC_GLOBAL,"title",strlen(Message)+1,Message);
      sprintf(Message,"EPIC Model Version %4.2f: specified fields are extracted from the running model at timestep interval -itextract",
                       grid.data_version);
      nc_put_att_text(nc_id,NC_GLOBAL,"source",strlen(Message)+1,Message);
    }
    else {
      sprintf(Message,"Explicit Planetary Isentropic-Coordinate (EPIC) Atmospheric Model, Version %4.2f: runable data file",grid.data_version);
      nc_put_att_text(nc_id,NC_GLOBAL,"title",strlen(Message)+1,Message);
      sprintf(Message,"EPIC Model Version %4.2f: model state (prognostic variables and tendencies)",
                       grid.data_version);
      nc_put_att_text(nc_id,NC_GLOBAL,"source",strlen(Message)+1,Message);
    }
    sprintf(Message,"Developed at the Comparative Planetology Laboratory (CPL), Louisville, KY");
    nc_put_att_text(nc_id,NC_GLOBAL,"institution",strlen(Message)+1,Message);
   /* 
    * NOTE: Regarding time stamping the  global attribute "history," gcc 3.3.5 has a bug that generates a segmentation fault
    *       with time() and ctime(), so we are currently not including a timestamp.
    */   
    sprintf(Message,"This file produced by the EPIC Atmospheric model");
    nc_put_att_text(nc_id,NC_GLOBAL,"history",strlen(Message)+1,Message);
    sprintf(Message,"The EPIC Model is downloadable as open source from the NASA Planetary Data System (PDS) Atmospheres Node\n");
    nc_put_att_text(nc_id,NC_GLOBAL,"references",strlen(Message)+1,Message);
    
    WRITEF(&grid.data_version,grid_data_version,1);

    WRITETIME(&var.start_time,var_start_time);

    WRITEI(&planet->index,planet_index,1);
    WRITEC(planet->name,planet_name,32);
    WRITEC(planet->type,planet_type,16);
    WRITEC(planet->orbital_epoch,planet_orbital_epoch,8);
    WRITEF(&planet->re,planet_re,1);
    WRITEF(&planet->rp,planet_rp,1);
    WRITEF(&planet->omega_sidereal,planet_omega_sidereal,1);
    WRITEF(&planet->omega_synodic,planet_omega_synodic,1);
    WRITEF(&planet->cp,planet_cp,1);
    WRITEF(&planet->rgas,planet_rgas,1);
    WRITEF(&planet->GM,planet_GM,1);
    WRITEF(&planet->J2,planet_J2,1);
    WRITEF(&planet->x_he,planet_x_he,1);
    WRITEF(&planet->x_h2,planet_x_h2,1);
    WRITEF(&planet->x_3,planet_x_3,1);
    WRITEF(&planet->a,planet_a,1);
    WRITEF(&planet->e,planet_e,1);
    WRITEF(&planet->i,planet_i,1);
    WRITEF(&planet->lon_ascending_node,planet_lon_ascending_node,1);
    WRITEF(&planet->lon_perihelion,planet_lon_perihelion,1);
    WRITEF(&planet->mean_lon,planet_mean_lon,1);
    WRITEF(&planet->orbit_period,planet_orbit_period,1);
    WRITEF(&planet->vernal_equinox_anomaly,planet_vernal_equinox_anomaly,1);
    WRITEF(&planet->kinvisc,planet_kinvisc,1);
    WRITEF(&planet->dynvisc,planet_dynvisc,1);
    WRITEF(&planet->k_a,planet_k_a,1);

    WRITEI(&grid.nk,grid_nk,1);
    WRITEI(&grid.nj,grid_nj,1);

    if (stretch_ni) {
      WRITEI(&stretch_ni,grid_ni,1);
    }
    else {
      WRITEI(&grid.ni,grid_ni,1);
    }

    WRITEI(&grid.jtp,grid_jtp,1);

    WRITEC(grid.geometry,grid_geometry,GEOM_STR);
    WRITEC(grid.uv_timestep_scheme,grid_uv_timestep_scheme,N_STR);
    WRITEC(grid.turbulence_scheme,grid_turbulence_scheme,N_STR);
    WRITEC(grid.stability_factor,grid_stability_factor,N_STR);

    WRITEC(var.hdry.advection_scheme,var_hdry_advection_scheme,N_STR);
    WRITEC(var.theta.advection_scheme,var_theta_advection_scheme,N_STR);
    WRITEC(var.fpara.advection_scheme,var_fpara_advection_scheme,N_STR);
    WRITEC(var.nu_turb.advection_scheme,var_nu_turb_advection_scheme,N_STR);
    WRITEC(var.species[FIRST_SPECIES].advection_scheme,var_species_advection_scheme,N_STR);

    WRITEF(&grid.globe_lonbot,grid_globe_lonbot,1);
    WRITEF(&grid.globe_lontop,grid_globe_lontop,1);
    WRITEF(&grid.globe_latbot,grid_globe_latbot,1);
    WRITEF(&grid.globe_lattop,grid_globe_lattop,1);
    WRITEC(grid.f_plane_map,grid_f_plane_map,GEOM_STR);       
    WRITEF(&grid.f_plane_lat0,grid_f_plane_lat0,1);
    WRITEF(&grid.f_plane_half_width,grid_f_plane_half_width,1);

    WRITEI(grid.wrap,grid_wrap,TOPDIM);
    WRITEI(grid.pad,grid_pad,TOPDIM);
    WRITEI(&grid.jlo,grid_jlo,1);
    WRITEI(&grid.jfirst,grid_jfirst,1);
    WRITEI(&grid.jlast,grid_jlast,1);
    WRITEI(&grid.k_sponge,grid_k_sponge,1);
    WRITEI(&grid.k_sigma,grid_k_sigma,1);
    WRITEF(&grid.du_vert,grid_du_vert,1);

    WRITEI(&grid.newt_cool_on,grid_newt_cool_on,1);
    WRITEI(&grid.microphysics_on,grid_microphysics_on,1);
    WRITEI(&grid.nmt_physics_on,grid_nmt_physics_on,1);
    WRITEI(var.on_list,var_on_list,LAST_INDEX-FIRST_INDEX+1);
    WRITEI(&var.extract_on,var_extract_on,1);
    WRITEI(var.extract_on_list,var_extract_on_list,LAST_INDEX-FIRST_INDEX+1);

    WRITEI(&var.ntp,var_ntp,1);
    WRITEI(&var.n_t_cool,var_n_t_cool,1);

    if (grid.nmt_physics_on == 1) {
      WRITEI(&nmt.nz,nmt_nz,1);
      WRITEI(&nmt.land,nmt_land,1);
      WRITEI(&nmt.dorad,nmt_dorad,1);
      WRITEI(&nmt.frad,nmt_frad,1);
      WRITEF(&nmt.radcool,nmt_radcool,1);
      WRITEF(&nmt.tpause,nmt_tpause,1);
      WRITEF(&nmt.radbrk,nmt_radbrk,1);
      WRITEF(&nmt.cvc,nmt_cvc,1);
      WRITEF(&nmt.cvs,nmt_cvs,1);
      WRITEF(&nmt.cvp,nmt_cvp,1);
      WRITEF(&nmt.cve,nmt_cve,1);
      WRITEF(&nmt.theslop,nmt_theslop,1);
      WRITEF(&nmt.pstiff,nmt_pstiff,1);
      WRITEF(&nmt.pscale,nmt_pscale,1);
      WRITEF(&nmt.pbltop,nmt_pbltop,1);
      WRITEF(&nmt.cdrag,nmt_cdrag,1);
      WRITEF(&nmt.wscale,nmt_wscale,1);
      WRITEF(&nmt.sst,nmt_sst,1);
      WRITEF(&nmt.lfrac,nmt_lfrac,1);
      WRITEF(&nmt.eflux0,nmt_eflux0,1);
      WRITEF(&nmt.cld,nmt_cld,1);
      WRITEF(&nmt.cfract,nmt_cfract,1);
      WRITEF(&nmt.dt,nmt_dt,1);
    }
  }

  if (portion == SIZE_DATA) {
    if (IAMNODE == 0) {
      nc_close(nc_id);
    }
    return;
  }

  /**********************************************
   *                                            *
   * Seam between SIZE_DATA and POST_SIZE_DATA. *
   *                                            *
   **********************************************/

  if (portion == POST_SIZE_DATA      ||
      portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA ||
      portion == ALL_DATA              ) {
    EPIC_FLOAT
      *buffer;

    nk = grid.nk;

    WRITEF(&grid.dlt,grid_dlt,1);
    WRITEF(&grid.dln,grid_dln,1);
    WRITEI(&grid.dt,grid_dt,1);

    WRITEF(&var.time_fp_bar,var_time_fp_bar,1);

    if (var.ntp > 0) {
      /*
       * Write reference temperature sounding profile.
       */
      buffer = fvector(0,var.ntp-1,dbmsname);
      for (i = 0; i < var.ntp; i++) {
        /*
         * Convert output pressure from Pa to hPa.
         */
        buffer[i] = var.pdat[i]/100.;
      }
      WRITEF(buffer,var_pdat,var.ntp);
      WRITEC("hPa",var_pdat_units,strlen("hPa")+1);
      WRITEF(var.tdat,var_tdat,var.ntp);
      WRITEC("K",var_tdat_units,strlen("K")+1);
      WRITEF(var.dtdat,var_dtdat,var.ntp);
      WRITEC("K",var_dtdat_units,strlen("K")+1);
      free_fvector(buffer,0,var.ntp-1,dbmsname);
    }

    if (var.n_t_cool > 0) {
      /*
       * Write Newtonian-cooling time constant profile.
       * Since float_triplet memory is striped, use a buffer.
       */
      buffer = fvector(0,var.n_t_cool-1,dbmsname);
      for (i = 0; i < var.n_t_cool; i++) {
        buffer[i] = exp(-var.t_cool_table[i].x)/100.;
      }
      WRITEF(buffer,var_t_cool_p,var.n_t_cool);
      WRITEC("hPa",var_t_cool_p_units,strlen("hPa")+1);
      for (i = 0; i < var.n_t_cool; i++) {
        buffer[i] = var.t_cool_table[i].y;
      }
      WRITEF(buffer,var_t_cool_tau,var.n_t_cool);
      WRITEC("s",var_t_cool_tau_units,strlen("s")+1);
      free_fvector(buffer,0,var.n_t_cool-1,dbmsname);
    }

    WRITED(grid.sigmatheta,grid_sigmatheta,2*(nk+1)+1);

    WRITEF(grid.p_ref,grid_p_ref,2*(nk+1)+1);
    WRITEF(grid.theta_ref,grid_theta_ref,2*(nk+1)+1);
    WRITEF(grid.h_min,grid_h_min,nk+1);

    WRITEF(&grid.sgth_bot,grid_sgth_bot,1);
    WRITEF(&grid.sgth_top,grid_sgth_top,1);
    /*
     * Double precision to improve calculation of diagnostic theta.
     */
    WRITED(&grid.zeta0,grid_zeta0,1);
    WRITED(&grid.zeta1,grid_zeta1,1);
    WRITED(&grid.hybrid_alpha,grid_hybrid_alpha,1);
    WRITED(&grid.sigma_sigma,grid_sigma_sigma,1);

    WRITEI(&grid.newt_cool_adjust,grid_newt_cool_adjust,1);
    WRITEI(&grid.diffusion_direction,grid_diffusion_direction,1);
    WRITEC(grid.eos,grid_eos,8);
    WRITEF(&grid.prandtl,grid_prandtl,1);
    WRITEF(&grid.tau_drag,grid_tau_drag,1);
    WRITEI(&grid.drag_v,grid_drag_v,1);
    WRITEI(&grid.drag_zonal_avg,grid_drag_zonal_avg,1);

    WRITEF(&grid.pbot,grid_pbot,1);
    WRITEF(&grid.press0,grid_press0,1);
    WRITEF(&grid.mont0,grid_mont0,1);
    WRITEF(&grid.topo_scale,grid_topo_scale,1);

    WRITEF(&grid.nudiv_nondim,grid_nudiv_nondim,1);
    WRITEF(grid.nu_nondim,grid_nu_nondim,MAX_NU_ORDER+1);
    WRITED(grid.nu,grid_nu,MAX_NU_ORDER+1);
  }

  if (portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA   ) {
    if (IAMNODE == 0) {
      nc_close(nc_id);
    }
    return;
  }

  /**********************************************************************
   *                                                                    *
   * Seam between HEADER_DATA/EXTRACT_HEADER_DATA and POST_HEADER_DATA. *
   *                                                                    *
   **********************************************************************/

  /*
   * Write non-array data that is time dependent.
   */
  WRITETIME(&var.model_time,var_model_time);

  if (stretch_ni) {
    grid.cfl_dt = cfl_dt(planet);
  }
  WRITEI(&grid.cfl_dt,grid_cfl_dt,1);

  WRITEI(&grid.aux_a,grid_aux_a,1);
  WRITEI(&grid.aux_b,grid_aux_b,1);
  WRITEI(&grid.aux_c,grid_aux_c,1);
  WRITEF(&grid.aux_fa,grid_aux_fa,1);
  WRITEF(&grid.aux_fb,grid_aux_fb,1);
  WRITEF(&grid.aux_fc,grid_aux_fc,1);

  if (IAMNODE == NODE0) {
    /*
     * Leave define mode for netCDF file:
     */
    nc_err = nc_enddef(nc_id);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"nc_enddef(), %s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
  }

  /*
   * Write time in days since var.start_time.
   */
  if (IAMNODE == NODE0) {
    t_index[ 0] = time_index;
    the_time[0] = TIME/86400.;

#if EPIC_PRECISION == DOUBLE_PRECISION
    nc_err = nc_put_var1_double(nc_id,var.info[0].coorid[NETCDF_T_INDEX],t_index,the_time);
#else
    nc_err = nc_put_var1_float(nc_id,var.info[0].coorid[NETCDF_T_INDEX],t_index,the_time);
#endif

    if (nc_err != NC_NOERR) {
      sprintf(Message,"t_index=%zu, %s",t_index[0],nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
  }

  /*
   * Write solar longitude, L_s [deg].
   */
  if (IAMNODE == NODE0) {
    t_index[0] = time_index;
    the_L_s[0] = L_s;

#if EPIC_PRECISION == DOUBLE_PRECISION
    nc_err = nc_put_var1_double(nc_id,var.l_s.info.id,t_index,the_L_s);
#else
    nc_err = nc_put_var1_float(nc_id,var.l_s.info.id,t_index,the_L_s);
#endif

    if (nc_err != NC_NOERR) {
      sprintf(Message,"Writing L_s, t_index=%zu, %s",t_index[0],nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
  }

  /* 
   * Loop over nodes to get all the information off one node before
   * moving on to the next. 
   */
  num_nodes = setup_write_array();

  /*
   * Load start and end vectors.
   */
  start[0] = 1;
  /*
   * NOTE: Set end[0] to grid.ni even in stretch_ni case.
   */
  end[0] = grid.ni;

  /*
   * Write 0 to grid.nk for all arrays. This puts the top and bottom interfaces into the file
   * for variables defined on the layer interfaces (which is the majority of the arrays).
   */
  start[2] = 0;
  end[  2] = grid.nk;

  start[3] = time_index;
  end[  3] = time_index;

  if (IAMNODE == NODE0) {
    fprintf(stdout,"  0%%");
    fflush(stdout);
  }

  for (node = 0; node < num_nodes; node++) {
    if (IAMNODE == NODE0) {
      fprintf(stdout,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)node/grid.we_num_nodes));
      fflush(stdout);
    }

    get_jlohi(node,num_nodes,start+1,end+1);

    /*
     * Write parameter arrays.
     */
    if (strcmp(planet->type,"terrestrial") == 0) {
      if (var.gz_surface.on) {
        if ((portion == EXTRACT_DATA && var.gz_surface.extract_on) ||
            (portion != EXTRACT_DATA)                                ) {
          write_array(node,TWODIM,start,end,stretch_ni,var.gz_surface.info[0].name,
                      var.gz_surface.info[0].index,var.gz_surface.value,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
    }

    if (var.u_spinup.on) {
      if ((portion == EXTRACT_DATA && var.u_spinup.extract_on) ||
          (portion != EXTRACT_DATA)                              ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.u_spinup.info[0].name,
                    var.u_spinup.info[0].index,var.u_spinup.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (portion == EXTRACT_DATA) {
      /*
       * Write selected diagnostic variables to extract file.
       */
      if (var.hdry3.on && var.hdry3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.hdry3.info[0].name,
                    var.hdry3.info[0].index,var.hdry3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.pdry3.on && var.pdry3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.pdry3.info[0].name,
                    var.pdry3.info[0].index,var.pdry3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.p2.on && var.p2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.p2.info[0].name,
                    var.p2.info[0].index,var.p2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      /*
       * P3 is written as if it is the prognostic variable, rather than HDRY,
       * so we handle HDRY here instead of P3.
       */
      if (var.hdry.on && var.hdry.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.hdry.info[0].name,
                    var.hdry.info[0].index,var.hdry.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.theta2.on && var.theta2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.theta2.info[0].name,
                    var.theta2.info[0].index,var.theta2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.h2.on && var.h2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.h2.info[0].name,
                    var.h2.info[0].index,var.h2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.h3.on && var.h3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.h3.info[0].name,
                    var.h3.info[0].index,var.h3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.t2.on && var.t2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.t2.info[0].name,
                    var.t2.info[0].index,var.t2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.t3.on && var.t3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.t3.info[0].name,
                    var.t3.info[0].index,var.t3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.rho2.on && var.rho2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.rho2.info[0].name,
                    var.rho2.info[0].index,var.rho2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.rho3.on && var.rho3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.rho3.info[0].name,
                    var.rho3.info[0].index,var.rho3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.exner3.on && var.exner3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.exner3.info[0].name,
                    var.exner3.info[0].index,var.exner3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.fgibb3.on && var.fgibb3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.fgibb3.info[0].name,
                    var.fgibb3.info[0].index,var.fgibb3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.gz3.on && var.gz3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.gz3.info[0].name,
                    var.gz3.info[0].index,var.gz3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.mont3.on && var.mont3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.mont3.info[0].name,
                    var.mont3.info[0].index,var.mont3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.heat3.on && var.heat3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.heat3.info[0].name,
                    var.heat3.info[0].index,var.heat3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.pv3.on && var.pv3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.pv3.info[0].name,
                    var.pv3.info[0].index,var.pv3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.ri2.on && var.ri2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.ri2.info[0].name,
                    var.ri2.info[0].index,var.ri2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.diffusion_coef_uv.on && var.diffusion_coef_uv.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.diffusion_coef_uv.info[0].name,
                    var.diffusion_coef_uv.info[0].index,var.diffusion_coef_uv.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.diffusion_coef_theta.on && var.diffusion_coef_theta.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.diffusion_coef_theta.info[0].name,
                    var.diffusion_coef_theta.info[0].index,var.diffusion_coef_theta.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.diffusion_coef_mass.on && var.diffusion_coef_mass.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.diffusion_coef_mass.info[0].name,
                    var.diffusion_coef_mass.info[0].index,var.diffusion_coef_mass.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.div_uv3.on && var.div_uv3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.div_uv3.info[0].name,
                    var.div_uv3.info[0].index,var.div_uv3.value,EPIC_FLOAT_ARRAY,nc_id);
      }

      /*
       * Collect together here output options that do not have permanent memory.
       * These variables must be calculated here.
       */
      if (var.vort3.extract_on) {
        if (!var.vort3.on) {
          /* Calculate variable here and store in BUFF3D memory. */
          var.vort3.value = buff3d;
          for (K = KLO; K <= KHI; K++) {
            vorticity(planet,ON_SIGMATHETA,RELATIVE,
                      var.u.value    +(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                      var.v.value    +(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                      NULL,
                      var.vort3.value+(K-Kshift)*Nelem2d);
          }
        }
        write_array(node,FOURDIM,start,end,stretch_ni,var.vort3.info[0].name,
                    var.vort3.info[0].index,var.vort3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.eddy_pv3.extract_on) {
        var.eddy_pv3.value = buff3d;
        for (K = KLO; K <= KHI; K++) {
          for (J = JLOPAD; J <= JHIPADPV; J++) {
            avg = 0.;
            for (I = ILO; I <= IHI; I++) {
              avg += PV3(K,J,I);
            }
            avg /= grid.ni;
            for (I = ILOPAD; I <= IHIPAD; I++) {
              EDDY_PV3(K,J,I) = PV3(K,J,I)-avg;
            }
          }
        }
        /* No need to apply bc_lateral() here. */
        write_array(node,FOURDIM,start,end,stretch_ni,var.eddy_pv3.info[0].name,
                    var.eddy_pv3.info[0].index,var.eddy_pv3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      /* End of collection of output options that do not have permanent memory. */

      if (var.w3.on && var.w3.extract_on) {
        /*
         * NOTE: W3 is one timestep behind all the other diagnostic variables at the point
         *       when the data are written to extract.nc, because of the way it is calculated.
         *       It is thus written to extract.nc here with its appropriately lagged time index.
         *       The main consequence is that the last frame of W3 in an extract.nc file
         *       is filled with a copy of the penultimate frame (to avoid a blank).
         */
        if (time_index > 0) {
          start[3] = end[3] = time_index-1;
          write_array(node,FOURDIM,start,end,stretch_ni,var.w3.info[0].name,
                      var.w3.info[0].index,var.w3.value,EPIC_FLOAT_ARRAY,nc_id);
          start[3] = end[3] = time_index;
        }
        /*
         * Write this lagged value into the last timeframe in extract.nc as a placeholder,
         * to avoid having garbage values in the last frame for W3.
         */
        start[3] = end[3] = time_index;
        write_array(node,FOURDIM,start,end,stretch_ni,var.w3.info[0].name,
                    var.w3.info[0].index,var.w3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.dzdt3.on && var.dzdt3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.dzdt3.info[0].name,
                    var.dzdt3.info[0].index,var.dzdt3.value,EPIC_FLOAT_ARRAY,nc_id);
      }

      if (grid.nmt_physics_on == 1) {
        if (var.dry_entropy.on && var.dry_entropy.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni, \
                      var.dry_entropy.info[0].name, \
                      var.dry_entropy.info[0].index,var.dry_entropy.value, \
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        if (var.moist_entropy.on && var.moist_entropy.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni, \
                      var.moist_entropy.info[0].name, \
                      var.moist_entropy.info[0].index,var.moist_entropy.value, \
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        if (var.sat_moist_entropy.on && var.sat_moist_entropy.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni, \
                      var.sat_moist_entropy.info[0].name, \
                      var.sat_moist_entropy.info[0].index, \
                      var.sat_moist_entropy.value, \
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        if (var.the_flux.on && var.the_flux.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni, \
                      var.the_flux.info[0].name, \
                      var.the_flux.info[0].index, \
                      var.the_flux.value, \
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        if (var.rt_flux.on && var.rt_flux.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni, \
                      var.rt_flux.info[0].name, \
                      var.rt_flux.info[0].index, \
                      var.rt_flux.value, \
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        if (var.u_flux.on && var.u_flux.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni, \
                      var.u_flux.info[0].name, \
                      var.u_flux.info[0].index, \
                      var.u_flux.value, \
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        if (var.v_flux.on && var.v_flux.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni, \
                      var.v_flux.info[0].name, \
                      var.v_flux.info[0].index, \
                      var.v_flux.value, \
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        if (var.convthrot.on && var.convthrot.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni, \
                      var.convthrot.info[0].name, \
                      var.convthrot.info[0].index, \
                      var.convthrot.value, \
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        if (var.fluxthrot.on && var.fluxthrot.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni, \
                      var.fluxthrot.info[0].name, \
                      var.fluxthrot.info[0].index, \
                      var.fluxthrot.value, \
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        if (var.rain_rate.on && var.rain_rate.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni, \
                      var.rain_rate.info[0].name, \
                      var.rain_rate.info[0].index, \
                      var.rain_rate.value, \
                      EPIC_FLOAT_ARRAY,nc_id);
        }
      }
    }

    /* 
     * Write prognostic variables, and their associated tendency and moment data as appropriate.
     */
    if (var.u.on) {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        if ((portion == EXTRACT_DATA && var.u.extract_on) ||
            (portion != EXTRACT_DATA                    )   ) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.u.info[0].name,
                      var.u.info[0].index,var.u.value,EPIC_FLOAT_ARRAY,nc_id);
        }
        if ((portion != VAR_DATA    ) &&
            (portion != EXTRACT_DATA)   ) {
          write_array(node,THREEDIM,start,end,stretch_ni,var.u.info_tend[0].name,
                      var.u.info_tend[0].index,var.u.tendency+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          write_array(node,THREEDIM,start,end,stretch_ni,var.u.info_tend[1].name,
                      var.u.info_tend[1].index,var.u.tendency+IT_MINUS2*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        if (portion == EXTRACT_DATA && var.u.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.u.info[0].name,
                      var.u.info[0].index,var.u.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
        else if (portion != EXTRACT_DATA) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.u.info[0].name,
                      var.u.info[0].index,var.u.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          write_array(node,FOURDIM,start,end,stretch_ni,var.u.info[1].name,
                      var.u.info[1].index,var.u.value+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }

    if (var.v.on) {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        if ((portion == EXTRACT_DATA && var.v.extract_on) ||
            (portion != EXTRACT_DATA                    )   ) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.v.info[0].name,
                      var.v.info[0].index,var.v.value,EPIC_FLOAT_ARRAY,nc_id);
        }
        if ((portion != VAR_DATA    ) &&
            (portion != EXTRACT_DATA)   ) {
          write_array(node,THREEDIM,start,end,stretch_ni,var.v.info_tend[0].name,
                      var.v.info_tend[0].index,var.v.tendency+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          write_array(node,THREEDIM,start,end,stretch_ni,var.v.info_tend[1].name,
                      var.v.info_tend[1].index,var.v.tendency+IT_MINUS2*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        if (portion == EXTRACT_DATA && var.v.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.v.info[0].name,
                      var.v.info[0].index,var.v.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
        else if (portion != EXTRACT_DATA) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.v.info[0].name,
                      var.v.info[0].index,var.v.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          write_array(node,FOURDIM,start,end,stretch_ni,var.v.info[1].name,
                      var.v.info[1].index,var.v.value+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }

    /*
     * We write the closely related diagnostic variable P3 instead of the prognostic variable HDRY, 
     * since pressure is more intuitive and HDRY does not contain the pressure data at the top of the model.
     */
    if (var.p3.on) {
      if ((portion == EXTRACT_DATA && var.p3.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.p3.info[0].name,
                    var.p3.info[0].index,var.p3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.theta.on) {
      if ((portion == EXTRACT_DATA && var.theta.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.theta.info[0].name,
                    var.theta.info[0].index,var.theta.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.fpara.on) {
      if ((portion == EXTRACT_DATA && var.fpara.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.fpara.info[0].name,
                    var.fpara.info[0].index,var.fpara.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    for (iq = 0; iq < grid.nq; iq++) {
      if ((portion == EXTRACT_DATA && var.species[grid.is[iq]].phase[grid.ip[iq]].extract_on) ||
          (portion != EXTRACT_DATA                                        )   ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.species[grid.is[iq]].phase[grid.ip[iq]].info[0].name,
                    var.species[grid.is[iq]].phase[grid.ip[iq]].info[0].index,
                    var.species[grid.is[iq]].phase[grid.ip[iq]].q,
                    EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.nu_turb.on) {
      if ((portion == EXTRACT_DATA && var.nu_turb.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.nu_turb.info[0].name,
                    var.nu_turb.info[0].index,var.nu_turb.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }
  }

  /* 
   * Close NODE0's data connection.
   */
  if (IAMNODE == NODE0) {
    nc_close(nc_id);
    fprintf(stdout,"\b\b\b\b%3d%%\n",100);
    fflush(stdout);
  }

  return;
}

/*======================= end of var_write() ==================================*/

/*====================== number_spots_in_file() =====================================*/

int number_spots_in_file( char *spots_file )
{
  int 
    nspots = 0;
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
    dbmsname[]="number_spots_in_file";

  spots = fopen(spots_file,"r");
  if (!spots) {
    sprintf(Message,"unable to open %s \n",spots_file);
    epic_error(dbmsname,Message);
  }
  while (nspots == 0) {
    fgets(buffer,FILE_STR,spots);
    char_pt = strchr(buffer,':');
    if (char_pt) {
      sscanf(char_pt+1,"%d",&nspots);
    }
  }
  fclose(spots);

  return nspots;
}

/*====================== end of number_spots_in_file() =====================================*/

/*====================== read_spots_file() =====================================*/

void read_spots_file( 
               char       *spots_file,
               EPIC_FLOAT *ampspot,
               EPIC_FLOAT *lonspot,
               EPIC_FLOAT *latspot,
               EPIC_FLOAT *pspot,
               EPIC_FLOAT *aspot,
               EPIC_FLOAT *bspot,
               EPIC_FLOAT *cspot_up,
               EPIC_FLOAT *cspot_down,
               int        adjust_amplitude )
{
  register int 
    ispot;
  int
    nspots=0;
  EPIC_FLOAT
    fspot,
    factor;
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
    dbmsname[]="read_spots_file";

  spots = fopen(spots_file,"r");
  while (nspots == 0) {
    fgets(buffer,FILE_STR,spots);
    char_pt = strchr(buffer,':');
    if (char_pt) {
      sscanf(char_pt+1,"%d",&nspots);
    }
  }

  fgets(buffer,FILE_STR,spots);
  for (ispot = 0; ispot < nspots; ispot++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(spots,"%lf %lf %lf %lf %lf %lf %lf %lf",
           lonspot+ispot,latspot+ispot,pspot+ispot,
           aspot+ispot,bspot+ispot,cspot_up+ispot,cspot_down+ispot,ampspot+ispot);
#else
    fscanf(spots,"%f %f %f %f %f %f %f %f",
           lonspot+ispot,latspot+ispot,pspot+ispot,
           aspot+ispot,bspot+ispot,cspot_up+ispot,cspot_down+ispot,ampspot+ispot);
#endif
    /* Convert pspot from hPa to Pa: */
    pspot[ispot] *= 100.;


    /* 
     * Convert ampspot to amp for mont.  The factor named 'factor' 
     * should make the maximum spot velocity close to the input 
     * ampspot[m/s] for a gaussian spot mont.
     */
    if (adjust_amplitude == ADJUST_SPOT_AMPLITUDE) {
      fspot           = 2.*planet->omega_sidereal*sin(latspot[ispot]*DEG);
      factor          = 1.166;
      //factor          = 0.65;
      ampspot[ispot] *= factor*bspot[ispot]*DEG*planet->re*fabs(fspot);
    }

    fgets(buffer,FILE_STR,spots);
  }
  fclose(spots);

}

/*====================== end of read_spots_file() =====================================*/

/*======================= lookup_netcdf() =====================================*/

/*
 * Read numbers and names of global attributes and variables 
 * contained in infile, which must be in netCDF format.
 * Returns NC_NOERR if no read error, otherwise returns nc_err.
 *
 * NOTE: ngatts,**gattname,num_progs,**varname should be declared static in the calling
 *       function, with their input values equal to the last call, in order to
 *       properly reallocate memory.
 */

int lookup_netcdf(char   *infile,
                  int    *nc_id,
                  int    *ngatts,
                  char ***gattname,
                  int    *num_vars,
                  char ***varname)
{
  int
    i,
    ndims,unlimdimid,
    nc_err = NC_NOERR;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lookup_netcdf";

  /*
   * Free previous memory:
   */
  for (i = 0; i < *ngatts; i++) {
    free((*gattname)[i]);
  }
  for (i = 0; i < *num_vars; i++) {
    free((*varname)[i]);
  }

  if (IAMNODE == NODE0) {
    /*
     * Decompress infile if necessary.
     */
    handle_file_compression(infile);
    nc_err = nc_open(infile,NC_NOWRITE,nc_id);
    if (nc_err == NC_NOERR) {
      nc_inq(*nc_id,&ndims,num_vars,ngatts,&unlimdimid);
    }
    else {
      return nc_err;
    }
  }

#if defined(EPIC_MPI)
  MPI_Bcast(&nc_err,1,MPI_INT,NODE0,para.comm);
#endif

#if defined(EPIC_MPI)
  MPI_Bcast(ngatts,   1,MPI_INT,NODE0,para.comm);
  MPI_Bcast(num_vars, 1,MPI_INT,NODE0,para.comm);
#endif

  /*
   * Reallocate memory for character arrays. 
   * Look up and store names.
   */
  *gattname = (char **)realloc(*gattname,(*ngatts)*sizeof(char *));
  for (i = 0; i < *ngatts; i++) {
    (*gattname)[i] = (char *)calloc(NC_MAX_NAME,sizeof(char));
    if (IAMNODE == NODE0) {
      nc_inq_attname(*nc_id,NC_GLOBAL,i,(*gattname)[i]);
    }

#if defined (EPIC_MPI)
    MPI_Bcast((*gattname)[i],NC_MAX_NAME,MPI_CHAR,NODE0,para.comm);
#endif 

  }

  *varname = (char **)realloc(*varname,(*num_vars)*sizeof(char *));
  for (i = 0; i < *num_vars; i++) {
    (*varname)[i] = (char *)calloc(NC_MAX_NAME,sizeof(char));
    if (IAMNODE == NODE0) {
      nc_inq_varname(*nc_id,i,(*varname)[i]);
    }

#if defined (EPIC_MPI)
    MPI_Bcast((*varname)[i],NC_MAX_NAME,MPI_CHAR,NODE0,para.comm);
#endif

  } 

  /* No errors. */
  return nc_err;
}

/*======================= end of lookup_netcdf() ==============================*/

/*======================= define_netcdf() =====================================*/
/*
 * Define variables and attributes for netCDF file.
 * We are following the CF conventions.
 */

#define DEFINE_NC_GRID(u) \
    var.u.info[0].dim = FOURDIM; \
    /* lon (I direction) */ \
    sprintf(coord_name,"lon_%s",var.u.info[0].name); \
    nc_err = nc_def_dim(nc_id,coord_name,idim,&var.u.info[0].dimid[NETCDF_I_INDEX]); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_def_var(nc_id,coord_name,float_type,ONEDIM, \
                        &var.u.info[0].dimid[NETCDF_I_INDEX],&var.u.info[0].coorid[NETCDF_I_INDEX]); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],"standard_name", \
                             strlen("longitude")+1,"longitude"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],"long_name", \
                             strlen("longitude")+1,"longitude"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],"units", \
                             strlen("degrees_east")+1,"degrees_east"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    /* lat (J direction) */ \
    sprintf(coord_name,"lat_%s",var.u.info[0].name); \
    nc_err = nc_def_dim(nc_id,coord_name,jdim,&var.u.info[0].dimid[NETCDF_J_INDEX]); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_def_var(nc_id,coord_name,float_type,ONEDIM, \
                        &var.u.info[0].dimid[NETCDF_J_INDEX],&var.u.info[0].coorid[NETCDF_J_INDEX]); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],"standard_name", \
                             strlen("latitude")+1,"latitude"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],"long_name", \
                             strlen("latitude (planetographic)")+1,"latitude (planetographic)"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],"units", \
                             strlen("degrees_north")+1,"degrees_north"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    /* sigmatheta (K direction) */ \
    sprintf(coord_name,"sigmatheta_%s",var.u.info[0].name); \
    nc_err = nc_def_dim(nc_id,coord_name,kdim,&var.u.info[0].dimid[NETCDF_K_INDEX]); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_def_var(nc_id,coord_name,NC_DOUBLE,ONEDIM, \
                        &var.u.info[0].dimid[NETCDF_K_INDEX],&var.u.info[0].coorid[NETCDF_K_INDEX]); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"long_name", \
                             strlen("hybrid sigma-theta vertical coordinate")+1,"hybrid sigma-theta vertical coordinate"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"positive", \
                             strlen("up")+1,"up"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"units", \
                             strlen("K")+1,"K"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    /* time */ \
    var.u.info[0].dimid[ NETCDF_T_INDEX] = var.info[0].dimid[NETCDF_T_INDEX]; \
    var.u.info[0].coorid[NETCDF_T_INDEX] = var.info[0].coorid[NETCDF_T_INDEX];

#define DEFINE_NC_VAR(p2,hdry,num,has_standard_name) \
    for (in = 0; in < num; in++) { \
      var.p2.info[in].dim = var.hdry.info[0].dim;     /* rhs index should be 0 */ \
      for (i = 0; i < var.p2.info[in].dim; i++) { \
        var.p2.info[in].dimid[i]  = var.hdry.info[0].dimid[i]; \
        var.p2.info[in].coorid[i] = var.hdry.info[0].coorid[i]; \
      } \
      nc_err = nc_def_var(nc_id,var.p2.info[in].name,float_type,var.p2.info[in].dim, \
                          var.p2.info[in].dimid,&var.p2.info[in].id); \
      if (nc_err != NC_NOERR) { \
        fprintf(stderr,"DEFINE_NC_VAR: %s\n",nc_strerror(nc_err)); \
      } \
      if (has_standard_name) { \
        nc_err = nc_put_att_text(nc_id,var.p2.info[in].id,"standard_name", \
                                 strlen(var.p2.info[in].standard_name)+1,var.p2.info[in].standard_name); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_VAR: %s\n",nc_strerror(nc_err)); \
        } \
      } \
      nc_err = nc_put_att_text(nc_id,var.p2.info[in].id,"long_name", \
                               strlen(var.p2.info[in].long_name)+1,var.p2.info[in].long_name); \
      if (nc_err != NC_NOERR) { \
        fprintf(stderr,"DEFINE_NC_VAR: %s\n",nc_strerror(nc_err)); \
      } \
      nc_err = nc_put_att_text(nc_id,var.p2.info[in].id,"units", \
                               strlen(var.p2.info[in].units)+1,var.p2.info[in].units); \
      if (nc_err != NC_NOERR) { \
        fprintf(stderr,"DEFINE_NC_VAR: %s\n",nc_strerror(nc_err)); \
      } \
    }

#define DEFINE_NC_JI(gz_surface,p3,has_standard_name) \
    nc_err = nc_def_var(nc_id,var.gz_surface.info[0].name,float_type,TWODIM, \
                        &var.p3.info[0].dimid[NETCDF_J_INDEX],&var.gz_surface.info[0].id); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_JI: %s\n",nc_strerror(nc_err)); \
    } \
    if (has_standard_name) { \
      nc_err = nc_put_att_text(nc_id,var.gz_surface.info[0].id,"standard_name", \
                               strlen(var.gz_surface.info[0].standard_name)+1,var.gz_surface.info[0].standard_name); \
      if (nc_err != NC_NOERR) { \
        fprintf(stderr,"DEFINE_NC_JI: %s\n",nc_strerror(nc_err)); \
      } \
    } \
    nc_err = nc_put_att_text(nc_id,var.gz_surface.info[0].id,"long_name", \
                             strlen(var.gz_surface.info[0].long_name)+1,var.gz_surface.info[0].long_name); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_JI: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.gz_surface.info[0].id,"units", \
                             strlen(var.gz_surface.info[0].units)+1,var.gz_surface.info[0].units); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_JI: %s\n",nc_strerror(nc_err)); \
    }


#define DEFINE_NC_TEND(u,num) \
    for (in = 0; in < num; in++) { \
      nc_err = nc_def_var(nc_id,var.u.info_tend[in].name,float_type,THREEDIM, \
                          &var.u.info[0].dimid[NETCDF_K_INDEX],&var.u.info_tend[in].id); \
      if (nc_err != NC_NOERR) { \
        fprintf(stderr,"DEFINE_NC_TEND: %s\n",nc_strerror(nc_err)); \
      } \
      nc_err = nc_put_att_text(nc_id,var.u.info_tend[in].id,"units", \
                               strlen(var.u.info_tend[in].units)+1,var.u.info_tend[in].units); \
      if (nc_err != NC_NOERR) { \
        fprintf(stderr,"DEFINE_NC_TEND: %s\n",nc_strerror(nc_err)); \
      } \
      nc_err = nc_put_att_text(nc_id,var.u.info_tend[in].id,"standard_name", \
                               strlen(var.u.info_tend[in].standard_name)+1,var.u.info_tend[in].standard_name); \
      if (nc_err != NC_NOERR) { \
        fprintf(stderr,"DEFINE_NC_TEND: %s\n",nc_strerror(nc_err)); \
      } \
      nc_err = nc_put_att_text(nc_id,var.u.info_tend[in].id,"long_name", \
                               strlen(var.u.info_tend[in].long_name)+1,var.u.info_tend[in].long_name); \
      if (nc_err != NC_NOERR) { \
        fprintf(stderr,"DEFINE_NC_TEND: %s\n",nc_strerror(nc_err)); \
      } \
    }

void define_netcdf(planetspec   *planet,
                   char         *outfile,
                   int           portion,
                   int           stretch_ni,
                   int           nc_id)
{
  int
    nc_err,
    kdim,jdim,idim,
    K,J,I,
    is,ip,
    in,i;
  size_t
    index[1],
    message_len;
  char
    coord_name[VAR_NM_SZ+8];
  nc_type
    float_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="define_netcdf";

  if (IAMNODE != NODE0) {
    return;
  }

  /*
   * Layers 0 to grid.nk.
   */
  kdim = grid.nk+1;

  jdim = grid.nj-grid.jlo+1;

  if (stretch_ni) {
    idim = stretch_ni;
  }
  else {
    idim = grid.ni;
  }

  /*
   * Specify floating-point type (precision).
   */
  if (EPIC_PRECISION == DOUBLE_PRECISION) {
    float_type = NC_DOUBLE;
  }
  else {
    float_type = NC_FLOAT;
  }

  /*
   * The variables all share the same time dimension, which is unlimited.
   *
   * time:
   */
  nc_def_dim(nc_id,"time",NC_UNLIMITED,&var.info[0].dimid[NETCDF_T_INDEX]); 
  nc_def_var(nc_id,"time",float_type,ONEDIM,
             &var.info[0].dimid[NETCDF_T_INDEX],&var.info[0].coorid[NETCDF_T_INDEX]);
  message_len = strftime(Message,N_STR,"days since %Y-%m-%d %H:%M:%S 0",gmtime(&var.start_time));
  nc_err = nc_put_att_text(nc_id,var.info[0].coorid[NETCDF_T_INDEX],"units",
                           message_len+1,Message);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  nc_err = nc_put_att_text(nc_id,var.info[0].coorid[NETCDF_T_INDEX],"standard_name",
                           strlen("time")+1,"time");
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  nc_err = nc_put_att_text(nc_id,var.info[0].coorid[NETCDF_T_INDEX],"long_name",
                           strlen("elapsed time")+1,"elapsed time");
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  nc_err = nc_put_att_text(nc_id,var.info[0].coorid[NETCDF_T_INDEX],"calendar",
                           strlen("julian")+1,"julian");
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /*
   * Define solar longitude variable, L_s, which is a function of time.
   * L_s = 0 corresponds to the planet's vernal equinox.
   */
  var.l_s.info.dim       = ONEDIM;
  var.l_s.info.dimid[0]  = var.info[0].dimid[NETCDF_T_INDEX];
  var.l_s.info.coorid[0] = var.info[0].coorid[NETCDF_T_INDEX];
  nc_err = nc_def_var(nc_id,var.l_s.info.name,float_type,var.l_s.info.dim,
                      var.l_s.info.dimid,&var.l_s.info.id);
  if (nc_err != NC_NOERR) {
    fprintf(stderr,"defining L_s: %s\n",nc_strerror(nc_err));
  }
  nc_err = nc_put_att_text(nc_id,var.l_s.info.id,"long_name",
                           strlen(var.l_s.info.long_name)+1,var.l_s.info.long_name);
  if (nc_err != NC_NOERR) {
    fprintf(stderr,"defining L_s: %s\n",nc_strerror(nc_err));
  }
  nc_err = nc_put_att_text(nc_id,var.l_s.info.id,"units",
                             strlen(var.l_s.info.units)+1,var.l_s.info.units);
  if (nc_err != NC_NOERR) {
    fprintf(stderr,"defining L_s: %s\n",nc_strerror(nc_err));
  }

  /*
   * The variables u, v, hdry, p3, and pv3 each reside on a different staggered
   * grid in the C-grid scheme, and we refer to these as the u-grid, v-grid, h-grid, etc.
   * Thus, the grids for these variables, followed by the variables themselves,
   * are defined before the other variables that share these grids.
   */
  DEFINE_NC_GRID(pv3);
  DEFINE_NC_GRID(u);
  DEFINE_NC_GRID(v);
  DEFINE_NC_GRID(hdry);
  DEFINE_NC_GRID(p3);

  /*
   * Define the base variables from which the grids derive their names.
   */
  if (var.u.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.u.extract_on) {
        DEFINE_NC_VAR(u,u,1,HAS_STANDARD_NAME);
      }
    } 
    else {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        DEFINE_NC_VAR(u,u,1,HAS_STANDARD_NAME);
        if (portion == ALL_DATA) {
          DEFINE_NC_TEND(u,2);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        DEFINE_NC_VAR(u,u,2,HAS_STANDARD_NAME);
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }
  }

  if (var.v.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.v.extract_on) {
        DEFINE_NC_VAR(v,v,1,HAS_STANDARD_NAME);
      }
    } 
    else {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        DEFINE_NC_VAR(v,v,1,HAS_STANDARD_NAME);
        if (portion == ALL_DATA) {
          DEFINE_NC_TEND(v,2);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        DEFINE_NC_VAR(v,v,2,HAS_STANDARD_NAME);
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }
  }

  if (var.p3.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.p3.extract_on) {
        DEFINE_NC_VAR(p3,p3,1,HAS_STANDARD_NAME);
      }
    }
    else {
      DEFINE_NC_VAR(p3,p3,1,HAS_STANDARD_NAME);
    }
  }

  if (var.pv3.on) {
    if (portion == EXTRACT_HEADER_DATA && var.pv3.extract_on) {
      DEFINE_NC_VAR(pv3,pv3,1,HAS_STANDARD_NAME);
    }
  }

  if (var.hdry.on) {
    if (portion == EXTRACT_HEADER_DATA && var.hdry.extract_on) {
      DEFINE_NC_VAR(hdry,hdry,1,NEEDS_STANDARD_NAME);
    }
  }

  /*
   * The remaining variables.
   */
  if (var.theta.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.theta.extract_on) {
        DEFINE_NC_VAR(theta,p3,1,HAS_STANDARD_NAME);
      }
    }
    else {
      DEFINE_NC_VAR(theta,p3,1,HAS_STANDARD_NAME);
    }
  }

  if (var.fpara.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.fpara.extract_on) {
        DEFINE_NC_VAR(fpara,p3,1,NEEDS_STANDARD_NAME);
      }
    }
    else {
      DEFINE_NC_VAR(fpara,p3,1,NEEDS_STANDARD_NAME);
    }
  }

  /*
   * NOTE: Not using grid.nq, grid.is[], grid.ip[] here because of
   *       the extract_on conditional.
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
        if (var.species[is].extract_on) {
          for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
            if (var.species[is].phase[ip].on) {
              DEFINE_NC_VAR(species[is].phase[ip],p3,1,NEEDS_STANDARD_NAME);
            }
          }
        }
      }
      else {
        for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
          if (var.species[is].phase[ip].on) {
            DEFINE_NC_VAR(species[is].phase[ip],p3,1,NEEDS_STANDARD_NAME);
          }
        }
      }
    }
  }

  if (var.nu_turb.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.nu_turb.extract_on) {
        DEFINE_NC_VAR(nu_turb,hdry,1,NEEDS_STANDARD_NAME);
      }
    }
    else {
      DEFINE_NC_VAR(nu_turb,hdry,1,NEEDS_STANDARD_NAME);
    }
  }

  if (portion == EXTRACT_HEADER_DATA) {
    if (var.hdry3.extract_on) {
      DEFINE_NC_VAR(hdry3,p3,1,NEEDS_STANDARD_NAME);
    }
    if (var.pdry3.extract_on) {
      DEFINE_NC_VAR(pdry3,p3,1,NEEDS_STANDARD_NAME);
    }
    if (var.p2.extract_on) {
      DEFINE_NC_VAR(p2,hdry,1,NEEDS_STANDARD_NAME);
    }
    if (var.theta2.extract_on) {
      DEFINE_NC_VAR(theta2,hdry,1,HAS_STANDARD_NAME);
    }
    if (var.h2.extract_on) {
      DEFINE_NC_VAR(h2,hdry,1,NEEDS_STANDARD_NAME);
    }
    if (var.h3.extract_on) {
      DEFINE_NC_VAR(h3,p3,1,NEEDS_STANDARD_NAME);
    }
    if (var.t2.extract_on) {
      DEFINE_NC_VAR(t2,hdry,1,HAS_STANDARD_NAME);
    }
    if (var.t3.extract_on) {
      DEFINE_NC_VAR(t3,p3,1,HAS_STANDARD_NAME);
    }
    if (var.rho2.extract_on) {
      DEFINE_NC_VAR(rho2,hdry,1,HAS_STANDARD_NAME);
    }
    if (var.rho3.extract_on) {
      DEFINE_NC_VAR(rho3,p3,1,HAS_STANDARD_NAME);
    }
    if (var.exner3.extract_on) {
      DEFINE_NC_VAR(exner3,p3,1,NEEDS_STANDARD_NAME);
    }
    if (var.fgibb3.extract_on) {
      DEFINE_NC_VAR(fgibb3,p3,1,NEEDS_STANDARD_NAME); 
    }
    if (var.gz3.extract_on) {
      DEFINE_NC_VAR(gz3,p3,1,HAS_STANDARD_NAME);
    }
    if (var.mont3.extract_on) {
      DEFINE_NC_VAR(mont3,p3,1,NEEDS_STANDARD_NAME);
    }
    if (var.ri2.extract_on) {
      DEFINE_NC_VAR(ri2,hdry,1,NEEDS_STANDARD_NAME);
    }
    if (var.heat3.extract_on) {
      DEFINE_NC_VAR(heat3,p3,1,NEEDS_STANDARD_NAME);
    }
    if (var.vort3.extract_on) {
      DEFINE_NC_VAR(vort3,pv3,1,HAS_STANDARD_NAME);
    }
    if (var.eddy_pv3.extract_on) {
      DEFINE_NC_VAR(eddy_pv3,pv3,1,NEEDS_STANDARD_NAME);
    }
    if (var.div_uv3.extract_on) {
      /*
       * NOTE: DIV_UV3 is lagged in leapfrog case via grid.it_uv_dis.
       */
      DEFINE_NC_VAR(div_uv3,p3,1,HAS_STANDARD_NAME);
    }
    if (var.w3.extract_on) {
      DEFINE_NC_VAR(w3,p3,1,NEEDS_STANDARD_NAME);
    }
    if (var.dzdt3.extract_on) {
      DEFINE_NC_VAR(dzdt3,p3,1,HAS_STANDARD_NAME);
    }
    if (var.diffusion_coef_uv.extract_on) {
      DEFINE_NC_VAR(diffusion_coef_uv,hdry,1,NEEDS_STANDARD_NAME);
    }
    if (var.diffusion_coef_theta.extract_on) {
      /* Carried on h-grid, even though THETA is carried on p-grid. */
      DEFINE_NC_VAR(diffusion_coef_theta,hdry,1,NEEDS_STANDARD_NAME);
    }
    if (var.diffusion_coef_mass.extract_on) {
      /* Carried on h-grid, even though Qs are carried on p-grid. */
      DEFINE_NC_VAR(diffusion_coef_mass,hdry,1,NEEDS_STANDARD_NAME);
    }

    if (grid.nmt_physics_on == 1) {
      if (var.dry_entropy.extract_on) {
        /* Carried on h-grid, even though Qs are carried on p-grid. */
        DEFINE_NC_VAR(dry_entropy,hdry,1,NEEDS_STANDARD_NAME);
      }
      if (var.moist_entropy.extract_on) {
        /* Carried on h-grid, even though Qs are carried on p-grid. */
        DEFINE_NC_VAR(moist_entropy,hdry,1,NEEDS_STANDARD_NAME);
      }
      if (var.sat_moist_entropy.extract_on) {
        /* Carried on h-grid, even though Qs are carried on p-grid. */
        DEFINE_NC_VAR(sat_moist_entropy,hdry,1,NEEDS_STANDARD_NAME);
      }
      if (var.the_flux.extract_on) {
        /* Carried on h-grid, even though Qs are carried on p-grid. */
        DEFINE_NC_VAR(the_flux,hdry,1,NEEDS_STANDARD_NAME);
      }
      if (var.rt_flux.extract_on) {
        /* Carried on h-grid, even though Qs are carried on p-grid. */
        DEFINE_NC_VAR(rt_flux,hdry,1,NEEDS_STANDARD_NAME);
      }
      if (var.u_flux.extract_on) {
        /* Carried on h-grid, even though Qs are carried on p-grid. */
        DEFINE_NC_VAR(u_flux,hdry,1,NEEDS_STANDARD_NAME);
      }
      if (var.v_flux.extract_on) {
        /* Carried on h-grid, even though Qs are carried on p-grid. */
        DEFINE_NC_VAR(v_flux,hdry,1,NEEDS_STANDARD_NAME);
      }
      if (var.convthrot.extract_on) {
        /* Carried on h-grid, even though Qs are carried on p-grid. */
        DEFINE_NC_VAR(convthrot,hdry,1,NEEDS_STANDARD_NAME);
      }
      if (var.fluxthrot.extract_on) {
        /* Carried on h-grid, even though Qs are carried on p-grid. */
        DEFINE_NC_VAR(fluxthrot,hdry,1,NEEDS_STANDARD_NAME);
      }
      if (var.rain_rate.extract_on) {
        /* Carried on h-grid, even though Qs are carried on p-grid. */
        DEFINE_NC_VAR(rain_rate,hdry,1,HAS_STANDARD_NAME);
      }
    }
  }

  /*
   * Parameter values.
   */
  if (var.u_spinup.on) {
    if ((portion == EXTRACT_HEADER_DATA && var.u_spinup.extract_on) ||
        (portion != EXTRACT_HEADER_DATA)                              ) {
      DEFINE_NC_VAR(u_spinup,u,1,NEEDS_STANDARD_NAME);
    }
  }

  if (strcmp(planet->type,"terrestrial") == 0) {
    if (var.gz_surface.on) {
      if ((portion == EXTRACT_HEADER_DATA && var.gz_surface.extract_on) ||
          (portion != EXTRACT_HEADER_DATA)                              ) {
        DEFINE_NC_JI(gz_surface,p3,HAS_STANDARD_NAME);
      }
    }
  }

  /*
   * Leave define mode:
   */
  nc_enddef(nc_id);

  /*
   * Assign values to staggered C-grid coordinates.
   */

  /*
   * longitude:
   */
  if (stretch_ni) {
    /*
     * Need to compute and output longitudes appropriate to stretch_ni.
     */
    EPIC_FLOAT
      *longitude;
    int
      ii;

    /*
     * Allocate memory.
     */
    longitude = fvector(0,2*(stretch_ni+1),dbmsname);

    /*
     * NOTE: The value of grid.dln should already be modified accordingly.
     */

    if (strcmp(grid.geometry,"globe") == 0) {
      longitude[0] = grid.globe_lonbot-grid.dln;
    }
    else {
      longitude[0] = -180.-grid.dln;
    }
    for (ii = 1; ii <= 2*(stretch_ni+1); ii++) {
      longitude[ii] = longitude[ii-1]+grid.dln*.5;
    }

    for (I = 1; I <= stretch_ni; I++) {
      index[0] = I-1;

#if EPIC_PRECISION == DOUBLE_PRECISION
      nc_put_var1_double(nc_id,var.pv3.info[0].coorid[NETCDF_I_INDEX],
                         index,&(longitude[2*I]));
      nc_put_var1_double(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],
                         index,&(longitude[2*I]));
      nc_put_var1_double(nc_id,var.v.info[0].coorid[NETCDF_I_INDEX],
                         index,&(longitude[2*I+1]));
      nc_put_var1_double(nc_id,var.hdry.info[0].coorid[NETCDF_I_INDEX],
                         index,&(longitude[2*I+1]));
      nc_put_var1_double(nc_id,var.p3.info[0].coorid[NETCDF_I_INDEX],
                         index,&(longitude[2*I+1]));
#else
      nc_put_var1_float(nc_id,var.pv3.info[0].coorid[NETCDF_I_INDEX],
                        index,&(longitude[2*I]));
      nc_put_var1_float(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],
                        index,&(longitude[2*I]));
      nc_put_var1_float(nc_id,var.v.info[0].coorid[NETCDF_I_INDEX],
                        index,&(longitude[2*I+1]));
      nc_put_var1_float(nc_id,var.hdry.info[0].coorid[NETCDF_I_INDEX],
                        index,&(longitude[2*I+1]));
      nc_put_var1_float(nc_id,var.p3.info[0].coorid[NETCDF_I_INDEX],
                        index,&(longitude[2*I+1]));
#endif
    }

    /*
     * Free allocated memory.
     */
    free_fvector(longitude,0,2*(stretch_ni+1),dbmsname);
  }
  else {
    for (I = 1; I <= grid.ni; I++) {
      index[0] = I-1;

#if EPIC_PRECISION == DOUBLE_PRECISION
      nc_put_var1_double(nc_id,var.pv3.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I]));
      nc_put_var1_double(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I]));
      nc_put_var1_double(nc_id,var.v.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I+1]));
      nc_put_var1_double(nc_id,var.hdry.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I+1]));
      nc_put_var1_double(nc_id,var.p3.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I+1]));
#else
      nc_put_var1_float(nc_id,var.pv3.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I]));
      nc_put_var1_float(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I]));
      nc_put_var1_float(nc_id,var.v.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I+1]));
      nc_put_var1_float(nc_id,var.hdry.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I+1]));
      nc_put_var1_float(nc_id,var.p3.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I+1]));
#endif
    }
  }

  /*
   * latitude:
   */
  for (J = grid.jlo; J <= grid.nj; J++) {
    index[0] = J-grid.jlo;

#if EPIC_PRECISION == DOUBLE_PRECISION
    nc_put_var1_double(nc_id,var.pv3.info[0].coorid[NETCDF_J_INDEX],
                       index,&(grid.lat[2*J]));
    nc_put_var1_double(nc_id,var.v.info[0].coorid[NETCDF_J_INDEX],
                       index,&(grid.lat[2*J]));
    nc_put_var1_double(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],
                       index,&(grid.lat[2*J+1]));
    nc_put_var1_double(nc_id,var.hdry.info[0].coorid[NETCDF_J_INDEX],
                       index,&(grid.lat[2*J+1]));
    nc_put_var1_double(nc_id,var.p3.info[0].coorid[NETCDF_J_INDEX],
                       index,&(grid.lat[2*J+1]));
#else
    nc_put_var1_float(nc_id,var.pv3.info[0].coorid[NETCDF_J_INDEX],
                      index,&(grid.lat[2*J]));
    nc_put_var1_float(nc_id,var.v.info[0].coorid[NETCDF_J_INDEX],
                      index,&(grid.lat[2*J]));
    nc_put_var1_float(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],
                      index,&(grid.lat[2*J+1]));
    nc_put_var1_float(nc_id,var.hdry.info[0].coorid[NETCDF_J_INDEX],
                      index,&(grid.lat[2*J+1]));
    nc_put_var1_float(nc_id,var.p3.info[0].coorid[NETCDF_J_INDEX],
                      index,&(grid.lat[2*J+1]));
#endif
  }

  /*
   * sigmatheta:
   */
  for (K = 0; K <= grid.nk; K++) {
    index[0] = K;
    nc_put_var1_double(nc_id,var.pv3.info[0].coorid[NETCDF_K_INDEX],
                       index,&(grid.sigmatheta[2*K+1]));
    nc_put_var1_double(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],
                       index,&(grid.sigmatheta[2*K+1]));
    nc_put_var1_double(nc_id,var.v.info[0].coorid[NETCDF_K_INDEX],
                       index,&(grid.sigmatheta[2*K+1]));
    nc_put_var1_double(nc_id,var.hdry.info[0].coorid[NETCDF_K_INDEX],
                       index,&(grid.sigmatheta[2*K]));
    nc_put_var1_double(nc_id,var.p3.info[0].coorid[NETCDF_K_INDEX],
                        index,&(grid.sigmatheta[2*K+1]));
  }

  /*
   * NOTE: The values for the time dimension are written when
   *       the variables are written.
   */

  /*
   * Put back into define mode before returning:
   */
  nc_err = nc_redef(nc_id);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"nc_redef(), %s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  return;
}

/*======================= end of define_netcdf() ==============================*/

/*======================= handle_file_compression() ===========================*/

/*
 * Removes any .gz suffix from file_name and  
 * decompresses the file if necessary.
 */
void handle_file_compression(char *file_name)
{
  char
    *ptr,
    alt_file_name[FILE_STR];
  FILE
    *fp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="handle_file_compression";

  ptr = strrchr(file_name,'.');
  if (ptr && strstr(ptr,".gz")) {
    /*
     * The file_name string has a .gz suffix.
     *
     * Check if the compressed file exists.
     */
    fp = fopen(file_name,"r");
    if (fp) {
      /*
       * The compressed file exists.
       * Decompress it.
       */
      fclose(fp);
      sprintf(Message,"gunzip %s",file_name);
      system(Message);
    }
    /*
     * Remove the .gz suffix from file_name and return.
     */
    *ptr = '\0';
    return;
  }
  else {
    /*
     * The file_name string does not have a .gz suffix.
     *
     * Check if the uncompressed file exists.
     */
    fp = fopen(file_name,"r");
    if (fp) {
      /*
       * The uncompressed file exists, so return.
       */
      fclose(fp);
      return;
    }
    else {
      /*
       * The uncompressed file does not exist.
       * Check if the compressed file exists.
       */
      sprintf(alt_file_name,"%s.gz",file_name);
      fp = fopen(alt_file_name,"r");
      if (fp) {
        /*
         * The compressed file exists.
         * Decompress the file and return.
         */
        sprintf(Message,"gunzip %s",alt_file_name);
        system(Message);
        fclose(fp);
        return;
      }
    }
  }

  return;
}

/*======================= end of handle_file_compression() ====================*/

/*======================= get_jlohi() =========================================*/

/*
 * Determines j range of specified node.
 * The value node = SETUP_GET_JLOHI signals that the data are being passed
 * in via jlo and jhi.
 */

void get_jlohi(int  node,
               int  num_nodes,
               int *jlo,
               int *jhi)
{
  int
    n;
  static int
    initialized = FALSE,
   *jlow,
   *jhigh;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_jlohi";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    jlow  = ivector(0,num_nodes-1,dbmsname);
    jhigh = ivector(0,num_nodes-1,dbmsname);

    initialized = TRUE;
  }

  if (node >= 0 && node < num_nodes) {
    *jlo = jlow[node];
    *jhi = jhigh[node];
  }
  else if (node == SETUP_GET_JLOHI) {
   /*
    * The information is passed in through the jlo, jhi arguments. 
    */
    for (n = 0; n < num_nodes; n++) {
      jlow[n]  = jlo[n];
      jhigh[n] = jhi[n];
    }
  }
  else {
    sprintf(Message,"node=%d,num_nodes=%d is invalid",node,num_nodes);
    epic_error(dbmsname,Message);
  } 

  return;
}

/*======================= end of get_jlohi() ==================================*/

/*======================= prompt_extract_on() =================================*/

#define PRINT_EXTRACT_ON(u,U_INDEX) \
    if (var.extract_on_list[U_INDEX] != -1) { \
      if (var.extract_on_list[U_INDEX] == TRUE) { \
        fprintf(stdout,"     %2d >< %s \n",U_INDEX,var.u.info[0].name); \
      } \
      else if (var.extract_on_list[U_INDEX] == FALSE) { \
        fprintf(stdout,"     %2d    %s \n",U_INDEX,var.u.info[0].name); \
      } \
      var.u.extract_on = var.extract_on_list[U_INDEX] = FALSE; \
    }

void prompt_extract_on(char *def_extract_str)
{
  int
    is,ip,index;
  char
    extract_str[N_STR]="",
    *ptr;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="prompt_extract_on";

  /*
   * Use def_extract_str to turn on extract_on_list.
   */
  var.extract_on = FALSE;
  for (index = FIRST_INDEX; index <= LAST_INDEX; index++) {
    var.extract_on_list[index] = FALSE;
  }

  for (ptr = def_extract_str; strcmp(ptr,"\0") != 0; ptr++) {
    if (isdigit(*ptr)) {
      /* Read in integer. */
      sscanf(ptr,"%d",&index);

      if (index >= FIRST_SPECIES && index <= LAST_SPECIES) {
        if (var.species[is].on != TRUE) {
          /*
           * Remove from def_extract_str.
           */
          sprintf(ptr," ");
        }
        else {
          var.extract_on             = TRUE;
          var.extract_on_list[index] = TRUE;
        }
      }
      else {
        var.extract_on             = TRUE;
        var.extract_on_list[index] = TRUE;
      }

      /* Skip to end of integer just read. */
      while (isdigit(*ptr)) {
        ptr++;
      }
      ptr--;
    }
  }

  /*
   * Print list of variable indices and names, and then reset extract_on to FALSE.
   */
  fprintf(stdout,"\n");
  fprintf(stdout,"   index   variable\n");
  PRINT_EXTRACT_ON(u,        U_INDEX);
  PRINT_EXTRACT_ON(v,        V_INDEX);
  PRINT_EXTRACT_ON(hdry,  HDRY_INDEX);
  PRINT_EXTRACT_ON(theta,THETA_INDEX);
  if (var.nu_turb.on) PRINT_EXTRACT_ON(nu_turb,NU_TURB_INDEX);
  if (var.fpara.on  ) PRINT_EXTRACT_ON(fpara,    FPARA_INDEX);
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      PRINT_EXTRACT_ON(species[is],is);
    }
  }

  PRINT_EXTRACT_ON(hdry3,HDRY3_INDEX);
  PRINT_EXTRACT_ON(pdry3,PDRY3_INDEX);
  PRINT_EXTRACT_ON(p2,P2_INDEX);
  PRINT_EXTRACT_ON(p3,P3_INDEX);
  PRINT_EXTRACT_ON(theta2,THETA2_INDEX);
  PRINT_EXTRACT_ON(h2,H2_INDEX);
  PRINT_EXTRACT_ON(h3,H3_INDEX);
  PRINT_EXTRACT_ON(t2,T2_INDEX);
  PRINT_EXTRACT_ON(t3,T3_INDEX);
  PRINT_EXTRACT_ON(rho2,RHO2_INDEX);
  PRINT_EXTRACT_ON(rho3,RHO3_INDEX);
  PRINT_EXTRACT_ON(exner3,EXNER3_INDEX);

  if (var.fgibb3.on) {PRINT_EXTRACT_ON(fgibb3,FGIBB3_INDEX);}

  PRINT_EXTRACT_ON(gz3,GZ3_INDEX);
  PRINT_EXTRACT_ON(mont3,MONT3_INDEX);
  PRINT_EXTRACT_ON(heat3,HEAT3_INDEX);
  PRINT_EXTRACT_ON(pv3,PV3_INDEX);
  PRINT_EXTRACT_ON(eddy_pv3,EDDY_PV3_INDEX);
  PRINT_EXTRACT_ON(ri2,RI2_INDEX);
  PRINT_EXTRACT_ON(vort3,VORT3_INDEX);
  PRINT_EXTRACT_ON(div_uv3,DIV_UV3_INDEX);
  PRINT_EXTRACT_ON(w3,W3_INDEX);
  PRINT_EXTRACT_ON(dzdt3,DZDT3_INDEX);
  PRINT_EXTRACT_ON(diffusion_coef_uv,DIFFUSION_COEF_UV_INDEX);
  PRINT_EXTRACT_ON(diffusion_coef_theta,DIFFUSION_COEF_THETA_INDEX);
  PRINT_EXTRACT_ON(diffusion_coef_mass,DIFFUSION_COEF_MASS_INDEX);

  if (var.u_spinup.on) {PRINT_EXTRACT_ON(u_spinup,U_SPINUP_INDEX);}

  PRINT_EXTRACT_ON(gz_surface,GZ_SURFACE_INDEX);

  if (grid.nmt_physics_on) {
    PRINT_EXTRACT_ON(dry_entropy,DRY_ENTROPY_INDEX);
    PRINT_EXTRACT_ON(moist_entropy,MOIST_ENTROPY_INDEX);
    PRINT_EXTRACT_ON(sat_moist_entropy,SAT_MOIST_ENTROPY_INDEX);
    PRINT_EXTRACT_ON(the_flux,THE_FLUX_INDEX);
    PRINT_EXTRACT_ON(rt_flux,RT_FLUX_INDEX);
    PRINT_EXTRACT_ON(u_flux,U_FLUX_INDEX);
    PRINT_EXTRACT_ON(v_flux,V_FLUX_INDEX);
    PRINT_EXTRACT_ON(convthrot,CONVTHROT_INDEX);
    PRINT_EXTRACT_ON(fluxthrot,FLUXTHROT_INDEX);
    PRINT_EXTRACT_ON(rain_rate,RAIN_RATE_INDEX);
  }

  fprintf(stdout,"\n");
  sprintf(Message,"On one line, input the indices of variables to be included in extract.nc.\n");
  input_string(Message,def_extract_str,extract_str);

  for (ptr = extract_str; strcmp(ptr,"\0") != 0; ptr++) {
    if (isdigit(*ptr)) {
      /* Read in integer. */
      sscanf(ptr,"%d",&index);

      var.extract_on             = TRUE;
      var.extract_on_list[index] = TRUE;

      /* Turn on extract_on. */
      switch(index) {
        case U_INDEX:
          var.u.extract_on = TRUE;
        break;
        case V_INDEX:
          var.v.extract_on = TRUE;
        break;
        case HDRY_INDEX:
          var.hdry.extract_on = TRUE;
        break;
        case THETA_INDEX:
          var.theta.extract_on = TRUE;
        break;
        case FPARA_INDEX:
          var.fpara.extract_on = TRUE;
        break;
        case NU_TURB_INDEX:
          var.nu_turb.extract_on = TRUE;
        break;
        case HDRY3_INDEX:
          var.hdry3.extract_on = TRUE;
        break;
        case PDRY3_INDEX:
          var.pdry3.extract_on = TRUE;
        break;
        case P2_INDEX:
          var.p2.extract_on = TRUE;
        break;
        case P3_INDEX:
          var.p3.extract_on = TRUE;
        break;
        case THETA2_INDEX:
          var.theta2.extract_on = TRUE;
        break;
        case H2_INDEX:
          var.h2.extract_on = TRUE;
        break;
        case H3_INDEX:
          var.h3.extract_on = TRUE;
        break;
        case T2_INDEX:
          var.t2.extract_on = TRUE;
        break;
        case T3_INDEX:
          var.t3.extract_on = TRUE;
        break;
        case RHO2_INDEX:
          var.rho2.extract_on = TRUE;
        break;
        case RHO3_INDEX:
          var.rho3.extract_on = TRUE;
        break;
        case EXNER3_INDEX:
          var.exner3.extract_on = TRUE;
        break;
        case FGIBB3_INDEX:
          var.fgibb3.extract_on = TRUE;
        break;
        case GZ3_INDEX:
          var.gz3.extract_on = TRUE;
        break;
        case MONT3_INDEX:
          var.mont3.extract_on = TRUE;
        break;
        case HEAT3_INDEX:
          var.heat3.extract_on = TRUE;
        break;
        case PV3_INDEX:
          var.pv3.extract_on = TRUE;
        break;
        case EDDY_PV3_INDEX:
          var.eddy_pv3.extract_on = TRUE;
        break;
        case RI2_INDEX:
          var.ri2.extract_on = TRUE;
        break;
        case VORT3_INDEX:
          var.vort3.extract_on = TRUE;
        break;
        case DIV_UV3_INDEX:
          var.div_uv3.extract_on = TRUE;
        break;
        case W3_INDEX:
          var.w3.extract_on = TRUE;
        break;
        case DZDT3_INDEX:
          var.dzdt3.extract_on = TRUE;
        break;
        case DIFFUSION_COEF_UV_INDEX:
          var.diffusion_coef_uv.extract_on = TRUE;
        break;
        case DIFFUSION_COEF_THETA_INDEX:
          var.diffusion_coef_theta.extract_on = TRUE;
        break;
        case DIFFUSION_COEF_MASS_INDEX:
          var.diffusion_coef_mass.extract_on = TRUE;
        break;
        case U_SPINUP_INDEX:
          var.u_spinup.extract_on = TRUE;
        break;
        case GZ_SURFACE_INDEX:
          var.gz_surface.extract_on = TRUE;
        break;
        /*
         * nmt_physics
         */
        case DRY_ENTROPY_INDEX:
          var.dry_entropy.extract_on = TRUE;
        break;
        case MOIST_ENTROPY_INDEX:
          var.moist_entropy.extract_on = TRUE;
        break;
        case SAT_MOIST_ENTROPY_INDEX:
          var.sat_moist_entropy.extract_on = TRUE;
        break;
        case THE_FLUX_INDEX:
          var.the_flux.extract_on = TRUE;
        break;
        case RT_FLUX_INDEX:
          var.rt_flux.extract_on = TRUE;
        break;
        case U_FLUX_INDEX:
          var.u_flux.extract_on = TRUE;
        break;
        case V_FLUX_INDEX:
          var.v_flux.extract_on = TRUE;
        break;
        case CONVTHROT_INDEX:
          var.convthrot.extract_on = TRUE;
        break;
        case FLUXTHROT_INDEX:
          var.fluxthrot.extract_on = TRUE;
        break;
        case RAIN_RATE_INDEX:
          var.rain_rate.extract_on = TRUE; 
        break;
        default:
          if (index < FIRST_SPECIES || index > LAST_SPECIES) {
            sprintf(Message,"unrecognized index=%d",index);
            epic_error(dbmsname,Message);
          }
          var.species[index].extract_on = TRUE;
          for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
            if (var.species[index].phase[ip].on == 1) {
              var.species[index].phase[ip].extract_on = TRUE;
            }
          }
        break;
      } /* end of switch */

      /* Skip to end of integer just read. */
      while (isdigit(*ptr)) {
        ptr++;
      }
      ptr--;
    }
  }

  return;
}

/*======================= end of prompt_extract_on() ==========================*/

/*======================= prompt_species_on() =================================*/

/*
 * Returns the number of optional species turned on.
 */

#define PRINT_SPECIES_ON(index) \
    if (var.on_list[index] != -1) { \
      if (var.on_list[index] == 1) { \
        fprintf(stdout,"  %2d  ><  %s \n",index,var.species[index].info[0].name); \
      } \
      else if (var.on_list[index] == 0) { \
        fprintf(stdout,"  %2d      %s \n",index,var.species[index].info[0].name); \
      } \
      var.species[index].on = var.on_list[index] = FALSE; \
    }

int prompt_species_on(planetspec *planet,
                      char       *def_species_str)
{
  int
    is,ip,index,
    count = 0;
  char
    species_str[N_STR],
    *ptr;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="prompt_species_on";

  /*
   * Prune optional-variable list as appropriate for each atmosphere.
   * Set on_list to -1 to take off the printed list.
   */
  for (index = FPARA_INDEX; index <= LAST_SPECIES; index++) {
    var.on_list[index]    = -1;
    var.species[index].on = 0;
  }

  if (strcmp(planet->name,"venus") == 0) {
    /*
     * Venus case.
     */
    ;
  }
  else if (strcmp(planet->name,"venus_llr05") == 0) {
    /*
     * Venus_LLR05 benchmark case.
     */
    ;
  }
  else if (strcmp(planet->name,"earth") == 0) {
    /*
     * Earth case.
     */
    var.on_list[H_2O_INDEX] = 0;
    var.on_list[CO_2_INDEX] = 0;
    var.on_list[O_3_INDEX ] = 0;
  }
  else if (strcmp(planet->name,"held_suarez") == 0) {
    /*
     * Held-Suarez benchmark case.
     */
    ;
  }
  else if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * Gas-giant case.
     */
    var.on_list[FPARA_INDEX ] = 0;
    var.on_list[H_2O_INDEX  ] = 0;
    var.on_list[NH_3_INDEX  ] = 0;
    var.on_list[H_2S_INDEX  ] = 0;
    var.on_list[CH_4_INDEX  ] = 0;
    var.on_list[C_2H_2_INDEX] = 0;
    var.on_list[C_2H_6_INDEX] = 0;
    var.on_list[NH_4SH_INDEX] = 0;
  }
  else if (strcmp(planet->name,"titan") == 0) {
    /*
     * Titan case.
     */
    var.on_list[CH_4_INDEX  ] = 0;
    var.on_list[C_2H_2_INDEX] = 0;
    var.on_list[C_2H_6_INDEX] = 0;
  }
  else {
    /*
     * Default to have them all appear on the list.
     */
    for (index = FPARA_INDEX; index <= LAST_SPECIES; index++) {
      var.on_list[index] = 0;
    }
  }

  /*
   * Use def_species_str to turn on var.on_list[index].
   */
  for (ptr = def_species_str; strcmp(ptr,"\0") != 0; ptr++) {
    if (isdigit(*ptr)) {
      /* Read in integer. */
      sscanf(ptr,"%d",&index);

      if (index == FPARA_INDEX) {
        var.fpara.on = var.on_list[index] = TRUE;
      }
      else if (index >= FIRST_SPECIES && index <= LAST_SPECIES) {
        var.species[index].on = var.on_list[index] = TRUE;
      }

      /* Skip to end of integer just read. */
      while (isdigit(*ptr)) {
        ptr++;
      }
      ptr--;
    }
  }

  /*
   * Print list of variable indices and names, and then reset to FALSE.
   */
  fprintf(stdout," Optional prognostic variables:\n");

  if (var.fpara.on == 1) {
    fprintf(stdout,"  %2d  ><  %s \n",FPARA_INDEX,var.fpara.info[0].name);
  }
  else if (planet->x_h2 > 0.) {
    fprintf(stdout,"  %2d      %s \n",FPARA_INDEX,var.fpara.info[0].name);
  }
  var.fpara.on = var.on_list[FPARA_INDEX] = FALSE;

  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    PRINT_SPECIES_ON(is);
  }

  fprintf(stdout,"\n");
  sprintf(Message,"On one line, input the indices of optional prognostic variables to be turned on [none => optional variables off].\n");
  input_string(Message,def_species_str,species_str);

  if (strcmp(species_str,"\n") == 0) {
    return count;
  }

  for (ptr = species_str; strcmp(ptr,"\0") != 0; ptr++) {
    if (isdigit(*ptr)) {
      /* Read in integer. */
      sscanf(ptr,"%d",&index);

      /* Turn on species. */
      if (index == FPARA_INDEX) {
        var.fpara.on = var.on_list[FPARA_INDEX] = TRUE;
        count++;
      }
      else if (index < FIRST_SPECIES || index > LAST_SPECIES) {
        sprintf(Message,"unrecognized index=%d",index);
        epic_error(dbmsname,Message);
      }
      else {
        var.species[index].on = var.on_list[index] = TRUE;
        count++;
      }

      /* Skip to end of integer just read. */
      while (isdigit(*ptr)) {
        ptr++;
      }
      ptr--;
    }
  }
  return count;
}

/*======================= end of prompt_species_on() ==========================*/

/*======================= bcast_char() ========================================*/

void bcast_char(int   node,
                char *str,
                int   num)
{

#if defined(EPIC_MPI)
  MPI_Bcast(str,num,MPI_CHAR,node,para.comm);
#endif


  return;
}

/*======================= end of bcast_char() =================================*/

/*======================= bcast_int() =========================================*/

void bcast_int(int  node,
               int *val,
               int  num)
{

#if defined(EPIC_MPI)
  MPI_Bcast(val,num,MPI_INT,node,para.comm);
#endif

  return;
}

/*======================= end of bcast_int() ==================================*/

/*======================= bcast_float() =======================================*/

void bcast_float(int         node,
                 EPIC_FLOAT *val,
                 int         num)
{

#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif

#if defined(EPIC_MPI)
  MPI_Bcast(val,num,float_type,node,para.comm);
#endif

  return;
}

/*======================= end of bcast_float() ================================*/

/*======================= bcast_double() ======================================*/

void bcast_double(int     node,
                  double *val,
                  int     num)
{

#if defined(EPIC_MPI)
  MPI_Bcast(val,num,MPI_DOUBLE,node,para.comm);
#endif

  return;
}

/*======================= end of bcast_double() ===============================*/

/*====================== read_spacing_file() ==================================*/
/*
 * Read layer spacing data from the file def->layer_spacing_dat.
 * The following lines without an asterisk illustrate the format of this file;
 * one can cut cut and paste these lines to make a template.
   nk         10
   K        p[hPa]
   0          1.0
   1          2.6
   2          6.4
   3         14.2
   4         35.1
   5        118.1
   6        375.6
   7       2416.9
   8       3624.2 
   9       4344.0
  10       6000.0
 */   
void read_spacing_file(init_defaultspec *def,
                       int               mode)
{
  register int
    K,kk;
  char
    token[32],
   *filename;
  EPIC_FLOAT
   *p_ref;
  FILE
    *spacing_file;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_spacing_file";

  filename = def->layer_spacing_dat;

  spacing_file = fopen(filename,"r");
  if (!spacing_file) {
    sprintf(Message,"cannot open %s",filename);
    epic_error(dbmsname,Message); 
  }

  fscanf(spacing_file,"%s %d",token,&def->nk);
  grid.nk = def->nk;

  /*
   * Allocate memory.
   */
  p_ref = fvector(0,2*(grid.nk+1),dbmsname);

  fscanf(spacing_file,"%s %s",token,token);

#if EPIC_PRECISION == DOUBLE_PRECISION
  for (K = 0; K <= grid.nk; K++) {
    fscanf(spacing_file,"%s %lf",token,p_ref+(2*K+1));
  }
#else
  for (K = 0; K <= grid.nk; K++) {
    fscanf(spacing_file,"%s %f",token,p_ref+(2*K+1));
  }
#endif

  /*
   * Fill in layer values.
   */
  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;
    p_ref[kk] = onto_kk(planet,P2_INDEX,p_ref[kk-1],p_ref[kk+1],kk,JLO,ILO);
  }
  p_ref[      0] = p_ref[1]*p_ref[1]/p_ref[2];
  p_ref[2*KHI+2] = p_ref[2*KHI+1]*p_ref[2*KHI+1]/p_ref[2*KHI];

  /*
   * Convert from hPa to Pa.
   */
  for (kk = 0; kk <= 2*(grid.nk+1); kk++) {
    p_ref[kk] *= 100.;
  }

  def->ptop = p_ref[1];
  def->pbot = grid.pbot = p_ref[2*grid.nk+1];

  if (mode == ALL_DATA) {
    for (kk = 0; kk <= 2*(grid.nk+1); kk++) {
      grid.p_ref[kk] = p_ref[kk];
    }
  }

  fclose(spacing_file);
  /*
   * Free allocated memory.
   */
  free_fvector(p_ref,0,2*(grid.nk+1),dbmsname);

  return;
}

/*====================== end of read_spacing_file() ===========================*/

/*======================= read_t_vs_p() =======================================*/

/*
 * NOTE: This function should only be called during initialization, not when the 
 *       model is running, and hence does not need to be MPI ready.
 */

int read_t_vs_p(planetspec *planet,
                int         portion)
{
  char
    infile[FILE_STR],
    header[N_STR];
  int
    nn,ntp;
  EPIC_FLOAT
    p1,t1,dt1;
  FILE
    *t_vs_p;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_t_vs_p";


  /* Look in local directory first.*/
  sprintf(infile,"./t_vs_p.%s",planet->name);
  t_vs_p = fopen(infile,"r");
  if (!t_vs_p) {
    sprintf(infile,EPIC4_PATH"/data/%s/t_vs_p.%s",planet->name,planet->name);
    t_vs_p = fopen(infile,"r");
    if (!t_vs_p) {
      sprintf(Message,"Cannot open file %s",infile);
      epic_error(dbmsname,Message);
    }
  }

  /* Skip over 6-line header: */
  for (nn = 0; nn < 6; nn++) {
    fgets(header,N_STR,t_vs_p);
  }
  /* input number of data points */
  fscanf(t_vs_p,"%d",&ntp); 
  if (portion == SIZE_DATA) {
    fclose(t_vs_p);
    return ntp;
  }

  /* 
   * Store in order of increasing sigmatheta.
   */
  for (nn = ntp-1; nn >= 0; nn--) { 

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(t_vs_p,"%lf %lf %lf",&p1,&t1,&dt1);
#else
    fscanf(t_vs_p,"%f %f %f",&p1,&t1,&dt1);
#endif

    /* convert from hPa to Pa */
    var.pdat[ nn] = 100.*p1; 
    var.tdat[ nn] = t1;
    var.dtdat[nn] = dt1;
  }

  fclose(t_vs_p);

  return ntp;
}

/*======================= end of read_t_vs_p() ================================*/

/*======================= read_meridional_plane() =============================*/

/*
 * Reads data from a text file arrayed in a pressure-latitude meridional plane.
 * An example of the format is the file epic/data/jupiter/temps.jupiter.HCS.
 *
 * The input argument infile is the name of the file, including its path, 
 * and portion is SIZE_DATA or POST_SIZE_DATA.
 *
 * The output arguments np and nlat are the number of pressure and latitude levels
 * in the requested data set.  
 * The output arguments p and lat are the respective 1D dimension values, where latitude
 * increases with index j and pressure decreases with index k.
 * The output argument value is the 2D array of data values associated with varname,
 * such that value[k][j] refers to pressure level k and latitude j.
 */

void read_meridional_plane(planetspec *planet,
                           char       *infile,
                           int         portion,
                           int        *np,
                           int        *nlat,
                           EPIC_FLOAT *p,
                           EPIC_FLOAT *lat,
                           EPIC_FLOAT *value)
{
  int
    k,j;
  char
    header[N_STR];
  FILE
    *input;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_meridional_plane";

  input = fopen(infile,"r");
  if (!input) {
    sprintf(Message,"error opening %s",infile);
    epic_error(dbmsname,Message);
  }
  /* Skip over 6-line header: */
  for (k = 0; k < 6; k++) {
    fgets(header,N_STR,input);
  }
  /* input number of data points */
  fscanf(input,"%d %d",np,nlat);
  if (portion == SIZE_DATA) {
    fclose(input);
    return;
  }

  for (j = 0; j < *nlat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(input,"%lf",lat+j);
#else
    fscanf(input,"%f", lat+j);
#endif

  }
  for (k = 0; k < *np; k++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(input,"%lf",p+k);
#else
    fscanf(input,"%f", p+k);
#endif

    /* Convert hPa to Pa. */
    p[k] *= 100.;

    for (j = 0; j < *nlat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(input,"%lf",value+(j+k*(*nlat)));
#else
    fscanf(input,"%f", value+(j+k*(*nlat)));
#endif

    }
  }

  fclose(input);
  return;
}

/*======================= end of read_meridional_plane() ======================*/
 
/*======================= input_float() =======================================*/

/* 
 * Read in floating-point datum, or set to default if input is a return ('\n').
 *
 * C.Santori, T.Dowling 
 */

EPIC_FLOAT input_float(char       *prompt,
                       EPIC_FLOAT  def) 
{
  char  
    c,
    buffer[N_STR];
  int   
    len;
  EPIC_FLOAT 
    ans;

#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif

  if (IAMNODE == NODE0) {
    fprintf(stdout,"%s[%g]: ",prompt,def);
    for (len = 0; (c = getchar()) != '\n' && len < N_STR; len++) {
      buffer[len]=c;
    }
    buffer[len] = '\0';
    if (len == 0) {
      ans = def;
    }
    else {

#if EPIC_PRECISION == DOUBLE_PRECISION
      sscanf(buffer,"%lf",&ans);
#else
      sscanf(buffer,"%f",&ans);
#endif

    }
  }

#if defined(EPIC_MPI)
   MPI_Bcast(&ans,1,float_type,NODE0,para.comm);
#endif

  return ans;
}

/*====================== end input_float() ====================================*/

/*====================== input_int() ========================================*/

/* 
 * Read in int, or set to default if input is a return ('\n').
 * C.Santori, T.Dowling 
 */

int input_int(char *prompt,
              int   def) 
{
  char  
    c,
    buffer[N_STR];
  int 
    ans,
    len;

  if (IAMNODE == NODE0) {
    fprintf(stdout,"%s[%d]: ",prompt,def);
    for (len = 0; (c = getchar()) != '\n' && len < N_STR; len++) {
      buffer[len]=c;
    }
    buffer[len] = '\0';
    if (len == 0) {
      ans = def;
    }
    else {
      sscanf(buffer,"%d",&ans);
    }
  }

#if defined(EPIC_MPI)
  MPI_Bcast(&ans,1,MPI_INT,NODE0,para.comm);
#endif

  return ans;
}

/*====================== end input_int() ====================================*/

/*====================== input_string() =====================================*/

/* 
 * Read in a string, or set to default if input is a return ('\n').
 *
 * C.Santori, T.Dowling 
 */

void input_string(char *prompt, 
                  char *def, 
                  char *ans) 
{
  char 
    c,
    buffer[N_STR];
  int  
    len;

  if (IAMNODE == NODE0) {
    fprintf(stdout,"%s[%s]: ",prompt,def);
    for (len = 0; (c = getchar()) != '\n' && len < N_STR; len++) {
      buffer[len]=c;
    }
    buffer[len] = '\0';
    if (len == 0) {
      strcpy(ans,def);
    }
    else {
      strcpy(ans,buffer);
      strcpy(def,buffer);
    }
  }

#if defined(EPIC_MPI)
  MPI_Bcast(ans,N_STR,MPI_CHAR,NODE0,para.comm);
#endif

}

/*====================== end input_string() ==================================*/

/*====================== print_model_description() ===========================*/
/*
 * Print to stdout a select listing of model parameters.
 */

void print_model_description(planetspec *planet)
{
  int
    is,index,
    dt_cfl;
  char
    *ptr;
  double  
    max_nu[MAX_NU_ORDER+1];  /* NOTE: declared as double, not EPIC_FLOAT */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="print_model_description";

  /* Calculate CFL dt: */
  dt_cfl = cfl_dt(planet);

  /*
   * Compute hyperviscosity coefficients.
   */
  set_max_nu(max_nu);

  fprintf(stdout,"\n");
  fprintf(stdout,"                     Version: %4.2f\n",grid.data_version);
#if EPIC_PRECISION == DOUBLE_PRECISION
  fprintf(stdout,"              Floating point: double precision\n");
#else
  fprintf(stdout,"              Floating point: single precision\n");
#endif
  /* 
   * Capitalize for printing.
   */
  ptr  = planet->name;
  *ptr = (char)toupper((int)*ptr);
  /*
   * If last character is a digit, capitalize all characters after first '_'. 
   * Otherwise, capitalize first character
   * after each '_', which handles cases like "Held_Suarez".
   */
  ptr = strchr(planet->name,'\0');
  ptr--;
  if (isdigit((int)*ptr)) {
    ptr = strchr(planet->name,'_');
    while (*ptr) {
      ptr++;
      *ptr = (char)toupper((int)*ptr);
    }
  }
  else {
    ptr = planet->name;
    while (*ptr) {
      if (*ptr == '_') {
        ptr++;
        *ptr = (char)toupper((int)*ptr);
      }
      ptr++;
    }
  }
  fprintf(stdout,"                      System: %s \n",planet->name);

  /* 
   * Return to lower case.
   */
  ptr = planet->name;
  while (*ptr) {
    *ptr = (char)tolower((int)*ptr);
    ptr++;
  }

  fprintf(stdout,"                    Geometry: %s \n",grid.geometry);
  fprintf(stdout,"                   Lon range: %-.1f to %-.1f deg\n",
         grid.lon[2*1],grid.lon[2*(grid.ni+1)]);
  fprintf(stdout,"                   Lat range: %-.1f to %-.1f deg \n",
         grid.lat[2*(grid.jlo)],grid.lat[2*(grid.nj+1)]);
  fprintf(stdout,"           Momentum timestep: %s\n",
                 grid.uv_timestep_scheme);
  fprintf(stdout,"              Mass advection: %s \n",
                 var.hdry.advection_scheme);
  fprintf(stdout,"             Theta advection: %s \n",
                 var.theta.advection_scheme);
  fprintf(stdout,"           Turbulence scheme: %s \n",
                 grid.turbulence_scheme);
  if (strcmp(grid.turbulence_scheme,"off") != 0) {
    fprintf(stdout,"   Ri-based stability factor: %s \n",
                   grid.stability_factor);
  }
  strftime(Message,N_STR,"days since %Y-%m-%d %H:%M:%S (UTC)",gmtime(&var.start_time));
  fprintf(stdout,"                        Time: %g %s \n",TIME/86400.,Message);

  if (fcmp(L_s,0.) == 0 || fcmp(L_s,360.) == 0) {
    fprintf(stdout,"                         L_s: %.1f deg (vernal equinox)\n",L_s);
  }
  else if (fcmp(L_s,90.) < 0) {
    fprintf(stdout,"                         L_s: %.1f deg (northern spring)\n",L_s);
  }
  else if (fcmp(L_s,90.) == 0) {
    fprintf(stdout,"                         L_s: %.1f deg (summer solstice)\n",L_s);
  }
  else if (fcmp(L_s,180.) < 0) {
    fprintf(stdout,"                         L_s: %.1f deg (northern summer)\n",L_s);
  }
  else if (fcmp(L_s,180.) == 0) {
    fprintf(stdout,"                         L_s: %.1f deg (autumnal equinox)\n",L_s);
  }
  else if (fcmp(L_s,270.) < 0) {
    fprintf(stdout,"                         L_s: %.1f deg (northern autumn)\n",L_s);
  }
  else if (fcmp(L_s,270.) == 0) {
    fprintf(stdout,"                         L_s: %.1f deg (winter solstice)\n",L_s);
  }
  else {
    fprintf(stdout,"                         L_s: %.1f deg (northern winter)\n",L_s);
  }

  /*
   * This CFL estimate is not currently a reliable guide for a 
   * numerically stable timestep, so we do not print it.
   */
  /****
  fprintf(stdout,"                    Timestep: dt = %d s, CFL dt = %d s\n",
                  grid.dt,dt_cfl);
   ****/
  fprintf(stdout,"                    Timestep: dt = %d s\n",grid.dt);


  if (grid.ni == 1) {
    if (grid.nj == grid.jlo) {
      /* 1D model. */
      fprintf(stdout,"                        Size: nk = %d, 1D \n",
              grid.nk);
    }
    else {
      /* 2D model. */
      fprintf(stdout,"                        Size: nk = %d, nj = %d, 2D \n",
              grid.nk,grid.nj);
    }
  }
  else {
    if (grid.nj == grid.jlo) {
      /* 2D model. */
      fprintf(stdout,"                        Size: nk = %d, ni = %d, 2D \n",
              grid.nk,grid.ni);
    }
    else {
      /* 3D model. */
      fprintf(stdout,"                        Size: nk = %d, nj = %d, ni = %d \n",
              grid.nk,grid.nj,grid.ni);
    }
  }
  /* 
   * Print species list: 
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      fprintf(stdout,"%28s: included \n",var.species[is].info[0].name);
    }
  }
  /*
   * Print chemical information:
   */
  if (var.fpara.on) {
    fprintf(stdout,"        para conversion time: t_p(1bar) = %8.1e s \n",
            var.time_fp_bar/1.e+5);
  }

  if (grid.newt_cool_on == TRUE) {
    if (grid.newt_cool_adjust) {
      fprintf(stdout,"           Newtonian cooling: on (layer average adjusted to zero)\n");
    }
    else {
      fprintf(stdout,"           Newtonian cooling: on (no adjustment to layer average)\n");
    }
    if (grid.prandtl == 0.) {
      fprintf(stdout,"               Rayleigh drag: off \n");
    } 
    else {
      if (grid.drag_v) {
        fprintf(stdout,"               Rayleigh drag: on, tau_rad/tau_drag = %4.2f \n",grid.prandtl);
      }
      else {
        fprintf(stdout,"               Rayleigh drag: on (v drag off), tau_rad/tau_drag = %4.2f \n",grid.prandtl);
      }
      if (grid.drag_zonal_avg) {
        fprintf(stdout," Explicit Rayleigh drag term: zonal avg of (u,v) \n");
      }
      else {
        fprintf(stdout," Explicit Rayleigh drag term: (u,v) \n");
      }
    }
  }
  else {
    fprintf(stdout,"           Newtonian cooling: off \n");
    if (grid.tau_drag >= 1.e+20) {
      fprintf(stdout,"               Rayleigh drag: off \n");
    }
    else {
      if (grid.drag_v) {
        fprintf(stdout,"               Rayleigh drag: on, tau_drag = %8.2e s \n",
                       grid.tau_drag);

      }
      else {
        fprintf(stdout,"               Rayleigh drag: on (v drag off), tau_drag = %8.2e s \n",
                       grid.tau_drag);
      }
      if (grid.drag_zonal_avg) {
        fprintf(stdout," Explicit Rayleigh drag term: zonal avg of (u,v) \n");
      }
      else {
        fprintf(stdout," Explicit Rayleigh drag term: (u,v) \n");
      }
    }
  }


  if (strcmp(planet->name,"earth") == 0) {
    if (!grid.nmt_physics_on) {
      fprintf(stdout,"                 NMT physics: off\n");
    }
    else {
      fprintf(stdout,"                 NMT physics: on\n");
    }
  }
  else {
    if (!grid.microphysics_on) {
      fprintf(stdout,"                Microphysics: off\n");
    }
    else {
      fprintf(stdout,"                Microphysics: on\n");
    }
  }

  if (grid.diffusion_direction == HORIZONTAL_AND_VERTICAL) {
    fprintf(stdout,"        Horizontal diffusion: on\n");
    fprintf(stdout,"          Vertical diffusion: on\n");
  }
  else if (grid.diffusion_direction == JUST_HORIZONTAL) {
    fprintf(stdout,"        Horizontal diffusion: on\n");
    fprintf(stdout,"          Vertical diffusion: off\n");
  }
  else if (grid.diffusion_direction == JUST_VERTICAL) {
    fprintf(stdout,"        Horizontal diffusion: off\n");
    fprintf(stdout,"          Vertical diffusion: on\n");
  }
  else {
    fprintf(stdout,"        Horizontal diffusion: off\n");
    fprintf(stdout,"          Vertical diffusion: off\n");
  }

  if (grid.nudiv_nondim > 0.) {
    fprintf(stdout,"          Divergence damping: %3.2g (%9.3e m^2/s)\n",
                    grid.nudiv_nondim,grid.nudiv_nondim*max_nu[2]);
  }
  else {
    fprintf(stdout,"          Divergence damping: off \n");
  }
  for (index = 4; index <= MAX_NU_ORDER; index += 2) {
    if (grid.nu[index] == 0.) {
      fprintf(stdout,"                       nu[%d]: off \n",index);
    }
    else {
      fprintf(stdout,"                       nu[%d]: %.3g (%9.3e m^%d/s)\n",
              index,grid.nu[index]/max_nu[index],grid.nu[index],index);
    }
  }
    
  if (grid.k_sponge == 0) {
    fprintf(stdout,"                    k_sponge: off\n");
  }
  else {
    fprintf(stdout,"                    k_sponge: %d\n",grid.k_sponge);
  }

  fprintf(stdout,"           Hybrid coordinate: k_sigma = %d, sigma_sigma = %.3f\n",grid.k_sigma,grid.sigma_sigma);

  if (strcmp(planet->type,"terrestrial") == 0) {
    if (grid.topo_scale == 1.) {
      fprintf(stdout,"           Topo_scale factor: off \n");
    }
    else {
      fprintf(stdout,"           Topo_scale factor: %-4.1f \n",grid.topo_scale);
    }
  }

  return;
}

/*====================== end print_model_description() =======================*/


/*====================== print_zonal_info() ==================================*/

void print_zonal_info(planetspec *planet)
{
  int
    K,J,I;
  EPIC_FLOAT
    zeta,zf,zfy,fy,
    *buffji;
  FILE
    *u_dat;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="print_zonal_info";

  if (IAMNODE != NODE0) {
    /* NOTE: not set up for MPI. */
    sprintf(Message,"not set up for MPI");
    epic_error(dbmsname,Message);
  }

  /*
   * Allocate memory.
   */
  buffji = fvector(0,Nelem2d-1,dbmsname);

  K = grid.nk;
  I = ILO;

  /*
   * Calculate absolute vorticity on pv-grid.
   */
  vorticity(planet,ON_SIGMATHETA,ABSOLUTE,
            var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
            var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
            NULL,
            buffji);

  u_dat = fopen("zonal_wind.dat","w");
  fprintf(u_dat," %s,  nk, nj, ni = %2d, %2d, %2d \n",
                 planet->name,grid.nk,grid.nj,grid.ni);
  fprintf(u_dat," zonal wind data at K,I = %2d, %2d \n",K,I);
  fprintf(u_dat,"      f,zeta,z+f units: 1.e-4 s^-1\n");
  fprintf(u_dat," df/dy,d(z+f)/dy units: 1.e-12 m^-1 s^-1 \n\n");
  fprintf(u_dat," lat(deg)  u(m/s)      f       zeta"
                "       zeta+f     df/dy     d(zeta+f)/dy \n");

  for (J = JHI; J >= JLO; J--) {
    /* Average vorticity onto u-grid. */
    zeta = .5*(BUFFJI(J,I)-grid.f[2*J]+BUFFJI(J+1,I)-grid.f[2*(J+1)]);
    zf   = .5*(BUFFJI(J,I)+BUFFJI(J+1,I));
    zfy  = grid.n[2*J+1]*(BUFFJI(J+1,I)-BUFFJI(J,I));
    fy   = grid.n[2*J+1]*(grid.f[2*J+2]-grid.f[2*J]);
    /* change units on zeta terms */
    /* change units on beta=df/dy terms */
    zeta *= 1.e+4;
    zf   *= 1.e+4;
    zfy  *= 1.e+12;
    fy   *= 1.e+12;
    fprintf(u_dat," %5.1f   %7.2f  %10.3e %10.3e %10.3e %10.3e %10.3e \n",
                  grid.lat[2*J+1],U(grid.it_uv,grid.nk,J,I),
                  grid.f[2*J+1]*1.e+4,zeta,zf,fy,zfy);
  }

  fclose(u_dat);

  /*
   * Free allocated memory.
   */
  free_fvector(buffji,0,Nelem2d-1,dbmsname);

  return;
}

/*====================== end of print_zonal_info() ===========================*/

/*====================== print_vertical_column() =============================*/

void print_vertical_column(planetspec *planet,
                           int         J,
                           int         I,
                           char       *filename)
{
  int
    K;
  EPIC_FLOAT
    pressure,theta,temperature,
    brunt2;
  FILE
    *vert_dat;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="print_vertical_column";

  if (IAMNODE != 0) {
    return;
  }

  vert_dat = fopen(filename,"w");
  fprintf(stdout,"\n Column data at lat = %.1f, lon = %.1f: \n",
                 grid.lat[2*J+1],grid.lon[2*I+1]);
  fprintf(stdout,"\n     K       sgth[K]     press[hPa]      theta[K]  N2[1/s^2]  \n");
  fprintf(vert_dat,"  Vertical profile for lat=%.1f lon=%.1f \n",
                    grid.lat[2*J+1],grid.lon[2*I+1]);
  fprintf(vert_dat,"   K     sgth[K]  press[hPa]      temp[K]  theta[K]  N2[1/s^2]  U[m/s]\n");
  K = 0;
  fprintf(stdout,"-- %4.1f -- %9.1f -- %12.5g -- %9.1f ----------- \n",
         (EPIC_FLOAT)K+0.5,grid.sigmatheta[2*K+1],P3(K,J,I)/100.,THETA(K,J,I));
  brunt2 = get_brunt2(planet,2*K+1,J,I);
  fprintf(vert_dat," %4.1f  %9.1f %13.6e  %9.1f %9.1f  %9.6f %7.2f\n",
                   (EPIC_FLOAT)K+.5,grid.sigmatheta[2*K+1],P3(K,J,I)/100.,T3(K,J,I),THETA(K,J,I),brunt2,U(grid.it_uv,K,J,I));
  for (K = 1; K <= grid.nk; K++) {
    brunt2 = get_brunt2(planet,2*K,J,I);
    fprintf(stdout,"   %4.1f    %9.1f                    %9.1f  %9.6f \n", 
                   (EPIC_FLOAT)K,grid.sigmatheta[2*K],THETA2(K,J,I),brunt2);

    brunt2 = get_brunt2(planet,2*K+1,J,I);
    fprintf(vert_dat," %4.1f  %9.1f %13.6e  %9.1f %9.1f  %9.6f %7.2f\n",
                     (EPIC_FLOAT)K+.5,grid.sigmatheta[2*K+1],P3(K,J,I)/100.,T3(K,J,I),THETA(K,J,I),brunt2,U(grid.it_uv,K,J,I));
    fprintf(stdout,"-- %4.1f -- %9.1f -- %12.5g -- %9.1f ----------- \n",
           (EPIC_FLOAT)K+0.5,grid.sigmatheta[2*K+1],P3(K,J,I)/100.,THETA(K,J,I));
  }

  fclose(vert_dat);

  return;
}

/*====================== end of print_vertical_column() ======================*/

/*====================== node0_barrier_print() ===============================*/

/*
 * Aaron Herrnstein.
 *
 * Routine prints a message using Processor 0.  MPI_Barrier() calls are placed
 * before and after the print statement.  This can be beneficial when tracking
 * down MPI bugs.
 */
void node0_barrier_print(char *calling_function,
                         char *Message, 
                         int  display_time)
{
  time_t rawtime;

#if defined(EPIC_MPI)
  MPI_Barrier(para.comm); 
#endif

  if (IAMNODE==NODE0) {
    /*
     * Print Message using processor 0
     */
    if (display_time) {
      time( &rawtime );
      fprintf(stdout, "%s(),  %s  at %s\n", calling_function, Message, ctime(&rawtime) );
    }
    else {
      fprintf(stdout, "%s(),  %s\n", calling_function, Message );
    }
    fflush(stdout);
  }

#if defined(EPIC_MPI)
  /* The barrier looks redundant here, but is necessary for proper io. */
  MPI_Barrier(para.comm);
#endif

}

/*====================== end of node0_barrier_print() ========================*/

/*====================== scdswap() ===========================================*/

/*
 *  Byte swapping for 8, 4 and 2 byte quantities.  
 */
#include <stdio.h>

void scdswap(char *arr, 
             int   nbytes, 
             int   cnt)
{
  char 
    buf[4];
  register int 
    nb=nbytes;
  register char 
   *parr, *pend;

  pend = arr+nb*cnt;

  switch(nb) {
    case 1:
    break;
    case 2:  
      for (parr=arr; parr<pend; parr+=nb) {
        buf[0]  = parr[0];
        parr[0] = parr[1];
        parr[1] = buf[0];
      }
    break;
    case 4:  
      for (parr=arr; parr<pend; parr+=nb) {
        buf[0]  = parr[0];
        buf[1]  = parr[1];
        parr[0] = parr[3];
        parr[1] = parr[2];
        parr[2] = buf[1];
        parr[3] = buf[0];
      }
    break;
    case 8:  
      for (parr=arr; parr<pend; parr+=nb) {
        buf[0]  = parr[0];
        buf[1]  = parr[1];
        buf[2]  = parr[2];
        buf[3]  = parr[3];
        parr[0] = parr[7];
        parr[1] = parr[6];
        parr[2] = parr[5];
        parr[3] = parr[4];
        parr[4] = buf[3];
        parr[5] = buf[2];
        parr[6] = buf[1];
        parr[7] = buf[0];
      }
    break;
    default: 
      fprintf(stderr," Bad length to scdswap()\n");
      exit(99);
    break;
  } 

  return;
} 

/*======================= end of scdswap() ===================================*/

/*======================= epic_error() =======================================*/
/*
 * Prints machine name, calling node rank, timestep, calling function name, 
 * and Message to stderr, then aborts.
 */
void epic_error(char *calling_function,
                char *Message)
{
  fprintf(stderr,"\n** EPIC Error: node=%d, timestep=%lu, %s(): %s\n",
                 IAMNODE,grid.itime,calling_function,Message);
  fflush(stderr);

#if defined(EPIC_MPI)
  MPI_Abort(MPI_COMM_WORLD,1);
#endif

  exit(1);
}

/*======================= end of epic_error() ================================*/

/*======================= epic_warning() =====================================*/
/*
 * Prints machine name, calling node rank, calling function name, 
 * and Message to stderr (does not exit).
 */
void epic_warning(char *calling_function,
                  char *Message)
{
  fprintf(stderr,"EPIC Warning: node=%d, timestep=%lu, %s(): %s\n",
                 IAMNODE,grid.itime,calling_function,Message);
  fflush(stderr);

  return;
}

/*======================= end of epic_warning() ==============================*/


/*======================= declare_copyright() ================================*/

void declare_copyright(void)
{
  static int
    initialized=FALSE;

  if (!initialized) {
    fprintf(stderr,"\n");
    fprintf(stderr," EPIC Model, Copyright (C) 1998-2009 Timothy E. Dowling \n");                                                                                         
    fprintf(stderr," This program is free software; you can redistribute it and/or \n");  
    fprintf(stderr," modify it under the terms of the GNU General Public License.  \n");    
    fprintf(stderr," This program is distributed WITHOUT ANY WARRANTY.             \n"); 
    fprintf(stderr,"\n");

    initialized = TRUE;
  } 
                                                                 
  return;
}

/*======================= end of declare_copyright() =========================*/

/* * * * * * * * * * * * *  end of epic_funcs_io.c  * * * * * * * * * * * * * */




