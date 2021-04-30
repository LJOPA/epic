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

 Explicit Planetary Isentropic-Coordinate 
      (EPIC) Atmospheric Model
            Version 4.30

  Required: 
    Unix operating system

  Included third-party, GPL-licence open software:
    netCDF: self-describing file format

  Optional: 
    MPI (Message Passing Interface), for running EPIC on parallel
    computers. Free versions of MPI are available, including openmpi
    and mpich, and commercial versions are available.
    
    A host of software packages will input the self-describing
    netcdf files used by EPIC (epic.nc and extract.nc).  They tend
    to have an Earth-continent grid superimposed by default, but this
    can easily be turned off.
    We recommend:
      Panoply:
        free from http://www.giss.nasa.gov/tools/panoply/
      IDV (Integrated Data Viewer): 
        free from http://www.unidata.ucar.edu/software/idv/

    Sample tools for Matlab, IDL, and Mathematica are included in
    $EPIC4_PATH/tools. We welcome additions to this collection. 

 Topics covered:
   A. INSTALLING
   B. TROUBLESHOOTING
   C. INITIALIZING
   D. CHANGING
   E. RUNNING FROM COMMAND-LINE
    
=========================================================================
=========================================================================

A. INSTALLING 

  1.  To decompress and de-archive the file epic.tar.gz:

      a) If a directory named epic already exists, delete it (or move it):
           rm -r epic
      b) Type 
           gunzip epic.tar.gz 
           tar xpf epic.tar
           gzip epic.tar

      This will produce a directory tree similar to:

epic/           Top EPIC-model directory.
  bin/          Executable files.
  data/         Data for the planets.
  help/         Help files for executables.
  include/      Header (.h) files.
  License.txt   Software license.
  netcdf/       Source code for Network Common Data Format.
  notes/        Brief remarks on version history, etc.
  README.txt    This file.
  src/
    mpg/        Grid decomposition MPI functions.
    mpi/        MPI (parallel-processor) specific code.
    shared/     Platform independent code.
    single/     Single-processor specific code.
  tmp/          A location for temporary files.
  tools/        Analysis tools that work with external software.
  util/         Unix utility programs such as shell scripts.

  2. To archive and compress the model:

      a) Clear any object code by typing
           cd $EPIC4_PATH/src
           make clear

      b) Move up to the directory containing epic and type
           touch epic
           tar cf epic.tar epic
           gzip epic.tar

     Note: usage of the environment variable $EPIC4_PATH
     is enabled after completing step 3.


  3. In your shell resource file (such as .cshrc or .bashrc), add the
     following lines, and then edit them for your environment:

     For a .cshrc or .tcshrc file:

setenv EPIC4_PATH        ~/epic
setenv MPI_TYPE         openmpi
setenv MACHINE_TYPE     `$EPIC4_PATH/util/get_machine_type.sh`
setenv EPIC_CC          gcc
setenv EPIC_CFLAG       -g
setenv EPIC_PRECISION   8

      For a .bashrc file:

export EPIC4_PATH=~/epic
export MPI_TYPE=openmpi
export MACHINE_TYPE=`$EPIC4_PATH/util/get_machine_type.sh`
export EPIC_CC=gcc
export EPIC_CFLAG=-g
export EPIC_PRECISION=8

      NOTE: It is important that the quotes around
            $EPIC4_PATH/util/get_machine_type.sh 
            be backwards single quotes on both ends (`...`),
            such that the output of get_machine_type.sh gets written
            into MACHINE_TYPE.

      EPIC4_PATH is the unix path where the EPIC model is kept.
      It is distinct from "EPIC_PATH" without the "4", which should be used for
      the odler Version 3 of the model (the pure isentropic-coordinate version).
      Keeping these distinct facilitates compiling different versions of
      the EPIC model on the same computer.

      MPI_TYPE
            Assuming the syntax for .cshrc, the choices are:

            No MPI: setenv MPI_TYPE none
          Open MPI: setenv MPI_TYPE openmpi
               LAM: setenv MPI_TYPE LAM
             mpich: setenv MPI_TYPE mpich

      MACHINE_TYPE is set by the shell script get_machine_type.sh.

      EPIC_CC is the name of the C compiler to be used.

           NOTE: Use the most ANSI-compatible compiler available.
                 For example, on the SunOS use c89, not the default cc,
                 or install and use the most recent gcc from Gnu.

      EPIC_CFLAG has the following recognized cases
           "-O" for optimizing
           "-g" for debugging and extra checks of the validity of values
           "-p" for profiling
           "-pO" for profiling and optimizing
           In each case, an attempt is made in the top makefile to use
           the best flags depending on the platform. 

           Alternatively, you can elect to not set EPIC_CFLAG ("unsetenv EPIC_CFLAG"),
           or set it to an unrecognized case, and then pass your own flags via the 
           environment variables CFLAGS and LDFLAGS, 
           to the C compiler and linker, respectively. This is easy
           to do and advanced users are welcome to do it.

           NOTE: gcc version 2.96 with -O2 optimization, which is what "EPIC_CFLAG -O" sets up,
                 sometimes yields spurious nans (not-a-numbers) when running EPIC, and should 
                 therefore not be used. No such problem is known to occur with "EPIC_CFLAG -g" 
                 for gcc 2.96.
                 
      EPIC_PRECISION is set to "4" or "8" to compile the model with 
           single-precision or double-precision floating-point variables, 
           respectively. Double precision is recommended for long runs.

    In addition, for MPI you may need to add lines to your path like the following:
      /usr/local/mpi/lib/hpux/ch_shmem   (that is, ../lib/$(ARCH)/$(COMM))
      /usr/local/mpi/bin

    These changes will take effect when you next log in. 
    You can make the changes take effect now by typing:
      source ~/.cshrc

  4.  To compile the model:

      a) Change directories to the source directory by typing
           cd $EPIC4_PATH/src 

      b) If desired, clear previously compiled object code by typing
           make clear
         If only a few changes have been made to the source code, then
         skipping this step will avoid unnecessary compilations.
         To be on the safe side, after modifications to header files like 
         $EPIC4_PATH/include/epic.h, one should use make clear.

      c) Compile the model by typing
           make install 

      If the compilation is successful, the executables will be  
      found in $EPIC4_PATH/bin.

  5.  The IDL tools in $EPIC4_PATH/tools/IDL use the environoment variable
      IDL_EPIC4_PATH, which should be set to be the working directory for
      IDL plots and movies.

=========================================================================
=========================================================================

B. TROUBLESHOOTING

 1. If the netCDF make fails, look at the file 
    $EPIC4_PATH/netcdf/src/INSTALL for examples of environment variables
    that work for various platforms.  Use a working combination for your
    platform (ideally ones with "-O" flags for optimization) in your 
    ~/.cshrc file, type "source ~/.cshrc" and then recompile.

 2. A "semget" or other semaphore error can usually be cleared up by using 
    ipcs and ipcrm, or for example mpich's sbin/cleanipcs script, which 
    calls these functions.

=========================================================================
=========================================================================

C. INITIALIZING

   The information needed to start the EPIC model is contained in
   a single epic.nc file. The program "initial" generates this
   information.  

=========================================================================
=========================================================================

D. CHANGING

   If you wish to change only a few parameters in an existing
   epic.nc file, use the program "change," which reads an epic.nc file,
   prompts for changes, and then writes a new file.
   You will probably want to customize the prompts in 
   $EPIC4_PATH/src/shared/epic_change.c. 

   To learn about other uses of the change program try the 
   "-h" or "-help" command-line flag.

=========================================================================
=========================================================================

E. RUNNING FROM COMMAND-LINE

   Typing the name of the epic program with the command-line flag
   -h or -help prints the options on the screen.  For example, type:

%epic -h

   to get a list of command-line flags and their meanings.
      
   To integrate the file epic.nc forward 10000 timesteps as a batch job,
   saving twice and backing up every 1000, with error messages to be 
   written to the file ``epig.log'', type:

%epic -itrun 10000 -itsave 5000 -itback 1000 epic.nc >& epic.log &

   To write an extract.nc file every 100 steps, first choose the variables
   to extract using initial or change, and then run the model with:

%epic -itrun 10000 -itsave 5000 -itback 1000 -itextract 100 epic.nc >& epic.log &

   To append to an existing extract.nc file called extract2.nc, type

%epic -itrun 10000 -itsave 5000 -itback 1000 -itextract 100 -append extract2.nc epic.nc >& epic.log &

=========================================================================
=========================================================================

