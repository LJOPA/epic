#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# *                                                                 *
# * Copyright (C) 1998-2009 Timothy E. Dowling                      *
# *                                                                 *
# * This program is free software; you can redistribute it and/or   *
# * modify it under the terms of the GNU General Public License     *
# * as published by the Free Software Foundation; either version 2  *
# * of the License, or (at your option) any later version.          *
# * A copy of this License is in the file:                          *
# *   $EPIC4_PATH/License.txt                                       *
# *                                                                 *
# * This program is distributed in the hope that it will be useful, *
# * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
# *                                                                 *
# * You should have received a copy of the GNU General Public       *
# * License along with this program; if not, write to the Free      *
# * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     *
# * Boston, MA 02110-1301, USA.                                     *
# *                                                                 *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#
#  The top Makefile for the EPIC model.
#
#  T. Dowling
#
#  NOTE: Each logical line after a target is executed as a separate
#        shell, therefore, to propagate variables use semicolons and
#        backslashes to string commands together into one line.
#
#  NOTE: Do not put spaces around the "=" in variable assignments.
#

SHELL = /bin/sh

all:
	@if [ $(MACHINE_TYPE) = LINUX ] ; then \
	  export CC=/usr/bin/gcc ; \
	  export FC="" ; \
	  export CXX="" ; \
	elif [ $(MACHINE_TYPE) = Darwin ] ; then \
	  export CC=/usr/bin/cc ; \
	  export FC="" ; \
	  export CXX="" ; \
	fi ; \
	cd $(EPIC4_PATH)/netcdf/src ; \
	if [ -f ./config.status ] ; then \
	  echo "netCDF is built" ; \
	else \
	  echo "Configuring netCDF (Network Common Data Format)..." ; \
	  ./configure --prefix=$(EPIC4_PATH)/netcdf > configure.log ; \
	  echo "Making netCDF..." ; \
	  make ; \
	  echo "Installing netCDF..." ; \
	  make install ; \
	fi ;
	@if [ "$(EPIC_CFLAG)" = "-g" ] ; then \
	  export CFLAGS="-g" ; \
	  export LDFLAGS="-g" ; \
	  export EPIC_CHECK="1" ; \
	elif [ "$(EPIC_CFLAG)" = "-O_opteron" ] ; then \
	  export CFLAGS="-march=opteron -O2 -pipe" ; \
	  export LDFLAGS="" ; \
	  export EPIC_CHECK="0" ; \
	elif [ "$(EPIC_CFLAG)" = "-O_CPL" ] ; then \
	  export CFLAGS="-O3 -march=athlon-xp" ; \
	  export LDFLAGS="" ; \
	  export EPIC_CHECK="0" ; \
	elif [ "$(EPIC_CFLAG)" = "-O_CPL2" ] ; then \
	  export CFLAGS="-O3" ; \
	  export LDFLAGS="" ; \
	  export EPIC_CHECK="0" ; \
	elif [ "$(EPIC_CFLAG)" = "-O_tdmac" ] ; then \
	  export CFLAGS="-O3" ; \
	  export LDFLAGS="" ; \
	  export EPIC_CHECK="0" ; \
	elif [ "$(EPIC_CFLAG)" = "-O" ] ; then \
	  if [ $(EPIC_CC) = gcc ] ; then \
	    export CFLAGS="-g -O2" ; \
	    export LDFLAGS="-g" ; \
	    export EPIC_CHECK="0" ; \
	  elif [ $(MACHINE_TYPE) = hpux ] ; then \
	    export CFLAGS="-O +Onolimit" ; \
	    export EPIC_CHECK="0" ; \
	  elif [ $(MACHINE_TYPE) = sp2 ] ; then \
	    export CFLAGS="-O -qmaxmem=-1" ; \
	    export EPIC_CHECK="0" ; \
	  else \
	    export CFLAGS="-O2" ; \
	    export EPIC_CHECK="0" ; \
	  fi ; \
	elif [ "$(EPIC_CFLAG)" = "-p" ] ; then \
	  if [ $(EPIC_CC) = gcc ] ; then \
	    export CFLAGS="-pg" ; \
	    export LDFLAGS="-pg" ; \
	    export EPIC_CHECK="0" ; \
	  elif [ $(MACHINE_TYPE) = hpux ] ; then \
	    export CFLAGS="-p" ; \
	    export LDFLAGS="-p" ; \
	    export EPIC_CHECK="0" ; \
	  elif [ $(MACHINE_TYPE) = sp2 ] ; then \
	    export CFLAGS="-p" ; \
	    export LDFLAGS="-p" ; \
	    export EPIC_CHECK="0" ; \
	  else \
	    export CFLAGS="-pg" ; \
	    export LDFLAGS="-pg" ; \
	    export EPIC_CHECK="0" ; \
	  fi ; \
	elif [ "$(EPIC_CFLAG)" = "-pO" ] ; then \
	  if [ $(EPIC_CC) = gcc ] ; then \
	    export CFLAGS="-pg -O2" ; \
	    export LDFLAGS="-pg" ; \
	    export EPIC_CHECK="0" ; \
	  elif [ $(MACHINE_TYPE) = hpux ] ; then \
	    export CFLAGS="-p -O +Onolimit" ; \
	    export LDFLAGS="-p" ; \
	    export EPIC_CHECK="0" ; \
	  elif [ $(MACHINE_TYPE) = sp2 ] ; then \
	    export CFLAGS="-p -O -qmaxmem=-1" ; \
	    export LDFLAGS="-p" ; \
	    export EPIC_CHECK="0" ; \
	  else \
	    export CFLAGS="-pg -O2" ; \
	    export LDFLAGS="-pg" ; \
	    export EPIC_CHECK="0" ; \
	  fi ; \
	else \
	  echo "" ; \
	  echo "No recognized EPIC_CFLAG; defaulting to environment-variable compiler and linker flags" ; \
	  echo "Setting EPIC_CHECK=0" ; \
	  echo "" ; \
	  export EPIC_CHECK="0" ; \
	fi ; \
	if [ $(MPI_TYPE) != none ] ; then \
	  cd $(EPIC4_PATH)/src/mpg ; \
	  make ; \
	  cd $(EPIC4_PATH)/src/mpi ; \
	  make ; \
	else \
	  echo MPI_TYPE = $(MPI_TYPE) ; \
	fi ; \
	cd $(EPIC4_PATH)/src/single ; \
	make

clean:
	@cd $(EPIC4_PATH)/netcdf/src ; \
	make clean ; \
	cd $(EPIC4_PATH)/src/mpi ; \
	make clean ; \
	cd $(EPIC4_PATH)/src/single ; \
	make clean ; \
	cd $(EPIC4_PATH)/src/mpg ; \
	make clean

distclean:
	@cd $(EPIC4_PATH)/src ; \
	rm -f *BAK *bak *~ *.a *.log ; \
	cd $(EPIC4_PATH)/netcdf/src ; \
	make distclean ; \
	cd $(EPIC4_PATH)/netcdf ; \
	rm -r -f bin include lib share ; \
	cd $(EPIC4_PATH)/src/mpi ; \
	make clear ; \
	cd $(EPIC4_PATH)/src/single ; \
	make clear ; \
	cd $(EPIC4_PATH)/src/mpg ; \
	make clear ; \
	cd $(EPIC4_PATH) ; \
	rm -f *BAK *bak *~ *% *.log ; \
	cd $(EPIC4_PATH)/tmp ; \
	rm -f *; \
	cd $(EPIC4_PATH)/help ; \
	rm -f *BAK *bak *~ *% *.log ; \
	cd $(EPIC4_PATH)/netcdf/bin ; \
	rm -f *.nc *.gz dump ; \
	cd $(EPIC4_PATH)/tools/idl ; \
	rm -f *BAK *bak *~ *% *.ps *.dat *.nc *.pca *.gif *.log; \
	cd $(EPIC4_PATH)/tools/matlab ; \
	DISTCLEAN_EPICncVis ; \
	cd $(EPIC4_PATH)/util ; \
	rm -f *BAK *bak *~ *% *.log ; \
	cd $(EPIC4_PATH)/include ; \
	rm -f *BAK *bak *~ *% *.log; \
	cd $(EPIC4_PATH)/src/shared ; \
	rm -f core *o *BAK *bak *~ *% *.log ; \
	cd $(EPIC4_PATH)/src/shared/rt ; \
	rm -f core *o *BAK *bak *~ *% *.log ; \
	cd $(EPIC4_PATH)/src/shared/microphysics ; \
	rm -f core *o *BAK *bak *~ *% *.log ; \
	cd $(EPIC4_PATH)/bin ; \
	rm -f *

clear: distclean

install: all
	@if [ $(MPI_TYPE) != none ] ; then \
	  cd $(EPIC4_PATH)/src/mpg ; \
	  make install ; \
	  cd $(EPIC4_PATH)/src/mpi ; \
	  make install; \
	else \
	  echo MPI_TYPE = $(MPI_TYPE) ; \
	fi ; \
	cd $(EPIC4_PATH)/src/single ; \
	make install

depend:
	@if [ $(MPI_TYPE) != none ] ; then \
	  cd $(EPIC4_PATH)/src/mpi ; \
	  make depend ; \
	else \
	  echo MPI_TYPE = $(MPI_TYPE) ; \
	fi ; \
	cd $(EPIC4_PATH)/src/single ; \
	make depend
