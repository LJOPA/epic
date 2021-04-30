#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# *                                                                 *
# * Copyright (C) 2008-2009 Thomas K. Greathouse                    *
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

RT_HOME     = $(EPIC4_PATH)/src/shared/rt
F77         = g77
LINKER      = g77
AR          = ar crl
RANLIB      = ranlib

PROGS = librt \
        sat_real8_hr

LIBSRCS     = gravsub.f pscale.f tominterplin.f tominterp.f mixrshay.f expint8.f \
              b.f locate.f column.f diurnaltest8.f uv8.f midtomexpb.f \
              fir.f tomsvisabsorp.f nir8_hitran.f callbruno.f bdiurn2.f \
              vis18.f vis28.f vis38.f nir8.f nir8_ir.f tomsfastest_new.f \
              exptablesub.f polin2.f irwin_rewrite.f polint.f irwin.f \
              shadow8.f numtostring5.f 

LIBOBJS     = gravsub.o pscale.o tominterplin.o tominterp.o mixrshay.o expint8.o \
              b.o locate.o column.o diurnaltest8.o uv8.o midtomexpb.o \
              fir.o tomsvisabsorp.o nir8_hitran.o callbruno.o bdiurn2.o \
              vis18.o vis28.o vis38.o nir8.o nir8_ir.o tomsfastest_new.o \
              exptablesub.o polin2.o irwin_rewrite.o polint.o irwin.o \
              shadow8.o numtostring5.o 

SATSRCS    = sat_real8_hr.f

SATOBJS    = sat_real8_hr.o

LBLSRCS    = lbyl_real8_hr.f

LBLOBJS    = lbyl_real8_hr.o

SRCS = $(SATSRCS) $(LBLSRCS)

INCLUDE_DIR =  -I$(RT_HOME) \
               -I$(EPIC4_PATH)/include

LD_FLAGS = $(LDFLAGS) $(XLDFLAGS)

LIB_PATH    = -L$(RT_HOME)

LIB_LIST    = -lrt

LIBS        = $(LIB_PATH) $(LIB_LIST)

.f.o:
	$(F77) $(F_FLAGS) -c $*.f

	     
F_FLAGS	    =  $(FFLAGS) $(INCLUDE_DIR)


all: $(PROGS)

librt: librt.a
	$(RANLIB) librt.a

librt.a: $(LIBOBJS)
	$(AR) librt.a $(LIBOBJS)

sat_real8_hr: $(LIBOBJS) $(SATOBJS)
	$(LINKER) $(LD_FLAGS) $(SATOBJS) $(LIBS) -o $@

lbyl_real8_hr: $(LIBOBJS) $(LBLOBJS)
	$(LINKER) $(LD_FLAGS) $(LBLOBJS) $(LIBS) -o $@

install: all

clear:
	/bin/rm -f *.o *.a *~ *.BAK

clean:
	/bin/rm -f *.o *.a *~ *.BAK

depend:
	makedepend -- $(F_FLAGS) -- $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.
