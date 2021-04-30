#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# *                                                                 *
# * Copyright (C) 1998 Timothy E. Dowling                           *
# *                                                                 *
# * This program is free software; you can redistribute it and/or   *
# * modify it under the terms of the GNU General Public License     *
# * as published by the Free Software Foundation; either version 2  *
# * of the License, or (at your option) any later version.          *
# * A copy of this License is in the file:                          *
# *   $EPIC_PATH/GNU_General_Public_License.                        *
# *                                                                 *
# * This program is distributed in the hope that it will be useful, *
# * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
# *                                                                 *
# * You should have received a copy of the GNU General Public       *
# * License along with this program; if not, write to the Free      *
# * Software Foundation, Inc., 59 Temple Place - Suite 330,         *
# * Boston, MA  02111-1307, USA.                                    *
# *                                                                 *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#
# Makefile for $EPIC_PATH/data/jupiter/gen_t_vs_p.c.
#

PROGS = gen_t_vs_p 

T_VS_PSRCS = $(EPIC_PATH)/data/jupiter/gen_t_vs_p.c

T_VS_POBJS = $(EPIC_PATH)/data/jupiter/gen_t_vs_p.o \
             $(EPIC_PATH)/src/shared/epic_funcs_util_se.o

HFILES = $(EPIC_PATH)/include/epic.h \
         $(EPIC_PATH)/include/avs_epic.h

SRCS = $(T_VS_PSRCS)  

# EPIC_SINGLE indicates that the computer is a single processor:

DEFINES = -DEPIC_SINGLE \
          -DEPIC_PATH=\"$(EPIC_PATH)\" 

INCL = -I. \
       -I/usr/include \
       -I$(EPIC_PATH)/include

#
# Environment variable EPIC_CFLAG is -O for optimize, -g for debug
#
CFLAGS = $(INCL) $(DEFINES) 

CC = $(EPIC_CC)
FC = $(EPIC_FC)

LFLAGS = -g

.SUFFIXES:
.SUFFIXES: .o .c  
.c.o: $(HFILES)
	if [ $(EPIC_CFLAG) = -g ] ; then \
	  ${CC} -g $(CFLAGS) -c $*.c -o $*.o ; \
	else \
	  ${CC} -O $(CFLAGS) -c $*.c -o $*.o ; \
	fi

all: $(PROGS)

gen_t_vs_p: $(T_VS_POBJS) $(HFILES)
	$(CC) $(LFLAGS) $(T_VS_POBJS) -o $@.$(MACHINE_TYPE) -lm -lc ;

clear:
	@rm -f core *o *BAK *bak *~ 

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.
