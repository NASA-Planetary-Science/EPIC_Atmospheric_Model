#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# *                                                                 *
# * Copyright (C) 1998-2012 Timothy E. Dowling                      *
# *                                                                 *
# * This program is free software; you can redistribute it and/or   *
# * modify it under the terms of the GNU General Public License     *
# * as published by the Free Software Foundation; either version 2  *
# * of the License, or (at your option) any later version.          *
# * A copy of this License is in the file:                          *
# *   $EPIC_PATH/License.txt                                       *
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
#  This is the single-processor Makefile for the EPIC model.
#  It also compiles initial and change.
#
#  The object-code suffix _se.o denotes single-epic.
#

PROGS = initial \
        change \
        epic

INITSRCS = $(EPIC_PATH)/src/core/epic_globals.c \
           $(EPIC_PATH)/src/core/epic_initial.c \
           $(EPIC_PATH)/src/core/epic_init_eigen.c \
           $(EPIC_PATH)/src/core/epic_funcs_io.c \
           $(EPIC_PATH)/src/core/epic_funcs_util.c \
           $(EPIC_PATH)/src/core/epic_funcs_init.c \
           $(EPIC_PATH)/src/core/epic_funcs_astron.c \
           $(EPIC_PATH)/src/core/epic_funcs_thermo.c \
           $(EPIC_PATH)/src/core/epic_funcs_diag.c \
           $(EPIC_PATH)/src/single/s_epic_functions.c \
           $(EPIC_PATH)/src/core/epic_timestep.c \
           $(EPIC_PATH)/src/core/epic_sensible_heating.c \
           $(EPIC_PATH)/src/core/epic_flux.c \
           $(EPIC_PATH)/src/core/epic_adjust.c \
           $(EPIC_PATH)/src/chemistry/epic_funcs_chem.c \
           $(EPIC_PATH)/src/clouds/epic_microphysics.c \
           $(EPIC_PATH)/src/clouds/epic_microphysics_funcs.c \
           $(EPIC_PATH)/src/subgrid/epic_subgrid.c \
           $(EPIC_PATH)/src/rt/rt_heating.c

INITOBJS = $(EPIC_PATH)/src/core/epic_globals_se.o \
           $(EPIC_PATH)/src/core/epic_initial_se.o \
           $(EPIC_PATH)/src/core/epic_init_eigen_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_io_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_util_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_init_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_astron_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_thermo_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_diag_se.o \
           $(EPIC_PATH)/src/single/s_epic_functions_se.o \
           $(EPIC_PATH)/src/core/epic_timestep_se.o \
           $(EPIC_PATH)/src/core/epic_sensible_heating_se.o \
           $(EPIC_PATH)/src/core/epic_flux_se.o \
           $(EPIC_PATH)/src/core/epic_adjust_se.o \
           $(EPIC_PATH)/src/chemistry/epic_funcs_chem_se.o \
           $(EPIC_PATH)/src/clouds/epic_microphysics_se.o \
           $(EPIC_PATH)/src/clouds/epic_microphysics_funcs_se.o \
           $(EPIC_PATH)/src/subgrid/epic_subgrid_se.o \
           $(EPIC_PATH)/src/rt/rt_heating_se.o
     

CHANGESRCS = $(EPIC_PATH)/src/core/epic_globals.c \
             $(EPIC_PATH)/src/core/epic_change.c \
             $(EPIC_PATH)/src/core/epic_init_eigen.c \
             $(EPIC_PATH)/src/core/epic_funcs_io.c \
             $(EPIC_PATH)/src/core/epic_funcs_util.c \
             $(EPIC_PATH)/src/core/epic_funcs_init.c \
             $(EPIC_PATH)/src/core/epic_funcs_astron.c \
             $(EPIC_PATH)/src/core/epic_funcs_thermo.c \
             $(EPIC_PATH)/src/core/epic_funcs_diag.c \
             $(EPIC_PATH)/src/single/s_epic_functions.c \
             $(EPIC_PATH)/src/core/epic_timestep.c \
             $(EPIC_PATH)/src/core/epic_sensible_heating.c \
             $(EPIC_PATH)/src/core/epic_flux.c \
             $(EPIC_PATH)/src/core/epic_adjust.c \
             $(EPIC_PATH)/src/chemistry/epic_funcs_chem.c \
             $(EPIC_PATH)/src/clouds/epic_microphysics.c \
             $(EPIC_PATH)/src/clouds/epic_microphysics_funcs.c \
             $(EPIC_PATH)/src/subgrid/epic_subgrid.c \
             $(EPIC_PATH)/src/rt/rt_heating.c

CHANGEOBJS = $(EPIC_PATH)/src/core/epic_globals_se.o \
             $(EPIC_PATH)/src/core/epic_change_se.o \
             $(EPIC_PATH)/src/core/epic_init_eigen_se.o \
             $(EPIC_PATH)/src/core/epic_funcs_io_se.o \
             $(EPIC_PATH)/src/core/epic_funcs_util_se.o \
             $(EPIC_PATH)/src/core/epic_funcs_init_se.o \
             $(EPIC_PATH)/src/core/epic_funcs_astron_se.o \
             $(EPIC_PATH)/src/core/epic_funcs_thermo_se.o \
             $(EPIC_PATH)/src/core/epic_funcs_diag_se.o \
             $(EPIC_PATH)/src/single/s_epic_functions_se.o \
             $(EPIC_PATH)/src/core/epic_timestep_se.o \
             $(EPIC_PATH)/src/core/epic_sensible_heating_se.o \
             $(EPIC_PATH)/src/core/epic_flux_se.o \
             $(EPIC_PATH)/src/core/epic_adjust_se.o \
             $(EPIC_PATH)/src/chemistry/epic_funcs_chem_se.o \
             $(EPIC_PATH)/src/clouds/epic_microphysics_se.o \
             $(EPIC_PATH)/src/clouds/epic_microphysics_funcs_se.o \
             $(EPIC_PATH)/src/subgrid/epic_subgrid_se.o \
             $(EPIC_PATH)/src/rt/rt_heating_se.o

EPIC_SRCS = $(EPIC_PATH)/src/core/epic_globals.c \
            $(EPIC_PATH)/src/single/s_epic_functions.c \
	    $(EPIC_PATH)/src/core/epic_main.c \
            $(EPIC_PATH)/src/core/epic_timestep.c \
            $(EPIC_PATH)/src/core/epic_sensible_heating.c \
            $(EPIC_PATH)/src/core/epic_flux.c \
            $(EPIC_PATH)/src/core/epic_funcs_io.c \
            $(EPIC_PATH)/src/core/epic_funcs_util.c \
            $(EPIC_PATH)/src/core/epic_funcs_astron.c \
            $(EPIC_PATH)/src/core/epic_funcs_thermo.c \
            $(EPIC_PATH)/src/core/epic_funcs_diag.c \
            $(EPIC_PATH)/src/core/epic_adjust.c \
            $(EPIC_PATH)/src/chemistry/epic_funcs_chem.c \
            $(EPIC_PATH)/src/clouds/epic_microphysics.c \
            $(EPIC_PATH)/src/clouds/epic_microphysics_funcs.c \
            $(EPIC_PATH)/src/subgrid/epic_subgrid.c \
            $(EPIC_PATH)/src/rt/rt_heating.c

EPICOBJS = $(EPIC_PATH)/src/core/epic_globals_se.o \
           $(EPIC_PATH)/src/single/s_epic_functions_se.o \
	   $(EPIC_PATH)/src/core/epic_main_se.o \
           $(EPIC_PATH)/src/core/epic_timestep_se.o \
           $(EPIC_PATH)/src/core/epic_sensible_heating_se.o \
           $(EPIC_PATH)/src/core/epic_flux_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_io_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_util_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_astron_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_thermo_se.o \
           $(EPIC_PATH)/src/core/epic_funcs_diag_se.o \
           $(EPIC_PATH)/src/core/epic_adjust_se.o \
           $(EPIC_PATH)/src/chemistry/epic_funcs_chem_se.o \
           $(EPIC_PATH)/src/clouds/epic_microphysics_se.o \
           $(EPIC_PATH)/src/clouds/epic_microphysics_funcs_se.o \
           $(EPIC_PATH)/src/subgrid/epic_subgrid_se.o \
           $(EPIC_PATH)/src/rt/rt_heating_se.o

HFILES = $(EPIC_PATH)/include/epic.h \
         $(EPIC_PATH)/include/epic_funcs_util.h \
         $(EPIC_PATH)/include/epic_microphysics.h \
         $(EPIC_PATH)/include/epic_subgrid.h \
         $(EPIC_PATH)/include/rt_heating.h

SRCS = $(INITSRCS) $(CHANGESRCS) $(EPIC_SRCS) 

# EPIC_SINGLE indicates that the computer is a single processor:

DEFINES = -DEPIC_SINGLE \
          -DEPIC_PATH=\"$(EPIC_PATH)\" \
          -D$(MACHINE_TYPE) \
          -DEPIC_PRECISION=$(EPIC_PRECISION) \
          -DEPIC_CHECK=$(EPIC_CHECK)

INCL = -I. \
       -I/usr/include \
       -I$(EPIC_PATH)/netcdf/include \
       -I$(EPIC_PATH)/include

C_FLAGS = $(CFLAGS) $(INCL) $(DEFINES) 

CC     = $(EPIC_CC)
LINKER = $(EPIC_CC)

LD_FLAGS = $(LDFLAGS) $(XLDFLAGS)

LIB_PATH = -L$(EPIC_PATH)/netcdf/lib \
           -L$(EPIC_PATH)/src/rt/cdisort \
           -L$(EPIC_PATH)/lib

LIB_LIST = -lnetcdf \
           -lcdisort \
           -lm

LIBS     = $(LIB_PATH) $(LIB_LIST)

.SUFFIXES:
.SUFFIXES: _se.o .c .f 
.c_se.o: $(HFILES)
	${CC} $(C_FLAGS) -c $*.c -o $*_se.o ; 

all: $(PROGS)

initial: $(INITOBJS) $(HFILES)
	$(LINKER) $(LD_FLAGS) $(INITOBJS) $(LIBS) -o $@

change: $(CHANGEOBJS) $(HFILES)
	$(LINKER) $(LD_FLAGS) $(CHANGEOBJS) $(LIBS) -o $@

epic: $(EPICOBJS) $(HFILES)
	$(LINKER) $(LD_FLAGS) $(EPICOBJS) $(LIBS) -o $@

clean:
	@rm -f core *o *BAK *bak *~ *% *.log ._*

clear:
	@rm -f core *o *BAK *bak *~ *% *.log ._*

install: all
	@for name in $(PROGS); do \
	  mv $$name $(EPIC_PATH)/bin ; \
	done

depend:
	makedepend -- $(C_FLAGS) -- $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.
