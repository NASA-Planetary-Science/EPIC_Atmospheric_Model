#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# *                                                                 *
# * Copyright (C) 1998-2019 Timothy E. Dowling                      *
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
#  NOTE: For Ubuntu linux, change 'SHELL = /bin/sh' to 'SHELL = /bin/bash'.
#

SHELL = /bin/sh

all:
	@export CC=$(EPIC_CC) ; \
	export FC="" ; \
	export CXX="" ; \
	cd $(EPIC_PATH)/netcdf/src ; \
	if [ -f ./config.status ] ; then \
	  echo "netCDF is built" ; \
	else \
	  echo "Configuring netCDF (Network Common Data Format)..." ; \
	  ./configure --prefix=$(EPIC_PATH)/netcdf >& configure.log ; \
	  echo "Making netCDF..." ; \
	  make ; \
	  echo "Installing netCDF..." ; \
	  make install ; \
	  cd .. ; \
	  if [ -d lib64 ] ; then \
	    ln -s lib64 lib ; \
	  fi ; \
	fi ;
	@if [ "$(EPIC_CFLAG)" = "-g" ] ; then \
	  export CFLAGS="-g" ; \
	  export LDFLAGS="-g" ; \
	  export EPIC_CHECK="1" ; \
	elif [ "$(EPIC_CFLAG)" = "-O_opteron" ] ; then \
	  export CFLAGS="-march=opteron -g -O2 -pipe" ; \
	  export LDFLAGS="" ; \
	  export EPIC_CHECK="0" ; \
	elif [ "$(EPIC_CFLAG)" = "-O_CPL" ] ; then \
	  export CFLAGS="-g -O3 -march=athlon-xp" ; \
	  export LDFLAGS="" ; \
	  export EPIC_CHECK="0" ; \
	elif [ "$(EPIC_CFLAG)" = "-O_CPL2" ] ; then \
	  export CFLAGS="-g -O3" ; \
	  export LDFLAGS="" ; \
	  export EPIC_CHECK="0" ; \
	elif [ "$(EPIC_CFLAG)" = "-O_tdmac" ] ; then \
	  export CFLAGS="-g -O3" ; \
	  export LDFLAGS="" ; \
	  export EPIC_CHECK="0" ; \
	elif [ "$(EPIC_CFLAG)" = "-O" ] ; then \
	  if [ $(EPIC_CC) = clang ] ; then \
	    export CFLAGS="-g -O2" ; \
	    export LDFLAGS="-g" ; \
	    export EPIC_CHECK="0" ; \
	  elif [ $(MACHINE_TYPE) = hpux ] ; then \
	    export CFLAGS="-g -O +Onolimit" ; \
	    export EPIC_CHECK="0" ; \
	  elif [ $(MACHINE_TYPE) = sp2 ] ; then \
	    export CFLAGS="-g -O -qmaxmem=-1" ; \
	    export EPIC_CHECK="0" ; \
	  else \
	    export CFLAGS="-g -O2" ; \
	    export EPIC_CHECK="0" ; \
	  fi ; \
	elif [ "$(EPIC_CFLAG)" = "-p" ] ; then \
	  if [ $(EPIC_CC) = clang ] ; then \
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
	  if [ $(EPIC_CC) = clang ] ; then \
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
	cd $(EPIC_PATH)/src/rt/cdisort ; \
	make -f Makefile_cdisort ; \
	if [ $(MPI_TYPE) != none ] ; then \
	  cd $(EPIC_PATH)/src/mpi/mpg ; \
	  make ; \
	  cd $(EPIC_PATH)/src/mpi ; \
	  make ; \
	else \
	  echo MPI_TYPE = $(MPI_TYPE) ; \
	fi ; \
	cd $(EPIC_PATH)/src/single ; \
	make

clean:
	@cd $(EPIC_PATH)/netcdf/src ; \
	make clean ; \
	cd $(EPIC_PATH)/src/rt/cdisort ; \
	make -f Makefile_cdisort clean ; \
	cd $(EPIC_PATH)/src/mpi ; \
	make clean ; \
	cd $(EPIC_PATH)/src/single ; \
	make clean ; \
	cd $(EPIC_PATH)/src/mpi/mpg ; \
	make clean

distclean:
	@cd $(EPIC_PATH)/src ; \
	rm -f *BAK *bak *~ *.a *.log ._*; \
	cd $(EPIC_PATH)/lib ; \
	rm -f *.a ._*; \
	cd $(EPIC_PATH)/netcdf/src ; \
	make distclean ; \
	cd $(EPIC_PATH)/netcdf ; \
	rm -r -f bin include lib share ._*; \
	cd $(EPIC_PATH)/src/rt/cdisort ; \
	make -f Makefile_cdisort clear ; \
	cd $(EPIC_PATH)/src/mpi ; \
	make clear ; \
	cd $(EPIC_PATH)/src/single ; \
	make clear ; \
	cd $(EPIC_PATH)/src/mpi/mpg ; \
	make clear ; \
	cd $(EPIC_PATH) ; \
	rm -f *BAK *bak *~ *% *.log ._*; \
	cd $(EPIC_PATH)/tmp ; \
	rm -f *; \
	cd $(EPIC_PATH)/help ; \
	rm -f *BAK *bak *~ *% *.log ._*; \
	cd $(EPIC_PATH)/tools/idl ; \
	rm -f *BAK *bak *~ *% *.ps *.dat *.nc *.pca *.gif *.log ._*; \
	cd $(EPIC_PATH)/tools/matlab ; \
	DISTCLEAN_EPICncVis ; \
	cd $(EPIC_PATH)/util ; \
	rm -f *BAK *bak *~ *% *.log ._*; \
	cd $(EPIC_PATH)/include ; \
	rm -f *BAK *bak *~ *% *.log ._*; \
	cd $(EPIC_PATH)/src/core ; \
	rm -f core *o *BAK *bak *~ *% *.log ._*; \
	cd $(EPIC_PATH)/src/subgrid ; \
	rm -f core *o *BAK *bak *~ *% *.log ._*; \
	cd $(EPIC_PATH)/src/clouds ; \
	rm -f core *o *BAK *bak *~ *% *.log ._*; \
	cd $(EPIC_PATH)/src/rt ; \
	rm -f core *o *BAK *bak *~ *% *.log ._*; \
	cd $(EPIC_PATH)/src/chemistry ; \
	rm -f core *o *BAK *bak *~ *% *.log ._*; \
	cd $(EPIC_PATH)/bin ; \
	rm -f * ; \
	cd $(EPIC_PATH)/data ; \
	rm -f ._* ;

clear: distclean

install: all
	@cd $(EPIC_PATH)/src/rt/cdisort ; \
	make -f Makefile_cdisort install ; \
	if [ $(MPI_TYPE) != none ] ; then \
	  cd $(EPIC_PATH)/src/mpi/mpg ; \
	  make install ; \
	  cd $(EPIC_PATH)/src/mpi ; \
	  make install; \
	else \
	  echo MPI_TYPE = $(MPI_TYPE) ; \
	fi ; \
	cd $(EPIC_PATH)/src/single ; \
	make install

