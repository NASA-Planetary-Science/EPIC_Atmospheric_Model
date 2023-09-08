# EPIC Atmospheric Model
Explicit Planetary Isentropic/Isobaric-Coordinate (EPIC) Atmospheric Model <br>
Multi-planet, hydrostatic global circulation model with isentropic, isobaric, or hybrid vertical coordinate

Required: Unix operating system

Included third-party, GPL-licence open software: netCDF, a self-describing file format
    
A host of software packages will input the self-describing netcdf files used by EPIC (epic.nc and extract.nc).
They tend to have an Earth-continent grid superimposed by default, but this can easily be turned off. For viewing .nc files, consider the package Panoply, which is free from NASA GISS: http://www.giss.nasa.gov/tools/panoply/

Optional: MPI (Message Passing Interface), for running EPIC on parallel computers. <br>
Free versions of MPI are available, including openmpi and mpich, and commercial versions are available.

Sample tools for Matlab, IDL, and Mathematica are included in $EPIC_PATH/tools. Additions are welcome. 
    
A. INSTALLING 

When source code is not managed (obtained) via GitHub:
  
1. To decompress and de-archive the file epic.tar.gz:

a) If a directory named epic already exists, delete it (or move it): <br>
```
rm -r -f epic
```
b) Type <br>
```
gunzip epic.tar.gz
tar xpf epic.tar
gzip epic.tar
```
This will produce a directory tree similar to:
```
%cd $EPIC_PATH
%ls -l
```
```
epic/           Top EPIC-model directory <br>
  bin/          Executable files <br>
  data/         Data files (planetary data, chemical data, etc.) <br>
  fftw/         Fast Fourier Transform code for parallel processors <br>
  help/         Help files for executables <br>
  include/      Header (.h) files <br>
  lib/          software libraries and/or symbolic links to libraries
  License.txt   Software license
  netcdf/       Source code for Network Common Data Format
  notes/        Brief remarks on version history, etc.
  README.txt    This file.
  src/
    attic/      A place to store valuable old code.
    chemistry/  Subroutines related to atmospheric chemistry
    clouds/     Cloud microphysics code
    core/       Dynamical core code
    mpi/        Parallel-processor specific code
    rt/         Radiative transfer code
    single/     Single-processor specific code
    subgrid/    Turbulence and other subgrid code
  tmp/          A location for temporary files
  tools/        Analysis tools that work with external software
  util/         Unix utility programs such as shell scripts
  ```
 To archive and compress the model:

  a) Clear any object code by typing
  ```
  cd $EPIC_PATH/src
  make clear
  ```
  b) Move up to the directory containing epic and type
  ```
  touch epic
  tar cf epic.tar epic
  gzip epic.tar
  ```
  Note: usage of the environment variable $EPIC_PATH is enabled after completing step 3.

  In your shell resource file (such as .cshrc), add the following environment-variable lines, and then edit them for your environment. A description of each is made below.

  For a .cshrc or .tcshrc file:
  ```
  setenv CC                  clang
  setenv EPIC_PATH           ~/epic
  setenv MPI_TYPE            openmpi
  setenv EPIC_CC             $CC
  setenv OMPI_CC             $CC
  setenv EPIC_CFLAG          -O
  setenv EPIC_PRECISION      8
  setenv MACHINE_TYPE        `$EPIC_PATH/util/get_machine_type.sh`
  ```
  For a .bashrc file:
  ```
  export CC=clang
  export EPIC_PATH=~/epic
  export MPI_TYPE=openmpi
  export EPIC_CC=${CC}
  export OMPI_CC=${CC}
  export EPIC_CFLAG=-O
  export EPIC_PRECISION=8
  export MACHINE_TYPE=`$EPIC_PATH/util/get_machine_type.sh`
  ```
  It is important that the quotes around $EPIC_PATH/util/get_machine_type.sh be backwards single quotes on both ends, \`...\`, such that the output of get_machine_type.sh is written into MACHINE_TYPE.

  EPIC_PATH is the unix path where the EPIC model is kept.

  MPI_TYPE <br>
  Assuming the syntax for .cshrc, the choices are:
  ```
    No MPI: none
  Open MPI: openmpi
       LAM: LAM
     mpich: mpich
  ```
  (Louisville tends to use openmpi and consequently it may be more reliable for EPIC.)

  EPIC_CC is the name of the C compiler to be used for EPIC source files. <br>
  OMPI_CC is the C compiler for Open MPI.
  ```
  NOTE: Use the most ANSI-compatible compiler available. For example, on the SunOS use c89, not the default cc, or install and use the most recent gcc from Gnu. For Intel-chip machines, icc is often the best choice. Louisville likes clang because of its static analyzer and helpful error messages.
  ```
  EPIC_CFLAG has the following recognized cases <br>
  ```
  "-O" for optimizing
  "-g" for debugging and extra checks of the validity of values
  "-p" for profiling
  "-pO" for profiling and optimizing
  ```
  In each case, an attempt is made in the top makefile to use the best flags depending on the platform.

  Alternatively, you can elect to not set EPIC_CFLAG ("unsetenv EPIC_CFLAG"), or set it to an unrecognized case, and then pass your own flags to the C compiler and linker via the environment variables CFLAGS and LDFLAGS, respectively.
                 
  EPIC_PRECISION is set to "4" or "8" to compile the model with single-precision or double-precision floating-point variables, respectively. Double precision is recommended.

  MACHINE_TYPE is set by the shell script get_machine_type.sh.

  In addition, for MPI you may need to add lines to your path like the following: <br>
```
/usr/local/mpi/lib/hpux/ch_shmem   (that is, ../lib/$(ARCH)/$(COMM))
/usr/local/mpi/bin
```
These changes will take effect when you next log in. You can make the changes take effect now by typing: <br>
```
source ~/.cshrc
```

4.  To compile the model:

a) Change directories to the source directory by typing <br>
```
cd $EPIC_PATH/src
```
b) If desired, clear previously compiled object code by typing <br>
```
make clear
```
If only a few changes have been made to the source code, then skipping this step will avoid unnecessary compilations. To be on the safe side, after modifications to header files like <br>
```
$EPIC_PATH/include/epic.h
```
one should use make clear.

c) Compile the model by typing
```
make install
```
If the compilation is successful, the executables will be found in $EPIC_PATH/bin.

NOTE: For Ubuntu linux, in the top Makefile in epic/src, change 'SHELL = /bin/sh' to 'SHELL = /bin/bash'.

6.  The IDL tools in $EPIC_PATH/tools/IDL use the environoment variable IDL_EPIC_PATH, which should be set to be the working directory for IDL plots and movies.

B. TROUBLESHOOTING

1. If the netCDF make fails, look at the file <br>
```
$EPIC_PATH/netcdf/src/INSTALL
```
for examples of environment variables that work for various platforms.  Use a working combination for your platform (ideally ones with "-O" flags for optimization) in your ~/.cshrc file, type "source ~/.cshrc" and then recompile.

 2. A "semget" or other semaphore error can usually be cleared up by using ipcs and ipcrm, or for example mpich's sbin/cleanipcs script.

C. INITIALIZING

The information needed to start the EPIC model is contained in a single epic.nc file. The program "initial" generates this information.

When using the isentropic vertical coordinate, occasionally an input T(p) profile has small glitches in it that cause the temperature lapse rate to be slightly superadiabatic, which results in a "dsgth < 0" type error from initial.  The program t_vs_p_smoother, which can be found in epic/tools/C, applies a weak diffusion to the temperature profile in an EPIC t_vs_p file, and this may clear up the issue. The strength of the diffusion is controlled by mudat in the code.

D. CHANGING

If you wish to change only a few parameters in an existing epic.nc file, use the program "change," which reads an epic.nc file, prompts for changes, and then writes a new file. You will probably want to customize the prompts in
```
$EPIC_PATH/src/shared/epic_change.c
```

To learn about other uses of the change program, try the "-h" or "-help" command-line flag.

E. RUNNING FROM COMMAND-LINE

Typing the name of the epic program with the command-line flag -h or -help prints the options on the screen.  For example, type:
```
%epic -h
```
to get a list of command-line flags and their meanings.
      
To integrate the file epic.nc forward 10000 timesteps, saving twice and backing up every 1000, with error messages to be written to the file ``epic.log'', type:
```
%epic -itrun 10000 -itsave 5000 -itback 1000 epic.nc >& epic.log &
```

To write to an extract.nc file every 100 steps, first choose the variables to extract using initial or change, and then run the model with:
```
%epic -itrun 10000 -itsave 5000 -itback 1000 -itextract 100 epic.nc >& epic.log &
```

To append to an existing extract.nc file called extract2.nc, type
```
%epic -itrun 10000 -itsave 5000 -itback 1000 -itextract 100 -append extract2.nc epic.nc >& epic.log &
```
