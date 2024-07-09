# EcoSIM

A biogeochemical modeling library spins off the ecosys model.

## Download

git clone --recursive git@github.com:jinyun1tang/EcoSIM.git

## Building

make config CC=icc CXX=icpc FC=ifort 

make install CC=icc CXX=icpc FC=ifort

Find the executable under ./local/bin/ecosim.f90.x

make test CC=icc CXX=icpc FC=ifort 

will do regression test to check if the code is downloaded properly.  (Hower, because compiler difference, tests may fail.)

For gnu compilers, change icc, icpc and ifort into gcc, g++ and gfortran, respectively.

## netcdf library

It is recommended you have the netcdf library installed on your system for C, C++ and fortran. 
For mac, this can be done through either macport of homebrew.

If a system provided netcdf libraray is provided. Use the following 

make config CC=icx CXX=icpc FC=ifort netcdfsys=1

to config, and use

make install CC=icx CXX=icpc FC=ifort netcdfsys=1

to build.

There is another approach when your system does not have netcdf support. Assume using
gnu compiler, first do

make stage CC=gcc CXX=g++ FC=gfortran 

then

make config CC=gcc CXX=g++ FC=gfortran netcdfsys=1

and 

make install CC=gcc CXX=g++ FC=gfortran netcdfsys=1


## more about build

read the Makefile.

## Examples
Examples can be found in examples/run_dir.
To run the example, for instance, dryland, goes to dryland, then type

../../../local/bin/ecosim.f90.x dryland.namelist

and output can be found at dryland_maize_outputs.



