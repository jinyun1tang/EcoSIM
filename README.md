# EcoSIM

A biogeochemical modeling library spins off the ecosys model.

## Download

git clone --recursive git@github.com:jinyun1tang/EcoSIM.git

## Building

make config CC=icc CXX=icpc FC=ifort F90=1

make install CC=icc CXX=icpc FC=ifort F90=1

Find the executable under ./local/bin/ecosim.f90.x

make test CC=icc CXX=icpc FC=ifort F90=1

will check if the code is downloaded properly.  (Hower, because compiler difference, tests may fail.)

For gnu compilers, change icc, icpc and ifort into gcc, g++ and gfortran, respectively.

## Examples
Examples can be found in examples/run_dir.
To run the example, for instance, dryland, goes to dryland, then type

../../../local/bin/ecosim.f90.x dryland.namelist

and output can be found at dryland_maize_outputs.
