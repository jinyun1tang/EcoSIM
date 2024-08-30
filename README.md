# EcoSIM

A biogeochemical modeling library spins off the ecosys model.

## Download

git clone --recursive git@github.com:jinyun1tang/EcoSIM.git

## Building

use the script to build the code:

sh build_EcoSIM.sh

you can also run as a bash script:

./build_EcoSIM.sh

This should build and link all the required tpls

The executable will be found in ./build/bin/ecosim.f90.x

This has been tested with gcc on multiple systems, but not intel compilers

## additional options

Opening the build_EcoSIM.sh script there are some opttional parameters you can set:

debug=0
mpi=0
shared=0
verbose=0
sanitize=0
regression_test=0

precision="double"
prefix=""
systype=""

#Leave empty to just use the environment variable compiler
CC=""
CXX=""
FC=""

You can also set options via the command line, for example:

./build_EcoSIM.sh CC=/path/to/cc CXX=/path/to/cxx FC=/path/to/fortran --debug --regression_test

run

./build_EcoSIM.sh --help 

for a full list of optional arguments
