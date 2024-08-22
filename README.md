# EcoSIM

A biogeochemical modeling library spins off the ecosys model.

## Download

git clone --recursive git@github.com:jinyun1tang/EcoSIM.git

## Building

use the shell script to build the code:

sh build_EcoSIM.sh

This should build and link all the required tpls

The executable will be found in ./build/bin/ecosim.f90.x

This has been tested with gcc on multiple systems, have not tried intel yet.

## additional options

Opening the build_EcoSIM.sh script there are some opttional parameters you can set:

debug=0
mpi=0
shared=0
verbose=0
sanitize=0
travis=0

precision="double"
prefix=""
systype=""


#Leave empty to just use the environment variable compiler
CC=""
CXX=""
FC=""

You shouldn't have to modify any of these, but they are there if you need them.



