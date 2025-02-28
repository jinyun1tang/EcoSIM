# EcoSIM

A biogeochemical modeling library spins off the ecosys model.

## Download

git clone --recursive git@github.com:jinyun1tang/EcoSIM.git

## Building

Use the script in the root directory to build the code:

sh build_EcoSIM.sh

or run it as a bash script:

./build_EcoSIM.sh

This should do everything including building and linking all the required tpls, and configuring, building and linking EcoSIM itself.

The executable will be found in ./build/bin/ecosim.f90.x

This has been tested with gcc on multiple systems, but not intel compilers. See the .github/workflow/ecosim-ci.yml file for examples. You will need cmake to use the build_EcoSIM.sh script.  

## Additional options

The build_EcoSIM.sh script lets you set some opttional parameters:

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

Run

./build_EcoSIM.sh --help 

for a full list of optional arguments
