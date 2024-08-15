#!/bin/bash

# ############################################################################ #
#
# EcoSIM build script
#
# ############################################################################ #

# A couple of helper functions:

function status_message()
{
  #Printing status messages in a nice format
  local GREEN='32m'
  if [ "${no_color}" -eq "${TRUE}" ]; then
    echo "[`date`] $1"
  else
    echo -e "[`date`]\033[$GREEN $1\033[m"
  fi
}

############# EDIT THESE! ##################


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

#This is a little confusing, but we have to move into the build dir
#and then point cmake to the top-level CMakeLists file which will
#be a directory up from the build dir hence setting
#ecosim_source_dir to ../
ecosim_source_dir='../'
ecosim_build_dir='./build/'

############## END EDIT ####################

cmake_binary=$(which cmake)

if [ "$shared" -eq 1 ]; then
    BUILDDIR="${BUILDDIR}-shared"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DBUILD_SHARED_LIBS=ON"
fi

if [ "$precision" = "double" ]; then
    BUILDDIR="${BUILDDIR}-double"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DECOSIM_PRECISION=double"
else
    BUILDDIR="${BUILDDIR}-${precision}"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DECOSIM_PRECISION=${precision}"
fi

if [ "$debug" -eq 0 ]; then
    BUILDDIR="${BUILDDIR}-Release"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_BUILD_TYPE=Release"
else
    BUILDDIR="${BUILDDIR}-Debug"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_BUILD_TYPE=Debug"
fi

if [ -n "$prefix"]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_INSTALL_PREFIX:PATH=${prefix}"
fi

if [ "$systype" == "Darwin" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DAPPLE=1"
elif [ "$systype" == "Linux" ]; then
    #Do I want this?
    #CONFIG_FLAGS="${CONFIG_FLAGS} -DLINUX=1"
    echo "We're on linux"
fi

if [ "$sanitize" -eq 1 ]; then
    BUILDDIR="${BUILDDIR}-AddressSanitizer"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DADDRESS_SANITIZER=1"
fi

echo "ATS_ECOSIM: ${ATS_ECOSIM}"


if [ -n $ATS_ECOSIM ]; then
    echo "Building ATS-EcoSIM"
    CONFIG_FLAGS="${CONFIG_FLAGS} -DATS_ECOSIM=1" 
    #Having this automatically set to use mpi compilers
    if [ "$mpi" -eq 0 ]; then
    	mpi=1
    fi
fi

#if mpi is set use the mpi compilers, if not use default
if [ "$mpi" -eq 1 ]; then
  # Use MPI versions of the compilers
  if [ -n "$MPICC" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_C_COMPILER=${MPICC}"
  fi

  if [ -n "$MPICXX" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_CXX_COMPILER=${MPICXX}"
  fi

  if [ -n "$MPIF90" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_Fortran_COMPILER=${MPIF90}"
  fi
else
  # Use non-MPI versions of the compilers
  if [ -n "$CC" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_C_COMPILER=${CC}"
  fi

  if [ -n "$CXX" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_CXX_COMPILER=${CXX}"
  fi

  if [ -n "$FC" ]; then
    CONFIG_FLAGS="${CONFIG_FLAGS} -DCMAKE_Fortran_COMPILER=${FC}"
  fi
fi


#Shouldn't need this anymore
#CONFIG_FLAGS="${CONFIG_FLAGS} -DATS_ECOSIM=$ATS_ECOSIM"

echo "Config Flags: "
echo "${CONFIG_FLAGS}"

cmd_configure="${cmake_binary} \
  ${CONFIG_FLAGS}
  ${ecosim_source_dir}"

# Check if the build exists
if [ ! -d "$ecosim_build_dir" ]; then
    echo "Directory $ecosim_build_dir does not exist. Creating..."
    mkdir -p "$ecosim_build_dir"
else
    echo "Directory $ecosim_build_dir already exists."
fi

# Configure the EcoSIM build
# Note: many of the options that used to be in the CMakeLists.txt file
# have been moved here to remove redundancies

#Do I need this?
#ATS_ECOSIM=$ATS_ECOSIM

echo "cmd_configure: $cmd_configure"
echo "building in: $ecosim_build_dir"
#run the configure command
#cd ${ecosim_build_dir}
#${cmd_configure}
cd build
${cmd_configure}

#This does the build
make -j ${parallel_jobs}
#Does the install
make install
