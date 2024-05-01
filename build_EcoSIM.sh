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
#cmake
cmake_binary='which cmake'

#Compilers 
#GNU compilers should work, intel compilers are trickier
compiler_c=gcc
compiler_cxx=g++
compiler_fc=gfortran

#Compiler flags
build_c_flags=
build_cxx_flags=
build_fort_flags=
build_link_flags=

#Install prefix
ecosim_install_prefix=
build_type=
structured=
unstructured=
shared=

#This is a little confusing, but we have to move into the build dir
#and then point cmake to the top-level CMakeLists file which will
#be a directory up from the build dir hence setting
#ecosim_source_dir to ../
ecosim_source_dir='../'
ecosim_build_dir='./build/'

############## END EDIT ####################

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
cmd_configure="${cmake_binary} \
    -DCMAKE_C_FLAGS:STRING="${build_c_flags}" \
    -DCMAKE_CXX_FLAGS:STRING="${build_cxx_flags}" \
    -DCMAKE_Fortran_FLAGS:STRING="${build_fort_flags}" \
    -DCMAKE_EXE_LINKER_FLAGS:STRING="${build_link_flags}" \
    -DCMAKE_INSTALL_PREFIX:STRING=${ecosim_install_prefix} \
    -DCMAKE_BUILD_TYPE:STRING=${build_type} \
    -DENABLE_Structured:BOOL=${structured} \
    -DENABLE_Unstructured:BOOL=${unstructured} \
    -DENABLE_OpenMP:BOOL=${openmp} \
    -DBUILD_SHARED_LIBS:BOOL=${shared} \
    ${ecosim_source_dir}"

#run the configure command
cd ${ecosim_build_dir}
${cmd_configure}

#This does the build
make -j ${parallel_jobs}
if [ $? -ne 0 ]; then
    error_message "Failed to build EcoSIM"
    exit_now 50
fi
status_message "EcoSIM build complete"

#Does the install
make install
if [ $? -ne 0 ]; then
  error_message "Failed to install EcoSIM"
  exit_now 50
fi
status_message "EcoSIM install complete"