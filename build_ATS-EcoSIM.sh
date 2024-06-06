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
#cmake_binary=$(which cmake)

#Compilers 
#GNU compilers should work, intel compilers are trickier
#compiler_c=$(which gcc)
#compiler_cxx=$(which g++)
#compiler_fc=$(which gfortran)

#Compiler flags
#build_c_flags='-std=c99 -Wall -Wpedantic -Wextra -Wno-sign-compare -Wno-unused-parameter -Wno-unused-function -Wno-cast-qual -Wno-unused-but-set-variable -Wno-overlength-strings -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-discarded-qualifiers -Wno-sign-conversion -Wno-maybe-uninitialized'
#build_cxx_flags=' '
#build_fort_flags='-W -Wall -Wno-unused-variable -Wno-unused-parameter -Wno-unused-function -Wuninitialized -Wno-unused-dummy-argument'
#build_link_flags=' '

#Install prefix
#ecosim_install_prefix=local
#build_type=Debug
#structured=0
#unstructured=0
#shared=0
#openmp=0

#This is a little confusing, but we have to move into the build dir
#and then point cmake to the top-level CMakeLists file which will
#be a directory up from the build dir hence setting
#ecosim_source_dir to ../
ecosim_source_dir='../'
ecosim_build_dir='./build/'

############## END EDIT ####################

cmake_binary='which cmake'

cmd_configure="${cmake_binary} \
      -DCMAKE_C_FLAGS:STRING="${build_c_flags}" \
      -DATS_ECOSIM=$ATS_ECOSIM
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

ATS_ECOSIM=$ATS_ECOSIM
echo $ATS_ECOSIM

echo "cmd_configure: $cmd_configure"
echo "building in: $ecosim_build_dir"
#run the configure command
#cd ${ecosim_build_dir}
#${cmd_configure}
cd build
pwd
cmake ../ -DATS_ECOSIM=$ATS_ECOSIM

#This does the build
make -j ${parallel_jobs}
#if [ $? -ne 0 ]; then
#    error_message "Failed to build EcoSIM"
#    exit_now 50
#fi
#status_message "EcoSIM build complete"

#Does the install
make install
#if [ $? -ne 0 ]; then
#  error_message "Failed to install EcoSIM"
#  exit_now 50
#fi
#status_message "EcoSIM install complete"
