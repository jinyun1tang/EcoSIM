# NOTE: this is the old CMake file for compiling with a UNIX Makefile
# Minimum CMake version, determined by hdf5.
cmake_minimum_required (VERSION 3.4)

# Version numbers.
set (ECOSIM_MAJOR_VERSION 0)
set (ECOSIM_MINOR_VERSION 1)
set (ECOSIM_PATCH_VERSION 0)
set (ECOSIM_VERSION "${ECOSIM_MAJOR_VERSION}.${ECOSIM_MINOR_VERSION}.${ECOSIM_PATCH_VERSION}")

message("We are in the EcoSIM top-level CMakeLists")

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Options for building EcoSIM. These come from the xSDK compliance rules.
#ioption(USE_XSDK_DEFAULTS "Set to use xSDK defaults for options [ON]." ON)
#option(CMAKE_INSTALL_PREFIX "Sets installation prefix [/usr/local].")
#option(XSDK_ENABLE_DEBUG "Enables Debug mode builds [OFF]." OFF)
#option(BUILD_SHARED_LIBS "Builds shared libraries [OFF]." OFF)
#option(XSDK_WITH_NETCDF "Enables support for netcdf [OFF]." ON)
#option(TPL_NETCDF_LIBRARIES "List of absolute paths to netcdf link libraries [].")
#option(TPL_NETCDF_INCLUDE_DIRS "List of absolute paths to netcdf include directories [].")

message(STATUS "ATS_ECOSIM_BUILD: ${ATS_ECOSIM_BUILD}")
if (ATS_ECOSIM_BUILD)
  message(STATUS "Building EcoSIM within ATS")
endif()

message(STATUS "Environment variables:")
foreach (env_var IN LISTS ENV)
    message(STATUS "${env_var}=${$ENV{env_var}}")
endforeach()

if (NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX local)
endif()

# set the number of cores to something reasonable so we don't over
# subscribe a laptop or exceed hpc login node restrictions. could
# switch to something more sophisticated like polymec if necessary.
set(NUMBER_OF_CORES 4)
math(EXPR NUM_BUILD_THREADS "${NUMBER_OF_CORES} + 1")

include(set_up_platform)
include(set_up_compilers)

# Make sure compilers are set. This must be done before enabling languages.
if (NOT CMAKE_C_COMPILER)
  if (NOT $ENV{CC} STREQUAL "")
    set(CMAKE_C_COMPILER $ENV{CC})
  else()
    set(CMAKE_C_COMPILER cc)
  endif()
endif()
if (NOT CMAKE_C_FLAGS)
  set(CMAKE_C_FLAGS $ENV{CFLAGS})
endif()
if (NOT CMAKE_CXX_COMPILER)
  if (NOT $ENV{CXX} STREQUAL "")
    set(CMAKE_CXX_COMPILER $ENV{CXX})
  else()
    set(CMAKE_CXX_COMPILER c++)
  endif()
endif()
if (NOT CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS $ENV{CXX_FLAGS})
endif()
if (NOT CMAKE_Fortran_COMPILER)
  if (NOT $ENV{FC} STREQUAL "")
    set(CMAKE_Fortran_COMPILER $ENV{FC})
  else()
    set(CMAKE_Fortran_COMPILER gfortran)
  endif()
endif()
if (NOT CMAKE_Fortran_FLAGS)
  set(CMAKE_Fortran_FLAGS $ENV{FCFLAGS})
endif()

set_up_platform()
enable_language(C)
enable_language(CXX)
enable_language(Fortran)


# We declare the project here.
project (ecosim)

message("-- C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})")
message("-- CXX compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID})")
message("-- Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")

set_up_compilers()

if (BGC)
  set (ECOSIM_BGC 1)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DECOSIM_BGC")
else()
  set (ECOSIM_BGC 0)
endif()

if (ECOSIM)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DECOSIM")
endif()

# Figure out the system type.
set(ECOSIM_HAVE_BOOL 1) # All reasonable C99 compilers have this now.
if (APPLE EQUAL 1)
  set(SYS_FLAGS "-DAPPLE=1")
  set(DYLIB_SUFFIX "dylib")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
else ()
  if (LINUX EQUAL 1)
    set(SYS_FLAGS "-DLINUX=1")
    set(DYLIB_SUFFIX "so")
  else()
    if (WIN32 EQUAL 1)
      set(ECOSIM_HAVE_BOOL 0) # MS doesn't have reasonable C compilers.
      set(SYS_FLAGS "-DWINDOWS=1")
      set(DYLIB_SUFFIX "dll")
    endif()
  endif ()
endif ()

# Here we make sure CMake-installed binaries use the correct runpath, and
# that the path is not stripped during installation.
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_BINARY_DIR}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Basic libraries to be linked in.
# added the standard mathematical library for C and C++
set(ECOSIM_LIBRARIES "")

set(ECOSIM_TPLS m)
if (${NEED_LAPACK})
  # NEED_LAPACK is set by set_up_platform()
  include(FindBLAS)
  include(FindLAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  if (${LAPACK_LIBRARY_DIR})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LAPACK_LIBRARY_DIR}")
  endif()
  if (${BLAS_LIBRARY_DIR})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${BLAS_LIBRARY_DIR}")
  endif()
  set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};${LAPACK_LIBRARIES};${BLAS_LIBRARIES})
endif()

set(ECOSIM_NEED_PFUNIT 0)
set(ECOSIM_NEED_PETSC 0)
set(ECOSIM_INCLUDE_DIRS "" CACHE LIST "Include directory of EcoSIM")

message("TPL_NETCDF_LIBRARIES: ${TPL_NETCDF_LIBRARIES}")
message("TPL_NETCDF_INCLUDE_DIRS: ${TPL_NETCDF_INCLUDE_DIRS}")
message("ECOSIM_LIBRARIES: ${ECOSIM_LIBRARIES}")

message("XSDK_WITH_NETCDF: ${XSDK_WITH_NETCDF}")
message("MACHINE: ${MACHINE}")

#if (XSDK_WITH_NETCDF AND MACHINE MATCHES "supported")
#  if ((TPL_NETCDF_LIBRARIES MATCHES "lib") AND (TPL_NETCDF_INCLUDE_DIRS MATCHES "include"))
#     message("ECOSIM has netcdf support from the system")
#     message("TPL_NETCDF_LIBRARIES ${TPL_NETCDF_LIBRARIES}")
#     message("TPL_NETCDF_INCLUDE_DIRS ${TPL_NETCDF_INCLUDE_DIRS}")
#     list(APPEND ECOSIM_TPLS ${TPL_NETCDF_LIBRARIES})
#     list(APPEND ECOSIM_INCLUDE_DIRS ${TPL_NETCDF_INCLUDE_DIRS})
#     set(ECOSIM_HAVE_NETCDF 1)
#  else()
#     message("ECOSIM does not have netcdf support from the system")
#     set(ECOSIM_HAVE_NETCDF 0)
#  endif()
#endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/ecosim.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/ecosim.cmake"
  @ONLY
)

# Other third-party libraries.
add_subdirectory(3rd-partylibs)

if (NOT ECOSIM_HAVE_NETCDF)
  message(("TPL_NETCDF_LIBRARIES: ${TPL_NETCDF_LIBRARIES}"))
  list(APPEND ECOSIM_TPLS ${TPL_NETCDF_LIBRARIES})
  list(APPEND ECOSIM_INCLUDE_DIRS ${TPL_NETCDF_INCLUDE_DIRS})
endif()
# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.

include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")

list(APPEND ECOSIM_INCLUDE_DIRS "${PROJECT_BINARY_DIR}/include")
include_directories(${ECOSIM_INCLUDE_DIRS})

enable_testing()

# Source code itself.
include_directories("${PROJECT_SOURCE_DIR}")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not specified, defaulting to Debug")
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

add_subdirectory(f90src)
add_subdirectory(drivers)
