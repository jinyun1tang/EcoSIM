# Install script for directory: /home/hutx/EcoSIM/f77src/ecosim_datatype

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/hutx/EcoSIM/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ecosim/datatype" TYPE FILE FILES
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk10.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk11a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk11b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk12a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk12b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk13a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk13b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk13c.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk13d.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk14.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk15a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk15b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk16.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk17.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk18a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk18b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk19a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk19b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk19c.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk19d.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk1cp.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk1cr.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk1g.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk1n.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk1p.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk1s.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk1u.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk20a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk20b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk20c.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk20d.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk20e.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk20f.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk21a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk21b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk22a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk22b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk22c.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk2a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk2b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk2c.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk3.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk5.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk6.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk8a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk8b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk9a.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk9b.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blk9c.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blkc.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/blktest.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/filec.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/files.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/parameters.h"
    "/home/hutx/EcoSIM/f77src/ecosim_datatype/solutepar.h"
    )
endif()

