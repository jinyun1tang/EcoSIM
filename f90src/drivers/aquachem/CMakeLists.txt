file(GLOB AQUACHEM_DRIVER_SOURCES "*.F90")
#set(AQUACHEM_DRIVER_SOURCES
#  AquaSaltChemMod.F90
#	AquachemMod.F90
#)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/geochem/box_chem)
include_directories(${CMAKE_BINARY_DIR}/f90src/drivers/boxshared)
include_directories(${CMAKE_BINARY_DIR}/f90src/APIs/)

#include(add_ecosim_library)
#add_ecosim_library(aquachem_driver ${AQUACHEM_DRIVER_SOURCES})

#set(ECOSIM_LIBRARIES aquachem_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)

include(add_ecosim_executable)
add_ecosim_executable(aquachem.x  ${AQUACHEM_DRIVER_SOURCES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS aquachem.x DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
