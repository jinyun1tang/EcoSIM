file(GLOB AQUACHEM_DRIVER_SOURCES "*.F90")
#set(AQUACHEM_DRIVER_SOURCES
#  AquaSaltChemMod.F90
#	AquachemMod.F90
#)

include(add_ecosim_executable)
add_ecosim_executable(aquachem.x  ${AQUACHEM_DRIVER_SOURCES})

target_include_directories(aquachem.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/
  ${CMAKE_BINARY_DIR}/f90src/geochem/box_chem
  ${CMAKE_BINARY_DIR}/drivers/boxshared
  ${CMAKE_BINARY_DIR}/f90src/APIs/)

#include(add_ecosim_library)
#add_ecosim_library(aquachem_driver ${AQUACHEM_DRIVER_SOURCES})

#set(ECOSIM_LIBRARIES aquachem_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS aquachem.x DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
