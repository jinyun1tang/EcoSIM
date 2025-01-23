file(GLOB AQUACHEM_DRIVER_SOURCES "*.F90")
#set(AQUACHEM_DRIVER_SOURCES
#  AquaSaltChemMod.F90
#	AquachemMod.F90
#)

include(add_ecosim_executable)
add_ecosim_executable(aquachem.x  ${AQUACHEM_DRIVER_SOURCES})

target_include_directories(aquachem.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/Geochem/Box_chem
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
