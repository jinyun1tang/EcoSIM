set(GEOCHEM_SOURCES
  SoluteChemDataType.F90
  SaltChemEquilibriaMod.F90
  ChemEquilibriaMod.F90
  InitSoluteMod.F90
  SoluteParMod.F90
  SoluteMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/SoilHydroTherm/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/)
include(add_ecosim_library)
add_ecosim_library(geochem ${GEOCHEM_SOURCES})

set(ECOSIM_LIBRARIES geochem;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS geochem DESTINATION ${CMAKE_BINARY_DIR}/lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/geochem)
endif()
