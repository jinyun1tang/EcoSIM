file(GLOB BOXSHARED_SOURCES "*.F90")
#set(BOXGEOCHEM_SOURCES
#  SoluteChemDataType.F90
#  SaltChemEquilibriaMod.F90
#  ChemEquilibriaMod.F90  
#  InitSoluteMod.F90
#  SoluteParMod.F90
#)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/mesh/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelpars/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/)

include(add_ecosim_library)
add_ecosim_library(boxshared ${BOXSHARED_SOURCES})

set(ECOSIM_LIBRARIES boxshared;${ECOSIM_LIBRARIES} PARENT_SCOPE)
set(ECOSIM_LIBRARIES boxshared;${ECOSIM_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS boxshared DESTINATION ${CMAKE_BINARY_DIR}/lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
