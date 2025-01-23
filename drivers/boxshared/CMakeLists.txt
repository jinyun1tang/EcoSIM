file(GLOB BOXSHARED_SOURCES "*.F90")
#set(BOXGEOCHEM_SOURCES
#  SoluteChemDataType.F90
#  SaltChemEquilibriaMod.F90
#  ChemEquilibriaMod.F90  
#  InitSoluteMod.F90
#  SoluteParMod.F90
#)
include_directories(${CMAKE_BINARY_DIR}/f90src/Utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/Minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/Mesh/)
include_directories(${CMAKE_BINARY_DIR}/f90src/Modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/Modelpars/)
include_directories(${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/IOutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/Plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/)

include(add_ecosim_library)
add_ecosim_library(boxshared ${BOXSHARED_SOURCES})

set(ECOSIM_LIBRARIES boxshared;${ECOSIM_LIBRARIES} PARENT_SCOPE)
set(ECOSIM_LIBRARIES boxshared;${ECOSIM_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS boxshared DESTINATION ${CMAKE_BINARY_DIR}/lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
