file(GLOB BOXGEOCHEM_SOURCES "*.F90")
#set(BOXGEOCHEM_SOURCES
#  SoluteChemDataType.F90
#  SaltChemEquilibriaMod.F90
#  ChemEquilibriaMod.F90  
#  InitSoluteMod.F90
#  SoluteParMod.F90
#)
include(add_ecosim_library)
add_ecosim_library(box_chem ${BOXGEOCHEM_SOURCES})

target_include_directories(box_chem  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/
  ${CMAKE_BINARY_DIR}/f90src/plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/)


set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};box_chem PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};box_chem)

install(TARGETS box_chem DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS box_chem DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
