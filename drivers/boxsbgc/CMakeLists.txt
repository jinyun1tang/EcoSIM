file(GLOB BOXSBGC_DRIVER_SOURCES "*.F90")
#set(BOXSBGC_DRIVER_SOURCES
#  ForcTypeMod.F90
#  MicIDMod.F90
#  ChemMod.F90
#  batchsbgc.F90
#)
include(add_ecosim_executable)
add_ecosim_executable(boxsbgc.x ${BOXSBGC_DRIVER_SOURCES})

target_include_directories(boxsbgc.x  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/IOutils/
  ${CMAKE_BINARY_DIR}/f90src/Geochem/Box_chem
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Box_Micmodel
  ${CMAKE_BINARY_DIR}/drivers/boxshared)
#include(add_ecosim_library)
#add_ecosim_library(boxsbgc_driver ${BOXSBGC_DRIVER_SOURCES})

#set(ECOSIM_LIBRARIES boxsbgc_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)

#link_libraries(${ECOSIM_LIBRARIES};${ECOSIM_TPLS})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS boxsbgc.x DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
