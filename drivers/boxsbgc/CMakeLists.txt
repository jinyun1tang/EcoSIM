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
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/
  ${CMAKE_BINARY_DIR}/f90src/geochem/box_chem
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model
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
