file(GLOB BOXSBGC_DRIVER_SOURCES "*.F90")
#set(BOXSBGC_DRIVER_SOURCES
#  ForcTypeMod.F90
#  MicIDMod.F90
#  ChemMod.F90
#  batchsbgc.F90
#)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelpars/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/geochem/box_chem)
include_directories(${CMAKE_BINARY_DIR}/f90src/drivers/boxshared)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model)

#include(add_ecosim_library)
#add_ecosim_library(boxsbgc_driver ${BOXSBGC_DRIVER_SOURCES})

#set(ECOSIM_LIBRARIES boxsbgc_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)

#link_libraries(${ECOSIM_LIBRARIES};${ECOSIM_TPLS})
include(add_ecosim_executable)
add_ecosim_executable(boxsbgc.x ${BOXSBGC_DRIVER_SOURCES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS boxsbgc.x DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
