set(GEOCHEML_SOURCES
  SoluteMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelpars/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/SoilHydroTherm/)
include_directories(${CMAKE_BINARY_DIR}/f90src/geochem/box_chem/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/)
include(add_ecosim_library)
add_ecosim_library(layers_chem ${GEOCHEML_SOURCES})

set(ECOSIM_LIBRARIES layers_chem;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS layers_chem DESTINATION ${CMAKE_BINARY_DIR}/lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/geochem/layers_chem)
endif()
