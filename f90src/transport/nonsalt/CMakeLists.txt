set(NONSALT_MOD_SOURCES
  TransfrDataMod.F90
  SurfaceFluxMod.F90
  InsideTranspMod.F90
  BoundaryTranspMod.F90
	TrnsfrMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/mesh/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelpars/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/SoilHydroTherm/)
include_directories(${CMAKE_BINARY_DIR}/f90src/APIs/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/layers_model)
include(add_ecosim_library)
add_ecosim_library(nonsalt ${NONSALT_MOD_SOURCES})

set(ECOSIM_LIBRARIES nonsalt;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS nonsalt DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/nonsalt)
endif()
