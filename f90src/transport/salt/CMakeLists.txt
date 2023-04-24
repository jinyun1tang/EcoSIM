set(SALT_MOD_SOURCES
  TrnsfrsDataMod.F90
  IngridTranspMod.F90
	TrnsfrsMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/mesh/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/SoilHydroTherm/)
include_directories(${CMAKE_BINARY_DIR}/f90src/APIs/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/layers_model)
include(add_ecosim_library)
add_ecosim_library(salt ${SALT_MOD_SOURCES})

set(ECOSIM_LIBRARIES salt;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS salt DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/salt)
endif()
