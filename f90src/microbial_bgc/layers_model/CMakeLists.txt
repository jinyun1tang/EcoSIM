set(LAYERS_SOURCES
  InitSOMBGCMod.F90
  NitroDisturbMod.F90
	NitrosMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model)
include(add_ecosim_library)
add_ecosim_library(layers_model ${LAYERS_SOURCES})

set(ECOSIM_LIBRARIES layers_model;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS layers_model DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/layers_model)
endif()
