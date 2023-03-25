set(ECOSIM_MESH_SOURCES
  GridConsts.F90
  GridMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include(add_ecosim_library)
add_ecosim_library(mesh ${ECOSIM_MESH_SOURCES})

set(ECOSIM_LIBRARIES mesh;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/mesh)
endif()
