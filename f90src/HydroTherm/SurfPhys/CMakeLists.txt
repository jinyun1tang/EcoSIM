set(SURFPHYS_SOURCES
  SurfPhysData.F90
  SurfPhysMod.F90
  SurfPhysAPI.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/mesh/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelpars/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/)
include_directories(${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/)
include(add_ecosim_library)
add_ecosim_library(SurfPhys ${SURFPHYS_SOURCES})

set(ECOSIM_LIBRARIES SurfPhys;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS SurfPhys DESTINATION ${CMAKE_BINARY_DIR}/lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/geochem/SurfPhys)
endif()
