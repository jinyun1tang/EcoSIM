set(SOILPHYS_SOURCES
  SoilHydroParaMod.F90
  WatsubDataMod.F90
  WatsubMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/mesh/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelpars/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/)
include_directories(${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/)
include_directories(${CMAKE_BINARY_DIR}/f90src/HydroTherm/SurfPhys/)
include(add_ecosim_library)
add_ecosim_library(SoilPhys ${SOILPHYS_SOURCES})

set(ECOSIM_LIBRARIES SoilPhys;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
#  install(TARGETS SoilPhys DESTINATION ${CMAKE_BINARY_DIR}/lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
