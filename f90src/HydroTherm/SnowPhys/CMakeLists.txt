set(SNOWPHYS_SOURCES
 SnowPhysData.F90
 SnowPhysMod.F90
 SnowBalanceMod.F90
 SnowTransportMod.F90
)
include(add_ecosim_library)
add_ecosim_library(SnowPhys ${SNOWPHYS_SOURCES})

target_include_directories(SnowPhys  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SnowPhys PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SnowPhys)
install(TARGETS SnowPhys DESTINATION ${CMAKE_BINARY_DIR}/lib)


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS SnowPhys DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
