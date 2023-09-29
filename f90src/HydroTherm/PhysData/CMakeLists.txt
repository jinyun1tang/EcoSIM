set(PHYSDATA_SOURCES
 PhysPars.F90
 HydroThermData.F90
 SoilPhysParaMod.F90
)
include(add_ecosim_library)
add_ecosim_library(PhysData ${PHYSDATA_SOURCES})

target_include_directories(PhysData  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};PhysData PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};PhysData)

install(TARGETS PhysData DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS PhysData DESTINATION lib)  
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
