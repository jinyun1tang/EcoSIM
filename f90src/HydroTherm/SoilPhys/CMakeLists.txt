set(SOILPHYS_SOURCES
  SoilHydroParaMod.F90
  WatsubDataMod.F90
  WatsubMod.F90
)

include(add_ecosim_library)
add_ecosim_library(SoilPhys ${SOILPHYS_SOURCES})

target_include_directories(SoilPhys  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/disturbances/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SurfPhys/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SoilPhys PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};SoilPhys)
install(TARGETS SoilPhys DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS SoilPhys DESTINATION lib)  
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
