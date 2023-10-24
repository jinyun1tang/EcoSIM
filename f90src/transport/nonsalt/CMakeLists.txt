set(NONSALT_MOD_SOURCES
  TranspNoSaltDataMod.F90
  SurfaceFluxMod.F90
  InsideTranspMod.F90
  BoundaryTranspMod.F90
	TranspNoSaltMod.F90
)
include(add_ecosim_library)
add_ecosim_library(nonsalt ${NONSALT_MOD_SOURCES})

target_include_directories(nonsalt  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/
  ${CMAKE_BINARY_DIR}/f90src/SoilHydroTherm/
  ${CMAKE_BINARY_DIR}/f90src/APIs/
  ${CMAKE_BINARY_DIR}/f90src/plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/layers_model)


set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};nonsalt PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};nonsalt)
install(TARGETS nonsalt DESTINATION ${CMAKE_BINARY_DIR}/lib)  

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS nonsalt DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
