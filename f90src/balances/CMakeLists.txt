set(BALANCES_MOD_SOURCES
  TFlxTypeMod.F90
  LateralTranspMod.F90
  SoilLayerDynMod.F90
  SnowBalMod.F90
  RunoffBalMod.F90
  ErosionBalMod.F90
	RedistMod.F90
  TillageMixMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/SoilHydroTherm/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/layers_model)
include(add_ecosim_library)
add_ecosim_library(balances ${BALANCES_MOD_SOURCES})

set(ECOSIM_LIBRARIES balances;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS balances DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/balances)
endif()
