set(BALANCES_MOD_SOURCES
  TFlxTypeMod.F90
  LateralTranspMod.F90
  SoilLayerDynMod.F90
  RunoffBalMod.F90
  ErosionBalMod.F90
	RedistMod.F90
  TillageMixMod.F90
  BalancesMod.F90
)
include(add_ecosim_library)
add_ecosim_library(balances ${BALANCES_MOD_SOURCES})

target_include_directories(balances  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SurfPhys/
  ${CMAKE_BINARY_DIR}/f90src/disturbances/  
  ${CMAKE_BINARY_DIR}/f90src/plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/layers_model)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};balances PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};balances)
install(TARGETS balances DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS balances DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
