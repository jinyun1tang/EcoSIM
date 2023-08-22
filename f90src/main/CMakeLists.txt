set(ECOSIM_MAIN_SOURCES
  InitEcoSIM.F90
  EcoSIMDesctruct.F90
)
include(add_ecosim_library)
add_ecosim_library(ecosim_main ${ECOSIM_MAIN_SOURCES})

target_include_directories(ecosim_main  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/APIData/
  ${CMAKE_BINARY_DIR}/f90src/plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/disturbances/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SoilPhys
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/layers_model
  ${CMAKE_BINARY_DIR}/f90src/modelforc/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_mods/
  ${CMAKE_BINARY_DIR}/f90src/balances/
  ${CMAKE_BINARY_DIR}/f90src/transport/nonsalt
  ${CMAKE_BINARY_DIR}/f90src/transport/salt
  ${CMAKE_BINARY_DIR}/f90src/APIs/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ecosim_main PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ecosim_main)
install(TARGETS ecosim_main DESTINATION ${CMAKE_BINARY_DIR}/lib)  

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosim_main DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
