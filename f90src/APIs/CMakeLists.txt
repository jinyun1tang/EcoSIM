set(APIS_SOURCES
	PlantAPI.F90
  PlantMod.F90
  MicBGCAPI.F90
  GeochemAPI.F90
  SurfPhysAPI.F90  
)

include(add_ecosim_library)
add_ecosim_library(APIs ${APIS_SOURCES})

target_include_directories(APIs PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/APIData/
  ${CMAKE_BINARY_DIR}/f90src/plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/layers_model
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/
  ${CMAKE_BINARY_DIR}/f90src/disturbances/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_mods/
  ${CMAKE_BINARY_DIR}/f90src/balances/
  ${CMAKE_BINARY_DIR}/f90src/transport/nonsalt
  ${CMAKE_BINARY_DIR}/f90src/transport/salt
  ${CMAKE_BINARY_DIR}/f90src/geochem/box_chem
  ${CMAKE_BINARY_DIR}/f90src/geochem/layers_chem)



set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};APIs PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};APIs)

install(TARGETS APIs DESTINATION ${CMAKE_BINARY_DIR}/lib)
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS APIs DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
