set(DISTURB_MOD_SOURCES
	ErosionMod.F90
  PerturbationMod.F90
  PlantDisturbMod.F90
  SoilDisturbMod.F90
)
include(add_ecosim_library)
add_ecosim_library(disturb_mods ${DISTURB_MOD_SOURCES})

target_include_directories(disturb_mods  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};disturb_mods PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};disturb_mods)
install(TARGETS disturb_mods DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS disturb_mods DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
