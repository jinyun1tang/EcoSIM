set(MODELFORC_MOD_SOURCES
  YearMod.F90
	DayMod.F90
	Hour1Mod.F90
	WthrMod.F90
)
include(add_ecosim_library)
add_ecosim_library(modelforc ${MODELFORC_MOD_SOURCES})

target_include_directories(modelforc  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/
  ${CMAKE_BINARY_DIR}/f90src/ATSUtils/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SoilPhys
  ${CMAKE_BINARY_DIR}/f90src/APIs/
  ${CMAKE_BINARY_DIR}/f90src/plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/geochem/box_chem
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/layers_model)


set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelforc PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelforc)
install(TARGETS modelforc DESTINATION ${CMAKE_BINARY_DIR}/lib)  

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS modelforc DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
