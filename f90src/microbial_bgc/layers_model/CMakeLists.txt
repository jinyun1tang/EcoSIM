set(LAYERS_SOURCES
  InitSOMBGCMod.F90
	SoilBGCNLayMod.F90
)
include(add_ecosim_library)
add_ecosim_library(layers_model ${LAYERS_SOURCES})

target_include_directories(layers_model  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/disturbances/
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};layers_model PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};layers_model)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS layers_model DESTINATION ${CMAKE_BINARY_DIR}/lib)  
  install(TARGETS layers_model DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
