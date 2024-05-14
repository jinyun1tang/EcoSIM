set(BOX_MODEL_SOURCES
  MicrobMathFuncMod.F90
  MicAutoCPLXMod.F90
  MicBGCMod.F90
  MicFluxTypeMod.F90
  MicForcTypeMod.F90
  MicStateTraitTypeMod.F90
  MicrobeDiagTypes.F90
)

include(add_ecosim_library)
add_ecosim_library(box_model ${BOX_MODEL_SOURCES})

target_include_directories(box_model  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};box_model PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};box_model)
install(TARGETS box_model DESTINATION ${CMAKE_BINARY_DIR}/lib)  

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS box_model DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
