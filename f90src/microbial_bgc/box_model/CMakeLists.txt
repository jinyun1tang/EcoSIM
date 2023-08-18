set(BOX_MODEL_SOURCES
  MicAutoCPLXMod.F90
  MicBGCMod.F90
  MicFluxTypeMod.F90
  MicForcTypeMod.F90
  MicStateTraitTypeMod.F90
  NitroDiagTypes.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelpars/)
include(add_ecosim_library)
add_ecosim_library(box_model ${BOX_MODEL_SOURCES})

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};box_model PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};box_model)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS box_model DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
