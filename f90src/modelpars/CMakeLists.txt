set(MODELPARS_SOURCES
  NitroPars.F90
  GrosubPars.F90
  TracerPropMod.F90
  MicBGCPars.F90
  EcoSiMParDataMod.F90
  ChemTracerParsMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include(add_ecosim_library)
add_ecosim_library(modelpars ${MODELPARS_SOURCES})

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelpars PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelpars)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS modelpars DESTINATION ${CMAKE_BINARY_DIR}/lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
