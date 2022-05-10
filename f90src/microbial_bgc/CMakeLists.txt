set(MICROBIAL_BGC_SOURCES
  InitSOMBGCMod.F90
  MicBGCPars.F90
  NitroDiagTypes.F90
  Nitro1LayerMod.F90
  NitroDisturbMod.F90
	NitroMod.F90
  NitroPars.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include(add_ecosim_library)
add_ecosim_library(microbial_bgc ${MICROBIAL_BGC_SOURCES})

set(ECOSIM_LIBRARIES microbial_bgc;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS microbial_bgc DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/microbial_bgc)
endif()
