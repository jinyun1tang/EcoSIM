set(MODELPARS_SOURCES
  TracerPropMod.F90
  ChemTracerParsMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include(add_ecosim_library)
add_ecosim_library(modelpars ${MODELPARS_SOURCES})

set(ECOSIM_LIBRARIES modelpars;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS modelpars DESTINATION ${CMAKE_BINARY_DIR}/lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/modelpars)
endif()
