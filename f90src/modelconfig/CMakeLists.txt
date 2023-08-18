set(ECOSIM_CONFG_SOURCES
  EcoSIMSolverPar.F90
  EcoSIMCtrlMod.F90
  EcoSIMConfig.F90
  ElmIDMod.F90
  TracerIDMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include(add_ecosim_library)
add_ecosim_library(modelconfig ${ECOSIM_CONFG_SOURCES})


set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelconfig PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};modelconfig)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS modelconfig DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
