set(ECOSIM_CONFG_SOURCES
  EcoSIMSolverPar.F90
  EcoSIMConfig.F90
  GridConsts.F90
  ElmIDMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include(add_ecosim_library)
add_ecosim_library(modelconfig ${ECOSIM_CONFG_SOURCES})


set(ECOSIM_LIBRARIES modelconfig;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/config)
endif()
