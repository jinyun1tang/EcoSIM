set(ECOSIM_UTIL_SOURCES
  data_kind_mod.F90
  MiniMathMod.F90
  abortutils.F90
  TestMod.F90
  clock.c
  getfilename.c
  timings.F90
  fileUtil.F90
  EcoSIMSolverPar.F90
  EcoSimConst.F90
  EcoSIMConfig.F90
)


include(add_ecosim_library)
add_ecosim_library(util ${ECOSIM_UTIL_SOURCES})


set(ECOSIM_LIBRARIES util;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/util)
endif()
