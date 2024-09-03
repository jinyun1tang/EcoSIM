set(ECOSIM_DRIVER_SOURCES
  EcoSIMAPI.F90
  ecosim.F90
)
include(add_ecosim_executable)
add_ecosim_executable(ecosim.f90.x  ${ECOSIM_DRIVER_SOURCES})
link_directories(())
target_include_directories(ecosim.f90.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/
  ${CMAKE_BINARY_DIR}/f90src/microbial_bgc/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SoilPhys
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys
  ${CMAKE_BINARY_DIR}/f90src/plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/geochem/layers_chem
  ${CMAKE_BINARY_DIR}/f90src/disturbances/
  ${CMAKE_BINARY_DIR}/f90src/modelforc/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_mods/
  ${CMAKE_BINARY_DIR}/f90src/transport/nonsalt
  ${CMAKE_BINARY_DIR}/f90src/transport/salt
  ${CMAKE_BINARY_DIR}/f90src/balances/
  ${CMAKE_BINARY_DIR}/f90src/APIs/
  ${CMAKE_BINARY_DIR}/f90src/main/)

message("ECOSIM_LIBRARIES: ${ECOSIM_LIBRARIES}")
message("ECOSIM_TPLS: ${ECOSIM_TPLS}")

#include(add_ecosim_library)
#add_ecosim_library(ecosiml ${ECOSIM_LIBRARIES})
#if (BUILD_SHARED_LIBS)
#else()
#  message("ECOSIM_LIBRARIES: ${ECOSIM_LIBRARIES}")
#  add_library(ecosiml INTERFACE)
#  foreach(lib  ${ECOSIM_LIBRARIES})
#    message("lib: ${lib}")
#    target_link_libraries(ecosiml INTERFACE ${lib})
#  endforeach()
#endif()


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosim.f90.x DESTINATION bin)
  install(TARGETS ecosim.f90.x DESTINATION ${CMAKE_BINARY_DIR}/bin)
  if (BUILD_SHARED_LIBS)
#    install(TARGETS ecosiml LIBRARY DESTINATION lib/)
  else()   
#    install(TARGETS ecosiml EXPORT ecosimlTargets)
  endif()
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
