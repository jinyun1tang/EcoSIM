set(ECOSIM_DRIVER_SOURCES
	readnamelist.F90
	soil.F90
	regressiontest.F90
)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/SoilHydroTherm/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/geochem/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_mods/)
include_directories(${CMAKE_BINARY_DIR}/f90src/chemflux/nonsalt)
include_directories(${CMAKE_BINARY_DIR}/f90src/chemflux/salt)
include_directories(${CMAKE_BINARY_DIR}/f90src/balances/)
include_directories(${CMAKE_BINARY_DIR}/f90src/APIs/)
include_directories(${CMAKE_BINARY_DIR}/f90src/main/)

include(add_ecosim_library)
add_ecosim_library(ecosim_driver ${ECOSIM_DRIVER_SOURCES})

set(ECOSIM_LIBRARIES ecosim_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)
set(ECOSIM_LIBRARIES ecosim_driver;${ECOSIM_LIBRARIES})

link_libraries(${ECOSIM_LIBRARIES};${ECOSIM_TPLS})
include(add_ecosim_executable)
add_ecosim_executable(ecosim.f90.x main.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosim.f90.x DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/ecosim/driver/main)
endif()
