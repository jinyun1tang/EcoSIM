set(BOXSBGC_DRIVER_SOURCES
	batchsbgc.F90
)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/APIs/)

include(add_ecosim_library)
add_ecosim_library(boxsbgc_driver ${BOXSBGC_SOURCES})

set(ECOSIM_LIBRARIES boxsbgc_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)
set(ECOSIM_LIBRARIES boxsbgc_driver;${ECOSIM_LIBRARIES})

link_libraries(${ECOSIM_LIBRARIES};${ECOSIM_TPLS})
include(add_ecosim_executable)
add_ecosim_executable(ecosim.f90.x main.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS boxsbgc.x DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/ecosim/drivers/boxsbgc)
endif()
