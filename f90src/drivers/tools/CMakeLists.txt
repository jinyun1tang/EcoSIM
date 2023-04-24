set(TOOLS_SOURCES
	ClimTransformer.F90
  ManagementReader.F90
  ClimReader.F90
  etimerTest.F90
  restartTest.F90
  NamelistTest.F90
  HFileTest.F90
)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/mesh/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)

include(add_ecosim_library)
add_ecosim_library(tools ${TOOLS_SOURCES})

set(ECOSIM_LIBRARIES tools;${ECOSIM_LIBRARIES} PARENT_SCOPE)
set(ECOSIM_LIBRARIES tools;${ECOSIM_LIBRARIES})

link_libraries(${ECOSIM_LIBRARIES};${ECOSIM_TPLS})
include(add_ecosim_executable)
add_ecosim_executable(ClimTransformer.x ClimTransformer.F90)
add_ecosim_executable(ManagementReader.x ManagementReader.F90)
add_ecosim_executable(ClimReader.x ClimReader.F90)
add_ecosim_executable(etimerTest.x etimerTest.F90)
add_ecosim_executable(NamelistTest.x NamelistTest.F90)
add_ecosim_executable(HFileTest.x HFileTest.F90)
add_ecosim_executable(restartTest.x restartTest.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ClimTransformer.x DESTINATION bin)
  install(TARGETS ManagementReader.x DESTINATION bin)
  install(TARGETS ClimReader.x DESTINATION bin)
  install(TARGETS etimerTest.x DESTINATION bin)
  install(TARGETS NamelistTest.x DESTINATION bin)
  install(TARGETS HFileTest.x DESTINATION bin)
  install(TARGETS restartTest.x DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/ecosim/drivers/tools)
endif()
