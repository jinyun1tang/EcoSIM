set(ECOSIM_DRIVER_SOURCES
	readnamelist.F90
	soil.F90
	regressiontest.F90
	main.F90
)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_mods/)

include(add_ecosim_library)
add_ecosim_library(ecosim_driver ${ECOSIM_DRIVER_SOURCES})

set(ECOSIM_LIBRARIES ecosim_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)

link_libraries(${ECOSIM_LIBRARIES};${ECOSIM_TPLS})
include(add_ecosim_executable)
add_ecosim_executable(ecosim.f90.x main.F90 readnamelist.F90 regressiontest.F90 soil.F90)

install(TARGETS ecosim.f90.x DESTINATION bin)
