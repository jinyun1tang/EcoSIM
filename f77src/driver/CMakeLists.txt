set(ECOSYS_DRIVER_SOURCES
	readnamelist.f
	soil.f
	regressiontest.f
	main.f
)

include_directories(${ECOSYS_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f77src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f77src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f77src/ecosys_src/)

include(add_ecosim_library)
add_ecosim_library(ecosys_driver ${ECOSYS_DRIVER_SOURCES})

set(ECOSIM_LIBRARIES ecosys_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)

link_libraries(${ECOSIM_LIBRARIES};${ECOSIM_TPLS})
include(add_ecosim_executable)
add_ecosim_executable(ecosys.x main.f readnamelist.f regressiontest.f soil.f)

install(TARGETS ecosys.x DESTINATION bin)
