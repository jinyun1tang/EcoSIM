set(ECOSIM_DRIVER_SOURCES
	readnamelist.f
	soil.f
	regressiontest.f
	main.f
)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f77src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f77src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f77src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f77src/ecosim_mods/)

include(add_ecosim_library)
add_ecosim_library(ecosim_driver ${ECOSIM_DRIVER_SOURCES})

set(ECOSIM_LIBRARIES ecosim_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)

link_libraries(${ECOSIM_LIBRARIES};${ECOSIM_TPLS})
include(add_ecosim_executable)
add_ecosim_executable(ecosim.x main.f readnamelist.f regressiontest.f soil.f)

install(TARGETS ecosim.x DESTINATION bin)
