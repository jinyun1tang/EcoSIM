set(IOUTILS_SOURCES 
	BLOCKDATA001.f
	foutp.f
	fouts.f
	outpd.f
	outph.f
	outsd.f
	outsh.f
	readi.f
	readq.f
	reads.f
	routp.f
	routq.f
	routs.f
	split.f
	splitc.f
	woutp.f
	woutq.f
	wouts.f
	splitp.c
	splits.c
	)
include_directories(${CMAKE_BINARY_DIR}/f77src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f77src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f77src/ioutils/ecosim_init/)
include_directories(${CMAKE_BINARY_DIR}/f77src/ecosim_init/)
include(add_ecosim_library)
add_ecosim_library(ioutils ${IOUTILS_SOURCES})

set(ECOSIM_LIBRARIES ioutils;${ECOSIM_LIBRARIES} PARENT_SCOPE)

message(STATUS "lib = ${ECOSIM_LIBRARIES}")
message(STATUS "iou = ${IOUTILS_SOURCES}")

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ioutils DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosys/ioutils)
endif()
