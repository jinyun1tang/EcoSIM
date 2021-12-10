set(ECOSYS_SOURCES
	day.f
	erosion.f
	exec.f
	extract.f
	hfunc.f
	hour1.f
	redist.f
	starte.f
	startq.f
	starts.f
	stomate.f
	trnsfr.f
	trnsfrs.f
	uptake.f
	visual.f
	watsub.f
	wthr.f
	)
include_directories(${CMAKE_BINARY_DIR}/f77src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f77src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f77src/ioutils/)
include(add_ecosim_library)
add_ecosim_library(ecosys_subs ${ECOSYS_SOURCES})

set(ECOSIM_LIBRARIES ecosys_subs;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosys_subs DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosys/subs)
endif()
