set(IOUTILS_SOURCES
	Hist1Mod.F90
	readimod.F90
	readqmod.F90
	readsmod.F90
	routp.F90
	routq.F90
	routs.F90
	split.F90
	splitc.F90
	woutp.F90
	woutq.F90
	wouts.F90
	splitp.c
	splits.c
  ncdio_pio.F90
  histMod.F90
	)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_mods/)
include(add_ecosim_library)
add_ecosim_library(ioutils ${IOUTILS_SOURCES})

set(ECOSIM_LIBRARIES ioutils;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ioutils DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/ioutils)
endif()
