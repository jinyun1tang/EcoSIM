set(PLANT_BGC_SOURCES
  InitVegBGC.F90
	UptakesMod.F90
	StartqMod.F90
  PlantDisturbMod.F90
  PlantDisturbsMod.F90
	ExtractsMod.F90
	GrosubsMod.F90
	StomatesMod.F90
  CanopyCondsMod.F90
  GrosubPars.F90
  UptakePars.F90
  StartqsMod.F90
  HfuncsMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/APIData/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include(add_ecosim_library)
add_ecosim_library(plant_bgc ${PLANT_BGC_SOURCES})

set(ECOSIM_LIBRARIES plant_bgc;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS plant_bgc DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/plant_bgc)
endif()
