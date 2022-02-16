set(ECOSIM_MOD_SOURCES
	DayMod.F90
	ErosionMod.F90
	ExecMod.F90
	ExtractMod.F90
	HfuncMod.F90
	Hour1Mod.F90
	RedistMod.F90
	SoluteMod.F90
	StarteMod.F90
	StartqMod.F90
	StartsMod.F90
	StomateMod.F90
	TrnsfrMod.F90
	TrnsfrsMod.F90
	UptakeMod.F90
	VisualMod.F90
	WatsubMod.F90
	WthrMod.F90
  InitEcoSIM.F90
  EcoSIMDesctruct.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/)
include(add_ecosim_library)
add_ecosim_library(ecosim_mods ${ECOSIM_MOD_SOURCES})

set(ECOSIM_LIBRARIES ecosim_mods;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosim_mods DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/mods)
endif()
