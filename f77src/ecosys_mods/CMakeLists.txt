configure_file(solutepar.h  solutepar.h COPYONLY)
set(ECOSYSMOD_SOURCES
	DayMod.f
	ErosionMod.f
	ExecMod.f
	ExtractMod.f
	GrosubMod.f
	HfuncMod.f
	Hour1Mod.f
	NitroMod.f
	RedistMod.f
	SoluteMod.f
	StarteMod.f
	StartqMod.f
	StartsMod.f
	StomateMod.f
	TrnsfrMod.f
	TrnsfrsMod.f
	UptakeMod.f
	VisualMod.f
	WatsubMod.f
	WthrMod.f
)
include_directories(${CMAKE_BINARY_DIR}/f77src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f77src/ioutils/)
include(add_ecosim_library)
add_ecosim_library(ecosys_mods ${ECOSYSMOD_SOURCES})

set(ECOSIM_LIBRARIES ecosys_mods;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosys_mods DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosys/mods)
endif()
