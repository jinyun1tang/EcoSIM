set(ECOSIM_MOD_SOURCES
    InitAllocMod.F90
	ExecMod.F90
	StarteMod.F90
	StartsMod.F90
	VisualMod.F90

)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/mesh/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelpars/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys)
include_directories(${CMAKE_BINARY_DIR}/f90src/HydroTherm/SoilPhys)
include_directories(${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData)
include_directories(${CMAKE_BINARY_DIR}/f90src/APIData/)
include_directories(${CMAKE_BINARY_DIR}/f90src/APIs/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/geochem/box_chem)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/layers_model)
include(add_ecosim_library)
add_ecosim_library(ecosim_mods ${ECOSIM_MOD_SOURCES})

set(ECOSIM_LIBRARIES ecosim_mods;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosim_mods DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/mods)
endif()
