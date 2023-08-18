set(APIS_SOURCES
	PlantAPI.F90
  MicBGCAPI.F90
  GeochemAPI.F90
  SurfPhysAPI.F90  
)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
include_directories(${CMAKE_BINARY_DIR}/f90src/mesh/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelpars/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/APIData/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/box_model)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/layers_model)
include_directories(${CMAKE_BINARY_DIR}/f90src/SoilHydroTherm/)
include_directories(${CMAKE_BINARY_DIR}/f90src/disturbances/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_mods/)
include_directories(${CMAKE_BINARY_DIR}/f90src/balances/)
include_directories(${CMAKE_BINARY_DIR}/f90src/transport/nonsalt)
include_directories(${CMAKE_BINARY_DIR}/f90src/transport/salt)
include_directories(${CMAKE_BINARY_DIR}/f90src/geochem/box_chem)
include_directories(${CMAKE_BINARY_DIR}/f90src/geochem/layers_chem)

include(add_ecosim_library)
add_ecosim_library(APIs ${APIS_SOURCES})

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};APIs PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};APIs)

#add_ecosim_library(APIs ${APIS_SOURCES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS APIs DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
