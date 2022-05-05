set(APIS_SOURCES
	PlantAPI.F90
)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/APIData/)
include_directories(${CMAKE_BINARY_DIR}/f90src/microbial_bgc/)
#include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_mods/)
include_directories(${CMAKE_BINARY_DIR}/f90src/plant_bgc/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_mods/)

include(add_ecosim_library)
add_ecosim_library(APIs ${APIS_SOURCES})

set(ECOSIM_LIBRARIES APIs;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS APIs DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/APIs)
endif()
