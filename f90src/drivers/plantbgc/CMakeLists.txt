set(PLANT_DRIVER_SOURCES
  PlantMod.F90
)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)

include(add_ecosim_library)
add_ecosim_library(plant_driver ${PLANT_DRIVER_SOURCES})

set(ECOSIM_LIBRARIES plant_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)
set(ECOSIM_LIBRARIES plant_driver;${ECOSIM_LIBRARIES})

link_libraries(${ECOSIM_LIBRARIES};${ECOSIM_TPLS})
include(add_ecosim_executable)
add_ecosim_executable(plant.x plantdriver.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS plant.x DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
