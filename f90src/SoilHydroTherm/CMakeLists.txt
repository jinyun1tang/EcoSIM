set(SOIL_HT_SOURCES
  WatsubPars.F90
  WatsubMod.F90
)
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)
include(add_ecosim_library)
add_ecosim_library(soil_ht ${SOIL_HT_SOURCES})

set(ECOSIM_LIBRARIES soil_ht;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS soil_ht DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/soil_ht)
endif()
