set(APIDATA_SOURCES
  PlantAPIData.F90
)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)

include(add_ecosim_library)
add_ecosim_library(APIData ${APIDATA_SOURCES})

set(ECOSIM_LIBRARIES APIData;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS APIData DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/APIData)
endif()
