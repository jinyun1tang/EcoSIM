set(MOCK_DRIVER_SOURCES
  MockMod.F90
  mockdriver.F90
)

include_directories(${ECOSIM_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/ioutils/)

#include(add_ecosim_library)
#add_ecosim_library(mock_driver ${MOCK_DRIVER_SOURCES})

#set(ECOSIM_LIBRARIES mock_driver;${ECOSIM_LIBRARIES} PARENT_SCOPE)

#link_libraries(${ECOSIM_LIBRARIES};${ECOSIM_TPLS})
include(add_ecosim_executable)
add_ecosim_executable(mock.x ${MOCK_DRIVER_SOURCES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mock.x DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
