set(ECOSIM_INIT_SOURCES
	StartqMod.f
	StartsMod.f
)

include_directories(${CMAKE_BINARY_DIR}/f77src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f77src/ecosim_datatype/)
add_ecosim_library(ecosim_init ${ECOSIM_INIT_SOURCES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosim_init DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosys/ecosim_init)
endif()

