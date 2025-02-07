set(TOOLS_SOURCES
	ATSEcoSIM_test.F90
)
include(add_ecosim_library)
include(add_ecosim_executable)

add_ecosim_executable(ATSEcoSIM_test.x ATSEcoSIM_test.F90)

target_include_directories(ATSEcoSIM_test.x PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/ATSUtils/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  )

#target_include_directories(ClimTransformer.x PRIVATE
#  ${CMAKE_BINARY_DIR}/f90src/utils/
#  ${CMAKE_BINARY_DIR}/f90src/mesh/
#  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
#  ${CMAKE_BINARY_DIR}/f90src/minimath/
#  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
#  ${CMAKE_BINARY_DIR}/f90src/ioutils/)


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ATSEcoSIM_test.x DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
