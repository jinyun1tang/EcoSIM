set(MODELDIAGS_MOD_SOURCES
  BalancesMod.F90
)
include(add_ecosim_library)
add_ecosim_library(ModelDiags ${MODELDIAGS_MOD_SOURCES})

target_include_directories(ModelDiags  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ModelDiags PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ModelDiags)
install(TARGETS ModelDiags DESTINATION ${CMAKE_BINARY_DIR}/lib)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ModelDiags DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
