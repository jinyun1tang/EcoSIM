set(ECOSIM_MESH_SOURCES
  GridConsts.F90
  GridMod.F90
)

#include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
#include_directories(${CMAKE_BINARY_DIR}/f90src/minimath/)
#include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)


include(add_ecosim_library)
add_ecosim_library(mesh ${ECOSIM_MESH_SOURCES})
#add_library(mesh OBJECT ${ECOSIM_MESH_SOURCES})
#install(TARGETS mesh OBJECTS DESTINATION lib)  

target_include_directories(mesh PRIVATE 
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};mesh PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};mesh)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mesh DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
