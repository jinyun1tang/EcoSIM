set(IOUTILS_SOURCES
    ReadManagementMod.F90
    ForcWriterMod.F90
	ClimReadMod.F90
    PlantInfoMod.F90
	HistDataType.F90
    HistFileMod.F90
    restUtilMod.F90
    RestartMod.F90
	readimod.F90
	readsmod.F90
    bhistMod.F90
)
include(add_ecosim_library)
add_ecosim_library(ioutils ${IOUTILS_SOURCES})
#add_library(ioutils OBJECT ${IOUTILS_SOURCES})
#install(TARGETS ioutils OBJECTS DESTINATION lib)

target_include_directories(ioutils  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/Utils/
  ${CMAKE_BINARY_DIR}/f90src/Minimath/
  ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/DebugTools/   
  ${CMAKE_BINARY_DIR}/f90src/Mesh/
  ${CMAKE_BINARY_DIR}/f90src/Modelpars/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/Plant_bgc/
  ${CMAKE_BINARY_DIR}/f90src/Microbial_bgc/Layers_Micmodel/  
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SoilPhys/
  ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
  ${CMAKE_BINARY_DIR}/f90src/Disturbances/
  ${CMAKE_BINARY_DIR}/f90src/Ecosim_mods/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ioutils PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ioutils)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ioutils DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
