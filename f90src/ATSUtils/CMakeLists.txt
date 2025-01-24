set(ATSUTILS_MOD_SOURCES
  SharedDataMod.F90
  ATSUtilsMod.F90
  ATSCPLMod.F90
  ATSEcoSIMInitMod.F90
  ATSEcoSIMAdvanceMod.F90
  BGC_containers.F90
  c_f_interface_module.F90
)
include(add_ecosim_library)
add_ecosim_library(ATSUtils_mods ${ATSUTILS_MOD_SOURCES})

target_include_directories(ATSUtils_mods  PRIVATE
   ${CMAKE_BINARY_DIR}/f90src/Utils/
   ${CMAKE_BINARY_DIR}/f90src/Minimath/
   ${CMAKE_BINARY_DIR}/f90src/Modelconfig/
   ${CMAKE_BINARY_DIR}/f90src/Mesh/
   ${CMAKE_BINARY_DIR}/f90src/Modelpars/
   ${CMAKE_BINARY_DIR}/f90src/Balances/   
   ${CMAKE_BINARY_DIR}/f90src/Ecosim_datatype/
   ${CMAKE_BINARY_DIR}/f90src/HydroTherm/PhysData/
   ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SoilPhys/
   ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SurfPhys/
   ${CMAKE_BINARY_DIR}/f90src/HydroTherm/SnowPhys/
   ${CMAKE_BINARY_DIR}/f90src/APIData/
   ${CMAKE_BINARY_DIR}/f90src/Ecosim_mods/)


set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ATSUtils_mods PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};ATSUtils_mods)

install(TARGETS ATSUtils_mods DESTINATION ${CMAKE_BINARY_DIR}/lib)
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ATSUtils_mods DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
