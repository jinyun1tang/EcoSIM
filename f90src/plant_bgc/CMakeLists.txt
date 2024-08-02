set(PLANT_BGC_SOURCES
  PlantMathFuncMod.F90
  PlantDisturbByFireMod.F90
  PlantDisturbByTillageMod.F90 
  PlantDisturbByGrazingMod.F90
  InitVegBGC.F90
  RootGasMod.F90
  PlantBalMod.F90
  NutUptakeMod.F90
	UptakesMod.F90
  NoduleBGCMod.F90
  PlantNonstElmDynMod.F90
  PlantBranchMod.F90
  RootMod.F90
	StartqMod.F90
  PhotoSynsMod.F90
  LitterFallMod.F90
  PlantDisturbsMod.F90
	ExtractsMod.F90
	GrosubsMod.F90
	StomatesMod.F90
  CanopyCondsMod.F90
  UptakePars.F90
  InitPlantMod.F90
  PlantPhenolMod.F90
)
include(add_ecosim_library)
add_ecosim_library(plant_bgc ${PLANT_BGC_SOURCES})

target_include_directories(plant_bgc  PRIVATE
  ${CMAKE_BINARY_DIR}/f90src/utils/
  ${CMAKE_BINARY_DIR}/f90src/minimath/
  ${CMAKE_BINARY_DIR}/f90src/modelconfig/
  ${CMAKE_BINARY_DIR}/f90src/mesh/
  ${CMAKE_BINARY_DIR}/f90src/ecosim_datatype/
  ${CMAKE_BINARY_DIR}/f90src/APIData/
  ${CMAKE_BINARY_DIR}/f90src/modelpars/
  ${CMAKE_BINARY_DIR}/f90src/ioutils/)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};plant_bgc PARENT_SCOPE)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};plant_bgc)
install(TARGETS plant_bgc DESTINATION ${CMAKE_BINARY_DIR}/lib)  

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS plant_bgc DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
