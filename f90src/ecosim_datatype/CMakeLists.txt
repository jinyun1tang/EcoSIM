set(ECOSIM_DATATYPE_SOURCES
  GridDataType.F90
  MicrobialDataType.F90
  SOMDataType.F90
  SoilPhysDataType.F90
  FlagDataType.F90
  ChemTranspDataType.F90
  PlantTraitDataType.F90
  FertilizerDataType.F90
  CanopyRadDataType.F90
  EcoSIMCtrlDataType.F90
  SoilHeatDataType.F90
  SoilWaterDataType.F90
  SoilBGCDataType.F90
  PlantDataRateType.F90
  ClimForcDataType.F90
  LandSurfDataType.F90
  SnowDataType.F90
  PlantMngmtDataType.F90
  SurfLitterDataType.F90
  SurfSoilDataType.F90
  CanopyDataType.F90
  RootDataType.F90
  EcoSimSumDataType.F90
  EcosimBGCFluxType.F90
  AqueChemDatatype.F90
  EcoSIMHistMod.F90
  SoilPropertyDataType.F90
  IrrigationDataType.F90
  SedimentDataType.F90
)

include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)
include_directories(${CMAKE_BINARY_DIR}/f90src/modelconfig/)
include_directories(${CMAKE_BINARY_DIR}/f90src/mesh/)
include(add_ecosim_library)
add_ecosim_library(ecosim_datatype ${ECOSIM_DATATYPE_SOURCES})

set(ECOSIM_LIBRARIES ecosim_datatype;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
#  install(TARGETS ecosim_datatype DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
