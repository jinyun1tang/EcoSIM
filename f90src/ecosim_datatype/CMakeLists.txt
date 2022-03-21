set(ECOSIM_DATATYPE_SOURCES
  GridDataType.F90
  MicrobialDataType.F90
  SOMDataType.F90
  SoilPhysDataType.F90
  FlagDataType.F90
  SoilChemDataType.F90
  PhenologyDataType.F90
  FertilizerDataType.F90
  VegDataType.F90
  EcoSIMCtrlDataType.F90
  SoilHeatDatatype.F90
  SoilWaterDataType.F90
  SoilBGCDataType.F90
  PlantDataStateType.F90
  PlantDataRateType.F90
  PlantDataCharType.F90
  ClimForcDataType.F90
  LandSurfDataType.F90
  SnowDataType.F90
  PlantMngmtDataType.F90
  SurfLitterDataType.F90
  SurfSoilDataType.F90
  CanopyDataType.F90
  RootDataType.F90
  EcoSimSumDataType.F90
  EcosysBGCFluxType.F90
  AqueChemDatatype.F90
  EcoSIMHistMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/f90src/utils/)

configure_file(solutepar.h  solutepar.h COPYONLY)


include(add_ecosim_library)
add_ecosim_library(ecosim_datatype ${ECOSIM_DATATYPE_SOURCES})

set(ECOSIM_LIBRARIES ecosim_datatype;${ECOSIM_LIBRARIES} PARENT_SCOPE)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosim_datatype DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/ecosim/datatype)
endif()
