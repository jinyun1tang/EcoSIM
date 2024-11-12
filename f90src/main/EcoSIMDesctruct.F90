module EcoSIMDesctruct

  implicit none

  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: DestructEcoSIM
  contains

  subroutine DestructEcoSIM
  use EcoSimSumDataType   , only : DestructEcoSimSum
  use FlagDataType        , only : DestructFlagData
  use WatsubDataMod       , only : DestructWatSubData
  use TranspSaltMod       , only : DestructTranspSalt
  use ErosionMod          , only : DestructErosion
  use SoilHeatDatatype    , only : DestructSoilHeatData
  use SoilPhysDataType    , only : DestructSoilPhysData
  use SurfSoilDataType    , only : DestructSurfSoilData
  use GridDataType        , only : DestructGridData
  use RootDataType        , only : DestructRootData
  use SnowDataType        , only : DestructSnowData
  use PlantTraitDataType  , only : DestructPlantTraits
  use LandSurfDataType    , only : DestructLandSurfData
  use EcoSIMCtrlDataType  , only : DestructEcoSIMCtrlData
  use EcoSIMCtrlMod       , only : pft_nfid,salt_model, plant_model
  use ncdio_pio           , only : ncd_pio_closefile
  use EcosimBGCFluxType   , only : DestructEcosimBGCFluxData
  use EcoSIMHistMod       , only : DestructEcoSIMHistData
  use CanopyDataType      , only : DestructCanopyData
  use ClimForcDataType    , only : DestructClimForcData
  use AqueChemDatatype    , only : DestructAquaChem
  use MicrobialDataType   , only : DestructMicrobialData
  use ChemTranspDataType  , only : DestructChemTranspData
  use SOMDataType         , only : DestructSOMData
  use FertilizerDataType  , only : DestructFertilizerData
  use CanopyRadDataType   , only : DestructCanopyRad
  use SurfLitterDataType  , only : DestructSurfLitter
  use SoilPropertyDataType, only : DestructSoilProperty
  use PlantDataRateType   , only : DestructPlantRates
  use IrrigationDataType  , only : DestructIrrigation
  use SoilBGCDataType     , only : DestructSoilBGCData
  use SedimentDataType    , only : DestructSedimentData
  use SoilWaterDataType   , only : DestructSoilWater
  use PlantAPIData        , only : DestructPlantAPIData
  use PlantMgmtDataType  , only : DestructPlantMngmtData
  use InitSOMBGCMOD       , only : DestructSOMBGC
  use TranspNoSaltMod           , only : DestructTranspNoSalt
  use SnowPhysData        , only : DestructSnowPhysData
  use HydroThermData      , only : DestructHydroThermData
  use BalanceCheckDataType, only : DestructBalanceCheckData
  use PerturbationMod     , only : destructSoilWarming
  implicit none

  call DestructMicrobialData

  call DestructChemTranspData

  call DestructSOMData

  call DestructFertilizerData

  call DestructCanopyRad

  call DestructSurfLitter

  call DestructIrrigation

  call DestructSedimentData

  call DestructSoilWater

  call DestructPlantAPIData

  call DestructPlantMngmtData

  call  DestructPlantRates

  call DestructAquaChem

  call DestructSOMBGC

  call DestructFlagData

  if(salt_model)then
    call DestructTranspSalt
  else
    call DestructTranspNoSalt
  endif
  call DestructEcoSIMCtrlData

  call DestructCanopyData

  call DestructEcoSIMHistData

  call DestructEcosimBGCFluxData

  call DestructLandSurfData

  call DestructGridData

  call DestructPlantTraits

  call DestructRootData

  call DestructSnowData

  call DestructSoilHeatData

  call DestructSoilPhysData

  call DestructSurfSoilData

  call DestructErosion

  call DestructSnowPhysData

  call DestructWatSubData

  call DestructHydroThermData

  call DestructEcoSimSum

  if(plant_model)call ncd_pio_closefile(pft_nfid)
 
  call DestructBalanceCheckData

  call destructSoilWarming()
  end subroutine DestructEcoSIM

end module EcoSIMDesctruct
