module EcoSIMDesctruct

  implicit none

  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: DestructEcoSIM
  contains

  subroutine DestructEcoSIM
  use EcoSimSumDataType   , only : DestructEcoSimSum
  use FlagDataType        , only : DestructFlagData,ISALTG
  use WatsubDataMod       , only : DestructWatSubData
  use TrnsfrsMod          , only : DestructTrnsfrs
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
  use EcoSIMCtrlMod       , only : pft_nfid,salt_model
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
  use PlantMngmtDataType  , only : DestructPlantMngmtData
  use InitSOMBGCMOD       , only : DestructSOMBGC
  use TrnsfrMod           , only : DestructTrnsfr

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
    call DestructTrnsfrs
  else
    call DestructTrnsfr
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

  call DestructWatSubData

  call DestructEcoSimSum

  call ncd_pio_closefile(pft_nfid)
  


  end subroutine DestructEcoSIM

end module EcoSIMDesctruct
