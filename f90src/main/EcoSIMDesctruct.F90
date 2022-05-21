module EcoSIMDesctruct

  implicit none


  public :: DestructEcoSIM
  contains

  subroutine DestructEcoSIM
  use AqueChemDatatype  , only : DestructAquaChem
  use MicrobialDataType , only : DestructMicrobialData
  use ChemTranspDataType  , only : DestructChemTranspData
  use SOMDataType       , only : DestructSOMData
  use FertilizerDataType, only : DestructFertilizerData
  use CanopyRadDataType       , only : DestructCanopyRad
  use SurfLitterDataType, only : DestructSurfLitter
  use SoilPropertyDataType, only : DestructSoilProperty
  use PlantDataRateType, only : DestructPlantRates
  use IrrigationDataType, only : DestructIrrigation
  use SoilBGCDataType, only : DestructSoilBGCData
  use SedimentDataType, only : DestructSedimentData
  use SoilWaterDataType, only : DestructSoilWater
  use PlantAPIData  , only : DestructPlantAPIData
  use PlantMngmtDataType, only : DestructPlantMngmtData
  use InitSOMBGCMOD, only : DestructSOMBGC
  use RedistMod, only : DestructRedist
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

  call DestructRedist
  end subroutine DestructEcoSIM

end module EcoSIMDesctruct
