module EcoSIMDesctruct

  implicit none


  public :: DestructEcoSIM
  contains

  subroutine DestructEcoSIM
  use MicrobialDataType , only : DestructMicrobialData
  use ChemTranspDataType  , only : DestructChemTranspData
  use SOMDataType       , only : DestructSOMData
  use FertilizerDataType, only : DestructFertilizerData
  use CanopyRadDataType       , only : DestructCanopyRad
  use SurfLitterDataType, only : DestructSurfLitter
  use SoilPropertyDataType, only : DestructSoilProperty
  use IrrigationDataType, only : DestructIrrigation
  use SoilBGCDataType, only : DestructSoilBGCData
  use SedimentDataType, only : DestructSedimentData
  use SoilWaterDataType, only : DestructSoilWater
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
  end subroutine DestructEcoSIM

end module EcoSIMDesctruct
