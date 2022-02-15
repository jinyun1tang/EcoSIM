module EcoSIMDesctruct

  implicit none


  public :: DestructEcoSIM
  contains

  subroutine DestructEcoSIM
  use MicrobialDataType, only : DestructMicrobialData
  use SoilChemDataType , only : DestructSoilChemData
  use SOMDataType      , only : DestructSOMData
  use FertilizerDataType, only : DestructFertilizerData
  implicit none

  call DestructMicrobialData

  call DestructSoilChemData

  call DestructSOMData

  call DestructFertilizerData

  end subroutine DestructEcoSIM

end module EcoSIMDesctruct
