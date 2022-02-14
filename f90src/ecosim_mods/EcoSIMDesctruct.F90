module EcoSIMDesctruct

  implicit none


  public :: DestructEcoSIM
  contains

  subroutine DestructEcoSIM
  use MicrobialDataType, only : DestructMicrobialData
  use SoilChemDataType , only : DestructSoilChemData
  use SOMDataType      , only : DestructSOMData
  implicit none

  call DestructMicrobialData

  call DestructSoilChemData

  call DestructSOMData
  end subroutine DestructEcoSIM

end module EcoSIMDesctruct
