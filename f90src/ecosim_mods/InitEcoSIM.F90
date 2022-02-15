module InitEcoSIM

  implicit none
  private
  public ::   InitModules

  contains

  subroutine InitModules(nmicbguilds)
  use NitroMod         , only : InitNitro
  use RedistMod        , only : InitRedist
  use MicrobialDataType, only : InitMicrobialData
  use SOMDataType      , only : InitSOMData
  use SoilChemDataType , only : InitSoilChemData
  use FertilizerDataType, only : InitFertilizerData
  implicit  none
  integer, intent(in) :: nmicbguilds

  call InitMicrobialData(nmicbguilds)

  call InitSoilChemData

  call InitSOMData

  call InitFertilizerData

  call initNitro

  call InitRedist

  end subroutine InitModules
end module InitEcoSIM
