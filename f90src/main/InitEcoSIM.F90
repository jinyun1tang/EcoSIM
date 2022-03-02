module InitEcoSIM

  implicit none
  private
  public ::   InitModules

  contains

  subroutine InitModules(nmicbguilds)
  use NitroMod          , only : InitNitro
  use RedistMod         , only : InitRedist
  use MicrobialDataType , only : InitMicrobialData
  use SOMDataType       , only : InitSOMData
  use SoilChemDataType  , only : InitSoilChemData
  use FertilizerDataType, only : InitFertilizerData
  use VegDataType       , only : InitVegData
  use GrosubMod         , only : InitGrosub
  use WatsubMod         , only : initWatsub
  use PlantDisturbMod   , only : InitPlantDisturbance
  use PhenologyDataType , only : InitPhenologyData
  use UptakeMod         , only : InitUptake
  implicit  none
  integer               , intent(in) :: nmicbguilds   !number of microbial guilds per group

! begin_execution
  call InitVegData

  call InitMicrobialData(nmicbguilds)

  call InitSoilChemData

  call InitSOMData

  call InitFertilizerData

  call InitPhenologyData

  call InitPlantDisturbance

  call InitGrosub

  call InitUptake

  call InitWatsub

  call initNitro

  call InitRedist

  end subroutine InitModules
end module InitEcoSIM
