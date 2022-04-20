module InitEcoSIM

  implicit none
  private
  public ::   InitModules

  contains

  subroutine InitModules(nmicbguilds)
  use NitroMod            , only : InitNitro
  use RedistMod           , only : InitRedist
  use MicrobialDataType   , only : InitMicrobialData
  use SOMDataType         , only : InitSOMData
  use ChemTranspDataType  , only : InitChemTranspData
  use FertilizerDataType  , only : InitFertilizerData
  use CanopyRadDataType   , only : InitCanopyRad
  use GrosubMod           , only : InitGrosub
  use WatsubMod           , only : initWatsub
  use PlantDisturbMod     , only : InitPlantDisturbance
  use PlantTraitDataType   , only : InitPlantTraits
  use UptakesMod           , only : InitUptake
  use Hour1Mod            , only : InitHour1
  use SoilPropertyDataType, only : InitSoilProperty
  use SurfLitterDataType  , only : InitSurfLitter
  use IrrigationDataType  , only : InitIrrigation
  use SoilBGCDataType     , only : InitSoilBGCData
  use SedimentDataType    , only : InitSedimentData
  use SoilWaterDataType   , only : InitSoilWater
  use PlantAPIData        , only : InitPlantAPIData
  use PlantMngmtDataType  , only : InitPlantMngmtData
  use GridConsts
  implicit  none
  integer                 , intent(in) :: nmicbguilds   !number of microbial guilds per group

! begin_execution
  call InitCanopyRad

  call InitPlantMngmtData

  call InitSoilProperty

  call InitSurfLitter

  call InitSedimentData

  call InitSoilWater

  call InitIrrigation

  call InitPlantAPIData(JZ,JC,JP,JSA,jcplx1,JLI,JNODS)

  call InitMicrobialData(nmicbguilds)

  call InitChemTranspData

  call InitSoilBGCData

  call InitSOMData

  call InitFertilizerData

  call InitPlantTraits

  call InitPlantDisturbance

  call InitHour1

  call InitGrosub

  call InitUptake

  call InitWatsub

  call initNitro

  call InitRedist

  end subroutine InitModules
end module InitEcoSIM
