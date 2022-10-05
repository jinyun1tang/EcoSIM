module InitEcoSIM

  implicit none
  private
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: InitModules
  public :: InitModules2
  contains

  subroutine InitModules(nmicbguilds)

  use SnowDataType        , only : InitSnowData
  use ErosionMod          , only : InitErosion
  use FlagDataType        , only : InitFlagData
  use SurfSoilDataType    , only : InitSurfSoilData
  use SoilPhysDataType    , only : InitSoilPhysData
  use SoilHeatDatatype    , only : InitSoilHeatData
  use GridDataType        , only : InitGridData
  use RootDataType        , only : InitRootData
  use LandSurfDataType    , only : InitLandSurfData
  use EcosimBGCFluxType   , only : InitEcosimBGCFluxData
  use CanopyDataType      , only : InitCanopyData
  use EcoSIMCtrlDataType  , only : InitEcoSIMCtrlData
  use EcoSIMHistMod       , only : InitEcoSIMHistData
  use ClimForcDataType    , only : InitClimForcData
  use NitrosMod            , only : InitNitro
  use RedistMod           , only : InitRedist
  use MicrobialDataType   , only : InitMicrobialData
  use SOMDataType         , only : InitSOMData
  use ChemTranspDataType  , only : InitChemTranspData
  use FertilizerDataType  , only : InitFertilizerData
  use CanopyRadDataType   , only : InitCanopyRad
  use GrosubsMod           , only : InitGrosub
  use WatsubMod           , only : initWatsub
  use AqueChemDatatype    , only : initaquachem
  use PlantDataRateType   , only : InitPlantRates
  use PlantDisturbsMod     , only : InitPlantDisturbance
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
  use InitSOMBGCMod       , only : InitSOMBGC
  use GridConsts
  use EcoSIMConfig        , only : jcplx1 => jcplx1c
  implicit  none
  integer                 , intent(in) :: nmicbguilds   !number of microbial guilds per group

! begin_execution
  call InitSOMBGC(nmicbguilds)

  call InitGridData

  call InitLandSurfData

  call InitEcoSIMCtrlData

  call InitCanopyData

  call InitCanopyRad

  call InitAquaChem

  call InitPlantMngmtData

  call InitPlantRates

  call InitSoilProperty

  call InitSurfLitter

  call InitSedimentData

  call InitSoilWater

  call InitIrrigation

  call InitPlantAPIData(JZ,JC,JP,JSA,jcplx1,JLI,JLA,JNODS)

  call InitMicrobialData

  call InitChemTranspData

  call InitSoilBGCData

  call InitSOMData

  call InitFertilizerData

  call InitPlantTraits

  call InitFlagData

  call InitPlantDisturbance

  call InitHour1

  call InitGrosub

  call InitUptake

  call InitWatsub

  call initNitro

  call InitRedist

  call InitRootData

  call InitClimForcData

  call InitEcosimBGCFluxData

  call InitSnowData

  call InitEcoSIMHistData

  call InitSoilHeatData

  call InitSurfSoilData

  call InitSoilPhysData

  call InitErosion

  end subroutine InitModules

!------------------------------------------------------------------------------------------

  subroutine InitModules2

  use FlagDataType , only : ISALTG
  use TrnsfrsMod   , only : InitTrnsfrs
  use TrnsfrMod    , only : InitTrnsfr

  implicit none


  if(ISALTG/=0)then
    call InitTrnsfrs
  else
    call InitTrnsfr
  endif


  end subroutine InitModules2
end module InitEcoSIM
