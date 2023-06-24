module InitEcoSIM
!
!! DESCRIPTION
! initialize the data structure for EcoSIM
  use EcoSIMCtrlMod
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: InitModules
  public :: InitModules2
  contains

  subroutine InitModules(nmicbguilds)

  use EcoSiMParDataMod    , only : micpar,pltpar
  use EcoSimSumDataType   , only : InitEcoSimSum
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
  use GrosubsMod          , only : InitGrosub
  use WatsubMod           , only : initWatsub
  use AqueChemDatatype    , only : initaquachem
  use PlantDataRateType   , only : InitPlantRates
  use PlantDisturbsMod     , only : InitPlantDisturbance
  use PlantTraitDataType   , only : InitPlantTraits
  use UptakesMod           , only : InitUptake
  use Hour1Mod             , only : InitHour1
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
  use HistDataType        , only : hist_ecosim
  use EcoSIMConfig        , only : jcplx1 => jcplx1c
  use TracerIDMod         , only : InitTracerIDs
  use SnowPhysData        , only : InitSnowPhysData
  use HydroThermData      , only : InitHydroThermData
  implicit  none
  
  integer                 , intent(in) :: nmicbguilds   !number of microbial guilds per group

! begin_execution

  call InitSOMBGC(nmicbguilds)

  call InitPlantMorphSize()

  call InitGrosub(jpstgs,JRS)

  call InitGridData

  call InitTracerIDs(salt_model)

  call InitLandSurfData

  call InitEcoSIMCtrlData

  call InitCanopyData

  call InitCanopyRad

  call InitAquaChem

  call InitPlantMngmtData

  call InitPlantRates(micpar%n_pltlitrk,pltpar%jroots)

  call InitSoilProperty

  call InitSurfLitter(micpar%n_litrsfk)

  call InitSedimentData

  call InitSoilWater

  call InitIrrigation

  call InitPlantAPIData()

  call InitMicrobialData

  call InitChemTranspData

  call InitSoilBGCData(pltpar%n_pltlitrk)

  call InitSOMData(micpar%n_litrsfk)

  call InitFertilizerData

  call InitPlantTraits(pltpar%n_pltlitrk)

  call InitFlagData

  call InitPlantDisturbance

  call InitHour1(micpar%n_litrsfk)

  call InitUptake

  call InitWatsub

  call initNitro

  call InitRedist

  call InitRootData(pltpar%jroots)

  call InitClimForcData

  call InitEcosimBGCFluxData

  call InitSnowData

  call InitEcoSIMHistData

  call InitSoilHeatData

  call InitSurfSoilData

  call InitHydroThermData

  call InitSnowPhysData
  
  call InitSoilPhysData
  
  call InitErosion

  call InitEcoSimSum

  call InitModules2
    
  call hist_ecosim%Init(bounds)

  end subroutine InitModules

!------------------------------------------------------------------------------------------

  subroutine InitModules2

  use FlagDataType , only : ISALTG
  use TrnsfrsMod   , only : InitTrnsfrs
  use TrnsfrMod    , only : InitTrnsfr

  implicit none


  if(salt_model)then
    call InitTrnsfrs
  else
    call InitTrnsfr
  endif


  end subroutine InitModules2

!------------------------------------------------------------------------------------------
  subroutine InitPlantMorphSize()
  use EcoSiMParDataMod, only : pltpar,micpar
  use GridConsts
  implicit none

  pltpar%JZ1    = JZ
  pltpar%JC1    = JC
  pltpar%JP1    = JP
  pltpar%JLA1   = JLA
  pltpar%JSA1   = JSA
  pltpar%JLI1   = JLI
  pltpar%JNODS1 = JNODS
  pltpar%iprotein =micpar%iprotein
  pltpar%icarbhyro=micpar%icarbhyro
  pltpar%icellulos=micpar%icellulos
  pltpar%ilignin  =micpar%ilignin
  pltpar%k_woody_litr=micpar%k_woody_litr
  pltpar%k_fine_litr=micpar%k_fine_litr
  !the following variable should be consistent with the soil bgc model
  pltpar%jcplx= micpar%jcplx
  pltpar%n_pltlitrk=micpar%n_pltlitrk
  pltpar%jsken  = micpar%jsken
  pltpar%Jlitgrp= 5     !number of liter groups
  pltpar%JBR    = 10    !number of branches
  JBR=pltpar%JBR
  Jlitgrp=pltpar%Jlitgrp
  JRS=pltpar%JRS


  end subroutine InitPlantMorphSize
end module InitEcoSIM
