module InitAllocMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod

implicit none
  character(len=*),private, parameter :: mod_filename = __FILE__

  public :: InitAlloc
  contains
  subroutine InitAlloc(NOMicrobeGuilds)
  !
  !DESCRIPTION
  !allocate memeory for ecosim data

  use EcoSiMParDataMod    , only : micpar,pltpar
  use EcoSimSumDataType   , only : InitEcoSimSum
  use SnowDataType        , only : InitSnowData
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
  use MicrobialDataType   , only : InitMicrobialData
  use SOMDataType         , only : InitSOMData
  use ChemTranspDataType  , only : InitChemTranspData
  use FertilizerDataType  , only : InitFertilizerData
  use CanopyRadDataType   , only : InitCanopyRad
  use GrosubsMod          , only : InitGrosub

  use AqueChemDatatype    , only : initaquachem
  use PlantDataRateType   , only : InitPlantRates

  use PlantTraitDataType   , only : InitPlantTraits
  use SoilPropertyDataType, only : InitSoilProperty
  use SurfLitterDataType  , only : InitSurfLitter
  use IrrigationDataType  , only : InitIrrigation
  use SoilBGCDataType     , only : InitSoilBGCData
  use SedimentDataType    , only : InitSedimentData
  use SoilWaterDataType   , only : InitSoilWater
  use PlantAPIData        , only : InitPlantAPIData
  use PlantMngmtDataType  , only : InitPlantMngmtData
  use InitSOMBGCMod       , only : InitSOMBGC

  use EcoSIMConfig        , only : jcplx1 => jcplx1c
  use TracerIDMod         , only : InitTracerIDs
  use SnowPhysData        , only : InitSnowPhysData
  use HydroThermData      , only : InitHydroThermData
  use GridConsts  

  implicit none
  integer                 , intent(in) :: NOMicrobeGuilds   !number of microbial guilds per group
! begin_execution

  call InitSOMBGC(NOMicrobeGuilds)

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

  call InitEcoSimSum

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

  end subroutine InitAlloc
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
end module InitAllocMod