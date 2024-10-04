module InitAllocMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod

implicit none
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

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
  use InitSoluteMod       , only : InitSoluteProperty
  use AqueChemDatatype    , only : initaquachem
  use PlantDataRateType   , only : InitPlantRates
  use EcoSIMCtrlMod       , only : plant_model
  use PlantTraitDataType   , only : InitPlantTraits
  use SoilPropertyDataType, only : InitSoilProperty
  use SurfLitterDataType  , only : InitSurfLitter
  use SurfPhysData        , only : InitSurfPhysData
  use IrrigationDataType  , only : InitIrrigation
  use SoilBGCDataType     , only : InitSoilBGCData
  use SedimentDataType    , only : InitSedimentData
  use SoilWaterDataType   , only : InitSoilWater
  use PlantAPIData        , only : InitPlantAPIData
  use PlantMgmtDataType  , only : InitPlantMngmtData
  use InitSOMBGCMod       , only : InitSOMBGC
  use EcoSIMConfig        , only : jcplx1 => jcplxcm1
  use TracerIDMod         , only : InitTracerIDs
  use SnowPhysData        , only : InitSnowPhysData
  use HydroThermData      , only : InitHydroThermData
  use PerturbationMod     , only : InitSoilWarming
  use GridConsts
  use WatsubMod           , only : InitWatsub
  use BalanceCheckDataType, only : InitBalanceCheckData
  implicit none
  integer                 , intent(in) :: NOMicrobeGuilds   !number of microbial guilds per group
! begin_execution

  call InitSOMBGC(NOMicrobeGuilds)

  call InitPlantMorphSize()

  if(plant_model)call InitGrosub(NumGrowthStages,MaxNumRootAxes)

  call InitGridData

  call InitTracerIDs(salt_model)

  call InitLandSurfData

  call InitEcoSIMCtrlData

  call InitCanopyData

  call InitCanopyRad

  call InitAquaChem

  call InitPlantMngmtData

  call InitPlantRates(micpar%NumOfPlantLitrCmplxs,pltpar%jroots)

  call InitSoilProperty

  call InitWatsub

  call InitSurfLitter(micpar%NumOfLitrCmplxs)

  call InitSedimentData

  call InitSoilWater

  call InitIrrigation

  call InitPlantAPIData()

  call InitMicrobialData

  call InitChemTranspData(salt_model)

  call InitSoilBGCData(pltpar%NumOfPlantLitrCmplxs)

  call InitSOMData(micpar%NumOfLitrCmplxs)

  call InitFertilizerData

  call InitPlantTraits(pltpar%NumOfPlantLitrCmplxs)

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

  if(salt_model)call InitSoluteProperty

  call InitSoilWarming
  
  call InitBalanceCheckData

  end subroutine InitAlloc
!------------------------------------------------------------------------------------------
  subroutine InitPlantMorphSize()
  use EcoSiMParDataMod, only : pltpar,micpar
  use GridConsts
  implicit none

  pltpar%JZ1    = JZ
  pltpar%NumOfCanopyLayers1    = NumOfCanopyLayers
  pltpar%JP1    = JP
  pltpar%NumOfLeafAzimuthSectors   = NumOfLeafAzimuthSectors
  pltpar%NumOfSkyAzimuSects1   = NumOfSkyAzimuSects
  pltpar%NumOfLeafZenithSectors1   = NumOfLeafZenithSectors
  pltpar%MaxNodesPerBranch1 = MaxNodesPerBranch
  pltpar%iprotein =micpar%iprotein
  pltpar%icarbhyro=micpar%icarbhyro
  pltpar%icellulos=micpar%icellulos
  pltpar%ilignin  =micpar%ilignin
  pltpar%k_woody_litr=micpar%k_woody_litr
  pltpar%k_fine_litr=micpar%k_fine_litr
  !the following variable should be consistent with the soil bgc model
  pltpar%jcplx= micpar%jcplx
  pltpar%NumOfPlantLitrCmplxs=micpar%NumOfPlantLitrCmplxs
  pltpar%jsken  = micpar%jsken
  pltpar%NumLitterGroups= 5     !number of liter groups
  pltpar%MaxNumBranches    = 10    !number of branches
  pltpar%MaxNumRootAxes=10
  pltpar%NumGrowthStages=10
  pltpar%NumOfPlantMorphUnits=7

  MaxNumBranches=pltpar%MaxNumBranches
  NumLitterGroups=pltpar%NumLitterGroups
  MaxNumRootAxes=pltpar%MaxNumRootAxes
  NumOfPlantMorphUnits=pltpar%NumOfPlantMorphUnits
  NumGrowthStages=pltpar%NumGrowthStages
  end subroutine InitPlantMorphSize
end module InitAllocMod
