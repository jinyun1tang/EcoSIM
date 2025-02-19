module PlantAPI
!
! interface to integrate the plant model
  use data_kind_mod   , only : r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only : micpar, pltpar
  use SoilPhysDataType, only : SurfAlbedo_col
  use MiniMathMod, only : AZMAX1
  use EcoSIMSolverPar
  use EcoSIMHistMod
  use SnowDataType
  use TracerIDMod
  use SurfLitterDataType
  use LandSurfDataType
  use SoilPropertyDataType
  use ChemTranspDataType
  use EcoSimSumDataType
  use SoilHeatDataType
  use SOMDataType
  use ClimForcDataType
  use EcoSIMCtrlDataType
  use GridDataType
  use RootDataType
  use SoilWaterDataType
  use CanopyDataType
  use PlantDataRateType
  use PlantTraitDataType
  use CanopyRadDataType
  use FlagDataType
  use EcosimBGCFluxType
  use FertilizerDataType
  use SoilBGCDataType
  use PlantMgmtDataType
  use PlantAPIData
implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: PlantAPISend,PlantAPICanMSend
  public :: PlantAPIRecv,PlantAPICanMRecv

  contains

!------------------------------------------------------------------------------------------

  subroutine PlantAPIRecv(I,J,NY,NX)
  !
  !DESCRIPTION
  !
  use EcoSIMConfig, only : jsken=>jskenc,jcplx => jcplxc
  use PlantAPIData, only : plt_rad
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: NB,NR,NZ,K,L,M,N,I1,NE

  I1=I+1;if(I1>DazCurrYear)I1=1
  NumActivePlants(NY,NX)                              = plt_site%NumActivePlants
  PlantPopu_col(NY,NX)                                = plt_site%PlantPopu_col
  ECO_ER_col(NY,NX)                                   = plt_bgcr%ECO_ER_col
  Eco_NBP_CumYr_col(NY,NX)                            = plt_bgcr%Eco_NBP_CumYr_col
  Air_Heat_Latent_store_col(NY,NX)                    = plt_ew%Air_Heat_Latent_store_col
  Air_Heat_Sens_store_col(NY,NX)                      = plt_ew%Air_Heat_Sens_store_col
  Eco_AutoR_CumYr_col(NY,NX)                          = plt_bgcr%Eco_AutoR_CumYr_col
  LitrFallStrutElms_col(1:NumPlantChemElms,NY,NX)     = plt_bgcr%LitrFallStrutElms_col(1:NumPlantChemElms)
  EcoHavstElmnt_CumYr_col(1:NumPlantChemElms,NY,NX)   = plt_distb%EcoHavstElmnt_CumYr_col(1:NumPlantChemElms)
  WatHeldOnCanopy_col(NY,NX)                          = plt_ew%WatHeldOnCanopy_col
  Eco_Heat_Sens_col(NY,NX)                            = plt_ew%Eco_Heat_Sens_col
  StandingDeadStrutElms_col(1:NumPlantChemElms,NY,NX) = plt_biom%StandingDeadStrutElms_col(1:NumPlantChemElms)
  H2OLoss_CumYr_col(NY,NX)                            = plt_ew%H2OLoss_CumYr_col
  StemArea_col(NY,NX)                              = plt_morph%StemArea_col
  HeatCanopy2Dist_col(NY,NX)                       = plt_ew%HeatCanopy2Dist_col
  HeatCanopy2Dist_col(NY,NX)                       = plt_ew%HeatCanopy2Dist_col
  CanopyLeafArea_col(NY,NX)                        = plt_morph%CanopyLeafArea_col
  Eco_NetRad_col(NY,NX)                            = plt_rad%Eco_NetRad_col
  Eco_Heat_Latent_col(NY,NX)                       = plt_ew%Eco_Heat_Latent_col
  Eco_Heat_GrndSurf_col(NY,NX)                     = plt_ew%Eco_Heat_GrndSurf_col
  QVegET_col(NY,NX)                                  = plt_ew%QVegET_col
  LWRadCanG(NY,NX)                                 = plt_ew%LWRadCanG
  VapXAir2Canopy_col(NY,NX)                        = plt_ew%VapXAir2Canopy_col
  HeatFlx2Canopy_col(NY,NX)                        = plt_ew%HeatFlx2Canopy_col
  CanopyWat_col(NY,NX)                             = plt_ew%CanopyWat_col
  CanopyHeatStor_col(NY,NX)                        = plt_ew%CanopyHeatStor_col
  TRootGasLossDisturb_pft(idg_beg:idg_NH3,NY,NX) = plt_rbgc%TRootGasLossDisturb_pft(idg_beg:idg_NH3)
  Canopy_NEE_col(NY,NX)                            = plt_bgcr%Canopy_NEE_col
  TPlantRootH2OUptake_col(NY,NX)                   = plt_ew%TPlantRootH2OUptake_col
  FERT(17:19,I1,NY,NX) = plt_distb%FERT(17:19)
  FERT(3,I1,NY,NX)                                                       = plt_distb%FERT(3)
  IYTYP(2,I1,NY,NX)                                                      = plt_distb%IYTYP
  FracRootStalkElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)  = plt_allom%FracRootStalkElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)
  FracRootElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)       = plt_allom%FracRootElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)
  FracShootLeafElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)  = plt_allom%FracShootLeafElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)
  FracShootStalkElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs) = plt_allom%FracShootStalkElmAlloc2Litr(1:NumPlantChemElms,1:NumOfPlantLitrCmplxs)
  QH2OLoss_lnds                                                          = plt_site%QH2OLoss_lnds

  DO L=1,NumOfCanopyLayers
    tCanLeafC_cl(L,NY,NX)        = plt_biom%tCanLeafC_cl(L)
    CanopyStemAareZ_col(L,NY,NX) = plt_morph%CanopyStemAareZ_col(L)
    CanopyLeafAareZ_col(L,NY,NX) = plt_morph%CanopyLeafAareZ_col(L)
  ENDDO
  
  DO L=NU(NY,NX),NL(NY,NX)
    
    DO K=1,jcplx
      DO NE=1,NumPlantChemElms
        REcoDOMProd_vr(NE,K,L,NY,NX)=plt_bgcr%REcoDOMProd_vr(NE,K,L)
      ENDDO
    ENDDO
    
    RO2UptkSoilM_vr(1:NPH,L,NY,NX)=plt_rbgc%RO2UptkSoilM_vr(1:NPH,L)    
  ENDDO

  DO L=0,NL(NY,NX)
    REcoH2PO4DmndBand_vr(L,NY,NX)      = plt_bgcr%REcoH2PO4DmndBand_vr(L)
    REcoH1PO4DmndBand_vr(L,NY,NX)      = plt_bgcr%REcoH1PO4DmndBand_vr(L)
    REcoNO3DmndBand_vr(L,NY,NX)        = plt_bgcr%REcoNO3DmndBand_vr(L)
    REcoNH4DmndBand_vr(L,NY,NX)        = plt_bgcr%REcoNH4DmndBand_vr(L)
    REcoH1PO4DmndSoil_vr(L,NY,NX)      = plt_bgcr%REcoH1PO4DmndSoil_vr(L)
    REcoH2PO4DmndSoil_vr(L,NY,NX)      = plt_bgcr%REcoH2PO4DmndSoil_vr(L)
    REcoNO3DmndSoil_vr(L,NY,NX)        = plt_bgcr%REcoNO3DmndSoil_vr(L)
    REcoNH4DmndSoil_vr(L,NY,NX)        = plt_bgcr%REcoNH4DmndSoil_vr(L)
    REcoO2DmndResp_vr(L,NY,NX)         = plt_bgcr%REcoO2DmndResp_vr(L)
    THeatLossRoot2Soil_vr(L,NY,NX)        = plt_ew%THeatLossRoot2Soil_vr(L)
    TPlantRootH2OLoss_vr(L,NY,NX) = plt_ew%TPlantRootH2OLoss_vr(L)
    DO  K=1,micpar%NumOfPlantLitrCmplxs
      DO  M=1,jsken
        DO NE=1,NumPlantChemElms        
          LitrfalStrutElms_vr(NE,M,K,L,NY,NX)=plt_bgcr%LitrfalStrutElms_vr(NE,M,K,L)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO L=1,NL(NY,NX)
    DO NE=1,NumPlantChemElms
      RootMassElm_vr(NE,L,NY,NX)=  sum(plt_biom%RootMassElm_pvr(NE,L,1:NP0(NY,NX)))
    ENDDO
    totRootLenDens_vr(L,NY,NX)                      = plt_morph%totRootLenDens_vr(L)
    trcg_root_vr(idg_beg:idg_NH3,L,NY,NX)         = plt_rbgc%trcg_root_vr(idg_beg:idg_NH3,L)
    trcg_air2root_flx_vr(idg_beg:idg_NH3,L,NY,NX) = plt_rbgc%trcg_air2root_flx_vr(idg_beg:idg_NH3,L)
    tRootCO2Emis2Root_vr(L,NY,NX)                        = plt_bgcr%tRootCO2Emis2Root_vr(L)
    tRO2MicrbUptk_vr(L,NY,NX)                       = plt_bgcr%tRO2MicrbUptk_vr(L)
    
    trcs_plant_uptake_vr(ids_beg:ids_end,L,NY,NX) =plt_rbgc%trcs_plant_uptake_vr(ids_beg:ids_end,L)
    DO  K=1,jcplx
      tRootMycoExud2Soil_vr(1:NumPlantChemElms,K,L,NY,NX)=plt_bgcr%tRootMycoExud2Soil_vr(1:NumPlantChemElms,K,L)
    ENDDO
  ENDDO

  DO NZ=1,NP0(NY,NX)
    Eco_GPP_CumYr_col(NY,NX)                            = Eco_GPP_CumYr_col(NY,NX)+plt_bgcr%GrossCO2Fix_pft(NZ)
    PARTS_brch(1:pltpar%NumOfPlantMorphUnits,1:pltpar%MaxNumBranches,NZ,NY,NX)= &
      plt_morph%PARTS_brch(1:pltpar%NumOfPlantMorphUnits,1:pltpar%MaxNumBranches,NZ)
    QdewCanopy_CumYr_pft(NZ,NY,NX)                              = QdewCanopy_CumYr_pft(NZ,NY,NX)+plt_ew%QdewCanopy_pft(NZ)  
    RootUptk_N_CumYr_pft(NZ,NY,NX)                              = plt_rbgc%RootUptk_N_CumYr_pft(NZ)
    RootUptk_P_CumYr_pft(NZ,NY,NX)                              = plt_rbgc%RootUptk_P_CumYr_pft(NZ)
    RootElms_pft(1:NumPlantChemElms,NZ,NY,NX)                   = plt_biom%RootElms_pft(1:NumPlantChemElms,NZ)
    ElmBalanceCum_pft(1:NumPlantChemElms,NZ,NY,NX)              = plt_site%ElmBalanceCum_pft(1:NumPlantChemElms,NZ)
    CanopyNonstElms_pft(1:NumPlantChemElms,NZ,NY,NX)            = plt_biom%CanopyNonstElms_pft(1:NumPlantChemElms,NZ)
    CanopyNodulElms_pft(1:NumPlantChemElms,NZ,NY,NX)            = plt_biom%CanopyNodulElms_pft(1:NumPlantChemElms,NZ)
    CanopyNodulNonstElms_pft(1:NumPlantChemElms,NZ,NY,NX)       = plt_biom%CanopyNodulNonstElms_pft(1:NumPlantChemElms,NZ)
    CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ,NY,NX)         = plt_biom%CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ)
    LitrfalStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)           = plt_bgcr%LitrfalStrutElms_pft(1:NumPlantChemElms,NZ)
    EcoHavstElmnt_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)        = plt_distb%EcoHavstElmnt_CumYr_pft(1:NumPlantChemElms,NZ)
    NetCumElmntFlx2Plant_pft(1:NumPlantChemElms,NZ,NY,NX)       = plt_pheno%NetCumElmntFlx2Plant_pft(1:NumPlantChemElms,NZ)
    SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX) = plt_bgcr%SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ)
    LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)     = plt_bgcr%LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ)
    EcoHavstElmntCum_pft(1:NumPlantChemElms,NZ,NY,NX)           = plt_distb%EcoHavstElmntCum_pft(1:NumPlantChemElms,NZ)
    PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)         = plt_rbgc%PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ)
    RootMycoExudElms_pft(1:NumPlantChemElms,NZ,NY,NX)           = plt_rbgc%RootMycoExudElms_pft(1:NumPlantChemElms,NZ)
    StandDeadStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)         = plt_biom%StandDeadStrutElms_pft(1:NumPlantChemElms,NZ)
    SeasonalNonstElms_pft(1:NumPlantChemElms,NZ,NY,NX)          = plt_biom%SeasonalNonstElms_pft(1:NumPlantChemElms,NZ)
    ShootStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)             = plt_biom%ShootStrutElms_pft(1:NumPlantChemElms,NZ)
    LeafStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)              = plt_biom%LeafStrutElms_pft(1:NumPlantChemElms,NZ)
    PetoleStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)            = plt_biom%PetoleStrutElms_pft(1:NumPlantChemElms,NZ)
    StalkStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)             = plt_biom%StalkStrutElms_pft(1:NumPlantChemElms,NZ)
    StalkRsrvElms_pft(1:NumPlantChemElms,NZ,NY,NX)              = plt_biom%StalkRsrvElms_pft(1:NumPlantChemElms,NZ)
    HuskStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)              = plt_biom%HuskStrutElms_pft(1:NumPlantChemElms,NZ)
    EarStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)               = plt_biom%EarStrutElms_pft(1:NumPlantChemElms,NZ)
    GrainStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)             = plt_biom%GrainStrutElms_pft(1:NumPlantChemElms,NZ)
    RootStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)              = plt_biom%RootStrutElms_pft(1:NumPlantChemElms,NZ)
    NodulStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)             = plt_biom%NodulStrutElms_pft(1:NumPlantChemElms,NZ)
    PlantRootSoilElmNetX_pft(1:NumPlantChemElms,NZ,NY,NX)       = plt_rbgc%PlantRootSoilElmNetX_pft(1:NumPlantChemElms,NZ)
    CanopyLeafArea_pft(NZ,NY,NX)                                = plt_morph%CanopyLeafArea_pft(NZ)
    CanopyMassC_pft(NZ,NY,NX)           = plt_biom%CanopyMassC_pft(NZ)
    CanopyStemArea_pft(NZ,NY,NX)        = plt_morph%CanopyStemArea_pft(NZ)
    NoduleNonstructCconc_pft(NZ,NY,NX)  = plt_biom%NoduleNonstructCconc_pft(NZ)
    CO2NetFix_pft(NZ,NY,NX)             = plt_bgcr%CO2NetFix_pft(NZ)
    CanopyGasCO2_pft(NZ,NY,NX)          = plt_photo%CanopyGasCO2_pft(NZ)
    LeafIntracellularCO2_pft(NZ,NY,NX)  = plt_photo%LeafIntracellularCO2_pft(NZ)
    aquCO2Intraleaf_pft(NZ,NY,NX)       = plt_photo%aquCO2Intraleaf_pft(NZ)
    ClumpFactor_pft(NZ,NY,NX)           = plt_morph%ClumpFactor_pft(NZ)
    rCNNonstRemob_pft(NZ,NY,NX)         = plt_allom%rCNNonstRemob_pft(NZ)
    rCPNonstRemob_pft(NZ,NY,NX)         = plt_allom%rCPNonstRemob_pft(NZ)
    RootFracRemobilizableBiom(NZ,NY,NX) = plt_allom%RootFracRemobilizableBiom(NZ)
    CNRTS_pft(NZ,NY,NX)                 = plt_allom%CNRTS_pft(NZ)
    CPRTS_pft(NZ,NY,NX)                 = plt_allom%CPRTS_pft(NZ)
    ETCanopy_CumYr_pft(NZ,NY,NX)        = plt_ew%ETCanopy_CumYr_pft(NZ)
    GrossCO2Fix_pft(NZ,NY,NX)           = plt_bgcr%GrossCO2Fix_pft(NZ)
    GrossCO2Fix_CumYr_pft(NZ,NY,NX)     = GrossCO2Fix_CumYr_pft(NZ,NY,NX)+GrossCO2Fix_pft(NZ,NY,NX)
    ChillHours_pft(NZ,NY,NX)            = plt_photo%ChillHours_pft(NZ)
    DiffCO2Atmos2Intracel_pft(NZ,NY,NX) = plt_photo%DiffCO2Atmos2Intracel_pft(NZ)
    DeltaTKC_pft(NZ,NY,NX)              = plt_ew%DeltaTKC_pft(NZ)
    ENGYX_pft(NZ,NY,NX)                 = plt_ew%ENGYX_pft(NZ)
    Transpiration_pft(NZ,NY,NX)         = plt_ew%Transpiration_pft(NZ)
    VapXAir2Canopy_pft(NZ,NY,NX)        = plt_ew%VapXAir2Canopy_pft(NZ)
    EvapTransLHeat_pft(NZ,NY,NX)         = plt_ew%EvapTransLHeat_pft(NZ)
    AirConc_pft(NZ,NY,NX)               = plt_photo%AirConc_pft(NZ)
    FracPARads2Canopy_pft(NZ,NY,NX)    = plt_rad%FracPARads2Canopy_pft(NZ)
    FracGroth2Node_pft(NZ,NY,NX)        = plt_allom%FracGroth2Node_pft(NZ)
    CanopySeedNum_pft(NZ,NY,NX)         = plt_morph%CanopySeedNum_pft(NZ)
    HypoctoHeight_pft(NZ,NY,NX)         = plt_morph%HypoctoHeight_pft(NZ)
    HighTempLimitSeed_pft(NZ,NY,NX)     = plt_pheno%HighTempLimitSeed_pft(NZ)
    HeatStorCanopy_pft(NZ,NY,NX)        = plt_ew%HeatStorCanopy_pft(NZ)
    canopy_growth_pft(NZ,NY,NX)         = plt_rbgc%canopy_growth_pft(NZ)
    CanopyHeight4WatUptake_pft(NZ,NY,NX)      = plt_morph%CanopyHeight4WatUptake_pft(NZ)
    IsPlantActive_pft(NZ,NY,NX)         = plt_pheno%IsPlantActive_pft(NZ)
    iPlantState_pft(NZ,NY,NX)           = plt_pheno%iPlantState_pft(NZ)
    iPlantShootState_pft(NZ,NY,NX)      = plt_pheno%iPlantShootState_pft(NZ)
    iPlantRootState_pft(NZ,NY,NX)       = plt_pheno%iPlantRootState_pft(NZ)
    iDayPlanting_pft(NZ,NY,NX)          = plt_distb%iDayPlanting_pft(NZ)
    iYearPlanting_pft(NZ,NY,NX)         = plt_distb%iYearPlanting_pft(NZ)
    doInitPlant_pft(NZ,NY,NX)           = plt_pheno%doInitPlant_pft(NZ)
    iDayPlantHarvest_pft(NZ,NY,NX)      = plt_distb%iDayPlantHarvest_pft(NZ)
    iYearPlantHarvest_pft(NZ,NY,NX)     = plt_distb%iYearPlantHarvest_pft(NZ)
    NIXBotRootLayer_pft(NZ,NY,NX)       = plt_morph%NIXBotRootLayer_pft(NZ)
    BranchNumber_pft(NZ,NY,NX)          = plt_morph%BranchNumber_pft(NZ)
    NumOfBranches_pft(NZ,NY,NX)         = plt_morph%NumOfBranches_pft(NZ)
    NumRootAxes_pft(NZ,NY,NX)           = plt_morph%NumRootAxes_pft(NZ)
    MaxSoiL4Root_pft(NZ,NY,NX)          = plt_morph%MaxSoiL4Root_pft(NZ)
    NGTopRootLayer_pft(NZ,NY,NX)        = plt_morph%NGTopRootLayer_pft(NZ)
    MainBranchNum_pft(NZ,NY,NX)         = plt_morph%MainBranchNum_pft(NZ)
    NumCogrowthNode_pft(NZ,NY,NX)        = plt_morph%NumCogrowthNode_pft(NZ)
    O2L(NZ,NY,NX)                                      = plt_photo%O2L(NZ)
    O2I(NZ,NY,NX)                                      = plt_photo%O2I(NZ)
    TempOffset_pft(NZ,NY,NX)                           = plt_pheno%TempOffset_pft(NZ)
    PlantO2Stress_pft(NZ,NY,NX)                        = plt_pheno%PlantO2Stress_pft(NZ)
    PPX_pft(NZ,NY,NX)                                  = plt_site%PPX_pft(NZ)
    PlantPopulation_pft(NZ,NY,NX)                      = plt_site%PlantPopulation_pft(NZ)
    PPI_pft(NZ,NY,NX)                                  = plt_site%PPI_pft(NZ)
    PSICanopy_pft(NZ,NY,NX)                            = plt_ew%PSICanopy_pft(NZ)
    PSICanopyOsmo_pft(NZ,NY,NX)                        = plt_ew%PSICanopyOsmo_pft(NZ)
    PSICanopyTurg_pft(NZ,NY,NX)                        = plt_ew%PSICanopyTurg_pft(NZ)
    PSICanPDailyMin(NZ,NY,NX)                          = plt_ew%PSICanPDailyMin(NZ)
    CO2FixCL_pft(NZ,NY,NX)                             = plt_rbgc%CO2FixCL_pft(NZ)
    CO2FixLL_pft(NZ,NY,NX)                             = plt_rbgc%CO2FixLL_pft(NZ)
    RootGasLossDisturb_pft(idg_beg:idg_NH3,NZ,NY,NX) = plt_bgcr%RootGasLossDisturb_pft(idg_beg:idg_NH3,NZ)
    MinCanPStomaResistH2O_pft(NZ,NY,NX)                = plt_photo%MinCanPStomaResistH2O_pft(NZ)
    H2OCuticleResist_pft(NZ,NY,NX)                     = plt_photo%H2OCuticleResist_pft(NZ)
    CO2CuticleResist_pft(NZ,NY,NX)                     = plt_photo%CO2CuticleResist_pft(NZ)
    NH3Dep2Can_pft(NZ,NY,NX)                           = plt_bgcr%NH3Dep2Can_pft(NZ)
    RadNet2Canopy_pft(NZ,NY,NX)                        = plt_rad%RadNet2Canopy_pft(NZ)
    ReistanceCanopy_pft(NZ,NY,NX)                                      = plt_ew%ReistanceCanopy_pft(NZ)
    CanPStomaResistH2O_pft(NZ,NY,NX)                   = plt_photo%CanPStomaResistH2O_pft(NZ)
    CanopyBndlResist_pft(NZ,NY,NX)                     = plt_photo%CanopyBndlResist_pft(NZ)
    PlantinDepz_pft(NZ,NY,NX)                          = plt_morph%PlantinDepz_pft(NZ)
    CO2Solubility_pft(NZ,NY,NX)                        = plt_photo%CO2Solubility_pft(NZ)
    LeafO2Solubility_pft(NZ,NY,NX)                     = plt_photo%LeafO2Solubility_pft(NZ)
    SeedTempSens_pft(NZ,NY,NX)                         = plt_pheno%SeedTempSens_pft(NZ)
    SeedVolumeMean_pft(NZ,NY,NX)                       = plt_morph%SeedVolumeMean_pft(NZ)
    SeedMeanLen_pft(NZ,NY,NX)                          = plt_morph%SeedMeanLen_pft(NZ)
    SeedAreaMean_pft(NZ,NY,NX)                         = plt_morph%SeedAreaMean_pft(NZ)
    SeedDepth_pft(NZ,NY,NX)                            = plt_morph%SeedDepth_pft(NZ)
    HeatXAir2PCan_pft(NZ,NY,NX)                            = plt_ew%HeatXAir2PCan_pft(NZ)
    GrossResp_pft(NZ,NY,NX)                            = plt_bgcr%GrossResp_pft(NZ)
    GrossRespC_CumYr_pft(NZ,NY,NX)                     = GrossRespC_CumYr_pft(NZ,NY,NX)+GrossResp_pft(NZ,NY,NX)
    CanopyRespC_CumYr_pft(NZ,NY,NX)                    = plt_bgcr%CanopyRespC_CumYr_pft(NZ)
    TC4LeafOut_pft(NZ,NY,NX)                           = plt_pheno%TC4LeafOut_pft(NZ)
    TC4LeafOff_pft(NZ,NY,NX)                           = plt_pheno%TC4LeafOff_pft(NZ)
    NH3Emis_CumYr_pft(NZ,NY,NX)                        = plt_bgcr%NH3Emis_CumYr_pft(NZ)
    NodulInfectElms_pft(1:NumPlantChemElms,NZ,NY,NX) = plt_bgcr%NodulInfectElms_pft(1:NumPlantChemElms,NZ)
    NodulInfectElmsCum_pft(1:NumPlantChemElms,NZ,NY,NX)=plt_bgcr%NodulInfectElmsCum_pft(1:NumPlantChemElms,NZ)
    PlantN2Fix_CumYr_pft(NZ,NY,NX)             = plt_bgcr%PlantN2Fix_CumYr_pft(NZ)
    TKC_pft(NZ,NY,NX)                              = plt_ew%TKC_pft(NZ)

    TdegCCanopy_pft(NZ,NY,NX)               = plt_ew%TdegCCanopy_pft(NZ)
    LWRadCanopy_pft(NZ,NY,NX)                  = plt_rad%LWRadCanopy_pft(NZ)
    TKCanopy_pft(NZ,NY,NX)                     = plt_ew%TKCanopy_pft(NZ)
    TKGroth_pft(NZ,NY,NX)                      = plt_pheno%TKGroth_pft(NZ)
    TCGroth_pft(NZ,NY,NX)                      = plt_pheno%TCGroth_pft(NZ)
    fTCanopyGroth_pft(NZ,NY,NX)                = plt_pheno%fTCanopyGroth_pft(NZ)
    RootN2Fix_pft(NZ,NY,NX)                    = plt_rbgc%RootN2Fix_pft(NZ)
    RootNH4Uptake_pft(NZ,NY,NX)                = plt_rbgc%RootNH4Uptake_pft(NZ)
    RootNO3Uptake_pft(NZ,NY,NX)                = plt_rbgc%RootNO3Uptake_pft(NZ)
    RootH2PO4Uptake_pft(NZ,NY,NX)              = plt_rbgc%RootH2PO4Uptake_pft(NZ)
    RootHPO4Uptake_pft(NZ,NY,NX)               = plt_rbgc%RootHPO4Uptake_pft(NZ)
    WatHeldOnCanopy_pft(NZ,NY,NX)              = plt_ew%WatHeldOnCanopy_pft(NZ)
    CanopyBiomWater_pft(NZ,NY,NX)                  = plt_ew%CanopyBiomWater_pft(NZ)
    CO2ByFire_CumYr_pft(NZ,NY,NX)              = plt_distb%CO2ByFire_CumYr_pft(NZ)
    CH4ByFire_CumYr_pft(NZ,NY,NX)              = plt_distb%CH4ByFire_CumYr_pft(NZ)
    O2ByFire_CumYr_pft(NZ,NY,NX)               = plt_distb%O2ByFire_CumYr_pft(NZ)
    NH3byFire_CumYr_pft(NZ,NY,NX)              = plt_distb%NH3byFire_CumYr_pft(NZ)
    N2ObyFire_CumYr_pft(NZ,NY,NX)              = plt_distb%N2ObyFire_CumYr_pft(NZ)
    PO4byFire_CumYr_pft(NZ,NY,NX)              = plt_distb%PO4byFire_CumYr_pft(NZ)
    VHeatCapCanopy_pft(NZ,NY,NX)                 = plt_ew%VHeatCapCanopy_pft(NZ)
    HoursTooLowPsiCan_pft(NZ,NY,NX)            = plt_pheno%HoursTooLowPsiCan_pft(NZ)
    SeedCPlanted_pft(NZ,NY,NX)                 = plt_biom%SeedCPlanted_pft(NZ)
    CanopyStalkC_pft(NZ,NY,NX)                 = plt_biom%CanopyStalkC_pft(NZ)
    CanopyLeafShethC_pft(NZ,NY,NX)             = plt_biom%CanopyLeafShethC_pft(NZ)
    RootBiomCPerPlant_pft(NZ,NY,NX)            = plt_biom%RootBiomCPerPlant_pft(NZ)
    Km4LeafaqCO2_pft(NZ,NY,NX)                 = plt_photo%Km4LeafaqCO2_pft(NZ)
    Km4RubiscoCarboxy_pft(NZ,NY,NX)            = plt_photo%Km4RubiscoCarboxy_pft(NZ)
    CanopyHeight_pft(NZ,NY,NX)                 = plt_morph%CanopyHeight_pft(NZ)
    NetPrimProduct_pft(NZ,NY,NX)               = plt_bgcr%NetPrimProduct_pft(NZ)
    ZERO4Groth_pft(NZ,NY,NX)                   = plt_biom%ZERO4Groth_pft(NZ)
    ZERO4Uptk_pft(NZ,NY,NX)                    = plt_rbgc%ZERO4Uptk_pft(NZ)
    ZERO4LeafVar_pft(NZ,NY,NX)                 = plt_biom%ZERO4LeafVar_pft(NZ)
    iPlantThermoAdaptZone_pft(NZ,NY,NX)        = plt_pheno%iPlantThermoAdaptZone_pft(NZ)
    FracCanopyHeightCut_pft(NZ,I,NY,NX)        = plt_distb%FracCanopyHeightCut_pft(NZ)
    iHarvstType_pft(NZ,I,NY,NX)                = plt_distb%iHarvstType_pft(NZ)
    jHarvst_pft(NZ,I,NY,NX)                    = plt_distb%jHarvst_pft(NZ)
    THIN_pft(NZ,I,NY,NX)                       = plt_distb%THIN_pft(NZ)
    ShootElms_pft(1:NumPlantChemElms,NZ,NY,NX) = plt_biom%ShootElms_pft(1:NumPlantChemElms,NZ)
    RootElms_pft(1:NumPlantChemElms,NZ,NY,NX)  = plt_biom%RootElms_pft(1:NumPlantChemElms,NZ)
    CanopyGrosRCO2_pft(NZ,NY,NX)               = plt_bgcr%CanopyGrosRCO2_pft(NZ)
    DO L=1,NL(NY,NX)
      RootNodulStrutElms_rpvr(1:NumPlantChemElms,L,NZ,NY,NX) = plt_biom%RootNodulStrutElms_rpvr(1:NumPlantChemElms,L,NZ)
      RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ,NY,NX) = plt_biom%RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ)
      RootN2Fix_pvr(L,NZ,NY,NX)                             = plt_bgcr%RootN2Fix_pvr(L,NZ)
      fTgrowRootP_vr(L,NZ,NY,NX)                            = plt_pheno%fTgrowRootP_vr(L,NZ)
    ENDDO
    DO L=1,NumOfCanopyLayers
      CanopyLeafAreaZ_pft(L,NZ,NY,NX) = plt_morph%CanopyLeafAreaZ_pft(L,NZ)
      CanopyLeafCLyr_pft(L,NZ,NY,NX)  = plt_biom%CanopyLeafCLyr_pft(L,NZ)
      CanopyStemAreaZ_pft(L,NZ,NY,NX) = plt_morph%CanopyStemAreaZ_pft(L,NZ)
    ENDDO

    DO L=0,NL(NY,NX)
      DO K=1,micpar%NumOfPlantLitrCmplxs
        DO M=1,jsken
          LitrfalStrutElms_pvr(1:NumPlantChemElms,M,K,L,NZ,NY,NX)=plt_bgcr%LitrfalStrutElms_pvr(1:NumPlantChemElms,M,K,L,NZ)
        ENDDO
      ENDDO
    ENDDO

    
    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO NE=1,NumPlantChemElms
        CanopyNonstElms_brch(NE,NB,NZ,NY,NX) = plt_biom%CanopyNonstElms_brch(NE,NB,NZ)
        EarStrutElms_brch(NE,NB,NZ,NY,NX)    = plt_biom%EarStrutElms_brch(NE,NB,NZ)
      ENDDO
      ShootC4NonstC_brch(NB,NZ,NY,NX)=plt_biom%ShootC4NonstC_brch(NB,NZ)      
    ENDDO


    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO NE=1,NumPlantChemElms
        CanopyNodulNonstElms_brch(NE,NB,NZ,NY,NX) = plt_biom%CanopyNodulNonstElms_brch(NE,NB,NZ)
        ShootStrutElms_brch(NE,NB,NZ,NY,NX)       = plt_biom%ShootStrutElms_brch(NE,NB,NZ)
        PetoleStrutElms_brch(NE,NB,NZ,NY,NX)      = plt_biom%PetoleStrutElms_brch(NE,NB,NZ)
        StalkStrutElms_brch(NE,NB,NZ,NY,NX)       = plt_biom%StalkStrutElms_brch(NE,NB,NZ)
        LeafPetoNonstElmConc_brch(NE,NB,NZ,NY,NX) = plt_biom%LeafPetoNonstElmConc_brch(NE,NB,NZ)
        LeafStrutElms_brch(NE,NB,NZ,NY,NX)        = plt_biom%LeafStrutElms_brch(NE,NB,NZ)
        StalkRsrvElms_brch(NE,NB,NZ,NY,NX)        = plt_biom%StalkRsrvElms_brch(NE,NB,NZ)
        HuskStrutElms_brch(NE,NB,NZ,NY,NX)        = plt_biom%HuskStrutElms_brch(NE,NB,NZ)
        GrainStrutElms_brch(NE,NB,NZ,NY,NX)       = plt_biom%GrainStrutElms_brch(NE,NB,NZ)
        CanopyNodulStrutElms_brch(NE,NB,NZ,NY,NX) = plt_biom%CanopyNodulStrutElms_brch(NE,NB,NZ)
        PetioleChemElmRemob_brch(NE,NB,NZ,NY,NX)  = plt_biom%PetioleChemElmRemob_brch(NE,NB,NZ)
      ENDDO
    ENDDO
    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO L=1,NumOfCanopyLayers
        CanopyStalkArea_lbrch(L,NB,NZ,NY,NX)=plt_morph%CanopyStalkArea_lbrch(L,NB,NZ)
      ENDDO
      Hours2LeafOut_brch(NB,NZ,NY,NX) = plt_pheno%Hours2LeafOut_brch(NB,NZ)
      LeafAreaLive_brch(NB,NZ,NY,NX)  = plt_morph%LeafAreaLive_brch(NB,NZ)
      LeafAreaDying_brch(NB,NZ,NY,NX) = plt_morph%LeafAreaDying_brch(NB,NZ)

      HourFailGrainFill_brch(NB,NZ,NY,NX)                         = plt_pheno%HourFailGrainFill_brch(NB,NZ)
      RubiscoActivity_brch(NB,NZ,NY,NX)                           = plt_photo%RubiscoActivity_brch(NB,NZ)
      C4PhotosynDowreg_brch(NB,NZ,NY,NX)                          = plt_photo%C4PhotosynDowreg_brch(NB,NZ)
      HoursDoingRemob_brch(NB,NZ,NY,NX)                           = plt_pheno%HoursDoingRemob_brch(NB,NZ)
      MatureGroup_brch(NB,NZ,NY,NX)                               = plt_pheno%MatureGroup_brch(NB,NZ)
      NodeNumNormByMatgrp_brch(NB,NZ,NY,NX)                       = plt_pheno%NodeNumNormByMatgrp_brch(NB,NZ)
      ReprodNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX)                = plt_pheno%ReprodNodeNumNormByMatrgrp_brch(NB,NZ)
      PotentialSeedSites_brch(NB,NZ,NY,NX)                        = plt_morph%PotentialSeedSites_brch(NB,NZ)
      SeedNumSet_brch(NB,NZ,NY,NX)                                = plt_morph%SeedNumSet_brch(NB,NZ)
      GrainSeedBiomCMean_brch(NB,NZ,NY,NX)                        = plt_allom%GrainSeedBiomCMean_brch(NB,NZ)
      CanPBranchHeight(NB,NZ,NY,NX)                               = plt_morph%CanPBranchHeight(NB,NZ)
      doRemobilization_brch(NB,NZ,NY,NX)                          = plt_pheno%doRemobilization_brch(NB,NZ)
      iPlantBranchState_brch(NB,NZ,NY,NX)                         = plt_pheno%iPlantBranchState_brch(NB,NZ)
      doPlantLeaveOff_brch(NB,NZ,NY,NX)                           = plt_pheno%doPlantLeaveOff_brch(NB,NZ)
      doPlantLeafOut_brch(NB,NZ,NY,NX)                            = plt_pheno%doPlantLeafOut_brch(NB,NZ)
      doInitLeafOut_brch(NB,NZ,NY,NX)                             = plt_pheno%doInitLeafOut_brch(NB,NZ)
      doSenescence_brch(NB,NZ,NY,NX)                              = plt_pheno%doSenescence_brch(NB,NZ)
      Prep4Literfall_brch(NB,NZ,NY,NX)                            = plt_pheno%Prep4Literfall_brch(NB,NZ)
      Hours4LiterfalAftMature_brch(NB,NZ,NY,NX)                   = plt_pheno%Hours4LiterfalAftMature_brch(NB,NZ)
      KHiestGroLeafNode_brch(NB,NZ,NY,NX)                         = plt_pheno%KHiestGroLeafNode_brch(NB,NZ)
      KLeafNumber_brch(NB,NZ,NY,NX)                               = plt_morph%KLeafNumber_brch(NB,NZ)
      KMinNumLeaf4GroAlloc_brch(NB,NZ,NY,NX)                      = plt_morph%KMinNumLeaf4GroAlloc_brch(NB,NZ)
      KLowestGroLeafNode_brch(NB,NZ,NY,NX)                        = plt_pheno%KLowestGroLeafNode_brch(NB,NZ)
      BranchNumber_brch(NB,NZ,NY,NX)                              = plt_morph%BranchNumber_brch(NB,NZ)
      ShootNodeNum_brch(NB,NZ,NY,NX)                              = plt_morph%ShootNodeNum_brch(NB,NZ)
      NodeNum2InitFloral_brch(NB,NZ,NY,NX)                        = plt_morph%NodeNum2InitFloral_brch(NB,NZ)
      NodeNumberAtAnthesis_brch(NB,NZ,NY,NX)                      = plt_morph%NodeNumberAtAnthesis_brch(NB,NZ)
      LeafElmntRemobFlx_brch(1:NumPlantChemElms,NB,NZ,NY,NX)      = plt_pheno%LeafElmntRemobFlx_brch(1:NumPlantChemElms,NB,NZ)
      PetioleChemElmRemobFlx_brch(1:NumPlantChemElms,NB,NZ,NY,NX) = plt_pheno%PetioleChemElmRemobFlx_brch(1:NumPlantChemElms,NB,NZ)

      TotalNodeNumNormByMatgrp_brch(NB,NZ,NY,NX)               = plt_pheno%TotalNodeNumNormByMatgrp_brch(NB,NZ)
      TotReproNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX)           = plt_pheno%TotReproNodeNumNormByMatrgrp_brch(NB,NZ)
      Hours4LenthenPhotoPeriod_brch(NB,NZ,NY,NX)               = plt_pheno%Hours4LenthenPhotoPeriod_brch(NB,NZ)
      Hours4ShortenPhotoPeriod_brch(NB,NZ,NY,NX)               = plt_pheno%Hours4ShortenPhotoPeriod_brch(NB,NZ)
      Hours4Leafout_brch(NB,NZ,NY,NX)                          = plt_pheno%Hours4Leafout_brch(NB,NZ)
      Hours4LeafOff_brch(NB,NZ,NY,NX)                          = plt_pheno%Hours4LeafOff_brch(NB,NZ)
      NumOfLeaves_brch(NB,NZ,NY,NX)                            = plt_morph%NumOfLeaves_brch(NB,NZ)
      LeafNumberAtFloralInit_brch(NB,NZ,NY,NX)                 = plt_pheno%LeafNumberAtFloralInit_brch(NB,NZ)
      LeafPetolBiomassC_brch(NB,NZ,NY,NX)                      = plt_biom%LeafPetolBiomassC_brch(NB,NZ)
      dReproNodeNumNormByMatG_brch(NB,NZ,NY,NX)                = plt_pheno%dReproNodeNumNormByMatG_brch(NB,NZ)
      LeafChemElmRemob_brch(1:NumPlantChemElms,NB,NZ,NY,NX)    = plt_biom%LeafChemElmRemob_brch(1:NumPlantChemElms,NB,NZ)
      SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ,NY,NX) = plt_biom%SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)
      StalkLiveBiomassC_brch(NB,NZ,NY,NX)                          = plt_biom%StalkLiveBiomassC_brch(NB,NZ)

      DO K=0,MaxNodesPerBranch
        LeafNodeArea_brch(K,NB,NZ,NY,NX)                          = plt_morph%LeafNodeArea_brch(K,NB,NZ)
        InternodeHeightDead_brch(K,NB,NZ,NY,NX)                  = plt_morph%InternodeHeightDead_brch(K,NB,NZ)
        LiveInterNodeHight_brch(K,NB,NZ,NY,NX)                    = plt_morph%LiveInterNodeHight_brch(K,NB,NZ)
        PetoleLensNode_brch(K,NB,NZ,NY,NX)                        = plt_morph%PetoleLensNode_brch(K,NB,NZ)
        InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX) = plt_biom%InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ)
        LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)      = plt_biom%LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)
        LeafProteinCNode_brch(K,NB,NZ,NY,NX)                      = plt_biom%LeafProteinCNode_brch(K,NB,NZ)
        PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)   = plt_biom%PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)
        PetoleProteinCNode_brch(K,NB,NZ,NY,NX)                    = plt_biom%PetoleProteinCNode_brch(K,NB,NZ)
      ENDDO
      DO  L=1,NumOfCanopyLayers
        DO N=1,NumOfLeafZenithSectors
          StemAreaZsec_brch(N,L,NB,NZ,NY,NX)=plt_morph%StemAreaZsec_brch(N,L,NB,NZ)
        ENDDO
      ENDDO
      DO K=0,MaxNodesPerBranch
        DO  L=1,NumOfCanopyLayers
          CanopyLeafArea_lpft(L,K,NB,NZ,NY,NX)                         = plt_morph%CanopyLeafArea_lpft(L,K,NB,NZ)
          LeafElmsByLayerNode_brch(1:NumPlantChemElms,L,K,NB,NZ,NY,NX) = plt_biom%LeafElmsByLayerNode_brch(1:NumPlantChemElms,L,K,NB,NZ)
        ENDDO
      ENDDO
      DO M=1,pltpar%NumGrowthStages
        iPlantCalendar_brch(M,NB,NZ,NY,NX)=plt_pheno%iPlantCalendar_brch(M,NB,NZ)
      ENDDO
      DO K=1,MaxNodesPerBranch
        DO  L=1,NumOfCanopyLayers
          DO N=1,NumOfLeafZenithSectors
            LeafAreaZsec_brch(N,L,K,NB,NZ,NY,NX)  = plt_morph%LeafAreaZsec_brch(N,L,K,NB,NZ)
            LeafAUnshaded_zsec(N,L,K,NB,NZ,NY,NX) = plt_photo%LeafAUnshaded_zsec(N,L,K,NB,NZ)
          ENDDO
        ENDDO

        CPOOL3_node(K,NB,NZ,NY,NX)                   = plt_photo%CPOOL3_node(K,NB,NZ)
        CPOOL4_node(K,NB,NZ,NY,NX)                   = plt_photo%CPOOL4_node(K,NB,NZ)
        CMassCO2BundleSheath_node(K,NB,NZ,NY,NX)     = plt_photo%CMassCO2BundleSheath_node(K,NB,NZ)
        CO2CompenPoint_node(K,NB,NZ,NY,NX)           = plt_photo%CO2CompenPoint_node(K,NB,NZ)
        RubiscoCarboxyEff_node(K,NB,NZ,NY,NX)        = plt_photo%RubiscoCarboxyEff_node(K,NB,NZ)
        C4CarboxyEff_node(K,NB,NZ,NY,NX)             = plt_photo%C4CarboxyEff_node(K,NB,NZ)
        LigthSatCarboxyRate_node(K,NB,NZ,NY,NX)      = plt_photo%LigthSatCarboxyRate_node(K,NB,NZ)
        LigthSatC4CarboxyRate_node(K,NB,NZ,NY,NX)    = plt_photo%LigthSatC4CarboxyRate_node(K,NB,NZ)
        NutrientCtrlonC4Carboxy_node(K,NB,NZ,NY,NX)  = plt_photo%NutrientCtrlonC4Carboxy_node(K,NB,NZ)
        CMassHCO3BundleSheath_node(K,NB,NZ,NY,NX)    = plt_photo%CMassHCO3BundleSheath_node(K,NB,NZ)
        Vmax4RubiscoCarboxy_pft(K,NB,NZ,NY,NX)       = plt_photo%Vmax4RubiscoCarboxy_pft(K,NB,NZ)
        CO2lmtRubiscoCarboxyRate_node(K,NB,NZ,NY,NX) = plt_photo%CO2lmtRubiscoCarboxyRate_node(K,NB,NZ)
        Vmax4PEPCarboxy_pft(K,NB,NZ,NY,NX)           = plt_photo%Vmax4PEPCarboxy_pft(K,NB,NZ)
        CO2lmtPEPCarboxyRate_node(K,NB,NZ,NY,NX)     = plt_photo%CO2lmtPEPCarboxyRate_node(K,NB,NZ)
      ENDDO
    ENDDO
    DO M=1,jsken
      StandDeadKCompElms_pft(1:NumPlantChemElms,M,NZ,NY,NX)=plt_biom%StandDeadKCompElms_pft(1:NumPlantChemElms,M,NZ)
    ENDDO

    DO  L=1,NL(NY,NX)

      DO K=1,jcplx      
        DOM_vr(idom_doc:idom_dop,K,L,NY,NX)=plt_soilchem%DOM_vr(idom_doc:idom_dop,K,L)
      ENDDO

      DO N=1,pltpar%jroots
        fRootGrowPSISense_pvr(N,L,NZ,NY,NX)                        = plt_pheno%fRootGrowPSISense_pvr(N,L,NZ)
        RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX)    = plt_biom%RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ)
        RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX) = plt_biom%RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ)
        RootProteinConc_rpvr(N,L,NZ,NY,NX)                         = plt_biom%RootProteinConc_rpvr(N,L,NZ)
        trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX)            = plt_rbgc%trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)
        trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX)            = plt_rbgc%trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)
        PSIRoot_pvr(N,L,NZ,NY,NX)                                  = plt_ew%PSIRoot_pvr(N,L,NZ)
        PSIRootOSMO_vr(N,L,NZ,NY,NX)                               = plt_ew%PSIRootOSMO_vr(N,L,NZ)
        PSIRootTurg_vr(N,L,NZ,NY,NX)                               = plt_ew%PSIRootTurg_vr(N,L,NZ)
        Root1stXNumL_pvr(N,L,NZ,NY,NX)                             = plt_morph%Root1stXNumL_pvr(N,L,NZ)
        Root2ndXNum_pvr(N,L,NZ,NY,NX)                              = plt_morph%Root2ndXNum_pvr(N,L,NZ)
        RootLenPerPlant_pvr(N,L,NZ,NY,NX)                          = plt_morph%RootLenPerPlant_pvr(N,L,NZ)
        RootLenDensPerPlant_pvr(N,L,NZ,NY,NX)                      = plt_morph%RootLenDensPerPlant_pvr(N,L,NZ)
        RootPoreVol_pvr(N,L,NZ,NY,NX)                              = plt_morph%RootPoreVol_pvr(N,L,NZ)
        RootVH2O_pvr(N,L,NZ,NY,NX)                                 = plt_morph%RootVH2O_pvr(N,L,NZ)
        Root1stRadius_pvr(N,L,NZ,NY,NX)                            = plt_morph%Root1stRadius_pvr(N,L,NZ)
        Root2ndRadius_pvr(N,L,NZ,NY,NX)                            = plt_morph%Root2ndRadius_pvr(N,L,NZ)
        RootAreaPerPlant_pvr(N,L,NZ,NY,NX)                         = plt_morph%RootAreaPerPlant_pvr(N,L,NZ)
        Root2ndAveLen_pvr(N,L,NZ,NY,NX)                            = plt_morph%Root2ndAveLen_pvr(N,L,NZ)
        RootRespPotent_pvr(N,L,NZ,NY,NX)                           = plt_rbgc%RootRespPotent_pvr(N,L,NZ)
        RootCO2EmisPot_pvr(N,L,NZ,NY,NX)                           = plt_rbgc%RootCO2EmisPot_pvr(N,L,NZ)
        RootCO2Autor_pvr(N,L,NZ,NY,NX)                             = plt_rbgc%RootCO2Autor_pvr(N,L,NZ)
        RootCO2Emis_pvr(N,L,NZ,NY,NX)                              = plt_rbgc%RootCO2Emis_pvr(N,L,NZ)
        RootO2Uptk_pvr(N,L,NZ,NY,NX)                               = plt_rbgc%RootO2Uptk_pvr(N,L,NZ)
        RootGasConductance_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX)     = plt_rbgc%RootGasConductance_pvr(idg_beg:idg_NH3,N,L,NZ)
        RootUptkSoiSol_vr(idg_CO2,N,L,NZ,NY,NX)                    = plt_rbgc%RootUptkSoiSol_vr(idg_CO2,N,L,NZ)
        RootUptkSoiSol_vr(idg_O2,N,L,NZ,NY,NX)         = plt_rbgc%RootUptkSoiSol_vr(idg_O2,N,L,NZ)
        RootUptkSoiSol_vr(idg_CH4,N,L,NZ,NY,NX)        = plt_rbgc%RootUptkSoiSol_vr(idg_CH4,N,L,NZ)
        RootUptkSoiSol_vr(idg_N2O,N,L,NZ,NY,NX)        = plt_rbgc%RootUptkSoiSol_vr(idg_N2O,N,L,NZ)
        RootUptkSoiSol_vr(idg_NH3,N,L,NZ,NY,NX)        = plt_rbgc%RootUptkSoiSol_vr(idg_NH3,N,L,NZ)
        RootUptkSoiSol_vr(idg_NH3B,N,L,NZ,NY,NX)       = plt_rbgc%RootUptkSoiSol_vr(idg_NH3B,N,L,NZ)
        RootUptkSoiSol_vr(idg_H2,N,L,NZ,NY,NX)         = plt_rbgc%RootUptkSoiSol_vr(idg_H2,N,L,NZ)
        trcg_air2root_flx_pvr(idg_CO2,N,L,NZ,NY,NX)    = plt_rbgc%trcg_air2root_flx_pvr(idg_CO2,N,L,NZ)
        trcg_air2root_flx_pvr(idg_O2,N,L,NZ,NY,NX)     = plt_rbgc%trcg_air2root_flx_pvr(idg_O2,N,L,NZ)
        trcg_air2root_flx_pvr(idg_CH4,N,L,NZ,NY,NX)    = plt_rbgc%trcg_air2root_flx_pvr(idg_CH4,N,L,NZ)
        trcg_air2root_flx_pvr(idg_N2O,N,L,NZ,NY,NX)    = plt_rbgc%trcg_air2root_flx_pvr(idg_N2O,N,L,NZ)
        trcg_air2root_flx_pvr(idg_NH3,N,L,NZ,NY,NX)    = plt_rbgc%trcg_air2root_flx_pvr(idg_NH3,N,L,NZ)
        trcg_air2root_flx_pvr(idg_H2,N,L,NZ,NY,NX)     = plt_rbgc%trcg_air2root_flx_pvr(idg_H2,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_CO2,N,L,NZ,NY,NX) = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_CO2,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_O2,N,L,NZ,NY,NX)  = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_O2,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_CH4,N,L,NZ,NY,NX) = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_CH4,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_N2O,N,L,NZ,NY,NX) = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_N2O,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_NH3,N,L,NZ,NY,NX) = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_NH3,N,L,NZ)
        trcg_Root_gas2aqu_flx_vr(idg_H2,N,L,NZ,NY,NX)  = plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_H2,N,L,NZ)
        RootNH4DmndSoil_pvr(N,L,NZ,NY,NX)              = plt_rbgc%RootNH4DmndSoil_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_NH4,N,L,NZ,NY,NX)        = plt_rbgc%RootNutUptake_pvr(ids_NH4,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_NH4,N,L,NZ,NY,NX)    = plt_rbgc%RootOUlmNutUptake_pvr(ids_NH4,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_NH4,N,L,NZ,NY,NX)    = plt_rbgc%RootCUlmNutUptake_pvr(ids_NH4,N,L,NZ)
        RootNH4DmndBand_pvr(N,L,NZ,NY,NX)              = plt_rbgc%RootNH4DmndBand_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_NH4B,N,L,NZ,NY,NX)       = plt_rbgc%RootNutUptake_pvr(ids_NH4B,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_NH4B,N,L,NZ,NY,NX)   = plt_rbgc%RootOUlmNutUptake_pvr(ids_NH4B,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_NH4B,N,L,NZ,NY,NX)   = plt_rbgc%RootCUlmNutUptake_pvr(ids_NH4B,N,L,NZ)
        RootNO3DmndSoil_pvr(N,L,NZ,NY,NX)              = plt_rbgc%RootNO3DmndSoil_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_NO3,N,L,NZ,NY,NX)        = plt_rbgc%RootNutUptake_pvr(ids_NO3,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_NO3,N,L,NZ,NY,NX)    = plt_rbgc%RootOUlmNutUptake_pvr(ids_NO3,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_NO3,N,L,NZ,NY,NX)    = plt_rbgc%RootCUlmNutUptake_pvr(ids_NO3,N,L,NZ)
        RootNO3DmndBand_pvr(N,L,NZ,NY,NX)              = plt_rbgc%RootNO3DmndBand_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_NO3B,N,L,NZ,NY,NX)       = plt_rbgc%RootNutUptake_pvr(ids_NO3B,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_NO3B,N,L,NZ,NY,NX)   = plt_rbgc%RootOUlmNutUptake_pvr(ids_NO3B,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_NO3B,N,L,NZ,NY,NX)   = plt_rbgc%RootCUlmNutUptake_pvr(ids_NO3B,N,L,NZ)
        RootH2PO4DmndSoil_pvr(N,L,NZ,NY,NX)            = plt_rbgc%RootH2PO4DmndSoil_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_H2PO4,N,L,NZ,NY,NX)      = plt_rbgc%RootNutUptake_pvr(ids_H2PO4,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_H2PO4,N,L,NZ,NY,NX)  = plt_rbgc%RootOUlmNutUptake_pvr(ids_H2PO4,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_H2PO4,N,L,NZ,NY,NX)  = plt_rbgc%RootCUlmNutUptake_pvr(ids_H2PO4,N,L,NZ)
        RootH2PO4DmndBand_pvr(N,L,NZ,NY,NX)            = plt_rbgc%RootH2PO4DmndBand_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_H2PO4B,N,L,NZ,NY,NX)     = plt_rbgc%RootNutUptake_pvr(ids_H2PO4B,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ,NY,NX) = plt_rbgc%RootOUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ,NY,NX) = plt_rbgc%RootCUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ)
        RootH1PO4DmndSoil_pvr(N,L,NZ,NY,NX)            = plt_rbgc%RootH1PO4DmndSoil_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_H1PO4,N,L,NZ,NY,NX)      = plt_rbgc%RootNutUptake_pvr(ids_H1PO4,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_H1PO4,N,L,NZ,NY,NX)  = plt_rbgc%RootOUlmNutUptake_pvr(ids_H1PO4,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_H1PO4,N,L,NZ,NY,NX)  = plt_rbgc%RootCUlmNutUptake_pvr(ids_H1PO4,N,L,NZ)
        RootH1PO4DmndBand_pvr(N,L,NZ,NY,NX)            = plt_rbgc%RootH1PO4DmndBand_pvr(N,L,NZ)
        RootNutUptake_pvr(ids_H1PO4B,N,L,NZ,NY,NX)     = plt_rbgc%RootNutUptake_pvr(ids_H1PO4B,N,L,NZ)
        RootOUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ,NY,NX) = plt_rbgc%RootOUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ)
        RootCUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ,NY,NX) = plt_rbgc%RootCUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ)
        RootO2Dmnd4Resp_pvr(N,L,NZ,NY,NX)              = plt_rbgc%RootO2Dmnd4Resp_pvr(N,L,NZ)
        AllPlantRootH2OLoss_vr(N,L,NZ,NY,NX)         = plt_ew%AllPlantRootH2OLoss_vr(N,L,NZ)
        RootMycoActiveBiomC_pvr(N,L,NZ,NY,NX)          = plt_biom%RootMycoActiveBiomC_pvr(N,L,NZ)
        PopuRootMycoC_pvr(N,L,NZ,NY,NX)                = AZMAX1(plt_biom%PopuRootMycoC_pvr(N,L,NZ))
        RootProteinC_pvr(N,L,NZ,NY,NX)                 = plt_biom%RootProteinC_pvr(N,L,NZ)
        RAutoRootO2Limter_rpvr(N,L,NZ,NY,NX)            = plt_rbgc%RAutoRootO2Limter_rpvr(N,L,NZ)
        RootCO2Autor_vr(L,NY,NX)                       = RootCO2Autor_vr(L,NY,NX)+RootCO2Autor_pvr(N,L,NZ,NY,NX)
      ENDDO
    ENDDO

    DO NR=1,NumOfCanopyLayers
      NIXBotRootLayer_rpft(NR,NZ,NY,NX)=plt_morph%NIXBotRootLayer_rpft(NR,NZ)
      DO N=1,MY(NZ,NY,NX)
        RootMyco1stElm_raxs(1:NumPlantChemElms,N,NR,NZ,NY,NX) = plt_biom%RootMyco1stElm_raxs(1:NumPlantChemElms,N,NR,NZ)
        Root1stDepz_pft(N,NR,NZ,NY,NX)                        = plt_morph%Root1stDepz_pft(N,NR,NZ)
      ENDDO
      DO L=1,NL(NY,NX)
        DO N=1,MY(NZ,NY,NX)
          RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ,NY,NX) = plt_biom%RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ)
          RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ,NY,NX) = plt_biom%RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ)
          Root1stLen_rpvr(N,L,NR,NZ,NY,NX)                              = plt_morph%Root1stLen_rpvr(N,L,NR,NZ)
          Root2ndLen_rpvr(N,L,NR,NZ,NY,NX)                               = plt_morph%Root2ndLen_rpvr(N,L,NR,NZ)
          Root2ndXNum_rpvr(N,L,NR,NZ,NY,NX)                             = plt_morph%Root2ndXNum_rpvr(N,L,NR,NZ)
        ENDDO
      ENDDO
    ENDDO

    FracBiomHarvsted(1:2,1:4,NZ,I,NY,NX)=plt_distb%FracBiomHarvsted(1:2,1:4,NZ)

    DO L=NU(NY,NX),MaxSoiL4Root_pft(NZ,NY,NX)
      DO K=1,jcplx
        DO N=1,MY(NZ,NY,NX)
          DO NE=1,NumPlantChemElms
            RootMycoExudEUptk_pvr(NE,N,K,L,NZ,NY,NX)=plt_rbgc%RootMycoExudEUptk_pvr(NE,N,K,L,NZ)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO M=1,jsken
      DO N=0,pltpar%NumLitterGroups
        DO NE=1,NumPlantChemElms        
          ElmAllocmat4Litr(NE,N,M,NZ,NY,NX)=plt_soilchem%ElmAllocmat4Litr(NE,N,M,NZ)
        enddo
      enddo
    ENDDO
    !the following needs to be updated 
    Root1stMaxRadius_pft(2,NZ,NY,NX) = plt_morph%Root1stMaxRadius_pft(2,NZ)
    Root2ndMaxRadius_pft(2,NZ,NY,NX) = plt_morph%Root2ndMaxRadius_pft(2,NZ)
    RootPorosity_pft(2,NZ,NY,NX)     = plt_morph%RootPorosity_pft(2,NZ)
    VmaxNH4Root_pft(2,NZ,NY,NX)      = plt_rbgc%VmaxNH4Root_pft(2,NZ)
    KmNH4Root_pft(2,NZ,NY,NX)        = plt_rbgc%KmNH4Root_pft(2,NZ)
    CMinNH4Root_pft(2,NZ,NY,NX)      = plt_rbgc%CMinNH4Root_pft(2,NZ)
    VmaxNO3Root_pft(2,NZ,NY,NX)      = plt_rbgc%VmaxNO3Root_pft(2,NZ)
    KmNO3Root_pft(2,NZ,NY,NX)        = plt_rbgc%KmNO3Root_pft(2,NZ)
    CminNO3Root_pft(2,NZ,NY,NX)      = plt_rbgc%CminNO3Root_pft(2,NZ)
    VmaxPO4Root_pft(2,NZ,NY,NX)      = plt_rbgc%VmaxPO4Root_pft(2,NZ)
    KmPO4Root_pft(2,NZ,NY,NX)        = plt_rbgc%KmPO4Root_pft(2,NZ)
    CMinPO4Root_pft(2,NZ,NY,NX)      = plt_rbgc%CMinPO4Root_pft(2,NZ)
    RootRadialResist_pft(2,NZ,NY,NX) = plt_morph%RootRadialResist_pft(2,NZ)
    RootAxialResist_pft(2,NZ,NY,NX)  = plt_morph%RootAxialResist_pft(2,NZ)

    DO N=1,MY(NZ,NY,NX)
      RootMycoNonstElms_pft(1:NumPlantChemElms,N,NZ,NY,NX) = plt_biom%RootMycoNonstElms_pft(1:NumPlantChemElms,N,NZ)
      RootPoreTortu4Gas(N,NZ,NY,NX)                        = plt_morph%RootPoreTortu4Gas(N,NZ)
      RootRaidus_rpft(N,NZ,NY,NX)                          = plt_morph%RootRaidus_rpft(N,NZ)
      RootVolPerMassC_pft(N,NZ,NY,NX)                      = plt_morph%RootVolPerMassC_pft(N,NZ)
      Root1stSpecLen_pft(N,NZ,NY,NX)                       = plt_morph%Root1stSpecLen_pft(N,NZ)
      Root2ndSpecLen_pft(N,NZ,NY,NX)                       = plt_morph%Root2ndSpecLen_pft(N,NZ)
      Root1stMaxRadius1_pft(N,NZ,NY,NX)                    = plt_morph%Root1stMaxRadius1_pft(N,NZ)
      Root2ndMaxRadius1_pft(N,NZ,NY,NX)                    = plt_morph%Root2ndMaxRadius1_pft(N,NZ)
      Root1stXSecArea_pft(N,NZ,NY,NX)                      = plt_morph%Root1stXSecArea_pft(N,NZ)
      Root2ndXSecArea_pft(N,NZ,NY,NX)                      = plt_morph%Root2ndXSecArea_pft(N,NZ)
    enDDO
  ENDDO

  end subroutine PlantAPIRecv


!------------------------------------------------------------------------------------------

  subroutine PlantAPISend(I,J,NY,NX)
  !
  !DESCRIPTION
  !Send data to plant model
  use PlantAPIData, only : plt_rad
  use EcoSIMConfig, only : jsken=>jskenc,jcplx => jcplxc
  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer :: K,L,M,N,NB,NZ,NR,I1,NE

  plt_site%DazCurrYear=DazCurrYear
  I1=I+1;if(I1>DazCurrYear)I1=1  
  plt_site%ZERO                       = ZERO
  plt_site%ZERO2                      = ZERO2
  plt_site%ALAT                       = ALAT(NY,NX)
  plt_site%ATCA                       = ATCA(NY,NX)
  plt_morph%LeafStalkArea_col         = LeafStalkArea_col(NY,NX)
  plt_morph%CanopyLeafArea_col        = CanopyLeafArea_col(NY,NX)
  plt_site%ALT                        = ALT(NY,NX)
  plt_site%CCO2EI                     = CCO2EI(NY,NX)
  plt_site%CO2EI                      = CO2EI(NY,NX)
  plt_bgcr%NetCO2Flx2Canopy_col       = NetCO2Flx2Canopy_col(NY,NX)
  plt_site%CO2E                       = CO2E_col(NY,NX)
  plt_site%AtmGasc(idg_beg:idg_NH3) = AtmGasCgperm3(idg_beg:idg_NH3,NY,NX)
  plt_site%DayLenthPrev               = DayLenthPrev_col(NY,NX)
  plt_site%DayLenthCurrent            = DayLensCurr_col(NY,NX)
  plt_ew%SnowDepth                    = SnowDepth_col(NY,NX)
  plt_site%DayLenthMax                = DayLenthMax(NY,NX)
  plt_site%KoppenClimZone             = KoppenClimZone_col(NY,NX)
  plt_site%iYearCurrent               = iYearCurrent
  plt_site%NL                         = NL(NY,NX)
  plt_site%NP0                        = NP0(NY,NX)
  plt_site%MaxNumRootLays             = MaxNumRootLays(NY,NX)
  plt_site%NP                         = NP(NY,NX)
  plt_site%NU                         = NUM(NY,NX)
  plt_site%OXYE                       = OXYE_col(NY,NX)
  plt_ew%AbvCanopyBndlResist_col      = AbvCanopyBndlResist_col(NY,NX)
  plt_ew%RIB                          = RIB_col(NY,NX)
  plt_rad%SineSunInclAnglNxtHour_col  = SineSunInclAnglNxtHour_col(NY,NX)
  plt_rad%SineSunInclAngle_col        = SineSunInclAngle_col(NY,NX)
  plt_ew%TKSnow                       = TKSnow_snvr(1,NY,NX)  !surface layer snow temperature
  plt_ew%TairK                        = TairK_col(NY,NX)
  plt_rad%LWRadGrnd                   = LWRadGrnd(NY,NX)
  plt_rad%LWRadSky_col                = LWRadSky_col(NY,NX)
  plt_ew%VPA                          = VPA_col(NY,NX)
  plt_distb%XCORP                     = XTillCorp_col(NY,NX)
  plt_site%SolarNoonHour_col          = SolarNoonHour_col(NY,NX)
  plt_site%ZEROS2                     = ZEROS2(NY,NX)
  plt_site%ZEROS                      = ZEROS(NY,NX)
  plt_ew%RoughHeight                  = RoughHeight_col(NY,NX)
  plt_morph%CanopyHeight_col          = CanopyHeight_col(NY,NX)
  plt_ew%ZERO4PlantDisplace_col       = ZERO4PlantDisplace_col(NY,NX)
  plt_distb%DCORP                     = DepzCorp_col(I,NY,NX)
  plt_distb%iSoilDisturbType_col      = iSoilDisturbType_col(I,NY,NX)
  plt_morph%CanopyHeightZ_col(0)      = CanopyHeightZ_col(0,NY,NX)
  DO  L=1,NumOfCanopyLayers
    plt_morph%CanopyHeightZ_col(L) = CanopyHeightZ_col(L,NY,NX)
    plt_rad%TAU_DirRadTransm(L)    = TAU_DirRadTransm(L,NY,NX)
    plt_rad%TAU_RadThru(L)         = TAU_RadThru(L,NY,NX)
  ENDDO
  plt_rad%TAU_DirRadTransm(NumOfCanopyLayers+1) = TAU_DirRadTransm(NumOfCanopyLayers+1,NY,NX)
  plt_rad%TAU_RadThru(NumOfCanopyLayers+1)      = TAU_RadThru(NumOfCanopyLayers+1,NY,NX)

  DO N=1,NumOfLeafZenithSectors
    plt_rad%SineLeafAngle(N)=SineLeafAngle(N)
  ENDDO

  DO L=1,NL(NY,NX)
    plt_soilchem%HydroCondMicP4RootUptake_vr(L) = HydroCondMicP4RootUptake_vr(L,NY,NX)
    plt_soilchem%GasDifc_vr(idg_beg:idg_end,L)  = GasDifc_vr(idg_beg:idg_end,L,NY,NX)
    plt_soilchem%SoilResit4RootPentrate_vr(L)   = SoilResit4RootPentrate_vr(L,NY,NX)
    plt_site%CumSoilThickMidL_vr(L)                        = CumSoilThickMidL_vr(L,NY,NX)
  ENDDO
  plt_allom%FracHour4LeafoffRemob(:) = FracHour4LeafoffRemob(:)
  DO L=1,NL(NY,NX)
    plt_soilchem%trcg_gasml_vr(idg_CO2,L) = trcg_gasml_vr(idg_CO2,L,NY,NX)
    plt_soilchem%trcg_gasml_vr(idg_O2,L)  = trcg_gasml_vr(idg_O2,L,NY,NX)
  ENDDO

  DO L=0,NL(NY,NX)
    plt_site%AREA3(L)                                 = AREA(3,L,NY,NX)
    plt_soilchem%SoilBulkDensity_vr(L)                = SoilBulkDensity_vr(L,NY,NX)
    plt_soilchem%trc_solcl_vr(ids_beg:ids_end,L)      = trc_solcl_vr(ids_beg:ids_end,L,NY,NX)
    plt_soilchem%SoluteDifusvty_vr(ids_beg:ids_end,L) = SoluteDifusvty_vr(ids_beg:ids_end,L,NY,NX)
    plt_soilchem%trcg_gascl_vr(idg_beg:idg_NH3,L)     = trcg_gascl_vr(idg_beg:idg_NH3,L,NY,NX)
    plt_soilchem%CSoilOrgM_vr(ielmc,L)                = CSoilOrgM_vr(ielmc,L,NY,NX)
    plt_site%CumSoilThickness_vr(L)                   = CumSoilThickness_vr(L,NY,NX)
    plt_site%FracSoiAsMicP_vr(L)                      = FracSoiAsMicP_vr(L,NY,NX)
    plt_soilchem%trcs_solml_vr(ids_beg:ids_end,L)     = trcs_solml_vr(ids_beg:ids_end,L,NY,NX)
    plt_soilchem%GasSolbility_vr(idg_beg:idg_NH3,L) = GasSolbility_vr(idg_beg:idg_NH3,L,NY,NX)

    plt_ew%ElvAdjstedSoilH2OPSIMPa_vr(L)      = ElvAdjstedSoilH2OPSIMPa_vr(L,NY,NX)
    plt_bgcr%RH2PO4EcoDmndSoilPrev_vr(L) = RH2PO4EcoDmndSoilPrev_vr(L,NY,NX)
    plt_bgcr%RH2PO4EcoDmndBandPrev_vr(L) = RH2PO4EcoDmndBandPrev_vr(L,NY,NX)
    plt_bgcr%RH1PO4EcoDmndSoilPrev_vr(L) = RH1PO4EcoDmndSoilPrev_vr(L,NY,NX)
    plt_bgcr%RH1PO4EcoDmndBandPrev_vr(L) = RH1PO4EcoDmndBandPrev_vr(L,NY,NX)
    plt_bgcr%RNO3EcoDmndSoilPrev_vr(L)   = RNO3EcoDmndSoilPrev_vr(L,NY,NX)
    plt_bgcr%RNH4EcoDmndSoilPrev_vr(L)   = RNH4EcoDmndSoilPrev_vr(L,NY,NX)
    plt_bgcr%RNH4EcoDmndBandPrev_vr(L)   = RNH4EcoDmndBandPrev_vr(L,NY,NX)
    plt_bgcr%RNO3EcoDmndBandPrev_vr(L)   = RNO3EcoDmndBandPrev_vr(L,NY,NX)
    plt_bgcr%RGasFlxPrev_vr(idg_beg:idg_NH3,L) = RGasFlxPrev_vr(idg_beg:idg_NH3,L,NY,NX)
    plt_bgcr%RO2AquaSourcePrev_vr(L)           = RO2AquaSourcePrev_vr(L,NY,NX)
    plt_bgcr%RO2EcoDmndPrev_vr(L)              = RO2EcoDmndPrev_vr(L,NY,NX)
    plt_ew%TKS_vr(L)                           = TKS_vr(L,NY,NX)

    plt_soilchem%THETW_vr(L)               = THETW_vr(L,NY,NX)
    plt_soilchem%SoilWatAirDry_vr(L)               = SoilWatAirDry_vr(L,NY,NX)
    plt_soilchem%TScal4Difsvity_vr(L)      = TScal4Difsvity_vr(L,NY,NX)
    plt_soilchem%VLSoilPoreMicP_vr(L)      = VLSoilPoreMicP_vr(L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_H1PO4B,L) = trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_NO3,L)    = trcs_VLN_vr(ids_NO3,L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_H1PO4,L)  = trcs_VLN_vr(ids_H1PO4,L,NY,NX)
    plt_soilchem%VLSoilMicP_vr(L)          = VLSoilMicP_vr(L,NY,NX)
    plt_soilchem%VLiceMicP_vr(L)           = VLiceMicP_vr(L,NY,NX)
    plt_soilchem%VLWatMicP_vr(L)           = VLWatMicP_vr(L,NY,NX)
    plt_soilchem%VLMicP_vr(L)              = VLMicP_vr(L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_NO3B,L)   = trcs_VLN_vr(ids_NO3B,L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_NH4,L)    = trcs_VLN_vr(ids_NH4,L,NY,NX)
    plt_soilchem%trcs_VLN_vr(ids_NH4B,L)   = trcs_VLN_vr(ids_NH4B,L,NY,NX)

    plt_site%DLYR3(L)     =DLYR_3D(3,L,NY,NX)
    DO K=1,jcplx
      plt_soilchem%FracBulkSOMC_vr(K,L)          = FracBulkSOMC_vr(K,L,NY,NX)
      plt_soilchem%DOM_vr(idom_doc:idom_dop,K,L) = DOM_vr(idom_doc:idom_dop,K,L,NY,NX)
    ENDDO
  ENDDO

  DO NZ=1,NP0(NY,NX)
!plant properties begin
    plt_photo%iPlantPhotosynthesisType(NZ)        = iPlantPhotosynthesisType(NZ,NY,NX)
    plt_pheno%iPlantRootProfile_pft(NZ)           = iPlantRootProfile_pft(NZ,NY,NX)
    plt_pheno%iPlantPhenolPattern_pft(NZ)         = iPlantPhenolPattern_pft(NZ,NY,NX)
    plt_pheno%iPlantDevelopPattern_pft(NZ)        = iPlantDevelopPattern_pft(NZ,NY,NX)
    plt_pheno%iPlantPhenolType_pft(NZ)            = iPlantPhenolType_pft(NZ,NY,NX)
    plt_pheno%iPlantPhotoperiodType_pft(NZ)       = iPlantPhotoperiodType_pft(NZ,NY,NX)
    plt_pheno%iPlantTurnoverPattern_pft(NZ)       = iPlantTurnoverPattern_pft(NZ,NY,NX)
    plt_pheno%PlantInitThermoAdaptZone(NZ)        = PlantInitThermoAdaptZone(NZ,NY,NX)
    plt_morph%iPlantGrainType_pft(NZ)             = iPlantGrainType_pft(NZ,NY,NX)
    plt_morph%iPlantNfixType_pft(NZ)              = iPlantNfixType_pft(NZ,NY,NX)
    plt_morph%MY(NZ)                              = MY(NZ,NY,NX)
    plt_photo%VmaxRubCarboxyRef_pft(NZ)           = VmaxRubCarboxyRef_pft(NZ,NY,NX)
    plt_photo%VmaxRubOxyRef_pft(NZ)               = VmaxRubOxyRef_pft(NZ,NY,NX)
    plt_photo%VmaxPEPCarboxyRef_pft(NZ)           = VmaxPEPCarboxyRef_pft(NZ,NY,NX)
    plt_photo%XKCO2(NZ)                           = XKCO2(NZ,NY,NX)
    plt_photo%XKO2(NZ)                            = XKO2(NZ,NY,NX)
    plt_photo%Km4PEPCarboxy_pft(NZ)               = Km4PEPCarboxy_pft(NZ,NY,NX)
    plt_photo%LeafRuBPConc_pft(NZ)                = LeafRuBPConc_pft(NZ,NY,NX)
    plt_photo%FracLeafProtinAsPEPCarboxyl_pft(NZ) = FracLeafProtinAsPEPCarboxyl_pft(NZ,NY,NX)
    plt_photo%SpecChloryfilAct_pft(NZ)            = SpecChloryfilAct_pft(NZ,NY,NX)
    plt_photo%LeafC3ChlorofilConc_pft(NZ)         = LeafC3ChlorofilConc_pft(NZ,NY,NX)
    plt_photo%LeafC4ChlorofilConc_pft(NZ)         = LeafC4ChlorofilConc_pft(NZ,NY,NX)
    plt_photo%CanPCi2CaRatio(NZ)                  = CanPCi2CaRatio(NZ,NY,NX)

    plt_pheno%RefNodeInitRate_pft(NZ)        = RefNodeInitRate_pft(NZ,NY,NX)
    plt_pheno%RefLeafAppearRate_pft(NZ)      = RefLeafAppearRate_pft(NZ,NY,NX)
    plt_pheno%TCChill4Seed_pft(NZ)           = TCChill4Seed_pft(NZ,NY,NX)
    plt_morph%rLen2WidthLeaf_pft(NZ)         = rLen2WidthLeaf_pft(NZ,NY,NX)
    plt_pheno%MinNonstC2InitBranch_pft(NZ)   = MinNonstC2InitBranch_pft(NZ,NY,NX)
    plt_morph%ShootNodeNumAtPlanting_pft(NZ) = ShootNodeNumAtPlanting_pft(NZ,NY,NX)
    plt_pheno%CriticPhotoPeriod_pft(NZ)      = CriticPhotoPeriod_pft(NZ,NY,NX)
    plt_pheno%PhotoPeriodSens_pft(NZ)        = PhotoPeriodSens_pft(NZ,NY,NX)
    plt_morph%SLA1_pft(NZ)                   = SLA1_pft(NZ,NY,NX)
    plt_morph%PetoLen2Mass_pft(NZ)           = PetoLen2Mass_pft(NZ,NY,NX)
    plt_morph%NodeLenPergC(NZ)               = NodeLenPergC(NZ,NY,NX)
    DO  N=1,NumOfLeafZenithSectors
      plt_morph%CLASS(N,NZ)=CLASS(N,NZ,NY,NX)
    ENDDO
    plt_morph%ClumpFactorInit_pft(NZ)     = ClumpFactorInit_pft(NZ,NY,NX)
    plt_morph%SineBranchAngle_pft(NZ)     = SineBranchAngle_pft(NZ,NY,NX)
    plt_morph%SinePetioleAngle_pft(NZ)    = SinePetioleAngle_pft(NZ,NY,NX)
    plt_morph%MaxPotentSeedNumber_pft(NZ) = MaxPotentSeedNumber_pft(NZ,NY,NX)
    plt_morph%MaxSeedNumPerSite_pft(NZ)   = MaxSeedNumPerSite_pft(NZ,NY,NX)
    plt_morph%MaxSeedCMass_pft(NZ)        = MaxSeedCMass_pft(NZ,NY,NX)
    plt_morph%SeedCMass_pft(NZ)           = SeedCMass_pft(NZ,NY,NX)
    plt_pheno%GrainFillRate25C_pft(NZ)    = GrainFillRate25C_pft(NZ,NY,NX)
    plt_biom%StandingDeadInitC_pft(NZ)    = StandingDeadInitC_pft(NZ,NY,NX)

    !initial root values
    DO N=1,MY(NZ,NY,NX)
      plt_morph%Root1stMaxRadius_pft(N,NZ) = Root1stMaxRadius_pft(N,NZ,NY,NX)
      plt_morph%Root2ndMaxRadius_pft(N,NZ) = Root2ndMaxRadius_pft(N,NZ,NY,NX)
      plt_morph%RootPorosity_pft(N,NZ)     = RootPorosity_pft(N,NZ,NY,NX)
      plt_morph%RootRadialResist_pft(N,NZ) = RootRadialResist_pft(N,NZ,NY,NX)
      plt_morph%RootAxialResist_pft(N,NZ)  = RootAxialResist_pft(N,NZ,NY,NX)
      plt_rbgc%VmaxNH4Root_pft(N,NZ)       = VmaxNH4Root_pft(N,NZ,NY,NX)
      plt_rbgc%KmNH4Root_pft(N,NZ)         = KmNH4Root_pft(N,NZ,NY,NX)
      plt_rbgc%CMinNH4Root_pft(N,NZ)       = CMinNH4Root_pft(N,NZ,NY,NX)
      plt_rbgc%VmaxNO3Root_pft(N,NZ)       = VmaxNO3Root_pft(N,NZ,NY,NX)
      plt_rbgc%KmNO3Root_pft(N,NZ)         = KmNO3Root_pft(N,NZ,NY,NX)
      plt_rbgc%CminNO3Root_pft(N,NZ)       = CminNO3Root_pft(N,NZ,NY,NX)
      plt_rbgc%VmaxPO4Root_pft(N,NZ)       = VmaxPO4Root_pft(N,NZ,NY,NX)
      plt_rbgc%KmPO4Root_pft(N,NZ)         = KmPO4Root_pft(N,NZ,NY,NX)
      plt_rbgc%CMinPO4Root_pft(N,NZ)       = CMinPO4Root_pft(N,NZ,NY,NX)  
    ENDDO
    plt_pheno%MinNonstC2InitRoot_pft(NZ)        = MinNonstC2InitRoot_pft(NZ,NY,NX)
    plt_pheno%ShutRutNonstElmntConducts_pft(NZ) = ShutRutNonstElmntConducts_pft(NZ,NY,NX)
    plt_morph%RootBranchFreq_pft(NZ)            = RootBranchFreq_pft(NZ,NY,NX)
    plt_ew%CanOsmoPsi0pt_pft(NZ)                = CanOsmoPsi0pt_pft(NZ,NY,NX)
    plt_photo%RCS(NZ)                           = RCS(NZ,NY,NX)
    plt_photo%CuticleResist_pft(NZ)             = CuticleResist_pft(NZ,NY,NX)

    plt_allom%LeafBiomGrowthYield(NZ)    = LeafBiomGrowthYield(NZ,NY,NX)
    plt_allom%PetioleBiomGrowthYield(NZ) = PetioleBiomGrowthYield(NZ,NY,NX)
    plt_allom%StalkBiomGrowthYield(NZ)   = StalkBiomGrowthYield(NZ,NY,NX)
    plt_allom%ReserveBiomGrowthYield(NZ) = ReserveBiomGrowthYield(NZ,NY,NX)
    plt_allom%HuskBiomGrowthYield(NZ)    = HuskBiomGrowthYield(NZ,NY,NX)
    plt_allom%EarBiomGrowthYield(NZ)     = EarBiomGrowthYield(NZ,NY,NX)
    plt_allom%GrainBiomGrowthYield(NZ)   = GrainBiomGrowthYield(NZ,NY,NX)
    plt_allom%RootBiomGrosYld_pft(NZ)    = RootBiomGrosYld_pft(NZ,NY,NX)
    plt_allom%NoduGrowthYield_pft(NZ)    = NoduGrowthYield_pft(NZ,NY,NX)
    plt_allom%CNLF(NZ)                   = CNLF(NZ,NY,NX)
    plt_allom%CNSHE(NZ)                  = CNSHE(NZ,NY,NX)
    plt_allom%rNCStalk_pft(NZ)           = rNCStalk_pft(NZ,NY,NX)
    plt_allom%rNCReserve_pft(NZ)         = rNCReserve_pft(NZ,NY,NX)
    plt_allom%rNCHusk_pft(NZ)            = rNCHusk_pft(NZ,NY,NX)
    plt_allom%rNCEar_pft(NZ)             = rNCEar_pft(NZ,NY,NX)
    plt_allom%CNGR(NZ)                   = CNGR(NZ,NY,NX)
    plt_allom%RootrNC_pft(NZ)            = RootrNC_pft(NZ,NY,NX)
    plt_allom%NodulerNC_pft(NZ)          = NodulerNC_pft(NZ,NY,NX)
    plt_allom%CPLF(NZ)                   = CPLF(NZ,NY,NX)
    plt_allom%CPSHE(NZ)                  = CPSHE(NZ,NY,NX)
    plt_allom%rPCStalk_pft(NZ)           = rPCStalk_pft(NZ,NY,NX)
    plt_allom%rPCReserve_pft(NZ)         = rPCReserve_pft(NZ,NY,NX)
    plt_allom%rPCHusk_pft(NZ)            = rPCHusk_pft(NZ,NY,NX)
    plt_allom%rPCEar_pft(NZ)             = rPCEar_pft(NZ,NY,NX)
    plt_allom%CPGR(NZ)                   = CPGR(NZ,NY,NX)
    plt_allom%RootrPC_pft(NZ)            = RootrPC_pft(NZ,NY,NX)
    plt_allom%NodulerPC_pft(NZ)          = NodulerPC_pft(NZ,NY,NX)

!plant properties end

    plt_morph%LeafStalkArea_pft(NZ)   = LeafStalkArea_pft(NZ,NY,NX)
    plt_distb%iPlantingYear_pft(NZ)   = iPlantingYear_pft(NZ,NY,NX)
    plt_distb%iPlantingDay_pft(NZ)    = iPlantingDay_pft(NZ,NY,NX)
    plt_distb%iHarvestYear_pft(NZ)    = iHarvestYear_pft(NZ,NY,NX)
    plt_rad%RadPARbyCanopy_pft(NZ)    = RadPARbyCanopy_pft(NZ,NY,NX)
    plt_rad%RadSWbyCanopy_pft(NZ)     = RadSWbyCanopy_pft(NZ,NY,NX)
    plt_ew%PrecIntcptByCanopy_pft(NZ) = PrecIntcptByCanopy_pft(NZ,NY,NX)

    plt_site%PPatSeeding_pft(NZ)          = PPatSeeding_pft(NZ,NY,NX)
    plt_distb%iHarvestDay_pft(NZ)         = iHarvestDay_pft(NZ,NY,NX)
    plt_morph%ClumpFactorNow_pft(NZ)      = ClumpFactorNow_pft(NZ,NY,NX)
    plt_site%DATAP(NZ)                    = DATAP(NZ,NY,NX)
    plt_pheno%MatureGroup_pft(NZ)         = MatureGroup_pft(NZ,NY,NX)
    plt_biom%AvgCanopyBiomC2Graze_pft(NZ) = AvgCanopyBiomC2Graze_pft(NZ,NY,NX)

    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      plt_pheno%HourReq4LeafOut_brch(NB,NZ)=HourReq4LeafOut_brch(NB,NZ,NY,NX)
      plt_pheno%HourReq4LeafOff_brch(NB,NZ)=HourReq4LeafOff_brch(NB,NZ,NY,NX)
    ENDDO

    DO L=1,NumOfCanopyLayers
      DO  M=1,NumOfSkyAzimuthSects
        DO  N=1,NumOfLeafZenithSectors
          plt_rad%RadPAR_zsec(N,M,L,NZ)    = RadPAR_zsec(N,M,L,NZ,NY,NX)
          plt_rad%RadDifPAR_zsec(N,M,L,NZ) = RadDifPAR_zsec(N,M,L,NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO L=1,NL(NY,NX)
    DO M=1,NPH
      plt_site%VLWatMicPM_vr(M,L)           = VLWatMicPM_vr(M,L,NY,NX)
      plt_site%VLsoiAirPM_vr(M,L)           = VLsoiAirPM_vr(M,L,NY,NX)
      plt_site%TortMicPM_vr(M,L)            = TortMicPM_vr(M,L,NY,NX)
      plt_site%FILMM_vr(M,L)                    = FILMM_vr(M,L,NY,NX)
      plt_soilchem%DiffusivitySolutEffM_vr(M,L) = DiffusivitySolutEffM_vr(M,L,NY,NX)
    ENDDO
  ENDDO

! sent variables also modified
  plt_site%NumActivePlants                               = NumActivePlants(NY,NX)
  plt_site%QH2OLoss_lnds                                 = QH2OLoss_lnds
  plt_site%PlantPopu_col                                 = PlantPopu_col(NY,NX)
  plt_bgcr%ECO_ER_col                                    = ECO_ER_col(NY,NX)
  plt_biom%StandingDeadStrutElms_col(1:NumPlantChemElms) = StandingDeadStrutElms_col(1:NumPlantChemElms,NY,NX)
  plt_bgcr%LitrFallStrutElms_col(1:NumPlantChemElms)     = LitrFallStrutElms_col(1:NumPlantChemElms,NY,NX)
  plt_morph%StemArea_col                                 = StemArea_col(NY,NX)
  plt_ew%Eco_Heat_Sens_col                               = Eco_Heat_Sens_col(NY,NX)
  plt_ew%WatHeldOnCanopy_col                             = WatHeldOnCanopy_col(NY,NX)
  plt_bgcr%Eco_NBP_CumYr_col                             = Eco_NBP_CumYr_col(NY,NX)
  plt_ew%Air_Heat_Latent_store_col                       = Air_Heat_Latent_store_col(NY,NX)
  plt_ew%Air_Heat_Sens_store_col                         = Air_Heat_Sens_store_col(NY,NX)
  plt_bgcr%Eco_AutoR_CumYr_col                           = Eco_AutoR_CumYr_col(NY,NX)
  plt_ew%H2OLoss_CumYr_col                               = H2OLoss_CumYr_col(NY,NX)
  plt_distb%EcoHavstElmnt_CumYr_col(1:NumPlantChemElms)  = EcoHavstElmnt_CumYr_col(1:NumPlantChemElms,NY,NX)
  plt_rad%Eco_NetRad_col                                 = Eco_NetRad_col(NY,NX)
  plt_ew%VapXAir2Canopy_col                              = VapXAir2Canopy_col(NY,NX)
  plt_ew%Eco_Heat_Latent_col                             = Eco_Heat_Latent_col(NY,NX)
  plt_rbgc%TRootGasLossDisturb_pft(idg_beg:idg_NH3)    = TRootGasLossDisturb_pft(idg_beg:idg_NH3,NY,NX)
  plt_ew%Eco_Heat_GrndSurf_col                           = Eco_Heat_GrndSurf_col(NY,NX)
  plt_ew%QVegET_col                                        = QVegET_col(NY,NX)
  plt_ew%HeatFlx2Canopy_col                              = HeatFlx2Canopy_col(NY,NX)
  plt_ew%LWRadCanG                                       = LWRadCanG(NY,NX)
  plt_ew%CanopyWat_col                                   = CanopyWat_col(NY,NX)
  plt_ew%CanopyHeatStor_col                              = CanopyHeatStor_col(NY,NX)
  plt_bgcr%Canopy_NEE_col                                = Canopy_NEE_col(NY,NX)
  plt_distb%FERT(1:20)                                   = FERT(1:20,I1,NY,NX)
  plt_ew%HeatCanopy2Dist_col                             = HeatCanopy2Dist_col(NY,NX)
  plt_ew%HeatCanopy2Dist_col                             = HeatCanopy2Dist_col(NY,NX)
  DO  L=1,NumOfCanopyLayers
    plt_morph%CanopyStemAareZ_col(L) = CanopyStemAareZ_col(L,NY,NX)
    plt_biom%tCanLeafC_cl(L)         = tCanLeafC_cl(L,NY,NX)
    plt_morph%CanopyLeafAareZ_col(L) = CanopyLeafAareZ_col(L,NY,NX)
  ENDDO

  DO L=0,NL(NY,NX)
    DO K=1,jcplx
      DO NE=1,NumPlantChemElms
        plt_bgcr%REcoDOMProd_vr(NE,K,L)=REcoDOMProd_vr(NE,K,L,NY,NX)
      ENDDO
    ENDDO
  ENDDO

  DO L=0,NL(NY,NX)
    plt_bgcr%REcoH2PO4DmndBand_vr(L) = REcoH2PO4DmndBand_vr(L,NY,NX)
    plt_bgcr%REcoH1PO4DmndBand_vr(L) = REcoH1PO4DmndBand_vr(L,NY,NX)
    plt_bgcr%REcoNO3DmndBand_vr(L)   = REcoNO3DmndBand_vr(L,NY,NX)
    plt_bgcr%REcoNH4DmndBand_vr(L)   = REcoNH4DmndBand_vr(L,NY,NX)
    plt_bgcr%REcoH1PO4DmndSoil_vr(L) = REcoH1PO4DmndSoil_vr(L,NY,NX)
    plt_bgcr%REcoH2PO4DmndSoil_vr(L) = REcoH2PO4DmndSoil_vr(L,NY,NX)
    plt_bgcr%REcoNO3DmndSoil_vr(L)   = REcoNO3DmndSoil_vr(L,NY,NX)
    plt_bgcr%REcoNH4DmndSoil_vr(L)   = REcoNH4DmndSoil_vr(L,NY,NX)
    plt_bgcr%REcoO2DmndResp_vr(L)    = REcoO2DmndResp_vr(L,NY,NX)
    plt_ew%THeatLossRoot2Soil_vr(L)     = THeatLossRoot2Soil_vr(L,NY,NX)
    plt_ew%TPlantRootH2OLoss_vr(L) = TPlantRootH2OLoss_vr(L,NY,NX)
    DO  K=1,micpar%NumOfPlantLitrCmplxs
      DO  M=1,jsken
        DO NE=1,NumPlantChemElms        
          plt_bgcr%LitrfalStrutElms_vr(NE,M,K,L)=LitrfalStrutElms_vr(NE,M,K,L,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO L=1,NL(NY,NX)
    plt_morph%totRootLenDens_vr(L)                   = totRootLenDens_vr(L,NY,NX)
    plt_rbgc%trcg_root_vr(idg_beg:idg_NH3,L)         = trcg_root_vr(idg_beg:idg_NH3,L,NY,NX)
    plt_rbgc%trcg_air2root_flx_vr(idg_beg:idg_NH3,L) = trcg_air2root_flx_vr(idg_beg:idg_NH3,L,NY,NX)
    plt_bgcr%tRootCO2Emis2Root_vr(L)                 = tRootCO2Emis2Root_vr(L,NY,NX)
    plt_bgcr%tRO2MicrbUptk_vr(L)                     = tRO2MicrbUptk_vr(L,NY,NX)
    DO  K=1,jcplx
      plt_bgcr%tRootMycoExud2Soil_vr(1:NumPlantChemElms,K,L)=tRootMycoExud2Soil_vr(1:NumPlantChemElms,K,L,NY,NX)
    ENDDO
  ENDDO


  DO NZ=1,NP0(NY,NX)
    plt_biom%RootElms_pft(1:NumPlantChemElms,NZ)                   = RootElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%ShootElms_pft(1:NumPlantChemElms,NZ)                  = ShootElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%LeafStrutElms_pft(1:NumPlantChemElms,NZ)              = LeafStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%EarStrutElms_pft(1:NumPlantChemElms,NZ)               = EarStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%RootStrutElms_pft(1:NumPlantChemElms,NZ)              = RootStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%HuskStrutElms_pft(1:NumPlantChemElms,NZ)              = HuskStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%StalkRsrvElms_pft(1:NumPlantChemElms,NZ)              = StalkRsrvElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%StalkStrutElms_pft(1:NumPlantChemElms,NZ)             = StalkStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%PetoleStrutElms_pft(1:NumPlantChemElms,NZ)            = PetoleStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%GrainStrutElms_pft(1:NumPlantChemElms,NZ)             = GrainStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_site%ElmBalanceCum_pft(1:NumPlantChemElms,NZ)              = ElmBalanceCum_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ)         = CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%CanopyNonstElms_pft(1:NumPlantChemElms,NZ)            = CanopyNonstElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%CanopyNodulElms_pft(1:NumPlantChemElms,NZ)            = CanopyNodulElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%CanopyNodulNonstElms_pft(1:NumPlantChemElms,NZ)       = CanopyNodulNonstElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_bgcr%LitrfalStrutElms_pft(1:NumPlantChemElms,NZ)           = LitrfalStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_distb%EcoHavstElmnt_CumYr_pft(1:NumPlantChemElms,NZ)       = EcoHavstElmnt_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_pheno%NetCumElmntFlx2Plant_pft(1:NumPlantChemElms,NZ)      = NetCumElmntFlx2Plant_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_bgcr%SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ) = SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_bgcr%LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ)     = LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_rbgc%PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ)         = PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_rbgc%RootUptk_N_CumYr_pft(NZ)                              = RootUptk_N_CumYr_pft(NZ,NY,NX)
    plt_rbgc%RootUptk_P_CumYr_pft(NZ)                              = RootUptk_P_CumYr_pft(NZ,NY,NX)
    plt_distb%EcoHavstElmntCum_pft(1:NumPlantChemElms,NZ)          = EcoHavstElmntCum_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_rbgc%RootMycoExudElms_pft(1:NumPlantChemElms,NZ)           = RootMycoExudElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%SeasonalNonstElms_pft(1:NumPlantChemElms,NZ)          = SeasonalNonstElms_pft(1:NumPlantChemElms,NZ,NY,NX)

    plt_biom%ShootStrutElms_pft(1:NumPlantChemElms,NZ)             = ShootStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%StandDeadStrutElms_pft(1:NumPlantChemElms,NZ)         = StandDeadStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_biom%NodulStrutElms_pft(1:NumPlantChemElms,NZ)             = NodulStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)

    plt_ew%TKCanopy_pft(NZ)            = TKCanopy_pft(NZ,NY,NX)
    plt_photo%LeafO2Solubility_pft(NZ) = LeafO2Solubility_pft(NZ,NY,NX)
    plt_ew%PSICanPDailyMin(NZ)         = PSICanPDailyMin(NZ,NY,NX)
    plt_ew%HeatStorCanopy_pft(NZ)      = HeatStorCanopy_pft(NZ,NY,NX)
    plt_photo%CanopyGasCO2_pft(NZ)     = CanopyGasCO2_pft(NZ,NY,NX)
    plt_photo%O2L(NZ)                  = O2L(NZ,NY,NX)

    plt_photo%aquCO2Intraleaf_pft(NZ)      = aquCO2Intraleaf_pft(NZ,NY,NX)
    plt_distb%FracBiomHarvsted(1:2,1:4,NZ) = FracBiomHarvsted(1:2,1:4,NZ,I,NY,NX)

    plt_pheno%iPlantState_pft(NZ)    = iPlantState_pft(NZ,NY,NX)
    plt_distb%iYearPlanting_pft(NZ)  = iYearPlanting_pft(NZ,NY,NX)
    plt_morph%NumCogrowthNode_pft(NZ) = NumCogrowthNode_pft(NZ,NY,NX)

    plt_pheno%iPlantShootState_pft(NZ) = iPlantShootState_pft(NZ,NY,NX)
    plt_pheno%iPlantRootState_pft(NZ)  = iPlantRootState_pft(NZ,NY,NX)
    plt_allom%CNRTS_pft(NZ)            = CNRTS_pft(NZ,NY,NX)
    plt_allom%CPRTS_pft(NZ)            = CPRTS_pft(NZ,NY,NX)

    plt_rbgc%ZERO4Uptk_pft(NZ)              = ZERO4Uptk_pft(NZ,NY,NX)
    plt_pheno%SeedTempSens_pft(NZ)          = SeedTempSens_pft(NZ,NY,NX)
    plt_rad%FracPARads2Canopy_pft(NZ)       = FracPARads2Canopy_pft(NZ,NY,NX)
    plt_bgcr%NH3Dep2Can_pft(NZ)             = NH3Dep2Can_pft(NZ,NY,NX)
    plt_ew%DeltaTKC_pft(NZ)                 = DeltaTKC_pft(NZ,NY,NX)
    plt_pheno%iPlantThermoAdaptZone_pft(NZ) = iPlantThermoAdaptZone_pft(NZ,NY,NX)
    plt_pheno%HighTempLimitSeed_pft(NZ)     = HighTempLimitSeed_pft(NZ,NY,NX)
    plt_ew%Transpiration_pft(NZ)            = Transpiration_pft(NZ,NY,NX)

    plt_biom%CanopyStalkC_pft(NZ) = CanopyStalkC_pft(NZ,NY,NX)

    plt_photo%ChillHours_pft(NZ) = ChillHours_pft(NZ,NY,NX)
    plt_ew%PSICanopyOsmo_pft(NZ) = PSICanopyOsmo_pft(NZ,NY,NX)

    plt_ew%TdegCCanopy_pft(NZ)           = TdegCCanopy_pft(NZ,NY,NX)
    plt_allom%RootFracRemobilizableBiom(NZ) = RootFracRemobilizableBiom(NZ,NY,NX)
    plt_photo%H2OCuticleResist_pft(NZ)      = H2OCuticleResist_pft(NZ,NY,NX)
    plt_biom%RootBiomCPerPlant_pft(NZ)      = RootBiomCPerPlant_pft(NZ,NY,NX)
    plt_morph%ClumpFactor_pft(NZ)           = ClumpFactor_pft(NZ,NY,NX)

    plt_site%PPI_pft(NZ)               = PPI_pft(NZ,NY,NX)
    plt_site%PPX_pft(NZ)               = PPX_pft(NZ,NY,NX)
    plt_distb%O2ByFire_CumYr_pft(NZ)   = O2ByFire_CumYr_pft(NZ,NY,NX)
    plt_ew%ENGYX_pft(NZ)               = ENGYX_pft(NZ,NY,NX)
    plt_bgcr%CanopyRespC_CumYr_pft(NZ) = CanopyRespC_CumYr_pft(NZ,NY,NX)
    plt_ew%ETCanopy_CumYr_pft(NZ)      = ETCanopy_CumYr_pft(NZ,NY,NX)

    plt_morph%BranchNumber_pft(NZ)    = BranchNumber_pft(NZ,NY,NX)
    plt_morph%NGTopRootLayer_pft(NZ)  = NGTopRootLayer_pft(NZ,NY,NX)
    plt_pheno%doInitPlant_pft(NZ)     = doInitPlant_pft(NZ,NY,NX)
    plt_morph%NIXBotRootLayer_pft(NZ) = NIXBotRootLayer_pft(NZ,NY,NX)
    plt_morph%NumRootAxes_pft(NZ)     = NumRootAxes_pft(NZ,NY,NX)
    plt_morph%MainBranchNum_pft(NZ)   = MainBranchNum_pft(NZ,NY,NX)
    plt_morph%NumOfBranches_pft(NZ)   = NumOfBranches_pft(NZ,NY,NX)

    plt_pheno%IsPlantActive_pft(NZ)     = IsPlantActive_pft(NZ,NY,NX)
    plt_distb%iDayPlanting_pft(NZ)      = iDayPlanting_pft(NZ,NY,NX)
    plt_distb%iDayPlantHarvest_pft(NZ)  = iDayPlantHarvest_pft(NZ,NY,NX)
    plt_distb%iYearPlantHarvest_pft(NZ) = iYearPlantHarvest_pft(NZ,NY,NX)

    plt_distb%FracCanopyHeightCut_pft(NZ) = FracCanopyHeightCut_pft(NZ,I,NY,NX)
    plt_distb%iHarvstType_pft(NZ)         = iHarvstType_pft(NZ,I,NY,NX)
    plt_distb%jHarvst_pft(NZ)             = jHarvst_pft(NZ,I,NY,NX)
    plt_distb%THIN_pft(NZ)                = THIN_pft(NZ,I,NY,NX)
    plt_morph%CanopyStemArea_pft(NZ)      = CanopyStemArea_pft(NZ,NY,NX)
    plt_morph%CanopyLeafArea_pft(NZ)      = CanopyLeafArea_pft(NZ,NY,NX)
    plt_biom%CanopyMassC_pft(NZ)          = CanopyMassC_pft(NZ,NY,NX)           
    plt_photo%O2I(NZ)                      = O2I(NZ,NY,NX)
    plt_photo%LeafIntracellularCO2_pft(NZ) = LeafIntracellularCO2_pft(NZ,NY,NX)

    plt_biom%NoduleNonstructCconc_pft(NZ) = NoduleNonstructCconc_pft(NZ,NY,NX)
    plt_bgcr%CO2NetFix_pft(NZ)            = CO2NetFix_pft(NZ,NY,NX)

    plt_allom%rCNNonstRemob_pft(NZ) = rCNNonstRemob_pft(NZ,NY,NX)
    plt_allom%rCPNonstRemob_pft(NZ) = rCPNonstRemob_pft(NZ,NY,NX)

    plt_photo%DiffCO2Atmos2Intracel_pft(NZ) = DiffCO2Atmos2Intracel_pft(NZ,NY,NX)
    plt_photo%AirConc_pft(NZ)               = AirConc_pft(NZ,NY,NX)
    plt_allom%FracGroth2Node_pft(NZ)        = FracGroth2Node_pft(NZ,NY,NX)

    plt_morph%HypoctoHeight_pft(NZ)                          = HypoctoHeight_pft(NZ,NY,NX)
    plt_rbgc%PlantRootSoilElmNetX_pft(1:NumPlantChemElms,NZ) = PlantRootSoilElmNetX_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_morph%CanopyHeight4WatUptake_pft(NZ)                       = CanopyHeight4WatUptake_pft(NZ,NY,NX)
    plt_morph%MaxSoiL4Root_pft(NZ)                           = MaxSoiL4Root_pft(NZ,NY,NX)
    plt_photo%CanopyBndlResist_pft(NZ)                       = CanopyBndlResist_pft(NZ,NY,NX)
    plt_photo%CanPStomaResistH2O_pft(NZ)                     = CanPStomaResistH2O_pft(NZ,NY,NX)
    plt_ew%TKC_pft(NZ)                                           = TKC_pft(NZ,NY,NX)
    plt_ew%HeatXAir2PCan_pft(NZ)                             = HeatXAir2PCan_pft(NZ,NY,NX)
    plt_rad%RadNet2Canopy_pft(NZ)                            = RadNet2Canopy_pft(NZ,NY,NX)
    plt_rad%LWRadCanopy_pft(NZ)                              = LWRadCanopy_pft(NZ,NY,NX)
    plt_ew%EvapTransLHeat_pft(NZ)                             = EvapTransLHeat_pft(NZ,NY,NX)
    plt_ew%VapXAir2Canopy_pft(NZ)                            = VapXAir2Canopy_pft(NZ,NY,NX)
    plt_photo%MinCanPStomaResistH2O_pft(NZ)                  = MinCanPStomaResistH2O_pft(NZ,NY,NX)
    plt_pheno%TempOffset_pft(NZ)                             = TempOffset_pft(NZ,NY,NX)

    plt_pheno%PlantO2Stress_pft(NZ)                       = PlantO2Stress_pft(NZ,NY,NX)
    plt_site%PlantPopulation_pft(NZ)                      = PlantPopulation_pft(NZ,NY,NX)
    plt_ew%PSICanopy_pft(NZ)                              = PSICanopy_pft(NZ,NY,NX)
    plt_ew%PSICanopyTurg_pft(NZ)                          = PSICanopyTurg_pft(NZ,NY,NX)
    plt_photo%CO2CuticleResist_pft(NZ)                    = CO2CuticleResist_pft(NZ,NY,NX)
    plt_bgcr%RootGasLossDisturb_pft(idg_beg:idg_NH3,NZ) = RootGasLossDisturb_pft(idg_beg:idg_NH3,NZ,NY,NX)
    plt_ew%ReistanceCanopy_pft(NZ)                        = ReistanceCanopy_pft(NZ,NY,NX)
    plt_photo%CO2Solubility_pft(NZ)                       = CO2Solubility_pft(NZ,NY,NX)
    plt_morph%SeedDepth_pft(NZ)                           = SeedDepth_pft(NZ,NY,NX)
    plt_morph%PlantinDepz_pft(NZ)                         = PlantinDepz_pft(NZ,NY,NX)
    plt_morph%SeedMeanLen_pft(NZ)                         = SeedMeanLen_pft(NZ,NY,NX)
    plt_morph%SeedVolumeMean_pft(NZ)                      = SeedVolumeMean_pft(NZ,NY,NX)
    plt_morph%SeedAreaMean_pft(NZ)                        = SeedAreaMean_pft(NZ,NY,NX)
    plt_pheno%TC4LeafOut_pft(NZ)                          = TC4LeafOut_pft(NZ,NY,NX)
    plt_pheno%TCGroth_pft(NZ)                             = TCGroth_pft(NZ,NY,NX)
    plt_pheno%TC4LeafOff_pft(NZ)                          = TC4LeafOff_pft(NZ,NY,NX)
    plt_pheno%TKGroth_pft(NZ)                             = TKGroth_pft(NZ,NY,NX)
    plt_pheno%fTCanopyGroth_pft(NZ)                       = fTCanopyGroth_pft(NZ,NY,NX)

    plt_photo%Km4LeafaqCO2_pft(NZ)                         = Km4LeafaqCO2_pft(NZ,NY,NX)
    plt_photo%Km4RubiscoCarboxy_pft(NZ)                    = Km4RubiscoCarboxy_pft(NZ,NY,NX)
    plt_bgcr%NH3Emis_CumYr_pft(NZ)                         = NH3Emis_CumYr_pft(NZ,NY,NX)
    plt_bgcr%NodulInfectElms_pft(1:NumPlantChemElms,NZ)    = NodulInfectElms_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_bgcr%NodulInfectElmsCum_pft(1:NumPlantChemElms,NZ) = NodulInfectElmsCum_pft(1:NumPlantChemElms,NZ,NY,NX)
    plt_bgcr%PlantN2Fix_CumYr_pft(NZ)                      = PlantN2Fix_CumYr_pft(NZ,NY,NX)
    plt_rbgc%RootN2Fix_pft(NZ)                             = RootN2Fix_pft(NZ,NY,NX)
    plt_rbgc%RootNO3Uptake_pft(NZ)                         = RootNO3Uptake_pft(NZ,NY,NX)
    plt_rbgc%RootNH4Uptake_pft(NZ)                         = RootNH4Uptake_pft(NZ,NY,NX)
    plt_rbgc%RootHPO4Uptake_pft(NZ)                        = RootHPO4Uptake_pft(NZ,NY,NX)
    plt_rbgc%RootH2PO4Uptake_pft(NZ)                       = RootH2PO4Uptake_pft(NZ,NY,NX)
    plt_ew%WatHeldOnCanopy_pft(NZ)                         = WatHeldOnCanopy_pft(NZ,NY,NX)
    plt_ew%VHeatCapCanopy_pft(NZ)                            = VHeatCapCanopy_pft(NZ,NY,NX)
    plt_distb%CH4ByFire_CumYr_pft(NZ)                      = CH4ByFire_CumYr_pft(NZ,NY,NX)
    plt_distb%CO2ByFire_CumYr_pft(NZ)                      = CO2ByFire_CumYr_pft(NZ,NY,NX)
    plt_distb%N2ObyFire_CumYr_pft(NZ)                      = N2ObyFire_CumYr_pft(NZ,NY,NX)
    plt_distb%NH3byFire_CumYr_pft(NZ)                      = NH3byFire_CumYr_pft(NZ,NY,NX)
    plt_distb%PO4byFire_CumYr_pft(NZ)                      = PO4byFire_CumYr_pft(NZ,NY,NX)
    plt_ew%CanopyBiomWater_pft(NZ)                             = CanopyBiomWater_pft(NZ,NY,NX)
    plt_pheno%HoursTooLowPsiCan_pft(NZ)                    = HoursTooLowPsiCan_pft(NZ,NY,NX)
    plt_biom%SeedCPlanted_pft(NZ)                          = SeedCPlanted_pft(NZ,NY,NX)
    plt_biom%CanopyLeafShethC_pft(NZ)                      = CanopyLeafShethC_pft(NZ,NY,NX)

    plt_biom%ZERO4LeafVar_pft(NZ)   = ZERO4LeafVar_pft(NZ,NY,NX)
    plt_biom%ZERO4Groth_pft(NZ)     = ZERO4Groth_pft(NZ,NY,NX)
    plt_morph%CanopyHeight_pft(NZ)  = CanopyHeight_pft(NZ,NY,NX)
    plt_bgcr%NetPrimProduct_pft(NZ) = NetPrimProduct_pft(NZ,NY,NX)


    DO L=1,NL(NY,NX)
      DO K=1,jcplx
        DO N=1,MY(NZ,NY,NX)
          DO NE=1,NumPlantChemElms
            plt_rbgc%RootMycoExudEUptk_pvr(NE,N,K,L,NZ)=RootMycoExudEUptk_pvr(NE,N,K,L,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    
    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO NE=1,NumPlantChemElms
        plt_biom%CanopyNonstElms_brch(NE,NB,NZ)      = CanopyNonstElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%CanopyNodulNonstElms_brch(NE,NB,NZ) = CanopyNodulNonstElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%ShootStrutElms_brch(NE,NB,NZ)       = ShootStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%LeafPetoNonstElmConc_brch(NE,NB,NZ) = LeafPetoNonstElmConc_brch(NE,NB,NZ,NY,NX)
        plt_biom%PetoleStrutElms_brch(NE,NB,NZ)      = PetoleStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%StalkStrutElms_brch(NE,NB,NZ)       = StalkStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%LeafStrutElms_brch(NE,NB,NZ)        = LeafStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%StalkRsrvElms_brch(NE,NB,NZ)        = StalkRsrvElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%HuskStrutElms_brch(NE,NB,NZ)        = HuskStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%GrainStrutElms_brch(NE,NB,NZ)       = GrainStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%EarStrutElms_brch(NE,NB,NZ)         = EarStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%CanopyNodulStrutElms_brch(NE,NB,NZ) = CanopyNodulStrutElms_brch(NE,NB,NZ,NY,NX)
        plt_biom%PetioleChemElmRemob_brch(NE,NB,NZ)  = PetioleChemElmRemob_brch(NE,NB,NZ,NY,NX)
      ENDDO
      plt_biom%ShootC4NonstC_brch(NB,NZ)=ShootC4NonstC_brch(NB,NZ,NY,NX)              
    ENDDO

    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      plt_photo%RubiscoActivity_brch(NB,NZ)  = RubiscoActivity_brch(NB,NZ,NY,NX)
      plt_photo%C4PhotosynDowreg_brch(NB,NZ) = C4PhotosynDowreg_brch(NB,NZ,NY,NX)
      plt_pheno%Hours2LeafOut_brch(NB,NZ)    = Hours2LeafOut_brch(NB,NZ,NY,NX)
      plt_morph%LeafAreaDying_brch(NB,NZ)    = LeafAreaDying_brch(NB,NZ,NY,NX)
      plt_morph%LeafAreaLive_brch(NB,NZ)     = LeafAreaLive_brch(NB,NZ,NY,NX)

      plt_pheno%dReproNodeNumNormByMatG_brch(NB,NZ)    = dReproNodeNumNormByMatG_brch(NB,NZ,NY,NX)
      plt_pheno%HourFailGrainFill_brch(NB,NZ)          = HourFailGrainFill_brch(NB,NZ,NY,NX)
      plt_pheno%HoursDoingRemob_brch(NB,NZ)            = HoursDoingRemob_brch(NB,NZ,NY,NX)
      plt_pheno%MatureGroup_brch(NB,NZ)                = MatureGroup_brch(NB,NZ,NY,NX)
      plt_pheno%NodeNumNormByMatgrp_brch(NB,NZ)        = NodeNumNormByMatgrp_brch(NB,NZ,NY,NX)
      plt_pheno%ReprodNodeNumNormByMatrgrp_brch(NB,NZ) = ReprodNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX)
      plt_morph%PotentialSeedSites_brch(NB,NZ)         = PotentialSeedSites_brch(NB,NZ,NY,NX)
      plt_morph%SeedNumSet_brch(NB,NZ)                 = SeedNumSet_brch(NB,NZ,NY,NX)
      plt_allom%GrainSeedBiomCMean_brch(NB,NZ)         = GrainSeedBiomCMean_brch(NB,NZ,NY,NX)
      plt_morph%CanPBranchHeight(NB,NZ)                = CanPBranchHeight(NB,NZ,NY,NX)
      plt_pheno%iPlantBranchState_brch(NB,NZ)          = iPlantBranchState_brch(NB,NZ,NY,NX)
      plt_pheno%doRemobilization_brch(NB,NZ)           = doRemobilization_brch(NB,NZ,NY,NX)
      plt_pheno%doPlantLeaveOff_brch(NB,NZ)            = doPlantLeaveOff_brch(NB,NZ,NY,NX)
      plt_pheno%doPlantLeafOut_brch(NB,NZ)             = doPlantLeafOut_brch(NB,NZ,NY,NX)
      plt_pheno%doInitLeafOut_brch(NB,NZ)              = doInitLeafOut_brch(NB,NZ,NY,NX)
      plt_pheno%doSenescence_brch(NB,NZ)               = doSenescence_brch(NB,NZ,NY,NX)
      plt_pheno%Prep4Literfall_brch(NB,NZ)             = Prep4Literfall_brch(NB,NZ,NY,NX)
      plt_pheno%Hours4LiterfalAftMature_brch(NB,NZ)    = Hours4LiterfalAftMature_brch(NB,NZ,NY,NX)
      plt_pheno%KHiestGroLeafNode_brch(NB,NZ)          = KHiestGroLeafNode_brch(NB,NZ,NY,NX)
      plt_pheno%KLowestGroLeafNode_brch(NB,NZ)         = KLowestGroLeafNode_brch(NB,NZ,NY,NX)
      plt_morph%BranchNumber_brch(NB,NZ)               = BranchNumber_brch(NB,NZ,NY,NX)
      plt_morph%ShootNodeNum_brch(NB,NZ)               = ShootNodeNum_brch(NB,NZ,NY,NX)
      plt_morph%NodeNum2InitFloral_brch(NB,NZ)                        = NodeNum2InitFloral_brch(NB,NZ,NY,NX)
      plt_morph%NodeNumberAtAnthesis_brch(NB,NZ)                      = NodeNumberAtAnthesis_brch(NB,NZ,NY,NX)
      plt_pheno%LeafElmntRemobFlx_brch(1:NumPlantChemElms,NB,NZ)      = LeafElmntRemobFlx_brch(1:NumPlantChemElms,NB,NZ,NY,NX)
      plt_pheno%PetioleChemElmRemobFlx_brch(1:NumPlantChemElms,NB,NZ) = PetioleChemElmRemobFlx_brch(1:NumPlantChemElms,NB,NZ,NY,NX)
      plt_pheno%TotalNodeNumNormByMatgrp_brch(NB,NZ)                  = TotalNodeNumNormByMatgrp_brch(NB,NZ,NY,NX)
      plt_pheno%TotReproNodeNumNormByMatrgrp_brch(NB,NZ)              = TotReproNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX)
      plt_pheno%LeafNumberAtFloralInit_brch(NB,NZ)                    = LeafNumberAtFloralInit_brch(NB,NZ,NY,NX)
      plt_morph%NumOfLeaves_brch(NB,NZ)                               = NumOfLeaves_brch(NB,NZ,NY,NX)
      plt_pheno%Hours4LenthenPhotoPeriod_brch(NB,NZ)                  = Hours4LenthenPhotoPeriod_brch(NB,NZ,NY,NX)
      plt_pheno%Hours4ShortenPhotoPeriod_brch(NB,NZ)                  = Hours4ShortenPhotoPeriod_brch(NB,NZ,NY,NX)
      plt_pheno%Hours4Leafout_brch(NB,NZ)                             = Hours4Leafout_brch(NB,NZ,NY,NX)
      plt_pheno%Hours4LeafOff_brch(NB,NZ)                             = Hours4LeafOff_brch(NB,NZ,NY,NX)
      plt_biom%LeafPetolBiomassC_brch(NB,NZ)                          = LeafPetolBiomassC_brch(NB,NZ,NY,NX)
      plt_biom%LeafChemElmRemob_brch(1:NumPlantChemElms,NB,NZ)        = LeafChemElmRemob_brch(1:NumPlantChemElms,NB,NZ,NY,NX)
      plt_biom%SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)     = SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ,NY,NX)
      plt_biom%StalkLiveBiomassC_brch(NB,NZ)                              = StalkLiveBiomassC_brch(NB,NZ,NY,NX)
      DO M=1,pltpar%NumGrowthStages
        plt_pheno%iPlantCalendar_brch(M,NB,NZ)=iPlantCalendar_brch(M,NB,NZ,NY,NX)
      ENDDO
      DO K=1,MaxNodesPerBranch
        plt_photo%CPOOL3_node(K,NB,NZ)                = CPOOL3_node(K,NB,NZ,NY,NX)
        plt_photo%CPOOL4_node(K,NB,NZ)                = CPOOL4_node(K,NB,NZ,NY,NX)
        plt_photo%CMassCO2BundleSheath_node(K,NB,NZ)  = CMassCO2BundleSheath_node(K,NB,NZ,NY,NX)
        plt_photo%CO2CompenPoint_node(K,NB,NZ)        = CO2CompenPoint_node(K,NB,NZ,NY,NX)
        plt_photo%RubiscoCarboxyEff_node(K,NB,NZ)     = RubiscoCarboxyEff_node(K,NB,NZ,NY,NX)
        plt_photo%CMassHCO3BundleSheath_node(K,NB,NZ) = CMassHCO3BundleSheath_node(K,NB,NZ,NY,NX)
      ENDDO
      DO K=0,MaxNodesPerBranch
        plt_morph%LeafNodeArea_brch(K,NB,NZ)                         = LeafNodeArea_brch(K,NB,NZ,NY,NX)
        plt_morph%InternodeHeightDead_brch(K,NB,NZ)                 = InternodeHeightDead_brch(K,NB,NZ,NY,NX)
        plt_morph%PetoleLensNode_brch(K,NB,NZ)                       = PetoleLensNode_brch(K,NB,NZ,NY,NX)
        plt_morph%LiveInterNodeHight_brch(K,NB,NZ)                   = LiveInterNodeHight_brch(K,NB,NZ,NY,NX)
        plt_biom%InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ) = InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)
        plt_biom%LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)      = LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)
        plt_biom%LeafProteinCNode_brch(K,NB,NZ)                      = LeafProteinCNode_brch(K,NB,NZ,NY,NX)
        plt_biom%PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)   = PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)
        plt_biom%PetoleProteinCNode_brch(K,NB,NZ)                    = PetoleProteinCNode_brch(K,NB,NZ,NY,NX)
      ENDDO

      DO K=0,MaxNodesPerBranch
        DO  L=1,NumOfCanopyLayers                    
          plt_morph%CanopyLeafArea_lpft(L,K,NB,NZ)                        = CanopyLeafArea_lpft(L,K,NB,NZ,NY,NX)
          plt_biom%LeafElmsByLayerNode_brch(1:NumPlantChemElms,L,K,NB,NZ) = LeafElmsByLayerNode_brch(1:NumPlantChemElms,L,K,NB,NZ,NY,NX)
        ENDDO
      ENDDO
      DO  L=1,NumOfCanopyLayers
        plt_morph%CanopyStalkArea_lbrch(L,NB,NZ)=CanopyStalkArea_lbrch(L,NB,NZ,NY,NX)
      ENDDO
    enddo

    DO L=1,NL(NY,NX)
      DO NE=1,NumPlantChemElms
        plt_biom%RootNodulStrutElms_rpvr(NE,L,NZ) =RootNodulStrutElms_rpvr(NE,L,NZ,NY,NX)
      ENDDO
    ENDDO
    DO L=1,NL(NY,NX)
      DO N=1,MY(NZ,NY,NX)
        plt_biom%RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ) = RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX)
        plt_rbgc%trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)         = trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX)
        plt_rbgc%trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)         = trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ,NY,NX)

        plt_biom%RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ) = RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX)
        plt_biom%RootProteinConc_rpvr(N,L,NZ)                         = RootProteinConc_rpvr(N,L,NZ,NY,NX)

        plt_ew%PSIRoot_pvr(N,L,NZ)                                = PSIRoot_pvr(N,L,NZ,NY,NX)
        plt_ew%PSIRootOSMO_vr(N,L,NZ)                             = PSIRootOSMO_vr(N,L,NZ,NY,NX)
        plt_ew%PSIRootTurg_vr(N,L,NZ)                             = PSIRootTurg_vr(N,L,NZ,NY,NX)
        plt_rbgc%RootRespPotent_pvr(N,L,NZ)                       = RootRespPotent_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootCO2EmisPot_pvr(N,L,NZ)       = RootCO2EmisPot_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootCO2Autor_pvr(N,L,NZ)         = RootCO2Autor_pvr(N,L,NZ,NY,NX)
        plt_morph%Root1stXNumL_pvr(N,L,NZ)        = Root1stXNumL_pvr(N,L,NZ,NY,NX)
        plt_morph%Root2ndXNum_pvr(N,L,NZ)         = Root2ndXNum_pvr(N,L,NZ,NY,NX)
        plt_morph%RootLenPerPlant_pvr(N,L,NZ)     = RootLenPerPlant_pvr(N,L,NZ,NY,NX)
        plt_morph%RootLenDensPerPlant_pvr(N,L,NZ) = RootLenDensPerPlant_pvr(N,L,NZ,NY,NX)
        plt_morph%RootPoreVol_pvr(N,L,NZ)         = RootPoreVol_pvr(N,L,NZ,NY,NX)
        plt_morph%RootVH2O_pvr(N,L,NZ)            = RootVH2O_pvr(N,L,NZ,NY,NX)
        plt_morph%Root1stRadius_pvr(N,L,NZ)       = Root1stRadius_pvr(N,L,NZ,NY,NX)
        plt_morph%Root2ndRadius_pvr(N,L,NZ)       = Root2ndRadius_pvr(N,L,NZ,NY,NX)
        plt_morph%RootAreaPerPlant_pvr(N,L,NZ)    = RootAreaPerPlant_pvr(N,L,NZ,NY,NX)
        plt_morph%Root2ndAveLen_pvr(N,L,NZ)       = Root2ndAveLen_pvr(N,L,NZ,NY,NX)

        plt_rbgc%RootO2Dmnd4Resp_pvr(N,L,NZ)     = RootO2Dmnd4Resp_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootNH4DmndSoil_pvr(N,L,NZ)     = RootNH4DmndSoil_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootNH4DmndBand_pvr(N,L,NZ)     = RootNH4DmndBand_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootNO3DmndSoil_pvr(N,L,NZ)     = RootNO3DmndSoil_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootNO3DmndBand_pvr(N,L,NZ)     = RootNO3DmndBand_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootH2PO4DmndSoil_pvr(N,L,NZ)   = RootH2PO4DmndSoil_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootH2PO4DmndBand_pvr(N,L,NZ)   = RootH2PO4DmndBand_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootH1PO4DmndSoil_pvr(N,L,NZ)   = RootH1PO4DmndSoil_pvr(N,L,NZ,NY,NX)
        plt_rbgc%RootH1PO4DmndBand_pvr(N,L,NZ)   = RootH1PO4DmndBand_pvr(N,L,NZ,NY,NX)
        plt_ew%AllPlantRootH2OLoss_vr(N,L,NZ)  = AllPlantRootH2OLoss_vr(N,L,NZ,NY,NX)
        plt_rbgc%RAutoRootO2Limter_rpvr(N,L,NZ)   = RAutoRootO2Limter_rpvr(N,L,NZ,NY,NX)
        plt_biom%RootMycoActiveBiomC_pvr(N,L,NZ) = RootMycoActiveBiomC_pvr(N,L,NZ,NY,NX)
        plt_biom%PopuRootMycoC_pvr(N,L,NZ)       = PopuRootMycoC_pvr(N,L,NZ,NY,NX)
        plt_biom%RootProteinC_pvr(N,L,NZ)        = RootProteinC_pvr(N,L,NZ,NY,NX)

      enddo
      plt_biom%RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ)=RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ,NY,NX)
    ENDDO
    DO L=1,NumOfCanopyLayers
      plt_morph%CanopyStemAreaZ_pft(L,NZ) = CanopyStemAreaZ_pft(L,NZ,NY,NX)
      plt_morph%CanopyLeafAreaZ_pft(L,NZ) = CanopyLeafAreaZ_pft(L,NZ,NY,NX)
      plt_biom%CanopyLeafCLyr_pft(L,NZ)   = CanopyLeafCLyr_pft(L,NZ,NY,NX)
    ENDDO
    DO N=1,MY(NZ,NY,NX)
      plt_morph%RootVolPerMassC_pft(N,NZ)   = RootVolPerMassC_pft(N,NZ,NY,NX)
      plt_morph%RootPoreTortu4Gas(N,NZ)     = RootPoreTortu4Gas(N,NZ,NY,NX)
      plt_morph%Root2ndXSecArea_pft(N,NZ)   = Root2ndXSecArea_pft(N,NZ,NY,NX)
      plt_morph%Root1stXSecArea_pft(N,NZ)   = Root1stXSecArea_pft(N,NZ,NY,NX)
      plt_morph%Root1stMaxRadius1_pft(N,NZ) = Root1stMaxRadius1_pft(N,NZ,NY,NX)
      plt_morph%Root2ndMaxRadius1_pft(N,NZ) = Root2ndMaxRadius1_pft(N,NZ,NY,NX)
      plt_morph%Root1stSpecLen_pft(N,NZ)    = Root1stSpecLen_pft(N,NZ,NY,NX)
      plt_morph%RootRaidus_rpft(N,NZ)       = RootRaidus_rpft(N,NZ,NY,NX)
      plt_morph%Root2ndSpecLen_pft(N,NZ)    = Root2ndSpecLen_pft(N,NZ,NY,NX)
    ENDDO
    DO NR=1,NumOfCanopyLayers
      plt_morph%NIXBotRootLayer_rpft(NR,NZ)=NIXBotRootLayer_rpft(NR,NZ,NY,NX)
      DO L=1,MaxNumRootLays(NY,NX)
        DO N=1,MY(NZ,NY,NX)
          plt_morph%Root1stLen_rpvr(N,L,NR,NZ)                             = Root1stLen_rpvr(N,L,NR,NZ,NY,NX)
          plt_morph%Root2ndLen_rpvr(N,L,NR,NZ)                              = Root2ndLen_rpvr(N,L,NR,NZ,NY,NX)
          plt_morph%Root2ndXNum_rpvr(N,L,NR,NZ)                            = Root2ndXNum_rpvr(N,L,NR,NZ,NY,NX)
          plt_biom%RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = &
            RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ,NY,NX)
          plt_biom%RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = &
            RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ,NY,NX)
        enddo
      enddo
      DO N=1,MY(NZ,NY,NX)
        plt_morph%Root1stDepz_pft(N,NR,NZ)                       = Root1stDepz_pft(N,NR,NZ,NY,NX)
        plt_biom%RootMyco1stElm_raxs(1:NumPlantChemElms,N,NR,NZ) = RootMyco1stElm_raxs(1:NumPlantChemElms,N,NR,NZ,NY,NX)
      enddo
    enddo
    
    DO M=1,jsken
      DO NE=1,NumPlantChemElms
        plt_biom%StandDeadKCompElms_pft(NE,M,NZ)=StandDeadKCompElms_pft(NE,M,NZ,NY,NX)
      ENDDO
    ENDDO
!!!!  LitrfalStrutElms_pvr in restart file?  
    DO L=0,MaxNumRootLays(NY,NX)
      DO K=1,micpar%NumOfPlantLitrCmplxs
        DO M=1,jsken
          plt_bgcr%LitrfalStrutElms_pvr(1:NumPlantChemElms,M,K,L,NZ)=LitrfalStrutElms_pvr(1:NumPlantChemElms,M,K,L,NZ,NY,NX)
        enddo
      enddo
    ENDDO

    DO M=1,jsken
      DO N=0,pltpar%NumLitterGroups
        DO NE=1,NumPlantChemElms        
          plt_soilchem%ElmAllocmat4Litr(NE,N,M,NZ)=ElmAllocmat4Litr(NE,N,M,NZ,NY,NX)
        enddo
      enddo
    ENDDO
!!!    
  ENDDO

  DO L=1,NL(NY,NX)
    DO M=1,NPH
      plt_rbgc%RO2UptkSoilM_vr(M,L)           = RO2UptkSoilM_vr(M,L,NY,NX)
      plt_soilchem%AirFilledSoilPoreM_vr(M,L) = AirFilledSoilPoreM_vr(M,L,NY,NX)
    ENDDO
  ENDDO
  end subroutine PlantAPISend

!------------------------------------------------------------------------------------------

  subroutine PlantAPICanMSend(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,N,M,NN,NZ,K,NB
!  Integers
  plt_site%KoppenClimZone        = KoppenClimZone_col(NY,NX)
  plt_site%NP                    = NP(NY,NX)
  plt_site%NU                    = NU(NY,NX)
  plt_morph%StemArea_col         = StemArea_col(NY,NX)
  plt_morph%CanopyLeafArea_col   = CanopyLeafArea_col(NY,NX)
  plt_site%ZEROS                 = ZEROS(NY,NX)
  plt_site%ZERO                  = ZERO
  plt_ew%SnowDepth               = SnowDepth_col(NY,NX)
  plt_ew%TairK                   = TairK_col(NY,NX)
  plt_morph%CanopyHeight_col     = CanopyHeight_col(NY,NX)
  plt_site%WindMesureHeight_col  = WindMesureHeight_col(NY,NX)
  plt_ew%ZERO4PlantDisplace_col  = ZERO4PlantDisplace_col(NY,NX)
  plt_site%WindSpeedAtm_col      = WindSpeedAtm_col(NY,NX)
  plt_ew%VLHeatCapSnowMin_col    = VLHeatCapSnowMin_col(NY,NX)
  plt_ew%VLHeatCapSurfSnow_col   = VLHeatCapSnow_snvr(1,NY,NX)
  plt_morph%CanopyHeightZ_col(0) = CanopyHeightZ_col(0,NY,NX)
  DO L=1,NumOfCanopyLayers
    plt_morph%CanopyStemAareZ_col(L) = CanopyStemAareZ_col(L,NY,NX)
    plt_morph%CanopyLeafAareZ_col(L) = CanopyLeafAareZ_col(L,NY,NX)
    plt_morph%CanopyHeightZ_col(L)   = CanopyHeightZ_col(L,NY,NX)
    plt_rad%TAU_DirRadTransm(L)      = TAU_DirRadTransm(L,NY,NX)
  ENDDO
  plt_rad%TAU_DirRadTransm(NumOfCanopyLayers+1)=TAU_DirRadTransm(NumOfCanopyLayers+1,NY,NX)

  DO L=0,NL(NY,NX)
    plt_site%AREA3(L)                 = AREA(3,L,NY,NX)
    plt_soilchem%VLSoilPoreMicP_vr(L) = VLSoilPoreMicP_vr(L,NY,NX)
    plt_soilchem%VLSoilMicP_vr(L)     = VLSoilMicP_vr(L,NY,NX)
    plt_soilchem%VLWatMicP_vr(L)      = VLWatMicP_vr(L,NY,NX)
  ENDDO
  plt_rad%SineSunInclAngle_col     = SineSunInclAngle_col(NY,NX)
  plt_site%SolarNoonHour_col       = SolarNoonHour_col(NY,NX)
  plt_morph%LeafStalkArea_col      = LeafStalkArea_col(NY,NX)
  plt_rad%GroundSurfAzimuth_col    = GroundSurfAzimuth_col(NY,NX)
  plt_rad%CosineGrndSlope_col      = CosineGrndSlope_col(NY,NX)
  plt_rad%SineGrndSlope_col        = SineGrndSlope_col(NY,NX)
  plt_rad%RadSWDiffus_col          = RadSWDiffus_col(NY,NX)
  plt_rad%RadPARDiffus_col         = RadPARDiffus_col(NY,NX)
  plt_rad%RadSWDirect_col          = RadSWDirect_col(NY,NX)
  plt_rad%RadPARDirect_col         = RadPARDirect_col(NY,NX)
  plt_site%SoilSurfRoughnesst0_col = SoilSurfRoughnesst0_col(NY,NX)
  plt_ew%VcumWatSnow_col           = VcumWatSnow_col(NY,NX)
  plt_ew%VcumIceSnow_col           = VcumIceSnow_col(NY,NX)
  plt_ew%VcumDrySnoWE_col          = VcumDrySnoWE_col(NY,NX)
  plt_rad%TotSineSkyAngles_grd     = TotSineSkyAngles_grd
  plt_rad%SoilAlbedo               = SoilAlbedo_col(NY,NX)
  plt_rad%SurfAlbedo_col           = SurfAlbedo_col(NY,NX)
  plt_site%ZEROS2                  = ZEROS2(NY,NX)
  plt_site%POROS1                  = POROS_vr(NU(NY,NX),NY,NX)
  DO NZ=1,NP(NY,NX)
    plt_morph%CanopyLeafArea_pft(NZ)   = CanopyLeafArea_pft(NZ,NY,NX)
    plt_morph%CanopyHeight_pft(NZ)     = CanopyHeight_pft(NZ,NY,NX)
    plt_morph%ClumpFactorNow_pft(NZ)   = ClumpFactorNow_pft(NZ,NY,NX)
    plt_rad%LeafSWabsorpty_pft(NZ)     = LeafSWabsorpty_pft(NZ,NY,NX)
    plt_rad%LeafPARabsorpty_pft(NZ)    = LeafPARabsorpty_pft(NZ,NY,NX)
    plt_rad%RadSWLeafTransmis_pft(NZ)  = RadSWLeafTransmis_pft(NZ,NY,NX)
    plt_rad%RadSWLeafAlbedo_pft(NZ)    = RadSWLeafAlbedo_pft(NZ,NY,NX)
    plt_rad%RadPARLeafTransmis_pft(NZ) = RadPARLeafTransmis_pft(NZ,NY,NX)
    plt_rad%CanopyPARalbedo_pft(NZ)    = CanopyPARalbedo_pft(NZ,NY,NX)
    plt_morph%NumOfBranches_pft(NZ)    = NumOfBranches_pft(NZ,NY,NX)
    plt_morph%ClumpFactor_pft(NZ)      = ClumpFactor_pft(NZ,NY,NX)

    DO NB=1,NumOfBranches_pft(NZ,NY,NX)
      DO K=0,MaxNodesPerBranch
        DO  L=1,NumOfCanopyLayers
          plt_morph%CanopyLeafArea_lpft(L,K,NB,NZ)=CanopyLeafArea_lpft(L,K,NB,NZ,NY,NX)
        ENDDO
      ENDDO
      DO  L=1,NumOfCanopyLayers
        plt_morph%CanopyStalkArea_lbrch(L,NB,NZ)=CanopyStalkArea_lbrch(L,NB,NZ,NY,NX)
      ENDDO
      DO K=1,MaxNodesPerBranch
        DO  L=1,NumOfCanopyLayers
          DO N=1,NumOfLeafZenithSectors
            plt_morph%LeafAreaZsec_brch(N,L,K,NB,NZ)=LeafAreaZsec_brch(N,L,K,NB,NZ,NY,NX)
          ENDDO
        ENDDO
      ENDDO
      DO  L=1,NumOfCanopyLayers
        DO N=1,NumOfLeafZenithSectors
          plt_morph%StemAreaZsec_brch(N,L,NB,NZ)=StemAreaZsec_brch(N,L,NB,NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO N=1,NumOfSkyAzimuthSects
    plt_rad%OMEGAG(N)=OMEGAG(N,NY,NX)
  ENDDO
  DO N=1,NumOfLeafZenithSectors
    plt_rad%CosineLeafAngle(N) = CosineLeafAngle(N)
    plt_rad%SineLeafAngle(N)   = SineLeafAngle(N)
  ENDDO
  DO NN=1,NumOfLeafAzimuthSectors
    DO M=1,NumOfLeafZenithSectors
      DO N=1,NumOfSkyAzimuthSects
        plt_rad%OMEGA(N,M,NN)             = OMEGA(N,M,NN)
        plt_rad%OMEGX(N,M,NN)             = OMEGX(N,M,NN)
        plt_rad%iScatteringDiffus(N,M,NN) = iScatteringDiffus(N,M,NN)
      ENDDO
    ENDDO
  ENDDO

  end subroutine PlantAPICanMSend

!------------------------------------------------------------------------------------------

  subroutine PlantAPICanMRecv(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  integer :: N,M,NN,L,NZ,K,NB

  ZERO4PlantDisplace_col(NY,NX)  = plt_ew%ZERO4PlantDisplace_col
  RoughHeight_col(NY,NX)         = plt_ew%RoughHeight
  AbvCanopyBndlResist_col(NY,NX) = plt_ew%AbvCanopyBndlResist_col
  RIB_col(NY,NX)                     = plt_ew%RIB
  
  CanopyHeight_col(NY,NX)    = plt_morph%CanopyHeight_col
  RadSWDirect_col(NY,NX)     = plt_rad%RadSWDirect_col
  RadSWDiffus_col(NY,NX)     = plt_rad%RadSWDiffus_col
  RadPARDirect_col(NY,NX)    = plt_rad%RadPARDirect_col
  RadPARDiffus_col(NY,NX)    = plt_rad%RadPARDiffus_col
  RadSWGrnd_col(NY,NX)       = plt_rad%RadSWGrnd_col
  FracSWRad2Grnd_col(NY,NX)  = plt_rad%FracSWRad2Grnd_col
  RadSWSolarBeam_col(NY,NX)  = plt_rad%RadSWSolarBeam_col
  RadPARSolarBeam_col(NY,NX) = plt_rad%RadPARSolarBeam_col
  
  DO L=0,NumOfCanopyLayers
    CanopyHeightZ_col(L,NY,NX)=plt_morph%CanopyHeightZ_col(L)
  ENDDO
  DO L=1,NumOfCanopyLayers
    TAU_DirRadTransm(L,NY,NX) = plt_rad%TAU_DirRadTransm(L)
    TAU_RadThru(L,NY,NX)      = plt_rad%TAU_RadThru(L)
  ENDDO
  LeafStalkArea_col(NY,NX)=plt_morph%LeafStalkArea_col
  DO NZ=1,NP(NY,NX)
    LeafStalkArea_pft(NZ,NY,NX)    =plt_morph%LeafStalkArea_pft(NZ)
    RadSWbyCanopy_pft(NZ,NY,NX)    =plt_rad%RadSWbyCanopy_pft(NZ)
    RadPARbyCanopy_pft(NZ,NY,NX)   =plt_rad%RadPARbyCanopy_pft(NZ)
    ClumpFactorNow_pft(NZ,NY,NX)   =plt_morph%ClumpFactorNow_pft(NZ)
    FracPARads2Canopy_pft(NZ,NY,NX)=plt_rad%FracPARads2Canopy_pft(NZ)
    StomatalStress_pft(NZ,NY,NX)   =plt_biom%StomatalStress_pft(NZ)
    Eco_RadSW_col(NY,NX)           =Eco_RadSW_col(NY,NX)+RadSWbyCanopy_pft(NZ,NY,NX)
    DO L=1,NumOfCanopyLayers
      DO M=1,NumOfSkyAzimuthSects
        DO  N=1,NumOfLeafZenithSectors
          RadDifPAR_zsec(N,M,L,NZ,NY,NX)=plt_rad%RadDifPAR_zsec(N,M,L,NZ)
          RadPAR_zsec(N,M,L,NZ,NY,NX)   =plt_rad%RadPAR_zsec(N,M,L,NZ)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  end subroutine PlantAPICanMRecv

end module PlantAPI
