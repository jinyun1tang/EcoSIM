module PlantTraitTableMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use GridConsts,    only: NumLeafZenithSectors
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  integer, target, allocatable :: iPlantPhotosynthesisType_tab(:)
  integer, target, allocatable :: iPlantRootProfile_tab(:)
  integer, target, allocatable :: iPlantPhenolPattern_tab(:)
  integer, target, allocatable :: iPlantDevelopPattern_tab(:)
  integer, target, allocatable :: iPlantNfixType_tab(:)
  integer, target, allocatable :: iPlantPhenolType_tab(:)
  integer, target, allocatable :: iPlantPhotoperiodType_tab(:)
  integer, target, allocatable :: iPlantTurnoverPattern_tab(:)
  integer, target, allocatable :: iPlantGrainType_tab(:)
  integer, target, allocatable :: Myco_tab(:)
  real(r8), target, allocatable :: PlantInitThermoAdaptZone_tab(:)
  real(r8), target, allocatable :: VmaxRubCarboxyRef_tab(:)
  real(r8), target, allocatable :: VmaxRubOxyRef_tab(:)
  real(r8), target, allocatable :: VmaxPEPCarboxyRef_tab(:)
  real(r8), target, allocatable :: XKCO2_tab(:)
  real(r8), target, allocatable :: XKO2_tab(:)
  real(r8), target, allocatable :: Km4PEPCarboxy_tab(:)
  real(r8), target, allocatable :: LeafRuBPConc_tab(:)
  real(r8), target, allocatable :: FracLeafProtAsPEPCarboxyl_tab(:)
  real(r8), target, allocatable :: SpecChloryfilAct_tab(:)
  real(r8), target, allocatable :: LeafC3ChlorofilConc_tab(:)
  real(r8), target, allocatable :: LeafC4ChlorofilConc_tab(:)
  real(r8), target, allocatable :: CanPCi2CaRatio_tab(:)
  real(r8), target, allocatable :: RadSWLeafAlbedo_tab(:)
  real(r8), target, allocatable :: CanopyPARalbedo_tab(:)
  real(r8), target, allocatable :: RadSWLeafTransmis_tab(:)
  real(r8), target, allocatable :: RadPARLeafTransmis_tab(:)
  real(r8), target, allocatable :: RefNodeInitRate_tab(:)
  real(r8), target, allocatable :: RefLeafAppearRate_tab(:)
  real(r8), target, allocatable :: TCChill4Seed_tab(:)
  real(r8), target, allocatable :: VRNLI_tab(:)
  real(r8), target, allocatable :: VRNXI_tab(:)
  real(r8), target, allocatable :: rLen2WidthLeaf_tab(:)
  real(r8), target, allocatable :: NonstCMinConc2InitBranch_tab(:)
  real(r8), target, allocatable :: GROUPX_tab(:)
  real(r8), target, allocatable :: ShootNodeNumAtPlanting_tab(:)
  real(r8), target, allocatable :: CriticPhotoPeriod_tab(:)
  real(r8), target, allocatable :: PhotoPeriodSens_tab(:)
  real(r8), target, allocatable :: SLA1_tab(:)
  real(r8), target, allocatable :: PetoLen2Mass_tab(:)
  real(r8), target, allocatable :: NodeLenPergC_tab(:)
  real(r8), target, allocatable :: LeafAngleClass_tab(:,:)
  real(r8), target, allocatable :: ClumpFactorInit_tab(:)
  real(r8), target, allocatable :: BranchAngle_tab(:)
  real(r8), target, allocatable :: PetioleAngle_tab(:)
  real(r8), target, allocatable :: MaxPotentSeedNumber_tab(:)
  real(r8), target, allocatable :: MaxSeedNumPerSite_tab(:)
  real(r8), target, allocatable :: SeedCMassMax_tab(:)
  real(r8), target, allocatable :: SeedCMass_tab(:)
  real(r8), target, allocatable :: GrainFillRate25C_tab(:)
  real(r8), target, allocatable :: StandingDeadInitC_tab(:)
  real(r8), target, allocatable :: Root1stMaxRadius_tab(:)
  real(r8), target, allocatable :: Root2ndMaxRadius_tab(:)
  real(r8), target, allocatable :: RootPorosity_tab(:)
  real(r8), target, allocatable :: MinNonstC2InitRoot_tab(:)
  real(r8), target, allocatable :: RootRadialResist_tab(:)
  real(r8), target, allocatable :: RootAxialResist_tab(:)
  real(r8), target, allocatable :: ShutRutNonstElmntConducts_tab(:)
  real(r8), target, allocatable :: RootBranchFreq_tab(:)
  real(r8), target, allocatable :: VmaxNH4Root_tab(:)
  real(r8), target, allocatable :: KmNH4Root_tab(:)
  real(r8), target, allocatable :: CMinNH4Root_tab(:)
  real(r8), target, allocatable :: VmaxNO3Root_tab(:)
  real(r8), target, allocatable :: KmNO3Root_tab(:)
  real(r8), target, allocatable :: CminNO3Root_tab(:)
  real(r8), target, allocatable :: VmaxPO4Root_tab(:)
  real(r8), target, allocatable :: KmPO4Root_tab(:)
  real(r8), target, allocatable :: CMinPO4Root_tab(:)
  real(r8), target, allocatable :: CanOsmoPsi0pt_tab(:)
  real(r8), target, allocatable :: RCS_tab(:)
  real(r8), target, allocatable :: CuticleResist_tab(:)
  real(r8), target, allocatable :: LeafBiomGrowthYld_tab(:)
  real(r8), target, allocatable :: PetioleBiomGrowthYld_tab(:)
  real(r8), target, allocatable :: StalkBiomGrowthYld_tab(:)
  real(r8), target, allocatable :: ReserveBiomGrowthYld_tab(:)
  real(r8), target, allocatable :: HuskBiomGrowthYld_tab(:)
  real(r8), target, allocatable :: EarBiomGrowthYld_tab(:)
  real(r8), target, allocatable :: GrainBiomGrowthYld_tab(:)
  real(r8), target, allocatable :: RootBiomGrosYld_tab(:)
  real(r8), target, allocatable :: NoduGrowthYield_tab(:)
  real(r8), target, allocatable :: rNCLeaf_tab(:)
  real(r8), target, allocatable :: rNCSheath_tab(:)
  real(r8), target, allocatable :: rNCStalk_tab(:)
  real(r8), target, allocatable :: rNCReserve_tab(:)
  real(r8), target, allocatable :: rNCHusk_tab(:)
  real(r8), target, allocatable :: rNCEar_tab(:)
  real(r8), target, allocatable :: rNCGrain_tab(:)
  real(r8), target, allocatable :: rNCRoot_tab(:)
  real(r8), target, allocatable :: rNCNodule_tab(:)
  real(r8), target, allocatable :: rPCLeaf_tab(:)
  real(r8), target, allocatable :: rPCSheath_tab(:)
  real(r8), target, allocatable :: rPCStalk_tab(:)
  real(r8), target, allocatable :: rPCReserve_tab(:)
  real(r8), target, allocatable :: rPCHusk_tab(:)
  real(r8), target, allocatable :: rPCEar_tab(:)
  real(r8), target, allocatable :: rPCGrain_tab(:)
  real(r8), target, allocatable :: rPCRootr_tab(:)
  real(r8), target, allocatable :: rPCNoduler_tab(:)
  character(len=10),allocatable :: pftss_tab(:)               !pft trait record name, including short pft-type name and koppen climate zone code
  character(len=40),allocatable :: pft_long_tab(:)            !long pft-type name
  character(len=4), allocatable :: pft_short_tab(:)           !short pft-type name
  character(len=2), allocatable :: koppen_clim_ncode_tab(:)   !numerical koppen climate code 
  character(len=3), allocatable :: koppen_clim_short_tab(:)   !short 3-letter of koppen climate zone code
  character(len=64),allocatable :: koppen_clim_long_tab(:)    !long Description of koppen climate zones

  contains

  Subroutine AllocPlantTraitTable(jroots,npft,nkopenclms,npfts)
  implicit none  
  integer, intent(in) :: jroots          !total root types
  integer, intent(in) :: npft            !total pft types, exclude koppen climate code
  integer, intent(in) :: nkopenclms      !total koppen climate code
  integer, intent(in) :: npfts           !total pft records, pft_short name + numerical koppen climate code

  allocate(iPlantPhotosynthesisType_tab(npfts));iPlantPhotosynthesisType_tab=0
  allocate(iPlantRootProfile_tab(npfts));iPlantRootProfile_tab=0
  allocate(iPlantPhenolPattern_tab(npfts));iPlantPhenolPattern_tab=0
  allocate(iPlantDevelopPattern_tab(npfts));iPlantDevelopPattern_tab=0
  allocate(iPlantNfixType_tab(npfts));iPlantNfixType_tab=0
  allocate(iPlantPhenolType_tab(npfts));iPlantPhenolType_tab=0
  allocate(iPlantPhotoperiodType_tab(npfts));iPlantPhotoperiodType_tab=0
  allocate(iPlantTurnoverPattern_tab(npfts));iPlantTurnoverPattern_tab=0
  allocate(iPlantGrainType_tab(npfts));iPlantGrainType_tab=0
  allocate(Myco_tab(npfts));Myco_tab=0
  allocate(PlantInitThermoAdaptZone_tab(npfts));PlantInitThermoAdaptZone_tab=0.0_r8
  allocate(VmaxRubCarboxyRef_tab(npfts));VmaxRubCarboxyRef_tab=0.0_r8
  allocate(VmaxRubOxyRef_tab(npfts));VmaxRubOxyRef_tab=0.0_r8
  allocate(VmaxPEPCarboxyRef_tab(npfts));VmaxPEPCarboxyRef_tab=0.0_r8
  allocate(XKCO2_tab(npfts));XKCO2_tab=0.0_r8
  allocate(XKO2_tab(npfts));XKO2_tab=0.0_r8
  allocate(Km4PEPCarboxy_tab(npfts));Km4PEPCarboxy_tab=0.0_r8
  allocate(LeafRuBPConc_tab(npfts));LeafRuBPConc_tab=0.0_r8
  allocate(FracLeafProtAsPEPCarboxyl_tab(npfts));FracLeafProtAsPEPCarboxyl_tab=0.0_r8
  allocate(SpecChloryfilAct_tab(npfts));SpecChloryfilAct_tab=0.0_r8
  allocate(LeafC3ChlorofilConc_tab(npfts));LeafC3ChlorofilConc_tab=0._r8
  allocate(LeafC4ChlorofilConc_tab(npfts));LeafC4ChlorofilConc_tab=0._r8
  allocate(CanPCi2CaRatio_tab(npfts));CanPCi2CaRatio_tab=0._r8
  allocate(RadSWLeafAlbedo_tab(npfts));RadSWLeafAlbedo_tab=0._r8
  allocate(CanopyPARalbedo_tab(npfts));CanopyPARalbedo_tab=0._r8
  allocate(RadSWLeafTransmis_tab(npfts));RadSWLeafTransmis_tab=0._r8
  allocate(RadPARLeafTransmis_tab(npfts));RadPARLeafTransmis_tab=0._r8
  allocate(RefNodeInitRate_tab(npfts));RefNodeInitRate_tab=0._r8
  allocate(RefLeafAppearRate_tab(npfts));RefLeafAppearRate_tab=0._r8
  allocate(TCChill4Seed_tab(npfts));TCChill4Seed_tab=0._r8
  allocate(VRNLI_tab(npfts));VRNLI_tab=0._r8
  allocate(VRNXI_tab(npfts));VRNXI_tab=0._r8
  allocate(rLen2WidthLeaf_tab(npfts));rLen2WidthLeaf_tab=0._r8
  allocate(NonstCMinConc2InitBranch_tab(npfts));NonstCMinConc2InitBranch_tab=0._r8
  allocate(GROUPX_tab(npfts));GROUPX_tab=0._r8
  allocate(ShootNodeNumAtPlanting_tab(npfts));ShootNodeNumAtPlanting_tab=0._r8
  allocate(CriticPhotoPeriod_tab(npfts));CriticPhotoPeriod_tab=0._r8
  allocate(PhotoPeriodSens_tab(npfts));PhotoPeriodSens_tab=0._r8
  allocate(SLA1_tab(npfts));SLA1_tab=0._r8
  allocate(PetoLen2Mass_tab(npfts));PetoLen2Mass_tab=0._r8
  allocate(NodeLenPergC_tab(npfts));NodeLenPergC_tab=0._r8
  allocate(LeafAngleClass_tab(1:NumLeafZenithSectors,npfts));LeafAngleClass_tab=0._r8
  allocate(ClumpFactorInit_tab(npfts));ClumpFactorInit_tab=0._r8
  allocate(BranchAngle_tab(npfts));BranchAngle_tab=0._r8
  allocate(PetioleAngle_tab(npfts));PetioleAngle_tab=0._r8
  allocate(MaxPotentSeedNumber_tab(npfts));MaxPotentSeedNumber_tab=0._r8
  allocate(MaxSeedNumPerSite_tab(npfts));MaxSeedNumPerSite_tab=0._r8
  allocate(SeedCMassMax_tab(npfts));SeedCMassMax_tab=0._r8
  allocate(SeedCMass_tab(npfts));SeedCMass_tab=0._r8
  allocate(GrainFillRate25C_tab(npfts));GrainFillRate25C_tab=0._r8
  allocate(StandingDeadInitC_tab(npfts));StandingDeadInitC_tab=0._r8
  allocate(Root1stMaxRadius_tab(npfts));Root1stMaxRadius_tab=0._r8
  allocate(Root2ndMaxRadius_tab(npfts));Root2ndMaxRadius_tab=0._r8
  allocate(RootPorosity_tab(npfts));RootPorosity_tab=0._r8
  allocate(MinNonstC2InitRoot_tab(npfts));MinNonstC2InitRoot_tab=0._r8
  allocate(RootRadialResist_tab(npfts));RootRadialResist_tab=0._r8
  allocate(RootAxialResist_tab(npfts));RootAxialResist_tab=0._r8
  allocate(ShutRutNonstElmntConducts_tab(npfts));ShutRutNonstElmntConducts_tab=0._r8
  allocate(RootBranchFreq_tab(npfts));RootBranchFreq_tab=0._r8
  allocate(VmaxNH4Root_tab(npfts));VmaxNH4Root_tab=0._r8
  allocate(KmNH4Root_tab(npfts));KmNH4Root_tab=0._r8
  allocate(CMinNH4Root_tab(npfts));CMinNH4Root_tab=0._r8
  allocate(VmaxNO3Root_tab(npfts));VmaxNO3Root_tab=0._r8
  allocate(KmNO3Root_tab(npfts));KmNO3Root_tab=0._r8
  allocate(CminNO3Root_tab(npfts));CminNO3Root_tab=0._r8
  allocate(VmaxPO4Root_tab(npfts));VmaxPO4Root_tab=0._r8
  allocate(KmPO4Root_tab(npfts));KmPO4Root_tab=0._r8
  allocate(CMinPO4Root_tab(npfts));CMinPO4Root_tab=0._r8
  allocate(CanOsmoPsi0pt_tab(npfts));CanOsmoPsi0pt_tab=0._r8
  allocate(RCS_tab(npfts));RCS_tab=0._r8
  allocate(CuticleResist_tab(npfts));CuticleResist_tab=0._r8
  allocate(LeafBiomGrowthYld_tab(npfts));LeafBiomGrowthYld_tab=0._r8
  allocate(PetioleBiomGrowthYld_tab(npfts));PetioleBiomGrowthYld_tab=0._r8
  allocate(StalkBiomGrowthYld_tab(npfts));StalkBiomGrowthYld_tab=0._r8
  allocate(ReserveBiomGrowthYld_tab(npfts));ReserveBiomGrowthYld_tab=0._r8
  allocate(HuskBiomGrowthYld_tab(npfts));HuskBiomGrowthYld_tab=0._r8
  allocate(EarBiomGrowthYld_tab(npfts));EarBiomGrowthYld_tab=0._r8
  allocate(GrainBiomGrowthYld_tab(npfts));GrainBiomGrowthYld_tab=0._r8
  allocate(RootBiomGrosYld_tab(npfts));RootBiomGrosYld_tab=0._r8
  allocate(NoduGrowthYield_tab(npfts));NoduGrowthYield_tab=0._r8
  allocate(rNCLeaf_tab(npfts));rNCLeaf_tab=0._r8
  allocate(rNCSheath_tab(npfts));rNCSheath_tab=0._r8
  allocate(rNCStalk_tab(npfts));rNCStalk_tab=0._r8
  allocate(rNCReserve_tab(npfts));rNCReserve_tab=0._r8
  allocate(rNCHusk_tab(npfts));rNCHusk_tab=0._r8
  allocate(rNCEar_tab(npfts));rNCEar_tab=0._r8
  allocate(rNCGrain_tab(npfts));rNCGrain_tab=0._r8
  allocate(rNCRoot_tab(npfts));rNCRoot_tab=0._r8
  allocate(rNCNodule_tab(npfts));rNCNodule_tab=0._r8
  allocate(rPCLeaf_tab(npfts));rPCLeaf_tab=0._r8
  allocate(rPCSheath_tab(npfts));rPCSheath_tab=0._r8
  allocate(rPCStalk_tab(npfts));rPCStalk_tab=0._r8
  allocate(rPCReserve_tab(npfts));rPCReserve_tab=0._r8
  allocate(rPCHusk_tab(npfts));rPCHusk_tab=0._r8
  allocate(rPCEar_tab(npfts));rPCEar_tab=0._r8
  allocate(rPCGrain_tab(npfts));rPCGrain_tab=0._r8
  allocate(rPCRootr_tab(npfts));rPCRootr_tab=0._r8
  allocate(rPCNoduler_tab(npfts));rPCNoduler_tab=0._r8
  allocate(pftss_tab(npfts))
  allocate(pft_long_tab(npft))
  allocate(pft_short_tab(npft))
  allocate(koppen_clim_ncode_tab(nkopenclms))
  allocate(koppen_clim_short_tab(nkopenclms))
  allocate(koppen_clim_long_tab(nkopenclms))
  
  end Subroutine AllocPlantTraitTable


!----------------------------------------------------------------------

  Subroutine DestructPlantTraitTable()
  use abortutils, only : destroy
  implicit none

  call destroy(iPlantPhotosynthesisType_tab)
  call destroy(iPlantRootProfile_tab)
  call destroy(iPlantPhenolPattern_tab)
  call destroy(iPlantDevelopPattern_tab)
  call destroy(iPlantNfixType_tab)
  call destroy(iPlantPhenolType_tab)
  call destroy(iPlantPhotoperiodType_tab)
  call destroy(iPlantTurnoverPattern_tab)
  call destroy(iPlantGrainType_tab)
  call destroy(Myco_tab)
  call destroy(PlantInitThermoAdaptZone_tab)
  call destroy(VmaxRubCarboxyRef_tab)
  call destroy(VmaxRubOxyRef_tab)
  call destroy(VmaxPEPCarboxyRef_tab)
  call destroy(XKCO2_tab)
  call destroy(XKO2_tab)
  call destroy(Km4PEPCarboxy_tab)
  call destroy(LeafRuBPConc_tab)
  call destroy(FracLeafProtAsPEPCarboxyl_tab)
  call destroy(SpecChloryfilAct_tab)
  call destroy(LeafC3ChlorofilConc_tab)
  call destroy(LeafC4ChlorofilConc_tab)
  call destroy(CanPCi2CaRatio_tab)
  call destroy(RadSWLeafAlbedo_tab)
  call destroy(CanopyPARalbedo_tab)
  call destroy(RadSWLeafTransmis_tab)
  call destroy(RadPARLeafTransmis_tab)
  call destroy(RefNodeInitRate_tab)
  call destroy(RefLeafAppearRate_tab)
  call destroy(TCChill4Seed_tab)
  call destroy(VRNLI_tab)
  call destroy(VRNXI_tab)
  call destroy(rLen2WidthLeaf_tab)
  call destroy(NonstCMinConc2InitBranch_tab)
  call destroy(GROUPX_tab)
  call destroy(ShootNodeNumAtPlanting_tab)
  call destroy(CriticPhotoPeriod_tab)
  call destroy(PhotoPeriodSens_tab)
  call destroy(SLA1_tab)
  call destroy(PetoLen2Mass_tab)
  call destroy(NodeLenPergC_tab)  
  call destroy(LeafAngleClass_tab)
  call destroy(ClumpFactorInit_tab)
  call destroy(BranchAngle_tab)
  call destroy(PetioleAngle_tab)
  call destroy(MaxPotentSeedNumber_tab)
  call destroy(MaxSeedNumPerSite_tab)
  call destroy(SeedCMassMax_tab)
  call destroy(SeedCMass_tab)
  call destroy(GrainFillRate25C_tab)
  call destroy(StandingDeadInitC_tab) 
  call destroy(Root1stMaxRadius_tab)
  call destroy(Root2ndMaxRadius_tab)
  call destroy(RootPorosity_tab)
  call destroy(MinNonstC2InitRoot_tab)
  call destroy(RootRadialResist_tab)
  call destroy(RootAxialResist_tab)
  call destroy(ShutRutNonstElmntConducts_tab)
  call destroy(RootBranchFreq_tab)  
  call destroy(VmaxNH4Root_tab)
  call destroy(KmNH4Root_tab)
  call destroy(CMinNH4Root_tab) 
  call destroy(VmaxNO3Root_tab)
  call destroy(KmNO3Root_tab)
  call destroy(CminNO3Root_tab)
  call destroy(VmaxPO4Root_tab)
  call destroy(KmPO4Root_tab)
  call destroy(CMinPO4Root_tab)
  call destroy(CanOsmoPsi0pt_tab)
  call destroy(RCS_tab)
  call destroy(CuticleResist_tab) 
  call destroy(LeafBiomGrowthYld_tab)
  call destroy(PetioleBiomGrowthYld_tab)
  call destroy(StalkBiomGrowthYld_tab)
  call destroy(ReserveBiomGrowthYld_tab)
  call destroy(HuskBiomGrowthYld_tab)
  call destroy(EarBiomGrowthYld_tab)
  call destroy(GrainBiomGrowthYld_tab)
  call destroy(RootBiomGrosYld_tab)
  call destroy(NoduGrowthYield_tab)
  call destroy(rNCLeaf_tab)
  call destroy(rNCSheath_tab)
  call destroy(rNCStalk_tab)
  call destroy(rNCReserve_tab)
  call destroy(rNCHusk_tab)
  call destroy(rNCEar_tab)
  call destroy(rNCGrain_tab)
  call destroy(rNCRoot_tab)
  call destroy(rNCNodule_tab)
  call destroy(rPCLeaf_tab)
  call destroy(rPCSheath_tab)
  call destroy(rPCStalk_tab)
  call destroy(rPCReserve_tab)
  call destroy(rPCHusk_tab)
  call destroy(rPCEar_tab)
  call destroy(rPCGrain_tab)
  call destroy(rPCRootr_tab)
  call destroy(rPCNoduler_tab)
  end Subroutine DestructPlantTraitTable
    
  end module PlantTraitTableMod
