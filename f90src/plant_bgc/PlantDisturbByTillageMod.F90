module PlantDisturbByTillageMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use PlantAPIData
  use ElmIDMod
  use PlantMathFuncMod
  implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: RemoveBiomByTillage
contains
!--------------------------------------------------------------------------------
  subroutine RemoveShootByTillage(I,J,NZ,XHVST,XHVST1)        
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: XHVST,XHVST1

  real(r8) :: WVPLT
  integer :: M,NB,NE,K,L
  real(r8) :: FDM,VOLWPX  
  associate(                                                             &
    jHarvst_pft                 => plt_distb%jHarvst_pft,                 &
    H2OLoss_CumYr_col           => plt_ew%H2OLoss_CumYr_col,              &
    inonstruct                  => pltpar%inonstruct,                     &
    istalk                      => pltpar%istalk,                         &
    inonfoliar                  => pltpar%inonfoliar,                     &
    icwood                      => pltpar%icwood,                         &
    ifoliar                     => pltpar%ifoliar,                        &
    CanopyWater_pft             => plt_ew%CanopyWater_pft,                &
    ElmAllocmat4Litr            => plt_soilchem%ElmAllocmat4Litr,         &
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft,        &
    SeasonalNonstElms_pft       => plt_biom%SeasonalNonstElms_pft,        &
    PSICanopy_pft               => plt_ew%PSICanopy_pft,                  &
    iPlantPhenolPattern_pft     => plt_pheno%iPlantPhenolPattern_pft,     &
    FracShootLeafElmAlloc2Litr  => plt_allom%FracShootLeafElmAlloc2Litr,  &
    FracShootStalkElmAlloc2Litr => plt_allom%FracShootStalkElmAlloc2Litr, &
    iPlantState_pft             => plt_pheno%iPlantState_pft,             &
    iPlantRootState_pft         => plt_pheno%iPlantRootState_pft,         &
    iPlantShootState_pft        => plt_pheno%iPlantShootState_pft,        &
    CMassHCO3BundleSheath_node  => plt_photo%CMassHCO3BundleSheath_node,  &
    CMassCO2BundleSheath_node   => plt_photo%CMassCO2BundleSheath_node,   &
    FracRootStalkElmAlloc2Litr  => plt_allom%FracRootStalkElmAlloc2Litr,  &
    CPOOL3_node                 => plt_photo%CPOOL3_node,                 &
    CPOOL4_node                 => plt_photo%CPOOL4_node,                 &
    QH2OLoss_lnds               => plt_site%QH2OLoss_lnds,                &
    k_fine_litr                 => pltpar%k_fine_litr,                    &
    k_woody_litr                => pltpar%k_woody_litr,                   &
    iYearCurrent                => plt_site%iYearCurrent,                 &
    iDayPlantHarvest_pft        => plt_distb%iDayPlantHarvest_pft,        &
    iYearPlantHarvest_pft       => plt_distb%iYearPlantHarvest_pft,       &
    LeafProteinCNode_brch       => plt_biom%LeafProteinCNode_brch,        &
    InternodeStrutElms_brch     => plt_biom%InternodeStrutElms_brch,      &
    CanopyLeafShethC_pft        => plt_biom%CanopyLeafShethC_pft,         &
    FracPARads2Canopy_pft      => plt_rad%FracPARads2Canopy_pft,        &
    CanopyStalkC_pft            => plt_biom%CanopyStalkC_pft,             &
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft,           &
    iPlantBranchState_brch      => plt_pheno%iPlantBranchState_brch,      &
    PlantPopulation_pft         => plt_site%PlantPopulation_pft,          &
    LitrfalStrutElms_pvr        => plt_bgcr%LitrfalStrutElms_pvr,         &
    CanopyNodulNonstElms_brch   => plt_biom%CanopyNodulNonstElms_brch,    &
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch,         &
    ShootC4NonstC_brch          => plt_biom%ShootC4NonstC_brch,           &
    VHeatCapCanP_pft            => plt_ew%VHeatCapCanP_pft,               &
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch,           &
    GrainStrutElms_brch         => plt_biom%GrainStrutElms_brch,          &
    EarStrutElms_brch           => plt_biom%EarStrutElms_brch,            &
    HuskStrutElms_brch          => plt_biom%HuskStrutElms_brch,           &
    StalkRsrvElms_brch          => plt_biom%StalkRsrvElms_brch,           &
    CanopyNodulStrutElms_brch   => plt_biom%CanopyNodulStrutElms_brch,    &
    GrainSeedBiomCMean_brch     => plt_allom%GrainSeedBiomCMean_brch,     &
    ShootStrutElms_brch         => plt_biom%ShootStrutElms_brch,          &
    StalkStrutElms_brch         => plt_biom%StalkStrutElms_brch,          &
    LeafElmsByLayerNode_brch    => plt_biom%LeafElmsByLayerNode_brch,     &
    PetoleStrutElms_brch        => plt_biom%PetoleStrutElms_brch,         &
    LeafElmntNode_brch          => plt_biom%LeafElmntNode_brch,           &
    LeafPetolBiomassC_brch      => plt_biom%LeafPetolBiomassC_brch,       &
    SenecStalkStrutElms_brch    => plt_biom%SenecStalkStrutElms_brch,     &
    StalkBiomassC_brch          => plt_biom%StalkBiomassC_brch,           &
    LeafAreaNode_brch           => plt_morph%LeafAreaNode_brch,           &
    LeafAreaLive_brch           => plt_morph%LeafAreaLive_brch,           &
    PotentialSeedSites_brch     => plt_morph%PotentialSeedSites_brch,     &
    SeedNumSet_brch             => plt_morph%SeedNumSet_brch,             &
    PetoleProteinCNode_brch     => plt_biom%PetoleProteinCNode_brch,      &
    PetioleElmntNode_brch       => plt_biom%PetioleElmntNode_brch,        &
    CanopyLeafArea_lpft         => plt_morph%CanopyLeafArea_lpft          &
  )

!     PPX,PP=PFT population per m2,grid cell
!     FracPARads2Canopy_pft=fraction of radiation received by each PFT canopy
!     VHeatCapCanP_pft=canopy heat capacity

  FracPARads2Canopy_pft(NZ)=FracPARads2Canopy_pft(NZ)*XHVST
  VHeatCapCanP_pft(NZ)=VHeatCapCanP_pft(NZ)*XHVST
  CanopyLeafShethC_pft(NZ)=0._r8
  CanopyStalkC_pft(NZ)=0._r8
!
!     TERMINATE BRANCHES IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     PP=PFT population
!
  D8975: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
      IF(PlantPopulation_pft(NZ).LE.0.0)then
        iPlantBranchState_brch(NB,NZ)=iDead
      endif
!
!     LitrFall FROM BRANCHES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     XHVST=fraction of PFT remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     ShootC4NonstC_brch=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
    D6380: DO M=1,jsken
      LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,M,k_fine_litr,0,NZ) &
        +XHVST1*(ElmAllocmat4Litr(ielmc,inonstruct,M,NZ)*(CanopyNonstElms_brch(ielmc,NB,NZ) &
        +CanopyNodulNonstElms_brch(ielmc,NB,NZ) &
        +ShootC4NonstC_brch(NB,NZ)+StalkRsrvElms_brch(ielmc,NB,NZ)) &
        +ElmAllocmat4Litr(ielmc,ifoliar,M,NZ)*(LeafStrutElms_brch(ielmc,NB,NZ)*FracShootStalkElmAlloc2Litr(ielmc,k_fine_litr) &
        +CanopyNodulStrutElms_brch(ielmc,NB,NZ)) &
        +ElmAllocmat4Litr(ielmc,inonfoliar,M,NZ)*(PetoleStrutElms_brch(ielmc,NB,NZ)*FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr) &
        +HuskStrutElms_brch(ielmc,NB,NZ)+EarStrutElms_brch(ielmc,NB,NZ)))

      DO NE=2,NumPlantChemElms
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
          *(ElmAllocmat4Litr(NE,inonstruct,M,NZ)*(CanopyNonstElms_brch(NE,NB,NZ)+CanopyNodulNonstElms_brch(NE,NB,NZ)&
          +StalkRsrvElms_brch(NE,NB,NZ)) &
          +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_fine_litr) &
          +CanopyNodulStrutElms_brch(NE,NB,NZ)) &
          +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*(PetoleStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_fine_litr) &
          +HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ)))
      ENDDO
    ENDDO D6380

    DO M=1,jsken
      DO NE=1,NumPlantChemElms
        LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+XHVST1 &
          *ElmAllocmat4Litr(NE,icwood,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootStalkElmAlloc2Litr(NE,k_woody_litr) &
          +PetoleStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_woody_litr))

        IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+XHVST1 &
            *ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
        ELSE
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
            *ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
        ENDIF
        LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+XHVST1 &
          *ElmAllocmat4Litr(NE,icwood,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

        LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
          *ElmAllocmat4Litr(NE,istalk,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
      ENDDO
    ENDDO
  !
  !     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
  !

    ShootC4NonstC_brch(NB,NZ)=ShootC4NonstC_brch(NB,NZ)*XHVST
    StalkBiomassC_brch(NB,NZ)=StalkBiomassC_brch(NB,NZ)*XHVST
    DO NE=1,NumPlantChemElms
      CanopyNonstElms_brch(NE,NB,NZ)      = CanopyNonstElms_brch(NE,NB,NZ)*XHVST
      CanopyNodulNonstElms_brch(NE,NB,NZ) = CanopyNodulNonstElms_brch(NE,NB,NZ)*XHVST
      ShootStrutElms_brch(NE,NB,NZ)       = ShootStrutElms_brch(NE,NB,NZ)*XHVST
      StalkRsrvElms_brch(NE,NB,NZ)        = StalkRsrvElms_brch(NE,NB,NZ)*XHVST
      HuskStrutElms_brch(NE,NB,NZ)        = HuskStrutElms_brch(NE,NB,NZ)*XHVST
      EarStrutElms_brch(NE,NB,NZ)         = EarStrutElms_brch(NE,NB,NZ)*XHVST
      GrainStrutElms_brch(NE,NB,NZ)       = GrainStrutElms_brch(NE,NB,NZ)*XHVST
      LeafStrutElms_brch(NE,NB,NZ)        = LeafStrutElms_brch(NE,NB,NZ)*XHVST
      CanopyNodulStrutElms_brch(NE,NB,NZ) = CanopyNodulStrutElms_brch(NE,NB,NZ)*XHVST
      PetoleStrutElms_brch(NE,NB,NZ)      = PetoleStrutElms_brch(NE,NB,NZ)*XHVST
      StalkStrutElms_brch(NE,NB,NZ)       = StalkStrutElms_brch(NE,NB,NZ)*XHVST
      SenecStalkStrutElms_brch(NE,NB,NZ)  = SenecStalkStrutElms_brch(NE,NB,NZ)*XHVST
    ENDDO

    PotentialSeedSites_brch(NB,NZ)=PotentialSeedSites_brch(NB,NZ)*XHVST
    SeedNumSet_brch(NB,NZ)=SeedNumSet_brch(NB,NZ)*XHVST
    GrainSeedBiomCMean_brch(NB,NZ)=GrainSeedBiomCMean_brch(NB,NZ)*XHVST
    LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)*XHVST
    LeafPetolBiomassC_brch(NB,NZ)=AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetoleStrutElms_brch(ielmc,NB,NZ))
    CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(NB,NZ)

    CanopyStalkC_pft(NZ)=CanopyStalkC_pft(NZ)+StalkBiomassC_brch(NB,NZ)
    D8970: DO K=0,MaxNodesPerBranch1
      IF(K.NE.0)THEN
        CPOOL3_node(K,NB,NZ)=CPOOL3_node(K,NB,NZ)*XHVST
        CPOOL4_node(K,NB,NZ)=CPOOL4_node(K,NB,NZ)*XHVST
        CMassCO2BundleSheath_node(K,NB,NZ)=CMassCO2BundleSheath_node(K,NB,NZ)*XHVST
        CMassHCO3BundleSheath_node(K,NB,NZ)=CMassHCO3BundleSheath_node(K,NB,NZ)*XHVST
      ENDIF
      LeafAreaNode_brch(K,NB,NZ)=LeafAreaNode_brch(K,NB,NZ)*XHVST

      LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)*XHVST
  !     PetoleLensNode_brch(K,NB,NZ)=PetoleLensNode_brch(K,NB,NZ)*XHVST

      PetoleProteinCNode_brch(K,NB,NZ)=PetoleProteinCNode_brch(K,NB,NZ)*XHVST
  !     LiveInterNodeHight_brch(K,NB,NZ)=LiveInterNodeHight_brch(K,NB,NZ)*XHVST
  !     InternodeHeightDying_brch(K,NB,NZ)=InternodeHeightDying_brch(K,NB,NZ)*XHVST
      DO NE=1,NumPlantChemElms
        InternodeStrutElms_brch(NE,K,NB,NZ) = InternodeStrutElms_brch(NE,K,NB,NZ)*XHVST
        LeafElmntNode_brch(NE,K,NB,NZ)      = LeafElmntNode_brch(NE,K,NB,NZ)*XHVST
        PetioleElmntNode_brch(NE,K,NB,NZ)   = PetioleElmntNode_brch(NE,K,NB,NZ)*XHVST
        DO L=1,NumOfCanopyLayers1
          LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)=LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)*XHVST
        ENDDO
      ENDDO
      D8965: DO L=1,NumOfCanopyLayers1
        CanopyLeafArea_lpft(L,K,NB,NZ)=CanopyLeafArea_lpft(L,K,NB,NZ)*XHVST
      ENDDO D8965
    ENDDO D8970
  ENDIF
  ENDDO D8975
!
!     PSICanopy_pft=canopy water potential
!     CanopyWater_pft=water volume in canopy
!     QH2OLoss_lnds,H2OLoss_CumYr_col=accumulated water loss for water balance calculation
!
  VOLWPX=CanopyWater_pft(NZ)
  WVPLT=AZMAX1(CanopyLeafShethC_pft(NZ)+CanopyStalkC_pft(NZ))

  FDM=get_FDM(PSICanopy_pft(NZ))
!        APSILT=ABS(PSICanopy_pft(NZ))
!        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)

  CanopyWater_pft(NZ)=ppmc*WVPLT/FDM
  QH2OLoss_lnds=QH2OLoss_lnds+VOLWPX-CanopyWater_pft(NZ)
  H2OLoss_CumYr_col=H2OLoss_CumYr_col+VOLWPX-CanopyWater_pft(NZ)
!
!     TERMINATE ROOTS IF TILLAGE IMPLEMENT 10 IS SELECTED
!
!     PP=PFT population
!     IDTHR,iPlantShootState_pft=PFT root,shoot living flag: 0=alive,1=dead
!     IDTH=PFT living flag: 0=alive,1=dead
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,and reseed
!     iDayPlantHarvest_pft,iYearPlantHarvest_pft=day,year of harvesting
!     iYearCurrent=current year
!
  IF(PlantPopulation_pft(NZ).LE.0.0_r8)THEN
    iPlantRootState_pft(NZ)   = iDead
    iPlantShootState_pft(NZ)  = iDead
    iPlantState_pft(NZ)       = iDead
    jHarvst_pft(NZ)           = jharvtyp_terminate
    iDayPlantHarvest_pft(NZ)  = I
    iYearPlantHarvest_pft(NZ) = iYearCurrent
  ENDIF
  end associate
  END subroutine RemoveShootByTillage


!------------------------------------------------------------------------------------------    
  subroutine RemoveBiomByTillage(I,J,NZ)
  !
  !Description
  !REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
  implicit none
  integer , intent(in) :: I,J,NZ
  integer :: L,K,M,N,NR,NE,NB,NTG
  real(r8) :: XHVST,XHVST1
  REAL(R8) :: APSILT
!     begin_execution
  associate(                                                              &
    iSoilDisturbType_col        => plt_distb%iSoilDisturbType_col,        &
    iPlantTurnoverPattern_pft   => plt_pheno%iPlantTurnoverPattern_pft,   &    
    iDayPlanting_pft            => plt_distb%iDayPlanting_pft,            &  
    iYearCurrent                => plt_site%iYearCurrent,                 &      
    iYearPlanting_pft           => plt_distb%iYearPlanting_pft,           &    
    iPlantRootProfile_pft       => plt_pheno%iPlantRootProfile_pft,       &  
    XCORP                       => plt_distb%XCORP,                       &    
    PPX_pft                     => plt_site%PPX_pft,                      &    
    PlantPopulation_pft         => plt_site%PlantPopulation_pft,          &      
    SolarNoonHour_col           => plt_site%SolarNoonHour_col             &    
  )
!     SolarNoonHour_col=hour of solar noon
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     iDayPlanting_pft,iYearPlanting_pft=day,year of planting
!     iYearCurrent=current year
!     iSoilDisturbType_col=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     XHVST=fraction of PFT remaining after disturbance
!
  IF(J.EQ.INT(SolarNoonHour_col) .AND. (iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR. (.not.is_plant_treelike(iPlantRootProfile_pft(NZ)))) &
    .AND. (I.NE.iDayPlanting_pft(NZ) .OR. iYearCurrent.NE.iYearPlanting_pft(NZ)))THEN

    IF(iSoilDisturbType_col.LE.10 .OR. NZ.NE.1)THEN
      IF(I.GT.iDayPlanting_pft(NZ) .OR. iYearCurrent.GT.iYearPlanting_pft(NZ))THEN
        XHVST=XCORP
        PPX_pft(NZ)=PPX_pft(NZ)*XHVST
        PlantPopulation_pft(NZ)=PlantPopulation_pft(NZ)*XHVST
        XHVST1=1._r8-XHVST
        
        call RemoveShootByTillage(I,J,NZ,XHVST,XHVST1)

!     LitrFall FROM ROOTS DURING TILLAGE
        call RemoveRootsByTillage(I,J,NZ,XHVST,XHVST1)

      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByTillage
!------------------------------------------------------------------------------------------
  subroutine RemoveRootsByTillage(I,J,NZ,XHVST,XHVST1)
  implicit none        
  integer, intent(in) :: I,J,NZ
  real(r8),intent(in) :: XHVST,XHVST1
  integer :: L,N,M,NE,NR,NTG

  associate(                                                            &
    MaxNumRootLays             => plt_site%MaxNumRootLays,              &
    MY                         => plt_morph%MY,                         &
    ElmAllocmat4Litr           => plt_soilchem%ElmAllocmat4Litr,        &
    trcg_rootml_pvr            => plt_rbgc%trcg_rootml_pvr,             &
    trcs_rootml_pvr            => plt_rbgc%trcs_rootml_pvr,             &
    RootMycoNonstElms_rpvr     => plt_biom%RootMycoNonstElms_rpvr,      &
    RootProteinC_pvr           => plt_biom%RootProteinC_pvr,            &
    PopuRootMycoC_pvr          => plt_biom% PopuRootMycoC_pvr,          &
    RootMycoActiveBiomC_pvr    => plt_biom%RootMycoActiveBiomC_pvr,     &
    RootMyco1stStrutElms_rpvr  => plt_biom%RootMyco1stStrutElms_rpvr,   &
    RootMyco1stElm_raxs        => plt_biom%RootMyco1stElm_raxs,         &
    SeasonalNonstElms_pft      => plt_biom%SeasonalNonstElms_pft,       &
    RootMyco2ndStrutElms_rpvr  => plt_biom%RootMyco2ndStrutElms_rpvr,   &
    RootNodulNonstElms_pvr     => plt_biom%RootNodulNonstElms_pvr,      &
    RootNodulStrutElms_pvr     => plt_biom%RootNodulStrutElms_pvr,      &
    FracRootStalkElmAlloc2Litr => plt_allom%FracRootStalkElmAlloc2Litr, &
    FracRootElmAlloc2Litr      => plt_allom%FracRootElmAlloc2Litr,      &
    QH2OLoss_lnds              => plt_site%QH2OLoss_lnds,               &
    NU                         => plt_site%NU,                          &
    k_fine_litr                => pltpar%k_fine_litr,                   &
    k_woody_litr               => pltpar%k_woody_litr,                  &
    inonstruct                 => pltpar%inonstruct,                    &
    iroot                      => pltpar%iroot,                         &
    icwood                     => pltpar%icwood,                        &
    RootGasLossDisturb_pft     => plt_bgcr%RootGasLossDisturb_pft,      &
    LitrfalStrutElms_pvr       => plt_bgcr%LitrfalStrutElms_pvr,        &
    RootCO2Autor_pvr           => plt_rbgc%RootCO2Autor_pvr,            &
    RootRespPotent_pvr         => plt_rbgc%RootRespPotent_pvr,          &
    RootCO2EmisPot_pvr         => plt_rbgc%RootCO2EmisPot_pvr,          &
    Root1stLen_rpvr            => plt_morph%Root1stLen_rpvr,            &
    RootVH2O_pvr               => plt_morph%RootVH2O_pvr,               &
    RootAreaPerPlant_pvr       => plt_morph%RootAreaPerPlant_pvr,       &
    RootPoreVol_pvr            => plt_morph%RootPoreVol_pvr,            &
    RootLenDensPerPlant_pvr    => plt_morph%RootLenDensPerPlant_pvr,    &
    Root1stXNumL_pvr           => plt_morph%Root1stXNumL_pvr,           &
    iPlantNfixType_pft         => plt_morph%iPlantNfixType_pft,         &
    RootLenPerPlant_pvr        => plt_morph%RootLenPerPlant_pvr,        &
    Root2ndXNum_rpvr           => plt_morph%Root2ndXNum_rpvr,           &
    Root2ndLen_pvr             => plt_morph%Root2ndLen_pvr,             &
    Root2ndXNum_pvr            => plt_morph%Root2ndXNum_pvr,            &
    NGTopRootLayer_pft         => plt_morph%NGTopRootLayer_pft,         &
    NumRootAxes_pft            => plt_morph%NumRootAxes_pft             &
  )
!
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     XHVST=fraction of PFT remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!
  
  DO NR=1,NumRootAxes_pft(NZ)
    DO N=1,MY(NZ)
      DO NE=1,NumPlantChemElms
        RootMyco1stElm_raxs(NE,N,NR,NZ)=RootMyco1stElm_raxs(NE,N,NR,NZ)*XHVST
      ENDDO
    ENDDO
  ENDDO

  D8980: DO L=NU,MaxNumRootLays
    D8985: DO N=1,MY(NZ)
      D6385: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+XHVST1 &
              *ElmAllocmat4Litr(NE,inonstruct,M,NZ)* RootMycoNonstElms_rpvr(NE,N,L,NZ)
          ENDDO

        DO NR=1,NumRootAxes_pft(NZ)
          DO NE=1,NumPlantChemElms
            LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,L,NZ)+XHVST1 &
              *ElmAllocmat4Litr(NE,icwood,M,NZ)*(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) &
              +RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))*FracRootElmAlloc2Litr(NE,k_woody_litr)

            LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+XHVST1 &
              *ElmAllocmat4Litr(NE,iroot,M,NZ)*(RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) &
              +RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ))*FracRootElmAlloc2Litr(NE,k_fine_litr)
          ENDDO
        ENDDO
      ENDDO D6385
!
!     RELEASE ROOT GAS CONTENTS DURING TILLAGE
!
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!

      DO NTG=idg_beg,idg_end-1
        RootGasLossDisturb_pft(NTG,NZ)=RootGasLossDisturb_pft(NTG,NZ)-XHVST1 &
          *(trcg_rootml_pvr(NTG,N,L,NZ)+trcs_rootml_pvr(NTG,N,L,NZ))
        trcg_rootml_pvr(NTG,N,L,NZ)=XHVST*trcg_rootml_pvr(NTG,N,L,NZ)
        trcs_rootml_pvr(NTG,N,L,NZ)=XHVST*trcs_rootml_pvr(NTG,N,L,NZ)
      ENDDO
!
!     ROOT STATE VARIABLES REMAINING AFTER TILLAGE
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     Root1stLen_rpvr,Root2ndLen_pvr=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RootMycoActiveBiomC_pvr, PopuRootMycoC_pvr=active,actual root C mass
!     RootProteinC_pvr=root protein C mass
!     RTN1,Root2ndXNum_pvr=number of primary,secondary root axes
!     RootLenDensPerPlant_pvr,RootLenPerPlant_pvr=root length density,root length per plant
!     RootVH2O_pvr,RootPoreVol_pvr=root or myco aqueous,gaseous volume
!     RootAreaPerPlant_pvr=root surface area per plant
!     RootRespPotent_pvr,RootCO2EmisPot_pvr,RootCO2Autor_pvr unlimited by O2,nonstructural C
!
      D8960: DO NR=1,NumRootAxes_pft(NZ)
        DO NE=1,NumPlantChemElms
          RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)*XHVST
          RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)=RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)*XHVST
        ENDDO
        Root1stLen_rpvr(N,L,NR,NZ)=Root1stLen_rpvr(N,L,NR,NZ)*XHVST
        Root2ndLen_pvr(N,L,NR,NZ)=Root2ndLen_pvr(N,L,NR,NZ)*XHVST
        Root2ndXNum_rpvr(N,L,NR,NZ)=Root2ndXNum_rpvr(N,L,NR,NZ)*XHVST
      ENDDO D8960

      DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)*XHVST
      ENDDO
      RootMycoActiveBiomC_pvr(N,L,NZ) = RootMycoActiveBiomC_pvr(N,L,NZ)*XHVST
      PopuRootMycoC_pvr(N,L,NZ)       = PopuRootMycoC_pvr(N,L,NZ)*XHVST
      RootProteinC_pvr(N,L,NZ)        = RootProteinC_pvr(N,L,NZ)*XHVST
      Root1stXNumL_pvr(N,L,NZ)        = Root1stXNumL_pvr(N,L,NZ)*XHVST
      Root2ndXNum_pvr(N,L,NZ)         = Root2ndXNum_pvr(N,L,NZ)*XHVST
      RootLenPerPlant_pvr(N,L,NZ)     = RootLenPerPlant_pvr(N,L,NZ)*XHVST
      RootLenDensPerPlant_pvr(N,L,NZ) = RootLenDensPerPlant_pvr(N,L,NZ)*XHVST
      RootPoreVol_pvr(N,L,NZ)         = RootPoreVol_pvr(N,L,NZ)*XHVST
      RootVH2O_pvr(N,L,NZ)            = RootVH2O_pvr(N,L,NZ)*XHVST
      RootAreaPerPlant_pvr(N,L,NZ)    = RootAreaPerPlant_pvr(N,L,NZ)*XHVST
      RootRespPotent_pvr(N,L,NZ)      = RootRespPotent_pvr(N,L,NZ)*XHVST
      RootCO2EmisPot_pvr(N,L,NZ)      = RootCO2EmisPot_pvr(N,L,NZ)*XHVST
      RootCO2Autor_pvr(N,L,NZ)        = RootCO2Autor_pvr(N,L,NZ)*XHVST
!
!     LitrFall AND STATE VARIABLES FOR NODULES DURING TILLAGE
!
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
      IF(is_plant_N2fix(iPlantNfixType_pft(NZ)).AND.N.EQ.ipltroot)THEN
        DO NE=1,NumPlantChemElms
          D6395: DO M=1,jsken
            LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,L,NZ)+&
              XHVST1*(ElmAllocmat4Litr(NE,iroot,M,NZ)*RootNodulStrutElms_pvr(NE,L,NZ) &
              +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*RootNodulNonstElms_pvr(NE,L,NZ))
          ENDDO D6395
          RootNodulStrutElms_pvr(NE,L,NZ)=RootNodulStrutElms_pvr(NE,L,NZ)*XHVST
          RootNodulNonstElms_pvr(NE,L,NZ)=RootNodulNonstElms_pvr(NE,L,NZ)*XHVST
        ENDDO
      ENDIF
    ENDDO D8985
  ENDDO D8980
!
!     LitrFall AND STATE VARIABLES FOR SEASONAL STORAGE RESERVES
!     DURING TILLAGE
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining after disturbance
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
  DO NE=1,NumPlantChemElms
    D6400: DO M=1,jsken
      LitrfalStrutElms_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ)=&
        LitrfalStrutElms_pvr(NE,M,k_woody_litr,NGTopRootLayer_pft(NZ),NZ) &
        +(XHVST1*ElmAllocmat4Litr(NE,inonstruct,M,NZ)*SeasonalNonstElms_pft(NE,NZ)) &
        *FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

      LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ)= &
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ) &
        +(XHVST1*ElmAllocmat4Litr(NE,inonstruct,M,NZ)*SeasonalNonstElms_pft(NE,NZ)) &
        *FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
    ENDDO D6400
    SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)*XHVST
  ENDDO
  end associate
  end subroutine RemoveRootsByTillage        

end module PlantDisturbByTillageMod