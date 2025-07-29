module PlantDisturbByTillageMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use DebugToolMod
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
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine RemoveShootByTillage(I,J,NZ,XHVST)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: XHVST
  
  character(len=*), parameter :: subname='RemoveShootByTillage'
  real(r8) :: XHVST1
  real(r8) :: WVPLT
  integer :: M,NB,NE,K,L
  real(r8) :: FDM,VOLWPX  
  associate(                                                               &
    inonstruct                  => pltpar%inonstruct                      ,& !input  :group id of plant nonstructural litter
    istalk                      => pltpar%istalk                          ,& !input  :group id of plant stalk litter group
    inonfoliar                  => pltpar%inonfoliar                      ,& !input  :group id of plant non-foliar litter group
    icwood                      => pltpar%icwood                          ,& !input  :group id of coarse woody litter
    ifoliar                     => pltpar%ifoliar                         ,& !input  :group id of plant foliar litter
    ElmAllocmat4Litr            => plt_soilchem%ElmAllocmat4Litr          ,& !input  :litter kinetic fraction, [-]
    iPlantPhenolType_pft        => plt_pheno%iPlantPhenolType_pft         ,& !input  :climate signal for phenological progress: none, temperature, water stress,[-]
    PSICanopy_pft               => plt_ew%PSICanopy_pft                   ,& !input  :canopy total water potential, [Mpa]
    iPlantPhenolPattern_pft     => plt_pheno%iPlantPhenolPattern_pft      ,& !input  :plant growth habit: annual or perennial,[-]
    FracShootLeafElmAlloc2Litr  => plt_allom%FracShootLeafElmAlloc2Litr   ,& !input  :woody element allocation, [-]
    FracShootPetolElmAlloc2Litr => plt_allom%FracShootPetolElmAlloc2Litr  ,& !input  :leaf element allocation,[-]
    FracWoodStalkElmAlloc2Litr  => plt_allom%FracWoodStalkElmAlloc2Litr   ,& !input  :woody element allocation,[-]
    k_fine_litr                 => pltpar%k_fine_litr                     ,& !input  :fine litter complex id
    k_woody_litr                => pltpar%k_woody_litr                    ,& !input  :woody litter complex id
    iYearCurrent                => plt_site%iYearCurrent                  ,& !input  :current year,[-]
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft            ,& !input  :number of branches,[-]
    PlantPopulation_pft         => plt_site%PlantPopulation_pft           ,& !input  :plant population, [d-2]
    H2OLoss_CumYr_col           => plt_ew%H2OLoss_CumYr_col               ,& !inoput :total subsurface water flux, [m3 d-2]
    CanopyBiomWater_pft         => plt_ew%CanopyBiomWater_pft             ,& !inoput :canopy water content, [m3 d-2]
    SeasonalNonstElms_pft       => plt_biom%SeasonalNonstElms_pft         ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    CMassHCO3BundleSheath_node  => plt_photo%CMassHCO3BundleSheath_node   ,& !inoput :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CMassCO2BundleSheath_node   => plt_photo%CMassCO2BundleSheath_node    ,& !inoput :bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
    CPOOL3_node                 => plt_photo%CPOOL3_node                  ,& !inoput :minimum sink strength for nonstructural C transfer, [g d-2]
    CPOOL4_node                 => plt_photo%CPOOL4_node                  ,& !inoput :leaf nonstructural C4 content in C4 photosynthesis, [g d-2]
    QH2OLoss_lnds               => plt_site%QH2OLoss_lnds                 ,& !inoput :total subsurface water loss flux over the landscape, [m3 d-2]
    LeafProteinCNode_brch       => plt_biom%LeafProteinCNode_brch         ,& !inoput :layer leaf protein C, [g d-2]
    InternodeStrutElms_brch     => plt_biom%InternodeStrutElms_brch       ,& !inoput :internode C, [g d-2]
    CanopyLeafShethC_pft        => plt_biom%CanopyLeafShethC_pft          ,& !inoput :canopy leaf + sheath C, [g d-2]
    FracPARads2Canopy_pft       => plt_rad%FracPARads2Canopy_pft          ,& !inoput :fraction of incoming PAR absorbed by canopy, [-]
    CanopySapwoodC_pft            => plt_biom%CanopySapwoodC_pft              ,& !inoput :canopy active stalk C, [g d-2]
    iPlantBranchState_brch      => plt_pheno%iPlantBranchState_brch       ,& !inoput :flag to detect branch death, [-]
    LitrfalStrutElms_pvr        => plt_bgcr%LitrfalStrutElms_pvr          ,& !inoput :plant LitrFall element, [g d-2 h-1]
    CanopyNodulNonstElms_brch   => plt_biom%CanopyNodulNonstElms_brch     ,& !inoput :branch nodule nonstructural element, [g d-2]
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch          ,& !inoput :branch nonstructural element, [g d-2]
    C4PhotoShootNonstC_brch     => plt_biom%C4PhotoShootNonstC_brch       ,& !inoput :branch shoot nonstrucal elelment, [g d-2]
    VHeatCapCanopy_pft          => plt_ew%VHeatCapCanopy_pft              ,& !inoput :canopy heat capacity, [MJ d-2 K-1]
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch            ,& !inoput :branch leaf structural element mass, [g d-2]
    GrainStrutElms_brch         => plt_biom%GrainStrutElms_brch           ,& !inoput :branch grain structural element mass, [g d-2]
    EarStrutElms_brch           => plt_biom%EarStrutElms_brch             ,& !inoput :branch ear structural chemical element mass, [g d-2]
    HuskStrutElms_brch          => plt_biom%HuskStrutElms_brch            ,& !inoput :branch husk structural element mass, [g d-2]
    StalkRsrvElms_brch          => plt_biom%StalkRsrvElms_brch            ,& !inoput :branch reserve element mass, [g d-2]
    CanopyNodulStrutElms_brch   => plt_biom%CanopyNodulStrutElms_brch     ,& !inoput :branch nodule structural element, [g d-2]
    GrainSeedBiomCMean_brch     => plt_allom%GrainSeedBiomCMean_brch      ,& !inoput :maximum grain C during grain fill, [g d-2]
    ShootStrutElms_brch         => plt_biom%ShootStrutElms_brch           ,& !inoput :branch shoot structural element mass, [g d-2]
    StalkStrutElms_brch         => plt_biom%StalkStrutElms_brch           ,& !inoput :branch stalk structural element mass, [g d-2]
    LeafElmsByLayerNode_brch    => plt_biom%LeafElmsByLayerNode_brch      ,& !inoput :layer leaf element, [g d-2]
    PetoleStrutElms_brch        => plt_biom%PetoleStrutElms_brch          ,& !inoput :branch sheath structural element, [g d-2]
    LeafElmntNode_brch          => plt_biom%LeafElmntNode_brch            ,& !inoput :leaf element, [g d-2]
    SenecStalkStrutElms_brch    => plt_biom%SenecStalkStrutElms_brch      ,& !inoput :branch stalk structural element, [g d-2]
    SapwoodBiomassC_brch      => plt_biom%SapwoodBiomassC_brch        ,& !inoput :branch live stalk C, [gC d-2]
    LeafArea_node           => plt_morph%LeafArea_node            ,& !inoput :leaf area, [m2 d-2]
    LeafAreaLive_brch           => plt_morph%LeafAreaLive_brch            ,& !inoput :branch leaf area, [m2 d-2]
    PotentialSeedSites_brch     => plt_morph%PotentialSeedSites_brch      ,& !inoput :branch potential grain number, [d-2]
    SeedNumSet_brch             => plt_morph%SeedNumSet_brch              ,& !inoput :branch grain number, [d-2]
    PetoleProteinCNode_brch     => plt_biom%PetoleProteinCNode_brch       ,& !inoput :layer sheath protein C, [g d-2]
    PetioleElmntNode_brch       => plt_biom%PetioleElmntNode_brch         ,& !inoput :sheath chemical element, [g d-2]
    CanopyLeafArea_lnode        => plt_morph%CanopyLeafArea_lnode         ,& !inoput :layer/node/branch leaf area, [m2 d-2]
    jHarvst_pft                 => plt_distb%jHarvst_pft                  ,& !output :flag for stand replacing disturbance,[-]
    iPlantState_pft             => plt_pheno%iPlantState_pft              ,& !output :flag for species death, [-]
    iPlantRootState_pft         => plt_pheno%iPlantRootState_pft          ,& !output :flag to detect root system death,[-]
    iPlantShootState_pft        => plt_pheno%iPlantShootState_pft         ,& !output :flag to detect canopy death,[-]
    iDayPlantHarvest_pft        => plt_distb%iDayPlantHarvest_pft         ,& !output :day of harvest,[-]
    iYearPlantHarvest_pft       => plt_distb%iYearPlantHarvest_pft        ,& !output :year of harvest,[-]
    LeafPetolBiomassC_brch      => plt_biom%LeafPetolBiomassC_brch         & !output :plant branch leaf + sheath C, [g d-2]
  )
  call PrintInfo('beg '//subname)
!     PPX,PP=PFT population per m2,grid cell
!     FracPARads2Canopy_pft=fraction of radiation received by each PFT canopy
!     VHeatCapCanopy_pft=canopy heat capacity

  XHVST1                    = 1._r8-XHVST
  FracPARads2Canopy_pft(NZ) = FracPARads2Canopy_pft(NZ)*XHVST
  VHeatCapCanopy_pft(NZ)    = VHeatCapCanopy_pft(NZ)*XHVST
  CanopyLeafShethC_pft(NZ)  = 0._r8
  CanopySapwoodC_pft(NZ)      = 0._r8
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
!     XHVST=fraction of PFT remaining after disturbance
!     C4PhotoShootNonstC_brch=total C4 nonstructural C in branch
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
        +C4PhotoShootNonstC_brch(NB,NZ)+StalkRsrvElms_brch(ielmc,NB,NZ)) &
        +ElmAllocmat4Litr(ielmc,ifoliar,M,NZ)*(LeafStrutElms_brch(ielmc,NB,NZ)*FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr) &
        +CanopyNodulStrutElms_brch(ielmc,NB,NZ)) &
        +ElmAllocmat4Litr(ielmc,inonfoliar,M,NZ)*(PetoleStrutElms_brch(ielmc,NB,NZ)*FracShootPetolElmAlloc2Litr(ielmc,k_fine_litr) &
        +HuskStrutElms_brch(ielmc,NB,NZ)+EarStrutElms_brch(ielmc,NB,NZ)))

      DO NE=2,NumPlantChemElms
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
          *(ElmAllocmat4Litr(NE,inonstruct,M,NZ)*(CanopyNonstElms_brch(NE,NB,NZ)+CanopyNodulNonstElms_brch(NE,NB,NZ)&
          +StalkRsrvElms_brch(NE,NB,NZ)) &
          +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_fine_litr) &
          +CanopyNodulStrutElms_brch(NE,NB,NZ)) &
          +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*(PetoleStrutElms_brch(NE,NB,NZ)*FracShootPetolElmAlloc2Litr(NE,k_fine_litr) &
          +HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ)))
      ENDDO
    ENDDO D6380

    DO M=1,jsken
      DO NE=1,NumPlantChemElms
        LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+XHVST1 &
          *ElmAllocmat4Litr(NE,icwood,M,NZ)*(LeafStrutElms_brch(NE,NB,NZ)*FracShootLeafElmAlloc2Litr(NE,k_woody_litr) &
          +PetoleStrutElms_brch(NE,NB,NZ)*FracShootPetolElmAlloc2Litr(NE,k_woody_litr))

        IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+XHVST1 &
            *ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
        ELSE
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
            *ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
        ENDIF
        LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+XHVST1 &
          *ElmAllocmat4Litr(NE,icwood,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_woody_litr)

        LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XHVST1 &
          *ElmAllocmat4Litr(NE,istalk,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)*FracWoodStalkElmAlloc2Litr(NE,k_fine_litr)
      ENDDO
    ENDDO
  !
  !     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
  !

    C4PhotoShootNonstC_brch(NB,NZ)     = C4PhotoShootNonstC_brch(NB,NZ)*XHVST
    SapwoodBiomassC_brch(NB,NZ) = SapwoodBiomassC_brch(NB,NZ)*XHVST
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

    PotentialSeedSites_brch(NB,NZ) = PotentialSeedSites_brch(NB,NZ)*XHVST
    SeedNumSet_brch(NB,NZ)         = SeedNumSet_brch(NB,NZ)*XHVST
    GrainSeedBiomCMean_brch(NB,NZ) = GrainSeedBiomCMean_brch(NB,NZ)*XHVST
    LeafAreaLive_brch(NB,NZ)       = LeafAreaLive_brch(NB,NZ)*XHVST
    LeafPetolBiomassC_brch(NB,NZ)  = AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetoleStrutElms_brch(ielmc,NB,NZ))
    CanopyLeafShethC_pft(NZ)       = CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(NB,NZ)
    CanopySapwoodC_pft(NZ)           = CanopySapwoodC_pft(NZ)+SapwoodBiomassC_brch(NB,NZ)
    D8970: DO K=0,MaxNodesPerBranch1
      IF(K.NE.0)THEN
        CPOOL3_node(K,NB,NZ)                = CPOOL3_node(K,NB,NZ)*XHVST
        CPOOL4_node(K,NB,NZ)                = CPOOL4_node(K,NB,NZ)*XHVST
        CMassCO2BundleSheath_node(K,NB,NZ)  = CMassCO2BundleSheath_node(K,NB,NZ)*XHVST
        CMassHCO3BundleSheath_node(K,NB,NZ) = CMassHCO3BundleSheath_node(K,NB,NZ)*XHVST
      ENDIF
      LeafArea_node(K,NB,NZ)=LeafArea_node(K,NB,NZ)*XHVST

      LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)*XHVST
  !     PetoleLensNode_brch(K,NB,NZ)=PetoleLensNode_brch(K,NB,NZ)*XHVST

      PetoleProteinCNode_brch(K,NB,NZ)=PetoleProteinCNode_brch(K,NB,NZ)*XHVST
  !     LiveInterNodeHight_brch(K,NB,NZ)=LiveInterNodeHight_brch(K,NB,NZ)*XHVST
  !     InternodeHeightDead_brch(K,NB,NZ)=InternodeHeightDead_brch(K,NB,NZ)*XHVST
      DO NE=1,NumPlantChemElms
        InternodeStrutElms_brch(NE,K,NB,NZ) = InternodeStrutElms_brch(NE,K,NB,NZ)*XHVST
        LeafElmntNode_brch(NE,K,NB,NZ)      = LeafElmntNode_brch(NE,K,NB,NZ)*XHVST
        PetioleElmntNode_brch(NE,K,NB,NZ)   = PetioleElmntNode_brch(NE,K,NB,NZ)*XHVST
        DO L=1,NumCanopyLayers1
          LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)=LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)*XHVST
        ENDDO
      ENDDO
      D8965: DO L=1,NumCanopyLayers1
        CanopyLeafArea_lnode(L,K,NB,NZ)=CanopyLeafArea_lnode(L,K,NB,NZ)*XHVST
      ENDDO D8965
    ENDDO D8970
  ENDIF
  ENDDO D8975
!
!     PSICanopy_pft=canopy water potential
!     CanopyBiomWater_pft=water volume in canopy
!     QH2OLoss_lnds,H2OLoss_CumYr_col=accumulated water loss for water balance calculation
!
  VOLWPX=CanopyBiomWater_pft(NZ)
  WVPLT=AZMAX1(CanopyLeafShethC_pft(NZ)+CanopySapwoodC_pft(NZ))

  FDM=get_FDM(PSICanopy_pft(NZ))
!        APSILT=ABS(PSICanopy_pft(NZ))
!        FDM=0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)

  CanopyBiomWater_pft(NZ)=ppmc*WVPLT/FDM
  QH2OLoss_lnds=QH2OLoss_lnds+VOLWPX-CanopyBiomWater_pft(NZ)
  H2OLoss_CumYr_col=H2OLoss_CumYr_col+VOLWPX-CanopyBiomWater_pft(NZ)
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
  call PrintInfo('end '//subname)
  end associate
  END subroutine RemoveShootByTillage

!----------------------------------------------------------------------------------------------------
  subroutine RemoveBiomByTillage(I,J,NZ)
  !
  !Description
  !REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
  implicit none
  integer , intent(in) :: I,J,NZ
  integer :: L,K,M,N,NR,NE,NB
  real(r8) :: XHVST
  REAL(R8) :: APSILT

!     begin_execution
  associate(                                                           &
    iSoilDisturbType_col      => plt_distb%iSoilDisturbType_col       ,& !input  :soil disturbance type, [-]
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft  ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    iDayPlanting_pft          => plt_distb%iDayPlanting_pft           ,& !input  :day of planting,[-]
    iYearCurrent              => plt_site%iYearCurrent                ,& !input  :current year,[-]
    iYearPlanting_pft         => plt_distb%iYearPlanting_pft          ,& !input  :year of planting,[-]
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft      ,& !input  :plant growth type (vascular, non-vascular),[-]
    XCORP                     => plt_distb%XCORP                      ,& !input  :factor for surface litter incorporation and soil mixing,[-]
    SolarNoonHour_col         => plt_site%SolarNoonHour_col           ,& !input  :time of solar noon, [h]
    PPX_pft                   => plt_site%PPX_pft                     ,& !inoput :plant population, [plants m-2]
    PlantPopulation_pft       => plt_site%PlantPopulation_pft          & !inoput :plant population, [d-2]
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
        XHVST                   = XCORP
        PPX_pft(NZ)             = PPX_pft(NZ)*XHVST
        PlantPopulation_pft(NZ) = PlantPopulation_pft(NZ)*XHVST
        
        call RemoveShootByTillage(I,J,NZ,XHVST)

!     LitrFall FROM ROOTS DURING TILLAGE
        call RemoveRootsByTillage(I,J,NZ,XHVST)

      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine RemoveBiomByTillage

!----------------------------------------------------------------------------------------------------
  subroutine RemoveRootsByTillage(I,J,NZ,XHVST)
  implicit none        
  integer, intent(in) :: I,J,NZ
  real(r8),intent(in) :: XHVST

  character(len=*), parameter :: subname='RemoveRootsByTillage'
  real(r8) :: XHVST1
  integer :: L,N,M,NE,NR,idg

  associate(                                                             &
    MaxNumRootLays             => plt_site%MaxNumRootLays               ,& !input  :maximum root layer number,[-]
    Myco_pft                   => plt_morph%Myco_pft                    ,& !input  :mycorrhizal type (no or yes),[-]
    ElmAllocmat4Litr           => plt_soilchem%ElmAllocmat4Litr         ,& !input  :litter kinetic fraction, [-]
    FracWoodStalkElmAlloc2Litr => plt_allom%FracWoodStalkElmAlloc2Litr  ,& !input  :woody element allocation,[-]
    FracRootElmAlloc2Litr      => plt_allom%FracRootElmAlloc2Litr       ,& !input  :C woody fraction in root,[-]
    NU                         => plt_site%NU                           ,& !input  :current soil surface layer number, [-]
    k_fine_litr                => pltpar%k_fine_litr                    ,& !input  :fine litter complex id
    k_woody_litr               => pltpar%k_woody_litr                   ,& !input  :woody litter complex id
    inonstruct                 => pltpar%inonstruct                     ,& !input  :group id of plant nonstructural litter
    iroot                      => pltpar%iroot                          ,& !input  :group id of plant root litter
    icwood                     => pltpar%icwood                         ,& !input  :group id of coarse woody litter
    iPlantNfixType_pft         => plt_morph%iPlantNfixType_pft          ,& !input  :N2 fixation type,[-]
    NGTopRootLayer_pft         => plt_morph%NGTopRootLayer_pft          ,& !input  :soil layer at planting depth, [-]
    NumRootAxes_pft            => plt_morph%NumRootAxes_pft             ,& !input  :root primary axis number,[-]
    trcg_rootml_pvr            => plt_rbgc%trcg_rootml_pvr              ,& !inoput :root gas content, [g d-2]
    trcs_rootml_pvr            => plt_rbgc%trcs_rootml_pvr              ,& !inoput :root aqueous content, [g d-2]
    RootMycoNonstElms_rpvr     => plt_biom%RootMycoNonstElms_rpvr       ,& !inoput :root layer nonstructural element, [g d-2]
    RootProteinC_pvr           => plt_biom%RootProteinC_pvr             ,& !inoput :root layer protein C, [gC d-2]
    PopuRootMycoC_pvr          => plt_biom% PopuRootMycoC_pvr           ,& !inoput :root layer C, [gC d-2]
    RootMycoActiveBiomC_pvr    => plt_biom%RootMycoActiveBiomC_pvr      ,& !inoput :root layer structural C, [gC d-2]
    RootMyco1stStrutElms_rpvr  => plt_biom%RootMyco1stStrutElms_rpvr    ,& !inoput :root layer element primary axes, [g d-2]
    RootMyco1stElm_raxs        => plt_biom%RootMyco1stElm_raxs          ,& !inoput :root C primary axes, [g d-2]
    SeasonalNonstElms_pft      => plt_biom%SeasonalNonstElms_pft        ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    RootMyco2ndStrutElms_rpvr  => plt_biom%RootMyco2ndStrutElms_rpvr    ,& !inoput :root layer element secondary axes, [g d-2]
    RootNodulNonstElms_rpvr    => plt_biom%RootNodulNonstElms_rpvr      ,& !inoput :root layer nonstructural element, [g d-2]
    RootNodulStrutElms_rpvr    => plt_biom%RootNodulStrutElms_rpvr      ,& !inoput :root layer nodule element, [g d-2]
    RootGasLossDisturb_pft     => plt_bgcr%RootGasLossDisturb_pft       ,& !inoput :gaseous flux fron root disturbance, [g d-2 h-1]
    LitrfalStrutElms_pvr       => plt_bgcr%LitrfalStrutElms_pvr         ,& !inoput :plant LitrFall element, [g d-2 h-1]
    RootCO2Autor_pvr           => plt_rbgc%RootCO2Autor_pvr             ,& !inoput :root respiration constrained by O2, [g d-2 h-1]
    RootRespPotent_pvr         => plt_rbgc%RootRespPotent_pvr           ,& !inoput :root respiration unconstrained by O2, [g d-2 h-1]
    RootCO2EmisPot_pvr         => plt_rbgc%RootCO2EmisPot_pvr           ,& !inoput :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    Root1stLen_rpvr            => plt_morph%Root1stLen_rpvr             ,& !inoput :root layer length primary axes, [m d-2]
    RootVH2O_pvr               => plt_morph%RootVH2O_pvr                ,& !inoput :root layer volume water, [m2 d-2]
    RootAreaPerPlant_pvr       => plt_morph%RootAreaPerPlant_pvr        ,& !inoput :root layer area per plant, [m p-1]
    RootPoreVol_rpvr            => plt_morph%RootPoreVol_rpvr             ,& !inoput :root layer volume air, [m2 d-2]
    RootLenDensPerPlant_pvr    => plt_morph%RootLenDensPerPlant_pvr     ,& !inoput :root layer length density, [m m-3]
    Root1stXNumL_rpvr           => plt_morph%Root1stXNumL_rpvr            ,& !inoput :root layer number primary axes, [d-2]
    RootLenPerPlant_pvr        => plt_morph%RootLenPerPlant_pvr         ,& !inoput :root layer length per plant, [m p-1]
    Root2ndXNum_rpvr           => plt_morph%Root2ndXNum_rpvr            ,& !inoput :root layer number secondary axes, [d-2]
    Root2ndLen_rpvr            => plt_morph%Root2ndLen_rpvr             ,& !inoput :root layer length secondary axes, [m d-2]
    Root2ndXNumL_rpvr            => plt_morph%Root2ndXNumL_rpvr              & !inoput :root layer number axes, [d-2]
  )
!
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from disturbance
!     XHVST=fraction of PFT remaining after disturbance
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!
  call PrintInfo('beg '//subname)
  XHVST1=1._r8-XHVST
  
  DO NR=1,NumRootAxes_pft(NZ)
    DO N=1,Myco_pft(NZ)
      DO NE=1,NumPlantChemElms
        RootMyco1stElm_raxs(NE,N,NR,NZ)=RootMyco1stElm_raxs(NE,N,NR,NZ)*XHVST
      ENDDO
    ENDDO
  ENDDO

  D8980: DO L=NU,MaxNumRootLays
    D8985: DO N=1,Myco_pft(NZ)
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

      DO idg=idg_beg,idg_NH3
        RootGasLossDisturb_pft(idg,NZ)=RootGasLossDisturb_pft(idg,NZ)-XHVST1 &
          *(trcg_rootml_pvr(idg,N,L,NZ)+trcs_rootml_pvr(idg,N,L,NZ))

        trcg_rootml_pvr(idg,N,L,NZ) = XHVST*trcg_rootml_pvr(idg,N,L,NZ)
        trcs_rootml_pvr(idg,N,L,NZ) = XHVST*trcs_rootml_pvr(idg,N,L,NZ)
      ENDDO
!
!     ROOT STATE VARIABLES REMAINING AFTER TILLAGE
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     Root1stLen_rpvr,Root2ndLen_rpvr=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RootMycoActiveBiomC_pvr, PopuRootMycoC_pvr=active,actual root C mass
!     RootProteinC_pvr=root protein C mass
!     RTN1,Root2ndXNumL_rpvr=number of primary,secondary root axes
!     RootLenDensPerPlant_pvr,RootLenPerPlant_pvr=root length density,root length per plant
!     RootVH2O_pvr,RootPoreVol_rpvr=root or myco aqueous,gaseous volume
!     RootAreaPerPlant_pvr=root surface area per plant
!     RootRespPotent_pvr,RootCO2EmisPot_pvr,RootCO2Autor_pvr unlimited by O2,nonstructural C
!
      D8960: DO NR=1,NumRootAxes_pft(NZ)
        DO NE=1,NumPlantChemElms
          RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ) = RootMyco1stStrutElms_rpvr(NE,N,L,NR,NZ)*XHVST
          RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ) = RootMyco2ndStrutElms_rpvr(NE,N,L,NR,NZ)*XHVST
        ENDDO
        Root1stLen_rpvr(N,L,NR,NZ)  = Root1stLen_rpvr(N,L,NR,NZ)*XHVST
        Root2ndLen_rpvr(N,L,NR,NZ)  = Root2ndLen_rpvr(N,L,NR,NZ)*XHVST
        Root2ndXNum_rpvr(N,L,NR,NZ) = Root2ndXNum_rpvr(N,L,NR,NZ)*XHVST
      ENDDO D8960

      DO NE=1,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)*XHVST
      ENDDO
      RootMycoActiveBiomC_pvr(N,L,NZ) = RootMycoActiveBiomC_pvr(N,L,NZ)*XHVST
      PopuRootMycoC_pvr(N,L,NZ)       = PopuRootMycoC_pvr(N,L,NZ)*XHVST
      RootProteinC_pvr(N,L,NZ)        = RootProteinC_pvr(N,L,NZ)*XHVST
      Root1stXNumL_rpvr(N,L,NZ)        = Root1stXNumL_rpvr(N,L,NZ)*XHVST
      Root2ndXNumL_rpvr(N,L,NZ)         = Root2ndXNumL_rpvr(N,L,NZ)*XHVST
      RootLenPerPlant_pvr(N,L,NZ)     = RootLenPerPlant_pvr(N,L,NZ)*XHVST
      RootLenDensPerPlant_pvr(N,L,NZ) = RootLenDensPerPlant_pvr(N,L,NZ)*XHVST
      RootPoreVol_rpvr(N,L,NZ)         = RootPoreVol_rpvr(N,L,NZ)*XHVST
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
              XHVST1*(ElmAllocmat4Litr(NE,iroot,M,NZ)*RootNodulStrutElms_rpvr(NE,L,NZ) &
              +ElmAllocmat4Litr(NE,inonstruct,M,NZ)*RootNodulNonstElms_rpvr(NE,L,NZ))
          ENDDO D6395
          RootNodulStrutElms_rpvr(NE,L,NZ) = RootNodulStrutElms_rpvr(NE,L,NZ)*XHVST
          RootNodulNonstElms_rpvr(NE,L,NZ) = RootNodulNonstElms_rpvr(NE,L,NZ)*XHVST
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
        *FracWoodStalkElmAlloc2Litr(NE,k_woody_litr)

      LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ)= &
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,NGTopRootLayer_pft(NZ),NZ) &
        +(XHVST1*ElmAllocmat4Litr(NE,inonstruct,M,NZ)*SeasonalNonstElms_pft(NE,NZ)) &
        *FracWoodStalkElmAlloc2Litr(NE,k_fine_litr)
    ENDDO D6400
    SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)*XHVST
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine RemoveRootsByTillage        
  ![tail]
end module PlantDisturbByTillageMod