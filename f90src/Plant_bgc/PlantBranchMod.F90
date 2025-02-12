module PlantBranchMod

! Description:
! module for plant biological transformations
  use minimathmod, only : isclose,safe_adb,AZMAX1
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod , only : etimer 
  use DebugToolMod
  use EcosimConst
  use GrosubPars
  use PlantAPIData
  use PhotoSynsMod
  use PlantMathFuncMod
  use NoduleBGCMod  
  use LitrFallMod
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: GrowOneBranch

  integer, parameter :: ibrch_leaf=1
  integer, parameter :: ibrch_petole=2
  integer, parameter :: ibrch_stalk=3
  integer, parameter :: ibrch_resrv=4
  integer, parameter :: ibrch_husk=5
  integer, parameter :: ibrch_ear=6
  integer, parameter :: ibrch_grain=7

  contains
!------------------------------------------------------------------------------------------

  subroutine GrowOneBranch(I,J,NB,NZ,TFN6_vr,CanopyHeight_copy,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,&
    CPRTW,TFN5,WaterStress4Groth,Stomata_Stress,WFNS,CanTurgPSIFun4Expans,PTRT,CanopyN2Fix_pft,BegRemoblize)
  implicit none
  integer, intent(in)  :: I,J,NB,NZ
  REAL(R8), INTENT(IN) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: CanopyHeight_copy(JP1)
  real(r8), intent(in) :: CNLFW,CPLFW
  real(r8), intent(in) :: CNSHW,CPSHW
  real(r8), intent(in) :: CNRTW,CPRTW
  real(r8), intent(in) :: TFN5,WaterStress4Groth
  real(r8), intent(in) :: Stomata_Stress
  real(r8), intent(in) :: WFNS,CanTurgPSIFun4Expans
  real(r8), intent(inout) :: CanopyN2Fix_pft(JP1)
  real(r8), intent(out) :: PTRT
  integer, intent(out) :: BegRemoblize
  real(r8) :: DMSHD
  integer  :: MinNodeNum  
  REAL(R8) :: CH2O3(pltpar%MaxNodesPerBranch1),CH2O4(pltpar%MaxNodesPerBranch1)
  integer  :: LRemob_brch
  REAL(R8) :: PART(pltpar%NumOfPlantMorphUnits)
  real(r8) :: ALLOCL
  REAL(R8) :: CNPG
  real(r8) :: CCE
  real(r8) :: CCC,CNC,CPC
  real(r8) :: EtoliationCoeff
  real(r8) :: Growth_brch(NumPlantChemElms,pltpar%NumOfPlantMorphUnits)
  real(r8) :: RCO2NonstC_brch
  REAL(R8) :: RCO2Maint_brch
  real(r8) :: RMxess_brch
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: SSL
  real(r8) :: CNLFM,CPLFM
  real(r8) :: CNLFX,CPLFX,CNSHX,CPSHX
  real(r8) :: NonstC4Groth_brch
  real(r8) :: RCO2NonstC4Nassim_brch
  real(r8) :: ShootStructE_brch(1:NumPlantChemElms)
  real(r8) :: DMLFB
  real(r8) :: DMSHB
  real(r8) :: CNLFB
  real(r8) :: CPLFB  
  real(r8) :: CNSHB,CPSHB

! begin_execution
  associate(                                                                &
    LeafPetolBiomassC_brch       =>  plt_biom%LeafPetolBiomassC_brch      , &
    LeafProteinCNode_brch        =>  plt_biom%LeafProteinCNode_brch       , &
    LeafPetoNonstElmConc_brch    =>  plt_biom%LeafPetoNonstElmConc_brch   , &
    LeafStrutElms_brch           =>  plt_biom%LeafStrutElms_brch          , &      
    PetoleStrutElms_brch         =>  plt_biom%PetoleStrutElms_brch        , &    
    iPlantPhotosynthesisType     =>  plt_photo%iPlantPhotosynthesisType   , &
    iPlantBranchState_brch       =>  plt_pheno%iPlantBranchState_brch     , &
    iPlantRootProfile_pft        =>  plt_pheno%iPlantRootProfile_pft      , &
    iPlantCalendar_brch          =>  plt_pheno%iPlantCalendar_brch        , &
    iPlantTurnoverPattern_pft    =>  plt_pheno%iPlantTurnoverPattern_pft  , &
    iPlantPhenolPattern_pft      =>  plt_pheno%iPlantPhenolPattern_pft    , &
    SineSunInclAnglNxtHour_col   =>  plt_rad%SineSunInclAnglNxtHour_col   , &
    PlantPopulation_pft          =>  plt_site%PlantPopulation_pft         , &
    ZERO                         =>  plt_site%ZERO                        , &
    canopy_growth_pft            =>  plt_rbgc%canopy_growth_pft           , &
    MainBranchNum_pft            =>  plt_morph%MainBranchNum_pft            &
  )

  IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN

    call CalcPartitionCoeff(I,J,NB,NZ,PART,PTRT,LRemob_brch,BegRemoblize)

    call UpdateBranchAllometry(I,J,NZ,NB,PART,CNLFW,CNRTW,CNSHW,CPLFW,&
      CPRTW,CPSHW,ShootStructE_brch,DMSHD,CNLFM,CPLFM,&
      CNSHX,CPSHX,CNLFX,CPLFX,DMLFB,DMSHB,CNLFB,CPLFB,CNSHB,CPSHB)
!
!   GROSS PRIMARY PRODUCTIVITY
!
    call UpdatePhotosynthates(I,J,NB,NZ,TFN6_vr,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX &
      ,CNLFX,CPLFX,ShootStructE_brch(ielmn),TFN5,WaterStress4Groth,Stomata_Stress,CanTurgPSIFun4Expans &
      ,CH2O3,CH2O4,CNPG,RCO2NonstC_brch,RCO2Maint_brch,RMxess_brch,NonstC4Groth_brch,RCO2NonstC4Nassim_brch)
!
!   TRANSFER OF C4 FIXATION PRODUCTS FROM NON-STRUCTURAL POOLS
!   IN MESOPHYLL TO THOSE IN BUNDLE SHEATH, DECARBOXYLATION
!   OF C4 FIXATION PRODUCTS IN BUNDLE SHEATH, LEAKAGE OF DECARBOXYLATION
!   PRODUCTS BACK TO MESOPHYLL IN C4 PLANTS
!
!   iPlantPhotosynthesisType=photosynthesis type:3=C3,4=C4
!
    IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo)THEN
      call C4PhotoProductTransfer(I,J,NZ,NB,CH2O3,CH2O4)
    ENDIF
    canopy_growth_pft(NZ) = canopy_growth_pft(NZ)+NonstC4Groth_brch

    call BranchBiomAllocate(I,J,NB,NZ,PART,NonstC4Groth_brch,&
      DMLFB,DMSHB,CNLFB,CPLFB,CNSHB,CPSHB,ZPLFD,CNPG,Growth_brch,EtoliationCoeff,MinNodeNum)

    CALL GrowLeavesOnBranch(I,J,NZ,NB,MinNodeNum,Growth_brch(:,ibrch_leaf),EtoliationCoeff,WFNS,ALLOCL)
!
!     DISTRIBUTE SHEATH OR PETIOLE GROWTH AMONG CURRENTLY GROWING NODES
!
    CALL GrowPetioleOnBranch(NZ,NB,MinNodeNum,Growth_brch(:,ibrch_petole),EtoliationCoeff,WFNS,ALLOCL)

!
!   DISTRIBUTE STALK GROWTH AMONG CURRENTLY GROWING NODES
!
    call GrowStalkOnBranch(NZ,NB,Growth_brch(:,ibrch_stalk),EtoliationCoeff)

!
    !   RECOVERY OF REMOBILIZABLE N,P DURING REMOBILIZATION DEPENDS
    !   ON SHOOT NON-STRUCTURAL C:N:P
    !
    !   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
    !   CNKI,CPKI=nonstruct N,P inhibition constant on growth (g N,P g-1 C)
    !   RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !   RCCZ,RCCY=min,max fractions for shoot C recycling
    !   RCCX,RCCQ=max fractions for shoot N,P recycling
    !   iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
    IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0.AND.LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZERO)THEN
      CCC=AZMAX1(AMIN1(1.0_r8 &
        ,LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ) &
        +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI) &
        ,LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ) &
        +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI)))
      CNC=AZMAX1(AMIN1(1.0_r8 &
        ,LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ) &
        +LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/CNKI)))
      CPC=AZMAX1(AMIN1(1.0_r8 &
        ,LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ) &
        +LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/CPKI)))
    ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
    ENDIF
    RCCC = RCCZ(iPlantRootProfile_pft(NZ))+CCC*RCCY(iPlantRootProfile_pft(NZ))
    RCCN = CNC*RCCX(iPlantRootProfile_pft(NZ))
    RCCP = CPC*RCCQ(iPlantRootProfile_pft(NZ))
!
!       WITHDRAW REMOBILIZABLE C,N,P FROM LOWEST NODE AFTER
!       MAXIMUM NODE NUMBER OF 25 IS REACHED
!
!       doSenescence_brch=PFT branch senescence flag
!       KHiestGroLeafNode_brch=integer of most recent leaf number
!       fTCanopyGroth_pft=temperature function for canopy growth
!       XRLA=rate of leaf appearance at 25 oC (h-1)
!       FSNC=fraction of lowest leaf to be remobilized
!
    call SenescenceBranch(NZ,NB,RCCC,RCCN,RCCP)
!
    call RemobilizeBranch(NZ,NB,BegRemoblize,LRemob_brch,RCCC,RCCN,RCCP,RMxess_brch)
!
!
!   DEATH IF MAIN STALK OF TREE DIES
!
!   iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+petiole,2=none,3=between 1,2
!   iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!   iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!   KMinGroingLeafNodeNum,KHiestGroLeafNode_brch=integer of lowest,highest leaf number currently growing

!
    IF(iPlantTurnoverPattern_pft(NZ).NE.0 .AND. is_plant_treelike(iPlantRootProfile_pft(NZ)) &
      .AND. iPlantBranchState_brch(MainBranchNum_pft(NZ),NZ).EQ.iDead)then
      iPlantBranchState_brch(NB,NZ)=iDead
    endif

!     REMOBILIZE EXCESS LEAF STRUCTURAL N,P
    call WithdrawBranchLeaves(I,J,NB,NZ,CNLFB,CPLFB)

    call AllocateLeaf2CanopyLayers(I,J,NB,NZ,CanopyHeight_copy)

!
    !     ALLOCATE LEAF AREA TO INCLINATION CLASSES ACCORDING TO
    !     DISTRIBUTION ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
    !
    !     SineSunInclAngle_col=sine of solar angle
    !     LeafAreaZsec_brch=leaf node surface area in canopy layer
    !     LeafNodeArea_brch,CanopyLeafArea_lpft=leaf node surface area in canopy layer
    !     ZC,DPTHS=canopy,snowpack height
    !     CLASS=leaf inclination class
    !
    IF(SineSunInclAnglNxtHour_col.GT.0.0_r8)THEN
      call LeafClassAllocation(NB,NZ)
    ENDIF

    call GrainFilling(I,J,NB,NZ,Growth_brch,Growth_brch(ielmc,ibrch_stalk))
!
    call ResetBranchPhenology(I,J,NB,NZ)
!   
    call BranchElmntTransfer(I,J,NB,NZ,BegRemoblize,WaterStress4Groth,CanTurgPSIFun4Expans)

!   CANOPY N2 FIXATION (CYANOBACTERIA)
!
    call CanopyNoduleBiochemistry(I,J,NZ,NB,TFN5,WaterStress4Groth,CanopyN2Fix_pft)

    LeafPetolBiomassC_brch(NB,NZ)=AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetoleStrutElms_brch(ielmc,NB,NZ))

  ENDIF
  end associate
  end subroutine GrowOneBranch
!------------------------------------------------------------------------------------------
  subroutine WithdrawBranchLeaves(I,J,NB,NZ,CNLFB,CPLFB)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in):: CNLFB,CPLFB
  integer :: KK,NE,K
  integer :: KMinGroingLeafNodeNum
  real(r8) :: ZPOOLD,XFRN1,XFRP1
  real(r8) :: CPOOLT
  real(r8) :: XFRE(1:NumPlantChemElms)
  real(r8) :: PPOOLD
  associate(                                                    &
    ZERO4Groth_pft         => plt_biom%ZERO4Groth_pft,          &
    KHiestGroLeafNode_brch => plt_pheno%KHiestGroLeafNode_brch, &
    LeafElmntNode_brch     => plt_biom%LeafElmntNode_brch,      &
    LeafStrutElms_brch     => plt_biom%LeafStrutElms_brch,      &
    LeafProteinCNode_brch  => plt_biom%LeafProteinCNode_brch,   &
    rCNNonstRemob_pft      => plt_allom%rCNNonstRemob_pft,      &
    rCPNonstRemob_pft      => plt_allom%rCPNonstRemob_pft,      &
    CanopyNonstElms_brch   => plt_biom%CanopyNonstElms_brch     &
  )

  KMinGroingLeafNodeNum=MAX(0,KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+1)
  D495: DO KK=KMinGroingLeafNodeNum,KHiestGroLeafNode_brch(NB,NZ)
    K=pMOD(KK,MaxNodesPerBranch1)
    IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.0.0_r8)THEN
      CPOOLT=LeafElmntNode_brch(ielmc,K,NB,NZ)+CanopyNonstElms_brch(ielmc,NB,NZ)
      IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
        ZPOOLD=LeafElmntNode_brch(ielmn,K,NB,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ) &
          -CanopyNonstElms_brch(ielmn,NB,NZ)*LeafElmntNode_brch(ielmc,K,NB,NZ)
        PPOOLD=LeafElmntNode_brch(ielmp,K,NB,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ) &
          -CanopyNonstElms_brch(ielmp,NB,NZ)*LeafElmntNode_brch(ielmc,K,NB,NZ)
          
        XFRN1=AZMAX1(AMIN1(1.0E-03_r8*ZPOOLD/CPOOLT,LeafElmntNode_brch(ielmn,K,NB,NZ) &
          -ZPLFM*CNLFB*LeafElmntNode_brch(ielmc,K,NB,NZ)))            
        XFRP1=AZMAX1(AMIN1(1.0E-03_r8*PPOOLD/CPOOLT,LeafElmntNode_brch(ielmp,K,NB,NZ) &
          -ZPLFM*CPLFB*LeafElmntNode_brch(ielmc,K,NB,NZ)))
          
        XFRE(ielmn)=AMAX1(XFRN1,10.0_r8*XFRP1)
        XFRE(ielmp)=AMAX1(XFRP1,0.10_r8*XFRN1)
        DO NE=2,NumPlantChemElms
          LeafElmntNode_brch(NE,K,NB,NZ) = LeafElmntNode_brch(NE,K,NB,NZ)-XFRE(NE)
          LeafStrutElms_brch(NE,NB,NZ)   = LeafStrutElms_brch(NE,NB,NZ)-XFRE(NE)
          CanopyNonstElms_brch(NE,NB,NZ) = CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE)
        ENDDO
        LeafProteinCNode_brch(K,NB,NZ)=AZMAX1(LeafProteinCNode_brch(K,NB,NZ)-&
          AMAX1(XFRE(ielmn)*rCNNonstRemob_pft(NZ),XFRE(ielmp)*rCPNonstRemob_pft(NZ)))
      ENDIF
    ENDIF
  ENDDO D495
  end associate
  end subroutine WithdrawBranchLeaves
!------------------------------------------------------------------------------------------

  subroutine BranchBiomAllocate(I,J,NB,NZ,PART,NonstC4Groth_brch,&
    DMLFB,DMSHB,CNLFB,CPLFB,CNSHB,CPSHB,ZPLFD,CNPG,Growth_brch,EtoliationCoeff,MinNodeNum)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: PART(pltpar%NumOfPlantMorphUnits)
  real(r8), intent(in) :: NonstC4Groth_brch                  !non-structural C used for branch growth
  real(r8), intent(in) :: DMLFB,DMSHB,CNLFB
  real(r8), intent(in) :: CPLFB,CNSHB,CPSHB,ZPLFD,CNPG
  real(r8), intent(out) :: Growth_brch(NumPlantChemElms,pltpar%NumOfPlantMorphUnits)
  integer , intent(out) :: MinNodeNum
  real(r8), intent(out) :: EtoliationCoeff

  character(len=*), parameter :: subname='BranchBiomAllocate'
  real(r8) :: CCE
  integer  :: NE
  
  associate(                                                         &
    StalkStrutElms_brch       => plt_biom%StalkStrutElms_brch,       &
    StalkRsrvElms_brch        => plt_biom%StalkRsrvElms_brch,        &
    HuskStrutElms_brch        => plt_biom%HuskStrutElms_brch,        &
    EarStrutElms_brch         => plt_biom%EarStrutElms_brch,         &
    PetoleStrutElms_brch      => plt_biom%PetoleStrutElms_brch,      &
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch,        &
    LeafPetoNonstElmConc_brch => plt_biom%LeafPetoNonstElmConc_brch, &
    LeafBiomGrowthYield       => plt_allom%LeafBiomGrowthYield,      &
    PetioleBiomGrowthYield    => plt_allom%PetioleBiomGrowthYield,   &
    EarBiomGrowthYield        => plt_allom%EarBiomGrowthYield,       &
    rNCEar_pft                => plt_allom%rNCEar_pft,               &
    rPCReserve_pft            => plt_allom%rPCReserve_pft,           &
    rPCHusk_pft               => plt_allom%rPCHusk_pft,              &
    rPCStalk_pft              => plt_allom%rPCStalk_pft,             &
    rNCStalk_pft              => plt_allom%rNCStalk_pft,             &
    StalkBiomGrowthYield      => plt_allom%StalkBiomGrowthYield,     &
    HuskBiomGrowthYield       => plt_allom%HuskBiomGrowthYield,      &
    rNCHusk_pft               => plt_allom%rNCHusk_pft,              &
    rPCEar_pft                => plt_allom%rPCEar_pft,               &
    ReserveBiomGrowthYield    => plt_allom%ReserveBiomGrowthYield,   &
    GrainBiomGrowthYield      => plt_allom%GrainBiomGrowthYield,     &
    rCNNonstRemob_pft         => plt_allom%rCNNonstRemob_pft,        &
    rCPNonstRemob_pft         => plt_allom%rCPNonstRemob_pft,        &
    SeedDepth_pft             => plt_morph%SeedDepth_pft,            &
    HypoctoHeight_pft         => plt_morph%HypoctoHeight_pft,        &
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft,        &
    rNCReserve_pft            => plt_allom%rNCReserve_pft            &
  )
!
!   C,N,P GROWTH OF LEAF, SHEATH OR PETIOLE, STALK,
!   STALK RESERVES, REPRODUCTIVE ORGANS, GRAIN
!

  call PrintInfo('beg '//subname)
  Growth_brch(ielmc,ibrch_leaf)   =NonstC4Groth_brch*PART(ibrch_leaf)  *DMLFB       !leaf
  Growth_brch(ielmc,ibrch_petole) =NonstC4Groth_brch*PART(ibrch_petole)*DMSHB       !petiole
  Growth_brch(ielmc,ibrch_stalk)  =NonstC4Groth_brch*PART(ibrch_stalk) *StalkBiomGrowthYield(NZ)  !stalk
  Growth_brch(ielmc,ibrch_resrv)  =NonstC4Groth_brch*PART(ibrch_resrv) *ReserveBiomGrowthYield(NZ)  !reserve
  Growth_brch(ielmc,ibrch_husk)   =NonstC4Groth_brch*PART(ibrch_husk)  *HuskBiomGrowthYield(NZ)  !husk
  Growth_brch(ielmc,ibrch_ear)    =NonstC4Groth_brch*PART(ibrch_ear)   *EarBiomGrowthYield(NZ)  !ear
  Growth_brch(ielmc,ibrch_grain)  =NonstC4Groth_brch*PART(ibrch_grain) *GrainBiomGrowthYield(NZ)    !grain

  Growth_brch(ielmn,ibrch_leaf)   =Growth_brch(ielmc,ibrch_leaf)  *CNLFB*(ZPLFM+ZPLFD*CNPG)
  Growth_brch(ielmn,ibrch_petole) =Growth_brch(ielmc,ibrch_petole)*CNSHB
  Growth_brch(ielmn,ibrch_stalk)  =Growth_brch(ielmc,ibrch_stalk) *rNCStalk_pft(NZ)
  Growth_brch(ielmn,ibrch_resrv)  =Growth_brch(ielmc,ibrch_resrv) *rNCReserve_pft(NZ)
  Growth_brch(ielmn,ibrch_husk)   =Growth_brch(ielmc,ibrch_husk)  *rNCHusk_pft(NZ)
  Growth_brch(ielmn,ibrch_ear)    =Growth_brch(ielmc,ibrch_ear)   *rNCEar_pft(NZ)
  Growth_brch(ielmn,ibrch_grain)  =Growth_brch(ielmc,ibrch_grain) *rNCReserve_pft(NZ)

  Growth_brch(ielmp,ibrch_leaf)  =Growth_brch(ielmc,ibrch_leaf)  *CPLFB*(ZPLFM+ZPLFD*CNPG)
  Growth_brch(ielmp,ibrch_petole)=Growth_brch(ielmc,ibrch_petole)*CPSHB
  Growth_brch(ielmp,ibrch_stalk) =Growth_brch(ielmc,ibrch_stalk) *rPCStalk_pft(NZ)
  Growth_brch(ielmp,ibrch_resrv) =Growth_brch(ielmc,ibrch_resrv) *rPCReserve_pft(NZ)
  Growth_brch(ielmp,ibrch_husk)  =Growth_brch(ielmc,ibrch_husk)  *rPCHusk_pft(NZ)
  Growth_brch(ielmp,ibrch_ear)   =Growth_brch(ielmc,ibrch_ear)   *rPCEar_pft(NZ)
  Growth_brch(ielmp,ibrch_grain) =Growth_brch(ielmc,ibrch_grain) *rPCReserve_pft(NZ)
  
  DO NE=1,NumPlantChemElms
    LeafStrutElms_brch(NE,NB,NZ)   = LeafStrutElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_leaf)
    PetoleStrutElms_brch(NE,NB,NZ) = PetoleStrutElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_petole)
    StalkStrutElms_brch(NE,NB,NZ)  = StalkStrutElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_stalk)
    StalkRsrvElms_brch(NE,NB,NZ)   = StalkRsrvElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_resrv)
    HuskStrutElms_brch(NE,NB,NZ)   = HuskStrutElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_husk)
    EarStrutElms_brch(NE,NB,NZ)    = EarStrutElms_brch(NE,NB,NZ)+Growth_brch(NE,ibrch_ear)
  ENDDO

!
!   ETOLIATION
!
!   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!   CNKI,CPKI=nonstruct N,P inhibition constant on growth (g N,P g-1 C)
!   ETOL=coefficient for etoliation effects on expansion,extension
!
  CCE=AMIN1(safe_adb(LeafPetoNonstElmConc_brch(ielmn,NB,NZ) &
    ,LeafPetoNonstElmConc_brch(ielmn,NB,NZ)+LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI) &
    ,safe_adb(LeafPetoNonstElmConc_brch(ielmp,NB,NZ) &
    ,LeafPetoNonstElmConc_brch(ielmp,NB,NZ)+LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI))

  EtoliationCoeff = 1.0_r8+CCE
!
!   DISTRIBUTE LEAF GROWTH AMONG CURRENTLY GROWING NODES
!
  IF(NB.EQ.MainBranchNum_pft(NZ) .AND. HypoctoHeight_pft(NZ).LE.SeedDepth_pft(NZ))THEN
    MinNodeNum=0
  ELSE
    MinNodeNum=1
  ENDIF
  call PrintInfo('end '//subname)
  end associate  
  end subroutine BranchBiomAllocate

!------------------------------------------------------------------------------------------
  subroutine UpdateBranchAllometry(I,J,NZ,NB,PART,CNLFW,CNRTW,CNSHW,CPLFW,&
    CPRTW,CPSHW,ShootStructE_brch,DMSHD,CNLFM,CPLFM,&
    CNSHX,CPSHX,CNLFX,CPLFX,DMLFB,DMSHB,CNLFB,CPLFB,CNSHB,CPSHB)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: CNLFW,CNRTW,CNSHW,CPLFW,CPRTW,CPSHW
  real(r8), intent(in) :: PART(pltpar%NumOfPlantMorphUnits)
  real(r8), intent(out) :: ShootStructE_brch(NumPlantChemElms)    
  real(r8), intent(out) :: DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX
  real(r8), intent(out) :: DMLFB,DMSHB,CNLFB,CPLFB,CNSHB,CPSHB

  character(len=*), parameter :: subname='UpdateBranchAllometry'
  real(r8)   :: DMSHT  
  associate(                                                    &
    LeafBiomGrowthYield    => plt_allom%LeafBiomGrowthYield,    &
    PetioleBiomGrowthYield => plt_allom%PetioleBiomGrowthYield, &
    EarBiomGrowthYield     => plt_allom%EarBiomGrowthYield,     &
    rNCEar_pft             => plt_allom%rNCEar_pft,             &
    rPCReserve_pft         => plt_allom%rPCReserve_pft,         &
    rPCHusk_pft            => plt_allom%rPCHusk_pft,            &
    rPCStalk_pft           => plt_allom%rPCStalk_pft,           &
    rNCStalk_pft           => plt_allom%rNCStalk_pft,           &
    rNCReserve_pft         => plt_allom%rNCReserve_pft,         &
    StalkBiomGrowthYield   => plt_allom%StalkBiomGrowthYield,   &
    HuskBiomGrowthYield    => plt_allom%HuskBiomGrowthYield,    &
    rNCHusk_pft            => plt_allom%rNCHusk_pft,            &
    rPCEar_pft             => plt_allom%rPCEar_pft,             &
    ReserveBiomGrowthYield => plt_allom%ReserveBiomGrowthYield, &
    GrainBiomGrowthYield   => plt_allom%GrainBiomGrowthYield,   &
    iPlantCalendar_brch    => plt_pheno%iPlantCalendar_brch,    &
    PetoleStrutElms_brch   => plt_biom%PetoleStrutElms_brch,    &
    GrainStrutElms_brch    => plt_biom%GrainStrutElms_brch,     &
    EarStrutElms_brch      => plt_biom%EarStrutElms_brch,       &
    HuskStrutElms_brch     => plt_biom%HuskStrutElms_brch,      &
    LeafStrutElms_brch     => plt_biom%LeafStrutElms_brch,      &
    StalkLiveBiomassC_brch => plt_biom%StalkLiveBiomassC_brch,  &
    RootBiomGrosYld_pft    => plt_allom%RootBiomGrosYld_pft     &
  )

  call PrintInfo('beg '//subname)
!   SHOOT COEFFICIENTS FOR GROWTH RESPIRATION AND N,P CONTENTS
!   FROM GROWTH YIELDS ENTERED IN 'READQ', AND FROM PARTITIONING
!   COEFFICIENTS ABOVE
!
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   ZPLFD=1.0_r8-ZPLFM
!
    IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0)THEN
      DMLFB = LeafBiomGrowthYield(NZ)
      DMSHB = PetioleBiomGrowthYield(NZ)
      CNLFB = CNLFW
      CNSHB = CNSHW
      CPLFB = CPLFW
      CPSHB = CPSHW
    ELSE
      DMLFB = RootBiomGrosYld_pft(NZ)
      DMSHB = RootBiomGrosYld_pft(NZ)
      CNLFB = CNRTW
      CNSHB = CNRTW
      CPLFB = CPRTW
      CPSHB = CPRTW
    ENDIF
    !part 1(leaf), (2) petiole, (3) stalk, (4) reserve, (5) husk, (6) ear, (7) grain
  DMSHT = PART(ibrch_leaf)*DMLFB+PART(ibrch_petole)*DMSHB+PART(ibrch_stalk)*StalkBiomGrowthYield(NZ) &
    +PART(ibrch_resrv)*ReserveBiomGrowthYield(NZ)+PART(ibrch_husk)*HuskBiomGrowthYield(NZ) &
    +PART(ibrch_ear)*EarBiomGrowthYield(NZ)+PART(ibrch_grain)*GrainBiomGrowthYield(NZ)

  DMSHD = 1.0_r8-DMSHT
  CNLFM = PART(ibrch_leaf)*DMLFB*ZPLFM*CNLFB
  CPLFM = PART(ibrch_leaf)*DMLFB*ZPLFM*CPLFB
  CNLFX = PART(ibrch_leaf)*DMLFB*ZPLFD*CNLFB
  CPLFX = PART(ibrch_leaf)*DMLFB*ZPLFD*CPLFB
  CNSHX = PART(ibrch_petole)*DMSHB*CNSHB &
    +PART(ibrch_stalk)*StalkBiomGrowthYield(NZ)*rNCStalk_pft(NZ) &
    +PART(ibrch_resrv)*ReserveBiomGrowthYield(NZ)*rNCReserve_pft(NZ) &
    +PART(ibrch_husk)*HuskBiomGrowthYield(NZ)*rNCHusk_pft(NZ) &
    +PART(ibrch_ear)*EarBiomGrowthYield(NZ)*rNCEar_pft(NZ) &
    +PART(ibrch_grain)*GrainBiomGrowthYield(NZ)*rNCReserve_pft(NZ)
  CPSHX = PART(ibrch_petole)*DMSHB*CPSHB &
    +PART(ibrch_stalk)*StalkBiomGrowthYield(NZ)*rPCStalk_pft(NZ) &
    +PART(ibrch_resrv)*ReserveBiomGrowthYield(NZ)*rPCReserve_pft(NZ) &
    +PART(ibrch_husk)*HuskBiomGrowthYield(NZ)*rPCHusk_pft(NZ) &
    +PART(ibrch_ear)*EarBiomGrowthYield(NZ)*rPCEar_pft(NZ) &
    +PART(ibrch_grain)*GrainBiomGrowthYield(NZ)*rPCReserve_pft(NZ)
!
!   TOTAL SHOOT STRUCTURAL N MASS FOR MAINTENANCE RESPIRATION
!
!   ShootStructE_brch(ielmn)=shoot structural N mass
!   WTLFBN,WTSHBN,WTHSBN,WTEARN,WTFRBN=leaf,petiole,husk,ear,grain N mass
!   rNCStalk_pft,StalkLiveBiomassC_brch=stalk N:C,sapwood mass
!   iPlantCalendar_brch(10=date of physiological maturity
!
  ShootStructE_brch(ielmc) = AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetoleStrutElms_brch(ielmc,NB,NZ)+StalkLiveBiomassC_brch(NB,NZ))
  ShootStructE_brch(ielmn) = AZMAX1(LeafStrutElms_brch(ielmn,NB,NZ)+PetoleStrutElms_brch(ielmn,NB,NZ)+StalkLiveBiomassC_brch(NB,NZ))*rNCStalk_pft(NZ)
  ShootStructE_brch(ielmp) = AZMAX1(LeafStrutElms_brch(ielmp,NB,NZ)+PetoleStrutElms_brch(ielmp,NB,NZ)+StalkLiveBiomassC_brch(NB,NZ))*rPCStalk_pft(NZ)

  IF(iPlantCalendar_brch(ipltcal_EndSeedFill,NB,NZ).EQ.0)THEN
    ShootStructE_brch(ielmn)=ShootStructE_brch(ielmn)+AZMAX1(HuskStrutElms_brch(ielmn,NB,NZ) &
      +EarStrutElms_brch(ielmn,NB,NZ)+GrainStrutElms_brch(ielmn,NB,NZ))

    ShootStructE_brch(ielmp)=ShootStructE_brch(ielmp)+AZMAX1(HuskStrutElms_brch(ielmp,NB,NZ) &
      +EarStrutElms_brch(ielmp,NB,NZ)+GrainStrutElms_brch(ielmp,NB,NZ))      
  ENDIF

  call PrintInfo('end '//subname)
  end associate
  end subroutine UpdateBranchAllometry  
!------------------------------------------------------------------------------------------
  subroutine CalcPartitionCoeff(I,J,NB,NZ,PART,PTRT,BegRemoblize,LRemob_brch)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  integer, intent(out) :: BegRemoblize,LRemob_brch
  real(r8), intent(out):: PART(pltpar%NumOfPlantMorphUnits),PTRT

  character(len=*), parameter :: subname='CalcPartitionCoeff'
  REAL(R8) :: ARLFI
  integer :: N,LP
  real(r8) :: FPARTL
  real(r8) :: PARTS
  real(r8) :: PARTX
  real(r8) :: TOTAL
  real(r8) :: PSILY(0:3)
  real(r8), parameter :: FPART1=1.00_r8
  real(r8), parameter :: FPART2=0.40_r8

  associate(                                                                          &
    TdegCCanopy_pft                   => plt_ew%TdegCCanopy_pft,                      &
    PSICanopy_pft                     => plt_ew%PSICanopy_pft,                        &
    StalkLiveBiomassC_brch            => plt_biom%StalkLiveBiomassC_brch,             &
    StalkRsrvElms_brch                => plt_biom%StalkRsrvElms_brch,                 &
    ZERO                              => plt_site%ZERO,                               &
    NU                                => plt_site%NU,                                 &
    AREA3                             => plt_site%AREA3,                              &
    FracHour4LeafoffRemob             => plt_allom%FracHour4LeafoffRemob,             &
    HourReq4LeafOff_brch              => plt_pheno%HourReq4LeafOff_brch,              &
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch,               &
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft,             &
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch, &
    iPlantDevelopPattern_pft          => plt_pheno%iPlantDevelopPattern_pft,          &
    HoursDoingRemob_brch              => plt_pheno%HoursDoingRemob_brch,              &
    iPlantTurnoverPattern_pft         => plt_pheno%iPlantTurnoverPattern_pft,         &
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch,                &
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch,     &
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft,           &
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft,              &
    TCChill4Seed_pft                  => plt_pheno%TCChill4Seed_pft,                  &
    NH3Dep2Can_brch                   => plt_rbgc%NH3Dep2Can_brch,                    &
    CO2NetFix_pft                     => plt_bgcr%CO2NetFix_pft,                      &
    NodeLenPergC                      => plt_morph%NodeLenPergC,                      &
    CanopyLeafArea_pft                => plt_morph%CanopyLeafArea_pft,                &
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft,                 &
    PARTS_brch                        => plt_morph%PARTS_brch                         &
  )

  call PrintInfo('beg '//subname)
  PSILY=real((/-200.0_r8,-2.0,-2.0,-2.0/),r8)

!     begin_execution

!
!     PARTITION GROWTH WITHIN EACH BRANCH FROM GROWTH STAGE
!     1=LEAF,2=SHEATH OR PETIOLE,3=STALK,4=RESERVE,
!     5,6=REPRODUCTIVE ORGANS,7=GRAIN
!
!     PART=organ partitioning fraction
!

  TOTAL=0._r8
  PART(1:NumOfPlantMorphUnits)=0._r8
!
!     IF BEFORE FLORAL INDUCTION
!
!     iPlantCalendar_brch(ipltcal_InitFloral,=floral initiation date
  
  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
    PART(ibrch_leaf)   = 0.725_r8
    PART(ibrch_petole) = 0.275_r8
!
!     IF BEFORE ANTHESIS
!
!     iPlantCalendar_brch(ipltcal_Anthesis,=start of anthesis and setting final seed number
!     TotalNodeNumNormByMatgrp_brch=total change in vegv node number normalized for maturity group
!
  ELSEIF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
    PART(ibrch_leaf)   = AMAX1(PART1X,0.725_r8-FPART1*TotalNodeNumNormByMatgrp_brch(NB,NZ))
    PART(ibrch_petole) = AMAX1(PART2X,0.275_r8-FPART2*TotalNodeNumNormByMatgrp_brch(NB,NZ))
    PARTS              = 1.0_r8-PART(ibrch_leaf)-PART(ibrch_petole)
    PART(ibrch_stalk)  = 0.60_r8*PARTS
    PART(ibrch_resrv)  = 0.30_r8*PARTS
    PARTX              = PARTS-PART(ibrch_stalk)-PART(ibrch_resrv)
    PART(ibrch_husk)   = 0.5_r8*PARTX
    PART(ibrch_ear)    = 0.5_r8*PARTX
!
!     IF BEFORE GRAIN FILLING, DETERMINATE OR INDETERMINATE
!
!     iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
!     iPlantDevelopPattern_pft=growth habit:0=determinate,1=indetermimate from PFT file
!     TotReproNodeNumNormByMatrgrp_brch=total change in reprv node number normalized for maturity group
!
  ELSEIF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
    IF(iPlantDevelopPattern_pft(NZ).EQ.ideterminate)THEN
      PART(ibrch_leaf)   = 0._r8
      PART(ibrch_petole) = 0._r8
    ELSE
      PART(ibrch_leaf)   = AMAX1(PART1X,(0.725_r8-FPART1)*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
      PART(ibrch_petole) = AMAX1(PART2X,(0.275_r8-FPART2)*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
    ENDIF
    PARTS             = 1.0_r8-PART(ibrch_leaf)-PART(ibrch_petole)
    PART(ibrch_stalk) = AZMAX1(0.60_r8*PARTS*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
    PART(ibrch_resrv) = AZMAX1(0.30_r8*PARTS*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
    PARTX             = PARTS-PART(ibrch_stalk)-PART(ibrch_resrv)
    PART(ibrch_husk)  = 0.5_r8*PARTX
    PART(ibrch_ear)   = 0.5_r8*PARTX
!
!     DURING GRAIN FILLING, DETERMINATE OR INDETERMINATE
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iPlantDevelopPattern_pft=growth habit:0=determinate,1=indetermimate
!
  ELSE
    IF(iPlantDevelopPattern_pft(NZ).EQ.ideterminate)THEN
      PART(ibrch_grain)=1.0_r8
    ELSE
      PART(ibrch_leaf)   = PART1X
      PART(ibrch_petole) = PART2X
      PARTS              = 1.0_r8-PART(ibrch_leaf)-PART(ibrch_petole)
      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
        PART(ibrch_stalk) = 0.125_r8*PARTS
        PART(ibrch_husk)  = 0.125_r8*PARTS
        PART(ibrch_ear)   = 0.125_r8*PARTS
        PART(ibrch_grain) = 0.625_r8*PARTS
      ELSE
        PART(ibrch_stalk) = 0.75_r8*PARTS
        PART(ibrch_grain) = 0.25_r8*PARTS
      ENDIF
    ENDIF
  ENDIF
!
!     IF AFTER GRAIN FILLING
!
!     iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantCalendar_brch(ipltcal_EndSeedFill,=physiological maturity date
!
  IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .AND. iPlantCalendar_brch(ipltcal_EndSeedFill,NB,NZ).NE.0)THEN
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
      PART(ibrch_resrv) = 0._r8
      PART(ibrch_stalk) = 0._r8
      PART(ibrch_grain) = 0._r8
    ELSE
      PART(ibrch_resrv) = PART(ibrch_resrv)+PART(ibrch_stalk)
      PART(ibrch_stalk) = 0._r8
      PART(ibrch_grain) = 0._r8
    ENDIF
  ENDIF
!
!     REDIRECT FROM STALK TO STALK RESERVES IF RESERVES BECOME LOW
!
!     WTRSVB,StalkLiveBiomassC_brch=stalk reserve,sapwood mass
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!
  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)THEN
    IF(StalkRsrvElms_brch(ielmc,NB,NZ).LT.XFRX*StalkLiveBiomassC_brch(NB,NZ))THEN
      D1020: DO N=1,NumOfPlantMorphUnits
        IF(N.NE.ibrch_resrv)THEN
          PART(ibrch_resrv) = PART(ibrch_resrv)+0.10_r8*PART(N)
          PART(N)           = PART(N)-0.10_r8*PART(N)
        ENDIF
      ENDDO D1020
!
!     REDIRECT FROM STALK RESERVES TO STALK IF RESERVES BECOME TOO LARGE
!
    ELSEIF(StalkRsrvElms_brch(ielmc,NB,NZ).GT.1.0_r8*StalkLiveBiomassC_brch(NB,NZ))THEN
      PART(ibrch_stalk) = PART(ibrch_stalk)+PART(ibrch_resrv)+PART(ibrch_grain)
      PART(ibrch_resrv) = 0._r8
      PART(ibrch_grain) = 0._r8
    ENDIF
  ENDIF
!
!     REDIRECT FROM LEAVES TO STALK IF LAI BECOMES TOO LARGE
!
!     CanopyLeafArea_pft=PFT leaf area
!
  ARLFI=CanopyLeafArea_pft(NZ)/AREA3(NU)
  IF(ARLFI.GT.5.0_r8)THEN
    FPARTL             = AZMAX1((10.0_r8-ARLFI)/5.0_r8)
    PART(ibrch_stalk)  = PART(ibrch_stalk)+(1.0_r8-FPARTL)*(PART(ibrch_leaf)+PART(ibrch_petole))
    PART(ibrch_leaf)   = FPARTL*PART(ibrch_leaf)
    PART(ibrch_petole) = FPARTL*PART(ibrch_petole)
  ENDIF

  IF(NB.EQ.MainBranchNum_pft(NZ))THEN
    PTRT=PART(ibrch_leaf)+PART(ibrch_petole)
  ENDIF
!
!     DECIDUOUS LEAF FALL AFTER GRAIN FILL IN DETERMINATES,
!     AFTER AUTUMNIZATION IN INDETERMINATES, OR AFTER SUSTAINED
!     WATER STRESS
!
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!     FracHour4LeafoffRemob=fraction of hours required for leafoff to initiate remobilization
!     iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date for setting final seed number
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     LRemob_brch,BegRemoblize=remobilization flags
!     FLGZ=control rate of remobilization
!
  IF((iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. &
     Hours4LeafOff_brch(NB,NZ).GE.FracHour4LeafoffRemob(iPlantPhenolType_pft(NZ))*HourReq4LeafOff_brch(NB,NZ)) &
    .OR. (iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).NE.0))THEN
    !set remobilization true
    BegRemoblize = itrue
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .OR. iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen)THEN
      !annual plant or evergreen perennial
      LRemob_brch                 = itrue
      HoursDoingRemob_brch(NB,NZ) = HoursDoingRemob_brch(NB,NZ)+1.0_r8
    ELSEIF((iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid .OR. &
      iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid) .AND. &
      TdegCCanopy_pft(NZ).LT.TCChill4Seed_pft(NZ))THEN
      LRemob_brch                 = itrue
      HoursDoingRemob_brch(NB,NZ) = HoursDoingRemob_brch(NB,NZ)+1.0_r8
    ELSEIF(iPlantPhenolType_pft(NZ).GE.2 .AND. PSICanopy_pft(NZ).LT.PSILY(iPlantRootProfile_pft(NZ)))THEN
      LRemob_brch                 = itrue
      HoursDoingRemob_brch(NB,NZ) = HoursDoingRemob_brch(NB,NZ)+1.0_r8
    ENDIF
    IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. iPlantPhenolType_pft(NZ).NE.iphenotyp_evgreen)THEN
      PART(ibrch_stalk)  = PART(ibrch_stalk)+0.5_r8*(PART(ibrch_leaf)+PART(ibrch_petole))
      PART(ibrch_resrv)  = PART(ibrch_resrv)+0.5_r8*(PART(ibrch_leaf)+PART(ibrch_petole))
      PART(ibrch_leaf)   = 0._r8
      PART(ibrch_petole) = 0._r8
    ENDIF
  ELSE
    BegRemoblize                = ifalse
    LRemob_brch                 = ifalse
    HoursDoingRemob_brch(NB,NZ) = 0._r8
  ENDIF
!
!     CHECK PARTITIONING COEFFICIENTS
!
  D1000: DO N=1,NumOfPlantMorphUnits
    IF(N.EQ.ibrch_stalk.AND.isclose(NodeLenPergC(NZ),0._r8))THEN
      PART(N)=0._r8
    ELSE
      PART(N)=AZMAX1(PART(N))
    ENDIF
    TOTAL=TOTAL+PART(N)
  ENDDO D1000
  
  IF(TOTAL.GT.ZERO)THEN
    D1010: DO N=1,NumOfPlantMorphUnits
      PART(N)=PART(N)/TOTAL
    ENDDO D1010
  ELSE
    D1015: DO N=1,NumOfPlantMorphUnits
      PART(N)=0._r8
    ENDDO D1015
  ENDIF

  PARTS_brch(:,NB,NZ)=PART  
  call PrintInfo('end '//subname)
  end associate
  end subroutine CalcPartitionCoeff

!------------------------------------------------------------------------------------------

  subroutine UpdatePhotosynthates(I,J,NB,NZ,TFN6_vr,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
    ShootStructN,TFN5,WaterStress4Groth,Stomata_Stress,CanTurgPSIFun4Expans,CH2O3,CH2O4,CNPG,&
    RCO2NonstC_brch,RCO2Maint_brch,RMxess_brch,NonstC4Groth_brch,RCO2NonstC4Nassim_brch)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: DMSHD
  real(r8), intent(in) :: CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX
  real(r8), intent(in) :: ShootStructN,TFN5,WaterStress4Groth
  real(r8), intent(in) :: Stomata_Stress,CanTurgPSIFun4Expans
  real(r8), intent(out) :: RCO2NonstC_brch
  real(r8), intent(out) :: RCO2Maint_brch
  real(r8), intent(out) :: CH2O3(pltpar%MaxNodesPerBranch1)
  real(r8), intent(out) :: CH2O4(pltpar%MaxNodesPerBranch1)
  REAL(R8), INTENT(OUT) :: CNPG
  real(r8), intent(out) :: RMxess_brch
  real(r8), intent(out) :: NonstC4Groth_brch
  real(r8), intent(out) :: RCO2NonstC4Nassim_brch

  character(len=*), parameter :: subname='UpdatePhotosynthates'
  real(r8) :: CO2F,CanopyNonstElm4Gros(NumPlantChemElms),CH2O
  real(r8) :: CH2OClm,CH2OLlm,dNonstCX
  integer :: K
  associate(                                                    &
    NH3Dep2Can_brch        =>  plt_rbgc%NH3Dep2Can_brch       , &
    CO2FixCL_pft           =>  plt_rbgc%CO2FixCL_pft          , &
    CO2FixLL_pft           =>  plt_rbgc%CO2FixLL_pft          , &        
    CanopyNonstElms_brch   =>  plt_biom%CanopyNonstElms_brch  , &
    RubiscoActivity_brch   =>  plt_photo%RubiscoActivity_brch , &    
    iPlantCalendar_brch    =>  plt_pheno%iPlantCalendar_brch    &
  )
! begin_execution

  call PrintInfo('beg '//subname)
! FDBK=N,P feedback inhibition on C3 CO2 fixation
!
  IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0)THEN
!  
    call ComputeGPP(I,J,NB,NZ,WaterStress4Groth,Stomata_Stress,CH2O3,CH2O4,CH2O,CO2F,CH2OClm,CH2OLlm)

!   SHOOT AUTOTROPHIC RESPIRATION AFTER EMERGENCE
!
    call ComputRAutoAfEmergence(I,J,NB,NZ,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
      CO2F,CH2O,TFN5,WaterStress4Groth,CanTurgPSIFun4Expans,ShootStructN,CanopyNonstElm4Gros,CNPG,&
      RCO2NonstC_brch,RCO2Maint_brch,RMxess_brch,&
      NonstC4Groth_brch,RCO2NonstC4Nassim_brch)

    CO2FixCL_pft(NZ) = CO2FixCL_pft(NZ)+CH2OClm   !carbon-dependent photosynthesis
    CO2FixLL_pft(NZ) = CO2FixLL_pft(NZ)+CH2OLlm   !light-dependent photosynthesis
!   SHOOT AUTOTROPHIC RESPIRATION BEFORE EMERGENCE
!
  ELSE
    CH2O  = 0._r8
    CH2O3 = 0._r8
    CH2O4 = 0._r8
    call ComputRAutoB4Emergence(I,J,NB,NZ,TFN6_vr,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,ShootStructN,&
      WaterStress4Groth,CanTurgPSIFun4Expans,CanopyNonstElm4Gros,CNPG,RCO2NonstC_brch,&
      RCO2Maint_brch,RMxess_brch,NonstC4Groth_brch,RCO2NonstC4Nassim_brch)
  ENDIF

!   REMOVE C,N,P USED IN MAINTENANCE + GROWTH REPIRATION AND GROWTH
!   FROM NON-STRUCTURAL POOLS
!
!   CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!   CH2O=total CH2O production
!   RCO2Maint_brch=maintenance respiration
!   RCO2NonstC_brch=respiration from non-structural C
!   NonstC4Groth_brch=total non-structural C used in growth and respiration
!   RCO2NonstC4Nassim_brch=respiration for N assimilation
!   CanopyNonstElm4Gros(ielmn),CanopyNonstElm4Gros(ielmp)=nonstructural N,P used in growth
!   NH3Dep2Can_brch=NH3 flux between atmosphere and branch from uptake.f
!   XFRE(ielmc),XFRE(ielmn),XFRE(ielmp)=branch-root layer C,N,P transfer
!
  dNonstCX                          = AMIN1(RCO2Maint_brch,RCO2NonstC_brch)+NonstC4Groth_brch+RCO2NonstC4Nassim_brch
  CanopyNonstElms_brch(ielmc,NB,NZ) = CanopyNonstElms_brch(ielmc,NB,NZ)+CH2O-dNonstCX
  CanopyNonstElms_brch(ielmn,NB,NZ) = CanopyNonstElms_brch(ielmn,NB,NZ)-CanopyNonstElm4Gros(ielmn)+NH3Dep2Can_brch(NB,NZ)
  CanopyNonstElms_brch(ielmp,NB,NZ) = CanopyNonstElms_brch(ielmp,NB,NZ)-CanopyNonstElm4Gros(ielmp)

  call PrintInfo('end '//subname)
  end associate
  end subroutine UpdatePhotosynthates
!------------------------------------------------------------------------------------------

  subroutine C4PhotoProductTransfer(I,J,NZ,NB,CH2O3,CH2O4)
  implicit none
  integer, intent(in) :: I,J,NZ,NB
  real(r8), intent(in) :: CH2O3(pltpar%MaxNodesPerBranch1),CH2O4(pltpar%MaxNodesPerBranch1)
  integer :: K
  real(r8) :: CPL4M,CCBS,CPL3K,aquCO2IntraLeafLeakFromBndsheth

!     begin_execution
  associate(                                                            &
    CO2NetFix_pft              => plt_bgcr%CO2NetFix_pft,               &
    CanopyGrosRCO2_pft         => plt_bgcr%CanopyGrosRCO2_pft,          &
    CanopyRespC_CumYr_pft      => plt_bgcr%CanopyRespC_CumYr_pft,       &
    Eco_AutoR_CumYr_col        => plt_bgcr%Eco_AutoR_CumYr_col,         &
    ECO_ER_col                 => plt_bgcr%ECO_ER_col,                  &
    LeafElmntNode_brch         => plt_biom%LeafElmntNode_brch,          &
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft,              &
    CMassHCO3BundleSheath_node => plt_photo%CMassHCO3BundleSheath_node, &
    CMassCO2BundleSheath_node  => plt_photo%CMassCO2BundleSheath_node,  &
    CPOOL3_node                => plt_photo%CPOOL3_node,                &
    CPOOL4_node                => plt_photo%CPOOL4_node,                &
    aquCO2Intraleaf_pft        => plt_photo%aquCO2Intraleaf_pft         &
  )
  D170: DO K=1,MaxNodesPerBranch1
    IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
!
!     MESOPHYLL TO BUNDLE SHEATH TRANSFER
!
!     WGLF=node leaf C mass
!     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
!     CH2O3,CH2O4=total CO2 fixation in bundle sheath,mesophyll
!     CPL4M=mesophyll to bundle sheath transfer of nonstructural C4
!     FBS,FMP=leaf water content in bundle sheath, mesophyll
!
      CPOOL3_node(K,NB,NZ)=CPOOL3_node(K,NB,NZ)-CH2O3(K)
      CPOOL4_node(K,NB,NZ)=CPOOL4_node(K,NB,NZ)+CH2O4(K)
      CPL4M=1.0_r8*(CPOOL4_node(K,NB,NZ)*LeafElmntNode_brch(ielmc,K,NB,NZ)*FBS &
        -CPOOL3_node(K,NB,NZ)*LeafElmntNode_brch(ielmc,K,NB,NZ)*FMP) &
        /(LeafElmntNode_brch(ielmc,K,NB,NZ)*(FBS+FMP))
      CPOOL4_node(K,NB,NZ)=CPOOL4_node(K,NB,NZ)-CPL4M
      CPOOL3_node(K,NB,NZ)=CPOOL3_node(K,NB,NZ)+CPL4M
!
!     BUNDLE SHEATH CO2 DECARBOXYLATION
!
!     CCBS=CO2 concn in bundle sheath (uM)
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!     WGLF=node leaf C mass
!     FBS,FMP=leaf water content in bundle sheath, mesophyll
!     CPL3K=bundle sheath CO2 decarboxylation
!     CO2KI=Ki for C3 leakage from bundle sheath to mesophyll in C4 (uM)
!     FCMassCO2BundleSheath_node,FCMassHCO3BundleSheath_node=partition decarboxylation to CO2,HCO3
!     CPOOL3_node=C4 nonstructural C mass in bundle sheath
!
      CCBS=AZMAX1(0.083E+09*CMassCO2BundleSheath_node(K,NB,NZ)/(LeafElmntNode_brch(ielmc,K,NB,NZ)*FBS))
      CPL3K=2.5E-02_r8*CPOOL3_node(K,NB,NZ)/(1.0_r8+CCBS/CO2KI)
      CPOOL3_node(K,NB,NZ)=CPOOL3_node(K,NB,NZ)-CPL3K
      CMassCO2BundleSheath_node(K,NB,NZ)=CMassCO2BundleSheath_node(K,NB,NZ)+FCMassCO2BundleSheath_node*CPL3K
      CMassHCO3BundleSheath_node(K,NB,NZ)=CMassHCO3BundleSheath_node(K,NB,NZ)+FCMassHCO3BundleSheath_node*CPL3K
!
!     BUNDLE SHEATH LEAKAGE
!
!     aquCO2IntraLeafLeakFromBndsheth=bundle sheath CO2 leakage, >0, leak
!     CCBS=CO2 concn in bundle sheath (uM)
!     aquCO2Intraleaf_pft=intercellular CO2 concentration (uM)
!     WGLF=node leaf C mass
!     FBS=leaf water content in bundle sheath
!     FCMassCO2BundleSheath_node,FCMassHCO3BundleSheath_node=partition decarboxylation to CO2,HCO3
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!
      aquCO2IntraLeafLeakFromBndsheth=5.0E-07_r8*(CCBS-aquCO2Intraleaf_pft(NZ))*&
        LeafElmntNode_brch(ielmc,K,NB,NZ)*FBS
      CMassCO2BundleSheath_node(K,NB,NZ)=CMassCO2BundleSheath_node(K,NB,NZ)-FCMassCO2BundleSheath_node*&
        aquCO2IntraLeafLeakFromBndsheth
      CMassHCO3BundleSheath_node(K,NB,NZ)=CMassHCO3BundleSheath_node(K,NB,NZ)-FCMassHCO3BundleSheath_node*&
        aquCO2IntraLeafLeakFromBndsheth

!
!     TOTAL C EXCHANGE
!
!     CanopyGrosRCO2_pft,CanopyRespC_CumYr_pft=total,above-ground PFT respiration
!     CO2NetFix_pft=PFT net CO2 fixation
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_CumYr_col=total autotrophic respiration
!     aquCO2IntraLeafLeakFromBndsheth=bundle sheath CO2 leakage
!
      CanopyGrosRCO2_pft(NZ)    = CanopyGrosRCO2_pft(NZ)-aquCO2IntraLeafLeakFromBndsheth
      CanopyRespC_CumYr_pft(NZ) = CanopyRespC_CumYr_pft(NZ)-aquCO2IntraLeafLeakFromBndsheth
      CO2NetFix_pft(NZ)         = CO2NetFix_pft(NZ)-aquCO2IntraLeafLeakFromBndsheth
      ECO_ER_col                = ECO_ER_col-aquCO2IntraLeafLeakFromBndsheth
      Eco_AutoR_CumYr_col       = Eco_AutoR_CumYr_col-aquCO2IntraLeafLeakFromBndsheth
    ENDIF
  ENDDO D170
  end associate
  end subroutine C4PhotoProductTransfer
!------------------------------------------------------------------------------------------

  subroutine RemobizeLeafNodes(NB,NZ,KN,RCCC,RCCN,RCCP,SenesFrac,SNCT,lgoto565)
  implicit none
  integer,intent(in) :: NB,NZ,KN
  real(r8), intent(in) :: RCCC,RCCN,RCCP
  real(r8), intent(in):: SenesFrac  
  real(r8), intent(inout) :: SNCT
  logical,intent(out):: lgoto565

  real(r8) :: FracSensResp4Leaf
  real(r8) :: FSNAS,FracNodeAsLeaf,FSNCS  
  real(r8) :: SensResp2Use,Frac2Remobl,FSNAL
  integer :: K,KK,NE,M
  real(r8) :: RCES(NumPlantChemElms),dleaf,dPetole,dmassE  
  real(r8) :: ElmntRemobFallingLeaf(NumPlantChemElms)

  !begin_execution
  associate(                                                              &
    KHiestGroLeafNode_brch      => plt_pheno%KHiestGroLeafNode_brch,      &
    PetioleElmntNode_brch       => plt_biom%PetioleElmntNode_brch,        &
    LitrfalStrutElms_pvr        => plt_bgcr%LitrfalStrutElms_pvr,         &
    CPOOL3_node                 => plt_photo%CPOOL3_node,                 &
    CPOOL4_node                 => plt_photo%CPOOL4_node,                 &
    inonfoliar                  => pltpar%inonfoliar,                     &
    icarbhyro                   => pltpar%icarbhyro,                      &
    icwood                      => pltpar%icwood,                         &
    ifoliar                     => pltpar%ifoliar,                        &
    FracShootLeafElmAlloc2Litr  => plt_allom%FracShootLeafElmAlloc2Litr,  &
    PetoleProteinCNode_brch     => plt_biom%PetoleProteinCNode_brch,      &
    LeafProteinCNode_brch       => plt_biom%LeafProteinCNode_brch,        &
    FracShootStalkElmAlloc2Litr => plt_allom%FracShootStalkElmAlloc2Litr, &
    LeafAreaLive_brch           => plt_morph%LeafAreaLive_brch,           &
    LeafNodeArea_brch           => plt_morph%LeafNodeArea_brch,           &
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch,           &
    LeafElmntNode_brch          => plt_biom%LeafElmntNode_brch,           &
    ElmAllocmat4Litr            => plt_soilchem%ElmAllocmat4Litr,         &
    PetoleStrutElms_brch        => plt_biom%PetoleStrutElms_brch,         &
    PetoleLensNode_brch         => plt_morph%PetoleLensNode_brch,         &
    rCNNonstRemob_pft           => plt_allom%rCNNonstRemob_pft,           &
    rCPNonstRemob_pft           => plt_allom%rCPNonstRemob_pft,           &
    k_fine_litr                 => pltpar%k_fine_litr,                    &
    k_woody_litr                => pltpar%k_woody_litr,                   &
    ZERO4LeafVar_pft            => plt_biom%ZERO4LeafVar_pft,             &
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft,               &
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch          &
  )  

  lgoto565=.false.
  D650: DO KK=KN,KHiestGroLeafNode_brch(NB,NZ)
    FracSensResp4Leaf=0._r8
    SensResp2Use=0._r8
    K=pMOD(KK,MaxNodesPerBranch1)
!
!       REMOBILIZATION OF LEAF C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
!
!       WGLF,PetioleElmntNode_brch=node leaf,petiole C mass
  !       SCNF,SCNSH=leaf,petiole senescence respiration
  !       ElmntRemobFallingLeaf(:)=remobilization of C,N,P from senescing leaf
  !       RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
  !       RCCZ,RCCY=min,max fractions for shoot C recycling
  !
    IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      FracNodeAsLeaf=LeafElmntNode_brch(ielmc,K,NB,NZ)/(LeafElmntNode_brch(ielmc,K,NB,NZ)+&
        PetioleElmntNode_brch(ielmc,K,NB,NZ))
      FracSensResp4Leaf=FracNodeAsLeaf*SNCT
      SensResp2Use=SNCT-FracSensResp4Leaf
      ElmntRemobFallingLeaf(ielmc)=RCCC*LeafElmntNode_brch(ielmc,K,NB,NZ)
      ElmntRemobFallingLeaf(ielmn)=LeafElmntNode_brch(ielmn,K,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
      ElmntRemobFallingLeaf(ielmp)=LeafElmntNode_brch(ielmp,K,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
!
  !         FRACTION OF CURRENT LEAF TO BE REMOBILIZED
  !
  !         Frac2Remobl,FSNAL=fraction of current leaf C,area to be remobilized
  !
      IF(ElmntRemobFallingLeaf(ielmc).GT.ZERO4Groth_pft(NZ))THEN
        Frac2Remobl=AZMAX1(AMIN1(1.0_r8,FracSensResp4Leaf/ElmntRemobFallingLeaf(ielmc)))
      ELSE
        Frac2Remobl=1.0_r8
      ENDIF
      FSNAL=Frac2Remobl
!
      !     NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
      !     TO FRACTIONS SET IN 'STARTQ'
      !
      !     CSNC,ZSNC,PSNC=literfall C,N,P
      !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
      !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
      !     Frac2Remobl=fraction of current leaf to be remobilized
      !     WGLF,WGLFN,WGLFP=node leaf C,N,P mass
      !     ElmntRemobFallingLeaf(:)=remobilization of C,N,P from senescing leaf
      !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
      !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
      !     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
      !
      
      D6310: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dleaf=AZMAX1(LeafElmntNode_brch(NE,K,NB,NZ)-ElmntRemobFallingLeaf(NE))*Frac2Remobl
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
            +ElmAllocmat4Litr(NE,icwood,M,NZ)*dleaf*FracShootStalkElmAlloc2Litr(NE,k_woody_litr)

          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*dleaf*FracShootStalkElmAlloc2Litr(NE,k_fine_litr)
        ENDDO
      ENDDO D6310
      IF(K.NE.0)THEN
        LitrfalStrutElms_pvr(ielmc,icarbhyro,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,icarbhyro,k_fine_litr,0,NZ) &
          +Frac2Remobl*AZMAX1(CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ))
        CPOOL3_node(K,NB,NZ)=CPOOL3_node(K,NB,NZ)*AZMAX1(1._r8-Frac2Remobl)
        CPOOL4_node(K,NB,NZ)=CPOOL4_node(K,NB,NZ)*AZMAX1(1._r8-Frac2Remobl)
      ENDIF
!
!     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
!
!     LeafAreaLive_brch=leaf area
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WGLF,WGLFN,WGLFP,LeafProteinCNode_brch=node leaf C,N,P,protein mass
!     Frac2Remobl=fraction of current leaf to be remobilized
!     CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios from startq.f
!
      LeafAreaLive_brch(NB,NZ)   = AZMAX1(LeafAreaLive_brch(NB,NZ)-FSNAL*LeafNodeArea_brch(K,NB,NZ))
      LeafNodeArea_brch(K,NB,NZ) = AZMAX1(LeafNodeArea_brch(K,NB,NZ)-FSNAL*LeafNodeArea_brch(K,NB,NZ))

      DO NE=1,NumPlantChemElms
        LeafStrutElms_brch(NE,NB,NZ)=AZMAX1(LeafStrutElms_brch(NE,NB,NZ)-Frac2Remobl*LeafElmntNode_brch(NE,K,NB,NZ))
        LeafElmntNode_brch(NE,K,NB,NZ)=AZMAX1(LeafElmntNode_brch(NE,K,NB,NZ)-Frac2Remobl*LeafElmntNode_brch(NE,K,NB,NZ))
      ENDDO
      LeafProteinCNode_brch(K,NB,NZ)=AZMAX1(LeafProteinCNode_brch(K,NB,NZ) &
        -Frac2Remobl*AMAX1(LeafElmntNode_brch(ielmn,K,NB,NZ)*rCNNonstRemob_pft(NZ) &
        ,LeafElmntNode_brch(ielmp,K,NB,NZ)*rCPNonstRemob_pft(NZ)))
!
!     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
!     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!     Frac2Remobl=fraction of current leaf C to be remobilized
!     ElmntRemobFallingLeaf(ielmc),ElmntRemobFallingLeaf(ielmn),ElmntRemobFallingLeaf(ielmp)=remobilization of C,N,P from senescing leaf
!     FracSensResp4Leaf,SNCT=remaining senescence respiration carried to next node
!
      CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)+Frac2Remobl*ElmntRemobFallingLeaf(ielmc)*SenesFrac
      DO NE=2,NumPlantChemElms
        CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+Frac2Remobl*ElmntRemobFallingLeaf(NE)
      ENDDO
      FracSensResp4Leaf = FracSensResp4Leaf-Frac2Remobl*ElmntRemobFallingLeaf(ielmc)
      SNCT              = SNCT-Frac2Remobl*ElmntRemobFallingLeaf(ielmc)
      IF(LeafStrutElms_brch(ielmc,NB,NZ).LE.ZERO4LeafVar_pft(NZ))THEN
        LeafStrutElms_brch(ielmc,NB,NZ) = 0._r8
        LeafAreaLive_brch(NB,NZ)        = 0._r8
      ENDIF
    !
      !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
      !
!        IF(FracSensResp4Leaf.LE.ZERO4Groth_pft(NZ))GO TO 564
      !
      !     OTHERWISE REMAINING C,N,P IN LEAF GOES TO LitrFall
      !
      !     CSNC,ZSNC,PSNC=literfall C,N,P
      !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
      !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
      !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
      !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
      !     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
      !     LeafAreaLive_brch=leaf area
      !     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
      !     WGLF,WGLFN,WGLFP,LeafProteinCNode_brch=node leaf C,N,P,protein mass
      !
    ELSE        
      D6315: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,icwood,M,NZ)*AZMAX1(LeafElmntNode_brch(NE,K,NB,NZ))*FracShootStalkElmAlloc2Litr(NE,k_woody_litr)
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*AZMAX1(LeafElmntNode_brch(NE,K,NB,NZ))*FracShootStalkElmAlloc2Litr(NE,k_fine_litr)
        ENDDO
      ENDDO D6315

      IF(K.NE.0)THEN
        LitrfalStrutElms_pvr(ielmc,icarbhyro,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(ielmc,icarbhyro,k_fine_litr,0,NZ)&
          +AZMAX1(CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ))
        CPOOL3_node(K,NB,NZ)=0._r8
        CPOOL4_node(K,NB,NZ)=0._r8
      ENDIF
      LeafAreaLive_brch(NB,NZ)=AZMAX1(LeafAreaLive_brch(NB,NZ)-LeafNodeArea_brch(K,NB,NZ))
      DO NE=1,NumPlantChemElms
        LeafStrutElms_brch(NE,NB,NZ)   = AZMAX1(LeafStrutElms_brch(NE,NB,NZ)-LeafElmntNode_brch(NE,K,NB,NZ))
        LeafElmntNode_brch(NE,K,NB,NZ) = 0._r8
      ENDDO
      LeafNodeArea_brch(K,NB,NZ)     = 0._r8
      LeafProteinCNode_brch(K,NB,NZ) = 0._r8
      IF(LeafStrutElms_brch(ielmc,NB,NZ).LE.ZERO4LeafVar_pft(NZ))THEN
        LeafStrutElms_brch(ielmc,NB,NZ) = 0._r8
        LeafAreaLive_brch(NB,NZ)        = 0._r8
      ENDIF
    ENDIF
!564   CONTINUE
!
  !     REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P DEPENDS ON
  !     NON-STRUCTURAL C:N:P
  !
  !     PetioleElmntNode_brch,WGSHN,WGSHP=node petiole C,N,P mass
  !     RCES(ielmc),RCES(ielmn),RCES(ielmp)=remobilization of C,N,P from senescing petiole
  !     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
  !
    IF(PetioleElmntNode_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      RCES(ielmc)=PetioleElmntNode_brch(ielmc,K,NB,NZ)*RCCC
      RCES(ielmn)=PetioleElmntNode_brch(ielmn,K,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
      RCES(ielmp)=PetioleElmntNode_brch(ielmp,K,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
!
    !     FRACTION OF REMOBILIZATION THAT CAN BE MET FROM CURRENT SHEATH
    !     OR PETIOLE
    !
    !     FSNCS,FSNAS=fraction of current petiole C,length to be remobilized
    !
      IF(RCES(ielmc).GT.ZERO4Groth_pft(NZ))THEN
        FSNCS=AZMAX1(AMIN1(1.0_r8,SensResp2Use/RCES(ielmc)))
      ELSE
        FSNCS=1.0_r8
      ENDIF
      FSNAS=FSNCS
  !
    !     NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
    !     TO FRACTIONS SET IN 'STARTQ'
    !
    !     CSNC,ZSNC,PSNC=literfall C,N,P
    !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
    !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
    !     FSNCS=fraction of current petiole to be remobilized
    !     PetioleElmntNode_brch,WGSHN,WGSHP=node petiole C,N,P mass
    !     RCES(ielmc),RCES(ielmn),RCES(ielmp)=remobilization of C,N,P from senescing petiole
    !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
    !     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
    !
      
      D6320: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dPetole=FSNCS*AZMAX1(PetioleElmntNode_brch(NE,K,NB,NZ)-RCES(NE))
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
            +ElmAllocmat4Litr(NE,icwood,M,NZ)*dPetole*FracShootLeafElmAlloc2Litr(NE,k_woody_litr)
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)&
            +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*dPetole*FracShootLeafElmAlloc2Litr(NE,k_fine_litr)
        ENDDO    
      ENDDO D6320
      DO NE=1,NumPlantChemElms
        dmassE=FSNCS*PetioleElmntNode_brch(NE,K,NB,NZ)
        PetoleStrutElms_brch(NE,NB,NZ)    = AZMAX1(PetoleStrutElms_brch(NE,NB,NZ)-dmassE)
        PetioleElmntNode_brch(NE,K,NB,NZ) = AZMAX1(PetioleElmntNode_brch(NE,K,NB,NZ)-dmassE)
      ENDDO
!
    !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
    !
    !     PetoleLensNode_brch=petiole length
    !     WTSHEB,WTLFBN,WTSHBP=branch petiole C,N,P mass
    !     PetioleElmntNode_brch,WGSHN,WGSHP,PetoleProteinCNode_brch=node petiole C,N,P,protein mass
    !     FSNCS=fraction of current petiole to be remobilized
    !     CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios from startq.f
    !
      PetoleLensNode_brch(K,NB,NZ)=PetoleLensNode_brch(K,NB,NZ)-FSNAS*PetoleLensNode_brch(K,NB,NZ)
      PetoleProteinCNode_brch(K,NB,NZ)=AZMAX1(PetoleProteinCNode_brch(K,NB,NZ) &
        -FSNCS*AMAX1(PetioleElmntNode_brch(ielmn,K,NB,NZ)*rCNNonstRemob_pft(NZ) &
        ,PetioleElmntNode_brch(ielmp,K,NB,NZ)*rCPNonstRemob_pft(NZ)))
!
!     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
!     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!     FSNCS=fraction of current petiole C to be remobilized
!     RCES(ielmc),RCES(ielmn),RCES(ielmp)=remobilization of C,N,P from senescing petiole
!     SensResp2Use,SNCT=remaining senescence respiration carried to next node
!
      CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)+FSNCS*RCES(ielmc)*SenesFrac
      DO NE=2,NumPlantChemElms
        CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+FSNCS*RCES(NE)
      ENDDO
      SensResp2Use=SensResp2Use-FSNCS*RCES(ielmc)
      SNCT=SNCT-FSNCS*RCES(ielmc)
      IF(PetoleStrutElms_brch(ielmc,NB,NZ).LE.ZERO4LeafVar_pft(NZ))THEN
        PetoleStrutElms_brch(ielmc,NB,NZ)=0._r8
      ENDIF
!
!     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
!
      lgoto565=(SensResp2Use.LE.ZERO4Groth_pft(NZ))
!
    !     OTHERWISE REMAINING C,N,P IN SHEATH OR PETIOLE GOES TO LitrFall
    !
    !     CSNC,ZSNC,PSNC=literfall C,N,P
    !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
    !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
    !     FWODB=C woody fraction in branch:0=woody,1=non-woody
    !     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
    !     PetoleLensNode_brch=petiole length
    !     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
    !     PetioleElmntNode_brch,WGSHN,WGSHP,PetoleProteinCNode_brch=node petiole C,N,P,protein mass
    !
    ELSE      
      D6325: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,icwood,M,NZ)*AZMAX1(PetioleElmntNode_brch(NE,K,NB,NZ))*FracShootLeafElmAlloc2Litr(NE,k_woody_litr)
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(PetioleElmntNode_brch(NE,K,NB,NZ))*FracShootLeafElmAlloc2Litr(NE,k_fine_litr)
        ENDDO    
      ENDDO D6325
      DO NE=1,NumPlantChemElms
        PetoleStrutElms_brch(NE,NB,NZ)    = AZMAX1(PetoleStrutElms_brch(NE,NB,NZ)-PetioleElmntNode_brch(NE,K,NB,NZ))
        PetioleElmntNode_brch(NE,K,NB,NZ) = 0._r8
      ENDDO
      PetoleLensNode_brch(K,NB,NZ)     = 0._r8
      PetoleProteinCNode_brch(K,NB,NZ) = 0._r8
      IF(PetoleStrutElms_brch(ielmc,NB,NZ).LE.ZERO4LeafVar_pft(NZ))THEN
        PetoleStrutElms_brch(ielmc,NB,NZ)=0._r8
      ENDIF
    ENDIF
  ENDDO D650
  end associate
  END subroutine RemobizeLeafNodes

!------------------------------------------------------------------------------------------
  subroutine RemobilizeLeafLayers(NumRemobLeafNodes,NB,NZ,RespSenesTot_brch,RCCC,RCCN,RCCP,SenesFrac)
  implicit none
  INTEGER, intent(in)  :: NB,NZ,NumRemobLeafNodes
  real(r8), intent(in) :: RespSenesTot_brch
  REAL(R8), INTENT(IN) :: RCCC,RCCN,RCCP
  real(r8), intent(inout):: SenesFrac
  integer :: N,M,K,KK,MXNOD,MNNOD,NE
  real(r8) :: FRCC
  real(r8) :: FSNCK
  real(r8) :: FSNCR
  real(r8) :: HTNODZ
  real(r8) :: totNumRemobLeafNodes_brch  
  real(r8) :: FStalkSenes
  real(r8) :: RCSC,RCSN,RCSP
  real(r8) :: RCEK(NumPlantChemElms)
  real(r8) :: RMxess_brch
  real(r8) :: RespSenesPhenol_brch
  real(r8) :: SNCT,dIntNode,dStalk,dNodeH
  integer :: KN  
  logical :: lgoto565
! begin_execution
  associate(                                                            &
    istalk                     => pltpar%istalk,                        &
    iroot                      => pltpar%iroot,                         &
    k_fine_litr                => pltpar%k_fine_litr,                   &
    k_woody_litr               => pltpar%k_woody_litr,                  &
    StalkRsrvElms_brch         => plt_biom%StalkRsrvElms_brch,          &
    StalkStrutElms_brch        => plt_biom%StalkStrutElms_brch,         &
    InternodeStrutElms_brch    => plt_biom%InternodeStrutElms_brch,     &
    StalkLiveBiomassC_brch     => plt_biom%StalkLiveBiomassC_brch,      &
    PetioleElmntNode_brch      => plt_biom%PetioleElmntNode_brch,       &
    SeasonalNonstElms_pft      => plt_biom%SeasonalNonstElms_pft,       &
    SenecStalkStrutElms_brch   => plt_biom%SenecStalkStrutElms_brch,    &
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft,              &
    PlantPopulation_pft        => plt_site%PlantPopulation_pft,         &
    FracRootStalkElmAlloc2Litr => plt_allom%FracRootStalkElmAlloc2Litr, &
    LitrfalStrutElms_pvr       => plt_bgcr%LitrfalStrutElms_pvr,        &
    ElmAllocmat4Litr           => plt_soilchem%ElmAllocmat4Litr,        &
    KLowestGroLeafNode_brch    => plt_pheno%KLowestGroLeafNode_brch,    &
    KHiestGroLeafNode_brch     => plt_pheno%KHiestGroLeafNode_brch,     &
    iPlantBranchState_brch     => plt_pheno%iPlantBranchState_brch,     &
    iPlantPhenolPattern_pft    => plt_pheno%iPlantPhenolPattern_pft,    &
    NumCogrowthNode_pft        => plt_morph%NumCogrowthNode_pft,        &
    LiveInterNodeHight_brch    => plt_morph%LiveInterNodeHight_brch,    &
    InternodeHeightDead_brch   => plt_morph%InternodeHeightDead_brch    &
  )
!     REMOBILIZATION AND LitrFall WHEN GROWTH RESPIRATION < 0
!     STARTING FROM LOWEST LEAFED NODE AND PROCEEDING UPWARDS
!
!     RespSenesTot_brch,SNCT=branch,node senescence respiration
!     KSNC=number of nodes undergoing remobilization
!
  totNumRemobLeafNodes_brch=NumRemobLeafNodes
  KN=MAX(0,KLowestGroLeafNode_brch(NB,NZ)-1)
  D575: DO N=1,NumRemobLeafNodes
    SNCT=RespSenesTot_brch/totNumRemobLeafNodes_brch

    call RemobizeLeafNodes(NB,NZ,KN,RCCC,RCCN,RCCP,SenesFrac,SNCT,lgoto565)    
    if(lgoto565)cycle
    
    KN=KN+1
    RMxess_brch=SNCT*(1.0_r8-SenesFrac)
!
!     REMOBILIZATION OF RESERVE C
!
!     WTRSVB=stalk reserve C mass
!     RMxess_brch=excess maintenance respiration
!
    IF(StalkRsrvElms_brch(ielmc,NB,NZ).GT.RMxess_brch)THEN
      StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)-RMxess_brch
      RMxess_brch=0._r8
      cycle
    ENDIF
!
!     REMOBILIZATION OF STALK C,N,P
!
!     FXFS=rate constant for remobilization of stalk C,N,P (h-1)
!     RespSenesPhenol_brch=phenologically-driven respiration senescence during late-season
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     WTSTKB,StalkLiveBiomassC_brch=stalk,sapwood C mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     MXNOD,MNNOD=max,min node number currently growing
!     NumCogrowthNode_pft=number of concurrently growing nodes
!     KHiestGroLeafNode_brch=integer of most recent leaf number
!
    RespSenesPhenol_brch=FXFS*RMxess_brch
    SNCT=RMxess_brch+RespSenesPhenol_brch
    IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. SNCT.GT.ZERO4Groth_pft(NZ) &
      .AND.StalkStrutElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      SenesFrac = RespSenesPhenol_brch/SNCT
      FRCC      = StalkLiveBiomassC_brch(NB,NZ)/StalkStrutElms_brch(ielmc,NB,NZ)
      RCSC      = RCCC*FRCC
      RCSN      = RCCN*FRCC
      RCSP      = RCCP*FRCC
      MXNOD     = KHiestGroLeafNode_brch(NB,NZ)
      MNNOD  = MAX(MIN(0,MAX(0,MXNOD-NumCogrowthNode_pft(NZ))),KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+2)
      MXNOD  = MAX(MXNOD,MNNOD)
      D1650: DO KK = MXNOD, MNNOD, -1
        K=pMOD(KK,MaxNodesPerBranch1)
    !
    !     REMOBILIZATION OF STALK C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
    !
    !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
    !     InternodeStrutElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
    !
        IF(InternodeStrutElms_brch(ielmc,K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
          RCEK(ielmc)=InternodeStrutElms_brch(ielmc,K,NB,NZ)*RCSC
          RCEK(ielmn)=InternodeStrutElms_brch(ielmn,K,NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
          RCEK(ielmp)=InternodeStrutElms_brch(ielmp,K,NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
    !
      !     FRACTION OF CURRENT NODE TO BE REMOBILIZED
      !
      !     FSNCS=fraction of lowest internode to be remobilized
!
          IF(RCEK(ielmc).GT.ZERO4Groth_pft(NZ))THEN
            FSNCK=AZMAX1(AMIN1(1.0,SNCT/RCEK(ielmc)))
          ELSE
            FSNCK=1.0_r8
          ENDIF
    !
      !     NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
      !     TO FRACTIONS SET IN 'STARTQ'
      !
      !     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
      !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
      !     FSNCK=fraction of lowest internode to be remobilized
      !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
      !     InternodeStrutElms_brch,WGNODN,WGNODP=senescing internode C,N,P mass
!
          
          D7310: DO M=1,jsken
            DO NE=1,NumPlantChemElms
              dIntNode=AZMAX1(InternodeStrutElms_brch(NE,K,NB,NZ)-RCEK(NE))*FSNCK
              LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
                +ElmAllocmat4Litr(NE,istalk,M,NZ)*dIntNode*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

              LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
                +ElmAllocmat4Litr(NE,istalk,M,NZ)*dIntNode*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
            ENDDO    
          ENDDO D7310

!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
      !
      !     FSNCK=fraction of lowest internode to be remobilized
      !     LiveInterNodeHight_brch,InternodeHeightDead_brch=living,senescing internode length
      !     WTSTKB,WTSTBN,WTSTBP,InternodeStrutElms_brch,WGNODN,WGNODP=C,N,P mass in senescing internode
      !
          DO NE=1,NumPlantChemElms
            dNodeH=FSNCK*InternodeStrutElms_brch(NE,K,NB,NZ)
            StalkStrutElms_brch(NE,NB,NZ)       = AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-dNodeH)
            InternodeStrutElms_brch(NE,K,NB,NZ) = AZMAX1(InternodeStrutElms_brch(NE,K,NB,NZ)-dNodeH)
          ENDDO
          dNodeH=FSNCK*InternodeHeightDead_brch(K,NB,NZ)
          LiveInterNodeHight_brch(K,NB,NZ)   = AZMAX1(LiveInterNodeHight_brch(K,NB,NZ)-dNodeH)
          InternodeHeightDead_brch(K,NB,NZ) = AZMAX1(InternodeHeightDead_brch(K,NB,NZ)-dNodeH)
!
      !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
      !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
      !
      !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
      !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
      !     FSNCK=fraction of lowest internode to be remobilized
      !     SNCT=remaining node senescence respiration
      !
          StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+FSNCK*RCEK(ielmc)*SenesFrac
          DO NE=2,NumPlantChemElms
            StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+FSNCK*RCEK(NE)
          ENDDO
          SNCT=SNCT-FSNCK*RCEK(ielmc)
    !
      !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
      !
          IF(SNCT.LE.ZERO4Groth_pft(NZ))then
            lgoto565=.true.
            exit
          else
            lgoto565=.false.
          endif
    !
    !       OTHERWISE REMAINING C,N,P IN NODE GOES TO LitrFall
    !
    !       CSNC,ZSNC,PSNC=literfall C,N,P
    !       CFOPC=fraction of plant litter allocated in nonstructural(0,*),
    !       WTSTKB,WTSTBN,WTSTBP,InternodeStrutElms_brch,WGNODN,WGNODP=C,N,P mass in senescing internode
    !       LiveInterNodeHight_brch,InternodeHeightDead_brch=living,senescing internode length
    !
        ELSE
          
          D7315: DO M=1,jsken
            DO NE=1,NumPlantChemElms
              LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
                +ElmAllocmat4Litr(NE,istalk,M,NZ)*AZMAX1(InternodeStrutElms_brch(NE,K,NB,NZ))*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

              LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
                +ElmAllocmat4Litr(NE,istalk,M,NZ)*AZMAX1(InternodeStrutElms_brch(NE,K,NB,NZ))*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
            ENDDO    
          ENDDO D7315
          DO NE=1,NumPlantChemElms
            StalkStrutElms_brch(NE,NB,NZ)=AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-InternodeStrutElms_brch(NE,K,NB,NZ))
            InternodeStrutElms_brch(NE,K,NB,NZ)=0._r8
          ENDDO
          LiveInterNodeHight_brch(K,NB,NZ)   = LiveInterNodeHight_brch(K,NB,NZ)-InternodeHeightDead_brch(K,NB,NZ)
          InternodeHeightDead_brch(K,NB,NZ) = 0._r8
        ENDIF
      ENDDO D1650
      if(lgoto565)cycle
!
    !   RESIDUAL STALK
    !
    !   RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
    !   WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
      IF(SenecStalkStrutElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
        RCEK(ielmc)=SenecStalkStrutElms_brch(ielmc,NB,NZ)*RCSC
        RCEK(ielmn)=SenecStalkStrutElms_brch(ielmn,NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
        RCEK(ielmp)=SenecStalkStrutElms_brch(ielmp,NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
    !
    !     FRACTION OF RESIDUAL STALK TO BE REMOBILIZED
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
    !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
        IF(RCEK(ielmc).GT.ZERO4Groth_pft(NZ))THEN
          FSNCR=AZMAX1(AMIN1(1.0,SNCT/RCEK(ielmc)))
        ELSE
          FSNCR=1.0_r8
        ENDIF
    !
    !     NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
    !     TO FRACTIONS SET IN 'STARTQ'
    !
        D8310: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            dStalk=FSNCR*AZMAX1(SenecStalkStrutElms_brch(NE,NB,NZ)-RCEK(NE))
            LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
              +ElmAllocmat4Litr(NE,istalk,M,NZ)*dStalk*FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
              +ElmAllocmat4Litr(NE,istalk,M,NZ)*dStalk*FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
          ENDDO    
        ENDDO D8310

    !
    !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in senescing stalk
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !     LiveInterNodeHight_brch,InternodeHeightDead_brch=living,senescing internode length
    !
        DO NE=1,NumPlantChemElms
          FStalkSenes=FSNCR*SenecStalkStrutElms_brch(NE,NB,NZ)
          StalkStrutElms_brch(NE,NB,NZ)=AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-FStalkSenes)
          SenecStalkStrutElms_brch(NE,NB,NZ)=AZMAX1(SenecStalkStrutElms_brch(NE,NB,NZ)-FStalkSenes)
        ENDDO
        HTNODZ=0._r8
        D8320: DO K=0,MaxNodesPerBranch1
          HTNODZ=AMAX1(HTNODZ,LiveInterNodeHight_brch(K,NB,NZ))
        ENDDO D8320

        HTNODZ=AZMAX1(HTNODZ-FSNCR*HTNODZ)
        D8325: DO K=0,MaxNodesPerBranch1
          LiveInterNodeHight_brch(K,NB,NZ)=AMIN1(HTNODZ,LiveInterNodeHight_brch(K,NB,NZ))
        ENDDO D8325
!
    !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
    !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
    !
    !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
    !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
    !     FSNCR=fraction of residual stalk to be remobilized
    !     SNCT=remaining node senescence respiration
    !
        StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+FSNCR*RCEK(ielmc)*SenesFrac
        DO NE=2,NumPlantChemElms
          StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+FSNCR*RCEK(NE)
        ENDDO
        SNCT=SNCT-FSNCR*RCEK(ielmc)
      ENDIF
  !
  !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
  !
      IF(SNCT.LE.ZERO4Groth_pft(NZ))cycle
  !
  !     OTHERWISE REMAINING C,N,P IN NODE GOES TO LitrFall
  !
  !     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
  !     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
  !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
  !
    ELSE      
      D8315: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
            +ElmAllocmat4Litr(NE,istalk,M,NZ)*AZMAX1(SenecStalkStrutElms_brch(NE,NB,NZ)) &
            *FracRootStalkElmAlloc2Litr(NE,k_woody_litr)

          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,istalk,M,NZ)*AZMAX1(SenecStalkStrutElms_brch(NE,NB,NZ)) &
            *FracRootStalkElmAlloc2Litr(NE,k_fine_litr)
        ENDDO    
      ENDDO D8315
      DO NE=1,NumPlantChemElms
        StalkStrutElms_brch(NE,NB,NZ)=AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-SenecStalkStrutElms_brch(NE,NB,NZ))
        SenecStalkStrutElms_brch(NE,NB,NZ)=0._r8
      ENDDO
!
!     REMOBILIZATION OF STORAGE C,N,P
!
!     WTRVC=storage C
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     RMxess_brch=remaining excess maintenance respiration
!
        RMxess_brch=SNCT/(1.0_r8+FXFS)
        IF(SeasonalNonstElms_pft(ielmc,NZ).GT.RMxess_brch)THEN
          SeasonalNonstElms_pft(ielmc,NZ) = SeasonalNonstElms_pft(ielmc,NZ)-RMxess_brch
          RMxess_brch                     = 0._r8
        ELSEIF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
          iPlantBranchState_brch(NB,NZ)=iDead
        ENDIF      
    ENDIF

565 CONTINUE
  ENDDO D575
  end associate
  end subroutine RemobilizeLeafLayers


!------------------------------------------------------------------------------------------

  subroutine AllocateLeaf2CanopyLayers(I,J,NB,NZ,CanopyHeight_copy)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NB,NZ
  real(r8), intent(in) :: CanopyHeight_copy(JP1)
  integer  :: LL,LU,L,K,k1,k2,KK,NE
  integer  :: KMinGroingLeafNodeNum,KLeafNumHighestGrowing
  integer  :: LNumHeightLeafTip,LNumHeightLeafBase
  integer  :: LNumHeightBranchTip,LNumHeightBranchBase,N
  real(r8) :: RadiusInnerStalk4Transp
  real(r8) :: YLeafElmntNode_brch(NumPlantChemElms)
  real(r8) :: YLeafArea_brch,LeafElevation,LeafLength
  real(r8) :: StalkSurfArea,StalkSectionArea
  real(r8) :: FRACL
  real(r8) :: HeightBranchBase
  real(r8) :: HeightStalk
  real(r8) :: HeightLeafNode,HeightLeafLow,HeightLeafTip
  real(r8) :: HeightLeafBase
  real(r8) :: StalkRadius
  real(r8) :: TotLeafElevation
! begin_execution
  associate(                                                          &
    LeafElmsByLayerNode_brch  => plt_biom%LeafElmsByLayerNode_brch,   &
    LeafElmntNode_brch        => plt_biom%LeafElmntNode_brch,         &
    CanopyLeafCLyr_pft        => plt_biom%CanopyLeafCLyr_pft,         &
    StalkLiveBiomassC_brch        => plt_biom%StalkLiveBiomassC_brch,         &
    StalkStrutElms_brch       => plt_biom%StalkStrutElms_brch,        &
    PetoleProteinCNode_brch   => plt_biom%PetoleProteinCNode_brch,    &
    FracGroth2Node_pft        => plt_allom%FracGroth2Node_pft,        &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,     &
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft,   &
    KHiestGroLeafNode_brch    => plt_pheno%KHiestGroLeafNode_brch,    &
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft, &
    KLowestGroLeafNode_brch   => plt_pheno%KLowestGroLeafNode_brch,   &
    rLen2WidthLeaf_pft        => plt_morph%rLen2WidthLeaf_pft,        &
    SeedDepth_pft             => plt_morph%SeedDepth_pft,             &
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft,         &
    CanopyLeafArea_lpft       => plt_morph%CanopyLeafArea_lpft,       &
    PetoleLensNode_brch       => plt_morph%PetoleLensNode_brch,       &
    CanopyStalkArea_lbrch     => plt_morph%CanopyStalkArea_lbrch,     &
    BranchNumber_brch         => plt_morph%BranchNumber_brch,         &
    LiveInterNodeHight_brch   => plt_morph%LiveInterNodeHight_brch,   &
    CanopyHeightZ_col         => plt_morph%CanopyHeightZ_col,         &
    CanopyHeight_pft          => plt_morph%CanopyHeight_pft,          &
    HypoctoHeight_pft         => plt_morph%HypoctoHeight_pft,         &
    CLASS                     => plt_morph%CLASS,                     &
    LeafNodeArea_brch         => plt_morph%LeafNodeArea_brch,         &
    CanopyLeafAreaZ_pft       => plt_morph%CanopyLeafAreaZ_pft,       &
    PlantPopulation_pft       => plt_site%PlantPopulation_pft,        &
    ZERO                      => plt_site%ZERO,                       &
    SineLeafAngle             => plt_rad%SineLeafAngle                &
  )
!   ALLOCATION OF LEAF AREA TO CANOPY LAYERS
!
!   HypoctoHeight_pft=hypocotyledon height
!   SeedDepth_pft=seeding depth
!   LeafNodeArea_brch=node leaf area
!   PetoleLensNode_brch=petiole length
!
  KLowestGroLeafNode_brch(NB,NZ)=0

  IF(HypoctoHeight_pft(NZ).LE.SeedDepth_pft(NZ) .AND. LeafNodeArea_brch(0,MainBranchNum_pft(NZ),NZ).GT.0.0_r8)THEN
    LeafLength            = SQRT(1.0E+02_r8*LeafNodeArea_brch(0,MainBranchNum_pft(NZ),NZ)/PlantPopulation_pft(NZ))
    HypoctoHeight_pft(NZ) = LeafLength+PetoleLensNode_brch(0,MainBranchNum_pft(NZ),NZ) &
      +LiveInterNodeHight_brch(0,MainBranchNum_pft(NZ),NZ)

  ENDIF
!
! IF CANOPY HAS EMERGED
!
  IF(HypoctoHeight_pft(NZ).GT.SeedDepth_pft(NZ))THEN
    D540: DO K=0,MaxNodesPerBranch1
      DO  L=1,NumOfCanopyLayers1
        CanopyLeafArea_lpft(L,K,NB,NZ)                         = 0._r8
        LeafElmsByLayerNode_brch(1:NumPlantChemElms,L,K,NB,NZ) = 0._r8
      enddo
    ENDDO D540
    D535: DO L=1,NumOfCanopyLayers1
      CanopyStalkArea_lbrch(L,NB,NZ)=0._r8
    ENDDO D535
!
!   BRANCH HEIGHT
!
!   iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+petiole,2=none,3=between 1,2
!   iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!   KLeafNumHighestGrowing,KLowestGroLeafNode_brch=integer of highest,lowest leaf number currently growing
!   LiveInterNodeHight_brch=internode length
!   HeightBranchBase=branch base height
!
    IF(iPlantTurnoverPattern_pft(NZ).NE.0 .AND. is_plant_treelike(iPlantRootProfile_pft(NZ)))THEN
      IF(NB.NE.MainBranchNum_pft(NZ))THEN
        KLeafNumHighestGrowing=MAX(1,KHiestGroLeafNode_brch(MainBranchNum_pft(NZ),NZ)-MaxNodesPerBranch1+1)
        IF(BranchNumber_brch(NB,NZ).GE.KLeafNumHighestGrowing)THEN
          K                = pMOD(BranchNumber_brch(NB,NZ),MaxNodesPerBranch1)
          HeightBranchBase = LiveInterNodeHight_brch(K,MainBranchNum_pft(NZ),NZ)
        ELSE
          HeightBranchBase=0._r8
        ENDIF
      ELSE
        HeightBranchBase=0._r8
      ENDIF
    ELSE
      HeightBranchBase=0._r8
    ENDIF
    KMinGroingLeafNodeNum=MAX(0,KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+1)
!
!   FOR ALL LEAFED NODES
!
    D560: DO KK=KMinGroingLeafNodeNum,KHiestGroLeafNode_brch(NB,NZ)
      K=pMOD(KK,MaxNodesPerBranch1)
      !
      !     HEIGHT OF STALK INTERNODE + SHEATH OR PETIOLE
      !     AND LENGTH OF LEAF
      !
      !     HeightStalk=stalk height
      !     LiveInterNodeHight_brch=internode length
      !     HeightLeafNode=leaf node height
      !     LeafNodeArea_brch=leaf node area
      !     PP=plant population
      !     FracGroth2Node_pft=scales node number for perennial vegetation (e.g. trees)
      !     LeafLength=leaf length
!
      HeightStalk    = HeightBranchBase+LiveInterNodeHight_brch(K,NB,NZ)
      HeightLeafNode = HeightStalk+PetoleLensNode_brch(K,NB,NZ)
      LeafLength     = AZMAX1(SQRT(rLen2WidthLeaf_pft(NZ)*AZMAX1(LeafNodeArea_brch(K,NB,NZ)) &
        /(PlantPopulation_pft(NZ)*FracGroth2Node_pft(NZ))))
      TotLeafElevation = 0._r8
      !
      !   ALLOCATE FRACTIONS OF LEAF IN EACH INCLINATION CLASS
      !     FROM HIGHEST TO LOWEST TO CANOPY LAYER
      !
      !     LeafElevation=leaf elevation
      !     CLASS=leaf inclination class
      !     LeafLength=leaf length
      !     ZC,CanopyHeight_copy=canopy height
      !     HeightLeafBase,HeightLeafTip=height of leaf base,tip
      !     ZL=height to bottom of each canopy layer
      !     LNumHeightLeafBase,LNumHeightLeafTip=layer number of leaf base,tip
      !     FRACL=leaf fraction in each layer
!
      D555: DO N=NumOfLeafZenithSectors1,1,-1
        LeafElevation = SineLeafAngle(N)*CLASS(N,NZ)*LeafLength
        HeightLeafLow = AMIN1(CanopyHeight_copy(NZ)+0.01_r8-LeafElevation,HeightLeafNode+TotLeafElevation)
        HeightLeafTip = AMIN1(CanopyHeight_copy(NZ)+0.01_r8,HeightLeafLow+LeafElevation)
        LU            = 0
        LL            = 0
        D550: DO L = NumOfCanopyLayers1, 1, -1
          IF(LU.EQ.1 .AND. LL.EQ.1)exit
          IF((HeightLeafTip.GT.CanopyHeightZ_col(L-1) .OR. CanopyHeightZ_col(L-1).LE.ZERO) .AND. LU.EQ.0)THEN
            LNumHeightLeafTip = MAX(1,L)
            LU                = 1
          ENDIF
          IF((HeightLeafLow.GT.CanopyHeightZ_col(L-1) .OR. CanopyHeightZ_col(L-1).LE.ZERO) .AND. LL.EQ.0)THEN
            LNumHeightLeafBase = MAX(1,L)
            LL                 = 1
          ENDIF
        ENDDO D550

        D570: DO L=LNumHeightLeafBase,LNumHeightLeafTip
          IF(LNumHeightLeafTip.EQ.LNumHeightLeafBase)THEN
            FRACL=CLASS(N,NZ)
          ELSEIF(HeightLeafTip.GT.HeightLeafLow .AND. CanopyHeightZ_col(L).GT.HeightLeafLow)THEN
            FRACL=CLASS(N,NZ)*(AMIN1(HeightLeafTip,CanopyHeightZ_col(L)) &
              -AMAX1(HeightLeafLow,CanopyHeightZ_col(L-1)))/(HeightLeafTip-HeightLeafLow)
          ELSE
            FRACL=CLASS(N,NZ)
          ENDIF
          YLeafArea_brch=FRACL*LeafNodeArea_brch(K,NB,NZ)
          DO NE=1,NumPlantChemElms
            YLeafElmntNode_brch(NE)=FRACL*LeafElmntNode_brch(NE,K,NB,NZ)
          ENDDO
!
    !     ACCUMULATE LAYER LEAF AREAS, C, N AND P CONTENTS
    !
    !     CanopyLeafArea_lpft=leaf node area in canopy layer
    !     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
    !     CanopyLeafAreaZ_pft,CanopyLeafCLyr_pft=total leaf area,C in canopy layer
    !     LiveInterNodeHight_brch=internode length
    !
          CanopyLeafArea_lpft(L,K,NB,NZ)=CanopyLeafArea_lpft(L,K,NB,NZ)+YLeafArea_brch
          
          if(abs(LeafNodeArea_brch(K,NB,NZ))>1.e10)stop
          DO NE=1,NumPlantChemElms
            LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)=LeafElmsByLayerNode_brch(NE,L,K,NB,NZ)+YLeafElmntNode_brch(NE)
          ENDDO
          CanopyLeafAreaZ_pft(L,NZ) = CanopyLeafAreaZ_pft(L,NZ)+YLeafArea_brch
          CanopyLeafCLyr_pft(L,NZ)  = CanopyLeafCLyr_pft(L,NZ)+YLeafElmntNode_brch(ielmc)
        ENDDO D570
        TotLeafElevation=TotLeafElevation+LeafElevation
        CanopyHeight_pft(NZ)=AMAX1(CanopyHeight_pft(NZ),HeightLeafTip)
      ENDDO D555

      IF(PetoleProteinCNode_brch(K,NB,NZ).GT.0.0_r8)THEN
        IF(KLowestGroLeafNode_brch(NB,NZ).EQ.0)KLowestGroLeafNode_brch(NB,NZ)=MIN(KK,KHiestGroLeafNode_brch(NB,NZ))
      ENDIF
    ENDDO D560

    IF(KLowestGroLeafNode_brch(NB,NZ).EQ.0)KLowestGroLeafNode_brch(NB,NZ)=KHiestGroLeafNode_brch(NB,NZ)
    K1=pMOD(KHiestGroLeafNode_brch(NB,NZ),MaxNodesPerBranch1)
    K2=pMOD(KHiestGroLeafNode_brch(NB,NZ)-1,MaxNodesPerBranch1)

    IF(isclose(LiveInterNodeHight_brch(K1,NB,NZ),0._r8))THEN
      LiveInterNodeHight_brch(K1,NB,NZ)=LiveInterNodeHight_brch(K2,NB,NZ)
    ENDIF
    HeightLeafBase=HeightBranchBase+AZMAX1(LiveInterNodeHight_brch(K1,NB,NZ))
!
  !     ALLOCATE STALK SURFACE AREA TO CANOPY LAYERS
  !
  !     LiveInterNodeHight_brch=internode length
  !     HeightLeafBase=leaf base height
  !     ZL=height to bottom of each canopy layer
  !     LNumHeightBranchBase,LNumHeightBranchTip=layer number of branch base,tip
  !     WTSTKB,StalkSurfArea=branch stalk mass,surface area
  !     FSTK=fraction of stalk area contributing to water,heat flow
  !     StalkMassDensity,SpecStalkVolume=stalk density (Mg m-3),specific volume (m3 g-1)
  !     StalkLiveBiomassC_brch=stalk sapwood mass
  !     FRACL=stalk fraction in each layer
  !     CanopyStalkArea_lbrch=total branch stalk surface area in each layer
  !
    IF(LiveInterNodeHight_brch(K1,NB,NZ).GT.0.0_r8)THEN
      LU=0
      LL=0
      D545: DO L=NumOfCanopyLayers1,1,-1
        IF(LU.EQ.1 .AND. LL.EQ.1)exit
        IF((HeightLeafBase.GT.CanopyHeightZ_col(L-1) .OR. CanopyHeightZ_col(L-1).LE.ZERO) .AND. LU.EQ.0)THEN
          LNumHeightBranchTip = MAX(1,L)
          LU                  = 1
        ENDIF
        IF((HeightBranchBase.GT.CanopyHeightZ_col(L-1).OR.CanopyHeightZ_col(L-1).LE.ZERO) .AND. LL.EQ.0)THEN
          LNumHeightBranchBase = MAX(1,L)
          LL                   = 1
        ENDIF
      ENDDO D545
      
      StalkRadius=SQRT(SpecStalkVolume*(AZMAX1(StalkStrutElms_brch(ielmc,NB,NZ))/PlantPopulation_pft(NZ)) &
        /(PICON*LiveInterNodeHight_brch(K1,NB,NZ)))
      StalkSurfArea = PICON*LiveInterNodeHight_brch(K1,NB,NZ)*PlantPopulation_pft(NZ)*StalkRadius

      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
        StalkLiveBiomassC_brch(NB,NZ)=StalkStrutElms_brch(ielmc,NB,NZ)
      ELSE        
        RadiusInnerStalk4Transp   = AMIN1(ZSTX,FSTK*StalkRadius)
        StalkSectionArea          = PICON*StalkRadius*(StalkRadius-RadiusInnerStalk4Transp)
        StalkLiveBiomassC_brch(NB,NZ) = StalkSectionArea*LiveInterNodeHight_brch(K1,NB,NZ)*PlantPopulation_pft(NZ)/SpecStalkVolume
      ENDIF

      D445: DO L=LNumHeightBranchBase,LNumHeightBranchTip
        IF(HeightLeafBase.GT.HeightBranchBase)THEN
          IF(HeightLeafBase.GT.CanopyHeightZ_col(L-1))THEN
            FRACL=(AMIN1(HeightLeafBase,CanopyHeightZ_col(L))-AMAX1(HeightBranchBase,CanopyHeightZ_col(L-1))) &
              /(HeightLeafBase-HeightBranchBase)
          ELSE
            FRACL=0._r8
          ENDIF
        ELSE
          FRACL=1.0_r8
        ENDIF
        CanopyStalkArea_lbrch(L,NB,NZ)=FRACL*StalkSurfArea
      ENDDO D445
    ELSE
      StalkLiveBiomassC_brch(NB,NZ)=0._r8
      D450: DO L=1,NumOfCanopyLayers1
        CanopyStalkArea_lbrch(L,NB,NZ)=0._r8
      ENDDO D450
    ENDIF
  ELSE
    StalkLiveBiomassC_brch(NB,NZ)=0._r8
    D455: DO L=1,NumOfCanopyLayers1
      CanopyStalkArea_lbrch(L,NB,NZ)=0._r8
    ENDDO D455
  ENDIF
  end associate
  end subroutine AllocateLeaf2CanopyLayers
!------------------------------------------------------------------------------------------
  subroutine LeafClassAllocation(NB,NZ)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8) :: dangle
  integer :: L,K,N
  ! begin_execution
  associate(                                                            &
    SineBranchAngle_pft         =>  plt_morph%SineBranchAngle_pft     , &
    LeafNodeArea_brch           =>  plt_morph%LeafNodeArea_brch       , &
    CanopyStalkArea_lbrch       =>  plt_morph%CanopyStalkArea_lbrch   , &
    CLASS                       =>  plt_morph%CLASS                   , &
    CanopyLeafArea_lpft         =>  plt_morph%CanopyLeafArea_lpft     , &
    StemAreaZsec_brch           =>  plt_morph%StemAreaZsec_brch       , &
    MainBranchNum_pft           =>  plt_morph%MainBranchNum_pft       , &
    LeafAreaZsec_brch           =>  plt_morph%LeafAreaZsec_brch         &
  )
  D900: DO K=1,MaxNodesPerBranch1
    DO  L=1,NumOfCanopyLayers1
      DO  N=1,NumOfLeafZenithSectors1
        LeafAreaZsec_brch(N,L,K,NB,NZ)=0._r8
      enddo
    enddo
  ENDDO D900
! ARLFXB=0._r8
! ARLFXL=0._r8
! SURFXX=0._r8
  D500: DO K=1,MaxNodesPerBranch1
!     ARLFXB=ARLFXB+LeafNodeArea_brch(K,NB,NZ)
    IF(LeafNodeArea_brch(K,NB,NZ).GT.0.0_r8)THEN
      D700: DO L=NumOfCanopyLayers1,1,-1
!       ARLFXL=ARLFXL+CanopyLeafArea_lpft(L,K,NB,NZ)
        D800: DO N=1,NumOfLeafZenithSectors1
          LeafAreaZsec_brch(N,L,K,NB,NZ)=AZMAX1(CLASS(N,NZ)*0.25_r8*CanopyLeafArea_lpft(L,K,NB,NZ))
  !       SURFXX=SURFXX+LeafAreaZsec_brch(N,L,K,NB,NZ)
        ENDDO D800
      ENDDO D700
    ENDIF
  ENDDO D500
!
! ALLOCATE STALK AREA TO INCLINATION CLASSES ACCORDING TO
! BRANCH ANGLE ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
!
! StemAreaZsec_brch=stalk surface area in canopy layer
! SineBranchAngle_pft=stem angle from horizontal
! CanopyStalkArea_lbrch=total branch stalk surface area in each layer
!
  D910: DO L=1,NumOfCanopyLayers1
    DO  N=1,NumOfLeafZenithSectors1
      StemAreaZsec_brch(N,L,NB,NZ)=0._r8
    enddo
  ENDDO D910

  IF(NB.EQ.MainBranchNum_pft(NZ))THEN
    N=NumOfLeafZenithSectors1
  ELSE
    dangle=PICON2h/real(NumOfLeafZenithSectors1,r8)
    N=MIN(NumOfLeafZenithSectors1,INT(ASIN(SineBranchAngle_pft(NZ))/dangle)+1)
  ENDIF
  D710: DO L=NumOfCanopyLayers1,1,-1
    StemAreaZsec_brch(N,L,NB,NZ)=CanopyStalkArea_lbrch(L,NB,NZ)/real(NumOfLeafZenithSectors1,r8)
  ENDDO D710
  end associate
  end subroutine LeafClassAllocation

!------------------------------------------------------------------------------------------

  subroutine GrainFilling(I,J,NB,NZ,Growth_brch,GROSTKC)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: Growth_brch(NumPlantChemElms),GROSTKC
  real(r8) :: ZPGRP,ZPGRN,ZNPGP,ZNPGN
  real(r8) :: MaxChemElmRsrv2Grain,ChemElmRsrv2Grain(NumPlantChemElms)
  real(r8) :: FGRNX
  real(r8) :: GRMXB
  real(r8) :: GROLM
  real(r8) :: GROLC
  real(r8) :: SeedSET
  integer :: NE
! begin_execution
  associate(                                                                &
    TdegCCanopy_pft           => plt_ew%TdegCCanopy_pft,              &
    LeafPetoNonstElmConc_brch    => plt_biom%LeafPetoNonstElmConc_brch,     &
    GrainStrutElms_brch          => plt_biom%GrainStrutElms_brch,           &
    StalkRsrvElms_brch           => plt_biom%StalkRsrvElms_brch,            &
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft,                &
    GrainSeedBiomCMean_brch      => plt_allom%GrainSeedBiomCMean_brch,      &
    CNGR                         => plt_allom%CNGR,                         &
    CPGR                         => plt_allom%CPGR,                         &
    fTgrowRootP_vr               => plt_pheno%fTgrowRootP_vr,               &
    fTCanopyGroth_pft            => plt_pheno%fTCanopyGroth_pft,            &
    iPlantPhenolType_pft         => plt_pheno%iPlantPhenolType_pft,         &
    Hours4LeafOff_brch           => plt_pheno%Hours4LeafOff_brch,           &
    HourFailGrainFill_brch       => plt_pheno%HourFailGrainFill_brch,       &
    HourReq4LeafOff_brch         => plt_pheno%HourReq4LeafOff_brch,         &
    GrainFillRate25C_pft         => plt_pheno%GrainFillRate25C_pft,         &
    SeedTempSens_pft             => plt_pheno%SeedTempSens_pft,             &
    iPlantPhenolPattern_pft      => plt_pheno%iPlantPhenolPattern_pft,      &
    HighTempLimitSeed_pft        => plt_pheno%HighTempLimitSeed_pft,        &
    dReproNodeNumNormByMatG_brch => plt_pheno%dReproNodeNumNormByMatG_brch, &
    iPlantCalendar_brch          => plt_pheno%iPlantCalendar_brch,          &
    TCChill4Seed_pft             => plt_pheno%TCChill4Seed_pft,             &
    PlantPopulation_pft          => plt_site%PlantPopulation_pft,           &
    MaxSeedNumPerSite_pft        => plt_morph%MaxSeedNumPerSite_pft,        &
    MaxSeedCMass_pft             => plt_morph%MaxSeedCMass_pft,             &
    MaxPotentSeedNumber_pft      => plt_morph%MaxPotentSeedNumber_pft,      &
    NGTopRootLayer_pft           => plt_morph%NGTopRootLayer_pft,           &
    PotentialSeedSites_brch      => plt_morph%PotentialSeedSites_brch,      &
    iPlantGrainType_pft          => plt_morph%iPlantGrainType_pft,          &
    MainBranchNum_pft            => plt_morph%MainBranchNum_pft,            &
    SeedNumSet_brch              => plt_morph%SeedNumSet_brch               &
  )
!
!   SET MAXIMUM GRAIN NUMBER FROM SHOOT MASS BEFORE ANTHESIS
!
!   iPlantCalendar_brch(ipltcal_Jointing,=start of stem elongation and setting max seed number
!   iPlantCalendar_brch(ipltcal_Anthesis,=start of anthesis and setting final seed number
!   PotentialSeedSites_brch=potential number of seed set sites
!   MaxPotentSeedNumber_pft=maximum potential seed number from PFT file
!     GROSTKC=stalk growth rate
!
  IF(iPlantCalendar_brch(ipltcal_Jointing,NB,NZ).NE.0 .AND. iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
    PotentialSeedSites_brch(NB,NZ)=PotentialSeedSites_brch(NB,NZ)+MaxPotentSeedNumber_pft(NZ)*AZMAX1(GROSTKC)
  ENDIF
!
!   SET FINAL GRAIN NUMBER FROM C,N,P NON-STRUCTURAL POOLS AFTER ANTHESIS
!
!   iPlantCalendar_brch(ipltcal_Anthesis,=start of anthesis and setting final seed number
!   iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
!   iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date setting for final seed number
!   iPlantCalendar_brch(ipltcal_SetSeedMass,=end of setting max seed size
!   SET=seed set limited by nonstructural C,N,P
!   CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!   TdegCCanopy_pft=canopy temperature
!   TCChill4Seed_pft=chilling temperature for CO2 fixation, seed loss (oC)
!   HTC=high temperature threshold for grain number loss
!   FGRNX=loss of seed set
!   SeedTempSens_pft=sensitivity to TdegCCanopy_pft> HTC,TdegCCanopy_pft< TCChill4Seed_pftfrom startq.f (seeds oC-1)
!   SeedNumSet_brch=seed set number
!   PotentialSeedSites_brch=potential number of seed set sites
!   MaxSeedNumPerSite_pft=maximum seed number per MaxPotentSeedNumber_pft from PFT file
!   dReproNodeNumNormByMatG_brch=change in reproductive node number normalized for maturity group
!
  IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).NE.0.AND.iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ).EQ.0)THEN
    SeedSET=AMIN1(LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ)+SETC) &
      ,LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ)+SETN) &
      ,LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ)+SETP))

    IF(TdegCCanopy_pft(NZ).LT.TCChill4Seed_pft(NZ))THEN
      IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
        FGRNX=SeedTempSens_pft(NZ)*(TCChill4Seed_pft(NZ)-TdegCCanopy_pft(NZ))
      ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN
        FGRNX=SeedTempSens_pft(NZ)*(TCChill4Seed_pft(NZ)-TdegCCanopy_pft(NZ))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSEIF(TdegCCanopy_pft(NZ).GT.HighTempLimitSeed_pft(NZ))THEN
      IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
        FGRNX=SeedTempSens_pft(NZ)*(TdegCCanopy_pft(NZ)-HighTempLimitSeed_pft(NZ))
      ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN
        FGRNX=SeedTempSens_pft(NZ)*(TdegCCanopy_pft(NZ)-HighTempLimitSeed_pft(NZ))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSE
      FGRNX=0._r8
    ENDIF
    IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).NE.0.AND.iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN
!     PotentialSeedSites_brch(NB,NZ)=PotentialSeedSites_brch(NB,NZ)*FGRNX
      SeedNumSet_brch(NB,NZ)=AMIN1(MaxSeedNumPerSite_pft(NZ)*PotentialSeedSites_brch(NB,NZ) &
        ,SeedNumSet_brch(NB,NZ)+MaxSeedNumPerSite_pft(NZ)*PotentialSeedSites_brch(NB,NZ) &
        *SeedSET*dReproNodeNumNormByMatG_brch(NB,NZ)-FGRNX*SeedNumSet_brch(NB,NZ))
    ENDIF
!

!     SET MAXIMUM GRAIN SIZE FROM C,N,P NON-STRUCTURAL POOLS AFTER ANTHESIS
!
!     GRMX=maximum individual seed size from PFT file (g)
!     dReproNodeNumNormByMatG_brch=change in reproductive node number normalized for maturity group
!     GrainSeedBiomCMean_brch=individual seed size
!     SET=seed set limited by nonstructural C,N,P
!
    IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0 &
      .AND.iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ).EQ.0)THEN
      GRMXB=MaxSeedCMass_pft(NZ)
      GrainSeedBiomCMean_brch(NB,NZ)=AMIN1(MaxSeedCMass_pft(NZ),GrainSeedBiomCMean_brch(NB,NZ) &
        +GRMXB*AMAX1(0.50_r8,SeedSET**0.25_r8)*dReproNodeNumNormByMatG_brch(NB,NZ))
    ENDIF
  ENDIF
!
!   GRAIN FILL BY TRANSLOCATION FROM STALK RESERVES
!   UNTIL GRAIN SINK (=FINAL GRAIN NUMBER X MAXIMUM
!   GRAIN SIZE) IS FILLED OR RESERVES ARE EXHAUSTED
!
!   iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
!   WTGRB=total seed C mass
!   GrainSeedBiomCMean_brch=individual seed size
!   SeedNumSet_brch=seed set number
!   GROLM=maximum grain fill rate
!   GrainFillRate25C_pft=grain filling rate at 25 oC from PFT file
!   fTCanopyGroth_pft=temperature function for canopy growth
!   fTgrowRootP_vr=temperature function for root growth
!
  IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0)THEN
    IF(GrainStrutElms_brch(ielmc,NB,NZ).GE.GrainSeedBiomCMean_brch(NB,NZ)*SeedNumSet_brch(NB,NZ))THEN
      GROLM=0._r8
    ELSEIF(iPlantGrainType_pft(NZ).EQ.igraintyp_abvgrnd)THEN
      GROLM=AZMAX1(GrainFillRate25C_pft(NZ)*SeedNumSet_brch(NB,NZ)*SQRT(fTCanopyGroth_pft(NZ)))
    ELSE
      GROLM=AZMAX1(GrainFillRate25C_pft(NZ)*SeedNumSet_brch(NB,NZ)*SQRT(fTgrowRootP_vr(NGTopRootLayer_pft(NZ),NZ)))
    ENDIF
!
!     GRAIN FILL RATE MAY BE CONSTRAINED BY HIGH GRAIN C:N OR C:P
!
!     WTGRB,WTGRBN,WTGRBP=total seed C,N,P mass
!     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
!     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
!     GROLM,GROLC=maximum,actual grain fill rate
!     MaxChemElmRsrv2Grain,ChemElmRsrv2Grain(ielmc)=maximum,actual C translocation rate from reserve to grain
!
    IF(GrainStrutElms_brch(ielmn,NB,NZ).LT.ZPGRM*CNGR(NZ)*GrainStrutElms_brch(ielmc,NB,NZ) &
      .OR.GrainStrutElms_brch(ielmp,NB,NZ).LT.ZPGRM*CPGR(NZ)*GrainStrutElms_brch(ielmc,NB,NZ))THEN
      GROLC=0._r8
    ELSE
      GROLC=GROLM
    ENDIF
    MaxChemElmRsrv2Grain=AMIN1(GROLM,StalkRsrvElms_brch(ielmc,NB,NZ))
    ChemElmRsrv2Grain(ielmc)=AMIN1(GROLC,StalkRsrvElms_brch(ielmc,NB,NZ))
!
!     GRAIN N OR P FILL RATE MAY BE LIMITED BY C:N OR C:P RATIOS
!     OF STALK RESERVES
!
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     ZNPGN,ZNPGP=effect of reserve N:C,P:C on grain fill N:C,P:C
!     SETN,SETP=Km for nonstructural N,P concn on seed set (g g-1)
!     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
!     ZPGRD=1.0_r8-ZPGRM
!     ZPGRN,ZPGRP=N:C,P:C ratios during grain fill
!     MaxChemElmRsrv2Grain,ChemElmRsrv2Grain(ielmc)=maximum,actual C translocation rate from reserve to grain
!     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
!     ChemElmRsrv2Grain(ielmn),ChemElmRsrv2Grain(ielmp)=N,P translocation rate from reserve to grain
!
    IF(StalkRsrvElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
      ZNPGN=StalkRsrvElms_brch(ielmn,NB,NZ)/(StalkRsrvElms_brch(ielmn,NB,NZ) &
        +SETN*StalkRsrvElms_brch(ielmc,NB,NZ))
      ZNPGP=StalkRsrvElms_brch(ielmp,NB,NZ)/(StalkRsrvElms_brch(ielmp,NB,NZ) &
        +SETP*StalkRsrvElms_brch(ielmc,NB,NZ))
      ZPGRN                    = ZPGRM+ZPGRD*AZMAX1(AMIN1(1.0_r8,ZNPGN))
      ZPGRP                    = ZPGRM+ZPGRD*AZMAX1(AMIN1(1.0_r8,ZNPGP))
      ChemElmRsrv2Grain(ielmn) = AMIN1(MaxChemElmRsrv2Grain*CNGR(NZ) &
        ,AZMAX1(StalkRsrvElms_brch(ielmn,NB,NZ)*ZPGRN) &
        ,(GrainStrutElms_brch(ielmc,NB,NZ)+ChemElmRsrv2Grain(ielmc))*CNGR(NZ) &
        -GrainStrutElms_brch(ielmn,NB,NZ))
      ChemElmRsrv2Grain(ielmp)=AMIN1(MaxChemElmRsrv2Grain*CPGR(NZ) &
        ,AZMAX1(StalkRsrvElms_brch(ielmp,NB,NZ)*ZPGRP) &
        ,(GrainStrutElms_brch(ielmc,NB,NZ)+ChemElmRsrv2Grain(ielmc))*CPGR(NZ) &
        -GrainStrutElms_brch(ielmp,NB,NZ))
    ELSE
      ChemElmRsrv2Grain(ielmn)=0._r8
      ChemElmRsrv2Grain(ielmp)=0._r8
    ENDIF
!
!     TRANSLOCATE C,N,P FROM STALK RESERVES TO GRAIN
!
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     GROGR=grain growth rate
!     ChemElmRsrv2Grain(ielmc),ChemElmRsrv2Grain(ielmn),ChemElmRsrv2Grain(ielmp)=C,N,P translocation rate from reserve to grain
!
    DO NE=1,NumPlantChemElms
      StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+Growth_brch(NE) &
        -ChemElmRsrv2Grain(NE)
      GrainStrutElms_brch(NE,NB,NZ)=GrainStrutElms_brch(NE,NB,NZ)+ChemElmRsrv2Grain(NE)
    ENDDO
  ELSE
    ChemElmRsrv2Grain(1:NumPlantChemElms)=0._r8
  ENDIF
!
!   SET DATE OF PHYSIOLOGICAL MATURITY WHEN GRAIN FILL
!   HAS STOPPED FOR SET PERIOD OF TIME
!
!   iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date setting for final seed number
!   ChemElmRsrv2Grain(ielmc)=C translocation rate from reserve to grain
!   PP=PFT population
!   HourFailGrainFill_brch=number of hours with no grain fill
!   Hours4PhyslMature=number of hours with no grain filling until physl maturity
!   iPlantCalendar_brch(ipltcal_EndSeedFill,=date of physiological maturity
!
  IF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).NE.0)THEN
    IF(ChemElmRsrv2Grain(ielmc).LE.1.0E-09_r8*PlantPopulation_pft(NZ))THEN
      HourFailGrainFill_brch(NB,NZ)=HourFailGrainFill_brch(NB,NZ)+1.0
    ELSE
      HourFailGrainFill_brch(NB,NZ)=0._r8
    ENDIF
    IF(HourFailGrainFill_brch(NB,NZ).GE.Hours4PhyslMature)THEN
      IF(iPlantCalendar_brch(ipltcal_EndSeedFill,NB,NZ).EQ.0)THEN
        iPlantCalendar_brch(ipltcal_EndSeedFill,NB,NZ)=I
      ENDIF
    ENDIF
!
!     TERMINATE ANNUALS AFTER GRAIN FILL
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     HourFailGrainFill_brch=number of hours with no grain fill
!     Hours4PhyslMature=number of hours with no grain filling until physiological maturity
!     Hours4SenesAftMature=number of hours after physiol maturity required for senescence
!     Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).NE.0)THEN
      IF(HourFailGrainFill_brch(NB,NZ).GT.Hours4PhyslMature+Hours4SenesAftMature(iPlantPhenolType_pft(NZ)))THEN
        Hours4LeafOff_brch(NB,NZ)=HourReq4LeafOff_brch(NB,NZ)+0.5_r8
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine GrainFilling
!------------------------------------------------------------------------------------------
  subroutine ResetNonAnnualBranch(I,J,NB,NZ)      
!  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NB,NZ
  integer :: K,M,NE

  associate(                                                                          &
    NodeNum2InitFloral_brch           => plt_morph%NodeNum2InitFloral_brch,           &
    NodeNumberAtAnthesis_brch         => plt_morph%NodeNumberAtAnthesis_brch,         &
    LeafNumberAtFloralInit_brch       => plt_pheno%LeafNumberAtFloralInit_brch,       &
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch,     &
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch, &
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft,                   &    
    icwood                            => pltpar%icwood,                               &  
    istalk                            => pltpar%istalk,                               &      
    ifoliar                           => pltpar%ifoliar,                              &
    inonfoliar                        => pltpar%inonfoliar,                           &        
    k_fine_litr                       => pltpar%k_fine_litr,                          &    
    k_woody_litr                      => pltpar%k_woody_litr,                         &
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft,                 &        
    MatureGroup_brch                  => plt_pheno%MatureGroup_brch,                  &
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch,               &
    HoursTooLowPsiCan_pft             => plt_pheno%HoursTooLowPsiCan_pft,             &
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft,           &    
    doPlantLeafOut_brch               => plt_pheno%doPlantLeafOut_brch,               &           
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch,                &
    iPlantTurnoverPattern_pft         => plt_pheno%iPlantTurnoverPattern_pft,         &        
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch,                 &     
    NumOfLeaves_brch                  => plt_morph%NumOfLeaves_brch,                  &    
    KLeafNumber_brch                  => plt_morph%KLeafNumber_brch,                  &    
    KHiestGroLeafNode_brch            => plt_pheno%KHiestGroLeafNode_brch,            &    
    LitrfalStrutElms_pvr              => plt_bgcr%LitrfalStrutElms_pvr,               &    
    HourFailGrainFill_brch            => plt_pheno%HourFailGrainFill_brch,            &    
    ElmAllocmat4Litr                  => plt_soilchem%ElmAllocmat4Litr,               &    
    LeafStrutElms_brch                => plt_biom%LeafStrutElms_brch,                 &    
    FracShootStalkElmAlloc2Litr       => plt_allom%FracShootStalkElmAlloc2Litr,       &    
    PetoleStrutElms_brch              => plt_biom%PetoleStrutElms_brch,               &    
    FracShootLeafElmAlloc2Litr        => plt_allom%FracShootLeafElmAlloc2Litr,        &   
    LeafAreaLive_brch                 => plt_morph%LeafAreaLive_brch,                 &   
    LeafElmntNode_brch                => plt_biom%LeafElmntNode_brch,                 &      
    LeafNodeArea_brch                 => plt_morph%LeafNodeArea_brch,                 &    
    PetioleElmntNode_brch             => plt_biom%PetioleElmntNode_brch,              &  
    PetoleProteinCNode_brch           => plt_biom%PetoleProteinCNode_brch,            &      
    PetoleLensNode_brch               => plt_morph%PetoleLensNode_brch,               &    
    LeafProteinCNode_brch             => plt_biom%LeafProteinCNode_brch,              &    
    HourReq4LeafOut_brch              => plt_pheno%HourReq4LeafOut_brch,              &
    ShootNodeNumAtPlanting_pft        => plt_morph%ShootNodeNumAtPlanting_pft,        &    
    HuskStrutElms_brch                => plt_biom%HuskStrutElms_brch,                 &        
    EarStrutElms_brch                 => plt_biom%EarStrutElms_brch,                  &    
    GrainStrutElms_brch               => plt_biom%GrainStrutElms_brch,                &    
    PotentialSeedSites_brch           => plt_morph%PotentialSeedSites_brch,           &    
    doInitLeafOut_brch                => plt_pheno%doInitLeafOut_brch       ,         &
    GrainSeedBiomCMean_brch           => plt_allom%GrainSeedBiomCMean_brch,           &    
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft,             &    
    StalkStrutElms_brch               => plt_biom%StalkStrutElms_brch,                &  
    SenecStalkStrutElms_brch          => plt_biom%SenecStalkStrutElms_brch,           &      
    InternodeStrutElms_brch           => plt_biom%InternodeStrutElms_brch,            &    
    LiveInterNodeHight_brch           => plt_morph%LiveInterNodeHight_brch,           &    
    InternodeHeightDead_brch         => plt_morph%InternodeHeightDead_brch,         &    
    doPlantLeaveOff_brch              => plt_pheno%doPlantLeaveOff_brch,              &    
    Prep4Literfall_brch               => plt_pheno%Prep4Literfall_brch,               &  
    Hours4LiterfalAftMature_brch      => plt_pheno%Hours4LiterfalAftMature_brch,      &      
    SeedNumSet_brch                   => plt_morph%SeedNumSet_brch                    &    
  )    
!    GROUP,MatureGroup_pft=node number required for floral initiation
!    NodeNum2InitFloral_brch=node number at floral initiation
!    NodeNumberAtAnthesis_brch=node number at flowering
!    VSTGX=leaf number on date of floral initiation
!    TotalNodeNumNormByMatgrp_brch=total change in vegve node number normalized for maturity group
!    TotReproNodeNumNormByMatrgrp_brch=total change in reprve node number normalized for maturity group
!    iPlantCalendar_brch(ipltcal_Emerge)=emergence date
!
  IF((doPlantLeafOut_brch(NB,NZ).EQ.iEnable .AND. iPlantPhenolPattern_pft(NZ).NE.iplt_annual) &
    .AND. (Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ)))THEN

!        IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
!          MatureGroup_brch(NB,NZ)=AZMAX1(MatureGroup_pft(NZ)-BranchNumber_brch(NB,NZ))
!        ELSE
      MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
!        ENDIF
    NodeNum2InitFloral_brch(NB,NZ)               = ShootNodeNum_brch(NB,NZ)
    NodeNumberAtAnthesis_brch(NB,NZ)             = 0._r8
    LeafNumberAtFloralInit_brch(NB,NZ)           = 0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ)         = 0._r8
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ)     = 0._r8
    iPlantCalendar_brch(ipltcal_Emerge,NB,NZ)    = I
    iPlantCalendar_brch(2:NumGrowthStages,NB,NZ) = 0

    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      HoursTooLowPsiCan_pft(NZ)=0._r8
    ENDIF
!
!   SPRING LEAF AND SHEATH RESET
!
!   doInitLeafOut_brch,doPlantLeafOut_brch=flag for initializing,enabling leafout
!   Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
!   PSTG=node number
!   NumOfLeaves_brch=number of leaves appeared
!     KHiestGroLeafNode_brch=integer of most recent leaf number currently growing
!     HourFailGrainFill_brch=number of hours with no grain fill
!     CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!     WT*,WT*N,WT*P=branch organ C,N,P mass
!     WG,WG*N,WG*P=node organ C,N,P mass
!     organ key:LF=leaf,SHE=petiole,STK=stalk,RSV=reserve
!     HSK=husk,EAR=ear,GR=grain,SHT=shoot
!     LeafAreaLive_brch,LeafNodeArea_brch=branch,node leaf area
!     InternodeStrutElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
!     InternodeHeightDead_brch,LiveInterNodeHight_brch=stalk height,stalk internode length
!     SeedNumSet_brch=seed set number
!     PotentialSeedSites_brch=potential number of seed set sites
!     GrainSeedBiomCMean_brch=individual seed size
!
    IF(doPlantLeafOut_brch(NB,NZ).EQ.iEnable .AND. iPlantPhenolPattern_pft(NZ).NE.iplt_annual &
      .AND. Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
      IF(iPlantTurnoverPattern_pft(NZ).EQ.0)THEN
        ShootNodeNum_brch(NB,NZ)      = ShootNodeNumAtPlanting_pft(NZ)
        NumOfLeaves_brch(NB,NZ)       = 0._r8
        KLeafNumber_brch(NB,NZ)       = 1
        KHiestGroLeafNode_brch(NB,NZ) = 1
        HourFailGrainFill_brch(NB,NZ) = 0._r8
        D5330: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,icwood,M,NZ)*AZMAX1(LeafStrutElms_brch(NE,NB,NZ))*FracShootStalkElmAlloc2Litr(NE,k_woody_litr) &
              +ElmAllocmat4Litr(NE,icwood,M,NZ)*AZMAX1(PetoleStrutElms_brch(NE,NB,NZ))*FracShootLeafElmAlloc2Litr(NE,k_woody_litr)

            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
              +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*AZMAX1(LeafStrutElms_brch(NE,NB,NZ))*FracShootStalkElmAlloc2Litr(NE,k_fine_litr) &
              +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(PetoleStrutElms_brch(NE,NB,NZ))*FracShootLeafElmAlloc2Litr(NE,k_fine_litr)
          ENDDO  
        ENDDO D5330
        LeafAreaLive_brch(NB,NZ)                       = 0._r8
        LeafStrutElms_brch(1:NumPlantChemElms,NB,NZ)   = 0._r8
        PetoleStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
        D5335: DO K=0,MaxNodesPerBranch1
          LeafNodeArea_brch(K,NB,NZ)                        = 0._r8
          PetoleLensNode_brch(K,NB,NZ)                      = 0._r8
          LeafProteinCNode_brch(K,NB,NZ)                    = 0._r8
          LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)    = 0._r8
          PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8
          PetoleProteinCNode_brch(K,NB,NZ)                  = 0._r8
        ENDDO D5335
      ENDIF
    ENDIF
!
!     RESIDUAL STALKS BECOME LitrFall IN GRASSES, SHRUBS AT
!     START OF SEASON
!
    IF((doPlantLeafOut_brch(NB,NZ).EQ.iEnable.AND.iPlantPhenolPattern_pft(NZ).NE.iplt_annual)&
      .AND.Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN

      D6245: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,inonfoliar,M,NZ) &
            *AZMAX1(HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ)+GrainStrutElms_brch(NE,NB,NZ))
        ENDDO
      ENDDO D6245
        
      HuskStrutElms_brch(1:NumPlantChemElms,NB,NZ)  = 0._r8
      EarStrutElms_brch(1:NumPlantChemElms,NB,NZ)   = 0._r8
      GrainStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
      
      PotentialSeedSites_brch(NB,NZ) = 0._r8
      SeedNumSet_brch(NB,NZ)         = 0._r8
      GrainSeedBiomCMean_brch(NB,NZ) = 0._r8
      IF(iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
        D6345: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)&
              +ElmAllocmat4Litr(NE,istalk,M,NZ)*AZMAX1(StalkStrutElms_brch(NE,NB,NZ))
          ENDDO
        ENDDO D6345
        StalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)      = 0._r8
        SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ) = 0._r8
        DO K=0,MaxNodesPerBranch1
          DO NE=1,NumPlantChemElms
            InternodeStrutElms_brch(NE,K,NB,NZ)=0._r8
          ENDDO
        ENDDO
        D6340: DO K=0,MaxNodesPerBranch1
          LiveInterNodeHight_brch(K,NB,NZ)   = 0._r8
          InternodeHeightDead_brch(K,NB,NZ) = 0._r8
        ENDDO D6340
      ENDIF
    ENDIF
  ENDIF
!
!   SPRING OR FALL FLAG RESET
!

  IF(doPlantLeafOut_brch(NB,NZ).EQ.iEnable .AND. Hours4Leafout_brch(NB,NZ) &
    .GE.HourReq4LeafOut_brch(NB,NZ))THEN
    doPlantLeafOut_brch(NB,NZ)          = iDisable
    doPlantLeaveOff_brch(NB,NZ)         = iEnable
    Prep4Literfall_brch(NB,NZ)          = ifalse
    Hours4LiterfalAftMature_brch(NB,NZ) = 0
    !doing leave off    
  ELSE
    doPlantLeafOut_brch(NB,NZ)          = iEnable
    doPlantLeaveOff_brch(NB,NZ)         = iDisable
    Prep4Literfall_brch(NB,NZ)          = itrue
    Hours4LiterfalAftMature_brch(NB,NZ) = 0
    doInitLeafOut_brch(NB,NZ)           = iEnable
  ENDIF
  end associate
  end subroutine ResetNonAnnualBranch
!------------------------------------------------------------------------------------------  
  subroutine ResetBranchPhenology(I,J,NB,NZ)

  implicit none
  integer, intent(in) :: I,J,nb,nz
  integer :: K,M,NE
  real(r8) :: FSNR1
  logical :: LeafOffTest,LeafOutTest
! begin_execution
  associate(                                                                          &
    istalk                            => pltpar%istalk,                               &
    iroot                             => pltpar%iroot,                                &
    k_fine_litr                       => pltpar%k_fine_litr,                          &    
    inonfoliar                        => pltpar%inonfoliar,                           &
    FracBiomHarvsted                  => plt_distb%FracBiomHarvsted,                  &
    iYearPlantHarvest_pft             => plt_distb%iYearPlantHarvest_pft,             &
    THIN_pft                          => plt_distb%THIN_pft,                          &
    FracCanopyHeightCut_pft           => plt_distb%FracCanopyHeightCut_pft,           &
    iHarvstType_pft                   => plt_distb%iHarvstType_pft,                   &
    jHarvst_pft                       => plt_distb%jHarvst_pft,                       &
    iDayPlanting_pft                  => plt_distb%iDayPlanting_pft,                  &
    iDayPlantHarvest_pft              => plt_distb%iDayPlantHarvest_pft,              &
    iYearPlanting_pft                 => plt_distb%iYearPlanting_pft,                 &
    GrainStrutElms_brch               => plt_biom%GrainStrutElms_brch,                &
    PopuRootMycoC_pvr                 => plt_biom% PopuRootMycoC_pvr,                 &
    SenecStalkStrutElms_brch          => plt_biom%SenecStalkStrutElms_brch,           &
    StalkStrutElms_brch               => plt_biom%StalkStrutElms_brch,                &
    RootMycoNonstElms_rpvr            => plt_biom%RootMycoNonstElms_rpvr,             &
    EarStrutElms_brch                 => plt_biom%EarStrutElms_brch,                  &
    HuskStrutElms_brch                => plt_biom%HuskStrutElms_brch,                 &
    SeasonalNonstElms_pft             => plt_biom%SeasonalNonstElms_pft,              &
    InternodeStrutElms_brch           => plt_biom%InternodeStrutElms_brch,            &
    GrainSeedBiomCMean_brch           => plt_allom%GrainSeedBiomCMean_brch,           &
    HourReq4LeafOff_brch              => plt_pheno%HourReq4LeafOff_brch,              &
    iPlantPhenolType_pft              => plt_pheno%iPlantPhenolType_pft,              &
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft,             &
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch,                &
    HourReq4LeafOut_brch              => plt_pheno%HourReq4LeafOut_brch,              &
    iPlantPhenolPattern_pft           => plt_pheno%iPlantPhenolPattern_pft,           &
    doPlantLeaveOff_brch              => plt_pheno%doPlantLeaveOff_brch,              &
    doInitPlant_pft                   => plt_pheno%doInitPlant_pft,                   &
    doPlantLeafOut_brch               => plt_pheno%doPlantLeafOut_brch,               &
    iPlantTurnoverPattern_pft         => plt_pheno%iPlantTurnoverPattern_pft,         &
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch,                &
    Prep4Literfall_brch               => plt_pheno%Prep4Literfall_brch,               &
    Hours4LiterfalAftMature_brch      => plt_pheno%Hours4LiterfalAftMature_brch,      &
    ElmAllocmat4Litr                  => plt_soilchem%ElmAllocmat4Litr,               &
    iYearCurrent                      => plt_site%iYearCurrent,                       &
    LitrfalStrutElms_pvr              => plt_bgcr%LitrfalStrutElms_pvr,               &
    MY                                => plt_morph%MY,                                &
    InternodeHeightDead_brch         => plt_morph%InternodeHeightDead_brch,         &
    LiveInterNodeHight_brch           => plt_morph%LiveInterNodeHight_brch,           &
    PotentialSeedSites_brch           => plt_morph%PotentialSeedSites_brch,           &
    MainBranchNum_pft                 => plt_morph%MainBranchNum_pft,                 &
    BranchNumber_brch                 => plt_morph%BranchNumber_brch,                 &
    SeedNumSet_brch                   => plt_morph%SeedNumSet_brch                    &
  )
!   RESET PHENOLOGY AT EMERGENCE ('Hours4Leafout_brch' > 'VRNL')
!   AND END OF SEASON ('Hours4LeafOff_brch' > 'VRNX')
!
!   iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!   iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!   doPlantLeafOut_brch=flag for enabling leafout:0=enable,1=disable
!   Hours4Leafout_brch,VRNL=leafout hours,hours required for leafout
!   doPlantLeaveOff_brch=flag for enabling leafoff:0=enable,1=disable
!   Hours4LeafOff_brch,VRNX=leafoff hours,hours required for leafoff
!
  IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual .OR. &
    (iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).GT.1))THEN

    LeafOutTest=(doPlantLeafOut_brch(NB,NZ).EQ.iEnable .AND. Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))
    LeafOffTest=(doPlantLeaveOff_brch(NB,NZ).EQ.iEnable .AND. Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ))
    IF(LeafOutTest .OR. LeafOffTest)THEN     
      !may include annual crop, which is set to evergreen
      CALL ResetNonAnnualBranch(I,J,NB,NZ)      
    ENDIF
  ENDIF
!
!   REPRODUCTIVE MATERIAL BECOMES LitrFall AT END OF SEASON
!

  IF(Prep4Literfall_brch(NB,NZ).EQ.itrue)THEN
    Hours4LiterfalAftMature_brch(NB,NZ)=Hours4LiterfalAftMature_brch(NB,NZ)+1
    IF(Hours4LiterfalAftMature_brch(NB,NZ).EQ.HoursReq4LiterfalAftMature)THEN
      Prep4Literfall_brch(NB,NZ)          = ifalse
      Hours4LiterfalAftMature_brch(NB,NZ) = 0
    ENDIF
    FSNR1=1.0_r8-FSNR

    D6330: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+ &
          FSNR*ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ))
        IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+ &
            FSNR*ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(GrainStrutElms_brch(NE,NB,NZ))
        ELSE
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+ &
            FSNR*ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*AZMAX1(GrainStrutElms_brch(NE,NB,NZ))
        ENDIF
      ENDDO
    ENDDO D6330

    DO NE=1,NumPlantChemElms  
      HuskStrutElms_brch(NE,NB,NZ)  = FSNR1*AZMAX1(HuskStrutElms_brch(NE,NB,NZ))
      EarStrutElms_brch(NE,NB,NZ)   = FSNR1*AZMAX1(EarStrutElms_brch(NE,NB,NZ))
      GrainStrutElms_brch(NE,NB,NZ) = FSNR1*AZMAX1(GrainStrutElms_brch(NE,NB,NZ))
    ENDDO
    PotentialSeedSites_brch(NB,NZ) = FSNR1*PotentialSeedSites_brch(NB,NZ)
    SeedNumSet_brch(NB,NZ)         = FSNR1*SeedNumSet_brch(NB,NZ)
    GrainSeedBiomCMean_brch(NB,NZ) = FSNR1*GrainSeedBiomCMean_brch(NB,NZ)
!
!     STALKS BECOME LitrFall IN GRASSES AT END OF SEASON
!
    IF((iPlantTurnoverPattern_pft(NZ).EQ.0 .OR. (.not.is_plant_treelike(iPlantRootProfile_pft(NZ)))) &
      .AND.iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN

      D6335: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)&
            +FSNR*ElmAllocmat4Litr(NE,istalk,M,NZ)*AZMAX1(StalkStrutElms_brch(NE,NB,NZ))
        ENDDO
      ENDDO D6335
      DO NE=1,NumPlantChemElms  
        StalkStrutElms_brch(NE,NB,NZ)=FSNR1*AZMAX1(StalkStrutElms_brch(NE,NB,NZ))
        SenecStalkStrutElms_brch(NE,NB,NZ)=FSNR1*SenecStalkStrutElms_brch(NE,NB,NZ)
      ENDDO
      DO K=0,MaxNodesPerBranch1
        DO NE=1,NumPlantChemElms
          InternodeStrutElms_brch(NE,K,NB,NZ)=FSNR1*InternodeStrutElms_brch(NE,K,NB,NZ)
        ENDDO
      ENDDO
      D2010: DO K=0,MaxNodesPerBranch1
    !     LiveInterNodeHight_brch(K,NB,NZ)=FSNR1*LiveInterNodeHight_brch(K,NB,NZ)
        InternodeHeightDead_brch(K,NB,NZ)=FSNR1*InternodeHeightDead_brch(K,NB,NZ)
      ENDDO D2010
    ENDIF

!
!   SELF-SEEDING ANNUALS IF COLD OR DROUGHT DECIDUOUS
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial
!     iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     iDayPlantHarvest_pft,iYearPlantHarvest_pft=day,year of harvesting
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     jHarvst_pft=terminate PFT:0=no,1=yes,2=yes,but reseed
!     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
!          iHarvstType_pft=3:reduction of clumping factor
!          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),iHarvstType_pft=5:fire
!     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
!          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
!     FracBiomHarvsted(1,1,FracBiomHarvsted(1,2,FracBiomHarvsted(1,3,FracBiomHarvsted(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     FracBiomHarvsted(2,1,FracBiomHarvsted(2,2,FracBiomHarvsted(2,3,FracBiomHarvsted(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosystem
!     iDayPlanting_pft,iYearPlanting_pft=day,year of planting
!     doInitPlant_pft=PFT initialization flag:0=no,1=yes
!
!     IF(J.EQ.INT(SolarNoonHour_col))THEN

    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      !deciduous annual plant      
      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).NE.0)THEN
        iDayPlantHarvest_pft(NZ)                           = I
        iYearPlantHarvest_pft(NZ)                          = iYearCurrent
        iHarvstType_pft(NZ)                                = iharvtyp_grain
        jHarvst_pft(NZ)                                    = jharvtyp_tmareseed
        FracCanopyHeightCut_pft(NZ)                        = 0._r8
        THIN_pft(NZ)                                       = 0._r8
        FracBiomHarvsted(ihav_pft,iplthvst_leaf,NZ)        = 1.0_r8
        FracBiomHarvsted(ihav_pft,iplthvst_finenonleaf,NZ) = 1.0_r8
        FracBiomHarvsted(ihav_pft,iplthvst_woody,NZ)       = 1.0_r8
        FracBiomHarvsted(ihav_pft,iplthvst_stdead,NZ)      = 1.0_r8
        FracBiomHarvsted(ihav_eco,iplthvst_leaf,NZ)        = 0._r8
        FracBiomHarvsted(ihav_eco,iplthvst_finenonleaf,NZ) = 1.0_r8
        FracBiomHarvsted(ihav_eco,iplthvst_woody,NZ)       = 0._r8
        FracBiomHarvsted(ihav_eco,iplthvst_stdead,NZ)      = 0._r8
        iDayPlanting_pft(NZ)                               = -1E+06
        iYearPlanting_pft(NZ)                              = -1E+06
        doInitPlant_pft(NZ)                                = itrue
!        write(134,*)I+J/24.,'self-seeding',NZ,doInitPlant_pft(NZ),iDayPlantHarvest_pft(NZ),iYearPlantHarvest_pft(NZ)
!        write(154,*)I+J/24.,'self-seeding',NZ,doInitPlant_pft(NZ),iDayPlantHarvest_pft(NZ),iYearPlantHarvest_pft(NZ)
      ENDIF
    ENDIF
!     ENDIF
  ENDIF
  end associate
  end subroutine ResetBranchPhenology
!------------------------------------------------------------------------------------------
  subroutine BranchElmntTransfer(I,J,NB,NZ,BegRemoblize,WaterStress4Groth,CanTurgPSIFun4Expans)
  use PlantNonstElmDynMod
  implicit none
  integer, intent(in) :: I,J,NB,NZ,BegRemoblize
  real(r8), intent(in) :: WaterStress4Groth
  real(r8), intent(in) :: CanTurgPSIFun4Expans
  integer :: NE
  real(r8) :: CPOOLD
  real(r8) :: FracCanopyCinStalk
  real(r8) :: ShootBiomC_brch
  real(r8) :: WVSTBX
  real(r8) :: WTRTTX,WTRSBX
  real(r8) :: WTRVCX
  real(r8) :: XFRE(1:NumPlantChemElms)
  logical :: PlantingChk,RemobChk,LeafOutChk,AnnualPlantChk
  ! begin_execution
  associate(                                                      &
    iDayPlanting_pft        => plt_distb%iDayPlanting_pft,        &
    iYearPlanting_pft       => plt_distb%iYearPlanting_pft,       &
    iYearCurrent            => plt_site%iYearCurrent,             &
    FracHour4LeafoffRemob   => plt_allom%FracHour4LeafoffRemob,   &
    RootElms_pft            => plt_biom%RootElms_pft,             &
    CanopyStalkC_pft        => plt_biom%CanopyStalkC_pft,         &
    CanopyNonstElms_brch    => plt_biom%CanopyNonstElms_brch,     &
    SeasonalNonstElms_pft   => plt_biom%SeasonalNonstElms_pft,    &
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft,           &
    StalkRsrvElms_brch      => plt_biom%StalkRsrvElms_brch,       &
    StalkLiveBiomassC_brch  => plt_biom%StalkLiveBiomassC_brch,   &
    iPlantCalendar_brch     => plt_pheno%iPlantCalendar_brch,     &
    fTCanopyGroth_pft       => plt_pheno%fTCanopyGroth_pft,       &
    HourReq4LeafOff_brch    => plt_pheno%HourReq4LeafOff_brch,    &
    Hours4Leafout_brch      => plt_pheno%Hours4Leafout_brch,      &
    HourReq4LeafOut_brch    => plt_pheno%HourReq4LeafOut_brch,    &
    doInitPlant_pft         => plt_pheno%doInitPlant_pft,         &
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft, &
    Hours4LeafOff_brch      => plt_pheno%Hours4LeafOff_brch,      &
    iPlantPhenolType_pft    => plt_pheno%iPlantPhenolType_pft,    &
    doInitLeafOut_brch      => plt_pheno%doInitLeafOut_brch,      &
    Hours2LeafOut_brch      => plt_pheno%Hours2LeafOut_brch,      &
    MainBranchNum_pft       => plt_morph%MainBranchNum_pft,       &
    MaxSoiL4Root_pft        => plt_morph%MaxSoiL4Root_pft         &
  )
!   TRANSFER C,N,P FROM SEASONAL STORAGE TO SHOOT AND ROOT
!   NON-STRUCTURAL C DURING SEED GERMINATION OR LEAFOUT
!
  PlantingChk=I.GE.iDayPlanting_pft(NZ) .AND. iYearCurrent.EQ.iYearPlanting_pft(NZ)
  RemobChk=Hours4LeafOff_brch(NB,NZ).LT.FracHour4LeafoffRemob(iPlantPhenolType_pft(NZ))*HourReq4LeafOff_brch(NB,NZ)  
  LeafOutChk=Hours4Leafout_brch(MainBranchNum_pft(NZ),NZ).GE.HourReq4LeafOut_brch(NB,NZ)
  AnnualPlantChk=iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. doInitPlant_pft(NZ).EQ.ifalse
  IF(AnnualPlantChk .OR.(PlantingChk.AND.RemobChk) .OR. (LeafOutChk.AND.RemobChk))THEN
!
  ! RESET TIME COUNTER
  !
  ! Hours2LeafOut_brch=hourly leafout counter
  ! doInitLeafOut_brch=flag for initializing leafout
  !
    IF(doInitLeafOut_brch(NB,NZ).EQ.iEnable)THEN
      Hours2LeafOut_brch(NB,NZ)=0._r8
      doInitLeafOut_brch(NB,NZ)=iDisable
    ENDIF
  !
  ! INCREMENT TIME COUNTER
  !
  ! iPlantPhotoperiodType_pft=photoperiod type:0=day neutral,1=short day,2=long day
  ! iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
  ! CriticPhotoPeriod_pft=critical photoperiod (h):<0=maximum daylength from site file
  ! PhotoPeriodSens_pft=photoperiod sensitivity (node h-1)
  ! DayLenthCurrent=daylength
  ! CanTurgPSIFun4Expans=expansion,extension function of canopy water potential
  ! fTCanopyGroth_pft=temperature function for canopy growth
  ! HourReq2InitSStor4LeafOut=number of hours required to initiate remobilization of storage C for leafout
  ! main branch leaf out
    
    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      call DevelopMainBranch(I,J,NB,NZ,CanTurgPSIFun4Expans)
    ENDIF
  !
  ! REDISTRIBUTE TRANFERRED C FROM MAIN STEM TO OTHER BRANCHES
  !
  ! ATRP=hourly leafout counter
  ! fTCanopyGroth_pft=temperature function for canopy growth
  ! HourReq2InitSStor4LeafOut=number of hours required for remobilization of storage C during leafout
  ! WaterStress4Groth=growth function of canopy water potential
  ! CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
  ! XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
  ! 
    IF(NB.NE.MainBranchNum_pft(NZ) &
      .AND. Hours2LeafOut_brch(NB,NZ).LE.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ)))THEN
      Hours2LeafOut_brch(NB,NZ)=Hours2LeafOut_brch(NB,NZ)+fTCanopyGroth_pft(NZ)*WaterStress4Groth
      DO NE=1,NumPlantChemElms
        XFRE(NE)=AZMAX1(0.05_r8*fTCanopyGroth_pft(NZ) &
          *(0.5_r8*(CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)+CanopyNonstElms_brch(NE,NB,NZ)) &
          -CanopyNonstElms_brch(NE,NB,NZ)))
        CanopyNonstElms_brch(NE,NB,NZ)                    = CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE)
        CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ) = CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)-XFRE(NE)
      ENDDO
    ENDIF
  ENDIF

!
! TRANSFER LEAF AND STALK NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
! IN PERENNIALS AFTER GRAIN FILL IN DETERMINATES, AFTER AUTUMNIZ'N
! IN INDETERMINATES, OR AFTER SUSTAINED WATER STRESS
!
!
  IF(BegRemoblize.EQ.itrue .AND. iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
    call SeasonStoreShootTransfer(I,J,NB,NZ)
  ENDIF
!
!   TRANSFER NON-STRUCTURAL C,N,P FROM LEAVES AND ROOTS TO RESERVES
!   IN STALKS DURING GRAIN FILL IN ANNUALS OR BETWEEN STALK RESERVES
!   AND LEAVES IN PERENNIALS ACCORDING TO CONCENTRATION DIFFERENCES
!
!   iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!   iPlantCalendar_brch(ipltcal_Jointing,=start of stem elongation and setting max seed number
!   iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date setting for final seed number
!   LeafPetolBiomassC_brch=leaf+petiole mass
!   StalkLiveBiomassC_brch=stalk sapwood mass

!
  IF((iPlantPhenolPattern_pft(NZ).EQ.iplt_annual  &
    .AND. iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).NE.0) &
    .OR. (iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial &
     .AND. iPlantCalendar_brch(ipltcal_Jointing,NB,NZ).NE.0))THEN

    call StalkRsrvShootNonstTransfer(I,J,NB,NZ)

    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.&
      iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).NE.0)THEN
      !stalk-root transfer
      call StalkRsrvRootNonstTransfer(I,J,NB,NZ)
    ENDIF

  ENDIF
!
!   REPLENISH BRANCH NON-STRUCTURAL POOL FROM
!   SEASONAL STORAGE POOL
!
!   StalkLiveBiomassC_brch,WVSTK=stalk,total stalk sapwood mass
!   WTRT=total root mass
!   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!   XFRX=maximum storage C content for remobiln from stalk,root reserves
!   XFRE(ielmc)=C transfer
!   Q: why are nitrogen and phosphorus not transferred?
  IF(StalkLiveBiomassC_brch(NB,NZ).GT.ZERO4Groth_pft(NZ)  &
    .AND. CanopyStalkC_pft(NZ).GT.ZERO4Groth_pft(NZ)  &
    .AND. RootElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ) &
    .AND. StalkRsrvElms_brch(ielmc,NB,NZ).LE.XFRX*StalkLiveBiomassC_brch(NB,NZ))THEN
    FracCanopyCinStalk              = StalkLiveBiomassC_brch(NB,NZ)/CanopyStalkC_pft(NZ)

    WVSTBX                          = StalkLiveBiomassC_brch(NB,NZ)
    WTRTTX                          = RootElms_pft(ielmc,NZ)*FracCanopyCinStalk
    ShootBiomC_brch                 = WVSTBX+WTRTTX
    WTRSBX                          = AZMAX1(StalkRsrvElms_brch(ielmc,NB,NZ))
    WTRVCX                          = AZMAX1(SeasonalNonstElms_pft(ielmc,NZ)*FracCanopyCinStalk)
    CPOOLD                          = (WTRVCX*WVSTBX-WTRSBX*WTRTTX)/ShootBiomC_brch

    XFRE(ielmc)                     = AZMAX1(XFRY*CPOOLD)
    StalkRsrvElms_brch(ielmc,NB,NZ) = StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
    SeasonalNonstElms_pft(ielmc,NZ) = SeasonalNonstElms_pft(ielmc,NZ)-XFRE(ielmc)
  ENDIF
  end associate
  end subroutine BranchElmntTransfer
!------------------------------------------------------------------------------------------

  subroutine DevelopMainBranch(I,J,NB,NZ,CanTurgPSIFun4Expans)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: CanTurgPSIFun4Expans
  real(r8)  :: TotPopuPlantRootC  
  real(r8) :: CH2OH
  real(r8) :: FXFC,FXFN  
  REAL(R8) :: cpoolt  
  real(r8) :: ATRPPD
  real(r8) :: NonstGradt  
  real(r8) :: PPDX,WFNSP
  real(r8) :: GFNX  
  real(r8) :: NonstElm2RootMyco(NumPlantChemElms)  
  real(r8) :: TotalRootNonstElms(1:NumPlantChemElms)  
  real(r8) :: ElmXferStore2Shoot(1:NumPlantChemElms)  
  real(r8) :: dTimeIncLeafout
  integer :: L,NE
  logical ::  PlantChk

  associate(                                                          &
    PhotoPeriodSens_pft       => plt_pheno%PhotoPeriodSens_pft,       &
    NU                        => plt_site%NU,                         &
    Hours2LeafOut_brch        => plt_pheno%Hours2LeafOut_brch,        &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,             &
    CriticPhotoPeriod_pft     => plt_pheno%CriticPhotoPeriod_pft,     &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,     &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,     &
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft,      &
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch,       &
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft,      &
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr,         &
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,        &
    fTCanopyGroth_pft         => plt_pheno%fTCanopyGroth_pft,         &
    iPlantPhenolPattern_pft   => plt_pheno%iPlantPhenolPattern_pft,   &
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft, &
    iPlantPhotoperiodType_pft => plt_pheno%iPlantPhotoperiodType_pft, &
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,          &
    DayLenthCurrent           => plt_site%DayLenthCurrent             &
  )  

  TotPopuPlantRootC         = 0._r8
  TotalRootNonstElms(ielmc) = 0._r8
  D4: DO L=NU,MaxSoiL4Root_pft(NZ)
    TotPopuPlantRootC         = TotPopuPlantRootC+AZMAX1(PopuRootMycoC_pvr(ipltroot,L,NZ))
    TotalRootNonstElms(ielmc) = TotalRootNonstElms(ielmc)+AZMAX1(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))
  ENDDO D4

  IF(iPlantPhotoperiodType_pft(NZ).EQ.iphotop_long &
    .AND.(iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid &
    .OR.iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid))THEN
    PPDX   = AZMAX1(CriticPhotoPeriod_pft(NZ)-PhotoPeriodSens_pft(NZ)-DayLenthCurrent)
    ATRPPD = EXP(-0.0_r8*PPDX)
  ELSE
    ATRPPD=1.0_r8
  ENDIF
      !
  IF(.not.is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
    WFNSP=CanTurgPSIFun4Expans
  ELSE
    WFNSP=1.0_r8
  ENDIF

  dTimeIncLeafout           = ATRPPD*fTCanopyGroth_pft(NZ)*WFNSP
  Hours2LeafOut_brch(NB,NZ) = Hours2LeafOut_brch(NB,NZ)+dTimeIncLeafout
  PlantChk                  = iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND. iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen
  IF(Hours2LeafOut_brch(NB,NZ).LE.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ)) .OR. plantChk)THEN

    IF(SeasonalNonstElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
      !CPOOLT:=total root nonst + branch nonst
      CPOOLT=TotalRootNonstElms(ielmc)+CanopyNonstElms_brch(ielmc,NB,NZ)
!
!       REMOBILIZE C FROM SEASONAL STORAGE AT FIRST-ORDER RATE
!       MODIFIED BY SOIL TEMPERATURE AT SEED DEPTH
!
!     GVMX=specific oxidation rate of storage C during leafout at 25 C
!     WTRVC=storage C
!     CH2OH=storage C oxidation rate during leafout
!     CPOOL,CPOOLR=non-structural C mass in branch,root
!     FXSH,FXRT=shoot-root partitioning of storage C during leafout
!     WTRTD=root C mass
!
      GFNX                              = GVMX(iPlantPhenolPattern_pft(NZ))*dTimeIncLeafout
      CH2OH                             = AZMAX1(GFNX*SeasonalNonstElms_pft(ielmc,NZ))
      SeasonalNonstElms_pft(ielmc,NZ)   = SeasonalNonstElms_pft(ielmc,NZ)-CH2OH
      CanopyNonstElms_brch(ielmc,NB,NZ) = CanopyNonstElms_brch(ielmc,NB,NZ)+CH2OH*FXSH(iPlantPhenolPattern_pft(NZ))
!         if root condition met
      IF(TotPopuPlantRootC.GT.ZERO4Groth_pft(NZ).AND.TotalRootNonstElms(ielmc).GT.ZERO4Groth_pft(NZ))THEN
        D50: DO L=NU,MaxSoiL4Root_pft(NZ)
          FXFC=AZMAX1(PopuRootMycoC_pvr(ipltroot,L,NZ))/TotPopuPlantRootC
          RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ) &
              +FXFC*CH2OH*FXRT(iPlantPhenolPattern_pft(NZ))
        ENDDO D50
      ELSE
          RootMycoNonstElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),NZ)=&
            RootMycoNonstElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),NZ)+CH2OH &
            *FXRT(iPlantPhenolPattern_pft(NZ))
      ENDIF

    ELSE
      CH2OH=0._r8
    ENDIF
  ELSE
    CH2OH=0._r8
  ENDIF
      !
      !     REMOBILIZE N,P FROM SEASONAL STORAGE AT FIRST-ORDER RATE
      !     MODIFIED BY SOIL TEMPERATURE AT SEED DEPTH
      !
      !     WTRVC,WTRVN,WTRVP=storage C,N,P
      !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
      !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
      !     NXferStore2Shoot,PXferStore2Shoot=N,P transfer from storage to shoot
      !     CH2OH=storage C oxidation rate during leafout
      !     FRSV=rate constant for remobiln of storage C,N,P during leafout C
      !     FXSH=shoot partitioning of storage C during leafout
      !
  IF(SeasonalNonstElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
      CPOOLT=AZMAX1(SeasonalNonstElms_pft(ielmc,NZ)+CanopyNonstElms_brch(ielmc,NB,NZ))          
      DO NE=2,NumPlantChemElms
        NonstGradt=(SeasonalNonstElms_pft(NE,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ)- &
          CanopyNonstElms_brch(NE,NB,NZ)*SeasonalNonstElms_pft(ielmc,NZ))/CPOOLT            
        ElmXferStore2Shoot(NE)=AZMAX1(FRSV(iPlantTurnoverPattern_pft(NZ))*NonstGradt)
      ENDDO

    ELSE
      DO NE=2,NumPlantChemElms
        ElmXferStore2Shoot(NE)=AZMAX1(FXSH(iPlantPhenolPattern_pft(NZ))*CH2OH*SeasonalNonstElms_pft(NE,NZ) &
          /SeasonalNonstElms_pft(ielmc,NZ))
      ENDDO  
    ENDIF
  ELSE
    !when there is no seasonal storage C
    DO NE=2,NumPlantChemElms
      ElmXferStore2Shoot(NE)=AZMAX1(FXSH(iPlantPhenolPattern_pft(NZ))*SeasonalNonstElms_pft(NE,NZ))
    ENDDO
  ENDIF
!
! ADD TO NON-STRUCTURAL POOLS IN ROOT
!
! CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
! WTRVC,WTRVN,WTRVP=storage C,N,P
! iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
! UPNH4R,UPPO4R=N,P transfer from storage to root
! FRSV=rate constant for remobiln of storage C,N,P during leafout
! FXRT=root partitioning of storage C during leafout
!

  DO NE=1,NumPlantChemElms
    TotalRootNonstElms(NE)=0._r8
    D3: DO L=NU,MaxSoiL4Root_pft(NZ)
      TotalRootNonstElms(NE)=TotalRootNonstElms(NE)+ &
        AZMAX1(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ))
    ENDDO D3
  ENDDO

  IF(SeasonalNonstElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
      CPOOLT=AMAX1(ZERO4Groth_pft(NZ),SeasonalNonstElms_pft(ielmc,NZ)+TotalRootNonstElms(ielmc))
      DO NE=2,NumPlantChemElms
        NonstGradt=(SeasonalNonstElms_pft(NE,NZ)*TotalRootNonstElms(ielmc)-&
          TotalRootNonstElms(NE)*SeasonalNonstElms_pft(ielmc,NZ))/CPOOLT            
        NonstElm2RootMyco(NE)=AZMAX1(FRSV(iPlantTurnoverPattern_pft(NZ))*NonstGradt)
      ENDDO
    ELSE
      DO NE=2,NumPlantChemElms
        NonstElm2RootMyco(NE)=AZMAX1(FXRT(iPlantPhenolPattern_pft(NZ))*CH2OH*&
          SeasonalNonstElms_pft(NE,NZ)/SeasonalNonstElms_pft(ielmc,NZ))
      ENDDO    
    ENDIF
  ELSE
    DO NE=2,NumPlantChemElms
      NonstElm2RootMyco(NE)=AZMAX1(FXRT(iPlantPhenolPattern_pft(NZ))*SeasonalNonstElms_pft(NE,NZ))
    ENDDO  
  ENDIF
  !
!     TRANSFER STORAGE FLUXES
!
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     NXferStore2Shoot,PXferStore2Shoot=N,P transfer from storage to shoot
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     UPNH4R,UPPO4R=N,P transfer from storage to root
!     FXFN=root layer allocation
!
  DO NE=2,NumPlantChemElms
    SeasonalNonstElms_pft(NE,NZ)   = SeasonalNonstElms_pft(NE,NZ)-ElmXferStore2Shoot(NE)-NonstElm2RootMyco(NE)
    CanopyNonstElms_brch(NE,NB,NZ) = CanopyNonstElms_brch(NE,NB,NZ)+ElmXferStore2Shoot(NE)
  ENDDO
  
  IF(TotPopuPlantRootC.GT.ZERO4Groth_pft(NZ).AND.TotalRootNonstElms(ielmc).GT.ZERO4Groth_pft(NZ))THEN
    D51: DO L=NU,MaxSoiL4Root_pft(NZ)
      FXFN=AZMAX1(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))/TotalRootNonstElms(ielmc)
      DO NE=2,NumPlantChemElms
        RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ) &
          +FXFN*NonstElm2RootMyco(NE)
      ENDDO  
    ENDDO D51
  ELSE
    DO NE=2,NumPlantChemElms
      RootMycoNonstElms_rpvr(NE,ipltroot,NGTopRootLayer_pft(NZ),NZ)=&
        RootMycoNonstElms_rpvr(NE,ipltroot,NGTopRootLayer_pft(NZ),NZ)+NonstElm2RootMyco(NE)   
    ENDDO            
  ENDIF
  end associate
  end subroutine DevelopMainBranch  
!------------------------------------------------------------------------------------------

  subroutine ComputRAutoAfEmergence(I,J,NB,NZ,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,CO2F,&
    CH2O,TFN5,WaterStress4Groth,CanTurgPSIFun4Expans,ShootStructN,CanopyNonstElm4Gros,CNPG,&
    RCO2NonstC_brch,RCO2Maint_brch,RMxess_brch,NonstC4Groth_brch,RCO2NonstC4Nassim_brch)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: DMSHD
  real(r8), intent(in) :: CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX
  real(r8), intent(in) :: ShootStructN,CO2F
  real(r8), intent(in) :: CH2O  !total CH2O production
  real(r8), intent(in) :: TFN5  !temperature function for canopy maintenance respiration
  real(r8), intent(in) :: WaterStress4Groth
  real(r8), intent(in) :: CanTurgPSIFun4Expans
  real(r8), intent(out) :: CanopyNonstElm4Gros(NumPlantChemElms)
  real(r8), INTENT(OUT) :: CNPG
  real(r8), intent(out) :: RCO2NonstC_brch
  real(r8), intent(out) :: RCO2Maint_brch
  real(r8), intent(out) :: RMxess_brch
  real(r8), intent(out) :: NonstC4Groth_brch
  real(r8), intent(out) :: RCO2NonstC4Nassim_brch
  real(r8) :: ZPOOLB
  real(r8) :: PPOOLB
  real(r8) :: RCO2X,RCO2Y,RGFNP
  real(r8) :: RgroCO2_ltd  !Nutient limited growth respiration
  real(r8) :: Rauto_brch,RCO2CM
  real(r8) :: CNG,CPG,RFN1,RFP1
! begin_execution
  associate(                                                         &
    CO2NetFix_pft             => plt_bgcr%CO2NetFix_pft,             &
    CanopyGrosRCO2_pft        => plt_bgcr%CanopyGrosRCO2_pft,        &
    Eco_AutoR_CumYr_col       => plt_bgcr%Eco_AutoR_CumYr_col,       &
    ECO_ER_col                => plt_bgcr%ECO_ER_col,                &
    GrossCO2Fix_pft           => plt_bgcr%GrossCO2Fix_pft,           &
    CanopyRespC_CumYr_pft     => plt_bgcr%CanopyRespC_CumYr_pft,     &
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch,      &
    LeafPetoNonstElmConc_brch => plt_biom%LeafPetoNonstElmConc_brch, &
    ZERO                      => plt_site%ZERO,                      &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,    &
    fTCanopyGroth_pft         => plt_pheno%fTCanopyGroth_pft,        &
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft,     &
    C4PhotosynDowreg_brch     => plt_photo%C4PhotosynDowreg_brch      &
  )
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
! CNKI,CPKI=nonstruct N,P inhibn constant on growth(g N,P g-1 C)
!
  IF(LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZERO)THEN
    CNG=LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI)
    CPG=LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI)
    CNPG=AMIN1(CNG,CPG)
  ELSE
    CNPG=1.0_r8
  ENDIF
!
! RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
! NON-STRUCTURAL C:N:P
!
! RCO2NonstC_brch=respiration from non-structural C
! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
! CPOOL=non-structural C mass
! fTCanopyGroth_pft=temperature function for canopy growth
! WaterStress4Groth=growth function of canopy water potential
! CNPG=N,P constraint on respiration
! C4PhotosynDowreg_brch=termination feedback inhibition on C3 CO2
!
  RCO2NonstC_brch=AZMAX1(VMXC*CanopyNonstElms_brch(ielmc,NB,NZ) &
    *fTCanopyGroth_pft(NZ))*CNPG*C4PhotosynDowreg_brch(NB,NZ)*WaterStress4Groth
!
! MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
!
! RCO2Maint_brch=maintenance respiration
! TFN5=
! ShootStructE_brch(ielmn)=shoot structural N mass
! iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
! iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
! WaterStress4Groth=growth function of canopy water potential
!
  RCO2Maint_brch=AZMAX1(RmSpecPlant*TFN5*ShootStructN)
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)) .OR. &
    iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    RCO2Maint_brch=RCO2Maint_brch*WaterStress4Groth
  ENDIF
  
!
! GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
! IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
! RCO2X=difference between non-structural C respn and mntc respn
! RCO2Y=growth respiration unlimited by N,P
! CanTurgPSIFun4Expans=expansion,extension function of canopy water potential
! RMxess_brch=excess maintenance respiration, drives remobilization & senescence
!
  RCO2X       = RCO2NonstC_brch-RCO2Maint_brch
  RCO2Y       = AZMAX1(RCO2X)*CanTurgPSIFun4Expans  
  RMxess_brch = AZMAX1(-RCO2X)
!  write(111,*)I+J/24.,RCO2NonstC_brch,RCO2Maint_brch,TFN5,ShootStructN,fTCanopyGroth_pft(NZ),WaterStress4Groth
!
! GROWTH RESPIRATION MAY BE LIMITED BY NON-STRUCTURAL N,P
! AVAILABLE FOR GROWTH
!
! RCO2Y,RgroCO2_ltd=growth respiration unlimited,limited by N,P
! CNLFX=diff between min and max leaf N prodn vs nonstruct C consumption
! CNSHX=N production vs nonstructural C consumption in rest of shoot
! ZPOOL,PPOOL=nonstructural N,P mass
! DMSHD=branch C respiration vs nonstructural C consumption
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
!
  IF(RCO2Y.GT.0.0_r8 .AND. (CNSHX.GT.0.0_r8 .OR. CNLFX.GT.0.0_r8))THEN
    ZPOOLB      = AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))
    PPOOLB      = AZMAX1(CanopyNonstElms_brch(ielmp,NB,NZ))
    RFN1=ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG)
    RFP1=PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG)
    RGFNP       = AMIN1(RFN1,RFP1)
    RgroCO2_ltd = AMIN1(RCO2Y,RGFNP)
  ELSE
    RgroCO2_ltd=0._r8
  ENDIF

!
! TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
! AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
! ENTERED IN 'READQ'
!
! NonstC4Groth_brch=total non-structural C used in growth and growth respiration
! RgroCO2_ltd=growth respiration limited by N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! CanopyNonstElm4Gros(ielmn),CanopyNonstElm4Gros(ielmp)=nonstructural N,P used in growth
! ZPOOL,PPOOL=nonstructural N,P mass
! CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! RCO2NonstC4Nassim_brch=respiration for N assimilation
! CH2O= 1/4 newly fixed C is used for N-assimilation (where is this number from?)
!
  
  NonstC4Groth_brch          = RgroCO2_ltd/DMSHD
!  write(116,*)I*1000+J,NonstC4Groth_brch,RCO2Y,CanTurgPSIFun4Expans,WaterStress4Groth,CanopyNonstElms_brch(ielmc,NB,NZ)  
  CanopyNonstElm4Gros(ielmn) = AZMAX1(AMIN1(CanopyNonstElms_brch(ielmn,NB,NZ),NonstC4Groth_brch*(CNSHX+CNLFM+CNLFX*CNPG)))
  CanopyNonstElm4Gros(ielmp) = AZMAX1(AMIN1(CanopyNonstElms_brch(ielmp,NB,NZ),NonstC4Groth_brch*(CPSHX+CPLFM+CPLFX*CNPG)))
  RCO2NonstC4Nassim_brch     = AZMAX1(1.70_r8*CanopyNonstElm4Gros(ielmn)-0.025_r8*CH2O)
!
! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
!
! Rauto_brch=total C respiration
! RCO2Maint_brch=maintenance respiration
! RCO2NonstC_brch=respiration from non-structural C
! RgroCO2_ltd=growth respiration limited by N,P
! RMxess_brch=excess maintenance respiration
! RCO2NonstC4Nassim_brch=respiration for N assimilation
! GrossCO2Fix_pft=total PFT CO2 fixation
! CO2F=total CO2 fixation
! CanopyGrosRCO2_pft,CanopyRespC_CumYr_pft=total,above-ground PFT respiration
! CO2NetFix_pft=PFT net CO2 fixation
! ECO_ER_col=ecosystem respiration
! Eco_AutoR_CumYr_col=total autotrophic respiration
!
  Rauto_brch                = AMIN1(RCO2Maint_brch,RCO2NonstC_brch)+RgroCO2_ltd+RMxess_brch+RCO2NonstC4Nassim_brch
  GrossCO2Fix_pft(NZ)       = GrossCO2Fix_pft(NZ)+CO2F
  CanopyGrosRCO2_pft(NZ)    = CanopyGrosRCO2_pft(NZ)-Rauto_brch
  CanopyRespC_CumYr_pft(NZ) = CanopyRespC_CumYr_pft(NZ)-Rauto_brch
  CO2NetFix_pft(NZ)         = CO2NetFix_pft(NZ)+CO2F-Rauto_brch
  ECO_ER_col                = ECO_ER_col-Rauto_brch
  Eco_AutoR_CumYr_col       = Eco_AutoR_CumYr_col-Rauto_brch

  end associate
  end subroutine ComputRAutoAfEmergence

!------------------------------------------------------------------------------------------

  subroutine ComputRAutoB4Emergence(I,J,NB,NZ,TFN6_vr,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,&
    CNLFX,CPLFX,ShootStructN,WaterStress4Groth,CanTurgPSIFun4Expans,CanopyNonstElm4Gros,CNPG,RCO2NonstC_brch,&
    RCO2Maint_brch,RMxess_brch,NonstC4Groth_brch,RCO2NonstC4Nassim_brch)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8),intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX
  real(r8), intent(in) :: ShootStructN
  real(r8), intent(in) :: WaterStress4Groth
  real(r8), intent(in) :: CanTurgPSIFun4Expans
  real(r8), intent(out) :: CanopyNonstElm4Gros(NumPlantChemElms)
  real(r8), intent(out) :: RCO2NonstC_brch
  real(r8), INTENT(OUT) :: CNPG,RCO2Maint_brch,RMxess_brch
  real(r8), intent(out) :: NonstC4Groth_brch    !nonstrucal C to drive biomass growth
  real(r8), intent(out) :: RCO2NonstC4Nassim_brch

  character(len=*), parameter :: subname='ComputRAutoB4Emergence'
  real(r8) :: ZPOOLB,ZADDBM,CGROSM
  real(r8) :: RCO2NonstC4NassimOUltd                !O2 unlimited respiration for N assimilation
  real(r8) :: FNP,RFN1,RFP1
  real(r8) :: PPOOLB
  real(r8) :: RCO2X,RCO2Y
  real(r8) :: RgroCO2_ltd   !oxygen & nutrient-limited growth respiration
  real(r8) :: Rauto_brch,RCO2CM
  real(r8) :: RCO2X_O2ulm,RCO2YM
  real(r8) :: RCO2GM
  real(r8) :: RCO2TM
  real(r8) :: SenesMaxByMaintDefcit,CNG,CPG
! begin_execution
  associate(                                                         &
    LeafPetoNonstElmConc_brch => plt_biom%LeafPetoNonstElmConc_brch, &
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch,      &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,    &
    fTgrowRootP_vr            => plt_pheno%fTgrowRootP_vr,           &
    RAutoRootO2Limter_rpvr    => plt_rbgc%RAutoRootO2Limter_rpvr,    &
    iPlantPhenolType_pft      => plt_pheno%iPlantPhenolType_pft,     &
    CO2NetFix_pft             => plt_bgcr%CO2NetFix_pft,             &
    RootRespPotent_pvr        => plt_rbgc%RootRespPotent_pvr,        &
    RootCO2EmisPot_pvr        => plt_rbgc%RootCO2EmisPot_pvr,        &
    RootCO2Autor_pvr          => plt_rbgc%RootCO2Autor_pvr,          &
    ZERO                      => plt_site%ZERO,                      &
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,       &
    C4PhotosynDowreg_brch     => plt_photo%C4PhotosynDowreg_brch     &
  )
!
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
  call PrintInfo('beg '//subname)
  IF(LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZERO)THEN
    CNG=LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI)
    CPG=LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI) 
    CNPG=AMIN1(CNG,CPG)
  ELSE
    CNPG=1.0_r8
  ENDIF
!
! RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
! NON-STRUCTURAL C:N:P, O2 UPTAKE
!
! RCO2CM,RCO2NonstC_brch=respiration from non-structural C unlimited,limited by O2
! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
! CPOOL=non-structural C mass
! fTgrowRootP_vr=temperature function for root growth
! WaterStress4Groth=growth function of canopy water potential
! CNPG=N,P constraint on respiration
! C4PhotosynDowreg_brch=termination feedback inhibition on C3 CO2
! RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
!
  RCO2CM=AZMAX1(VMXC*CanopyNonstElms_brch(ielmc,NB,NZ) &
    *fTgrowRootP_vr(NGTopRootLayer_pft(NZ),NZ))*WaterStress4Groth*CNPG*C4PhotosynDowreg_brch(NB,NZ)
  RCO2NonstC_brch=RCO2CM*RAutoRootO2Limter_rpvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)  
!
! MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
!
! RCO2Maint_brch=maintenance respiration
! TFN6_vr=temperature function for root maintenance respiration
! ShootStructE_brch(ielmn)=shoot structural N mass
! iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
! iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
! WaterStress4Groth=growth function of canopy water potential
!
  RCO2Maint_brch=AZMAX1(RmSpecPlant*TFN6_vr(NGTopRootLayer_pft(NZ))*ShootStructN)
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)) .OR. &
    iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    RCO2Maint_brch=RCO2Maint_brch*WaterStress4Groth
  ENDIF
!
! GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
! IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
! RCO2X_O2ulm,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
! RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
! CanTurgPSIFun4Expans=expansion,extension function of canopy water potential
! SenesMaxByMaintDefcit,RMxess_brch=excess maintenance respiration unltd,ltd by O2
!
  RCO2X_O2ulm           = RCO2CM-RCO2Maint_brch
  RCO2X                 = RCO2NonstC_brch-RCO2Maint_brch
  RCO2YM                = AZMAX1(RCO2X_O2ulm)*CanTurgPSIFun4Expans
  RCO2Y                 = AZMAX1(RCO2X)*CanTurgPSIFun4Expans
  SenesMaxByMaintDefcit = AZMAX1(-RCO2X_O2ulm)
  RMxess_brch           = AZMAX1(-RCO2X)
!  write(111,*)I+J/24.,RCO2NonstC_brch,RCO2Maint_brch,TFN6_vr(NGTopRootLayer_pft(NZ)),ShootStructN, &
!    fTgrowRootP_vr(NGTopRootLayer_pft(NZ),NZ),WaterStress4Groth,RAutoRootO2Limter_rpvr(ipltroot,NGTopRootLayer_pft(NZ),NZ),&
!    NGTopRootLayer_pft(NZ)  
!
! GROWTH RESPIRATION MAY BE LIMITED BY NON-STRUCTURAL N,P
! AVAILABLE FOR GROWTH
!
! RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
! CNLFX=diff between min and max leaf N prodn vs nonstruct C consumption
! CNSHX=N production vs nonstructural C consumption in rest of shoot
! ZPOOL,PPOOL=nonstructural N,P mass
! FNP=growth respiration limited by O2 and N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! RCO2GM,RgroCO2_ltd=growth respiration unltd,ltd by O2 and limited by N,P
! RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
!
  IF(CNSHX.GT.0.0_r8.OR.CNLFX.GT.0.0_r8)THEN
    ZPOOLB = AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))
    PPOOLB = AZMAX1(CanopyNonstElms_brch(ielmp,NB,NZ))
    RFN1=ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG)
    RFP1=PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG)
    FNP    = AMIN1(RFN1,RFP1)
    IF(RCO2YM.GT.0.0_r8)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
    ELSE
      RCO2GM=0._r8
    ENDIF
    IF(RCO2Y.GT.0.0_r8)THEN
      RgroCO2_ltd=AMIN1(RCO2Y,FNP*RAutoRootO2Limter_rpvr(ipltroot,NGTopRootLayer_pft(NZ),NZ))
!      write(123,*)I+J/24.,RCO2Y,RFN1*RAutoRootO2Limter_rpvr(ipltroot,NGTopRootLayer_pft(NZ),NZ),&
!        RFP1*RAutoRootO2Limter_rpvr(ipltroot,NGTopRootLayer_pft(NZ),NZ),DMSHD,&
!        RAutoRootO2Limter_rpvr(ipltroot,NGTopRootLayer_pft(NZ),NZ),RCO2CM
    ELSE
      RgroCO2_ltd=0._r8
    ENDIF

  ELSE
    RCO2GM=0._r8
    RgroCO2_ltd=0._r8
  ENDIF
!
! TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
! AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
! ENTERED IN 'READQ'
!
! CGROSM,NonstC4Groth_brch=total non-structural C used in growth and respn unltd,ltd by O2
! RCO2GM,RgroCO2_ltd=growth respiration unltd,ltd by O2 and limited by N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! ZADDBM,CanopyNonstElm4Gros(ielmn),CanopyNonstElm4Gros(ielmp)=nonstructural N,P unltd,ltd by O2 used in growth
! ZPOOL,PPOOL=nonstructural N,P mass
! CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! RCO2NonstC4NassimOUltd,RCO2NonstC4Nassim_brch=respiration for N assimilation unltd,ltd by O2
!
  CGROSM                     = RCO2GM/DMSHD
  NonstC4Groth_brch          = RgroCO2_ltd/DMSHD
  ZADDBM                     = AZMAX1(CGROSM*(CNSHX+CNLFM+CNLFX*CNPG))
  CanopyNonstElm4Gros(ielmn) = AZMAX1(NonstC4Groth_brch*(CNSHX+CNLFM+CNLFX*CNPG))
  CanopyNonstElm4Gros(ielmp) = AZMAX1(NonstC4Groth_brch*(CPSHX+CPLFM+CPLFX*CNPG))
  RCO2NonstC4NassimOUltd     = AZMAX1(1.70_r8*ZADDBM)
  RCO2NonstC4Nassim_brch     = AZMAX1(1.70_r8*CanopyNonstElm4Gros(ielmn))
!
! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
!
! RCO2TM,Rauto_brch=total C respiration unltd,ltd by O2
! RCO2Maint_brch=maintenance respiration
! RCO2GM,RgroCO2_ltd=growth respiration limited by N,P unltd,ltd by O2
! SenesMaxByMaintDefcit,RMxess_brch=excess maintenance respiration unltd,ltd by O2
! RCO2NonstC4NassimOUltd,RCO2NonstC4Nassim_brch=respiration for N assimilation unltd,ltd by O2
! RootCO2Autor_pvr=total root respiration
! RootRespPotent_pvr,RootCO2EmisPot_pvr=RootCO2Autor_pvr unltd by O2,nonstructural C
!
  RCO2TM                                                 = RCO2Maint_brch+RCO2GM+SenesMaxByMaintDefcit+RCO2NonstC4NassimOUltd
  Rauto_brch                                             = RCO2Maint_brch+RgroCO2_ltd+RMxess_brch+RCO2NonstC4Nassim_brch
  RootRespPotent_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ) = RootRespPotent_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+RCO2TM
  RootCO2EmisPot_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ) = RootCO2EmisPot_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+Rauto_brch
  RootCO2Autor_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)   = RootCO2Autor_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)-Rauto_brch

  call PrintInfo('end '//subname)
  end associate
  end subroutine ComputRAutoB4Emergence

!------------------------------------------------------------------------------------------

  subroutine GrowLeavesOnBranch(I,J,NZ,NB,MinNodeNum,GrowthLeaf,EtoliationCoeff,WFNS,ALLOCL)
  implicit none
  integer , intent(in) :: I,J  
  integer , intent(in) :: NB,NZ
  integer , intent(in) :: MinNodeNum
  real(r8), intent(in) :: GrowthLeaf(NumPlantChemElms)
  real(r8), intent(in) :: EtoliationCoeff
  real(r8), intent(in) :: WFNS
  REAL(R8), INTENT(OUT):: ALLOCL
  integer :: K,NE,KK
  integer :: MXNOD,MNNOD,KNOD
  real(r8) :: GNOD
  REAL(R8) :: GrowthSLA,LeafAreaGrowth,SpecAreaLeafGrowth
  REAL(R8) :: GrowthElms(NumPlantChemElms)

  associate(                                                    &
    NumCogrowthNode_pft     => plt_morph%NumCogrowthNode_pft,   &
    LeafNodeArea_brch      => plt_morph%LeafNodeArea_brch,      &
    LeafAreaLive_brch      => plt_morph%LeafAreaLive_brch,      &
    SLA1_pft               => plt_morph%SLA1_pft,               &
    LeafElmntNode_brch     => plt_biom%LeafElmntNode_brch,      &
    LeafProteinCNode_brch  => plt_biom%LeafProteinCNode_brch,   &
    KHiestGroLeafNode_brch => plt_pheno%KHiestGroLeafNode_brch, &
    ZERO4LeafVar_pft       => plt_biom%ZERO4LeafVar_pft,        &
    rCNNonstRemob_pft      => plt_allom%rCNNonstRemob_pft,      &
    rCPNonstRemob_pft      => plt_allom%rCPNonstRemob_pft,      &
    FracGroth2Node_pft     => plt_allom%FracGroth2Node_pft,     &
    PlantPopulation_pft    => plt_site%PlantPopulation_pft      &
  )
  
  IF(GrowthLeaf(ielmc).GT.0.0_r8)THEN
    MXNOD  = KHiestGroLeafNode_brch(NB,NZ)
    MNNOD  = MAX(MinNodeNum,MXNOD-NumCogrowthNode_pft(NZ)+1)
    MXNOD  = MAX(MXNOD,MNNOD)
    KNOD   = MXNOD-MNNOD+1
    GNOD   = KNOD           !number of growing nodes
    ALLOCL = 1.0_r8/GNOD
    DO NE     = 1, NumPlantChemElms
      GrowthElms(NE)=ALLOCL*GrowthLeaf(NE)
    ENDDO
    GrowthSLA=ALLOCL*FracGroth2Node_pft(NZ)*NumCogrowthNode_pft(NZ)
!
!     GROWTH AT EACH CURRENT NODE
!
!     WGLF,WGLFN,WGLFP,LeafProteinCNode_brch=node leaf C,N,P,protein mass
!     GRO,GrowthElms(ielmn),GrowthElms(ielmp)=leaf C,N,P growth at each node
!     CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios from startq.f
!
    D490: DO KK=MNNOD,MXNOD
      K=MOD(KK,MaxNodesPerBranch1)
      IF(K.EQ.0.AND.KK.NE.0)K=MaxNodesPerBranch1
        DO NE=1,NumPlantChemElms
          LeafElmntNode_brch(NE,K,NB,NZ)=LeafElmntNode_brch(NE,K,NB,NZ)+GrowthElms(NE)
        ENDDO
        LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)+ &
          AMIN1(GrowthElms(ielmn)*rCNNonstRemob_pft(NZ),GrowthElms(ielmp)*rCPNonstRemob_pft(NZ))
!
!         SPECIFIC LEAF AREA FUNCTION OF CURRENT LEAF MASS
!         AT EACH NODE
!
!         SpecAreaLeafGrowth=specific area of leaf growth
!         ETOL=coefficient for etoliation effects on expansion,extension
!         SLA1_pft=growth in leaf area vs mass from PFT file
!         SLA2=parameter for calculating leaf area expansion
!         WGLF=leaf C mass
!         PP=PFT population
!         GrowthSLA=allocation of leaf area growth to each node
!         WFNS=turgor expansion,extension function
!         LeafAreaGrowth,GRO=leaf area,mass growth
!         LeafAreaLive_brch,LeafNodeArea_brch=branch,node leaf area
!
      SpecAreaLeafGrowth=EtoliationCoeff*SLA1_pft(NZ)*(AMAX1(ZERO4LeafVar_pft(NZ) &
        ,LeafElmntNode_brch(ielmc,K,NB,NZ))/(PlantPopulation_pft(NZ)*GrowthSLA))**SLA2*WFNS
      LeafAreaGrowth             = GrowthElms(ielmc)*SpecAreaLeafGrowth
      LeafAreaLive_brch(NB,NZ)   = LeafAreaLive_brch(NB,NZ)+LeafAreaGrowth
      LeafNodeArea_brch(K,NB,NZ) = LeafNodeArea_brch(K,NB,NZ)+LeafAreaGrowth
    ENDDO D490
  ENDIF
  end associate  
  END subroutine GrowLeavesOnBranch
!------------------------------------------------------------------------------------------
  subroutine GrowPetioleOnBranch(NZ,NB,MinNodeNum,GrowthPetiole,EtoliationCoeff,WFNS,ALLOCL)
  implicit none
  integer, intent(in) :: NZ,MinNodeNum,NB
  real(r8), intent(in) :: GrowthPetiole(NumPlantChemElms)
  REAL(R8), INTENT(IN) :: EtoliationCoeff
  real(r8), intent(in) :: WFNS
  REAL(R8), INTENT(IN) :: ALLOCL
  integer :: MXNOD,MNNOD,NE,KK,K
  real(r8) :: GNOD,ALLOCS,GSSL,GROS
  REAL(R8) :: GrowthElms(NumPlantChemElms)
  REAL(R8) :: SSL
  
  associate(                                                     &
    NumCogrowthNode_pft     => plt_morph%NumCogrowthNode_pft,    &
    PetoleLensNode_brch     => plt_morph%PetoleLensNode_brch,    &
    PetoLen2Mass_pft        => plt_morph%PetoLen2Mass_pft,       &
    SinePetioleAngle_pft    => plt_morph%SinePetioleAngle_pft,   &
    PetioleElmntNode_brch   => plt_biom%PetioleElmntNode_brch,   &
    PetoleProteinCNode_brch => plt_biom%PetoleProteinCNode_brch, &
    LeafElmntNode_brch      => plt_biom%LeafElmntNode_brch,      &
    ZERO4LeafVar_pft        => plt_biom%ZERO4LeafVar_pft,        &
    rCNNonstRemob_pft       => plt_allom%rCNNonstRemob_pft,      &
    rCPNonstRemob_pft       => plt_allom%rCPNonstRemob_pft,      &
    FracGroth2Node_pft      => plt_allom%FracGroth2Node_pft,     &
    PlantPopulation_pft     => plt_site%PlantPopulation_pft,     &
    KHiestGroLeafNode_brch  => plt_pheno%KHiestGroLeafNode_brch  &
  )

  IF(GrowthPetiole(ielmc).GT.0.0_r8)THEN
    MXNOD=KHiestGroLeafNode_brch(NB,NZ)
    MNNOD=MAX(MinNodeNum,MXNOD-NumCogrowthNode_pft(NZ)+1)
    MXNOD=MAX(MXNOD,MNNOD)
    GNOD=MXNOD-MNNOD+1
    ALLOCS=1.0_r8/GNOD
    DO NE=1,NumPlantChemElms
      GrowthElms(NE)=ALLOCS*GrowthPetiole(NE)
    ENDDO
    GSSL=ALLOCL*FracGroth2Node_pft(NZ)*NumCogrowthNode_pft(NZ)
!
!       GROWTH AT EACH CURRENT NODE
!
!       PetioleElmntNode_brch,WGSHN,WGSHP,PetoleProteinCNode_brch=node petiole C,N,P,protein mass
!       GRO,GrowthElms(ielmn),GrowthElms(ielmp)=petiole C,N,P growth at each node
!       CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios from startq.f
!
    D505: DO KK=MNNOD,MXNOD
      K=pMOD(KK,MaxNodesPerBranch1)

      DO NE=1,NumPlantChemElms
        PetioleElmntNode_brch(NE,K,NB,NZ)=PetioleElmntNode_brch(NE,K,NB,NZ)+GrowthElms(NE)
      ENDDO
      PetoleProteinCNode_brch(K,NB,NZ)=PetoleProteinCNode_brch(K,NB,NZ) &
        +AMIN1(GrowthElms(ielmn)*rCNNonstRemob_pft(NZ),GrowthElms(ielmp)*rCPNonstRemob_pft(NZ))
!
!           SPECIFIC SHEATH OR PETIOLE LENGTH FUNCTION OF CURRENT MASS
!           AT EACH NODE
    !
    !   SSL=specific length of petiole growth
    !   ETOL=coefficient for etoliation effects on expansion,extension
    !   PetoLen2Mass_pft=growth in petiole length vs mass from PFT file
    !   SSL2=parameter for calculating petiole extension
    !   PetioleElmntNode_brch=petiole C mass
    !   PP=PFT population
    !   GSSL=allocation of petiole length growth to each node
    !   WFNS=turgor expansion,extension function
    !   GROS,GRO=petiole length,mass growth
    !   PetoleLensNode_brch=petiole length
!
      IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.0.0_r8)THEN
        SSL=EtoliationCoeff*PetoLen2Mass_pft(NZ)*(AMAX1(ZERO4LeafVar_pft(NZ) &
          ,PetioleElmntNode_brch(ielmc,K,NB,NZ))/(PlantPopulation_pft(NZ)*GSSL))**SSL2*WFNS
        GROS=GrowthElms(ielmc)/PlantPopulation_pft(NZ)*SSL
        PetoleLensNode_brch(K,NB,NZ)=PetoleLensNode_brch(K,NB,NZ)+GROS*SinePetioleAngle_pft(NZ)
      ENDIF
    ENDDO D505
  ENDIF
  end associate
  end subroutine GrowPetioleOnBranch  
!------------------------------------------------------------------------------------------

  subroutine GrowStalkOnBranch(NZ,NB,GrowthStalk,ETOL)
  implicit none
  integer, intent(in) :: NZ,NB
  real(r8), INTENT(IN) :: GrowthStalk(NumPlantChemElms)
  REAL(R8), INTENT(IN) :: ETOL
  
  character(len=*), parameter :: subname='GrowStalkOnBranch'  
  REAL(R8) :: GrowthElms(NumPlantChemElms)
  integer :: NN,MXNOD,MNNOD,K1,K2,NE,KX,KK
  REAL(R8) :: ALLOCN,GNOD,SpecLenStalkGrowth
  real(r8) :: StalkLenGrowth,StalkAveStrutC_brch

  associate(                                                        &
    NumCogrowthNode_pft      => plt_morph%NumCogrowthNode_pft,      &
    NodeLenPergC             => plt_morph%NodeLenPergC,             &
    SineBranchAngle_pft      => plt_morph%SineBranchAngle_pft,      &
    InternodeHeightDead_brch => plt_morph%InternodeHeightDead_brch, &
    LiveInterNodeHight_brch  => plt_morph%LiveInterNodeHight_brch,  &
    InternodeStrutElms_brch  => plt_biom%InternodeStrutElms_brch,   &
    StalkStrutElms_brch      => plt_biom%StalkStrutElms_brch,       &
    PlantPopulation_pft      => plt_site%PlantPopulation_pft,       &
    iPlantCalendar_brch      => plt_pheno%iPlantCalendar_brch,      &
    KHiestGroLeafNode_brch   => plt_pheno%KHiestGroLeafNode_brch    &
  )
  !plant hasnot emerged
  IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).EQ.0)THEN  
    NN=0
  !plant emerged  
  ELSE
    NN=1
  ENDIF

  MXNOD = KHiestGroLeafNode_brch(NB,NZ)
  MNNOD = MAX(MIN(NN,MAX(NN,MXNOD-NumCogrowthNode_pft(NZ))),KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+2)
  MXNOD = MAX(MXNOD,MNNOD)

  IF(GrowthStalk(ielmc).GT.0.0_r8)THEN
    GNOD=MXNOD-MNNOD+1
    ALLOCN=1.0_r8/GNOD
    DO NE=1,NumPlantChemElms
      GrowthElms(NE)=ALLOCN*GrowthStalk(NE)
    ENDDO
!
!     SPECIFIC INTERNODE LENGTH FUNCTION OF CURRENT STALK MASS
!     AT EACH NODE
!
!     SpecLenStalkGrowth=specific length of stalk growth
!     ETOL=coefficient for etoliation effects on expansion,extension
!     NodeLenPergC=growth in stalk length vs mass from PFT file
!     SNL2=parameter for calculating stalk extension
!     WTSKB=stalk C mass
!     PP=PFT population
!     StalkLenGrowth,GRO=stalk length,mass growth
!
    StalkAveStrutC_brch = StalkStrutElms_brch(ielmc,NB,NZ)/PlantPopulation_pft(NZ)
    SpecLenStalkGrowth  = ETOL*NodeLenPergC(NZ)*(StalkAveStrutC_brch)**SNL2
    StalkLenGrowth      = GrowthElms(ielmc)/PlantPopulation_pft(NZ)*SpecLenStalkGrowth
    KX                  = pMOD(MNNOD-1,MaxNodesPerBranch1)

!
!     GROWTH AT EACH CURRENT NODE
!
!     InternodeStrutElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
!     GRO,GrowthElms(ielmn),GrowthElms(ielmp)=stalk C,N,P growth at each node
!     InternodeHeightDead_brch,LiveInterNodeHight_brch=stalk height,stalk internode length
!     SineBranchAngle_pft=sine of stalk angle from horizontal from PFT file
!
    D510: DO KK=MNNOD,MXNOD
      K1 = pMOD(KK,MaxNodesPerBranch1)
      K2 = pMOD(KK-1,MaxNodesPerBranch1)
      DO NE = 1, NumPlantChemElms
        InternodeStrutElms_brch(NE,K1,NB,NZ)=InternodeStrutElms_brch(NE,K1,NB,NZ)+GrowthElms(NE)
      ENDDO
      InternodeHeightDead_brch(K1,NB,NZ)=InternodeHeightDead_brch(K1,NB,NZ)+StalkLenGrowth*SineBranchAngle_pft(NZ)
      IF(K1.NE.0)THEN
        LiveInterNodeHight_brch(K1,NB,NZ)=InternodeHeightDead_brch(K1,NB,NZ)+LiveInterNodeHight_brch(K2,NB,NZ)
      ELSE
        LiveInterNodeHight_brch(K1,NB,NZ)=InternodeHeightDead_brch(K1,NB,NZ)
      ENDIF
    ENDDO D510
  ENDIF
  END associate  
  end subroutine GrowStalkOnBranch

!------------------------------------------------------------------------------------------

  subroutine RemobilizeBranch(NZ,NB,BegRemoblize,LRemob_brch,RCCC,RCCN,RCCP,RMxess_brch)
  implicit none
  integer, intent(in) :: NZ,NB
  integer, intent(in) :: LRemob_brch,BegRemoblize
  real(r8), intent(in) :: RCCC,RCCN,RCCP
  real(r8), intent(inout) :: RMxess_brch
  integer :: NumRemobLeafNodes  
  real(r8) :: RespSenesTot_brch
  real(r8) :: SenesFrac  
  integer  :: NBZ(MaxNumBranches)
  real(r8) :: XFRE(1:NumPlantChemElms)
  REAL(R8) :: RCO2V,RespSenesPhenol_brch
  integer :: NBY,NBX,NBL
  integer  :: NBK,NE

  associate(                                                              &
    fTCanopyGroth_pft         =>  plt_pheno%fTCanopyGroth_pft           , &  
    iPlantTurnoverPattern_pft =>  plt_pheno%iPlantTurnoverPattern_pft   , &    
    HoursDoingRemob_brch      =>  plt_pheno%HoursDoingRemob_brch        , &
    KHiestGroLeafNode_brch    =>  plt_pheno%KHiestGroLeafNode_brch      , &    
    iPlantPhenolPattern_pft   =>  plt_pheno%iPlantPhenolPattern_pft     , &   
    iPlantRootProfile_pft     =>  plt_pheno%iPlantRootProfile_pft       , &     
    iPlantBranchState_brch    =>  plt_pheno%iPlantBranchState_brch      , &
    KLowestGroLeafNode_brch   =>  plt_pheno%KLowestGroLeafNode_brch     , &    
    NumOfBranches_pft         =>  plt_morph%NumOfBranches_pft           , &    
    MainBranchNum_pft         =>  plt_morph%MainBranchNum_pft           , &         
    BranchNumber_brch         =>  plt_morph%BranchNumber_brch           , &        
    CanopyNonstElms_brch      =>  plt_biom%CanopyNonstElms_brch         , &   
    LeafPetolBiomassC_brch    =>  plt_biom%LeafPetolBiomassC_brch       , &        
    StalkRsrvElms_brch        =>  plt_biom%StalkRsrvElms_brch           , &      
    ZERO4Groth_pft            =>  plt_biom%ZERO4Groth_pft                 &  
  ) 
  !     REMOBILIZATION OF STALK RESERVE C,N,P IF GROWTH RESPIRATION < 0

!
!     RMxess_brch=excess maintenance respiration
!     WTRSVB=stalk reserve C mass
!     RCO2V=remobilization of stalk reserve C
!     VMXC=rate constant for nonstructural C oxidation in respiration
!     fTCanopyGroth_pft=temperature function for canopy growth
!
  IF(BegRemoblize.EQ.ifalse)THEN
    IF(RMxess_brch.GT.0.0_r8 .AND. StalkRsrvElms_brch(ielmc,NB,NZ).GT.0.0_r8)THEN
      RCO2V                           = AMIN1(RMxess_brch,VMXC*StalkRsrvElms_brch(ielmc,NB,NZ)*fTCanopyGroth_pft(NZ))
      StalkRsrvElms_brch(ielmc,NB,NZ) = StalkRsrvElms_brch(ielmc,NB,NZ)-RCO2V
      RMxess_brch                     = RMxess_brch-RCO2V
    ENDIF
  ENDIF
!
!       TOTAL REMOBILIZATION = GROWTH RESPIRATION < 0 + DECIDUOUS LEAF
!       FALL DURING AUTUMN + REMOBILZATION DURING GRAIN FILL IN ANNUALS
!
!       iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!       LRemob_brch,BegRemoblize=remobilization flags
!       RespSenesPhenol_brch=phenologically-driven respiration senescence during late-season
!       RateConst4ShootSeaStoreNonstXfer=rate constant for plant-storage nonstructural C,N,P exchange
!       iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+petiole,2=none,3=between 1,2
!       LeafPetolBiomassC_brch=leaf+petiole mass
!       FLGZ=control rate of remobilization
!       Hours4FullSenes=number of hours until full senescence after physl maturity
!       RespSenesTot_brch=total senescence respiration
!       KHiestGroLeafNode_brch,KLowestGroLeafNode_brch=integer of highest,lowest leaf number currently growing
!       KSNC=number of nodes undergoing remobilization
!       SenesFrac=ratio of phenologically-driven vs total senescence respiration
!
  IF(BegRemoblize.EQ.itrue .AND. LRemob_brch.EQ.itrue .AND. iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
    RespSenesPhenol_brch=RateConst4ShootSeaStoreNonstXfer(iPlantTurnoverPattern_pft(NZ))*LeafPetolBiomassC_brch(NB,NZ) &
      *AMIN1(1.0_r8,HoursDoingRemob_brch(NB,NZ)/Hours4FullSenes)
  ELSE
    RespSenesPhenol_brch=0._r8
  ENDIF
  RespSenesTot_brch=RMxess_brch+RespSenesPhenol_brch

  IF(RespSenesTot_brch.GT.ZERO4Groth_pft(NZ))THEN
    SenesFrac         = RespSenesPhenol_brch/RespSenesTot_brch
    NumRemobLeafNodes = INT(0.5_r8*(KHiestGroLeafNode_brch(NB,NZ)-KLowestGroLeafNode_brch(NB,NZ)))+1
!
!     TRANSFER NON-STRUCTURAL C,N,P FROM BRANCHES TO MAIN STEM
!     IF MAIN STEM POOLS ARE DEPLETED
!
!     iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!     XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
!
    IF(iPlantTurnoverPattern_pft(NZ).NE.0 .AND. is_plant_treelike(iPlantRootProfile_pft(NZ)) &
      .AND.NB.EQ.MainBranchNum_pft(NZ) .AND. isclose(SenesFrac,0._r8))THEN
      NBY=0
      D584: DO NBL=1,NumOfBranches_pft(NZ)
        NBZ(NBL)=0
      ENDDO D584

      D586: DO NBL=1,NumOfBranches_pft(NZ)
        NBX=KHiestGroLeafNode_brch(NB,NZ)
        D585: DO NBK=1,NumOfBranches_pft(NZ)
          IF(iPlantBranchState_brch(NBK,NZ).EQ.iLive .AND. NBK.NE.MainBranchNum_pft(NZ) &
            .AND. BranchNumber_brch(NBK,NZ).LT.NBX .AND. BranchNumber_brch(NBK,NZ).GT.NBY)THEN
            NBZ(NBL) = NBK
            NBX      = BranchNumber_brch(NBK,NZ)
          ENDIF
        ENDDO D585
        IF(NBZ(NBL).NE.0)THEN
          NBY=BranchNumber_brch(NBZ(NBL),NZ)
        ENDIF
      ENDDO D586

      D580: DO NBL=1,NumOfBranches_pft(NZ)
        IF(NBZ(NBL).NE.0)THEN          
          IF(BranchNumber_brch(NBZ(NBL),NZ).NE.MainBranchNum_pft(NZ))THEN
            IF(CanopyNonstElms_brch(ielmc,NBZ(NBL),NZ).GT.0._r8)THEN
              XFRE(ielmc)=1.0E-02_r8*AMIN1(RespSenesTot_brch,CanopyNonstElms_brch(ielmc,NBZ(NBL),NZ))
              DO NE=2,NumPlantChemElms
                XFRE(NE)=XFRE(ielmc)*CanopyNonstElms_brch(NE,NBZ(NBL),NZ)/CanopyNonstElms_brch(ielmc,NBZ(NBL),NZ)
              ENDDO
            ELSE
              XFRE(ielmc)=0._r8
              DO NE=2,NumPlantChemElms
                XFRE(NE)=1.0E-02_r8*CanopyNonstElms_brch(NE,NBZ(NBL),NZ)
              ENDDO
            ENDIF
            DO NE=1,NumPlantChemElms
              CanopyNonstElms_brch(NE,NBZ(NBL),NZ)=CanopyNonstElms_brch(NE,NBZ(NBL),NZ)-XFRE(NE)
            ENDDO
            CanopyNonstElms_brch(ielmc,MainBranchNum_pft(NZ),NZ)=CanopyNonstElms_brch(ielmc,MainBranchNum_pft(NZ),NZ) &
              +XFRE(ielmc)*SenesFrac
            DO NE=2,NumPlantChemElms
              CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)=CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)+XFRE(NE)
            ENDDO
            RespSenesTot_brch=RespSenesTot_brch-XFRE(ielmc)
            IF(RespSenesTot_brch.LE.0.0_r8)exit
          ENDIF
        ENDIF
      ENDDO D580
    ENDIF
    IF(RespSenesTot_brch.GT.0.0_r8) &
      call RemobilizeLeafLayers(NumRemobLeafNodes,NB,nz,RespSenesTot_brch,RCCC,RCCN,RCCP,SenesFrac)
  ENDIF
  end associate
  end subroutine RemobilizeBranch
!------------------------------------------------------------------------------------------
  subroutine SenescenceBranch(NZ,NB,RCCC,RCCN,RCCP)
  implicit none
  integer, intent(in) :: NZ,NB
  real(r8), intent(in) :: RCCC,RCCN,RCCP

  character(len=*), parameter :: subname='SenescenceBranch'
  INTEGER :: K,KMinGroingLeafNodeNum,M,NE
  real(r8) :: FSNC
  real(r8) :: dRemoblE,FracRemobAsLeaf
  real(r8) :: FSNCS
  
  associate(                                                              &
    icwood                      => pltpar%icwood,                         &
    ifoliar                     => pltpar%ifoliar,                        &
    inonfoliar                  => pltpar%inonfoliar,                     &
    k_woody_litr                => pltpar%k_woody_litr,                   &
    k_fine_litr                 => pltpar%k_fine_litr,                    &
    ElmAllocmat4Litr            => plt_soilchem%ElmAllocmat4Litr,         &
    FracShootStalkElmAlloc2Litr => plt_allom%FracShootStalkElmAlloc2Litr, &
    FracShootLeafElmAlloc2Litr  => plt_allom%FracShootLeafElmAlloc2Litr,  &
    rCNNonstRemob_pft           => plt_allom%rCNNonstRemob_pft,           &
    rCPNonstRemob_pft           => plt_allom%rCPNonstRemob_pft,           &
    LitrfalStrutElms_pvr        => plt_bgcr%LitrfalStrutElms_pvr,         &
    SenecStalkStrutElms_brch    => plt_biom%SenecStalkStrutElms_brch,     &
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft,               &
    PetioleElmntNode_brch       => plt_biom%PetioleElmntNode_brch,        &
    LeafProteinCNode_brch       => plt_biom%LeafProteinCNode_brch,        &
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch,         &
    LeafStrutElms_brch          => plt_biom%LeafStrutElms_brch,           &
    PetoleProteinCNode_brch     => plt_biom%PetoleProteinCNode_brch,      &
    PetioleChemElmRemob_brch    => plt_biom%PetioleChemElmRemob_brch,     &
    InternodeStrutElms_brch     => plt_biom%InternodeStrutElms_brch,      &
    LeafChemElmRemob_brch       => plt_biom%LeafChemElmRemob_brch,        &
    LeafElmntNode_brch          => plt_biom%LeafElmntNode_brch,           &
    PetoleStrutElms_brch        => plt_biom%PetoleStrutElms_brch,         &
    CanPBranchHeight            => plt_morph%CanPBranchHeight,            &
    InternodeHeightDead_brch    => plt_morph%InternodeHeightDead_brch,    &
    LeafNodeArea_brch           => plt_morph%LeafNodeArea_brch,           &
    LeafAreaDying_brch          => plt_morph%LeafAreaDying_brch,          &
    LeafAreaLive_brch           => plt_morph%LeafAreaLive_brch,           &
    PetoleLensNode_brch         => plt_morph%PetoleLensNode_brch,         &
    doSenescence_brch           => plt_pheno%doSenescence_brch,           &
    doRemobilization_brch       => plt_pheno%doRemobilization_brch,       &
    PetioleChemElmRemobFlx_brch => plt_pheno%PetioleChemElmRemobFlx_brch, &
    RefLeafAppearRate_pft       => plt_pheno%RefLeafAppearRate_pft,       &
    KHiestGroLeafNode_brch      => plt_pheno%KHiestGroLeafNode_brch,      &
    LeafElmntRemobFlx_brch      => plt_pheno%LeafElmntRemobFlx_brch,      &
    fTCanopyGroth_pft           => plt_pheno%fTCanopyGroth_pft            &
  )    
  call PrintInfo('beg '//subname)
  IF(doSenescence_brch(NB,NZ).EQ.itrue)THEN
    KMinGroingLeafNodeNum=MAX(0,KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+1)
    IF(KMinGroingLeafNodeNum.GT.0)THEN
      K    = pMOD(KMinGroingLeafNodeNum,MaxNodesPerBranch1)
      FSNC = fTCanopyGroth_pft(NZ)*RefLeafAppearRate_pft(NZ)
!
      !   REMOBILIZATION OF LEAF C,N,P ALSO DEPENDS ON STRUCTURAL C:N:P
      !
      !   doRemobilization_brch=flag for remobilization
      !   WGLF,WGLFN,WGLFP=node leaf C,N,P mass
      !   LeafNodeArea_brch=node leaf area
      !   ElmntRemobFallingLeaf(ielmc)X,ElmntRemobFallingLeaf(ielmn)X,ElmntRemobFallingLeaf(ielmp)X=remobilization of C,N,P from senescing leaf
!
      IF(doRemobilization_brch(NB,NZ).EQ.itrue)THEN
        DO NE=1,NumPlantChemElms
          LeafChemElmRemob_brch(NE,NB,NZ)=AZMAX1(LeafElmntNode_brch(NE,K,NB,NZ))
        ENDDO
        LeafAreaDying_brch(NB,NZ)=AZMAX1(LeafNodeArea_brch(K,NB,NZ))
        IF(LeafChemElmRemob_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
          LeafElmntRemobFlx_brch(ielmc,NB,NZ) = LeafChemElmRemob_brch(ielmc,NB,NZ)*RCCC
          LeafElmntRemobFlx_brch(ielmn,NB,NZ) = LeafChemElmRemob_brch(ielmn,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
          LeafElmntRemobFlx_brch(ielmp,NB,NZ) = LeafChemElmRemob_brch(ielmp,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
        ELSE
          LeafElmntRemobFlx_brch(1:NumPlantChemElms,NB,NZ)=0._r8
        ENDIF
      ENDIF
!
  !       FRACTION OF CURRENT LEAF TO BE REMOBILIZED
  !
  !       FSNC,FSNCL=fraction of lowest leaf to be remobilized
!
      IF(FSNC*LeafChemElmRemob_brch(ielmc,NB,NZ).GT.LeafElmntNode_brch(ielmc,K,NB,NZ) &
        .AND.LeafChemElmRemob_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
        FracRemobAsLeaf = AZMAX1(LeafElmntNode_brch(ielmc,K,NB,NZ)/LeafChemElmRemob_brch(ielmc,NB,NZ))
      ELSE
        FracRemobAsLeaf=FSNC
      ENDIF
!
  !       NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
  !       TO FRACTIONS SET IN 'STARTQ'
  !
  !       CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
  !       CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
  !       FSNCL=fraction of lowest leaf to be remobilized
  !       ElmntRemobFallingLeaf(ielmc)X,ElmntRemobFallingLeaf(ielmn)X,ElmntRemobFallingLeaf(ielmp)X=remobilization of C,N,P from senescing leaf
  !       WGLFX,WGLFNX,WGLFPX=senescing leaf C,N,P mass
  !       FWODB=C woody fraction in other organs:0=woody,1=non-woody
  !       FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!
      DO NE=1,NumPlantChemElms
        D6300: DO M=1,jsken
          dRemoblE=FracRemobAsLeaf*AZMAX1(LeafChemElmRemob_brch(NE,NB,NZ)-LeafElmntRemobFlx_brch(NE,NB,NZ)) 
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,icwood,M,NZ)*dRemoblE*FracShootStalkElmAlloc2Litr(NE,k_woody_litr)
            
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,ifoliar,M,NZ)*dRemoblE*FracShootStalkElmAlloc2Litr(NE,k_fine_litr)
        ENDDO D6300
      ENDDO
!
!       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
!
!       FSNCL=fraction of lowest leaf to be remobilized
!       LeafAreaLive_brch,LeafAreaDying_brch=branch living,senescing leaf area
!       WTLFB,WTLFBN,WTLFBP,WGLFX,WGLFNX,WGLFPX=C,N,P mass in living,senescing leaf
!       LeafProteinCNode_brch=leaf protein mass
!       CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios from startq.f
!       CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!       ElmntRemobFallingLeaf(ielmc)X,ElmntRemobFallingLeaf(ielmn)X,ElmntRemobFallingLeaf(ielmp)X=remobilization of C,N,P from senescing leaf
!
      LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)-FracRemobAsLeaf*LeafAreaDying_brch(NB,NZ)
      DO NE=1,NumPlantChemElms
        LeafStrutElms_brch(NE,NB,NZ)   = LeafStrutElms_brch(NE,NB,NZ)-FracRemobAsLeaf*LeafChemElmRemob_brch(NE,NB,NZ)
        LeafElmntNode_brch(NE,K,NB,NZ) = LeafElmntNode_brch(NE,K,NB,NZ)-FracRemobAsLeaf*LeafChemElmRemob_brch(NE,NB,NZ)
        CanopyNonstElms_brch(NE,NB,NZ) = CanopyNonstElms_brch(NE,NB,NZ)+FracRemobAsLeaf*LeafElmntRemobFlx_brch(NE,NB,NZ)
      ENDDO
      LeafNodeArea_brch(K,NB,NZ)     = LeafNodeArea_brch(K,NB,NZ)-FracRemobAsLeaf*LeafAreaDying_brch(NB,NZ)
      LeafProteinCNode_brch(K,NB,NZ) = AZMAX1(LeafProteinCNode_brch(K,NB,NZ) &
        -FracRemobAsLeaf*AMAX1(LeafChemElmRemob_brch(ielmn,NB,NZ)*rCNNonstRemob_pft(NZ) &
        ,LeafChemElmRemob_brch(ielmp,NB,NZ)*rCPNonstRemob_pft(NZ)))
!
!       REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P ALSO DEPENDS ON
!       STRUCTURAL C:N:P
!
!       doRemobilization_brch=flag for remobilization
!       PetioleElmntNode_brch,WGSHN,WGSHP=node petiole C,N,P mass
!       CanPBranchHeight=petiole length
!       RCES(ielmc)X,RCES(ielmn)X,RCES(ielmp)X=remobilization of C,N,P from senescing petiole
!
      IF(doRemobilization_brch(NB,NZ).EQ.itrue)THEN
        DO NE=1,NumPlantChemElms
          PetioleChemElmRemob_brch(NE,NB,NZ)=AZMAX1(PetioleElmntNode_brch(NE,K,NB,NZ))
        ENDDO

        CanPBranchHeight(NB,NZ)=AZMAX1(PetoleLensNode_brch(K,NB,NZ))
        IF(PetioleChemElmRemob_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
          PetioleChemElmRemobFlx_brch(ielmc,NB,NZ) = RCCC*PetioleChemElmRemob_brch(ielmc,NB,NZ)
          PetioleChemElmRemobFlx_brch(ielmn,NB,NZ) = PetioleChemElmRemob_brch(ielmn,NB,NZ) &
            *(RCCN+(1.0_r8-RCCN)*PetioleChemElmRemobFlx_brch(ielmc,NB,NZ)/PetioleChemElmRemob_brch(ielmc,NB,NZ))
          PetioleChemElmRemobFlx_brch(ielmp,NB,NZ)=PetioleChemElmRemob_brch(ielmp,NB,NZ) &
            *(RCCP+(1.0_r8-RCCP)*PetioleChemElmRemobFlx_brch(ielmc,NB,NZ)/PetioleChemElmRemob_brch(ielmc,NB,NZ))
        ELSE
          PetioleChemElmRemobFlx_brch(1:NumPlantChemElms,NB,NZ)=0._r8
        ENDIF
        
        DO NE=1,NumPlantChemElms
          SenecStalkStrutElms_brch(NE,NB,NZ)=SenecStalkStrutElms_brch(NE,NB,NZ)+&
            InternodeStrutElms_brch(NE,K,NB,NZ)
        ENDDO
        InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8
        InternodeHeightDead_brch(K,NB,NZ)                  = 0._r8
      ENDIF
!
!       FRACTION OF CURRENT SHEATH TO BE REMOBILIZED
!
!       FSNCS=fraction of lowest petiole to be remobilized
!
      IF(FSNC*PetioleChemElmRemob_brch(ielmc,NB,NZ).GT.PetioleElmntNode_brch(ielmc,K,NB,NZ) &
        .AND.PetioleChemElmRemob_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
        FSNCS=AZMAX1(PetioleElmntNode_brch(ielmc,K,NB,NZ)/PetioleChemElmRemob_brch(ielmc,NB,NZ))
      ELSE
        FSNCS=FSNC
      ENDIF
!
!       NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
!       TO FRACTIONS SET IN 'STARTQ'
!
!       CSNC,ZSNC,PSNC=C,N,P LitrFall from senescence
!       CFOPC,CFOPN,CFOPC=fraction of LitrFall C,N,P allocated to litter components
!       FSNCS=fraction of lowest petiole to be remobilized
!       RCES(ielmc)X,RCES(ielmn)X,RCES(ielmp)X=remobilization of C,N,P from senescing petiole
!       WGSHX,WGSHNX,WGSHPX=senescing petiole C,N,P mass
!       FWODB=C woody fraction in other organs:0=woody,1=non-woody
!       FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
!
      
      D6305: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          dRemoblE                                     = FSNCS*AZMAX1(PetioleChemElmRemob_brch(NE,NB,NZ)-PetioleChemElmRemobFlx_brch(NE,NB,NZ))
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) = LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,icwood,M,NZ)*dRemoblE*FracShootLeafElmAlloc2Litr(NE,k_woody_litr)

          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +ElmAllocmat4Litr(NE,inonfoliar,M,NZ)*dRemoblE*FracShootLeafElmAlloc2Litr(NE,k_fine_litr)
        ENDDO    
      ENDDO D6305      
!
!       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
!
!       FSNCS=fraction of lowest petiole to be remobilized
!       PetoleLensNode_brch,CanPBranchHeight=living,senescing petiole length
!       WTSHB,WTSHBN,WTSHBP,WGSHEX,WGSHNX,WGSHPX=C,N,P mass in living,senescing petiole
!       PetoleProteinCNode_brch=petiole protein mass
!       CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios from startq.f
!       CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!       RCES(ielmc)X,RCES(ielmn)X,RCES(ielmp)X=remobilization of C,N,P from senescing petiole
!
      DO NE=1,NumPlantChemElms
        PetoleStrutElms_brch(NE,NB,NZ)    = PetoleStrutElms_brch(NE,NB,NZ)-FSNCS*PetioleChemElmRemob_brch(NE,NB,NZ)
        PetioleElmntNode_brch(NE,K,NB,NZ) = PetioleElmntNode_brch(NE,K,NB,NZ)-FSNCS*PetioleChemElmRemob_brch(NE,NB,NZ)
        CanopyNonstElms_brch(NE,NB,NZ)    = CanopyNonstElms_brch(NE,NB,NZ)+FSNCS*PetioleChemElmRemobFlx_brch(NE,NB,NZ)
      ENDDO
      PetoleLensNode_brch(K,NB,NZ)     = PetoleLensNode_brch(K,NB,NZ)-FSNCS*CanPBranchHeight(NB,NZ)
      PetoleProteinCNode_brch(K,NB,NZ) = AZMAX1(PetoleProteinCNode_brch(K,NB,NZ) &
        -FSNCS*AMAX1(PetioleChemElmRemob_brch(ielmn,NB,NZ)*rCNNonstRemob_pft(NZ) &
        ,PetioleChemElmRemob_brch(ielmp,NB,NZ)*rCPNonstRemob_pft(NZ)))

    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  END associate
  end subroutine SenescenceBranch  

end module PlantBranchMod
