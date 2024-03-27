module PlantBranchMod

! Description:
! module for plant biological transformations
  use minimathmod, only : isclose,safe_adb,AZMAX1
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod , only : etimer 
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
  integer, parameter :: ibrch_petiole=2
  integer, parameter :: ibrch_stalk=3
  integer, parameter :: ibrch_reserve=4
  integer, parameter :: ibrch_husk=5
  integer, parameter :: ibrch_ear=6
  integer, parameter :: ibrch_grain=7
  integer, save :: iter=0
  integer :: II,JJ
  contains
!------------------------------------------------------------------------------------------

  subroutine GrowOneBranch(I,J,NB,NZ,TFN6_vr,CanopyHeight_copy,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,TFN5,WFNG,&
    Stomata_Activity,WFNS,WFNSG,PTRT,CanopyN2Fix_pft,BegRemoblize)
  implicit none
  integer, intent(in)  :: I,J,NB,NZ
  REAL(R8), INTENT(IN) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: CanopyHeight_copy(JP1)
  real(r8), intent(in) :: CNLFW,CPLFW
  real(r8), intent(in) :: CNSHW,CPSHW
  real(r8), intent(in) :: CNRTW,CPRTW
  real(r8), intent(in) :: TFN5,WFNG
  real(r8), intent(in) :: Stomata_Activity
  real(r8), intent(in) :: WFNS,WFNSG
  real(r8), intent(inout) :: CanopyN2Fix_pft(JP1)
  real(r8), intent(out) :: PTRT
  integer, intent(out) :: BegRemoblize
  real(r8) :: DMSHD
  integer  :: K,KNOD,KK,K1,K2,KMinGroingLeafNodeNum
  integer  :: NN,M,N,NNOD1,NE
  real(r8) :: ZPOOLD,XFRN1,XFRP1
  REAL(R8) :: CH2O3(pltpar%MaxNodesPerBranch1),CH2O4(pltpar%MaxNodesPerBranch1)
  integer  :: IFLGY
  REAL(R8) :: PART(pltpar%NumOfPlantMorphUnits)
  real(r8) :: ALLOCL,ALLOCS
  REAL(R8) :: ALLOCN
  REAL(R8) :: CNPG
  real(r8) :: CCE
  real(r8) :: CCC,CNC,CPC
  REAL(R8) :: cpoolt
  real(r8) :: DMLFB
  real(r8) :: DMSHB
  real(r8) :: DMSHT
  real(r8) :: CNLFB
  real(r8) :: CPLFB
  real(r8) :: ETOL
  real(r8) :: GrowthGrain(NumPlantChemElms)
  real(r8) :: GrowthLeaf(NumPlantChemElms)
  real(r8) :: GrowthReserve(NumPlantChemElms)
  real(r8) :: GrowthHusk(NumPlantChemElms)
  real(r8) :: GrowthEar(NumPlantChemElms)
  real(r8) :: GrowthPetiole(NumPlantChemElms)
  real(r8) :: GrowthStalk(NumPlantChemElms)
  real(r8) :: GNOD
  REAL(R8) :: GrowthChemElmt(NumPlantChemElms)
  real(r8) :: GSSL,GROS
  real(r8) :: PPOOLD,RCO2C
  REAL(R8) :: RMNCS
  real(r8) :: RMxess_brch
  real(r8) :: RCCC,RCCN,RCCP
  real(r8) :: SSL
  real(r8) :: CNSHB,CPSHB
  real(r8) :: CNLFM,CPLFM
  real(r8) :: CNLFX,CPLFX,CNSHX,CPSHX
  real(r8) :: NonStructalC4Growth_brch
  real(r8) :: CNRDM,Rauto4Nassim_brch
  real(r8) :: ShootStructN
  real(r8) :: XFRE(1:NumPlantChemElms)
! begin_execution
  associate(                                                    &
    LeafBiomGrowthYield          =>  plt_allom%LeafBiomGrowthYield    , &
    PetioleBiomGrowthYield       =>  plt_allom%PetioleBiomGrowthYield   , &
    EarBiomGrowthYield           =>  plt_allom%EarBiomGrowthYield   , &
    rNCEar_pft                   =>  plt_allom%rNCEar_pft    , &
    rPCReserve_pft               =>  plt_allom%rPCReserve_pft    , &
    rPCHusk_pft                  =>  plt_allom%rPCHusk_pft    , &
    rPCStalk_pft                 =>  plt_allom%rPCStalk_pft    , &
    rNCStalk_pft                 =>  plt_allom%rNCStalk_pft   , &
    StalkBiomGrowthYield         =>  plt_allom%StalkBiomGrowthYield   , &
    HuskBiomGrowthYield          =>  plt_allom%HuskBiomGrowthYield   , &
    rNCHusk_pft                  =>  plt_allom%rNCHusk_pft    , &
    rPCEar_pft                   =>  plt_allom%rPCEar_pft    , &
    ReserveBiomGrowthYield       =>  plt_allom%ReserveBiomGrowthYield   , &
    GrainBiomGrowthYield         =>  plt_allom%GrainBiomGrowthYield    , &
    rCNNonstructRemob_pft        =>  plt_allom%rCNNonstructRemob_pft    , &
    rCPNonstructRemob_pft        =>  plt_allom%rCPNonstructRemob_pft     , &
    rNCReserve_pft               =>  plt_allom%rNCReserve_pft    , &
    RootBiomGrowthYield          =>  plt_allom%RootBiomGrowthYield     , &
    FNOD                         =>  plt_allom%FNOD     , &
    StalkStrutElms_brch          =>  plt_biom%StalkStrutElms_brch   , &
    StalkRsrvElms_brch           =>  plt_biom%StalkRsrvElms_brch   , &
    HuskStrutElms_brch           =>  plt_biom%HuskStrutElms_brch   , &
    EarStrutElms_brch            =>  plt_biom%EarStrutElms_brch   , &
    StalkBiomassC_brch           =>  plt_biom%StalkBiomassC_brch    , &
    CanopyNonstElms_brch         =>  plt_biom%CanopyNonstElms_brch    , &
    LeafPetolBiomassC_brch       =>  plt_biom%LeafPetolBiomassC_brch     , &
    PetoleStrutElms_brch         =>  plt_biom%PetoleStrutElms_brch  , &
    LeafStrutElms_brch           =>  plt_biom%LeafStrutElms_brch   , &
    LeafProteinCNode_brch        =>  plt_biom%LeafProteinCNode_brch      , &
    LeafElmntNode_brch           =>  plt_biom%LeafElmntNode_brch     , &
    LeafPetoNonstElmConc_brch    =>  plt_biom%LeafPetoNonstElmConc_brch    , &
    PetioleProteinCNode_brch     =>  plt_biom%PetioleProteinCNode_brch    , &
    GrainStrutElms_brch          =>  plt_biom%GrainStrutElms_brch    , &
    ZEROP                        =>  plt_biom%ZEROP     , &
    ZEROL                        =>  plt_biom%ZEROL     , &
    iPlantPhotosynthesisType     =>  plt_photo%iPlantPhotosynthesisType   , &
    iPlantBranchState_brch       =>  plt_pheno%iPlantBranchState_brch    , &
    KLowestGroLeafNode_brch      =>  plt_pheno%KLowestGroLeafNode_brch   , &    
    KHiestGroLeafNode_brch       =>  plt_pheno%KHiestGroLeafNode_brch   , &
    iPlantRootProfile_pft        =>  plt_pheno%iPlantRootProfile_pft   , &
    fTgrowCanP                   =>  plt_pheno%fTgrowCanP     , &
    HoursDoingRemob_brch         =>  plt_pheno%HoursDoingRemob_brch    , &
    iPlantCalendar_brch          =>  plt_pheno%iPlantCalendar_brch   , &
    iPlantTurnoverPattern_pft    =>  plt_pheno%iPlantTurnoverPattern_pft   , &
    iPlantPhenolPattern_pft      =>  plt_pheno%iPlantPhenolPattern_pft   , &
    TCelciusChill4Seed           =>  plt_pheno%TCelciusChill4Seed     , &
    SineSolarIncliAngleNextHour  =>  plt_rad%SineSolarIncliAngleNextHour      , &
    PlantPopulation_pft          =>  plt_site%PlantPopulation_pft        , &
    ZERO                         =>  plt_site%ZERO      , &
    PSICanopy_pft                =>  plt_ew%PSICanopy_pft      , &
    MainBranchNum_pft            =>  plt_morph%MainBranchNum_pft      , &
    HypoctoHeight_pft            =>  plt_morph%HypoctoHeight_pft   , &
    SeedDepth_pft                =>  plt_morph%SeedDepth_pft   , &
    InternodeHeightLive_brch     =>  plt_morph%InternodeHeightLive_brch  , &
    NumCogrowNode                =>  plt_morph%NumCogrowNode   , &
    PARTS_brch                   =>  plt_morph%PARTS_brch         &
  )

  II=I;JJ=J
  LeafPetolBiomassC_brch(NB,NZ)=AZMAX1(LeafStrutElms_brch(ielmc,NB,NZ)+PetoleStrutElms_brch(ielmc,NB,NZ))
!  write(101,*)'grow branch',NB,NZ,LeafPetolBiomassC_brch(NB,NZ)
  IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN

    call CalcPartitionCoeff(I,J,NB,NZ,PART,PTRT,IFLGY,BegRemoblize)

    PARTS_brch(:,NB,NZ)=PART
!
!   SHOOT COEFFICIENTS FOR GROWTH RESPIRATION AND N,P CONTENTS
!   FROM GROWTH YIELDS ENTERED IN 'READQ', AND FROM PARTITIONING
!   COEFFICIENTS ABOVE
!
!   DM*B=C production vs nonstructural C consumption
!   CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood content
!   *LF=leaf,*SHE=petiole,*STK=stalk,*RSV=stalk reserve,*HSK=husk
!   *EAR=ear,*GR=grain from PFT file,*SH=shoot
!   DMSHT=branch C production vs nonstructural C consumption
!   DMSHD=branch C respiration vs nonstructural C consumption
!   CN*M,CP*M=min N,P production vs nonstructural C consumption
!   CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
!   CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   ZPLFD=1.0_r8-ZPLFM
!
    IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0)THEN
      DMLFB=LeafBiomGrowthYield(NZ)
      DMSHB=PetioleBiomGrowthYield(NZ)
      CNLFB=CNLFW
      CNSHB=CNSHW
      CPLFB=CPLFW
      CPSHB=CPSHW
    ELSE
      DMLFB=RootBiomGrowthYield(NZ)
      DMSHB=RootBiomGrowthYield(NZ)
      CNLFB=CNRTW
      CNSHB=CNRTW
      CPLFB=CPRTW
      CPSHB=CPRTW
    ENDIF
    !part 1(leaf), (2) petiole, (3) stalk, (4) reserve, (5) husk, (6) ear, (7) grain
    DMSHT=PART(ibrch_leaf)*DMLFB+PART(ibrch_petiole)*DMSHB+PART(ibrch_stalk)*StalkBiomGrowthYield(NZ) &
      +PART(ibrch_reserve)*ReserveBiomGrowthYield(NZ)+PART(ibrch_husk)*HuskBiomGrowthYield(NZ) &
      +PART(ibrch_ear)*EarBiomGrowthYield(NZ)+PART(ibrch_grain)*GrainBiomGrowthYield(NZ)

    DMSHD=1.0_r8-DMSHT
    CNLFM=PART(ibrch_leaf)*DMLFB*ZPLFM*CNLFB
    CPLFM=PART(ibrch_leaf)*DMLFB*ZPLFM*CPLFB
    CNLFX=PART(ibrch_leaf)*DMLFB*ZPLFD*CNLFB
    CPLFX=PART(ibrch_leaf)*DMLFB*ZPLFD*CPLFB
    CNSHX=PART(ibrch_petiole)*DMSHB*CNSHB &
      +PART(ibrch_stalk)*StalkBiomGrowthYield(NZ)*rNCStalk_pft(NZ) &
      +PART(ibrch_reserve)*ReserveBiomGrowthYield(NZ)*rNCReserve_pft(NZ) &
      +PART(ibrch_husk)*HuskBiomGrowthYield(NZ)*rNCHusk_pft(NZ) &
      +PART(ibrch_ear)*EarBiomGrowthYield(NZ)*rNCEar_pft(NZ) &
      +PART(ibrch_grain)*GrainBiomGrowthYield(NZ)*rNCReserve_pft(NZ)
    CPSHX=PART(ibrch_petiole)*DMSHB*CPSHB &
      +PART(ibrch_stalk)*StalkBiomGrowthYield(NZ)*rPCStalk_pft(NZ) &
      +PART(ibrch_reserve)*ReserveBiomGrowthYield(NZ)*rPCReserve_pft(NZ) &
      +PART(ibrch_husk)*HuskBiomGrowthYield(NZ)*rPCHusk_pft(NZ) &
      +PART(ibrch_ear)*EarBiomGrowthYield(NZ)*rPCEar_pft(NZ) &
      +PART(ibrch_grain)*GrainBiomGrowthYield(NZ)*rPCReserve_pft(NZ)
!
!   TOTAL SHOOT STRUCTURAL N MASS FOR MAINTENANCE RESPIRATION
!
!   ShootStructN=shoot structural N mass
!   WTLFBN,WTSHBN,WTHSBN,WTEARN,WTFRBN=leaf,petiole,husk,ear,grain N mass
!   rNCStalk_pft,StalkBiomassC_brch=stalk N:C,sapwood mass
!   iPlantCalendar_brch(10=date of physiological maturity
!
    ShootStructN=AZMAX1(LeafStrutElms_brch(ielmn,NB,NZ)+PetoleStrutElms_brch(ielmn,NB,NZ) &
      +rNCStalk_pft(NZ)*StalkBiomassC_brch(NB,NZ))
    
    IF(iPlantCalendar_brch(ipltcal_EndSeedFill,NB,NZ).EQ.0)THEN
      ShootStructN=ShootStructN+AZMAX1(HuskStrutElms_brch(ielmn,NB,NZ) &
        +EarStrutElms_brch(ielmn,NB,NZ)+GrainStrutElms_brch(ielmn,NB,NZ))
    ENDIF
!
!   GROSS PRIMARY PRODUCTIVITY
!
    call UpdatePhotosynthates(I,J,NB,NZ,TFN6_vr,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX &
      ,CNLFX,CPLFX,ShootStructN,TFN5,WFNG,Stomata_Activity,WFNSG,CH2O3,CH2O4,CNPG &
      ,RCO2C,RMNCS,RMxess_brch,NonStructalC4Growth_brch,CNRDM,Rauto4Nassim_brch)
!
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
!
!   C,N,P GROWTH OF LEAF, SHEATH OR PETIOLE, STALK,
!   STALK RESERVES, REPRODUCTIVE ORGANS, GRAIN
!
!   GRO*,GRO*N,GRO*P=organ C,N,P growth rate
!   DM*=C production vs nonstructural C consumption
!   organ key:LF=leaf,SHE=petiole,STK=stalk,RSV=reserve
!   HSK=husk,EAR=ear,GR=grain,SHT=shoot
!   PART=organ partitioning fraction
!   NonStructalC4Growth_brch=total non-structural C used in growth and growth respiration
!   CN*,CP*=N:C,P:C ratios in plant organs
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   ZPLFD=1.0_r8-ZPLFM
!   CNPG=N,P constraint on growth respiration
!   WT*,WT*N,WT*P=organ C,N,P mass
!
    GrowthLeaf(ielmc)   =NonStructalC4Growth_brch*PART(ibrch_leaf)*DMLFB       !leaf
    GrowthPetiole(ielmc)=NonStructalC4Growth_brch*PART(ibrch_petiole)*DMSHB       !petiole
    GrowthStalk(ielmc)  =NonStructalC4Growth_brch*PART(ibrch_stalk)*StalkBiomGrowthYield(NZ)  !stalk
    GrowthReserve(ielmc)=NonStructalC4Growth_brch*PART(ibrch_reserve)*ReserveBiomGrowthYield(NZ)  !reserve
    GrowthHusk(ielmc)   =NonStructalC4Growth_brch*PART(ibrch_husk)*HuskBiomGrowthYield(NZ)  !husk
    GrowthEar(ielmc)    =NonStructalC4Growth_brch*PART(ibrch_ear)*EarBiomGrowthYield(NZ)  !ear
    GrowthGrain(ielmc)  =NonStructalC4Growth_brch*PART(ibrch_grain)*GrainBiomGrowthYield(NZ)    !grain

    GrowthLeaf(ielmn)   =GrowthLeaf(ielmc)*CNLFB*(ZPLFM+ZPLFD*CNPG)
    GrowthPetiole(ielmn)=GrowthPetiole(ielmc)*CNSHB
    GrowthStalk(ielmn)  =GrowthStalk(ielmc)*rNCStalk_pft(NZ)
    GrowthReserve(ielmn)=GrowthReserve(ielmc)*rNCReserve_pft(NZ)
    GrowthHusk(ielmn)   =GrowthHusk(ielmc)*rNCHusk_pft(NZ)
    GrowthEar(ielmn)    =GrowthEar(ielmc)*rNCEar_pft(NZ)
    GrowthGrain(ielmn)  =GrowthGrain(ielmc)*rNCReserve_pft(NZ)

    GrowthLeaf(ielmp)=GrowthLeaf(ielmc)*CPLFB*(ZPLFM+ZPLFD*CNPG)
    GrowthPetiole(ielmp)=GrowthPetiole(ielmc)*CPSHB
    GrowthStalk(ielmp)=GrowthStalk(ielmc)*rPCStalk_pft(NZ)
    GrowthReserve(ielmp)=GrowthReserve(ielmc)*rPCReserve_pft(NZ)
    GrowthHusk(ielmp)=GrowthHusk(ielmc)*rPCHusk_pft(NZ)
    GrowthEar(ielmp)=GrowthEar(ielmc)*rPCEar_pft(NZ)
    GrowthGrain(ielmp)=GrowthGrain(ielmc)*rPCReserve_pft(NZ)

    if(NZ==1)THEN
      write(185,'(I3,X,I4,9(X,F13.6))')NB,etimer%get_curr_year(),I+J/24.&
       ,NonStructalC4Growth_brch,PART(1:7)
!      GrowthStalk(NE),GrowthReserve(NE),GrowthHusk(NE),GrowthEar(NE),GrowthGrain(NE),NE=1,NumPlantChemElms)
!      WRITE(187,'(I3,8(X,F13.6))')NB,I+J/24.,(PART(NE),NE=1,7)
    else
      write(186,'(I3,X,I4,9(X,F13.6))')NB,etimer%get_curr_year(),I+J/24.&
       ,NonStructalC4Growth_brch,PART(1:7)
!      GrowthStalk(NE),GrowthReserve(NE),GrowthHusk(NE),GrowthEar(NE),GrowthGrain(NE),NE=1,NumPlantChemElms)
!      WRITE(188,'(I3,8(X,F13.6))')NB,I+J/24.,(PART(NE),NE=1,7)
    endif
    
    DO NE=1,NumPlantChemElms
      LeafStrutElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ)+GrowthLeaf(NE)
      PetoleStrutElms_brch(NE,NB,NZ)=PetoleStrutElms_brch(NE,NB,NZ)+GrowthPetiole(NE)
      StalkStrutElms_brch(NE,NB,NZ)=StalkStrutElms_brch(NE,NB,NZ)+GrowthStalk(NE)
      StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+GrowthReserve(NE)
      HuskStrutElms_brch(NE,NB,NZ)=HuskStrutElms_brch(NE,NB,NZ)+GrowthHusk(NE)
      EarStrutElms_brch(NE,NB,NZ)=EarStrutElms_brch(NE,NB,NZ)+GrowthEar(NE)
    ENDDO
    if(NZ==1)THEN
      WRITE(203,'(I3,7(X,F16.6))')NB,I+J/24.,(StalkRsrvElms_brch(NE,NB,NZ),GrowthReserve(NE),NE=1,NumPlantChemElms)
    ELSE
      WRITE(204,'(I3,7(X,F16.6))')NB,I+J/24.,(StalkRsrvElms_brch(NE,NB,NZ),GrowthReserve(NE),NE=1,NumPlantChemElms)
    ENDIF

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

    ETOL=1.0_r8+CCE
!
!   DISTRIBUTE LEAF GROWTH AMONG CURRENTLY GROWING NODES
!
    IF(NB.EQ.MainBranchNum_pft(NZ).AND.HypoctoHeight_pft(NZ).LE.SeedDepth_pft(NZ))THEN
      NNOD1=0
    ELSE
      NNOD1=1
    ENDIF
!    if(NZ==1)then
!      write(163,*)'ht',I,NB,MainBranchNum_pft(NZ),HypoctoHeight_pft(NZ),SeedDepth_pft(NZ)
!    else
!      write(164,*)'ht',I,NB,MainBranchNum_pft(NZ),HypoctoHeight_pft(NZ),SeedDepth_pft(NZ)
!    endif

    CALL GrowLeavesOnBranch(NZ,NB,NNOD1,GrowthLeaf,ETOL,WFNS,ALLOCL)
!
!     DISTRIBUTE SHEATH OR PETIOLE GROWTH AMONG CURRENTLY GROWING NODES
!
    CALL GrowPetioleOnBranch(NZ,NB,NNOD1,GrowthPetiole,ETOL,WFNS,ALLOCL)

!
!   DISTRIBUTE STALK GROWTH AMONG CURRENTLY GROWING NODES
!
    call GrowStalkOnBranch(NZ,NB,GrowthStalk,ETOL)

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
        ,safe_adb(LeafPetoNonstElmConc_brch(ielmn,NB,NZ),LeafPetoNonstElmConc_brch(ielmn,NB,NZ) &
        +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI) &
        ,safe_adb(LeafPetoNonstElmConc_brch(ielmp,NB,NZ),LeafPetoNonstElmConc_brch(ielmp,NB,NZ) &
        +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI)))
      CNC=AZMAX1(AMIN1(1.0_r8 &
        ,safe_adb(LeafPetoNonstElmConc_brch(ielmc,NB,NZ),LeafPetoNonstElmConc_brch(ielmc,NB,NZ) &
        +LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/CNKI)))
      CPC=AZMAX1(AMIN1(1.0_r8 &
        ,safe_adb(LeafPetoNonstElmConc_brch(ielmc,NB,NZ),LeafPetoNonstElmConc_brch(ielmc,NB,NZ) &
        +LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/CPKI)))
    ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
    ENDIF
    RCCC=RCCZ(iPlantRootProfile_pft(NZ))+CCC*RCCY(iPlantRootProfile_pft(NZ))
    RCCN=CNC*RCCX(iPlantRootProfile_pft(NZ))
    RCCP=CPC*RCCQ(iPlantRootProfile_pft(NZ))
!
!       WITHDRAW REMOBILIZABLE C,N,P FROM LOWEST NODE AFTER
!       MAXIMUM NODE NUMBER OF 25 IS REACHED
!
!       doSenescence_brch=PFT branch senescence flag
!       KHiestGroLeafNode_brch=integer of most recent leaf number
!       fTgrowCanP=temperature function for canopy growth
!       XRLA=rate of leaf appearance at 25 oC (h-1)
!       FSNC=fraction of lowest leaf to be remobilized
!
    call SenescenceBranch(NZ,NB,RCCC,RCCN,RCCP)
!
    call RemobilizeBranch(NZ,NB,BegRemoblize,IFLGY,RCCC,RCCN,RCCP,RMxess_brch)
!
!
!   DEATH IF MAIN STALK OF TREE DIES
!
!   iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+petiole,2=none,3=between 1,2
!   iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!   iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!   KMinGroingLeafNodeNum,KHiestGroLeafNode_brch=integer of lowest,highest leaf number currently growing
!   WGLF,WGLFN,WGLFP,LeafProteinCNode_brch=node leaf C,N,P,protein mass
!   CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!   ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
!   CNLFB,CPLFB=N:C,P:C ratios in leaf
!
    IF(iPlantTurnoverPattern_pft(NZ).NE.0.AND.is_plant_treelike(iPlantRootProfile_pft(NZ)) &
      .AND.iPlantBranchState_brch(MainBranchNum_pft(NZ),NZ).EQ.iDead)then
      iPlantBranchState_brch(NB,NZ)=iDead
    endif

!
!     REMOBILIZE EXCESS LEAF STRUCTURAL N,P
!
    KMinGroingLeafNodeNum=MAX(0,KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+1)
    D495: DO KK=KMinGroingLeafNodeNum,KHiestGroLeafNode_brch(NB,NZ)
      K=pMOD(KK,MaxNodesPerBranch1)
      IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.0.0_r8)THEN
        CPOOLT=LeafElmntNode_brch(ielmc,K,NB,NZ)+CanopyNonstElms_brch(ielmc,NB,NZ)
        IF(CPOOLT.GT.ZEROP(NZ))THEN
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
            LeafElmntNode_brch(NE,K,NB,NZ)=LeafElmntNode_brch(NE,K,NB,NZ)-XFRE(NE)
            LeafStrutElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ)-XFRE(NE)
            CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE)
            if(CanopyNonstElms_brch(NE,NB,NZ)<0._r8)then
              write(*,*)'435negative CanopyNonstElms_brch',NE,NB,NZ,CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE),XFRE(NE)
              stop
            endif
          ENDDO
          LeafProteinCNode_brch(K,NB,NZ)=AZMAX1(LeafProteinCNode_brch(K,NB,NZ)-&
            AMAX1(XFRE(ielmn)*rCNNonstructRemob_pft(NZ),XFRE(ielmp)*rCPNonstructRemob_pft(NZ)))
        ENDIF
      ENDIF
    ENDDO D495
!   KK inherits value from loop 495, is it right?

    call AllocateLeafToCanopyLayers(NB,NZ,CanopyHeight_copy)

!
    !     ALLOCATE LEAF AREA TO INCLINATION CLASSES ACCORDING TO
    !     DISTRIBUTION ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
    !
    !     SineSolarIncliAngle=sine of solar angle
    !     LeafAreaZsec_brch=leaf node surface area in canopy layer
    !     LeafAreaNode_brch,CanopyLeafAreaByLayer_pft=leaf node surface area in canopy layer
    !     ZC,DPTHS=canopy,snowpack height
    !     CLASS=leaf inclination class
    !
    IF(SineSolarIncliAngleNextHour.GT.0.0_r8)THEN
      call LeafClassAllocation(NB,NZ)
    ENDIF

    call GrainFilling(I,NB,NZ,GrowthGrain,GrowthStalk(ielmc))
!
    call PhenologyReset(I,NB,NZ)
!   
    print*,'branchbf',NZ,StalkRsrvElms_brch(1,1,NZ)
    call BranchElmntTransfer(I,J,NB,NZ,BegRemoblize,WFNG,WFNSG)

!   CANOPY N2 FIXATION (CYANOBACTERIA)
!
    print*,'nodulebf',NZ,StalkRsrvElms_brch(1,1,NZ)
    call CanopyNoduleBiochemistry(I,J,NZ,NB,TFN5,WFNG,CanopyN2Fix_pft)
  ENDIF
  end associate
  end subroutine GrowOneBranch

!------------------------------------------------------------------------------------------

  subroutine CalcPartitionCoeff(I,J,NB,NZ,part,PTRT,BegRemoblize,IFLGY)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  integer, intent(out) :: BegRemoblize,IFLGY
  real(r8), intent(out):: PART(pltpar%NumOfPlantMorphUnits),PTRT
  REAL(R8) :: ARLFI

  integer :: N,LP
  real(r8) :: FPARTL
  real(r8) :: PARTS
  real(r8) :: PARTX
  real(r8) :: TOTAL
  real(r8) :: PSILY(0:3)
  real(r8), parameter :: FPART1=1.00_r8
  real(r8), parameter :: FPART2=0.40_r8

  associate(                              &
    TCelciusCanopy_pft                 =>  plt_ew%TCelciusCanopy_pft         , &
    PSICanopy_pft                      =>  plt_ew%PSICanopy_pft       , &
    StalkBiomassC_brch                 =>  plt_biom%StalkBiomassC_brch     , &
    StalkRsrvElms_brch                =>  plt_biom%StalkRsrvElms_brch    , &
    ZERO                               =>  plt_site%ZERO       , &
    NU                                 =>  plt_site%NU         , &
    AREA3                              =>  plt_site%AREA3      , &
    FracHour4LeafoffRemob              =>  plt_allom%FracHour4LeafoffRemob      , &
    HourReq4LeafOff_brch               =>  plt_pheno%HourReq4LeafOff_brch     , &
    iPlantCalendar_brch                =>  plt_pheno%iPlantCalendar_brch    , &
    iPlantRootProfile_pft              =>  plt_pheno%iPlantRootProfile_pft    , &
    TotReproNodeNumNormByMatrgrp_brch  =>  plt_pheno%TotReproNodeNumNormByMatrgrp_brch    , &
    iPlantDevelopPattern_pft           =>  plt_pheno%iPlantDevelopPattern_pft    , &
    HoursDoingRemob_brch               =>  plt_pheno%HoursDoingRemob_brch     , &
    iPlantTurnoverPattern_pft          =>  plt_pheno%iPlantTurnoverPattern_pft    , &
    Hours4LeafOff_brch                 =>  plt_pheno%Hours4LeafOff_brch      , &
    TotalNodeNumNormByMatgrp_brch      =>  plt_pheno%TotalNodeNumNormByMatgrp_brch    , &
    iPlantPhenolPattern_pft         =>  plt_pheno%iPlantPhenolPattern_pft    , &
    iPlantPhenolType_pft               =>  plt_pheno%iPlantPhenolType_pft    , &
    TCelciusChill4Seed                 =>  plt_pheno%TCelciusChill4Seed      , &
    NH3Dep2Can_brch                       =>  plt_rbgc%NH3Dep2Can_brch      , &
    CO2NetFix_pft                      =>  plt_bgcr%CO2NetFix_pft       , &
    NodeLenPergC                       =>  plt_morph%NodeLenPergC      , &
    CanopyLeafArea_pft                 =>  plt_morph%CanopyLeafArea_pft     , &
    MainBranchNum_pft                  =>  plt_morph%MainBranchNum_pft         &
  )

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
!
!  if(NZ==1)THEN
!  write(195,'(I3,F13.6,10(X,I6)))')NB,I+J/24.,(iPlantCalendar_brch(LP,NB,NZ),LP=1,10)
!  ELSE
!  write(196,'(I3,F13.6,10(X,I6)))')NB,I+J/24.,(iPlantCalendar_brch(LP,NB,NZ),LP=1,10)
!  ENDIF
  
  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
    PART(ibrch_leaf)=0.725_r8
    PART(ibrch_petiole)=0.275_r8
!
!     IF BEFORE ANTHESIS
!
!     iPlantCalendar_brch(ipltcal_Anthesis,=start of anthesis and setting final seed number
!     TotalNodeNumNormByMatgrp_brch=total change in vegv node number normalized for maturity group
!
  ELSEIF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
    PART(ibrch_leaf)=AMAX1(PART1X,0.725_r8-FPART1*TotalNodeNumNormByMatgrp_brch(NB,NZ))
    PART(ibrch_petiole)=AMAX1(PART2X,0.275_r8-FPART2*TotalNodeNumNormByMatgrp_brch(NB,NZ))
    PARTS=1.0_r8-PART(ibrch_leaf)-PART(ibrch_petiole)
    PART(ibrch_stalk)=0.60_r8*PARTS
    PART(ibrch_reserve)=0.30_r8*PARTS
    PARTX=PARTS-PART(ibrch_stalk)-PART(ibrch_reserve)
    PART(ibrch_husk)=0.5_r8*PARTX
    PART(ibrch_ear)=0.5_r8*PARTX
!
!     IF BEFORE GRAIN FILLING, DETERMINATE OR INDETERMINATE
!
!     iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
!     iPlantDevelopPattern_pft=growth habit:0=determinate,1=indetermimate from PFT file
!     TotReproNodeNumNormByMatrgrp_brch=total change in reprv node number normalized for maturity group
!
  ELSEIF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
    IF(iPlantDevelopPattern_pft(NZ).EQ.ideterminate)THEN
      PART(ibrch_leaf)=0._r8
      PART(ibrch_petiole)=0._r8
    ELSE
      PART(ibrch_leaf)=AMAX1(PART1X,(0.725_r8-FPART1)*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
      PART(ibrch_petiole)=AMAX1(PART2X,(0.275_r8-FPART2)*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
    ENDIF
    PARTS=1.0_r8-PART(ibrch_leaf)-PART(ibrch_petiole)
    PART(ibrch_stalk)=AZMAX1(0.60_r8*PARTS*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
    PART(ibrch_reserve)=AZMAX1(0.30_r8*PARTS*(1.0_r8-TotReproNodeNumNormByMatrgrp_brch(NB,NZ)))
    PARTX=PARTS-PART(ibrch_stalk)-PART(ibrch_reserve)
    PART(ibrch_husk)=0.5_r8*PARTX
    PART(ibrch_ear)=0.5_r8*PARTX
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
      PART(ibrch_leaf)=PART1X
      PART(ibrch_petiole)=PART2X
      PARTS=1.0_r8-PART(ibrch_leaf)-PART(ibrch_petiole)
      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
        PART(ibrch_stalk)=0.125*PARTS
        PART(ibrch_husk)=0.125*PARTS
        PART(ibrch_ear)=0.125*PARTS
        PART(ibrch_grain)=0.625*PARTS
      ELSE
        PART(ibrch_stalk)=0.75*PARTS
        PART(ibrch_grain)=0.25*PARTS
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
      PART(ibrch_reserve)=0._r8
      PART(ibrch_stalk)=0._r8
      PART(ibrch_grain)=0._r8
    ELSE
      PART(ibrch_reserve)=PART(ibrch_reserve)+PART(ibrch_stalk)
      PART(ibrch_stalk)=0._r8
      PART(ibrch_grain)=0._r8
    ENDIF
  ENDIF
!
!     REDIRECT FROM STALK TO STALK RESERVES IF RESERVES BECOME LOW
!
!     WTRSVB,StalkBiomassC_brch=stalk reserve,sapwood mass
!     XFRX=maximum storage C content for remobiln from stalk,root reserves
!
  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)THEN
    IF(StalkRsrvElms_brch(ielmc,NB,NZ).LT.XFRX*StalkBiomassC_brch(NB,NZ))THEN
      D1020: DO N=1,NumOfPlantMorphUnits
        IF(N.NE.4)THEN
          PART(ibrch_reserve)=PART(ibrch_reserve)+0.10*PART(N)
          PART(N)=PART(N)-0.10*PART(N)
        ENDIF
      ENDDO D1020
!
!     REDIRECT FROM STALK RESERVES TO STALK IF RESERVES BECOME TOO LARGE
!
    ELSEIF(StalkRsrvElms_brch(ielmc,NB,NZ).GT.1.0*StalkBiomassC_brch(NB,NZ))THEN
      PART(ibrch_stalk)=PART(ibrch_stalk)+PART(ibrch_reserve)+PART(ibrch_grain)
      PART(ibrch_reserve)=0._r8
      PART(ibrch_grain)=0._r8
    ENDIF
  ENDIF
!
!     REDIRECT FROM LEAVES TO STALK IF LAI BECOMES TOO LARGE
!
!     CanopyLeafArea_pft=PFT leaf area
!
  ARLFI=CanopyLeafArea_pft(NZ)/AREA3(NU)
  IF(ARLFI.GT.5.0_r8)THEN
    FPARTL=AZMAX1((10.0_r8-ARLFI)/5.0_r8)
    PART(ibrch_stalk)=PART(ibrch_stalk)+(1.0_r8-FPARTL)*(PART(ibrch_leaf)+PART(ibrch_petiole))
    PART(ibrch_leaf)=FPARTL*PART(ibrch_leaf)
    PART(ibrch_petiole)=FPARTL*PART(ibrch_petiole)
  ENDIF
  IF(NB.EQ.MainBranchNum_pft(NZ))THEN
    PTRT=PART(ibrch_leaf)+PART(ibrch_petiole)
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
!     IFLGY,BegRemoblize=remobilization flags
!     FLGZ=control rate of remobilization
!
  IF((iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. &
     Hours4LeafOff_brch(NB,NZ).GE.FracHour4LeafoffRemob(iPlantPhenolType_pft(NZ))*HourReq4LeafOff_brch(NB,NZ)) &
    .OR.(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).NE.0))THEN
    !set remobilization true
    BegRemoblize=1
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.OR.iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen)THEN
      !annual plant or evergreen perennial
      IFLGY=itrue
      HoursDoingRemob_brch(NB,NZ)=HoursDoingRemob_brch(NB,NZ)+1.0_r8
    ELSEIF((iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid.OR.&
      iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid).AND.&
      TCelciusCanopy_pft(NZ).LT.TCelciusChill4Seed(NZ))THEN
      IFLGY=itrue
      HoursDoingRemob_brch(NB,NZ)=HoursDoingRemob_brch(NB,NZ)+1.0_r8
    ELSEIF(iPlantPhenolType_pft(NZ).GE.2 .AND. PSICanopy_pft(NZ).LT.PSILY(iPlantRootProfile_pft(NZ)))THEN
      IFLGY=itrue
      HoursDoingRemob_brch(NB,NZ)=HoursDoingRemob_brch(NB,NZ)+1.0_r8
    ENDIF
    IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.iphenotyp_evgreen)THEN
      PART(ibrch_stalk)=PART(ibrch_stalk)+0.5_r8*(PART(ibrch_leaf)+PART(ibrch_petiole))
      PART(ibrch_reserve)=PART(ibrch_reserve)+0.5_r8*(PART(ibrch_leaf)+PART(ibrch_petiole))
      PART(ibrch_leaf)=0._r8
      PART(ibrch_petiole)=0._r8
    ENDIF
  ELSE
    BegRemoblize=0
    IFLGY=0
    HoursDoingRemob_brch(NB,NZ)=0._r8
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
  end associate
  end subroutine CalcPartitionCoeff

!------------------------------------------------------------------------------------------

  subroutine UpdatePhotosynthates(I,J,NB,NZ,TFN6_vr,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
    ShootStructN,TFN5,WFNG,Stomata_Activity,WFNSG,CH2O3,CH2O4,CNPG,RCO2C,RMNCS,&
    RMxess_brch,NonStructalC4Growth_brch,CNRDM,Rauto4Nassim_brch)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: DMSHD
  real(r8), intent(in) :: CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,ShootStructN,TFN5,WFNG
  real(r8), intent(in) :: Stomata_Activity,WFNSG
  real(r8), intent(out) :: RCO2C
  real(r8), intent(out) :: RMNCS
  real(r8), intent(out) :: CH2O3(pltpar%MaxNodesPerBranch1)
  real(r8), intent(out) :: CH2O4(pltpar%MaxNodesPerBranch1)
  REAL(R8), INTENT(OUT) :: CNPG
  real(r8), intent(out) :: RMxess_brch
  real(r8), intent(out) :: NonStructalC4Growth_brch
  real(r8), intent(out) :: CNRDM,Rauto4Nassim_brch
  real(r8) :: CO2F,CanopyNonstElm4Gros(NumPlantChemElms),CH2O,CH2OClm,CH2OLlm
  integer :: K
  associate(                        &
    NH3Dep2Can_brch       =>  plt_rbgc%NH3Dep2Can_brch   , &
    CanopyNonstElms_brch  =>  plt_biom%CanopyNonstElms_brch  , &
    iPlantCalendar_brch   =>  plt_pheno%iPlantCalendar_brch   &
  )
! begin_execution
! FDBK=N,P feedback inhibition on C3 CO2 fixation
! SineSolarIncliAngle=sine of solar angle
! RADP=total PAR absorbed by canopy
! CanopyGasCO2_pft=canopy air CO2 concentration
!
!  write(101,*)'emerge cal',NB,NZ,iPlantCalendar_brch(ipltcal_Emerge,NB,NZ)
  IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).NE.0)THEN
!    write(101,*)'computegpp branch',NB,NZ
    call ComputeGPP(NB,NZ,WFNG,Stomata_Activity,CH2O3,CH2O4,CH2O,CO2F,CH2OClm,CH2OLlm)
!    if(NZ==1)then
!    write(101,'(I4,5(X,F13.6))')NB,I+J/24.,CO2F,CH2O,CH2OClm,CH2OLlm
!    else
!    write(102,'(I4,5(X,F13.6))')NB,I+J/24.,CO2F,CH2O,CH2OClm,CH2OLlm
!    endif
!    CH2O=max(CH2OClm,CH2OLlm);CO2F=CH2O
!
!   SHOOT AUTOTROPHIC RESPIRATION AFTER EMERGENCE
!
    call ComputRAutoAfEmergence(NB,NZ,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,&
      CO2F,CH2O,TFN5,WFNG,WFNSG,ShootStructN,CanopyNonstElm4Gros,CNPG,RCO2C,RMNCS,RMxess_brch,&
      NonStructalC4Growth_brch,Rauto4Nassim_brch)

!   SHOOT AUTOTROPHIC RESPIRATION BEFORE EMERGENCE
!
  ELSE
    CH2O=0._r8
    call ComputRAutoB4Emergence(I,NB,NZ,TFN6_vr,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,ShootStructN,&
      WFNG,WFNSG,CanopyNonstElm4Gros,CNPG,RCO2C,RMNCS,RMxess_brch,NonStructalC4Growth_brch,CNRDM,Rauto4Nassim_brch)
!    write(101,*)'ComputRAutoB4Emergence',NonStructalC4Growth_brch  
  ENDIF

!   REMOVE C,N,P USED IN MAINTENANCE + GROWTH REPIRATION AND GROWTH
!   FROM NON-STRUCTURAL POOLS
!
!   CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!   CH2O=total CH2O production
!   RMNCS=maintenance respiration
!   RCO2C=respiration from non-structural C
!   NonStructalC4Growth_brch=total non-structural C used in growth and respiration
!   Rauto4Nassim_brch=respiration for N assimilation
!   CanopyNonstElm4Gros(ielmn),CanopyNonstElm4Gros(ielmp)=nonstructural N,P used in growth
!   NH3Dep2Can_brch=NH3 flux between atmosphere and branch from uptake.f
!   XFRE(ielmc),XFRE(ielmn),XFRE(ielmp)=branch-root layer C,N,P transfer
!
  CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)+CH2O-&
    AMIN1(RMNCS,RCO2C)-NonStructalC4Growth_brch-Rauto4Nassim_brch
  CanopyNonstElms_brch(ielmn,NB,NZ)=CanopyNonstElms_brch(ielmn,NB,NZ)-CanopyNonstElm4Gros(ielmn)+NH3Dep2Can_brch(NB,NZ)
  CanopyNonstElms_brch(ielmp,NB,NZ)=CanopyNonstElms_brch(ielmp,NB,NZ)-CanopyNonstElm4Gros(ielmp)
  if(CanopyNonstElms_brch(ielmp,NB,NZ)<0._r8)then
    write(*,*)'806CanopyNonstElms_brch(ielmp,NB,NZ)',NB,CanopyNonstElms_brch(ielmp,NB,NZ)+CanopyNonstElm4Gros(ielmp),CanopyNonstElm4Gros(ielmp)
    stop
  endif
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
  associate(                              &
    CO2NetFix_pft                =>  plt_bgcr%CO2NetFix_pft       , &
    GrossResp_pft                =>  plt_bgcr%GrossResp_pft      , &
    CanopyPlusNoduRespC_pft      =>  plt_bgcr%CanopyPlusNoduRespC_pft      , &
    Eco_AutoR_col                =>  plt_bgcr%Eco_AutoR_col       , &
    ECO_ER_col                   =>  plt_bgcr%ECO_ER_col       , &
    LeafElmntNode_brch           =>  plt_biom%LeafElmntNode_brch      , &
    ZEROP                        =>  plt_biom%ZEROP      , &
    CMassHCO3BundleSheath_node   =>  plt_photo%CMassHCO3BundleSheath_node      , &
    CMassCO2BundleSheath_node    =>  plt_photo%CMassCO2BundleSheath_node      , &
    CPOOL3_node                  =>  plt_photo%CPOOL3_node    , &
    CPOOL4_node                  =>  plt_photo%CPOOL4_node    , &
    aquCO2Intraleaf_pft          =>  plt_photo%aquCO2Intraleaf_pft        &
  )
  D170: DO K=1,MaxNodesPerBranch1
    IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZEROP(NZ))THEN
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
      CPL3K=2.5E-02_r8*CPOOL3_node(K,NB,NZ)/(1.0+CCBS/CO2KI)
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
!     GrossResp_pft,CanopyPlusNoduRespC_pft=total,above-ground PFT respiration
!     CO2NetFix_pft=PFT net CO2 fixation
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_col=total autotrophic respiration
!     aquCO2IntraLeafLeakFromBndsheth=bundle sheath CO2 leakage
!
      GrossResp_pft(NZ)=GrossResp_pft(NZ)-aquCO2IntraLeafLeakFromBndsheth
      CanopyPlusNoduRespC_pft(NZ)=CanopyPlusNoduRespC_pft(NZ)-aquCO2IntraLeafLeakFromBndsheth
      CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-aquCO2IntraLeafLeakFromBndsheth
      ECO_ER_col=ECO_ER_col-aquCO2IntraLeafLeakFromBndsheth
      Eco_AutoR_col=Eco_AutoR_col-aquCO2IntraLeafLeakFromBndsheth
    ENDIF
  ENDDO D170
  end associate
  end subroutine C4PhotoProductTransfer
!------------------------------------------------------------------------------------------
  subroutine RemobilizeLeafLayers(NumRemobLeafNodes,NB,NZ,XKSNC,SNCX,RCCC,RCCN,RCCP,SNCF)
  implicit none
  INTEGER, intent(in)    :: NB,NZ,NumRemobLeafNodes
  real(r8), intent(in)   :: XKSNC,SNCX
  REAL(R8), INTENT(IN) :: RCCC,RCCN,RCCP
  real(r8), intent(inout):: SNCF
  integer :: N,M,K,KK,MXNOD,MNNOD,NE
  real(r8) :: FSNCL,FSNCS
  real(r8) :: FNCLF,FSNAL,FSNAS,FRCC
  real(r8) :: FSNCK
  real(r8) :: FSNCR
  real(r8) :: HTNODZ
  real(r8) :: ElmntRemobFallingLeaf(NumPlantChemElms)
  real(r8) :: RCES(NumPlantChemElms)
  real(r8) :: FStalkSenes
  real(r8) :: RCSC,RCSN,RCSP
  real(r8) :: RCEK(NumPlantChemElms)
  real(r8) :: RMxess_brch
  real(r8) :: SNCZ
  real(r8) :: SNCT
  real(r8) :: SNCLF
  real(r8) :: SNCSH
  integer :: KN  
! begin_execution
  associate(                             &
    ifoliar                         =>  pltpar%ifoliar      , &
    istalk                          =>  pltpar%istalk       , &
    iroot                           =>  pltpar%iroot        , &
    inonfoliar                      =>  pltpar%inonfoliar     , &
    icwood                          =>  pltpar%icwood       , &
    k_fine_litr                     =>  pltpar%k_fine_litr  , &
    k_woody_litr                    =>  pltpar%k_woody_litr , &
    StalkRsrvElms_brch              =>  plt_biom%StalkRsrvElms_brch    , &
    StalkStrutElms_brch             =>  plt_biom%StalkStrutElms_brch    , &
    InternodeStrutElms_brch         =>  plt_biom%InternodeStrutElms_brch     , &
    StalkBiomassC_brch              =>  plt_biom%StalkBiomassC_brch     , &
    PetioleElmntNode_brch           =>  plt_biom%PetioleElmntNode_brch     , &
    LeafElmntNode_brch              =>  plt_biom%LeafElmntNode_brch      , &
    LeafProteinCNode_brch           =>  plt_biom%LeafProteinCNode_brch       , &
    NonStrutElms_pft                =>  plt_biom%NonStrutElms_pft      , &
    LeafStrutElms_brch              =>  plt_biom%LeafStrutElms_brch    , &
    PetoleStrutElms_brch            =>  plt_biom%PetoleStrutElms_brch   , &
    SenecStalkStrutElms_brch        =>  plt_biom%SenecStalkStrutElms_brch    , &
    PetioleProteinCNode_brch        =>  plt_biom%PetioleProteinCNode_brch     , &
    ZEROP                           =>  plt_biom%ZEROP      , &
    ZEROL                           =>  plt_biom%ZEROL      , &
    CanopyNonstElms_brch            =>  plt_biom%CanopyNonstElms_brch     , &
    PlantPopulation_pft             =>  plt_site%PlantPopulation_pft        , &
    CPOOL3_node                     =>  plt_photo%CPOOL3_node    , &
    CPOOL4_node                     =>  plt_photo%CPOOL4_node    , &
    FWOODE                          =>  plt_allom%FWOODE    , &
    FWODBE                          =>  plt_allom%FWODBE    , &
    FWODLE                          =>  plt_allom%FWODLE    , &
    rCNNonstructRemob_pft           =>  plt_allom%rCNNonstructRemob_pft     , &
    rCPNonstructRemob_pft           =>  plt_allom%rCPNonstructRemob_pft      , &
    LitfalStrutElms_pvr             =>  plt_bgcr%LitfalStrutElms_pvr       , &
    CFOPE                           =>  plt_soilchem%CFOPE  , &
    icarbhyro                       =>  pltpar%icarbhyro    , &
    KLowestGroLeafNode_brch         =>  plt_pheno%KLowestGroLeafNode_brch, &
    KHiestGroLeafNode_brch          =>   plt_pheno%KHiestGroLeafNode_brch   , &
    iPlantBranchState_brch          =>   plt_pheno%iPlantBranchState_brch    , &
    iPlantPhenolPattern_pft         =>   plt_pheno%iPlantPhenolPattern_pft   , &
    LeafAreaLive_brch               =>   plt_morph%LeafAreaLive_brch    , &
    NumCogrowNode                   =>   plt_morph%NumCogrowNode    , &
    InternodeHeightLive_brch        =>   plt_morph%InternodeHeightLive_brch   , &
    PetioleLengthNode_brch          =>   plt_morph%PetioleLengthNode_brch    , &
    LeafAreaNode_brch               =>   plt_morph%LeafAreaNode_brch    , &
    InternodeHeightDying_brch       =>   plt_morph%InternodeHeightDying_brch     &
  )
!     REMOBILIZATION AND LitrFall WHEN GROWTH RESPIRATION < 0
!     STARTING FROM LOWEST LEAFED NODE AND PROCEEDING UPWARDS
!
!     SNCX,SNCT=branch,node senescence respiration
!     KSNC=number of nodes undergoing remobilization
!
  KN=MAX(0,KLowestGroLeafNode_brch(NB,NZ)-1)
  D575: DO N=1,NumRemobLeafNodes
    SNCT=SNCX/XKSNC
    D650: DO KK=KN,KHiestGroLeafNode_brch(NB,NZ)
      SNCLF=0._r8
      SNCSH=0._r8
      K=pMOD(KK,MaxNodesPerBranch1)
!
!       REMOBILIZATION OF LEAF C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
!
!       WGLF,PetioleElmntNode_brch=node leaf,petiole C mass
    !       SCNF,SCNSH=leaf,petiole senescence respiration
    !       ElmntRemobFallingLeaf(ielmc),ElmntRemobFallingLeaf(ielmn),ElmntRemobFallingLeaf(ielmp)=remobilization of C,N,P from senescing leaf
    !       RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
    !       RCCZ,RCCY=min,max fractions for shoot C recycling
    !
      IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.ZEROP(NZ))THEN
        FNCLF=LeafElmntNode_brch(ielmc,K,NB,NZ)/(LeafElmntNode_brch(ielmc,K,NB,NZ)+&
          PetioleElmntNode_brch(ielmc,K,NB,NZ))
        SNCLF=FNCLF*SNCT
        SNCSH=SNCT-SNCLF
        ElmntRemobFallingLeaf(ielmc)=RCCC*LeafElmntNode_brch(ielmc,K,NB,NZ)
        ElmntRemobFallingLeaf(ielmn)=LeafElmntNode_brch(ielmn,K,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
        ElmntRemobFallingLeaf(ielmp)=LeafElmntNode_brch(ielmp,K,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
!
    !         FRACTION OF CURRENT LEAF TO BE REMOBILIZED
    !
    !         FSNCL,FSNAL=fraction of current leaf C,area to be remobilized
    !
        IF(ElmntRemobFallingLeaf(ielmc).GT.ZEROP(NZ))THEN
          FSNCL=AZMAX1(AMIN1(1.0_r8,SNCLF/ElmntRemobFallingLeaf(ielmc)))
        ELSE
          FSNCL=1.0_r8
        ENDIF
        FSNAL=FSNCL
!
        !     NON-REMOBILIZABLE C,N,P BECOMES LitrFall ALLOCATED
        !     TO FRACTIONS SET IN 'STARTQ'
        !
        !     CSNC,ZSNC,PSNC=literfall C,N,P
        !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
        !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
        !     FSNCL=fraction of current leaf to be remobilized
        !     WGLF,WGLFN,WGLFP=node leaf C,N,P mass
        !     ElmntRemobFallingLeaf(ielmc),ElmntRemobFallingLeaf(ielmn),ElmntRemobFallingLeaf(ielmp)=remobilization of C,N,P from senescing leaf
        !     FWODB=C woody fraction in other organs:0=woody,1=non-woody
        !     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
        !     CPOOL3_node,CPOOL4_node=C4 nonstructural C mass in bundle sheath,mesophyll
        !
        
        D6310: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
              +CFOPE(NE,icwood,M,NZ)*FSNCL*(LeafElmntNode_brch(NE,K,NB,NZ) &
              -ElmntRemobFallingLeaf(NE))*FWODLE(NE,k_woody_litr)

            LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,ifoliar,M,NZ)*FSNCL*(LeafElmntNode_brch(NE,K,NB,NZ) &
              -ElmntRemobFallingLeaf(NE))*FWODLE(NE,k_fine_litr)
          ENDDO
        ENDDO D6310
        IF(K.NE.0)THEN
          LitfalStrutElms_pvr(ielmc,icarbhyro,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(ielmc,icarbhyro,k_fine_litr,0,NZ) &
            +FSNCL*(CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ))
          CPOOL3_node(K,NB,NZ)=CPOOL3_node(K,NB,NZ)-FSNCL*CPOOL3_node(K,NB,NZ)
          CPOOL4_node(K,NB,NZ)=CPOOL4_node(K,NB,NZ)-FSNCL*CPOOL4_node(K,NB,NZ)
        ENDIF
!
!     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
!
!     LeafAreaLive_brch=leaf area
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WGLF,WGLFN,WGLFP,LeafProteinCNode_brch=node leaf C,N,P,protein mass
!     FSNCL=fraction of current leaf to be remobilized
!     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!
        LeafAreaLive_brch(NB,NZ)=AZMAX1(LeafAreaLive_brch(NB,NZ)-FSNAL*LeafAreaNode_brch(K,NB,NZ))
        LeafAreaNode_brch(K,NB,NZ)=LeafAreaNode_brch(K,NB,NZ)-FSNAL*LeafAreaNode_brch(K,NB,NZ)

        DO NE=1,NumPlantChemElms
          LeafStrutElms_brch(NE,NB,NZ)=AZMAX1(LeafStrutElms_brch(NE,NB,NZ)-FSNCL*LeafElmntNode_brch(NE,K,NB,NZ))
          LeafElmntNode_brch(NE,K,NB,NZ)=LeafElmntNode_brch(NE,K,NB,NZ)-FSNCL*LeafElmntNode_brch(NE,K,NB,NZ)
        ENDDO
        LeafProteinCNode_brch(K,NB,NZ)=AZMAX1(LeafProteinCNode_brch(K,NB,NZ) &
          -FSNCL*AMAX1(LeafElmntNode_brch(ielmn,K,NB,NZ)*rCNNonstructRemob_pft(NZ) &
          ,LeafElmntNode_brch(ielmp,K,NB,NZ)*rCPNonstructRemob_pft(NZ)))
!
!     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
!     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!     FSNCL=fraction of current leaf C to be remobilized
!     ElmntRemobFallingLeaf(ielmc),ElmntRemobFallingLeaf(ielmn),ElmntRemobFallingLeaf(ielmp)=remobilization of C,N,P from senescing leaf
!     SNCLF,SNCT=remaining senescence respiration carried to next node
!
        CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)+FSNCL*ElmntRemobFallingLeaf(ielmc)*SNCF
        DO NE=2,NumPlantChemElms
          CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+FSNCL*ElmntRemobFallingLeaf(NE)
          if(CanopyNonstElms_brch(NE,NB,NZ)<0._r8)then
            write(*,*)'1077CanopyNonstElms_brch(NE,NB,NZ)',NE,NB,CanopyNonstElms_brch(NE,NB,NZ)&
              -FSNCL*ElmntRemobFallingLeaf(NE),FSNCL*ElmntRemobFallingLeaf(NE)
            stop
          endif
        ENDDO
        SNCLF=SNCLF-FSNCL*ElmntRemobFallingLeaf(ielmc)
        SNCT=SNCT-FSNCL*ElmntRemobFallingLeaf(ielmc)
        IF(LeafStrutElms_brch(ielmc,NB,NZ).LE.ZEROL(NZ))THEN
          LeafStrutElms_brch(ielmc,NB,NZ)=0._r8
          LeafAreaLive_brch(NB,NZ)=0._r8
        ENDIF
      !
        !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
        !
!        IF(SNCLF.LE.ZEROP(NZ))GO TO 564
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
            LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,icwood,M,NZ) &
              *LeafElmntNode_brch(NE,K,NB,NZ)*FWODLE(NE,k_woody_litr)

            LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,ifoliar,M,NZ) &
              *LeafElmntNode_brch(NE,K,NB,NZ)*FWODLE(NE,k_fine_litr)
          ENDDO
        ENDDO D6315
        IF(K.NE.0)THEN
          LitfalStrutElms_pvr(ielmc,icarbhyro,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(ielmc,icarbhyro,k_fine_litr,0,NZ)&
            +CPOOL3_node(K,NB,NZ)+CPOOL4_node(K,NB,NZ)
          CPOOL3_node(K,NB,NZ)=0._r8
          CPOOL4_node(K,NB,NZ)=0._r8
        ENDIF
        LeafAreaLive_brch(NB,NZ)=AZMAX1(LeafAreaLive_brch(NB,NZ)-LeafAreaNode_brch(K,NB,NZ))
        DO NE=1,NumPlantChemElms
          LeafStrutElms_brch(NE,NB,NZ)=AZMAX1(LeafStrutElms_brch(NE,NB,NZ)-LeafElmntNode_brch(NE,K,NB,NZ))
          LeafElmntNode_brch(NE,K,NB,NZ)=0._r8
        ENDDO
        LeafAreaNode_brch(K,NB,NZ)=0._r8
        LeafProteinCNode_brch(K,NB,NZ)=0._r8
        IF(LeafStrutElms_brch(ielmc,NB,NZ).LE.ZEROL(NZ))THEN
          LeafStrutElms_brch(ielmc,NB,NZ)=0._r8
          LeafAreaLive_brch(NB,NZ)=0._r8
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
      IF(PetioleElmntNode_brch(ielmc,K,NB,NZ).GT.ZEROP(NZ))THEN
        RCES(ielmc)=PetioleElmntNode_brch(ielmc,K,NB,NZ)*RCCC
        RCES(ielmn)=PetioleElmntNode_brch(ielmn,K,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
        RCES(ielmp)=PetioleElmntNode_brch(ielmp,K,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
!
      !     FRACTION OF REMOBILIZATION THAT CAN BE MET FROM CURRENT SHEATH
      !     OR PETIOLE
      !
      !     FSNCS,FSNAS=fraction of current petiole C,length to be remobilized
      !
        IF(RCES(ielmc).GT.ZEROP(NZ))THEN
          FSNCS=AZMAX1(AMIN1(1.0_r8,SNCSH/RCES(ielmc)))
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
            LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,icwood,M,NZ) &
              *FSNCS*(PetioleElmntNode_brch(NE,K,NB,NZ)-RCES(NE))*FWODBE(NE,k_woody_litr)
            LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,inonfoliar,M,NZ) &
              *FSNCS*(PetioleElmntNode_brch(NE,K,NB,NZ)-RCES(NE))*FWODBE(NE,k_fine_litr)
          ENDDO    
        ENDDO D6320
        DO NE=1,NumPlantChemElms
          PetoleStrutElms_brch(NE,NB,NZ)=AZMAX1(PetoleStrutElms_brch(NE,NB,NZ)-FSNCS*PetioleElmntNode_brch(NE,K,NB,NZ))
          PetioleElmntNode_brch(NE,K,NB,NZ)=PetioleElmntNode_brch(NE,K,NB,NZ)-FSNCS*PetioleElmntNode_brch(NE,K,NB,NZ)
        ENDDO
!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
      !
      !     PetioleLengthNode_brch=petiole length
      !     WTSHEB,WTLFBN,WTSHBP=branch petiole C,N,P mass
      !     PetioleElmntNode_brch,WGSHN,WGSHP,PetioleProteinCNode_brch=node petiole C,N,P,protein mass
      !     FSNCS=fraction of current petiole to be remobilized
      !     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
      !

        PetioleLengthNode_brch(K,NB,NZ)=PetioleLengthNode_brch(K,NB,NZ)-FSNAS*PetioleLengthNode_brch(K,NB,NZ)
        PetioleProteinCNode_brch(K,NB,NZ)=AZMAX1(PetioleProteinCNode_brch(K,NB,NZ) &
          -FSNCS*AMAX1(PetioleElmntNode_brch(ielmn,K,NB,NZ)*rCNNonstructRemob_pft(NZ) &
          ,PetioleElmntNode_brch(ielmp,K,NB,NZ)*rCPNonstructRemob_pft(NZ)))
!
!     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
  !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
  !
  !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
  !     FSNCS=fraction of current petiole C to be remobilized
  !     RCES(ielmc),RCES(ielmn),RCES(ielmp)=remobilization of C,N,P from senescing petiole
  !     SNCSH,SNCT=remaining senescence respiration carried to next node
  !
        CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)+FSNCS*RCES(ielmc)*SNCF
        DO NE=2,NumPlantChemElms
          CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+FSNCS*RCES(NE)
        ENDDO
        SNCSH=SNCSH-FSNCS*RCES(ielmc)
        SNCT=SNCT-FSNCS*RCES(ielmc)
        IF(PetoleStrutElms_brch(ielmc,NB,NZ).LE.ZEROL(NZ))THEN
          PetoleStrutElms_brch(ielmc,NB,NZ)=0._r8
        ENDIF
!
  !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
  !
        IF(SNCSH.LE.ZEROP(NZ))GO TO 565
  !
      !     OTHERWISE REMAINING C,N,P IN SHEATH OR PETIOLE GOES TO LitrFall
      !
      !     CSNC,ZSNC,PSNC=literfall C,N,P
      !     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
      !     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
      !     FWODB=C woody fraction in branch:0=woody,1=non-woody
      !     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
      !     PetioleLengthNode_brch=petiole length
      !     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
      !     PetioleElmntNode_brch,WGSHN,WGSHP,PetioleProteinCNode_brch=node petiole C,N,P,protein mass
      !
      ELSE
        
        D6325: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,icwood,M,NZ) &
              *PetioleElmntNode_brch(NE,K,NB,NZ)*FWODBE(NE,k_woody_litr)
            LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,inonfoliar,M,NZ) &
              *PetioleElmntNode_brch(NE,K,NB,NZ)*FWODBE(NE,k_fine_litr)
          ENDDO    
        ENDDO D6325
        DO NE=1,NumPlantChemElms
          PetoleStrutElms_brch(NE,NB,NZ)=AZMAX1(PetoleStrutElms_brch(NE,NB,NZ)-PetioleElmntNode_brch(NE,K,NB,NZ))
          PetioleElmntNode_brch(NE,K,NB,NZ)=0._r8
        ENDDO
        PetioleLengthNode_brch(K,NB,NZ)=0._r8
        PetioleProteinCNode_brch(K,NB,NZ)=0._r8
        IF(PetoleStrutElms_brch(ielmc,NB,NZ).LE.ZEROL(NZ))THEN
          PetoleStrutElms_brch(ielmc,NB,NZ)=0._r8
        ENDIF
      ENDIF
    ENDDO D650
    KN=KN+1
    RMxess_brch=SNCT*(1.0_r8-SNCF)
!
!     REMOBILIZATION OF RESERVE C
!
!     WTRSVB=stalk reserve C mass
!     RMxess_brch=excess maintenance respiration
!
    IF(StalkRsrvElms_brch(ielmc,NB,NZ).GT.RMxess_brch)THEN
      StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)-RMxess_brch
      write(*,*)'1273StalkRsrvElms_brch(ielmc,NB,NZ)',StalkRsrvElms_brch(ielmc,NB,NZ),RMxess_brch
      RMxess_brch=0._r8
      go to 565
    ENDIF
!
!     REMOBILIZATION OF STALK C,N,P
!
!     FXFS=rate constant for remobilization of stalk C,N,P (h-1)
!     SNCZ=phenologically-driven respiration senescence during late-season
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     WTSTKB,StalkBiomassC_brch=stalk,sapwood C mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     MXNOD,MNNOD=max,min node number currently growing
!     NumCogrowNode=number of concurrently growing nodes
!     KHiestGroLeafNode_brch=integer of most recent leaf number
!
    SNCZ=FXFS*RMxess_brch
    SNCT=RMxess_brch+SNCZ
    IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual .AND. SNCT.GT.ZEROP(NZ) &
      .AND.StalkStrutElms_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      SNCF=SNCZ/SNCT
      FRCC=StalkBiomassC_brch(NB,NZ)/StalkStrutElms_brch(ielmc,NB,NZ)
      RCSC=RCCC*FRCC
      RCSN=RCCN*FRCC
      RCSP=RCCP*FRCC
      MXNOD=KHiestGroLeafNode_brch(NB,NZ)
      MNNOD=MAX(MIN(0,MAX(0,MXNOD-NumCogrowNode(NZ))),KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+2)
      MXNOD=MAX(MXNOD,MNNOD)
      D1650: DO KK=MXNOD,MNNOD,-1
        K=pMOD(KK,MaxNodesPerBranch1)
    !
    !     REMOBILIZATION OF STALK C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
    !
    !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
    !     InternodeStrutElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
    !
        IF(InternodeStrutElms_brch(ielmc,K,NB,NZ).GT.ZEROP(NZ))THEN
          RCEK(ielmc)=InternodeStrutElms_brch(ielmc,K,NB,NZ)*RCSC
          RCEK(ielmn)=InternodeStrutElms_brch(ielmn,K,NB,NZ)*(RCSN+(1.0_r8-RCSN)*RCSC)
          RCEK(ielmp)=InternodeStrutElms_brch(ielmp,K,NB,NZ)*(RCSP+(1.0_r8-RCSP)*RCSC)
    !
      !     FRACTION OF CURRENT NODE TO BE REMOBILIZED
      !
      !     FSNCS=fraction of lowest internode to be remobilized
!
          IF(RCEK(ielmc).GT.ZEROP(NZ))THEN
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
              LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,istalk,M,NZ) &
                *FSNCK*(InternodeStrutElms_brch(NE,K,NB,NZ)-RCEK(NE))*FWOODE(NE,k_woody_litr)

              LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+CFOPE(NE,istalk,M,NZ) &
                *FSNCK*(InternodeStrutElms_brch(NE,K,NB,NZ)-RCEK(NE))*FWOODE(NE,k_fine_litr)
            ENDDO    
          ENDDO D7310

!
      !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
      !
      !     FSNCK=fraction of lowest internode to be remobilized
      !     InternodeHeightLive_brch,InternodeHeightDying_brch=living,senescing internode length
      !     WTSTKB,WTSTBN,WTSTBP,InternodeStrutElms_brch,WGNODN,WGNODP=C,N,P mass in senescing internode
      !
          DO NE=1,NumPlantChemElms
            StalkStrutElms_brch(NE,NB,NZ)=AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-FSNCK*InternodeStrutElms_brch(NE,K,NB,NZ))
            InternodeStrutElms_brch(NE,K,NB,NZ)=InternodeStrutElms_brch(NE,K,NB,NZ)-FSNCK*InternodeStrutElms_brch(NE,K,NB,NZ)
          ENDDO
          InternodeHeightLive_brch(K,NB,NZ)=InternodeHeightLive_brch(K,NB,NZ)-FSNCK*InternodeHeightDying_brch(K,NB,NZ)
          InternodeHeightDying_brch(K,NB,NZ)=InternodeHeightDying_brch(K,NB,NZ)-FSNCK*InternodeHeightDying_brch(K,NB,NZ)
!
      !     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
      !     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
      !
      !     RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
      !     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
      !     FSNCK=fraction of lowest internode to be remobilized
      !     SNCT=remaining node senescence respiration
      !
          write(*,*)'1367StalkRsrvElms_brch(ielmc,NB,NZ)',StalkRsrvElms_brch(ielmc,NB,NZ),FSNCK*RCEK(ielmc)*SNCF      
          StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+FSNCK*RCEK(ielmc)*SNCF
          DO NE=2,NumPlantChemElms
            StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+FSNCK*RCEK(NE)
          ENDDO
          SNCT=SNCT-FSNCK*RCEK(ielmc)
    !
      !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
      !
          IF(SNCT.LE.ZEROP(NZ))GO TO 565
    !
    !       OTHERWISE REMAINING C,N,P IN NODE GOES TO LitrFall
    !
    !       CSNC,ZSNC,PSNC=literfall C,N,P
    !       CFOPC=fraction of plant litter allocated in nonstructural(0,*),
    !       WTSTKB,WTSTBN,WTSTBP,InternodeStrutElms_brch,WGNODN,WGNODP=C,N,P mass in senescing internode
    !       InternodeHeightLive_brch,InternodeHeightDying_brch=living,senescing internode length
    !
        ELSE
          
          D7315: DO M=1,jsken
            DO NE=1,NumPlantChemElms
              LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
                +CFOPE(NE,istalk,M,NZ)*InternodeStrutElms_brch(NE,K,NB,NZ)*FWOODE(NE,k_woody_litr)

              LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
                +CFOPE(NE,istalk,M,NZ)*InternodeStrutElms_brch(NE,K,NB,NZ)*FWOODE(NE,k_fine_litr)
            ENDDO    
          ENDDO D7315
          DO NE=1,NumPlantChemElms
            StalkStrutElms_brch(NE,NB,NZ)=AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-InternodeStrutElms_brch(NE,K,NB,NZ))
            InternodeStrutElms_brch(NE,K,NB,NZ)=0._r8
          ENDDO
          InternodeHeightLive_brch(K,NB,NZ)=InternodeHeightLive_brch(K,NB,NZ)-InternodeHeightDying_brch(K,NB,NZ)
          InternodeHeightDying_brch(K,NB,NZ)=0._r8
        ENDIF
      ENDDO D1650
!
    !   RESIDUAL STALK
    !
    !   RCEK(ielmc),RCEK(ielmn),RCEK(ielmp)=remobilization of C,N,P from senescing internode
    !   WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !
      IF(SenecStalkStrutElms_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
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
        IF(RCEK(ielmc).GT.ZEROP(NZ))THEN
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
            LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
              +CFOPE(NE,istalk,M,NZ)*FSNCR*(SenecStalkStrutElms_brch(NE,NB,NZ)-RCEK(NE))*FWOODE(NE,k_woody_litr)

            LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
              +CFOPE(NE,istalk,M,NZ)*FSNCR*(SenecStalkStrutElms_brch(NE,NB,NZ)-RCEK(NE))*FWOODE(NE,k_fine_litr)
          ENDDO    
        ENDDO D8310

    !
    !     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
    !
    !     FSNCR=fraction of residual stalk to be remobilized
    !     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in senescing stalk
    !     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
    !     InternodeHeightLive_brch,InternodeHeightDying_brch=living,senescing internode length
    !
        DO NE=1,NumPlantChemElms
          FStalkSenes=FSNCR*SenecStalkStrutElms_brch(NE,NB,NZ)
          StalkStrutElms_brch(NE,NB,NZ)=AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-FStalkSenes)

          SenecStalkStrutElms_brch(NE,NB,NZ)=AZMAX1(SenecStalkStrutElms_brch(NE,NB,NZ)-FStalkSenes)
        ENDDO
        HTNODZ=0._r8
        D8320: DO K=0,MaxNodesPerBranch1
          HTNODZ=AMAX1(HTNODZ,InternodeHeightLive_brch(K,NB,NZ))
        ENDDO D8320

        HTNODZ=AZMAX1(HTNODZ-FSNCR*HTNODZ)
        D8325: DO K=0,MaxNodesPerBranch1
          InternodeHeightLive_brch(K,NB,NZ)=AMIN1(HTNODZ,InternodeHeightLive_brch(K,NB,NZ))
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
        StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+FSNCR*RCEK(ielmc)*SNCF
        DO NE=2,NumPlantChemElms
          StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+FSNCR*RCEK(NE)
        ENDDO
        SNCT=SNCT-FSNCR*RCEK(ielmc)
      ENDIF
  !
  !     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
  !
      IF(SNCT.LE.ZEROP(NZ))go to 565
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
          LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)&
            +CFOPE(NE,istalk,M,NZ)*SenecStalkStrutElms_brch(NE,NB,NZ)*FWOODE(NE,k_woody_litr)

          LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
            +CFOPE(NE,istalk,M,NZ)*SenecStalkStrutElms_brch(NE,NB,NZ)*FWOODE(NE,k_fine_litr)
        ENDDO    
      ENDDO D8315
      DO NE=1,NumPlantChemElms
        StalkStrutElms_brch(NE,NB,NZ)=AZMAX1(StalkStrutElms_brch(NE,NB,NZ)-SenecStalkStrutElms_brch(NE,NB,NZ))
        SenecStalkStrutElms_brch(NE,NB,NZ)=0._r8
      ENDDO
    ENDIF
!
!     REMOBILIZATION OF STORAGE C,N,P
!
!     WTRVC=storage C
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     RMxess_brch=remaining excess maintenance respiration
!
    RMxess_brch=SNCT/(1.0_r8+FXFS)
    write(*,*)'remob1503NonStrutElms_pft(ielmc,NZ)',NZ,NonStrutElms_pft(ielmc,NZ)
    IF(NonStrutElms_pft(ielmc,NZ).GT.RMxess_brch)THEN
      NonStrutElms_pft(ielmc,NZ)=NonStrutElms_pft(ielmc,NZ)-RMxess_brch
      RMxess_brch=0._r8
    ELSEIF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
      iPlantBranchState_brch(NB,NZ)=iDead
    ENDIF
565 CONTINUE
  ENDDO D575
  end associate
  end subroutine RemobilizeLeafLayers


!------------------------------------------------------------------------------------------

  subroutine AllocateLeafToCanopyLayers(NB,NZ,CanopyHeight_copy)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8), intent(in) :: CanopyHeight_copy(JP1)
  integer  :: LL,LU,L,K,k1,k2,KK,NE
  integer  :: KMinGroingLeafNodeNum,KLeafNumHighestGrowing
  integer  :: LNumHeightLeafTip,LNumHeightLeafBase
  integer  :: LHTBRU,LHTBRL,N
  real(r8) :: ZSTK
  real(r8) :: YLeafElmntNode_brch(NumPlantChemElms)
  real(r8) :: YLeafArea_brch,LeafElevation,LeafLength
  real(r8) :: ARSTKB,ASTV
  real(r8) :: FRACL
  real(r8) :: HeightBranchBase
  real(r8) :: HeightStalk
  real(r8) :: HeightLeafNode,HeightLeafLow,HeightLeafTip
  real(r8) :: HeightLeafBase
  real(r8) :: RSTK
  real(r8) :: TLGLF
! begin_execution
  associate(                            &
    LeafChemElmByLayerNode_brch  =>  plt_biom%LeafChemElmByLayerNode_brch  , &
    LeafElmntNode_brch           =>  plt_biom%LeafElmntNode_brch   , &
    CanopyLeafCLyr_pft           =>  plt_biom%CanopyLeafCLyr_pft   , &
    StalkBiomassC_brch           =>  plt_biom%StalkBiomassC_brch  , &
    StalkStrutElms_brch           =>  plt_biom%StalkStrutElms_brch , &
    PetioleProteinCNode_brch     =>  plt_biom%PetioleProteinCNode_brch  , &
    FNOD                         =>  plt_allom%FNOD   , &
    iPlantRootProfile_pft        =>  plt_pheno%iPlantRootProfile_pft , &
    iPlantPhenolPattern_pft   =>  plt_pheno%iPlantPhenolPattern_pft , &
    KHiestGroLeafNode_brch       =>  plt_pheno%KHiestGroLeafNode_brch , &
    iPlantTurnoverPattern_pft    =>  plt_pheno%iPlantTurnoverPattern_pft , &
    KLowestGroLeafNode_brch      =>  plt_pheno%KLowestGroLeafNode_brch , &
    WDLF                         =>  plt_morph%WDLF    , &
    SeedDepth_pft                =>  plt_morph%SeedDepth_pft   , &
    MainBranchNum_pft            =>  plt_morph%MainBranchNum_pft     , &
    CanopyLeafAreaByLayer_pft    =>  plt_morph%CanopyLeafAreaByLayer_pft   , &
    PetioleLengthNode_brch       =>  plt_morph%PetioleLengthNode_brch   , &
    CanopyStemArea_lbrch          =>  plt_morph%CanopyStemArea_lbrch   , &
    BranchNumber_brch            =>  plt_morph%BranchNumber_brch    , &
    InternodeHeightLive_brch     =>  plt_morph%InternodeHeightLive_brch  , &
    CanopyHeightz_col            =>  plt_morph%CanopyHeightz_col      , &
    CanopyHeight_pft             =>  plt_morph%CanopyHeight_pft     , &
    HypoctoHeight_pft            =>  plt_morph%HypoctoHeight_pft   , &
    CLASS                        =>  plt_morph%CLASS   , &
    LeafAreaNode_brch            =>  plt_morph%LeafAreaNode_brch   , &
    CanopyLeafALyr_pft           =>  plt_morph%CanopyLeafALyr_pft   , &
    PlantPopulation_pft          =>  plt_site%PlantPopulation_pft       , &
    ZERO                         =>  plt_site%ZERO     , &
    SineLeafAngle                =>  plt_rad%SineLeafAngle        &
  )
!   ALLOCATION OF LEAF AREA TO CANOPY LAYERS
!
!   HypoctoHeight_pft=hypocotyledon height
!   SeedDepth_pft=seeding depth
!   LeafAreaNode_brch=node leaf area
!   PetioleLengthNode_brch=petiole length
!
  KLowestGroLeafNode_brch(NB,NZ)=0
  IF(HypoctoHeight_pft(NZ).LE.SeedDepth_pft(NZ).AND. &
    LeafAreaNode_brch(0,MainBranchNum_pft(NZ),NZ).GT.0.0_r8)THEN
    LeafLength=SQRT(1.0E+02_r8*LeafAreaNode_brch(0,MainBranchNum_pft(NZ),NZ)/PlantPopulation_pft(NZ))
    HypoctoHeight_pft(NZ)=LeafLength+PetioleLengthNode_brch(0,MainBranchNum_pft(NZ),NZ) &
      +InternodeHeightLive_brch(0,MainBranchNum_pft(NZ),NZ)
  ENDIF
!
! IF CANOPY HAS EMERGED
!
  IF(HypoctoHeight_pft(NZ).GT.SeedDepth_pft(NZ))THEN
    D540: DO K=0,MaxNodesPerBranch1
      DO  L=1,NumOfCanopyLayers1
        CanopyLeafAreaByLayer_pft(L,K,NB,NZ)=0._r8
        LeafChemElmByLayerNode_brch(1:NumPlantChemElms,L,K,NB,NZ)=0._r8
      enddo
    ENDDO D540
    D535: DO L=1,NumOfCanopyLayers1
      CanopyStemArea_lbrch(L,NB,NZ)=0._r8
    ENDDO D535
!
!   BRANCH HEIGHT
!
!   iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+petiole,2=none,3=between 1,2
!   iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!   KLeafNumHighestGrowing,KLowestGroLeafNode_brch=integer of highest,lowest leaf number currently growing
!   InternodeHeightLive_brch=internode length
!   HeightBranchBase=branch base height
!
    IF(iPlantTurnoverPattern_pft(NZ).NE.0.AND.is_plant_treelike(iPlantRootProfile_pft(NZ)))THEN
      IF(NB.NE.MainBranchNum_pft(NZ))THEN
        KLeafNumHighestGrowing=MAX(1,KHiestGroLeafNode_brch(MainBranchNum_pft(NZ),NZ)-MaxNodesPerBranch1+1)
        IF(BranchNumber_brch(NB,NZ).GE.KLeafNumHighestGrowing)THEN
          K=pMOD(BranchNumber_brch(NB,NZ),MaxNodesPerBranch1)
          HeightBranchBase=InternodeHeightLive_brch(K,MainBranchNum_pft(NZ),NZ)
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
      !     InternodeHeightLive_brch=internode length
      !     HeightLeafNode=leaf node height
      !     LeafAreaNode_brch=leaf node area
      !     PP=plant population
      !     FNOD=scales node number for perennial vegetation (e.g. trees)
      !     LeafLength=leaf length
!
      HeightStalk=HeightBranchBase+InternodeHeightLive_brch(K,NB,NZ)
      HeightLeafNode=HeightStalk+PetioleLengthNode_brch(K,NB,NZ)
      LeafLength=AZMAX1(SQRT(WDLF(NZ)*AZMAX1(LeafAreaNode_brch(K,NB,NZ))/(PlantPopulation_pft(NZ)*FNOD(NZ))))
      TLGLF=0._r8
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
        LeafElevation=SineLeafAngle(N)*CLASS(N,NZ)*LeafLength
        HeightLeafLow=AMIN1(CanopyHeight_copy(NZ)+0.01_r8-LeafElevation,HeightLeafNode+TLGLF)
        HeightLeafTip=AMIN1(CanopyHeight_copy(NZ)+0.01_r8,HeightLeafLow+LeafElevation)
        LU=0
        LL=0
        D550: DO L=NumOfCanopyLayers1,1,-1
          IF(LU.EQ.1 .AND. LL.EQ.1)exit
          IF((HeightLeafTip.GT.CanopyHeightz_col(L-1).OR.CanopyHeightz_col(L-1).LE.ZERO).AND.LU.EQ.0)THEN
            LNumHeightLeafTip=MAX(1,L)
            LU=1
          ENDIF
          IF((HeightLeafLow.GT.CanopyHeightz_col(L-1).OR.CanopyHeightz_col(L-1).LE.ZERO).AND.LL.EQ.0)THEN
            LNumHeightLeafBase=MAX(1,L)
            LL=1
          ENDIF
        ENDDO D550
        D570: DO L=LNumHeightLeafBase,LNumHeightLeafTip
          IF(LNumHeightLeafTip.EQ.LNumHeightLeafBase)THEN
            FRACL=CLASS(N,NZ)
          ELSEIF(HeightLeafTip.GT.HeightLeafLow.AND.CanopyHeightz_col(L).GT.HeightLeafLow)THEN
            FRACL=CLASS(N,NZ)*(AMIN1(HeightLeafTip,CanopyHeightz_col(L)) &
              -AMAX1(HeightLeafLow,CanopyHeightz_col(L-1)))/(HeightLeafTip-HeightLeafLow)
          ELSE
            FRACL=CLASS(N,NZ)
          ENDIF
          YLeafArea_brch=FRACL*LeafAreaNode_brch(K,NB,NZ)
          DO NE=1,NumPlantChemElms
            YLeafElmntNode_brch(NE)=FRACL*LeafElmntNode_brch(NE,K,NB,NZ)
          ENDDO
!
    !     ACCUMULATE LAYER LEAF AREAS, C, N AND P CONTENTS
    !
    !     CanopyLeafAreaByLayer_pft=leaf node area in canopy layer
    !     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
    !     CanopyLeafALyr_pft,CanopyLeafCLyr_pft=total leaf area,C in canopy layer
    !     InternodeHeightLive_brch=internode length
    !
          CanopyLeafAreaByLayer_pft(L,K,NB,NZ)=CanopyLeafAreaByLayer_pft(L,K,NB,NZ)+YLeafArea_brch
          DO NE=1,NumPlantChemElms
            LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)=LeafChemElmByLayerNode_brch(NE,L,K,NB,NZ)+YLeafElmntNode_brch(NE)
          ENDDO
          CanopyLeafALyr_pft(L,NZ)=CanopyLeafALyr_pft(L,NZ)+YLeafArea_brch
          CanopyLeafCLyr_pft(L,NZ)=CanopyLeafCLyr_pft(L,NZ)+YLeafElmntNode_brch(ielmc)
        ENDDO D570
        TLGLF=TLGLF+LeafElevation
        CanopyHeight_pft(NZ)=AMAX1(CanopyHeight_pft(NZ),HeightLeafTip)
      ENDDO D555
      IF(PetioleProteinCNode_brch(K,NB,NZ).GT.0.0_r8)THEN
        IF(KLowestGroLeafNode_brch(NB,NZ).EQ.0)KLowestGroLeafNode_brch(NB,NZ)=MIN(KK,KHiestGroLeafNode_brch(NB,NZ))
      ENDIF
    ENDDO D560

    IF(KLowestGroLeafNode_brch(NB,NZ).EQ.0)KLowestGroLeafNode_brch(NB,NZ)=KHiestGroLeafNode_brch(NB,NZ)
    K1=pMOD(KHiestGroLeafNode_brch(NB,NZ),MaxNodesPerBranch1)
    K2=pMOD(KHiestGroLeafNode_brch(NB,NZ)-1,MaxNodesPerBranch1)
    IF(isclose(InternodeHeightLive_brch(K1,NB,NZ),0._r8))THEN
      InternodeHeightLive_brch(K1,NB,NZ)=InternodeHeightLive_brch(K2,NB,NZ)
    ENDIF
    HeightLeafBase=HeightBranchBase+AZMAX1(InternodeHeightLive_brch(K1,NB,NZ))
!
  !     ALLOCATE STALK SURFACE AREA TO CANOPY LAYERS
  !
  !     InternodeHeightLive_brch=internode length
  !     HeightLeafBase=leaf base height
  !     ZL=height to bottom of each canopy layer
  !     LHTBRL,LHTBRU=layer number of branch base,tip
  !     WTSTKB,ARSTKB=branch stalk mass,surface area
  !     FSTK=fraction of stalk area contributing to water,heat flow
  !     DSTK,VSTK=stalk density (Mg m-3),specific volume (m3 g-1)
  !     StalkBiomassC_brch=stalk sapwood mass
  !     FRACL=stalk fraction in each layer
  !     CanopyStemArea_lbrch=total branch stalk surface area in each layer
  !
  !     IF(NZ.EQ.1)THEN
  !     WRITE(*,6679)'K1',I,J,NZ,NB,K1,KHiestGroLeafNode_brch(NB,NZ)
  !    2,InternodeHeightLive_brch(K1,NB,NZ)
  !6679  FORMAT(A8,6I4,12E12.4)
  !     ENDIF

    IF(InternodeHeightLive_brch(K1,NB,NZ).GT.0.0_r8)THEN
      LU=0
      LL=0
      D545: DO L=NumOfCanopyLayers1,1,-1
        IF(LU.EQ.1 .AND. LL.EQ.1)exit
        IF((HeightLeafBase.GT.CanopyHeightz_col(L-1).OR.CanopyHeightz_col(L-1).LE.ZERO).AND.LU.EQ.0)THEN
          LHTBRU=MAX(1,L)
          LU=1
        ENDIF
        IF((HeightBranchBase.GT.CanopyHeightz_col(L-1).OR.CanopyHeightz_col(L-1).LE.ZERO).AND.LL.EQ.0)THEN
          LHTBRL=MAX(1,L)
          LL=1
        ENDIF
      ENDDO D545
      
      RSTK=SQRT(VSTK*(AZMAX1(StalkStrutElms_brch(ielmc,NB,NZ))/PlantPopulation_pft(NZ)) &
        /(PICON*InternodeHeightLive_brch(K1,NB,NZ)))
      ARSTKB=PICON*InternodeHeightLive_brch(K1,NB,NZ)*PlantPopulation_pft(NZ)*RSTK
      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
        StalkBiomassC_brch(NB,NZ)=StalkStrutElms_brch(ielmc,NB,NZ)
      ELSE
        ZSTK=AMIN1(ZSTX,FSTK*RSTK)
        ASTV=PICON*(2.0_r8*RSTK*ZSTK-ZSTK*ZSTK)
        StalkBiomassC_brch(NB,NZ)=ASTV/VSTK*InternodeHeightLive_brch(K1,NB,NZ)*PlantPopulation_pft(NZ)
      ENDIF
    !     IF(NZ.EQ.1)THEN
        !     WRITE(*,6677)'WVSTK',I,J,NX,NY,NZ,NB,K1,StalkBiomassC_brch(NB,NZ)
        !    2,ASTV,VSTK,InternodeHeightLive_brch(K1,NB,NZ),PlantPopulation_pft(NZ)
        !6677  FORMAT(A8,7I4,12E12.4)
        !     ENDIF
      D445: DO L=LHTBRL,LHTBRU
        IF(HeightLeafBase.GT.HeightBranchBase)THEN
          IF(HeightLeafBase.GT.CanopyHeightz_col(L-1))THEN
            FRACL=(AMIN1(HeightLeafBase,CanopyHeightz_col(L))-AMAX1(HeightBranchBase &
              ,CanopyHeightz_col(L-1)))/(HeightLeafBase-HeightBranchBase)
          ELSE
            FRACL=0._r8
          ENDIF
        ELSE
          FRACL=1.0_r8
        ENDIF
        CanopyStemArea_lbrch(L,NB,NZ)=FRACL*ARSTKB
      ENDDO D445
    ELSE
      StalkBiomassC_brch(NB,NZ)=0._r8
      D450: DO L=1,NumOfCanopyLayers1
        CanopyStemArea_lbrch(L,NB,NZ)=0._r8
      ENDDO D450
    ENDIF
  ELSE
    StalkBiomassC_brch(NB,NZ)=0._r8
    D455: DO L=1,NumOfCanopyLayers1
      CanopyStemArea_lbrch(L,NB,NZ)=0._r8
    ENDDO D455
  ENDIF
  end associate
  end subroutine AllocateLeafToCanopyLayers
!------------------------------------------------------------------------------------------
  subroutine LeafClassAllocation(NB,NZ)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8) :: dangle
  integer :: L,K,N
  ! begin_execution
  associate(                             &
    SineBranchAngle_pft         =>  plt_morph%SineBranchAngle_pft    , &
    LeafAreaNode_brch           =>  plt_morph%LeafAreaNode_brch    , &
    CanopyStemArea_lbrch         =>  plt_morph%CanopyStemArea_lbrch    , &
    CLASS                       =>  plt_morph%CLASS    , &
    CanopyLeafAreaByLayer_pft   =>  plt_morph%CanopyLeafAreaByLayer_pft    , &
    StemAreaZsec_brch           =>  plt_morph%StemAreaZsec_brch    , &
    MainBranchNum_pft           =>  plt_morph%MainBranchNum_pft      , &
    LeafAreaZsec_brch           =>  plt_morph%LeafAreaZsec_brch       &
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
!     ARLFXB=ARLFXB+LeafAreaNode_brch(K,NB,NZ)
    IF(LeafAreaNode_brch(K,NB,NZ).GT.0.0_r8)THEN
      D700: DO L=NumOfCanopyLayers1,1,-1
!       ARLFXL=ARLFXL+CanopyLeafAreaByLayer_pft(L,K,NB,NZ)
        D800: DO N=1,NumOfLeafZenithSectors1
          LeafAreaZsec_brch(N,L,K,NB,NZ)=AZMAX1(CLASS(N,NZ)*0.25_r8*CanopyLeafAreaByLayer_pft(L,K,NB,NZ))
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
! CanopyStemArea_lbrch=total branch stalk surface area in each layer
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
    StemAreaZsec_brch(N,L,NB,NZ)=CanopyStemArea_lbrch(L,NB,NZ)/real(NumOfLeafZenithSectors1,r8)
  ENDDO D710
  end associate
  end subroutine LeafClassAllocation

!------------------------------------------------------------------------------------------

  subroutine GrainFilling(I,NB,NZ,GrowthGrain,GROSTKC)
  implicit none
  integer, intent(in) :: I,NB,NZ
  real(r8), intent(in) :: GrowthGrain(NumPlantChemElms),GROSTKC
  real(r8) :: ZPGRP,ZPGRN,ZNPGP,ZNPGN
  real(r8) :: MaxChemElmRsrv2Grain,ChemElmRsrv2Grain(NumPlantChemElms)
  real(r8) :: FGRNX
  real(r8) :: GRMXB
  real(r8) :: GROLM
  real(r8) :: GROLC
  real(r8) :: SeedSET
  integer :: NE
! begin_execution
  associate(                              &
    TCelciusCanopy_pft                    =>  plt_ew%TCelciusCanopy_pft        , &
    LeafPetoNonstElmConc_brch             =>  plt_biom%LeafPetoNonstElmConc_brch    , &
    GrainStrutElms_brch                    =>  plt_biom%GrainStrutElms_brch    , &
    StalkRsrvElms_brch                   =>  plt_biom%StalkRsrvElms_brch   , &
    ZEROP                                 =>  plt_biom%ZEROP     , &
    GrainSeedBiomCMean_brch               =>  plt_allom%GrainSeedBiomCMean_brch    , &
    CNGR                                  =>  plt_allom%CNGR     , &
    CPGR                                  =>  plt_allom%CPGR     , &
    fTgrowRootP_vr                        =>  plt_pheno%fTgrowRootP_vr     , &
    fTgrowCanP                            =>  plt_pheno%fTgrowCanP     , &
    iPlantPhenolType_pft                  =>  plt_pheno%iPlantPhenolType_pft   , &
    Hours4LeafOff_brch                    =>  plt_pheno%Hours4LeafOff_brch     , &
    HourFailGrainFill_brch                =>  plt_pheno%HourFailGrainFill_brch     , &
    HourReq4LeafOff_brch                  =>  plt_pheno%HourReq4LeafOff_brch    , &
    GrainFillRate25C_pft                  =>  plt_pheno%GrainFillRate25C_pft    , &
    SSTX                                  =>  plt_pheno%SSTX     , &
    iPlantPhenolPattern_pft               =>  plt_pheno%iPlantPhenolPattern_pft   , &
    HighTCLimtSeed_pft                    =>  plt_pheno%HighTCLimtSeed_pft     , &
    dReproNodeNumNormByMatG_brch          =>  plt_pheno%dReproNodeNumNormByMatG_brch   , &
    iPlantCalendar_brch                   =>  plt_pheno%iPlantCalendar_brch   , &
    TCelciusChill4Seed                    =>  plt_pheno%TCelciusChill4Seed     , &
    PlantPopulation_pft                   =>  plt_site%PlantPopulation_pft        , &
    MaxSeedNumPerSite_pft                 =>  plt_morph%MaxSeedNumPerSite_pft     , &
    MaxSeedCMass                          =>  plt_morph%MaxSeedCMass    , &
    MaxPotentSeedNumber_pft               =>  plt_morph%MaxPotentSeedNumber_pft     , &
    NGTopRootLayer_pft                    =>  plt_morph%NGTopRootLayer_pft      , &
    PotentialSeedSites_brch               =>  plt_morph%PotentialSeedSites_brch    , &
    iPlantGrainType_pft                   =>  plt_morph%iPlantGrainType_pft   , &
    SeedNumSet_brch                    =>  plt_morph%SeedNumSet_brch      &
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
!   TCelciusCanopy_pft=canopy temperature
!   TCelciusChill4Seed=chilling temperature for CO2 fixation, seed loss (oC)
!   HTC=high temperature threshold for grain number loss
!   FGRNX=loss of seed set
!   SSTX=sensitivity to TCelciusCanopy_pft> HTC,TCelciusCanopy_pft< TCelciusChill4Seedfrom startq.f (seeds oC-1)
!   SeedNumSet_brch=seed set number
!   PotentialSeedSites_brch=potential number of seed set sites
!   MaxSeedNumPerSite_pft=maximum seed number per MaxPotentSeedNumber_pft from PFT file
!   dReproNodeNumNormByMatG_brch=change in reproductive node number normalized for maturity group
!
  IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).NE.0.AND.iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ).EQ.0)THEN
    SeedSET=AMIN1(LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ)+SETC) &
      ,LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ)+SETN) &
      ,LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ)+SETP))

    IF(TCelciusCanopy_pft(NZ).LT.TCelciusChill4Seed(NZ))THEN
      IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
        FGRNX=SSTX(NZ)*(TCelciusChill4Seed(NZ)-TCelciusCanopy_pft(NZ))
      ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN
        FGRNX=SSTX(NZ)*(TCelciusChill4Seed(NZ)-TCelciusCanopy_pft(NZ))
      ELSE
        FGRNX=0._r8
      ENDIF
    ELSEIF(TCelciusCanopy_pft(NZ).GT.HighTCLimtSeed_pft(NZ))THEN
      IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
        FGRNX=SSTX(NZ)*(TCelciusCanopy_pft(NZ)-HighTCLimtSeed_pft(NZ))
      ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN
        FGRNX=SSTX(NZ)*(TCelciusCanopy_pft(NZ)-HighTCLimtSeed_pft(NZ))
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
      GRMXB=MaxSeedCMass(NZ)
      GrainSeedBiomCMean_brch(NB,NZ)=AMIN1(MaxSeedCMass(NZ),GrainSeedBiomCMean_brch(NB,NZ) &
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
!   fTgrowCanP=temperature function for canopy growth
!   fTgrowRootP_vr=temperature function for root growth
!
  IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0)THEN
    IF(GrainStrutElms_brch(ielmc,NB,NZ).GE.GrainSeedBiomCMean_brch(NB,NZ)*SeedNumSet_brch(NB,NZ))THEN
      GROLM=0._r8
    ELSEIF(iPlantGrainType_pft(NZ).EQ.igraintyp_abvgrnd)THEN
      GROLM=AZMAX1(GrainFillRate25C_pft(NZ)*SeedNumSet_brch(NB,NZ)*SQRT(fTgrowCanP(NZ)))
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
    IF(GrainStrutElms_brch(ielmn,NB,NZ).LT.ZPGRM*CNGR(NZ) &
      *GrainStrutElms_brch(ielmc,NB,NZ).OR.GrainStrutElms_brch(ielmp,NB,NZ).LT.ZPGRM &
      *CPGR(NZ)*GrainStrutElms_brch(ielmc,NB,NZ))THEN
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
    IF(StalkRsrvElms_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      ZNPGN=StalkRsrvElms_brch(ielmn,NB,NZ)/(StalkRsrvElms_brch(ielmn,NB,NZ) &
        +SETN*StalkRsrvElms_brch(ielmc,NB,NZ))
      ZNPGP=StalkRsrvElms_brch(ielmp,NB,NZ)/(StalkRsrvElms_brch(ielmp,NB,NZ) &
        +SETP*StalkRsrvElms_brch(ielmc,NB,NZ))
      ZPGRN=ZPGRM+ZPGRD*AZMAX1(AMIN1(1.0_r8,ZNPGN))
      ZPGRP=ZPGRM+ZPGRD*AZMAX1(AMIN1(1.0_r8,ZNPGP))
      ChemElmRsrv2Grain(ielmn)=AMIN1(MaxChemElmRsrv2Grain*CNGR(NZ) &
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
      StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+GrowthGrain(NE) &
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
    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
      IF(HourFailGrainFill_brch(NB,NZ).GT.Hours4PhyslMature+Hours4SenesAftMature(iPlantPhenolType_pft(NZ)))THEN
        Hours4LeafOff_brch(NB,NZ)=HourReq4LeafOff_brch(NB,NZ)+0.5_r8
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine GrainFilling
!------------------------------------------------------------------------------------------

  subroutine PhenologyReset(I,NB,NZ)

  implicit none
  integer, intent(in) :: I,nb,nz
  integer :: K,M,NE
  real(r8) :: FSNR1
! begin_execution
  associate(                        &
    ifoliar                                =>  pltpar%ifoliar   , &
    istalk                                 =>  pltpar%istalk    , &
    iroot                                  =>  pltpar%iroot     , &
    inonfoliar                             =>  pltpar%inonfoliar  , &
    icwood                                 =>  pltpar%icwood    , &
    k_fine_litr                            =>  pltpar%k_fine_litr, &
    k_woody_litr                           =>  pltpar%k_woody_litr, &
    EHVST                                  =>  plt_distb%EHVST     , &
    iYearPlantHarvest_pft                  =>  plt_distb%iYearPlantHarvest_pft      , &
    THIN_pft                               =>  plt_distb%THIN_pft      , &
    HVST                                   =>  plt_distb%HVST      , &
    iHarvstType_pft                        =>  plt_distb%iHarvstType_pft     , &
    jHarvst_pft                            =>  plt_distb%jHarvst_pft     , &
    iDayPlanting_pft                       =>  plt_distb%iDayPlanting_pft     , &
    iDayPlantHarvest_pft                   =>  plt_distb%iDayPlantHarvest_pft     , &
    iYearPlanting_pft                      =>  plt_distb%iYearPlanting_pft      , &
    GrainStrutElms_brch                     =>  plt_biom%GrainStrutElms_brch     , &
    LeafElmntNode_brch                     =>  plt_biom%LeafElmntNode_brch      , &
    PopuRootMycoC_pvr                      =>  plt_biom% PopuRootMycoC_pvr     , &
    SenecStalkStrutElms_brch                =>  plt_biom%SenecStalkStrutElms_brch    , &
    StalkStrutElms_brch                     =>  plt_biom%StalkStrutElms_brch    , &
    RootMycoNonstElms_rpvr                   =>  plt_biom%RootMycoNonstElms_rpvr     , &
    PetoleStrutElms_brch                     =>  plt_biom%PetoleStrutElms_brch   , &
    EarStrutElms_brch                       =>  plt_biom%EarStrutElms_brch    , &
    LeafProteinCNode_brch                  =>  plt_biom%LeafProteinCNode_brch       , &
    PetioleProteinCNode_brch               =>  plt_biom%PetioleProteinCNode_brch     , &
    HuskStrutElms_brch                      =>  plt_biom%HuskStrutElms_brch    , &
    PetioleElmntNode_brch                  =>  plt_biom%PetioleElmntNode_brch     , &
    LeafStrutElms_brch                      =>  plt_biom%LeafStrutElms_brch    , &
    NonStrutElms_pft                      =>  plt_biom%NonStrutElms_pft      , &
    InternodeStrutElms_brch                  =>  plt_biom%InternodeStrutElms_brch     , &
    FWODLE                                 =>  plt_allom%FWODLE    , &
    FWODBE                                 =>  plt_allom%FWODBE    , &
    GrainSeedBiomCMean_brch                =>  plt_allom%GrainSeedBiomCMean_brch      , &
    HourReq4LeafOff_brch                   =>  plt_pheno%HourReq4LeafOff_brch     , &
    iPlantPhenolType_pft                   =>  plt_pheno%iPlantPhenolType_pft    , &
    iPlantRootProfile_pft                  =>  plt_pheno%iPlantRootProfile_pft    , &
    Hours4LeafOff_brch                     =>  plt_pheno%Hours4LeafOff_brch      , &
    HourReq4LeafOut_brch                   =>  plt_pheno%HourReq4LeafOut_brch     , &
    iPlantPhenolPattern_pft                =>  plt_pheno%iPlantPhenolPattern_pft    , &
    doPlantLeaveOff_brch                   =>  plt_pheno%doPlantLeaveOff_brch     , &
    doInitPlant_pft                        =>  plt_pheno%doInitPlant_pft     , &
    iPlantCalendar_brch                    =>  plt_pheno%iPlantCalendar_brch    , &
    doPlantLeafOut_brch                    =>  plt_pheno%doPlantLeafOut_brch     , &
    iPlantTurnoverPattern_pft              =>  plt_pheno%iPlantTurnoverPattern_pft    , &
    TotalNodeNumNormByMatgrp_brch          =>  plt_pheno%TotalNodeNumNormByMatgrp_brch    , &
    TotReproNodeNumNormByMatrgrp_brch      =>  plt_pheno%TotReproNodeNumNormByMatrgrp_brch    , &
    Hours4Leafout_brch                     =>  plt_pheno%Hours4Leafout_brch      , &
    MatureGroup_brch                       =>  plt_pheno%MatureGroup_brch    , &
    LeafNumberAtFloralInit_brch            =>  plt_pheno%LeafNumberAtFloralInit_brch    , &
    HourFailGrainFill_brch                 =>  plt_pheno%HourFailGrainFill_brch      , &
    Prep4Literfall_brch                    =>  plt_pheno%Prep4Literfall_brch     , &
    doInitLeafOut_brch                     =>  plt_pheno%doInitLeafOut_brch     , &
    MatureGroup_pft                        =>  plt_pheno%MatureGroup_pft   , &
    HoursCanopyPSITooLow                   =>  plt_pheno%HoursCanopyPSITooLow     , &
    Hours4LiterfalAftMature_brch           =>  plt_pheno%Hours4LiterfalAftMature_brch     , &
    KHiestGroLeafNode_brch                 =>  plt_pheno%KHiestGroLeafNode_brch    , &
    CFOPE                                  =>  plt_soilchem%CFOPE  , &
    iYearCurrent                           =>  plt_site%iYearCurrent       , &
    LitfalStrutElms_pvr                      =>  plt_bgcr%LitfalStrutElms_pvr       , &
    MY                                     =>  plt_morph%MY        , &
    NumOfLeaves_brch                       =>  plt_morph%NumOfLeaves_brch     , &
    XTLI                                   =>  plt_morph%XTLI      , &
    PetioleLengthNode_brch                 =>  plt_morph%PetioleLengthNode_brch     , &
    KLeafNumber_brch                       =>  plt_morph%KLeafNumber_brch    , &
    InternodeHeightDying_brch              =>  plt_morph%InternodeHeightDying_brch    , &
    InternodeHeightLive_brch               =>  plt_morph%InternodeHeightLive_brch    , &
    PotentialSeedSites_brch                =>  plt_morph%PotentialSeedSites_brch     , &
    NodeNumberAtAnthesis_brch              =>  plt_morph%NodeNumberAtAnthesis_brch     , &
    MainBranchNum_pft                      =>  plt_morph%MainBranchNum_pft       , &
    BranchNumber_brch                      =>  plt_morph%BranchNumber_brch      , &
    LeafAreaLive_brch                      =>  plt_morph%LeafAreaLive_brch     , &
    LeafAreaNode_brch                      =>  plt_morph%LeafAreaNode_brch     , &
    ShootNodeNumber_brch                   =>  plt_morph%ShootNodeNumber_brch     , &
    NodeNum2InitFloral_brch                =>  plt_morph%NodeNum2InitFloral_brch     , &
    SeedNumSet_brch                     =>  plt_morph%SeedNumSet_brch       &
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
  IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual.OR. &
    (iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).GT.1))THEN
!    write(192,*)'emerge plant',NB,doPlantLeafOut_brch(NB,NZ),Hours4Leafout_brch(NB,NZ),HourReq4LeafOut_brch(NB,NZ)
    IF((doPlantLeafOut_brch(NB,NZ).EQ.iEnable .AND. Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ)) &
      .OR.(doPlantLeaveOff_brch(NB,NZ).EQ.iEnable.AND.Hours4LeafOff_brch(NB,NZ).GE.HourReq4LeafOff_brch(NB,NZ)))THEN
      
!
  !    SPRING PHENOLOGY RESET
  !
  !    GROUP,MatureGroup_pft=node number required for floral initiation
  !    NodeNum2InitFloral_brch=node number at floral initiation
  !    NodeNumberAtAnthesis_brch=node number at flowering
  !    VSTGX=leaf number on date of floral initiation
  !    TotalNodeNumNormByMatgrp_brch=total change in vegve node number normalized for maturity group
  !    TotReproNodeNumNormByMatrgrp_brch=total change in reprve node number normalized for maturity group
  !    iPlantCalendar_brch(ipltcal_Emerge,=emergence date
!
      IF((doPlantLeafOut_brch(NB,NZ).EQ.iEnable.AND.iPlantPhenolPattern_pft(NZ).NE.iplt_annual) &
        .AND.(Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ)))THEN
        IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
          MatureGroup_brch(NB,NZ)=AZMAX1(MatureGroup_pft(NZ)-BranchNumber_brch(NB,NZ))
        ELSE
          MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
        ENDIF
        NodeNum2InitFloral_brch(NB,NZ)=ShootNodeNumber_brch(NB,NZ)
        NodeNumberAtAnthesis_brch(NB,NZ)=0._r8
        LeafNumberAtFloralInit_brch(NB,NZ)=0._r8
        TotalNodeNumNormByMatgrp_brch(NB,NZ)=0._r8
        TotReproNodeNumNormByMatrgrp_brch(NB,NZ)=0._r8
        iPlantCalendar_brch(ipltcal_Emerge,NB,NZ)=I
        D2005: DO M=2,NumGrowthStages
          iPlantCalendar_brch(M,NB,NZ)=0
        ENDDO D2005
        IF(NB.EQ.MainBranchNum_pft(NZ))THEN
          HoursCanopyPSITooLow(NZ)=0._r8
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
    !     LeafAreaLive_brch,LeafAreaNode_brch=branch,node leaf area
    !     InternodeStrutElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
    !     InternodeHeightDying_brch,InternodeHeightLive_brch=stalk height,stalk internode length
    !     SeedNumSet_brch=seed set number
    !     PotentialSeedSites_brch=potential number of seed set sites
    !     GrainSeedBiomCMean_brch=individual seed size
!
        IF(doPlantLeafOut_brch(NB,NZ).EQ.iEnable.AND.iPlantPhenolPattern_pft(NZ).NE.iplt_annual &
          .AND.Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
          IF(iPlantTurnoverPattern_pft(NZ).EQ.0)THEN
            ShootNodeNumber_brch(NB,NZ)=XTLI(NZ)
            NumOfLeaves_brch(NB,NZ)=0._r8
            KLeafNumber_brch(NB,NZ)=1
            KHiestGroLeafNode_brch(NB,NZ)=1
            HourFailGrainFill_brch(NB,NZ)=0._r8
            D5330: DO M=1,jsken
              DO NE=1,NumPlantChemElms
                LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
                  +CFOPE(NE,icwood,M,NZ)*LeafStrutElms_brch(NE,NB,NZ)*FWODLE(NE,k_woody_litr) &
                  +CFOPE(NE,icwood,M,NZ)*PetoleStrutElms_brch(NE,NB,NZ)*FWODBE(NE,k_woody_litr)

                LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
                  +CFOPE(NE,ifoliar,M,NZ)*LeafStrutElms_brch(NE,NB,NZ)*FWODLE(NE,k_fine_litr) &
                  +CFOPE(NE,inonfoliar,M,NZ)*PetoleStrutElms_brch(NE,NB,NZ)*FWODBE(NE,k_fine_litr)
              ENDDO  
            ENDDO D5330
            LeafAreaLive_brch(NB,NZ)=0._r8
            LeafStrutElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
            PetoleStrutElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
            D5335: DO K=0,MaxNodesPerBranch1
              LeafAreaNode_brch(K,NB,NZ)=0._r8
              PetioleLengthNode_brch(K,NB,NZ)=0._r8
              LeafProteinCNode_brch(K,NB,NZ)=0._r8
              LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)=0._r8
              PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)=0._r8
              PetioleProteinCNode_brch(K,NB,NZ)=0._r8
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
              LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
                +CFOPE(NE,inonfoliar,M,NZ) &
                *(HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ)+GrainStrutElms_brch(NE,NB,NZ))
            ENDDO
          ENDDO D6245
            
          HuskStrutElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
          EarStrutElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
          GrainStrutElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
          
          PotentialSeedSites_brch(NB,NZ)=0._r8
          SeedNumSet_brch(NB,NZ)=0._r8
          GrainSeedBiomCMean_brch(NB,NZ)=0._r8
          IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
            D6345: DO M=1,jsken
              DO NE=1,NumPlantChemElms
                LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)&
                  +CFOPE(NE,istalk,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)
              ENDDO
            ENDDO D6345
            StalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
            SenecStalkStrutElms_brch(1:NumPlantChemElms,NB,NZ)=0._r8
            DO K=0,MaxNodesPerBranch1
              DO NE=1,NumPlantChemElms
                InternodeStrutElms_brch(NE,K,NB,NZ)=0._r8
              ENDDO
            ENDDO
            D6340: DO K=0,MaxNodesPerBranch1
              InternodeHeightLive_brch(K,NB,NZ)=0._r8
              InternodeHeightDying_brch(K,NB,NZ)=0._r8
            ENDDO D6340
          ENDIF
        ENDIF
      ENDIF
!
  !   SPRING OR FALL FLAG RESET
  !
      IF(doPlantLeafOut_brch(NB,NZ).EQ.iEnable.AND.Hours4Leafout_brch(NB,NZ) &
        .GE.HourReq4LeafOut_brch(NB,NZ))THEN
        doPlantLeafOut_brch(NB,NZ)=iDisable
        doPlantLeaveOff_brch(NB,NZ)=iEnable
        Prep4Literfall_brch(NB,NZ)=ifalse
        Hours4LiterfalAftMature_brch(NB,NZ)=0
      ELSE
        !doing leave off
        doPlantLeafOut_brch(NB,NZ)=iEnable
        doPlantLeaveOff_brch(NB,NZ)=iDisable
        Prep4Literfall_brch(NB,NZ)=itrue
        Hours4LiterfalAftMature_brch(NB,NZ)=0
        doInitLeafOut_brch(NB,NZ)=iEnable
      ENDIF
    ENDIF
  ENDIF
!
!   REPRODUCTIVE MATERIAL BECOMES LitrFall AT END OF SEASON
!
  IF(Prep4Literfall_brch(NB,NZ).EQ.itrue)THEN
    Hours4LiterfalAftMature_brch(NB,NZ)=Hours4LiterfalAftMature_brch(NB,NZ)+1
    IF(Hours4LiterfalAftMature_brch(NB,NZ).EQ.HoursReq4LiterfalAftMature)THEN
      Prep4Literfall_brch(NB,NZ)=ifalse
      Hours4LiterfalAftMature_brch(NB,NZ)=0
    ENDIF
    FSNR1=1.0_r8-FSNR

    D6330: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+ &
          FSNR*CFOPE(NE,inonfoliar,M,NZ)*(HuskStrutElms_brch(NE,NB,NZ)+EarStrutElms_brch(NE,NB,NZ))
        IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
          NonStrutElms_pft(NE,NZ)=NonStrutElms_pft(NE,NZ)+ &
            FSNR*CFOPE(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
          print*,'plantbranch2375NonStrutElms_pft(NE,NZ)',NZ,NE,NonStrutElms_pft(NE,NZ)  
        ELSE
          LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+ &
            FSNR*CFOPE(NE,inonfoliar,M,NZ)*GrainStrutElms_brch(NE,NB,NZ)
        ENDIF
      ENDDO
    ENDDO D6330

    DO NE=1,NumPlantChemElms  
      HuskStrutElms_brch(NE,NB,NZ)=FSNR1*HuskStrutElms_brch(NE,NB,NZ)
      EarStrutElms_brch(NE,NB,NZ)=FSNR1*EarStrutElms_brch(NE,NB,NZ)
      GrainStrutElms_brch(NE,NB,NZ)=FSNR1*GrainStrutElms_brch(NE,NB,NZ)
    ENDDO
    PotentialSeedSites_brch(NB,NZ)=FSNR1*PotentialSeedSites_brch(NB,NZ)
    SeedNumSet_brch(NB,NZ)=FSNR1*SeedNumSet_brch(NB,NZ)
    GrainSeedBiomCMean_brch(NB,NZ)=FSNR1*GrainSeedBiomCMean_brch(NB,NZ)
!
!     STALKS BECOME LitrFall IN GRASSES AT END OF SEASON
!
    IF((iPlantTurnoverPattern_pft(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ)))) &
      .AND.iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN

      D6335: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)&
            +FSNR*CFOPE(NE,istalk,M,NZ)*StalkStrutElms_brch(NE,NB,NZ)
        ENDDO
      ENDDO D6335
      DO NE=1,NumPlantChemElms  
        StalkStrutElms_brch(NE,NB,NZ)=FSNR1*StalkStrutElms_brch(NE,NB,NZ)
        SenecStalkStrutElms_brch(NE,NB,NZ)=FSNR1*SenecStalkStrutElms_brch(NE,NB,NZ)
      ENDDO
      DO K=0,MaxNodesPerBranch1
        DO NE=1,NumPlantChemElms
          InternodeStrutElms_brch(NE,K,NB,NZ)=FSNR1*InternodeStrutElms_brch(NE,K,NB,NZ)
        ENDDO
      ENDDO
      D2010: DO K=0,MaxNodesPerBranch1
    !     InternodeHeightLive_brch(K,NB,NZ)=FSNR1*InternodeHeightLive_brch(K,NB,NZ)
        InternodeHeightDying_brch(K,NB,NZ)=FSNR1*InternodeHeightDying_brch(K,NB,NZ)
      ENDDO D2010
    ENDIF

!
!     SELF-SEEDING ANNUALS IF COLD OR DROUGHT DECIDUOUS
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
!     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from PFT
!     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
!           leaf,non-foliar,woody, standing dead removed from ecosystem
!     iDayPlanting_pft,iYearPlanting_pft=day,year of planting
!     doInitPlant_pft=PFT initialization flag:0=no,1=yes
!
!     IF(J.EQ.INT(SolarNoonHour_col))THEN

    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      !deciduous annual plant
      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenolType_pft(NZ).NE.0)THEN
        iDayPlantHarvest_pft(NZ)=I
        iYearPlantHarvest_pft(NZ)=iYearCurrent
        iHarvstType_pft(NZ)=1
        jHarvst_pft(NZ)=jharvtyp_tmareseed
        HVST(NZ)=0._r8
        THIN_pft(NZ)=0._r8
        EHVST(1,iplthvst_leaf,NZ)=1.0_r8
        EHVST(1,iplthvst_finenonleaf,NZ)=1.0_r8
        EHVST(1,iplthvst_woody,NZ)=1.0_r8
        EHVST(1,iplthvst_stdead,NZ)=1.0_r8
        EHVST(2,iplthvst_leaf,NZ)=0._r8
        EHVST(2,iplthvst_finenonleaf,NZ)=1.0_r8
        EHVST(2,iplthvst_woody,NZ)=0._r8
        EHVST(2,iplthvst_stdead,NZ)=0._r8
        iDayPlanting_pft(NZ)=-1E+06
        iYearPlanting_pft(NZ)=-1E+06
        doInitPlant_pft(NZ)=itrue
      ENDIF
    ENDIF
!     ENDIF
  ENDIF
  end associate
  end subroutine PhenologyReset
!------------------------------------------------------------------------------------------
  subroutine BranchElmntTransfer(I,J,NB,NZ,BegRemoblize,WFNG,WFNSG)
  implicit none
  integer, intent(in) :: I,J,NB,NZ,BegRemoblize
  real(r8), intent(in) :: WFNG,WFNSG
  integer :: L,NE
  real(r8) :: NonstGradt
  real(r8) :: XFRPX,XFRCX,XFRNX
  real(r8) :: ATRPPD
  REAL(R8) :: cpoolt
  real(r8) :: CH2OH
  real(r8) :: CWTRSV
  REAL(R8) :: CWTRSN,CWTRSP
  real(r8) :: CNR,CPR
  real(r8) :: CNL,CPL
  real(r8) :: CPOOLD
  real(r8) :: DATRP
  real(r8) :: FXFC,FXFN
  real(r8) :: FracCanopyCinStalk
  real(r8) :: GFNX
  real(r8) :: PPDX
  real(r8) :: TotalRootNonstElms(1:NumPlantChemElms)
  real(r8) :: ElmXferStore2Shoot(1:NumPlantChemElms)
  real(r8) :: TotPopuPlantRootC,WFNSP
  real(r8) :: ShootBiomC_brch,WTRTRX
  real(r8) :: WTPLTX,WVSTBX
  real(r8) :: WTRTTX,WTRSBX
  real(r8) :: WTRVCX
  real(r8) :: XFRE(1:NumPlantChemElms)
  real(r8) :: NonstElm2RootMyco(NumPlantChemElms)
  logical :: PlantingChk,RemobChk,LeafOutChk,PlantChk
  ! begin_execution
  associate(                          &
    iDayPlanting_pft                =>  plt_distb%iDayPlanting_pft  , &
    iYearPlanting_pft               =>  plt_distb%iYearPlanting_pft   , &
    iYearCurrent                    =>  plt_site%iYearCurrent    , &
    ZEROS2                          =>  plt_site%ZEROS2   , &
    DayLenthCurrent                 =>  plt_site%DayLenthCurrent    , &
    k_woody_litr                    =>  pltpar%k_woody_litr, &    
    NU                              =>  plt_site%NU      , &
    FracHour4LeafoffRemob           =>  plt_allom%FracHour4LeafoffRemob   , &
    FWODRE                          =>  plt_allom%FWODRE , &
    RootElms_pft                    =>  plt_biom%RootElms_pft   , &
     PopuRootMycoC_pvr              =>  plt_biom% PopuRootMycoC_pvr  , &
    RootMycoActiveBiomC_pvr         =>  plt_biom%RootMycoActiveBiomC_pvr  , &
    CanopyStalkC_pft                =>  plt_biom%CanopyStalkC_pft   , &
    CanopyNonstElms_brch            =>  plt_biom%CanopyNonstElms_brch  , &
     RootMycoNonstElms_rpvr         =>  plt_biom%RootMycoNonstElms_rpvr  , &
    NonStrutElms_pft                =>  plt_biom%NonStrutElms_pft   , &
    LeafPetoNonstElmConc_brch       =>  plt_biom%LeafPetoNonstElmConc_brch  , &
    LeafPetolBiomassC_brch          =>  plt_biom%LeafPetolBiomassC_brch   , &
    ZEROP                           =>  plt_biom%ZEROP   , &
    StalkRsrvElms_brch              =>  plt_biom%StalkRsrvElms_brch , &
    StalkBiomassC_brch              =>  plt_biom%StalkBiomassC_brch  , &
    VLSoilPoreMicP                  =>  plt_soilchem%VLSoilPoreMicP, &
    iPlantCalendar_brch             =>  plt_pheno%iPlantCalendar_brch , &
    fTgrowCanP                      =>  plt_pheno%fTgrowCanP   , &
    HourReq4LeafOff_brch            =>  plt_pheno%HourReq4LeafOff_brch  , &
    Hours4Leafout_brch              =>  plt_pheno%Hours4Leafout_brch   , &
    HourReq4LeafOut_brch            =>  plt_pheno%HourReq4LeafOut_brch  , &
    doInitPlant_pft                 =>  plt_pheno%doInitPlant_pft  , &
    PhotoPeriodSens_pft             =>  plt_pheno%PhotoPeriodSens_pft  , &
    iPlantPhenolPattern_pft         =>  plt_pheno%iPlantPhenolPattern_pft , &
    iPlantPhotoperiodType_pft       =>  plt_pheno%iPlantPhotoperiodType_pft , &
    CriticPhotoPeriod_pft           =>  plt_pheno%CriticPhotoPeriod_pft   , &
    Hours4LeafOff_brch              =>  plt_pheno%Hours4LeafOff_brch   , &
    iPlantTurnoverPattern_pft       =>  plt_pheno%iPlantTurnoverPattern_pft , &
    iPlantRootProfile_pft           =>  plt_pheno%iPlantRootProfile_pft , &
    iPlantPhenolType_pft            =>  plt_pheno%iPlantPhenolType_pft , &
    doInitLeafOut_brch              =>  plt_pheno%doInitLeafOut_brch  , &
    Hours2LeafOut_brch              =>  plt_pheno%Hours2LeafOut_brch  , &
    NGTopRootLayer_pft              =>  plt_morph%NGTopRootLayer_pft   , &
    MainBranchNum_pft               =>  plt_morph%MainBranchNum_pft    , &
    MaxSoiL4Root                    =>  plt_morph%MaxSoiL4Root      &
  )
!   TRANSFER C,N,P FROM SEASONAL STORAGE TO SHOOT AND ROOT
!   NON-STRUCTURAL C DURING SEED GERMINATION OR LEAFOUT
!
  PlantingChk=I.GE.iDayPlanting_pft(NZ).AND.iYearCurrent.EQ.iYearPlanting_pft(NZ)
  RemobChk=Hours4LeafOff_brch(NB,NZ).LT.FracHour4LeafoffRemob(iPlantPhenolType_pft(NZ))*HourReq4LeafOff_brch(NB,NZ)  
  LeafOutChk=Hours4Leafout_brch(MainBranchNum_pft(NZ),NZ).GE.HourReq4LeafOut_brch(NB,NZ)

  IF((iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.doInitPlant_pft(NZ).EQ.ifalse) &
    .OR.(PlantingChk.AND.RemobChk) .OR. (LeafOutChk.AND.RemobChk))THEN
    TotPopuPlantRootC=0._r8
    TotalRootNonstElms(ielmc)=0._r8
    D4: DO L=NU,MaxSoiL4Root(NZ)
      TotPopuPlantRootC=TotPopuPlantRootC+AZMAX1(PopuRootMycoC_pvr(ipltroot,L,NZ))
      TotalRootNonstElms(ielmc)=TotalRootNonstElms(ielmc)+&
        AZMAX1(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))
    ENDDO D4
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
  ! WFNSG=expansion,extension function of canopy water potential
  ! fTgrowCanP=temperature function for canopy growth
  ! HourReq2InitSStor4LeafOut=number of hours required to initiate remobilization of storage C for leafout
  ! main branch leaf out
    IF(NB.EQ.MainBranchNum_pft(NZ))THEN
      IF(iPlantPhotoperiodType_pft(NZ).EQ.iphotop_long &
        .AND.(iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldecid &
        .OR.iPlantPhenolType_pft(NZ).EQ.iphenotyp_coldroutdecid))THEN
        PPDX=AZMAX1(CriticPhotoPeriod_pft(NZ)-PhotoPeriodSens_pft(NZ)-DayLenthCurrent)
        ATRPPD=EXP(-0.0_r8*PPDX)
      ELSE
        ATRPPD=1.0_r8
      ENDIF
      IF(.not.is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
        WFNSP=WFNSG
      ELSE
        WFNSP=1.0_r8
      ENDIF
      DATRP=ATRPPD*fTgrowCanP(NZ)*WFNSP
      Hours2LeafOut_brch(NB,NZ)=Hours2LeafOut_brch(NB,NZ)+DATRP
      PlantChk=iPlantPhenolPattern_pft(NZ).EQ.iplt_annual .AND.iPlantPhenolType_pft(NZ).EQ.iphenotyp_evgreen
      IF(Hours2LeafOut_brch(NB,NZ).LE.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ)).OR.(plantChk))THEN
!        write(101,*)'NonStrutElms_pft LEAFOUT',NonStrutElms_pft(ielmc,NZ),NZ
        IF(NonStrutElms_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
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
          GFNX=GVMX(iPlantPhenolPattern_pft(NZ))*DATRP
          CH2OH=AZMAX1(GFNX*NonStrutElms_pft(ielmc,NZ))
          NonStrutElms_pft(ielmc,NZ)=NonStrutElms_pft(ielmc,NZ)-CH2OH
          CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)+CH2OH &
            *FXSH(iPlantPhenolPattern_pft(NZ))
          write(*,*)'2618NonStrutElms_pft(ielmc,NZ)',NZ,NonStrutElms_pft(ielmc,NZ)  
          IF(ABS(NonStrutElms_pft(ielmc,NZ))>1.E10)STOP
!          write(101,*)'LEAFOUTCanopyNonstElms_brch',CanopyNonstElms_brch(ielmc,NB,NZ),NB,NZ
!         if root condition met
          IF(TotPopuPlantRootC.GT.ZEROP(NZ).AND.TotalRootNonstElms(ielmc).GT.ZEROP(NZ))THEN
            D50: DO L=NU,MaxSoiL4Root(NZ)
              FXFC=AZMAX1(PopuRootMycoC_pvr(ipltroot,L,NZ))/TotPopuPlantRootC
              RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)&
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
      IF(NonStrutElms_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
        IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
          CPOOLT=AZMAX1(NonStrutElms_pft(ielmc,NZ)+CanopyNonstElms_brch(ielmc,NB,NZ))          
          DO NE=2,NumPlantChemElms
            NonstGradt=(NonStrutElms_pft(NE,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ)- &
              CanopyNonstElms_brch(NE,NB,NZ)*NonStrutElms_pft(ielmc,NZ))/CPOOLT            
            ElmXferStore2Shoot(NE)=AZMAX1(FRSV(iPlantTurnoverPattern_pft(NZ))*NonstGradt)
          ENDDO

        ELSE
          DO NE=2,NumPlantChemElms
            ElmXferStore2Shoot(NE)=AZMAX1(FXSH(iPlantPhenolPattern_pft(NZ))*CH2OH*NonStrutElms_pft(NE,NZ) &
              /NonStrutElms_pft(ielmc,NZ))
          ENDDO  
        ENDIF
      ELSE
        !when there is no seasonal storage C
        DO NE=2,NumPlantChemElms
          ElmXferStore2Shoot(NE)=AZMAX1(FXSH(iPlantPhenolPattern_pft(NZ))*NonStrutElms_pft(NE,NZ))
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
        D3: DO L=NU,MaxSoiL4Root(NZ)
          TotalRootNonstElms(NE)=TotalRootNonstElms(NE)+ &
            AZMAX1(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ))
        ENDDO D3
      ENDDO

      IF(NonStrutElms_pft(ielmc,NZ).GT.ZEROP(NZ))THEN
        IF(iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
          CPOOLT=AMAX1(ZEROP(NZ),NonStrutElms_pft(ielmc,NZ)+TotalRootNonstElms(ielmc))
          DO NE=2,NumPlantChemElms
            NonstGradt=(NonStrutElms_pft(NE,NZ)*TotalRootNonstElms(ielmc)-&
              TotalRootNonstElms(NE)*NonStrutElms_pft(ielmc,NZ))/CPOOLT            
            NonstElm2RootMyco(NE)=AZMAX1(FRSV(iPlantTurnoverPattern_pft(NZ))*NonstGradt)
          ENDDO

        ELSE
          DO NE=2,NumPlantChemElms
            NonstElm2RootMyco(NE)=AZMAX1(FXRT(iPlantPhenolPattern_pft(NZ))*CH2OH*&
              NonStrutElms_pft(NE,NZ)/NonStrutElms_pft(ielmc,NZ))
          ENDDO    
        ENDIF
      ELSE
        DO NE=2,NumPlantChemElms
          NonstElm2RootMyco(NE)=AZMAX1(FXRT(iPlantPhenolPattern_pft(NZ))*NonStrutElms_pft(NE,NZ))
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
        NonStrutElms_pft(NE,NZ)=NonStrutElms_pft(NE,NZ)-ElmXferStore2Shoot(NE)-NonstElm2RootMyco(NE)
        CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+ElmXferStore2Shoot(NE)
      ENDDO
      
      IF(TotPopuPlantRootC.GT.ZEROP(NZ).AND.TotalRootNonstElms(ielmc).GT.ZEROP(NZ))THEN
        D51: DO L=NU,MaxSoiL4Root(NZ)
          FXFN=AZMAX1(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))/TotalRootNonstElms(ielmc)
          DO NE=2,NumPlantChemElms
            RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)&
              +FXFN*NonstElm2RootMyco(NE)
            if(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)<0._r8)stop            
          ENDDO  
        ENDDO D51
      ELSE
        DO NE=2,NumPlantChemElms
          RootMycoNonstElms_rpvr(NE,ipltroot,NGTopRootLayer_pft(NZ),NZ)=&
            RootMycoNonstElms_rpvr(NE,ipltroot,NGTopRootLayer_pft(NZ),NZ)+NonstElm2RootMyco(NE)
        ENDDO            
      ENDIF
    ENDIF
  !
  ! REDISTRIBUTE TRANFERRED C FROM MAIN STEM TO OTHER BRANCHES
  !
  ! ATRP=hourly leafout counter
  ! fTgrowCanP=temperature function for canopy growth
  ! HourReq2InitSStor4LeafOut=number of hours required for remobilization of storage C during leafout
  ! WFNG=growth function of canopy water potential
  ! CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
  ! XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
  ! 
!    write(101,*)NB,MainBranchNum_pft(NZ),Hours2LeafOut_brch(NB,NZ).LE.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ))
    IF(NB.NE.MainBranchNum_pft(NZ) .AND. &
      Hours2LeafOut_brch(NB,NZ).LE.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ)))THEN
      Hours2LeafOut_brch(NB,NZ)=Hours2LeafOut_brch(NB,NZ)+fTgrowCanP(NZ)*WFNG
!      write(101,*)'Hours2LeafOut_brch',Hours2LeafOut_brch(NB,NZ),NB,NZ
      DO NE=1,NumPlantChemElms
        XFRE(NE)=AZMAX1(0.05_r8*fTgrowCanP(NZ) &
          *(0.5_r8*(CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)+CanopyNonstElms_brch(NE,NB,NZ)) &
          -CanopyNonstElms_brch(NE,NB,NZ)))
        CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE)
        CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)=CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)-XFRE(NE)
        if(CanopyNonstElms_brch(NE,NB,NZ)<0._r8)then
          write(*,*)'2755CanopyNonstElms_brch(NE,NB,NZ)',NE,NB,CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE),XFRE(NE)          
          stop
        endif
      ENDDO
    ENDIF
  ENDIF

!
! TRANSFER LEAF AND STALK NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
! IN PERENNIALS AFTER GRAIN FILL IN DETERMINATES, AFTER AUTUMNIZ'N
! IN INDETERMINATES, OR AFTER SUSTAINED WATER STRESS
!
! StalkBiomassC_brch=stalk sapwood mass
! WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
! FXFB=rate constant for plant-storage nonstructural C,N,P exchange
! XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
! BegRemoblize=remobilization flag
! StalkBiomassC_brch=stalk sapwood mass, it holds reserve biomass
! WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
! XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
! WTRVC,WTRVN,WTRVP=storage C,N,P
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!
  IF(BegRemoblize.EQ.itrue.AND.iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
    IF(StalkBiomassC_brch(NB,NZ).GT.ZEROP(NZ).AND. StalkRsrvElms_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      CWTRSV=AZMAX1(StalkRsrvElms_brch(ielmc,NB,NZ)/StalkBiomassC_brch(NB,NZ))
      CWTRSN=AZMAX1(StalkRsrvElms_brch(ielmn,NB,NZ)/StalkBiomassC_brch(NB,NZ))
      CWTRSP=AZMAX1(StalkRsrvElms_brch(ielmp,NB,NZ)/StalkBiomassC_brch(NB,NZ))
      CNR=CWTRSV/(CWTRSV+CWTRSN/CNKI)
      CPR=CWTRSV/(CWTRSV+CWTRSP/CPKI)
    ELSE
      CNR=0._r8
      CPR=0._r8
    ENDIF
    XFRCX=FXFB(iPlantTurnoverPattern_pft(NZ))*AZMAX1(StalkRsrvElms_brch(ielmc,NB,NZ))
    XFRNX=FXFB(iPlantTurnoverPattern_pft(NZ))*AZMAX1(StalkRsrvElms_brch(ielmn,NB,NZ))*(1.0_r8+CNR)
    XFRPX=FXFB(iPlantTurnoverPattern_pft(NZ))*AZMAX1(StalkRsrvElms_brch(ielmp,NB,NZ))*(1.0_r8+CPR)
    XFRE(ielmc)=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
    XFRE(ielmn)=AMIN1(XFRNX,XFRE(ielmc)*CNMX,XFRPX*CNMX/CPMN*0.5_r8)
    XFRE(ielmp)=AMIN1(XFRPX,XFRE(ielmc)*CPMX,XFRNX*CPMX/CNMN*0.5_r8)
    DO NE=1,NumPlantChemElms
      StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)-XFRE(NE)
      NonStrutElms_pft(NE,NZ)=NonStrutElms_pft(NE,NZ)+XFRE(NE)
      if(NE==1)print*,'2809StalkRsrvElms_brch(NE,NB,NZ)',II+JJ/24.,NE,StalkRsrvElms_brch(NE,NB,NZ),XFRE(NE)
    ENDDO
    IF(LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      CNL=LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ) &
        +LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/CNKI)
      CPL=LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ) &
        +LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/CPKI)
    ELSE
      CNL=0._r8
      CPL=0._r8
    ENDIF
    XFRCX=FXFB(iPlantTurnoverPattern_pft(NZ))*AZMAX1(CanopyNonstElms_brch(ielmc,NB,NZ))
    XFRNX=FXFB(iPlantTurnoverPattern_pft(NZ))*AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))*(1.0_r8+CNL)
    XFRPX=FXFB(iPlantTurnoverPattern_pft(NZ))*AZMAX1(CanopyNonstElms_brch(ielmp,NB,NZ))*(1.0_r8+CPL)
    XFRE(ielmc)=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
    XFRE(ielmn)=AMIN1(XFRNX,XFRE(ielmc)*CNMX,XFRPX*CNMX/CPMN*0.5_r8)
    XFRE(ielmp)=AMIN1(XFRPX,XFRE(ielmc)*CPMX,XFRNX*CPMX/CNMN*0.5_r8)
    DO NE=1,NumPlantChemElms
      CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE)
      NonStrutElms_pft(NE,NZ)=NonStrutElms_pft(NE,NZ)+XFRE(NE)
      if(CanopyNonstElms_brch(NE,NB,NZ)<0._r8)then
        write(*,*)'2821CanopyNonstElms_brch(NE,NB,NZ)',NE,NB,CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE),XFRE(NE)
        stop
      endif
    ENDDO
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
!   StalkBiomassC_brch=stalk sapwood mass
!   CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!   FXFY,FXFZ=rate constant for plant-reserve nonstructural C,N,P exchange
!   XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
!   CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
  IF((iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.&
    iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).NE.0) &
    .OR.(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial.AND.&
    iPlantCalendar_brch(ipltcal_Jointing,NB,NZ).NE.0))THEN
    ShootBiomC_brch=LeafPetolBiomassC_brch(NB,NZ)+StalkBiomassC_brch(NB,NZ)
    CPOOLT=CanopyNonstElms_brch(ielmc,NB,NZ)+StalkRsrvElms_brch(ielmc,NB,NZ)
    IF(ShootBiomC_brch.GT.ZEROP(NZ))THEN
      CPOOLD=(CanopyNonstElms_brch(ielmc,NB,NZ)*StalkBiomassC_brch(NB,NZ) &
        -StalkRsrvElms_brch(ielmc,NB,NZ)*LeafPetolBiomassC_brch(NB,NZ))/ShootBiomC_brch
      XFRE(ielmc)=FXFY(iPlantPhenolPattern_pft(NZ))*CPOOLD
      CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)-XFRE(ielmc)
      StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
    ENDIF
    IF(CPOOLT.GT.ZEROP(NZ))THEN
      DO NE=2,NumPlantChemElms
        NonstGradt=(CanopyNonstElms_brch(NE,NB,NZ)*StalkRsrvElms_brch(ielmc,NB,NZ) &
          -StalkRsrvElms_brch(NE,NB,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ))/CPOOLT
        XFRE(NE)=FXFZ(iPlantPhenolPattern_pft(NZ))*NonstGradt
        CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE)
        StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+XFRE(NE)
        if(CanopyNonstElms_brch(NE,NB,NZ)<0._r8)then
          write(*,*)'2869CanopyNonstElms_brch(NE,NB,NZ)',NE,NB,CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE),XFRE(NE)
          stop
        endif
      ENDDO
    ENDIF

    IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual.AND.&
      iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).NE.0)THEN
      !stalk-root transfer
      D2050: DO L=NU,MaxSoiL4Root(NZ)
        IF(VLSoilPoreMicP(L).GT.ZEROS2)THEN
          WTRTRX=AMAX1(ZEROP(NZ),RootMycoActiveBiomC_pvr(ipltroot,L,NZ)*FWODRE(ielmc,k_woody_litr))
          WTPLTX=WTRTRX+StalkBiomassC_brch(NB,NZ)
          IF(WTPLTX.GT.ZEROP(NZ))THEN
            CPOOLD=(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)*StalkBiomassC_brch(NB,NZ) &
              -StalkRsrvElms_brch(ielmc,NB,NZ)*WTRTRX)/WTPLTX
            XFRE(ielmc)=AZMAX1(FXFY(iPlantPhenolPattern_pft(NZ))*CPOOLD)
            RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)-XFRE(ielmc)
            StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
            CPOOLT=RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)+StalkRsrvElms_brch(ielmc,NB,NZ)
            write(*,*)'2904StalkRsrvElms_brch(ielmc,NB,NZ)',II+JJ/24.,StalkRsrvElms_brch(ielmc,NB,NZ),XFRE(ielmc)
            IF(CPOOLT.GT.ZEROP(NZ))THEN
              DO NE=2,NumPlantChemElms
                NonstGradt=(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)*StalkRsrvElms_brch(ielmc,NB,NZ)-&
                StalkRsrvElms_brch(NE,NB,NZ)*RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))/CPOOLT
                XFRE(NE)=AZMAX1(FXFZ(iPlantPhenolPattern_pft(NZ))*NonstGradt)              
                RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)-XFRE(NE)
                StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+XFRE(NE)
                if(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)<0._r8)then
                  write(*,*)'2898RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)',NE,&
                    RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)+XFRE(NE),XFRE(NE)
                  stop
                endif
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDDO D2050
    ENDIF
  ENDIF
!
!   REPLENISH BRANCH NON-STRUCTURAL POOL FROM
!   SEASONAL STORAGE POOL
!
!   StalkBiomassC_brch,WVSTK=stalk,total stalk sapwood mass
!   WTRT=total root mass
!   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!   XFRX=maximum storage C content for remobiln from stalk,root reserves
!   XFRE(ielmc)=C transfer
!   Q: why are nitrogen and phosphorus not transferred?
  IF(StalkBiomassC_brch(NB,NZ).GT.ZEROP(NZ).AND.CanopyStalkC_pft(NZ).GT.ZEROP(NZ) &
    .AND.RootElms_pft(ielmc,NZ).GT.ZEROP(NZ) &
    .AND.StalkRsrvElms_brch(ielmc,NB,NZ).LE.XFRX*StalkBiomassC_brch(NB,NZ))THEN
    FracCanopyCinStalk=StalkBiomassC_brch(NB,NZ)/CanopyStalkC_pft(NZ)
    WVSTBX=StalkBiomassC_brch(NB,NZ)
    WTRTTX=RootElms_pft(ielmc,NZ)*FracCanopyCinStalk
    ShootBiomC_brch=WVSTBX+WTRTTX
    WTRSBX=AZMAX1(StalkRsrvElms_brch(ielmc,NB,NZ))
    WTRVCX=AZMAX1(NonStrutElms_pft(ielmc,NZ)*FracCanopyCinStalk)
    CPOOLD=(WTRVCX*WVSTBX-WTRSBX*WTRTTX)/ShootBiomC_brch
    XFRE(ielmc)=AZMAX1(XFRY*CPOOLD)
    StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
    NonStrutElms_pft(ielmc,NZ)=NonStrutElms_pft(ielmc,NZ)-XFRE(ielmc)
    write(*,*)'plantbranch2934NonStrutElms_pft(ielmc,NZ)',NZ,NonStrutElms_pft(ielmc,NZ)
    if(abs(NonStrutElms_pft(ielmc,NZ))>1.e15)stop
  ENDIF
  end associate
  end subroutine BranchElmntTransfer
!------------------------------------------------------------------------------------------

  subroutine ComputRAutoAfEmergence(NB,NZ,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,CO2F,&
    CH2O,TFN5,WFNG,WFNSG,ShootStructN,CanopyNonstElm4Gros,CNPG,RCO2C,RMNCS,RMxess_brch,&
    NonStructalC4Growth_brch,Rauto4Nassim_brch)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8), intent(out) :: CanopyNonstElm4Gros(NumPlantChemElms)
  real(r8), INTENT(OUT) :: CNPG
  real(r8), intent(out) :: RCO2C,RMNCS,RMxess_brch,NonStructalC4Growth_brch,Rauto4Nassim_brch
  real(r8), intent(in) :: DMSHD
  real(r8), intent(in) :: CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX
  real(r8), intent(in) :: ShootStructN,CO2F,CH2O,TFN5
  real(r8), intent(in) :: WFNG,WFNSG
  real(r8) :: ZPOOLB
  real(r8) :: PPOOLB
  real(r8) :: RCO2X,RCO2Y
  real(r8) :: RgroCO2_ltd  !Nutient limited growth respiration
  real(r8) :: Rauto_brch,RCO2CM
! begin_execution
  associate(                             &
    CO2NetFix_pft                      =>  plt_bgcr%CO2NetFix_pft    , &
    GrossResp_pft                      =>  plt_bgcr%GrossResp_pft   , &
    Eco_AutoR_col                      =>  plt_bgcr%Eco_AutoR_col    , &
    ECO_ER_col                         =>  plt_bgcr%ECO_ER_col    , &
    Eco_GPP_col                        =>  plt_bgcr%Eco_GPP_col    , &
    GrossCO2Fix_pft                    =>  plt_bgcr%GrossCO2Fix_pft   , &
    CanopyPlusNoduRespC_pft            =>  plt_bgcr%CanopyPlusNoduRespC_pft   , &
    CanopyNonstElms_brch                  =>  plt_biom%CanopyNonstElms_brch  , &
    LeafPetoNonstElmConc_brch          =>  plt_biom%LeafPetoNonstElmConc_brch  , &
    ZERO                               =>  plt_site%ZERO     , &
    iPlantRootProfile_pft              =>  plt_pheno%iPlantRootProfile_pft , &
    fTgrowCanP                         =>  plt_pheno%fTgrowCanP   , &
    iPlantPhenolType_pft               =>  plt_pheno%iPlantPhenolType_pft , &
    C4PhotosynDowreg_brch              =>  plt_photo%C4PhotosynDowreg_brch    &
  )
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
! CNKI,CPKI=nonstruct N,P inhibn constant on growth(g N,P g-1 C)
!
  IF(LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZERO)THEN
    CNPG=AMIN1(LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI),LeafPetoNonstElmConc_brch(ielmp,NB,NZ) &
      /(LeafPetoNonstElmConc_brch(ielmp,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI))
  ELSE
    CNPG=1.0_r8
  ENDIF
!
! RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
! NON-STRUCTURAL C:N:P
!
! RCO2C=respiration from non-structural C
! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
! CPOOL=non-structural C mass
! fTgrowCanP=temperature function for canopy growth
! WFNG=growth function of canopy water potential
! CNPG=N,P constraint on respiration
! C4PhotosynDowreg_brch=termination feedback inhibition on C3 CO2
!
  RCO2C=AZMAX1(VMXC*CanopyNonstElms_brch(ielmc,NB,NZ) &
    *fTgrowCanP(NZ))*CNPG*C4PhotosynDowreg_brch(NB,NZ)*WFNG
!
! MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
!
! RMNCS=maintenance respiration
! TFN5=temperature function for canopy maintenance respiration
! ShootStructN=shoot structural N mass
! iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
! iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
! WFNG=growth function of canopy water potential
!
  RMNCS=AZMAX1(RmSpecPlant*TFN5*ShootStructN)
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)).OR.&
    iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    RMNCS=RMNCS*WFNG
  ENDIF
!
! GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
! IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
! RCO2X=difference between non-structural C respn and mntc respn
! RCO2Y=growth respiration unlimited by N,P
! WFNSG=expansion,extension function of canopy water potential
! RMxess_brch=excess maintenance respiration, drives remobilization & senescence
!
  RCO2X=RCO2C-RMNCS
  RCO2Y=AZMAX1(RCO2X)*WFNSG
  RMxess_brch=AZMAX1(-RCO2X)
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
  IF(RCO2Y.GT.0.0_r8.AND.(CNSHX.GT.0.0_r8.OR.CNLFX.GT.0.0_r8))THEN
    ZPOOLB=AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))
    PPOOLB=AZMAX1(CanopyNonstElms_brch(ielmp,NB,NZ))
    RgroCO2_ltd=AMIN1(RCO2Y,ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG) &
      ,PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))
  ELSE
    RgroCO2_ltd=0._r8
  ENDIF

!
! TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
! AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
! ENTERED IN 'READQ'
!
! NonStructalC4Growth_brch=total non-structural C used in growth and growth respiration
! RgroCO2_ltd=growth respiration limited by N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! CanopyNonstElm4Gros(ielmn),CanopyNonstElm4Gros(ielmp)=nonstructural N,P used in growth
! ZPOOL,PPOOL=nonstructural N,P mass
! CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! Rauto4Nassim_brch=respiration for N assimilation
! CH2O=total CH2O production, 1/4 newly fixed C is used for N-assimilation (where is this number from?)
!
  
  NonStructalC4Growth_brch=RgroCO2_ltd/DMSHD
  CanopyNonstElm4Gros(ielmn)=AZMAX1(AMIN1(CanopyNonstElms_brch(ielmn,NB,NZ),NonStructalC4Growth_brch*(CNSHX+CNLFM+CNLFX*CNPG)))
  CanopyNonstElm4Gros(ielmp)=AZMAX1(AMIN1(CanopyNonstElms_brch(ielmp,NB,NZ),NonStructalC4Growth_brch*(CPSHX+CPLFM+CPLFX*CNPG)))  
  Rauto4Nassim_brch=AZMAX1(1.70_r8*CanopyNonstElm4Gros(ielmn)-0.025_r8*CH2O)
  
!
! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
!
! Rauto_brch=total C respiration
! RMNCS=maintenance respiration
! RCO2C=respiration from non-structural C
! RgroCO2_ltd=growth respiration limited by N,P
! RMxess_brch=excess maintenance respiration
! Rauto4Nassim_brch=respiration for N assimilation
! GrossCO2Fix_pft=total PFT CO2 fixation
! CO2F=total CO2 fixation
! GrossResp_pft,CanopyPlusNoduRespC_pft=total,above-ground PFT respiration
! CO2NetFix_pft=PFT net CO2 fixation
! Eco_GPP_col=ecosystem GPP
! ECO_ER_col=ecosystem respiration
! Eco_AutoR_col=total autotrophic respiration
!
  Rauto_brch=AMIN1(RMNCS,RCO2C)+RgroCO2_ltd+RMxess_brch+Rauto4Nassim_brch
  GrossCO2Fix_pft(NZ)=GrossCO2Fix_pft(NZ)+CO2F
  GrossResp_pft(NZ)=GrossResp_pft(NZ)-Rauto_brch
  CanopyPlusNoduRespC_pft(NZ)=CanopyPlusNoduRespC_pft(NZ)-Rauto_brch
  CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)+CO2F-Rauto_brch
  Eco_GPP_col=Eco_GPP_col+CO2F
  ECO_ER_col=ECO_ER_col-Rauto_brch
  Eco_AutoR_col=Eco_AutoR_col-Rauto_brch

  end associate
  end subroutine ComputRAutoAfEmergence

!------------------------------------------------------------------------------------------

  subroutine ComputRAutoB4Emergence(I,NB,NZ,TFN6_vr,DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,&
    CNLFX,CPLFX,ShootStructN,WFNG,WFNSG,CanopyNonstElm4Gros,CNPG,RCO2C,RMNCS,RMxess_brch,&
    NonStructalC4Growth_brch,CNRDM,Rauto4Nassim_brch)
  implicit none
  integer, intent(in) :: I,NB,NZ
  real(r8),intent(in) :: TFN6_vr(JZ1)
  real(r8), intent(in) :: DMSHD,CNLFM,CPLFM,CNSHX,CPSHX,CNLFX,CPLFX,ShootStructN,WFNG
  real(r8), intent(in) :: WFNSG
  real(r8), intent(out) :: CanopyNonstElm4Gros(NumPlantChemElms),rco2c
  real(r8), INTENT(OUT) :: CNPG,RMNCS,RMxess_brch
  real(r8), intent(out) :: NonStructalC4Growth_brch,CNRDM,Rauto4Nassim_brch
  real(r8) :: ZPOOLB,ZADDBM,CGROSM
  real(r8) :: FNP
  real(r8) :: PPOOLB
  real(r8) :: RCO2X,RCO2Y
  real(r8) :: RgroCO2_ltd   !oxygen & nutrient-limited growth respiration
  real(r8) :: Rauto_brch,RCO2CM
  real(r8) :: RCO2X_O2ulm,RCO2YM
  real(r8) :: RCO2GM
  real(r8) :: RCO2TM
  real(r8) :: SNCRM
! begin_execution
  associate(                                                                            &
    LeafPetoNonstElmConc_brch          =>  plt_biom%LeafPetoNonstElmConc_brch   , &
    CanopyNonstElms_brch                =>  plt_biom%CanopyNonstElms_brch   , &
    iPlantRootProfile_pft              =>  plt_pheno%iPlantRootProfile_pft  , &
    fTgrowRootP_vr                     =>  plt_pheno%fTgrowRootP_vr    , &
    RAutoRootO2Limter_pvr              =>  plt_rbgc%RAutoRootO2Limter_pvr     , &
    iPlantPhenolType_pft               =>  plt_pheno%iPlantPhenolType_pft  , &
    CO2NetFix_pft                      =>  plt_bgcr%CO2NetFix_pft     , &
    RootRespPotent_pvr                 =>  plt_rbgc%RootRespPotent_pvr    , &
    RCO2N_pvr                          =>  plt_rbgc%RCO2N_pvr    , &
    RCO2A_pvr                          =>  plt_rbgc%RCO2A_pvr    , &
    ZERO                               =>  plt_site%ZERO      , &
    NGTopRootLayer_pft                 =>  plt_morph%NGTopRootLayer_pft     , &
    C4PhotosynDowreg_brch              =>  plt_photo%C4PhotosynDowreg_brch     &
  )
  iter=iter+1
!
! N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
!
! CNPG=N,P constraint on growth respiration
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch
! CNKI,CPKI=nonstructural N,P inhibition constant on growth
!
  IF(LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZERO)THEN
    CNPG=AMIN1(LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmn,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CNKI), &
      LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmp,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmc,NB,NZ)*CPKI))
  ELSE
    CNPG=1.0_r8
  ENDIF
!
! RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
! NON-STRUCTURAL C:N:P, O2 UPTAKE
!
! RCO2CM,RCO2C=respiration from non-structural C unlimited,limited by O2
! VMXC=rate constant for nonstructural C oxidation in respiration (h-1)
! CPOOL=non-structural C mass
! fTgrowRootP_vr=temperature function for root growth
! WFNG=growth function of canopy water potential
! CNPG=N,P constraint on respiration
! C4PhotosynDowreg_brch=termination feedback inhibition on C3 CO2
! RAutoRootO2Limter_pvr=constraint by O2 consumption on all root processes
!
  RCO2CM=AZMAX1(VMXC*CanopyNonstElms_brch(ielmc,NB,NZ) &
    *fTgrowRootP_vr(NGTopRootLayer_pft(NZ),NZ))*CNPG*WFNG*C4PhotosynDowreg_brch(NB,NZ)
!  write(104,*)'RCO2CM',iter,RCO2CM,CanopyNonstElms_brch(ielmc,NB,NZ),NB,NZ,NGTopRootLayer_pft(NZ)  
  RCO2C=RCO2CM*RAutoRootO2Limter_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)
!
! MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
!
! RMNCS=maintenance respiration
! TFN6_vr=temperature function for root maintenance respiration
! ShootStructN=shoot structural N mass
! iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
! iPlantPhenolType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
! WFNG=growth function of canopy water potential
!
  RMNCS=AZMAX1(RmSpecPlant*TFN6_vr(NGTopRootLayer_pft(NZ))*ShootStructN)
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)).OR.&
    iPlantPhenolType_pft(NZ).EQ.iphenotyp_drouhtdecidu)THEN
    RMNCS=RMNCS*WFNG
  ENDIF
!
! GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
! IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
! RCO2X_O2ulm,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
! RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by N,P
! WFNSG=expansion,extension function of canopy water potential
! SNCRM,RMxess_brch=excess maintenance respiration unltd,ltd by O2
!
  RCO2X_O2ulm=RCO2CM-RMNCS
  RCO2X=RCO2C-RMNCS
  RCO2YM=AZMAX1(RCO2X_O2ulm)*WFNSG
  RCO2Y=AZMAX1(RCO2X)*WFNSG
  SNCRM=AZMAX1(-RCO2X_O2ulm)
  RMxess_brch=AZMAX1(-RCO2X)
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
! RAutoRootO2Limter_pvr=constraint by O2 consumption on all root processes
!
  IF(CNSHX.GT.0.0_r8.OR.CNLFX.GT.0.0_r8)THEN
    ZPOOLB=AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))
    PPOOLB=AZMAX1(CanopyNonstElms_brch(ielmp,NB,NZ))
    FNP=AMIN1(ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG),PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))
    IF(RCO2YM.GT.0.0_r8)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
    ELSE
      RCO2GM=0._r8
    ENDIF
    IF(RCO2Y.GT.0.0_r8)THEN
      RgroCO2_ltd=AMIN1(RCO2Y,FNP*RAutoRootO2Limter_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ))
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
! CGROSM,NonStructalC4Growth_brch=total non-structural C used in growth and respn unltd,ltd by O2
! RCO2GM,RgroCO2_ltd=growth respiration unltd,ltd by O2 and limited by N,P
! DMSHD=branch C respiration vs nonstructural C consumption
! ZADDBM,CanopyNonstElm4Gros(ielmn),CanopyNonstElm4Gros(ielmp)=nonstructural N,P unltd,ltd by O2 used in growth
! ZPOOL,PPOOL=nonstructural N,P mass
! CNSHX,CPSHX=N,P production vs nonstructural C consumption in rest of shoot
! CNLFM,CPLFM=min leaf N,P production vs nonstructural C consumption
! CNLFX,CPLFX=diff between min and max leaf N,P prodn vs nonstruct C consumption
! CNPG=N,P constraint on growth respiration
! CNRDM,Rauto4Nassim_brch=respiration for N assimilation unltd,ltd by O2
!
  CGROSM=RCO2GM/DMSHD
  NonStructalC4Growth_brch=RgroCO2_ltd/DMSHD
  ZADDBM=AZMAX1(CGROSM*(CNSHX+CNLFM+CNLFX*CNPG))
  CanopyNonstElm4Gros(ielmn)=AZMAX1(NonStructalC4Growth_brch*(CNSHX+CNLFM+CNLFX*CNPG))
  CanopyNonstElm4Gros(ielmp)=AZMAX1(NonStructalC4Growth_brch*(CPSHX+CPLFM+CPLFX*CNPG))
  CNRDM=AZMAX1(1.70_r8*ZADDBM)
  Rauto4Nassim_brch=AZMAX1(1.70_r8*CanopyNonstElm4Gros(ielmn))
!
! TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
! ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
!
! RCO2TM,Rauto_brch=total C respiration unltd,ltd by O2
! RMNCS=maintenance respiration
! RCO2GM,RgroCO2_ltd=growth respiration limited by N,P unltd,ltd by O2
! SNCRM,RMxess_brch=excess maintenance respiration unltd,ltd by O2
! CNRDM,Rauto4Nassim_brch=respiration for N assimilation unltd,ltd by O2
! RCO2A_pvr=total root respiration
! RootRespPotent_pvr,RCO2N_pvr=RCO2A_pvr unltd by O2,nonstructural C
!
  RCO2TM=RMNCS+RCO2GM+SNCRM+CNRDM
  Rauto_brch=RMNCS+RgroCO2_ltd+RMxess_brch+Rauto4Nassim_brch
  RootRespPotent_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RootRespPotent_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+RCO2TM
  RCO2N_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RCO2N_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)+Rauto_brch
  RCO2A_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=RCO2A_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)-Rauto_brch

  end associate
  end subroutine ComputRAutoB4Emergence

!------------------------------------------------------------------------------------------

  subroutine GrowLeavesOnBranch(NZ,NB,NNOD1,GrowthLeaf,ETOL,WFNS,ALLOCL)
  implicit none
  integer , intent(in) :: NB,NZ
  integer , intent(in) :: NNOD1
  real(r8), intent(in) :: GrowthLeaf(NumPlantChemElms)
  real(r8), intent(in) :: ETOL
  real(r8), intent(in) :: WFNS
  REAL(R8), INTENT(OUT):: ALLOCL
  integer :: K,NE,KK
  integer :: MXNOD,MNNOD,KNOD
  real(r8) :: GNOD
  REAL(R8) :: GrowthSLA,LeafAreaGrowth,SpecAreaLeafGrowth
  REAL(R8) :: GrowthChemElmt(NumPlantChemElms)

  associate(                                                                &
    NumCogrowNode                =>   plt_morph%NumCogrowNode ,   &        
    LeafAreaNode_brch            =>   plt_morph%LeafAreaNode_brch    , &    
    LeafAreaLive_brch            =>   plt_morph%LeafAreaLive_brch   , &    
    SLA1                         =>   plt_morph%SLA1    , &
    LeafElmntNode_brch           =>   plt_biom%LeafElmntNode_brch     , &    
    LeafProteinCNode_brch        =>   plt_biom%LeafProteinCNode_brch      , &    
    KHiestGroLeafNode_brch       =>   plt_pheno%KHiestGroLeafNode_brch   , &        
    ZEROL                        =>   plt_biom%ZEROL     , &    
    rCNNonstructRemob_pft        =>   plt_allom%rCNNonstructRemob_pft    , &
    rCPNonstructRemob_pft        =>   plt_allom%rCPNonstructRemob_pft     , &    
    FNOD                         =>   plt_allom%FNOD     , &
    PlantPopulation_pft          =>   plt_site%PlantPopulation_pft          &
  )
!  write(101,*)'grow leave',NB,NZ,GrowthLeaf(ielmc)
  IF(GrowthLeaf(ielmc).GT.0.0_r8)THEN
    MXNOD=KHiestGroLeafNode_brch(NB,NZ)
    MNNOD=MAX(NNOD1,MXNOD-NumCogrowNode(NZ)+1)
    MXNOD=MAX(MXNOD,MNNOD)
    KNOD=MXNOD-MNNOD+1
    GNOD=KNOD
    ALLOCL=1.0_r8/GNOD
    DO NE=1,NumPlantChemElms
      GrowthChemElmt(NE)=ALLOCL*GrowthLeaf(NE)
    ENDDO
    GrowthSLA=ALLOCL*FNOD(NZ)*NumCogrowNode(NZ)
!
!     GROWTH AT EACH CURRENT NODE
!
!     WGLF,WGLFN,WGLFP,LeafProteinCNode_brch=node leaf C,N,P,protein mass
!     GRO,GrowthChemElmt(ielmn),GrowthChemElmt(ielmp)=leaf C,N,P growth at each node
!     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!
    D490: DO KK=MNNOD,MXNOD
      K=MOD(KK,MaxNodesPerBranch1)
      IF(K.EQ.0.AND.KK.NE.0)K=MaxNodesPerBranch1
        DO NE=1,NumPlantChemElms
          LeafElmntNode_brch(NE,K,NB,NZ)=LeafElmntNode_brch(NE,K,NB,NZ)+GrowthChemElmt(NE)
        ENDDO
        LeafProteinCNode_brch(K,NB,NZ)=LeafProteinCNode_brch(K,NB,NZ)+ &
          AMIN1(GrowthChemElmt(ielmn)*rCNNonstructRemob_pft(NZ),GrowthChemElmt(ielmp)*rCPNonstructRemob_pft(NZ))
!
!         SPECIFIC LEAF AREA FUNCTION OF CURRENT LEAF MASS
!         AT EACH NODE
!
!         SpecAreaLeafGrowth=specific area of leaf growth
!         ETOL=coefficient for etoliation effects on expansion,extension
!         SLA1=growth in leaf area vs mass from PFT file
!         SLA2=parameter for calculating leaf area expansion
!         WGLF=leaf C mass
!         PP=PFT population
!         GrowthSLA=allocation of leaf area growth to each node
!         WFNS=turgor expansion,extension function
!         LeafAreaGrowth,GRO=leaf area,mass growth
!         LeafAreaLive_brch,LeafAreaNode_brch=branch,node leaf area
!
      SpecAreaLeafGrowth=ETOL*SLA1(NZ)*(AMAX1(ZEROL(NZ) &
        ,LeafElmntNode_brch(ielmc,K,NB,NZ))/(PlantPopulation_pft(NZ)*GrowthSLA))**SLA2*WFNS
      LeafAreaGrowth=GrowthChemElmt(ielmc)*SpecAreaLeafGrowth
      LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)+LeafAreaGrowth
!      if(NZ==1)THEN
!      WRITE(171,*)'leafarea growth',NB,LeafAreaGrowth
!      ELSE
!      WRITE(172,*)'leafarea growth',NB,LeafAreaGrowth
!      ENDIF
      LeafAreaNode_brch(K,NB,NZ)=LeafAreaNode_brch(K,NB,NZ)+LeafAreaGrowth
    ENDDO D490
  ENDIF
  end associate  
  END subroutine GrowLeavesOnBranch      
!------------------------------------------------------------------------------------------
  subroutine GrowPetioleOnBranch(NZ,NB,NNOD1,GrowthPetiole,ETOL,WFNS,ALLOCL)
  implicit none
  integer, intent(in) :: NZ,NNOD1,NB
  real(r8), intent(in) :: GrowthPetiole(NumPlantChemElms)
  REAL(R8), INTENT(IN) :: ETOL
  real(r8), intent(in) :: WFNS
  REAL(R8), INTENT(IN) :: ALLOCL
  integer :: MXNOD,MNNOD,NE,KK,K
  real(r8) :: GNOD,ALLOCS,GSSL,GROS
  REAL(R8) :: GrowthChemElmt(NumPlantChemElms)
  REAL(R8) :: SSL
  
  associate(                           &
    NumCogrowNode                =>   plt_morph%NumCogrowNode ,   &      
    PetioleLengthNode_brch       =>   plt_morph%PetioleLengthNode_brch   , &    
    PetoLen2Mass_pft             =>   plt_morph%PetoLen2Mass_pft    , &    
    SinePetioleAngle_pft         =>   plt_morph%SinePetioleAngle_pft   , &
    PetioleElmntNode_brch        =>  plt_biom%PetioleElmntNode_brch    , &    
    PetioleProteinCNode_brch     =>  plt_biom%PetioleProteinCNode_brch    , & 
    LeafElmntNode_brch           =>  plt_biom%LeafElmntNode_brch     , &       
    ZEROL                        =>  plt_biom%ZEROL     , &        
    rCNNonstructRemob_pft        =>  plt_allom%rCNNonstructRemob_pft    , &
    rCPNonstructRemob_pft        =>  plt_allom%rCPNonstructRemob_pft     , &    
    FNOD                         =>  plt_allom%FNOD     , &
    PlantPopulation_pft          =>  plt_site%PlantPopulation_pft        , &    
    KHiestGroLeafNode_brch       =>  plt_pheno%KHiestGroLeafNode_brch     &        
  )

  IF(GrowthPetiole(ielmc).GT.0.0_r8)THEN
    MXNOD=KHiestGroLeafNode_brch(NB,NZ)
    MNNOD=MAX(NNOD1,MXNOD-NumCogrowNode(NZ)+1)
    MXNOD=MAX(MXNOD,MNNOD)
    GNOD=MXNOD-MNNOD+1
    ALLOCS=1.0_r8/GNOD
    DO NE=1,NumPlantChemElms
      GrowthChemElmt(NE)=ALLOCS*GrowthPetiole(NE)
    ENDDO
    GSSL=ALLOCL*FNOD(NZ)*NumCogrowNode(NZ)
!
!       GROWTH AT EACH CURRENT NODE
!
!       PetioleElmntNode_brch,WGSHN,WGSHP,PetioleProteinCNode_brch=node petiole C,N,P,protein mass
!       GRO,GrowthChemElmt(ielmn),GrowthChemElmt(ielmp)=petiole C,N,P growth at each node
!       CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!
    D505: DO KK=MNNOD,MXNOD
      K=pMOD(KK,MaxNodesPerBranch1)

      DO NE=1,NumPlantChemElms
        PetioleElmntNode_brch(NE,K,NB,NZ)=PetioleElmntNode_brch(NE,K,NB,NZ)+GrowthChemElmt(NE)
      ENDDO
      PetioleProteinCNode_brch(K,NB,NZ)=PetioleProteinCNode_brch(K,NB,NZ) &
        +AMIN1(GrowthChemElmt(ielmn)*rCNNonstructRemob_pft(NZ),GrowthChemElmt(ielmp)*rCPNonstructRemob_pft(NZ))
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
    !   PetioleLengthNode_brch=petiole length
!
      IF(LeafElmntNode_brch(ielmc,K,NB,NZ).GT.0.0_r8)THEN
        SSL=ETOL*PetoLen2Mass_pft(NZ)*(AMAX1(ZEROL(NZ) &
          ,PetioleElmntNode_brch(ielmc,K,NB,NZ))/(PlantPopulation_pft(NZ)*GSSL))**SSL2*WFNS
        GROS=GrowthChemElmt(ielmc)/PlantPopulation_pft(NZ)*SSL
        PetioleLengthNode_brch(K,NB,NZ)=PetioleLengthNode_brch(K,NB,NZ)+GROS*SinePetioleAngle_pft(NZ)
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
  REAL(R8) :: GrowthChemElmt(NumPlantChemElms)

  integer :: NN,MXNOD,MNNOD,K1,K2,NE,KX,KK
  REAL(R8) :: ALLOCN,GNOD,SpecLenStalkGrowth
  real(r8) :: StalkLenGrowth

  associate(                                                       &
    NumCogrowNode                =>   plt_morph%NumCogrowNode,    &  
    NodeLenPergC                 =>   plt_morph%NodeLenPergC    , &      
    SineBranchAngle_pft          =>   plt_morph%SineBranchAngle_pft   , &    
    InternodeHeightDying_brch    =>   plt_morph%InternodeHeightDying_brch  , &
    InternodeHeightLive_brch     =>   plt_morph%InternodeHeightLive_brch  , &    
    InternodeStrutElms_brch        =>   plt_biom%InternodeStrutElms_brch     , &    
    StalkStrutElms_brch           =>   plt_biom%StalkStrutElms_brch   , &
    PlantPopulation_pft          =>   plt_site%PlantPopulation_pft        , &          
    iPlantCalendar_brch          =>   plt_pheno%iPlantCalendar_brch   , &              
    KHiestGroLeafNode_brch       =>   plt_pheno%KHiestGroLeafNode_brch     &        
  )

  IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).EQ.0)THEN
    NN=0
  ELSE
    NN=1
  ENDIF

  MXNOD=KHiestGroLeafNode_brch(NB,NZ)
  MNNOD=MAX(MIN(NN,MAX(NN,MXNOD-NumCogrowNode(NZ))),&
    KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+2)
  MXNOD=MAX(MXNOD,MNNOD)

  IF(GrowthStalk(ielmc).GT.0.0_r8)THEN
    GNOD=MXNOD-MNNOD+1
    ALLOCN=1.0_r8/GNOD
    DO NE=1,NumPlantChemElms
      GrowthChemElmt(NE)=ALLOCN*GrowthStalk(NE)
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
    SpecLenStalkGrowth=ETOL*NodeLenPergC(NZ)*(StalkStrutElms_brch(ielmc,NB,NZ)/PlantPopulation_pft(NZ))**SNL2
    StalkLenGrowth=GrowthChemElmt(ielmc)/PlantPopulation_pft(NZ)*SpecLenStalkGrowth
    KX=pMOD(MNNOD-1,MaxNodesPerBranch1)

!
!     GROWTH AT EACH CURRENT NODE
!
!     InternodeStrutElms_brch,WGNODN,WGNODP=node stalk C,N,P mass
!     GRO,GrowthChemElmt(ielmn),GrowthChemElmt(ielmp)=stalk C,N,P growth at each node
!     InternodeHeightDying_brch,InternodeHeightLive_brch=stalk height,stalk internode length
!     SineBranchAngle_pft=sine of stalk angle from horizontal from PFT file
!
    D510: DO KK=MNNOD,MXNOD
      K1=pMOD(KK,MaxNodesPerBranch1)            
      K2=pMOD(KK-1,MaxNodesPerBranch1)
      DO NE=1,NumPlantChemElms
        InternodeStrutElms_brch(NE,K1,NB,NZ)=InternodeStrutElms_brch(NE,K1,NB,NZ)+GrowthChemElmt(NE)
      ENDDO
      InternodeHeightDying_brch(K1,NB,NZ)=InternodeHeightDying_brch(K1,NB,NZ) &
        +StalkLenGrowth*SineBranchAngle_pft(NZ)
      IF(K1.NE.0)THEN
        InternodeHeightLive_brch(K1,NB,NZ)=InternodeHeightDying_brch(K1,NB,NZ) &
          +InternodeHeightLive_brch(K2,NB,NZ)
      ELSE
        InternodeHeightLive_brch(K1,NB,NZ)=InternodeHeightDying_brch(K1,NB,NZ)
      ENDIF
    ENDDO D510
  ENDIF
  END associate  
  end subroutine GrowStalkOnBranch

!------------------------------------------------------------------------------------------

  subroutine RemobilizeBranch(NZ,NB,BegRemoblize,IFLGY,RCCC,RCCN,RCCP,RMxess_brch)
  implicit none
  integer, intent(in) :: NZ,NB
  integer, intent(in) :: IFLGY,BegRemoblize
  real(r8), intent(in) :: RCCC,RCCN,RCCP
  real(r8), intent(inout) :: RMxess_brch
  integer :: NumRemobLeafNodes  
  real(r8) :: XKSNC  
  real(r8) :: SNCX
  real(r8) :: SNCF  
  integer  :: NBZ(MaxNumBranches)
  real(r8) :: XFRE(1:NumPlantChemElms)
  REAL(R8) :: RCO2V,SNCZ
  integer :: NBY,NBX,NBL
  integer  :: NBK,NE

  associate(                                                          &
    fTgrowCanP                =>  plt_pheno%fTgrowCanP              , &  
    iPlantTurnoverPattern_pft =>  plt_pheno%iPlantTurnoverPattern_pft   , &    
    HoursDoingRemob_brch      =>  plt_pheno%HoursDoingRemob_brch    , &
    KHiestGroLeafNode_brch    =>  plt_pheno%KHiestGroLeafNode_brch   , &    
    iPlantPhenolPattern_pft=>  plt_pheno%iPlantPhenolPattern_pft   , &   
    iPlantRootProfile_pft     =>  plt_pheno%iPlantRootProfile_pft   , &     
    iPlantBranchState_brch    =>  plt_pheno%iPlantBranchState_brch    , &
    KLowestGroLeafNode_brch   =>  plt_pheno%KLowestGroLeafNode_brch   , &    
    NumOfBranches_pft         =>  plt_morph%NumOfBranches_pft     , &    
    MainBranchNum_pft         =>  plt_morph%MainBranchNum_pft      , &         
    BranchNumber_brch         =>  plt_morph%BranchNumber_brch    , &        
    CanopyNonstElms_brch         =>  plt_biom%CanopyNonstElms_brch    , &   
    LeafPetolBiomassC_brch    =>  plt_biom%LeafPetolBiomassC_brch     , &        
    StalkRsrvElms_brch       =>  plt_biom%StalkRsrvElms_brch   , &      
    ZEROP                     =>  plt_biom%ZEROP       &  
  ) 
  !     REMOBILIZATION OF STALK RESERVE C,N,P IF GROWTH RESPIRATION < 0

!
!     RMxess_brch=excess maintenance respiration
!     WTRSVB=stalk reserve C mass
!     RCO2V=remobilization of stalk reserve C
!     VMXC=rate constant for nonstructural C oxidation in respiration
!     fTgrowCanP=temperature function for canopy growth
!
  IF(BegRemoblize.EQ.ifalse)THEN
    IF(RMxess_brch.GT.0.0_r8.AND.StalkRsrvElms_brch(ielmc,NB,NZ).GT.0.0_r8)THEN
      RCO2V=AMIN1(RMxess_brch,VMXC*StalkRsrvElms_brch(ielmc,NB,NZ)*fTgrowCanP(NZ))
      StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)-RCO2V
      RMxess_brch=RMxess_brch-RCO2V
    ENDIF
  ENDIF
!
!       TOTAL REMOBILIZATION = GROWTH RESPIRATION < 0 + DECIDUOUS LEAF
!       FALL DURING AUTUMN + REMOBILZATION DURING GRAIN FILL IN ANNUALS
!
!       iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!       IFLGY,BegRemoblize=remobilization flags
!       SNCZ=phenologically-driven respiration senescence during late-season
!       FXFB=rate constant for plant-storage nonstructural C,N,P exchange
!       iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+petiole,2=none,3=between 1,2
!       LeafPetolBiomassC_brch=leaf+petiole mass
!       FLGZ=control rate of remobilization
!       Hours4FullSenescence=number of hours until full senescence after physl maturity
!       SNCX=total senescence respiration
!       KHiestGroLeafNode_brch,KLowestGroLeafNode_brch=integer of highest,lowest leaf number currently growing
!       KSNC=number of nodes undergoing remobilization
!       SNCF=ratio of phenologically-driven vs total senescence respiration
!
  IF(BegRemoblize.EQ.itrue .AND. IFLGY.EQ.itrue .AND. iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
    SNCZ=FXFB(iPlantTurnoverPattern_pft(NZ))*LeafPetolBiomassC_brch(NB,NZ) &
      *AMIN1(1.0_r8,HoursDoingRemob_brch(NB,NZ)/Hours4FullSenescence)
  ELSE
    SNCZ=0._r8
  ENDIF
  SNCX=RMxess_brch+SNCZ

  IF(SNCX.GT.ZEROP(NZ))THEN
    SNCF=SNCZ/SNCX
    NumRemobLeafNodes=INT(0.5_r8*(KHiestGroLeafNode_brch(NB,NZ)-KLowestGroLeafNode_brch(NB,NZ)))+1
    XKSNC=NumRemobLeafNodes
!
!     TRANSFER NON-STRUCTURAL C,N,P FROM BRANCHES TO MAIN STEM
!     IF MAIN STEM POOLS ARE DEPLETED
!
!     iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantRootProfile_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!     XFRE(ielmc),XFRE(ielmn),XFRE(ielmc)=nonstructural C,N,P transfer
!
    IF(iPlantTurnoverPattern_pft(NZ).NE.0.AND.is_plant_treelike(iPlantRootProfile_pft(NZ)) &
      .AND.NB.EQ.MainBranchNum_pft(NZ).AND.isclose(SNCF,0._r8))THEN
      NBY=0
      D584: DO NBL=1,NumOfBranches_pft(NZ)
        NBZ(NBL)=0
      ENDDO D584

      D586: DO NBL=1,NumOfBranches_pft(NZ)
        NBX=KHiestGroLeafNode_brch(NB,NZ)
        D585: DO NBK=1,NumOfBranches_pft(NZ)
          IF(iPlantBranchState_brch(NBK,NZ).EQ.iLive .AND. NBK.NE.MainBranchNum_pft(NZ) &
            .AND. BranchNumber_brch(NBK,NZ).LT.NBX .AND. BranchNumber_brch(NBK,NZ).GT.NBY)THEN
            NBZ(NBL)=NBK
            NBX=BranchNumber_brch(NBK,NZ)
          ENDIF
        ENDDO D585
        IF(NBZ(NBL).NE.0)THEN
          NBY=BranchNumber_brch(NBZ(NBL),NZ)
        ENDIF
      ENDDO D586

      D580: DO NBL=1,NumOfBranches_pft(NZ)
        IF(NBZ(NBL).NE.0)THEN
          !should KK be NBX?
          IF(BranchNumber_brch(NBZ(NBL),NZ).NE.MainBranchNum_pft(NZ))THEN
            IF(CanopyNonstElms_brch(ielmc,NBZ(NBL),NZ).GT.0._r8)THEN
              XFRE(ielmc)=1.0E-02_r8*AMIN1(SNCX,CanopyNonstElms_brch(ielmc,NBZ(NBL),NZ))
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
              if(CanopyNonstElms_brch(NE,NBZ(NBL),NZ)<0._r8)then
                write(*,*)'3677CanopyNonstElms_brch(NE,NBZ(NBL),NZ)',NE,NBZ(NBL),CanopyNonstElms_brch(NE,NBZ(NBL),NZ)+XFRE(NE)&
                  ,XFRE(NE)
                stop
              endif
            ENDDO
            CanopyNonstElms_brch(ielmc,MainBranchNum_pft(NZ),NZ)=CanopyNonstElms_brch(ielmc,MainBranchNum_pft(NZ),NZ) &
              +XFRE(ielmc)*SNCF
            DO NE=2,NumPlantChemElms
              CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)=CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)+XFRE(NE)
              if(CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)<0._r8)then
                write(*,*)'3688CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)',NE &
                  ,CanopyNonstElms_brch(NE,MainBranchNum_pft(NZ),NZ)-XFRE(NE),XFRE(NE)
                stop
              endif
            ENDDO
            SNCX=SNCX-XFRE(ielmc)
            IF(SNCX.LE.0.0_r8)exit
          ENDIF
        ENDIF
      ENDDO D580
    ENDIF
    IF(SNCX.GT.0.0_r8) &
      call RemobilizeLeafLayers(NumRemobLeafNodes,NB,nz,XKSNC,SNCX,RCCC,RCCN,RCCP,SNCF)
  ENDIF
  end associate
  end subroutine RemobilizeBranch
!------------------------------------------------------------------------------------------
  subroutine SenescenceBranch(NZ,NB,RCCC,RCCN,RCCP)
  implicit none
  integer, intent(in) :: NZ,NB
  real(r8), intent(in) :: RCCC,RCCN,RCCP
  INTEGER :: K,KMinGroingLeafNodeNum,M,NE
  real(r8) :: FSNC
  real(r8) :: FSNCL
  real(r8) :: FSNCS
  
  associate(                                                          &
    icwood                       =>  pltpar%icwood    , &  
    ifoliar                      =>  pltpar%ifoliar   , &        
    inonfoliar                   =>  pltpar%inonfoliar  , &    
    k_woody_litr                 =>  pltpar%k_woody_litr, &
    k_fine_litr                  =>  pltpar%k_fine_litr  , &    
    CFOPE                        =>  plt_soilchem%CFOPE , &    
    FWODLE                       =>  plt_allom%FWODLE   , &    
    FWODBE                       =>  plt_allom%FWODBE   , &    
    rCNNonstructRemob_pft        =>  plt_allom%rCNNonstructRemob_pft    , &
    rCPNonstructRemob_pft        =>  plt_allom%rCPNonstructRemob_pft     , &    
    LitfalStrutElms_pvr          =>  plt_bgcr%LitfalStrutElms_pvr      , &            
    SenecStalkStrutElms_brch      =>  plt_biom%SenecStalkStrutElms_brch   , &    
    ZEROP                        =>  plt_biom%ZEROP     , &    
    PetioleElmntNode_brch        =>  plt_biom%PetioleElmntNode_brch    , &    
    LeafProteinCNode_brch        =>  plt_biom%LeafProteinCNode_brch      , &    
    CanopyNonstElms_brch            =>  plt_biom%CanopyNonstElms_brch    , &    
    LeafStrutElms_brch            =>  plt_biom%LeafStrutElms_brch   , &    
    PetioleProteinCNode_brch     =>  plt_biom%PetioleProteinCNode_brch    , &
    PetioleChemElmRemob_brch     =>  plt_biom%PetioleChemElmRemob_brch   , &    
    InternodeStrutElms_brch        =>  plt_biom%InternodeStrutElms_brch    , &    
    LeafChemElmRemob_brch        =>  plt_biom%LeafChemElmRemob_brch    , & 
    LeafElmntNode_brch           =>  plt_biom%LeafElmntNode_brch     , & 
    PetoleStrutElms_brch           =>  plt_biom%PetoleStrutElms_brch  , &         
    CanPBranchHeight             =>  plt_morph%CanPBranchHeight  , &         
    InternodeHeightDying_brch    =>  plt_morph%InternodeHeightDying_brch  , &    
    LeafAreaNode_brch            =>  plt_morph%LeafAreaNode_brch   , &    
    LeafAreaDying_brch           =>  plt_morph%LeafAreaDying_brch   , &
    LeafAreaLive_brch            =>  plt_morph%LeafAreaLive_brch   , &
    PetioleLengthNode_brch       =>  plt_morph%PetioleLengthNode_brch   , &    
    doSenescence_brch            =>  plt_pheno%doSenescence_brch    , &    
    doRemobilization_brch        =>  plt_pheno%doRemobilization_brch    , &    
    PetioleChemElmRemobFlx_brch  =>  plt_pheno%PetioleChemElmRemobFlx_brch    , &    
    RefLeafAppearRate_pft        =>  plt_pheno%RefLeafAppearRate_pft    , &    
    KHiestGroLeafNode_brch       =>  plt_pheno%KHiestGroLeafNode_brch      , &    
    LeafElmntRemobFlx_brch       =>  plt_pheno%LeafElmntRemobFlx_brch    , &    
    fTgrowCanP                   =>  plt_pheno%fTgrowCanP                &  
  )    
  IF(doSenescence_brch(NB,NZ).EQ.itrue)THEN
    KMinGroingLeafNodeNum=MAX(0,KHiestGroLeafNode_brch(NB,NZ)-MaxNodesPerBranch1+1)
    IF(KMinGroingLeafNodeNum.GT.0)THEN
      K=pMOD(KMinGroingLeafNodeNum,MaxNodesPerBranch1)               
      FSNC=fTgrowCanP(NZ)*RefLeafAppearRate_pft(NZ)
!
      !   REMOBILIZATION OF LEAF C,N,P ALSO DEPENDS ON STRUCTURAL C:N:P
      !
      !   doRemobilization_brch=flag for remobilization
      !   WGLF,WGLFN,WGLFP=node leaf C,N,P mass
      !   LeafAreaNode_brch=node leaf area
      !   ElmntRemobFallingLeaf(ielmc)X,ElmntRemobFallingLeaf(ielmn)X,ElmntRemobFallingLeaf(ielmp)X=remobilization of C,N,P from senescing leaf
!
      IF(doRemobilization_brch(NB,NZ).EQ.itrue)THEN
        DO NE=1,NumPlantChemElms
          LeafChemElmRemob_brch(NE,NB,NZ)=AZMAX1(LeafElmntNode_brch(NE,K,NB,NZ))
        ENDDO
        LeafAreaDying_brch(NB,NZ)=AZMAX1(LeafAreaNode_brch(K,NB,NZ))
        IF(LeafChemElmRemob_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
          LeafElmntRemobFlx_brch(ielmc,NB,NZ)=LeafChemElmRemob_brch(ielmc,NB,NZ)*RCCC
          LeafElmntRemobFlx_brch(ielmn,NB,NZ)=LeafChemElmRemob_brch(ielmn,NB,NZ)*(RCCN+(1.0_r8-RCCN)*RCCC)
          LeafElmntRemobFlx_brch(ielmp,NB,NZ)=LeafChemElmRemob_brch(ielmp,NB,NZ)*(RCCP+(1.0_r8-RCCP)*RCCC)
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
        .AND.LeafChemElmRemob_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
        FSNCL=AZMAX1(LeafElmntNode_brch(ielmc,K,NB,NZ)/LeafChemElmRemob_brch(ielmc,NB,NZ))
      ELSE
        FSNCL=FSNC
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
          LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
            +CFOPE(NE,icwood,M,NZ)*FSNCL*(LeafChemElmRemob_brch(NE,NB,NZ)-LeafElmntRemobFlx_brch(NE,NB,NZ)) &
            *FWODLE(NE,k_woody_litr)
            
          LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +CFOPE(NE,ifoliar,M,NZ)*FSNCL*(LeafChemElmRemob_brch(NE,NB,NZ)-LeafElmntRemobFlx_brch(NE,NB,NZ)) &
            *FWODLE(NE,k_fine_litr)
        ENDDO D6300
      ENDDO
!
!       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
!
!       FSNCL=fraction of lowest leaf to be remobilized
!       LeafAreaLive_brch,LeafAreaDying_brch=branch living,senescing leaf area
!       WTLFB,WTLFBN,WTLFBP,WGLFX,WGLFNX,WGLFPX=C,N,P mass in living,senescing leaf
!       LeafProteinCNode_brch=leaf protein mass
!       CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!       CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!       ElmntRemobFallingLeaf(ielmc)X,ElmntRemobFallingLeaf(ielmn)X,ElmntRemobFallingLeaf(ielmp)X=remobilization of C,N,P from senescing leaf
!
      LeafAreaLive_brch(NB,NZ)=LeafAreaLive_brch(NB,NZ)-FSNCL*LeafAreaDying_brch(NB,NZ)
      DO NE=1,NumPlantChemElms
        LeafStrutElms_brch(NE,NB,NZ)=LeafStrutElms_brch(NE,NB,NZ)-FSNCL*LeafChemElmRemob_brch(NE,NB,NZ)
        LeafElmntNode_brch(NE,K,NB,NZ)=LeafElmntNode_brch(NE,K,NB,NZ)-FSNCL*LeafChemElmRemob_brch(NE,NB,NZ)
        CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+FSNCL*LeafElmntRemobFlx_brch(NE,NB,NZ)
        if(CanopyNonstElms_brch(NE,NB,NZ)<0._r8)then
          write(*,*)'3829CanopyNonstElms_brch(NE,NB,NZ)',NE,NB,CanopyNonstElms_brch(NE,NB,NZ),FSNCL*LeafElmntRemobFlx_brch(NE,NB,NZ)
          stop
        endif
      ENDDO
      LeafAreaNode_brch(K,NB,NZ)=LeafAreaNode_brch(K,NB,NZ)-FSNCL*LeafAreaDying_brch(NB,NZ)
      LeafProteinCNode_brch(K,NB,NZ)=AZMAX1(LeafProteinCNode_brch(K,NB,NZ) &
        -FSNCL*AMAX1(LeafChemElmRemob_brch(ielmn,NB,NZ)*rCNNonstructRemob_pft(NZ) &
        ,LeafChemElmRemob_brch(ielmp,NB,NZ)*rCPNonstructRemob_pft(NZ)))
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
        CanPBranchHeight(NB,NZ)=AZMAX1(PetioleLengthNode_brch(K,NB,NZ))
        IF(PetioleChemElmRemob_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
          PetioleChemElmRemobFlx_brch(ielmc,NB,NZ)=RCCC*PetioleChemElmRemob_brch(ielmc,NB,NZ)
          PetioleChemElmRemobFlx_brch(ielmn,NB,NZ)=PetioleChemElmRemob_brch(ielmn,NB,NZ) &
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
        InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ)=0._r8
        InternodeHeightDying_brch(K,NB,NZ)=0._r8
      ENDIF
!
!       FRACTION OF CURRENT SHEATH TO BE REMOBILIZED
!
!       FSNCS=fraction of lowest petiole to be remobilized
!
      IF(FSNC*PetioleChemElmRemob_brch(ielmc,NB,NZ).GT.PetioleElmntNode_brch(ielmc,K,NB,NZ) &
        .AND.PetioleChemElmRemob_brch(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
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
          LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ) &
            +CFOPE(NE,icwood,M,NZ)*FSNCS*(PetioleChemElmRemob_brch(NE,NB,NZ)&
            -PetioleChemElmRemobFlx_brch(NE,NB,NZ))*FWODBE(NE,k_woody_litr)

          LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ) &
            +CFOPE(NE,inonfoliar,M,NZ)*FSNCS*(PetioleChemElmRemob_brch(NE,NB,NZ)&
            -PetioleChemElmRemobFlx_brch(NE,NB,NZ))*FWODBE(NE,k_fine_litr)
        ENDDO    
      ENDDO D6305      
!
!       UPDATE STATE VARIABLES FOR REMOBILIZATION AND LitrFall
!
!       FSNCS=fraction of lowest petiole to be remobilized
!       PetioleLengthNode_brch,CanPBranchHeight=living,senescing petiole length
!       WTSHB,WTSHBN,WTSHBP,WGSHEX,WGSHNX,WGSHPX=C,N,P mass in living,senescing petiole
!       PetioleProteinCNode_brch=petiole protein mass
!       CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios from startq.f
!       CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
!       RCES(ielmc)X,RCES(ielmn)X,RCES(ielmp)X=remobilization of C,N,P from senescing petiole
!
      DO NE=1,NumPlantChemElms
        PetoleStrutElms_brch(NE,NB,NZ)=PetoleStrutElms_brch(NE,NB,NZ)-FSNCS &
          *PetioleChemElmRemob_brch(NE,NB,NZ)
        PetioleElmntNode_brch(NE,K,NB,NZ)=PetioleElmntNode_brch(NE,K,NB,NZ)-FSNCS &
          *PetioleChemElmRemob_brch(NE,NB,NZ)
        CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+FSNCS*PetioleChemElmRemobFlx_brch(NE,NB,NZ)
      ENDDO
      PetioleLengthNode_brch(K,NB,NZ)=PetioleLengthNode_brch(K,NB,NZ)-FSNCS*CanPBranchHeight(NB,NZ)

      PetioleProteinCNode_brch(K,NB,NZ)=AZMAX1(PetioleProteinCNode_brch(K,NB,NZ) &
        -FSNCS*AMAX1(PetioleChemElmRemob_brch(ielmn,NB,NZ)*rCNNonstructRemob_pft(NZ) &
        ,PetioleChemElmRemob_brch(ielmp,NB,NZ)*rCPNonstructRemob_pft(NZ)))

    ENDIF
  ENDIF
  END associate
  end subroutine SenescenceBranch  

end module PlantBranchMod
