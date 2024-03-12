module grosubsMod
!!
! Description:
! module for plant biological transformations
  use minimathmod, only : safe_adb,AZMAX1
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use EcoSiMParDataMod, only : pltpar
  use GrosubPars
  use PlantAPIData
  use PhotoSynsMod
  use PlantMathFuncMod
  use RootMod, only : RootBGCModel
  use LitterFallMod
  use PlantBranchMod
  implicit none

  private

! end_include_section

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
! DIMENSION VCO2(400,366,05)
!
!     RTSK=relative primary root sink strength 0.25=shallow,4.0=deep root profile
!     FXRN=rate constant for plant-bacteria nonstructl C,N,P exchange (h-1)
!     FXFB=rate constant for leaf-storage nonstructl C,N,P exchange (h-1)
!     FXFR=rate constant for root-storage nonstructl C,N,P exchange (h-1)
!     FPART=parameter for allocating nonstructural C to shoot organs
!     FXSH,FXRT=shoot-root partitioning of storage C during leafout
!     FRSV=rate constant for remobiln of storage C,N,P during leafout (h-1)
!     FXFY,FXFZ=rate constant for leaf-reserve nonstructural C,N,P exchange (h-1)
!     EFIRE=combustion  of N,P relative to C
!     PSILY=canopy water potential below which leafoff is induced (MPa)
!     HourFailGrainFill_brchY=number of hours after physiol maturity required for senescence
!     ATRPX=number of hours required to initiate remobilization of storage C for leafout
!     GVMX=specific oxidation rate of nonstructural C during leafout at 25 C
!     FracHour4LeafoffRemob=fraction of hours required for leafoff to initiate remobilization
!

  integer,private :: curday,curhour
  public :: grosubs
  public :: InitGrosub
  contains

  subroutine InitGrosub(NumGrowthStages,MaxNumRootAxes)

  implicit none
  integer, intent(out) :: NumGrowthStages,MaxNumRootAxes

  call InitVegPars(pltpar)
  NumGrowthStages = pltpar%NumGrowthStages
  MaxNumRootAxes = pltpar%MaxNumRootAxes


  end subroutine InitGrosub
!------------------------------------------------------------------------------------------

  subroutine grosubs(I,J)
!
!     THIS subroutine CALCULATES ALL PLANT BIOLOGICAL TRANSFORMATIONS
!
  use PlantDisturbsMod, only : RemoveBiomassByDisturbance
  implicit none
  integer, intent(in) :: I, J

  real(r8) :: CanopyHeight_copy(JP1)
  integer :: L,K,M
  integer :: NZ,NE
  real(r8) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
! begin_execution
  associate(                            &
    IsPlantActive_pft    => plt_pheno%IsPlantActive_pft   , &
    NP       => plt_site%NP       , &
    NP0      => plt_site%NP0      , &
    MaxNumRootLays       => plt_site%MaxNumRootLays       , &
    CO2NetFix_pft     => plt_bgcr%CO2NetFix_pft     , &
    LitterFallChemElmnt_pft    => plt_bgcr%LitterFallChemElmnt_pft    , &
    LitterFallChemElmnt_pvr     => plt_bgcr%LitterFallChemElmnt_pvr     , &
    CanopyHeight_pft      => plt_morph%CanopyHeight_pft       &
  )
!     TOTAL AGB FOR GRAZING IN LANDSCAPE SECTION
!
!
!     INITIALIZE SENESCENCE ARRAYS
!

  D9980: DO NZ=1,NP0
    D1: DO L=0,MaxNumRootLays
      DO K=1,pltpar%NumOfPlantLitrCmplxs
        DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitterFallChemElmnt_pvr(NE,M,K,L,NZ)=0._r8
          ENDDO
        ENDDO
      ENDDO
    ENDDO D1
    LitterFallChemElmnt_pft(1:NumPlantChemElms,NZ)=0._r8
    CO2NetFix_pft(NZ)=0._r8
    CanopyHeight_copy(NZ)=CanopyHeight_pft(NZ)
    CanopyHeight_pft(NZ)=0._r8
  ENDDO D9980
!
!     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
!
  D9985: DO NZ=1,NP

! IsPlantActive_pft= flag for living pft
    IF(IsPlantActive_pft(NZ).EQ.iPlantIsActive)THEN
      call GrowPlant(I,J,NZ,CanopyHeight_copy,ShootNonstructC_brch)
    ENDIF

!   HARVEST STANDING DEAD
    call RemoveBiomassByDisturbance(I,J,NZ,ShootNonstructC_brch)
  ENDDO D9985
!
! TRANSFORMATIONS IN LIVING OR DEAD PLANT POPULATIONS
  call LiveDeadTransformation(I,J)
  end associate
  END subroutine grosubs

!------------------------------------------------------------------------------------------

  subroutine LiveDeadTransformation(I,J)
  implicit none
  integer, intent(in) :: I,J

  integer :: L,K,NZ,M,NE,NB
  real(r8) :: XFRC,XFRN,XFRP,XFRE
!     begin_execution

  associate(                       &
    k_fine_litr=> pltpar%k_fine_litr ,&
    k_woody_litr=> pltpar%k_woody_litr,&
    iDayPlanting_pft   => plt_distb%iDayPlanting_pft   , &
    iYearPlanting_pft    => plt_distb%iYearPlanting_pft    , &
    EcoHavstElmntCum_pft  => plt_distb%EcoHavstElmntCum_pft  , &
    EcoHavstElmnt_pft   => plt_distb%EcoHavstElmnt_pft   , &
    PO4byFire_pft   => plt_distb%PO4byFire_pft   , &
    N2ObyFire_pft   => plt_distb%N2ObyFire_pft   , &
    CH4ByFire_pft   => plt_distb%CH4ByFire_pft   , &
    CO2ByFire_pft   => plt_distb%CO2ByFire_pft   , &
    NH3byFire_pft   => plt_distb%NH3byFire_pft   , &
    StandingDeadKCompChemElms_pft  => plt_biom%StandingDeadKCompChemElms_pft   , &
    RootElmnts_pft   => plt_biom%RootElmnts_pft    , &
    ShootChemElms_pft  => plt_biom%ShootChemElms_pft  , &
    NoduleChemElms_pft   => plt_biom%NoduleChemElms_pft    , &
    NonstructalElmnts_pft   => plt_biom%NonstructalElmnts_pft    , &
    StandingDeadChemElms_pft  => plt_biom%StandingDeadChemElms_pft   , &
    fTgrowCanP    => plt_pheno%fTgrowCanP    , &
    NetCumElmntFlx2Plant_pft   => plt_pheno%NetCumElmntFlx2Plant_pft   , &
    IsPlantActive_pft   => plt_pheno%IsPlantActive_pft   , &
    iPlantMorphologyType_pft  => plt_pheno%iPlantMorphologyType_pft  , &
    doInitPlant_pft   => plt_pheno%doInitPlant_pft   , &
    doPlantLeafOut_brch   => plt_pheno%doPlantLeafOut_brch   , &
    iPlantTurnoverPattern_pft  => plt_pheno%iPlantTurnoverPattern_pft  , &
    HourThreshold4LeafOut_brch   => plt_pheno%HourThreshold4LeafOut_brch   , &
    Hours4Leafout_brch    => plt_pheno%Hours4Leafout_brch    , &
    ElmntBalanceCum_pft    => plt_site%ElmntBalanceCum_pft     , &
    NP0     => plt_site%NP0      , &
    MaxNumRootLays      => plt_site%MaxNumRootLays       , &
    iYearCurrent    => plt_site%iYearCurrent     , &
    LitterFallChemElmnt_pvr    => plt_bgcr%LitterFallChemElmnt_pvr     , &
    NetPrimaryProductvity_pft    => plt_bgcr%NetPrimaryProductvity_pft     , &
    PlantN2FixCum_pft  => plt_bgcr%PlantN2FixCum_pft   , &
    LitterFallChemElmnt_pft   => plt_bgcr%LitterFallChemElmnt_pft    , &
    NH3EmiCum_pft   => plt_bgcr%NH3EmiCum_pft    , &
    LitrfallChemElms_pft   => plt_bgcr%LitrfallChemElms_pft    , &
    SurfLitrfallChemElms_pft   => plt_bgcr%SurfLitrfallChemElms_pft    , &
    GrossResp_pft   => plt_bgcr%GrossResp_pft    , &
    GrossCO2Fix_pft   => plt_bgcr%GrossCO2Fix_pft    , &
    PlantExudChemElmCum_pft  => plt_rbgc%PlantExudChemElmCum_pft   , &
    CumSoilThickness  => plt_site%CumSoilThickness   , &
    PlantinDepth  => plt_morph%PlantinDepth  , &
    NumOfBranches_pft     => plt_morph%NumOfBranches_pft       &
  )
  D9975: DO NZ=1,NP0
!
!     ACTIVATE DORMANT SEEDS
!
    D205: DO NB=1,NumOfBranches_pft(NZ)
      IF(doInitPlant_pft(NZ).EQ.itrue)THEN
        IF(doPlantLeafOut_brch(NB,NZ).EQ.iEnable.AND.Hours4Leafout_brch(NB,NZ).GE.HourThreshold4LeafOut_brch(NB,NZ))THEN
          iDayPlanting_pft(NZ)=I
          iYearPlanting_pft(NZ)=iYearCurrent
          PlantinDepth(NZ)=0.005_r8+CumSoilThickness(0)
          !mark plant as initialized
          doInitPlant_pft(NZ)=ifalse
        ENDIF
      ENDIF
    ENDDO D205
!
!     LITTERFALL FROM STANDING DEAD
!
!     XFRC,XFRN,XFRP=litterfall from standing dead
!     fTgrowCanP=temperature function for canopy growth
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     CSNC,ZSNC,PSNC=C,N,P litterfall
!
    DO NE=1,NumPlantChemElms
      D6235: DO M=1,jsken
        XFRE=1.5814E-05_r8*fTgrowCanP(NZ)*StandingDeadKCompChemElms_pft(NE,M,NZ)
        IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.iPlantMorphologyType_pft(NZ).LE.1)THEN
          LitterFallChemElmnt_pvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pvr(NE,M,k_fine_litr,0,NZ)+XFRE
        ELSE
          LitterFallChemElmnt_pvr(NE,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pvr(NE,M,k_woody_litr,0,NZ)+XFRE
        ENDIF
        StandingDeadKCompChemElms_pft(NE,M,NZ)=StandingDeadKCompChemElms_pft(NE,M,NZ)-XFRE
      ENDDO D6235
    ENDDO
!
!     ACCUMULATE TOTAL SURFACE, SUBSURFACE LITTERFALL
!
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P litterfall
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P litterfall
!     HCSNC,HZSNC,HPSNC=hourly C,N,P litterfall
!
    DO K=1,pltpar%NumOfPlantLitrCmplxs
      DO NE=1,NumPlantChemElms
        D6430: DO M=1,jsken
          SurfLitrfallChemElms_pft(NE,NZ)=SurfLitrfallChemElms_pft(NE,NZ)+LitterFallChemElmnt_pvr(NE,M,K,0,NZ)
          D8955: DO L=0,MaxNumRootLays
            LitterFallChemElmnt_pft(NE,NZ)=LitterFallChemElmnt_pft(NE,NZ)+LitterFallChemElmnt_pvr(NE,M,K,L,NZ)
            LitrfallChemElms_pft(NE,NZ)=LitrfallChemElms_pft(NE,NZ)+LitterFallChemElmnt_pvr(NE,M,K,L,NZ)
          ENDDO D8955
        ENDDO D6430
      enddo
    ENDDO
!
!     TOTAL STANDING DEAD
!
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
    DO NE=1,NumPlantChemElms
      StandingDeadChemElms_pft(NE,NZ)=sum(StandingDeadKCompChemElms_pft(NE,1:jsken,NZ))
    ENDDO

!
!     PLANT C BALANCE = TOTAL C STATE VARIABLES + TOTAL
!     AUTOTROPHIC RESPIRATION + TOTAL LITTERFALL - TOTAL EXUDATION
!     - TOTAL CO2 FIXATION
!
!     BALC=PFT C balance
!     WTSHT,WTRT,WTND,WTRVC,WTSTG=PFT shoot,root,bacteria,storage,standing dead C
!     NetPrimaryProductvity_pft=cumulative PFT NPP
!     TCSNC=cumulative PFT C litterfall
!     TCUPTK=cumulative PFT root-soil C exchange
!     RSETC=cumulative C balance from previous year
!     THVSTC=cumulative PFT C removed from ecosystem from previous year
!     HVSTC=total PFT C removed from ecosystem in current year
!     CO2ByFire_pft,CH4ByFire_pft=CO2,CH4 emission from disturbance
!
    NetPrimaryProductvity_pft(NZ)=GrossCO2Fix_pft(NZ)+GrossResp_pft(NZ)

    IF(IsPlantActive_pft(NZ).EQ.iPlantIsActive)THEN
      !check for living plant
      DO NE=1,NumPlantChemElms
        ElmntBalanceCum_pft(NE,NZ)=ShootChemElms_pft(NE,NZ)+RootElmnts_pft(NE,NZ)+NoduleChemElms_pft(NE,NZ) &
          +NonstructalElmnts_pft(NE,NZ)+LitrfallChemElms_pft(NE,NZ)-PlantExudChemElmCum_pft(NE,NZ) &
          -NetCumElmntFlx2Plant_pft(NE,NZ)+StandingDeadChemElms_pft(NE,NZ)+EcoHavstElmnt_pft(NE,NZ)+EcoHavstElmntCum_pft(NE,NZ)
      ENDDO
      ElmntBalanceCum_pft(ielmc,NZ)=ElmntBalanceCum_pft(ielmc,NZ)-NetPrimaryProductvity_pft(NZ)-CO2ByFire_pft(NZ)-CH4ByFire_pft(NZ)
!
!     PLANT N BALANCE = TOTAL N STATE VARIABLES + TOTAL N LITTERFALL
!     - TOTAL N UPTAKE FROM SOIL - TOTAL N ABSORPTION FROM ATMOSPHERE
!
!     BALN=PFT N balance
!     WTSHN,WTRTN,WTNDN,WTRVN,WTSTGN=PFT shoot,root,bacteria,storage,standing dead N
!     TZSNC=cumulative PFT N litterfall
!     TZUPTK=cumulative PFT root-soil N exchange
!     NH3EmiCum_pft=cumulative NH3 exchange
!     RSETN=cumulative N balance from previous year
!     THVSTN=cumulative PFT N removed from ecosystem from previous year
!     HVSTN=total PFT N removed from ecosystem in current year
!     NH3byFire_pft,N2ObyFire_pft=NH3,N2O emission from disturbance
!     PlantN2FixCum_pft=cumulative PFT N2 fixation
!
      ElmntBalanceCum_pft(ielmn,NZ)=ElmntBalanceCum_pft(ielmn,NZ)-NH3EmiCum_pft(NZ) &
        -NH3byFire_pft(NZ)-N2ObyFire_pft(NZ)-PlantN2FixCum_pft(NZ)
!
!     PLANT P BALANCE = TOTAL P STATE VARIABLES + TOTAL P LITTERFALL
!     - TOTAL P UPTAKE FROM SOIL
!
!     BALP=PFT N balance
!     WTSHP,WTRTP,WTNDP,WTRVP,WTSTGP=PFT shoot,root,bacteria,storage,standing dead P
!     TPSNC=cumulative PFT P litterfall
!     TPUPTK=cumulative PFT root-soil P exchange
!     RSETP=cumulative P balance from previous year
!     THVSTP=cumulative PFT P removed from ecosystem from previous year
!     HVSTP=total PFT P removed from ecosystem in current year
!     PO4byFire_pft=PO4 emission from disturbance
!
      ElmntBalanceCum_pft(ielmp,NZ)=ElmntBalanceCum_pft(ielmp,NZ)-PO4byFire_pft(NZ)
    ENDIF
  ENDDO D9975
  end associate
  end subroutine LiveDeadTransformation
!------------------------------------------------------------------------------------------

  subroutine GrowPlant(I,J,NZ,CanopyHeight_copy,ShootNonstructC_brch)
  use PlantDisturbsMod, only : RemoveBiomByManagement
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: CanopyHeight_copy(JP1)
  real(r8), intent(out) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)

  real(r8)  :: CanopyN2Fix_pft(JP1)
  integer  :: ICHK1(2,JZ1),IDTHRN,NB
  integer  :: NRX(2,JZ1),BegRemoblize
  real(r8) :: TFN6(JZ1)
  real(r8) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8) :: PTRT
  real(r8) :: RootAreaPopu
  real(r8) :: TFN5
  real(r8) :: WFNG
  real(r8) :: Stomata_Activity
  real(r8) :: WFNS,WFNSG
! begin_execution
  associate(                              &
    iPlantMorphologyType_pft => plt_pheno%iPlantMorphologyType_pft      , &
    iPlantRootState_pft => plt_pheno%iPlantRootState_pft      , &
    iPlantShootState_pft  => plt_pheno%iPlantShootState_pft       , &
    RootN2Fix_pft   => plt_rbgc%RootN2Fix_pft         , &
    RootH2PO4Uptake_pft  => plt_rbgc%RootH2PO4Uptake_pft        , &
    RootNH4Uptake_pft  => plt_rbgc%RootNH4Uptake_pft        , &
    RootHPO4Uptake_pft  => plt_rbgc%RootHPO4Uptake_pft        , &
    RootNO3Uptake_pft  => plt_rbgc%RootNO3Uptake_pft        , &
    PlantRootSoilChemNetX_pft => plt_rbgc%PlantRootSoilChemNetX_pft       , &
    RootExudChemElmnt_pft  => plt_rbgc%RootExudChemElmnt_pft        , &
    NumOfBranches_pft    => plt_morph%NumOfBranches_pft         , &
    NumRootAxes_pft   => plt_morph%NumRootAxes_pft          &
  )
  IF(iPlantShootState_pft(NZ).EQ.iLive.OR.iPlantRootState_pft(NZ).EQ.iLive)THEN
    CanopyN2Fix_pft(NZ)=0._r8
    BegRemoblize = 0
    call StagePlantForGrowth(I,J,NZ,ICHK1,NRX,TFN6,CNLFW,CPLFW,&
      CNSHW,CPSHW,CNRTW,CPRTW,RootAreaPopu,TFN5,WFNG,Stomata_Activity,WFNS,WFNSG)
!
!     CALCULATE GROWTH OF EACH BRANCH
!
!     WTLFB,WTSHEB,LeafPetolBiomassC_brch=leaf,petiole,leaf+petiole mass
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!
    DO  NB=1,NumOfBranches_pft(NZ)
      call GrowOneBranch(I,J,NB,NZ,TFN6,CanopyHeight_copy,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,&
        TFN5,WFNG,Stomata_Activity,WFNS,WFNSG,PTRT,CanopyN2Fix_pft,BegRemoblize)
    ENDDO
!
    call RootBGCModel(I,J,NZ,BegRemoblize,ICHK1,IDTHRN,NRX,PTRT,TFN6,CNRTW,CPRTW,RootAreaPopu)
!
    call ComputeTotalBiom(NZ,ShootNonstructC_brch)
  ELSE
    PlantRootSoilChemNetX_pft(1:NumPlantChemElms,NZ)=RootExudChemElmnt_pft(1:NumPlantChemElms,NZ)
    PlantRootSoilChemNetX_pft(ielmn,NZ)=PlantRootSoilChemNetX_pft(ielmn,NZ)+RootNH4Uptake_pft(NZ)+RootNO3Uptake_pft(NZ)+RootN2Fix_pft(NZ)
    PlantRootSoilChemNetX_pft(ielmp,NZ)=PlantRootSoilChemNetX_pft(ielmp,NZ)+RootH2PO4Uptake_pft(NZ)+RootHPO4Uptake_pft(NZ)
  ENDIF
!
  call RemoveBiomByManagement(I,J,NZ,ShootNonstructC_brch)
!
!     RESET DEAD BRANCHES
  call ResetDeadBranch(I,J,NZ,ShootNonstructC_brch)
!
  call AccumulateStates(I,J,NZ,CanopyN2Fix_pft)
  end associate
  end subroutine GrowPlant

!------------------------------------------------------------------------------------------
  subroutine StagePlantForGrowth(I,J,NZ,ICHK1,NRX,TFN6,CNLFW,CPLFW,CNSHW,&
    CPSHW,CNRTW,CPRTW,RootAreaPopu,TFN5,WFNG,Stomata_Activity,WFNS,WFNSG)
  integer, intent(in) :: I,J,NZ
  integer, intent(out):: ICHK1(jroots,JZ1)
  integer, intent(out):: NRX(jroots,JZ1)
  REAL(R8), INTENT(OUT):: TFN6(JZ1)
  REAL(R8), INTENT(OUT) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8), intent(out) :: RootAreaPopu,TFN5,WFNG,Stomata_Activity
  real(r8), intent(out) :: WFNS,WFNSG
  integer :: L,NR,N
  real(r8) :: ACTVM,RTK,STK,TKCM,TKSM
!     begin_execution

  associate(                            &
    TKC                     =>  plt_ew%TKC         , &
    TKS                     =>  plt_ew%TKS         , &
    PSICanopy_pft           =>  plt_ew%PSICanopy_pft      , &
    PSICanopyTurg_pft       =>  plt_ew%PSICanopyTurg_pft       , &
    PlantPopulation_pft     =>  plt_site%PlantPopulation_pft        , &
    NU                      =>  plt_site%NU        , &
    MaxNumRootLays          =>  plt_site%MaxNumRootLays        , &
    OFFST                   =>  plt_pheno%OFFST    , &
    ZEROP                   =>  plt_biom%ZEROP     , &
    CanopyStalkC_pft        =>  plt_biom%CanopyStalkC_pft     , &
    RootProteinC_pvr        =>  plt_biom%RootProteinC_pvr     , &
    RootBiomCPerPlant_pft   =>  plt_biom%RootBiomCPerPlant_pft     , &
    StalkChemElms_pft       =>  plt_biom%StalkChemElms_pft    , &
    RootElmnts_pft  =>  plt_biom%RootElmnts_pft     , &
    CanopyLeafCpft_lyr  =>  plt_biom%CanopyLeafCpft_lyr     , &
    iPlantTurnoverPattern_pft =>  plt_pheno%iPlantTurnoverPattern_pft   , &
    iPlantMorphologyType_pft =>  plt_pheno%iPlantMorphologyType_pft   , &
    RCO2A  =>  plt_rbgc%RCO2A     , &
    RootRespPotential_vr  =>  plt_rbgc%RootRespPotential_vr     , &
    RCO2N  =>  plt_rbgc%RCO2N     , &
    RootrNC_pft  =>  plt_allom%RootrNC_pft    , &
    RootrPC_pft  =>  plt_allom%RootrPC_pft    , &
    FracHour4LeafoffRemob   =>  plt_allom%FracHour4LeafoffRemob     , &
    FWODLE =>  plt_allom%FWODLE   , &
    FWODBE =>  plt_allom%FWODBE   , &
    FWOODE =>  plt_allom%FWOODE   , &
    FWODRE =>  plt_allom%FWODRE   , &
    CNLF   =>  plt_allom%CNLF     , &
    CPLF   =>  plt_allom%CPLF     , &
    CNSHE  =>  plt_allom%CNSHE    , &
    CPSHE  =>  plt_allom%CPSHE    , &
    rNCStalk_pft =>  plt_allom%rNCStalk_pft   , &
    rPCStalk_pft  =>  plt_allom%rPCStalk_pft    , &
    k_fine_litr=> pltpar%k_fine_litr,&
    k_woody_litr=> pltpar%k_woody_litr,&
    RCS    =>  plt_photo%RCS      , &
    PrimRootXNumL_pvr  =>  plt_morph%PrimRootXNumL_pvr    , &
    SecndRootXNum_pvr   =>  plt_morph%SecndRootXNum_pvr     , &
    MY                  =>  plt_morph%MY       , &
    CanopyLeafApft_lyr  =>  plt_morph%CanopyLeafApft_lyr    , &
    CanopyStemApft_lyr  =>  plt_morph%CanopyStemApft_lyr    , &
    NumRootAxes_pft     =>  plt_morph%NumRootAxes_pft       &
  )
  D2: DO L=1,NumOfCanopyLayers1
    CanopyLeafApft_lyr(L,NZ)=0._r8
    CanopyLeafCpft_lyr(L,NZ)=0._r8
    CanopyStemApft_lyr(L,NZ)=0._r8
  ENDDO D2
  D5: DO NR=1,NumRootAxes_pft(NZ)
    DO  N=1,MY(NZ)
      NRX(N,NR)=0
      ICHK1(N,NR)=0
    enddo
  ENDDO D5
  D9: DO N=1,MY(NZ)
    D6: DO L=NU,MaxNumRootLays
      RootProteinC_pvr(N,L,NZ)=0._r8
      PrimRootXNumL_pvr(N,L,NZ)=0._r8
      SecndRootXNum_pvr(N,L,NZ)=0._r8
      RootRespPotential_vr(N,L,NZ)=0._r8
      RCO2N(N,L,NZ)=0._r8
      RCO2A(N,L,NZ)=0._r8
    ENDDO D6
  ENDDO D9
!
!     iPlantTurnoverPattern_pft=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     WTSTK,WVSTK=stalk,sapwood mass
!     FWOOD,FWODB=C woody fraction in stalk,other organs:0=woody,1=non-woody
!     CN*,CP*=N:C,P:C ratios in plant organs from PFT files
!     CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood content
!     *LF=leaf,*SHE=petiole,*STK=stalk,*RT=root
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
!     FWOODN,FWOODP=N,P woody fraction in stalk:0=woody,1=non-woody
!
  IF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantMorphologyType_pft(NZ)))&
    .OR.StalkChemElms_pft(ielmc,NZ).LE.ZEROP(NZ))THEN
    FWODBE(ielmc,k_fine_litr)=1.0_r8
    FWOODE(ielmc,k_fine_litr)=1.0_r8
    FWODRE(ielmc,k_fine_litr)=1.0_r8
  ELSE
    FWODBE(ielmc,k_fine_litr)=1.0_r8
    FWOODE(ielmc,k_fine_litr)=SQRT(CanopyStalkC_pft(NZ)/StalkChemElms_pft(ielmc,NZ))
    FWODRE(ielmc,k_fine_litr)=SQRT(FRTX*CanopyStalkC_pft(NZ)/StalkChemElms_pft(ielmc,NZ))
  ENDIF

  FWODBE(ielmc,k_woody_litr)=1.0_r8-FWODBE(ielmc,k_fine_litr)
  FWOODE(ielmc,k_woody_litr)=1.0_r8-FWOODE(ielmc,k_fine_litr)
  FWODRE(ielmc,k_woody_litr)=1.0_r8-FWODRE(ielmc,k_fine_litr)

  CNLFW=FWODBE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FWODBE(ielmc,k_fine_litr)*CNLF(NZ)
  CPLFW=FWODBE(ielmc,k_woody_litr)*rPCStalk_pft(NZ)+FWODBE(ielmc,k_fine_litr)*CPLF(NZ)

  CNSHW=FWODBE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FWODBE(ielmc,k_fine_litr)*CNSHE(NZ)
  CPSHW=FWODBE(ielmc,k_woody_litr)*rPCStalk_pft(NZ)+FWODBE(ielmc,k_fine_litr)*CPSHE(NZ)

  CNRTW=FWODRE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FWODRE(ielmc,k_fine_litr)*RootrNC_pft(NZ)
  CPRTW=FWODRE(ielmc,k_woody_litr)*rPCStalk_pft(NZ)+FWODRE(ielmc,k_fine_litr)*RootrPC_pft(NZ)

  FWODLE(ielmc,1:NumOfPlantLitrCmplxs)=FWODBE(ielmc,1:NumOfPlantLitrCmplxs)

  FWODLE(ielmn,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNLFW
  FWODLE(ielmp,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*rPCStalk_pft(NZ)/CPLFW

  FWODBE(ielmn,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNSHW
  FWODBE(ielmp,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*rPCStalk_pft(NZ)/CPSHW

  FWOODE(ielmn,k_woody_litr)=FWOODE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNRTW
  FWOODE(ielmp,k_woody_litr)=FWOODE(ielmc,k_woody_litr)*rPCStalk_pft(NZ)/CPRTW

  FWODRE(ielmn,k_woody_litr)=FWODRE(ielmc,k_woody_litr)*RootrNC_pft(NZ)/CNRTW
  FWODRE(ielmp,k_woody_litr)=FWODRE(ielmc,k_woody_litr)*RootrPC_pft(NZ)/CPRTW

  FWODLE(ielmn,k_fine_litr)=1.0_r8-FWODLE(ielmn,k_woody_litr)
  FWODLE(ielmp,k_fine_litr)=1.0_r8-FWODLE(ielmp,k_woody_litr)
  FWODBE(ielmn,k_fine_litr)=1.0_r8-FWODBE(ielmn,k_woody_litr)
  FWODBE(ielmp,k_fine_litr)=1.0_r8-FWODBE(ielmp,k_woody_litr)
  FWOODE(ielmn,k_fine_litr)=1.0_r8-FWOODE(ielmn,k_woody_litr)
  FWOODE(ielmp,k_fine_litr)=1.0_r8-FWOODE(ielmp,k_woody_litr)
  FWODRE(ielmn,k_fine_litr)=1.0_r8-FWODRE(ielmn,k_woody_litr)
  FWODRE(ielmp,k_fine_litr)=1.0_r8-FWODRE(ielmp,k_woody_litr)
!
!     SHOOT AND ROOT TEMPERATURE FUNCTIONS FOR MAINTENANCE
!     RESPIRATION FROM TEMPERATURES WITH OFFSETS FOR THERMAL ADAPTATION
!
!     TKC,TKCM=canopy temperature,canopy temp used in Arrhenius eqn
!     TKS,TKSM=soil temperature,soil temp used in Arrhenius eqn
!     OFFST=shift in Arrhenius curve for thermal adaptation
!     TFN5,TFN6=temperature function for canopy,root mntc respn (25 oC =1)
!     8.3143,710.0=gas constant,enthalpy
!     62500,195000,232500=energy of activn,high,low temp inactivn(KJ mol-1)
!
  TKCM=TKC(NZ)+OFFST(NZ)

  TFN5=calc_plant_maint_tempf(TKCM)  
  D7: DO L=NU,MaxNumRootLays
    TKSM=TKS(L)+OFFST(NZ)
    TFN6(L)=calc_plant_maint_tempf(TKSM)  
  ENDDO D7
!
!     PRIMARY ROOT NUMBER
!
!     RootBiomCPerPlant_pft=root mass per plant used to calculate primary root number
!     WTRT,PP=root mass,PFT population
!     RootAreaPopu=multiplier for number of primary root axes
!
  RootBiomCPerPlant_pft(NZ)=AMAX1(0.999992087_r8*RootBiomCPerPlant_pft(NZ),&
    RootElmnts_pft(ielmc,NZ)/PlantPopulation_pft(NZ))
  RootAreaPopu=AMAX1(1.0_r8,RootBiomCPerPlant_pft(NZ)**0.667_r8)*PlantPopulation_pft(NZ)
!
!     WATER STRESS FUNCTIONS FOR EXPANSION AND GROWTH RESPIRATION
!     FROM CANOPY TURGOR
!
!     WFNS=turgor expansion,extension function
!     PSICanopyTurg_pft,PSIMin4OrganExtension=current,minimum canopy turgor potl for expansion,extension
!     Stomata_Activity=stomatal resistance function of canopy turgor
!     PSICanopy_pft=canopy water potential
!     WFNG=growth function of canopy water potential
!     WFNSG=expansion,extension function of canopy water potential
!
  WFNS=AMIN1(1.0_r8,AZMAX1(PSICanopyTurg_pft(NZ)-PSIMin4OrganExtension))

  IF(is_plant_bryophyte(iPlantMorphologyType_pft(NZ)))THEN
    !bryophyte, no turgor
    Stomata_Activity=1.0_r8
    WFNG=EXP(0.05_r8*PSICanopy_pft(NZ))
    WFNSG=WFNS**0.10_r8
  ELSE
    !others
    Stomata_Activity=EXP(RCS(NZ)*PSICanopyTurg_pft(NZ))
    WFNG=EXP(0.10_r8*PSICanopy_pft(NZ))
    WFNSG=WFNS**0.25_r8
  ENDIF
  end associate
  end subroutine StagePlantForGrowth
!------------------------------------------------------------------------------------------

  subroutine ComputeTotalBiom(NZ,ShootNonstructC_brch)

  integer, intent(in) :: NZ
  real(r8), intent(out) :: ShootNonstructC_brch(NumOfCanopyLayers1,JP1)
  integer :: L,K,N,NE,NB
!     begin_execution
  associate(                                 &
    LeafChemElms_brch    =>  plt_biom%LeafChemElms_brch    , &
    GrainChemElms_brch     =>  plt_biom%GrainChemElms_brch     , &
    RootMycoNonstructElmnt_vr     =>  plt_biom%RootMycoNonstructElmnt_vr     , &
    PopuPlantRootC_vr     =>  plt_biom% PopuPlantRootC_vr     , &
    ShootChemElmnt_brch    =>  plt_biom%ShootChemElmnt_brch    , &
    EarChemElms_brch    =>  plt_biom%EarChemElms_brch    , &
    StalkChemElms_brch    =>  plt_biom%StalkChemElms_brch    , &
    PetoleChemElmnt_brch   =>  plt_biom%PetoleChemElmnt_brch   , &
    HuskChemElms_brch    =>  plt_biom%HuskChemElms_brch    , &
    ReserveElmnts_brch    =>  plt_biom%ReserveElmnts_brch    , &
    NonstructElmnt_brch     =>  plt_biom%NonstructElmnt_brch     , &
    NU         =>  plt_site%NU         , &
    CPOOL3     =>  plt_photo%CPOOL3    , &
    CPOOL4     =>  plt_photo%CPOOL4    , &
    CMassHCO3BundleSheath_node       =>  plt_photo%CMassHCO3BundleSheath_node      , &
    CMassCO2BundleSheath_node       =>  plt_photo%CMassCO2BundleSheath_node      , &
    MY         =>  plt_morph%MY        , &
    MaxSoiL4Root         =>  plt_morph%MaxSoiL4Root        , &
    NumOfBranches_pft        =>  plt_morph%NumOfBranches_pft         &
  )
!     TOTAL C,N,P IN EACH BRANCH
!
!     ShootNonstructC_brch=total C4 nonstructural C in branch
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CMassCO2BundleSheath_node,CMassHCO3BundleSheath_node=aqueous CO2,HCO3-C mass in bundle sheath
!     CPOOL,ZPOOL,PPOOL=C3 non-structural C,N,P mass
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     ShootNonstructC_brch=total C4 nonstructural C in branch
!     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
!     ReserveElmnts_brch,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iPlantPhenologyType_pft=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
  DO NE=1,NumPlantChemElms
    DO NB=1,NumOfBranches_pft(NZ)
      ShootChemElmnt_brch(NE,NB,NZ)=LeafChemElms_brch(NE,NB,NZ) &
        +PetoleChemElmnt_brch(NE,NB,NZ)+StalkChemElms_brch(NE,NB,NZ)+ReserveElmnts_brch(NE,NB,NZ) &
        +HuskChemElms_brch(NE,NB,NZ)+EarChemElms_brch(NE,NB,NZ)+GrainChemElms_brch(NE,NB,NZ) &
        +NonstructElmnt_brch(NE,NB,NZ)
    ENDDO
  ENDDO

  D320: DO NB=1,NumOfBranches_pft(NZ)
    ShootNonstructC_brch(NB,NZ)=0._r8
    D325: DO K=1,MaxNodesPerBranch1
      ShootNonstructC_brch(NB,NZ)=ShootNonstructC_brch(NB,NZ)+CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ) &
        +CMassCO2BundleSheath_node(K,NB,NZ)+CMassHCO3BundleSheath_node(K,NB,NZ)
    ENDDO D325
    ShootChemElmnt_brch(ielmc,NB,NZ)=ShootChemElmnt_brch(ielmc,NB,NZ)+ShootNonstructC_brch(NB,NZ)
  ENDDO D320

!
!     TOTAL C,N,P IN ROOTS AND MYCORRHIZAE IN EACH SOIL LAYER
!
!     WTRTD=root C mass
!     CPOOLR=non-structural C mass in root
!     HCUPTK,HZUPTK,HPUPTK=net PFT root-soil C,N,P exchange
!     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
!     RootNH4Uptake_pft,RootNO3Uptake_pft,RootH2PO4Uptake_pft,RootHPO4Uptake_pft=PFT uptake of NH4,NO3,H2PO4,HPO4
!     RootN2Fix_pft=PFT N2 fixation
!
  D345: DO N=1,MY(NZ)
    DO  L=NU,MaxSoiL4Root(NZ)
      PopuPlantRootC_vr(N,L,NZ)= PopuPlantRootC_vr(N,L,NZ)+RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)
    enddo
  ENDDO D345
  end associate
  end subroutine ComputeTotalBiom
!------------------------------------------------------------------------------------------

  subroutine AccumulateStates(I,J,NZ,CanopyN2Fix_pft)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: CanopyN2Fix_pft(JP1)
  integer :: L,NR,N,NE,NB
!     begin_execution
  associate(                            &
    NonstructElmnt_brch   =>  plt_biom%NonstructElmnt_brch  , &
    RootMycoNonstructElmnt_vr   =>  plt_biom%RootMycoNonstructElmnt_vr  , &
    NoduleNonstructElmnt_brch   =>  plt_biom%NoduleNonstructElmnt_brch  , &
    NoduleChemElms_pft    =>  plt_biom%NoduleChemElms_pft   , &
    RootStructElmnt_pft   =>  plt_biom%RootStructElmnt_pft  , &
    RootElmnts_pft    =>  plt_biom%RootElmnts_pft   , &
    CanopyNoduleChemElmnt_brch   =>  plt_biom%CanopyNoduleChemElmnt_brch  , &
    ShootChemElmnt_brch  =>  plt_biom%ShootChemElmnt_brch , &
    StalkChemElms_brch  =>  plt_biom%StalkChemElms_brch , &
    HuskChemElms_brch  =>  plt_biom%HuskChemElms_brch , &
    ReserveElmnts_brch   =>  plt_biom%ReserveElmnts_brch  , &
    EarChemElms_brch  =>  plt_biom%EarChemElms_brch , &
    LeafPetolBiomassC_brch    =>  plt_biom%LeafPetolBiomassC_brch   , &
    GrainChemElms_brch   =>  plt_biom%GrainChemElms_brch  , &
    LeafChemElms_brch  =>  plt_biom%LeafChemElms_brch , &
    StalkBiomassC_brch   =>  plt_biom%StalkBiomassC_brch  , &
    PetoleChemElmnt_brch =>  plt_biom%PetoleChemElmnt_brch, &
    ShootChemElms_pft   =>  plt_biom%ShootChemElms_pft  , &
    LeafChemElms_pft    =>  plt_biom%LeafChemElms_pft   , &
    PetioleChemElms_pft   =>  plt_biom%PetioleChemElms_pft  , &
    StalkChemElms_pft   =>  plt_biom%StalkChemElms_pft  , &
    CanopyStalkC_pft    =>  plt_biom%CanopyStalkC_pft   , &
    ReserveChemElms_pft   =>  plt_biom%ReserveChemElms_pft  , &
    HuskChemElms_pft   =>  plt_biom%HuskChemElms_pft  , &
    EarChemElms_pft   =>  plt_biom%EarChemElms_pft  , &
    GrainChemElms_pft    =>  plt_biom%GrainChemElms_pft   , &
    CanopyLeafShethC_pft     =>  plt_biom%CanopyLeafShethC_pft    , &
    RootNodueChemElmnt_pvr   =>  plt_biom%RootNodueChemElmnt_pvr  , &
    CanopyNonstructElements_pft   =>  plt_biom%CanopyNonstructElements_pft  , &
    NoduleNonstructElmnt_pft   =>  plt_biom%NoduleNonstructElmnt_pft  , &
    Root1stStructChemElmnt_pvr   =>  plt_biom%Root1stStructChemElmnt_pvr  , &
    Root2ndStructChemElmnt_pvr   =>  plt_biom%Root2ndStructChemElmnt_pvr  , &
    RootNoduleNonstructElmnt_vr  =>  plt_biom%RootNoduleNonstructElmnt_vr , &
    PlantN2FixCum_pft   =>  plt_bgcr%PlantN2FixCum_pft  , &
    PlantExudChemElmCum_pft   =>  plt_rbgc%PlantExudChemElmCum_pft  , &
    RootHPO4Uptake_pft    =>  plt_rbgc%RootHPO4Uptake_pft   , &
    PlantRootSoilChemNetX_pft   =>  plt_rbgc%PlantRootSoilChemNetX_pft  , &
    RootN2Fix_pft     =>  plt_rbgc%RootN2Fix_pft    , &
    RootH2PO4Uptake_pft    =>  plt_rbgc%RootH2PO4Uptake_pft   , &
    RootNO3Uptake_pft    =>  plt_rbgc%RootNO3Uptake_pft   , &
    RootNH4Uptake_pft    =>  plt_rbgc%RootNH4Uptake_pft   , &
    RootExudChemElmnt_pft    =>  plt_rbgc%RootExudChemElmnt_pft   , &
    MaxNumRootLays       =>  plt_site%MaxNumRootLays      , &
    NU       =>  plt_site%NU      , &
    NumOfBranches_pft      =>  plt_morph%NumOfBranches_pft    , &
    MY       =>  plt_morph%MY     , &
    MaxSoiL4Root       =>  plt_morph%MaxSoiL4Root     , &
    NumRootAxes_pft     =>  plt_morph%NumRootAxes_pft   , &
    LeafAreaLive_brch    =>  plt_morph%LeafAreaLive_brch  , &
    CanopyStemA_pft    =>  plt_morph%CanopyStemA_pft  , &
    CanopyBranchStemApft_lyr    =>  plt_morph%CanopyBranchStemApft_lyr  , &
    SeedNumberSet_brch    =>  plt_morph%SeedNumberSet_brch  , &
    CanopyLeafArea_pft    =>  plt_morph%CanopyLeafArea_pft  , &
    CanopyStemApft_lyr    =>  plt_morph%CanopyStemApft_lyr  , &
    iPlantNfixType   =>  plt_morph%iPlantNfixType , &
    CanopySeedNumber_pft     =>  plt_morph%CanopySeedNumber_pft     &
  )
!
!     ACCUMULATE PFT STATE VARIABLES FROM BRANCH STATE VARIABLES
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
!     ReserveElmnts_brch,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     LeafAreaLive_brch=branch leaf area
!     CanopyBranchStemApft_lyr=total branch stalk surface area in each layer
!     SeedNumberSet_brch=seed set number
!
  DO NE=1,NumPlantChemElms
    CanopyNonstructElements_pft(NE,NZ)=sum(NonstructElmnt_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    ShootChemElms_pft(NE,NZ)=sum(ShootChemElmnt_brch(NE,1:NumOfBranches_pft(NZ),NZ))

    PetioleChemElms_pft(NE,NZ)=sum(PetoleChemElmnt_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    StalkChemElms_pft(NE,NZ)=sum(StalkChemElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    LeafChemElms_pft(NE,NZ)=sum(LeafChemElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    ReserveChemElms_pft(NE,NZ)=sum(ReserveElmnts_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    HuskChemElms_pft(NE,NZ)=sum(HuskChemElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    GrainChemElms_pft(NE,NZ)=sum(GrainChemElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    EarChemElms_pft(NE,NZ)=sum(EarChemElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    !root state variables
    RootElmnts_pft(NE,NZ)=sum(RootMycoNonstructElmnt_vr(NE,1:MY(NZ),NU:MaxNumRootLays,NZ))
    RootStructElmnt_pft(NE,NZ)=sum(Root1stStructChemElmnt_pvr(NE,1:MY(NZ),NU:MaxNumRootLays,1:NumRootAxes_pft(NZ),NZ)) &
      +sum(Root2ndStructChemElmnt_pvr(NE,1:MY(NZ),NU:MaxNumRootLays,1:NumRootAxes_pft(NZ),NZ))
    RootElmnts_pft(NE,NZ)=RootElmnts_pft(NE,NZ)+RootStructElmnt_pft(NE,NZ)
  ENDDO

  CanopyStalkC_pft(NZ)=sum(StalkBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopyLeafShethC_pft(NZ) =sum(LeafPetolBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopySeedNumber_pft(NZ) =sum(SeedNumberSet_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopyLeafArea_pft(NZ)=sum(LeafAreaLive_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopyStemA_pft(NZ)=sum(CanopyBranchStemApft_lyr(1:NumOfCanopyLayers1,1:NumOfBranches_pft(NZ),NZ))
  CanopyStemApft_lyr(1:NumOfCanopyLayers1,1:NumOfBranches_pft(NZ))=0._r8

  DO NB=1,NumOfBranches_pft(NZ)
    DO L=1,NumOfCanopyLayers1
      CanopyStemApft_lyr(L,NZ)=CanopyStemApft_lyr(L,NZ)+CanopyBranchStemApft_lyr(L,NB,NZ)
    ENDDO
  ENDDO

!
!     ACCUMULATE NODULE STATE VATIABLES FROM NODULE LAYER VARIABLES
!
!     iPlantNfixType=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!
  IF(is_plant_N2fix(iPlantNfixType(NZ)))THEN
    IF(is_canopy_N2fix(iPlantNfixType(NZ)))THEN
      DO NE=1,NumPlantChemElms
        D7950: DO NB=1,NumOfBranches_pft(NZ)
          NoduleNonstructElmnt_pft(NE,NZ)=NoduleNonstructElmnt_pft(NE,NZ)+NoduleNonstructElmnt_brch(NE,NB,NZ)
        ENDDO D7950
        NoduleChemElms_pft(NE,NZ)=sum(CanopyNoduleChemElmnt_brch(NE,1:NumOfBranches_pft(NZ),NZ))+ &
          sum(NoduleNonstructElmnt_brch(NE,1:NumOfBranches_pft(NZ),NZ))
      ENDDO
    ELSEIF(is_root_N2fix(iPlantNfixType(NZ)))THEN
      DO NE=1,NumPlantChemElms
        NoduleChemElms_pft(NE,NZ)=sum(RootNodueChemElmnt_pvr(NE,NU:MaxSoiL4Root(NZ),NZ))+&
          sum(RootNoduleNonstructElmnt_vr(NE,NU:MaxSoiL4Root(NZ),NZ))
      ENDDO
    ENDIF
  ENDIF
!
!     ACCUMULATE TOTAL SOIL-PLANT C,N,P EXCHANGE
!
!     HCUPTK,HZUPTK,HPUPTK=net PFT root-soil C,N,P exchange
!     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
!     RootNH4Uptake_pft,RootNO3Uptake_pft,RootH2PO4Uptake_pft,RootHPO4Uptake_pft=PFT uptake of NH4,NO3,H2PO4,HPO4
!     RootN2Fix_pft=PFT N2 fixation
!     TCUPTK,TZUPTK,TPUPTK=cumulative PFT root-soil C,N,P exchange
!     PlantN2FixCum_pft=cumulative PFT N2 fixation
!
  PlantRootSoilChemNetX_pft(1:NumPlantChemElms,NZ)=RootExudChemElmnt_pft(1:NumPlantChemElms,NZ)
  PlantRootSoilChemNetX_pft(ielmn,NZ)=PlantRootSoilChemNetX_pft(ielmn,NZ)+&
    RootNH4Uptake_pft(NZ)+RootNO3Uptake_pft(NZ)+RootN2Fix_pft(NZ)
  PlantRootSoilChemNetX_pft(ielmp,NZ)=PlantRootSoilChemNetX_pft(ielmp,NZ)+RootH2PO4Uptake_pft(NZ)+RootHPO4Uptake_pft(NZ)

  PlantExudChemElmCum_pft(1:NumPlantChemElms,NZ)=PlantExudChemElmCum_pft(1:NumPlantChemElms,NZ)+&
    RootExudChemElmnt_pft(1:NumPlantChemElms,NZ)
  PlantExudChemElmCum_pft(ielmn,NZ)=PlantExudChemElmCum_pft(ielmn,NZ)+ &
    RootNH4Uptake_pft(NZ)+RootNO3Uptake_pft(NZ)
  PlantExudChemElmCum_pft(ielmp,NZ)=PlantExudChemElmCum_pft(ielmp,NZ)+ &
    RootH2PO4Uptake_pft(NZ)+RootHPO4Uptake_pft(NZ)
  PlantN2FixCum_pft(NZ)=PlantN2FixCum_pft(NZ)+RootN2Fix_pft(NZ)+CanopyN2Fix_pft(NZ)
  end associate
  end subroutine AccumulateStates

end module grosubsMod
