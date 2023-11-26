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
!     FLG4Y=number of hours after physiol maturity required for senescence
!     ATRPX=number of hours required to initiate remobilization of storage C for leafout
!     GVMX=specific oxidation rate of nonstructural C during leafout at 25 C
!     FVRN=fraction of hours required for leafoff to initiate remobilization
!

  integer,private :: curday,curhour
  public :: grosubs
  public :: InitGrosub
  contains

  subroutine InitGrosub(NumGrowthStages,JRS)

  implicit none
  integer, intent(out) :: NumGrowthStages,JRS

  call InitVegPars(pltpar)
  NumGrowthStages = pltpar%NumGrowthStages
  jrs = pltpar%JRS


  end subroutine InitGrosub
!------------------------------------------------------------------------------------------

  subroutine grosubs(I,J)
!
!     THIS subroutine CALCULATES ALL PLANT BIOLOGICAL TRANSFORMATIONS
!
  use PlantDisturbsMod, only : RemoveBiomassByDisturbance
  implicit none
  integer, intent(in) :: I, J

  real(r8) :: ZCX(JP1)
  integer :: L,K,M
  integer :: NZ,NE
  real(r8) :: CPOOLK(NumOfCanopyLayers1,JP1)
! begin_execution
  associate(                            &
    IsPlantActive    => plt_pheno%IsPlantActive   , &
    NP       => plt_site%NP       , &
    NP0      => plt_site%NP0      , &
    NJ       => plt_site%NJ       , &
    CO2NetFix_pft     => plt_bgcr%CO2NetFix_pft     , &
    HESNC    => plt_bgcr%HESNC    , &
    ESNC     => plt_bgcr%ESNC     , &
    CanopyHeight       => plt_morph%CanopyHeight        &
  )
!     TOTAL AGB FOR GRAZING IN LANDSCAPE SECTION
!
!
!     INITIALIZE SENESCENCE ARRAYS
!

  D9980: DO NZ=1,NP0
    D1: DO L=0,NJ
      DO K=1,pltpar%NumOfPlantLitrCmplxs
        DO M=1,jsken
          DO NE=1,NumOfPlantChemElmnts
            ESNC(NE,M,K,L,NZ)=0._r8
          ENDDO
        ENDDO
      ENDDO
    ENDDO D1
    HESNC(1:NumOfPlantChemElmnts,NZ)=0._r8
    CO2NetFix_pft(NZ)=0._r8
    ZCX(NZ)=CanopyHeight(NZ)
    CanopyHeight(NZ)=0._r8
  ENDDO D9980
!
!     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
!
  D9985: DO NZ=1,NP

! IsPlantActive= flag for living pft
    IF(IsPlantActive(NZ).EQ.iPlantIsActive)THEN
      call GrowPlant(I,J,NZ,ZCX,CPOOLK)
    ENDIF

!   HARVEST STANDING DEAD
    call RemoveBiomassByDisturbance(I,J,NZ,CPOOLK)
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
    iDayPlanting   => plt_distb%iDayPlanting   , &
    iYearPlanting    => plt_distb%iYearPlanting    , &
    THVSTE  => plt_distb%THVSTE  , &
    HVSTE   => plt_distb%HVSTE   , &
    VPO4F   => plt_distb%VPO4F   , &
    VN2OF   => plt_distb%VN2OF   , &
    CH4ByFire_pft   => plt_distb%CH4ByFire_pft   , &
    CO2ByFire_pft   => plt_distb%CO2ByFire_pft   , &
    VNH3F   => plt_distb%VNH3F   , &
    StandingDeadKCompChemElmnts_pft  => plt_biom%StandingDeadKCompChemElmnts_pft   , &
    RootChemElmnts_pft   => plt_biom%RootChemElmnts_pft    , &
    ShootChemElmnts_pft  => plt_biom%ShootChemElmnts_pft  , &
    NoduleChemElmnts_pft   => plt_biom%NoduleChemElmnts_pft    , &
    NonstructalChemElmnts_pft   => plt_biom%NonstructalChemElmnts_pft    , &
    StandingDeadChemElmnts_pft  => plt_biom%StandingDeadChemElmnts_pft   , &
    fTgrowCanP    => plt_pheno%fTgrowCanP    , &
    RSETE   => plt_pheno%RSETE   , &
    IsPlantActive   => plt_pheno%IsPlantActive   , &
    iPlantMorphologyType  => plt_pheno%iPlantMorphologyType  , &
    doInitPlant   => plt_pheno%doInitPlant   , &
    doPlantLeafOut   => plt_pheno%doPlantLeafOut   , &
    iPlantTurnoverPattern  => plt_pheno%iPlantTurnoverPattern  , &
    HourThreshold4LeafOut   => plt_pheno%HourThreshold4LeafOut   , &
    Hours4Leafout    => plt_pheno%Hours4Leafout    , &
    BALE    => plt_site%BALE     , &
    NP0     => plt_site%NP0      , &
    NJ      => plt_site%NJ       , &
    iYearCurrent    => plt_site%iYearCurrent     , &
    ESNC    => plt_bgcr%ESNC     , &
    NetPrimaryProductvity_pft    => plt_bgcr%NetPrimaryProductvity_pft     , &
    TZUPFX  => plt_bgcr%TZUPFX   , &
    HESNC   => plt_bgcr%HESNC    , &
    TNH3C   => plt_bgcr%TNH3C    , &
    LitrfallChemElmnts_pft   => plt_bgcr%LitrfallChemElmnts_pft    , &
    SurfLitrfallChemElmnts_pft   => plt_bgcr%SurfLitrfallChemElmnts_pft    , &
    GrossResp_pft   => plt_bgcr%GrossResp_pft    , &
    GrossCO2Fix_pft   => plt_bgcr%GrossCO2Fix_pft    , &
    PlantExudChemElmnts_pft  => plt_rbgc%PlantExudChemElmnts_pft   , &
    CumSoilThickness  => plt_site%CumSoilThickness   , &
    PlantinDepth  => plt_morph%PlantinDepth  , &
    NumOfBranches_pft     => plt_morph%NumOfBranches_pft       &
  )
  D9975: DO NZ=1,NP0
!
!     ACTIVATE DORMANT SEEDS
!
    D205: DO NB=1,NumOfBranches_pft(NZ)
      IF(doInitPlant(NZ).EQ.itrue)THEN
        IF(doPlantLeafOut(NB,NZ).EQ.iEnable.AND.Hours4Leafout(NB,NZ).GE.HourThreshold4LeafOut(NB,NZ))THEN
          iDayPlanting(NZ)=I
          iYearPlanting(NZ)=iYearCurrent
          PlantinDepth(NZ)=0.005_r8+CumSoilThickness(0)
          !mark plant as initialized
          doInitPlant(NZ)=ifalse
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
    DO NE=1,NumOfPlantChemElmnts
      D6235: DO M=1,jsken
        XFRE=1.5814E-05_r8*fTgrowCanP(NZ)*StandingDeadKCompChemElmnts_pft(NE,M,NZ)
        IF(iPlantTurnoverPattern(NZ).EQ.0.OR.iPlantMorphologyType(NZ).LE.1)THEN
          ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+XFRE
        ELSE
          ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ)+XFRE
        ENDIF
        StandingDeadKCompChemElmnts_pft(NE,M,NZ)=StandingDeadKCompChemElmnts_pft(NE,M,NZ)-XFRE
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
      DO NE=1,NumOfPlantChemElmnts
        D6430: DO M=1,jsken
          SurfLitrfallChemElmnts_pft(NE,NZ)=SurfLitrfallChemElmnts_pft(NE,NZ)+ESNC(NE,M,K,0,NZ)
          D8955: DO L=0,NJ
            HESNC(NE,NZ)=HESNC(NE,NZ)+ESNC(NE,M,K,L,NZ)
            LitrfallChemElmnts_pft(NE,NZ)=LitrfallChemElmnts_pft(NE,NZ)+ESNC(NE,M,K,L,NZ)
          ENDDO D8955
        ENDDO D6430
      enddo
    ENDDO
!
!     TOTAL STANDING DEAD
!
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
    DO NE=1,NumOfPlantChemElmnts
      StandingDeadChemElmnts_pft(NE,NZ)=sum(StandingDeadKCompChemElmnts_pft(NE,1:jsken,NZ))
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

    IF(IsPlantActive(NZ).EQ.iPlantIsActive)THEN
      !check for living plant
      DO NE=1,NumOfPlantChemElmnts
        BALE(NE,NZ)=ShootChemElmnts_pft(NE,NZ)+RootChemElmnts_pft(NE,NZ)+NoduleChemElmnts_pft(NE,NZ) &
          +NonstructalChemElmnts_pft(NE,NZ)+LitrfallChemElmnts_pft(NE,NZ)-PlantExudChemElmnts_pft(NE,NZ) &
          -RSETE(NE,NZ)+StandingDeadChemElmnts_pft(NE,NZ)+HVSTE(NE,NZ)+THVSTE(NE,NZ)
      ENDDO
      BALE(ielmc,NZ)=BALE(ielmc,NZ)-NetPrimaryProductvity_pft(NZ)-CO2ByFire_pft(NZ)-CH4ByFire_pft(NZ)
!
!     PLANT N BALANCE = TOTAL N STATE VARIABLES + TOTAL N LITTERFALL
!     - TOTAL N UPTAKE FROM SOIL - TOTAL N ABSORPTION FROM ATMOSPHERE
!
!     BALN=PFT N balance
!     WTSHN,WTRTN,WTNDN,WTRVN,WTSTGN=PFT shoot,root,bacteria,storage,standing dead N
!     TZSNC=cumulative PFT N litterfall
!     TZUPTK=cumulative PFT root-soil N exchange
!     TNH3C=cumulative NH3 exchange
!     RSETN=cumulative N balance from previous year
!     THVSTN=cumulative PFT N removed from ecosystem from previous year
!     HVSTN=total PFT N removed from ecosystem in current year
!     VNH3F,VN2OF=NH3,N2O emission from disturbance
!     TZUPFX=cumulative PFT N2 fixation
!
      BALE(ielmn,NZ)=BALE(ielmn,NZ)-TNH3C(NZ)-VNH3F(NZ)-VN2OF(NZ)-TZUPFX(NZ)
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
!     VPO4F=PO4 emission from disturbance
!
      BALE(ielmp,NZ)=BALE(ielmp,NZ)-VPO4F(NZ)
    ENDIF
  ENDDO D9975
  end associate
  end subroutine LiveDeadTransformation
!------------------------------------------------------------------------------------------

  subroutine GrowPlant(I,J,NZ,ZCX,CPOOLK)
  use PlantDisturbsMod, only : RemoveBiomByManagement
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: ZCX(JP1)
  real(r8), intent(out) :: CPOOLK(NumOfCanopyLayers1,JP1)

  real(r8)  :: UPNFC(JP1)
  integer  :: ICHK1(2,JZ1),IDTHRN,NB
  integer  :: NRX(2,JZ1),IFLGZ
  real(r8) :: TFN6(JZ1)
  real(r8) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8) :: PTRT
  real(r8) :: XRTN1
  real(r8) :: TFN5
  real(r8) :: WFNG
  real(r8) :: Stomata_Activity
  real(r8) :: WFNS,WFNSG
! begin_execution
  associate(                              &
    iPlantMorphologyType => plt_pheno%iPlantMorphologyType      , &
    iPlantRootState => plt_pheno%iPlantRootState      , &
    iPlantShootState  => plt_pheno%iPlantShootState       , &
    UPNF   => plt_rbgc%UPNF         , &
    UPH2P  => plt_rbgc%UPH2P        , &
    UPNH4  => plt_rbgc%UPNH4        , &
    UPH1P  => plt_rbgc%UPH1P        , &
    UPNO3  => plt_rbgc%UPNO3        , &
    HEUPTK => plt_rbgc%HEUPTK       , &
    UPOME  => plt_rbgc%UPOME        , &
    NumOfBranches_pft    => plt_morph%NumOfBranches_pft         , &
    NumRootAxes_pft   => plt_morph%NumRootAxes_pft          &
  )
  IF(iPlantShootState(NZ).EQ.0.OR.iPlantRootState(NZ).EQ.iLive)THEN
    UPNFC(NZ)=0._r8
    IFLGZ = 0
    call StagePlantForGrowth(I,J,NZ,ICHK1,NRX,TFN6,CNLFW,CPLFW,&
      CNSHW,CPSHW,CNRTW,CPRTW,XRTN1,TFN5,WFNG,Stomata_Activity,WFNS,WFNSG)
!
!     CALCULATE GROWTH OF EACH BRANCH
!
!     WTLFB,WTSHEB,LeafPetioleBiomassC_brch=leaf,petiole,leaf+petiole mass
!     iPlantBranchState=branch living flag: 0=alive,1=dead
!
    DO  NB=1,NumOfBranches_pft(NZ)
      call GrowOneBranch(I,J,NB,NZ,TFN6,ZCX,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,&
        TFN5,WFNG,Stomata_Activity,WFNS,WFNSG,PTRT,UPNFC,IFLGZ)
    ENDDO
!
    call RootBGCModel(I,J,NZ,IFLGZ,ICHK1,IDTHRN,NRX,PTRT,TFN6,CNRTW,CPRTW,XRTN1)
!
    call ComputeTotalBiom(NZ,CPOOLK)
  ELSE
    HEUPTK(1:NumOfPlantChemElmnts,NZ)=UPOME(1:NumOfPlantChemElmnts,NZ)
    HEUPTK(ielmn,NZ)=HEUPTK(ielmn,NZ)+UPNH4(NZ)+UPNO3(NZ)+UPNF(NZ)
    HEUPTK(ielmp,NZ)=HEUPTK(ielmp,NZ)+UPH2P(NZ)+UPH1P(NZ)
  ENDIF
!
  call RemoveBiomByManagement(I,J,NZ,CPOOLK)
!
!     RESET DEAD BRANCHES
  call ResetDeadBranch(I,J,NZ,CPOOLK)
!
  call AccumulateStates(I,J,NZ,UPNFC)
  end associate
  end subroutine GrowPlant

!------------------------------------------------------------------------------------------
  subroutine StagePlantForGrowth(I,J,NZ,ICHK1,NRX,TFN6,CNLFW,CPLFW,CNSHW,&
    CPSHW,CNRTW,CPRTW,XRTN1,TFN5,WFNG,Stomata_Activity,WFNS,WFNSG)
  integer, intent(in) :: I,J,NZ
  integer, intent(out):: ICHK1(2,JZ1)
  integer, intent(out):: NRX(2,JZ1)
  REAL(R8), INTENT(OUT):: TFN6(JZ1)
  REAL(R8), INTENT(OUT) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,XRTN1,TFN5,WFNG,Stomata_Activity
  real(r8), intent(out) :: WFNS,WFNSG
  integer :: L,NR,N
  real(r8) :: ACTVM,RTK,STK,TKCM,TKSM
!     begin_execution

  associate(                            &
    TKC    =>  plt_ew%TKC         , &
    TKS    =>  plt_ew%TKS         , &
    PSICanP  =>  plt_ew%PSICanP       , &
    PSICanPTurg  =>  plt_ew%PSICanPTurg       , &
    pftPlantPopulation     =>  plt_site%pftPlantPopulation        , &
    NU     =>  plt_site%NU        , &
    NJ     =>  plt_site%NJ        , &
    OFFST  =>  plt_pheno%OFFST    , &
    ZEROP  =>  plt_biom%ZEROP     , &
    CanPStalkC  =>  plt_biom%CanPStalkC     , &
    WSRTL  =>  plt_biom%WSRTL     , &
    WTRTA  =>  plt_biom%WTRTA     , &
    StalkChemElmnts =>  plt_biom%StalkChemElmnts    , &
    RootChemElmnts_pft  =>  plt_biom%RootChemElmnts_pft     , &
    CanopyLeafCpft_lyr  =>  plt_biom%CanopyLeafCpft_lyr     , &
    iPlantTurnoverPattern =>  plt_pheno%iPlantTurnoverPattern   , &
    iPlantMorphologyType =>  plt_pheno%iPlantMorphologyType   , &
    RCO2A  =>  plt_rbgc%RCO2A     , &
    RCO2M  =>  plt_rbgc%RCO2M     , &
    RCO2N  =>  plt_rbgc%RCO2N     , &
    RootrNC_pft  =>  plt_allom%RootrNC_pft    , &
    CPRT   =>  plt_allom%CPRT     , &
    FVRN   =>  plt_allom%FVRN     , &
    FWODLE =>  plt_allom%FWODLE   , &
    FWODBE =>  plt_allom%FWODBE   , &
    FWOODE =>  plt_allom%FWOODE   , &
    FWODRE =>  plt_allom%FWODRE   , &
    CNLF   =>  plt_allom%CNLF     , &
    CPLF   =>  plt_allom%CPLF     , &
    CNSHE  =>  plt_allom%CNSHE    , &
    CPSHE  =>  plt_allom%CPSHE    , &
    rNCStalk_pft =>  plt_allom%rNCStalk_pft   , &
    CPSTK  =>  plt_allom%CPSTK    , &
    k_fine_litr=> pltpar%k_fine_litr,&
    k_woody_litr=> pltpar%k_woody_litr,&
    RCS    =>  plt_photo%RCS      , &
    PrimRootXNumL  =>  plt_morph%PrimRootXNumL    , &
    SecndRootXNumL   =>  plt_morph%SecndRootXNumL     , &
    MY     =>  plt_morph%MY       , &
    CanopyLeafApft_lyr  =>  plt_morph%CanopyLeafApft_lyr    , &
    CanopyStemApft_lyr  =>  plt_morph%CanopyStemApft_lyr    , &
    NumRootAxes_pft   =>  plt_morph%NumRootAxes_pft       &
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
    D6: DO L=NU,NJ
      WSRTL(N,L,NZ)=0._r8
      PrimRootXNumL(N,L,NZ)=0._r8
      SecndRootXNumL(N,L,NZ)=0._r8
      RCO2M(N,L,NZ)=0._r8
      RCO2N(N,L,NZ)=0._r8
      RCO2A(N,L,NZ)=0._r8
    ENDDO D6
  ENDDO D9
!
!     iPlantTurnoverPattern=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     WTSTK,WVSTK=stalk,sapwood mass
!     FWOOD,FWODB=C woody fraction in stalk,other organs:0=woody,1=non-woody
!     CN*,CP*=N:C,P:C ratios in plant organs from PFT files
!     CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood content
!     *LF=leaf,*SHE=petiole,*STK=stalk,*RT=root
!     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
!     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
!     FWOODN,FWOODP=N,P woody fraction in stalk:0=woody,1=non-woody
!
  IF(iPlantTurnoverPattern(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType(NZ)))&
    .OR.StalkChemElmnts(ielmc,NZ).LE.ZEROP(NZ))THEN
    FWODBE(ielmc,k_fine_litr)=1.0_r8
    FWOODE(ielmc,k_fine_litr)=1.0_r8
    FWODRE(ielmc,k_fine_litr)=1.0_r8
  ELSE
    FWODBE(ielmc,k_fine_litr)=1.0_r8
    FWOODE(ielmc,k_fine_litr)=SQRT(CanPStalkC(NZ)/StalkChemElmnts(ielmc,NZ))
    FWODRE(ielmc,k_fine_litr)=SQRT(FRTX*CanPStalkC(NZ)/StalkChemElmnts(ielmc,NZ))
  ENDIF

  FWODBE(ielmc,k_woody_litr)=1.0_r8-FWODBE(ielmc,k_fine_litr)
  FWOODE(ielmc,k_woody_litr)=1.0_r8-FWOODE(ielmc,k_fine_litr)
  FWODRE(ielmc,k_woody_litr)=1.0_r8-FWODRE(ielmc,k_fine_litr)

  CNLFW=FWODBE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FWODBE(ielmc,k_fine_litr)*CNLF(NZ)
  CPLFW=FWODBE(ielmc,k_woody_litr)*CPSTK(NZ)+FWODBE(ielmc,k_fine_litr)*CPLF(NZ)

  CNSHW=FWODBE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FWODBE(ielmc,k_fine_litr)*CNSHE(NZ)
  CPSHW=FWODBE(ielmc,k_woody_litr)*CPSTK(NZ)+FWODBE(ielmc,k_fine_litr)*CPSHE(NZ)

  CNRTW=FWODRE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)+FWODRE(ielmc,k_fine_litr)*RootrNC_pft(NZ)
  CPRTW=FWODRE(ielmc,k_woody_litr)*CPSTK(NZ)+FWODRE(ielmc,k_fine_litr)*CPRT(NZ)

  FWODLE(ielmc,1:NumOfPlantLitrCmplxs)=FWODBE(ielmc,1:NumOfPlantLitrCmplxs)

  FWODLE(ielmn,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNLFW
  FWODLE(ielmp,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*CPSTK(NZ)/CPLFW

  FWODBE(ielmn,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNSHW
  FWODBE(ielmp,k_woody_litr)=FWODBE(ielmc,k_woody_litr)*CPSTK(NZ)/CPSHW

  FWOODE(ielmn,k_woody_litr)=FWOODE(ielmc,k_woody_litr)*rNCStalk_pft(NZ)/CNRTW
  FWOODE(ielmp,k_woody_litr)=FWOODE(ielmc,k_woody_litr)*CPSTK(NZ)/CPRTW

  FWODRE(ielmn,k_woody_litr)=FWODRE(ielmc,k_woody_litr)*RootrNC_pft(NZ)/CNRTW
  FWODRE(ielmp,k_woody_litr)=FWODRE(ielmc,k_woody_litr)*CPRT(NZ)/CPRTW

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
  RTK=RGAS*TKCM
  STK=710.0_r8*TKCM
  ACTVM=1._r8+EXP((195000._r8-STK)/RTK)+EXP((STK-232500._r8)/RTK)
  TFN5=EXP(25.214_r8-62500._r8/RTK)/ACTVM
  D7: DO L=NU,NJ
    TKSM=TKS(L)+OFFST(NZ)
    RTK=RGAS*TKSM
    STK=710.0_r8*TKSM
    ACTVM=1+EXP((195000._r8-STK)/RTK)+EXP((STK-232500._r8)/RTK)
    TFN6(L)=EXP(25.214_r8-62500._r8/RTK)/ACTVM
  ENDDO D7
!
!     PRIMARY ROOT NUMBER
!
!     WTRTA=root mass per plant used to calculate primary root number
!     WTRT,PP=root mass,PFT population
!     XRTN1=multiplier for number of primary root axes
!
  WTRTA(NZ)=AMAX1(0.999992087_r8*WTRTA(NZ),RootChemElmnts_pft(ielmc,NZ)/pftPlantPopulation(NZ))
  XRTN1=AMAX1(1.0_r8,WTRTA(NZ)**0.667_r8)*pftPlantPopulation(NZ)
!
!     WATER STRESS FUNCTIONS FOR EXPANSION AND GROWTH RESPIRATION
!     FROM CANOPY TURGOR
!
!     WFNS=turgor expansion,extension function
!     PSICanPTurg,PSILM=current,minimum canopy turgor potl for expansion,extension
!     Stomata_Activity=stomatal resistance function of canopy turgor
!     PSICanP=canopy water potential
!     WFNG=growth function of canopy water potential
!     WFNSG=expansion,extension function of canopy water potential
!
  WFNS=AMIN1(1.0_r8,AZMAX1(PSICanPTurg(NZ)-PSILM))

  IF(is_plant_bryophyte(iPlantMorphologyType(NZ)))THEN
    !bryophyte, no turgor
    Stomata_Activity=1.0_r8
    WFNG=EXP(0.05_r8*PSICanP(NZ))
    WFNSG=WFNS**0.10_r8
  ELSE
    !others
    Stomata_Activity=EXP(RCS(NZ)*PSICanPTurg(NZ))
    WFNG=EXP(0.10_r8*PSICanP(NZ))
    WFNSG=WFNS**0.25_r8
  ENDIF
  end associate
  end subroutine StagePlantForGrowth
!------------------------------------------------------------------------------------------

  subroutine ComputeTotalBiom(NZ,CPOOLK)

  integer, intent(in) :: NZ
  real(r8), intent(out) :: CPOOLK(NumOfCanopyLayers1,JP1)
  integer :: L,K,N,NE,NB
!     begin_execution
  associate(                                 &
    LeafChemElmnts_brch    =>  plt_biom%LeafChemElmnts_brch    , &
    GrainChemElmnts_brch     =>  plt_biom%GrainChemElmnts_brch     , &
     RootMycoNonstructElmnt_vr     =>  plt_biom% RootMycoNonstructElmnt_vr     , &
     PopuPlantRootC_vr     =>  plt_biom% PopuPlantRootC_vr     , &
    WTSHTBE    =>  plt_biom%WTSHTBE    , &
    EarChemElmnts_brch    =>  plt_biom%EarChemElmnts_brch    , &
    StalkChemElmnts_brch    =>  plt_biom%StalkChemElmnts_brch    , &
    PetioleChemElmnts_brch   =>  plt_biom%PetioleChemElmnts_brch   , &
    HuskChemElmnts_brch    =>  plt_biom%HuskChemElmnts_brch    , &
    ReserveChemElmnts_brch    =>  plt_biom%ReserveChemElmnts_brch    , &
    NonstructElmnt_brch     =>  plt_biom%NonstructElmnt_brch     , &
    NU         =>  plt_site%NU         , &
    CPOOL3     =>  plt_photo%CPOOL3    , &
    CPOOL4     =>  plt_photo%CPOOL4    , &
    HCOB       =>  plt_photo%HCOB      , &
    CO2B       =>  plt_photo%CO2B      , &
    MY         =>  plt_morph%MY        , &
    NI         =>  plt_morph%NI        , &
    NumOfBranches_pft        =>  plt_morph%NumOfBranches_pft         &
  )
!     TOTAL C,N,P IN EACH BRANCH
!
!     CPOOLK=total C4 nonstructural C in branch
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     CPOOL,ZPOOL,PPOOL=C3 non-structural C,N,P mass
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CPOOLK=total C4 nonstructural C in branch
!     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
!     ReserveChemElmnts_brch,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     FWODB=C woody fraction in other organs:0=woody,1=non-woody
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     iPlantPhenologyPattern=growth habit:0=annual,1=perennial from PFT file
!     iPlantPhenologyType=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
  DO NE=1,NumOfPlantChemElmnts
    DO NB=1,NumOfBranches_pft(NZ)
      WTSHTBE(NE,NB,NZ)=LeafChemElmnts_brch(NE,NB,NZ) &
        +PetioleChemElmnts_brch(NE,NB,NZ)+StalkChemElmnts_brch(NE,NB,NZ)+ReserveChemElmnts_brch(NE,NB,NZ) &
        +HuskChemElmnts_brch(NE,NB,NZ)+EarChemElmnts_brch(NE,NB,NZ)+GrainChemElmnts_brch(NE,NB,NZ) &
        +NonstructElmnt_brch(NE,NB,NZ)
    ENDDO
  ENDDO

  D320: DO NB=1,NumOfBranches_pft(NZ)
    CPOOLK(NB,NZ)=0._r8
    D325: DO K=1,MaxNodesPerBranch1
      CPOOLK(NB,NZ)=CPOOLK(NB,NZ)+CPOOL3(K,NB,NZ)+CPOOL4(K,NB,NZ) &
        +CO2B(K,NB,NZ)+HCOB(K,NB,NZ)
    ENDDO D325
    WTSHTBE(ielmc,NB,NZ)=WTSHTBE(ielmc,NB,NZ)+CPOOLK(NB,NZ)
  ENDDO D320

!
!     TOTAL C,N,P IN ROOTS AND MYCORRHIZAE IN EACH SOIL LAYER
!
!     WTRTD=root C mass
!     CPOOLR=non-structural C mass in root
!     HCUPTK,HZUPTK,HPUPTK=net PFT root-soil C,N,P exchange
!     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
!     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
!     UPNF=PFT N2 fixation
!
  D345: DO N=1,MY(NZ)
    DO  L=NU,NI(NZ)
       PopuPlantRootC_vr(N,L,NZ)= PopuPlantRootC_vr(N,L,NZ)+ RootMycoNonstructElmnt_vr(ielmc,N,L,NZ)
    enddo
  ENDDO D345
  end associate
  end subroutine ComputeTotalBiom
!------------------------------------------------------------------------------------------

  subroutine AccumulateStates(I,J,NZ,UPNFC)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: UPNFC(JP1)
  integer :: L,NR,N,NE,NB
!     begin_execution
  associate(                            &
    NonstructElmnt_brch   =>  plt_biom%NonstructElmnt_brch  , &
     RootMycoNonstructElmnt_vr   =>  plt_biom% RootMycoNonstructElmnt_vr  , &
    NoduleNonstructElmnt_brch   =>  plt_biom%NoduleNonstructElmnt_brch  , &
    NoduleChemElmnts_pft    =>  plt_biom%NoduleChemElmnts_pft   , &
    WTRTSE   =>  plt_biom%WTRTSE  , &
    RootChemElmnts_pft    =>  plt_biom%RootChemElmnts_pft   , &
    WTNDBE   =>  plt_biom%WTNDBE  , &
    WTSHTBE  =>  plt_biom%WTSHTBE , &
    StalkChemElmnts_brch  =>  plt_biom%StalkChemElmnts_brch , &
    HuskChemElmnts_brch  =>  plt_biom%HuskChemElmnts_brch , &
    ReserveChemElmnts_brch   =>  plt_biom%ReserveChemElmnts_brch  , &
    EarChemElmnts_brch  =>  plt_biom%EarChemElmnts_brch , &
    LeafPetioleBiomassC_brch    =>  plt_biom%LeafPetioleBiomassC_brch   , &
    GrainChemElmnts_brch   =>  plt_biom%GrainChemElmnts_brch  , &
    LeafChemElmnts_brch  =>  plt_biom%LeafChemElmnts_brch , &
    StalkBiomassC_brch   =>  plt_biom%StalkBiomassC_brch  , &
    PetioleChemElmnts_brch =>  plt_biom%PetioleChemElmnts_brch, &
    ShootChemElmnts_pft   =>  plt_biom%ShootChemElmnts_pft  , &
    LeafChemElmnts    =>  plt_biom%LeafChemElmnts   , &
    SheathChemElmnts   =>  plt_biom%SheathChemElmnts  , &
    StalkChemElmnts   =>  plt_biom%StalkChemElmnts  , &
    CanPStalkC    =>  plt_biom%CanPStalkC   , &
    CanopyReserveChemElmnts   =>  plt_biom%CanopyReserveChemElmnts  , &
    HuskChemElmnts   =>  plt_biom%HuskChemElmnts  , &
    EarChemElmnts   =>  plt_biom%EarChemElmnts  , &
    GrainChemElmnts    =>  plt_biom%GrainChemElmnts   , &
    CanopyLeafShethC_pft     =>  plt_biom%CanopyLeafShethC_pft    , &
    WTNDLE   =>  plt_biom%WTNDLE  , &
    CanopyNonstructElements_pft   =>  plt_biom%CanopyNonstructElements_pft  , &
    NoduleNonstructElmnt_pft   =>  plt_biom%NoduleNonstructElmnt_pft  , &
    WTRT1E   =>  plt_biom%WTRT1E  , &
    WTRT2E   =>  plt_biom%WTRT2E  , &
    RootNoduleNonstructElmnt_vr  =>  plt_biom%RootNoduleNonstructElmnt_vr , &
    TZUPFX   =>  plt_bgcr%TZUPFX  , &
    PlantExudChemElmnts_pft   =>  plt_rbgc%PlantExudChemElmnts_pft  , &
    UPH1P    =>  plt_rbgc%UPH1P   , &
    HEUPTK   =>  plt_rbgc%HEUPTK  , &
    UPNF     =>  plt_rbgc%UPNF    , &
    UPH2P    =>  plt_rbgc%UPH2P   , &
    UPNO3    =>  plt_rbgc%UPNO3   , &
    UPNH4    =>  plt_rbgc%UPNH4   , &
    UPOME    =>  plt_rbgc%UPOME   , &
    NJ       =>  plt_site%NJ      , &
    NU       =>  plt_site%NU      , &
    NumOfBranches_pft      =>  plt_morph%NumOfBranches_pft    , &
    MY       =>  plt_morph%MY     , &
    NI       =>  plt_morph%NI     , &
    NumRootAxes_pft     =>  plt_morph%NumRootAxes_pft   , &
    CanopyBranchLeafA_pft    =>  plt_morph%CanopyBranchLeafA_pft  , &
    CanopyStemA_pft    =>  plt_morph%CanopyStemA_pft  , &
    CanopyBranchStemApft_lyr    =>  plt_morph%CanopyBranchStemApft_lyr  , &
    GRNOB    =>  plt_morph%GRNOB  , &
    CanopyLeafA_pft    =>  plt_morph%CanopyLeafA_pft  , &
    CanopyStemApft_lyr    =>  plt_morph%CanopyStemApft_lyr  , &
    iPlantNfixType   =>  plt_morph%iPlantNfixType , &
    GRNO     =>  plt_morph%GRNO     &
  )
!
!     ACCUMULATE PFT STATE VARIABLES FROM BRANCH STATE VARIABLES
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
!     ReserveChemElmnts_brch,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     CanopyBranchLeafA_pft=branch leaf area
!     CanopyBranchStemApft_lyr=total branch stalk surface area in each layer
!     GRNOB=seed set number
!
  DO NE=1,NumOfPlantChemElmnts
    CanopyNonstructElements_pft(NE,NZ)=sum(NonstructElmnt_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    ShootChemElmnts_pft(NE,NZ)=sum(WTSHTBE(NE,1:NumOfBranches_pft(NZ),NZ))
    SheathChemElmnts(NE,NZ)=sum(PetioleChemElmnts_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    StalkChemElmnts(NE,NZ)=sum(StalkChemElmnts_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    LeafChemElmnts(NE,NZ)=sum(LeafChemElmnts_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    CanopyReserveChemElmnts(NE,NZ)=sum(ReserveChemElmnts_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    HuskChemElmnts(NE,NZ)=sum(HuskChemElmnts_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    GrainChemElmnts(NE,NZ)=sum(GrainChemElmnts_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    EarChemElmnts(NE,NZ)=sum(EarChemElmnts_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    !root state variables
    RootChemElmnts_pft(NE,NZ)=sum( RootMycoNonstructElmnt_vr(NE,1:MY(NZ),NU:NJ,NZ))
    WTRTSE(NE,NZ)=sum(WTRT1E(NE,1:MY(NZ),NU:NJ,1:NumRootAxes_pft(NZ),NZ)) &
      +sum(WTRT2E(NE,1:MY(NZ),NU:NJ,1:NumRootAxes_pft(NZ),NZ))
    RootChemElmnts_pft(NE,NZ)=RootChemElmnts_pft(NE,NZ)+WTRTSE(NE,NZ)
  ENDDO

  CanPStalkC(NZ)=sum(StalkBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopyLeafShethC_pft(NZ) =sum(LeafPetioleBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))
  GRNO(NZ) =sum(GRNOB(1:NumOfBranches_pft(NZ),NZ))
  CanopyLeafA_pft(NZ)=sum(CanopyBranchLeafA_pft(1:NumOfBranches_pft(NZ),NZ))
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
  IF(iPlantNfixType(NZ).NE.0)THEN
    IF(iPlantNfixType(NZ).GE.4)THEN
      DO NE=1,NumOfPlantChemElmnts
        D7950: DO NB=1,NumOfBranches_pft(NZ)
          NoduleNonstructElmnt_pft(NE,NZ)=NoduleNonstructElmnt_pft(NE,NZ)+NoduleNonstructElmnt_brch(NE,NB,NZ)
        ENDDO D7950
        NoduleChemElmnts_pft(NE,NZ)=sum(WTNDBE(NE,1:NumOfBranches_pft(NZ),NZ))+sum(NoduleNonstructElmnt_brch(NE,1:NumOfBranches_pft(NZ),NZ))
      ENDDO
    ELSEIF(iPlantNfixType(NZ).GE.1.AND.iPlantNfixType(NZ).LE.3)THEN
      DO NE=1,NumOfPlantChemElmnts
        NoduleChemElmnts_pft(NE,NZ)=sum(WTNDLE(NE,NU:NI(NZ),NZ))+sum(RootNoduleNonstructElmnt_vr(NE,NU:NI(NZ),NZ))
      ENDDO
    ENDIF
  ENDIF
!
!     ACCUMULATE TOTAL SOIL-PLANT C,N,P EXCHANGE
!
!     HCUPTK,HZUPTK,HPUPTK=net PFT root-soil C,N,P exchange
!     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
!     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
!     UPNF=PFT N2 fixation
!     TCUPTK,TZUPTK,TPUPTK=cumulative PFT root-soil C,N,P exchange
!     TZUPFX=cumulative PFT N2 fixation
!
  HEUPTK(1:NumOfPlantChemElmnts,NZ)=UPOME(1:NumOfPlantChemElmnts,NZ)
  HEUPTK(ielmn,NZ)=HEUPTK(ielmn,NZ)+UPNH4(NZ)+UPNO3(NZ)+UPNF(NZ)
  HEUPTK(ielmp,NZ)=HEUPTK(ielmp,NZ)+UPH2P(NZ)+UPH1P(NZ)

  PlantExudChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=PlantExudChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)+UPOME(1:NumOfPlantChemElmnts,NZ)
  PlantExudChemElmnts_pft(ielmn,NZ)=PlantExudChemElmnts_pft(ielmn,NZ)+UPNH4(NZ)+UPNO3(NZ)
  PlantExudChemElmnts_pft(ielmp,NZ)=PlantExudChemElmnts_pft(ielmp,NZ)+UPH2P(NZ)+UPH1P(NZ)
  TZUPFX(NZ)=TZUPFX(NZ)+UPNF(NZ)+UPNFC(NZ)
  end associate
  end subroutine AccumulateStates

end module grosubsMod
