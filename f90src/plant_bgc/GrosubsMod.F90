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
  use LitrFallMod
  use PlantBranchMod
  use PlantBalMod
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
!     HourReq2InitSStor4LeafOut=number of hours required to initiate remobilization of storage C for leafout
!     GVMX=specific oxidation rate of nonstructural C during leafout at 25 C
!     FracHour4LeafoffRemob=fraction of hours required for leafoff to initiate remobilization
!

  integer  :: curday,curhour
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
  real(r8) :: ShootC4NonstC_brch(NumOfCanopyLayers1,JP1)
! begin_execution
  associate(                                                      &
    NP0                     =>  plt_site%NP0                    , &  
    IsPlantActive_pft       =>  plt_pheno%IsPlantActive_pft     , &
    CanopyHeight_pft        =>  plt_morph%CanopyHeight_pft      , &        
    NP                      =>  plt_site%NP                       &
  )
!     TOTAL AGB FOR GRAZING IN LANDSCAPE SECTION
!
!
!     INITIALIZE SENESCENCE ARRAYS
  DO NZ=1,NP0
    CanopyHeight_copy(NZ)=CanopyHeight_pft(NZ)
    CanopyHeight_pft(NZ)=0._r8
  ENDDO  
!
!     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
!
  D9985: DO NZ=1,NP

! IsPlantActive_pft= flag for living pft
    IF(IsPlantActive_pft(NZ).EQ.iActive)THEN
      call GrowPlant(I,J,NZ,CanopyHeight_copy)
    ENDIF

    call SumPlantBiom(I,J,NZ,'bfdistb')
!   HARVEST STANDING DEAD
    call RemoveBiomassByDisturbance(I,J,NZ)
    
    call SumPlantBiom(I,J,NZ,'bflvdeadtrns')
  ENDDO D9985
!
! TRANSFORMATIONS IN LIVING OR DEAD PLANT POPULATIONS
  
  call LiveDeadTransformation(I,J)
  DO NZ=1,NP
    call SumPlantBiom(I,J,NZ,'exgrosubs')
  ENDDO
  DO NZ=1,2
    write(111,*)''
    write(112,*)''
  ENDDO
  end associate
  END subroutine grosubs

!------------------------------------------------------------------------------------------

  subroutine LiveDeadTransformation(I,J)
  implicit none
  integer, intent(in) :: I,J

  integer :: L,K,NZ,M,NE,NB
  real(r8) :: XFRC,XFRN,XFRP,XFRE
!     begin_execution

  associate(                                                               &
    k_fine_litr                    => pltpar%k_fine_litr                 , &
    k_woody_litr                   => pltpar%k_woody_litr                , &
    iDayPlanting_pft               => plt_distb%iDayPlanting_pft         , &
    iYearPlanting_pft              => plt_distb%iYearPlanting_pft        , &
    EcoHavstElmntCum_pft           => plt_distb%EcoHavstElmntCum_pft     , &
    EcoHavstElmnt_pft              => plt_distb%EcoHavstElmnt_pft        , &
    PO4byFire_pft                  => plt_distb%PO4byFire_pft            , &
    N2ObyFire_pft                  => plt_distb%N2ObyFire_pft            , &
    CH4ByFire_pft                  => plt_distb%CH4ByFire_pft            , &
    CO2ByFire_pft                  => plt_distb%CO2ByFire_pft            , &
    NH3byFire_pft                  => plt_distb%NH3byFire_pft            , &
    StandDeadKCompElms_pft         => plt_biom%StandDeadKCompElms_pft    , &
    RootElms_pft                   => plt_biom%RootElms_pft              , &
    ShootStrutElms_pft             => plt_biom%ShootStrutElms_pft        , &
    NodulStrutElms_pft             => plt_biom%NodulStrutElms_pft        , &
    SeasonalNonstElms_pft          => plt_biom%SeasonalNonstElms_pft     , &
    StandDeadStrutElms_pft         => plt_biom%StandDeadStrutElms_pft    , &
    fTgrowCanP                     => plt_pheno%fTgrowCanP               , &
    NetCumElmntFlx2Plant_pft       => plt_pheno%NetCumElmntFlx2Plant_pft , &
    IsPlantActive_pft              => plt_pheno%IsPlantActive_pft        , &
    iPlantRootProfile_pft          => plt_pheno%iPlantRootProfile_pft    , &
    doInitPlant_pft                => plt_pheno%doInitPlant_pft          , &
    doPlantLeafOut_brch            => plt_pheno%doPlantLeafOut_brch      , &
    iPlantTurnoverPattern_pft      => plt_pheno%iPlantTurnoverPattern_pft, &
    HourReq4LeafOut_brch           => plt_pheno%HourReq4LeafOut_brch     , &
    Hours4Leafout_brch             => plt_pheno%Hours4Leafout_brch       , &
    ElmBalanceCum_pft              => plt_site%ElmBalanceCum_pft         , &
    NP0                            => plt_site%NP0                       , &
    MaxNumRootLays                 => plt_site%MaxNumRootLays            , &
    iYearCurrent                   => plt_site%iYearCurrent              , &
    LitrfalStrutElms_pvr           => plt_bgcr%LitrfalStrutElms_pvr      , &
    NetPrimProduct_pft             => plt_bgcr%NetPrimProduct_pft        , &
    PlantN2FixCum_pft              => plt_bgcr%PlantN2FixCum_pft         , &
    LitrfalStrutElms_pft           => plt_bgcr%LitrfalStrutElms_pft      , &
    NH3EmiCum_pft                  => plt_bgcr%NH3EmiCum_pft             , &
    LitrfalStrutElmsCum_pft        => plt_bgcr%LitrfalStrutElmsCum_pft     , &
    SurfLitrfalStrutElmsCum_pft    => plt_bgcr%SurfLitrfalStrutElmsCum_pft , &
    GrossResp_pft                  => plt_bgcr%GrossResp_pft               , &
    GrossCO2Fix_pft                => plt_bgcr%GrossCO2Fix_pft             , &
    PlantExudChemElmCum_pft        => plt_rbgc%PlantExudChemElmCum_pft     , &
    CumSoilThickness               => plt_site%CumSoilThickness            , &
    PlantinDepth                   => plt_morph%PlantinDepth               , &
    NumOfBranches_pft              => plt_morph%NumOfBranches_pft            &
  )
  D9975: DO NZ=1,NP0
!
!     ACTIVATE DORMANT SEEDS
!
    D205: DO NB=1,NumOfBranches_pft(NZ)
      IF(doInitPlant_pft(NZ).EQ.itrue)THEN
        IF(doPlantLeafOut_brch(NB,NZ).EQ.iEnable.AND.Hours4Leafout_brch(NB,NZ).GE.HourReq4LeafOut_brch(NB,NZ))THEN
          iDayPlanting_pft(NZ)=I
          iYearPlanting_pft(NZ)=iYearCurrent
          PlantinDepth(NZ)=0.005_r8+CumSoilThickness(0)
          !mark plant as initialized
          doInitPlant_pft(NZ)=ifalse
        ENDIF
      ENDIF
    ENDDO D205
!
!     LitrFall FROM STANDING DEAD
!
!     XFRC,XFRN,XFRP=LitrFall from standing dead
!     fTgrowCanP=temperature function for canopy growth
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     CSNC,ZSNC,PSNC=C,N,P LitrFall
!
    
    D6235: DO M=1,jsken
      DO NE=1,NumPlantChemElms
        XFRE=1.5814E-05_r8*fTgrowCanP(NZ)*StandDeadKCompElms_pft(NE,M,NZ)
        IF(iPlantTurnoverPattern_pft(NZ).EQ.0.OR.iPlantRootProfile_pft(NZ).LE.1)THEN
          LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_fine_litr,0,NZ)+XFRE
        ELSE
          LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)=LitrfalStrutElms_pvr(NE,M,k_woody_litr,0,NZ)+XFRE
        ENDIF
        StandDeadKCompElms_pft(NE,M,NZ)=StandDeadKCompElms_pft(NE,M,NZ)-XFRE
      ENDDO
    ENDDO D6235
!
!     ACCUMULATE TOTAL SURFACE, SUBSURFACE LitrFall
!
!     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P LitrFall
!     TCSNC,TZSNC,TPSNC=cumulative C,N,P LitrFall
!     HCSNC,HZSNC,HPSNC=hourly C,N,P LitrFall
!   diagnose surface literfall 
    DO K=1,pltpar%NumOfPlantLitrCmplxs      
      D6431: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          SurfLitrfalStrutElmsCum_pft(NE,NZ)=SurfLitrfalStrutElmsCum_pft(NE,NZ)+LitrfalStrutElms_pvr(NE,M,K,0,NZ)
        ENDDO
      ENDDO D6431  
    ENDDO
!   the following purposely starts from layer 0.
    LitrfalStrutElms_pft(1:NumPlantChemElms,NZ)=0._r8
    D8955: DO L=0,MaxNumRootLays
      DO K=1,pltpar%NumOfPlantLitrCmplxs      
        D6430: DO M=1,jsken
          DO NE=1,NumPlantChemElms
            LitrfalStrutElms_pft(NE,NZ)=LitrfalStrutElms_pft(NE,NZ)+LitrfalStrutElms_pvr(NE,M,K,L,NZ)
          enddo         
        ENDDO D6430      
      ENDDO
    ENDDO D8955    

    DO NE=1,NumPlantChemElms
      LitrfalStrutElmsCum_pft(NE,NZ)=LitrfalStrutElmsCum_pft(NE,NZ)+LitrfalStrutElms_pft(NE,NZ)
    ENDDO

!
!     PLANT C BALANCE = TOTAL C STATE VARIABLES + TOTAL
!     AUTOTROPHIC RESPIRATION + TOTAL LitrFall - TOTAL EXUDATION
!     - TOTAL CO2 FIXATION
!
!     BALC=PFT C balance
!     WTSHT,WTRT,WTND,WTRVC,WTSTG=PFT shoot,root,bacteria,storage,standing dead C
!     NetPrimProduct_pft=cumulative PFT NPP
!     TCSNC=cumulative PFT C LitrFall
!     TCUPTK=cumulative PFT root-soil C exchange
!     RSETC=cumulative C balance from previous year
!     THVSTC=cumulative PFT C removed from ecosystem from previous year
!     HVSTC=total PFT C removed from ecosystem in current year
!     CO2ByFire_pft,CH4ByFire_pft=CO2,CH4 emission from disturbance
!
    NetPrimProduct_pft(NZ)=GrossCO2Fix_pft(NZ)+GrossResp_pft(NZ)

    IF(IsPlantActive_pft(NZ).EQ.iActive)THEN
      !check for living plant
      DO NE=1,NumPlantChemElms
        ElmBalanceCum_pft(NE,NZ)=ShootStrutElms_pft(NE,NZ)+RootElms_pft(NE,NZ)+NodulStrutElms_pft(NE,NZ) &
          +SeasonalNonstElms_pft(NE,NZ)+LitrfalStrutElmsCum_pft(NE,NZ)-PlantExudChemElmCum_pft(NE,NZ) &
          -NetCumElmntFlx2Plant_pft(NE,NZ)+StandDeadStrutElms_pft(NE,NZ)&
          +EcoHavstElmnt_pft(NE,NZ)+EcoHavstElmntCum_pft(NE,NZ)
      ENDDO
      ElmBalanceCum_pft(ielmc,NZ)=ElmBalanceCum_pft(ielmc,NZ)-NetPrimProduct_pft(NZ) &
        -CO2ByFire_pft(NZ)-CH4ByFire_pft(NZ)
!
!     PLANT N BALANCE = TOTAL N STATE VARIABLES + TOTAL N LitrFall
!     - TOTAL N UPTAKE FROM SOIL - TOTAL N ABSORPTION FROM ATMOSPHERE
!
!     BALN=PFT N balance
!     WTSHN,WTRTN,WTNDN,WTRVN,WTSTGN=PFT shoot,root,bacteria,storage,standing dead N
!     TZSNC=cumulative PFT N LitrFall
!     TZUPTK=cumulative PFT root-soil N exchange
!     NH3EmiCum_pft=cumulative NH3 exchange
!     RSETN=cumulative N balance from previous year
!     THVSTN=cumulative PFT N removed from ecosystem from previous year
!     HVSTN=total PFT N removed from ecosystem in current year
!     NH3byFire_pft,N2ObyFire_pft=NH3,N2O emission from disturbance
!     PlantN2FixCum_pft=cumulative PFT N2 fixation
!
      ElmBalanceCum_pft(ielmn,NZ)=ElmBalanceCum_pft(ielmn,NZ)-NH3EmiCum_pft(NZ) &
        -NH3byFire_pft(NZ)-N2ObyFire_pft(NZ)-PlantN2FixCum_pft(NZ)
!
!     PLANT P BALANCE = TOTAL P STATE VARIABLES + TOTAL P LitrFall
!     - TOTAL P UPTAKE FROM SOIL
!
!     BALP=PFT N balance
!     WTSHP,WTRTP,WTNDP,WTRVP,WTSTGP=PFT shoot,root,bacteria,storage,standing dead P
!     TPSNC=cumulative PFT P LitrFall
!     TPUPTK=cumulative PFT root-soil P exchange
!     RSETP=cumulative P balance from previous year
!     THVSTP=cumulative PFT P removed from ecosystem from previous year
!     HVSTP=total PFT P removed from ecosystem in current year
!     PO4byFire_pft=PO4 emission from disturbance
!
      ElmBalanceCum_pft(ielmp,NZ)=ElmBalanceCum_pft(ielmp,NZ)-PO4byFire_pft(NZ)
    ENDIF
  ENDDO D9975
  end associate
  end subroutine LiveDeadTransformation
!------------------------------------------------------------------------------------------

  subroutine GrowPlant(I,J,NZ,CanopyHeight_copy)
  !
  !Description
  !plant growth
  use PlantDisturbsMod, only : RemoveBiomByManagement
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: CanopyHeight_copy(JP1)

  real(r8)  :: CanopyN2Fix_pft(JP1)
  integer  :: NB
  integer  :: BegRemoblize
  real(r8) :: TFN6_vr(JZ1)
  real(r8) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8) :: PTRT
  real(r8) :: RootAreaPopu
  real(r8) :: TFN5
  real(r8) :: WFNG
  real(r8) :: Stomata_Activity
  real(r8) :: WFNS,WFNSG
! begin_execution
  associate(                                                            &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft      , &
    iPlantRootState_pft       => plt_pheno%iPlantRootState_pft        , &
    iPlantShootState_pft      => plt_pheno%iPlantShootState_pft       , &
    RootN2Fix_pft             => plt_rbgc%RootN2Fix_pft               , &
    RootH2PO4Uptake_pft       => plt_rbgc%RootH2PO4Uptake_pft         , &
    RootNH4Uptake_pft         => plt_rbgc%RootNH4Uptake_pft           , &
    RootHPO4Uptake_pft        => plt_rbgc%RootHPO4Uptake_pft          , &
    RootNO3Uptake_pft         => plt_rbgc%RootNO3Uptake_pft           , &
    PlantRootSoilElmNetX_pft  => plt_rbgc%PlantRootSoilElmNetX_pft    , &
    RootMycoExudElms_pft      => plt_rbgc%RootMycoExudElms_pft        , &
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft          , &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft              &
  )
  IF(iPlantShootState_pft(NZ).EQ.iLive .OR. iPlantRootState_pft(NZ).EQ.iLive)THEN
    CanopyN2Fix_pft(NZ)=0._r8
    BegRemoblize = 0

!    call SumPlantBiom(I,J,NZ,'bfstageplant')
    
    call StagePlantForGrowth(I,J,NZ,TFN6_vr,CNLFW,CPLFW,&
      CNSHW,CPSHW,CNRTW,CPRTW,RootAreaPopu,TFN5,WFNG,Stomata_Activity,WFNS,WFNSG)
!
!     CALCULATE GROWTH OF EACH BRANCH
!
!     WTLFB,WTSHEB,LeafPetolBiomassC_brch=leaf,petiole,leaf+petiole mass
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!
    call SumPlantBiom(I,J,NZ,'bfgrowpbrch')
    
    DO  NB=1,NumOfBranches_pft(NZ)
      call GrowOneBranch(I,J,NB,NZ,TFN6_vr,CanopyHeight_copy,CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW,&
        TFN5,WFNG,Stomata_Activity,WFNS,WFNSG,PTRT,CanopyN2Fix_pft,BegRemoblize)
    ENDDO
!
    call SumPlantBiom(I,J,NZ,'bfRootBGCM')
    call RootBGCModel(I,J,NZ,BegRemoblize,PTRT,TFN6_vr,CNRTW,CPRTW,RootAreaPopu)
!
    call ComputeTotalBiom(I,J,NZ)
  ENDIF
!
  call SumPlantBiom(I,J,NZ,'bfrmbiom')
  call RemoveBiomByManagement(I,J,NZ)
!
!     RESET DEAD BRANCHES
  call SumPlantBiom(I,J,NZ,'bfresetdead')
  call ResetDeadBranch(I,J,NZ)
!  
  call AccumulateStates(I,J,NZ,CanopyN2Fix_pft)
  end associate
  end subroutine GrowPlant

!------------------------------------------------------------------------------------------
  subroutine StagePlantForGrowth(I,J,NZ,TFN6_vr,CNLFW,CPLFW,CNSHW,&
    CPSHW,CNRTW,CPRTW,RootAreaPopu,TFN5,WFNG,Stomata_Activity,WFNS,WFNSG)
  integer, intent(in) :: I,J,NZ
  REAL(R8), INTENT(OUT):: TFN6_vr(JZ1)
  REAL(R8), INTENT(OUT) :: CNLFW,CPLFW,CNSHW,CPSHW,CNRTW,CPRTW
  real(r8), intent(out) :: RootAreaPopu,TFN5,WFNG,Stomata_Activity
  real(r8), intent(out) :: WFNS,WFNSG
  integer :: L,NR,N,NE
  real(r8) :: ACTVM,RTK,STK,TKCM,TKSM
!     begin_execution

  associate(                            &
    TKC                       =>  plt_ew%TKC         , &
    TKS                       =>  plt_ew%TKS         , &
    PSICanopy_pft             =>  plt_ew%PSICanopy_pft      , &
    PSICanopyTurg_pft         =>  plt_ew%PSICanopyTurg_pft       , &
    PlantPopulation_pft       =>  plt_site%PlantPopulation_pft        , &
    NU                        =>  plt_site%NU        , &
    MaxNumRootLays            =>  plt_site%MaxNumRootLays        , &
    OFFST                     =>  plt_pheno%OFFST    , &
    ZEROP                     =>  plt_biom%ZEROP     , &
    CanopyStalkC_pft          =>  plt_biom%CanopyStalkC_pft     , &
    RootProteinC_pvr          =>  plt_biom%RootProteinC_pvr     , &
    RootBiomCPerPlant_pft     =>  plt_biom%RootBiomCPerPlant_pft     , &
    StalkStrutElms_pft        =>  plt_biom%StalkStrutElms_pft    , &
    RootElms_pft              =>  plt_biom%RootElms_pft     , &
    CanopyLeafCLyr_pft        =>  plt_biom%CanopyLeafCLyr_pft     , &
    iPlantTurnoverPattern_pft =>  plt_pheno%iPlantTurnoverPattern_pft   , &
    iPlantRootProfile_pft     =>  plt_pheno%iPlantRootProfile_pft   , &
    RCO2A_pvr                 =>  plt_rbgc%RCO2A_pvr     , &
    RootRespPotent_pvr        =>  plt_rbgc%RootRespPotent_pvr     , &
    RCO2N_pvr                 =>  plt_rbgc%RCO2N_pvr     , &
    RootrNC_pft               =>  plt_allom%RootrNC_pft    , &
    RootrPC_pft               =>  plt_allom%RootrPC_pft    , &
    FracHour4LeafoffRemob     =>  plt_allom%FracHour4LeafoffRemob     , &
    FWODLE                    =>  plt_allom%FWODLE   , &
    FWODBE                    =>  plt_allom%FWODBE   , &
    FWOODE                    =>  plt_allom%FWOODE   , &
    FWODRE                    =>  plt_allom%FWODRE   , &
    CNLF                      =>  plt_allom%CNLF     , &
    CPLF                      =>  plt_allom%CPLF     , &
    CNSHE                     =>  plt_allom%CNSHE    , &
    CPSHE                     =>  plt_allom%CPSHE    , &
    rNCStalk_pft              =>  plt_allom%rNCStalk_pft   , &
    rPCStalk_pft              =>  plt_allom%rPCStalk_pft    , &
    k_fine_litr               =>  pltpar%k_fine_litr,&
    k_woody_litr              =>  pltpar%k_woody_litr,&
    RCS                       =>  plt_photo%RCS      , &
    Root1stXNumL_pvr          =>  plt_morph%Root1stXNumL_pvr    , &
    Root2ndXNum_pvr           =>  plt_morph%Root2ndXNum_pvr     , &
    MY                        =>  plt_morph%MY       , &
    CanopyLeafALyr_pft        =>  plt_morph%CanopyLeafALyr_pft    , &
    CanopyStemArea_lpft       =>  plt_morph%CanopyStemArea_lpft    , &
    NumRootAxes_pft           =>  plt_morph%NumRootAxes_pft       &
  )
  D2: DO L=1,NumOfCanopyLayers1
    CanopyLeafALyr_pft(L,NZ)=0._r8
    CanopyLeafCLyr_pft(L,NZ)=0._r8
    CanopyStemArea_lpft(L,NZ)=0._r8
  ENDDO D2

  D9: DO N=1,MY(NZ)
    D6: DO L=NU,MaxNumRootLays
      RootProteinC_pvr(N,L,NZ)=0._r8
      Root1stXNumL_pvr(N,L,NZ)=0._r8
      Root2ndXNum_pvr(N,L,NZ)=0._r8
      RootRespPotent_pvr(N,L,NZ)=0._r8
      RCO2N_pvr(N,L,NZ)=0._r8
      RCO2A_pvr(N,L,NZ)=0._r8
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
    .OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ)))&
    .OR.StalkStrutElms_pft(ielmc,NZ).LE.ZEROP(NZ))THEN
    FWODBE(ielmc,k_fine_litr)=1.0_r8
    FWOODE(ielmc,k_fine_litr)=1.0_r8
    FWODRE(ielmc,k_fine_litr)=1.0_r8
  ELSE
    FWODBE(ielmc,k_fine_litr)=1.0_r8
    FWOODE(ielmc,k_fine_litr)=SQRT(CanopyStalkC_pft(NZ)/StalkStrutElms_pft(ielmc,NZ))
    FWODRE(ielmc,k_fine_litr)=SQRT(FRTX*CanopyStalkC_pft(NZ)/StalkStrutElms_pft(ielmc,NZ))
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

  DO NE=2,NumPlantChemElms
    FWODLE(NE,k_fine_litr)=1.0_r8-FWODLE(NE,k_woody_litr)
    FWODBE(NE,k_fine_litr)=1.0_r8-FWODBE(NE,k_woody_litr)
    FWOODE(NE,k_fine_litr)=1.0_r8-FWOODE(NE,k_woody_litr)
    FWODRE(NE,k_fine_litr)=1.0_r8-FWODRE(NE,k_woody_litr)
  ENDDO
!
!     SHOOT AND ROOT TEMPERATURE FUNCTIONS FOR MAINTENANCE
!     RESPIRATION FROM TEMPERATURES WITH OFFSETS FOR THERMAL ADAPTATION
!
!     TKC,TKCM=canopy temperature,canopy temp used in Arrhenius eqn
!     TKS,TKSM=soil temperature,soil temp used in Arrhenius eqn
!     OFFST=shift in Arrhenius curve for thermal adaptation
!     TFN5,TFN6_vr=temperature function for canopy,root mntc respn (25 oC =1)
!     8.3143,710.0=gas constant,enthalpy
!     62500,195000,232500=energy of activn,high,low temp inactivn(KJ mol-1)
!
  TKCM=TKC(NZ)+OFFST(NZ)

  TFN5=calc_plant_maint_tempf(TKCM)  
  D7: DO L=NU,MaxNumRootLays
    TKSM=TKS(L)+OFFST(NZ)
    TFN6_vr(L)=calc_plant_maint_tempf(TKSM)  
  ENDDO D7
!
!     PRIMARY ROOT NUMBER
!
!     RootBiomCPerPlant_pft=root mass per plant used to calculate primary root number
!     WTRT,PP=root mass,PFT population
!     RootAreaPopu=multiplier for number of primary root axes
!
  RootBiomCPerPlant_pft(NZ)=AMAX1(0.999992087_r8*RootBiomCPerPlant_pft(NZ),&
    RootElms_pft(ielmc,NZ)/PlantPopulation_pft(NZ))
  RootAreaPopu=AMAX1(1.0_r8,RootBiomCPerPlant_pft(NZ)**0.667_r8)*PlantPopulation_pft(NZ)
!
!     WATER STRESS FUNCTIONS FOR EXPANSION AND GROWTH RESPIRATION
!     FROM CANOPY TURGOR
!
!     WFNS=turgor expansion,extension function
!     PSICanopyTurg_pft,PSIMin4OrganExtens=current,minimum canopy turgor potl for expansion,extension
!     Stomata_Activity=stomatal resistance function of canopy turgor
!     PSICanopy_pft=canopy water potential
!     WFNG=growth function of canopy water potential
!     WFNSG=expansion,extension function of canopy water potential
!
  WFNS=AMIN1(1.0_r8,AZMAX1(PSICanopyTurg_pft(NZ)-PSIMin4OrganExtens))

  IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
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

  subroutine ComputeTotalBiom(I,J,NZ)

  integer, intent(in) :: I,J,NZ
  
  integer :: L,N
!     begin_execution
  associate(                                                             &
    NU                             =>  plt_site%NU                     , &  
    MY                             =>  plt_morph%MY                    , &        
    MaxSoiL4Root                   =>  plt_morph%MaxSoiL4Root          , &    
    NumRootAxes_pft                =>  plt_morph%NumRootAxes_pft       , &
    MaxNumRootLays                 =>  plt_site%MaxNumRootLays         , &        
    RootMycoNonstElms_rpvr         =>  plt_biom%RootMycoNonstElms_rpvr , &
    PopuRootMycoC_pvr              =>  plt_biom% PopuRootMycoC_pvr       &    
  )

  call SumPlantBiom(I,J,NZ,'computotb')
! add the nonstrucal components
  D3451: DO N=1,MY(NZ)
    DO  L=NU,MaxSoiL4Root(NZ)
      PopuRootMycoC_pvr(N,L,NZ)=PopuRootMycoC_pvr(N,L,NZ)+RootMycoNonstElms_rpvr(ielmc,N,L,NZ)
    enddo
  ENDDO D3451

  end associate

  end subroutine ComputeTotalBiom
!------------------------------------------------------------------------------------------

  subroutine AccumulateStates(I,J,NZ,CanopyN2Fix_pft)
  
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: CanopyN2Fix_pft(JP1)
  integer :: L,NR,N,NE,NB
  real(r8) :: root1st, root2nd
!     begin_execution
  associate(                                                                &
    CanopyNonstElms_brch          =>  plt_biom%CanopyNonstElms_brch      , &
    RootMycoNonstElms_rpvr        =>  plt_biom%RootMycoNonstElms_rpvr    , &
    CanopyNodulNonstElms_brch     =>  plt_biom%CanopyNodulNonstElms_brch , &
    CanopyNodulElms_pft           =>  plt_biom%CanopyNodulElms_pft       , &
    RootNodulElms_pft             =>  plt_biom%RootNodulElms_pft         , &
    RootStrutElms_pft             =>  plt_biom%RootStrutElms_pft         , &
    RootElms_pft                  =>  plt_biom%RootElms_pft              , &
    CanopyNodulStrutElms_brch     =>  plt_biom%CanopyNodulStrutElms_brch , &
    ShootStrutElms_brch           =>  plt_biom%ShootStrutElms_brch       , &
    StalkStrutElms_brch           =>  plt_biom%StalkStrutElms_brch       , &
    HuskStrutElms_brch            =>  plt_biom%HuskStrutElms_brch        , &
    StalkRsrvElms_brch            =>  plt_biom%StalkRsrvElms_brch        , &
    EarStrutElms_brch             =>  plt_biom%EarStrutElms_brch         , &
    LeafPetolBiomassC_brch        =>  plt_biom%LeafPetolBiomassC_brch    , &
    GrainStrutElms_brch           =>  plt_biom%GrainStrutElms_brch       , &
    LeafStrutElms_brch            =>  plt_biom%LeafStrutElms_brch        , &
    StalkBiomassC_brch            =>  plt_biom%StalkBiomassC_brch        , &
    PetoleStrutElms_brch          =>  plt_biom%PetoleStrutElms_brch      , &
    ShootStrutElms_pft            =>  plt_biom%ShootStrutElms_pft        , &
    LeafStrutElms_pft             =>  plt_biom%LeafStrutElms_pft         , &
    PetioleStrutElms_pft          =>  plt_biom%PetioleStrutElms_pft      , &
    StalkStrutElms_pft            =>  plt_biom%StalkStrutElms_pft        , &
    CanopyStalkC_pft              =>  plt_biom%CanopyStalkC_pft          , &
    StalkRsrvElms_pft             =>  plt_biom%StalkRsrvElms_pft         , &
    HuskStrutElms_pft             =>  plt_biom%HuskStrutElms_pft         , &
    EarStrutElms_pft              =>  plt_biom%EarStrutElms_pft          , &
    GrainStrutElms_pft            =>  plt_biom%GrainStrutElms_pft        , &
    CanopyLeafShethC_pft          =>  plt_biom%CanopyLeafShethC_pft      , &
    RootNodulStrutElms_pvr        =>  plt_biom%RootNodulStrutElms_pvr    , &
    CanopyNonstElms_pft           =>  plt_biom%CanopyNonstElms_pft  , &
    CanopyNodulNonstElms_pft      =>  plt_biom%CanopyNodulNonstElms_pft  , &
    RootMyco1stStrutElms_rpvr     =>  plt_biom%RootMyco1stStrutElms_rpvr  , &
    RootMyco2ndStrutElms_rpvr     =>  plt_biom%RootMyco2ndStrutElms_rpvr  , &
    RootNodulNonstElms_pvr        =>  plt_biom%RootNodulNonstElms_pvr , &
    PlantN2FixCum_pft             =>  plt_bgcr%PlantN2FixCum_pft  , &
    PlantExudChemElmCum_pft       =>  plt_rbgc%PlantExudChemElmCum_pft  , &
    RootHPO4Uptake_pft            =>  plt_rbgc%RootHPO4Uptake_pft   , &
    PlantRootSoilElmNetX_pft      =>  plt_rbgc%PlantRootSoilElmNetX_pft  , &
    RootN2Fix_pft                 =>  plt_rbgc%RootN2Fix_pft    , &
    RootH2PO4Uptake_pft           =>  plt_rbgc%RootH2PO4Uptake_pft   , &
    RootNO3Uptake_pft             =>  plt_rbgc%RootNO3Uptake_pft   , &
    RootNH4Uptake_pft             =>  plt_rbgc%RootNH4Uptake_pft   , &
    RootMycoExudElms_pft          =>  plt_rbgc%RootMycoExudElms_pft   , &
    MaxNumRootLays                =>  plt_site%MaxNumRootLays      , &
    NU                            =>  plt_site%NU      , &
    NumOfBranches_pft             =>  plt_morph%NumOfBranches_pft    , &
    MY                            =>  plt_morph%MY     , &
    MaxSoiL4Root                  =>  plt_morph%MaxSoiL4Root     , &
    NumRootAxes_pft               =>  plt_morph%NumRootAxes_pft   , &
    LeafAreaLive_brch             =>  plt_morph%LeafAreaLive_brch  , &
    CanopyStemArea_pft            =>  plt_morph%CanopyStemArea_pft  , &
    CanopyStemArea_lbrch          =>  plt_morph%CanopyStemArea_lbrch  , &
    SeedNumSet_brch               =>  plt_morph%SeedNumSet_brch  , &
    CanopyLeafArea_pft            =>  plt_morph%CanopyLeafArea_pft  , &
    CanopyStemArea_lpft           =>  plt_morph%CanopyStemArea_lpft  , &
    iPlantNfixType                =>  plt_morph%iPlantNfixType , &
    CanopySeedNum_pft             =>  plt_morph%CanopySeedNum_pft     &
  )
!
!     ACCUMULATE PFT STATE VARIABLES FROM BRANCH STATE VARIABLES
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
!     StalkRsrvElms_brch,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
!     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
!     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     LeafAreaLive_brch=branch leaf area
!     CanopyStemArea_lbrch=total branch stalk surface area in each layer
!     SeedNumSet_brch=seed set number
!
  call SumPlantBiom(I,J,NZ,'Accumstates')
  
  DO NE=1,NumPlantChemElms
    CanopyNonstElms_pft(NE,NZ)=sum(CanopyNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    ShootStrutElms_pft(NE,NZ)=sum(ShootStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    PetioleStrutElms_pft(NE,NZ)=sum(PetoleStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    StalkStrutElms_pft(NE,NZ)=sum(StalkStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    LeafStrutElms_pft(NE,NZ)=sum(LeafStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    StalkRsrvElms_pft(NE,NZ)=sum(StalkRsrvElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    if(StalkRsrvElms_pft(NE,NZ)>1.e20)then
      print*,NZ,NE,I+J/24.,StalkRsrvElms_pft(NE,NZ),NumOfBranches_pft(NZ)
      print*,'755',StalkRsrvElms_brch(NE,1,NZ)
      stop
    endif
    HuskStrutElms_pft(NE,NZ)=sum(HuskStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    GrainStrutElms_pft(NE,NZ)=sum(GrainStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    EarStrutElms_pft(NE,NZ)=sum(EarStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    !root state variables
    !sum structural biomass

  ENDDO

  CanopyStalkC_pft(NZ)=sum(StalkBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopyLeafShethC_pft(NZ) =sum(LeafPetolBiomassC_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopySeedNum_pft(NZ) =sum(SeedNumSet_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopyLeafArea_pft(NZ)=sum(LeafAreaLive_brch(1:NumOfBranches_pft(NZ),NZ))
  CanopyStemArea_pft(NZ)=sum(CanopyStemArea_lbrch(1:NumOfCanopyLayers1,1:NumOfBranches_pft(NZ),NZ))
  CanopyStemArea_lpft(1:NumOfCanopyLayers1,1:NZ)=0._r8

  DO NB=1,NumOfBranches_pft(NZ)
    DO L=1,NumOfCanopyLayers1
      CanopyStemArea_lpft(L,NZ)=CanopyStemArea_lpft(L,NZ)+CanopyStemArea_lbrch(L,NB,NZ)
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
          CanopyNodulNonstElms_pft(NE,NZ)=CanopyNodulNonstElms_pft(NE,NZ)+CanopyNodulNonstElms_brch(NE,NB,NZ)
        ENDDO D7950
        CanopyNodulElms_pft(NE,NZ)=sum(CanopyNodulStrutElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))+ &
          sum(CanopyNodulNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
      ENDDO
    ELSEIF(is_root_N2fix(iPlantNfixType(NZ)))THEN
      DO NE=1,NumPlantChemElms
        RootNodulElms_pft(NE,NZ)=sum(RootNodulStrutElms_pvr(NE,NU:MaxSoiL4Root(NZ),NZ))+&
          sum(RootNodulNonstElms_pvr(NE,NU:MaxSoiL4Root(NZ),NZ))
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
  PlantRootSoilElmNetX_pft(1:NumPlantChemElms,NZ)=RootMycoExudElms_pft(1:NumPlantChemElms,NZ)

  PlantRootSoilElmNetX_pft(ielmn,NZ)=PlantRootSoilElmNetX_pft(ielmn,NZ)+&
    RootNH4Uptake_pft(NZ)+RootNO3Uptake_pft(NZ)+RootN2Fix_pft(NZ)
  PlantRootSoilElmNetX_pft(ielmp,NZ)=PlantRootSoilElmNetX_pft(ielmp,NZ)+&
    RootH2PO4Uptake_pft(NZ)+RootHPO4Uptake_pft(NZ)

  PlantExudChemElmCum_pft(1:NumPlantChemElms,NZ)=PlantExudChemElmCum_pft(1:NumPlantChemElms,NZ)+&
    RootMycoExudElms_pft(1:NumPlantChemElms,NZ)
  PlantExudChemElmCum_pft(ielmn,NZ)=PlantExudChemElmCum_pft(ielmn,NZ)+ &
    RootNH4Uptake_pft(NZ)+RootNO3Uptake_pft(NZ)
  PlantExudChemElmCum_pft(ielmp,NZ)=PlantExudChemElmCum_pft(ielmp,NZ)+ &
    RootH2PO4Uptake_pft(NZ)+RootHPO4Uptake_pft(NZ)

  PlantN2FixCum_pft(NZ)=PlantN2FixCum_pft(NZ)+RootN2Fix_pft(NZ)+CanopyN2Fix_pft(NZ)
  end associate
  end subroutine AccumulateStates


end module grosubsMod
