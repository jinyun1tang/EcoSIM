module HfuncsMod
!!
! Description:
! code to do plant phenology

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use PlantAPIData
  use minimathmod, only : AZMAX1
  use StartqsMod, only : StartPlants
  use PlantMathFuncMod
  implicit none

  private


  character(len=*), parameter :: mod_filename = &
  __FILE__

!
! PSIMin4LeafExpansion=minimum canopy turgor potential for leaf expansion (MPa)
! PSILX=minimum canopy water potential for leafout of drought-deciduous PFT (MPa)
! PSILY=minimum canopy water potential for leafoff of drought-deciduous PFT (MPa)
! GSTGG,GSTGR=normalized growth stage durations for vegetative,reproductive phenology
! NBX=maximum branch number for PFT defined by iPlantTurnoverPatternin PFT file
! MaxHour4LeafOutOff=maximum hours for leafout,leafoff
!
  real(r8), PARAMETER :: PSIMin4LeafExpansion=0.1_r8
  real(r8), parameter :: PSILX=-0.2_r8
  real(r8) ,PARAMETER :: GSTGG=2.00_r8
  real(r8), PARAMETER :: GSTGR=0.667_r8
  real(r8), PARAMETER :: MaxHour4LeafOutOff=3600.0_r8
  real(r8), parameter :: PSILY(0:3)=real((/-200.0,-2.0,-2.0,-2.0/),r8)
  integer , parameter :: NBX(0:3)=(/5,1,1,1/)

  public :: hfuncs
  contains

  subroutine hfuncs(I,J)
!
!     THIS subroutine CALCULATES PLANT PHENOLOGY
!
  implicit none

  integer, intent(in) :: I,J
  INTEGER :: NB, NZ, LeafNumberCurrentGrowing

! begin_execution
  associate(                           &
    PSICanopy_pft  =>  plt_ew%PSICanopy_pft    , &
    Hours4ShortenPhotoPeriod    =>  plt_pheno%Hours4ShortenPhotoPeriod   , &
    doInitPlant   =>  plt_pheno%doInitPlant  , &
    doRemobilization_brch   =>  plt_pheno%doRemobilization_brch  , &
    LeafNumberAtFloralInit_brch  =>  plt_pheno%LeafNumberAtFloralInit_brch , &
    Hours4LenthenPhotoPeriod    =>  plt_pheno%Hours4LenthenPhotoPeriod   , &
    IsPlantActive   =>  plt_pheno%IsPlantActive  , &
    iPlantBranchState   =>  plt_pheno%iPlantBranchState  , &
    iPlantCalendar_brch  =>  plt_pheno%iPlantCalendar_brch , &
    HoursCanopyPSITooLow    =>  plt_pheno%HoursCanopyPSITooLow   , &
    KLeafNodeNumber  =>  plt_pheno%KLeafNodeNumber , &
    iPlantMorphologyType  =>  plt_pheno%iPlantMorphologyType , &
    DayLenthCurrent    =>  plt_site%DayLenthCurrent    , &
    DATAP   =>  plt_site%DATAP   , &
    PPT     =>  plt_site%PPT     , &
    DayLenthPrev    =>  plt_site%DayLenthPrev    , &
    NP      =>  plt_site%NP      , &
    PlantPopulation_pft      =>  plt_site%PlantPopulation_pft      , &
    KLeafNumber_brch  =>  plt_morph%KLeafNumber_brch , &
    NumOfLeaves_brch   =>  plt_morph%NumOfLeaves_brch  , &
    NumOfBranches_pft     =>  plt_morph%NumOfBranches_pft    , &
    NumOfMainBranch_pft     =>  plt_morph%NumOfMainBranch_pft      &
  )
  D9985: DO NZ=1,NP

    IF(DATAP(NZ).NE.'NO')THEN
!
!     PPT=total biome population
!
      PPT=PPT+PlantPopulation_pft(NZ)
!
!         SET CROP FLAG ACCORDINGTopRootLayer_pftTO PLANTING, HARVEST DATES, DEATH,
!         1 = ALIVE, 0 = NOT ALIVE
!         DATAP=PFT file name
!
      call set_flags(I,J,NZ)
!
!         INITIALIZE VARIABLES IN ACTIVE PFT
!
      IF(IsPlantActive(NZ).EQ.iPlantIsActive)THEN

        call stage_phenology_vars(I,J,NZ)

        call root_shoot_branching(I,J,NZ)
!
!           THE REST OF THE subroutine MODELS THE PHENOLOGY OF EACH BRANCH
!
!           doInitLeafOut_brch,doPlantLeafOut_brch=flags for initializing leafout,leafoff
!           Hours4Leafout=leafout hours
!
        IF(iPlantCalendar_brch(ipltcal_Emerge,NumOfMainBranch_pft(NZ),NZ).NE.0.OR.doInitPlant(NZ).EQ.itrue)THEN
          
          D2010: DO NB=1,NumOfBranches_pft(NZ)
            IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
              call living_branch_phenology(I,J,NB,nz)
            ENDIF
!
!               KLeafNodeNumber=integer of most recent leaf number currently growing
!
            LeafNumberCurrentGrowing=KLeafNodeNumber(NB,NZ)
            IF(LeafNumberAtFloralInit_brch(NB,NZ).LE.ppmc)THEN
              KLeafNodeNumber(NB,NZ)=INT(NumOfLeaves_brch(NB,NZ))+1
            ELSE
              KLeafNodeNumber(NB,NZ)=INT(AMIN1(NumOfLeaves_brch(NB,NZ),LeafNumberAtFloralInit_brch(NB,NZ)))+1
            ENDIF
            KLeafNumber_brch(NB,NZ)=MIN(MaxNodesPerBranch1-1,KLeafNodeNumber(NB,NZ))
            IF(KLeafNodeNumber(NB,NZ).GT.LeafNumberCurrentGrowing)THEN
              doRemobilization_brch(NB,NZ)=itrue
            ELSE
              doRemobilization_brch(NB,NZ)=ifalse
            ENDIF
!
!               PHENOLOGY
!
!               DayLenthPrev,DLYN=daylength of previous,current day
!               Hours4LenthenPhotoPeriod,Hours4ShortenPhotoPeriod=hourly counter for lengthening,shortening photoperiods
!
            IF(iPlantBranchState(NB,NZ).EQ.iLive.OR.doInitPlant(NZ).EQ.itrue)THEN
              IF(DayLenthCurrent.GE.DayLenthPrev)THEN
                Hours4LenthenPhotoPeriod(NB,NZ)=Hours4LenthenPhotoPeriod(NB,NZ)+1.0_r8
                Hours4ShortenPhotoPeriod(NB,NZ)=0.0_r8
              ELSE
                Hours4LenthenPhotoPeriod(NB,NZ)=0.0_r8
                Hours4ShortenPhotoPeriod(NB,NZ)=Hours4ShortenPhotoPeriod(NB,NZ)+1.0_r8
              ENDIF

              call branch_specific_phenology(I,J,NB,NZ)

            ENDIF
          ENDDO D2010
!
!             WATER STRESS INDICATOR
!
!             PSICanopy_pft=canopy total water potential
!             PSILY=minimum canopy water potential for leafoff
!             HoursCanopyPSITooLow=number of hours PSICanopy_pft(< PSILY (for output only)
!
          IF(PSICanopy_pft(NZ).LT.PSILY(iPlantMorphologyType(NZ)))THEN
            HoursCanopyPSITooLow(NZ)=HoursCanopyPSITooLow(NZ)+1.0_r8
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDDO D9985
  RETURN
  end associate
  END subroutine hfuncs
!------------------------------------------------------------------------------------------

  subroutine set_flags(I,J,NZ)
  use EcoSIMCtrlDataType, only : iYearCurrent

  implicit none
  integer, intent(in) :: I,J,NZ

  INTEGER :: L

! begin_execution
  associate(                        &
    iYearPlanting    =>  plt_distb%iYearPlanting    , &
    iYearPlantHarvest    =>  plt_distb%iYearPlantHarvest    , &  !year of harvest
    iDayPlanting   =>  plt_distb%iDayPlanting   , &  !day of planting
    iDayPlantHarvest   =>  plt_distb%iDayPlantHarvest   , &  !day of harvest
    DATAP   =>  plt_site%DATAP    , &
    IDATA   =>  plt_site%IDATA    , &
    iYearCurrent    =>  plt_site%iYearCurrent     , &
    NumActivePlants   =>  plt_site%NumActivePlants    , &
    Eco_NBP_col    =>  plt_bgcr%Eco_NBP_col     , &
    SeedCPlanted_pft   =>  plt_biom%SeedCPlanted_pft    , &
    iPlantState   =>  plt_pheno%iPlantState   , &
    IsPlantActive   =>  plt_pheno%IsPlantActive     &
  )
  !first hour of day
  IF(J.EQ.1)THEN
    IF(iDayPlanting(NZ).LE.iDayPlantHarvest(NZ) &
      .OR.iYearPlanting(NZ).LT.iYearPlantHarvest(NZ))THEN
      !planting is feasible
      IF(I.GE.iDayPlanting(NZ).OR.iYearCurrent.GT.iYearPlanting(NZ))THEN
        !planted 
        IF(I.GT.iDayPlantHarvest(NZ).AND.iYearCurrent.GE.iYearPlantHarvest(NZ) &
          .AND.iPlantState(NZ).EQ.iDead)THEN
          !post harvest
          IsPlantActive(NZ)=iPlantIsDormant
        ELSE
          IF(I.EQ.iDayPlanting(NZ).AND.iYearCurrent.EQ.iYearPlanting(NZ))THEN
            !planting day
            IsPlantActive(NZ)=iPlantIsDormant
            iPlantState(NZ)=iLive
            CALL StartPlants(NZ,NZ)
            Eco_NBP_col=Eco_NBP_col+SeedCPlanted_pft(NZ)
          ENDIF
          !the living plant has actual properties set
          IF(DATAP(NZ).NE.'NO'.AND.iPlantState(NZ).EQ.iLive)then
            IsPlantActive(NZ)=iPlantIsActive
            print*,'plant is activated',NZ
          endif
        ENDIF
      ELSE
        IsPlantActive(NZ)=iPlantIsDormant
        print*,'plant is dormant',NZ
      ENDIF
    ELSE
      IF((I.LT.iDayPlanting(NZ).AND.I.GT.iDayPlantHarvest(NZ) &
        .AND.iYearCurrent.GE.iYearPlantHarvest(NZ).AND.iPlantState(NZ).EQ.iDead) &
        .OR.(I.LT.iDayPlanting(NZ).AND.iYearPlanting(NZ) &
        .GT.iYearPlantHarvest(NZ)))THEN
        IsPlantActive(NZ)=iPlantIsDormant
      ELSE
        IF(I.EQ.iDayPlanting(NZ).AND.iYearCurrent.EQ.iYearPlanting(NZ))THEN
          IsPlantActive(NZ)=iPlantIsDormant
          iPlantState(NZ)=iLive
          CALL StartPlants(NZ,NZ)
          Eco_NBP_col=Eco_NBP_col+SeedCPlanted_pft(NZ)
        ENDIF
        IF(DATAP(NZ).NE.'NO'.AND.iPlantState(NZ).EQ.iLive)then
          IsPlantActive(NZ)=iPlantIsActive
          print*,'plant is activated',NZ
        endif
      ENDIF
    ENDIF
    NumActivePlants=NumActivePlants+IsPlantActive(NZ)
    print*,'NumActivePlants=',NumActivePlants
  ENDIF
  
  end associate
  end subroutine set_flags
!------------------------------------------------------------------------------------------

  subroutine root_shoot_branching(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: NB
! begin_execution
  associate(                            &
    CanopyNonstructElementConc_pft  =>   plt_biom%CanopyNonstructElementConc_pft  , &
    NonstructalChemElmnts_pft   =>   plt_biom%NonstructalChemElmnts_pft   , &
    MatureGroup_brch  =>   plt_pheno%MatureGroup_brch , &
    iPlantCalendar_brch  =>   plt_pheno%iPlantCalendar_brch , &
    doInitPlant   =>   plt_pheno%doInitPlant  , &
    iPlantRootState_pft   =>   plt_pheno%iPlantRootState_pft  , &
    iPlantPhenologyPattern_pft  =>   plt_pheno%iPlantPhenologyPattern_pft , &
    iPlantBranchState   =>   plt_pheno%iPlantBranchState  , &
    iPlantShootState_pft   =>   plt_pheno%iPlantShootState_pft  , &
    iPlantTurnoverPattern  =>   plt_pheno%iPlantTurnoverPattern , &
    MinNonstructuralC4InitRoot     =>   plt_pheno%MinNonstructuralC4InitRoot    , &
    MatureGroup_pft =>   plt_pheno%MatureGroup_pft, &
    MinNonstructalC4InitBranch     =>   plt_pheno%MinNonstructalC4InitBranch    , &
    PlantPopulation_pft      =>   plt_site%PlantPopulation_pft      , &
    Hours4Leafout    =>   plt_pheno%Hours4Leafout   , &
    PSIRootTurg_vr   =>   plt_ew%PSIRootTurg_vr     , &
    FNOD    =>   plt_allom%FNOD   , &
    NumRootAxes_pft    =>   plt_morph%NumRootAxes_pft   , &
    NumOfMainBranch_pft     =>   plt_morph%NumOfMainBranch_pft    , &
    NumOfBranches_pft     =>   plt_morph%NumOfBranches_pft    , &
    NumConCurrentGrowinNode   =>   plt_morph%NumConCurrentGrowinNode  , &
    BranchNumber_pft     =>   plt_morph%BranchNumber_pft    , &
    BranchNumber_brch    =>   plt_morph%BranchNumber_brch   , &
    NGTopRootLayer_pft     =>   plt_morph%NGTopRootLayer_pft    , &
    XTLI    =>   plt_morph%XTLI   , &
    ShootNodeNumber   =>   plt_morph%ShootNodeNumber    &
  )

!
! ADD BRANCH TO SHOOT IF PLANT GROWTH STAGE, SHOOT NON-STRUCTURAL
! CONCENTRATION PERMIT
!
! doInitPlant=PFT initialization flag:0=no,1=yes
! PSIRootTurg_vr=root turgor potential
! iPlantPhenologyPattern_pft=growth habit from PFT file
! iPlantCalendar_brch(ipltcal_InitFloral,=floral initiation date
! NumOfBranches_pft=primary root axis number
! WTRVC=nonstructural C storage
! PB=nonstructural C concentration needed for branching
! iPlantBranchState=branch life flag:0=living,1=dead
! PSTG=node number
! FNOD=scales node number for perennial vegetation (e.g. trees)
! NumConCurrentGrowinNode=number of concurrently growing nodes
! XTLI,GROUP=node number at planting,floral initiation
! IBTYP: setup for phenologically-driven above-ground turnover


  IF(doInitPlant(NZ).EQ.ifalse)THEN
    !plant initialized
    IF(J.EQ.1.AND.PlantPopulation_pft(NZ).GT.0.0_r8)THEN
      !first hour of the day, population > 0
      IF(PSIRootTurg_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ).GT.PSIMin4LeafExpansion)THEN
        !root is hydraulically active
        IF(iPlantPhenologyPattern_pft(NZ).NE.iplt_annual.OR. &
          iPlantCalendar_brch(ipltcal_InitFloral,NumOfMainBranch_pft(NZ),NZ).EQ.0)THEN
          !perennial plant or flower not initiated for annual plant 
          IF((NumOfBranches_pft(NZ).EQ.0.AND.NonstructalChemElmnts_pft(ielmc,NZ).GT.0.0_r8) &
            .OR.(CanopyNonstructElementConc_pft(ielmc,NZ).GT.MinNonstructalC4InitBranch(NZ) &
            .AND.MinNonstructalC4InitBranch(NZ).GT.0.0_r8))THEN

            D120: DO NB=1,NumOfCanopyLayers1
              IF(iPlantBranchState(NB,NZ).EQ.iDead)THEN
                IF(NB.EQ.NumOfMainBranch_pft(NZ).OR.ShootNodeNumber(NumOfMainBranch_pft(NZ),NZ) &
                  .GT.BranchNumber_pft(NZ)+NumConCurrentGrowinNode(NZ)/FNOD(NZ)+XTLI(NZ))THEN
                  !initiate a new branch
                  BranchNumber_pft(NZ)=BranchNumber_pft(NZ)+1
                  NumOfBranches_pft(NZ)=MIN(NBX(iPlantTurnoverPattern(NZ)),MAX(NB,NumOfBranches_pft(NZ)))
                  BranchNumber_brch(NB,NZ)=BranchNumber_pft(NZ)-1
                  iPlantShootState_pft(NZ)=iLive
                  iPlantBranchState(NB,NZ)=iLive
                  Hours4Leafout(NB,NZ)=0.0_r8
                  IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual)THEN
                    MatureGroup_brch(NB,NZ)=AZMAX1(MatureGroup_pft(NZ)-BranchNumber_brch(NB,NZ))
                  ELSE
                    MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
                  ENDIF
                  exit
                ENDIF
              ENDIF
            ENDDO D120
          ENDIF
        ENDIF
      ENDIF
!
!     ADD AXIS TO ROOT IF PLANT GROWTH STAGE, ROOT NON-STRUCTURAL C
!     CONCENTRATION PERMIT
!
!     PR=nonstructural C concentration needed for root branching
!     FNOD: parameter for allocation of growth to nodes
!     XLI: number of nodes in seed
!     PSTG: node number
!     NumOfMainBranch_pft: number of main branch
!     CanopyNonstructElementConc_pft: canopy nonstructural element concentration
!     PSIRootTurg_vr: root turgor pressure
!     NonstructalChemElmnts_pft: non-structural carbon

      IF(PSIRootTurg_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ).GT.PSIMin4LeafExpansion)THEN
        IF(NumRootAxes_pft(NZ).EQ.0 .OR.ShootNodeNumber(NumOfMainBranch_pft(NZ),NZ).GT.NumRootAxes_pft(NZ)/FNOD(NZ)+XTLI(NZ))THEN
          IF((NumRootAxes_pft(NZ).EQ.0 .AND. NonstructalChemElmnts_pft(ielmc,NZ).GT.0.0_r8) &
            .OR.(CanopyNonstructElementConc_pft(ielmc,NZ).GT.MinNonstructuralC4InitRoot(NZ) & 
            .AND.MinNonstructuralC4InitRoot(NZ).GT.0.0_r8))THEN
            NumRootAxes_pft(NZ)=MIN(NumOfCanopyLayers1,NumRootAxes_pft(NZ)+1)
            iPlantRootState_pft(NZ)=iLive
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine root_shoot_branching
!------------------------------------------------------------------------------------------

  subroutine stage_phenology_vars(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: NB,N,L,NE,BranchNumber_pftX
  real(r8):: ARLSP
  associate(                           &
    CanopyLeafShethC_pft   =>  plt_biom%CanopyLeafShethC_pft     , &
    LeafPetioleBiomassC_brch  =>  plt_biom%LeafPetioleBiomassC_brch    , &
    ShootChemElmnts_pft =>  plt_biom%ShootChemElmnts_pft   , &
    NoduleNonstructCconc_pft =>  plt_biom%NoduleNonstructCconc_pft   , &
    LeafPetioNonstructElmntConc_brch =>  plt_biom%LeafPetioNonstructElmntConc_brch   , &
    CanopyNonstructElementConc_pft =>  plt_biom%CanopyNonstructElementConc_pft   , &
    CanopyNonstructElements_pft =>  plt_biom%CanopyNonstructElements_pft   , &
    NoduleNonstructElmnt_brch =>  plt_biom%NoduleNonstructElmnt_brch   , &
     RootMycoNonstructElmnt_vr =>  plt_biom%RootMycoNonstructElmnt_vr   , &
    NonstructElmnt_brch =>  plt_biom%NonstructElmnt_brch   , &
    RootNonstructElementConcpft_vr =>  plt_biom%RootNonstructElementConcpft_vr   , &
    ZEROL  =>  plt_biom%ZEROL    , &
    ZEROP  =>  plt_biom%ZEROP    , &
    RootStructBiomC_vr =>  plt_biom%RootStructBiomC_vr   , &
    NoduleNonstructElmnt_pft =>  plt_biom%NoduleNonstructElmnt_pft   , &
    WatByPCanopy  =>  plt_ew%WatByPCanopy      , &
    VHeatCapCanP  =>  plt_ew%VHeatCapCanP      , &
    NU     =>  plt_site%NU       , &
    iPlantBranchState  =>  plt_pheno%iPlantBranchState   , &
    iPlantCalendar_brch =>  plt_pheno%iPlantCalendar_brch  , &
    NumOfMainBranch_pft    =>  plt_morph%NumOfMainBranch_pft     , &
    PrimRootDepth  =>  plt_morph%PrimRootDepth   , &
    MY     =>  plt_morph%MY      , &
    CanopyLeafArea_pft  =>  plt_morph%CanopyLeafArea_pft   , &
    NGTopRootLayer_pft    =>  plt_morph%NGTopRootLayer_pft     , &
    NIXBotRootLayer_pft   =>  plt_morph%NIXBotRootLayer_pft     , &
    NumOfBranches_pft    =>  plt_morph%NumOfBranches_pft     , &
    BranchNumber_brch   =>  plt_morph%BranchNumber_brch    , &
    HypoctoylHeight  =>  plt_morph%HypoctoylHeight   , &
    SeedinDepth  =>  plt_morph%SeedinDepth   , &
    CanopyStemA_pft  =>  plt_morph%CanopyStemA_pft   , &
    NI     =>  plt_morph%NI        &
  )
  plt_bgcr%RootGasLoss_disturb(idg_beg:idg_end-1,NZ)=0.0_r8
  CanopyNonstructElements_pft(1:NumOfPlantChemElmnts,NZ)=0.0_r8
  NI(NZ)=NIXBotRootLayer_pft(NZ)
  NGTopRootLayer_pft(NZ)=MIN(NI(NZ),MAX(NGTopRootLayer_pft(NZ),NU))
  NumOfMainBranch_pft(NZ)=1
  BranchNumber_pftX=1.0E+06_r8
!
! TOTAL PLANT NON-STRUCTURAL C, N, P
!
! CPOOL*,ZPOOL*,PPOOL*=non-structl C,N,P in branch(NB),canopy(g)
! CPOLN*,ZPOLN*,PPOLN*=non-structl C,N,P in branch,canopy nodules (g)
! NumOfMainBranch_pft=main branch number
!
  DO NE=1,NumOfPlantChemElmnts
    D140: DO NB=1,NumOfBranches_pft(NZ)
      IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
        CanopyNonstructElements_pft(NE,NZ)=CanopyNonstructElements_pft(NE,NZ)+NonstructElmnt_brch(NE,NB,NZ)
        NoduleNonstructElmnt_pft(NE,NZ)=NoduleNonstructElmnt_pft(NE,NZ)+NoduleNonstructElmnt_brch(NE,NB,NZ)
      ENDIF
    ENDDO D140
  ENDDO

  !find main branch number, which is the most recent live branch
  DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
      print*,'BranchNumber_brch(NB,NZ).LT.BranchNumber_pftX',BranchNumber_brch(NB,NZ),BranchNumber_pftX
      IF(BranchNumber_brch(NB,NZ).LT.BranchNumber_pftX)THEN
        NumOfMainBranch_pft(NZ)=NB
        BranchNumber_pftX=BranchNumber_brch(NB,NZ)
      ENDIF
    ENDIF
  ENDDO
!
! NON-STRUCTURAL C, N, P CONCENTRATIONS IN ROOT
!
! WTRTL=root mas(g)
! CPOOLR,ZPOOLR,PPOOLR=non-structl C,N,P in root(1),myco(2)(g)
! CCPOLR,CZPOLR,CPPOLR=non-structl C,N,P concn in root(1),myco(2)(g g-1)
!
  D180: DO N=1,MY(NZ)
    D160: DO L=NU,NI(NZ)
      IF(RootStructBiomC_vr(N,L,NZ).GT.ZEROL(NZ))THEN
        DO NE=1,NumOfPlantChemElmnts
          RootNonstructElementConcpft_vr(NE,N,L,NZ)=AZMAX1(RootMycoNonstructElmnt_vr(NE,N,L,NZ)/RootStructBiomC_vr(N,L,NZ))
        ENDDO
      ELSE
        DO NE=1,NumOfPlantChemElmnts
          RootNonstructElementConcpft_vr(NE,N,L,NZ)=1.0_r8
        ENDDO
      ENDIF
    ENDDO D160
  ENDDO D180
!
! NON-STRUCTURAL C, N, P CONCENTRATIONS IN SHOOT
!
! CCPOLP,CZPOLP,CPPOLP=nonstructural C,N,P concn in canopy(g g-1)
! NoduleNonstructCconc_pft=nonstructural C concentration in canopy nodules
! CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
!
  IF(CanopyLeafShethC_pft(NZ).GT.ZEROL(NZ))THEN
    DO NE=1,NumOfPlantChemElmnts
      CanopyNonstructElementConc_pft(NE,NZ)=AZMAX1(AMIN1(1.0_r8,CanopyNonstructElements_pft(NE,NZ)/CanopyLeafShethC_pft(NZ)))
    ENDDO
    NoduleNonstructCconc_pft(NZ)=AZMAX1(AMIN1(1.0_r8,NoduleNonstructElmnt_pft(ielmc,NZ)/CanopyLeafShethC_pft(NZ)))
  ELSE
    CanopyNonstructElementConc_pft(1:NumOfPlantChemElmnts,NZ)=1.0_r8
    NoduleNonstructCconc_pft(NZ)=1.0_r8
  ENDIF
  DO NE=1,NumOfPlantChemElmnts
    D190: DO NB=1,NumOfBranches_pft(NZ)
      IF(LeafPetioleBiomassC_brch(NB,NZ).GT.ZEROP(NZ))THEN
        LeafPetioNonstructElmntConc_brch(NE,NB,NZ)=AZMAX1(NonstructElmnt_brch(NE,NB,NZ)/LeafPetioleBiomassC_brch(NB,NZ))
      ELSE
        LeafPetioNonstructElmntConc_brch(NE,NB,NZ)=1.0_r8
      ENDIF
    ENDDO D190
  ENDDO
!
! EMERGENCE DATE FROM COTYLEDON HEIGHT, LEAF AREA, ROOT DEPTH
!
! iPlantCalendar_brch(ipltcal_Emerge,=emergence date
! CanopyLeafArea_pft,CanopyStemA_pft=leaf,stalk areas
! HypoctoylHeight=hypocotyledon height
! SeedinDepth=seeding depth
! PrimRootDepth=primary root depth
! VHeatCapCanP,WTSHT,WatByPCanopy=canopy heat capacity,mass,water content
!
  IF(iPlantCalendar_brch(ipltcal_Emerge,NumOfMainBranch_pft(NZ),NZ).EQ.0)THEN
    ARLSP=CanopyLeafArea_pft(NZ)+CanopyStemA_pft(NZ)
    IF((HypoctoylHeight(NZ).GT.SeedinDepth(NZ)).AND.(ARLSP.GT.ZEROL(NZ)) &
      .AND.(PrimRootDepth(1,1,NZ).GT.SeedinDepth(NZ)+ppmc))THEN
      iPlantCalendar_brch(ipltcal_Emerge,NumOfMainBranch_pft(NZ),NZ)=I
      VHeatCapCanP(NZ)=cpw*(ShootChemElmnts_pft(ielmc,NZ)*10.0E-06_r8+WatByPCanopy(NZ))
      write(101,*)'emergence',I,NumOfMainBranch_pft(NZ),NZ
    ENDIF
  ENDIF
  end associate
  end subroutine stage_phenology_vars
!------------------------------------------------------------------------------------------

  subroutine branch_specific_phenology(I,J,NB,NZ)
  implicit none
  integer, intent(in) :: I,J,NB,NZ

  associate(                           &
    PSICanopy_pft  =>  plt_ew%PSICanopy_pft    , &
    PSICanPTurg   =>  plt_ew%PSICanPTurg     , &
    DayLenthCurrent    =>  plt_site%DayLenthCurrent    , &
    DayLenthPrev    =>  plt_site%DayLenthPrev    , &
    DayLenthMax    =>  plt_site%DayLenthMax    , &
    ALAT    =>  plt_site%ALAT    , &
    doPlantLeafOut_brch   =>  plt_pheno%doPlantLeafOut_brch  , &
    iPlantCalendar_brch  =>  plt_pheno%iPlantCalendar_brch , &
    iPlantMorphologyType  =>  plt_pheno%iPlantMorphologyType , &
    Hours4ShortenPhotoPeriod    =>  plt_pheno%Hours4ShortenPhotoPeriod   , &
    Hours4Leafout    =>  plt_pheno%Hours4Leafout   , &
    HourThreshold4LeafOut_brch   =>  plt_pheno%HourThreshold4LeafOut_brch  , &
    iPlantPhenologyType_pft  =>  plt_pheno%iPlantPhenologyType_pft , &
    TCelciusChill4Leaf    =>  plt_pheno%TCelciusChill4Leaf   , &
    Hours4LenthenPhotoPeriod    =>  plt_pheno%Hours4LenthenPhotoPeriod   , &
    TCG     =>  plt_pheno%TCG    , &
    TCelciusChill4Seed    =>  plt_pheno%TCelciusChill4Seed   , &
    TCelcius4LeafOffHarden     =>  plt_pheno%TCelcius4LeafOffHarden    , &
    Hours4LeafOff    =>  plt_pheno%Hours4LeafOff   , &
    HourThreshold4LeafOff   =>  plt_pheno%HourThreshold4LeafOff  , &
    doPlantLeaveOff_brch   =>  plt_pheno%doPlantLeaveOff_brch    &
  )
!
! CALCULATE EVERGREEN PHENOLOGY DURINGTopRootLayer_pftLENGTHENINGTopRootLayer_pftPHOTOPERIODS
!
! iPlantPhenologyType_pft=phenology type from PFT file
! DayLenthPrev,DLYN=daylength of previous,current day
! Hours4Leafout,Hours4LeafOff=leafout,leafoff hours
! Hours4LenthenPhotoPeriod=hourly counter for lengthening photoperiods
! doPlantLeaveOff_brch=flag for enabling leafoff:0=enable,1=disable
! ALAT=latitude
!
  IF(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_evgreen)THEN
    IF(DayLenthCurrent.GE.DayLenthPrev)THEN
      Hours4Leafout(NB,NZ)=Hours4LenthenPhotoPeriod(NB,NZ)
      IF(Hours4Leafout(NB,NZ).GE.HourThreshold4LeafOut_brch(NB,NZ) &
        .OR.(ALAT.GT.0.0_r8.AND.I.EQ.173) &
        .OR.(ALAT.LT.0.0_r8.AND.I.EQ.355))THEN
        Hours4LeafOff(NB,NZ)=0.0_r8
        doPlantLeaveOff_brch(NB,NZ)=iEnable
      ENDIF
    ENDIF
!
!   CALCULATE EVERGREEN PHENOLOGY DURINGTopRootLayer_pftSHORTENINGTopRootLayer_pftPHOTOPERIODS
!
!   Hours4Leafout,Hours4LeafOff=leafout,leafoff hours
!   Hours4ShortenPhotoPeriod=hourly counter for shortening photoperiods
!   doPlantLeafOut_brch=flag for enabling leafout:0=enable,1=disable
!   ALAT=latitude
!
    IF(DayLenthCurrent.LT.DayLenthPrev)THEN
      Hours4LeafOff(NB,NZ)=Hours4ShortenPhotoPeriod(NB,NZ)
      IF(Hours4LeafOff(NB,NZ).GE.HourThreshold4LeafOff(NB,NZ) &
        .OR.(ALAT.GT.0.0_r8.AND.I.EQ.355) &
        .OR.(ALAT.LT.0.0_r8.AND.I.EQ.173))THEN
        Hours4Leafout(NB,NZ)=0.0_r8
        doPlantLeafOut_brch(NB,NZ)=iEnable
      ENDIF
    ENDIF
!
!   CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS ABOVE
!   SPECIFIED TEMPERATURE DURINGTopRootLayer_pftLENGTHENINGTopRootLayer_pftPHOTOPERIODS
!
!   iPlantPhenologyType_pft=phenology type from PFT file
!   DayLenthPrev,DLYN=daylength of previous,current day
!   Hours4Leafout,VRNL=leafout hours,hours required for leafout
!   Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!   doPlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
!   TCG,TCZ,TCelciusChill4Seed=canopy temp,leafout threshold temp,chilling temp
!   ALAT=latitude
!   iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
!
  ELSEIF(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldecidu)THEN
    IF((DayLenthCurrent.GE.DayLenthPrev.OR.&
      (DayLenthCurrent.LT.DayLenthPrev.AND.Hours4LeafOff(NB,NZ).LT.HourThreshold4LeafOff(NB,NZ))) &
      .AND.doPlantLeafOut_brch(NB,NZ).EQ.iEnable)THEN
      IF(TCG(NZ).GE.TCelciusChill4Leaf(NZ))THEN
        Hours4Leafout(NB,NZ)=Hours4Leafout(NB,NZ)+1.0_r8
      ENDIF
      IF(Hours4Leafout(NB,NZ).LT.HourThreshold4LeafOut_brch(NB,NZ))THEN
        IF(TCG(NZ).LT.TCelciusChill4Seed(NZ))THEN
          Hours4Leafout(NB,NZ)=AZMAX1(Hours4Leafout(NB,NZ)-1.0)
        ENDIF
      ENDIF
      IF(Hours4Leafout(NB,NZ).GE.HourThreshold4LeafOut_brch(NB,NZ) &
        .OR.(ALAT.GT.0.0.AND.I.EQ.173) &
        .OR.(ALAT.LT.0.0.AND.I.EQ.355))THEN
        Hours4LeafOff(NB,NZ)=0.0
      ENDIF
    ENDIF

    IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0.OR. &
      (DayLenthCurrent.LT.DayLenthPrev.AND.DayLenthCurrent.LT.12.0_r8))THEN
      doPlantLeaveOff_brch(NB,NZ)=iEnable
    ENDIF
!
!     CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS BELOW
!     SPECIFIED TEMPERATURE DURINGTopRootLayer_pftSHORTENINGTopRootLayer_pftPHOTOPERIODS
!
!     DayLenthPrev,DLYN=daylength of previous,current day
!     Hours4Leafout,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!     doPlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,TCelciusChill4Seed=canopy temp,leafout threshold temp,chilling temp
!     ALAT=latitude
!     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
!
    IF(DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeaveOff_brch(NB,NZ).EQ.iEnable &
      .AND.iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)THEN
      IF(TCG(NZ).LE.TCelcius4LeafOffHarden(NZ))THEN
        Hours4LeafOff(NB,NZ)=Hours4LeafOff(NB,NZ)+1.0_r8
      ENDIF
      IF(Hours4LeafOff(NB,NZ).GE.HourThreshold4LeafOff(NB,NZ).AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable)THEN
        Hours4Leafout(NB,NZ)=0.0
        doPlantLeafOut_brch(NB,NZ)=iEnable
      ENDIF
    ENDIF

!
!     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS
!     ABOVE SPECIFIED WATER POTENTIAL DURINGTopRootLayer_pftDORMANCY
!
!     iPlantPhenologyType_pft=phenology type from PFT file
!     Hours4Leafout,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff=leafoff hours
!     doPlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
!     PSICanopy_pft=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
!
  ELSEIF(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_drouhtdecidu.OR. &
    iPlantPhenologyType_pft(NZ).EQ.4.OR.iPlantPhenologyType_pft(NZ).EQ.5)THEN
    !type 4 and 5 are place holders for subtropical/tropical evergreen
    !
    IF(doPlantLeafOut_brch(NB,NZ).EQ.iEnable)THEN
      IF(PSICanopy_pft(NZ).GE.PSILX)THEN
        Hours4Leafout(NB,NZ)=Hours4Leafout(NB,NZ)+1.0_r8
      ENDIF
      IF(Hours4Leafout(NB,NZ).LT.HourThreshold4LeafOut_brch(NB,NZ))THEN
        IF(PSICanopy_pft(NZ).LT.PSILY(iPlantMorphologyType(NZ)))THEN
          Hours4Leafout(NB,NZ)=AZMAX1(Hours4Leafout(NB,NZ)-12.0)
        ENDIF
      ENDIF
      IF(Hours4Leafout(NB,NZ).GE.HourThreshold4LeafOut_brch(NB,NZ))THEN
        Hours4LeafOff(NB,NZ)=0.0_r8
        IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iEnable
      ENDIF
    ENDIF
    IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iEnable
!
!     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayer_pftHOURS
!     BELOW SPECIFIED WATER POTENTIAL DURINGTopRootLayer_pftGROWINGTopRootLayer_pftSEASON
!
!     Hours4Leafout=leafout hours,hours required for leafout
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!     doPlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
!     PSICanopy_pft=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
!     Hours4LenthenPhotoPeriod,Hours4ShortenPhotoPeriod=hourly counter for lengthening,shortening photoperiods
!     MaxHour4LeafOutOff=maximum hours for leafout,leafoff
!
      IF(doPlantLeafOut_brch(NB,NZ).EQ.iDisable.AND.doPlantLeaveOff_brch(NB,NZ).EQ.iEnable)THEN
        IF(PSICanopy_pft(NZ).LT.PSILY(iPlantMorphologyType(NZ)))THEN
          Hours4LeafOff(NB,NZ)=Hours4LeafOff(NB,NZ)+1.0_r8
        ENDIF
        IF(iPlantPhenologyType_pft(NZ).EQ.4)THEN
          IF(Hours4ShortenPhotoPeriod(NB,NZ).GT.MaxHour4LeafOutOff)THEN
            Hours4LeafOff(NB,NZ)=Hours4ShortenPhotoPeriod(NB,NZ)
          ENDIF
        ELSEIF(iPlantPhenologyType_pft(NZ).EQ.5)THEN
          IF(Hours4LenthenPhotoPeriod(NB,NZ).GT.MaxHour4LeafOutOff)THEN
            Hours4LeafOff(NB,NZ)=Hours4LenthenPhotoPeriod(NB,NZ)
          ENDIF
        ENDIF
        IF(Hours4LeafOff(NB,NZ).GE.HourThreshold4LeafOff(NB,NZ).AND. &
          doPlantLeafOut_brch(NB,NZ).EQ.iDisable)THEN
          Hours4Leafout(NB,NZ)=0.0_r8
          doPlantLeafOut_brch(NB,NZ)=iEnable
        ENDIF
      ENDIF
!
!     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
!     HOURS ABOVE SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
!     LENGTHENINGTopRootLayer_pftPHOTOPERIODS
!
!     iPlantPhenologyType_pft=phenology type from PFT file
!     DayLenthPrev,DLYN=daylength of previous,current day
!     Hours4Leafout,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!     PSICanopy_pft=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     doPlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,TCelciusChill4Seed=canopy temp,leafout threshold temp,chilling temp
!     ALAT=latitude
!     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
!
    ELSEIF(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldroutdecidu)THEN
      IF((DayLenthCurrent.GE.DayLenthPrev.OR.DayLenthCurrent.GE.DayLenthMax-2.0_r8) &
        .AND.doPlantLeafOut_brch(NB,NZ).EQ.iEnable)THEN
        IF(TCG(NZ).GE.TCelciusChill4Leaf(NZ).AND.PSICanPTurg(NZ).GT.PSIMin4LeafExpansion)THEN
          Hours4Leafout(NB,NZ)=Hours4Leafout(NB,NZ)+1.0_r8
        ENDIF
        IF(Hours4Leafout(NB,NZ).LT.HourThreshold4LeafOut_brch(NB,NZ))THEN
          IF(TCG(NZ).LT.TCelciusChill4Seed(NZ).OR.PSICanPTurg(NZ).LT.PSIMin4LeafExpansion)THEN
            Hours4Leafout(NB,NZ)=AZMAX1(Hours4Leafout(NB,NZ)-1.5_r8)
          ENDIF
        ENDIF
        IF(Hours4Leafout(NB,NZ).GE.HourThreshold4LeafOut_brch(NB,NZ))THEN
          Hours4LeafOff(NB,NZ)=0.0_r8
          IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iEnable
        ENDIF
      ENDIF
      IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff_brch(NB,NZ)=iEnable
!
!     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
!     HOURS BELOW SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
!     SHORTENINGTopRootLayer_pftPHOTOPERIODS
!
!     DayLenthPrev,DLYN=daylength of previous,current day
!     Hours4Leafout,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!     doPlantLeafOut_brch,doPlantLeaveOff_brch=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,TCelciusChill4Seed=canopy temp,leafout threshold temp,chilling temp
!     PSICanopy_pft=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
!
      IF((DayLenthCurrent.LT.DayLenthPrev.OR.DayLenthCurrent.LT.24.0_r8-DayLenthMax+2.0_r8) &
        .AND.doPlantLeaveOff_brch(NB,NZ).EQ.iEnable)THEN
        IF(TCG(NZ).LE.TCelcius4LeafOffHarden(NZ) &
          .OR.PSICanopy_pft(NZ).LT.PSILY(iPlantMorphologyType(NZ)))THEN
          Hours4LeafOff(NB,NZ)=Hours4LeafOff(NB,NZ)+1.0_r8
        ENDIF
        IF(Hours4LeafOff(NB,NZ).GE.HourThreshold4LeafOff(NB,NZ) &
          .AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable)THEN
          Hours4Leafout(NB,NZ)=0.0_r8
          doPlantLeafOut_brch(NB,NZ)=iEnable
        ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine branch_specific_phenology
!------------------------------------------------------------------------------------------

  subroutine living_branch_phenology(I,J,NB,NZ)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NB  !plant branch id
  integer, intent(in) :: NZ  !plant species id

  real(r8) :: TFNP,WFNG
  real(r8) :: ACTV,OFNG
  real(r8) :: PPD
  real(r8) :: RTK
  real(r8) :: STK,TKCO
  real(r8) :: NodeInitRate,LeafAppearRate

! begin_execution
  associate(                        &
    iYearPlanting    =>  plt_distb%iYearPlanting    , &
    iDayPlanting   =>  plt_distb%iDayPlanting   , &
    SnowDepth   =>  plt_ew%SnowDepth      , &
    PSICanopy_pft  =>  plt_ew%PSICanopy_pft     , &
    ZERO    =>  plt_site%ZERO     , &
    DayLenthCurrent    =>  plt_site%DayLenthCurrent     , &
    iYearCurrent    =>  plt_site%iYearCurrent     , &
    DayLenthPrev    =>  plt_site%DayLenthPrev     , &
    TKG     =>  plt_pheno%TKG     , &
    iPlantCalendar_brch  =>  plt_pheno%iPlantCalendar_brch  , &
    RefLeafAppearRate    =>  plt_pheno%RefLeafAppearRate    , &
    ReprodNodeNumNormByMatrgrp_brch   =>  plt_pheno%ReprodNodeNumNormByMatrgrp_brch   , &
    OFFST   =>  plt_pheno%OFFST   , &
    iPlantPhenologyPattern_pft  =>  plt_pheno%iPlantPhenologyPattern_pft  , &
    doPlantLeafOut_brch   =>  plt_pheno%doPlantLeafOut_brch   , &
    NumOfLeaves_brch   =>  plt_morph%NumOfLeaves_brch   , &
    NodeNumNormByMatgrp_brch   =>  plt_pheno%NodeNumNormByMatgrp_brch   , &
    HourReprodNodeNumNormByMatrgrp_brch  =>  plt_pheno%HourReprodNodeNumNormByMatrgrp_brch  , &
    doSenescence_brch   =>  plt_pheno%doSenescence_brch   , &
    LeafNumberAtFloralInit_brch  =>  plt_pheno%LeafNumberAtFloralInit_brch  , &
    iPlantDevelopPattern  =>  plt_pheno%iPlantDevelopPattern  , &
    TotalReprodNodeNumNormByMatrgrp_brch  =>  plt_pheno%TotalReprodNodeNumNormByMatrgrp_brch  , &
    CriticalPhotoPeriod_pft    =>  plt_pheno%CriticalPhotoPeriod_pft    , &
    HourlyNodeNumNormByMatgrp_brch  =>  plt_pheno%HourlyNodeNumNormByMatgrp_brch  , &
    PhotoPeriodSens_pft   =>  plt_pheno%PhotoPeriodSens_pft   , &
    TotalNodeNumNormByMatgrp_brch  =>  plt_pheno%TotalNodeNumNormByMatgrp_brch  , &
    MatureGroup_brch  =>  plt_pheno%MatureGroup_brch  , &
    iPlantPhenologyType_pft  =>  plt_pheno%iPlantPhenologyType_pft  , &
    HourThreshold4LeafOut_brch   =>  plt_pheno%HourThreshold4LeafOut_brch   , &
    iPlantPhotoperiodType  =>  plt_pheno%iPlantPhotoperiodType  , &
    doInitLeafOut_brch   =>  plt_pheno%doInitLeafOut_brch   , &
    MatureGroup_pft =>  plt_pheno%MatureGroup_pft , &
    PlantO2Stress    =>  plt_pheno%PlantO2Stress    , &
    Hours4Leafout    =>  plt_pheno%Hours4Leafout    , &
    Hours4LeafOff    =>  plt_pheno%Hours4LeafOff    , &
    HourThreshold4LeafOff   =>  plt_pheno%HourThreshold4LeafOff   , &
    RefNodeInitRate   =>  plt_pheno%RefNodeInitRate   , &
    CanopyHeight_pft     =>  plt_morph%CanopyHeight_pft     , &
    ShootNodeNumber   =>   plt_morph%ShootNodeNumber  , &
    NodeNumberToInitFloral   =>   plt_morph%NodeNumberToInitFloral  , &
    NodeNumberAtAnthesis   =>   plt_morph%NodeNumberAtAnthesis  , &
    NumOfMainBranch_pft     =>   plt_morph%NumOfMainBranch_pft      &
  )
  IF(iPlantCalendar_brch(ipltcal_Emerge,NB,NZ).EQ.0)THEN
    !plant emergence
    iPlantCalendar_brch(ipltcal_Emerge,NB,NZ)=I
    doInitLeafOut_brch(NB,NZ)=iDisable
    doPlantLeafOut_brch(NB,NZ)=iEnable
    Hours4Leafout(NB,NZ)=0.5_r8*Hours4Leafout(NumOfMainBranch_pft(NZ),NZ)
  ENDIF
!
! CALCULATE NODE INITIATION AND LEAF APPEARANCE RATES
! FROM TEMPERATURE FUNCTION CALCULATED IN 'UPTAKE' AND
! RATES AT 25C ENTERED IN 'READQ' EXCEPT WHEN DORMANT
!
! iPlantPhenologyType_pft=phenology type from PFT file, 0=evergreen,1=cold deciduous, 2=drought deciduous,3=1+2
! Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
! TKG,TKCO=canopy temperature,canopy temp used in Arrhenius eqn
! OFFST=shift in Arrhenius curve for thermal adaptation
! TFNP=temperature function for phenology (25 oC =1 )
! 8.313,710.0=gas constant,enthalpy
! 60000,197500,218500=energy of activn,high,low temp inactivn(KJ mol-1)
! NodeInitRate,LeafAppearRate=rates of node initiation,leaf appearance
! XRNI,XRLA=rate of node initiation,leaf appearance at 25 oC (h-1)
!
  IF(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_evgreen.OR.Hours4LeafOff(NB,NZ).LT.HourThreshold4LeafOff(NB,NZ))THEN
    TKCO=TKG(NZ)+OFFST(NZ)
    TFNP=calc_leave_grow_tempf(TKCO)    
    NodeInitRate=AZMAX1(RefNodeInitRate(NZ)*TFNP)
    LeafAppearRate=AZMAX1(RefLeafAppearRate(NZ)*TFNP)
!
!   NODE INITIATION AND LEAF APPEARANCE RATES SLOWED BY LOW TURGOR
!
!   PSICanPTurg=leaf turgor potential
!   WFNG=water stress effect on phenology
!   only annual plants depends on moisture
    IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual)THEN
      IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
        WFNG=EXP(0.025_r8*PSICanopy_pft(NZ))
        NodeInitRate=NodeInitRate*WFNG
        LeafAppearRate=LeafAppearRate*WFNG
      ENDIF
      IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
        OFNG=SQRT(PlantO2Stress(NZ))
        NodeInitRate=NodeInitRate*OFNG
        LeafAppearRate=LeafAppearRate*OFNG
      ENDIF
    ENDIF
!
!   ACCUMULATE NODE INITIATION AND LEAF APPEARANCE RATES
!   INTO TOTAL NUMBER OF NODES AND LEAVES
!
!   PSTG,NumOfLeaves_brch=number of nodes initiated,leaves appeared
!
    ShootNodeNumber(NB,NZ)=ShootNodeNumber(NB,NZ)+NodeInitRate
    NumOfLeaves_brch(NB,NZ)=NumOfLeaves_brch(NB,NZ)+LeafAppearRate
!
!   USE TOTAL NUMBER OF NODES TO CALCULATE PROGRESSION THROUGH
!   VEGETATIVE AND REPRODUCTIVE GROWTH STAGES. THIS PROGRESSION
!   IS USED TO SET START AND END DATES FOR GROWTH STAGES BELOW
!
!   NodeNumNormByMatgrp_brch=vegetative node number normalized for maturity group
!   MatureGroup_pft=node number required for floral initiation
!   HourlyNodeNumNormByMatgrp_brch,TotalNodeNumNormByMatgrp_brch=hourly,total change in NodeNumNormByMatgrp_brch
!   ReprodNodeNumNormByMatrgrp_brch=reproductive node number normalized for maturity group
!   NodeNumberAtAnthesis=node number at flowering
!   HourReprodNodeNumNormByMatrgrp_brch,TotalReprodNodeNumNormByMatrgrp_brch=hourly,total change in ReprodNodeNumNormByMatrgrp_brch
!   doSenescence_brch=PFT senescence flag
!
    IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).NE.0)THEN
      NodeNumNormByMatgrp_brch(NB,NZ)=(ShootNodeNumber(NB,NZ)-NodeNumberToInitFloral(NB,NZ))/MatureGroup_pft(NZ)
      HourlyNodeNumNormByMatgrp_brch(NB,NZ)=NodeInitRate/(MatureGroup_pft(NZ)*GSTGG)
      TotalNodeNumNormByMatgrp_brch(NB,NZ)=TotalNodeNumNormByMatgrp_brch(NB,NZ)+HourlyNodeNumNormByMatgrp_brch(NB,NZ)
    ENDIF
    IF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).NE.0)THEN
      ReprodNodeNumNormByMatrgrp_brch(NB,NZ)=(ShootNodeNumber(NB,NZ)-NodeNumberAtAnthesis(NB,NZ))/MatureGroup_pft(NZ)
      HourReprodNodeNumNormByMatrgrp_brch(NB,NZ)=NodeInitRate/(MatureGroup_pft(NZ)*GSTGR)
      TotalReprodNodeNumNormByMatrgrp_brch(NB,NZ)=TotalReprodNodeNumNormByMatrgrp_brch(NB,NZ) &
        +HourReprodNodeNumNormByMatrgrp_brch(NB,NZ)
    ENDIF
    doSenescence_brch(NB,NZ)=itrue
  ELSE
    doSenescence_brch(NB,NZ)=ifalse
  ENDIF
!
! REPRODUCTIVE GROWTH STAGES ADVANCE WHEN THRESHOLD NUMBER
! OF NODES HAVE BEEN INITIATED. FIRST DETERMINE PHOTOPERIOD
! AND TEMPERATURE EFFECTS ON FINAL VEG NODE NUMBER FROM
! NUMBER OF INITIATED NODES
!
! PSTG=number of nodes initiated
! NodeNumberToInitFloral=node number at floral initiation
! GROUP=node number required for floral initiation
! Hours4Leafout,VRNL=leafout hours,hours required for leafout
! DayLenthPrev,DLYN=daylength of previous,current day
! iPlantPhenologyPattern_pft=growth habit from PFT file
! iPlantPhenologyType_pft=phenology type from PFT file
! ZC,SnowDepth=canopy height,snowpack depth
!
  IF(iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
    IF(ShootNodeNumber(NB,NZ).GT.MatureGroup_brch(NB,NZ)+NodeNumberToInitFloral(NB,NZ) &
      .AND.((Hours4Leafout(NB,NZ).GE.HourThreshold4LeafOut_brch(NB,NZ)) &
      .OR.(I.GE.iDayPlanting(NZ).AND.iYearCurrent.EQ.iYearPlanting(NZ) &
      .AND.DayLenthCurrent.GT.DayLenthPrev)) &
      .OR.(((iPlantPhenologyPattern_pft(NZ).EQ.iplt_perennial.AND.(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldecidu &
      .OR.iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldroutdecidu)) &
      .OR.(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenologyType_pft(NZ).EQ.iphenotyp_evgreen)) &
      .AND.CanopyHeight_pft(NZ).GE.SnowDepth-ZERO &
      .AND.DayLenthCurrent.LT.DayLenthPrev))THEN
!
!     FINAL VEGETATIVE NODE NUMBER DEPENDS ON PHOTOPERIOD FROM 'DAY'
!     AND ON MATURITY GROUP, CRITICAL PHOTOPERIOD AND PHOTOPERIOD
!     SENSITIVITY ENTERED IN 'READQ'
!
!     iPlantPhotoperiodType=photoperiod type from PFT file
!     PPD=photoperiod sensitivity
!     CriticalPhotoPeriod_pft=critical photoperiod from PFT file
!     iPlantCalendar_brch(ipltcal_InitFloral,=date of floral initiation
!     VSTGX=node number on date of floral initiation
!
      IF(iPlantPhotoperiodType(NZ).EQ.iphotop_neutral)THEN
        PPD=0.0_r8
      ELSE
        PPD=AZMAX1(CriticalPhotoPeriod_pft(NZ)-DayLenthCurrent)
        IF(iPlantPhotoperiodType(NZ).EQ.iphotop_short.AND.DayLenthCurrent.GE.DayLenthPrev)PPD=0.0_r8
      ENDIF

      IF(iPlantPhotoperiodType(NZ).EQ.iphotop_neutral &
        .OR.(iPlantPhotoperiodType(NZ).EQ.iphotop_short.AND.PPD.GT.PhotoPeriodSens_pft(NZ)) &
        .OR.(iPlantPhotoperiodType(NZ).EQ.iphotop_long.AND.PPD.LT.PhotoPeriodSens_pft(NZ)) &
        .OR.(((iPlantPhenologyPattern_pft(NZ).EQ.iplt_perennial.AND.(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldecidu &
        .OR.iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldroutdecidu)) &
        .OR.(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual.AND.iPlantPhenologyType_pft(NZ).EQ.iphenotyp_evgreen)) &
        .AND.CanopyHeight_pft(NZ).GE.SnowDepth-ZERO &
        .AND.DayLenthCurrent.LT.DayLenthPrev))THEN
        write(101,*)'plant init floral',I,NB,NZ
        iPlantCalendar_brch(ipltcal_InitFloral,NB,NZ)=I
        NodeNumberToInitFloral(NB,NZ)=ShootNodeNumber(NB,NZ)
        IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual.AND.iPlantDevelopPattern(NZ).EQ.ideterminate)THEN
          LeafNumberAtFloralInit_brch(NB,NZ)=ShootNodeNumber(NB,NZ)
        ENDIF
      ENDIF
    ENDIF
!
!   STEM ELONGATION
!
!   NodeNumNormByMatgrp_brch=vegetative node number normalized for maturity group
!   GSTGG=normalized growth stage durations for vegetative phenology
!   iPlantCalendar_brch(ipltcal_Jointing,=start of stem elongation and setting max seed number
!
  ELSEIF(iPlantCalendar_brch(ipltcal_Jointing,NB,NZ).EQ.0)THEN
    IF(NodeNumNormByMatgrp_brch(NB,NZ).GT.0.25_r8*GSTGG &
      .OR.((iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldecidu.OR.iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldroutdecidu) &
      .AND.iPlantPhenologyPattern_pft(NZ).NE.iplt_annual.AND.iPlantPhotoperiodType(NZ).NE.1 &
      .AND.DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ)) &
      .OR.(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_drouhtdecidu.AND.iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual) &
      .AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ))THEN
      iPlantCalendar_brch(ipltcal_Jointing,NB,NZ)=I
      write(101,*)'plant jointing',I,NB,NZ
    ENDIF
!
!   iPlantCalendar_brch(ipltcal_Elongation,=mid stem elongation
!
  ELSEIF(iPlantCalendar_brch(ipltcal_Elongation,NB,NZ).EQ.0)THEN
    IF(NodeNumNormByMatgrp_brch(NB,NZ).GT.0.50_r8*GSTGG &
      .OR.((iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldecidu.OR.iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldroutdecidu) &
      .AND.iPlantPhenologyPattern_pft(NZ).NE.iplt_annual.AND.iPlantPhotoperiodType(NZ).NE.1 &
      .AND.DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ)) &
      .OR.(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_drouhtdecidu.AND.iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual) &
      .AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ))THEN
      iPlantCalendar_brch(ipltcal_Elongation,NB,NZ)=I
      write(101,*)'plant elongation',I,NB,NZ
      IF(iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual.AND.iPlantDevelopPattern(NZ).NE.ideterminate)THEN
        LeafNumberAtFloralInit_brch(NB,NZ)=ShootNodeNumber(NB,NZ)
      ENDIF
    ENDIF
!
!   iPlantCalendar_brch(ipltcal_Heading,=end of stem elongation and setting max seed number
!
  ELSEIF(iPlantCalendar_brch(ipltcal_Heading,NB,NZ).EQ.0)THEN
    IF(NodeNumNormByMatgrp_brch(NB,NZ).GT.1.00_r8*GSTGG &
      .OR.((iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldecidu.OR.iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldroutdecidu) &
      .AND.iPlantPhenologyPattern_pft(NZ).NE.iplt_annual.AND.iPlantPhotoperiodType(NZ).NE.1 &
      .AND.DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ)) &
      .OR.(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_drouhtdecidu.AND.iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual) &
      .AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ))THEN
      iPlantCalendar_brch(ipltcal_Heading,NB,NZ)=I
      write(101,*)'plant heading',I,NB,NZ
    ENDIF
!
!   ANTHESIS OCCURS WHEN THE NUMBER OF LEAVES THAT HAVE APPEARED
!   EQUALS THE NUMBER OF NODES INITIATED WHEN THE FINAL VEGETATIVE
!   NODE NUMBER WAS SET ABOVE
!
!   iPlantCalendar_brch(ipltcal_Anthesis,=start of anthesis and setting final seed number
!   NumOfLeaves_brch=number of leaves appeared
!   NodeNumberToInitFloral=node number at floral initiation
!   ISTYP,IWTYP,iPlantPhotoperiodType=growth habit,phenology,photoperiod type from PFT file
!   doPlantLeafOut_brch=flag for enabling leafout:0=enable,1=disable
!   Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!   DayLenthPrev,DLYN=daylength of previous,current day
!   NodeNumberAtAnthesis=number of nodes at anthesis
!
  ELSEIF(iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
    IF((NumOfLeaves_brch(NB,NZ).GT.NodeNumberToInitFloral(NB,NZ)) &
      .OR.(iPlantPhenologyPattern_pft(NZ).NE.iplt_annual.AND.iPlantCalendar_brch(ipltcal_Heading,NB,NZ).NE.0) &
      .OR.((iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldecidu.OR.iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldroutdecidu) &
      .AND.iPlantPhenologyPattern_pft(NZ).NE.iplt_annual.AND.iPlantPhotoperiodType(NZ).NE.1 &
      .AND.DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ)) &
      .OR.(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_drouhtdecidu.AND.iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual) &
      .AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ))THEN
      IF(NB.EQ.NumOfMainBranch_pft(NZ).OR.iPlantCalendar_brch(ipltcal_Anthesis,NumOfMainBranch_pft(NZ),NZ).NE.0)THEN
        iPlantCalendar_brch(ipltcal_Anthesis,NB,NZ)=I
        NodeNumberAtAnthesis(NB,NZ)=ShootNodeNumber(NB,NZ)
        write(101,*)'plant anthesis',I,NB,NZ
      ENDIF
    ENDIF
!
!   START GRAIN FILL PERIOD
!
!   iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
!   ReprodNodeNumNormByMatrgrp_brch=reproductive node number normalized for maturity group
!   GSTGR=normalized growth stage durations for reproductive phenology
!
!
  ELSEIF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
    IF(ReprodNodeNumNormByMatrgrp_brch(NB,NZ).GT.0.50_r8*GSTGR &
      .OR.((iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldecidu.OR.iPlantPhenologyType_pft(NZ).EQ.iphenotyp_coldroutdecidu) &
      .AND.iPlantPhenologyPattern_pft(NZ).NE.iplt_annual.AND.iPlantPhotoperiodType(NZ).NE.1 &
      .AND.DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ)) &
      .OR.(iPlantPhenologyType_pft(NZ).EQ.iphenotyp_drouhtdecidu.AND.iPlantPhenologyPattern_pft(NZ).EQ.iplt_annual) &
      .AND.doPlantLeafOut_brch(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ))THEN
        iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ)=I
        write(101,*)'plant begin seed fill',I,NB,NZ
    ENDIF
!
!   END SEED NUMBER SET PERIOD
!
!   iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date setting for final seed number
!
  ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN
    IF(ReprodNodeNumNormByMatrgrp_brch(NB,NZ).GT.1.00_r8*GSTGR)THEN
      iPlantCalendar_brch(ipltcal_SetSeedNumber,NB,NZ)=I
      write(101,*)'plant set seed number',I,NB,NZ
    ENDIF
!
!   END SEED SIZE SET PERIOD
!
!   iPlantCalendar_brch(ipltcal_SetSeedMass,=end of setting max seed size
!
  ELSEIF(iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ).EQ.0)THEN
    IF(ReprodNodeNumNormByMatrgrp_brch(NB,NZ).GT.1.50_r8*GSTGR)THEN
      iPlantCalendar_brch(ipltcal_SetSeedMass,NB,NZ)=I
      write(101,*)'plant set seed mass',I,NB,NZ
    ENDIF
  ENDIF
  end associate
  end subroutine living_branch_phenology

end module HfuncsMod
