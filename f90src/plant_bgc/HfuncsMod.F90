module HfuncsMod
!!
! Description:
! code to do plant phenology

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use PlantAPIData
  use minimathmod, only : AZMAX1
  use StartqsMod, only : StartPlants
  implicit none

  private


  character(len=*), parameter :: mod_filename = &
  __FILE__

  real(r8) :: RLA,STK,TKCO,TFNP
  real(r8) :: WFNG
  integer :: KVSTGX,NB,NBTX,N

!
! PSILM=minimum canopy turgor potential for leaf expansion (MPa)
! PSILX=minimum canopy water potential for leafout of drought-deciduous PFT (MPa)
! PSILY=minimum canopy water potential for leafoff of drought-deciduous PFT (MPa
! GSTGG,GSTGR=normalized growth stage durations for vegetative,reproductive phenology
! NBX=maximum branch number for PFT defined by IBTYP in PFT file
! VRNE=maximum hours for leafout,leafoff
!
  real(r8), PARAMETER :: PSILM=0.1_r8
  real(r8), parameter :: PSILX=-0.2_r8
  real(r8) ,PARAMETER :: GSTGG=2.00_r8
  real(r8), PARAMETER :: GSTGR=0.667_r8
  real(r8), PARAMETER :: VRNE=3600.0_r8
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
  INTEGER :: NZ

! begin_execution
  associate(                           &
    PSICanP   =>  plt_ew%PSICanP     , &
    VRNZ    =>  plt_pheno%VRNZ   , &
    doInitPlant   =>  plt_pheno%doInitPlant  , &
    doPlantRemobilization   =>  plt_pheno%doPlantRemobilization  , &
    VSTGX   =>  plt_pheno%VSTGX  , &
    VRNY    =>  plt_pheno%VRNY   , &
    IsPlantActive   =>  plt_pheno%IsPlantActive  , &
    iPlantBranchState   =>  plt_pheno%iPlantBranchState  , &
    iPlantCalendar   =>  plt_pheno%iPlantCalendar  , &
    WSTR    =>  plt_pheno%WSTR   , &
    KVSTG   =>  plt_pheno%KVSTG  , &
    IGTYP   =>  plt_pheno%IGTYP  , &
    DayLenthCurrent    =>  plt_site%DayLenthCurrent    , &
    DATAP   =>  plt_site%DATAP   , &
    PPT     =>  plt_site%PPT     , &
    DayLenthPrev    =>  plt_site%DayLenthPrev    , &
    NP      =>  plt_site%NP      , &
    pftPlantPopulation      =>  plt_site%pftPlantPopulation      , &
    KLEAF   =>  plt_morph%KLEAF  , &
    VSTG    =>  plt_morph%VSTG   , &
    NumOfBranches_pft     =>  plt_morph%NumOfBranches_pft    , &
    NB1     =>  plt_morph%NB1      &
  )
  D9985: DO NZ=1,NP

    IF(DATAP(NZ).NE.'NO')THEN
!
!     PPT=total biome population
!
      PPT=PPT+pftPlantPopulation(NZ)
!
!         SET CROP FLAG ACCORDINGTopRootLayerTO PLANTING, HARVEST DATES, DEATH,
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
!           doInitLeafOut,doPlantLeafOut=flags for initializing leafout,leafoff
!           VRNS=leafout hours
!
        IF(iPlantCalendar(ipltcal_Emerge,NB1(NZ),NZ).NE.0.OR.doInitPlant(NZ).EQ.itrue)THEN
          D2010: DO NB=1,NumOfBranches_pft(NZ)
            IF(iPlantBranchState(NB,NZ).EQ.iliving_branch)THEN
              call living_branch_phenology(I,J,NB,nz)
            ENDIF
!
!               KVSTG=integer of most recent leaf number currently growing
!
            KVSTGX=KVSTG(NB,NZ)
            IF(VSTGX(NB,NZ).LE.ppmc)THEN
              KVSTG(NB,NZ)=INT(VSTG(NB,NZ))+1
            ELSE
              KVSTG(NB,NZ)=INT(AMIN1(VSTG(NB,NZ),VSTGX(NB,NZ)))+1
            ENDIF
            KLEAF(NB,NZ)=MIN(MaxCanopyNodes1-1,KVSTG(NB,NZ))
            IF(KVSTG(NB,NZ).GT.KVSTGX)THEN
              doPlantRemobilization(NB,NZ)=itrue
            ELSE
              doPlantRemobilization(NB,NZ)=ifalse
            ENDIF
!
!               PHENOLOGY
!
!               DayLenthPrev,DLYN=daylength of previous,current day
!               VRNY,VRNZ=hourly counter for lengthening,shortening photoperiods
!
            IF(iPlantBranchState(NB,NZ).EQ.iliving_branch.OR.doInitPlant(NZ).EQ.itrue)THEN
              IF(DayLenthCurrent.GE.DayLenthPrev)THEN
                VRNY(NB,NZ)=VRNY(NB,NZ)+1.0_r8
                VRNZ(NB,NZ)=0.0_r8
              ELSE
                VRNY(NB,NZ)=0.0_r8
                VRNZ(NB,NZ)=VRNZ(NB,NZ)+1.0_r8
              ENDIF

              call pft_specific_phenology(I,J,NZ)

            ENDIF
          ENDDO D2010
!
!             WATER STRESS INDICATOR
!
!             PSICanP=canopy total water potential
!             PSILY=minimum canopy water potential for leafoff
!             WSTR=number of hours PSICanP < PSILY (for output only)
!
          IF(PSICanP(NZ).LT.PSILY(IGTYP(NZ)))THEN
            WSTR(NZ)=WSTR(NZ)+1.0_r8
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
  ENDIF
  end associate
  end subroutine set_flags
!------------------------------------------------------------------------------------------

  subroutine root_shoot_branching(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ

! begin_execution
  associate(                            &
    CanopyNonstructElementConc_pft  =>   plt_biom%CanopyNonstructElementConc_pft  , &
    WTRVE   =>   plt_biom%WTRVE   , &
    GROUP   =>   plt_pheno%GROUP  , &
    iPlantCalendar   =>   plt_pheno%iPlantCalendar  , &
    doInitPlant   =>   plt_pheno%doInitPlant  , &
    iPlantRootState   =>   plt_pheno%iPlantRootState  , &
    ISTYP   =>   plt_pheno%ISTYP  , &
    iPlantBranchState   =>   plt_pheno%iPlantBranchState  , &
    iPlantShootState   =>   plt_pheno%iPlantShootState  , &
    IBTYP   =>   plt_pheno%IBTYP  , &
    PR      =>   plt_pheno%PR     , &
    GROUPI  =>   plt_pheno%GROUPI , &
    PB      =>   plt_pheno%PB     , &
    pftPlantPopulation      =>   plt_site%pftPlantPopulation      , &
    VRNS    =>   plt_pheno%VRNS   , &
    PSIRootTurg   =>   plt_ew%PSIRootTurg     , &
    FNOD    =>   plt_allom%FNOD   , &
    NRT     =>   plt_morph%NRT    , &
    NB1     =>   plt_morph%NB1    , &
    NumOfBranches_pft     =>   plt_morph%NumOfBranches_pft    , &
    NNOD    =>   plt_morph%NNOD   , &
    NBT     =>   plt_morph%NBT    , &
    BranchNumber_brchpft    =>   plt_morph%BranchNumber_brchpft   , &
    NGTopRootLayer     =>   plt_morph%NGTopRootLayer    , &
    XTLI    =>   plt_morph%XTLI   , &
    PSTG    =>   plt_morph%PSTG     &
  )

!
! ADD BRANCH TO SHOOT IF PLANT GROWTH STAGE, SHOOT NON-STRUCTURAL
! CONCENTRATION PERMIT
!
! doInitPlant=PFT initialization flag:0=no,1=yes
! PSIRootTurg=root turgor potential
! ISTYP=growth habit from PFT file
! iPlantCalendar(ipltcal_InitFloral,=floral initiation date
! NumOfBranches_pft=primary root axis number
! WTRVC=nonstructural C storage
! PB=nonstructural C concentration needed for branching
! iPlantBranchState=branch life flag:0=living,1=dead
! PSTG=node number
! FNOD=scales node number for perennial vegetation (e.g. trees)
! NNOD=number of concurrently growing nodes
! XTLI,GROUP=node number at planting,floral initiation
! IBTYP: setup for phenologically-driven above-ground turnover


  IF(doInitPlant(NZ).EQ.ifalse)THEN
    !plant initialized
    IF(J.EQ.1.AND.pftPlantPopulation(NZ).GT.0.0_r8)THEN
      IF(PSIRootTurg(ipltroot,NGTopRootLayer(NZ),NZ).GT.PSILM)THEN
        IF(ISTYP(NZ).NE.iplt_annual.OR.iPlantCalendar(ipltcal_InitFloral,NB1(NZ),NZ).EQ.0)THEN
          IF((NumOfBranches_pft(NZ).EQ.0.AND.WTRVE(ielmc,NZ).GT.0.0_r8) &
            .OR.(CanopyNonstructElementConc_pft(ielmc,NZ).GT.PB(NZ).AND.PB(NZ).GT.0.0_r8))THEN
            D120: DO NB=1,NumOfCanopyLayers1
              IF(iPlantBranchState(NB,NZ).EQ.iDead)THEN
                IF(NB.EQ.NB1(NZ).OR.PSTG(NB1(NZ),NZ).GT.NBT(NZ)+NNOD(NZ)/FNOD(NZ)+XTLI(NZ))THEN
                  NBT(NZ)=NBT(NZ)+1
                  NumOfBranches_pft(NZ)=MIN(NBX(IBTYP(NZ)),MAX(NB,NumOfBranches_pft(NZ)))
                  BranchNumber_brchpft(NB,NZ)=NBT(NZ)-1
                  iPlantShootState(NZ)=iLive
                  iPlantBranchState(NB,NZ)=iLive
                  VRNS(NB,NZ)=0.0_r8
                  IF(ISTYP(NZ).EQ.iplt_annual)THEN
                    GROUP(NB,NZ)=AZMAX1(GROUPI(NZ)-BranchNumber_brchpft(NB,NZ))
                  ELSE
                    GROUP(NB,NZ)=GROUPI(NZ)
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
!     NB1: number of main branch
!     CanopyNonstructElementConc_pft: canopy nonstructural element concentration
!     PSIRootTurg: root turgor pressure
!     WTRVE: non-structural carbon

      IF(PSIRootTurg(ipltroot,NGTopRootLayer(NZ),NZ).GT.PSILM)THEN
        IF(NRT(NZ).EQ.0 .OR. PSTG(NB1(NZ),NZ).GT.NRT(NZ)/FNOD(NZ)+XTLI(NZ))THEN
          IF((NRT(NZ).EQ.0 .AND. WTRVE(ielmc,NZ).GT.0.0_r8) &
            .OR.(CanopyNonstructElementConc_pft(ielmc,NZ).GT.PR(NZ) & 
            .AND.PR(NZ).GT.0.0_r8))THEN
            NRT(NZ)=MIN(NumOfCanopyLayers1,NRT(NZ)+1)
            iPlantRootState(NZ)=iLive
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

  integer :: NB,N,L,NE
  real(r8):: ARLSP
  associate(                           &
    CanopyLeafShethC_pft   =>  plt_biom%CanopyLeafShethC_pft     , &
    CanPBLeafShethC  =>  plt_biom%CanPBLeafShethC    , &
    CanPShootElmMass =>  plt_biom%CanPShootElmMass   , &
    NoduleNonstructCconc_pft =>  plt_biom%NoduleNonstructCconc_pft   , &
    CEPOLB =>  plt_biom%CEPOLB   , &
    CanopyNonstructElementConc_pft =>  plt_biom%CanopyNonstructElementConc_pft   , &
    CanopyNonstructElements_pft =>  plt_biom%CanopyNonstructElements_pft   , &
    EPOLNB =>  plt_biom%EPOLNB   , &
    EPOOLR =>  plt_biom%EPOOLR   , &
    EPOOL  =>  plt_biom%EPOOL    , &
    RootNonstructElementConcpft_vr =>  plt_biom%RootNonstructElementConcpft_vr   , &
    ZEROL  =>  plt_biom%ZEROL    , &
    ZEROP  =>  plt_biom%ZEROP    , &
    WTRTL  =>  plt_biom%WTRTL    , &
    EPOLNP =>  plt_biom%EPOLNP   , &
    WatByPCan  =>  plt_ew%WatByPCan      , &
    VHeatCapCanP  =>  plt_ew%VHeatCapCanP      , &
    NU     =>  plt_site%NU       , &
    iPlantBranchState  =>  plt_pheno%iPlantBranchState   , &
    iPlantCalendar  =>  plt_pheno%iPlantCalendar   , &
    NB1    =>  plt_morph%NB1     , &
    PrimRootDepth  =>  plt_morph%PrimRootDepth   , &
    MY     =>  plt_morph%MY      , &
    CanopyLeafA_pft  =>  plt_morph%CanopyLeafA_pft   , &
    NGTopRootLayer    =>  plt_morph%NGTopRootLayer     , &
    NIXBotRootLayer   =>  plt_morph%NIXBotRootLayer     , &
    NumOfBranches_pft    =>  plt_morph%NumOfBranches_pft     , &
    BranchNumber_brchpft   =>  plt_morph%BranchNumber_brchpft    , &
    HypoctoylHeight  =>  plt_morph%HypoctoylHeight   , &
    SeedinDepth  =>  plt_morph%SeedinDepth   , &
    CanopyStemA_pft  =>  plt_morph%CanopyStemA_pft   , &
    NI     =>  plt_morph%NI        &
  )
  plt_bgcr%RootGasLoss_disturb(idg_beg:idg_end-1,NZ)=0.0_r8
  CanopyNonstructElements_pft(1:NumOfPlantChemElements,NZ)=0.0_r8
  NI(NZ)=NIXBotRootLayer(NZ)
  NGTopRootLayer(NZ)=MIN(NI(NZ),MAX(NGTopRootLayer(NZ),NU))
  NB1(NZ)=1
  NBTX=1.0E+06_r8
!
! TOTAL PLANT NON-STRUCTURAL C, N, P
!
! CPOOL*,ZPOOL*,PPOOL*=non-structl C,N,P in branch(NB),canopy(g)
! CPOLN*,ZPOLN*,PPOLN*=non-structl C,N,P in branch,canopy nodules (g)
! NB1=main branch number
!
  DO NE=1,NumOfPlantChemElements
    D140: DO NB=1,NumOfBranches_pft(NZ)
      IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
        CanopyNonstructElements_pft(NE,NZ)=CanopyNonstructElements_pft(NE,NZ)+EPOOL(NE,NB,NZ)
        EPOLNP(NE,NZ)=EPOLNP(NE,NZ)+EPOLNB(NE,NB,NZ)
      ENDIF
    ENDDO D140
  ENDDO

  DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState(NB,NZ).EQ.iLive)THEN
      IF(BranchNumber_brchpft(NB,NZ).LT.NBTX)THEN
        NB1(NZ)=NB
        NBTX=BranchNumber_brchpft(NB,NZ)
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
      IF(WTRTL(N,L,NZ).GT.ZEROL(NZ))THEN
        DO NE=1,NumOfPlantChemElements
          RootNonstructElementConcpft_vr(NE,N,L,NZ)=AZMAX1(EPOOLR(NE,N,L,NZ)/WTRTL(N,L,NZ))
        ENDDO
      ELSE
        DO NE=1,NumOfPlantChemElements
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
    DO NE=1,NumOfPlantChemElements
      CanopyNonstructElementConc_pft(NE,NZ)=AZMAX1(AMIN1(1.0_r8,CanopyNonstructElements_pft(NE,NZ)/CanopyLeafShethC_pft(NZ)))
    ENDDO
    NoduleNonstructCconc_pft(NZ)=AZMAX1(AMIN1(1.0_r8,EPOLNP(ielmc,NZ)/CanopyLeafShethC_pft(NZ)))
  ELSE
    CanopyNonstructElementConc_pft(1:NumOfPlantChemElements,NZ)=1.0_r8
    NoduleNonstructCconc_pft(NZ)=1.0_r8
  ENDIF
  DO NE=1,NumOfPlantChemElements
    D190: DO NB=1,NumOfBranches_pft(NZ)
      IF(CanPBLeafShethC(NB,NZ).GT.ZEROP(NZ))THEN
        CEPOLB(NE,NB,NZ)=AZMAX1(EPOOL(NE,NB,NZ)/CanPBLeafShethC(NB,NZ))
      ELSE
        CEPOLB(NE,NB,NZ)=1.0_r8
      ENDIF
    ENDDO D190
  ENDDO
!
! EMERGENCE DATE FROM COTYLEDON HEIGHT, LEAF AREA, ROOT DEPTH
!
! iPlantCalendar(ipltcal_Emerge,=emergence date
! CanopyLeafA_pft,CanopyStemA_pft=leaf,stalk areas
! HypoctoylHeight=hypocotyledon height
! SeedinDepth=seeding depth
! PrimRootDepth=primary root depth
! VHeatCapCanP,WTSHT,WatByPCan=canopy heat capacity,mass,water content
!
  IF(iPlantCalendar(ipltcal_Emerge,NB1(NZ),NZ).EQ.0)THEN
    ARLSP=CanopyLeafA_pft(NZ)+CanopyStemA_pft(NZ)
    IF((HypoctoylHeight(NZ).GT.SeedinDepth(NZ)).AND.(ARLSP.GT.ZEROL(NZ)) &
      .AND.(PrimRootDepth(1,1,NZ).GT.SeedinDepth(NZ)+ppmc))THEN
      iPlantCalendar(ipltcal_Emerge,NB1(NZ),NZ)=I
      VHeatCapCanP(NZ)=cpw*(CanPShootElmMass(ielmc,NZ)*10.0E-06_r8+WatByPCan(NZ))
    ENDIF
  ENDIF
  end associate
  end subroutine stage_phenology_vars
!------------------------------------------------------------------------------------------

  subroutine pft_specific_phenology(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ

  associate(                           &
    PSICanP   =>  plt_ew%PSICanP     , &
    PSICanPTurg   =>  plt_ew%PSICanPTurg     , &
    DayLenthCurrent    =>  plt_site%DayLenthCurrent    , &
    DayLenthPrev    =>  plt_site%DayLenthPrev    , &
    DayLenthMax    =>  plt_site%DayLenthMax    , &
    ALAT    =>  plt_site%ALAT    , &
    doPlantLeafOut   =>  plt_pheno%doPlantLeafOut  , &
    iPlantCalendar   =>  plt_pheno%iPlantCalendar  , &
    IGTYP   =>  plt_pheno%IGTYP  , &
    VRNZ    =>  plt_pheno%VRNZ   , &
    VRNS    =>  plt_pheno%VRNS   , &
    HourThreshold4LeafOut   =>  plt_pheno%HourThreshold4LeafOut  , &
    IWTYP   =>  plt_pheno%IWTYP  , &
    TCZ     =>  plt_pheno%TCZ    , &
    VRNY    =>  plt_pheno%VRNY   , &
    TCG     =>  plt_pheno%TCG    , &
    CTC     =>  plt_pheno%CTC    , &
    TCX     =>  plt_pheno%TCX    , &
    Hours4LeafOff    =>  plt_pheno%Hours4LeafOff   , &
    HourThreshold4LeafOff   =>  plt_pheno%HourThreshold4LeafOff  , &
    doPlantLeaveOff   =>  plt_pheno%doPlantLeaveOff    &
  )
!
! CALCULATE EVERGREEN PHENOLOGY DURINGTopRootLayerLENGTHENINGTopRootLayerPHOTOPERIODS
!
! IWTYP=phenology type from PFT file
! DayLenthPrev,DLYN=daylength of previous,current day
! VRNS,Hours4LeafOff=leafout,leafoff hours
! VRNY=hourly counter for lengthening photoperiods
! doPlantLeaveOff=flag for enabling leafoff:0=enable,1=disable
! ALAT=latitude
!
  IF(IWTYP(NZ).EQ.0)THEN
    IF(DayLenthCurrent.GE.DayLenthPrev)THEN
      VRNS(NB,NZ)=VRNY(NB,NZ)
      IF(VRNS(NB,NZ).GE.HourThreshold4LeafOut(NB,NZ) &
        .OR.(ALAT.GT.0.0_r8.AND.I.EQ.173) &
        .OR.(ALAT.LT.0.0_r8.AND.I.EQ.355))THEN
        Hours4LeafOff(NB,NZ)=0.0
        doPlantLeaveOff(NB,NZ)=iEnable
      ENDIF
    ENDIF
!
!   CALCULATE EVERGREEN PHENOLOGY DURINGTopRootLayerSHORTENINGTopRootLayerPHOTOPERIODS
!
!   VRNS,Hours4LeafOff=leafout,leafoff hours
!   VRNZ=hourly counter for shortening photoperiods
!   doPlantLeafOut=flag for enabling leafout:0=enable,1=disable
!   ALAT=latitude
!
    IF(DayLenthCurrent.LT.DayLenthPrev)THEN
      Hours4LeafOff(NB,NZ)=VRNZ(NB,NZ)
      IF(Hours4LeafOff(NB,NZ).GE.HourThreshold4LeafOff(NB,NZ) &
        .OR.(ALAT.GT.0.0_r8.AND.I.EQ.355) &
        .OR.(ALAT.LT.0.0_r8.AND.I.EQ.173))THEN
        VRNS(NB,NZ)=0.0
        doPlantLeafOut(NB,NZ)=iEnable
      ENDIF
    ENDIF
!
!   CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayerHOURS ABOVE
!   SPECIFIED TEMPERATURE DURINGTopRootLayerLENGTHENINGTopRootLayerPHOTOPERIODS
!
!   IWTYP=phenology type from PFT file
!   DayLenthPrev,DLYN=daylength of previous,current day
!   VRNS,VRNL=leafout hours,hours required for leafout
!   Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!   doPlantLeafOut,doPlantLeaveOff=flag for enabling leafout,leafoff:0=enable,1=disable
!   TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!   ALAT=latitude
!   iPlantCalendar(ipltcal_InitFloral,=date of floral initiation
!
  ELSEIF(IWTYP(NZ).EQ.1)THEN
    IF((DayLenthCurrent.GE.DayLenthPrev.OR.(DayLenthCurrent.LT.DayLenthPrev.AND.Hours4LeafOff(NB,NZ).LT.HourThreshold4LeafOff(NB,NZ))) &
      .AND.doPlantLeafOut(NB,NZ).EQ.iEnable)THEN
      IF(TCG(NZ).GE.TCZ(NZ))THEN
        VRNS(NB,NZ)=VRNS(NB,NZ)+1.0
      ENDIF
      IF(VRNS(NB,NZ).LT.HourThreshold4LeafOut(NB,NZ))THEN
        IF(TCG(NZ).LT.CTC(NZ))THEN
          VRNS(NB,NZ)=AZMAX1(VRNS(NB,NZ)-1.0)
        ENDIF
      ENDIF
      IF(VRNS(NB,NZ).GE.HourThreshold4LeafOut(NB,NZ) &
        .OR.(ALAT.GT.0.0.AND.I.EQ.173) &
        .OR.(ALAT.LT.0.0.AND.I.EQ.355))THEN
        Hours4LeafOff(NB,NZ)=0.0
      ENDIF
    ENDIF

    IF(iPlantCalendar(ipltcal_InitFloral,NB,NZ).NE.0.OR. &
      (DayLenthCurrent.LT.DayLenthPrev.AND.DayLenthCurrent.LT.12.0_r8))THEN
      doPlantLeaveOff(NB,NZ)=iEnable
    ENDIF
!
!     CALCULATE WINTER DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayerHOURS BELOW
!     SPECIFIED TEMPERATURE DURINGTopRootLayerSHORTENINGTopRootLayerPHOTOPERIODS
!
!     DayLenthPrev,DLYN=daylength of previous,current day
!     VRNS,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!     doPlantLeafOut,doPlantLeaveOff=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!     ALAT=latitude
!     iPlantCalendar(ipltcal_InitFloral,=date of floral initiation
!
    IF(DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeaveOff(NB,NZ).EQ.iEnable &
      .AND.iPlantCalendar(ipltcal_InitFloral,NB,NZ).NE.0)THEN
      IF(TCG(NZ).LE.TCX(NZ))THEN
        Hours4LeafOff(NB,NZ)=Hours4LeafOff(NB,NZ)+1.0_r8
      ENDIF
      IF(Hours4LeafOff(NB,NZ).GE.HourThreshold4LeafOff(NB,NZ).AND.doPlantLeafOut(NB,NZ).EQ.iDisable)THEN
        VRNS(NB,NZ)=0.0
        doPlantLeafOut(NB,NZ)=iEnable
      ENDIF
    ENDIF

!
!     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayerHOURS
!     ABOVE SPECIFIED WATER POTENTIAL DURINGTopRootLayerDORMANCY
!
!     IWTYP=phenology type from PFT file
!     VRNS,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff=leafoff hours
!     doPlantLeafOut,doPlantLeaveOff=flag for enabling leafout,leafoff:0=enable,1=disable
!     PSICanP=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     iPlantCalendar(ipltcal_InitFloral,=date of floral initiation
!
  ELSEIF(IWTYP(NZ).EQ.2.OR.IWTYP(NZ).EQ.4.OR.IWTYP(NZ).EQ.5)THEN
    IF(doPlantLeafOut(NB,NZ).EQ.iEnable)THEN
      IF(PSICanP(NZ).GE.PSILX)THEN
        VRNS(NB,NZ)=VRNS(NB,NZ)+1.0_r8
      ENDIF
      IF(VRNS(NB,NZ).LT.HourThreshold4LeafOut(NB,NZ))THEN
        IF(PSICanP(NZ).LT.PSILY(IGTYP(NZ)))THEN
          VRNS(NB,NZ)=AZMAX1(VRNS(NB,NZ)-12.0)
        ENDIF
      ENDIF
      IF(VRNS(NB,NZ).GE.HourThreshold4LeafOut(NB,NZ))THEN
        Hours4LeafOff(NB,NZ)=0.0_r8
        IF(iPlantCalendar(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff(NB,NZ)=iEnable
      ENDIF
    ENDIF
    IF(iPlantCalendar(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff(NB,NZ)=iEnable
!
!     CALCULATE DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATINGTopRootLayerHOURS
!     BELOW SPECIFIED WATER POTENTIAL DURINGTopRootLayerGROWINGTopRootLayerSEASON
!
!     VRNS=leafout hours,hours required for leafout
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!     doPlantLeafOut,doPlantLeaveOff=flag for enabling leafout,leafoff:0=enable,1=disable
!     PSICanP=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     iPlantCalendar(ipltcal_InitFloral,=date of floral initiation
!     VRNY,VRNZ=hourly counter for lengthening,shortening photoperiods
!     VRNE=maximum hours for leafout,leafoff
!
      IF(doPlantLeafOut(NB,NZ).EQ.iDisable.AND.doPlantLeaveOff(NB,NZ).EQ.iEnable)THEN
        IF(PSICanP(NZ).LT.PSILY(IGTYP(NZ)))THEN
          Hours4LeafOff(NB,NZ)=Hours4LeafOff(NB,NZ)+1.0_r8
        ENDIF
        IF(IWTYP(NZ).EQ.4)THEN
          IF(VRNZ(NB,NZ).GT.VRNE)THEN
            Hours4LeafOff(NB,NZ)=VRNZ(NB,NZ)
          ENDIF
        ELSEIF(IWTYP(NZ).EQ.5)THEN
          IF(VRNY(NB,NZ).GT.VRNE)THEN
            Hours4LeafOff(NB,NZ)=VRNY(NB,NZ)
          ENDIF
        ENDIF
        IF(Hours4LeafOff(NB,NZ).GE.HourThreshold4LeafOff(NB,NZ).AND.doPlantLeafOut(NB,NZ).EQ.iDisable)THEN
          VRNS(NB,NZ)=0.0_r8
          doPlantLeafOut(NB,NZ)=iEnable
        ENDIF
      ENDIF
!
!     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
!     HOURS ABOVE SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
!     LENGTHENINGTopRootLayerPHOTOPERIODS
!
!     IWTYP=phenology type from PFT file
!     DayLenthPrev,DLYN=daylength of previous,current day
!     VRNS,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!     PSICanP=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     doPlantLeafOut,doPlantLeaveOff=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!     ALAT=latitude
!     iPlantCalendar(ipltcal_InitFloral,=date of floral initiation
!
    ELSEIF(IWTYP(NZ).EQ.3)THEN
      IF((DayLenthCurrent.GE.DayLenthPrev.OR.DayLenthCurrent.GE.DayLenthMax-2.0_r8).AND.doPlantLeafOut(NB,NZ).EQ.iEnable)THEN
        IF(TCG(NZ).GE.TCZ(NZ).AND.PSICanPTurg(NZ).GT.PSILM)THEN
          VRNS(NB,NZ)=VRNS(NB,NZ)+1.0_r8
        ENDIF
        IF(VRNS(NB,NZ).LT.HourThreshold4LeafOut(NB,NZ))THEN
          IF(TCG(NZ).LT.CTC(NZ).OR.PSICanPTurg(NZ).LT.PSILM)THEN
            VRNS(NB,NZ)=AZMAX1(VRNS(NB,NZ)-1.5_r8)
          ENDIF
        ENDIF
        IF(VRNS(NB,NZ).GE.HourThreshold4LeafOut(NB,NZ))THEN
          Hours4LeafOff(NB,NZ)=0.0_r8
          IF(iPlantCalendar(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff(NB,NZ)=iEnable
        ENDIF
      ENDIF
      IF(iPlantCalendar(ipltcal_InitFloral,NB,NZ).NE.0)doPlantLeaveOff(NB,NZ)=iEnable
!
!     CALCULATE WINTER AND DROUGHT DECIDUOUS PHENOLOGY BY ACCUMULATING
!     HOURS BELOW SPECIFIED TEMPERATURE OR WATER POTENTIAL DURING
!     SHORTENINGTopRootLayerPHOTOPERIODS
!
!     DayLenthPrev,DLYN=daylength of previous,current day
!     VRNS,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!     doPlantLeafOut,doPlantLeaveOff=flag for enabling leafout,leafoff:0=enable,1=disable
!     TCG,TCZ,CTC=canopy temp,leafout threshold temp,chilling temp
!     PSICanP=canopy total water potential
!     PSILX,PSILY=minimum canopy water potential for leafout,leafoff
!     ALAT=latitude
!     iPlantCalendar(ipltcal_InitFloral,=date of floral initiation
!
      IF((DayLenthCurrent.LT.DayLenthPrev.OR.DayLenthCurrent.LT.24.0_r8-DayLenthMax+2.0_r8) &
        .AND.doPlantLeaveOff(NB,NZ).EQ.iEnable)THEN
        IF(TCG(NZ).LE.TCX(NZ).OR.PSICanP(NZ).LT.PSILY(IGTYP(NZ)))THEN
          Hours4LeafOff(NB,NZ)=Hours4LeafOff(NB,NZ)+1.0_r8
        ENDIF
        IF(Hours4LeafOff(NB,NZ).GE.HourThreshold4LeafOff(NB,NZ).AND.doPlantLeafOut(NB,NZ).EQ.iDisable)THEN
          VRNS(NB,NZ)=0.0
          doPlantLeafOut(NB,NZ)=iEnable
        ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine pft_specific_phenology
!------------------------------------------------------------------------------------------

  subroutine living_branch_phenology(I,J,NB,NZ)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NB  !plant branch id
  integer, intent(in) :: NZ  !plant species id

  real(r8) :: ACTV,OFNG
  real(r8) :: PPD
  real(r8) :: RTK
  real(r8) :: RNI

! begin_execution
  associate(                        &
    iYearPlanting    =>  plt_distb%iYearPlanting    , &
    iDayPlanting   =>  plt_distb%iDayPlanting   , &
    SnowDepth   =>  plt_ew%SnowDepth      , &
    PSICanP   =>  plt_ew%PSICanP      , &
    ZERO    =>  plt_site%ZERO     , &
    DayLenthCurrent    =>  plt_site%DayLenthCurrent     , &
    iYearCurrent    =>  plt_site%iYearCurrent     , &
    DayLenthPrev    =>  plt_site%DayLenthPrev     , &
    TKG     =>  plt_pheno%TKG     , &
    iPlantCalendar   =>  plt_pheno%iPlantCalendar   , &
    XRLA    =>  plt_pheno%XRLA    , &
    GSTGF   =>  plt_pheno%GSTGF   , &
    OFFST   =>  plt_pheno%OFFST   , &
    ISTYP   =>  plt_pheno%ISTYP   , &
    doPlantLeafOut   =>  plt_pheno%doPlantLeafOut   , &
    VSTG    =>  plt_morph%VSTG    , &
    GSTGI   =>  plt_pheno%GSTGI   , &
    DGSTGF  =>  plt_pheno%DGSTGF  , &
    doPlantSenescence   =>  plt_pheno%doPlantSenescence   , &
    VSTGX   =>  plt_pheno%VSTGX   , &
    IDTYP   =>  plt_pheno%IDTYP   , &
    TGSTGF  =>  plt_pheno%TGSTGF  , &
    XDL     =>  plt_pheno%XDL     , &
    DGSTGI  =>  plt_pheno%DGSTGI  , &
    XPPD    =>  plt_pheno%XPPD    , &
    TGSTGI  =>  plt_pheno%TGSTGI  , &
    GROUP   =>  plt_pheno%GROUP   , &
    IWTYP   =>  plt_pheno%IWTYP   , &
    HourThreshold4LeafOut   =>  plt_pheno%HourThreshold4LeafOut   , &
    IPTYP   =>  plt_pheno%IPTYP   , &
    doInitLeafOut   =>  plt_pheno%doInitLeafOut   , &
    GROUPI  =>  plt_pheno%GROUPI  , &
    PlantO2Stress    =>  plt_pheno%PlantO2Stress    , &
    VRNS    =>  plt_pheno%VRNS    , &
    Hours4LeafOff    =>  plt_pheno%Hours4LeafOff    , &
    HourThreshold4LeafOff   =>  plt_pheno%HourThreshold4LeafOff   , &
    XRNI    =>  plt_pheno%XRNI    , &
    CanopyHeight      =>  plt_morph%CanopyHeight      , &
    PSTG    =>   plt_morph%PSTG   , &
    PSTGI   =>   plt_morph%PSTGI  , &
    PSTGF   =>   plt_morph%PSTGF  , &
    NB1     =>   plt_morph%NB1      &
  )
  IF(iPlantCalendar(ipltcal_Emerge,NB,NZ).EQ.0)THEN
    !plant emergence
    iPlantCalendar(ipltcal_Emerge,NB,NZ)=I
    doInitLeafOut(NB,NZ)=iDisable
    doPlantLeafOut(NB,NZ)=iEnable
    VRNS(NB,NZ)=0.5_r8*VRNS(NB1(NZ),NZ)
  ENDIF
!
! CALCULATE NODE INITIATION AND LEAF APPEARANCE RATES
! FROM TEMPERATURE FUNCTION CALCULATED IN 'UPTAKE' AND
! RATES AT 25C ENTERED IN 'READQ' EXCEPT WHEN DORMANT
!
! IWTYP=phenology type from PFT file, 0=evergreen,1=cold deciduous, 2=drought deciduous,3=1+2
! Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
! TKG,TKCO=canopy temperature,canopy temp used in Arrhenius eqn
! OFFST=shift in Arrhenius curve for thermal adaptation
! TFNP=temperature function for phenology (25 oC =1 )
! 8.313,710.0=gas constant,enthalpy
! 60000,197500,218500=energy of activn,high,low temp inactivn(KJ mol-1)
! RNI,RLA=rates of node initiation,leaf appearance
! XRNI,XRLA=rate of node initiation,leaf appearance at 25 oC (h-1)
!
  IF(IWTYP(NZ).EQ.0.OR.Hours4LeafOff(NB,NZ).LT.HourThreshold4LeafOff(NB,NZ))THEN
    TKCO=TKG(NZ)+OFFST(NZ)
    RTK=RGAS*TKCO
    STK=710.0_r8*TKCO
    ACTV=1+EXP((197500_r8-STK)/RTK)+EXP((STK-218500._r8)/RTK)
    TFNP=EXP(24.269_r8-60000._r8/RTK)/ACTV
    
    RNI=AZMAX1(XRNI(NZ)*TFNP)
    RLA=AZMAX1(XRLA(NZ)*TFNP)
!
!   NODE INITIATION AND LEAF APPEARANCE RATES SLOWED BY LOW TURGOR
!
!   PSICanPTurg=leaf turgor potential
!   WFNG=water stress effect on phenology
!
    IF(ISTYP(NZ).EQ.iplt_annual)THEN
      IF(iPlantCalendar(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
        WFNG=EXP(0.025_r8*PSICanP(NZ))
        RNI=RNI*WFNG
        RLA=RLA*WFNG
      ENDIF
      IF(iPlantCalendar(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
        OFNG=SQRT(PlantO2Stress(NZ))
        RNI=RNI*OFNG
        RLA=RLA*OFNG
      ENDIF
    ENDIF
!
!   ACCUMULATE NODE INITIATION AND LEAF APPEARANCE RATES
!   INTO TOTAL NUMBER OF NODES AND LEAVES
!
!   PSTG,VSTG=number of nodes initiated,leaves appeared
!
    PSTG(NB,NZ)=PSTG(NB,NZ)+RNI
    VSTG(NB,NZ)=VSTG(NB,NZ)+RLA
!
!   USE TOTAL NUMBER OF NODES TO CALCULATE PROGRESSION THROUGH
!   VEGETATIVE AND REPRODUCTIVE GROWTH STAGES. THIS PROGRESSION
!   IS USED TO SET START AND END DATES FOR GROWTH STAGES BELOW
!
!   GSTGI=vegetative node number normalized for maturity group
!   PSTGI=node number at floral initiation
!   GROUPI=node number required for floral initiation
!   DGSTGI,TGSTGI=hourly,total change in GSTGI
!   GSTGF=reproductive node number normalized for maturity group
!   PSTGF=node number at flowering
!   DGSTGF,TGSTGF=hourly,total change in GSTGF
!   doPlantSenescence=PFT senescence flag
!
    IF(iPlantCalendar(ipltcal_InitFloral,NB,NZ).NE.0)THEN
      GSTGI(NB,NZ)=(PSTG(NB,NZ)-PSTGI(NB,NZ))/GROUPI(NZ)
      DGSTGI(NB,NZ)=RNI/(GROUPI(NZ)*GSTGG)
      TGSTGI(NB,NZ)=TGSTGI(NB,NZ)+DGSTGI(NB,NZ)
    ENDIF
    IF(iPlantCalendar(ipltcal_Anthesis,NB,NZ).NE.0)THEN
      GSTGF(NB,NZ)=(PSTG(NB,NZ)-PSTGF(NB,NZ))/GROUPI(NZ)
      DGSTGF(NB,NZ)=RNI/(GROUPI(NZ)*GSTGR)
      TGSTGF(NB,NZ)=TGSTGF(NB,NZ)+DGSTGF(NB,NZ)
    ENDIF
    doPlantSenescence(NB,NZ)=itrue
  ELSE
    doPlantSenescence(NB,NZ)=ifalse
  ENDIF
!
! REPRODUCTIVE GROWTH STAGES ADVANCE WHEN THRESHOLD NUMBER
! OF NODES HAVE BEEN INITIATED. FIRST DETERMINE PHOTOPERIOD
! AND TEMPERATURE EFFECTS ON FINAL VEG NODE NUMBER FROM
! NUMBER OF INITIATED NODES
!
! PSTG=number of nodes initiated
! PSTGI=node number at floral initiation
! GROUP=node number required for floral initiation
! VRNS,VRNL=leafout hours,hours required for leafout
! DayLenthPrev,DLYN=daylength of previous,current day
! ISTYP=growth habit from PFT file
! IWTYP=phenology type from PFT file
! ZC,SnowDepth=canopy height,snowpack depth
!
  IF(iPlantCalendar(ipltcal_InitFloral,NB,NZ).EQ.0)THEN
    IF(PSTG(NB,NZ).GT.GROUP(NB,NZ)+PSTGI(NB,NZ) &
      .AND.((VRNS(NB,NZ).GE.HourThreshold4LeafOut(NB,NZ)) &
      .OR.(I.GE.iDayPlanting(NZ).AND.iYearCurrent.EQ.iYearPlanting(NZ) &
      .AND.DayLenthCurrent.GT.DayLenthPrev)) &
      .OR.(((ISTYP(NZ).EQ.iplt_preanu.AND.(IWTYP(NZ).EQ.1 &
      .OR.IWTYP(NZ).EQ.3)) &
      .OR.(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).EQ.0)) &
      .AND.CanopyHeight(NZ).GE.SnowDepth-ZERO &
      .AND.DayLenthCurrent.LT.DayLenthPrev))THEN
!
!     FINAL VEGETATIVE NODE NUMBER DEPENDS ON PHOTOPERIOD FROM 'DAY'
!     AND ON MATURITY GROUP, CRITICAL PHOTOPERIOD AND PHOTOPERIOD
!     SENSITIVITY ENTERED IN 'READQ'
!
!     IPTYP=photoperiod type from PFT file
!     PPD=photoperiod sensitivity
!     XDL=critical photoperiod from PFT file
!     iPlantCalendar(ipltcal_InitFloral,=date of floral initiation
!     VSTGX=node number on date of floral initiation
!
      IF(IPTYP(NZ).EQ.0)THEN
        PPD=0.0_r8
      ELSE
        PPD=AZMAX1(XDL(NZ)-DayLenthCurrent)
        IF(IPTYP(NZ).EQ.1.AND.DayLenthCurrent.GE.DayLenthPrev)PPD=0.0_r8
      ENDIF

      IF(IPTYP(NZ).EQ.0 &
        .OR.(IPTYP(NZ).EQ.1.AND.PPD.GT.XPPD(NZ)) &
        .OR.(IPTYP(NZ).EQ.2.AND.PPD.LT.XPPD(NZ)) &
        .OR.(((ISTYP(NZ).EQ.iplt_preanu.AND.(IWTYP(NZ).EQ.1 &
        .OR.IWTYP(NZ).EQ.3)) &
        .OR.(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).EQ.0)) &
        .AND.CanopyHeight(NZ).GE.SnowDepth-ZERO &
        .AND.DayLenthCurrent.LT.DayLenthPrev))THEN
        iPlantCalendar(ipltcal_InitFloral,NB,NZ)=I
        PSTGI(NB,NZ)=PSTG(NB,NZ)
        IF(ISTYP(NZ).EQ.iplt_annual.AND.IDTYP(NZ).EQ.0)THEN
          VSTGX(NB,NZ)=PSTG(NB,NZ)
        ENDIF
      ENDIF
    ENDIF
!
!   STEM ELONGATION
!
!   GSTGI=vegetative node number normalized for maturity group
!   GSTGG=normalized growth stage durations for vegetative phenology
!   iPlantCalendar(ipltcal_Jointing,=start of stem elongation and setting max seed number
!
  ELSEIF(iPlantCalendar(ipltcal_Jointing,NB,NZ).EQ.0)THEN
    IF(GSTGI(NB,NZ).GT.0.25_r8*GSTGG &
      .OR.((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.ISTYP(NZ).NE.iplt_annual.AND.IPTYP(NZ).NE.1 &
      .AND.DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeafOut(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ)) &
      .OR.(IWTYP(NZ).EQ.2.AND.ISTYP(NZ).EQ.iplt_annual) &
      .AND.doPlantLeafOut(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ))THEN
      iPlantCalendar(ipltcal_Jointing,NB,NZ)=I
    ENDIF
!
!   iPlantCalendar(ipltcal_Elongation,=mid stem elongation
!
  ELSEIF(iPlantCalendar(ipltcal_Elongation,NB,NZ).EQ.0)THEN
    IF(GSTGI(NB,NZ).GT.0.50*GSTGG &
      .OR.((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.ISTYP(NZ).NE.iplt_annual.AND.IPTYP(NZ).NE.1 &
      .AND.DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeafOut(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ)) &
      .OR.(IWTYP(NZ).EQ.2.AND.ISTYP(NZ).EQ.iplt_annual) &
      .AND.doPlantLeafOut(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ))THEN
      iPlantCalendar(ipltcal_Elongation,NB,NZ)=I
      IF(ISTYP(NZ).EQ.iplt_annual.AND.IDTYP(NZ).NE.0)THEN
        VSTGX(NB,NZ)=PSTG(NB,NZ)
      ENDIF
    ENDIF
!
!   iPlantCalendar(ipltcal_Heading,=end of stem elongation and setting max seed number
!
  ELSEIF(iPlantCalendar(ipltcal_Heading,NB,NZ).EQ.0)THEN
    IF(GSTGI(NB,NZ).GT.1.00*GSTGG &
      .OR.((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.ISTYP(NZ).NE.iplt_annual.AND.IPTYP(NZ).NE.1 &
      .AND.DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeafOut(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ)) &
      .OR.(IWTYP(NZ).EQ.2.AND.ISTYP(NZ).EQ.iplt_annual) &
      .AND.doPlantLeafOut(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ))THEN
      iPlantCalendar(ipltcal_Heading,NB,NZ)=I
    ENDIF
!
!   ANTHESIS OCCURS WHEN THE NUMBER OF LEAVES THAT HAVE APPEARED
!   EQUALS THE NUMBER OF NODES INITIATED WHEN THE FINAL VEGETATIVE
!   NODE NUMBER WAS SET ABOVE
!
!   iPlantCalendar(ipltcal_Anthesis,=start of anthesis and setting final seed number
!   VSTG=number of leaves appeared
!   PSTGI=node number at floral initiation
!   ISTYP,IWTYP,IPTYP=growth habit,phenology,photoperiod type from PFT file
!   doPlantLeafOut=flag for enabling leafout:0=enable,1=disable
!   Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!   DayLenthPrev,DLYN=daylength of previous,current day
!   PSTGF=number of nodes at anthesis
!
  ELSEIF(iPlantCalendar(ipltcal_Anthesis,NB,NZ).EQ.0)THEN
    IF((VSTG(NB,NZ).GT.PSTGI(NB,NZ)) &
      .OR.(ISTYP(NZ).NE.iplt_annual.AND.iPlantCalendar(ipltcal_Heading,NB,NZ).NE.0) &
      .OR.((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.ISTYP(NZ).NE.iplt_annual.AND.IPTYP(NZ).NE.1 &
      .AND.DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeafOut(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ)) &
      .OR.(IWTYP(NZ).EQ.2.AND.ISTYP(NZ).EQ.iplt_annual) &
      .AND.doPlantLeafOut(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ))THEN
      IF(NB.EQ.NB1(NZ).OR.iPlantCalendar(ipltcal_Anthesis,NB1(NZ),NZ).NE.0)THEN
        iPlantCalendar(ipltcal_Anthesis,NB,NZ)=I
        PSTGF(NB,NZ)=PSTG(NB,NZ)
      ENDIF
    ENDIF
!
!   START GRAIN FILL PERIOD
!
!   iPlantCalendar(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
!   GSTGF=reproductive node number normalized for maturity group
!   GSTGR=normalized growth stage durations for reproductive phenology
!
!
  ELSEIF(iPlantCalendar(ipltcal_BeginSeedFill,NB,NZ).EQ.0)THEN
    IF(GSTGF(NB,NZ).GT.0.50*GSTGR &
      .OR.((IWTYP(NZ).EQ.1.OR.IWTYP(NZ).EQ.3) &
      .AND.ISTYP(NZ).NE.iplt_annual.AND.IPTYP(NZ).NE.1 &
      .AND.DayLenthCurrent.LT.DayLenthPrev.AND.doPlantLeafOut(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ)) &
      .OR.(IWTYP(NZ).EQ.2.AND.ISTYP(NZ).EQ.iplt_annual) &
      .AND.doPlantLeafOut(NB,NZ).EQ.iDisable &
      .AND.Hours4LeafOff(NB,NZ).GT.HourThreshold4LeafOff(NB,NZ))THEN
        iPlantCalendar(ipltcal_BeginSeedFill,NB,NZ)=I
    ENDIF
!
!   END SEED NUMBER SET PERIOD
!
!   iPlantCalendar(ipltcal_SetSeedNumber,=end date setting for final seed number
!
  ELSEIF(iPlantCalendar(ipltcal_SetSeedNumber,NB,NZ).EQ.0)THEN
    IF(GSTGF(NB,NZ).GT.1.00*GSTGR)THEN
      iPlantCalendar(ipltcal_SetSeedNumber,NB,NZ)=I
    ENDIF
!
!   END SEED SIZE SET PERIOD
!
!   iPlantCalendar(ipltcal_SetSeedMass,=end of setting max seed size
!
  ELSEIF(iPlantCalendar(ipltcal_SetSeedMass,NB,NZ).EQ.0)THEN
    IF(GSTGF(NB,NZ).GT.1.50*GSTGR)THEN
      iPlantCalendar(ipltcal_SetSeedMass,NB,NZ)=I
    ENDIF
  ENDIF
  end associate
  end subroutine living_branch_phenology

end module HfuncsMod
