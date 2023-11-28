module LitterFallMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: ResetDeadBranch
  contains
!------------------------------------------------------------------------------------------

  subroutine ResetDeadBranch(I,J,NZ,CPOOLK)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8),intent(inout) :: CPOOLK(NumOfCanopyLayers1,JP1)
  integer :: IDTHY

!     begin_execution
  associate(                                 &
    iYearPlantHarvest       =>   plt_distb%iYearPlantHarvest     , &
    JHVST      =>   plt_distb%JHVST    , &
    iDayPlantHarvest      =>   plt_distb%iDayPlantHarvest    , &
    UVOLO      =>   plt_ew%UVOLO       , &
    CanWatP      =>   plt_ew%CanWatP       , &
    RootChemElmnts_pft      =>   plt_biom%RootChemElmnts_pft     , &
    NonstructalChemElmnts_pft      =>   plt_biom%NonstructalChemElmnts_pft     , &
    iPlantPhenologyPattern     =>   plt_pheno%iPlantPhenologyPattern   , &
    iPlantRootState     =>   plt_pheno%iPlantRootState   , &
    iPlantShootState      =>   plt_pheno%iPlantShootState    , &
    doInitPlant      =>   plt_pheno%doInitPlant    , &
    HoursCanopyPSITooLow       =>   plt_pheno%HoursCanopyPSITooLow     , &
    iPlantCalendar      =>   plt_pheno%iPlantCalendar    , &
    pftPlantPopulation         =>   plt_site%pftPlantPopulation        , &
    iYearCurrent       =>   plt_site%iYearCurrent      , &
    VOLWOU     =>   plt_site%VOLWOU    , &
    ZNOON      =>   plt_site%ZNOON     , &
    HypoctoylHeight      =>   plt_morph%HypoctoylHeight    , &
    NumOfBranches_pft        =>   plt_morph%NumOfBranches_pft      , &
    NumOfMainBranch_pft        =>   plt_morph%NumOfMainBranch_pft      , &
    BranchNumber_pft        =>   plt_morph%BranchNumber_pft        &
  )
!
!     ZNOON=hour of solar noon
!     iPlantCalendar(ipltcal_Emerge,=emergence date
!     iPlantPhenologyPattern=growth habit:0=annual,1=perennial from PFT file
!     iDayPlantHarvest,iYearPlantHarvest=day,year of harvesting
!     iYearCurrent=current year
!     iPlantBranchState=branch living flag: 0=alive,1=dead
!     GROUP=node number required for floral initiation
!     NodeNumberToInitFloral=node number at floral initiation
!     NodeNumberAtAnthesis=node number at flowering
!     NumOfLeaves_brch=number of leaves appeared
!     KLeafNodeNumber=integer of most recent leaf number currently growing
!     VSTGX=leaf number on date of floral initiation
!     TGSTGI=total change in vegve node number normalized for maturity group
!     TGSTGF=total change in reprve node number normalized for maturity group
!     FLG4=number of hours with no grain fill
!     doInitLeafOut=flag for initializing leafout
!     Hours4Leafout,VRNL=leafout hours,hours required for leafout
!     Hours4LeafOff,VRNX=leafoff hours,hours required for leafoff
!     HourCounter4LeafOut_brch=hourly leafout counter
!     RubiscoActivity_brpft,FDBKX=N,P feedback inhibition on C3 CO2 fixation
!     doInitLeafOut,doPlantLeafOut=flag for initializing,enabling leafout
!     doPlantLeaveOff=flag for enabling leafoff:0=enable,1=disable
!     IFLGQ=current hours after physl maturity until start of litterfall
!
  IF(J.EQ.INT(ZNOON).AND.iPlantCalendar(ipltcal_Emerge,NumOfMainBranch_pft(NZ),NZ).NE.0 &
    .AND.(iPlantPhenologyPattern(NZ).NE.iplt_annual.OR.(I.GE.iDayPlantHarvest(NZ) &
    .AND.iYearCurrent.GE.iYearPlantHarvest(NZ))))THEN
    IDTHY=0
!
!     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
!
    call LiterfallFromDeadBranches(I,J,NZ,IDTHY,CPOOLK)

    IF(IDTHY.EQ.NumOfBranches_pft(NZ))THEN
      iPlantShootState(NZ)=iDead
      BranchNumber_pft(NZ)=0
      HoursCanopyPSITooLow(NZ)=0._r8
      IF(doInitPlant(NZ).EQ.itrue)THEN
        NumOfBranches_pft(NZ)=1
      ELSE
        NumOfBranches_pft(NZ)=0
      ENDIF
      HypoctoylHeight(NZ)=0._r8
      VOLWOU=VOLWOU+CanWatP(NZ)
      UVOLO=UVOLO+CanWatP(NZ)
      CanWatP(NZ)=0._r8
!
!     RESET LIVING FLAGS
!
!     WTRVC,WTRT=PFT storage,root C
!     iPlantPhenologyPattern=growth habit:0=annual,1=perennial
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     PP=PFT population
!     iPlantShootState,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!
      IF(NonstructalChemElmnts_pft(ielmc,NZ).LT.1.0E-04_r8*RootChemElmnts_pft(ielmc,NZ).AND.iPlantPhenologyPattern(NZ).NE.iplt_annual)then
        iPlantRootState(NZ)=iDead
      endif
      IF(iPlantPhenologyPattern(NZ).EQ.iplt_annual)then
        iPlantRootState(NZ)=iDead
      endif
      IF(JHVST(NZ).NE.ihv_noaction)then
        iPlantRootState(NZ)=iDead
      endif
      IF(pftPlantPopulation(NZ).LE.0.0_r8)then
        iPlantRootState(NZ)=iDead
      endif
      IF(iPlantRootState(NZ).EQ.iDead)then
        iPlantShootState(NZ)=iDead
      endif
    ENDIF
!
!     DEAD ROOTS
!
!
!     LITTERFALL FROM DEAD ROOTS
!
    call LiterfallFromDeadRoots(I,J,NZ)
!
    call LiterfallFromRootShootStorage(I,J,NZ,CPOOLK)
  ENDIF
  end associate
  end subroutine ResetDeadBranch
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromRootShootStorage(I,J,NZ,CPOOLK)

  use EcoSIMCtrlDataType, only : iYearCurrent
  implicit none
  integer, intent(in) :: I,J,NZ
  REAL(R8),INTENT(INOUT) :: CPOOLK(NumOfCanopyLayers1,JP1)
  integer :: L,M,NR,NB,N,NE
!     begin_execution
  associate(                            &
    JHVST      =>   plt_distb%JHVST   , &
    iYearPlanting       =>   plt_distb%iYearPlanting    , &
    iDayPlanting      =>   plt_distb%iDayPlanting   , &
    EarChemElmnts_brch    =>   plt_biom%EarChemElmnts_brch  , &
    NonstructElmnt_brch     =>   plt_biom%NonstructElmnt_brch   , &
    StalkChemElmnts_brch    =>   plt_biom%StalkChemElmnts_brch  , &
    HuskChemElmnts_brch    =>   plt_biom%HuskChemElmnts_brch  , &
    PetioleChemElmnts_brch   =>   plt_biom%PetioleChemElmnts_brch , &
    WTNDBE     =>   plt_biom%WTNDBE   , &
    ReserveChemElmnts_brch    =>   plt_biom%ReserveChemElmnts_brch  , &
    NoduleNonstructElmnt_brch     =>   plt_biom%NoduleNonstructElmnt_brch   , &
    NonstructalChemElmnts_pft      =>   plt_biom%NonstructalChemElmnts_pft    , &
    GrainChemElmnts_brch     =>   plt_biom%GrainChemElmnts_brch   , &
    WTRT1E     =>   plt_biom%WTRT1E   , &
    LeafChemElmnts_brch    =>   plt_biom%LeafChemElmnts_brch  , &
     RootMycoNonstructElmnt_vr     =>   plt_biom% RootMycoNonstructElmnt_vr   , &
    StandingDeadKCompChemElmnts_pft     =>   plt_biom%StandingDeadKCompChemElmnts_pft   , &
    WTRT2E     =>   plt_biom%WTRT2E   , &
    FWODLE     =>   plt_allom%FWODLE  , &
    FWODBE     =>   plt_allom%FWODBE  , &
    FWOODE     =>   plt_allom%FWOODE  , &
    FWODRE     =>   plt_allom%FWODRE  , &
    doInitPlant      =>   plt_pheno%doInitPlant   , &
    iPlantPhenologyPattern     =>   plt_pheno%iPlantPhenologyPattern  , &
    iPlantPhenologyType     =>   plt_pheno%iPlantPhenologyType  , &
    iPlantRootState     =>   plt_pheno%iPlantRootState  , &
    iPlantMorphologyType     =>   plt_pheno%iPlantMorphologyType  , &
    iPlantTurnoverPattern     =>   plt_pheno%iPlantTurnoverPattern  , &
    iPlantShootState      =>   plt_pheno%iPlantShootState   , &
    icwood     =>   pltpar%icwood     , &
    ifoliar    =>   pltpar%ifoliar    , &
    k_fine_litr=>   pltpar%k_fine_litr, &
    k_woody_litr=>  pltpar%k_woody_litr, &
    LitterFallChemElmnt_pftvr       =>   plt_bgcr%LitterFallChemElmnt_pftvr     , &
    infoliar   =>   pltpar%infoliar   , &
    istalk     =>   pltpar%istalk     , &
    iroot      =>   pltpar%iroot      , &
    instruct   =>   pltpar%instruct   , &
    LYRC       =>   plt_site%LYRC     , &
    NJ         =>   plt_site%NJ       , &
    NU         =>   plt_site%NU       , &
    IDATA      =>   plt_site%IDATA    , &
    CFOPE      =>   plt_soilchem%CFOPE, &
    MY         =>   plt_morph%MY      , &
    NGTopRootLayer         =>   plt_morph%NGTopRootLayer      , &
    NumOfBranches_pft        =>   plt_morph%NumOfBranches_pft     , &
    NumRootAxes_pft       =>   plt_morph%NumRootAxes_pft      &
  )
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE
!     RESERVES FROM SHOOT AT DEATH
!
!     iPlantShootState,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!     doInitPlant=PFT initialization flag:0=no,1=yes
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CPOOLK=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
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
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
  IF(iPlantShootState(NZ).EQ.iDead.AND.iPlantRootState(NZ).EQ.iDead)THEN
    !both plant shoots and roots are dead
    IF(doInitPlant(NZ).EQ.ifalse)THEN
      D6425: DO M=1,jsken
        D8825: DO NB=1,NumOfBranches_pft(NZ)

          LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(ielmc,M,k_fine_litr,0,NZ) &
            +CFOPE(ielmc,instruct,M,NZ)*CPOOLK(NB,NZ)

          DO NE=1,NumOfPlantChemElmnts
            LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(LeafChemElmnts_brch(NE,NB,NZ)*FWODLE(NE,k_woody_litr) &
              +PetioleChemElmnts_brch(NE,NB,NZ)*FWODBE(NE,k_woody_litr))

            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,instruct,M,NZ)*(NonstructElmnt_brch(NE,NB,NZ)+NoduleNonstructElmnt_brch(NE,NB,NZ) &
              +ReserveChemElmnts_brch(NE,NB,NZ)) &
              +CFOPE(NE,ifoliar,M,NZ)*(LeafChemElmnts_brch(NE,NB,NZ)*FWODLE(NE,k_fine_litr) &
              +WTNDBE(NE,NB,NZ)) &
              +CFOPE(NE,infoliar,M,NZ)*(PetioleChemElmnts_brch(NE,NB,NZ)*FWODBE(NE,k_fine_litr) &
              +HuskChemElmnts_brch(NE,NB,NZ)+EarChemElmnts_brch(NE,NB,NZ))

            IF(iPlantPhenologyPattern(NZ).EQ.iplt_annual.AND.iPlantPhenologyType(NZ).NE.0)THEN
              NonstructalChemElmnts_pft(NE,NZ)=NonstructalChemElmnts_pft(NE,NZ)+CFOPE(NE,infoliar,M,NZ)*GrainChemElmnts_brch(NE,NB,NZ)
            ELSE
              LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,infoliar,M,NZ)*GrainChemElmnts_brch(NE,NB,NZ)
            ENDIF
            IF(iPlantTurnoverPattern(NZ).EQ.0.OR.iPlantMorphologyType(NZ).LE.1)THEN
              !all above ground
              LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,istalk,M,NZ)*StalkChemElmnts_brch(NE,NB,NZ)
            ELSE
              StandingDeadKCompChemElmnts_pft(NE,M,NZ)=StandingDeadKCompChemElmnts_pft(NE,M,NZ)+CFOPE(NE,icwood,M,NZ)*StalkChemElmnts_brch(NE,NB,NZ)
            ENDIF
          ENDDO
        ENDDO D8825
!
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE
!     RESERVES FROM ROOT AND STORGE AT DEATH
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!
        DO NE=1,NumOfPlantChemElmnts
          D6415: DO L=NU,NJ
            DO N=1,MY(NZ)
              LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ) &
                +CFOPE(NE,instruct,M,NZ)* RootMycoNonstructElmnt_vr(NE,N,L,NZ)
              DO NR=1,NumRootAxes_pft(NZ)
                LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)+CFOPE(NE,icwood,M,NZ) &
                  *(WTRT1E(NE,N,L,NR,NZ)+WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)

                LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
                  *(WTRT1E(NE,N,L,NR,NZ)+WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_fine_litr)
              ENDDO
            ENDDO
          ENDDO D6415
          LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ) &
            +CFOPE(NE,instruct,M,NZ)*NonstructalChemElmnts_pft(NE,NZ)*FWOODE(NE,k_woody_litr)

          LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ) &
            +CFOPE(NE,instruct,M,NZ)*NonstructalChemElmnts_pft(NE,NZ)*FWOODE(NE,k_fine_litr)
        ENDDO
      ENDDO D6425
!
      call ResetBranchRootStates(NZ,CPOOLK)
    ENDIF
!
!     RESEED DEAD PERENNIALS
!
!     iPlantPhenologyPattern=growth habit:0=annual,1=perennial from PFT file
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     LYRC=number of days in current year
!     iDayPlanting,iYearPlanting=day,year of planting
!
    IF(iPlantPhenologyPattern(NZ).NE.iplt_annual.AND.JHVST(NZ).EQ.ihv_noaction)THEN
      IF(I.LT.LYRC)THEN
        iDayPlanting(NZ)=I+1
        iYearPlanting(NZ)=iYearCurrent
      ELSE
        iDayPlanting(NZ)=1
        iYearPlanting(NZ)=iYearCurrent+1
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine LiterfallFromRootShootStorage
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromDeadRoots(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: L,M,NR,N,NE,NTG
!     begin_execution
  associate(                                &
    WTRT2E    =>   plt_biom%WTRT2E    , &
    RTWT1E    =>   plt_biom%RTWT1E    , &
    WTRT1E    =>   plt_biom%WTRT1E    , &
     RootMycoNonstructElmnt_vr    =>   plt_biom% RootMycoNonstructElmnt_vr    , &
     PopuPlantRootC_vr    =>   plt_biom% PopuPlantRootC_vr    , &
    RootStructBiomC_vr    =>   plt_biom%RootStructBiomC_vr    , &
    WSRTL     =>   plt_biom%WSRTL     , &
    RootNoduleNonstructElmnt_vr   =>   plt_biom%RootNoduleNonstructElmnt_vr   , &
    WTNDLE    =>   plt_biom%WTNDLE    , &
    FWODRE    =>   plt_allom%FWODRE   , &
    iPlantRootState    =>   plt_pheno%iPlantRootState   , &
    CFOPE     =>   plt_soilchem%CFOPE , &
    trcg_rootml      =>   plt_rbgc%trcg_rootml  , &
    trcs_rootml => plt_rbgc%trcs_rootml, &
    NJ        =>   plt_site%NJ        , &
    NU        =>   plt_site%NU        , &
    RootGasLoss_disturb     =>   plt_bgcr%RootGasLoss_disturb     , &
    LitterFallChemElmnt_pftvr      =>   plt_bgcr%LitterFallChemElmnt_pftvr      , &
    icwood    =>   pltpar%icwood      , &
    iroot     =>   pltpar%iroot       , &
    instruct   =>   pltpar%instruct   , &
    k_fine_litr=>  pltpar%k_fine_litr , &
    k_woody_litr=> pltpar%k_woody_litr, &
    RootLenDensNLP     =>   plt_morph%RootLenDensNLP    , &
    RTVLW     =>   plt_morph%RTVLW    , &
    RTARP     =>   plt_morph%RTARP    , &
    RTVLP     =>   plt_morph%RTVLP    , &
    SecndRootXNumL      =>   plt_morph%SecndRootXNumL     , &
    PrimRootXNumL     =>   plt_morph%PrimRootXNumL    , &
    PrimRootDepth    =>   plt_morph%PrimRootDepth   , &
    AveSecndRootLen     =>   plt_morph%AveSecndRootLen    , &
    PrimRootRadius     =>   plt_morph%PrimRootRadius    , &
    SecndRootRadius     =>   plt_morph%SecndRootRadius    , &
    MaxPrimRootRadius    =>   plt_morph%MaxPrimRootRadius   , &
    MaxSecndRootRadius    =>   plt_morph%MaxSecndRootRadius   , &
    RootLenPerP     =>   plt_morph%RootLenPerP    , &
    PrimRootLen     =>   plt_morph%PrimRootLen    , &
    SecndRootLen     =>   plt_morph%SecndRootLen    , &
    iPlantNfixType    =>   plt_morph%iPlantNfixType   , &
    RTN2      =>   plt_morph%RTN2     , &
    MY        =>   plt_morph%MY       , &
    NIXBotRootLayer      =>   plt_morph%NIXBotRootLayer     , &
    NINR      =>   plt_morph%NINR     , &
    SeedinDepth     =>   plt_morph%SeedinDepth    , &
    NGTopRootLayer        =>   plt_morph%NGTopRootLayer       , &
    NumRootAxes_pft      =>   plt_morph%NumRootAxes_pft       &
  )
!     IDTHR=PFT root living flag: 0=alive,1=dead
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
!     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
!

  IF(iPlantRootState(NZ).EQ.iDead)THEN
    !add root to litterfall
    D8900: DO N=1,MY(NZ)
      D8895: DO L=NU,NJ
        DO NE=1,NumOfPlantChemElmnts
          D6410: DO M=1,jsken
            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ) &
              +CFOPE(NE,instruct,M,NZ)* RootMycoNonstructElmnt_vr(NE,N,L,NZ)
            DO  NR=1,NumRootAxes_pft(NZ)
              LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,L,NZ)+CFOPE(NE,icwood,M,NZ) &
                *(WTRT1E(NE,N,L,NR,NZ)+WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)

              LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
                *(WTRT1E(NE,N,L,NR,NZ)+WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_fine_litr)
            enddo
          ENDDO D6410
        ENDDO
!
!     RELEASE GAS CONTENTS OF DEAD ROOTS
!
        DO NTG=idg_beg,idg_end-1
          RootGasLoss_disturb(NTG,NZ)=RootGasLoss_disturb(NTG,NZ)-trcg_rootml(NTG,N,L,NZ)-trcs_rootml(NTG,N,L,NZ)
        ENDDO
        trcg_rootml(idg_beg:idg_end-1,N,L,NZ)=0._r8
        trcs_rootml(idg_beg:idg_end-1,N,L,NZ)=0._r8
!
!     RESET STATE VARIABLES OF DEAD ROOTS
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     PrimRootLen,SecndRootLen=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     RootStructBiomC_vr,WTRTD=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,SecndRootXNumL=number of primary,secondary root axes
!     RootLenDensNLP,RootLenPerP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     AveSecndRootLen=average secondary root length
!
        D8870: DO NR=1,NumRootAxes_pft(NZ)
          WTRT1E(1:NumOfPlantChemElmnts,N,L,NR,NZ)=0._r8
          WTRT2E(1:NumOfPlantChemElmnts,N,L,NR,NZ)=0._r8
          RTWT1E(1:NumOfPlantChemElmnts,N,NR,NZ)=0._r8
          PrimRootLen(N,L,NR,NZ)=0._r8
          SecndRootLen(N,L,NR,NZ)=0._r8
          RTN2(N,L,NR,NZ)=0._r8
        ENDDO D8870
         RootMycoNonstructElmnt_vr(:,N,L,NZ)=0._r8
        RootStructBiomC_vr(N,L,NZ)=0._r8
         PopuPlantRootC_vr(N,L,NZ)=0._r8
        WSRTL(N,L,NZ)=0._r8
        PrimRootXNumL(N,L,NZ)=0._r8
        SecndRootXNumL(N,L,NZ)=0._r8
        RootLenPerP(N,L,NZ)=0._r8
        RootLenDensNLP(N,L,NZ)=0._r8
        RTVLP(N,L,NZ)=0._r8
        RTVLW(N,L,NZ)=0._r8
        PrimRootRadius(N,L,NZ)=MaxPrimRootRadius(N,NZ)
        SecndRootRadius(N,L,NZ)=MaxSecndRootRadius(N,NZ)
        RTARP(N,L,NZ)=0._r8
        AveSecndRootLen(N,L,NZ)=MinAve2ndRootLen
!
!     LITTERFALL AND STATE VARIABLES FROM DEAD NODULES
!
!     iPlantNfixType=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
        IF(iPlantNfixType(NZ).NE.0.AND.N.EQ.1)THEN
          D6420: DO M=1,jsken
            DO NE=1,NumOfPlantChemElmnts            
              LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
                *WTNDLE(NE,L,NZ)+CFOPE(NE,instruct,M,NZ)*RootNoduleNonstructElmnt_vr(NE,L,NZ)
            ENDDO
          ENDDO D6420
          WTNDLE(1:NumOfPlantChemElmnts,L,NZ)=0._r8
          RootNoduleNonstructElmnt_vr(1:NumOfPlantChemElmnts,L,NZ)=0._r8          
        ENDIF
      ENDDO D8895
    ENDDO D8900
!
!     RESET DEPTH VARIABLES OF DEAD ROOTS
!
!     NINR=deepest root layer
!     RTDP1=primary root depth from soil surface
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!
    D8795: DO NR=1,NumRootAxes_pft(NZ)
      NINR(NR,NZ)=NGTopRootLayer(NZ)
      D8790: DO N=1,MY(NZ)
        PrimRootDepth(N,NR,NZ)=SeedinDepth(NZ)
        RTWT1E(1:NumOfPlantChemElmnts,N,NR,NZ)=0._r8
      ENDDO D8790
    ENDDO D8795
    NIXBotRootLayer(NZ)=NGTopRootLayer(NZ)
    NumRootAxes_pft(NZ)=0
  ENDIF
  end associate
  end subroutine LiterfallFromDeadRoots
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromDeadBranches(I,J,NZ,IDTHY,CPOOLK)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: IDTHY
  real(r8), intent(inout) :: CPOOLK(NumOfCanopyLayers1,JP1)
  integer :: M,NE,NB
!     begin_execution
  associate(                                &
    IHVST     =>  plt_distb%IHVST     , &
    FWODBE    =>  plt_allom%FWODBE    , &
    FWODLE    =>  plt_allom%FWODLE    , &
    NoduleNonstructElmnt_brch    =>  plt_biom%NoduleNonstructElmnt_brch     , &
    StalkChemElmnts_brch   =>  plt_biom%StalkChemElmnts_brch    , &
    HuskChemElmnts_brch   =>  plt_biom%HuskChemElmnts_brch    , &
    PetioleChemElmnts_brch  =>  plt_biom%PetioleChemElmnts_brch   , &
    WTNDBE    =>  plt_biom%WTNDBE     , &
    LeafChemElmnts_brch   =>  plt_biom%LeafChemElmnts_brch    , &
    ZEROP     =>  plt_biom%ZEROP      , &
    EarChemElmnts_brch   =>  plt_biom%EarChemElmnts_brch    , &
    NonstructElmnt_brch    =>  plt_biom%NonstructElmnt_brch     , &
    GrainChemElmnts_brch    =>  plt_biom%GrainChemElmnts_brch     , &
    ReserveChemElmnts_brch   =>  plt_biom%ReserveChemElmnts_brch    , &
    StandingDeadKCompChemElmnts_pft    =>  plt_biom%StandingDeadKCompChemElmnts_pft     , &
    NonstructalChemElmnts_pft     =>  plt_biom%NonstructalChemElmnts_pft      , &
    CFOPE     =>  plt_soilchem%CFOPE  , &
    iPlantBranchState     =>  plt_pheno%iPlantBranchState     , &
    MatureGroup_brch    =>  plt_pheno%MatureGroup_brch    , &
    VSTGX     =>  plt_pheno%VSTGX     , &
    KLeafNodeNumber    =>  plt_pheno%KLeafNodeNumber    , &
    TGSTGI    =>  plt_pheno%TGSTGI    , &
    TGSTGF    =>  plt_pheno%TGSTGF    , &
    Hours4Leafout      =>  plt_pheno%Hours4Leafout      , &
    Hours4LeafOff      =>  plt_pheno%Hours4LeafOff      , &
    Hours4LenthenPhotoPeriod      =>  plt_pheno%Hours4LenthenPhotoPeriod      , &
    Hours4ShortenPhotoPeriod      =>  plt_pheno%Hours4ShortenPhotoPeriod      , &
    HourCounter4LeafOut_brch      =>  plt_pheno%HourCounter4LeafOut_brch      , &
    FLG4      =>  plt_pheno%FLG4      , &
    doInitLeafOut     =>  plt_pheno%doInitLeafOut     , &
    doPlantLeafOut     =>  plt_pheno%doPlantLeafOut     , &
    IFLGR     =>  plt_pheno%IFLGR     , &
    IFLGQ     =>  plt_pheno%IFLGQ     , &
    MatureGroup_pft   =>  plt_pheno%MatureGroup_pft   , &
    doPlantLeaveOff     =>  plt_pheno%doPlantLeaveOff     , &
    iPlantCalendar     =>  plt_pheno%iPlantCalendar     , &
    iPlantTurnoverPattern    =>  plt_pheno%iPlantTurnoverPattern    , &
    iPlantMorphologyType    =>  plt_pheno%iPlantMorphologyType    , &
    iPlantPhenologyType    =>  plt_pheno%iPlantPhenologyType    , &
    iPlantPhenologyPattern    =>  plt_pheno%iPlantPhenologyPattern    , &
    LitterFallChemElmnt_pftvr      =>  plt_bgcr%LitterFallChemElmnt_pftvr       , &
    NumOfBranches_pft       =>  plt_morph%NumOfBranches_pft       , &
    NodeNumberToInitFloral     =>  plt_morph%NodeNumberToInitFloral     , &
    ShootNodeNumber     =>  plt_morph%ShootNodeNumber     , &
    BranchNumber_brch      =>  plt_morph%BranchNumber_brch      , &
    NodeNumberAtAnthesis     =>  plt_morph%NodeNumberAtAnthesis     , &
    XTLI      =>  plt_morph%XTLI      , &
    NumOfLeaves_brch     =>  plt_morph%NumOfLeaves_brch     , &
    KLEAF     =>  plt_morph%KLEAF     , &
    istalk    =>  pltpar%istalk       , &
    k_fine_litr=>   pltpar%k_fine_litr, &
    k_woody_litr=> pltpar%k_woody_litr, &
    instruct   =>   pltpar%instruct   , &
    icwood     =>   pltpar%icwood     , &
    infoliar   =>  pltpar%infoliar    , &
    ifoliar    => pltpar%ifoliar      , &
    RubiscoActivity_brpft      =>  plt_photo%RubiscoActivity_brpft      , &
    FDBKX     =>  plt_photo%FDBKX       &
  )
  D8845: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState(NB,NZ).EQ.iDead)THEN
      MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
      ShootNodeNumber(NB,NZ)=XTLI(NZ)
      NodeNumberToInitFloral(NB,NZ)=ShootNodeNumber(NB,NZ)
      NodeNumberAtAnthesis(NB,NZ)=0._r8
      NumOfLeaves_brch(NB,NZ)=0._r8
      VSTGX(NB,NZ)=0._r8
      KLEAF(NB,NZ)=1
      KLeafNodeNumber(NB,NZ)=1
      TGSTGI(NB,NZ)=0._r8
      TGSTGF(NB,NZ)=0._r8
      Hours4Leafout(NB,NZ)=0._r8
      Hours4LeafOff(NB,NZ)=0._r8
      Hours4LenthenPhotoPeriod(NB,NZ)=0._r8
      Hours4ShortenPhotoPeriod(NB,NZ)=0._r8
      HourCounter4LeafOut_brch(NB,NZ)=0._r8
      FLG4(NB,NZ)=0._r8
      RubiscoActivity_brpft(NB,NZ)=1.0_r8
      FDBKX(NB,NZ)=1.0_r8
      doInitLeafOut(NB,NZ)=0
      doPlantLeafOut(NB,NZ)=iDisable
      doPlantLeaveOff(NB,NZ)=iEnable
      IFLGR(NB,NZ)=0
      IFLGQ(NB,NZ)=0
      BranchNumber_brch(NB,NZ)=0
      D8850: DO M=1,pltpar%NumGrowthStages
        iPlantCalendar(M,NB,NZ)=0
      ENDDO D8850
!
!     LITTERFALL FROM DEAD BRANCHES
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
!     WTEARB,WTEABN,WTEABP=ear C,N,P mass
!     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     iPlantPhenologyPattern=growth habit:0=annual,1=perennial from PFT file
!     iPlantPhenologyType=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     iPlantTurnoverPattern=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     iPlantMorphologyType=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
      D6405: DO M=1,jsken
        DO NE=1,NumOfPlantChemElmnts        
          LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
            +CFOPE(NE,instruct,M,NZ)*NoduleNonstructElmnt_brch(NE,NB,NZ) &
            +CFOPE(NE,ifoliar,M,NZ)*(LeafChemElmnts_brch(NE,NB,NZ)*FWODLE(NE,k_fine_litr) &
            +WTNDBE(NE,NB,NZ)) &
            +CFOPE(NE,infoliar,M,NZ)*(PetioleChemElmnts_brch(NE,NB,NZ)*FWODBE(NE,k_fine_litr) &
            +HuskChemElmnts_brch(NE,NB,NZ)+EarChemElmnts_brch(NE,NB,NZ))

          LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_woody_litr,0,NZ) &
            +CFOPE(NE,icwood,M,NZ)*(LeafChemElmnts_brch(NE,NB,NZ)*FWODBE(NE,k_woody_litr) &
            +PetioleChemElmnts_brch(NE,NB,NZ)*FWODBE(NE,k_woody_litr))

          IF(iPlantPhenologyPattern(NZ).EQ.iplt_annual.AND.iPlantPhenologyType(NZ).NE.0)THEN
            NonstructalChemElmnts_pft(NE,NZ)=NonstructalChemElmnts_pft(NE,NZ)+CFOPE(NE,infoliar,M,NZ)*GrainChemElmnts_brch(NE,NB,NZ)
          ELSE
            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,infoliar,M,NZ)*GrainChemElmnts_brch(NE,NB,NZ)
          ENDIF
          IF(iPlantTurnoverPattern(NZ).EQ.0.OR.iPlantMorphologyType(NZ).LE.1)THEN
            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,istalk,M,NZ)*StalkChemElmnts_brch(NE,NB,NZ)
          ELSE
            StandingDeadKCompChemElmnts_pft(NE,M,NZ)=StandingDeadKCompChemElmnts_pft(NE,M,NZ)+CFOPE(NE,icwood,M,NZ)*StalkChemElmnts_brch(NE,NB,NZ)
          ENDIF
        ENDDO
      ENDDO D6405

!
!     RECOVER NON-STRUCTURAL C,N,P FROM BRANCH TO
!     SEASONAL STORAGE RESERVES
!
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOOLK=total C4 nonstructural C in branch
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!
      NonstructalChemElmnts_pft(ielmc,NZ)=NonstructalChemElmnts_pft(ielmc,NZ)+CPOOLK(NB,NZ)
      DO NE=1,NumOfPlantChemElmnts
        NonstructalChemElmnts_pft(NE,NZ)=NonstructalChemElmnts_pft(NE,NZ)+NonstructElmnt_brch(NE,NB,NZ)
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          D6406: DO M=1,jsken
            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,instruct,M,NZ)*ReserveChemElmnts_brch(NE,NB,NZ)
          ENDDO D6406
        ELSE
          NonstructalChemElmnts_pft(NE,NZ)=NonstructalChemElmnts_pft(NE,NZ)+ReserveChemElmnts_brch(NE,NB,NZ)
        ENDIF
      ENDDO
!
      call ResetDeadRootStates(NB,NZ,CPOOLK)

      IDTHY=IDTHY+1
    ENDIF
  ENDDO D8845
  end associate
  end subroutine LiterfallFromDeadBranches

!------------------------------------------------------------------------------------------

  subroutine ResetBranchRootStates(NZ,CPOOLK)
  implicit none
  integer, intent(in) :: NZ
  real(r8),INTENT(OUT) :: CPOOLK(NumOfCanopyLayers1,JP1)
  integer :: L,NR,N,NE,NB
!     begin_execution
  associate(                           &
    NonstructElmnt_brch => plt_biom%NonstructElmnt_brch        , &
    NoduleNonstructElmnt_brch => plt_biom%NoduleNonstructElmnt_brch        , &
    WTSHTBE=> plt_biom%WTSHTBE       , &
    LeafChemElmnts_brch=> plt_biom%LeafChemElmnts_brch       , &
    ReserveChemElmnts_brch=> plt_biom%ReserveChemElmnts_brch       , &
    PetioleChemElmnts_brch=> plt_biom%PetioleChemElmnts_brch      , &
    HuskChemElmnts_brch=> plt_biom%HuskChemElmnts_brch       , &
    EarChemElmnts_brch=> plt_biom%EarChemElmnts_brch       , &
    GrainChemElmnts_brch => plt_biom%GrainChemElmnts_brch        , &
    StalkBiomassC_brch => plt_biom%StalkBiomassC_brch        , &
     RootMycoNonstructElmnt_vr => plt_biom% RootMycoNonstructElmnt_vr        , &
    WTSTXBE=> plt_biom%WTSTXBE       , &
    LeafPetioleBiomassC_brch  => plt_biom%LeafPetioleBiomassC_brch         , &
    NonstructalChemElmnts_pft  => plt_biom%NonstructalChemElmnts_pft         , &
    WTRT1E => plt_biom%WTRT1E        , &
    WTRT2E => plt_biom%WTRT2E        , &
    WTNDBE => plt_biom%WTNDBE        , &
    StalkChemElmnts_brch=> plt_biom%StalkChemElmnts_brch       , &
    RTWT1E => plt_biom%RTWT1E        , &
    NJ     => plt_site%NJ            , &
    NU     => plt_site%NU            , &
    iPlantState  => plt_pheno%iPlantState        , &
    SecndRootLen  => plt_morph%SecndRootLen        , &
    RTN2   => plt_morph%RTN2         , &
    PrimRootLen  => plt_morph%PrimRootLen        , &
    MY     => plt_morph%MY           , &
    NumOfBranches_pft    => plt_morph%NumOfBranches_pft          , &
    NumRootAxes_pft   => plt_morph%NumRootAxes_pft           &
  )
!     RESET BRANCH STATE VARIABLES
!
  DO NE=1,NumOfPlantChemElmnts
    DO NB=1,NumOfBranches_pft(NZ)
      NonstructElmnt_brch(NE,NB,NZ)=0._r8
      NoduleNonstructElmnt_brch(NE,NB,NZ)=0._r8
      WTSHTBE(NE,NB,NZ)=0._r8
      LeafChemElmnts_brch(NE,NB,NZ)=0._r8
      PetioleChemElmnts_brch(NE,NB,NZ)=0._r8
      StalkChemElmnts_brch(NE,NB,NZ)=0._r8
      ReserveChemElmnts_brch(NE,NB,NZ)=0._r8
    ENDDO
  ENDDO
  D8835: DO NB=1,NumOfBranches_pft(NZ)
    CPOOLK(NB,NZ)=0._r8
    StalkBiomassC_brch(NB,NZ)=0._r8
    WTNDBE(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
    HuskChemElmnts_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
    EarChemElmnts_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
    GrainChemElmnts_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
    LeafPetioleBiomassC_brch(NB,NZ)=0._r8
    WTSTXBE(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  ENDDO D8835
!
!     RESET ROOT STATE VARIABLES
!
  D6416: DO L=NU,NJ
    DO  N=1,MY(NZ)
       RootMycoNonstructElmnt_vr(1:NumOfPlantChemElmnts,N,L,NZ)=0._r8
      DO  NR=1,NumRootAxes_pft(NZ)
        WTRT1E(1:NumOfPlantChemElmnts,N,L,NR,NZ)=0._r8
        WTRT2E(1:NumOfPlantChemElmnts,N,L,NR,NZ)=0._r8
        RTWT1E(1:NumOfPlantChemElmnts,N,NR,NZ)=0._r8
        PrimRootLen(N,L,NR,NZ)=0._r8
        SecndRootLen(N,L,NR,NZ)=0._r8
        RTN2(N,L,NR,NZ)=0._r8
      enddo
    enddo
  ENDDO D6416
  NonstructalChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  iPlantState(NZ)=1
  end associate
  end subroutine ResetBranchRootStates

!------------------------------------------------------------------------------------------

  subroutine ResetDeadRootStates(NB,NZ,CPOOLK)
!     RESET STATE VARIABLES FROM DEAD BRANCHES
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8),intent(inout) :: CPOOLK(NumOfCanopyLayers1,JP1)
  integer :: L,K,N
!     begin_execution
  associate(                          &
    NonstructElmnt_brch   => plt_biom%NonstructElmnt_brch     , &
    NoduleNonstructElmnt_brch   => plt_biom%NoduleNonstructElmnt_brch     , &
    WTSHTBE  => plt_biom%WTSHTBE    , &
    LeafChemElmnts_brch  => plt_biom%LeafChemElmnts_brch    , &
    WTNDBE   => plt_biom%WTNDBE     , &
    StalkChemElmnts_brch  => plt_biom%StalkChemElmnts_brch    , &
    LeafChemElmntNode_brch    => plt_biom%LeafChemElmntNode_brch      , &
    ReserveChemElmnts_brch  => plt_biom%ReserveChemElmnts_brch    , &
    WTSTXBE  => plt_biom%WTSTXBE    , &
    StalkBiomassC_brch   => plt_biom%StalkBiomassC_brch     , &
    GrainChemElmnts_brch   => plt_biom%GrainChemElmnts_brch     , &
    LeafPetioleBiomassC_brch    => plt_biom%LeafPetioleBiomassC_brch      , &
    HuskChemElmnts_brch  => plt_biom%HuskChemElmnts_brch    , &
    EarChemElmnts_brch  => plt_biom%EarChemElmnts_brch    , &
    PetioleProteinCNode_brch   => plt_biom%PetioleProteinCNode_brch     , &
    LeafProteinCNode_brch     => plt_biom%LeafProteinCNode_brch       , &
    InternodeChemElmnt_brch   => plt_biom%InternodeChemElmnt_brch     , &
    PetioleElmntNode_brch   => plt_biom%PetioleElmntNode_brch     , &
    WGLFLE   => plt_biom%WGLFLE     , &
    CanopyLeafCpft_lyr    => plt_biom%CanopyLeafCpft_lyr      , &
    PetioleChemElmnts_brch => plt_biom%PetioleChemElmnts_brch   , &
    CPOOL3   => plt_photo%CPOOL3    , &
    CPOOL4   => plt_photo%CPOOL4    , &
    CO2B     => plt_photo%CO2B      , &
    HCOB     => plt_photo%HCOB      , &
    GRWTB    => plt_allom%GRWTB     , &
    GRNXB    => plt_morph%GRNXB     , &
    GRNOB    => plt_morph%GRNOB     , &
    LeafAreaLive_brch    => plt_morph%LeafAreaLive_brch     , &
    CanopyBranchStemApft_lyr    => plt_morph%CanopyBranchStemApft_lyr     , &
    NumOfBranches_pft      => plt_morph%NumOfBranches_pft       , &
    LeafAreaNode_brch    => plt_morph%LeafAreaNode_brch     , &
    PetioleLengthNode_brch    => plt_morph%PetioleLengthNode_brch     , &
    InternodeHeightDying_brch   => plt_morph%InternodeHeightDying_brch    , &
    LeafA_lyrnodbrchpft     => plt_morph%LeafA_lyrnodbrchpft      , &
    InternodeHeightLive_brch   => plt_morph%InternodeHeightLive_brch    , &
    CanopyLeafApft_lyr    => plt_morph%CanopyLeafApft_lyr     , &
    StemA_lyrnodbrchpft    => plt_morph%StemA_lyrnodbrchpft     , &
    CanPLNBLA    => plt_morph%CanPLNBLA       &
  )
!
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
!     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
!     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
!     WTEARB,WTEABN,WTEABP=ear C,N,P mass
!     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!     iPlantPhenologyPattern=growth habit:0=annual,1=perennial from PFT file
!     GRNOB=seed set number
!     GRNXB=potential number of seed set sites
!     GRWTB=individual seed size
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     LeafProteinCNode_brch=leaf protein mass
!     LeafAreaLive_brch=branch leaf area
!     WGLF,WGLFN,WGLFP,LeafProteinCNode_brch=node leaf C,N,P,protein mass
!     PetioleElmntNode_brch,WGSHN,WGSHP,PetioleProteinCNode_brch=node petiole C,N,P,protein mass
!     InternodeChemElmnt_brch,WGNODN,WGNODP=node stalk C,N,P mass
!
  CPOOLK(NB,NZ)=0._r8
  NonstructElmnt_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  NoduleNonstructElmnt_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  WTSHTBE(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  LeafChemElmnts_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  WTNDBE(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  PetioleChemElmnts_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  StalkChemElmnts_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  ReserveChemElmnts_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  HuskChemElmnts_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  EarChemElmnts_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  GrainChemElmnts_brch(1:NumOfPlantChemElmnts,NB,NZ)=0._r8
  StalkBiomassC_brch(NB,NZ)=0._r8
  LeafPetioleBiomassC_brch(NB,NZ)=0._r8
  GRNXB(NB,NZ)=0._r8
  GRNOB(NB,NZ)=0._r8
  GRWTB(NB,NZ)=0._r8
  LeafAreaLive_brch(NB,NZ)=0._r8
  WTSTXBE(1:NumOfPlantChemElmnts,NB,NZ)=0._r8

  D8855: DO K=0,MaxNodesPerBranch1
    IF(K.NE.0)THEN
      CPOOL3(K,NB,NZ)=0._r8
      CPOOL4(K,NB,NZ)=0._r8
      CO2B(K,NB,NZ)=0._r8
      HCOB(K,NB,NZ)=0._r8
    ENDIF
    LeafAreaNode_brch(K,NB,NZ)=0._r8
    InternodeHeightLive_brch(K,NB,NZ)=0._r8
    InternodeHeightDying_brch(K,NB,NZ)=0._r8
    PetioleLengthNode_brch(K,NB,NZ)=0._r8
    LeafProteinCNode_brch(K,NB,NZ)=0._r8
    PetioleProteinCNode_brch(K,NB,NZ)=0._r8
    LeafChemElmntNode_brch(1:NumOfPlantChemElmnts,K,NB,NZ)=0._r8
    PetioleElmntNode_brch(1:NumOfPlantChemElmnts,K,NB,NZ)=0._r8
    InternodeChemElmnt_brch(1:NumOfPlantChemElmnts,K,NB,NZ)=0._r8
    D8865: DO L=1,NumOfCanopyLayers1
      CanopyLeafApft_lyr(L,NZ)=CanopyLeafApft_lyr(L,NZ)-CanPLNBLA(L,K,NB,NZ)
      CanopyLeafCpft_lyr(L,NZ)=CanopyLeafCpft_lyr(L,NZ)-WGLFLE(ielmc,L,K,NB,NZ)
      CanPLNBLA(L,K,NB,NZ)=0._r8
      WGLFLE(1:NumOfPlantChemElmnts,L,K,NB,NZ)=0._r8
      IF(K.NE.0)THEN
        D8860: DO N=1,JLI1
          LeafA_lyrnodbrchpft(N,L,K,NB,NZ)=0._r8
        ENDDO D8860
      ENDIF
    ENDDO D8865
  ENDDO D8855
  D8875: DO L=1,NumOfCanopyLayers1
    CanopyBranchStemApft_lyr(L,NB,NZ)=0._r8
    DO  N=1,JLI1
      StemA_lyrnodbrchpft(N,L,NB,NZ)=0._r8
    enddo
  ENDDO D8875
  end associate
  end subroutine ResetDeadRootStates


end module LitterFallMod
