module LitterFallMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
implicit none
  private
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: ResetDeadBranch
  contains
!------------------------------------------------------------------------------------------

  subroutine ResetDeadBranch(I,J,NZ,CPOOLK)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8),intent(inout) :: CPOOLK(JC1,JP1)
  integer :: IDTHY

!     begin_execution
  associate(                                 &
    IYRH       =>   plt_distb%IYRH     , &
    JHVST      =>   plt_distb%JHVST    , &
    IDAYH      =>   plt_distb%IDAYH    , &
    UVOLO      =>   plt_ew%UVOLO       , &
    CanWatP      =>   plt_ew%CanWatP       , &
    WTRTE      =>   plt_biom%WTRTE     , &
    WTRVE      =>   plt_biom%WTRVE     , &
    ISTYP      =>   plt_pheno%ISTYP    , &
    IDTHR      =>   plt_pheno%IDTHR    , &
    IDTHP      =>   plt_pheno%IDTHP    , &
    IFLGI      =>   plt_pheno%IFLGI    , &
    WSTR       =>   plt_pheno%WSTR     , &
    IDAY       =>   plt_pheno%IDAY     , &
    pftPlantPopulation         =>   plt_site%pftPlantPopulation        , &
    IYRC       =>   plt_site%IYRC      , &
    VOLWOU     =>   plt_site%VOLWOU    , &
    ZNOON      =>   plt_site%ZNOON     , &
    HypoctoylHeight      =>   plt_morph%HypoctoylHeight    , &
    NBR        =>   plt_morph%NBR      , &
    NB1        =>   plt_morph%NB1      , &
    NBT        =>   plt_morph%NBT        &
  )
!
!     ZNOON=hour of solar noon
!     IDAY(1,=emergence date
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IDAYH,IYRH=day,year of harvesting
!     IYRC=current year
!     IDTHB=branch living flag: 0=alive,1=dead
!     GROUP=node number required for floral initiation
!     PSTGI=node number at floral initiation
!     PSTGF=node number at flowering
!     VSTG=number of leaves appeared
!     KVSTG=integer of most recent leaf number currently growing
!     VSTGX=leaf number on date of floral initiation
!     TGSTGI=total change in vegve node number normalized for maturity group
!     TGSTGF=total change in reprve node number normalized for maturity group
!     FLG4=number of hours with no grain fill
!     IFLGA=flag for initializing leafout
!     VRNS,VRNL=leafout hours,hours required for leafout
!     VRNF,VRNX=leafoff hours,hours required for leafoff
!     ATRP=hourly leafout counter
!     FDBK,FDBKX=N,P feedback inhibition on C3 CO2 fixation
!     IFLGA,IFLGE=flag for initializing,enabling leafout
!     IFLGF=flag for enabling leafoff:0=enable,1=disable
!     IFLGQ=current hours after physl maturity until start of litterfall
!
  IF(J.EQ.INT(ZNOON).AND.IDAY(1,NB1(NZ),NZ).NE.0 &
    .AND.(ISTYP(NZ).NE.iplt_annual.OR.(I.GE.IDAYH(NZ) &
    .AND.IYRC.GE.IYRH(NZ))))THEN
    IDTHY=0
!
!     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
!
    call LiterfallFromDeadBranches(I,J,NZ,IDTHY,CPOOLK)

    IF(IDTHY.EQ.NBR(NZ))THEN
      IDTHP(NZ)=ibrdead
      NBT(NZ)=0
      WSTR(NZ)=0._r8
      IF(IFLGI(NZ).EQ.1)THEN
        NBR(NZ)=1
      ELSE
        NBR(NZ)=0
      ENDIF
      HypoctoylHeight(NZ)=0._r8
      VOLWOU=VOLWOU+CanWatP(NZ)
      UVOLO=UVOLO+CanWatP(NZ)
      CanWatP(NZ)=0._r8
!
!     RESET LIVING FLAGS
!
!     WTRVC,WTRT=PFT storage,root C
!     ISTYP=growth habit:0=annual,1=perennial
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     PP=PFT population
!     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!
      IF(WTRVE(ielmc,NZ).LT.1.0E-04_r8*WTRTE(ielmc,NZ).AND.ISTYP(NZ).NE.iplt_annual)then
        IDTHR(NZ)=ibrdead
      endif  
      IF(ISTYP(NZ).EQ.iplt_annual)then
        IDTHR(NZ)=ibrdead
      endif  
      IF(JHVST(NZ).NE.ihv_noaction)then
        IDTHR(NZ)=ibrdead
      endif  
      IF(pftPlantPopulation(NZ).LE.0.0)then
        IDTHR(NZ)=ibrdead
      endif  
      IF(IDTHR(NZ).EQ.ibrdead)then
        IDTHP(NZ)=ibrdead
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
  
  use EcoSIMCtrlDataType, only : iyear_cur
  implicit none
  integer, intent(in) :: I,J,NZ
  REAL(R8),INTENT(INOUT) :: CPOOLK(JC1,JP1)
  integer :: L,M,NR,NB,N,NE
!     begin_execution
  associate(                            &
    JHVST      =>   plt_distb%JHVST   , &
    IYR0       =>   plt_distb%IYR0    , &
    IDAY0      =>   plt_distb%IDAY0   , &
    WTEARBE    =>   plt_biom%WTEARBE  , &
    EPOOL      =>   plt_biom%EPOOL    , &
    WTSTKBE    =>   plt_biom%WTSTKBE  , &
    WTHSKBE    =>   plt_biom%WTHSKBE  , &
    WTSHEBE    =>   plt_biom%WTSHEBE  , &
    WTNDBE     =>   plt_biom%WTNDBE   , &
    WTRSVBE    =>   plt_biom%WTRSVBE  , &
    EPOLNB     =>   plt_biom%EPOLNB   , &
    WTRVE      =>   plt_biom%WTRVE    , &
    WTGRBE     =>   plt_biom%WTGRBE   , &
    WTRT1E     =>   plt_biom%WTRT1E   , &
    WTLFBE     =>   plt_biom%WTLFBE   , &
    EPOOLR     =>   plt_biom%EPOOLR   , &
    WTSTDE     =>   plt_biom%WTSTDE   , &
    WTRT2E     =>   plt_biom%WTRT2E   , &
    FWODLE     =>   plt_allom%FWODLE  , &
    FWODBE     =>   plt_allom%FWODBE  , &
    FWOODE     =>   plt_allom%FWOODE  , &
    FWODRE     =>   plt_allom%FWODRE  , &
    IFLGI      =>   plt_pheno%IFLGI   , &
    ISTYP      =>   plt_pheno%ISTYP   , &
    IWTYP      =>   plt_pheno%IWTYP   , &
    IDTHR      =>   plt_pheno%IDTHR   , &
    IGTYP      =>   plt_pheno%IGTYP   , &
    IBTYP      =>   plt_pheno%IBTYP   , &
    IDTHP      =>   plt_pheno%IDTHP   , &
    icwood     =>   pltpar%icwood     , &
    ifoliar    =>   pltpar%ifoliar    , &
    k_fine_litr=>   pltpar%k_fine_litr, &
    k_woody_litr=>  pltpar%k_woody_litr, &
    ESNC       =>   plt_bgcr%ESNC     , &
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
    NBR        =>   plt_morph%NBR     , &
    NRT        =>   plt_morph%NRT       &
  )
!     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE
!     RESERVES FROM SHOOT AT DEATH
!
!     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!     IFLGI=PFT initialization flag:0=no,1=yes
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
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     WTRVC,WTRVN,WTRVP=storage C,N,P
!     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
!
  IF(IDTHP(NZ).EQ.1.AND.IDTHR(NZ).EQ.1)THEN
    IF(IFLGI(NZ).EQ.0)THEN
      D6425: DO M=1,jsken
        D8825: DO NB=1,NBR(NZ)

          ESNC(ielmc,M,k_fine_litr,0,NZ)=ESNC(ielmc,M,k_fine_litr,0,NZ) &
            +CFOPE(ielmc,instruct,M,NZ)*CPOOLK(NB,NZ)

          DO NE=1,npelms
            ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ) &
              +CFOPE(NE,icwood,M,NZ)*(WTLFBE(NE,NB,NZ)*FWODLE(NE,k_woody_litr) &
              +WTSHEBE(NE,NB,NZ)*FWODBE(NE,k_woody_litr))

            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,instruct,M,NZ)*(EPOOL(NE,NB,NZ)+EPOLNB(NE,NB,NZ) &
              +WTRSVBE(NE,NB,NZ)) &
              +CFOPE(NE,ifoliar,M,NZ)*(WTLFBE(NE,NB,NZ)*FWODLE(NE,k_fine_litr) &
              +WTNDBE(NE,NB,NZ)) &
              +CFOPE(NE,infoliar,M,NZ)*(WTSHEBE(NE,NB,NZ)*FWODBE(NE,k_fine_litr) &
              +WTHSKBE(NE,NB,NZ)+WTEARBE(NE,NB,NZ))

            IF(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).NE.0)THEN
              WTRVE(NE,NZ)=WTRVE(NE,NZ)+CFOPE(NE,infoliar,M,NZ)*WTGRBE(NE,NB,NZ)
            ELSE
              ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,infoliar,M,NZ)*WTGRBE(NE,NB,NZ)
            ENDIF
            IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
!all above ground
              ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ)+CFOPE(NE,istalk,M,NZ)*WTSTKBE(NE,NB,NZ)
            ELSE
              WTSTDE(NE,M,NZ)=WTSTDE(NE,M,NZ)+CFOPE(NE,icwood,M,NZ)*WTSTKBE(NE,NB,NZ)
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
        DO NE=1,npelms
          D6415: DO L=NU,NJ
            DO N=1,MY(NZ)
              ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ) &
                +CFOPE(NE,instruct,M,NZ)*EPOOLR(NE,N,L,NZ)
              DO NR=1,NRT(NZ)
                ESNC(NE,M,k_woody_litr,L,NZ)=ESNC(NE,M,k_woody_litr,L,NZ)+CFOPE(NE,icwood,M,NZ) &
                  *(WTRT1E(NE,N,L,NR,NZ)+WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)

                ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
                  *(WTRT1E(NE,N,L,NR,NZ)+WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_fine_litr)
              ENDDO
            ENDDO
          ENDDO D6415
          ESNC(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ)=ESNC(NE,M,k_woody_litr,NGTopRootLayer(NZ),NZ) &
            +CFOPE(NE,instruct,M,NZ)*WTRVE(NE,NZ)*FWOODE(NE,k_woody_litr)

          ESNC(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ)=ESNC(NE,M,k_fine_litr,NGTopRootLayer(NZ),NZ) &
            +CFOPE(NE,instruct,M,NZ)*WTRVE(NE,NZ)*FWOODE(NE,k_fine_litr)
        ENDDO
      ENDDO D6425
!
      call ResetBranchRootStates(NZ,CPOOLK)
    ENDIF
!
!     RESEED DEAD PERENNIALS
!
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     LYRC=number of days in current year
!     IDAY0,IYR0=day,year of planting
!
    IF(ISTYP(NZ).NE.iplt_annual.AND.JHVST(NZ).EQ.ihv_noaction)THEN
      IF(I.LT.LYRC)THEN
        IDAY0(NZ)=I+1
        IYR0(NZ)=iyear_cur
      ELSE
        IDAY0(NZ)=1
        IYR0(NZ)=iyear_cur+1
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
    EPOOLR    =>   plt_biom%EPOOLR    , &
    RootCPZR    =>   plt_biom%RootCPZR    , &
    WTRTL     =>   plt_biom%WTRTL     , &
    WSRTL     =>   plt_biom%WSRTL     , &
    EPOOLN    =>   plt_biom%EPOOLN    , &
    WTNDLE    =>   plt_biom%WTNDLE    , &
    FWODRE    =>   plt_allom%FWODRE   , &
    IDTHR     =>   plt_pheno%IDTHR    , &
    CFOPE     =>   plt_soilchem%CFOPE , &
    trcg_rootml      =>   plt_rbgc%trcg_rootml  , &
    trcs_rootml => plt_rbgc%trcs_rootml, &
    NJ        =>   plt_site%NJ        , &
    NU        =>   plt_site%NU        , &
    RootGasLoss_disturb     =>   plt_bgcr%RootGasLoss_disturb     , &
    ESNC      =>   plt_bgcr%ESNC      , &
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
    INTYP     =>   plt_morph%INTYP    , &
    RTN2      =>   plt_morph%RTN2     , &
    MY        =>   plt_morph%MY       , &
    NIXBotRootLayer      =>   plt_morph%NIXBotRootLayer     , &
    NINR      =>   plt_morph%NINR     , &
    SeedinDepth     =>   plt_morph%SeedinDepth    , &
    NGTopRootLayer        =>   plt_morph%NGTopRootLayer       , &
    NRT       =>   plt_morph%NRT        &
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

  IF(IDTHR(NZ).EQ.1)THEN
    D8900: DO N=1,MY(NZ)

      D8895: DO L=NU,NJ
        DO NE=1,npelms
          D6410: DO M=1,jsken
            ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ) &
              +CFOPE(NE,instruct,M,NZ)*EPOOLR(NE,N,L,NZ)
            DO  NR=1,NRT(NZ)
              ESNC(NE,M,k_woody_litr,L,NZ)=ESNC(NE,M,k_woody_litr,L,NZ)+CFOPE(NE,icwood,M,NZ) &
                *(WTRT1E(NE,N,L,NR,NZ)+WTRT2E(NE,N,L,NR,NZ))*FWODRE(NE,k_woody_litr)

              ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
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
!     WTRTL,WTRTD=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,SecndRootXNumL=number of primary,secondary root axes
!     RootLenDensNLP,RootLenPerP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     AveSecndRootLen=average secondary root length
!
        D8870: DO NR=1,NRT(NZ)
          WTRT1E(1:npelms,N,L,NR,NZ)=0._r8
          WTRT2E(1:npelms,N,L,NR,NZ)=0._r8
          RTWT1E(1:npelms,N,NR,NZ)=0._r8
          PrimRootLen(N,L,NR,NZ)=0._r8
          SecndRootLen(N,L,NR,NZ)=0._r8
          RTN2(N,L,NR,NZ)=0._r8
        ENDDO D8870
        EPOOLR(:,N,L,NZ)=0._r8
        WTRTL(N,L,NZ)=0._r8
        RootCPZR(N,L,NZ)=0._r8
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
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
        IF(INTYP(NZ).NE.0.AND.N.EQ.1)THEN
          D6420: DO M=1,jsken
            DO NE=1,npelms            
              ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)+CFOPE(NE,iroot,M,NZ) &
                *WTNDLE(NE,L,NZ)+CFOPE(NE,instruct,M,NZ)*EPOOLN(NE,L,NZ)
            ENDDO    
          ENDDO D6420
          WTNDLE(1:npelms,L,NZ)=0._r8
          EPOOLN(1:npelms,L,NZ)=0._r8          
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
    D8795: DO NR=1,NRT(NZ)
      NINR(NR,NZ)=NGTopRootLayer(NZ)
      D8790: DO N=1,MY(NZ)
        PrimRootDepth(N,NR,NZ)=SeedinDepth(NZ)
        RTWT1E(1:npelms,N,NR,NZ)=0._r8
      ENDDO D8790
    ENDDO D8795
    NIXBotRootLayer(NZ)=NGTopRootLayer(NZ)
    NRT(NZ)=0
  ENDIF
  end associate
  end subroutine LiterfallFromDeadRoots
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromDeadBranches(I,J,NZ,IDTHY,CPOOLK)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer, intent(inout) :: IDTHY
  real(r8), intent(inout) :: CPOOLK(JC1,JP1)
  integer :: M,NE,NB
!     begin_execution
  associate(                                &
    IHVST     =>  plt_distb%IHVST     , &
    FWODBE    =>  plt_allom%FWODBE    , &
    FWODLE    =>  plt_allom%FWODLE    , &
    EPOLNB    =>  plt_biom%EPOLNB     , &
    WTSTKBE   =>  plt_biom%WTSTKBE    , &
    WTHSKBE   =>  plt_biom%WTHSKBE    , &
    WTSHEBE   =>  plt_biom%WTSHEBE    , &
    WTNDBE    =>  plt_biom%WTNDBE     , &
    WTLFBE    =>  plt_biom%WTLFBE     , &
    ZEROP     =>  plt_biom%ZEROP      , &
    WTEARBE   =>  plt_biom%WTEARBE    , &
    EPOOL     =>  plt_biom%EPOOL      , &
    WTGRBE    =>  plt_biom%WTGRBE     , &
    WTRSVBE   =>  plt_biom%WTRSVBE    , &
    WTSTDE    =>  plt_biom%WTSTDE     , &
    WTRVE     =>  plt_biom%WTRVE      , &
    CFOPE     =>  plt_soilchem%CFOPE  , &
    IDTHB     =>  plt_pheno%IDTHB     , &
    GROUP     =>  plt_pheno%GROUP     , &
    VSTGX     =>  plt_pheno%VSTGX     , &
    KVSTG     =>  plt_pheno%KVSTG     , &
    TGSTGI    =>  plt_pheno%TGSTGI    , &
    TGSTGF    =>  plt_pheno%TGSTGF    , &
    VRNS      =>  plt_pheno%VRNS      , &
    VRNF      =>  plt_pheno%VRNF      , &
    VRNY      =>  plt_pheno%VRNY      , &
    VRNZ      =>  plt_pheno%VRNZ      , &
    ATRP      =>  plt_pheno%ATRP      , &
    FLG4      =>  plt_pheno%FLG4      , &
    IFLGA     =>  plt_pheno%IFLGA     , &
    IFLGE     =>  plt_pheno%IFLGE     , &
    IFLGR     =>  plt_pheno%IFLGR     , &
    IFLGQ     =>  plt_pheno%IFLGQ     , &
    GROUPI    =>  plt_pheno%GROUPI    , &
    IFLGF     =>  plt_pheno%IFLGF     , &
    IDAY      =>  plt_pheno%IDAY      , &
    IBTYP     =>  plt_pheno%IBTYP     , &
    IGTYP     =>  plt_pheno%IGTYP     , &
    IWTYP     =>  plt_pheno%IWTYP     , &
    ISTYP     =>  plt_pheno%ISTYP     , &
    ESNC      =>  plt_bgcr%ESNC       , &
    NBR       =>  plt_morph%NBR       , &
    PSTGI     =>  plt_morph%PSTGI     , &
    PSTG      =>  plt_morph%PSTG      , &
    NBTB      =>  plt_morph%NBTB      , &
    PSTGF     =>  plt_morph%PSTGF     , &
    XTLI      =>  plt_morph%XTLI      , &
    VSTG      =>  plt_morph%VSTG      , &
    KLEAF     =>  plt_morph%KLEAF     , &
    istalk    =>  pltpar%istalk       , &
    k_fine_litr=>   pltpar%k_fine_litr, &
    k_woody_litr=> pltpar%k_woody_litr, &
    instruct   =>   pltpar%instruct   , &
    icwood     =>   pltpar%icwood     , &
    infoliar   =>  pltpar%infoliar    , &
    ifoliar    => pltpar%ifoliar      , &
    FDBK      =>  plt_photo%FDBK      , &
    FDBKX     =>  plt_photo%FDBKX       &
  )
  D8845: DO NB=1,NBR(NZ)
    IF(IDTHB(NB,NZ).EQ.ibrdead)THEN
      GROUP(NB,NZ)=GROUPI(NZ)
      PSTG(NB,NZ)=XTLI(NZ)
      PSTGI(NB,NZ)=PSTG(NB,NZ)
      PSTGF(NB,NZ)=0._r8
      VSTG(NB,NZ)=0._r8
      VSTGX(NB,NZ)=0._r8
      KLEAF(NB,NZ)=1
      KVSTG(NB,NZ)=1
      TGSTGI(NB,NZ)=0._r8
      TGSTGF(NB,NZ)=0._r8
      VRNS(NB,NZ)=0._r8
      VRNF(NB,NZ)=0._r8
      VRNY(NB,NZ)=0._r8
      VRNZ(NB,NZ)=0._r8
      ATRP(NB,NZ)=0._r8
      FLG4(NB,NZ)=0._r8
      FDBK(NB,NZ)=1.0_r8
      FDBKX(NB,NZ)=1.0_r8
      IFLGA(NB,NZ)=0
      IFLGE(NB,NZ)=1
      IFLGF(NB,NZ)=0
      IFLGR(NB,NZ)=0
      IFLGQ(NB,NZ)=0
      NBTB(NB,NZ)=0
      D8850: DO M=1,pltpar%jpstgs
        IDAY(M,NB,NZ)=0
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
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IWTYP=phenology type:0=evergreen,1=cold decid,2=drought decid,3=1+2
!     IBTYP=turnover:0=all abve-grd,1=all leaf+petiole,2=none,3=between 1,2
!     IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!
      D6405: DO M=1,jsken
        DO NE=1,npelms        
          ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ) &
            +CFOPE(NE,instruct,M,NZ)*EPOLNB(NE,NB,NZ) &
            +CFOPE(NE,ifoliar,M,NZ)*(WTLFBE(NE,NB,NZ)*FWODLE(NE,k_fine_litr) &
            +WTNDBE(NE,NB,NZ)) &
            +CFOPE(NE,infoliar,M,NZ)*(WTSHEBE(NE,NB,NZ)*FWODBE(NE,k_fine_litr) &
            +WTHSKBE(NE,NB,NZ)+WTEARBE(NE,NB,NZ))

          ESNC(NE,M,k_woody_litr,0,NZ)=ESNC(NE,M,k_woody_litr,0,NZ) &
            +CFOPE(NE,icwood,M,NZ)*(WTLFBE(NE,NB,NZ)*FWODBE(NE,k_woody_litr) &
            +WTSHEBE(NE,NB,NZ)*FWODBE(NE,k_woody_litr))

          IF(ISTYP(NZ).EQ.iplt_annual.AND.IWTYP(NZ).NE.0)THEN
            WTRVE(NE,NZ)=WTRVE(NE,NZ)+CFOPE(NE,infoliar,M,NZ)*WTGRBE(NE,NB,NZ)
          ELSE
            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,infoliar,M,NZ)*WTGRBE(NE,NB,NZ)
          ENDIF
          IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,istalk,M,NZ)*WTSTKBE(NE,NB,NZ)
          ELSE
            WTSTDE(NE,M,NZ)=WTSTDE(NE,M,NZ)+CFOPE(NE,icwood,M,NZ)*WTSTKBE(NE,NB,NZ)
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
      WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)+CPOOLK(NB,NZ)
      DO NE=1,npelms
        WTRVE(NE,NZ)=WTRVE(NE,NZ)+EPOOL(NE,NB,NZ)
        IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
          D6406: DO M=1,jsken
            ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ) &
              +CFOPE(NE,instruct,M,NZ)*WTRSVBE(NE,NB,NZ)
          ENDDO D6406
        ELSE
          WTRVE(NE,NZ)=WTRVE(NE,NZ)+WTRSVBE(NE,NB,NZ)
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
  real(r8),INTENT(OUT) :: CPOOLK(JC1,JP1)
  integer :: L,NR,N,NE,NB
!     begin_execution
  associate(                           &
    EPOOL  => plt_biom%EPOOL         , &
    EPOLNB => plt_biom%EPOLNB        , &
    WTSHTBE=> plt_biom%WTSHTBE       , &
    WTLFBE => plt_biom%WTLFBE        , &
    WTRSVBE=> plt_biom%WTRSVBE       , &
    WTSHEBE=> plt_biom%WTSHEBE       , &
    WTHSKBE=> plt_biom%WTHSKBE       , &
    WTEARBE=> plt_biom%WTEARBE       , &
    WTGRBE => plt_biom%WTGRBE        , &
    CanPBStalkC => plt_biom%CanPBStalkC        , &
    EPOOLR => plt_biom%EPOOLR        , &
    WTSTXBE=> plt_biom%WTSTXBE       , &
    CanPBLeafShethC  => plt_biom%CanPBLeafShethC         , &
    WTRVE  => plt_biom%WTRVE         , &
    WTRT1E => plt_biom%WTRT1E        , &
    WTRT2E => plt_biom%WTRT2E        , &
    WTNDBE => plt_biom%WTNDBE        , &
    WTSTKBE=> plt_biom%WTSTKBE       , &
    RTWT1E => plt_biom%RTWT1E        , &
    NJ     => plt_site%NJ            , &
    NU     => plt_site%NU            , &
    IDTH   => plt_pheno%IDTH         , &
    SecndRootLen  => plt_morph%SecndRootLen        , &
    RTN2   => plt_morph%RTN2         , &
    PrimRootLen  => plt_morph%PrimRootLen        , &
    MY     => plt_morph%MY           , &
    NBR    => plt_morph%NBR          , &
    NRT    => plt_morph%NRT            &
  )
!     RESET BRANCH STATE VARIABLES
!
  DO NE=1,npelms
    DO NB=1,NBR(NZ)
      EPOOL(NE,NB,NZ)=0._r8
      EPOLNB(NE,NB,NZ)=0._r8
      WTSHTBE(NE,NB,NZ)=0._r8
      WTLFBE(NE,NB,NZ)=0._r8
      WTSHEBE(NE,NB,NZ)=0._r8
      WTSTKBE(NE,NB,NZ)=0._r8
      WTRSVBE(NE,NB,NZ)=0._r8
    ENDDO
  ENDDO
  D8835: DO NB=1,NBR(NZ)
    CPOOLK(NB,NZ)=0._r8
    CanPBStalkC(NB,NZ)=0._r8
    WTNDBE(1:npelms,NB,NZ)=0._r8
    WTHSKBE(1:npelms,NB,NZ)=0._r8
    WTEARBE(1:npelms,NB,NZ)=0._r8
    WTGRBE(1:npelms,NB,NZ)=0._r8
    CanPBLeafShethC(NB,NZ)=0._r8
    WTSTXBE(1:npelms,NB,NZ)=0._r8
  ENDDO D8835
!
!     RESET ROOT STATE VARIABLES
!
  D6416: DO L=NU,NJ
    DO  N=1,MY(NZ)
      EPOOLR(1:npelms,N,L,NZ)=0._r8
      DO  NR=1,NRT(NZ)
        WTRT1E(1:npelms,N,L,NR,NZ)=0._r8
        WTRT2E(1:npelms,N,L,NR,NZ)=0._r8
        RTWT1E(1:npelms,N,NR,NZ)=0._r8
        PrimRootLen(N,L,NR,NZ)=0._r8
        SecndRootLen(N,L,NR,NZ)=0._r8
        RTN2(N,L,NR,NZ)=0._r8
      enddo
    enddo
  ENDDO D6416
  WTRVE(1:npelms,NZ)=0._r8
  IDTH(NZ)=1
  end associate
  end subroutine ResetBranchRootStates

!------------------------------------------------------------------------------------------

  subroutine ResetDeadRootStates(NB,NZ,CPOOLK)
!     RESET STATE VARIABLES FROM DEAD BRANCHES
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8),intent(inout) :: CPOOLK(JC1,JP1)
  integer :: L,K,N
!     begin_execution
  associate(                          &
    EPOOL    => plt_biom%EPOOL      , &
    EPOLNB   => plt_biom%EPOLNB     , &
    WTSHTBE  => plt_biom%WTSHTBE    , &
    WTLFBE   => plt_biom%WTLFBE     , &
    WTNDBE   => plt_biom%WTNDBE     , &
    WTSTKBE  => plt_biom%WTSTKBE    , &
    WGLFE    => plt_biom%WGLFE      , &
    WTRSVBE  => plt_biom%WTRSVBE    , &
    WTSTXBE  => plt_biom%WTSTXBE    , &
    CanPBStalkC   => plt_biom%CanPBStalkC     , &
    WTGRBE   => plt_biom%WTGRBE     , &
    CanPBLeafShethC    => plt_biom%CanPBLeafShethC      , &
    WTHSKBE  => plt_biom%WTHSKBE    , &
    WTEARBE  => plt_biom%WTEARBE    , &
    WSSHE    => plt_biom%WSSHE      , &
    WSLF     => plt_biom%WSLF       , &
    WGNODE   => plt_biom%WGNODE     , &
    WGSHE    => plt_biom%WGSHE      , &
    WGLFLE   => plt_biom%WGLFLE     , &
    WGLFV    => plt_biom%WGLFV      , &
    WTSHEBE  => plt_biom%WTSHEBE    , &
    CPOOL3   => plt_photo%CPOOL3    , &
    CPOOL4   => plt_photo%CPOOL4    , &
    CO2B     => plt_photo%CO2B      , &
    HCOB     => plt_photo%HCOB      , &
    GRWTB    => plt_allom%GRWTB     , &
    GRNXB    => plt_morph%GRNXB     , &
    GRNOB    => plt_morph%GRNOB     , &
    CanPBLA    => plt_morph%CanPBLA     , &
    CanPLBSA    => plt_morph%CanPLBSA     , &
    NBR      => plt_morph%NBR       , &
    ARLF1    => plt_morph%ARLF1     , &
    HTSHE    => plt_morph%HTSHE     , &
    HTNODX   => plt_morph%HTNODX    , &
    SURF     => plt_morph%SURF      , &
    HTNODE   => plt_morph%HTNODE    , &
    ARLFV    => plt_morph%ARLFV     , &
    SURFB    => plt_morph%SURFB     , &
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
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     GRNOB=seed set number
!     GRNXB=potential number of seed set sites
!     GRWTB=individual seed size
!     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
!     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
!     WSLF=leaf protein mass
!     CanPBLA=branch leaf area
!     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
!     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
!
  CPOOLK(NB,NZ)=0._r8
  EPOOL(1:npelms,NB,NZ)=0._r8
  EPOLNB(1:npelms,NB,NZ)=0._r8
  WTSHTBE(1:npelms,NB,NZ)=0._r8
  WTLFBE(1:npelms,NB,NZ)=0._r8
  WTNDBE(1:npelms,NB,NZ)=0._r8
  WTSHEBE(1:npelms,NB,NZ)=0._r8
  WTSTKBE(1:npelms,NB,NZ)=0._r8
  WTRSVBE(1:npelms,NB,NZ)=0._r8
  WTHSKBE(1:npelms,NB,NZ)=0._r8
  WTEARBE(1:npelms,NB,NZ)=0._r8
  WTGRBE(1:npelms,NB,NZ)=0._r8
  CanPBStalkC(NB,NZ)=0._r8
  CanPBLeafShethC(NB,NZ)=0._r8
  GRNXB(NB,NZ)=0._r8
  GRNOB(NB,NZ)=0._r8
  GRWTB(NB,NZ)=0._r8
  CanPBLA(NB,NZ)=0._r8
  WTSTXBE(1:npelms,NB,NZ)=0._r8

  D8855: DO K=0,JNODS1
    IF(K.NE.0)THEN
      CPOOL3(K,NB,NZ)=0._r8
      CPOOL4(K,NB,NZ)=0._r8
      CO2B(K,NB,NZ)=0._r8
      HCOB(K,NB,NZ)=0._r8
    ENDIF
    ARLF1(K,NB,NZ)=0._r8
    HTNODE(K,NB,NZ)=0._r8
    HTNODX(K,NB,NZ)=0._r8
    HTSHE(K,NB,NZ)=0._r8
    WSLF(K,NB,NZ)=0._r8
    WSSHE(K,NB,NZ)=0._r8
    WGLFE(1:npelms,K,NB,NZ)=0._r8
    WGSHE(1:npelms,K,NB,NZ)=0._r8
    WGNODE(1:npelms,K,NB,NZ)=0._r8
    D8865: DO L=1,JC1
      ARLFV(L,NZ)=ARLFV(L,NZ)-CanPLNBLA(L,K,NB,NZ)
      WGLFV(L,NZ)=WGLFV(L,NZ)-WGLFLE(ielmc,L,K,NB,NZ)
      CanPLNBLA(L,K,NB,NZ)=0._r8
      WGLFLE(1:npelms,L,K,NB,NZ)=0._r8
      IF(K.NE.0)THEN
        D8860: DO N=1,JLI1
          SURF(N,L,K,NB,NZ)=0._r8
        ENDDO D8860
      ENDIF
    ENDDO D8865
  ENDDO D8855
  D8875: DO L=1,JC1
    CanPLBSA(L,NB,NZ)=0._r8
    DO  N=1,JLI1
      SURFB(N,L,NB,NZ)=0._r8
    enddo
  ENDDO D8875
  end associate
  end subroutine ResetDeadRootStates


end module LitterFallMod
