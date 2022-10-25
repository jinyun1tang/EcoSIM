module LitterFallMod

  use data_kind_mod, only : r8 => SHR_KIND_R8
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
    VOLWP      =>   plt_ew%VOLWP       , &
    WTRTE      =>   plt_biom%WTRTE     , &
    WTRVE      =>   plt_biom%WTRVE     , &
    ISTYP      =>   plt_pheno%ISTYP    , &
    IDTHR      =>   plt_pheno%IDTHR    , &
    IDTHP      =>   plt_pheno%IDTHP    , &
    IFLGI      =>   plt_pheno%IFLGI    , &
    WSTR       =>   plt_pheno%WSTR     , &
    IDAY       =>   plt_pheno%IDAY     , &
    PP         =>   plt_site%PP        , &
    IYRC       =>   plt_site%IYRC      , &
    VOLWOU     =>   plt_site%VOLWOU    , &
    ZNOON      =>   plt_site%ZNOON     , &
    HTCTL      =>   plt_morph%HTCTL    , &
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
    .AND.(ISTYP(NZ).NE.0.OR.(I.GE.IDAYH(NZ) &
    .AND.IYRC.GE.IYRH(NZ))))THEN
    IDTHY=0
!
!     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
!
    call LiterfallFromDeadBranches(I,J,NZ,IDTHY,CPOOLK)

    IF(IDTHY.EQ.NBR(NZ))THEN
      IDTHP(NZ)=1
      NBT(NZ)=0
      WSTR(NZ)=0._r8
      IF(IFLGI(NZ).EQ.1)THEN
        NBR(NZ)=1
      ELSE
        NBR(NZ)=0
      ENDIF
      HTCTL(NZ)=0._r8
      VOLWOU=VOLWOU+VOLWP(NZ)
      UVOLO=UVOLO+VOLWP(NZ)
      VOLWP(NZ)=0._r8
!
!     RESET LIVING FLAGS
!
!     WTRVC,WTRT=PFT storage,root C
!     ISTYP=growth habit:0=annual,1=perennial
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     PP=PFT population
!     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!
      IF(WTRVE(ielmc,NZ).LT.1.0E-04*WTRTE(ielmc,NZ) &
        .AND.ISTYP(NZ).NE.0)IDTHR(NZ)=1
      IF(ISTYP(NZ).EQ.0)IDTHR(NZ)=1
      IF(JHVST(NZ).NE.0)IDTHR(NZ)=1
      IF(PP(NZ).LE.0.0)IDTHR(NZ)=1
      IF(IDTHR(NZ).EQ.1)IDTHP(NZ)=1
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
  implicit none
  integer, intent(in) :: I,J,NZ
  REAL(R8),INTENT(INOUT) :: CPOOLK(JC1,JP1)
  integer :: L,M,NR,NB,N
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
    FWODLN     =>   plt_allom%FWODLN  , &
    FWODLP     =>   plt_allom%FWODLP  , &
    FWODB      =>   plt_allom%FWODB   , &
    FWODSP     =>   plt_allom%FWODSP  , &
    FWODSN     =>   plt_allom%FWODSN  , &
    FWOOD      =>   plt_allom%FWOOD   , &
    FWODR      =>   plt_allom%FWODR   , &
    FWODRN     =>   plt_allom%FWODRN  , &
    FWODRP     =>   plt_allom%FWODRP  , &
    FWOODN     =>   plt_allom%FWOODN  , &
    FWOODP     =>   plt_allom%FWOODP  , &
    IFLGI      =>   plt_pheno%IFLGI   , &
    ISTYP      =>   plt_pheno%ISTYP   , &
    IWTYP      =>   plt_pheno%IWTYP   , &
    IDTHR      =>   plt_pheno%IDTHR   , &
    IGTYP      =>   plt_pheno%IGTYP   , &
    IBTYP      =>   plt_pheno%IBTYP   , &
    IDTHP      =>   plt_pheno%IDTHP   , &
    ESNC       =>   plt_bgcr%ESNC     , &
    LYRC       =>   plt_site%LYRC     , &
    NJ         =>   plt_site%NJ       , &
    NU         =>   plt_site%NU       , &
    IDATA      =>   plt_site%IDATA    , &
    CFOPP      =>   plt_soilchem%CFOPP, &
    CFOPN      =>   plt_soilchem%CFOPN, &
    CFOPC      =>   plt_soilchem%CFOPC, &
    MY         =>   plt_morph%MY      , &
    NG         =>   plt_morph%NG      , &
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
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ) &
            +CFOPC(0,M,NZ)*(EPOOL(NB,ielmc,NZ)+EPOLNB(NB,ielmc,NZ) &
            +CPOOLK(NB,NZ)+WTRSVBE(NB,ielmc,NZ)) &
            +CFOPC(1,M,NZ)*(WTLFBE(NB,ielmc,NZ)*FWODB(1) &
            +WTNDBE(NB,ielmc,NZ)) &
            +CFOPC(2,M,NZ)*(WTSHEBE(NB,ielmc,NZ)*FWODB(1) &
            +WTHSKBE(NB,ielmc,NZ)+WTEARBE(NB,ielmc,NZ))
          ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ) &
            +CFOPC(5,M,NZ)*(WTLFBE(NB,ielmc,NZ)*FWODB(0) &
            +WTSHEBE(NB,ielmc,NZ)*FWODB(0))
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ) &
            +CFOPN(0,M,NZ)*(EPOOL(NB,ielmn,NZ)+EPOLNB(NB,ielmn,NZ) &
            +WTRSVBE(NB,ielmn,NZ)) &
            +CFOPN(1,M,NZ)*(WTLFBE(NB,ielmn,NZ)*FWODLN(1) &
            +WTNDBE(NB,ielmn,NZ)) &
            +CFOPN(2,M,NZ)*(WTSHEBE(NB,ielmn,NZ)*FWODSN(1) &
            +WTHSKBE(NB,ielmn,NZ)+WTEARBE(NB,ielmn,NZ))
          ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ) &
            +CFOPN(5,M,NZ)*(WTLFBE(NB,ielmn,NZ)*FWODLN(0) &
            +WTSHEBE(NB,ielmn,NZ)*FWODSN(0))
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ) &
            +CFOPP(0,M,NZ)*(EPOOL(NB,ielmp,NZ)+EPOLNB(NB,ielmp,NZ) &
            +WTRSVBE(NB,ielmp,NZ)) &
            +CFOPP(1,M,NZ)*(WTLFBE(NB,ielmp,NZ)*FWODLP(1) &
            +WTNDBE(NB,ielmp,NZ)) &
            +CFOPP(2,M,NZ)*(WTSHEBE(NB,ielmp,NZ)*FWODSP(1) &
            +WTHSKBE(NB,ielmp,NZ)+WTEARBE(NB,ielmp,NZ))
          ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ) &
            +CFOPP(5,M,NZ)*(WTLFBE(NB,ielmp,NZ)*FWODLP(0) &
            +WTSHEBE(NB,ielmp,NZ)*FWODSP(0))
          IF(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).NE.0)THEN
            WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)+CFOPC(2,M,NZ)*WTGRBE(NB,ielmc,NZ)
            WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)+CFOPN(2,M,NZ)*WTGRBE(NB,ielmn,NZ)
            WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)+CFOPP(2,M,NZ)*WTGRBE(NB,ielmp,NZ)
          ELSE
            ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(2,M,NZ)*WTGRBE(NB,ielmc,NZ)
            ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(2,M,NZ)*WTGRBE(NB,ielmn,NZ)
            ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(2,M,NZ)*WTGRBE(NB,ielmp,NZ)
          ENDIF
          IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
            ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(3,M,NZ)*WTSTKBE(NB,ielmc,NZ)
            ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(3,M,NZ)*WTSTKBE(NB,ielmn,NZ)
            ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(3,M,NZ)*WTSTKBE(NB,ielmp,NZ)
          ELSE
            WTSTDE(M,ielmc,NZ)=WTSTDE(M,ielmc,NZ)+CFOPC(5,M,NZ)*WTSTKBE(NB,ielmc,NZ)
            WTSTDE(M,ielmn,NZ)=WTSTDE(M,ielmn,NZ)+CFOPN(5,M,NZ)*WTSTKBE(NB,ielmn,NZ)
            WTSTDE(M,ielmp,NZ)=WTSTDE(M,ielmp,NZ)+CFOPP(5,M,NZ)*WTSTKBE(NB,ielmp,NZ)
          ENDIF
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
        D6415: DO L=NU,NJ
          DO N=1,MY(NZ)
            ESNC(M,ielmc,1,L,NZ)=ESNC(M,ielmc,1,L,NZ)+CFOPC(0,M,NZ)*EPOOLR(ielmc,N,L,NZ)
            ESNC(M,ielmn,1,L,NZ)=ESNC(M,ielmn,1,L,NZ)+CFOPN(0,M,NZ)*EPOOLR(ielmn,N,L,NZ)
            ESNC(M,ielmp,1,L,NZ)=ESNC(M,ielmp,1,L,NZ)+CFOPP(0,M,NZ)*EPOOLR(ielmp,N,L,NZ)
            DO NR=1,NRT(NZ)
              ESNC(M,ielmc,0,L,NZ)=ESNC(M,ielmc,0,L,NZ)+CFOPC(5,M,NZ) &
                *(WTRT1E(ielmc,N,L,NR,NZ)+WTRT2E(ielmc,N,L,NR,NZ))*FWODR(0)
              ESNC(M,ielmn,0,L,NZ)=ESNC(M,ielmn,0,L,NZ)+CFOPN(5,M,NZ) &
                *(WTRT1E(ielmn,N,L,NR,NZ)+WTRT2E(ielmn,N,L,NR,NZ))*FWODRN(0)
              ESNC(M,ielmp,0,L,NZ)=ESNC(M,ielmp,0,L,NZ)+CFOPP(5,M,NZ) &
                *(WTRT1E(ielmp,N,L,NR,NZ)+WTRT2E(ielmp,N,L,NR,NZ))*FWODRP(0)
              ESNC(M,ielmc,1,L,NZ)=ESNC(M,ielmc,1,L,NZ)+CFOPC(4,M,NZ) &
                *(WTRT1E(ielmc,N,L,NR,NZ)+WTRT2E(ielmc,N,L,NR,NZ))*FWODR(1)
              ESNC(M,ielmn,1,L,NZ)=ESNC(M,ielmn,1,L,NZ)+CFOPN(4,M,NZ) &
                *(WTRT1E(ielmn,N,L,NR,NZ)+WTRT2E(ielmn,N,L,NR,NZ))*FWODRN(1)
              ESNC(M,ielmp,1,L,NZ)=ESNC(M,ielmp,1,L,NZ)+CFOPP(4,M,NZ) &
                *(WTRT1E(ielmp,N,L,NR,NZ)+WTRT2E(ielmp,N,L,NR,NZ))*FWODRP(1)
            ENDDO
          ENDDO
        ENDDO D6415
        ESNC(M,ielmc,0,NG(NZ),NZ)=ESNC(M,ielmc,0,NG(NZ),NZ) &
          +CFOPC(0,M,NZ)*WTRVE(ielmc,NZ)*FWOOD(0)
        ESNC(M,ielmn,0,NG(NZ),NZ)=ESNC(M,ielmn,0,NG(NZ),NZ) &
          +CFOPN(0,M,NZ)*WTRVE(ielmn,NZ)*FWOODN(0)
        ESNC(M,ielmp,0,NG(NZ),NZ)=ESNC(M,ielmp,0,NG(NZ),NZ) &
          +CFOPP(0,M,NZ)*WTRVE(ielmp,NZ)*FWOODP(0)
        ESNC(M,ielmc,1,NG(NZ),NZ)=ESNC(M,ielmc,1,NG(NZ),NZ) &
          +CFOPC(0,M,NZ)*WTRVE(ielmc,NZ)*FWOOD(1)
        ESNC(M,ielmn,1,NG(NZ),NZ)=ESNC(M,ielmn,1,NG(NZ),NZ) &
          +CFOPN(0,M,NZ)*WTRVE(ielmn,NZ)*FWOODN(1)
        ESNC(M,ielmp,1,NG(NZ),NZ)=ESNC(M,ielmp,1,NG(NZ),NZ) &
          +CFOPP(0,M,NZ)*WTRVE(ielmp,NZ)*FWOODP(1)
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
    IF(ISTYP(NZ).NE.0.AND.JHVST(NZ).EQ.0)THEN
      IF(I.LT.LYRC)THEN
        IDAY0(NZ)=I+1
        IYR0(NZ)=IDATA(3)
      ELSE
        IDAY0(NZ)=1
        IYR0(NZ)=IDATA(3)+1
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine LiterfallFromRootShootStorage
!------------------------------------------------------------------------------------------

  subroutine LiterfallFromDeadRoots(I,J,NZ)
  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: L,M,NR,N
!     begin_execution
  associate(                                &
    WTRT2E    =>   plt_biom%WTRT2E    , &
    RTWT1     =>   plt_biom%RTWT1     , &
    RTWT1N    =>   plt_biom%RTWT1N    , &
    RTWT1P    =>   plt_biom%RTWT1P    , &
    WTRT1E    =>   plt_biom%WTRT1E    , &
    EPOOLR    =>   plt_biom%EPOOLR    , &
    WTRTD     =>   plt_biom%WTRTD     , &
    WTRTL     =>   plt_biom%WTRTL     , &
    WSRTL     =>   plt_biom%WSRTL     , &
    EPOOLN    =>   plt_biom%EPOOLN    , &
    WTNDLE    =>   plt_biom%WTNDLE    , &
    FWODR     =>   plt_allom%FWODR    , &
    FWODRN    =>   plt_allom%FWODRN   , &
    FWODRP    =>   plt_allom%FWODRP   , &
    IDTHR     =>   plt_pheno%IDTHR    , &
    CFOPC     =>   plt_soilchem%CFOPC , &
    CFOPN     =>   plt_soilchem%CFOPN , &
    CFOPP     =>   plt_soilchem%CFOPP , &
    CO2P      =>   plt_rbgc%CO2P  , &
    OXYP      =>   plt_rbgc%OXYP  , &
    CO2A      =>   plt_rbgc%CO2A  , &
    OXYA      =>   plt_rbgc%OXYA  , &
    CH4P      =>   plt_rbgc%CH4P  , &
    CH4A      =>   plt_rbgc%CH4A  , &
    H2GP      =>   plt_rbgc%H2GP  , &
    H2GA      =>   plt_rbgc%H2GA  , &
    Z2OP      =>   plt_rbgc%Z2OP      , &
    Z2OA      =>   plt_rbgc%Z2OA      , &
    ZH3P      =>   plt_rbgc%ZH3P      , &
    ZH3A      =>   plt_rbgc%ZH3A      , &
    NJ        =>   plt_site%NJ        , &
    NU        =>   plt_site%NU        , &
    RH2GZ     =>   plt_bgcr%RH2GZ     , &
    RNH3Z     =>   plt_bgcr%RNH3Z     , &
    RN2OZ     =>   plt_bgcr%RN2OZ     , &
    RCH4Z     =>   plt_bgcr%RCH4Z     , &
    ROXYZ     =>   plt_bgcr%ROXYZ     , &
    RCO2Z     =>   plt_bgcr%RCO2Z     , &
    ESNC      =>   plt_bgcr%ESNC      , &
    RTDNP     =>   plt_morph%RTDNP    , &
    RTVLW     =>   plt_morph%RTVLW    , &
    RTARP     =>   plt_morph%RTARP    , &
    RTVLP     =>   plt_morph%RTVLP    , &
    RTNL      =>   plt_morph%RTNL     , &
    RTN1      =>   plt_morph%RTN1     , &
    RTDP1     =>   plt_morph%RTDP1    , &
    RTLGA     =>   plt_morph%RTLGA    , &
    RRAD1     =>   plt_morph%RRAD1    , &
    RRAD2     =>   plt_morph%RRAD2    , &
    RRAD1M    =>   plt_morph%RRAD1M   , &
    RRAD2M    =>   plt_morph%RRAD2M   , &
    RTLGP     =>   plt_morph%RTLGP    , &
    RTLG1     =>   plt_morph%RTLG1    , &
    RTLG2     =>   plt_morph%RTLG2    , &
    INTYP     =>   plt_morph%INTYP    , &
    RTN2      =>   plt_morph%RTN2     , &
    MY        =>   plt_morph%MY       , &
    NIX       =>   plt_morph%NIX      , &
    NINR      =>   plt_morph%NINR     , &
    SDPTH     =>   plt_morph%SDPTH    , &
    NG        =>   plt_morph%NG       , &
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
    DO 8900 N=1,MY(NZ)
      DO 8895 L=NU,NJ
        D6410: DO M=1,jsken
          ESNC(M,ielmc,1,L,NZ)=ESNC(M,ielmc,1,L,NZ)+CFOPC(0,M,NZ)*EPOOLR(ielmc,N,L,NZ)
          ESNC(M,ielmn,1,L,NZ)=ESNC(M,ielmn,1,L,NZ)+CFOPN(0,M,NZ)*EPOOLR(ielmn,N,L,NZ)
          ESNC(M,ielmp,1,L,NZ)=ESNC(M,ielmp,1,L,NZ)+CFOPP(0,M,NZ)*EPOOLR(ielmp,N,L,NZ)
          DO  NR=1,NRT(NZ)
            ESNC(M,ielmc,0,L,NZ)=ESNC(M,ielmc,0,L,NZ)+CFOPC(5,M,NZ) &
              *(WTRT1E(ielmc,N,L,NR,NZ)+WTRT2E(ielmc,N,L,NR,NZ))*FWODR(0)
            ESNC(M,ielmn,0,L,NZ)=ESNC(M,ielmn,0,L,NZ)+CFOPN(5,M,NZ) &
              *(WTRT1E(ielmn,N,L,NR,NZ)+WTRT2E(ielmn,N,L,NR,NZ))*FWODRN(0)
            ESNC(M,ielmp,0,L,NZ)=ESNC(M,ielmp,0,L,NZ)+CFOPP(5,M,NZ) &
              *(WTRT1E(ielmp,N,L,NR,NZ)+WTRT2E(ielmp,N,L,NR,NZ))*FWODRP(0)
            ESNC(M,ielmc,1,L,NZ)=ESNC(M,ielmc,1,L,NZ)+CFOPC(4,M,NZ) &
              *(WTRT1E(ielmc,N,L,NR,NZ)+WTRT2E(ielmc,N,L,NR,NZ))*FWODR(1)
            ESNC(M,ielmn,1,L,NZ)=ESNC(M,ielmn,1,L,NZ)+CFOPN(4,M,NZ) &
              *(WTRT1E(ielmn,N,L,NR,NZ)+WTRT2E(ielmn,N,L,NR,NZ))*FWODRN(1)
            ESNC(M,ielmp,1,L,NZ)=ESNC(M,ielmp,1,L,NZ)+CFOPP(4,M,NZ) &
              *(WTRT1E(ielmp,N,L,NR,NZ)+WTRT2E(ielmp,N,L,NR,NZ))*FWODRP(1)
          enddo
        ENDDO D6410
!
!     RELEASE GAS CONTENTS OF DEAD ROOTS
!
        RCO2Z(NZ)=RCO2Z(NZ)-CO2A(N,L,NZ)-CO2P(N,L,NZ)
        ROXYZ(NZ)=ROXYZ(NZ)-OXYA(N,L,NZ)-OXYP(N,L,NZ)
        RCH4Z(NZ)=RCH4Z(NZ)-CH4A(N,L,NZ)-CH4P(N,L,NZ)
        RN2OZ(NZ)=RN2OZ(NZ)-Z2OA(N,L,NZ)-Z2OP(N,L,NZ)
        RNH3Z(NZ)=RNH3Z(NZ)-ZH3A(N,L,NZ)-ZH3P(N,L,NZ)
        RH2GZ(NZ)=RH2GZ(NZ)-H2GA(N,L,NZ)-H2GP(N,L,NZ)
        CO2A(N,L,NZ)=0._r8
        OXYA(N,L,NZ)=0._r8
        CH4A(N,L,NZ)=0._r8
        Z2OA(N,L,NZ)=0._r8
        ZH3A(N,L,NZ)=0._r8
        H2GA(N,L,NZ)=0._r8
        CO2P(N,L,NZ)=0._r8
        OXYP(N,L,NZ)=0._r8
        CH4P(N,L,NZ)=0._r8
        Z2OP(N,L,NZ)=0._r8
        ZH3P(N,L,NZ)=0._r8
        H2GP(N,L,NZ)=0._r8
!
!     RESET STATE VARIABLES OF DEAD ROOTS
!
!     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
!     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!     RTLG1,RTLG2=primary,secondary root length
!     RTN2=number of secondary root axes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!     WTRTL,WTRTD=active,actual root C mass
!     WSRTL=root protein C mass
!     RTN1,RTNL=number of primary,secondary root axes
!     RTDNP,RTLGP=root length density,root length per plant
!     RTVLW,RTVLP=root or myco aqueous,gaseous volume
!     RTARP=root surface area per plant
!     RTLGA=average secondary root length
!
        D8870: DO NR=1,NRT(NZ)
          WTRT1E(ielmc,N,L,NR,NZ)=0._r8
          WTRT1E(ielmn,N,L,NR,NZ)=0._r8
          WTRT1E(ielmp,N,L,NR,NZ)=0._r8
          WTRT2E(ielmc,N,L,NR,NZ)=0._r8
          WTRT2E(ielmn,N,L,NR,NZ)=0._r8
          WTRT2E(ielmp,N,L,NR,NZ)=0._r8
          RTWT1(N,NR,NZ)=0._r8
          RTWT1N(N,NR,NZ)=0._r8
          RTWT1P(N,NR,NZ)=0._r8
          RTLG1(N,L,NR,NZ)=0._r8
          RTLG2(N,L,NR,NZ)=0._r8
          RTN2(N,L,NR,NZ)=0._r8
        ENDDO D8870
        EPOOLR(:,N,L,NZ)=0._r8
        WTRTL(N,L,NZ)=0._r8
        WTRTD(N,L,NZ)=0._r8
        WSRTL(N,L,NZ)=0._r8
        RTN1(N,L,NZ)=0._r8
        RTNL(N,L,NZ)=0._r8
        RTLGP(N,L,NZ)=0._r8
        RTDNP(N,L,NZ)=0._r8
        RTVLP(N,L,NZ)=0._r8
        RTVLW(N,L,NZ)=0._r8
        RRAD1(N,L,NZ)=RRAD1M(N,NZ)
        RRAD2(N,L,NZ)=RRAD2M(N,NZ)
        RTARP(N,L,NZ)=0._r8
        RTLGA(N,L,NZ)=RTLGAX
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
            ESNC(M,ielmc,1,L,NZ)=ESNC(M,ielmc,1,L,NZ)+CFOPC(4,M,NZ) &
              *WTNDLE(L,ielmc,NZ)+CFOPC(0,M,NZ)*EPOOLN(L,ielmc,NZ)
            ESNC(M,ielmn,1,L,NZ)=ESNC(M,ielmn,1,L,NZ)+CFOPN(4,M,NZ) &
              *WTNDLE(L,ielmn,NZ)+CFOPN(0,M,NZ)*EPOOLN(L,ielmn,NZ)
            ESNC(M,ielmp,1,L,NZ)=ESNC(M,ielmp,1,L,NZ)+CFOPP(4,M,NZ) &
              *WTNDLE(L,ielmp,NZ)+CFOPP(0,M,NZ)*EPOOLN(L,ielmp,NZ)
          ENDDO D6420
          WTNDLE(L,ielmc,NZ)=0._r8
          WTNDLE(L,ielmn,NZ)=0._r8
          WTNDLE(L,ielmp,NZ)=0._r8
          EPOOLN(L,ielmc,NZ)=0._r8
          EPOOLN(L,ielmn,NZ)=0._r8
          EPOOLN(L,ielmp,NZ)=0._r8
        ENDIF
8895  CONTINUE
8900  CONTINUE
!
!     RESET DEPTH VARIABLES OF DEAD ROOTS
!
!     NINR=deepest root layer
!     RTDP1=primary root depth from soil surface
!     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
!
    DO 8795 NR=1,NRT(NZ)
      NINR(NR,NZ)=NG(NZ)
      DO 8790 N=1,MY(NZ)
        RTDP1(N,NR,NZ)=SDPTH(NZ)
        RTWT1(N,NR,NZ)=0._r8
        RTWT1N(N,NR,NZ)=0._r8
        RTWT1P(N,NR,NZ)=0._r8
8790  CONTINUE
8795  CONTINUE
    NIX(NZ)=NG(NZ)
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
  integer :: M,NB
!     begin_execution
  associate(                                &
    IHVST     =>  plt_distb%IHVST     , &
    FWODB     =>  plt_allom%FWODB     , &
    FWODSP    =>  plt_allom%FWODSP    , &
    FWODSN    =>  plt_allom%FWODSN    , &
    FWODLN    =>  plt_allom%FWODLN    , &
    FWODLP    =>  plt_allom%FWODLP    , &
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
    CFOPC     =>  plt_soilchem%CFOPC  , &
    CFOPN     =>  plt_soilchem%CFOPN  , &
    CFOPP     =>  plt_soilchem%CFOPP  , &
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
    FDBK      =>  plt_photo%FDBK      , &
    FDBKX     =>  plt_photo%FDBKX       &
  )
  DO 8845 NB=1,NBR(NZ)
    IF(IDTHB(NB,NZ).EQ.1)THEN
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
      DO 8850 M=1,10
        IDAY(M,NB,NZ)=0
8850  CONTINUE
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
        ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ) &
          +CFOPC(0,M,NZ)*EPOLNB(NB,ielmc,NZ) &
          +CFOPC(1,M,NZ)*(WTLFBE(NB,ielmc,NZ)*FWODB(1) &
          +WTNDBE(NB,ielmc,NZ)) &
          +CFOPC(2,M,NZ)*(WTSHEBE(NB,ielmc,NZ)*FWODB(1) &
          +WTHSKBE(NB,ielmc,NZ)+WTEARBE(NB,ielmc,NZ))
        ESNC(M,ielmc,0,0,NZ)=ESNC(M,ielmc,0,0,NZ) &
          +CFOPC(5,M,NZ)*(WTLFBE(NB,ielmc,NZ)*FWODB(0) &
          +WTSHEBE(NB,ielmc,NZ)*FWODB(0))
        ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ) &
          +CFOPN(0,M,NZ)*EPOLNB(NB,ielmn,NZ) &
          +CFOPN(1,M,NZ)*(WTLFBE(NB,ielmn,NZ)*FWODLN(1) &
          +WTNDBE(NB,ielmn,NZ)) &
          +CFOPN(2,M,NZ)*(WTSHEBE(NB,ielmn,NZ)*FWODSN(1) &
          +WTHSKBE(NB,ielmn,NZ)+WTEARBE(NB,ielmn,NZ))
        ESNC(M,ielmn,0,0,NZ)=ESNC(M,ielmn,0,0,NZ) &
          +CFOPN(5,M,NZ)*(WTLFBE(NB,ielmn,NZ)*FWODLN(0) &
          +WTSHEBE(NB,ielmn,NZ)*FWODSN(0))
        ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ) &
          +CFOPP(0,M,NZ)*EPOLNB(NB,ielmp,NZ) &
          +CFOPP(1,M,NZ)*(WTLFBE(NB,ielmp,NZ)*FWODLP(1) &
          +WTNDBE(NB,ielmp,NZ)) &
          +CFOPP(2,M,NZ)*(WTSHEBE(NB,ielmp,NZ)*FWODSP(1) &
          +WTHSKBE(NB,ielmp,NZ)+WTEARBE(NB,ielmp,NZ))
        ESNC(M,ielmp,0,0,NZ)=ESNC(M,ielmp,0,0,NZ) &
          +CFOPP(5,M,NZ)*(WTLFBE(NB,ielmp,NZ)*FWODLP(0) &
          +WTSHEBE(NB,ielmp,NZ)*FWODSP(0))
        IF(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).NE.0)THEN
          WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)+CFOPC(2,M,NZ)*WTGRBE(NB,ielmc,NZ)
          WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)+CFOPN(2,M,NZ)*WTGRBE(NB,ielmn,NZ)
          WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)+CFOPP(2,M,NZ)*WTGRBE(NB,ielmp,NZ)
        ELSE
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(2,M,NZ)*WTGRBE(NB,ielmc,NZ)
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(2,M,NZ)*WTGRBE(NB,ielmn,NZ)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(2,M,NZ)*WTGRBE(NB,ielmp,NZ)
        ENDIF
        IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ)+CFOPC(3,M,NZ)*WTSTKBE(NB,ielmc,NZ)
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ)+CFOPN(3,M,NZ)*WTSTKBE(NB,ielmn,NZ)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ)+CFOPP(3,M,NZ)*WTSTKBE(NB,ielmp,NZ)
        ELSE
          WTSTDE(M,ielmc,NZ)=WTSTDE(M,ielmc,NZ)+CFOPC(5,M,NZ)*WTSTKBE(NB,ielmc,NZ)
          WTSTDE(M,ielmn,NZ)=WTSTDE(M,ielmn,NZ)+CFOPN(5,M,NZ)*WTSTKBE(NB,ielmn,NZ)
          WTSTDE(M,ielmp,NZ)=WTSTDE(M,ielmp,NZ)+CFOPP(5,M,NZ)*WTSTKBE(NB,ielmp,NZ)
        ENDIF
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
      WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)+EPOOL(NB,ielmc,NZ)+CPOOLK(NB,NZ)
      WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)+EPOOL(NB,ielmn,NZ)
      WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)+EPOOL(NB,ielmp,NZ)
      IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
        D6406: DO M=1,jsken
          ESNC(M,ielmc,1,0,NZ)=ESNC(M,ielmc,1,0,NZ) &
            +CFOPC(0,M,NZ)*WTRSVBE(NB,ielmc,NZ)
          ESNC(M,ielmn,1,0,NZ)=ESNC(M,ielmn,1,0,NZ) &
            +CFOPN(0,M,NZ)*WTRSVBE(NB,ielmn,NZ)
          ESNC(M,ielmp,1,0,NZ)=ESNC(M,ielmp,1,0,NZ) &
            +CFOPP(0,M,NZ)*WTRSVBE(NB,ielmp,NZ)
        ENDDO D6406
      ELSE
        WTRVE(ielmc,NZ)=WTRVE(ielmc,NZ)+WTRSVBE(NB,ielmc,NZ)
        WTRVE(ielmn,NZ)=WTRVE(ielmn,NZ)+WTRSVBE(NB,ielmn,NZ)
        WTRVE(ielmp,NZ)=WTRVE(ielmp,NZ)+WTRSVBE(NB,ielmp,NZ)
      ENDIF
!
      call ResetDeadRootStates(NB,NZ,CPOOLK)

      IDTHY=IDTHY+1
    ENDIF
8845  CONTINUE
  end associate
  end subroutine LiterfallFromDeadBranches

!------------------------------------------------------------------------------------------

  subroutine ResetBranchRootStates(NZ,CPOOLK)
  implicit none
  integer, intent(in) :: NZ
  real(r8),INTENT(OUT) :: CPOOLK(JC1,JP1)
  integer :: L,NR,N,NB,NE
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
    WVSTKB => plt_biom%WVSTKB        , &
    EPOOLR => plt_biom%EPOOLR        , &
    WTSTXB => plt_biom%WTSTXB        , &
    WTSTXN => plt_biom%WTSTXN        , &
    WTSTXP => plt_biom%WTSTXP        , &
    WTLSB  => plt_biom%WTLSB         , &
    WTRVE  => plt_biom%WTRVE         , &
    WTRT1E => plt_biom%WTRT1E        , &
    WTRT2E => plt_biom%WTRT2E        , &
    WTNDBE => plt_biom%WTNDBE        , &
    WTSTKBE=> plt_biom%WTSTKBE       , &
    RTWT1  => plt_biom%RTWT1         , &
    RTWT1N => plt_biom%RTWT1N        , &
    RTWT1P => plt_biom%RTWT1P        , &
    NJ     => plt_site%NJ            , &
    NU     => plt_site%NU            , &
    IDTH   => plt_pheno%IDTH         , &
    RTLG2  => plt_morph%RTLG2        , &
    RTN2   => plt_morph%RTN2         , &
    RTLG1  => plt_morph%RTLG1        , &
    MY     => plt_morph%MY           , &
    NBR    => plt_morph%NBR          , &
    NRT    => plt_morph%NRT            &
  )
!     RESET BRANCH STATE VARIABLES
!
  DO NE=1,npelms
    DO NB=1,NBR(NZ)
      EPOOL(NB,NE,NZ)=0._r8
      EPOLNB(NB,NE,NZ)=0._r8
      WTSHTBE(NB,NE,NZ)=0._r8
      WTLFBE(NB,NE,NZ)=0._r8
      WTSHEBE(NB,NE,NZ)=0._r8
      WTSTKBE(NB,NE,NZ)=0._r8
      WTRSVBE(NB,NE,NZ)=0._r8
    ENDDO
  ENDDO
  D8835: DO NB=1,NBR(NZ)

    CPOOLK(NB,NZ)=0._r8
    WTNDBE(NB,ielmc,NZ)=0._r8
    WVSTKB(NB,NZ)=0._r8
    WTHSKBE(NB,ielmc,NZ)=0._r8
    WTEARBE(NB,ielmc,NZ)=0._r8
    WTGRBE(NB,ielmc,NZ)=0._r8
    WTLSB(NB,NZ)=0._r8
    WTNDBE(NB,ielmn,NZ)=0._r8
    WTHSKBE(NB,ielmn,NZ)=0._r8
    WTEARBE(NB,ielmn,NZ)=0._r8
    WTGRBE(NB,ielmn,NZ)=0._r8
    WTNDBE(NB,ielmp,NZ)=0._r8
    WTHSKBE(NB,ielmp,NZ)=0._r8
    WTEARBE(NB,ielmp,NZ)=0._r8
    WTGRBE(NB,ielmp,NZ)=0._r8
    WTSTXB(NB,NZ)=0._r8
    WTSTXN(NB,NZ)=0._r8
    WTSTXP(NB,NZ)=0._r8
  ENDDO D8835
!
!     RESET ROOT STATE VARIABLES
!
  D6416: DO L=NU,NJ
    DO  N=1,MY(NZ)
      EPOOLR(:,N,L,NZ)=0._r8
      DO  NR=1,NRT(NZ)
        WTRT1E(ielmc,N,L,NR,NZ)=0._r8
        WTRT1E(ielmn,N,L,NR,NZ)=0._r8
        WTRT1E(ielmp,N,L,NR,NZ)=0._r8
        WTRT2E(ielmc,N,L,NR,NZ)=0._r8
        WTRT2E(ielmn,N,L,NR,NZ)=0._r8
        WTRT2E(ielmp,N,L,NR,NZ)=0._r8
        RTWT1(N,NR,NZ)=0._r8
        RTWT1N(N,NR,NZ)=0._r8
        RTWT1P(N,NR,NZ)=0._r8
        RTLG1(N,L,NR,NZ)=0._r8
        RTLG2(N,L,NR,NZ)=0._r8
        RTN2(N,L,NR,NZ)=0._r8
      enddo
    enddo
  ENDDO D6416
  WTRVE(ielmc,NZ)=0._r8
  WTRVE(ielmn,NZ)=0._r8
  WTRVE(ielmp,NZ)=0._r8
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
    WGLF     => plt_biom%WGLF       , &
    WTRSVBE  => plt_biom%WTRSVBE    , &
    WTSTXN   => plt_biom%WTSTXN     , &
    WTSTXP   => plt_biom%WTSTXP     , &
    WTSTXB   => plt_biom%WTSTXB     , &
    WVSTKB   => plt_biom%WVSTKB     , &
    WTGRBE   => plt_biom%WTGRBE     , &
    WTLSB    => plt_biom%WTLSB      , &
    WTHSKBE  => plt_biom%WTHSKBE    , &
    WTEARBE  => plt_biom%WTEARBE    , &
    WSSHE    => plt_biom%WSSHE      , &
    WGLFN    => plt_biom%WGLFN      , &
    WSLF     => plt_biom%WSLF       , &
    WGLFP    => plt_biom%WGLFP      , &
    WGLFLP   => plt_biom%WGLFLP     , &
    WGNODE   => plt_biom%WGNODE     , &
    WGNODN   => plt_biom%WGNODN     , &
    WGNODP   => plt_biom%WGNODP     , &
    WGSHP    => plt_biom%WGSHP      , &
    WGSHN    => plt_biom%WGSHN      , &
    WGSHE    => plt_biom%WGSHE      , &
    WGLFLN   => plt_biom%WGLFLN     , &
    WGLFL    => plt_biom%WGLFL      , &
    WGLFV    => plt_biom%WGLFV      , &
    WTSHEBE  => plt_biom%WTSHEBE    , &
    CPOOL3   => plt_photo%CPOOL3    , &
    CPOOL4   => plt_photo%CPOOL4    , &
    CO2B     => plt_photo%CO2B      , &
    HCOB     => plt_photo%HCOB      , &
    GRWTB    => plt_allom%GRWTB     , &
    GRNXB    => plt_morph%GRNXB     , &
    GRNOB    => plt_morph%GRNOB     , &
    ARLFB    => plt_morph%ARLFB     , &
    ARSTK    => plt_morph%ARSTK     , &
    NBR      => plt_morph%NBR       , &
    ARLF1    => plt_morph%ARLF1     , &
    HTSHE    => plt_morph%HTSHE     , &
    HTNODX   => plt_morph%HTNODX    , &
    SURF     => plt_morph%SURF      , &
    HTNODE   => plt_morph%HTNODE    , &
    ARLFV    => plt_morph%ARLFV     , &
    SURFB    => plt_morph%SURFB     , &
    ARLFL    => plt_morph%ARLFL       &
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
!     ARLFB=branch leaf area
!     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
!     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
!     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
!
  EPOOL(NB,ielmc,NZ)=0._r8
  CPOOLK(NB,NZ)=0._r8
  EPOOL(NB,ielmn,NZ)=0._r8
  EPOOL(NB,ielmp,NZ)=0._r8
  EPOLNB(NB,ielmc,NZ)=0._r8
  EPOLNB(NB,ielmn,NZ)=0._r8
  EPOLNB(NB,ielmp,NZ)=0._r8
  WTSHTBE(NB,ielmc,NZ)=0._r8
  WTLFBE(NB,ielmc,NZ)=0._r8
  WTNDBE(NB,ielmc,NZ)=0._r8
  WTSHEBE(NB,ielmc,NZ)=0._r8
  WTSTKBE(NB,ielmc,NZ)=0._r8
  WVSTKB(NB,NZ)=0._r8
  WTRSVBE(NB,ielmc,NZ)=0._r8
  WTHSKBE(NB,ielmc,NZ)=0._r8
  WTEARBE(NB,ielmc,NZ)=0._r8
  WTGRBE(NB,ielmc,NZ)=0._r8
  WTLSB(NB,NZ)=0._r8
  WTSHTBE(NB,ielmn,NZ)=0._r8
  WTLFBE(NB,ielmn,NZ)=0._r8
  WTNDBE(NB,ielmn,NZ)=0._r8
  WTSHEBE(NB,ielmn,NZ)=0._r8
  WTSTKBE(NB,ielmn,NZ)=0._r8
  WTRSVBE(NB,ielmn,NZ)=0._r8
  WTHSKBE(NB,ielmn,NZ)=0._r8
  WTEARBE(NB,ielmn,NZ)=0._r8
  WTGRBE(NB,ielmn,NZ)=0._r8
  WTSHTBE(NB,ielmp,NZ)=0._r8
  WTLFBE(NB,ielmp,NZ)=0._r8
  WTNDBE(NB,ielmp,NZ)=0._r8
  WTSHEBE(NB,ielmp,NZ)=0._r8
  WTSTKBE(NB,ielmp,NZ)=0._r8
  WTRSVBE(NB,ielmp,NZ)=0._r8
  WTHSKBE(NB,ielmp,NZ)=0._r8
  WTEARBE(NB,ielmp,NZ)=0._r8
  WTGRBE(NB,ielmp,NZ)=0._r8
  GRNXB(NB,NZ)=0._r8
  GRNOB(NB,NZ)=0._r8
  GRWTB(NB,NZ)=0._r8
  ARLFB(NB,NZ)=0._r8
  WTSTXB(NB,NZ)=0._r8
  WTSTXN(NB,NZ)=0._r8
  WTSTXP(NB,NZ)=0._r8
  DO 8855 K=0,JNODS1
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
    WGLF(K,NB,NZ)=0._r8
    WSLF(K,NB,NZ)=0._r8
    WGLFN(K,NB,NZ)=0._r8
    WGLFP(K,NB,NZ)=0._r8
    WGSHE(K,NB,NZ)=0._r8
    WSSHE(K,NB,NZ)=0._r8
    WGSHN(K,NB,NZ)=0._r8
    WGSHP(K,NB,NZ)=0._r8
    WGNODE(K,NB,NZ)=0._r8
    WGNODN(K,NB,NZ)=0._r8
    WGNODP(K,NB,NZ)=0._r8
    DO 8865 L=1,JC1
      ARLFV(L,NZ)=ARLFV(L,NZ)-ARLFL(L,K,NB,NZ)
      WGLFV(L,NZ)=WGLFV(L,NZ)-WGLFL(L,K,NB,NZ)
      ARLFL(L,K,NB,NZ)=0._r8
      WGLFL(L,K,NB,NZ)=0._r8
      WGLFLN(L,K,NB,NZ)=0._r8
      WGLFLP(L,K,NB,NZ)=0._r8
      IF(K.NE.0)THEN
        DO 8860 N=1,JLI1
          SURF(N,L,K,NB,NZ)=0._r8
8860    CONTINUE
      ENDIF
8865  CONTINUE
8855  CONTINUE
  DO 8875 L=1,JC1
    ARSTK(L,NB,NZ)=0._r8
    DO  N=1,JLI1
      SURFB(N,L,NB,NZ)=0._r8
    enddo
8875  CONTINUE
  end associate
  end subroutine ResetDeadRootStates


end module LitterFallMod
