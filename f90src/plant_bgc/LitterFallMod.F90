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
      IF(WTRVE(NZ,ielmc).LT.1.0E-04*WTRTE(NZ,ielmc) &
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
  associate(                                &
    JHVST      =>   plt_distb%JHVST   , &
    IYR0       =>   plt_distb%IYR0    , &
    IDAY0      =>   plt_distb%IDAY0   , &
    WTEARB     =>   plt_biom%WTEARB   , &
    WTEABN     =>   plt_biom%WTEABN   , &
    WTEABP     =>   plt_biom%WTEABP   , &
    CPOOL      =>   plt_biom%CPOOL    , &
    ZPOLNB     =>   plt_biom%ZPOLNB   , &
    PPOOL      =>   plt_biom%PPOOL    , &
    WTSTBP     =>   plt_biom%WTSTBP   , &
    ZPOOL      =>   plt_biom%ZPOOL    , &
    PPOLNB     =>   plt_biom%PPOLNB   , &
    WTSTKB     =>   plt_biom%WTSTKB   , &
    WTSTBN     =>   plt_biom%WTSTBN   , &
    WTRSBP     =>   plt_biom%WTRSBP   , &
    WTGRBN     =>   plt_biom%WTGRBN   , &
    WTGRBP     =>   plt_biom%WTGRBP   , &
    WTHSKB     =>   plt_biom%WTHSKB   , &
    WTSHEB     =>   plt_biom%WTSHEB   , &
    WTNDB      =>   plt_biom%WTNDB    , &
    WTSHBP     =>   plt_biom%WTSHBP   , &
    WTRSVB     =>   plt_biom%WTRSVB   , &
    WTLFBN     =>   plt_biom%WTLFBN   , &
    CPOLNB     =>   plt_biom%CPOLNB   , &
    WTRSBN     =>   plt_biom%WTRSBN   , &
    WTLFBP     =>   plt_biom%WTLFBP   , &
    WTNDBP     =>   plt_biom%WTNDBP   , &
    WTRVE      =>   plt_biom%WTRVE    , &
    WTGRB      =>   plt_biom%WTGRB    , &
    WTRT1      =>   plt_biom%WTRT1    , &
    WTRT1N     =>   plt_biom%WTRT1N   , &
    WTRT1P     =>   plt_biom%WTRT1P   , &
    WTLFB      =>   plt_biom%WTLFB    , &
    EPOOLR     =>   plt_biom%EPOOLR   , &
    WTNDBN     =>   plt_biom%WTNDBN   , &
    WTSHBN     =>   plt_biom%WTSHBN   , &
    WTHSBP     =>   plt_biom%WTHSBP   , &
    WTHSBN     =>   plt_biom%WTHSBN   , &
    WTSTDE     =>   plt_biom%WTSTDE   , &
    WTRT2      =>   plt_biom%WTRT2    , &
    WTRT2N     =>   plt_biom%WTRT2N   , &
    WTRT2P     =>   plt_biom%WTRT2P   , &
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
          ESNC(M,1,0,NZ,ielmc)=ESNC(M,1,0,NZ,ielmc) &
            +CFOPC(0,M,NZ)*(CPOOL(NB,NZ)+CPOLNB(NB,NZ) &
            +CPOOLK(NB,NZ)+WTRSVB(NB,NZ)) &
            +CFOPC(1,M,NZ)*(WTLFB(NB,NZ)*FWODB(1) &
            +WTNDB(NB,NZ)) &
            +CFOPC(2,M,NZ)*(WTSHEB(NB,NZ)*FWODB(1) &
            +WTHSKB(NB,NZ)+WTEARB(NB,NZ))
          ESNC(M,0,0,NZ,ielmc)=ESNC(M,0,0,NZ,ielmc) &
            +CFOPC(5,M,NZ)*(WTLFB(NB,NZ)*FWODB(0) &
            +WTSHEB(NB,NZ)*FWODB(0))
          ESNC(M,1,0,NZ,ielmn)=ESNC(M,1,0,NZ,ielmn) &
            +CFOPN(0,M,NZ)*(ZPOOL(NB,NZ)+ZPOLNB(NB,NZ) &
            +WTRSBN(NB,NZ)) &
            +CFOPN(1,M,NZ)*(WTLFBN(NB,NZ)*FWODLN(1) &
            +WTNDBN(NB,NZ)) &
            +CFOPN(2,M,NZ)*(WTSHBN(NB,NZ)*FWODSN(1) &
            +WTHSBN(NB,NZ)+WTEABN(NB,NZ))
          ESNC(M,0,0,NZ,ielmn)=ESNC(M,0,0,NZ,ielmn) &
            +CFOPN(5,M,NZ)*(WTLFBN(NB,NZ)*FWODLN(0) &
            +WTSHBN(NB,NZ)*FWODSN(0))
          ESNC(M,1,0,NZ,ielmp)=ESNC(M,1,0,NZ,ielmp) &
            +CFOPP(0,M,NZ)*(PPOOL(NB,NZ)+PPOLNB(NB,NZ) &
            +WTRSBP(NB,NZ)) &
            +CFOPP(1,M,NZ)*(WTLFBP(NB,NZ)*FWODLP(1) &
            +WTNDBP(NB,NZ)) &
            +CFOPP(2,M,NZ)*(WTSHBP(NB,NZ)*FWODSP(1) &
            +WTHSBP(NB,NZ)+WTEABP(NB,NZ))
          ESNC(M,0,0,NZ,ielmp)=ESNC(M,0,0,NZ,ielmp) &
            +CFOPP(5,M,NZ)*(WTLFBP(NB,NZ)*FWODLP(0) &
            +WTSHBP(NB,NZ)*FWODSP(0))
          IF(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).NE.0)THEN
            WTRVE(NZ,ielmc)=WTRVE(NZ,ielmc)+CFOPC(2,M,NZ)*WTGRB(NB,NZ)
            WTRVE(NZ,ielmn)=WTRVE(NZ,ielmn)+CFOPN(2,M,NZ)*WTGRBN(NB,NZ)
            WTRVE(NZ,ielmp)=WTRVE(NZ,ielmp)+CFOPP(2,M,NZ)*WTGRBP(NB,NZ)
          ELSE
            ESNC(M,1,0,NZ,ielmc)=ESNC(M,1,0,NZ,ielmc)+CFOPC(2,M,NZ)*WTGRB(NB,NZ)
            ESNC(M,1,0,NZ,ielmn)=ESNC(M,1,0,NZ,ielmn)+CFOPN(2,M,NZ)*WTGRBN(NB,NZ)
            ESNC(M,1,0,NZ,ielmp)=ESNC(M,1,0,NZ,ielmp)+CFOPP(2,M,NZ)*WTGRBP(NB,NZ)
          ENDIF
          IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
            ESNC(M,1,0,NZ,ielmc)=ESNC(M,1,0,NZ,ielmc)+CFOPC(3,M,NZ)*WTSTKB(NB,NZ)
            ESNC(M,1,0,NZ,ielmn)=ESNC(M,1,0,NZ,ielmn)+CFOPN(3,M,NZ)*WTSTBN(NB,NZ)
            ESNC(M,1,0,NZ,ielmp)=ESNC(M,1,0,NZ,ielmp)+CFOPP(3,M,NZ)*WTSTBP(NB,NZ)
          ELSE
            WTSTDE(M,NZ,ielmc)=WTSTDE(M,NZ,ielmc)+CFOPC(5,M,NZ)*WTSTKB(NB,NZ)
            WTSTDE(M,NZ,ielmn)=WTSTDE(M,NZ,ielmn)+CFOPN(5,M,NZ)*WTSTBN(NB,NZ)
            WTSTDE(M,NZ,ielmp)=WTSTDE(M,NZ,ielmp)+CFOPP(5,M,NZ)*WTSTBP(NB,NZ)
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
            ESNC(M,1,L,NZ,ielmc)=ESNC(M,1,L,NZ,ielmc)+CFOPC(0,M,NZ)*EPOOLR(N,L,NZ,ielmc)
            ESNC(M,1,L,NZ,ielmn)=ESNC(M,1,L,NZ,ielmn)+CFOPN(0,M,NZ)*EPOOLR(N,L,NZ,ielmn)
            ESNC(M,1,L,NZ,ielmp)=ESNC(M,1,L,NZ,ielmp)+CFOPP(0,M,NZ)*EPOOLR(N,L,NZ,ielmp)
            DO NR=1,NRT(NZ)
              ESNC(M,0,L,NZ,ielmc)=ESNC(M,0,L,NZ,ielmc)+CFOPC(5,M,NZ) &
                *(WTRT1(N,L,NR,NZ)+WTRT2(N,L,NR,NZ))*FWODR(0)
              ESNC(M,0,L,NZ,ielmn)=ESNC(M,0,L,NZ,ielmn)+CFOPN(5,M,NZ) &
                *(WTRT1N(N,L,NR,NZ)+WTRT2N(N,L,NR,NZ))*FWODRN(0)
              ESNC(M,0,L,NZ,ielmp)=ESNC(M,0,L,NZ,ielmp)+CFOPP(5,M,NZ) &
                *(WTRT1P(N,L,NR,NZ)+WTRT2P(N,L,NR,NZ))*FWODRP(0)
              ESNC(M,1,L,NZ,ielmc)=ESNC(M,1,L,NZ,ielmc)+CFOPC(4,M,NZ) &
                *(WTRT1(N,L,NR,NZ)+WTRT2(N,L,NR,NZ))*FWODR(1)
              ESNC(M,1,L,NZ,ielmn)=ESNC(M,1,L,NZ,ielmn)+CFOPN(4,M,NZ) &
                *(WTRT1N(N,L,NR,NZ)+WTRT2N(N,L,NR,NZ))*FWODRN(1)
              ESNC(M,1,L,NZ,ielmp)=ESNC(M,1,L,NZ,ielmp)+CFOPP(4,M,NZ) &
                *(WTRT1P(N,L,NR,NZ)+WTRT2P(N,L,NR,NZ))*FWODRP(1)
            ENDDO
          ENDDO
        ENDDO D6415
        ESNC(M,0,NG(NZ),NZ,ielmc)=ESNC(M,0,NG(NZ),NZ,ielmc) &
          +CFOPC(0,M,NZ)*WTRVE(NZ,ielmc)*FWOOD(0)
        ESNC(M,0,NG(NZ),NZ,ielmn)=ESNC(M,0,NG(NZ),NZ,ielmn) &
          +CFOPN(0,M,NZ)*WTRVE(NZ,ielmn)*FWOODN(0)
        ESNC(M,0,NG(NZ),NZ,ielmp)=ESNC(M,0,NG(NZ),NZ,ielmp) &
          +CFOPP(0,M,NZ)*WTRVE(NZ,ielmp)*FWOODP(0)
        ESNC(M,1,NG(NZ),NZ,ielmc)=ESNC(M,1,NG(NZ),NZ,ielmc) &
          +CFOPC(0,M,NZ)*WTRVE(NZ,ielmc)*FWOOD(1)
        ESNC(M,1,NG(NZ),NZ,ielmn)=ESNC(M,1,NG(NZ),NZ,ielmn) &
          +CFOPN(0,M,NZ)*WTRVE(NZ,ielmn)*FWOODN(1)
        ESNC(M,1,NG(NZ),NZ,ielmp)=ESNC(M,1,NG(NZ),NZ,ielmp) &
          +CFOPP(0,M,NZ)*WTRVE(NZ,ielmp)*FWOODP(1)
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
    WTRT2     =>   plt_biom%WTRT2     , &
    WTRT2N    =>   plt_biom%WTRT2N    , &
    WTRT2P    =>   plt_biom%WTRT2P    , &
    RTWT1     =>   plt_biom%RTWT1     , &
    RTWT1N    =>   plt_biom%RTWT1N    , &
    RTWT1P    =>   plt_biom%RTWT1P    , &
    WTRT1     =>   plt_biom%WTRT1     , &
    WTRT1N    =>   plt_biom%WTRT1N    , &
    WTRT1P    =>   plt_biom%WTRT1P    , &
    EPOOLR    =>   plt_biom%EPOOLR    , &
    WTRTD     =>   plt_biom%WTRTD     , &
    ZPOOLN    =>   plt_biom%ZPOOLN    , &
    WTRTL     =>   plt_biom%WTRTL     , &
    WSRTL     =>   plt_biom%WSRTL     , &
    CPOOLN    =>   plt_biom%CPOOLN    , &
    WTNDL     =>   plt_biom%WTNDL     , &
    WTNDLN    =>   plt_biom%WTNDLN    , &
    WTNDLP    =>   plt_biom%WTNDLP    , &
    PPOOLN    =>   plt_biom%PPOOLN    , &
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
          ESNC(M,1,L,NZ,ielmc)=ESNC(M,1,L,NZ,ielmc)+CFOPC(0,M,NZ)*EPOOLR(N,L,NZ,ielmc)
          ESNC(M,1,L,NZ,ielmn)=ESNC(M,1,L,NZ,ielmn)+CFOPN(0,M,NZ)*EPOOLR(N,L,NZ,ielmn)
          ESNC(M,1,L,NZ,ielmp)=ESNC(M,1,L,NZ,ielmp)+CFOPP(0,M,NZ)*EPOOLR(N,L,NZ,ielmp)
          DO  NR=1,NRT(NZ)
            ESNC(M,0,L,NZ,ielmc)=ESNC(M,0,L,NZ,ielmc)+CFOPC(5,M,NZ) &
              *(WTRT1(N,L,NR,NZ)+WTRT2(N,L,NR,NZ))*FWODR(0)
            ESNC(M,0,L,NZ,ielmn)=ESNC(M,0,L,NZ,ielmn)+CFOPN(5,M,NZ) &
              *(WTRT1N(N,L,NR,NZ)+WTRT2N(N,L,NR,NZ))*FWODRN(0)
            ESNC(M,0,L,NZ,ielmp)=ESNC(M,0,L,NZ,ielmp)+CFOPP(5,M,NZ) &
              *(WTRT1P(N,L,NR,NZ)+WTRT2P(N,L,NR,NZ))*FWODRP(0)
            ESNC(M,1,L,NZ,ielmc)=ESNC(M,1,L,NZ,ielmc)+CFOPC(4,M,NZ) &
              *(WTRT1(N,L,NR,NZ)+WTRT2(N,L,NR,NZ))*FWODR(1)
            ESNC(M,1,L,NZ,ielmn)=ESNC(M,1,L,NZ,ielmn)+CFOPN(4,M,NZ) &
              *(WTRT1N(N,L,NR,NZ)+WTRT2N(N,L,NR,NZ))*FWODRN(1)
            ESNC(M,1,L,NZ,ielmp)=ESNC(M,1,L,NZ,ielmp)+CFOPP(4,M,NZ) &
              *(WTRT1P(N,L,NR,NZ)+WTRT2P(N,L,NR,NZ))*FWODRP(1)
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
          WTRT1(N,L,NR,NZ)=0._r8
          WTRT1N(N,L,NR,NZ)=0._r8
          WTRT1P(N,L,NR,NZ)=0._r8
          WTRT2(N,L,NR,NZ)=0._r8
          WTRT2N(N,L,NR,NZ)=0._r8
          WTRT2P(N,L,NR,NZ)=0._r8
          RTWT1(N,NR,NZ)=0._r8
          RTWT1N(N,NR,NZ)=0._r8
          RTWT1P(N,NR,NZ)=0._r8
          RTLG1(N,L,NR,NZ)=0._r8
          RTLG2(N,L,NR,NZ)=0._r8
          RTN2(N,L,NR,NZ)=0._r8
        ENDDO D8870
        EPOOLR(N,L,NZ,:)=0._r8
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
            ESNC(M,1,L,NZ,ielmc)=ESNC(M,1,L,NZ,ielmc)+CFOPC(4,M,NZ) &
              *WTNDL(L,NZ)+CFOPC(0,M,NZ)*CPOOLN(L,NZ)
            ESNC(M,1,L,NZ,ielmn)=ESNC(M,1,L,NZ,ielmn)+CFOPN(4,M,NZ) &
              *WTNDLN(L,NZ)+CFOPN(0,M,NZ)*ZPOOLN(L,NZ)
            ESNC(M,1,L,NZ,ielmp)=ESNC(M,1,L,NZ,ielmp)+CFOPP(4,M,NZ) &
              *WTNDLP(L,NZ)+CFOPP(0,M,NZ)*PPOOLN(L,NZ)
          ENDDO D6420
          WTNDL(L,NZ)=0._r8
          WTNDLN(L,NZ)=0._r8
          WTNDLP(L,NZ)=0._r8
          CPOOLN(L,NZ)=0._r8
          ZPOOLN(L,NZ)=0._r8
          PPOOLN(L,NZ)=0._r8
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
    WTEABN    =>  plt_biom%WTEABN     , &
    CPOLNB    =>  plt_biom%CPOLNB     , &
    ZPOLNB    =>  plt_biom%ZPOLNB     , &
    PPOLNB    =>  plt_biom%PPOLNB     , &
    WTSTKB    =>  plt_biom%WTSTKB     , &
    WTHSKB    =>  plt_biom%WTHSKB     , &
    WTSHEB    =>  plt_biom%WTHSKB     , &
    WTHSBN    =>  plt_biom%WTHSBN     , &
    WTHSBP    =>  plt_biom%WTHSBP     , &
    WTSHBP    =>  plt_biom%WTSHBP     , &
    WTNDB     =>  plt_biom%WTNDB      , &
    WTLFB     =>  plt_biom%WTLFB      , &
    ZEROP     =>  plt_biom%ZEROP      , &
    WTLFBN    =>  plt_biom%WTLFBN     , &
    WTLFBP    =>  plt_biom%WTLFBP     , &
    WTNDBN    =>  plt_biom%WTNDBN     , &
    WTNDBP    =>  plt_biom%WTNDBP     , &
    WTSHBN    =>  plt_biom%WTSHBN     , &
    PPOOL     =>  plt_biom%PPOOL      , &
    WTEARB    =>  plt_biom%WTEARB     , &
    WTEABP    =>  plt_biom%WTEABP     , &
    WTSTBN    =>  plt_biom%WTSTBN     , &
    WTSTBP    =>  plt_biom%WTSTBP     , &
    ZPOOL     =>  plt_biom%ZPOOL      , &
    CPOOL     =>  plt_biom%CPOOL      , &
    WTGRBP    =>  plt_biom%WTGRBP     , &
    WTGRB     =>  plt_biom%WTGRB      , &
    WTRSVB    =>  plt_biom%WTRSVB     , &
    WTRSBN    =>  plt_biom%WTRSBN     , &
    WTRSBP    =>  plt_biom%WTRSBP     , &
    WTSTDE    =>  plt_biom%WTSTDE     , &
    WTGRBN    =>  plt_biom%WTGRBN     , &
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
        ESNC(M,1,0,NZ,ielmc)=ESNC(M,1,0,NZ,ielmc) &
          +CFOPC(0,M,NZ)*CPOLNB(NB,NZ) &
          +CFOPC(1,M,NZ)*(WTLFB(NB,NZ)*FWODB(1) &
          +WTNDB(NB,NZ)) &
          +CFOPC(2,M,NZ)*(WTSHEB(NB,NZ)*FWODB(1) &
          +WTHSKB(NB,NZ)+WTEARB(NB,NZ))
        ESNC(M,0,0,NZ,ielmc)=ESNC(M,0,0,NZ,ielmc) &
          +CFOPC(5,M,NZ)*(WTLFB(NB,NZ)*FWODB(0) &
          +WTSHEB(NB,NZ)*FWODB(0))
        ESNC(M,1,0,NZ,ielmn)=ESNC(M,1,0,NZ,ielmn) &
          +CFOPN(0,M,NZ)*ZPOLNB(NB,NZ) &
          +CFOPN(1,M,NZ)*(WTLFBN(NB,NZ)*FWODLN(1) &
          +WTNDBN(NB,NZ)) &
          +CFOPN(2,M,NZ)*(WTSHBN(NB,NZ)*FWODSN(1) &
          +WTHSBN(NB,NZ)+WTEABN(NB,NZ))
        ESNC(M,0,0,NZ,ielmn)=ESNC(M,0,0,NZ,ielmn) &
          +CFOPN(5,M,NZ)*(WTLFBN(NB,NZ)*FWODLN(0) &
          +WTSHBN(NB,NZ)*FWODSN(0))
        ESNC(M,1,0,NZ,ielmp)=ESNC(M,1,0,NZ,ielmp) &
          +CFOPP(0,M,NZ)*PPOLNB(NB,NZ) &
          +CFOPP(1,M,NZ)*(WTLFBP(NB,NZ)*FWODLP(1) &
          +WTNDBP(NB,NZ)) &
          +CFOPP(2,M,NZ)*(WTSHBP(NB,NZ)*FWODSP(1) &
          +WTHSBP(NB,NZ)+WTEABP(NB,NZ))
        ESNC(M,0,0,NZ,ielmp)=ESNC(M,0,0,NZ,ielmp) &
          +CFOPP(5,M,NZ)*(WTLFBP(NB,NZ)*FWODLP(0) &
          +WTSHBP(NB,NZ)*FWODSP(0))
        IF(ISTYP(NZ).EQ.0.AND.IWTYP(NZ).NE.0)THEN
          WTRVE(NZ,ielmc)=WTRVE(NZ,ielmc)+CFOPC(2,M,NZ)*WTGRB(NB,NZ)
          WTRVE(NZ,ielmn)=WTRVE(NZ,ielmn)+CFOPN(2,M,NZ)*WTGRBN(NB,NZ)
          WTRVE(NZ,ielmp)=WTRVE(NZ,ielmp)+CFOPP(2,M,NZ)*WTGRBP(NB,NZ)
        ELSE
          ESNC(M,1,0,NZ,ielmc)=ESNC(M,1,0,NZ,ielmc)+CFOPC(2,M,NZ)*WTGRB(NB,NZ)
          ESNC(M,1,0,NZ,ielmn)=ESNC(M,1,0,NZ,ielmn)+CFOPN(2,M,NZ)*WTGRBN(NB,NZ)
          ESNC(M,1,0,NZ,ielmp)=ESNC(M,1,0,NZ,ielmp)+CFOPP(2,M,NZ)*WTGRBP(NB,NZ)
        ENDIF
        IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
          ESNC(M,1,0,NZ,ielmc)=ESNC(M,1,0,NZ,ielmc)+CFOPC(3,M,NZ)*WTSTKB(NB,NZ)
          ESNC(M,1,0,NZ,ielmn)=ESNC(M,1,0,NZ,ielmn)+CFOPN(3,M,NZ)*WTSTBN(NB,NZ)
          ESNC(M,1,0,NZ,ielmp)=ESNC(M,1,0,NZ,ielmp)+CFOPP(3,M,NZ)*WTSTBP(NB,NZ)
        ELSE
          WTSTDE(M,NZ,ielmc)=WTSTDE(M,NZ,ielmc)+CFOPC(5,M,NZ)*WTSTKB(NB,NZ)
          WTSTDE(M,NZ,ielmn)=WTSTDE(M,NZ,ielmn)+CFOPN(5,M,NZ)*WTSTBN(NB,NZ)
          WTSTDE(M,NZ,ielmp)=WTSTDE(M,NZ,ielmp)+CFOPP(5,M,NZ)*WTSTBP(NB,NZ)
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
      WTRVE(NZ,ielmc)=WTRVE(NZ,ielmc)+CPOOL(NB,NZ)+CPOOLK(NB,NZ)
      WTRVE(NZ,ielmn)=WTRVE(NZ,ielmn)+ZPOOL(NB,NZ)
      WTRVE(NZ,ielmp)=WTRVE(NZ,ielmp)+PPOOL(NB,NZ)
      IF(IHVST(NZ).NE.4.AND.IHVST(NZ).NE.6)THEN
        D6406: DO M=1,jsken
          ESNC(M,1,0,NZ,ielmc)=ESNC(M,1,0,NZ,ielmc) &
            +CFOPC(0,M,NZ)*WTRSVB(NB,NZ)
          ESNC(M,1,0,NZ,ielmn)=ESNC(M,1,0,NZ,ielmn) &
            +CFOPN(0,M,NZ)*WTRSBN(NB,NZ)
          ESNC(M,1,0,NZ,ielmp)=ESNC(M,1,0,NZ,ielmp) &
            +CFOPP(0,M,NZ)*WTRSBP(NB,NZ)
        ENDDO D6406
      ELSE
        WTRVE(NZ,ielmc)=WTRVE(NZ,ielmc)+WTRSVB(NB,NZ)
        WTRVE(NZ,ielmn)=WTRVE(NZ,ielmn)+WTRSBN(NB,NZ)
        WTRVE(NZ,ielmp)=WTRVE(NZ,ielmp)+WTRSBP(NB,NZ)
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
  integer :: L,NR,N,NB
!     begin_execution
  associate(                               &
    CPOOL  => plt_biom%CPOOL         , &
    ZPOOL  => plt_biom%ZPOOL         , &
    CPOLNB => plt_biom%CPOLNB        , &
    ZPOLNB => plt_biom%ZPOLNB        , &
    WTSHTB => plt_biom%WTSHTB        , &
    WTLFB  => plt_biom%WTLFB         , &
    WTRSVB => plt_biom%WTRSVB        , &
    WTSHEB => plt_biom%WTSHEB        , &
    WTHSKB => plt_biom%WTHSKB        , &
    WTEARB => plt_biom%WTEARB        , &
    WTGRB  => plt_biom%WTGRB         , &
    WVSTKB => plt_biom%WVSTKB        , &
    WTRSBN => plt_biom%WTRSBN        , &
    WTHSBN => plt_biom%WTHSBN        , &
    EPOOLR => plt_biom%EPOOLR        , &
    WTEABN => plt_biom%WTEABN        , &
    WTGRBN => plt_biom%WTGRBN        , &
    WTSHTP => plt_biom%WTSHTP        , &
    WTLFBP => plt_biom%WTLFBP        , &
    WTNDBP => plt_biom%WTNDBP        , &
    WTSHBP => plt_biom%WTSHBP        , &
    WTSTBP => plt_biom%WTSTBP        , &
    WTRSBP => plt_biom%WTRSBP        , &
    WTHSBP => plt_biom%WTHSBP        , &
    WTEABP => plt_biom%WTEABP        , &
    WTGRBP => plt_biom%WTGRBP        , &
    WTSTXB => plt_biom%WTSTXB        , &
    WTSTXN => plt_biom%WTSTXN        , &
    WTSTXP => plt_biom%WTSTXP        , &
    WTSTBN => plt_biom%WTSTBN        , &
    WTLSB  => plt_biom%WTLSB         , &
    WTSHBN => plt_biom%WTSHBN        , &
    WTRVE  => plt_biom%WTRVE         , &
    WTLFBN => plt_biom%WTLFBN        , &
    WTSHTN => plt_biom%WTSHTN        , &
    WTRT1  => plt_biom%WTRT1         , &
    WTRT1N => plt_biom%WTRT1N        , &
    WTRT1P => plt_biom%WTRT1P        , &
    WTNDBN => plt_biom%WTNDBN        , &
    WTRT2  => plt_biom%WTRT2         , &
    WTRT2N => plt_biom%WTRT2N        , &
    WTRT2P => plt_biom%WTRT2P        , &
    WTNDB  => plt_biom%WTNDB         , &
    WTSTKB => plt_biom%WTSTKB        , &
    RTWT1  => plt_biom%RTWT1         , &
    RTWT1N => plt_biom%RTWT1N        , &
    RTWT1P => plt_biom%RTWT1P        , &
    PPOLNB => plt_biom%PPOLNB        , &
    PPOOL  => plt_biom%PPOOL         , &
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
  DO 8835 NB=1,NBR(NZ)
    CPOOL(NB,NZ)=0._r8
    CPOOLK(NB,NZ)=0._r8
    ZPOOL(NB,NZ)=0._r8
    PPOOL(NB,NZ)=0._r8
    CPOLNB(NB,NZ)=0._r8
    ZPOLNB(NB,NZ)=0._r8
    PPOLNB(NB,NZ)=0._r8
    WTSHTB(NB,NZ)=0._r8
    WTLFB(NB,NZ)=0._r8
    WTNDB(NB,NZ)=0._r8
    WTSHEB(NB,NZ)=0._r8
    WTSTKB(NB,NZ)=0._r8
    WVSTKB(NB,NZ)=0._r8
    WTRSVB(NB,NZ)=0._r8
    WTHSKB(NB,NZ)=0._r8
    WTEARB(NB,NZ)=0._r8
    WTGRB(NB,NZ)=0._r8
    WTLSB(NB,NZ)=0._r8
    WTSHTN(NB,NZ)=0._r8
    WTLFBN(NB,NZ)=0._r8
    WTNDBN(NB,NZ)=0._r8
    WTSHBN(NB,NZ)=0._r8
    WTSTBN(NB,NZ)=0._r8
    WTRSBN(NB,NZ)=0._r8
    WTHSBN(NB,NZ)=0._r8
    WTEABN(NB,NZ)=0._r8
    WTGRBN(NB,NZ)=0._r8
    WTSHTP(NB,NZ)=0._r8
    WTLFBP(NB,NZ)=0._r8
    WTNDBP(NB,NZ)=0._r8
    WTSHBP(NB,NZ)=0._r8
    WTSTBP(NB,NZ)=0._r8
    WTRSBP(NB,NZ)=0._r8
    WTHSBP(NB,NZ)=0._r8
    WTEABP(NB,NZ)=0._r8
    WTGRBP(NB,NZ)=0._r8
    WTSTXB(NB,NZ)=0._r8
    WTSTXN(NB,NZ)=0._r8
    WTSTXP(NB,NZ)=0._r8
8835  CONTINUE
!
!     RESET ROOT STATE VARIABLES
!
  D6416: DO L=NU,NJ
    DO  N=1,MY(NZ)
      EPOOLR(N,L,NZ,:)=0._r8
      DO  NR=1,NRT(NZ)
        WTRT1(N,L,NR,NZ)=0._r8
        WTRT1N(N,L,NR,NZ)=0._r8
        WTRT1P(N,L,NR,NZ)=0._r8
        WTRT2(N,L,NR,NZ)=0._r8
        WTRT2N(N,L,NR,NZ)=0._r8
        WTRT2P(N,L,NR,NZ)=0._r8
        RTWT1(N,NR,NZ)=0._r8
        RTWT1N(N,NR,NZ)=0._r8
        RTWT1P(N,NR,NZ)=0._r8
        RTLG1(N,L,NR,NZ)=0._r8
        RTLG2(N,L,NR,NZ)=0._r8
        RTN2(N,L,NR,NZ)=0._r8
      enddo
    enddo
  ENDDO D6416
  WTRVE(NZ,ielmc)=0._r8
  WTRVE(NZ,ielmn)=0._r8
  WTRVE(NZ,ielmp)=0._r8
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
  associate(                              &
    CPOOL    => plt_biom%CPOOL      , &
    ZPOOL    => plt_biom%ZPOOL      , &
    PPOOL    => plt_biom%PPOOL      , &
    CPOLNB   => plt_biom%CPOLNB     , &
    PPOLNB   => plt_biom%PPOLNB     , &
    WTSHTB   => plt_biom%WTSHTB     , &
    WTLFB    => plt_biom%WTLFB      , &
    WTNDB    => plt_biom%WTNDB      , &
    WTSTKB   => plt_biom%WTSTKB     , &
    WGLF     => plt_biom%WGLF       , &
    WTRSVB   => plt_biom%WTRSVB     , &
    WTRSBP   => plt_biom%WTRSBP     , &
    WTHSBP   => plt_biom%WTHSBP     , &
    WTSTXN   => plt_biom%WTSTXN     , &
    WTSTXP   => plt_biom%WTSTXP     , &
    WTSTXB   => plt_biom%WTSTXB     , &
    WTSHBP   => plt_biom%WTSHBP     , &
    WTEABP   => plt_biom%WTEABP     , &
    WTGRBP   => plt_biom%WTGRBP     , &
    WTLFBP   => plt_biom%WTLFBP     , &
    WTNDBP   => plt_biom%WTNDBP     , &
    WTSHTN   => plt_biom%WTSHTN     , &
    WTSHBN   => plt_biom%WTSHBN     , &
    WTSTBP   => plt_biom%WTSTBP     , &
    WTNDBN   => plt_biom%WTNDBN     , &
    WTSTBN   => plt_biom%WTSTBN     , &
    WTSHTP   => plt_biom%WTSHTP     , &
    WTEABN   => plt_biom%WTEABN     , &
    WTGRBN   => plt_biom%WTGRBN     , &
    WTHSBN   => plt_biom%WTHSBN     , &
    WTRSBN   => plt_biom%WTRSBN     , &
    WTLFBN   => plt_biom%WTLFBN     , &
    WVSTKB   => plt_biom%WVSTKB     , &
    WTGRB    => plt_biom%WTGRB      , &
    WTLSB    => plt_biom%WTLSB      , &
    WTHSKB   => plt_biom%WTHSKB     , &
    WTEARB   => plt_biom%WTEARB     , &
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
    WTSHEB   => plt_biom%WTSHEB     , &
    ZPOLNB   => plt_biom%ZPOLNB     , &
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
  CPOOL(NB,NZ)=0._r8
  CPOOLK(NB,NZ)=0._r8
  ZPOOL(NB,NZ)=0._r8
  PPOOL(NB,NZ)=0._r8
  CPOLNB(NB,NZ)=0._r8
  ZPOLNB(NB,NZ)=0._r8
  PPOLNB(NB,NZ)=0._r8
  WTSHTB(NB,NZ)=0._r8
  WTLFB(NB,NZ)=0._r8
  WTNDB(NB,NZ)=0._r8
  WTSHEB(NB,NZ)=0._r8
  WTSTKB(NB,NZ)=0._r8
  WVSTKB(NB,NZ)=0._r8
  WTRSVB(NB,NZ)=0._r8
  WTHSKB(NB,NZ)=0._r8
  WTEARB(NB,NZ)=0._r8
  WTGRB(NB,NZ)=0._r8
  WTLSB(NB,NZ)=0._r8
  WTSHTN(NB,NZ)=0._r8
  WTLFBN(NB,NZ)=0._r8
  WTNDBN(NB,NZ)=0._r8
  WTSHBN(NB,NZ)=0._r8
  WTSTBN(NB,NZ)=0._r8
  WTRSBN(NB,NZ)=0._r8
  WTHSBN(NB,NZ)=0._r8
  WTEABN(NB,NZ)=0._r8
  WTGRBN(NB,NZ)=0._r8
  WTSHTP(NB,NZ)=0._r8
  WTLFBP(NB,NZ)=0._r8
  WTNDBP(NB,NZ)=0._r8
  WTSHBP(NB,NZ)=0._r8
  WTSTBP(NB,NZ)=0._r8
  WTRSBP(NB,NZ)=0._r8
  WTHSBP(NB,NZ)=0._r8
  WTEABP(NB,NZ)=0._r8
  WTGRBP(NB,NZ)=0._r8
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
