module LitterFallMod

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
implicit none
  private
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
    IYRHs1       =>   plt_distb%IYRHs1     , &
    JHVSTs1      =>   plt_distb%JHVSTs1    , &
    IDAYHs1      =>   plt_distb%IDAYHs1    , &
    UVOLOs1      =>   plt_ew%UVOLOs1       , &
    VOLWPs1      =>   plt_ew%VOLWPs1       , &
    WTRTs1       =>   plt_biom%WTRTs1      , &
    WTRVCs1      =>   plt_biom%WTRVCs1     , &
    ISTYPs1      =>   plt_pheno%ISTYPs1    , &
    IDTHRs1      =>   plt_pheno%IDTHRs1    , &
    IDTHPs1      =>   plt_pheno%IDTHPs1    , &
    IFLGIs1      =>   plt_pheno%IFLGIs1    , &
    WSTRs1       =>   plt_pheno%WSTRs1     , &
    IDAYs1       =>   plt_pheno%IDAYs1     , &
    PPs1         =>   plt_site%PPs1        , &
    IYRCs1       =>   plt_site%IYRCs1      , &
    VOLWOUs1     =>   plt_site%VOLWOUs1    , &
    ZNOONs1      =>   plt_site%ZNOONs1     , &
    HTCTLs1      =>   plt_morph%HTCTLs1    , &
    NBRs1        =>   plt_morph%NBRs1      , &
    NB1s1        =>   plt_morph%NB1s1      , &
    NBTs1        =>   plt_morph%NBTs1        &
  )
!
!     ZNOON=hour of solar noon
!     IDAYs1(1,=emergence date
!     ISTYP=growth habit:0=annual,1=perennial from PFT file
!     IDAYH,IYRH=day,year of harvesting
!     IYRCs1=current year
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
  IF(J.EQ.INT(ZNOONs1).AND.IDAYs1(1,NB1s1(NZ),NZ).NE.0 &
    .AND.(ISTYPs1(NZ).NE.0.OR.(I.GE.IDAYHs1(NZ) &
    .AND.IYRCs1.GE.IYRHs1(NZ))))THEN
    IDTHY=0
!
!     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
!
    call LiterfallFromDeadBranches(I,J,NZ,IDTHY,CPOOLK)

    IF(IDTHY.EQ.NBRs1(NZ))THEN
      IDTHPs1(NZ)=1
      NBTs1(NZ)=0
      WSTRs1(NZ)=0._r8
      IF(IFLGIs1(NZ).EQ.1)THEN
        NBRs1(NZ)=1
      ELSE
        NBRs1(NZ)=0
      ENDIF
      HTCTLs1(NZ)=0._r8
      VOLWOUs1=VOLWOUs1+VOLWPs1(NZ)
      UVOLOs1=UVOLOs1+VOLWPs1(NZ)
      VOLWPs1(NZ)=0._r8
!
!     RESET LIVING FLAGS
!
!     WTRVC,WTRT=PFT storage,root C
!     ISTYP=growth habit:0=annual,1=perennial
!     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
!     PP=PFT population
!     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
!
      IF(WTRVCs1(NZ).LT.1.0E-04*WTRTs1(NZ) &
        .AND.ISTYPs1(NZ).NE.0)IDTHRs1(NZ)=1
      IF(ISTYPs1(NZ).EQ.0)IDTHRs1(NZ)=1
      IF(JHVSTs1(NZ).NE.0)IDTHRs1(NZ)=1
      IF(PPs1(NZ).LE.0.0)IDTHRs1(NZ)=1
      IF(IDTHRs1(NZ).EQ.1)IDTHPs1(NZ)=1
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
    JHVSTs1      =>   plt_distb%JHVSTs1   , &
    IYR0s1       =>   plt_distb%IYR0s1    , &
    IDAY0s1      =>   plt_distb%IDAY0s1   , &
    WTEARBs1     =>   plt_biom%WTEARBs1   , &
    WTEABNs1     =>   plt_biom%WTEABNs1   , &
    WTEABPs1     =>   plt_biom%WTEABPs1   , &
    CPOOLs1      =>   plt_biom%CPOOLs1    , &
    ZPOLNBs1     =>   plt_biom%ZPOLNBs1   , &
    PPOOLs1      =>   plt_biom%PPOOLs1    , &
    WTSTBPs1     =>   plt_biom%WTSTBPs1   , &
    ZPOOLs1      =>   plt_biom%ZPOOLs1    , &
    PPOLNBs1     =>   plt_biom%PPOLNBs1   , &
    WTSTKBs1     =>   plt_biom%WTSTKBs1   , &
    WTSTBNs1     =>   plt_biom%WTSTBNs1   , &
    WTRSBPs1     =>   plt_biom%WTRSBPs1   , &
    WTGRBNs1     =>   plt_biom%WTGRBNs1   , &
    WTRVPs1      =>   plt_biom%WTRVPs1    , &
    WTGRBPs1     =>   plt_biom%WTGRBPs1   , &
    WTHSKBs1     =>   plt_biom%WTHSKBs1   , &
    WTSHEBs1     =>   plt_biom%WTSHEBs1   , &
    WTNDBs1      =>   plt_biom%WTNDBs1    , &
    WTSHBPs1     =>   plt_biom%WTSHBPs1   , &
    WTRSVBs1     =>   plt_biom%WTRSVBs1   , &
    WTLFBNs1     =>   plt_biom%WTLFBNs1   , &
    CPOLNBs1     =>   plt_biom%CPOLNBs1   , &
    WTRSBNs1     =>   plt_biom%WTRSBNs1   , &
    WTLFBPs1     =>   plt_biom%WTLFBPs1   , &
    WTRVNs1      =>   plt_biom%WTRVNs1    , &
    WTNDBPs1     =>   plt_biom%WTNDBPs1   , &
    WTRVCs1      =>   plt_biom%WTRVCs1    , &
    WTGRBs1      =>   plt_biom%WTGRBs1    , &
    WTRT1s1      =>   plt_biom%WTRT1s1    , &
    WTRT1Ns1     =>   plt_biom%WTRT1Ns1   , &
    WTRT1Ps1     =>   plt_biom%WTRT1Ps1   , &
    WTLFBs1      =>   plt_biom%WTLFBs1    , &
    CPOOLRs1     =>   plt_biom%CPOOLRs1   , &
    ZPOOLRs1     =>   plt_biom%ZPOOLRs1   , &
    PPOOLRs1     =>   plt_biom%PPOOLRs1   , &
    WTSTDNs1     =>   plt_biom%WTSTDNs1   , &
    WTSTDPs1     =>   plt_biom%WTSTDPs1   , &
    WTNDBNs1     =>   plt_biom%WTNDBNs1   , &
    WTSHBNs1     =>   plt_biom%WTSHBNs1   , &
    WTHSBPs1     =>   plt_biom%WTHSBPs1   , &
    WTHSBNs1     =>   plt_biom%WTHSBNs1   , &
    WTSTDGs1     =>   plt_biom%WTSTDGs1   , &
    WTRT2s1      =>   plt_biom%WTRT2s1    , &
    WTRT2Ns1     =>   plt_biom%WTRT2Ns1   , &
    WTRT2Ps1     =>   plt_biom%WTRT2Ps1   , &
    FWODLNs1     =>   plt_allom%FWODLNs1  , &
    FWODLPs1     =>   plt_allom%FWODLPs1  , &
    FWODBs1      =>   plt_allom%FWODBs1   , &
    FWODSPs1     =>   plt_allom%FWODSPs1  , &
    FWODSNs1     =>   plt_allom%FWODSNs1  , &
    FWOODs1      =>   plt_allom%FWOODs1   , &
    FWODRs1      =>   plt_allom%FWODRs1   , &
    FWODRNs1     =>   plt_allom%FWODRNs1  , &
    FWODRPs1     =>   plt_allom%FWODRPs1  , &
    FWOODNs1     =>   plt_allom%FWOODNs1  , &
    FWOODPs1     =>   plt_allom%FWOODPs1  , &
    IFLGIs1      =>   plt_pheno%IFLGIs1   , &
    ISTYPs1      =>   plt_pheno%ISTYPs1   , &
    IWTYPs1      =>   plt_pheno%IWTYPs1   , &
    IDTHRs1      =>   plt_pheno%IDTHRs1   , &
    IGTYPs1      =>   plt_pheno%IGTYPs1   , &
    IBTYPs1      =>   plt_pheno%IBTYPs1   , &
    IDTHPs1      =>   plt_pheno%IDTHPs1   , &
    CSNCs1       =>   plt_bgcr%CSNCs1     , &
    ZSNCs1       =>   plt_bgcr%ZSNCs1     , &
    PSNCs1       =>   plt_bgcr%PSNCs1     , &
    LYRCs1       =>   plt_site%LYRCs1     , &
    NJs1         =>   plt_site%NJs1       , &
    NUs1         =>   plt_site%NUs1       , &
    IDATAs1      =>   plt_site%IDATAs1    , &
    CFOPPs1      =>   plt_soilchem%CFOPPs1, &
    CFOPNs1      =>   plt_soilchem%CFOPNs1, &
    CFOPCs1      =>   plt_soilchem%CFOPCs1, &
    MYs1         =>   plt_morph%MYs1      , &
    NGs1         =>   plt_morph%NGs1      , &
    NBRs1        =>   plt_morph%NBRs1     , &
    NRTs1        =>   plt_morph%NRTs1       &
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
  IF(IDTHPs1(NZ).EQ.1.AND.IDTHRs1(NZ).EQ.1)THEN
    IF(IFLGIs1(NZ).EQ.0)THEN
      DO 6425 M=1,4
        DO 8825 NB=1,NBRs1(NZ)
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
            +CFOPCs1(0,M,NZ)*(CPOOLs1(NB,NZ)+CPOLNBs1(NB,NZ) &
            +CPOOLK(NB,NZ)+WTRSVBs1(NB,NZ)) &
            +CFOPCs1(1,M,NZ)*(WTLFBs1(NB,NZ)*FWODBs1(1) &
            +WTNDBs1(NB,NZ)) &
            +CFOPCs1(2,M,NZ)*(WTSHEBs1(NB,NZ)*FWODBs1(1) &
            +WTHSKBs1(NB,NZ)+WTEARBs1(NB,NZ))
          CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ) &
            +CFOPCs1(5,M,NZ)*(WTLFBs1(NB,NZ)*FWODBs1(0) &
            +WTSHEBs1(NB,NZ)*FWODBs1(0))
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
            +CFOPNs1(0,M,NZ)*(ZPOOLs1(NB,NZ)+ZPOLNBs1(NB,NZ) &
            +WTRSBNs1(NB,NZ)) &
            +CFOPNs1(1,M,NZ)*(WTLFBNs1(NB,NZ)*FWODLNs1(1) &
            +WTNDBNs1(NB,NZ)) &
            +CFOPNs1(2,M,NZ)*(WTSHBNs1(NB,NZ)*FWODSNs1(1) &
            +WTHSBNs1(NB,NZ)+WTEABNs1(NB,NZ))
          ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ) &
            +CFOPNs1(5,M,NZ)*(WTLFBNs1(NB,NZ)*FWODLNs1(0) &
            +WTSHBNs1(NB,NZ)*FWODSNs1(0))
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
            +CFOPPs1(0,M,NZ)*(PPOOLs1(NB,NZ)+PPOLNBs1(NB,NZ) &
            +WTRSBPs1(NB,NZ)) &
            +CFOPPs1(1,M,NZ)*(WTLFBPs1(NB,NZ)*FWODLPs1(1) &
            +WTNDBPs1(NB,NZ)) &
            +CFOPPs1(2,M,NZ)*(WTSHBPs1(NB,NZ)*FWODSPs1(1) &
            +WTHSBPs1(NB,NZ)+WTEABPs1(NB,NZ))
          PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ) &
            +CFOPPs1(5,M,NZ)*(WTLFBPs1(NB,NZ)*FWODLPs1(0) &
            +WTSHBPs1(NB,NZ)*FWODSPs1(0))
          IF(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).NE.0)THEN
            WTRVCs1(NZ)=WTRVCs1(NZ)+CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
            WTRVNs1(NZ)=WTRVNs1(NZ)+CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
            WTRVPs1(NZ)=WTRVPs1(NZ)+CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
          ELSE
            CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
            ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
            PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
          ENDIF
          IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
            CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(3,M,NZ)*WTSTKBs1(NB,NZ)
            ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(3,M,NZ)*WTSTBNs1(NB,NZ)
            PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(3,M,NZ)*WTSTBPs1(NB,NZ)
          ELSE
            WTSTDGs1(M,NZ)=WTSTDGs1(M,NZ)+CFOPCs1(5,M,NZ)*WTSTKBs1(NB,NZ)
            WTSTDNs1(M,NZ)=WTSTDNs1(M,NZ)+CFOPNs1(5,M,NZ)*WTSTBNs1(NB,NZ)
            WTSTDPs1(M,NZ)=WTSTDPs1(M,NZ)+CFOPPs1(5,M,NZ)*WTSTBPs1(NB,NZ)
          ENDIF
8825    CONTINUE
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
        DO 6415 L=NUs1,NJs1
          DO N=1,MYs1(NZ)
            CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(0,M,NZ)*CPOOLRs1(N,L,NZ)
            ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(0,M,NZ)*ZPOOLRs1(N,L,NZ)
            PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(0,M,NZ)*PPOOLRs1(N,L,NZ)
            DO NR=1,NRTs1(NZ)
              CSNCs1(M,0,L,NZ)=CSNCs1(M,0,L,NZ)+CFOPCs1(5,M,NZ) &
                *(WTRT1s1(N,L,NR,NZ)+WTRT2s1(N,L,NR,NZ))*FWODRs1(0)
              ZSNCs1(M,0,L,NZ)=ZSNCs1(M,0,L,NZ)+CFOPNs1(5,M,NZ) &
                *(WTRT1Ns1(N,L,NR,NZ)+WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(0)
              PSNCs1(M,0,L,NZ)=PSNCs1(M,0,L,NZ)+CFOPPs1(5,M,NZ) &
                *(WTRT1Ps1(N,L,NR,NZ)+WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(0)
              CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ) &
                *(WTRT1s1(N,L,NR,NZ)+WTRT2s1(N,L,NR,NZ))*FWODRs1(1)
              ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ) &
                *(WTRT1Ns1(N,L,NR,NZ)+WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(1)
              PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ) &
                *(WTRT1Ps1(N,L,NR,NZ)+WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(1)
            ENDDO
          ENDDO
6415    CONTINUE
        CSNCs1(M,0,NGs1(NZ),NZ)=CSNCs1(M,0,NGs1(NZ),NZ) &
          +CFOPCs1(0,M,NZ)*WTRVCs1(NZ)*FWOODs1(0)
        ZSNCs1(M,0,NGs1(NZ),NZ)=ZSNCs1(M,0,NGs1(NZ),NZ) &
          +CFOPNs1(0,M,NZ)*WTRVNs1(NZ)*FWOODNs1(0)
        PSNCs1(M,0,NGs1(NZ),NZ)=PSNCs1(M,0,NGs1(NZ),NZ) &
          +CFOPPs1(0,M,NZ)*WTRVPs1(NZ)*FWOODPs1(0)
        CSNCs1(M,1,NGs1(NZ),NZ)=CSNCs1(M,1,NGs1(NZ),NZ) &
          +CFOPCs1(0,M,NZ)*WTRVCs1(NZ)*FWOODs1(1)
        ZSNCs1(M,1,NGs1(NZ),NZ)=ZSNCs1(M,1,NGs1(NZ),NZ) &
          +CFOPNs1(0,M,NZ)*WTRVNs1(NZ)*FWOODNs1(1)
        PSNCs1(M,1,NGs1(NZ),NZ)=PSNCs1(M,1,NGs1(NZ),NZ) &
          +CFOPPs1(0,M,NZ)*WTRVPs1(NZ)*FWOODPs1(1)
6425  CONTINUE
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
    IF(ISTYPs1(NZ).NE.0.AND.JHVSTs1(NZ).EQ.0)THEN
      IF(I.LT.LYRCs1)THEN
        IDAY0s1(NZ)=I+1
        IYR0s1(NZ)=IDATAs1(3)
      ELSE
        IDAY0s1(NZ)=1
        IYR0s1(NZ)=IDATAs1(3)+1
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
    WTRT2s1     =>   plt_biom%WTRT2s1     , &
    WTRT2Ns1    =>   plt_biom%WTRT2Ns1    , &
    WTRT2Ps1    =>   plt_biom%WTRT2Ps1    , &
    RTWT1s1     =>   plt_biom%RTWT1s1     , &
    RTWT1Ns1    =>   plt_biom%RTWT1Ns1    , &
    RTWT1Ps1    =>   plt_biom%RTWT1Ps1    , &
    WTRT1s1     =>   plt_biom%WTRT1s1     , &
    WTRT1Ns1    =>   plt_biom%WTRT1Ns1    , &
    WTRT1Ps1    =>   plt_biom%WTRT1Ps1    , &
    CPOOLRs1    =>   plt_biom%CPOOLRs1    , &
    ZPOOLRs1    =>   plt_biom%ZPOOLRs1    , &
    PPOOLRs1    =>   plt_biom%PPOOLRs1    , &
    WTRTDs1     =>   plt_biom%WTRTDs1     , &
    ZPOOLNs1    =>   plt_biom%ZPOOLNs1    , &
    WTRTLs1     =>   plt_biom%WTRTLs1     , &
    WSRTLs1     =>   plt_biom%WSRTLs1     , &
    CPOOLNs1    =>   plt_biom%CPOOLNs1    , &
    WTNDLs1     =>   plt_biom%WTNDLs1     , &
    WTNDLNs1    =>   plt_biom%WTNDLNs1    , &
    WTNDLPs1    =>   plt_biom%WTNDLPs1    , &
    PPOOLNs1    =>   plt_biom%PPOOLNs1    , &
    FWODRs1     =>   plt_allom%FWODRs1    , &
    FWODRNs1    =>   plt_allom%FWODRNs1   , &
    FWODRPs1    =>   plt_allom%FWODRPs1   , &
    IDTHRs1     =>   plt_pheno%IDTHRs1    , &
    CFOPCs1     =>   plt_soilchem%CFOPCs1 , &
    CFOPNs1     =>   plt_soilchem%CFOPNs1 , &
    CFOPPs1     =>   plt_soilchem%CFOPPs1 , &
    CO2Ps1      =>   plt_rbgc%CO2Ps1  , &
    OXYPs1      =>   plt_rbgc%OXYPs1  , &
    CO2As1      =>   plt_rbgc%CO2As1  , &
    OXYAs1      =>   plt_rbgc%OXYAs1  , &
    CH4Ps1      =>   plt_rbgc%CH4Ps1  , &
    CH4As1      =>   plt_rbgc%CH4As1  , &
    H2GPs1      =>   plt_rbgc%H2GPs1  , &
    H2GAs1      =>   plt_rbgc%H2GAs1  , &
    Z2OPs1      =>   plt_rbgc%Z2OPs1      , &
    Z2OAs1      =>   plt_rbgc%Z2OAs1      , &
    ZH3Ps1      =>   plt_rbgc%ZH3Ps1      , &
    ZH3As1      =>   plt_rbgc%ZH3As1      , &
    NJs1        =>   plt_site%NJs1        , &
    NUs1        =>   plt_site%NUs1        , &
    RH2GZs1     =>   plt_bgcr%RH2GZs1     , &
    RNH3Zs1     =>   plt_bgcr%RNH3Zs1     , &
    RN2OZs1     =>   plt_bgcr%RN2OZs1     , &
    RCH4Zs1     =>   plt_bgcr%RCH4Zs1     , &
    ROXYZs1     =>   plt_bgcr%ROXYZs1     , &
    RCO2Zs1     =>   plt_bgcr%RCO2Zs1     , &
    CSNCs1      =>   plt_bgcr%CSNCs1      , &
    ZSNCs1      =>   plt_bgcr%ZSNCs1      , &
    PSNCs1      =>   plt_bgcr%PSNCs1      , &
    RTDNPs1     =>   plt_morph%RTDNPs1    , &
    RTVLWs1     =>   plt_morph%RTVLWs1    , &
    RTARPs1     =>   plt_morph%RTARPs1    , &
    RTVLPs1     =>   plt_morph%RTVLPs1    , &
    RTNLs1      =>   plt_morph%RTNLs1     , &
    RTN1s1      =>   plt_morph%RTN1s1     , &
    RTDP1s1     =>   plt_morph%RTDP1s1    , &
    RTLGAs1     =>   plt_morph%RTLGAs1    , &
    RRAD1s1     =>   plt_morph%RRAD1s1    , &
    RRAD2s1     =>   plt_morph%RRAD2s1    , &
    RRAD1Ms1    =>   plt_morph%RRAD1Ms1   , &
    RRAD2Ms1    =>   plt_morph%RRAD2Ms1   , &
    RTLGPs1     =>   plt_morph%RTLGPs1    , &
    RTLG1s1     =>   plt_morph%RTLG1s1    , &
    RTLG2s1     =>   plt_morph%RTLG2s1    , &
    INTYPs1     =>   plt_morph%INTYPs1    , &
    RTN2s1      =>   plt_morph%RTN2s1     , &
    MYs1        =>   plt_morph%MYs1       , &
    NIXs1       =>   plt_morph%NIXs1      , &
    NINRs1      =>   plt_morph%NINRs1     , &
    SDPTHs1     =>   plt_morph%SDPTHs1    , &
    NGs1        =>   plt_morph%NGs1       , &
    NRTs1       =>   plt_morph%NRTs1        &
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

  IF(IDTHRs1(NZ).EQ.1)THEN
    DO 8900 N=1,MYs1(NZ)
      DO 8895 L=NUs1,NJs1
        DO 6410 M=1,4
          CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(0,M,NZ)*CPOOLRs1(N,L,NZ)
          ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(0,M,NZ)*ZPOOLRs1(N,L,NZ)
          PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(0,M,NZ)*PPOOLRs1(N,L,NZ)
          DO  NR=1,NRTs1(NZ)
            CSNCs1(M,0,L,NZ)=CSNCs1(M,0,L,NZ)+CFOPCs1(5,M,NZ) &
              *(WTRT1s1(N,L,NR,NZ)+WTRT2s1(N,L,NR,NZ))*FWODRs1(0)
            ZSNCs1(M,0,L,NZ)=ZSNCs1(M,0,L,NZ)+CFOPNs1(5,M,NZ) &
              *(WTRT1Ns1(N,L,NR,NZ)+WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(0)
            PSNCs1(M,0,L,NZ)=PSNCs1(M,0,L,NZ)+CFOPPs1(5,M,NZ) &
              *(WTRT1Ps1(N,L,NR,NZ)+WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(0)
            CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ) &
              *(WTRT1s1(N,L,NR,NZ)+WTRT2s1(N,L,NR,NZ))*FWODRs1(1)
            ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ) &
              *(WTRT1Ns1(N,L,NR,NZ)+WTRT2Ns1(N,L,NR,NZ))*FWODRNs1(1)
            PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ) &
              *(WTRT1Ps1(N,L,NR,NZ)+WTRT2Ps1(N,L,NR,NZ))*FWODRPs1(1)
          enddo
6410    CONTINUE
!
!     RELEASE GAS CONTENTS OF DEAD ROOTS
!
        RCO2Zs1(NZ)=RCO2Zs1(NZ)-CO2As1(N,L,NZ)-CO2Ps1(N,L,NZ)
        ROXYZs1(NZ)=ROXYZs1(NZ)-OXYAs1(N,L,NZ)-OXYPs1(N,L,NZ)
        RCH4Zs1(NZ)=RCH4Zs1(NZ)-CH4As1(N,L,NZ)-CH4Ps1(N,L,NZ)
        RN2OZs1(NZ)=RN2OZs1(NZ)-Z2OAs1(N,L,NZ)-Z2OPs1(N,L,NZ)
        RNH3Zs1(NZ)=RNH3Zs1(NZ)-ZH3As1(N,L,NZ)-ZH3Ps1(N,L,NZ)
        RH2GZs1(NZ)=RH2GZs1(NZ)-H2GAs1(N,L,NZ)-H2GPs1(N,L,NZ)
        CO2As1(N,L,NZ)=0._r8
        OXYAs1(N,L,NZ)=0._r8
        CH4As1(N,L,NZ)=0._r8
        Z2OAs1(N,L,NZ)=0._r8
        ZH3As1(N,L,NZ)=0._r8
        H2GAs1(N,L,NZ)=0._r8
        CO2Ps1(N,L,NZ)=0._r8
        OXYPs1(N,L,NZ)=0._r8
        CH4Ps1(N,L,NZ)=0._r8
        Z2OPs1(N,L,NZ)=0._r8
        ZH3Ps1(N,L,NZ)=0._r8
        H2GPs1(N,L,NZ)=0._r8
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
        DO 8870 NR=1,NRTs1(NZ)
          WTRT1s1(N,L,NR,NZ)=0._r8
          WTRT1Ns1(N,L,NR,NZ)=0._r8
          WTRT1Ps1(N,L,NR,NZ)=0._r8
          WTRT2s1(N,L,NR,NZ)=0._r8
          WTRT2Ns1(N,L,NR,NZ)=0._r8
          WTRT2Ps1(N,L,NR,NZ)=0._r8
          RTWT1s1(N,NR,NZ)=0._r8
          RTWT1Ns1(N,NR,NZ)=0._r8
          RTWT1Ps1(N,NR,NZ)=0._r8
          RTLG1s1(N,L,NR,NZ)=0._r8
          RTLG2s1(N,L,NR,NZ)=0._r8
          RTN2s1(N,L,NR,NZ)=0._r8
8870    CONTINUE
        CPOOLRs1(N,L,NZ)=0._r8
        ZPOOLRs1(N,L,NZ)=0._r8
        PPOOLRs1(N,L,NZ)=0._r8
        WTRTLs1(N,L,NZ)=0._r8
        WTRTDs1(N,L,NZ)=0._r8
        WSRTLs1(N,L,NZ)=0._r8
        RTN1s1(N,L,NZ)=0._r8
        RTNLs1(N,L,NZ)=0._r8
        RTLGPs1(N,L,NZ)=0._r8
        RTDNPs1(N,L,NZ)=0._r8
        RTVLPs1(N,L,NZ)=0._r8
        RTVLWs1(N,L,NZ)=0._r8
        RRAD1s1(N,L,NZ)=RRAD1Ms1(N,NZ)
        RRAD2s1(N,L,NZ)=RRAD2Ms1(N,NZ)
        RTARPs1(N,L,NZ)=0._r8
        RTLGAs1(N,L,NZ)=RTLGAX
!
!     LITTERFALL AND STATE VARIABLES FROM DEAD NODULES
!
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
        IF(INTYPs1(NZ).NE.0.AND.N.EQ.1)THEN
          DO 6420 M=1,4
            CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ) &
              *WTNDLs1(L,NZ)+CFOPCs1(0,M,NZ)*CPOOLNs1(L,NZ)
            ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ) &
              *WTNDLNs1(L,NZ)+CFOPNs1(0,M,NZ)*ZPOOLNs1(L,NZ)
            PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ) &
              *WTNDLPs1(L,NZ)+CFOPPs1(0,M,NZ)*PPOOLNs1(L,NZ)
6420      CONTINUE
          WTNDLs1(L,NZ)=0._r8
          WTNDLNs1(L,NZ)=0._r8
          WTNDLPs1(L,NZ)=0._r8
          CPOOLNs1(L,NZ)=0._r8
          ZPOOLNs1(L,NZ)=0._r8
          PPOOLNs1(L,NZ)=0._r8
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
    DO 8795 NR=1,NRTs1(NZ)
      NINRs1(NR,NZ)=NGs1(NZ)
      DO 8790 N=1,MYs1(NZ)
        RTDP1s1(N,NR,NZ)=SDPTHs1(NZ)
        RTWT1s1(N,NR,NZ)=0._r8
        RTWT1Ns1(N,NR,NZ)=0._r8
        RTWT1Ps1(N,NR,NZ)=0._r8
8790  CONTINUE
8795  CONTINUE
    NIXs1(NZ)=NGs1(NZ)
    NRTs1(NZ)=0
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
    IHVSTs1     =>  plt_distb%IHVSTs1     , &
    FWODBs1     =>  plt_allom%FWODBs1     , &
    FWODSPs1    =>  plt_allom%FWODSPs1    , &
    FWODSNs1    =>  plt_allom%FWODSNs1    , &
    FWODLNs1    =>  plt_allom%FWODLNs1    , &
    FWODLPs1    =>  plt_allom%FWODLPs1    , &
    WTEABNs1    =>  plt_biom%WTEABNs1     , &
    CPOLNBs1    =>  plt_biom%CPOLNBs1     , &
    ZPOLNBs1    =>  plt_biom%ZPOLNBs1     , &
    PPOLNBs1    =>  plt_biom%PPOLNBs1     , &
    WTSTKBs1    =>  plt_biom%WTSTKBs1     , &
    WTHSKBs1    =>  plt_biom%WTHSKBs1     , &
    WTSHEBs1    =>  plt_biom%WTHSKBs1     , &
    WTHSBNs1    =>  plt_biom%WTHSBNs1     , &
    WTHSBPs1    =>  plt_biom%WTHSBPs1     , &
    WTSHBPs1    =>  plt_biom%WTSHBPs1     , &
    WTNDBs1     =>  plt_biom%WTNDBs1      , &
    WTLFBs1     =>  plt_biom%WTLFBs1      , &
    ZEROPs1     =>  plt_biom%ZEROPs1      , &
    WTLFBNs1    =>  plt_biom%WTLFBNs1     , &
    WTLFBPs1    =>  plt_biom%WTLFBPs1     , &
    WTNDBNs1    =>  plt_biom%WTNDBNs1     , &
    WTNDBPs1    =>  plt_biom%WTNDBPs1     , &
    WTSHBNs1    =>  plt_biom%WTSHBNs1     , &
    PPOOLs1     =>  plt_biom%PPOOLs1      , &
    WTEARBs1    =>  plt_biom%WTEARBs1     , &
    WTEABPs1    =>  plt_biom%WTEABPs1     , &
    WTSTBNs1    =>  plt_biom%WTSTBNs1     , &
    WTSTBPs1    =>  plt_biom%WTSTBPs1     , &
    ZPOOLs1     =>  plt_biom%ZPOOLs1      , &
    CPOOLs1     =>  plt_biom%CPOOLs1      , &
    WTRVPs1     =>  plt_biom%WTRVPs1      , &
    WTGRBPs1    =>  plt_biom%WTGRBPs1     , &
    WTGRBs1     =>  plt_biom%WTGRBs1      , &
    WTRSVBs1    =>  plt_biom%WTRSVBs1     , &
    WTRSBNs1    =>  plt_biom%WTRSBNs1     , &
    WTRSBPs1    =>  plt_biom%WTRSBPs1     , &
    WTSTDGs1    =>  plt_biom%WTSTDGs1     , &
    WTSTDNs1    =>  plt_biom%WTSTDNs1     , &
    WTSTDPs1    =>  plt_biom%WTSTDPs1     , &
    WTGRBNs1    =>  plt_biom%WTGRBNs1     , &
    WTRVCs1     =>  plt_biom%WTRVCs1      , &
    WTRVNs1     =>  plt_biom%WTRVNs1      , &
    CFOPCs1     =>  plt_soilchem%CFOPCs1  , &
    CFOPNs1     =>  plt_soilchem%CFOPNs1  , &
    CFOPPs1     =>  plt_soilchem%CFOPPs1  , &
    IDTHBs1     =>  plt_pheno%IDTHBs1     , &
    GROUPs1     =>  plt_pheno%GROUPs1     , &
    VSTGXs1     =>  plt_pheno%VSTGXs1     , &
    KVSTGs1     =>  plt_pheno%KVSTGs1     , &
    TGSTGIs1    =>  plt_pheno%TGSTGIs1    , &
    TGSTGFs1    =>  plt_pheno%TGSTGFs1    , &
    VRNSs1      =>  plt_pheno%VRNSs1      , &
    VRNFs1      =>  plt_pheno%VRNFs1      , &
    VRNYs1      =>  plt_pheno%VRNYs1      , &
    VRNZs1      =>  plt_pheno%VRNZs1      , &
    ATRPs1      =>  plt_pheno%ATRPs1      , &
    FLG4s1      =>  plt_pheno%FLG4s1      , &
    IFLGAs1     =>  plt_pheno%IFLGAs1     , &
    IFLGEs1     =>  plt_pheno%IFLGEs1     , &
    IFLGRs1     =>  plt_pheno%IFLGRs1     , &
    IFLGQs1     =>  plt_pheno%IFLGQs1     , &
    GROUPIs1    =>  plt_pheno%GROUPIs1    , &
    IFLGFs1     =>  plt_pheno%IFLGFs1     , &
    IDAYs1      =>  plt_pheno%IDAYs1      , &
    IBTYPs1     =>  plt_pheno%IBTYPs1     , &
    IGTYPs1     =>  plt_pheno%IGTYPs1     , &
    IWTYPs1     =>  plt_pheno%IWTYPs1     , &
    ISTYPs1     =>  plt_pheno%ISTYPs1     , &
    CSNCs1      =>  plt_bgcr%CSNCs1       , &
    ZSNCs1      =>  plt_bgcr%ZSNCs1       , &
    PSNCs1      =>  plt_bgcr%PSNCs1       , &
    NBRs1       =>  plt_morph%NBRs1       , &
    PSTGIs1     =>  plt_morph%PSTGIs1     , &
    PSTGs1      =>  plt_morph%PSTGs1      , &
    NBTBs1      =>  plt_morph%NBTBs1      , &
    PSTGFs1     =>  plt_morph%PSTGFs1     , &
    XTLIs1      =>  plt_morph%XTLIs1      , &
    VSTGs1      =>  plt_morph%VSTGs1      , &
    KLEAFs1     =>  plt_morph%KLEAFs1     , &
    FDBKs1      =>  plt_photo%FDBKs1      , &
    FDBKXs1     =>  plt_photo%FDBKXs1       &
  )
  DO 8845 NB=1,NBRs1(NZ)
    IF(IDTHBs1(NB,NZ).EQ.1)THEN
      GROUPs1(NB,NZ)=GROUPIs1(NZ)
      PSTGs1(NB,NZ)=XTLIs1(NZ)
      PSTGIs1(NB,NZ)=PSTGs1(NB,NZ)
      PSTGFs1(NB,NZ)=0._r8
      VSTGs1(NB,NZ)=0._r8
      VSTGXs1(NB,NZ)=0._r8
      KLEAFs1(NB,NZ)=1
      KVSTGs1(NB,NZ)=1
      TGSTGIs1(NB,NZ)=0._r8
      TGSTGFs1(NB,NZ)=0._r8
      VRNSs1(NB,NZ)=0._r8
      VRNFs1(NB,NZ)=0._r8
      VRNYs1(NB,NZ)=0._r8
      VRNZs1(NB,NZ)=0._r8
      ATRPs1(NB,NZ)=0._r8
      FLG4s1(NB,NZ)=0._r8
      FDBKs1(NB,NZ)=1.0_r8
      FDBKXs1(NB,NZ)=1.0_r8
      IFLGAs1(NB,NZ)=0
      IFLGEs1(NB,NZ)=1
      IFLGFs1(NB,NZ)=0
      IFLGRs1(NB,NZ)=0
      IFLGQs1(NB,NZ)=0
      NBTBs1(NB,NZ)=0
      DO 8850 M=1,10
        IDAYs1(M,NB,NZ)=0
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
      DO 6405 M=1,4
        CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
          +CFOPCs1(0,M,NZ)*CPOLNBs1(NB,NZ) &
          +CFOPCs1(1,M,NZ)*(WTLFBs1(NB,NZ)*FWODBs1(1) &
          +WTNDBs1(NB,NZ)) &
          +CFOPCs1(2,M,NZ)*(WTSHEBs1(NB,NZ)*FWODBs1(1) &
          +WTHSKBs1(NB,NZ)+WTEARBs1(NB,NZ))
        CSNCs1(M,0,0,NZ)=CSNCs1(M,0,0,NZ) &
          +CFOPCs1(5,M,NZ)*(WTLFBs1(NB,NZ)*FWODBs1(0) &
          +WTSHEBs1(NB,NZ)*FWODBs1(0))
        ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
          +CFOPNs1(0,M,NZ)*ZPOLNBs1(NB,NZ) &
          +CFOPNs1(1,M,NZ)*(WTLFBNs1(NB,NZ)*FWODLNs1(1) &
          +WTNDBNs1(NB,NZ)) &
          +CFOPNs1(2,M,NZ)*(WTSHBNs1(NB,NZ)*FWODSNs1(1) &
          +WTHSBNs1(NB,NZ)+WTEABNs1(NB,NZ))
        ZSNCs1(M,0,0,NZ)=ZSNCs1(M,0,0,NZ) &
          +CFOPNs1(5,M,NZ)*(WTLFBNs1(NB,NZ)*FWODLNs1(0) &
          +WTSHBNs1(NB,NZ)*FWODSNs1(0))
        PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
          +CFOPPs1(0,M,NZ)*PPOLNBs1(NB,NZ) &
          +CFOPPs1(1,M,NZ)*(WTLFBPs1(NB,NZ)*FWODLPs1(1) &
          +WTNDBPs1(NB,NZ)) &
          +CFOPPs1(2,M,NZ)*(WTSHBPs1(NB,NZ)*FWODSPs1(1) &
          +WTHSBPs1(NB,NZ)+WTEABPs1(NB,NZ))
        PSNCs1(M,0,0,NZ)=PSNCs1(M,0,0,NZ) &
          +CFOPPs1(5,M,NZ)*(WTLFBPs1(NB,NZ)*FWODLPs1(0) &
          +WTSHBPs1(NB,NZ)*FWODSPs1(0))
        IF(ISTYPs1(NZ).EQ.0.AND.IWTYPs1(NZ).NE.0)THEN
          WTRVCs1(NZ)=WTRVCs1(NZ)+CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
          WTRVNs1(NZ)=WTRVNs1(NZ)+CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
          WTRVPs1(NZ)=WTRVPs1(NZ)+CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
        ELSE
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(2,M,NZ)*WTGRBs1(NB,NZ)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(2,M,NZ)*WTGRBNs1(NB,NZ)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(2,M,NZ)*WTGRBPs1(NB,NZ)
        ENDIF
        IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(3,M,NZ)*WTSTKBs1(NB,NZ)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(3,M,NZ)*WTSTBNs1(NB,NZ)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(3,M,NZ)*WTSTBPs1(NB,NZ)
        ELSE
          WTSTDGs1(M,NZ)=WTSTDGs1(M,NZ)+CFOPCs1(5,M,NZ)*WTSTKBs1(NB,NZ)
          WTSTDNs1(M,NZ)=WTSTDNs1(M,NZ)+CFOPNs1(5,M,NZ)*WTSTBNs1(NB,NZ)
          WTSTDPs1(M,NZ)=WTSTDPs1(M,NZ)+CFOPPs1(5,M,NZ)*WTSTBPs1(NB,NZ)
        ENDIF
6405  CONTINUE
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
      WTRVCs1(NZ)=WTRVCs1(NZ)+CPOOLs1(NB,NZ)+CPOOLK(NB,NZ)
      WTRVNs1(NZ)=WTRVNs1(NZ)+ZPOOLs1(NB,NZ)
      WTRVPs1(NZ)=WTRVPs1(NZ)+PPOOLs1(NB,NZ)
      IF(IHVSTs1(NZ).NE.4.AND.IHVSTs1(NZ).NE.6)THEN
        DO 6406 M=1,4
          CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ) &
            +CFOPCs1(0,M,NZ)*WTRSVBs1(NB,NZ)
          ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ) &
            +CFOPNs1(0,M,NZ)*WTRSBNs1(NB,NZ)
          PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ) &
            +CFOPPs1(0,M,NZ)*WTRSBPs1(NB,NZ)
6406    CONTINUE
      ELSE
        WTRVCs1(NZ)=WTRVCs1(NZ)+WTRSVBs1(NB,NZ)
        WTRVNs1(NZ)=WTRVNs1(NZ)+WTRSBNs1(NB,NZ)
        WTRVPs1(NZ)=WTRVPs1(NZ)+WTRSBPs1(NB,NZ)
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
    CPOOLs1  => plt_biom%CPOOLs1         , &
    ZPOOLs1  => plt_biom%ZPOOLs1         , &
    CPOLNBs1 => plt_biom%CPOLNBs1        , &
    ZPOLNBs1 => plt_biom%ZPOLNBs1        , &
    WTSHTBs1 => plt_biom%WTSHTBs1        , &
    WTLFBs1  => plt_biom%WTLFBs1         , &
    WTRSVBs1 => plt_biom%WTRSVBs1        , &
    WTSHEBs1 => plt_biom%WTSHEBs1        , &
    WTHSKBs1 => plt_biom%WTHSKBs1        , &
    WTEARBs1 => plt_biom%WTEARBs1        , &
    WTGRBs1  => plt_biom%WTGRBs1         , &
    WVSTKBs1 => plt_biom%WVSTKBs1        , &
    WTRSBNs1 => plt_biom%WTRSBNs1        , &
    WTHSBNs1 => plt_biom%WTHSBNs1        , &
    CPOOLRs1 => plt_biom%CPOOLRs1        , &
    ZPOOLRs1 => plt_biom%ZPOOLRs1        , &
    PPOOLRs1 => plt_biom%PPOOLRs1        , &
    WTEABNs1 => plt_biom%WTEABNs1        , &
    WTGRBNs1 => plt_biom%WTGRBNs1        , &
    WTSHTPs1 => plt_biom%WTSHTPs1        , &
    WTLFBPs1 => plt_biom%WTLFBPs1        , &
    WTNDBPs1 => plt_biom%WTNDBPs1        , &
    WTSHBPs1 => plt_biom%WTSHBPs1        , &
    WTSTBPs1 => plt_biom%WTSTBPs1        , &
    WTRSBPs1 => plt_biom%WTRSBPs1        , &
    WTHSBPs1 => plt_biom%WTHSBPs1        , &
    WTEABPs1 => plt_biom%WTEABPs1        , &
    WTGRBPs1 => plt_biom%WTGRBPs1        , &
    WTSTXBs1 => plt_biom%WTSTXBs1        , &
    WTSTXNs1 => plt_biom%WTSTXNs1        , &
    WTSTXPs1 => plt_biom%WTSTXPs1        , &
    WTSTBNs1 => plt_biom%WTSTBNs1        , &
    WTLSBs1  => plt_biom%WTLSBs1         , &
    WTSHBNs1 => plt_biom%WTSHBNs1        , &
    WTRVCs1  => plt_biom%WTRVCs1         , &
    WTRVNs1  => plt_biom%WTRVNs1         , &
    WTRVPs1  => plt_biom%WTRVPs1         , &
    WTLFBNs1 => plt_biom%WTLFBNs1        , &
    WTSHTNs1 => plt_biom%WTSHTNs1        , &
    WTRT1s1  => plt_biom%WTRT1s1         , &
    WTRT1Ns1 => plt_biom%WTRT1Ns1        , &
    WTRT1Ps1 => plt_biom%WTRT1Ps1        , &
    WTNDBNs1 => plt_biom%WTNDBNs1        , &
    WTRT2s1  => plt_biom%WTRT2s1         , &
    WTRT2Ns1 => plt_biom%WTRT2Ns1        , &
    WTRT2Ps1 => plt_biom%WTRT2Ps1        , &
    WTNDBs1  => plt_biom%WTNDBs1         , &
    WTSTKBs1 => plt_biom%WTSTKBs1        , &
    RTWT1s1  => plt_biom%RTWT1s1         , &
    RTWT1Ns1 => plt_biom%RTWT1Ns1        , &
    RTWT1Ps1 => plt_biom%RTWT1Ps1        , &
    PPOLNBs1 => plt_biom%PPOLNBs1        , &
    PPOOLs1  => plt_biom%PPOOLs1         , &
    NJs1     => plt_site%NJs1            , &
    NUs1     => plt_site%NUs1            , &
    IDTHs1   => plt_pheno%IDTHs1         , &
    RTLG2s1  => plt_morph%RTLG2s1        , &
    RTN2s1   => plt_morph%RTN2s1         , &
    RTLG1s1  => plt_morph%RTLG1s1        , &
    MYs1     => plt_morph%MYs1           , &
    NBRs1    => plt_morph%NBRs1          , &
    NRTs1    => plt_morph%NRTs1            &
  )
!     RESET BRANCH STATE VARIABLES
!
  DO 8835 NB=1,NBRs1(NZ)
    CPOOLs1(NB,NZ)=0._r8
    CPOOLK(NB,NZ)=0._r8
    ZPOOLs1(NB,NZ)=0._r8
    PPOOLs1(NB,NZ)=0._r8
    CPOLNBs1(NB,NZ)=0._r8
    ZPOLNBs1(NB,NZ)=0._r8
    PPOLNBs1(NB,NZ)=0._r8
    WTSHTBs1(NB,NZ)=0._r8
    WTLFBs1(NB,NZ)=0._r8
    WTNDBs1(NB,NZ)=0._r8
    WTSHEBs1(NB,NZ)=0._r8
    WTSTKBs1(NB,NZ)=0._r8
    WVSTKBs1(NB,NZ)=0._r8
    WTRSVBs1(NB,NZ)=0._r8
    WTHSKBs1(NB,NZ)=0._r8
    WTEARBs1(NB,NZ)=0._r8
    WTGRBs1(NB,NZ)=0._r8
    WTLSBs1(NB,NZ)=0._r8
    WTSHTNs1(NB,NZ)=0._r8
    WTLFBNs1(NB,NZ)=0._r8
    WTNDBNs1(NB,NZ)=0._r8
    WTSHBNs1(NB,NZ)=0._r8
    WTSTBNs1(NB,NZ)=0._r8
    WTRSBNs1(NB,NZ)=0._r8
    WTHSBNs1(NB,NZ)=0._r8
    WTEABNs1(NB,NZ)=0._r8
    WTGRBNs1(NB,NZ)=0._r8
    WTSHTPs1(NB,NZ)=0._r8
    WTLFBPs1(NB,NZ)=0._r8
    WTNDBPs1(NB,NZ)=0._r8
    WTSHBPs1(NB,NZ)=0._r8
    WTSTBPs1(NB,NZ)=0._r8
    WTRSBPs1(NB,NZ)=0._r8
    WTHSBPs1(NB,NZ)=0._r8
    WTEABPs1(NB,NZ)=0._r8
    WTGRBPs1(NB,NZ)=0._r8
    WTSTXBs1(NB,NZ)=0._r8
    WTSTXNs1(NB,NZ)=0._r8
    WTSTXPs1(NB,NZ)=0._r8
8835  CONTINUE
!
!     RESET ROOT STATE VARIABLES
!
  DO 6416 L=NUs1,NJs1
    DO  N=1,MYs1(NZ)
      CPOOLRs1(N,L,NZ)=0._r8
      ZPOOLRs1(N,L,NZ)=0._r8
      PPOOLRs1(N,L,NZ)=0._r8
      DO  NR=1,NRTs1(NZ)
        WTRT1s1(N,L,NR,NZ)=0._r8
        WTRT1Ns1(N,L,NR,NZ)=0._r8
        WTRT1Ps1(N,L,NR,NZ)=0._r8
        WTRT2s1(N,L,NR,NZ)=0._r8
        WTRT2Ns1(N,L,NR,NZ)=0._r8
        WTRT2Ps1(N,L,NR,NZ)=0._r8
        RTWT1s1(N,NR,NZ)=0._r8
        RTWT1Ns1(N,NR,NZ)=0._r8
        RTWT1Ps1(N,NR,NZ)=0._r8
        RTLG1s1(N,L,NR,NZ)=0._r8
        RTLG2s1(N,L,NR,NZ)=0._r8
        RTN2s1(N,L,NR,NZ)=0._r8
      enddo
    enddo
6416  CONTINUE
  WTRVCs1(NZ)=0._r8
  WTRVNs1(NZ)=0._r8
  WTRVPs1(NZ)=0._r8
  IDTHs1(NZ)=1
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
    CPOOLs1    => plt_biom%CPOOLs1      , &
    ZPOOLs1    => plt_biom%ZPOOLs1      , &
    PPOOLs1    => plt_biom%PPOOLs1      , &
    CPOLNBs1   => plt_biom%CPOLNBs1     , &
    PPOLNBs1   => plt_biom%PPOLNBs1     , &
    WTSHTBs1   => plt_biom%WTSHTBs1     , &
    WTLFBs1    => plt_biom%WTLFBs1      , &
    WTNDBs1    => plt_biom%WTNDBs1      , &
    WTSTKBs1   => plt_biom%WTSTKBs1     , &
    WGLFs1     => plt_biom%WGLFs1       , &
    WTRSVBs1   => plt_biom%WTRSVBs1     , &
    WTRSBPs1   => plt_biom%WTRSBPs1     , &
    WTHSBPs1   => plt_biom%WTHSBPs1     , &
    WTSTXNs1   => plt_biom%WTSTXNs1     , &
    WTSTXPs1   => plt_biom%WTSTXPs1     , &
    WTSTXBs1   => plt_biom%WTSTXBs1     , &
    WTSHBPs1   => plt_biom%WTSHBPs1     , &
    WTEABPs1   => plt_biom%WTEABPs1     , &
    WTGRBPs1   => plt_biom%WTGRBPs1     , &
    WTLFBPs1   => plt_biom%WTLFBPs1     , &
    WTNDBPs1   => plt_biom%WTNDBPs1     , &
    WTSHTNs1   => plt_biom%WTSHTNs1     , &
    WTSHBNs1   => plt_biom%WTSHBNs1     , &
    WTSTBPs1   => plt_biom%WTSTBPs1     , &
    WTNDBNs1   => plt_biom%WTNDBNs1     , &
    WTSTBNs1   => plt_biom%WTSTBNs1     , &
    WTSHTPs1   => plt_biom%WTSHTPs1     , &
    WTEABNs1   => plt_biom%WTEABNs1     , &
    WTGRBNs1   => plt_biom%WTGRBNs1     , &
    WTHSBNs1   => plt_biom%WTHSBNs1     , &
    WTRSBNs1   => plt_biom%WTRSBNs1     , &
    WTLFBNs1   => plt_biom%WTLFBNs1     , &
    WVSTKBs1   => plt_biom%WVSTKBs1     , &
    WTGRBs1    => plt_biom%WTGRBs1      , &
    WTLSBs1    => plt_biom%WTLSBs1      , &
    WTHSKBs1   => plt_biom%WTHSKBs1     , &
    WTEARBs1   => plt_biom%WTEARBs1     , &
    WSSHEs1    => plt_biom%WSSHEs1      , &
    WGLFNs1    => plt_biom%WGLFNs1      , &
    WSLFs1     => plt_biom%WSLFs1       , &
    WGLFPs1    => plt_biom%WGLFPs1      , &
    WGLFLPs1   => plt_biom%WGLFLPs1     , &
    WGNODEs1   => plt_biom%WGNODEs1     , &
    WGNODNs1   => plt_biom%WGNODNs1     , &
    WGNODPs1   => plt_biom%WGNODPs1     , &
    WGSHPs1    => plt_biom%WGSHPs1      , &
    WGSHNs1    => plt_biom%WGSHNs1      , &
    WGSHEs1    => plt_biom%WGSHEs1      , &
    WGLFLNs1   => plt_biom%WGLFLNs1     , &
    WGLFLs1    => plt_biom%WGLFLs1      , &
    WGLFVs1    => plt_biom%WGLFVs1      , &
    WTSHEBs1   => plt_biom%WTSHEBs1     , &
    ZPOLNBs1   => plt_biom%ZPOLNBs1     , &
    CPOOL3s1   => plt_photo%CPOOL3s1    , &
    CPOOL4s1   => plt_photo%CPOOL4s1    , &
    CO2Bs1     => plt_photo%CO2Bs1      , &
    HCOBs1     => plt_photo%HCOBs1      , &
    GRWTBs1    => plt_allom%GRWTBs1     , &
    GRNXBs1    => plt_morph%GRNXBs1     , &
    GRNOBs1    => plt_morph%GRNOBs1     , &
    ARLFBs1    => plt_morph%ARLFBs1     , &
    ARSTKs1    => plt_morph%ARSTKs1     , &
    NBRs1      => plt_morph%NBRs1       , &
    ARLF1s1    => plt_morph%ARLF1s1     , &
    HTSHEs1    => plt_morph%HTSHEs1     , &
    HTNODXs1   => plt_morph%HTNODXs1    , &
    SURFs1     => plt_morph%SURFs1      , &
    HTNODEs1   => plt_morph%HTNODEs1    , &
    ARLFVs1    => plt_morph%ARLFVs1     , &
    SURFBs1    => plt_morph%SURFBs1     , &
    ARLFLs1    => plt_morph%ARLFLs1       &
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
  CPOOLs1(NB,NZ)=0._r8
  CPOOLK(NB,NZ)=0._r8
  ZPOOLs1(NB,NZ)=0._r8
  PPOOLs1(NB,NZ)=0._r8
  CPOLNBs1(NB,NZ)=0._r8
  ZPOLNBs1(NB,NZ)=0._r8
  PPOLNBs1(NB,NZ)=0._r8
  WTSHTBs1(NB,NZ)=0._r8
  WTLFBs1(NB,NZ)=0._r8
  WTNDBs1(NB,NZ)=0._r8
  WTSHEBs1(NB,NZ)=0._r8
  WTSTKBs1(NB,NZ)=0._r8
  WVSTKBs1(NB,NZ)=0._r8
  WTRSVBs1(NB,NZ)=0._r8
  WTHSKBs1(NB,NZ)=0._r8
  WTEARBs1(NB,NZ)=0._r8
  WTGRBs1(NB,NZ)=0._r8
  WTLSBs1(NB,NZ)=0._r8
  WTSHTNs1(NB,NZ)=0._r8
  WTLFBNs1(NB,NZ)=0._r8
  WTNDBNs1(NB,NZ)=0._r8
  WTSHBNs1(NB,NZ)=0._r8
  WTSTBNs1(NB,NZ)=0._r8
  WTRSBNs1(NB,NZ)=0._r8
  WTHSBNs1(NB,NZ)=0._r8
  WTEABNs1(NB,NZ)=0._r8
  WTGRBNs1(NB,NZ)=0._r8
  WTSHTPs1(NB,NZ)=0._r8
  WTLFBPs1(NB,NZ)=0._r8
  WTNDBPs1(NB,NZ)=0._r8
  WTSHBPs1(NB,NZ)=0._r8
  WTSTBPs1(NB,NZ)=0._r8
  WTRSBPs1(NB,NZ)=0._r8
  WTHSBPs1(NB,NZ)=0._r8
  WTEABPs1(NB,NZ)=0._r8
  WTGRBPs1(NB,NZ)=0._r8
  GRNXBs1(NB,NZ)=0._r8
  GRNOBs1(NB,NZ)=0._r8
  GRWTBs1(NB,NZ)=0._r8
  ARLFBs1(NB,NZ)=0._r8
  WTSTXBs1(NB,NZ)=0._r8
  WTSTXNs1(NB,NZ)=0._r8
  WTSTXPs1(NB,NZ)=0._r8
  DO 8855 K=0,JNODS1
    IF(K.NE.0)THEN
      CPOOL3s1(K,NB,NZ)=0._r8
      CPOOL4s1(K,NB,NZ)=0._r8
      CO2Bs1(K,NB,NZ)=0._r8
      HCOBs1(K,NB,NZ)=0._r8
    ENDIF
    ARLF1s1(K,NB,NZ)=0._r8
    HTNODEs1(K,NB,NZ)=0._r8
    HTNODXs1(K,NB,NZ)=0._r8
    HTSHEs1(K,NB,NZ)=0._r8
    WGLFs1(K,NB,NZ)=0._r8
    WSLFs1(K,NB,NZ)=0._r8
    WGLFNs1(K,NB,NZ)=0._r8
    WGLFPs1(K,NB,NZ)=0._r8
    WGSHEs1(K,NB,NZ)=0._r8
    WSSHEs1(K,NB,NZ)=0._r8
    WGSHNs1(K,NB,NZ)=0._r8
    WGSHPs1(K,NB,NZ)=0._r8
    WGNODEs1(K,NB,NZ)=0._r8
    WGNODNs1(K,NB,NZ)=0._r8
    WGNODPs1(K,NB,NZ)=0._r8
    DO 8865 L=1,JC1
      ARLFVs1(L,NZ)=ARLFVs1(L,NZ)-ARLFLs1(L,K,NB,NZ)
      WGLFVs1(L,NZ)=WGLFVs1(L,NZ)-WGLFLs1(L,K,NB,NZ)
      ARLFLs1(L,K,NB,NZ)=0._r8
      WGLFLs1(L,K,NB,NZ)=0._r8
      WGLFLNs1(L,K,NB,NZ)=0._r8
      WGLFLPs1(L,K,NB,NZ)=0._r8
      IF(K.NE.0)THEN
        DO 8860 N=1,JLI1
          SURFs1(N,L,K,NB,NZ)=0._r8
8860    CONTINUE
      ENDIF
8865  CONTINUE
8855  CONTINUE
  DO 8875 L=1,JC1
    ARSTKs1(L,NB,NZ)=0._r8
    DO  N=1,JLI1
      SURFBs1(N,L,NB,NZ)=0._r8
    enddo
8875  CONTINUE
  end associate
  end subroutine ResetDeadRootStates


end module LitterFallMod
