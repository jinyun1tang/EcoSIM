module StartqsMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use EcoSIMConfig
  use PlantAPIData
  implicit none

  private

  public :: startqs
  contains

  SUBROUTINE startqs(NZ1Q,NZ2Q)
!
!     THIS SUBROUTINE INITIALIZES ALL PLANT VARIABLES
!
  implicit none
  integer, intent(in) :: NZ1Q,NZ2Q

  integer :: K,L,M,NZ,NZ2X
!     begin_execution

  associate(                            &
    PPs1      => plt_site%PPs1        , &
    NUs1      => plt_site%NUs1        , &
    NPs1      => plt_site%NPs1        , &
    NLs1      => plt_site%NLs1        , &
    ZEROs1    => plt_site%ZEROs1      , &
    AREA3s1   => plt_site%AREA3s1     , &
    ZEROQs1   => plt_rbgc%ZEROQs1     , &
    ZEROPs1   => plt_biom%ZEROPs1     , &
    ZEROLs1   => plt_biom%ZEROLs1     , &
    IFLGCs1   => plt_pheno%IFLGCs1      &
  )
!
!     INITIALIZE SHOOT GROWTH VARIABLES
!
!     IFLGC=PFT flag:0=not active,1=active
!     IYR0,IDAY0,IYRH,IDAYH=year,day of planting,arvesting
!     PPI,PPX=initial,current population (m-2)
!     CF,CFI=current,initial clumping factor
!     RSMH=cuticular resistance to water (h m-1)
!     RCMX=cuticular resistance to CO2 (s m-1)
!     CNWS,CPWS=protein:N,protein:P ratios
!     CWSRT=maximum root protein concentration (g g-1)
!     O2I=intercellular O2 concentration in C3,C4 PFT (umol mol-1)
!

      NZ2X=MIN(NZ2Q,NPs1)
      DO 9985 NZ=NZ1Q,NZ2X
        IF(IFLGCs1(NZ).EQ.0)THEN

          call InitShootGrowth(NZ)

          call PlantLitterFraction(NZ)

          call PFTThermalAcclimation(NZ)

          call InitDimensionsandUptake(NZ)

          call InitPlantPhenoMorphoBio(NZ)

          call InitMassBalance(NZ)

          call InitPlantHeatandWater(NZ)

          call InitRootMychorMorphoBio(NZ)

          call InitSeedMorphoBio(NZ)
        ENDIF
        ZEROPs1(NZ)=ZEROs1*PPs1(NZ)
        ZEROQs1(NZ)=ZEROs1*PPs1(NZ)/AREA3s1(NUs1)
        ZEROLs1(NZ)=ZEROs1*PPs1(NZ)*1.0E+06_r8
9985  CONTINUE
!
!     FILL OUT UNUSED ARRAYS
!
      DO 9986 NZ=NPs1+1,5
        plt_bgcr%TCSN0s1(NZ)=0._r8
        plt_bgcr%TZSN0s1(NZ)=0._r8
        plt_bgcr%TPSN0s1(NZ)=0._r8
        plt_bgcr%TCSNCs1(NZ)=0._r8
        plt_bgcr%TZSNCs1(NZ)=0._r8
        plt_bgcr%TPSNCs1(NZ)=0._r8
        plt_biom%WTSTGs1(NZ)=0._r8
        plt_biom%WTSTGNs1(NZ)=0._r8
        plt_biom%WTSTGPs1(NZ)=0._r8
        DO 6401 L=1,NLs1
          DO  K=0,1
            DO  M=1,jsken
              plt_bgcr%CSNCs1(M,K,L,NZ)=0._r8
              plt_bgcr%ZSNCs1(M,K,L,NZ)=0._r8
              plt_bgcr%PSNCs1(M,K,L,NZ)=0._r8
            enddo
          enddo
6401    CONTINUE
9986  CONTINUE
  RETURN
  end associate
  END subroutine startqs
!------------------------------------------------------------------------------------------

  subroutine InitShootGrowth(NZ)

  implicit none
  integer, intent(in) :: NZ
  associate(                         &
    IDAYXs1  =>  plt_distb%IDAYXs1 , &
    IYRYs1   =>  plt_distb%IYRYs1  , &
    IYRXs1   =>  plt_distb%IYRXs1  , &
    IYR0s1   =>  plt_distb%IYR0s1  , &
    IYRHs1   =>  plt_distb%IYRHs1  , &
    IDAYHs1  =>  plt_distb%IDAYHs1 , &
    IDAY0s1  =>  plt_distb%IDAY0s1 , &
    IDAYYs1  =>  plt_distb%IDAYYs1 , &
    RSMXs1   =>  plt_photo%RSMXs1  , &
    PPIs1    =>  plt_site%PPIs1    , &
    PPXs1    =>  plt_site%PPXs1    , &
    PPZs1    =>  plt_site%PPZs1    , &
    CWSRTs1  =>  plt_allom%CWSRTs1 , &
    CNWSs1   =>  plt_allom%CNWSs1  , &
    CPWSs1   =>  plt_allom%CPWSs1  , &
    CNRTs1   =>  plt_allom%CNRTs1  , &
    CPRTs1   =>  plt_allom%CPRTs1  , &
    O2Is1    =>  plt_photo%O2Is1   , &
    RCMXs1   =>  plt_photo%RCMXs1  , &
    ICTYPs1  =>  plt_photo%ICTYPs1 , &
    RSMHs1   =>  plt_photo%RSMHs1  , &
    CFIs1    =>  plt_morph%CFIs1   , &
    CFs1     =>  plt_morph%CFs1    , &
    NRTs1    =>  plt_morph%NRTs1     &
  )
  IYR0s1(NZ)=IYRXs1(NZ)
  IDAY0s1(NZ)=IDAYXs1(NZ)
  IYRHs1(NZ)=IYRYs1(NZ)
  IDAYHs1(NZ)=IDAYYs1(NZ)
  PPIs1(NZ)=PPZs1(NZ)
  PPXs1(NZ)=PPIs1(NZ)
  CFs1(NZ)=CFIs1(NZ)

  RSMHs1(NZ)=RSMXs1(NZ)/3600.0_r8
  RCMXs1(NZ)=RSMXs1(NZ)*1.56_r8
  CNWSs1(NZ)=2.5_r8
  CPWSs1(NZ)=25.0_r8
  CWSRTs1(NZ)=AMIN1(CNRTs1(NZ)*CNWSs1(NZ),CPRTs1(NZ)*CPWSs1(NZ))
  IF(ICTYPs1(NZ).EQ.3)THEN
    O2Is1(NZ)=2.10E+05_r8
  ELSE
    O2Is1(NZ)=3.96E+05_r8
  ENDIF
  end associate
  end subroutine InitShootGrowth
!------------------------------------------------------------------------------------------

  subroutine PlantLitterFraction(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: N,M
  real(r8) :: CNOPC(4),CPOPC(4)
  REAL(R8) :: CNOPCT,CPOPCT

  associate(                             &
    XRLAs1   =>  plt_pheno%XRLAs1      , &
    IBTYPs1  =>  plt_pheno%IBTYPs1     , &
    IGTYPs1  =>  plt_pheno%IGTYPs1     , &
    GROUPIs1 =>  plt_pheno%GROUPIs1    , &
    CFOPCs1  =>  plt_soilchem%CFOPCs1  , &
    CFOPNs1  =>  plt_soilchem%CFOPNs1  , &
    CFOPPs1  =>  plt_soilchem%CFOPPs1  , &
    FNODs1   =>  plt_allom%FNODs1      , &
    INTYPs1  =>  plt_morph%INTYPs1     , &
    NNODs1   =>  plt_morph%NNODs1        &
  )
!
!     FRACTIONS OF PLANT LITTER ALLOCATED TO KINETIC COMPONENTS
!     PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
!
!     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
!     foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), coarse woody (5,*)
!
!     NONSTRUCTURAL
!
  CFOPCs1(0,1,NZ)=0.0_r8
  CFOPCs1(0,2,NZ)=0.67_r8
  CFOPCs1(0,3,NZ)=0.33_r8
  CFOPCs1(0,4,NZ)=0.0_r8
!
!     NON-VASCULAR (E.G. MOSSES)
!
  IF(IGTYPs1(NZ).EQ.0)THEN
    CFOPCs1(1,1,NZ)=0.07_r8
    CFOPCs1(1,2,NZ)=0.25_r8
    CFOPCs1(1,3,NZ)=0.30_r8
    CFOPCs1(1,4,NZ)=0.38_r8
    CFOPCs1(2,1,NZ)=0.07_r8
    CFOPCs1(2,2,NZ)=0.25_r8
    CFOPCs1(2,3,NZ)=0.30_r8
    CFOPCs1(2,4,NZ)=0.38_r8
!
!     LEGUMES
!
  ELSEIF(INTYPs1(NZ).NE.0)THEN
    CFOPCs1(1,1,NZ)=0.16_r8
    CFOPCs1(1,2,NZ)=0.38_r8
    CFOPCs1(1,3,NZ)=0.34_r8
    CFOPCs1(1,4,NZ)=0.12_r8
    CFOPCs1(2,1,NZ)=0.07_r8
    CFOPCs1(2,2,NZ)=0.41_r8
    CFOPCs1(2,3,NZ)=0.37_r8
    CFOPCs1(2,4,NZ)=0.15_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
    CFOPCs1(1,1,NZ)=0.08_r8
    CFOPCs1(1,2,NZ)=0.41_r8
    CFOPCs1(1,3,NZ)=0.36_r8
    CFOPCs1(1,4,NZ)=0.15_r8
    CFOPCs1(2,1,NZ)=0.07_r8
    CFOPCs1(2,2,NZ)=0.41_r8
    CFOPCs1(2,3,NZ)=0.36_r8
    CFOPCs1(2,4,NZ)=0.16_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(IBTYPs1(NZ).EQ.1.OR.IBTYPs1(NZ).EQ.3)THEN
    CFOPCs1(1,1,NZ)=0.07_r8
    CFOPCs1(1,2,NZ)=0.34_r8
    CFOPCs1(1,3,NZ)=0.36_r8
    CFOPCs1(1,4,NZ)=0.23_r8
    CFOPCs1(2,1,NZ)=0.0_r8
    CFOPCs1(2,2,NZ)=0.045_r8
    CFOPCs1(2,3,NZ)=0.660_r8
    CFOPCs1(2,4,NZ)=0.295_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPCs1(1,1,NZ)=0.07_r8
    CFOPCs1(1,2,NZ)=0.25_r8
    CFOPCs1(1,3,NZ)=0.38_r8
    CFOPCs1(1,4,NZ)=0.30_r8
    CFOPCs1(2,1,NZ)=0.0_r8
    CFOPCs1(2,2,NZ)=0.045_r8
    CFOPCs1(2,3,NZ)=0.660_r8
    CFOPCs1(2,4,NZ)=0.295_r8
  ENDIF
!
!     FRACTIONS OF WOODY LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     NON-VASCULAR
!
  IF(IGTYPs1(NZ).EQ.0)THEN
    CFOPCs1(3,1,NZ)=0.07_r8
    CFOPCs1(3,2,NZ)=0.25_r8
    CFOPCs1(3,3,NZ)=0.30_r8
    CFOPCs1(3,4,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
    CFOPCs1(3,1,NZ)=0.03_r8
    CFOPCs1(3,2,NZ)=0.25_r8
    CFOPCs1(3,3,NZ)=0.57_r8
    CFOPCs1(3,4,NZ)=0.15_r8
!
!     DECIDUOUS AND CONIFEROUS TREES
!
  ELSE
    CFOPCs1(3,1,NZ)=0.0_r8
    CFOPCs1(3,2,NZ)=0.045_r8
    CFOPCs1(3,3,NZ)=0.660_r8
    CFOPCs1(3,4,NZ)=0.295_r8
  ENDIF
!
!     FRACTIONS OF FINE ROOT LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN PC&E 25:601-608
!
!     NON-VASCULAR
!
  IF(IGTYPs1(NZ).EQ.0)THEN
    CFOPCs1(4,1,NZ)=0.07_r8
    CFOPCs1(4,2,NZ)=0.25_r8
    CFOPCs1(4,3,NZ)=0.30_r8
    CFOPCs1(4,4,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
    CFOPCs1(4,1,NZ)=0.057_r8
    CFOPCs1(4,2,NZ)=0.263_r8
    CFOPCs1(4,3,NZ)=0.542_r8
    CFOPCs1(4,4,NZ)=0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(IBTYPs1(NZ).EQ.1.OR.IBTYPs1(NZ).EQ.3)THEN
    CFOPCs1(4,1,NZ)=0.059_r8
    CFOPCs1(4,2,NZ)=0.308_r8
    CFOPCs1(4,3,NZ)=0.464_r8
    CFOPCs1(4,4,NZ)=0.169_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPCs1(4,1,NZ)=0.059_r8
    CFOPCs1(4,2,NZ)=0.308_r8
    CFOPCs1(4,3,NZ)=0.464_r8
    CFOPCs1(4,4,NZ)=0.169_r8
  ENDIF
!
!     COARSE WOODY LITTER FROM BOLES AND ROOTS
!
  CFOPCs1(5,1,NZ)=0.00_r8
  CFOPCs1(5,2,NZ)=0.045_r8
  CFOPCs1(5,3,NZ)=0.660_r8
  CFOPCs1(5,4,NZ)=0.295_r8
!
!     INITIALIZE C-N AND C-P RATIOS IN PLANT LITTER
!
!     CNOPC,CPOPC=fractions to allocate N,P to kinetic components
!     CFOPN,CFOPP=distribution of litter N,P to kinetic components
!
  CNOPC(1)=0.020_r8
  CNOPC(2)=0.010_r8
  CNOPC(3)=0.010_r8
  CNOPC(4)=0.020_r8
  CPOPC(1)=0.0020_r8
  CPOPC(2)=0.0010_r8
  CPOPC(3)=0.0010_r8
  CPOPC(4)=0.0020_r8
  DO 110 N=0,5
    CNOPCT=0.0_r8
    CPOPCT=0.0_r8
    DO 100 M=1,4
      CNOPCT=CNOPCT+CFOPCs1(N,M,NZ)*CNOPC(M)
      CPOPCT=CPOPCT+CFOPCs1(N,M,NZ)*CPOPC(M)
100 CONTINUE
    DO 105 M=1,4
      CFOPNs1(N,M,NZ)=CFOPCs1(N,M,NZ)*CNOPC(M)/CNOPCT
      CFOPPs1(N,M,NZ)=CFOPCs1(N,M,NZ)*CPOPC(M)/CPOPCT
105 CONTINUE
110 CONTINUE
!
!     CONCURRENT NODE GROWTH
!
!     FNOD=scales node number for perennial vegetation (e.g. trees)
!     NNOD=number of concurrently growing nodes
!
  IF(IBTYPs1(NZ).EQ.0.OR.IGTYPs1(NZ).LE.1)THEN
    FNODs1(NZ)=1.0
    IF(GROUPIs1(NZ).LE.10)THEN
      NNODs1(NZ)=3
    ELSEIF(GROUPIs1(NZ).LE.15)THEN
      NNODs1(NZ)=4
    ELSE
      NNODs1(NZ)=5
    ENDIF
  ELSE
    FNODs1(NZ)=AMAX1(1.0,0.04/XRLAs1(NZ))
    NNODs1(NZ)=24
  ENDIF
  end associate
  end subroutine PlantLitterFraction
!------------------------------------------------------------------------------------------

  subroutine PFTThermalAcclimation(NZ)

  implicit none
  integer, intent(in) :: NZ
  real(r8), parameter :: TCZD = 5.0_r8        !basal value for threshold temperature for spring leafout/dehardening	oC
  real(r8), parameter :: TCXD = 12.0_r8       !basal value for threshold temperature for autumn leafoff/hardening	oC

  associate(                            &
    DATAPs1    =>  plt_site%DATAPs1   , &
    ICTYPs1    =>  plt_photo%ICTYPs1  , &
    HTCs1      =>  plt_pheno%HTCs1    , &
    TCXs1      =>  plt_pheno%TCXs1    , &
    TCZs1      =>  plt_pheno%TCZs1    , &
    OFFSTs1    =>  plt_pheno%OFFSTs1  , &
    ZTYPIs1    =>  plt_pheno%ZTYPIs1  , &
    ZTYPs1     =>  plt_pheno%ZTYPs1   , &
    SSTXs1     =>  plt_pheno%SSTXs1     &
  )
!
!     PFT THERMAL ACCLIMATION
!
!     ZTYP,ZTYPI=dynamic,initial thermal adaptation zone from PFT file
!     OFFST=shift in Arrhenius curve for thermal adaptation (oC)
!     TCZ,TCX=threshold temperature for leafout,leafoff
!     HTC=high temperature threshold for grain number loss (oC)
!     SSTX=sensitivity to HTC (seeds oC-1 above HTC)
!
  ZTYPs1(NZ)=ZTYPIs1(NZ)
  OFFSTs1(NZ)=2.667_r8*(2.5_r8-ZTYPs1(NZ))
  TCZs1(NZ)=TCZD-OFFSTs1(NZ)
  TCXs1(NZ)=AMIN1(15.0_r8,TCXD-OFFSTs1(NZ))
  IF(ICTYPs1(NZ).EQ.3)THEN
    IF(DATAPs1(NZ)(1:4).EQ.'soyb')THEN
      HTCs1(NZ)=30.0_r8+3.0_r8*ZTYPs1(NZ)
      SSTXs1(NZ)=0.002_r8
    ELSE
      HTCs1(NZ)=27.0_r8+3.0_r8*ZTYPs1(NZ)
      SSTXs1(NZ)=0.002_r8
    ENDIF
  ELSE
    HTCs1(NZ)=27.0_r8+3.0_r8*ZTYPs1(NZ)
    SSTXs1(NZ)=0.005_r8
  ENDIF
  end associate
  end subroutine PFTThermalAcclimation
!------------------------------------------------------------------------------------------

  subroutine InitDimensionsandUptake(NZ)

  implicit none
  integer, intent(in) :: NZ
  INTEGER :: L,N,NR
  associate(                             &
    CNRTSs1    =>  plt_allom%CNRTSs1   , &
    CPRTSs1    =>  plt_allom%CPRTSs1   , &
    DMRTs1     =>  plt_allom%DMRTs1    , &
    CNRTs1     =>  plt_allom%CNRTs1    , &
    CPRTs1     =>  plt_allom%CPRTs1    , &
    UPMNPOs1   =>  plt_rbgc%UPMNPOs1   , &
    UPKMPOs1   =>  plt_rbgc%UPKMPOs1   , &
    UPMXPOs1   =>  plt_rbgc%UPMXPOs1   , &
    UPMNZOs1   =>  plt_rbgc%UPMNZOs1   , &
    UPKMZOs1   =>  plt_rbgc%UPKMZOs1   , &
    UPMXZOs1   =>  plt_rbgc%UPMXZOs1   , &
    UPMNZHs1   =>  plt_rbgc%UPMNZHs1   , &
    UPMXZHs1   =>  plt_rbgc%UPMXZHs1   , &
    UPKMZHs1   =>  plt_rbgc%UPKMZHs1   , &
    CDPTHZs1   =>  plt_site%CDPTHZs1   , &
    NUs1       =>  plt_site%NUs1       , &
    NLs1       =>  plt_site%NLs1       , &
    NGs1       =>  plt_morph%NGs1      , &
    SDPTHs1    =>  plt_morph%SDPTHs1   , &
    SDARs1     =>  plt_morph%SDARs1    , &
    GRDMs1     =>  plt_morph%GRDMs1    , &
    SDPTHIs1   =>  plt_morph%SDPTHIs1  , &
    RRAD2Xs1   =>  plt_morph%RRAD2Xs1  , &
    RTAR2Xs1   =>  plt_morph%RTAR2Xs1  , &
    RRAD1Xs1   =>  plt_morph%RRAD1Xs1  , &
    RTAR1Xs1   =>  plt_morph%RTAR1Xs1  , &
    RRAD2Ms1   =>  plt_morph%RRAD2Ms1  , &
    RRAD1Ms1   =>  plt_morph%RRAD1Ms1  , &
    RTLG2Xs1   =>  plt_morph%RTLG2Xs1  , &
    NIXs1      =>  plt_morph%NIXs1     , &
    RSRRs1     =>  plt_morph%RSRRs1    , &
    RTLG1Xs1   =>  plt_morph%RTLG1Xs1  , &
    PORTs1     =>  plt_morph%PORTs1    , &
    PORTXs1    =>  plt_morph%PORTXs1   , &
    RRADPs1    =>  plt_morph%RRADPs1   , &
    DMVLs1     =>  plt_morph%DMVLs1    , &
    RSRAs1     =>  plt_morph%RSRAs1    , &
    NINRs1     =>  plt_morph%NINRs1    , &
    SDVLs1     =>  plt_morph%SDVLs1    , &
    SDLGs1     =>  plt_morph%SDLGs1      &
  )
!
!     SEED CHARACTERISTICS
!
!     SDVL,SDLG,SDAR=seed volume(m3),length(m),AREA3s1(NUs1)(m2)
!     GRDM=seed C mass (g) from PFT file
!
  SDVLs1(NZ)=GRDMs1(NZ)*5.0E-06
  SDLGs1(NZ)=2.0*(0.75*SDVLs1(NZ)/PICON)**0.33
  SDARs1(NZ)=4.0*PICON*(SDLGs1(NZ)/2.0)**2
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) DIMENSIONS, UPTAKE PARAMETERS
!
!     SDPTH=seeding depth(m) from PFT management file
!     CDPTHZ=depth to soil layer bottom from surface(m)
!     NG,NIX,NINR=seeding,upper,lower rooting layer
!     CNRTS,CPRTS=N,P root growth yield
!     RRAD1M,RRAD2M=maximum primary,secondary mycorrhizal radius (m)
!     PORT=mycorrhizal porosity
!     UPMXZH,UPKMZH,UPMNZH=NH4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake(g m-2 h-1),Km(uM), min concn (uM)
!     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     RSRR,RSRA=radial,axial root resistivity (m2 MPa-1 h-1)
!
  SDPTHs1(NZ)=SDPTHIs1(NZ)
  DO 9795 L=NUs1,NLs1
    IF(SDPTHs1(NZ).GE.CDPTHZs1(L-1).AND.SDPTHs1(NZ).LT.CDPTHZs1(L))THEN
      NGs1(NZ)=L
      NIXs1(NZ)=L
      DO 9790 NR=1,JC1
        NINRs1(NR,NZ)=L
9790  CONTINUE
    ENDIF
9795  CONTINUE
  CNRTSs1(NZ)=CNRTs1(NZ)*DMRTs1(NZ)
  CPRTSs1(NZ)=CPRTs1(NZ)*DMRTs1(NZ)
  RRAD1Ms1(2,NZ)=5.0E-06_r8
  RRAD2Ms1(2,NZ)=5.0E-06_r8
  PORTs1(2,NZ)=PORTs1(1,NZ)
  UPMXZHs1(2,NZ)=UPMXZHs1(1,NZ)
  UPKMZHs1(2,NZ)=UPKMZHs1(1,NZ)
  UPMNZHs1(2,NZ)=UPMNZHs1(1,NZ)
  UPMXZOs1(2,NZ)=UPMXZOs1(1,NZ)
  UPKMZOs1(2,NZ)=UPKMZOs1(1,NZ)
  UPMNZOs1(2,NZ)=UPMNZOs1(1,NZ)
  UPMXPOs1(2,NZ)=UPMXPOs1(1,NZ)
  UPKMPOs1(2,NZ)=UPKMPOs1(1,NZ)
  UPMNPOs1(2,NZ)=UPMNPOs1(1,NZ)
  RSRRs1(2,NZ)=1.0E+04_r8
  RSRAs1(2,NZ)=1.0E+12_r8
!
!     PORTX=tortuosity for gas transport
!     RRADP=path length for radial diffusion within root (m)
!     DMVL=volume:C ratio (m3 g-1)
!     RTLG1X,RTLG2X=specific primary,secondary root length (m g-1)
!     RTAR1X,RTAR2X=specific primary,secondary root area (m2 g-1)
!
  DO 500 N=1,2
    PORTXs1(N,NZ)=PORTs1(N,NZ)**1.33_r8
    RRADPs1(N,NZ)=LOG(1.0_r8/SQRT(AMAX1(0.01_r8,PORTs1(N,NZ))))
    DMVLs1(N,NZ)=1.0E-06_r8/(0.05_r8*(1.0-PORTs1(N,NZ)))
    RTLG1Xs1(N,NZ)=DMVLs1(N,NZ)/(PICON*RRAD1Ms1(N,NZ)**2)
    RTLG2Xs1(N,NZ)=DMVLs1(N,NZ)/(PICON*RRAD2Ms1(N,NZ)**2)
    RRAD1Xs1(N,NZ)=RRAD1Ms1(N,NZ)
!    2*SQRT(0.25*(1.0-PORT(N,NZ)))
    RRAD2Xs1(N,NZ)=RRAD2Ms1(N,NZ)
!    2*SQRT(0.25*(1.0-PORT(N,NZ)))
    RTAR1Xs1(N,NZ)=PICON*RRAD1Xs1(N,NZ)**2
    RTAR2Xs1(N,NZ)=PICON*RRAD2Xs1(N,NZ)**2
500 CONTINUE
  end associate
  end subroutine InitDimensionsandUptake
!------------------------------------------------------------------------------------------

  subroutine InitPlantPhenoMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: K,L,M,N,NB
  associate(                           &
    NUs1      =>  plt_site%NUs1      , &
    PPXs1     =>  plt_site%PPXs1     , &
    PPs1      =>  plt_site%PPs1      , &
    ALATs1    =>  plt_site%ALATs1    , &
    AREA3s1   =>  plt_site%AREA3s1   , &
    IGTYPs1   =>  plt_pheno%IGTYPs1  , &
    VSTGXs1   =>  plt_pheno%VSTGXs1  , &
    GROUPs1   =>  plt_pheno%GROUPs1  , &
    KVSTGs1   =>  plt_pheno%KVSTGs1  , &
    KVSTGNs1  =>  plt_pheno%KVSTGNs1 , &
    GSTGIs1   =>  plt_pheno%GSTGIs1  , &
    GSTGFs1   =>  plt_pheno%GSTGFs1  , &
    FLG4s1    =>  plt_pheno%FLG4s1   , &
    FLGZs1    =>  plt_pheno%FLGZs1   , &
    VRNYs1    =>  plt_pheno%VRNYs1   , &
    VRNFs1    =>  plt_pheno%VRNFs1   , &
    IDTHBs1   =>  plt_pheno%IDTHBs1  , &
    VRNZs1    =>  plt_pheno%VRNZs1   , &
    IDAYs1    =>  plt_pheno%IDAYs1   , &
    VRNSs1    =>  plt_pheno%VRNSs1   , &
    ATRPs1    =>  plt_pheno%ATRPs1   , &
    WSTRs1    =>  plt_pheno%WSTRs1   , &
    GROUPIs1  =>  plt_pheno%GROUPIs1 , &
    RCCSXs1   =>  plt_pheno%RCCSXs1  , &
    TGSTGFs1  =>  plt_pheno%TGSTGFs1 , &
    TGSTGIs1  =>  plt_pheno%TGSTGIs1 , &
    FDBKXs1   =>  plt_photo%FDBKXs1  , &
    CPOOL3s1  =>  plt_photo%CPOOL3s1 , &
    CPOOL4s1  =>  plt_photo%CPOOL4s1 , &
    CHILLs1   =>  plt_photo%CHILLs1  , &
    FDBKs1    =>  plt_photo%FDBKs1   , &
    HCOBs1    =>  plt_photo%HCOBs1   , &
    CO2Bs1    =>  plt_photo%CO2Bs1   , &
    VSTGs1    =>  plt_morph%VSTGs1   , &
    ZCs1      =>  plt_morph%ZCs1     , &
    KLEAFs1   =>  plt_morph%KLEAFs1  , &
    XTLIs1    =>  plt_morph%XTLIs1   , &
    NBTs1     =>  plt_morph%NBTs1    , &
    PSTGs1    =>  plt_morph%PSTGs1   , &
    ARLFPs1   =>  plt_morph%ARLFPs1  , &
    ARSTKs1   =>  plt_morph%ARSTKs1  , &
    SURFBs1   =>  plt_morph%SURFBs1  , &
    ARSTPs1   =>  plt_morph%ARSTPs1  , &
    PSTGFs1   =>  plt_morph%PSTGFs1  , &
    GRNXBs1   =>  plt_morph%GRNXBs1  , &
    HTNODEs1  =>  plt_morph%HTNODEs1 , &
    GRNOBs1   =>  plt_morph%GRNOBs1  , &
    ARLFBs1   =>  plt_morph%ARLFBs1  , &
    ARLFZs1   =>  plt_morph%ARLFZs1  , &
    ARLFLs1   =>  plt_morph%ARLFLs1  , &
    ARLFVs1   =>  plt_morph%ARLFVs1  , &
    ARSTVs1   =>  plt_morph%ARSTVs1  , &
    SURFs1    =>  plt_morph%SURFs1   , &
    ARLF1s1   =>  plt_morph%ARLF1s1  , &
    HTSHEXs1  =>  plt_morph%HTSHEXs1 , &
    HTCTLs1   =>  plt_morph%HTCTLs1  , &
    NBTBs1    =>  plt_morph%NBTBs1   , &
    PSTGIs1   =>  plt_morph%PSTGIs1  , &
    KLEAFXs1  =>  plt_morph%KLEAFXs1 , &
    NBRs1     =>  plt_morph%NBRs1      &
  )
!
!     INITIALIZE PLANT PHENOLOGY
!
!     PP=population (grid cell-1)
!
  PPs1(NZ)=PPXs1(NZ)*AREA3s1(NUs1)
  plt_pheno%IFLGIs1(NZ)=0
  plt_pheno%IDTHPs1(NZ)=0
  plt_pheno%IDTHRs1(NZ)=0
  NBTs1(NZ)=0
  NBRs1(NZ)=0
  HTCTLs1(NZ)=0._r8
  ZCs1(NZ)=0._r8
  DO 10 NB=1,10
    plt_pheno%IFLGAs1(NB,NZ)=0
    plt_pheno%IFLGEs1(NB,NZ)=0
    plt_pheno%IFLGFs1(NB,NZ)=0
    plt_pheno%IFLGRs1(NB,NZ)=0
    plt_pheno%IFLGQs1(NB,NZ)=0
    GROUPs1(NB,NZ)=GROUPIs1(NZ)
    PSTGs1(NB,NZ)=XTLIs1(NZ)
    PSTGIs1(NB,NZ)=PSTGs1(NB,NZ)
    PSTGFs1(NB,NZ)=0._r8
    VSTGs1(NB,NZ)=0._r8
    VSTGXs1(NB,NZ)=0._r8
    KLEAFs1(NB,NZ)=1
    KLEAFXs1(NB,NZ)=1
    KVSTGs1(NB,NZ)=1
    KVSTGNs1(NB,NZ)=0
    GSTGIs1(NB,NZ)=0._r8
    GSTGFs1(NB,NZ)=0._r8
    TGSTGIs1(NB,NZ)=0._r8
    TGSTGFs1(NB,NZ)=0._r8
    VRNYs1(NB,NZ)=0._r8
    VRNZs1(NB,NZ)=0._r8
    VRNSs1(NB,NZ)=VRNYs1(NB,NZ)
    VRNFs1(NB,NZ)=VRNZs1(NB,NZ)
    ATRPs1(NB,NZ)=0._r8
    FDBKs1(NB,NZ)=1.0
    FDBKXs1(NB,NZ)=1.0
    FLG4s1(NB,NZ)=0
    FLGZs1(NB,NZ)=0
    NBTBs1(NB,NZ)=0
    plt_pheno%IDTHBs1(NB,NZ)=1
    DO 15 M=1,10
      IDAYs1(M,NB,NZ)=0
15  CONTINUE
10  CONTINUE
!
!     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
!
  WSTRs1(NZ)=0._r8
  CHILLs1(NZ)=0._r8
  DO 25 NB=1,10
    plt_biom%CPOOLs1(NB,NZ)=0._r8
    plt_biom%ZPOOLs1(NB,NZ)=0._r8
    plt_biom%PPOOLs1(NB,NZ)=0._r8
    plt_biom%CPOLNBs1(NB,NZ)=0._r8
    plt_biom%ZPOLNBs1(NB,NZ)=0._r8
    plt_biom%PPOLNBs1(NB,NZ)=0._r8
    plt_biom%WTSHTBs1(NB,NZ)=0._r8
    plt_biom%WTLFBs1(NB,NZ)=0._r8
    plt_biom%WTNDBs1(NB,NZ)=0._r8
    plt_biom%WTSHEBs1(NB,NZ)=0._r8
    plt_biom%WTSTKBs1(NB,NZ)=0._r8
    plt_biom%WVSTKBs1(NB,NZ)=0._r8
    plt_biom%WTRSVBs1(NB,NZ)=0._r8
    plt_biom%WTHSKBs1(NB,NZ)=0._r8
    plt_biom%WTEARBs1(NB,NZ)=0._r8
    plt_biom%WTGRBs1(NB,NZ)=0._r8
    plt_biom%WTLSBs1(NB,NZ)=0._r8
    plt_biom%WTSHTNs1(NB,NZ)=0._r8
    plt_biom%WTLFBNs1(NB,NZ)=0._r8
    plt_biom%WTNDBNs1(NB,NZ)=0._r8
    plt_biom%WTSHBNs1(NB,NZ)=0._r8
    plt_biom%WTSTBNs1(NB,NZ)=0._r8
    plt_biom%WTRSBNs1(NB,NZ)=0._r8
    plt_biom%WTHSBNs1(NB,NZ)=0._r8
    plt_biom%WTEABNs1(NB,NZ)=0._r8
    plt_biom%WTGRBNs1(NB,NZ)=0._r8
    plt_biom%WTSHTPs1(NB,NZ)=0._r8
    plt_biom%WTLFBPs1(NB,NZ)=0._r8
    plt_biom%WTNDBPs1(NB,NZ)=0._r8
    plt_biom%WTSHBPs1(NB,NZ)=0._r8
    plt_biom%WTSTBPs1(NB,NZ)=0._r8
    plt_biom%WTRSBPs1(NB,NZ)=0._r8
    plt_biom%WTHSBPs1(NB,NZ)=0._r8
    plt_biom%WTEABPs1(NB,NZ)=0._r8
    plt_biom%WTGRBPs1(NB,NZ)=0._r8
    GRNXBs1(NB,NZ)=0._r8
    GRNOBs1(NB,NZ)=0._r8
    plt_allom%GRWTBs1(NB,NZ)=0._r8
    ARLFBs1(NB,NZ)=0._r8
    plt_rbgc%RNH3Bs1(NB,NZ)=0._r8
    plt_pheno%RCZLXs1(NB,NZ)=0._r8
    plt_pheno%RCPLXs1(NB,NZ)=0._r8
    plt_pheno%RCCLXs1(NB,NZ)=0._r8
    plt_biom%WGLFXs1(NB,NZ)=0._r8
    plt_biom%WGLFNXs1(NB,NZ)=0._r8
    plt_biom%WGLFPXs1(NB,NZ)=0._r8
    ARLFZs1(NB,NZ)=0._r8
    plt_pheno%RCZSXs1(NB,NZ)=0._r8
    plt_pheno%RCPSXs1(NB,NZ)=0._r8
    plt_pheno%RCCSXs1(NB,NZ)=0._r8
    plt_biom%WTSTXBs1(NB,NZ)=0._r8
    plt_biom%WTSTXNs1(NB,NZ)=0._r8
    plt_biom%WTSTXPs1(NB,NZ)=0._r8
    plt_biom%WGSHEXs1(NB,NZ)=0._r8
    plt_biom%WGSHNXs1(NB,NZ)=0._r8
    plt_biom%WGSHPXs1(NB,NZ)=0._r8
    HTSHEXs1(NB,NZ)=0._r8
    DO 5 L=1,JC1
      ARSTKs1(L,NB,NZ)=0._r8
      DO N=1,JLI1
        SURFBs1(N,L,NB,NZ)=0._r8
      enddo
5   CONTINUE
    DO K=0,JNODS1
      ARLF1s1(K,NB,NZ)=0._r8
      HTNODEs1(K,NB,NZ)=0._r8
      plt_morph%HTNODXs1(K,NB,NZ)=0._r8
      plt_morph%HTSHEs1(K,NB,NZ)=0._r8
      plt_biom%WGLFs1(K,NB,NZ)=0._r8
      plt_biom%WSLFs1(K,NB,NZ)=0._r8
      plt_biom%WGLFNs1(K,NB,NZ)=0._r8
      plt_biom%WGLFPs1(K,NB,NZ)=0._r8
      plt_biom%WGSHEs1(K,NB,NZ)=0._r8
      plt_biom%WSSHEs1(K,NB,NZ)=0._r8
      plt_biom%WGSHNs1(K,NB,NZ)=0._r8
      plt_biom%WGSHPs1(K,NB,NZ)=0._r8
      plt_biom%WGNODEs1(K,NB,NZ)=0._r8
      plt_biom%WGNODNs1(K,NB,NZ)=0._r8
      plt_biom%WGNODPs1(K,NB,NZ)=0._r8
      DO 55 L=1,JC1
        ARLFLs1(L,K,NB,NZ)=0._r8
        plt_biom%WGLFLs1(L,K,NB,NZ)=0._r8
        plt_biom%WGLFLNs1(L,K,NB,NZ)=0._r8
        plt_biom%WGLFLPs1(L,K,NB,NZ)=0._r8
55    CONTINUE
      IF(K.NE.0)THEN
        CPOOL3s1(K,NB,NZ)=0._r8
        CO2Bs1(K,NB,NZ)=0._r8
        HCOBs1(K,NB,NZ)=0._r8
        CPOOL4s1(K,NB,NZ)=0._r8
        DO 45 L=1,JC1
          DO N=1,JLI1
            SURFs1(N,L,K,NB,NZ)=0._r8
          enddo
45      CONTINUE
      ENDIF
    enddo
25  CONTINUE
  DO 35 L=1,JC1
    ARLFVs1(L,NZ)=0._r8
    plt_biom%WGLFVs1(L,NZ)=0._r8
    ARSTVs1(L,NZ)=0._r8
35  CONTINUE
  plt_biom%CPOOLPs1(NZ)=0._r8
  plt_biom%ZPOOLPs1(NZ)=0._r8
  plt_biom%PPOOLPs1(NZ)=0._r8
  plt_biom%CCPOLPs1(NZ)=0._r8
  plt_biom%CCPLNPs1(NZ)=0._r8
  plt_biom%CZPOLPs1(NZ)=0._r8
  plt_biom%CPPOLPs1(NZ)=0._r8
  plt_biom%WTSHTs1(NZ)=0._r8
  plt_biom%WTLFs1(NZ)=0._r8
  plt_biom%WTSHEs1(NZ)=0._r8
  plt_biom%WTSTKs1(NZ)=0._r8
  plt_biom%WVSTKs1(NZ)=0._r8
  plt_biom%WTRSVs1(NZ)=0._r8
  plt_biom%WTHSKs1(NZ)=0._r8
  plt_biom%WTEARs1(NZ)=0._r8
  plt_biom%WTGRs1(NZ)=0._r8
  plt_biom%WTRTts1(NZ)=0._r8
  plt_biom%WTRTSs1(NZ)=0._r8
  plt_biom%WTNDs1(NZ)=0._r8
  plt_biom%WTLSs1(NZ)=0._r8
  plt_biom%WTSHNs1(NZ)=0._r8
  plt_biom%WTLFNs1(NZ)=0._r8
  plt_biom%WTSHENs1(NZ)=0._r8
  plt_biom%WTSTKNs1(NZ)=0._r8
  plt_biom%WTRSVNs1(NZ)=0._r8
  plt_biom%WTHSKNs1(NZ)=0._r8
  plt_biom%WTEARNs1(NZ)=0._r8
  plt_biom%WTGRNNs1(NZ)=0._r8
  plt_biom%WTNDNs1(NZ)=0._r8
  plt_biom%WTSHPs1(NZ)=0._r8
  plt_biom%WTLFPs1(NZ)=0._r8
  plt_biom%WTSHEPs1(NZ)=0._r8
  plt_biom%WTSTKPs1(NZ)=0._r8
  plt_biom%WTRSVPs1(NZ)=0._r8
  plt_biom%WTHSKPs1(NZ)=0._r8
  plt_biom%WTEARPs1(NZ)=0._r8
  plt_biom%WTGRNPs1(NZ)=0._r8
  plt_biom%WTNDPs1(NZ)=0._r8
  ARLFPs1(NZ)=0._r8
  plt_biom%WTRTAs1(NZ)=0._r8
  ARSTPs1(NZ)=0._r8
  end associate
  end subroutine InitPlantPhenoMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitMassBalance(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: M
  real(r8) :: WTSTDX

  associate(                           &
    NUs1     => plt_site%NUs1        , &
    AREA3s1  => plt_site%AREA3s1     , &
    CFOPCs1  => plt_soilchem%CFOPCs1 , &
    CFOPNs1  => plt_soilchem%CFOPNs1 , &
    CFOPPs1  => plt_soilchem%CFOPPs1 , &
    CTRANs1  => plt_ew%CTRANs1       , &
    WTSTDGs1 => plt_biom%WTSTDGs1    , &
    WTSTDNs1 => plt_biom%WTSTDNs1    , &
    WTSTDPs1 => plt_biom%WTSTDPs1    , &
    WTSTGs1  => plt_biom%WTSTGs1     , &
    WTSTGPs1 => plt_biom%WTSTGPs1    , &
    WTSTDIs1 => plt_biom%WTSTDIs1    , &
    WTSTGNs1 => plt_biom%WTSTGNs1    , &
    WTRTAs1  => plt_biom%WTRTAs1     , &
    CNSTKs1  => plt_allom%CNSTKs1    , &
    CPSTKs1  => plt_allom%CPSTKs1    , &
    TCUPTKs1 => plt_rbgc%TCUPTKs1    , &
    TZUPTKs1 => plt_rbgc%TZUPTKs1    , &
    TPUPTKs1 => plt_rbgc%TPUPTKs1    , &
    TNH3Cs1  => plt_bgcr%TNH3Cs1     , &
    RNH3Cs1  => plt_bgcr%RNH3Cs1     , &
    TPSN0s1  => plt_bgcr%TPSN0s1     , &
    TZSN0s1  => plt_bgcr%TZSN0s1     , &
    TCSN0s1  => plt_bgcr%TCSN0s1     , &
    CARBNs1  => plt_bgcr%CARBNs1     , &
    TCO2As1  => plt_bgcr%TCO2As1     , &
    TCO2Ts1  => plt_bgcr%TCO2Ts1     , &
    TZUPFXs1 => plt_bgcr%TZUPFXs1    , &
    TZSNCs1  => plt_bgcr%TZSNCs1     , &
    TPSNCs1  => plt_bgcr%TPSNCs1     , &
    TCSNCs1  => plt_bgcr%TCSNCs1     , &
    ARSTPs1  => plt_morph%ARSTPs1    , &
    RSETCs1  => plt_pheno%RSETCs1    , &
    RSETNs1  => plt_pheno%RSETNs1    , &
    RSETPs1  => plt_pheno%RSETPs1      &

  )
!
!     INITIALIZE MASS BALANCE CHECKS
!
  IF(.not.is_restart_run.AND.is_first_year)THEN
    CARBNs1(NZ)=0._r8
    TCSN0s1(NZ)=0._r8
    TZSN0s1(NZ)=0._r8
    TPSN0s1(NZ)=0._r8
    TCO2Ts1(NZ)=0._r8
    TCO2As1(NZ)=0._r8
    TCUPTKs1(NZ)=0._r8
    TCSNCs1(NZ)=0._r8
    TZUPTKs1(NZ)=0._r8
    TZSNCs1(NZ)=0._r8
    TPUPTKs1(NZ)=0._r8
    TPSNCs1(NZ)=0._r8
    TZUPFXs1(NZ)=0._r8
    RNH3Cs1(NZ)=0._r8
    TNH3Cs1(NZ)=0._r8
    plt_distb%VCO2Fs1(NZ)=0._r8
    plt_distb%VCH4Fs1(NZ)=0._r8
    plt_distb%VOXYFs1(NZ)=0._r8
    plt_distb%VNH3Fs1(NZ)=0._r8
    plt_distb%VN2OFs1(NZ)=0._r8
    plt_distb%VPO4Fs1(NZ)=0._r8
    plt_distb%THVSTCs1(NZ)=0._r8
    plt_distb%THVSTNs1(NZ)=0._r8
    plt_distb%THVSTPs1(NZ)=0._r8
    plt_distb%HVSTCs1(NZ)=0._r8
    plt_distb%HVSTNs1(NZ)=0._r8
    plt_distb%HVSTPs1(NZ)=0._r8
    RSETCs1(NZ)=0._r8
    RSETNs1(NZ)=0._r8
    RSETPs1(NZ)=0._r8
    CTRANs1(NZ)=0._r8
    WTSTGs1(NZ)=0._r8
    WTSTGNs1(NZ)=0._r8
    WTSTGPs1(NZ)=0._r8
    WTSTDX=WTSTDIs1(NZ)*AREA3s1(NUs1)
    DO 155 M=1,4
      WTSTDGs1(M,NZ)=WTSTDX*CFOPCs1(5,M,NZ)
      WTSTDNs1(M,NZ)=WTSTDX*CNSTKs1(NZ)*CFOPNs1(5,M,NZ)
      WTSTDPs1(M,NZ)=WTSTDX*CPSTKs1(NZ)*CFOPPs1(5,M,NZ)
      WTSTGs1(NZ)=WTSTGs1(NZ)+WTSTDGs1(M,NZ)
      WTSTGNs1(NZ)=WTSTGNs1(NZ)+WTSTDNs1(M,NZ)
      WTSTGPs1(NZ)=WTSTGPs1(NZ)+WTSTDPs1(M,NZ)
155 CONTINUE
  ENDIF
  end associate
  end subroutine InitMassBalance
!------------------------------------------------------------------------------------------

  subroutine InitPlantHeatandWater(NZ)

  implicit none
  integer, intent(in) :: NZ
  associate(                          &
    ATCAs1     =>  plt_site%ATCAs1  , &
    OSMOs1     =>  plt_ew%OSMOs1    , &
    TKCs1      =>  plt_ew%TKCs1     , &
    EPs1       =>  plt_ew%EPs1      , &
    VHCPCs1    =>  plt_ew%VHCPCs1   , &
    PSILTs1    =>  plt_ew%PSILTs1   , &
    PSILGs1    =>  plt_ew%PSILGs1   , &
    PSILOs1    =>  plt_ew%PSILOs1   , &
    ENGYXs1    =>  plt_ew%ENGYXs1   , &
    DTKCs1     =>  plt_ew%DTKCs1    , &
    TCCs1      =>  plt_ew%TCCs1     , &
    TKGs1      =>  plt_pheno%TKGs1  , &
    TCGs1      =>  plt_pheno%TCGs1  , &
    TFN3s1     =>  plt_pheno%TFN3s1 , &
    WTSHTs1    =>  plt_biom%WTSHTs1 , &
    FRADPs1    =>  plt_rad%FRADPs1    &
  )
!
!     INITIALIZE PLANT HEAT AND WATER STATUS
!
!     VHCPC=canopy heat capacity (MJ m-3 K-1)
!     TCC,TKC=canopy temperature for growth (oC,K)
!     TCG,TKG=canopy temperature for phenology (oC,K)
!     PSILT,PSILO,PSILG=canopy total,osmotic,turgor water potl(MPa)
!
  VHCPCs1(NZ)=cpw*WTSHTs1(NZ)*10.0E-06
  ENGYXs1(NZ)=0._r8
  DTKCs1(NZ)=0._r8
  TCCs1(NZ)=ATCAs1
  TKCs1(NZ)=TCCs1(NZ)+TC2K
  TCGs1(NZ)=TCCs1(NZ)
  TKGs1(NZ)=TCGs1(NZ)+TC2K
  TFN3s1(NZ)=1.0
  PSILTs1(NZ)=-1.0E-03
  PSILOs1(NZ)=OSMOs1(NZ)+PSILTs1(NZ)
  PSILGs1(NZ)=AMAX1(0.0,PSILTs1(NZ)-PSILOs1(NZ))
  EPs1(NZ)=0._r8
  FRADPs1(NZ)=0._r8
  end associate
  end subroutine InitPlantHeatandWater
!------------------------------------------------------------------------------------------

  subroutine InitRootMychorMorphoBio(NZ)
  implicit none
  integer, intent(in) :: NZ
  integer :: K,L,M,N,NR
  REAL(R8) :: CCO2A
  REAL(R8) :: CCO2P
  REAL(R8) :: COXYA
  REAL(R8) :: COXYP
  associate(                             &
    OSMOs1     =>  plt_ew%OSMOs1       , &
    PSILTs1    =>  plt_ew%PSILTs1      , &
    PSIRTs1    =>  plt_ew%PSIRTs1      , &
    PSIRGs1    =>  plt_ew%PSIRGs1      , &
    PSIROs1    =>  plt_ew%PSIROs1      , &
    RUPNFs1    =>  plt_bgcr%RUPNFs1    , &
    CO2Ps1     =>  plt_rbgc%CO2Ps1     , &
    CO2As1     =>  plt_rbgc%CO2As1     , &
    OXYEs1     =>  plt_site%OXYEs1     , &
    COXYEs1    =>  plt_site%COXYEs1    , &
    CO2EIs1    =>  plt_site%CO2EIs1    , &
    NLs1       =>  plt_site%NLs1       , &
    ATCAs1     =>  plt_site%ATCAs1     , &
    CCO2EIs1   =>  plt_site%CCO2EIs1   , &
    CWSRTs1    =>  plt_allom%CWSRTs1   , &
    CWSRTLs1   =>  plt_biom%CWSRTLs1   , &
    WSRTLs1    =>  plt_biom%WSRTLs1    , &
    SDPTHs1    =>  plt_morph%SDPTHs1   , &
    RTVLWs1    =>  plt_morph%RTVLWs1   , &
    RTVLPs1    =>  plt_morph%RTVLPs1   , &
    RRAD2s1    =>  plt_morph%RRAD2s1   , &
    RRAD1s1    =>  plt_morph%RRAD1s1   , &
    RRAD1Ms1   =>  plt_morph%RRAD1Ms1  , &
    RRAD2Ms1   =>  plt_morph%RRAD2Ms1  , &
    NRTs1      =>  plt_morph%NRTs1       &
  )
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
!
!     PSIRT,PSIRO,PSIRG=root,myco total,osmotic,turgor water potl(MPa)
!     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
!     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
!
  NRTs1(NZ)=0
  plt_rbgc%UPNH4s1(NZ)=0._r8
  plt_rbgc%UPNO3s1(NZ)=0._r8
  plt_rbgc%UPH2Ps1(NZ)=0._r8
  plt_rbgc%UPH1Ps1(NZ)=0._r8
  plt_rbgc%UPNFs1(NZ)=0._r8
  DO 40 N=1,2
    DO 20 L=1,NLs1
      plt_ew%UPWTRs1(N,L,NZ)=0._r8
      PSIRTs1(N,L,NZ)=-0.01
      PSIROs1(N,L,NZ)=OSMOs1(NZ)+PSIRTs1(N,L,NZ)
      PSIRGs1(N,L,NZ)=AMAX1(0.0_r8,PSIRTs1(N,L,NZ)-PSIROs1(N,L,NZ))
      plt_biom%CPOOLRs1(N,L,NZ)=0._r8
      plt_biom%ZPOOLRs1(N,L,NZ)=0._r8
      plt_biom%PPOOLRs1(N,L,NZ)=0._r8
      plt_biom%CCPOLRs1(N,L,NZ)=0._r8
      plt_biom%CZPOLRs1(N,L,NZ)=0._r8
      plt_biom%CPPOLRs1(N,L,NZ)=0._r8
      CWSRTLs1(N,L,NZ)=CWSRTs1(NZ)
      plt_biom%WTRTLs1(N,L,NZ)=0._r8
      plt_biom%WTRTDs1(N,L,NZ)=0._r8
      WSRTLs1(N,L,NZ)=0._r8
      plt_morph%RTN1s1(N,L,NZ)=0._r8
      plt_morph%RTNLs1(N,L,NZ)=0._r8
      plt_morph%RTLGPs1(N,L,NZ)=0._r8
      plt_morph%RTDNPs1(N,L,NZ)=0._r8
      RTVLPs1(N,L,NZ)=0._r8
      RTVLWs1(N,L,NZ)=0._r8
      RRAD1s1(N,L,NZ)=RRAD1Ms1(N,NZ)
      RRAD2s1(N,L,NZ)=RRAD2Ms1(N,NZ)
      plt_morph%RTARPs1(N,L,NZ)=0._r8
      plt_morph%RTLGAs1(N,L,NZ)=1.0E-03
      plt_rbgc%RUPNH4s1(N,L,NZ)=0._r8
      plt_rbgc%RUPNO3s1(N,L,NZ)=0._r8
      plt_rbgc%RUPH2Ps1(N,L,NZ)=0._r8
      plt_rbgc%RUPH1Ps1(N,L,NZ)=0._r8
      plt_rbgc%RUPNHBs1(N,L,NZ)=0._r8
      plt_rbgc%RUPNOBs1(N,L,NZ)=0._r8
      plt_rbgc%RUPH2Bs1(N,L,NZ)=0._r8
      plt_rbgc%RUPH1Bs1(N,L,NZ)=0._r8
      plt_rbgc%ROXYPs1(N,L,NZ)=0._r8
      plt_rbgc%RUNNHPs1(N,L,NZ)=0._r8
      plt_rbgc%RUNNBPs1(N,L,NZ)=0._r8
      plt_rbgc%RUNNOPs1(N,L,NZ)=0._r8
      plt_rbgc%RUNNXPs1(N,L,NZ)=0._r8
      plt_rbgc%RUPP2Ps1(N,L,NZ)=0._r8
      plt_rbgc%RUPP1Ps1(N,L,NZ)=0._r8
      plt_rbgc%RUPP2Bs1(N,L,NZ)=0._r8
      plt_rbgc%RUPP1Bs1(N,L,NZ)=0._r8
      CCO2A=CCO2EIs1
      CCO2P=0.030*EXP(-2.621_r8-0.0317_r8*ATCAs1)*CO2EIs1
      CO2As1(N,L,NZ)=CCO2A*RTVLPs1(N,L,NZ)
      CO2Ps1(N,L,NZ)=CCO2P*RTVLWs1(N,L,NZ)
      plt_rbgc%RCOFLAs1(N,L,NZ)=0._r8
      plt_rbgc%RCODFAs1(N,L,NZ)=0._r8
      plt_rbgc%RCO2Ss1(N,L,NZ)=0._r8
      plt_rbgc%RCO2Ps1(N,L,NZ)=0._r8
      COXYA=COXYEs1
      COXYP=0.032*EXP(-6.175_r8-0.0211_r8*ATCAs1)*OXYEs1
      plt_rbgc%OXYAs1(N,L,NZ)=COXYA*RTVLPs1(N,L,NZ)
      plt_rbgc%OXYPs1(N,L,NZ)=COXYP*RTVLWs1(N,L,NZ)
      plt_rbgc%CH4As1(N,L,NZ)=0._r8
      plt_rbgc%CH4Ps1(N,L,NZ)=0._r8
      plt_rbgc%Z2OAs1(N,L,NZ)=0._r8
      plt_rbgc%Z2OPs1(N,L,NZ)=0._r8
      plt_rbgc%ZH3As1(N,L,NZ)=0._r8
      plt_rbgc%ZH3Ps1(N,L,NZ)=0._r8
      plt_rbgc%H2GAs1(N,L,NZ)=0._r8
      plt_rbgc%H2GPs1(N,L,NZ)=0._r8
      plt_rbgc%WFRs1(N,L,NZ)=1.0
      DO 30 NR=1,JC1
        plt_morph%RTN2s1(N,L,NR,NZ)=0._r8
        plt_morph%RTLG1s1(N,L,NR,NZ)=0._r8
        plt_morph%RTLG2s1(N,L,NR,NZ)=0._r8
        plt_morph%RTDP1s1(N,NR,NZ)=SDPTHs1(NZ)
        plt_biom%WTRT1s1(N,L,NR,NZ)=0._r8
        plt_biom%WTRT1Ns1(N,L,NR,NZ)=0._r8
        plt_biom%WTRT1Ps1(N,L,NR,NZ)=0._r8
        plt_biom%WTRT2s1(N,L,NR,NZ)=0._r8
        plt_biom%WTRT2Ns1(N,L,NR,NZ)=0._r8
        plt_biom%WTRT2Ps1(N,L,NR,NZ)=0._r8
        plt_biom%RTWT1s1(N,NR,NZ)=0._r8
        plt_biom%RTWT1Ns1(N,NR,NZ)=0._r8
        plt_biom%RTWT1Ps1(N,NR,NZ)=0._r8
30    CONTINUE
      IF(N.EQ.1)THEN
        DO 6400 K=0,1
          DO  M=1,jsken
            plt_bgcr%CSNCs1(M,K,L,NZ)=0._r8
            plt_bgcr%ZSNCs1(M,K,L,NZ)=0._r8
            plt_bgcr%PSNCs1(M,K,L,NZ)=0._r8
          enddo
6400    CONTINUE
        plt_biom%CPOOLNs1(L,NZ)=0._r8
        plt_biom%ZPOOLNs1(L,NZ)=0._r8
        plt_biom%PPOOLNs1(L,NZ)=0._r8
        plt_biom%WTNDLs1(L,NZ)=0._r8
        plt_biom%WTNDLNs1(L,NZ)=0._r8
        plt_biom%WTNDLPs1(L,NZ)=0._r8
        RUPNFs1(L,NZ)=0._r8
      ENDIF
20  CONTINUE
40  CONTINUE

  plt_rbgc%RUPNH4s1(1:2,NLs1+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPNHBs1(1:2,NLs1+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPH2Ps1(1:2,NLs1+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPH2Bs1(1:2,NLs1+1:JZ1,NZ)=0._r8
  plt_morph%RTDNPs1(1:2,NLs1+1:JZ1,NZ)=0._r8
  end associate
  end subroutine InitRootMychorMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitSeedMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  REAL(R8) :: FDM

  associate(                             &
    PPs1       =>   plt_site%PPs1      , &
    PSILTs1    =>   plt_ew%PSILTs1     , &
    VOLWCs1    =>   plt_ew%VOLWCs1     , &
    VOLWPs1    =>   plt_ew%VOLWPs1     , &
    CWSRTs1    =>   plt_allom%CWSRTs1  , &
    CNGRs1     =>   plt_allom%CNGRs1   , &
    CPGRs1     =>   plt_allom%CPGRs1   , &
    RTWT1s1    =>   plt_biom%RTWT1s1   , &
    RTWT1Ns1   =>   plt_biom%RTWT1Ns1  , &
    RTWT1Ps1   =>   plt_biom%RTWT1Ps1  , &
    WTRVXs1    =>   plt_biom%WTRVXs1   , &
    WTRT1s1    =>   plt_biom%WTRT1s1   , &
    WTRT1Ns1   =>   plt_biom%WTRT1Ns1  , &
    WTRT1Ps1   =>   plt_biom%WTRT1Ps1  , &
    WTLFBNs1   =>   plt_biom%WTLFBNs1  , &
    WTLSBs1    =>   plt_biom%WTLSBs1   , &
    WTLSs1     =>   plt_biom%WTLSs1    , &
    WTRTDs1    =>   plt_biom%WTRTDs1   , &
    WTSHEBs1   =>   plt_biom%WTSHEBs1  , &
    WTRTLs1    =>   plt_biom%WTRTLs1   , &
    WSRTLs1    =>   plt_biom%WSRTLs1   , &
    CPOOLRs1   =>   plt_biom%CPOOLRs1  , &
    ZPOOLRs1   =>   plt_biom%ZPOOLRs1  , &
    PPOOLRs1   =>   plt_biom%PPOOLRs1  , &
    CPOOLs1    =>   plt_biom%CPOOLs1   , &
    ZPOOLs1    =>   plt_biom%ZPOOLs1   , &
    PPOOLs1    =>   plt_biom%PPOOLs1   , &
    WTLFBs1    =>   plt_biom%WTLFBs1   , &
    WTLFBPs1   =>   plt_biom%WTLFBPs1  , &
    WTRVCs1    =>   plt_biom%WTRVCs1   , &
    WTRVNs1    =>   plt_biom%WTRVNs1   , &
    WTRVPs1    =>   plt_biom%WTRVPs1   , &
    GRDMs1     =>   plt_morph%GRDMs1   , &
    RTDP1s1    =>   plt_morph%RTDP1s1  , &
    NGs1       =>   plt_morph%NGs1       &
  )
!
!     INITIALIZE SEED MORPHOLOGY AND BIOMASS
!
!     WTRVC,WTRVN,WTRVP=C,N,P in storage reserves (g)
!     WTLFB,WTLFBN,WTLFBP=C,N,P in leaves (g)
!     WTLSB=C in leaves+petioles (g)
!     FDM-dry matter fraction (g DM C g FM C-1)
!     VOLWP,VOLWC=water volume in,on canopy (m3)
!     CPOOL,ZPOOL,PPOOL=C,N,P in canopy nonstructural pools (g)
!     WTRT1,WTRT1N,WTRT1P=C,N,P in primary root layer (g)
!     RTWT1,RTWT1N,RTWT1P=total C,N,P in primary root (g)
!     WTRTL,WTRTD=total root C mass (g)
!     WSRTL=total root protein C mass (g)
!     CPOOLR,ZPOOLR,PPOOLR=C,N,P in root,myco nonstructural pools (g)
!
  WTRVXs1(NZ)=GRDMs1(NZ)*PPs1(NZ)
  WTRVCs1(NZ)=WTRVXs1(NZ)
  WTRVNs1(NZ)=CNGRs1(NZ)*WTRVCs1(NZ)
  WTRVPs1(NZ)=CPGRs1(NZ)*WTRVCs1(NZ)
  WTLFBNs1(1,NZ)=CNGRs1(NZ)*WTLFBs1(1,NZ)
  WTLFBPs1(1,NZ)=CPGRs1(NZ)*WTLFBs1(1,NZ)
  WTLSBs1(1,NZ)=WTLFBs1(1,NZ)+WTSHEBs1(1,NZ)
  WTLSs1(NZ)=WTLSs1(NZ)+WTLSBs1(1,NZ)
  FDM=AMIN1(1.0,0.16-0.045*PSILTs1(NZ))
  VOLWPs1(NZ)=1.0E-06*WTLSs1(NZ)/FDM
  VOLWCs1(NZ)=0._r8
  ZPOOLs1(1,NZ)=CNGRs1(NZ)*CPOOLs1(1,NZ)
  PPOOLs1(1,NZ)=CPGRs1(NZ)*CPOOLs1(1,NZ)
  WTRT1Ns1(1,NGs1(NZ),1,NZ)=CNGRs1(NZ)*WTRT1s1(1,NGs1(NZ),1,NZ)
  WTRT1Ps1(1,NGs1(NZ),1,NZ)=CPGRs1(NZ)*WTRT1s1(1,NGs1(NZ),1,NZ)
  RTWT1Ns1(1,1,NZ)=CNGRs1(NZ)*RTWT1s1(1,1,NZ)
  RTWT1Ps1(1,1,NZ)=CPGRs1(NZ)*RTWT1s1(1,1,NZ)
  WTRTLs1(1,NGs1(NZ),NZ)=WTRT1s1(1,NGs1(NZ),1,NZ)
  WTRTDs1(1,NGs1(NZ),NZ)=WTRT1s1(1,NGs1(NZ),1,NZ)
  WSRTLs1(1,NGs1(NZ),NZ)=WTRTLs1(1,NGs1(NZ),NZ)*CWSRTs1(NZ)
  ZPOOLRs1(1,NGs1(NZ),NZ)=CNGRs1(NZ)*CPOOLRs1(1,NGs1(NZ),NZ)
  PPOOLRs1(1,NGs1(NZ),NZ)=CPGRs1(NZ)*CPOOLRs1(1,NGs1(NZ),NZ)

  end associate
  end subroutine InitSeedMorphoBio

  end module StartqsMod
