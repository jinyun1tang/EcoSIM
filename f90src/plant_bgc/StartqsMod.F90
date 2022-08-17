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
        TCSN0s1(NZ)=0._r8
        TZSN0s1(NZ)=0._r8
        TPSN0s1(NZ)=0._r8
        TCSNCs1(NZ)=0._r8
        TZSNCs1(NZ)=0._r8
        TPSNCs1(NZ)=0._r8
        WTSTGs1(NZ)=0._r8
        WTSTGNs1(NZ)=0._r8
        WTSTGPs1(NZ)=0._r8
        DO 6401 L=1,NLs1
          DO  K=0,1
            DO  M=1,4
              CSNCs1(M,K,L,NZ)=0._r8
              ZSNCs1(M,K,L,NZ)=0._r8
              PSNCs1(M,K,L,NZ)=0._r8
            enddo
          enddo
6401    CONTINUE
9986  CONTINUE
  RETURN
  END subroutine startqs
!------------------------------------------------------------------------------------------

  subroutine InitShootGrowth(NZ)

  implicit none
  integer, intent(in) :: NZ
  associate(                         &
    O2Is1    =>  plt_photo%O2Is1   , &
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

  associate(                         &
    NNODs1   =>  plt_morph%NNODs1    &
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
  end subroutine PFTThermalAcclimation
!------------------------------------------------------------------------------------------

  subroutine InitDimensionsandUptake(NZ)

  implicit none
  integer, intent(in) :: NZ
  INTEGER :: L,N,NR
  associate(                             &
    SDPTHs1    =>  plt_morph%SDPTHs1   , &
    SDARs1     =>  plt_morph%SDARs1    , &
    SDPTHIs1   =>  plt_morph%SDPTHIs1  , &
    NIXs1      =>  plt_morph%NIXs1     , &
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
    FDBKXs1   =>  plt_photo%FDBKXs1  , &
    FDBKs1    =>  plt_photo%FDBKs1   , &
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
    NBTBs1    =>  plt_morph%NBTBs1   , &
    HTNODXs1  =>  plt_morph%HTNODXs1 , &
    HTSHEs1   =>  plt_morph%HTSHEs1  , &
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
  IFLGIs1(NZ)=0
  IDTHPs1(NZ)=0
  IDTHRs1(NZ)=0
  NBTs1(NZ)=0
  NBRs1(NZ)=0
  HTCTLs1(NZ)=0._r8
  ZCs1(NZ)=0._r8
  DO 10 NB=1,10
    IFLGAs1(NB,NZ)=0
    IFLGEs1(NB,NZ)=0
    IFLGFs1(NB,NZ)=0
    IFLGRs1(NB,NZ)=0
    IFLGQs1(NB,NZ)=0
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
    IDTHBs1(NB,NZ)=1
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
    CPOOLs1(NB,NZ)=0._r8
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
    RNH3Bs1(NB,NZ)=0._r8
    RCZLXs1(NB,NZ)=0._r8
    RCPLXs1(NB,NZ)=0._r8
    RCCLXs1(NB,NZ)=0._r8
    WGLFXs1(NB,NZ)=0._r8
    WGLFNXs1(NB,NZ)=0._r8
    WGLFPXs1(NB,NZ)=0._r8
    ARLFZs1(NB,NZ)=0._r8
    RCZSXs1(NB,NZ)=0._r8
    RCPSXs1(NB,NZ)=0._r8
    RCCSXs1(NB,NZ)=0._r8
    WTSTXBs1(NB,NZ)=0._r8
    WTSTXNs1(NB,NZ)=0._r8
    WTSTXPs1(NB,NZ)=0._r8
    WGSHEXs1(NB,NZ)=0._r8
    WGSHNXs1(NB,NZ)=0._r8
    WGSHPXs1(NB,NZ)=0._r8
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
      DO 55 L=1,JC1
        ARLFLs1(L,K,NB,NZ)=0._r8
        WGLFLs1(L,K,NB,NZ)=0._r8
        WGLFLNs1(L,K,NB,NZ)=0._r8
        WGLFLPs1(L,K,NB,NZ)=0._r8
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
    WGLFVs1(L,NZ)=0._r8
    ARSTVs1(L,NZ)=0._r8
35  CONTINUE
  CPOOLPs1(NZ)=0._r8
  ZPOOLPs1(NZ)=0._r8
  PPOOLPs1(NZ)=0._r8
  CCPOLPs1(NZ)=0._r8
  CCPLNPs1(NZ)=0._r8
  CZPOLPs1(NZ)=0._r8
  CPPOLPs1(NZ)=0._r8
  WTSHTs1(NZ)=0._r8
  WTLFs1(NZ)=0._r8
  WTSHEs1(NZ)=0._r8
  WTSTKs1(NZ)=0._r8
  WVSTKs1(NZ)=0._r8
  WTRSVs1(NZ)=0._r8
  WTHSKs1(NZ)=0._r8
  WTEARs1(NZ)=0._r8
  WTGRs1(NZ)=0._r8
  WTRTts1(NZ)=0._r8
  WTRTSs1(NZ)=0._r8
  WTNDs1(NZ)=0._r8
  WTLSs1(NZ)=0._r8
  WTSHNs1(NZ)=0._r8
  WTLFNs1(NZ)=0._r8
  WTSHENs1(NZ)=0._r8
  WTSTKNs1(NZ)=0._r8
  WTRSVNs1(NZ)=0._r8
  WTHSKNs1(NZ)=0._r8
  WTEARNs1(NZ)=0._r8
  WTGRNNs1(NZ)=0._r8
  WTNDNs1(NZ)=0._r8
  WTSHPs1(NZ)=0._r8
  WTLFPs1(NZ)=0._r8
  WTSHEPs1(NZ)=0._r8
  WTSTKPs1(NZ)=0._r8
  WTRSVPs1(NZ)=0._r8
  WTHSKPs1(NZ)=0._r8
  WTEARPs1(NZ)=0._r8
  WTGRNPs1(NZ)=0._r8
  WTNDPs1(NZ)=0._r8
  ARLFPs1(NZ)=0._r8
  WTRTAs1(NZ)=0._r8
  ARSTPs1(NZ)=0._r8
  end associate
  end subroutine InitPlantPhenoMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitMassBalance(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: M
  real(r8) :: WTSTDX
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
    VCO2Fs1(NZ)=0._r8
    VCH4Fs1(NZ)=0._r8
    VOXYFs1(NZ)=0._r8
    VNH3Fs1(NZ)=0._r8
    VN2OFs1(NZ)=0._r8
    VPO4Fs1(NZ)=0._r8
    THVSTCs1(NZ)=0._r8
    THVSTNs1(NZ)=0._r8
    THVSTPs1(NZ)=0._r8
    HVSTCs1(NZ)=0._r8
    HVSTNs1(NZ)=0._r8
    HVSTPs1(NZ)=0._r8
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
  end subroutine InitMassBalance
!------------------------------------------------------------------------------------------

  subroutine InitPlantHeatandWater(NZ)

  implicit none
  integer, intent(in) :: NZ
  associate(                          &
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
    SDPTHs1    =>  plt_morph%SDPTHs1   , &
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
  UPNH4s1(NZ)=0._r8
  UPNO3s1(NZ)=0._r8
  UPH2Ps1(NZ)=0._r8
  UPH1Ps1(NZ)=0._r8
  UPNFs1(NZ)=0._r8
  DO 40 N=1,2
    DO 20 L=1,NLs1
      UPWTRs1(N,L,NZ)=0._r8
      PSIRTs1(N,L,NZ)=-0.01
      PSIROs1(N,L,NZ)=OSMOs1(NZ)+PSIRTs1(N,L,NZ)
      PSIRGs1(N,L,NZ)=AMAX1(0.0_r8,PSIRTs1(N,L,NZ)-PSIROs1(N,L,NZ))
      CPOOLRs1(N,L,NZ)=0._r8
      ZPOOLRs1(N,L,NZ)=0._r8
      PPOOLRs1(N,L,NZ)=0._r8
      CCPOLRs1(N,L,NZ)=0._r8
      CZPOLRs1(N,L,NZ)=0._r8
      CPPOLRs1(N,L,NZ)=0._r8
      CWSRTLs1(N,L,NZ)=CWSRTs1(NZ)
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
      RTLGAs1(N,L,NZ)=1.0E-03
      RUPNH4s1(N,L,NZ)=0._r8
      RUPNO3s1(N,L,NZ)=0._r8
      RUPH2Ps1(N,L,NZ)=0._r8
      RUPH1Ps1(N,L,NZ)=0._r8
      RUPNHBs1(N,L,NZ)=0._r8
      RUPNOBs1(N,L,NZ)=0._r8
      RUPH2Bs1(N,L,NZ)=0._r8
      RUPH1Bs1(N,L,NZ)=0._r8
      ROXYPs1(N,L,NZ)=0._r8
      RUNNHPs1(N,L,NZ)=0._r8
      RUNNBPs1(N,L,NZ)=0._r8
      RUNNOPs1(N,L,NZ)=0._r8
      RUNNXPs1(N,L,NZ)=0._r8
      RUPP2Ps1(N,L,NZ)=0._r8
      RUPP1Ps1(N,L,NZ)=0._r8
      RUPP2Bs1(N,L,NZ)=0._r8
      RUPP1Bs1(N,L,NZ)=0._r8
      CCO2A=CCO2EIs1
      CCO2P=0.030*EXP(-2.621_r8-0.0317_r8*ATCAs1)*CO2EIs1
      CO2As1(N,L,NZ)=CCO2A*RTVLPs1(N,L,NZ)
      CO2Ps1(N,L,NZ)=CCO2P*RTVLWs1(N,L,NZ)
      RCOFLAs1(N,L,NZ)=0._r8
      RCODFAs1(N,L,NZ)=0._r8
      RCO2Ss1(N,L,NZ)=0._r8
      RCO2Ps1(N,L,NZ)=0._r8
      COXYA=COXYEs1
      COXYP=0.032*EXP(-6.175_r8-0.0211_r8*ATCAs1)*OXYEs1
      OXYAs1(N,L,NZ)=COXYA*RTVLPs1(N,L,NZ)
      OXYPs1(N,L,NZ)=COXYP*RTVLWs1(N,L,NZ)
      CH4As1(N,L,NZ)=0._r8
      CH4Ps1(N,L,NZ)=0._r8
      Z2OAs1(N,L,NZ)=0._r8
      Z2OPs1(N,L,NZ)=0._r8
      ZH3As1(N,L,NZ)=0._r8
      ZH3Ps1(N,L,NZ)=0._r8
      H2GAs1(N,L,NZ)=0._r8
      H2GPs1(N,L,NZ)=0._r8
      WFRs1(N,L,NZ)=1.0
      DO 30 NR=1,JC1
        RTN2s1(N,L,NR,NZ)=0._r8
        RTLG1s1(N,L,NR,NZ)=0._r8
        WTRT1s1(N,L,NR,NZ)=0._r8
        WTRT1Ns1(N,L,NR,NZ)=0._r8
        WTRT1Ps1(N,L,NR,NZ)=0._r8
        RTLG2s1(N,L,NR,NZ)=0._r8
        WTRT2s1(N,L,NR,NZ)=0._r8
        WTRT2Ns1(N,L,NR,NZ)=0._r8
        WTRT2Ps1(N,L,NR,NZ)=0._r8
        RTDP1s1(N,NR,NZ)=SDPTHs1(NZ)
        RTWT1s1(N,NR,NZ)=0._r8
        RTWT1Ns1(N,NR,NZ)=0._r8
        RTWT1Ps1(N,NR,NZ)=0._r8
30    CONTINUE
      IF(N.EQ.1)THEN
        DO 6400 K=0,1
          DO  M=1,4
            CSNCs1(M,K,L,NZ)=0._r8
            ZSNCs1(M,K,L,NZ)=0._r8
            PSNCs1(M,K,L,NZ)=0._r8
          enddo
6400    CONTINUE
        CPOOLNs1(L,NZ)=0._r8
        ZPOOLNs1(L,NZ)=0._r8
        PPOOLNs1(L,NZ)=0._r8
        WTNDLs1(L,NZ)=0._r8
        WTNDLNs1(L,NZ)=0._r8
        WTNDLPs1(L,NZ)=0._r8
        RUPNFs1(L,NZ)=0._r8
      ENDIF
20  CONTINUE
40  CONTINUE

  RUPNH4s1(1:2,NLs1+1:JZ1,NZ)=0._r8
  RUPNHBs1(1:2,NLs1+1:JZ1,NZ)=0._r8
  RUPH2Ps1(1:2,NLs1+1:JZ1,NZ)=0._r8
  RUPH2Bs1(1:2,NLs1+1:JZ1,NZ)=0._r8
  RTDNPs1(1:2,NLs1+1:JZ1,NZ)=0._r8
  end associate
  end subroutine InitRootMychorMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitSeedMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  REAL(R8) :: FDM
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
  end subroutine InitSeedMorphoBio

  end module StartqsMod
