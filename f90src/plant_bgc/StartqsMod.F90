module StartqsMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use minimathmod, only : AZMAX1
  use EcoSIMConfig
  use PlantAPIData
  use MicBGCPars, only : micpar
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = __FILE__
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
    PP      => plt_site%PP        , &
    NU      => plt_site%NU        , &
    NP      => plt_site%NP        , &
    NL      => plt_site%NL        , &
    ZERO    => plt_site%ZERO      , &
    AREA3   => plt_site%AREA3     , &
    ZEROQ   => plt_rbgc%ZEROQ     , &
    ZEROP   => plt_biom%ZEROP     , &
    ZEROL   => plt_biom%ZEROL     , &
    IFLGC   => plt_pheno%IFLGC      &
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

      NZ2X=MIN(NZ2Q,NP)
      D9985: DO NZ=NZ1Q,NZ2X

        IF(IFLGC(NZ).EQ.0)THEN

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
        ZEROP(NZ)=ZERO*PP(NZ)
        ZEROQ(NZ)=ZERO*PP(NZ)/AREA3(NU)
        ZEROL(NZ)=ZERO*PP(NZ)*1.0E+06_r8
      ENDDO D9985
!
!     FILL OUT UNUSED ARRAYS
!
      D9986: DO NZ=NP+1,JP1
        plt_bgcr%TESN0(1:npelms,NZ)=0._r8
        plt_bgcr%TESNC(1:npelms,NZ)=0._r8
        plt_biom%WTSTGE(1:npelms,NZ)=0._r8
        D6401: DO L=1,NL
          DO  K=0,micpar%n_pltlitrk
            plt_bgcr%ESNC(1:jsken,1:npelms,K,L,NZ)=0._r8
          enddo
        ENDDO D6401
      ENDDO D9986
  RETURN
  end associate
  END subroutine startqs
!------------------------------------------------------------------------------------------

  subroutine InitShootGrowth(NZ)

  implicit none
  integer, intent(in) :: NZ
  associate(                         &
    IDAYX  =>  plt_distb%IDAYX , &
    IYRY   =>  plt_distb%IYRY  , &
    IYRX   =>  plt_distb%IYRX  , &
    IYR0   =>  plt_distb%IYR0  , &
    IYRH   =>  plt_distb%IYRH  , &
    IDAYH  =>  plt_distb%IDAYH , &
    IDAY0  =>  plt_distb%IDAY0 , &
    IDAYY  =>  plt_distb%IDAYY , &
    RSMX   =>  plt_photo%RSMX  , &
    PPI    =>  plt_site%PPI    , &
    PPX    =>  plt_site%PPX    , &
    PPZ    =>  plt_site%PPZ    , &
    CWSRT  =>  plt_allom%CWSRT , &
    CNWS   =>  plt_allom%CNWS  , &
    CPWS   =>  plt_allom%CPWS  , &
    CNRT   =>  plt_allom%CNRT  , &
    CPRT   =>  plt_allom%CPRT  , &
    O2I    =>  plt_photo%O2I   , &
    RCMX   =>  plt_photo%RCMX  , &
    ICTYP  =>  plt_photo%ICTYP , &
    RSMH   =>  plt_photo%RSMH  , &
    CFI    =>  plt_morph%CFI   , &
    CF     =>  plt_morph%CF    , &
    NRT    =>  plt_morph%NRT     &
  )
  IYR0(NZ)=IYRX(NZ)
  IDAY0(NZ)=IDAYX(NZ)
  IYRH(NZ)=IYRY(NZ)
  IDAYH(NZ)=IDAYY(NZ)
  PPI(NZ)=PPZ(NZ)
  PPX(NZ)=PPI(NZ)
  CF(NZ)=CFI(NZ)

  RSMH(NZ)=RSMX(NZ)/3600.0_r8
  RCMX(NZ)=RSMX(NZ)*1.56_r8
  CNWS(NZ)=2.5_r8
  CPWS(NZ)=25.0_r8
  CWSRT(NZ)=AMIN1(CNRT(NZ)*CNWS(NZ),CPRT(NZ)*CPWS(NZ))
  IF(ICTYP(NZ).EQ.3)THEN
    O2I(NZ)=2.10E+05_r8
  ELSE
    O2I(NZ)=3.96E+05_r8
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
    XRLA   =>  plt_pheno%XRLA      , &
    IBTYP  =>  plt_pheno%IBTYP     , &
    IGTYP  =>  plt_pheno%IGTYP     , &
    GROUPI =>  plt_pheno%GROUPI    , &
    CFOPC  =>  plt_soilchem%CFOPC  , &
    CFOPN  =>  plt_soilchem%CFOPN  , &
    CFOPP  =>  plt_soilchem%CFOPP  , &
    FNOD   =>  plt_allom%FNOD      , &
    INTYP  =>  plt_morph%INTYP     , &
    NNOD   =>  plt_morph%NNOD        &
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
  CFOPC(0,1,NZ)=0.0_r8
  CFOPC(0,2,NZ)=0.67_r8
  CFOPC(0,3,NZ)=0.33_r8
  CFOPC(0,4,NZ)=0.0_r8
!
!     NON-VASCULAR (E.G. MOSSES)
!
  IF(IGTYP(NZ).EQ.0)THEN
    CFOPC(1,1,NZ)=0.07_r8
    CFOPC(1,2,NZ)=0.25_r8
    CFOPC(1,3,NZ)=0.30_r8
    CFOPC(1,4,NZ)=0.38_r8
    CFOPC(2,1,NZ)=0.07_r8
    CFOPC(2,2,NZ)=0.25_r8
    CFOPC(2,3,NZ)=0.30_r8
    CFOPC(2,4,NZ)=0.38_r8
!
!     LEGUMES
!
  ELSEIF(INTYP(NZ).NE.0)THEN
    CFOPC(1,1,NZ)=0.16_r8
    CFOPC(1,2,NZ)=0.38_r8
    CFOPC(1,3,NZ)=0.34_r8
    CFOPC(1,4,NZ)=0.12_r8
    CFOPC(2,1,NZ)=0.07_r8
    CFOPC(2,2,NZ)=0.41_r8
    CFOPC(2,3,NZ)=0.37_r8
    CFOPC(2,4,NZ)=0.15_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
    CFOPC(1,1,NZ)=0.08_r8
    CFOPC(1,2,NZ)=0.41_r8
    CFOPC(1,3,NZ)=0.36_r8
    CFOPC(1,4,NZ)=0.15_r8
    CFOPC(2,1,NZ)=0.07_r8
    CFOPC(2,2,NZ)=0.41_r8
    CFOPC(2,3,NZ)=0.36_r8
    CFOPC(2,4,NZ)=0.16_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(IBTYP(NZ).EQ.1.OR.IBTYP(NZ).EQ.3)THEN
    CFOPC(1,1,NZ)=0.07_r8
    CFOPC(1,2,NZ)=0.34_r8
    CFOPC(1,3,NZ)=0.36_r8
    CFOPC(1,4,NZ)=0.23_r8
    CFOPC(2,1,NZ)=0.0_r8
    CFOPC(2,2,NZ)=0.045_r8
    CFOPC(2,3,NZ)=0.660_r8
    CFOPC(2,4,NZ)=0.295_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPC(1,1,NZ)=0.07_r8
    CFOPC(1,2,NZ)=0.25_r8
    CFOPC(1,3,NZ)=0.38_r8
    CFOPC(1,4,NZ)=0.30_r8
    CFOPC(2,1,NZ)=0.0_r8
    CFOPC(2,2,NZ)=0.045_r8
    CFOPC(2,3,NZ)=0.660_r8
    CFOPC(2,4,NZ)=0.295_r8
  ENDIF
!
!     FRACTIONS OF WOODY LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     NON-VASCULAR
!
  IF(IGTYP(NZ).EQ.0)THEN
    CFOPC(3,1,NZ)=0.07_r8
    CFOPC(3,2,NZ)=0.25_r8
    CFOPC(3,3,NZ)=0.30_r8
    CFOPC(3,4,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
    CFOPC(3,1,NZ)=0.03_r8
    CFOPC(3,2,NZ)=0.25_r8
    CFOPC(3,3,NZ)=0.57_r8
    CFOPC(3,4,NZ)=0.15_r8
!
!     DECIDUOUS AND CONIFEROUS TREES
!
  ELSE
    CFOPC(3,1,NZ)=0.0_r8
    CFOPC(3,2,NZ)=0.045_r8
    CFOPC(3,3,NZ)=0.660_r8
    CFOPC(3,4,NZ)=0.295_r8
  ENDIF
!
!     FRACTIONS OF FINE ROOT LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN PC&E 25:601-608
!
!     NON-VASCULAR
!
  IF(IGTYP(NZ).EQ.0)THEN
    CFOPC(4,1,NZ)=0.07_r8
    CFOPC(4,2,NZ)=0.25_r8
    CFOPC(4,3,NZ)=0.30_r8
    CFOPC(4,4,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
    CFOPC(4,1,NZ)=0.057_r8
    CFOPC(4,2,NZ)=0.263_r8
    CFOPC(4,3,NZ)=0.542_r8
    CFOPC(4,4,NZ)=0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(IBTYP(NZ).EQ.1.OR.IBTYP(NZ).EQ.3)THEN
    CFOPC(4,1,NZ)=0.059_r8
    CFOPC(4,2,NZ)=0.308_r8
    CFOPC(4,3,NZ)=0.464_r8
    CFOPC(4,4,NZ)=0.169_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPC(4,1,NZ)=0.059_r8
    CFOPC(4,2,NZ)=0.308_r8
    CFOPC(4,3,NZ)=0.464_r8
    CFOPC(4,4,NZ)=0.169_r8
  ENDIF
!
!     COARSE WOODY LITTER FROM BOLES AND ROOTS
!
  CFOPC(5,1,NZ)=0.00_r8
  CFOPC(5,2,NZ)=0.045_r8
  CFOPC(5,3,NZ)=0.660_r8
  CFOPC(5,4,NZ)=0.295_r8
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
      CNOPCT=CNOPCT+CFOPC(N,M,NZ)*CNOPC(M)
      CPOPCT=CPOPCT+CFOPC(N,M,NZ)*CPOPC(M)
100 CONTINUE
    DO 105 M=1,4
      CFOPN(N,M,NZ)=CFOPC(N,M,NZ)*CNOPC(M)/CNOPCT
      CFOPP(N,M,NZ)=CFOPC(N,M,NZ)*CPOPC(M)/CPOPCT
105 CONTINUE
110 CONTINUE
!
!     CONCURRENT NODE GROWTH
!
!     FNOD=scales node number for perennial vegetation (e.g. trees)
!     NNOD=number of concurrently growing nodes
!
  IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
    FNOD(NZ)=1.0
    IF(GROUPI(NZ).LE.10)THEN
      NNOD(NZ)=3
    ELSEIF(GROUPI(NZ).LE.15)THEN
      NNOD(NZ)=4
    ELSE
      NNOD(NZ)=5
    ENDIF
  ELSE
    FNOD(NZ)=AMAX1(1.0,0.04/XRLA(NZ))
    NNOD(NZ)=24
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
    DATAP    =>  plt_site%DATAP   , &
    ICTYP    =>  plt_photo%ICTYP  , &
    HTC      =>  plt_pheno%HTC    , &
    TCX      =>  plt_pheno%TCX    , &
    TCZ      =>  plt_pheno%TCZ    , &
    OFFST    =>  plt_pheno%OFFST  , &
    ZTYPI    =>  plt_pheno%ZTYPI  , &
    ZTYP     =>  plt_pheno%ZTYP   , &
    SSTX     =>  plt_pheno%SSTX     &
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
  ZTYP(NZ)=ZTYPI(NZ)
  OFFST(NZ)=2.667_r8*(2.5_r8-ZTYP(NZ))
  TCZ(NZ)=TCZD-OFFST(NZ)
  TCX(NZ)=AMIN1(15.0_r8,TCXD-OFFST(NZ))
  IF(ICTYP(NZ).EQ.3)THEN
    IF(DATAP(NZ)(1:4).EQ.'soyb')THEN
      HTC(NZ)=30.0_r8+3.0_r8*ZTYP(NZ)
      SSTX(NZ)=0.002_r8
    ELSE
      HTC(NZ)=27.0_r8+3.0_r8*ZTYP(NZ)
      SSTX(NZ)=0.002_r8
    ENDIF
  ELSE
    HTC(NZ)=27.0_r8+3.0_r8*ZTYP(NZ)
    SSTX(NZ)=0.005_r8
  ENDIF
  end associate
  end subroutine PFTThermalAcclimation
!------------------------------------------------------------------------------------------

  subroutine InitDimensionsandUptake(NZ)

  implicit none
  integer, intent(in) :: NZ
  INTEGER :: L,N,NR
  associate(                             &
    CNRTS    =>  plt_allom%CNRTS   , &
    CPRTS    =>  plt_allom%CPRTS   , &
    DMRT     =>  plt_allom%DMRT    , &
    CNRT     =>  plt_allom%CNRT    , &
    CPRT     =>  plt_allom%CPRT    , &
    UPMNPO   =>  plt_rbgc%UPMNPO   , &
    UPKMPO   =>  plt_rbgc%UPKMPO   , &
    UPMXPO   =>  plt_rbgc%UPMXPO   , &
    UPMNZO   =>  plt_rbgc%UPMNZO   , &
    UPKMZO   =>  plt_rbgc%UPKMZO   , &
    UPMXZO   =>  plt_rbgc%UPMXZO   , &
    UPMNZH   =>  plt_rbgc%UPMNZH   , &
    UPMXZH   =>  plt_rbgc%UPMXZH   , &
    UPKMZH   =>  plt_rbgc%UPKMZH   , &
    CDPTHZ   =>  plt_site%CDPTHZ   , &
    NU       =>  plt_site%NU       , &
    NL       =>  plt_site%NL       , &
    NG       =>  plt_morph%NG      , &
    SDPTH    =>  plt_morph%SDPTH   , &
    SDAR     =>  plt_morph%SDAR    , &
    GRDM     =>  plt_morph%GRDM    , &
    SDPTHI   =>  plt_morph%SDPTHI  , &
    RRAD2X   =>  plt_morph%RRAD2X  , &
    RTAR2X   =>  plt_morph%RTAR2X  , &
    RRAD1X   =>  plt_morph%RRAD1X  , &
    RTAR1X   =>  plt_morph%RTAR1X  , &
    RRAD2M   =>  plt_morph%RRAD2M  , &
    RRAD1M   =>  plt_morph%RRAD1M  , &
    RTLG2X   =>  plt_morph%RTLG2X  , &
    NIX      =>  plt_morph%NIX     , &
    RSRR     =>  plt_morph%RSRR    , &
    RTLG1X   =>  plt_morph%RTLG1X  , &
    PORT     =>  plt_morph%PORT    , &
    PORTX    =>  plt_morph%PORTX   , &
    RRADP    =>  plt_morph%RRADP   , &
    DMVL     =>  plt_morph%DMVL    , &
    RSRA     =>  plt_morph%RSRA    , &
    NINR     =>  plt_morph%NINR    , &
    SDVL     =>  plt_morph%SDVL    , &
    SDLG     =>  plt_morph%SDLG      &
  )
!
!     SEED CHARACTERISTICS
!
!     SDVL,SDLG,SDAR=seed volume(m3),length(m),AREA3(NU)(m2)
!     GRDM=seed C mass (g) from PFT file
!
  SDVL(NZ)=GRDM(NZ)*5.0E-06_r8
  SDLG(NZ)=2.0_r8*(0.75*SDVL(NZ)/PICON)**0.33_r8
  SDAR(NZ)=4.0_r8*PICON*(SDLG(NZ)/2.0_r8)**2_r8
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
  SDPTH(NZ)=SDPTHI(NZ)
  D9795: DO L=NU,NL
    IF(SDPTH(NZ).GE.CDPTHZ(L-1).AND.SDPTH(NZ).LT.CDPTHZ(L))THEN
      NG(NZ)=L
      NIX(NZ)=L
      D9790: DO NR=1,JC1
        NINR(NR,NZ)=L
      ENDDO D9790
    ENDIF
  ENDDO D9795
  CNRTS(NZ)=CNRT(NZ)*DMRT(NZ)
  CPRTS(NZ)=CPRT(NZ)*DMRT(NZ)
  RRAD1M(2,NZ)=5.0E-06_r8
  RRAD2M(2,NZ)=5.0E-06_r8
  PORT(2,NZ)=PORT(1,NZ)
  UPMXZH(2,NZ)=UPMXZH(1,NZ)
  UPKMZH(2,NZ)=UPKMZH(1,NZ)
  UPMNZH(2,NZ)=UPMNZH(1,NZ)
  UPMXZO(2,NZ)=UPMXZO(1,NZ)
  UPKMZO(2,NZ)=UPKMZO(1,NZ)
  UPMNZO(2,NZ)=UPMNZO(1,NZ)
  UPMXPO(2,NZ)=UPMXPO(1,NZ)
  UPKMPO(2,NZ)=UPKMPO(1,NZ)
  UPMNPO(2,NZ)=UPMNPO(1,NZ)
  RSRR(2,NZ)=1.0E+04_r8
  RSRA(2,NZ)=1.0E+12_r8
!
!     PORTX=tortuosity for gas transport
!     RRADP=path length for radial diffusion within root (m)
!     DMVL=volume:C ratio (m3 g-1)
!     RTLG1X,RTLG2X=specific primary,secondary root length (m g-1)
!     RTAR1X,RTAR2X=specific primary,secondary root area (m2 g-1)
!
  DO 500 N=1,2
    PORTX(N,NZ)=PORT(N,NZ)**1.33_r8
    RRADP(N,NZ)=LOG(1.0_r8/SQRT(AMAX1(0.01_r8,PORT(N,NZ))))
    DMVL(N,NZ)=ppmc/(0.05_r8*(1.0-PORT(N,NZ)))
    RTLG1X(N,NZ)=DMVL(N,NZ)/(PICON*RRAD1M(N,NZ)**2)
    RTLG2X(N,NZ)=DMVL(N,NZ)/(PICON*RRAD2M(N,NZ)**2)
    RRAD1X(N,NZ)=RRAD1M(N,NZ)
!    2*SQRT(0.25*(1.0-PORT(N,NZ)))
    RRAD2X(N,NZ)=RRAD2M(N,NZ)
!    2*SQRT(0.25*(1.0-PORT(N,NZ)))
    RTAR1X(N,NZ)=PICON*RRAD1X(N,NZ)**2
    RTAR2X(N,NZ)=PICON*RRAD2X(N,NZ)**2
500 CONTINUE
  end associate
  end subroutine InitDimensionsandUptake
!------------------------------------------------------------------------------------------

  subroutine InitPlantPhenoMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: K,L,M,N,NB
  associate(                           &
    NU      =>  plt_site%NU      , &
    PPX     =>  plt_site%PPX     , &
    PP      =>  plt_site%PP      , &
    ALAT    =>  plt_site%ALAT    , &
    AREA3   =>  plt_site%AREA3   , &
    IGTYP   =>  plt_pheno%IGTYP  , &
    VSTGX   =>  plt_pheno%VSTGX  , &
    GROUP   =>  plt_pheno%GROUP  , &
    KVSTG   =>  plt_pheno%KVSTG  , &
    KVSTGN  =>  plt_pheno%KVSTGN , &
    GSTGI   =>  plt_pheno%GSTGI  , &
    GSTGF   =>  plt_pheno%GSTGF  , &
    FLG4    =>  plt_pheno%FLG4   , &
    FLGZ    =>  plt_pheno%FLGZ   , &
    VRNY    =>  plt_pheno%VRNY   , &
    VRNF    =>  plt_pheno%VRNF   , &
    IDTHB   =>  plt_pheno%IDTHB  , &
    VRNZ    =>  plt_pheno%VRNZ   , &
    IDAY    =>  plt_pheno%IDAY   , &
    VRNS    =>  plt_pheno%VRNS   , &
    ATRP    =>  plt_pheno%ATRP   , &
    WSTR    =>  plt_pheno%WSTR   , &
    GROUPI  =>  plt_pheno%GROUPI , &
    RCESX   =>  plt_pheno%RCESX  , &
    TGSTGF  =>  plt_pheno%TGSTGF , &
    TGSTGI  =>  plt_pheno%TGSTGI , &
    FDBKX   =>  plt_photo%FDBKX  , &
    CPOOL3  =>  plt_photo%CPOOL3 , &
    CPOOL4  =>  plt_photo%CPOOL4 , &
    CHILL   =>  plt_photo%CHILL  , &
    FDBK    =>  plt_photo%FDBK   , &
    HCOB    =>  plt_photo%HCOB   , &
    CO2B    =>  plt_photo%CO2B   , &
    VSTG    =>  plt_morph%VSTG   , &
    ZC      =>  plt_morph%ZC     , &
    KLEAF   =>  plt_morph%KLEAF  , &
    XTLI    =>  plt_morph%XTLI   , &
    NBT     =>  plt_morph%NBT    , &
    PSTG    =>  plt_morph%PSTG   , &
    ARLFP   =>  plt_morph%ARLFP  , &
    ARSTK   =>  plt_morph%ARSTK  , &
    SURFB   =>  plt_morph%SURFB  , &
    ARSTP   =>  plt_morph%ARSTP  , &
    PSTGF   =>  plt_morph%PSTGF  , &
    GRNXB   =>  plt_morph%GRNXB  , &
    HTNODE  =>  plt_morph%HTNODE , &
    GRNOB   =>  plt_morph%GRNOB  , &
    ARLFB   =>  plt_morph%ARLFB  , &
    ARLFZ   =>  plt_morph%ARLFZ  , &
    ARLFL   =>  plt_morph%ARLFL  , &
    ARLFV   =>  plt_morph%ARLFV  , &
    ARSTV   =>  plt_morph%ARSTV  , &
    SURF    =>  plt_morph%SURF   , &
    ARLF1   =>  plt_morph%ARLF1  , &
    HTSHEX  =>  plt_morph%HTSHEX , &
    HTCTL   =>  plt_morph%HTCTL  , &
    NBTB    =>  plt_morph%NBTB   , &
    PSTGI   =>  plt_morph%PSTGI  , &
    KLEAFX  =>  plt_morph%KLEAFX , &
    NBR     =>  plt_morph%NBR      &
  )
!
!     INITIALIZE PLANT PHENOLOGY
!
!     PP=population (grid cell-1)
!
  PP(NZ)=PPX(NZ)*AREA3(NU)
  plt_pheno%IFLGI(NZ)=0
  plt_pheno%IDTHP(NZ)=0
  plt_pheno%IDTHR(NZ)=0
  NBT(NZ)=0
  NBR(NZ)=0
  HTCTL(NZ)=0._r8
  ZC(NZ)=0._r8
  D10: DO NB=1,JC1
    plt_pheno%IFLGA(NB,NZ)=0
    plt_pheno%IFLGE(NB,NZ)=0
    plt_pheno%IFLGF(NB,NZ)=0
    plt_pheno%IFLGR(NB,NZ)=0
    plt_pheno%IFLGQ(NB,NZ)=0
    GROUP(NB,NZ)=GROUPI(NZ)
    PSTG(NB,NZ)=XTLI(NZ)
    PSTGI(NB,NZ)=PSTG(NB,NZ)
    PSTGF(NB,NZ)=0._r8
    VSTG(NB,NZ)=0._r8
    VSTGX(NB,NZ)=0._r8
    KLEAF(NB,NZ)=1
    KLEAFX(NB,NZ)=1
    KVSTG(NB,NZ)=1
    KVSTGN(NB,NZ)=0
    GSTGI(NB,NZ)=0._r8
    GSTGF(NB,NZ)=0._r8
    TGSTGI(NB,NZ)=0._r8
    TGSTGF(NB,NZ)=0._r8
    VRNY(NB,NZ)=0._r8
    VRNZ(NB,NZ)=0._r8
    VRNS(NB,NZ)=VRNY(NB,NZ)
    VRNF(NB,NZ)=VRNZ(NB,NZ)
    ATRP(NB,NZ)=0._r8
    FDBK(NB,NZ)=1.0
    FDBKX(NB,NZ)=1.0
    FLG4(NB,NZ)=0
    FLGZ(NB,NZ)=0
    NBTB(NB,NZ)=0
    plt_pheno%IDTHB(NB,NZ)=1
    D15: DO M=1,10
      IDAY(M,NB,NZ)=0
    ENDDO D15
  ENDDO D10
!
!     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
!
  WSTR(NZ)=0._r8
  CHILL(NZ)=0._r8
  plt_biom%EPOOL(1:JC1,1:npelms,NZ)=0._r8
  plt_biom%EPOLNB(1:JC1,1:npelms,NZ)=0._r8
  plt_biom%WTSHEBE(1:JC1,1:npelms,NZ)=0._r8
  plt_biom%WTSHTBE(1:JC1,1:npelms,NZ)=0._r8
  plt_biom%WTSTKBE(1:JC1,1:npelms,NZ)=0._r8
  plt_biom%WTLFBE(1:JC1,1:npelms,NZ)=0._r8
  plt_biom%WTRSVBE(1:JC1,1:npelms,NZ)=0._r8
  plt_biom%WTHSKBE(1:JC1,1:npelms,NZ)=0._r8
  plt_biom%WTGRBE(1:JC1,1:npelms,NZ)=0._r8
  plt_biom%WTEARBE(1:JC1,1:npelms,NZ)=0._r8
  plt_biom%WTNDBE(1:JC1,1:npelms,NZ)=0._r8
  D25: DO NB=1,JC1
    plt_biom%WVSTKB(NB,NZ)=0._r8
    plt_biom%WTLSB(NB,NZ)=0._r8
    GRNXB(NB,NZ)=0._r8
    GRNOB(NB,NZ)=0._r8
    plt_allom%GRWTB(NB,NZ)=0._r8
    ARLFB(NB,NZ)=0._r8
    plt_rbgc%RNH3B(NB,NZ)=0._r8
    plt_pheno%RCELX(NB,1:npelms,NZ)=0._r8
    plt_biom%WGLFEX(NB,1:npelms,NZ)=0._r8
    ARLFZ(NB,NZ)=0._r8
    plt_pheno%RCESX(NB,1:npelms,NZ)=0._r8
    plt_biom%WTSTXB(NB,NZ)=0._r8
    plt_biom%WTSTXN(NB,NZ)=0._r8
    plt_biom%WTSTXP(NB,NZ)=0._r8
    plt_biom%WGSHEXE(NB,1:npelms,NZ)=0._r8
    HTSHEX(NB,NZ)=0._r8
    D5: DO L=1,JC1
      ARSTK(L,NB,NZ)=0._r8
      DO N=1,JLI1
        SURFB(N,L,NB,NZ)=0._r8
      enddo
    ENDDO D5
    DO K=0,JNODS1
      ARLF1(K,NB,NZ)=0._r8
      HTNODE(K,NB,NZ)=0._r8
      plt_morph%HTNODX(K,NB,NZ)=0._r8
      plt_morph%HTSHE(K,NB,NZ)=0._r8
      plt_biom%WGLFE(K,NB,1:npelms,NZ)=0._r8
      plt_biom%WSLF(K,NB,NZ)=0._r8
      plt_biom%WGSHE(K,NB,1:npelms,NZ)=0._r8
      plt_biom%WSSHE(K,NB,NZ)=0._r8
      plt_biom%WGNODE(K,NB,NZ)=0._r8
      plt_biom%WGNODN(K,NB,NZ)=0._r8
      plt_biom%WGNODP(K,NB,NZ)=0._r8
      D55: DO L=1,JC1
        ARLFL(L,K,NB,NZ)=0._r8
        plt_biom%WGLFL(L,K,NB,NZ)=0._r8
        plt_biom%WGLFLN(L,K,NB,NZ)=0._r8
        plt_biom%WGLFLP(L,K,NB,NZ)=0._r8
      ENDDO D55
      IF(K.NE.0)THEN
        CPOOL3(K,NB,NZ)=0._r8
        CO2B(K,NB,NZ)=0._r8
        HCOB(K,NB,NZ)=0._r8
        CPOOL4(K,NB,NZ)=0._r8
        D45: DO L=1,JC1
          DO N=1,JLI1
            SURF(N,L,K,NB,NZ)=0._r8
          enddo
        ENDDO D45
      ENDIF
    enddo
  ENDDO D25
  D35: DO L=1,JC1
    ARLFV(L,NZ)=0._r8
    plt_biom%WGLFV(L,NZ)=0._r8
    ARSTV(L,NZ)=0._r8
  ENDDO D35
  plt_biom%EPOOLP(1:npelms,NZ)=0._r8
  plt_biom%CEPOLP(1:npelms,NZ)=0._r8
  plt_biom%CCPLNP(NZ)=0._r8
  plt_biom%WTSHTE(1:npelms,NZ)=0._r8
  plt_biom%WTLFE(1:npelms,NZ)=0._r8
  plt_biom%WTSHEE(1:npelms,NZ)=0._r8
  plt_biom%WTSTKE(1:npelms,NZ)=0._r8
  plt_biom%WVSTK(NZ)=0._r8
  plt_biom%WTRSVE(1:npelms,NZ)=0._r8
  plt_biom%WTHSKE(1:npelms,NZ)=0._r8
  plt_biom%WTEARE(1:npelms,NZ)=0._r8
  plt_biom%WTGRE(1:npelms,NZ)=0._r8
  plt_biom%WTRTSE(1:npelms,NZ)=0._r8
  plt_biom%WTNDE(1:npelms,NZ)=0._r8
  plt_biom%WTLS(NZ)=0._r8
  ARLFP(NZ)=0._r8
  plt_biom%WTRTA(NZ)=0._r8
  ARSTP(NZ)=0._r8
  end associate
  end subroutine InitPlantPhenoMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitMassBalance(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: M
  real(r8) :: WTSTDX

  associate(                           &
    NU     => plt_site%NU        , &
    AREA3  => plt_site%AREA3     , &
    CFOPC  => plt_soilchem%CFOPC , &
    CFOPN  => plt_soilchem%CFOPN , &
    CFOPP  => plt_soilchem%CFOPP , &
    CTRAN  => plt_ew%CTRAN       , &
    WTSTDE => plt_biom%WTSTDE    , &
    WTSTGE => plt_biom%WTSTGE    , &
    WTSTDI => plt_biom%WTSTDI    , &
    WTRTA  => plt_biom%WTRTA     , &
    CNSTK  => plt_allom%CNSTK    , &
    CPSTK  => plt_allom%CPSTK    , &
    TEUPTK => plt_rbgc%TEUPTK    , &
    TNH3C  => plt_bgcr%TNH3C     , &
    RNH3C  => plt_bgcr%RNH3C     , &
    TESN0  => plt_bgcr%TESN0     , &
    CARBN  => plt_bgcr%CARBN     , &
    TCO2A  => plt_bgcr%TCO2A     , &
    TCO2T  => plt_bgcr%TCO2T     , &
    TZUPFX => plt_bgcr%TZUPFX    , &
    TESNC  => plt_bgcr%TESNC     , &
    ARSTP  => plt_morph%ARSTP    , &
    RSETE  => plt_pheno%RSETE      &

  )
!
!     INITIALIZE MASS BALANCE CHECKS
!
  IF(.not.is_restart_run.AND.is_first_year)THEN
    CARBN(NZ)=0._r8
    TESN0(1:npelms,NZ)=0._r8
    TCO2T(NZ)=0._r8
    TCO2A(NZ)=0._r8
    TEUPTK(1:npelms,NZ)=0._r8
    TESNC(1:npelms,NZ)=0._r8

    TZUPFX(NZ)=0._r8
    RNH3C(NZ)=0._r8
    TNH3C(NZ)=0._r8
    plt_distb%VCO2F(NZ)=0._r8
    plt_distb%VCH4F(NZ)=0._r8
    plt_distb%VOXYF(NZ)=0._r8
    plt_distb%VNH3F(NZ)=0._r8
    plt_distb%VN2OF(NZ)=0._r8
    plt_distb%VPO4F(NZ)=0._r8
    plt_distb%THVSTE(1:npelms,NZ)=0._r8
    plt_distb%HVSTE(1:npelms,NZ)=0._r8
    RSETE(1:npelms,NZ)=0._r8
    CTRAN(NZ)=0._r8
    WTSTGE(1:npelms,NZ)=0._r8
    WTSTDX=WTSTDI(NZ)*AREA3(NU)
    D155: DO M=1,jsken
      WTSTDE(M,ielmc,NZ)=WTSTDX*CFOPC(5,M,NZ)
      WTSTDE(M,ielmn,NZ)=WTSTDX*CNSTK(NZ)*CFOPN(5,M,NZ)
      WTSTDE(M,ielmp,NZ)=WTSTDX*CPSTK(NZ)*CFOPP(5,M,NZ)
      WTSTGE(ielmc,NZ)=WTSTGE(ielmc,NZ)+WTSTDE(M,ielmc,NZ)
      WTSTGE(ielmn,NZ)=WTSTGE(ielmn,NZ)+WTSTDE(M,ielmn,NZ)
      WTSTGE(ielmp,NZ)=WTSTGE(ielmp,NZ)+WTSTDE(M,ielmp,NZ)
    ENDDO D155
  ENDIF
  end associate
  end subroutine InitMassBalance
!------------------------------------------------------------------------------------------

  subroutine InitPlantHeatandWater(NZ)

  implicit none
  integer, intent(in) :: NZ
  associate(                          &
    ATCA     =>  plt_site%ATCA  , &
    OSMO     =>  plt_ew%OSMO    , &
    TKC      =>  plt_ew%TKC     , &
    EP       =>  plt_ew%EP      , &
    VHCPC    =>  plt_ew%VHCPC   , &
    PSILT    =>  plt_ew%PSILT   , &
    PSILG    =>  plt_ew%PSILG   , &
    PSILO    =>  plt_ew%PSILO   , &
    ENGYX    =>  plt_ew%ENGYX   , &
    DTKC     =>  plt_ew%DTKC    , &
    TCC      =>  plt_ew%TCC     , &
    TKG      =>  plt_pheno%TKG  , &
    TCG      =>  plt_pheno%TCG  , &
    TFN3     =>  plt_pheno%TFN3 , &
    WTSHTE   =>  plt_biom%WTSHTE , &
    FRADP    =>  plt_rad%FRADP    &
  )
!
!     INITIALIZE PLANT HEAT AND WATER STATUS
!
!     VHCPC=canopy heat capacity (MJ m-3 K-1)
!     TCC,TKC=canopy temperature for growth (oC,K)
!     TCG,TKG=canopy temperature for phenology (oC,K)
!     PSILT,PSILO,PSILG=canopy total,osmotic,turgor water potl(MPa)
!
  VHCPC(NZ)=cpw*WTSHTE(ielmc,NZ)*10.0E-06
  ENGYX(NZ)=0._r8
  DTKC(NZ)=0._r8
  TCC(NZ)=ATCA
  TKC(NZ)=TCC(NZ)+TC2K
  TCG(NZ)=TCC(NZ)
  TKG(NZ)=TCG(NZ)+TC2K
  TFN3(NZ)=1.0
  PSILT(NZ)=-1.0E-03
  PSILO(NZ)=OSMO(NZ)+PSILT(NZ)
  PSILG(NZ)=AZMAX1(PSILT(NZ)-PSILO(NZ))
  EP(NZ)=0._r8
  FRADP(NZ)=0._r8
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
    OSMO     =>  plt_ew%OSMO       , &
    PSILT    =>  plt_ew%PSILT      , &
    PSIRT    =>  plt_ew%PSIRT      , &
    PSIRG    =>  plt_ew%PSIRG      , &
    PSIRO    =>  plt_ew%PSIRO      , &
    RUPNF    =>  plt_bgcr%RUPNF    , &
    CO2P     =>  plt_rbgc%CO2P     , &
    CO2A     =>  plt_rbgc%CO2A     , &
    OXYE     =>  plt_site%OXYE     , &
    COXYE    =>  plt_site%COXYE    , &
    CO2EI    =>  plt_site%CO2EI    , &
    NL       =>  plt_site%NL       , &
    ATCA     =>  plt_site%ATCA     , &
    CCO2EI   =>  plt_site%CCO2EI   , &
    CWSRT    =>  plt_allom%CWSRT   , &
    CWSRTL   =>  plt_biom%CWSRTL   , &
    WSRTL    =>  plt_biom%WSRTL    , &
    SDPTH    =>  plt_morph%SDPTH   , &
    RTVLW    =>  plt_morph%RTVLW   , &
    RTVLP    =>  plt_morph%RTVLP   , &
    RRAD2    =>  plt_morph%RRAD2   , &
    RRAD1    =>  plt_morph%RRAD1   , &
    RRAD1M   =>  plt_morph%RRAD1M  , &
    RRAD2M   =>  plt_morph%RRAD2M  , &
    NRT      =>  plt_morph%NRT       &
  )
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
!
!     PSIRT,PSIRO,PSIRG=root,myco total,osmotic,turgor water potl(MPa)
!     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
!     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
!
  NRT(NZ)=0
  plt_rbgc%UPNH4(NZ)=0._r8
  plt_rbgc%UPNO3(NZ)=0._r8
  plt_rbgc%UPH2P(NZ)=0._r8
  plt_rbgc%UPH1P(NZ)=0._r8
  plt_rbgc%UPNF(NZ)=0._r8
  D40: DO N=1,2
    D20: DO L=1,NL
      plt_ew%UPWTR(N,L,NZ)=0._r8
      PSIRT(N,L,NZ)=-0.01
      PSIRO(N,L,NZ)=OSMO(NZ)+PSIRT(N,L,NZ)
      PSIRG(N,L,NZ)=AZMAX1(PSIRT(N,L,NZ)-PSIRO(N,L,NZ))
      plt_biom%EPOOLR(1:npelms,N,L,NZ)=0._r8
      plt_biom%CCPOLR(N,L,NZ)=0._r8
      plt_biom%CZPOLR(N,L,NZ)=0._r8
      plt_biom%CPPOLR(N,L,NZ)=0._r8
      CWSRTL(N,L,NZ)=CWSRT(NZ)
      plt_biom%WTRTL(N,L,NZ)=0._r8
      plt_biom%WTRTD(N,L,NZ)=0._r8
      WSRTL(N,L,NZ)=0._r8
      plt_morph%RTN1(N,L,NZ)=0._r8
      plt_morph%RTNL(N,L,NZ)=0._r8
      plt_morph%RTLGP(N,L,NZ)=0._r8
      plt_morph%RTDNP(N,L,NZ)=0._r8
      RTVLP(N,L,NZ)=0._r8
      RTVLW(N,L,NZ)=0._r8
      RRAD1(N,L,NZ)=RRAD1M(N,NZ)
      RRAD2(N,L,NZ)=RRAD2M(N,NZ)
      plt_morph%RTARP(N,L,NZ)=0._r8
      plt_morph%RTLGA(N,L,NZ)=1.0E-03
      plt_rbgc%RUPNH4(N,L,NZ)=0._r8
      plt_rbgc%RUPNO3(N,L,NZ)=0._r8
      plt_rbgc%RUPH2P(N,L,NZ)=0._r8
      plt_rbgc%RUPH1P(N,L,NZ)=0._r8
      plt_rbgc%RUPNHB(N,L,NZ)=0._r8
      plt_rbgc%RUPNOB(N,L,NZ)=0._r8
      plt_rbgc%RUPH2B(N,L,NZ)=0._r8
      plt_rbgc%RUPH1B(N,L,NZ)=0._r8
      plt_rbgc%ROXYP(N,L,NZ)=0._r8
      plt_rbgc%RUNNHP(N,L,NZ)=0._r8
      plt_rbgc%RUNNBP(N,L,NZ)=0._r8
      plt_rbgc%RUNNOP(N,L,NZ)=0._r8
      plt_rbgc%RUNNXP(N,L,NZ)=0._r8
      plt_rbgc%RUPP2P(N,L,NZ)=0._r8
      plt_rbgc%RUPP1P(N,L,NZ)=0._r8
      plt_rbgc%RUPP2B(N,L,NZ)=0._r8
      plt_rbgc%RUPP1B(N,L,NZ)=0._r8
      CCO2A=CCO2EI
      CCO2P=0.030*EXP(-2.621_r8-0.0317_r8*ATCA)*CO2EI
      CO2A(N,L,NZ)=CCO2A*RTVLP(N,L,NZ)
      CO2P(N,L,NZ)=CCO2P*RTVLW(N,L,NZ)
      plt_rbgc%RCOFLA(N,L,NZ)=0._r8
      plt_rbgc%RCODFA(N,L,NZ)=0._r8
      plt_rbgc%RCO2S(N,L,NZ)=0._r8
      plt_rbgc%RCO2P(N,L,NZ)=0._r8
      COXYA=COXYE
      COXYP=0.032*EXP(-6.175_r8-0.0211_r8*ATCA)*OXYE
      plt_rbgc%OXYA(N,L,NZ)=COXYA*RTVLP(N,L,NZ)
      plt_rbgc%OXYP(N,L,NZ)=COXYP*RTVLW(N,L,NZ)
      plt_rbgc%CH4A(N,L,NZ)=0._r8
      plt_rbgc%CH4P(N,L,NZ)=0._r8
      plt_rbgc%Z2OA(N,L,NZ)=0._r8
      plt_rbgc%Z2OP(N,L,NZ)=0._r8
      plt_rbgc%ZH3A(N,L,NZ)=0._r8
      plt_rbgc%ZH3P(N,L,NZ)=0._r8
      plt_rbgc%H2GA(N,L,NZ)=0._r8
      plt_rbgc%H2GP(N,L,NZ)=0._r8
      plt_rbgc%WFR(N,L,NZ)=1.0
      D30: DO NR=1,JC1
        plt_morph%RTN2(N,L,NR,NZ)=0._r8
        plt_morph%RTLG1(N,L,NR,NZ)=0._r8
        plt_morph%RTLG2(N,L,NR,NZ)=0._r8
        plt_morph%RTDP1(N,NR,NZ)=SDPTH(NZ)
        plt_biom%WTRT1E(ielmc,N,L,NR,NZ)=0._r8
        plt_biom%WTRT1E(ielmn,N,L,NR,NZ)=0._r8
        plt_biom%WTRT1E(ielmp,N,L,NR,NZ)=0._r8
        plt_biom%WTRT2E(ielmc,N,L,NR,NZ)=0._r8
        plt_biom%WTRT2E(ielmn,N,L,NR,NZ)=0._r8
        plt_biom%WTRT2E(ielmp,N,L,NR,NZ)=0._r8
        plt_biom%RTWT1(N,NR,NZ)=0._r8
        plt_biom%RTWT1N(N,NR,NZ)=0._r8
        plt_biom%RTWT1P(N,NR,NZ)=0._r8
      ENDDO D30
      IF(N.EQ.1)THEN
        D6400: DO K=0,micpar%k_non_woody_litr
          plt_bgcr%ESNC(1:jsken,1:npelms,K,L,NZ)=0._r8
        ENDDO D6400
        plt_biom%EPOOLN(L,1:npelms,NZ)=0._r8
        plt_biom%WTNDLE(L,1:npelms,NZ)=0._r8
        RUPNF(L,NZ)=0._r8
      ENDIF
    ENDDO D20
  ENDDO D40

  plt_rbgc%RUPNH4(1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPNHB(1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPH2P(1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPH2B(1:2,NL+1:JZ1,NZ)=0._r8
  plt_morph%RTDNP(1:2,NL+1:JZ1,NZ)=0._r8
  end associate
  end subroutine InitRootMychorMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitSeedMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  REAL(R8) :: FDM

  associate(                             &
    PP       =>   plt_site%PP      , &
    PSILT    =>   plt_ew%PSILT     , &
    VOLWC    =>   plt_ew%VOLWC     , &
    VOLWP    =>   plt_ew%VOLWP     , &
    CWSRT    =>   plt_allom%CWSRT  , &
    CNGR     =>   plt_allom%CNGR   , &
    CPGR     =>   plt_allom%CPGR   , &
    RTWT1    =>   plt_biom%RTWT1   , &
    RTWT1N   =>   plt_biom%RTWT1N  , &
    RTWT1P   =>   plt_biom%RTWT1P  , &
    WTRVX    =>   plt_biom%WTRVX   , &
    WTRT1E   =>   plt_biom%WTRT1E  , &
    WTLSB    =>   plt_biom%WTLSB   , &
    WTLS     =>   plt_biom%WTLS    , &
    WTRTD    =>   plt_biom%WTRTD   , &
    WTSHEBE  =>   plt_biom%WTSHEBE , &
    WTRTL    =>   plt_biom%WTRTL   , &
    WSRTL    =>   plt_biom%WSRTL   , &
    EPOOLR   =>   plt_biom%EPOOLR  , &
    EPOOL    =>   plt_biom%EPOOL   , &
    WTLFBE   =>   plt_biom%WTLFBE  , &
    WTRVE    =>   plt_biom%WTRVE   , &
    GRDM     =>   plt_morph%GRDM   , &
    RTDP1    =>   plt_morph%RTDP1  , &
    NG       =>   plt_morph%NG       &
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
  WTRVX(NZ)=GRDM(NZ)*PP(NZ)
  WTRVE(ielmc,NZ)=WTRVX(NZ)
  WTRVE(ielmn,NZ)=CNGR(NZ)*WTRVE(ielmc,NZ)
  WTRVE(ielmp,NZ)=CPGR(NZ)*WTRVE(ielmc,NZ)
  WTLFBE(1,ielmn,NZ)=CNGR(NZ)*WTLFBE(1,ielmc,NZ)
  WTLFBE(1,ielmp,NZ)=CPGR(NZ)*WTLFBE(1,ielmc,NZ)
  WTLSB(1,NZ)=WTLFBE(1,ielmc,NZ)+WTSHEBE(1,ielmc,NZ)
  WTLS(NZ)=WTLS(NZ)+WTLSB(1,NZ)
  FDM=AMIN1(1.0_r8,0.16_r8-0.045_r8*PSILT(NZ))
  VOLWP(NZ)=ppmc*WTLS(NZ)/FDM
  VOLWC(NZ)=0._r8
  EPOOL(1,ielmn,NZ)=CNGR(NZ)*EPOOL(1,ielmc,NZ)
  EPOOL(1,ielmp,NZ)=CPGR(NZ)*EPOOL(1,ielmc,NZ)
  WTRT1E(ielmn,1,NG(NZ),1,NZ)=CNGR(NZ)*WTRT1E(ielmc,1,NG(NZ),1,NZ)
  WTRT1E(ielmp,1,NG(NZ),1,NZ)=CPGR(NZ)*WTRT1E(ielmc,1,NG(NZ),1,NZ)
  RTWT1N(1,1,NZ)=CNGR(NZ)*RTWT1(1,1,NZ)
  RTWT1P(1,1,NZ)=CPGR(NZ)*RTWT1(1,1,NZ)
  WTRTL(1,NG(NZ),NZ)=WTRT1E(ielmc,1,NG(NZ),1,NZ)
  WTRTD(1,NG(NZ),NZ)=WTRT1E(ielmc,1,NG(NZ),1,NZ)
  WSRTL(1,NG(NZ),NZ)=WTRTL(1,NG(NZ),NZ)*CWSRT(NZ)
  EPOOLR(ielmn,1,NG(NZ),NZ)=CNGR(NZ)*EPOOLR(ielmc,1,NG(NZ),NZ)
  EPOOLR(ielmp,1,NG(NZ),NZ)=CPGR(NZ)*EPOOLR(ielmc,1,NG(NZ),NZ)

  end associate
  end subroutine InitSeedMorphoBio

  end module StartqsMod
