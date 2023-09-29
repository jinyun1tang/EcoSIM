module StartqsMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use minimathmod, only : AZMAX1
  use EcoSIMConfig
  use PlantAPIData
  use TracerIDMod
  use EcoSiMParDataMod, only : pltpar
  use UnitMod, only : units
  use PlantMathFuncMod
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
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
    pftPlantPopulation      => plt_site%pftPlantPopulation        , &
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
!     CF,ClumpFactort0=current,initial clumping factor
!     MaxCanPStomaResistH2O=cuticular resistance to water (h m-1)
!     RCMX=cuticular resistance to CO2 (s m-1)
!     CNWS,CPWS=protein:N,protein:P ratios
!     RootFracRemobilizableBiom=maximum root protein concentration (g g-1)
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
        ZEROP(NZ)=ZERO*pftPlantPopulation(NZ)
        ZEROQ(NZ)=ZERO*pftPlantPopulation(NZ)/AREA3(NU)
        ZEROL(NZ)=ZERO*pftPlantPopulation(NZ)*1.0E+06_r8
      ENDDO D9985
!
!     FILL OUT UNUSED ARRAYS
!
      D9986: DO NZ=NP+1,JP1
        plt_bgcr%TESN0(1:NumOfPlantChemElements,NZ)=0._r8
        plt_bgcr%TESNC(1:NumOfPlantChemElements,NZ)=0._r8
        plt_biom%WTSTGE(1:NumOfPlantChemElements,NZ)=0._r8
        D6401: DO L=1,NL
          DO  K=1,pltpar%NumOfPlantLitrCmplxs
            plt_bgcr%ESNC(1:NumOfPlantChemElements,1:jsken,K,L,NZ)=0._r8
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
    RootFracRemobilizableBiom  =>  plt_allom%RootFracRemobilizableBiom , &
    CNWS   =>  plt_allom%CNWS  , &
    CPWS   =>  plt_allom%CPWS  , &
    CNRT   =>  plt_allom%CNRT  , &
    CPRT   =>  plt_allom%CPRT  , &
    O2I    =>  plt_photo%O2I   , &
    RCMX   =>  plt_photo%RCMX  , &
    ICTYP  =>  plt_photo%ICTYP , &
    MaxCanPStomaResistH2O   =>  plt_photo%MaxCanPStomaResistH2O  , &
    ClumpFactort0    =>  plt_morph%ClumpFactort0   , &
    ClumpFactor    =>  plt_morph%ClumpFactor   , &
    NRT    =>  plt_morph%NRT     &
  )
  IYR0(NZ)=IYRX(NZ)
  IDAY0(NZ)=IDAYX(NZ)
  IYRH(NZ)=IYRY(NZ)
  IDAYH(NZ)=IDAYY(NZ)
  PPI(NZ)=PPZ(NZ)
  PPX(NZ)=PPI(NZ)
  ClumpFactor(NZ)=ClumpFactort0(NZ)

  MaxCanPStomaResistH2O(NZ)=RSMX(NZ)/3600.0_r8
  RCMX(NZ)=RSMX(NZ)*1.56_r8
  CNWS(NZ)=2.5_r8
  CPWS(NZ)=25.0_r8
  RootFracRemobilizableBiom(NZ)=AMIN1(CNRT(NZ)*CNWS(NZ),CPRT(NZ)*CPWS(NZ))
  IF(ICTYP(NZ).EQ.ic3_photo)THEN
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

  associate(                         &
    instruct => pltpar%instruct, &
    ifoliar  => pltpar%ifoliar , &
    infoliar => pltpar%infoliar, &
    istalk   => pltpar%istalk  , &
    iroot    => pltpar%iroot   , &
    icwood   => pltpar%icwood  , &
    iprotein  => pltpar%iprotein  ,&
    icarbhyro => pltpar%icarbhyro ,&
    icellulos => pltpar%icellulos ,&
    ilignin  =>  pltpar%ilignin   , &
    Jlitgrp  => pltpar%Jlitgrp , &
    XRLA   =>  plt_pheno%XRLA      , &
    IBTYP  =>  plt_pheno%IBTYP     , &
    IGTYP  =>  plt_pheno%IGTYP     , &
    GROUPI =>  plt_pheno%GROUPI    , &
    CFOPE  =>  plt_soilchem%CFOPE  , &
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
  CFOPE(ielmc,instruct,iprotein,NZ)=0.0_r8
  CFOPE(ielmc,instruct,icarbhyro,NZ)=0.67_r8
  CFOPE(ielmc,instruct,icellulos,NZ)=0.33_r8
  CFOPE(ielmc,instruct,ilignin,NZ)=0.0_r8
!
!     NON-VASCULAR (E.G. MOSSES)
!
  IF(IGTYP(NZ).EQ.0)THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ)=0.30_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ)=0.38_r8

    CFOPE(ielmc,infoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,infoliar,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,infoliar,icellulos,NZ)=0.30_r8
    CFOPE(ielmc,infoliar,ilignin,NZ)=0.38_r8
!
!     LEGUMES
!
  ELSEIF(INTYP(NZ).NE.0)THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ)=0.16_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ)=0.38_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ)=0.34_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ)=0.12_r8

    CFOPE(ielmc,infoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,infoliar,icarbhyro,NZ)=0.41_r8
    CFOPE(ielmc,infoliar,icellulos,NZ)=0.37_r8
    CFOPE(ielmc,infoliar,ilignin,NZ)=0.15_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ)=0.08_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ)=0.41_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ)=0.36_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ)=0.15_r8

    CFOPE(ielmc,infoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,infoliar,icarbhyro,NZ)=0.41_r8
    CFOPE(ielmc,infoliar,icellulos,NZ)=0.36_r8
    CFOPE(ielmc,infoliar,ilignin,NZ)=0.16_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(IBTYP(NZ).EQ.1.OR.IBTYP(NZ).GE.3)THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ)=0.34_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ)=0.36_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ)=0.23_r8

    CFOPE(ielmc,infoliar,iprotein,NZ)=0.0_r8
    CFOPE(ielmc,infoliar,icarbhyro,NZ)=0.045_r8
    CFOPE(ielmc,infoliar,icellulos,NZ)=0.660_r8
    CFOPE(ielmc,infoliar,ilignin,NZ)=0.295_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPE(ielmc,ifoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ)=0.38_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ)=0.30_r8

    CFOPE(ielmc,infoliar,iprotein,NZ)=0.0_r8
    CFOPE(ielmc,infoliar,icarbhyro,NZ)=0.045_r8
    CFOPE(ielmc,infoliar,icellulos,NZ)=0.660_r8
    CFOPE(ielmc,infoliar,ilignin,NZ)=0.295_r8
  ENDIF
!
!     FRACTIONS OF WOODY LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     NON-VASCULAR
!
  IF(IGTYP(NZ).EQ.0)THEN
    CFOPE(ielmc,istalk,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,istalk,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,istalk,icellulos,NZ)=0.30_r8
    CFOPE(ielmc,istalk,ilignin,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
    CFOPE(ielmc,istalk,iprotein,NZ)=0.03_r8
    CFOPE(ielmc,istalk,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,istalk,icellulos,NZ)=0.57_r8
    CFOPE(ielmc,istalk,ilignin,NZ)=0.15_r8
!
!     DECIDUOUS AND CONIFEROUS TREES
!
  ELSE
    CFOPE(ielmc,istalk,iprotein,NZ)=0.0_r8
    CFOPE(ielmc,istalk,icarbhyro,NZ)=0.045_r8
    CFOPE(ielmc,istalk,icellulos,NZ)=0.660_r8
    CFOPE(ielmc,istalk,ilignin,NZ)=0.295_r8
  ENDIF
!
!     FRACTIONS OF FINE ROOT LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN PC&E 25:601-608
!
!     NON-VASCULAR
!
  IF(IGTYP(NZ).EQ.0)THEN
    CFOPE(ielmc,iroot,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,iroot,icellulos,NZ)=0.30_r8
    CFOPE(ielmc,iroot,ilignin,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
    CFOPE(ielmc,iroot,iprotein,NZ)=0.057_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ)=0.263_r8
    CFOPE(ielmc,iroot,icellulos,NZ)=0.542_r8
    CFOPE(ielmc,iroot,ilignin,NZ)=0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(IBTYP(NZ).EQ.1.OR.IBTYP(NZ).GE.3)THEN
    CFOPE(ielmc,iroot,iprotein,NZ)=0.059_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ)=0.308_r8
    CFOPE(ielmc,iroot,icellulos,NZ)=0.464_r8
    CFOPE(ielmc,iroot,ilignin,NZ)=0.169_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPE(ielmc,iroot,iprotein,NZ)=0.059_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ)=0.308_r8
    CFOPE(ielmc,iroot,icellulos,NZ)=0.464_r8
    CFOPE(ielmc,iroot,ilignin,NZ)=0.169_r8
  ENDIF
!
!     COARSE WOODY LITTER FROM BOLES AND ROOTS
!
  CFOPE(ielmc,icwood,iprotein,NZ)=0.00_r8
  CFOPE(ielmc,icwood,icarbhyro,NZ)=0.045_r8
  CFOPE(ielmc,icwood,icellulos,NZ)=0.660_r8
  CFOPE(ielmc,icwood,ilignin,NZ)=0.295_r8
!
!     INITIALIZE C-N AND C-P RATIOS IN PLANT LITTER
!
!     CNOPC,CPOPC=fractions to allocate N,P to kinetic components
!     CFOPN,CFOPP=distribution of litter N,P to kinetic components
!
  CNOPC(iprotein)=0.020_r8
  CNOPC(icarbhyro)=0.010_r8
  CNOPC(icellulos)=0.010_r8
  CNOPC(ilignin)=0.020_r8

  CPOPC(iprotein)=0.0020_r8
  CPOPC(icarbhyro)=0.0010_r8
  CPOPC(icellulos)=0.0010_r8
  CPOPC(ilignin)=0.0020_r8

  D110: DO N=0,Jlitgrp
    CNOPCT=0.0_r8
    CPOPCT=0.0_r8
    D100: DO M=1,jsken
      CNOPCT=CNOPCT+CFOPE(ielmc,N,M,NZ)*CNOPC(M)
      CPOPCT=CPOPCT+CFOPE(ielmc,N,M,NZ)*CPOPC(M)
    ENDDO D100
    D105: DO M=1,jsken
      CFOPE(ielmn,N,M,NZ)=CFOPE(ielmc,N,M,NZ)*CNOPC(M)/CNOPCT
      CFOPE(ielmp,N,M,NZ)=CFOPE(ielmc,N,M,NZ)*CPOPC(M)/CPOPCT
    ENDDO D105
  ENDDO D110
!
!     CONCURRENT NODE GROWTH
!
!     FNOD=scales node number for perennial vegetation (e.g. trees)
!     NNOD=number of concurrently growing nodes
!
  IF(IBTYP(NZ).EQ.0.OR.IGTYP(NZ).LE.1)THEN
    FNOD(NZ)=1.0_r8
    IF(GROUPI(NZ).LE.10)THEN
      NNOD(NZ)=3
    ELSEIF(GROUPI(NZ).LE.15)THEN
      NNOD(NZ)=4
    ELSE
      NNOD(NZ)=5
    ENDIF
  ELSE
    FNOD(NZ)=AMAX1(1.0_r8,0.04_r8/XRLA(NZ))
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
  IF(ICTYP(NZ).EQ.ic3_photo)THEN
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
    BiomGrowthYieldRoot     =>  plt_allom%BiomGrowthYieldRoot    , &
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
    CumSoilThickness   =>  plt_site%CumSoilThickness   , &
    NU       =>  plt_site%NU       , &
    NL       =>  plt_site%NL       , &
    NGTopRootLayer       =>  plt_morph%NGTopRootLayer      , &
    SeedinDepth    =>  plt_morph%SeedinDepth   , &
    SeedArea     =>  plt_morph%SeedArea    , &
    SeedCMass     =>  plt_morph%SeedCMass    , &
    PlantinDepth   =>  plt_morph%PlantinDepth  , &
    MaxSecndRootRadius1   =>  plt_morph%MaxSecndRootRadius1  , &
    SecndRootXSecArea   =>  plt_morph%SecndRootXSecArea  , &
    MaxPrimRootRadius1   =>  plt_morph%MaxPrimRootRadius1  , &
    PrimRootXSecArea   =>  plt_morph%PrimRootXSecArea  , &
    MaxSecndRootRadius   =>  plt_morph%MaxSecndRootRadius  , &
    MaxPrimRootRadius   =>  plt_morph%MaxPrimRootRadius  , &
    SecndRootSpecLen   =>  plt_morph%SecndRootSpecLen  , &
    NIXBotRootLayer     =>  plt_morph%NIXBotRootLayer    , &
    RSRR     =>  plt_morph%RSRR    , &
    PrimRootSpecLen   =>  plt_morph%PrimRootSpecLen  , &
    RootPorosity    =>  plt_morph%RootPorosity   , &
    RootPoreTortu4Gas    =>  plt_morph%RootPoreTortu4Gas   , &
    RRADP    =>  plt_morph%RRADP   , &
    DMVL     =>  plt_morph%DMVL    , &
    RSRA     =>  plt_morph%RSRA    , &
    NINR     =>  plt_morph%NINR    , &
    SeedVolume     =>  plt_morph%SeedVolume    , &
    SeedLength     =>  plt_morph%SeedLength      &
  )
!
!     SEED CHARACTERISTICS
!
!     SeedVolume,SeedLength,SeedArea=seed volume(m3),length(m),AREA3(NU)(m2)
!     SeedCMass=seed C mass (g) from PFT file
!
  call calc_seed_geometry(SeedCMass(NZ),SeedVolume(NZ),SeedLength(NZ),SeedArea(NZ))
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) DIMENSIONS, UPTAKE PARAMETERS
!
!     SeedinDepth=seeding depth(m) from PFT management file
!     CumSoilThickness=depth to soil layer bottom from surface(m)
!     NG,NIX,NINR=seeding,upper,lower rooting layer
!     CNRTS,CPRTS=N,P root growth yield
!     MaxPrimRootRadius,MaxSecndRootRadius=maximum primary,secondary mycorrhizal radius (m)
!     PORT=mycorrhizal porosity
!     UPMXZH,UPKMZH,UPMNZH=NH4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake(g m-2 h-1),Km(uM), min concn (uM)
!     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     RSRR,RSRA=radial,axial root resistivity (m2 MPa-1 h-1)
!
  SeedinDepth(NZ)=PlantinDepth(NZ)
  D9795: DO L=NU,NL
    IF(SeedinDepth(NZ).GE.CumSoilThickness(L-1).AND.SeedinDepth(NZ).LT.CumSoilThickness(L))THEN
      NGTopRootLayer(NZ)=L
      NIXBotRootLayer(NZ)=L
      D9790: DO NR=1,pltpar%JRS
        NINR(NR,NZ)=L
      ENDDO D9790
    ENDIF
  ENDDO D9795
  CNRTS(NZ)=CNRT(NZ)*BiomGrowthYieldRoot(NZ)
  CPRTS(NZ)=CPRT(NZ)*BiomGrowthYieldRoot(NZ)
  MaxPrimRootRadius(2,NZ)=5.0E-06_r8
  MaxSecndRootRadius(2,NZ)=5.0E-06_r8
  RootPorosity(2,NZ)=RootPorosity(1,NZ)
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
!     RootPoreTortu4Gas=tortuosity for gas transport
!     RRADP=path length for radial diffusion within root (m)
!     DMVL=volume:C ratio (m3 g-1)
!     PrimRootSpecLen,SecndRootSpecLen=specific primary,secondary root length (m g-1)
!     PrimRootXSecArea,SecndRootXSecArea=specific primary,secondary root area (m2 g-1)
!
  D500: DO N=1,2
    RootPoreTortu4Gas(N,NZ)=RootPorosity(N,NZ)**1.33_r8
    RRADP(N,NZ)=LOG(1.0_r8/SQRT(AMAX1(0.01_r8,RootPorosity(N,NZ))))
    DMVL(N,NZ)=ppmc/(0.05_r8*(1.0_r8-RootPorosity(N,NZ)))
    PrimRootSpecLen(N,NZ)=DMVL(N,NZ)/(PICON*MaxPrimRootRadius(N,NZ)**2._r8)
    SecndRootSpecLen(N,NZ)=DMVL(N,NZ)/(PICON*MaxSecndRootRadius(N,NZ)**2._r8)
    MaxPrimRootRadius1(N,NZ)=MaxPrimRootRadius(N,NZ)
!    2*SQRT(0.25*(1.0-RootPorosity(N,NZ)))
    MaxSecndRootRadius1(N,NZ)=MaxSecndRootRadius(N,NZ)
!    2*SQRT(0.25*(1.0-RootPorosity(N,NZ)))
    PrimRootXSecArea(N,NZ)=PICON*MaxPrimRootRadius1(N,NZ)**2
    SecndRootXSecArea(N,NZ)=PICON*MaxSecndRootRadius1(N,NZ)**2
  ENDDO D500
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
    pftPlantPopulation      =>  plt_site%pftPlantPopulation      , &
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
    HourCounter4LeafOut_brch    =>  plt_pheno%HourCounter4LeafOut_brch   , &
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
    CanopyHeight      =>  plt_morph%CanopyHeight     , &
    KLEAF   =>  plt_morph%KLEAF  , &
    XTLI    =>  plt_morph%XTLI   , &
    NBT     =>  plt_morph%NBT    , &
    PSTG    =>  plt_morph%PSTG   , &
    CanopyLeafA_pft   =>  plt_morph%CanopyLeafA_pft  , &
    CanopyBranchStemApft_lyr   =>  plt_morph%CanopyBranchStemApft_lyr  , &
    SURFB   =>  plt_morph%SURFB  , &
    CanPSA   =>  plt_morph%CanPSA  , &
    PSTGF   =>  plt_morph%PSTGF  , &
    GRNXB   =>  plt_morph%GRNXB  , &
    HTNODE  =>  plt_morph%HTNODE , &
    GRNOB   =>  plt_morph%GRNOB  , &
    CanopyBranchLeafA_pft   =>  plt_morph%CanopyBranchLeafA_pft  , &
    ARLFZ   =>  plt_morph%ARLFZ  , &
    CanPLNBLA   =>  plt_morph%CanPLNBLA  , &
    CanopyLeafApft_lyr   =>  plt_morph%CanopyLeafApft_lyr  , &
    CanopyStemApft_lyr   =>  plt_morph%CanopyStemApft_lyr  , &
    SURF    =>  plt_morph%SURF   , &
    ARLF1   =>  plt_morph%ARLF1  , &
    CanPBranchHeight  =>  plt_morph%CanPBranchHeight , &
    HypoctoylHeight   =>  plt_morph%HypoctoylHeight  , &
    BranchNumber_brchpft    =>  plt_morph%BranchNumber_brchpft   , &
    PSTGI   =>  plt_morph%PSTGI  , &
    KLEAFX  =>  plt_morph%KLEAFX , &
    NumOfBranches_pft     =>  plt_morph%NumOfBranches_pft      &
  )
!
!     INITIALIZE PLANT PHENOLOGY
!
!     PP=population (grid cell-1)
!
  pftPlantPopulation(NZ)=PPX(NZ)*AREA3(NU)
  plt_pheno%IFLGI(NZ)=0
  plt_pheno%IDTHP(NZ)=0
  plt_pheno%IDTHR(NZ)=0
  NBT(NZ)=0
  NumOfBranches_pft(NZ)=0
  HypoctoylHeight(NZ)=0._r8
  CanopyHeight(NZ)=0._r8
  D10: DO NB=1,JBR
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
    HourCounter4LeafOut_brch(NB,NZ)=0._r8
    FDBK(NB,NZ)=1.0
    FDBKX(NB,NZ)=1.0
    FLG4(NB,NZ)=0
    FLGZ(NB,NZ)=0
    BranchNumber_brchpft(NB,NZ)=0
    plt_pheno%IDTHB(NB,NZ)=ibrdead
    D15: DO M=1,pltpar%jpstgs
      IDAY(M,NB,NZ)=0
    ENDDO D15
  ENDDO D10
!
!     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
!
  WSTR(NZ)=0._r8
  CHILL(NZ)=0._r8
  plt_biom%EPOOL(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%EPOLNB(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WTSHEBE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WTSHTBE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WTSTKBE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WTLFBE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WTRSVBE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WTHSKBE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WTGRBE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WTEARBE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WTNDBE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_pheno%RCELX(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WGLFEX(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_pheno%RCESX(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WGSHEXE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8
  plt_biom%WTSTXBE(1:NumOfPlantChemElements,1:JBR,NZ)=0._r8

  D25: DO NB=1,JBR
    plt_biom%CanPBStalkC(NB,NZ)=0._r8
    plt_biom%CanPBLeafShethC(NB,NZ)=0._r8
    GRNXB(NB,NZ)=0._r8
    GRNOB(NB,NZ)=0._r8
    plt_allom%GRWTB(NB,NZ)=0._r8
    CanopyBranchLeafA_pft(NB,NZ)=0._r8
    plt_rbgc%RNH3B(NB,NZ)=0._r8
    ARLFZ(NB,NZ)=0._r8
    CanPBranchHeight(NB,NZ)=0._r8
    D5: DO L=1,NumOfCanopyLayers1
      CanopyBranchStemApft_lyr(L,NB,NZ)=0._r8
      DO N=1,JLI1
        SURFB(N,L,NB,NZ)=0._r8
      enddo
    ENDDO D5
    DO K=0,JNODS1
      ARLF1(K,NB,NZ)=0._r8
      HTNODE(K,NB,NZ)=0._r8
      plt_morph%HTNODX(K,NB,NZ)=0._r8
      plt_morph%CanPSheathHeight(K,NB,NZ)=0._r8
      plt_biom%WSLF(K,NB,NZ)=0._r8
      plt_biom%WSSHE(K,NB,NZ)=0._r8
      plt_biom%WGLFE(1:NumOfPlantChemElements,K,NB,NZ)=0._r8
      plt_biom%WGNODE(1:NumOfPlantChemElements,K,NB,NZ)=0._r8
      plt_biom%WGSHE(1:NumOfPlantChemElements,K,NB,NZ)=0._r8


      D55: DO L=1,NumOfCanopyLayers1
        CanPLNBLA(L,K,NB,NZ)=0._r8
        plt_biom%WGLFLE(1:NumOfPlantChemElements,L,K,NB,NZ)=0._r8
      ENDDO D55
      IF(K.NE.0)THEN
        CPOOL3(K,NB,NZ)=0._r8
        CO2B(K,NB,NZ)=0._r8
        HCOB(K,NB,NZ)=0._r8
        CPOOL4(K,NB,NZ)=0._r8
        D45: DO L=1,NumOfCanopyLayers1
          DO N=1,JLI1
            SURF(N,L,K,NB,NZ)=0._r8
          enddo
        ENDDO D45
      ENDIF
    enddo
  ENDDO D25
  D35: DO L=1,NumOfCanopyLayers1
    CanopyLeafApft_lyr(L,NZ)=0._r8
    plt_biom%CanopyLeafCpft_lyr(L,NZ)=0._r8
    CanopyStemApft_lyr(L,NZ)=0._r8
  ENDDO D35
  plt_biom%CanopyNonstructElements_pft(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%CanopyNonstructElementConc_pft(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%NoduleNonstructCconc_pft(NZ)=0._r8
  plt_biom%CanPShootElmMass(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%WTLFE(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%WTSHEE(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%WTSTKE(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%CanPStalkC(NZ)=0._r8
  plt_biom%WTRSVE(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%WTHSKE(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%WTEARE(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%WTGRE(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%WTRTSE(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%WTNDE(1:NumOfPlantChemElements,NZ)=0._r8
  plt_biom%CanopyLeafShethC_pft(NZ)=0._r8
  CanopyLeafA_pft(NZ)=0._r8
  plt_biom%WTRTA(NZ)=0._r8
  CanPSA(NZ)=0._r8
  end associate
  end subroutine InitPlantPhenoMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitMassBalance(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: M,NE
  real(r8) :: WTSTDX

  associate(                           &
    NU     => plt_site%NU        , &
    AREA3  => plt_site%AREA3     , &
    CFOPE  => plt_soilchem%CFOPE , &
    ETCanP  => plt_ew%ETCanP       , &
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
    CanPSA  => plt_morph%CanPSA    , &
    icwood => pltpar%icwood      , &
    RSETE  => plt_pheno%RSETE      &

  )
!
!     INITIALIZE MASS BALANCE CHECKS
!
  IF(.not.is_restart().AND.is_first_year)THEN
    CARBN(NZ)=0._r8
    TESN0(1:NumOfPlantChemElements,NZ)=0._r8
    TCO2T(NZ)=0._r8
    TCO2A(NZ)=0._r8
    TEUPTK(1:NumOfPlantChemElements,NZ)=0._r8
    TESNC(1:NumOfPlantChemElements,NZ)=0._r8

    TZUPFX(NZ)=0._r8
    RNH3C(NZ)=0._r8
    TNH3C(NZ)=0._r8
    plt_distb%VCO2F(NZ)=0._r8
    plt_distb%VCH4F(NZ)=0._r8
    plt_distb%VOXYF(NZ)=0._r8
    plt_distb%VNH3F(NZ)=0._r8
    plt_distb%VN2OF(NZ)=0._r8
    plt_distb%VPO4F(NZ)=0._r8
    plt_distb%THVSTE(1:NumOfPlantChemElements,NZ)=0._r8
    plt_distb%HVSTE(1:NumOfPlantChemElements,NZ)=0._r8
    RSETE(1:NumOfPlantChemElements,NZ)=0._r8
    ETCanP(NZ)=0._r8
    WTSTGE(1:NumOfPlantChemElements,NZ)=0._r8
    WTSTDX=WTSTDI(NZ)*AREA3(NU)
    D155: DO M=1,jsken
      WTSTDE(ielmc,M,NZ)=WTSTDX*CFOPE(ielmc,icwood,M,NZ)
      WTSTDE(ielmn,M,NZ)=WTSTDX*CNSTK(NZ)*CFOPE(ielmn,icwood,M,NZ)
      WTSTDE(ielmp,M,NZ)=WTSTDX*CPSTK(NZ)*CFOPE(ielmp,icwood,M,NZ)
    ENDDO D155
    DO NE=1,NumOfPlantChemElements
      WTSTGE(NE,NZ)=WTSTGE(NE,NZ)+sum(WTSTDE(NE,1:jsken,NZ))
    ENDDO
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
    PTrans      =>  plt_ew%PTrans     , &
    VHeatCapCanP    =>  plt_ew%VHeatCapCanP   , &
    PSICanP    =>  plt_ew%PSICanP   , &
    PSICanPTurg    =>  plt_ew%PSICanPTurg   , &
    PSICanPOsmo    =>  plt_ew%PSICanPOsmo   , &
    ENGYX    =>  plt_ew%ENGYX   , &
    DTKC     =>  plt_ew%DTKC    , &
    TCC      =>  plt_ew%TCC     , &
    TKG      =>  plt_pheno%TKG  , &
    TCG      =>  plt_pheno%TCG  , &
    fTgrowCanP     =>  plt_pheno%fTgrowCanP , &
    CanPShootElmMass   =>  plt_biom%CanPShootElmMass , &
    FracPARByCanP    =>  plt_rad%FracPARByCanP    &
  )
!
!     INITIALIZE PLANT HEAT AND WATER STATUS
!
!     VHeatCapCanP=canopy heat capacity (MJ m-3 K-1)
!     TCC,TKC=canopy temperature for growth (oC,K)
!     TCG,TKG=canopy temperature for phenology (oC,K)
!     PSICanP,PSICanPOsmo,PSICanPTurg=canopy total,osmotic,turgor water potl(MPa)
!
  VHeatCapCanP(NZ)=cpw*CanPShootElmMass(ielmc,NZ)*10.0E-06
  ENGYX(NZ)=0._r8
  DTKC(NZ)=0._r8
  TCC(NZ)=ATCA
  TKC(NZ)=units%Celcius2Kelvin(TCC(NZ))
  TCG(NZ)=TCC(NZ)
  TKG(NZ)=units%Celcius2Kelvin(TCG(NZ))
  fTgrowCanP(NZ)=1.0
  PSICanP(NZ)=-1.0E-03
  PSICanPOsmo(NZ)=OSMO(NZ)+PSICanP(NZ)
  PSICanPTurg(NZ)=AZMAX1(PSICanP(NZ)-PSICanPOsmo(NZ))
  PTrans(NZ)=0._r8
  FracPARByCanP(NZ)=0._r8
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
    PSICanP    =>  plt_ew%PSICanP      , &
    PSIRoot    =>  plt_ew%PSIRoot      , &
    PSIRootTurg    =>  plt_ew%PSIRootTurg      , &
    PSIRootOSMO    =>  plt_ew%PSIRootOSMO      , &
    RUPNF    =>  plt_bgcr%RUPNF    , &
    trcg_rootml     =>  plt_rbgc%trcg_rootml     , &
    trcs_rootml     =>  plt_rbgc%trcs_rootml     , &
    OXYE     =>  plt_site%OXYE     , &
    COXYE    =>  plt_site%COXYE    , &
    CO2EI    =>  plt_site%CO2EI    , &
    NL       =>  plt_site%NL       , &
    ATCA     =>  plt_site%ATCA     , &
    CCO2EI   =>  plt_site%CCO2EI   , &
    RootFracRemobilizableBiom    =>  plt_allom%RootFracRemobilizableBiom   , &
    RootProteinConc_pftvr   =>  plt_biom%RootProteinConc_pftvr   , &
    WSRTL    =>  plt_biom%WSRTL    , &
    SeedinDepth    =>  plt_morph%SeedinDepth   , &
    RTVLW    =>  plt_morph%RTVLW   , &
    RTVLP    =>  plt_morph%RTVLP   , &
    SecndRootRadius    =>  plt_morph%SecndRootRadius   , &
    PrimRootRadius    =>  plt_morph%PrimRootRadius   , &
    MaxPrimRootRadius   =>  plt_morph%MaxPrimRootRadius  , &
    MaxSecndRootRadius   =>  plt_morph%MaxSecndRootRadius  , &
    NRT      =>  plt_morph%NRT       &
  )
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
!
!     PSIRoot,PSIRootOSMO,PSIRootTurg=root,myco total,osmotic,turgor water potl(MPa)
!     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
!     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
!
  NRT(NZ)=0
  plt_rbgc%UPNH4(NZ)=0._r8
  plt_rbgc%UPNO3(NZ)=0._r8
  plt_rbgc%UPH2P(NZ)=0._r8
  plt_rbgc%UPH1P(NZ)=0._r8
  plt_rbgc%UPNF(NZ)=0._r8
  D40: DO N=1,pltpar%jroots
    D20: DO L=1,NL
      plt_ew%AllPlantRootH2OUptake_vr(N,L,NZ)=0._r8
      PSIRoot(N,L,NZ)=-0.01
      PSIRootOSMO(N,L,NZ)=OSMO(NZ)+PSIRoot(N,L,NZ)
      PSIRootTurg(N,L,NZ)=AZMAX1(PSIRoot(N,L,NZ)-PSIRootOSMO(N,L,NZ))
      plt_biom%EPOOLR(1:NumOfPlantChemElements,N,L,NZ)=0._r8
      plt_biom%RootNonstructElementConcpft_vr(1:NumOfPlantChemElements,N,L,NZ)=0._r8
      RootProteinConc_pftvr(N,L,NZ)=RootFracRemobilizableBiom(NZ)
      plt_biom%WTRTL(N,L,NZ)=0._r8
      plt_biom%PopPlantRootC_vr(N,L,NZ)=0._r8
      WSRTL(N,L,NZ)=0._r8
      plt_morph%PrimRootXNumL(N,L,NZ)=0._r8
      plt_morph%SecndRootXNumL(N,L,NZ)=0._r8
      plt_morph%RootLenPerP(N,L,NZ)=0._r8
      plt_morph%RootLenDensNLP(N,L,NZ)=0._r8
      RTVLP(N,L,NZ)=0._r8
      RTVLW(N,L,NZ)=0._r8
      PrimRootRadius(N,L,NZ)=MaxPrimRootRadius(N,NZ)
      SecndRootRadius(N,L,NZ)=MaxSecndRootRadius(N,NZ)
      plt_morph%RTARP(N,L,NZ)=0._r8
      plt_morph%AveSecndRootLen(N,L,NZ)=1.0E-03
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
      trcg_rootml(idg_CO2,N,L,NZ)=CCO2A*RTVLP(N,L,NZ)
      trcs_rootml(idg_CO2,N,L,NZ)=CCO2P*RTVLW(N,L,NZ)
      plt_rbgc%trcg_RFLA(idg_CO2,N,L,NZ)=0._r8
      plt_rbgc%trcg_RDFA(idg_CO2,N,L,NZ)=0._r8
      plt_rbgc%RCO2S(N,L,NZ)=0._r8
      plt_rbgc%RCO2P(N,L,NZ)=0._r8
      COXYA=COXYE
      COXYP=0.032_r8*EXP(-6.175_r8-0.0211_r8*ATCA)*OXYE
      plt_rbgc%trcg_rootml(idg_O2,N,L,NZ)=COXYA*RTVLP(N,L,NZ)
      plt_rbgc%trcs_rootml(idg_O2,N,L,NZ)=COXYP*RTVLW(N,L,NZ)
      plt_rbgc%trcg_rootml(idg_beg:idg_end-1,N,L,NZ)=0._r8
      plt_rbgc%trcs_rootml(idg_beg:idg_end-1,N,L,NZ)=0._r8
      plt_rbgc%WFR(N,L,NZ)=1.0
      D30: DO NR=1,JRS
        plt_morph%RTN2(N,L,NR,NZ)=0._r8
        plt_morph%PrimRootLen(N,L,NR,NZ)=0._r8
        plt_morph%SecndRootLen(N,L,NR,NZ)=0._r8
        plt_morph%PrimRootDepth(N,NR,NZ)=SeedinDepth(NZ)
        plt_biom%WTRT1E(1:NumOfPlantChemElements,N,L,NR,NZ)=0._r8
        plt_biom%WTRT2E(1:NumOfPlantChemElements,N,L,NR,NZ)=0._r8
        plt_biom%RTWT1E(1:NumOfPlantChemElements,N,NR,NZ)=0._r8
      ENDDO D30
      IF(N.EQ.1)THEN
        D6400: DO K=1,pltpar%NumOfPlantLitrCmplxs
          plt_bgcr%ESNC(1:NumOfPlantChemElements,1:jsken,K,L,NZ)=0._r8
        ENDDO D6400
        plt_biom%EPOOLN(1:NumOfPlantChemElements,L,NZ)=0._r8
        plt_biom%WTNDLE(1:NumOfPlantChemElements,L,NZ)=0._r8
        RUPNF(L,NZ)=0._r8
      ENDIF
    ENDDO D20
  ENDDO D40

  plt_rbgc%RUPNH4(1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPNHB(1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPH2P(1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPH2B(1:2,NL+1:JZ1,NZ)=0._r8
  plt_morph%RootLenDensNLP(1:2,NL+1:JZ1,NZ)=0._r8
  end associate
  end subroutine InitRootMychorMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitSeedMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  REAL(R8) :: FDM

  associate(                             &
    pftPlantPopulation       =>   plt_site%pftPlantPopulation      , &
    PSICanP    =>   plt_ew%PSICanP     , &
    WatByPCan    =>   plt_ew%WatByPCan     , &
    CanWatP    =>   plt_ew%CanWatP     , &
    RootFracRemobilizableBiom    =>   plt_allom%RootFracRemobilizableBiom  , &
    CNGR     =>   plt_allom%CNGR   , &
    CPGR     =>   plt_allom%CPGR   , &
    RTWT1E   =>   plt_biom%RTWT1E  , &
    WTRVX    =>   plt_biom%WTRVX   , &
    WTRT1E   =>   plt_biom%WTRT1E  , &
    CanPBLeafShethC    =>   plt_biom%CanPBLeafShethC   , &
    CanopyLeafShethC_pft     =>   plt_biom%CanopyLeafShethC_pft    , &
    PopPlantRootC_vr    =>   plt_biom%PopPlantRootC_vr   , &
    WTSHEBE  =>   plt_biom%WTSHEBE , &
    WTRTL    =>   plt_biom%WTRTL   , &
    WSRTL    =>   plt_biom%WSRTL   , &
    EPOOLR   =>   plt_biom%EPOOLR  , &
    EPOOL    =>   plt_biom%EPOOL   , &
    WTLFBE   =>   plt_biom%WTLFBE  , &
    WTRVE    =>   plt_biom%WTRVE   , &
    SeedCMass     =>   plt_morph%SeedCMass   , &
    PrimRootDepth    =>   plt_morph%PrimRootDepth  , &
    NGTopRootLayer      =>   plt_morph%NGTopRootLayer      &
  )
!
!     INITIALIZE SEED MORPHOLOGY AND BIOMASS
!
!     WTRVC,WTRVN,WTRVP=C,N,P in storage reserves (g)
!     WTLFB,WTLFBN,WTLFBP=C,N,P in leaves (g)
!     CanPBLeafShethC=C in leaves+petioles (g)
!     FDM-dry matter fraction (g DM C g FM C-1)
!     CanWatP,WatByPCan=water volume in,on canopy (m3)
!     CPOOL,ZPOOL,PPOOL=C,N,P in canopy nonstructural pools (g)
!     WTRT1,WTRT1N,WTRT1P=C,N,P in primary root layer (g)
!     RTWT1,RTWT1N,RTWT1P=total C,N,P in primary root (g)
!     WTRTL,PopPlantRootC_vr=total root C mass (g)
!     WSRTL=total root protein C mass (g)
!     CPOOLR,ZPOOLR,PPOOLR=C,N,P in root,myco nonstructural pools (g)
!
  WTRVX(NZ)=SeedCMass(NZ)*pftPlantPopulation(NZ)
  WTRVE(ielmc,NZ)=WTRVX(NZ)
  WTRVE(ielmn,NZ)=CNGR(NZ)*WTRVE(ielmc,NZ)
  WTRVE(ielmp,NZ)=CPGR(NZ)*WTRVE(ielmc,NZ)
  WTLFBE(ielmn,1,NZ)=CNGR(NZ)*WTLFBE(ielmc,1,NZ)
  WTLFBE(ielmp,1,NZ)=CPGR(NZ)*WTLFBE(ielmc,1,NZ)
  CanPBLeafShethC(1,NZ)=WTLFBE(ielmc,1,NZ)+WTSHEBE(ielmc,1,NZ)
  CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+CanPBLeafShethC(1,NZ)
  FDM=AMIN1(1.0_r8,0.16_r8-0.045_r8*PSICanP(NZ))
  CanWatP(NZ)=ppmc*CanopyLeafShethC_pft(NZ)/FDM
  WatByPCan(NZ)=0._r8
  EPOOL(ielmn,1,NZ)=CNGR(NZ)*EPOOL(ielmc,1,NZ)
  EPOOL(ielmp,1,NZ)=CPGR(NZ)*EPOOL(ielmc,1,NZ)
  WTRT1E(ielmn,ipltroot,NGTopRootLayer(NZ),1,NZ)=CNGR(NZ)*WTRT1E(ielmc,ipltroot,NGTopRootLayer(NZ),1,NZ)
  WTRT1E(ielmp,ipltroot,NGTopRootLayer(NZ),1,NZ)=CPGR(NZ)*WTRT1E(ielmc,ipltroot,NGTopRootLayer(NZ),1,NZ)
  RTWT1E(ielmn,1,1,NZ)=CNGR(NZ)*RTWT1E(ielmc,1,1,NZ)
  RTWT1E(ielmp,1,1,NZ)=CPGR(NZ)*RTWT1E(ielmc,1,1,NZ)
  WTRTL(ipltroot,NGTopRootLayer(NZ),NZ)=WTRT1E(ielmc,ipltroot,NGTopRootLayer(NZ),1,NZ)
  PopPlantRootC_vr(ipltroot,NGTopRootLayer(NZ),NZ)=WTRT1E(ielmc,ipltroot,NGTopRootLayer(NZ),1,NZ)
  WSRTL(1,NGTopRootLayer(NZ),NZ)=WTRTL(ipltroot,NGTopRootLayer(NZ),NZ)*RootFracRemobilizableBiom(NZ)
  EPOOLR(ielmn,1,NGTopRootLayer(NZ),NZ)=CNGR(NZ)*EPOOLR(ielmc,1,NGTopRootLayer(NZ),NZ)
  EPOOLR(ielmp,1,NGTopRootLayer(NZ),NZ)=CPGR(NZ)*EPOOLR(ielmc,1,NGTopRootLayer(NZ),NZ)

  end associate
  end subroutine InitSeedMorphoBio

  end module StartqsMod
