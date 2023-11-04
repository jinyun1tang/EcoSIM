module StartqMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GridConsts
  use FlagDataType
  use EcosimConst
  use TracerIDMod
  use EcoSIMCtrlDataType
  use minimathmod, only : AZMAX1
  use PlantDataRateType
  use ClimForcDataType
  use PlantTraitDataType
  use PlantMngmtDataType
  use CanopyDataType
  use CanopyRadDataType
  use RootDataType
  use EcoSIMHistMod
  use GridDataType
  use EcoSIMConfig
  use GrosubPars
  use UnitMod, only : units
  use EcoSiMParDataMod, only : pltpar
  use PlantMathFuncMod
  implicit none

  private

  public :: startq
  contains

  SUBROUTINE startq(NHWQ,NHEQ,NVNQ,NVSQ,NZ1Q,NZ2Q)
!
!     THIS SUBROUTINE INITIALIZES ALL PLANT VARIABLES
!
  implicit none
  integer, intent(in) :: NHWQ,NHEQ,NVNQ,NVSQ
  integer, intent(in) :: NZ1Q,NZ2Q

  integer :: NY,NX,K,L,M,NZ,NZ2X
!     begin_execution
!
!     INITIALIZE SHOOT GROWTH VARIABLES
!
!     IsPlantActive=PFT flag:0=not active,1=active
!     iYearPlanting,iDayPlanting,iYearPlantHarvest,iDayPlantHarvest=year,day of planting,arvesting
!     PPI,PPX=initial,current population (m-2)
!     CF,ClumpFactort0=current,initial clumping factor
!     MaxCanPStomaResistH2O=cuticular resistance to water (h m-1)
!     RCMX=cuticular resistance to CO2 (s m-1)
!     CNWS,CPWS=protein:N,protein:P ratios
!    RootFracRemobilizableBiom=maximum root protein concentration (g g-1)
!     O2I=intercellular O2 concentration in C3,C4 PFT (umol mol-1)
!

  D9995: DO NX=NHWQ,NHEQ
    D9990: DO NY=NVNQ,NVSQ
      NZ2X=MIN(NZ2Q,NP(NY,NX))
      D9985: DO NZ=NZ1Q,NZ2X
        IF(IsPlantActive(NZ,NY,NX).EQ.iPlantIsDormant)THEN

          call InitShootGrowth(NZ,NY,NX)

          call PlantLitterFractions(NZ,NY,NX)

          call PFTThermalAcclimation(NZ,NY,NX)

          call InitDimensionsandUptake(NZ,NY,NX)

          call InitPlantPhenoMorphoBio(NZ,NY,NX)

          call InitMassBalance(NZ,NY,NX)

          call InitPlantHeatandWater(NZ,NY,NX)

          call InitRootMychorMorphoBio(NZ,NY,NX)

          call InitSeedMorphoBio(NZ,NY,NX)
        ENDIF
        ZEROP(NZ,NY,NX)=ZERO*pftPlantPopulation(NZ,NY,NX)
        ZEROQ(NZ,NY,NX)=ZERO*pftPlantPopulation(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        ZEROL(NZ,NY,NX)=ZERO*pftPlantPopulation(NZ,NY,NX)*1.0E+06
      ENDDO D9985
!
!     FILL OUT UNUSED ARRAYS
!
      D9986: DO NZ=NP(NY,NX)+1,JP
        TESN0(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
        TESNC(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
        WTSTGE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
        D6401: DO L=1,NL(NY,NX)
          DO  K=1,pltpar%NumOfPlantLitrCmplxs
            DO  M=1,jskenc
              ESNC(1:NumOfPlantChemElements,M,K,L,NZ,NY,NX)=0._r8
            enddo
          enddo
        ENDDO D6401
      ENDDO D9986
    ENDDO D9990
  ENDDO D9995
  RETURN
  END subroutine startq
!------------------------------------------------------------------------------------------

  subroutine InitShootGrowth(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX

  iYearPlanting(NZ,NY,NX)=IYRX(NZ,NY,NX)   !planting year
  iDayPlanting(NZ,NY,NX)=IDAYX(NZ,NY,NX) !planting day
  iYearPlantHarvest(NZ,NY,NX)=IYRY(NZ,NY,NX)
  iDayPlantHarvest(NZ,NY,NX)=IDAYY(NZ,NY,NX)
  PPI(NZ,NY,NX)=PPZ(NZ,NY,NX)
  PPX(NZ,NY,NX)=PPI(NZ,NY,NX)
  ClumpFactor(NZ,NY,NX)=ClumpFactort0(NZ,NY,NX)       !clumping factor
  
  MaxCanPStomaResistH2O(NZ,NY,NX)=RSMX(NZ,NY,NX)/3600.0_r8
  RCMX(NZ,NY,NX)=RSMX(NZ,NY,NX)*1.56_r8
  CNWS(NZ,NY,NX)=2.5_r8
  CPWS(NZ,NY,NX)=25.0_r8
 RootFracRemobilizableBiom(NZ,NY,NX)=AMIN1(CNRT(NZ,NY,NX)*CNWS(NZ,NY,NX),CPRT(NZ,NY,NX)*CPWS(NZ,NY,NX))
  IF(ICTYP(NZ,NY,NX).EQ.ic3_photo)THEN
    O2I(NZ,NY,NX)=2.10E+05_r8
  ELSE
    O2I(NZ,NY,NX)=3.96E+05_r8
  ENDIF
  end subroutine InitShootGrowth
!------------------------------------------------------------------------------------------

  subroutine PlantLitterFractions(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer :: N,M
  real(r8) :: CNOPC(4),CPOPC(4)
  REAL(R8) :: CNOPCT,CPOPCT

  associate(                       &
    iprotein  => pltpar%iprotein  ,&
    icarbhyro => pltpar%icarbhyro ,&
    icellulos => pltpar%icellulos ,&
    ilignin  =>  pltpar%ilignin   , &
    Jlitgrp  => pltpar%Jlitgrp , &
    instruct => pltpar%instruct, &
    ifoliar  => pltpar%ifoliar , &
    infoliar => pltpar%infoliar, &
    istalk   => pltpar%istalk  , &
    iroot    => pltpar%iroot   , &
    icwood   => pltpar%icwood    &
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
  CFOPE(ielmc,instruct,iprotein,NZ,NY,NX)=0.0_r8
  CFOPE(ielmc,instruct,icarbhyro,NZ,NY,NX)=0.67_r8
  CFOPE(ielmc,instruct,icellulos,NZ,NY,NX)=0.33_r8
  CFOPE(ielmc,instruct,ilignin,NZ,NY,NX)=0.0_r8
!
!     NON-VASCULAR (E.G. MOSSES)
!
  IF(IGTYP(NZ,NY,NX).EQ.0)THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ,NY,NX)=0.07_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ,NY,NX)=0.25_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ,NY,NX)=0.30_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ,NY,NX)=0.38_r8

    CFOPE(ielmc,infoliar,iprotein,NZ,NY,NX)=0.07_r8
    CFOPE(ielmc,infoliar,icarbhyro,NZ,NY,NX)=0.25_r8
    CFOPE(ielmc,infoliar,icellulos,NZ,NY,NX)=0.30_r8
    CFOPE(ielmc,infoliar,ilignin,NZ,NY,NX)=0.38_r8
!
!     LEGUMES
!
  ELSEIF(INTYP(NZ,NY,NX).NE.0)THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ,NY,NX)=0.16_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ,NY,NX)=0.38_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ,NY,NX)=0.34_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ,NY,NX)=0.12_r8

    CFOPE(ielmc,infoliar,iprotein,NZ,NY,NX)=0.07_r8
    CFOPE(ielmc,infoliar,icarbhyro,NZ,NY,NX)=0.41_r8
    CFOPE(ielmc,infoliar,icellulos,NZ,NY,NX)=0.37_r8
    CFOPE(ielmc,infoliar,ilignin,NZ,NY,NX)=0.15_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ,NY,NX)=0.08_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ,NY,NX)=0.41_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ,NY,NX)=0.36_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ,NY,NX)=0.15_r8

    CFOPE(ielmc,infoliar,iprotein,NZ,NY,NX)=0.07_r8
    CFOPE(ielmc,infoliar,icarbhyro,NZ,NY,NX)=0.41_r8
    CFOPE(ielmc,infoliar,icellulos,NZ,NY,NX)=0.36_r8
    CFOPE(ielmc,infoliar,ilignin,NZ,NY,NX)=0.16_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(IBTYP(NZ,NY,NX).EQ.1.OR.IBTYP(NZ,NY,NX).GE.3)THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ,NY,NX)=0.07_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ,NY,NX)=0.34_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ,NY,NX)=0.36_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ,NY,NX)=0.23_r8

    CFOPE(ielmc,infoliar,iprotein,NZ,NY,NX)=0.0_r8
    CFOPE(ielmc,infoliar,icarbhyro,NZ,NY,NX)=0.045_r8
    CFOPE(ielmc,infoliar,icellulos,NZ,NY,NX)=0.660_r8
    CFOPE(ielmc,infoliar,ilignin,NZ,NY,NX)=0.295_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPE(ielmc,ifoliar,iprotein,NZ,NY,NX)=0.07_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ,NY,NX)=0.25_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ,NY,NX)=0.38_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ,NY,NX)=0.30_r8

    CFOPE(ielmc,infoliar,iprotein,NZ,NY,NX)=0.0_r8
    CFOPE(ielmc,infoliar,icarbhyro,NZ,NY,NX)=0.045_r8
    CFOPE(ielmc,infoliar,icellulos,NZ,NY,NX)=0.660_r8
    CFOPE(ielmc,infoliar,ilignin,NZ,NY,NX)=0.295_r8
  ENDIF
!
!     FRACTIONS OF WOODY LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     NON-VASCULAR
!
  IF(IGTYP(NZ,NY,NX).EQ.0)THEN
    CFOPE(ielmc,istalk,iprotein,NZ,NY,NX)=0.07_r8
    CFOPE(ielmc,istalk,icarbhyro,NZ,NY,NX)=0.25_r8
    CFOPE(ielmc,istalk,icellulos,NZ,NY,NX)=0.30_r8
    CFOPE(ielmc,istalk,ilignin,NZ,NY,NX)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
    CFOPE(ielmc,istalk,iprotein,NZ,NY,NX)=0.03_r8
    CFOPE(ielmc,istalk,icarbhyro,NZ,NY,NX)=0.25_r8
    CFOPE(ielmc,istalk,icellulos,NZ,NY,NX)=0.57_r8
    CFOPE(ielmc,istalk,ilignin,NZ,NY,NX)=0.15_r8
!
!     DECIDUOUS AND CONIFEROUS TREES
!
  ELSE
    CFOPE(ielmc,istalk,iprotein,NZ,NY,NX)=0.0_r8
    CFOPE(ielmc,istalk,icarbhyro,NZ,NY,NX)=0.045_r8
    CFOPE(ielmc,istalk,icellulos,NZ,NY,NX)=0.660_r8
    CFOPE(ielmc,istalk,ilignin,NZ,NY,NX)=0.295_r8
  ENDIF
!
!     FRACTIONS OF FINE ROOT LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN PC&E 25:601-608
!
!     NON-VASCULAR
!
  IF(IGTYP(NZ,NY,NX).EQ.0)THEN
    CFOPE(ielmc,iroot,iprotein,NZ,NY,NX)=0.07_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ,NY,NX)=0.25_r8
    CFOPE(ielmc,iroot,icellulos,NZ,NY,NX)=0.30_r8
    CFOPE(ielmc,iroot,ilignin,NZ,NY,NX)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
    CFOPE(ielmc,iroot,iprotein,NZ,NY,NX)=0.057_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ,NY,NX)=0.263_r8
    CFOPE(ielmc,iroot,icellulos,NZ,NY,NX)=0.542_r8
    CFOPE(ielmc,iroot,ilignin,NZ,NY,NX)=0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(IBTYP(NZ,NY,NX).EQ.1.OR.IBTYP(NZ,NY,NX).GE.3)THEN
    CFOPE(ielmc,iroot,iprotein,NZ,NY,NX)=0.059_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ,NY,NX)=0.308_r8
    CFOPE(ielmc,iroot,icellulos,NZ,NY,NX)=0.464_r8
    CFOPE(ielmc,iroot,ilignin,NZ,NY,NX)=0.169_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPE(ielmc,iroot,iprotein,NZ,NY,NX)=0.059_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ,NY,NX)=0.308_r8
    CFOPE(ielmc,iroot,icellulos,NZ,NY,NX)=0.464_r8
    CFOPE(ielmc,iroot,ilignin,NZ,NY,NX)=0.169_r8
  ENDIF
!
!     COARSE WOODY LITTER FROM BOLES AND ROOTS
!
  CFOPE(ielmc,icwood,iprotein,NZ,NY,NX)=0.00_r8
  CFOPE(ielmc,icwood,icarbhyro,NZ,NY,NX)=0.045_r8
  CFOPE(ielmc,icwood,icellulos,NZ,NY,NX)=0.660_r8
  CFOPE(ielmc,icwood,ilignin,NZ,NY,NX)=0.295_r8
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
    D100: DO M=1,jskenc
      CNOPCT=CNOPCT+CFOPE(ielmc,N,M,NZ,NY,NX)*CNOPC(M)
      CPOPCT=CPOPCT+CFOPE(ielmc,N,M,NZ,NY,NX)*CPOPC(M)
    ENDDO D100
    D105: DO M=1,jskenc
      CFOPE(ielmn,N,M,NZ,NY,NX)=CFOPE(ielmc,N,M,NZ,NY,NX)*CNOPC(M)/CNOPCT
      CFOPE(ielmp,N,M,NZ,NY,NX)=CFOPE(ielmc,N,M,NZ,NY,NX)*CPOPC(M)/CPOPCT
    ENDDO D105
  ENDDO D110
!
!     CONCURRENT NODE GROWTH
!
!     FNOD=scales node number for perennial vegetation (e.g. trees)
!     NNOD=number of concurrently growing nodes
!
  IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
! deciduous or shallow root
    FNOD(NZ,NY,NX)=1.0_r8
!
    IF(GROUPI(NZ,NY,NX).LE.10)THEN
      NNOD(NZ,NY,NX)=3
    ELSEIF(GROUPI(NZ,NY,NX).LE.15)THEN
      NNOD(NZ,NY,NX)=4
    ELSE
      NNOD(NZ,NY,NX)=5
    ENDIF
  ELSE
    FNOD(NZ,NY,NX)=AMAX1(1.0_r8,0.04_r8/XRLA(NZ,NY,NX))
    NNOD(NZ,NY,NX)=24
  ENDIF
  end associate
  end subroutine PlantLitterFractions
!------------------------------------------------------------------------------------------

  subroutine PFTThermalAcclimation(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
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
  ZTYP(NZ,NY,NX)=ZTYPI(NZ,NY,NX)
  OFFST(NZ,NY,NX)=2.667*(2.5-ZTYP(NZ,NY,NX))
  TCZ(NZ,NY,NX)=TCZD-OFFST(NZ,NY,NX)
  TCX(NZ,NY,NX)=AMIN1(15.0,TCXD-OFFST(NZ,NY,NX))
  IF(ICTYP(NZ,NY,NX).EQ.3)THEN
    IF(DATAP(NZ,NY,NX)(1:4).EQ.'soyb')THEN
      HTC(NZ,NY,NX)=30.0_r8+3.0_r8*ZTYP(NZ,NY,NX)
      SSTX(NZ,NY,NX)=0.002_r8
    ELSE
      HTC(NZ,NY,NX)=27.0_r8+3.0_r8*ZTYP(NZ,NY,NX)
      SSTX(NZ,NY,NX)=0.002_r8
    ENDIF
  ELSE
    HTC(NZ,NY,NX)=27.0_r8+3.0_r8*ZTYP(NZ,NY,NX)
    SSTX(NZ,NY,NX)=0.005_r8
  ENDIF
  end subroutine PFTThermalAcclimation
!------------------------------------------------------------------------------------------

  subroutine InitDimensionsandUptake(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  INTEGER :: L,N,NR
!
!     SEED CHARACTERISTICS
!
  call calc_seed_geometry(SeedCMass(NZ,NY,NX),SeedVolume(NZ,NY,NX),&
    SeedLength(NZ,NY,NX),SeedArea(NZ,NY,NX))

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
  SeedinDepth(NZ,NY,NX)=PlantinDepth(NZ,NY,NX)
  D9795: DO L=NU(NY,NX),NL(NY,NX)
    IF(SeedinDepth(NZ,NY,NX).GE.CumSoilThickness(L-1,NY,NX) &
      .AND.SeedinDepth(NZ,NY,NX).LT.CumSoilThickness(L,NY,NX))THEN
      NGTopRootLayer(NZ,NY,NX)=L
      NIXBotRootLayer(NZ,NY,NX)=L
      D9790: DO NR=1,pltpar%JRS
        NINR(NR,NZ,NY,NX)=L
      ENDDO D9790
    ENDIF
  ENDDO D9795  
  CNRTS(NZ,NY,NX)=CNRT(NZ,NY,NX)*BiomGrowthYieldRoot(NZ,NY,NX)
  CPRTS(NZ,NY,NX)=CPRT(NZ,NY,NX)*BiomGrowthYieldRoot(NZ,NY,NX)
  MaxPrimRootRadius(2,NZ,NY,NX)=5.0E-06
  MaxSecndRootRadius(2,NZ,NY,NX)=5.0E-06
  RootPorosity(2,NZ,NY,NX)=RootPorosity(1,NZ,NY,NX)
  UPMXZH(2,NZ,NY,NX)=UPMXZH(1,NZ,NY,NX)
  UPKMZH(2,NZ,NY,NX)=UPKMZH(1,NZ,NY,NX)
  UPMNZH(2,NZ,NY,NX)=UPMNZH(1,NZ,NY,NX)
  UPMXZO(2,NZ,NY,NX)=UPMXZO(1,NZ,NY,NX)
  UPKMZO(2,NZ,NY,NX)=UPKMZO(1,NZ,NY,NX)
  UPMNZO(2,NZ,NY,NX)=UPMNZO(1,NZ,NY,NX)
  UPMXPO(2,NZ,NY,NX)=UPMXPO(1,NZ,NY,NX)
  UPKMPO(2,NZ,NY,NX)=UPKMPO(1,NZ,NY,NX)
  UPMNPO(2,NZ,NY,NX)=UPMNPO(1,NZ,NY,NX)
  RSRR(2,NZ,NY,NX)=1.0E+04
  RSRA(2,NZ,NY,NX)=1.0E+12
!
!     RootPoreTortu4Gas=tortuosity for gas transport
!     RRADP=path length for radial diffusion within root (m)
!     DMVL=volume:C ratio (m3 g-1)
!     PrimRootSpecLen,SecndRootSpecLen=specific primary,secondary root length (m g-1)
!     PrimRootXSecArea,SecndRootXSecArea=specific primary,secondary root area (m2 g-1)
!
  D500: DO N=1,2
    RootPoreTortu4Gas(N,NZ,NY,NX)=RootPorosity(N,NZ,NY,NX)**1.33_r8
    RRADP(N,NZ,NY,NX)=LOG(1.0_r8/SQRT(AMAX1(0.01_r8,RootPorosity(N,NZ,NY,NX))))
    DMVL(N,NZ,NY,NX)=ppmc/(0.05_r8*(1.0-RootPorosity(N,NZ,NY,NX)))
    PrimRootSpecLen(N,NZ,NY,NX)=DMVL(N,NZ,NY,NX)/(PICON*MaxPrimRootRadius(N,NZ,NY,NX)**2)
    SecndRootSpecLen(N,NZ,NY,NX)=DMVL(N,NZ,NY,NX)/(PICON*MaxSecndRootRadius(N,NZ,NY,NX)**2)
    MaxPrimRootRadius1(N,NZ,NY,NX)=MaxPrimRootRadius(N,NZ,NY,NX)
!    2*SQRT(0.25*(1.0-RootPorosity(N,NZ,NY,NX)))
    MaxSecndRootRadius1(N,NZ,NY,NX)=MaxSecndRootRadius(N,NZ,NY,NX)
!    2*SQRT(0.25*(1.0-RootPorosity(N,NZ,NY,NX)))
    PrimRootXSecArea(N,NZ,NY,NX)=PICON*MaxPrimRootRadius1(N,NZ,NY,NX)**2._r8
    SecndRootXSecArea(N,NZ,NY,NX)=PICON*MaxSecndRootRadius1(N,NZ,NY,NX)**2._r8
  ENDDO D500 
  end subroutine InitDimensionsandUptake
!------------------------------------------------------------------------------------------

  subroutine InitPlantPhenoMorphoBio(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  integer :: K,L,M,N,NB
!
!     INITIALIZE PLANT PHENOLOGY
!
!     PP=population (grid cell-1)
!
  pftPlantPopulation(NZ,NY,NX)=PPX(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
  doInitPlant(NZ,NY,NX)=ifalse
  iPlantShootState(NZ,NY,NX)=iDead
  iPlantRootState(NZ,NY,NX)=iDead
  NBT(NZ,NY,NX)=0
  NumOfBranches_pft(NZ,NY,NX)=0
  HypoctoylHeight(NZ,NY,NX)=0._r8
  CanopyHeight(NZ,NY,NX)=0._r8
  D10: DO NB=1,JBR
    doInitLeafOut(NB,NZ,NY,NX)=0
    doPlantLeafOut(NB,NZ,NY,NX)=iDisable
    doPlantLeaveOff(NB,NZ,NY,NX)=iDisable
    IFLGR(NB,NZ,NY,NX)=0
    IFLGQ(NB,NZ,NY,NX)=0
    GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
    PSTG(NB,NZ,NY,NX)=XTLI(NZ,NY,NX)
    PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
    PSTGF(NB,NZ,NY,NX)=0._r8
    VSTG(NB,NZ,NY,NX)=0._r8
    VSTGX(NB,NZ,NY,NX)=0._r8
    KLEAF(NB,NZ,NY,NX)=1
    KLEAFX(NB,NZ,NY,NX)=1
    KVSTG(NB,NZ,NY,NX)=1
    KVSTGN(NB,NZ,NY,NX)=0
    GSTGI(NB,NZ,NY,NX)=0._r8
    GSTGF(NB,NZ,NY,NX)=0._r8
    TGSTGI(NB,NZ,NY,NX)=0._r8
    TGSTGF(NB,NZ,NY,NX)=0._r8
    VRNY(NB,NZ,NY,NX)=0._r8
    VRNZ(NB,NZ,NY,NX)=0._r8
    VRNS(NB,NZ,NY,NX)=VRNY(NB,NZ,NY,NX)
    Hours4LeafOff(NB,NZ,NY,NX)=VRNZ(NB,NZ,NY,NX)
    HourCounter4LeafOut_brch(NB,NZ,NY,NX)=0._r8
    RubiscoActivity_brpft(NB,NZ,NY,NX)=1.0
    FDBKX(NB,NZ,NY,NX)=1.0
    FLG4(NB,NZ,NY,NX)=0
    FLGZ(NB,NZ,NY,NX)=0
    BranchNumber_brchpft(NB,NZ,NY,NX)=0
    iPlantBranchState(NB,NZ,NY,NX)=1
    D15: DO M=1,NumGrowthStages
      iPlantCalendar(M,NB,NZ,NY,NX)=0
    ENDDO D15
  ENDDO D10
!
!     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
!
  WSTR(NZ,NY,NX)=0._r8
  CHILL(NZ,NY,NX)=0._r8
  EPOOL(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  EPOLNB(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WTSHTBE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WTSHEBE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WTSTKBE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WTLFBE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WTRSVBE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WTHSKBE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WTGRBE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WTEARBE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WTNDBE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  RCELX(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  RCESX(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WGSHEXE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8  
  WTSTXBE(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  WGLFEX(1:NumOfPlantChemElements,1:JBR,NZ,NY,NX)=0._r8
  
  D25: DO NB=1,JBR
    CanPBStalkC(NB,NZ,NY,NX)=0._r8
    CanPBLeafShethC(NB,NZ,NY,NX)=0._r8
    GRNXB(NB,NZ,NY,NX)=0._r8
    GRNOB(NB,NZ,NY,NX)=0._r8
    GRWTB(NB,NZ,NY,NX)=0._r8
    CanopyBranchLeafA_pft(NB,NZ,NY,NX)=0._r8
    RNH3B(NB,NZ,NY,NX)=0._r8
    ARLFZ(NB,NZ,NY,NX)=0._r8
    CanPBranchHeight(NB,NZ,NY,NX)=0._r8
    
    D5: DO L=1,JC
      CanopyBranchStemApft_lyr(L,NB,NZ,NY,NX)=0._r8
      DO N=1,JLI
        StemA_lyrnodbrchpft(N,L,NB,NZ,NY,NX)=0._r8
      enddo
    ENDDO D5
    DO K=0,MaxCanopyNodes
      ARLF(K,NB,NZ,NY,NX)=0._r8
      HTNODE(K,NB,NZ,NY,NX)=0._r8
      HTNODX(K,NB,NZ,NY,NX)=0._r8
      CanPSheathHeight(K,NB,NZ,NY,NX)=0._r8
      WGLFE(1:NumOfPlantChemElements,K,NB,NZ,NY,NX)=0._r8
      WGSHE(1:NumOfPlantChemElements,K,NB,NZ,NY,NX)=0._r8
      WGNODE(1:NumOfPlantChemElements,K,NB,NZ,NY,NX)=0._r8
      WSLF(K,NB,NZ,NY,NX)=0._r8
      WSSHE(K,NB,NZ,NY,NX)=0._r8

      D55: DO L=1,JC
        CanPLNBLA(L,K,NB,NZ,NY,NX)=0._r8
        WGLFLE(1:NumOfPlantChemElements,L,K,NB,NZ,NY,NX)=0._r8
      ENDDO D55
      IF(K.NE.0)THEN
        CPOOL3(K,NB,NZ,NY,NX)=0._r8
        CO2B(K,NB,NZ,NY,NX)=0._r8
        HCOB(K,NB,NZ,NY,NX)=0._r8
        CPOOL4(K,NB,NZ,NY,NX)=0._r8
        D45: DO L=1,JC
          DO N=1,JLI
            LeafA_lyrnodbrchpft(N,L,K,NB,NZ,NY,NX)=0._r8
          enddo
        ENDDO D45
      ENDIF
    enddo
  ENDDO D25
  D35: DO L=1,JC
    CanopyLeafApft_lyr(L,NZ,NY,NX)=0._r8
    CanopyLeafCpft_lyr(L,NZ,NY,NX)=0._r8
    CanopyStemApft_lyr(L,NZ,NY,NX)=0._r8
  ENDDO D35
  CanopyNonstructElements_pft(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  CanopyNonstructElementConc_pft(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  NoduleNonstructCconc_pft(NZ,NY,NX)=0._r8
  CanPShootElmMass(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  WTLFE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  WTSHEE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  WTSTKE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  CanPStalkC(NZ,NY,NX)=0._r8
  WTRSVE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  WTHSKE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  WTEARE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  WTGRE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  WTRTE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  WTRTSE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  WTNDE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
  CanopyLeafShethC_pft(NZ,NY,NX)=0._r8

  CanopyLeafA_pft(NZ,NY,NX)=0._r8
  WTRTA(NZ,NY,NX)=0._r8
  CanopyStemA_pft(NZ,NY,NX)=0._r8
  end subroutine InitPlantPhenoMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitMassBalance(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  integer :: M
  real(r8) :: WTSTDX
  associate(                 &
    icwood => pltpar%icwood  &
  )
!
!     INITIALIZE MASS BALANCE CHECKS
!
  IF(.not.is_restart().AND.is_first_year)THEN
    CARBN(NZ,NY,NX)=0._r8
    TESN0(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
    TCO2T(NZ,NY,NX)=0._r8
    TCO2A(NZ,NY,NX)=0._r8
    TEUPTK(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
    TESNC(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
    TZUPFX(NZ,NY,NX)=0._r8
    RNH3C(NZ,NY,NX)=0._r8
    TNH3C(NZ,NY,NX)=0._r8
    VCO2F(NZ,NY,NX)=0._r8
    VCH4F(NZ,NY,NX)=0._r8
    VOXYF(NZ,NY,NX)=0._r8
    VNH3F(NZ,NY,NX)=0._r8
    VN2OF(NZ,NY,NX)=0._r8
    VPO4F(NZ,NY,NX)=0._r8
    THVSTE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
    HVSTE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
    RSETE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
    ETCanP(NZ,NY,NX)=0._r8
    WTSTGE(1:NumOfPlantChemElements,NZ,NY,NX)=0._r8
    WTSTDX=WTSTDI(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
    D155: DO M=1,jskenc
      WTSTDE(ielmc,M,NZ,NY,NX)=WTSTDX*CFOPE(ielmc,icwood,M,NZ,NY,NX)
      WTSTDE(ielmn,M,NZ,NY,NX)=WTSTDX*CNSTK(NZ,NY,NX)*CFOPE(ielmn,icwood,M,NZ,NY,NX)
      WTSTDE(ielmp,M,NZ,NY,NX)=WTSTDX*CPSTK(NZ,NY,NX)*CFOPE(ielmp,icwood,M,NZ,NY,NX)
      WTSTGE(ielmc,NZ,NY,NX)=WTSTGE(ielmc,NZ,NY,NX)+WTSTDE(ielmc,M,NZ,NY,NX)
      WTSTGE(ielmn,NZ,NY,NX)=WTSTGE(ielmn,NZ,NY,NX)+WTSTDE(ielmn,M,NZ,NY,NX)
      WTSTGE(ielmp,NZ,NY,NX)=WTSTGE(ielmp,NZ,NY,NX)+WTSTDE(ielmp,M,NZ,NY,NX)
    ENDDO D155
  ENDIF
  end associate
  end subroutine InitMassBalance
!------------------------------------------------------------------------------------------

  subroutine InitPlantHeatandWater(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
!
!     INITIALIZE PLANT HEAT AND WATER STATUS
!
!     VHeatCapCanP=canopy heat capacity (MJ m-3 K-1)
!     TCC,TKC=canopy temperature for growth (oC,K)
!     TCG,TKG=canopy temperature for phenology (oC,K)
!     PSICanP,PSICanPOsmo,PSICanPTurg=canopy total,osmotic,turgor water potl(MPa)
!
  VHeatCapCanP(NZ,NY,NX)=cpw*CanPShootElmMass(ielmc,NZ,NY,NX)*10.0E-06
  ENGYX(NZ,NY,NX)=0._r8
  DTKC(NZ,NY,NX)=0._r8
  TCC(NZ,NY,NX)=ATCA(NY,NX)
  TKC(NZ,NY,NX)=units%Celcius2Kelvin(TCC(NZ,NY,NX))
  TCG(NZ,NY,NX)=TCC(NZ,NY,NX)
  TKG(NZ,NY,NX)=units%Celcius2Kelvin(TCG(NZ,NY,NX))
  fTgrowCanP(NZ,NY,NX)=1.0
  PSICanP(NZ,NY,NX)=-1.0E-03
  PSICanPOsmo(NZ,NY,NX)=OSMO(NZ,NY,NX)+PSICanP(NZ,NY,NX)
  PSICanPTurg(NZ,NY,NX)=AZMAX1(PSICanP(NZ,NY,NX)-PSICanPOsmo(NZ,NY,NX))
  PTrans(NZ,NY,NX)=0._r8
  FracPARByCanP(NZ,NY,NX)=0._r8
  end subroutine InitPlantHeatandWater
!------------------------------------------------------------------------------------------

  subroutine InitRootMychorMorphoBio(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ, NY, NX
  integer :: K,L,M,N,NR
  REAL(R8) :: CCO2A
  REAL(R8) :: CCO2P
  REAL(R8) :: COXYA
  REAL(R8) :: COXYP
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
!
!     PSIRoot,PSIRootOSMO,PSIRootTurg=root,myco total,osmotic,turgor water potl(MPa)
!     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
!     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
!
  NRT(NZ,NY,NX)=0
  UPNH4(NZ,NY,NX)=0._r8
  UPNO3(NZ,NY,NX)=0._r8
  UPH2P(NZ,NY,NX)=0._r8
  UPH1P(NZ,NY,NX)=0._r8
  UPNF(NZ,NY,NX)=0._r8
  D40: DO N=1,pltpar%jroots
    D20: DO L=1,NL(NY,NX)
      AllPlantRootH2OUptake_vr(N,L,NZ,NY,NX)=0._r8
      PSIRoot(N,L,NZ,NY,NX)=-0.01
      PSIRootOSMO(N,L,NZ,NY,NX)=OSMO(NZ,NY,NX)+PSIRoot(N,L,NZ,NY,NX)
      PSIRootTurg(N,L,NZ,NY,NX)=AZMAX1(PSIRoot(N,L,NZ,NY,NX)-PSIRootOSMO(N,L,NZ,NY,NX))
      EPOOLR(1:NumOfPlantChemElements,N,L,NZ,NY,NX)=0._r8
      RootNonstructElementConcpft_vr(1:NumOfPlantChemElements,N,L,NZ,NY,NX)=0._r8
      RootProteinConc_pftvr(N,L,NZ,NY,NX)=RootFracRemobilizableBiom(NZ,NY,NX)
      WTRTL(N,L,NZ,NY,NX)=0._r8
      PopPlantRootC_vr(N,L,NZ,NY,NX)=0._r8
      WSRTL(N,L,NZ,NY,NX)=0._r8
      PrimRootXNumL(N,L,NZ,NY,NX)=0._r8
      SecndRootXNumL(N,L,NZ,NY,NX)=0._r8
      RootLenPerP(N,L,NZ,NY,NX)=0._r8
      RootLenDensNLP(N,L,NZ,NY,NX)=0._r8
      RTVLP(N,L,NZ,NY,NX)=0._r8
      RTVLW(N,L,NZ,NY,NX)=0._r8
      PrimRootRadius(N,L,NZ,NY,NX)=MaxPrimRootRadius(N,NZ,NY,NX)
      SecndRootRadius(N,L,NZ,NY,NX)=MaxSecndRootRadius(N,NZ,NY,NX)
      RTARP(N,L,NZ,NY,NX)=0._r8
      AveSecndRootLen(N,L,NZ,NY,NX)=1.0E-03
      RUPNH4(N,L,NZ,NY,NX)=0._r8
      RUPNO3(N,L,NZ,NY,NX)=0._r8
      RUPH2P(N,L,NZ,NY,NX)=0._r8
      RUPH1P(N,L,NZ,NY,NX)=0._r8
      RUPNHB(N,L,NZ,NY,NX)=0._r8
      RUPNOB(N,L,NZ,NY,NX)=0._r8
      RUPH2B(N,L,NZ,NY,NX)=0._r8
      RUPH1B(N,L,NZ,NY,NX)=0._r8
      ROXYP(N,L,NZ,NY,NX)=0._r8
      RUNNHP(N,L,NZ,NY,NX)=0._r8
      RUNNBP(N,L,NZ,NY,NX)=0._r8
      RUNNOP(N,L,NZ,NY,NX)=0._r8
      RUNNXP(N,L,NZ,NY,NX)=0._r8
      RUPP2P(N,L,NZ,NY,NX)=0._r8
      RUPP1P(N,L,NZ,NY,NX)=0._r8
      RUPP2B(N,L,NZ,NY,NX)=0._r8
      RUPP1B(N,L,NZ,NY,NX)=0._r8
      CCO2A=CCO2EI(NY,NX)
      CCO2P=0.030_r8*EXP(-2.621_r8-0.0317_r8*ATCA(NY,NX))*CO2EI(NY,NX)
      trcg_rootml(idg_CO2,N,L,NZ,NY,NX)=CCO2A*RTVLP(N,L,NZ,NY,NX)
      trcs_rootml(idg_CO2,N,L,NZ,NY,NX)=CCO2P*RTVLW(N,L,NZ,NY,NX)
      trcg_air2root_flx_pft_vr(idg_CO2,N,L,NZ,NY,NX)=0._r8
      trcg_Root_DisEvap_flx_vr(idg_CO2,N,L,NZ,NY,NX)=0._r8
      RCO2S(N,L,NZ,NY,NX)=0._r8
      RCO2P(N,L,NZ,NY,NX)=0._r8
      COXYA=AtmGgms(idg_O2,NY,NX)
      COXYP=0.032_r8*EXP(-6.175_r8-0.0211_r8*ATCA(NY,NX))*OXYE(NY,NX)
      trcg_rootml(idg_beg:idg_end-1,N,L,NZ,NY,NX)=0._r8
      trcs_rootml(idg_beg:idg_end-1,N,L,NZ,NY,NX)=0._r8
      trcg_rootml(idg_O2,N,L,NZ,NY,NX)=COXYA*RTVLP(N,L,NZ,NY,NX)
      trcs_rootml(idg_O2,N,L,NZ,NY,NX)=COXYP*RTVLW(N,L,NZ,NY,NX)

      WFR(N,L,NZ,NY,NX)=1.0
      D30: DO NR=1,JRS
        RTN2(N,L,NR,NZ,NY,NX)=0._r8
        PrimRootLen(N,L,NR,NZ,NY,NX)=0._r8
        WTRT1E(1:NumOfPlantChemElements,N,L,NR,NZ,NY,NX)=0._r8
        SecndRootLen(N,L,NR,NZ,NY,NX)=0._r8
        WTRT2E(1:NumOfPlantChemElements,N,L,NR,NZ,NY,NX)=0._r8
        PrimRootDepth(N,NR,NZ,NY,NX)=SeedinDepth(NZ,NY,NX)
        RTWT1E(1:NumOfPlantChemElements,N,NR,NZ,NY,NX)=0._r8
      ENDDO D30
      IF(N.EQ.1)THEN
        D6400: DO K=1,pltpar%NumOfPlantLitrCmplxs
          DO  M=1,jskenc
            ESNC(1:NumOfPlantChemElements,M,K,L,NZ,NY,NX)=0._r8
          enddo
        ENDDO D6400
        EPOOLN(1:NumOfPlantChemElements,L,NZ,NY,NX)=0._r8
        WTNDLE(1:NumOfPlantChemElements,L,NZ,NY,NX)=0._r8
        RUPNF(L,NZ,NY,NX)=0._r8
      ENDIF
    ENDDO D20
  ENDDO D40

  RUPNH4(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)=0._r8
  RUPNHB(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)=0._r8
  RUPH2P(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)=0._r8
  RUPH2B(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)=0._r8
  RootLenDensNLP(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)=0._r8
  end subroutine InitRootMychorMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitSeedMorphoBio(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  REAL(R8) :: FDM
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
!     WTRTL,WTRTD=total root C mass (g)
!     WSRTL=total root protein C mass (g)
!     CPOOLR,ZPOOLR,PPOOLR=C,N,P in root,myco nonstructural pools (g)
!
  SeedCPlanted_pft(NZ,NY,NX)=SeedCMass(NZ,NY,NX)*pftPlantPopulation(NZ,NY,NX)
  WTRVE(ielmc,NZ,NY,NX)=SeedCPlanted_pft(NZ,NY,NX)
  WTRVE(ielmn,NZ,NY,NX)=CNGR(NZ,NY,NX)*WTRVE(ielmc,NZ,NY,NX)
  WTRVE(ielmp,NZ,NY,NX)=CPGR(NZ,NY,NX)*WTRVE(ielmc,NZ,NY,NX)
  WTLFBE(ielmn,1,NZ,NY,NX)=CNGR(NZ,NY,NX)*WTLFBE(ielmc,1,NZ,NY,NX)
  WTLFBE(ielmp,1,NZ,NY,NX)=CPGR(NZ,NY,NX)*WTLFBE(ielmc,1,NZ,NY,NX)
  CanPBLeafShethC(1,NZ,NY,NX)=WTLFBE(ielmc,1,NZ,NY,NX)+WTSHEBE(ielmc,1,NZ,NY,NX)
  CanopyLeafShethC_pft(NZ,NY,NX)=CanopyLeafShethC_pft(NZ,NY,NX)+CanPBLeafShethC(1,NZ,NY,NX)  
  FDM=AMIN1(1.0_r8,0.16_r8-0.045_r8*PSICanP(NZ,NY,NX))
  CanWatP(NZ,NY,NX)=ppmc*CanopyLeafShethC_pft(NZ,NY,NX)/FDM
  WatByPCan(NZ,NY,NX)=0._r8
  EPOOL(ielmn,1,NZ,NY,NX)=CNGR(NZ,NY,NX)*EPOOL(ielmc,1,NZ,NY,NX)
  EPOOL(ielmp,1,NZ,NY,NX)=CPGR(NZ,NY,NX)*EPOOL(ielmc,1,NZ,NY,NX)
  WTRT1E(ielmn,ipltroot,NGTopRootLayer(NZ,NY,NX),1,NZ,NY,NX)=CNGR(NZ,NY,NX)*WTRT1E(ielmc,ipltroot,NGTopRootLayer(NZ,NY,NX),1,NZ,NY,NX)
  WTRT1E(ielmp,ipltroot,NGTopRootLayer(NZ,NY,NX),1,NZ,NY,NX)=CPGR(NZ,NY,NX)*WTRT1E(ielmc,ipltroot,NGTopRootLayer(NZ,NY,NX),1,NZ,NY,NX)
  RTWT1E(ielmn,1,1,NZ,NY,NX)=CNGR(NZ,NY,NX)*RTWT1E(ielmc,1,1,NZ,NY,NX)
  RTWT1E(ielmp,1,1,NZ,NY,NX)=CPGR(NZ,NY,NX)*RTWT1E(ielmc,1,1,NZ,NY,NX)
  WTRTL(ipltroot,NGTopRootLayer(NZ,NY,NX),NZ,NY,NX)=WTRT1E(ielmc,ipltroot,NGTopRootLayer(NZ,NY,NX),1,NZ,NY,NX)
  PopPlantRootC_vr(ipltroot,NGTopRootLayer(NZ,NY,NX),NZ,NY,NX)=WTRT1E(ielmc,ipltroot,NGTopRootLayer(NZ,NY,NX),1,NZ,NY,NX)
  WSRTL(1,NGTopRootLayer(NZ,NY,NX),NZ,NY,NX)=WTRTL(ipltroot,NGTopRootLayer(NZ,NY,NX),NZ,NY,NX)*RootFracRemobilizableBiom(NZ,NY,NX)
  EPOOLR(ielmn,1,NGTopRootLayer(NZ,NY,NX),NZ,NY,NX)=CNGR(NZ,NY,NX)*EPOOLR(ielmc,1,NGTopRootLayer(NZ,NY,NX),NZ,NY,NX)
  EPOOLR(ielmp,1,NGTopRootLayer(NZ,NY,NX),NZ,NY,NX)=CPGR(NZ,NY,NX)*EPOOLR(ielmc,1,NGTopRootLayer(NZ,NY,NX),NZ,NY,NX)
  end subroutine InitSeedMorphoBio

  end module StartqMod
