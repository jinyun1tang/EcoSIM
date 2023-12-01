module StartqsMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use minimathmod, only : AZMAX1
  use EcoSIMConfig
  use PlantAPIData
  use TracerIDMod
  use EcoSiMParDataMod, only : pltpar
  use PlantMathFuncMod
  use UnitMod, only : units
  use PlantMathFuncMod
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: StartPlants
  contains

  SUBROUTINE StartPlants(NZ1Q,NZ2Q)
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
    IsPlantActive   => plt_pheno%IsPlantActive      &
  )
!
!     INITIALIZE SHOOT GROWTH VARIABLES
!
!     IsPlantActive=PFT flag:0=not active,1=active
!     iYearPlanting,iDayPlanting,iYearPlantHarvest,iDayPlantHarvest=year,day of planting,arvesting
!     PPI,PPX=initial,current population (m-2)
!     CF,ClumpFactort0=current,initial clumping factor
!     MaxCanPStomaResistH2O=cuticular resistance to water (h m-1)
!     CO2CuticleResist_pft=cuticular resistance to CO2 (s m-1)
!     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios
!     RootFracRemobilizableBiom=maximum root protein concentration (g g-1)
!     O2I=intercellular O2 concentration in C3,C4 PFT (umol mol-1)
!

      NZ2X=MIN(NZ2Q,NP)
      D9985: DO NZ=NZ1Q,NZ2X

        IF(IsPlantActive(NZ).EQ.iPlantIsDormant)THEN

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
      ENDDO D9985

      DO NZ=NZ1Q,NZ2X
        ZEROP(NZ)=ZERO*pftPlantPopulation(NZ)
        ZEROQ(NZ)=ZERO*pftPlantPopulation(NZ)/AREA3(NU)
        ZEROL(NZ)=ZERO*pftPlantPopulation(NZ)*1.0E+06_r8
      ENDDO  
!
!     FILL OUT UNUSED ARRAYS
!
      D9986: DO NZ=NP+1,JP1
        plt_bgcr%SurfLitrfallChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
        plt_bgcr%LitrfallChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
        plt_biom%StandingDeadChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
        D6401: DO L=1,NL
          DO  K=1,pltpar%NumOfPlantLitrCmplxs
            plt_bgcr%LitterFallChemElmnt_pftvr(1:NumOfPlantChemElmnts,1:jsken,K,L,NZ)=0._r8
          enddo
        ENDDO D6401
      ENDDO D9986
  RETURN
  end associate
  END subroutine StartPlants
!------------------------------------------------------------------------------------------

  subroutine InitShootGrowth(NZ)

  implicit none
  integer, intent(in) :: NZ
  associate(                         &
    IDAYX  =>  plt_distb%IDAYX , &
    IYRY   =>  plt_distb%IYRY  , &
    IYRX   =>  plt_distb%IYRX  , &
    iYearPlanting   =>  plt_distb%iYearPlanting  , &
    iYearPlantHarvest   =>  plt_distb%iYearPlantHarvest  , &
    iDayPlantHarvest  =>  plt_distb%iDayPlantHarvest , &
    iDayPlanting  =>  plt_distb%iDayPlanting , &
    IDAYY  =>  plt_distb%IDAYY , &
    RSMX   =>  plt_photo%RSMX  , &
    PPI    =>  plt_site%PPI    , &
    PPX    =>  plt_site%PPX    , &
    PPZ    =>  plt_site%PPZ    , &
    RootFracRemobilizableBiom  =>  plt_allom%RootFracRemobilizableBiom , &
    rCNNonstructRemob_pft  =>  plt_allom%rCNNonstructRemob_pft , &
    rCPNonstructRemob_pft   =>  plt_allom%rCPNonstructRemob_pft  , &
    RootrNC_pft  =>  plt_allom%RootrNC_pft , &
    RootrPC_pft  =>  plt_allom%RootrPC_pft , &
    O2I    =>  plt_photo%O2I   , &
    CO2CuticleResist_pft   =>  plt_photo%CO2CuticleResist_pft  , &
    iPlantPhotosynthesisType =>  plt_photo%iPlantPhotosynthesisType, &
    MaxCanPStomaResistH2O   =>  plt_photo%MaxCanPStomaResistH2O  , &
    ClumpFactort0    =>  plt_morph%ClumpFactort0   , &
    ClumpFactor    =>  plt_morph%ClumpFactor   , &
    NumRootAxes_pft   =>  plt_morph%NumRootAxes_pft    &
  )
  iYearPlanting(NZ)=IYRX(NZ)
  iDayPlanting(NZ)=IDAYX(NZ)
  iYearPlantHarvest(NZ)=IYRY(NZ)
  iDayPlantHarvest(NZ)=IDAYY(NZ)
  PPI(NZ)=PPZ(NZ)
  PPX(NZ)=PPI(NZ)
  ClumpFactor(NZ)=ClumpFactort0(NZ)

  MaxCanPStomaResistH2O(NZ)=RSMX(NZ)/3600.0_r8
  CO2CuticleResist_pft(NZ)=RSMX(NZ)*1.56_r8
  rCNNonstructRemob_pft(NZ)=2.5_r8
  rCPNonstructRemob_pft(NZ)=25.0_r8
  RootFracRemobilizableBiom(NZ)=AMIN1(RootrNC_pft(NZ)*rCNNonstructRemob_pft(NZ),RootrPC_pft(NZ)*rCPNonstructRemob_pft(NZ))
  IF(iPlantPhotosynthesisType(NZ).EQ.ic3_photo)THEN
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
    inonstruct => pltpar%inonstruct, &
    ifoliar  => pltpar%ifoliar , &
    inonfoliar => pltpar%inonfoliar, &
    istalk   => pltpar%istalk  , &
    iroot    => pltpar%iroot   , &
    icwood   => pltpar%icwood  , &
    iprotein  => pltpar%iprotein  ,&
    icarbhyro => pltpar%icarbhyro ,&
    icellulos => pltpar%icellulos ,&
    ilignin  =>  pltpar%ilignin   , &
    NumLitterGroups  => pltpar%NumLitterGroups , &
    RefLeafAppearRate   =>  plt_pheno%RefLeafAppearRate      , &
    iPlantTurnoverPattern =>  plt_pheno%iPlantTurnoverPattern    , &
    iPlantMorphologyType =>  plt_pheno%iPlantMorphologyType    , &
    MatureGroup_pft=>  plt_pheno%MatureGroup_pft   , &
    CFOPE  =>  plt_soilchem%CFOPE  , &
    FNOD   =>  plt_allom%FNOD      , &
    iPlantNfixType =>  plt_morph%iPlantNfixType    , &
    NumConCurrentGrowinNode  =>  plt_morph%NumConCurrentGrowinNode       &
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
  CFOPE(ielmc,inonstruct,iprotein,NZ)=0.0_r8
  CFOPE(ielmc,inonstruct,icarbhyro,NZ)=0.67_r8
  CFOPE(ielmc,inonstruct,icellulos,NZ)=0.33_r8
  CFOPE(ielmc,inonstruct,ilignin,NZ)=0.0_r8
!
!     NON-VASCULAR (E.G. MOSSES)
!
  IF(is_plant_bryophyte(iPlantMorphologyType(NZ)))THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ)=0.30_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ)=0.38_r8

    CFOPE(ielmc,inonfoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,inonfoliar,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,inonfoliar,icellulos,NZ)=0.30_r8
    CFOPE(ielmc,inonfoliar,ilignin,NZ)=0.38_r8
!
!     LEGUMES
!
  ELSEIF(is_plant_N2fix(iPlantNfixType(NZ)))THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ)=0.16_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ)=0.38_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ)=0.34_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ)=0.12_r8

    CFOPE(ielmc,inonfoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,inonfoliar,icarbhyro,NZ)=0.41_r8
    CFOPE(ielmc,inonfoliar,icellulos,NZ)=0.37_r8
    CFOPE(ielmc,inonfoliar,ilignin,NZ)=0.15_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType(NZ))))THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ)=0.08_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ)=0.41_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ)=0.36_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ)=0.15_r8

    CFOPE(ielmc,inonfoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,inonfoliar,icarbhyro,NZ)=0.41_r8
    CFOPE(ielmc,inonfoliar,icellulos,NZ)=0.36_r8
    CFOPE(ielmc,inonfoliar,ilignin,NZ)=0.16_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(iPlantTurnoverPattern(NZ).EQ.1.OR.iPlantTurnoverPattern(NZ).GE.3)THEN
    CFOPE(ielmc,ifoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ)=0.34_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ)=0.36_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ)=0.23_r8

    CFOPE(ielmc,inonfoliar,iprotein,NZ)=0.0_r8
    CFOPE(ielmc,inonfoliar,icarbhyro,NZ)=0.045_r8
    CFOPE(ielmc,inonfoliar,icellulos,NZ)=0.660_r8
    CFOPE(ielmc,inonfoliar,ilignin,NZ)=0.295_r8
!
!     CONIFEROUS TREES
!
  ELSE
    CFOPE(ielmc,ifoliar,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,ifoliar,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,ifoliar,icellulos,NZ)=0.38_r8
    CFOPE(ielmc,ifoliar,ilignin,NZ)=0.30_r8

    CFOPE(ielmc,inonfoliar,iprotein,NZ)=0.0_r8
    CFOPE(ielmc,inonfoliar,icarbhyro,NZ)=0.045_r8
    CFOPE(ielmc,inonfoliar,icellulos,NZ)=0.660_r8
    CFOPE(ielmc,inonfoliar,ilignin,NZ)=0.295_r8
  ENDIF
!
!     FRACTIONS OF WOODY LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     NON-VASCULAR
!
  IF(is_plant_bryophyte(iPlantMorphologyType(NZ)))THEN
    CFOPE(ielmc,istalk,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,istalk,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,istalk,icellulos,NZ)=0.30_r8
    CFOPE(ielmc,istalk,ilignin,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType(NZ))))THEN
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
  IF(is_plant_bryophyte(iPlantMorphologyType(NZ)))THEN
    CFOPE(ielmc,iroot,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,iroot,icellulos,NZ)=0.30_r8
    CFOPE(ielmc,iroot,ilignin,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType(NZ))))THEN
    CFOPE(ielmc,iroot,iprotein,NZ)=0.057_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ)=0.263_r8
    CFOPE(ielmc,iroot,icellulos,NZ)=0.542_r8
    CFOPE(ielmc,iroot,ilignin,NZ)=0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(iPlantTurnoverPattern(NZ).EQ.1.OR.iPlantTurnoverPattern(NZ).GE.3)THEN
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

  D110: DO N=0,NumLitterGroups
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
!     NumConCurrentGrowinNode=number of concurrently growing nodes
!
  IF(iPlantTurnoverPattern(NZ).EQ.0.OR.(.not.is_plant_treelike(iPlantMorphologyType(NZ))))THEN
    FNOD(NZ)=1.0_r8
    IF(MatureGroup_pft(NZ).LE.10)THEN
      NumConCurrentGrowinNode(NZ)=3
    ELSEIF(MatureGroup_pft(NZ).LE.15)THEN
      NumConCurrentGrowinNode(NZ)=4
    ELSE
      NumConCurrentGrowinNode(NZ)=5
    ENDIF
  ELSE
    FNOD(NZ)=AMAX1(1.0_r8,0.04_r8/RefLeafAppearRate(NZ))
    NumConCurrentGrowinNode(NZ)=24
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
    iPlantPhotosynthesisType   =>  plt_photo%iPlantPhotosynthesisType , &
    HTC      =>  plt_pheno%HTC    , &
    TCelcius4LeafOffHarden      =>  plt_pheno%TCelcius4LeafOffHarden    , &
    TCelciusChill4Leaf     =>  plt_pheno%TCelciusChill4Leaf   , &
    OFFST    =>  plt_pheno%OFFST  , &
    iPlantInitThermoAdaptZone   =>  plt_pheno%iPlantInitThermoAdaptZone , &
    iPlantThermoAdaptZone    =>  plt_pheno%iPlantThermoAdaptZone  , &
    SSTX     =>  plt_pheno%SSTX     &
  )
!
!     PFT THERMAL ACCLIMATION
!
!     ZTYP,iPlantInitThermoAdaptZone=dynamic,initial thermal adaptation zone from PFT file
!     OFFST=shift in Arrhenius curve for thermal adaptation (oC)
!     TCZ,TCelcius4LeafOffHarden=threshold temperature for leafout,leafoff
!     HTC=high temperature threshold for grain number loss (oC)
!     SSTX=sensitivity to HTC (seeds oC-1 above HTC)
!
  iPlantThermoAdaptZone(NZ)=iPlantInitThermoAdaptZone(NZ)
  OFFST(NZ)=2.667_r8*(2.5_r8-iPlantThermoAdaptZone(NZ))
  TCelciusChill4Leaf(NZ)=TCZD-OFFST(NZ)
  TCelcius4LeafOffHarden(NZ)=AMIN1(15.0_r8,TCXD-OFFST(NZ))
  IF(iPlantPhotosynthesisType(NZ).EQ.ic3_photo)THEN
    IF(DATAP(NZ)(1:4).EQ.'soyb')THEN
      HTC(NZ)=30.0_r8+3.0_r8*iPlantThermoAdaptZone(NZ)
      SSTX(NZ)=0.002_r8
    ELSE
      HTC(NZ)=27.0_r8+3.0_r8*iPlantThermoAdaptZone(NZ)
      SSTX(NZ)=0.002_r8
    ENDIF
  ELSE
    HTC(NZ)=27.0_r8+3.0_r8*iPlantThermoAdaptZone(NZ)
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
    RootBiomGrowthYield     =>  plt_allom%RootBiomGrowthYield    , &
    RootrNC_pft    =>  plt_allom%RootrNC_pft   , &
    RootrPC_pft    =>  plt_allom%RootrPC_pft   , &
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
    NGTopRootLayer_pft       =>  plt_morph%NGTopRootLayer_pft      , &
    SeedinDepth    =>  plt_morph%SeedinDepth   , &
    SeedAreaMean_pft     =>  plt_morph%SeedAreaMean_pft    , &
    SeedCMass     =>  plt_morph%SeedCMass    , &
    PlantinDepth   =>  plt_morph%PlantinDepth  , &
    MaxSecndRootRadius1   =>  plt_morph%MaxSecndRootRadius1  , &
    SecndRootXSecArea   =>  plt_morph%SecndRootXSecArea  , &
    MaxPrimRootRadius1   =>  plt_morph%MaxPrimRootRadius1  , &
    PrimRootXSecArea   =>  plt_morph%PrimRootXSecArea  , &
    MaxSecndRootRadius   =>  plt_morph%MaxSecndRootRadius  , &
    MaxPrimRootRadius   =>  plt_morph%MaxPrimRootRadius  , &
    SecndRootSpecLen   =>  plt_morph%SecndRootSpecLen  , &
    NIXBotRootLayer_pft     =>  plt_morph%NIXBotRootLayer_pft    , &
    RSRR     =>  plt_morph%RSRR    , &
    PrimRootSpecLen   =>  plt_morph%PrimRootSpecLen  , &
    RootPorosity    =>  plt_morph%RootPorosity   , &
    RootPoreTortu4Gas    =>  plt_morph%RootPoreTortu4Gas   , &
    RRADP    =>  plt_morph%RRADP   , &
    RootVolPerMassC_pft     =>  plt_morph%RootVolPerMassC_pft    , &
    RSRA     =>  plt_morph%RSRA    , &
    NIXBotRootLayer_rpft     =>  plt_morph%NIXBotRootLayer_rpft    , &
    SeedVolumeMean_pft     =>  plt_morph%SeedVolumeMean_pft    , &
    SeedLengthMean_pft     =>  plt_morph%SeedLengthMean_pft      &
  )
!
!     SEED CHARACTERISTICS
!
!     SeedVolumeMean_pft,SeedLengthMean_pft,SeedAreaMean_pft=seed volume(m3),length(m),AREA3(NU)(m2)
!     SeedCMass=seed C mass (g) from PFT file
!
  call calc_seed_geometry(SeedCMass(NZ),SeedVolumeMean_pft(NZ),SeedLengthMean_pft(NZ),SeedAreaMean_pft(NZ))
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) DIMENSIONS, UPTAKE PARAMETERS
!
!     SeedinDepth=seeding depth(m) from PFT management file
!     CumSoilThickness=depth to soil layer bottom from surface(m)
!     NG,NIX,NIXBotRootLayer_rpft=seeding,upper,lower rooting layer
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
      NGTopRootLayer_pft(NZ)=L
      NIXBotRootLayer_pft(NZ)=L
      D9790: DO NR=1,pltpar%JRS
        NIXBotRootLayer_rpft(NR,NZ)=L
      ENDDO D9790
    ENDIF
  ENDDO D9795
  CNRTS(NZ)=RootrNC_pft(NZ)*RootBiomGrowthYield(NZ)
  CPRTS(NZ)=RootrPC_pft(NZ)*RootBiomGrowthYield(NZ)
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
!     RootVolPerMassC_pft=volume:C ratio (m3 g-1)
!     PrimRootSpecLen,SecndRootSpecLen=specific primary,secondary root length (m g-1)
!     PrimRootXSecArea,SecndRootXSecArea=specific primary,secondary root area (m2 g-1)
!
  D500: DO N=1,2
    RootPoreTortu4Gas(N,NZ)=RootPorosity(N,NZ)**1.33_r8
    RRADP(N,NZ)=LOG(1.0_r8/SQRT(AMAX1(0.01_r8,RootPorosity(N,NZ))))
    RootVolPerMassC_pft(N,NZ)=ppmc/(0.05_r8*(1.0_r8-RootPorosity(N,NZ)))
    PrimRootSpecLen(N,NZ)=RootVolPerMassC_pft(N,NZ)/(PICON*MaxPrimRootRadius(N,NZ)**2._r8)
    SecndRootSpecLen(N,NZ)=RootVolPerMassC_pft(N,NZ)/(PICON*MaxSecndRootRadius(N,NZ)**2._r8)
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
    iPlantMorphologyType  =>  plt_pheno%iPlantMorphologyType , &
    LeafNumberAtFloralInit_brch  =>  plt_pheno%LeafNumberAtFloralInit_brch , &
    MatureGroup_brch  =>  plt_pheno%MatureGroup_brch , &
    KLeafNodeNumber  =>  plt_pheno%KLeafNodeNumber , &
    KLeafNumLowestGrowing_pft  =>  plt_pheno%KLeafNumLowestGrowing_pft , &
    NodeNumNormByMatgrp_brch   =>  plt_pheno%NodeNumNormByMatgrp_brch  , &
    ReprodNodeNumNormByMatrgrp_brch   =>  plt_pheno%ReprodNodeNumNormByMatrgrp_brch  , &
    HourWithNoGrainFill_brch    =>  plt_pheno%HourWithNoGrainFill_brch   , &
    HoursDoingRemob_brch   =>  plt_pheno%HoursDoingRemob_brch  , &
    Hours4LenthenPhotoPeriod    =>  plt_pheno%Hours4LenthenPhotoPeriod   , &
    Hours4LeafOff    =>  plt_pheno%Hours4LeafOff   , &
    iPlantBranchState   =>  plt_pheno%iPlantBranchState  , &
    Hours4ShortenPhotoPeriod    =>  plt_pheno%Hours4ShortenPhotoPeriod   , &
    iPlantCalendar   =>  plt_pheno%iPlantCalendar  , &
    Hours4Leafout    =>  plt_pheno%Hours4Leafout   , &
    HourCounter4LeafOut_brch    =>  plt_pheno%HourCounter4LeafOut_brch   , &
    HoursCanopyPSITooLow    =>  plt_pheno%HoursCanopyPSITooLow   , &
    MatureGroup_pft =>  plt_pheno%MatureGroup_pft, &
    PetioleChemElmntRemobFlx_brch   =>  plt_pheno%PetioleChemElmntRemobFlx_brch  , &
    TotalReprodNodeNumNormByMatrgrp_brch  =>  plt_pheno%TotalReprodNodeNumNormByMatrgrp_brch , &
    TotalNodeNumNormByMatgrp_brch  =>  plt_pheno%TotalNodeNumNormByMatgrp_brch , &
    C4PhotosynDowreg_brch   =>  plt_photo%C4PhotosynDowreg_brch  , &
    CPOOL3  =>  plt_photo%CPOOL3 , &
    CPOOL4  =>  plt_photo%CPOOL4 , &
    CHILL   =>  plt_photo%CHILL  , &
    RubiscoActivity_brpft    =>  plt_photo%RubiscoActivity_brpft   , &
    CMassHCO3BundleSheath_node    =>  plt_photo%CMassHCO3BundleSheath_node   , &
    CMassCO2BundleSheath_node    =>  plt_photo%CMassCO2BundleSheath_node   , &
    NumOfLeaves_brch   =>  plt_morph%NumOfLeaves_brch  , &
    CanopyHeight      =>  plt_morph%CanopyHeight     , &
    KLeafNumber_brch  =>  plt_morph%KLeafNumber_brch , &
    XTLI    =>  plt_morph%XTLI   , &
    BranchNumber_pft     =>  plt_morph%BranchNumber_pft    , &
    ShootNodeNumber   =>  plt_morph%ShootNodeNumber  , &
    CanopyLeafA_pft   =>  plt_morph%CanopyLeafA_pft  , &
    CanopyBranchStemApft_lyr   =>  plt_morph%CanopyBranchStemApft_lyr  , &
    StemA_lyrnodbrchpft   =>  plt_morph%StemA_lyrnodbrchpft  , &
    CanopyStemA_pft   =>  plt_morph%CanopyStemA_pft  , &
    NodeNumberAtAnthesis   =>  plt_morph%NodeNumberAtAnthesis  , &
    PotentialSeedSites_brch   =>  plt_morph%PotentialSeedSites_brch  , &
    InternodeHeightLive_brch  =>  plt_morph%InternodeHeightLive_brch , &
    SeedNumberSet_brch   =>  plt_morph%SeedNumberSet_brch  , &
    LeafAreaLive_brch   =>  plt_morph%LeafAreaLive_brch  , &
    LeafAreaDying_brch   =>  plt_morph%LeafAreaDying_brch  , &
    CanopyLeafAreaByLayer_pft   =>  plt_morph%CanopyLeafAreaByLayer_pft  , &
    CanopyLeafApft_lyr   =>  plt_morph%CanopyLeafApft_lyr  , &
    CanopyStemApft_lyr   =>  plt_morph%CanopyStemApft_lyr  , &
    LeafA_lyrnodbrchpft    =>  plt_morph%LeafA_lyrnodbrchpft   , &
    LeafAreaNode_brch   =>  plt_morph%LeafAreaNode_brch  , &
    CanPBranchHeight  =>  plt_morph%CanPBranchHeight , &
    HypoctoylHeight   =>  plt_morph%HypoctoylHeight  , &
    BranchNumber_brch    =>  plt_morph%BranchNumber_brch   , &
    NodeNumberToInitFloral   =>  plt_morph%NodeNumberToInitFloral  , &
    KLEAFX  =>  plt_morph%KLEAFX , &
    NumOfBranches_pft     =>  plt_morph%NumOfBranches_pft      &
  )
!
!     INITIALIZE PLANT PHENOLOGY
!
!     PP=population (grid cell-1)
!
  pftPlantPopulation(NZ)=PPX(NZ)*AREA3(NU)
  plt_pheno%doInitPlant(NZ)=ifalse
  plt_pheno%iPlantShootState(NZ)=iDead
  plt_pheno%iPlantRootState(NZ)=iDead
  BranchNumber_pft(NZ)=0
  NumOfBranches_pft(NZ)=0
  HypoctoylHeight(NZ)=0._r8
  CanopyHeight(NZ)=0._r8
  D10: DO NB=1,MaxNumBranches
    plt_pheno%doInitLeafOut_brch(NB,NZ)=0
    plt_pheno%doPlantLeafOut_brch(NB,NZ)=iDisable
    plt_pheno%doPlantLeaveOff_brch(NB,NZ)=iDisable
    plt_pheno%Prep4Literfall_brch(NB,NZ)=ifalse
    plt_pheno%Hours4LiterfalAftMature_brch(NB,NZ)=0
    MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
    ShootNodeNumber(NB,NZ)=XTLI(NZ)
    NodeNumberToInitFloral(NB,NZ)=ShootNodeNumber(NB,NZ)
    NodeNumberAtAnthesis(NB,NZ)=0._r8
    NumOfLeaves_brch(NB,NZ)=0._r8
    LeafNumberAtFloralInit_brch(NB,NZ)=0._r8
    KLeafNumber_brch(NB,NZ)=1
    KLEAFX(NB,NZ)=1
    KLeafNodeNumber(NB,NZ)=1
    KLeafNumLowestGrowing_pft(NB,NZ)=0
    NodeNumNormByMatgrp_brch(NB,NZ)=0._r8
    ReprodNodeNumNormByMatrgrp_brch(NB,NZ)=0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ)=0._r8
    TotalReprodNodeNumNormByMatrgrp_brch(NB,NZ)=0._r8
    Hours4LenthenPhotoPeriod(NB,NZ)=0._r8
    Hours4ShortenPhotoPeriod(NB,NZ)=0._r8
    Hours4Leafout(NB,NZ)=Hours4LenthenPhotoPeriod(NB,NZ)
    Hours4LeafOff(NB,NZ)=Hours4ShortenPhotoPeriod(NB,NZ)
    HourCounter4LeafOut_brch(NB,NZ)=0._r8
    RubiscoActivity_brpft(NB,NZ)=1.0
    C4PhotosynDowreg_brch(NB,NZ)=1.0
    HourWithNoGrainFill_brch(NB,NZ)=0
    HoursDoingRemob_brch(NB,NZ)=0
    BranchNumber_brch(NB,NZ)=0
    plt_pheno%iPlantBranchState(NB,NZ)=iDead
    D15: DO M=1,pltpar%NumGrowthStages
      iPlantCalendar(M,NB,NZ)=0
    ENDDO D15
  ENDDO D10
!
!     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
!
  HoursCanopyPSITooLow(NZ)=0._r8
  CHILL(NZ)=0._r8
  plt_biom%NonstructElmnt_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%NoduleNonstructElmnt_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%PetioleChemElmnts_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%ShootChemElmnt_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%StalkChemElmnts_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%LeafChemElmnts_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%ReserveChemElmnts_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%HuskChemElmnts_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%GrainChemElmnts_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%EarChemElmnts_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%CanopyNoduleChemElmnt_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_pheno%LeafElmntRemobFlx_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%LeafChemElmntRemob_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_pheno%PetioleChemElmntRemobFlx_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%PetioleChemElmntRemob_brch(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8
  plt_biom%BranchStalkChemElmnts_pft_pft(1:NumOfPlantChemElmnts,1:MaxNumBranches,NZ)=0._r8

  D25: DO NB=1,MaxNumBranches
    plt_biom%StalkBiomassC_brch(NB,NZ)=0._r8
    plt_biom%LeafPetioleBiomassC_brch(NB,NZ)=0._r8
    PotentialSeedSites_brch(NB,NZ)=0._r8
    SeedNumberSet_brch(NB,NZ)=0._r8
    plt_allom%GrainSeedBiomCMean_brch(NB,NZ)=0._r8
    LeafAreaLive_brch(NB,NZ)=0._r8
    plt_rbgc%RNH3B(NB,NZ)=0._r8
    LeafAreaDying_brch(NB,NZ)=0._r8
    CanPBranchHeight(NB,NZ)=0._r8
    D5: DO L=1,NumOfCanopyLayers1
      CanopyBranchStemApft_lyr(L,NB,NZ)=0._r8
      DO N=1,JLI1
        StemA_lyrnodbrchpft(N,L,NB,NZ)=0._r8
      enddo
    ENDDO D5
    DO K=0,MaxNodesPerBranch1
      LeafAreaNode_brch(K,NB,NZ)=0._r8
      InternodeHeightLive_brch(K,NB,NZ)=0._r8
      plt_morph%InternodeHeightDying_brch(K,NB,NZ)=0._r8
      plt_morph%PetioleLengthNode_brch(K,NB,NZ)=0._r8
      plt_biom%LeafProteinCNode_brch(K,NB,NZ)=0._r8
      plt_biom%PetioleProteinCNode_brch(K,NB,NZ)=0._r8
      plt_biom%LeafChemElmntNode_brch(1:NumOfPlantChemElmnts,K,NB,NZ)=0._r8
      plt_biom%InternodeChemElmnt_brch(1:NumOfPlantChemElmnts,K,NB,NZ)=0._r8
      plt_biom%PetioleElmntNode_brch(1:NumOfPlantChemElmnts,K,NB,NZ)=0._r8


      D55: DO L=1,NumOfCanopyLayers1
        CanopyLeafAreaByLayer_pft(L,K,NB,NZ)=0._r8
        plt_biom%LeafChemElmntByLayer_pft(1:NumOfPlantChemElmnts,L,K,NB,NZ)=0._r8
      ENDDO D55
      IF(K.NE.0)THEN
        CPOOL3(K,NB,NZ)=0._r8
        CMassCO2BundleSheath_node(K,NB,NZ)=0._r8
        CMassHCO3BundleSheath_node(K,NB,NZ)=0._r8
        CPOOL4(K,NB,NZ)=0._r8
        D45: DO L=1,NumOfCanopyLayers1
          DO N=1,JLI1
            LeafA_lyrnodbrchpft(N,L,K,NB,NZ)=0._r8
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
  plt_biom%CanopyNonstructElements_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%CanopyNonstructElementConc_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%NoduleNonstructCconc_pft(NZ)=0._r8
  plt_biom%ShootChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%LeafChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%PetioleChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%StalkChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%CanopyStalkC_pft(NZ)=0._r8
  plt_biom%ReserveChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%HuskChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%EarChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%GrainChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%RootStructChemElmnt_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%NoduleChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
  plt_biom%CanopyLeafShethC_pft(NZ)=0._r8
  CanopyLeafA_pft(NZ)=0._r8
  plt_biom%RootBiomCPerPlant_pft(NZ)=0._r8
  CanopyStemA_pft(NZ)=0._r8
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
    StandingDeadKCompChemElmnts_pft => plt_biom%StandingDeadKCompChemElmnts_pft    , &
    StandingDeadChemElmnts_pft => plt_biom%StandingDeadChemElmnts_pft    , &
    StandingDeadInitC_pft => plt_biom%StandingDeadInitC_pft    , &
    RootBiomCPerPlant_pft  => plt_biom%RootBiomCPerPlant_pft     , &
    rNCStalk_pft => plt_allom%rNCStalk_pft   , &
    rPCStalk_pft  => plt_allom%rPCStalk_pft    , &
    PlantExudChemElmntCum_pft => plt_rbgc%PlantExudChemElmntCum_pft    , &
    NH3EmiCum_pft  => plt_bgcr%NH3EmiCum_pft     , &
    RNH3C  => plt_bgcr%RNH3C     , &
    SurfLitrfallChemElmnts_pft  => plt_bgcr%SurfLitrfallChemElmnts_pft     , &
    GrossCO2Fix_pft  => plt_bgcr%GrossCO2Fix_pft     , &
    CanopyPlusNoduRespC_pft  => plt_bgcr%CanopyPlusNoduRespC_pft     , &
    GrossResp_pft  => plt_bgcr%GrossResp_pft     , &
    PlantN2FixCum_pft => plt_bgcr%PlantN2FixCum_pft    , &
    LitrfallChemElmnts_pft  => plt_bgcr%LitrfallChemElmnts_pft     , &
    CanopyStemA_pft  => plt_morph%CanopyStemA_pft    , &
    icwood => pltpar%icwood      , &
    RSETE  => plt_pheno%RSETE      &

  )
!
!     INITIALIZE MASS BALANCE CHECKS
!
  IF(.not.is_restart().AND.is_first_year)THEN
    GrossCO2Fix_pft(NZ)=0._r8
    SurfLitrfallChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
    GrossResp_pft(NZ)=0._r8
    CanopyPlusNoduRespC_pft(NZ)=0._r8
    PlantExudChemElmntCum_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
    LitrfallChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8

    PlantN2FixCum_pft(NZ)=0._r8
    RNH3C(NZ)=0._r8
    NH3EmiCum_pft(NZ)=0._r8
    plt_distb%CO2ByFire_pft(NZ)=0._r8
    plt_distb%CH4ByFire_pft(NZ)=0._r8
    plt_distb%VOXYF(NZ)=0._r8
    plt_distb%VNH3F(NZ)=0._r8
    plt_distb%VN2OF(NZ)=0._r8
    plt_distb%VPO4F(NZ)=0._r8
    plt_distb%THVSTE(1:NumOfPlantChemElmnts,NZ)=0._r8
    plt_distb%HVSTE(1:NumOfPlantChemElmnts,NZ)=0._r8
    RSETE(1:NumOfPlantChemElmnts,NZ)=0._r8
    ETCanP(NZ)=0._r8
    StandingDeadChemElmnts_pft(1:NumOfPlantChemElmnts,NZ)=0._r8
    WTSTDX=StandingDeadInitC_pft(NZ)*AREA3(NU)
    D155: DO M=1,jsken
      StandingDeadKCompChemElmnts_pft(ielmc,M,NZ)=WTSTDX*CFOPE(ielmc,icwood,M,NZ)
      StandingDeadKCompChemElmnts_pft(ielmn,M,NZ)=WTSTDX*rNCStalk_pft(NZ)*CFOPE(ielmn,icwood,M,NZ)
      StandingDeadKCompChemElmnts_pft(ielmp,M,NZ)=WTSTDX*rPCStalk_pft(NZ)*CFOPE(ielmp,icwood,M,NZ)
    ENDDO D155
    DO NE=1,NumOfPlantChemElmnts
      StandingDeadChemElmnts_pft(NE,NZ)=StandingDeadChemElmnts_pft(NE,NZ)+sum(StandingDeadKCompChemElmnts_pft(NE,1:jsken,NZ))
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
    TCelciusCanopy     =>  plt_ew%TCelciusCanopy    , &
    TKG      =>  plt_pheno%TKG  , &
    TCG      =>  plt_pheno%TCG  , &
    fTgrowCanP     =>  plt_pheno%fTgrowCanP , &
    ShootChemElmnts_pft   =>  plt_biom%ShootChemElmnts_pft , &
    FracPARbyCanopy_pft    =>  plt_rad%FracPARbyCanopy_pft    &
  )
!
!     INITIALIZE PLANT HEAT AND WATER STATUS
!
!     VHeatCapCanP=canopy heat capacity (MJ m-3 K-1)
!     TCelciusCanopy,TKC=canopy temperature for growth (oC,K)
!     TCG,TKG=canopy temperature for phenology (oC,K)
!     PSICanP,PSICanPOsmo,PSICanPTurg=canopy total,osmotic,turgor water potl(MPa)
!
  VHeatCapCanP(NZ)=cpw*ShootChemElmnts_pft(ielmc,NZ)*10.0E-06
  ENGYX(NZ)=0._r8
  DTKC(NZ)=0._r8
  TCelciusCanopy(NZ)=ATCA
  TKC(NZ)=units%Celcius2Kelvin(TCelciusCanopy(NZ))
  TCG(NZ)=TCelciusCanopy(NZ)
  TKG(NZ)=units%Celcius2Kelvin(TCG(NZ))
  fTgrowCanP(NZ)=1.0
  PSICanP(NZ)=-1.0E-03
  PSICanPOsmo(NZ)=OSMO(NZ)+PSICanP(NZ)
  PSICanPTurg(NZ)=AZMAX1(PSICanP(NZ)-PSICanPOsmo(NZ))
  PTrans(NZ)=0._r8
  FracPARbyCanopy_pft(NZ)=0._r8
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
    RootN2Fix_pvr    =>  plt_bgcr%RootN2Fix_pvr    , &
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
    RootProteinC_pvr    =>  plt_biom%RootProteinC_pvr    , &
    SeedinDepth    =>  plt_morph%SeedinDepth   , &
    RootVH2O_vr   =>  plt_morph%RootVH2O_vr  , &
    RootVolume_vr    =>  plt_morph%RootVolume_vr   , &
    SecndRootRadius    =>  plt_morph%SecndRootRadius   , &
    PrimRootRadius    =>  plt_morph%PrimRootRadius   , &
    MaxPrimRootRadius   =>  plt_morph%MaxPrimRootRadius  , &
    MaxSecndRootRadius   =>  plt_morph%MaxSecndRootRadius  , &
    NumRootAxes_pft     =>  plt_morph%NumRootAxes_pft      &
  )
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
!
!     PSIRoot,PSIRootOSMO,PSIRootTurg=root,myco total,osmotic,turgor water potl(MPa)
!     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
!     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
!
  NumRootAxes_pft(NZ)=0
  plt_rbgc%RootNH4Uptake_pft(NZ)=0._r8
  plt_rbgc%RootNO3Uptake_pft(NZ)=0._r8
  plt_rbgc%RootH2PO4Uptake_pft(NZ)=0._r8
  plt_rbgc%RootHPO4Uptake_pft(NZ)=0._r8
  plt_rbgc%RootN2Fix_pft(NZ)=0._r8
  D40: DO N=1,pltpar%jroots
    D20: DO L=1,NL
      plt_ew%AllPlantRootH2OUptake_vr(N,L,NZ)=0._r8
      PSIRoot(N,L,NZ)=-0.01
      PSIRootOSMO(N,L,NZ)=OSMO(NZ)+PSIRoot(N,L,NZ)
      PSIRootTurg(N,L,NZ)=AZMAX1(PSIRoot(N,L,NZ)-PSIRootOSMO(N,L,NZ))
      plt_biom%RootMycoNonstructElmnt_vr(1:NumOfPlantChemElmnts,N,L,NZ)=0._r8
      plt_biom%RootNonstructElementConcpft_vr(1:NumOfPlantChemElmnts,N,L,NZ)=0._r8
      RootProteinConc_pftvr(N,L,NZ)=RootFracRemobilizableBiom(NZ)
      plt_biom%RootStructBiomC_vr(N,L,NZ)=0._r8
      plt_biom% PopuPlantRootC_vr(N,L,NZ)=0._r8
      RootProteinC_pvr(N,L,NZ)=0._r8
      plt_morph%PrimRootXNumL(N,L,NZ)=0._r8
      plt_morph%SecndRootXNum_pvr(N,L,NZ)=0._r8
      plt_morph%RootLenPerPopu_pvr(N,L,NZ)=0._r8
      plt_morph%RootLenthDensPerPopu_pvr(N,L,NZ)=0._r8
      RootVolume_vr(N,L,NZ)=0._r8
      RootVH2O_vr(N,L,NZ)=0._r8
      PrimRootRadius(N,L,NZ)=MaxPrimRootRadius(N,NZ)
      SecndRootRadius(N,L,NZ)=MaxSecndRootRadius(N,NZ)
      plt_morph%RootAreaPerPlant_vr(N,L,NZ)=0._r8
      plt_morph%AveSecndRootLen(N,L,NZ)=1.0E-03
      plt_rbgc%RootNH4Uptake_pvr(N,L,NZ)=0._r8
      plt_rbgc%RootNO3Uptake_pvr(N,L,NZ)=0._r8
      plt_rbgc%RootH2PO4Uptake_pvr(N,L,NZ)=0._r8
      plt_rbgc%RootHPO4Uptake_pvr(N,L,NZ)=0._r8
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
      trcg_rootml(idg_CO2,N,L,NZ)=CCO2A*RootVolume_vr(N,L,NZ)
      trcs_rootml(idg_CO2,N,L,NZ)=CCO2P*RootVH2O_vr(N,L,NZ)
      plt_rbgc%trcg_air2root_flx_pft_vr(idg_CO2,N,L,NZ)=0._r8
      plt_rbgc%trcg_Root_DisEvap_flx_vr(idg_CO2,N,L,NZ)=0._r8
      plt_rbgc%RCO2S(N,L,NZ)=0._r8
      plt_rbgc%RCO2P(N,L,NZ)=0._r8
      COXYA=COXYE
      COXYP=0.032_r8*EXP(-6.175_r8-0.0211_r8*ATCA)*OXYE
      plt_rbgc%trcg_rootml(idg_O2,N,L,NZ)=COXYA*RootVolume_vr(N,L,NZ)
      plt_rbgc%trcs_rootml(idg_O2,N,L,NZ)=COXYP*RootVH2O_vr(N,L,NZ)
      plt_rbgc%trcg_rootml(idg_beg:idg_end-1,N,L,NZ)=0._r8
      plt_rbgc%trcs_rootml(idg_beg:idg_end-1,N,L,NZ)=0._r8
      plt_rbgc%RootAutoRO2Limiter_pvr(N,L,NZ)=1.0
      D30: DO NR=1,JRS
        plt_morph%SecndRootXNum_rpvr(N,L,NR,NZ)=0._r8
        plt_morph%PrimRootLen(N,L,NR,NZ)=0._r8
        plt_morph%SecndRootLen(N,L,NR,NZ)=0._r8
        plt_morph%PrimRootDepth(N,NR,NZ)=SeedinDepth(NZ)
        plt_biom%Root1stStructChemElmnt_pvr(1:NumOfPlantChemElmnts,N,L,NR,NZ)=0._r8
        plt_biom%Root2ndStructChemElmnt_pvr(1:NumOfPlantChemElmnts,N,L,NR,NZ)=0._r8
        plt_biom%Root1stChemElmnt(1:NumOfPlantChemElmnts,N,NR,NZ)=0._r8
      ENDDO D30
      IF(N.EQ.1)THEN
        D6400: DO K=1,pltpar%NumOfPlantLitrCmplxs
          plt_bgcr%LitterFallChemElmnt_pftvr(1:NumOfPlantChemElmnts,1:jsken,K,L,NZ)=0._r8
        ENDDO D6400
        plt_biom%RootNoduleNonstructElmnt_vr(1:NumOfPlantChemElmnts,L,NZ)=0._r8
        plt_biom%RootNodueChemElmnt_pvr(1:NumOfPlantChemElmnts,L,NZ)=0._r8
        RootN2Fix_pvr(L,NZ)=0._r8
      ENDIF
    ENDDO D20
  ENDDO D40

  plt_rbgc%RootNH4Uptake_pvr(1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPNHB(1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RootH2PO4Uptake_pvr(1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RUPH2B(1:2,NL+1:JZ1,NZ)=0._r8
  plt_morph%RootLenthDensPerPopu_pvr(1:2,NL+1:JZ1,NZ)=0._r8
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
    WatByPCanopy    =>   plt_ew%WatByPCanopy     , &
    CanWatP    =>   plt_ew%CanWatP     , &
    RootFracRemobilizableBiom    =>   plt_allom%RootFracRemobilizableBiom  , &
    CNGR     =>   plt_allom%CNGR   , &
    CPGR     =>   plt_allom%CPGR   , &
    Root1stChemElmnt   =>   plt_biom%Root1stChemElmnt  , &
    SeedCPlanted_pft    =>   plt_biom%SeedCPlanted_pft   , &
    Root1stStructChemElmnt_pvr   =>   plt_biom%Root1stStructChemElmnt_pvr  , &
    LeafPetioleBiomassC_brch    =>   plt_biom%LeafPetioleBiomassC_brch   , &
    CanopyLeafShethC_pft     =>   plt_biom%CanopyLeafShethC_pft    , &
     PopuPlantRootC_vr    =>   plt_biom% PopuPlantRootC_vr   , &
    PetioleChemElmnts_brch =>   plt_biom%PetioleChemElmnts_brch, &
    RootStructBiomC_vr   =>   plt_biom%RootStructBiomC_vr  , &
    RootProteinC_pvr    =>   plt_biom%RootProteinC_pvr   , &
     RootMycoNonstructElmnt_vr   =>   plt_biom%RootMycoNonstructElmnt_vr  , &
    NonstructElmnt_brch   =>   plt_biom%NonstructElmnt_brch  , &
    LeafChemElmnts_brch  =>   plt_biom%LeafChemElmnts_brch , &
    NonstructalChemElmnts_pft    =>   plt_biom%NonstructalChemElmnts_pft   , &
    SeedCMass     =>   plt_morph%SeedCMass   , &
    PrimRootDepth    =>   plt_morph%PrimRootDepth  , &
    NGTopRootLayer_pft      =>   plt_morph%NGTopRootLayer_pft      &
  )
!
!     INITIALIZE SEED MORPHOLOGY AND BIOMASS
!
!     WTRVC,WTRVN,WTRVP=C,N,P in storage reserves (g)
!     WTLFB,WTLFBN,WTLFBP=C,N,P in leaves (g)
!     LeafPetioleBiomassC_brch=C in leaves+petioles (g)
!     FDM-dry matter fraction (g DM C g FM C-1)
!     CanWatP,WatByPCanopy=water volume in,on canopy (m3)
!     CPOOL,ZPOOL,PPOOL=C,N,P in canopy nonstructural pools (g)
!     WTRT1,WTRT1N,WTRT1P=C,N,P in primary root layer (g)
!     RTWT1,RTWT1N,RTWT1P=total C,N,P in primary root (g)
!     RootStructBiomC_vr, PopuPlantRootC_vr=total root C mass (g)
!     RootProteinC_pvr=total root protein C mass (g)
!     CPOOLR,ZPOOLR,PPOOLR=C,N,P in root,myco nonstructural pools (g)
!
  SeedCPlanted_pft(NZ)=SeedCMass(NZ)*pftPlantPopulation(NZ)
  NonstructalChemElmnts_pft(ielmc,NZ)=SeedCPlanted_pft(NZ)
  NonstructalChemElmnts_pft(ielmn,NZ)=CNGR(NZ)*NonstructalChemElmnts_pft(ielmc,NZ)
  NonstructalChemElmnts_pft(ielmp,NZ)=CPGR(NZ)*NonstructalChemElmnts_pft(ielmc,NZ)
  LeafChemElmnts_brch(ielmn,1,NZ)=CNGR(NZ)*LeafChemElmnts_brch(ielmc,1,NZ)
  LeafChemElmnts_brch(ielmp,1,NZ)=CPGR(NZ)*LeafChemElmnts_brch(ielmc,1,NZ)
  LeafPetioleBiomassC_brch(1,NZ)=LeafChemElmnts_brch(ielmc,1,NZ)+PetioleChemElmnts_brch(ielmc,1,NZ)
  CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetioleBiomassC_brch(1,NZ)
  FDM=AMIN1(1.0_r8,0.16_r8-0.045_r8*PSICanP(NZ))
  CanWatP(NZ)=ppmc*CanopyLeafShethC_pft(NZ)/FDM
  WatByPCanopy(NZ)=0._r8
  NonstructElmnt_brch(ielmn,1,NZ)=CNGR(NZ)*NonstructElmnt_brch(ielmc,1,NZ)
  NonstructElmnt_brch(ielmp,1,NZ)=CPGR(NZ)*NonstructElmnt_brch(ielmc,1,NZ)
  Root1stStructChemElmnt_pvr(ielmn,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)=CNGR(NZ)*Root1stStructChemElmnt_pvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  Root1stStructChemElmnt_pvr(ielmp,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)=CPGR(NZ)*Root1stStructChemElmnt_pvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  Root1stChemElmnt(ielmn,1,1,NZ)=CNGR(NZ)*Root1stChemElmnt(ielmc,1,1,NZ)
  Root1stChemElmnt(ielmp,1,1,NZ)=CPGR(NZ)*Root1stChemElmnt(ielmc,1,1,NZ)
  RootStructBiomC_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=Root1stStructChemElmnt_pvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
   PopuPlantRootC_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=Root1stStructChemElmnt_pvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  RootProteinC_pvr(1,NGTopRootLayer_pft(NZ),NZ)=RootStructBiomC_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ)*RootFracRemobilizableBiom(NZ)
   RootMycoNonstructElmnt_vr(ielmn,1,NGTopRootLayer_pft(NZ),NZ)=CNGR(NZ)* RootMycoNonstructElmnt_vr(ielmc,1,NGTopRootLayer_pft(NZ),NZ)
   RootMycoNonstructElmnt_vr(ielmp,1,NGTopRootLayer_pft(NZ),NZ)=CPGR(NZ)* RootMycoNonstructElmnt_vr(ielmc,1,NGTopRootLayer_pft(NZ),NZ)

  end associate
  end subroutine InitSeedMorphoBio

  end module StartqsMod
