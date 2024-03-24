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
    PlantPopulation_pft      => plt_site%PlantPopulation_pft        , &
    NU                       => plt_site%NU        , &
    NP                       => plt_site%NP        , &
    NL                       => plt_site%NL        , &
    ZERO                     => plt_site%ZERO      , &
    AREA3                    => plt_site%AREA3     , &
    ZEROQ                    => plt_rbgc%ZEROQ     , &
    ZEROP                    => plt_biom%ZEROP     , &
    ZEROL                    => plt_biom%ZEROL     , &
    IsPlantActive_pft        => plt_pheno%IsPlantActive_pft      &
  )
!
!     INITIALIZE SHOOT GROWTH VARIABLES
!
!     IsPlantActive_pft=PFT flag:0=not active,1=active
!     iYearPlanting_pft,iDayPlanting_pft,iYearPlantHarvest_pft,iDayPlantHarvest_pft=year,day of planting,arvesting
!     PPI,PPX=initial,current population (m-2)
!     CF,ClumpFactorInit_pft=current,initial clumping factor
!     MaxCanPStomaResistH2O_pft=cuticular resistance to water (h m-1)
!     CO2CuticleResist_pft=cuticular resistance to CO2 (s m-1)
!     CNWS,rCPNonstructRemob_pft=protein:N,protein:P ratios
!     RootFracRemobilizableBiom=maximum root protein concentration (g g-1)
!     O2I=intercellular O2 concentration in C3,C4 PFT (umol mol-1)
!

      NZ2X=MIN(NZ2Q,NP)
      D9985: DO NZ=NZ1Q,NZ2X

        IF(IsPlantActive_pft(NZ).EQ.iPlantIsDormant)THEN

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
        ZEROP(NZ)=ZERO*PlantPopulation_pft(NZ)
        ZEROQ(NZ)=ZERO*PlantPopulation_pft(NZ)/AREA3(NU)
        ZEROL(NZ)=ZERO*PlantPopulation_pft(NZ)*1.0E+06_r8
      ENDDO  
!
!     FILL OUT UNUSED ARRAYS
!
      D9986: DO NZ=NP+1,JP1
        plt_bgcr%SurfLitrfallChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
        plt_bgcr%LitrfallChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
        plt_biom%StandingDeadChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
        D6401: DO L=1,NL
          DO  K=1,pltpar%NumOfPlantLitrCmplxs
            plt_bgcr%LitfalChemElm_pvr(1:NumPlantChemElms,1:jsken,K,L,NZ)=0._r8
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
    IDAYX                      =>  plt_distb%IDAYX , &
    IYRY                       =>  plt_distb%IYRY  , &
    IYRX                       =>  plt_distb%IYRX  , &
    iYearPlanting_pft          =>  plt_distb%iYearPlanting_pft  , &
    iYearPlantHarvest_pft      =>  plt_distb%iYearPlantHarvest_pft  , &
    iDayPlantHarvest_pft       =>  plt_distb%iDayPlantHarvest_pft , &
    iDayPlanting_pft           =>  plt_distb%iDayPlanting_pft , &
    IDAYY                      =>  plt_distb%IDAYY , &
    RSMX                       =>  plt_photo%RSMX  , &
    PPI                        =>  plt_site%PPI    , &
    PPX                        =>  plt_site%PPX    , &
    PPZ                        =>  plt_site%PPZ    , &
    RootFracRemobilizableBiom  =>  plt_allom%RootFracRemobilizableBiom , &
    rCNNonstructRemob_pft      =>  plt_allom%rCNNonstructRemob_pft , &
    rCPNonstructRemob_pft      =>  plt_allom%rCPNonstructRemob_pft  , &
    RootrNC_pft                =>  plt_allom%RootrNC_pft , &
    RootrPC_pft                =>  plt_allom%RootrPC_pft , &
    O2I                        =>  plt_photo%O2I   , &
    CO2CuticleResist_pft       =>  plt_photo%CO2CuticleResist_pft  , &
    iPlantPhotosynthesisType   =>  plt_photo%iPlantPhotosynthesisType, &
    MaxCanPStomaResistH2O_pft  =>  plt_photo%MaxCanPStomaResistH2O_pft  , &
    ClumpFactorInit_pft        =>  plt_morph%ClumpFactorInit_pft   , &
    ClumpFactor                =>  plt_morph%ClumpFactor   , &
    NumRootAxes_pft            =>  plt_morph%NumRootAxes_pft    &
  )
  iYearPlanting_pft(NZ)=IYRX(NZ)
  iDayPlanting_pft(NZ)=IDAYX(NZ)
  iYearPlantHarvest_pft(NZ)=IYRY(NZ)
  iDayPlantHarvest_pft(NZ)=IDAYY(NZ)
  PPI(NZ)=PPZ(NZ)
  PPX(NZ)=PPI(NZ)
  ClumpFactor(NZ)=ClumpFactorInit_pft(NZ)

  MaxCanPStomaResistH2O_pft(NZ)=RSMX(NZ)/3600.0_r8
  CO2CuticleResist_pft(NZ)=RSMX(NZ)*1.56_r8
  rCNNonstructRemob_pft(NZ)=2.5_r8
  rCPNonstructRemob_pft(NZ)=25.0_r8
  RootFracRemobilizableBiom(NZ)=AMIN1(RootrNC_pft(NZ)*rCNNonstructRemob_pft(NZ),&
    RootrPC_pft(NZ)*rCPNonstructRemob_pft(NZ))
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
    inonstruct                => pltpar%inonstruct, &
    ifoliar                   => pltpar%ifoliar , &
    inonfoliar                => pltpar%inonfoliar, &
    istalk                    => pltpar%istalk  , &
    iroot                     => pltpar%iroot   , &
    icwood                    => pltpar%icwood  , &
    iprotein                  => pltpar%iprotein  ,&
    icarbhyro                 => pltpar%icarbhyro ,&
    icellulos                 => pltpar%icellulos ,&
    ilignin                   =>  pltpar%ilignin   , &
    NumLitterGroups           => pltpar%NumLitterGroups , &
    RefLeafAppearRate_pft     =>  plt_pheno%RefLeafAppearRate_pft      , &
    iPlantTurnoverPattern_pft =>  plt_pheno%iPlantTurnoverPattern_pft    , &
    iPlantRootProfile_pft  =>  plt_pheno%iPlantRootProfile_pft    , &
    MatureGroup_pft           =>  plt_pheno%MatureGroup_pft   , &
    CFOPE                     =>  plt_soilchem%CFOPE  , &
    FNOD                      =>  plt_allom%FNOD      , &
    iPlantNfixType            =>  plt_morph%iPlantNfixType    , &
    NumCogrowNode   =>  plt_morph%NumCogrowNode       &
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
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
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
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
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
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.1 &
    .OR.iPlantTurnoverPattern_pft(NZ).GE.3)THEN
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
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
    CFOPE(ielmc,istalk,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,istalk,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,istalk,icellulos,NZ)=0.30_r8
    CFOPE(ielmc,istalk,ilignin,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
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
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
    CFOPE(ielmc,iroot,iprotein,NZ)=0.07_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ)=0.25_r8
    CFOPE(ielmc,iroot,icellulos,NZ)=0.30_r8
    CFOPE(ielmc,iroot,ilignin,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
    CFOPE(ielmc,iroot,iprotein,NZ)=0.057_r8
    CFOPE(ielmc,iroot,icarbhyro,NZ)=0.263_r8
    CFOPE(ielmc,iroot,icellulos,NZ)=0.542_r8
    CFOPE(ielmc,iroot,ilignin,NZ)=0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.1&
    .OR.iPlantTurnoverPattern_pft(NZ).GE.3)THEN
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
!     NumCogrowNode=number of concurrently growing nodes
!
  IF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
    FNOD(NZ)=1.0_r8
    IF(MatureGroup_pft(NZ).LE.10)THEN
      NumCogrowNode(NZ)=3
    ELSEIF(MatureGroup_pft(NZ).LE.15)THEN
      NumCogrowNode(NZ)=4
    ELSE
      NumCogrowNode(NZ)=5
    ENDIF
  ELSE
    FNOD(NZ)=AMAX1(1.0_r8,0.04_r8/RefLeafAppearRate_pft(NZ))
    NumCogrowNode(NZ)=24
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
    DATAP                      =>  plt_site%DATAP   , &
    iPlantPhotosynthesisType   =>  plt_photo%iPlantPhotosynthesisType , &
    HighTCLimtSeed_pft         =>  plt_pheno%HighTCLimtSeed_pft   , &
    TCelcius4LeafOffHarden_pft =>  plt_pheno%TCelcius4LeafOffHarden_pft    , &
    TCelsChill4Leaf_pft        =>  plt_pheno%TCelsChill4Leaf_pft   , &
    OFFST                      =>  plt_pheno%OFFST  , &
    PlantInitThermoAdaptZone   =>  plt_pheno%PlantInitThermoAdaptZone , &
    iPlantThermoAdaptZone      =>  plt_pheno%iPlantThermoAdaptZone  , &
    SSTX                       =>  plt_pheno%SSTX     &
  )
!
!     PFT THERMAL ACCLIMATION
!
!     ZTYP,PlantInitThermoAdaptZone=dynamic,initial thermal adaptation zone from PFT file
!     OFFST=shift in Arrhenius curve for thermal adaptation (oC)
!     TCZ,TCelcius4LeafOffHarden_pft=threshold temperature for leafout,leafoff
!     HTC=high temperature threshold for grain number loss (oC)
!     SSTX=sensitivity to HTC (seeds oC-1 above HTC)
!
  iPlantThermoAdaptZone(NZ)=PlantInitThermoAdaptZone(NZ)
  OFFST(NZ)=2.667_r8*(2.5_r8-iPlantThermoAdaptZone(NZ))
  TCelsChill4Leaf_pft(NZ)=TCZD-OFFST(NZ)
  TCelcius4LeafOffHarden_pft(NZ)=AMIN1(15.0_r8,TCXD-OFFST(NZ))
  IF(iPlantPhotosynthesisType(NZ).EQ.ic3_photo)THEN
    IF(DATAP(NZ)(1:4).EQ.'soyb')THEN
      HighTCLimtSeed_pft(NZ)=30.0_r8+3.0_r8*iPlantThermoAdaptZone(NZ)
      SSTX(NZ)=0.002_r8
    ELSE
      HighTCLimtSeed_pft(NZ)=27.0_r8+3.0_r8*iPlantThermoAdaptZone(NZ)
      SSTX(NZ)=0.002_r8
    ENDIF
  ELSE
    HighTCLimtSeed_pft(NZ)=27.0_r8+3.0_r8*iPlantThermoAdaptZone(NZ)
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
    CNRTS                     =>  plt_allom%CNRTS   , &
    CPRTS                     =>  plt_allom%CPRTS   , &
    RootBiomGrowthYield       =>  plt_allom%RootBiomGrowthYield    , &
    RootrNC_pft               =>  plt_allom%RootrNC_pft   , &
    RootrPC_pft               =>  plt_allom%RootrPC_pft   , &
    CMinPO4Root_pft           =>  plt_rbgc%CMinPO4Root_pft   , &
    KmPO4Root_pft             =>  plt_rbgc%KmPO4Root_pft   , &
    VmaxPO4Root_pft           =>  plt_rbgc%VmaxPO4Root_pft   , &
    CminNO3Root_pft           =>  plt_rbgc%CminNO3Root_pft   , &
    KmNO3Root_pft             =>  plt_rbgc%KmNO3Root_pft   , &
    VmaxNO3Root_pft           =>  plt_rbgc%VmaxNO3Root_pft   , &
    CMinNH4Root_pft           =>  plt_rbgc%CMinNH4Root_pft   , &
    VmaxNH4Root_pft           =>  plt_rbgc%VmaxNH4Root_pft   , &
    KmNH4Root_pft             =>  plt_rbgc%KmNH4Root_pft   , &
    CumSoilThickness          =>  plt_site%CumSoilThickness   , &
    NU                        =>  plt_site%NU       , &
    NL                        =>  plt_site%NL       , &
    NGTopRootLayer_pft        =>  plt_morph%NGTopRootLayer_pft      , &
    SeedDepth_pft             =>  plt_morph%SeedDepth_pft   , &
    SeedAreaMean_pft          =>  plt_morph%SeedAreaMean_pft    , &
    SeedCMass                 =>  plt_morph%SeedCMass    , &
    PlantinDepth              =>  plt_morph%PlantinDepth  , &
    Max2ndRootRadius_pft1         =>  plt_morph%Max2ndRootRadius_pft1  , &
    Root2ndXSecArea_pft         =>  plt_morph%Root2ndXSecArea_pft  , &
    Max1stRootRadius_pft1         =>  plt_morph%Max1stRootRadius_pft1  , &
    Root1stXSecArea_pft          =>  plt_morph%Root1stXSecArea_pft  , &
    Max2ndRootRadius_pft          =>  plt_morph%Max2ndRootRadius_pft  , &
    Max1stRootRadius_pft          =>  plt_morph%Max1stRootRadius_pft  , &
    Root2ndSpecLen_pft          =>  plt_morph%Root2ndSpecLen_pft  , &
    NIXBotRootLayer_pft       =>  plt_morph%NIXBotRootLayer_pft    , &
    RSRR                      =>  plt_morph%RSRR    , &
    Root1stSpecLen_pft           =>  plt_morph%Root1stSpecLen_pft  , &
    RootPorosity_pft              =>  plt_morph%RootPorosity_pft   , &
    RootPoreTortu4Gas         =>  plt_morph%RootPoreTortu4Gas   , &
    RootRaidus_rpft           =>  plt_morph%RootRaidus_rpft   , &
    RootVolPerMassC_pft       =>  plt_morph%RootVolPerMassC_pft    , &
    RSRA                      =>  plt_morph%RSRA    , &
    NIXBotRootLayer_rpft      =>  plt_morph%NIXBotRootLayer_rpft    , &
    SeedVolumeMean_pft        =>  plt_morph%SeedVolumeMean_pft    , &
    SeedMeanLen_pft        =>  plt_morph%SeedMeanLen_pft      &
  )
!
!     SEED CHARACTERISTICS
!
!     SeedVolumeMean_pft,SeedMeanLen_pft,SeedAreaMean_pft=seed volume(m3),length(m),AREA3(NU)(m2)
!     SeedCMass=seed C mass (g) from PFT file
!
  call calc_seed_geometry(SeedCMass(NZ),SeedVolumeMean_pft(NZ),SeedMeanLen_pft(NZ),SeedAreaMean_pft(NZ))
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) DIMENSIONS, UPTAKE PARAMETERS
!
!     SeedDepth_pft=seeding depth(m) from PFT management file
!     CumSoilThickness=depth to soil layer bottom from surface(m)
!     NG,NIX,NIXBotRootLayer_rpft=seeding,upper,lower rooting layer
!     CNRTS,CPRTS=N,P root growth yield
!     Max1stRootRadius_pft,Max2ndRootRadius_pft=maximum primary,secondary mycorrhizal radius (m)
!     PORT=mycorrhizal porosity
!     VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake(g m-2 h-1),Km(uM), min concn (uM)
!     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     RSRR,RSRA=radial,axial root resistivity (m2 MPa-1 h-1)
!
  SeedDepth_pft(NZ)=PlantinDepth(NZ)
  D9795: DO L=NU,NL
    IF(SeedDepth_pft(NZ).GE.CumSoilThickness(L-1) &
      .AND.SeedDepth_pft(NZ).LT.CumSoilThickness(L))THEN
      NGTopRootLayer_pft(NZ)=L
      NIXBotRootLayer_pft(NZ)=L
      D9790: DO NR=1,pltpar%MaxNumRootAxes
        NIXBotRootLayer_rpft(NR,NZ)=L
      ENDDO D9790
    ENDIF
  ENDDO D9795
  CNRTS(NZ)=RootrNC_pft(NZ)*RootBiomGrowthYield(NZ)
  CPRTS(NZ)=RootrPC_pft(NZ)*RootBiomGrowthYield(NZ)
  Max1stRootRadius_pft(2,NZ)=5.0E-06_r8
  Max2ndRootRadius_pft(2,NZ)=5.0E-06_r8
  RootPorosity_pft(2,NZ)=RootPorosity_pft(1,NZ)
  VmaxNH4Root_pft(2,NZ)=VmaxNH4Root_pft(1,NZ)
  KmNH4Root_pft(2,NZ)=KmNH4Root_pft(1,NZ)
  CMinNH4Root_pft(2,NZ)=CMinNH4Root_pft(1,NZ)
  VmaxNO3Root_pft(2,NZ)=VmaxNO3Root_pft(1,NZ)
  KmNO3Root_pft(2,NZ)=KmNO3Root_pft(1,NZ)
  CminNO3Root_pft(2,NZ)=CminNO3Root_pft(1,NZ)
  VmaxPO4Root_pft(2,NZ)=VmaxPO4Root_pft(1,NZ)
  KmPO4Root_pft(2,NZ)=KmPO4Root_pft(1,NZ)
  CMinPO4Root_pft(2,NZ)=CMinPO4Root_pft(1,NZ)
  RSRR(2,NZ)=1.0E+04_r8
  RSRA(2,NZ)=1.0E+12_r8
!
!     RootPoreTortu4Gas=tortuosity for gas transport
!     RootRaidus_rpft=path length for radial diffusion within root (m)
!     RootVolPerMassC_pft=volume:C ratio (m3 g-1)
!     Root1stSpecLen_pft,Root2ndSpecLen_pft=specific primary,secondary root length (m g-1)
!     Root1stXSecArea_pft,Root2ndXSecArea_pft=specific primary,secondary root area (m2 g-1)
!
  D500: DO N=1,2
    RootPoreTortu4Gas(N,NZ)=RootPorosity_pft(N,NZ)**1.33_r8
    RootRaidus_rpft(N,NZ)=LOG(1.0_r8/SQRT(AMAX1(0.01_r8,RootPorosity_pft(N,NZ))))
    RootVolPerMassC_pft(N,NZ)=ppmc/(0.05_r8*(1.0_r8-RootPorosity_pft(N,NZ)))
    Root1stSpecLen_pft(N,NZ)=RootVolPerMassC_pft(N,NZ)/(PICON*Max1stRootRadius_pft(N,NZ)**2._r8)
    Root2ndSpecLen_pft(N,NZ)=RootVolPerMassC_pft(N,NZ)/(PICON*Max2ndRootRadius_pft(N,NZ)**2._r8)
    Max1stRootRadius_pft1(N,NZ)=Max1stRootRadius_pft(N,NZ)
!    2*SQRT(0.25*(1.0-RootPorosity_pft(N,NZ)))
    Max2ndRootRadius_pft1(N,NZ)=Max2ndRootRadius_pft(N,NZ)
!    2*SQRT(0.25*(1.0-RootPorosity_pft(N,NZ)))
    Root1stXSecArea_pft(N,NZ)=PICON*Max1stRootRadius_pft1(N,NZ)**2
    Root2ndXSecArea_pft(N,NZ)=PICON*Max2ndRootRadius_pft1(N,NZ)**2
  ENDDO D500
  end associate
  end subroutine InitDimensionsandUptake
!------------------------------------------------------------------------------------------

  subroutine InitPlantPhenoMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: K,L,M,N,NB
  associate(                           &
    NU                                =>  plt_site%NU      , &
    PPX                               =>  plt_site%PPX     , &
    PlantPopulation_pft               =>  plt_site%PlantPopulation_pft      , &
    ALAT                              =>  plt_site%ALAT    , &
    AREA3                             =>  plt_site%AREA3   , &
    iPlantRootProfile_pft          =>  plt_pheno%iPlantRootProfile_pft , &
    LeafNumberAtFloralInit_brch       =>  plt_pheno%LeafNumberAtFloralInit_brch , &
    MatureGroup_brch                  =>  plt_pheno%MatureGroup_brch , &
    KHiestGroLeafNode_brch                   =>  plt_pheno%KHiestGroLeafNode_brch , &
    KLowestGroLeafNode_brch         =>  plt_pheno%KLowestGroLeafNode_brch , &
    NodeNumNormByMatgrp_brch          =>  plt_pheno%NodeNumNormByMatgrp_brch  , &
    ReprodNodeNumNormByMatrgrp_brch   =>  plt_pheno%ReprodNodeNumNormByMatrgrp_brch  , &
    HourFailGrainFill_brch            =>  plt_pheno%HourFailGrainFill_brch   , &
    HoursDoingRemob_brch              =>  plt_pheno%HoursDoingRemob_brch  , &
    Hours4LenthenPhotoPeriod_brch     =>  plt_pheno%Hours4LenthenPhotoPeriod_brch   , &
    Hours4LeafOff_brch                =>  plt_pheno%Hours4LeafOff_brch   , &
    iPlantBranchState_brch            =>  plt_pheno%iPlantBranchState_brch  , &
    Hours4ShortenPhotoPeriod_brch     =>  plt_pheno%Hours4ShortenPhotoPeriod_brch   , &
    iPlantCalendar_brch               =>  plt_pheno%iPlantCalendar_brch , &
    Hours4Leafout_brch                =>  plt_pheno%Hours4Leafout_brch   , &
    HourCount2LeafOut_brch          =>  plt_pheno%HourCount2LeafOut_brch   , &
    HoursCanopyPSITooLow              =>  plt_pheno%HoursCanopyPSITooLow   , &
    MatureGroup_pft                   =>  plt_pheno%MatureGroup_pft, &
    PetioleChemElmRemobFlx_brch       =>  plt_pheno%PetioleChemElmRemobFlx_brch  , &
    TotReproNodeNumNormByMatrgrp_brch =>  plt_pheno%TotReproNodeNumNormByMatrgrp_brch , &
    TotalNodeNumNormByMatgrp_brch     =>  plt_pheno%TotalNodeNumNormByMatgrp_brch , &
    C4PhotosynDowreg_brch             =>  plt_photo%C4PhotosynDowreg_brch  , &
    CPOOL3                            =>  plt_photo%CPOOL3 , &
    CPOOL4                            =>  plt_photo%CPOOL4 , &
    CHILL                             =>  plt_photo%CHILL  , &
    RubiscoActivity_brch              =>  plt_photo%RubiscoActivity_brch   , &
    CMassHCO3BundleSheath_node        =>  plt_photo%CMassHCO3BundleSheath_node   , &
    CMassCO2BundleSheath_node         =>  plt_photo%CMassCO2BundleSheath_node   , &
    NumOfLeaves_brch                  =>  plt_morph%NumOfLeaves_brch  , &
    CanopyHeight_pft                  =>  plt_morph%CanopyHeight_pft    , &
    KLeafNumber_brch                  =>  plt_morph%KLeafNumber_brch , &
    XTLI                              =>  plt_morph%XTLI   , &
    BranchNumber_pft                  =>  plt_morph%BranchNumber_pft    , &
    ShootNodeNumber_brch              =>  plt_morph%ShootNodeNumber_brch  , &
    CanopyLeafArea_pft                =>  plt_morph%CanopyLeafArea_pft  , &
    CanopyStemALyr_brch          =>  plt_morph%CanopyStemALyr_brch  , &
    StemAreaZsec_brch                 =>  plt_morph%StemAreaZsec_brch  , &
    CanopyStemA_pft                   =>  plt_morph%CanopyStemA_pft  , &
    NodeNumberAtAnthesis_brch         =>  plt_morph%NodeNumberAtAnthesis_brch  , &
    PotentialSeedSites_brch           =>  plt_morph%PotentialSeedSites_brch  , &
    InternodeHeightLive_brch          =>  plt_morph%InternodeHeightLive_brch , &
    SeedNumSet_brch                =>  plt_morph%SeedNumSet_brch  , &
    LeafAreaLive_brch                 =>  plt_morph%LeafAreaLive_brch  , &
    LeafAreaDying_brch                =>  plt_morph%LeafAreaDying_brch  , &
    CanopyLeafAreaByLayer_pft         =>  plt_morph%CanopyLeafAreaByLayer_pft  , &
    CanopyLeafALyr_pft                =>  plt_morph%CanopyLeafALyr_pft  , &
    CanopyStemApft_lyr                =>  plt_morph%CanopyStemApft_lyr  , &
    LeafAreaZsec_brch                 =>  plt_morph%LeafAreaZsec_brch   , &
    LeafAreaNode_brch                 =>  plt_morph%LeafAreaNode_brch  , &
    CanPBranchHeight                  =>  plt_morph%CanPBranchHeight , &
    HypoctoHeight_pft                 =>  plt_morph%HypoctoHeight_pft  , &
    BranchNumber_brch                 =>  plt_morph%BranchNumber_brch   , &
    NodeNum2InitFloral_brch       =>  plt_morph%NodeNum2InitFloral_brch  , &
    KLEAFX                            =>  plt_morph%KLEAFX , &
    NumOfBranches_pft                 =>  plt_morph%NumOfBranches_pft      &
  )
!
!     INITIALIZE PLANT PHENOLOGY
!
!     PP=population (grid cell-1)
!
  PlantPopulation_pft(NZ)=PPX(NZ)*AREA3(NU)
  plt_pheno%doInitPlant_pft(NZ)=ifalse
  plt_pheno%iPlantShootState_pft(NZ)=iDead
  plt_pheno%iPlantRootState_pft(NZ)=iDead
  BranchNumber_pft(NZ)=0
  NumOfBranches_pft(NZ)=0
  HypoctoHeight_pft(NZ)=0._r8
  CanopyHeight_pft(NZ)=0._r8
  D10: DO NB=1,MaxNumBranches
    plt_pheno%doInitLeafOut_brch(NB,NZ)=iEnable
    plt_pheno%doPlantLeafOut_brch(NB,NZ)=iDisable
    plt_pheno%doPlantLeaveOff_brch(NB,NZ)=iDisable
    plt_pheno%Prep4Literfall_brch(NB,NZ)=ifalse
    plt_pheno%Hours4LiterfalAftMature_brch(NB,NZ)=0
    MatureGroup_brch(NB,NZ)=MatureGroup_pft(NZ)
    ShootNodeNumber_brch(NB,NZ)=XTLI(NZ)
    NodeNum2InitFloral_brch(NB,NZ)=ShootNodeNumber_brch(NB,NZ)
    NodeNumberAtAnthesis_brch(NB,NZ)=0._r8
    NumOfLeaves_brch(NB,NZ)=0._r8
    LeafNumberAtFloralInit_brch(NB,NZ)=0._r8
    KLeafNumber_brch(NB,NZ)=1
    KLEAFX(NB,NZ)=1
    KHiestGroLeafNode_brch(NB,NZ)=1
    KLowestGroLeafNode_brch(NB,NZ)=0
    NodeNumNormByMatgrp_brch(NB,NZ)=0._r8
    ReprodNodeNumNormByMatrgrp_brch(NB,NZ)=0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ)=0._r8
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ)=0._r8
    Hours4LenthenPhotoPeriod_brch(NB,NZ)=0._r8
    Hours4ShortenPhotoPeriod_brch(NB,NZ)=0._r8
    Hours4Leafout_brch(NB,NZ)=Hours4LenthenPhotoPeriod_brch(NB,NZ)
    Hours4LeafOff_brch(NB,NZ)=Hours4ShortenPhotoPeriod_brch(NB,NZ)
    HourCount2LeafOut_brch(NB,NZ)=0._r8
    RubiscoActivity_brch(NB,NZ)=1.0
    C4PhotosynDowreg_brch(NB,NZ)=1.0
    HourFailGrainFill_brch(NB,NZ)=0
    HoursDoingRemob_brch(NB,NZ)=0
    BranchNumber_brch(NB,NZ)=0
    plt_pheno%iPlantBranchState_brch(NB,NZ)=iDead
    D15: DO M=1,pltpar%NumGrowthStages
      iPlantCalendar_brch(M,NB,NZ)=0
    ENDDO D15
  ENDDO D10
!
!     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
!
  HoursCanopyPSITooLow(NZ)=0._r8
  CHILL(NZ)=0._r8
  plt_biom%NonstructElm_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%NodulNonstElm_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%PetoleChemElm_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%ShootChemElm_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%StalkChemElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%LeafChemElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%StalkRsrvElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%HuskChemElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%GrainChemElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%EarChemElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%CanopyNodulChemElm_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_pheno%LeafElmntRemobFlx_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%LeafChemElmRemob_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_pheno%PetioleChemElmRemobFlx_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%PetioleChemElmRemob_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8
  plt_biom%SenecStalkChemElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)=0._r8

  D25: DO NB=1,MaxNumBranches
    plt_biom%StalkBiomassC_brch(NB,NZ)=0._r8
    plt_biom%LeafPetolBiomassC_brch(NB,NZ)=0._r8
    PotentialSeedSites_brch(NB,NZ)=0._r8
    SeedNumSet_brch(NB,NZ)=0._r8
    plt_allom%GrainSeedBiomCMean_brch(NB,NZ)=0._r8
    LeafAreaLive_brch(NB,NZ)=0._r8
    plt_rbgc%NH3Dep2_brch(NB,NZ)=0._r8
    LeafAreaDying_brch(NB,NZ)=0._r8
    CanPBranchHeight(NB,NZ)=0._r8
    D5: DO L=1,NumOfCanopyLayers1
      CanopyStemALyr_brch(L,NB,NZ)=0._r8
      DO N=1,NumOfLeafZenithSectors1
        StemAreaZsec_brch(N,L,NB,NZ)=0._r8
      enddo
    ENDDO D5
    DO K=0,MaxNodesPerBranch1
      LeafAreaNode_brch(K,NB,NZ)=0._r8
      InternodeHeightLive_brch(K,NB,NZ)=0._r8
      plt_morph%InternodeHeightDying_brch(K,NB,NZ)=0._r8
      plt_morph%PetioleLengthNode_brch(K,NB,NZ)=0._r8
      plt_biom%LeafProteinCNode_brch(K,NB,NZ)=0._r8
      plt_biom%PetioleProteinCNode_brch(K,NB,NZ)=0._r8
      plt_biom%LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)=0._r8
      plt_biom%InternodeChemElm_brch(1:NumPlantChemElms,K,NB,NZ)=0._r8
      plt_biom%PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)=0._r8


      D55: DO L=1,NumOfCanopyLayers1
        CanopyLeafAreaByLayer_pft(L,K,NB,NZ)=0._r8
        plt_biom%LeafChemElmByLayerNode_brch(1:NumPlantChemElms,L,K,NB,NZ)=0._r8
      ENDDO D55
      IF(K.NE.0)THEN
        CPOOL3(K,NB,NZ)=0._r8
        CMassCO2BundleSheath_node(K,NB,NZ)=0._r8
        CMassHCO3BundleSheath_node(K,NB,NZ)=0._r8
        CPOOL4(K,NB,NZ)=0._r8
        D45: DO L=1,NumOfCanopyLayers1
          DO N=1,NumOfLeafZenithSectors1
            LeafAreaZsec_brch(N,L,K,NB,NZ)=0._r8
          enddo
        ENDDO D45
      ENDIF
    enddo
  ENDDO D25
  D35: DO L=1,NumOfCanopyLayers1
    CanopyLeafALyr_pft(L,NZ)=0._r8
    plt_biom%CanopyLeafCLyr_pft(L,NZ)=0._r8
    CanopyStemApft_lyr(L,NZ)=0._r8
  ENDDO D35
  plt_biom%CanopyNonstructElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%CanopyNonstructElmConc_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%NoduleNonstructCconc_pft(NZ)=0._r8
  plt_biom%ShootChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%LeafChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%PetioleChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%StalkChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%CanopyStalkC_pft(NZ)=0._r8
  plt_biom%StalkRsrvElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%HuskChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%EarChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%GrainChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%RootStructElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%NoduleChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
  plt_biom%CanopyLeafShethC_pft(NZ)=0._r8
  CanopyLeafArea_pft(NZ)=0._r8
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
    NU                            => plt_site%NU        , &
    AREA3                         => plt_site%AREA3     , &
    CFOPE                         => plt_soilchem%CFOPE , &
    ETCanopy_pft                  => plt_ew%ETCanopy_pft       , &
    StandDeadKCompElms_pft => plt_biom%StandDeadKCompElms_pft    , &
    StandingDeadChemElms_pft      => plt_biom%StandingDeadChemElms_pft    , &
    StandingDeadInitC_pft         => plt_biom%StandingDeadInitC_pft    , &
    RootBiomCPerPlant_pft         => plt_biom%RootBiomCPerPlant_pft     , &
    rNCStalk_pft                  => plt_allom%rNCStalk_pft   , &
    rPCStalk_pft                  => plt_allom%rPCStalk_pft    , &
    PlantExudChemElmCum_pft       => plt_rbgc%PlantExudChemElmCum_pft    , &
    NH3EmiCum_pft                 => plt_bgcr%NH3EmiCum_pft     , &
    NH3Dep2Can_pft                => plt_bgcr%NH3Dep2Can_pft     , &
    SurfLitrfallChemElms_pft      => plt_bgcr%SurfLitrfallChemElms_pft     , &
    GrossCO2Fix_pft               => plt_bgcr%GrossCO2Fix_pft     , &
    CanopyPlusNoduRespC_pft       => plt_bgcr%CanopyPlusNoduRespC_pft     , &
    GrossResp_pft                 => plt_bgcr%GrossResp_pft     , &
    PlantN2FixCum_pft             => plt_bgcr%PlantN2FixCum_pft    , &
    LitrfallChemElms_pft          => plt_bgcr%LitrfallChemElms_pft     , &
    CanopyStemA_pft               => plt_morph%CanopyStemA_pft    , &
    icwood                        => pltpar%icwood      , &
    NetCumElmntFlx2Plant_pft      => plt_pheno%NetCumElmntFlx2Plant_pft      &

  )
!
!     INITIALIZE MASS BALANCE CHECKS
!
  IF(.not.is_restart().AND.is_first_year)THEN
    GrossCO2Fix_pft(NZ)=0._r8
    SurfLitrfallChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
    GrossResp_pft(NZ)=0._r8
    CanopyPlusNoduRespC_pft(NZ)=0._r8
    PlantExudChemElmCum_pft(1:NumPlantChemElms,NZ)=0._r8
    LitrfallChemElms_pft(1:NumPlantChemElms,NZ)=0._r8

    PlantN2FixCum_pft(NZ)=0._r8
    NH3Dep2Can_pft(NZ)=0._r8
    NH3EmiCum_pft(NZ)=0._r8
    plt_distb%CO2ByFire_pft(NZ)=0._r8
    plt_distb%CH4ByFire_pft(NZ)=0._r8
    plt_distb%O2ByFire_pft(NZ)=0._r8
    plt_distb%NH3byFire_pft(NZ)=0._r8
    plt_distb%N2ObyFire_pft(NZ)=0._r8
    plt_distb%PO4byFire_pft(NZ)=0._r8
    plt_distb%EcoHavstElmntCum_pft(1:NumPlantChemElms,NZ)=0._r8
    plt_distb%EcoHavstElmnt_pft(1:NumPlantChemElms,NZ)=0._r8
    NetCumElmntFlx2Plant_pft(1:NumPlantChemElms,NZ)=0._r8
    ETCanopy_pft(NZ)=0._r8
    StandingDeadChemElms_pft(1:NumPlantChemElms,NZ)=0._r8
    WTSTDX=StandingDeadInitC_pft(NZ)*AREA3(NU)
    D155: DO M=1,jsken
      StandDeadKCompElms_pft(ielmc,M,NZ)=WTSTDX*CFOPE(ielmc,icwood,M,NZ)
      StandDeadKCompElms_pft(ielmn,M,NZ)=WTSTDX*rNCStalk_pft(NZ)*CFOPE(ielmn,icwood,M,NZ)
      StandDeadKCompElms_pft(ielmp,M,NZ)=WTSTDX*rPCStalk_pft(NZ)*CFOPE(ielmp,icwood,M,NZ)
    ENDDO D155
    DO NE=1,NumPlantChemElms
      StandingDeadChemElms_pft(NE,NZ)=StandingDeadChemElms_pft(NE,NZ)+&
        sum(StandDeadKCompElms_pft(NE,1:jsken,NZ))
    ENDDO
  ENDIF
  end associate
  end subroutine InitMassBalance
!------------------------------------------------------------------------------------------

  subroutine InitPlantHeatandWater(NZ)

  implicit none
  integer, intent(in) :: NZ
  associate(                          &
    ATCA                     =>  plt_site%ATCA  , &
    OSMO                     =>  plt_ew%OSMO    , &
    TKC                      =>  plt_ew%TKC     , &
    Transpiration_pft        =>  plt_ew%Transpiration_pft    , &
    VHeatCapCanP             =>  plt_ew%VHeatCapCanP   , &
    PSICanopy_pft            =>  plt_ew%PSICanopy_pft  , &
    PSICanopyTurg_pft        =>  plt_ew%PSICanopyTurg_pft   , &
    PSICanopyOsmo_pft        =>  plt_ew%PSICanopyOsmo_pft   , &
    ENGYX                    =>  plt_ew%ENGYX   , &
    DTKC                     =>  plt_ew%DTKC    , &
    TCelciusCanopy_pft       =>  plt_ew%TCelciusCanopy_pft    , &
    TKG                      =>  plt_pheno%TKG  , &
    TCG                      =>  plt_pheno%TCG  , &
    fTgrowCanP               =>  plt_pheno%fTgrowCanP , &
    ShootChemElms_pft        =>  plt_biom%ShootChemElms_pft , &
    FracRadPARbyCanopy_pft   =>  plt_rad%FracRadPARbyCanopy_pft    &
  )
!
!     INITIALIZE PLANT HEAT AND WATER STATUS
!
!     VHeatCapCanP=canopy heat capacity (MJ m-3 K-1)
!     TCelciusCanopy_pft,TKC=canopy temperature for growth (oC,K)
!     TCG,TKG=canopy temperature for phenology (oC,K)
!     PSICanopy_pft,PSICanopyOsmo_pft,PSICanopyTurg_pft=canopy total,osmotic,turgor water potl(MPa)
!
  VHeatCapCanP(NZ)=cpw*ShootChemElms_pft(ielmc,NZ)*10.0E-06
  ENGYX(NZ)=0._r8
  DTKC(NZ)=0._r8
  TCelciusCanopy_pft(NZ)=ATCA
  TKC(NZ)=units%Celcius2Kelvin(TCelciusCanopy_pft(NZ))
  TCG(NZ)=TCelciusCanopy_pft(NZ)
  TKG(NZ)=units%Celcius2Kelvin(TCG(NZ))
  fTgrowCanP(NZ)=1.0
  PSICanopy_pft(NZ)=-1.0E-03
  PSICanopyOsmo_pft(NZ)=OSMO(NZ)+PSICanopy_pft(NZ)
  PSICanopyTurg_pft(NZ)=AZMAX1(PSICanopy_pft(NZ)-PSICanopyOsmo_pft(NZ))
  Transpiration_pft(NZ)=0._r8
  FracRadPARbyCanopy_pft(NZ)=0._r8
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
    OSMO                         =>  plt_ew%OSMO       , &
    PSICanopy_pft                =>  plt_ew%PSICanopy_pft     , &
    PSIRoot_pvr                   =>  plt_ew%PSIRoot_pvr     , &
    PSIRootTurg_vr               =>  plt_ew%PSIRootTurg_vr      , &
    PSIRootOSMO_vr               =>  plt_ew%PSIRootOSMO_vr      , &
    RootN2Fix_pvr                =>  plt_bgcr%RootN2Fix_pvr    , &
    trcg_rootml_pvr               =>  plt_rbgc%trcg_rootml_pvr     , &
    trcs_rootml_pvr               =>  plt_rbgc%trcs_rootml_pvr     , &
    OXYE                         =>  plt_site%OXYE     , &
    COXYE                        =>  plt_site%COXYE    , &
    CO2EI                        =>  plt_site%CO2EI    , &
    NL                           =>  plt_site%NL       , &
    ATCA                         =>  plt_site%ATCA     , &
    CCO2EI                       =>  plt_site%CCO2EI   , &
    RootFracRemobilizableBiom    =>  plt_allom%RootFracRemobilizableBiom   , &
    RootProteinConc_pvr          =>  plt_biom%RootProteinConc_pvr   , &
    RootProteinC_pvr             =>  plt_biom%RootProteinC_pvr    , &
    SeedDepth_pft                =>  plt_morph%SeedDepth_pft   , &
    RootVH2O_pvr                  =>  plt_morph%RootVH2O_pvr  , &
    RootPoreVol_pvr                =>  plt_morph%RootPoreVol_pvr   , &
    Radius2ndRoot_pvr          =>  plt_morph%Radius2ndRoot_pvr   , &
    Root1stRadius_pvr           =>  plt_morph%Root1stRadius_pvr   , &
    Max1stRootRadius_pft             =>  plt_morph%Max1stRootRadius_pft  , &
    Max2ndRootRadius_pft             =>  plt_morph%Max2ndRootRadius_pft  , &
    NumRootAxes_pft              =>  plt_morph%NumRootAxes_pft      &
  )
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
!
!     PSIRoot_pvr,PSIRootOSMO_vr,PSIRootTurg_vr=root,myco total,osmotic,turgor water potl(MPa)
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
      PSIRoot_pvr(N,L,NZ)=-0.01_r8
      PSIRootOSMO_vr(N,L,NZ)=OSMO(NZ)+PSIRoot_pvr(N,L,NZ)
      PSIRootTurg_vr(N,L,NZ)=AZMAX1(PSIRoot_pvr(N,L,NZ)-PSIRootOSMO_vr(N,L,NZ))
      plt_biom%RootMycoNonstElm_pvr(1:NumPlantChemElms,N,L,NZ)=0._r8
      plt_biom%RootNonstructElmConc_pvr(1:NumPlantChemElms,N,L,NZ)=0._r8
      RootProteinConc_pvr(N,L,NZ)=RootFracRemobilizableBiom(NZ)
      plt_biom%RootMycoActiveBiomC_pvr(N,L,NZ)=0._r8
      plt_biom% PopuPlantRootC_vr(N,L,NZ)=0._r8
      RootProteinC_pvr(N,L,NZ)=0._r8
      plt_morph%Root1stXNumL_pvr(N,L,NZ)=0._r8
      plt_morph%Root2ndXNum_pvr(N,L,NZ)=0._r8
      plt_morph%RootLenPerPlant_pvr(N,L,NZ)=0._r8
      plt_morph%RootLenDensPerPlant_pvr(N,L,NZ)=0._r8
      RootPoreVol_pvr(N,L,NZ)=0._r8
      RootVH2O_pvr(N,L,NZ)=0._r8
      Root1stRadius_pvr(N,L,NZ)=Max1stRootRadius_pft(N,NZ)
      Radius2ndRoot_pvr(N,L,NZ)=Max2ndRootRadius_pft(N,NZ)
      plt_morph%RootAreaPerPlant_pvr(N,L,NZ)=0._r8
      plt_morph%AveLen2ndRoot_pvr(N,L,NZ)=1.0E-03
      plt_rbgc%RootNutUptake_pvr(ids_NH4,N,L,NZ)=0._r8
      plt_rbgc%RootNutUptake_pvr(ids_NO3,N,L,NZ)=0._r8
      plt_rbgc%RootNutUptake_pvr(ids_H2PO4,N,L,NZ)=0._r8
      plt_rbgc%RootNutUptake_pvr(ids_H1PO4,N,L,NZ)=0._r8
      plt_rbgc%RootNutUptake_pvr(ids_NH4B,N,L,NZ)=0._r8
      plt_rbgc%RootNutUptake_pvr(ids_NO3B,N,L,NZ)=0._r8
      plt_rbgc%RootNutUptake_pvr(ids_H2PO4B,N,L,NZ)=0._r8
      plt_rbgc%RootNutUptake_pvr(ids_H1PO4B,N,L,NZ)=0._r8
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
      trcg_rootml_pvr(idg_CO2,N,L,NZ)=CCO2A*RootPoreVol_pvr(N,L,NZ)
      trcs_rootml_pvr(idg_CO2,N,L,NZ)=CCO2P*RootVH2O_pvr(N,L,NZ)
      plt_rbgc%trcg_air2root_flx__pvr(idg_CO2,N,L,NZ)=0._r8
      plt_rbgc%trcg_Root_DisEvap_flx_vr(idg_CO2,N,L,NZ)=0._r8
      plt_rbgc%RUPGasSol_vr(idg_CO2,N,L,NZ)=0._r8
      plt_rbgc%RCO2P(N,L,NZ)=0._r8
      COXYA=COXYE
      COXYP=0.032_r8*EXP(-6.175_r8-0.0211_r8*ATCA)*OXYE
      plt_rbgc%trcg_rootml_pvr(idg_O2,N,L,NZ)=COXYA*RootPoreVol_pvr(N,L,NZ)
      plt_rbgc%trcs_rootml_pvr(idg_O2,N,L,NZ)=COXYP*RootVH2O_pvr(N,L,NZ)
      plt_rbgc%trcg_rootml_pvr(idg_beg:idg_end-1,N,L,NZ)=0._r8
      plt_rbgc%trcs_rootml_pvr(idg_beg:idg_end-1,N,L,NZ)=0._r8
      plt_rbgc%RootAutoRO2Limiter_pvr(N,L,NZ)=1.0
      D30: DO NR=1,MaxNumRootAxes
        plt_morph%Root2ndXNum_rpvr(N,L,NR,NZ)=0._r8
        plt_morph%Root1stLen_rpvr(N,L,NR,NZ)=0._r8
        plt_morph%Root2ndLen_pvr(N,L,NR,NZ)=0._r8
        plt_morph%Root1stDepz_pft(N,NR,NZ)=SeedDepth_pft(NZ)
        plt_biom%Root1stStructElm_rpvr(1:NumPlantChemElms,N,L,NR,NZ)=0._r8
        plt_biom%Root2ndStructElm_pvr(1:NumPlantChemElms,N,L,NR,NZ)=0._r8
        plt_biom%Root1stElm_raxs(1:NumPlantChemElms,N,NR,NZ)=0._r8
      ENDDO D30
      IF(N.EQ.1)THEN
        D6400: DO K=1,pltpar%NumOfPlantLitrCmplxs
          plt_bgcr%LitfalChemElm_pvr(1:NumPlantChemElms,1:jsken,K,L,NZ)=0._r8
        ENDDO D6400
        plt_biom%RootNodulNonstElm_pvr(1:NumPlantChemElms,L,NZ)=0._r8
        plt_biom%RootNodulElm_pvr(1:NumPlantChemElms,L,NZ)=0._r8
        RootN2Fix_pvr(L,NZ)=0._r8
      ENDIF
    ENDDO D20
  ENDDO D40

  plt_rbgc%RootNutUptake_pvr(ids_NH4,1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RootNutUptake_pvr(ids_NH4B,1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RootNutUptake_pvr(ids_H2PO4,1:2,NL+1:JZ1,NZ)=0._r8
  plt_rbgc%RootNutUptake_pvr(ids_H2PO4B,1:2,NL+1:JZ1,NZ)=0._r8
  plt_morph%RootLenDensPerPlant_pvr(1:2,NL+1:JZ1,NZ)=0._r8
  end associate
  end subroutine InitRootMychorMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitSeedMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  REAL(R8) :: FDM

  associate(                             &
    PlantPopulation_pft          =>   plt_site%PlantPopulation_pft      , &
    PSICanopy_pft                =>   plt_ew%PSICanopy_pft    , &
    WatByPCanopy                 =>   plt_ew%WatByPCanopy     , &
    CanopyWater_pft              =>   plt_ew%CanopyWater_pft     , &
    RootFracRemobilizableBiom    =>   plt_allom%RootFracRemobilizableBiom  , &
    CNGR                         =>   plt_allom%CNGR   , &
    CPGR                         =>   plt_allom%CPGR   , &
    Root1stElm_raxs              =>   plt_biom%Root1stElm_raxs  , &
    SeedCPlanted_pft             =>   plt_biom%SeedCPlanted_pft   , &
    Root1stStructElm_rpvr        =>   plt_biom%Root1stStructElm_rpvr  , &
    LeafPetolBiomassC_brch       =>   plt_biom%LeafPetolBiomassC_brch   , &
    CanopyLeafShethC_pft         =>   plt_biom%CanopyLeafShethC_pft    , &
     PopuPlantRootC_vr           =>   plt_biom% PopuPlantRootC_vr   , &
    PetoleChemElm_brch           =>   plt_biom%PetoleChemElm_brch, &
    RootMycoActiveBiomC_pvr          =>   plt_biom%RootMycoActiveBiomC_pvr  , &
    RootProteinC_pvr             =>   plt_biom%RootProteinC_pvr   , &
     RootMycoNonstElm_pvr        =>   plt_biom%RootMycoNonstElm_pvr  , &
    NonstructElm_brch            =>   plt_biom%NonstructElm_brch  , &
    LeafChemElms_brch            =>   plt_biom%LeafChemElms_brch , &
    NonstructElms_pft            =>   plt_biom%NonstructElms_pft   , &
    SeedCMass                    =>   plt_morph%SeedCMass   , &
    Root1stDepz_pft              =>   plt_morph%Root1stDepz_pft  , &
    NGTopRootLayer_pft           =>   plt_morph%NGTopRootLayer_pft      &
  )
!
!     INITIALIZE SEED MORPHOLOGY AND BIOMASS
!
!     WTRVC,WTRVN,WTRVP=C,N,P in storage reserves (g)
!     WTLFB,WTLFBN,WTLFBP=C,N,P in leaves (g)
!     LeafPetolBiomassC_brch=C in leaves+petioles (g)
!     FDM-dry matter fraction (g DM C g FM C-1)
!     CanopyWater_pft,WatByPCanopy=water volume in,on canopy (m3)
!     CPOOL,ZPOOL,PPOOL=C,N,P in canopy nonstructural pools (g)
!     WTRT1,WTRT1N,WTRT1P=C,N,P in primary root layer (g)
!     RTWT1,RTWT1N,RTWT1P=total C,N,P in primary root (g)
!     RootMycoActiveBiomC_pvr, PopuPlantRootC_vr=total root C mass (g)
!     RootProteinC_pvr=total root protein C mass (g)
!     CPOOLR,ZPOOLR,PPOOLR=C,N,P in root,myco nonstructural pools (g)
!
  SeedCPlanted_pft(NZ)=SeedCMass(NZ)*PlantPopulation_pft(NZ)
  NonstructElms_pft(ielmc,NZ)=SeedCPlanted_pft(NZ)
  NonstructElms_pft(ielmn,NZ)=CNGR(NZ)*NonstructElms_pft(ielmc,NZ)
  NonstructElms_pft(ielmp,NZ)=CPGR(NZ)*NonstructElms_pft(ielmc,NZ)
  LeafChemElms_brch(ielmn,1,NZ)=CNGR(NZ)*LeafChemElms_brch(ielmc,1,NZ)
  LeafChemElms_brch(ielmp,1,NZ)=CPGR(NZ)*LeafChemElms_brch(ielmc,1,NZ)
  LeafPetolBiomassC_brch(1,NZ)=LeafChemElms_brch(ielmc,1,NZ)+PetoleChemElm_brch(ielmc,1,NZ)
  CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(1,NZ)
  FDM=AMIN1(1.0_r8,0.16_r8-0.045_r8*PSICanopy_pft(NZ))
  CanopyWater_pft(NZ)=ppmc*CanopyLeafShethC_pft(NZ)/FDM
  WatByPCanopy(NZ)=0._r8
  NonstructElm_brch(ielmn,1,NZ)=CNGR(NZ)*NonstructElm_brch(ielmc,1,NZ)
  NonstructElm_brch(ielmp,1,NZ)=CPGR(NZ)*NonstructElm_brch(ielmc,1,NZ)
  Root1stStructElm_rpvr(ielmn,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)=CNGR(NZ) &
    *Root1stStructElm_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  Root1stStructElm_rpvr(ielmp,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)=CPGR(NZ) &
    *Root1stStructElm_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  Root1stElm_raxs(ielmn,1,1,NZ)=CNGR(NZ)*Root1stElm_raxs(ielmc,1,1,NZ)
  Root1stElm_raxs(ielmp,1,1,NZ)=CPGR(NZ)*Root1stElm_raxs(ielmc,1,1,NZ)
  RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=&
    Root1stStructElm_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  PopuPlantRootC_vr(ipltroot,NGTopRootLayer_pft(NZ),NZ)= &
    Root1stStructElm_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  RootProteinC_pvr(1,NGTopRootLayer_pft(NZ),NZ)=RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)&
    *RootFracRemobilizableBiom(NZ)
  RootMycoNonstElm_pvr(ielmn,1,NGTopRootLayer_pft(NZ),NZ)=CNGR(NZ)&
    *RootMycoNonstElm_pvr(ielmc,1,NGTopRootLayer_pft(NZ),NZ)
  RootMycoNonstElm_pvr(ielmp,1,NGTopRootLayer_pft(NZ),NZ)=CPGR(NZ) &
    *RootMycoNonstElm_pvr(ielmc,1,NGTopRootLayer_pft(NZ),NZ)

  end associate
  end subroutine InitSeedMorphoBio

  end module StartqsMod
