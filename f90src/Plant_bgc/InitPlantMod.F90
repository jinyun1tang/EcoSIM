module InitPlantMod
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use UnitMod,          only: units
  use minimathmod,      only: AZMAX1
  use EcoSiMParDataMod, only: pltpar
  use EcosimConst
  use EcoSIMConfig
  use PlantAPIData
  use TracerIDMod
  use PlantMathFuncMod
  use GrosubPars

  use PlantMathFuncMod
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: StartPlants
  public :: InitPlantPhenoMorphoBio
  contains

  SUBROUTINE StartPlants(NZ1Q,NZ2Q)
!
!     THIS SUBROUTINE INITIALIZES ALL PLANT VARIABLES
!
  implicit none
  integer, intent(in) :: NZ1Q,NZ2Q

  integer :: K,L,M,NZ,NZ2X
!     begin_execution

  associate(                                             &
    AREA3               => plt_site%AREA3,               &  !input
    IsPlantActive_pft   => plt_pheno%IsPlantActive_pft,  &  !input
    NL                  => plt_site%NL,                  &  !input
    NP                  => plt_site%NP,                  &  !input
    NU                  => plt_site%NU,                  &  !input
    PlantPopulation_pft => plt_site%PlantPopulation_pft, &  !input
    ZERO                => plt_site%ZERO,                &  !input
    ZERO4Groth_pft      => plt_biom%ZERO4Groth_pft,      &  !output
    ZERO4LeafVar_pft    => plt_biom%ZERO4LeafVar_pft,    &  !output
    ZERO4Uptk_pft       => plt_rbgc%ZERO4Uptk_pft        &  !output
  )
!
!     INITIALIZE SHOOT GROWTH VARIABLES
!
!     IsPlantActive_pft=PFT flag:0=not active,1=active
!     iYearPlanting_pft,iDayPlanting_pft,iYearPlantHarvest_pft,iDayPlantHarvest_pft=year,day of planting,arvesting
!     PPI,PPX=initial,current population (m-2)
!     CF,ClumpFactorInit_pft=current,initial clumping factor
!     H2OCuticleResist_pft=cuticular resistance to water (h m-1)
!     CO2CuticleResist_pft=cuticular resistance to CO2 (s m-1)
!     CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios
!     RootFracRemobilizableBiom=maximum root protein concentration (g g-1)
!     O2I=intercellular O2 concentration in C3,C4 PFT (umol mol-1)
!

      NZ2X=MIN(NZ2Q,NP)
      D9985: DO NZ=NZ1Q,NZ2X

        IF(IsPlantActive_pft(NZ).EQ.iDormant)THEN

          call InitShootGrowth(NZ)

          call PlantLitterFraction(NZ)

          call PFTThermalAcclimation(NZ)

          call InitDimensionsandUptake(NZ)

          call InitPlantPhenoMorphoBio(NZ)

          call InitMassBalance(NZ)

          call InitRootMychorMorphoBio(NZ)

          call InitSeedMorphoBio(NZ)

          call InitPlantHeatWater(NZ)

        ENDIF
      ENDDO D9985

      DO NZ=NZ1Q,NZ2X
        ZERO4Groth_pft(NZ)   = ZERO*PlantPopulation_pft(NZ)
        ZERO4Uptk_pft(NZ)    = ZERO*PlantPopulation_pft(NZ)/AREA3(NU)
        ZERO4LeafVar_pft(NZ) = ZERO*PlantPopulation_pft(NZ)*1.0E+06_r8
      ENDDO  
!
!     FILL OUT UNUSED ARRAYS
!
      D9986: DO NZ=NP+1,JP1
        plt_bgcr%SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ) = 0._r8
        plt_bgcr%LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ)     = 0._r8
        plt_biom%StandDeadStrutElms_pft(1:NumPlantChemElms,NZ)         = 0._r8
        D6401: DO L=1,NL
          DO  K=1,pltpar%NumOfPlantLitrCmplxs
            plt_bgcr%LitrfalStrutElms_pvr(1:NumPlantChemElms,1:jsken,K,L,NZ)=0._r8
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
  associate(                                                          &
    ClumpFactorInit_pft       => plt_morph%ClumpFactorInit_pft,       &  !input
    CuticleResist_pft         => plt_photo%CuticleResist_pft,         &  !input
    PPatSeeding_pft           => plt_site%PPatSeeding_pft,            &  !input
    RootrNC_pft               => plt_allom%RootrNC_pft,               &  !input
    RootrPC_pft               => plt_allom%RootrPC_pft,               &  !input
    iHarvestDay_pft           => plt_distb%iHarvestDay_pft,           &  !input
    iHarvestYear_pft          => plt_distb%iHarvestYear_pft,          &  !input
    iPlantPhotosynthesisType  => plt_photo%iPlantPhotosynthesisType,  &  !input
    iPlantingDay_pft          => plt_distb%iPlantingDay_pft,          &  !input
    iPlantingYear_pft         => plt_distb%iPlantingYear_pft,         &  !input
    PPI_pft                   => plt_site%PPI_pft,                    &  !inoput
    rCNNonstRemob_pft         => plt_allom%rCNNonstRemob_pft,         &  !inoput
    rCPNonstRemob_pft         => plt_allom%rCPNonstRemob_pft,         &  !inoput
    CO2CuticleResist_pft      => plt_photo%CO2CuticleResist_pft,      &  !output
    ClumpFactor_pft           => plt_morph%ClumpFactor_pft,           &  !output
    H2OCuticleResist_pft      => plt_photo%H2OCuticleResist_pft,      &  !output
    O2I_pft                   => plt_photo%O2I_pft,                   &  !output
    PPX_pft                   => plt_site%PPX_pft,                    &  !output
    RootFracRemobilizableBiom => plt_allom%RootFracRemobilizableBiom, &  !output
    iDayPlantHarvest_pft      => plt_distb%iDayPlantHarvest_pft,      &  !output
    iDayPlanting_pft          => plt_distb%iDayPlanting_pft,          &  !output
    iYearPlantHarvest_pft     => plt_distb%iYearPlantHarvest_pft,     &  !output
    iYearPlanting_pft         => plt_distb%iYearPlanting_pft          &  !output
  )
  iYearPlanting_pft(NZ)     = iPlantingYear_pft(NZ)
  iDayPlanting_pft(NZ)      = iPlantingDay_pft(NZ)
  iYearPlantHarvest_pft(NZ) = iHarvestYear_pft(NZ)
  iDayPlantHarvest_pft(NZ)  = iHarvestDay_pft(NZ)
  PPI_pft(NZ)               = PPatSeeding_pft(NZ)
  PPX_pft(NZ)               = PPI_pft(NZ)
  ClumpFactor_pft(NZ)       = ClumpFactorInit_pft(NZ)

  H2OCuticleResist_pft(NZ)      = CuticleResist_pft(NZ)/3600.0_r8
  CO2CuticleResist_pft(NZ)      = CuticleResist_pft(NZ)*1.56_r8    !1.56 = sqrt(44./18.)
  rCNNonstRemob_pft(NZ)         = 2.5_r8
  rCPNonstRemob_pft(NZ)         = 25.0_r8
  RootFracRemobilizableBiom(NZ) = AMIN1(RootrNC_pft(NZ)*rCNNonstRemob_pft(NZ),&
    RootrPC_pft(NZ)*rCPNonstRemob_pft(NZ))
  IF(iPlantPhotosynthesisType(NZ).EQ.ic3_photo)THEN
    O2I_pft(NZ)=2.10E+05_r8
  ELSE
    O2I_pft(NZ)=3.96E+05_r8
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

  associate(                                                          &
    MatureGroup_pft           => plt_pheno%MatureGroup_pft,           &  !input
    NumLitterGroups           => pltpar%NumLitterGroups,              &  !input
    RefLeafAppearRate_pft     => plt_pheno%RefLeafAppearRate_pft,     &  !input
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft,        &  !input
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,     &  !input
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft, &  !input
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr,       &  !inoput
    FracGroth2Node_pft        => plt_allom%FracGroth2Node_pft,        &  !output
    NumCogrowthNode_pft       => plt_morph%NumCogrowthNode_pft,       &  !output
    icarbhyro                 => pltpar%icarbhyro,                    &  !output
    icellulos                 => pltpar%icellulos,                    &  !output
    icwood                    => pltpar%icwood,                       &  !output
    ifoliar                   => pltpar%ifoliar,                      &  !output
    ilignin                   => pltpar%ilignin,                      &  !output
    inonfoliar                => pltpar%inonfoliar,                   &  !output
    inonstruct                => pltpar%inonstruct,                   &  !output
    iprotein                  => pltpar%iprotein,                     &  !output
    iroot                     => pltpar%iroot,                        &  !output
    istalk                    => pltpar%istalk                        &  !output
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
  ElmAllocmat4Litr(ielmc,inonstruct,iprotein,NZ)  = 0.0_r8
  ElmAllocmat4Litr(ielmc,inonstruct,icarbhyro,NZ) = 0.67_r8
  ElmAllocmat4Litr(ielmc,inonstruct,icellulos,NZ) = 0.33_r8
  ElmAllocmat4Litr(ielmc,inonstruct,ilignin,NZ)   = 0.0_r8
!
!     NON-VASCULAR (E.G. MOSSES)
!
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
    ElmAllocmat4Litr(ielmc,ifoliar,iprotein,NZ)  = 0.07_r8
    ElmAllocmat4Litr(ielmc,ifoliar,icarbhyro,NZ) = 0.25_r8
    ElmAllocmat4Litr(ielmc,ifoliar,icellulos,NZ) = 0.30_r8
    ElmAllocmat4Litr(ielmc,ifoliar,ilignin,NZ)   = 0.38_r8

    ElmAllocmat4Litr(ielmc,inonfoliar,iprotein,NZ)  = 0.07_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,icarbhyro,NZ) = 0.25_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,icellulos,NZ) = 0.30_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,ilignin,NZ)   = 0.38_r8
!
!     LEGUMES
!
  ELSEIF(is_plant_N2fix(iPlantNfixType_pft(NZ)))THEN
    ElmAllocmat4Litr(ielmc,ifoliar,iprotein,NZ)  = 0.16_r8
    ElmAllocmat4Litr(ielmc,ifoliar,icarbhyro,NZ) = 0.38_r8
    ElmAllocmat4Litr(ielmc,ifoliar,icellulos,NZ) = 0.34_r8
    ElmAllocmat4Litr(ielmc,ifoliar,ilignin,NZ)   = 0.12_r8

    ElmAllocmat4Litr(ielmc,inonfoliar,iprotein,NZ)  = 0.07_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,icarbhyro,NZ) = 0.41_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,icellulos,NZ) = 0.37_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,ilignin,NZ)   = 0.15_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
    ElmAllocmat4Litr(ielmc,ifoliar,iprotein,NZ)  = 0.08_r8
    ElmAllocmat4Litr(ielmc,ifoliar,icarbhyro,NZ) = 0.41_r8
    ElmAllocmat4Litr(ielmc,ifoliar,icellulos,NZ) = 0.36_r8
    ElmAllocmat4Litr(ielmc,ifoliar,ilignin,NZ)   = 0.15_r8

    ElmAllocmat4Litr(ielmc,inonfoliar,iprotein,NZ)  = 0.07_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,icarbhyro,NZ) = 0.41_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,icellulos,NZ) = 0.36_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,ilignin,NZ)   = 0.16_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.1 &
    .OR.iPlantTurnoverPattern_pft(NZ).GE.3)THEN
    ElmAllocmat4Litr(ielmc,ifoliar,iprotein,NZ)  = 0.07_r8
    ElmAllocmat4Litr(ielmc,ifoliar,icarbhyro,NZ) = 0.34_r8
    ElmAllocmat4Litr(ielmc,ifoliar,icellulos,NZ) = 0.36_r8
    ElmAllocmat4Litr(ielmc,ifoliar,ilignin,NZ)   = 0.23_r8

    ElmAllocmat4Litr(ielmc,inonfoliar,iprotein,NZ)  = 0.0_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,icarbhyro,NZ) = 0.045_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,icellulos,NZ) = 0.660_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,ilignin,NZ)   = 0.295_r8
!
!     CONIFEROUS TREES
!
  ELSE
    ElmAllocmat4Litr(ielmc,ifoliar,iprotein,NZ)  = 0.07_r8
    ElmAllocmat4Litr(ielmc,ifoliar,icarbhyro,NZ) = 0.25_r8
    ElmAllocmat4Litr(ielmc,ifoliar,icellulos,NZ) = 0.38_r8
    ElmAllocmat4Litr(ielmc,ifoliar,ilignin,NZ)   = 0.30_r8

    ElmAllocmat4Litr(ielmc,inonfoliar,iprotein,NZ)  = 0.0_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,icarbhyro,NZ) = 0.045_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,icellulos,NZ) = 0.660_r8
    ElmAllocmat4Litr(ielmc,inonfoliar,ilignin,NZ)   = 0.295_r8
  ENDIF
!
!     FRACTIONS OF WOODY LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     NON-VASCULAR
!
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
    ElmAllocmat4Litr(ielmc,istalk,iprotein,NZ)  = 0.07_r8
    ElmAllocmat4Litr(ielmc,istalk,icarbhyro,NZ) = 0.25_r8
    ElmAllocmat4Litr(ielmc,istalk,icellulos,NZ) = 0.30_r8
    ElmAllocmat4Litr(ielmc,istalk,ilignin,NZ)   = 0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
    ElmAllocmat4Litr(ielmc,istalk,iprotein,NZ)  = 0.03_r8
    ElmAllocmat4Litr(ielmc,istalk,icarbhyro,NZ) = 0.25_r8
    ElmAllocmat4Litr(ielmc,istalk,icellulos,NZ) = 0.57_r8
    ElmAllocmat4Litr(ielmc,istalk,ilignin,NZ)   = 0.15_r8
!
!     DECIDUOUS AND CONIFEROUS TREES
!
  ELSE
    ElmAllocmat4Litr(ielmc,istalk,iprotein,NZ)  = 0.0_r8
    ElmAllocmat4Litr(ielmc,istalk,icarbhyro,NZ) = 0.045_r8
    ElmAllocmat4Litr(ielmc,istalk,icellulos,NZ) = 0.660_r8
    ElmAllocmat4Litr(ielmc,istalk,ilignin,NZ)   = 0.295_r8
  ENDIF
!
!     FRACTIONS OF FINE ROOT LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN PC&E 25:601-608
!
!     NON-VASCULAR
!
  IF(is_root_shallow(iPlantRootProfile_pft(NZ)))THEN
    ElmAllocmat4Litr(ielmc,iroot,iprotein,NZ)=0.07_r8
    ElmAllocmat4Litr(ielmc,iroot,icarbhyro,NZ)=0.25_r8
    ElmAllocmat4Litr(ielmc,iroot,icellulos,NZ)=0.30_r8
    ElmAllocmat4Litr(ielmc,iroot,ilignin,NZ)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
    ElmAllocmat4Litr(ielmc,iroot,iprotein,NZ)  = 0.057_r8
    ElmAllocmat4Litr(ielmc,iroot,icarbhyro,NZ) = 0.263_r8
    ElmAllocmat4Litr(ielmc,iroot,icellulos,NZ) = 0.542_r8
    ElmAllocmat4Litr(ielmc,iroot,ilignin,NZ)   = 0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ).EQ.1&
    .OR.iPlantTurnoverPattern_pft(NZ).GE.3)THEN
    ElmAllocmat4Litr(ielmc,iroot,iprotein,NZ)  = 0.059_r8
    ElmAllocmat4Litr(ielmc,iroot,icarbhyro,NZ) = 0.308_r8
    ElmAllocmat4Litr(ielmc,iroot,icellulos,NZ) = 0.464_r8
    ElmAllocmat4Litr(ielmc,iroot,ilignin,NZ)   = 0.169_r8
!
!     CONIFEROUS TREES
!
  ELSE
    ElmAllocmat4Litr(ielmc,iroot,iprotein,NZ)  = 0.059_r8
    ElmAllocmat4Litr(ielmc,iroot,icarbhyro,NZ) = 0.308_r8
    ElmAllocmat4Litr(ielmc,iroot,icellulos,NZ) = 0.464_r8
    ElmAllocmat4Litr(ielmc,iroot,ilignin,NZ)   = 0.169_r8
  ENDIF
!
!     COARSE WOODY LITTER FROM BOLES AND ROOTS
!
  ElmAllocmat4Litr(ielmc,icwood,iprotein,NZ)  = 0.00_r8
  ElmAllocmat4Litr(ielmc,icwood,icarbhyro,NZ) = 0.045_r8
  ElmAllocmat4Litr(ielmc,icwood,icellulos,NZ) = 0.660_r8
  ElmAllocmat4Litr(ielmc,icwood,ilignin,NZ)   = 0.295_r8
!
!     INITIALIZE C-N AND C-P RATIOS IN PLANT LITTER
!
!     CNOPC,CPOPC=fractions to allocate N,P to kinetic components
!     CFOPN,CFOPP=distribution of litter N,P to kinetic components
!
  CNOPC(iprotein)  = 0.020_r8
  CNOPC(icarbhyro) = 0.010_r8
  CNOPC(icellulos) = 0.010_r8
  CNOPC(ilignin)   = 0.020_r8

  CPOPC(iprotein)  = 0.0020_r8
  CPOPC(icarbhyro) = 0.0010_r8
  CPOPC(icellulos) = 0.0010_r8
  CPOPC(ilignin)   = 0.0020_r8

  D110: DO N=0,NumLitterGroups
    CNOPCT=0.0_r8
    CPOPCT=0.0_r8
    D100: DO M=1,jsken
      CNOPCT=CNOPCT+ElmAllocmat4Litr(ielmc,N,M,NZ)*CNOPC(M)
      CPOPCT=CPOPCT+ElmAllocmat4Litr(ielmc,N,M,NZ)*CPOPC(M)
    ENDDO D100
    D105: DO M=1,jsken
      ElmAllocmat4Litr(ielmn,N,M,NZ)=ElmAllocmat4Litr(ielmc,N,M,NZ)*CNOPC(M)/CNOPCT
      ElmAllocmat4Litr(ielmp,N,M,NZ)=ElmAllocmat4Litr(ielmc,N,M,NZ)*CPOPC(M)/CPOPCT
    ENDDO D105
  ENDDO D110
!
!     CONCURRENT NODE GROWTH
!
!     FracGroth2Node_pft=scales node number for perennial vegetation (e.g. trees)
!     NumCogrowthNode_pft=number of concurrently growing nodes
!     iPlantTurnoverPattern_pft=turnover:0=all aboveground,1=all leaf+petiole,2=none,3=between 1,2!
  IF(iPlantTurnoverPattern_pft(NZ).EQ.0 &
    .OR.(.not.is_plant_treelike(iPlantRootProfile_pft(NZ))))THEN
    FracGroth2Node_pft(NZ)=1.0_r8
    IF(MatureGroup_pft(NZ).LE.10)THEN
      NumCogrowthNode_pft(NZ)=3
    ELSEIF(MatureGroup_pft(NZ).LE.15)THEN
      NumCogrowthNode_pft(NZ)=4
    ELSE
      NumCogrowthNode_pft(NZ)=5
    ENDIF
  ELSE
    !not grasslike plant
    FracGroth2Node_pft(NZ)=AMAX1(1.0_r8,0.04_r8/RefLeafAppearRate_pft(NZ))
    NumCogrowthNode_pft(NZ)=24
  ENDIF
  end associate
  end subroutine PlantLitterFraction
!------------------------------------------------------------------------------------------

  subroutine PFTThermalAcclimation(NZ)

  implicit none
  integer, intent(in) :: NZ
  real(r8), parameter :: TCZD = 5.0_r8        !basal value for threshold temperature for spring leafout/dehardening	oC
  real(r8), parameter :: TCXD = 12.0_r8       !basal value for threshold temperature for autumn leafoff/hardening	oC

  associate(                                                          &
    DATAP                     => plt_site%DATAP,                      &  !input
    PlantInitThermoAdaptZone  => plt_pheno%PlantInitThermoAdaptZone,  &  !input
    iPlantPhotosynthesisType  => plt_photo%iPlantPhotosynthesisType,  &  !input
    TempOffset_pft            => plt_pheno%TempOffset_pft,            &  !inoput
    iPlantThermoAdaptZone_pft => plt_pheno%iPlantThermoAdaptZone_pft, &  !inoput
    HighTempLimitSeed_pft     => plt_pheno%HighTempLimitSeed_pft,     &  !output
    SeedTempSens_pft          => plt_pheno%SeedTempSens_pft,          &  !output
    TC4LeafOff_pft            => plt_pheno%TC4LeafOff_pft,            &  !output
    TC4LeafOut_pft            => plt_pheno%TC4LeafOut_pft             &  !output
  )
!
!     PFT THERMAL ACCLIMATION
!
!     ZTYP,PlantInitThermoAdaptZone=dynamic,initial thermal adaptation zone from PFT file
!     TempOffset_pft=shift in Arrhenius curve for thermal adaptation (oC)
!     TCZ,TC4LeafOff_pft=threshold temperature for leafout,leafoff
!     HTC=high temperature threshold for grain number loss (oC)
!     SeedTempSens_pft=sensitivity to HTC (seeds oC-1 above HTC)
!
  iPlantThermoAdaptZone_pft(NZ) = PlantInitThermoAdaptZone(NZ)
  TempOffset_pft(NZ)            = 2.667_r8*(2.5_r8-iPlantThermoAdaptZone_pft(NZ))
  TC4LeafOut_pft(NZ)            = TCZD-TempOffset_pft(NZ)
  TC4LeafOff_pft(NZ)            = AMIN1(15.0_r8,TCXD-TempOffset_pft(NZ))
  IF(iPlantPhotosynthesisType(NZ).EQ.ic3_photo)THEN
    IF(DATAP(NZ)(1:4).EQ.'soyb')THEN
      HighTempLimitSeed_pft(NZ)=30.0_r8+3.0_r8*iPlantThermoAdaptZone_pft(NZ)
      SeedTempSens_pft(NZ)=0.002_r8
    ELSE
      HighTempLimitSeed_pft(NZ)=27.0_r8+3.0_r8*iPlantThermoAdaptZone_pft(NZ)
      SeedTempSens_pft(NZ)=0.002_r8
    ENDIF
  ELSE
    HighTempLimitSeed_pft(NZ)=27.0_r8+3.0_r8*iPlantThermoAdaptZone_pft(NZ)
    SeedTempSens_pft(NZ)=0.005_r8
  ENDIF
  end associate
  end subroutine PFTThermalAcclimation
!------------------------------------------------------------------------------------------

  subroutine InitDimensionsandUptake(NZ)

  implicit none
  integer, intent(in) :: NZ
  INTEGER :: L,N,NR
  associate(                                                  &
    CumSoilThickness_vr   => plt_site%CumSoilThickness_vr,    &  !input
    NL                    => plt_site%NL,                     &  !input
    NU                    => plt_site%NU,                     &  !input
    PlantinDepz_pft       => plt_morph%PlantinDepz_pft,       &  !input
    RootBiomGrosYld_pft   => plt_allom%RootBiomGrosYld_pft,   &  !input
    RootrNC_pft           => plt_allom%RootrNC_pft,           &  !input
    RootrPC_pft           => plt_allom%RootrPC_pft,           &  !input
    SeedAreaMean_pft      => plt_morph%SeedAreaMean_pft,      &  !input
    SeedCMass_pft         => plt_morph%SeedCMass_pft,         &  !input
    SeedMeanLen_pft       => plt_morph%SeedMeanLen_pft,       &  !input
    SeedVolumeMean_pft    => plt_morph%SeedVolumeMean_pft,    &  !input
    CMinNH4Root_pft       => plt_rbgc%CMinNH4Root_pft,        &  !inoput
    CMinPO4Root_pft       => plt_rbgc%CMinPO4Root_pft,        &  !inoput
    CminNO3Root_pft       => plt_rbgc%CminNO3Root_pft,        &  !inoput
    KmNH4Root_pft         => plt_rbgc%KmNH4Root_pft,          &  !inoput
    KmNO3Root_pft         => plt_rbgc%KmNO3Root_pft,          &  !inoput
    KmPO4Root_pft         => plt_rbgc%KmPO4Root_pft,          &  !inoput
    Root1stMaxRadius1_pft => plt_morph%Root1stMaxRadius1_pft, &  !inoput
    Root1stMaxRadius_pft  => plt_morph%Root1stMaxRadius_pft,  &  !inoput
    Root2ndMaxRadius1_pft => plt_morph%Root2ndMaxRadius1_pft, &  !inoput
    Root2ndMaxRadius_pft  => plt_morph%Root2ndMaxRadius_pft,  &  !inoput
    RootPorosity_pft      => plt_morph%RootPorosity_pft,      &  !inoput
    RootVolPerMassC_pft   => plt_morph%RootVolPerMassC_pft,   &  !inoput
    SeedDepth_pft         => plt_morph%SeedDepth_pft,         &  !inoput
    VmaxNH4Root_pft       => plt_rbgc%VmaxNH4Root_pft,        &  !inoput
    VmaxNO3Root_pft       => plt_rbgc%VmaxNO3Root_pft,        &  !inoput
    VmaxPO4Root_pft       => plt_rbgc%VmaxPO4Root_pft,        &  !inoput
    CNRTS_pft             => plt_allom%CNRTS_pft,             &  !output
    CPRTS_pft             => plt_allom%CPRTS_pft,             &  !output
    NGTopRootLayer_pft    => plt_morph%NGTopRootLayer_pft,    &  !output
    NIXBotRootLayer_pft   => plt_morph%NIXBotRootLayer_pft,   &  !output
    NIXBotRootLayer_rpft  => plt_morph%NIXBotRootLayer_rpft,  &  !output
    Root1stSpecLen_pft    => plt_morph%Root1stSpecLen_pft,    &  !output
    Root1stXSecArea_pft   => plt_morph%Root1stXSecArea_pft,   &  !output
    Root2ndSpecLen_pft    => plt_morph%Root2ndSpecLen_pft,    &  !output
    Root2ndXSecArea_pft   => plt_morph%Root2ndXSecArea_pft,   &  !output
    RootAxialResist_pft   => plt_morph%RootAxialResist_pft,   &  !output
    RootPoreTortu4Gas_pft => plt_morph%RootPoreTortu4Gas_pft, &  !output
    RootRadialResist_pft  => plt_morph%RootRadialResist_pft,  &  !output
    RootRaidus_rpft       => plt_morph%RootRaidus_rpft        &  !output
  )
!
!     SEED CHARACTERISTICS
!
!     SeedVolumeMean_pft,SeedMeanLen_pft,SeedAreaMean_pft=seed volume(m3),length(m),AREA3(NU)(m2)
!     SeedCMass=seed C mass (g) from PFT file
!
  call calc_seed_geometry(SeedCMass_pft(NZ),SeedVolumeMean_pft(NZ),SeedMeanLen_pft(NZ),SeedAreaMean_pft(NZ))
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) DIMENSIONS, UPTAKE PARAMETERS
!
!     SeedDepth_pft=seeding depth(m) from PFT management file
!     CumSoilThickness_vr=depth to soil layer bottom from surface(m)
!     NG,NIX,NIXBotRootLayer_rpft=seeding,upper,lower rooting layer
!     CNRTS_pft,CPRTS_pft=N,P root growth yield
!     Root1stMaxRadius_pft,Root2ndMaxRadius_pft=maximum primary,secondary mycorrhizal radius (m)
!     PORT=mycorrhizal porosity
!     VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake(g m-2 h-1),Km(uM), min concn (uM)
!     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake(g m-2 h-1),Km(uM),min concn (uM)
!     RootRadialResist_pft,RootAxialResist_pft=radial,axial root resistivity (m2 MPa-1 h-1)
!
  SeedDepth_pft(NZ)=PlantinDepz_pft(NZ)
  D9795: DO L=NU,NL
    IF(SeedDepth_pft(NZ).GE.CumSoilThickness_vr(L-1) &
      .AND.SeedDepth_pft(NZ).LT.CumSoilThickness_vr(L))THEN
      NGTopRootLayer_pft(NZ)  = L
      NIXBotRootLayer_pft(NZ) = L
      D9790: DO NR=1,pltpar%MaxNumRootAxes
        NIXBotRootLayer_rpft(NR,NZ)=L
      ENDDO D9790
    ENDIF
  ENDDO D9795
  CNRTS_pft(NZ)              = RootrNC_pft(NZ)*RootBiomGrosYld_pft(NZ)
  CPRTS_pft(NZ)              = RootrPC_pft(NZ)*RootBiomGrosYld_pft(NZ)
  Root1stMaxRadius_pft(2,NZ) = 5.0E-06_r8
  Root2ndMaxRadius_pft(2,NZ) = 5.0E-06_r8
  RootPorosity_pft(2,NZ)     = RootPorosity_pft(1,NZ)
  VmaxNH4Root_pft(2,NZ)      = VmaxNH4Root_pft(1,NZ)
  KmNH4Root_pft(2,NZ)        = KmNH4Root_pft(1,NZ)
  CMinNH4Root_pft(2,NZ)      = CMinNH4Root_pft(1,NZ)
  VmaxNO3Root_pft(2,NZ)      = VmaxNO3Root_pft(1,NZ)
  KmNO3Root_pft(2,NZ)        = KmNO3Root_pft(1,NZ)
  CminNO3Root_pft(2,NZ)      = CminNO3Root_pft(1,NZ)
  VmaxPO4Root_pft(2,NZ)      = VmaxPO4Root_pft(1,NZ)
  KmPO4Root_pft(2,NZ)        = KmPO4Root_pft(1,NZ)
  CMinPO4Root_pft(2,NZ)      = CMinPO4Root_pft(1,NZ)
  RootRadialResist_pft(2,NZ) = 1.0E+04_r8
  RootAxialResist_pft(2,NZ)  = 1.0E+12_r8
!
!     RootPoreTortu4Gas_pft=tortuosity for gas transport
!     RootRaidus_rpft=path length for radial diffusion within root (m)
!     RootVolPerMassC_pft=volume:C ratio (m3 g-1)
!     Root1stSpecLen_pft,Root2ndSpecLen_pft=specific primary,secondary root length (m g-1)
!     Root1stXSecArea_pft,Root2ndXSecArea_pft=specific primary,secondary root area (m2 g-1)
!
  D500: DO N=1,2
    RootPoreTortu4Gas_pft(N,NZ)     = RootPorosity_pft(N,NZ)**1.33_r8
    RootRaidus_rpft(N,NZ)       = LOG(1.0_r8/SQRT(AMAX1(0.01_r8,RootPorosity_pft(N,NZ))))
    RootVolPerMassC_pft(N,NZ)   = ppmc/(0.05_r8*(1.0_r8-RootPorosity_pft(N,NZ)))
    Root1stSpecLen_pft(N,NZ)    = RootVolPerMassC_pft(N,NZ)/(PICON*Root1stMaxRadius_pft(N,NZ)**2._r8)
    Root2ndSpecLen_pft(N,NZ)    = RootVolPerMassC_pft(N,NZ)/(PICON*Root2ndMaxRadius_pft(N,NZ)**2._r8)
    Root1stMaxRadius1_pft(N,NZ) = Root1stMaxRadius_pft(N,NZ)
!    2*SQRT(0.25*(1.0-RootPorosity_pft(N,NZ)))
    Root2ndMaxRadius1_pft(N,NZ)=Root2ndMaxRadius_pft(N,NZ)
!    2*SQRT(0.25*(1.0-RootPorosity_pft(N,NZ)))
    Root1stXSecArea_pft(N,NZ)=PICON*Root1stMaxRadius1_pft(N,NZ)**2
    Root2ndXSecArea_pft(N,NZ)=PICON*Root2ndMaxRadius1_pft(N,NZ)**2
  ENDDO D500
  end associate
  end subroutine InitDimensionsandUptake
!------------------------------------------------------------------------------------------

  subroutine InitPlantPhenoMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: K,L,M,N,NB
  associate(                                                                          &
    AREA3                             => plt_site%AREA3,                              &  !input
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft,                   &  !input
    NU                                => plt_site%NU,                                 &  !input
    PPX_pft                           => plt_site%PPX_pft,                            &  !input
    PetioleChemElmRemobFlx_brch       => plt_pheno%PetioleChemElmRemobFlx_brch,       &  !input
    ShootNodeNumAtPlanting_pft        => plt_morph%ShootNodeNumAtPlanting_pft,        &  !input
    iPlantBranchState_brch            => plt_pheno%iPlantBranchState_brch,            &  !input
    Hours4LenthenPhotoPeriod_brch     => plt_pheno%Hours4LenthenPhotoPeriod_brch,     &  !inoput
    Hours4ShortenPhotoPeriod_brch     => plt_pheno%Hours4ShortenPhotoPeriod_brch,     &  !inoput
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch,                 &  !inoput
    BranchNumber_brch                 => plt_morph%BranchNumber_brch,                 &  !output
    BranchNumber_pft                  => plt_morph%BranchNumber_pft,                  &  !output
    C4PhotosynDowreg_brch             => plt_photo%C4PhotosynDowreg_brch,             &  !output
    CMassCO2BundleSheath_node         => plt_photo%CMassCO2BundleSheath_node,         &  !output
    CMassHCO3BundleSheath_node        => plt_photo%CMassHCO3BundleSheath_node,        &  !output
    CPOOL3_node                       => plt_photo%CPOOL3_node,                       &  !output
    CPOOL4_node                       => plt_photo%CPOOL4_node,                       &  !output
    CanPBranchHeight                  => plt_morph%CanPBranchHeight,                  &  !output
    CanopyHeight_pft                  => plt_morph%CanopyHeight_pft,                  &  !output
    CanopyLeafAreaZ_pft               => plt_morph%CanopyLeafAreaZ_pft,               &  !output
    CanopyLeafArea_lnode              => plt_morph%CanopyLeafArea_lnode,              &  !output
    CanopyLeafArea_pft                => plt_morph%CanopyLeafArea_pft,                &  !output
    CanopyStalkArea_lbrch             => plt_morph%CanopyStalkArea_lbrch,             &  !output
    CanopyStemAreaZ_pft               => plt_morph%CanopyStemAreaZ_pft,               &  !output
    CanopyStemArea_pft                => plt_morph%CanopyStemArea_pft,                &  !output
    ChillHours_pft                    => plt_photo%ChillHours_pft,                    &  !output
    HourFailGrainFill_brch            => plt_pheno%HourFailGrainFill_brch,            &  !output
    Hours2LeafOut_brch                => plt_pheno%Hours2LeafOut_brch,                &  !output
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch,                &  !output
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch,                &  !output
    HoursDoingRemob_brch              => plt_pheno%HoursDoingRemob_brch,              &  !output
    HoursTooLowPsiCan_pft             => plt_pheno%HoursTooLowPsiCan_pft,             &  !output
    HypoctoHeight_pft                 => plt_morph%HypoctoHeight_pft,                 &  !output
    KHiestGroLeafNode_brch            => plt_pheno%KHiestGroLeafNode_brch,            &  !output
    KLeafNumber_brch                  => plt_morph%KLeafNumber_brch,                  &  !output
    KLowestGroLeafNode_brch           => plt_pheno%KLowestGroLeafNode_brch,           &  !output
    KMinNumLeaf4GroAlloc_brch         => plt_morph%KMinNumLeaf4GroAlloc_brch,         &  !output
    LeafAreaDying_brch                => plt_morph%LeafAreaDying_brch,                &  !output
    LeafAreaLive_brch                 => plt_morph%LeafAreaLive_brch,                 &  !output
    LeafAreaZsec_brch                 => plt_morph%LeafAreaZsec_brch,                 &  !output
    LeafNodeArea_brch                 => plt_morph%LeafNodeArea_brch,                 &  !output
    LeafNumberAtFloralInit_brch       => plt_pheno%LeafNumberAtFloralInit_brch,       &  !output
    LiveInterNodeHight_brch           => plt_morph%LiveInterNodeHight_brch,           &  !output
    MatureGroup_brch                  => plt_pheno%MatureGroup_brch,                  &  !output
    NodeNum2InitFloral_brch           => plt_morph%NodeNum2InitFloral_brch,           &  !output
    NodeNumNormByMatgrp_brch          => plt_pheno%NodeNumNormByMatgrp_brch,          &  !output
    NodeNumberAtAnthesis_brch         => plt_morph%NodeNumberAtAnthesis_brch,         &  !output
    NumOfBranches_pft                 => plt_morph%NumOfBranches_pft,                 &  !output
    NumOfLeaves_brch                  => plt_morph%NumOfLeaves_brch,                  &  !output
    PlantPopulation_pft               => plt_site%PlantPopulation_pft,                &  !output
    PotentialSeedSites_brch           => plt_morph%PotentialSeedSites_brch,           &  !output
    ReprodNodeNumNormByMatrgrp_brch   => plt_pheno%ReprodNodeNumNormByMatrgrp_brch,   &  !output
    RubiscoActivity_brch              => plt_photo%RubiscoActivity_brch,              &  !output
    SeedNumSet_brch                   => plt_morph%SeedNumSet_brch,                   &  !output
    StemAreaZsec_brch                 => plt_morph%StemAreaZsec_brch,                 &  !output
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch, &  !output
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch,     &  !output
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch                &  !output
  )
!
!     INITIALIZE PLANT PHENOLOGY
!
!     PP=population (grid cell-1)
!
  PlantPopulation_pft(NZ)=PPX_pft(NZ)*AREA3(NU)
  plt_pheno%doInitPlant_pft(NZ)      = ifalse
  plt_pheno%iPlantShootState_pft(NZ) = iLive
  plt_pheno%iPlantRootState_pft(NZ)  = iLive
  BranchNumber_pft(NZ)               = 0
  NumOfBranches_pft(NZ)              = 0
  HypoctoHeight_pft(NZ)              = 0._r8
  CanopyHeight_pft(NZ)               = 0._r8
  D10: DO NB=1,MaxNumBranches
    plt_pheno%doInitLeafOut_brch(NB,NZ)           = iEnable
    plt_pheno%doPlantLeafOut_brch(NB,NZ)          = iEnable
    plt_pheno%doPlantLeaveOff_brch(NB,NZ)         = iEnable
    plt_pheno%Prep4Literfall_brch(NB,NZ)          = ifalse
    plt_pheno%Hours4LiterfalAftMature_brch(NB,NZ) = 0
    MatureGroup_brch(NB,NZ)                       = MatureGroup_pft(NZ)
    ShootNodeNum_brch(NB,NZ)                      = ShootNodeNumAtPlanting_pft(NZ)
    NodeNum2InitFloral_brch(NB,NZ)                = ShootNodeNum_brch(NB,NZ)
    NodeNumberAtAnthesis_brch(NB,NZ)              = 0._r8
    NumOfLeaves_brch(NB,NZ)                       = 0._r8
    LeafNumberAtFloralInit_brch(NB,NZ)            = 0._r8
    KLeafNumber_brch(NB,NZ)                       = 1
    KMinNumLeaf4GroAlloc_brch(NB,NZ)              = 1
    KHiestGroLeafNode_brch(NB,NZ)                 = 1
    KLowestGroLeafNode_brch(NB,NZ)                = 0
    NodeNumNormByMatgrp_brch(NB,NZ)               = 0._r8
    ReprodNodeNumNormByMatrgrp_brch(NB,NZ)        = 0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ)          = 0._r8
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ)      = 0._r8
    Hours4LenthenPhotoPeriod_brch(NB,NZ)          = 0._r8
    Hours4ShortenPhotoPeriod_brch(NB,NZ)          = 0._r8
    Hours4Leafout_brch(NB,NZ)                     = Hours4LenthenPhotoPeriod_brch(NB,NZ)
    Hours4LeafOff_brch(NB,NZ)                     = Hours4ShortenPhotoPeriod_brch(NB,NZ)
    Hours2LeafOut_brch(NB,NZ)                     = 0._r8
    RubiscoActivity_brch(NB,NZ)                   = 1.0
    C4PhotosynDowreg_brch(NB,NZ)                  = 1.0
    HourFailGrainFill_brch(NB,NZ)                 = 0
    HoursDoingRemob_brch(NB,NZ)                   = 0
    BranchNumber_brch(NB,NZ)                      = 0
    plt_pheno%iPlantBranchState_brch(NB,NZ)       = iDead
    D15: DO M=1,pltpar%NumGrowthStages
      iPlantCalendar_brch(M,NB,NZ)=0
    ENDDO D15
  ENDDO D10
!
!     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
!
  HoursTooLowPsiCan_pft(NZ)=0._r8
  ChillHours_pft(NZ)=0._r8
  plt_biom%CanopyNonstElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)         = 0._r8
  plt_biom%CanopyNodulNonstElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)    = 0._r8
  plt_biom%PetoleStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)         = 0._r8
  plt_biom%ShootStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)          = 0._r8
  plt_biom%StalkStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)          = 0._r8
  plt_biom%LeafStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)           = 0._r8
  plt_biom%StalkRsrvElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)           = 0._r8
  plt_biom%HuskStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)           = 0._r8
  plt_biom%GrainStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)          = 0._r8
  plt_biom%EarStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)            = 0._r8
  plt_biom%CanopyNodulStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)    = 0._r8
  plt_pheno%LeafElmntRemobFlx_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)      = 0._r8
  plt_biom%LeafChemElmRemob_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)        = 0._r8
  plt_pheno%PetioleChemElmRemobFlx_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ) = 0._r8
  plt_biom%PetioleChemElmRemob_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)     = 0._r8
  plt_biom%SenecStalkStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ)     = 0._r8

  D25: DO NB=1,MaxNumBranches
    plt_biom%StalkLiveBiomassC_brch(NB,NZ)   = 0._r8
    plt_biom%LeafPetolBiomassC_brch(NB,NZ)   = 0._r8
    PotentialSeedSites_brch(NB,NZ)           = 0._r8
    SeedNumSet_brch(NB,NZ)                   = 0._r8
    plt_allom%GrainSeedBiomCMean_brch(NB,NZ) = 0._r8
    LeafAreaLive_brch(NB,NZ)                 = 0._r8
    plt_rbgc%NH3Dep2Can_brch(NB,NZ)          = 0._r8
    LeafAreaDying_brch(NB,NZ)                = 0._r8
    CanPBranchHeight(NB,NZ)                  = 0._r8
    D5: DO L=1,NumOfCanopyLayers1
      CanopyStalkArea_lbrch(L,NB,NZ)=0._r8
      DO N=1,NumOfLeafZenithSectors1
        StemAreaZsec_brch(N,L,NB,NZ)=0._r8
      enddo
    ENDDO D5

    DO K=0,MaxNodesPerBranch1
      LeafNodeArea_brch(K,NB,NZ)                                   = 0._r8
      LiveInterNodeHight_brch(K,NB,NZ)                             = 0._r8
      plt_morph%InternodeHeightDead_brch(K,NB,NZ)                  = 0._r8
      plt_morph%PetoleLensNode_brch(K,NB,NZ)                       = 0._r8
      plt_biom%LeafProteinCNode_brch(K,NB,NZ)                      = 0._r8
      plt_biom%PetoleProteinCNode_brch(K,NB,NZ)                    = 0._r8
      plt_biom%LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)      = 0._r8
      plt_biom%InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8
      plt_biom%PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)   = 0._r8

      D55: DO L=1,NumOfCanopyLayers1
        CanopyLeafArea_lnode(L,K,NB,NZ)=0._r8
        plt_biom%LeafElmsByLayerNode_brch(1:NumPlantChemElms,L,K,NB,NZ)=0._r8
      ENDDO D55

      IF(K.NE.0)THEN
        CPOOL3_node(K,NB,NZ)                = 0._r8
        CMassCO2BundleSheath_node(K,NB,NZ)  = 0._r8
        CMassHCO3BundleSheath_node(K,NB,NZ) = 0._r8
        CPOOL4_node(K,NB,NZ)                = 0._r8
        D45: DO L=1,NumOfCanopyLayers1
          DO N=1,NumOfLeafZenithSectors1
            LeafAreaZsec_brch(N,L,K,NB,NZ)=0._r8
          enddo
        ENDDO D45
      ENDIF
    enddo
  ENDDO D25
  
  D35: DO L=1,NumOfCanopyLayers1
    CanopyLeafAreaZ_pft(L,NZ)         = 0._r8
    plt_biom%CanopyLeafCLyr_pft(L,NZ) = 0._r8
    CanopyStemAreaZ_pft(L,NZ)         = 0._r8
  ENDDO D35
  plt_biom%CanopyNonstElms_pft(1:NumPlantChemElms,NZ)    = 0._r8
  plt_biom%CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ) = 0._r8
  plt_biom%NoduleNonstructCconc_pft(NZ)                  = 0._r8
  plt_biom%ShootStrutElms_pft(1:NumPlantChemElms,NZ)     = 0._r8
  plt_biom%LeafStrutElms_pft(1:NumPlantChemElms,NZ)      = 0._r8
  plt_biom%PetoleStrutElms_pft(1:NumPlantChemElms,NZ)    = 0._r8
  plt_biom%StalkStrutElms_pft(1:NumPlantChemElms,NZ)     = 0._r8
  plt_biom%CanopyStalkC_pft(NZ)                          = 0._r8
  plt_biom%StalkRsrvElms_pft(1:NumPlantChemElms,NZ)      = 0._r8
  plt_biom%HuskStrutElms_pft(1:NumPlantChemElms,NZ)      = 0._r8
  plt_biom%EarStrutElms_pft(1:NumPlantChemElms,NZ)       = 0._r8
  plt_biom%GrainStrutElms_pft(1:NumPlantChemElms,NZ)     = 0._r8
  plt_biom%RootStrutElms_pft(1:NumPlantChemElms,NZ)      = 0._r8
  plt_biom%NodulStrutElms_pft(1:NumPlantChemElms,NZ)     = 0._r8
  plt_biom%CanopyLeafShethC_pft(NZ)                      = 0._r8
  CanopyLeafArea_pft(NZ)                                 = 0._r8
  plt_biom%RootBiomCPerPlant_pft(NZ)                     = 0._r8
  CanopyStemArea_pft(NZ)                                 = 0._r8
  end associate
  end subroutine InitPlantPhenoMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitMassBalance(NZ)

  implicit none
  integer, intent(in) :: NZ
  integer :: M,NE
  real(r8) :: WTSTDX

  associate(                                                                   &
    AREA3                          => plt_site%AREA3,                          &  !input
    ElmAllocmat4Litr               => plt_soilchem%ElmAllocmat4Litr,           &  !input
    NU                             => plt_site%NU,                             &  !input
    StandingDeadInitC_pft          => plt_biom%StandingDeadInitC_pft,          &  !input
    icwood                         => pltpar%icwood,                           &  !input
    rNCStalk_pft                   => plt_allom%rNCStalk_pft,                  &  !input
    rPCStalk_pft                   => plt_allom%rPCStalk_pft,                  &  !input
    StandDeadKCompElms_pft         => plt_biom%StandDeadKCompElms_pft,         &  !inoput
    StandDeadStrutElms_pft         => plt_biom%StandDeadStrutElms_pft,         &  !inoput
    CanopyRespC_CumYr_pft          => plt_bgcr%CanopyRespC_CumYr_pft,          &  !output
    ETCanopy_CumYr_pft             => plt_ew%ETCanopy_CumYr_pft,               &  !output
    GrossCO2Fix_pft                => plt_bgcr%GrossCO2Fix_pft,                &  !output
    GrossResp_pft                  => plt_bgcr%GrossResp_pft,                  &  !output
    LitrfalStrutElms_CumYr_pft     => plt_bgcr%LitrfalStrutElms_CumYr_pft,     &  !output
    NH3Dep2Can_pft                 => plt_bgcr%NH3Dep2Can_pft,                 &  !output
    NH3Emis_CumYr_pft              => plt_bgcr%NH3Emis_CumYr_pft,              &  !output
    NetCumElmntFlx2Plant_pft       => plt_pheno%NetCumElmntFlx2Plant_pft,      &  !output
    PlantExudElm_CumYr_pft         => plt_rbgc%PlantExudElm_CumYr_pft,         &  !output
    PlantN2Fix_CumYr_pft           => plt_bgcr%PlantN2Fix_CumYr_pft,           &  !output
    RootUptk_N_CumYr_pft           => plt_rbgc%RootUptk_N_CumYr_pft,           &  !output
    RootUptk_P_CumYr_pft           => plt_rbgc%RootUptk_P_CumYr_pft,           &  !output
    SurfLitrfalStrutElms_CumYr_pft => plt_bgcr%SurfLitrfalStrutElms_CumYr_pft  &  !output
  )
!
!     INITIALIZE MASS BALANCE CHECKS
!
  IF(.not.is_restart().AND.is_first_year)THEN
    GrossCO2Fix_pft(NZ)                                   = 0._r8
    SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ) = 0._r8
    GrossResp_pft(NZ)                                     = 0._r8
    CanopyRespC_CumYr_pft(NZ)                             = 0._r8
    PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ)         = 0._r8
    RootUptk_N_CumYr_pft(NZ)                              = 0._r8
    RootUptk_P_CumYr_pft(NZ)                              = 0._r8
    LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ)     = 0._r8

    PlantN2Fix_CumYr_pft(NZ)                                 = 0._r8
    NH3Dep2Can_pft(NZ)                                       = 0._r8
    NH3Emis_CumYr_pft(NZ)                                    = 0._r8
    plt_distb%CO2ByFire_CumYr_pft(NZ)                        = 0._r8
    plt_distb%CH4ByFire_CumYr_pft(NZ)                        = 0._r8
    plt_distb%O2ByFire_CumYr_pft(NZ)                         = 0._r8
    plt_distb%NH3byFire_CumYr_pft(NZ)                        = 0._r8
    plt_distb%N2ObyFire_CumYr_pft(NZ)                        = 0._r8
    plt_distb%PO4byFire_CumYr_pft(NZ)                        = 0._r8
    plt_distb%EcoHavstElmntCum_pft(1:NumPlantChemElms,NZ)    = 0._r8
    plt_distb%EcoHavstElmnt_CumYr_pft(1:NumPlantChemElms,NZ) = 0._r8
    NetCumElmntFlx2Plant_pft(1:NumPlantChemElms,NZ)          = 0._r8
    StandDeadStrutElms_pft(1:NumPlantChemElms,NZ)            = 0._r8
    ETCanopy_CumYr_pft(NZ)                                   = 0._r8
    WTSTDX                                                   = StandingDeadInitC_pft(NZ)*AREA3(NU)
    D155: DO M=1,jsken
      StandDeadKCompElms_pft(ielmc,M,NZ)=WTSTDX*ElmAllocmat4Litr(ielmc,icwood,M,NZ)
      StandDeadKCompElms_pft(ielmn,M,NZ)=WTSTDX*rNCStalk_pft(NZ)*ElmAllocmat4Litr(ielmn,icwood,M,NZ)
      StandDeadKCompElms_pft(ielmp,M,NZ)=WTSTDX*rPCStalk_pft(NZ)*ElmAllocmat4Litr(ielmp,icwood,M,NZ)
    ENDDO D155
    DO NE=1,NumPlantChemElms
      StandDeadStrutElms_pft(NE,NZ)=StandDeadStrutElms_pft(NE,NZ)+&
        sum(StandDeadKCompElms_pft(NE,1:jsken,NZ))
    ENDDO
  ENDIF
  end associate
  end subroutine InitMassBalance
!------------------------------------------------------------------------------------------

  subroutine InitPlantHeatWater(NZ)

  implicit none
  integer, intent(in) :: NZ
  REAL(R8) :: FDM

  associate(                                                &
    ATCA                  => plt_site%ATCA,                 &  !input
    CanOsmoPsi0pt_pft     => plt_ew%CanOsmoPsi0pt_pft,      &  !input
    CanopyLeafShethC_pft  => plt_biom%CanopyLeafShethC_pft, &  !input
    ShootStrutElms_pft    => plt_biom%ShootStrutElms_pft,   &  !input
    CanopyBiomWater_pft   => plt_ew%CanopyBiomWater_pft,    &  !inoput
    HeatCanopy2Dist_col   => plt_ew%HeatCanopy2Dist_col,    &  !inoput
    PSICanopyOsmo_pft     => plt_ew%PSICanopyOsmo_pft,      &  !inoput
    PSICanopy_pft         => plt_ew%PSICanopy_pft,          &  !inoput
    QCanopyWat2Dist_col   => plt_ew%QCanopyWat2Dist_col,    &  !inoput
    TCGroth_pft           => plt_pheno%TCGroth_pft,         &  !inoput
    TKC_pft               => plt_ew%TKC_pft,                &  !inoput
    TdegCCanopy_pft       => plt_ew%TdegCCanopy_pft,        &  !inoput
    VHeatCapCanopy_pft    => plt_ew%VHeatCapCanopy_pft,     &  !inoput
    DeltaTKC_pft          => plt_ew%DeltaTKC_pft,           &  !output
    ENGYX_pft             => plt_ew%ENGYX_pft,              &  !output
    FracPARads2Canopy_pft => plt_rad%FracPARads2Canopy_pft, &  !output
    PSICanopyTurg_pft     => plt_ew%PSICanopyTurg_pft,      &  !output
    TKGroth_pft           => plt_pheno%TKGroth_pft,         &  !output
    Transpiration_pft     => plt_ew%Transpiration_pft,      &  !output
    fTCanopyGroth_pft     => plt_pheno%fTCanopyGroth_pft    &  !output
  )
!
!     INITIALIZE PLANT HEAT AND WATER STATUS
!
!     VHeatCapCanopy_pft=canopy heat capacity (MJ m-3 K-1)
!     TdegCCanopy_pft,TKC=canopy temperature for growth (oC,K)
!     TCGroth_pft,TKGroth_pft=canopy temperature for phenology (oC,K)
!     PSICanopy_pft,PSICanopyOsmo_pft,PSICanopyTurg_pft=canopy total,osmotic,turgor water potl(MPa)
!
  VHeatCapCanopy_pft(NZ)    = cpw*ShootStrutElms_pft(ielmc,NZ)*10.0E-06
  ENGYX_pft(NZ)             = 0._r8
  DeltaTKC_pft(NZ)          = 0._r8
  TdegCCanopy_pft(NZ)       = ATCA
  TKC_pft(NZ)               = units%Celcius2Kelvin(TdegCCanopy_pft(NZ))
  TCGroth_pft(NZ)           = TdegCCanopy_pft(NZ)
  TKGroth_pft(NZ)           = units%Celcius2Kelvin(TCGroth_pft(NZ))
  fTCanopyGroth_pft(NZ)     = 1.0
  PSICanopy_pft(NZ)         = -1.0E-03
  PSICanopyOsmo_pft(NZ)     = CanOsmoPsi0pt_pft(NZ)+PSICanopy_pft(NZ)
  PSICanopyTurg_pft(NZ)     = AZMAX1(PSICanopy_pft(NZ)-PSICanopyOsmo_pft(NZ))
  Transpiration_pft(NZ)     = 0._r8
  FracPARads2Canopy_pft(NZ) = 0._r8
  FDM                       = get_FDM(PSICanopy_pft(NZ))
  CanopyBiomWater_pft(NZ)   = ppmc*CanopyLeafShethC_pft(NZ)/FDM
  VHeatCapCanopy_pft(NZ)    = cpw*(ShootStrutElms_pft(ielmc,NZ)*SpecStalkVolume+CanopyBiomWater_pft(NZ))
  QCanopyWat2Dist_col       = QCanopyWat2Dist_col-CanopyBiomWater_pft(NZ)
  HeatCanopy2Dist_col       = HeatCanopy2Dist_col-VHeatCapCanopy_pft(NZ)*TKC_pft(NZ)

  end associate
  end subroutine InitPlantHeatWater
!------------------------------------------------------------------------------------------

  subroutine InitRootMychorMorphoBio(NZ)
  implicit none
  integer, intent(in) :: NZ
  integer :: K,L,M,N,NR
  REAL(R8) :: CCO2A
  REAL(R8) :: CCO2P
  REAL(R8) :: COXYA
  REAL(R8) :: COXYP
  associate(                                                          &
    ATCA                      => plt_site%ATCA,                       &  !input
    CCO2EI_gperm3                    => plt_site%CCO2EI_gperm3,                     &  !input
    CO2EI                     => plt_site%CO2EI,                      &  !input
    COXYE                     => plt_site%COXYE,                      &  !input
    CanOsmoPsi0pt_pft         => plt_ew%CanOsmoPsi0pt_pft,            &  !input
    NL                        => plt_site%NL,                         &  !input
    OXYE                      => plt_site%OXYE,                       &  !input
    Root1stMaxRadius_pft      => plt_morph%Root1stMaxRadius_pft,      &  !input
    Root2ndMaxRadius_pft      => plt_morph%Root2ndMaxRadius_pft,      &  !input
    RootFracRemobilizableBiom => plt_allom%RootFracRemobilizableBiom, &  !input
    SeedDepth_pft             => plt_morph%SeedDepth_pft,             &  !input
    PSIRootOSMO_vr            => plt_ew%PSIRootOSMO_vr,               &  !inoput
    PSIRoot_pvr               => plt_ew%PSIRoot_pvr,                  &  !inoput
    RootPoreVol_pvr           => plt_morph%RootPoreVol_pvr,           &  !inoput
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr,              &  !inoput
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr,            &  !inoput
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr,            &  !inoput
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft,           &  !output
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr,               &  !output
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr,         &  !output
    Root2ndRadius_pvr         => plt_morph%Root2ndRadius_pvr,         &  !output
    RootN2Fix_pvr             => plt_bgcr%RootN2Fix_pvr,              &  !output
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr,           &  !output
    RootProteinConc_rpvr      => plt_biom%RootProteinConc_rpvr        &  !output
  )
!
!     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
!
!     PSIRoot_pvr,PSIRootOSMO_vr,PSIRootTurg_vr=root,myco total,osmotic,turgor water potl(MPa)
!     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
!     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
!
  NumRootAxes_pft(NZ)=0
  plt_rbgc%RootNH4Uptake_pft(NZ)   = 0._r8
  plt_rbgc%RootNO3Uptake_pft(NZ)   = 0._r8
  plt_rbgc%RootH2PO4Uptake_pft(NZ) = 0._r8
  plt_rbgc%RootHPO4Uptake_pft(NZ)  = 0._r8
  plt_rbgc%RootN2Fix_pft(NZ)       = 0._r8
  DO NR=1,MaxNumRootAxes
    DO N=1,pltpar%jroots
      plt_morph%Root1stDepz_pft(N,NR,NZ)=SeedDepth_pft(NZ)
      plt_biom%RootMyco1stElm_raxs(1:NumPlantChemElms,N,NR,NZ)=0._r8  
    ENDDO
  ENDDO      
  D40: DO N=1,pltpar%jroots
    D20: DO L=1,NL
      plt_ew%AllPlantRootH2OLoss_pvr(N,L,NZ)                        = 0._r8
      PSIRoot_pvr(N,L,NZ)                                          = -0.01_r8
      PSIRootOSMO_vr(N,L,NZ)                                       = CanOsmoPsi0pt_pft(NZ)+PSIRoot_pvr(N,L,NZ)
      PSIRootTurg_vr(N,L,NZ)                                       = AZMAX1(PSIRoot_pvr(N,L,NZ)-PSIRootOSMO_vr(N,L,NZ))
      plt_biom%RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ)   = 0._r8
      plt_biom%RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ) = 0._r8
      RootProteinConc_rpvr(N,L,NZ)                                  = RootFracRemobilizableBiom(NZ)
      plt_biom%RootMycoActiveBiomC_pvr(N,L,NZ)                     = 0._r8
      plt_biom% PopuRootMycoC_pvr(N,L,NZ)=0._r8
      RootProteinC_pvr(N,L,NZ)                                 = 0._r8
      plt_morph%Root1stXNumL_pvr(N,L,NZ)                       = 0._r8
      plt_morph%Root2ndXNum_pvr(N,L,NZ)                        = 0._r8
      plt_morph%RootLenPerPlant_pvr(N,L,NZ)                    = 0._r8
      plt_morph%RootLenDensPerPlant_pvr(N,L,NZ)                = 0._r8
      RootPoreVol_pvr(N,L,NZ)                                  = 0._r8
      RootVH2O_pvr(N,L,NZ)                                     = 0._r8
      Root1stRadius_pvr(N,L,NZ)                                = Root1stMaxRadius_pft(N,NZ)
      Root2ndRadius_pvr(N,L,NZ)                                = Root2ndMaxRadius_pft(N,NZ)
      plt_morph%RootAreaPerPlant_pvr(N,L,NZ)                   = 0._r8
      plt_morph%Root2ndMeanLens_pvr(N,L,NZ)                    = 1.0E-03
      plt_rbgc%RootNutUptake_pvr(ids_NH4B:ids_nuts_end,N,L,NZ) = 0._r8
      plt_rbgc%RootO2Dmnd4Resp_pvr(N,L,NZ)                     = 0._r8
      plt_rbgc%RootNH4DmndSoil_pvr(N,L,NZ)                     = 0._r8
      plt_rbgc%RootNH4DmndBand_pvr(N,L,NZ)                     = 0._r8
      plt_rbgc%RootNO3DmndSoil_pvr(N,L,NZ)                     = 0._r8
      plt_rbgc%RootNO3DmndBand_pvr(N,L,NZ)                     = 0._r8
      plt_rbgc%RootH2PO4DmndSoil_pvr(N,L,NZ)                   = 0._r8
      plt_rbgc%RootH1PO4DmndSoil_pvr(N,L,NZ)                   = 0._r8
      plt_rbgc%RootH2PO4DmndBand_pvr(N,L,NZ)                   = 0._r8
      plt_rbgc%RootH1PO4DmndBand_pvr(N,L,NZ)                   = 0._r8
      plt_rbgc%trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)         = 0._r8
      plt_rbgc%trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)         = 0._r8
      CCO2A                                             = CCO2EI_gperm3
      CCO2P                                             = 0.030*EXP(-2.621_r8-0.0317_r8*ATCA)*CO2EI
      trcg_rootml_pvr(idg_CO2,N,L,NZ)                   = CCO2A*RootPoreVol_pvr(N,L,NZ)
      trcs_rootml_pvr(idg_CO2,N,L,NZ)                   = CCO2P*RootVH2O_pvr(N,L,NZ)
      plt_rbgc%trcg_air2root_flx_pvr(idg_CO2,N,L,NZ)    = 0._r8
      plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_CO2,N,L,NZ) = 0._r8
      plt_rbgc%RootUptkSoiSol_pvr(idg_CO2,N,L,NZ)       = 0._r8
      plt_rbgc%RCO2Emis2Root_pvr(N,L,NZ)                = 0._r8
      COXYA                                             = COXYE
      COXYP                                             = 0.032_r8*EXP(-6.175_r8-0.0211_r8*ATCA)*OXYE
      plt_rbgc%trcg_rootml_pvr(idg_O2,N,L,NZ)=COXYA*RootPoreVol_pvr(N,L,NZ)
      plt_rbgc%trcs_rootml_pvr(idg_O2,N,L,NZ)=COXYP*RootVH2O_pvr(N,L,NZ)
      plt_rbgc%RAutoRootO2Limter_rpvr(N,L,NZ)=1.0
      D30: DO NR=1,MaxNumRootAxes
        plt_morph%Root2ndXNum_rpvr(N,L,NR,NZ)                            = 0._r8
        plt_morph%Root1stLen_rpvr(N,L,NR,NZ)                             = 0._r8
        plt_morph%Root2ndLen_rpvr(N,L,NR,NZ)                              = 0._r8
        plt_biom%RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
        plt_biom%RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ) = 0._r8
      ENDDO D30
      IF(N.EQ.ipltroot)THEN
        D6400: DO K=1,pltpar%NumOfPlantLitrCmplxs
          plt_bgcr%LitrfalStrutElms_pvr(1:NumPlantChemElms,1:jsken,K,L,NZ)=0._r8
        ENDDO D6400
        plt_biom%RootNodulNonstElms_rpvr(1:NumPlantChemElms,L,NZ)=0._r8
        plt_biom%RootNodulStrutElms_rpvr(1:NumPlantChemElms,L,NZ)=0._r8
        RootN2Fix_pvr(L,NZ)=0._r8
      ENDIF
    ENDDO D20
  ENDDO D40

  plt_morph%RootLenDensPerPlant_pvr(1:2,NL+1:JZ1,NZ)=0._r8
  end associate
  end subroutine InitRootMychorMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitSeedMorphoBio(NZ)

  implicit none
  integer, intent(in) :: NZ

  associate(                                                          &
    CNGR_pft                  => plt_allom%CNGR_pft,                  &  !input
    CPGR_pft                  => plt_allom%CPGR_pft,                  &  !input
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,        &  !input
    PetoleStrutElms_brch      => plt_biom%PetoleStrutElms_brch,       &  !input
    PlantPopulation_pft       => plt_site%PlantPopulation_pft,        &  !input
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr,         &  !input
    RootFracRemobilizableBiom => plt_allom%RootFracRemobilizableBiom, &  !input
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr,  &  !input
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr,    &  !input
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,     &  !input
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr,           &  !input
    SeedCMass_pft             => plt_morph%SeedCMass_pft,             &  !input
    CanopyLeafShethC_pft      => plt_biom%CanopyLeafShethC_pft,       &  !inoput
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch,       &  !inoput
    LeafPetolBiomassC_brch    => plt_biom%LeafPetolBiomassC_brch,     &  !inoput
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch,         &  !inoput
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs,        &  !inoput
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft,      &  !inoput
    SeedCPlanted_pft          => plt_biom%SeedCPlanted_pft,           &  !inoput
    WatHeldOnCanopy_pft       => plt_ew%WatHeldOnCanopy_pft           &  !output
  )
!
!     INITIALIZE SEED MORPHOLOGY AND BIOMASS
!
!     WTRVC,WTRVN,WTRVP=C,N,P in storage reserves (g)
!     WTLFB,WTLFBN,WTLFBP=C,N,P in leaves (g)
!     LeafPetolBiomassC_brch=C in leaves+petioles (g)
!     FDM-dry matter fraction (g DM C g FM C-1)
!     CanopyBiomWater_pft,WatHeldOnCanopy_pft=water volume in,on canopy (m3)
!     CPOOL,ZPOOL,PPOOL=C,N,P in canopy nonstructural pools (g)
!     WTRT1,WTRT1N,WTRT1P=C,N,P in primary root layer (g)
!     RTWT1,RTWT1N,RTWT1P=total C,N,P in primary root (g)
!     RootMycoActiveBiomC_pvr, PopuRootMycoC_pvr=total root C mass (g)
!     RootProteinC_pvr=total root protein C mass (g)
!     CPOOLR,ZPOOLR,PPOOLR=C,N,P in root,myco nonstructural pools (g)
!
  SeedCPlanted_pft(NZ)            = SeedCMass_pft(NZ)*PlantPopulation_pft(NZ)
  SeasonalNonstElms_pft(ielmc,NZ) = SeedCPlanted_pft(NZ)
  SeasonalNonstElms_pft(ielmn,NZ) = CNGR_pft(NZ)*SeasonalNonstElms_pft(ielmc,NZ)
  SeasonalNonstElms_pft(ielmp,NZ) = CPGR_pft(NZ)*SeasonalNonstElms_pft(ielmc,NZ)
  LeafStrutElms_brch(ielmn,1,NZ)  = CNGR_pft(NZ)*LeafStrutElms_brch(ielmc,1,NZ)
  LeafStrutElms_brch(ielmp,1,NZ)  = CPGR_pft(NZ)*LeafStrutElms_brch(ielmc,1,NZ)
  LeafPetolBiomassC_brch(1,NZ)    = LeafStrutElms_brch(ielmc,1,NZ)+PetoleStrutElms_brch(ielmc,1,NZ)
  CanopyLeafShethC_pft(NZ)        = CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(1,NZ)

  WatHeldOnCanopy_pft(NZ)         = 0._r8
  CanopyNonstElms_brch(ielmn,1,NZ) = CNGR_pft(NZ)*CanopyNonstElms_brch(ielmc,1,NZ)
  CanopyNonstElms_brch(ielmp,1,NZ) = CPGR_pft(NZ)*CanopyNonstElms_brch(ielmc,1,NZ)
  RootMyco1stStrutElms_rpvr(ielmn,ipltroot,NGTopRootLayer_pft(NZ),1,NZ) = CNGR_pft(NZ) &
    *RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  RootMyco1stStrutElms_rpvr(ielmp,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)=CPGR_pft(NZ) &
    *RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  RootMyco1stElm_raxs(ielmn,1,1,NZ)=CNGR_pft(NZ)*RootMyco1stElm_raxs(ielmc,1,1,NZ)
  RootMyco1stElm_raxs(ielmp,1,1,NZ)=CPGR_pft(NZ)*RootMyco1stElm_raxs(ielmc,1,1,NZ)
  RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=&
    RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  PopuRootMycoC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)= &
    RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  RootProteinC_pvr(1,NGTopRootLayer_pft(NZ),NZ)=RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)&
    *RootFracRemobilizableBiom(NZ)
  RootMycoNonstElms_rpvr(ielmn,1,NGTopRootLayer_pft(NZ),NZ)=CNGR_pft(NZ)&
    *RootMycoNonstElms_rpvr(ielmc,1,NGTopRootLayer_pft(NZ),NZ)
  RootMycoNonstElms_rpvr(ielmp,1,NGTopRootLayer_pft(NZ),NZ)=CPGR_pft(NZ) &
    *RootMycoNonstElms_rpvr(ielmc,1,NGTopRootLayer_pft(NZ),NZ)

  end associate
  end subroutine InitSeedMorphoBio

  end module InitPlantMod
