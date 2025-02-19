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
    PlantPopulation_pft => plt_site%PlantPopulation_pft, &
    NU                  => plt_site%NU,                  &
    NP                  => plt_site%NP,                  &
    NL                  => plt_site%NL,                  &
    ZERO                => plt_site%ZERO,                &
    AREA3               => plt_site%AREA3,               &
    ZERO4Uptk_pft       => plt_rbgc%ZERO4Uptk_pft,       &
    ZERO4Groth_pft      => plt_biom%ZERO4Groth_pft,      &
    ZERO4LeafVar_pft    => plt_biom%ZERO4LeafVar_pft,    &
    IsPlantActive_pft   => plt_pheno%IsPlantActive_pft   &
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
        ZERO4Groth_pft(NZ)=ZERO*PlantPopulation_pft(NZ)
        ZERO4Uptk_pft(NZ)=ZERO*PlantPopulation_pft(NZ)/AREA3(NU)
        ZERO4LeafVar_pft(NZ)=ZERO*PlantPopulation_pft(NZ)*1.0E+06_r8  
      ENDDO  
!
!     FILL OUT UNUSED ARRAYS
!
      D9986: DO NZ=NP+1,JP1
        plt_bgcr%SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ)=0._r8
        plt_bgcr%LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ)=0._r8
        plt_biom%StandDeadStrutElms_pft(1:NumPlantChemElms,NZ)=0._r8
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
    iPlantingDay_pft          => plt_distb%iPlantingDay_pft,          &
    iHarvestYear_pft          => plt_distb%iHarvestYear_pft,          &
    iPlantingYear_pft         => plt_distb%iPlantingYear_pft,         &
    iYearPlanting_pft         => plt_distb%iYearPlanting_pft,         &
    iYearPlantHarvest_pft     => plt_distb%iYearPlantHarvest_pft,     &
    iDayPlantHarvest_pft      => plt_distb%iDayPlantHarvest_pft,      &
    iDayPlanting_pft          => plt_distb%iDayPlanting_pft,          &
    iHarvestDay_pft           => plt_distb%iHarvestDay_pft,           &
    CuticleResist_pft         => plt_photo%CuticleResist_pft,         &
    PPI_pft                   => plt_site%PPI_pft,                    &
    PPX_pft                   => plt_site%PPX_pft,                    &
    PPatSeeding_pft           => plt_site%PPatSeeding_pft,            &
    RootFracRemobilizableBiom => plt_allom%RootFracRemobilizableBiom, &
    rCNNonstRemob_pft         => plt_allom%rCNNonstRemob_pft,         &
    rCPNonstRemob_pft         => plt_allom%rCPNonstRemob_pft,         &
    RootrNC_pft               => plt_allom%RootrNC_pft,               &
    RootrPC_pft               => plt_allom%RootrPC_pft,               &
    O2I                       => plt_photo%O2I,                       &
    CO2CuticleResist_pft      => plt_photo%CO2CuticleResist_pft,      &
    iPlantPhotosynthesisType  => plt_photo%iPlantPhotosynthesisType,  &
    H2OCuticleResist_pft      => plt_photo%H2OCuticleResist_pft,      &
    ClumpFactorInit_pft       => plt_morph%ClumpFactorInit_pft,       &
    ClumpFactor_pft           => plt_morph%ClumpFactor_pft,           &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft            &
  )
  iYearPlanting_pft(NZ)=iPlantingYear_pft(NZ)
  iDayPlanting_pft(NZ)=iPlantingDay_pft(NZ)
  iYearPlantHarvest_pft(NZ)=iHarvestYear_pft(NZ)
  iDayPlantHarvest_pft(NZ)=iHarvestDay_pft(NZ)
  PPI_pft(NZ)=PPatSeeding_pft(NZ)
  PPX_pft(NZ)=PPI_pft(NZ)
  ClumpFactor_pft(NZ)=ClumpFactorInit_pft(NZ)

  H2OCuticleResist_pft(NZ)=CuticleResist_pft(NZ)/3600.0_r8        
  CO2CuticleResist_pft(NZ)=CuticleResist_pft(NZ)*1.56_r8    !1.56=sqrt(44./18.)
  rCNNonstRemob_pft(NZ)=2.5_r8
  rCPNonstRemob_pft(NZ)=25.0_r8
  RootFracRemobilizableBiom(NZ)=AMIN1(RootrNC_pft(NZ)*rCNNonstRemob_pft(NZ),&
    RootrPC_pft(NZ)*rCPNonstRemob_pft(NZ))
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

  associate(                                                          &
    inonstruct                => pltpar%inonstruct,                   &
    ifoliar                   => pltpar%ifoliar,                      &
    inonfoliar                => pltpar%inonfoliar,                   &
    istalk                    => pltpar%istalk,                       &
    iroot                     => pltpar%iroot,                        &
    icwood                    => pltpar%icwood,                       &
    iprotein                  => pltpar%iprotein,                     &
    icarbhyro                 => pltpar%icarbhyro,                    &
    icellulos                 => pltpar%icellulos,                    &
    ilignin                   => pltpar%ilignin,                      &
    NumLitterGroups           => pltpar%NumLitterGroups,              &
    RefLeafAppearRate_pft     => plt_pheno%RefLeafAppearRate_pft,     &
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft, &
    iPlantRootProfile_pft     => plt_pheno%iPlantRootProfile_pft,     &
    MatureGroup_pft           => plt_pheno%MatureGroup_pft,           &
    ElmAllocmat4Litr          => plt_soilchem%ElmAllocmat4Litr,       &
    FracGroth2Node_pft        => plt_allom%FracGroth2Node_pft,        &
    iPlantNfixType_pft        => plt_morph%iPlantNfixType_pft,        &
    NumCogrowthNode_pft        => plt_morph%NumCogrowthNode_pft         &
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

  associate(                                                        &
    DATAP                    => plt_site%DATAP,                     &
    iPlantPhotosynthesisType => plt_photo%iPlantPhotosynthesisType, &
    HighTempLimitSeed_pft    => plt_pheno%HighTempLimitSeed_pft,    &
    TC4LeafOff_pft           => plt_pheno%TC4LeafOff_pft,           &
    TC4LeafOut_pft           => plt_pheno%TC4LeafOut_pft,           &
    TempOffset_pft           => plt_pheno%TempOffset_pft,           &
    PlantInitThermoAdaptZone => plt_pheno%PlantInitThermoAdaptZone, &
    iPlantThermoAdaptZone_pft    => plt_pheno%iPlantThermoAdaptZone_pft,    &
    SeedTempSens_pft         => plt_pheno%SeedTempSens_pft          &
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
  iPlantThermoAdaptZone_pft(NZ)=PlantInitThermoAdaptZone(NZ)
  TempOffset_pft(NZ)=2.667_r8*(2.5_r8-iPlantThermoAdaptZone_pft(NZ))
  TC4LeafOut_pft(NZ)=TCZD-TempOffset_pft(NZ)
  TC4LeafOff_pft(NZ)=AMIN1(15.0_r8,TCXD-TempOffset_pft(NZ))
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
    CNRTS_pft             => plt_allom%CNRTS_pft,             &
    CPRTS_pft             => plt_allom%CPRTS_pft,             &
    RootBiomGrosYld_pft   => plt_allom%RootBiomGrosYld_pft,   &
    RootrNC_pft           => plt_allom%RootrNC_pft,           &
    RootrPC_pft           => plt_allom%RootrPC_pft,           &
    CMinPO4Root_pft       => plt_rbgc%CMinPO4Root_pft,        &
    KmPO4Root_pft         => plt_rbgc%KmPO4Root_pft,          &
    VmaxPO4Root_pft       => plt_rbgc%VmaxPO4Root_pft,        &
    CminNO3Root_pft       => plt_rbgc%CminNO3Root_pft,        &
    KmNO3Root_pft         => plt_rbgc%KmNO3Root_pft,          &
    VmaxNO3Root_pft       => plt_rbgc%VmaxNO3Root_pft,        &
    CMinNH4Root_pft       => plt_rbgc%CMinNH4Root_pft,        &
    VmaxNH4Root_pft       => plt_rbgc%VmaxNH4Root_pft,        &
    KmNH4Root_pft         => plt_rbgc%KmNH4Root_pft,          &
    CumSoilThickness_vr   => plt_site%CumSoilThickness_vr,    &
    NU                    => plt_site%NU,                     &
    NL                    => plt_site%NL,                     &
    NGTopRootLayer_pft    => plt_morph%NGTopRootLayer_pft,    &
    SeedDepth_pft         => plt_morph%SeedDepth_pft,         &
    SeedAreaMean_pft      => plt_morph%SeedAreaMean_pft,      &
    SeedCMass_pft         => plt_morph%SeedCMass_pft,         &
    PlantinDepz_pft       => plt_morph%PlantinDepz_pft,       &
    Root2ndMaxRadius1_pft => plt_morph%Root2ndMaxRadius1_pft, &
    Root2ndXSecArea_pft   => plt_morph%Root2ndXSecArea_pft,   &
    Root1stMaxRadius1_pft => plt_morph%Root1stMaxRadius1_pft, &
    Root1stXSecArea_pft   => plt_morph%Root1stXSecArea_pft,   &
    Root2ndMaxRadius_pft  => plt_morph%Root2ndMaxRadius_pft,  &
    Root1stMaxRadius_pft  => plt_morph%Root1stMaxRadius_pft,  &
    Root2ndSpecLen_pft    => plt_morph%Root2ndSpecLen_pft,    &
    NIXBotRootLayer_pft   => plt_morph%NIXBotRootLayer_pft,   &
    RootRadialResist_pft => plt_morph%RootRadialResist_pft, &
    Root1stSpecLen_pft    => plt_morph%Root1stSpecLen_pft,    &
    RootPorosity_pft      => plt_morph%RootPorosity_pft,      &
    RootPoreTortu4Gas     => plt_morph%RootPoreTortu4Gas,     &
    RootRaidus_rpft       => plt_morph%RootRaidus_rpft,       &
    RootVolPerMassC_pft   => plt_morph%RootVolPerMassC_pft,   &
    RootAxialResist_pft  => plt_morph%RootAxialResist_pft,  &
    NIXBotRootLayer_rpft  => plt_morph%NIXBotRootLayer_rpft,  &
    SeedVolumeMean_pft    => plt_morph%SeedVolumeMean_pft,    &
    SeedMeanLen_pft       => plt_morph%SeedMeanLen_pft        &
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
      NGTopRootLayer_pft(NZ)=L
      NIXBotRootLayer_pft(NZ)=L
      D9790: DO NR=1,pltpar%MaxNumRootAxes
        NIXBotRootLayer_rpft(NR,NZ)=L
      ENDDO D9790
    ENDIF
  ENDDO D9795
  CNRTS_pft(NZ)=RootrNC_pft(NZ)*RootBiomGrosYld_pft(NZ)
  CPRTS_pft(NZ)=RootrPC_pft(NZ)*RootBiomGrosYld_pft(NZ)
  Root1stMaxRadius_pft(2,NZ)=5.0E-06_r8
  Root2ndMaxRadius_pft(2,NZ)=5.0E-06_r8
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
  RootRadialResist_pft(2,NZ)=1.0E+04_r8
  RootAxialResist_pft(2,NZ)=1.0E+12_r8
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
    Root1stSpecLen_pft(N,NZ)=RootVolPerMassC_pft(N,NZ)/(PICON*Root1stMaxRadius_pft(N,NZ)**2._r8)
    Root2ndSpecLen_pft(N,NZ)=RootVolPerMassC_pft(N,NZ)/(PICON*Root2ndMaxRadius_pft(N,NZ)**2._r8)
    Root1stMaxRadius1_pft(N,NZ)=Root1stMaxRadius_pft(N,NZ)
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
    NU                                => plt_site%NU,                                 &
    PPX_pft                           => plt_site%PPX_pft,                            &
    PlantPopulation_pft               => plt_site%PlantPopulation_pft,                &
    ALAT                              => plt_site%ALAT,                               &
    AREA3                             => plt_site%AREA3,                              &
    iPlantRootProfile_pft             => plt_pheno%iPlantRootProfile_pft,             &
    LeafNumberAtFloralInit_brch       => plt_pheno%LeafNumberAtFloralInit_brch,       &
    MatureGroup_brch                  => plt_pheno%MatureGroup_brch,                  &
    KHiestGroLeafNode_brch            => plt_pheno%KHiestGroLeafNode_brch,            &
    KLowestGroLeafNode_brch           => plt_pheno%KLowestGroLeafNode_brch,           &
    NodeNumNormByMatgrp_brch          => plt_pheno%NodeNumNormByMatgrp_brch,          &
    ReprodNodeNumNormByMatrgrp_brch   => plt_pheno%ReprodNodeNumNormByMatrgrp_brch,   &
    HourFailGrainFill_brch            => plt_pheno%HourFailGrainFill_brch,            &
    HoursDoingRemob_brch              => plt_pheno%HoursDoingRemob_brch,              &
    Hours4LenthenPhotoPeriod_brch     => plt_pheno%Hours4LenthenPhotoPeriod_brch,     &
    Hours4LeafOff_brch                => plt_pheno%Hours4LeafOff_brch,                &
    iPlantBranchState_brch            => plt_pheno%iPlantBranchState_brch,            &
    Hours4ShortenPhotoPeriod_brch     => plt_pheno%Hours4ShortenPhotoPeriod_brch,     &
    iPlantCalendar_brch               => plt_pheno%iPlantCalendar_brch,               &
    Hours4Leafout_brch                => plt_pheno%Hours4Leafout_brch,                &
    Hours2LeafOut_brch                => plt_pheno%Hours2LeafOut_brch,                &
    HoursTooLowPsiCan_pft             => plt_pheno%HoursTooLowPsiCan_pft,             &
    MatureGroup_pft                   => plt_pheno%MatureGroup_pft,                   &
    PetioleChemElmRemobFlx_brch       => plt_pheno%PetioleChemElmRemobFlx_brch,       &
    TotReproNodeNumNormByMatrgrp_brch => plt_pheno%TotReproNodeNumNormByMatrgrp_brch, &
    TotalNodeNumNormByMatgrp_brch     => plt_pheno%TotalNodeNumNormByMatgrp_brch,     &
    C4PhotosynDowreg_brch             => plt_photo%C4PhotosynDowreg_brch,             &
    CPOOL3_node                       => plt_photo%CPOOL3_node,                       &
    CPOOL4_node                       => plt_photo%CPOOL4_node,                       &
    ChillHours_pft                    => plt_photo%ChillHours_pft,                    &
    RubiscoActivity_brch              => plt_photo%RubiscoActivity_brch,              &
    CMassHCO3BundleSheath_node        => plt_photo%CMassHCO3BundleSheath_node,        &
    CMassCO2BundleSheath_node         => plt_photo%CMassCO2BundleSheath_node,         &
    NumOfLeaves_brch                  => plt_morph%NumOfLeaves_brch,                  &
    CanopyHeight_pft                  => plt_morph%CanopyHeight_pft,                  &
    KLeafNumber_brch                  => plt_morph%KLeafNumber_brch,                  &
    ShootNodeNumAtPlanting_pft        => plt_morph%ShootNodeNumAtPlanting_pft,        &
    BranchNumber_pft                  => plt_morph%BranchNumber_pft,                  &
    ShootNodeNum_brch                 => plt_morph%ShootNodeNum_brch,                 &
    CanopyLeafArea_pft                => plt_morph%CanopyLeafArea_pft,                &
    CanopyStalkArea_lbrch             => plt_morph%CanopyStalkArea_lbrch,             &
    StemAreaZsec_brch                 => plt_morph%StemAreaZsec_brch,                 &
    CanopyStemArea_pft                => plt_morph%CanopyStemArea_pft,                &
    NodeNumberAtAnthesis_brch         => plt_morph%NodeNumberAtAnthesis_brch,         &
    PotentialSeedSites_brch           => plt_morph%PotentialSeedSites_brch,           &
    LiveInterNodeHight_brch           => plt_morph%LiveInterNodeHight_brch,           &
    SeedNumSet_brch                   => plt_morph%SeedNumSet_brch,                   &
    LeafAreaLive_brch                 => plt_morph%LeafAreaLive_brch,                 &
    LeafAreaDying_brch                => plt_morph%LeafAreaDying_brch,                &
    CanopyLeafArea_lpft               => plt_morph%CanopyLeafArea_lpft,               &
    CanopyLeafAreaZ_pft               => plt_morph%CanopyLeafAreaZ_pft,               &
    CanopyStemAreaZ_pft               => plt_morph%CanopyStemAreaZ_pft,               &
    LeafAreaZsec_brch                 => plt_morph%LeafAreaZsec_brch,                 &
    LeafNodeArea_brch                 => plt_morph%LeafNodeArea_brch,                 &
    CanPBranchHeight                  => plt_morph%CanPBranchHeight,                  &
    HypoctoHeight_pft                 => plt_morph%HypoctoHeight_pft,                 &
    BranchNumber_brch                 => plt_morph%BranchNumber_brch,                 &
    NodeNum2InitFloral_brch           => plt_morph%NodeNum2InitFloral_brch,           &
    KMinNumLeaf4GroAlloc_brch         => plt_morph%KMinNumLeaf4GroAlloc_brch,         &
    NumOfBranches_pft                 => plt_morph%NumOfBranches_pft                  &
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
    plt_biom%StalkLiveBiomassC_brch(NB,NZ)       = 0._r8
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
      plt_morph%InternodeHeightDead_brch(K,NB,NZ)                 = 0._r8
      plt_morph%PetoleLensNode_brch(K,NB,NZ)                       = 0._r8
      plt_biom%LeafProteinCNode_brch(K,NB,NZ)                      = 0._r8
      plt_biom%PetoleProteinCNode_brch(K,NB,NZ)                    = 0._r8
      plt_biom%LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)      = 0._r8
      plt_biom%InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ) = 0._r8
      plt_biom%PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ)   = 0._r8

      D55: DO L=1,NumOfCanopyLayers1
        CanopyLeafArea_lpft(L,K,NB,NZ)=0._r8
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
    NU                             => plt_site%NU,                             &
    AREA3                          => plt_site%AREA3,                          &
    ElmAllocmat4Litr               => plt_soilchem%ElmAllocmat4Litr,           &
    ETCanopy_CumYr_pft             => plt_ew%ETCanopy_CumYr_pft,               &
    StandDeadKCompElms_pft         => plt_biom%StandDeadKCompElms_pft,         &
    StandDeadStrutElms_pft         => plt_biom%StandDeadStrutElms_pft,         &
    StandingDeadInitC_pft          => plt_biom%StandingDeadInitC_pft,          &
    RootBiomCPerPlant_pft          => plt_biom%RootBiomCPerPlant_pft,          &
    rNCStalk_pft                   => plt_allom%rNCStalk_pft,                  &
    rPCStalk_pft                   => plt_allom%rPCStalk_pft,                  &
    PlantExudElm_CumYr_pft         => plt_rbgc%PlantExudElm_CumYr_pft,         &
    RootUptk_N_CumYr_pft           => plt_rbgc%RootUptk_N_CumYr_pft,           &
    RootUptk_P_CumYr_pft           => plt_rbgc%RootUptk_P_CumYr_pft,           &
    NH3Emis_CumYr_pft              => plt_bgcr%NH3Emis_CumYr_pft,              &
    NH3Dep2Can_pft                 => plt_bgcr%NH3Dep2Can_pft,                 &
    SurfLitrfalStrutElms_CumYr_pft => plt_bgcr%SurfLitrfalStrutElms_CumYr_pft, &
    GrossCO2Fix_pft                => plt_bgcr%GrossCO2Fix_pft,                &
    CanopyRespC_CumYr_pft          => plt_bgcr%CanopyRespC_CumYr_pft,          &
    GrossResp_pft                  => plt_bgcr%GrossResp_pft,                  &
    PlantN2Fix_CumYr_pft           => plt_bgcr%PlantN2Fix_CumYr_pft,           &
    LitrfalStrutElms_CumYr_pft     => plt_bgcr%LitrfalStrutElms_CumYr_pft,     &
    CanopyStemArea_pft             => plt_morph%CanopyStemArea_pft,            &
    icwood                         => pltpar%icwood,                           &
    NetCumElmntFlx2Plant_pft       => plt_pheno%NetCumElmntFlx2Plant_pft       &
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
    ATCA                  => plt_site%ATCA,                 &
    CanOsmoPsi0pt_pft     => plt_ew%CanOsmoPsi0pt_pft,      &
    HeatCanopy2Dist_col   => plt_ew%HeatCanopy2Dist_col,    &
    TKC_pft               => plt_ew%TKC_pft,                &
    Transpiration_pft     => plt_ew%Transpiration_pft,      &
    VHeatCapCanopy_pft    => plt_ew%VHeatCapCanopy_pft,     &
    PSICanopy_pft         => plt_ew%PSICanopy_pft,          &
    PSICanopyTurg_pft     => plt_ew%PSICanopyTurg_pft,      &
    PSICanopyOsmo_pft     => plt_ew%PSICanopyOsmo_pft,      &
    ENGYX_pft             => plt_ew%ENGYX_pft,              &
    DeltaTKC_pft          => plt_ew%DeltaTKC_pft,           &
    TdegCCanopy_pft    => plt_ew%TdegCCanopy_pft,     &
    TKGroth_pft           => plt_pheno%TKGroth_pft,         &
    TCGroth_pft           => plt_pheno%TCGroth_pft,         &
    CanopyBiomWater_pft       => plt_ew%CanopyBiomWater_pft,        &
    QCanopyWat2Dist_col   => plt_ew%QCanopyWat2Dist_col,    &
    fTCanopyGroth_pft     => plt_pheno%fTCanopyGroth_pft,   &
    CanopyLeafShethC_pft  => plt_biom%CanopyLeafShethC_pft, &
    ShootStrutElms_pft    => plt_biom%ShootStrutElms_pft,   &
    FracPARads2Canopy_pft => plt_rad%FracPARads2Canopy_pft  &
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
  TdegCCanopy_pft(NZ)    = ATCA
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
  CanopyBiomWater_pft(NZ)       = ppmc*CanopyLeafShethC_pft(NZ)/FDM
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
    CanOsmoPsi0pt_pft         => plt_ew%CanOsmoPsi0pt_pft,            &
    PSICanopy_pft             => plt_ew%PSICanopy_pft,                &
    PSIRoot_pvr               => plt_ew%PSIRoot_pvr,                  &
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr,               &
    PSIRootOSMO_vr            => plt_ew%PSIRootOSMO_vr,               &
    RootN2Fix_pvr             => plt_bgcr%RootN2Fix_pvr,              &
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr,            &
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr,            &
    OXYE                      => plt_site%OXYE,                       &
    COXYE                     => plt_site%COXYE,                      &
    CO2EI                     => plt_site%CO2EI,                      &
    NL                        => plt_site%NL,                         &
    ATCA                      => plt_site%ATCA,                       &
    CCO2EI                    => plt_site%CCO2EI,                     &
    RootFracRemobilizableBiom => plt_allom%RootFracRemobilizableBiom, &
    RootProteinConc_rpvr       => plt_biom%RootProteinConc_rpvr,        &
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr,           &
    SeedDepth_pft             => plt_morph%SeedDepth_pft,             &
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr,              &
    RootPoreVol_pvr           => plt_morph%RootPoreVol_pvr,           &
    Root2ndRadius_pvr         => plt_morph%Root2ndRadius_pvr,         &
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr,         &
    Root1stMaxRadius_pft      => plt_morph%Root1stMaxRadius_pft,      &
    Root2ndMaxRadius_pft      => plt_morph%Root2ndMaxRadius_pft,      &
    NumRootAxes_pft           => plt_morph%NumRootAxes_pft            &
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
      plt_ew%AllPlantRootH2OLoss_vr(N,L,NZ)                      = 0._r8
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
      plt_morph%Root2ndAveLen_pvr(N,L,NZ)                      = 1.0E-03
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
      plt_rbgc%trcg_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)       = 0._r8
      plt_rbgc%trcs_rootml_pvr(idg_beg:idg_NH3,N,L,NZ)       = 0._r8
      CCO2A                                                    = CCO2EI
      CCO2P                                             = 0.030*EXP(-2.621_r8-0.0317_r8*ATCA)*CO2EI
      trcg_rootml_pvr(idg_CO2,N,L,NZ)                   = CCO2A*RootPoreVol_pvr(N,L,NZ)
      trcs_rootml_pvr(idg_CO2,N,L,NZ)                   = CCO2P*RootVH2O_pvr(N,L,NZ)
      plt_rbgc%trcg_air2root_flx_pvr(idg_CO2,N,L,NZ)   = 0._r8
      plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_CO2,N,L,NZ) = 0._r8
      plt_rbgc%RootUptkSoiSol_vr(idg_CO2,N,L,NZ)             = 0._r8
      plt_rbgc%RootCO2Emis_pvr(N,L,NZ)                  = 0._r8
      COXYA                                             = COXYE
      COXYP=0.032_r8*EXP(-6.175_r8-0.0211_r8*ATCA)*OXYE
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
    PlantPopulation_pft       => plt_site%PlantPopulation_pft,        &
    PSICanopy_pft             => plt_ew%PSICanopy_pft,                &
    WatHeldOnCanopy_pft       => plt_ew%WatHeldOnCanopy_pft,          &
    RootFracRemobilizableBiom => plt_allom%RootFracRemobilizableBiom, &
    CNGR                      => plt_allom%CNGR,                      &
    CPGR                      => plt_allom%CPGR,                      &
    RootMyco1stElm_raxs       => plt_biom%RootMyco1stElm_raxs,        &
    SeedCPlanted_pft          => plt_biom%SeedCPlanted_pft,           &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr,  &
    LeafPetolBiomassC_brch    => plt_biom%LeafPetolBiomassC_brch,     &
    CanopyLeafShethC_pft      => plt_biom%CanopyLeafShethC_pft,       &
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr,         &
    PetoleStrutElms_brch      => plt_biom%PetoleStrutElms_brch,       &
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr,    &
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr,           &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,     &
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch,       &
    LeafStrutElms_brch        => plt_biom%LeafStrutElms_brch,         &
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft,      &
    SeedCMass_pft             => plt_morph%SeedCMass_pft,             &
    Root1stDepz_pft           => plt_morph%Root1stDepz_pft,           &
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft         &
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
  SeasonalNonstElms_pft(ielmn,NZ) = CNGR(NZ)*SeasonalNonstElms_pft(ielmc,NZ)
  SeasonalNonstElms_pft(ielmp,NZ) = CPGR(NZ)*SeasonalNonstElms_pft(ielmc,NZ)
  LeafStrutElms_brch(ielmn,1,NZ)  = CNGR(NZ)*LeafStrutElms_brch(ielmc,1,NZ)
  LeafStrutElms_brch(ielmp,1,NZ)  = CPGR(NZ)*LeafStrutElms_brch(ielmc,1,NZ)
  LeafPetolBiomassC_brch(1,NZ)    = LeafStrutElms_brch(ielmc,1,NZ)+PetoleStrutElms_brch(ielmc,1,NZ)
  CanopyLeafShethC_pft(NZ)        = CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(1,NZ)

  WatHeldOnCanopy_pft(NZ)         = 0._r8
  CanopyNonstElms_brch(ielmn,1,NZ) = CNGR(NZ)*CanopyNonstElms_brch(ielmc,1,NZ)
  CanopyNonstElms_brch(ielmp,1,NZ) = CPGR(NZ)*CanopyNonstElms_brch(ielmc,1,NZ)
  RootMyco1stStrutElms_rpvr(ielmn,ipltroot,NGTopRootLayer_pft(NZ),1,NZ) = CNGR(NZ) &
    *RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  RootMyco1stStrutElms_rpvr(ielmp,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)=CPGR(NZ) &
    *RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  RootMyco1stElm_raxs(ielmn,1,1,NZ)=CNGR(NZ)*RootMyco1stElm_raxs(ielmc,1,1,NZ)
  RootMyco1stElm_raxs(ielmp,1,1,NZ)=CPGR(NZ)*RootMyco1stElm_raxs(ielmc,1,1,NZ)
  RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)=&
    RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  PopuRootMycoC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)= &
    RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ),1,NZ)
  RootProteinC_pvr(1,NGTopRootLayer_pft(NZ),NZ)=RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ),NZ)&
    *RootFracRemobilizableBiom(NZ)
  RootMycoNonstElms_rpvr(ielmn,1,NGTopRootLayer_pft(NZ),NZ)=CNGR(NZ)&
    *RootMycoNonstElms_rpvr(ielmc,1,NGTopRootLayer_pft(NZ),NZ)
  RootMycoNonstElms_rpvr(ielmp,1,NGTopRootLayer_pft(NZ),NZ)=CPGR(NZ) &
    *RootMycoNonstElms_rpvr(ielmc,1,NGTopRootLayer_pft(NZ),NZ)

  end associate
  end subroutine InitSeedMorphoBio

  end module InitPlantMod
