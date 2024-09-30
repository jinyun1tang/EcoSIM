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
  use PlantMgmtDataType
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
  character(len=*),private, parameter :: mod_filename =&
   __FILE__
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
!     IsPlantActive_pft=PFT flag:0=not active,1=active
!     iYearPlanting_pft,iDayPlanting_pft,iYearPlantHarvest_pft,iDayPlantHarvest_pft=year,day of planting,arvesting
!     PPI,PPX=initial,current population (m-2)
!     CF,ClumpFactorInit_pft=current,initial clumping factor
!     H2OCuticleResist_pft=cuticular resistance to water (h m-1)
!     CO2CuticleResist_pft=cuticular resistance to CO2 (s m-1)
!     CNWS,rCPNonstRemob_pft=protein:N,protein:P ratios
!    RootFracRemobilizableBiom=maximum root protein concentration (g g-1)
!     O2I=intercellular O2 concentration in C3,C4 PFT (umol mol-1)
!

  D9995: DO NX=NHWQ,NHEQ
    D9990: DO NY=NVNQ,NVSQ
      NZ2X=MIN(NZ2Q,NP(NY,NX))
      D9985: DO NZ=NZ1Q,NZ2X
        IF(IsPlantActive_pft(NZ,NY,NX).EQ.iDormant)THEN

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
        ZERO4Groth_pft(NZ,NY,NX)   = ZERO*PlantPopulation_pft(NZ,NY,NX)
        ZERO4Uptk_pft(NZ,NY,NX)    = ZERO*PlantPopulation_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        ZERO4LeafVar_pft(NZ,NY,NX) = ZERO*PlantPopulation_pft(NZ,NY,NX)*1.0E+06_r8
      ENDDO D9985
!
!     FILL OUT UNUSED ARRAYS
!
      D9986: DO NZ=NP(NY,NX)+1,JP
        SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX) = 0._r8
        LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)     = 0._r8
        StandDeadStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)      = 0._r8
        D6401: DO L=1,NL(NY,NX)
          DO  K=1,pltpar%NumOfPlantLitrCmplxs
            DO  M=1,jskenc
              LitrfalStrutElms_pvr(1:NumPlantChemElms,M,K,L,NZ,NY,NX)=0._r8
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

  iYearPlanting_pft(NZ,NY,NX)=iPlantingYear_pft(NZ,NY,NX)   !planting year
  iDayPlanting_pft(NZ,NY,NX)=iPlantingDay_pft(NZ,NY,NX) !planting day
  iYearPlantHarvest_pft(NZ,NY,NX)=iHarvestYear_pft(NZ,NY,NX)
  iDayPlantHarvest_pft(NZ,NY,NX)=iHarvestDay_pft(NZ,NY,NX)
  PPI_pft(NZ,NY,NX)=PPatSeeding_pft(NZ,NY,NX)
  PPX_pft(NZ,NY,NX)=PPI_pft(NZ,NY,NX)
  ClumpFactor_pft(NZ,NY,NX)=ClumpFactorInit_pft(NZ,NY,NX)       !clumping factor
  
  H2OCuticleResist_pft(NZ,NY,NX)=CuticleResist_pft(NZ,NY,NX)/3600.0_r8
  CO2CuticleResist_pft(NZ,NY,NX)=CuticleResist_pft(NZ,NY,NX)*1.56_r8
  rCNNonstRemob_pft(NZ,NY,NX)=2.5_r8
  rCPNonstRemob_pft(NZ,NY,NX)=25.0_r8
  RootFracRemobilizableBiom(NZ,NY,NX)=AMIN1(RootrNC_pft(NZ,NY,NX)*rCNNonstRemob_pft(NZ,NY,NX)&
    ,RootrPC_pft(NZ,NY,NX)*rCPNonstRemob_pft(NZ,NY,NX))
  IF(iPlantPhotosynthesisType(NZ,NY,NX).EQ.ic3_photo)THEN
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
    iprotein         => pltpar%iprotein  ,&
    icarbhyro        => pltpar%icarbhyro ,&
    icellulos        => pltpar%icellulos ,&
    ilignin          =>  pltpar%ilignin   , &
    NumLitterGroups  => pltpar%NumLitterGroups , &
    inonstruct       => pltpar%inonstruct, &
    ifoliar          => pltpar%ifoliar , &
    inonfoliar       => pltpar%inonfoliar, &
    istalk           => pltpar%istalk  , &
    iroot            => pltpar%iroot   , &
    icwood           => pltpar%icwood    &
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
  ElmAllocMat4Litr(ielmc,inonstruct,iprotein,NZ,NY,NX)=0.0_r8
  ElmAllocMat4Litr(ielmc,inonstruct,icarbhyro,NZ,NY,NX)=0.67_r8
  ElmAllocMat4Litr(ielmc,inonstruct,icellulos,NZ,NY,NX)=0.33_r8
  ElmAllocMat4Litr(ielmc,inonstruct,ilignin,NZ,NY,NX)=0.0_r8
!
!     NON-VASCULAR (E.G. MOSSES)
!
  IF(is_root_shallow(iPlantRootProfile_pft(NZ,NY,NX)))THEN
    ElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ,NY,NX)=0.07_r8
    ElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ,NY,NX)=0.25_r8
    ElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ,NY,NX)=0.30_r8
    ElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ,NY,NX)=0.38_r8

    ElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ,NY,NX)=0.07_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ,NY,NX)=0.25_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ,NY,NX)=0.30_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ,NY,NX)=0.38_r8
!
!     LEGUMES
!
  ELSEIF(iPlantNfixType_pft(NZ,NY,NX).NE.0)THEN
    ElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ,NY,NX)=0.16_r8
    ElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ,NY,NX)=0.38_r8
    ElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ,NY,NX)=0.34_r8
    ElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ,NY,NX)=0.12_r8

    ElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ,NY,NX)=0.07_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ,NY,NX)=0.41_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ,NY,NX)=0.37_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ,NY,NX)=0.15_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.0 .OR. &
    (.not.is_plant_treelike(iPlantRootProfile_pft(NZ,NY,NX))))THEN
    ElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ,NY,NX)=0.08_r8
    ElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ,NY,NX)=0.41_r8
    ElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ,NY,NX)=0.36_r8
    ElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ,NY,NX)=0.15_r8

    ElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ,NY,NX)=0.07_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ,NY,NX)=0.41_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ,NY,NX)=0.36_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ,NY,NX)=0.16_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.1 .OR. iPlantTurnoverPattern_pft(NZ,NY,NX).GE.3)THEN
    ElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ,NY,NX)=0.07_r8
    ElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ,NY,NX)=0.34_r8
    ElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ,NY,NX)=0.36_r8
    ElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ,NY,NX)=0.23_r8

    ElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ,NY,NX)=0.0_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ,NY,NX)=0.045_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ,NY,NX)=0.660_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ,NY,NX)=0.295_r8
!
!     CONIFEROUS TREES
!
  ELSE
    ElmAllocMat4Litr(ielmc,ifoliar,iprotein,NZ,NY,NX)=0.07_r8
    ElmAllocMat4Litr(ielmc,ifoliar,icarbhyro,NZ,NY,NX)=0.25_r8
    ElmAllocMat4Litr(ielmc,ifoliar,icellulos,NZ,NY,NX)=0.38_r8
    ElmAllocMat4Litr(ielmc,ifoliar,ilignin,NZ,NY,NX)=0.30_r8

    ElmAllocMat4Litr(ielmc,inonfoliar,iprotein,NZ,NY,NX)=0.0_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,icarbhyro,NZ,NY,NX)=0.045_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,icellulos,NZ,NY,NX)=0.660_r8
    ElmAllocMat4Litr(ielmc,inonfoliar,ilignin,NZ,NY,NX)=0.295_r8
  ENDIF
!
!     FRACTIONS OF WOODY LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     NON-VASCULAR
!
  IF(is_root_shallow(iPlantRootProfile_pft(NZ,NY,NX)))THEN
    ElmAllocMat4Litr(ielmc,istalk,iprotein,NZ,NY,NX)=0.07_r8
    ElmAllocMat4Litr(ielmc,istalk,icarbhyro,NZ,NY,NX)=0.25_r8
    ElmAllocMat4Litr(ielmc,istalk,icellulos,NZ,NY,NX)=0.30_r8
    ElmAllocMat4Litr(ielmc,istalk,ilignin,NZ,NY,NX)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.0 .OR. &
    (.not.is_plant_treelike(iPlantRootProfile_pft(NZ,NY,NX))))THEN
    ElmAllocMat4Litr(ielmc,istalk,iprotein,NZ,NY,NX)=0.03_r8
    ElmAllocMat4Litr(ielmc,istalk,icarbhyro,NZ,NY,NX)=0.25_r8
    ElmAllocMat4Litr(ielmc,istalk,icellulos,NZ,NY,NX)=0.57_r8
    ElmAllocMat4Litr(ielmc,istalk,ilignin,NZ,NY,NX)=0.15_r8
!
!     DECIDUOUS AND CONIFEROUS TREES
!
  ELSE
    ElmAllocMat4Litr(ielmc,istalk,iprotein,NZ,NY,NX)=0.0_r8
    ElmAllocMat4Litr(ielmc,istalk,icarbhyro,NZ,NY,NX)=0.045_r8
    ElmAllocMat4Litr(ielmc,istalk,icellulos,NZ,NY,NX)=0.660_r8
    ElmAllocMat4Litr(ielmc,istalk,ilignin,NZ,NY,NX)=0.295_r8
  ENDIF
!
!     FRACTIONS OF FINE ROOT LITTER ALLOCATED TO
!     PROTEIN, CH2O, CELLULOSE, LIGNIN PC&E 25:601-608
!
!     NON-VASCULAR
!
  IF(is_root_shallow(iPlantRootProfile_pft(NZ,NY,NX)))THEN
    ElmAllocMat4Litr(ielmc,iroot,iprotein,NZ,NY,NX)=0.07_r8
    ElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ,NY,NX)=0.25_r8
    ElmAllocMat4Litr(ielmc,iroot,icellulos,NZ,NY,NX)=0.30_r8
    ElmAllocMat4Litr(ielmc,iroot,ilignin,NZ,NY,NX)=0.38_r8
!
!     ANNUALS, GRASSES, SHRUBS
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.0 .OR. &
    (.not.is_plant_treelike(iPlantRootProfile_pft(NZ,NY,NX))))THEN
    ElmAllocMat4Litr(ielmc,iroot,iprotein,NZ,NY,NX)=0.057_r8
    ElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ,NY,NX)=0.263_r8
    ElmAllocMat4Litr(ielmc,iroot,icellulos,NZ,NY,NX)=0.542_r8
    ElmAllocMat4Litr(ielmc,iroot,ilignin,NZ,NY,NX)=0.138_r8
!
!     DECIDUOUS TREES
!
  ELSEIF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.1 .OR. iPlantTurnoverPattern_pft(NZ,NY,NX).GE.3)THEN
    ElmAllocMat4Litr(ielmc,iroot,iprotein,NZ,NY,NX)=0.059_r8
    ElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ,NY,NX)=0.308_r8
    ElmAllocMat4Litr(ielmc,iroot,icellulos,NZ,NY,NX)=0.464_r8
    ElmAllocMat4Litr(ielmc,iroot,ilignin,NZ,NY,NX)=0.169_r8
!
!     CONIFEROUS TREES
!
  ELSE
    ElmAllocMat4Litr(ielmc,iroot,iprotein,NZ,NY,NX)=0.059_r8
    ElmAllocMat4Litr(ielmc,iroot,icarbhyro,NZ,NY,NX)=0.308_r8
    ElmAllocMat4Litr(ielmc,iroot,icellulos,NZ,NY,NX)=0.464_r8
    ElmAllocMat4Litr(ielmc,iroot,ilignin,NZ,NY,NX)=0.169_r8
  ENDIF
!
!     COARSE WOODY LITTER FROM BOLES AND ROOTS
!
  ElmAllocMat4Litr(ielmc,icwood,iprotein,NZ,NY,NX)=0.00_r8
  ElmAllocMat4Litr(ielmc,icwood,icarbhyro,NZ,NY,NX)=0.045_r8
  ElmAllocMat4Litr(ielmc,icwood,icellulos,NZ,NY,NX)=0.660_r8
  ElmAllocMat4Litr(ielmc,icwood,ilignin,NZ,NY,NX)=0.295_r8
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
    D100: DO M=1,jskenc
      CNOPCT=CNOPCT+ElmAllocMat4Litr(ielmc,N,M,NZ,NY,NX)*CNOPC(M)
      CPOPCT=CPOPCT+ElmAllocMat4Litr(ielmc,N,M,NZ,NY,NX)*CPOPC(M)
    ENDDO D100
    D105: DO M=1,jskenc
      ElmAllocMat4Litr(ielmn,N,M,NZ,NY,NX)=ElmAllocMat4Litr(ielmc,N,M,NZ,NY,NX)*CNOPC(M)/CNOPCT
      ElmAllocMat4Litr(ielmp,N,M,NZ,NY,NX)=ElmAllocMat4Litr(ielmc,N,M,NZ,NY,NX)*CPOPC(M)/CPOPCT
    ENDDO D105
  ENDDO D110
!
!     CONCURRENT NODE GROWTH
!
!     FracGroth2Node_pft=scales node number for perennial vegetation (e.g. trees)
!     NumCogrothNode_pft=number of concurrently growing nodes
!
  IF(iPlantTurnoverPattern_pft(NZ,NY,NX).EQ.0 .OR. &
    (.not.is_plant_treelike(iPlantRootProfile_pft(NZ,NY,NX))))THEN
    ! deciduous or shallow root
    FracGroth2Node_pft(NZ,NY,NX)=1.0_r8
!
    IF(MatureGroup_pft(NZ,NY,NX).LE.10)THEN
      NumCogrothNode_pft(NZ,NY,NX)=3
    ELSEIF(MatureGroup_pft(NZ,NY,NX).LE.15)THEN
      NumCogrothNode_pft(NZ,NY,NX)=4
    ELSE
      NumCogrothNode_pft(NZ,NY,NX)=5
    ENDIF
  ELSE
    FracGroth2Node_pft(NZ,NY,NX)=AMAX1(1.0_r8,0.04_r8/RefLeafAppearRate_pft(NZ,NY,NX))
    NumCogrothNode_pft(NZ,NY,NX)=24
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
!     ZTYP,PlantInitThermoAdaptZone=dynamic,initial thermal adaptation zone from PFT file
!     TempOffset_pft=shift in Arrhenius curve for thermal adaptation (oC)
!     TCZ,TC4LeafOff_pft=threshold temperature for leafout,leafoff
!     HTC=high temperature threshold for grain number loss (oC)
!     SeedTempSens_pft=sensitivity to HTC (seeds oC-1 above HTC)
!
  iPlantThermoAdaptZone_pft(NZ,NY,NX)=PlantInitThermoAdaptZone(NZ,NY,NX)
  TempOffset_pft(NZ,NY,NX)=2.667*(2.5-iPlantThermoAdaptZone_pft(NZ,NY,NX))
  TC4LeafOut_pft(NZ,NY,NX)=TCZD-TempOffset_pft(NZ,NY,NX)
  TC4LeafOff_pft(NZ,NY,NX)=AMIN1(15.0,TCXD-TempOffset_pft(NZ,NY,NX))
  IF(iPlantPhotosynthesisType(NZ,NY,NX).EQ.3)THEN
    IF(DATAP(NZ,NY,NX)(1:4).EQ.'soyb')THEN
      HighTempLimitSeed_pft(NZ,NY,NX)=30.0_r8+3.0_r8*iPlantThermoAdaptZone_pft(NZ,NY,NX)
      SeedTempSens_pft(NZ,NY,NX)=0.002_r8
    ELSE
      HighTempLimitSeed_pft(NZ,NY,NX)=27.0_r8+3.0_r8*iPlantThermoAdaptZone_pft(NZ,NY,NX)
      SeedTempSens_pft(NZ,NY,NX)=0.002_r8
    ENDIF
  ELSE
    HighTempLimitSeed_pft(NZ,NY,NX)=27.0_r8+3.0_r8*iPlantThermoAdaptZone_pft(NZ,NY,NX)
    SeedTempSens_pft(NZ,NY,NX)=0.005_r8
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
  call calc_seed_geometry(SeedCMass_pft(NZ,NY,NX),SeedVolumeMean_pft(NZ,NY,NX),&
    SeedMeanLen_pft(NZ,NY,NX),SeedAreaMean_pft(NZ,NY,NX))

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
!     RSRR,RSRA=radial,axial root resistivity (m2 MPa-1 h-1)
!
  SeedDepth_pft(NZ,NY,NX)=PlantinDepz_pft(NZ,NY,NX)
  D9795: DO L=NU(NY,NX),NL(NY,NX)
    IF(SeedDepth_pft(NZ,NY,NX).GE.CumSoilThickness_vr(L-1,NY,NX) &
      .AND.SeedDepth_pft(NZ,NY,NX).LT.CumSoilThickness_vr(L,NY,NX))THEN
      !find the seeding layer
      NGTopRootLayer_pft(NZ,NY,NX)=L
      NIXBotRootLayer_pft(NZ,NY,NX)=L
      D9790: DO NR=1,pltpar%MaxNumRootAxes
        NIXBotRootLayer_rpft(NR,NZ,NY,NX)=L
      ENDDO D9790
    ENDIF
  ENDDO D9795  
  CNRTS_pft(NZ,NY,NX)=RootrNC_pft(NZ,NY,NX)*RootBiomGrosYld_pft(NZ,NY,NX)
  CPRTS_pft(NZ,NY,NX)=RootrPC_pft(NZ,NY,NX)*RootBiomGrosYld_pft(NZ,NY,NX)
  Root1stMaxRadius_pft(imycorrhz,NZ,NY,NX)=5.0E-06
  Root2ndMaxRadius_pft(imycorrhz,NZ,NY,NX)=5.0E-06
  RootPorosity_pft(imycorrhz,NZ,NY,NX)=RootPorosity_pft(1,NZ,NY,NX)
  VmaxNH4Root_pft(imycorrhz,NZ,NY,NX)=VmaxNH4Root_pft(1,NZ,NY,NX)
  KmNH4Root_pft(imycorrhz,NZ,NY,NX)=KmNH4Root_pft(1,NZ,NY,NX)
  CMinNH4Root_pft(imycorrhz,NZ,NY,NX)=CMinNH4Root_pft(1,NZ,NY,NX)
  VmaxNO3Root_pft(imycorrhz,NZ,NY,NX)=VmaxNO3Root_pft(1,NZ,NY,NX)
  KmNO3Root_pft(imycorrhz,NZ,NY,NX)=KmNO3Root_pft(1,NZ,NY,NX)
  CminNO3Root_pft(imycorrhz,NZ,NY,NX)=CminNO3Root_pft(1,NZ,NY,NX)
  VmaxPO4Root_pft(imycorrhz,NZ,NY,NX)=VmaxPO4Root_pft(1,NZ,NY,NX)
  KmPO4Root_pft(imycorrhz,NZ,NY,NX)=KmPO4Root_pft(1,NZ,NY,NX)
  CMinPO4Root_pft(imycorrhz,NZ,NY,NX)=CMinPO4Root_pft(1,NZ,NY,NX)
  RootRadialResist_pft(imycorrhz,NZ,NY,NX)=1.0E+04
  RootAxialResist_pft(imycorrhz,NZ,NY,NX)=1.0E+12
!
!     RootPoreTortu4Gas=tortuosity for gas transport
!     RootRaidus_rpft=path length for radial diffusion within root (m)
!     RootVolPerMassC_pft=volume:C ratio (m3 g-1)
!     Root1stSpecLen_pft,Root2ndSpecLen_pft=specific primary,secondary root length (m g-1)
!     Root1stXSecArea_pft,Root2ndXSecArea_pft=specific primary,secondary root area (m2 g-1)
!
  D500: DO N=1,2
    RootPoreTortu4Gas(N,NZ,NY,NX)=RootPorosity_pft(N,NZ,NY,NX)**1.33_r8
    RootRaidus_rpft(N,NZ,NY,NX)=LOG(1.0_r8/SQRT(AMAX1(0.01_r8,RootPorosity_pft(N,NZ,NY,NX))))
    RootVolPerMassC_pft(N,NZ,NY,NX)=ppmc/(0.05_r8*(1.0-RootPorosity_pft(N,NZ,NY,NX)))
    Root1stSpecLen_pft(N,NZ,NY,NX)=RootVolPerMassC_pft(N,NZ,NY,NX)/(PICON*Root1stMaxRadius_pft(N,NZ,NY,NX)**2)
    Root2ndSpecLen_pft(N,NZ,NY,NX)=RootVolPerMassC_pft(N,NZ,NY,NX)/(PICON*Root2ndMaxRadius_pft(N,NZ,NY,NX)**2)
    Root1stMaxRadius1_pft(N,NZ,NY,NX)=Root1stMaxRadius_pft(N,NZ,NY,NX)
!    2*SQRT(0.25*(1.0-RootPorosity_pft(N,NZ,NY,NX)))
    Root2ndMaxRadius1_pft(N,NZ,NY,NX)=Root2ndMaxRadius_pft(N,NZ,NY,NX)
!    2*SQRT(0.25*(1.0-RootPorosity_pft(N,NZ,NY,NX)))
    Root1stXSecArea_pft(N,NZ,NY,NX)=PICON*Root1stMaxRadius1_pft(N,NZ,NY,NX)**2._r8
    Root2ndXSecArea_pft(N,NZ,NY,NX)=PICON*Root2ndMaxRadius1_pft(N,NZ,NY,NX)**2._r8
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
  PlantPopulation_pft(NZ,NY,NX)  = PPX_pft(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
  doInitPlant_pft(NZ,NY,NX)      = ifalse
  iPlantShootState_pft(NZ,NY,NX) = iLive
  iPlantRootState_pft(NZ,NY,NX)  = iLive
  BranchNumber_pft(NZ,NY,NX)     = 0
  NumOfBranches_pft(NZ,NY,NX)    = 0
  HypoctoHeight_pft(NZ,NY,NX)    = 0._r8
  CanopyHeight_pft(NZ,NY,NX)     = 0._r8
  D10: DO NB=1,MaxNumBranches
    doInitLeafOut_brch(NB,NZ,NY,NX)                = iEnable
    doPlantLeafOut_brch(NB,NZ,NY,NX)               = iEnable
    doPlantLeaveOff_brch(NB,NZ,NY,NX)              = iEnable
    Prep4Literfall_brch(NB,NZ,NY,NX)               = ifalse
    Hours4LiterfalAftMature_brch(NB,NZ,NY,NX)      = 0
    MatureGroup_brch(NB,NZ,NY,NX)                  = MatureGroup_pft(NZ,NY,NX)
    ShootNodeNum_brch(NB,NZ,NY,NX)                 = ShootNodeNumAtPlanting_pft(NZ,NY,NX)
    NodeNum2InitFloral_brch(NB,NZ,NY,NX)           = ShootNodeNum_brch(NB,NZ,NY,NX)
    NodeNumberAtAnthesis_brch(NB,NZ,NY,NX)         = 0._r8
    NumOfLeaves_brch(NB,NZ,NY,NX)                  = 0._r8
    LeafNumberAtFloralInit_brch(NB,NZ,NY,NX)       = 0._r8
    KLeafNumber_brch(NB,NZ,NY,NX)                  = 1
    KMinNumLeaf4GroAlloc_brch(NB,NZ,NY,NX)         = 1
    KHiestGroLeafNode_brch(NB,NZ,NY,NX)            = 1
    KLowestGroLeafNode_brch(NB,NZ,NY,NX)           = 0
    NodeNumNormByMatgrp_brch(NB,NZ,NY,NX)          = 0._r8
    ReprodNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX)   = 0._r8
    TotalNodeNumNormByMatgrp_brch(NB,NZ,NY,NX)     = 0._r8
    TotReproNodeNumNormByMatrgrp_brch(NB,NZ,NY,NX) = 0._r8
    Hours4LenthenPhotoPeriod_brch(NB,NZ,NY,NX)     = 0._r8
    Hours4ShortenPhotoPeriod_brch(NB,NZ,NY,NX)     = 0._r8
    Hours4Leafout_brch(NB,NZ,NY,NX)                = Hours4LenthenPhotoPeriod_brch(NB,NZ,NY,NX)
    Hours4LeafOff_brch(NB,NZ,NY,NX)                = Hours4ShortenPhotoPeriod_brch(NB,NZ,NY,NX)
    Hours2LeafOut_brch(NB,NZ,NY,NX)                = 0._r8
    RubiscoActivity_brch(NB,NZ,NY,NX)              = 1.0_r8
    C4PhotosynDowreg_brch(NB,NZ,NY,NX)             = 1.0_r8
    HourFailGrainFill_brch(NB,NZ,NY,NX)            = 0
    HoursDoingRemob_brch(NB,NZ,NY,NX)              = 0
    BranchNumber_brch(NB,NZ,NY,NX)                 = 0
    iPlantBranchState_brch(NB,NZ,NY,NX)            = iDead
    D15: DO M=1,NumGrowthStages
      iPlantCalendar_brch(M,NB,NZ,NY,NX)=0
    ENDDO D15
  ENDDO D10
!
!     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
!
  HoursTooLowPsiCan_pft(NZ,NY,NX)=0._r8
  ChillHours_pft(NZ,NY,NX)=0._r8
  CanopyNonstElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  CanopyNodulNonstElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  ShootStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  PetoleStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  StalkStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  LeafStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  StalkRsrvElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  HuskStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  GrainStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  EarStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  CanopyNodulStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  LeafElmntRemobFlx_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  PetioleChemElmRemobFlx_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  PetioleChemElmRemob_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8  
  SenecStalkStrutElms_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  LeafChemElmRemob_brch(1:NumPlantChemElms,1:MaxNumBranches,NZ,NY,NX)=0._r8
  
  D25: DO NB=1,MaxNumBranches
    StalkBiomassC_brch(NB,NZ,NY,NX)=0._r8
    LeafPetolBiomassC_brch(NB,NZ,NY,NX)=0._r8
    PotentialSeedSites_brch(NB,NZ,NY,NX)=0._r8
    SeedNumSet_brch(NB,NZ,NY,NX)=0._r8
    GrainSeedBiomCMean_brch(NB,NZ,NY,NX)=0._r8
    LeafAreaLive_brch(NB,NZ,NY,NX)=0._r8
    NH3Dep2Can_brch(NB,NZ,NY,NX)=0._r8
    LeafAreaDying_brch(NB,NZ,NY,NX)=0._r8
    CanPBranchHeight(NB,NZ,NY,NX)=0._r8
    
    D5: DO L=1,NumOfCanopyLayers
      CanopyStalkArea_lbrch(L,NB,NZ,NY,NX)=0._r8
      DO N=1,NumOfLeafZenithSectors
        StemAreaZsec_brch(N,L,NB,NZ,NY,NX)=0._r8
      enddo
    ENDDO D5

    DO K=0,MaxNodesPerBranch
      LeafAreaNode_brch(K,NB,NZ,NY,NX)                          = 0._r8
      LiveInterNodeHight_brch(K,NB,NZ,NY,NX)                    = 0._r8
      InternodeHeightDying_brch(K,NB,NZ,NY,NX)                  = 0._r8
      PetoleLensNode_brch(K,NB,NZ,NY,NX)                        = 0._r8
      LeafElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)      = 0._r8
      PetioleElmntNode_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX)   = 0._r8
      InternodeStrutElms_brch(1:NumPlantChemElms,K,NB,NZ,NY,NX) = 0._r8
      LeafProteinCNode_brch(K,NB,NZ,NY,NX)                      = 0._r8
      PetoleProteinCNode_brch(K,NB,NZ,NY,NX)                    = 0._r8

      D55: DO L=1,NumOfCanopyLayers
        CanopyLeafArea_lpft(L,K,NB,NZ,NY,NX)=0._r8
        LeafElmsByLayerNode_brch(1:NumPlantChemElms,L,K,NB,NZ,NY,NX)=0._r8
      ENDDO D55

      IF(K.NE.0)THEN
        CPOOL3_node(K,NB,NZ,NY,NX)=0._r8
        CMassCO2BundleSheath_node(K,NB,NZ,NY,NX)=0._r8
        CMassHCO3BundleSheath_node(K,NB,NZ,NY,NX)=0._r8
        CPOOL4_node(K,NB,NZ,NY,NX)=0._r8
        D45: DO L=1,NumOfCanopyLayers
          DO N=1,NumOfLeafZenithSectors
            LeafAreaZsec_brch(N,L,K,NB,NZ,NY,NX)=0._r8
          enddo
        ENDDO D45
      ENDIF
    enddo
  ENDDO D25
  D35: DO L=1,NumOfCanopyLayers
    CanopyLeafAreaZ_pft(L,NZ,NY,NX)=0._r8
    CanopyLeafCLyr_pft(L,NZ,NY,NX)=0._r8
    CanopyStemAreaZ_pft(L,NZ,NY,NX)=0._r8
  ENDDO D35
  CanopyNonstElms_pft(1:NumPlantChemElms,NZ,NY,NX)    = 0._r8
  CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ,NY,NX) = 0._r8
  NoduleNonstructCconc_pft(NZ,NY,NX)                  = 0._r8
  ShootStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)     = 0._r8
  LeafStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)      = 0._r8
  PetoleStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)    = 0._r8
  StalkStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)     = 0._r8
  CanopyStalkC_pft(NZ,NY,NX)                          = 0._r8
  StalkRsrvElms_pft(1:NumPlantChemElms,NZ,NY,NX)      = 0._r8
  HuskStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)      = 0._r8
  EarStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)       = 0._r8
  GrainStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)     = 0._r8
  RootElms_pft(1:NumPlantChemElms,NZ,NY,NX)           = 0._r8
  RootStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)      = 0._r8
  NodulStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)     = 0._r8
  CanopyLeafShethC_pft(NZ,NY,NX)                      = 0._r8

  CanopyLeafArea_pft(NZ,NY,NX)    = 0._r8
  RootBiomCPerPlant_pft(NZ,NY,NX) = 0._r8
  CanopyStemArea_pft(NZ,NY,NX)    = 0._r8

  end subroutine InitPlantPhenoMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitMassBalance(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  integer :: M,NE
  real(r8) :: WTSTDX
  associate(                 &
    icwood => pltpar%icwood  &
  )
!
!     INITIALIZE MASS BALANCE CHECKS
!
  IF(.not.is_restart().AND.is_first_year)THEN
    GrossCO2Fix_pft(NZ,NY,NX)                                   = 0._r8
    SurfLitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX) = 0._r8
    GrossResp_pft(NZ,NY,NX)                                     = 0._r8
    CanopyRespC_CumYr_pft(NZ,NY,NX)                             = 0._r8
    PlantExudElm_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)         = 0._r8
    RootUptk_N_CumYr_pft(NZ,NY,NX)                              = 0._r8
    RootUptk_P_CumYr_pft(NZ,NY,NX)                              = 0._r8
    LitrfalStrutElms_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)     = 0._r8
    PlantN2Fix_CumYr_pft(NZ,NY,NX)                              = 0._r8
    NH3Dep2Can_pft(NZ,NY,NX)                                    = 0._r8
    NH3Emis_CumYr_pft(NZ,NY,NX)                                 = 0._r8
    CO2ByFire_CumYr_pft(NZ,NY,NX)                               = 0._r8
    CH4ByFire_CumYr_pft(NZ,NY,NX)                               = 0._r8
    O2ByFire_CumYr_pft(NZ,NY,NX)                                = 0._r8
    NH3byFire_CumYr_pft(NZ,NY,NX)                               = 0._r8
    N2ObyFire_CumYr_pft(NZ,NY,NX)                               = 0._r8
    PO4byFire_CumYr_pft(NZ,NY,NX)                               = 0._r8
    EcoHavstElmntCum_pft(1:NumPlantChemElms,NZ,NY,NX)           = 0._r8
    EcoHavstElmnt_CumYr_pft(1:NumPlantChemElms,NZ,NY,NX)        = 0._r8
    NetCumElmntFlx2Plant_pft(1:NumPlantChemElms,NZ,NY,NX)       = 0._r8
    ETCanopy_CumYr_pft(NZ,NY,NX)                                = 0._r8
    StandDeadStrutElms_pft(1:NumPlantChemElms,NZ,NY,NX)         = 0._r8
    WTSTDX                                                      = StandingDeadInitC_pft(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
    D155: DO M=1,jskenc
      StandDeadKCompElms_pft(ielmc,M,NZ,NY,NX)=WTSTDX*ElmAllocMat4Litr(ielmc,icwood,M,NZ,NY,NX)
      StandDeadKCompElms_pft(ielmn,M,NZ,NY,NX)=WTSTDX*rNCStalk_pft(NZ,NY,NX)*ElmAllocMat4Litr(ielmn,icwood,M,NZ,NY,NX)
      StandDeadKCompElms_pft(ielmp,M,NZ,NY,NX)=WTSTDX*rPCStalk_pft(NZ,NY,NX)*ElmAllocMat4Litr(ielmp,icwood,M,NZ,NY,NX)
      
      DO NE=1,NumPlantChemElms
      StandDeadStrutElms_pft(NE,NZ,NY,NX)=StandDeadStrutElms_pft(NE,NZ,NY,NX) &
        +StandDeadKCompElms_pft(NE,M,NZ,NY,NX)
      ENDDO  
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
!     VHeatCapCanP_pft=canopy heat capacity (MJ m-3 K-1)
!     TCelciusCanopy_pft,TKC=canopy temperature for growth (oC,K)
!     TCGroth_pft,TKGroth_pft=canopy temperature for phenology (oC,K)
!     PSICanopy_pft,PSICanopyOsmo_pft,PSICanopyTurg_pft=canopy total,osmotic,turgor water potl(MPa)
!
  VHeatCapCanP_pft(NZ,NY,NX)       = cpw*ShootStrutElms_pft(ielmc,NZ,NY,NX)*10.0E-06
  ENGYX_pft(NZ,NY,NX)              = 0._r8
  DeltaTKC_pft(NZ,NY,NX)           = 0._r8
  TCelciusCanopy_pft(NZ,NY,NX)     = ATCA(NY,NX)
  TKC(NZ,NY,NX)                    = units%Celcius2Kelvin(TCelciusCanopy_pft(NZ,NY,NX))
  TCGroth_pft(NZ,NY,NX)            = TCelciusCanopy_pft(NZ,NY,NX)
  TKGroth_pft(NZ,NY,NX)            = units%Celcius2Kelvin(TCGroth_pft(NZ,NY,NX))
  fTCanopyGroth_pft(NZ,NY,NX)      = 1.0_r8
  PSICanopy_pft(NZ,NY,NX)          = -1.0E-03_r8
  PSICanopyOsmo_pft(NZ,NY,NX)      = CanOsmoPsi0pt_pft(NZ,NY,NX)+PSICanopy_pft(NZ,NY,NX)
  PSICanopyTurg_pft(NZ,NY,NX)      = AZMAX1(PSICanopy_pft(NZ,NY,NX)-PSICanopyOsmo_pft(NZ,NY,NX))
  Transpiration_pft(NZ,NY,NX)      = 0._r8
  FracPARRadbyCanopy_pft(NZ,NY,NX) = 0._r8
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
!     PSIRoot_pvr,PSIRootOSMO_vr,PSIRootTurg_vr=root,myco total,osmotic,turgor water potl(MPa)
!     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
!     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
!
  NumRootAxes_pft(NZ,NY,NX)=0
  RootNH4Uptake_pft(NZ,NY,NX)=0._r8
  RootNO3Uptake_pft(NZ,NY,NX)=0._r8
  RootH2PO4Uptake_pft(NZ,NY,NX)=0._r8
  RootHPO4Uptake_pft(NZ,NY,NX)=0._r8
  RootN2Fix_pft(NZ,NY,NX)=0._r8

  DO NR=1,MaxNumRootAxes
    DO N=1,pltpar%jroots
      Root1stDepz_pft(N,NR,NZ,NY,NX)=SeedDepth_pft(NZ,NY,NX)
      RootMyco1stElm_raxs(1:NumPlantChemElms,N,NR,NZ,NY,NX)=0._r8
    ENDDO
  ENDDO  
  D40: DO N=1,pltpar%jroots
    D20: DO L=1,NL(NY,NX)
      AllPlantRootH2OUptake_vr(N,L,NZ,NY,NX)=0._r8
      PSIRoot_pvr(N,L,NZ,NY,NX)=-0.01_r8
      PSIRootOSMO_vr(N,L,NZ,NY,NX)=CanOsmoPsi0pt_pft(NZ,NY,NX)+PSIRoot_pvr(N,L,NZ,NY,NX)
      PSIRootTurg_vr(N,L,NZ,NY,NX)=AZMAX1(PSIRoot_pvr(N,L,NZ,NY,NX)-PSIRootOSMO_vr(N,L,NZ,NY,NX))
      RootMycoNonstElms_rpvr(1:NumPlantChemElms,N,L,NZ,NY,NX)=0._r8
      RootNonstructElmConc_pvr(1:NumPlantChemElms,N,L,NZ,NY,NX)=0._r8
      RootProteinConc_pvr(N,L,NZ,NY,NX)=RootFracRemobilizableBiom(NZ,NY,NX)
      RootMycoActiveBiomC_pvr(N,L,NZ,NY,NX)=0._r8
      PopuRootMycoC_pvr(N,L,NZ,NY,NX)=0._r8
      RootProteinC_pvr(N,L,NZ,NY,NX)=0._r8
      Root1stXNumL_pvr(N,L,NZ,NY,NX)=0._r8
      Root2ndXNum_pvr(N,L,NZ,NY,NX)=0._r8
      RootLenPerPlant_pvr(N,L,NZ,NY,NX)=0._r8
      RootLenDensPerPlant_pvr(N,L,NZ,NY,NX)=0._r8
      RootPoreVol_pvr(N,L,NZ,NY,NX)=0._r8
      RootVH2O_pvr(N,L,NZ,NY,NX)=0._r8
      Root1stRadius_pvr(N,L,NZ,NY,NX)=Root1stMaxRadius_pft(N,NZ,NY,NX)
      Root2ndRadius_pvr(N,L,NZ,NY,NX)=Root2ndMaxRadius_pft(N,NZ,NY,NX)
      RootAreaPerPlant_pvr(N,L,NZ,NY,NX)=0._r8
      Root2ndAveLen_pvr(N,L,NZ,NY,NX)=1.0E-03
      RootNutUptake_pvr(ids_NH4,N,L,NZ,NY,NX)=0._r8
      RootNutUptake_pvr(ids_NO3,N,L,NZ,NY,NX)=0._r8
      RootNutUptake_pvr(ids_H2PO4,N,L,NZ,NY,NX)=0._r8
      RootNutUptake_pvr(ids_H1PO4,N,L,NZ,NY,NX)=0._r8
      RootNutUptake_pvr(ids_NH4B,N,L,NZ,NY,NX)=0._r8
      RootNutUptake_pvr(ids_NO3B,N,L,NZ,NY,NX)=0._r8
      RootNutUptake_pvr(ids_H2PO4B,N,L,NZ,NY,NX)=0._r8
      RootNutUptake_pvr(ids_H1PO4B,N,L,NZ,NY,NX)=0._r8
      RootO2Dmnd4Resp_pvr(N,L,NZ,NY,NX)=0._r8
      RootNH4DmndSoil_pvr(N,L,NZ,NY,NX)=0._r8
      RootNH4DmndBand_pvr(N,L,NZ,NY,NX)=0._r8
      RootNO3DmndSoil_pvr(N,L,NZ,NY,NX)=0._r8
      RootNO3DmndBand_pvr(N,L,NZ,NY,NX)=0._r8
      RootH2PO4DmndSoil_pvr(N,L,NZ,NY,NX)=0._r8
      RootH1PO4DmndSoil_pvr(N,L,NZ,NY,NX)=0._r8
      RootH2PO4DmndBand_pvr(N,L,NZ,NY,NX)=0._r8
      RootH1PO4DmndBand_pvr(N,L,NZ,NY,NX)=0._r8
      CCO2A=CCO2EI(NY,NX)
      CCO2P=0.030_r8*EXP(-2.621_r8-0.0317_r8*ATCA(NY,NX))*CO2EI(NY,NX)
      trcg_rootml_pvr(idg_CO2,N,L,NZ,NY,NX)=CCO2A*RootPoreVol_pvr(N,L,NZ,NY,NX)
      trcs_rootml_pvr(idg_CO2,N,L,NZ,NY,NX)=CCO2P*RootVH2O_pvr(N,L,NZ,NY,NX)
      trcg_air2root_flx__pvr(idg_CO2,N,L,NZ,NY,NX)=0._r8
      trcg_Root_DisEvap_flx_vr(idg_CO2,N,L,NZ,NY,NX)=0._r8
      RUPGasSol_vr(idg_CO2,N,L,NZ,NY,NX)=0._r8
      RootCO2Emis_pvr(N,L,NZ,NY,NX)=0._r8
      COXYA=AtmGasCgperm3(idg_O2,NY,NX)
      COXYP=0.032_r8*EXP(-6.175_r8-0.0211_r8*ATCA(NY,NX))*OXYE(NY,NX)
      trcg_rootml_pvr(idg_beg:idg_end-1,N,L,NZ,NY,NX)=0._r8
      trcs_rootml_pvr(idg_beg:idg_end-1,N,L,NZ,NY,NX)=0._r8
      trcg_rootml_pvr(idg_O2,N,L,NZ,NY,NX)=COXYA*RootPoreVol_pvr(N,L,NZ,NY,NX)
      trcs_rootml_pvr(idg_O2,N,L,NZ,NY,NX)=COXYP*RootVH2O_pvr(N,L,NZ,NY,NX)
      RAutoRootO2Limter_pvr(N,L,NZ,NY,NX)=1.0_r8

      D30: DO NR=1,MaxNumRootAxes
        Root2ndXNum_rpvr(N,L,NR,NZ,NY,NX)=0._r8
        Root1stLen_rpvr(N,L,NR,NZ,NY,NX)=0._r8
        RootMyco1stStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ,NY,NX)=0._r8
        Root2ndLen_pvr(N,L,NR,NZ,NY,NX)=0._r8
        RootMyco2ndStrutElms_rpvr(1:NumPlantChemElms,N,L,NR,NZ,NY,NX)=0._r8
      ENDDO D30

      IF(N.EQ.1)THEN
        D6400: DO K=1,pltpar%NumOfPlantLitrCmplxs
          DO  M=1,jskenc
            LitrfalStrutElms_pvr(1:NumPlantChemElms,M,K,L,NZ,NY,NX)=0._r8
          enddo
        ENDDO D6400
        RootNodulNonstElms_pvr(1:NumPlantChemElms,L,NZ,NY,NX)=0._r8
        RootNodulStrutElms_pvr(1:NumPlantChemElms,L,NZ,NY,NX)=0._r8
        RootN2Fix_pvr(L,NZ,NY,NX)=0._r8
      ENDIF
    ENDDO D20
  ENDDO D40

  RootNutUptake_pvr(ids_NH4,1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)    = 0._r8
  RootNutUptake_pvr(ids_NH4B,1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)   = 0._r8
  RootNutUptake_pvr(ids_H2PO4,1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)  = 0._r8
  RootNutUptake_pvr(ids_H2PO4B,1:2,NL(NY,NX)+1:JZ,NZ,NY,NX) = 0._r8
  RootLenDensPerPlant_pvr(1:2,NL(NY,NX)+1:JZ,NZ,NY,NX)      = 0._r8
  end subroutine InitRootMychorMorphoBio
!------------------------------------------------------------------------------------------

  subroutine InitSeedMorphoBio(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ, NY, NX
  REAL(R8) :: FDM
  INTEGER :: NR,L,NB,NE
!
!     INITIALIZE SEED MORPHOLOGY AND BIOMASS
!
!     WTRVC,WTRVN,WTRVP=C,N,P in storage reserves (g)
!     WTLFB,WTLFBN,WTLFBP=C,N,P in leaves (g)
!     LeafPetolBiomassC_brch=C in leaves+petioles (g)
!     FDM-dry matter fraction (g DM C g FM C-1)
!     CanopyWater_pft,WatByPCanopy_pft=water volume in,on canopy (m3)
!     CPOOL,ZPOOL,PPOOL=C,N,P in canopy nonstructural pools (g)
!     WTRT1,WTRT1N,WTRT1P=C,N,P in primary root layer (g)
!     RTWT1,RTWT1N,RTWT1P=total C,N,P in primary root (g)
!     RootMycoActiveBiomC_pvr,WTRTD=total root C mass (g)
!     RootProteinC_pvr=total root protein C mass (g)
!     CPOOLR,ZPOOLR,PPOOLR=C,N,P in root,myco nonstructural pools (g)
!
  SeedCPlanted_pft(NZ,NY,NX)            = SeedCMass_pft(NZ,NY,NX)*PlantPopulation_pft(NZ,NY,NX)
  SeasonalNonstElms_pft(ielmc,NZ,NY,NX) = SeedCPlanted_pft(NZ,NY,NX)
  SeasonalNonstElms_pft(ielmn,NZ,NY,NX) = CNGR(NZ,NY,NX)*SeasonalNonstElms_pft(ielmc,NZ,NY,NX)
  SeasonalNonstElms_pft(ielmp,NZ,NY,NX) = CPGR(NZ,NY,NX)*SeasonalNonstElms_pft(ielmc,NZ,NY,NX)
  LeafStrutElms_brch(ielmn,1,NZ,NY,NX)  = CNGR(NZ,NY,NX)*LeafStrutElms_brch(ielmc,1,NZ,NY,NX)
  LeafStrutElms_brch(ielmp,1,NZ,NY,NX)  = CPGR(NZ,NY,NX)*LeafStrutElms_brch(ielmc,1,NZ,NY,NX)
  LeafPetolBiomassC_brch(1,NZ,NY,NX)    = LeafStrutElms_brch(ielmc,1,NZ,NY,NX)+PetoleStrutElms_brch(ielmc,1,NZ,NY,NX)
  CanopyLeafShethC_pft(NZ,NY,NX)        = CanopyLeafShethC_pft(NZ,NY,NX)+LeafPetolBiomassC_brch(1,NZ,NY,NX)
    
  FDM=AMIN1(1.0_r8,0.16_r8-0.045_r8*PSICanopy_pft(NZ,NY,NX))
  CanopyWater_pft(NZ,NY,NX)=ppmc*CanopyLeafShethC_pft(NZ,NY,NX)/FDM
  WatByPCanopy_pft(NZ,NY,NX)=0._r8
  CanopyNonstElms_brch(ielmn,1,NZ,NY,NX)=CNGR(NZ,NY,NX)*CanopyNonstElms_brch(ielmc,1,NZ,NY,NX)
  CanopyNonstElms_brch(ielmp,1,NZ,NY,NX)=CPGR(NZ,NY,NX)*CanopyNonstElms_brch(ielmc,1,NZ,NY,NX)
  
  RootMyco1stStrutElms_rpvr(ielmn,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)= &
    CNGR(NZ,NY,NX)*RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)
  RootMyco1stStrutElms_rpvr(ielmp,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)= &
    CPGR(NZ,NY,NX)*RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)
  RootMyco1stElm_raxs(ielmn,ipltroot,1,NZ,NY,NX)=CNGR(NZ,NY,NX)*RootMyco1stElm_raxs(ielmc,ipltroot,1,NZ,NY,NX)
  RootMyco1stElm_raxs(ielmp,ipltroot,1,NZ,NY,NX)=CPGR(NZ,NY,NX)*RootMyco1stElm_raxs(ielmc,ipltroot,1,NZ,NY,NX)
  RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)= &
    RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)
  PopuRootMycoC_pvr(ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)= &
    RootMyco1stStrutElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),1,NZ,NY,NX)
  RootProteinC_pvr(ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)= &
    RootMycoActiveBiomC_pvr(ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)*RootFracRemobilizableBiom(NZ,NY,NX)
  RootMycoNonstElms_rpvr(ielmn,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)= CNGR(NZ,NY,NX)&
    *RootMycoNonstElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)
  RootMycoNonstElms_rpvr(ielmp,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)=CPGR(NZ,NY,NX) &
    *RootMycoNonstElms_rpvr(ielmc,ipltroot,NGTopRootLayer_pft(NZ,NY,NX),NZ,NY,NX)


  !shoots
  NB=1
  DO NE=1,NumPlantChemElms
    ShootStrutElms_brch(NE,NB,NZ,NY,NX)=LeafStrutElms_brch(NE,NB,NZ,NY,NX) &
      +PetoleStrutElms_brch(NE,NB,NZ,NY,NX)+StalkStrutElms_brch(NE,NB,NZ,NY,NX) &
      +StalkRsrvElms_brch(NE,NB,NZ,NY,NX)+HuskStrutElms_brch(NE,NB,NZ,NY,NX)&
      +EarStrutElms_brch(NE,NB,NZ,NY,NX)+GrainStrutElms_brch(NE,NB,NZ,NY,NX) 
        
    ShootElms_brch(NE,NB,NZ,NY,NX)=ShootStrutElms_brch(NE,NB,NZ,NY,NX)+CanopyNonstElms_brch(NE,NB,NZ,NY,NX)
  ENDDO
  
  
  DO NE=1,NumPlantChemElms
    ShootElms_pft(NE,NZ,NY,NX)=ShootElms_brch(NE,1,NZ,NY,NX)+SeasonalNonstElms_pft(NE,NZ,NY,NX) 
  ENDDO

  RootElms_pft(1:NumPlantChemElms,NZ,NY,NX)=0._r8
  !roots
  DO NR=1,MaxNumRootAxes
    DO L=NU(NY,NX),NGTopRootLayer_pft(NZ,NY,NX)
      DO NE=1,NumPlantChemElms
        !add reserve to struct
        RootElms_pft(NE,NZ,NY,NX)=RootElms_pft(NE,NZ,NY,NX)+RootMyco1stStrutElms_rpvr(NE,ipltroot,L,NR,NZ,NY,NX)&
          +RootMyco2ndStrutElms_rpvr(NE,ipltroot,L,NR,NZ,NY,NX)
      ENDDO
    ENDDO    
  ENDDO
  
  DO L=NU(NY,NX),NGTopRootLayer_pft(NZ,NY,NX)
    DO NE=1,NumPlantChemElms
        !add reserve to struct
      RootElms_pft(NE,NZ,NY,NX)=RootElms_pft(NE,NZ,NY,NX)+RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ,NY,NX)
    ENDDO
  ENDDO    

  end subroutine InitSeedMorphoBio

  end module StartqMod
