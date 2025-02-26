module UptakesMod
  use data_kind_mod , only : r8 => DAT_KIND_R8
  use data_const_mod, only : GravAcceleration=>DAT_CONST_G
  use StomatesMod   , only : StomatalDynamics
  use EcoSIMCtrlMod , only : etimer, lverb  
  use UnitMod       , only : units  
  use GrosubPars
  use minimathmod
  use DebugToolMod
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use NutUptakeMod
  use PlantAPIData
  use PlantMathFuncMod
  use ElmIDMod, only : itrue
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public :: RootUptakes
  public :: InitUptake

  real(r8), parameter :: mGravAccelerat=1.e-3_r8*GravAcceleration  !gravitational constant devided by 1000, MPa/m.
  contains

  subroutine InitUptake

  implicit none

  call InitUptakePars

  end subroutine InitUptake
!------------------------------------------------------------------------------------------

  subroutine RootUptakes(I,J)
!
!     THIS subroutine CALCULATES EXCHANGES OF ENERGY, C, N AND P
!     BETWEEN THE CANOPY AND THE ATMOSPHERE AND BETWEEN ROOTS AND THE SOIL
!
  implicit none
  integer, intent(in) :: I, J
  
  character(len=*), parameter :: subname='RootUptakes'
  integer :: NN,N,NZ,K,L
  real(r8) :: FracGrndByPFT
  real(r8) :: VHeatCapCanopyAir
  real(r8) :: CanopyMassC                !Canopy C mass, g/m2
  real(r8) :: TotalSoilPSIMPa_vr(JZ1)    !total soil matric pressure, removing elevation adjustment [MPa]
  real(r8) :: PathLen_pvr(jroots,JZ1)
  real(r8) :: FineRootRadius_rvr(jroots,JZ1)
  real(r8) :: RootResist(jroots,JZ1)
  real(r8) :: RootResistSoi(jroots,JZ1)
  real(r8) :: RootResistPrimary(jroots,JZ1)
  real(r8) :: RootResist2ndary(jroots,JZ1)
  real(r8) :: SoiH2OResist(jroots,JZ1)
  real(r8) :: SoilRootResistance_rvr(jroots,JZ1),AllRootC_vr(JZ1)
  real(r8) :: FracPRoot4Uptake_pvr(jroots,JZ1,JP1)
  real(r8) :: MinFracPRoot4Uptake_pvr(jroots,JZ1,JP1)
  real(r8) :: FracPrimRootOccupiedLay_pvr(JZ1,JP1)
  real(r8) :: RootAreaDivRadius_vr(jroots,JZ1)
  real(r8) :: AirMicPore4Fill_vr(JZ1)  !
  real(r8) :: WatAvail4Uptake_vr(JZ1)
  real(r8) :: TKCX,CNDT,PrecpHeatbyCanopy
  real(r8) :: VHeatCapCanopyPrev_pft  !canopy heat capacity, J/K
  real(r8) :: DIFF
  real(r8) :: PSIGravCanopyHeight               !gravitation potential at effective canopy height [MPa]
  real(r8) :: cumPRootH2OUptake,HeatEvapSens
  real(r8) :: CumPlantHeatLoss2Soil
  real(r8) :: FDMP
  logical  :: HydroActivePlant
  logical  :: SoiLayerHasRoot_rvr(jroots,JZ1)
!     begin_execution
  associate(                                                       &
    CanopyBndlResist_pft      => plt_photo%CanopyBndlResist_pft,   &
    CanopyBiomWater_pft       => plt_ew%CanopyBiomWater_pft,       &
    TKCanopy_pft              => plt_ew%TKCanopy_pft,              &
    PrecIntcptByCanopy_pft    => plt_ew%PrecIntcptByCanopy_pft,    &
    EvapTransLHeat_pft        => plt_ew%EvapTransLHeat_pft,        &
    TKC_pft                   => plt_ew%TKC_pft,                   &
    PSICanopy_pft             => plt_ew%PSICanopy_pft,             &
    HeatXAir2PCan_pft         => plt_ew%HeatXAir2PCan_pft,         &
    Air_Heat_Sens_store_col   => plt_ew%Air_Heat_Sens_store_col,   &
    VapXAir2Canopy_pft        => plt_ew%VapXAir2Canopy_pft,        &
    TairK                     => plt_ew%TairK,                     &
    Air_Heat_Latent_store_col => plt_ew%Air_Heat_Latent_store_col, &
    Transpiration_pft         => plt_ew%Transpiration_pft,         &
    WatHeldOnCanopy_pft       => plt_ew%WatHeldOnCanopy_pft,       &
    CumSoilThickness_vr       => plt_site%CumSoilThickness_vr,     &
    PlantPopu_col             => plt_site%PlantPopu_col,           &
    NP                        => plt_site%NP,                      &
    PlantPopulation_pft       => plt_site%PlantPopulation_pft,     &
    NU                        => plt_site%NU,                      &
    AREA3                     => plt_site%AREA3,                   &
    ZEROS                     => plt_site%ZEROS,                   &
    CanopyMassC_pft           => plt_biom%CanopyMassC_pft,         &
    PopuRootMycoC_pvr         => plt_biom% PopuRootMycoC_pvr,      &
    ZERO4LeafVar_pft          => plt_biom%ZERO4LeafVar_pft,        &
    CanopyLeafShethC_pft      => plt_biom%CanopyLeafShethC_pft,    &
    CanopyStalkC_pft          => plt_biom%CanopyStalkC_pft,        &
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch,    &
    IsPlantActive_pft         => plt_pheno%IsPlantActive_pft,      &
    CanopyLeafArea_col        => plt_morph%CanopyLeafArea_col,     &
    Root1stDepz_pft           => plt_morph%Root1stDepz_pft,        &
    LeafStalkArea_col         => plt_morph%LeafStalkArea_col,      &
    LeafStalkArea_pft         => plt_morph%LeafStalkArea_pft,      &
    SeedDepth_pft             => plt_morph%SeedDepth_pft,          &
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft,      &
    FracPARads2Canopy_pft     => plt_rad%FracPARads2Canopy_pft     &
  )

  call PrintInfo('beg '//subname)

  call PrepH2ONutrientUptake(TotalSoilPSIMPa_vr,AllRootC_vr,AirMicPore4Fill_vr,WatAvail4Uptake_vr)
!
!     IF PLANT SPECIES EXISTS

  DO NZ=1,NP

    IF(IsPlantActive_pft(NZ).EQ.iActive .AND. PlantPopulation_pft(NZ).GT.0.0_r8)THEN
      
      call UpdateCanopyProperty(NZ)

!     STOMATE=solve for minimum canopy stomatal resistance
      
      CALL StomatalDynamics(I,J,NZ)
!
      if(lverb)write(*,*)'CALCULATE VARIABLES USED IN ROOT UPTAKE OF WATER AND NUTRIENTS'
      call UpdateRootProperty(NZ,PathLen_pvr,FineRootRadius_rvr,AllRootC_vr,FracPRoot4Uptake_pvr,&
        MinFracPRoot4Uptake_pvr,FracPrimRootOccupiedLay_pvr,RootAreaDivRadius_vr)
!
!     CALCULATE CANOPY WATER STATUS FROM CONVERGENCE SOLUTION FOR
!     TRANSPIRATION - ROOT WATER UPTAKE = CHANGE IN CANOPY WATER CONTENT
!
      CanopyMassC            = AZMAX1(CanopyLeafShethC_pft(NZ)+CanopyStalkC_pft(NZ))
      VHeatCapCanopyPrev_pft = cpw*(CanopyMassC*SpecStalkVolume+WatHeldOnCanopy_pft(NZ)+CanopyBiomWater_pft(NZ))

      HydroActivePlant=(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).NE.0)             &  !plant emerged
        .AND.(LeafStalkArea_pft(NZ).GT.ZERO4LeafVar_pft(NZ).AND.FracPARads2Canopy_pft(NZ).GT.0.0_r8)   &  !active canopy
        .AND.(Root1stDepz_pft(ipltroot,1,NZ).GT.SeedDepth_pft(NZ)+CumSoilThickness_vr(0))              &  !active root
        .and. CanopyMassC>0._r8

      IF(HydroActivePlant)THEN  
!
        call CalcResistance(NZ,PathLen_pvr,FineRootRadius_rvr,RootAreaDivRadius_vr,RootResist,RootResistSoi,&
          RootResistPrimary,RootResist2ndary,SoiH2OResist,SoilRootResistance_rvr,CNDT,PSIGravCanopyHeight,SoiLayerHasRoot_rvr)
!
!       INITIALIZE CANOPY WATER POTENTIAL, OTHER VARIABLES USED IN ENERGY
!       BALANCE THAT DON'T NEED TO BE RECALCULATED DURING CONVERGENCE
!
        PSICanopy_pft(NZ)      = AMIN1(-ppmc,0.667_r8*PSICanopy_pft(NZ))
        Transpiration_pft(NZ)  = 0.0_r8
        VapXAir2Canopy_pft(NZ) = 0.0_r8
        PrecpHeatbyCanopy      = PrecIntcptByCanopy_pft(NZ)*cpw*TairK

        IF(LeafStalkArea_col.GT.ZEROS)THEN
          !the grid has significant canopy (leaf+steam) area
          FracGrndByPFT=LeafStalkArea_pft(NZ)/LeafStalkArea_col*AMIN1(1.0_r8,0.5_r8*CanopyLeafArea_col/AREA3(NU))
        ELSEIF(PlantPopu_col.GT.ZEROS)THEN
          !total population is > 0
          FracGrndByPFT=PlantPopulation_pft(NZ)/PlantPopu_col
        ELSE
          FracGrndByPFT=1.0_r8/NP
        ENDIF

        TKCX   = TKC_pft(NZ)
!
!     CONVERGENCE SOLUTION
!
        NN=CanopyEnergyH2OIteration(I,J,NZ,FracGrndByPFT,CanopyMassC,&
          TotalSoilPSIMPa_vr,VHeatCapCanopyAir,DIFF,cumPRootH2OUptake,CumPlantHeatLoss2Soil,&
          HeatEvapSens,FDMP,SoilRootResistance_rvr,FracPRoot4Uptake_pvr,AirMicPore4Fill_vr,&
          WatAvail4Uptake_vr,TKCX,CNDT,VHeatCapCanopyPrev_pft,PrecpHeatbyCanopy,PSIGravCanopyHeight,SoiLayerHasRoot_rvr)
!
!     IF CONVERGENCE NOT ACHIEVED (RARE), SET DEFAULT
!     TEMPERATURES, ENERGY FLUXES, WATER POTENTIALS, RESISTANCES
!
        call HandlingDivergence(I,J,NN,NZ,TotalSoilPSIMPa_vr,DIFF,FDMP)

        call UpdatePlantWaterVars(NZ,VHeatCapCanopyAir,TotalSoilPSIMPa_vr,RootResist,SoiH2OResist,SoilRootResistance_rvr,&
          TKCX,VHeatCapCanopyPrev_pft,PrecpHeatbyCanopy,cumPRootH2OUptake,CumPlantHeatLoss2Soil,HeatEvapSens,SoiLayerHasRoot_rvr)
!
!     DEFAULT VALUES IF PLANT SPECIES DOES NOT EXIST
!
      ELSE
        call HandleBareSoil(NZ,TotalSoilPSIMPa_vr,FDMP)
      ENDIF

      Air_Heat_Latent_store_col = Air_Heat_Latent_store_col+EvapTransLHeat_pft(NZ)*CanopyBndlResist_pft(NZ)
      Air_Heat_Sens_store_col   = Air_Heat_Sens_store_col+HeatXAir2PCan_pft(NZ)*CanopyBndlResist_pft(NZ)
      
      call SetCanopyGrowthFuncs(I,J,NZ)
    
      call PlantNutientO2Uptake(I,J,NZ,FDMP,PathLen_pvr,FineRootRadius_rvr,FracPRoot4Uptake_pvr,&
        MinFracPRoot4Uptake_pvr,FracPrimRootOccupiedLay_pvr,RootAreaDivRadius_vr)

    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
  RETURN
  end associate
  END subroutine RootUptakes
!------------------------------------------------------------------------

  subroutine PrepH2ONutrientUptake(TotalSoilPSIMPa_vr,AllRootC_vr,AirMicPore4Fill_vr,WatAvail4Uptake_vr)
!
!     prepare for uptake calculation
  implicit none
  real(r8), intent(out) :: TotalSoilPSIMPa_vr(JZ1)   !total soil matric pressure [MPa]
  real(r8), intent(out) :: AllRootC_vr(JZ1)
  real(r8), intent(out) :: AirMicPore4Fill_vr(JZ1)
  real(r8), intent(out) :: WatAvail4Uptake_vr(JZ1)

  character(len=*), parameter :: subname='PrepH2ONutrientUptake'  
  integer :: NZ, L, N
  real(r8) :: ARLSC

  associate(                                                         &
    ZERO                       => plt_site%ZERO,                     &
    NP                         => plt_site%NP,                       &
    MaxNumRootLays             => plt_site%MaxNumRootLays,           &
    VLWatMicPM_vr              => plt_site%VLWatMicPM_vr,            &
    ALT                        => plt_site%ALT,                      &
    NU                         => plt_site%NU,                       &
    NP0                        => plt_site%NP0,                      &
    ElvAdjstedSoilH2OPSIMPa_vr => plt_ew%ElvAdjstedSoilH2OPSIMPa_vr, &
    PopuRootMycoC_pvr          => plt_biom% PopuRootMycoC_pvr,       &
    VLMicP_vr                  => plt_soilchem%VLMicP_vr,            &
    VLiceMicP_vr               => plt_soilchem%VLiceMicP_vr,         &
    SoilWatAirDry_vr           => plt_soilchem%SoilWatAirDry_vr,     &
    SoilBulkDensity_vr         => plt_soilchem%SoilBulkDensity_vr,   &
    VLWatMicP_vr               => plt_soilchem%VLWatMicP_vr,         &
    VLSoilMicP_vr              => plt_soilchem%VLSoilMicP_vr,        &
    CanopyStemArea_pft         => plt_morph%CanopyStemArea_pft,      &
    MY                         => plt_morph%MY,                      &
    CanopyLeafArea_pft         => plt_morph%CanopyLeafArea_pft,      &
    RadNet2Canopy_pft          => plt_rad%RadNet2Canopy_pft,         &
    LWRadCanopy_pft            => plt_rad%LWRadCanopy_pft            &
  )
!
!     RESET TOTAL UPTAKE ARRAYS
!
!     CanopyLeafArea_pft,CanopyStemArea_pft=leaf,stalk areas
!
  call PrintInfo('beg '//subname)

  ARLSC=0.0_r8
  D9984: DO NZ=1,NP0
!     TKC_pft(NZ)=TairK+DeltaTKC_pft(NZ)
!     TdegCCanopy_pft(NZ)=TKC_pft(NZ)-TC2K
    ARLSC   = ARLSC+CanopyLeafArea_pft(NZ)+CanopyStemArea_pft(NZ)
    RadNet2Canopy_pft(NZ)                                = 0.0_r8
    plt_ew%EvapTransLHeat_pft(NZ)                         = 0.0_r8
    plt_ew%HeatXAir2PCan_pft(NZ)                         = 0.0_r8
    plt_ew%HeatStorCanopy_pft(NZ)                        = 0.0_r8
    LWRadCanopy_pft(NZ)                                  = 0.0_r8
    plt_ew%Transpiration_pft(NZ)                         = 0.0_r8
    plt_ew%VapXAir2Canopy_pft(NZ)                        = 0.0_r8
    plt_rbgc%RootMycoExudElms_pft(1:NumPlantChemElms,NZ) = 0.0_r8
    plt_rbgc%RootNH4Uptake_pft(NZ)                       = 0.0_r8
    plt_rbgc%RootNO3Uptake_pft(NZ)                       = 0.0_r8
    plt_rbgc%RootH2PO4Uptake_pft(NZ)                     = 0.0_r8
    plt_rbgc%RootHPO4Uptake_pft(NZ)                      = 0.0_r8
    plt_rbgc%RootN2Fix_pft(NZ)                           = 0.0_r8
!
!     RESET UPTAKE ARRAYS
!
    DO  L=NU,MaxNumRootLays
      DO  N=1,MY(NZ)
        plt_ew%AllPlantRootH2OLoss_vr(N,L,NZ)                   = 0.0_r8
        plt_rbgc%RootCO2Emis_pvr(N,L,NZ)                          = 0.0_r8
        plt_rbgc%RootO2Uptk_pvr(N,L,NZ)                           = 0.0_r8
        plt_rbgc%RootUptkSoiSol_vr(idg_beg:idg_end,N,L,NZ)        = 0.0_r8
        plt_rbgc%trcg_air2root_flx_pvr(idg_beg:idg_NH3,N,L,NZ)    = 0.0_r8
        plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_beg:idg_NH3,N,L,NZ) = 0.0_r8
      enddo
    enddo
  ENDDO D9984
!
! NPH is the last iteration from solving for soil heat-moisture hydrothermal dynamics 
  D9000: DO L=NU,MaxNumRootLays
    !remove elevation dependence
    TotalSoilPSIMPa_vr(L)=ElvAdjstedSoilH2OPSIMPa_vr(L)-mGravAccelerat*ALT
    IF(SoilBulkDensity_vr(L).GT.ZERO)THEN
      WatAvail4Uptake_vr(L) = VLWatMicPM_vr(NPH,L)-SoilWatAirDry_vr(L)*VLSoilMicP_vr(L)     !maximum amount of water for uptake
      AirMicPore4Fill_vr(L) = AZMAX1(VLMicP_vr(L)-VLWatMicP_vr(L)-VLiceMicP_vr(L))   !air volume in soil
    ELSE
      WatAvail4Uptake_vr(L) = VLWatMicPM_vr(NPH,L)
      AirMicPore4Fill_vr(L) = 0.0_r8
    ENDIF
    AllRootC_vr(L)=0.0_r8
    D9005: DO NZ=1,NP
      DO  N=1,MY(NZ)
      !turn off to include dead plants as well
!     IF(IsPlantActive_pft(NZ).EQ.iActive .AND. PlantPopulation_pft(NZ).GT.0.0)THEN
        AllRootC_vr(L)=AllRootC_vr(L)+AZMAX1(PopuRootMycoC_pvr(N,L,NZ))
!     ENDIF
      enddo
    ENDDO D9005
  ENDDO D9000

  call PrintInfo('end '//subname)
  end associate
  end subroutine PrepH2ONutrientUptake
!------------------------------------------------------------------------
  subroutine UpdateCanopyProperty(NZ)
!
!     update canopy characterization
  implicit none
  integer, intent(in) :: NZ

  character(len=*), parameter :: subname='UpdateCanopyProperty'
  real(r8) :: ALFZ       !shape parameter for windspeed attenuation in canopy
  real(r8) :: TFRADP     !total fraction of PAR on canopy
  real(r8) :: BndlCanopyReistance_pft(JP1)
  integer :: NB,K,L,N,NZZ

  associate(                                                          &
    ZERO                      => plt_site%ZERO,                       &
    NP                        => plt_site%NP,                         &
    KoppenClimZone            => plt_site%KoppenClimZone,             &
    AbvCanopyBndlResist_col   => plt_ew%AbvCanopyBndlResist_col,      & !input:
    TairK                     => plt_ew%TairK,                        & !input:
    DeltaTKC_pft              => plt_ew%DeltaTKC_pft,                 & !input:
    ZERO4PlantDisplace_col    => plt_ew%ZERO4PlantDisplace_col,       & !input:
    RoughHeight               => plt_ew%RoughHeight,                  &
    ReistanceCanopy_pft       => plt_ew%ReistanceCanopy_pft,          &
    TKCanopy_pft              => plt_ew%TKCanopy_pft,                 &
    LeafProteinCNode_brch     => plt_biom%LeafProteinCNode_brch,      &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,             &
    FracPARads2Canopy_pft     => plt_rad%FracPARads2Canopy_pft,       &
    LeafAUnshaded_zsec        => plt_photo%LeafAUnshaded_zsec,        &
    LeafNodeArea_brch         => plt_morph%LeafNodeArea_brch,         &
    KMinNumLeaf4GroAlloc_brch => plt_morph%KMinNumLeaf4GroAlloc_brch, &
    CanopyHeight_pft          => plt_morph%CanopyHeight_pft,          & !input:
    ClumpFactorNow_pft        => plt_morph%ClumpFactorNow_pft,        &
    LeafStalkArea_pft         => plt_morph%LeafStalkArea_pft,         &
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft,         &
    LeafAreaZsec_brch         => plt_morph%LeafAreaZsec_brch,         &
    CanopyHeight_col          => plt_morph%CanopyHeight_col           & !input:
  )
!
!     APPLY CLUMPING FACTOR TO LEAF SURFACE AREA DEFINED BY
!     INCLINATION N, LAYER L, NODE K, BRANCH NB, SPECIES NZ,
!     N-S POSITION NY, E-W POSITION NX(AZIMUTH M ASSUMED UNIFORM)
!
  call PrintInfo('beg '//subname)
  D500: DO NB=1,NumOfBranches_pft(NZ)
    D550: DO K=1,MaxNodesPerBranch1
!
!     NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
!
!     LeafNodeArea_brch=leaf area
!     LeafProteinCNode_brch=leaf protein content
!     LeafAUnshaded_zsec,LeafAreaZsec_brch=unself-shaded,total leaf surface area
!     ClumpFactorNow_pft=clumping factor from PFT file
!
      IF(LeafNodeArea_brch(K,NB,NZ).GT.ZERO4Groth_pft(NZ) .AND. LeafProteinCNode_brch(K,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
        KMinNumLeaf4GroAlloc_brch(NB,NZ)=K
      ENDIF

      D600: DO L=NumOfCanopyLayers1,1,-1
        D650: DO N=1,NumOfLeafZenithSectors1
          LeafAUnshaded_zsec(N,L,K,NB,NZ)=LeafAreaZsec_brch(N,L,K,NB,NZ)*ClumpFactorNow_pft(NZ)
        ENDDO D650
      ENDDO D600
    ENDDO D550
  ENDDO D500
!
  if(lverb)write(*,*)'   CANOPY HEIGHT FROM HEIGHT OF MAXIMUM LEAF LAYER',LeafStalkArea_pft(NZ)
!
!     DIFFUSIVE RESISTANCE OF OTHER TALLER CANOPIES TO HEAT AND VAPOR
!     TRANSFER OF CURRENT CANOPY ADDED TO BOUNDARY LAYER RESISTANCE
!     OF TALLEST CANOPY CALCULATED IN 'HOUR1'
!
!     KoppenClimZone=Koppen climate zone
!     ZC,ZT,RoughHeight=PFT canopy,grid biome,surface roughness height
!     FracPARads2Canopy_pft=fraction of radiation received by each PFT canopy
!     AbvCanopyBndlResist_col=biome canopy isothermal boundary layer resistance
!     BndlCanopyReistance_pft,ReistanceCanopy_pft=additional,total PFT canopy isothermal blr
!
  IF(LeafStalkArea_pft(NZ).GT.0.0_r8)THEN
    IF(KoppenClimZone.GE.0)THEN
      TFRADP=0.0_r8
      D700: DO NZZ=1,NP
        IF(CanopyHeight_pft(NZZ).GT.CanopyHeight_pft(NZ)+RoughHeight)THEN
          TFRADP=TFRADP+FracPARads2Canopy_pft(NZZ)
        ENDIF
      ENDDO D700
      ALFZ=2.0_r8*TFRADP
      IF(AbvCanopyBndlResist_col.GT.ZERO .AND. CanopyHeight_col.GT.ZERO .AND. ALFZ.GT.ZERO)THEN
        BndlCanopyReistance_pft(NZ)=AMIN1(RACX,AZMAX1(CanopyHeight_col*EXP(ALFZ) &
          /(ALFZ/AbvCanopyBndlResist_col)*(EXP(-ALFZ*CanopyHeight_pft(NZ)/CanopyHeight_col) &
          -EXP(-ALFZ*(ZERO4PlantDisplace_col+RoughHeight)/CanopyHeight_col))))
      ELSE
        BndlCanopyReistance_pft(NZ)=0.0_r8
      ENDIF
    ELSE
      BndlCanopyReistance_pft(NZ)=0.0_r8
    ENDIF
  ELSE
    BndlCanopyReistance_pft(NZ)=RACX
  ENDIF
  ReistanceCanopy_pft(NZ)=AbvCanopyBndlResist_col+BndlCanopyReistance_pft(NZ)
!
!     INITIALIZE CANOPY TEMPERATURE WITH CURRENT AIR TEMPERATURE AND
!     LAST HOUR'S CANOPY-AIR TEMPERATURE DIFFERENCE, AND CALL A
!     subroutine TO CALCULATE MINIMUM CANOPY STOMATAL RESISTANCE
!     FOR SUBSEQUENT USE IN ENERGY EXCHANGE CALCULATIONS
!
!     TKCanopy_pft=initial estimate of canopy temperature TKC
!     TairK=current air temperature
!     DeltaTKC_pft=TKC-TairK from previous hour
!
  TKCanopy_pft(NZ)=TairK+DeltaTKC_pft(NZ)

  call PrintInfo('end '//subname)
  end associate
  end subroutine UpdateCanopyProperty
!------------------------------------------------------------------------
  subroutine UpdateRootProperty(NZ,PathLen_pvr,FineRootRadius_rvr,AllRootC_vr,&
    FracPRoot4Uptake_pvr,MinFracPRoot4Uptake_pvr,FracPrimRootOccupiedLay_pvr,RootAreaDivRadius_vr)
!
!     update root characterization

  implicit none
  integer , intent(in) :: NZ
  real(r8), intent(in) :: AllRootC_vr(JZ1)
  real(r8), intent(out) :: PathLen_pvr(jroots,JZ1)
  real(r8), intent(out) :: FineRootRadius_rvr(jroots,JZ1)
  real(r8), intent(out) :: FracPRoot4Uptake_pvr(jroots,JZ1,JP1)
  real(r8), intent(out) :: MinFracPRoot4Uptake_pvr(jroots,JZ1,JP1)
  real(r8), intent(out) :: FracPrimRootOccupiedLay_pvr(JZ1,JP1)
  real(r8), intent(out) :: RootAreaDivRadius_vr(jroots,JZ1)
  real(r8) :: RootDepZ     !primary root depth
  real(r8) :: RTDPX        !Root occupied thickness in current layer
  integer :: N,L,NR

  associate(                                                      &
    PopuRootMycoC_pvr       => plt_biom% PopuRootMycoC_pvr,       &
    FracSoiAsMicP_vr        => plt_site%FracSoiAsMicP_vr,         &
    CumSoilThickness_vr     => plt_site%CumSoilThickness_vr,      & !in: depth from soil surface to layer bottom
    DLYR3                   => plt_site%DLYR3,                    &
    ZEROS                   => plt_site%ZEROS,                    &
    ZERO                    => plt_site%ZERO,                     &
    PlantPopulation_pft     => plt_site%PlantPopulation_pft,      &
    NU                      => plt_site%NU,                       &
    NumRootAxes_pft         => plt_morph%NumRootAxes_pft,         &
    MY                      => plt_morph%MY,                      &
    RootLenDensPerPlant_pvr => plt_morph%RootLenDensPerPlant_pvr, &
    Root2ndMaxRadius1_pft   => plt_morph%Root2ndMaxRadius1_pft,   &
    HypoctoHeight_pft       => plt_morph%HypoctoHeight_pft,       & !in: hypocotyledon height
    RootLenPerPlant_pvr     => plt_morph%RootLenPerPlant_pvr,     &
    Root1stDepz_pft         => plt_morph%Root1stDepz_pft,         &
    RootPorosity_pft        => plt_morph%RootPorosity_pft,        &
    RootVH2O_pvr            => plt_morph%RootVH2O_pvr,            &
    Root2ndMaxRadius_pft    => plt_morph%Root2ndMaxRadius_pft,    &
    SeedDepth_pft           => plt_morph%SeedDepth_pft,           & !in: seeding depth
    MaxSoiL4Root_pft        => plt_morph%MaxSoiL4Root_pft         &
  )
  
  D2000: DO N=1,MY(NZ)

    DO  L=NU,MaxSoiL4Root_pft(NZ)
      IF(N.EQ.ipltroot)THEN
        !obtain plant rooting depth
        RootDepZ=0.0_r8
        D2005: DO NR=1,NumRootAxes_pft(NZ)
          RootDepZ=AMAX1(RootDepZ,Root1stDepz_pft(ipltroot,NR,NZ))
        ENDDO D2005

        IF(L.EQ.NU)THEN
          FracPrimRootOccupiedLay_pvr(L,NZ)=1.0_r8
        ELSE
          IF(DLYR3(L).GT.ZERO)THEN
            RTDPX = AZMAX1(RootDepZ-CumSoilThickness_vr(L-1))
            RTDPX = AZMAX1(AMIN1(DLYR3(L),RTDPX)-AZMAX1(SeedDepth_pft(NZ)-CumSoilThickness_vr(L-1)-HypoctoHeight_pft(NZ)))
            FracPrimRootOccupiedLay_pvr(L,NZ) = RTDPX/DLYR3(L)
          ELSE
            FracPrimRootOccupiedLay_pvr(L,NZ)=0.0_r8
          ENDIF
        ENDIF

      ENDIF
      IF(AllRootC_vr(L).GT.ZEROS)THEN
        FracPRoot4Uptake_pvr(N,L,NZ)=PopuRootMycoC_pvr(N,L,NZ)/AllRootC_vr(L)
      ELSE
        FracPRoot4Uptake_pvr(N,L,NZ)=1.0_r8
      ENDIF
      MinFracPRoot4Uptake_pvr(N,L,NZ)=AMIN1(FMN,AMAX1(FMN,FracPRoot4Uptake_pvr(N,L,NZ)))

      IF(RootLenDensPerPlant_pvr(N,L,NZ).GT.ZERO .AND. FracPrimRootOccupiedLay_pvr(L,NZ).GT.ZERO)THEN
        FineRootRadius_rvr(N,L)=AMAX1(Root2ndMaxRadius1_pft(N,NZ),SQRT((RootVH2O_pvr(N,L,NZ) &
          /(1.0_r8-RootPorosity_pft(N,NZ)))/(PICON*PlantPopulation_pft(NZ)*RootLenPerPlant_pvr(N,L,NZ))))
        PathLen_pvr(N,L)=AMAX1(1.001_r8*FineRootRadius_rvr(N,L) &
          ,1.0_r8/(SQRT(PICON*(RootLenDensPerPlant_pvr(N,L,NZ)/FracPrimRootOccupiedLay_pvr(L,NZ))/FracSoiAsMicP_vr(L))))
        RootAreaDivRadius_vr(N,L)=TwoPiCON*RootLenPerPlant_pvr(N,L,NZ)/FracPrimRootOccupiedLay_pvr(L,NZ)
      ELSE
        FineRootRadius_rvr(N,L)   = Root2ndMaxRadius_pft(N,NZ)
        PathLen_pvr(N,L)          = 1.001_r8*FineRootRadius_rvr(N,L)
        RootAreaDivRadius_vr(N,L) = TwoPiCON*RootLenPerPlant_pvr(N,L,NZ)
      ENDIF
    enddo
  ENDDO D2000
  end associate
  end subroutine UpdateRootProperty
!------------------------------------------------------------------------

  subroutine HandlingDivergence(I,J,NN,NZ,TotalSoilPSIMPa_vr,DIFF,FDMP)

  implicit none
  integer  , intent(in) :: NN, I, J
  integer  , intent(in) :: NZ
  REAL(R8) , INTENT(IN) :: TotalSoilPSIMPa_vr(JZ1)
  real(r8) , intent(in) :: DIFF   !divergence check
  real(r8), intent(out):: FDMP
  real(r8) :: APSILT
  real(r8) :: CCPOLT
  real(r8) :: FTHRM
  real(r8) :: OSWT,Stomata_Stress

  integer :: N,L
! begin_execution
  associate(                                                         &
   ReistanceCanopy_pft       => plt_ew%ReistanceCanopy_pft,          &
   DeltaTKC_pft              => plt_ew%DeltaTKC_pft,                 &
   CanOsmoPsi0pt_pft         => plt_ew%CanOsmoPsi0pt_pft,            &
   PSICanopyOsmo_pft         => plt_ew%PSICanopyOsmo_pft,            &
   TKC_pft                   => plt_ew%TKC_pft,                      &
   TairK                     => plt_ew%TairK,                        &
   TKS_vr                    => plt_ew%TKS_vr,                       &
   PSIRootOSMO_vr            => plt_ew%PSIRootOSMO_vr,               &
   TdegCCanopy_pft        => plt_ew%TdegCCanopy_pft,           &
   AllPlantRootH2OLoss_vr  => plt_ew%AllPlantRootH2OLoss_vr,     &
   PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr,               &
   PSICanopyTurg_pft         => plt_ew%PSICanopyTurg_pft,            &
   PSIRoot_pvr               => plt_ew%PSIRoot_pvr,                  &
   VHeatCapCanopy_pft        => plt_ew%VHeatCapCanopy_pft,           &
   PSICanopy_pft             => plt_ew%PSICanopy_pft,                &
   NU                        => plt_site%NU,                         &
   AREA3                     => plt_site%AREA3,                      &
   iYearCurrent              => plt_site%iYearCurrent,               &
   MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,          &
   NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,        &
   MY                        => plt_morph%MY,                        &
   RootNonstructElmConc_rpvr => plt_biom%RootNonstructElmConc_rpvr,  &
   CanopyNonstElmConc_pft    => plt_biom%CanopyNonstElmConc_pft,     &
   ShootStrutElms_pft        => plt_biom%ShootStrutElms_pft,         &
   CanopyBndlResist_pft      => plt_photo%CanopyBndlResist_pft,      &
   MinCanPStomaResistH2O_pft => plt_photo%MinCanPStomaResistH2O_pft, &
   CanPStomaResistH2O_pft    => plt_photo%CanPStomaResistH2O_pft,    &
   H2OCuticleResist_pft      => plt_photo%H2OCuticleResist_pft,      &
   RCS                       => plt_photo%RCS,                       &
   FracPARads2Canopy_pft     => plt_rad%FracPARads2Canopy_pft,       &
   LWRadCanopy_pft           => plt_rad%LWRadCanopy_pft              &
  )
  IF(NN.GE.MaxIterNum)THEN
    WRITE(*,9999)iYearCurrent,I,J,NZ
9999  FORMAT('CONVERGENCE FOR WATER UPTAKE NOT ACHIEVED ON   ',6I4)

    IF(DIFF.GT.0.5_r8)THEN
      plt_rad%RadNet2Canopy_pft(NZ) = 0.0_r8
      plt_ew%EvapTransLHeat_pft(NZ)  = 0.0_r8
      plt_ew%HeatXAir2PCan_pft(NZ)  = 0.0_r8
      plt_ew%HeatStorCanopy_pft(NZ) = 0.0_r8
      plt_ew%VapXAir2Canopy_pft(NZ) = 0.0_r8
      plt_ew%Transpiration_pft(NZ)  = 0.0_r8
      TKC_pft(NZ)                   = TairK+DeltaTKC_pft(NZ)
      TdegCCanopy_pft(NZ)        = units%Kelvin2Celcius(TKC_pft(NZ))
      FTHRM                         = EMMC*stefboltz_const*FracPARads2Canopy_pft(NZ)*AREA3(NU)
      LWRadCanopy_pft(NZ)           = FTHRM*TKC_pft(NZ)**4._r8
      PSICanopy_pft(NZ)             = TotalSoilPSIMPa_vr(NGTopRootLayer_pft(NZ))
      CCPOLT                        = CanopyNonstElmConc_pft(ielmc,NZ)+CanopyNonstElmConc_pft(ielmn,NZ)&
        +CanopyNonstElmConc_pft(ielmp,NZ)

      CALL update_osmo_turg_pressure(PSICanopy_pft(NZ),CCPOLT,CanOsmoPsi0pt_pft(NZ),TKC_pft(NZ) &
        ,PSICanopyOsmo_pft(NZ),PSICanopyTurg_pft(NZ),FDMP)

      CanopyBndlResist_pft(NZ)   = ReistanceCanopy_pft(NZ)
      VHeatCapCanopy_pft(NZ)     = cpw*(ShootStrutElms_pft(ielmc,NZ)*10.0E-06_r8)
      DeltaTKC_pft(NZ)           = 0.0_r8
      Stomata_Stress             = EXP(RCS(NZ)*PSICanopyTurg_pft(NZ))
      CanPStomaResistH2O_pft(NZ) = MinCanPStomaResistH2O_pft(NZ) &
        +(H2OCuticleResist_pft(NZ)-MinCanPStomaResistH2O_pft(NZ))*Stomata_Stress

      D4290: DO N=1,MY(NZ)
        DO  L=NU,MaxSoiL4Root_pft(NZ)
          PSIRoot_pvr(N,L,NZ) = TotalSoilPSIMPa_vr(L)
          CCPOLT              = sum(RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ))
          CALL update_osmo_turg_pressure(PSIRoot_pvr(N,L,NZ),CCPOLT,CanOsmoPsi0pt_pft(NZ),TKS_vr(L),&
            PSIRootOSMO_vr(N,L,NZ),PSIRootTurg_vr(N,L,NZ))

          AllPlantRootH2OLoss_vr(N,L,NZ)=0.0_r8
      enddo
      ENDDO D4290
    ENDIF
  ENDIF
  end associate
  end subroutine HandlingDivergence
!------------------------------------------------------------------------------
  function CanopyEnergyH2OIteration(I,J,NZ,FracGrndByPFT,CanopyMassC,TotalSoilPSIMPa_vr,&
    VHeatCapCanopyAir,DIFF,cumPRootH2OUptake,CumPlantHeatLoss2Soil,HeatEvapSens,FDMP,&
    SoilRootResistance_rvr,FracPRoot4Uptake_pvr,AirMicPore4Fill_vr,WatAvail4Uptake_vr,TKCX,CNDT,&
    VHeatCapCanopyPrev_pft,PrecpHeatbyCanopy,PSIGravCanopyHeight,SoiLayerHasRoot_rvr) result(NN)

  implicit none
  integer  , intent(in) :: I, J
  integer  , intent(in) :: NZ
  real(r8) , intent(in) :: FracGrndByPFT,CanopyMassC
  real(r8) , intent(in) :: TotalSoilPSIMPa_vr(JZ1)
  real(r8) , intent(in) :: SoilRootResistance_rvr(jroots,JZ1)
  real(r8) , intent(in) :: FracPRoot4Uptake_pvr(jroots,JZ1,JP1)
  real(r8) , intent(in) :: AirMicPore4Fill_vr(JZ1),WatAvail4Uptake_vr(JZ1)
  real(r8) , intent(in) :: TKCX
  real(r8) , intent(in) :: CNDT
  real(r8) , intent(in) :: VHeatCapCanopyPrev_pft    !canopy heat capacity at previous time step, MJ/K
  real(r8) , intent(in) :: PrecpHeatbyCanopy         !heat added to canopy by precipitation [MJ]
  real(r8) , intent(in) :: PSIGravCanopyHeight       !Graviational water potential at effective canopy height [MPa]
  logical  , intent(in) :: SoiLayerHasRoot_rvr(jroots,JZ1)
  real(r8) , intent(out):: VHeatCapCanopyAir,DIFF
  real(r8) , intent(out):: cumPRootH2OUptake
  real(r8) , intent(out):: CumPlantHeatLoss2Soil  
  real(r8) , intent(out):: HeatEvapSens   !sensible heat due to evaporation/condensation  [MJ]
  real(r8) , intent(out):: FDMP
  real(r8) :: APSILT
  real(r8) :: CCPOLT
  real(r8) :: DTHS1
  real(r8) :: DIFFZ,DIFFU,DPSI
  real(r8) :: EX,PTransPre
  real(r8) :: FTHRM
  real(r8) :: LWRad2Canopy   !long wave radiation (from sky + ground) on canopy [MJ]
  real(r8) :: HeatAdd2Can
  real(r8) :: EffGrndAreaByPFT4H2O,EffGrndAreaByPFT4Heat,PSICanPPre,EvapConductCanopy
  real(r8) :: PSILC               !canopy height adjusted canopy water potential, [MPa]
  real(r8) :: RSSUX,RSSU,RA1,RSSZ
  real(r8) :: TKC1,TKCY
  real(r8) :: cumPRootH2OUptakePre
  real(r8) :: SymplasmicWatPrev,SymplasmicWat  !previous, current biomass-bounded water
  real(r8) :: VPC
  real(r8) :: XC,Stomata_Stress
  real(r8) :: RichardsNO
  integer  :: IC
  logical :: LIterationExit
  real(r8) :: DPSI_old
  character(len=64) :: fmt
!     return variables
  integer :: NN
!     local variables
  integer :: N,L
!     begin_execution

  associate(                                                          &
    NU                        => plt_site%NU,                         &
    AREA3                     => plt_site%AREA3,                      &
    QdewCanopy_pft            => plt_ew%QdewCanopy_pft,               &
    DeltaTKC_pft              => plt_ew%DeltaTKC_pft,                 &
    TdegCCanopy_pft        => plt_ew%TdegCCanopy_pft,           &
    PSICanopyOsmo_pft         => plt_ew%PSICanopyOsmo_pft,            &
    CanOsmoPsi0pt_pft         => plt_ew%CanOsmoPsi0pt_pft,            &
    TKC_pft                   => plt_ew%TKC_pft,                      &
    ReistanceCanopy_pft       => plt_ew%ReistanceCanopy_pft,          &
    VPA                       => plt_ew%VPA,                          &  !air vapor pressure
    TKS_vr                    => plt_ew%TKS_vr,                       &
    TairK                     => plt_ew%TairK,                        &
    RIB                       => plt_ew%RIB,                          &
    Transpiration_pft         => plt_ew%Transpiration_pft,            &  !ton H2O/d2
    PSICanopy_pft             => plt_ew%PSICanopy_pft,                &  !canopy total water potential [MPa]
    CanopyBiomWater_pft       => plt_ew%CanopyBiomWater_pft,          &  !water held in canopy biomass [m3 d-2]
    AllPlantRootH2OLoss_vr    => plt_ew%AllPlantRootH2OLoss_vr,       &
    TKCanopy_pft              => plt_ew%TKCanopy_pft,                 &
    VapXAir2Canopy_pft        => plt_ew%VapXAir2Canopy_pft,           &
    WatHeldOnCanopy_pft       => plt_ew%WatHeldOnCanopy_pft,          &
    PSICanopyTurg_pft         => plt_ew%PSICanopyTurg_pft,            &
    VHeatCapCanopy_pft        => plt_ew%VHeatCapCanopy_pft,           &
    PrecIntcptByCanopy_pft    => plt_ew%PrecIntcptByCanopy_pft,       &
    EvapTransLHeat_pft        => plt_ew%EvapTransLHeat_pft,           &  !out: latent heat due to evaporation/condensation [MJ]
    ZERO4LeafVar_pft          => plt_biom%ZERO4LeafVar_pft,           &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,             &
    CanopyNonstElmConc_pft    => plt_biom%CanopyNonstElmConc_pft,     &
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,          &
    MY                        => plt_morph%MY,                        &
    MinCanPStomaResistH2O_pft => plt_photo%MinCanPStomaResistH2O_pft, &
    RCS                       => plt_photo%RCS,                       &
    CanopyBndlResist_pft      => plt_photo%CanopyBndlResist_pft,      &
    H2OCuticleResist_pft      => plt_photo%H2OCuticleResist_pft,      &
    CanPStomaResistH2O_pft    => plt_photo%CanPStomaResistH2O_pft,    &
    RadNet2Canopy_pft         => plt_rad%RadNet2Canopy_pft,           &
    RadSWbyCanopy_pft         => plt_rad%RadSWbyCanopy_pft,           &
    LWRadSky_col              => plt_rad%LWRadSky_col,                &
    FracPARads2Canopy_pft     => plt_rad%FracPARads2Canopy_pft,       &
    LWRadCanopy_pft           => plt_rad%LWRadCanopy_pft,             &
    LWRadGrnd                 => plt_rad%LWRadGrnd                    &
  )
  
  !CCPOLT: total nonstructural canopy C,N,P concentration
  !FTHRM:coefficient for LW emitted by canopy  
  !LWRad2Canopy:long-wave absorbed by canopy
  FTHRM        = EMMC*stefboltz_const*FracPARads2Canopy_pft(NZ)*AREA3(NU)
  LWRad2Canopy = (LWRadSky_col+LWRadGrnd)*FracPARads2Canopy_pft(NZ)
  CCPOLT       = CanopyNonstElmConc_pft(ielmc,NZ)+CanopyNonstElmConc_pft(ielmn,NZ) &
    +CanopyNonstElmConc_pft(ielmp,NZ)
! ReistanceCanopy_pft          = canopy isothermal boundary later resistance

  cumPRootH2OUptake     = 0.0_r8
  EffGrndAreaByPFT4H2O  = FracGrndByPFT*AREA3(NU)                    !area coverd by pft m2
  EffGrndAreaByPFT4Heat = FracGrndByPFT*AREA3(NU)*1.25E-03_r8        !1.25e-3 is heat capacity of air [MJ/(m3 K)], assuming air volume is not affected by plant mass
  RA1                   = ReistanceCanopy_pft(NZ)

  IC                    = 0
  XC                    = 0.5_r8   !learning/updating rate for canopy temperature, not used now
  LIterationExit        = .false.
  PSICanPPre            = 0.0_r8
  PTransPre             = 0.0_r8
  cumPRootH2OUptakePre  = 0.0_r8
  SymplasmicWatPrev     = 0.0_r8
  DPSI_old              = 0._r8
  CumPlantHeatLoss2Soil = 0._r8

  D4000: DO NN=1,MaxIterNum
!
!     NET RADIATION FROM ABSORBED SW AND NET LW
!
!     LWRadCanopy_pft=LW emitted by canopy
!     DTHS1=net LW absorbed by canopy
!     RadSWbyCanopy_pft=total SW absorbed by canopy
!     RadNet2Canopy_pft=net SW+LW absorbed by canopy
!
    TKC1                  = TKCanopy_pft(NZ)
    LWRadCanopy_pft(NZ)   = FTHRM*TKC1**4._r8              !long wave radiation
    DTHS1                 = LWRad2Canopy-LWRadCanopy_pft(NZ)*2.0_r8        !net long wave radiation to canopy
    RadNet2Canopy_pft(NZ) = RadSWbyCanopy_pft(NZ)+DTHS1  !total radiation to canopy
!
!     BOUNDARY LAYER RESISTANCE FROM RICHARDSON NUMBER
!
!     RI=Ricardson's number
!     RA=canopy boundary layer resistance
!     EvapConductCanopy,VHeatCapCanopyAir=canopy latent,sensible heat conductance
!
    RichardsNO               = RichardsonNumber(RIB,TairK,TKC1)
    CanopyBndlResist_pft(NZ) = AMAX1(MinCanopyBndlResist_pft,0.9_r8*RA1 &
      ,AMIN1(1.1_r8*RA1,ReistanceCanopy_pft(NZ)/(1.0_r8-10.0_r8*RichardsNO)))

    RA1               = CanopyBndlResist_pft(NZ)
    EvapConductCanopy = EffGrndAreaByPFT4H2O/CanopyBndlResist_pft(NZ)  !m2/(h/m)                 = m3/h
    VHeatCapCanopyAir = EffGrndAreaByPFT4Heat/CanopyBndlResist_pft(NZ) !canopy air heat capacity,
!
!     CANOPY WATER AND OSMOTIC POTENTIALS
!
!     CanOsmoPsi0pt_pft=osmotic potential at PSICanopy_pft=0 from PFT file
!     PSICanopyOsmo_pft,PSICanopyTurg_pft=canopy osmotic,turgor water potential
!
     call update_osmo_turg_pressure(PSICanopy_pft(NZ),CCPOLT,CanOsmoPsi0pt_pft(NZ),TKC1 &
       ,PSICanopyOsmo_pft(NZ),PSICanopyTurg_pft(NZ),FDMP)
!
!     CANOPY STOMATAL RESISTANCE
!
!     RCS=shape parameter for CanPStomaResistH2O_pftvs PSICanopyTurg_pft from PFT file
!     RC=canopy stomatal resistance
!     MinCanPStomaResistH2O_pft=minimum CanPStomaResistH2O_pftat PSICanopy_pft=0 from stomate.f
!     CuticleResist_pft=cuticular resistance from PFT file
!
    Stomata_Stress             = EXP(RCS(NZ)*PSICanopyTurg_pft(NZ))
    CanPStomaResistH2O_pft(NZ) = MinCanPStomaResistH2O_pft(NZ)+Stomata_Stress &
      *(H2OCuticleResist_pft(NZ)-MinCanPStomaResistH2O_pft(NZ))
!
!     CANOPY VAPOR PRESSURE AND EVAPORATION OF INTERCEPTED WATER
!     OR TRANSPIRATION OF UPTAKEN WATER
!
!     VPC,VPA=vapor pressure inside canopy, in atmosphere
!     TKC1=canopy temperature
!     EX=canopy-atmosphere water flux
!     RA,RZ=canopy boundary layer,surface resistance
!     VapXAir2Canopy_pft=water flux to,from air to canopy surfaces, [ton H2O/(h*m2)]
!     WatHeldOnCanopy_pft=water volume on canopy surfaces
!     EP=water flux to,from inside canopy
!     EvapTransLHeat_pft=canopy latent heat flux
!     HeatEvapSens=convective heat flux from EvapTransLHeat_pft, <0, into air
!     VAP=latent heat of evaporation
!     EvapConductCanopy=aerodynamic conductance

    VPC = vapsat(tkc1)*EXP(18.0_r8*AMAX1(PSICanopy_pft(NZ),-5000._r8)/(RGASC*TKC1))  !ton H2O/m3
    EX  = EvapConductCanopy*(VPA-VPC)   !air to canopy water vap flux, [m/h]*[ton/m3] = [ton H2O/(h*m2)]
    
    !Dew condensation to canopy >0 to canopy
    IF(EX.GT.0.0_r8)THEN     
      VapXAir2Canopy_pft(NZ) = EX*CanopyBndlResist_pft(NZ)/(CanopyBndlResist_pft(NZ)+RZ)
      QdewCanopy_pft(NZ)     = VapXAir2Canopy_pft(NZ)
      EX                     = 0.0_r8
      HeatEvapSens           = VapXAir2Canopy_pft(NZ)*cpw*TairK               !enthalpy of condensed water add to canopy, MJ/(h)
    !canopy lose water, and there is canopy-held water  
    ELSEIF(EX.LE.0.0_r8)THEN 
      if(WatHeldOnCanopy_pft(NZ).GT.0.0_r8)THEN
        !evaporation, and there is water stored in canopy
        !VapXAir2Canopy_pft <0._r8, off canopy as evaporation, cannot be more than WatHeldOnCanopy_pft(NZ)
        VapXAir2Canopy_pft(NZ) = AMAX1(EX*CanopyBndlResist_pft(NZ)/(CanopyBndlResist_pft(NZ)+RZ),-WatHeldOnCanopy_pft(NZ))
        EX                     = EX-VapXAir2Canopy_pft(NZ)                !demand for transpiration
        HeatEvapSens           = VapXAir2Canopy_pft(NZ)*cpw*TKC1          !enthalpy of evaporated water leaving canopy
      ELSE
        VapXAir2Canopy_pft(NZ) =0._r8
        HeatEvapSens = 0._r8  
      ENDIF
    ENDIF

    !Transpiration_pft<0 means canopy lose water through transpiration
    !EvapTransLHeat_pft:latent heat flux, negative means into atmosphere

    Transpiration_pft(NZ)  = EX*CanopyBndlResist_pft(NZ)/(CanopyBndlResist_pft(NZ)+CanPStomaResistH2O_pft(NZ))
    EvapTransLHeat_pft(NZ) = (Transpiration_pft(NZ)+VapXAir2Canopy_pft(NZ))*EvapLHTC
    HeatEvapSens           = HeatEvapSens + Transpiration_pft(NZ)*cpw*TKC1
!
!     SENSIBLE + STORAGE HEAT FROM RN, LE AND CONVECTIVE HEAT FLUXES
!
!     HeatSenAddStore=initial estimate of sensible+storage heat flux
!     PrecpHeatbyCanopy=convective heat flux from precip to canopy
!
    HeatAdd2Can=RadNet2Canopy_pft(NZ)+EvapTransLHeat_pft(NZ)+HeatEvapSens+PrecpHeatbyCanopy-CumPlantHeatLoss2Soil
!
!     SOLVE FOR CANOPY TEMPERATURE CAUSED BY SENSIBLE + STORAGE HEAT
!
!     VHeatCapCanopy_pft=canopy heat capacity
!     TKCY=equilibrium canopy temperature for HeatSenAddStore
!
!   VHeatCapCanopyPrev_pft= canopy heat capacity, including biomass, held and bounded water
!   assuming transpiration does not change water binded by the canopy, i.e. CanopyBiomWater_pft(NZ)
!   canopy bound water is changing
    VHeatCapCanopy_pft(NZ)=VHeatCapCanopyPrev_pft+cpw*(PrecIntcptByCanopy_pft(NZ)+VapXAir2Canopy_pft(NZ) &
       +Transpiration_pft(NZ)-cumPRootH2OUptake)
!   new canopy temperature 
    TKCY=(TKCX*VHeatCapCanopyPrev_pft+HeatAdd2Can+TairK*VHeatCapCanopyAir)/(VHeatCapCanopy_pft(NZ)+VHeatCapCanopyAir)

    !limit canopy temperature to no different from air for more than 10 K
    TKCY=AMIN1(TairK+10.0_r8,AMAX1(TairK-10.0_r8,TKCY))
!
!     RESET CANOPY TEMPERATURE FOR NEXT ITERATION
!
!     XC,IC=magnitude,direction of change in canopy temp for next cycle
!
    IF((IC.EQ.0 .AND. TKCY.GT.TKC1) .OR. (IC.EQ.1 .AND. TKCY.LT.TKC1))THEN
      XC=0.5_r8*XC
    ENDIF
    !0.1 is the learning/updating rate
    TKCanopy_pft(NZ)=TKC1+0.1_r8*(TKCY-TKC1)

    IF(TKCY.GT.TKC1)THEN
      !warming
      IC=1
    ELSE
      !cooling or no change
      IC=0
    ENDIF

!
!     IF CONVERGENCE CRITERION IS MET OR ON EVERY TENTH ITERATION,
!     PROCEED TO WATER BALANCE
!
!
    IF(ABS(TKCY-TKC1).LT.0.05_r8 .OR. (NN/10)*10.EQ.NN)THEN
      cumPRootH2OUptake = 0.0_r8
      PSILC             = PSICanopy_pft(NZ)+PSIGravCanopyHeight
!
!     ROOT WATER UPTAKE FROM SOIL-CANOPY WATER POTENTIALS,
!     SOIL + ROOT HYDRAULIC RESISTANCES
!
!     SoiLayerHasRoot_rvr=rooted layer flag
!     AllPlantRootH2OLoss_vr=root water uptake from soil layer > 0
!     WatAvail4Uptake_vr,AirMicPore4Fill_vr=water volume available for uptake,air volume
!     FracPRoot4Uptake_pvr=PFT fraction of biome root mass
!     PSILC=height corrected canopy water potential 
!     SoilRootResistance_rvr=total soil+root resistance
!     cumPRootH2OUptake=total root water uptake from soil 
!     CumPlantHeatLoss2Soil (<0.), add heat to canopy 

      CumPlantHeatLoss2Soil=0._r8
      D4200: DO N=1,MY(NZ)
        D4201: DO L=NU,MaxSoiL4Root_pft(NZ)
          IF(SoiLayerHasRoot_rvr(N,L))THEN
            !<0 active uptake
            AllPlantRootH2OLoss_vr(N,L,NZ)=AMAX1(AZMIN1(-WatAvail4Uptake_vr(L)*FracPRoot4Uptake_pvr(N,L,NZ)), &
              AMIN1((PSILC-TotalSoilPSIMPa_vr(L))/SoilRootResistance_rvr(N,L), &
              AirMicPore4Fill_vr(L)*FracPRoot4Uptake_pvr(N,L,NZ)))

            !plant/myco lose water to soil > 0
            IF(AllPlantRootH2OLoss_vr(N,L,NZ).GT.0.0_r8)THEN              
              !why multiply 0.1 here? I don't know
              AllPlantRootH2OLoss_vr(N,L,NZ)=0.1_r8*AllPlantRootH2OLoss_vr(N,L,NZ)

              !plant moves heat from canopy to soil,
              CumPlantHeatLoss2Soil=CumPlantHeatLoss2Soil+cpw*AllPlantRootH2OLoss_vr(N,L,NZ)*TKC1

            !plant/myco gains water from soil < 0
            else  
              CumPlantHeatLoss2Soil=CumPlantHeatLoss2Soil+cpw*AllPlantRootH2OLoss_vr(N,L,NZ)*TKS_vr(L)
            ENDIF
            cumPRootH2OUptake=cumPRootH2OUptake+AllPlantRootH2OLoss_vr(N,L,NZ)
          ELSE
            AllPlantRootH2OLoss_vr(N,L,NZ)=0.0_r8
          ENDIF
        enddo D4201
      ENDDO D4200
      !turn it off at the moment
!
!     TEST TRANSPIRATION - ROOT WATER UPTAKE VS. CHANGE IN CANOPY
!     WATER STORAGE
!
!     SymplasmicWat  = total water can be held in biomass
!     CanopyBiomWater_pft= canopy water content held in biomass
!     DIFFZ= extra canopy water binding capacity
!     DIFFU= change to canopy binded water
!     DIFFU-DIFFZ=residual of canopy binded water,
!     >0 gain binded water (increase leaf P, more positive), <0 lose binded water (more uptake, decrease leaf water P) 
!     DIFF=normalized difference between DIFFZ and DIFFU
!     5.0E-03=acceptance criterion for DIFF
!     RSSZ=change in canopy water potl vs change in canopy water cnt
!     RSSU=change in canopy water potl vs change in transpiration
      
      SymplasmicWat = ppmc*CanopyMassC/FDMP              
      DIFFZ         = SymplasmicWat-CanopyBiomWater_pft(NZ)    !biomass water deficit
      DIFFU         = Transpiration_pft(NZ)-cumPRootH2OUptake  !root water uptake excess

      !ideally, the difference between DIFFZ and DIFFU should be as small as possible      
      IF(.not.isclose(cumPRootH2OUptake,0.0_r8))THEN
        DIFF=ABS((DIFFU-DIFFZ)/cumPRootH2OUptake)
      ELSE
        DIFF=ABS(DIFFU-DIFFZ)/SymplasmicWat
      ENDIF

      !the relative difference is small enough
      IF(DIFF.LT.5.0E-03_r8)THEN
        IF(LIterationExit)EXIT
        LIterationExit=.true.
        CALL StomatalDynamics(I,J,NZ)
        CYCLE
      ENDIF
      
      IF(ABS(SymplasmicWat-SymplasmicWatPrev).GT.ZERO4Groth_pft(NZ))THEN
        RSSZ=ABS((PSICanopy_pft(NZ)-PSICanPPre)/(SymplasmicWat-SymplasmicWatPrev))
      ELSEIF(CNDT.GT.ZERO4Groth_pft(NZ))THEN
        RSSZ=1.0_r8/CNDT   !resistance
      ELSE
        RSSZ=ZERO4LeafVar_pft(NZ)
      ENDIF
      
      IF(ABS(Transpiration_pft(NZ)-PTransPre).GT.ZERO4Groth_pft(NZ))THEN
        RSSUX=ABS((PSICanopy_pft(NZ)-PSICanPPre)/(Transpiration_pft(NZ)-PTransPre))
        IF(CNDT.GT.ZERO4Groth_pft(NZ))THEN
          RSSU=AMIN1(1.0_r8/CNDT,RSSUX)  !resistance for uptake
        ELSE
          RSSU=RSSUX
        ENDIF
      ELSEIF(ABS(cumPRootH2OUptake-cumPRootH2OUptakePre).GT.ZERO4Groth_pft(NZ))THEN
        RSSUX=ABS((PSICanopy_pft(NZ)-PSICanPPre)/(cumPRootH2OUptake-cumPRootH2OUptakePre))
        IF(CNDT.GT.ZERO4Groth_pft(NZ))THEN
          RSSU=AMIN1(1.0_r8/CNDT,RSSUX)  !resistance
        ELSE
          RSSU=RSSUX
        ENDIF
      ELSEIF(CNDT.GT.ZERO4Groth_pft(NZ))THEN
        RSSU=1.0_r8/CNDT
      ELSE
        RSSU=ZERO4LeafVar_pft(NZ)
      ENDIF
!
!     CHANGE IN CANOPY WATER POTENTIAL REQUIRED TO BRING AGREEMENT
!     BETWEEN TRANSPIRATION - ROOT WATER UPTAKE AND CHANGE IN CANOPY
!     WATER STORAGE
!
!     DPSI=change in PSICanopy_pft(for next convergence cycle
!     1.0E-03=acceptance criterion for DPSI
!
      DPSI=AMIN1(AMIN1(RSSZ,RSSU)*(DIFFU-DIFFZ),ABS(PSICanopy_pft(NZ)))

!     IF CONVERGENCE CRITERION IS MET THEN FINISH,
!     OTHERWISE START NEXT ITERATION WITH CANOPY WATER POTENTIAL
!     TRANSPIRATION, UPTAKE AND WATER CONTENT FROM CURRENT ITERATION
!
      IF((NN.GE.30 .AND. ABS(DPSI).LT.1.0E-03_r8) .OR. NN.GE.MaxIterNum)then
        IF(LIterationExit)EXIT
        LIterationExit=.true.
        CALL StomatalDynamics(I,J,NZ)      
      ELSE
        !prepare for next iteration
        PSICanPPre           = PSICanopy_pft(NZ)
        PTransPre            = Transpiration_pft(NZ)
        cumPRootH2OUptakePre = cumPRootH2OUptake
        SymplasmicWatPrev    = SymplasmicWat

        PSICanopy_pft(NZ) = AZMIN1(PSICanopy_pft(NZ)+0.1_r8*DPSI)
        DPSI_old          = DPSI
        XC                = 0.50_r8!
      ENDIF
    ENDIF
  ENDDO D4000

!
!     FINAL CANOPY TEMPERATURE, DIFFERENCE WITH AIR TEMPERATURE
!
!     TKC=final estimate of canopy temperature TKCanopy_pft
!     TairK=current air temperature
!     DeltaTKC_pft=TKC-TairK for next hour
!
  TKC_pft(NZ)            = TKCanopy_pft(NZ)
  TdegCCanopy_pft(NZ) = units%Kelvin2Celcius(TKC_pft(NZ))
  DeltaTKC_pft(NZ)       = TKC_pft(NZ)-TairK

  end associate
  end function CanopyEnergyH2OIteration
!------------------------------------------------------------------------
  subroutine CalcResistance(NZ,PathLen_pvr,FineRootRadius_rvr,RootAreaDivRadius_vr,&
    RootResist,RootResistSoi,RootResistPrimary,RootResist2ndary,SoiH2OResist,&
    SoilRootResistance_rvr,CNDT,PSIGravCanopyHeight,SoiLayerHasRoot_rvr)

  implicit none
  integer, intent(in)   :: NZ
  real(r8), intent(in)  :: PathLen_pvr(jroots,JZ1),FineRootRadius_rvr(jroots,JZ1)
  real(r8), intent(in)  :: RootAreaDivRadius_vr(jroots,JZ1)
  real(r8), intent(out) :: RootResist(jroots,JZ1),RootResistSoi(jroots,JZ1)
  real(r8), intent(out) :: RootResistPrimary(jroots,JZ1),RootResist2ndary(jroots,JZ1)
  real(r8), intent(out) :: SoiH2OResist(jroots,JZ1)
  real(r8), intent(out) :: SoilRootResistance_rvr(jroots,JZ1)   !added soil and root resistance
  real(r8), intent(out) :: CNDT                   !total root conductance
  real(r8), intent(out) :: PSIGravCanopyHeight    !gravimetric water potential at CanopyHeight4WatUptake_pft, [MPa]
  logical , intent(out) :: SoiLayerHasRoot_rvr(jroots,JZ1)

  character(len=*), parameter :: subname='CalcResistance'
  real(r8) :: FRADW,FRAD1,FRAD2
  real(r8) :: RSSL,RTAR2
  integer :: N, L
  associate(                                                                 &
    CumSoilThickMidL_vr         => plt_site%CumSoilThickMidL_vr,             &
    PlantPopulation_pft         => plt_site%PlantPopulation_pft,             &
    VLWatMicPM_vr               => plt_site%VLWatMicPM_vr,                   &
    ZERO                        => plt_site%ZERO,                            &
    ZEROS2                      => plt_site%ZEROS2,                          &
    NU                          => plt_site%NU,                              &
    PSICanopy_pft               => plt_ew%PSICanopy_pft,                     &
    ZERO4Groth_pft              => plt_biom%ZERO4Groth_pft,                  &
    THETW_vr                    => plt_soilchem%THETW_vr,                    &
    VLMicP_vr                   => plt_soilchem%VLMicP_vr,                   &
    HydroCondMicP4RootUptake_vr => plt_soilchem%HydroCondMicP4RootUptake_vr, &
    VLSoilPoreMicP_vr           => plt_soilchem%VLSoilPoreMicP_vr,           &
    Root2ndXNum_pvr             => plt_morph%Root2ndXNum_pvr,                &
    Root1stXNumL_pvr            => plt_morph%Root1stXNumL_pvr,               &
    Root2ndMaxRadius_pft        => plt_morph%Root2ndMaxRadius_pft,           &
    RootAxialResist_pft         => plt_morph%RootAxialResist_pft,            &
    Root2ndRadius_pvr           => plt_morph%Root2ndRadius_pvr,              &
    RootLenDensPerPlant_pvr     => plt_morph%RootLenDensPerPlant_pvr,        &
    CanopyHeight4WatUptake_pft  => plt_morph%CanopyHeight4WatUptake_pft,     &  !effective canopy height for water uptake calculation
    RootLenPerPlant_pvr         => plt_morph%RootLenPerPlant_pvr,            &
    Root2ndAveLen_pvr           => plt_morph%Root2ndAveLen_pvr,              &
    Root1stRadius_pvr           => plt_morph%Root1stRadius_pvr,              &
    MY                          => plt_morph%MY,                             &
    RootRadialResist_pft        => plt_morph%RootRadialResist_pft,           &
    CanopyHeight_pft            => plt_morph%CanopyHeight_pft,               &  !canopy height [m]
    NGTopRootLayer_pft          => plt_morph%NGTopRootLayer_pft,             &
    MaxSoiL4Root_pft            => plt_morph%MaxSoiL4Root_pft                &
  )

  !     GRAVIMETRIC WATER POTENTIAL FROM CANOPY HEIGHT
  !
  !     CanopyHeight4WatUptake_pft=canopy height for water uptake
  !     FRADW=conducting elements of stalk relative to those of primary root
  !     PSICanopy_pft=canopy total water potential
  !     EMODW=wood modulus of elasticity (MPa)
  call PrintInfo('beg '//subname)

  CNDT                           = 0.0_r8
  CanopyHeight4WatUptake_pft(NZ) = 0.80_r8*CanopyHeight_pft(NZ)
  PSIGravCanopyHeight            = mGravAccelerat*CanopyHeight4WatUptake_pft(NZ)
  FRADW                          = 1.0E+04_r8*(AMAX1(0.5_r8,1.0_r8+PSICanopy_pft(NZ)/EMODW))**4._r8
!
  !     SOIL AND ROOT HYDRAULIC RESISTANCES TO ROOT WATER UPTAKE
  !
  !      VLSoilPoreMicP_vr,VLWatMicPM,THETW=soil,water volume,content
  !     RootLenDensPerPlant_pvr,RootLenPerPlant_pvr=root length density,root length per plant
  !     HydroCondMicP4RootUptake_vr=soil hydraulic conductivity for root uptake
  !     Root1stXNumL_pvr,Root2ndXNum_pvr=number of root,myco primary,secondary axes
  !     SoiLayerHasRoot_rvr:1=rooted,0=not rooted
  !     N:1=root,2=mycorrhizae
!
  D3880: DO N=1,MY(NZ)
    DO  L=NU,MaxSoiL4Root_pft(NZ)

      SoiLayerHasRoot_rvr(N,L)=VLSoilPoreMicP_vr(L).GT.ZEROS2        &  
        .AND. VLWatMicPM_vr(NPH,L).GT.ZEROS2                         &
        .AND. RootLenDensPerPlant_pvr(N,L,NZ).GT.ZERO                &
        .AND. HydroCondMicP4RootUptake_vr(L).GT.ZERO                 &
        .AND. Root1stXNumL_pvr(ipltroot,L,NZ).GT.ZERO4Groth_pft(NZ)  &
        .AND. Root2ndXNum_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ)          &
        .AND. THETW_vr(L).GT.ZERO
      if(SoiLayerHasRoot_rvr(N,L))THEN
        
        !
        !     SOIL HYDRAULIC RESISTANCE FROM RADIAL UPTAKE GEOMETRY
        !     AND SOIL HYDRAULIC CONDUCTIVITY
        !
        !     SoiH2OResist=soil hydraulic resistance
        !     PP=plant population
        !     PathLen_pvr=path length of water and nutrient uptake
        !     FineRootRadius_rvr,RootAreaDivRadius_vr=root radius,surface/radius area
        !
        RSSL              = (LOG(PathLen_pvr(N,L)/FineRootRadius_rvr(N,L))/RootAreaDivRadius_vr(N,L))/PlantPopulation_pft(NZ)
        SoiH2OResist(N,L) = RSSL/HydroCondMicP4RootUptake_vr(L)
        !
        !     RADIAL ROOT RESISTANCE FROM ROOT AREA AND RADIAL RESISTIVITY
        !     ENTERED IN 'READQ'
        !
        !     Root2ndRadius_pvr=secondary root radius
        !     RootLenPerPlant_pvr=root length per plant
        !     RootResistSoi=radial resistance
        !     RootRadialResist_pft=radial resistivity from PFT file
        !     VLMicP,VLWatMicPM=soil micropore,water volume
        !
        RTAR2              = TwoPiCON*Root2ndRadius_pvr(N,L,NZ)*RootLenPerPlant_pvr(N,L,NZ)*PlantPopulation_pft(NZ)
        RootResistSoi(N,L) = RootRadialResist_pft(N,NZ)/RTAR2*VLMicP_vr(L)/VLWatMicPM_vr(NPH,L)
!
        !     ROOT AXIAL RESISTANCE FROM RADII AND LENGTHS OF PRIMARY AND
        !     SECONDARY ROOTS AND FROM AXIAL RESISTIVITY ENTERED IN 'READQ'
        !
        !     FRAD1,FRAD2=primary,secondary root radius relative to maximum
        !     secondary radius from PFT file Root2ndMaxRadius_pft at which RootAxialResist_pft is defined
        !     Root1stRadius_pvr,Root2ndRadius_pvr=primary,secondary root radius
        !     RootAxialResist_pft=axial resistivity from PFT file
        !     DPTHZ=depth of primary root from surface
        !     RootResistPrimary,RootResist2ndary=axial resistance of primary,secondary roots
        !     Root2ndAveLen_pvr=average secondary root length
        !     Root1stXNumL_pvr,Root2ndXNum_pvr=number of primary,secondary axes
!
        FRAD1                  = (Root1stRadius_pvr(N,L,NZ)/Root2ndMaxRadius_pft(N,NZ))**4._r8
        RootResistPrimary(N,L) = RootAxialResist_pft(N,NZ)*CumSoilThickMidL_vr(L)/(FRAD1*Root1stXNumL_pvr(ipltroot,L,NZ)) &
          +RootAxialResist_pft(ipltroot,NZ)*CanopyHeight4WatUptake_pft(NZ)/(FRADW*Root1stXNumL_pvr(ipltroot,L,NZ))

        FRAD2                = (Root2ndRadius_pvr(N,L,NZ)/Root2ndMaxRadius_pft(N,NZ))**4._r8
        RootResist2ndary(N,L) = RootAxialResist_pft(N,NZ)*Root2ndAveLen_pvr(N,L,NZ)/(FRAD2*Root2ndXNum_pvr(N,L,NZ))
        !
        !     TOTAL ROOT RESISTANCE = SOIL + RADIAL + AXIAL
        !
        !     RootResist=root radial+axial resistance
        !     SoilRootResistance_rvr=total soil+root resistance
        !     CNDT=total soil+root conductance for all layers
        ! assuming all roots work in parallel

        RootResist(N,L)       = RootResistSoi(N,L)+RootResistPrimary(N,L)+RootResist2ndary(N,L)
        SoilRootResistance_rvr(N,L) = SoiH2OResist(N,L)+RootResist(N,L)
        CNDT                  = CNDT+1.0_r8/SoilRootResistance_rvr(N,L)

      ENDIF
    enddo
  ENDDO D3880
  call PrintInfo('end '//subname)
  end associate
  end subroutine CalcResistance
!------------------------------------------------------------------------

  subroutine HandleBareSoil(NZ,TotalSoilPSIMPa_vr,FDMP)
  !
  !plant has not emerged yet
  use data_const_mod, only : spval => DAT_CONST_SPVAL
  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: TotalSoilPSIMPa_vr(JZ1)
  real(r8), intent(out):: FDMP
  integer :: N,L
  real(r8) :: APSILT
  real(r8) :: CCPOLT
  real(r8) :: FTHRM
  real(r8) :: OSWT,Stomata_Stress

! begin_execution
  associate(                                                          &
    QdewCanopy_pft            => plt_ew%QdewCanopy_pft,               &
    TdegCCanopy_pft        => plt_ew%TdegCCanopy_pft,           &
    CanOsmoPsi0pt_pft         => plt_ew%CanOsmoPsi0pt_pft,            &
    TKSnow                    => plt_ew%TKSnow,                       &
    TKC_pft                   => plt_ew%TKC_pft,                      &
    TKS_vr                    => plt_ew%TKS_vr,                       &
    TairK                     => plt_ew%TairK,                        &
    SnowDepth                 => plt_ew%SnowDepth,                    &
    DeltaTKC_pft              => plt_ew%DeltaTKC_pft,                 &
    ReistanceCanopy_pft       => plt_ew%ReistanceCanopy_pft,          &
    VHeatCapCanopy_pft        => plt_ew%VHeatCapCanopy_pft,           &
    PSICanopyOsmo_pft         => plt_ew%PSICanopyOsmo_pft,            &
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr,               &
    AllPlantRootH2OLoss_vr  => plt_ew%AllPlantRootH2OLoss_vr,     &
    PSICanopyTurg_pft         => plt_ew%PSICanopyTurg_pft,            &
    PSICanopy_pft             => plt_ew%PSICanopy_pft,                &
    PSIRootOSMO_vr            => plt_ew%PSIRootOSMO_vr,               &
    PSIRoot_pvr               => plt_ew%PSIRoot_pvr,                  &
    Transpiration_pft         => plt_ew%Transpiration_pft,            &
    HeatStorCanopy_pft        => plt_ew%HeatStorCanopy_pft,           &
    EvapTransLHeat_pft        => plt_ew%EvapTransLHeat_pft,           &
    HeatXAir2PCan_pft         => plt_ew%HeatXAir2PCan_pft,            &
    VapXAir2Canopy_pft        => plt_ew%VapXAir2Canopy_pft,           &
    NU                        => plt_site%NU,                         &
    ZERO                      => plt_site%ZERO,                       &
    AREA3                     => plt_site%AREA3,                      &
    CanopyNonstElmConc_pft    => plt_biom%CanopyNonstElmConc_pft,     &
    RootNonstructElmConc_rpvr => plt_biom%RootNonstructElmConc_rpvr,  &
    ShootStrutElms_pft        => plt_biom%ShootStrutElms_pft,         &
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,          &
    CanopyHeight_pft          => plt_morph%CanopyHeight_pft,          &
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft,        &
    MY                        => plt_morph%MY,                        &
    RCS                       => plt_photo%RCS,                       &
    CanopyBndlResist_pft      => plt_photo%CanopyBndlResist_pft,      &
    CanPStomaResistH2O_pft    => plt_photo%CanPStomaResistH2O_pft,    &
    MinCanPStomaResistH2O_pft => plt_photo%MinCanPStomaResistH2O_pft, &
    H2OCuticleResist_pft      => plt_photo%H2OCuticleResist_pft,      &
    FracPARads2Canopy_pft     => plt_rad%FracPARads2Canopy_pft,       &
    LWRadCanopy_pft           => plt_rad%LWRadCanopy_pft,             &
    RadNet2Canopy_pft         => plt_rad%RadNet2Canopy_pft            &
  )
  RadNet2Canopy_pft(NZ)  = 0.0_r8
  EvapTransLHeat_pft(NZ)  = 0.0_r8
  HeatXAir2PCan_pft(NZ)  = 0.0_r8
  HeatStorCanopy_pft(NZ) = 0.0_r8
  VapXAir2Canopy_pft(NZ) = 0.0_r8
  Transpiration_pft(NZ)  = 0.0_r8
  QdewCanopy_pft(NZ)     = 0.0_r8
  IF(CanopyHeight_pft(NZ).GE.SnowDepth-ZERO .or. TKSnow==spval)THEN
    TKC_pft(NZ)=TairK
  ELSE
    TKC_pft(NZ)=TKSnow
  ENDIF
  TdegCCanopy_pft(NZ) = units%Kelvin2Celcius(TKC_pft(NZ))
  FTHRM                  = EMMC*stefboltz_const*FracPARads2Canopy_pft(NZ)*AREA3(NU)
  LWRadCanopy_pft(NZ)    = FTHRM*TKC_pft(NZ)**4._r8
  PSICanopy_pft(NZ)      = TotalSoilPSIMPa_vr(NGTopRootLayer_pft(NZ))
  CCPOLT                 = sum(CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ))

  call update_osmo_turg_pressure(PSICanopy_pft(NZ),CCPOLT,CanOsmoPsi0pt_pft(NZ),TKC_pft(NZ)&
    ,PSICanopyOsmo_pft(NZ),PSICanopyTurg_pft(NZ),FDMP)

  Stomata_Stress             = EXP(RCS(NZ)*PSICanopyTurg_pft(NZ))
  CanPStomaResistH2O_pft(NZ) = MinCanPStomaResistH2O_pft(NZ) &
    +(H2OCuticleResist_pft(NZ)-MinCanPStomaResistH2O_pft(NZ))*Stomata_Stress
  CanopyBndlResist_pft(NZ) = ReistanceCanopy_pft(NZ)
  VHeatCapCanopy_pft(NZ)     = cpw*(ShootStrutElms_pft(ielmc,NZ)*10.0E-06_r8)
  DeltaTKC_pft(NZ)         = 0.0_r8

  DO N=1,MY(NZ)
    DO  L=NU,MaxSoiL4Root_pft(NZ)
      PSIRoot_pvr(N,L,NZ)              = TotalSoilPSIMPa_vr(L)
      CCPOLT                           = sum(RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ))
      AllPlantRootH2OLoss_vr(N,L,NZ) = 0.0_r8

      call update_osmo_turg_pressure(PSIRoot_pvr(N,L,NZ),CCPOLT,CanOsmoPsi0pt_pft(NZ),TKS_vr(L),&
        PSIRootOSMO_vr(N,L,NZ),PSIRootTurg_vr(N,L,NZ))

    enddo
  ENDDO
  end associate
  end subroutine HandleBareSoil
!------------------------------------------------------------------------

  subroutine UpdatePlantWaterVars(NZ,VHeatCapCanopyAir,TotalSoilPSIMPa_vr,RootResist,SoiH2OResist,SoilRootResistance_rvr,&
    TKCX,VHeatCapCanopyPrev_pft,PrecpHeatbyCanopy,cumPRootH2OUptake,CumPlantHeatLoss2Soil,HeatEvapSens,SoiLayerHasRoot_rvr)
  !
  !Description
  !Update canopy heat and water states after the numerical iterations.
  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: VHeatCapCanopyAir
  real(r8), intent(in) :: TotalSoilPSIMPa_vr(JZ1)    !Elevation adjusted soil water potential [MPa]
  real(r8), intent(in) :: RootResist(jroots,JZ1)
  real(r8), intent(in) :: SoiH2OResist(jroots,JZ1)
  real(r8), intent(in) :: SoilRootResistance_rvr(jroots,JZ1)
  real(r8), intent(in) :: TKCX,VHeatCapCanopyPrev_pft
  real(r8), intent(in) :: PrecpHeatbyCanopy
  real(r8), intent(in) :: cumPRootH2OUptake
  real(r8), intent(in) :: HeatEvapSens                     !sensible heat associated with evaporation [MJ/h]
  real(r8), intent(in) :: CumPlantHeatLoss2Soil
  logical , intent(in) :: SoiLayerHasRoot_rvr(jroots,JZ1)  !indicator of root prescence

  character(len=*), parameter :: subname='UpdatePlantWaterVars'
  real(r8) :: CCPOLT,CanopyMassC
  real(r8) :: OSWT
  integer :: N,L
  associate(                                                         &
    NU                        => plt_site%NU,                        &
    CanOsmoPsi0pt_pft         => plt_ew%CanOsmoPsi0pt_pft,           &
    TairK                     => plt_ew%TairK,                       &
    TKC_pft                   => plt_ew%TKC_pft,                     &
    TKS_vr                    => plt_ew%TKS_vr,                      &
    TKCanopy_pft              => plt_ew%TKCanopy_pft,                &
    WatHeldOnCanopy_pft       => plt_ew%WatHeldOnCanopy_pft,         &
    CanopyBiomWater_pft       => plt_ew%CanopyBiomWater_pft,         &
    VHeatCapCanopy_pft        => plt_ew%VHeatCapCanopy_pft,          &
    HeatStorCanopy_pft        => plt_ew%HeatStorCanopy_pft,          & !loss of heat from canopy
    PrecIntcptByCanopy_pft    => plt_ew%PrecIntcptByCanopy_pft,      &
    Transpiration_pft         => plt_ew%Transpiration_pft,           &
    CanopyLeafShethC_pft      => plt_biom%CanopyLeafShethC_pft,      &
    CanopyStalkC_pft          => plt_biom%CanopyStalkC_pft,          &
    VapXAir2Canopy_pft        => plt_ew%VapXAir2Canopy_pft,          &
    PSIRoot_pvr               => plt_ew%PSIRoot_pvr,                 &
    HeatXAir2PCan_pft         => plt_ew%HeatXAir2PCan_pft,           &
    PSIRootOSMO_vr            => plt_ew%PSIRootOSMO_vr,              &
    PSIRootTurg_vr            => plt_ew%PSIRootTurg_vr,              &
    PSICanopy_pft             => plt_ew%PSICanopy_pft,               &
    RootNonstructElmConc_rpvr => plt_biom%RootNonstructElmConc_rpvr, &
    MY                        => plt_morph%MY,                       &
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft          &
  )
  !
  !     CANOPY SURFACE WATER STORAGE, SENSIBLE AND STORAGE HEAT FLUXES
  !     (NOT EXPLICITLY CALCULATED IN CONVERGENCE SOLUTION)
  !
  !     VOLWP,WatHeldOnCanopy_pft=water volume in canopy,on canopy surfaces
  !     HeatXAir2PCan_pft=canopy sensible,storage heat fluxes
  !     HeatStorCanopy_pft=Radiation-latent heat-sens, i.e. (>0) residual heat from air to canopy
  !     VHeatCapCanopyPrev_pft,VHeatCapCanopy_pft=previous,current canopy heat capacity
  !     VHeatCapCanopyAir=canopy sensible heat conductance
  !     HeatEvapSens=convective heat flux from latent heat flux (>0)
  !     PrecpHeatbyCanopy=convective heat flux from precip to canopy
  !     cumPRootH2OUptake < 0., add water to canopy 
  !     Transpiration_pft < 0., lost from canopy
  CanopyBiomWater_pft(NZ)     = CanopyBiomWater_pft(NZ)+Transpiration_pft(NZ)-cumPRootH2OUptake
  CanopyMassC             = AZMAX1(CanopyLeafShethC_pft(NZ)+CanopyStalkC_pft(NZ))
  WatHeldOnCanopy_pft(NZ) = WatHeldOnCanopy_pft(NZ)+PrecIntcptByCanopy_pft(NZ)+VapXAir2Canopy_pft(NZ)
  VHeatCapCanopy_pft(NZ)  = cpw*(CanopyMassC*SpecStalkVolume+WatHeldOnCanopy_pft(NZ)+CanopyBiomWater_pft(NZ))
  HeatXAir2PCan_pft(NZ)   = VHeatCapCanopyAir*(TairK-TKCanopy_pft(NZ))
  HeatStorCanopy_pft(NZ)  = TKCX*VHeatCapCanopyPrev_pft-TKCanopy_pft(NZ)*VHeatCapCanopy_pft(NZ)+HeatEvapSens+PrecpHeatbyCanopy
  !
  !     ROOT TOTAL, OSMOTIC AND TURGOR WATER POTENTIALS
  !
  !     PSIRoot_pvr,PSICanopy_pft=root,canopy total water potential
  !     TotalSoilPSIMPa_vr=total soil water potential PSIST adjusted for surf elevn
  !     SoiH2OResist,SoilRootResistance_rvr,RootResist=soil,soil+root,root radial+axial resistance
  !     PSIRootOSMO_vr,PSIRootTurg_vr=root osmotic,turgor water potential
  !     CanOsmoPsi0pt_pft=osmotic potential at PSIRoot_pvr=0 from PFT file

  !compute root pressure assuming zero water storage capacity in root, cf. Eq. (36) Grant (1998), Ecological modelling.
  D4505: DO N=1,MY(NZ)
    D4510: DO L=NU,MaxSoiL4Root_pft(NZ)
      IF(SoiLayerHasRoot_rvr(N,L))THEN
        PSIRoot_pvr(N,L,NZ)=AZMIN1((TotalSoilPSIMPa_vr(L)*RootResist(N,L) &
          +PSICanopy_pft(NZ)*SoiH2OResist(N,L))/SoilRootResistance_rvr(N,L))
      ELSE
        PSIRoot_pvr(N,L,NZ)=TotalSoilPSIMPa_vr(L)
      ENDIF           
      !obtain total reserve
      CCPOLT=sum(RootNonstructElmConc_rpvr(1:NumPlantChemElms,N,L,NZ))

      CALL update_osmo_turg_pressure(PSIRoot_pvr(N,L,NZ),CCPOLT,CanOsmoPsi0pt_pft(NZ),TKS_vr(L),&
        PSIRootOSMO_vr(N,L,NZ),PSIRootTurg_vr(N,L,NZ))

    ENDDO D4510
  ENDDO D4505
  call PrintInfo('end '//subname)
  end associate
  end subroutine UpdatePlantWaterVars
!------------------------------------------------------------------------

  subroutine SetCanopyGrowthFuncs(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ

  character(len=*), parameter :: subname='SetCanopyGrowthFuncs'
  real(r8) :: ACTV,RTK,STK,TKGO,TKSO
  integer :: L
  associate(                                              &
    TdegCCanopy_pft  => plt_ew%TdegCCanopy_pft,     &
    TairK               => plt_ew%TairK,                  &
    TKC_pft             => plt_ew%TKC_pft,                &
    TKS_vr              => plt_ew%TKS_vr,                 &
    PSICanPDailyMin     => plt_ew%PSICanPDailyMin,        &
    PSICanopy_pft       => plt_ew%PSICanopy_pft,          &
    NU                  => plt_site%NU,                   &
    ChillHours_pft      => plt_photo%ChillHours_pft,      &
    TempOffset_pft      => plt_pheno%TempOffset_pft,      &
    TCChill4Seed_pft    => plt_pheno%TCChill4Seed_pft,    &
    fTgrowRootP_vr      => plt_pheno%fTgrowRootP_vr,      &
    TCGroth_pft         => plt_pheno%TCGroth_pft,         &
    TKGroth_pft         => plt_pheno%TKGroth_pft,         &
    iPlantCalendar_brch => plt_pheno%iPlantCalendar_brch, &
    fTCanopyGroth_pft   => plt_pheno%fTCanopyGroth_pft,   &
    MaxSoiL4Root_pft    => plt_morph%MaxSoiL4Root_pft,    &
    MainBranchNum_pft   => plt_morph%MainBranchNum_pft    &
  )
  !
  !     SET CANOPY GROWTH TEMPERATURE FROM SOIL SURFACE
  !     OR CANOPY TEMPERATURE DEPENDING ON GROWTH STAGE
  !
  call PrintInfo('beg '//subname)
  IF(iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).EQ.0)THEN
    !before seed emergence
    TKGroth_pft(NZ)=TKS_vr(NU)
    !     ELSEIF((iPlantTurnoverPattern_pft(NZ).EQ.0.OR.iPlantRootProfile_pft(NZ).LE.1)
    !    2.AND.iPlantCalendar_brch(ipltcal_InitFloral,MainBranchNum_pft(NZ),NZ).EQ.0)THEN
    !     TKGroth_pft(NZ)=TKS_vr(NU)
  ELSE
    TKGroth_pft(NZ)=TKC_pft(NZ)
  ENDIF
  TCGroth_pft(NZ)=units%Kelvin2Celcius(TKGroth_pft(NZ))
  !
  !     ARRHENIUS FUNCTION FOR CANOPY AND ROOT GROWTH WITH OFFSET
  !     FOR ZONE OF THERMAL ADAPTATION ENTERED IN 'READQ'
  !
  !     TKGroth_pft,TKGO=canopy temperature,canopy temp used in Arrhenius eqn
  !     TKS_vr,TKSO=soil temperature,soil temp used in Arrhenius eqn
  !     TempOffset_pft=shift in Arrhenius curve for thermal adaptation
  !     fTCanopyGroth_pft,fTgrowRootP_vr=temperature function for canopy,root growth (25 oC =1)
  !     RGAS,710.0=gas constant,enthalpy
  !     62500,197500,222500=energy of activn,high,low temp inactivn(KJ mol-1)
  !     PSICanPDailyMin=minimum daily canopy water potential
  !
  TKGO                  = TKGroth_pft(NZ)+TempOffset_pft(NZ)
  fTCanopyGroth_pft(NZ) = calc_canopy_grow_tempf(TKGO)

  D100: DO L=NU,MaxSoiL4Root_pft(NZ)
    TKSO                 = TKS_vr(L)+TempOffset_pft(NZ)
    fTgrowRootP_vr(L,NZ) = calc_root_grow_tempf(TKSO)
  ENDDO D100
  PSICanPDailyMin(NZ)=AMIN1(PSICanPDailyMin(NZ),PSICanopy_pft(NZ))
  !
  !     DIURNAL CHILLING
  !
  !     TCChill4Seed_pft=chilling temperature from PFT file
  !     CHILL=accumulated chilling hours used to limit CO2 fixn in stomate.f
  !
  IF(TdegCCanopy_pft(NZ).LT.TCChill4Seed_pft(NZ))THEN
    ChillHours_pft(NZ)=AMIN1(24.0_r8,ChillHours_pft(NZ)+1.0_r8)
  ELSE
    ChillHours_pft(NZ)=AZMAX1(ChillHours_pft(NZ)-1.0_r8)
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine SetCanopyGrowthFuncs

end module UptakesMod
