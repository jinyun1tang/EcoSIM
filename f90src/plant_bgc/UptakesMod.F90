module UptakesMod
  use data_kind_mod , only : r8 => DAT_KIND_R8
  use data_const_mod, only : GravAcceleration=>DAT_CONST_G
  use StomatesMod   , only : stomates
  use minimathmod
  use UnitMod       , only : units
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use NutUptakeMod
  use PlantAPIData
  use PlantMathFuncMod

  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public :: uptakes
  public :: InitUptake

  real(r8), parameter :: mGravAccelerat=1.e-3_r8*GravAcceleration  !gravitational constant devided by 1000.
  contains

  subroutine InitUptake

  implicit none

  call InitUptakePars

  end subroutine InitUptake
!------------------------------------------------------------------------------------------

  subroutine uptakes(I,J)
!
!     THIS subroutine CALCULATES EXCHANGES OF ENERGY, C, N AND P
!     BETWEEN THE CANOPY AND THE ATMOSPHERE AND BETWEEN ROOTS AND THE SOIL
!
  implicit none
  integer, intent(in) :: I, J

  integer :: NN,N,NZ,K,L
  real(r8) :: FracGrndByPFT
  real(r8) :: PopPlantO2Uptake,PopPlantO2Demand,HeatSensConductCanP
  real(r8) :: CanPMassC
  real(r8) :: ElvAdjstedtSoiPSIMPa(JZ1)
  real(r8) :: PATH(2,JZ1)
  real(r8) :: FineRootRadius(2,JZ1),RootResist(2,JZ1)
  real(r8) :: RootResistSoi(2,JZ1),RootResistRadial(2,JZ1)
  real(r8) :: RootResistAxial(2,JZ1),SoiH2OResist(2,JZ1)
  real(r8) :: SoiAddRootResist(2,JZ1),AllRootC_vr(JZ1)
  real(r8) :: FracPRoot4Uptake(2,JZ1,JP1),MinFracPRoot4Uptake(2,JZ1,JP1)
  real(r8) :: FracSoiLayByPrimRoot(JZ1,JP1),RootAreaDivRadius(2,JZ1)
  real(r8) :: AirPoreAvail4Fill(JZ1),WatAvail4Uptake(JZ1)
  real(r8) :: TKCX,CNDT,PrecHeatIntcptByCanP1,VHCPX
  real(r8) :: DIFF,PSILH,cumPRootH2OUptake,VFLXC
  real(r8) :: FDMP
  integer :: LayrHasRoot(2,JZ1)
!     begin_execution
  associate(                         &
    CanPbndlResist    => plt_photo%CanPbndlResist    , &
    CanWatP  => plt_ew%CanWatP     , &
    TCelciusCanopy   => plt_ew%TCelciusCanopy      , &
    TKCZ   => plt_ew%TKCZ      , &
    PrecIntcptByCanP   => plt_ew%PrecIntcptByCanP      , &
    EvapTransHeatP  => plt_ew%EvapTransHeatP     , &
    TKC    => plt_ew%TKC       , &
    PSICanP  => plt_ew%PSICanP     , &
    HeatXAir2PCan  => plt_ew%HeatXAir2PCan     , &
    Canopy_Heat_Sens_col   => plt_ew%Canopy_Heat_Sens_col      , &
    DTKC   => plt_ew%DTKC      , &
    VapXAir2PCan  => plt_ew%VapXAir2PCan     , &
    TairK    => plt_ew%TairK       , &
    Canopy_Heat_Latent_col   => plt_ew%Canopy_Heat_Latent_col      , &
    PTrans    => plt_ew%PTrans       , &
    WatByPCan  => plt_ew%WatByPCan     , &
    CumSoilThickness => plt_site%CumSoilThickness  , &
    PPT    => plt_site%PPT     , &
    NP     => plt_site%NP      , &
    pftPlantPopulation     => plt_site%pftPlantPopulation      , &
    NU     => plt_site%NU      , &
    AREA3  => plt_site%AREA3   , &
    ZEROS  => plt_site%ZEROS   , &
    PopPlantRootC_vr  => plt_biom%PopPlantRootC_vr   , &
    ZEROL  => plt_biom%ZEROL   , &
    ZEROP  => plt_biom%ZEROP   , &
    CanopyLeafShethC_pft   => plt_biom%CanopyLeafShethC_pft    , &
    CanPStalkC  => plt_biom%CanPStalkC   , &
    iPlantCalendar  => plt_pheno%iPlantCalendar  , &
    PlantO2Stress   => plt_pheno%PlantO2Stress   , &
    IsPlantActive  => plt_pheno%IsPlantActive  , &
    CanopyLA_grd  => plt_morph%CanopyLA_grd  , &
    PrimRootDepth  => plt_morph%PrimRootDepth  , &
    CanopyArea_grid  => plt_morph%CanopyArea_grid  , &
    CanopyArea_pft  => plt_morph%CanopyArea_pft  , &
    SeedinDepth  => plt_morph%SeedinDepth  , &
    NumOfMainBranch_pft    => plt_morph%NumOfMainBranch_pft    , &
    FracPARByCanP  => plt_rad%FracPARByCanP      &
  )

  call PrepH2ONutrientUptake(ElvAdjstedtSoiPSIMPa,AllRootC_vr,AirPoreAvail4Fill,WatAvail4Uptake)
!

!     IF PLANT SPECIES EXISTS
!
  DO NZ=1,NP
    PopPlantO2Uptake=0.0_r8
    PopPlantO2Demand=0.0_r8

    IF(IsPlantActive(NZ).EQ.iPlantIsActive.AND.pftPlantPopulation(NZ).GT.0.0_r8)THEN

      call UpdateCanopyProperty(NZ)

!     STOMATE=solve for minimum canopy stomatal resistance
      CALL STOMATEs(I,J,NZ)
!
!     CALCULATE VARIABLES USED IN ROOT UPTAKE OF WATER AND NUTRIENTS
      call UpdateRootProperty(NZ,PATH,FineRootRadius,AllRootC_vr,FracPRoot4Uptake,MinFracPRoot4Uptake,FracSoiLayByPrimRoot,RootAreaDivRadius)
!
!     CALCULATE CANOPY WATER STATUS FROM CONVERGENCE SOLUTION FOR
!     TRANSPIRATION - ROOT WATER UPTAKE = CHANGE IN CANOPY WATER CONTENT
!
!     (AG: - originally this line had a N0B1 here )
      IF((iPlantCalendar(ipltcal_Emerge,NumOfMainBranch_pft(NZ),NZ).NE.0).AND.(CanopyArea_pft(NZ).GT.ZEROL(NZ) &
        .AND.FracPARByCanP(NZ).GT.0.0_r8).AND.(PrimRootDepth(1,1,NZ).GT.SeedinDepth(NZ)+CumSoilThickness(0)))THEN
        !leaf area > 0, absorped par>0, and rooting depth > seeding depth
!
        call CalcResistance(NZ,PATH,FineRootRadius,RootAreaDivRadius,RootResist,RootResistSoi,&
          RootResistRadial,RootResistAxial,SoiH2OResist,SoiAddRootResist,CNDT,PSILH,LayrHasRoot)
!
!     INITIALIZE CANOPY WATER POTENTIAL, OTHER VARIABLES USED IN ENERGY
!     BALANCE THAT DON'T NEED TO BE RECALCULATED DURING CONVERGENCE
!
!     PSICanP=initial estimate of total canopy water potential
!     PrecIntcptByCanP,PrecHeatIntcptByCanP1=convective water,heat flux from precip to canopy
!     FTHRM,LWRad2CanP=LW emitted,absorbed by canopy
!     FracGrndByPFT=fraction of ground area AREA occupied by PFT
!     CCPOLT=total nonstructural canopy C,N,P concentration
!     CCPOLP,CZPOLP,CPPOLP=nonstructural C,N,P concn in canopy
!     OSWT=molar mass of CCPOLT
!     TKCX=intermediate estimate of TKC used in convergence solution
!     CanPMassC=leaf+petiole+stalk mass
!     VSTK=specific stalk volume
!     VHCPX=canopy heat capacity
!     CanWatP,WatByPCan=water volume in canopy,on canopy surfaces
!
        PSICanP(NZ)=AMIN1(-ppmc,0.667_r8*PSICanP(NZ))
        PTrans(NZ)=0.0_r8
        VapXAir2PCan(NZ)=0.0_r8
        PrecHeatIntcptByCanP1=PrecIntcptByCanP(NZ)*cpw*TairK

        IF(CanopyArea_grid.GT.ZEROS)THEN
          !the grid has significant canopy (leaf+steam) area
          FracGrndByPFT=CanopyArea_pft(NZ)/CanopyArea_grid*AMIN1(1.0_r8,0.5_r8*CanopyLA_grd/AREA3(NU))
        ELSEIF(PPT.GT.ZEROS)THEN
          !total population is > 0
          FracGrndByPFT=pftPlantPopulation(NZ)/PPT
        ELSE
          FracGrndByPFT=1.0_r8/NP
        ENDIF

        TKCX=TKC(NZ)
        !Canopy C mass, g/m2
        CanPMassC=AZMAX1(CanopyLeafShethC_pft(NZ)+CanPStalkC(NZ))
        !canopy heat capacity J/K
        VHCPX=cpw*(CanPMassC*VSTK+WatByPCan(NZ)+CanWatP(NZ))
!
!     CONVERGENCE SOLUTION
!
        NN=CanopyEnergyH2OIteration(I,J,NZ,FracGrndByPFT,CanPMassC,&
          ElvAdjstedtSoiPSIMPa,HeatSensConductCanP,DIFF,cumPRootH2OUptake,&
          VFLXC,FDMP,SoiAddRootResist,FracPRoot4Uptake,AirPoreAvail4Fill,&
          WatAvail4Uptake,TKCX,CNDT,VHCPX,PrecHeatIntcptByCanP1,PSILH,LayrHasRoot)
!
!     FINAL CANOPY TEMPERATURE, DIFFERENCE WITH AIR TEMPERATURE
!
!     TKC=final estimate of canopy temperature TKCZ
!     TairK=current air temperature
!     DTKC=TKC-TairK for next hour
!
        TKC(NZ)=TKCZ(NZ)
        TCelciusCanopy(NZ)=units%Kelvin2Celcius(TKC(NZ))
        DTKC(NZ)=TKC(NZ)-TairK
!
!     IF CONVERGENCE NOT ACHIEVED (RARE), SET DEFAULT
!     TEMPERATURES, ENERGY FLUXES, WATER POTENTIALS, RESISTANCES
!
        call HandlingDivergence(I,J,NN,NZ,ElvAdjstedtSoiPSIMPa,DIFF,FDMP)

        call UpdateCanopyWater(NZ,HeatSensConductCanP,ElvAdjstedtSoiPSIMPa,RootResist,SoiH2OResist,&
          SoiAddRootResist,TKCX,VHCPX,PrecHeatIntcptByCanP1,cumPRootH2OUptake,VFLXC,LayrHasRoot)
!
!     DEFAULT VALUES IF PLANT SPECIES DOES NOT EXIST
!
      ELSE
        call HandleBareSoil(NZ,ElvAdjstedtSoiPSIMPa,FDMP)
      ENDIF

      call SetCanopyGrowthFuncs(NZ)

      call PopPlantNutientO2Uptake(NZ,FDMP,PopPlantO2Uptake,PopPlantO2Demand,&
        PATH,FineRootRadius,FracPRoot4Uptake,MinFracPRoot4Uptake,FracSoiLayByPrimRoot,RootAreaDivRadius)

      Canopy_Heat_Latent_col=Canopy_Heat_Latent_col+EvapTransHeatP(NZ)*CanPbndlResist(NZ)
      Canopy_Heat_Sens_col=Canopy_Heat_Sens_col+HeatXAir2PCan(NZ)*CanPbndlResist(NZ)
      IF(PopPlantO2Demand.GT.ZEROP(NZ))THEN
        PlantO2Stress(NZ)=PopPlantO2Uptake/PopPlantO2Demand
      ELSE
        PlantO2Stress(NZ)=0.0_r8
      ENDIF
    ENDIF
  ENDDO
  RETURN
  end associate
  END subroutine uptakes
!------------------------------------------------------------------------

  subroutine PrepH2ONutrientUptake(ElvAdjstedtSoiPSIMPa,AllRootC_vr,AirPoreAvail4Fill,WatAvail4Uptake)
!
!     prepare for uptake calculation
  implicit none
  real(r8), intent(out) :: ElvAdjstedtSoiPSIMPa(JZ1)
  real(r8), intent(out) :: AllRootC_vr(JZ1),AirPoreAvail4Fill(JZ1),WatAvail4Uptake(JZ1)
  integer :: NZ, L, N
  real(r8) :: ARLSC

  associate(                          &
    ZERO   => plt_site%ZERO     , &
    NP     => plt_site%NP       , &
    NJ     => plt_site%NJ       , &
    VLWatMicPM  => plt_site%VLWatMicPM    , &
    ALT    => plt_site%ALT      , &
    NU     => plt_site%NU       , &
    NP0    => plt_site%NP0      , &
    TotalSoilH2OPSIMPa  => plt_ew%TotalSoilH2OPSIMPa      , &
    PopPlantRootC_vr  => plt_biom%PopPlantRootC_vr    , &
    VLMicP   => plt_soilchem%VLMicP   , &
    VLiceMicP   => plt_soilchem%VLiceMicP , &
    THETY  => plt_soilchem%THETY, &
    SoiBulkDensity   => plt_soilchem%SoiBulkDensity , &
    VLWatMicP   => plt_soilchem%VLWatMicP , &
    VLSoilMicP   => plt_soilchem%VLSoilMicP , &
    CanopyStemA_pft  => plt_morph%CanopyStemA_pft   , &
    MY     => plt_morph%MY      , &
    CanopyLeafA_pft  => plt_morph%CanopyLeafA_pft   , &
    RadNet2CanP   => plt_rad%RadNet2CanP      , &
    LWRadCanP  => plt_rad%LWRadCanP       &
  )
!
!     RESET TOTAL UPTAKE ARRAYS
!
!     CanopyLeafA_pft,CanopyStemA_pft=leaf,stalk areas
!

  ARLSC=0.0_r8
  D9984: DO NZ=1,NP0
!     TKC(NZ)=TairK+DTKC(NZ)
!     TCelciusCanopy(NZ)=TKC(NZ)-TC2K
    ARLSC=ARLSC+CanopyLeafA_pft(NZ)+CanopyStemA_pft(NZ)
    RadNet2CanP(NZ)=0.0_r8
    plt_ew%EvapTransHeatP(NZ)=0.0_r8
    plt_ew%HeatXAir2PCan(NZ)=0.0_r8
    plt_ew%HeatStorCanP(NZ)=0.0_r8
    LWRadCanP(NZ)=0.0_r8
    plt_ew%PTrans(NZ)=0.0_r8
    plt_ew%VapXAir2PCan(NZ)=0.0_r8
    plt_rbgc%UPOME(1:NumOfPlantChemElements,NZ)=0.0_r8
    plt_rbgc%UPNH4(NZ)=0.0_r8
    plt_rbgc%UPNO3(NZ)=0.0_r8
    plt_rbgc%UPH2P(NZ)=0.0_r8
    plt_rbgc%UPH1P(NZ)=0.0_r8
    plt_rbgc%UPNF(NZ)=0.0_r8
!
!     RESET UPTAKE ARRAYS
!
    DO  L=NU,NJ
      DO  N=1,MY(NZ)
        plt_ew%AllPlantRootH2OUptake_vr(N,L,NZ)=0.0_r8
        plt_rbgc%RCO2P(N,L,NZ)=0.0_r8
        plt_rbgc%RUPOXP(N,L,NZ)=0.0_r8
        plt_rbgc%RCO2S(N,L,NZ)=0.0_r8
        plt_rbgc%RUPOXS(N,L,NZ)=0.0_r8
        plt_rbgc%RUPCHS(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN2S(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN3S(N,L,NZ)=0.0_r8
        plt_rbgc%RUPN3B(N,L,NZ)=0.0_r8
        plt_rbgc%RUPHGS(N,L,NZ)=0.0_r8
        plt_rbgc%trcg_air2root_flx_pft_vr(idg_beg:idg_end-1,N,L,NZ)=0.0_r8
        plt_rbgc%trcg_Root_DisEvap_flx_vr(idg_beg:idg_end-1,N,L,NZ)=0.0_r8
      enddo
    enddo
  ENDDO D9984
!
!     ElvAdjstedtSoiPSIMPa=total soil water potential PSIST adjusted for surf elevn
!     ALT=surface elevation
!     WatAvail4Uptake,VLWatMicPM=water volume available for uptake,total water volume
!     THETY,VLSoilPoreMicP=hygroscopic SWC,soil volume
!     AirPoreAvail4Fill=air volume, used to fill by water coming out from plant root/mycorrhizae
!     AllRootC_vr=total biome root mass
!
  D9000: DO L=NU,NJ
    ElvAdjstedtSoiPSIMPa(L)=TotalSoilH2OPSIMPa(L)-mGravAccelerat*ALT
    IF(SoiBulkDensity(L).GT.ZERO)THEN
      WatAvail4Uptake(L)=VLWatMicPM(NPH,L)-THETY(L)*VLSoilMicP(L)     !maximum amount of water for uptake
      AirPoreAvail4Fill(L)=AZMAX1(VLMicP(L)-VLWatMicP(L)-VLiceMicP(L))  !air volume
    ELSE
      WatAvail4Uptake(L)=VLWatMicPM(NPH,L)
      AirPoreAvail4Fill(L)=0.0_r8
    ENDIF
    AllRootC_vr(L)=0.0_r8
    D9005: DO NZ=1,NP
      DO  N=1,MY(NZ)
!     IF(IsPlantActive(NZ).EQ.iPlantIsActive.AND.pftPlantPopulation(NZ).GT.0.0)THEN
      AllRootC_vr(L)=AllRootC_vr(L)+AZMAX1(PopPlantRootC_vr(N,L,NZ))
!     ENDIF
      enddo
    ENDDO D9005
  ENDDO D9000
  end associate
  end subroutine PrepH2ONutrientUptake
!------------------------------------------------------------------------
  subroutine UpdateCanopyProperty(NZ)
!
!     update canopy characterization
  implicit none
  integer, intent(in) :: NZ
  real(r8) :: ALFZ
  real(r8) :: TFRADP,RACZ(JP1)
  integer :: NB,K,L,N,NZZ

  associate(                          &
    ZERO   =>  plt_site%ZERO    , &
    NP     =>  plt_site%NP      , &
    IETYP  =>  plt_site%IETYP   , &
    BndlResistAboveCanG   =>  plt_ew%BndlResistAboveCanG      , &
    TairK    =>  plt_ew%TairK       , &
    DTKC   =>  plt_ew%DTKC      , &
    ZeroPlanDisp     =>  plt_ew%ZeroPlanDisp        , &
    RoughHeight     =>  plt_ew%RoughHeight        , &
    RAZ    =>  plt_ew%RAZ       , &
    TKCZ   =>  plt_ew%TKCZ      , &
    WSLF   =>  plt_biom%WSLF    , &
    ZEROP  =>  plt_biom%ZEROP   , &
    FracPARByCanP  =>  plt_rad%FracPARByCanP    , &
    LeafAUnshaded_seclyrnodbrpft  =>  plt_photo%LeafAUnshaded_seclyrnodbrpft  , &
    ARLF1  =>  plt_morph%ARLF1  , &
    KLEAFX =>  plt_morph%KLEAFX , &
    CanopyHeight     =>  plt_morph%CanopyHeight     , &
    ClumpFactort    =>  plt_morph%ClumpFactort    , &
    CanopyArea_pft  =>  plt_morph%CanopyArea_pft  , &
    NumOfBranches_pft    =>  plt_morph%NumOfBranches_pft    , &
    LeafA_lyrnodbrchpft   =>  plt_morph%LeafA_lyrnodbrchpft   , &
    GridMaxCanopyHeight     =>  plt_morph%GridMaxCanopyHeight       &
  )
!
!     APPLY CLUMPING FACTOR TO LEAF SURFACE AREA DEFINED BY
!     INCLINATION N, LAYER L, NODE K, BRANCH NB, SPECIES NZ,
!     N-S POSITION NY, E-W POSITION NX(AZIMUTH M ASSUMED UNIFORM)
!
  D500: DO NB=1,NumOfBranches_pft(NZ)
    D550: DO K=1,MaxNodesPerBranch1
!
!     NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
!
!     ARLF=leaf area
!     WSLF=leaf protein content
!     LeafAUnshaded_seclyrnodbrpft,LeafA_lyrnodbrchpft=unself-shaded,total leaf surface area
!     ClumpFactort=clumping factor from PFT file
!
      IF(ARLF1(K,NB,NZ).GT.ZEROP(NZ).AND.WSLF(K,NB,NZ).GT.ZEROP(NZ))THEN
        KLEAFX(NB,NZ)=K
      ENDIF
      D600: DO L=NumOfCanopyLayers1,1,-1
        D650: DO N=1,JLI1
          LeafAUnshaded_seclyrnodbrpft(N,L,K,NB,NZ)=LeafA_lyrnodbrchpft(N,L,K,NB,NZ)*ClumpFactort(NZ)
        ENDDO D650
      ENDDO D600
    ENDDO D550
  ENDDO D500
!
!     CANOPY HEIGHT FROM HEIGHT OF MAXIMUM LEAF LAYER
!
!     DIFFUSIVE RESISTANCE OF OTHER TALLER CANOPIES TO HEAT AND VAPOR
!     TRANSFER OF CURRENT CANOPY ADDED TO BOUNDARY LAYER RESISTANCE
!     OF TALLEST CANOPY CALCULATED IN 'HOUR1'
!
!     IETYP=Koppen climate zone
!     ZC,ZT,RoughHeight=PFT canopy,grid biome,surface roughness height
!     FracPARByCanP=fraction of radiation received by each PFT canopy
!     ALFZ=shape parameter for windspeed attenuation in canopy
!     BndlResistAboveCanG=biome canopy isothermal boundary layer resistance
!     RACZ,RAZ=additional,total PFT canopy isothermal blr
!
  IF(CanopyArea_pft(NZ).GT.0.0_r8)THEN
    IF(IETYP.GE.0)THEN
      TFRADP=0.0_r8
      D700: DO NZZ=1,NP
        IF(CanopyHeight(NZZ).GT.CanopyHeight(NZ)+RoughHeight)THEN
          TFRADP=TFRADP+FracPARByCanP(NZZ)
        ENDIF
      ENDDO D700
      ALFZ=2.0_r8*TFRADP
      IF(BndlResistAboveCanG.GT.ZERO.AND.GridMaxCanopyHeight.GT.ZERO.AND.ALFZ.GT.ZERO)THEN
        RACZ(NZ)=AMIN1(RACX,AZMAX1(GridMaxCanopyHeight*EXP(ALFZ) &
          /(ALFZ/BndlResistAboveCanG)*(EXP(-ALFZ*CanopyHeight(NZ)/GridMaxCanopyHeight) &
          -EXP(-ALFZ*(ZeroPlanDisp+RoughHeight)/GridMaxCanopyHeight))))
      ELSE
        RACZ(NZ)=0.0_r8
      ENDIF
    ELSE
      RACZ(NZ)=0.0_r8
    ENDIF
  ELSE
    RACZ(NZ)=RACX
  ENDIF
  RAZ(NZ)=BndlResistAboveCanG+RACZ(NZ)
!
!     INITIALIZE CANOPY TEMPERATURE WITH CURRENT AIR TEMPERATURE AND
!     LAST HOUR'S CANOPY-AIR TEMPERATURE DIFFERENCE, AND CALL A
!     subroutine TO CALCULATE MINIMUM CANOPY STOMATAL RESISTANCE
!     FOR SUBSEQUENT USE IN ENERGY EXCHANGE CALCULATIONS
!
!     TKCZ=initial estimate of canopy temperature TKC
!     TairK=current air temperature
!     DTKC=TKC-TairK from previous hour
!
  TKCZ(NZ)=TairK+DTKC(NZ)
  end associate
  end subroutine UpdateCanopyProperty
!------------------------------------------------------------------------
  subroutine UpdateRootProperty(NZ,PATH,FineRootRadius,AllRootC_vr,&
    FracPRoot4Uptake,MinFracPRoot4Uptake,FracSoiLayByPrimRoot,RootAreaDivRadius)
!
!     update root characterization

  implicit none
  integer , intent(in) :: NZ
  real(r8), intent(in) :: AllRootC_vr(JZ1)
  real(r8), intent(out) :: PATH(2,JZ1),FineRootRadius(2,JZ1)
  real(r8), intent(out) :: FracPRoot4Uptake(2,JZ1,JP1),MinFracPRoot4Uptake(2,JZ1,JP1)
  real(r8), intent(out) :: FracSoiLayByPrimRoot(JZ1,JP1),RootAreaDivRadius(2,JZ1)
  real(r8) :: RootDepZ,RTDPX
  integer :: N,L,NR

  associate(                           &
    PopPlantRootC_vr  =>  plt_biom%PopPlantRootC_vr    , &
    FracSoiAsMicP   =>  plt_site%FracSoiAsMicP     , &
    CumSoilThickness =>  plt_site%CumSoilThickness   , &
    DLYR3  =>  plt_site%DLYR3    , &
    ZEROS  =>  plt_site%ZEROS    , &
    ZERO   =>  plt_site%ZERO     , &
    pftPlantPopulation     =>  plt_site%pftPlantPopulation       , &
    NU     =>  plt_site%NU       , &
    NRT    =>  plt_morph%NRT     , &
    MY     =>  plt_morph%MY      , &
    RootLenDensNLP  =>  plt_morph%RootLenDensNLP   , &
    MaxSecndRootRadius1 =>  plt_morph%MaxSecndRootRadius1  , &
    HypoctoylHeight  =>  plt_morph%HypoctoylHeight   , &
    RootLenPerP  =>  plt_morph%RootLenPerP   , &
    PrimRootDepth  =>  plt_morph%PrimRootDepth   , &
    RootPorosity   =>  plt_morph%RootPorosity    , &
    RTVLW  =>  plt_morph%RTVLW   , &
    MaxSecndRootRadius =>  plt_morph%MaxSecndRootRadius  , &
    SeedinDepth  =>  plt_morph%SeedinDepth   , &
    NI     =>  plt_morph%NI        &
  )
!     RootDepZ,PrimRootDepth=primary root depth
!     FracSoiLayByPrimRoot=fraction of each soil layer with primary root
!     DLYR=layer thickness
!     CumSoilThickness=depth from soil surface to layer bottom
!     SeedinDepth=seeding depth
!     HypoctoylHeight=hypocotyledon height
!     FracPRoot4Uptake=PFT fraction of biome root mass
!     RootLenDensNLP,RootLenPerP=root length density,root length per plant
!     FineRootRadius=root radius
!     RootPorosity=root porosity
!     FracSoiAsMicP=micropore fraction excluding macropore,rock
!     RTVLW=root aqueous volume
!     PP=plant population
!     PATH=path length of water and nutrient uptake
!     RootAreaDivRadius=root surface area/radius for uptake,diffusivity
!
  D2000: DO N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      IF(N.EQ.ipltroot)THEN
        RootDepZ=0.0_r8
        D2005: DO NR=1,NRT(NZ)
          RootDepZ=AMAX1(RootDepZ,PrimRootDepth(ipltroot,NR,NZ))
        ENDDO D2005
        IF(L.EQ.NU)THEN
          FracSoiLayByPrimRoot(L,NZ)=1.0_r8
        ELSE
          IF(DLYR3(L).GT.ZERO)THEN
            RTDPX=AZMAX1(RootDepZ-CumSoilThickness(L-1))
            RTDPX=AZMAX1(AMIN1(DLYR3(L),RTDPX)-AZMAX1(SeedinDepth(NZ)-CumSoilThickness(L-1)-HypoctoylHeight(NZ)))
            FracSoiLayByPrimRoot(L,NZ)=RTDPX/DLYR3(L)
          ELSE
            FracSoiLayByPrimRoot(L,NZ)=0.0_r8
          ENDIF
        ENDIF

      ENDIF
      IF(AllRootC_vr(L).GT.ZEROS)THEN
        FracPRoot4Uptake(N,L,NZ)=AZMAX1(PopPlantRootC_vr(N,L,NZ))/AllRootC_vr(L)
      ELSE
        FracPRoot4Uptake(N,L,NZ)=1.0_r8
      ENDIF
      MinFracPRoot4Uptake(N,L,NZ)=FMN*FracPRoot4Uptake(N,L,NZ)
      IF(RootLenDensNLP(N,L,NZ).GT.ZERO.AND.FracSoiLayByPrimRoot(L,NZ).GT.ZERO)THEN
        FineRootRadius(N,L)=AMAX1(MaxSecndRootRadius1(N,NZ),SQRT((RTVLW(N,L,NZ) &
          /(1.0_r8-RootPorosity(N,NZ)))/(PICON*pftPlantPopulation(NZ)*RootLenPerP(N,L,NZ))))
        PATH(N,L)=AMAX1(1.001_r8*FineRootRadius(N,L),1.0_r8/(SQRT(PICON*(RootLenDensNLP(N,L,NZ)/FracSoiLayByPrimRoot(L,NZ))/FracSoiAsMicP(L))))
        RootAreaDivRadius(N,L)=TwoPiCON*RootLenPerP(N,L,NZ)/FracSoiLayByPrimRoot(L,NZ)
      ELSE
        FineRootRadius(N,L)=MaxSecndRootRadius(N,NZ)
        PATH(N,L)=1.001_r8*FineRootRadius(N,L)
        RootAreaDivRadius(N,L)=TwoPiCON*RootLenPerP(N,L,NZ)
      ENDIF
    enddo
  ENDDO D2000
  end associate
  end subroutine UpdateRootProperty
!------------------------------------------------------------------------

  subroutine HandlingDivergence(I,J,NN,NZ,ElvAdjstedtSoiPSIMPa,DIFF,FDMP)

  implicit none
  integer  , intent(in) :: NN, I, J
  integer  , intent(in) :: NZ
  REAL(R8) ,INTENT(IN) :: ElvAdjstedtSoiPSIMPa(JZ1),DIFF
  real(r8), intent(out):: FDMP
  real(r8) :: APSILT
  real(r8) :: CCPOLT
  real(r8) :: FTHRM,FDMR
  real(r8) :: OSWT,Stomata_Activity

  integer :: N,L
! begin_execution
  associate(                         &
   RAZ     => plt_ew%RAZ       , &
   DTKC    => plt_ew%DTKC      , &
   OSMO    => plt_ew%OSMO      , &
   PSICanPOsmo   => plt_ew%PSICanPOsmo     , &
   TKC     => plt_ew%TKC       , &
   TairK     => plt_ew%TairK       , &
   TKS     => plt_ew%TKS       , &
   PSIRootOSMO   => plt_ew%PSIRootOSMO     , &
   TCelciusCanopy    => plt_ew%TCelciusCanopy      , &
   AllPlantRootH2OUptake_vr   => plt_ew%AllPlantRootH2OUptake_vr     , &
   PSIRootTurg   => plt_ew%PSIRootTurg     , &
   PSICanPTurg   => plt_ew%PSICanPTurg     , &
   PSIRoot   => plt_ew%PSIRoot     , &
   VHeatCapCanP   => plt_ew%VHeatCapCanP     , &
   PSICanP   => plt_ew%PSICanP     , &
   NU      => plt_site%NU      , &
   AREA3   => plt_site%AREA3   , &
   iYearCurrent    => plt_site%iYearCurrent    , &
   NI      => plt_morph%NI     , &
   NGTopRootLayer     => plt_morph%NGTopRootLayer    , &
   MY      => plt_morph%MY     , &
   RootNonstructElementConcpft_vr  => plt_biom%RootNonstructElementConcpft_vr  , &
   CanopyNonstructElementConc_pft  => plt_biom%CanopyNonstructElementConc_pft  , &
   ShootChemElmnts_pft  => plt_biom%ShootChemElmnts_pft  , &
   CanPbndlResist     => plt_photo%CanPbndlResist    , &
   MinCanPStomaResistH2O    => plt_photo%MinCanPStomaResistH2O   , &
   CanPStomaResistH2O     => plt_photo%CanPStomaResistH2O    , &
   MaxCanPStomaResistH2O    => plt_photo%MaxCanPStomaResistH2O   , &
   RCS     => plt_photo%RCS    , &
   FracPARByCanP   => plt_rad%FracPARByCanP    , &
   LWRadCanP   => plt_rad%LWRadCanP      &
  )
  IF(NN.GE.MaxIterNum)THEN
    WRITE(*,9999)iYearCurrent,I,J,NZ
9999  FORMAT('CONVERGENCE FOR WATER UPTAKE NOT ACHIEVED ON   ',6I4)

    IF(DIFF.GT.0.5_r8)THEN
      plt_rad%RadNet2CanP(NZ)=0.0_r8
      plt_ew%EvapTransHeatP(NZ)=0.0_r8
      plt_ew%HeatXAir2PCan(NZ)=0.0_r8
      plt_ew%HeatStorCanP(NZ)=0.0_r8
      plt_ew%VapXAir2PCan(NZ)=0.0_r8
      plt_ew%PTrans(NZ)=0.0_r8
      TKC(NZ)=TairK+DTKC(NZ)
      TCelciusCanopy(NZ)=units%Kelvin2Celcius(TKC(NZ))
      FTHRM=EMMC*2.04E-10_r8*FracPARByCanP(NZ)*AREA3(NU)
      LWRadCanP(NZ)=FTHRM*TKC(NZ)**4._r8
      PSICanP(NZ)=ElvAdjstedtSoiPSIMPa(NGTopRootLayer(NZ))      
      CCPOLT=CanopyNonstructElementConc_pft(ielmc,NZ)+CanopyNonstructElementConc_pft(ielmn,NZ)+CanopyNonstructElementConc_pft(ielmp,NZ)

      CALL update_osmo_turg_pressure(PSICanP(NZ),CCPOLT,OSMO(NZ),TKC(NZ),PSICanPOsmo(NZ),PSICanPTurg(NZ),FDMP)

      Stomata_Activity=EXP(RCS(NZ)*PSICanPTurg(NZ))
      CanPStomaResistH2O(NZ)=MinCanPStomaResistH2O(NZ)+(MaxCanPStomaResistH2O(NZ)-MinCanPStomaResistH2O(NZ))*Stomata_Activity
      CanPbndlResist(NZ)=RAZ(NZ)
      VHeatCapCanP(NZ)=cpw*(ShootChemElmnts_pft(ielmc,NZ)*10.0E-06_r8)
      DTKC(NZ)=0.0_r8
      D4290: DO N=1,MY(NZ)
        DO  L=NU,NI(NZ)
          PSIRoot(N,L,NZ)=ElvAdjstedtSoiPSIMPa(L)
          CCPOLT=sum(RootNonstructElementConcpft_vr(1:NumOfPlantChemElements,N,L,NZ))
          CALL update_osmo_turg_pressure(PSIRoot(N,L,NZ),CCPOLT,OSMO(NZ),TKS(L),&
            PSIRootOSMO(N,L,NZ),PSIRootTurg(N,L,NZ))

          AllPlantRootH2OUptake_vr(N,L,NZ)=0.0_r8
      enddo
      ENDDO D4290
    ENDIF
  ENDIF
  end associate
  end subroutine HandlingDivergence
!------------------------------------------------------------------------------
  function CanopyEnergyH2OIteration(I,J,NZ,FracGrndByPFT,CanPMassC,&
    ElvAdjstedtSoiPSIMPa,HeatSensConductCanP,DIFF,cumPRootH2OUptake,VFLXC,FDMP,&
    SoiAddRootResist,FracPRoot4Uptake,AirPoreAvail4Fill,WatAvail4Uptake,TKCX,CNDT,&
    VHCPX,PrecHeatIntcptByCanP1,PSILH,LayrHasRoot) result(NN)
  implicit none
  integer  , intent(in) :: I, J
  integer  , intent(in) :: NZ
  real(r8) , intent(in) :: FracGrndByPFT,CanPMassC,ElvAdjstedtSoiPSIMPa(JZ1)
  real(r8) , intent(in) :: SoiAddRootResist(2,JZ1),FracPRoot4Uptake(2,JZ1,JP1),AirPoreAvail4Fill(JZ1),WatAvail4Uptake(JZ1)
  real(r8) , intent(in) :: TKCX,CNDT,VHCPX,PrecHeatIntcptByCanP1,PSILH
  integer  , intent(in) :: LayrHasRoot(2,JZ1)
  real(r8) , intent(out):: HeatSensConductCanP,DIFF,cumPRootH2OUptake,VFLXC,FDMP
  real(r8) :: APSILT
  real(r8) :: CCPOLT
  real(r8) :: DTHS1
  real(r8) :: DIFFZ,DIFFU,DPSI
  real(r8) :: EX,PTransPre
  real(r8) :: FTHRM,FDMR
  real(r8) :: LWRad2CanP
  real(r8) :: HeatSenAddStore
  real(r8) :: EffGrndAreaByPFT4H2O,EffGrndAreaByPFT4Heat,PSICanPPre,HeatLatentConductCanP
  real(r8) :: PSILC
  real(r8) :: RSSUX,RSSU,RA1,RSSZ
  real(r8) :: TKC1,TKCY
  real(r8) :: cumPRootH2OUptakePre
  real(r8) :: VOLWPX,VPC,SymplasmicWat
  real(r8) :: XC,Stomata_Activity
  real(r8) :: RichardsNO
  integer  :: IC,ICHK
!     return variables
  integer :: NN
!     local variables
  integer :: N,L
!     begin_execution

  associate(                          &
    NU      => plt_site%NU      , &
    AREA3   => plt_site%AREA3   , &
    PSICanPOsmo   => plt_ew%PSICanPOsmo     , &
    OSMO    => plt_ew%OSMO      , &
    RAZ     => plt_ew%RAZ       , &
    VPA     => plt_ew%VPA       , &
    TairK     => plt_ew%TairK       , &
    RIB     => plt_ew%RIB       , &
    PTrans     => plt_ew%PTrans       , &
    PSICanP   => plt_ew%PSICanP     , &
    CanWatP   => plt_ew%CanWatP     , &
    AllPlantRootH2OUptake_vr   => plt_ew%AllPlantRootH2OUptake_vr     , &
    TKCZ    => plt_ew%TKCZ      , &
    VapXAir2PCan   => plt_ew%VapXAir2PCan     , &
    VHeatCapCanP   => plt_ew%VHeatCapCanP     , &
    PSICanPTurg   => plt_ew%PSICanPTurg     , &
    WatByPCan   => plt_ew%WatByPCan     , &
    PrecIntcptByCanP    => plt_ew%PrecIntcptByCanP      , &
    EvapTransHeatP   => plt_ew%EvapTransHeatP     , &
    ZEROL   => plt_biom%ZEROL   , &
    ZEROP   => plt_biom%ZEROP   , &
    CanopyNonstructElementConc_pft  => plt_biom%CanopyNonstructElementConc_pft  , &
    NI      => plt_morph%NI     , &
    MY      => plt_morph%MY     , &
    MinCanPStomaResistH2O    => plt_photo%MinCanPStomaResistH2O   , &
    RCS     => plt_photo%RCS    , &
    CanPbndlResist     => plt_photo%CanPbndlResist    , &
    MaxCanPStomaResistH2O    => plt_photo%MaxCanPStomaResistH2O   , &
    CanPStomaResistH2O     => plt_photo%CanPStomaResistH2O    , &
    RadNet2CanP    => plt_rad%RadNet2CanP     , &
    SWRadByCanP    => plt_rad%SWRadByCanP     , &
    LWRadSky     => plt_rad%LWRadSky      , &
    FracPARByCanP   => plt_rad%FracPARByCanP    , &
    LWRadCanP   => plt_rad%LWRadCanP    , &
    LWRadGrnd  => plt_rad%LWRadGrnd     &
  )
  CCPOLT=CanopyNonstructElementConc_pft(ielmc,NZ)+CanopyNonstructElementConc_pft(ielmn,NZ)+CanopyNonstructElementConc_pft(ielmp,NZ)
  FTHRM=EMMC*2.04E-10_r8*FracPARByCanP(NZ)*AREA3(NU)
  LWRad2CanP=(LWRadSky+LWRadGrnd)*FracPARByCanP(NZ)
!     RAZ=canopy isothermal boundary later resistance

  cumPRootH2OUptake=0.0_r8
  EffGrndAreaByPFT4H2O=FracGrndByPFT*AREA3(NU)            !m2
  EffGrndAreaByPFT4Heat=FracGrndByPFT*AREA3(NU)*1.25E-03_r8            !
  RA1=RAZ(NZ)

  IC=0
  XC=0.5_r8   !learning/updating rate for canopy temperature, not used now
  ICHK=0
  PSICanPPre=0.0_r8
  PTransPre=0.0_r8
  cumPRootH2OUptakePre=0.0_r8
  VOLWPX=0.0_r8

  DO NN=1,MaxIterNum
!
!     NET RADIATION FROM ABSORBED SW AND NET LW
!
!     LWRadCanP=LW emitted by canopy
!     DTHS1=net LW absorbed by canopy
!     SWRadByCanP=total SW absorbed by canopy
!     RadNet2CanP=net SW+LW absorbed by canopy
!
    TKC1=TKCZ(NZ)
    LWRadCanP(NZ)=FTHRM*TKC1**4._r8   !long wave radiation
    DTHS1=LWRad2CanP-LWRadCanP(NZ)*2.0_r8  !canopy long wave radiation goes both upper and down
    RadNet2CanP(NZ)=SWRadByCanP(NZ)+DTHS1
!
!     BOUNDARY LAYER RESISTANCE FROM RICHARDSON NUMBER
!
!     RI=Ricardson's number
!     RA=canopy boundary layer resistance
!     HeatLatentConductCanP,HeatSensConductCanP=canopy latent,sensible heat conductance
!
    RichardsNO=RichardsonNumber(RIB,TairK,TKC1)

    CanPbndlResist(NZ)=AMAX1(MinCanPbndlResist,0.9_r8*RA1,AMIN1(1.1_r8*RA1,RAZ(NZ)/(1.0_r8-10.0_r8*RichardsNO)))

    RA1=CanPbndlResist(NZ)
    HeatLatentConductCanP=EffGrndAreaByPFT4H2O/CanPbndlResist(NZ)  !m2/(h/m)=m/h
    HeatSensConductCanP=EffGrndAreaByPFT4Heat/CanPbndlResist(NZ)
!
!     CANOPY WATER AND OSMOTIC POTENTIALS
!
!     PSICanP=canopy total water potential
!     OSMO=osmotic potential at PSICanP=0 from PFT file
!     PSICanPOsmo,PSICanPTurg=canopy osmotic,turgor water potential
!
     call update_osmo_turg_pressure(PSICanP(NZ),CCPOLT,OSMO(NZ),TKC1,PSICanPOsmo(NZ),PSICanPTurg(NZ),FDMP)
!
!     CANOPY STOMATAL RESISTANCE
!
!     RCS=shape parameter for CanPStomaResistH2Ovs PSICanPTurg from PFT file
!     RC=canopy stomatal resistance
!     MinCanPStomaResistH2O=minimum CanPStomaResistH2Oat PSICanP=0 from stomate.f
!     RSMX=cuticular resistance from PFT file
!
    Stomata_Activity=EXP(RCS(NZ)*PSICanPTurg(NZ))

    CanPStomaResistH2O(NZ)=MinCanPStomaResistH2O(NZ)+Stomata_Activity &
      *(MaxCanPStomaResistH2O(NZ)-MinCanPStomaResistH2O(NZ))
!
!     CANOPY VAPOR PRESSURE AND EVAPORATION OF INTERCEPTED WATER
!     OR TRANSPIRATION OF UPTAKEN WATER
!
!     VPC,VPA=vapor pressure inside canopy, in atmosphere
!     TKC1=canopy temperature
!     PSICanP=canopy total water potential
!     EX=canopy-atmosphere water flux
!     RA,RZ=canopy boundary layer,surface resistance
!     VapXAir2PCan=water flux to,from air to canopy surfaces
!     WatByPCan=water volume on canopy surfaces
!     EP=water flux to,from inside canopy
!     EvapTransHeatP=canopy latent heat flux
!     VFLXC=convective heat flux from EvapTransHeatP
!     VAP=latent heat of evaporation
!     HeatLatentConductCanP=aerodynamic conductance

    VPC=vapsat(tkc1)*EXP(18.0_r8*PSICanP(NZ)/(RGAS*TKC1))
    EX=HeatLatentConductCanP*(VPA-VPC)   !air to canopy water vap flux, [m/h]*[ton/m3]
    IF(EX.GT.0.0_r8)THEN
      !condensation  > 0._r8, to canopy
      VapXAir2PCan(NZ)=EX*CanPbndlResist(NZ)/(CanPbndlResist(NZ)+RZ)
      EX=0.0_r8
      VFLXC=VapXAir2PCan(NZ)*cpw*TairK                !enthalpy of evaporated water, MJ/(h*m2)
    ELSEIF(EX.LE.0.0_r8.AND.WatByPCan(NZ).GT.0.0_r8)THEN
      !evaporation, and there is water stored in canopy
      !<0._r8, off canopy, cannot be more than WatByPCan(NZ)
      VapXAir2PCan(NZ)=AMAX1(EX*CanPbndlResist(NZ)/(CanPbndlResist(NZ)+RZ),-WatByPCan(NZ))
      EX=EX-VapXAir2PCan(NZ)                !extra water needs transpiration
      VFLXC=VapXAir2PCan(NZ)*cpw*TKC1       !enthalpy of evaporated water
    ENDIF
    !PTrans <0 means transpiration into atmosphere
    PTrans(NZ)=EX*CanPbndlResist(NZ)/(CanPbndlResist(NZ)+CanPStomaResistH2O(NZ))
    EvapTransHeatP(NZ)=(PTrans(NZ)+VapXAir2PCan(NZ))*EvapLHTC   !latent heat flux, negative means into atmosphere

!
!     SENSIBLE + STORAGE HEAT FROM RN, LE AND CONVECTIVE HEAT FLUXES
!
!     HeatSenAddStore=initial estimate of sensible+storage heat flux
!     PrecHeatIntcptByCanP1=convective heat flux from precip to canopy
!
    HeatSenAddStore=RadNet2CanP(NZ)+EvapTransHeatP(NZ)+VFLXC+PrecHeatIntcptByCanP1
!
!     SOLVE FOR CANOPY TEMPERATURE CAUSED BY SENSIBLE + STORAGE HEAT
!
!     VHeatCapCanP=canopy heat capacity
!     TKCY=equilibrium canopy temperature for HeatSenAddStore
!
!   VHCPX= canopy heat capacity
    VHeatCapCanP(NZ)=VHCPX+cpw*(VapXAir2PCan(NZ)+PrecIntcptByCanP(NZ))
    TKCY=(TKCX*VHCPX+TairK*HeatSensConductCanP+HeatSenAddStore)/(VHeatCapCanP(NZ)+HeatSensConductCanP)

    !canopy temperature not differ from air more more than 10 K
    TKCY=AMIN1(TairK+10.0_r8,AMAX1(TairK-10.0_r8,TKCY))
!
!     RESET CANOPY TEMPERATURE FOR NEXT ITERATION
!
!     XC,IC=magnitude,direction of change in canopy temp for next cycle
!
    IF((IC.EQ.0.AND.TKCY.GT.TKC1).OR.(IC.EQ.1.AND.TKCY.LT.TKC1))THEN
      XC=0.5_r8*XC
    ENDIF
    !0.1 is the learning/updating rate
    TKCZ(NZ)=TKC1+0.1_r8*(TKCY-TKC1)
    IF(TKCY.GT.TKC1)THEN
      !warming
      IC=1
    ELSE
      !cooling
      IC=0
    ENDIF

!
!     IF CONVERGENCE CRITERION IS MET OR ON EVERY TENTH ITERATION,
!     PROCEED TO WATER BALANCE
!
!     PSILC=canopy water potential adjusted for canopy height
!
    IF(ABS(TKCY-TKC1).LT.0.05_r8.OR.(NN/10)*10.EQ.NN)THEN
      cumPRootH2OUptake=0.0_r8
      PSILC=PSICanP(NZ)-PSILH
!
!     ROOT WATER UPTAKE FROM SOIL-CANOPY WATER POTENTIALS,
!     SOIL + ROOT HYDRAULIC RESISTANCES
!
!     LayrHasRoot=rooted layer flag
!     AllPlantRootH2OUptake_vr=root water uptake from soil layer > 0
!     WatAvail4Uptake,AirPoreAvail4Fill=water volume available for uptake,air volume
!     FracPRoot4Uptake=PFT fraction of biome root mass
!     PSILC=canopy water potential adjusted for canopy height
!     ElvAdjstedtSoiPSIMPa=total soil water potential PSIST adjusted for surf elevn
!     SoiAddRootResist=total soil+root resistance
!     cumPRootH2OUptake=total water uptake from soil profile
!
      DO N=1,MY(NZ)
        DO  L=NU,NI(NZ)
          IF(LayrHasRoot(N,L).EQ.1)THEN
            AllPlantRootH2OUptake_vr(N,L,NZ)=AMAX1(AZMIN1(-WatAvail4Uptake(L)*FracPRoot4Uptake(N,L,NZ)) &
              ,AMIN1((PSILC-ElvAdjstedtSoiPSIMPa(L))/SoiAddRootResist(N,L),AirPoreAvail4Fill(L)*FracPRoot4Uptake(N,L,NZ)))

            IF(AllPlantRootH2OUptake_vr(N,L,NZ).GT.0.0_r8)THEN
              !plant/myco lose water to soil
              AllPlantRootH2OUptake_vr(N,L,NZ)=0.1_r8*AllPlantRootH2OUptake_vr(N,L,NZ)
            ENDIF
            cumPRootH2OUptake=cumPRootH2OUptake+AllPlantRootH2OUptake_vr(N,L,NZ)
          ELSE
            AllPlantRootH2OUptake_vr(N,L,NZ)=0.0_r8
          ENDIF
        enddo
      ENDDO
!
!     TEST TRANSPIRATION - ROOT WATER UPTAKE VS. CHANGE IN CANOPY
!     WATER STORAGE
!
!     SymplasmicWat,CanWatP=canopy water content
!     DIFFZ,DIFFU=change in canopy water content,transpiration-uptake
!     DIFF=normalized difference between DIFFZ and DIFFU
!     5.0E-03=acceptance criterion for DIFF
!     RSSZ=change in canopy water potl vs change in canopy water cnt
!     RSSU=change in canopy water potl vs change in transpiration

      SymplasmicWat=ppmc*CanPMassC/FDMP
      DIFFZ=SymplasmicWat-CanWatP(NZ)
      DIFFU=PTrans(NZ)-cumPRootH2OUptake
      IF(.not.isclose(cumPRootH2OUptake,0.0_r8))THEN
        DIFF=ABS((DIFFU-DIFFZ)/cumPRootH2OUptake)
      ELSE
        DIFF=ABS((DIFFU-DIFFZ)/SymplasmicWat)
      ENDIF
      IF(DIFF.LT.5.0E-03_r8)THEN
        IF(ICHK.EQ.1)EXIT
        ICHK=1
        CALL STOMATEs(I,J,NZ)
        CYCLE
      ENDIF
      IF(ABS(SymplasmicWat-VOLWPX).GT.ZEROP(NZ))THEN
        RSSZ=ABS((PSICanP(NZ)-PSICanPPre)/(SymplasmicWat-VOLWPX))
      ELSEIF(CNDT.GT.ZEROP(NZ))THEN
        RSSZ=1.0_r8/CNDT
      ELSE
        RSSZ=ZEROL(NZ)
      ENDIF
      IF(ABS(PTrans(NZ)-PTransPre).GT.ZEROP(NZ))THEN
        RSSUX=ABS((PSICanP(NZ)-PSICanPPre)/(PTrans(NZ)-PTransPre))
        IF(CNDT.GT.ZEROP(NZ))THEN
          RSSU=AMIN1(1.0_r8/CNDT,RSSUX)
        ELSE
          RSSU=RSSUX
        ENDIF
      ELSEIF(ABS(cumPRootH2OUptake-cumPRootH2OUptakePre).GT.ZEROP(NZ))THEN
        RSSUX=ABS((PSICanP(NZ)-PSICanPPre)/(cumPRootH2OUptake-cumPRootH2OUptakePre))
        IF(CNDT.GT.ZEROP(NZ))THEN
          RSSU=AMIN1(1.0_r8/CNDT,RSSUX)
        ELSE
          RSSU=RSSUX
        ENDIF
      ELSEIF(CNDT.GT.ZEROP(NZ))THEN
        RSSU=1.0_r8/CNDT
      ELSE
        RSSU=ZEROL(NZ)
      ENDIF
!
!     CHANGE IN CANOPY WATER POTENTIAL REQUIRED TO BRING AGREEMENT
!     BETWEEN TRANSPIRATION - ROOT WATER UPTAKE AND CHANGE IN CANOPY
!     WATER STORAGE
!
!     DPSI=change in PSICanP for next convergence cycle
!     1.0E-03=acceptance criterion for DPSI
!
      DPSI=AMIN1(AMIN1(RSSZ,RSSU)*(DIFFU-DIFFZ),ABS(PSICanP(NZ)))

!     IF CONVERGENCE CRITERION IS MET THEN FINISH,
!     OTHERWISE START NEXT ITERATION WITH CANOPY WATER POTENTIAL
!     TRANSPIRATION, UPTAKE AND WATER CONTENT FROM CURRENT ITERATION
!
      IF(.not.((NN.GE.30.AND.ABS(DPSI).LT.1.0E-03_r8).OR.NN.GE.MaxIterNum))then
        !make copies for next iteration
        PSICanPPre=PSICanP(NZ)
        PTransPre=PTrans(NZ)
        cumPRootH2OUptakePre=cumPRootH2OUptake
        VOLWPX=SymplasmicWat
        PSICanP(NZ)=AZMIN1(PSICanP(NZ)+0.5_r8*DPSI)

        XC=0.50_r8
        cycle
!
!     RESET MIN STOMATAL RESISTANCE IN STOMATE.F BEFORE FINAL ITERATION
!
      ELSE
        IF(ICHK.EQ.1)EXIT
        ICHK=1
        CALL STOMATEs(I,J,NZ)
      ENDIF
    ENDIF
  ENDDO
  end associate
  end function CanopyEnergyH2OIteration
!------------------------------------------------------------------------
  subroutine CalcResistance(NZ,PATH,FineRootRadius,RootAreaDivRadius,&
    RootResist,RootResistSoi,RootResistRadial,RootResistAxial,SoiH2OResist,&
    SoiAddRootResist,CNDT,PSILH,LayrHasRoot)

  implicit none
  integer, intent(in)   :: NZ
  real(r8), intent(in)  :: PATH(2,JZ1),FineRootRadius(2,JZ1),RootAreaDivRadius(2,JZ1)
  real(r8), intent(out) :: RootResist(2,JZ1),RootResistSoi(2,JZ1)
  real(r8), intent(out) :: RootResistRadial(2,JZ1),RootResistAxial(2,JZ1)
  real(r8), intent(out) :: SoiH2OResist(2,JZ1),SoiAddRootResist(2,JZ1),CNDT
  real(r8), intent(out) :: PSILH
  integer, intent(out) :: LayrHasRoot(2,JZ1)
  real(r8) :: FRADW,FRAD1,FRAD2
  real(r8) :: RSSL,RTAR2
  integer :: N, L
  associate(                          &
    DPTHZ  => plt_site%DPTHZ                            , &
    pftPlantPopulation     => plt_site%pftPlantPopulation       , &
    VLWatMicPM  => plt_site%VLWatMicPM    , &
    ZERO   => plt_site%ZERO     , &
    ZEROS2 => plt_site%ZEROS2   , &
    NU     => plt_site%NU       , &
    PSICanP  => plt_ew%PSICanP      , &
    ZEROP  => plt_biom%ZEROP    , &
    THETW  => plt_soilchem%THETW, &
    VLMicP   => plt_soilchem%VLMicP , &
    HydroCondMicP4RootUptake   => plt_soilchem%HydroCondMicP4RootUptake , &
    VLSoilPoreMicP   => plt_soilchem%VLSoilPoreMicP , &
    SecndRootXNumL   => plt_morph%SecndRootXNumL    , &
    PrimRootXNumL  => plt_morph%PrimRootXNumL   , &
    MaxSecndRootRadius => plt_morph%MaxSecndRootRadius  , &
    RSRA   => plt_morph%RSRA    , &
    SecndRootRadius  => plt_morph%SecndRootRadius   , &
    RootLenDensNLP  => plt_morph%RootLenDensNLP   , &
    CanPHeight4WatUptake  => plt_morph%CanPHeight4WatUptake   , &
    RootLenPerP  => plt_morph%RootLenPerP   , &
    AveSecndRootLen  => plt_morph%AveSecndRootLen   , &
    PrimRootRadius  => plt_morph%PrimRootRadius   , &
    MY     => plt_morph%MY      , &
    RSRR   => plt_morph%RSRR    , &
    CanopyHeight     => plt_morph%CanopyHeight      , &
    NGTopRootLayer    => plt_morph%NGTopRootLayer     , &
    NI     => plt_morph%NI        &
  )

  !     GRAVIMETRIC WATER POTENTIAL FROM CANOPY HEIGHT
  !
  !     CanPHeight4WatUptake=canopy height for water uptake
  !     PSILH=gravimetric water potential at CanPHeight4WatUptake
  !     FRADW=conducting elements of stalk relative to those of primary root
  !     PSICanP=canopy total water potential
  !     EMODW=wood modulus of elasticity (MPa)
!
  CNDT=0.0_r8
  CanPHeight4WatUptake(NZ)=0.80_r8*CanopyHeight(NZ)
  PSILH=-mGravAccelerat*CanPHeight4WatUptake(NZ)
  FRADW=1.0E+04_r8*(AMAX1(0.5_r8,1.0_r8+PSICanP(NZ)/EMODW))**4._r8
!
  !     SOIL AND ROOT HYDRAULIC RESISTANCES TO ROOT WATER UPTAKE
  !
  !      VLSoilPoreMicP,VLWatMicPM,THETW=soil,water volume,content
  !     RootLenDensNLP,RootLenPerP=root length density,root length per plant
  !     HydroCondMicP4RootUptake=soil hydraulic conductivity for root uptake
  !     PrimRootXNumL,SecndRootXNumL=number of root,myco primary,secondary axes
  !     LayrHasRoot:1=rooted,0=not rooted
  !     N:1=root,2=mycorrhizae
!
  D3880: DO N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      IF(VLSoilPoreMicP(L).GT.ZEROS2 &
        .AND.VLWatMicPM(NPH,L).GT.ZEROS2 &
        .AND.RootLenDensNLP(N,L,NZ).GT.ZERO &
        .AND.HydroCondMicP4RootUptake(L).GT.ZERO &
        .AND.PrimRootXNumL(ipltroot,L,NZ).GT.ZEROP(NZ) &
        .AND.SecndRootXNumL(N,L,NZ).GT.ZEROP(NZ) &
        .AND.THETW(L).GT.ZERO)THEN
        LayrHasRoot(N,L)=1
        !
        !     SOIL HYDRAULIC RESISTANCE FROM RADIAL UPTAKE GEOMETRY
        !     AND SOIL HYDRAULIC CONDUCTIVITY
        !
        !     SoiH2OResist=soil hydraulic resistance
        !     PP=plant population
        !     PATH=path length of water and nutrient uptake
        !     FineRootRadius,RootAreaDivRadius=root radius,surface/radius area
        !
        RSSL=(LOG(PATH(N,L)/FineRootRadius(N,L))/RootAreaDivRadius(N,L))/pftPlantPopulation(NZ)
        SoiH2OResist(N,L)=RSSL/HydroCondMicP4RootUptake(L)
        !
        !     RADIAL ROOT RESISTANCE FROM ROOT AREA AND RADIAL RESISTIVITY
        !     ENTERED IN 'READQ'
        !
        !     SecndRootRadius=secondary root radius
        !     RootLenPerP=root length per plant
        !     RootResistSoi=radial resistance
        !     RSRR=radial resistivity from PFT file
        !     VLMicP,VLWatMicPM=soil micropore,water volume
        !
        RTAR2=TwoPiCON*SecndRootRadius(N,L,NZ)*RootLenPerP(N,L,NZ)*pftPlantPopulation(NZ)
        RootResistSoi(N,L)=RSRR(N,NZ)/RTAR2*VLMicP(L)/VLWatMicPM(NPH,L)
!
        !     ROOT AXIAL RESISTANCE FROM RADII AND LENGTHS OF PRIMARY AND
        !     SECONDARY ROOTS AND FROM AXIAL RESISTIVITY ENTERED IN 'READQ'
        !
        !     FRAD1,FRAD2=primary,secondary root radius relative to maximum
        !     secondary radius from PFT file MaxSecndRootRadius at which RSRA is defined
        !     PrimRootRadius,SecndRootRadius=primary,secondary root radius
        !     RSRA=axial resistivity from PFT file
        !     DPTHZ=depth of primary root from surface
        !     RootResistRadial,RootResistAxial=axial resistance of primary,secondary roots
        !     AveSecndRootLen=average secondary root length
        !     PrimRootXNumL,SecndRootXNumL=number of primary,secondary axes
!
        FRAD1=(PrimRootRadius(N,L,NZ)/MaxSecndRootRadius(N,NZ))**4._r8
        RootResistRadial(N,L)=RSRA(N,NZ)*DPTHZ(L)/(FRAD1*PrimRootXNumL(ipltroot,L,NZ)) &
          +RSRA(ipltroot,NZ)*CanPHeight4WatUptake(NZ)/(FRADW*PrimRootXNumL(ipltroot,L,NZ))
        FRAD2=(SecndRootRadius(N,L,NZ)/MaxSecndRootRadius(N,NZ))**4._r8
        RootResistAxial(N,L)=RSRA(N,NZ)*AveSecndRootLen(N,L,NZ)/(FRAD2*SecndRootXNumL(N,L,NZ))
      ELSE
        LayrHasRoot(N,L)=0
      ENDIF
    enddo
  ENDDO D3880

  D3890: DO N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      IF(LayrHasRoot(N,L).EQ.1)THEN
        !
        !     TOTAL ROOT RESISTANCE = SOIL + RADIAL + AXIAL
        !
        !     RootResist=root radial+axial resistance
        !     SoiAddRootResist=total soil+root resistance
        !     CNDT=total soil+root conductance for all layers
        !
        RootResist(N,L)=RootResistSoi(N,L)+RootResistRadial(N,L)+RootResistAxial(N,L)
        SoiAddRootResist(N,L)=SoiH2OResist(N,L)+RootResist(N,L)
        CNDT=CNDT+1.0/SoiAddRootResist(N,L)

      ENDIF
    enddo
  ENDDO D3890
  end associate
  end subroutine CalcResistance
!------------------------------------------------------------------------

  subroutine HandleBareSoil(NZ,ElvAdjstedtSoiPSIMPa,FDMP)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: ElvAdjstedtSoiPSIMPa(JZ1)
  real(r8), intent(out):: FDMP
  integer :: N,L
  real(r8) :: APSILT
  real(r8) :: CCPOLT
  real(r8) :: FTHRM,FDMR
  real(r8) :: OSWT,Stomata_Activity

! begin_execution
  associate(                         &
    TCelciusCanopy   =>  plt_ew%TCelciusCanopy      , &
    OSMO   =>  plt_ew%OSMO      , &
    TKSnow    =>  plt_ew%TKSnow       , &
    TKC    =>  plt_ew%TKC       , &
    TKS    =>  plt_ew%TKS       , &
    TairK    =>  plt_ew%TairK       , &
    SnowDepth  =>  plt_ew%SnowDepth     , &
    DTKC   =>  plt_ew%DTKC      , &
    RAZ    =>  plt_ew%RAZ       , &
    VHeatCapCanP  =>  plt_ew%VHeatCapCanP     , &
    PSICanPOsmo  =>  plt_ew%PSICanPOsmo     , &
    PSIRootTurg  =>  plt_ew%PSIRootTurg     , &
    AllPlantRootH2OUptake_vr  =>  plt_ew%AllPlantRootH2OUptake_vr     , &
    PSICanPTurg  =>  plt_ew%PSICanPTurg     , &
    PSICanP  =>  plt_ew%PSICanP     , &
    PSIRootOSMO  =>  plt_ew%PSIRootOSMO     , &
    PSIRoot  =>  plt_ew%PSIRoot     , &
    PTrans    =>  plt_ew%PTrans       , &
    HeatStorCanP  =>  plt_ew%HeatStorCanP     , &
    EvapTransHeatP  =>  plt_ew%EvapTransHeatP     , &
    HeatXAir2PCan  =>  plt_ew%HeatXAir2PCan     , &
    VapXAir2PCan  =>  plt_ew%VapXAir2PCan     , &
    NU     =>  plt_site%NU      , &
    ZERO   =>  plt_site%ZERO    , &
    AREA3  =>  plt_site%AREA3   , &
    CanopyNonstructElementConc_pft =>  plt_biom%CanopyNonstructElementConc_pft  , &
    RootNonstructElementConcpft_vr =>  plt_biom%RootNonstructElementConcpft_vr  , &
    ShootChemElmnts_pft =>  plt_biom%ShootChemElmnts_pft  , &
    NI     =>  plt_morph%NI     , &
    CanopyHeight     =>  plt_morph%CanopyHeight     , &
    NGTopRootLayer    =>  plt_morph%NGTopRootLayer    , &
    MY     =>  plt_morph%MY     , &
    RCS    =>  plt_photo%RCS    , &
    CanPbndlResist    =>  plt_photo%CanPbndlResist    , &
    CanPStomaResistH2O    =>  plt_photo%CanPStomaResistH2O    , &
    MinCanPStomaResistH2O   =>  plt_photo%MinCanPStomaResistH2O   , &
    MaxCanPStomaResistH2O   =>  plt_photo%MaxCanPStomaResistH2O   , &
    FracPARByCanP  =>  plt_rad%FracPARByCanP    , &
    LWRadCanP  =>  plt_rad%LWRadCanP    , &
    RadNet2CanP   =>  plt_rad%RadNet2CanP       &
  )
  RadNet2CanP(NZ)=0.0_r8
  EvapTransHeatP(NZ)=0.0_r8
  HeatXAir2PCan(NZ)=0.0_r8
  HeatStorCanP(NZ)=0.0_r8
  VapXAir2PCan(NZ)=0.0_r8
  PTrans(NZ)=0.0_r8
  IF(CanopyHeight(NZ).GE.SnowDepth-ZERO)THEN
    TKC(NZ)=TairK
  ELSE
    TKC(NZ)=TKSnow
  ENDIF
  TCelciusCanopy(NZ)=units%Kelvin2Celcius(TKC(NZ))
  FTHRM=EMMC*2.04E-10_r8*FracPARByCanP(NZ)*AREA3(NU)
  LWRadCanP(NZ)=FTHRM*TKC(NZ)**4._r8
  PSICanP(NZ)=ElvAdjstedtSoiPSIMPa(NGTopRootLayer(NZ))
  
  CCPOLT=sum(CanopyNonstructElementConc_pft(1:NumOfPlantChemElements,NZ))

  call update_osmo_turg_pressure(PSICanP(NZ),CCPOLT,OSMO(NZ),TKC(NZ),PSICanPOsmo(NZ),PSICanPTurg(NZ),FDMP)

  Stomata_Activity=EXP(RCS(NZ)*PSICanPTurg(NZ))
  CanPStomaResistH2O(NZ)=MinCanPStomaResistH2O(NZ)+(MaxCanPStomaResistH2O(NZ)-MinCanPStomaResistH2O(NZ))*Stomata_Activity
  CanPbndlResist(NZ)=RAZ(NZ)
  VHeatCapCanP(NZ)=cpw*(ShootChemElmnts_pft(ielmc,NZ)*10.0E-06_r8)
  DTKC(NZ)=0.0_r8

  DO N=1,MY(NZ)
    DO  L=NU,NI(NZ)
      PSIRoot(N,L,NZ)=ElvAdjstedtSoiPSIMPa(L)      
      CCPOLT=sum(RootNonstructElementConcpft_vr(1:NumOfPlantChemElements,N,L,NZ))

      call update_osmo_turg_pressure(PSIRoot(N,L,NZ),CCPOLT,OSMO(NZ),TKS(L),&
        PSIRootOSMO(N,L,NZ),PSIRootTurg(N,L,NZ))

      AllPlantRootH2OUptake_vr(N,L,NZ)=0.0_r8
    enddo
  ENDDO
  end associate
  end subroutine HandleBareSoil
!------------------------------------------------------------------------

  subroutine UpdateCanopyWater(NZ,HeatSensConductCanP,ElvAdjstedtSoiPSIMPa,RootResist,SoiH2OResist,&
    SoiAddRootResist,TKCX,VHCPX,PrecHeatIntcptByCanP1,cumPRootH2OUptake,VFLXC,LayrHasRoot)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: HeatSensConductCanP,ElvAdjstedtSoiPSIMPa(JZ1),RootResist(2,JZ1)
  real(r8), intent(in) :: SoiH2OResist(2,JZ1),SoiAddRootResist(2,JZ1)
  real(r8), intent(in) :: TKCX,VHCPX,PrecHeatIntcptByCanP1,cumPRootH2OUptake,VFLXC
  integer , intent(in) :: LayrHasRoot(2,JZ1)
  real(r8) :: CCPOLT
  real(r8) :: FDMR
  real(r8) :: OSWT
  integer :: N,L
  associate(                           &
    NU       => plt_site%NU      , &
    OSMO     => plt_ew%OSMO      , &
    TairK      => plt_ew%TairK       , &
    TKC      => plt_ew%TKC       , &
    TKS      => plt_ew%TKS       , &
    TKCZ     => plt_ew%TKCZ      , &
    WatByPCan    => plt_ew%WatByPCan     , &
    CanWatP   => plt_ew%CanWatP    , &
    VHeatCapCanP    => plt_ew%VHeatCapCanP     , &
    HeatStorCanP    => plt_ew%HeatStorCanP     , &
    PrecIntcptByCanP     => plt_ew%PrecIntcptByCanP      , &
    PTrans      => plt_ew%PTrans       , &
    VapXAir2PCan    => plt_ew%VapXAir2PCan     , &
    PSIRoot    => plt_ew%PSIRoot     , &
    HeatXAir2PCan    => plt_ew%HeatXAir2PCan     , &
    PSIRootOSMO    => plt_ew%PSIRootOSMO     , &
    PSIRootTurg    => plt_ew%PSIRootTurg     , &
    PSICanP    => plt_ew%PSICanP     , &
    RootNonstructElementConcpft_vr   => plt_biom%RootNonstructElementConcpft_vr  , &
    MY       => plt_morph%MY     , &
    NI       => plt_morph%NI       &
  )
  !
  !     CANOPY SURFACE WATER STORAGE, SENSIBLE AND STORAGE HEAT FLUXES
  !     (NOT EXPLICITLY CALCULATED IN CONVERGENCE SOLUTION)
  !
  !     VOLWP,WatByPCan=water volume in canopy,on canopy surfaces
  !     HeatXAir2PCan,HeatStorCanP=canopy sensible,storage heat fluxes
  !     VHCPX,VHeatCapCanP=previous,current canopy heat capacity
  !     HeatSensConductCanP=canopy sensible heat conductance
  !     VFLXC=convective heat flux from latent heat flux
  !     PrecHeatIntcptByCanP1=convective heat flux from precip to canopy
  !
  CanWatP(NZ)=CanWatP(NZ)+PTrans(NZ)-cumPRootH2OUptake
  WatByPCan(NZ)=WatByPCan(NZ)+PrecIntcptByCanP(NZ)+VapXAir2PCan(NZ)
  HeatXAir2PCan(NZ)=HeatSensConductCanP*(TairK-TKCZ(NZ))
  HeatStorCanP(NZ)=TKCX*VHCPX-TKCZ(NZ)*VHeatCapCanP(NZ)+VFLXC+PrecHeatIntcptByCanP1
  !
  !     ROOT TOTAL, OSMOTIC AND TURGOR WATER POTENTIALS
  !
  !     PSIRoot,PSICanP=root,canopy total water potential
  !     ElvAdjstedtSoiPSIMPa=total soil water potential PSIST adjusted for surf elevn
  !     SoiH2OResist,SoiAddRootResist,RootResist=soil,soil+root,root radial+axial resistance
  !     PSIRootOSMO,PSIRootTurg=root osmotic,turgor water potential
  !     FDMR=dry matter content
  !     OSMO=osmotic potential at PSIRoot=0 from PFT file
!
  !
  D4505: DO N=1,MY(NZ)
    D4510: DO L=NU,NI(NZ)
      IF(LayrHasRoot(N,L).EQ.1)THEN
        PSIRoot(N,L,NZ)=AZMIN1((ElvAdjstedtSoiPSIMPa(L)*RootResist(N,L) &
          +PSICanP(NZ)*SoiH2OResist(N,L))/SoiAddRootResist(N,L))
      ELSE
        PSIRoot(N,L,NZ)=ElvAdjstedtSoiPSIMPa(L)
      ENDIF           
      CCPOLT=sum(RootNonstructElementConcpft_vr(1:NumOfPlantChemElements,N,L,NZ))

      CALL update_osmo_turg_pressure(PSIRoot(N,L,NZ),CCPOLT,OSMO(NZ),TKS(L),&
        PSIRootOSMO(N,L,NZ),PSIRootTurg(N,L,NZ))

    ENDDO D4510
  ENDDO D4505
  end associate
  end subroutine UpdateCanopyWater
!------------------------------------------------------------------------

  subroutine SetCanopyGrowthFuncs(NZ)

  implicit none
  integer, intent(in) :: NZ
  real(r8) :: ACTV,RTK,STK,TKGO,TKSO
  integer :: L
  associate(                          &
    TCelciusCanopy   =>  plt_ew%TCelciusCanopy      , &
    TKC    =>  plt_ew%TKC       , &
    TKS    =>  plt_ew%TKS       , &
    PSICanPDailyMin  =>  plt_ew%PSICanPDailyMin     , &
    PSICanP  =>  plt_ew%PSICanP     , &
    NU     =>  plt_site%NU      , &
    CHILL  =>  plt_photo%CHILL  , &
    OFFST  =>  plt_pheno%OFFST  , &
    TCelciusChill4Seed   =>  plt_pheno%TCelciusChill4Seed   , &
    fTgrowRootP   =>  plt_pheno%fTgrowRootP   , &
    TCG    =>  plt_pheno%TCG    , &
    TKG    =>  plt_pheno%TKG    , &
    iPlantCalendar  =>  plt_pheno%iPlantCalendar  , &
    fTgrowCanP   =>  plt_pheno%fTgrowCanP   , &
    NI     =>  plt_morph%NI     , &
    NumOfMainBranch_pft    =>  plt_morph%NumOfMainBranch_pft      &
  )
  !
  !     SET CANOPY GROWTH TEMPERATURE FROM SOIL SURFACE
  !     OR CANOPY TEMPERATURE DEPENDING ON GROWTH STAGE
  !
  IF(iPlantCalendar(ipltcal_Emerge,NumOfMainBranch_pft(NZ),NZ).EQ.0)THEN
    TKG(NZ)=TKS(NU)
    !     ELSEIF((iPlantTurnoverPattern(NZ).EQ.0.OR.iPlantMorphologyType(NZ).LE.1)
    !    2.AND.iPlantCalendar(ipltcal_InitFloral,NumOfMainBranch_pft(NZ),NZ).EQ.0)THEN
    !     TKG(NZ)=TKS(NU)
  ELSE
    TKG(NZ)=TKC(NZ)
  ENDIF
  TCG(NZ)=units%Kelvin2Celcius(TKG(NZ))
  !
  !     ARRHENIUS FUNCTION FOR CANOPY AND ROOT GROWTH WITH OFFSET
  !     FOR ZONE OF THERMAL ADAPTATION ENTERED IN 'READQ'
  !
  !     TKG,TKGO=canopy temperature,canopy temp used in Arrhenius eqn
  !     TKS,TKSO=soil temperature,soil temp used in Arrhenius eqn
  !     OFFST=shift in Arrhenius curve for thermal adaptation
  !     fTgrowCanP,fTgrowRootP=temperature function for canopy,root growth (25 oC =1)
  !     RGAS,710.0=gas constant,enthalpy
  !     62500,197500,222500=energy of activn,high,low temp inactivn(KJ mol-1)
  !     PSICanPDailyMin=minimum daily canopy water potential
  !
  TKGO=TKG(NZ)+OFFST(NZ)
  fTgrowCanP(NZ)=calc_canopy_grow_tempf(TKGO)

  D100: DO L=NU,NI(NZ)
    TKSO=TKS(L)+OFFST(NZ)
    fTgrowRootP(L,NZ)=calc_root_grow_tempf(TKSO)
  ENDDO D100
  PSICanPDailyMin(NZ)=AMIN1(PSICanPDailyMin(NZ),PSICanP(NZ))
  !
  !     DIURNAL CHILLING
  !
  !     TCelciusChill4Seed=chilling temperature from PFT file
  !     CHILL=accumulated chilling hours used to limit CO2 fixn in stomate.f
  !
  IF(TCelciusCanopy(NZ).LT.TCelciusChill4Seed(NZ))THEN
    CHILL(NZ)=AMIN1(24.0_r8,CHILL(NZ)+1.0_r8)
  ELSE
    CHILL(NZ)=AZMAX1(CHILL(NZ)-1.0_r8)
  ENDIF
  end associate
  end subroutine SetCanopyGrowthFuncs

end module UptakesMod
