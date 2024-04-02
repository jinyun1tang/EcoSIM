module StartsMod
!!
! Description:
! code to initalize soil variables

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : padr, print_info,check_bool
  use minimathMod, only : isclose, AZMAX1,AZMIN1
  use EcosimConst
  use TracerIDMod
  use SnowDataType
  use MicrobialDataType
  use EcoSIMSolverPar
  use SOMDataType
  use EcoSIMCtrlMod
  use ChemTranspDataType
  use FertilizerDataType
  use SnowPhysMod, only : InitSnowLayers
  use InitSOMBGCMod, only : InitSOMConsts,InitSOMProfile,InitSOMVars
  use CanopyRadDataType
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use LandSurfDataType
  use InitVegBGC, only : InitIrradianceGeometry
  use PlantTraitDataType
  use PlantDataRateType
  use SurfLitterDataType
  use CanopyDataType
  use SurfSoilDataType
  use SoilBGCDataType
  use PlantMngmtDataType
  use EcoSimSumDataType
  use RootDataType
  use EcosimBGCFluxType
  use SoilPropertyDataType
  use IrrigationDataType
  use SedimentDataType
  use GridDataType
  use MiniFuncMod
  use UnitMod, only : units
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  !
  !
  !     BulkDensLitR=dry bulk density of woody(0),fine(1),manure(2) litter
  !     FORGC=minimum SOC for organic soil (g Mg-1)
  !      FVLWB,FCH4F=maximum SWC,CH4 emission fraction for combustion
  !     PSIHY=hygroscopic water potential (MPa)
  !     FCI,WPI=FC,WP for water retention by ice (MPa)
  !     CDPTHSI=depth to bottom of snowpack layers
  !     POROQ=Penman Water Linear Reduction tortuosity used in gas flux calculations
  !

  public :: starts
  public :: set_ecosim_solver
  public :: startsim
  contains

  SUBROUTINE starts(NHW,NHE,NVN,NVS)
  !
  !     THIS SUBROUTINE INITIALIZES ALL SOIL VARIABLES

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX,L,NGL
  REAL(R8) :: ALTY
  real(r8) :: ALTZG
  real(r8) :: tPBOT
  real(r8) :: LandScape1stSoiLayDepth
  real(r8) :: YSIN(NumOfSkyAzimuSects),YCOS(NumOfSkyAzimuSects)
  real(r8) :: SkyAzimuthAngle(NumOfSkyAzimuSects)
! begin_execution


  !  Initialize controlling parameters
  call InitControlParms
  !
  !  IRRADIANCE INTERCEPTION GEOMETRY, plant model
  call InitIrradianceGeometry(YSIN,YCOS,SkyAzimuthAngle)
  !  CALCULATE ELEVATION OF EACH GRID CELL
  !
  call InitGridElevation(NHW,NHE,NVN,NVS,YSIN,YCOS,SkyAzimuthAngle,ALTY)
  !
  !     INITIALIZE ACCUMULATORS AND MASS BALANCE CHECKS
  !     OF EACH GRID CELL
  !
  call InitAccumulators()

! this assumes the whole landscape is a grid
  ALTZG=0.0_r8
  LandScape1stSoiLayDepth=0.0_r8

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      !
      !     MINIMUM SURFACE ELEVATION IN LANDSCAPE
      !
      !     ALT=surface elevation relative to maximum
      !     ALTZG=minimum surface elevation in landscape
      !CDPTH: depth to bottom of soil layer
      ALT(NY,NX)=ALT(NY,NX)-ALTY
      IF(NX.EQ.NHW.AND.NY.EQ.NVN)THEN
        ALTZG=ALT(NY,NX)
      ELSE
        ALTZG=MIN(ALTZG,ALT(NY,NX))
      ENDIF
      !
      LandScape1stSoiLayDepth=AMAX1(LandScape1stSoiLayDepth,CumDepth2LayerBottom(NU(NY,NX),NY,NX)) !topsoil layer depth
!
!     INITIALIZE ATMOSPHERE VARIABLES
!
!     C*E=atmospheric concentration (g m-3)
!     *E=atmospheric concentration from readi.f (umol mol-1)
!     CO2=CO2,CH4=CH4,OXY=O2,Z2G=N2,Z2O=N2O,NH3=NH3,H2G=H2
!
      tPBOT=PBOT(NY,NX)/1.01325E+02_r8
      CCO2EI(NY,NX)=CO2EI(NY,NX)*5.36E-04_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_CO2,NY,NX)=CO2E(NY,NX)*5.36E-04_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_CH4,NY,NX)=CH4E(NY,NX)*5.36E-04_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_O2,NY,NX)=OXYE(NY,NX)*1.43E-03_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_N2,NY,NX)=Z2GE(NY,NX)*1.25E-03_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_N2O,NY,NX)=Z2OE(NY,NX)*1.25E-03_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_NH3,NY,NX)=ZNH3E(NY,NX)*6.25E-04_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_H2,NY,NX)=H2GE(NY,NX)*8.92E-05_r8*Tref/TairKClimMean(NY,NX)*tPBOT
!
!     MICROBIAL THERMAL ADAPTATION
!
!     OFFSET=shift in Arrhenius curve used in nitro.f (oC)
!     ATCS=mean annual soil temperature (OC)
!
      OFFSET(NY,NX)=fOFFSET(ATCS(NY,NX))

!
!     INITIALIZE WATER POTENTIAL VARIABLES FOR SOIL LAYERS
!
!     LOGPSIMX,LOGPSIMN,LOGPSIAtSat=log water potential at FC,WP,POROS
!     LOGPSIMXD,LOGPSIMND=LOGPSIMX-LOGPSIAtSat,LOGPSIMN-LOGPSIMX
!
      LOGPSIAtSat(NY,NX)=LOG(-PSIPS)
      LOGPSIFLD(NY,NX)=LOG(-PSIAtFldCapacity(NY,NX))
      LOGPSIMN(NY,NX)=LOG(-PSIAtWiltPoint(NY,NX))
      LOGPSIMXD(NY,NX)=LOGPSIFLD(NY,NX)-LOGPSIAtSat(NY,NX)
      LOGPSIMND(NY,NX)=LOGPSIMN(NY,NX)-LOGPSIFLD(NY,NX)
!     SoilMicPMassLayer(0,NY,NX)=0.0_r8
!
!     DISTRIBUTION OF OM AMONG FRACTIONS OF DIFFERING
!     BIOLOGICAL ACTIVITY
!
      call InitHGrid(NY,NX)

      call InitLayerDepths(NY,NX)
    ! ActiveLayDepth=active layer depth (m)
      ActiveLayDepth(NY,NX)=9999.0_r8
!
!     INITIALIZE SNOWPACK LAYERS
      call InitSnowLayers(NY,NX)

!     VHeatCapLitR,VHCPNX=minimum heat capacities for solving
!      surface litter,soil layer water and heat fluxes
      VHeatCapLitR(NY,NX)=VLHeatCapLitRMin*AREA(3,NU(NY,NX),NY,NX)
      VHCPNX(NY,NX)=VLHeatCapSoiMin*AREA(3,NU(NY,NX),NY,NX)

!
!     SURFACE WATER STORAGE AND LOWER HEAT SINK

!      SoilHeatSrcDepth=depth at which soil heat sink-source calculated
!     TCNDG=assumed thermal conductivity below lower soil boundary
!     (MJ m-1 K-1 h-1)
!     TKSD=deep source/sink temperature from geothermal flux(K)

      SoilHeatSrcDepth(NY,NX)=AMAX1(10.0_r8,CumDepth2LayerBottom(NL(NY,NX),NY,NX)+1.0_r8)
      TCS(0,NY,NX)=ATCS(NY,NX)
      TKS(0,NY,NX)=ATKS(NY,NX)
      TKSD(NY,NX)=ATKS(NY,NX)+2.052E-04_r8*SoilHeatSrcDepth(NY,NX)/TCNDG
!

    ENDDO
  ENDDO

!
!     INITIALIZE COMMUNITY CANOPY
!
  CanopyHeight_col(:,:)=0.0_r8
  CanopyHeightZ_col(0:NumOfCanopyLayers,:,:)=0.0_r8
  CanopyLeafAareZ_col(1:NumOfCanopyLayers,:,:)=0.0_r8
  CanopyStemAareZ_col(1:NumOfCanopyLayers,:,:)=0.0_r8
  WGLFT(1:NumOfCanopyLayers,:,:)=0.0_r8

!
  call InitSoilVars(NHW,NHE,NVN,NVS,ALTZG,LandScape1stSoiLayDepth)

  RETURN
  END subroutine starts
!------------------------------------------------------------------------------------------
  subroutine InitSoilVars(NHW,NHE,NVN,NVS,ALTZG,LandScape1stSoiLayDepth)
  !     N3,N2,N1=L,NY,NX of source grid cell
  !     N6,N5,N4=L,NY,NX of destination grid cell
  !      ALTZG=minimum surface elevation in landscape
  !     DTBLI,DTBLDI=depth of natural,artificial water table
  !     WaterTBLSlope=slope of natural water table relative to landscape surface
  !     in geography, slope =rise/run
  !     ExtWaterTablet0,DTBLD=depth of natural,artificial water table adjusted for elevn
  !     DepthInternalWTBL=depth to internal water table
  !     DIST=distance between adjacent layers:1=EW,2=NS,3=vertical(m)
  !     XDPTH=x-section area/distance in solute flux calculations (m2/m)
  !     DISP=dispersivity parameter in solute flux calculations (m2 h-1)
  !
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: ALTZG,LandScape1stSoiLayDepth
  integer :: NY,NX,L,N
  integer :: N1,N2,N3,N4,N5,N6


  ! begin_execution

!     INITIALIZE SEDIMENT LOAD IN EROSION MODEL
!
  IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
    SED(:,:)=0.0_r8
  ENDIF
!
!     INITIALIZE GRID CELL DIMENSIONS
! altz: topographic altitude
! ALT: grid altitude
! DTBLI: external water table depth, before applying the altitude correction
! WaterTBLSlope: slope of water table relative to surface slope
! ExtWaterTablet0: external water table depth
! DTBLDI: depth of artificial water table
! DepthInternalWTBL: internal water table depth
! DTBLD: artifical water table depth, before applying the altitude correction
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      ALTZ(NY,NX)=ALTZG
      IF(SoiBulkDensity(NU(NY,NX),NY,NX).GT.0.0_r8)THEN
        ExtWaterTablet0(NY,NX)=DTBLI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))*(1.0_r8-WaterTBLSlope(NY,NX))
        DTBLD(NY,NX)=AZMAX1(DTBLDI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))*(1.0_r8-WaterTBLSlope(NY,NX)))
      ELSE
        ExtWaterTablet0(NY,NX)=0.0_r8
        DTBLD(NY,NX)=0.0_r8
      ENDIF
      DepthInternalWTBL(NY,NX)=ExtWaterTablet0(NY,NX)
    ENDDO
  ENDDO

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      DO L=1,NL(NY,NX)
        N1=NX
        N2=NY
        N3=L
        DO N=FlowDirIndicator(N2,N1),3
          IF(N.EQ.1)THEN
            ! in direction x, west-east
            IF(NX.EQ.NHE)THEN
              cycle
            ELSE
              N4=NX+1
              N5=NY
              N6=L
            ENDIF
          ELSEIF(N.EQ.2)THEN
            ! in direction y, north-south
            IF(NY.EQ.NVS)THEN
              cycle
            ELSE
              N4=NX
              N5=NY+1
              N6=L
            ENDIF
          ELSEIF(N.EQ.3)THEN
            ! in vertical, up-down
            IF(L.EQ.NL(NY,NX))THEN
              cycle
            ELSE
              N4=NX
              N5=NY
              N6=L+1
            ENDIF
          ENDIF
          DIST(N,N6,N5,N4)=0.5_r8*(DLYR(N,N3,N2,N1)+DLYR(N,N6,N5,N4))
          XDPTH(N,N6,N5,N4)=AREA(N,N3,N2,N1)/DIST(N,N6,N5,N4)
          !1.07 is a scaling parameter for dispersion calculation, reference?
          DISP(N,N6,N5,N4)=0.20_r8*DIST(N,N6,N5,N4)**1.07_r8
        ENDDO

        IF(L.EQ.NU(NY,NX))THEN
          DIST(3,N3,N2,N1)=0.5_r8*DLYR(3,N3,N2,N1)
          XDPTH(3,N3,N2,N1)=AREA(3,N3,N2,N1)/DIST(3,N3,N2,N1)
          DISP(3,N3,N2,N1)=0.20_r8*DIST(3,N3,N2,N1)**1.07_r8
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      !
      !     ALLOCATE LITTER,SOC TO WOODY,NON-WOODY,MANURE,POC AND HUMUS
      !
      call InitSoilProfile(NY,NX,LandScape1stSoiLayDepth)
      !
      !     SURFACE LITTER HEAT CAPACITY
      !
      SoilMicPMassLayerMn(NY,NX)=AZMAX1(SAND(NU(NY,NX),NY,NX)+SILT(NU(NY,NX),NY,NX)+CLAY(NU(NY,NX),NY,NX))
      VHeatCapacity(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VLWatMicP(0,NY,NX)+cpi*VLiceMicP(0,NY,NX)
      VHeatCapacitySoilM(0,NY,NX)=0.0_r8
      VLMicPt0(0,NY,NX)=0.0_r8
    ENDDO
  ENDDO
  end subroutine InitSoilVars
!------------------------------------------------------------------------------------------
  subroutine InitSoilProfile(NY,NX,LandScape1stSoiLayDepth)

  implicit none
  integer, intent(in) :: NY,NX
  REAL(R8),INTENT(IN) :: LandScape1stSoiLayDepth
  integer  :: L,M,K,N,KK,NN,NGL
  real(r8) :: CORGCM,HCX,TORGC
  real(r8) :: CORGL,TORGLL,FCX
  REAL(R8) :: PTDS
  real(r8) :: TORGM
  real(r8) :: VSAND
  real(r8) :: TORGL(JZ)
  real(r8) :: VMINL
  real(r8) :: VORGC
! begin_execution
! RSC,RSC,RSP=C,N,P in fine(1),woody(0),manure(2) litter (g m-2)
! CORGC,CORGR,CORGN,CORGP=SOC,POC,SON,SOP (g Mg-1)
! BKVL=soil mass (Mg)
!

!
!     INITIALIZE SOM FROM ORGANIC INPUTS IN SOIL FILE FROM 'READS'
!
!     CORGC,CORGR,CORGN,CORGP=SOC,POC,SON,SOP (g Mg-1)
!
  TORGC=0.0_r8
  D1190: DO L=NU(NY,NX),NL(NY,NX)
    !     CORGCZ=CORGC(L,NY,NX)
    !     CORGRZ=CORGR(L,NY,NX)
    !     CORGNZ=CORGN(L,NY,NX)
    !     CORGPZ=CORGP(L,NY,NX)
    !
    !     ALLOCATE SOC TO POC(3) AND HUMUS(4)
    !
    !     CORGCX(3)=CORGRZ
    !     CORGCX(4)=AZMAX1(CORGCZ-CORGCX(3))
    !     CORGNX(3)=AMIN1(CNRH(3)*CORGCX(3),CORGNZ)
    !     CORGNX(4)=AZMAX1(CORGNZ-CORGNX(3))
    !     CORGPX(3)=AMIN1(CPRH(3)*CORGCX(3),CORGPZ)
    !     CORGPX(4)=AZMAX1(CORGPZ-CORGPX(3))
    CORGL=AZMAX1(CORGC(L,NY,NX)-CORGR(L,NY,NX))
    TORGL(L)=TORGC+CORGL*SoilMicPMassLayer(L,NY,NX)/AREA(3,L,NY,NX)*0.5_r8
    TORGC=TORGC+CORGL*SoilMicPMassLayer(L,NY,NX)/AREA(3,L,NY,NX)
  ENDDO D1190
!
!     PARAMETERS TO ALLOCATE HUMUS TO LESS OR MORE RECALCITRANT FRACTIONS
!
!     TORGL=accumulated humus down to soil layer (g m-2)
!     TORGM=TORGL used to calculate allocation (g m-2)
!     HCX=shape parameter for depth effect on allocation
!
  TORGM=AMAX1(2.0E+03_r8,AMIN1(5.0E+03_r8,0.25_r8*TORGL(MaxNumRootLays(NY,NX))))
  IF(TORGM.GT.ZERO)THEN
    HCX=LOG(0.5_r8)/TORGM
  ELSE
    HCX=0.0_r8
  ENDIF

  D1200: DO L=0,NL(NY,NX)
    !
    if(L==0)then
      TORGLL=0.0_r8
    else
      TORGLL=TORGL(L)
    endif

    call InitSOMProfile(L,NY,NX,HCX,LandScape1stSoiLayDepth,TORGLL,CORGCM,FCX)
    !
    !     LAYER WATER, ICE, AIR CONTENTS
    !
    !     PSISE,PSISA=water potential at saturation,air entry (MPa)
    !     PTDS=particle density (Mg m-3)
    !     POROS=total porosity
    !     VOLA,VOLAH=micropore,macropore volume
    !     THW,THI=initial soil water,ice content
    !     VOLW,VOLWH=micropore,macropore water volume(m3)
    !     VOLI,VOLIH=micropore,macropore ice volume(m3)
    !     VOLP=total air volume (m3)
    !
    PSISE(L,NY,NX)=PSIPS
    PSISoilAirEntry(L,NY,NX)=-1.5E-03_r8
    ROXYF(L,NY,NX)=0.0_r8
    RCO2F(L,NY,NX)=0.0_r8
    ROXYL(L,NY,NX)=0.0_r8
    RCH4F(L,NY,NX)=0.0_r8
    RCH4L(L,NY,NX)=0.0_r8

    IF(L.GT.0)THEN
      IF(SoiBulkDensity(L,NY,NX).GT.ZERO)THEN
        !it is a soil layer
        !compute particle density
        PTDS=ppmc*(1.30_r8*CORGCM+2.66_r8*(1.0E+06_r8-CORGCM))
        POROS(L,NY,NX)=1.0_r8-(SoiBulkDensity(L,NY,NX)/PTDS)
      ELSE
        !for ponding water
        PTDS=0.0_r8
        POROS(L,NY,NX)=1.0_r8
      ENDIF
      POROSI(L,NY,NX)=POROS(L,NY,NX)*FracSoiAsMicP(L,NY,NX)
      VLMicP(L,NY,NX)=POROS(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
      VLMicPt0(L,NY,NX)=VLMicP(L,NY,NX)
      VLMacP(L,NY,NX)=SoilFracAsMacP(L,NY,NX)*VGeomLayert0(L,NY,NX)
      !
      !     LAYER HEAT CONTENTS
      !
      !     SAND,SILT,CLAY=sand,silt,clay mass (Mg)
      !     VORGC,VMINL,VSAND=volumetric fractions of SOC,non-sand,sand
      !     VHCM,VHCP=volumetric dry,wet soil heat capacity (MJ m-3 K-1)
      !     TKS,TCS=soil temperature (oC,K)
      !     THETW,THETI,THETP=micropore water,ice,air concentration (m3 m-3)
!
      SAND(L,NY,NX)=CSAND(L,NY,NX)*SoilMicPMassLayer(L,NY,NX)
      SILT(L,NY,NX)=CSILT(L,NY,NX)*SoilMicPMassLayer(L,NY,NX)
      CLAY(L,NY,NX)=CCLAY(L,NY,NX)*SoilMicPMassLayer(L,NY,NX)
      IF(SoiBulkDensity(L,NY,NX).GT.ZERO)THEN
        ! PTDS=particle density (Mg m-3)
        ! soil volumetric heat capacity
        VORGC=CORGCM*ppmc*SoiBulkDensity(L,NY,NX)/PTDS
        VMINL=(CSILT(L,NY,NX)+CCLAY(L,NY,NX))*SoiBulkDensity(L,NY,NX)/PTDS
        VSAND=CSAND(L,NY,NX)*SoiBulkDensity(L,NY,NX)/PTDS
        VHeatCapacitySoilM(L,NY,NX)=((2.496_r8*VORGC+2.385_r8*VMINL+2.128_r8*VSAND) &
          *FracSoiAsMicP(L,NY,NX)+2.128_r8*ROCK(L,NY,NX))*VGeomLayer(L,NY,NX)
      ELSE
        VHeatCapacitySoilM(L,NY,NX)=0.0_r8
      ENDIF
!
      !     INITIAL SOIL WATER AND ICE CONTENTS
!
      IF(ISOIL(isoi_fc,L,NY,NX).EQ.0.AND.ISOIL(isoi_wp,L,NY,NX).EQ.0)THEN
      ! field capacity and wilting point are read from input
        IF(THW(L,NY,NX).GT.1.0_r8)THEN
          THETW(L,NY,NX)=POROS(L,NY,NX)
        ELSEIF(isclose(THW(L,NY,NX),1.0_r8))THEN
          THETW(L,NY,NX)=FieldCapacity(L,NY,NX)
        ELSEIF(isclose(THW(L,NY,NX),0.0_r8))THEN
          THETW(L,NY,NX)=WiltPoint(L,NY,NX)
        ELSEIF(THW(L,NY,NX).LT.0.0)THEN
          THETW(L,NY,NX)=0.0_r8
        ELSE
          THETW(L,NY,NX)=THW(L,NY,NX)
        ENDIF
        IF(THI(L,NY,NX).GT.1.0_r8)THEN
          THETI(L,NY,NX)=AZMAX1(AMIN1(POROS(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
        ELSEIF(isclose(THI(L,NY,NX),1.0_r8))THEN
          THETI(L,NY,NX)=AZMAX1(AMIN1(FieldCapacity(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
        ELSEIF(isclose(THI(L,NY,NX),0.0_r8))THEN
          THETI(L,NY,NX)=AZMAX1(AMIN1(WiltPoint(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
        ELSEIF(THI(L,NY,NX).LT.0.0_r8)THEN
          THETI(L,NY,NX)=0.0_r8
        ELSE
          THETI(L,NY,NX)=THI(L,NY,NX)
        ENDIF
        VLWatMicP(L,NY,NX)=THETW(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
        VLWatMicPX(L,NY,NX)=VLWatMicP(L,NY,NX)
        VLWatMacP(L,NY,NX)=THETW(L,NY,NX)*VLMacP(L,NY,NX)
        VLiceMicP(L,NY,NX)=THETI(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
        VLiceMacP(L,NY,NX)=THETI(L,NY,NX)*VLMacP(L,NY,NX)
!       total air-filled porosity, micropores + macropores
        VLsoiAirP(L,NY,NX)=AZMAX1(VLMicP(L,NY,NX)-VLWatMicP(L,NY,NX)-VLiceMicP(L,NY,NX)) &
          +AZMAX1(VLMacP(L,NY,NX)-VLWatMacP(L,NY,NX)-VLiceMacP(L,NY,NX))
        VHeatCapacity(L,NY,NX)=VHeatCapacitySoilM(L,NY,NX)+cpw*(VLWatMicP(L,NY,NX) &
          +VLWatMacP(L,NY,NX))+cpi*(VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX))
        THETWZ(L,NY,NX)=THETW(L,NY,NX)
        THETIZ(L,NY,NX)=THETI(L,NY,NX)
      ENDIF
    ENDIF
    TKS(L,NY,NX)=ATKS(NY,NX)
    TCS(L,NY,NX)=ATCS(NY,NX)
    !
    !     INITIALIZE SOM VARIABLES
    call InitSOMVars(L,NY,NX,FCX)
    !
  ENDDO D1200
  POROSI(0,NY,NX)=1._r8  !this is added for numerical fixing
  !
  !  INITIALIZE FERTILIZER ARRAYS
  call initFertArrays(NY,NX)

  end subroutine InitSoilProfile
!------------------------------------------------------------------------------------------

  subroutine initFertArrays(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L2

! begin_execution
  L2=NL(NY,NX)
  FertN_soil(ifertn_beg:ifertn_end,0:L2,NY,NX)=0._r8
  FertN_band(ifertnb_beg:ifertnb_end,1:L2,NY,NX)=0._r8
  trcs_VLN_vr(ids_NH4,0:L2,NY,NX)=1.0_r8
  trcs_VLN_vr(idg_NH3,0:L2,NY,NX)=trcs_VLN_vr(ids_NH4,0:L2,NY,NX)
  trcs_VLN_vr(ids_NO3,0:L2,NY,NX)=1.0_r8
  trcs_VLN_vr(ids_NO2,0:L2,NY,NX)=trcs_VLN_vr(ids_NO3,0:L2,NY,NX)
  trcs_VLN_vr(ids_H1PO4,0:L2,NY,NX)=1.0_r8
  trcs_VLN_vr(ids_H2PO4,0:L2,NY,NX)=trcs_VLN_vr(ids_H1PO4,0:L2,NY,NX)
  trcs_VLN_vr(ids_NH4B,0:L2,NY,NX)=0.0_r8
  trcs_VLN_vr(idg_NH3B,0:L2,NY,NX)= trcs_VLN_vr(ids_NH4B,0:L2,NY,NX)
  trcs_VLN_vr(ids_NO3B,0:L2,NY,NX)=0.0_r8
  trcs_VLN_vr(ids_NO2B,0:L2,NY,NX)=trcs_VLN_vr(ids_NO3B,0:L2,NY,NX)
  trcs_VLN_vr(ids_H1PO4B,0:L2,NY,NX)=0.0_r8
  trcs_VLN_vr(ids_H2PO4B,0:L2,NY,NX)=trcs_VLN_vr(ids_H1PO4B,0:L2,NY,NX)

  ROXYX(0:L2,NY,NX)=0.0_r8
  RNH4X(0:L2,NY,NX)=0.0_r8
  RNO3X(0:L2,NY,NX)=0.0_r8
  RNO2X(0:L2,NY,NX)=0.0_r8
  RN2OX(0:L2,NY,NX)=0.0_r8
  RPO4X(0:L2,NY,NX)=0.0_r8
  RP14X(0:L2,NY,NX)=0.0_r8
  RVMXC(0:L2,NY,NX)=0.0_r8
  RNHBX(0:L2,NY,NX)=0.0_r8
  RN3BX(0:L2,NY,NX)=0.0_r8
  RN2BX(0:L2,NY,NX)=0.0_r8
  RPOBX(0:L2,NY,NX)=0.0_r8
  RP1BX(0:L2,NY,NX)=0.0_r8
  RVMBC(0:L2,NY,NX)=0.0_r8
  ZNHUI(0:L2,NY,NX)=0.0_r8
  ZNHU0(0:L2,NY,NX)=0.0_r8
  ZNFNI(0:L2,NY,NX)=0.0_r8
  ZNFN0(0:L2,NY,NX)=0.0_r8

  WDNHB(1:L2,NY,NX)=0.0_r8
  DPNHB(1:L2,NY,NX)=0.0_r8
  WDNOB(1:L2,NY,NX)=0.0_r8
  DPNOB(1:L2,NY,NX)=0.0_r8
  WDPOB(1:L2,NY,NX)=0.0_r8
  DPPOB(1:L2,NY,NX)=0.0_r8
  COCU(1:jcplx,1:L2,NY,NX)=0.0_r8
  CONU(1:jcplx,1:L2,NY,NX)=0.0_r8
  COPU(1:jcplx,1:L2,NY,NX)=0.0_r8
  COAU(1:jcplx,1:L2,NY,NX)=0.0_r8

  end subroutine initFertArrays

!------------------------------------------------------------------------------------------
  subroutine InitGridElevation(NHW,NHE,NVN,NVS,YSIN,YCOS,SkyAzimuthAngle,ALTY)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: YSIN(NumOfSkyAzimuSects),YCOS(NumOfSkyAzimuSects),SkyAzimuthAngle(NumOfSkyAzimuSects)
  REAL(R8),INTENT(OUT):: ALTY
  integer :: NY,NX,N,NN
  REAL(R8) :: DGAZI
  real(r8) :: SineGrndSurfAzimuth_col(JY,JX),CosineGrndSurfAzimuth_col(JY,JX)  !diagnostic

! begin_execution
! GroundSurfAzimuth_col=ground surface azimuth, aspect in radians
! SineGrndSlope_col,CosineGrndSlope_col=sine,cosine of ground surface
! OMEGAG=incident sky angle at ground surface
! SLOPE=sine of ground surface slope in (0)aspect, (1)EW,(2)NS directions
! ALT=ground surface elevation
! ALTY=maximum surface elevation in landscape
! XGridRunoffFlag=runoff boundary flags:0=not possible,1=possible
! ASP=aspect angle in degree
  ALTY=0.0_r8
  write(*,1112)'NY','NX','east','west','south','north','altitude','Dist(m):E-W','Dist(m):N-S',&
    'aspect(o)','slope(o)','slope0','slope-east','slope-north','SineGrndSlope_col','CosineGrndSlope_col','SineGrndSurfAzimuth_col'
1112    FORMAT(2A4,4A6,20A12)
  D9985: DO NX=NHW,NHE
    D9980: DO NY=NVN,NVS
      ZEROS(NY,NX)=ZERO*DH(NY,NX)*DV(NY,NX)
      ZEROS2(NY,NX)=ZERO2*DH(NY,NX)*DV(NY,NX)
!     compute slopes
      GroundSurfAzimuth_col(NY,NX)=ASP(NY,NX)*RadianPerDegree   !radian
      SineGrndSurfAzimuth_col(NY,NX)=ABS(SIN(GroundSurfAzimuth_col(NY,NX)))
      CosineGrndSurfAzimuth_col(NY,NX)=ABS(COS(GroundSurfAzimuth_col(NY,NX)))
      SLOPE(0,NY,NX)=AMAX1(1.745E-04_r8,SIN(SL(NY,NX)*RadianPerDegree))  !small slope approximation

      IF(ASP(NY,NX).GE.0.0_r8.AND.ASP(NY,NX).LT.90.0_r8)THEN
      ! along the northeast
        SLOPE(1,NY,NX)=-SLOPE(0,NY,NX)*COS(ASP(NY,NX)*RadianPerDegree)    !to south (thus -)
        SLOPE(2,NY,NX)=SLOPE(0,NY,NX)*SIN(ASP(NY,NX)*RadianPerDegree)     !to east
        XGridRunoffFlag(1,1,NY,NX)=.true.    !east
        XGridRunoffFlag(2,1,NY,NX)=.false.
        XGridRunoffFlag(1,2,NY,NX)=.false.
        XGridRunoffFlag(2,2,NY,NX)=.true.    !north
      ELSEIF(ASP(NY,NX).GE.90.0_r8.AND.ASP(NY,NX).LT.180.0_r8)THEN
      ! along the southeast
        SLOPE(1,NY,NX)=SLOPE(0,NY,NX)*SIN((ASP(NY,NX)-90.0_r8)*RadianPerDegree)   !to south
        SLOPE(2,NY,NX)=SLOPE(0,NY,NX)*COS((ASP(NY,NX)-90.0_r8)*RadianPerDegree)   !to east
        XGridRunoffFlag(1,1,NY,NX)=.false.
        XGridRunoffFlag(2,1,NY,NX)=.true.   !west to east
        XGridRunoffFlag(1,2,NY,NX)=.false.
        XGridRunoffFlag(2,2,NY,NX)=.true.   !north to south
      ELSEIF(ASP(NY,NX).GE.180.0_r8.AND.ASP(NY,NX).LT.270.0_r8)THEN
      !along the southwest
        SLOPE(1,NY,NX)=SLOPE(0,NY,NX)*COS((ASP(NY,NX)-180.0_r8)*RadianPerDegree)    !to south
        SLOPE(2,NY,NX)=-SLOPE(0,NY,NX)*SIN((ASP(NY,NX)-180.0_r8)*RadianPerDegree)   !to east (thus -)
        XGridRunoffFlag(1,1,NY,NX)=.false.
        XGridRunoffFlag(2,1,NY,NX)=.true.  !west
        XGridRunoffFlag(1,2,NY,NX)=.true.  !south
        XGridRunoffFlag(2,2,NY,NX)=.false.
      ELSEIF(ASP(NY,NX).GE.270.0_r8.AND.ASP(NY,NX).LE.360.0_r8)THEN
      ! along the northwest
        SLOPE(1,NY,NX)=-SLOPE(0,NY,NX)*SIN((ASP(NY,NX)-270.0_r8)*RadianPerDegree)   !to north (thus -)
        SLOPE(2,NY,NX)=-SLOPE(0,NY,NX)*COS((ASP(NY,NX)-270.0_r8)*RadianPerDegree)   !to west (thus -)
        XGridRunoffFlag(1,1,NY,NX)=.true.  !east
        XGridRunoffFlag(2,1,NY,NX)=.false.
        XGridRunoffFlag(1,2,NY,NX)=.true.  !south
        XGridRunoffFlag(2,2,NY,NX)=.false.
      ENDIF
      SLOPE(3,NY,NX)=-1.0_r8

      IF(.not.isclose(SLOPE(1,NY,NX),0.0_r8).OR.(.not.isclose(SLOPE(2,NY,NX),0.0_r8)))THEN
        FSLOPE(1,NY,NX)=ABS(SLOPE(1,NY,NX))/(ABS(SLOPE(1,NY,NX))+ABS(SLOPE(2,NY,NX)))  !
        FSLOPE(2,NY,NX)=ABS(SLOPE(2,NY,NX))/(ABS(SLOPE(1,NY,NX))+ABS(SLOPE(2,NY,NX)))
      ELSE
        FSLOPE(1,NY,NX)=0.5_r8
        FSLOPE(2,NY,NX)=0.5_r8
      ENDIF

!    compute incident sky angle at ground surface
      SineGrndSlope_col(NY,NX)=SLOPE(0,NY,NX)    !this is exact
      CosineGrndSlope_col(NY,NX)=SQRT(1.0_r8-SineGrndSlope_col(NY,NX)**2._r8)
      D240: DO N=1,NumOfSkyAzimuSects
        DGAZI=COS(GroundSurfAzimuth_col(NY,NX)-SkyAzimuthAngle(N))
        OMEGAG(N,NY,NX)=AZMAX1(AMIN1(1.0_r8,CosineGrndSlope_col(NY,NX)*YSIN(N)+ &
          SineGrndSlope_col(NY,NX)*YCOS(N)*DGAZI))
      ENDDO D240
!     compute ground surface elevation
!     DH, length in e-w direction
      IF(NX.EQ.NHW)THEN
        IF(NY.EQ.NVN)THEN
          !(west, north) corner
          ALT(NY,NX)=0.5_r8*DH(NY,NX)*SLOPE(1,NY,NX)+0.5_r8*DV(NY,NX)*SLOPE(2,NY,NX)
        ELSE
          !west boundary
          ALT(NY,NX)=ALT(NY-1,NX) &
            +1.0_r8*DH(NY,NX)*SLOPE(1,NY,NX) &
            +0.5_r8*DV(NY-1,NX)*(SLOPE(2,NY-1,NX)) &
            +0.5_r8*DV(NY,NX)*SLOPE(2,NY,NX)
        ENDIF
      ELSE
        IF(NY.EQ.NVN)THEN
          !north boundary
          ALT(NY,NX)=ALT(NY,NX-1) &
            +0.5_r8*DH(NY,NX-1)*SLOPE(1,NY,NX-1) &
            +0.5_r8*DH(NY,NX)*SLOPE(1,NY,NX) &
            +0.5_r8*DV(NY,NX-1)*SLOPE(2,NY,NX-1) &
            +0.5_r8*DV(NY,NX)*SLOPE(2,NY,NX)
        ELSE
          ALT(NY,NX)=(ALT(NY,NX-1) &
            +0.5_r8*DH(NY,NX-1)*SLOPE(1,NY,NX-1) &
            +0.5_r8*DH(NY,NX)*SLOPE(1,NY,NX) &
            +ALT(NY-1,NX) &
            +0.5_r8*DV(NY-1,NX)*SLOPE(2,NY-1,NX) &
            +0.5_r8*DV(NY,N)*SLOPE(2,NY,NX))/2.0
        ENDIF
      ENDIF

      IF(NX.EQ.NHW.AND.NY.EQ.NVN)THEN
        !(west,north) corner
        ALTY=ALT(NY,NX)
      ELSE
        ALTY=MAX(ALTY,ALT(NY,NX))
      ENDIF
      WRITE(*,1111)NX,NY,((XGridRunoffFlag(NN,N,NY,NX),NN=1,2),N=1,2) &
        ,ALT(NY,NX),DH(NY,NX),DV(NY,NX),ASP(NY,NX),SL(NY,NX) &
        ,SLOPE(0,NY,NX),SLOPE(1,NY,NX),SLOPE(2,NY,NX) &
        ,SineGrndSlope_col(NY,NX),CosineGrndSurfAzimuth_col(NY,NX),SineGrndSurfAzimuth_col(NY,NX)
1111  FORMAT(2I4,4L6,20E12.4)
    ENDDO D9980
  ENDDO D9985
  end subroutine InitGridElevation
!------------------------------------------------------------------------------------------
  subroutine InitControlParms
  implicit none
  !     begin_execution
  !
  !     NPH=no. of cycles per hour for water, heat and solute flux calculns
  !     NPT=number of cycles per water iteration for gas flux calculations
  !     NPG=number of cycles per hour for gas flux calculations
  !     NPR,NPS=number of cycles NPH-1 for litter, snowpack flux calculns
  !     THETX=minimum air-filled porosity for gas flux calculations
  !     THETPI,DENSICE=ice porosity,density
  !     BulkDensLitR=surface litter bulk density, Mg m-3

  BulkDensLitR=(/0.0333_r8,0.0167_r8,0.0167_r8/)

  call InitSOMConsts
  !     NDIM=1
  !     IF(NHE.GT.NHW)NDIM=NDIM+1
  !     IF(NVS.GT.NVN)NDIM=NDIM+1
  !     XDIM=1.0/NDIM
  ZERO=1.0E-15_r8
  ZERO2=1.0E-08_r8
  TAREA=0.0_r8  !land scape area
  !
  !     INITIALIZE MASS BALANCE CHECKS
  !
  CRAIN=0.0_r8
  HEATIN=0.0_r8
  CO2GIN=0.0_r8
  OXYGIN=0.0_r8
  H2GIN=0.0_r8
  TZIN=0.0_r8
  ZN2GIN=0.0_r8
  TPIN=0.0_r8
  TORGF=0.0_r8
  TORGN=0.0_r8
  TORGP=0.0_r8
  VOLWOU=0.0_r8
  CEVAP=0.0_r8
  CRUN=0.0_r8
  HEATOU=0.0_r8
  OXYGOU=0.0_r8
  H2GOU=0.0_r8
  TSEDOU=0.0_r8
  TOMOU(ielmc)=0.0_r8
  TOMOU(ielmn)=0.0_r8
  TOMOU(ielmp)=0.0_r8
  XESN(ielmc)=0.0_r8
  XESN(ielmn)=0.0_r8
  XESN(ielmp)=0.0_r8
  TIONIN=0.0_r8
  TIONOU=0.0_r8
  end subroutine InitControlParms
!------------------------------------------------------------------------------------------
  subroutine InitAccumulators()
  implicit none

!     begin_execution
  TDTPX(:,:,:)=0.0_r8
  TDTPN(:,:,:)=0.0_r8
  TDRAD(:,:,:)=1.0_r8
  TDWND(:,:,:)=1.0_r8
  TDHUM(:,:,:)=1.0_r8
  TDPRC(:,:,:)=1.0_r8
  TDIRI(:,:,:)=1.0_r8
  TDCN4(:,:,:)=1.0_r8
  TDCNO(:,:,:)=1.0_r8

  IUTYP(:,:)=0
  IFNHB(:,:)=0
  IFNOB(:,:)=0
  IFPOB(:,:)=0
  IFLGS(:,:)=1
  NumActivePlants(:,:)=0
  ATCA(:,:)=ATCAI(:,:)
  ATCS(:,:)=ATCAI(:,:)
  TairKClimMean(:,:)=units%Celcius2Kelvin(ATCA)
  ATKS(:,:)=units%Celcius2Kelvin(ATCS)
  URAIN(:,:)=0.0_r8

  CO2byFire_col(:,:)=0.0_r8
  CH4byFire_col(:,:)=0.0_r8
  O2byFire_col(:,:)=0.0_r8
  N2ObyFire_col(:,:)=0.0_r8
  NH3byFire_col(:,:)=0.0_r8
  PO4byFire_col(:,:)=0.0_r8
  AmendCFlx_col(:,:)=0.0_r8
  FertNFlx_col(:,:)=0.0_r8
  FerPFlx_col(:,:)=0.0_r8
  UVOLO(:,:)=0.0_r8
  UEVAP(:,:)=0.0_r8
  URUN(:,:)=0.0_r8
  USEDOU(:,:)=0.0_r8
  UCOP(:,:)=0.0_r8
  HydroSufDOCFlx_col(:,:)=0.0_r8
  HydroSubsDOCFlx_col(:,:)=0.0_r8
  HydroSufDONFlx_col(:,:)=0.0_r8
  HydroSubsDONFlx_col(:,:)=0.0_r8
  HydroSufDOPFlx_col(:,:)=0.0_r8
  HydroSubsDOPFlx_col(:,:)=0.0_r8
  HydroSufDICFlx_col(:,:)=0.0_r8
  HydroSubsDICFlx_col(:,:)=0.0_r8
  HydroSufDINFlx_col(:,:)=0.0_r8
  HydroSubsDINFlx_col(:,:)=0.0_r8
  HydroSufDIPFlx_col(:,:)=0.0_r8
  HydroSubsDIPFlx_col(:,:)=0.0_r8
  HydroIonFlx_col(:,:)=0.0_r8
  LiterfalOrgM_col(ielmc,:,:)=0.0_r8
  LiterfalOrgM_col(ielmn,:,:)=0.0_r8
  LiterfalOrgM_col(ielmp,:,:)=0.0_r8
  UDRAIN(:,:)=0.0_r8
  ZDRAIN(:,:)=0.0_r8
  PDRAIN(:,:)=0.0_r8
  DPNH4(:,:)=0.0_r8
  DPNO3(:,:)=0.0_r8
  DPPO4(:,:)=0.0_r8
  trc_solml_vr(idg_O2,0,:,:)=0.0_r8
  FracSWRad2Grnd(:,:)=1.0_r8
  LWRadBySurf(:,:)=0.0_r8
  LWRadCanG(:,:)=0.0_r8
  Eco_NetRad_col(:,:)=0.0_r8
  Eco_Heat_Latent_col(:,:)=0.0_r8
  Eco_Heat_Sens_col(:,:)=0.0_r8
  Eco_Heat_Grnd_col(:,:)=0.0_r8
  Canopy_Heat_Latent_col(:,:)=0.0_r8
  Canopy_Heat_Sens_col(:,:)=0.0_r8
  TLEX(:,:)=0.0_r8
  TSHX(:,:)=0.0_r8
  Eco_NEE_col(:,:)=0.0_r8
  CanH2OHeldVg(:,:)=0.0_r8
  CanopyLeafArea_col(:,:)=0.0_r8
  StemArea_col(:,:)=0.0_r8
  PrecIntcptByCanG(:,:)=0.0_r8
  PPT(:,:)=0.0_r8
  DayLenthCurrent(:,:)=12.0_r8
  SurfAlbedo_col(:,:)=SoilAlbedo(:,:)
  EcoHavstElmnt_col(:,:,:)=0.0_r8
  EnergyImpact4Erosion(:,:)=0.0_r8
  end subroutine InitAccumulators

!------------------------------------------------------------------------------------------
  subroutine InitHGrid(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  DO  L=0,NL(NY,NX)
    DLYRI(1,L,NY,NX)=DH(NY,NX)    !east-west direction
    DLYRI(2,L,NY,NX)=DV(NY,NX)    !north-south direction
    DLYR(1,L,NY,NX)=DLYRI(1,L,NY,NX)
    DLYR(2,L,NY,NX)=DLYRI(2,L,NY,NX)
    AREA(3,L,NY,NX)=DLYR(1,L,NY,NX)*DLYR(2,L,NY,NX)
  ENDDO
  end subroutine InitHGrid
!------------------------------------------------------------------------------------------
  subroutine InitLayerDepths(NY,NX)
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer, intent(in) :: NY, NX
  integer :: L,K
  real(r8) :: VLitR0
  associate(                              &
    NumOfLitrCmplxs    => micpar%NumOfLitrCmplxs    , &
    k_woody_litr => micpar%k_woody_litr , &
    k_fine_litr  => micpar%k_fine_litr  , &
    k_manure     => micpar%k_manure       &
  )

!     begin_execution
  DO  L=0,NL(NY,NX)
!
! LAYER DEPTHS AND THEIR PHYSICAL PROPERTIES

! surface litter:L=0,soil layer:L>0
! DLYR,AREA=layer thickness,x-sectional area:1=EW,2=NS,3=vertical
! ORGC=organic C content
! VOLT,VOLX=volume including,excluding macropores+rock
! BKVL=mass
! CDPTH,DPTH=depth to bottom,midpoint
!
    IF(L.EQ.0)THEN
      ! surface litter residue layer
      TAREA=TAREA+AREA(3,L,NY,NX)
      CumSoilThickness(L,NY,NX)=0.0_r8
      ORGC(L,NY,NX)=SUM(RSC(1:NumOfLitrCmplxs,L,NY,NX))*AREA(3,L,NY,NX)
      ORGCX(L,NY,NX)=ORGC(L,NY,NX)
      VLitR0=0._r8
      DO K=1,NumOfLitrCmplxs
        VLitR0=VLitR0+RSC(K,L,NY,NX)/BulkDensLitR(K)
      ENDDO
      !volume of litter layer
      VLitR(NY,NX)=VLitR0*ppmc*AREA(3,L,NY,NX)
      VGeomLayer(L,NY,NX)=VLitR(NY,NX)
      VLSoilPoreMicP_vr(L,NY,NX)=VGeomLayer(L,NY,NX)
      VLSoilMicP(L,NY,NX)=VLSoilPoreMicP_vr(L,NY,NX)
      VGeomLayert0(L,NY,NX)=VGeomLayer(L,NY,NX)
      SoilMicPMassLayer(L,NY,NX)=MWC2Soil*ORGC(L,NY,NX)  !mass of soil layer, Mg/d2
      !thickness of litter layer 
      DLYRI(3,L,NY,NX)=VLSoilPoreMicP_vr(L,NY,NX)/AREA(3,L,NY,NX)  
      DLYR(3,L,NY,NX)=DLYRI(3,L,NY,NX)
    ELSE
!     if it is a standing water, no macropore fraction
!     DPTH=depth of layer middle
!     CumSoilThickness=soil thickness from surface to bottom of layer L, [m]
!     FracSoiAsMicP=micropore fraction
!     DPTHZ=depth to middle of soil layer from  surface of grid cell [m]
!     VOLT=total volume
!     VOLX=total micropore volume
      IF(SoiBulkDensityt0(L,NY,NX).LE.ZERO)SoilFracAsMacP(L,NY,NX)=0.0_r8
!     thickness:=bottom depth-upper depth
      DLYRI(3,L,NY,NX)=(CumDepth2LayerBottom(L,NY,NX)-CumDepth2LayerBottom(L-1,NY,NX))
      call check_bool(DLYRI(3,L,NY,NX)<0._r8,'negative soil layer thickness',&
        __LINE__,mod_filename)
      DLYR(3,L,NY,NX)=DLYRI(3,L,NY,NX)
      SoiDepthMidLay(L,NY,NX)=0.5_r8*(CumDepth2LayerBottom(L,NY,NX)+CumDepth2LayerBottom(L-1,NY,NX))
      CumSoilThickness(L,NY,NX)=CumDepth2LayerBottom(L,NY,NX)-CumDepth2LayerBottom(NU(NY,NX),NY,NX)+DLYR(3,NU(NY,NX),NY,NX)
      DPTHZ(L,NY,NX)=0.5_r8*(CumSoilThickness(L,NY,NX)+CumSoilThickness(L-1,NY,NX))
      VGeomLayer(L,NY,NX)=amax1(AREA(3,L,NY,NX)*DLYR(3,L,NY,NX),1.e-8_r8)
      VLSoilPoreMicP_vr(L,NY,NX)=VGeomLayer(L,NY,NX)*FracSoiAsMicP(L,NY,NX)
      VLSoilMicP(L,NY,NX)=VLSoilPoreMicP_vr(L,NY,NX)
      VGeomLayert0(L,NY,NX)=VGeomLayer(L,NY,NX)
!     bulk density is defined only for soil with micropores
!     bulk soil mass evaluated as micropore volume
      SoilMicPMassLayer(L,NY,NX)=SoiBulkDensity(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
      RTDNT(L,NY,NX)=0.0_r8
    ENDIF
    AREA(1,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(2,L,NY,NX)
    AREA(2,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(1,L,NY,NX)
  ENDDO
  CumDepth2LayerBottom(0,NY,NX)=CumDepth2LayerBottom(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX)
  CumSoilDeptht0(NY,NX)=CumDepth2LayerBottom(0,NY,NX)
  AREA(3,NL(NY,NX)+1:JZ,NY,NX)=DLYR(1,NL(NY,NX),NY,NX)*DLYR(2,NL(NY,NX),NY,NX)

  end associate
  end subroutine InitLayerDepths

!------------------------------------------------------------------------------------------

  subroutine set_ecosim_solver(NPXS1,NPYS1,NCYC_LITR,NCYC_SNOW)

  implicit none
  integer, intent(in) :: NPXS1,NPYS1
  integer, intent(in) :: NCYC_LITR,NCYC_SNOW
  !     begin_execution
  real(r8) :: XNPV

  NPX=NPXS1   !number of cycles per hour for water,heat,solute flux calcns
  NPY=NPYS1   !number of cycles per NPX for gas flux calcns

  NPH=NPX
  NPT=NPY
  NPG=NPH*NPT

  NPR=NCYC_LITR     !sub-cycles of litter
  NPS=NCYC_SNOW     !sub-cycles of snow iteration
  dts_HeatWatTP=1.0_r8/NPH
  dt_GasCyc=1.0_r8/NPT
  dts_gas=1.0_r8/NPG
  XNPR=1.0_r8/NPR
  XNPS=1.0_r8/NPS
  dts_snohttp=dts_HeatWatTP*XNPS
  dts_litrhtwtp=dts_HeatWatTP*XNPR
  dts_litrvapht=dts_litrhtwtp*XNPS
  XNPV=XNPR*XNPS
  XNPD=600.0_r8*dts_gas
  dts_wat=AMIN1(1.0_r8,5.0_r8*dts_HeatWatTP)  !adjust/recompute the time step for water/heat update, no greater than 1 hour
  dts_sno=dts_wat*XNPS
  XNPB=dts_wat*XNPR
  dt_watvap=dts_wat*XNPV

  end subroutine set_ecosim_solver

!------------------------------------------------------------------------------------------

  subroutine startsim(NHW,NHE,NVN,NVS)
  use SoilHydroParaMod, only : ComputeSoilHydroPars
  use SoilPhysParaMod, only : SetDeepSoil
  implicit none
  integer,intent(in) :: NHW,NHE,NVN,NVS
  REAL(R8) :: ALTY
  real(r8) :: ALTZG
  real(r8) :: tPBOT
  integer :: NY,NX,NM
  real(r8) :: LandScape1stSoiLayDepth
  real(r8) :: YSIN(NumOfSkyAzimuSects),YCOS(NumOfSkyAzimuSects),SkyAzimuthAngle(NumOfSkyAzimuSects)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      NM=MaxNumRootLays(NY,NX)+1
      call ComputeSoilHydroPars(NY,NX,NU(NY,NX),NM)
      call SetDeepSoil(NY,NX,NM,JZ)
    enddo
  enddo
  !  Initialize controlling parameters
  call InitControlParms

  !  IRRADIANCE INTERCEPTION GEOMETRY, plant model
  call InitIrradianceGeometry(YSIN,YCOS,SkyAzimuthAngle)
  !  CALCULATE ELEVATION OF EACH GRID CELL
  !
  call InitGridElevation(NHW,NHE,NVN,NVS,YSIN,YCOS,SkyAzimuthAngle,ALTY)


  call InitAccumulators()

  ALTZG=0.0_r8
  LandScape1stSoiLayDepth=0.0_r8   !pay attention to how it is set for many-grid simulations

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      !     MINIMUM SURFACE ELEVATION IN LANDSCAPE
      !
      !     ALT=surface elevation relative to maximum
      !     ALTZG=minimum surface elevation in landscape
      !CDPTH: depth to bottom of soil layer
      ALT(NY,NX)=ALT(NY,NX)-ALTY
      IF(NX.EQ.NHW.AND.NY.EQ.NVN)THEN
        ALTZG=ALT(NY,NX)
      ELSE
        ALTZG=MIN(ALTZG,ALT(NY,NX))
      ENDIF
      !
      LandScape1stSoiLayDepth=AMAX1(LandScape1stSoiLayDepth,CumDepth2LayerBottom(NU(NY,NX),NY,NX)) !topsoil layer depth
!
!     INITIALIZE ATMOSPHERE VARIABLES
!
!     C*E=atmospheric concentration (g m-3)
!     *E=atmospheric concentration from readi.f (umol mol-1)
!     CO2=CO2,CH4=CH4,OXY=O2,Z2G=N2,Z2O=N2O,NH3=NH3,H2G=H2
!     TairKClimMean=mean annual air temperature (K)
!
      tPBOT=1._r8
      CCO2EI(NY,NX)=CO2EI(NY,NX)*5.36E-04_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_CO2,NY,NX)=CO2E(NY,NX)*5.36E-04_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_CH4,NY,NX)=CH4E(NY,NX)*5.36E-04_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_O2,NY,NX)=OXYE(NY,NX)*1.43E-03_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_N2,NY,NX)=Z2GE(NY,NX)*1.25E-03_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_N2O,NY,NX)=Z2OE(NY,NX)*1.25E-03_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_NH3,NY,NX)=ZNH3E(NY,NX)*6.25E-04_r8*Tref/TairKClimMean(NY,NX)*tPBOT
      AtmGasCgperm3(idg_H2,NY,NX)=H2GE(NY,NX)*8.92E-05_r8*Tref/TairKClimMean(NY,NX)*tPBOT
!
!     MICROBIAL THERMAL ADAPTATION
!
!     OFFSET=shift in Arrhenius curve used in nitro.f (oC)
!     ATCS=mean annual soil temperature (OC)
!
      OFFSET(NY,NX)=fOFFSET(ATCS(NY,NX))

!
!     INITIALIZE WATER POTENTIAL VARIABLES FOR SOIL LAYERS
!
!     LOGPSIMX,PSIMN,LOGPSIAtSat=log water potential at FC,WP,POROS
!     LOGPSIMXD,LOGPSIMND=LOGPSIMX-LOGPSIAtSat,PSIMN-LOGPSIMX
!
      LOGPSIAtSat(NY,NX)=LOG(-PSIPS)
      LOGPSIFLD(NY,NX)=LOG(-PSIAtFldCapacity(NY,NX))
      LOGPSIMN(NY,NX)=LOG(-PSIAtWiltPoint(NY,NX))
      LOGPSIMXD(NY,NX)=LOGPSIFLD(NY,NX)-LOGPSIAtSat(NY,NX)
      LOGPSIMND(NY,NX)=LOGPSIMN(NY,NX)-LOGPSIFLD(NY,NX)
!     SoilMicPMassLayer(0,NY,NX)=0.0_r8
!
!     DISTRIBUTION OF OM AMONG FRACTIONS OF DIFFERING
!     BIOLOGICAL ACTIVITY
!

      call InitLayerDepths(NY,NX)
    ! ActiveLayDepth=active layer depth (m)
      ActiveLayDepth(NY,NX)=9999.0_r8
!
!     INITIALIZE SNOWPACK LAYERS
      call InitSnowLayers(NY,NX)

!     VHeatCapLitR,VHCPNX=minimum heat capacities for solving
!      surface litter,soil layer water and heat fluxes
      VHeatCapLitR(NY,NX)=VLHeatCapLitRMin*AREA(3,NU(NY,NX),NY,NX)
      VHCPNX(NY,NX)=VLHeatCapSoiMin*AREA(3,NU(NY,NX),NY,NX)

!
!     SURFACE WATER STORAGE AND LOWER HEAT SINK

!      SoilHeatSrcDepth=depth at which soil heat sink-source calculated
!     TCNDG=assumed thermal conductivity below lower soil boundary
!     (MJ m-1 K-1 h-1)
!     TKSD=deep source/sink temperature from geothermal flux(K)

      SoilHeatSrcDepth(NY,NX)=AMAX1(10.0_r8,CumDepth2LayerBottom(NL(NY,NX),NY,NX)+1.0_r8)
      TCS(0,NY,NX)=ATCS(NY,NX)
      TKS(0,NY,NX)=ATKS(NY,NX)
      TKSD(NY,NX)=ATKS(NY,NX)+2.052E-04_r8*SoilHeatSrcDepth(NY,NX)/TCNDG
!

    ENDDO
  ENDDO

!
!     INITIALIZE COMMUNITY CANOPY
!
  CanopyHeight_col(:,:)=0.0_r8
  CanopyHeightZ_col(0:NumOfCanopyLayers,:,:)=0.0_r8
  CanopyLeafAareZ_col(1:NumOfCanopyLayers,:,:)=0.0_r8
  CanopyStemAareZ_col(1:NumOfCanopyLayers,:,:)=0.0_r8
  WGLFT(1:NumOfCanopyLayers,:,:)=0.0_r8

!
  call InitSoilVars(NHW,NHE,NVN,NVS,ALTZG,LandScape1stSoiLayDepth)
  end subroutine startsim

end module StartsMod
