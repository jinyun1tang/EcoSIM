module StartsMod
!!
! Description:
! code to initalize soil variables

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : padr, print_info,check_bool
  use minimathMod, only : test_aeqb, test_aneb, AZMAX1,AZMIN1
  use EcosimConst
  use TracerIDMod
  use SnowDataType  
  use MicrobialDataType
  use EcoSIMSolverPar
  use SOMDataType
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
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = __FILE__
  !
  !
  !     BKRS=dry bulk density of woody(0),fine(1),manure(2) litter
  !     FORGC=minimum SOC for organic soil (g Mg-1)
  !      FVLWB,FCH4F=maximum SWC,CH4 emission fraction for combustion
  !     PSIHY=hygroscopic water potential (MPa)
  !     FCI,WPI=FC,WP for water retention by ice (MPa)
  !     CDPTHSI=depth to bottom of snowpack layers
  !     POROQ=Penman Water Linear Reduction tortuosity used in gas flux calculations
  !

  public :: starts
  public :: set_ecosim_solver
  contains

  SUBROUTINE starts(NHW,NHE,NVN,NVS)
  !
  !     THIS SUBROUTINE INITIALIZES ALL SOIL VARIABLES

  use InitVegBGC, only : InitIrradianceGeometry
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX,L,NGL
  REAL(R8) :: ALTY
  real(r8) :: ALTZG
  real(r8) :: tPBOT
  real(r8) :: CDPTHG
  real(r8) :: YSIN(JSA),YCOS(JSA),YAZI(JSA)
! begin_execution


  !  Initialize controlling parameters
  call InitControlParms
  !
  !     IRRADIANCE INTERCEPTION GEOMETRY
  call InitIrradianceGeometry(YSIN,YCOS,YAZI)
  !     CALCULATE ELEVATION OF EACH GRID CELL
  !
  call InitGridElevation(NHW,NHE,NVN,NVS,YSIN,YCOS,YAZI,ALTY)
  !
  !     INITIALIZE ACCUMULATORS AND MASS BALANCE CHECKS
  !     OF EACH GRID CELL
  !
  call InitAccumulators()

  ALTZG=0.0_r8
  CDPTHG=0.0_r8
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
      CDPTHG=AMAX1(CDPTHG,CDPTH(NU(NY,NX),NY,NX)) !topsoil layer depth
!
!     INITIALIZE ATMOSPHERE VARIABLES
!
!     C*E=atmospheric concentration (g m-3)
!     *E=atmospheric concentration from readi.f (umol mol-1)
!     CO2=CO2,CH4=CH4,OXY=O2,Z2G=N2,Z2O=N2O,NH3=NH3,H2G=H2
!     ATKA=mean annual air temperature (K)
!
      tPBOT=PBOT(NY,NX)/1.01325E+02_r8
      CCO2EI(NY,NX)=CO2EI(NY,NX)*5.36E-04_r8*Tref/ATKA(NY,NX)*tPBOT
      AtmGgms(idg_CO2,NY,NX)=CO2E(NY,NX)*5.36E-04_r8*Tref/ATKA(NY,NX)*tPBOT
      AtmGgms(idg_CH4,NY,NX)=CH4E(NY,NX)*5.36E-04_r8*Tref/ATKA(NY,NX)*tPBOT
      AtmGgms(idg_O2,NY,NX)=OXYE(NY,NX)*1.43E-03_r8*Tref/ATKA(NY,NX)*tPBOT
      AtmGgms(idg_N2,NY,NX)=Z2GE(NY,NX)*1.25E-03_r8*Tref/ATKA(NY,NX)*tPBOT
      AtmGgms(idg_N2O,NY,NX)=Z2OE(NY,NX)*1.25E-03_r8*Tref/ATKA(NY,NX)*tPBOT
      AtmGgms(idg_NH3,NY,NX)=ZNH3E(NY,NX)*6.25E-04_r8*Tref/ATKA(NY,NX)*tPBOT
      AtmGgms(idg_H2,NY,NX)=H2GE(NY,NX)*8.92E-05_r8*Tref/ATKA(NY,NX)*tPBOT
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
!     PSIMX,PSIMN,PSIMS=log water potential at FC,WP,POROS
!     PSISD,PSIMD=PSIMX-PSIMS,PSIMN-PSIMX
!
      PSIMS(NY,NX)=LOG(-PSIPS)
      PSIMX(NY,NX)=LOG(-PSIFC(NY,NX))
      PSIMN(NY,NX)=LOG(-PSIWP(NY,NX))
      PSISD(NY,NX)=PSIMX(NY,NX)-PSIMS(NY,NX)
      PSIMD(NY,NX)=PSIMN(NY,NX)-PSIMX(NY,NX)
!     BKVL(0,NY,NX)=0.0_r8
!
!     DISTRIBUTION OF OM AMONG FRACTIONS OF DIFFERING
!     BIOLOGICAL ACTIVITY
!
      call InitLayerDepths(NY,NX)
    ! DPTHA=active layer depth (m)
      DPTHA(NY,NX)=9999.0_r8
!
!     INITIALIZE SNOWPACK LAYERS
      call InitSnowLayers(NY,NX)

!     VHCPRX,VHCPNX=minimum heat capacities for solving
!      surface litter,soil layer water and heat fluxes
      VHCPRX(NY,NX)=VHCPRMin*AREA(3,NU(NY,NX),NY,NX)
      VHCPNX(NY,NX)=VHCPNMin*AREA(3,NU(NY,NX),NY,NX)

!
!     SURFACE WATER STORAGE AND LOWER HEAT SINK

!      DPTHSK=depth at which soil heat sink-source calculated
!     TCNDG=assumed thermal conductivity below lower soil boundary
!     (MJ m-1 K-1 h-1)
!     TKSD=deep source/sink temperature from geothermal flux(K)
      
      DPTHSK(NY,NX)=AMAX1(10.0_r8,CDPTH(NL(NY,NX),NY,NX)+1.0_r8)
      TCS(0,NY,NX)=ATCS(NY,NX)
      TKS(0,NY,NX)=ATKS(NY,NX)
      TKSD(NY,NX)=ATKS(NY,NX)+2.052E-04_r8*DPTHSK(NY,NX)/TCNDG
!

    ENDDO
  ENDDO

!
!     INITIALIZE COMMUNITY CANOPY
!
  ZT(:,:)=0.0_r8
  ZL(0:JC,:,:)=0.0_r8
  ARLFT(1:JC,:,:)=0.0_r8
  ARSTT(1:JC,:,:)=0.0_r8
  WGLFT(1:JC,:,:)=0.0_r8

!
  call InitSoilVars(NHW,NHE,NVN,NVS,ALTZG,CDPTHG)

  RETURN
  END subroutine starts
!------------------------------------------------------------------------------------------
  subroutine InitSoilVars(NHW,NHE,NVN,NVS,ALTZG,CDPTHG)
  !     N3,N2,N1=L,NY,NX of source grid cell
  !     N6,N5,N4=L,NY,NX of destination grid cell
  !      ALTZG=minimum surface elevation in landscape
  !     DTBLI,DTBLDI=depth of natural,artificial water table
  !     DTBLG=slope of natural water table relative to landscape surface
  !     in geography, slope =rise/run
  !     DTBLZ,DTBLD=depth of natural,artificial water table adjusted for elevn
  !     DPTHT=depth to internal water table
  !     DIST=distance between adjacent layers:1=EW,2=NS,3=vertical(m)
  !     XDPTH=x-section area/distance in solute flux calculations (m2/m)
  !     DISP=dispersivity parameter in solute flux calculations (m2 h-1)
  !
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: ALTZG,CDPTHG
  integer :: NY,NX,L,N
  integer :: N1,N2,N3,N4,N5,N6


  ! begin_execution

!     INITIALIZE SEDIMENT LOAD IN EROSION MODEL
!
  IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
    SED(:,:)=0.0_r8
  ENDIF
!
!     INITIALIZE GRID CELL DIMENSIONS
! altz: topographic altitude
! ALT: grid altitude
! DTBLI: external water table depth, before applying the altitude correction
! DTBLG: slope of water table relative to surface slope
! DTBLZ: external water table depth
! DTBLDI: depth of artificial water table
! DPTHT: internal water table depth
! DTBLD: artifical water table depth, before applying the altitude correction
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      ALTZ(NY,NX)=ALTZG
      IF(BKDS(NU(NY,NX),NY,NX).GT.0.0_r8)THEN
        DTBLZ(NY,NX)=DTBLI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))*(1.0_r8-DTBLG(NY,NX))
        DTBLD(NY,NX)=AZMAX1(DTBLDI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))*(1.0_r8-DTBLG(NY,NX)))
      ELSE
        DTBLZ(NY,NX)=0.0_r8
        DTBLD(NY,NX)=0.0_r8
      ENDIF
      DPTHT(NY,NX)=DTBLZ(NY,NX)
    ENDDO
  ENDDO

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      DO L=1,NL(NY,NX)
        N1=NX
        N2=NY
        N3=L
        DO N=NCN(N2,N1),3
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
      call InitSoilProfile(NY,NX,CDPTHG)
      !
      !     SURFACE LITTER HEAT CAPACITY
      !
      BKVLNM(NY,NX)=AZMAX1(SAND(NU(NY,NX),NY,NX)+SILT(NU(NY,NX),NY,NX)+CLAY(NU(NY,NX),NY,NX))
      VHCP(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VOLW(0,NY,NX)+cpi*VOLI(0,NY,NX)
      VHCM(0,NY,NX)=0.0_r8
      VOLAI(0,NY,NX)=0.0_r8
    ENDDO
  ENDDO
  end subroutine InitSoilVars
!------------------------------------------------------------------------------------------
  subroutine InitSoilProfile(NY,NX,CDPTHG)

  implicit none
  integer, intent(in) :: NY,NX
  REAL(R8),INTENT(IN) :: CDPTHG
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
    TORGL(L)=TORGC+CORGL*BKVL(L,NY,NX)/AREA(3,L,NY,NX)*0.5
    TORGC=TORGC+CORGL*BKVL(L,NY,NX)/AREA(3,L,NY,NX)
  ENDDO D1190
!
!     PARAMETERS TO ALLOCATE HUMUS TO LESS OR MORE RECALCITRANT FRACTIONS
!
!     TORGL=accumulated humus down to soil layer (g m-2)
!     TORGM=TORGL used to calculate allocation (g m-2)
!     HCX=shape parameter for depth effect on allocation
!
  TORGM=AMAX1(2.0E+03_r8,AMIN1(5.0E+03_r8,0.25_r8*TORGL(NJ(NY,NX))))
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

    call InitSOMProfile(L,NY,NX,HCX,CDPTHG,TORGLL,CORGCM,FCX)
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
    PSISA(L,NY,NX)=-1.5E-03_r8
    ROXYF(L,NY,NX)=0.0_r8
    RCO2F(L,NY,NX)=0.0_r8
    ROXYL(L,NY,NX)=0.0_r8
    RCH4F(L,NY,NX)=0.0_r8
    RCH4L(L,NY,NX)=0.0_r8

    IF(L.GT.0)THEN
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
        !it is a soil layer
        !compute particle density
        PTDS=ppmc*(1.30_r8*CORGCM+2.66_r8*(1.0E+06_r8-CORGCM))
        POROS(L,NY,NX)=1.0_r8-(BKDS(L,NY,NX)/PTDS)
      ELSE
        !for ponding water
        PTDS=0.0_r8
        POROS(L,NY,NX)=1.0_r8
      ENDIF
      POROSI(L,NY,NX)=POROS(L,NY,NX)*FMPR(L,NY,NX)
      VOLA(L,NY,NX)=POROS(L,NY,NX)*VOLX(L,NY,NX)
      VOLAI(L,NY,NX)=VOLA(L,NY,NX)
      VOLAH(L,NY,NX)=FHOL(L,NY,NX)*VOLTI(L,NY,NX)
      !
      !     LAYER HEAT CONTENTS
      !
      !     SAND,SILT,CLAY=sand,silt,clay mass (Mg)
      !     VORGC,VMINL,VSAND=volumetric fractions of SOC,non-sand,sand
      !     VHCM,VHCP=volumetric dry,wet soil heat capacity (MJ m-3 K-1)
      !     TKS,TCS=soil temperature (oC,K)
      !     THETW,THETI,THETP=micropore water,ice,air concentration (m3 m-3)
!
      SAND(L,NY,NX)=CSAND(L,NY,NX)*BKVL(L,NY,NX)
      SILT(L,NY,NX)=CSILT(L,NY,NX)*BKVL(L,NY,NX)
      CLAY(L,NY,NX)=CCLAY(L,NY,NX)*BKVL(L,NY,NX)
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
        ! PTDS=particle density (Mg m-3)
        ! soil volumetric heat capacity
        VORGC=CORGCM*ppmc*BKDS(L,NY,NX)/PTDS
        VMINL=(CSILT(L,NY,NX)+CCLAY(L,NY,NX))*BKDS(L,NY,NX)/PTDS
        VSAND=CSAND(L,NY,NX)*BKDS(L,NY,NX)/PTDS
        VHCM(L,NY,NX)=((2.496_r8*VORGC+2.385_r8*VMINL+2.128_r8*VSAND) &
          *FMPR(L,NY,NX)+2.128_r8*ROCK(L,NY,NX))*VOLT(L,NY,NX)
      ELSE
        VHCM(L,NY,NX)=0.0_r8
      ENDIF
!
      !     INITIAL SOIL WATER AND ICE CONTENTS
!
      IF(ISOIL(isoi_fc,L,NY,NX).EQ.0.AND.ISOIL(isoi_wp,L,NY,NX).EQ.0)THEN
      ! field capacity and wilting point are read from input
        IF(THW(L,NY,NX).GT.1.0)THEN
          THETW(L,NY,NX)=POROS(L,NY,NX)
        ELSEIF(test_aeqb(THW(L,NY,NX),1.0_r8))THEN
          THETW(L,NY,NX)=FC(L,NY,NX)
        ELSEIF(test_aeqb(THW(L,NY,NX),0.0_r8))THEN
          THETW(L,NY,NX)=WP(L,NY,NX)
        ELSEIF(THW(L,NY,NX).LT.0.0)THEN
          THETW(L,NY,NX)=0.0_r8
        ELSE
          THETW(L,NY,NX)=THW(L,NY,NX)
        ENDIF
        IF(THI(L,NY,NX).GT.1.0_r8)THEN
          THETI(L,NY,NX)=AZMAX1(AMIN1(POROS(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
        ELSEIF(test_aeqb(THI(L,NY,NX),1.0_r8))THEN
          THETI(L,NY,NX)=AZMAX1(AMIN1(FC(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
        ELSEIF(test_aeqb(THI(L,NY,NX),0.0_r8))THEN
          THETI(L,NY,NX)=AZMAX1(AMIN1(WP(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
        ELSEIF(THI(L,NY,NX).LT.0.0_r8)THEN
          THETI(L,NY,NX)=0.0_r8
        ELSE
          THETI(L,NY,NX)=THI(L,NY,NX)
        ENDIF
        VOLW(L,NY,NX)=THETW(L,NY,NX)*VOLX(L,NY,NX)
        VOLWX(L,NY,NX)=VOLW(L,NY,NX)
        VOLWH(L,NY,NX)=THETW(L,NY,NX)*VOLAH(L,NY,NX)
        VOLI(L,NY,NX)=THETI(L,NY,NX)*VOLX(L,NY,NX)
        VOLIH(L,NY,NX)=THETI(L,NY,NX)*VOLAH(L,NY,NX)
!       VOLP: total air-filled porosity, micropores + macropores
        VOLP(L,NY,NX)=AZMAX1(VOLA(L,NY,NX)-VOLW(L,NY,NX) &
          -VOLI(L,NY,NX))+AZMAX1(VOLAH(L,NY,NX)-VOLWH(L,NY,NX) &
          -VOLIH(L,NY,NX))
        VHCP(L,NY,NX)=VHCM(L,NY,NX)+cpw*(VOLW(L,NY,NX) &
          +VOLWH(L,NY,NX))+cpi*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
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
  trcs_VLN(ids_NH4,0:L2,NY,NX)=1.0_r8
  trcs_VLN(idg_NH3,0:L2,NY,NX)=trcs_VLN(ids_NH4,0:L2,NY,NX)
  trcs_VLN(ids_NO3,0:L2,NY,NX)=1.0_r8
  trcs_VLN(ids_NO2,0:L2,NY,NX)=trcs_VLN(ids_NO3,0:L2,NY,NX)
  trcs_VLN(ids_H1PO4,0:L2,NY,NX)=1.0_r8
  trcs_VLN(ids_H2PO4,0:L2,NY,NX)=trcs_VLN(ids_H1PO4,0:L2,NY,NX)
  trcs_VLN(ids_NH4B,0:L2,NY,NX)=0.0_r8
  trcs_VLN(idg_NH3B,0:L2,NY,NX)= trcs_VLN(ids_NH4B,0:L2,NY,NX)
  trcs_VLN(ids_NO3B,0:L2,NY,NX)=0.0_r8
  trcs_VLN(ids_NO2B,0:L2,NY,NX)=trcs_VLN(ids_NO3B,0:L2,NY,NX)
  trcs_VLN(ids_H1PO4B,0:L2,NY,NX)=0.0_r8
  trcs_VLN(ids_H2PO4B,0:L2,NY,NX)=trcs_VLN(ids_H1PO4B,0:L2,NY,NX)

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
  subroutine InitGridElevation(NHW,NHE,NVN,NVS,YSIN,YCOS,YAZI,ALTY)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: YSIN(JSA),YCOS(JSA),YAZI(JSA)
  REAL(R8),INTENT(OUT):: ALTY
  integer :: NY,NX,N,NN
  REAL(R8) :: DGAZI
  real(r8) :: GSINA(JY,JX),GCOSA(JY,JX)  !diagnostic

! begin_execution
! GAZI=ground surface azimuth, aspect in radians
! GSIN,GCOS=sine,cosine of ground surface
! OMEGAG=incident sky angle at ground surface
! SLOPE=sine of ground surface slope in (0)aspect, (1)EW,(2)NS directions
! ALT=ground surface elevation
! ALTY=maximum surface elevation in landscape
! IRCHG=runoff boundary flags:0=not possible,1=possible
! ASP=aspect angle in degree
  ALTY=0.0_r8
  D9985: DO NX=NHW,NHE
    D9980: DO NY=NVN,NVS
      ZEROS(NY,NX)=ZERO*DH(NY,NX)*DV(NY,NX)
      ZEROS2(NY,NX)=ZERO2*DH(NY,NX)*DV(NY,NX)
!     compute slopes
      GAZI(NY,NX)=ASP(NY,NX)/RDN   !radian
      GSINA(NY,NX)=ABS(SIN(GAZI(NY,NX)))
      GCOSA(NY,NX)=ABS(COS(GAZI(NY,NX)))
      SLOPE(0,NY,NX)=AMAX1(1.745E-04_r8,SIN(SL(NY,NX)/RDN))  !small slope approximation

      IF(ASP(NY,NX).GE.0.0_r8.AND.ASP(NY,NX).LT.90.0_r8)THEN
      ! along the northeast
        SLOPE(1,NY,NX)=-SLOPE(0,NY,NX)*COS(ASP(NY,NX)/RDN)    !to east
        SLOPE(2,NY,NX)=SLOPE(0,NY,NX)*SIN(ASP(NY,NX)/RDN)     !to north
        IRCHG(1,1,NY,NX)=1    !east
        IRCHG(2,1,NY,NX)=0
        IRCHG(1,2,NY,NX)=0
        IRCHG(2,2,NY,NX)=1    !north
      ELSEIF(ASP(NY,NX).GE.90.0_r8.AND.ASP(NY,NX).LT.180.0_r8)THEN
      ! along the northwest
        SLOPE(1,NY,NX)=SLOPE(0,NY,NX)*SIN((ASP(NY,NX)-90.0_r8)/RDN)   !to west
        SLOPE(2,NY,NX)=SLOPE(0,NY,NX)*COS((ASP(NY,NX)-90.0_r8)/RDN)   !to north
        IRCHG(1,1,NY,NX)=0
        IRCHG(2,1,NY,NX)=1   !west
        IRCHG(1,2,NY,NX)=0
        IRCHG(2,2,NY,NX)=1   !north
      ELSEIF(ASP(NY,NX).GE.180.0_r8.AND.ASP(NY,NX).LT.270.0_r8)THEN
      !along the southwest
        SLOPE(1,NY,NX)=SLOPE(0,NY,NX)*COS((ASP(NY,NX)-180.0_r8)/RDN)    !to west
        SLOPE(2,NY,NX)=-SLOPE(0,NY,NX)*SIN((ASP(NY,NX)-180.0_r8)/RDN)   !to south
        IRCHG(1,1,NY,NX)=0
        IRCHG(2,1,NY,NX)=1  !west
        IRCHG(1,2,NY,NX)=1  !south
        IRCHG(2,2,NY,NX)=0
      ELSEIF(ASP(NY,NX).GE.270.0_r8.AND.ASP(NY,NX).LE.360.0_r8)THEN
      ! along the southeast
        SLOPE(1,NY,NX)=-SLOPE(0,NY,NX)*SIN((ASP(NY,NX)-270.0_r8)/RDN)   !to east
        SLOPE(2,NY,NX)=-SLOPE(0,NY,NX)*COS((ASP(NY,NX)-270.0_r8)/RDN)   !to south
        IRCHG(1,1,NY,NX)=1  !east
        IRCHG(2,1,NY,NX)=0
        IRCHG(1,2,NY,NX)=1  !south
        IRCHG(2,2,NY,NX)=0
      ENDIF
      SLOPE(3,NY,NX)=-1.0_r8

      IF(test_aneb(SLOPE(1,NY,NX),0.0_r8).OR.test_aneb(SLOPE(2,NY,NX),0.0_r8))THEN
        FSLOPE(1,NY,NX)=ABS(SLOPE(1,NY,NX))/(ABS(SLOPE(1,NY,NX))+ABS(SLOPE(2,NY,NX)))  !
        FSLOPE(2,NY,NX)=ABS(SLOPE(2,NY,NX))/(ABS(SLOPE(1,NY,NX))+ABS(SLOPE(2,NY,NX)))
      ELSE
        FSLOPE(1,NY,NX)=0.5_r8
        FSLOPE(2,NY,NX)=0.5_r8
      ENDIF

!    compute incident sky angle at ground surface
      GSIN(NY,NX)=SLOPE(0,NY,NX)    !this is exact
      GCOS(NY,NX)=SQRT(1.0_r8-GSIN(NY,NX)**2._r8)
      D240: DO N=1,JSA
        DGAZI=COS(GAZI(NY,NX)-YAZI(N))
        OMEGAG(N,NY,NX)=AZMAX1(AMIN1(1.0_r8,GCOS(NY,NX)*YSIN(N)+GSIN(NY,NX)*YCOS(N)*DGAZI))
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
      WRITE(*,1111)'ALT',NX,NY,((IRCHG(NN,N,NY,NX),NN=1,2),N=1,2) &
        ,ALT(NY,NX),DH(NY,NX),DV(NY,NX),ASP(NY,NX),SL(NY,NX) &
        ,SLOPE(0,NY,NX),SLOPE(1,NY,NX),SLOPE(2,NY,NX) &
        ,GSIN(NY,NX),GCOSA(NY,NX),GSINA(NY,NX)
1111  FORMAT(A8,6I4,20E12.4)
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
  !     THETPI,DENSI=ice porosity,density
  !     BKRS=surface litter bulk density, Mg m-3

  BKRS=(/0.0333_r8,0.0167_r8,0.0167_r8/)

  call InitSOMConsts
  !     NDIM=1
  !     IF(NHE.GT.NHW)NDIM=NDIM+1
  !     IF(NVS.GT.NVN)NDIM=NDIM+1
  !     XDIM=1.0/NDIM
  ZERO=1.0E-15_r8
  ZERO2=1.0E-08_r8
  TAREA=0.0_r8
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
  TCOU=0.0_r8
  TZOU=0.0_r8
  TPOU=0.0_r8
  XCSN=0.0_r8
  XZSN=0.0_r8
  XPSN=0.0_r8
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
  TDCO2(:,:,:)=1.0_r8
  TDCN4(:,:,:)=1.0_r8
  TDCNO(:,:,:)=1.0_r8

  IUTYP(:,:)=0
  IFNHB(:,:)=0
  IFNOB(:,:)=0
  IFPOB(:,:)=0
  IFLGS(:,:)=1
  IFLGT(:,:)=0
  ATCA(:,:)=ATCAI(:,:)
  ATCS(:,:)=ATCAI(:,:)
  ATKA(:,:)=ATCA(:,:)+TC2K
  ATKS(:,:)=ATCS(:,:)+TC2K
  URAIN(:,:)=0.0_r8
  UCO2G(:,:)=0.0_r8
  UCH4G(:,:)=0.0_r8
  UOXYG(:,:)=0.0_r8
  UN2GG(:,:)=0.0_r8
  UN2OG(:,:)=0.0_r8
  UNH3G(:,:)=0.0_r8
  UN2GS(:,:)=0.0_r8
  UCO2F(:,:)=0.0_r8
  UCH4F(:,:)=0.0_r8
  UOXYF(:,:)=0.0_r8
  UN2OF(:,:)=0.0_r8
  UNH3F(:,:)=0.0_r8
  UPO4F(:,:)=0.0_r8
  UORGF(:,:)=0.0_r8
  UFERTN(:,:)=0.0_r8
  UFERTP(:,:)=0.0_r8
  UVOLO(:,:)=0.0_r8
  UEVAP(:,:)=0.0_r8
  URUN(:,:)=0.0_r8
  USEDOU(:,:)=0.0_r8
  UCOP(:,:)=0.0_r8
  UDOCQ(:,:)=0.0_r8
  UDOCD(:,:)=0.0_r8
  UDONQ(:,:)=0.0_r8
  UDOND(:,:)=0.0_r8
  UDOPQ(:,:)=0.0_r8
  UDOPD(:,:)=0.0_r8
  UDICQ(:,:)=0.0_r8
  UDICD(:,:)=0.0_r8
  UDINQ(:,:)=0.0_r8
  UDIND(:,:)=0.0_r8
  UDIPQ(:,:)=0.0_r8
  UDIPD(:,:)=0.0_r8
  UIONOU(:,:)=0.0_r8
  UXCSN(:,:)=0.0_r8
  UXZSN(:,:)=0.0_r8
  UXPSN(:,:)=0.0_r8
  UDRAIN(:,:)=0.0_r8
  ZDRAIN(:,:)=0.0_r8
  PDRAIN(:,:)=0.0_r8
  DPNH4(:,:)=0.0_r8
  DPNO3(:,:)=0.0_r8
  DPPO4(:,:)=0.0_r8
  trc_solml(idg_O2,0,:,:)=0.0_r8
  FRADG(:,:)=1.0_r8
  THRMG(:,:)=0.0_r8
  THRMC(:,:)=0.0_r8
  TRN(:,:)=0.0_r8
  TLE(:,:)=0.0_r8
  TSH(:,:)=0.0_r8
  TGH(:,:)=0.0_r8
  TLEC(:,:)=0.0_r8
  TSHC(:,:)=0.0_r8
  TLEX(:,:)=0.0_r8
  TSHX(:,:)=0.0_r8
  TCNET(:,:)=0.0_r8
  TVOLWC(:,:)=0.0_r8
  ARLFC(:,:)=0.0_r8
  ARSTC(:,:)=0.0_r8
  TFLWC(:,:)=0.0_r8
  PPT(:,:)=0.0_r8
  DYLN(:,:)=12.0_r8
  ALBX(:,:)=ALBS(:,:)
  XHVSTE(:,:,:)=0.0_r8
  ENGYP(:,:)=0.0_r8
  end subroutine InitAccumulators
!------------------------------------------------------------------------------------------
  subroutine InitLayerDepths(NY,NX)
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer, intent(in) :: NY, NX
  integer :: L,K
  real(r8) :: VOLR0
  associate(                              &
    n_litrsfk    => micpar%n_litrsfk    , &
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
    DLYRI(1,L,NY,NX)=DH(NY,NX)
    DLYRI(2,L,NY,NX)=DV(NY,NX)
    DLYR(1,L,NY,NX)=DLYRI(1,L,NY,NX)
    DLYR(2,L,NY,NX)=DLYRI(2,L,NY,NX)
    AREA(3,L,NY,NX)=DLYR(1,L,NY,NX)*DLYR(2,L,NY,NX)
    IF(L.EQ.0)THEN
      ! surface litter residue layer
      TAREA=TAREA+AREA(3,L,NY,NX)
      CDPTHZ(L,NY,NX)=0.0_r8
      ORGC(L,NY,NX)=SUM(RSC(1:n_litrsfk,L,NY,NX))*AREA(3,L,NY,NX)
      ORGCX(L,NY,NX)=ORGC(L,NY,NX)
      VOLR0=0._r8
      DO K=1,n_litrsfk
        VOLR0=VOLR0+RSC(K,L,NY,NX)/BKRS(K)
      ENDDO
      VOLR(NY,NX)=VOLR0*ppmc*AREA(3,L,NY,NX)
      VOLT(L,NY,NX)=VOLR(NY,NX)
      VOLX(L,NY,NX)=VOLT(L,NY,NX)
      VOLY(L,NY,NX)=VOLX(L,NY,NX)
      VOLTI(L,NY,NX)=VOLT(L,NY,NX)
      BKVL(L,NY,NX)=MWC2Soil*ORGC(L,NY,NX)  !mass of soil layer, Mg/d2
      DLYRI(3,L,NY,NX)=VOLX(L,NY,NX)/AREA(3,L,NY,NX)
      DLYR(3,L,NY,NX)=DLYRI(3,L,NY,NX)
    ELSE
!     if it is a standing water, no macropore fraction
      IF(BKDSI(L,NY,NX).LE.ZERO)FHOL(L,NY,NX)=0.0
!     thickness:=bottom depth-upper depth
      DLYRI(3,L,NY,NX)=(CDPTH(L,NY,NX)-CDPTH(L-1,NY,NX))
      call check_bool(DLYRI(3,L,NY,NX)<0._r8,'negative soil layer thickness',&
        __LINE__,mod_filename)
      !FMPR: micropore fraction
      DLYR(3,L,NY,NX)=DLYRI(3,L,NY,NX)
!     DPTH: depth of layer middle
      DPTH(L,NY,NX)=0.5_r8*(CDPTH(L,NY,NX)+CDPTH(L-1,NY,NX))
!     CDPTHZ: soil thickness from surface to bottom of layer L, [m]
      CDPTHZ(L,NY,NX)=CDPTH(L,NY,NX)-CDPTH(NU(NY,NX),NY,NX)+DLYR(3,NU(NY,NX),NY,NX)
!     DPTHZ: depth to middle of soil layer from  surface of grid cell [m]
      DPTHZ(L,NY,NX)=0.5_r8*(CDPTHZ(L,NY,NX)+CDPTHZ(L-1,NY,NX))
!     VOLT: total volume
      VOLT(L,NY,NX)=amax1(AREA(3,L,NY,NX)*DLYR(3,L,NY,NX),1.e-8_r8)
!     VOLX: total micropore volume
      VOLX(L,NY,NX)=VOLT(L,NY,NX)*FMPR(L,NY,NX)
      VOLY(L,NY,NX)=VOLX(L,NY,NX)
      VOLTI(L,NY,NX)=VOLT(L,NY,NX)
!     bulk density is defined only for soil with micropores      
!     bulk soil mass evaluated as micropore volume
      BKVL(L,NY,NX)=BKDS(L,NY,NX)*VOLX(L,NY,NX)
      RTDNT(L,NY,NX)=0.0_r8
    ENDIF
    AREA(1,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(2,L,NY,NX)
    AREA(2,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(1,L,NY,NX)
  ENDDO
  CDPTH(0,NY,NX)=CDPTH(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX)
  CDPTHI(NY,NX)=CDPTH(0,NY,NX)
  AREA(3,NL(NY,NX)+1:JZ,NY,NX)=DLYR(1,NL(NY,NX),NY,NX)*DLYR(2,NL(NY,NX),NY,NX)

  end associate
  end subroutine InitLayerDepths
    
!------------------------------------------------------------------------------------------
  
  subroutine set_ecosim_solver(NPXS1,NPYS1)
  
  implicit none
  integer, intent(in) :: NPXS1,NPYS1
  !     begin_execution
  real(r8) :: XNPV
  
  NPX=NPXS1   !number of cycles per hour for water,heat,solute flux calcns
  NPY=NPYS1   !number of cycles per NPX for gas flux calcns
    
  NPH=NPX
  NPT=NPY
  NPG=NPH*NPT

  NPR=30
  NPS=10
  XNPH=1.0_r8/NPH
  XNPT=1.0_r8/NPT
  XNPG=1.0_r8/NPG
  XNPR=1.0_r8/NPR
  XNPS=1.0_r8/NPS
  XNPY=XNPH*XNPS
  XNPZ=XNPH*XNPR
  XNPQ=XNPZ*XNPS
  XNPV=XNPR*XNPS
  XNPD=600.0_r8*XNPG
  XNPX=AMIN1(1.0_r8,20.0_r8*XNPH)
  XNPA=XNPX*XNPS
  XNPB=XNPX*XNPR
  XNPC=XNPX*XNPV
  
  end subroutine set_ecosim_solver

end module StartsMod
