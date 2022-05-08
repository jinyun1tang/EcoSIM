module StartsMod
!!
! Description:
! code to initalize soil variables

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : padr, print_info,check_bool
  use minimathMod, only : test_aeqb, test_aneb
  use EcosimConst
  use MicrobialDataType
  use EcoSIMSolverPar
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use InitSOMBGC
  use CanopyRadDataType
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use LandSurfDataType
  use SnowDataType
  use PlantTraitDataType
  use PlantDataRateType
  use SurfLitterDataType
  use CanopyDataType
  use SurfSoilDataType
  use SoilBGCDataType
  use PlantMngmtDataType
  use EcoSimSumDataType
  use RootDataType
  use EcosysBGCFluxType
  use SoilPropertyDataType
  use IrrigationDataType
  use SedimentDataType
  use GridDataType
  implicit none

  private

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
  character(len=*), private, parameter :: mod_filename = __FILE__

  public :: starts
  contains

  SUBROUTINE starts(NHW,NHE,NVN,NVS)
  !
  !     THIS SUBROUTINE INITIALIZES ALL SOIL VARIABLES
  !
  use InitSOMBGC, only : InitMicbStoichiometry
  use InitVegBGC, only : InitIrradianceGeometry
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX,L,NGL
  REAL(R8) :: ALTY
  real(r8) :: ALTZG
  real(r8) :: CDPTHG
  real(r8) :: YSIN(JSA),YCOS(JSA),YAZI(JSA)
! begin_execution


  !  Initialize controlling parameters
  call InitControlParameters
  !
  !     IRRADIANCE INTERCEPTION GEOMETRY
  call InitIrradianceGeometry(YSIN,YCOS,YAZI)
  !
  !     INITIALIZE C-N AND C-P RATIOS OF RESIDUE AND SOIL
  !
  call InitMicbStoichiometry
  !
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
  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS

      !
      !     MINIMUM SURFACE ELEVATION IN LANDSCAPE
      !
      !     ALT=surface elevation relative to maximum
      !     ALTZG=minimum surface elevation in landscape
      !
      ALT(NY,NX)=ALT(NY,NX)-ALTY
      IF(NX.EQ.NHW.AND.NY.EQ.NVN)THEN
        ALTZG=ALT(NY,NX)
      ELSE
        ALTZG=MIN(ALTZG,ALT(NY,NX))
      ENDIF
      CDPTHG=AMAX1(CDPTHG,CDPTH(NU(NY,NX),NY,NX))
!
!     INITIALIZE ATMOSPHERE VARIABLES
!
!     C*E=atmospheric concentration (g m-3)
!     *E=atmospheric concentration from readi.f (umol mol-1)
!     CO2=CO2,CH4=CH4,OXY=O2,Z2G=N2,Z2O=N2O,NH3=NH3,H2G=H2
!     ATKA=mean annual air temperature (K)
!
      CCO2EI(NY,NX)=CO2EI(NY,NX)*5.36E-04_r8*Tref/ATKA(NY,NX)
      CCO2E(NY,NX)=CO2E(NY,NX)*5.36E-04_r8*Tref/ATKA(NY,NX)
      CCH4E(NY,NX)=CH4E(NY,NX)*5.36E-04_r8*Tref/ATKA(NY,NX)
      COXYE(NY,NX)=OXYE(NY,NX)*1.43E-03_r8*Tref/ATKA(NY,NX)
      CZ2GE(NY,NX)=Z2GE(NY,NX)*1.25E-03_r8*Tref/ATKA(NY,NX)
      CZ2OE(NY,NX)=Z2OE(NY,NX)*1.25E-03_r8*Tref/ATKA(NY,NX)
      CNH3E(NY,NX)=ZNH3E(NY,NX)*6.25E-04_r8*Tref/ATKA(NY,NX)
      CH2GE(NY,NX)=H2GE(NY,NX)*8.92E-05_r8*Tref/ATKA(NY,NX)
!
!     MICROBIAL THERMAL ADAPTATION
!
!     OFFSET=shift in Arrhenius curve used in nitro.f (oC)
!     ATCS=mean annual soil temperature (OC)
!
      OFFSET(NY,NX)=0.333_r8*(12.5_r8-AMAX1(0.0_r8,AMIN1(25.0_r8,ATCS(NY,NX))))
      !     WRITE(*,2222)'OFFSET',OFFSET(NY,NX),ATCS(NY,NX)
!2222  FORMAT(A8,2E12.4)
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
!
!     INITIALIZE SNOWPACK LAYERS
      call InitSnowLayers(NY,NX)
!
!     SURFACE WATER STORAGE AND LOWER HEAT SINK
!
!     VHCPWX,VHCPRX,VHCPNX=minimum heat capacities for solving
!      snowpack,surface litter,soil layer water and heat fluxes
!      DPTHSK=depth at which soil heat sink-source calculated
!     TCNDG=assumed thermal conductivity below lower soil boundary
!     (MJ m-1 K-1 h-1)
!     TKSD=deep source/sink temperature from geothermal flux(K)
!
      VHCPWX(NY,NX)=VHCPWMin*AREA(3,NU(NY,NX),NY,NX)
      VHCPRX(NY,NX)=VHCPRMin*AREA(3,NU(NY,NX),NY,NX)
      VHCPNX(NY,NX)=VHCPNMin*AREA(3,NU(NY,NX),NY,NX)
      DPTHSK(NY,NX)=AMAX1(10.0_r8,CDPTH(NL(NY,NX),NY,NX)+1.0_r8)
      TCS(0,NY,NX)=ATCS(NY,NX)
      TKS(0,NY,NX)=ATKS(NY,NX)
      TKSD(NY,NX)=ATKS(NY,NX)+2.052E-04*DPTHSK(NY,NX)/TCNDG
!

9990  CONTINUE
9995  CONTINUE

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
  DO 9895 NX=NHW,NHE
    DO 9890 NY=NVN,NVS
      ALTZ(NY,NX)=ALTZG
      IF(BKDS(NU(NY,NX),NY,NX).GT.0.0_r8)THEN
        DTBLZ(NY,NX)=DTBLI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX)) &
          *(1.0_r8-DTBLG(NY,NX))
        DTBLD(NY,NX)=AMAX1(0.0_r8,DTBLDI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX)) &
          *(1.0_r8-DTBLG(NY,NX)))
      ELSE
        DTBLZ(NY,NX)=0.0_r8
        DTBLD(NY,NX)=0.0_r8
      ENDIF
      DPTHT(NY,NX)=DTBLZ(NY,NX)
      DO 4400 L=1,NL(NY,NX)
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
          DIST(N,N6,N5,N4)=0.5*(DLYR(N,N3,N2,N1)+DLYR(N,N6,N5,N4))
          XDPTH(N,N6,N5,N4)=AREA(N,N3,N2,N1)/DIST(N,N6,N5,N4)
!1.07 is a scaling parameter for dispersion calculation, reference?
          DISP(N,N6,N5,N4)=0.20*DIST(N,N6,N5,N4)**1.07_r8
        ENDDO
        IF(L.EQ.NU(NY,NX))THEN
          DIST(3,N3,N2,N1)=0.5*DLYR(3,N3,N2,N1)
          XDPTH(3,N3,N2,N1)=AREA(3,N3,N2,N1)/DIST(3,N3,N2,N1)
          DISP(3,N3,N2,N1)=0.20*DIST(3,N3,N2,N1)**1.07_r8
        ENDIF
4400  CONTINUE

      !
      !     ALLOCATE LITTER,SOC TO WOODY,NON-WOODY,MANURE,POC AND HUMUS
      !
      call InitSoilProfile(NY,NX,CDPTHG)
      !
      !     SURFACE LITTER HEAT CAPACITY
      !
      BKVLNM(NY,NX)=AMAX1(0.0,SAND(NU(NY,NX),NY,NX) &
        +SILT(NU(NY,NX),NY,NX)+CLAY(NU(NY,NX),NY,NX))
      VHCP(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VOLW(0,NY,NX) &
        +cpi*VOLI(0,NY,NX)
      VHCM(0,NY,NX)=0.0_r8
      VOLAI(0,NY,NX)=0.0_r8
9890  CONTINUE
9895  CONTINUE
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
  DO 1190 L=NU(NY,NX),NL(NY,NX)
    !     CORGCZ=CORGC(L,NY,NX)
    !     CORGRZ=CORGR(L,NY,NX)
    !     CORGNZ=CORGN(L,NY,NX)
    !     CORGPZ=CORGP(L,NY,NX)
    !
    !     ALLOCATE SOC TO POC(3) AND HUMUS(4)
    !
    !     CORGCX(3)=CORGRZ
    !     CORGCX(4)=AMAX1(0.0,CORGCZ-CORGCX(3))
    !     CORGNX(3)=AMIN1(CNRH(3)*CORGCX(3),CORGNZ)
    !     CORGNX(4)=AMAX1(0.0,CORGNZ-CORGNX(3))
    !     CORGPX(3)=AMIN1(CPRH(3)*CORGCX(3),CORGPZ)
    !     CORGPX(4)=AMAX1(0.0,CORGPZ-CORGPX(3))
    CORGL=AMAX1(0.0,CORGC(L,NY,NX)-CORGR(L,NY,NX))
    TORGL(L)=TORGC+CORGL*BKVL(L,NY,NX)/AREA(3,L,NY,NX)*0.5
    TORGC=TORGC+CORGL*BKVL(L,NY,NX)/AREA(3,L,NY,NX)
1190  CONTINUE
!
!     PARAMETERS TO ALLOCATE HUMUS TO LESS OR MORE RECALCITRANT FRACTIONS
!
!     TORGL=accumulated humus down to soil layer (g m-2)
!     TORGM=TORGL used to calculate allocation (g m-2)
!     HCX=shape parameter for depth effect on allocation
!
  TORGM=AMAX1(2.0E+03,AMIN1(5.0E+03,0.25*TORGL(NJ(NY,NX))))
  IF(TORGM.GT.ZERO)THEN
    HCX=LOG(0.5)/TORGM
  ELSE
    HCX=0.0_r8
  ENDIF

  DO 1200 L=0,NL(NY,NX)
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
        PTDS=1.0E-06_r8*(1.30*CORGCM+2.66_r8*(1.0E+06_r8-CORGCM))
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
        VORGC=CORGCM*1.0E-06*BKDS(L,NY,NX)/PTDS
        VMINL=(CSILT(L,NY,NX)+CCLAY(L,NY,NX))*BKDS(L,NY,NX)/PTDS
        VSAND=CSAND(L,NY,NX)*BKDS(L,NY,NX)/PTDS
        VHCM(L,NY,NX)=((2.496*VORGC+2.385*VMINL+2.128*VSAND) &
          *FMPR(L,NY,NX)+2.128*ROCK(L,NY,NX))*VOLT(L,NY,NX)
      ELSE
        VHCM(L,NY,NX)=0.0_r8
      ENDIF
!
      !     INITIAL SOIL WATER AND ICE CONTENTS
!
      IF(ISOIL(1,L,NY,NX).EQ.0.AND.ISOIL(2,L,NY,NX).EQ.0)THEN
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
        IF(THI(L,NY,NX).GT.1.0)THEN
          THETI(L,NY,NX)=AMAX1(0.0,AMIN1(POROS(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
        ELSEIF(test_aeqb(THI(L,NY,NX),1.0_r8))THEN
          THETI(L,NY,NX)=AMAX1(0.0,AMIN1(FC(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
        ELSEIF(test_aeqb(THI(L,NY,NX),0.0_r8))THEN
          THETI(L,NY,NX)=AMAX1(0.0,AMIN1(WP(L,NY,NX),POROS(L,NY,NX)-THW(L,NY,NX)))
        ELSEIF(THI(L,NY,NX).LT.0.0)THEN
          THETI(L,NY,NX)=0.0_r8
        ELSE
          THETI(L,NY,NX)=THI(L,NY,NX)
        ENDIF
        VOLW(L,NY,NX)=THETW(L,NY,NX)*VOLX(L,NY,NX)
        VOLWX(L,NY,NX)=VOLW(L,NY,NX)
        VOLWH(L,NY,NX)=THETW(L,NY,NX)*VOLAH(L,NY,NX)
        VOLI(L,NY,NX)=THETI(L,NY,NX)*VOLX(L,NY,NX)
        VOLIH(L,NY,NX)=THETI(L,NY,NX)*VOLAH(L,NY,NX)
        VOLP(L,NY,NX)=AMAX1(0.0,VOLA(L,NY,NX)-VOLW(L,NY,NX) &
          -VOLI(L,NY,NX))+AMAX1(0.0,VOLAH(L,NY,NX)-VOLWH(L,NY,NX) &
          -VOLIH(L,NY,NX))
        VHCP(L,NY,NX)=VHCM(L,NY,NX)+cpw*(VOLW(L,NY,NX) &
          +VOLWH(L,NY,NX))+cpi*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
        THETWZ(L,NY,NX)=THETW(L,NY,NX)
        THETIZ(L,NY,NX)=THETI(L,NY,NX)
          !     WRITE(*,2425)'VOLWS',NX,NY,L
          !    2,VOLW(L,NY,NX),THETW(L,NY,NX),VOLI(L,NY,NX),THETI(L,NY,NX)
          !    3,VOLX(L,NY,NX),POROS(L,NY,NX),TKS(L,NY,NX),VHCP(L,NY,NX)
!2425  FORMAT(A8,3I4,20E12.4)
      ENDIF
    ENDIF
    TKS(L,NY,NX)=ATKS(NY,NX)
    TCS(L,NY,NX)=ATCS(NY,NX)
    !
    !     INITIALIZE SOM VARIABLES
    call InitSOMVars(L,NY,NX,FCX)
    !
1200  CONTINUE
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
  ZNH4FA(0:L2,NY,NX)=0.0_r8
  ZNH3FA(0:L2,NY,NX)=0.0_r8
  ZNHUFA(0:L2,NY,NX)=0.0_r8
  ZNO3FA(0:L2,NY,NX)=0.0_r8

  VLNH4(0:L2,NY,NX)=1.0
  VLNO3(0:L2,NY,NX)=1.0
  VLPO4(0:L2,NY,NX)=1.0
  VLNHB(0:L2,NY,NX)=0.0_r8
  VLNOB(0:L2,NY,NX)=0.0_r8
  VLPOB(0:L2,NY,NX)=0.0_r8
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

  ZNH4FB(1:L2,NY,NX)=0.0_r8
  ZNH3FB(1:L2,NY,NX)=0.0_r8
  ZNHUFB(1:L2,NY,NX)=0.0_r8
  ZNO3FB(1:L2,NY,NX)=0.0_r8
  WDNHB(1:L2,NY,NX)=0.0_r8
  DPNHB(1:L2,NY,NX)=0.0_r8
  WDNOB(1:L2,NY,NX)=0.0_r8
  DPNOB(1:L2,NY,NX)=0.0_r8
  WDPOB(1:L2,NY,NX)=0.0_r8
  DPPOB(1:L2,NY,NX)=0.0_r8
  COCU(0:jcplx1,1:L2,NY,NX)=0.0_r8
  CONU(0:jcplx1,1:L2,NY,NX)=0.0_r8
  COPU(0:jcplx1,1:L2,NY,NX)=0.0_r8
  COAU(0:jcplx1,1:L2,NY,NX)=0.0_r8

  end subroutine initFertArrays
!------------------------------------------------------------------------------------------
  subroutine InitSnowLayers(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  real(r8) :: DLYRSI
  real(r8) :: VOLSWI
  real(r8), parameter :: CDPTHSI(JS)=(/0.05_r8,0.15_r8,0.30_r8,0.60_r8,1.00_r8/)
  integer :: L
! begin_execution
!
! CDPTHS=depth to bottom
! DENS0=snow density (Mg m-3)
! VOLSS,VOLWS,VOLIS,VOLS=snow,water,ice,total snowpack volume(m3)
! DPTHA=active layer depth (m)
! CDPTHSI=depth to bottom of snowpack layers
! DLYRS=snowpack layer thickness (m)
! VOLSSL,VOLWSL,VOLISL,VOLSL=snow,water,ice,total layer volume(m3)
! DENSS=layer density (Mg m-3)
! TKW,TCW=later temperature K,oC
! VHCPW=layer volumetric heat capacity (MJ m-3 K-1)
!
  CDPTHS(0,NY,NX)=0.0_r8
  DENS0(NY,NX)=0.10
  VOLSS(NY,NX)=DPTHS(NY,NX)*DENS0(NY,NX)*DH(NY,NX)*DV(NY,NX)
  VOLWS(NY,NX)=0.0_r8
  VOLIS(NY,NX)=0.0_r8
  VOLS(NY,NX)=VOLSS(NY,NX)/DENS0(NY,NX)+VOLWS(NY,NX)+VOLIS(NY,NX)
  DPTHA(NY,NX)=9999.0
  VOLSWI=0.0_r8
  DO 9580 L=1,JS
    IF(L.EQ.1)THEN
      DLYRSI=CDPTHSI(L)
      DLYRS(L,NY,NX)=AMIN1(DLYRSI,DPTHS(NY,NX))
    ELSE
      DLYRSI=CDPTHSI(L)-CDPTHSI(L-1)
      DLYRS(L,NY,NX)=AMIN1(DLYRSI,AMAX1(0.0,DPTHS(NY,NX)-CDPTHSI(L-1)))
    ENDIF
    VOLSSL(L,NY,NX)=DLYRS(L,NY,NX)*DENS0(NY,NX)*DH(NY,NX)*DV(NY,NX)
    VOLWSL(L,NY,NX)=0.0_r8
    VOLISL(L,NY,NX)=0.0_r8
    IF(L.EQ.1)THEN
      VOLSWI=VOLSWI+0.5*(VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX)+VOLISL(L,NY,NX)*DENSI)
    ELSE
      VOLSWI=VOLSWI+0.5*(VOLSSL(L-1,NY,NX)+VOLWSL(L-1,NY,NX) &
        +VOLISL(L-1,NY,NX)*DENSI+VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX) &
        +VOLISL(L,NY,NX)*DENSI)
    ENDIF
    DENSS(L,NY,NX)=DENS0(NY,NX)
    VOLSL(L,NY,NX)=VOLSSL(L,NY,NX)/DENSS(L,NY,NX)+VOLWSL(L,NY,NX)+VOLISL(L,NY,NX)
    VOLSI(L,NY,NX)=DLYRSI*DH(NY,NX)*DV(NY,NX)
    CDPTHS(L,NY,NX)=CDPTHS(L-1,NY,NX)+DLYRS(L,NY,NX)
    TKW(L,NY,NX)=AMIN1(Tref,ATKA(NY,NX))
    TCW(L,NY,NX)=AMIN1(0.0,ATCA(NY,NX))
    VHCPW(L,NY,NX)=cps*VOLSSL(L,NY,NX)+cpw*VOLWSL(L,NY,NX)+cpi*VOLISL(L,NY,NX)
9580  CONTINUE
  end subroutine InitSnowLayers

!------------------------------------------------------------------------------------------
  subroutine InitGridElevation(NHW,NHE,NVN,NVS,YSIN,YCOS,YAZI,ALTY)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: YSIN(JSA),YCOS(JSA),YAZI(JSA)
  REAL(R8),INTENT(OUT):: ALTY
  integer :: NY,NX,N,NN
  REAL(R8) :: DGAZI
  real(r8) :: GSINA(JY,JX),GCOSA(JY,JX)

! begin_execution
! GAZI=ground surface azimuth
! GSIN,GCOS=sine,cosine of ground surface
! OMEGAG=incident sky aNGLe at ground surface
! SLOPE=sine of ground surface slope in (0)aspect, (1)EW,(2)NS directions
! ALT=ground surface elevation
! ALTY=maximum surface elevation in landscape
! IRCHG=runoff boundary flags:0=not possible,1=possible
!
  ALTY=0.0
  DO 9985 NX=NHW,NHE
    DO 9980 NY=NVN,NVS
      ZEROS(NY,NX)=ZERO*DH(NY,NX)*DV(NY,NX)
      ZEROS2(NY,NX)=ZERO2*DH(NY,NX)*DV(NY,NX)
!     compute slopes
      GAZI(NY,NX)=ASP(NY,NX)/RDN
      GSINA(NY,NX)=ABS(SIN(GAZI(NY,NX)))
      GCOSA(NY,NX)=ABS(COS(GAZI(NY,NX)))
      SLOPE(0,NY,NX)=AMAX1(1.745E-04,SIN(SL(NY,NX)/RDN))
      IF(ASP(NY,NX).GE.0.0.AND.ASP(NY,NX).LT.90.0)THEN
        SLOPE(1,NY,NX)=-SLOPE(0,NY,NX)*COS(ASP(NY,NX)/RDN)
        SLOPE(2,NY,NX)=SLOPE(0,NY,NX)*SIN(ASP(NY,NX)/RDN)
        IRCHG(1,1,NY,NX)=1
        IRCHG(2,1,NY,NX)=0
        IRCHG(1,2,NY,NX)=0
        IRCHG(2,2,NY,NX)=1
      ELSEIF(ASP(NY,NX).GE.90.0.AND.ASP(NY,NX).LT.180.0)THEN
        SLOPE(1,NY,NX)=SLOPE(0,NY,NX)*SIN((ASP(NY,NX)-90.0)/RDN)
        SLOPE(2,NY,NX)=SLOPE(0,NY,NX)*COS((ASP(NY,NX)-90.0)/RDN)
        IRCHG(1,1,NY,NX)=0
        IRCHG(2,1,NY,NX)=1
        IRCHG(1,2,NY,NX)=0
        IRCHG(2,2,NY,NX)=1
      ELSEIF(ASP(NY,NX).GE.180.0.AND.ASP(NY,NX).LT.270.0)THEN
        SLOPE(1,NY,NX)=SLOPE(0,NY,NX)*COS((ASP(NY,NX)-180.0)/RDN)
        SLOPE(2,NY,NX)=-SLOPE(0,NY,NX)*SIN((ASP(NY,NX)-180.0)/RDN)
        IRCHG(1,1,NY,NX)=0
        IRCHG(2,1,NY,NX)=1
        IRCHG(1,2,NY,NX)=1
        IRCHG(2,2,NY,NX)=0
      ELSEIF(ASP(NY,NX).GE.270.0.AND.ASP(NY,NX).LE.360.0)THEN
        SLOPE(1,NY,NX)=-SLOPE(0,NY,NX)*SIN((ASP(NY,NX)-270.0)/RDN)
        SLOPE(2,NY,NX)=-SLOPE(0,NY,NX)*COS((ASP(NY,NX)-270.0)/RDN)
        IRCHG(1,1,NY,NX)=1
        IRCHG(2,1,NY,NX)=0
        IRCHG(1,2,NY,NX)=1
        IRCHG(2,2,NY,NX)=0
      ENDIF
      SLOPE(3,NY,NX)=-1.0
      IF(test_aneb(SLOPE(1,NY,NX),0.0_r8).OR.test_aneb(SLOPE(2,NY,NX),0.0_r8))THEN
        FSLOPE(1,NY,NX)=ABS(SLOPE(1,NY,NX))/(ABS(SLOPE(1,NY,NX))+ABS(SLOPE(2,NY,NX)))
        FSLOPE(2,NY,NX)=ABS(SLOPE(2,NY,NX))/(ABS(SLOPE(1,NY,NX))+ABS(SLOPE(2,NY,NX)))
      ELSE
        FSLOPE(1,NY,NX)=0.5
        FSLOPE(2,NY,NX)=0.5
      ENDIF
!    compute incident sky aNGLe at ground surface
      GSIN(NY,NX)=SLOPE(0,NY,NX)
      GCOS(NY,NX)=SQRT(1.0-GSIN(NY,NX)**2)
      DO 240 N=1,JSA
        DGAZI=COS(GAZI(NY,NX)-YAZI(N))
        OMEGAG(N,NY,NX)=AMAX1(0.0,AMIN1(1.0,GCOS(NY,NX)*YSIN(N)+GSIN(NY,NX)*YCOS(N)*DGAZI))
240   CONTINUE
!     compute ground surface elevation
      IF(NX.EQ.NHW)THEN
        IF(NY.EQ.NVN)THEN
          ALT(NY,NX)=0.5*DH(NY,NX)*SLOPE(1,NY,NX)+0.5*DV(NY,NX)*SLOPE(2,NY,NX)
        ELSE
          ALT(NY,NX)=ALT(NY-1,NX) &
            +1.0*DH(NY,NX)*SLOPE(1,NY,NX) &
            +0.5*DV(NY-1,NX)*(SLOPE(2,NY-1,NX)) &
            +0.5*DV(NY,NX)*SLOPE(2,NY,NX)
        ENDIF
      ELSE
        IF(NY.EQ.NVN)THEN
          ALT(NY,NX)=ALT(NY,NX-1) &
            +0.5*DH(NY,NX-1)*SLOPE(1,NY,NX-1) &
            +0.5*DH(NY,NX)*SLOPE(1,NY,NX) &
            +0.5*DV(NY,NX-1)*SLOPE(2,NY,NX-1) &
            +0.5*DV(NY,NX)*SLOPE(2,NY,NX)
        ELSE
          ALT(NY,NX)=(ALT(NY,NX-1) &
            +0.5*DH(NY,NX-1)*SLOPE(1,NY,NX-1) &
            +0.5*DH(NY,NX)*SLOPE(1,NY,NX) &
            +ALT(NY-1,NX) &
            +0.5*DV(NY-1,NX)*SLOPE(2,NY-1,NX) &
            +0.5*DV(NY,N)*SLOPE(2,NY,NX))/2.0
        ENDIF
      ENDIF
      IF(NX.EQ.NHW.AND.NY.EQ.NVN)THEN
        ALTY=ALT(NY,NX)
      ELSE
        ALTY=MAX(ALTY,ALT(NY,NX))
      ENDIF
      WRITE(*,1111)'ALT',NX,NY,((IRCHG(NN,N,NY,NX),NN=1,2),N=1,2) &
        ,ALT(NY,NX),DH(NY,NX),DV(NY,NX),ASP(NY,NX),SL(NY,NX) &
        ,SLOPE(0,NY,NX),SLOPE(1,NY,NX),SLOPE(2,NY,NX) &
        ,GSIN(NY,NX),GCOSA(NY,NX),GSINA(NY,NX)
1111  FORMAT(A8,6I4,20E12.4)
9980  CONTINUE
9985  CONTINUE
  end subroutine InitGridElevation
!------------------------------------------------------------------------------------------
  subroutine InitControlParameters
  implicit none
  !     begin_execution
  real(r8) :: XNPV
  !
  !     NPH=no. of cycles h-1 for water, heat and solute flux calculns
  !     NPT=number of cycles NPH-1 for gas flux calculations
  !     NPG=number of cycles h-1 for gas flux calculations
  !     NPR,NPS=number of cycles NPH-1 for litter,snowpack flux calculns
  !     THETX=minimum air-filled porosity for gas flux calculations
  !     THETPI,DENSI=ice porosity,density
  !
  BKRS=(/0.0333_r8,0.0167_r8,0.0167_r8/)

  call InitSOMConsts

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
  XNPD=600.0*XNPG
  XNPX=AMIN1(1.0_r8,20.0_r8*XNPH)
  XNPA=XNPX*XNPS
  XNPB=XNPX*XNPR
  XNPC=XNPX*XNPV
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
  end subroutine InitControlParameters
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
  OXYS(0,:,:)=0.0_r8
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
  XHVSTC(:,:)=0.0_r8
  XHVSTN(:,:)=0.0_r8
  XHVSTP(:,:)=0.0_r8
  ENGYP(:,:)=0.0_r8
  end subroutine InitAccumulators
!------------------------------------------------------------------------------------------
  subroutine InitLayerDepths(NY,NX)

  implicit none
  integer, intent(in) :: NY, NX
  integer :: L

!     begin_execution
  DO 1195 L=0,NL(NY,NX)
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
      TAREA=TAREA+AREA(3,L,NY,NX)
      CDPTHZ(L,NY,NX)=0.0
      ORGC(L,NY,NX)=(RSC(0,L,NY,NX)+RSC(1,L,NY,NX)+RSC(2,L,NY,NX))*AREA(3,L,NY,NX)
      ORGCX(L,NY,NX)=ORGC(L,NY,NX)
      VOLR(NY,NX)=(RSC(0,L,NY,NX)*1.0E-06/BKRS(0) &
        +RSC(1,L,NY,NX)*1.0E-06/BKRS(1)+RSC(2,L,NY,NX)*1.0E-06/BKRS(2)) &
        *AREA(3,L,NY,NX)
      VOLT(L,NY,NX)=VOLR(NY,NX)
      VOLX(L,NY,NX)=VOLT(L,NY,NX)
      VOLY(L,NY,NX)=VOLX(L,NY,NX)
      VOLTI(L,NY,NX)=VOLT(L,NY,NX)
      BKVL(L,NY,NX)=1.82E-06*ORGC(L,NY,NX)
      DLYRI(3,L,NY,NX)=VOLX(L,NY,NX)/AREA(3,L,NY,NX)
      DLYR(3,L,NY,NX)=DLYRI(3,L,NY,NX)
    ELSE
      IF(BKDSI(L,NY,NX).LE.ZERO)FHOL(L,NY,NX)=0.0
      DLYRI(3,L,NY,NX)=(CDPTH(L,NY,NX)-CDPTH(L-1,NY,NX))
      call check_bool(DLYRI(3,L,NY,NX)<0._r8,'negative soil layer thickness',&
        __LINE__,mod_filename)
!      if(abs(DLYRI(3,L,NY,NX))<1.e-10_r8)then
!      print*,'L,NX,NY=',L,NX,NY
!      print*,'DLYRI(3,L,NY,NX)=',DLYRI(3,L,NY,NX)
!      call print_info('CDPTH(L,NY,NX)',(/padr('CDPTH(L)',12),
!      2padr('CDPTH(L-1)',12)/),
!      3(/CDPTH(L,NY,NX),CDPTH(L-1,NY,NX)/))
!      endif
      DLYR(3,L,NY,NX)=DLYRI(3,L,NY,NX)
      DPTH(L,NY,NX)=0.5*(CDPTH(L,NY,NX)+CDPTH(L-1,NY,NX))
      CDPTHZ(L,NY,NX)=CDPTH(L,NY,NX)-CDPTH(NU(NY,NX),NY,NX)+DLYR(3,NU(NY,NX),NY,NX)
      DPTHZ(L,NY,NX)=0.5*(CDPTHZ(L,NY,NX)+CDPTHZ(L-1,NY,NX))
      VOLT(L,NY,NX)=amax1(AREA(3,L,NY,NX)*DLYR(3,L,NY,NX),1.e-8_r8)
!      if(abs(DLYR(3,L,NY,NX))<1.e-10_r8)then
!      print*,'L,NX,NY=',L,NX,NY
!      print*,'VOLT(L,NY,NX)=',VOLT(L,NY,NX)
!      call print_info('DLYR(3,L,NY,NX)==0.',(/padr('DLYR(3)',10)/),
!     2(/DLYR(3,L,NY,NX)/))
!      endif
      VOLX(L,NY,NX)=VOLT(L,NY,NX)*FMPR(L,NY,NX)
      VOLY(L,NY,NX)=VOLX(L,NY,NX)
      VOLTI(L,NY,NX)=VOLT(L,NY,NX)
      BKVL(L,NY,NX)=BKDS(L,NY,NX)*VOLX(L,NY,NX)
      RTDNT(L,NY,NX)=0.0
    ENDIF
    AREA(1,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(2,L,NY,NX)
    AREA(2,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(1,L,NY,NX)
1195  CONTINUE
  CDPTH(0,NY,NX)=CDPTH(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX)
  CDPTHI(NY,NX)=CDPTH(0,NY,NX)
  AREA(3,NL(NY,NX)+1:JZ,NY,NX)=DLYR(1,NL(NY,NX),NY,NX)*DLYR(2,NL(NY,NX),NY,NX)
  end subroutine InitLayerDepths
end module StartsMod
