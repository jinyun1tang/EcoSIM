module StartsMod
!!
! Description:
! code to initalize soil variables

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : padr, print_info
  use minimathMod, only : test_aeqb, test_aneb
  use EcosimConst
  use MicrobialDataType
  use SOMDataType
  use SoilChemDataType
  use InitSOMBGC
  implicit none

  private
  include "parameters.h"
  include "blkc.h"
  include "blk2a.h"
  include "blk2b.h"
  include "blk2c.h"
  include "blk5.h"
  include "blk8a.h"
  include "blk8b.h"
  include "blk11a.h"
  include "blk11b.h"
  include "blk13a.h"
  include "blk13b.h"
  include "blk13c.h"
  include "blk16.h"
  include "blk18a.h"
  include "blk18b.h"

  real(r8) :: ALTY,ALTZG,CDPTHG
  real(r8) :: DAZI,DGAZI,DLYRSI
  real(r8) :: OMEGY,OMEGZ
  real(r8) :: PTDS,TORGM
  real(r8) :: VOLSWI,VORGC,VMINL,VSAND,YAGL
  real(r8) :: ZAGL

  real(r8) :: YSIN(4),YCOS(4),YAZI(4),ZAZI(4) &
    ,GSINA(JY,JX),GCOSA(JY,JX),ALTX(JV,JH),CDPTHSI(JS) &
    ,TORGL(JZ)
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
  real(r8), PARAMETER :: PSIPS=-0.5E-03_r8,RDN=57.29577951_r8

  public :: starts
  contains

  SUBROUTINE starts(NHW,NHE,NVN,NVS)
  !
  !     THIS SUBROUTINE INITIALIZES ALL SOIL VARIABLES
  !
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX,L,NGL
  !     begin_execution


  !     Initialize controling parameters
  call InitControlParameters
  !
  !     IRRADIANCE INTERCEPTION GEOMETRY
  call InitIrradianceGeometry
  !
  !     INITIALIZE C-N AND C-P RATIOS OF RESIDUE AND SOIL
  !
  call InitCNPRatios
  !
  !     CALCULATE ELEVATION OF EACH GRID CELL
  !
  call InitGridElevation(NHW,NHE,NVN,NVS)
  !
  !     INITIALIZE ACCUMULATORS AND MASS BALANCE CHECKS
  !     OF EACH GRID CELL
  !
  ALTZG=0.0_r8
  CDPTHG=0.0_r8
  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS

      call InitAccumulators(NY,NX)
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
      TCNDG=8.1E-03_r8
      TKSD(NY,NX)=ATKS(NY,NX)+2.052E-04*DPTHSK(NY,NX)/TCNDG
      !
      !     INITIALIZE COMMUNITY CANOPY
!
      ZT(NY,NX)=0.0_r8
      ZL(0,NY,NX)=0.0_r8
      DO 1925 L=1,JC
        ZL(L,NY,NX)=0.0_r8
        ARLFT(L,NY,NX)=0.0_r8
        ARSTT(L,NY,NX)=0.0_r8
        WGLFT(L,NY,NX)=0.0_r8
1925  CONTINUE
!
!     INITIALIZE SEDIMENT LOAD IN EROSION MODEL
!
      IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
        SED(NY,NX)=0.0_r8
      ENDIF
9990  CONTINUE
9995  CONTINUE
!
!     INITIALIZE GRID CELL DIMENSIONS
!
  call InitSoilVars(NHW,NHE,NVN,NVS)

  RETURN
  END subroutine starts
!------------------------------------------------------------------------------------------
  subroutine InitSoilVars(NHW,NHE,NVN,NVS)
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

  integer :: NY,NX,L,N
  integer :: N1,N2,N3,N4,N5,N6


  ! begin_execution
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
        DO 4320 N=NCN(N2,N1),3
          IF(N.EQ.1)THEN
            IF(NX.EQ.NHE)THEN
              GO TO 4320
            ELSE
              N4=NX+1
              N5=NY
              N6=L
            ENDIF
          ELSEIF(N.EQ.2)THEN
            IF(NY.EQ.NVS)THEN
              GO TO 4320
            ELSE
              N4=NX
              N5=NY+1
              N6=L
            ENDIF
          ELSEIF(N.EQ.3)THEN
            IF(L.EQ.NL(NY,NX))THEN
              GO TO 4320
            ELSE
              N4=NX
              N5=NY
              N6=L+1
            ENDIF
          ENDIF
          DIST(N,N6,N5,N4)=0.5*(DLYR(N,N3,N2,N1)+DLYR(N,N6,N5,N4))
          XDPTH(N,N6,N5,N4)=AREA(N,N3,N2,N1)/DIST(N,N6,N5,N4)
          DISP(N,N6,N5,N4)=0.20*DIST(N,N6,N5,N4)**1.07
4320    CONTINUE
        IF(L.EQ.NU(NY,NX))THEN
          DIST(3,N3,N2,N1)=0.5*DLYR(3,N3,N2,N1)
          XDPTH(3,N3,N2,N1)=AREA(3,N3,N2,N1)/DIST(3,N3,N2,N1)
          DISP(3,N3,N2,N1)=0.20*DIST(3,N3,N2,N1)**1.07
        ENDIF
4400  CONTINUE

      !
      !     ALLOCATE LITTER,SOC TO WOODY,NON-WOODY,MANURE,POC AND HUMUS
      !
      call InitSoilProfile(NY,NX)
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
  subroutine InitIrradianceGeometry

  implicit none
  integer :: L,M,N
  !     begin_execution
  !     ZSIN,ZCOS=sine,cosine of leaf inclination class
  !     ZAZI=leaf azimuth class
  !     YAZI,YSIN,YCOS=sky azimuth,sine,cosine of sky azimuth
  !     OMEGA,OMEGX=incident aNGLe of diffuse radn at leaf,horizontal surface
  !     IALBY:1=backscattering,2=forward scattering of sky radiation
  !

  ZSIN(1)=0.195
  ZSIN(2)=0.556
  ZSIN(3)=0.831
  ZSIN(4)=0.981
  ZCOS(1)=0.981
  ZCOS(2)=0.831
  ZCOS(3)=0.556
  ZCOS(4)=0.195
  DO 205 L=1,4
    ZAZI(L)=(L-0.5)*3.1416/4.0
205   CONTINUE
  DO 230 N=1,4
    YAZI(N)=3.1416*(2*N-1)/4.0
    YAGL=3.1416/4.0
    YSIN(N)=SIN(YAGL)
    YCOS(N)=COS(YAGL)
    TYSIN=TYSIN+YSIN(N)
    DO 225 L=1,4
      DAZI=COS(ZAZI(L)-YAZI(N))
      DO  M=1,4
        OMEGY=ZCOS(M)*YSIN(N)+ZSIN(M)*YCOS(N)*DAZI
        OMEGA(N,M,L)=ABS(OMEGY)
        OMEGX(N,M,L)=OMEGA(N,M,L)/YSIN(N)
        IF(ZCOS(M).GT.YSIN(N))THEN
          OMEGZ=ACOS(OMEGY)
        ELSE
          OMEGZ=-ACOS(OMEGY)
        ENDIF
        IF(OMEGZ.GT.-1.5708)THEN
          ZAGL=YAGL+2.0*OMEGZ
        ELSE
          ZAGL=YAGL-2.0*(3.1416+OMEGZ)
        ENDIF
        IF(ZAGL.GT.0.0.AND.ZAGL.LT.3.1416)THEN
          IALBY(N,M,L)=1
        ELSE
          IALBY(N,M,L)=2
        ENDIF
      ENDDO
225 CONTINUE
230 CONTINUE
  end subroutine InitIrradianceGeometry
!------------------------------------------------------------------------------------------
  subroutine InitSoilProfile(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  integer  :: L,M,K,N,KK,NN,NGL
  real(r8) :: CORGCM,HCX,TORGC
  real(r8) :: CORGL,TORGLL
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

    call InitSOMProfile(L,NY,NX,HCX,CDPTHG,TORGLL,CORGCM)
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
    call InitSOMVars(L,NY,NX)

    !
    !     INITIALIZE FERTILIZER ARRAYS
    !
    ZNH4FA(L,NY,NX)=0.0_r8
    ZNH3FA(L,NY,NX)=0.0_r8
    ZNHUFA(L,NY,NX)=0.0_r8
    ZNO3FA(L,NY,NX)=0.0_r8
    IF(L.GT.0)THEN
      ZNH4FB(L,NY,NX)=0.0_r8
      ZNH3FB(L,NY,NX)=0.0_r8
      ZNHUFB(L,NY,NX)=0.0_r8
      ZNO3FB(L,NY,NX)=0.0_r8
      WDNHB(L,NY,NX)=0.0_r8
      DPNHB(L,NY,NX)=0.0_r8
      WDNOB(L,NY,NX)=0.0_r8
      DPNOB(L,NY,NX)=0.0_r8
      WDPOB(L,NY,NX)=0.0_r8
      DPPOB(L,NY,NX)=0.0_r8
    ENDIF
    VLNH4(L,NY,NX)=1.0
    VLNO3(L,NY,NX)=1.0
    VLPO4(L,NY,NX)=1.0
    VLNHB(L,NY,NX)=0.0_r8
    VLNOB(L,NY,NX)=0.0_r8
    VLPOB(L,NY,NX)=0.0_r8
    ROXYX(L,NY,NX)=0.0_r8
    RNH4X(L,NY,NX)=0.0_r8
    RNO3X(L,NY,NX)=0.0_r8
    RNO2X(L,NY,NX)=0.0_r8
    RN2OX(L,NY,NX)=0.0_r8
    RPO4X(L,NY,NX)=0.0_r8
    RP14X(L,NY,NX)=0.0_r8
    RVMXC(L,NY,NX)=0.0_r8
    RNHBX(L,NY,NX)=0.0_r8
    RN3BX(L,NY,NX)=0.0_r8
    RN2BX(L,NY,NX)=0.0_r8
    RPOBX(L,NY,NX)=0.0_r8
    RP1BX(L,NY,NX)=0.0_r8
    RVMBC(L,NY,NX)=0.0_r8
    DO 1250 K=0,4
      IF(L.GT.0)THEN
        COCU(K,L,NY,NX)=0.0_r8
        CONU(K,L,NY,NX)=0.0_r8
        COPU(K,L,NY,NX)=0.0_r8
        COAU(K,L,NY,NX)=0.0_r8
      ENDIF
1250  CONTINUE
    ZNHUI(L,NY,NX)=0.0_r8
    ZNHU0(L,NY,NX)=0.0_r8
    ZNFNI(L,NY,NX)=0.0_r8
    ZNFN0(L,NY,NX)=0.0_r8
1200  CONTINUE
  end subroutine InitSoilProfile

!------------------------------------------------------------------------------------------
  subroutine InitSnowLayers(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX
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
  subroutine InitCNPRatios

  integer :: K,N,NGL
! begin_execution
! CNOFC,CPOFC=fractions to allocate N,P to kinetic components
! CNOMC,CPOMC=maximum N:C and P:C ratios in microbial biomass

  CNOFC(1,0)=0.005_r8
  CNOFC(2,0)=0.005_r8
  CNOFC(3,0)=0.005_r8
  CNOFC(4,0)=0.020_r8
  CPOFC(1,0)=0.0005_r8
  CPOFC(2,0)=0.0005_r8
  CPOFC(3,0)=0.0005_r8
  CPOFC(4,0)=0.0020_r8
  CNOFC(1,1)=0.020_r8
  CNOFC(2,1)=0.020_r8
  CNOFC(3,1)=0.020_r8
  CNOFC(4,1)=0.020_r8
  CPOFC(1,1)=0.0020_r8
  CPOFC(2,1)=0.0020_r8
  CPOFC(3,1)=0.0020_r8
  CPOFC(4,1)=0.0020_r8
  CNOFC(1,2)=0.020_r8
  CNOFC(2,2)=0.020_r8
  CNOFC(3,2)=0.020_r8
  CNOFC(4,2)=0.020_r8
  CPOFC(1,2)=0.0020_r8
  CPOFC(2,2)=0.0020_r8
  CPOFC(3,2)=0.0020_r8
  CPOFC(4,2)=0.0020_r8
  FL(1)=0.55_r8
  FL(2)=0.45_r8
  DO 95 K=0,5
    DO  N=1,7
      IF(K.LE.4.AND.N.EQ.3)THEN
        DO NGL=1,JG
          CNOMC(1,NGL,N,K)=0.15
          CNOMC(2,NGL,N,K)=0.09
          CPOMC(1,NGL,N,K)=0.015
          CPOMC(2,NGL,N,K)=0.009
        ENDDO
      ELSE
        do NGL=1,JG
          CNOMC(1,NGL,N,K)=0.225
          CNOMC(2,NGL,N,K)=0.135
          CPOMC(1,NGL,N,K)=0.0225
          CPOMC(2,NGL,N,K)=0.0135
        enddo
      ENDIF
      do NGL=1,JG
        CNOMC(3,NGL,N,K)=FL(1)*CNOMC(1,NGL,N,K)+FL(2)*CNOMC(2,NGL,N,K)
        CPOMC(3,NGL,N,K)=FL(1)*CPOMC(1,NGL,N,K)+FL(2)*CPOMC(2,NGL,N,K)
      enddo
     enddo
95  CONTINUE
    end subroutine InitCNPRatios
!------------------------------------------------------------------------------------------
  subroutine InitGridElevation(NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX,N,NN
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
      GSIN(NY,NX)=SLOPE(0,NY,NX)
      GCOS(NY,NX)=SQRT(1.0-GSIN(NY,NX)**2)
      DO 240 N=1,4
        DGAZI=COS(GAZI(NY,NX)-YAZI(N))
        OMEGAG(N,NY,NX)=AMAX1(0.0,AMIN1(1.0,GCOS(NY,NX)*YSIN(N)+GSIN(NY,NX)*YCOS(N)*DGAZI))
240   CONTINUE
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
  !
  !     NPH=no. of cycles h-1 for water, heat and solute flux calculns
  !     NPT=number of cycles NPH-1 for gas flux calculations
  !     NPG=number of cycles h-1 for gas flux calculations
  !     NPR,NPS=number of cycles NPH-1 for litter,snowpack flux calculns
  !     THETX=minimum air-filled porosity for gas flux calculations
  !     THETPI,DENSI=ice porosity,density
  !
  BKRS=(/0.0333_r8,0.0167_r8,0.0167_r8/)

  FORGC=0.1E+06_r8
  FVLWB=1.0_r8
  FCH4F=0.01_r8
  PSIHY=-2500.0_r8
  FCI=0.05_r8
  WPI=0.025_r8
  CDPTHSI=(/0.05_r8,0.15_r8,0.30_r8,0.60_r8,1.00_r8/)
  POROQ=0.66_r8

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
  THETX=1.0E-03_r8
  THETPI=0.00_r8
  DENSI=0.92_r8-THETPI
  DENSJ=1.0_r8-DENSI
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
  VAP=2465.0_r8   !kJ/kg
  VAPS=2834.0_r8
  OXKM=0.080_r8
  TYSIN=0.0_r8
  end subroutine InitControlParameters
!------------------------------------------------------------------------------------------
  subroutine InitAccumulators(NY,NX)
  implicit none
  integer, intent(in) :: NY, NX
  integer :: N
!     begin_execution
  DO 600 N=1,12
    TDTPX(NY,NX,N)=0.0
    TDTPN(NY,NX,N)=0.0
    TDRAD(NY,NX,N)=1.0
    TDWND(NY,NX,N)=1.0
    TDHUM(NY,NX,N)=1.0
    TDPRC(NY,NX,N)=1.0
    TDIRI(NY,NX,N)=1.0
    TDCO2(NY,NX,N)=1.0
    TDCN4(NY,NX,N)=1.0
    TDCNO(NY,NX,N)=1.0
600 CONTINUE
  IUTYP(NY,NX)=0
  IFNHB(NY,NX)=0
  IFNOB(NY,NX)=0
  IFPOB(NY,NX)=0
  IFLGS(NY,NX)=1
  IFLGT(NY,NX)=0
  ATCA(NY,NX)=ATCAI(NY,NX)
  ATCS(NY,NX)=ATCAI(NY,NX)
  ATKA(NY,NX)=ATCA(NY,NX)+TC2K
  ATKS(NY,NX)=ATCS(NY,NX)+TC2K
  URAIN(NY,NX)=0.0
  UCO2G(NY,NX)=0.0
  UCH4G(NY,NX)=0.0
  UOXYG(NY,NX)=0.0
  UN2GG(NY,NX)=0.0
  UN2OG(NY,NX)=0.0
  UNH3G(NY,NX)=0.0
  UN2GS(NY,NX)=0.0
  UCO2F(NY,NX)=0.0
  UCH4F(NY,NX)=0.0
  UOXYF(NY,NX)=0.0
  UN2OF(NY,NX)=0.0
  UNH3F(NY,NX)=0.0
  UPO4F(NY,NX)=0.0
  UORGF(NY,NX)=0.0
  UFERTN(NY,NX)=0.0
  UFERTP(NY,NX)=0.0
  UVOLO(NY,NX)=0.0
  UEVAP(NY,NX)=0.0
  URUN(NY,NX)=0.0
  USEDOU(NY,NX)=0.0
  UCOP(NY,NX)=0.0
  UDOCQ(NY,NX)=0.0
  UDOCD(NY,NX)=0.0
  UDONQ(NY,NX)=0.0
  UDOND(NY,NX)=0.0
  UDOPQ(NY,NX)=0.0
  UDOPD(NY,NX)=0.0
  UDICQ(NY,NX)=0.0
  UDICD(NY,NX)=0.0
  UDINQ(NY,NX)=0.0
  UDIND(NY,NX)=0.0
  UDIPQ(NY,NX)=0.0
  UDIPD(NY,NX)=0.0
  UIONOU(NY,NX)=0.0
  UXCSN(NY,NX)=0.0
  UXZSN(NY,NX)=0.0
  UXPSN(NY,NX)=0.0
  UDRAIN(NY,NX)=0.0
  ZDRAIN(NY,NX)=0.0
  PDRAIN(NY,NX)=0.0
  DPNH4(NY,NX)=0.0
  DPNO3(NY,NX)=0.0
  DPPO4(NY,NX)=0.0
  OXYS(0,NY,NX)=0.0
  FRADG(NY,NX)=1.0
  THRMG(NY,NX)=0.0
  THRMC(NY,NX)=0.0
  TRN(NY,NX)=0.0
  TLE(NY,NX)=0.0
  TSH(NY,NX)=0.0
  TGH(NY,NX)=0.0
  TLEC(NY,NX)=0.0
  TSHC(NY,NX)=0.0
  TLEX(NY,NX)=0.0
  TSHX(NY,NX)=0.0
  TCNET(NY,NX)=0.0
  TVOLWC(NY,NX)=0.0
  ARLFC(NY,NX)=0.0
  ARSTC(NY,NX)=0.0
  TFLWC(NY,NX)=0.0
  PPT(NY,NX)=0.0
  DYLN(NY,NX)=12.0
  ALBX(NY,NX)=ALBS(NY,NX)
  XHVSTC(NY,NX)=0.0
  XHVSTN(NY,NX)=0.0
  XHVSTP(NY,NX)=0.0
  ENGYP(NY,NX)=0.0
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
  end subroutine InitLayerDepths
end module StartsMod
