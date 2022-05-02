module CanopyCondMod

  use data_kind_mod, only : r8 => SHR_KIND_R8
  USE EcoSIMCtrlDataType
  use PlantTraitDataType
  use LandSurfDataType
  use SnowDataType
  use EcoSimConst
  use GridConsts
  use ClimForcDataType
  use GridDataType
  use FlagDataType
  use CanopyDataType
  use SoilPhysDataType
  use SoilPropertyDataType
  use CanopyRadDataType
  use SoilWaterDataType
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__

  public :: CanopyConditionModel

  real(r8), parameter :: RAM=2.78E-03_r8    !minimum boundary layer resistance (h m-1)
  real(r8), parameter :: ALBRW=0.1_r8       !stalk albedo for shortwave
  real(r8), parameter :: ALBPW=0.1_r8       !stalk albedo for PAR
  real(r8), parameter :: ABSRW=1.0_r8-ALBRW
  real(r8), parameter :: ABSPW=1.0_r8-ALBPW
  real(r8), parameter :: CFW=0.5_r8         !stalk clumping factor,FORGW=minimum SOC or organic soil (g Mg-1)

  contains
  subroutine CanopyConditionModel(I,J,NY,NX,DPTH0)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: DPTH0(JY,JX) !water+ice depth at surface

  call MultiLayerSurfaceRadiation(I,J,NY,NX,DPTH0)

  call DivideCanopyLayerByLAI(NY,NX)

  call CalcBoundaryLayerProperties(NY,NX,DPTH0)

  end subroutine CanopyConditionModel

!------------------------------------------------------------------------------------------

  subroutine CalcBoundaryLayerProperties(NY,NX,DPTH0)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: DPTH0(JY,JX)
  real(r8) :: ARLSC
  real(r8) :: ARLSG
  real(r8) :: ZX,ZY,ZE
  REAL(R8) :: ZZ
!     begin_execution
!     CANOPY ZERO PLANE AND ROUGHNESS HEIGHTS
!
!     ARLFC,ARSTC=leaf,stalk area of combined canopy
!     DPTHS,DPTH0=snowpack,surface water depths
!     ZT,ZD,ZR=canopy,zero plane displacement,roughness height
!     ZZ=reference height for wind speed
!
  ARLSC=ARLFC(NY,NX)+ARSTC(NY,NX)
  IF(ARLSC.GT.ZEROS(NY,NX) &
    .AND.ZT(NY,NX).GE.DPTHS(NY,NX)-ZERO &
    .AND.ZT(NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
    ARLSG=ARLSC/AREA(3,NU(NY,NX),NY,NX)
    ZX=EXP(-0.5*ARLSG)
    ZY=1.0-ZX
    ZD(NY,NX)=ZT(NY,NX)*AMAX1(0.0,1.0-2.0/ARLSG*ZY)
    ZE=ZT(NY,NX)*AMAX1(0.05,ZX*ZY)
  ELSE
    ZD(NY,NX)=0.0
    ZE=0.0
  ENDIF
  IF(IFLGW.EQ.1)THEN
    ZZ=Z0(NY,NX)+ZT(NY,NX)
  ELSE
    ZZ=AMAX1(Z0(NY,NX),ZD(NY,NX)+2.0)
  ENDIF
  IF(IETYP(NY,NX).GE.0)THEN
    IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
      ZR(NY,NX)=AMAX1(0.001,ZE,ZW)
    ELSE
      ZR(NY,NX)=AMAX1(0.001,ZE,ZS(NY,NX))
    ENDIF
!
!     CANOPY ISOTHERMAL BOUNDARY LAYER RESISTANCE
!
!     RAB,RAM=biome canopy,minimum isothermal boundary layer resistance
!     UA=wind speed
!     RIB=canopy isothermal Richardson number
!
    RAB(NY,NX)=AMAX1(RAM,(LOG((ZZ-ZD(NY,NX))/ZR(NY,NX)))**2/(0.168*UA(NY,NX)))
    RIB(NY,NX)=1.27E+08*(ZZ-ZR(NY,NX))/(UA(NY,NX)**2*TKA(NY,NX))
  ELSE
    RAB(NY,NX)=RAM
    RIB(NY,NX)=0.0
  ENDIF
  end subroutine CalcBoundaryLayerProperties

!------------------------------------------------------------------------------------------

  subroutine DivideCanopyLayerByLAI(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  real(r8) :: ZL1(0:JZ,JY,JX)
  real(r8) :: ART,ARL
  real(r8) :: ARX
  real(r8) :: DZL
  integer :: NZ,L
  !     begin_execution
  !
  !     DIVISION OF CANOPY INTO JC LAYERS WITH EQUAL LAI
  !
  !     ZT,ZC=heights of combined canopy,PFT canopy
  !     ZL=height to bottom of each canopy layer
  !     ARLFC,ARSTC=leaf,stalk area of combined canopy
  !     ARLFT,ARSTT=leaf,stalk area of combined canopy layer
  !
  ZT(NY,NX)=0.0
  DO 9685 NZ=1,NP(NY,NX)
    ZT(NY,NX)=AMAX1(ZT(NY,NX),ZC(NZ,NY,NX))
9685  CONTINUE
  ZL(JC,NY,NX)=ZT(NY,NX)+0.01
  ZL1(JC,NY,NX)=ZL(JC,NY,NX)
  ZL1(0,NY,NX)=0.0
  ART=(ARLFC(NY,NX)+ARSTC(NY,NX))/JC
  IF(ART.GT.ZEROS(NY,NX))THEN
    DO 2765 L=JC,2,-1
      ARL=ARLFT(L,NY,NX)+ARSTT(L,NY,NX)
      IF(ARL.GT.1.01*ART)THEN
        DZL=ZL(L,NY,NX)-ZL(L-1,NY,NX)
        ZL1(L-1,NY,NX)=ZL(L-1,NY,NX)+0.5*AMIN1(1.0,(ARL-ART)/ARL)*DZL
      ELSEIF(ARL.LT.0.99*ART)THEN
        ARX=ARLFT(L-1,NY,NX)+ARSTT(L-1,NY,NX)
        DZL=ZL(L-1,NY,NX)-ZL(L-2,NY,NX)
        IF(ARX.GT.ZEROS(NY,NX))THEN
          ZL1(L-1,NY,NX)=ZL(L-1,NY,NX)-0.5*AMIN1(1.0,(ART-ARL)/ARX)*DZL
        ENDIF
      ELSE
        ZL1(L-1,NY,NX)=ZL(L-1,NY,NX)
      ENDIF
      !     IF(J.EQ.12)THEN
      !     WRITE(*,3233)'ZL',I,J,NX,NY,L,ZL1(L,NY,NX),ZL1(L-1,NY,NX)
      !    3,ZL(L,NY,NX),ZL(L-1,NY,NX),ART,ARL,ARX
      !    2,DZL,ARLFC(NY,NX),ARSTC(NY,NX),ARLFT(L,NY,NX),ARSTT(L,NY,NX)
      !    3,ARLFT(L-1,NY,NX),ARSTT(L-1,NY,NX)
      !    3,ZL(JC,NY,NX),ZT(NY,NX),ZC(1,NY,NX)
      !3233  FORMAT(A8,5I4,30E12.4)
      !     ENDIF
2765  CONTINUE
    DO 2770 L=JC,2,-1
      ZL(L-1,NY,NX)=ZL1(L-1,NY,NX)
!     ZL(L-1,NY,NX)=AMAX1(0.0,AMIN1(ZL(L,NY,NX)-1.0E-06
!    2,ZL(L-1,NY,NX)))
2770  CONTINUE
  ENDIF
  end subroutine DivideCanopyLayerByLAI

!------------------------------------------------------------------------------------------

  subroutine MultiLayerSurfaceRadiation(I,J,NY,NX,DPTH0)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(in) :: DPTH0(JY,JX)
  integer :: NB,NZ,L,K,M,N,NN
  integer :: IALBS(JLI,JSA)
  real(r8) :: TAUY(0:JC+1)
  real(r8) :: PARDIR(JLI,JSA,JP,JY,JX)
  real(r8) :: PARDIW(JLI,JSA,JP,JY,JX)
  real(r8) :: RABSL(0:JC+1)
  real(r8) :: RABPL(0:JC+1)
  real(r8) :: RAFSL(0:JC+1)
  real(r8) :: RAFPL(0:JC+1)
  real(r8) :: RADSL(JP,JY,JX)
  real(r8) :: RADPL(JP,JY,JX)
  real(r8) :: RAYSW(JP,JY,JX)
  real(r8) :: RAYPW(JP,JY,JX)
  real(r8) :: RADSW(JP,JY,JX)
  real(r8) :: RADPW(JP,JY,JX)
  real(r8) :: RADWA(JP,JY,JX)
  real(r8) :: RADSA(JP,JY,JX)
  real(r8) :: RAPSA(JP,JY,JX)
  REAL(R8) :: RAPWA(JP,JY,JX)
  real(r8) :: RADS1(JP,JY,JX),RADS2(JP,JY,JX)
  real(r8) :: RADP1(JP,JY,JX),RADP2(JP,JY,JX)
  real(r8) :: RADQ1(JP,JY,JX),RADQ2(JP,JY,JX)
  real(r8) :: RADW1(JP,JY,JX),RADW2(JP,JY,JX)
  real(r8) :: RAYSL(JP,JY,JX),RAYPL(JP,JY,JX)
  real(r8) :: RAYS1(JP,JY,JX),RAYS2(JP,JY,JX)
  real(r8) :: RAYP1(JP,JY,JX),RAYP2(JP,JY,JX)
  real(r8) :: RAYW1(JP,JY,JX),RAYW2(JP,JY,JX)
  real(r8) :: RAYQ1(JP,JY,JX),RAYQ2(JP,JY,JX)
  real(r8) :: BETA(JLI,JSA)                         !sine of direct solar radiation on leaf surface, [-]
  real(r8) :: BETX(JLI,JSA)                         !sine of direct solar radiation on leaf surface/sine of direct solar radiation, [-]
  REAL(R8) :: RDNDIR(JLI,JSA,JP,JY,JX),RDNDIW(JLI,JSA,JP,JY,JX)
  real(r8) :: TSURF(JLI,JZ,JP,JY,JX),TSURFB(JLI,JZ,JP,JY,JX)
  real(r8) :: TRADC(JY,JX),TRAPC(JY,JX),TRADG(JY,JX),TRAPG(JY,JX)
  real(r8) :: ALBW,ALBG
  real(r8) :: BETAG,BETY,BETZ
  real(r8) :: DGAZI,DAZI
  real(r8) :: FSNOW,FRADPT
  real(r8) :: RADYL,RADYN,RADYW
  real(r8) :: RADST
  real(r8) :: RADWT,RADPT,RADQT
  real(r8) :: RADSG,RADYG
  real(r8) :: RAPYL,RAPYN,RAPYW
  real(r8) :: RAPSG,RAPYG,RAPG
  real(r8) :: STOPS, STOPX, STOPY, STOPZ, STOPSZ, STOPYZ
  real(r8) :: RA1ST,RA1WT,RA1PT,RA1QT,RA2ST
  real(r8) :: RA2WT,RA2PT,RA2QT
  REAL(R8) :: RASG,SAGL
  real(r8) :: THETW1
  real(r8) :: TSURFX,TSURFY,TSURFZ,TSURFS
  real(r8) :: TSURWY,TSURWZ,TSURWS,TSURWX
  real(r8) :: XTAUS,XTAUY,XAREA
  real(r8) :: YAREA
  REAL(R8) :: ZAZI,ZAGL
  !     begin_execution
  !     MULTILAYER CANOPY INTERECEPTION OF DIRECT AND DIFFUSE RADIATION
  !     IN SW AND VISIBLE BANDS BY INCLINATION N, AZIMUTH M, LAYER L,
  !     NODE K, BRANCH NB, PFT NZ
  !
  !     ARLFS,ARLSS=leaf+stalk area of combined,each PFT canopy
  !     ZL=height to bottom of canopy layer
  !     DPTHS,DPTH0=snowpack,surface water depths
  !     ARLFL,ARSTK=leaf,stalk areas of PFT
  !     RAD,RAP=vertical direct+diffuse SW,PAR
  !     RADS,RADY,RAPS,RAPY=solar beam direct,diffuse SW,PAR
  !     SSIN,TYSIN=sine of solar,sky angles
  !     RADC,RADP=total SW,PAR absorbed by canopy
  !     CFX=clumping factor for self-shading
  !
  ARLSS(NY,NX)=0.0
  DO 1135 NZ=1,NP(NY,NX)
    ARLFS(NZ,NY,NX)=0.0
    DO  NB=1,NBR(NZ,NY,NX)
      DO  L=1,JC
        IF(ZL(L-1,NY,NX).GE.DPTHS(NY,NX)-ZERO &
          .AND.ZL(L-1,NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
          DO 1130 K=1,JNODS
            ARLFS(NZ,NY,NX)=ARLFS(NZ,NY,NX)+ARLFL(L,K,NB,NZ,NY,NX)
            ARLSS(NY,NX)=ARLSS(NY,NX)+ARLFL(L,K,NB,NZ,NY,NX)
1130      CONTINUE
          ARLFS(NZ,NY,NX)=ARLFS(NZ,NY,NX)+ARSTK(L,NB,NZ,NY,NX)
          ARLSS(NY,NX)=ARLSS(NY,NX)+ARSTK(L,NB,NZ,NY,NX)
        ENDIF
      enddo
    enddo
1135  CONTINUE
  IF(SSIN(NY,NX).GT.ZERO)THEN
    RAD(NY,NX)=RADS(NY,NX)*SSIN(NY,NX)+RADY(NY,NX)*TYSIN
    RAP(NY,NX)=RAPS(NY,NX)*SSIN(NY,NX)+RAPY(NY,NX)*TYSIN
  ELSE
    RADS(NY,NX)=0.0
    RADY(NY,NX)=0.0
    RAPS(NY,NX)=0.0
    RAPY(NY,NX)=0.0
    RAD(NY,NX)=0.0
    RAP(NY,NX)=0.0
  ENDIF
  TRADC(NY,NX)=0.0
  TRAPC(NY,NX)=0.0
  DO 1025 NZ=1,NP(NY,NX)
    RADC(NZ,NY,NX)=0.0
    RADP(NZ,NY,NX)=0.0
    CFX(NZ,NY,NX)=CF(NZ,NY,NX)*(1.0-0.025 &
      *ARLFP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX))
1025  CONTINUE
  !
  !     ANGLE BETWEEN SUN AND GROUND SURFACE
  !
  !     SAZI,SCOS=solar azimuth,cosine of solar angle
  !     BETAG=incident solar angle at ground surface
  !     GCOS,GSIN=cos,sin of ground surface
  !     ZNOON=hour of solar noon from weather file
  !
  IF(SSIN(NY,NX).GT.ZERO)THEN
    SAZI=0.2618*(ZNOON(NY,NX)-J)+4.7124
    SCOS=SQRT(1.0-SSIN(NY,NX)**2)
    DGAZI=COS(GAZI(NY,NX)-SAZI)
    BETAG=AMAX1(0.0,AMIN1(1.0,GCOS(NY,NX)*SSIN(NY,NX) &
      +GSIN(NY,NX)*SCOS*DGAZI))
    IF(ARLSS(NY,NX).GT.0.0)THEN
      SAGL=ASIN(SSIN(NY,NX))
      !
      !     ABSORBED RADIATION FROM OPTICAL PROPERTIES ENTERED IN 'READS'
      !
      !     RADSA,RADWA,RAPSA,RAPWA=SW,PAR absorbed at leaf,stalk surface
      !     perpendicular to incoming radiation
      !
      DO 1050 NZ=1,NP(NY,NX)
        RADSA(NZ,NY,NX)=RADS(NY,NX)*ABSR(NZ,NY,NX)
        RADWA(NZ,NY,NX)=RADS(NY,NX)*ABSRW
        RAPSA(NZ,NY,NX)=RAPS(NY,NX)*ABSP(NZ,NY,NX)
        RAPWA(NZ,NY,NX)=RAPS(NY,NX)*ABSPW
1050  CONTINUE
      !
      !     ANGLES BETWEEN SUN OR SKY ZONES AND FOLIAR SURFACES
      !
      !     ZAZI=leaf azimuth
      !     BETA,BETX=incident angle of direct radiation at leaf,horizontal surface
      !     ZAGL=determines forward vs backscattering
      !     IALBS=flag for forward vs backscattering
      !
      DO 1100 M=1,JSA
        ZAZI=SAZI+(M-0.5)*PICON/real(M,r8)
        DAZI=COS(ZAZI-SAZI)
        DO  N=1,JLI
          BETY=ZCOS(N)*SSIN(NY,NX)+ZSIN(N)*SCOS*DAZI
          BETA(N,M)=ABS(BETY)
          BETX(N,M)=BETA(N,M)/SSIN(NY,NX)
          IF(ZCOS(N).GT.SSIN(NY,NX))THEN
            BETZ=ACOS(BETY)
          ELSE
            BETZ=-ACOS(BETY)
          ENDIF
          IF(BETZ.GT.-PICON2)THEN
            ZAGL=SAGL+2.0*BETZ
          ELSE
            ZAGL=SAGL-2.0*(PICON+BETZ)
          ENDIF
          IF(ZAGL.GT.0.0.AND.ZAGL.LT.PICON)THEN
            IALBS(N,M)=1
          ELSE
            IALBS(N,M)=2
          ENDIF
!
          !     INTENSITY OF ABSORBED DIRECT RADIATION AT LEAF SURFACES
          !
          !     RDNDIR,RDNDIW,PARDIR,PARDIW=SW,PAR flux absorbed by leaf,stalk surfaces
          !     PAR,PARDIF=direct,diffuse PAR flux
          !     RADYL,RAPYL=solar beam diffuse SW,PAR flux
          !     RAFYL,RAFPL=forward scattered diffuse SW,PAR flux
          !     TAUS,TAUY=fraction of direct,diffuse radiation transmitted
!
          DO  NZ=1,NP(NY,NX)
            RDNDIR(N,M,NZ,NY,NX)=RADSA(NZ,NY,NX)*ABS(BETA(N,M))
            RDNDIW(N,M,NZ,NY,NX)=RADWA(NZ,NY,NX)*ABS(BETA(N,M))
            PARDIR(N,M,NZ,NY,NX)=RAPSA(NZ,NY,NX)*ABS(BETA(N,M))
            PARDIW(N,M,NZ,NY,NX)=RAPWA(NZ,NY,NX)*ABS(BETA(N,M))
            DO L=1,JC
              PARDIF(N,M,L,NZ,NY,NX)=0.0
              PAR(N,M,L,NZ,NY,NX)=PARDIR(N,M,NZ,NY,NX)
            enddo
          enddo
        enddo
1100  CONTINUE
      XAREA=1.00/AREA(3,NU(NY,NX),NY,NX)
      YAREA=0.25/AREA(3,NU(NY,NX),NY,NX)
      RADYL=RADY(NY,NX)
      RAPYL=RAPY(NY,NX)
      TAUS(JC+1,NY,NX)=1.0
      TAUY(JC+1)=1.0
      RAFSL(JC+1)=0.0
      RAFPL(JC+1)=0.0
      STOPS=0.0
      !
      !     RESET ARRAYS OF SUNLIT AND SHADED LEAF AREAS IN DIFFERENT
      !     LAYERS AND ANGLE CLASSES
      !
      !     TSURF,TSURFB,SURF,SURFB=leaf,stalk total,PFT surface area
!
      DO 1150 NZ=1,NP(NY,NX)
        DO  L=1,JC
          DO  N=1,JLI
            TSURF(N,L,NZ,NY,NX)=0.0
            TSURFB(N,L,NZ,NY,NX)=0.0
          enddo
        enddo
1150  CONTINUE
      DO 1200 NZ=1,NP(NY,NX)
        DO  NB=1,NBR(NZ,NY,NX)
          DO  L=1,JC
            IF(ZL(L-1,NY,NX).GT.DPTHS(NY,NX)-ZERO &
             .AND.ZL(L-1,NY,NX).GT.DPTH0(NY,NX)-ZERO)THEN
              DO 1205 N=1,JLI
                DO 1210 K=1,JNODS
                  TSURF(N,L,NZ,NY,NX)=TSURF(N,L,NZ,NY,NX)+SURF(N,L,K,NB,NZ,NY,NX)
1210            CONTINUE
                TSURFB(N,L,NZ,NY,NX)=TSURFB(N,L,NZ,NY,NX)+SURFB(N,L,NB,NZ,NY,NX)
1205          CONTINUE
            ENDIF
          enddo
        enddo
1200  CONTINUE
      !
      !     CALCULATE ABSORPTION, REFLECTION AND TRANSMISSION OF DIRECT AND
      !     DIFFUSE DOWNWARD TOTAL AND VISIBLE RADIATION BY EACH SPECIES
      !     NZ IN EACH LAYER L
      !
      !     RAFYL,RAFPL=forward scattered diffuse SW,PAR
      !     RABYL,RABPL=backscattered diffuse SW,PAR
      !     RADYL,RAPYL=solar beam diffuse SW,PAR
      !     STOPY,STOPSZ,STOPYZ=fraction of direct,diffuse radiation intercepted
      !
      DO 1800 L=JC,1,-1
        IF(ZL(L-1,NY,NX).GE.DPTHS(NY,NX)-ZERO &
          .AND.ZL(L-1,NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
          RADYL=RADYL*TAUY(L+1)+RAFSL(L+1)
          RAPYL=RAPYL*TAUY(L+1)+RAFPL(L+1)
          RAFSL(L)=0.0
          RAFPL(L)=0.0
          RABSL(L)=0.0
          RABPL(L)=0.0
          STOPY=0.0
          STOPSZ=0.0
          STOPYZ=0.0
    !
          !     RESET ACCUMULATORS OB ABSORBED, REFLECTED AND TRANSMITTED RADIATION
          !
          !     RADSL,RADSW,RADPL,RADPW=direct atmosph SW,PAR absbd by leaf,stalk surf
          !     RAYSL,RAYSW,RAYPL,RAYPW=diffuse atmosph SW,PAR absbd by leaf,stalk surf
          !     RADS1,RADW1,RADP1,RADQ1=backscattered direct SW,PAR absbd by leaf,stalk surf
          !     RAYS1,RAYW1,RAYP1,RAYQ1=backscattered diffuse SW,PAR absbd by leaf,stalk surf
          !     RADS2,RADW2,RADP2,RADQ2=fwd scattered direct SW,PAR absbd by leaf,stalk surf
          !     RAYS2,RAYW2,RAYP2,RAYQ2=fwd scattered diffuse SW,PAR absbd by leaf,stalk surf
    !
          DO 1500 NZ=1,NP(NY,NX)
            RADSL(NZ,NY,NX)=0.0
            RADSW(NZ,NY,NX)=0.0
            RADPL(NZ,NY,NX)=0.0
            RADPW(NZ,NY,NX)=0.0
            RAYSL(NZ,NY,NX)=0.0
            RAYSW(NZ,NY,NX)=0.0
            RAYPL(NZ,NY,NX)=0.0
            RAYPW(NZ,NY,NX)=0.0
            RADS1(NZ,NY,NX)=0.0
            RADW1(NZ,NY,NX)=0.0
            RADP1(NZ,NY,NX)=0.0
            RADQ1(NZ,NY,NX)=0.0
            RAYS1(NZ,NY,NX)=0.0
            RAYW1(NZ,NY,NX)=0.0
            RAYP1(NZ,NY,NX)=0.0
            RAYQ1(NZ,NY,NX)=0.0
            RADS2(NZ,NY,NX)=0.0
            RADW2(NZ,NY,NX)=0.0
            RADP2(NZ,NY,NX)=0.0
            RADQ2(NZ,NY,NX)=0.0
            RAYS2(NZ,NY,NX)=0.0
            RAYW2(NZ,NY,NX)=0.0
            RAYP2(NZ,NY,NX)=0.0
            RAYQ2(NZ,NY,NX)=0.0
    !
      !     LEAF SURFACE AREA IN EACH INCLINATION CLASS N, AZIMUTH CLASS M,
      !     LAYER L AND SPECIES NZ
      !
      !     TSURFY=unself-shaded leaf area
      !     TSURFZ=unself-shaded leaf area m-2 in each azimuth class
      !     TSURFS=TSURFY with shading from canopy layers above
      !     TSURFX=TSURFS m-2
      !     TSURWY=unself-shaded stalk area
      !     TSURWZ=unself-shaded stalk area m-2 in each azimuth class
      !     TSURWS=TSURWY with shading from canopy layers above
      !     TSURWX=TSURWS m-2
      !
            DO 1600 N=1,JLI
              TSURFY=TSURF(N,L,NZ,NY,NX)*CFX(NZ,NY,NX)
              TSURFZ=TSURFY*YAREA
              TSURFS=TSURFY*TAUS(L+1,NY,NX)
              TSURFX=TSURFS*XAREA
              TSURWY=TSURFB(N,L,NZ,NY,NX)*CFW
              TSURWZ=TSURWY*YAREA
              TSURWS=TSURWY*TAUS(L+1,NY,NX)
              TSURWX=TSURWS*XAREA
              !
              !     ABSORPTION OF DIRECT RADIATION BY SUNLIT LEAF SURFACES
              !
              !     STOPZ=accumulated horizontal area of intercepted direct radiation
              !
              DO 1700 M=1,JSA
                RADSL(NZ,NY,NX)=RADSL(NZ,NY,NX)+TSURFS*RDNDIR(N,M,NZ,NY,NX)
                RADSW(NZ,NY,NX)=RADSW(NZ,NY,NX)+TSURWS*RDNDIW(N,M,NZ,NY,NX)
                RADPL(NZ,NY,NX)=RADPL(NZ,NY,NX)+TSURFS*PARDIR(N,M,NZ,NY,NX)
                RADPW(NZ,NY,NX)=RADPW(NZ,NY,NX)+TSURWS*PARDIW(N,M,NZ,NY,NX)
                STOPSZ=STOPSZ+(TSURFX+TSURWX)*BETX(N,M)
!
          !     BACKSCATTERING OF REFLECTED DIRECT RADIATION
          !
                IF(IALBS(N,M).EQ.1)THEN
                  RADS1(NZ,NY,NX)=RADS1(NZ,NY,NX)+TSURFS*RDNDIR(N,M,NZ,NY,NX)
                  RADW1(NZ,NY,NX)=RADW1(NZ,NY,NX)+TSURWS*RDNDIW(N,M,NZ,NY,NX)
                  RADP1(NZ,NY,NX)=RADP1(NZ,NY,NX)+TSURFS*PARDIR(N,M,NZ,NY,NX)
                  RADQ1(NZ,NY,NX)=RADQ1(NZ,NY,NX)+TSURWS*PARDIW(N,M,NZ,NY,NX)
                  !
                  ! FORWARD SCATTERING OF REFLECTED DIRECT RADIATION
                  !
                ELSE
                  RADS2(NZ,NY,NX)=RADS2(NZ,NY,NX)+TSURFS*RDNDIR(N,M,NZ,NY,NX)
                  RADW2(NZ,NY,NX)=RADW2(NZ,NY,NX)+TSURWS*RDNDIW(N,M,NZ,NY,NX)
                  RADP2(NZ,NY,NX)=RADP2(NZ,NY,NX)+TSURFS*PARDIR(N,M,NZ,NY,NX)
                  RADQ2(NZ,NY,NX)=RADQ2(NZ,NY,NX)+TSURWS*PARDIW(N,M,NZ,NY,NX)
                ENDIF
    !
                !     INTENSITY OF ABSORBED DIFFUSE RADIATION AT LEAF SURFACES
                !
                !     RADYN,RADYW,RAPYN,RAPYW=diffuse SW,PAR flux absorbed by leaf,stalk surf
                !     OMEGA,OMEGX=incident angle of diffuse radn at leaf,horizontal surface
!
                DO 1750 NN=1,JLA
                  RADYN=RADYL*OMEGA(M,N,NN)*ABSR(NZ,NY,NX)
                  RADYW=RADYL*OMEGA(M,N,NN)*ABSRW
                  RAPYN=RAPYL*OMEGA(M,N,NN)*ABSP(NZ,NY,NX)
                  RAPYW=RAPYL*OMEGA(M,N,NN)*ABSPW
                  PARDIF(N,M,L,NZ,NY,NX)=PARDIF(N,M,L,NZ,NY,NX)+RAPYN
                  PAR(N,M,L,NZ,NY,NX)=PAR(N,M,L,NZ,NY,NX)+RAPYN
!
                  !     ABSORPTION OF DIFFUSE RADIATION BY SHADED LEAF SURFACES
                  !
                  !     STOPYZ=accumulated horizontal area of intercepted diffuse radiation
                  !
                  RAYSL(NZ,NY,NX)=RAYSL(NZ,NY,NX)+TSURFY*RADYN
                  RAYSW(NZ,NY,NX)=RAYSW(NZ,NY,NX)+TSURWY*RADYW
                  RAYPL(NZ,NY,NX)=RAYPL(NZ,NY,NX)+TSURFY*RAPYN
                  RAYPW(NZ,NY,NX)=RAYPW(NZ,NY,NX)+TSURWY*RAPYW
                  STOPYZ=STOPYZ+(TSURFZ+TSURWZ)*OMEGX(M,N,NN)
    !
                  !     BACKSCATTERING OF REFLECTED DIFFUSE RADIATION
                  !
                  IF(IALBY(M,N,NN).EQ.1)THEN
                    RAYS1(NZ,NY,NX)=RAYS1(NZ,NY,NX)+TSURFY*RADYN
                    RAYW1(NZ,NY,NX)=RAYW1(NZ,NY,NX)+TSURWY*RADYW
                    RAYP1(NZ,NY,NX)=RAYP1(NZ,NY,NX)+TSURFY*RAPYN
                    RAYQ1(NZ,NY,NX)=RAYQ1(NZ,NY,NX)+TSURWY*RAPYW
                    !
                    !     FORWARD SCATTERING OF REFLECTED DIFFUSE RADIATION
                    !
                  ELSE
                    RAYS2(NZ,NY,NX)=RAYS2(NZ,NY,NX)+TSURFY*RADYN
                    RAYW2(NZ,NY,NX)=RAYW2(NZ,NY,NX)+TSURWY*RADYW
                    RAYP2(NZ,NY,NX)=RAYP2(NZ,NY,NX)+TSURFY*RAPYN
                    RAYQ2(NZ,NY,NX)=RAYQ2(NZ,NY,NX)+TSURWY*RAPYW
                  ENDIF
1750            CONTINUE
1700          CONTINUE
1600        CONTINUE
1500      CONTINUE
          !
          !     ACCUMULATED INTERCEPTION BY CANOPY LAYER
          !
          !     XTAUS=interception of direct radiation in current layer
          !     STOPZ=accumulated interception of direct radiation from topmost layer
          !     TAUS=transmission of direct radiation to next lower layer
!
          IF(STOPS+STOPSZ.GT.1.0)THEN
            IF(STOPSZ.GT.ZERO)THEN
              XTAUS=(1.0-STOPS)/((1.0-STOPS)-(1.0-STOPS-STOPSZ))
            ELSE
              XTAUS=0.0
            ENDIF
            TAUS(L+1,NY,NX)=TAUS(L+1,NY,NX)*XTAUS
            STOPSZ=STOPSZ*XTAUS
            DO 1510 NZ=1,NP(NY,NX)
              RADSL(NZ,NY,NX)=RADSL(NZ,NY,NX)*XTAUS
              RADSW(NZ,NY,NX)=RADSW(NZ,NY,NX)*XTAUS
              RADPL(NZ,NY,NX)=RADPL(NZ,NY,NX)*XTAUS
              RADPW(NZ,NY,NX)=RADPW(NZ,NY,NX)*XTAUS
              RADS1(NZ,NY,NX)=RADS1(NZ,NY,NX)*XTAUS
              RADW1(NZ,NY,NX)=RADW1(NZ,NY,NX)*XTAUS
              RADP1(NZ,NY,NX)=RADP1(NZ,NY,NX)*XTAUS
              RADQ1(NZ,NY,NX)=RADQ1(NZ,NY,NX)*XTAUS
              RADS2(NZ,NY,NX)=RADS2(NZ,NY,NX)*XTAUS
              RADW2(NZ,NY,NX)=RADW2(NZ,NY,NX)*XTAUS
              RADP2(NZ,NY,NX)=RADP2(NZ,NY,NX)*XTAUS
              RADQ2(NZ,NY,NX)=RADQ2(NZ,NY,NX)*XTAUS
1510        CONTINUE
          ENDIF
!
          !     XTAUY=interception of diffuse radiation in current layer
          !     STOPYZ=accumulated interception of diffuse radiation from topmost layer
          !     TAUY=transmission of diffuse radiation to next lower layer
!
          IF(STOPY+STOPYZ.GT.1.0)THEN
            XTAUY=(1.0-STOPY)/((1.0-STOPY)-(1.0-STOPY-STOPYZ))
            TAUY(L+1)=TAUY(L+1)*XTAUY
            STOPYZ=STOPYZ*XTAUY
            DO 1520 NZ=1,NP(NY,NX)
              RAYSL(NZ,NY,NX)=RAYSL(NZ,NY,NX)*XTAUY
              RAYSW(NZ,NY,NX)=RAYSW(NZ,NY,NX)*XTAUY
              RAYPL(NZ,NY,NX)=RAYPL(NZ,NY,NX)*XTAUY
              RAYPW(NZ,NY,NX)=RAYPW(NZ,NY,NX)*XTAUY
              RAYS1(NZ,NY,NX)=RAYS1(NZ,NY,NX)*XTAUY
              RAYW1(NZ,NY,NX)=RAYW1(NZ,NY,NX)*XTAUY
              RAYP1(NZ,NY,NX)=RAYP1(NZ,NY,NX)*XTAUY
              RAYQ1(NZ,NY,NX)=RAYQ1(NZ,NY,NX)*XTAUY
              RAYS2(NZ,NY,NX)=RAYS2(NZ,NY,NX)*XTAUY
              RAYW2(NZ,NY,NX)=RAYW2(NZ,NY,NX)*XTAUY
              RAYP2(NZ,NY,NX)=RAYP2(NZ,NY,NX)*XTAUY
              RAYQ2(NZ,NY,NX)=RAYQ2(NZ,NY,NX)*XTAUY
              DO 1730 N=1,JLI
                DO  M=1,JSA
                  PARDIF(N,M,L,NZ,NY,NX)=PARDIF(N,M,L,NZ,NY,NX)*XTAUY
                  PAR(N,M,L,NZ,NY,NX)=PARDIR(N,M,NZ,NY,NX)+PARDIF(N,M,L,NZ,NY,NX)
                enddo
1730          CONTINUE
1520        CONTINUE
          ENDIF
          !
          !     TOTAL RADIATION ABSORBED, REFLECTED AND TRANSMITTED BY ALL PFTs
          !
          !     RADST,RADWT,RADPT,RADQT=total atmospheric SW,PAR absorbed by leaf,stalk
          !     RA1ST,RA1WT,RA1PT,RA1QT=total backscattered SW,PAR absd by leaf,stalk
          !     RA2ST,RA2WT,RA2PT,RA2QT=total fwd scattered SW,PAR absd by leaf,stalk
          !     RAFSL,RAFPL=total fwd scattered SW,PAR to next layer
          !     RABSL,RABPL=total back scattered SW,PAR to next layer
          !     RADC,TRADC,RADP,TRADP=total atmospheric SW,PAR absbd by each,all PFT
          !     STOPS,STOPY=accumulated interception of direct,diffuse radiation
          !     TAUS,TAUY=transmission of direct,diffuse radiation to next lower layer
          !
          DO 1530 NZ=1,NP(NY,NX)
            RADST=RADSL(NZ,NY,NX)+RAYSL(NZ,NY,NX)
            RADWT=RADSW(NZ,NY,NX)+RAYSW(NZ,NY,NX)
            RADPT=RADPL(NZ,NY,NX)+RAYPL(NZ,NY,NX)
            RADQT=RADPW(NZ,NY,NX)+RAYPW(NZ,NY,NX)
            RA1ST=RADS1(NZ,NY,NX)+RAYS1(NZ,NY,NX)
            RA1WT=RADW1(NZ,NY,NX)+RAYW1(NZ,NY,NX)
            RA1PT=RADP1(NZ,NY,NX)+RAYP1(NZ,NY,NX)
            RA1QT=RADP1(NZ,NY,NX)+RAYQ1(NZ,NY,NX)
            RA2ST=RADS2(NZ,NY,NX)+RAYS2(NZ,NY,NX)
            RA2WT=RADW2(NZ,NY,NX)+RAYW2(NZ,NY,NX)
            RA2PT=RADP2(NZ,NY,NX)+RAYP2(NZ,NY,NX)
            RA2QT=RADQ2(NZ,NY,NX)+RAYQ2(NZ,NY,NX)
            RAFSL(L)=RAFSL(L)+(RADST*TAUR(NZ,NY,NX) &
              +RA2ST*ALBR(NZ,NY,NX)+RA2WT*ALBRW)*YAREA
            RAFPL(L)=RAFPL(L)+(RADPT*TAUP(NZ,NY,NX) &
              +RA2PT*ALBP(NZ,NY,NX)+RA2QT*ALBPW)*YAREA
            RABSL(L)=RABSL(L)+(RA1ST*ALBR(NZ,NY,NX)+RA1WT*ALBRW)*YAREA
            RABPL(L)=RABPL(L)+(RA1PT*ALBP(NZ,NY,NX)+RA1QT*ALBPW)*YAREA
            RADC(NZ,NY,NX)=RADC(NZ,NY,NX)+RADST+RADWT
            RADP(NZ,NY,NX)=RADP(NZ,NY,NX)+RADPT+RADQT
            TRADC(NY,NX)=TRADC(NY,NX)+RADST+RADWT
            TRAPC(NY,NX)=TRAPC(NY,NX)+RADPT+RADQT
1530      CONTINUE
          STOPS=STOPS+STOPSZ
          STOPY=STOPY+STOPYZ
          TAUS(L,NY,NX)=1.0-STOPS
          TAU0(L,NY,NX)=1.0-TAUS(L,NY,NX)
          TAUY(L)=1.0-STOPY
        ELSE
          RAFSL(L)=RAFSL(L+1)
          RAFPL(L)=RAFPL(L+1)
          TAUS(L,NY,NX)=TAUS(L+1,NY,NX)
          TAU0(L,NY,NX)=1.0-TAUS(L,NY,NX)
          TAUY(L)=TAUY(L+1)
        ENDIF
1800  CONTINUE
      !
      !     DIRECT AND DIFFUSE RADIATION ABSORBED AT GROUND SURFACE
      !
      !     RADSG,RADYG,RAPSG,RAPYG=direct,diffuse SW,PAR at horizl ground surface
      !     RADS,RAPS =solar beam direct SW,PAR flux
      !     TAUS,TAUY=transmission of direct,diffuse radiation below canopy
      !     RADYL,RAPYL=solar beam diffuse SW,PAR flux
      !     RASG,RAPG=total SW,PAR at ground surface
      !     BETAG,OMEGAG=incident solar,sky angle at ground surface
!
      RADSG=RADS(NY,NX)*TAUS(1,NY,NX)
      RADYG=RADYL*TAUY(1)+RAFSL(1)
      RAPSG=RAPS(NY,NX)*TAUS(1,NY,NX)
      RAPYG=RAPYL*TAUY(1)+RAFPL(1)
      RASG=ABS(BETAG)*RADSG
      RAPG=ABS(BETAG)*RAPSG
      DO 20 N=1,JSA
        RASG=RASG+ABS(OMEGAG(N,NY,NX))*RADYG
        RAPG=RAPG+ABS(OMEGAG(N,NY,NX))*RAPYG
20    CONTINUE
      RADG(NY,NX)=RASG*AREA(3,NU(NY,NX),NY,NX)
!
      !     RADIATION REFLECTED FROM GROUND SURFACE
      !
      !     VHCPW,VHCPWX=current,minimum snowpack heat capacity
      !     ALBW,VOLSS,VOLWS,VOLIS=snowpack albedo,snow,water,ice volume
      !     ALBG,ALBS,FSNOW=ground,soil albedo,snow cover fraction
      !     THETW1=soil surface water content
      !     RABSL,RADPL=SW,PAR backscatter from ground surface
      !     TRADG,TRAPG=SW,PAR absorbed by ground surface
!
      IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
        ALBW=(0.80*VOLSS(NY,NX)+0.30*VOLIS(NY,NX)+0.06*VOLWS(NY,NX)) &
          /(VOLSS(NY,NX)+VOLIS(NY,NX)+VOLWS(NY,NX))
        FSNOW=AMIN1((DPTHS(NY,NX)/0.07)**2,1.0)
        ALBG=FSNOW*ALBW+(1.0-FSNOW)*ALBS(NY,NX)
      ELSE
        IF(VOLX(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
          THETW1=AMIN1(POROS(NU(NY,NX),NY,NX) &
            ,VOLW(NU(NY,NX),NY,NX)/VOLY(NU(NY,NX),NY,NX))
        ELSE
          THETW1=0.0
        ENDIF
        ALBG=AMIN1(ALBX(NY,NX),ALBS(NY,NX) &
          +AMAX1(0.0,ALBX(NY,NX)-THETW1))
      ENDIF
      RABSL(0)=RASG*ALBG*0.25
      RABPL(0)=RAPG*ALBG*0.25
      TRADG(NY,NX)=(1.0-ALBG)*RASG*AREA(3,NU(NY,NX),NY,NX)
      TRAPG(NY,NX)=(1.0-ALBG)*RAPG*AREA(3,NU(NY,NX),NY,NX)
!
      !     ADD RADIATION FROM SCATTERING THROUGH CANOPY LAYERS
      !
      !     RABSL,RABPL=total backscattered SW,PAR to next layer
      !     RAFSL,RAFPL=total fwd scattered SW,PAR to next layer
      !     RADYN,RADYW,RAPYN,RAPYW=leaf,stalk SW,PAR absbd fwd+back flux
      !     RAYSL,RAYSW,RAYPL,RAYPW=total leaf,stalk SW,PAR absbd fwd+back
      !     RADC,TRADC,RADP,TRADP=total SW,PAR absbd by each,all PFT
!
      RADYL=0.0
      RAPYL=0.0
      TAUY(0)=1.0
      RAFSL(0)=0.0
      RAFPL(0)=0.0
      DO 2800 L=1,JC
        IF(ZL(L-1,NY,NX).GE.DPTHS(NY,NX)-ZERO &
         .AND.ZL(L-1,NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
          RADYL=RADYL*TAUY(L-1)+RAFSL(L-1)+RABSL(L-1)
          RAPYL=RAPYL*TAUY(L-1)+RAFPL(L-1)+RABPL(L-1)
          RAFSL(L)=0.0
          RAFPL(L)=0.0
          DO 2500 NZ=1,NP(NY,NX)
            RAYSL(NZ,NY,NX)=0.0
            RAYSW(NZ,NY,NX)=0.0
            RAYPL(NZ,NY,NX)=0.0
            RAYPW(NZ,NY,NX)=0.0
            DO 2600 N=1,JLI
              TSURFY=TSURF(N,L,NZ,NY,NX)*CFX(NZ,NY,NX)
              TSURWY=TSURFB(N,L,NZ,NY,NX)*CFW
              DO 2700 M=1,JSA
                DO 2750 NN=1,JLA
                  RADYN=RADYL*OMEGA(M,N,NN)*ABSR(NZ,NY,NX)
                  RADYW=RADYL*OMEGA(M,N,NN)*ABSRW
                  RAPYN=RAPYL*OMEGA(M,N,NN)*ABSP(NZ,NY,NX)
                  RAPYW=RAPYL*OMEGA(M,N,NN)*ABSPW
                  PARDIF(N,M,L,NZ,NY,NX)=PARDIF(N,M,L,NZ,NY,NX)+RAPYN
                  PAR(N,M,L,NZ,NY,NX)=PAR(N,M,L,NZ,NY,NX)+RAPYN
                  RAYSL(NZ,NY,NX)=RAYSL(NZ,NY,NX)+TSURFY*RADYN
                  RAYSW(NZ,NY,NX)=RAYSW(NZ,NY,NX)+TSURWY*RADYW
                  RAYPL(NZ,NY,NX)=RAYPL(NZ,NY,NX)+TSURFY*RAPYN
                  RAYPW(NZ,NY,NX)=RAYPW(NZ,NY,NX)+TSURWY*RAPYW
2750            CONTINUE
2700          CONTINUE
2600        CONTINUE
            RAFSL(L)=RAFSL(L)+RAYSL(NZ,NY,NX)*TAUR(NZ,NY,NX)*YAREA
            RAFPL(L)=RAFPL(L)+RAYPL(NZ,NY,NX)*TAUP(NZ,NY,NX)*YAREA
            RADC(NZ,NY,NX)=RADC(NZ,NY,NX)+RAYSL(NZ,NY,NX)+RAYSW(NZ,NY,NX)
            RADP(NZ,NY,NX)=RADP(NZ,NY,NX)+RAYPL(NZ,NY,NX)+RAYPW(NZ,NY,NX)
            TRADC(NY,NX)=TRADC(NY,NX)+RAYSL(NZ,NY,NX)+RAYSW(NZ,NY,NX)
            TRAPC(NY,NX)=TRAPC(NY,NX)+RAYPL(NZ,NY,NX)+RAYPW(NZ,NY,NX)
2500      CONTINUE
        ELSE
          RAFSL(L)=RAFSL(L-1)
          RAFPL(L)=RAFPL(L-1)
          RABSL(L)=RABSL(L-1)
          RABPL(L)=RABPL(L-1)
        ENDIF
2800  CONTINUE
!
      !     RADIATION AT GROUND SURFACE IF NO CANOPY
!
    ELSE
      RASG=ABS(BETAG)*RADS(NY,NX)
      DO 120 N=1,JSA
        RASG=RASG+ABS(OMEGAG(N,NY,NX))*RADY(NY,NX)
120   CONTINUE
      RADG(NY,NX)=RASG*AREA(3,NU(NY,NX),NY,NX)
      DO 135 NZ=1,NP(NY,NX)
        RADC(NZ,NY,NX)=0.0
        RADP(NZ,NY,NX)=0.0
135   CONTINUE
    ENDIF
!
    !     IF NO RADIATION
!
  ELSE
    RADG(NY,NX)=0.0
    DO 125 NZ=1,NP(NY,NX)
      RADC(NZ,NY,NX)=0.0
      RADP(NZ,NY,NX)=0.0
125 CONTINUE
  ENDIF
  !
  !     CANOPY AND GROUND SKY FRACTIONS USED FOR BOUNDARY LAYER CALCULNS
  !
  !     FRADG=fraction of radiation received by ground surface
  !     FRADP=fraction of radiation received by each PFT canopy
  !     ARLSS,ARLFS=leaf+stalk area of all PFTs,each PFT
  !
  FRADG(NY,NX)=1.0
  IF(ARLSS(NY,NX).GT.ZEROS(NY,NX))THEN
    FRADPT=1.0-EXP(-0.65*ARLSS(NY,NX)/AREA(3,NU(NY,NX),NY,NX))
    DO 145 NZ=1,NP(NY,NX)
      FRADP(NZ,NY,NX)=FRADPT*ARLFS(NZ,NY,NX)/ARLSS(NY,NX)
      FRADG(NY,NX)=FRADG(NY,NX)-FRADP(NZ,NY,NX)
145 CONTINUE
  ELSE
    FRADG(NY,NX)=1.0
    DO 146 NZ=1,NP(NY,NX)
      FRADP(NZ,NY,NX)=0.0
146 CONTINUE
  ENDIF
  end subroutine MultiLayerSurfaceRadiation
end module CanopyCondMod
