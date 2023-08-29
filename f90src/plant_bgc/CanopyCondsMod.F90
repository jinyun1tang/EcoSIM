module CanopyCondsMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSimConst
  use EcoSIMConfig
  use PlantAPIData
  use minimathmod, only : AZMAX1,isnan
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

  integer :: curday
  contains
  subroutine CanopyConditionModel(I,J,DPTH0)
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(in) :: DPTH0 !water+ice depth at surface

  curday=I
  call MultiLayerSurfaceRadiation(I,J,DPTH0)

  call DivideCanopyLayerByLAI()

  call CalcBoundaryLayerProperties(DPTH0)

  end subroutine CanopyConditionModel

!------------------------------------------------------------------------------------------

  subroutine CalcBoundaryLayerProperties(DPTH0)
  implicit none
  real(r8), intent(in) :: DPTH0
  real(r8) :: ARLSC
  real(r8) :: ARLSG
  real(r8) :: ZX,ZY,ZE
  REAL(R8) :: ZZ
!     begin_execution
  associate(                          &
    UA      => plt_site%UA      , &
    ZERO    => plt_site%ZERO    , &
    AREA3   => plt_site%AREA3   , &
    IETYP   => plt_site%IETYP   , &
    SoiSurfRoughnesst0      => plt_site%SoiSurfRoughnesst0      , &
    WindMesHeight      => plt_site%WindMesHeight      , &
    ZEROS   => plt_site%ZEROS   , &
    NU      => plt_site%NU      , &
    BndlResistAboveCanG     => plt_ew%BndlResistAboveCanG       , &
    ZeroPlanDisp      => plt_ew%ZeroPlanDisp        , &
    RoughHeight      => plt_ew%RoughHeight        , &
    RIB     => plt_ew%RIB       , &
    TairK     => plt_ew%TairK       , &
    VLHeatCapSnowMN  => plt_ew%VLHeatCapSnowMN    , &
    SnowDepth   => plt_ew%SnowDepth     , &
    VHCPW1  => plt_ew%VHCPW1    , &
    GridMaxCanopyHeight      => plt_morph%GridMaxCanopyHeight     , &
    IRTYP   => plt_morph%IRTYP  , &
    StemAreag   => plt_morph%StemAreag  , &
    CanGLA   => plt_morph%CanGLA    &
  )
!     CANOPY ZERO PLANE AND ROUGHNESS HEIGHTS
!
!     CanGLA,StemAreag=leaf,stalk area of combined canopy
!     SnowDepth,DPTH0=snowpack,surface water depths
!     ZT,ZeroPlanDisp,RoughHeight=canopy,zero plane displacement,roughness height
!     ZZ=reference height for wind speed
!
  ARLSC=CanGLA+StemAreag
  IF(ARLSC.GT.ZEROS.AND.GridMaxCanopyHeight.GE.SnowDepth-ZERO.AND.GridMaxCanopyHeight.GE.DPTH0-ZERO)THEN
    ARLSG=ARLSC/AREA3(NU)
    ZX=EXP(-0.5_r8*ARLSG)
    ZY=1.0_r8-ZX
    ZeroPlanDisp=GridMaxCanopyHeight*AZMAX1(1.0_r8-2.0_r8/ARLSG*ZY)
    ZE=GridMaxCanopyHeight*AMAX1(0.05_r8,ZX*ZY)
  ELSE
    ZeroPlanDisp=0.0_r8
    ZE=0.0_r8
  ENDIF
  IF(IFLGW.EQ.1)THEN
    ZZ=WindMesHeight+GridMaxCanopyHeight
  ELSE
    ZZ=AMAX1(WindMesHeight,ZeroPlanDisp+2.0_r8)
  ENDIF
  IF(IETYP.GE.0)THEN
    IF(VHCPW1.GT.VLHeatCapSnowMN)THEN
      RoughHeight=AMAX1(0.001,ZE,ZW)
    ELSE
      RoughHeight=AMAX1(0.001,ZE,SoiSurfRoughnesst0)
    ENDIF
!
!     CANOPY ISOTHERMAL BOUNDARY LAYER RESISTANCE
!
!     BndlResistAboveCanG,RAM=biome canopy,minimum isothermal boundary layer resistance
!     UA=wind speed
!     RIB=canopy isothermal Richardson number
!
    BndlResistAboveCanG=AMAX1(RAM,(LOG((ZZ-ZeroPlanDisp)/RoughHeight))**2._r8/(0.168_r8*UA))
    RIB=1.27E+08_r8*(ZZ-RoughHeight)/(UA**2*TairK)
  ELSE
    BndlResistAboveCanG=RAM
    RIB=0.0_r8
  ENDIF
  end associate
  end subroutine CalcBoundaryLayerProperties

!------------------------------------------------------------------------------------------

  subroutine DivideCanopyLayerByLAI()
  implicit none

  real(r8) :: ZL1(0:JZ1)
  real(r8) :: ART,ARL
  real(r8) :: ARX
  real(r8) :: DZL
  integer :: NZ,L
  !     begin_execution
  associate(                      &
    NP      => plt_site%NP      , &
    ZEROS   => plt_site%ZEROS   , &
    GridMaxCanopyHeight      => plt_morph%GridMaxCanopyHeight     , &
    CanopyHeightz      => plt_morph%CanopyHeightz     , &
    CanopyHeight      => plt_morph%CanopyHeight     , &
    ARSTT   => plt_morph%ARSTT  , &
    ARLFT   => plt_morph%ARLFT  , &
    StemAreag   => plt_morph%StemAreag  , &
    CanGLA   => plt_morph%CanGLA    &
  )
  !
  !     DIVISION OF CANOPY INTO JC LAYERS WITH EQUAL LAI
  !
  !     ZT,ZC=heights of combined canopy,PFT canopy
  !     ZL=height to bottom of each canopy layer
  !     CanGLA,StemAreag=leaf,stalk area of combined canopy
  !     ARLFT,ARSTT=leaf,stalk area of combined canopy layer
  !
  GridMaxCanopyHeight=0.0
  D9685: DO NZ=1,NP
    GridMaxCanopyHeight=AMAX1(GridMaxCanopyHeight,CanopyHeight(NZ))
  ENDDO D9685
  CanopyHeightz(JC1)=GridMaxCanopyHeight+0.01_r8
  ZL1(JC1)=CanopyHeightz(JC1)
  ZL1(0)=0.0
  ART=(CanGLA+StemAreag)/JC1
  IF(ART.GT.ZEROS)THEN
    DO 2765 L=JC1,2,-1
      ARL=ARLFT(L)+ARSTT(L)
      IF(ARL.GT.1.01*ART)THEN
        DZL=CanopyHeightz(L)-CanopyHeightz(L-1)
        ZL1(L-1)=CanopyHeightz(L-1)+0.5*AMIN1(1.0,(ARL-ART)/ARL)*DZL
      ELSEIF(ARL.LT.0.99*ART)THEN
        ARX=ARLFT(L-1)+ARSTT(L-1)
        DZL=CanopyHeightz(L-1)-CanopyHeightz(L-2)
        IF(ARX.GT.ZEROS)THEN
          ZL1(L-1)=CanopyHeightz(L-1)-0.5_r8*AMIN1(1.0_r8,(ART-ARL)/ARX)*DZL
        ELSE
          ZL1(L-1)=CanopyHeightz(L-1)
        ENDIF
      ELSE
        ZL1(L-1)=CanopyHeightz(L-1)
      ENDIF
2765  CONTINUE
    D2770: DO L=JC1,2,-1
      CanopyHeightz(L-1)=ZL1(L-1)
!     CanopyHeightz(L-1)=AZMAX1(AMIN1(CanopyHeightz(L)-ppmc
!    2,CanopyHeightz(L-1)))
    ENDDO D2770
  ENDIF
  end associate
  end subroutine DivideCanopyLayerByLAI

!------------------------------------------------------------------------------------------

  subroutine MultiLayerSurfaceRadiation(I,J,DPTH0)
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(in) :: DPTH0
  integer :: NB,NZ,L,K,M,N,NN
  integer :: IALBS(JLI1,JSA1)
  real(r8) :: TAUY(0:JC1+1)
  real(r8) :: PARDIR(JLI1,JSA1,JP1)
  real(r8) :: PARDIW(JLI1,JSA1,JP1)
  real(r8) :: RABSL(0:JC1+1)
  real(r8) :: RABPL(0:JC1+1)
  real(r8) :: RAFSL(0:JC1+1)
  real(r8) :: RAFPL(0:JC1+1)
  real(r8) :: RADSL(JP1)
  real(r8) :: RADPL(JP1)
  real(r8) :: RAYSW(JP1)
  real(r8) :: RAYPW(JP1)
  real(r8) :: RADSW(JP1)
  real(r8) :: RADPW(JP1)
  real(r8) :: RADWA(JP1)
  real(r8) :: RADSA(JP1)
  real(r8) :: RAPSA(JP1)
  REAL(R8) :: RAPWA(JP1)
  real(r8) :: RADS1(JP1),RADS2(JP1)
  real(r8) :: RADP1(JP1),RADP2(JP1)
  real(r8) :: RADQ1(JP1),RADQ2(JP1)
  real(r8) :: RADW1(JP1),RADW2(JP1)
  real(r8) :: RAYSL(JP1),RAYPL(JP1)
  real(r8) :: RAYS1(JP1),RAYS2(JP1)
  real(r8) :: RAYP1(JP1),RAYP2(JP1)
  real(r8) :: RAYW1(JP1),RAYW2(JP1)
  real(r8) :: RAYQ1(JP1),RAYQ2(JP1)
  real(r8) :: BETA(JLI1,JSA1)                         !sine of direct solar radiation on leaf surface, [-]
  real(r8) :: BETX(JLI1,JSA1)                         !sine of direct solar radiation on leaf surface/sine of direct solar radiation, [-]
  REAL(R8) :: RDNDIR(JLI1,JSA1,JP1),RDNDIW(JLI1,JSA1,JP1)
  real(r8) :: TSURF(JLI1,JZ1,JP1),TSURFB(JLI1,JZ1,JP1)
  real(r8) :: TRADC,TRAPC,TRADG,TRAPG
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
  real(r8) :: SAZI,SCOS
  !     begin_execution
  !     MULTILAYER CANOPY INTERECEPTION OF DIRECT AND DIFFUSE RADIATION
  !     IN SW AND VISIBLE BANDS BY INCLINATION N, AZIMUTH M, LAYER L,
  !     NODE K, BRANCH NB, PFT NZ
  !
  !     CanPA,CanGA=leaf+stalk area of combined,each PFT canopy
  !     ZL=height to bottom of canopy layer
  !     SnowDepth,DPTH0=snowpack,surface water depths
  !     CanPLNBLA,CanPLBSA=leaf,stalk areas of PFT
  !     RAD,RAP=vertical direct+diffuse SW,PAR
  !     RADS,RADY,RAPS,RAPY=solar beam direct,diffuse SW,PAR
  !     SSIN,TYSIN=sine of solar,sky angles
  !     SWRadByCanP,RADP=total SW,PAR absorbed by canopy
  !     CFX=clumping factor for self-shading
  !
  associate(                       &
    VHCPW1  => plt_ew%VHCPW1 , &
    VcumIceSnow   => plt_ew%VcumIceSnow  , &
    VcumDrySnoWE   => plt_ew%VcumDrySnoWE  , &
    VcumWatSnow   => plt_ew%VcumWatSnow  , &
    VLHeatCapSnowMN  => plt_ew%VLHeatCapSnowMN , &
    SnowDepth   => plt_ew%SnowDepth  , &
    ALBR    => plt_rad%ALBR  , &
    ALBS    => plt_rad%ALBS  , &
    ALBX    => plt_rad%ALBX  , &
    ALBP    => plt_rad%ALBP  , &
    FRADG   => plt_rad%FRADG , &
    GAZI    => plt_rad%GAZI  , &
    GCOS    => plt_rad%GCOS  , &
    GSIN    => plt_rad%GSIN  , &
    IALBY   => plt_rad%IALBY , &
    OMEGA   => plt_rad%OMEGA , &
    FracPARByCanP   => plt_rad%FracPARByCanP , &
    PAR     => plt_rad%PAR   , &
    PARDIF  => plt_rad%PARDIF, &
    OMEGAG  => plt_rad%OMEGAG   , &
    OMEGX   => plt_rad%OMEGX    , &
    RAD0    => plt_rad%RAD0     , &
    RADG    => plt_rad%RADG     , &
    RADS    => plt_rad%RADS     , &
    RAP0    => plt_rad%RAP0     , &
    SWRadByCanP    => plt_rad%SWRadByCanP     , &
    RAPY    => plt_rad%RAPY     , &
    RAPS    => plt_rad%RAPS     , &
    RADY    => plt_rad%RADY     , &
    SSIN    => plt_rad%SSIN     , &
    TYSIN   => plt_rad%TYSIN    , &
    TAU0    => plt_rad%TAU0     , &
    TAUS    => plt_rad%TAUS     , &
    ZSIN    => plt_rad%ZSIN     , &
    ABSP    => plt_rad%ABSP     , &
    ABSR    => plt_rad%ABSR     , &
    TAUP    => plt_rad%TAUP     , &
    TAUR    => plt_rad%TAUR     , &
    PARByCanP    => plt_rad%PARByCanP     , &
    ZCOS    => plt_rad%ZCOS     , &
    ZEROS   => plt_site%ZEROS   , &
    NU      => plt_site%NU      , &
    AREA3   => plt_site%AREA3   , &
    NP      => plt_site%NP      , &
    ZERO    => plt_site%ZERO    , &
    ZEROS2  => plt_site%ZEROS2  , &
    POROS1  => plt_site%POROS1  , &
    ZNOON   => plt_site%ZNOON   , &
    VLSoilPoreMicP    => plt_soilchem%VLSoilPoreMicP, &
    VLSoilMicP    => plt_soilchem%VLSoilMicP, &
    VLWatMicP    => plt_soilchem%VLWatMicP, &
    CFX     => plt_morph%CFX    , &
    CanopyHeightz      => plt_morph%CanopyHeightz     , &
    NBR     => plt_morph%NBR    , &
    SURF    => plt_morph%SURF   , &
    SURFB   => plt_morph%SURFB  , &
    CanPLA   => plt_morph%CanPLA  , &
    CanPLNBLA   => plt_morph%CanPLNBLA  , &
    CanPLBSA   => plt_morph%CanPLBSA  , &
    CanPA   => plt_morph%CanPA  , &
    CF      => plt_morph%CF     , &
    CanGA   => plt_morph%CanGA    &
  )
  CanGA=0.0
  D1135: DO NZ=1,NP
    CanPA(NZ)=0.0
    DO  NB=1,NBR(NZ)
      DO  L=1,JC1

        IF(CanopyHeightz(L-1).GE.SnowDepth-ZERO.AND.CanopyHeightz(L-1).GE.DPTH0-ZERO)THEN
          !above snow depth and above water/ice surface
          D1130: DO K=1,JNODS1
            CanPA(NZ)=CanPA(NZ)+CanPLNBLA(L,K,NB,NZ)
            CanGA=CanGA+CanPLNBLA(L,K,NB,NZ)
          ENDDO D1130
          !add stem/stalk area
          CanPA(NZ)=CanPA(NZ)+CanPLBSA(L,NB,NZ)
          CanGA=CanGA+CanPLBSA(L,NB,NZ)
        ENDIF
      enddo
    enddo
  ENDDO D1135
  IF(SSIN.GT.ZERO)THEN
    RAD0=RADS*SSIN+RADY*TYSIN
    RAP0=RAPS*SSIN+RAPY*TYSIN
  ELSE
    RADS=0.0_r8
    RADY=0.0_r8
    RAPS=0.0_r8
    RAPY=0.0_r8
    RAD0=0.0_r8
    RAP0=0.0_r8
  ENDIF
  TRADC=0.0_r8
  TRAPC=0.0_r8
  D1025: DO NZ=1,NP
    SWRadByCanP(NZ)=0.0_r8
    PARByCanP(NZ)=0.0_r8
    CFX(NZ)=CF(NZ)*(1.0_r8-0.025_r8*CanPLA(NZ)/AREA3(NU))
  ENDDO D1025
  !
  !     ANGLE BETWEEN SUN AND GROUND SURFACE
  !
  !     SAZI,SCOS=solar azimuth,cosine of solar angle
  !     BETAG=incident solar angle at ground surface
  !     GCOS,GSIN=cos,sin of ground surface
  !     ZNOON=hour of solar noon from weather file
  !     0.2618=pi/12 (hrs)
  IF(SSIN.GT.ZERO)THEN
    SAZI=0.2618_r8*(ZNOON-J)+4.7124_r8
    SCOS=SQRT(1.0_r8-SSIN**2._r8)
    DGAZI=COS(GAZI-SAZI)
    BETAG=AZMAX1(AMIN1(1.0_r8,GCOS*SSIN+GSIN*SCOS*DGAZI))
    IF(CanGA.GT.0.0_r8)THEN
      SAGL=ASIN(SSIN)
      !
      !     ABSORBED RADIATION FROM OPTICAL PROPERTIES ENTERED IN 'READS'
      !
      !     RADSA,RADWA,RAPSA,RAPWA=SW,PAR absorbed at leaf,stalk surface
      !     perpendicular to incoming radiation
      !
      D1050: DO NZ=1,NP
        RADSA(NZ)=RADS*ABSR(NZ)
        RADWA(NZ)=RADS*ABSRW
        RAPSA(NZ)=RAPS*ABSP(NZ)
        RAPWA(NZ)=RAPS*ABSPW
      ENDDO D1050
      !
      !     ANGLES BETWEEN SUN OR SKY ZONES AND FOLIAR SURFACES
      !
      !     ZAZI=leaf azimuth
      !     BETA,BETX=incident angle of direct radiation at leaf,horizontal surface
      !     ZAGL=determines forward vs backscattering
      !     IALBS=flag for forward vs backscattering
      !
      D1100: DO M=1,JSA1
        ZAZI=SAZI+(M-0.5_r8)*PICON/real(M,r8)
        DAZI=COS(ZAZI-SAZI)
        DO  N=1,JLI1
          BETY=ZCOS(N)*SSIN+ZSIN(N)*SCOS*DAZI
          BETA(N,M)=ABS(BETY)
          BETX(N,M)=BETA(N,M)/SSIN
          IF(ZCOS(N).GT.SSIN)THEN
            BETZ=ACOS(BETY)
          ELSE
            BETZ=-ACOS(BETY)
          ENDIF
          IF(BETZ.GT.-PICON2h)THEN
            ZAGL=SAGL+2.0_r8*BETZ
          ELSE
            ZAGL=SAGL-2.0_r8*(PICON+BETZ)
          ENDIF
          IF(ZAGL.GT.0.0_r8.AND.ZAGL.LT.PICON)THEN
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
          DO  NZ=1,NP
            RDNDIR(N,M,NZ)=RADSA(NZ)*ABS(BETA(N,M))
            RDNDIW(N,M,NZ)=RADWA(NZ)*ABS(BETA(N,M))
            PARDIR(N,M,NZ)=RAPSA(NZ)*ABS(BETA(N,M))
            PARDIW(N,M,NZ)=RAPWA(NZ)*ABS(BETA(N,M))
            DO L=1,JC1
              PARDIF(N,M,L,NZ)=0.0_r8
              PAR(N,M,L,NZ)=PARDIR(N,M,NZ)
            enddo
          enddo
        enddo
      ENDDO D1100
      XAREA=1.00_r8/AREA3(NU)
      YAREA=0.25_r8/AREA3(NU)
      RADYL=RADY
      RAPYL=RAPY
      TAUS(JC1+1)=1.0_r8
      TAUY(JC1+1)=1.0_r8
      RAFSL(JC1+1)=0.0_r8
      RAFPL(JC1+1)=0.0_r8
      STOPS=0.0_r8
      !
      !     RESET ARRAYS OF SUNLIT AND SHADED LEAF AREAS IN DIFFERENT
      !     LAYERS AND ANGLE CLASSES
      !
      !     TSURF,TSURFB,SURF,SURFB=leaf,stalk total,PFT surface area
!
      D1150: DO NZ=1,NP
        DO  L=1,JC1
          DO  N=1,JLI1
            TSURF(N,L,NZ)=0.0_r8
            TSURFB(N,L,NZ)=0.0_r8
          enddo
        enddo
      ENDDO D1150
      D1200: DO NZ=1,NP
        DO  NB=1,NBR(NZ)
          DO  L=1,JC1
            IF(CanopyHeightz(L-1).GT.SnowDepth-ZERO.AND.CanopyHeightz(L-1).GT.DPTH0-ZERO)THEN
              D1205: DO N=1,JLI1
                D1210: DO K=1,JNODS1
                  TSURF(N,L,NZ)=TSURF(N,L,NZ)+SURF(N,L,K,NB,NZ)
                ENDDO D1210
                TSURFB(N,L,NZ)=TSURFB(N,L,NZ)+SURFB(N,L,NB,NZ)
              ENDDO D1205
            ENDIF
          enddo
        enddo
      ENDDO D1200
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
      D1800: DO L=JC1,1,-1
        IF(CanopyHeightz(L-1).GE.SnowDepth-ZERO.AND.CanopyHeightz(L-1).GE.DPTH0-ZERO)THEN
          RADYL=RADYL*TAUY(L+1)+RAFSL(L+1)
          RAPYL=RAPYL*TAUY(L+1)+RAFPL(L+1)
          RAFSL(L)=0.0_r8
          RAFPL(L)=0.0_r8
          RABSL(L)=0.0_r8
          RABPL(L)=0.0_r8
          STOPY=0.0_r8
          STOPSZ=0.0_r8
          STOPYZ=0.0_r8
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
          D1500: DO NZ=1,NP
            RADSL(NZ)=0.0_r8
            RADSW(NZ)=0.0_r8
            RADPL(NZ)=0.0_r8
            RADPW(NZ)=0.0_r8
            RAYSL(NZ)=0.0_r8
            RAYSW(NZ)=0.0_r8
            RAYPL(NZ)=0.0_r8
            RAYPW(NZ)=0.0_r8
            RADS1(NZ)=0.0_r8
            RADW1(NZ)=0.0_r8
            RADP1(NZ)=0.0_r8
            RADQ1(NZ)=0.0_r8
            RAYS1(NZ)=0.0_r8
            RAYW1(NZ)=0.0_r8
            RAYP1(NZ)=0.0_r8
            RAYQ1(NZ)=0.0_r8
            RADS2(NZ)=0.0_r8
            RADW2(NZ)=0.0_r8
            RADP2(NZ)=0.0_r8
            RADQ2(NZ)=0.0_r8
            RAYS2(NZ)=0.0_r8
            RAYW2(NZ)=0.0_r8
            RAYP2(NZ)=0.0_r8
            RAYQ2(NZ)=0.0_r8
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
            D1600: DO N=1,JLI1
              TSURFY=TSURF(N,L,NZ)*CFX(NZ)
              TSURFZ=TSURFY*YAREA
              TSURFS=TSURFY*TAUS(L+1)
              TSURFX=TSURFS*XAREA
              TSURWY=TSURFB(N,L,NZ)*CFW
              TSURWZ=TSURWY*YAREA
              TSURWS=TSURWY*TAUS(L+1)
              TSURWX=TSURWS*XAREA
              !
              !     ABSORPTION OF DIRECT RADIATION BY SUNLIT LEAF SURFACES
              !
              !     STOPZ=accumulated horizontal area of intercepted direct radiation
              !
              D1700: DO M=1,JSA1
                RADSL(NZ)=RADSL(NZ)+TSURFS*RDNDIR(N,M,NZ)
                RADSW(NZ)=RADSW(NZ)+TSURWS*RDNDIW(N,M,NZ)
                RADPL(NZ)=RADPL(NZ)+TSURFS*PARDIR(N,M,NZ)
                RADPW(NZ)=RADPW(NZ)+TSURWS*PARDIW(N,M,NZ)
                STOPSZ=STOPSZ+(TSURFX+TSURWX)*BETX(N,M)
!
          !     BACKSCATTERING OF REFLECTED DIRECT RADIATION
          !
                IF(IALBS(N,M).EQ.1)THEN
                  RADS1(NZ)=RADS1(NZ)+TSURFS*RDNDIR(N,M,NZ)
                  RADW1(NZ)=RADW1(NZ)+TSURWS*RDNDIW(N,M,NZ)
                  RADP1(NZ)=RADP1(NZ)+TSURFS*PARDIR(N,M,NZ)
                  RADQ1(NZ)=RADQ1(NZ)+TSURWS*PARDIW(N,M,NZ)
                  !
                  ! FORWARD SCATTERING OF REFLECTED DIRECT RADIATION
                  !
                ELSE
                  RADS2(NZ)=RADS2(NZ)+TSURFS*RDNDIR(N,M,NZ)
                  RADW2(NZ)=RADW2(NZ)+TSURWS*RDNDIW(N,M,NZ)
                  RADP2(NZ)=RADP2(NZ)+TSURFS*PARDIR(N,M,NZ)
                  RADQ2(NZ)=RADQ2(NZ)+TSURWS*PARDIW(N,M,NZ)
                ENDIF
    !
                !     INTENSITY OF ABSORBED DIFFUSE RADIATION AT LEAF SURFACES
                !
                !     RADYN,RADYW,RAPYN,RAPYW=diffuse SW,PAR flux absorbed by leaf,stalk surf
                !     OMEGA,OMEGX=incident angle of diffuse radn at leaf,horizontal surface
!
                D1750: DO NN=1,JLA1
                  RADYN=RADYL*OMEGA(M,N,NN)*ABSR(NZ)
                  RADYW=RADYL*OMEGA(M,N,NN)*ABSRW
                  RAPYN=RAPYL*OMEGA(M,N,NN)*ABSP(NZ)
                  RAPYW=RAPYL*OMEGA(M,N,NN)*ABSPW
                  PARDIF(N,M,L,NZ)=PARDIF(N,M,L,NZ)+RAPYN
                  PAR(N,M,L,NZ)=PAR(N,M,L,NZ)+RAPYN
!
                  !     ABSORPTION OF DIFFUSE RADIATION BY SHADED LEAF SURFACES
                  !
                  !     STOPYZ=accumulated horizontal area of intercepted diffuse radiation
                  !
                  RAYSL(NZ)=RAYSL(NZ)+TSURFY*RADYN
                  RAYSW(NZ)=RAYSW(NZ)+TSURWY*RADYW
                  RAYPL(NZ)=RAYPL(NZ)+TSURFY*RAPYN
                  RAYPW(NZ)=RAYPW(NZ)+TSURWY*RAPYW
                  STOPYZ=STOPYZ+(TSURFZ+TSURWZ)*OMEGX(M,N,NN)
    !
                  !     BACKSCATTERING OF REFLECTED DIFFUSE RADIATION
                  !
                  IF(IALBY(M,N,NN).EQ.1)THEN
                    RAYS1(NZ)=RAYS1(NZ)+TSURFY*RADYN
                    RAYW1(NZ)=RAYW1(NZ)+TSURWY*RADYW
                    RAYP1(NZ)=RAYP1(NZ)+TSURFY*RAPYN
                    RAYQ1(NZ)=RAYQ1(NZ)+TSURWY*RAPYW
                    !
                    !     FORWARD SCATTERING OF REFLECTED DIFFUSE RADIATION
                    !
                  ELSE
                    RAYS2(NZ)=RAYS2(NZ)+TSURFY*RADYN
                    RAYW2(NZ)=RAYW2(NZ)+TSURWY*RADYW
                    RAYP2(NZ)=RAYP2(NZ)+TSURFY*RAPYN
                    RAYQ2(NZ)=RAYQ2(NZ)+TSURWY*RAPYW
                  ENDIF
                ENDDO D1750
              ENDDO D1700
            ENDDO D1600
          ENDDO D1500
          !
          !     ACCUMULATED INTERCEPTION BY CANOPY LAYER
          !
          !     XTAUS=interception of direct radiation in current layer
          !     STOPZ=accumulated interception of direct radiation from topmost layer
          !     TAUS=transmission of direct radiation to next lower layer
!
          IF(STOPS+STOPSZ.GT.1.0_r8)THEN
            IF(STOPSZ.GT.ZERO)THEN
              XTAUS=(1.0_r8-STOPS)/((1.0_r8-STOPS)-(1.0_r8-STOPS-STOPSZ))
            ELSE
              XTAUS=0.0_r8
            ENDIF
            TAUS(L+1)=TAUS(L+1)*XTAUS
            STOPSZ=STOPSZ*XTAUS
            D1510: DO NZ=1,NP
              RADSL(NZ)=RADSL(NZ)*XTAUS
              RADSW(NZ)=RADSW(NZ)*XTAUS
              RADPL(NZ)=RADPL(NZ)*XTAUS
              RADPW(NZ)=RADPW(NZ)*XTAUS
              RADS1(NZ)=RADS1(NZ)*XTAUS
              RADW1(NZ)=RADW1(NZ)*XTAUS
              RADP1(NZ)=RADP1(NZ)*XTAUS
              RADQ1(NZ)=RADQ1(NZ)*XTAUS
              RADS2(NZ)=RADS2(NZ)*XTAUS
              RADW2(NZ)=RADW2(NZ)*XTAUS
              RADP2(NZ)=RADP2(NZ)*XTAUS
              RADQ2(NZ)=RADQ2(NZ)*XTAUS
            ENDDO D1510
          ENDIF
!
          !     XTAUY=interception of diffuse radiation in current layer
          !     STOPYZ=accumulated interception of diffuse radiation from topmost layer
          !     TAUY=transmission of diffuse radiation to next lower layer
!
          IF(STOPY+STOPYZ.GT.1.0)THEN
            XTAUY=(1.0_r8-STOPY)/((1.0_r8-STOPY)-(1.0_r8-STOPY-STOPYZ))
            TAUY(L+1)=TAUY(L+1)*XTAUY
            STOPYZ=STOPYZ*XTAUY
            D1520: DO NZ=1,NP
              RAYSL(NZ)=RAYSL(NZ)*XTAUY
              RAYSW(NZ)=RAYSW(NZ)*XTAUY
              RAYPL(NZ)=RAYPL(NZ)*XTAUY
              RAYPW(NZ)=RAYPW(NZ)*XTAUY
              RAYS1(NZ)=RAYS1(NZ)*XTAUY
              RAYW1(NZ)=RAYW1(NZ)*XTAUY
              RAYP1(NZ)=RAYP1(NZ)*XTAUY
              RAYQ1(NZ)=RAYQ1(NZ)*XTAUY
              RAYS2(NZ)=RAYS2(NZ)*XTAUY
              RAYW2(NZ)=RAYW2(NZ)*XTAUY
              RAYP2(NZ)=RAYP2(NZ)*XTAUY
              RAYQ2(NZ)=RAYQ2(NZ)*XTAUY
              D1730: DO N=1,JLI1
                DO  M=1,JSA1
                  PARDIF(N,M,L,NZ)=PARDIF(N,M,L,NZ)*XTAUY
                  PAR(N,M,L,NZ)=PARDIR(N,M,NZ)+PARDIF(N,M,L,NZ)
                enddo
              ENDDO D1730
            ENDDO D1520
          ENDIF
          !
          !     TOTAL RADIATION ABSORBED, REFLECTED AND TRANSMITTED BY ALL PFTs
          !
          !     RADST,RADWT,RADPT,RADQT=total atmospheric SW,PAR absorbed by leaf,stalk
          !     RA1ST,RA1WT,RA1PT,RA1QT=total backscattered SW,PAR absd by leaf,stalk
          !     RA2ST,RA2WT,RA2PT,RA2QT=total fwd scattered SW,PAR absd by leaf,stalk
          !     RAFSL,RAFPL=total fwd scattered SW,PAR to next layer
          !     RABSL,RABPL=total back scattered SW,PAR to next layer
          !     SWRadByCanP,TRADC,RADP,TRADP=total atmospheric SW,PAR absbd by each,all PFT
          !     STOPS,STOPY=accumulated interception of direct,diffuse radiation
          !     TAUS,TAUY=transmission of direct,diffuse radiation to next lower layer
          !
          D1530: DO NZ=1,NP
            RADST=RADSL(NZ)+RAYSL(NZ)
            RADWT=RADSW(NZ)+RAYSW(NZ)
            RADPT=RADPL(NZ)+RAYPL(NZ)
            RADQT=RADPW(NZ)+RAYPW(NZ)
            RA1ST=RADS1(NZ)+RAYS1(NZ)
            RA1WT=RADW1(NZ)+RAYW1(NZ)
            RA1PT=RADP1(NZ)+RAYP1(NZ)
            RA1QT=RADP1(NZ)+RAYQ1(NZ)
            RA2ST=RADS2(NZ)+RAYS2(NZ)
            RA2WT=RADW2(NZ)+RAYW2(NZ)
            RA2PT=RADP2(NZ)+RAYP2(NZ)
            RA2QT=RADQ2(NZ)+RAYQ2(NZ)
            RAFSL(L)=RAFSL(L)+(RADST*TAUR(NZ)+RA2ST*ALBR(NZ)+RA2WT*ALBRW)*YAREA
            RAFPL(L)=RAFPL(L)+(RADPT*TAUP(NZ)+RA2PT*ALBP(NZ)+RA2QT*ALBPW)*YAREA
            RABSL(L)=RABSL(L)+(RA1ST*ALBR(NZ)+RA1WT*ALBRW)*YAREA
            RABPL(L)=RABPL(L)+(RA1PT*ALBP(NZ)+RA1QT*ALBPW)*YAREA
            SWRadByCanP(NZ)=SWRadByCanP(NZ)+RADST+RADWT
            PARByCanP(NZ)=PARByCanP(NZ)+RADPT+RADQT
            TRADC=TRADC+RADST+RADWT
            TRAPC=TRAPC+RADPT+RADQT
          ENDDO D1530
          STOPS=STOPS+STOPSZ
          STOPY=STOPY+STOPYZ
          TAUS(L)=1.0_r8-STOPS
          TAU0(L)=1.0_r8-TAUS(L)
          TAUY(L)=1.0_r8-STOPY
        ELSE
          RAFSL(L)=RAFSL(L+1)
          RAFPL(L)=RAFPL(L+1)
          TAUS(L)=TAUS(L+1)
          TAU0(L)=1.0_r8-TAUS(L)
          TAUY(L)=TAUY(L+1)
        ENDIF
      ENDDO D1800
      !
      !     DIRECT AND DIFFUSE RADIATION ABSORBED AT GROUND SURFACE
      !
      !     RADSG,RADYG,RAPSG,RAPYG=direct,diffuse SW,PAR at horiCanopyHeightz ground surface
      !     RADS,RAPS =solar beam direct SW,PAR flux
      !     TAUS,TAUY=transmission of direct,diffuse radiation below canopy
      !     RADYL,RAPYL=solar beam diffuse SW,PAR flux
      !     RASG,RAPG=total SW,PAR at ground surface
      !     BETAG,OMEGAG=incident solar,sky angle at ground surface
!
      RADSG=RADS*TAUS(1)
      RADYG=RADYL*TAUY(1)+RAFSL(1)
      RAPSG=RAPS*TAUS(1)
      RAPYG=RAPYL*TAUY(1)+RAFPL(1)
      RASG=ABS(BETAG)*RADSG
      RAPG=ABS(BETAG)*RAPSG
      D20: DO N=1,JSA1
        RASG=RASG+ABS(OMEGAG(N))*RADYG
        RAPG=RAPG+ABS(OMEGAG(N))*RAPYG
      ENDDO D20
      RADG=RASG*AREA3(NU)
!
      !     RADIATION REFLECTED FROM GROUND SURFACE
      !
      !     VHCPW,VLHeatCapSnowMN=current,minimum snowpack heat capacity
      !     ALBW,VcumDrySnoWE,VcumWatSnow,VcumIceSnow=snowpack albedo,snow,water,ice volume
      !     ALBG,ALBS,FSNOW=ground,soil albedo,snow cover fraction
      !     THETW1=soil surface water content
      !     RABSL,RADPL=SW,PAR backscatter from ground surface
      !     TRADG,TRAPG=SW,PAR absorbed by ground surface
!
      IF(VHCPW1.GT.VLHeatCapSnowMN)THEN
        ALBW=(0.80_r8*VcumDrySnoWE+0.30_r8*VcumIceSnow+0.06_r8*VcumWatSnow)/(VcumDrySnoWE+VcumIceSnow+VcumWatSnow)
        !the following partition differs from that used in the surface physics module
        FSNOW=AMIN1((SnowDepth/0.07_r8)**2._r8,1.0_r8)
        ALBG=FSNOW*ALBW+(1.0_r8-FSNOW)*ALBS
      ELSE
        IF(VLSoilPoreMicP(NU).GT.ZEROS2)THEN
          THETW1=AMIN1(POROS1,VLWatMicP(NU)/VLSoilMicP(NU))
        ELSE
          THETW1=0.0_r8
        ENDIF
        ALBG=AMIN1(ALBX,ALBS+AZMAX1(ALBX-THETW1))
      ENDIF
      RABSL(0)=RASG*ALBG*0.25_r8
      RABPL(0)=RAPG*ALBG*0.25_r8
      TRADG=(1.0_r8-ALBG)*RASG*AREA3(NU)
      TRAPG=(1.0_r8-ALBG)*RAPG*AREA3(NU)
!
      !     ADD RADIATION FROM SCATTERING THROUGH CANOPY LAYERS
      !
      !     RABSL,RABPL=total backscattered SW,PAR to next layer
      !     RAFSL,RAFPL=total fwd scattered SW,PAR to next layer
      !     RADYN,RADYW,RAPYN,RAPYW=leaf,stalk SW,PAR absbd fwd+back flux
      !     RAYSL,RAYSW,RAYPL,RAYPW=total leaf,stalk SW,PAR absbd fwd+back
      !     SWRadByCanP,TRADC,PARByCanP,TRADP=total SW,PAR absbd by each,all PFT
!
      RADYL=0.0_r8
      RAPYL=0.0_r8
      TAUY(0)=1.0_r8
      RAFSL(0)=0.0_r8
      RAFPL(0)=0.0_r8
      D2800: DO L=1,JC1
        IF(CanopyHeightz(L-1).GE.SnowDepth-ZERO.AND.CanopyHeightz(L-1).GE.DPTH0-ZERO)THEN
          RADYL=RADYL*TAUY(L-1)+RAFSL(L-1)+RABSL(L-1)
          RAPYL=RAPYL*TAUY(L-1)+RAFPL(L-1)+RABPL(L-1)
          RAFSL(L)=0.0
          RAFPL(L)=0.0_r8
          D2500: DO NZ=1,NP
            RAYSL(NZ)=0.0_r8
            RAYSW(NZ)=0.0_r8
            RAYPL(NZ)=0.0_r8
            RAYPW(NZ)=0.0_r8
            D2600: DO N=1,JLI1
              TSURFY=TSURF(N,L,NZ)*CFX(NZ)
              TSURWY=TSURFB(N,L,NZ)*CFW
              D2700: DO M=1,JSA1
                D2750: DO NN=1,JLA1
                  RADYN=RADYL*OMEGA(M,N,NN)*ABSR(NZ)
                  RADYW=RADYL*OMEGA(M,N,NN)*ABSRW
                  RAPYN=RAPYL*OMEGA(M,N,NN)*ABSP(NZ)
                  RAPYW=RAPYL*OMEGA(M,N,NN)*ABSPW
                  PARDIF(N,M,L,NZ)=PARDIF(N,M,L,NZ)+RAPYN
                  PAR(N,M,L,NZ)=PAR(N,M,L,NZ)+RAPYN
                  RAYSL(NZ)=RAYSL(NZ)+TSURFY*RADYN
                  RAYSW(NZ)=RAYSW(NZ)+TSURWY*RADYW
                  RAYPL(NZ)=RAYPL(NZ)+TSURFY*RAPYN
                  RAYPW(NZ)=RAYPW(NZ)+TSURWY*RAPYW
                ENDDO D2750
              ENDDO D2700
            ENDDO D2600
            RAFSL(L)=RAFSL(L)+RAYSL(NZ)*TAUR(NZ)*YAREA
            RAFPL(L)=RAFPL(L)+RAYPL(NZ)*TAUP(NZ)*YAREA
            SWRadByCanP(NZ)=SWRadByCanP(NZ)+RAYSL(NZ)+RAYSW(NZ)
            PARByCanP(NZ)=PARByCanP(NZ)+RAYPL(NZ)+RAYPW(NZ)
            TRADC=TRADC+RAYSL(NZ)+RAYSW(NZ)
            TRAPC=TRAPC+RAYPL(NZ)+RAYPW(NZ)
          ENDDO D2500
        ELSE
          RAFSL(L)=RAFSL(L-1)
          RAFPL(L)=RAFPL(L-1)
          RABSL(L)=RABSL(L-1)
          RABPL(L)=RABPL(L-1)
        ENDIF
      ENDDO D2800
!
      !     RADIATION AT GROUND SURFACE IF NO CANOPY
!
    ELSE
      RASG=ABS(BETAG)*RADS
      D120: DO N=1,JSA1
        RASG=RASG+ABS(OMEGAG(N))*RADY
      ENDDO D120
      RADG=RASG*AREA3(NU)
      D135: DO NZ=1,NP
        SWRadByCanP(NZ)=0.0_r8
        PARByCanP(NZ)=0.0_r8
      ENDDO D135
    ENDIF
!
    !     IF NO RADIATION
!
  ELSE
    RADG=0.0_r8
    D125: DO NZ=1,NP
      SWRadByCanP(NZ)=0.0_r8
      PARByCanP(NZ)=0.0_r8
    ENDDO D125
  ENDIF
  !
  !     CANOPY AND GROUND SKY FRACTIONS USED FOR BOUNDARY LAYER CALCULNS
  !
  !     FRADG=fraction of radiation received by ground surface
  !     FRADP=fraction of radiation received by each PFT canopy
  !     CanGA,CanPA=leaf+stalk area of all PFTs,each PFT
  !
  FRADG=1.0_r8
  IF(CanGA.GT.ZEROS)THEN
    FRADPT=1.0_r8-EXP(-0.65_r8*CanGA/AREA3(NU))
    D145: DO NZ=1,NP
      FracPARByCanP(NZ)=FRADPT*CanPA(NZ)/CanGA
      FRADG=FRADG-FracPARByCanP(NZ)
    ENDDO D145
  ELSE
    FRADG=1.0_r8
    D146: DO NZ=1,NP
      FracPARByCanP(NZ)=0.0_r8
    ENDDO D146
  ENDIF
  end associate
  end subroutine MultiLayerSurfaceRadiation
end module CanopyCondsMod
