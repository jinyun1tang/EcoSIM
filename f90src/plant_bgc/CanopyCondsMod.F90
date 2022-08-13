module CanopyCondsMod

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcoSimConst
  use EcoSIMConfig
  use PlantAPIData
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
  subroutine CanopyConditionModel(I,J,DPTH0s1)
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(in) :: DPTH0s1 !water+ice depth at surface

  curday=I
  call MultiLayerSurfaceRadiation(I,J,DPTH0s1)

  call DivideCanopyLayerByLAI()

  call CalcBoundaryLayerProperties(DPTH0s1)

  end subroutine CanopyConditionModel

!------------------------------------------------------------------------------------------

  subroutine CalcBoundaryLayerProperties(DPTH0s1)
  implicit none
  real(r8), intent(in) :: DPTH0s1
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
  ARLSC=ARLFCs1+ARSTCs1
  IF(ARLSC.GT.ZEROSs1.AND.ZTs1.GE.DPTHSs1-ZEROs1 &
    .AND.ZTs1.GE.DPTH0s1-ZEROs1)THEN
    ARLSG=ARLSC/AREA3s1(NUs1)
    ZX=EXP(-0.5*ARLSG)
    ZY=1.0-ZX
    ZDs1=ZTs1*AMAX1(0.0_r8,1.0_r8-2.0_r8/ARLSG*ZY)
    ZE=ZTs1*AMAX1(0.05_r8,ZX*ZY)
  ELSE
    ZDs1=0.0_r8
    ZE=0.0_r8
  ENDIF
  IF(IFLGW.EQ.1)THEN
    ZZ=Z0s1+ZTs1
  ELSE
    ZZ=AMAX1(Z0s1,ZDs1+2.0_r8)
  ENDIF
  IF(IETYPs1.GE.0)THEN
    IF(VHCPW1s1.GT.VHCPWXs1)THEN
      ZRs1=AMAX1(0.001,ZE,ZW)
    ELSE
      ZRs1=AMAX1(0.001,ZE,ZSs1)
    ENDIF
!
!     CANOPY ISOTHERMAL BOUNDARY LAYER RESISTANCE
!
!     RAB,RAM=biome canopy,minimum isothermal boundary layer resistance
!     UA=wind speed
!     RIB=canopy isothermal Richardson number
!
    RABs1=AMAX1(RAM,(LOG((ZZ-ZDs1)/ZRs1))**2/(0.168*UAs1))
    RIBs1=1.27E+08*(ZZ-ZRs1)/(UAs1**2*TKAs1)
  ELSE
    RABs1=RAM
    RIBs1=0.0
  ENDIF
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
  !
  !     DIVISION OF CANOPY INTO JC LAYERS WITH EQUAL LAI
  !
  !     ZT,ZC=heights of combined canopy,PFT canopy
  !     ZL=height to bottom of each canopy layer
  !     ARLFC,ARSTC=leaf,stalk area of combined canopy
  !     ARLFT,ARSTT=leaf,stalk area of combined canopy layer
  !
  ZTs1=0.0
  DO 9685 NZ=1,NPs1
    ZTs1=AMAX1(ZTs1,ZCs1(NZ))
9685  CONTINUE
  ZLs1(JC1)=ZTs1+0.01
  ZL1(JC1)=ZLs1(JC1)
  ZL1(0)=0.0
  ART=(ARLFCs1+ARSTCs1)/JC1
  IF(ART.GT.ZEROSs1)THEN
    DO 2765 L=JC1,2,-1
      ARL=ARLFTs1(L)+ARSTTs1(L)
      IF(ARL.GT.1.01*ART)THEN
        DZL=ZLs1(L)-ZLs1(L-1)
        ZL1(L-1)=ZLs1(L-1)+0.5*AMIN1(1.0,(ARL-ART)/ARL)*DZL
      ELSEIF(ARL.LT.0.99*ART)THEN
        ARX=ARLFTs1(L-1)+ARSTTs1(L-1)
        DZL=ZLs1(L-1)-ZLs1(L-2)
        IF(ARX.GT.ZEROSs1)THEN
          ZL1(L-1)=ZLs1(L-1)-0.5_r8*AMIN1(1.0_r8,(ART-ARL)/ARX)*DZL
        ELSE
          ZL1(L-1)=ZLs1(L-1)
        ENDIF
      ELSE
        ZL1(L-1)=ZLs1(L-1)
      ENDIF
2765  CONTINUE
    DO 2770 L=JC1,2,-1
      ZLs1(L-1)=ZL1(L-1)
!     ZLs1(L-1)=AMAX1(0.0,AMIN1(ZLs1(L)-1.0E-06
!    2,ZLs1(L-1)))
2770  CONTINUE
  ENDIF
  end subroutine DivideCanopyLayerByLAI

!------------------------------------------------------------------------------------------

  subroutine MultiLayerSurfaceRadiation(I,J,DPTH0s1)
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(in) :: DPTH0s1
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
  associate(                     &
    ALBRs1  => plt_rad%ALBRs1  , &
    ALBSs1  => plt_rad%ALBSs1  , &
    ALBXs1  => plt_rad%ALBXs1  , &
    ALBPs1  => plt_rad%ALBPs1  , &
    FRADGs1 => plt_rad%FRADGs1 , &
    GAZIs1  => plt_rad%GAZIs1  , &
    GCOSs1  => plt_rad%GCOSs1  , &
    GSINs1  => plt_rad%GSINs1  , &
    IALBYs1 => plt_rad%IALBYs1 , &
    OMEGAs1 => plt_rad%OMEGAs1 , &
    OMEGAGs1=> plt_rad%OMEGAGs1, &
    OMEGXs1 => plt_rad%OMEGXs1 , &
    RAD0s1  => plt_rad%RAD0s1  , &
    RADGs1  => plt_rad%RADGs1  , &
    RADSs1  => plt_rad%RADSs1  , &
    RAP0s1  => plt_rad%RAP0s1  , &
    RADCs1  => plt_rad%RADCs1  , &
    RAPYs1  => plt_rad%RAPYs1  , &
    RAPSs1  => plt_rad%RAPSs1  , &
    RADYs1  => plt_rad%RADYs1  , &
    RAD1s1  => plt_rad%RAD1s1  , &
    TYSINs1 => plt_rad%TYSINs1 , &
    TAU0s1  => plt_rad%TAU0s1  , &
    TAUSs1  => plt_rad%TAUSs1  , &
    ZSINs1  => plt_rad%ZSINs1    &
  )
  ARLSSs1=0.0
  DO 1135 NZ=1,NPs1
    ARLFSs1(NZ)=0.0
    DO  NB=1,NBRs1(NZ)
      DO  L=1,JC1
        if(ZLs1(L-1)/=ZLs1(L-1))&
        print*,L,ZLs1(L),ZLs1(L-1)
        IF(ZLs1(L-1).GE.DPTHSs1-ZEROs1 &
          .AND.ZLs1(L-1).GE.DPTH0s1-ZEROs1)THEN
          DO 1130 K=1,JNODS1
            ARLFSs1(NZ)=ARLFSs1(NZ)+ARLFLs1(L,K,NB,NZ)
            ARLSSs1=ARLSSs1+ARLFLs1(L,K,NB,NZ)
1130      CONTINUE
          ARLFSs1(NZ)=ARLFSs1(NZ)+ARSTKs1(L,NB,NZ)
          ARLSSs1=ARLSSs1+ARSTKs1(L,NB,NZ)
        ENDIF
      enddo
    enddo
1135  CONTINUE
  IF(SSINs1.GT.ZEROs1)THEN
    RAD0s1=RADSs1*SSINs1+RADYs1*TYSINs1
    RAP0s1=RAPSs1*SSINs1+RAPYs1*TYSINs1
  ELSE
    RADSs1=0.0
    RADYs1=0.0
    RAPSs1=0.0
    RAPYs1=0.0
    RAD0s1=0.0
    RAP0s1=0.0
  ENDIF
  TRADC=0.0
  TRAPC=0.0
  DO 1025 NZ=1,NPs1
    RADCs1(NZ)=0.0
    RADPs1(NZ)=0.0
    CFXs1(NZ)=CFs1(NZ)*(1.0_r8-0.025_r8*ARLFPs1(NZ)/AREA3s1(NUs1))
1025  CONTINUE
  !
  !     ANGLE BETWEEN SUN AND GROUND SURFACE
  !
  !     SAZI,SCOS=solar azimuth,cosine of solar angle
  !     BETAG=incident solar angle at ground surface
  !     GCOS,GSIN=cos,sin of ground surface
  !     ZNOON=hour of solar noon from weather file
  !     0.2618=pi/12 (hrs)
  IF(SSINs1.GT.ZEROs1)THEN
    SAZI=0.2618_r8*(ZNOONs1-J)+4.7124_r8
    SCOS=SQRT(1.0-SSINs1**2)
    DGAZI=COS(GAZIs1-SAZI)
    BETAG=AMAX1(0.0,AMIN1(1.0,GCOSs1*SSINs1 &
      +GSINs1*SCOS*DGAZI))
    IF(ARLSSs1.GT.0.0)THEN
      SAGL=ASIN(SSINs1)
      !
      !     ABSORBED RADIATION FROM OPTICAL PROPERTIES ENTERED IN 'READS'
      !
      !     RADSA,RADWA,RAPSA,RAPWA=SW,PAR absorbed at leaf,stalk surface
      !     perpendicular to incoming radiation
      !
      DO 1050 NZ=1,NPs1
        RADSA(NZ)=RADSs1*ABSRs1(NZ)
        RADWA(NZ)=RADSs1*ABSRW
        RAPSA(NZ)=RAPSs1*ABSPs1(NZ)
        RAPWA(NZ)=RAPSs1*ABSPW
1050  CONTINUE
      !
      !     ANGLES BETWEEN SUN OR SKY ZONES AND FOLIAR SURFACES
      !
      !     ZAZI=leaf azimuth
      !     BETA,BETX=incident angle of direct radiation at leaf,horizontal surface
      !     ZAGL=determines forward vs backscattering
      !     IALBS=flag for forward vs backscattering
      !
      DO 1100 M=1,JSA1
        ZAZI=SAZI+(M-0.5)*PICON/real(M,r8)
        DAZI=COS(ZAZI-SAZI)
        DO  N=1,JLI1
          BETY=ZCOSs1(N)*SSINs1+ZSINs1(N)*SCOS*DAZI
          BETA(N,M)=ABS(BETY)
          BETX(N,M)=BETA(N,M)/SSINs1
          IF(ZCOSs1(N).GT.SSINs1)THEN
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
          DO  NZ=1,NPs1
            RDNDIR(N,M,NZ)=RADSA(NZ)*ABS(BETA(N,M))
            RDNDIW(N,M,NZ)=RADWA(NZ)*ABS(BETA(N,M))
            PARDIR(N,M,NZ)=RAPSA(NZ)*ABS(BETA(N,M))
            PARDIW(N,M,NZ)=RAPWA(NZ)*ABS(BETA(N,M))
            DO L=1,JC1
              PARDIFs1(N,M,L,NZ)=0.0
              PARs1(N,M,L,NZ)=PARDIR(N,M,NZ)
            enddo
          enddo
        enddo
1100  CONTINUE
      XAREA=1.00/AREA3s1(NUs1)
      YAREA=0.25/AREA3s1(NUs1)
      RADYL=RADYs1
      RAPYL=RAPYs1
      TAUSs1(JC1+1)=1.0
      TAUY(JC1+1)=1.0
      RAFSL(JC1+1)=0.0
      RAFPL(JC1+1)=0.0
      STOPS=0.0
      !
      !     RESET ARRAYS OF SUNLIT AND SHADED LEAF AREAS IN DIFFERENT
      !     LAYERS AND ANGLE CLASSES
      !
      !     TSURF,TSURFB,SURF,SURFB=leaf,stalk total,PFT surface area
!
      DO 1150 NZ=1,NPs1
        DO  L=1,JC1
          DO  N=1,JLI1
            TSURF(N,L,NZ)=0.0
            TSURFB(N,L,NZ)=0.0
          enddo
        enddo
1150  CONTINUE
      DO 1200 NZ=1,NPs1
        DO  NB=1,NBRs1(NZ)
          DO  L=1,JC1
            IF(ZLs1(L-1).GT.DPTHSs1-ZEROs1.AND.ZLs1(L-1).GT.DPTH0s1-ZEROs1)THEN
              DO 1205 N=1,JLI1
                DO 1210 K=1,JNODS1
                  TSURF(N,L,NZ)=TSURF(N,L,NZ)+SURFs1(N,L,K,NB,NZ)
1210            CONTINUE
                TSURFB(N,L,NZ)=TSURFB(N,L,NZ)+SURFBs1(N,L,NB,NZ)
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
      DO 1800 L=JC1,1,-1
        IF(ZLs1(L-1).GE.DPTHSs1-ZEROs1 &
          .AND.ZLs1(L-1).GE.DPTH0s1-ZEROs1)THEN
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
          DO 1500 NZ=1,NPs1
            RADSL(NZ)=0.0
            RADSW(NZ)=0.0
            RADPL(NZ)=0.0
            RADPW(NZ)=0.0
            RAYSL(NZ)=0.0
            RAYSW(NZ)=0.0
            RAYPL(NZ)=0.0
            RAYPW(NZ)=0.0
            RADS1(NZ)=0.0
            RADW1(NZ)=0.0
            RADP1(NZ)=0.0
            RADQ1(NZ)=0.0
            RAYS1(NZ)=0.0
            RAYW1(NZ)=0.0
            RAYP1(NZ)=0.0
            RAYQ1(NZ)=0.0
            RADS2(NZ)=0.0
            RADW2(NZ)=0.0
            RADP2(NZ)=0.0
            RADQ2(NZ)=0.0
            RAYS2(NZ)=0.0
            RAYW2(NZ)=0.0
            RAYP2(NZ)=0.0
            RAYQ2(NZ)=0.0
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
            DO 1600 N=1,JLI1
              TSURFY=TSURF(N,L,NZ)*CFXs1(NZ)
              TSURFZ=TSURFY*YAREA
              TSURFS=TSURFY*TAUSs1(L+1)
              TSURFX=TSURFS*XAREA
              TSURWY=TSURFB(N,L,NZ)*CFW
              TSURWZ=TSURWY*YAREA
              TSURWS=TSURWY*TAUSs1(L+1)
              TSURWX=TSURWS*XAREA
              !
              !     ABSORPTION OF DIRECT RADIATION BY SUNLIT LEAF SURFACES
              !
              !     STOPZ=accumulated horizontal area of intercepted direct radiation
              !
              DO 1700 M=1,JSA1
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
                DO 1750 NN=1,JLA1
                  RADYN=RADYL*OMEGAs1(M,N,NN)*ABSRs1(NZ)
                  RADYW=RADYL*OMEGAs1(M,N,NN)*ABSRW
                  RAPYN=RAPYL*OMEGAs1(M,N,NN)*ABSPs1(NZ)
                  RAPYW=RAPYL*OMEGAs1(M,N,NN)*ABSPW
                  PARDIFs1(N,M,L,NZ)=PARDIFs1(N,M,L,NZ)+RAPYN
                  PARs1(N,M,L,NZ)=PARs1(N,M,L,NZ)+RAPYN
!
                  !     ABSORPTION OF DIFFUSE RADIATION BY SHADED LEAF SURFACES
                  !
                  !     STOPYZ=accumulated horizontal area of intercepted diffuse radiation
                  !
                  RAYSL(NZ)=RAYSL(NZ)+TSURFY*RADYN
                  RAYSW(NZ)=RAYSW(NZ)+TSURWY*RADYW
                  RAYPL(NZ)=RAYPL(NZ)+TSURFY*RAPYN
                  RAYPW(NZ)=RAYPW(NZ)+TSURWY*RAPYW
                  STOPYZ=STOPYZ+(TSURFZ+TSURWZ)*OMEGXs1(M,N,NN)
    !
                  !     BACKSCATTERING OF REFLECTED DIFFUSE RADIATION
                  !
                  IF(IALBYs1(M,N,NN).EQ.1)THEN
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
            IF(STOPSZ.GT.ZEROs1)THEN
              XTAUS=(1.0-STOPS)/((1.0-STOPS)-(1.0-STOPS-STOPSZ))
            ELSE
              XTAUS=0.0
            ENDIF
            TAUSs1(L+1)=TAUSs1(L+1)*XTAUS
            STOPSZ=STOPSZ*XTAUS
            DO 1510 NZ=1,NPs1
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
            DO 1520 NZ=1,NPs1
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
              DO 1730 N=1,JLI1
                DO  M=1,JSA1
                  PARDIFs1(N,M,L,NZ)=PARDIFs1(N,M,L,NZ)*XTAUY
                  PARs1(N,M,L,NZ)=PARDIR(N,M,NZ)+PARDIFs1(N,M,L,NZ)
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
          DO 1530 NZ=1,NPs1
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
            RAFSL(L)=RAFSL(L)+(RADST*TAURs1(NZ) &
              +RA2ST*ALBRs1(NZ)+RA2WT*ALBRW)*YAREA
            RAFPL(L)=RAFPL(L)+(RADPT*TAUPs1(NZ) &
              +RA2PT*ALBPs1(NZ)+RA2QT*ALBPW)*YAREA
            RABSL(L)=RABSL(L)+(RA1ST*ALBRs1(NZ)+RA1WT*ALBRW)*YAREA
            RABPL(L)=RABPL(L)+(RA1PT*ALBPs1(NZ)+RA1QT*ALBPW)*YAREA
            RADCs1(NZ)=RADCs1(NZ)+RADST+RADWT
            RADPs1(NZ)=RADPs1(NZ)+RADPT+RADQT
            TRADC=TRADC+RADST+RADWT
            TRAPC=TRAPC+RADPT+RADQT
1530      CONTINUE
          STOPS=STOPS+STOPSZ
          STOPY=STOPY+STOPYZ
          TAUSs1(L)=1.0-STOPS
          TAU0s1(L)=1.0-TAUSs1(L)
          TAUY(L)=1.0-STOPY
        ELSE
          RAFSL(L)=RAFSL(L+1)
          RAFPL(L)=RAFPL(L+1)
          TAUSs1(L)=TAUSs1(L+1)
          TAU0s1(L)=1.0-TAUSs1(L)
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
      RADSG=RADSs1*TAUSs1(1)
      RADYG=RADYL*TAUY(1)+RAFSL(1)
      RAPSG=RAPSs1*TAUSs1(1)
      RAPYG=RAPYL*TAUY(1)+RAFPL(1)
      RASG=ABS(BETAG)*RADSG
      RAPG=ABS(BETAG)*RAPSG
      DO 20 N=1,JSA1
        RASG=RASG+ABS(OMEGAGs1(N))*RADYG
        RAPG=RAPG+ABS(OMEGAGs1(N))*RAPYG
20    CONTINUE
      RADGs1=RASG*AREA3s1(NUs1)
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
      IF(VHCPW1s1.GT.VHCPWXs1)THEN
        ALBW=(0.80_r8*VOLSSs1+0.30_r8*VOLISs1+0.06_r8*VOLWSs1) &
          /(VOLSSs1+VOLISs1+VOLWSs1)
        FSNOW=AMIN1((DPTHSs1/0.07)**2,1.0_r8)
        ALBG=FSNOW*ALBW+(1.0_r8-FSNOW)*ALBSs1
      ELSE
        IF(VOLXs1(NUs1).GT.ZEROS2s1)THEN
          THETW1=AMIN1(POROS1s1,VOLWs1(NUs1)/VOLYs1(NUs1))
        ELSE
          THETW1=0.0
        ENDIF
        ALBG=AMIN1(ALBXs1,ALBSs1+AMAX1(0.0_r8,ALBXs1-THETW1))
      ENDIF
      RABSL(0)=RASG*ALBG*0.25_r8
      RABPL(0)=RAPG*ALBG*0.25_r8
      TRADG=(1.0_r8-ALBG)*RASG*AREA3s1(NUs1)
      TRAPG=(1.0_r8-ALBG)*RAPG*AREA3s1(NUs1)
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
      DO 2800 L=1,JC1
        IF(ZLs1(L-1).GE.DPTHSs1-ZEROs1 &
         .AND.ZLs1(L-1).GE.DPTH0s1-ZEROs1)THEN
          RADYL=RADYL*TAUY(L-1)+RAFSL(L-1)+RABSL(L-1)
          RAPYL=RAPYL*TAUY(L-1)+RAFPL(L-1)+RABPL(L-1)
          RAFSL(L)=0.0
          RAFPL(L)=0.0
          DO 2500 NZ=1,NPs1
            RAYSL(NZ)=0.0
            RAYSW(NZ)=0.0
            RAYPL(NZ)=0.0
            RAYPW(NZ)=0.0
            DO 2600 N=1,JLI1
              TSURFY=TSURF(N,L,NZ)*CFXs1(NZ)
              TSURWY=TSURFB(N,L,NZ)*CFW
              DO 2700 M=1,JSA1
                DO 2750 NN=1,JLA1
                  RADYN=RADYL*OMEGAs1(M,N,NN)*ABSRs1(NZ)
                  RADYW=RADYL*OMEGAs1(M,N,NN)*ABSRW
                  RAPYN=RAPYL*OMEGAs1(M,N,NN)*ABSPs1(NZ)
                  RAPYW=RAPYL*OMEGAs1(M,N,NN)*ABSPW
                  PARDIFs1(N,M,L,NZ)=PARDIFs1(N,M,L,NZ)+RAPYN
                  PARs1(N,M,L,NZ)=PARs1(N,M,L,NZ)+RAPYN
                  RAYSL(NZ)=RAYSL(NZ)+TSURFY*RADYN
                  RAYSW(NZ)=RAYSW(NZ)+TSURWY*RADYW
                  RAYPL(NZ)=RAYPL(NZ)+TSURFY*RAPYN
                  RAYPW(NZ)=RAYPW(NZ)+TSURWY*RAPYW
2750            CONTINUE
2700          CONTINUE
2600        CONTINUE
            RAFSL(L)=RAFSL(L)+RAYSL(NZ)*TAURs1(NZ)*YAREA
            RAFPL(L)=RAFPL(L)+RAYPL(NZ)*TAUPs1(NZ)*YAREA
            RADCs1(NZ)=RADCs1(NZ)+RAYSL(NZ)+RAYSW(NZ)
            RADPs1(NZ)=RADPs1(NZ)+RAYPL(NZ)+RAYPW(NZ)
            TRADC=TRADC+RAYSL(NZ)+RAYSW(NZ)
            TRAPC=TRAPC+RAYPL(NZ)+RAYPW(NZ)
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
      RASG=ABS(BETAG)*RADSs1
      DO 120 N=1,JSA1
        RASG=RASG+ABS(OMEGAGs1(N))*RADYs1
120   CONTINUE
      RADGs1=RASG*AREA3s1(NUs1)
      DO 135 NZ=1,NPs1
        RADCs1(NZ)=0.0
        RADPs1(NZ)=0.0
135   CONTINUE
    ENDIF
!
    !     IF NO RADIATION
!
  ELSE
    RADGs1=0.0
    DO 125 NZ=1,NPs1
      RADCs1(NZ)=0.0
      RADPs1(NZ)=0.0
125 CONTINUE
  ENDIF
  !
  !     CANOPY AND GROUND SKY FRACTIONS USED FOR BOUNDARY LAYER CALCULNS
  !
  !     FRADG=fraction of radiation received by ground surface
  !     FRADP=fraction of radiation received by each PFT canopy
  !     ARLSS,ARLFS=leaf+stalk area of all PFTs,each PFT
  !
  FRADGs1=1.0
  IF(ARLSSs1.GT.ZEROSs1)THEN
    FRADPT=1.0-EXP(-0.65*ARLSSs1/AREA3s1(NUs1))
    DO 145 NZ=1,NPs1
      FRADPs1(NZ)=FRADPT*ARLFSs1(NZ)/ARLSSs1
      FRADGs1=FRADGs1-FRADPs1(NZ)
145 CONTINUE
  ELSE
    FRADGs1=1.0
    DO 146 NZ=1,NPs1
      FRADPs1(NZ)=0.0
146 CONTINUE
  ENDIF
  end associate
  end subroutine MultiLayerSurfaceRadiation
end module CanopyCondsMod
