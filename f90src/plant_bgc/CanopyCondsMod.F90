module CanopyCondsMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSimConst
  use EcoSIMConfig
  use PlantAPIData
  use minimathmod, only : AZMAX1,isnan
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: CanopyConditionModel

  real(r8), parameter :: RAM=2.78E-03_r8    !minimum boundary layer resistance (h m-1)
  real(r8), parameter :: StalkAlbedo4SWRad=0.1_r8                       !stalk albedo for shortwave
  real(r8), parameter :: StalkAlbedo4PARRad=0.1_r8                      !stalk albedo for PAR
  real(r8), parameter :: StalkAbsorpty4SWRad=1.0_r8-StalkAlbedo4SWRad   !stalk absorpt for 
  real(r8), parameter :: StalkAbsorpty4PARRad=1.0_r8-StalkAlbedo4PARRad
  real(r8), parameter :: StalkClumpFactor=0.5_r8                                     !stalk clumping factor,FORGW=minimum SOC or organic soil (g Mg-1)

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
    WindSpeedAtm      => plt_site%WindSpeedAtm      , &
    ZERO    => plt_site%ZERO    , &
    AREA3   => plt_site%AREA3   , &
    KoppenClimZone  => plt_site%KoppenClimZone  , &
    SoiSurfRoughnesst0      => plt_site%SoiSurfRoughnesst0      , &
    WindMesHeight      => plt_site%WindMesHeight      , &
    ZEROS   => plt_site%ZEROS   , &
    NU      => plt_site%NU      , &
    BndlResistAboveCanG     => plt_ew%BndlResistAboveCanG       , &
    ZeroPlanDisp      => plt_ew%ZeroPlanDisp        , &
    RoughHeight      => plt_ew%RoughHeight        , &
    RIB     => plt_ew%RIB       , &
    TairK     => plt_ew%TairK       , &
    VLHeatCapSnowMin  => plt_ew%VLHeatCapSnowMin    , &
    SnowDepth   => plt_ew%SnowDepth     , &
    VLHeatCapSurfSnow  => plt_ew%VLHeatCapSurfSnow    , &
    MaxCanopyHeight_grd      => plt_morph%MaxCanopyHeight_grd     , &
    StemArea_grd   => plt_morph%StemArea_grd  , &
    CanopyLeafArea_grd   => plt_morph%CanopyLeafArea_grd    &
  )
!     CANOPY ZERO PLANE AND ROUGHNESS HEIGHTS
!
!     CanopyLeafArea_grd,StemArea_grd=leaf,stalk area of combined canopy
!     SnowDepth,DPTH0=snowpack,surface water depths
!     ZT,ZeroPlanDisp,RoughHeight=canopy,zero plane displacement,roughness height
!     ZZ=reference height for wind speed
!
  ARLSC=CanopyLeafArea_grd+StemArea_grd
  IF(ARLSC.GT.ZEROS.AND.MaxCanopyHeight_grd.GE.SnowDepth-ZERO.AND.MaxCanopyHeight_grd.GE.DPTH0-ZERO)THEN
    ARLSG=ARLSC/AREA3(NU)
    ZX=EXP(-0.5_r8*ARLSG)
    ZY=1.0_r8-ZX
    ZeroPlanDisp=MaxCanopyHeight_grd*AZMAX1(1.0_r8-2.0_r8/ARLSG*ZY)
    ZE=MaxCanopyHeight_grd*AMAX1(0.05_r8,ZX*ZY)
  ELSE
    ZeroPlanDisp=0.0_r8
    ZE=0.0_r8
  ENDIF
  IF(IFLGW.EQ.1)THEN
    ZZ=WindMesHeight+MaxCanopyHeight_grd
  ELSE
    ZZ=AMAX1(WindMesHeight,ZeroPlanDisp+2.0_r8)
  ENDIF

  IF(KoppenClimZone.GE.0)THEN
    IF(VLHeatCapSurfSnow.GT.VLHeatCapSnowMin)THEN
      RoughHeight=AMAX1(0.001_r8,ZE,ZW)
    ELSE
      RoughHeight=AMAX1(0.001_r8,ZE,SoiSurfRoughnesst0)
    ENDIF
!
!     CANOPY ISOTHERMAL BOUNDARY LAYER RESISTANCE
!
!     BndlResistAboveCanG,RAM=biome canopy,minimum isothermal boundary layer resistance
!     WindSpeedAtm=wind speed
!     RIB=canopy isothermal Richardson number
!
    BndlResistAboveCanG=AMAX1(RAM,(LOG((ZZ-ZeroPlanDisp)/RoughHeight))**2._r8/(0.168_r8*WindSpeedAtm))
    RIB=1.27E+08_r8*(ZZ-RoughHeight)/(WindSpeedAtm**2._r8*TairK)
  ELSE
    BndlResistAboveCanG=RAM
    RIB=0.0_r8
  ENDIF
  end associate
  end subroutine CalcBoundaryLayerProperties

!------------------------------------------------------------------------------------------

  subroutine DivideCanopyLayerByLAI()
  implicit none

  real(r8) :: ZL1(0:NumOfCanopyLayers1)
  real(r8) :: ART,ARL
  real(r8) :: ARX
  real(r8) :: DZL
  integer :: NZ,L
  !     begin_execution
  associate(                      &
    NP      => plt_site%NP      , &
    ZEROS   => plt_site%ZEROS   , &
    MaxCanopyHeight_grd      => plt_morph%MaxCanopyHeight_grd     , &
    CanopyHeightz_col      => plt_morph%CanopyHeightz_col     , &
    CanopyHeight_pft     => plt_morph%CanopyHeight_pft    , &
    CanopyStemA_lyr   => plt_morph%CanopyStemA_lyr  , &
    CanopyLAgrid_lyr   => plt_morph%CanopyLAgrid_lyr  , &
    StemArea_grd   => plt_morph%StemArea_grd  , &
    CanopyLeafArea_grd   => plt_morph%CanopyLeafArea_grd    &
  )
  !
  !     DIVISION OF CANOPY INTO JC LAYERS WITH EQUAL LAI
  !
  !     ZT,ZC=heights of combined canopy,PFT canopy
  !     ZL=height to bottom of each canopy layer
  !     CanopyLeafArea_grd,StemArea_grd=leaf,stalk area of combined canopy
  !     CanopyLAgrid_lyr,CanopyStemA_lyr=leaf,stalk area of combined canopy layer
  !
  MaxCanopyHeight_grd=0.0
  D9685: DO NZ=1,NP
    MaxCanopyHeight_grd=AMAX1(MaxCanopyHeight_grd,CanopyHeight_pft(NZ))
  ENDDO D9685  
  CanopyHeightz_col(NumOfCanopyLayers1)=MaxCanopyHeight_grd+0.01_r8
  ZL1(NumOfCanopyLayers1)=CanopyHeightz_col(NumOfCanopyLayers1)
  ZL1(0)=0.0_r8
  ART=(CanopyLeafArea_grd+StemArea_grd)/NumOfCanopyLayers1

  IF(ART.GT.ZEROS)THEN
    D2765: DO L=NumOfCanopyLayers1,2,-1
      ARL=CanopyLAgrid_lyr(L)+CanopyStemA_lyr(L)
      IF(ARL.GT.1.01_r8*ART)THEN
        DZL=CanopyHeightz_col(L)-CanopyHeightz_col(L-1)
        ZL1(L-1)=CanopyHeightz_col(L-1)+0.5_r8*AMIN1(1.0,(ARL-ART)/ARL)*DZL
      ELSEIF(ARL.LT.0.99_r8*ART)THEN
        ARX=CanopyLAgrid_lyr(L-1)+CanopyStemA_lyr(L-1)
        DZL=CanopyHeightz_col(L-1)-CanopyHeightz_col(L-2)
        IF(ARX.GT.ZEROS)THEN
          ZL1(L-1)=CanopyHeightz_col(L-1)-0.5_r8*AMIN1(1.0_r8,(ART-ARL)/ARX)*DZL
        ELSE
          ZL1(L-1)=CanopyHeightz_col(L-1)
        ENDIF
      ELSE
        ZL1(L-1)=CanopyHeightz_col(L-1)
      ENDIF
    ENDDO D2765
    D2770: DO L=NumOfCanopyLayers1,2,-1
      CanopyHeightz_col(L-1)=ZL1(L-1)
!     CanopyHeightz_col(L-1)=AZMAX1(AMIN1(CanopyHeightz_col(L)-ppmc
!    2,CanopyHeightz_col(L-1)))
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
  integer :: IALBS(NumOfLeafZenithSectors1,NumOfSkyAzimuthSectors1)
  real(r8) :: TAUY(0:NumOfCanopyLayers1+1)
  real(r8) :: RadPARbyLeafSurf_zsec(NumOfLeafZenithSectors1,NumOfSkyAzimuthSectors1,JP1)
  real(r8) :: RadPARbyStalkSurf_zsec(NumOfLeafZenithSectors1,NumOfSkyAzimuthSectors1,JP1)
  real(r8) :: baksRadSW2NextL(0:NumOfCanopyLayers1+1)
  real(r8) :: baksRadPAR2NextL(0:NumOfCanopyLayers1+1)
  real(r8) :: fwdsRadSW2NextL(0:NumOfCanopyLayers1+1)
  real(r8) :: fwdsRadPAR2NextL(0:NumOfCanopyLayers1+1)
  real(r8) :: DirectSWbyLeafL_pft(JP1)
  real(r8) :: DirectPARbyLeafL_pft(JP1)
  real(r8) :: DiffusSWbyStalkL_pft(JP1)
  real(r8) :: DiffusPARbyStalkL_pft(JP1)
  real(r8) :: DirectSWbyStalkL_pft(JP1)
  real(r8) :: DirectPARbyStalkL_pft(JP1)
  real(r8) :: RadSWatStalkSurf_pft(JP1)
  real(r8) :: RadSWatLeafSurf_pft(JP1)
  real(r8) :: RadPARatLeafSurf_pft(JP1)
  REAL(R8) :: RadPARatStalkSurf_pft(JP1)
  real(r8) :: baksDirectSWbyLeafL_pft(JP1),fwdsDirectSWbyLeafL_pft(JP1)
  real(r8) :: baksDirectPARbyLeafL_pft(JP1),fwdsDirectPARbyLeafL_pft(JP1)
  real(r8) :: baksDirectPARbyStalkL_pft(JP1),fwdsDirectPARbyStalkL_pft(JP1)
  real(r8) :: baksDirectSWbyStalkL_pft(JP1),fwdsDirectSWbyStalkL_pft(JP1)
  real(r8) :: DiffusSWbyLeafL_pft(JP1),DiffusPARbyLeafL_pft(JP1)
  real(r8) :: baksDiffusSWbyLeafL_pft(JP1),fwdsDiffusSWbyLeafL_pft(JP1)
  real(r8) :: baksDiffusPARbyLeafL_pft(JP1),fwdsDiffusPARbyLeafL_pft(JP1)
  real(r8) :: baksDiffusSWbyStalkL_pft(JP1),fwdsDiffusSWbyStalkL_pft(JP1)
  real(r8) :: baksDiffusPARbyStalkL_pft(JP1),fwdsDiffusPARbyStalkL_pft(JP1)
  real(r8) :: BETA(NumOfLeafZenithSectors1,NumOfSkyAzimuthSectors1)                         !sine of direct solar radiation on leaf surface, [-]
  real(r8) :: BETX(NumOfLeafZenithSectors1,NumOfSkyAzimuthSectors1)                         !sine of direct solar radiation on leaf surface/sine of direct solar radiation, [-]
  REAL(R8) :: RadSWbyLeafSurf_zsec(NumOfLeafZenithSectors1,NumOfSkyAzimuthSectors1,JP1)
  real(r8) :: RadSWbyStalkSurf_zsec(NumOfLeafZenithSectors1,NumOfSkyAzimuthSectors1,JP1)
  real(r8) :: CanopyLeafArea_zsec(NumOfLeafZenithSectors1,NumOfCanopyLayers1,JP1)
  real(r8) :: CanopyStemArea_zsec(NumOfLeafZenithSectors1,NumOfCanopyLayers1,JP1)
  real(r8) :: TRADC,TRAPC,TRADG,TRAPG
  real(r8) :: SnowpackAlbedo,GrndAlbedo
  real(r8) :: GrndIncidSolarAngle,BETY,BETZ
  real(r8) :: DGAZI,DAZI
  real(r8) :: FracGrndBySnow,FRadPARbyLeafT
  real(r8) :: RadSWDiffusL,diffusSWLeafAbsorptAzclass,diffusSWStalkAbsorptAzclass
  real(r8) :: RadSWbyLeafT
  real(r8) :: RadSWbyStalkT,RadPARbyLeafT,RadPARbyStalkT
  real(r8) :: RADSG,RADYG
  real(r8) :: PARDiffusL,diffusPARLeafAbsorptAzclass,diffusPARStalkAbsorptAzclass
  real(r8) :: RAPSG,RAPYG,RadPAR_Grnd
  real(r8) :: fracDirectRadAbsorbtCum, fracDiffusRadAbsorbtCum,STOPZ
  real(r8) :: fracDirectRadAbsorbt, fracDiffusRadAbsorbt
  real(r8) :: baksRadSWbyLeafT,baksRadSWbyStalkT,baksRadPARbyLeafT,baksRadPARbyStalkT,fwdsRadSWbyLeafT
  real(r8) :: fwdsRadSWbyStalkT,fwdsRadPARbyLeafT,fwdsRadPARbyStalkT
  REAL(R8) :: RadSW_Grnd,SolarAngle
  real(r8) :: THETW1
  real(r8) :: TSURFX,UnselfShadeLeafArea,UnselfShadeLeafAreaAzclass,TSURFS
  real(r8) :: UnselfShadeStalkArea,UnselfShadeStalkAreaAzclass,TSURWS,TSURWX
  real(r8) :: XTAU_RadCapt,XTAUY,XAREA
  real(r8) :: YAREA
  REAL(R8) :: LeafAzimuthAngle,ZAGL
  real(r8) :: SolarAzimuthAngle,CosineSolarAngle
  !     begin_execution
  !     MULTILAYER CANOPY INTERECEPTION OF DIRECT AND DIFFUSE RADIATION
  !     IN SW AND VISIBLE BANDS BY INCLINATION N, AZIMUTH M, LAYER L,
  !     NODE K, BRANCH NB, PFT NZ
  !
  !     CanopyArea_pft,CanopyArea_grd=leaf+stalk area of combined,each PFT canopy
  !     ZL=height to bottom of canopy layer
  !     SnowDepth,DPTH0=snowpack,surface water depths
  !     CanopyLeafAreaByLayer_pft,CanopyBranchStemApft_lyr=leaf,stalk areas of PFT
  !     RAD,RadPARSolarBeam_col=vertical direct+diffuse SW,PAR
  !     RADS,RADY,RAPS,PARDiffus_col=solar beam direct,diffuse SW,PAR
  !     SineSolarAngle,TotSineSkyAngles_grd=sine of solar,sky angles
  !     RadSWbyCanopy_pft,RADP=total SW,PAR absorbed by canopy
  !     ClumpFactorCurrent_pft=clumping factor for self-shading
  !
  associate(                       &
    VLHeatCapSurfSnow  => plt_ew%VLHeatCapSurfSnow , &
    VcumIceSnow   => plt_ew%VcumIceSnow  , &
    VcumDrySnoWE   => plt_ew%VcumDrySnoWE  , &
    VcumWatSnow   => plt_ew%VcumWatSnow  , &
    VLHeatCapSnowMin  => plt_ew%VLHeatCapSnowMin , &
    SnowDepth   => plt_ew%SnowDepth  , &
    CanopySWAlbedo_pft    => plt_rad%CanopySWAlbedo_pft  , &
    SoilAlbedo    => plt_rad%SoilAlbedo  , &
    SurfAlbedo_col    => plt_rad%SurfAlbedo_col  , &
    CanopyPARalbedo_pft    => plt_rad%CanopyPARalbedo_pft  , &
    FracSWRad2Grnd   => plt_rad%FracSWRad2Grnd , &
    GroundSurfAzimuth_col   => plt_rad%GroundSurfAzimuth_col , &
    CosineGrndSlope_col    => plt_rad%CosineGrndSlope_col  , &
    SineGrndSlope_col    => plt_rad%SineGrndSlope_col  , &
    IALBY   => plt_rad%IALBY , &
    OMEGA   => plt_rad%OMEGA , &
    FracRadPARbyCanopy_pft   => plt_rad%FracRadPARbyCanopy_pft , &
    PARDirect_zsec    => plt_rad%PARDirect_zsec  , &
    PARDiffus_zsec  => plt_rad%PARDiffus_zsec, &
    OMEGAG  => plt_rad%OMEGAG   , &
    OMEGX   => plt_rad%OMEGX    , &
    RadSWSolarBeam_col    => plt_rad%RadSWSolarBeam_col     , &  !MJ/hour/m2
    SWRadOnGrnd    => plt_rad%SWRadOnGrnd     , &
    RadSWDirect_col   => plt_rad%RadSWDirect_col    , &  !MJ/hour/m2
    RadPARSolarBeam_col    => plt_rad%RadPARSolarBeam_col     , &
    RadSWbyCanopy_pft    => plt_rad%RadSWbyCanopy_pft     , &
    PARDiffus_col   => plt_rad%PARDiffus_col    , & !umol /m2/s
    PARDirect_col   => plt_rad%PARDirect_col    , & !umol /m2/s
    RadSWDiffus_col   => plt_rad%RadSWDiffus_col    , &
    SineSolarAngle    => plt_rad%SineSolarAngle     , &
    TotSineSkyAngles_grd   => plt_rad%TotSineSkyAngles_grd    , &
    TAU_RadThru    => plt_rad%TAU_RadThru     , &
    TAU_RadCapt    => plt_rad%TAU_RadCapt     , &
    SineLeafAngle    => plt_rad%SineLeafAngle     , &
    CanopyPARabsorpty_pft    => plt_rad%CanopyPARabsorpty_pft     , &
    CanopySWabsorpty_pft    => plt_rad%CanopySWabsorpty_pft     , &
    TAUP    => plt_rad%TAUP     , &
    TAUR    => plt_rad%TAUR     , &
    RadPARbyCanopy_pft    => plt_rad%RadPARbyCanopy_pft     , &
    CosineLeafAngle    => plt_rad%CosineLeafAngle     , &
    ZEROS   => plt_site%ZEROS   , &
    NU      => plt_site%NU      , &
    AREA3   => plt_site%AREA3   , &
    NP      => plt_site%NP      , &
    ZERO    => plt_site%ZERO    , &
    ZEROS2  => plt_site%ZEROS2  , &
    POROS1  => plt_site%POROS1  , &
    SolarNoonHour_col  => plt_site%SolarNoonHour_col  , &
    VLSoilPoreMicP    => plt_soilchem%VLSoilPoreMicP, &
    VLSoilMicP    => plt_soilchem%VLSoilMicP, &
    VLWatMicP    => plt_soilchem%VLWatMicP, &
    ClumpFactorCurrent_pft     => plt_morph%ClumpFactorCurrent_pft    , &
    CanopyHeightz_col      => plt_morph%CanopyHeightz_col     , &
    NumOfBranches_pft     => plt_morph%NumOfBranches_pft    , &
    LeafAreaZsec_brch   => plt_morph%LeafAreaZsec_brch  , &
    StemAreaZsec_brch   => plt_morph%StemAreaZsec_brch  , &
    CanopyLeafArea_pft   => plt_morph%CanopyLeafArea_pft  , &
    CanopyLeafAreaByLayer_pft   => plt_morph%CanopyLeafAreaByLayer_pft  , &
    CanopyBranchStemApft_lyr   => plt_morph%CanopyBranchStemApft_lyr  , &
    CanopyArea_pft   => plt_morph%CanopyArea_pft  , &
    ClumpFactor     => plt_morph%ClumpFactor    , &
    CanopyArea_grd   => plt_morph%CanopyArea_grd    &
  )
  CanopyArea_grd=0.0_r8
  D1135: DO NZ=1,NP
    CanopyArea_pft(NZ)=0.0_r8
    DO  NB=1,NumOfBranches_pft(NZ)
      DO  L=1,NumOfCanopyLayers1
        
        IF(CanopyHeightz_col(L-1).GE.SnowDepth-ZERO.AND.CanopyHeightz_col(L-1).GE.DPTH0-ZERO)THEN
          !above snow depth and above water/ice surface
          D1130: DO K=1,MaxNodesPerBranch1
            CanopyArea_pft(NZ)=CanopyArea_pft(NZ)+CanopyLeafAreaByLayer_pft(L,K,NB,NZ)
            CanopyArea_grd=CanopyArea_grd+CanopyLeafAreaByLayer_pft(L,K,NB,NZ)
          ENDDO D1130
          !add stem/stalk area
          CanopyArea_pft(NZ)=CanopyArea_pft(NZ)+CanopyBranchStemApft_lyr(L,NB,NZ)
          CanopyArea_grd=CanopyArea_grd+CanopyBranchStemApft_lyr(L,NB,NZ)
        ENDIF
      enddo
    enddo
  ENDDO D1135
  IF(SineSolarAngle.GT.ZERO)THEN
    RadSWSolarBeam_col=RadSWDirect_col*SineSolarAngle+RadSWDiffus_col*TotSineSkyAngles_grd
    RadPARSolarBeam_col=PARDirect_col*SineSolarAngle+PARDiffus_col*TotSineSkyAngles_grd
  ELSE
    RadSWDirect_col=0.0_r8
    RadSWDiffus_col=0.0_r8
    PARDirect_col=0.0_r8
    PARDiffus_col=0.0_r8
    RadSWSolarBeam_col=0.0_r8
    RadPARSolarBeam_col=0.0_r8
  ENDIF
  TRADC=0.0_r8
  TRAPC=0.0_r8
  D1025: DO NZ=1,NP
    RadSWbyCanopy_pft(NZ)=0.0_r8
    RadPARbyCanopy_pft(NZ)=0.0_r8
    ClumpFactorCurrent_pft(NZ)=ClumpFactor(NZ)*(1.0_r8-0.025_r8*CanopyLeafArea_pft(NZ)/AREA3(NU))
  ENDDO D1025
  !
  !     ANGLE BETWEEN SUN AND GROUND SURFACE
  !
  !     SolarAzimuthAngle,CosineSolarAngle=solar azimuth,cosine of solar angle, radian, east be zero
  !     GrndIncidSolarAngle=incident solar angle at ground surface
  !     CosineGrndSlope_col,SineGrndSlope_col=cos,sin of ground surface
  !     SolarNoonHour_col=hour of solar noon from weather file
  !     0.2618=pi/12 (hrs), 4.7124=1.5*pi=270 degree
  IF(SineSolarAngle.GT.ZERO)THEN
    SolarAzimuthAngle=0.2618_r8*(SolarNoonHour_col-J)+4.7124_r8
    CosineSolarAngle=SQRT(1.0_r8-SineSolarAngle**2._r8)
    DGAZI=COS(GroundSurfAzimuth_col-SolarAzimuthAngle)
    GrndIncidSolarAngle=AZMAX1(AMIN1(1.0_r8,CosineGrndSlope_col*SineSolarAngle+&
      SineGrndSlope_col*CosineSolarAngle*DGAZI))
    IF(CanopyArea_grd.GT.0.0_r8)THEN
      SolarAngle=ASIN(SineSolarAngle)
      !
      !     ABSORBED RADIATION FROM OPTICAL PROPERTIES ENTERED IN 'READS'
      !
      !     RadSWatLeafSurf_pft,RadSWatStalkSurf_pft,RadPARatLeafSurf_pft,
      !     RadPARatStalkSurf_pft=SW,PAR absorbed at leaf,stalk surface
      !     perpendicular to incoming radiation
      !
      D1050: DO NZ=1,NP
        RadSWatLeafSurf_pft(NZ)=RadSWDirect_col*CanopySWabsorpty_pft(NZ)
        RadSWatStalkSurf_pft(NZ)=RadSWDirect_col*StalkAbsorpty4SWRad
        RadPARatLeafSurf_pft(NZ)=PARDirect_col*CanopyPARabsorpty_pft(NZ)
        RadPARatStalkSurf_pft(NZ)=PARDirect_col*StalkAbsorpty4PARRad
      ENDDO D1050
      !
      !     ANGLES BETWEEN SUN OR SKY ZONES AND FOLIAR SURFACES
      !
      !     LeafAzimuthAngle=leaf azimuth
      !     BETA,BETX=incident angle of direct radiation at leaf,horizontal surface
      !     ZAGL=determines forward vs backscattering
      !     IALBS=flag for forward vs backscattering
      !
      D1100: DO M=1,NumOfSkyAzimuthSectors1
        LeafAzimuthAngle=SolarAzimuthAngle+(M-0.5_r8)*PICON/real(NumOfSkyAzimuthSectors1,r8)
        DAZI=COS(LeafAzimuthAngle-SolarAzimuthAngle)
        DO  N=1,NumOfLeafZenithSectors1
          BETY=CosineLeafAngle(N)*SineSolarAngle+SineLeafAngle(N)*CosineSolarAngle*DAZI
          BETA(N,M)=ABS(BETY)
          BETX(N,M)=BETA(N,M)/SineSolarAngle
          IF(CosineLeafAngle(N).GT.SineSolarAngle)THEN
            BETZ=ACOS(BETY)
          ELSE
            BETZ=-ACOS(BETY)
          ENDIF
          IF(BETZ.GT.-PICON2h)THEN
            ZAGL=SolarAngle+2.0_r8*BETZ
          ELSE
            ZAGL=SolarAngle-2.0_r8*(PICON+BETZ)
          ENDIF
          IF(ZAGL.GT.0.0_r8.AND.ZAGL.LT.PICON)THEN
            IALBS(N,M)=1
          ELSE
            IALBS(N,M)=2
          ENDIF
!
          !     INTENSITY OF ABSORBED DIRECT RADIATION AT LEAF SURFACES
          !
          !     RadSWbyLeafSurf_zsec,RadSWbyStalkSurf_zsec,RadPARbyLeafSurf_zsec,RadPARbyStalkSurf_zsec=SW,PAR flux absorbed by leaf,stalk surfaces
          !     PAR,PARDiffus_zsec=direct,diffuse PAR flux
          !     RadSWDiffusL,PARDiffusL=solar beam diffuse SW,PAR flux
          !     RAFYL,fwdsRadPAR2NextL=forward scattered diffuse SW,PAR flux
          !     TAU_RadCapt,TAUY=fraction of direct,diffuse radiation transmitted
!
          DO  NZ=1,NP
            RadSWbyLeafSurf_zsec(N,M,NZ)=RadSWatLeafSurf_pft(NZ)*ABS(BETA(N,M))
            RadSWbyStalkSurf_zsec(N,M,NZ)=RadSWatStalkSurf_pft(NZ)*ABS(BETA(N,M))
            RadPARbyLeafSurf_zsec(N,M,NZ)=RadPARatLeafSurf_pft(NZ)*ABS(BETA(N,M))
            RadPARbyStalkSurf_zsec(N,M,NZ)=RadPARatStalkSurf_pft(NZ)*ABS(BETA(N,M))
            DO L=1,NumOfCanopyLayers1
              PARDiffus_zsec(N,M,L,NZ)=0.0_r8
              PARDirect_zsec(N,M,L,NZ)=RadPARbyLeafSurf_zsec(N,M,NZ)
            enddo
          enddo
        enddo
      ENDDO D1100
      XAREA=1.00_r8/AREA3(NU)
      YAREA=1.00_r8/(AREA3(NU)*REAL(NumOfLeafZenithSectors1,R8))
      RadSWDiffusL=RadSWDiffus_col
      PARDiffusL=PARDiffus_col
      TAU_RadCapt(NumOfCanopyLayers1+1)=1.0_r8
      TAUY(NumOfCanopyLayers1+1)=1.0_r8
      fwdsRadSW2NextL(NumOfCanopyLayers1+1)=0.0_r8
      fwdsRadPAR2NextL(NumOfCanopyLayers1+1)=0.0_r8
      fracDirectRadAbsorbtCum=0.0_r8
      !
      !     RESET ARRAYS OF SUNLIT AND SHADED LEAF AREAS IN DIFFERENT
      !     LAYERS AND ANGLE CLASSES
      !
      !     TSURF,CanopyStemArea_zsec,SURF,StemAreaZsec_brch=leaf,stalk total,PFT surface area
!
      D1150: DO NZ=1,NP
        DO  L=1,NumOfCanopyLayers1
          DO  N=1,NumOfLeafZenithSectors1
            CanopyLeafArea_zsec(N,L,NZ)=0.0_r8
            CanopyStemArea_zsec(N,L,NZ)=0.0_r8
          enddo
        enddo
      ENDDO D1150
      D1200: DO NZ=1,NP
        DO  NB=1,NumOfBranches_pft(NZ)
          DO  L=1,NumOfCanopyLayers1
            IF(CanopyHeightz_col(L-1).GT.SnowDepth-ZERO.AND.CanopyHeightz_col(L-1).GT.DPTH0-ZERO)THEN
              D1205: DO N=1,NumOfLeafZenithSectors1
                D1210: DO K=1,MaxNodesPerBranch1
                  CanopyLeafArea_zsec(N,L,NZ)=CanopyLeafArea_zsec(N,L,NZ)+LeafAreaZsec_brch(N,L,K,NB,NZ)
                ENDDO D1210
                CanopyStemArea_zsec(N,L,NZ)=CanopyStemArea_zsec(N,L,NZ)+StemAreaZsec_brch(N,L,NB,NZ)
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
      D1800: DO L=NumOfCanopyLayers1,1,-1
        IF(CanopyHeightz_col(L-1).GE.SnowDepth-ZERO.AND.CanopyHeightz_col(L-1).GE.DPTH0-ZERO)THEN
          RadSWDiffusL=RadSWDiffusL*TAUY(L+1)+fwdsRadSW2NextL(L+1)
          PARDiffusL=PARDiffusL*TAUY(L+1)+fwdsRadPAR2NextL(L+1)
          fwdsRadSW2NextL(L)=0.0_r8
          fwdsRadPAR2NextL(L)=0.0_r8
          baksRadSW2NextL(L)=0.0_r8
          baksRadPAR2NextL(L)=0.0_r8
          fracDiffusRadAbsorbtCum=0.0_r8
          fracDirectRadAbsorbt=0.0_r8
          fracDiffusRadAbsorbt=0.0_r8
    !
    !     RESET ACCUMULATORS OB ABSORBED, REFLECTED AND TRANSMITTED RADIATION
    !
    !
          D1500: DO NZ=1,NP
            DirectSWbyLeafL_pft(NZ)=0.0_r8
            DirectSWbyStalkL_pft(NZ)=0.0_r8
            DirectPARbyLeafL_pft(NZ)=0.0_r8
            DirectPARbyStalkL_pft(NZ)=0.0_r8
            DiffusSWbyLeafL_pft(NZ)=0.0_r8
            DiffusSWbyStalkL_pft(NZ)=0.0_r8
            DiffusPARbyLeafL_pft(NZ)=0.0_r8
            DiffusPARbyStalkL_pft(NZ)=0.0_r8
            baksDirectSWbyLeafL_pft(NZ)=0.0_r8
            baksDirectSWbyStalkL_pft(NZ)=0.0_r8
            baksDirectPARbyLeafL_pft(NZ)=0.0_r8
            baksDirectPARbyStalkL_pft(NZ)=0.0_r8
            baksDiffusSWbyLeafL_pft(NZ)=0.0_r8
            baksDiffusSWbyStalkL_pft(NZ)=0.0_r8
            baksDiffusPARbyLeafL_pft(NZ)=0.0_r8
            baksDiffusPARbyStalkL_pft(NZ)=0.0_r8
            fwdsDirectSWbyLeafL_pft(NZ)=0.0_r8
            fwdsDirectSWbyStalkL_pft(NZ)=0.0_r8
            fwdsDirectPARbyLeafL_pft(NZ)=0.0_r8
            fwdsDirectPARbyStalkL_pft(NZ)=0.0_r8
            fwdsDiffusSWbyLeafL_pft(NZ)=0.0_r8
            fwdsDiffusSWbyStalkL_pft(NZ)=0.0_r8
            fwdsDiffusPARbyLeafL_pft(NZ)=0.0_r8
            fwdsDiffusPARbyStalkL_pft(NZ)=0.0_r8
    !
      !     LEAF SURFACE AREA IN EACH INCLINATION CLASS N, AZIMUTH CLASS M,
      !     LAYER L AND SPECIES NZ
      !
      !     UnselfShadeLeafArea=unself-shaded leaf area
      !     UnselfShadeLeafAreaAzclass=unself-shaded leaf area m-2 in each azimuth class
      !     TSURFS=UnselfShadeLeafArea with shading from canopy layers above
      !     TSURFX=TSURFS m-2
      !     UnselfShadeStalkArea=unself-shaded stalk area
      !     UnselfShadeStalkAreaAzclass=unself-shaded stalk area m-2 in each azimuth class
      !     TSURWS=UnselfShadeStalkArea with shading from canopy layers above
      !     TSURWX=TSURWS m-2
      !
            D1600: DO N=1,NumOfLeafZenithSectors1
              UnselfShadeLeafArea=CanopyLeafArea_zsec(N,L,NZ)*ClumpFactorCurrent_pft(NZ)
              UnselfShadeLeafAreaAzclass=UnselfShadeLeafArea*YAREA
              TSURFS=UnselfShadeLeafArea*TAU_RadCapt(L+1)
              TSURFX=TSURFS*XAREA
              UnselfShadeStalkArea=CanopyStemArea_zsec(N,L,NZ)*StalkClumpFactor
              UnselfShadeStalkAreaAzclass=UnselfShadeStalkArea*YAREA
              TSURWS=UnselfShadeStalkArea*TAU_RadCapt(L+1)
              TSURWX=TSURWS*XAREA
              !
              !     ABSORPTION OF DIRECT RADIATION BY SUNLIT LEAF SURFACES
              !
              !     STOPZ=accumulated horizontal area of intercepted direct radiation
              !
              D1700: DO M=1,NumOfSkyAzimuthSectors1
                DirectSWbyLeafL_pft(NZ)=DirectSWbyLeafL_pft(NZ)+TSURFS*RadSWbyLeafSurf_zsec(N,M,NZ)
                DirectSWbyStalkL_pft(NZ)=DirectSWbyStalkL_pft(NZ)+TSURWS*RadSWbyStalkSurf_zsec(N,M,NZ)
                DirectPARbyLeafL_pft(NZ)=DirectPARbyLeafL_pft(NZ)+TSURFS*RadPARbyLeafSurf_zsec(N,M,NZ)
                DirectPARbyStalkL_pft(NZ)=DirectPARbyStalkL_pft(NZ)+TSURWS*RadPARbyStalkSurf_zsec(N,M,NZ)
                fracDirectRadAbsorbt=fracDirectRadAbsorbt+(TSURFX+TSURWX)*BETX(N,M)
!
          !     BACKSCATTERING OF REFLECTED DIRECT RADIATION
          !
                IF(IALBS(N,M).EQ.1)THEN
                  baksDirectSWbyLeafL_pft(NZ)=baksDirectSWbyLeafL_pft(NZ)+TSURFS*RadSWbyLeafSurf_zsec(N,M,NZ)
                  baksDirectSWbyStalkL_pft(NZ)=baksDirectSWbyStalkL_pft(NZ)+TSURWS*RadSWbyStalkSurf_zsec(N,M,NZ)
                  baksDirectPARbyLeafL_pft(NZ)=baksDirectPARbyLeafL_pft(NZ)+TSURFS*RadPARbyLeafSurf_zsec(N,M,NZ)
                  baksDirectPARbyStalkL_pft(NZ)=baksDirectPARbyStalkL_pft(NZ)+TSURWS*RadPARbyStalkSurf_zsec(N,M,NZ)
                  !
                  ! FORWARD SCATTERING OF REFLECTED DIRECT RADIATION
                  !
                ELSE
                  fwdsDirectSWbyLeafL_pft(NZ)=fwdsDirectSWbyLeafL_pft(NZ)+TSURFS*RadSWbyLeafSurf_zsec(N,M,NZ)
                  fwdsDirectSWbyStalkL_pft(NZ)=fwdsDirectSWbyStalkL_pft(NZ)+TSURWS*RadSWbyStalkSurf_zsec(N,M,NZ)
                  fwdsDirectPARbyLeafL_pft(NZ)=fwdsDirectPARbyLeafL_pft(NZ)+TSURFS*RadPARbyLeafSurf_zsec(N,M,NZ)
                  fwdsDirectPARbyStalkL_pft(NZ)=fwdsDirectPARbyStalkL_pft(NZ)+TSURWS*RadPARbyStalkSurf_zsec(N,M,NZ)
                ENDIF
    !
                !     INTENSITY OF ABSORBED DIFFUSE RADIATION AT LEAF SURFACES
                !
                !     diffusSWLeafAbsorptAzclass,diffusSWStalkAbsorptAzclass,diffusPARLeafAbsorptAzclass,diffusPARStalkAbsorptAzclass=diffuse SW,PAR flux absorbed by leaf,stalk surf
                !     OMEGA,OMEGX=incident angle of diffuse radn at leaf,horizontal surface
!
                D1750: DO NN=1,NumOfLeafAzimuthSectors1
                  diffusSWLeafAbsorptAzclass=RadSWDiffusL*OMEGA(M,N,NN)*CanopySWabsorpty_pft(NZ)
                  diffusSWStalkAbsorptAzclass=RadSWDiffusL*OMEGA(M,N,NN)*StalkAbsorpty4SWRad
                  diffusPARLeafAbsorptAzclass=PARDiffusL*OMEGA(M,N,NN)*CanopyPARabsorpty_pft(NZ)
                  diffusPARStalkAbsorptAzclass=PARDiffusL*OMEGA(M,N,NN)*StalkAbsorpty4PARRad
                  PARDiffus_zsec(N,M,L,NZ)=PARDiffus_zsec(N,M,L,NZ)+diffusPARLeafAbsorptAzclass
                  PARDirect_zsec(N,M,L,NZ)=PARDirect_zsec(N,M,L,NZ)+diffusPARLeafAbsorptAzclass
!
                  !     ABSORPTION OF DIFFUSE RADIATION BY SHADED LEAF SURFACES
                  !
                  !     fracDiffusRadAbsorbt=accumulated horizontal area of intercepted diffuse radiation
                  !
                  DiffusSWbyLeafL_pft(NZ)=DiffusSWbyLeafL_pft(NZ)+UnselfShadeLeafArea*diffusSWLeafAbsorptAzclass
                  DiffusSWbyStalkL_pft(NZ)=DiffusSWbyStalkL_pft(NZ)+UnselfShadeStalkArea*diffusSWStalkAbsorptAzclass
                  DiffusPARbyLeafL_pft(NZ)=DiffusPARbyLeafL_pft(NZ)+UnselfShadeLeafArea*diffusPARLeafAbsorptAzclass
                  DiffusPARbyStalkL_pft(NZ)=DiffusPARbyStalkL_pft(NZ)+UnselfShadeStalkArea*diffusPARStalkAbsorptAzclass
                  fracDiffusRadAbsorbt=fracDiffusRadAbsorbt+(UnselfShadeLeafAreaAzclass+UnselfShadeStalkAreaAzclass)*OMEGX(M,N,NN)
    !
                  !     BACKSCATTERING OF REFLECTED DIFFUSE RADIATION
                  !
                  IF(IALBY(M,N,NN).EQ.1)THEN
                    baksDiffusSWbyLeafL_pft(NZ)=baksDiffusSWbyLeafL_pft(NZ)+UnselfShadeLeafArea*diffusSWLeafAbsorptAzclass
                    baksDiffusSWbyStalkL_pft(NZ)=baksDiffusSWbyStalkL_pft(NZ)+UnselfShadeStalkArea*diffusSWStalkAbsorptAzclass
                    baksDiffusPARbyLeafL_pft(NZ)=baksDiffusPARbyLeafL_pft(NZ)+UnselfShadeLeafArea*diffusPARLeafAbsorptAzclass
                    baksDiffusPARbyStalkL_pft(NZ)=baksDiffusPARbyStalkL_pft(NZ)+UnselfShadeStalkArea*diffusPARStalkAbsorptAzclass
                    !
                    !     FORWARD SCATTERING OF REFLECTED DIFFUSE RADIATION
                    !
                  ELSE
                    fwdsDiffusSWbyLeafL_pft(NZ)=fwdsDiffusSWbyLeafL_pft(NZ)+UnselfShadeLeafArea*diffusSWLeafAbsorptAzclass
                    fwdsDiffusSWbyStalkL_pft(NZ)=fwdsDiffusSWbyStalkL_pft(NZ)+UnselfShadeStalkArea*diffusSWStalkAbsorptAzclass
                    fwdsDiffusPARbyLeafL_pft(NZ)=fwdsDiffusPARbyLeafL_pft(NZ)+UnselfShadeLeafArea*diffusPARLeafAbsorptAzclass
                    fwdsDiffusPARbyStalkL_pft(NZ)=fwdsDiffusPARbyStalkL_pft(NZ)+UnselfShadeStalkArea*diffusPARStalkAbsorptAzclass
                  ENDIF
                ENDDO D1750
              ENDDO D1700
            ENDDO D1600
          ENDDO D1500
          !
          !     ACCUMULATED INTERCEPTION BY CANOPY LAYER
          !
          !     XTAU_RadCapt=interception of direct radiation in current layer
          !     STOPZ=accumulated interception of direct radiation from topmost layer
          !     TAU_RadCapt=transmission of direct radiation to next lower layer
!
          IF(fracDirectRadAbsorbtCum+fracDirectRadAbsorbt.GT.1.0_r8)THEN
            IF(fracDirectRadAbsorbt.GT.ZERO)THEN
              XTAU_RadCapt=(1.0_r8-fracDirectRadAbsorbtCum)/((1.0_r8-fracDirectRadAbsorbtCum)-&
                (1.0_r8-fracDirectRadAbsorbtCum-fracDirectRadAbsorbt))
            ELSE
              XTAU_RadCapt=0.0_r8
            ENDIF
            TAU_RadCapt(L+1)=TAU_RadCapt(L+1)*XTAU_RadCapt
            fracDirectRadAbsorbt=fracDirectRadAbsorbt*XTAU_RadCapt
            D1510: DO NZ=1,NP
              DirectSWbyLeafL_pft(NZ)=DirectSWbyLeafL_pft(NZ)*XTAU_RadCapt
              DirectSWbyStalkL_pft(NZ)=DirectSWbyStalkL_pft(NZ)*XTAU_RadCapt
              DirectPARbyLeafL_pft(NZ)=DirectPARbyLeafL_pft(NZ)*XTAU_RadCapt
              DirectPARbyStalkL_pft(NZ)=DirectPARbyStalkL_pft(NZ)*XTAU_RadCapt
              baksDirectSWbyLeafL_pft(NZ)=baksDirectSWbyLeafL_pft(NZ)*XTAU_RadCapt
              baksDirectSWbyStalkL_pft(NZ)=baksDirectSWbyStalkL_pft(NZ)*XTAU_RadCapt
              baksDirectPARbyLeafL_pft(NZ)=baksDirectPARbyLeafL_pft(NZ)*XTAU_RadCapt
              baksDirectPARbyStalkL_pft(NZ)=baksDirectPARbyStalkL_pft(NZ)*XTAU_RadCapt
              fwdsDirectSWbyLeafL_pft(NZ)=fwdsDirectSWbyLeafL_pft(NZ)*XTAU_RadCapt
              fwdsDirectSWbyStalkL_pft(NZ)=fwdsDirectSWbyStalkL_pft(NZ)*XTAU_RadCapt
              fwdsDirectPARbyLeafL_pft(NZ)=fwdsDirectPARbyLeafL_pft(NZ)*XTAU_RadCapt
              fwdsDirectPARbyStalkL_pft(NZ)=fwdsDirectPARbyStalkL_pft(NZ)*XTAU_RadCapt
            ENDDO D1510
          ENDIF
!
          !     XTAUY=interception of diffuse radiation in current layer
          !     fracDiffusRadAbsorbt=accumulated interception of diffuse radiation from topmost layer
          !     TAUY=transmission of diffuse radiation to next lower layer
!
          IF(fracDiffusRadAbsorbtCum+fracDiffusRadAbsorbt.GT.1.0_r8)THEN
            XTAUY=(1.0_r8-fracDiffusRadAbsorbtCum)/((1.0_r8-fracDiffusRadAbsorbtCum)-&
              (1.0_r8-fracDiffusRadAbsorbtCum-fracDiffusRadAbsorbt))
            TAUY(L+1)=TAUY(L+1)*XTAUY
            fracDiffusRadAbsorbt=fracDiffusRadAbsorbt*XTAUY
            D1520: DO NZ=1,NP
              DiffusSWbyLeafL_pft(NZ)=DiffusSWbyLeafL_pft(NZ)*XTAUY
              DiffusSWbyStalkL_pft(NZ)=DiffusSWbyStalkL_pft(NZ)*XTAUY
              DiffusPARbyLeafL_pft(NZ)=DiffusPARbyLeafL_pft(NZ)*XTAUY
              DiffusPARbyStalkL_pft(NZ)=DiffusPARbyStalkL_pft(NZ)*XTAUY
              baksDiffusSWbyLeafL_pft(NZ)=baksDiffusSWbyLeafL_pft(NZ)*XTAUY
              baksDiffusSWbyStalkL_pft(NZ)=baksDiffusSWbyStalkL_pft(NZ)*XTAUY
              baksDiffusPARbyLeafL_pft(NZ)=baksDiffusPARbyLeafL_pft(NZ)*XTAUY
              baksDiffusPARbyStalkL_pft(NZ)=baksDiffusPARbyStalkL_pft(NZ)*XTAUY
              fwdsDiffusSWbyLeafL_pft(NZ)=fwdsDiffusSWbyLeafL_pft(NZ)*XTAUY
              fwdsDiffusSWbyStalkL_pft(NZ)=fwdsDiffusSWbyStalkL_pft(NZ)*XTAUY
              fwdsDiffusPARbyLeafL_pft(NZ)=fwdsDiffusPARbyLeafL_pft(NZ)*XTAUY
              fwdsDiffusPARbyStalkL_pft(NZ)=fwdsDiffusPARbyStalkL_pft(NZ)*XTAUY
              D1730: DO N=1,NumOfLeafZenithSectors1
                DO  M=1,NumOfSkyAzimuthSectors1
                  PARDiffus_zsec(N,M,L,NZ)=PARDiffus_zsec(N,M,L,NZ)*XTAUY
                  PARDirect_zsec(N,M,L,NZ)=RadPARbyLeafSurf_zsec(N,M,NZ)+PARDiffus_zsec(N,M,L,NZ)
                enddo
              ENDDO D1730
            ENDDO D1520
          ENDIF
          !
          !     TOTAL RADIATION ABSORBED, REFLECTED AND TRANSMITTED BY ALL PFTs
          !
          !     RadSWbyCanopy_pft,TRADC,RADP,TRADP=total atmospheric SW,PAR absbd by each,all PFT
          !     fracDirectRadAbsorbtCum,fracDiffusRadAbsorbtCum=accumulated interception of direct,diffuse radiation
          !     TAU_RadCapt,TAUY=transmission of direct,diffuse radiation to next lower layer
          !
          D1530: DO NZ=1,NP
            RadSWbyLeafT=DirectSWbyLeafL_pft(NZ)+DiffusSWbyLeafL_pft(NZ)
            RadSWbyStalkT=DirectSWbyStalkL_pft(NZ)+DiffusSWbyStalkL_pft(NZ)
            RadPARbyLeafT=DirectPARbyLeafL_pft(NZ)+DiffusPARbyLeafL_pft(NZ)
            RadPARbyStalkT=DirectPARbyStalkL_pft(NZ)+DiffusPARbyStalkL_pft(NZ)
            baksRadSWbyLeafT=baksDirectSWbyLeafL_pft(NZ)+baksDiffusSWbyLeafL_pft(NZ)
            baksRadSWbyStalkT=baksDirectSWbyStalkL_pft(NZ)+baksDiffusSWbyStalkL_pft(NZ)
            baksRadPARbyLeafT=baksDirectPARbyLeafL_pft(NZ)+baksDiffusPARbyLeafL_pft(NZ)
            baksRadPARbyStalkT=baksDirectPARbyStalkL_pft(NZ)+baksDiffusPARbyStalkL_pft(NZ)
            fwdsRadSWbyLeafT=fwdsDirectSWbyLeafL_pft(NZ)+fwdsDiffusSWbyLeafL_pft(NZ)
            fwdsRadSWbyStalkT=fwdsDirectSWbyStalkL_pft(NZ)+fwdsDiffusSWbyStalkL_pft(NZ)
            fwdsRadPARbyLeafT=fwdsDirectPARbyLeafL_pft(NZ)+fwdsDiffusPARbyLeafL_pft(NZ)
            fwdsRadPARbyStalkT=fwdsDirectPARbyStalkL_pft(NZ)+fwdsDiffusPARbyStalkL_pft(NZ)
            fwdsRadSW2NextL(L)=fwdsRadSW2NextL(L)+(RadSWbyLeafT*TAUR(NZ)+fwdsRadSWbyLeafT*CanopySWAlbedo_pft(NZ)+&
              fwdsRadSWbyStalkT*StalkAlbedo4SWRad)*YAREA
            fwdsRadPAR2NextL(L)=fwdsRadPAR2NextL(L)+(RadPARbyLeafT*TAUP(NZ)+fwdsRadPARbyLeafT*CanopyPARalbedo_pft(NZ)+&
              fwdsRadPARbyStalkT*StalkAlbedo4PARRad)*YAREA
            baksRadSW2NextL(L)=baksRadSW2NextL(L)+(baksRadSWbyLeafT*CanopySWAlbedo_pft(NZ)+ &
              baksRadSWbyStalkT*StalkAlbedo4SWRad)*YAREA
            baksRadPAR2NextL(L)=baksRadPAR2NextL(L)+(baksRadPARbyLeafT*CanopyPARalbedo_pft(NZ)+&
              baksRadPARbyStalkT*StalkAlbedo4PARRad)*YAREA
            RadSWbyCanopy_pft(NZ)=RadSWbyCanopy_pft(NZ)+RadSWbyLeafT+RadSWbyStalkT
            RadPARbyCanopy_pft(NZ)=RadPARbyCanopy_pft(NZ)+RadPARbyLeafT+RadPARbyStalkT
            TRADC=TRADC+RadSWbyLeafT+RadSWbyStalkT
            TRAPC=TRAPC+RadPARbyLeafT+RadPARbyStalkT
          ENDDO D1530
          fracDirectRadAbsorbtCum=fracDirectRadAbsorbtCum+fracDirectRadAbsorbt
          fracDiffusRadAbsorbtCum=fracDiffusRadAbsorbtCum+fracDiffusRadAbsorbt
          TAU_RadCapt(L)=1.0_r8-fracDirectRadAbsorbtCum
          TAU_RadThru(L)=1.0_r8-TAU_RadCapt(L)
          TAUY(L)=1.0_r8-fracDiffusRadAbsorbtCum
        ELSE
          fwdsRadSW2NextL(L)=fwdsRadSW2NextL(L+1)
          fwdsRadPAR2NextL(L)=fwdsRadPAR2NextL(L+1)
          TAU_RadCapt(L)=TAU_RadCapt(L+1)
          TAU_RadThru(L)=1.0_r8-TAU_RadCapt(L)
          TAUY(L)=TAUY(L+1)
        ENDIF
      ENDDO D1800
      !
      !     DIRECT AND DIFFUSE RADIATION ABSORBED AT GROUND SURFACE
      !
      !     RADSG,RADYG,RAPSG,RAPYG=direct,diffuse SW,PAR at horiCanopyHeightz_col ground surface
      !     RADS,PARDirect_col=solar beam direct SW,PAR flux
      !     TAU_RadCapt,TAUY=transmission of direct,diffuse radiation below canopy
      !     RadSWDiffusL,PARDiffusL=solar beam diffuse SW,PAR flux
      !     RadSW_Grnd,RadPAR_Grnd=total SW,PAR at ground surface
      !     GrndIncidSolarAngle,OMEGAG=incident solar,sky angle at ground surface
!
      RADSG=RadSWDirect_col*TAU_RadCapt(1)
      RADYG=RadSWDiffusL*TAUY(1)+fwdsRadSW2NextL(1)
      RAPSG=PARDirect_col*TAU_RadCapt(1)
      RAPYG=PARDiffusL*TAUY(1)+fwdsRadPAR2NextL(1)
      RadSW_Grnd=ABS(GrndIncidSolarAngle)*RADSG
      RadPAR_Grnd=ABS(GrndIncidSolarAngle)*RAPSG
      D20: DO N=1,NumOfSkyAzimuthSectors1
        RadSW_Grnd=RadSW_Grnd+ABS(OMEGAG(N))*RADYG
        RadPAR_Grnd=RadPAR_Grnd+ABS(OMEGAG(N))*RAPYG
      ENDDO D20 
      SWRadOnGrnd=RadSW_Grnd*AREA3(NU)
!
      !     RADIATION REFLECTED FROM GROUND SURFACE
      !
      !     VHCPW,VLHeatCapSnowMin=current,minimum snowpack heat capacity
      !     SnowpackAlbedo,VcumDrySnoWE,VcumWatSnow,VcumIceSnow=snowpack albedo,snow,water,ice volume
      !     GrndAlbedo,SoilAlbedo,FracGrndBySnow=ground,soil albedo,snow cover fraction
      !     THETW1=soil surface water content
      !     baksRadSW2NextL,DirectPARbyLeafL_pft=SW,PAR backscatter from ground surface
      !     TRADG,TRAPG=SW,PAR absorbed by ground surface
!
      IF(VLHeatCapSurfSnow.GT.VLHeatCapSnowMin)THEN
        SnowpackAlbedo=(0.80_r8*VcumDrySnoWE+0.30_r8*VcumIceSnow+0.06_r8*VcumWatSnow) &
          /(VcumDrySnoWE+VcumIceSnow+VcumWatSnow)
        !the following partition differs from that used in the surface physics module  
        FracGrndBySnow=AMIN1((SnowDepth/0.07_r8)**2._r8,1.0_r8)
        GrndAlbedo=FracGrndBySnow*SnowpackAlbedo+(1.0_r8-FracGrndBySnow)*SoilAlbedo
      ELSE
        IF(VLSoilPoreMicP(NU).GT.ZEROS2)THEN
          THETW1=AMIN1(POROS1,VLWatMicP(NU)/VLSoilMicP(NU))
        ELSE
          THETW1=0.0_r8
        ENDIF
        GrndAlbedo=AMIN1(SurfAlbedo_col,SoilAlbedo+AZMAX1(SurfAlbedo_col-THETW1))
      ENDIF
      baksRadSW2NextL(0)=RadSW_Grnd*GrndAlbedo*0.25_r8
      baksRadPAR2NextL(0)=RadPAR_Grnd*GrndAlbedo*0.25_r8
      TRADG=(1.0_r8-GrndAlbedo)*RadSW_Grnd*AREA3(NU)
      TRAPG=(1.0_r8-GrndAlbedo)*RadPAR_Grnd*AREA3(NU)
!
      !     ADD RADIATION FROM SCATTERING THROUGH CANOPY LAYERS
      !
      !     baksRadSW2NextL,baksRadPAR2NextL=total backscattered SW,PAR to next layer
      !     fwdsRadSW2NextL,fwdsRadPAR2NextL=total fwd scattered SW,PAR to next layer
      !     diffusSWLeafAbsorptAzclass,diffusSWStalkAbsorptAzclass,diffusPARLeafAbsorptAzclass,diffusPARStalkAbsorptAzclass=leaf,stalk SW,PAR absbd fwd+back flux
      !     DiffusSWbyLeafL_pft,DiffusSWbyStalkL_pft,DiffusPARbyLeafL_pft,DiffusPARbyStalkL_pft=total leaf,stalk SW,PAR absbd fwd+back
      !     RadSWbyCanopy_pft,TRADC,RadPARbyCanopy_pft,TRADP=total SW,PAR absbd by each,all PFT
!
      RadSWDiffusL=0.0_r8
      PARDiffusL=0.0_r8
      TAUY(0)=1.0_r8
      fwdsRadSW2NextL(0)=0.0_r8
      fwdsRadPAR2NextL(0)=0.0_r8
      D2800: DO L=1,NumOfCanopyLayers1
        !the following line shuts off radiation when it is below water 
        IF(CanopyHeightz_col(L-1).GE.SnowDepth-ZERO.AND.CanopyHeightz_col(L-1).GE.DPTH0-ZERO)THEN
          RadSWDiffusL=RadSWDiffusL*TAUY(L-1)+fwdsRadSW2NextL(L-1)+baksRadSW2NextL(L-1)
          PARDiffusL=PARDiffusL*TAUY(L-1)+fwdsRadPAR2NextL(L-1)+baksRadPAR2NextL(L-1)
          fwdsRadSW2NextL(L)=0.0
          fwdsRadPAR2NextL(L)=0.0_r8
          D2500: DO NZ=1,NP
            DiffusSWbyLeafL_pft(NZ)=0.0_r8
            DiffusSWbyStalkL_pft(NZ)=0.0_r8
            DiffusPARbyLeafL_pft(NZ)=0.0_r8
            DiffusPARbyStalkL_pft(NZ)=0.0_r8
            D2600: DO N=1,NumOfLeafZenithSectors1
              UnselfShadeLeafArea=CanopyLeafArea_zsec(N,L,NZ)*ClumpFactorCurrent_pft(NZ)
              UnselfShadeStalkArea=CanopyStemArea_zsec(N,L,NZ)*StalkClumpFactor
              D2700: DO M=1,NumOfSkyAzimuthSectors1
                D2750: DO NN=1,NumOfLeafAzimuthSectors1
                  diffusSWLeafAbsorptAzclass=RadSWDiffusL*OMEGA(M,N,NN)*CanopySWabsorpty_pft(NZ)
                  diffusSWStalkAbsorptAzclass=RadSWDiffusL*OMEGA(M,N,NN)*StalkAbsorpty4SWRad
                  diffusPARLeafAbsorptAzclass=PARDiffusL*OMEGA(M,N,NN)*CanopyPARabsorpty_pft(NZ)
                  diffusPARStalkAbsorptAzclass=PARDiffusL*OMEGA(M,N,NN)*StalkAbsorpty4PARRad
                  PARDiffus_zsec(N,M,L,NZ)=PARDiffus_zsec(N,M,L,NZ)+diffusPARLeafAbsorptAzclass
                  PARDirect_zsec(N,M,L,NZ)=PARDirect_zsec(N,M,L,NZ)+diffusPARLeafAbsorptAzclass
                  DiffusSWbyLeafL_pft(NZ)=DiffusSWbyLeafL_pft(NZ)+UnselfShadeLeafArea*diffusSWLeafAbsorptAzclass
                  DiffusSWbyStalkL_pft(NZ)=DiffusSWbyStalkL_pft(NZ)+UnselfShadeStalkArea*diffusSWStalkAbsorptAzclass
                  DiffusPARbyLeafL_pft(NZ)=DiffusPARbyLeafL_pft(NZ)+UnselfShadeLeafArea*diffusPARLeafAbsorptAzclass
                  DiffusPARbyStalkL_pft(NZ)=DiffusPARbyStalkL_pft(NZ)+UnselfShadeStalkArea*diffusPARStalkAbsorptAzclass
                ENDDO D2750
              ENDDO D2700
            ENDDO D2600
            fwdsRadSW2NextL(L)=fwdsRadSW2NextL(L)+DiffusSWbyLeafL_pft(NZ)*TAUR(NZ)*YAREA
            fwdsRadPAR2NextL(L)=fwdsRadPAR2NextL(L)+DiffusPARbyLeafL_pft(NZ)*TAUP(NZ)*YAREA
            RadSWbyCanopy_pft(NZ)=RadSWbyCanopy_pft(NZ)+DiffusSWbyLeafL_pft(NZ)+DiffusSWbyStalkL_pft(NZ)
            RadPARbyCanopy_pft(NZ)=RadPARbyCanopy_pft(NZ)+DiffusPARbyLeafL_pft(NZ)+DiffusPARbyStalkL_pft(NZ)
            TRADC=TRADC+DiffusSWbyLeafL_pft(NZ)+DiffusSWbyStalkL_pft(NZ)
            TRAPC=TRAPC+DiffusPARbyLeafL_pft(NZ)+DiffusPARbyStalkL_pft(NZ)
          ENDDO D2500
        ELSE
          fwdsRadSW2NextL(L)=fwdsRadSW2NextL(L-1)
          fwdsRadPAR2NextL(L)=fwdsRadPAR2NextL(L-1)
          baksRadSW2NextL(L)=baksRadSW2NextL(L-1)
          baksRadPAR2NextL(L)=baksRadPAR2NextL(L-1)
        ENDIF
      ENDDO D2800
!
      !     RADIATION AT GROUND SURFACE IF NO CANOPY
!
    ELSE
      RadSW_Grnd=ABS(GrndIncidSolarAngle)*RadSWDirect_col
      D120: DO N=1,NumOfSkyAzimuthSectors1
        RadSW_Grnd=RadSW_Grnd+ABS(OMEGAG(N))*RadSWDiffus_col
      ENDDO D120
      SWRadOnGrnd=RadSW_Grnd*AREA3(NU)
      D135: DO NZ=1,NP
        RadSWbyCanopy_pft(NZ)=0.0_r8
        RadPARbyCanopy_pft(NZ)=0.0_r8
      ENDDO D135
    ENDIF
!
    !     IF NO RADIATION
!
  ELSE
    SWRadOnGrnd=0.0_r8
    D125: DO NZ=1,NP
      RadSWbyCanopy_pft(NZ)=0.0_r8
      RadPARbyCanopy_pft(NZ)=0.0_r8
    ENDDO D125
  ENDIF
  !
  !     CANOPY AND GROUND SKY FRACTIONS USED FOR BOUNDARY LAYER CALCULNS
  !
  !     FracSWRad2Grnd=fraction of radiation received by ground surface
  !     FRADP=fraction of radiation received by each PFT canopy
  !     CanopyArea_grd,CanopyArea_pft=leaf+stalk area of all PFTs,each PFT
  !
  FracSWRad2Grnd=1.0_r8
  IF(CanopyArea_grd.GT.ZEROS)THEN
    FRadPARbyLeafT=1.0_r8-EXP(-0.65_r8*CanopyArea_grd/AREA3(NU))
    D145: DO NZ=1,NP
      FracRadPARbyCanopy_pft(NZ)=FRadPARbyLeafT*CanopyArea_pft(NZ)/CanopyArea_grd
      FracSWRad2Grnd=FracSWRad2Grnd-FracRadPARbyCanopy_pft(NZ)
    ENDDO D145
  ELSE
    FracSWRad2Grnd=1.0_r8
    D146: DO NZ=1,NP
      FracRadPARbyCanopy_pft(NZ)=0.0_r8
    ENDDO D146
  ENDIF
  end associate
  end subroutine MultiLayerSurfaceRadiation
end module CanopyCondsMod
