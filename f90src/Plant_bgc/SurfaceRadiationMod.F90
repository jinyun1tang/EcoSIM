module SurfaceRadiationMod

  use data_kind_mod,      only: r8 => DAT_KIND_R8
  use minimathmod,        only: AZMAX1,   isnan
  use PlantBGCPars,         only: iforward, ibackward
  use PrescribePhenolMod, only: SetCanopyProfile
  use EcoSIMCtrlMod,      only: ldo_sp_mode,ldo_radiation_test
  use EcoSimConst
  use EcoSIMConfig
  use PlantAPIData
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: CanopyConditionModel

  real(r8), parameter :: RAM                = 2.78E-03_r8    !minimum boundary layer resistance (h m-1)
  real(r8), parameter :: RadSWStalkAlbedo   = 0.1_r8                       !stalk albedo for shortwave
  real(r8), parameter :: StalkAlbedo4PARRad = 0.1_r8                       !stalk albedo for PAR
  real(r8), parameter :: StalkSWAbsorpty    = 1.0_r8-RadSWStalkAlbedo      !stalk absorpt for shortwave radiation
  real(r8), parameter :: StalkPARAbsorpty   = 1.0_r8-StalkAlbedo4PARRad
  real(r8), parameter :: StalkClumpFactor   = 0.5_r8           !stalk clumping factor, FORGW = minimum SOC or organic soil (g Mg-1)

  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine CanopyConditionModel(I,J,DepthSurfWatIce)
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(in) :: DepthSurfWatIce !water+ice depth at surface

  real(r8) :: LeafAreaZsec_lpft(NumLeafZenithSectors1,NumCanopyLayers1,JP1)
  real(r8) :: StemAreaZsec_lpft(NumLeafZenithSectors1,NumCanopyLayers1,JP1)    

  if(ldo_sp_mode)then
    !do prescribed phenolgoy mode
    call SetCanopyProfile(I,J,LeafAreaZsec_lpft,StemAreaZsec_lpft)
  else
    call DivideCanopyAreaByHeight(I,J)

    call SummaryCanopyAREA(I,J,DepthSurfWatIce,LeafAreaZsec_lpft,StemAreaZsec_lpft)
  endif
   
  call SurfaceRadiation(I,J,DepthSurfWatIce,LeafAreaZsec_lpft,StemAreaZsec_lpft)

  call CalcBoundaryLayerProperties(DepthSurfWatIce)

  end subroutine CanopyConditionModel


!----------------------------------------------------------------------------------------------------
  subroutine CalcBoundaryLayerProperties(DepthSurfWatIce)
  implicit none
  real(r8), intent(in) :: DepthSurfWatIce
  real(r8) :: ARLSC
  real(r8) :: ARLSG
  real(r8) :: ZX,ZY,ZE
  REAL(R8) :: ZZ
!     begin_execution
  associate(                                                      &
    WindSpeedAtm_col        => plt_site%WindSpeedAtm_col         ,& !input  :wind speed, [m h-1]
    WindMesureHeight_col    => plt_site%WindMesureHeight_col     ,& !input  :wind speed measurement height, [m]
    ZERO                    => plt_site%ZERO                     ,& !input  :threshold zero for numerical stability, [-]
    AREA3                   => plt_site%AREA3                    ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    KoppenClimZone          => plt_site%KoppenClimZone           ,& !input  :Koppen climate zone for the grid,[-]
    SoilSurfRoughnesst0_col => plt_site%SoilSurfRoughnesst0_col  ,& !input  :initial soil surface roughness height, [m]
    ZEROS                   => plt_site%ZEROS                    ,& !input  :threshold zero for numerical stability,[-]
    NU                      => plt_site%NU                       ,& !input  :current soil surface layer number, [-]
    TairK                   => plt_ew%TairK                      ,& !input  :air temperature, [K]
    VLHeatCapSnowMin_col    => plt_ew%VLHeatCapSnowMin_col       ,& !input  :minimum snowpack heat capacity, [MJ d-2 K-1]
    SnowDepth               => plt_ew%SnowDepth                  ,& !input  :snowpack depth, [m]
    VLHeatCapSurfSnow_col   => plt_ew%VLHeatCapSurfSnow_col      ,& !input  :snowpack heat capacity, [MJ m-3 K-1]
    CanopyHeight_col        => plt_morph%CanopyHeight_col        ,& !input  :canopy height , [m]
    StemArea_col            => plt_morph%StemArea_col            ,& !input  :grid canopy stem area, [m2 d-2]
    CanopyLeafArea_col      => plt_morph%CanopyLeafArea_col      ,& !input  :grid canopy leaf area, [m2 d-2]
    AbvCanopyBndlResist_col => plt_ew%AbvCanopyBndlResist_col    ,& !output :isothermal boundary layer resistance, [h m-1]
    ZERO4PlantDisplace_col  => plt_ew%ZERO4PlantDisplace_col     ,& !output :zero plane displacement height, [m]
    RoughHeight             => plt_ew%RoughHeight                ,& !output :canopy surface roughness height, [m]
    RIB                     => plt_ew%RIB                         & !output :Richardson number for calculating boundary layer resistance, [-]
  )
!     CANOPY ZERO PLANE AND ROUGHNESS HEIGHTS
!
!     CanopyLeafArea_col,StemArea_col=leaf,stalk area of combined canopy
!     SnowDepth,DepthSurfWatIce=snowpack,surface water depths
!     ZT,ZERO4PlantDisplace_col,RoughHeight=canopy,zero plane displacement,roughness height
!     ZZ=reference height for wind speed
!
  ARLSC=CanopyLeafArea_col+StemArea_col
  IF(ARLSC.GT.ZEROS .AND. CanopyHeight_col.GE.SnowDepth-ZERO .AND. CanopyHeight_col.GE.DepthSurfWatIce-ZERO)THEN
    ARLSG                  = ARLSC/AREA3(NU)
    ZX                     = EXP(-0.5_r8*ARLSG)
    ZY                     = 1.0_r8-ZX
    ZERO4PlantDisplace_col = CanopyHeight_col*AZMAX1(1.0_r8-2.0_r8/ARLSG*ZY)
    ZE                     = CanopyHeight_col*AMAX1(0.05_r8,ZX*ZY)
  ELSE
    ZERO4PlantDisplace_col = 0.0_r8
    ZE                     = 0.0_r8
  ENDIF
  IF(IFLGW.EQ.1)THEN
    ZZ=WindMesureHeight_col+CanopyHeight_col
  ELSE
    ZZ=AMAX1(WindMesureHeight_col,ZERO4PlantDisplace_col+2.0_r8)
  ENDIF
  
  IF(KoppenClimZone.GE.0)THEN
    IF(VLHeatCapSurfSnow_col.GT.VLHeatCapSnowMin_col)THEN
      RoughHeight=AMAX1(0.001_r8,ZE,ZW)
    ELSE
      RoughHeight=AMAX1(0.001_r8,ZE,SoilSurfRoughnesst0_col)
    ENDIF
!
!     CANOPY ISOTHERMAL BOUNDARY LAYER RESISTANCE
!
!     AbvCanopyBndlResist_col,RAM=biome canopy,minimum isothermal boundary layer resistance
!     WindSpeedAtm_col=wind speed
!     RIB=canopy isothermal Richardson number
!   here 0.168 = cvonkarman**2
    AbvCanopyBndlResist_col = AMAX1(RAM,(LOG((ZZ-ZERO4PlantDisplace_col)/RoughHeight))**2._r8/(0.168_r8*WindSpeedAtm_col))
    RIB                     = 1.27E+08_r8*(ZZ-RoughHeight)/(WindSpeedAtm_col**2._r8*TairK)

  ELSE
    AbvCanopyBndlResist_col = RAM
    RIB                     = 0.0_r8
    RoughHeight             = 0._r8
  ENDIF

  end associate
  end subroutine CalcBoundaryLayerProperties

!----------------------------------------------------------------------------------------------------
  subroutine DivideCanopyAreaByHeight(I,J)
  implicit none
  integer, intent(in) :: I,J

  real(r8) :: ZL1(0:NumCanopyLayers1)
  real(r8) :: AreaInterval,AreaL
  real(r8) :: ARX  !interval canopy area: leaf+stem
  real(r8) :: DZL  !canopy interval height 
  integer :: NZ,L,K,NB,N
  !     begin_execution
  associate(                                               &
    NP                  => plt_site%NP                    ,& !input  :current number of plant species,[-]
    ZEROS               => plt_site%ZEROS                 ,& !input  :threshold zero for numerical stability,[-]
    CanopyHeight_pft    => plt_morph%CanopyHeight_pft     ,& !input  :canopy height, [m]
    CanopyStemAareZ_col => plt_morph%CanopyStemAareZ_col  ,& !input  :total stem area, [m2 d-2]
    CanopyLeafAareZ_col => plt_morph%CanopyLeafAareZ_col  ,& !input  :total leaf area, [m2 d-2]
    StemArea_col        => plt_morph%StemArea_col         ,& !input  :grid canopy stem area, [m2 d-2]
    CanopyLeafArea_col  => plt_morph%CanopyLeafArea_col   ,& !input  :grid canopy leaf area, [m2 d-2]
    CanopyHeight_col    => plt_morph%CanopyHeight_col     ,& !inoput :canopy height , [m]
    CanopyHeightZ_col   => plt_morph%CanopyHeightZ_col     & !output :canopy layer height, [m]
  )
  !
  !     DIVISION OF CANOPY INTO NumCanopyLayers LAYERS WITH EQUAL LAI
  !
  !     ZT,ZC=heights of combined canopy,PFT canopy
  !     ZL=height to bottom of each canopy layer
  !     CanopyLeafArea_col,StemArea_col=leaf,stalk area of combined canopy
  !     CanopyLeafAareZ_col,CanopyStemAareZ_col=leaf,stalk area of combined canopy layer
  !
  if(I==6.and.J==22.and..false.)then
    !set up input
    CanopyHeight_pft(1)=0.2
    plt_morph%CanopyLeafArea_pft(1)=0.1
    CanopyLeafAareZ_col(1)=0.1
    CanopyStemAareZ_col(1)=0.05
  endif  
  CanopyHeight_col=0.0_r8
  D9685: DO NZ=1,NP
    CanopyHeight_col=AMAX1(CanopyHeight_col,CanopyHeight_pft(NZ))
  ENDDO D9685  
  CanopyHeightZ_col(NumCanopyLayers1) = CanopyHeight_col+0.01_r8
  ZL1(NumCanopyLayers1)               = CanopyHeightZ_col(NumCanopyLayers1)
  ZL1(0)                                = 0.0_r8

  !divide total are into NumCanopyLayers1, from top to bottom
  AreaInterval=(CanopyLeafArea_col+StemArea_col)/NumCanopyLayers1  
  IF(AreaInterval.GT.ZEROS)THEN
    D2765: DO L=NumCanopyLayers1,2,-1
      AreaL=CanopyLeafAareZ_col(L)+CanopyStemAareZ_col(L)

      !greater than the mean leaf area or area-interval
      IF(AreaL.GT.1.01_r8*AreaInterval)THEN
        DZL      = CanopyHeightZ_col(L)-CanopyHeightZ_col(L-1)
        ZL1(L-1) = CanopyHeightZ_col(L-1)+0.5_r8*AMIN1(1.0_r8,(AreaL-AreaInterval)/AreaL)*DZL
      ELSEIF(AreaL.LT.0.99_r8*AreaInterval)THEN
        ARX = CanopyLeafAareZ_col(L-1)+CanopyStemAareZ_col(L-1)
        DZL = CanopyHeightZ_col(L-1)-CanopyHeightZ_col(L-2)
        !layer L-1 has significant leaf+stem (canopy) area
        IF(ARX.GT.ZEROS)THEN
          ZL1(L-1)=CanopyHeightZ_col(L-1)-0.5_r8*AMIN1(1.0_r8,(AreaInterval-AreaL)/ARX)*DZL
        ELSE
          ZL1(L-1)=CanopyHeightZ_col(L-1)
        ENDIF
      ELSE
        ZL1(L-1)=CanopyHeightZ_col(L-1)
      ENDIF
    ENDDO D2765

    D2770: DO L=NumCanopyLayers1,2,-1
      CanopyHeightZ_col(L-1)=ZL1(L-1)
    ENDDO D2770
  ENDIF
  if(I==1 .and. J==18 .and. .false.)then  
    write(456,*)'single',J, plt_rad%SineSunInclAngle_col  
    write(456,*)(CanopyHeightZ_col(L),L=1,NumCanopyLayers1)
  endif

  end associate
  end subroutine DivideCanopyAreaByHeight

!----------------------------------------------------------------------------------------------------
  subroutine SummaryCanopyAREA(I,J,DepthSurfWatIce,LeafAreaZsec_lpft,StemAreaZsec_lpft)
  !
  !Description:
  !Summarize canopy leaf and steam area
  implicit none
  integer , intent(in) :: I,J
  real(r8), intent(in) :: DepthSurfWatIce   !surface water/ice thickness above soil surface  
  integer :: NZ,NB,L,K,N
  real(r8), intent(out) :: LeafAreaZsec_lpft(NumLeafZenithSectors1,NumCanopyLayers1,JP1)    !leaf area in different angle sector
  real(r8), intent(out) :: StemAreaZsec_lpft(NumLeafZenithSectors1,NumCanopyLayers1,JP1)    !stem area in different angle sector

  associate(                                                   &
    NP                    => plt_site%NP                      ,& !input  :current number of plant species,[-]
    ZERO                  => plt_site%ZERO                    ,& !input  :threshold zero for numerical stability, [-]
    SnowDepth             => plt_ew%SnowDepth                 ,& !input  :snowpack depth, [m]
    NumOfBranches_pft     => plt_morph%NumOfBranches_pft      ,& !input  :number of branches,[-]
    LeafAreaZsec_brch     => plt_morph%LeafAreaZsec_brch      ,& !input  :leaf surface area, [m2 d-2]
    CanopyLeafArea_lnode  => plt_morph%CanopyLeafArea_lnode   ,& !input  :layer/node/branch leaf area, [m2 d-2]
    StemAreaZsec_brch     => plt_morph%StemAreaZsec_brch      ,& !input  :stem surface area, [m2 d-2]
    CanopyStalkArea_lbrch => plt_morph%CanopyStalkArea_lbrch  ,& !input  :plant canopy layer branch stem area, [m2 d-2]
    CanopyHeightZ_col     => plt_morph%CanopyHeightZ_col      ,& !input  :canopy layer height, [m]
    LeafStalkArea_pft     => plt_morph%LeafStalkArea_pft      ,& !output :plant leaf+stem/stalk area, [m2 d-2]
    LeafStalkArea_col     => plt_morph%LeafStalkArea_col       & !output :stalk area of combined, each PFT canopy,[m^2 d-2]
  )
  
  if(I==6.and.J==22.and..false.)then
  CanopyLeafArea_lnode(1,1,1,1)=0.1
  CanopyStalkArea_lbrch(1,1,1)=0.05
  LeafAreaZsec_brch(:,1,1,1,1)=(/0.025,0.025,0.025,0.025/)
  StemAreaZsec_brch(:,1,1,1)=(/0.0125,0.0125,0.0125,0.0125/)
  endif
  LeafStalkArea_col=0.0_r8
  D1135: DO NZ=1,NP

    LeafStalkArea_pft(NZ)=0.0_r8
    DO  NB=1,NumOfBranches_pft(NZ)
      DO  L=1,NumCanopyLayers1    
        !above snow and water    
        IF(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO .AND. CanopyHeightZ_col(L-1).GE.DepthSurfWatIce-ZERO)THEN
          !above snow depth and above water/ice surface
          !add all nodes over a branch
          D1130: DO K=1,MaxNodesPerBranch1
            LeafStalkArea_pft(NZ) = LeafStalkArea_pft(NZ)+CanopyLeafArea_lnode(L,K,NB,NZ)
            LeafStalkArea_col     = LeafStalkArea_col+CanopyLeafArea_lnode(L,K,NB,NZ)
            
          ENDDO D1130
          !add stem/stalk area
          LeafStalkArea_pft(NZ) = LeafStalkArea_pft(NZ)+CanopyStalkArea_lbrch(L,NB,NZ)
          LeafStalkArea_col     = LeafStalkArea_col+CanopyStalkArea_lbrch(L,NB,NZ)
        ENDIF
      enddo
    enddo
    if(LeafStalkArea_pft(NZ)>1.e10)then
      write(*,*)'canopy',LeafStalkArea_pft(NZ)
      stop
    endif  
  ENDDO D1135

  !summarize branches into different angle/depth classes
  D1150: DO NZ=1,NP
    DO  L=1,NumCanopyLayers1
      DO  N=1,NumLeafZenithSectors1
        LeafAreaZsec_lpft(N,L,NZ)=0.0_r8
        StemAreaZsec_lpft(N,L,NZ)=0.0_r8
      enddo
    enddo
  ENDDO D1150

  D1200: DO NZ=1,NP
    DO  NB=1,NumOfBranches_pft(NZ)
      DO  L=1,NumCanopyLayers1
        IF(CanopyHeightZ_col(L-1).GT.SnowDepth-ZERO .AND. CanopyHeightZ_col(L-1).GT.DepthSurfWatIce-ZERO)THEN
          D1205: DO N=1,NumLeafZenithSectors1
            D1210: DO K=1,MaxNodesPerBranch1
              LeafAreaZsec_lpft(N,L,NZ)=LeafAreaZsec_lpft(N,L,NZ)+LeafAreaZsec_brch(N,L,K,NB,NZ)
            ENDDO D1210
            StemAreaZsec_lpft(N,L,NZ)=StemAreaZsec_lpft(N,L,NZ)+StemAreaZsec_brch(N,L,NB,NZ)
          ENDDO D1205
        ENDIF
      enddo
    enddo
  ENDDO D1200
  end associate
  end subroutine SummaryCanopyArea

!----------------------------------------------------------------------------------------------------
  subroutine SurfaceRadiation(I,J,DepthSurfWatIce,LeafAreaZsec_lpft,StemAreaZsec_lpft)
  !
  !
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(in) :: DepthSurfWatIce   !surface water/ice thickness above soil surface    
  real(r8), intent(in) :: LeafAreaZsec_lpft(NumLeafZenithSectors1,NumCanopyLayers1,JP1)
  real(r8), intent(in) :: StemAreaZsec_lpft(NumLeafZenithSectors1,NumCanopyLayers1,JP1)  

  real(r8) :: DGAZI
  real(r8) :: SolarAzimuthAngle,CosineSunInclAngle,GrndIncidSolarAngle 
  integer  :: NZ,N
  real(r8) :: FRadPARbyLeafT,RadSW_Grnd
  !     begin_execution
  associate(                                                 &
    ZEROS                 => plt_site%ZEROS                 ,& !input  :threshold zero for numerical stability,[-]
    ZERO                  => plt_site%ZERO                  ,& !input  :threshold zero for numerical stability, [-]
    NP                    => plt_site%NP                    ,& !input  :current number of plant species,[-]
    NU                    => plt_site%NU                    ,& !input  :current soil surface layer number, [-]
    AREA3                 => plt_site%AREA3                 ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    CosineGrndSlope_col   => plt_rad%CosineGrndSlope_col    ,& !input  :cosine of slope, [-]
    SineGrndSlope_col     => plt_rad%SineGrndSlope_col      ,& !input  :sine of slope, [-]
    GroundSurfAzimuth_col => plt_rad%GroundSurfAzimuth_col  ,& !input  :azimuth of slope, [-]
    TotSineSkyAngles_grd  => plt_rad%TotSineSkyAngles_grd   ,& !input  :sine of sky angles,[-]
    OMEGAG                => plt_rad%OMEGAG                 ,& !input  :sine of solar beam on leaf surface, [-]
    SolarNoonHour_col     => plt_site%SolarNoonHour_col     ,& !input  :time of solar noon, [h]
    SineSunInclAngle_col  => plt_rad%SineSunInclAngle_col   ,& !input  :sine of solar angle, [-]
    CanopyLeafArea_pft    => plt_morph%CanopyLeafArea_pft   ,& !input  :plant canopy leaf area, [m2 d-2]
    LeafStalkArea_pft     => plt_morph%LeafStalkArea_pft    ,& !input  :plant leaf+stem/stalk area, [m2 d-2]
    ClumpFactor_pft       => plt_morph%ClumpFactor_pft      ,& !input  :clumping factor for self-shading in canopy layer, [-]
    LeafStalkArea_col     => plt_morph%LeafStalkArea_col    ,& !input  :stalk area of combined, each PFT canopy,[m^2 d-2]
    RadSWDirect_col       => plt_rad%RadSWDirect_col        ,& !inoput :direct shortwave radiation, [W m-2]
    RadSWDiffus_col       => plt_rad%RadSWDiffus_col        ,& !inoput :diffuse shortwave radiation, [W m-2]
    RadPARDirect_col      => plt_rad%RadPARDirect_col       ,& !inoput :direct PAR, [umol m-2 s-1]
    RadPARDiffus_col      => plt_rad%RadPARDiffus_col       ,& !inoput :diffuse PAR, [umol m-2 s-1]
    FracSWRad2Grnd_col    => plt_rad%FracSWRad2Grnd_col     ,& !inoput :fraction of radiation intercepted by ground surface, [-]
    RadSWbyCanopy_pft     => plt_rad%RadSWbyCanopy_pft      ,& !output :canopy absorbed shortwave radiation, [MJ d-2 h-1]
    RadSWSolarBeam_col    => plt_rad%RadSWSolarBeam_col     ,& !output :shortwave radiation in solar beam, [MJ m-2 h-1]
    RadPARSolarBeam_col   => plt_rad%RadPARSolarBeam_col    ,& !output :PAR radiation in solar beam, [umol m-2 s-1]
    RadPARbyCanopy_pft    => plt_rad%RadPARbyCanopy_pft     ,& !output :canopy absorbed PAR, [umol m-2 s-1]
    FracPARads2Canopy_pft => plt_rad%FracPARads2Canopy_pft  ,& !output :fraction of incoming PAR absorbed by canopy, [-]
    RadSWGrnd_col         => plt_rad%RadSWGrnd_col          ,& !output :radiation intercepted by ground surface, [MJ m-2 h-1]
    ClumpFactorNow_pft    => plt_morph%ClumpFactorNow_pft    & !output :clumping factor for self-shading in canopy layer at current LAI, [-]
  )
  !     MULTILAYER CANOPY INTERECEPTION OF DIRECT AND DIFFUSE RADIATION
  !     IN SW AND VISIBLE BANDS BY INCLINATION N, AZIMUTH M, LAYER L,
  !     NODE K, BRANCH NB, PFT NZ
  !
  !     LeafStalkArea_pft,LeafStalkArea_col=leaf+stalk area of combined,each PFT canopy
  !     ZL=height to bottom of canopy layer
  !     SnowDepth,DepthSurfWatIce=snowpack,surface water depths
  !     CanopyLeafArea_lnode,CanopyStalkArea_lbrch=leaf,stalk areas of PFT
  !     RAD,RadPARSolarBeam_col=vertical direct+diffuse SW,PAR
  !     RADS,RADY,RAPS,RadPARDiffus_col=solar beam direct,diffuse SW,PAR
  !     SineSunInclAngle_col,TotSineSkyAngles_grd=sine of solar,sky angles
  !     RadSWbyCanopy_pft,RADP=total SW,PAR absorbed by canopy
  !     ClumpFactorNow_pft=clumping factor for self-shading
  !

  IF(SineSunInclAngle_col.GT.ZERO)THEN
    RadSWSolarBeam_col =RadSWDirect_col*SineSunInclAngle_col+RadSWDiffus_col*TotSineSkyAngles_grd
    RadPARSolarBeam_col=RadPARDirect_col*SineSunInclAngle_col+RadPARDiffus_col*TotSineSkyAngles_grd
  ELSE
    RadSWDirect_col     = 0.0_r8
    RadSWDiffus_col     = 0.0_r8
    RadPARDirect_col    = 0.0_r8
    RadPARDiffus_col    = 0.0_r8
    RadSWSolarBeam_col  = 0.0_r8
    RadPARSolarBeam_col = 0.0_r8
  ENDIF
  if(I==6.and.J==22.and..false.)then
  write(456,*)RadSWDirect_col,RadSWDiffus_col,SineSunInclAngle_col,TotSineSkyAngles_grd
  write(456,*)'sw, par=',RadSWSolarBeam_col,RadPARSolarBeam_col,ClumpFactor_pft(1)
  endif
  !compute leaf clumping factor
  D1025: DO NZ=1,NP
    RadSWbyCanopy_pft(NZ)  = 0.0_r8
    RadPARbyCanopy_pft(NZ) = 0.0_r8
    ClumpFactorNow_pft(NZ) = ClumpFactor_pft(NZ)*(1.0_r8-0.025_r8*CanopyLeafArea_pft(NZ)/AREA3(NU))
  ENDDO D1025
  !
  !     ANGLE BETWEEN SUN AND GROUND SURFACE
  !
  !     CosineSunInclAngle=solar azimuth,cosine of solar angle, radian, east be zero
  !     SolarAzimuthAngle uses the mathematical definition from Campbell and Norman, 1998
  !     GrndIncidSolarAngle=incident solar angle at ground surface
  !     CosineGrndSlope_col,SineGrndSlope_col=cos,sin of ground surface
  !     SolarNoonHour_col=hour of solar noon from weather file
  !     0.2618=pi/12 (hrs), 4.7124=1.5*pi=270 degree
  IF(SineSunInclAngle_col.GT.ZERO)THEN
    !the solar azimuth angle is computed according to north hemisphere,
    SolarAzimuthAngle   = PICON*(SolarNoonHour_col-J)/12._r8+1.5_r8*PICON
    CosineSunInclAngle  = SQRT(1.0_r8-SineSunInclAngle_col**2._r8)
    DGAZI               = COS(GroundSurfAzimuth_col-SolarAzimuthAngle)
    GrndIncidSolarAngle = AZMAX1(AMIN1(1.0_r8,CosineGrndSlope_col*SineSunInclAngle_col+&
      SineGrndSlope_col*CosineSunInclAngle*DGAZI))

    !when there is canopy
    IF(LeafStalkArea_col.GT.0.0_r8)THEN

      if(ldo_radiation_test)then
        RadSW_Grnd=ABS(GrndIncidSolarAngle)*RadSWDirect_col
        D121: DO N=1,NumOfSkyAzimuthSects1
          RadSW_Grnd=RadSW_Grnd+ABS(OMEGAG(N))*RadSWDiffus_col
        ENDDO D121
      else
        call MultiCanLayerRadiation(I,J,DepthSurfWatIce,LeafAreaZsec_lpft,StemAreaZsec_lpft,&
          SolarAzimuthAngle,CosineSunInclAngle,GrndIncidSolarAngle,RadSW_Grnd)
      endif
      !     RADIATION AT GROUND SURFACE IF NO CANOPY      
    ELSE
      !plug in lake radiation below
      RadSW_Grnd=ABS(GrndIncidSolarAngle)*RadSWDirect_col
      D120: DO N=1,NumOfSkyAzimuthSects1
        RadSW_Grnd=RadSW_Grnd+ABS(OMEGAG(N))*RadSWDiffus_col
      ENDDO D120

      D135: DO NZ=1,NP
        RadSWbyCanopy_pft(NZ)  = 0.0_r8
        RadPARbyCanopy_pft(NZ) = 0.0_r8
      ENDDO D135
    ENDIF
    RadSWGrnd_col=RadSW_Grnd*AREA3(NU)    
    !
    !     IF NO RADIATION
    !
  ELSE
    RadSWGrnd_col=0.0_r8
    D125: DO NZ=1,NP
      RadSWbyCanopy_pft(NZ)  = 0.0_r8
      RadPARbyCanopy_pft(NZ) = 0.0_r8
    ENDDO D125
  ENDIF
  
  if(I==6.and.J==22.and..false.)then
  write(456,*)'csw,cpar=',RadSWbyCanopy_pft(1),RadPARbyCanopy_pft(1)
  endif  
  !
  !     CANOPY AND GROUND SKY FRACTIONS USED FOR BOUNDARY LAYER CALCULNS
  !
  !     FracSWRad2Grnd_col=fraction of radiation received by ground surface
  !     FRADP=fraction of radiation received by each PFT canopy
  !     LeafStalkArea_col,LeafStalkArea_pft=leaf+stalk area of all PFTs,each PFT
  !
  FracSWRad2Grnd_col=1.0_r8
  IF(LeafStalkArea_col.GT.ZEROS)THEN
    !Beer's law
    FRadPARbyLeafT=1.0_r8-EXP(-0.65_r8*LeafStalkArea_col/AREA3(NU))
    D145: DO NZ=1,NP
      FracPARads2Canopy_pft(NZ) = FRadPARbyLeafT*LeafStalkArea_pft(NZ)/LeafStalkArea_col
      if(.not.ldo_radiation_test)FracSWRad2Grnd_col = FracSWRad2Grnd_col-FracPARads2Canopy_pft(NZ)
    ENDDO D145
  ELSE
    FracSWRad2Grnd_col=1.0_r8
    D146: DO NZ=1,NP
      FracPARads2Canopy_pft(NZ)=0.0_r8
    ENDDO D146
  ENDIF
  if(I==6.and.J==22.and..false.)then
  write(456,*)'fpar',FracPARads2Canopy_pft(1)
  endif
  end associate
  end subroutine SurfaceRadiation

!----------------------------------------------------------------------------------------------------
  subroutine MultiCanLayerRadiation(I,J,DepthSurfWatIce,LeafAreaZsec_lpft,StemAreaZsec_lpft,&
    SolarAzimuthAngle,CosineSunInclAngle,GrndIncidSolarAngle,RadSW_Grnd)
  !
  !Description:
  ! Model multiple canopy layer radiation using bidirectional reflectance distribution function.
  !Ref: Grant et al., (1989), AFM, SIMULATION OF CANOPY PHOTOSYNTHESIS IN MAIZE AND SOYBEAN.
  !
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(in) :: DepthSurfWatIce   !surface water/ice thickness above soil surface    
  real(r8), intent(in) :: LeafAreaZsec_lpft(NumLeafZenithSectors1,NumCanopyLayers1,JP1)
  real(r8), intent(in) :: StemAreaZsec_lpft(NumLeafZenithSectors1,NumCanopyLayers1,JP1)  
  real(r8), intent(in) :: SolarAzimuthAngle
  real(r8), intent(in) :: CosineSunInclAngle  
  real(r8), intent(in) :: GrndIncidSolarAngle
  real(r8), intent(out):: RadSW_Grnd
  integer :: NB,NZ,L,K,M,N,NN
  integer :: iScatteringDirect(NumLeafZenithSectors1,NumOfSkyAzimuthSects1)
  real(r8) :: TAU_DifuseRTransmit(0:NumCanopyLayers1+1)
  real(r8) :: RadPARDirLeafSurf_zsec(NumLeafZenithSectors1,NumOfSkyAzimuthSects1,JP1)
  real(r8) :: RadPARDirStalkSurf_zsec(NumLeafZenithSectors1,NumOfSkyAzimuthSects1,JP1)
  real(r8) :: RadSWBakScat2NextL(0:NumCanopyLayers1+1)
  real(r8) :: RadPARBakScat2NextL(0:NumCanopyLayers1+1)
  real(r8) :: RadSWFwdScat2NextL(0:NumCanopyLayers1+1)
  real(r8) :: RadPARFwdScat2NextL(0:NumCanopyLayers1+1)
  real(r8) :: RadDirSWbyLeaf_pft(JP1)
  real(r8) :: RadDirPARbyLeaf_pft(JP1)
  real(r8) :: RadDifSWbyStalk_pft(JP1)
  real(r8) :: RadDifPARbyStalk_pft(JP1)
  real(r8) :: RadDirSWbyStalk_pft(JP1)
  real(r8) :: RadDirPARbyStalk_pft(JP1)
  real(r8) :: RadSWbyStalkSurf_pft(JP1)
  real(r8) :: RadSWbyLeafSurf_pft(JP1)
  real(r8) :: RadPARbyLeafSurf_pft(JP1)
  REAL(R8) :: RadPARbyStalkSurf_pft(JP1)
  real(r8) :: bakScatRadDirSWbyLeaf_pft(JP1),fwdScatRadDirSWbyLeaf_pft(JP1)
  real(r8) :: bakScatRadDirPARbyLeaf_pft(JP1),fwdScatRadDirPARbyLeaf_pft(JP1)
  real(r8) :: bakScatRadDirPARbyStalk_pft(JP1),fwdScatRadDirPARbyStalk_pft(JP1)
  real(r8) :: bakScatRadDirSWbyStalk_pft(JP1),fwdScatRadDirSWbyStalk_pft(JP1)
  real(r8) :: RadDifSWbyLeaf_pft(JP1),RadDifPARbyLeaf_pft(JP1)
  real(r8) :: bakScatRadDifSWbyLeaf_pft(JP1),fwdScatRadDifSWbyLeaf_pft(JP1)
  real(r8) :: bakScatRadDifPARbyLeaf_pft(JP1),fwdScatRadDifPARbyLeaf_pft(JP1)
  real(r8) :: bakScatRadDifSWbyStalk_pft(JP1),fwdScatRadDifSWbyStalk_pft(JP1)
  real(r8) :: bakScatRadDifPARbyStalk_pft(JP1),fwdScatRadDifPARbyStalk_pft(JP1)
  real(r8) :: BETA(NumLeafZenithSectors1,NumOfSkyAzimuthSects1)                         !sine of direct solar radiation on leaf surface, [-]
  real(r8) :: BETX(NumLeafZenithSectors1,NumOfSkyAzimuthSects1)                         !sine of direct solar radiation on leaf surface/sine of direct solar radiation, [-]
  REAL(R8) :: RadSWbyLeafSurf_zsec(NumLeafZenithSectors1,NumOfSkyAzimuthSects1,JP1)     !shortwave radiation on leave in each canopy sector
  real(r8) :: RadSWbyStalkSurf_zsec(NumLeafZenithSectors1,NumOfSkyAzimuthSects1,JP1)    !shortwave radiation on stalk in each canopy sector

  real(r8) :: SnowpackAlbedo,GrndAlbedo
  real(r8) :: BETY,BETZ
  real(r8) :: DAZI
  real(r8) :: FracGrndBySnow
  real(r8) :: RadSWDiffusL,diffusSWLeafAbsorptAzclass,diffusSWStalkAbsorptAzclass
  real(r8) :: RadSWbyLeafT
  real(r8) :: RadSWbyStalkT,RadPARbyLeafT,RadPARbyStalkT
  real(r8) :: RADSG,RADYG
  real(r8) :: RadPARDiffusL,diffusPARLeafAbsorptAzclass,diffusPARStalkAbsorptAzclass
  real(r8) :: RAPSG,RAPYG,RadPAR_Grnd
  real(r8) :: FracDirRadAbsorbtCum, FracDifRadAbsorbtCum,STOPZ
  real(r8) :: FracDirRadAbsorbt, FracDifRadAbsorbt
  real(r8) :: bakScatRadSWbyLeafT,bakScatRadSWbyStalkT,bakScatRadPARbyLeafT,bakScatRadPARbyStalkT,fwdScatRadSWbyLeafT
  real(r8) :: fwdScatRadSWbyStalkT,fwdScatRadPARbyLeafT,fwdScatRadPARbyStalkT
  REAL(R8) :: SolarAngle
  real(r8) :: THETW1
  real(r8) :: LeafIntceptArea,UnselfShadeLeafArea,UnselfShadeLeafAreaAzclass,TSurfLeaf
  real(r8) :: UnselfShadeStalkArea,UnselfShadeStalkAreaAzclass,TSurfStalk,StalkIntceptArea
  real(r8) :: XTAUS,XTAUY,XAREA
  real(r8) :: YAREA
  REAL(R8) :: LeafAzimuthAngle,ZAGL


  associate(                                                    &
    VLHeatCapSurfSnow_col  => plt_ew%VLHeatCapSurfSnow_col     ,& !input  :snowpack heat capacity, [MJ m-3 K-1]
    VcumIceSnow_col        => plt_ew%VcumIceSnow_col           ,& !input  :ice volume in snowpack, [m3 d-2]
    VcumDrySnoWE_col       => plt_ew%VcumDrySnoWE_col          ,& !input  :snow volume in snowpack (water equivalent), [m3 d-2]
    VcumWatSnow_col        => plt_ew%VcumWatSnow_col           ,& !input  :water volume in snowpack, [m3 d-2]
    VLHeatCapSnowMin_col   => plt_ew%VLHeatCapSnowMin_col      ,& !input  :minimum snowpack heat capacity, [MJ d-2 K-1]
    SnowDepth              => plt_ew%SnowDepth                 ,& !input  :snowpack depth, [m]
    RadSWLeafAlbedo_pft    => plt_rad%RadSWLeafAlbedo_pft      ,& !input  :canopy shortwave albedo, [-]
    SoilAlbedo             => plt_rad%SoilAlbedo               ,& !input  :soil albedo,[-]
    SurfAlbedo_col         => plt_rad%SurfAlbedo_col           ,& !input  :Surface albedo,[-]
    CanopyPARalbedo_pft    => plt_rad%CanopyPARalbedo_pft      ,& !input  :canopy PAR albedo, [-]
    iScatteringDiffus      => plt_rad%iScatteringDiffus        ,& !input  :flag for calculating backscattering of radiation in canopy,[-]
    OMEGA                  => plt_rad%OMEGA                    ,& !input  :sine of indirect sky radiation on leaf surface,[-]
    OMEGAG                 => plt_rad%OMEGAG                   ,& !input  :sine of solar beam on leaf surface, [-]
    OMEGX                  => plt_rad%OMEGX                    ,& !input  :sine of indirect sky radiation on leaf surface/sine of indirect sky radiation,[-]
    RadSWDirect_col        => plt_rad%RadSWDirect_col          ,& !input  :direct shortwave radiation, [W m-2]
    RadSWDiffus_col        => plt_rad%RadSWDiffus_col          ,& !input  :diffuse shortwave radiation, [W m-2]
    RadPARDiffus_col       => plt_rad%RadPARDiffus_col         ,& !input  :diffuse PAR, [umol m-2 s-1]
    RadPARDirect_col       => plt_rad%RadPARDirect_col         ,& !input  :direct PAR, [umol m-2 s-1]
    SineSunInclAngle_col   => plt_rad%SineSunInclAngle_col     ,& !input  :sine of solar angle, [-]
    SineLeafAngle          => plt_rad%SineLeafAngle            ,& !input  :sine of leaf angle,[-]
    LeafPARabsorpty_pft    => plt_rad%LeafPARabsorpty_pft      ,& !input  :canopy PAR absorptivity,[-]
    LeafSWabsorpty_pft     => plt_rad%LeafSWabsorpty_pft       ,& !input  :canopy shortwave absorptivity, [-]
    RadPARLeafTransmis_pft => plt_rad%RadPARLeafTransmis_pft   ,& !input  :canopy PAR transmissivity, [-]
    RadSWLeafTransmis_pft  => plt_rad%RadSWLeafTransmis_pft    ,& !input  :canopy shortwave transmissivity, [-]
    CosineLeafAngle        => plt_rad%CosineLeafAngle          ,& !input  :cosine of leaf angle,[-]
    NU                     => plt_site%NU                      ,& !input  :current soil surface layer number, [-]
    AREA3                  => plt_site%AREA3                   ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    NP                     => plt_site%NP                      ,& !input  :current number of plant species,[-]
    ZERO                   => plt_site%ZERO                    ,& !input  :threshold zero for numerical stability, [-]
    ZEROS2                 => plt_site%ZEROS2                  ,& !input  :threshold zero for numerical stability,[-]
    POROS1                 => plt_site%POROS1                  ,& !input  :top layer soil porosity, [m3 m-3]
    VLSoilPoreMicP_vr      => plt_soilchem%VLSoilPoreMicP_vr   ,& !input  :volume of soil layer, [m3 d-2]
    VLSoilMicP_vr          => plt_soilchem%VLSoilMicP_vr       ,& !input  :total micropore volume in layer, [m3 d-2]
    VLWatMicP_vr           => plt_soilchem%VLWatMicP_vr        ,& !input  :soil micropore water content, [m3 d-2]
    ClumpFactorNow_pft     => plt_morph%ClumpFactorNow_pft     ,& !input  :clumping factor for self-shading in canopy layer at current LAI, [-]
    CanopyHeightZ_col      => plt_morph%CanopyHeightZ_col      ,& !input  :canopy layer height, [m]
    RadPAR_zsec            => plt_rad%RadPAR_zsec              ,& !inoput :direct incoming PAR, [umol m-2 s-1]
    RadDifPAR_zsec         => plt_rad%RadDifPAR_zsec           ,& !inoput :diffuse incoming PAR, [umol m-2 s-1]
    RadSWbyCanopy_pft      => plt_rad%RadSWbyCanopy_pft        ,& !inoput :canopy absorbed shortwave radiation, [MJ d-2 h-1]
    TAU_DirectRTransmit    => plt_rad%TAU_DirectRTransmit      ,& !inoput :fraction of radiation intercepted by canopy layer, [-]
    RadPARbyCanopy_pft     => plt_rad%RadPARbyCanopy_pft       ,& !inoput :canopy absorbed PAR, [umol m-2 s-1]
    TAU_RadThru            => plt_rad%TAU_RadThru               & !output :fraction of radiation transmitted by canopy layer, [-]
  )
  
  SolarAngle=ASIN(SineSunInclAngle_col)
  !
  !     ABSORBED RADIATION FROM OPTICAL PROPERTIES ENTERED IN 'READS'
  !
  !     RadSWbyLeafSurf_pft,RadSWbyStalkSurf_pft,RadPARbyLeafSurf_pft,
  !     RadPARbyStalkSurf_pft=SW,PAR absorbed at leaf,stalk surface
  !     perpendicular to incoming radiation
  !
  !Distribute SW and PAR over all pfts in the grid
  D1050: DO NZ=1,NP
    RadSWbyLeafSurf_pft(NZ)   = RadSWDirect_col*LeafSWabsorpty_pft(NZ)
    RadSWbyStalkSurf_pft(NZ)  = RadSWDirect_col*StalkSWAbsorpty
    RadPARbyLeafSurf_pft(NZ)  = RadPARDirect_col*LeafPARabsorpty_pft(NZ)
    RadPARbyStalkSurf_pft(NZ) = RadPARDirect_col*StalkPARAbsorpty
  ENDDO D1050
!  write(3333,*)I*100+J,'LeafSWabsorpty_pft(NZ)',RadSWDirect_col,LeafSWabsorpty_pft(NZ),StalkSWAbsorpty
  !distribute radiation into different leaf/canopy sector
  !     ANGLES BETWEEN SUN OR SKY ZONES AND FOLIAR SURFACES
  !
  !     LeafAzimuthAngle=leaf azimuth
  !     BETA,BETX=incident angle of direct radiation at leaf,horizontal surface
  !     ZAGL=determines forward vs backscattering
  !     iScatteringDirect=flag for forward vs backscattering
  !
  D1100: DO M=1,NumOfSkyAzimuthSects1
    LeafAzimuthAngle = SolarAzimuthAngle+(M-0.5_r8)*PICON/real(NumOfSkyAzimuthSects1,r8)
    DAZI             = COS(LeafAzimuthAngle-SolarAzimuthAngle)
    DO N = 1, NumLeafZenithSectors1
      BETY      = CosineLeafAngle(N)*SineSunInclAngle_col+SineLeafAngle(N)*CosineSunInclAngle*DAZI
      BETA(N,M) = ABS(BETY)
      BETX(N,M) = BETA(N,M)/SineSunInclAngle_col

      !compute incident angle BETZ, eq. (12)/(13)
      IF(CosineLeafAngle(N).GT.SineSunInclAngle_col)THEN
        BETZ=ACOS(BETY)
      ELSE
        BETZ=-ACOS(BETY)
      ENDIF

      IF(BETZ.GT.-PICON2h)THEN
        ZAGL=SolarAngle+2.0_r8*BETZ
      ELSE
        ZAGL=SolarAngle-2.0_r8*(PICON+BETZ)
      ENDIF

      IF(ZAGL.GT.0.0_r8 .AND. ZAGL.LT.PICON)THEN
        iScatteringDirect(N,M)=ibackward
      ELSE
        iScatteringDirect(N,M)=iforward
      ENDIF
!
        !     INTENSITY OF ABSORBED DIRECT RADIATION AT LEAF SURFACES
        !
        !     RadSWbyLeafSurf_zsec,RadSWbyStalkSurf_zsec,
        !     RadPARDirLeafSurf_zsec,RadPARDirStalkSurf_zsec=SW,PAR flux absorbed by leaf,stalk surfaces
        !     PAR,RadDifPAR_zsec=direct,diffuse PAR flux
        !     RadSWDiffusL,RadPARDiffusL=solar beam diffuse SW,PAR flux
        !     RAFYL,RadPARFwdScat2NextL=forward scattered diffuse SW,PAR flux
        !     TAU_DirectRTransmit,TAU_DifuseRTransmit=fraction of direct,diffuse radiation transmitted
!
      DO  NZ=1,NP
        RadSWbyLeafSurf_zsec(N,M,NZ)    = RadSWbyLeafSurf_pft(NZ)*ABS(BETA(N,M))
        RadSWbyStalkSurf_zsec(N,M,NZ)   = RadSWbyStalkSurf_pft(NZ)*ABS(BETA(N,M))
        RadPARDirLeafSurf_zsec(N,M,NZ)  = RadPARbyLeafSurf_pft(NZ)*ABS(BETA(N,M))
        RadPARDirStalkSurf_zsec(N,M,NZ) = RadPARbyStalkSurf_pft(NZ)*ABS(BETA(N,M))

        DO L=1,NumCanopyLayers1
          RadDifPAR_zsec(N,M,L,NZ) = 0.0_r8
          RadPAR_zsec(N,M,L,NZ)    = RadPARDirLeafSurf_zsec(N,M,NZ)
        enddo
      enddo
    enddo
  ENDDO D1100

  XAREA                                     = 1.00_r8/AREA3(NU)
  YAREA                                     = 1.00_r8/(AREA3(NU)*REAL(NumLeafZenithSectors1,R8))
  RadSWDiffusL                              = RadSWDiffus_col
  RadPARDiffusL                             = RadPARDiffus_col
  TAU_DirectRTransmit(NumCanopyLayers1+1) = 1.0_r8
  TAU_DifuseRTransmit(NumCanopyLayers1+1) = 1.0_r8
  RadSWFwdScat2NextL(NumCanopyLayers1+1)  = 0.0_r8
  RadPARFwdScat2NextL(NumCanopyLayers1+1) = 0.0_r8
  FracDirRadAbsorbtCum                      = 0.0_r8
  FracDifRadAbsorbtCum                      = 0.0_r8
  !distribute radiation from top of canopy to ground
  !it is only applicable to upland vegetation
  !
  !     RESET ARRAYS OF SUNLIT AND SHADED LEAF AREAS IN DIFFERENT
  !     LAYERS AND ANGLE CLASSES
  !
  !     TSURF,StemAreaZsec_lpft,SURF,StemAreaZsec_brch=leaf,stalk total,PFT surface area
  !
  !
  !     CALCULATE ABSORPTION, REFLECTION AND TRANSMISSION OF DIRECT AND
  !     DIFFUSE DOWNWARD TOTAL AND VISIBLE RADIATION BY EACH SPECIES
  !     NZ IN EACH LAYER L
  !
!  write(3355,*)I*100+J,'leaf',((LeafAreaZsec_lpft(N,L,1),N=1,NumLeafZenithSectors1),L=NumCanopyLayers1,1,-1)
!  write(3355,*)I*100+J,'stalk',((StemAreaZsec_lpft(N,L,1),N=1,NumLeafZenithSectors1),L=NumCanopyLayers1,1,-1)
  D1800: DO L=NumCanopyLayers1,1,-1
    !next layer is above snow, and above water
    IF(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO .AND. CanopyHeightZ_col(L-1).GE.DepthSurfWatIce-ZERO)THEN
      RadSWDiffusL           = RadSWDiffusL *TAU_DifuseRTransmit(L+1)+RadSWFwdScat2NextL(L+1)
      RadPARDiffusL          = RadPARDiffusL*TAU_DifuseRTransmit(L+1)+RadPARFwdScat2NextL(L+1)
      RadSWFwdScat2NextL(L)  = 0.0_r8
      RadPARFwdScat2NextL(L) = 0.0_r8
      RadSWBakScat2NextL(L)  = 0.0_r8
      RadPARBakScat2NextL(L) = 0.0_r8
      FracDirRadAbsorbt      = 0.0_r8
      FracDifRadAbsorbt      = 0.0_r8
      !
      !  RESET ACCUMULATORS OB ABSORBED, REFLECTED AND TRANSMITTED RADIATION
      !
      !
      D1500: DO NZ=1,NP
        RadDirSWbyLeaf_pft(NZ)   = 0.0_r8
        RadDirSWbyStalk_pft(NZ)  = 0.0_r8
        RadDirPARbyLeaf_pft(NZ)  = 0.0_r8
        RadDirPARbyStalk_pft(NZ) = 0.0_r8

        RadDifSWbyLeaf_pft(NZ)   = 0.0_r8
        RadDifSWbyStalk_pft(NZ)  = 0.0_r8
        RadDifPARbyLeaf_pft(NZ)  = 0.0_r8
        RadDifPARbyStalk_pft(NZ) = 0.0_r8

        bakScatRadDirSWbyLeaf_pft(NZ)   = 0.0_r8
        bakScatRadDirSWbyStalk_pft(NZ)  = 0.0_r8
        bakScatRadDirPARbyLeaf_pft(NZ)  = 0.0_r8
        bakScatRadDirPARbyStalk_pft(NZ) = 0.0_r8

        bakScatRadDifSWbyLeaf_pft(NZ)   = 0.0_r8
        bakScatRadDifSWbyStalk_pft(NZ)  = 0.0_r8
        bakScatRadDifPARbyLeaf_pft(NZ)  = 0.0_r8
        bakScatRadDifPARbyStalk_pft(NZ) = 0.0_r8

        fwdScatRadDirSWbyLeaf_pft(NZ)   = 0.0_r8
        fwdScatRadDirSWbyStalk_pft(NZ)  = 0.0_r8
        fwdScatRadDirPARbyLeaf_pft(NZ)  = 0.0_r8
        fwdScatRadDirPARbyStalk_pft(NZ) = 0.0_r8

        fwdScatRadDifSWbyLeaf_pft(NZ)   = 0.0_r8
        fwdScatRadDifSWbyStalk_pft(NZ)  = 0.0_r8
        fwdScatRadDifPARbyLeaf_pft(NZ)  = 0.0_r8
        fwdScatRadDifPARbyStalk_pft(NZ) = 0.0_r8

        !     LEAF SURFACE AREA IN EACH INCLINATION CLASS N, AZIMUTH CLASS M,
        !     LAYER L AND SPECIES NZ
        !
        !     UnselfShadeLeafArea=unself-shaded leaf area
        !     UnselfShadeLeafAreaAzclass=unself-shaded leaf area m-2 in each azimuth class
        !     TSurfLeaf=UnselfShadeLeafArea with shading from canopy layers above
        !     LeafIntceptArea=TSurfLeaf m-2
        !     UnselfShadeStalkArea=unself-shaded stalk area
        !     UnselfShadeStalkAreaAzclass=unself-shaded stalk area m-2 in each azimuth class
        !     TSurfStalk=UnselfShadeStalkArea with shading from canopy layers above
        !     StalkIntceptArea=TSurfStalk m-2
        !
        D1600: DO N=1,NumLeafZenithSectors1          
          UnselfShadeLeafArea        = LeafAreaZsec_lpft(N,L,NZ)*ClumpFactorNow_pft(NZ)
          UnselfShadeLeafAreaAzclass = UnselfShadeLeafArea*YAREA
          TSurfLeaf                  = UnselfShadeLeafArea*TAU_DirectRTransmit(L+1)
          LeafIntceptArea            = TSurfLeaf*XAREA

          UnselfShadeStalkArea        = StemAreaZsec_lpft(N,L,NZ)*StalkClumpFactor
          UnselfShadeStalkAreaAzclass = UnselfShadeStalkArea*YAREA
          TSurfStalk                  = UnselfShadeStalkArea*TAU_DirectRTransmit(L+1)
          StalkIntceptArea            = TSurfStalk*XAREA

          !
          !     ABSORPTION OF DIRECT RADIATION BY SUNLIT LEAF SURFACES
          !
          !     STOPZ=accumulated horizontal area of intercepted direct radiation
          !
          D1700: DO M=1,NumOfSkyAzimuthSects1
            RadDirSWbyLeaf_pft(NZ)   = RadDirSWbyLeaf_pft(NZ)+TSurfLeaf*RadSWbyLeafSurf_zsec(N,M,NZ)
            RadDirSWbyStalk_pft(NZ)  = RadDirSWbyStalk_pft(NZ)+TSurfStalk*RadSWbyStalkSurf_zsec(N,M,NZ)
            RadDirPARbyLeaf_pft(NZ)  = RadDirPARbyLeaf_pft(NZ)+TSurfLeaf*RadPARDirLeafSurf_zsec(N,M,NZ)
            RadDirPARbyStalk_pft(NZ) = RadDirPARbyStalk_pft(NZ)+TSurfStalk*RadPARDirStalkSurf_zsec(N,M,NZ)
            FracDirRadAbsorbt        = FracDirRadAbsorbt+(LeafIntceptArea+StalkIntceptArea)*BETX(N,M)

            !     BACKSCATTERING OF REFLECTED DIRECT RADIATION
      
            IF(iScatteringDirect(N,M).EQ.ibackward)THEN
              bakScatRadDirSWbyLeaf_pft(NZ)   = bakScatRadDirSWbyLeaf_pft(NZ)+TSurfLeaf*RadSWbyLeafSurf_zsec(N,M,NZ)
              bakScatRadDirSWbyStalk_pft(NZ)  = bakScatRadDirSWbyStalk_pft(NZ)+TSurfStalk*RadSWbyStalkSurf_zsec(N,M,NZ)
              bakScatRadDirPARbyLeaf_pft(NZ)  = bakScatRadDirPARbyLeaf_pft(NZ)+TSurfLeaf*RadPARDirLeafSurf_zsec(N,M,NZ)
              bakScatRadDirPARbyStalk_pft(NZ) = bakScatRadDirPARbyStalk_pft(NZ)+TSurfStalk*RadPARDirStalkSurf_zsec(N,M,NZ)
              !
              ! FORWARD SCATTERING OF REFLECTED DIRECT RADIATION
              !
            ELSE
              fwdScatRadDirSWbyLeaf_pft(NZ)   = fwdScatRadDirSWbyLeaf_pft(NZ)+TSurfLeaf*RadSWbyLeafSurf_zsec(N,M,NZ)
              fwdScatRadDirSWbyStalk_pft(NZ)  = fwdScatRadDirSWbyStalk_pft(NZ)+TSurfStalk*RadSWbyStalkSurf_zsec(N,M,NZ)
              fwdScatRadDirPARbyLeaf_pft(NZ)  = fwdScatRadDirPARbyLeaf_pft(NZ)+TSurfLeaf*RadPARDirLeafSurf_zsec(N,M,NZ)
              fwdScatRadDirPARbyStalk_pft(NZ) = fwdScatRadDirPARbyStalk_pft(NZ)+TSurfStalk*RadPARDirStalkSurf_zsec(N,M,NZ)
            ENDIF

            !     INTENSITY OF ABSORBED DIFFUSE RADIATION AT LEAF SURFACES
            !
            !     diffusSWLeafAbsorptAzclass,diffusSWStalkAbsorptAzclass,diffusPARLeafAbsorptAzclass,diffusPARStalkAbsorptAzclass=diffuse SW,PAR flux absorbed by leaf,stalk surf
            !     OMEGA,OMEGX=incident angle of diffuse radn at leaf,horizontal surface

            D1750: DO NN=1,NumOfLeafAzimuthSectors1
              diffusSWLeafAbsorptAzclass   = RadSWDiffusL*OMEGA(M,N,NN)*LeafSWabsorpty_pft(NZ)
              diffusSWStalkAbsorptAzclass  = RadSWDiffusL*OMEGA(M,N,NN)*StalkSWAbsorpty
              diffusPARLeafAbsorptAzclass  = RadPARDiffusL*OMEGA(M,N,NN)*LeafPARabsorpty_pft(NZ)
              diffusPARStalkAbsorptAzclass = RadPARDiffusL*OMEGA(M,N,NN)*StalkPARAbsorpty
              RadDifPAR_zsec(N,M,L,NZ)     = RadDifPAR_zsec(N,M,L,NZ)+diffusPARLeafAbsorptAzclass
              RadPAR_zsec(N,M,L,NZ)        = RadPAR_zsec(N,M,L,NZ)+diffusPARLeafAbsorptAzclass

              !     ABSORPTION OF DIFFUSE RADIATION BY SHADED LEAF SURFACES
              !
              !     FracDifRadAbsorbt=accumulated horizontal area of intercepted diffuse radiation
              !
              RadDifSWbyLeaf_pft(NZ)   = RadDifSWbyLeaf_pft(NZ)+UnselfShadeLeafArea*diffusSWLeafAbsorptAzclass
              RadDifSWbyStalk_pft(NZ)  = RadDifSWbyStalk_pft(NZ)+UnselfShadeStalkArea*diffusSWStalkAbsorptAzclass
              RadDifPARbyLeaf_pft(NZ)  = RadDifPARbyLeaf_pft(NZ)+UnselfShadeLeafArea*diffusPARLeafAbsorptAzclass
              RadDifPARbyStalk_pft(NZ) = RadDifPARbyStalk_pft(NZ)+UnselfShadeStalkArea*diffusPARStalkAbsorptAzclass
              FracDifRadAbsorbt        = FracDifRadAbsorbt+(UnselfShadeLeafAreaAzclass+UnselfShadeStalkAreaAzclass)*OMEGX(M,N,NN)

              !     BACKSCATTERING OF REFLECTED DIFFUSE RADIATION
              !
              IF(iScatteringDiffus(M,N,NN).EQ.ibackward)THEN
                bakScatRadDifSWbyLeaf_pft(NZ)   = bakScatRadDifSWbyLeaf_pft(NZ)+UnselfShadeLeafArea*diffusSWLeafAbsorptAzclass
                bakScatRadDifSWbyStalk_pft(NZ)  = bakScatRadDifSWbyStalk_pft(NZ)+UnselfShadeStalkArea*diffusSWStalkAbsorptAzclass
                bakScatRadDifPARbyLeaf_pft(NZ)  = bakScatRadDifPARbyLeaf_pft(NZ)+UnselfShadeLeafArea*diffusPARLeafAbsorptAzclass
                bakScatRadDifPARbyStalk_pft(NZ) = bakScatRadDifPARbyStalk_pft(NZ)+UnselfShadeStalkArea*diffusPARStalkAbsorptAzclass
                !
                !     FORWARD SCATTERING OF REFLECTED DIFFUSE RADIATION
                !
              ELSE
                fwdScatRadDifSWbyLeaf_pft(NZ)   = fwdScatRadDifSWbyLeaf_pft(NZ)+UnselfShadeLeafArea*diffusSWLeafAbsorptAzclass
                fwdScatRadDifSWbyStalk_pft(NZ)  = fwdScatRadDifSWbyStalk_pft(NZ)+UnselfShadeStalkArea*diffusSWStalkAbsorptAzclass
                fwdScatRadDifPARbyLeaf_pft(NZ)  = fwdScatRadDifPARbyLeaf_pft(NZ)+UnselfShadeLeafArea*diffusPARLeafAbsorptAzclass
                fwdScatRadDifPARbyStalk_pft(NZ) = fwdScatRadDifPARbyStalk_pft(NZ)+UnselfShadeStalkArea*diffusPARStalkAbsorptAzclass
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
      !     TAU_DirectRTransmit=transmission of direct radiation to next lower layer
      !
      IF(FracDirRadAbsorbtCum+FracDirRadAbsorbt.GT.1.0_r8)THEN
        IF(FracDirRadAbsorbt.GT.ZERO)THEN
          XTAUS=(1.0_r8-FracDirRadAbsorbtCum)/((1.0_r8-FracDirRadAbsorbtCum)-&
            (1.0_r8-FracDirRadAbsorbtCum-FracDirRadAbsorbt))
        ELSE
          XTAUS=0.0_r8
        ENDIF
        TAU_DirectRTransmit(L+1) = TAU_DirectRTransmit(L+1)*XTAUS
        FracDirRadAbsorbt        = FracDirRadAbsorbt*XTAUS

        D1510: DO NZ=1,NP
          RadDirSWbyLeaf_pft(NZ)   = RadDirSWbyLeaf_pft(NZ)*XTAUS
          RadDirSWbyStalk_pft(NZ)  = RadDirSWbyStalk_pft(NZ)*XTAUS
          RadDirPARbyLeaf_pft(NZ)  = RadDirPARbyLeaf_pft(NZ)*XTAUS
          RadDirPARbyStalk_pft(NZ) = RadDirPARbyStalk_pft(NZ)*XTAUS

          bakScatRadDirSWbyLeaf_pft(NZ)   = bakScatRadDirSWbyLeaf_pft(NZ)*XTAUS
          bakScatRadDirSWbyStalk_pft(NZ)  = bakScatRadDirSWbyStalk_pft(NZ)*XTAUS
          bakScatRadDirPARbyLeaf_pft(NZ)  = bakScatRadDirPARbyLeaf_pft(NZ)*XTAUS
          bakScatRadDirPARbyStalk_pft(NZ) = bakScatRadDirPARbyStalk_pft(NZ)*XTAUS

          fwdScatRadDirSWbyLeaf_pft(NZ)   = fwdScatRadDirSWbyLeaf_pft(NZ)*XTAUS
          fwdScatRadDirSWbyStalk_pft(NZ)  = fwdScatRadDirSWbyStalk_pft(NZ)*XTAUS
          fwdScatRadDirPARbyLeaf_pft(NZ)  = fwdScatRadDirPARbyLeaf_pft(NZ)*XTAUS
          fwdScatRadDirPARbyStalk_pft(NZ) = fwdScatRadDirPARbyStalk_pft(NZ)*XTAUS
        ENDDO D1510
      ENDIF
      !
      !     XTAUY=interception of diffuse radiation in current layer
      !     FracDifRadAbsorbt=accumulated interception of diffuse radiation from topmost layer
      !     TAU_DifuseRTransmit=transmission of diffuse radiation to next lower layer
      !
      IF(FracDifRadAbsorbtCum+FracDifRadAbsorbt.GT.1.0_r8)THEN
        XTAUY=(1.0_r8-FracDifRadAbsorbtCum)/((1.0_r8-FracDifRadAbsorbtCum)-&
          (1.0_r8-FracDifRadAbsorbtCum-FracDifRadAbsorbt))
        TAU_DifuseRTransmit(L+1) = TAU_DifuseRTransmit(L+1)*XTAUY
        FracDifRadAbsorbt        = FracDifRadAbsorbt*XTAUY
        D1520: DO NZ=1,NP
          RadDifSWbyLeaf_pft(NZ)   = RadDifSWbyLeaf_pft(NZ)*XTAUY
          RadDifSWbyStalk_pft(NZ)  = RadDifSWbyStalk_pft(NZ)*XTAUY
          RadDifPARbyLeaf_pft(NZ)  = RadDifPARbyLeaf_pft(NZ)*XTAUY
          RadDifPARbyStalk_pft(NZ) = RadDifPARbyStalk_pft(NZ)*XTAUY

          bakScatRadDifSWbyLeaf_pft(NZ)   = bakScatRadDifSWbyLeaf_pft(NZ)*XTAUY
          bakScatRadDifSWbyStalk_pft(NZ)  = bakScatRadDifSWbyStalk_pft(NZ)*XTAUY
          bakScatRadDifPARbyLeaf_pft(NZ)  = bakScatRadDifPARbyLeaf_pft(NZ)*XTAUY
          bakScatRadDifPARbyStalk_pft(NZ) = bakScatRadDifPARbyStalk_pft(NZ)*XTAUY

          fwdScatRadDifSWbyLeaf_pft(NZ)   = fwdScatRadDifSWbyLeaf_pft(NZ)*XTAUY
          fwdScatRadDifSWbyStalk_pft(NZ)  = fwdScatRadDifSWbyStalk_pft(NZ)*XTAUY
          fwdScatRadDifPARbyLeaf_pft(NZ)  = fwdScatRadDifPARbyLeaf_pft(NZ)*XTAUY
          fwdScatRadDifPARbyStalk_pft(NZ) = fwdScatRadDifPARbyStalk_pft(NZ)*XTAUY

          D1730: DO N=1,NumLeafZenithSectors1
            DO  M=1,NumOfSkyAzimuthSects1
              RadDifPAR_zsec(N,M,L,NZ) = RadDifPAR_zsec(N,M,L,NZ)*XTAUY
              RadPAR_zsec(N,M,L,NZ)    = RadPARDirLeafSurf_zsec(N,M,NZ)+RadDifPAR_zsec(N,M,L,NZ)
            enddo
          ENDDO D1730
        ENDDO D1520
      ENDIF
      !
      !     TOTAL RADIATION ABSORBED, REFLECTED AND TRANSMITTED BY ALL PFTs
      !
      !     RadSWbyCanopy_pft,TRADC,RADP,TRADP=total atmospheric SW,PAR absbd by each,all PFT
      !     FracDirRadAbsorbtCum,FracDifRadAbsorbtCum=accumulated interception of direct,diffuse radiation
      !     TAU_DirectRTransmit,TAU_DifuseRTransmit=transmission of direct,diffuse radiation to next lower layer
      !
      D1530: DO NZ=1,NP
        RadSWbyLeafT   = RadDirSWbyLeaf_pft(NZ)+RadDifSWbyLeaf_pft(NZ)
        RadSWbyStalkT  = RadDirSWbyStalk_pft(NZ)+RadDifSWbyStalk_pft(NZ)
        RadPARbyLeafT  = RadDirPARbyLeaf_pft(NZ)+RadDifPARbyLeaf_pft(NZ)
        RadPARbyStalkT = RadDirPARbyStalk_pft(NZ)+RadDifPARbyStalk_pft(NZ)

        bakScatRadSWbyLeafT   = bakScatRadDirSWbyLeaf_pft(NZ)+bakScatRadDifSWbyLeaf_pft(NZ)
        bakScatRadSWbyStalkT  = bakScatRadDirSWbyStalk_pft(NZ)+bakScatRadDifSWbyStalk_pft(NZ)
        bakScatRadPARbyLeafT  = bakScatRadDirPARbyLeaf_pft(NZ)+bakScatRadDifPARbyLeaf_pft(NZ)
        bakScatRadPARbyStalkT = bakScatRadDirPARbyStalk_pft(NZ)+bakScatRadDifPARbyStalk_pft(NZ)

        fwdScatRadSWbyLeafT   = fwdScatRadDirSWbyLeaf_pft(NZ)+fwdScatRadDifSWbyLeaf_pft(NZ)
        fwdScatRadSWbyStalkT  = fwdScatRadDirSWbyStalk_pft(NZ)+fwdScatRadDifSWbyStalk_pft(NZ)
        fwdScatRadPARbyLeafT  = fwdScatRadDirPARbyLeaf_pft(NZ)+fwdScatRadDifPARbyLeaf_pft(NZ)
        fwdScatRadPARbyStalkT = fwdScatRadDirPARbyStalk_pft(NZ)+fwdScatRadDifPARbyStalk_pft(NZ)

        RadSWFwdScat2NextL(L)=RadSWFwdScat2NextL(L)+(RadSWbyLeafT*RadSWLeafTransmis_pft(NZ) &
          +fwdScatRadSWbyLeafT*RadSWLeafAlbedo_pft(NZ) &
          +fwdScatRadSWbyStalkT*RadSWStalkAlbedo)*YAREA
        RadPARFwdScat2NextL(L)=RadPARFwdScat2NextL(L)+(RadPARbyLeafT*RadPARLeafTransmis_pft(NZ) &
          +fwdScatRadPARbyLeafT*CanopyPARalbedo_pft(NZ) &
          +fwdScatRadPARbyStalkT*StalkAlbedo4PARRad)*YAREA
        RadSWBakScat2NextL(L)=RadSWBakScat2NextL(L)+(bakScatRadSWbyLeafT*RadSWLeafAlbedo_pft(NZ)+ &
          bakScatRadSWbyStalkT*RadSWStalkAlbedo)*YAREA
        RadPARBakScat2NextL(L)=RadPARBakScat2NextL(L)+(bakScatRadPARbyLeafT*CanopyPARalbedo_pft(NZ)+&
          bakScatRadPARbyStalkT*StalkAlbedo4PARRad)*YAREA
          
        !accumulate shortwave radiation on canopy 
        RadSWbyCanopy_pft(NZ)  = RadSWbyCanopy_pft(NZ)+RadSWbyLeafT+RadSWbyStalkT
        RadPARbyCanopy_pft(NZ) = RadPARbyCanopy_pft(NZ)+RadPARbyLeafT+RadPARbyStalkT
      ENDDO D1530
      FracDirRadAbsorbtCum   = FracDirRadAbsorbtCum+FracDirRadAbsorbt
      FracDifRadAbsorbtCum   = FracDifRadAbsorbtCum+FracDifRadAbsorbt
      TAU_DirectRTransmit(L) = 1.0_r8-FracDirRadAbsorbtCum
      TAU_RadThru(L)         = 1.0_r8-TAU_DirectRTransmit(L)
      TAU_DifuseRTransmit(L) = 1.0_r8-FracDifRadAbsorbtCum
    ELSE
      RadSWFwdScat2NextL(L)  = RadSWFwdScat2NextL(L+1)
      RadPARFwdScat2NextL(L) = RadPARFwdScat2NextL(L+1)
      TAU_DirectRTransmit(L) = TAU_DirectRTransmit(L+1)
      TAU_RadThru(L)         = 1.0_r8-TAU_DirectRTransmit(L)
      TAU_DifuseRTransmit(L) = TAU_DifuseRTransmit(L+1)
    ENDIF

  ENDDO D1800

  !
  !     DIRECT AND DIFFUSE RADIATION ABSORBED AT GROUND SURFACE
  !
  !     RADSG,RADYG,RAPSG,RAPYG=direct,diffuse SW,PAR at horiCanopyHeightZ_col ground surface
  !     RADS,RadPARDirect_col=solar beam direct SW,PAR flux
  !     TAU_DirectRTransmit,TAU_DifuseRTransmit=transmission of direct,diffuse radiation below canopy
  !     RadSWDiffusL,RadPARDiffusL=solar beam diffuse SW,PAR flux
  !     RadSW_Grnd,RadPAR_Grnd=total SW,PAR at ground surface
  !     GrndIncidSolarAngle,OMEGAG=incident solar,sky angle at ground surface
  !

  RADSG = RadSWDirect_col*TAU_DirectRTransmit(1)
  RADYG = RadSWDiffusL*TAU_DifuseRTransmit(1)+RadSWFwdScat2NextL(1)
  RAPSG = RadPARDirect_col*TAU_DirectRTransmit(1)
  RAPYG = RadPARDiffusL*TAU_DifuseRTransmit(1)+RadPARFwdScat2NextL(1)

  RadSW_Grnd  = ABS(GrndIncidSolarAngle)*RADSG
  RadPAR_Grnd = ABS(GrndIncidSolarAngle)*RAPSG
!  write(3333,*)I*100+J,RadSW_Grnd,RADSG,RadSWDirect_col,'0RadSWDirect_col',TAU_DirectRTransmit  

  D20: DO N=1,NumOfSkyAzimuthSects1
    RadSW_Grnd  = RadSW_Grnd+ABS(OMEGAG(N))*RADYG
    RadPAR_Grnd = RadPAR_Grnd+ABS(OMEGAG(N))*RAPYG
  ENDDO D20 
!  write(3333,*)I*100+J,RadSW_Grnd,RADYG,'1RadSWDirect_col',TAU_DifuseRTransmit(1),RadSWFwdScat2NextL(1)
!  if(I>=16)stop

  !
  !     RADIATION REFLECTED FROM GROUND SURFACE
  !
  !     VHCPW,VLHeatCapSnowMin_col=current,minimum snowpack heat capacity
  !     SnowpackAlbedo,VcumDrySnoWE_col,VcumWatSnow_col,VcumIceSnow_col=snowpack albedo,snow,water,ice volume
  !     GrndAlbedo,SoilAlbedo,FracGrndBySnow=ground,soil albedo,snow cover fraction
  !     THETW1=soil surface water content
  !     RadSWBakScat2NextL,RadDirPARbyLeaf_pft=SW,PAR backscatter from ground surface
  !     TRADG,TRAPG=SW,PAR absorbed by ground surface
  !
  IF(VLHeatCapSurfSnow_col.GT.VLHeatCapSnowMin_col)THEN
    SnowpackAlbedo=(0.85_r8*VcumDrySnoWE_col+0.30_r8*VcumIceSnow_col+0.06_r8*VcumWatSnow_col) &
      /(VcumDrySnoWE_col+VcumIceSnow_col+VcumWatSnow_col)
    !the following partition differs from that used in the surface physics module  
    FracGrndBySnow = AMIN1((SnowDepth/0.07_r8)**2._r8,1.0_r8)
    GrndAlbedo     = FracGrndBySnow*SnowpackAlbedo+(1.0_r8-FracGrndBySnow)*SoilAlbedo
  ELSE
    IF(VLSoilPoreMicP_vr(NU).GT.ZEROS2)THEN
      THETW1=AMIN1(POROS1,VLWatMicP_vr(NU)/VLSoilMicP_vr(NU))
    ELSE
      THETW1=0.0_r8
    ENDIF
    GrndAlbedo=AMIN1(SurfAlbedo_col,SoilAlbedo+AZMAX1(SurfAlbedo_col-THETW1))
  ENDIF
  RadSWBakScat2NextL(0)  = RadSW_Grnd*GrndAlbedo*0.25_r8
  RadPARBakScat2NextL(0) = RadPAR_Grnd*GrndAlbedo*0.25_r8
  !
  !     ADD RADIATION FROM SCATTERING THROUGH CANOPY LAYERS
  !
  !     RadSWBakScat2NextL,RadPARBakScat2NextL=total backscattered SW,PAR to next layer
  !     RadSWFwdScat2NextL,RadPARFwdScat2NextL=total fwd scattered SW,PAR to next layer
  !     diffusSWLeafAbsorptAzclass,diffusSWStalkAbsorptAzclass,diffusPARLeafAbsorptAzclass,diffusPARStalkAbsorptAzclass=leaf,stalk SW,PAR absbd fwd+back flux
  !     RadDifSWbyLeaf_pft,RadDifSWbyStalk_pft,RadDifPARbyLeaf_pft,RadDifPARbyStalk_pft=total leaf,stalk SW,PAR absbd fwd+back
  !     RadSWbyCanopy_pft,TRADC,RadPARbyCanopy_pft,TRADP=total SW,PAR absbd by each,all PFT
  !
  RadSWDiffusL           = 0.0_r8
  RadPARDiffusL          = 0.0_r8
  TAU_DifuseRTransmit(0) = 1.0_r8
  RadSWFwdScat2NextL(0)  = 0.0_r8
  RadPARFwdScat2NextL(0) = 0.0_r8

  D2800: DO L=1,NumCanopyLayers1
    !the following line shuts off radiation when it is below water 
    IF(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO .AND. CanopyHeightZ_col(L-1).GE.DepthSurfWatIce-ZERO)THEN
      RadSWDiffusL           = RadSWDiffusL*TAU_DifuseRTransmit(L-1)+RadSWFwdScat2NextL(L-1)+RadSWBakScat2NextL(L-1)
      RadPARDiffusL          = RadPARDiffusL*TAU_DifuseRTransmit(L-1)+RadPARFwdScat2NextL(L-1)+RadPARBakScat2NextL(L-1)
      RadSWFwdScat2NextL(L)  = 0.0
      RadPARFwdScat2NextL(L) = 0.0_r8
      D2500: DO NZ=1,NP
        RadDifSWbyLeaf_pft(NZ)   = 0.0_r8
        RadDifSWbyStalk_pft(NZ)  = 0.0_r8
        RadDifPARbyLeaf_pft(NZ)  = 0.0_r8
        RadDifPARbyStalk_pft(NZ) = 0.0_r8

        D2600: DO N=1,NumLeafZenithSectors1
          UnselfShadeLeafArea  = LeafAreaZsec_lpft(N,L,NZ)*ClumpFactorNow_pft(NZ)
          UnselfShadeStalkArea = StemAreaZsec_lpft(N,L,NZ)*StalkClumpFactor
          D2700: DO M=1,NumOfSkyAzimuthSects1
            D2750: DO NN=1,NumOfLeafAzimuthSectors1
              diffusSWLeafAbsorptAzclass   = RadSWDiffusL*OMEGA(M,N,NN)*LeafSWabsorpty_pft(NZ)
              diffusSWStalkAbsorptAzclass  = RadSWDiffusL*OMEGA(M,N,NN)*StalkSWAbsorpty
              diffusPARLeafAbsorptAzclass  = RadPARDiffusL*OMEGA(M,N,NN)*LeafPARabsorpty_pft(NZ)
              diffusPARStalkAbsorptAzclass = RadPARDiffusL*OMEGA(M,N,NN)*StalkPARAbsorpty

              RadDifPAR_zsec(N,M,L,NZ) = RadDifPAR_zsec(N,M,L,NZ)+diffusPARLeafAbsorptAzclass
              RadPAR_zsec(N,M,L,NZ)    = RadPAR_zsec(N,M,L,NZ)+diffusPARLeafAbsorptAzclass
              RadDifSWbyLeaf_pft(NZ)   = RadDifSWbyLeaf_pft(NZ)+UnselfShadeLeafArea*diffusSWLeafAbsorptAzclass
              RadDifSWbyStalk_pft(NZ)  = RadDifSWbyStalk_pft(NZ)+UnselfShadeStalkArea*diffusSWStalkAbsorptAzclass
              RadDifPARbyLeaf_pft(NZ)  = RadDifPARbyLeaf_pft(NZ)+UnselfShadeLeafArea*diffusPARLeafAbsorptAzclass
              RadDifPARbyStalk_pft(NZ) = RadDifPARbyStalk_pft(NZ)+UnselfShadeStalkArea*diffusPARStalkAbsorptAzclass
            ENDDO D2750
          ENDDO D2700
        ENDDO D2600
        RadSWFwdScat2NextL(L)  = RadSWFwdScat2NextL(L)+RadDifSWbyLeaf_pft(NZ)*RadSWLeafTransmis_pft(NZ)*YAREA
        RadPARFwdScat2NextL(L) = RadPARFwdScat2NextL(L)+RadDifPARbyLeaf_pft(NZ)*RadPARLeafTransmis_pft(NZ)*YAREA
        RadSWbyCanopy_pft(NZ)  = RadSWbyCanopy_pft(NZ)+RadDifSWbyLeaf_pft(NZ)+RadDifSWbyStalk_pft(NZ)
        RadPARbyCanopy_pft(NZ) = RadPARbyCanopy_pft(NZ)+RadDifPARbyLeaf_pft(NZ)+RadDifPARbyStalk_pft(NZ)
      ENDDO D2500
    ELSE
      RadSWFwdScat2NextL(L)  = RadSWFwdScat2NextL(L-1)
      RadPARFwdScat2NextL(L) = RadPARFwdScat2NextL(L-1)
      RadSWBakScat2NextL(L)  = RadSWBakScat2NextL(L-1)
      RadPARBakScat2NextL(L) = RadPARBakScat2NextL(L-1)
    ENDIF
  ENDDO D2800        
  if(I==6.and.J==22.and..false.)then
  write(456,*)'swfwd','pafwd','swbak','pabak'
  write(456,*)(RadSWFwdScat2NextL(L),L=0,NumCanopyLayers1+1)
  write(456,*)(RadPARFwdScat2NextL(L),L=0,NumCanopyLayers1+1)
  write(456,*)(RadSWBakScat2NextL(L),L=0,NumCanopyLayers1+1)
  write(456,*)(RadPARBakScat2NextL(L),L=0,NumCanopyLayers1+1)
  write(456,*)'parz'
  DO L=1,NumCanopyLayers1  
  write(456,*)((RadPAR_zsec(N,M,L,1),N=1,NumLeafZenithSectors1),M=1,NumOfSkyAzimuthSects1)
  write(456,*)((RadDifPAR_zsec(N,M,L,1),N=1,NumLeafZenithSectors1),M=1,NumOfSkyAzimuthSects1)
  ENDDO
  endif
  end associate
  end subroutine MultiCanLayerRadiation
  ![tail]
end module SurfaceRadiationMod
