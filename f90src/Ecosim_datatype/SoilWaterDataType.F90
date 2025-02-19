module SoilWaterDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8),target,allocatable ::  TXGridSurfRunoff_2DH(:,:)                 !
  real(r8),target,allocatable ::  THeatXGridBySurfRunoff_2DH(:,:)           !
  logical, target,allocatable ::  iPondFlag_col(:,:)                        !
  integer, target,allocatable ::  iPondBotLev_col(:,:)                      !Bottom level ID
  real(r8),target,allocatable ::  ThetaAir_vr(:,:,:)                        !air concentration [m3 m-3]
  real(r8),target,allocatable ::  VLsoiAirP_vr(:,:,:)                       !soil air content [m3 d-2]
  real(r8),target,allocatable ::  THETW_vr(:,:,:)                           !volumetric water content [m3 m-3]
  real(r8),target,allocatable ::  THETI_vr(:,:,:)                           !volumetric ice content [m3 m-3]
  real(r8),target,allocatable ::  ThetaH2OZ_vr(:,:,:)                       !volumetric moblize water [m3 m-3]
  real(r8),target,allocatable ::  ThetaICEZ_vr(:,:,:)                       !volumetric mobile ice [m3 m-3]
  real(r8),target,allocatable ::  VLWatMicP_vr(:,:,:)                       !soil micropore water content [m3 d-2]
  real(r8),target,allocatable ::  VLiceMicP_vr(:,:,:)                       !soil micropore ice content   [m3 d-2]
  real(r8),target,allocatable ::  VLWatMacP_vr(:,:,:)                       !soil macropore water content [m3 d-2]
  real(r8),target,allocatable ::  PSISoilMatricP_vr(:,:,:)                  !soil micropore matric water potential [MPa]
  real(r8),target,allocatable ::  ElvAdjstedSoilH2OPSIMPa_vr(:,:,:)         !elevation adjusted total soil micropore total water potential [MPa]
  real(r8),target,allocatable ::  VLWatMicPX_vr(:,:,:)                      !soil micropore water content before wetting front [m3 d-2]
  real(r8),target,allocatable ::  FWatExMacP2MicP_vr(:,:,:)                 !soil macropore - micropore water transfer [m3 d-2 h-1]
  real(r8),target,allocatable ::  VLiceMacP_vr(:,:,:)                       !soil macropore ice content [m3 d-2]
  real(r8),target,allocatable ::  VLWatMicPM_vr(:,:,:,:)                    !soil micropore water content, [m3 d-2]
  real(r8),target,allocatable ::  VLWatMacPM_vr(:,:,:,:)                    !soil macropore water content, [m3 d-2]
  real(r8),target,allocatable ::  VLsoiAirPM_vr(:,:,:,:)                    !soil air content, [m3 d-2]
  real(r8),target,allocatable ::  FILMM_vr(:,:,:,:)                         !soil water film thickness , [m]
  real(r8),target,allocatable ::  WaterTBLSlope_col(:,:)                    !slope of water table relative to surface slope, [-]
  real(r8),target,allocatable ::  WtblDepzTile_col(:,:)                     !depth of artificial water table
  real(r8),target,allocatable ::  TileWaterTable_col(:,:)                   !artificial water table depth, [m]
  real(r8),target,allocatable ::  DTBLD(:,:)                                !depth of artificial water table adjusted for elevation
  real(r8),target,allocatable ::  DepzIntWTBL_col(:,:)                      !internal water table depth, [m]
  real(r8),target,allocatable ::  ExtWaterTablet0_col(:,:)                  !initial external water table depth, elevation corrected [m]
  real(r8),target,allocatable ::  ExtWaterTable_col(:,:)                    !current external water table depth, elevation corrected [m]
  real(r8),target,allocatable ::  NatWtblDepz_col(:,:)                      !external water table depth, [m]
  real(r8),target,allocatable ::  EnergyImpact4ErosionM(:,:,:)              !total energy impact for erosion
  real(r8),target,allocatable ::  XVLMobileWaterLitRM(:,:,:)                !excess water+ice
  real(r8),target,allocatable ::  XVLMobileWatMicPM(:,:,:)                  !excess water
  real(r8),target,allocatable ::  XVLiceMicPM(:,:,:)                        !excess ice
  real(r8),target,allocatable ::  HydroCond_3D(:,:,:,:,:)                   !hydraulic conductivity at different moisture levels
  real(r8),target,allocatable ::  HydroCondMacP_vr(:,:,:)                   !macropore hydraulic conductivity, [m MPa-1 h-1]
  real(r8),target,allocatable ::  HydroCondMicP4RootUptake_vr(:,:,:)        !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]
  real(r8),target,allocatable ::  SurfRunoffWatFluxM_2DH(:,:,:)             !runoff water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  RunoffVelocityM_col(:,:,:)                !runoff velocity, [m t-1]
  integer,target,allocatable ::  IFLBM(:,:,:,:,:)                           !flag for directional surface runoff
  logical,target,allocatable ::  XGridRunoffFlag(:,:,:,:)                   !enables or disables boundary water flux depending on aspect, [-]
  integer,target,allocatable ::  IFLBH(:,:,:,:)                             !flag for directional runoff, related to IFLBM
  real(r8),target,allocatable ::  RechargNorthSubSurf(:,:)                  !northern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargEastSubSurf(:,:)                   !eastern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargSouthSubSurf(:,:)                  !southern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargWestSubSurf(:,:)                   !western subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargRateNorthWTBL(:,:)                 !northern subsurface boundary water flux rate constant, [h-1]
  real(r8),target,allocatable ::  RechargRateEastWTBL(:,:)                  !eastern subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RechargRateSouthWTBL(:,:)                 !southern subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RechargRateWestWTBL(:,:)                  !western subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RechargNorthSurf(:,:)                     !northern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargEastSurf(:,:)                      !eastern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargSouthSurf(:,:)                     !southern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargWestSurf(:,:)                      !western surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargBottom_col(:,:)                    !lower subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  WaterFlow2MicPM_3D(:,:,:,:,:)             !micropore water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  WaterFlow2MacPM_3D(:,:,:,:,:)             !macropore water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  ReductVLsoiAirPM_vr(:,:,:,:)              !change in soil air volume for layer from last to current iteration, [g d-2 t-1] >0, shrink
  real(r8),target,allocatable ::  FWatExMacP2MicPM_vr(:,:,:,:)              !soil macropore - micropore water transfer, [g d-2 t-1]
  real(r8),target,allocatable ::  WatFlowSno2MicPM_col(:,:,:)               !meltwater flux into soil micropores
  real(r8),target,allocatable ::  WatFlowSno2MacPM_col(:,:,:)               !meltwater flux into soil macropores
  real(r8),target,allocatable ::  AirFilledSoilPoreM_vr(:,:,:,:)            !air-filled soil porosity, [m3 m-3]
  real(r8),target,allocatable ::  TortMicPM_vr(:,:,:,:)                     !soil tortuosity, []
  real(r8),target,allocatable ::  TortMacPM_vr(:,:,:,:)                     !macropore tortuosity, []
  real(r8),target,allocatable ::  DiffusivitySolutEffM_vr(:,:,:,:)          !coefficient for dissolution - volatilization, []
  real(r8),target,allocatable ::  SoilResit4RootPentrate_vr(:,:,:)          !soil hydraulic resistance, [MPa h m-2]
  real(r8),target,allocatable ::  PSISE_vr(:,:,:)                           !soil water potential at saturation, [Mpa]
  real(r8),target,allocatable ::  PSISoilAirEntry(:,:,:)                    !soil water potential at air entry, [Mpa]
  real(r8),target,allocatable ::  PSISoilOsmotic_vr(:,:,:)                  !osmotic soil water potential , [Mpa]
  real(r8),target,allocatable ::  PSIGrav_vr(:,:,:)                         !gravimetric soil water potential , [Mpa]
  real(r8),target,allocatable ::  SoilWatAirDry_vr(:,:,:)                   !air-dry water content, [m3 m-3]
  real(r8),target,allocatable ::  ThetaSat_vr(:,:,:)                        !micropore class water content
  real(r8),target,allocatable ::  WaterFlowSoiMicPX_3D(:,:,:,:)             !unsaturated water flux , [m3 d-2 h-1]
  real(r8),target,allocatable ::  EvapoTransp_col(:,:)                      !evapotranspiration
  real(r8),target,allocatable ::  QEvap_CumYr_col(:,:)                      !total evaporation, [m3 d-2]
  real(r8),target,allocatable ::  QRain_CumYr_col(:,:)                      !total precipitation, [m3 d-2]
  real(r8),target,allocatable ::  Qrunoff_CumYr_col(:,:)                    !total surface runoff, [m3 d-2]
  real(r8),target,allocatable ::  WatMass_col(:,:)                          !total soil water content, [m3 d-2]
  real(r8),target,allocatable ::  H2OLoss_CumYr_col(:,:)                    !total subsurface water flux, [m3 d-2]
  real(r8),target,allocatable ::  QDrain_col(:,:)                           !total water drainage below root zone, [m3 d-2]
  real(r8),target,allocatable ::  XGridSurfRunoff_2DH(:,:,:,:)              !soil surface runoff water, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HeatXGridBySurfRunoff_2DH(:,:,:,:)        !soil surface runoff heat, [MJ d-2 h-1]
  real(r8),target,allocatable ::  QRunSurf_col(:,:)                         !runoff from surface water, [m3 d-2 h-1]
  real(r8),target,allocatable ::  QDischar_col(:,:)                         !water discharge, [m3 d-2 h-1]
  real(r8),target,allocatable ::  QflxSurfRunoffM_2DH(:,:,:,:,:)            !surface runoff in iteration M,
  real(r8),target,allocatable ::  Qinflx2Soil_col(:,:)                      !infiltration into soil [m3 d-2 h-1]
  real(r8),target,allocatable ::  SoilWatMassBeg_col(:,:)                   !soil water mass at the begnining of time step
  real(r8),target,allocatable ::  SoilWatMassEnd_col(:,:)                   !soil water mass at the end of time step
  real(r8),target,allocatable ::  Rain2Soil_col(:,:)                        !water flow into soil due to precipitation (+ surface irrigation), [m3 H2O/d2/h]
  real(r8),target,allocatable ::  QdewCanopy_CumYr_pft(:,:,:)
  real(r8),target,allocatable :: QSnoWatXfer2Soil_col(:,:)
  real(r8),target,allocatable :: QSnoIceXfer2Soil_col(:,:)
  real(r8),target,allocatable :: PrecipAtm2LandSurf_col(:,:)                !precipiation from atmosphere to land surface 
  real(r8),target,allocatable :: RainPrecThrufall_col(:,:)                  !precipitation through canopy [m3 H2O d-2 h-1]
  real(r8),target,allocatable :: RainPrec2Sno_col(:,:)                      !rainfall to snow [m3 H2O d-2 h-1]
  real(r8),target,allocatable :: Rain2ExposedSurf_col(:,:)                  !rainfall to exposed surface [m3 H2O d-2 h-1]
  real(r8),target,allocatable :: QWatIntLaterFlow_col(:,:)                  !Internal lateral flow between grids [m3 H2O d-2 h-1]
  private :: InitAllocate
  contains

  subroutine InitSoilWater
  implicit none

  call InitAllocate
  end subroutine InitSoilWater

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none

  allocate(iPondBotLev_col(JY,JX)); iPondBotLev_col=0
  allocate(iPondFlag_col(JY,JX)); iPondFlag_col =.false.
  allocate(QWatIntLaterFlow_col(JY,JX)); QWatIntLaterFlow_col=0._r8
  allocate(Rain2ExposedSurf_col(JY,JX)); Rain2ExposedSurf_col=0._r8
  allocate(RainPrec2Sno_col(JY,JX)); RainPrec2Sno_col = 0._r8
  allocate(RainPrecThrufall_col(JY,JX)); RainPrecThrufall_col=0._r8
  allocate(PrecipAtm2LandSurf_col(JY,JX)); PrecipAtm2LandSurf_col=0._r8
  allocate(QSnoIceXfer2Soil_col(JY,JX)); QSnoIceXfer2Soil_col=0._r8
  allocate(QSnoWatXfer2Soil_col(JY,JX)); QSnoWatXfer2Soil_col=0._r8
  allocate(Rain2Soil_col(JY,JX));  Rain2Soil_col=0._r8
  allocate(SoilWatMassBeg_col(JY,JX)); SoilWatMassBeg_col=0._r8
  allocate(SoilWatMassEnd_col(JY,JX)); SoilWatMassEnd_col=0._r8
  allocate(QdewCanopy_CumYr_pft(JZ,JY,JX)); QdewCanopy_CumYr_pft=0._r8
  allocate(EvapoTransp_col(JY,JX)); EvapoTransp_col=0._r8
  allocate(Qinflx2Soil_col(JY,JX)); Qinflx2Soil_col=0._r8
  allocate(ThetaAir_vr(0:JZ,JY,JX));  ThetaAir_vr=0._r8
  allocate(VLsoiAirP_vr(0:JZ,JY,JX));   VLsoiAirP_vr=0._r8
  allocate(THETW_vr(0:JZ,JY,JX));  THETW_vr=0._r8
  allocate(THETI_vr(0:JZ,JY,JX));  THETI_vr=0._r8
  allocate(ThetaH2OZ_vr(0:JZ,JY,JX)); ThetaH2OZ_vr=0._r8
  allocate(ThetaICEZ_vr(0:JZ,JY,JX)); ThetaICEZ_vr=0._r8
  allocate(VLWatMicP_vr(0:JZ,JY,JX));   VLWatMicP_vr=0._r8
  allocate(VLiceMicP_vr(0:JZ,JY,JX));   VLiceMicP_vr=0._r8
  allocate(VLWatMacP_vr(JZ,JY,JX));    VLWatMacP_vr=0._r8
  allocate(PSISoilMatricP_vr(0:JZ,JY,JX));  PSISoilMatricP_vr=0._r8
  allocate(ElvAdjstedSoilH2OPSIMPa_vr(0:JZ,JY,JX));  ElvAdjstedSoilH2OPSIMPa_vr=0._r8
  allocate(VLWatMicPX_vr(0:JZ,JY,JX));  VLWatMicPX_vr=0._r8
  allocate(FWatExMacP2MicP_vr(JZ,JY,JX));     FWatExMacP2MicP_vr=0._r8
  allocate(VLiceMacP_vr(JZ,JY,JX));    VLiceMacP_vr=0._r8
  allocate(VLWatMicPM_vr(60,0:JZ,JY,JX));VLWatMicPM_vr=0._r8
  allocate(VLWatMacPM_vr(60,JZ,JY,JX));VLWatMacPM_vr=0._r8
  allocate(VLsoiAirPM_vr(60,0:JZ,JY,JX));VLsoiAirPM_vr=0._r8
  allocate(FILMM_vr(60,0:JZ,JY,JX));FILMM_vr=0._r8
  allocate(WaterTBLSlope_col(JY,JX));       WaterTBLSlope_col=0._r8
  allocate(WtblDepzTile_col(JY,JX));      WtblDepzTile_col=0._r8
  allocate(TileWaterTable_col(JY,JX));       TileWaterTable_col=0._r8
  allocate(DTBLD(JY,JX));       DTBLD=0._r8
  allocate(DepzIntWTBL_col(JY,JX));       DepzIntWTBL_col=0._r8
  allocate(ExtWaterTablet0_col(JY,JX));       ExtWaterTablet0_col=0._r8
  allocate(ExtWaterTable_col(JY,JX));       ExtWaterTable_col=0._r8
  allocate(NatWtblDepz_col(JY,JX));       NatWtblDepz_col=0._r8
  allocate(EnergyImpact4ErosionM(60,JY,JX));   EnergyImpact4ErosionM=0._r8
  allocate(XVLMobileWaterLitRM(60,JY,JX));   XVLMobileWaterLitRM=0._r8
  allocate(XVLMobileWatMicPM(60,JY,JX));   XVLMobileWatMicPM=0._r8
  allocate(XVLiceMicPM(60,JY,JX));   XVLiceMicPM=0._r8
  allocate(HydroCond_3D(3,100,0:JZ,JY,JX));HydroCond_3D=0._r8
  allocate(HydroCondMacP_vr(JZ,JY,JX));     HydroCondMacP_vr=0._r8
  allocate(HydroCondMicP4RootUptake_vr(JZ,JY,JX));     HydroCondMicP4RootUptake_vr=0._r8
  allocate(SurfRunoffWatFluxM_2DH(60,JV,JH));      SurfRunoffWatFluxM_2DH=0._r8
  allocate(RunoffVelocityM_col(60,JY,JX));      RunoffVelocityM_col=0._r8
  allocate(IFLBM(60,2,2,JY,JX));IFLBM=0
  allocate(XGridRunoffFlag(2,2,JY,JX));   XGridRunoffFlag=.false.
  allocate(IFLBH(2,2,JY,JX));   IFLBH=0
  allocate(RechargNorthSubSurf(JY,JX));      RechargNorthSubSurf=0._r8
  allocate(RechargEastSubSurf(JY,JX));      RechargEastSubSurf=0._r8
  allocate(RechargSouthSubSurf(JY,JX));      RechargSouthSubSurf=0._r8
  allocate(RechargWestSubSurf(JY,JX));      RechargWestSubSurf=0._r8
  allocate(RechargRateNorthWTBL(JY,JX));      RechargRateNorthWTBL=0._r8
  allocate(RechargRateEastWTBL(JY,JX));      RechargRateEastWTBL=0._r8
  allocate(RechargRateSouthWTBL(JY,JX));      RechargRateSouthWTBL=0._r8
  allocate(RechargRateWestWTBL(JY,JX));      RechargRateWestWTBL=0._r8
  allocate(RechargNorthSurf(JY,JX));       RechargNorthSurf=0._r8
  allocate(RechargEastSurf(JY,JX));       RechargEastSurf=0._r8
  allocate(RechargSouthSurf(JY,JX));       RechargSouthSurf=0._r8
  allocate(RechargWestSurf(JY,JX));       RechargWestSurf=0._r8
  allocate(RechargBottom_col(JY,JX));       RechargBottom_col=0._r8
  allocate(WaterFlow2MicPM_3D(60,3,JD,JV,JH));WaterFlow2MicPM_3D=0._r8
  allocate(WaterFlow2MacPM_3D(60,3,JD,JV,JH));WaterFlow2MacPM_3D=0._r8
  allocate(ReductVLsoiAirPM_vr(60,JZ,JY,JX));  ReductVLsoiAirPM_vr=0._r8
  allocate(FWatExMacP2MicPM_vr(60,JZ,JY,JX)); FWatExMacP2MicPM_vr=0._r8
  allocate(WatFlowSno2MicPM_col(60,JY,JX));    WatFlowSno2MicPM_col=0._r8
  allocate(WatFlowSno2MacPM_col(60,JY,JX));    WatFlowSno2MacPM_col=0._r8
  allocate(AirFilledSoilPoreM_vr(60,0:JZ,JY,JX));AirFilledSoilPoreM_vr=0._r8
  allocate(TortMicPM_vr(60,0:JZ,JY,JX));TortMicPM_vr=0._r8
  allocate(TortMacPM_vr(60,JZ,JY,JX)); TortMacPM_vr=0._r8
  allocate(DiffusivitySolutEffM_vr(60,0:JZ,JY,JX));DiffusivitySolutEffM_vr=0._r8
  allocate(SoilResit4RootPentrate_vr(JZ,JY,JX));     SoilResit4RootPentrate_vr=0._r8
  allocate(PSISE_vr(0:JZ,JY,JX));  PSISE_vr=0._r8
  allocate(PSISoilAirEntry(0:JZ,JY,JX));  PSISoilAirEntry=0._r8
  allocate(PSISoilOsmotic_vr(0:JZ,JY,JX));  PSISoilOsmotic_vr=0._r8
  allocate(PSIGrav_vr(0:JZ,JY,JX));  PSIGrav_vr=0._r8
  allocate(SoilWatAirDry_vr(0:JZ,JY,JX));  SoilWatAirDry_vr=0._r8
  allocate(ThetaSat_vr(0:JZ,JY,JX));  ThetaSat_vr=0._r8
  allocate(WaterFlowSoiMicPX_3D(3,JD,JV,JH));   WaterFlowSoiMicPX_3D=0._r8
  allocate(QEvap_CumYr_col(JY,JX));       QEvap_CumYr_col=0._r8
  allocate(QRain_CumYr_col(JY,JX));       QRain_CumYr_col=0._r8
  allocate(Qrunoff_CumYr_col(JY,JX));        Qrunoff_CumYr_col=0._r8
  allocate(WatMass_col(JY,JX));       WatMass_col=0._r8
  allocate(H2OLoss_CumYr_col(JY,JX));       H2OLoss_CumYr_col=0._r8
  allocate(QDrain_col(JY,JX));      QDrain_col=0._r8
  allocate(XGridSurfRunoff_2DH(2,2,JV,JH));      XGridSurfRunoff_2DH=0._r8
  allocate(HeatXGridBySurfRunoff_2DH(2,2,JV,JH));     HeatXGridBySurfRunoff_2DH=0._r8
  allocate(QRunSurf_col(JY,JX));        QRunSurf_col=0._r8
  allocate(QDischar_col(JY,JX));       QDischar_col=0._r8
  allocate(QflxSurfRunoffM_2DH(60,2,2,JV,JH)); QflxSurfRunoffM_2DH=0._r8
  allocate(TXGridSurfRunoff_2DH(JY,JX));         TXGridSurfRunoff_2DH=0._r8
  allocate(THeatXGridBySurfRunoff_2DH(JY,JX));        THeatXGridBySurfRunoff_2DH=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSoilWater
  use abortutils, only : destroy
  implicit none

  call destroy(iPondBotLev_col)
  call destroy(iPondFlag_col)
  call destroy(QWatIntLaterFlow_col)
  call destroy(Rain2ExposedSurf_col)
  call destroy(RainPrec2Sno_col)
  call destroy(RainPrecThrufall_col)
  call destroy(PrecipAtm2LandSurf_col)
  call destroy(QSnoIceXfer2Soil_col)
  call destroy(QSnoWatXfer2Soil_col)
  call destroy(Rain2Soil_col)
  call destroy(SoilWatMassBeg_col)
  call destroy(SoilWatMassEnd_col)
  call destroy(QdewCanopy_CumYr_pft)
  call destroy(ThetaAir_vr)
  call destroy(VLsoiAirP_vr)
  call destroy(THETW_vr)
  call destroy(THETI_vr)
  call destroy(ThetaH2OZ_vr)
  call destroy(ThetaICEZ_vr)
  call destroy(VLWatMicP_vr)
  call destroy(VLiceMicP_vr)
  call destroy(VLWatMacP_vr)
  call destroy(PSISoilMatricP_vr)
  call destroy(ElvAdjstedSoilH2OPSIMPa_vr)
  call destroy(VLWatMicPX_vr)
  call destroy(FWatExMacP2MicP_vr)
  call destroy(VLiceMacP_vr)
  call destroy(VLWatMicPM_vr)
  call destroy(VLWatMacPM_vr)
  call destroy(VLsoiAirPM_vr)
  call destroy(FILMM_vr)
  call destroy(WaterTBLSlope_col)
  call destroy(WtblDepzTile_col)
  call destroy(TileWaterTable_col)
  call destroy(DTBLD)
  call destroy(DepzIntWTBL_col)
  call destroy(ExtWaterTablet0_col)
  call destroy(ExtWaterTable_col)
  call destroy(NatWtblDepz_col)
  call destroy(EnergyImpact4ErosionM)
  call destroy(XVLMobileWaterLitRM)
  call destroy(XVLMobileWatMicPM)
  call destroy(XVLiceMicPM)
  call destroy(HydroCond_3D)
  call destroy(HydroCondMacP_vr)
  call destroy(HydroCondMicP4RootUptake_vr)
  call destroy(SurfRunoffWatFluxM_2DH)
  call destroy(RunoffVelocityM_col)
  call destroy(IFLBM)
  call destroy(XGridRunoffFlag)
  call destroy(IFLBH)
  call destroy(RechargNorthSubSurf)
  call destroy(RechargEastSubSurf)
  call destroy(RechargSouthSubSurf)
  call destroy(RechargWestSubSurf)
  call destroy(RechargRateNorthWTBL)
  call destroy(RechargRateEastWTBL)
  call destroy(RechargRateSouthWTBL)
  call destroy(RechargRateWestWTBL)
  call destroy(RechargNorthSurf)
  call destroy(RechargEastSurf)
  call destroy(RechargSouthSurf)
  call destroy(RechargWestSurf)
  call destroy(RechargBottom_col)
  call destroy(WaterFlow2MicPM_3D)
  call destroy(WaterFlow2MacPM_3D)
  call destroy(ReductVLsoiAirPM_vr)
  call destroy(FWatExMacP2MicPM_vr)
  call destroy(WatFlowSno2MicPM_col)
  call destroy(WatFlowSno2MacPM_col)
  call destroy(AirFilledSoilPoreM_vr)
  call destroy(TortMicPM_vr)
  call destroy(TortMacPM_vr)
  call destroy(DiffusivitySolutEffM_vr)
  call destroy(SoilResit4RootPentrate_vr)
  call destroy(PSISE_vr)
  call destroy(PSISoilAirEntry)
  call destroy(PSISoilOsmotic_vr)
  call destroy(PSIGrav_vr)
  call destroy(SoilWatAirDry_vr)
  call destroy(ThetaSat_vr)
  call destroy(WaterFlowSoiMicPX_3D)
  call destroy(QEvap_CumYr_col)
  call destroy(EvapoTransp_col)
  call destroy(QRain_CumYr_col)
  call destroy(Qrunoff_CumYr_col)
  call destroy(WatMass_col)
  call destroy(H2OLoss_CumYr_col)
  call destroy(QDrain_col)
  call destroy(XGridSurfRunoff_2DH)
  call destroy(HeatXGridBySurfRunoff_2DH)
  call destroy(QRunSurf_col)
  call destroy(QDischar_col)
  call destroy(QflxSurfRunoffM_2DH)
  call destroy(Qinflx2Soil_col)
  call destroy(TXGridSurfRunoff_2DH)
  call destroy(THeatXGridBySurfRunoff_2DH)

  end subroutine DestructSoilWater

end module SoilWaterDataType
