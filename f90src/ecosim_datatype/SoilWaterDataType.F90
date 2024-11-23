module SoilWaterDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8),target,allocatable ::  TXGridSurfRunoff_2DH(:,:)                           !
  real(r8),target,allocatable ::  THeatXGridBySurfRunoff_2DH(:,:)                          !

  real(r8),target,allocatable ::  ThetaAir_vr(:,:,:)                      !air concentration [m3 m-3]
  real(r8),target,allocatable ::  VLsoiAirP_vr(:,:,:)                       !soil air content [m3 d-2]
  real(r8),target,allocatable ::  THETW_vr(:,:,:)                      !volumetric water content [m3 m-3]
  real(r8),target,allocatable ::  THETI_vr(:,:,:)                      !volumetric ice content [m3 m-3]
  real(r8),target,allocatable ::  ThetaH2OZ_vr(:,:,:)                     !volumetric moblize water [m3 m-3]
  real(r8),target,allocatable ::  ThetaICEZ_vr(:,:,:)                     !volumetric mobile ice [m3 m-3]
  real(r8),target,allocatable ::  VLWatMicP_vr(:,:,:)                       !soil micropore water content [m3 d-2]
  real(r8),target,allocatable ::  VLiceMicP_vr(:,:,:)                       !soil micropore ice content   [m3 d-2]
  real(r8),target,allocatable ::  VLWatMacP_vr(:,:,:)                      !soil macropore water content [m3 d-2]
  real(r8),target,allocatable ::  PSISoilMatricP_vr(:,:,:)             !soil micropore matric water potential [MPa]
  real(r8),target,allocatable ::  TotalSoilH2OPSIMPa_vr(:,:,:)                      !soil micropore total water potential [MPa]
  real(r8),target,allocatable ::  VLWatMicPX_vr(:,:,:)                      !soil micropore water content before wetting front [m3 d-2]
  real(r8),target,allocatable ::  FWatExMacP2MicP(:,:,:)                       !soil macropore - micropore water transfer [m3 d-2 h-1]
  real(r8),target,allocatable ::  VLiceMacP_vr(:,:,:)                      !soil macropore ice content [m3 d-2]
  real(r8),target,allocatable ::  VLWatMicPM_vr(:,:,:,:)                    !soil micropore water content, [m3 d-2]
  real(r8),target,allocatable ::  VLWatMacPM(:,:,:,:)                   !soil macropore water content, [m3 d-2]
  real(r8),target,allocatable ::  VLsoiAirPM(:,:,:,:)                    !soil air content, [m3 d-2]
  real(r8),target,allocatable ::  FILM(:,:,:,:)                     !soil water film thickness , [m]
  real(r8),target,allocatable ::  WaterTBLSlope(:,:)                !slope of water table relative to surface slope, [-]
  real(r8),target,allocatable ::  WtblDepzTile_col(:,:)                       !depth of artificial water table
  real(r8),target,allocatable ::  DTBLY(:,:)                        !artificial water table depth, [m]
  real(r8),target,allocatable ::  DTBLD(:,:)                        !depth of artificial water table adjusted for elevation
  real(r8),target,allocatable ::  DepthInternalWTBL(:,:)            !internal water table depth, [m]
  real(r8),target,allocatable ::  ExtWaterTablet0(:,:)              !initial external water table depth, elevation corrected [m]
  real(r8),target,allocatable ::  ExtWaterTable_col(:,:)                !current external water table depth, elevation corrected [m]
  real(r8),target,allocatable ::  NatWtblDepz_col(:,:)                        !external water table depth, [m]
  real(r8),target,allocatable ::  EnergyImpact4ErosionM(:,:,:)                     !total energy impact for erosion
  real(r8),target,allocatable ::  XVLMobileWaterLitRM(:,:,:)                     !excess water+ice
  real(r8),target,allocatable ::  XVLMobileWatMicPM(:,:,:)                     !excess water
  real(r8),target,allocatable ::  XVLiceMicPM(:,:,:)                     !excess ice
  real(r8),target,allocatable ::  HydroCond_3D(:,:,:,:,:)                   !hydraulic conductivity at different moisture levels
  real(r8),target,allocatable ::  HydroCondMacP_vr(:,:,:)                       !macropore hydraulic conductivity, [m MPa-1 h-1]
  real(r8),target,allocatable ::  HydroCondMicP4RootUptake_vr(:,:,:)                       !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]
  real(r8),target,allocatable ::  WatFlux4ErosionM_2DH(:,:,:)                        !runoff water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  RunoffVelocity(:,:,:)                        !runoff velocity, [m t-1]
  integer,target,allocatable ::  IFLBM(:,:,:,:,:)                   !flag for directional surface runoff
  logical,target,allocatable ::  XGridRunoffFlag(:,:,:,:)           !enables or disables boundary water flux depending on aspect, [-]
  integer,target,allocatable ::  IFLBH(:,:,:,:)                     !flag for directional runoff, related to IFLBM
  real(r8),target,allocatable ::  RechargNorthSubSurf(:,:)                       !northern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargEastSubSurf(:,:)                       !eastern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargSouthSubSurf(:,:)                       !southern subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargWestSubSurf(:,:)                       !western subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargRateNorthWTBL(:,:)                       !northern subsurface boundary water flux rate constant, [h-1]
  real(r8),target,allocatable ::  RechargRateEastWTBL(:,:)                       !eastern subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RechargRateSouthWTBL(:,:)                       !southern subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RechargRateWestWTBL(:,:)                       !western subsurface boundary water flux  rate constant, [h-1]
  real(r8),target,allocatable ::  RechargNorthSurf(:,:)                        !northern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargEastSurf(:,:)                        !eastern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargSouthSurf(:,:)                        !southern surface boundary water flux , [-]
  real(r8),target,allocatable ::  RechargWestSurf(:,:)                        !western surface boundary water flux , [-]
  real(r8),target,allocatable ::  RCHGD(:,:)                        !lower subsurface boundary water flux , [-]
  real(r8),target,allocatable ::  WaterFlow2MicPM_3D(:,:,:,:,:)                   !micropore water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  WaterFlow2MacPM_3D(:,:,:,:,:)                  !macropore water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  ReductVLsoiAirPM(:,:,:,:)                     !change in soil air volume for layer from last to current iteration, [g d-2 t-1] >0, shrink
  real(r8),target,allocatable ::  FWatExMacP2MicPM(:,:,:,:)                    !soil macropore - micropore water transfer, [g d-2 t-1]
  real(r8),target,allocatable ::  WatFlowSno2MicPM(:,:,:)                      !meltwater flux into soil micropores
  real(r8),target,allocatable ::  WatFlowSno2MacPM(:,:,:)                      !meltwater flux into soil macropores
  real(r8),target,allocatable ::  THETPM(:,:,:,:)                   !soil air-filled porosity, [m3 m-3]
  real(r8),target,allocatable ::  TortMicPM_vr(:,:,:,:)                     !soil tortuosity, []
  real(r8),target,allocatable ::  TortMacPM(:,:,:,:)                    !macropore tortuosity, []
  real(r8),target,allocatable ::  DiffusivitySolutEff(:,:,:,:)                     !coefficient for dissolution - volatilization, []
  real(r8),target,allocatable ::  SoilResit4RootPentrate_vr(:,:,:)                       !soil hydraulic resistance, [MPa h m-2]
  real(r8),target,allocatable ::  PSISE_vr(:,:,:)                      !soil water potential at saturation, [Mpa]
  real(r8),target,allocatable ::  PSISoilAirEntry(:,:,:)                      !soil water potential at air entry, [Mpa]
  real(r8),target,allocatable ::  PSISoilOsmotic_vr(:,:,:)                      !osmotic soil water potential , [Mpa]
  real(r8),target,allocatable ::  PSIGrav_vr(:,:,:)                      !gravimetric soil water potential , [Mpa]
  real(r8),target,allocatable ::  THETY_vr(:,:,:)                      !air-dry water content, [m3 m-3]
  real(r8),target,allocatable ::  ThetaSat_vr(:,:,:)                      !micropore class water content
  real(r8),target,allocatable ::  WaterFlowSoiMicPX(:,:,:,:)                     !unsaturated water flux , [m3 d-2 h-1]
  real(r8),target,allocatable ::  EvapoTransp_col(:,:)              !evapotranspiration
  real(r8),target,allocatable ::  QEvap_CumYr_col(:,:)                        !total evaporation, [m3 d-2]
  real(r8),target,allocatable ::  QRain_CumYr_col(:,:)                        !total precipitation, [m3 d-2]
  real(r8),target,allocatable ::  Qrunoff_CumYr_col(:,:)                         !total surface runoff, [m3 d-2]
  real(r8),target,allocatable ::  WatMass_col(:,:)                !total soil water content, [m3 d-2]
  real(r8),target,allocatable ::  H2OLoss_CumYr_col(:,:)                        !total subsurface water flux, [m3 d-2]
  real(r8),target,allocatable ::  QDrain_col(:,:)                       !total water drainage below root zone, [m3 d-2]
  real(r8),target,allocatable ::  XGridSurfRunoff_2DH(:,:,:,:)                       !soil surface runoff water, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HeatXGridBySurfRunoff_2DH(:,:,:,:)                      !soil surface runoff heat, [MJ d-2 h-1]
  real(r8),target,allocatable ::  QRunSurf_col(:,:)                         !runoff from surface water, [m3 d-2 h-1]
  real(r8),target,allocatable ::  QDischar_col(:,:)                !water discharge, [m3 d-2 h-1]
  real(r8),target,allocatable ::  QflxSurfRunoffM_2DH(:,:,:,:,:)        !surface runoff,
  real(r8),target,allocatable ::  Qinflx2Soil_col(:,:)
  real(r8),target,allocatable :: QdewCanopy_CumYr_pft(:,:,:)
  private :: InitAllocate
  contains

  subroutine InitSoilWater
  implicit none

  call InitAllocate
  end subroutine InitSoilWater

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
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
  allocate(TotalSoilH2OPSIMPa_vr(0:JZ,JY,JX));  TotalSoilH2OPSIMPa_vr=0._r8
  allocate(VLWatMicPX_vr(0:JZ,JY,JX));  VLWatMicPX_vr=0._r8
  allocate(FWatExMacP2MicP(JZ,JY,JX));     FWatExMacP2MicP=0._r8
  allocate(VLiceMacP_vr(JZ,JY,JX));    VLiceMacP_vr=0._r8
  allocate(VLWatMicPM_vr(60,0:JZ,JY,JX));VLWatMicPM_vr=0._r8
  allocate(VLWatMacPM(60,JZ,JY,JX));VLWatMacPM=0._r8
  allocate(VLsoiAirPM(60,0:JZ,JY,JX));VLsoiAirPM=0._r8
  allocate(FILM(60,0:JZ,JY,JX));FILM=0._r8
  allocate(WaterTBLSlope(JY,JX));       WaterTBLSlope=0._r8
  allocate(WtblDepzTile_col(JY,JX));      WtblDepzTile_col=0._r8
  allocate(DTBLY(JY,JX));       DTBLY=0._r8
  allocate(DTBLD(JY,JX));       DTBLD=0._r8
  allocate(DepthInternalWTBL(JY,JX));       DepthInternalWTBL=0._r8
  allocate(ExtWaterTablet0(JY,JX));       ExtWaterTablet0=0._r8
  allocate(ExtWaterTable_col(JY,JX));       ExtWaterTable_col=0._r8
  allocate(NatWtblDepz_col(JY,JX));       NatWtblDepz_col=0._r8
  allocate(EnergyImpact4ErosionM(60,JY,JX));   EnergyImpact4ErosionM=0._r8
  allocate(XVLMobileWaterLitRM(60,JY,JX));   XVLMobileWaterLitRM=0._r8
  allocate(XVLMobileWatMicPM(60,JY,JX));   XVLMobileWatMicPM=0._r8
  allocate(XVLiceMicPM(60,JY,JX));   XVLiceMicPM=0._r8
  allocate(HydroCond_3D(3,100,0:JZ,JY,JX));HydroCond_3D=0._r8
  allocate(HydroCondMacP_vr(JZ,JY,JX));     HydroCondMacP_vr=0._r8
  allocate(HydroCondMicP4RootUptake_vr(JZ,JY,JX));     HydroCondMicP4RootUptake_vr=0._r8
  allocate(WatFlux4ErosionM_2DH(60,JV,JH));      WatFlux4ErosionM_2DH=0._r8
  allocate(RunoffVelocity(60,JY,JX));      RunoffVelocity=0._r8
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
  allocate(RCHGD(JY,JX));       RCHGD=0._r8
  allocate(WaterFlow2MicPM_3D(60,3,JD,JV,JH));WaterFlow2MicPM_3D=0._r8
  allocate(WaterFlow2MacPM_3D(60,3,JD,JV,JH));WaterFlow2MacPM_3D=0._r8
  allocate(ReductVLsoiAirPM(60,JZ,JY,JX));  ReductVLsoiAirPM=0._r8
  allocate(FWatExMacP2MicPM(60,JZ,JY,JX)); FWatExMacP2MicPM=0._r8
  allocate(WatFlowSno2MicPM(60,JY,JX));    WatFlowSno2MicPM=0._r8
  allocate(WatFlowSno2MacPM(60,JY,JX));    WatFlowSno2MacPM=0._r8
  allocate(THETPM(60,0:JZ,JY,JX));THETPM=0._r8
  allocate(TortMicPM_vr(60,0:JZ,JY,JX));TortMicPM_vr=0._r8
  allocate(TortMacPM(60,JZ,JY,JX)); TortMacPM=0._r8
  allocate(DiffusivitySolutEff(60,0:JZ,JY,JX));DiffusivitySolutEff=0._r8
  allocate(SoilResit4RootPentrate_vr(JZ,JY,JX));     SoilResit4RootPentrate_vr=0._r8
  allocate(PSISE_vr(0:JZ,JY,JX));  PSISE_vr=0._r8
  allocate(PSISoilAirEntry(0:JZ,JY,JX));  PSISoilAirEntry=0._r8
  allocate(PSISoilOsmotic_vr(0:JZ,JY,JX));  PSISoilOsmotic_vr=0._r8
  allocate(PSIGrav_vr(0:JZ,JY,JX));  PSIGrav_vr=0._r8
  allocate(THETY_vr(0:JZ,JY,JX));  THETY_vr=0._r8
  allocate(ThetaSat_vr(0:JZ,JY,JX));  ThetaSat_vr=0._r8
  allocate(WaterFlowSoiMicPX(3,JD,JV,JH));   WaterFlowSoiMicPX=0._r8
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
  call destroy(TotalSoilH2OPSIMPa_vr)
  call destroy(VLWatMicPX_vr)
  call destroy(FWatExMacP2MicP)
  call destroy(VLiceMacP_vr)
  call destroy(VLWatMicPM_vr)
  call destroy(VLWatMacPM)
  call destroy(VLsoiAirPM)
  call destroy(FILM)
  call destroy(WaterTBLSlope)
  call destroy(WtblDepzTile_col)
  call destroy(DTBLY)
  call destroy(DTBLD)
  call destroy(DepthInternalWTBL)
  call destroy(ExtWaterTablet0)
  call destroy(ExtWaterTable_col)
  call destroy(NatWtblDepz_col)
  call destroy(EnergyImpact4ErosionM)
  call destroy(XVLMobileWaterLitRM)
  call destroy(XVLMobileWatMicPM)
  call destroy(XVLiceMicPM)
  call destroy(HydroCond_3D)
  call destroy(HydroCondMacP_vr)
  call destroy(HydroCondMicP4RootUptake_vr)
  call destroy(WatFlux4ErosionM_2DH)
  call destroy(RunoffVelocity)
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
  call destroy(RCHGD)
  call destroy(WaterFlow2MicPM_3D)
  call destroy(WaterFlow2MacPM_3D)
  call destroy(ReductVLsoiAirPM)
  call destroy(FWatExMacP2MicPM)
  call destroy(WatFlowSno2MicPM)
  call destroy(WatFlowSno2MacPM)
  call destroy(THETPM)
  call destroy(TortMicPM_vr)
  call destroy(TortMacPM)
  call destroy(DiffusivitySolutEff)
  call destroy(SoilResit4RootPentrate_vr)
  call destroy(PSISE_vr)
  call destroy(PSISoilAirEntry)
  call destroy(PSISoilOsmotic_vr)
  call destroy(PSIGrav_vr)
  call destroy(THETY_vr)
  call destroy(ThetaSat_vr)
  call destroy(WaterFlowSoiMicPX)
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
