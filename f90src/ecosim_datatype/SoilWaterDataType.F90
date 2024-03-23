module SoilWaterDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  THETP(:,:,:)                      !air concentration [m3 m-3]
  real(r8),target,allocatable ::  VLsoiAirP(:,:,:)                       !soil air content [m3 d-2]
  real(r8),target,allocatable ::  THETW(:,:,:)                      !volumetric water content [m3 m-3]
  real(r8),target,allocatable ::  THETI(:,:,:)                      !volumetric ice content [m3 m-3]
  real(r8),target,allocatable ::  THETWZ(:,:,:)                     !volumetric moblize water [m3 m-3]
  real(r8),target,allocatable ::  THETIZ(:,:,:)                     !volumetric mobile ice [m3 m-3]
  real(r8),target,allocatable ::  VLWatMicP(:,:,:)                       !soil micropore water content [m3 d-2]
  real(r8),target,allocatable ::  VLiceMicP(:,:,:)                       !soil micropore ice content   [m3 d-2]
  real(r8),target,allocatable ::  VLWatMacP(:,:,:)                      !soil macropore water content [m3 d-2]
  real(r8),target,allocatable ::  PSISoilMatricP(:,:,:)             !soil micropore matric water potential [MPa]
  real(r8),target,allocatable ::  TotalSoilH2OPSIMPa(:,:,:)                      !soil micropore total water potential [MPa]
  real(r8),target,allocatable ::  VLWatMicPX(:,:,:)                      !soil micropore water content before wetting front [m3 d-2]
  real(r8),target,allocatable ::  FWatExMacP2MicP(:,:,:)                       !soil macropore - micropore water transfer [m3 d-2 h-1]
  real(r8),target,allocatable ::  VLiceMacP(:,:,:)                      !soil macropore ice content [m3 d-2]
  real(r8),target,allocatable ::  VLWatMicPM(:,:,:,:)                    !soil micropore water content, [m3 d-2]
  real(r8),target,allocatable ::  VLWatMacPM(:,:,:,:)                   !soil macropore water content, [m3 d-2]
  real(r8),target,allocatable ::  VLsoiAirPM(:,:,:,:)                    !soil air content, [m3 d-2]
  real(r8),target,allocatable ::  FILM(:,:,:,:)                     !soil water film thickness , [m]
  real(r8),target,allocatable ::  WaterTBLSlope(:,:)                !slope of water table relative to surface slope, [-]
  real(r8),target,allocatable ::  DTBLDI(:,:)                       !depth of artificial water table
  real(r8),target,allocatable ::  DTBLY(:,:)                        !artificial water table depth, [m]
  real(r8),target,allocatable ::  DTBLD(:,:)                        !depth of artificial water table adjusted for elevation
  real(r8),target,allocatable ::  DepthInternalWTBL(:,:)            !internal water table depth, [m]
  real(r8),target,allocatable ::  ExtWaterTablet0(:,:)              !initial external water table depth, elevation corrected [m]
  real(r8),target,allocatable ::  ExtWaterTable(:,:)                !current external water table depth, elevation corrected [m]
  real(r8),target,allocatable ::  DTBLI(:,:)                        !external water table depth, [m]
  real(r8),target,allocatable ::  EnergyImpact4ErosionM(:,:,:)                     !total energy impact for erosion
  real(r8),target,allocatable ::  XVLMobileWaterLitRM(:,:,:)                     !excess water+ice
  real(r8),target,allocatable ::  XVLMobileWatMicPM(:,:,:)                     !excess water
  real(r8),target,allocatable ::  XVLiceMicPM(:,:,:)                     !excess ice
  real(r8),target,allocatable ::  HydroCond3D(:,:,:,:,:)                   !hydraulic conductivity at different moisture levels
  real(r8),target,allocatable ::  HydroCondMacP(:,:,:)                       !macropore hydraulic conductivity, [m MPa-1 h-1]
  real(r8),target,allocatable ::  HydroCondMicP4RootUptake(:,:,:)                       !soil micropore hydraulic conductivity for root water uptake [m MPa-1 h-1]
  real(r8),target,allocatable ::  WatFlux4ErosionM(:,:,:)                        !runoff water flux, [m3 d-2 t-1]
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
  real(r8),target,allocatable ::  WaterFlow2MicPM(:,:,:,:,:)                   !micropore water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  WaterFlow2MacPM(:,:,:,:,:)                  !macropore water flux, [m3 d-2 t-1]
  real(r8),target,allocatable ::  ReductVLsoiAirPM(:,:,:,:)                     !change in soil air volume for layer from last to current iteration, [g d-2 t-1] >0, shrink
  real(r8),target,allocatable ::  FWatExMacP2MicPM(:,:,:,:)                    !soil macropore - micropore water transfer, [g d-2 t-1]
  real(r8),target,allocatable ::  WatFlowSno2MicPM(:,:,:)                      !meltwater flux into soil micropores
  real(r8),target,allocatable ::  WatFlowSno2MacPM(:,:,:)                      !meltwater flux into soil macropores
  real(r8),target,allocatable ::  THETPM(:,:,:,:)                   !soil air-filled porosity, [m3 m-3]
  real(r8),target,allocatable ::  TortMicPM(:,:,:,:)                     !soil tortuosity, []
  real(r8),target,allocatable ::  TortMacPM(:,:,:,:)                    !macropore tortuosity, []
  real(r8),target,allocatable ::  DiffusivitySolutEff(:,:,:,:)                     !coefficient for dissolution - volatilization, []
  real(r8),target,allocatable ::  SoilResit4RootPentrate_vr(:,:,:)                       !soil hydraulic resistance, [MPa h m-2]
  real(r8),target,allocatable ::  PSISE(:,:,:)                      !soil water potential at saturation, [Mpa]
  real(r8),target,allocatable ::  PSISoilAirEntry(:,:,:)                      !soil water potential at air entry, [Mpa]
  real(r8),target,allocatable ::  PSISoilOsmotic(:,:,:)                      !osmotic soil water potential , [Mpa]
  real(r8),target,allocatable ::  PSIGrav(:,:,:)                      !gravimetric soil water potential , [Mpa]
  real(r8),target,allocatable ::  THETY(:,:,:)                      !air-dry water content, [m3 m-3]
  real(r8),target,allocatable ::  Theta_sat(:,:,:)                      !micropore class water content
  real(r8),target,allocatable ::  WaterFlowSoiMicPX(:,:,:,:)                     !unsaturated water flux , [m3 d-2 h-1]
  real(r8),target,allocatable ::  UEVAP(:,:)                        !total evaporation, [m3 d-2]
  real(r8),target,allocatable ::  URAIN(:,:)                        !total precipitation, [m3 d-2]
  real(r8),target,allocatable ::  URUN(:,:)                         !total surface runoff, [m3 d-2]
  real(r8),target,allocatable ::  UVLWatMicP(:,:)                        !total soil water content, [m3 d-2]
  real(r8),target,allocatable ::  UVOLO(:,:)                        !total subsurface water flux, [m3 d-2]
  real(r8),target,allocatable ::  UDRAIN(:,:)                       !total water drainage below root zone, [m3 d-2]
  real(r8),target,allocatable ::  Wat2GridBySurfRunoff(:,:,:,:)                       !soil surface runoff water, [m3 d-2 h-1]
  real(r8),target,allocatable ::  Heat2GridBySurfRunoff(:,:,:,:)                      !soil surface runoff heat, [MJ d-2 h-1]
  real(r8),target,allocatable ::  WQRH(:,:)                         !runoff from surface water, [m3 d-2 h-1]
  real(r8),target,allocatable ::  FWatDischarge(:,:)                !water discharge, [m3 d-2 h-1]
  real(r8),target,allocatable ::  QflxSurfRunoffM(:,:,:,:,:)        !surface runoff,
  real(r8),target,allocatable ::  Qinflx2Soil_col(:,:)
  private :: InitAllocate
  contains

  subroutine InitSoilWater
  implicit none

  call InitAllocate
  end subroutine InitSoilWater

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(Qinflx2Soil_col(JY,JX)); Qinflx2Soil_col=0._r8
  allocate(THETP(0:JZ,JY,JX));  THETP=0._r8
  allocate(VLsoiAirP(0:JZ,JY,JX));   VLsoiAirP=0._r8
  allocate(THETW(0:JZ,JY,JX));  THETW=0._r8
  allocate(THETI(0:JZ,JY,JX));  THETI=0._r8
  allocate(THETWZ(0:JZ,JY,JX)); THETWZ=0._r8
  allocate(THETIZ(0:JZ,JY,JX)); THETIZ=0._r8
  allocate(VLWatMicP(0:JZ,JY,JX));   VLWatMicP=0._r8
  allocate(VLiceMicP(0:JZ,JY,JX));   VLiceMicP=0._r8
  allocate(VLWatMacP(JZ,JY,JX));    VLWatMacP=0._r8
  allocate(PSISoilMatricP(0:JZ,JY,JX));  PSISoilMatricP=0._r8
  allocate(TotalSoilH2OPSIMPa(0:JZ,JY,JX));  TotalSoilH2OPSIMPa=0._r8
  allocate(VLWatMicPX(0:JZ,JY,JX));  VLWatMicPX=0._r8
  allocate(FWatExMacP2MicP(JZ,JY,JX));     FWatExMacP2MicP=0._r8
  allocate(VLiceMacP(JZ,JY,JX));    VLiceMacP=0._r8
  allocate(VLWatMicPM(60,0:JZ,JY,JX));VLWatMicPM=0._r8
  allocate(VLWatMacPM(60,JZ,JY,JX));VLWatMacPM=0._r8
  allocate(VLsoiAirPM(60,0:JZ,JY,JX));VLsoiAirPM=0._r8
  allocate(FILM(60,0:JZ,JY,JX));FILM=0._r8
  allocate(WaterTBLSlope(JY,JX));       WaterTBLSlope=0._r8
  allocate(DTBLDI(JY,JX));      DTBLDI=0._r8
  allocate(DTBLY(JY,JX));       DTBLY=0._r8
  allocate(DTBLD(JY,JX));       DTBLD=0._r8
  allocate(DepthInternalWTBL(JY,JX));       DepthInternalWTBL=0._r8
  allocate(ExtWaterTablet0(JY,JX));       ExtWaterTablet0=0._r8
  allocate(ExtWaterTable(JY,JX));       ExtWaterTable=0._r8
  allocate(DTBLI(JY,JX));       DTBLI=0._r8
  allocate(EnergyImpact4ErosionM(60,JY,JX));   EnergyImpact4ErosionM=0._r8
  allocate(XVLMobileWaterLitRM(60,JY,JX));   XVLMobileWaterLitRM=0._r8
  allocate(XVLMobileWatMicPM(60,JY,JX));   XVLMobileWatMicPM=0._r8
  allocate(XVLiceMicPM(60,JY,JX));   XVLiceMicPM=0._r8
  allocate(HydroCond3D(3,100,0:JZ,JY,JX));HydroCond3D=0._r8
  allocate(HydroCondMacP(JZ,JY,JX));     HydroCondMacP=0._r8
  allocate(HydroCondMicP4RootUptake(JZ,JY,JX));     HydroCondMicP4RootUptake=0._r8
  allocate(WatFlux4ErosionM(60,JV,JH));      WatFlux4ErosionM=0._r8
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
  allocate(WaterFlow2MicPM(60,3,JD,JV,JH));WaterFlow2MicPM=0._r8
  allocate(WaterFlow2MacPM(60,3,JD,JV,JH));WaterFlow2MacPM=0._r8
  allocate(ReductVLsoiAirPM(60,JZ,JY,JX));  ReductVLsoiAirPM=0._r8
  allocate(FWatExMacP2MicPM(60,JZ,JY,JX)); FWatExMacP2MicPM=0._r8
  allocate(WatFlowSno2MicPM(60,JY,JX));    WatFlowSno2MicPM=0._r8
  allocate(WatFlowSno2MacPM(60,JY,JX));    WatFlowSno2MacPM=0._r8
  allocate(THETPM(60,0:JZ,JY,JX));THETPM=0._r8
  allocate(TortMicPM(60,0:JZ,JY,JX));TortMicPM=0._r8
  allocate(TortMacPM(60,JZ,JY,JX)); TortMacPM=0._r8
  allocate(DiffusivitySolutEff(60,0:JZ,JY,JX));DiffusivitySolutEff=0._r8
  allocate(SoilResit4RootPentrate_vr(JZ,JY,JX));     SoilResit4RootPentrate_vr=0._r8
  allocate(PSISE(0:JZ,JY,JX));  PSISE=0._r8
  allocate(PSISoilAirEntry(0:JZ,JY,JX));  PSISoilAirEntry=0._r8
  allocate(PSISoilOsmotic(0:JZ,JY,JX));  PSISoilOsmotic=0._r8
  allocate(PSIGrav(0:JZ,JY,JX));  PSIGrav=0._r8
  allocate(THETY(0:JZ,JY,JX));  THETY=0._r8
  allocate(Theta_sat(0:JZ,JY,JX));  Theta_sat=0._r8
  allocate(WaterFlowSoiMicPX(3,JD,JV,JH));   WaterFlowSoiMicPX=0._r8
  allocate(UEVAP(JY,JX));       UEVAP=0._r8
  allocate(URAIN(JY,JX));       URAIN=0._r8
  allocate(URUN(JY,JX));        URUN=0._r8
  allocate(UVLWatMicP(JY,JX));       UVLWatMicP=0._r8
  allocate(UVOLO(JY,JX));       UVOLO=0._r8
  allocate(UDRAIN(JY,JX));      UDRAIN=0._r8
  allocate(Wat2GridBySurfRunoff(2,2,JV,JH));      Wat2GridBySurfRunoff=0._r8
  allocate(Heat2GridBySurfRunoff(2,2,JV,JH));     Heat2GridBySurfRunoff=0._r8
  allocate(WQRH(JY,JX));        WQRH=0._r8
  allocate(FWatDischarge(JY,JX));       FWatDischarge=0._r8
  allocate(QflxSurfRunoffM(60,2,2,JV,JH)); QflxSurfRunoffM=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSoilWater
  use abortutils, only : destroy
  implicit none
  call destroy(THETP)
  call destroy(VLsoiAirP)
  call destroy(THETW)
  call destroy(THETI)
  call destroy(THETWZ)
  call destroy(THETIZ)
  call destroy(VLWatMicP)
  call destroy(VLiceMicP)
  call destroy(VLWatMacP)
  call destroy(PSISoilMatricP)
  call destroy(TotalSoilH2OPSIMPa)
  call destroy(VLWatMicPX)
  call destroy(FWatExMacP2MicP)
  call destroy(VLiceMacP)
  call destroy(VLWatMicPM)
  call destroy(VLWatMacPM)
  call destroy(VLsoiAirPM)
  call destroy(FILM)
  call destroy(WaterTBLSlope)
  call destroy(DTBLDI)
  call destroy(DTBLY)
  call destroy(DTBLD)
  call destroy(DepthInternalWTBL)
  call destroy(ExtWaterTablet0)
  call destroy(ExtWaterTable)
  call destroy(DTBLI)
  call destroy(EnergyImpact4ErosionM)
  call destroy(XVLMobileWaterLitRM)
  call destroy(XVLMobileWatMicPM)
  call destroy(XVLiceMicPM)
  call destroy(HydroCond3D)
  call destroy(HydroCondMacP)
  call destroy(HydroCondMicP4RootUptake)
  call destroy(WatFlux4ErosionM)
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
  call destroy(WaterFlow2MicPM)
  call destroy(WaterFlow2MacPM)
  call destroy(ReductVLsoiAirPM)
  call destroy(FWatExMacP2MicPM)
  call destroy(WatFlowSno2MicPM)
  call destroy(WatFlowSno2MacPM)
  call destroy(THETPM)
  call destroy(TortMicPM)
  call destroy(TortMacPM)
  call destroy(DiffusivitySolutEff)
  call destroy(SoilResit4RootPentrate_vr)
  call destroy(PSISE)
  call destroy(PSISoilAirEntry)
  call destroy(PSISoilOsmotic)
  call destroy(PSIGrav)
  call destroy(THETY)
  call destroy(Theta_sat)
  call destroy(WaterFlowSoiMicPX)
  call destroy(UEVAP)
  call destroy(URAIN)
  call destroy(URUN)
  call destroy(UVLWatMicP)
  call destroy(UVOLO)
  call destroy(UDRAIN)
  call destroy(Wat2GridBySurfRunoff)
  call destroy(Heat2GridBySurfRunoff)
  call destroy(WQRH)
  call destroy(FWatDischarge)
  call destroy(QflxSurfRunoffM)
  call destroy(Qinflx2Soil_col)
  end subroutine DestructSoilWater

end module SoilWaterDataType
