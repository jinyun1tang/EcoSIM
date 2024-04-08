module ATSEcoSIMInitMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilWaterDataType
  use SharedDataMod
  use GridDataType
  use SOMDataType
  USE SoilPhysDataType
  use LandSurfDataType
  use EcoSIMCtrlMod
  use HydroThermData, only : PSISM1, TKSoi1, VLHeatCapacity, &
      SoilFracAsMicP, VLWatMicP1, VLiceMicP1 !need the only as some vars
  use CanopyDataType, only: SWRadOnGrnd
  use ClimForcDataType, only : LWRadSky, TairK, &
      VPA, WindSpeedAtm, RainH
  use SoilPropertyDataType
implicit none
  character(len=*), private, parameter :: mod_filename=&
  __FILE__
  public :: Init_EcoSIM_Soil
  contains

  subroutine Init_EcoSIM_Soil(NYS)
  use EcoSimConst
  use GridMod           , only : SetMeshATS
  use InitAllocMod
  use StartsMod, only : startsim, set_ecosim_solver
  implicit none
  integer :: NY,NX,L,NHW,NHE,NVN,NVS
  integer, intent(in) :: NYS
  real(r8) :: YSIN(NumOfSkyAzimuSects),YCOS(NumOfSkyAzimuSects),SkyAzimuthAngle(NumOfSkyAzimuSects)

  NHW=1;NHE=1;NVN=1;NVS=NYS

  call SetMeshATS(NHW,NVN,NHE,NVS)
  call set_ecosim_solver(1, 1, 1, 1)
  call InitAlloc(NOMicrobeGuilds=1)

  !Setting some variables by hand as they are set to bad
  !Values by default in InitAlloc, should this change?
  FlowDirIndicator = 3 !Basically a switch, setting to 3 removes lateral flow
  MaxNumRootLays = 1 !Is the number of layers down the roots go
  NX=1
  ATS_cpl_mode=.true.

  do NY=1,NYS
    NU(NY,NX)=a_NU(NY)
    NL(NY,NX)=a_NL(NY)
    AREA(3,0,NY,NX)=a_AREA3(0,NY)
    AREA(3,NU(NY,NX),NY,NX)=a_AREA3(0,NY)
    ASP(NY,NX)=a_ASP(NY)
    !TairKClimMean(NY,NX)=a_ATKA(NY)
    !CO2E(NY,NX)=atm_co2
    !CH4E(NY,NX)=atm_ch4
    !OXYE(NY,NX)=atm_o2
    !Z2GE(NY,NX)=atm_n2
    !Z2OE(NY,NX)=atm_n2o
    !ZNH3E(NY,NX)=atm_nh3
    !H2GE(NY,NX)=atm_H2
    TairK(NY,NX)=tairc(NY) !it's already in K??
    VPA(NY,NX) = vpair(NY)
    WindSpeedAtm(NY,NX) = uwind(NY)  
    DO L=NU(NY,NX),NL(NY,NX)
      TKSoi1(L,NY,NX) = a_TEMP(L,NY)
      CumDepth2LayerBottom(L,NY,NX)=a_CumDepth2LayerBottom(L,NY)
      POROS(L,NY,NX)=a_PORO(L,NY)
      SoiBulkDensityt0(L,NY,NX)=a_BKDSI(L,NY)
      CORGC(L,NY,NX)=a_CORGC(L,NY)
      CORGN(L,NY,NX)=a_CORGN(L,NY)
      CORGP(L,NY,NX)=a_CORGP(L,NY)
    ENDDO
  ENDDO

  PSIAtFldCapacity = pressure_at_field_capacity
  PSIAtWiltPoint = pressure_at_wilting_point

  call startsim(NHW,NHE,NVN,NVS)

  end subroutine Init_EcoSIM_Soil


end module ATSEcoSIMInitMod
