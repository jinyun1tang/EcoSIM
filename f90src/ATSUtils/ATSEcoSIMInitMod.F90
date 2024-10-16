module ATSEcoSIMInitMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilWaterDataType
  use SharedDataMod
  use GridDataType
  use SOMDataType
  USE SoilPhysDataType
  use LandSurfDataType
  use EcoSIMCtrlMod
  use HydroThermData, only : PSISM1_vr, TKSoi1_vr, VHeatCapacity1_vr, &
      SoilFracAsMicP_vr, VLWatMicP1_vr, VLiceMicP1_vr !need the only as some vars
  use CanopyDataType, only: RadSWGrnd_col
  use ClimForcDataType, only : LWRadSky, TairK_col, &
      VPA, WindSpeedAtm_col, RainH
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
  real(r8) :: DORGC(JZ),DVLiceMicP(JZ)
  real(r8) :: TXCO2(JY,JX),DORGE(JY,JX)
  real(r8) :: VOLISO,VOLPT,VOLTT
  real(r8) :: TFLWT

  NHW=1;NHE=1;NVN=1;NVS=NYS
  !Setting some variables
  !That ecosim needs to recognize that it is running in coupled mode
  !with ATS and to turn off features unsupported in the coupler
  ATS_cpl_mode=.true.
  plant_model=.false.
  microbial_model=.false.
  soichem_model=.false.
  snowRedist_model=.false.
  disp_planttrait=.false.
  disp_modelconfig=.false.

  !Calling some setup functions
  call SetMeshATS(NHW,NVN,NHE,NVS)
  call set_ecosim_solver(1, 1, 1, 1)
  call InitAlloc(NOMicrobeGuilds=1)

  !setting a few variables 
  FlowDirIndicator = 3 !Basically a switch, setting to 3 removes lateral flow
  MaxNumRootLays = 1 !Is the number of layers down the roots go
  NX=1

  VOLISO=0.0_r8
  TFLWT=0.0_r8
  VOLPT=0.0_r8
  VOLTT=0.0_r8
 
  do NY=1,NYS
    TXCO2(NY,NX)=0.0_r8
    DORGE(NY,NX)=0.0_r8
    QRunSurf_col(NY,NX)=0.0_r8
    NU(NY,NX)=a_NU(NY)
    NL(NY,NX)=a_NL(NY)
    AREA(3,0,NY,NX)=a_AREA3(0,NY)
    AREA(3,NU(NY,NX),NY,NX)=a_AREA3(0,NY)
    ASP_col(NY,NX)=a_ASP(NY)
    !TairKClimMean(NY,NX)=a_ATKA(NY)
    !CO2E_col(NY,NX)=atm_co2
    !CH4E_col(NY,NX)=atm_ch4
    !OXYE(NY,NX)=atm_o2
    !Z2GE(NY,NX)=atm_n2
    !Z2OE(NY,NX)=atm_n2o
    !ZNH3E(NY,NX)=atm_nh3
    !H2GE(NY,NX)=atm_H2
    TairK_col(NY,NX)=tairc(NY) !it's already in K??
    !convert VPA from ATS units (Pa) to EcoSIM (MPa)
    VPA(NY,NX) = vpair(NY)/1.0e6_r8
    !convert WindSpeedAtm from ATS units (m s^-1) to EcoSIM (m h^-1)
    WindSpeedAtm_col(NY,NX)  = uwind(NY)*3600.0_r8

    DO L=NU(NY,NX),NL(NY,NX)
      TKSoi1_vr(L,NY,NX) = a_TEMP(L,NY)
      CumDepz2LayerBot_vr(L,NY,NX)=a_CumDepz2LayerBot_vr(L,NY)
      POROS_vr(L,NY,NX)=a_PORO(L,NY)
      SoiBulkDensityt0_vr(L,NY,NX)=a_BKDSI(L,NY)
      CSoilOrgM_vr(ielmc,L,NY,NX)=a_CORGC(L,NY)
      CSoilOrgM_vr(ielmn,L,NY,NX)=a_CORGN(L,NY)
      CSoilOrgM_vr(ielmp,L,NY,NX)=a_CORGP(L,NY)
    ENDDO
  ENDDO

  PSIAtFldCapacity = pressure_at_field_capacity
  PSIAtWiltPoint = pressure_at_wilting_point

  call startsim(NHW,NHE,NVN,NVS)

  end subroutine Init_EcoSIM_Soil


end module ATSEcoSIMInitMod
