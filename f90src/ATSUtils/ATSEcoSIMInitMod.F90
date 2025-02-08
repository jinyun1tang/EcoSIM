module ATSEcoSIMInitMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilWaterDataType
  use SharedDataMod
  use GridDataType
  use SOMDataType
  USE SoilPhysDataType
  use LandSurfDataType
  use EcoSIMCtrlMod
  use HydroThermData, only : PSISM1_vr, TKSoil1_vr, VHeatCapacity1_vr, &
      SoilFracAsMicP_vr, VLWatMicP1_vr, VLiceMicP1_vr, FracSoiPAsWat_vr, &
      FracSoiPAsIce_vr, AirFilledSoilPore_vr!need the only as some vars
  use CanopyDataType, only: RadSWGrnd_col
  use ClimForcDataType, only : LWRadSky_col, TairK_col, &
      VPA_col, WindSpeedAtm_col, RainH, VPK_col
  use SoilPropertyDataType
  use SurfLitterDataType
  use EcoSIMConfig
  use MiniMathMod
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
  real(r8) :: YSIN(NumOfSkyAzimuthSects),YCOS(NumOfSkyAzimuthSects),SkyAzimuthAngle(NumOfSkyAzimuthSects)
  real(r8) :: DORGC(JZ),DVLiceMicP(JZ)
  real(r8) :: TXCO2(JY,JX),DORGE(JY,JX)
  real(r8) :: VOLISO,VOLPT,VOLTT
  real(r8) :: TFLWT
  real(r8) :: VPS(JY,JX)

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
  column_mode=.true.

  !Calling some setup functions
  call SetMeshATS(NHW,NVN,NHE,NVS)
  call set_ecosim_solver(30, 10, 20, 20)
  call InitAlloc(NOMicrobeGuilds=1)

  !setting a few variables 
  FlowDirIndicator = 3 !Basically a switch, setting to 3 removes lateral flow
  MaxNumRootLays   = 1 !Is the number of layers down the roots go
  NX               = 1

  VOLISO = 0.0_r8
  TFLWT  = 0.0_r8
  VOLPT  = 0.0_r8
  VOLTT  = 0.0_r8
 
  do NY=1,NYS
    TXCO2(NY,NX)            = 0.0_r8
    DORGE(NY,NX)            = 0.0_r8
    QRunSurf_col(NY,NX)     = 0.0_r8
    NU(NY,NX)               = a_NU(NY)
    NL(NY,NX)               = a_NL(NY)
    !AREA(3,0,NY,NX)         = a_AREA3(0,NY)
    !AREA(3,NU(NY,NX),NY,NX) = a_AREA3(0,NY)
    ASP_col(NY,NX)          = a_ASP(NY)
    !TairKClimMean(NY,NX)   = a_ATKA(NY)
    !CO2E_col(NY,NX)=atm_co2
    !CH4E_col(NY,NX)=atm_ch4
    !OXYE_col(NY,NX)=atm_o2
    !Z2GE_col(NY,NX)=atm_n2
    !Z2OE_col(NY,NX)=atm_n2o
    !ZNH3E_col(NY,NX)=atm_nh3
    !H2GE_col(NY,NX)=atm_H2
    TairK_col(NY,NX)=tairc(NY) !it's already in K??
    !convert VPA from ATS units (Pa) to EcoSIM (MPa)
    !VPA(NY,NX) = vpair(NY)/1.0e6_r8
    !convert WindSpeedAtm from ATS units (m s^-1) to EcoSIM (m h^-1)
    !VPS(NY,NX)              = vapsat0(TairK_col(NY,NX))*EXP(-ALTI(NY,NX)/7272.0_r8)
    VPK_col(NY,NX)          = vpair(NY)/1.0e3 !vapor pressure in kPa
    !VPK_col(NY,NX)          = AMIN1(VPK_col(NY,NX),VPS(NY,NX))
    VPA_col(NY,NX)              = VPK_col(NY,NX)*2.173E-03_r8/TairK_col(NY,NX)
    
    WindSpeedAtm_col(NY,NX)  = uwind(NY)*3600.0_r8

    !Need to check if litter area is set or not
    if(VGeomLayer_vr(0,NY,NX).EQ.0.0)then
      VGeomLayer_vr(0,NY,NX) = 0.1
    endif
    POROS0_col(NY,NX) = 0.5
    

    DO L=NU(NY,NX),NL(NY,NX)
      TKSoil1_vr(L,NY,NX) = a_TEMP(L,NY)
      CumDepz2LayBottom_vr(L,NY,NX)=a_CumDepz2LayBottom_vr(L,NY)
      POROS_vr(L,NY,NX)=a_PORO(L,NY)
      AREA(3,L,NY,NX)=a_AREA3(L,NY)
      !write(*,*) "AREA(3,L,NY,NX) = ", AREA(3,L,NY,NX), ", a_AREA3(L,NY) = ", a_AREA3(L,NY)
      SoiBulkDensityt0_vr(L,NY,NX)=a_BKDSI(L,NY)
      SoilBulkDensity_vr(L,NY,NX)=a_BKDSI(L,NY)
      SoilFracAsMicP_vr(L,NY,NX) = 1.0
      CSoilOrgM_vr(ielmc,L,NY,NX)=a_CORGC(L,NY)
      CSoilOrgM_vr(ielmn,L,NY,NX)=a_CORGN(L,NY)
      CSoilOrgM_vr(ielmp,L,NY,NX)=a_CORGP(L,NY)
      
      DH(NY,NX) = 1.0
      DV(NY,NX) = 1.0

    ENDDO
  ENDDO

  PSIAtFldCapacity = pressure_at_field_capacity
  PSIAtWiltPoint = pressure_at_wilting_point

  call startsim(NHW,NHE,NVN,NVS)

  do NY=1,NYS
    DO L=NU(NY,NX),NL(NY,NX)
      AREA(3,L,NY,NX)=a_AREA3(L,NY)
    ENDDO
  ENDDO

  end subroutine Init_EcoSIM_Soil


end module ATSEcoSIMInitMod
