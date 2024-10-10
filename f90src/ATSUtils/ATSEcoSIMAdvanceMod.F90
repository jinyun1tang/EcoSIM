module ATSEcoSIMAdvanceMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilWaterDataType
  use SharedDataMod
  use GridDataType
  use GridConsts
  use SOMDataType
  USE SoilPhysDataType
  use LandSurfDataType
  use CanopyDataType, only: RadSWGrnd_col
  !use PlantAPIData, only: CO2E, CH4E, OXYE, Z2GE, Z2OE, ZNH3E, &
  !    H2GE
  use ClimForcDataType, only : LWRadSky, TairK_col, &
      VPA, WindSpeedAtm_col, RainH  
  use SoilPropertyDataType
  use HydroThermData, only : PSISM1_vr, TKSoi1_vr, VHeatCapacity1_vr, &
      SoilFracAsMicP_vr, VLWatMicP1_vr, VLiceMicP1_vr !need the only as some vars are double defined
  use EcoSIMSolverPar, only : NPH, dts_HeatWatTP
implicit none
  character(len=*), private, parameter :: mod_filename=&
  __FILE__
  public :: RunEcoSIMSurfaceBalance
  contains

  subroutine RunEcoSIMSurfaceBalance(NYS)
  use EcoSimConst
  use GridMod           , only : SetMeshATS
  use SurfPhysMod       , only : RunSurfacePhysModel, StageSurfacePhysModel
  use StartsMod         , only : set_ecosim_solver
  implicit none
  integer :: NY,NX,L,NHW,NHE,NVN,NVS, I, J, M, heat_vec_size
  integer, intent(in) :: NYS
  real(r8) :: YSIN(NumOfSkyAzimuSects),YCOS(NumOfSkyAzimuSects),SkyAzimuthAngle(NumOfSkyAzimuSects)
  real(r8) :: ResistanceLitRLay(JY,JX)
  real(r8) :: KSatReductByRainKineticEnergy(JY,JX)
  real(r8) :: HeatFluxAir2Soi(JY,JX)
  real(r8) :: TopLayWatVol(JY,JX)
  real(r8) :: Qinfl2MicP(JY,JX)
  real(r8) :: HInfl2Soil(JY,JX)

  NHW=1;NHE=1;NVN=1;NVS=NYS

  call SetMeshATS(NHW,NVN,NHE,NVS)

  NX=1

  do NY=1,NYS
    NU(NY,NX)=a_NU(NY)
    NL(NY,NX)=a_NL(NY)
    a_AREA3(0,NY) = 1.0_r8
    AREA(3,0,NY,NX)=a_AREA3(0,NY)
    AREA(3,NU(NY,NX),NY,NX)=a_AREA3(0,NY)
    AREA(3,2,NY,NX)=a_AREA3(0,NY)


    ASP_col(NY,NX)=a_ASP(NY)
    !TairKClimMean(NY,NX)=a_ATKA(NY)
    !CO2E(NY,NX)=atm_co2
    !CH4E(NY,NX)=atm_ch4
    !OXYE(NY,NX)=atm_o2
    !Z2GE(NY,NX)=atm_n2
    !Z2OE(NY,NX)=atm_n2o
    !ZNH3E(NY,NX)=atm_nh3
    !H2GE(NY,NX)=atm_H2
    TairK_col(NY,NX)=tairc(NY)
    !convert VPA from ATS units (Pa) to EcoSIM (MPa)
    VPA(NY,NX) = vpair(NY)/1.0e6_r8
    !convert WindSpeedAtm_col from ATS units (m s^-1) to EcoSIM (m h^-1)
    WindSpeedAtm_col(NY,NX) = uwind(NY)*3600.0_r8
    !converting radiation units from ATS (W m^-2) to EcoSIM (MJ m^-2 h^-1)
    RadSWGrnd_col(NY,NX) = swrad(NY)*0.0036_r8
    LWRadSky(NY,NX) = sunrad(NY)*0.0036_r8
    RainH(NY,NX) = p_rain(NY)
    DO L=NU(NY,NX),NL(NY,NX)
      CumDepz2LayerBot_vr(L,NY,NX)=a_CumDepz2LayerBot_vr(L,NY)
      !Convert Bulk Density from ATS (kg m^-3) to EcoSIM (Mg m^-3)
      SoiBulkDensityt0_vr(L,NY,NX) = a_BKDSI(L,NY)/1.0e3_r8
      CSoilOrgM_vr(ielmc,L,NY,NX)  = a_CORGC(L,NY)
      CSoilOrgM_vr(ielmn,L,NY,NX)  = a_CORGN(L,NY)
      CSoilOrgM_vr(ielmp,L,NY,NX)  = a_CORGP(L,NY)
      VLWatMicP1_vr(L,NY,NX)       = a_WC(L,NY)
      VLiceMicP1_vr(L,NY,NX)       = 0.0
      TKSoi1_vr(L,NY,NX)              = a_TEMP(L,NY)
      VHeatCapacity1_vr(L,NY,NX)   = heat_capacity
      SoilFracAsMicP_vr(L,NY,NX)   = 1.0
      PSISM1_vr(L,NY,NX)           = a_MATP(L,NY)
      POROS_vr(L,NY,NX)            = a_PORO(L,NY)
      !AREA3(L,NY,NX)              = a_AREA3(L,NY)
   ENDDO
   POROS_vr(0,NY,NX) = POROS_vr(1,NY,NX)
  ENDDO

  PSIAtFldCapacity = pressure_at_field_capacity
  PSIAtWiltPoint = pressure_at_wilting_point

  call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)

  DO M=1,NPH
    call RunSurfacePhysModel(I,J,M,NHE,NHW,NVS,NVN,ResistanceLitRLay,&
      KSatReductByRainKineticEnergy,TopLayWatVol,HeatFluxAir2Soi,Qinfl2MicP,Hinfl2Soil)
  ENDDO
  
  write(*,*) "Heat and water souces: "
  write(*,*) "Hinfl2Soil = ", Hinfl2Soil, " MJ"
  write(*,*) "HeatFluxAir2Soi = ", HeatFluxAir2Soi
  write(*,*) "Qinfl2MicP = ", Qinfl2MicP 
  write(*,*) "Timestep in EcoSIM: ", dts_HeatWatTP, " hr"
  DO NY=1,NYS
    write(*,*) "ATS true E source: ", Hinfl2Soil(NY,1) / (dts_HeatWatTP*3600._r8)
    !for every column send the top layer to the transfer var
    surf_e_source(NY) = Hinfl2Soil(NY,1) / (dts_HeatWatTP*3600._r8)
    surf_w_source(NY) = Qinfl2MicP(NY,1)
    write(*,*) "After conversion ", surf_e_source(NY) , " MJ/s" 
  ENDDO

  end subroutine RunEcoSIMSurfaceBalance

end module ATSEcoSIMAdvanceMod
