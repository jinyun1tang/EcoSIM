module ATSEcoSIMAdvanceMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilWaterDataType
  use SharedDataMod
  use GridDataType
  use GridConsts
  use SOMDataType
  USE SoilPhysDataType
  use SnowDataType
  use LandSurfDataType
  use CanopyDataType, only: RadSWGrnd_col
  !use PlantAPIData, only: CO2E, CH4E, OXYE, Z2GE, Z2OE, ZNH3E, &
  !    H2GE
  use ClimForcDataType 
  use SoilPropertyDataType
  use HydroThermData, only : PSISM1_vr, TKSoi1_vr, VHeatCapacity1_vr, &
      SoilFracAsMicP_vr, VLWatMicP1_vr, VLiceMicP1_vr, VLairMicP1_vr, &
      FracSoiPAsWat_vr, FracSoiPAsIce_vr, FracSoiPAsAir_vr !need the only as some vars are double defined
  use EcoSIMSolverPar, only : NPH, dts_HeatWatTP
  use UnitMod    , only : units
  use EcoSIMCtrlDataType
  use MiniMathMod
implicit none
  character(len=*), private, parameter :: mod_filename=&
  __FILE__
  public :: RunEcoSIMSurfaceBalance
  contains

  subroutine RunEcoSIMSurfaceBalance(NYS)
  use EcoSimConst
  use GridMod           , only : SetMeshATS
  use SurfPhysMod       , only : RunSurfacePhysModel, StageSurfacePhysModel, & 
      SetHourlyAccumulatorsATS
  use StartsMod         , only : set_ecosim_solver
  use SnowBalanceMod    , only : SnowMassUpdate
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
  !real(r8) :: SnowDepth_col(JY,JX)
  real(r8) :: PrecAsRain(JY,JX)
  real(r8) :: PrecAsSnow(JY,JX)
  real(r8) :: VLTSoiPore
  real(r8), PARAMETER :: TSNOW=-0.25_r8  !oC, threshold temperature for snowfall

  NHW=1;NHE=1;NVN=1;NVS=NYS

  call SetMeshATS(NHW,NVN,NHE,NVS)

  NX=1

  do NY=1, NYS
    call SetHourlyAccumulatorsATS(NY,NX)
  enddo

  do NY=1,NYS
    NU(NY,NX)=a_NU(NY)
    NL(NY,NX)=a_NL(NY)
    !a_AREA3(0,NY) = 1.0_r8
    !AREA(3,0,NY,NX)=a_AREA3(1,NY)
    !AREA(3,NU(NY,NX),NY,NX)=a_AREA3(1,NY)
    !AREA(3,2,NY,NX)=a_AREA3(1,NY)

    ASP_col(NY,NX)=a_ASP(NY)
    !TairKClimMean(NY,NX)=a_ATKA(NY)
    !CO2E_col(NY,NX)=atm_co2
    !CH4E_col(NY,NX)=atm_ch4
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
    !converting precipitation units from m s^-1 to mm hr^-1
    RainH(NY,NX) = p_rain(NY)*1000.0_r8*3600.0_r8
    TCA_col(NY,NX) = units%Kelvin2Celcius(TairK_col(NY,NX))
    DO L=NU(NY,NX),NL(NY,NX)
      CumDepz2LayerBot_vr(L,NY,NX)=a_CumDepz2LayerBot_vr(L,NY)
      AREA(3,L,NY,NX)=a_AREA3(L,NY)
      !write(*,*) "AREA(3,L,NY,NX) = ", AREA(3,L,NY,NX), ", a_AREA3(L,NY) = ", a_AREA3(L,NY)
      !Convert Bulk Density from ATS (kg m^-3) to EcoSIM (Mg m^-3)
      SoiBulkDensity_vr(L,NY,NX)=a_BKDSI(L,NY)/1.0e3_r8
      CSoilOrgM_vr(ielmc,L,NY,NX)=a_CORGC(L,NY)
      CSoilOrgM_vr(ielmn,L,NY,NX)=a_CORGN(L,NY)
      CSoilOrgM_vr(ielmp,L,NY,NX)=a_CORGP(L,NY)
      VLWatMicP1_vr(L,NY,NX)=a_WC(L,NY)
      VLiceMicP1_vr(L,NY,NX)=0.0
      TKSoi1_vr(L,NY,NX) = a_TEMP(L,NY)
      VHeatCapacity1_vr(L,NY,NX) = heat_capacity
      SoilFracAsMicP_vr(L,NY,NX) = 1.0
      PSISM1_vr(L,NY,NX) = a_MATP(L,NY)
      POROS_vr(L,NY,NX) = a_PORO(L,NY)
      !AREA3(L,NY,NX) = a_AREA3(L,NY)
      !Setting the soil moisture for snow-soil heat exchange
      VLTSoiPore = VLSoilMicP_vr(L,NY,NX)
      IF(VLTSoiPore.GT.ZEROS2(NY,NX))THEN
        !fraction as water
        FracSoiPAsWat_vr(L,NY,NX)=AZMAX1t(VLWatMicP1_vr(L,NY,NX)/VLTSoiPore)
        !fraction as ice
        FracSoiPAsIce_vr(L,NY,NX)=AZMAX1t(VLiceMicP1_vr(L,NY,NX)/VLTSoiPore)
        !fraction as air
        FracSoiPAsAir_vr(L,NY,NX)=AZMAX1t(VLairMicP1_vr(L,NY,NX)/VLTSoiPore)
      ELSE
        FracSoiPAsWat_vr(L,NY,NX)=POROS_vr(L,NY,NX)
        FracSoiPAsIce_vr(L,NY,NX)=0.0_r8
        FracSoiPAsAir_vr(L,NY,NX)=0.0_r8
      ENDIF
   ENDDO
    IF(TCA_col(NY,NX).GT.TSNOW)THEN
      PrecAsRain(NY,NX)=RAINH(NY,NX)
      PrecAsSnow(NY,NX)=0.0_r8
    ELSE
      PrecAsRain(NY,NX)=0.0_r8
      PrecAsSnow(NY,NX)=RAINH(NY,NX)
    ENDIF
    RainFalPrec(NY,NX)=PrecAsRain(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
    SnoFalPrec_col(NY,NX)=PrecAsSnow(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
    POROS_vr(0,NY,NX) = POROS_vr(1,NY,NX)
  ENDDO

  PSIAtFldCapacity = pressure_at_field_capacity
  PSIAtWiltPoint = pressure_at_wilting_point

  !write(*,*) "before stage: SnoFalPrec(NY,NX) = ", SnoFalPrec(1,1), " AREA(3,NU(NY,NX),NY,NX) = ", AREA(3,1,1,1)

  call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)

  !call UpdateSoilMoistureFromATS(I,J,NHW,NHE,NVN,NVS)

  !write(*,*) "after stage: SnoFalPrec(NY,NX) = ", SnoFalPrec(1,1), " AREA(3,NU(NY,NX),NY,NX) = ", AREA(3,1,1,1)

  DO M=1,NPH
    call RunSurfacePhysModel(I,J,M,NHE,NHW,NVS,NVN,ResistanceLitRLay,&
      KSatReductByRainKineticEnergy,TopLayWatVol,HeatFluxAir2Soi,Qinfl2MicP,Hinfl2Soil)
  ENDDO
 
  !write(*,*) "after run: SnoFalPrec(NY,NX) = ", SnoFalPrec(1,1), " AREA(3,NU(NY,NX),NY,NX) = ", AREA(3,1,1,1)

  !write(*,*) "Heat and water souces: "
  DO NY=1,NYS
    !for every column send the top layer to the transfer var
    !Convert heat and water flux from EcoSIM units (flux/ hr)
    !to ATS units (flux / s)
    surf_e_source(NY) = Hinfl2Soil(NY,1) / (dts_HeatWatTP*3600._r8)
    surf_w_source(NY) = Qinfl2MicP(NY,1) / (dts_HeatWatTP*3600._r8)
    !write(*,*) "After conversion ", surf_e_source(NY) , " MJ/s" 
    !write(*,*) "Water conversion ", surf_w_source(NY) , " m/s"
    surf_snow_depth(NY) = SnowDepth_col(NY,1)
  ENDDO

  write(*,*) "snow_depth = ", surf_snow_depth(1) 
  write(*,*) "Q_e = ", surf_e_source(1) , " MJ/s" 
  write(*,*) "Q_w ", surf_w_source(1) , " m/s"

  !open(unit=10, file="fluxes.txt", status="unknown", position="append")
  !write(10,*) "prec = ", RAINH(1,1), " m/s", " Q_e = ", surf_e_source(1), " MJ/s", " Q_w = ", surf_w_source(1), " m/s"
  !close(10)

  end subroutine RunEcoSIMSurfaceBalance

  subroutine UpdateSoilMoistureFromATS(I,J,NHW,NHE,NVN,NVS)
  implicit none
  
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX, L    
  real(r8) :: VLTSoiPore

  !Modified to remove all MacP1 as it should all be zero
  DO NY=1,NYS
    DO L=NU(NY,NX),NL(NY,NX)
      !VLTSoiPore = VLSoilMicP_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)
      VLTSoiPore = VLSoilMicP_vr(L,NY,NX)
      IF(VLTSoiPore.GT.ZEROS2(NY,NX))THEN
        !fraction as water
        FracSoiPAsWat_vr(L,NY,NX)=AZMAX1t(VLWatMicP1_vr(L,NY,NX)/VLTSoiPore)
        !fraction as ice
        FracSoiPAsIce_vr(L,NY,NX)=AZMAX1t(VLiceMicP1_vr(L,NY,NX)/VLTSoiPore)
        !fraction as air
        FracSoiPAsAir_vr(L,NY,NX)=AZMAX1t(VLairMicP1_vr(L,NY,NX)/VLTSoiPore)
      ELSE
        FracSoiPAsWat_vr(L,NY,NX)=POROS_vr(L,NY,NX)
        FracSoiPAsIce_vr(L,NY,NX)=0.0_r8
        FracSoiPAsAir_vr(L,NY,NX)=0.0_r8
      ENDIF
      !write(*,*) "NY, NX ", NY, NX
      !write(*,*) "FracSoiPAsWat_vr(L,NY,NX): ", FracSoiPAsWat_vr(L,NY,NX)
    ENDDO
  ENDDO
  end subroutine UpdateSoilMoistureFromATS

end module ATSEcoSIMAdvanceMod
