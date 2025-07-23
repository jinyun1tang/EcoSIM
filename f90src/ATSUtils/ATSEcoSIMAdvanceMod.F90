module ATSEcoSIMAdvanceMod
  !
  !Description
  
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
  use ClimForcDataType, only : LWRadSky_col, TairK_col, &
      VPA_col, WindSpeedAtm_col, RainH, VPK_col  
  use SoilPropertyDataType
  use HydroThermData, only : PSISM1_vr, TKSoil1_vr, VHeatCapacity1_vr, &
      SoilFracAsMicP_vr, VLWatMicP1_vr, VLiceMicP1_vr, FracSoiPAsWat_vr, &
      FracSoiPAsIce_vr, FracAirFilledSoilPore_vr, VLairMicP1_vr!need the only as some vars
  use EcoSIMSolverPar, only : NPH, dts_HeatWatTP
  use UnitMod    , only : units
  use EcoSIMCtrlDataType
  use MiniMathMod
  use ClimForcDataType

implicit none
  character(len=*), private, parameter :: mod_filename=&
  __FILE__
  public :: RunEcoSIMSurfaceBalance
  contains

  ! Function to check for NaN in an array
  function is_nan(x) result(mask)
    real(r8), intent(in) :: x(:)
    logical, dimension(size(x)) :: mask
    integer :: i

    !allocate(mask(size(x)))
    do i = 1, size(x)
      mask(i) = (x(i) /= x(i))  ! NaN is the only value that is not equal to itself
    end do
  end function is_nan

  subroutine RunEcoSIMSurfaceBalance(NYS)
  !
  use EcoSimConst
  use GridMod           , only : SetMeshATS
  use SurfPhysMod       , only : RunSurfacePhysModelM, StageSurfacePhysModel, UpdateSurfaceAtM, &
      SetHourlyAccumulatorsATS
  use SnowBalanceMod    , only : SnowMassUpdate
  use StartsMod         , only : set_ecosim_solver
  use SnowBalanceMod    , only : SnowMassUpdate, SnowpackLayering
  implicit none
  integer, intent(in) :: NYS  !Number of columns?

  integer :: NY,NX,L,NHW,NHE,NVN,NVS, I, J, M, heat_vec_size, NPH_Test 
  real(r8) :: Wat_next
  real(r8) :: YSIN(NumOfSkyAzimuthSects),YCOS(NumOfSkyAzimuthSects),SkyAzimuthAngle(NumOfSkyAzimuthSects)
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
  real(r8) :: VPS(JY,JX)
  real(r8) :: EMM
  real(r8), PARAMETER :: TSNOW=-0.25_r8  !oC, threshold temperature for snowfall
  real(r8) :: Qinfl2MicPM(JY,JX)
  real(r8) :: Hinfl2SoilM(JY,JX)
  real(r8) :: VLWat_test(JZ,JY,JX)  

  NHW=1;NHE=1;NVN=1;NVS=NYS
  I=1;J=1
  NPH_Test=1
  NX=1

  call SetMeshATS(NHW,NVN,NHE,NVS)

  NX=1

  do NY=1, NYS
    call SetHourlyAccumulatorsATS(NY,NX)
  enddo

  do NY=1,NYS
    NU_col(NY,NX)               = a_NU(NY)
    NL_col(NY,NX)               = a_NL(NY)

    ASP_col(NY,NX)=a_ASP(NY)
    !TairKClimMean_col(NY,NX) = a_ATKA(NY)
    !CO2E_col(NY,NX)      = atm_co2
    !CH4E_col(NY,NX)      = atm_ch4
    !OXYE_col(NY,NX)      = atm_o2
    !Z2GE_col(NY,NX)      = atm_n2
    !Z2OE_col(NY,NX)      = atm_n2o
    !ZNH3E_col(NY,NX)     = atm_nh3
    !H2GE_col(NY,NX)      = atm_H2
    TairK_col(NY,NX)      = tairc(NY)
    !convert VPA from ATS units (Pa) to EcoSIM (MPa)
    !VPA(NY,NX) = vpair(NY)/1.0e6_r8

    !VPS(NY,NX)              = vapsat0(TairK_col(NY,NX))*EXP(-ALTI_col(NY,NX)/7272.0_r8)
    VPK_col(NY,NX)          = vpair(NY)/1.0e3 !vapor pressure in kPa
    !VPK_col(NY,NX)          = AMIN1(VPK_col(NY,NX),VPS(NY,NX))
    VPA_col(NY,NX)              = VPK_col(NY,NX)*2.173E-03_r8/TairK_col(NY,NX)
    !convert WindSpeedAtm_col from ATS units (m s^-1) to EcoSIM (m h^-1)
    WindSpeedAtm_col(NY,NX) = uwind(NY)*3600.0_r8
    !converting radiation units from ATS (W m^-2) to EcoSIM (MJ m^-2 h^-1)
    RadSWGrnd_col(NY,NX) = 0.0

    !EMM = 2.445 !There is a more elaborate calcuation of sky emissivity but I don't think we'll need that yet
    EMM = 0.684
    SkyLonwRad_col(NY,NX) = EMM*stefboltz_const*TairK_col(NY,NX)**4._r8
    LWRadSky_col(NY,NX) = SkyLonwRad_col(NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    TCA_col(NY,NX) = units%Kelvin2Celcius(TairK_col(NY,NX))
    DO L=NU_col(NY,NX),NL_col(NY,NX)
      CumDepz2LayBottom_vr(L,NY,NX) = a_CumDepz2LayBottom_vr(L,NY)
      !Convert Bulk Density from ATS (kg m^-3) to EcoSIM (Mg m^-3)
      SoiBulkDensityt0_vr(L,NY,NX) = a_BKDSI(L,NY)/1.0e3_r8
      CSoilOrgM_vr(ielmc,L,NY,NX)  = a_CORGC(L,NY)
      CSoilOrgM_vr(ielmn,L,NY,NX)  = a_CORGN(L,NY)
      CSoilOrgM_vr(ielmp,L,NY,NX)  = a_CORGP(L,NY)
      !Convert ATS units (mols) to EcoSIM units (mH2O)
      VLWatMicP1_vr(L,NY,NX)       = a_WC(L,NY)/(a_LDENS(L,NY)*AREA_3D(3,NU_col(NY,NX),NY,NX))
      !VLWatMicP1_vr(L,NY,NX)      = a_WC(L,NY)/(a_LDENS(L,NY))
      VLiceMicP1_vr(L,NY,NX)       = 0.0
      TKSoil1_vr(L,NY,NX)          = a_TEMP(L,NY)
      VHeatCapacity1_vr(L,NY,NX)   = heat_capacity
      SoilFracAsMicP_vr(L,NY,NX)   = 1.0
      !Convert Matric Pressure from ATS [Pa] to EcoSIM [MPa]
      PSISM1_vr(L,NY,NX)           = a_MATP(L,NY)/1.0e6
      POROS_vr(L,NY,NX)            = a_PORO(L,NY)
      !AREA3(L,NY,NX)              = a_AREA3(L,NY)
      VLTSoiPore = VLSoilMicP_vr(L,NY,NX)
      IF(VLTSoiPore.GT.ZEROS2(NY,NX))THEN
        !fraction as water
        FracSoiPAsWat_vr(L,NY,NX)=AZMAX1t(VLWatMicP1_vr(L,NY,NX)/VLTSoiPore)
        !fraction as ice
        FracSoiPAsIce_vr(L,NY,NX)=AZMAX1t(VLiceMicP1_vr(L,NY,NX)/VLTSoiPore)
        !fraction as air
        FracAirFilledSoilPore_vr(L,NY,NX)=AZMAX1t(VLairMicP1_vr(L,NY,NX)/VLTSoiPore)
      ELSE
        FracSoiPAsWat_vr(L,NY,NX)=POROS_vr(L,NY,NX)
        FracSoiPAsIce_vr(L,NY,NX)=0.0_r8
        FracAirFilledSoilPore_vr(L,NY,NX)=0.0_r8
      ENDIF    
    ENDDO

    !Setting the snow, if passed as total precipitation do full temp calc, else
    !Just set snow and rain to their given values
    !5/16 - I don't think we need to account for area here at all
    !       I am switching out PrecAsRain and just directly substituting
    !       with RainFalPrec
    if(p_bool)then
      IF(TCA_col(NY,NX).GT.TSNOW)THEN
        RainFalPrec_col(NY,NX)=p_total(NY)*3600.0 !convert from m/s to m/hr
        SnoFalPrec_col(NY,NX)=0.0_r8
      ELSE
        RainFalPrec_col(NY,NX)=0.0_r8
        SnoFalPrec_col(NY,NX)=p_total(NY)*3600.0 !convert to m SWE/s to m SWE/hr
      ENDIF
    else
      RainFalPrec_col(NY,NX)=p_rain(NY)*3600.0 !convert from m/s to m/hr
      SnoFalPrec_col(NY,NX)=p_snow(NY)*3600.0 !convert from m SWE/s to mSWE/hr 
    endif
    !Set Prec equal to variables for after Irrigation and canopy processing
    !since that is not included yet
    PrecRainAndIrrig_col = RainFalPrec_col
    RainPrecThrufall_col = RainFalPrec_col

    !RainFalPrec_col(NY,NX)=PrecAsRain(NY,NX)*AREA(3,NU_col(NY,NX),NY,NX)
    !SnoFalPrec_col(NY,NX)=PrecAsSnow(NY,NX)*AREA(3,NU_col(NY,NX),NY,NX)
    POROS_vr(0,NY,NX) = 1.0
  ENDDO

  !write(*,*) "(ATSEcoSIMAdvance) RainFalPrec_col: ", RainFalPrec_col(1,1), "m/s, PrecAsSnow: " , SnoFalPrec_col(1,1), " m/s" 
  PSIAtFldCapacity_col = pressure_at_field_capacity
  PSIAtWiltPoint_col = pressure_at_wilting_point

  call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)

  !Actually I update this just in the loop above??
  !call UpdateSoilMoistureFromATS(I,J,NHW,NHE,NVN,NVS)

  VHeatCapacity1_vr(0,1,1) = 0.0

  !These arrays sometimes have junk in them when initialized so we have to zero
  !them out before iteration
  do NY=1, NYS
    Qinfl2MicPM(NY,NX)=0.0
    Hinfl2SoilM(NY,NX)=0.0
    Qinfl2MicP(NY,NX)=0.0
    Hinfl2Soil(NY,NX)=0.0
  enddo

  !This does the subcycling of the land surface model
  DO M=1,NPH

    call RunSurfacePhysModelM(I,J,M,NHE,NHW,NVS,NVN,ResistanceLitRLay,&
      KSatReductByRainKineticEnergy,TopLayWatVol,HeatFluxAir2Soi,Qinfl2MicPM,Hinfl2SoilM)


    Qinfl2MicP = Qinfl2MicP+Qinfl2MicPM
    Hinfl2Soil = Hinfl2Soil+Hinfl2SoilM
    !also update state variables for iteration M 
    call UpdateSurfaceAtM(I,J,M,NHW,NHE,NVN,NVS)

  ENDDO
  do NY=1,NYS
    call SnowMassUpdate(I,J,NY,NX,Qinfl2MicPM(NY,NX),Hinfl2SoilM(NY,NX))
  ENDDO

  !write(*,*) "After SnowMassUpdate computation of Qinfl2MicP: "
  ! Check for NaN in surf_w_source
  if (any(is_nan(Qinfl2MicP(:,1)))) then
    write(*,*) "NaN found in Winfl2MicP at indices:", pack([(i, i=1,NYS)], is_nan(Qinfl2MicP(:,1)))
  end if

  DO NY=1,NYS
    !for every column send the top layer to the transfer var
    !Convert heat and water flux flux from the subcycled value
    !to ATS units (flux / s)
    write(*,*) "NY: ", NY, " Hinfl2Soil: ", Hinfl2Soil(NY,1), " Qinfl2MicP: ", Qinfl2MicP(NY,1)
    surf_e_source(NY) = Hinfl2Soil(NY,1) / (dts_HeatWatTP)
    surf_w_source(NY) = Qinfl2MicP(NY,1) / (dts_HeatWatTP)
    surf_snow_depth(NY) = SnowDepth_col(NY,1)
  ENDDO
 
  !write(*,*) "After setting surf_w_source to Qinfl2MicP: "
  ! Check for NaN in surf_w_source
  if (any(is_nan(surf_w_source))) then
    write(*,*) "NaN found in surf_w_source at indices:", pack([(i, i=1,NYS)], is_nan(surf_w_source))
  end if

  !Compute potential water loss(or gain) before next EcoSIM run
  !Wat_next = VLWatMicP1_vr(1,1,1) - Qinfl2MicP(NY,1) / (dts_HeatWatTP)
  !write(*,*) "(End EcoSIM Advance) Total Water Volume in top layer: ", VLWatMicP1_vr(1,1,1), " m, Q_w: ", surf_w_source(1)
  !write(*,*) "After an hour of this flux the water content should be: ", Wat_next

  end subroutine RunEcoSIMSurfaceBalance

end module ATSEcoSIMAdvanceMod
