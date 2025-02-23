module ClimForcDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod, only : warming_exp
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8) :: RMAX       !maximum hourly radiation,	[MJ m-2 h-1]
  real(r8) :: TAVG1      !parameter to calculate hourly  air temperature from daily value	[oC]
  real(r8) :: TAVG2      !parameter to calculate hourly  air temperature from daily value	[oC]
  real(r8) :: TAVG3      !parameter to calculate hourly  air temperature from daily value	[oC]
  real(r8) :: AMP1       !parameter to calculate hourly  air temperature from daily value	[oC]
  real(r8) :: AMP2       !parameter to calculate hourly  air temperature from daily value	[oC]
  real(r8) :: AMP3       !parameter to calculate hourly  air temperature from daily value	[oC]
  real(r8) :: VAVG1      !parameter to calculate hourly  vapor pressure from daily value	[kPa]
  real(r8) :: VAVG2      !parameter to calculate hourly  vapor pressure from daily value	[kPa]
  real(r8) :: VAVG3      !parameter to calculate hourly  vapor pressure from daily value	[kPa]
  real(r8) :: VMP1       !parameter to calculate hourly  vapor pressure from daily value	[kPa]
  real(r8) :: VMP2       !parameter to calculate hourly  vapor pressure from daily value	[kPa]
  real(r8) :: VMP3       !parameter to calculate hourly  vapor pressure from daily value	[kPa]
  real(r8) :: SAZI                              !solar azimuth of solar angle
  real(r8) :: SCOS                              !cosine of solar angle
  real(r8) :: DOY                               !day of year

  real(r8) :: atm_co2_mon(12)
  real(r8) :: atm_ch4_mon(12)
  real(r8) :: atm_n2o_mon(12)
  real(r8) :: TMPX(366)                         !maximum daily air temperature, [oC]
  real(r8) :: TMPN(366)                         !minimum daily air temperature, [oC]
  real(r8) :: SRAD(366)                         !daily solar radiation, [MJ m-2 d-1]
  real(r8) :: RAIN(366)                         !daily precipitation, [mm d-1 ]
  real(r8) :: WIND(366)                         !daily wind travel, [m d-1]
  real(r8) :: DWPT(2,366)                       !daily dewpoint temperature, [oC]

  real(r8) :: TMP_hrly(24,366)                      !hourly air temperature, [oC]
  real(r8) :: SWRad_hrly(24,366)                !hourly solar radiation, [MJ m-2 h-1]
  real(r8) :: RAINH(24,366)                     !hourly precipitation, [mm h-1]
  real(r8) :: WINDH(24,366)                     !hourly wind speed, [m h-1]
  real(r8) :: DWPTH(24,366)                     !hourly dewpoint temperature, [oC]
  real(r8) :: RadLWClm(24,366)                     !longwave radiation (MJ m-2 h-1)
  real(r8) :: PBOT_hrly(24,366)
  real(r8) :: DRAD(12)                          !change factor for radiation, [-]
  real(r8) :: DTMPX(12)                         !change factor for maximum temperature, [-]
  real(r8) :: DTMPN(12)                         !change factor for minimum temperature, [-]
  real(r8) :: DHUM(12)                          !change factor for humidity, [-]
  real(r8) :: DPREC(12)                         !change factor for precipitation, [-]
  real(r8) :: DWIND(12)                         !change factor for wind speed, [-]
  real(r8) :: DCN4R(12)                         !change factor for NH4 in precipitation, [-]
  real(r8) :: DCNOR(12)                         !change factor for NO3 in precipitation, [-]

  real(r8),target,allocatable ::  Eco_RadSW_col(:,:)       !shortwave radiation absorbed by the ecosystem [MJ/h]
  real(r8),target,allocatable ::  TKS_ref_vr(:,:,:,:)      !reference tempeature profile from control run to warming experiment [K]
  real(r8),target,allocatable ::  WDPTHD(:,:,:)                      !
  real(r8),target,allocatable ::  TDTPX(:,:,:)                       !accumulated change  for maximum temperature, [-]
  real(r8),target,allocatable ::  TDTPN(:,:,:)                       !accumulated change  for minimum temperature, [-]
  real(r8),target,allocatable ::  TDRAD(:,:,:)                       !accumulated change  for radiation, [-]
  real(r8),target,allocatable ::  TDHUM(:,:,:)                       !accumulated change  for humidity, [-]
  real(r8),target,allocatable ::  TDPRC(:,:,:)                       !accumulated change  for precipitation, [-]
  real(r8),target,allocatable ::  TDWND(:,:,:)                       !accumulated change  for wind speed, [-]
  real(r8),target,allocatable ::  TDCN4(:,:,:)                       !accumulated change  for NH4 in precipitation, [-]
  real(r8),target,allocatable ::  TDCNO(:,:,:)                       !accumulated change  for NO3 in precipitation, [-]
  real(r8),target,allocatable ::  TCA_col(:,:)                           !air temperature, [oC]
  real(r8),target,allocatable ::  TairK_col(:,:)                         !air temperature, [K]
  real(r8),target,allocatable ::  WindSpeedAtm_col(:,:)                  !wind speed, [m h-1]
  real(r8),target,allocatable ::  VPA_col(:,:)                           !atmospheric vapor concentration, [m3 m-3]
  real(r8),target,allocatable ::  VPK_col(:,:)                           !atmospheric vapor pressure, [kPa]
  real(r8),target,allocatable ::  PBOT_col(:,:)                          !atmospheric pressure [kPa]
  real(r8),target,allocatable ::  DayLensCurr_col(:,:)                          !daylength, [h]
  real(r8),target,allocatable ::  DayLenthPrev_col(:,:)                          !daylength of previous day, [h]
  real(r8),target,allocatable ::  DayLenthMax(:,:)                          !maximum daylength, [h]
  real(r8),target,allocatable ::  OMEGAG(:,:,:)                      !sine of solar beam on leaf surface, [-]
  real(r8),target,allocatable ::  LWRadSky_col(:,:)                  !sky longwave radiation , [MJ/h]
  real(r8),target,allocatable ::  TRAD(:,:)                          !total daily solar radiation, [MJ d-1]
  real(r8),target,allocatable ::  TAMX(:,:)                          !daily maximum air temperature , [oC]
  real(r8),target,allocatable ::  TAMN(:,:)                          !daily minimum air temperature , [oC]
  real(r8),target,allocatable ::  HUDX(:,:)                          !daily maximum vapor pressure , [kPa]
  real(r8),target,allocatable ::  HUDN(:,:)                          !daily minimum vapor pressure , [kPa]
  real(r8),target,allocatable ::  TWIND(:,:)                         !total daily wind travel, [m d-1]
!  real(r8),target,allocatable ::  PrecDaily_col(:,:)                !total daily precipitation, [m d-1]
  real(r8),target,allocatable ::  SkyLonwRad_col(:,:)                !sky longwave radiation , [MJ m-2 h-1]
  real(r8),target,allocatable ::  TempOffset_col(:,:)                !TempOffset_col for calculating temperature in Arrhenius curves, [oC]
  real(r8),target,allocatable ::  PrecDirect2Grnd_col(:,:)                     !direct precipitation at ground surface used to calculate soil erosion, [m h-1]
  real(r8),target,allocatable ::  PrecIndirect2Grnd_col(:,:)                         !indirect precipitation at ground surface used to calculate soil erosion, [m h-1]
  real(r8),target,allocatable ::  CO2EI(:,:)                         !initial atmospheric CO2 concentration, [umol mol-1]
  real(r8),target,allocatable ::  CCO2EI(:,:)                        !initial atmospheric CO2 concentration, [gC m-3]

  real(r8),target,allocatable ::  AtmGasCgperm3(:,:,:)               !atmospheric gas concentration in g m-3
  real(r8),target,allocatable ::  AtmGmms(:,:,:)                     !atmospheric gas concentration in umol mol-1
  real(r8),target,allocatable ::  OXYE_col(:,:)                          !atmospheric O2 concentration, [umol mol-1]
  real(r8),target,allocatable ::  Z2OE_col(:,:)                          !atmospheric N2O concentration, [umol mol-1]
  real(r8),target,allocatable ::  Z2GE_col(:,:)                      !atmospheric N2 concentration, [umol mol-1]
  real(r8),target,allocatable ::  ZNH3E_col(:,:)                     !atmospheric NH3 concentration, [umol mol-1]
  real(r8),target,allocatable ::  CH4E_col(:,:)                      !atmospheric CH4 concentration, [umol mol-1]
  real(r8),target,allocatable ::  H2GE_col(:,:)                      !atmospheric H2 concentration, [umol mol-1]
  real(r8),target,allocatable ::  CO2E_col(:,:)                      !atmospheric CO2 concentration, [umol mol-1]
  real(r8),target,allocatable ::  ARGE_col(:,:)                      !atmospheric AR concentration, [umol mol-1]
  real(r8),target,allocatable ::  SolarNoonHour_col(:,:)             !time of solar noon, [h]
  real(r8),target,allocatable ::  RadSWDirect_col(:,:)               !direct shortwave radiation, [W m-2]
  real(r8),target,allocatable ::  RadSWDiffus_col(:,:)               !diffuse shortwave radiation, [W m-2]
  real(r8),target,allocatable ::  RadPARDirect_col(:,:)              !direct PAR, [umol m-2 s-1]
  real(r8),target,allocatable ::  RadPARDiffus_col(:,:)              !diffuse PAR, [umol m-2 s-1]
  real(r8),target,allocatable ::  SineSunInclAngle_col(:,:)          !sine of solar angle, [-]
  real(r8),target,allocatable ::  SineSunInclAnglNxtHour_col(:,:)    !sine of solar angle next hour, [-]
  real(r8),target,allocatable ::  TLEX_col(:,:)                          !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8),target,allocatable ::  TSHX_col(:,:)                          !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8),target,allocatable ::  Air_Heat_Latent_store_col(:,:)                          !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8),target,allocatable ::  Air_Heat_Sens_store_col(:,:)                          !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8),target,allocatable ::  SoilHeatSrcDepth_col(:,:)                        !depth of soil heat sink/source, [m]
  real(r8),target,allocatable ::  TKSD(:,:)                          !temperature of soil heat sink/source, [oC]
  real(r8),target,allocatable ::  ATCAI(:,:)                         !initial mean annual air temperature, [oC]
  real(r8),target,allocatable ::  RadSWSolarBeam_col(:,:)                           !shortwave radiation in solar beam, [MJ m-2 h-1]
  real(r8),target,allocatable ::  RadPARSolarBeam_col(:,:)                           !PAR radiation in solar beam, [umol m-2 s-1]
  real(r8),target,allocatable ::  ATCA(:,:)                          !mean annual air temperature, [oC]
  real(r8),target,allocatable ::  ATCS(:,:)                          !mean annual soil temperature, [oC]
  real(r8),target,allocatable ::  TairKClimMean(:,:)                 !mean annual air temperature, [K]
  real(r8),target,allocatable ::  ATKS(:,:)                          !mean annual soil temperature, [K]
  real(r8),target,allocatable ::  RainFalPrec_col(:,:)                   !rainfall, [m3 d-2 h-1]
  real(r8),target,allocatable ::  SnoFalPrec_col(:,:)                    !snowfall, [m3 d-2 h-1]
  real(r8),target,allocatable ::  Irrigation_col(:,:)                 !irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  PrecAtm_col(:,:)                         !rainfall + snowfall, [m3 d-2 h-1]
  real(r8),target,allocatable ::  PrecRainAndIrrig_col(:,:)                         !rainfall + irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  EnergyImpact4Erosion(:,:)                         !cumulative rainfall energy impact on soil surface
  real(r8),target,allocatable ::  pH_rain_col(:,:)                   !precipitation pH, [-]
  real(r8),target,allocatable ::  CN4RI(:,:)                         !precipitation initial NH4 concentration, [mol N m-3]
  real(r8),target,allocatable ::  CNORI(:,:)                         !precipitation initial NO3 concentration, [mol N m-3]
  real(r8),target,allocatable ::  NH4_rain_mole_conc(:,:)            !precipitation  NH4 concentration, [mol m-3]
  real(r8),target,allocatable ::  NO3_rain_mole_conc(:,:)            !precipitation  NO3 concentration, [mol m-3]
  real(r8),target,allocatable ::  H2PO4_rain_mole_conc(:,:)          !precipitation  H2PO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CALR(:,:)                          !precipitation  Al concentration, [mol m-3]
  real(r8),target,allocatable ::  CFER(:,:)                          !precipitation  Fe concentration, [mol m-3]
  real(r8),target,allocatable ::  CHYR(:,:)                          !precipitation  H concentration, [mol m-3]
  real(r8),target,allocatable ::  CCAR(:,:)                          !precipitation  Ca concentration, [mol m-3]
  real(r8),target,allocatable ::  CMGR(:,:)                          !precipitation  Mg concentration, [mol m-3]
  real(r8),target,allocatable ::  CNAR(:,:)                          !precipitation  Na concentration, [mol m-3]
  real(r8),target,allocatable ::  CKAR(:,:)                          !precipitation  K concentration, [mol m-3]
  real(r8),target,allocatable ::  COHR(:,:)                          !precipitation  OH concentration, [mol m-3]
  real(r8),target,allocatable ::  CSOR(:,:)                          !precipitation  SO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CCLR(:,:)                          !precipitation  Cl concentration, [mol m-3]
  real(r8),target,allocatable ::  CC3R(:,:)                          !precipitation  CO3 concentration, [mol m-3]
  real(r8),target,allocatable ::  CHCR(:,:)                          !precipitation  HCO3 concentration, [mol m-3]
  real(r8),target,allocatable ::  trcg_rain_mole_conc_col(:,:,:)     !precipitation volatile concentration, [mol m-3]
  real(r8),target,allocatable ::  CAL1R(:,:)                         !precipitation  AlOH concentration, [mol m-3]
  real(r8),target,allocatable ::  CAL2R(:,:)                         !precipitation  AlOH2 concentration, [mol m-3]
  real(r8),target,allocatable ::  CAL3R(:,:)                         !precipitation  AlOH3 concentration, [mol m-3]
  real(r8),target,allocatable ::  CAL4R(:,:)                         !precipitation  AlOH4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CALSR(:,:)                         !precipitation  AlSO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CFE1R(:,:)                         !precipitation  FeOH concentration, [mol m-3]
  real(r8),target,allocatable ::  CFE2R(:,:)                         !precipitation  FeOH2 concentration, [mol m-3]
  real(r8),target,allocatable ::  CFE3R(:,:)                         !precipitation  FeOH3 concentration, [mol m-3]
  real(r8),target,allocatable ::  CFE4R(:,:)                         !precipitation  FeOH4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CFESR(:,:)                         !precipitation  FeSO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CCAOR(:,:)                         !precipitation  CaOH concentration, [mol m-3]
  real(r8),target,allocatable ::  CCACR(:,:)                         !precipitation  CaCO3 concentration, [mol m-3]
  real(r8),target,allocatable ::  CCAHR(:,:)                         !precipitation  CaHCO3 concentration, [mol m-3]
  real(r8),target,allocatable ::  CCASR(:,:)                         !precipitation  CaSO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CMGOR(:,:)                         !precipitation  MgOH concentration, [mol m-3]
  real(r8),target,allocatable ::  CMGCR(:,:)                         !precipitation  MgCO3 concentration, [mol m-3]
  real(r8),target,allocatable ::  CMGHR(:,:)                         !precipitation  MgHCO3 concentration, [mol m-3]
  real(r8),target,allocatable ::  CMGSR(:,:)                         !precipitation  MgSO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CNACR(:,:)                         !precipitation  NaCO3 concentration, [mol m-3]
  real(r8),target,allocatable ::  CNASR(:,:)                         !precipitation  NaSO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CKASR(:,:)                         !precipitation  K concentration, [mol m-3]
  real(r8),target,allocatable ::  CH0PR(:,:)                         !precipitation  PO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  HPO4_rain_mole_conc(:,:)                !precipitation  HPO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CH3PR(:,:)                         !precipitation  H3PO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CF1PR(:,:)                         !precipitation  FeHPO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CF2PR(:,:)                         !precipitation  FeH2PO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CC0PR(:,:)                         !precipitation  CaPO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CC1PR(:,:)                         !precipitation  CaHPO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  CC2PR(:,:)                         !precipitation  CaH4P2O8 concentration, [mol m-3]
  real(r8),target,allocatable ::  CM1PR(:,:)                         !precipitation  MgHPO4 concentration, [mol m-3]
  real(r8),target,allocatable ::  GDD_col(:,:)    !growing degree day with base temperature at oC
  real(r8),target,allocatable ::  PrecHeat_col(:,:)    !precipitation heat to surface [MJ/d2/h]
  real(r8),target,allocatable ::  RainLitr_col(:,:)  !water from aboveground falling litter
  contains
!----------------------------------------------------------------------

  subroutine InitClimForcData()
  use TracerIDMod
  implicit none

  if(len(trim(warming_exp))>10)then
    allocate(TKS_ref_vr(8784,JZ,JY,JX));TKS_ref_vr=0._r8
  endif
  allocate(Eco_RadSW_col(JY,JX)); Eco_RadSW_col=0._r8
  allocate(GDD_col(JY,JX)); GDD_col=0._r8
  allocate(WDPTHD(366,JY,JX));  WDPTHD=0._r8
  allocate(TDTPX(12,JY,JX));    TDTPX=0._r8
  allocate(TDTPN(12,JY,JX));    TDTPN=0._r8
  allocate(TDRAD(12,JY,JX));    TDRAD=0._r8
  allocate(TDHUM(12,JY,JX));    TDHUM=0._r8
  allocate(TDPRC(12,JY,JX));    TDPRC=0._r8
  allocate(TDWND(12,JY,JX));    TDWND=0._r8
  allocate(TDCN4(12,JY,JX));    TDCN4=0._r8
  allocate(TDCNO(12,JY,JX));    TDCNO=0._r8
  allocate(PrecHeat_col(JY,JX)); PrecHeat_col=0._r8
  allocate(RainLitr_col(JY,JX)); RainLitr_col=0._r8
  allocate(TCA_col(JY,JX));         TCA_col=0._r8
  allocate(TairK_col(JY,JX));         TairK_col=0._r8
  allocate(WindSpeedAtm_col(JY,JX));          WindSpeedAtm_col=0._r8
  allocate(VPA_col(JY,JX));         VPA_col=0._r8
  allocate(VPK_col(JY,JX));         VPK_col=0._r8
  allocate(PBOT_col(JY,JX));        PBOT_col=1.01325E+02_r8
  allocate(DayLensCurr_col(JY,JX));        DayLensCurr_col=0._r8
  allocate(DayLenthPrev_col(JY,JX));        DayLenthPrev_col=0._r8
  allocate(DayLenthMax(JY,JX));        DayLenthMax=0._r8
  allocate(OMEGAG(NumOfSkyAzimuthSects,JY,JX));  OMEGAG=0._r8
  allocate(LWRadSky_col(JY,JX));         LWRadSky_col=0._r8
  allocate(TRAD(JY,JX));        TRAD=0._r8
  allocate(TAMX(JY,JX));        TAMX=0._r8
  allocate(TAMN(JY,JX));        TAMN=0._r8
  allocate(HUDX(JY,JX));        HUDX=0._r8
  allocate(HUDN(JY,JX));        HUDN=0._r8
  allocate(TWIND(JY,JX));       TWIND=0._r8
!  allocate(PrecDaily_col(JY,JX));        PrecDaily_col=0._r8
  allocate(SkyLonwRad_col(JY,JX));        SkyLonwRad_col=0._r8
  allocate(TempOffset_col(JY,JX));      TempOffset_col=0._r8
  allocate(PrecDirect2Grnd_col(JY,JX));       PrecDirect2Grnd_col=0._r8
  allocate(PrecIndirect2Grnd_col(JY,JX));       PrecIndirect2Grnd_col=0._r8
  allocate(CO2EI(JY,JX));       CO2EI=0._r8
  allocate(CCO2EI(JY,JX));      CCO2EI=0._r8

  allocate(AtmGasCgperm3(idg_beg:idg_end,JY,JX)); AtmGasCgperm3=0._r8
  allocate(AtmGmms(idg_beg:idg_end,JY,JX)); AtmGmms=0._r8

  allocate(OXYE_col(JY,JX));        OXYE_col=0._r8
  allocate(Z2OE_col(JY,JX));        Z2OE_col=0._r8
  allocate(Z2GE_col(JY,JX));        Z2GE_col=0._r8
  allocate(ZNH3E_col(JY,JX));       ZNH3E_col=0._r8
  allocate(CH4E_col(JY,JX));        CH4E_col=0._r8
  allocate(H2GE_col(JY,JX));        H2GE_col=0._r8
  allocate(CO2E_col(JY,JX));        CO2E_col=0._r8
  allocate(ARGE_col(JY,JX));        ARGE_col=0._r8
  allocate(SolarNoonHour_col(JY,JX));       SolarNoonHour_col=0._r8
  allocate(RadSWDirect_col(JY,JX));        RadSWDirect_col=0._r8
  allocate(RadSWDiffus_col(JY,JX));        RadSWDiffus_col=0._r8
  allocate(RadPARDirect_col(JY,JX));        RadPARDirect_col=0._r8
  allocate(RadPARDiffus_col(JY,JX));        RadPARDiffus_col=0._r8
  allocate(SineSunInclAngle_col(JY,JX));        SineSunInclAngle_col=0._r8
  allocate(SineSunInclAnglNxtHour_col(JY,JX));       SineSunInclAnglNxtHour_col=0._r8
  allocate(TLEX_col(JY,JX));        TLEX_col=0._r8
  allocate(TSHX_col(JY,JX));        TSHX_col=0._r8
  allocate(Air_Heat_Latent_store_col(JY,JX));        Air_Heat_Latent_store_col=0._r8
  allocate(Air_Heat_Sens_store_col(JY,JX));        Air_Heat_Sens_store_col=0._r8
  allocate(SoilHeatSrcDepth_col(JY,JX));      SoilHeatSrcDepth_col=0._r8
  allocate(TKSD(JY,JX));        TKSD=0._r8
  allocate(ATCAI(JY,JX));       ATCAI=0._r8
  allocate(RadSWSolarBeam_col(JY,JX));         RadSWSolarBeam_col=0._r8
  allocate(RadPARSolarBeam_col(JY,JX));         RadPARSolarBeam_col=0._r8
  allocate(ATCA(JY,JX));        ATCA=0._r8
  allocate(ATCS(JY,JX));        ATCS=0._r8
  allocate(TairKClimMean(JY,JX));        TairKClimMean=0._r8
  allocate(ATKS(JY,JX));        ATKS=0._r8
  allocate(RainFalPrec_col(JY,JX));       RainFalPrec_col=0._r8
  allocate(Irrigation_col(JY,JX));    Irrigation_col=0._r8
  allocate(SnoFalPrec_col(JY,JX));       SnoFalPrec_col=0._r8
  allocate(PrecAtm_col(JY,JX));       PrecAtm_col=0._r8
  allocate(PrecRainAndIrrig_col(JY,JX));       PrecRainAndIrrig_col=0._r8
  allocate(EnergyImpact4Erosion(JY,JX));       EnergyImpact4Erosion=0._r8
  allocate(pH_rain_col(JY,JX));         pH_rain_col=0._r8
  allocate(CN4RI(JY,JX));       CN4RI=0._r8
  allocate(CNORI(JY,JX));       CNORI=0._r8
  allocate(NH4_rain_mole_conc(JY,JX));        NH4_rain_mole_conc=0._r8
  allocate(NO3_rain_mole_conc(JY,JX));        NO3_rain_mole_conc=0._r8
  allocate(H2PO4_rain_mole_conc(JY,JX));        H2PO4_rain_mole_conc=0._r8
  allocate(CALR(JY,JX));        CALR=0._r8
  allocate(CFER(JY,JX));        CFER=0._r8
  allocate(CHYR(JY,JX));        CHYR=0._r8
  allocate(CCAR(JY,JX));        CCAR=0._r8
  allocate(CMGR(JY,JX));        CMGR=0._r8
  allocate(CNAR(JY,JX));        CNAR=0._r8
  allocate(CKAR(JY,JX));        CKAR=0._r8
  allocate(COHR(JY,JX));        COHR=0._r8
  allocate(CSOR(JY,JX));        CSOR=0._r8
  allocate(CCLR(JY,JX));        CCLR=0._r8
  allocate(CC3R(JY,JX));        CC3R=0._r8
  allocate(CHCR(JY,JX));        CHCR=0._r8
  allocate(trcg_rain_mole_conc_col(idg_beg:idg_NH3,JY,JX));        trcg_rain_mole_conc_col=0._r8
  allocate(CAL1R(JY,JX));       CAL1R=0._r8
  allocate(CAL2R(JY,JX));       CAL2R=0._r8
  allocate(CAL3R(JY,JX));       CAL3R=0._r8
  allocate(CAL4R(JY,JX));       CAL4R=0._r8
  allocate(CALSR(JY,JX));       CALSR=0._r8
  allocate(CFE1R(JY,JX));       CFE1R=0._r8
  allocate(CFE2R(JY,JX));       CFE2R=0._r8
  allocate(CFE3R(JY,JX));       CFE3R=0._r8
  allocate(CFE4R(JY,JX));       CFE4R=0._r8
  allocate(CFESR(JY,JX));       CFESR=0._r8
  allocate(CCAOR(JY,JX));       CCAOR=0._r8
  allocate(CCACR(JY,JX));       CCACR=0._r8
  allocate(CCAHR(JY,JX));       CCAHR=0._r8
  allocate(CCASR(JY,JX));       CCASR=0._r8
  allocate(CMGOR(JY,JX));       CMGOR=0._r8
  allocate(CMGCR(JY,JX));       CMGCR=0._r8
  allocate(CMGHR(JY,JX));       CMGHR=0._r8
  allocate(CMGSR(JY,JX));       CMGSR=0._r8
  allocate(CNACR(JY,JX));       CNACR=0._r8
  allocate(CNASR(JY,JX));       CNASR=0._r8
  allocate(CKASR(JY,JX));       CKASR=0._r8
  allocate(CH0PR(JY,JX));       CH0PR=0._r8
  allocate(HPO4_rain_mole_conc(JY,JX));       HPO4_rain_mole_conc=0._r8
  allocate(CH3PR(JY,JX));       CH3PR=0._r8
  allocate(CF1PR(JY,JX));       CF1PR=0._r8
  allocate(CF2PR(JY,JX));       CF2PR=0._r8
  allocate(CC0PR(JY,JX));       CC0PR=0._r8
  allocate(CC1PR(JY,JX));       CC1PR=0._r8
  allocate(CC2PR(JY,JX));       CC2PR=0._r8
  allocate(CM1PR(JY,JX));       CM1PR=0._r8

  end subroutine InitClimForcData

!----------------------------------------------------------------------
  subroutine DestructClimForcData
  use abortutils, only : destroy
  implicit none
  call destroy(PrecHeat_col)
  call destroy(RainLitr_col)
  call destroy(WDPTHD)
  call destroy(TDTPX)
  call destroy(TDTPN)
  call destroy(TDRAD)
  call destroy(TDHUM)
  call destroy(TDPRC)
  call destroy(TDWND)
  call destroy(TDCN4)
  call destroy(TDCNO)
  call destroy(TCA_col)
  call destroy(TairK_col)
  call destroy(WindSpeedAtm_col)
  call destroy(VPA_col)
  call destroy(VPK_col)
  call destroy(PBOT_col)
  call destroy(DayLensCurr_col)
  call destroy(DayLenthPrev_col)
  call destroy(DayLenthMax)
  call destroy(OMEGAG)
  call destroy(LWRadSky_col)
  call destroy(TRAD)
  call destroy(TAMX)
  call destroy(TAMN)
  call destroy(HUDX)
  call destroy(HUDN)
  call destroy(TWIND)
  call destroy(ARGE_col)
!  call destroy(PrecDaily_col)
  call destroy(SkyLonwRad_col)
  call destroy(TempOffset_col)
  call destroy(PrecDirect2Grnd_col)
  call destroy(PrecIndirect2Grnd_col)
  call destroy(CO2EI)
  call destroy(CCO2EI)

  call destroy(AtmGasCgperm3)
  call destroy(AtmGmms)
  call destroy(TKS_ref_vr)
  call destroy(Eco_RadSW_col)
  call destroy(OXYE_col)
  call destroy(Z2OE_col)
  call destroy(Z2GE_col)
  call destroy(ZNH3E_col)
  call destroy(CH4E_col)
  call destroy(H2GE_col)

  call destroy(SolarNoonHour_col)
  call destroy(CO2E_col)
  call destroy(RadSWDirect_col)
  call destroy(RadSWDiffus_col)
  call destroy(RadPARDirect_col)
  call destroy(RadPARDiffus_col)
  call destroy(SineSunInclAngle_col)
  call destroy(SineSunInclAnglNxtHour_col)
  call destroy(TLEX_col)
  call destroy(TSHX_col)
  call destroy(Air_Heat_Latent_store_col)
  call destroy(Air_Heat_Sens_store_col)
  call destroy(SoilHeatSrcDepth_col)
  call destroy(TKSD)
  call destroy(ATCAI)
  call destroy(RadSWSolarBeam_col)
  call destroy(RadPARSolarBeam_col)
  call destroy(ATCA)
  call destroy(ATCS)
  call destroy(TairKClimMean)
  call destroy(ATKS)
  call destroy(RainFalPrec_col)
  call destroy(Irrigation_col)
  call destroy(SnoFalPrec_col)
  call destroy(PrecAtm_col)
  call destroy(PrecRainAndIrrig_col)
  call destroy(EnergyImpact4Erosion)
  call destroy(pH_rain_col)
  call destroy(CN4RI)
  call destroy(CNORI)
  call destroy(NH4_rain_mole_conc)
  call destroy(NO3_rain_mole_conc)
  call destroy(H2PO4_rain_mole_conc)
  call destroy(CALR)
  call destroy(CFER)
  call destroy(CHYR)
  call destroy(CCAR)
  call destroy(CMGR)
  call destroy(CNAR)
  call destroy(CKAR)
  call destroy(COHR)
  call destroy(CSOR)
  call destroy(CCLR)
  call destroy(CC3R)
  call destroy(CHCR)
  call destroy(trcg_rain_mole_conc_col)
  call destroy(CAL1R)
  call destroy(CAL2R)
  call destroy(CAL3R)
  call destroy(CAL4R)
  call destroy(CALSR)
  call destroy(CFE1R)
  call destroy(CFE2R)
  call destroy(CFE3R)
  call destroy(CFE4R)
  call destroy(CFESR)
  call destroy(CCAOR)
  call destroy(CCACR)
  call destroy(CCAHR)
  call destroy(CCASR)
  call destroy(CMGOR)
  call destroy(CMGCR)
  call destroy(CMGHR)
  call destroy(CMGSR)
  call destroy(CNACR)
  call destroy(CNASR)
  call destroy(CKASR)
  call destroy(CH0PR)
  call destroy(HPO4_rain_mole_conc)
  call destroy(CH3PR)
  call destroy(CF1PR)
  call destroy(CF2PR)
  call destroy(CC0PR)
  call destroy(CC1PR)
  call destroy(CC2PR)
  call destroy(CM1PR)
  end subroutine DestructClimForcData
end module ClimForcDataType
