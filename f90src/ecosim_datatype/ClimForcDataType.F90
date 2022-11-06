module ClimForcDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: DECLIN
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

  real(r8) :: TMPX(366)                         !maximum daily air temperature, [oC]
  real(r8) :: TMPN(366)                         !minimum daily air temperature, [oC]
  real(r8) :: SRAD(366)                         !daily solar radiation, [MJ m-2 d-1]
  real(r8) :: RAIN(366)                         !daily precipitation, [mm d-1 ]
  real(r8) :: WIND(366)                         !daily wind travel, [m d-1]
  real(r8) :: DWPT(2,366)                       !daily dewpoint temperature, [oC]
  real(r8) :: TMPH(24,366)                      !hourly air temperature, [oC]
  real(r8) :: SRADH(24,366)                     !hourly solar radiation, [MJ m-2 h-1]
  real(r8) :: RAINH(24,366)                     !hourly precipitation, [mm h-1]
  real(r8) :: WINDH(24,366)                     !hourly wind speed, [m h-1]
  real(r8) :: DWPTH(24,366)                     !hourly dewpoint temperature, [oC]
  real(r8) :: XRADH(24,366)                     !longwave radiation (MJ m-2 h-1)

  real(r8) :: DRAD(12)                          !change factor for radiation, [-]
  real(r8) :: DTMPX(12)                         !change factor for maximum temperature, [-]
  real(r8) :: DTMPN(12)                         !change factor for minimum temperature, [-]
  real(r8) :: DHUM(12)                          !change factor for humidity, [-]
  real(r8) :: DPREC(12)                         !change factor for precipitation, [-]
  real(r8) :: DWIND(12)                         !change factor for wind speed, [-]
  real(r8) :: DCO2E(12)                         !change factor for atmospheric CO2 concentration, [-]
  real(r8) :: DCN4R(12)                         !change factor for NH4 in precipitation, [-]
  real(r8) :: DCNOR(12)                         !change factor for NO3 in precipitation, [-]

  real(r8),pointer ::  WDPTHD(:,:,:)                      !
  real(r8),pointer ::  TDTPX(:,:,:)                       !accumulated change  for maximum temperature, [-]
  real(r8),pointer ::  TDTPN(:,:,:)                       !accumulated change  for minimum temperature, [-]
  real(r8),pointer ::  TDRAD(:,:,:)                       !accumulated change  for radiation, [-]
  real(r8),pointer ::  TDHUM(:,:,:)                       !accumulated change  for humidity, [-]
  real(r8),pointer ::  TDPRC(:,:,:)                       !accumulated change  for precipitation, [-]
  real(r8),pointer ::  TDWND(:,:,:)                       !accumulated change  for wind speed, [-]
  real(r8),pointer ::  TDCO2(:,:,:)                       !accumulated change  for atmospheric CO2 concentration, [-]
  real(r8),pointer ::  TDCN4(:,:,:)                       !accumulated change  for NH4 in precipitation, [-]
  real(r8),pointer ::  TDCNO(:,:,:)                       !accumulated change  for NO3 in precipitation, [-]
  real(r8),pointer ::  TCA(:,:)                           !air temperature, [oC]
  real(r8),pointer ::  TKA(:,:)                           !air temperature, [K]
  real(r8),pointer ::  UA(:,:)                            !wind speed, [m h-1]
  real(r8),pointer ::  VPA(:,:)                           !vapor concentration, [m3 m-3]
  real(r8),pointer ::  VPK(:,:)                           !vapor pressure, [kPa]
  real(r8),pointer ::  DYLN(:,:)                          !daylength, [h]
  real(r8),pointer ::  DYLX(:,:)                          !daylength of previous day, [h]
  real(r8),pointer ::  DYLM(:,:)                          !maximum daylength, [h]
  real(r8),pointer ::  OMEGAG(:,:,:)                      !sine of solar beam on leaf surface, [-]
  real(r8),pointer ::  THS(:,:)                           !sky longwave radiation , [MJ d-2 h-1]
  real(r8),pointer ::  TRAD(:,:)                          !total daily solar radiation, [MJ d-1]
  real(r8),pointer ::  TAMX(:,:)                          !daily maximum air temperature , [oC]
  real(r8),pointer ::  TAMN(:,:)                          !daily minimum air temperature , [oC]
  real(r8),pointer ::  HUDX(:,:)                          !daily maximum vapor pressure , [kPa]
  real(r8),pointer ::  HUDN(:,:)                          !daily minimum vapor pressure , [kPa]
  real(r8),pointer ::  TWIND(:,:)                         !total daily wind travel, [m d-1]
  real(r8),pointer ::  TRAI(:,:)                          !total daily precipitation, [m d-1]
  real(r8),pointer ::  THSX(:,:)                          !sky longwave radiation , [MJ m-2 h-1]
  real(r8),pointer ::  OFFSET(:,:)                        !offset for calculating temperature in Arrhenius curves, [oC]
  real(r8),pointer ::  PRECD(:,:)                         !precipitation at ground surface used to calculate soil erosion, [m h-1]
  real(r8),pointer ::  PRECB(:,:)                         !precipitation at ground surface used to calculate soil erosion, [m h-1]
  real(r8),pointer ::  CO2EI(:,:)                         !initial atmospheric CO2 concentration, [umol mol-1]
  real(r8),pointer ::  CCO2EI(:,:)                        !initial atmospheric CO2 concentration, [g m-3]

  real(r8),pointer ::  AtmGgms(:,:,:)                     !atmospheric gas concentration in g m-3
  real(r8),pointer ::  AtmGmms(:,:,:)                     !atmospheric gas concentration in umol mol-1
  real(r8),pointer ::  OXYE(:,:)                          !atmospheric O2 concentration, [umol mol-1]
  real(r8),pointer ::  Z2OE(:,:)                          !atmospheric N2O concentration, [umol mol-1]
  real(r8),pointer ::  Z2GE(:,:)                          !atmospheric N2 concentration, [umol mol-1]
  real(r8),pointer ::  ZNH3E(:,:)                         !atmospheric NH3 concentration, [umol mol-1]
  real(r8),pointer ::  CH4E(:,:)                          !atmospheric CH4 concentration, [umol mol-1]
  real(r8),pointer ::  H2GE(:,:)                          !atmospheric H2 concentration, [umol mol-1]
  real(r8),pointer ::  CO2E(:,:)                          !atmospheric CO2 concentration, [umol mol-1]

  real(r8),pointer ::  ZNOON(:,:)                         !time of solar noon, [h]
  real(r8),pointer ::  RADS(:,:)                          !direct shortwave radiation, [W m-2]
  real(r8),pointer ::  RADY(:,:)                          !diffuse shortwave radiation, [W m-2]
  real(r8),pointer ::  RAPS(:,:)                          !direct PAR, [umol m-2 s-1]
  real(r8),pointer ::  RAPY(:,:)                          !diffuse PAR, [umol m-2 s-1]
  real(r8),pointer ::  SSIN(:,:)                          !sine of solar angle, [-]
  real(r8),pointer ::  SSINN(:,:)                         !sine of solar angle next hour, [-]
  real(r8),pointer ::  TLEX(:,:)                          !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8),pointer ::  TSHX(:,:)                          !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8),pointer ::  TLEC(:,:)                          !total latent heat flux x boundary layer resistance, [MJ m-1]
  real(r8),pointer ::  TSHC(:,:)                          !total sensible heat flux x boundary layer resistance, [MJ m-1]
  real(r8),pointer ::  DPTHSK(:,:)                        !depth of soil heat sink/source, [m]
  real(r8),pointer ::  TKSD(:,:)                          !temperature of soil heat sink/source, [oC]
  real(r8),pointer ::  ATCAI(:,:)                         !initial mean annual air temperature, [oC]
  real(r8),pointer ::  RAD(:,:)                           !shortwave radiation in solar beam, [MJ m-2 h-1]
  real(r8),pointer ::  RAP(:,:)                           !PAR radiation in solar beam, [umol m-2 s-1]
  real(r8),pointer ::  ATCA(:,:)                          !mean annual air temperature, [oC]
  real(r8),pointer ::  ATCS(:,:)                          !mean annual soil temperature, [oC]
  real(r8),pointer ::  ATKA(:,:)                          !mean annual air temperature, [K]
  real(r8),pointer ::  ATKS(:,:)                          !mean annual soil temperature, [K]
  real(r8),pointer ::  PRECR(:,:)                         !rainfall, [m3 d-2 h-1]
  real(r8),pointer ::  PRECW(:,:)                         !snowfall, [m3 d-2 h-1]
  real(r8),pointer ::  PRECQ(:,:)                         !rainfall + snowfall, [m3 d-2 h-1]
  real(r8),pointer ::  PRECA(:,:)                         !rainfall + irrigation, [m3 d-2 h-1]
  real(r8),pointer ::  ENGYP(:,:)                         !cumulative rainfall energy impact on soil surface
  real(r8),pointer ::  PHR(:,:)                           !precipitation pH, [-]
  real(r8),pointer ::  CN4RI(:,:)                         !precipitation initial NH4 concentration, [g m-3]
  real(r8),pointer ::  CNORI(:,:)                         !precipitation initial NO3 concentration, [g m-3]
  real(r8),pointer ::  CN4R(:,:)                          !precipitation  NH4 concentration, [g m-3]
  real(r8),pointer ::  CN3R(:,:)                          !precipitation  NH3 concentration, [g m-3]
  real(r8),pointer ::  CNOR(:,:)                          !precipitation  NO3 concentration, [g m-3]
  real(r8),pointer ::  CPOR(:,:)                          !precipitation  H2PO4 concentration, [g m-3]
  real(r8),pointer ::  CALR(:,:)                          !precipitation  Al concentration, [g m-3]
  real(r8),pointer ::  CFER(:,:)                          !precipitation  Fe concentration, [g m-3]
  real(r8),pointer ::  CHYR(:,:)                          !precipitation  H concentration, [g m-3]
  real(r8),pointer ::  CCAR(:,:)                          !precipitation  Ca concentration, [g m-3]
  real(r8),pointer ::  CMGR(:,:)                          !precipitation  Mg concentration, [g m-3]
  real(r8),pointer ::  CNAR(:,:)                          !precipitation  Na concentration, [g m-3]
  real(r8),pointer ::  CKAR(:,:)                          !precipitation  K concentration, [g m-3]
  real(r8),pointer ::  COHR(:,:)                          !precipitation  OH concentration, [g m-3]
  real(r8),pointer ::  CSOR(:,:)                          !precipitation  SO4 concentration, [g m-3]
  real(r8),pointer ::  CCLR(:,:)                          !precipitation  Cl concentration, [g m-3]
  real(r8),pointer ::  CC3R(:,:)                          !precipitation  CO3 concentration, [g m-3]
  real(r8),pointer ::  CHCR(:,:)                          !precipitation  HCO3 concentration, [g m-3]
  real(r8),pointer ::  CCHR(:,:)                          !precipitation  CH4 concentration, [g m-3]
  real(r8),pointer ::  CAL1R(:,:)                         !precipitation  AlOH concentration, [g m-3]
  real(r8),pointer ::  CAL2R(:,:)                         !precipitation  AlOH2 concentration, [g m-3]
  real(r8),pointer ::  CAL3R(:,:)                         !precipitation  AlOH3 concentration, [g m-3]
  real(r8),pointer ::  CAL4R(:,:)                         !precipitation  AlOH4 concentration, [g m-3]
  real(r8),pointer ::  CALSR(:,:)                         !precipitation  AlSO4 concentration, [g m-3]
  real(r8),pointer ::  CFE1R(:,:)                         !precipitation  FeOH concentration, [g m-3]
  real(r8),pointer ::  CFE2R(:,:)                         !precipitation  FeOH2 concentration, [g m-3]
  real(r8),pointer ::  CFE3R(:,:)                         !precipitation  FeOH3 concentration, [g m-3]
  real(r8),pointer ::  CFE4R(:,:)                         !precipitation  FeOH4 concentration, [g m-3]
  real(r8),pointer ::  CFESR(:,:)                         !precipitation  FeSO4 concentration, [g m-3]
  real(r8),pointer ::  CCAOR(:,:)                         !precipitation  CaOH concentration, [g m-3]
  real(r8),pointer ::  CCACR(:,:)                         !precipitation  CaCO3 concentration, [g m-3]
  real(r8),pointer ::  CCAHR(:,:)                         !precipitation  CaHCO3 concentration, [g m-3]
  real(r8),pointer ::  CCASR(:,:)                         !precipitation  CaSO4 concentration, [g m-3]
  real(r8),pointer ::  CMGOR(:,:)                         !precipitation  MgOH concentration, [g m-3]
  real(r8),pointer ::  CMGCR(:,:)                         !precipitation  MgCO3 concentration, [g m-3]
  real(r8),pointer ::  CMGHR(:,:)                         !precipitation  MgHCO3 concentration, [g m-3]
  real(r8),pointer ::  CMGSR(:,:)                         !precipitation  MgSO4 concentration, [g m-3]
  real(r8),pointer ::  CNACR(:,:)                         !precipitation  NaCO3 concentration, [g m-3]
  real(r8),pointer ::  CNASR(:,:)                         !precipitation  NaSO4 concentration, [g m-3]
  real(r8),pointer ::  CKASR(:,:)                         !precipitation  K concentration, [g m-3]
  real(r8),pointer ::  CH0PR(:,:)                         !precipitation  PO4 concentration, [g m-3]
  real(r8),pointer ::  CH1PR(:,:)                         !precipitation  HPO4 concentration, [g m-3]
  real(r8),pointer ::  CH3PR(:,:)                         !precipitation  H3PO4 concentration, [g m-3]
  real(r8),pointer ::  CF1PR(:,:)                         !precipitation  FeHPO4 concentration, [g m-3]
  real(r8),pointer ::  CF2PR(:,:)                         !precipitation  FeH2PO4 concentration, [g m-3]
  real(r8),pointer ::  CC0PR(:,:)                         !precipitation  CaPO4 concentration, [g m-3]
  real(r8),pointer ::  CC1PR(:,:)                         !precipitation  CaHPO4 concentration, [g m-3]
  real(r8),pointer ::  CC2PR(:,:)                         !precipitation  CaH2PO4 concentration, [g m-3]
  real(r8),pointer ::  CM1PR(:,:)                         !precipitation  MgHPO4 concentration, [g m-3]
  real(r8),pointer ::  CCOR(:,:)                          !precipitation  CO2 concentration, [g m-3]
  real(r8),pointer ::  COXR(:,:)                          !precipitation  O2 concentration, [g m-3]
  real(r8),pointer ::  CNNR(:,:)                          !precipitation  N2 concentration, [g m-3]
  real(r8),pointer ::  CN2R(:,:)                          !precipitation  N2O concentration, [g m-3]
  contains
!----------------------------------------------------------------------

  subroutine InitClimForcData
  use TracerIDMod
  implicit none
  allocate(WDPTHD(366,JY,JX));  WDPTHD=0._r8
  allocate(TDTPX(JY,JX,12));    TDTPX=0._r8
  allocate(TDTPN(JY,JX,12));    TDTPN=0._r8
  allocate(TDRAD(JY,JX,12));    TDRAD=0._r8
  allocate(TDHUM(JY,JX,12));    TDHUM=0._r8
  allocate(TDPRC(JY,JX,12));    TDPRC=0._r8
  allocate(TDWND(JY,JX,12));    TDWND=0._r8
  allocate(TDCO2(JY,JX,12));    TDCO2=0._r8
  allocate(TDCN4(JY,JX,12));    TDCN4=0._r8
  allocate(TDCNO(JY,JX,12));    TDCNO=0._r8
  allocate(TCA(JY,JX));         TCA=0._r8
  allocate(TKA(JY,JX));         TKA=0._r8
  allocate(UA(JY,JX));          UA=0._r8
  allocate(VPA(JY,JX));         VPA=0._r8
  allocate(VPK(JY,JX));         VPK=0._r8
  allocate(DYLN(JY,JX));        DYLN=0._r8
  allocate(DYLX(JY,JX));        DYLX=0._r8
  allocate(DYLM(JY,JX));        DYLM=0._r8
  allocate(OMEGAG(JSA,JY,JX));  OMEGAG=0._r8
  allocate(THS(JY,JX));         THS=0._r8
  allocate(TRAD(JY,JX));        TRAD=0._r8
  allocate(TAMX(JY,JX));        TAMX=0._r8
  allocate(TAMN(JY,JX));        TAMN=0._r8
  allocate(HUDX(JY,JX));        HUDX=0._r8
  allocate(HUDN(JY,JX));        HUDN=0._r8
  allocate(TWIND(JY,JX));       TWIND=0._r8
  allocate(TRAI(JY,JX));        TRAI=0._r8
  allocate(THSX(JY,JX));        THSX=0._r8
  allocate(OFFSET(JY,JX));      OFFSET=0._r8
  allocate(PRECD(JY,JX));       PRECD=0._r8
  allocate(PRECB(JY,JX));       PRECB=0._r8
  allocate(CO2EI(JY,JX));       CO2EI=0._r8
  allocate(CCO2EI(JY,JX));      CCO2EI=0._r8

  allocate(AtmGgms(idg_beg:idg_end,JY,JX)); AtmGgms=0._r8
  allocate(AtmGmms(idg_beg:idg_end,JY,JX)); AtmGmms=0._r8

  allocate(OXYE(JY,JX));        OXYE=0._r8
  allocate(Z2OE(JY,JX));        Z2OE=0._r8
  allocate(Z2GE(JY,JX));        Z2GE=0._r8
  allocate(ZNH3E(JY,JX));       ZNH3E=0._r8
  allocate(CH4E(JY,JX));        CH4E=0._r8
  allocate(H2GE(JY,JX));        H2GE=0._r8
  allocate(CO2E(JY,JX));        CO2E=0._r8

  allocate(ZNOON(JY,JX));       ZNOON=0._r8
  allocate(RADS(JY,JX));        RADS=0._r8
  allocate(RADY(JY,JX));        RADY=0._r8
  allocate(RAPS(JY,JX));        RAPS=0._r8
  allocate(RAPY(JY,JX));        RAPY=0._r8
  allocate(SSIN(JY,JX));        SSIN=0._r8
  allocate(SSINN(JY,JX));       SSINN=0._r8
  allocate(TLEX(JY,JX));        TLEX=0._r8
  allocate(TSHX(JY,JX));        TSHX=0._r8
  allocate(TLEC(JY,JX));        TLEC=0._r8
  allocate(TSHC(JY,JX));        TSHC=0._r8
  allocate(DPTHSK(JY,JX));      DPTHSK=0._r8
  allocate(TKSD(JY,JX));        TKSD=0._r8
  allocate(ATCAI(JY,JX));       ATCAI=0._r8
  allocate(RAD(JY,JX));         RAD=0._r8
  allocate(RAP(JY,JX));         RAP=0._r8
  allocate(ATCA(JY,JX));        ATCA=0._r8
  allocate(ATCS(JY,JX));        ATCS=0._r8
  allocate(ATKA(JY,JX));        ATKA=0._r8
  allocate(ATKS(JY,JX));        ATKS=0._r8
  allocate(PRECR(JY,JX));       PRECR=0._r8
  allocate(PRECW(JY,JX));       PRECW=0._r8
  allocate(PRECQ(JY,JX));       PRECQ=0._r8
  allocate(PRECA(JY,JX));       PRECA=0._r8
  allocate(ENGYP(JY,JX));       ENGYP=0._r8
  allocate(PHR(JY,JX));         PHR=0._r8
  allocate(CN4RI(JY,JX));       CN4RI=0._r8
  allocate(CNORI(JY,JX));       CNORI=0._r8
  allocate(CN4R(JY,JX));        CN4R=0._r8
  allocate(CN3R(JY,JX));        CN3R=0._r8
  allocate(CNOR(JY,JX));        CNOR=0._r8
  allocate(CPOR(JY,JX));        CPOR=0._r8
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
  allocate(CCHR(JY,JX));        CCHR=0._r8
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
  allocate(CH1PR(JY,JX));       CH1PR=0._r8
  allocate(CH3PR(JY,JX));       CH3PR=0._r8
  allocate(CF1PR(JY,JX));       CF1PR=0._r8
  allocate(CF2PR(JY,JX));       CF2PR=0._r8
  allocate(CC0PR(JY,JX));       CC0PR=0._r8
  allocate(CC1PR(JY,JX));       CC1PR=0._r8
  allocate(CC2PR(JY,JX));       CC2PR=0._r8
  allocate(CM1PR(JY,JX));       CM1PR=0._r8
  allocate(CCOR(JY,JX));        CCOR=0._r8
  allocate(COXR(JY,JX));        COXR=0._r8
  allocate(CNNR(JY,JX));        CNNR=0._r8
  allocate(CN2R(JY,JX));        CN2R=0._r8
  end subroutine InitClimForcData

!----------------------------------------------------------------------
  subroutine DestructClimForcData
  use abortutils, only : destroy
  implicit none
  call destroy(WDPTHD)
  call destroy(TDTPX)
  call destroy(TDTPN)
  call destroy(TDRAD)
  call destroy(TDHUM)
  call destroy(TDPRC)
  call destroy(TDWND)
  call destroy(TDCO2)
  call destroy(TDCN4)
  call destroy(TDCNO)
  call destroy(TCA)
  call destroy(TKA)
  call destroy(UA)
  call destroy(VPA)
  call destroy(VPK)
  call destroy(DYLN)
  call destroy(DYLX)
  call destroy(DYLM)
  call destroy(OMEGAG)
  call destroy(THS)
  call destroy(TRAD)
  call destroy(TAMX)
  call destroy(TAMN)
  call destroy(HUDX)
  call destroy(HUDN)
  call destroy(TWIND)
  call destroy(TRAI)
  call destroy(THSX)
  call destroy(OFFSET)
  call destroy(PRECD)
  call destroy(PRECB)
  call destroy(CO2EI)
  call destroy(CCO2EI)

  call destroy(AtmGgms)
  call destroy(AtmGmms)

  call destroy(OXYE)
  call destroy(Z2OE)
  call destroy(Z2GE)
  call destroy(ZNH3E)
  call destroy(CH4E)
  call destroy(H2GE)

  call destroy(ZNOON)
  call destroy(CO2E)
  call destroy(RADS)
  call destroy(RADY)
  call destroy(RAPS)
  call destroy(RAPY)
  call destroy(SSIN)
  call destroy(SSINN)
  call destroy(TLEX)
  call destroy(TSHX)
  call destroy(TLEC)
  call destroy(TSHC)
  call destroy(DPTHSK)
  call destroy(TKSD)
  call destroy(ATCAI)
  call destroy(RAD)
  call destroy(RAP)
  call destroy(ATCA)
  call destroy(ATCS)
  call destroy(ATKA)
  call destroy(ATKS)
  call destroy(PRECR)
  call destroy(PRECW)
  call destroy(PRECQ)
  call destroy(PRECA)
  call destroy(ENGYP)
  call destroy(PHR)
  call destroy(CN4RI)
  call destroy(CNORI)
  call destroy(CN4R)
  call destroy(CN3R)
  call destroy(CNOR)
  call destroy(CPOR)
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
  call destroy(CCHR)
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
  call destroy(CH1PR)
  call destroy(CH3PR)
  call destroy(CF1PR)
  call destroy(CF2PR)
  call destroy(CC0PR)
  call destroy(CC1PR)
  call destroy(CC2PR)
  call destroy(CM1PR)
  call destroy(CCOR)
  call destroy(COXR)
  call destroy(CNNR)
  call destroy(CN2R)
  end subroutine DestructClimForcData
end module ClimForcDataType
