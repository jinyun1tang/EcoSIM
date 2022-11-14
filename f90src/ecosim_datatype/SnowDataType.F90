module SnowDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMCtrlMod, only : salt_model
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),target, allocatable ::  VHCPWM(:,:,:,:)                    !volumetric heat capacity of snowpack
  real(r8),target, allocatable ::  FLQWM(:,:,:,:)                     !snowpack water flux
  real(r8),target, allocatable ::  QSM(:,:,:,:)                       !runoff snow flux, [m3 d-2 t-1]
  REAL(R8),target, allocatable ::  ALBS(:,:)                          !snowpack albedo
  real(r8),target, allocatable ::  DENS0(:,:)                         !snowpack density, [Mg m-3]
  real(r8),target, allocatable ::  TCW(:,:,:)                         !snow temperature, [oC]
  real(r8),target, allocatable ::  TKW(:,:,:)                         !snow temperature, [K]
  real(r8),target, allocatable ::  VHCPW(:,:,:)                       !snowpack heat capacity, [MJ m-3 K-1]
  real(r8),target, allocatable ::  VOLSSL(:,:,:)                      !snow water equivalent volume in snowpack layer
  real(r8),target, allocatable ::  VOLWSL(:,:,:)                      !snow water volume in snowpack layer
  real(r8),target, allocatable ::  VOLISL(:,:,:)                      !snow ice volume in snowpack layer
  real(r8),target, allocatable ::  VOLSL(:,:,:)                       !snow volume in snowpack layer
  real(r8),target, allocatable ::  DENSS(:,:,:)                       !snowpack density, [Mg m-3]
  real(r8),target, allocatable ::  DLYRS(:,:,:)                       !snowpack layer depth
  real(r8),target, allocatable ::  XFLWW(:,:,:)                       !hourly snow water transfer
  real(r8),target, allocatable ::  XFLWS(:,:,:)                       !hourly snow transfer
  real(r8),target, allocatable ::  XFLWI(:,:,:)                       !hourly snow ice transfer
  real(r8),target, allocatable ::  XHFLWW(:,:,:)                      !hourly convective heat flux from water transfer
  real(r8),target, allocatable ::  XWFLXS(:,:,:)                      !hourly convective heat flux from snow transfer
  real(r8),target, allocatable ::  XWFLXI(:,:,:)                      !hourly convective heat flux from ice transfer
  real(r8),target, allocatable ::  CDPTHS(:,:,:)                      !cumulative depth to bottom of snowpack layer
  real(r8),target, allocatable ::  VOLSI(:,:,:)                       !Initial snowpack volume, [m3 d-2]
  real(r8),target, allocatable ::  DPTHS(:,:)                         !snowpack depth, [m]
  real(r8),target, allocatable ::  VOLSS(:,:)                         !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8),target, allocatable ::  VOLWS(:,:)                         !water volume in snowpack, [m3 d-2]
  real(r8),target, allocatable ::  VOLIS(:,:)                         !ice volume in snowpack, [m3 d-2]
  real(r8),target, allocatable ::  VOLS(:,:)                          !snowpack volume, [m3 d-2]
  real(r8),target, allocatable ::  VHCPWX(:,:)                        !snowpack heat capacity from previous time step, [MJ d-2 K-1]
  real(r8),target, allocatable ::  FLSW(:,:,:)                        !water from snowpack to soil micropores
  real(r8),target, allocatable ::  FLSWH(:,:,:)                       !water from snowpack to soil macropores
  real(r8),target, allocatable ::  HFLSW(:,:,:)                       !convective heat from snowpack to soil
  real(r8),target, allocatable ::  FLSWR(:,:,:)                       !water flux from snowpack to litter
  real(r8),target, allocatable ::  HFLSWR(:,:,:)                      !convective heat flux from snowpack to litter
  real(r8),target, allocatable ::  QS(:,:,:)                          !snowpack runoff snow, [m3 d-2 h-1]
  real(r8),target, allocatable ::  QW(:,:,:)                          !snowpack runoff water, [m3 d-2 h-1]
  real(r8),target, allocatable ::  QI(:,:,:)                          !snowpack runoff ice, [m3 d-2 h-1]
  real(r8),target, allocatable ::  HQS(:,:,:)                         !snowpack runoff heat, [MJ d-2 h-1]
  real(r8),target, allocatable ::  XCOQSS(:,:,:)                      !snowpack runoff CO2 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  XCHQSS(:,:,:)                      !snowpack runoff CH4 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  XOXQSS(:,:,:)                      !snowpack runoff O2 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  XNGQSS(:,:,:)                      !snowpack runoff N2 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  XN2QSS(:,:,:)                      !snowpack runoff N2O flux, [g d-2 h-1]
  real(r8),target, allocatable ::  XN4QSS(:,:,:)                      !snowpack runoff NH4 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  XN3QSS(:,:,:)                      !snowpack runoff NH3 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  XNOQSS(:,:,:)                      !snowpack runoff NO3 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  XP4QSS(:,:,:)                      !snowpack runoff PO4 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  XP1QSS(:,:,:)                      !snowpack runoff HPO4 flux, [g d-2 h-1]

  real(r8),target, allocatable ::  trcg_solsml(:,:,:,:)               ! snowpack dual phase disolved tracers
  real(r8),target, allocatable ::  trcn_solsml(:,:,:,:)               ! snowpack nutrient dissolved tracers
  real(r8),target, allocatable ::  trcs_solsml(:,:,:,:)               ! snowpack salt dissolved tracers

  real(r8),target, allocatable ::  CO2W(:,:,:)                        !snowpack CO2, [mol d-2]
  real(r8),target, allocatable ::  CH4W(:,:,:)                        !snowpack CH4, [mol d-2]
  real(r8),target, allocatable ::  OXYW(:,:,:)                        !snowpack O2, [mol d-2]
  real(r8),target, allocatable ::  ZN2W(:,:,:)                        !snowpack N2O, [mol d-2]
  real(r8),target, allocatable ::  ZNGW(:,:,:)                        !snowpack N2, [mol d-2]
  real(r8),target, allocatable ::  ZN3W(:,:,:)                        !snowpack NH3, [mol d-2]

  real(r8),target, allocatable ::  ZN4W(:,:,:)                        !snowpack NH4, [mol d-2]
  real(r8),target, allocatable ::  ZNOW(:,:,:)                        !snowpack NO3, [mol d-2]
  real(r8),target, allocatable ::  Z1PW(:,:,:)                        !snowpack HPO4,[mol d-2]
  real(r8),target, allocatable ::  ZHPW(:,:,:)                        !snowpack H2PO4, [mol d-2]

  real(r8),target, allocatable ::  ZALW(:,:,:)                        !snowpack Al, [mol d-2]
  real(r8),target, allocatable ::  ZFEW(:,:,:)                        !snowpack Fe, [mol d-2]
  real(r8),target, allocatable ::  ZHYW(:,:,:)                        !snowpack H, [mol d-2]
  real(r8),target, allocatable ::  ZCAW(:,:,:)                        !snowpack Ca, [mol d-2]
  real(r8),target, allocatable ::  ZMGW(:,:,:)                        !snowpack Mg, [mol d-2]
  real(r8),target, allocatable ::  ZNAW(:,:,:)                        !snowpack Na, [mol d-2]
  real(r8),target, allocatable ::  ZKAW(:,:,:)                        !snowpack K, [mol d-2]
  real(r8),target, allocatable ::  ZOHW(:,:,:)                        !snowpack OH, [mol d-2]
  real(r8),target, allocatable ::  ZSO4W(:,:,:)                       !snowpack SO4, [mol d-2]
  real(r8),target, allocatable ::  ZCLW(:,:,:)                        !snowpack Cl, [mol d-2]
  real(r8),target, allocatable ::  ZCO3W(:,:,:)                       !snowpack CO3, [mol d-2]
  real(r8),target, allocatable ::  ZHCO3W(:,:,:)                      !snowpack HCO3, [mol d-2]
  real(r8),target, allocatable ::  ZALH1W(:,:,:)                      !snowpack AlOH, [mol d-2]
  real(r8),target, allocatable ::  ZALH2W(:,:,:)                      !snowpack AlOH2, [mol d-2]
  real(r8),target, allocatable ::  ZALH3W(:,:,:)                      !snowpack AlOH3, [mol d-2]
  real(r8),target, allocatable ::  ZALH4W(:,:,:)                      !snowpack AlOH4, [mol d-2]
  real(r8),target, allocatable ::  ZALSW(:,:,:)                       !snowpack AlSO4, [mol d-2]
  real(r8),target, allocatable ::  ZFEH1W(:,:,:)                      !snowpack FeOH, [mol d-2]
  real(r8),target, allocatable ::  ZFEH2W(:,:,:)                      !snowpack FeOH2, [mol d-2]
  real(r8),target, allocatable ::  ZFEH3W(:,:,:)                      !snowpack FeOH3, [mol d-2]
  real(r8),target, allocatable ::  ZFEH4W(:,:,:)                      !snowpack F3OH4, [mol d-2]
  real(r8),target, allocatable ::  ZFESW(:,:,:)                       !snowpack FeSO4, [mol d-2]
  real(r8),target, allocatable ::  ZCAOW(:,:,:)                       !snowpack CaOH2, [mol d-2]
  real(r8),target, allocatable ::  ZCACW(:,:,:)                       !snowpack CaCO3, [mol d-2]
  real(r8),target, allocatable ::  ZCAHW(:,:,:)                       !snowpack CaHCO3, [mol d-2]
  real(r8),target, allocatable ::  ZCASW(:,:,:)                       !snowpack CaSO4, [mol d-2]
  real(r8),target, allocatable ::  ZMGOW(:,:,:)                       !snowpack MgOH2, [mol d-2]
  real(r8),target, allocatable ::  ZMGCW(:,:,:)                       !snowpack MgCO3, [mol d-2]
  real(r8),target, allocatable ::  ZMGHW(:,:,:)                       !snowpack MgHCO3, [mol d-2]
  real(r8),target, allocatable ::  ZMGSW(:,:,:)                       !snowpack MgSO4, [mol d-2]
  real(r8),target, allocatable ::  ZNACW(:,:,:)                       !snowpack NaCO3, [mol d-2]
  real(r8),target, allocatable ::  ZNASW(:,:,:)                       !snowpack NaSO4, [mol d-2]
  real(r8),target, allocatable ::  ZKASW(:,:,:)                       !snowpack KSO4, [mol d-2]
  real(r8),target, allocatable ::  H0PO4W(:,:,:)                      !snowpack PO4, [mol d-2]
  real(r8),target, allocatable ::  H3PO4W(:,:,:)                      !snowpack H3PO4, [mol d-2]
  real(r8),target, allocatable ::  ZFE1PW(:,:,:)                      !snowpack FeHPO4, [mol d-2]
  real(r8),target, allocatable ::  ZFE2PW(:,:,:)                      !snowpack FeH2PO4, [mol d-2]
  real(r8),target, allocatable ::  ZCA0PW(:,:,:)                      !snowpack CaPO4, [mol d-2]
  real(r8),target, allocatable ::  ZCA1PW(:,:,:)                      !snowpack CaHPO4, [mol d-2]
  real(r8),target, allocatable ::  ZCA2PW(:,:,:)                      !snowpack CaH2PO4, [mol d-2]
  real(r8),target, allocatable ::  ZMG1PW(:,:,:)                      !snowpack MgHPO4, [mol d-2]

  real(r8),target, allocatable ::  XQSAL(:,:,:)                       !total Al in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSFE(:,:,:)                       !total Fe in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSHY(:,:,:)                       !total H in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSCA(:,:,:)                       !total Ca in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSMG(:,:,:)                       !total Mg in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSNA(:,:,:)                       !total Na in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSKA(:,:,:)                       !total K in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSOH(:,:,:)                       !total OH in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSSO(:,:,:)                       !total SO4 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSCL(:,:,:)                       !total Cl in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSC3(:,:,:)                       !total CO3 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSHC(:,:,:)                       !total HCO3 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSAL1(:,:,:)                      !total AlOH in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSAL2(:,:,:)                      !total AlOH2 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSAL3(:,:,:)                      !total AlOH3 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSAL4(:,:,:)                      !total AlOH4 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSALS(:,:,:)                      !total AlSO4 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSFE1(:,:,:)                      !total FeOH in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSFE2(:,:,:)                      !total FeOH2 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSFE3(:,:,:)                      !total FeOH3 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSFE4(:,:,:)                      !total FeOH4 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSFES(:,:,:)                      !total FeSO4 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSCAO(:,:,:)                      !total CaOH in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSCAC(:,:,:)                      !total CaCO3 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSCAH(:,:,:)                      !total CaHCO3 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSCAS(:,:,:)                      !total CaSO4 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSMGO(:,:,:)                      !total MgOH in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSMGC(:,:,:)                      !total MgCO3 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSMGH(:,:,:)                      !total MgHCO3 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSMGS(:,:,:)                      !total MgSO4 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSNAC(:,:,:)                      !total NaCO3 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSNAS(:,:,:)                      !total NaSO4 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSKAS(:,:,:)                      !total KSO4 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSH0P(:,:,:)                      !total PO4 in snow drif, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSH1P(:,:,:)                      !total HPO4 in snow drift , [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSH3P(:,:,:)                      !total H3PO4 in snow drift , [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSF1P(:,:,:)                      !total FeHPO4 in snow drift , [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSF2P(:,:,:)                      !total FeH2PO4 in snow drift , [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSC0P(:,:,:)                      !total CaPO4 in snow drift , [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSC1P(:,:,:)                      !total CaHPO4 in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSC2P(:,:,:)                      !total CaH2PO4 in snow drift , [mol d-2 h-1]
  real(r8),target, allocatable ::  XQSM1P(:,:,:)                      !total MgHPO4 in snow drift , [mol d-2 h-1]
!----------------------------------------------------------------------

contains
  subroutine InitSnowData

  implicit none
  allocate(VHCPWM(60,JS,JY,JX));VHCPWM=0._r8
  allocate(FLQWM(60,JS,JY,JX)); FLQWM=0._r8
  allocate(QSM(60,2,JV,JH));    QSM=0._r8
  allocate(ALBS(JY,JX));        ALBS=0._r8
  allocate(DENS0(JY,JX));       DENS0=0._r8
  allocate(TCW(JS,JY,JX));      TCW=0._r8
  allocate(TKW(JS,JY,JX));      TKW=0._r8
  allocate(VHCPW(JS,JY,JX));    VHCPW=0._r8
  allocate(VOLSSL(JS,JY,JX));   VOLSSL=0._r8
  allocate(VOLWSL(JS,JY,JX));   VOLWSL=0._r8
  allocate(VOLISL(JS,JY,JX));   VOLISL=0._r8
  allocate(VOLSL(JS,JY,JX));    VOLSL=0._r8
  allocate(DENSS(JS,JY,JX));    DENSS=0._r8
  allocate(DLYRS(JS,JY,JX));    DLYRS=0._r8
  allocate(XFLWW(JS,JY,JX));    XFLWW=0._r8
  allocate(XFLWS(JS,JY,JX));    XFLWS=0._r8
  allocate(XFLWI(JS,JY,JX));    XFLWI=0._r8
  allocate(XHFLWW(JS,JY,JX));   XHFLWW=0._r8
  allocate(XWFLXS(JS,JY,JX));   XWFLXS=0._r8
  allocate(XWFLXI(JS,JY,JX));   XWFLXI=0._r8
  allocate(CDPTHS(0:JS,JY,JX)); CDPTHS=0._r8
  allocate(VOLSI(JS,JY,JX));    VOLSI=0._r8
  allocate(DPTHS(JY,JX));       DPTHS=0._r8
  allocate(VOLSS(JY,JX));       VOLSS=0._r8
  allocate(VOLWS(JY,JX));       VOLWS=0._r8
  allocate(VOLIS(JY,JX));       VOLIS=0._r8
  allocate(VOLS(JY,JX));        VOLS=0._r8
  allocate(VHCPWX(JY,JX));      VHCPWX=0._r8
  allocate(FLSW(JS,JY,JX));     FLSW=0._r8
  allocate(FLSWH(JS,JY,JX));    FLSWH=0._r8
  allocate(HFLSW(JS,JY,JX));    HFLSW=0._r8
  allocate(FLSWR(JS,JY,JX));    FLSWR=0._r8
  allocate(HFLSWR(JS,JY,JX));   HFLSWR=0._r8
  allocate(QS(2,JV,JH));        QS=0._r8
  allocate(QW(2,JV,JH));        QW=0._r8
  allocate(QI(2,JV,JH));        QI=0._r8
  allocate(HQS(2,JV,JH));       HQS=0._r8
  allocate(XCOQSS(2,JV,JH));    XCOQSS=0._r8
  allocate(XCHQSS(2,JV,JH));    XCHQSS=0._r8
  allocate(XOXQSS(2,JV,JH));    XOXQSS=0._r8
  allocate(XNGQSS(2,JV,JH));    XNGQSS=0._r8
  allocate(XN2QSS(2,JV,JH));    XN2QSS=0._r8
  allocate(XN4QSS(2,JV,JH));    XN4QSS=0._r8
  allocate(XN3QSS(2,JV,JH));    XN3QSS=0._r8
  allocate(XNOQSS(2,JV,JH));    XNOQSS=0._r8
  allocate(XP4QSS(2,JV,JH));    XP4QSS=0._r8
  allocate(XP1QSS(2,JV,JH));    XP1QSS=0._r8

! exclude NH3B
  allocate(trcg_solsml(idg_beg:idg_end-1,JS,JY,JX));trcg_solsml=0._r8
  allocate(trcn_solsml(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_solsml=0._r8
  if(salt_model)then
    allocate(trcs_solsml(idsa_beg:idsa_end,JS,JY,JX)); trcs_solsml=0._r8
  endif

  allocate(CO2W(JS,JY,JX));     CO2W=0._r8
  allocate(CH4W(JS,JY,JX));     CH4W=0._r8
  allocate(OXYW(JS,JY,JX));     OXYW=0._r8
  allocate(ZN2W(JS,JY,JX));     ZN2W=0._r8
  allocate(ZNGW(JS,JY,JX));     ZNGW=0._r8
  allocate(ZN4W(JS,JY,JX));     ZN4W=0._r8
  allocate(ZN3W(JS,JY,JX));     ZN3W=0._r8
  allocate(ZNOW(JS,JY,JX));     ZNOW=0._r8
  allocate(Z1PW(JS,JY,JX));     Z1PW=0._r8
  allocate(ZHPW(JS,JY,JX));     ZHPW=0._r8
  if(salt_model)then
    allocate(ZALW(JS,JY,JX));     ZALW=0._r8
    allocate(ZFEW(JS,JY,JX));     ZFEW=0._r8
    allocate(ZHYW(JS,JY,JX));     ZHYW=0._r8
    allocate(ZCAW(JS,JY,JX));     ZCAW=0._r8
    allocate(ZMGW(JS,JY,JX));     ZMGW=0._r8
    allocate(ZNAW(JS,JY,JX));     ZNAW=0._r8
    allocate(ZKAW(JS,JY,JX));     ZKAW=0._r8
    allocate(ZOHW(JS,JY,JX));     ZOHW=0._r8
    allocate(ZSO4W(JS,JY,JX));    ZSO4W=0._r8
    allocate(ZCLW(JS,JY,JX));     ZCLW=0._r8
    allocate(ZCO3W(JS,JY,JX));    ZCO3W=0._r8
    allocate(ZHCO3W(JS,JY,JX));   ZHCO3W=0._r8
    allocate(ZALH1W(JS,JY,JX));   ZALH1W=0._r8
    allocate(ZALH2W(JS,JY,JX));   ZALH2W=0._r8
    allocate(ZALH3W(JS,JY,JX));   ZALH3W=0._r8
    allocate(ZALH4W(JS,JY,JX));   ZALH4W=0._r8
    allocate(ZALSW(JS,JY,JX));    ZALSW=0._r8
    allocate(ZFEH1W(JS,JY,JX));   ZFEH1W=0._r8
    allocate(ZFEH2W(JS,JY,JX));   ZFEH2W=0._r8
    allocate(ZFEH3W(JS,JY,JX));   ZFEH3W=0._r8
    allocate(ZFEH4W(JS,JY,JX));   ZFEH4W=0._r8
    allocate(ZFESW(JS,JY,JX));    ZFESW=0._r8
    allocate(ZCAOW(JS,JY,JX));    ZCAOW=0._r8
    allocate(ZCACW(JS,JY,JX));    ZCACW=0._r8
    allocate(ZCAHW(JS,JY,JX));    ZCAHW=0._r8
    allocate(ZCASW(JS,JY,JX));    ZCASW=0._r8
    allocate(ZMGOW(JS,JY,JX));    ZMGOW=0._r8
    allocate(ZMGCW(JS,JY,JX));    ZMGCW=0._r8
    allocate(ZMGHW(JS,JY,JX));    ZMGHW=0._r8
    allocate(ZMGSW(JS,JY,JX));    ZMGSW=0._r8
    allocate(ZNACW(JS,JY,JX));    ZNACW=0._r8
    allocate(ZNASW(JS,JY,JX));    ZNASW=0._r8
    allocate(ZKASW(JS,JY,JX));    ZKASW=0._r8
    allocate(H0PO4W(JS,JY,JX));   H0PO4W=0._r8
    allocate(H3PO4W(JS,JY,JX));   H3PO4W=0._r8
    allocate(ZFE1PW(JS,JY,JX));   ZFE1PW=0._r8
    allocate(ZFE2PW(JS,JY,JX));   ZFE2PW=0._r8
    allocate(ZCA0PW(JS,JY,JX));   ZCA0PW=0._r8
    allocate(ZCA1PW(JS,JY,JX));   ZCA1PW=0._r8
    allocate(ZCA2PW(JS,JY,JX));   ZCA2PW=0._r8
    allocate(ZMG1PW(JS,JY,JX));   ZMG1PW=0._r8
  endif
  allocate(XQSAL(2,JV,JH));     XQSAL=0._r8
  allocate(XQSFE(2,JV,JH));     XQSFE=0._r8
  allocate(XQSHY(2,JV,JH));     XQSHY=0._r8
  allocate(XQSCA(2,JV,JH));     XQSCA=0._r8
  allocate(XQSMG(2,JV,JH));     XQSMG=0._r8
  allocate(XQSNA(2,JV,JH));     XQSNA=0._r8
  allocate(XQSKA(2,JV,JH));     XQSKA=0._r8
  allocate(XQSOH(2,JV,JH));     XQSOH=0._r8
  allocate(XQSSO(2,JV,JH));     XQSSO=0._r8
  allocate(XQSCL(2,JV,JH));     XQSCL=0._r8
  allocate(XQSC3(2,JV,JH));     XQSC3=0._r8
  allocate(XQSHC(2,JV,JH));     XQSHC=0._r8
  allocate(XQSAL1(2,JV,JH));    XQSAL1=0._r8
  allocate(XQSAL2(2,JV,JH));    XQSAL2=0._r8
  allocate(XQSAL3(2,JV,JH));    XQSAL3=0._r8
  allocate(XQSAL4(2,JV,JH));    XQSAL4=0._r8
  allocate(XQSALS(2,JV,JH));    XQSALS=0._r8
  allocate(XQSFE1(2,JV,JH));    XQSFE1=0._r8
  allocate(XQSFE2(2,JV,JH));    XQSFE2=0._r8
  allocate(XQSFE3(2,JV,JH));    XQSFE3=0._r8
  allocate(XQSFE4(2,JV,JH));    XQSFE4=0._r8
  allocate(XQSFES(2,JV,JH));    XQSFES=0._r8
  allocate(XQSCAO(2,JV,JH));    XQSCAO=0._r8
  allocate(XQSCAC(2,JV,JH));    XQSCAC=0._r8
  allocate(XQSCAH(2,JV,JH));    XQSCAH=0._r8
  allocate(XQSCAS(2,JV,JH));    XQSCAS=0._r8
  allocate(XQSMGO(2,JV,JH));    XQSMGO=0._r8
  allocate(XQSMGC(2,JV,JH));    XQSMGC=0._r8
  allocate(XQSMGH(2,JV,JH));    XQSMGH=0._r8
  allocate(XQSMGS(2,JV,JH));    XQSMGS=0._r8
  allocate(XQSNAC(2,JV,JH));    XQSNAC=0._r8
  allocate(XQSNAS(2,JV,JH));    XQSNAS=0._r8
  allocate(XQSKAS(2,JV,JH));    XQSKAS=0._r8
  allocate(XQSH0P(2,JV,JH));    XQSH0P=0._r8
  allocate(XQSH1P(2,JV,JH));    XQSH1P=0._r8
  allocate(XQSH3P(2,JV,JH));    XQSH3P=0._r8
  allocate(XQSF1P(2,JV,JH));    XQSF1P=0._r8
  allocate(XQSF2P(2,JV,JH));    XQSF2P=0._r8
  allocate(XQSC0P(2,JV,JH));    XQSC0P=0._r8
  allocate(XQSC1P(2,JV,JH));    XQSC1P=0._r8
  allocate(XQSC2P(2,JV,JH));    XQSC2P=0._r8
  allocate(XQSM1P(2,JV,JH));    XQSM1P=0._r8
  end subroutine InitSnowData

!----------------------------------------------------------------------
  subroutine DestructSnowData
  use abortutils, only : destroy
  implicit none
  call destroy(VHCPWM)
  call destroy(FLQWM)
  call destroy(QSM)
  call destroy(ALBS)
  call destroy(DENS0)
  call destroy(TCW)
  call destroy(TKW)
  call destroy(VHCPW)
  call destroy(VOLSSL)
  call destroy(VOLWSL)
  call destroy(VOLISL)
  call destroy(VOLSL)
  call destroy(DENSS)
  call destroy(DLYRS)
  call destroy(XFLWW)
  call destroy(XFLWS)
  call destroy(XFLWI)
  call destroy(XHFLWW)
  call destroy(XWFLXS)
  call destroy(XWFLXI)
  call destroy(CDPTHS)
  call destroy(VOLSI)
  call destroy(DPTHS)
  call destroy(VOLSS)
  call destroy(VOLWS)
  call destroy(VOLIS)
  call destroy(VOLS)
  call destroy(VHCPWX)
  call destroy(FLSW)
  call destroy(FLSWH)
  call destroy(HFLSW)
  call destroy(FLSWR)
  call destroy(HFLSWR)
  call destroy(QS)
  call destroy(QW)
  call destroy(QI)
  call destroy(HQS)
  call destroy(XCOQSS)
  call destroy(XCHQSS)
  call destroy(XOXQSS)
  call destroy(XNGQSS)
  call destroy(XN2QSS)
  call destroy(XN4QSS)
  call destroy(XN3QSS)
  call destroy(XNOQSS)
  call destroy(XP4QSS)
  call destroy(XP1QSS)
  call destroy(CO2W)
  call destroy(CH4W)
  call destroy(OXYW)
  call destroy(ZN2W)
  call destroy(ZNGW)
  call destroy(ZN4W)
  call destroy(ZN3W)
  call destroy(ZNOW)
  call destroy(Z1PW)
  call destroy(ZHPW)

  call destroy(ZALW)
  call destroy(ZFEW)
  call destroy(ZHYW)
  call destroy(ZCAW)
  call destroy(ZMGW)
  call destroy(ZNAW)
  call destroy(ZKAW)
  call destroy(ZOHW)
  call destroy(ZSO4W)
  call destroy(ZCLW)
  call destroy(ZCO3W)
  call destroy(ZHCO3W)
  call destroy(ZALH1W)
  call destroy(ZALH2W)
  call destroy(ZALH3W)
  call destroy(ZALH4W)
  call destroy(ZALSW)
  call destroy(ZFEH1W)
  call destroy(ZFEH2W)
  call destroy(ZFEH3W)
  call destroy(ZFEH4W)
  call destroy(ZFESW)
  call destroy(ZCAOW)
  call destroy(ZCACW)
  call destroy(ZCAHW)
  call destroy(ZCASW)
  call destroy(ZMGOW)
  call destroy(ZMGCW)
  call destroy(ZMGHW)
  call destroy(ZMGSW)
  call destroy(ZNACW)
  call destroy(ZNASW)
  call destroy(ZKASW)
  call destroy(H0PO4W)
  call destroy(H3PO4W)
  call destroy(ZFE1PW)
  call destroy(ZFE2PW)
  call destroy(ZCA0PW)
  call destroy(ZCA1PW)
  call destroy(ZCA2PW)
  call destroy(ZMG1PW)
  
  call destroy(XQSAL)
  call destroy(XQSFE)
  call destroy(XQSHY)
  call destroy(XQSCA)
  call destroy(XQSMG)
  call destroy(XQSNA)
  call destroy(XQSKA)
  call destroy(XQSOH)
  call destroy(XQSSO)
  call destroy(XQSCL)
  call destroy(XQSC3)
  call destroy(XQSHC)
  call destroy(XQSAL1)
  call destroy(XQSAL2)
  call destroy(XQSAL3)
  call destroy(XQSAL4)
  call destroy(XQSALS)
  call destroy(XQSFE1)
  call destroy(XQSFE2)
  call destroy(XQSFE3)
  call destroy(XQSFE4)
  call destroy(XQSFES)
  call destroy(XQSCAO)
  call destroy(XQSCAC)
  call destroy(XQSCAH)
  call destroy(XQSCAS)
  call destroy(XQSMGO)
  call destroy(XQSMGC)
  call destroy(XQSMGH)
  call destroy(XQSMGS)
  call destroy(XQSNAC)
  call destroy(XQSNAS)
  call destroy(XQSKAS)
  call destroy(XQSH0P)
  call destroy(XQSH1P)
  call destroy(XQSH3P)
  call destroy(XQSF1P)
  call destroy(XQSF2P)
  call destroy(XQSC0P)
  call destroy(XQSC1P)
  call destroy(XQSC2P)
  call destroy(XQSM1P)
  end subroutine DestructSnowData

end module SnowDataType
