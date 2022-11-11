module SoilBGCDataType

!
! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use ElmIDMod    , only : npelms
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken => jskenc
implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),pointer ::  CNH4(:,:,:)                       !soil NH4 content, [mg kg-1]
  real(r8),pointer ::  CNO3(:,:,:)                       !soil NO3 content, [mg kg-1]
  real(r8),pointer ::  CPO4(:,:,:)                       !soil PO4 content, [mg kg-1]

  real(r8),pointer :: CPO4B(:,:,:)          !PO4 concentration band micropore	[g m-3]
  real(r8),pointer :: trc_soHml(:,:,:,:)                 !solute mass in macropore [g d-2]
  real(r8),pointer :: trc_solml(:,:,:,:)                 !solute mass in micropore [g d-2]
  real(r8),pointer :: trc_solcl(:,:,:,:)                 !solute concentration in micropre [g m-3]
  real(r8),pointer :: trc_gascl(:,:,:,:)                 !gaseous concentation [g m-3]

  real(r8),pointer ::  ZNFNI(:,:,:)                      !current nitrification inhibition activity
  real(r8),pointer ::  ZNFN0(:,:,:)                      !initial nitrification inhibition activity
  real(r8),pointer ::  ZNHUI(:,:,:)                      !current inhibition activity
  real(r8),pointer ::  ZNHU0(:,:,:)                      !urea hydrolysis inhibition activity

  real(r8),pointer :: CNH3G(:,:,:)          !gaseous NH3 concentration	[g m-3]
  real(r8),pointer :: CZ2GG(:,:,:)          !gaseous N2 concentration	[g m-3]
  real(r8),pointer :: CZ2OG(:,:,:)          !gaseous N2O concentration	[g m-3]
  real(r8),pointer :: COXYG(:,:,:)          !gaseous O2 concentration	[g m-3]
  real(r8),pointer :: CCH4G(:,:,:)          !gaseous CH4 concentration	[g m-3]

  real(r8),pointer :: CCO2G(:,:,:)          !gaseous CO2 concentration	[g m-3]

  real(r8),pointer :: CH1P4(:,:,:)          !aqueous H1PO4 concentration non-band [g m-3]


  real(r8),pointer :: CH2GG(:,:,:)          !gaseous H2 concentration	[g m-3]

  real(r8),pointer :: trc_gasml(:,:,:,:)     !layer mass of gases [g d-2]

  real(r8),pointer :: CPO4S(:,:,:)          !PO4 concentration non-band micropore	[g m-3]
  real(r8),pointer :: PH(:,:,:)             !soil pH
  real(r8),pointer :: CEC(:,:,:)            !soil cation exchange capacity	[cmol kg-1]
  real(r8),pointer :: AEC(:,:,:)            !soil anion exchange capacity	[cmol kg-1]

  real(r8),pointer ::  ROXSK(:,:,:,:)                     !total O2 sink, [g d-2 t-1]
  real(r8),pointer ::  HCO2G(:,:)                         !soil CO2 flux, [g d-2 h-1]
  real(r8),pointer ::  HCH4G(:,:)                         !soil CH4 flux, [g d-2 h-1]
  real(r8),pointer ::  HOXYG(:,:)                         !soil O2 flux, [g d-2 h-1]
  real(r8),pointer ::  HN2OG(:,:)                         !soil N2O flux, [g d-2 h-1]
  real(r8),pointer ::  HNH3G(:,:)                         !soil NH3 flux, [g d-2 h-1]
  real(r8),pointer ::  UORGF(:,:)                         !total C amendment, [g d-2]
  real(r8),pointer ::  UFERTN(:,:)                        !total fertilizer N amendment, [g d-2]
  real(r8),pointer ::  UFERTP(:,:)                        !total fertilizer P amendment, [g d-2]
  real(r8),pointer ::  UDOCQ(:,:)                         !total surface DOC flux, [g d-2]
  real(r8),pointer ::  UDOCD(:,:)                         !total subsurface DOC flux, [g d-2]
  real(r8),pointer ::  UXCSN(:,:)                         !total litterfall C, [g d-2]
  real(r8),pointer ::  UXZSN(:,:)                         !total litterfall N, [g d-2]
  real(r8),pointer ::  UXPSN(:,:)                         !total litterfall P, [g d-2]
  real(r8),pointer ::  UDONQ(:,:)                         !total surface DON flux, [g d-2]
  real(r8),pointer ::  UDOND(:,:)                         !total subsurface DON flux, [g d-2]
  real(r8),pointer ::  UDOPQ(:,:)                         !total surface DOP flux, [g d-2]
  real(r8),pointer ::  UDOPD(:,:)                         !total subsurface DOP flux, [g d-2]
  real(r8),pointer ::  UPP4(:,:)                          !total soil precipited P, [g d-2]
  real(r8),pointer ::  UN2GS(:,:)                         !total N2 fixation, [g d-2]
  real(r8),pointer ::  UH2GG(:,:)                         !total H2 flux, []
  real(r8),pointer ::  HN2GG(:,:)                         !soil N2 flux, [g d-2 h-1]
  real(r8),pointer ::  UN2GG(:,:)                         !total soil N2 flux, [g d-2]
  real(r8),pointer ::  UCO2G(:,:)                         !total soil CO2 flux, [g d-2]
  real(r8),pointer ::  UCH4G(:,:)                         !total soil CH4 flux, [g d-2]
  real(r8),pointer ::  UOXYG(:,:)                         !total soil O2 flux, [g d-2]
  real(r8),pointer ::  UNH3G(:,:)                         !total soil NH3 flux, [g d-2]
  real(r8),pointer ::  UN2OG(:,:)                         !total soil N2O flux, [g d-2]
  real(r8),pointer ::  UCOP(:,:)                          !total soil autotrophic respiration, [g d-2]
  real(r8),pointer ::  USEDOU(:,:)                        !total sediment subsurface flux, [Mg d-2]
  real(r8),pointer ::  UDICQ(:,:)                         !total surface DIC flux, [g d-2]
  real(r8),pointer ::  UDICD(:,:)                         !total subsurface DIC flux, [g d-2]
  real(r8),pointer ::  UDINQ(:,:)                         !total surface DIN flux, [g d-2]
  real(r8),pointer ::  UDIND(:,:)                         !total subsurface DIN flux, [g d-2]
  real(r8),pointer ::  UDIPQ(:,:)                         !total surface DIP flux, [g d-2]
  real(r8),pointer ::  UDIPD(:,:)                         !total subsurface DIP flux, [g d-2]
  real(r8),pointer ::  WTSTGT(:,:)                        !total standing dead C, [g d-2]
  real(r8),pointer ::  ZDRAIN(:,:)                        !total N drainage below root zone, [g d-2]
  real(r8),pointer ::  PDRAIN(:,:)                        !total P drainage below root zone, [g d-2]
  real(r8),pointer ::  UION(:,:)                          !total soil ion content, [mol d-2]
  real(r8),pointer ::  UIONOU(:,:)                        !total subsurface ion flux, [mol d-2]
  real(r8),pointer ::  XNO2S(:,:,:)                       !total NO2 exchange, [g d-2 h-1]
  real(r8),pointer ::  RUPOXO(:,:,:)                      !microbial O2 uptake, [g d-2 h-1]
  real(r8),pointer ::  RCO2O(:,:,:)                       !microbial net CO2 exchange, [g d-2 h-1]
  real(r8),pointer ::  RCH4O(:,:,:)                       !microbial net CH4 exchange, [g d-2 h-1]
  real(r8),pointer ::  RH2GO(:,:,:)                       !microbial net H2 exchange, [g d-2 h-1]
  real(r8),pointer ::  RN2G(:,:,:)                        !microbial net N2 exchange, [g d-2 h-1]
  real(r8),pointer ::  RN2O(:,:,:)                        !microbial net N2O exchange, [g d-2 h-1]
  real(r8),pointer ::  XNO2B(:,:,:)                       !net microbial NO2 exchange band, [g d-2 h-1]
  real(r8),pointer ::  XNH4B(:,:,:)                       !net microbial NH4 exchange band, [g d-2 h-1]
  real(r8),pointer ::  XNO3B(:,:,:)                       !net microbial NO3 exchange band, [g d-2 h-1]
  real(r8),pointer ::  XH2BS(:,:,:)                       !net microbial PO4 exchange band, [g d-2 h-1]
  real(r8),pointer ::  XH1BS(:,:,:)                       !net microbial HPO4 exchange band, [g d-2 h-1]
  real(r8),pointer ::  XN2GS(:,:,:)                       !net microbial N2 exchange, [g d-2 h-1]
  real(r8),pointer ::  XH1PS(:,:,:)                       !net microibal HPO4 exchange non-band, [g d-2 h-1]
  real(r8),pointer ::  XH2PS(:,:,:)                       !net microbial PO4 exchange nonband, [g d-2 h-1]
  real(r8),pointer ::  XNH4S(:,:,:)                       !net microbial NH4 exchange non-band, [g d-2 h-1]
  real(r8),pointer ::  XNO3S(:,:,:)                       !net microbial NO3 exchange non-band, [g d-2 h-1]
  real(r8),pointer ::  XOQCS(:,:,:,:)                     !net microbial DOC flux, [g d-2 h-1]
  real(r8),pointer ::  XOQNS(:,:,:,:)                     !net microbial DON flux, [g d-2 h-1]
  real(r8),pointer ::  XOQPS(:,:,:,:)                     !net microbial DOP flux, [g d-2 h-1]
  real(r8),pointer ::  XOQAS(:,:,:,:)                     !net microbial acetate flux, [g d-2 h-1]
  real(r8),pointer ::  TOQCK(:,:,:)                       !total respiration of DOC+DOA in soil layer
  real(r8),pointer ::  VOLQ(:,:,:)                        !soil water volume occupied by microial biomass, [m3 m-3]
  real(r8),pointer ::  TFNQ(:,:,:)                        !constraints of temperature and water potential on microbial activity, []
  real(r8),pointer ::  ESNT(:,:,:,:,:,:)                    !total litterfall C, [g d-2 h-1]
  real(r8),pointer ::  VLNH4(:,:,:)                       !NH4 non-band volume fracrion, []
  real(r8),pointer ::  VLNHB(:,:,:)                       !NH4 band volume fracrion, []
  real(r8),pointer ::  VLNO3(:,:,:)                       !NO3 non-band volume fracrion, []
  real(r8),pointer ::  VLNOB(:,:,:)                       !NO3 band volume fracrion, []
  real(r8),pointer ::  VLPO4(:,:,:)                       !PO4 non-band volume fracrion, []
  real(r8),pointer ::  VLPOB(:,:,:)                       !PO4 band volume fracrion, []
  real(r8),pointer ::  WDNHB(:,:,:)                       !width of NH4 band, [m]
  real(r8),pointer ::  DPNHB(:,:,:)                       !depth of NH4 band, [m]
  real(r8),pointer ::  WDNOB(:,:,:)                       !width of NO3 band, [m]
  real(r8),pointer ::  DPNOB(:,:,:)                       !depth of NO4 band, [m]
  real(r8),pointer ::  WDPOB(:,:,:)                       !width of PO4 band, [m]
  real(r8),pointer ::  DPPOB(:,:,:)                       !depth of PO4 band, [m]
  real(r8),pointer ::  DPNH4(:,:)                         !total depth of NH4 band, [m]
  real(r8),pointer ::  DPNO3(:,:)                         !total depth of NO3 band, [m]
  real(r8),pointer ::  DPPO4(:,:)                         !total depth of PO4 band, [m]
  real(r8),pointer ::  RVMXC(:,:,:)                       !total chemodenitrification N2O uptake non-band unconstrained by N2O, [g d-2 h-1]
  real(r8),pointer ::  RVMBC(:,:,:)                       !total chemodenitrification N2O uptake band unconstrained by N2O, [g d-2 h-1]
  real(r8),pointer ::  XCODFR(:,:)                        !soil surface CO2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),pointer ::  XCHDFR(:,:)                        !soil surface CH4 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),pointer ::  XOXDFR(:,:)                        !soil surface O2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),pointer ::  XNGDFR(:,:)                        !soil surface N2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),pointer ::  XN2DFR(:,:)                        !soil surface N2O dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),pointer ::  XN3DFR(:,:)                        !soil surface NH3 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),pointer ::  XHGDFR(:,:)                        !soil surface H2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),pointer ::  XCOBBL(:,:,:)                      !CO2 bubbling, [g d-2 h-1]
  real(r8),pointer ::  XCHBBL(:,:,:)                      !CH4 bubbling, [g d-2 h-1]
  real(r8),pointer ::  XOXBBL(:,:,:)                      !O2 bubbling, [g d-2 h-1]
  real(r8),pointer ::  XNGBBL(:,:,:)                      !N2 bubbling, [g d-2 h-1]
  real(r8),pointer ::  XN2BBL(:,:,:)                      !N2O bubbling, [g d-2 h-1]
  real(r8),pointer ::  XN3BBL(:,:,:)                      !NH3 bubbling non-band, [g d-2 h-1]
  real(r8),pointer ::  XNBBBL(:,:,:)                      !NH3 bubbling band, [g d-2 h-1]
  real(r8),pointer ::  XHGBBL(:,:,:)                      !H2 bubbling, [g d-2 h-1]
  real(r8),pointer ::  XZHYS(:,:,:)                       !total H+ production
  real(r8),pointer ::  FLW(:,:,:,:)                       !water flux micropore, [m3 d-2 h-1]
  real(r8),pointer ::  FLWH(:,:,:,:)                      !water flux macropore, [m3 d-2 h-1]
  real(r8),pointer ::  HFLW(:,:,:,:)                      !convective heat flux micropore, [MJ d-2 h-1]
  real(r8),pointer ::  XCOFLS(:,:,:,:)                    !aqueous CO2 flux micropore, [g d-2 h-1]
  real(r8),pointer ::  XCHFLS(:,:,:,:)                    !aqueous CH4 flux micropore, [g d-2 h-1]
  real(r8),pointer ::  XOXFLS(:,:,:,:)                    !aqueous O2 flux micropore, [g d-2 h-1]
  real(r8),pointer ::  XNGFLS(:,:,:,:)                    !aqueous N2 flux micropore, [g d-2 h-1]
  real(r8),pointer ::  XN2FLS(:,:,:,:)                    !aqueous N2O flux micropore, [g d-2 h-1]
  real(r8),pointer ::  XHGFLS(:,:,:,:)                    !aqueous H2 flux micropore, [g d-2 h-1]
  real(r8),pointer ::  XN4FLW(:,:,:,:)                    !aqueous NH4 flux non-band micropore, [g d-2 h-1]
  real(r8),pointer ::  XN3FLW(:,:,:,:)                    !aqueous NH3 flux non-band micropore, [g d-2 h-1]
  real(r8),pointer ::  XNOFLW(:,:,:,:)                    !aqueous NO3 flux non-band micropore, [g d-2 h-1]
  real(r8),pointer ::  XH2PFS(:,:,:,:)                    !aqueous PO4 flux non-band micropore, [g d-2 h-1]
  real(r8),pointer ::  XNXFLS(:,:,:,:)                    !aqueous NO2 flux non-band micropore, [g d-2 h-1]
  real(r8),pointer ::  XN4FLB(:,:,:,:)                    !aqueous NH4 flux band micropore, [g d-2 h-1]
  real(r8),pointer ::  XN3FLB(:,:,:,:)                    !aqueous NH3 flux band micropore, [g d-2 h-1]
  real(r8),pointer ::  XNOFLB(:,:,:,:)                    !aqueous NO3 flux band micropore, [g d-2 h-1]
  real(r8),pointer ::  XH2BFB(:,:,:,:)                    !aqueous PO4 flux band micropore, [g d-2 h-1]
  real(r8),pointer ::  XNXFLB(:,:,:,:)                    !aqueous NO2 flux band micropore, [g d-2 h-1]
  real(r8),pointer ::  XOCFLS(:,:,:,:,:)                  !DOC flux micropore, [g d-2 h-1]
  real(r8),pointer ::  XONFLS(:,:,:,:,:)                  !DON flux micropore, [g d-2 h-1]

  real(r8),pointer ::  XOAFLS(:,:,:,:,:)                  !aqueous acetate flux, [g d-2 h-1]

  real(r8),pointer ::  R3GasADTFlx(:,:,:,:,:)             !3D gaseous fluxes, [g d-2 h-1]
  real(r8),pointer ::  XOPFLS(:,:,:,:,:)                  !DOP flux micropore, [g d-2 h-1]
  real(r8),pointer ::  XH1PFS(:,:,:,:)                    !total HPO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8),pointer ::  XH1BFB(:,:,:,:)                    !total HPO4 in micropore water flux band, [mol d-2 h-1]
  real(r8),pointer ::  XCOFHS(:,:,:,:)                    !aqueous CO2 flux macropore, [g d-2 h-1]
  real(r8),pointer ::  XCHFHS(:,:,:,:)                    !aqueous CH4 flux macropore, [g d-2 h-1]
  real(r8),pointer ::  XOXFHS(:,:,:,:)                    !aqueous O2 flux macropore, [g d-2 h-1]
  real(r8),pointer ::  XNGFHS(:,:,:,:)                    !aqueous N2 flux macropore, [g d-2 h-1]
  real(r8),pointer ::  XN2FHS(:,:,:,:)                    !aqueous N2O flux macropore, [g d-2 h-1]
  real(r8),pointer ::  XHGFHS(:,:,:,:)                    !aqueous H2 flux macropore, [g d-2 h-1]
  real(r8),pointer ::  XN4FHW(:,:,:,:)                    !aqueous NH4 flux non-band macropore, [g d-2 h-1]
  real(r8),pointer ::  XN3FHW(:,:,:,:)                    !aqueous NH3 flux non-band macropore, [g d-2 h-1]
  real(r8),pointer ::  XNOFHW(:,:,:,:)                    !aqueous NO3 flux non-band macropore, [g d-2 h-1]
  real(r8),pointer ::  XH2PHS(:,:,:,:)                    !aqueous PO4 flux non-band macropore, [g d-2 h-1]
  real(r8),pointer ::  XNXFHS(:,:,:,:)                    !aqueous NO2 flux non-band macropore, [g d-2 h-1]
  real(r8),pointer ::  XN4FHB(:,:,:,:)                    !aqueous NH4 flux band macropore, [g d-2 h-1]
  real(r8),pointer ::  XN3FHB(:,:,:,:)                    !aqueous NH3 flux band macropore, [g d-2 h-1]
  real(r8),pointer ::  XNOFHB(:,:,:,:)                    !aqueous NO3 flux band macropore, [g d-2 h-1]
  real(r8),pointer ::  XNXFHB(:,:,:,:)                    !aqueous PO4 flux band macropore, [g d-2 h-1]
  real(r8),pointer ::  XH2BHB(:,:,:,:)                    !aqueous NO2 flux band macropore, [g d-2 h-1]
  real(r8),pointer ::  XOCFHS(:,:,:,:,:)                  !DOC flux macropore, [g d-2 h-1]
  real(r8),pointer ::  XONFHS(:,:,:,:,:)                  !DON flux macropore, [g d-2 h-1]
  real(r8),pointer ::  XOPFHS(:,:,:,:,:)                  !DOP flux macropore, [g d-2 h-1]
  real(r8),pointer ::  XOAFHS(:,:,:,:,:)                  !acetate flux macropore, [g d-2 h-1]
  real(r8),pointer ::  XH1PHS(:,:,:,:)                    !total HPO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8),pointer ::  XH1BHB(:,:,:,:)                    !total HPO4 in macropore water flux band, [mol d-2 h-1]

  private :: InitAllocate
  contains

  subroutine InitSoilBGCData(n_pltlitrk)

  implicit none
  integer, intent(in) :: n_pltlitrk

  call InitAllocate(n_pltlitrk)
  end subroutine InitSoilBGCData

!------------------------------------------------------------------------------------------

  subroutine InitAllocate(n_pltlitrk)
  implicit none
  integer, intent(in) :: n_pltlitrk

  allocate(CNH4(JZ,JY,JX));     CNH4=0._r8
  allocate(CNO3(JZ,JY,JX));     CNO3=0._r8
  allocate(CPO4(JZ,JY,JX));     CPO4=0._r8

  allocate(trc_gasml(idg_beg:idg_end,JZ,JY,JX)); trc_gasml=0._r8
  allocate(trc_soHml(ids_beg:ids_end,0:JZ,JY,JX)); trc_soHml=0._r8
  allocate(trc_solml(ids_beg:ids_end,0:JZ,JY,JX)); trc_solml=0._r8
  allocate(trc_solcl(ids_beg:ids_end,0:JZ,JY,JX)); trc_solcl=0._r8
  allocate(trc_gascl(idg_beg:idg_end,0:JZ,JY,JX)); trc_gascl=0._r8

  allocate(ZNFNI(0:JZ,JY,JX));  ZNFNI=0._r8
  allocate(ZNFN0(0:JZ,JY,JX));  ZNFN0=0._r8
  allocate(ZNHUI(0:JZ,JY,JX));  ZNHUI=0._r8
  allocate(ZNHU0(0:JZ,JY,JX));  ZNHU0=0._r8
  allocate(CPO4B(0:JZ,JY,JX));CPO4B(0:JZ,JY,JX)=0._r8
  allocate(CNH3G(0:JZ,JY,JX));CNH3G(0:JZ,JY,JX)=0._r8
  allocate(CZ2GG(0:JZ,JY,JX));CZ2GG(0:JZ,JY,JX)=0._r8
  allocate(CZ2OG(0:JZ,JY,JX));CZ2OG(0:JZ,JY,JX)=0._r8


  allocate(COXYG(0:JZ,JY,JX));COXYG(0:JZ,JY,JX)=0._r8
  allocate(CCH4G(0:JZ,JY,JX));CCH4G(0:JZ,JY,JX)=0._r8

  allocate(CCO2G(0:JZ,JY,JX));CCO2G(0:JZ,JY,JX)=0._r8

  allocate(CH1P4(0:JZ,JY,JX));CH1P4(0:JZ,JY,JX)=0._r8

  allocate(CH2GG(0:JZ,JY,JX));CH2GG(0:JZ,JY,JX)=0._r8
  allocate(PH(0:JZ,JY,JX));PH(0:JZ,JY,JX)=0._r8
  allocate(CEC(JZ,JY,JX));CEC(JZ,JY,JX)=0._r8
  allocate(AEC(JZ,JY,JX));AEC(JZ,JY,JX)=0._r8

  allocate(ROXSK(60,0:JZ,JY,JX));ROXSK=0._r8
  allocate(HCO2G(JY,JX));       HCO2G=0._r8
  allocate(HCH4G(JY,JX));       HCH4G=0._r8
  allocate(HOXYG(JY,JX));       HOXYG=0._r8
  allocate(HN2OG(JY,JX));       HN2OG=0._r8
  allocate(HNH3G(JY,JX));       HNH3G=0._r8
  allocate(UORGF(JY,JX));       UORGF=0._r8
  allocate(UFERTN(JY,JX));      UFERTN=0._r8
  allocate(UFERTP(JY,JX));      UFERTP=0._r8
  allocate(UDOCQ(JY,JX));       UDOCQ=0._r8
  allocate(UDOCD(JY,JX));       UDOCD=0._r8
  allocate(UXCSN(JY,JX));       UXCSN=0._r8
  allocate(UXZSN(JY,JX));       UXZSN=0._r8
  allocate(UXPSN(JY,JX));       UXPSN=0._r8
  allocate(UDONQ(JY,JX));       UDONQ=0._r8
  allocate(UDOND(JY,JX));       UDOND=0._r8
  allocate(UDOPQ(JY,JX));       UDOPQ=0._r8
  allocate(UDOPD(JY,JX));       UDOPD=0._r8
  allocate(UPP4(JY,JX));        UPP4=0._r8
  allocate(UN2GS(JY,JX));       UN2GS=0._r8
  allocate(UH2GG(JY,JX));       UH2GG=0._r8
  allocate(HN2GG(JY,JX));       HN2GG=0._r8
  allocate(UN2GG(JY,JX));       UN2GG=0._r8
  allocate(UCO2G(JY,JX));       UCO2G=0._r8
  allocate(UCH4G(JY,JX));       UCH4G=0._r8
  allocate(UOXYG(JY,JX));       UOXYG=0._r8
  allocate(UNH3G(JY,JX));       UNH3G=0._r8
  allocate(UN2OG(JY,JX));       UN2OG=0._r8
  allocate(UCOP(JY,JX));        UCOP=0._r8
  allocate(USEDOU(JY,JX));      USEDOU=0._r8
  allocate(UDICQ(JY,JX));       UDICQ=0._r8
  allocate(UDICD(JY,JX));       UDICD=0._r8
  allocate(UDINQ(JY,JX));       UDINQ=0._r8
  allocate(UDIND(JY,JX));       UDIND=0._r8
  allocate(UDIPQ(JY,JX));       UDIPQ=0._r8
  allocate(UDIPD(JY,JX));       UDIPD=0._r8
  allocate(WTSTGT(JY,JX));      WTSTGT=0._r8
  allocate(ZDRAIN(JY,JX));      ZDRAIN=0._r8
  allocate(PDRAIN(JY,JX));      PDRAIN=0._r8
  allocate(UION(JY,JX));        UION=0._r8
  allocate(UIONOU(JY,JX));      UIONOU=0._r8
  allocate(XNO2S(0:JZ,JY,JX));  XNO2S=0._r8
  allocate(RUPOXO(0:JZ,JY,JX)); RUPOXO=0._r8
  allocate(RCO2O(0:JZ,JY,JX));  RCO2O=0._r8
  allocate(RCH4O(0:JZ,JY,JX));  RCH4O=0._r8
  allocate(RH2GO(0:JZ,JY,JX));  RH2GO=0._r8
  allocate(RN2G(0:JZ,JY,JX));   RN2G=0._r8
  allocate(RN2O(0:JZ,JY,JX));   RN2O=0._r8
  allocate(XNO2B(0:JZ,JY,JX));  XNO2B=0._r8
  allocate(XNH4B(0:JZ,JY,JX));  XNH4B=0._r8
  allocate(XNO3B(0:JZ,JY,JX));  XNO3B=0._r8
  allocate(XH2BS(0:JZ,JY,JX));  XH2BS=0._r8
  allocate(XH1BS(0:JZ,JY,JX));  XH1BS=0._r8
  allocate(XN2GS(0:JZ,JY,JX));  XN2GS=0._r8
  allocate(XH1PS(0:JZ,JY,JX));  XH1PS=0._r8
  allocate(XH2PS(0:JZ,JY,JX));  XH2PS=0._r8
  allocate(XNH4S(0:JZ,JY,JX));  XNH4S=0._r8
  allocate(XNO3S(0:JZ,JY,JX));  XNO3S=0._r8
  allocate(XOQCS(1:jcplx,0:JZ,JY,JX));XOQCS=0._r8
  allocate(XOQNS(1:jcplx,0:JZ,JY,JX));XOQNS=0._r8
  allocate(XOQPS(1:jcplx,0:JZ,JY,JX));XOQPS=0._r8
  allocate(XOQAS(1:jcplx,0:JZ,JY,JX));XOQAS=0._r8
  allocate(TOQCK(0:JZ,JY,JX));  TOQCK=0._r8
  allocate(VOLQ(0:JZ,JY,JX));   VOLQ=0._r8
  allocate(TFNQ(0:JZ,JY,JX));   TFNQ=0._r8
  allocate(ESNT(jsken,npelms,1:n_pltlitrk,0:JZ,JY,JX));ESNT=0._r8
  allocate(VLNH4(0:JZ,JY,JX));  VLNH4=0._r8
  allocate(VLNHB(0:JZ,JY,JX));  VLNHB=0._r8
  allocate(VLNO3(0:JZ,JY,JX));  VLNO3=0._r8
  allocate(VLNOB(0:JZ,JY,JX));  VLNOB=0._r8
  allocate(VLPO4(0:JZ,JY,JX));  VLPO4=0._r8
  allocate(VLPOB(0:JZ,JY,JX));  VLPOB=0._r8
  allocate(WDNHB(JZ,JY,JX));    WDNHB=0._r8
  allocate(DPNHB(JZ,JY,JX));    DPNHB=0._r8
  allocate(WDNOB(JZ,JY,JX));    WDNOB=0._r8
  allocate(DPNOB(JZ,JY,JX));    DPNOB=0._r8
  allocate(WDPOB(JZ,JY,JX));    WDPOB=0._r8
  allocate(DPPOB(JZ,JY,JX));    DPPOB=0._r8
  allocate(DPNH4(JY,JX));       DPNH4=0._r8
  allocate(DPNO3(JY,JX));       DPNO3=0._r8
  allocate(DPPO4(JY,JX));       DPPO4=0._r8
  allocate(RVMXC(0:JZ,JY,JX));  RVMXC=0._r8
  allocate(RVMBC(0:JZ,JY,JX));  RVMBC=0._r8
  allocate(XCODFR(JY,JX));      XCODFR=0._r8
  allocate(XCHDFR(JY,JX));      XCHDFR=0._r8
  allocate(XOXDFR(JY,JX));      XOXDFR=0._r8
  allocate(XNGDFR(JY,JX));      XNGDFR=0._r8
  allocate(XN2DFR(JY,JX));      XN2DFR=0._r8
  allocate(XN3DFR(JY,JX));      XN3DFR=0._r8
  allocate(XHGDFR(JY,JX));      XHGDFR=0._r8
  allocate(XCOBBL(JZ,JY,JX));   XCOBBL=0._r8
  allocate(XCHBBL(JZ,JY,JX));   XCHBBL=0._r8
  allocate(XOXBBL(JZ,JY,JX));   XOXBBL=0._r8
  allocate(XNGBBL(JZ,JY,JX));   XNGBBL=0._r8
  allocate(XN2BBL(JZ,JY,JX));   XN2BBL=0._r8
  allocate(XN3BBL(JZ,JY,JX));   XN3BBL=0._r8
  allocate(XNBBBL(JZ,JY,JX));   XNBBBL=0._r8
  allocate(XHGBBL(JZ,JY,JX));   XHGBBL=0._r8
  allocate(XZHYS(0:JZ,JY,JX));  XZHYS=0._r8
  allocate(FLW(3,JD,JV,JH));    FLW=0._r8
  allocate(FLWH(3,JD,JV,JH));   FLWH=0._r8
  allocate(HFLW(3,JD,JV,JH));   HFLW=0._r8
  allocate(XCOFLS(3,0:JD,JV,JH));XCOFLS=0._r8
  allocate(XCHFLS(3,0:JD,JV,JH));XCHFLS=0._r8
  allocate(XOXFLS(3,0:JD,JV,JH));XOXFLS=0._r8
  allocate(XNGFLS(3,0:JD,JV,JH));XNGFLS=0._r8
  allocate(XN2FLS(3,0:JD,JV,JH));XN2FLS=0._r8
  allocate(XHGFLS(3,0:JD,JV,JH));XHGFLS=0._r8
  allocate(XN4FLW(3,0:JD,JV,JH));XN4FLW=0._r8
  allocate(XN3FLW(3,0:JD,JV,JH));XN3FLW=0._r8
  allocate(XNOFLW(3,0:JD,JV,JH));XNOFLW=0._r8
  allocate(XH2PFS(3,0:JD,JV,JH));XH2PFS=0._r8
  allocate(XNXFLS(3,0:JD,JV,JH));XNXFLS=0._r8
  allocate(XN4FLB(3,JD,JV,JH)); XN4FLB=0._r8
  allocate(XN3FLB(3,JD,JV,JH)); XN3FLB=0._r8
  allocate(XNOFLB(3,JD,JV,JH)); XNOFLB=0._r8
  allocate(XH2BFB(3,JD,JV,JH)); XH2BFB=0._r8
  allocate(XNXFLB(3,JD,JV,JH)); XNXFLB=0._r8
  allocate(XOCFLS(1:jcplx,3,0:JD,JV,JH));XOCFLS=0._r8
  allocate(XONFLS(1:jcplx,3,0:JD,JV,JH));XONFLS=0._r8
  allocate(XOAFLS(1:jcplx,3,0:JD,JV,JH));XOAFLS=0._r8
  allocate(R3GasADTFlx(idg_beg:idg_end,3,JD,JV,JH));R3GasADTFlx=0._r8
  allocate(XOPFLS(1:jcplx,3,0:JD,JV,JH));XOPFLS=0._r8
  allocate(CPO4S(JZ,JY,JX));CPO4S(JZ,JY,JX)=0._r8
  allocate(XH1PFS(3,0:JD,JV,JH));XH1PFS=0._r8
  allocate(XH1BFB(3,0:JD,JV,JH));XH1BFB=0._r8
  allocate(XCOFHS(3,JD,JV,JH)); XCOFHS=0._r8
  allocate(XCHFHS(3,JD,JV,JH)); XCHFHS=0._r8
  allocate(XOXFHS(3,JD,JV,JH)); XOXFHS=0._r8
  allocate(XNGFHS(3,JD,JV,JH)); XNGFHS=0._r8
  allocate(XN2FHS(3,JD,JV,JH)); XN2FHS=0._r8
  allocate(XHGFHS(3,JD,JV,JH)); XHGFHS=0._r8
  allocate(XN4FHW(3,JD,JV,JH)); XN4FHW=0._r8
  allocate(XN3FHW(3,JD,JV,JH)); XN3FHW=0._r8
  allocate(XNOFHW(3,JD,JV,JH)); XNOFHW=0._r8
  allocate(XH2PHS(3,JD,JV,JH)); XH2PHS=0._r8
  allocate(XNXFHS(3,JD,JV,JH)); XNXFHS=0._r8
  allocate(XN4FHB(3,JD,JV,JH)); XN4FHB=0._r8
  allocate(XN3FHB(3,JD,JV,JH)); XN3FHB=0._r8
  allocate(XNOFHB(3,JD,JV,JH)); XNOFHB=0._r8
  allocate(XNXFHB(3,0:JD,JV,JH));XNXFHB=0._r8
  allocate(XH2BHB(3,JD,JV,JH)); XH2BHB=0._r8
  allocate(XOCFHS(1:jcplx,3,JD,JV,JH));XOCFHS=0._r8
  allocate(XONFHS(1:jcplx,3,JD,JV,JH));XONFHS=0._r8
  allocate(XOPFHS(1:jcplx,3,JD,JV,JH));XOPFHS=0._r8
  allocate(XOAFHS(1:jcplx,3,JD,JV,JH));XOAFHS=0._r8
  allocate(XH1PHS(3,JD,JV,JH)); XH1PHS=0._r8
  allocate(XH1BHB(3,JD,JV,JH)); XH1BHB=0._r8

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructSoilBGCData
  use abortutils, only : destroy

  implicit none
  call destroy(CNH4)
  call destroy(CNO3)
  call destroy(CPO4)

  call destroy(trc_gasml)
  call destroy(CPO4B)

  call destroy(trc_solml)
  call destroy(trc_soHml)

  call destroy(ZNFNI)
  call destroy(ZNFN0)
  call destroy(ZNHUI)
  call destroy(ZNHU0)
  call destroy(CNH3G)
  call destroy(CZ2GG)
  call destroy(CZ2OG)

  call destroy(COXYG)
  call destroy(CCH4G)

  call destroy(CCO2G)

  call destroy(CH1P4)


  call destroy(CH2GG)
  call destroy(PH)
  call destroy(CEC)
  call destroy(AEC)
  call destroy(CPO4S)

  call destroy(ROXSK)
  call destroy(HCO2G)
  call destroy(HCH4G)
  call destroy(HOXYG)
  call destroy(HN2OG)
  call destroy(HNH3G)
  call destroy(UORGF)
  call destroy(UFERTN)
  call destroy(UFERTP)
  call destroy(UDOCQ)
  call destroy(UDOCD)
  call destroy(UXCSN)
  call destroy(UXZSN)
  call destroy(UXPSN)
  call destroy(UDONQ)
  call destroy(UDOND)
  call destroy(UDOPQ)
  call destroy(UDOPD)
  call destroy(UPP4)
  call destroy(UN2GS)
  call destroy(UH2GG)
  call destroy(HN2GG)
  call destroy(UN2GG)
  call destroy(UCO2G)
  call destroy(UCH4G)
  call destroy(UOXYG)
  call destroy(UNH3G)
  call destroy(UN2OG)
  call destroy(UCOP)
  call destroy(USEDOU)
  call destroy(UDICQ)
  call destroy(UDICD)
  call destroy(UDINQ)
  call destroy(UDIND)
  call destroy(UDIPQ)
  call destroy(UDIPD)
  call destroy(WTSTGT)
  call destroy(ZDRAIN)
  call destroy(PDRAIN)
  call destroy(UION)
  call destroy(UIONOU)
  call destroy(XNO2S)
  call destroy(RUPOXO)
  call destroy(RCO2O)
  call destroy(RCH4O)
  call destroy(RH2GO)
  call destroy(RN2G)
  call destroy(RN2O)
  call destroy(XNO2B)
  call destroy(XNH4B)
  call destroy(XNO3B)
  call destroy(XH2BS)
  call destroy(XH1BS)
  call destroy(XN2GS)
  call destroy(XH1PS)
  call destroy(XH2PS)
  call destroy(XNH4S)
  call destroy(XNO3S)
  call destroy(XOQCS)
  call destroy(XOQNS)
  call destroy(XOQPS)
  call destroy(XOQAS)
  call destroy(TOQCK)
  call destroy(VOLQ)
  call destroy(TFNQ)
  call destroy(ESNT)
  call destroy(VLNH4)
  call destroy(VLNHB)
  call destroy(VLNO3)
  call destroy(VLNOB)
  call destroy(VLPO4)
  call destroy(VLPOB)
  call destroy(WDNHB)
  call destroy(DPNHB)
  call destroy(WDNOB)
  call destroy(DPNOB)
  call destroy(WDPOB)
  call destroy(DPPOB)
  call destroy(DPNH4)
  call destroy(DPNO3)
  call destroy(DPPO4)
  call destroy(RVMXC)
  call destroy(RVMBC)
  call destroy(XCODFR)
  call destroy(XCHDFR)
  call destroy(XOXDFR)
  call destroy(XNGDFR)
  call destroy(XN2DFR)
  call destroy(XN3DFR)
  call destroy(XHGDFR)
  call destroy(XCOBBL)
  call destroy(XCHBBL)
  call destroy(XOXBBL)
  call destroy(XNGBBL)
  call destroy(XN2BBL)
  call destroy(XN3BBL)
  call destroy(XNBBBL)
  call destroy(XHGBBL)
  call destroy(XZHYS)
  call destroy(FLW)
  call destroy(FLWH)
  call destroy(HFLW)
  call destroy(XCOFLS)
  call destroy(XCHFLS)
  call destroy(XOXFLS)
  call destroy(XNGFLS)
  call destroy(XN2FLS)
  call destroy(XHGFLS)
  call destroy(XN4FLW)
  call destroy(XN3FLW)
  call destroy(XNOFLW)
  call destroy(XH2PFS)
  call destroy(XNXFLS)
  call destroy(XN4FLB)
  call destroy(XN3FLB)
  call destroy(XNOFLB)
  call destroy(XH2BFB)
  call destroy(XNXFLB)
  call destroy(XOCFLS)
  call destroy(XONFLS)

  call destroy(XOAFLS)
  call destroy(XOPFLS)
  call destroy(XH1PFS)
  call destroy(XH1BFB)
  call destroy(XCOFHS)
  call destroy(XCHFHS)
  call destroy(XOXFHS)
  call destroy(XNGFHS)
  call destroy(XN2FHS)
  call destroy(XHGFHS)
  call destroy(XN4FHW)
  call destroy(XN3FHW)
  call destroy(XNOFHW)
  call destroy(XH2PHS)
  call destroy(XNXFHS)
  call destroy(XN4FHB)
  call destroy(XN3FHB)
  call destroy(XNOFHB)
  call destroy(XNXFHB)
  call destroy(XH2BHB)
  call destroy(XOCFHS)
  call destroy(XONFHS)
  call destroy(XOPFHS)
  call destroy(XOAFHS)
  call destroy(XH1PHS)
  call destroy(XH1BHB)

  end subroutine DestructSoilBGCData

end module SoilBGCDataType
