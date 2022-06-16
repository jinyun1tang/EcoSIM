module SoilBGCDataType

!
! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable ::  CNH4(:,:,:)                       !soil NH4 content, [mg kg-1]
  real(r8),allocatable ::  CNO3(:,:,:)                       !soil NO3 content, [mg kg-1]
  real(r8),allocatable ::  CPO4(:,:,:)                       !soil PO4 content, [mg kg-1]
  real(r8),allocatable ::  H2PO4(:,:,:)                      !PO4 non-band micropore, [g d-2]
  real(r8),allocatable ::  ZNH4B(:,:,:)                      !NH4 band micropore, [g d-2]
  real(r8),allocatable ::  ZNH3B(:,:,:)                      !NH3 band micropore, [g d-2]
  real(r8),allocatable ::  ZNO3B(:,:,:)                      !NO3 band micropore, [g d-2]
  real(r8),allocatable ::  H2POB(:,:,:)                      !PO4 band micropore, [g d-2]
  real(r8),allocatable ::  ZNO2S(:,:,:)                      !NO2  non-band micropore, [g d-2]
  real(r8),allocatable ::  ZNH3G(:,:,:)                      !gaseous NH3, [g d-2]
  real(r8),allocatable ::  Z2GG(:,:,:)                       !gaseous N2, [g d-2]
  real(r8),allocatable ::  Z2GS(:,:,:)                       !aqueous N2 micropore, [g d-2]
  real(r8),allocatable ::  Z2OG(:,:,:)                       !gaseous N2O, [g d-2]
  real(r8),allocatable ::  Z2OS(:,:,:)                       !aqueous N2O micropore, [g d-2]
  real(r8),allocatable ::  ZNH4SH(:,:,:)                     !NH4 non-band macropore, [g d-2]
  real(r8),allocatable ::  ZNH3SH(:,:,:)                     !NH3 non-band macropore, [g d-2]
  real(r8),allocatable ::  ZNO3SH(:,:,:)                     !NO3 non-band macropore, [g d-2]
  real(r8),allocatable ::  H2PO4H(:,:,:)                     !PO4 non-band macropore, [g d-2]
  real(r8),allocatable ::  ZNH4BH(:,:,:)                     !NH4 band macropore, [g d-2]
  real(r8),allocatable ::  ZNH3BH(:,:,:)                     !NH3 band macropore, [g d-2]
  real(r8),allocatable ::  ZNO3BH(:,:,:)                     !NO3 band macropore, [g d-2]
  real(r8),allocatable ::  H2POBH(:,:,:)                     !PO4 band macropore, [g d-2]
  real(r8),allocatable ::  ZNO2SH(:,:,:)                     !NO2  non-band macropore, [g d-2]
  real(r8),allocatable ::  Z2GSH(:,:,:)                      !aqueous N2 macropore, [g d-2]
  real(r8),allocatable ::  Z2OSH(:,:,:)                      !aqueous N2O macropore, [g d-2]
  real(r8),allocatable ::  ZNO2BH(:,:,:)                     !NO2 band macropore, [g d-2]
  real(r8),allocatable ::  ZNO2B(:,:,:)                      !NO2  band micropore, [g d-2]
  real(r8),allocatable ::  ZNH4S(:,:,:)                      !NH4 non-band micropore, [g d-2]
  real(r8),allocatable ::  ZNH3S(:,:,:)                      !NH3 non-band micropore, [g d-2]
  real(r8),allocatable ::  ZNO3S(:,:,:)                      !NO3 non-band micropore, [g d-2]
  real(r8),allocatable ::  ZNFNI(:,:,:)                      !current nitrification inhibition activity
  real(r8),allocatable ::  ZNFN0(:,:,:)                      !initial nitrification inhibition activity
  real(r8),allocatable ::  ZNHUI(:,:,:)                      !current inhibition activity
  real(r8),allocatable ::  ZNHU0(:,:,:)                      !urea hydrolysis inhibition activity
  real(r8),allocatable :: CZ2GS(:,:,:)          !aqueous N2 concentration micropore	[g m-3]
  real(r8),allocatable :: CNH4S(:,:,:)          !NH4 concentration non-band micropore	[g m-3]
  real(r8),allocatable :: CNH3S(:,:,:)          !NH3 concentration non-band micropore	[g m-3]
  real(r8),allocatable :: CNO3S(:,:,:)          !NO3 concentration non-band micropore	[g m-3]
  real(r8),allocatable :: CPO4S(:,:,:)          !PO4 concentration non-band micropore	[g m-3]
  real(r8),allocatable :: CNH4B(:,:,:)          !NH4 concentration band micropore	[g m-3]
  real(r8),allocatable :: CNH3B(:,:,:)          !NH3 concentration band micropore	[g m-3]
  real(r8),allocatable :: CNO3B(:,:,:)          !NO3 concentration band micropore	[g m-3]
  real(r8),allocatable :: CPO4B(:,:,:)          !PO4 concentration band micropore	[g m-3]
  real(r8),allocatable :: CNO2S(:,:,:)          !NO2 concentration non-band micropore	[g m-3]
  real(r8),allocatable :: CNH3G(:,:,:)          !gaseous NH3 concentration	[g m-3]
  real(r8),allocatable :: CZ2GG(:,:,:)          !gaseous N2 concentration	[g m-3]
  real(r8),allocatable :: CZ2OG(:,:,:)          !gaseous N2O concentration	[g m-3]
  real(r8),allocatable :: CZ2OS(:,:,:)          !aqueous N2O concentration micropore	[g m-3]
  real(r8),allocatable :: COXYG(:,:,:)          !gaseous O2 concentration	[g m-3]
  real(r8),allocatable :: CCH4G(:,:,:)          !gaseous CH4 concentration	[g m-3]
  real(r8),allocatable :: COXYS(:,:,:)          !aqueous O2 concentration micropore	[g m-3]
  real(r8),allocatable :: CCO2G(:,:,:)          !gaseous CO2 concentration	[g m-3]
  real(r8),allocatable :: CCO2S(:,:,:)          !aqueous CO2 concentration micropore	[g m-3]
  real(r8),allocatable :: CCH4S(:,:,:)          !aqueous CH4 concentration micropore	[g m-3]
  real(r8),allocatable :: CH1P4(:,:,:)          !aqueous H1PO4 concentration non-band [g m-3]
  real(r8),allocatable :: CH1P4B(:,:,:)         !aqueous H1PO4 concentration band [g m-3]
  real(r8),allocatable :: CNO2B(:,:,:)          !aqueous HNO2 concentration band [g m-3]
  real(r8),allocatable :: CH2GS(:,:,:)          !aqueous H2 concentration	[g m-3]
  real(r8),allocatable :: CH2P4(:,:,:)          !aqueous PO4 concentration non-band	[g m-3]
  real(r8),allocatable :: CH2P4B(:,:,:)         !aqueous PO4 concentration band	[g m-3]
  real(r8),allocatable :: CH2GG(:,:,:)          !gaseous H2 concentration	[g m-3]
  real(r8),allocatable :: OXYS(:,:,:)           !aqueous O2  micropore	[g d-2]
  real(r8),allocatable :: OXYSH(:,:,:)          !aqueous O2 macropore	g [d-2]
  real(r8),allocatable :: OXYG(:,:,:)           !gaseous O2 	[g d-2]
  real(r8),allocatable :: CO2G(:,:,:)           !gaseous CO2	[g d-2]
  real(r8),allocatable :: CO2S(:,:,:)           !aqueous CO2  micropore	[g d-2]
  real(r8),allocatable :: CO2SH(:,:,:)          !aqueous CO2  macropore	[g d-2]
  real(r8),allocatable :: CH4G(:,:,:)           !gaseous CH4	[g d-2]
  real(r8),allocatable :: CH4S(:,:,:)           !aqueous CO2  micropore	[g d-2]
  real(r8),allocatable :: CH4SH(:,:,:)          !aqueous CO2  macropore	[g d-2]
  real(r8),allocatable :: H2GS(:,:,:)           !aqueous H2 	[g d-2]
  real(r8),allocatable :: H2GSH(:,:,:)          !aqueous H2 macropore	[g d-2]
  real(r8),allocatable :: H2GG(:,:,:)           !gaseous H2	[g d-2]
  real(r8),allocatable :: PH(:,:,:)             !soil pH
  real(r8),allocatable :: CEC(:,:,:)            !soil cation exchange capacity	[cmol kg-1]
  real(r8),allocatable :: AEC(:,:,:)            !soil anion exchange capacity	[cmol kg-1]

  real(r8),allocatable ::  H1PO4(:,:,:)                       !soil aqueous HPO4 content micropore non-band, [mol d-2]
  real(r8),allocatable ::  H1POB(:,:,:)                       !soil aqueous HPO4 content micropore band, [mol d-2]
  real(r8),allocatable ::  H1PO4H(:,:,:)                      !soil aqueous HPO4 content non-band macropore, [mol d-2]
  real(r8),allocatable ::  H1POBH(:,:,:)                      !soil aqueous HPO4 content band macropore, [mol d-2]
  real(r8),allocatable ::  ROXSK(:,:,:,:)                     !total O2 sink, [g d-2 t-1]
  real(r8),allocatable ::  HCO2G(:,:)                         !soil CO2 flux, [g d-2 h-1]
  real(r8),allocatable ::  HCH4G(:,:)                         !soil CH4 flux, [g d-2 h-1]
  real(r8),allocatable ::  HOXYG(:,:)                         !soil O2 flux, [g d-2 h-1]
  real(r8),allocatable ::  HN2OG(:,:)                         !soil N2O flux, [g d-2 h-1]
  real(r8),allocatable ::  HNH3G(:,:)                         !soil NH3 flux, [g d-2 h-1]
  real(r8),allocatable ::  UORGF(:,:)                         !total C amendment, [g d-2]
  real(r8),allocatable ::  UFERTN(:,:)                        !total fertilizer N amendment, [g d-2]
  real(r8),allocatable ::  UFERTP(:,:)                        !total fertilizer P amendment, [g d-2]
  real(r8),allocatable ::  UDOCQ(:,:)                         !total surface DOC flux, [g d-2]
  real(r8),allocatable ::  UDOCD(:,:)                         !total subsurface DOC flux, [g d-2]
  real(r8),allocatable ::  UXCSN(:,:)                         !total litterfall C, [g d-2]
  real(r8),allocatable ::  UXZSN(:,:)                         !total litterfall N, [g d-2]
  real(r8),allocatable ::  UXPSN(:,:)                         !total litterfall P, [g d-2]
  real(r8),allocatable ::  UDONQ(:,:)                         !total surface DON flux, [g d-2]
  real(r8),allocatable ::  UDOND(:,:)                         !total subsurface DON flux, [g d-2]
  real(r8),allocatable ::  UDOPQ(:,:)                         !total surface DOP flux, [g d-2]
  real(r8),allocatable ::  UDOPD(:,:)                         !total subsurface DOP flux, [g d-2]
  real(r8),allocatable ::  UPP4(:,:)                          !total soil precipited P, [g d-2]
  real(r8),allocatable ::  UN2GS(:,:)                         !total N2 fixation, [g d-2]
  real(r8),allocatable ::  UH2GG(:,:)                         !total H2 flux, []
  real(r8),allocatable ::  HN2GG(:,:)                         !soil N2 flux, [g d-2 h-1]
  real(r8),allocatable ::  UN2GG(:,:)                         !total soil N2 flux, [g d-2]
  real(r8),allocatable ::  UCO2G(:,:)                         !total soil CO2 flux, [g d-2]
  real(r8),allocatable ::  UCH4G(:,:)                         !total soil CH4 flux, [g d-2]
  real(r8),allocatable ::  UOXYG(:,:)                         !total soil O2 flux, [g d-2]
  real(r8),allocatable ::  UNH3G(:,:)                         !total soil NH3 flux, [g d-2]
  real(r8),allocatable ::  UN2OG(:,:)                         !total soil N2O flux, [g d-2]
  real(r8),allocatable ::  UCOP(:,:)                          !total soil autotrophic respiration, [g d-2]
  real(r8),allocatable ::  USEDOU(:,:)                        !total sediment subsurface flux, [Mg d-2]
  real(r8),allocatable ::  UDICQ(:,:)                         !total surface DIC flux, [g d-2]
  real(r8),allocatable ::  UDICD(:,:)                         !total subsurface DIC flux, [g d-2]
  real(r8),allocatable ::  UDINQ(:,:)                         !total surface DIN flux, [g d-2]
  real(r8),allocatable ::  UDIND(:,:)                         !total subsurface DIN flux, [g d-2]
  real(r8),allocatable ::  UDIPQ(:,:)                         !total surface DIP flux, [g d-2]
  real(r8),allocatable ::  UDIPD(:,:)                         !total subsurface DIP flux, [g d-2]
  real(r8),allocatable ::  WTSTGT(:,:)                        !total standing dead C, [g d-2]
  real(r8),allocatable ::  ZDRAIN(:,:)                        !total N drainage below root zone, [g d-2]
  real(r8),allocatable ::  PDRAIN(:,:)                        !total P drainage below root zone, [g d-2]
  real(r8),allocatable ::  UION(:,:)                          !total soil ion content, [mol d-2]
  real(r8),allocatable ::  UIONOU(:,:)                        !total subsurface ion flux, [mol d-2]
  real(r8),allocatable ::  XNO2S(:,:,:)                       !total NO2 exchange, [g d-2 h-1]
  real(r8),allocatable ::  RUPOXO(:,:,:)                      !microbial O2 uptake, [g d-2 h-1]
  real(r8),allocatable ::  RCO2O(:,:,:)                       !microbial net CO2 exchange, [g d-2 h-1]
  real(r8),allocatable ::  RCH4O(:,:,:)                       !microbial net CH4 exchange, [g d-2 h-1]
  real(r8),allocatable ::  RH2GO(:,:,:)                       !microbial net H2 exchange, [g d-2 h-1]
  real(r8),allocatable ::  RN2G(:,:,:)                        !microbial net N2 exchange, [g d-2 h-1]
  real(r8),allocatable ::  RN2O(:,:,:)                        !microbial net N2O exchange, [g d-2 h-1]
  real(r8),allocatable ::  XNO2B(:,:,:)                       !net microbial NO2 exchange band, [g d-2 h-1]
  real(r8),allocatable ::  XNH4B(:,:,:)                       !net microbial NH4 exchange band, [g d-2 h-1]
  real(r8),allocatable ::  XNO3B(:,:,:)                       !net microbial NO3 exchange band, [g d-2 h-1]
  real(r8),allocatable ::  XH2BS(:,:,:)                       !net microbial PO4 exchange band, [g d-2 h-1]
  real(r8),allocatable ::  XH1BS(:,:,:)                       !net microbial HPO4 exchange band, [g d-2 h-1]
  real(r8),allocatable ::  XN2GS(:,:,:)                       !net microbial N2 exchange, [g d-2 h-1]
  real(r8),allocatable ::  XH1PS(:,:,:)                       !net microibal HPO4 exchange non-band, [g d-2 h-1]
  real(r8),allocatable ::  XH2PS(:,:,:)                       !net microbial PO4 exchange nonband, [g d-2 h-1]
  real(r8),allocatable ::  XNH4S(:,:,:)                       !net microbial NH4 exchange non-band, [g d-2 h-1]
  real(r8),allocatable ::  XNO3S(:,:,:)                       !net microbial NO3 exchange non-band, [g d-2 h-1]
  real(r8),allocatable ::  XOQCS(:,:,:,:)                     !net microbial DOC flux, [g d-2 h-1]
  real(r8),allocatable ::  XOQNS(:,:,:,:)                     !net microbial DON flux, [g d-2 h-1]
  real(r8),allocatable ::  XOQPS(:,:,:,:)                     !net microbial DOP flux, [g d-2 h-1]
  real(r8),allocatable ::  XOQAS(:,:,:,:)                     !net microbial acetate flux, [g d-2 h-1]
  real(r8),allocatable ::  TOQCK(:,:,:)                       !total respiration of DOC+DOA in soil layer
  real(r8),allocatable ::  VOLQ(:,:,:)                        !soil water volume occupied by microial biomass, [m3 m-3]
  real(r8),allocatable ::  TFNQ(:,:,:)                        !constraints of temperature and water potential on microbial activity, []
  real(r8),allocatable ::  CSNT(:,:,:,:,:)                    !total litterfall C, [g d-2 h-1]
  real(r8),allocatable ::  ZSNT(:,:,:,:,:)                    !total litterfall N, [g d-2 h-1]
  real(r8),allocatable ::  PSNT(:,:,:,:,:)                    !total litterfall P, [g d-2 h-1]
  real(r8),allocatable ::  VLNH4(:,:,:)                       !NH4 non-band volume fracrion, []
  real(r8),allocatable ::  VLNHB(:,:,:)                       !NH4 band volume fracrion, []
  real(r8),allocatable ::  VLNO3(:,:,:)                       !NO3 non-band volume fracrion, []
  real(r8),allocatable ::  VLNOB(:,:,:)                       !NO3 band volume fracrion, []
  real(r8),allocatable ::  VLPO4(:,:,:)                       !PO4 non-band volume fracrion, []
  real(r8),allocatable ::  VLPOB(:,:,:)                       !PO4 band volume fracrion, []
  real(r8),allocatable ::  WDNHB(:,:,:)                       !width of NH4 band, [m]
  real(r8),allocatable ::  DPNHB(:,:,:)                       !depth of NH4 band, [m]
  real(r8),allocatable ::  WDNOB(:,:,:)                       !width of NO3 band, [m]
  real(r8),allocatable ::  DPNOB(:,:,:)                       !depth of NO4 band, [m]
  real(r8),allocatable ::  WDPOB(:,:,:)                       !width of PO4 band, [m]
  real(r8),allocatable ::  DPPOB(:,:,:)                       !depth of PO4 band, [m]
  real(r8),allocatable ::  DPNH4(:,:)                         !total depth of NH4 band, [m]
  real(r8),allocatable ::  DPNO3(:,:)                         !total depth of NO3 band, [m]
  real(r8),allocatable ::  DPPO4(:,:)                         !total depth of PO4 band, [m]
  real(r8),allocatable ::  RVMXC(:,:,:)                       !total chemodenitrification N2O uptake non-band unconstrained by N2O, [g d-2 h-1]
  real(r8),allocatable ::  RVMBC(:,:,:)                       !total chemodenitrification N2O uptake band unconstrained by N2O, [g d-2 h-1]
  real(r8),allocatable ::  XCODFR(:,:)                        !soil surface CO2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),allocatable ::  XCHDFR(:,:)                        !soil surface CH4 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),allocatable ::  XOXDFR(:,:)                        !soil surface O2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),allocatable ::  XNGDFR(:,:)                        !soil surface N2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),allocatable ::  XN2DFR(:,:)                        !soil surface N2O dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),allocatable ::  XN3DFR(:,:)                        !soil surface NH3 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),allocatable ::  XHGDFR(:,:)                        !soil surface H2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),allocatable ::  XCOBBL(:,:,:)                      !CO2 bubbling, [g d-2 h-1]
  real(r8),allocatable ::  XCHBBL(:,:,:)                      !CH4 bubbling, [g d-2 h-1]
  real(r8),allocatable ::  XOXBBL(:,:,:)                      !O2 bubbling, [g d-2 h-1]
  real(r8),allocatable ::  XNGBBL(:,:,:)                      !N2 bubbling, [g d-2 h-1]
  real(r8),allocatable ::  XN2BBL(:,:,:)                      !N2O bubbling, [g d-2 h-1]
  real(r8),allocatable ::  XN3BBL(:,:,:)                      !NH3 bubbling non-band, [g d-2 h-1]
  real(r8),allocatable ::  XNBBBL(:,:,:)                      !NH3 bubbling band, [g d-2 h-1]
  real(r8),allocatable ::  XHGBBL(:,:,:)                      !H2 bubbling, [g d-2 h-1]
  real(r8),allocatable ::  XZHYS(:,:,:)                       !total H+ production
  real(r8),allocatable ::  FLW(:,:,:,:)                       !water flux micropore, [m3 d-2 h-1]
  real(r8),allocatable ::  FLWH(:,:,:,:)                      !water flux macropore, [m3 d-2 h-1]
  real(r8),allocatable ::  HFLW(:,:,:,:)                      !convective heat flux micropore, [MJ d-2 h-1]
  real(r8),allocatable ::  XCOFLS(:,:,:,:)                    !aqueous CO2 flux micropore, [g d-2 h-1]
  real(r8),allocatable ::  XCHFLS(:,:,:,:)                    !aqueous CH4 flux micropore, [g d-2 h-1]
  real(r8),allocatable ::  XOXFLS(:,:,:,:)                    !aqueous O2 flux micropore, [g d-2 h-1]
  real(r8),allocatable ::  XNGFLS(:,:,:,:)                    !aqueous N2 flux micropore, [g d-2 h-1]
  real(r8),allocatable ::  XN2FLS(:,:,:,:)                    !aqueous N2O flux micropore, [g d-2 h-1]
  real(r8),allocatable ::  XHGFLS(:,:,:,:)                    !aqueous H2 flux micropore, [g d-2 h-1]
  real(r8),allocatable ::  XN4FLW(:,:,:,:)                    !aqueous NH4 flux non-band micropore, [g d-2 h-1]
  real(r8),allocatable ::  XN3FLW(:,:,:,:)                    !aqueous NH3 flux non-band micropore, [g d-2 h-1]
  real(r8),allocatable ::  XNOFLW(:,:,:,:)                    !aqueous NO3 flux non-band micropore, [g d-2 h-1]
  real(r8),allocatable ::  XH2PFS(:,:,:,:)                    !aqueous PO4 flux non-band micropore, [g d-2 h-1]
  real(r8),allocatable ::  XNXFLS(:,:,:,:)                    !aqueous NO2 flux non-band micropore, [g d-2 h-1]
  real(r8),allocatable ::  XN4FLB(:,:,:,:)                    !aqueous NH4 flux band micropore, [g d-2 h-1]
  real(r8),allocatable ::  XN3FLB(:,:,:,:)                    !aqueous NH3 flux band micropore, [g d-2 h-1]
  real(r8),allocatable ::  XNOFLB(:,:,:,:)                    !aqueous NO3 flux band micropore, [g d-2 h-1]
  real(r8),allocatable ::  XH2BFB(:,:,:,:)                    !aqueous PO4 flux band micropore, [g d-2 h-1]
  real(r8),allocatable ::  XNXFLB(:,:,:,:)                    !aqueous NO2 flux band micropore, [g d-2 h-1]
  real(r8),allocatable ::  XOCFLS(:,:,:,:,:)                  !DOC flux micropore, [g d-2 h-1]
  real(r8),allocatable ::  XONFLS(:,:,:,:,:)                  !DON flux micropore, [g d-2 h-1]
  real(r8),allocatable ::  XCHFLG(:,:,:,:)                    !gaseous CH4 flux, [g d-2 h-1]
  real(r8),allocatable ::  XOAFLS(:,:,:,:,:)                  !aqueous acetate flux, [g d-2 h-1]
  real(r8),allocatable ::  XCOFLG(:,:,:,:)                    !gaseous CO2 flux, [g d-2 h-1]
  real(r8),allocatable ::  XOPFLS(:,:,:,:,:)                  !DOP flux micropore, [g d-2 h-1]
  real(r8),allocatable ::  XOXFLG(:,:,:,:)                    !gaseous O24 flux, [g d-2 h-1]
  real(r8),allocatable ::  XNGFLG(:,:,:,:)                    !gaseous N2 flux, [g d-2 h-1]
  real(r8),allocatable ::  XHGFLG(:,:,:,:)                    !gaseous H2 flux, [g d-2 h-1]
  real(r8),allocatable ::  XN2FLG(:,:,:,:)                    !gaseous N2O flux, [g d-2 h-1]
  real(r8),allocatable ::  XN3FLG(:,:,:,:)                    !gaseous NH3 flux, [g d-2 h-1]
  real(r8),allocatable ::  XH1PFS(:,:,:,:)                    !total HPO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XH1BFB(:,:,:,:)                    !total HPO4 in micropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XCOFHS(:,:,:,:)                    !aqueous CO2 flux macropore, [g d-2 h-1]
  real(r8),allocatable ::  XCHFHS(:,:,:,:)                    !aqueous CH4 flux macropore, [g d-2 h-1]
  real(r8),allocatable ::  XOXFHS(:,:,:,:)                    !aqueous O2 flux macropore, [g d-2 h-1]
  real(r8),allocatable ::  XNGFHS(:,:,:,:)                    !aqueous N2 flux macropore, [g d-2 h-1]
  real(r8),allocatable ::  XN2FHS(:,:,:,:)                    !aqueous N2O flux macropore, [g d-2 h-1]
  real(r8),allocatable ::  XHGFHS(:,:,:,:)                    !aqueous H2 flux macropore, [g d-2 h-1]
  real(r8),allocatable ::  XN4FHW(:,:,:,:)                    !aqueous NH4 flux non-band macropore, [g d-2 h-1]
  real(r8),allocatable ::  XN3FHW(:,:,:,:)                    !aqueous NH3 flux non-band macropore, [g d-2 h-1]
  real(r8),allocatable ::  XNOFHW(:,:,:,:)                    !aqueous NO3 flux non-band macropore, [g d-2 h-1]
  real(r8),allocatable ::  XH2PHS(:,:,:,:)                    !aqueous PO4 flux non-band macropore, [g d-2 h-1]
  real(r8),allocatable ::  XNXFHS(:,:,:,:)                    !aqueous NO2 flux non-band macropore, [g d-2 h-1]
  real(r8),allocatable ::  XN4FHB(:,:,:,:)                    !aqueous NH4 flux band macropore, [g d-2 h-1]
  real(r8),allocatable ::  XN3FHB(:,:,:,:)                    !aqueous NH3 flux band macropore, [g d-2 h-1]
  real(r8),allocatable ::  XNOFHB(:,:,:,:)                    !aqueous NO3 flux band macropore, [g d-2 h-1]
  real(r8),allocatable ::  XNXFHB(:,:,:,:)                    !aqueous PO4 flux band macropore, [g d-2 h-1]
  real(r8),allocatable ::  XH2BHB(:,:,:,:)                    !aqueous NO2 flux band macropore, [g d-2 h-1]
  real(r8),allocatable ::  XOCFHS(:,:,:,:,:)                  !DOC flux macropore, [g d-2 h-1]
  real(r8),allocatable ::  XONFHS(:,:,:,:,:)                  !DON flux macropore, [g d-2 h-1]
  real(r8),allocatable ::  XOPFHS(:,:,:,:,:)                  !DOP flux macropore, [g d-2 h-1]
  real(r8),allocatable ::  XOAFHS(:,:,:,:,:)                  !acetate flux macropore, [g d-2 h-1]
  real(r8),allocatable ::  XH1PHS(:,:,:,:)                    !total HPO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XH1BHB(:,:,:,:)                    !total HPO4 in macropore water flux band, [mol d-2 h-1]

  private :: InitAllocate
  contains

  subroutine InitSoilBGCData

  implicit none

  call InitAllocate
  end subroutine InitSoilBGCData

!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none
  allocate(CNH4(JZ,JY,JX));     CNH4=0._r8
  allocate(CNO3(JZ,JY,JX));     CNO3=0._r8
  allocate(CPO4(JZ,JY,JX));     CPO4=0._r8
  allocate(H2PO4(0:JZ,JY,JX));  H2PO4=0._r8
  allocate(ZNH4B(0:JZ,JY,JX));  ZNH4B=0._r8
  allocate(ZNH3B(0:JZ,JY,JX));  ZNH3B=0._r8
  allocate(ZNO3B(0:JZ,JY,JX));  ZNO3B=0._r8
  allocate(H2POB(0:JZ,JY,JX));  H2POB=0._r8
  allocate(ZNO2S(0:JZ,JY,JX));  ZNO2S=0._r8
  allocate(ZNH3G(JZ,JY,JX));    ZNH3G=0._r8
  allocate(Z2GG(JZ,JY,JX));     Z2GG=0._r8
  allocate(Z2GS(0:JZ,JY,JX));   Z2GS=0._r8
  allocate(Z2OG(JZ,JY,JX));     Z2OG=0._r8
  allocate(Z2OS(0:JZ,JY,JX));   Z2OS=0._r8
  allocate(ZNH4SH(JZ,JY,JX));   ZNH4SH=0._r8
  allocate(ZNH3SH(JZ,JY,JX));   ZNH3SH=0._r8
  allocate(ZNO3SH(JZ,JY,JX));   ZNO3SH=0._r8
  allocate(H2PO4H(JZ,JY,JX));   H2PO4H=0._r8
  allocate(ZNH4BH(JZ,JY,JX));   ZNH4BH=0._r8
  allocate(ZNH3BH(JZ,JY,JX));   ZNH3BH=0._r8
  allocate(ZNO3BH(JZ,JY,JX));   ZNO3BH=0._r8
  allocate(H2POBH(JZ,JY,JX));   H2POBH=0._r8
  allocate(ZNO2SH(JZ,JY,JX));   ZNO2SH=0._r8
  allocate(Z2GSH(JZ,JY,JX));    Z2GSH=0._r8
  allocate(Z2OSH(JZ,JY,JX));    Z2OSH=0._r8
  allocate(ZNO2BH(JZ,JY,JX));   ZNO2BH=0._r8
  allocate(ZNO2B(0:JZ,JY,JX));  ZNO2B=0._r8
  allocate(ZNH4S(0:JZ,JY,JX));  ZNH4S=0._r8
  allocate(ZNH3S(0:JZ,JY,JX));  ZNH3S=0._r8
  allocate(ZNO3S(0:JZ,JY,JX));  ZNO3S=0._r8
  allocate(ZNFNI(0:JZ,JY,JX));  ZNFNI=0._r8
  allocate(ZNFN0(0:JZ,JY,JX));  ZNFN0=0._r8
  allocate(ZNHUI(0:JZ,JY,JX));  ZNHUI=0._r8
  allocate(ZNHU0(0:JZ,JY,JX));  ZNHU0=0._r8
  allocate(CZ2GS(0:JZ,JY,JX));CZ2GS(0:JZ,JY,JX)=0._r8
  allocate(CNH4S(0:JZ,JY,JX));CNH4S(0:JZ,JY,JX)=0._r8
  allocate(CNH3S(0:JZ,JY,JX));CNH3S(0:JZ,JY,JX)=0._r8
  allocate(CNO3S(0:JZ,JY,JX));CNO3S(0:JZ,JY,JX)=0._r8
  allocate(CPO4S(JZ,JY,JX));CPO4S(JZ,JY,JX)=0._r8
  allocate(CNH4B(0:JZ,JY,JX));CNH4B(0:JZ,JY,JX)=0._r8
  allocate(CNH3B(0:JZ,JY,JX));CNH3B(0:JZ,JY,JX)=0._r8
  allocate(CNO3B(0:JZ,JY,JX));CNO3B(0:JZ,JY,JX)=0._r8
  allocate(CPO4B(0:JZ,JY,JX));CPO4B(0:JZ,JY,JX)=0._r8
  allocate(CNO2S(0:JZ,JY,JX));CNO2S(0:JZ,JY,JX)=0._r8
  allocate(CNH3G(0:JZ,JY,JX));CNH3G(0:JZ,JY,JX)=0._r8
  allocate(CZ2GG(0:JZ,JY,JX));CZ2GG(0:JZ,JY,JX)=0._r8
  allocate(CZ2OG(0:JZ,JY,JX));CZ2OG(0:JZ,JY,JX)=0._r8
  allocate(CZ2OS(0:JZ,JY,JX));CZ2OS(0:JZ,JY,JX)=0._r8
  allocate(OXYG(JZ,JY,JX));OXYG(JZ,JY,JX)=0._r8
  allocate(OXYS(0:JZ,JY,JX));OXYS(0:JZ,JY,JX)=0._r8
  allocate(OXYSH(JZ,JY,JX));OXYSH(JZ,JY,JX)=0._r8
  allocate(CO2G(JZ,JY,JX));CO2G(JZ,JY,JX)=0._r8
  allocate(CO2S(0:JZ,JY,JX));CO2S(0:JZ,JY,JX)=0._r8
  allocate(CO2SH(JZ,JY,JX));CO2SH(JZ,JY,JX)=0._r8
  allocate(CH4G(JZ,JY,JX));CH4G(JZ,JY,JX)=0._r8
  allocate(CH4S(0:JZ,JY,JX));CH4S(0:JZ,JY,JX)=0._r8
  allocate(CH4SH(JZ,JY,JX));CH4SH(JZ,JY,JX)=0._r8
  allocate(COXYG(0:JZ,JY,JX));COXYG(0:JZ,JY,JX)=0._r8
  allocate(CCH4G(0:JZ,JY,JX));CCH4G(0:JZ,JY,JX)=0._r8
  allocate(COXYS(0:JZ,JY,JX));COXYS(0:JZ,JY,JX)=0._r8
  allocate(CCO2G(0:JZ,JY,JX));CCO2G(0:JZ,JY,JX)=0._r8
  allocate(CCO2S(0:JZ,JY,JX));CCO2S(0:JZ,JY,JX)=0._r8
  allocate(CCH4S(0:JZ,JY,JX));CCH4S(0:JZ,JY,JX)=0._r8
  allocate(CH1P4(0:JZ,JY,JX));CH1P4(0:JZ,JY,JX)=0._r8
  allocate(CH1P4B(0:JZ,JY,JX));CH1P4B(0:JZ,JY,JX)=0._r8
  allocate(CNO2B(0:JZ,JY,JX));CNO2B(0:JZ,JY,JX)=0._r8
  allocate(H2GS(0:JZ,JY,JX));H2GS(0:JZ,JY,JX)=0._r8
  allocate(CH2GS(0:JZ,JY,JX));CH2GS(0:JZ,JY,JX)=0._r8
  allocate(CH2P4(0:JZ,JY,JX));CH2P4(0:JZ,JY,JX)=0._r8
  allocate(CH2P4B(0:JZ,JY,JX));CH2P4B(0:JZ,JY,JX)=0._r8
  allocate(H2GSH(JZ,JY,JX));H2GSH(JZ,JY,JX)=0._r8
  allocate(H2GG(JZ,JY,JX));H2GG(JZ,JY,JX)=0._r8
  allocate(CH2GG(0:JZ,JY,JX));CH2GG(0:JZ,JY,JX)=0._r8
  allocate(PH(0:JZ,JY,JX));PH(0:JZ,JY,JX)=0._r8
  allocate(CEC(JZ,JY,JX));CEC(JZ,JY,JX)=0._r8
  allocate(AEC(JZ,JY,JX));AEC(JZ,JY,JX)=0._r8
  allocate(H1PO4(0:JZ,JY,JX));  H1PO4=0._r8
  allocate(H1POB(0:JZ,JY,JX));  H1POB=0._r8
  allocate(H1PO4H(JZ,JY,JX));   H1PO4H=0._r8
  allocate(H1POBH(JZ,JY,JX));   H1POBH=0._r8
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
  allocate(XOQCS(0:jcplx1,0:JZ,JY,JX));XOQCS=0._r8
  allocate(XOQNS(0:jcplx1,0:JZ,JY,JX));XOQNS=0._r8
  allocate(XOQPS(0:jcplx1,0:JZ,JY,JX));XOQPS=0._r8
  allocate(XOQAS(0:jcplx1,0:JZ,JY,JX));XOQAS=0._r8
  allocate(TOQCK(0:JZ,JY,JX));  TOQCK=0._r8
  allocate(VOLQ(0:JZ,JY,JX));   VOLQ=0._r8
  allocate(TFNQ(0:JZ,JY,JX));   TFNQ=0._r8
  allocate(CSNT(4,0:1,0:JZ,JY,JX));CSNT=0._r8
  allocate(ZSNT(4,0:1,0:JZ,JY,JX));ZSNT=0._r8
  allocate(PSNT(4,0:1,0:JZ,JY,JX));PSNT=0._r8
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
  allocate(XOCFLS(0:jcplx1,3,0:JD,JV,JH));XOCFLS=0._r8
  allocate(XONFLS(0:jcplx1,3,0:JD,JV,JH));XONFLS=0._r8
  allocate(XCHFLG(3,0:JD,JV,JH));XCHFLG=0._r8
  allocate(XOAFLS(0:jcplx1,3,0:JD,JV,JH));XOAFLS=0._r8
  allocate(XCOFLG(3,JD,JV,JH)); XCOFLG=0._r8
  allocate(XOPFLS(0:jcplx1,3,0:JD,JV,JH));XOPFLS=0._r8
  allocate(XOXFLG(3,JD,JV,JH)); XOXFLG=0._r8
  allocate(XNGFLG(3,JD,JV,JH)); XNGFLG=0._r8
  allocate(XHGFLG(3,JD,JV,JH)); XHGFLG=0._r8
  allocate(XN2FLG(3,JD,JV,JH)); XN2FLG=0._r8
  allocate(XN3FLG(3,JD,JV,JH)); XN3FLG=0._r8
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
  allocate(XOCFHS(0:jcplx1,3,JD,JV,JH));XOCFHS=0._r8
  allocate(XONFHS(0:jcplx1,3,JD,JV,JH));XONFHS=0._r8
  allocate(XOPFHS(0:jcplx1,3,JD,JV,JH));XOPFHS=0._r8
  allocate(XOAFHS(0:jcplx1,3,JD,JV,JH));XOAFHS=0._r8
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
  call destroy(H2PO4)
  call destroy(ZNH4B)
  call destroy(ZNH3B)
  call destroy(ZNO3B)
  call destroy(H2POB)
  call destroy(ZNO2S)
  call destroy(ZNH3G)
  call destroy(Z2GG)
  call destroy(Z2GS)
  call destroy(Z2OG)
  call destroy(Z2OS)
  call destroy(ZNH4SH)
  call destroy(ZNH3SH)
  call destroy(ZNO3SH)
  call destroy(H2PO4H)
  call destroy(ZNH4BH)
  call destroy(ZNH3BH)
  call destroy(ZNO3BH)
  call destroy(H2POBH)
  call destroy(ZNO2SH)
  call destroy(Z2GSH)
  call destroy(Z2OSH)
  call destroy(ZNO2BH)
  call destroy(ZNO2B)
  call destroy(ZNH4S)
  call destroy(ZNH3S)
  call destroy(ZNO3S)
  call destroy(ZNFNI)
  call destroy(ZNFN0)
  call destroy(ZNHUI)
  call destroy(ZNHU0)
  call destroy(CZ2GS)
  call destroy(CNH4S)
  call destroy(CNH3S)
  call destroy(CNO3S)
  call destroy(CPO4S)
  call destroy(CNH4B)
  call destroy(CNH3B)
  call destroy(CNO3B)
  call destroy(CPO4B)
  call destroy(CNO2S)
  call destroy(CNH3G)
  call destroy(CZ2GG)
  call destroy(CZ2OG)
  call destroy(CZ2OS)
  call destroy(OXYG)
  call destroy(OXYS)
  call destroy(OXYSH)
  call destroy(CO2G)
  call destroy(CO2S)
  call destroy(CO2SH)
  call destroy(CH4G)
  call destroy(CH4S)
  call destroy(CH4SH)
  call destroy(COXYG)
  call destroy(CCH4G)
  call destroy(COXYS)
  call destroy(CCO2G)
  call destroy(CCO2S)
  call destroy(CCH4S)
  call destroy(CH1P4)
  call destroy(CH1P4B)
  call destroy(CNO2B)
  call destroy(H2GS)
  call destroy(CH2GS)
  call destroy(CH2P4)
  call destroy(CH2P4B)
  call destroy(H2GSH)
  call destroy(H2GG)
  call destroy(CH2GG)
  call destroy(PH)
  call destroy(CEC)
  call destroy(AEC)

  call destroy(H1PO4)
  call destroy(H1POB)
  call destroy(H1PO4H)
  call destroy(H1POBH)
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
  call destroy(CSNT)
  call destroy(ZSNT)
  call destroy(PSNT)
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
  call destroy(XCHFLG)
  call destroy(XOAFLS)
  call destroy(XCOFLG)
  call destroy(XOPFLS)
  call destroy(XOXFLG)
  call destroy(XNGFLG)
  call destroy(XHGFLG)
  call destroy(XN2FLG)
  call destroy(XN3FLG)
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
