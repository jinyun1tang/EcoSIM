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

  real(r8) :: H1PO4(0:JZ,JY,JX)                 !soil aqueous HPO4 content micropore non-band, [mol d-2]
  real(r8) :: H1POB(0:JZ,JY,JX)                 !soil aqueous HPO4 content micropore band, [mol d-2]
  real(r8) :: H1PO4H(JZ,JY,JX)                  !soil aqueous HPO4 content non-band macropore, [mol d-2]
  real(r8) :: H1POBH(JZ,JY,JX)                  !soil aqueous HPO4 content band macropore, [mol d-2]

  real(r8) :: ROXSK(60,0:JZ,JY,JX)              !total O2 sink, [g d-2 t-1]
  real(r8) :: HCO2G(JY,JX)                      !soil CO2 flux, [g d-2 h-1]
  real(r8) :: HCH4G(JY,JX)                      !soil CH4 flux, [g d-2 h-1]
  real(r8) :: HOXYG(JY,JX)                      !soil O2 flux, [g d-2 h-1]
  real(r8) :: HN2OG(JY,JX)                      !soil N2O flux, [g d-2 h-1]
  real(r8) :: HNH3G(JY,JX)                      !soil NH3 flux, [g d-2 h-1]
  real(r8) :: UORGF(JY,JX)                      !total C amendment, [g d-2]
  real(r8) :: UFERTN(JY,JX)                     !total fertilizer N amendment, [g d-2]
  real(r8) :: UFERTP(JY,JX)                     !total fertilizer P amendment, [g d-2]
  real(r8) :: UDOCQ(JY,JX)                      !total surface DOC flux, [g d-2]
  real(r8) :: UDOCD(JY,JX)                      !total subsurface DOC flux, [g d-2]
  real(r8) :: UXCSN(JY,JX)                      !total litterfall C, [g d-2]
  real(r8) :: UXZSN(JY,JX)                      !total litterfall N, [g d-2]
  real(r8) :: UXPSN(JY,JX)                      !total litterfall P, [g d-2]
  real(r8) :: UDONQ(JY,JX)                      !total surface DON flux, [g d-2]
  real(r8) :: UDOND(JY,JX)                      !total subsurface DON flux, [g d-2]
  real(r8) :: UDOPQ(JY,JX)                      !total surface DOP flux, [g d-2]
  real(r8) :: UDOPD(JY,JX)                      !total subsurface DOP flux, [g d-2]
  real(r8) :: UPP4(JY,JX)                       !total soil precipited P, [g d-2]
  real(r8) :: UN2GS(JY,JX)                      !total N2 fixation, [g d-2]
  real(r8) :: UH2GG(JY,JX)                      !total H2 flux, []
  real(r8) :: HN2GG(JY,JX)                      !soil N2 flux, [g d-2 h-1]
  real(r8) :: UN2GG(JY,JX)                      !total soil N2 flux, [g d-2]
  real(r8) :: UCO2G(JY,JX)                      !total soil CO2 flux, [g d-2]
  real(r8) :: UCH4G(JY,JX)                      !total soil CH4 flux, [g d-2]
  real(r8) :: UOXYG(JY,JX)                      !total soil O2 flux, [g d-2]
  real(r8) :: UNH3G(JY,JX)                      !total soil NH3 flux, [g d-2]
  real(r8) :: UN2OG(JY,JX)                      !total soil N2O flux, [g d-2]
  real(r8) :: UCOP(JY,JX)                       !total soil autotrophic respiration, [g d-2]
  real(r8) :: USEDOU(JY,JX)                     !total sediment subsurface flux, [Mg d-2]
  real(r8) :: UDICQ(JY,JX)                      !total surface DIC flux, [g d-2]
  real(r8) :: UDICD(JY,JX)                      !total subsurface DIC flux, [g d-2]
  real(r8) :: UDINQ(JY,JX)                      !total surface DIN flux, [g d-2]
  real(r8) :: UDIND(JY,JX)                      !total subsurface DIN flux, [g d-2]
  real(r8) :: UDIPQ(JY,JX)                      !total surface DIP flux, [g d-2]
  real(r8) :: UDIPD(JY,JX)                      !total subsurface DIP flux, [g d-2]
  real(r8) :: WTSTGT(JY,JX)                     !total standing dead C, [g d-2]

  real(r8) :: ZDRAIN(JY,JX)                     !total N drainage below root zone, [g d-2]
  real(r8) :: PDRAIN(JY,JX)                     !total P drainage below root zone, [g d-2]
  real(r8) :: UION(JY,JX)                       !total soil ion content, [mol d-2]
  real(r8) :: UIONOU(JY,JX)                     !total subsurface ion flux, [mol d-2]

! biological and chemical rates

  real(r8) :: XNO2S(0:JZ,JY,JX)                 !total NO2 exchange, [g d-2 h-1]
  real(r8) :: RUPOXO(0:JZ,JY,JX)                !microbial O2 uptake, [g d-2 h-1]
  real(r8) :: RCO2O(0:JZ,JY,JX)                 !microbial net CO2 exchange, [g d-2 h-1]
  real(r8) :: RCH4O(0:JZ,JY,JX)                 !microbial net CH4 exchange, [g d-2 h-1]
  real(r8) :: RH2GO(0:JZ,JY,JX)                 !microbial net H2 exchange, [g d-2 h-1]
  real(r8) :: RN2G(0:JZ,JY,JX)                  !microbial net N2 exchange, [g d-2 h-1]
  real(r8) :: RN2O(0:JZ,JY,JX)                  !microbial net N2O exchange, [g d-2 h-1]
  real(r8) :: XNO2B(0:JZ,JY,JX)                 !net microbial NO2 exchange band, [g d-2 h-1]
  real(r8) :: XNH4B(0:JZ,JY,JX)                 !net microbial NH4 exchange band, [g d-2 h-1]
  real(r8) :: XNO3B(0:JZ,JY,JX)                 !net microbial NO3 exchange band, [g d-2 h-1]
  real(r8) :: XH2BS(0:JZ,JY,JX)                 !net microbial PO4 exchange band, [g d-2 h-1]
  real(r8) :: XH1BS(0:JZ,JY,JX)                 !net microbial HPO4 exchange band, [g d-2 h-1]
  real(r8) :: XN2GS(0:JZ,JY,JX)                 !net microbial N2 exchange, [g d-2 h-1]
  real(r8) :: XH1PS(0:JZ,JY,JX)                 !net microibal HPO4 exchange non-band, [g d-2 h-1]
  real(r8) :: XH2PS(0:JZ,JY,JX)                 !net microbial PO4 exchange nonband, [g d-2 h-1]
  real(r8) :: XNH4S(0:JZ,JY,JX)                 !net microbial NH4 exchange non-band, [g d-2 h-1]
  real(r8) :: XNO3S(0:JZ,JY,JX)                 !net microbial NO3 exchange non-band, [g d-2 h-1]
  real(r8) :: XOQCS(0:jcplx1,0:JZ,JY,JX)             !net microbial DOC flux, [g d-2 h-1]
  real(r8) :: XOQNS(0:jcplx1,0:JZ,JY,JX)             !net microbial DON flux, [g d-2 h-1]
  real(r8) :: XOQPS(0:jcplx1,0:JZ,JY,JX)             !net microbial DOP flux, [g d-2 h-1]
  real(r8) :: XOQAS(0:jcplx1,0:JZ,JY,JX)             !net microbial acetate flux, [g d-2 h-1]
  real(r8) :: TOQCK(0:JZ,JY,JX)                 !total respiration of DOC+DOA in soil layer
  real(r8) :: VOLQ(0:JZ,JY,JX)                  !soil water volume occupied by microial biomass, [m3 m-3]
  real(r8) :: TFNQ(0:JZ,JY,JX)                  !constraints of temperature and water potential on microbial activity, []
  real(r8) :: CSNT(4,0:1,0:JZ,JY,JX)            !total litterfall C, [g d-2 h-1]
  real(r8) :: ZSNT(4,0:1,0:JZ,JY,JX)            !total litterfall N, [g d-2 h-1]
  real(r8) :: PSNT(4,0:1,0:JZ,JY,JX)            !total litterfall P, [g d-2 h-1]

  real(r8) :: VLNH4(0:JZ,JY,JX)                 !NH4 non-band volume fracrion, []
  real(r8) :: VLNHB(0:JZ,JY,JX)                 !NH4 band volume fracrion, []
  real(r8) :: VLNO3(0:JZ,JY,JX)                 !NO3 non-band volume fracrion, []
  real(r8) :: VLNOB(0:JZ,JY,JX)                 !NO3 band volume fracrion, []
  real(r8) :: VLPO4(0:JZ,JY,JX)                 !PO4 non-band volume fracrion, []
  real(r8) :: VLPOB(0:JZ,JY,JX)                 !PO4 band volume fracrion, []
  real(r8) :: WDNHB(JZ,JY,JX)                   !width of NH4 band, [m]
  real(r8) :: DPNHB(JZ,JY,JX)                   !depth of NH4 band, [m]
  real(r8) :: WDNOB(JZ,JY,JX)                   !width of NO3 band, [m]
  real(r8) :: DPNOB(JZ,JY,JX)                   !depth of NO4 band, [m]
  real(r8) :: WDPOB(JZ,JY,JX)                   !width of PO4 band, [m]
  real(r8) :: DPPOB(JZ,JY,JX)                   !depth of PO4 band, [m]
  real(r8) :: DPNH4(JY,JX)                      !total depth of NH4 band, [m]
  real(r8) :: DPNO3(JY,JX)                      !total depth of NO3 band, [m]
  real(r8) :: DPPO4(JY,JX)                      !total depth of PO4 band, [m]
  real(r8) :: RVMXC(0:JZ,JY,JX)                 !total chemodenitrification N2O uptake non-band unconstrained by N2O, [g d-2 h-1]
  real(r8) :: RVMBC(0:JZ,JY,JX)                 !total chemodenitrification N2O uptake band unconstrained by N2O, [g d-2 h-1]
  real(r8) :: XCODFR(JY,JX)                     !soil surface CO2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8) :: XCHDFR(JY,JX)                     !soil surface CH4 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8) :: XOXDFR(JY,JX)                     !soil surface O2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8) :: XNGDFR(JY,JX)                     !soil surface N2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8) :: XN2DFR(JY,JX)                     !soil surface N2O dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8) :: XN3DFR(JY,JX)                     !soil surface NH3 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8) :: XHGDFR(JY,JX)                     !soil surface H2 dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8) :: XCOBBL(JZ,JY,JX)                  !CO2 bubbling, [g d-2 h-1]
  real(r8) :: XCHBBL(JZ,JY,JX)                  !CH4 bubbling, [g d-2 h-1]
  real(r8) :: XOXBBL(JZ,JY,JX)                  !O2 bubbling, [g d-2 h-1]
  real(r8) :: XNGBBL(JZ,JY,JX)                  !N2 bubbling, [g d-2 h-1]
  real(r8) :: XN2BBL(JZ,JY,JX)                  !N2O bubbling, [g d-2 h-1]
  real(r8) :: XN3BBL(JZ,JY,JX)                  !NH3 bubbling non-band, [g d-2 h-1]
  real(r8) :: XNBBBL(JZ,JY,JX)                  !NH3 bubbling band, [g d-2 h-1]
  real(r8) :: XHGBBL(JZ,JY,JX)                  !H2 bubbling, [g d-2 h-1]
  real(r8) :: XZHYS(0:JZ,JY,JX)                 !total H+ production
  real(r8) :: FLW(3,JD,JV,JH)                   !water flux micropore, [m3 d-2 h-1]
  real(r8) :: FLWH(3,JD,JV,JH)                  !water flux macropore, [m3 d-2 h-1]
  real(r8) :: HFLW(3,JD,JV,JH)                  !convective heat flux micropore, [MJ d-2 h-1]
  real(r8) :: XCOFLS(3,0:JD,JV,JH)              !aqueous CO2 flux micropore, [g d-2 h-1]
  real(r8) :: XCHFLS(3,0:JD,JV,JH)              !aqueous CH4 flux micropore, [g d-2 h-1]
  real(r8) :: XOXFLS(3,0:JD,JV,JH)              !aqueous O2 flux micropore, [g d-2 h-1]
  real(r8) :: XNGFLS(3,0:JD,JV,JH)              !aqueous N2 flux micropore, [g d-2 h-1]
  real(r8) :: XN2FLS(3,0:JD,JV,JH)              !aqueous N2O flux micropore, [g d-2 h-1]
  real(r8) :: XHGFLS(3,0:JD,JV,JH)              !aqueous H2 flux micropore, [g d-2 h-1]
  real(r8) :: XN4FLW(3,0:JD,JV,JH)              !aqueous NH4 flux non-band micropore, [g d-2 h-1]
  real(r8) :: XN3FLW(3,0:JD,JV,JH)              !aqueous NH3 flux non-band micropore, [g d-2 h-1]
  real(r8) :: XNOFLW(3,0:JD,JV,JH)              !aqueous NO3 flux non-band micropore, [g d-2 h-1]
  real(r8) :: XH2PFS(3,0:JD,JV,JH)              !aqueous PO4 flux non-band micropore, [g d-2 h-1]
  real(r8) :: XNXFLS(3,0:JD,JV,JH)              !aqueous NO2 flux non-band micropore, [g d-2 h-1]
  real(r8) :: XN4FLB(3,JD,JV,JH)                !aqueous NH4 flux band micropore, [g d-2 h-1]
  real(r8) :: XN3FLB(3,JD,JV,JH)                !aqueous NH3 flux band micropore, [g d-2 h-1]
  real(r8) :: XNOFLB(3,JD,JV,JH)                !aqueous NO3 flux band micropore, [g d-2 h-1]
  real(r8) :: XH2BFB(3,JD,JV,JH)                !aqueous PO4 flux band micropore, [g d-2 h-1]
  real(r8) :: XNXFLB(3,JD,JV,JH)                !aqueous NO2 flux band micropore, [g d-2 h-1]
  real(r8) :: XOCFLS(0:jcplx1,3,0:JD,JV,JH)          !DOC flux micropore, [g d-2 h-1]
  real(r8) :: XONFLS(0:jcplx1,3,0:JD,JV,JH)          !DON flux micropore, [g d-2 h-1]
  real(r8) :: XCHFLG(3,0:JD,JV,JH)              !gaseous CH4 flux, [g d-2 h-1]
  real(r8) :: XOAFLS(0:jcplx1,3,0:JD,JV,JH)          !aqueous acetate flux, [g d-2 h-1]
  real(r8) :: XCOFLG(3,JD,JV,JH)                !gaseous CO2 flux, [g d-2 h-1]
  real(r8) :: XOPFLS(0:jcplx1,3,0:JD,JV,JH)          !DOP flux micropore, [g d-2 h-1]
  real(r8) :: XOXFLG(3,JD,JV,JH)                !gaseous O24 flux, [g d-2 h-1]
  real(r8) :: XNGFLG(3,JD,JV,JH)                !gaseous N2 flux, [g d-2 h-1]
  real(r8) :: XHGFLG(3,JD,JV,JH)                !gaseous H2 flux, [g d-2 h-1]
  real(r8) :: XN2FLG(3,JD,JV,JH)                !gaseous N2O flux, [g d-2 h-1]
  real(r8) :: XN3FLG(3,JD,JV,JH)                !gaseous NH3 flux, [g d-2 h-1]
  real(r8) :: XH1PFS(3,0:JD,JV,JH)              !total HPO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8) :: XH1BFB(3,0:JD,JV,JH)              !total HPO4 in micropore water flux band, [mol d-2 h-1]
  real(r8) :: XCOFHS(3,JD,JV,JH)                !aqueous CO2 flux macropore, [g d-2 h-1]
  real(r8) :: XCHFHS(3,JD,JV,JH)                !aqueous CH4 flux macropore, [g d-2 h-1]
  real(r8) :: XOXFHS(3,JD,JV,JH)                !aqueous O2 flux macropore, [g d-2 h-1]
  real(r8) :: XNGFHS(3,JD,JV,JH)                !aqueous N2 flux macropore, [g d-2 h-1]
  real(r8) :: XN2FHS(3,JD,JV,JH)                !aqueous N2O flux macropore, [g d-2 h-1]
  real(r8) :: XHGFHS(3,JD,JV,JH)                !aqueous H2 flux macropore, [g d-2 h-1]
  real(r8) :: XN4FHW(3,JD,JV,JH)                !aqueous NH4 flux non-band macropore, [g d-2 h-1]
  real(r8) :: XN3FHW(3,JD,JV,JH)                !aqueous NH3 flux non-band macropore, [g d-2 h-1]
  real(r8) :: XNOFHW(3,JD,JV,JH)                !aqueous NO3 flux non-band macropore, [g d-2 h-1]
  real(r8) :: XH2PHS(3,JD,JV,JH)                !aqueous PO4 flux non-band macropore, [g d-2 h-1]
  real(r8) :: XNXFHS(3,JD,JV,JH)                !aqueous NO2 flux non-band macropore, [g d-2 h-1]
  real(r8) :: XN4FHB(3,JD,JV,JH)                !aqueous NH4 flux band macropore, [g d-2 h-1]
  real(r8) :: XN3FHB(3,JD,JV,JH)                !aqueous NH3 flux band macropore, [g d-2 h-1]
  real(r8) :: XNOFHB(3,JD,JV,JH)                !aqueous NO3 flux band macropore, [g d-2 h-1]
  real(r8) :: XNXFHB(3,0:JD,JV,JH)              !aqueous PO4 flux band macropore, [g d-2 h-1]
  real(r8) :: XH2BHB(3,JD,JV,JH)                !aqueous NO2 flux band macropore, [g d-2 h-1]
  real(r8) :: XOCFHS(0:jcplx1,3,JD,JV,JH)            !DOC flux macropore, [g d-2 h-1]
  real(r8) :: XONFHS(0:jcplx1,3,JD,JV,JH)            !DON flux macropore, [g d-2 h-1]
  real(r8) :: XOPFHS(0:jcplx1,3,JD,JV,JH)            !DOP flux macropore, [g d-2 h-1]
  real(r8) :: XOAFHS(0:jcplx1,3,JD,JV,JH)            !acetate flux macropore, [g d-2 h-1]
  real(r8) :: XH1PHS(3,JD,JV,JH)                !total HPO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8) :: XH1BHB(3,JD,JV,JH)                !total HPO4 in macropore water flux band, [mol d-2 h-1]

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
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructSoilBGCData

  implicit none
  if (allocated(CNH4))     deallocate(CNH4)
  if (allocated(CNO3))     deallocate(CNO3)
  if (allocated(CPO4))     deallocate(CPO4)
  if (allocated(H2PO4))    deallocate(H2PO4)
  if (allocated(ZNH4B))    deallocate(ZNH4B)
  if (allocated(ZNH3B))    deallocate(ZNH3B)
  if (allocated(ZNO3B))    deallocate(ZNO3B)
  if (allocated(H2POB))    deallocate(H2POB)
  if (allocated(ZNO2S))    deallocate(ZNO2S)
  if (allocated(ZNH3G))    deallocate(ZNH3G)
  if (allocated(Z2GG))     deallocate(Z2GG)
  if (allocated(Z2GS))     deallocate(Z2GS)
  if (allocated(Z2OG))     deallocate(Z2OG)
  if (allocated(Z2OS))     deallocate(Z2OS)
  if (allocated(ZNH4SH))   deallocate(ZNH4SH)
  if (allocated(ZNH3SH))   deallocate(ZNH3SH)
  if (allocated(ZNO3SH))   deallocate(ZNO3SH)
  if (allocated(H2PO4H))   deallocate(H2PO4H)
  if (allocated(ZNH4BH))   deallocate(ZNH4BH)
  if (allocated(ZNH3BH))   deallocate(ZNH3BH)
  if (allocated(ZNO3BH))   deallocate(ZNO3BH)
  if (allocated(H2POBH))   deallocate(H2POBH)
  if (allocated(ZNO2SH))   deallocate(ZNO2SH)
  if (allocated(Z2GSH))    deallocate(Z2GSH)
  if (allocated(Z2OSH))    deallocate(Z2OSH)
  if (allocated(ZNO2BH))   deallocate(ZNO2BH)
  if (allocated(ZNO2B))    deallocate(ZNO2B)
  if (allocated(ZNH4S))    deallocate(ZNH4S)
  if (allocated(ZNH3S))    deallocate(ZNH3S)
  if (allocated(ZNO3S))    deallocate(ZNO3S)
  if (allocated(ZNFNI))    deallocate(ZNFNI)
  if (allocated(ZNFN0))    deallocate(ZNFN0)
  if (allocated(ZNHUI))    deallocate(ZNHUI)
  if (allocated(ZNHU0))    deallocate(ZNHU0)
  if(allocated(CZ2GS))deallocate(CZ2GS)
  if(allocated(CNH4S))deallocate(CNH4S)
  if(allocated(CNH3S))deallocate(CNH3S)
  if(allocated(CNO3S))deallocate(CNO3S)
  if(allocated(CPO4S))deallocate(CPO4S)
  if(allocated(CNH4B))deallocate(CNH4B)
  if(allocated(CNH3B))deallocate(CNH3B)
  if(allocated(CNO3B))deallocate(CNO3B)
  if(allocated(CPO4B))deallocate(CPO4B)
  if(allocated(CNO2S))deallocate(CNO2S)
  if(allocated(CNH3G))deallocate(CNH3G)
  if(allocated(CZ2GG))deallocate(CZ2GG)
  if(allocated(CZ2OG))deallocate(CZ2OG)
  if(allocated(CZ2OS))deallocate(CZ2OS)
  if(allocated(OXYG))deallocate(OXYG)
  if(allocated(OXYS))deallocate(OXYS)
  if(allocated(OXYSH))deallocate(OXYSH)
  if(allocated(CO2G))deallocate(CO2G)
  if(allocated(CO2S))deallocate(CO2S)
  if(allocated(CO2SH))deallocate(CO2SH)
  if(allocated(CH4G))deallocate(CH4G)
  if(allocated(CH4S))deallocate(CH4S)
  if(allocated(CH4SH))deallocate(CH4SH)
  if(allocated(COXYG))deallocate(COXYG)
  if(allocated(CCH4G))deallocate(CCH4G)
  if(allocated(COXYS))deallocate(COXYS)
  if(allocated(CCO2G))deallocate(CCO2G)
  if(allocated(CCO2S))deallocate(CCO2S)
  if(allocated(CCH4S))deallocate(CCH4S)
  if(allocated(CH1P4))deallocate(CH1P4)
  if(allocated(CH1P4B))deallocate(CH1P4B)
  if(allocated(CNO2B))deallocate(CNO2B)
  if(allocated(H2GS))deallocate(H2GS)
  if(allocated(CH2GS))deallocate(CH2GS)
  if(allocated(CH2P4))deallocate(CH2P4)
  if(allocated(CH2P4B))deallocate(CH2P4B)
  if(allocated(H2GSH))deallocate(H2GSH)
  if(allocated(H2GG))deallocate(H2GG)
  if(allocated(CH2GG))deallocate(CH2GG)
  if(allocated(PH))deallocate(PH)
  if(allocated(CEC))deallocate(CEC)
  if(allocated(AEC))deallocate(AEC)
  end subroutine DestructSoilBGCData

end module SoilBGCDataType
