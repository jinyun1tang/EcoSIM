module SoilBGCDataType

!
! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod    , only : NumOfPlantChemElements
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken => jskenc
implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  CNH4(:,:,:)                       !soil NH4 content, [mg kg-1]
  real(r8),target,allocatable ::  CNO3(:,:,:)                       !soil NO3 content, [mg kg-1]
  real(r8),target,allocatable ::  CPO4(:,:,:)                       !soil PO4 content, [mg kg-1]

  real(r8),target,allocatable :: CPO4B(:,:,:)                       !PO4 concentration band micropore	[g m-3]
  real(r8),target,allocatable :: CPO4S(:,:,:)                       !PO4 concentration non-band micropore	[g m-3]

  real(r8),target,allocatable :: trc_soHml(:,:,:,:)                 !solute mass in macropore [g d-2]
  real(r8),target,allocatable :: trc_solml(:,:,:,:)                 !solute mass in micropore [g d-2]
  real(r8),target,allocatable :: trc_solcl(:,:,:,:)                 !solute concentration in micropre [g m-3]
  real(r8),target,allocatable :: trc_gascl(:,:,:,:)                 !gaseous concentation [g m-3]

  real(r8),target,allocatable ::  ZNFNI(:,:,:)                      !current nitrification inhibition activity
  real(r8),target,allocatable ::  ZNFN0(:,:,:)                      !initial nitrification inhibition activity
  real(r8),target,allocatable ::  ZNHUI(:,:,:)                      !current inhibition activity
  real(r8),target,allocatable ::  ZNHU0(:,:,:)                      !urea hydrolysis inhibition activity

  real(r8),target,allocatable :: trc_gasml(:,:,:,:)     !layer mass of gases [g d-2]

  real(r8),target,allocatable :: PH(:,:,:)             !soil pH
  real(r8),target,allocatable :: CEC(:,:,:)            !soil cation exchange capacity	[cmol kg-1]
  real(r8),target,allocatable :: AEC(:,:,:)            !soil anion exchange capacity	[cmol kg-1]

  real(r8),target,allocatable ::  ROXSK(:,:,:,:)                     !total O2 sink, [g d-2 t-1]
  real(r8),target,allocatable ::  HCO2G(:,:)                         !soil CO2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  HCH4G(:,:)                         !soil CH4 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  HOXYG(:,:)                         !soil O2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  HN2OG(:,:)                         !soil N2O flux, [g d-2 h-1]
  real(r8),target,allocatable ::  HNH3G(:,:)                         !soil NH3 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  UORGF(:,:)                         !total C amendment, [g d-2]
  real(r8),target,allocatable ::  UFERTN(:,:)                        !total fertilizer N amendment, [g d-2]
  real(r8),target,allocatable ::  UFERTP(:,:)                        !total fertilizer P amendment, [g d-2]
  real(r8),target,allocatable ::  UDOCQ(:,:)                         !total surface DOC flux, [g d-2]
  real(r8),target,allocatable ::  UDOCD(:,:)                         !total subsurface DOC flux, [g d-2]
  real(r8),target,allocatable ::  UXCSN(:,:)                         !total litterfall C, [g d-2]
  real(r8),target,allocatable ::  UXZSN(:,:)                         !total litterfall N, [g d-2]
  real(r8),target,allocatable ::  UXPSN(:,:)                         !total litterfall P, [g d-2]
  real(r8),target,allocatable ::  UDONQ(:,:)                         !total surface DON flux, [g d-2]
  real(r8),target,allocatable ::  UDOND(:,:)                         !total subsurface DON flux, [g d-2]
  real(r8),target,allocatable ::  UDOPQ(:,:)                         !total surface DOP flux, [g d-2]
  real(r8),target,allocatable ::  UDOPD(:,:)                         !total subsurface DOP flux, [g d-2]
  real(r8),target,allocatable ::  UPP4(:,:)                          !total soil precipited P, [g d-2]
  real(r8),target,allocatable ::  UN2GS(:,:)                         !total N2 fixation, [g d-2]
  real(r8),target,allocatable ::  UH2GG(:,:)                         !total H2 flux, []
  real(r8),target,allocatable ::  HN2GG(:,:)                         !soil N2 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  UN2GG(:,:)                         !total soil N2 flux, [g d-2]
  real(r8),target,allocatable ::  UCO2G(:,:)                         !total soil CO2 flux, [g d-2]
  real(r8),target,allocatable ::  UCH4G(:,:)                         !total soil CH4 flux, [g d-2]
  real(r8),target,allocatable ::  UOXYG(:,:)                         !total soil O2 flux, [g d-2]
  real(r8),target,allocatable ::  UNH3G(:,:)                         !total soil NH3 flux, [g d-2]
  real(r8),target,allocatable ::  UN2OG(:,:)                         !total soil N2O flux, [g d-2]
  real(r8),target,allocatable ::  UCOP(:,:)                          !total soil autotrophic respiration, [g d-2]
  real(r8),target,allocatable ::  USEDOU(:,:)                        !total sediment subsurface flux, [Mg d-2]
  real(r8),target,allocatable ::  UDICQ(:,:)                         !total surface DIC flux, [g d-2]
  real(r8),target,allocatable ::  UDICD(:,:)                         !total subsurface DIC flux, [g d-2]
  real(r8),target,allocatable ::  UDINQ(:,:)                         !total surface DIN flux, [g d-2]
  real(r8),target,allocatable ::  UDIND(:,:)                         !total subsurface DIN flux, [g d-2]
  real(r8),target,allocatable ::  UDIPQ(:,:)                         !total surface DIP flux, [g d-2]
  real(r8),target,allocatable ::  UDIPD(:,:)                         !total subsurface DIP flux, [g d-2]
  real(r8),target,allocatable ::  WTSTGET(:,:,:)                        !total standing dead C, [g d-2]
  real(r8),target,allocatable ::  ZDRAIN(:,:)                        !total N drainage below root zone, [g d-2]
  real(r8),target,allocatable ::  PDRAIN(:,:)                        !total P drainage below root zone, [g d-2]
  real(r8),target,allocatable ::  UION(:,:)                          !total soil ion content, [mol d-2]
  real(r8),target,allocatable ::  UIONOU(:,:)                        !total subsurface ion flux, [mol d-2]
  real(r8),target,allocatable ::  XNO2S(:,:,:)                       !total NO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  RUPOXO(:,:,:)                      !microbial O2 uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RCO2O(:,:,:)                       !microbial net CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  RCH4O(:,:,:)                       !microbial net CH4 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  RH2GO(:,:,:)                       !microbial net H2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  RN2G(:,:,:)                        !microbial net N2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  RN2O(:,:,:)                        !microbial net N2O exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  XNO2B(:,:,:)                       !net microbial NO2 exchange band, [g d-2 h-1]
  real(r8),target,allocatable ::  XNH4B(:,:,:)                       !net microbial NH4 exchange band, [g d-2 h-1]
  real(r8),target,allocatable ::  XNO3B(:,:,:)                       !net microbial NO3 exchange band, [g d-2 h-1]
  real(r8),target,allocatable ::  XH2BS(:,:,:)                       !net microbial PO4 exchange band, [g d-2 h-1]
  real(r8),target,allocatable ::  XH1BS(:,:,:)                       !net microbial HPO4 exchange band, [g d-2 h-1]
  real(r8),target,allocatable ::  XN2GS(:,:,:)                       !net microbial N2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  XH1PS(:,:,:)                       !net microibal HPO4 exchange non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  XH2PS(:,:,:)                       !net microbial PO4 exchange nonband, [g d-2 h-1]
  real(r8),target,allocatable ::  XNH4S(:,:,:)                       !net microbial NH4 exchange non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  XNO3S(:,:,:)                       !net microbial NO3 exchange non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  XOQCS(:,:,:,:)                     !net microbial DOC flux, [g d-2 h-1]
  real(r8),target,allocatable ::  XOQNS(:,:,:,:)                     !net microbial DON flux, [g d-2 h-1]
  real(r8),target,allocatable ::  XOQPS(:,:,:,:)                     !net microbial DOP flux, [g d-2 h-1]
  real(r8),target,allocatable ::  XOQAS(:,:,:,:)                     !net microbial acetate flux, [g d-2 h-1]
  real(r8),target,allocatable ::  TOQCK(:,:,:)                       !total respiration of DOC+DOA in soil layer
  real(r8),target,allocatable ::  VOLQ(:,:,:)                        !soil water volume occupied by microial biomass, [m3 m-3]
  real(r8),target,allocatable ::  TFNQ(:,:,:)                        !constraints of temperature and water potential on microbial activity, []
  real(r8),target,allocatable ::  LitrfalChemElemnts_vr(:,:,:,:,:,:)                    !total litterfall C, [g d-2 h-1]
  real(r8),target,allocatable :: trcs_VLN(:,:,:,:)

  real(r8),target,allocatable ::  VLNHB(:,:,:)                       !NH4 band volume fracrion, []
  real(r8),target,allocatable ::  VLNOB(:,:,:)                       !NO3 band volume fracrion, []
  real(r8),target,allocatable ::  VLPO4(:,:,:)                       !PO4 non-band volume fracrion, []
  real(r8),target,allocatable ::  VLPOB(:,:,:)                       !PO4 band volume fracrion, []
  real(r8),target,allocatable ::  WDNHB(:,:,:)                       !width of NH4 band, [m]
  real(r8),target,allocatable ::  DPNHB(:,:,:)                       !depth of NH4 band, [m]
  real(r8),target,allocatable ::  WDNOB(:,:,:)                       !width of NO3 band, [m]
  real(r8),target,allocatable ::  DPNOB(:,:,:)                       !depth of NO4 band, [m]
  real(r8),target,allocatable ::  WDPOB(:,:,:)                       !width of PO4 band, [m]
  real(r8),target,allocatable ::  DPPOB(:,:,:)                       !depth of PO4 band, [m]
  real(r8),target,allocatable ::  DPNH4(:,:)                         !total depth of NH4 band, [m]
  real(r8),target,allocatable ::  DPNO3(:,:)                         !total depth of NO3 band, [m]
  real(r8),target,allocatable ::  DPPO4(:,:)                         !total depth of PO4 band, [m]
  real(r8),target,allocatable ::  RVMXC(:,:,:)                       !total chemodenitrification N2O uptake non-band unconstrained by N2O, [g d-2 h-1]
  real(r8),target,allocatable ::  RVMBC(:,:,:)                       !total chemodenitrification N2O uptake band unconstrained by N2O, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_XDFR(:,:,:)                   !soil surface gas dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_XBLL(:,:,:,:)                      !CO2 bubbling, [g d-2 h-1]
  real(r8),target,allocatable ::  XZHYS(:,:,:)                       !total H+ production
  real(r8),target,allocatable ::  WaterFlowSoiMicP(:,:,:,:)                       !water flux micropore, [m3 d-2 h-1]
  real(r8),target,allocatable ::  WaterFlowMacP(:,:,:,:)                      !water flux macropore, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HeatFlow2Soil(:,:,:,:)                      !convective heat flux micropore, [MJ d-2 h-1]

  real(r8),target,allocatable ::  trcs_XFLS(:,:,:,:,:)
  real(r8),target,allocatable ::  XOCFLS(:,:,:,:,:)                  !DOC flux micropore, [g d-2 h-1]
  real(r8),target,allocatable ::  XONFLS(:,:,:,:,:)                  !DON flux micropore, [g d-2 h-1]

  real(r8),target,allocatable ::  XOAFLS(:,:,:,:,:)                  !aqueous acetate flux, [g d-2 h-1]

  real(r8),target,allocatable ::  trcs_XFHS(:,:,:,:,:)
  real(r8),target,allocatable ::  R3GasADTFlx(:,:,:,:,:)             !3D gaseous fluxes, [g d-2 h-1]
  real(r8),target,allocatable ::  XOPFLS(:,:,:,:,:)                  !DOP flux micropore, [g d-2 h-1]
  real(r8),target,allocatable ::  XOCFHS(:,:,:,:,:)                  !DOC flux macropore, [g d-2 h-1]
  real(r8),target,allocatable ::  XONFHS(:,:,:,:,:)                  !DON flux macropore, [g d-2 h-1]
  real(r8),target,allocatable ::  XOPFHS(:,:,:,:,:)                  !DOP flux macropore, [g d-2 h-1]
  real(r8),target,allocatable ::  XOAFHS(:,:,:,:,:)                  !acetate flux macropore, [g d-2 h-1]

  private :: InitAllocate
  contains

  subroutine InitSoilBGCData(NumOfPlantLitrCmplxs)

  implicit none
  integer, intent(in) :: NumOfPlantLitrCmplxs

  call InitAllocate(NumOfPlantLitrCmplxs)
  end subroutine InitSoilBGCData

!------------------------------------------------------------------------------------------

  subroutine InitAllocate(NumOfPlantLitrCmplxs)
  implicit none
  integer, intent(in) :: NumOfPlantLitrCmplxs

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
  allocate(WTSTGET(NumOfPlantChemElements,JY,JX));      WTSTGET=0._r8
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
  allocate(LitrfalChemElemnts_vr(NumOfPlantChemElements,jsken,1:NumOfPlantLitrCmplxs,0:JZ,JY,JX));LitrfalChemElemnts_vr=0._r8
  allocate(trcs_VLN(ids_beg:ids_end,0:JZ,JY,JX));trcs_VLN=1._r8

  allocate(VLNHB(0:JZ,JY,JX));  VLNHB=0._r8

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
  allocate(trcg_XDFR(idg_beg:idg_end-1,JY,JX));      trcg_XDFR=0._r8
  allocate(trcg_XBLL(idg_beg:idg_end,JZ,JY,JX));  trcg_XBLL=0._r8
  allocate(XZHYS(0:JZ,JY,JX));  XZHYS=0._r8
  allocate(WaterFlowSoiMicP(3,JD,JV,JH));    WaterFlowSoiMicP=0._r8
  allocate(WaterFlowMacP(3,JD,JV,JH));   WaterFlowMacP=0._r8
  allocate(HeatFlow2Soil(3,JD,JV,JH));   HeatFlow2Soil=0._r8

  allocate(trcs_XFLS(ids_beg:ids_end,3,0:JD,JV,JH));trcs_XFLS=0._r8
  allocate(XOCFLS(1:jcplx,3,0:JD,JV,JH));XOCFLS=0._r8
  allocate(XONFLS(1:jcplx,3,0:JD,JV,JH));XONFLS=0._r8
  allocate(XOAFLS(1:jcplx,3,0:JD,JV,JH));XOAFLS=0._r8
  allocate(R3GasADTFlx(idg_beg:idg_end,3,JD,JV,JH));R3GasADTFlx=0._r8
  allocate(trcs_XFHS(ids_beg:ids_end,3,0:JD,JV,JH));trcs_XFHS=0._r8
  allocate(XOPFLS(1:jcplx,3,0:JD,JV,JH));XOPFLS=0._r8
  allocate(CPO4S(JZ,JY,JX));CPO4S(JZ,JY,JX)=0._r8
  allocate(XOCFHS(1:jcplx,3,JD,JV,JH));XOCFHS=0._r8
  allocate(XONFHS(1:jcplx,3,JD,JV,JH));XONFHS=0._r8
  allocate(XOPFHS(1:jcplx,3,JD,JV,JH));XOPFHS=0._r8
  allocate(XOAFHS(1:jcplx,3,JD,JV,JH));XOAFHS=0._r8

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
  call destroy(WTSTGET)
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
  call destroy(LitrfalChemElemnts_vr)

  call destroy(trcs_VLN)
  call destroy(trcg_XBLL)
  call destroy(VLNHB)
  call destroy(trcg_XDFR)

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
  call destroy(XZHYS)
  call destroy(WaterFlowSoiMicP)
  call destroy(WaterFlowMacP)
  call destroy(HeatFlow2Soil)
  call destroy(XOCFLS)
  call destroy(XONFLS)

  call destroy(XOAFLS)
  call destroy(XOPFLS)
  call destroy(XOCFHS)
  call destroy(XONFHS)
  call destroy(XOPFHS)
  call destroy(XOAFHS)
  end subroutine DestructSoilBGCData

end module SoilBGCDataType
