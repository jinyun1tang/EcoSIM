module SoilBGCDataType

!
! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod    , only : NumOfPlantChemElmnts
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
  real(r8),target,allocatable ::  SurfGasFlx(:,:,:)                  !soil gas flux, [g d-2 h-1]
  real(r8),target,allocatable ::  UORGF(:,:)                         !total C amendment, [g d-2]
  real(r8),target,allocatable ::  UFERTN(:,:)                        !total fertilizer N amendment, [g d-2]
  real(r8),target,allocatable ::  UFERTP(:,:)                        !total fertilizer P amendment, [g d-2]
  real(r8),target,allocatable ::  HDOCQ(:,:)                         !total surface DOC flux, [g d-2]
  real(r8),target,allocatable ::  HDOCD(:,:)                         !total subsurface DOC flux, [g d-2]
  real(r8),target,allocatable ::  UXCSN(:,:)                         !total litterfall C, [g d-2]
  real(r8),target,allocatable ::  UXZSN(:,:)                         !total litterfall N, [g d-2]
  real(r8),target,allocatable ::  UXPSN(:,:)                         !total litterfall P, [g d-2]
  real(r8),target,allocatable ::  UDONQ(:,:)                         !total surface DON flux, [g d-2]
  real(r8),target,allocatable ::  HDOND(:,:)                         !total subsurface DON flux, [g d-2]
  real(r8),target,allocatable ::  UDOPQ(:,:)                         !total surface DOP flux, [g d-2]
  real(r8),target,allocatable ::  HDOPD(:,:)                         !total subsurface DOP flux, [g d-2]
  real(r8),target,allocatable ::  UPP4(:,:)                          !total soil precipited P, [g d-2]
  real(r8),target,allocatable ::  UCOP(:,:)                          !total soil autotrophic respiration, [g d-2]
  real(r8),target,allocatable ::  USEDOU(:,:)                        !total sediment subsurface flux, [Mg d-2]
  real(r8),target,allocatable ::  HDICQ(:,:)                         !total surface DIC flux, [g d-2]
  real(r8),target,allocatable ::  HDICD(:,:)                         !total subsurface DIC flux, [g d-2]
  real(r8),target,allocatable ::  UDINQ(:,:)                         !total surface DIN flux, [g d-2]
  real(r8),target,allocatable ::  HDIND(:,:)                         !total subsurface DIN flux, [g d-2]
  real(r8),target,allocatable ::  UDIPQ(:,:)                         !total surface DIP flux, [g d-2]
  real(r8),target,allocatable ::  HDIPD(:,:)                         !total subsurface DIP flux, [g d-2]
  real(r8),target,allocatable ::  StandingDeadChemElmnt_col(:,:,:)                        !total standing dead C, [g d-2]
  real(r8),target,allocatable ::  ZDRAIN(:,:)                        !total N drainage below root zone, [g d-2]
  real(r8),target,allocatable ::  PDRAIN(:,:)                        !total P drainage below root zone, [g d-2]
  real(r8),target,allocatable ::  UION(:,:)                          !total soil ion content, [mol d-2]
  real(r8),target,allocatable ::  UIONOU(:,:)                        !total subsurface ion flux, [mol d-2]
  real(r8),target,allocatable ::  RNutMicbTransf_vr(:,:,:,:)         !total nutrient exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_RMicbTransf_vr(:,:,:,:)       !microbial gases transformation, [g d-2 h-1]
  real(r8),target,allocatable ::  Micb_N2Fixation_vr(:,:,:)                       !net microbial N2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  RDOM_micb_flx(:,:,:,:,:)                     !net microbial DOC flux, [g d-2 h-1]
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
  real(r8),target,allocatable ::  trcg_surf_disevap_flx(:,:,:)                   !soil surface gas dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_ebu_flx_vr(:,:,:,:)                      !CO2 bubbling, [g d-2 h-1]
  real(r8),target,allocatable ::  XZHYS(:,:,:)                       !total H+ production
  real(r8),target,allocatable ::  WaterFlowSoiMicP(:,:,:,:)                       !water flux micropore, [m3 d-2 h-1]
  real(r8),target,allocatable ::  WaterFlowMacP(:,:,:,:)                      !water flux macropore, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HeatFlow2Soil(:,:,:,:)                      !convective heat flux micropore, [MJ d-2 h-1]

  real(r8),target,allocatable ::  trcs_3DTransp2MicP(:,:,:,:,:)
  real(r8),target,allocatable ::  DOM_3DMicp_Transp_flx(:,:,:,:,:,:)                  !DOC flux micropore, [g d-2 h-1]

  real(r8),target,allocatable ::  trcs_3DTransp2MacP(:,:,:,:,:)
  real(r8),target,allocatable ::  Gas_3DAdvDif_Flx_vr(:,:,:,:,:)             !3D gaseous fluxes, [g d-2 h-1]
  real(r8),target,allocatable ::  DOM_3DMacp_Transp_flx(:,:,:,:,:,:)                  !DOC flux macropore, [g d-2 h-1]

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
  allocate(SurfGasFlx(idg_beg:idg_NH3,JY,JX));  SurfGasFlx=0._r8
  allocate(UORGF(JY,JX));       UORGF=0._r8
  allocate(UFERTN(JY,JX));      UFERTN=0._r8
  allocate(UFERTP(JY,JX));      UFERTP=0._r8
  allocate(HDOCQ(JY,JX));       HDOCQ=0._r8
  allocate(HDOCD(JY,JX));       HDOCD=0._r8
  allocate(UXCSN(JY,JX));       UXCSN=0._r8
  allocate(UXZSN(JY,JX));       UXZSN=0._r8
  allocate(UXPSN(JY,JX));       UXPSN=0._r8
  allocate(UDONQ(JY,JX));       UDONQ=0._r8
  allocate(HDOND(JY,JX));       HDOND=0._r8
  allocate(UDOPQ(JY,JX));       UDOPQ=0._r8
  allocate(HDOPD(JY,JX));       HDOPD=0._r8
  allocate(UPP4(JY,JX));        UPP4=0._r8

  allocate(UCOP(JY,JX));        UCOP=0._r8
  allocate(USEDOU(JY,JX));      USEDOU=0._r8
  allocate(HDICQ(JY,JX));       HDICQ=0._r8
  allocate(HDICD(JY,JX));       HDICD=0._r8
  allocate(UDINQ(JY,JX));       UDINQ=0._r8
  allocate(HDIND(JY,JX));       HDIND=0._r8
  allocate(UDIPQ(JY,JX));       UDIPQ=0._r8
  allocate(HDIPD(JY,JX));       HDIPD=0._r8
  allocate(StandingDeadChemElmnt_col(NumOfPlantChemElmnts,JY,JX));      StandingDeadChemElmnt_col=0._r8
  allocate(ZDRAIN(JY,JX));      ZDRAIN=0._r8
  allocate(PDRAIN(JY,JX));      PDRAIN=0._r8
  allocate(UION(JY,JX));        UION=0._r8
  allocate(UIONOU(JY,JX));      UIONOU=0._r8
  allocate(RNutMicbTransf_vr(ids_NH4B:ids_nuts_end,0:JZ,JY,JX)); RNutMicbTransf_vr=0._r8
  allocate(trcg_RMicbTransf_vr(idg_beg:idg_NH3-1,0:JZ,JY,JX)); trcg_RMicbTransf_vr=0._r8
  allocate(Micb_N2Fixation_vr(0:JZ,JY,JX));  Micb_N2Fixation_vr=0._r8

  allocate(RDOM_micb_flx(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));RDOM_micb_flx=0._r8
  allocate(TOQCK(0:JZ,JY,JX));  TOQCK=0._r8
  allocate(VOLQ(0:JZ,JY,JX));   VOLQ=0._r8
  allocate(TFNQ(0:JZ,JY,JX));   TFNQ=0._r8
  allocate(LitrfalChemElemnts_vr(NumOfPlantChemElmnts,jsken,1:NumOfPlantLitrCmplxs,0:JZ,JY,JX));LitrfalChemElemnts_vr=0._r8
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
  allocate(trcg_surf_disevap_flx(idg_beg:idg_end-1,JY,JX));      trcg_surf_disevap_flx=0._r8
  allocate(trcg_ebu_flx_vr(idg_beg:idg_end,JZ,JY,JX));  trcg_ebu_flx_vr=0._r8
  allocate(XZHYS(0:JZ,JY,JX));  XZHYS=0._r8
  allocate(WaterFlowSoiMicP(3,JD,JV,JH));    WaterFlowSoiMicP=0._r8
  allocate(WaterFlowMacP(3,JD,JV,JH));   WaterFlowMacP=0._r8
  allocate(HeatFlow2Soil(3,JD,JV,JH));   HeatFlow2Soil=0._r8

  allocate(trcs_3DTransp2MicP(ids_beg:ids_end,3,0:JD,JV,JH));trcs_3DTransp2MicP=0._r8
  allocate(DOM_3DMicp_Transp_flx(idom_beg:idom_end,1:jcplx,3,0:JD,JV,JH));DOM_3DMicp_Transp_flx=0._r8
  allocate(Gas_3DAdvDif_Flx_vr(idg_beg:idg_end,3,JD,JV,JH));Gas_3DAdvDif_Flx_vr=0._r8
  allocate(trcs_3DTransp2MacP(ids_beg:ids_end,3,0:JD,JV,JH));trcs_3DTransp2MacP=0._r8
  allocate(CPO4S(JZ,JY,JX));CPO4S(JZ,JY,JX)=0._r8
  allocate(DOM_3DMacp_Transp_flx(idom_beg:idom_end,1:jcplx,3,JD,JV,JH));DOM_3DMacp_Transp_flx=0._r8

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
  call destroy(UORGF)
  call destroy(UFERTN)
  call destroy(UFERTP)
  call destroy(HDOCQ)
  call destroy(HDOCD)
  call destroy(UXCSN)
  call destroy(UXZSN)
  call destroy(UXPSN)
  call destroy(UDONQ)
  call destroy(HDOND)
  call destroy(UDOPQ)
  call destroy(HDOPD)
  call destroy(UPP4)
  call destroy(UCOP)
  call destroy(USEDOU)
  call destroy(HDICQ)
  call destroy(HDICD)
  call destroy(UDINQ)
  call destroy(HDIND)
  call destroy(UDIPQ)
  call destroy(HDIPD)
  call destroy(StandingDeadChemElmnt_col)
  call destroy(ZDRAIN)
  call destroy(PDRAIN)
  call destroy(UION)
  call destroy(UIONOU)
  call destroy(Micb_N2Fixation_vr)
  call destroy(RNutMicbTransf_vr)
  call destroy(RDOM_micb_flx)
  call destroy(TOQCK)
  call destroy(VOLQ)
  call destroy(TFNQ)
  call destroy(LitrfalChemElemnts_vr)

  call destroy(trcs_VLN)
  call destroy(trcg_ebu_flx_vr)
  call destroy(VLNHB)
  call destroy(trcg_surf_disevap_flx)

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
  call destroy(DOM_3DMicp_Transp_flx)
  call destroy(DOM_3DMacp_Transp_flx)
  call destroy(trcg_RMicbTransf_vr)
  end subroutine DestructSoilBGCData

end module SoilBGCDataType
