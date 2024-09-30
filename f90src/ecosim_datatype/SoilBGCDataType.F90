module SoilBGCDataType

!
! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod    , only : NumPlantChemElms
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

  real(r8),target,allocatable :: trc_soHml_vr(:,:,:,:)                 !solute mass in macropore [g d-2]
  real(r8),target,allocatable :: trc_solml_vr(:,:,:,:)                 !solute mass in micropore [g d-2]
  real(r8),target,allocatable :: trc_solcl_vr(:,:,:,:)                 !solute concentration in micropre [g m-3]
  real(r8),target,allocatable :: trc_gascl_vr(:,:,:,:)                 !gaseous concentation [g m-3]

  real(r8),target,allocatable ::  ZNFNI(:,:,:)                      !current nitrification inhibition activity
  real(r8),target,allocatable ::  ZNFN0(:,:,:)                      !initial nitrification inhibition activity
  real(r8),target,allocatable ::  ZNHUI(:,:,:)                      !current inhibition activity
  real(r8),target,allocatable ::  ZNHU0(:,:,:)                      !urea hydrolysis inhibition activity

  real(r8),target,allocatable :: trc_gasml_vr(:,:,:,:)     !layer mass of gases [g d-2]

  real(r8),target,allocatable :: PH(:,:,:)             !soil pH
  real(r8),target,allocatable :: CEC(:,:,:)            !soil cation exchange capacity	[cmol kg-1]
  real(r8),target,allocatable :: AEC(:,:,:)            !soil anion exchange capacity	[cmol kg-1]

  real(r8),target,allocatable ::  RO2UptkSoilM_vr(:,:,:,:)                     !total O2 sink, [g d-2 t-1]
  real(r8),target,allocatable ::  SurfGasFlx_col(:,:,:)                  !soil gas flux, [g d-2 h-1]
  real(r8),target,allocatable ::  AmendCFlx_CumYr_col(:,:)                         !total C amendment, [g d-2]
  real(r8),target,allocatable ::  FertNFlx_CumYr_col(:,:)                        !total fertilizer N amendment, [g d-2]
  real(r8),target,allocatable ::  FerPFlx_CumYr_col(:,:)                        !total fertilizer P amendment, [g d-2]
  real(r8),target,allocatable ::  HydroSufDOCFlx_col(:,:)                         !total surface DOC flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDOCFlx_col(:,:)                         !total subsurface DOC flux, [g d-2]
  real(r8),target,allocatable ::  LiterfalOrgM_col(:,:,:)                         !total LitrFall C, [g d-2]
  real(r8),target,allocatable ::  HydroSufDONFlx_CumYr_col(:,:)                         !total surface DON flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDONFlx_col(:,:)                         !total subsurface DON flux, [g d-2]
  real(r8),target,allocatable ::  HydroSufDOPFlx_CumYr_col(:,:)                         !total surface DOP flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDOPFlx_col(:,:)                         !total subsurface DOP flux, [g d-2]
  real(r8),target,allocatable ::  tXPO4_col(:,:)                          !total soil precipited P, [g d-2]
  real(r8),target,allocatable ::  RootResp_CumYr_col(:,:)                          !total soil autotrophic respiration, [g d-2]
  real(r8),target,allocatable ::  SedmErossLoss_CumYr_col(:,:)                        !total sediment subsurface flux, [Mg d-2]
  real(r8),target,allocatable ::  HydroSufDICFlx_col(:,:)                         !total surface DIC flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDICFlx_col(:,:)                         !total subsurface DIC flux, [g d-2]
  real(r8),target,allocatable ::  HydroSufDINFlx_CumYr_col(:,:)                         !total surface DIN flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDINFlx_col(:,:)                         !total subsurface DIN flux, [g d-2]
  real(r8),target,allocatable ::  HydroSufDIPFlx_CumYr_col(:,:)                         !total surface DIP flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDIPFlx_col(:,:)                         !total subsurface DIP flux, [g d-2]
  real(r8),target,allocatable ::  StandingDeadStrutElms_col(:,:,:)                        !total standing dead C, [g d-2]
  real(r8),target,allocatable ::  ZDRAIN(:,:)                        !total N drainage below root zone, [g d-2]
  real(r8),target,allocatable ::  PDRAIN(:,:)                        !total P drainage below root zone, [g d-2]
  real(r8),target,allocatable ::  UION(:,:)                          !total soil ion content, [mol d-2]
  real(r8),target,allocatable ::  HydroIonFlx_CumYr_col(:,:)                        !total subsurface ion flux, [mol d-2]
  real(r8),target,allocatable ::  RNutMicbTransf_vr(:,:,:,:)         !total nutrient exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_RMicbTransf_vr(:,:,:,:)       !microbial gases transformation, [g d-2 h-1]
  real(r8),target,allocatable ::  Micb_N2Fixation_vr(:,:,:)                       !net microbial N2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoDOMProd_vr(:,:,:,:,:)          !net plant+microbial DOC flux, >0 into soil [g d-2 h-1]
  real(r8),target,allocatable ::  RDOMMicProd_vr(:,:,:,:,:)          !microbial dom flux, > 0 into soil [g d-2 h-1]
  real(r8),target,allocatable ::  TMicHeterAct_vr(:,:,:)                       !total respiration of DOC+DOA in soil layer
  real(r8),target,allocatable ::  VWatMicrobAct_vr(:,:,:)                        !soil water volume occupied by microial biomass, [m3 m-3]
  real(r8),target,allocatable ::  TSens4MicbGrwoth_vr(:,:,:)                        !constraints of temperature and water potential on microbial activity, []
  real(r8),target,allocatable ::  LitrfalStrutElms_vr(:,:,:,:,:,:)                    !total LitrFall C, [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_VLN_vr(:,:,:,:)
  real(r8),target,allocatable ::  tRDOE2Die_col(:,:,:)
  real(r8),target,allocatable ::  VLNHB(:,:,:)                       !NH4 band volume fracrion, []
  real(r8),target,allocatable ::  VLNOB(:,:,:)                       !NO3 band volume fracrion, []
  real(r8),target,allocatable ::  VLPO4(:,:,:)                       !PO4 non-band volume fracrion, []
  real(r8),target,allocatable ::  VLPOB(:,:,:)                       !PO4 band volume fracrion, []
  real(r8),target,allocatable ::  BandWidthNH4_vr(:,:,:)                       !width of NH4 band, [m]
  real(r8),target,allocatable ::  BandThicknessNH4_vr(:,:,:)                       !depth of NH4 band, [m]
  real(r8),target,allocatable ::  BandWidthNO3_vr(:,:,:)                       !width of NO3 band, [m]
  real(r8),target,allocatable ::  BandThicknessNO3_vr(:,:,:)                       !depth of NO4 band, [m]
  real(r8),target,allocatable ::  BandWidthPO4_vr(:,:,:)                       !width of PO4 band, [m]
  real(r8),target,allocatable ::  BandThicknessPO4_vr(:,:,:)                       !depth of PO4 band, [m]
  real(r8),target,allocatable ::  BandDepthNH4_col(:,:)                         !total depth of NH4 band, [m]
  real(r8),target,allocatable ::  BandDepthNO3_col(:,:)                         !total depth of NO3 band, [m]
  real(r8),target,allocatable ::  BandDepthPO4_col(:,:)                         !total depth of PO4 band, [m]
  real(r8),target,allocatable ::  RNO2DmndSoilChemo_vr(:,:,:)                       !total chemodenitrification N2O uptake non-band unconstrained by N2O, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2DmndBandChemo_vr(:,:,:)                       !total chemodenitrification N2O uptake band unconstrained by N2O, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_surf_disevap_flx(:,:,:)                   !soil surface gas dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_ebu_flx_vr(:,:,:,:)                      !gas bubbling, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_ebu_flx_col(:,:,:)
  real(r8),target,allocatable ::  trcg_pltroot_flx_col(:,:,:)
  real(r8),target,allocatable ::  XZHYS(:,:,:)                       !total H+ production
  real(r8),target,allocatable ::  WaterFlowSoiMicP_3D(:,:,:,:)                       !water flux micropore, [m3 d-2 h-1]
  real(r8),target,allocatable ::  WaterFlowMacP_3D(:,:,:,:)                      !water flux macropore, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HeatFlow2Soil_3D(:,:,:,:)                      !convective heat flux micropore, [MJ d-2 h-1]

  real(r8),target,allocatable ::  trcs_Transp2MicP_3D(:,:,:,:,:)
  real(r8),target,allocatable ::  DOM_MicpTransp_3D(:,:,:,:,:,:)                  !DOC flux micropore, [g d-2 h-1]

  real(r8),target,allocatable ::  trcs_Transp2MacP_3D(:,:,:,:,:)
  real(r8),target,allocatable ::  Gas_3DAdvDif_Flx_vr(:,:,:,:,:)             !3D gaseous fluxes, [g d-2 h-1]
  real(r8),target,allocatable ::  DOM_3DMacp_Transp_flx(:,:,:,:,:,:)                  !DOC flux macropore, [g d-2 h-1]

  real(r8),target,allocatable :: RCH4ProdHydrog_vr(:,:,:)
  real(r8),target,allocatable :: RCH4ProdAcetcl_vr(:,:,:)
  real(r8),target,allocatable :: RCH4Oxi_aero_vr(:,:,:)
  real(r8),target,allocatable :: RFerment_vr(:,:,:)
  real(r8),target,allocatable :: RNH3oxi_vr(:,:,:)
  real(r8),target,allocatable :: RN2ODeniProd_vr(:,:,:)    !denitrification N2O production
  real(r8),target,allocatable :: RN2ONitProd_vr(:,:,:)
  real(r8),target,allocatable :: RN2OChemoProd_vr(:,:,:)    !chemo N2O production
  real(r8),target,allocatable :: RN2ORedux_vr(:,:,:)     !N2O reduction into N2
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

  allocate(trc_gasml_vr(idg_beg:idg_end,JZ,JY,JX)); trc_gasml_vr=0._r8
  allocate(trc_soHml_vr(ids_beg:ids_end,0:JZ,JY,JX)); trc_soHml_vr=0._r8
  allocate(trc_solml_vr(ids_beg:ids_end,0:JZ,JY,JX)); trc_solml_vr=0._r8
  allocate(trc_solcl_vr(ids_beg:ids_end,0:JZ,JY,JX)); trc_solcl_vr=0._r8
  allocate(trc_gascl_vr(idg_beg:idg_end,0:JZ,JY,JX)); trc_gascl_vr=0._r8
  allocate(tRDOE2Die_col(1:NumPlantChemElms,JY,JX)); tRDOE2Die_col=0._r8

  allocate(ZNFNI(0:JZ,JY,JX));  ZNFNI=0._r8
  allocate(ZNFN0(0:JZ,JY,JX));  ZNFN0=0._r8
  allocate(ZNHUI(0:JZ,JY,JX));  ZNHUI=0._r8
  allocate(ZNHU0(0:JZ,JY,JX));  ZNHU0=0._r8
  allocate(CPO4B(0:JZ,JY,JX));CPO4B(0:JZ,JY,JX)=0._r8

  allocate(RCH4ProdHydrog_vr(0:JZ,JY,JX)); RCH4ProdHydrog_vr=0._r8
  allocate(RCH4ProdAcetcl_vr(0:JZ,JY,JX)); RCH4ProdAcetcl_vr=0._r8
  allocate(RCH4Oxi_aero_vr(0:JZ,JY,JX)); RCH4Oxi_aero_vr=0._r8
  allocate(RFerment_vr(0:JZ,JY,JX)); RFerment_vr=0._r8
  allocate(RNH3oxi_vr(0:JZ,JY,JX)); RNH3oxi_vr=0._r8
  allocate(RN2ODeniProd_vr(0:JZ,JY,JX)); RN2ODeniProd_vr=0._r8
  allocate(RN2ONitProd_vr(0:JZ,JY,JX)); RN2ONitProd_vr=0._r8
  allocate(RN2OChemoProd_vr(0:JZ,JY,JX)); RN2OChemoProd_vr=0._r8
  allocate(RN2ORedux_vr(0:JZ,JY,JX));RN2ORedux_vr=0._r8
  allocate(PH(0:JZ,JY,JX));PH(0:JZ,JY,JX)=0._r8
  allocate(CEC(JZ,JY,JX));CEC(JZ,JY,JX)=0._r8
  allocate(AEC(JZ,JY,JX));AEC(JZ,JY,JX)=0._r8

  allocate(RO2UptkSoilM_vr(60,0:JZ,JY,JX));RO2UptkSoilM_vr=0._r8
  allocate(SurfGasFlx_col(idg_beg:idg_NH3,JY,JX));  SurfGasFlx_col=0._r8
  allocate(AmendCFlx_CumYr_col(JY,JX));       AmendCFlx_CumYr_col=0._r8
  allocate(FertNFlx_CumYr_col(JY,JX));      FertNFlx_CumYr_col=0._r8
  allocate(FerPFlx_CumYr_col(JY,JX));      FerPFlx_CumYr_col=0._r8
  allocate(HydroSufDOCFlx_col(JY,JX));       HydroSufDOCFlx_col=0._r8
  allocate(HydroSubsDOCFlx_col(JY,JX));       HydroSubsDOCFlx_col=0._r8
  allocate(LiterfalOrgM_col(NumPlantChemElms,JY,JX));       LiterfalOrgM_col=0._r8
  allocate(HydroSufDONFlx_CumYr_col(JY,JX));       HydroSufDONFlx_CumYr_col=0._r8
  allocate(HydroSubsDONFlx_col(JY,JX));       HydroSubsDONFlx_col=0._r8
  allocate(HydroSufDOPFlx_CumYr_col(JY,JX));       HydroSufDOPFlx_CumYr_col=0._r8
  allocate(HydroSubsDOPFlx_col(JY,JX));       HydroSubsDOPFlx_col=0._r8
  allocate(tXPO4_col(JY,JX));        tXPO4_col=0._r8

  allocate(RootResp_CumYr_col(JY,JX));        RootResp_CumYr_col=0._r8
  allocate(SedmErossLoss_CumYr_col(JY,JX));      SedmErossLoss_CumYr_col=0._r8
  allocate(HydroSufDICFlx_col(JY,JX));       HydroSufDICFlx_col=0._r8
  allocate(HydroSubsDICFlx_col(JY,JX));       HydroSubsDICFlx_col=0._r8
  allocate(HydroSufDINFlx_CumYr_col(JY,JX));       HydroSufDINFlx_CumYr_col=0._r8
  allocate(HydroSubsDINFlx_col(JY,JX));       HydroSubsDINFlx_col=0._r8
  allocate(HydroSufDIPFlx_CumYr_col(JY,JX));       HydroSufDIPFlx_CumYr_col=0._r8
  allocate(HydroSubsDIPFlx_col(JY,JX));       HydroSubsDIPFlx_col=0._r8
  allocate(StandingDeadStrutElms_col(NumPlantChemElms,JY,JX));      StandingDeadStrutElms_col=0._r8
  allocate(ZDRAIN(JY,JX));      ZDRAIN=0._r8
  allocate(PDRAIN(JY,JX));      PDRAIN=0._r8
  allocate(UION(JY,JX));        UION=0._r8
  allocate(HydroIonFlx_CumYr_col(JY,JX));      HydroIonFlx_CumYr_col=0._r8
  allocate(RNutMicbTransf_vr(ids_NH4B:ids_nuts_end,0:JZ,JY,JX)); RNutMicbTransf_vr=0._r8
  allocate(trcg_RMicbTransf_vr(idg_beg:idg_NH3-1,0:JZ,JY,JX)); trcg_RMicbTransf_vr=0._r8
  allocate(Micb_N2Fixation_vr(0:JZ,JY,JX));  Micb_N2Fixation_vr=0._r8

  allocate(REcoDOMProd_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));REcoDOMProd_vr=0._r8
  allocate(RDOMMicProd_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));RDOMMicProd_vr=0._r8
  allocate(TMicHeterAct_vr(0:JZ,JY,JX));  TMicHeterAct_vr=0._r8
  allocate(VWatMicrobAct_vr(0:JZ,JY,JX));   VWatMicrobAct_vr=0._r8
  allocate(TSens4MicbGrwoth_vr(0:JZ,JY,JX));   TSens4MicbGrwoth_vr=0._r8
  allocate(LitrfalStrutElms_vr(NumPlantChemElms,jsken,1:NumOfPlantLitrCmplxs,0:JZ,JY,JX));LitrfalStrutElms_vr=0._r8
  allocate(trcs_VLN_vr(ids_beg:ids_end,0:JZ,JY,JX));trcs_VLN_vr=1._r8

  allocate(VLNHB(0:JZ,JY,JX));  VLNHB=0._r8

  allocate(VLNOB(0:JZ,JY,JX));  VLNOB=0._r8
  allocate(VLPO4(0:JZ,JY,JX));  VLPO4=0._r8
  allocate(VLPOB(0:JZ,JY,JX));  VLPOB=0._r8
  allocate(BandWidthNH4_vr(JZ,JY,JX));    BandWidthNH4_vr=0._r8
  allocate(BandThicknessNH4_vr(JZ,JY,JX));    BandThicknessNH4_vr=0._r8
  allocate(BandWidthNO3_vr(JZ,JY,JX));    BandWidthNO3_vr=0._r8
  allocate(BandThicknessNO3_vr(JZ,JY,JX));    BandThicknessNO3_vr=0._r8
  allocate(BandWidthPO4_vr(JZ,JY,JX));    BandWidthPO4_vr=0._r8
  allocate(BandThicknessPO4_vr(JZ,JY,JX));    BandThicknessPO4_vr=0._r8
  allocate(BandDepthNH4_col(JY,JX));       BandDepthNH4_col=0._r8
  allocate(BandDepthNO3_col(JY,JX));       BandDepthNO3_col=0._r8
  allocate(BandDepthPO4_col(JY,JX));       BandDepthPO4_col=0._r8
  allocate(RNO2DmndSoilChemo_vr(0:JZ,JY,JX));  RNO2DmndSoilChemo_vr=0._r8
  allocate(RNO2DmndBandChemo_vr(0:JZ,JY,JX));  RNO2DmndBandChemo_vr=0._r8
  allocate(trcg_surf_disevap_flx(idg_beg:idg_end-1,JY,JX));      trcg_surf_disevap_flx=0._r8
  allocate(trcg_ebu_flx_vr(idg_beg:idg_end,JZ,JY,JX));  trcg_ebu_flx_vr=0._r8
  allocate(trcg_ebu_flx_col(idg_beg:idg_NH3,JY,JX)); trcg_ebu_flx_col=0._r8
  allocate(trcg_pltroot_flx_col(idg_beg:idg_NH3,JY,JX)); trcg_pltroot_flx_col=0._r8

  allocate(XZHYS(0:JZ,JY,JX));  XZHYS=0._r8
  allocate(WaterFlowSoiMicP_3D(3,JD,JV,JH));    WaterFlowSoiMicP_3D=0._r8
  allocate(WaterFlowMacP_3D(3,JD,JV,JH));   WaterFlowMacP_3D=0._r8
  allocate(HeatFlow2Soil_3D(3,JD,JV,JH));   HeatFlow2Soil_3D=0._r8

  allocate(trcs_Transp2MicP_3D(ids_beg:ids_end,3,0:JD,JV,JH));trcs_Transp2MicP_3D=0._r8
  allocate(DOM_MicpTransp_3D(idom_beg:idom_end,1:jcplx,3,0:JD,JV,JH));DOM_MicpTransp_3D=0._r8
  allocate(Gas_3DAdvDif_Flx_vr(idg_beg:idg_end,3,JD,JV,JH));Gas_3DAdvDif_Flx_vr=0._r8
  allocate(trcs_Transp2MacP_3D(ids_beg:ids_end,3,0:JD,JV,JH));trcs_Transp2MacP_3D=0._r8
  allocate(CPO4S(JZ,JY,JX));CPO4S(JZ,JY,JX)=0._r8
  allocate(DOM_3DMacp_Transp_flx(idom_beg:idom_end,1:jcplx,3,JD,JV,JH));DOM_3DMacp_Transp_flx=0._r8

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructSoilBGCData
  use abortutils, only : destroy

  implicit none
  
  call destroy(trcg_ebu_flx_col)
  call destroy(trcg_pltroot_flx_col)

  call destroy(tRDOE2Die_col)
  call destroy(CNH4)
  call destroy(CNO3)
  call destroy(CPO4)

  call destroy(RCH4ProdAcetcl_vr)
  call destroy(RCH4ProdHydrog_vr)
  call destroy(RCH4Oxi_aero_vr)
  call destroy(RFerment_vr)
  call destroy(RNH3oxi_vr)
  call destroy(RN2ONitProd_vr)  
  call destroy(RN2ODeniProd_vr)
  call destroy(RN2OChemoProd_vr)
  call destroy(RN2ORedux_vr)
  call destroy(trc_gasml_vr)
  call destroy(CPO4B)

  call destroy(trc_solml_vr)
  call destroy(trc_soHml_vr)

  call destroy(ZNFNI)
  call destroy(ZNFN0)
  call destroy(ZNHUI)
  call destroy(ZNHU0)

  call destroy(PH)
  call destroy(CEC)
  call destroy(AEC)
  call destroy(CPO4S)
  call destroy(RO2UptkSoilM_vr)
  call destroy(AmendCFlx_CumYr_col)
  call destroy(FertNFlx_CumYr_col)
  call destroy(FerPFlx_CumYr_col)
  call destroy(HydroSufDOCFlx_col)
  call destroy(HydroSubsDOCFlx_col)
  call destroy(LiterfalOrgM_col)
  call destroy(HydroSufDONFlx_CumYr_col)
  call destroy(HydroSubsDONFlx_col)
  call destroy(HydroSufDOPFlx_CumYr_col)
  call destroy(HydroSubsDOPFlx_col)
  call destroy(tXPO4_col)
  call destroy(RootResp_CumYr_col)
  call destroy(SedmErossLoss_CumYr_col)
  call destroy(HydroSufDICFlx_col)
  call destroy(HydroSubsDICFlx_col)
  call destroy(HydroSufDINFlx_CumYr_col)
  call destroy(HydroSubsDINFlx_col)
  call destroy(HydroSufDIPFlx_CumYr_col)
  call destroy(HydroSubsDIPFlx_col)
  call destroy(StandingDeadStrutElms_col)
  call destroy(ZDRAIN)
  call destroy(PDRAIN)
  call destroy(UION)
  call destroy(HydroIonFlx_CumYr_col)
  call destroy(Micb_N2Fixation_vr)
  call destroy(RNutMicbTransf_vr)
  call destroy(REcoDOMProd_vr)
  call destroy(RDOMMicProd_vr)
  call destroy(TMicHeterAct_vr)
  call destroy(VWatMicrobAct_vr)
  call destroy(TSens4MicbGrwoth_vr)
  call destroy(LitrfalStrutElms_vr)

  call destroy(trcs_VLN_vr)
  call destroy(trcg_ebu_flx_vr)
  call destroy(VLNHB)
  call destroy(trcg_surf_disevap_flx)

  call destroy(VLNOB)
  call destroy(VLPO4)
  call destroy(VLPOB)
  call destroy(BandWidthNH4_vr)
  call destroy(BandThicknessNH4_vr)
  call destroy(BandWidthNO3_vr)
  call destroy(BandThicknessNO3_vr)
  call destroy(BandWidthPO4_vr)
  call destroy(BandThicknessPO4_vr)
  call destroy(BandDepthNH4_col)
  call destroy(BandDepthNO3_col)
  call destroy(BandDepthPO4_col)
  call destroy(RNO2DmndSoilChemo_vr)
  call destroy(RNO2DmndBandChemo_vr)
  call destroy(XZHYS)
  call destroy(WaterFlowSoiMicP_3D)
  call destroy(WaterFlowMacP_3D)
  call destroy(HeatFlow2Soil_3D)
  call destroy(DOM_MicpTransp_3D)
  call destroy(DOM_3DMacp_Transp_flx)
  call destroy(trcg_RMicbTransf_vr)
  end subroutine DestructSoilBGCData

end module SoilBGCDataType
