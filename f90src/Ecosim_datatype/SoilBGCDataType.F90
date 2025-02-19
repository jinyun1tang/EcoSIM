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

  real(r8),target,allocatable ::  CNH4_vr(:,:,:)                        !soil NH4 content, [mg kg-1]
  real(r8),target,allocatable ::  CNO3_vr(:,:,:)                        !soil NO3 content, [mg kg-1]
  real(r8),target,allocatable ::  CPO4_vr(:,:,:)                        !soil PO4 content, [mg kg-1]

  real(r8),target,allocatable :: CPO4B_vr(:,:,:)                        !PO4 concentration band micropore	[g m-3]
  real(r8),target,allocatable :: CPO4S_vr(:,:,:)                        !PO4 concentration non-band micropore	[g m-3]

  real(r8),target,allocatable :: trcs_soHml_vr(:,:,:,:)               !solute mass in macropore [g d-2]
  real(r8),target,allocatable :: trcs_solml_vr(:,:,:,:)               !solute mass in micropore [g d-2]
  real(r8),target,allocatable :: trc_solcl_vr(:,:,:,:)               !solute concentration in micropre [g m-3]
  real(r8),target,allocatable :: trcg_gascl_vr(:,:,:,:)               !gaseous concentation in micropore [g m-3]
  real(r8),target,allocatable :: tRHydlySOM_vr(:,:,:,:)              !solid SOM hydrolysis rate [g/m2/hr]
  real(r8),target,allocatable :: tRHydlyBioReSOM_vr(:,:,:,:)         !microbial residual hydrolysis rate [g/m2/hr]
  real(r8),target,allocatable :: tRHydlySoprtOM_vr(:,:,:,:)          !sorbed OM hydrolysis rate [g/m2/hr]

  real(r8),target,allocatable ::  ZNFNI_vr(:,:,:)                       !current nitrification inhibition activity
  real(r8),target,allocatable ::  ZNFN0_vr(:,:,:)                       !initial nitrification inhibition activity
  real(r8),target,allocatable ::  ZNHUI_vr(:,:,:)                       !current inhibition activity
  real(r8),target,allocatable ::  ZNHU0_vr(:,:,:)                       !urea hydrolysis inhibition activity
  real(r8),target,allocatable ::  trcg_soilMass_col(:,:,:)           !column integrated volatile tracer mass in soil at the moment [g d-2]
  real(r8),target,allocatable ::  trcg_soilMass_beg_col(:,:,:)           !column integrated volatile tracer mass in soil at the moment [g d-2]  
  real(r8),target,allocatable ::  trcg_gasml_vr(:,:,:,:)             !layer mass of gases in micropores [g d-2]
  real(r8),target,allocatable ::  trcg_TotalMass_beg_col(:,:,:)      !column integrated volatile tracer mass at the begining of time step [g d-2]
  real(r8),target,allocatable ::  trcg_TotalMass_col(:,:,:)          !column integrated volatile tracer mass at the moment [g d-2]
  real(r8),target,allocatable ::  PH_vr(:,:,:)                       !soil pH
  real(r8),target,allocatable ::  CEC_vr(:,:,:)                      !soil cation exchange capacity	[cmol kg-1]
  real(r8),target,allocatable ::  AEC_vr(:,:,:)                      !soil anion exchange capacity	[cmol kg-1]
  real(r8),target,allocatable ::  TempSensDecomp_vr(:,:,:)           !temperature dependense of microbial activity
  real(r8),target,allocatable ::  MoistSensDecomp_vr(:,:,:)          !moisture dependence of microbial activity
  real(r8),target,allocatable ::  SurfGasDifFlx_col(:,:,:)           !surface gas flux in advection+diffusion [g d-2 h-1]
  real(r8),target,allocatable ::  RO2UptkSoilM_vr(:,:,:,:)           !total O2 sink, [g d-2 h-1]
  real(r8),target,allocatable ::  SurfGasEmisFlx_col(:,:,:)          !surface gas flux, including diffusion, ebullition, wet deposition and plant transp [g d-2 h-1]
  real(r8),target,allocatable ::  GasHydroLossFlx_col(:,:,:)         !hydrological loss of volatile tracers [g d-2 h-1]
  real(r8),target,allocatable ::  AmendCFlx_CumYr_col(:,:)           !total C amendment, [g d-2]
  real(r8),target,allocatable ::  FertNFlx_CumYr_col(:,:)            !total fertilizer N amendment, [g d-2]
  real(r8),target,allocatable ::  FerPFlx_CumYr_col(:,:)             !total fertilizer P amendment, [g d-2]
  real(r8),target,allocatable ::  HydroSufDOCFlx_col(:,:)            !total surface DOC flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDOCFlx_col(:,:)           !total subsurface DOC flux, [g d-2]
  real(r8),target,allocatable ::  LiterfalOrgM_col(:,:,:)            !total LitrFall C, [g d-2]
  real(r8),target,allocatable ::  HydroSufDONFlx_CumYr_col(:,:)      !total surface DON flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDONFlx_col(:,:)           !total subsurface DON flux, [g d-2]
  real(r8),target,allocatable ::  HydroSufDOPFlx_CumYr_col(:,:)      !total surface DOP flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDOPFlx_col(:,:)           !total subsurface DOP flux, [g d-2]
  real(r8),target,allocatable ::  tXPO4_col(:,:)                     !total soil precipited P, [g d-2]
  real(r8),target,allocatable ::  RootResp_CumYr_col(:,:)            !total soil autotrophic respiration, [g d-2]
  real(r8),target,allocatable ::  SedmErossLoss_CumYr_col(:,:)       !total sediment subsurface flux, [Mg d-2]
  real(r8),target,allocatable ::  HydroSufDICFlx_col(:,:)            !total surface DIC flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDICFlx_col(:,:)           !total subsurface DIC flux, [g d-2]
  real(r8),target,allocatable ::  HydroSufDINFlx_CumYr_col(:,:)      !total surface DIN flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDINFlx_col(:,:)           !total subsurface DIN flux, [g d-2]
  real(r8),target,allocatable ::  HydroSufDIPFlx_CumYr_col(:,:)      !total surface DIP flux, [g d-2]
  real(r8),target,allocatable ::  HydroSubsDIPFlx_col(:,:)           !total subsurface DIP flux, [g d-2]
  real(r8),target,allocatable ::  StandingDeadStrutElms_col(:,:,:)   !total standing dead C, [g d-2]
  real(r8),target,allocatable ::  ZDRAIN_col(:,:)                    !total N drainage below root zone, [g d-2]
  real(r8),target,allocatable ::  PDRAIN_col(:,:)                    !total P drainage below root zone, [g d-2]
  real(r8),target,allocatable ::  UION_col(:,:)                      !total soil ion content, [mol d-2]
  real(r8),target,allocatable ::  HydroIonFlx_CumYr_col(:,:)         !total subsurface ion flux, [mol d-2]
  real(r8),target,allocatable ::  RNut_MicbRelease_vr(:,:,:,:)         !total nutrient exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_RMicbUptake_vr(:,:,:,:)       !microbial gases transformation, [g d-2 h-1]
  real(r8),target,allocatable ::  Micb_N2Fixation_vr(:,:,:)          !net microbial N2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoDOMProd_vr(:,:,:,:,:)          !net plant+microbial DOC flux, >0 into soil [g d-2 h-1]
  real(r8),target,allocatable ::  RDOMMicProd_vr(:,:,:,:,:)          !microbial dom flux, > 0 into soil [g d-2 h-1]
  real(r8),target,allocatable ::  TMicHeterActivity_vr(:,:,:)        !total respiration of DOC+DOA in soil layer
  real(r8),target,allocatable ::  VWatMicrobAct_vr(:,:,:)            !soil water volume occupied by microial biomass, [m3 m-3]
  real(r8),target,allocatable ::  TSens4MicbGrwoth_vr(:,:,:)         !constraints of temperature and water potential on microbial activity, []
  real(r8),target,allocatable ::  LitrfalStrutElms_vr(:,:,:,:,:,:)   !total LitrFall C, [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_VLN_vr(:,:,:,:)               !effective volume fraction of nutrient solutes [0-1]
  real(r8),target,allocatable ::  tRDIM2DOM_col(:,:,:)               !conversion flux from DIM into DOM [g d-2 h-1]
  real(r8),target,allocatable ::  Gas_NetProd_col(:,:,:)             !net production of gas [g d-2 h-1]
  real(r8),target,allocatable ::  OxyDecompLimiter_vr(:,:,:)         !decomposer oxygen limitation
  real(r8),target,allocatable ::  RO2DecompUptk_vr(:,:,:)            !decompoer oxygen uptake rate
  real(r8),target,allocatable ::  BandWidthNH4_vr(:,:,:)             !width of NH4 band, [m]
  real(r8),target,allocatable ::  BandThicknessNH4_vr(:,:,:)         !depth of NH4 band, [m]
  real(r8),target,allocatable ::  BandWidthNO3_vr(:,:,:)             !width of NO3 band, [m]
  real(r8),target,allocatable ::  BandThicknessNO3_vr(:,:,:)         !depth of NO4 band, [m]
  real(r8),target,allocatable ::  BandWidthPO4_vr(:,:,:)             !width of PO4 band, [m]
  real(r8),target,allocatable ::  BandThicknessPO4_vr(:,:,:)         !depth of PO4 band, [m]
  real(r8),target,allocatable ::  BandDepthNH4_col(:,:)              !total depth of NH4 band, [m]
  real(r8),target,allocatable ::  BandDepthNO3_col(:,:)              !total depth of NO3 band, [m]
  real(r8),target,allocatable ::  BandDepthPO4_col(:,:)              !total depth of PO4 band, [m]
  real(r8),target,allocatable ::  RNO2DmndSoilChemo_vr(:,:,:)        !total chemodenitrification N2O uptake non-band unconstrained by N2O, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2DmndBandChemo_vr(:,:,:)        !total chemodenitrification N2O uptake band unconstrained by N2O, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_DisolEvap_Atm2Litr_flx(:,:,:) !soil surface gas dissolution (+ve) - volatilization (-ve), [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_ebu_flx_vr(:,:,:,:)           !<0., active gas bubbling, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_ebu_flx_col(:,:,:)            !vertically integrated ebullition flux [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_pltroot_flx_col(:,:,:)        !plant-aided gas transport flux [g d-2 h-1]
  real(r8),target,allocatable ::  RProd_Hp_vr(:,:,:)                 !total H+ production 
  real(r8),target,allocatable ::  WaterFlowSoiMicP_3D(:,:,:,:)       !water flux micropore, [m3 d-2 h-1]
  real(r8),target,allocatable ::  WaterFlowSoiMacP_3D(:,:,:,:)       !water flux macropore, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HeatFlow2Soil_3D(:,:,:,:)          !convective heat flux micropore, [MJ d-2 h-1]
  real(r8),target,allocatable ::  Gas_Prod_TP_cumRes_col(:,:,:)      !Cumulative difference in gas belowground production and surface flux [g d-2]
  real(r8),target,allocatable ::  trcs_TransptMicP_3D(:,:,:,:,:)     !tracer solute transport in micropore [g d-2 h-1]
  real(r8),target,allocatable ::  DOM_MicpTransp_3D(:,:,:,:,:,:)     !DOC flux micropore, [g d-2 h-1]
  real(r8),target,allocatable ::  Gas_WetDeposition_col(:,:,:)       !wet gas deposition due to irrigation and rainfall [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_TransptMacP_3D(:,:,:,:,:)     !tracer solute transport in macropore [g d-2 h-1]
  real(r8),target,allocatable ::  Gas_AdvDif_Flx_3D(:,:,:,:,:)       !3D gaseous fluxes, [g d-2 h-1]
  real(r8),target,allocatable ::  DOM_Macp_Transp_flx_3D(:,:,:,:,:,:) !DOC flux macropore, [g d-2 h-1]
  real(r8),target,allocatable ::  Soil_Gas_pressure_vr(:,:,:)         !soil gas pressure, [Pa]
  real(r8),target,allocatable ::  CO2_Gas_Frac_vr(:,:,:)              !volumetric concentation of gaseous CO2 [ppmv]
  real(r8),target,allocatable ::  O2_Gas_Frac_vr(:,:,:)              !volumetric concentation of gaseous O2 [ppmv]
  real(r8),target,allocatable ::  Ar_Gas_frac_vr(:,:,:)               !volumetric concentation of Ar gas  [ppmv]
  real(r8),target,allocatable ::  CH4_gas_frac_vr(:,:,:)              !volumetric concentation of CH4 gas [ppmv]
  real(r8),target,allocatable :: RCH4ProdHydrog_vr(:,:,:)             !Hydrogenotrophic CH4 production rate [gC d-2 h-1]
  real(r8),target,allocatable :: RCH4ProdAcetcl_vr(:,:,:)             !Acetoclastic CH4 production rate [gC d-2 h-1]
  real(r8),target,allocatable :: RCH4Oxi_aero_vr(:,:,:)               !Aerobic CH4 oxidation rate [gC d-2 h-1]
  real(r8),target,allocatable :: RFerment_vr(:,:,:)                   !Fermentation rate [gC d-2 h-1]
  real(r8),target,allocatable :: RNH3oxi_vr(:,:,:)                    !NH3 oxidation rate [gN d-2 h-1]
  real(r8),target,allocatable :: RN2ODeniProd_vr(:,:,:)              !denitrification N2O production [gN d-2 h-1]
  real(r8),target,allocatable :: RN2ONitProd_vr(:,:,:)               !Nitrification N2O produciton rate [gN d-2 h-1]
  real(r8),target,allocatable :: RN2OChemoProd_vr(:,:,:)             !chemo N2O production [gN d-2 h-1]
  real(r8),target,allocatable :: RN2ORedux_vr(:,:,:)                 !N2O reduction into N2  [gN d-2 h-1]
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

  allocate(CNH4_vr(JZ,JY,JX));     CNH4_vr=0._r8
  allocate(CNO3_vr(JZ,JY,JX));     CNO3_vr=0._r8
  allocate(CPO4_vr(JZ,JY,JX));     CPO4_vr=0._r8

  allocate(trcg_soilMass_col(idg_beg:idg_end,JY,JX)); trcg_soilMass_col=0._r8
  allocate(trcg_soilMass_beg_col(idg_beg:idg_end,JY,JX)); trcg_soilMass_beg_col=0._r8
  allocate(Gas_Prod_TP_cumRes_col(idg_beg:idg_NH3,JY,JX)); Gas_Prod_TP_cumRes_col=0._r8
  allocate(trcg_TotalMass_beg_col(idg_beg:idg_end,JY,JX)); trcg_TotalMass_beg_col=0._r8
  allocate(trcg_TotalMass_col(idg_beg:idg_end,JY,JX)); trcg_TotalMass_col=0._r8
  allocate(trcg_gasml_vr(idg_beg:idg_NH3,JZ,JY,JX)); trcg_gasml_vr=0._r8
  allocate(trcs_soHml_vr(ids_beg:ids_end,JZ,JY,JX)); trcs_soHml_vr=0._r8
  allocate(trcs_solml_vr(ids_beg:ids_end,0:JZ,JY,JX)); trcs_solml_vr=0._r8
  allocate(trc_solcl_vr(ids_beg:ids_end,0:JZ,JY,JX)); trc_solcl_vr=0._r8
  allocate(trcg_gascl_vr(idg_beg:idg_NH3,0:JZ,JY,JX)); trcg_gascl_vr=0._r8
  allocate(tRDIM2DOM_col(1:NumPlantChemElms,JY,JX)); tRDIM2DOM_col=0._r8
  allocate(Gas_WetDeposition_col(idg_beg:idg_NH3,JY,JX)); Gas_WetDeposition_col=0._r8
  allocate(tRHydlySOM_vr(1:NumPlantChemElms,0:JZ,JY,JX)); tRHydlySOM_vr=0._r8
  allocate(tRHydlyBioReSOM_vr(1:NumPlantChemElms,0:JZ,JY,JX));tRHydlyBioReSOM_vr=0._r8
  allocate(tRHydlySoprtOM_vr(1:NumPlantChemElms,0:JZ,JY,JX));      tRHydlySoprtOM_vr=0._r8

  allocate(ZNFNI_vr(0:JZ,JY,JX));  ZNFNI_vr=0._r8
  allocate(ZNFN0_vr(0:JZ,JY,JX));  ZNFN0_vr=0._r8
  allocate(ZNHUI_vr(0:JZ,JY,JX));  ZNHUI_vr  =0._r8
  allocate(ZNHU0_vr(0:JZ,JY,JX));  ZNHU0_vr=0._r8
  allocate(CPO4B_vr(0:JZ,JY,JX));CPO4B_vr(0:JZ,JY,JX)=0._r8
  allocate(O2_Gas_Frac_vr(1:JZ,JY,JX)) ; O2_Gas_Frac_vr = 0._r8
  allocate(CO2_Gas_Frac_vr(1:JZ,JY,JX)) ; CO2_Gas_Frac_vr = 0._r8
  allocate(CH4_Gas_Frac_vr(1:JZ,JY,JX)) ; CH4_Gas_Frac_vr = 0._r8  
  allocate(Ar_Gas_Frac_vr(1:JZ,JY,JX)) ; Ar_Gas_Frac_vr = 0._r8  
  allocate(Soil_Gas_pressure_vr(1:JZ,JY,JX)); Soil_Gas_pressure_vr=0._r8  
  allocate(RCH4ProdHydrog_vr(0:JZ,JY,JX)); RCH4ProdHydrog_vr=0._r8
  allocate(RCH4ProdAcetcl_vr(0:JZ,JY,JX)); RCH4ProdAcetcl_vr=0._r8
  allocate(RCH4Oxi_aero_vr(0:JZ,JY,JX)); RCH4Oxi_aero_vr=0._r8
  allocate(RFerment_vr(0:JZ,JY,JX)); RFerment_vr=0._r8
  allocate(RNH3oxi_vr(0:JZ,JY,JX)); RNH3oxi_vr=0._r8
  allocate(RN2ODeniProd_vr(0:JZ,JY,JX)); RN2ODeniProd_vr=0._r8
  allocate(RN2ONitProd_vr(0:JZ,JY,JX)); RN2ONitProd_vr=0._r8
  allocate(RN2OChemoProd_vr(0:JZ,JY,JX)); RN2OChemoProd_vr=0._r8
  allocate(RN2ORedux_vr(0:JZ,JY,JX));RN2ORedux_vr=0._r8
  allocate(PH_vr(0:JZ,JY,JX));PH_vr(0:JZ,JY,JX)=0._r8
  allocate(CEC_vr(JZ,JY,JX));CEC_vr(JZ,JY,JX)=0._r8
  allocate(AEC_vr(JZ,JY,JX));AEC_vr(JZ,JY,JX)=0._r8

  allocate(RO2UptkSoilM_vr(60,0:JZ,JY,JX));RO2UptkSoilM_vr=0._r8
  allocate(GasHydroLossFlx_col(idg_beg:idg_end,JY,JX)); GasHydroLossFlx_col=0._r8
  allocate(SurfGasEmisFlx_col(idg_beg:idg_NH3,JY,JX));  SurfGasEmisFlx_col=0._r8
  allocate(SurfGasDifFlx_col(idg_beg:idg_NH3,JY,JX)); SurfGasDifFlx_col=0._r8
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
  allocate(Gas_NetProd_col(idg_beg:idg_NH3,JY,JX)); Gas_NetProd_col=0._r8
  allocate(RootResp_CumYr_col(JY,JX));        RootResp_CumYr_col=0._r8
  allocate(SedmErossLoss_CumYr_col(JY,JX));      SedmErossLoss_CumYr_col=0._r8
  allocate(HydroSufDICFlx_col(JY,JX));       HydroSufDICFlx_col=0._r8
  allocate(HydroSubsDICFlx_col(JY,JX));       HydroSubsDICFlx_col=0._r8
  allocate(HydroSufDINFlx_CumYr_col(JY,JX));       HydroSufDINFlx_CumYr_col=0._r8
  allocate(HydroSubsDINFlx_col(JY,JX));       HydroSubsDINFlx_col=0._r8
  allocate(HydroSufDIPFlx_CumYr_col(JY,JX));       HydroSufDIPFlx_CumYr_col=0._r8
  allocate(HydroSubsDIPFlx_col(JY,JX));       HydroSubsDIPFlx_col=0._r8
  allocate(StandingDeadStrutElms_col(NumPlantChemElms,JY,JX));      StandingDeadStrutElms_col=0._r8
  allocate(ZDRAIN_col(JY,JX));      ZDRAIN_col=0._r8
  allocate(PDRAIN_col(JY,JX));      PDRAIN_col=0._r8
  allocate(UION_col(JY,JX));        UION_col=0._r8
  allocate(HydroIonFlx_CumYr_col(JY,JX));      HydroIonFlx_CumYr_col=0._r8
  allocate(RNut_MicbRelease_vr(ids_NH4B:ids_nuts_end,0:JZ,JY,JX)); RNut_MicbRelease_vr=0._r8
  allocate(trcs_RMicbUptake_vr(idg_beg:idg_NH3-1,0:JZ,JY,JX)); trcs_RMicbUptake_vr=0._r8
  allocate(Micb_N2Fixation_vr(0:JZ,JY,JX));  Micb_N2Fixation_vr=0._r8

  allocate(REcoDOMProd_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));REcoDOMProd_vr=0._r8
  allocate(RDOMMicProd_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));RDOMMicProd_vr=0._r8
  allocate(TMicHeterActivity_vr(0:JZ,JY,JX));  TMicHeterActivity_vr=0._r8
  allocate(VWatMicrobAct_vr(0:JZ,JY,JX));   VWatMicrobAct_vr=0._r8
  allocate(TSens4MicbGrwoth_vr(0:JZ,JY,JX));   TSens4MicbGrwoth_vr=0._r8
  allocate(LitrfalStrutElms_vr(NumPlantChemElms,jsken,1:NumOfPlantLitrCmplxs,0:JZ,JY,JX));LitrfalStrutElms_vr=0._r8
  allocate(trcs_VLN_vr(ids_beg:ids_end,0:JZ,JY,JX));trcs_VLN_vr=1._r8

  allocate(BandWidthNH4_vr(JZ,JY,JX));    BandWidthNH4_vr=0._r8
  allocate(OxyDecompLimiter_vr(0:JZ,JY,JX)); OxyDecompLimiter_vr=0._r8
  allocate(RO2DecompUptk_vr(0:JZ,JY,JX)); RO2DecompUptk_vr=0._r8
  allocate(BandThicknessNH4_vr(JZ,JY,JX));    BandThicknessNH4_vr=0._r8
  allocate(BandWidthNO3_vr(JZ,JY,JX));    BandWidthNO3_vr=0._r8
  allocate(TempSensDecomp_vr(0:JZ,JY,JX)); TempSensDecomp_vr=0._r8
  allocate(MoistSensDecomp_vr(0:JZ,JY,JX)); MoistSensDecomp_vr=0._r8
  allocate(BandThicknessNO3_vr(JZ,JY,JX));    BandThicknessNO3_vr=0._r8
  allocate(BandWidthPO4_vr(JZ,JY,JX));    BandWidthPO4_vr=0._r8
  allocate(BandThicknessPO4_vr(JZ,JY,JX));    BandThicknessPO4_vr=0._r8
  allocate(BandDepthNH4_col(JY,JX));       BandDepthNH4_col=0._r8
  allocate(BandDepthNO3_col(JY,JX));       BandDepthNO3_col=0._r8
  allocate(BandDepthPO4_col(JY,JX));       BandDepthPO4_col=0._r8
  allocate(RNO2DmndSoilChemo_vr(0:JZ,JY,JX));  RNO2DmndSoilChemo_vr=0._r8
  allocate(RNO2DmndBandChemo_vr(0:JZ,JY,JX));  RNO2DmndBandChemo_vr=0._r8
  allocate(trcg_DisolEvap_Atm2Litr_flx(idg_beg:idg_NH3,JY,JX));      trcg_DisolEvap_Atm2Litr_flx=0._r8
  allocate(trcg_ebu_flx_vr(idg_beg:idg_end,JZ,JY,JX));  trcg_ebu_flx_vr=0._r8
  allocate(trcg_ebu_flx_col(idg_beg:idg_NH3,JY,JX)); trcg_ebu_flx_col=0._r8
  allocate(trcg_pltroot_flx_col(idg_beg:idg_NH3,JY,JX)); trcg_pltroot_flx_col=0._r8

  allocate(RProd_Hp_vr(0:JZ,JY,JX));  RProd_Hp_vr=0._r8
  allocate(WaterFlowSoiMicP_3D(3,JD,JV,JH));    WaterFlowSoiMicP_3D=0._r8
  allocate(WaterFlowSoiMacP_3D(3,JD,JV,JH));   WaterFlowSoiMacP_3D=0._r8
  allocate(HeatFlow2Soil_3D(3,JD,JV,JH));   HeatFlow2Soil_3D=0._r8

  allocate(trcs_TransptMicP_3D(ids_beg:ids_end,3,0:JD,JV,JH));trcs_TransptMicP_3D=0._r8
  allocate(DOM_MicpTransp_3D(idom_beg:idom_end,1:jcplx,3,0:JD,JV,JH));DOM_MicpTransp_3D=0._r8
  allocate(Gas_AdvDif_Flx_3D(idg_beg:idg_end,3,JD,JV,JH));Gas_AdvDif_Flx_3D=0._r8
  allocate(trcs_TransptMacP_3D(ids_beg:ids_end,3,0:JD,JV,JH));trcs_TransptMacP_3D=0._r8
  allocate(CPO4S_vr(JZ,JY,JX));CPO4S_vr(JZ,JY,JX)=0._r8
  allocate(DOM_Macp_Transp_flx_3D(idom_beg:idom_end,1:jcplx,3,JD,JV,JH));DOM_Macp_Transp_flx_3D=0._r8

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructSoilBGCData
  use abortutils, only : destroy

  implicit none
  call destroy(O2_Gas_Frac_vr)
  call destroy(CO2_Gas_Frac_vr)
  call destroy(CH4_Gas_Frac_vr)  
  call destroy(Ar_Gas_Frac_vr)
  call destroy(Gas_NetProd_col)
  call destroy(Gas_WetDeposition_col)
  call destroy(Gas_Prod_TP_cumRes_col)
  call destroy(trcg_ebu_flx_col)
  call destroy(trcg_pltroot_flx_col)
  call destroy(tRDIM2DOM_col)
  call destroy(CNH4_vr)
  call destroy(CNO3_vr)
  call destroy(CPO4_vr)
  call destroy(trc_solcl_vr)
  call destroy(tRHydlySOM_vr)
  call destroy(tRHydlyBioReSOM_vr)
  call destroy(tRHydlySoprtOM_vr)
  call destroy(trcg_gascl_vr)
  call destroy(RCH4ProdAcetcl_vr)
  call destroy(RCH4ProdHydrog_vr)
  call destroy(RCH4Oxi_aero_vr)
  call destroy(RFerment_vr)
  call destroy(RNH3oxi_vr)
  call destroy(RN2ONitProd_vr)  
  call destroy(RN2ODeniProd_vr)
  call destroy(RN2OChemoProd_vr)
  call destroy(RN2ORedux_vr)
  call destroy(trcg_gasml_vr)
  call destroy(CPO4B_vr)
  call destroy(OxyDecompLimiter_vr)
  call destroy(RO2DecompUptk_vr)
  call destroy(trcs_solml_vr)
  call destroy(trcs_soHml_vr)

  call destroy(ZNFNI_vr)
  call destroy(ZNFN0_vr)
  call destroy(ZNHUI_vr)
  call destroy(ZNHU0_vr)
  call destroy(Soil_Gas_pressure_vr)
  call destroy(PH_vr)
  call destroy(CEC_vr)
  call destroy(AEC_vr)
  call destroy(CPO4S_vr)
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
  call destroy(ZDRAIN_col)
  call destroy(PDRAIN_col)
  call destroy(UION_col)
  call destroy(HydroIonFlx_CumYr_col)
  call destroy(Micb_N2Fixation_vr)
  call destroy(RNut_MicbRelease_vr)
  call destroy(REcoDOMProd_vr)
  call destroy(RDOMMicProd_vr)
  call destroy(TMicHeterActivity_vr)
  call destroy(VWatMicrobAct_vr)
  call destroy(TSens4MicbGrwoth_vr)
  call destroy(LitrfalStrutElms_vr)
  call destroy(trcg_TotalMass_beg_col)
  call destroy(trcg_TotalMass_col)
  call destroy(SurfGasDifFlx_col)
  call destroy(SurfGasEmisFlx_col)
  call destroy(GasHydroLossFlx_col)
  call destroy(trcs_VLN_vr)
  call destroy(trcg_ebu_flx_vr)
  call destroy(trcg_DisolEvap_Atm2Litr_flx)
  call destroy(trcg_soilMass_beg_col)
  call destroy(trcg_soilMass_col)
  call destroy(BandWidthNH4_vr)
  call destroy(BandThicknessNH4_vr)
  call destroy(BandWidthNO3_vr)
  call destroy(TempSensDecomp_vr)
  call destroy(MoistSensDecomp_vr)
  call destroy(BandThicknessNO3_vr)
  call destroy(BandWidthPO4_vr)
  call destroy(BandThicknessPO4_vr)
  call destroy(BandDepthNH4_col)
  call destroy(BandDepthNO3_col)
  call destroy(BandDepthPO4_col)
  call destroy(RNO2DmndSoilChemo_vr)
  call destroy(RNO2DmndBandChemo_vr)
  call destroy(RProd_Hp_vr)
  call destroy(WaterFlowSoiMicP_3D)
  call destroy(WaterFlowSoiMacP_3D)
  call destroy(HeatFlow2Soil_3D)
  call destroy(DOM_MicpTransp_3D)
  call destroy(DOM_Macp_Transp_flx_3D)
  call destroy(trcs_RMicbUptake_vr)
  end subroutine DestructSoilBGCData

end module SoilBGCDataType
