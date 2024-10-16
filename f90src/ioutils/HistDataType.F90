module HistDataType
!
! this module is an intermediate step to support ascii output
! when output is done with netcdf, no id is needed.
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use data_const_mod,   only: spval  => DAT_CONST_SPVAL, ispval => DAT_CONST_ISPVAL
  use SoilBGCNLayMod,   only: SumMicbGroup,              sumDOML, sumMicBiomLayL
  use UnitMod,          only: units
  use MiniMathMod,      only: safe_adb,                  AZMAX1
  use EcoSiMParDataMod, only: pltpar,                    micpar
  use GridConsts
  use GridMod
  use HistFileMod
  use ElmIDMod
  use GridDataType
  use EcoSIMCtrlDataType
  use EcosimConst
  use EcoSIMHistMod
  use FlagDataType
  use PlantTraitDataType
  use PlantDataRateType
  use ClimForcDataType
  use CanopyDataType
  use RootDataType
  use SOMDataType
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use SnowDataType
  use ChemTranspDataType
  use PlantMgmtDataType
  use EcosimBGCFluxType
  use SoilPropertyDataType
  use SoilBGCDataType
  use AqueChemDatatype
  use SurfSoilDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  type, public :: histdata_type
  real(r8),pointer   :: h1D_tFIRE_CO2_col(:)     !CO2byFire_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tFIRE_CH4_col(:)     !CH4byFire_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_cNH4_LITR_col(:)     !(trc_solml_vr(ids_NH4,0,NY,NX)+14.0*trcx_solml_vr(idx_NH4,0,NY,NX))/VLSoilMicPMass_vr(0,NY,NX)
  real(r8),pointer   :: h1D_cNO3_LITR_col(:)      !(trc_solml_vr(ids_NO3,0,NY,NX)+trc_solml_vr(ids_NO2,0,NY,NX))/VLSoilMicPMass_vr(0,NY,NX)                            
  real(r8),pointer   :: h1D_ECO_HVST_N_col(:)     !EcoHavstElmnt_CumYr_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_NET_N_MIN_col(:)      !-NetNH4Mineralize_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tLITR_P_col(:)    !tLitrOM_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_HUMUS_C_col(:)        !tHumOM_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_HUMUS_N_col(:)        !tHumOM_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_HUMUS_P_col(:)        !tHumOM_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_AMENDED_P_col(:)       !FerPFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tLITRf_P_FLX_col(:)   !LiterfalOrgM_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tEXCH_PO4_col(:)       !tHxPO4_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX), exchangeable 
  real(r8),pointer   :: h1D_SUR_DOP_FLX_col(:)    !HydroSufDOPFlx_CumYr_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_SUB_DOP_FLX_col(:)    !HydroSubsDOPFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_SUR_DIP_FLX_col(:)    !HydroSufDIPFlx_CumYr_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_SUB_DIP_FLX_col(:)    !HydroSubsDIPFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_HeatFlx2Grnd_col(:)   !

  real(r8),pointer   :: h1D_Qinfl2soi_col(:)      !
  real(r8),pointer   :: h1D_Qdrain_col(:)          !drainage
  real(r8),pointer   :: h1D_tPRECIP_P_col(:)       !tXPO4_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tMICRO_P_col(:)        !tMicBiome_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tSoilOrgC_col(:)
  real(r8),pointer   :: h1D_tSoilOrgN_col(:)
  real(r8),pointer   :: h1D_tSoilOrgP_col(:)
  real(r8),pointer   :: h1D_PO4_FIRE_col(:)       !PO4byFire_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_cPO4_LITR_col(:)      !trc_solml_vr(ids_H2PO4,0,NY,NX)/VLSoilMicPMass_vr(0,NY,NX)
  real(r8),pointer   :: h1D_cEXCH_P_LITR_col(:)     !31.0*(trcx_solml_vr(idx_HPO4,0,NY,NX)+trcx_solml_vr(idx_H2PO4,0,NY,NX))/VLSoilMicPMass_vr(0,NY,NX)
  real(r8),pointer   :: h1D_ECO_HVST_P_col(:)     !EcoHavstElmnt_CumYr_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_NET_P_MIN_col(:)      !-NetPO4Mineralize_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tSALT_DISCHG_FLX_col(:)    !HydroIonFlx_CumYr_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_PSI_SURF_col(:)       !PSISM(0,NY,NX)
  real(r8),pointer   :: h1D_SURF_ELEV_col(:)      !-CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)+DLYR(3,0,NY,NX)
  real(r8),pointer   :: h1D_tLITR_N_col(:)      !tLitrOM_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h2D_BotDEPZ_vr(:,:)
  real(r8),pointer   :: h1D_tRAD_col(:)
  real(r8),pointer   :: h1D_tHeatUptk_col(:)
  real(r8),pointer   :: h1D_AMENDED_N_col(:)       !FertNFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tLITRf_N_FLX_col(:)  !LiterfalOrgM_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tLITRf_C_FLX_col(:)  !LiterfalOrgM_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tNH4X_col(:)           !tNH4_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX), total NH3+NH4 content
  real(r8),pointer   :: h1D_tNO3_col(:)           !tNO3_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX), total NO3+NO2 content
  real(r8),pointer   :: h1D_SUR_DON_FLX_col(:)    !HydroSufDONFlx_CumYr_col(NY,NX)/TAREA, daily flux
  real(r8),pointer   :: h1D_SUB_DON_FLX_col(:)    !HydroSubsDONFlx_col(NY,NX)/TAREA, daily flux
  real(r8),pointer   :: h1D_SUR_DIN_FLX_col(:)    !HydroSufDINFlx_CumYr_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_SUB_DIN_FLX_col(:)    !HydroSubsDINFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_tMICRO_N_col(:)        !tMicBiome_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_TEMP_LITR_col(:)      !TCS(0,NY,NX)
  real(r8),pointer   :: h1D_TEMP_SNOW_col(:)      !TCSnow_snvr(1,NY,NX)
  real(r8),pointer   :: h1D_tLITR_C_col(:)      !tLitrOM_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_AMENDED_C_col(:)      !AmendCFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tMICRO_C_col(:)        !tMicBiome_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_OMC_LITR_col(:)       !SoilOrgM_vr(ielmc,0,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total residual C
  real(r8),pointer   :: h1D_OMN_LITR_col(:)       !SoilOrgM_vr(ielmn,0,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total residual N
  real(r8),pointer   :: h1D_OMP_LITR_col(:)       !SoilOrgM_vr(ielmp,0,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total residual P
  real(r8),pointer   :: h1D_SUR_DOC_FLX_col(:)    !HydroSufDOCFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_SUB_DOC_FLX_col(:)    !HydroSubsDOCFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_SUR_DIC_FLX_col(:)    !HydroSufDICFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_SUB_DIC_FLX_col(:)    !HydroSubsDICFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: h1D_ATM_CO2_col(:)        !CO2E_col(NY,NX)
  real(r8),pointer   :: h1D_ATM_CH4_col(:)        !CH4E
  real(r8),pointer   :: h1D_NBP_col(:)            !Eco_NBP_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ECO_HVST_C_col(:)     !EcoHavstElmnt_CumYr_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ECO_LAI_col(:)        !CanopyLeafArea_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_Eco_GPP_CumYr_col(:)        !Eco_GPP_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ECO_RA_col(:)         !Eco_AutoR_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_Eco_NPP_CumYr_col(:)        !Eco_NPP_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_Eco_HR_CumYr_col(:)         !Eco_HR_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tDIC_col(:)        !DIC_mass_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX), total soil DIC
  real(r8),pointer   :: h1D_tSTANDING_DEAD_C_col(:)       !StandingDeadStrutElms_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tSTANDING_DEAD_N_col(:)       !StandingDeadStrutElms_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_tSTANDING_DEAD_P_col(:)       !StandingDeadStrutElms_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)    
  real(r8),pointer   :: h1D_tPRECN_col(:)          !1000.0_r8*QRain_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ECO_ET_col(:)             !1000.0_r8*QEvap_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_N2O_LITR_col(:)       !trc_solcl_vr(idg_N2O,0,NY,NX)
  real(r8),pointer   :: h1D_NH3_LITR_col(:)       !trc_solcl_vr(idg_NH3,0,NY,NX)
  real(r8),pointer   :: h1D_SOL_RADN_col(:)       !RAD(NY,NX)*277.8, W m-2
  real(r8),pointer   :: h1D_AIR_TEMP_col(:)       !TCA_col(NY,NX)
  real(r8),pointer   :: h1D_PATM_col(:)           !atmospheric pressure
  real(r8),pointer   :: h1D_HUM_col(:)            !VPK_col(NY,NX)
  real(r8),pointer   :: h1D_WIND_col(:)           !WindSpeedAtm_col(NY,NX)/secs1hour
  real(r8),pointer   :: h1D_PREC_col(:)           !(RainFalPrec(NY,NX)+SnoFalPrec_col(NY,NX))*1000.0/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_SOIL_RN_col(:)        !HeatByRadiation_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX) 
  real(r8),pointer   :: h1D_SOIL_LE_col(:)        !HeatEvapAir2Surf_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_SOIL_H_col(:)         !HeatSensAir2Surf_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_SOIL_G_col(:)         !-(HeatNet2Surf_col(NY,NX)-HeatSensVapAir2Surf_col(NY,NX))*MJ2W/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ECO_RN_col(:)         !Eco_NetRad_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ECO_LE_col(:)         !Eco_Heat_Latent_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_Eco_HeatSen_col(:)    !Eco_Heat_Sens_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ECO_Heat2G_col(:)          !Eco_Heat_Grnd_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_O2_LITR_col(:)       !trc_solcl_vr(idg_O2,0,NY,NX)
  real(r8),pointer   :: h1D_MIN_LWP_ptc(:)       !PSICanPDailyMin(NZ,NY,NX), minimum daily canopy water potential, [MPa]
  real(r8),pointer   :: h1D_SOIL_CO2_FLX_col(:)  !SurfGasFlx_col(idg_CO2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815, umol m-2 s-1, 1.e6/(12*3600)=23.14815
  real(r8),pointer   :: h1D_ECO_CO2_FLX_col(:)   !Eco_NEE_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815
  real(r8),pointer   :: h1D_CH4_FLX_col(:)       !SurfGasFlx_col(idg_CH4,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815, umol m-2 s-1, 1.e6/(12*3600)=23.14815
  real(r8),pointer   :: h1D_CH4_EBU_flx_col(:)
  real(r8),pointer   :: h1D_CH4_PLTROOT_flx_col(:)
  real(r8),pointer   :: h1D_O2_FLX_col(:)        !SurfGasFlx_col(idg_O2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*8.68056,  umol m-2 s-1, 1.e6/(32*3600)=8.68056
  real(r8),pointer   :: h1D_CO2_LITR_col(:)      !trc_solcl_vr(idg_CO2,0,NY,NX)
  real(r8),pointer   :: h1D_EVAPN_col(:)          !VapXAir2GSurf_col(NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_CANET_col(:)         !canopy evaportranspiration
  real(r8),pointer   :: h1D_RUNOFF_FLX_col(:)         !-QRunSurf_col(NY,NX)*1000.0/TAREA, 
  real(r8),pointer   :: h1D_SEDIMENT_FLX_col(:)       !SedmErossLoss_CumYr_col(NY,NX)*1000.0/TAREA, soil mass 
  real(r8),pointer   :: h1D_tSWC_col(:)        !WatMass_col(NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX), volumetric soil water content
  real(r8),pointer   :: h1D_tHeat_col(:) 
  real(r8),pointer   :: h1D_QDISCHG_FLX_col(:)         !QDischar_col(NY,NX)*1000.0/TAREA
  real(r8),pointer   :: h1D_HeatDISCHG_FLX_col(:)
  real(r8),pointer   :: h1D_SNOWPACK_col(:)       !AZMAX1((VOLSS(NY,NX)+VcumIceSnow_col(NY,NX)*DENSICE+VOLWS(NY,NX))*1000.0/AREA(3,NU(NY,NX),NY,NX))
  real(r8),pointer   :: h1D_SURF_WTR_col(:)       !ThetaH2OZ_vr(0,NY,NX)
  real(r8),pointer   :: h1D_SURF_ICE_col(:)       !ThetaICEZ_vr(0,NY,NX)
  real(r8),pointer   :: h1D_ACTV_LYR_col(:)       !-(ActiveLayDepth(NY,NX)-CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX))
  real(r8),pointer   :: h1D_WTR_TBL_col(:)        !-(DepthInternalWTBL(NY,NX)-CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX))
  real(r8),pointer   :: h1D_sN2O_FLX_col(:)        !SurfGasFlx_col(idg_N2O,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_sN2G_FLX_col(:)        !SurfGasFlx_col(idg_N2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_sNH3_FLX_col(:)        !SurfGasFlx_col(idg_NH3,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_frcPARabs_ptc(:)      !fraction of PAR absorbed
  real(r8),pointer   :: h1D_PAR_CAN_ptc(:)        !PAR absorbed by Canopy, umol /m2/s
  real(r8),pointer   :: h1D_PAR_col(:)            !incoming PAR, umol/s
  real(r8),pointer   :: h1D_Plant_C_ptc(:)        !whole plant C  
  real(r8),pointer   :: h1D_Plant_N_ptc(:)        !whole plant N  
  real(r8),pointer   :: h1D_Plant_P_ptc(:)        !whole plant P  
  real(r8),pointer   :: h1D_stomatal_stress_ptc(:)
  real(r8),pointer   :: h1D_CANDew_ptc(:)
  real(r8),pointer   :: h1D_VHeatCap_litr_col(:)
  real(r8),pointer   :: h1D_LEAF_PC_ptc(:)       !(LeafStrutElms_pft(ielmp,NZ,NY,NX)+CanopyNonstElms_pft(ielmp,NZ,NY,NX))/(LeafStrutElms_pft(ielmc,NZ,NY,NX)+CanopyNonstElms_pft(ielmc,NZ,NY,NX)),mass based CP ratio of leaf
  real(r8),pointer   :: h2D_tSOC_vr(:,:)        !SoilOrgM_vr(ielmc,1:JZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total soil C
  real(r8),pointer   :: h2D_litrC_vr(:,:)
  real(r8),pointer   :: h2D_litrN_vr(:,:)
  real(r8),pointer   :: h2D_litrP_vr(:,:)
  real(r8),pointer   :: h2D_tSON_vr(:,:)
  real(r8),pointer   :: h2D_tSOP_vr(:,:)  
  real(r8),pointer   :: h2D_NO3_vr(:,:)
  real(r8),pointer   :: h2D_NH4_vr(:,:)
  real(r8),pointer   :: h2D_VHeatCap_vr(:,:)  
  real(r8),pointer   :: h2D_DOC_vr(:,:)
  real(r8),pointer   :: h2D_DON_vr(:,:)
  real(r8),pointer   :: h2D_DOP_vr(:,:)
  real(r8),pointer   :: h2D_acetate_vr(:,:)
  real(r8),pointer   :: h1D_DOC_litr_col(:)
  real(r8),pointer   :: h1D_DON_litr_col(:)
  real(r8),pointer   :: h1D_DOP_litr_col(:)
  real(r8),pointer   :: h1D_acetate_litr_col(:)
  real(r8),pointer   :: h1D_CAN_RN_ptc(:)        !277.8*RadNet2Canopy_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), W m-2
  real(r8),pointer   :: h1D_CAN_LE_ptc(:)        !277.8*EvapTransHeat_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_CAN_H_ptc(:)         !277.8*HeatXAir2PCan_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_CAN_G_ptc(:)         !277.8*HeatStorCanopy_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_CAN_TEMPC_ptc(:)      !TCelciusCanopy_pft(NZ,NY,NX)
  real(r8),pointer   :: h1D_CAN_TEMPFN_ptc(:)       !fTCanopyGroth_pft(NZ,NY,NX), canopy temperature growth function/stress
  real(r8),pointer   :: h1D_CAN_CO2_FLX_ptc(:)   !CO2NetFix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.148, umol m-2 s-1
  real(r8),pointer   :: h1D_CAN_GPP_ptc(:)       !GrossCO2Fix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX),  gross CO2 fixation, gC m-2/hr
  real(r8),pointer   :: h1D_CAN_RA_ptc(:)        !CanopyRespC_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total autotrophic respiration
  real(r8),pointer   :: h1D_cTNC_ptc(:)          !CanopyNonstElmConc_pft(ielmc,NZ,NY,NX), canopy nonstructural C concentration, 
  real(r8),pointer   :: h1D_cTNN_ptc(:)          !CanopyNonstElmConc_pft(ielmn,NZ,NY,NX)
  real(r8),pointer   :: h1D_cTNP_ptc(:)          !CanopyNonstElmConc_pft(ielmp,NZ,NY,NX)
  real(r8),pointer   :: h1D_STOML_RSC_CO2_ptc(:) !CanPStomaResistH2O_pft(NZ,NY,NX)*1.56*secs1hour, s m-1, for CO2
  real(r8),pointer   :: h1D_BLYR_RSC_CO2_ptc(:)  !CanopyBndlResist_pft(NZ,NY,NX)*1.34*secs1hour, s m-1, for CO2
  real(r8),pointer   :: h1D_CAN_CO2_ptc(:)       !CanopyGasCO2_pft(NZ,NY,NX), umol mol-1
  real(r8),pointer   :: h1D_LAI_ptc(:)           !LeafStalkArea_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant leaf area, include stalk
  real(r8),pointer   :: h1D_PSI_CAN_ptc(:)       !PSICanopy_pft(NZ,NY,NX), canopy total water potential , MPa
  real(r8),pointer   :: h1D_TURG_CAN_ptc(:)      !PSICanopyTurg_pft(NZ,NY,NX), canopy turgor water potential, MPa
  real(r8),pointer   :: h1D_STOM_RSC_H2O_ptc(:)  !CanPStomaResistH2O_pft(NZ,NY,NX)*secs1hour, s m-1, for H2O
  real(r8),pointer   :: h1D_BLYR_RSC_H2O_ptc(:)  !CanopyBndlResist_pft(NZ,NY,NX)*secs1hour, s m-1, for H2O
  real(r8),pointer   :: h1D_TRANSPN_ptc(:)       !Transpiration_pft(NZ,NY,NX)*1000.0_r8/AREA(3,NU(NY,NX),NY,NX), canopy transpiration mm H2O/m2/h
  real(r8),pointer   :: h1D_NH4_UPTK_FLX_ptc(:)      !RootNH4Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_NO3_UPTK_FLX_ptc(:)      !RootNO3Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_N2_FIXN_FLX_ptc(:)       !RootN2Fix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_cNH3_FLX_ptc(:)       !NH3Dep2Can_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_PO4_UPTK_FLX_ptc(:)      !RootH2PO4Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h2D_Root1stStrutC_ptc(:,:) !primary root structural C biomass
  real(r8),pointer   :: h2D_Root1stStrutN_ptc(:,:)
  real(r8),pointer   :: h2D_Root1stStrutP_ptc(:,:)
  real(r8),pointer   :: h2D_Root2ndStrutC_ptc(:,:)
  real(r8),pointer   :: h2D_Root2ndStrutN_ptc(:,:)
  real(r8),pointer   :: h2D_Root2ndStrutP_ptc(:,:)
  real(r8),pointer   :: h1D_SHOOT_C_ptc(:)       !ShootStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_SHOOTST_C_ptc(:)       
  real(r8),pointer   :: h1D_SHOOTST_N_ptc(:)       
  real(r8),pointer   :: h1D_SHOOTST_P_ptc(:)             
  real(r8),pointer   :: h1D_LEAF_C_ptc(:)        !LeafStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_Petole_C_ptc(:)        !PetoleStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), canopy sheath element
  real(r8),pointer   :: h1D_STALK_C_ptc(:)       !StalkStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_RESERVE_C_ptc(:)     !StalkRsrvElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_HUSK_C_ptc(:)        !(HuskStrutElms_pft(ielmc,NZ,NY,NX)+EarStrutElms_pft(ielmc,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_GRAIN_C_ptc(:)       !GrainStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ROOT_NONSTC_ptc(:)
  real(r8),pointer   :: h1D_ROOT_NONSTN_ptc(:)
  real(r8),pointer   :: h1D_ROOT_NONSTP_ptc(:)
  real(r8),pointer   :: h1D_ROOT_C_ptc(:)        !RootElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ROOTST_C_ptc(:)     
  real(r8),pointer   :: h1D_ROOTST_N_ptc(:) 
  real(r8),pointer   :: h1D_ROOTST_P_ptc(:)   
  real(r8),pointer   :: h1D_NODULE_C_ptc(:)      !NodulStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), nodule
  real(r8),pointer   :: h1D_STORED_C_ptc(:)      !NonStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_GRAIN_NO_ptc(:)      !CanopySeedNum_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_LAIb_ptc(:)          !CanopyLeafArea_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total branch leaf area
  real(r8),pointer   :: h1D_EXUD_CumYr_C_FLX_ptc(:)        !PlantExudElm_CumYr_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_LITRf_C_FLX_ptc(:)       !LitrfalStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_LITRf_P_FLX_ptc(:)       !LitrfalStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_SURF_LITRf_C_FLX_ptc(:)  !SurfLitrfalStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_AUTO_RESP_FLX_ptc(:)     !GrossResp_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_HVST_C_FLX_ptc(:)        !EcoHavstElmnt_CumYr_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_PLANT_BALANCE_C_ptc(:)     !ElmBalanceCum_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_STANDING_DEAD_C_ptc(:)    !StandDeadStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_FIREp_CO2_FLX_ptc(:)     !CO2ByFire_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant CO2 from fire
  real(r8),pointer   :: h1D_FIREp_CH4_FLX_ptc(:)     !CH4ByFire_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_NPP_ptc(:)           !NetPrimProduct_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_CAN_HT_ptc(:)        !CanopyHeight_pft(NZ,NY,NX), canopy height, m
  real(r8),pointer   :: h1D_POPN_ptc(:)          !PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant population
  real(r8),pointer   :: h1D_tTRANSPN_ptc(:)      !-ETCanopy_CumYr_pft(NZ,NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX), total transpiration
  real(r8),pointer   :: h1D_WTR_STRESS_ptc(:)    !HoursTooLowPsiCan_pft(NZ,NY,NX)
  real(r8),pointer   :: h1D_OXY_STRESS_ptc(:)    !OSTR(NZ,NY,NX)
  real(r8),pointer   :: h1D_SHOOT_N_ptc(:)       !ShootStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_LEAF_N_ptc(:)        !LeafStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_Petole_N_ptc(:)        !PetoleStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_STALK_N_ptc(:)       !StalkStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_RESERVE_N_ptc(:)     !StalkRsrvElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_HUSK_N_ptc(:)        !(HuskStrutElms_pft(ielmn,NZ,NY,NX)+EarStrutElms_pft(ielmn,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_GRAIN_N_ptc(:)       !GrainStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ROOT_N_ptc(:)        !RootElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_NODULE_N_ptc(:)         !NodulStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_STORED_N_ptc(:)      !NonStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_EXUD_N_FLX_ptc(:)        !PlantExudElm_CumYr_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_Uptk_N_Flx_ptc(:)
  real(r8),pointer   :: h1D_Uptk_P_Flx_ptc(:)
  real(r8),pointer   :: h1D_LITRf_N_FLX_ptc(:)       !LitrfalStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total plant LitrFall N
  real(r8),pointer   :: h1D_TL_N_FIXED_FLX_ptc(:)    !PlantN2Fix_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total plant N2 fixation
  real(r8),pointer   :: h1D_HVST_N_FLX_ptc(:)        !EcoHavstElmnt_CumYr_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_NH3can_FLX_ptc(:)    !NH3Emis_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_PLANT_BALANCE_N_ptc(:)     !ElmBalanceCum_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_STANDING_DEAD_N_ptc(:)    !StandDeadStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_FIREp_N_FLX_ptc(:)        !NH3byFire_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant N emission from fire
  real(r8),pointer   :: h1D_SURF_LITRf_N_FLX_ptc(:)   !SurfLitrfalStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), surface litter fall
  real(r8),pointer   :: h1D_SHOOT_P_ptc(:)       !ShootStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_LEAF_P_ptc(:)        !LeafStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_Petole_P_ptc(:)        !PetoleStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_STALK_P_ptc(:)       !StalkStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_RESERVE_P_ptc(:)     !StalkRsrvElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_HUSK_P_ptc(:)        !(HuskStrutElms_pft(ielmp,NZ,NY,NX)+EarStrutElms_pft(ielmp,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_GRAIN_P_ptc(:)       !GrainStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_ROOT_P_ptc(:)        !RootElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_NODULE_P_ptc(:)         !NodulStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_STORED_P_ptc(:)      !NonStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_EXUD_P_FLX_ptc(:)        !PlantExudElm_CumYr_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_LITTERf_P_ptc(:)     !LitrfalStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_HVST_P_FLX_ptc(:)        !EcoHavstElmnt_CumYr_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_PLANT_BALANCE_P_ptc(:)     !ElmBalanceCum_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_STANDING_DEAD_P_ptc(:)    !StandDeadStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_FIREp_P_FLX_ptc(:)        !PO4byFire_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_SURF_LITRf_P_FLX_ptc(:)  !SurfLitrfalStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: h1D_BRANCH_NO_ptc(:)     !NumOfBranches_pft(NZ,NY,NX)
  real(r8),pointer   :: h1D_SHOOT_NONSTC_ptc(:)   !CanopyNonstElms_pft(ielmc,NY,NX)
  real(r8),pointer   :: h1D_SHOOT_NONSTN_ptc(:)   !
  real(r8),pointer   :: h1D_SHOOT_NONSTP_ptc(:)   !
  real(r8),pointer   :: h1D_CFIX_lmtf_pft(:)
  real(r8),pointer   :: h1D_MainBranchNO_pft(:)
  real(r8),pointer   :: h1D_LEAF_NC_ptc(:)       !(LeafStrutElms_pft(ielmn,NZ,NY,NX)+CanopyNonstElms_pft(ielmn,NZ,NY,NX))/(LeafStrutElms_pft(ielmc,NZ,NY,NX)+CanopyNonstElms_pft(ielmc,NZ,NY,NX)),mass based CN ratio of leaf  
  real(r8),pointer    :: h1D_Growth_Stage_ptc(:)    !plant development stage, integer, 0-10, planting, emergence, floral_init, jointing, 
                                      !elongation, heading, anthesis, seed_fill, see_no_set, seed_mass_set, end_seed_fill
  real(r8),pointer   :: h2D_LEAF_NODE_NO_ptc(:,:)       !NumOfLeaves_brch(MainBranchNum_pft(NZ,NY,NX),NZ,NY,NX), leaf NO
  real(r8),pointer   :: h2D_RUB_ACTVN_ptc(:,:)     !RubiscoActivity_brch(MainBranchNum_pft(NZ,NY,NX),NZ,NY,NX), branch down-regulation of CO2 fixation
  real(r8),pointer   :: h3D_PARTS_ptc(:,:,:)       !

  real(r8),pointer   :: h2D_AeroHrBactC_vr(:,:)     !aerobic heterotropic bacteria
  real(r8),pointer   :: h2D_AeroHrFungC_vr(:,:)   !aerobic heterotropic fungi
  real(r8),pointer   :: h2D_faculDenitC_vr(:,:)   !facultative denitrifier
  real(r8),pointer   :: h2D_fermentorC_vr(:,:)  !fermentor
  real(r8),pointer   :: h2D_acetometgC_vr(:,:)  !acetogenic methanogen
  real(r8),pointer   :: h2D_aeroN2fixC_vr(:,:)  !aerobic N2 fixer
  real(r8),pointer   :: h2D_anaeN2FixC_vr(:,:)  !anaerobic N2 fixer
  real(r8),pointer   :: h2D_NH3OxiBactC_vr(:,:)
  real(r8),pointer   :: h2D_NO2OxiBactC_vr(:,:)
  real(r8),pointer   :: h2D_CH4AeroOxiC_vr(:,:)
  real(r8),pointer   :: h2D_H2MethogenC_vr(:,:)

  real(r8),pointer   :: h2D_RCH4ProdHydrog_vr(:,:)
  real(r8),pointer   :: h2D_RCH4ProdAcetcl_vr(:,:)
  real(r8),pointer   :: h2D_RCH4Oxi_aero_vr(:,:)
  real(r8),pointer   :: h2D_RFerment_vr(:,:)
  real(r8),pointer   :: h2D_nh3oxi_vr(:,:)
  real(r8),pointer   :: h2D_n2oprod_vr(:,:)
  real(r8),pointer   :: h2D_RootMassC_vr(:,:)
  real(r8),pointer   :: h2D_RootMassN_vr(:,:)
  real(r8),pointer   :: h2D_RootMassP_vr(:,:)
  real(r8),pointer   :: h1D_RCH4ProdHydrog_litr_col(:)
  real(r8),pointer   :: h1D_RCH4ProdAcetcl_litr_col(:)  
  real(r8),pointer   :: h1D_RCH4Oxi_aero_litr_col(:)
  real(r8),pointer   :: h1D_RFermen_litr_col(:)
  real(r8),pointer   :: h1D_NH3oxi_litr_col(:)
  real(r8),pointer   :: h1D_N2oprod_litr_col(:)

  real(r8),pointer   :: h2D_AeroHrBactN_vr(:,:)     !aerobic heterotropic bacteria
  real(r8),pointer   :: h2D_AeroHrFungN_vr(:,:)   !aerobic heterotropic fungi
  real(r8),pointer   :: h2D_faculDenitN_vr(:,:)   !facultative denitrifier
  real(r8),pointer   :: h2D_fermentorN_vr(:,:)  !fermentor
  real(r8),pointer   :: h2D_acetometgN_vr(:,:)  !acetogenic methanogen
  real(r8),pointer   :: h2D_aeroN2fixN_vr(:,:)  !aerobic N2 fixer
  real(r8),pointer   :: h2D_anaeN2FixN_vr(:,:)  !anaerobic N2 fixer
  real(r8),pointer   :: h2D_NH3OxiBactN_vr(:,:)
  real(r8),pointer   :: h2D_NO2OxiBactN_vr(:,:)
  real(r8),pointer   :: h2D_CH4AeroOxiN_vr(:,:)
  real(r8),pointer   :: h2D_H2MethogenN_vr(:,:)

  real(r8),pointer   :: h2D_AeroHrBactP_vr(:,:)     !aerobic heterotropic bacteria
  real(r8),pointer   :: h2D_AeroHrFungP_vr(:,:)   !aerobic heterotropic fungi
  real(r8),pointer   :: h2D_faculDenitP_vr(:,:)   !facultative denitrifier
  real(r8),pointer   :: h2D_fermentorP_vr(:,:)  !fermentor
  real(r8),pointer   :: h2D_acetometgP_vr(:,:)  !acetogenic methanogen
  real(r8),pointer   :: h2D_aeroN2fixP_vr(:,:)  !aerobic N2 fixer
  real(r8),pointer   :: h2D_anaeN2FixP_vr(:,:)  !anaerobic N2 fixer
  real(r8),pointer   :: h2D_NH3OxiBactP_vr(:,:)
  real(r8),pointer   :: h2D_NO2OxiBactP_vr(:,:)
  real(r8),pointer   :: h2D_CH4AeroOxiP_vr(:,:)
  real(r8),pointer   :: h2D_H2MethogenP_vr(:,:)
  real(r8),pointer   :: h2D_MicroBiomeE_litr_col(:,:)   !total microbial biomass in litr
  real(r8),pointer   :: h2D_AeroHrBactE_litr_col(:,:)     !aerobic heterotropic bacteria
  real(r8),pointer   :: h2D_AeroHrFungE_litr_col(:,:)   !aerobic heterotropic fungi
  real(r8),pointer   :: h2D_faculDenitE_litr_col(:,:)   !facultative denitrifier
  real(r8),pointer   :: h2D_fermentorE_litr_col(:,:)  !fermentor
  real(r8),pointer   :: h2D_acetometgE_litr_col(:,:)  !acetogenic methanogen
  real(r8),pointer   :: h2D_aeroN2fixE_litr_col(:,:)  !aerobic N2 fixer
  real(r8),pointer   :: h2D_anaeN2FixE_litr_col(:,:)  !anaerobic N2 fixer
  real(r8),pointer   :: h2D_NH3OxiBactE_litr_col(:,:)
  real(r8),pointer   :: h2D_NO2OxiBactE_litr_col(:,:)
  real(r8),pointer   :: h2D_CH4AeroOxiE_litr_col(:,:)
  real(r8),pointer   :: h2D_H2MethogenE_litr_col(:,:)

  real(r8),pointer   :: h2D_RNITRIF_vr(:,:)
  real(r8),pointer   :: h2D_CO2_vr(:,:)        !trc_solcl_vr(idg_CO2,1:JZ,NY,NX)
  real(r8),pointer   :: h2D_CH4_vr(:,:)        !trc_solcl_vr(idg_CH4,1:JZ,NY,NX)
  real(r8),pointer   :: h2D_O2_vr(:,:)         !trc_solcl_vr(idg_O2,1:JZ,NY,NX)
  real(r8),pointer   :: h2D_N2O_vr(:,:)         !trc_solcl_vr(idg_N2O,1:JZ,NY,NX)
  real(r8),pointer   :: h2D_NH3_vr(:,:)         !trc_solcl_vr(idg_NH3,1:JZ,NY,NX)
  real(r8),pointer   :: h2D_TEMP_vr(:,:)        !TCS(1:JZ,NY,NX)
  real(r8),pointer   :: h2D_HeatUptk_vr(:,:)    !Heat uptake by root  
  real(r8),pointer   :: h2D_HeatFlow_vr(:,:)    !
  real(r8),pointer   :: h2D_VSPore_vr(:,:)
  real(r8),pointer   :: h2D_VSM_vr    (:,:)       !ThetaH2OZ_vr(1:JZ,NY,NX)
  real(r8),pointer   :: h2D_VSICE_vr    (:,:)         !ThetaICEZ_vr(1:JZ,NY,NX)  
  real(r8),pointer   :: h2D_PSI_vr(:,:)         !PSISM(1:JZ,NY,NX)+PSISO(1:JZ,NY,NX)
  real(r8),pointer   :: h2D_PsiO_vr(:,:)
  real(r8),pointer   :: h2D_cNH4t_vr(:,:)       !(trc_solml_vr(ids_NH4,1:JZ,NY,NX)+trc_solml_vr(ids_NH4B,1:JZ,NY,NX) &
                                                                  !+14.0*(trcx_solml_vr(idx_NH4,1:JZ,NY,NX)+trcx_solml_vr(idx_NH4B,1:JZ,NY,NX)))/VLSoilMicPMass_vr(1:JZ,NY,NX)
  real(r8),pointer   :: h2D_RootH2OUP_vr(:,:)   !root water uptake flux                                 
  real(r8),pointer   :: h2D_cNO3t_vr(:,:)       !(trc_solml_vr(ids_NO3,1:JZ,NY,NX)+trc_solml_vr(ids_NO3B,1:JZ,NY,NX) &
                                                                  !+trc_solml_vr(ids_NO2,1,NY,NX)+trc_solml_vr(ids_NO2B,1,NY,NX))/VLSoilMicPMass_vr(1,NY,NX)
  real(r8),pointer   :: h2D_cPO4_vr(:,:)        !(trc_solml_vr(ids_H1PO4,1:JZ,NY,NX)+trc_solml_vr(ids_H1PO4B,1,NY,NX)+trc_solml_vr(ids_H2PO4,1,NY,NX)+trc_solml_vr(ids_H2PO4B,1,NY,NX))/VLWatMicP_vr(1,NY,NX)
  real(r8),pointer   :: h2D_cEXCH_P_vr(:,:)     !31.0*(trcx_solml_vr(idx_HPO4,1,NY,NX)+trcx_solml_vr(idx_H2PO4,1,NY,NX)+trcx_solml_vr(idx_HPO4B,1,NY,NX)+trcx_solml_vr(idx_H2PO4B,1,NY,NX))/VLSoilMicPMass_vr(1,NY,NX)
  real(r8),pointer   :: h2D_ECND_vr(:,:)         !ECND_vr(1:JZ,NY,NX)

  real(r8),pointer   :: h2D_PSI_RT_vr_ptc(:,:)     !PSIRoot_pvr(1,1:JZ,NZ,NY,NX), root total water potential , MPa
  real(r8),pointer   :: h2D_prtUP_NH4_vr_ptc(:,:)     !(RootNutUptake_pvr(ids_NH4,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NH4,2,1:JZ,NZ,NY,NX) &
                                                                   !+RootNutUptake_pvr(ids_NH4B,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NH4B,2,1:JZ,NZ,NY,NX))/AREA(3,1,NY,NX)
  real(r8),pointer   :: h2D_prtUP_NO3_vr_ptc(:,:)     !(RootNutUptake_pvr(ids_NO3,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NO3,2,1:JZ,NZ,NY,NX) &
                                                                   !+RootNutUptake_pvr(ids_NO3B,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NO3B,2,1:JZ,NZ,NY,NX))/AREA(3,1,NY,NX)
  real(r8),pointer   :: h2D_prtUP_PO4_vr_ptc(:,:)     !(RootNutUptake_pvr(ids_H2PO4,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_H2PO4,2,1:JZ,NZ,NY,NX) &
                                                                   !+RootNutUptake_pvr(ids_H2PO4B,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_H2PO4B,2,1:JZ,NZ,NY,NX))/AREA(3,1,NY,NX)
  real(r8),pointer   :: h2D_DNS_RT_vr_ptc(:,:)     !RootLenDensPerPlant_pvr(1,1:JZ,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  contains
    procedure, public :: Init  => init_hist_data
    procedure, public :: hist_update
  end type histdata_type

  type(histdata_type),public :: hist_ecosim
  contains

  subroutine init_hist_data(this,bounds)

  implicit none
  class(histdata_type) :: this
  type(bounds_type), intent(in) :: bounds
  integer :: beg_col, end_col
  integer :: beg_ptc, end_ptc
  real(r8), pointer :: data1d_ptr(:)
  integer , pointer :: idata1d_ptr(:)
  real(r8), pointer :: data2d_ptr(:,:)
  integer :: nbr
  character(len=32) :: fieldname

  beg_col=1;end_col=bounds%ncols
  beg_ptc=1;end_ptc=bounds%npfts
  allocate(this%h1D_tFIRE_CO2_col(beg_col:end_col)) ;this%h1D_tFIRE_CO2_col(:)=spval    
  allocate(this%h1D_tFIRE_CH4_col(beg_col:end_col)) ;this%h1D_tFIRE_CH4_col(:)=spval    
  allocate(this%h1D_cNH4_LITR_col(beg_col:end_col)) ;this%h1D_cNH4_LITR_col(:)=spval    
  allocate(this%h1D_cNO3_LITR_col(beg_col:end_col)) ;this%h1D_cNO3_LITR_col(:)=spval    
  allocate(this%h1D_ECO_HVST_C_col(beg_col:end_col));this%h1D_ECO_HVST_C_col(:)=spval   
  allocate(this%h1D_ECO_HVST_N_col(beg_col:end_col));this%h1D_ECO_HVST_N_col(:)=spval   
  allocate(this%h1D_ECO_HVST_P_col(beg_col:end_col));this%h1D_ECO_HVST_P_col(:)=spval   
  allocate(this%h1D_NET_N_MIN_col(beg_col:end_col)) ;this%h1D_NET_N_MIN_col(:) =spval   
  allocate(this%h1D_tLITR_P_col(beg_col:end_col)); this%h1D_tLITR_P_col(:)=spval   
  allocate(this%h1D_HUMUS_C_col(beg_col:end_col))  ;this%h1D_HUMUS_C_col(:)=spval 
  allocate(this%h1D_HUMUS_N_col(beg_col:end_col))  ;this%h1D_HUMUS_N_col(:)=spval 
  allocate(this%h1D_HUMUS_P_col(beg_col:end_col))  ;this%h1D_HUMUS_P_col(:)=spval 
  allocate(this%h1D_AMENDED_P_col(beg_col:end_col));this%h1D_AMENDED_P_col(:)=spval 
  allocate(this%h1D_tLITRf_C_FLX_col(beg_col:end_col))  ;this%h1D_tLITRf_C_FLX_col(:)=spval 
  allocate(this%h1D_tLITRf_N_FLX_col(beg_col:end_col))  ;this%h1D_tLITRf_N_FLX_col(:)=spval 
  allocate(this%h1D_tLITRf_P_FLX_col(beg_col:end_col))  ;this%h1D_tLITRf_P_FLX_col(:)=spval 
  allocate(this%h1D_tEXCH_PO4_col(beg_col:end_col))     ;this%h1D_tEXCH_PO4_col(:)=spval 
  allocate(this%h1D_SUR_DOP_FLX_col(beg_col:end_col))   ;this%h1D_SUR_DOP_FLX_col(:)=spval 
  allocate(this%h1D_SUB_DOP_FLX_col(beg_col:end_col))   ;this%h1D_SUB_DOP_FLX_col(:)=spval 
  allocate(this%h1D_SUR_DIP_FLX_col(beg_col:end_col))   ;this%h1D_SUR_DIP_FLX_col(:)=spval
  allocate(this%h1D_SUB_DIP_FLX_col(beg_col:end_col))   ;this%h1D_SUB_DIP_FLX_col(:)=spval
  allocate(this%h1D_HeatFlx2Grnd_col(beg_col:end_col))     ;this%h1D_HeatFlx2Grnd_col(:)=spval
  allocate(this%h1D_Qinfl2soi_col(beg_col:end_col))     ;this%h1D_Qinfl2soi_col(:)=spval
  allocate(this%h1D_Qdrain_col(beg_col:end_col))       ; this%h1D_Qdrain_col(:)=spval
  allocate(this%h1D_tSALT_DISCHG_FLX_col(beg_col:end_col)) ;this%h1D_tSALT_DISCHG_FLX_col(:)=spval
  allocate(this%h1D_SUR_DON_FLX_col(beg_col:end_col))    ;this%h1D_SUR_DON_FLX_col(:)=spval
  allocate(this%h1D_SUB_DON_FLX_col(beg_col:end_col))    ;this%h1D_SUB_DON_FLX_col(:)=spval
  allocate(this%h1D_SUR_DIN_FLX_col(beg_col:end_col))    ;this%h1D_SUR_DIN_FLX_col(:)=spval
  allocate(this%h1D_SUB_DIN_FLX_col(beg_col:end_col))    ;this%h1D_SUB_DIN_FLX_col(:)=spval 
  allocate(this%h1D_SUR_DOC_FLX_col(beg_col:end_col))    ;this%h1D_SUR_DOC_FLX_col(:)=spval
  allocate(this%h1D_SUB_DOC_FLX_col(beg_col:end_col))    ;this%h1D_SUB_DOC_FLX_col(:)=spval 
  allocate(this%h1D_SUR_DIC_FLX_col(beg_col:end_col))    ;this%h1D_SUR_DIC_FLX_col(:)=spval 
  allocate(this%h1D_SUB_DIC_FLX_col(beg_col:end_col))    ;this%h1D_SUB_DIC_FLX_col(:)=spval 

  allocate(this%h1D_tPRECIP_P_col(beg_col:end_col))      ;this%h1D_tPRECIP_P_col(:)=spval 
  allocate(this%h1D_tMICRO_P_col(beg_col:end_col))       ;this%h1D_tMICRO_P_col(:)=spval
  allocate(this%h1D_PO4_FIRE_col(beg_col:end_col))       ;this%h1D_PO4_FIRE_col(:)=spval
  allocate(this%h1D_cPO4_LITR_col(beg_col:end_col))      ;this%h1D_cPO4_LITR_col(:)=spval
  allocate(this%h1D_cEXCH_P_LITR_col(beg_col:end_col))   ;this%h1D_cEXCH_P_LITR_col(:)=spval
  allocate(this%h1D_NET_P_MIN_col(beg_col:end_col))      ;this%h1D_NET_P_MIN_col(:)=spval

  allocate(this%h1D_PSI_SURF_col(beg_col:end_col))       ;this%h1D_PSI_SURF_col(:)=spval
  allocate(this%h1D_SURF_ELEV_col(beg_col:end_col))      ;this%h1D_SURF_ELEV_col(:)=spval
  allocate(this%h1D_tLITR_C_col(beg_col:end_col));this%h1D_tLITR_C_col(:)=spval
  allocate(this%h1D_tLITR_N_col(beg_col:end_col));this%h1D_tLITR_N_col(:)=spval
  allocate(this%h1D_tRAD_col(beg_col:end_col)); this%h1D_tRAD_col(:)=spval
  allocate(this%h2D_BotDEPZ_vr(beg_col:end_col,1:JZ)); this%h2D_BotDEPZ_vr(:,:)=spval
  allocate(this%h1D_AMENDED_N_col(beg_col:end_col))       ;this%h1D_AMENDED_N_col(:)=spval
  allocate(this%h1D_tHeatUptk_col(beg_col:end_col))       ;this%h1D_tHeatUptk_col(:)=spval
  allocate(this%h1D_tNH4X_col(beg_col:end_col))           ;this%h1D_tNH4X_col(:)=spval
  allocate(this%h1D_tNO3_col(beg_col:end_col))            ;this%h1D_tNO3_col(:)=spval
  allocate(this%h1D_tMICRO_N_col(beg_col:end_col))        ;this%h1D_tMICRO_N_col(:)=spval
  allocate(this%h1D_TEMP_LITR_col(beg_col:end_col))       ;this%h1D_TEMP_LITR_col(:)=spval
  allocate(this%h1D_TEMP_SNOW_col(beg_col:end_col))       ;this%h1D_TEMP_SNOW_col(:)=spval

  allocate(this%h1D_AMENDED_C_col(beg_col:end_col))       ;this%h1D_AMENDED_C_col(:)=spval
  allocate(this%h1D_tMICRO_C_col(beg_col:end_col))        ;this%h1D_tMICRO_C_col(:)=spval
  allocate(this%h1D_tSoilOrgC_col(beg_col:end_col))       ;this%h1D_tSoilOrgC_col(:)=spval
  allocate(this%h1D_tSoilOrgN_col(beg_col:end_col))       ;this%h1D_tSoilOrgN_col(:)=spval
  allocate(this%h1D_tSoilOrgP_col(beg_col:end_col))       ;this%h1D_tSoilOrgP_col(:)=spval
  allocate(this%h1D_OMC_LITR_col(beg_col:end_col))        ;this%h1D_OMC_LITR_col(:)=spval
  allocate(this%h1D_OMN_LITR_col(beg_col:end_col))        ;this%h1D_OMN_LITR_col(:)=spval
  allocate(this%h1D_OMP_LITR_col(beg_col:end_col))        ;this%h1D_OMP_LITR_col(:)=spval

  allocate(this%h1D_ATM_CO2_col(beg_col:end_col))         ;this%h1D_ATM_CO2_col(:)=spval
  allocate(this%h1D_ATM_CH4_col(beg_col:end_col))         ;this%h1D_ATM_CH4_col(:)=spval
  allocate(this%h1D_NBP_col(beg_col:end_col))             ;this%h1D_NBP_col(:)=spval

  allocate(this%h1D_ECO_LAI_col(beg_col:end_col))         ;this%h1D_ECO_LAI_col(:)=spval
  allocate(this%h1D_Eco_GPP_CumYr_col(beg_col:end_col))         ;this%h1D_Eco_GPP_CumYr_col(:)=spval
  allocate(this%h1D_ECO_RA_col(beg_col:end_col))          ;this%h1D_ECO_RA_col(:)=spval
  allocate(this%h1D_Eco_NPP_CumYr_col(beg_col:end_col))         ;this%h1D_Eco_NPP_CumYr_col(:)=spval
  allocate(this%h1D_Eco_HR_CumYr_col(beg_col:end_col))          ;this%h1D_Eco_HR_CumYr_col(:)=spval
  allocate(this%h1D_tDIC_col(beg_col:end_col))            ;this%h1D_tDIC_col=spval
  allocate(this%h1D_tSTANDING_DEAD_C_col(beg_col:end_col));this%h1D_tSTANDING_DEAD_C_col=spval 
  allocate(this%h1D_tSTANDING_DEAD_N_col(beg_col:end_col));this%h1D_tSTANDING_DEAD_N_col=spval  
  allocate(this%h1D_tSTANDING_DEAD_P_col(beg_col:end_col));this%h1D_tSTANDING_DEAD_P_col=spval  
  allocate(this%h1D_tPRECN_col(beg_col:end_col))          ;this%h1D_tPRECN_col(:)=spval
  allocate(this%h1D_ECO_ET_col(beg_col:end_col))              ;this%h1D_ECO_ET_col(:)=spval
  allocate(this%h1D_N2O_LITR_col(beg_col:end_col))        ;this%h1D_N2O_LITR_col(:)=spval
  allocate(this%h1D_NH3_LITR_col(beg_col:end_col))        ;this%h1D_NH3_LITR_col(:)=spval
  allocate(this%h1D_SOL_RADN_col(beg_col:end_col))        ;this%h1D_SOL_RADN_col(:)=spval
  allocate(this%h1D_AIR_TEMP_col(beg_col:end_col))        ;this%h1D_AIR_TEMP_col(:)=spval
  allocate(this%h1D_PATM_col(beg_col:end_col))            ;this%h1D_PATM_col(:)=spval
  allocate(this%h1D_HUM_col(beg_col:end_col))             ;this%h1D_HUM_col(:)=spval
  allocate(this%h1D_WIND_col(beg_col:end_col))            ;this%h1D_WIND_col(:)=spval
  allocate(this%h1D_PREC_col(beg_col:end_col))            ;this%h1D_PREC_col(:)=spval
  allocate(this%h1D_SOIL_RN_col(beg_col:end_col))         ;this%h1D_SOIL_RN_col(:)=spval
  allocate(this%h1D_SOIL_LE_col(beg_col:end_col))         ;this%h1D_SOIL_LE_col(:)=spval
  allocate(this%h1D_SOIL_H_col(beg_col:end_col))          ;this%h1D_SOIL_H_col(:)=spval
  allocate(this%h1D_SOIL_G_col(beg_col:end_col))          ;this%h1D_SOIL_G_col(:)=spval
  allocate(this%h1D_ECO_RN_col(beg_col:end_col))          ;this%h1D_ECO_RN_col(:)=spval
  allocate(this%h1D_ECO_LE_col(beg_col:end_col))          ;this%h1D_ECO_LE_col(:)=spval
  allocate(this%h1D_Eco_HeatSen_col(beg_col:end_col))        ;this%h1D_Eco_HeatSen_col(:)=spval
  allocate(this%h1D_ECO_Heat2G_col(beg_col:end_col))           ;this%h1D_ECO_Heat2G_col(:)=spval
  allocate(this%h1D_O2_LITR_col(beg_col:end_col))         ;this%h1D_O2_LITR_col(:)=spval
  allocate(this%h1D_MIN_LWP_ptc(beg_ptc:end_ptc))         ;this%h1D_MIN_LWP_ptc(:)=spval
  allocate(this%h1D_SOIL_CO2_FLX_col(beg_col:end_col))    ;this%h1D_SOIL_CO2_FLX_col(:)=spval
  allocate(this%h1D_ECO_CO2_FLX_col(beg_col:end_col))     ;this%h1D_ECO_CO2_FLX_col(:)=spval
  allocate(this%h1D_CH4_FLX_col(beg_col:end_col))         ;this%h1D_CH4_FLX_col(:)=spval
  allocate(this%h1D_CH4_EBU_flx_col(beg_col:end_col))     ;this%h1D_CH4_EBU_flx_col(:)=spval
  allocate(this%h1D_CH4_PLTROOT_flx_col(beg_col:end_col)) ;this%h1D_CH4_PLTROOT_flx_col(:)=spval
  allocate(this%h1D_O2_FLX_col(beg_col:end_col))          ;this%h1D_O2_FLX_col(:)=spval
  allocate(this%h1D_CO2_LITR_col(beg_col:end_col))        ;this%h1D_CO2_LITR_col(:)=spval
  allocate(this%h1D_EVAPN_col(beg_col:end_col))           ;this%h1D_EVAPN_col(:)=spval
  allocate(this%h1D_CANET_col(beg_col:end_col))           ;this%h1D_CANET_col(:)=spval
  allocate(this%h1D_PAR_col(beg_col:end_col))             ;this%h1D_PAR_col(:)=spval
  allocate(this%h1D_tSWC_col(beg_col:end_col))            ;this%h1D_tSWC_col(:)=spval
  allocate(this%h1D_tHeat_col(beg_col:end_col))           ;this%h1D_tHeat_col(:)=spval
  allocate(this%h1D_SNOWPACK_col(beg_col:end_col))        ;this%h1D_SNOWPACK_col(:)=spval
  allocate(this%h1D_SURF_WTR_col(beg_col:end_col))        ;this%h1D_SURF_WTR_col(:)=spval
  allocate(this%h1D_SURF_ICE_col(beg_col:end_col))        ;this%h1D_SURF_ICE_col(:)=spval
  allocate(this%h1D_ACTV_LYR_col(beg_col:end_col))        ;this%h1D_ACTV_LYR_col(:)=spval
  allocate(this%h1D_WTR_TBL_col(beg_col:end_col))         ;this%h1D_WTR_TBL_col(:)=spval
  allocate(this%h1D_sN2O_FLX_col(beg_col:end_col))        ;this%h1D_sN2O_FLX_col(:)=spval
  allocate(this%h1D_sN2G_FLX_col(beg_col:end_col))        ;this%h1D_sN2G_FLX_col(:)=spval
  allocate(this%h1D_sNH3_FLX_col(beg_col:end_col))        ;this%h1D_sNH3_FLX_col(:)=spval
  allocate(this%h1D_VHeatCap_litr_col(beg_col:end_col))   ;this%h1D_VHeatCap_litr_col(:)=spval
  allocate(this%h1D_RUNOFF_FLX_col(beg_col:end_col))      ;this%h1D_RUNOFF_FLX_col(:)=spval
  allocate(this%h1D_SEDIMENT_FLX_col(beg_col:end_col))    ;this%h1D_SEDIMENT_FLX_col(:)=spval
  allocate(this%h1D_QDISCHG_FLX_col(beg_col:end_col))      ;this%h1D_QDISCHG_FLX_col(:)=spval
  allocate(this%h1D_HeatDISCHG_FLX_col(beg_col:end_col))  ; this%h1D_HeatDISCHG_FLX_col(:)=spval
  allocate(this%h1D_LEAF_PC_ptc(beg_ptc:end_ptc))         ;this%h1D_LEAF_PC_ptc(:)=spval
  allocate(this%h1D_CAN_RN_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_RN_ptc(:)=spval
  allocate(this%h1D_CAN_LE_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_LE_ptc(:)=spval
  allocate(this%h1D_CAN_H_ptc(beg_ptc:end_ptc))           ;this%h1D_CAN_H_ptc(:)=spval
  allocate(this%h1D_CAN_G_ptc(beg_ptc:end_ptc))           ;this%h1D_CAN_G_ptc(:)=spval
  allocate(this%h1D_CAN_TEMPC_ptc(beg_ptc:end_ptc))        ;this%h1D_CAN_TEMPC_ptc(:)=spval
  allocate(this%h1D_CAN_TEMPFN_ptc(beg_ptc:end_ptc))      ;this%h1D_CAN_TEMPFN_ptc(:)=spval
  allocate(this%h1D_CAN_CO2_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_CAN_CO2_FLX_ptc(:)=spval
  allocate(this%h1D_CAN_GPP_ptc(beg_ptc:end_ptc))         ;this%h1D_CAN_GPP_ptc(:)=spval
  allocate(this%h1D_CAN_RA_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_RA_ptc(:)=spval
  allocate(this%h1D_cTNC_ptc(beg_ptc:end_ptc))            ;this%h1D_cTNC_ptc(:)=spval
  allocate(this%h1D_cTNN_ptc(beg_ptc:end_ptc))            ;this%h1D_cTNN_ptc(:)=spval
  allocate(this%h1D_cTNP_ptc(beg_ptc:end_ptc))            ;this%h1D_cTNP_ptc(:)=spval
  allocate(this%h1D_STOML_RSC_CO2_ptc(beg_ptc:end_ptc))   ;this%h1D_STOML_RSC_CO2_ptc(:)=spval
  allocate(this%h1D_BLYR_RSC_CO2_ptc(beg_ptc:end_ptc))    ;this%h1D_BLYR_RSC_CO2_ptc(:)=spval
  allocate(this%h1D_CAN_CO2_ptc(beg_ptc:end_ptc))         ;this%h1D_CAN_CO2_ptc(:)=spval
  allocate(this%h1D_LAI_ptc(beg_ptc:end_ptc))             ;this%h1D_LAI_ptc(:)=spval
  allocate(this%h1D_PSI_CAN_ptc(beg_ptc:end_ptc))         ;this%h1D_PSI_CAN_ptc(:)=spval
  allocate(this%h1D_TURG_CAN_ptc(beg_ptc:end_ptc))        ;this%h1D_TURG_CAN_ptc(:)=spval
  allocate(this%h1D_STOM_RSC_H2O_ptc(beg_ptc:end_ptc))    ;this%h1D_STOM_RSC_H2O_ptc(:)=spval
  allocate(this%h1D_BLYR_RSC_H2O_ptc(beg_ptc:end_ptc))    ;this%h1D_BLYR_RSC_H2O_ptc(:)=spval
  allocate(this%h1D_TRANSPN_ptc(beg_ptc:end_ptc))         ;this%h1D_TRANSPN_ptc(:)=spval
  allocate(this%h1D_NH4_UPTK_FLX_ptc(beg_ptc:end_ptc))    ;this%h1D_NH4_UPTK_FLX_ptc(:)=spval
  allocate(this%h1D_NO3_UPTK_FLX_ptc(beg_ptc:end_ptc))    ;this%h1D_NO3_UPTK_FLX_ptc(:)=spval
  allocate(this%h1D_N2_FIXN_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_N2_FIXN_FLX_ptc(:)=spval
  allocate(this%h1D_cNH3_FLX_ptc(beg_ptc:end_ptc))        ;this%h1D_cNH3_FLX_ptc(:)=spval
  allocate(this%h1D_PO4_UPTK_FLX_ptc(beg_ptc:end_ptc))    ;this%h1D_PO4_UPTK_FLX_ptc(:)=spval
  allocate(this%h1D_SHOOT_C_ptc(beg_ptc:end_ptc))         ;this%h1D_SHOOT_C_ptc(:)=spval
  allocate(this%h1D_SHOOTST_C_ptc(beg_ptc:end_ptc))       ;this%h1D_SHOOTST_C_ptc(:)=spval
  allocate(this%h1D_SHOOTST_N_ptc(beg_ptc:end_ptc))       ;this%h1D_SHOOTST_N_ptc(:)=spval
  allocate(this%h1D_SHOOTST_P_ptc(beg_ptc:end_ptc))       ;this%h1D_SHOOTST_P_ptc(:)=spval
  allocate(this%h1D_Plant_C_ptc(beg_ptc:end_ptc))         ;this%h1D_Plant_C_ptc(:)=spval
  allocate(this%h1D_frcPARabs_ptc(beg_ptc:end_ptc))       ;this%h1D_frcPARabs_ptc(:)=spval
  allocate(this%h1D_PAR_CAN_ptc(beg_ptc:end_ptc))         ;this%h1D_PAR_CAN_ptc(:)=spval
  allocate(this%h1D_LEAF_C_ptc(beg_ptc:end_ptc))          ;this%h1D_LEAF_C_ptc(:)=spval
  allocate(this%h1D_Petole_C_ptc(beg_ptc:end_ptc))        ;this%h1D_Petole_C_ptc(:)=spval
  allocate(this%h1D_STALK_C_ptc(beg_ptc:end_ptc))         ;this%h1D_STALK_C_ptc(:)=spval
  allocate(this%h1D_RESERVE_C_ptc(beg_ptc:end_ptc))       ;this%h1D_RESERVE_C_ptc(:)=spval
  allocate(this%h1D_HUSK_C_ptc(beg_ptc:end_ptc))          ;this%h1D_HUSK_C_ptc(:)=spval
  allocate(this%h1D_GRAIN_C_ptc(beg_ptc:end_ptc))         ;this%h1D_GRAIN_C_ptc(:)=spval
  allocate(this%h1D_ROOT_C_ptc(beg_ptc:end_ptc))          ;this%h1D_ROOT_C_ptc(:)=spval
  allocate(this%h1D_ROOTST_C_ptc(beg_ptc:end_ptc))        ;this%h1D_ROOTST_C_ptc(:)=spval  
  allocate(this%h1D_ROOTST_N_ptc(beg_ptc:end_ptc))        ;this%h1D_ROOTST_N_ptc(:)=spval  
  allocate(this%h1D_ROOTST_P_ptc(beg_ptc:end_ptc))        ;this%h1D_ROOTST_P_ptc(:)=spval      
  allocate(this%h1D_NODULE_C_ptc(beg_ptc:end_ptc))        ;this%h1D_NODULE_C_ptc(:)=spval
  allocate(this%h1D_STORED_C_ptc(beg_ptc:end_ptc))        ;this%h1D_STORED_C_ptc(:)=spval
  allocate(this%h1D_GRAIN_NO_ptc(beg_ptc:end_ptc))        ;this%h1D_GRAIN_NO_ptc(:)=spval
  allocate(this%h1D_LAIb_ptc(beg_ptc:end_ptc))            ;this%h1D_LAIb_ptc(:)=spval
  allocate(this%h1D_EXUD_CumYr_C_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_EXUD_CumYr_C_FLX_ptc(:)=spval
  allocate(this%h1D_LITRf_C_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_LITRf_C_FLX_ptc=spval
  allocate(this%h1D_LITRf_P_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_LITRf_P_FLX_ptc=spval
  allocate(this%h1D_SURF_LITRf_C_FLX_ptc(beg_ptc:end_ptc));this%h1D_SURF_LITRf_C_FLX_ptc(:)=spval
  allocate(this%h1D_AUTO_RESP_FLX_ptc(beg_ptc:end_ptc))   ;this%h1D_AUTO_RESP_FLX_ptc(:)=spval
  allocate(this%h1D_HVST_C_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_HVST_C_FLX_ptc(:)=spval
  allocate(this%h1D_STANDING_DEAD_C_ptc(beg_ptc:end_ptc)) ;this%h1D_STANDING_DEAD_C_ptc(:)=spval
  allocate(this%h1D_FIREp_CO2_FLX_ptc(beg_ptc:end_ptc))   ;this%h1D_FIREp_CO2_FLX_ptc(:)=spval
  allocate(this%h1D_FIREp_CH4_FLX_ptc(beg_ptc:end_ptc))   ;this%h1D_FIREp_CH4_FLX_ptc(:)=spval
  allocate(this%h1D_NPP_ptc(beg_ptc:end_ptc))             ;this%h1D_NPP_ptc(:)=spval
  allocate(this%h1D_CAN_HT_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_HT_ptc(:)=spval
  allocate(this%h1D_POPN_ptc(beg_ptc:end_ptc))            ;this%h1D_POPN_ptc(:)=spval
  allocate(this%h1D_tTRANSPN_ptc(beg_ptc:end_ptc))        ;this%h1D_tTRANSPN_ptc(:)=spval
  allocate(this%h1D_WTR_STRESS_ptc(beg_ptc:end_ptc))      ;this%h1D_WTR_STRESS_ptc(:)=spval
  allocate(this%h1D_OXY_STRESS_ptc(beg_ptc:end_ptc))      ;this%h1D_OXY_STRESS_ptc(:)=spval
  allocate(this%h1D_SHOOT_N_ptc(beg_ptc:end_ptc))         ;this%h1D_SHOOT_N_ptc(:)=spval
  allocate(this%h1D_Plant_N_ptc(beg_ptc:end_ptc))         ;this%h1D_Plant_N_ptc(:)=spval
  allocate(this%h1D_LEAF_N_ptc(beg_ptc:end_ptc))          ;this%h1D_LEAF_N_ptc(:)=spval
  allocate(this%h1D_Petole_N_ptc(beg_ptc:end_ptc))       ;this%h1D_Petole_N_ptc(:)=spval
  allocate(this%h1D_STALK_N_ptc(beg_ptc:end_ptc))         ;this%h1D_STALK_N_ptc(:)=spval
  allocate(this%h1D_RESERVE_N_ptc(beg_ptc:end_ptc))       ;this%h1D_RESERVE_N_ptc(:)=spval
  allocate(this%h1D_HUSK_N_ptc(beg_ptc:end_ptc))          ;this%h1D_HUSK_N_ptc(:)=spval
  allocate(this%h1D_GRAIN_N_ptc(beg_ptc:end_ptc))         ;this%h1D_GRAIN_N_ptc(:)=spval
  allocate(this%h1D_ROOT_N_ptc(beg_ptc:end_ptc))          ;this%h1D_ROOT_N_ptc(:)=spval
  allocate(this%h1D_NODULE_N_ptc(beg_ptc:end_ptc))        ;this%h1D_NODULE_N_ptc(:)=spval
  allocate(this%h1D_STORED_N_ptc(beg_ptc:end_ptc))        ;this%h1D_STORED_N_ptc(:)=spval
  allocate(this%h1D_EXUD_N_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_EXUD_N_FLX_ptc(:)=spval
  allocate(this%h1D_Uptk_N_Flx_ptc(beg_ptc:end_ptc))      ;this%h1D_Uptk_N_Flx_ptc(:)=spval
  allocate(this%h1D_Uptk_P_Flx_ptc(beg_ptc:end_ptc))      ;this%h1D_Uptk_P_Flx_ptc(:)=spval
  allocate(this%h1D_LITRf_N_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_LITRf_N_FLX_ptc(:)=spval
  allocate(this%h1D_TL_N_FIXED_FLX_ptc(beg_ptc:end_ptc))  ;this%h1D_TL_N_FIXED_FLX_ptc(:)=spval
  allocate(this%h1D_HVST_N_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_HVST_N_FLX_ptc(:)=spval
  allocate(this%h1D_NH3can_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_NH3can_FLX_ptc(:)=spval
  allocate(this%h1D_PLANT_BALANCE_C_ptc(beg_ptc:end_ptc)) ;this%h1D_PLANT_BALANCE_C_ptc(:)=spval
  allocate(this%h1D_PLANT_BALANCE_N_ptc(beg_ptc:end_ptc)) ;this%h1D_PLANT_BALANCE_N_ptc(:)=spval
  allocate(this%h1D_PLANT_BALANCE_P_ptc(beg_ptc:end_ptc)) ;this%h1D_PLANT_BALANCE_P_ptc(:)=spval
  allocate(this%h1D_STANDING_DEAD_N_ptc(beg_ptc:end_ptc)) ;this%h1D_STANDING_DEAD_N_ptc(:)=spval
  allocate(this%h1D_FIREp_N_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_FIREp_N_FLX_ptc(:)=spval
  allocate(this%h1D_SURF_LITRf_N_FLX_ptc(beg_ptc:end_ptc));this%h1D_SURF_LITRf_N_FLX_ptc(:)=spval
  allocate(this%h1D_SHOOT_P_ptc(beg_ptc:end_ptc))         ;this%h1D_SHOOT_P_ptc(:)=spval
  allocate(this%h1D_Plant_P_ptc(beg_ptc:end_ptc))         ;this%h1D_Plant_P_ptc(:)=spval
  allocate(this%h1D_stomatal_stress_ptc(beg_ptc:end_ptc)) ;this%h1D_stomatal_stress_ptc(:)=spval
  allocate(this%h1D_CANDew_ptc(beg_ptc:end_ptc)); this%h1D_CANDew_ptc(:)=spval
  allocate(this%h1D_LEAF_P_ptc(beg_ptc:end_ptc))          ;this%h1D_LEAF_P_ptc(:)=spval
  allocate(this%h1D_Petole_P_ptc(beg_ptc:end_ptc))        ;this%h1D_Petole_P_ptc(:)=spval
  allocate(this%h1D_STALK_P_ptc(beg_ptc:end_ptc))         ;this%h1D_STALK_P_ptc(:)=spval
  allocate(this%h1D_RESERVE_P_ptc(beg_ptc:end_ptc))       ;this%h1D_RESERVE_P_ptc(:)=spval
  allocate(this%h1D_HUSK_P_ptc(beg_ptc:end_ptc))          ;this%h1D_HUSK_P_ptc(:)=spval
  allocate(this%h1D_GRAIN_P_ptc(beg_ptc:end_ptc))         ;this%h1D_GRAIN_P_ptc(:)=spval
  allocate(this%h1D_ROOT_P_ptc(beg_ptc:end_ptc))          ;this%h1D_ROOT_P_ptc(:)=spval
  allocate(this%h1D_NODULE_P_ptc(beg_ptc:end_ptc))        ;this%h1D_NODULE_P_ptc(:)=spval
  allocate(this%h1D_STORED_P_ptc(beg_ptc:end_ptc))        ;this%h1D_STORED_P_ptc(:)=spval
  allocate(this%h1D_EXUD_P_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_EXUD_P_FLX_ptc(:)=spval
  allocate(this%h1D_LITTERf_P_ptc(beg_ptc:end_ptc))       ;this%h1D_LITTERf_P_ptc(:)=spval
  allocate(this%h1D_HVST_P_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_HVST_P_FLX_ptc(:)=spval
  allocate(this%h1D_STANDING_DEAD_P_ptc(beg_ptc:end_ptc)) ;this%h1D_STANDING_DEAD_P_ptc(:)=spval
  allocate(this%h1D_FIREp_P_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_FIREp_P_FLX_ptc(:)=spval
  allocate(this%h1D_SURF_LITRf_P_FLX_ptc(beg_ptc:end_ptc));this%h1D_SURF_LITRf_P_FLX_ptc(:)=spval
  allocate(this%h1D_BRANCH_NO_ptc(beg_ptc:end_ptc))       ;this%h1D_BRANCH_NO_ptc(:)=spval
  allocate(this%h1D_Growth_Stage_ptc(beg_ptc:end_ptc))    ;this%h1D_Growth_Stage_ptc(:)=spval
  allocate(this%h1D_LEAF_NC_ptc(beg_ptc:end_ptc));        ;this%h1D_LEAF_NC_ptc(:)=spval
  allocate(this%h1D_SHOOT_NONSTC_ptc(beg_ptc:end_ptc))     ;this%h1D_SHOOT_NONSTC_ptc(:)=spval
  allocate(this%h1D_SHOOT_NONSTN_ptc(beg_ptc:end_ptc))     ;this%h1D_SHOOT_NONSTN_ptc(:)=spval
  allocate(this%h1D_SHOOT_NONSTP_ptc(beg_ptc:end_ptc))     ;this%h1D_SHOOT_NONSTP_ptc(:)=spval
  allocate(this%h1D_CFIX_lmtf_pft(beg_ptc:end_ptc))        ;this%h1D_CFIX_lmtf_pft(:)=spval
  allocate(this%h1D_MainBranchNO_pft(beg_ptc:end_ptc))     ;this%h1D_MainBranchNO_pft(:)=ispval
  allocate(this%h1D_ROOT_NONSTC_ptc(beg_ptc:end_ptc))     ;this%h1D_ROOT_NONSTC_ptc(:)=spval
  allocate(this%h1D_ROOT_NONSTN_ptc(beg_ptc:end_ptc))     ;this%h1D_ROOT_NONSTN_ptc(:)=spval
  allocate(this%h1D_ROOT_NONSTP_ptc(beg_ptc:end_ptc))     ;this%h1D_ROOT_NONSTP_ptc(:)=spval
  allocate(this%h2D_tSOC_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSOC_vr(:,:)=spval
  allocate(this%h2D_tSON_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSON_vr(:,:)=spval
  allocate(this%h2D_tSOP_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSOP_vr(:,:)=spval
  allocate(this%h2D_litrC_vr(beg_col:end_col,1:JZ))   ;this%h2D_litrC_vr(:,:)=spval
  allocate(this%h2D_litrN_vr(beg_col:end_col,1:JZ))   ;this%h2D_litrN_vr(:,:)=spval
  allocate(this%h2D_litrP_vr(beg_col:end_col,1:JZ))   ;this%h2D_litrP_vr(:,:)=spval
  allocate(this%h2D_VHeatCap_vr(beg_col:end_col,1:JZ));this%h2D_VHeatCap_vr(:,:)=spval
  allocate(this%h2D_LEAF_NODE_NO_ptc(beg_ptc:end_ptc,1:MaxNumBranches));this%h2D_LEAF_NODE_NO_ptc(:,:)=spval
  allocate(this%h2D_RUB_ACTVN_ptc(beg_ptc:end_ptc,1:MaxNumBranches));  this%h2D_RUB_ACTVN_ptc(:,:)=spval
  allocate(this%h2D_RNITRIF_vr(beg_col:end_col,1:JZ))    ;this%h2D_RNITRIF_vr(:,:)=spval
  allocate(this%h2D_CO2_vr(beg_col:end_col,1:JZ))        ;this%h2D_CO2_vr(:,:)=spval
  allocate(this%h2D_CH4_vr(beg_col:end_col,1:JZ))        ;this%h2D_CH4_vr(:,:)=spval
  allocate(this%h2D_O2_vr(beg_col:end_col,1:JZ))         ;this%h2D_O2_vr(:,:)=spval
  allocate(this%h2D_N2O_vr(beg_col:end_col,1:JZ))        ;this%h2D_N2O_vr(:,:)=spval
  allocate(this%h2D_NH3_vr(beg_col:end_col,1:JZ))        ;this%h2D_NH3_vr(:,:)=spval
  allocate(this%h2D_TEMP_vr(beg_col:end_col,1:JZ))       ;this%h2D_TEMP_vr(:,:)=spval
  allocate(this%h2D_RootMassC_vr(beg_col:end_col,1:JZ)); this%h2D_RootMassC_vr(:,:)=spval
  allocate(this%h2D_RootMassN_vr(beg_col:end_col,1:JZ)); this%h2D_RootMassN_vr(:,:)=spval
  allocate(this%h2D_RootMassP_vr(beg_col:end_col,1:JZ)); this%h2D_RootMassP_vr(:,:)=spval

  allocate(this%h2D_DOC_vr(beg_col:end_col,1:JZ)); this%h2D_DOC_vr(:,:)=spval
  allocate(this%h2D_DON_vr(beg_col:end_col,1:JZ)); this%h2D_DON_vr(:,:)=spval
  allocate(this%h2D_DOP_vr(beg_col:end_col,1:JZ)); this%h2D_DOP_vr(:,:)=spval
  allocate(this%h2D_acetate_vr(beg_col:end_col,1:JZ)); this%h2D_acetate_vr(:,:)=spval
  allocate(this%h1D_DOC_litr_col(beg_col:end_col));  this%h1D_DOC_litr_col(:)=spval
  allocate(this%h1D_DON_litr_col(beg_col:end_col)); this%h1D_DON_litr_col(:)=spval
  allocate(this%h1D_DOP_litr_col(beg_col:end_col)); this%h1D_DOP_litr_col(:)=spval
  allocate(this%h1D_acetate_litr_col(beg_col:end_col)); this%h1D_acetate_litr_col(:)=spval

  allocate(this%h2D_AeroHrBactC_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrBactC_vr(:,:)=spval
  allocate(this%h2D_AeroHrFungC_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrFungC_vr(:,:)=spval
  allocate(this%h2D_faculDenitC_vr(beg_col:end_col,1:JZ)); this%h2D_faculDenitC_vr(:,:)=spval
  allocate(this%h2D_fermentorC_vr(beg_col:end_col,1:JZ));  this%h2D_fermentorC_vr(:,:)=spval
  allocate(this%h2D_acetometgC_vr(beg_col:end_col,1:JZ));  this%h2D_acetometgC_vr(:,:)=spval
  allocate(this%h2D_aeroN2fixC_vr(beg_col:end_col,1:JZ));  this%h2D_aeroN2fixC_vr(:,:)=spval
  allocate(this%h2D_anaeN2FixC_vr(beg_col:end_col,1:JZ));  this%h2D_anaeN2FixC_vr(:,:)=spval
  allocate(this%h2D_NH3OxiBactC_vr(beg_col:end_col,1:JZ)); this%h2D_NH3OxiBactC_vr(:,:)=spval
  allocate(this%h2D_NO2OxiBactC_vr(beg_col:end_col,1:JZ)); this%h2D_NO2OxiBactC_vr(:,:)=spval
  allocate(this%h2D_CH4AeroOxiC_vr(beg_col:end_col,1:JZ)); this%h2D_CH4AeroOxiC_vr(:,:)=spval
  allocate(this%h2D_H2MethogenC_vr(beg_col:end_col,1:JZ)); this%h2D_H2MethogenC_vr(:,:)=spval

  allocate(this%h2D_RCH4ProdHydrog_vr(beg_col:end_col,1:JZ));   this%h2D_RCH4ProdHydrog_vr(:,:)=spval
  allocate(this%h2D_RCH4ProdAcetcl_vr(beg_col:end_col,1:JZ));   this%h2D_RCH4ProdAcetcl_vr(:,:)=spval
  allocate(this%h2D_RCH4Oxi_aero_vr(beg_col:end_col,1:JZ)); this%h2D_RCH4Oxi_aero_vr(:,:)=spval
  allocate(this%h2D_RFerment_vr(beg_col:end_col,1:JZ)); this%h2D_RFerment_vr(:,:)=spval
  allocate(this%h2D_nh3oxi_vr(beg_col:end_col,1:JZ));  this%h2D_nh3oxi_vr(:,:)=spval
  allocate(this%h2D_n2oprod_vr(beg_col:end_col,1:JZ));  this%h2D_n2oprod_vr(:,:)=spval

  allocate(this%h1D_RCH4ProdHydrog_litr_col(beg_col:end_col));  this%h1D_RCH4ProdHydrog_litr_col(:)=spval
  allocate(this%h1D_RCH4ProdAcetcl_litr_col(beg_col:end_col));  this%h1D_RCH4ProdAcetcl_litr_col(:)=spval
  allocate(this%h1D_RCH4Oxi_aero_litr_col(beg_col:end_col)); this%h1D_RCH4Oxi_aero_litr_col(:)=spval
  allocate(this%h1D_RFermen_litr_col(beg_col:end_col));  this%h1D_RFermen_litr_col(:)=spval
  allocate(this%h1D_nh3oxi_litr_col(beg_col:end_col)); this%h1D_nh3oxi_litr_col(:)=spval
  allocate(this%h1D_n2oprod_litr_col(beg_col:end_col));  this%h1D_n2oprod_litr_col(:)=spval

  allocate(this%h2D_AeroHrBactN_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrBactN_vr(:,:)=spval
  allocate(this%h2D_AeroHrFungN_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrFungN_vr(:,:)=spval
  allocate(this%h2D_faculDenitN_vr(beg_col:end_col,1:JZ)); this%h2D_faculDenitN_vr(:,:)=spval
  allocate(this%h2D_fermentorN_vr(beg_col:end_col,1:JZ));  this%h2D_fermentorN_vr(:,:)=spval
  allocate(this%h2D_acetometgN_vr(beg_col:end_col,1:JZ));  this%h2D_acetometgN_vr(:,:)=spval
  allocate(this%h2D_aeroN2fixN_vr(beg_col:end_col,1:JZ));  this%h2D_aeroN2fixN_vr(:,:)=spval
  allocate(this%h2D_anaeN2FixN_vr(beg_col:end_col,1:JZ));  this%h2D_anaeN2FixN_vr(:,:)=spval
  allocate(this%h2D_NH3OxiBactN_vr(beg_col:end_col,1:JZ)); this%h2D_NH3OxiBactN_vr(:,:)=spval
  allocate(this%h2D_NO2OxiBactN_vr(beg_col:end_col,1:JZ)); this%h2D_NO2OxiBactN_vr(:,:)=spval
  allocate(this%h2D_CH4AeroOxiN_vr(beg_col:end_col,1:JZ)); this%h2D_CH4AeroOxiN_vr(:,:)=spval
  allocate(this%h2D_H2MethogenN_vr(beg_col:end_col,1:JZ)); this%h2D_H2MethogenN_vr(:,:)=spval

  allocate(this%h2D_AeroHrBactP_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrBactP_vr(:,:)=spval
  allocate(this%h2D_AeroHrFungP_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrFungP_vr(:,:)=spval
  allocate(this%h2D_faculDenitP_vr(beg_col:end_col,1:JZ)); this%h2D_faculDenitP_vr(:,:)=spval
  allocate(this%h2D_fermentorP_vr(beg_col:end_col,1:JZ));  this%h2D_fermentorP_vr(:,:)=spval
  allocate(this%h2D_acetometgP_vr(beg_col:end_col,1:JZ));  this%h2D_acetometgP_vr(:,:)=spval
  allocate(this%h2D_aeroN2fixP_vr(beg_col:end_col,1:JZ));  this%h2D_aeroN2fixP_vr(:,:)=spval
  allocate(this%h2D_anaeN2FixP_vr(beg_col:end_col,1:JZ));  this%h2D_anaeN2FixP_vr(:,:)=spval
  allocate(this%h2D_NH3OxiBactP_vr(beg_col:end_col,1:JZ)); this%h2D_NH3OxiBactP_vr(:,:)=spval
  allocate(this%h2D_NO2OxiBactP_vr(beg_col:end_col,1:JZ)); this%h2D_NO2OxiBactP_vr(:,:)=spval
  allocate(this%h2D_CH4AeroOxiP_vr(beg_col:end_col,1:JZ)); this%h2D_CH4AeroOxiP_vr(:,:)=spval
  allocate(this%h2D_H2MethogenP_vr(beg_col:end_col,1:JZ)); this%h2D_H2MethogenP_vr(:,:)=spval

  allocate(this%h2D_MicroBiomeE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_MicroBiomeE_litr_col(:,:)=spval
  allocate(this%h2D_AeroHrBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_AeroHrBactE_litr_col(:,:)=spval
  allocate(this%h2D_AeroHrFungE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_AeroHrFungE_litr_col(:,:)=spval
  allocate(this%h2D_faculDenitE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_faculDenitE_litr_col(:,:)=spval
  allocate(this%h2D_fermentorE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_fermentorE_litr_col(:,:)=spval
  allocate(this%h2D_acetometgE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_acetometgE_litr_col(:,:)=spval
  allocate(this%h2D_aeroN2fixE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_aeroN2fixE_litr_col(:,:)=spval
  allocate(this%h2D_anaeN2FixE_litr_col(beg_col:end_col,1:NumPlantChemElms)); this%h2D_anaeN2FixE_litr_col(:,:)=spval
  allocate(this%h2D_NH3OxiBactE_litr_col(beg_col:end_col,1:NumPlantChemElms));this%h2D_NH3OxiBactE_litr_col(:,:)=spval
  allocate(this%h2D_NO2OxiBactE_litr_col(beg_col:end_col,1:NumPlantChemElms));this%h2D_NO2OxiBactE_litr_col(:,:)=spval
  allocate(this%h2D_CH4AeroOxiE_litr_col(beg_col:end_col,1:NumPlantChemElms));this%h2D_CH4AeroOxiE_litr_col(:,:)=spval
  allocate(this%h2D_H2MethogenE_litr_col(beg_col:end_col,1:NumPlantChemElms));this%h2D_H2MethogenE_litr_col(:,:)=spval


  allocate(this%h2D_HeatFlow_vr(beg_col:end_col,1:JZ))   ;this%h2D_HeatFlow_vr(:,:)=spval
  allocate(this%h2D_HeatUptk_vr(beg_col:end_col,1:JZ))   ;this%h2D_HeatUptk_vr(:,:)=spval
  allocate(this%h2D_VSM_vr    (beg_col:end_col,1:JZ))     ;this%h2D_VSM_vr    (:,:)=spval
  allocate(this%h2D_VSPore_vr (beg_col:end_col,1:JZ))     ;this%h2D_VSPore_vr (:,:)=spval
  allocate(this%h2D_VSICE_vr    (beg_col:end_col,1:JZ))       ;this%h2D_VSICE_vr    (:,:)=spval
  allocate(this%h2D_PSI_vr(beg_col:end_col,1:JZ))        ;this%h2D_PSI_vr(:,:)=spval
  allocate(this%h2D_PsiO_vr(beg_col:end_col,1:JZ)) ; this%h2D_PsiO_vr(:,:)=spval
  allocate(this%h2D_RootH2OUP_vr(beg_col:end_col,1:JZ))  ;this%h2D_RootH2OUP_vr(:,:)=spval
  allocate(this%h2D_cNH4t_vr(beg_col:end_col,1:JZ))      ;this%h2D_cNH4t_vr(:,:)=spval
  allocate(this%h2D_cNO3t_vr(beg_col:end_col,1:JZ))      ;this%h2D_cNO3t_vr(:,:)=spval
                                                             
  allocate(this%h2D_cPO4_vr(beg_col:end_col,1:JZ))       ;this%h2D_cPO4_vr(:,:)=spval
  allocate(this%h2D_cEXCH_P_vr(beg_col:end_col,1:JZ))    ;this%h2D_cEXCH_P_vr(:,:)=spval
  allocate(this%h2D_ECND_vr(beg_col:end_col,1:JZ))       ;this%h2D_ECND_vr(:,:)=spval
  allocate(this%h2D_PSI_RT_vr_ptc(beg_ptc:end_ptc,1:JZ))     ;this%h2D_PSI_RT_vr_ptc(:,:)=spval
  allocate(this%h2D_prtUP_NH4_vr_ptc(beg_ptc:end_ptc,1:JZ))  ;this%h2D_prtUP_NH4_vr_ptc(:,:)=spval                                                              
  allocate(this%h2D_prtUP_NO3_vr_ptc(beg_ptc:end_ptc,1:JZ))  ;this%h2D_prtUP_NO3_vr_ptc(:,:)=spval                                                              
  allocate(this%h2D_prtUP_PO4_vr_ptc(beg_ptc:end_ptc,1:JZ))  ;this%h2D_prtUP_PO4_vr_ptc(:,:)=spval                                                              
  allocate(this%h2D_DNS_RT_vr_ptc(beg_ptc:end_ptc,1:JZ))     ;this%h2D_DNS_RT_vr_ptc(:,:)=spval
  allocate(this%h2D_Root1stStrutC_ptc(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root1stStrutC_ptc=spval
  allocate(this%h2D_Root1stStrutN_ptc(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root1stStrutN_ptc=spval
  allocate(this%h2D_Root1stStrutP_ptc(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root1stStrutP_ptc=spval
  allocate(this%h2D_Root2ndStrutC_ptc(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root2ndStrutC_ptc=spval
  allocate(this%h2D_Root2ndStrutN_ptc(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root2ndStrutN_ptc=spval
  allocate(this%h2D_Root2ndStrutP_ptc(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root2ndStrutP_ptc=spval

  allocate(this%h3D_PARTS_ptc(beg_ptc:end_ptc,1:NumOfPlantMorphUnits,1:MaxNumBranches));this%h3D_PARTS_ptc(:,:,:)=spval
  !-----------------------------------------------------------------------
  ! initialize history fields 
  !--------------------------------------------------------------------
  data1d_ptr => this%h1D_tFIRE_CO2_col(beg_col:end_col) 
  call hist_addfld1d(fname='tFIRE_CO2',units='gC m-2',avgflag='A',&
    long_name='cumulative CO2 flux from fire',ptr_col=data1d_ptr)        

  data1d_ptr => this%h1D_tFIRE_CH4_col(beg_col:end_col)  
  call hist_addfld1d(fname='tFIRE_CH4',units='',avgflag='A', &
    long_name='cumulative CH4 flux from fire', ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_cNH4_LITR_col(beg_col:end_col) 
  call hist_addfld1d(fname='cNH4_LITR',units='gN NH4/g litter',avgflag='A', &
    long_name='NH4 concentration in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_cNO3_LITR_col(beg_col:end_col)        
  call hist_addfld1d(fname='cNO3_LITR',units='gN NO3/g litter',avgflag='A',&
    long_name='NO3 concentration in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_HVST_C_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_HVST_C',units='gC/m2',avgflag='A',&
    long_name='Harvested C',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_HVST_N_col(beg_col:end_col)     
  call hist_addfld1d(fname='ECO_HVST_N',units='gN/m2',avgflag='A',&
    long_name='Harvested N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_HVST_P_col(beg_col:end_col)   
  call hist_addfld1d(fname='ECO_HVST_P',units='gP/m2',avgflag='A',&
    long_name='Harvested P',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_NET_N_MIN_col(beg_col:end_col)  
  call hist_addfld1d(fname='NET_N_MIN',units='gN/m2',avgflag='A',&
    long_name='Cumulative net microbial NH4 mineralization (<0 immobilization)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tLITR_C_col(beg_col:end_col)   
  call hist_addfld1d(fname='tLITR_C',units='gC/m2',avgflag='A',&
    long_name='column integrated total litter C',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tRAD_col(beg_col:end_col)
  call hist_addfld1d(fname='RADN',units='MJ/m2/hr',avgflag='A',&
    long_name='total incoming solar radiation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tLITR_N_col(beg_col:end_col)      
  call hist_addfld1d(fname='tLITR_N',units='gN/m2',avgflag='A',&
    long_name='column integrated total litter N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tLITR_P_col(beg_col:end_col)  
  call hist_addfld1d(fname='tLITR_P',units='gP/m2',avgflag='A',&
    long_name='column integrated total litter P',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HUMUS_C_col(beg_col:end_col) 
  call hist_addfld1d(fname='HUMUS_C',units='gC/m2',avgflag='A',&
    long_name='colum integrated humus C',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HUMUS_N_col(beg_col:end_col)         
  call hist_addfld1d(fname='HUMUS_N',units='gN/m2',avgflag='A',&
    long_name='colum integrated humus N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HUMUS_P_col(beg_col:end_col)       
  call hist_addfld1d(fname='HUMUS_P',units='gP/m2',avgflag='A',&
    long_name='colum integrated humus P',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_AMENDED_C_col(beg_col:end_col)   
  call hist_addfld1d(fname='AMENDED_C',units='gC/m2',avgflag='A',&
    long_name='total fertilizer C amendment',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_AMENDED_N_col(beg_col:end_col)   
  call hist_addfld1d(fname='AMENDED_N',units='gN/m2',avgflag='A',&
    long_name='total fertilizer N amendment',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_AMENDED_P_col(beg_col:end_col)        
  call hist_addfld1d(fname='AMENDED_P',units='gP/m2',avgflag='A',&
    long_name='total fertilizer P amendment',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tLITRf_C_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='tLITRf_C',units='gC/m2/hr',avgflag='A',&
    long_name='total LitrFall C',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tLITRf_N_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='tLITRf_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total LitrFall N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tLITRf_P_FLX_col(beg_col:end_col) 
  call hist_addfld1d(fname='tLITRf_P',units='gP/m2/hr',avgflag='A',&
    long_name='total LitrFall P',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tEXCH_PO4_col(beg_col:end_col)      
  call hist_addfld1d(fname='tEXCH_PO4',units='gP/m2',avgflag='A',&
    long_name='total bioavailable (exchangeable) mineral P: H2PO4+HPO4',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUR_DOC_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUR_DOC_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='total surface DOC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUR_DON_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUR_DON_FLX',units='gN/m2',avgflag='A',&
    long_name='cumulative surface DON flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUR_DOP_FLX_col(beg_col:end_col) 
  call hist_addfld1d(fname='SUR_DOP_FLX',units='gP/m2',avgflag='A',&
    long_name='Cumulative surface DOP flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUB_DOC_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUB_DOC_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='total subsurface DOC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUB_DON_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUB_DON_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total subsurface DON flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUB_DOP_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUB_DOP_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total subsurface DOP flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUR_DIC_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUR_DIC_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='total surface DIC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUR_DIN_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUR_DIN_FLX',units='gN/m2',avgflag='A',&
    long_name='cumulative total surface DIN flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUR_DIP_FLX_col(beg_col:end_col)     
  call hist_addfld1d(fname='SUR_DIP_FLX',units='gP/m2',avgflag='A',&
    long_name='Cumulative surface DIP flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUB_DIC_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUB_DIC_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='total subsurface DIC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUB_DIN_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUB_DIN_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='landscape total subsurface DIN flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUB_DIP_FLX_col(beg_col:end_col)     
  call hist_addfld1d(fname='SUB_DIP_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total subsurface DIP flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HeatFlx2Grnd_col(beg_col:end_col)
  call hist_addfld1d(fname='HeatFlx2Grnd_col',units='MJ/m2/hr',avgflag='A',&
    long_name='Heat flux into the ground',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Qinfl2soi_col(beg_col:end_col)
  call hist_addfld1d(fname='Qinfl2soi_col',units='mm H2O/hr',avgflag='A',&
    long_name='Water flux into the ground',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Qdrain_col(beg_col:end_col)
  call hist_addfld1d(fname='Qdrain_col',units='mm H2O/hr',avgflag='A',&
    long_name='Drainage water flux out of the soil column',ptr_col=data1d_ptr)      


  IF(salt_model)THEN
    data1d_ptr => this%h1D_tSALT_DISCHG_FLX_col(beg_col:end_col) 
    call hist_addfld1d(fname='tSALT_DISCHG_FLX',units='mol/m2/hr',avgflag='A',&
      long_name='total subsurface ion flux',ptr_col=data1d_ptr)      
  endif

  data1d_ptr => this%h1D_tPRECIP_P_col(beg_col:end_col)  
  call hist_addfld1d(fname='tPRECIP_P',units='gP/m2',avgflag='A',&
    long_name='column integrated total soil precipited P',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tSoilOrgC_col(beg_col:end_col)    
  call hist_addfld1d(fname='tSoilOrgC',units='gC/m2',avgflag='A', &
    long_name='total soil organic C',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tSoilOrgN_col(beg_col:end_col)    
  call hist_addfld1d(fname='tSoilOrgN',units='gN/m2',avgflag='A', &
    long_name='total soil organic N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tSoilOrgP_col(beg_col:end_col)    
  call hist_addfld1d(fname='tSoilOrgP',units='gP/m2',avgflag='A', &
    long_name='total soil organic P',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tMICRO_C_col(beg_col:end_col)    
  call hist_addfld1d(fname='tMICRO_C',units='gC/m2',avgflag='A', &
    long_name='total micriobial C',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tMICRO_N_col(beg_col:end_col)    
  call hist_addfld1d(fname='tMICRO_N',units='gN/m2',avgflag='A', &
    long_name='total micriobial N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tMICRO_P_col(beg_col:end_col)    
  call hist_addfld1d(fname='tMICRO_P',units='gP/m2',avgflag='A', &
    long_name='total micriobial P',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_PO4_FIRE_col(beg_col:end_col)  
  call hist_addfld1d(fname='PO4_FIRE',units='gP/m2',avgflag='A',&
    long_name='cumulative PO4 flux from fire',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_cPO4_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='cPO4_LITR',units='gP/g litr',avgflag='A',&
    long_name='PO4 concentration in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_cEXCH_P_LITR_col(beg_col:end_col)     
  call hist_addfld1d(fname='cEXCH_P_LITR',units='gP/g litr',avgflag='A',&
    long_name='concentration of exchangeable inorganic P in litterr',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_NET_P_MIN_col(beg_col:end_col)    
  call hist_addfld1d(fname='NET_P_MIN',units='gP/m2',avgflag='A',&
    long_name='Cumulative net microbial P mineralization (<0 immobilization)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HUM_col(beg_col:end_col)       
  call hist_addfld1d(fname='HMAX_AIR',units='kPa',avgflag='X',&
    long_name='daily maximum vapor pressure',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HUM_col(beg_col:end_col)    
  call hist_addfld1d(fname='HMIN_AIR',units='kPa',avgflag='M',&
    long_name='daily maximum vapor pressure',ptr_col=data1d_ptr)      
    
  data1d_ptr => this%h1D_PSI_SURF_col(beg_col:end_col)    
  call hist_addfld1d(fname='PSI_SURF',units='MPa',avgflag='A',&
    long_name='soil micropore matric water potential',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SURF_ELEV_col(beg_col:end_col)  
  call hist_addfld1d(fname='SURF_ELEV',units='m',avgflag='A',&
    long_name='Surface elevation, including litter layer',ptr_col=data1d_ptr)      
    
  data1d_ptr => this%h1D_tNH4X_col(beg_col:end_col)     
  call hist_addfld1d(fname='tNH4',units='gN/m2',avgflag='A', &
    long_name='column integrated NH4+NH3',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tNO3_col(beg_col:end_col)          
  call hist_addfld1d(fname='tNO3',units='gN/m2',avgflag='A',&
    long_name='Column integrated NO3+NO2 content',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_TEMP_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='TEMP_LITR',units='oC',avgflag='A',&
    long_name='Litter layer temperature',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_TEMP_SNOW_col(beg_col:end_col)    
  call hist_addfld1d(fname='TEMP_SNOW',units='oC',avgflag='A',&
    long_name='First snow layer temperature',ptr_col=data1d_ptr)      
    
  data1d_ptr => this%h1D_OMC_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='Surf_LitrC',units='gC/m2',avgflag='A',&
    long_name='total litter residual C, including microbes',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_OMN_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='Surf_LitrN',units='gN/m2',avgflag='A',&
    long_name='total litter residual N, including microbes',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_OMP_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='Surf_LitrP',units='gP/m2',avgflag='A',&
    long_name='total litter residual P, including microbes',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ATM_CO2_col(beg_col:end_col)   
  call hist_addfld1d(fname='ATM_CO2',units='umol/mol',avgflag='A',&
    long_name='Atmospheric CO2 concentration',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ATM_CH4_col(beg_col:end_col)   
  call hist_addfld1d(fname='ATM_CH4',units='umol/mol',avgflag='A',&
    long_name='Atmospheric CH4 concentration',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_NBP_col(beg_col:end_col)   
  call hist_addfld1d(fname='NBP',units='gC/m2',avgflag='A',&
    long_name='Cumulative net biosphere productivity',ptr_col=data1d_ptr)      
    
  data1d_ptr => this%h1D_ECO_LAI_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_LAI',units='m2/m2',avgflag='A',&
    long_name='ecosystem LAI',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Eco_GPP_CumYr_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_GPP',units='gC/m2',avgflag='A',&
    long_name='cumulative ecosystem GPP',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_RA_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_RA',units='gC/m2',avgflag='A',&
    long_name='cumulative ecosystem autotrophic respiration',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Eco_NPP_CumYr_col(beg_col:end_col)      
  call hist_addfld1d(fname='ECO_NPP',units='gC/m2',avgflag='A',&
    long_name='cumulative ecosystem NPP',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Eco_HR_CumYr_col(beg_col:end_col)      
  call hist_addfld1d(fname='ECO_RH',units='gC/m2',avgflag='A',&
    long_name='cumulative ecosystem heterotrophic respiration',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tDIC_col(beg_col:end_col)       
  call hist_addfld1d(fname='tDIC',units='gC/m2',avgflag='A',&
    long_name='column integrated total soil DIC: CO2+CH4',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tSTANDING_DEAD_C_col(beg_col:end_col)     
  call hist_addfld1d(fname='tSTANDING_DEAD_C',units='gC/m2',avgflag='A',&
    long_name='total standing dead C',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tSTANDING_DEAD_N_col(beg_col:end_col)     
  call hist_addfld1d(fname='tSTANDING_DEAD_N',units='gN/m2',avgflag='A',&
    long_name='total standing dead N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tSTANDING_DEAD_P_col(beg_col:end_col)     
  call hist_addfld1d(fname='tSTANDING_DEAD_P',units='gP/m2',avgflag='A',&
    long_name='total standing dead P',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tPRECN_col(beg_col:end_col)      
  call hist_addfld1d(fname='tPRECN',units='mm/m2',avgflag='A',&
    long_name='cumulative precipitation, including irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_ET_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_ET',units='mm H2O/m2',avgflag='A',&
    long_name='cumulative total evapotranspiration',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_N2O_LITR_col(beg_col:end_col)      
  call hist_addfld1d(fname='N2O_LITR',units='g/m3',avgflag='A',&
    long_name='N2O solute concentration in soil micropres',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_NH3_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='NH3_LITR',units='g/m3',avgflag='A',&
    long_name='NH3 solute concentration in soil micropres',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SOL_RADN_col(beg_col:end_col)      
  call hist_addfld1d(fname='SOL_RADN',units='MJ/m2/day',avgflag='A',&
    long_name='shortwave radiation in solar beam',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_AIR_TEMP_col(beg_col:end_col)      
  call hist_addfld1d(fname='AIR_TEMP',units='oC',avgflag='A',&
    long_name='air temperature',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_PATM_col(beg_col:end_col)      
  call hist_addfld1d(fname='PATM',units='kPa',avgflag='A',&
    long_name='atmospheric pressure',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HUM_col(beg_col:end_col)    
  call hist_addfld1d(fname='HUM',units='kPa',avgflag='A',&
    long_name='vapor pressure',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_WIND_col(beg_col:end_col)   
  call hist_addfld1d(fname='WIND',units='m/s',avgflag='A',&
    long_name='wind speed',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_PREC_col(beg_col:end_col)   
  call hist_addfld1d(fname='PREC',units='mm H2O/m2/hr',avgflag='A',&
    long_name='Total precipitation, excluding irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SOIL_RN_col(beg_col:end_col)   
  call hist_addfld1d(fname='SOIL_RN',units='W/m2',avgflag='A',&
    long_name='total net radiation at ground surface',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SOIL_LE_col(beg_col:end_col)      
  call hist_addfld1d(fname='SOIL_LE',units='W/m2',avgflag='A',&
    long_name='Latent heat flux into ground surface',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SOIL_H_col(beg_col:end_col)      
  call hist_addfld1d(fname='SOIL_H',units='W/m2',avgflag='A',&
    long_name='Sensible heat flux into ground surface',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SOIL_G_col(beg_col:end_col)  
  call hist_addfld1d(fname='SOIL_G',units='W/m2',avgflag='A',&
    long_name='total heat flux out of ground surface (>0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_RN_col(beg_col:end_col)     
  call hist_addfld1d(fname='ECO_Radnet',units='W/m2',avgflag='A',&
    long_name='ecosystem net radiation (>0 into ecosystem)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_LE_col(beg_col:end_col)        
  call hist_addfld1d(fname='ECO_LE',units='W/m2',avgflag='A',&
    long_name='ecosystem latent heat flux (>0 into surface)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Eco_HeatSen_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_HeatS',units='W/m2',avgflag='A',&
    long_name='ecosystem sensible heat flux (>0 into surface)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_Heat2G_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_Heat2G',units='W/m2',avgflag='A',&
    long_name='ecosystem storage heat flux (>0 into surface)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_O2_LITR_col(beg_col:end_col)      
  call hist_addfld1d(fname='O2_LITR',units='g/m3',avgflag='A',&
    long_name='O2 solute concentration in litter layer',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_MIN_LWP_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='MIN_LWP',units='MPa',avgflag='A',&
    long_name='minimum daily canopy water potential',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SOIL_CO2_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SOIL_CO2_FLX',units='umol C/m2/s',avgflag='A',&
    long_name='soil CO2 flux (< 0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_CO2_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_NEE_CO2',units='umol C/m2/s',avgflag='A',&
    long_name='ecosystem net CO2 exchange',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CH4_FLX_col(beg_col:end_col)     
  call hist_addfld1d(fname='CH4_FLX',units='umol C/m2/s',avgflag='A',&
    long_name='soil CH4 flux (<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CH4_EBU_flx_col(beg_col:end_col)     
  call hist_addfld1d(fname='CH4_EBU_FLX',units='umol C/m2/s',avgflag='A',&
    long_name='soil CH4 ebullition flux (<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CH4_PLTROOT_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4_PLTROOT_FLX',units='umol C/m2/s',avgflag='A',&
    long_name='soil CH4 flux through plants(<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_O2_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='O2_FLX',units='umol O2/m2/s',avgflag='A',&
    long_name='soil O2 flux (<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CO2_LITR_col(beg_col:end_col)      
  call hist_addfld1d(fname='CO2_LITR',units='gC/m3',avgflag='A',&
    long_name='CO2 solute concentration in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_EVAPN_col(beg_col:end_col)      
  call hist_addfld1d(fname='EVAPN',units='mm H2O/m2/hr',avgflag='A',&
    long_name='total ground surface evaporation(<0 int atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CANET_col(beg_col:end_col)
  call hist_addfld1d(fname='CANET',units='mm H2O/m2/hr',avgflag='A',&
    long_name='total canopy evapotranspiration(<0 int atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tSWC_col(beg_col:end_col)  
  call hist_addfld1d(fname='tSWC',units='mmH2O/m2',avgflag='A', &
    long_name='column integrated water content (include snow)',ptr_col=data1d_ptr)        

  data1d_ptr => this%h1D_tHeat_col(beg_col:end_col)
  call hist_addfld1d(fname='tSHeat',units='MJ/m2',avgflag='A', &
    long_name='column integrated heat content (include snow)',ptr_col=data1d_ptr)        

  data1d_ptr => this%h1D_SNOWPACK_col(beg_col:end_col)      
  call hist_addfld1d(fname='SNOWPACK',units='mmH2O/m2',&
    avgflag='A',long_name='total water equivalent snow',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SURF_WTR_col(beg_col:end_col)   
  call hist_addfld1d(fname='SURF_WTR',units='m3/m3',avgflag='A',&
    long_name='Volumetric water content in surface litter layer',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SURF_ICE_col(beg_col:end_col)    
  call hist_addfld1d(fname='SURF_ICE',units='m3/m3',avgflag='A',&
    long_name='Volumetric ice content in surface litter layer',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ACTV_LYR_col(beg_col:end_col)    
  call hist_addfld1d(fname='ACTV_LYR',units='m',avgflag='A',&
    long_name='active layer depth',ptr_col=data1d_ptr)        

  data1d_ptr => this%h1D_WTR_TBL_col(beg_col:end_col)      
  call hist_addfld1d(fname='WTR_TBL',units='m',avgflag='A',&
    long_name='internal water table depth',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_sN2O_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='sN2O_FLX',units='g/m2/hr',avgflag='A',&
    long_name='soil N2O flux (<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_PAR_col(beg_col:end_col)
  call hist_addfld1d(fname='PAR',units='umol m-2 s-1',avgflag='A',&
    long_name='Direct plus diffusive incoming photosynthetic photon flux density',&
    ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_sN2G_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='sN2G_FLX',units='g/m2/hr',&
    avgflag='A',long_name='soil N2 flux (<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_sNH3_FLX_col(beg_col:end_col)       
  call hist_addfld1d(fname='sNH3_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='soil NH3 flux (<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_VHeatCap_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='vHeatCap_litr',units='MJ/m3/K',avgflag='A',&
    long_name='surface litter heat capacity',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_RUNOFF_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='RUNOFF_FLX',units='mmH2O/m2/hr',avgflag='A',&
    long_name='landscape runoff from surface water',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SEDIMENT_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='SEDIMENT_FLX',units='kg/m2/hr',avgflag='A',&
    long_name='total sediment subsurface flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_QDISCHG_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='Qdischarge_flx',units='mmH2O/m2/hr',avgflag='A',&
    long_name='grid water discharge',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HeatDISCHG_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='HeatDischarge_flx',units='MJ/m2/hr',avgflag='A',&
    long_name='grid heat flux through discharge',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_LEAF_PC_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='LEAF_rPC',units='gP/gC',avgflag='I',&
    long_name='mass based leaf PC ratio',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CAN_RN_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CAN_RN',units='W/m2',avgflag='A',&
    long_name='canopy net radiation',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CAN_LE_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CAN_LE',units='W/m2',avgflag='A',&
    long_name='canopy latent heat flux (<0 to ATM)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CAN_H_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='CAN_H',units='W/m2',avgflag='A',&
    long_name='canopy sensible heat flux',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CAN_G_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='CAN_G',units='W/m2',avgflag='A',&
    long_name='canopy storage heat flux',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CAN_TEMPC_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_TEMPC',units='oC',avgflag='A',&
    long_name='canopy temperature',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CAN_TEMPFN_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='TEMP_FN',units='none',avgflag='A',&
    long_name='canopy temperature growth function/stress,[0-1]',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CAN_CO2_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_CO2_FLX',units='umol C/m2/s',avgflag='A',&
    long_name='canopy net CO2 exchange',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CAN_GPP_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_GPP',units='gC/m2/hr',avgflag='A',&
    long_name='plant gross CO2 fixation',ptr_patch=data1d_ptr)      
  
  data1d_ptr => this%h1D_CAN_RA_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_RA',units='gC/m2/hr',avgflag='A',&
    long_name='total autotrophic respiration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_cTNC_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='cTNC',units='gC/gC',avgflag='A',&
    long_name='canopy nonstructural C concentration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_cTNN_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='cTNN',units='gN/gC',avgflag='A',&
    long_name='canopy nonstructural N concentration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_cTNP_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='cTNP',units='gP/gC',avgflag='A',&
    long_name='canopy nonstructural P concentration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STOML_RSC_CO2_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='STOML_RSC_CO2',units='s/m',avgflag='A',&
    long_name='canopy stomatal resistance for CO2',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_BLYR_RSC_CO2_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='BLYR_RSC_CO2',units='s/m',avgflag='A',&
    long_name='canopy boundary layer resistance for CO2',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CAN_CO2_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CAN_CO2',units='umol/mol',avgflag='A',&
    long_name='canopy gaesous CO2 concentration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_LAI_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='LAIstk',units='m2/m2',avgflag='A',&
    long_name='whole plant leaf area, including stalk',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_PSI_CAN_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='PSI_CAN',units='MPa',avgflag='A',&
    long_name='canopy total water potential',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_TURG_CAN_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='TURG_CAN',units='MPa',avgflag='A',&
    long_name='canopy turgor water potential',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STOM_RSC_H2O_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='STOM_RSC',units='s/m',avgflag='A',&
    long_name='canopy stomatal resistance for H2O',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_BLYR_RSC_H2O_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='BLYR_RSC_H2O',units='s/m',avgflag='A',&
    long_name='canopy boundary layer resistance for H2O',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_TRANSPN_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='QvTransp',units='mmH2O/m2/h',avgflag='A',&
    long_name='canopy transpiration (<0 into atmosphere)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_NH4_UPTK_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='UPTK_NH4_FLX',units='gN/m2/hr',&
    avgflag='A',long_name='total root uptake of NH4',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_NO3_UPTK_FLX_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='UPTK_NO3_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total root uptake of NO3',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_N2_FIXN_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='N2_FIXN_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total root N2 fixation',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_cNH3_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='cNH3_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='*canopy NH3 flux',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_PO4_UPTK_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='UPTK_PO4_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total root uptake of PO4',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SHOOT_C_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='SHOOT_C',units='gC/m2',avgflag='A',&
    long_name='plant shoot C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SHOOTST_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOTST_C',units='gC/m2',avgflag='A',&
    long_name='plant shoot structural C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SHOOTST_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOTST_N',units='gN/m2',avgflag='A',&
    long_name='plant shoot structural N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SHOOTST_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOTST_P',units='gP/m2',avgflag='A',&
    long_name='plant shoot structural P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_frcPARabs_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='frcPARabs',units='none',avgflag='A',&
    long_name='fraction of PAR absorbed by plant',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_PAR_CAN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CANOPY_PAR',units='umol m-2 s-1',avgflag='A',&
    long_name='PAR absorbed by plant canopy',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_Plant_C_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='PLANT_C',units='gC/m2',avgflag='A',&
    long_name='plant C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_LEAF_C_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='LEAF_C',units='gC/m2',avgflag='A',&
    long_name='canopy leaf C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_Petole_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='PETIOLE_C',units='gC/m2',avgflag='A',&
    long_name='canopy sheath C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STALK_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='STALK_C',units='gC/m2',avgflag='A',&
    long_name='canopy stalk C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_RESERVE_C_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='RESERVE_C',units='gC/m2',avgflag='A',&
    long_name='canopy reserve C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_HUSK_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='HUSK_C',units='gC/m2',avgflag='A',&
    long_name='canopy husk C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_GRAIN_C_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='GRAIN_C',units='gC/m2',avgflag='A',&
    long_name='canopy grain C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_ROOT_C_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='ROOT_C',units='gC/m2',avgflag='A',&
    long_name='plant root C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_ROOTST_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ROOTST_C',units='gC/m2',avgflag='A',&
    long_name='plant root structural C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_ROOTST_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ROOTST_N',units='gN/m2',avgflag='A',&
    long_name='plant root structural N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_ROOTST_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ROOTST_P',units='gP/m2',avgflag='A',&
    long_name='plant root structural P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_NODULE_C_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='NODULE_C',units='gC/m2',avgflag='A',&
    long_name='root total nodule C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STORED_C_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='SSTORED_C',units='gC/m2',avgflag='A',&
    long_name='plant seasonal storage of nonstructural C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_ROOT_NONSTC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ROOT_NONSTC',units='gC/m2',avgflag='A',&
    long_name='plant root storage of nonstructural C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_ROOT_NONSTN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ROOT_NONSTN',units='gN/m2',avgflag='A',&
    long_name='plant root storage of nonstructural N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_ROOT_NONSTP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ROOT_NONSTP',units='gP/m2',avgflag='A',&
    long_name='plant root storage of nonstructural P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SHOOT_NONSTC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_NONSTC',units='gC/m2',avgflag='A',&
    long_name='plant leaf storage of nonstructural C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SHOOT_NONSTN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_NONSTN',units='gN/m2',avgflag='A',&
    long_name='plant leaf storage of nonstructural N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SHOOT_NONSTP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_NONSTP',units='gP/m2',avgflag='A',&
    long_name='plant leaf storage of nonstructural P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CFIX_lmtf_pft(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CFIX_rC2L_pft',units='gP/m2',avgflag='A',&
    long_name='plant rubisco to light limiting ratio (>1 light-limited)',ptr_patch=data1d_ptr)      


  data1d_ptr => this%h1D_GRAIN_NO_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='GRAIN_NO',units='1/m2',avgflag='A',&
    long_name='canopy grain number',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_LAIb_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='LAI_xstk',units='m2/m2',avgflag='A',&
    long_name='whole plant leaf area, exclude stalk',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_EXUD_CumYr_C_FLX_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='EXUD_CumYr_C_FLX',units='gC/m2',avgflag='A',&
    long_name='Cumulative root organic C uptake (<0 exudation into soil)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_LITRf_C_FLX_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='LITRf_C',units='gC/m2/hr',avgflag='A',&
    long_name='total plant LitrFall C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SURF_LITRf_C_FLX_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='SURF_LITRf_C_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='Cumulative plant LitrFall C to the soil surface',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_AUTO_RESP_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='AUTO_RESP',units='gC/m2/hr',avgflag='A',&
    long_name='whole plant autotrophic respiration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_HVST_C_FLX_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='HVST_C_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='plant C harvest',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_PLANT_BALANCE_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='PLANT_BALANCE_C',units='gC/m2',avgflag='A',&
    long_name='plant C balance?',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STANDING_DEAD_C_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='STANDING_DEAD_C',units='gC/m2',avgflag='A',&
    long_name='pft Standing dead C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_FIREp_CO2_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='FIREp_CO2_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='plant CO2 from fire',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_FIREp_CH4_FLX_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='FIREp_CH4_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='plant CH4 emission from fire',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_NPP_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='NPP_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant net primary productivity',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CAN_HT_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_HT',units='m',avgflag='A',&
    long_name='Canopy height',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_POPN_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='POPN',units='1/m2',avgflag='A',&
    long_name='Plant population',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_tTRANSPN_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='tTRANSPN',units='mmH2O/m2/hr',avgflag='A',&
    long_name='Cumulative canopy evapotranspiration (>0 into atmosphere)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_WTR_STRESS_ptc(beg_ptc:end_ptc)    !HoursTooLowPsiCan_pft(NZ,NY,NX)
  call hist_addfld1d(fname='WTR_STRESS',units='hr',avgflag='A',&
    long_name='canopy plant water stress indicator: number of ' &
    //'hours PSICanopy_pft(< PSILY',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_OXY_STRESS_ptc(beg_ptc:end_ptc)    !OSTR(NZ,NY,NX)
  call hist_addfld1d(fname='OXY_STRESS',units='none',avgflag='A',&
    long_name='plant O2 stress indicator [0->1 decrease]',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SHOOT_N_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='SHOOT_N',units='gN/m2',avgflag='A',&
    long_name='plant shoot N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_Plant_N_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='PLANT_N',units='gN/m2',avgflag='A',&
    long_name='plant N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_LEAF_N_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='LEAF_N',units='gN/m2',avgflag='A',&
    long_name='Canopy leaf N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_Petole_N_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='Petiole_N',units='gN/m2',avgflag='A',&
    long_name='canopy sheath N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STALK_N_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='STALK_N',units='gN/m2',avgflag='A',&
    long_name='Canopy stalk N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_RESERVE_N_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='RESERVE_N',units='gN/m2',avgflag='A',&
    long_name='Canopy reserve N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_HUSK_N_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='HUSK_N',units='gN/m2',avgflag='A',&
    long_name='Canopy husk N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_GRAIN_N_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='GRAIN_N',units='gN/m2',avgflag='A',&
    long_name='Canopy grain C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_ROOT_N_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='ROOT_N',units='gN/m2',avgflag='A',&
    long_name='Root nitrogen',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_NODULE_N_ptc(beg_ptc:end_ptc)        
  call hist_addfld1d(fname='NODULE_N',units='gN/m2',avgflag='A',&
    long_name='root total nodule N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STORED_N_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='SSTORED_N',units='gN/m2',avgflag='A',&
    long_name='plant seasonal storage of nonstructural N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_EXUD_N_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='EXUD_CumYr_N_FLX',units='gN/m2',avgflag='A',&
    long_name='Cumulative Root organic N uptake (<0 exudation into soil)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_Uptk_N_Flx_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Uptk_N_CumYr_FLX',units='gN/m2',avgflag='A',&
    long_name='Cumulative Root N uptake (<0 exudation to soil)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_Uptk_P_Flx_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Uptk_P_CumYr_FLX',units='gP/m2',avgflag='A',&
    long_name='Cumulative Root P uptake (<0 exudation to soil)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_LITRf_N_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='LITRf_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total plant LitrFall N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_TL_N_FIXED_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='TL_N_FIXED_FLX',units='gN/m2',avgflag='A',&
    long_name='cumulative plant N2 fixation',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_HVST_N_FLX_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='HVST_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='plant N harvest',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_NH3can_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='NH3can_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total canopy NH3 flux',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_PLANT_BALANCE_N_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='PLANT_BALANCE_N',units='gC/m2',avgflag='A',&
    long_name='plant N balance?',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STANDING_DEAD_N_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='STANDING_DEAD_N',units='gN/m2',avgflag='A',&
    long_name='pft standing dead N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_FIREp_N_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='FIREp_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='plant N emission from fire',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SURF_LITRf_N_FLX_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='SURF_LITRf_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total surface LitrFall N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SHOOT_P_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='SHOOT_P',units='gP/m2',avgflag='A',&
    long_name='plant shoot P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_Plant_P_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='PLANT_P',units='gP/m2',avgflag='A',&
    long_name='plant P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_stomatal_stress_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='stomatal_stress',units='none',avgflag='A',&
    long_name='stomatal stress from root turogr [0-1 stress]',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CANDew_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CANOPY_DEW',units='mm H2O/m2',avgflag='A',&
    long_name='Cumulative canopy dew deposition',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_LEAF_P_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='LEAF_P',units='gP/m2',avgflag='A',&
    long_name='Canopy leaf P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_Petole_P_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='Petiole_P',units='gP/m2',avgflag='A',&
    long_name='canopy sheath P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STALK_P_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='STALK_P',units='gP/m2',avgflag='A',&
    long_name='Plant stalk P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_RESERVE_P_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='RESERVE_P',units='gP/m2',avgflag='A',&
    long_name='Plant reserve P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_HUSK_P_ptc(beg_ptc:end_ptc)        
  call hist_addfld1d(fname='HUSK_P',units='gP/m2',avgflag='A',&
    long_name='Husk P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_GRAIN_P_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='GRAIN_P',units='gP/m2',avgflag='A',&
    long_name='canopy grain P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_ROOT_P_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='ROOT_P',units='gP/m2',avgflag='A',&
    long_name='plant root P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_NODULE_P_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='NODULE_P',units='gP/m2',avgflag='A',&
    long_name='root total nodule P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STORED_P_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='SSTORED_P',units='gP/m2',avgflag='A',&
    long_name='plant seasonal storage of nonstructural P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_EXUD_P_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='EXUD_CumYr_P_FLX',units='gP/m2',avgflag='A',&
    long_name='Cumulative root organic P uptake (<0 exudation into soil)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_LITRf_P_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='LITRf_P_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total plant LitrFall P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_HVST_P_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='HVST_P_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='Plant P harvest',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_PLANT_BALANCE_P_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='PLANT_BALANCE_P',units='gP/m2',avgflag='A',&
    long_name='plant P balance?',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STANDING_DEAD_P_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='STANDING_DEAD_P',units='gP/m2',avgflag='A',&
    long_name='pft Standing dead P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_FIREp_P_FLX_ptc(beg_ptc:end_ptc)             
  call hist_addfld1d(fname='FIREp_P_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='plant PO4 emission from fire',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_SURF_LITRf_P_FLX_ptc(beg_ptc:end_ptc)         !SurfLitrfalStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  call hist_addfld1d(fname='SURF_LITRf_P_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='plant LitrFall P to the soil surface',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_BRANCH_NO_ptc(beg_ptc:end_ptc)            !NumOfBranches_pft(NZ,NY,NX)
  call hist_addfld1d(fname='BRANCH_NO',units='none',avgflag='I',&
    long_name='Plant branch number',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_Growth_Stage_ptc(beg_ptc:end_ptc)  !plant development stage, integer, 0-10, planting, emergence, floral_init, jointing, 
                                                               !elongation, heading, anthesis, seed_fill, see_no_set, seed_mass_set, end_seed_fill
  call hist_addfld1d(fname='Growth_Stage',units='none',avgflag='I',&
    long_name='plant development stage, integer, 0-planting, 1-emergence, 2-floral_init, 3-jointing,'// &
    '4-elongation, 5-heading, 6-anthesis, 7-seed_fill, 8-see_no_set, 9-seed_mass_set, 10-end_seed_fill',&
    ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_LEAF_NC_ptc(beg_ptc:end_ptc)            
  call hist_addfld1d(fname='LEAF_rNC',units='gN/gC',avgflag='A',&
    long_name='mass based plant leaf NC ratio',ptr_patch=data1d_ptr)      

    data1d_ptr => this%h1D_MainBranchNO_pft(beg_ptc:end_ptc)            
  call hist_addfld1d(fname='MainBranchNO',units='-',avgflag='A',&
    long_name='Main branch number',ptr_patch=data1d_ptr)      

  data2d_ptr => this%h2D_litrC_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='litrC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved litter C',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_litrN_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='litrN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved litter N',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_litrP_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='litrP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved litter P',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_tSOC_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tSOC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total soil organic C (everything organic)',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_tSON_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tSON_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total soil organic N (everything organic)',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_BotDEPZ_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='BOTDEPZ_vr',units='m',type2d='levsoi',avgflag='A',&
    long_name='Bottom depth of soil layer',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_tSOP_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tSOP_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total soil organic P (everything organic)',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_VHeatCap_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='VHeatCap_vr',units='MJ/m3/K',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved Volumetric heat capacity',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_LEAF_NODE_NO_ptc(beg_ptc:end_ptc,1:MaxNumBranches)        !NumOfLeaves_brch(MainBranchNum_pft(NZ,NY,NX),NZ,NY,NX), leaf NO
  call hist_addfld2d(fname='LEAF_NODE_NO',units='none',type2d='nbranches',avgflag='I',&
    long_name='Leaf number',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_RUB_ACTVN_ptc(beg_ptc:end_ptc,1:MaxNumBranches)      !RubiscoActivity_brch(MainBranchNum_pft(NZ,NY,NX),NZ,NY,NX), branch down-regulation of CO2 fixation
  call hist_addfld2d(fname='RUB_ACTVN',units='none',type2d='nbranches',avgflag='A',&
    long_name='branch rubisco activity for CO2 fixation, 0-1',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_RNITRIF_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RNITRIF_vr',units='gN/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Nitrification rate from NH3 oxidation in each soil layer',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_CO2_vr(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_CO2,1:JZ,NY,NX)
  call hist_addfld2d(fname='CO2_vr',units='gC/m3 water',type2d='levsoi',avgflag='A',&
    long_name='solute concentration of CO2 in soil micropre water',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_CH4_vr(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_CH4,1:JZ,NY,NX)
  call hist_addfld2d(fname='CH4_vr',units='gC/m3 water',type2d='levsoi',avgflag='A',&
    long_name='solute concentration of CH4 in soil micropre water',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_O2_vr(beg_col:end_col,1:JZ)           !trc_solcl_vr(idg_O2,1:JZ,NY,NX)
  call hist_addfld2d(fname='O2_vr',units='g/m3 water',type2d='levsoi',avgflag='A',&
    long_name='solute concentration of O2 in soil micropre water',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_N2O_vr(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_N2O,1:JZ,NY,NX)
  call hist_addfld2d(fname='N2O_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='solute concentration of N2O in soil micropre water',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_NH3_vr(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_NH3,1:JZ,NY,NX)
  call hist_addfld2d(fname='NH3_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='solute concentration of NH3 in soil micropre water',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_TEMP_vr(beg_col:end_col,1:JZ)         !TCS(1:JZ,NY,NX)
  call hist_addfld2d(fname='TEMP_vr',units='oC',type2d='levsoi',avgflag='A',&
    long_name='soil temperature profile',ptr_col=data2d_ptr)      
!-----

  data2d_ptr =>  this%h2D_RootMassC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Root C density profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_RootMassN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Root N density profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_RootMassP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Root P density profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_DOC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='DOC_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='DOC profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_DON_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='DON_vr',units='gN/m2',type2d='levsoi',avgflag='A',&
    long_name='DON profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_DOP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='DOP_vr',units='gP/m2',type2d='levsoi',avgflag='A',&
    long_name='DOP profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_acetate_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='acetate_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Acetate profile',ptr_col=data2d_ptr)      

  data1d_ptr => this%h1D_DOC_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='DOC_litr',units='gC/m2',avgflag='A',&
    long_name='DOC in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_DON_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='DON_litr',units='gN/m2',avgflag='A',&
    long_name='DON in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_DOP_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='DOP_litr',units='gP/m2',avgflag='A',&
    long_name='DOP in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_acetate_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='Acetate_litr',units='gC/m2',avgflag='A',&
    long_name='Acetate in litter',ptr_col=data1d_ptr)      

!------
  data2d_ptr =>  this%h2D_AeroHrBactC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrBacterC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic bacteria C profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_AeroHrFungC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrFungiC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic fungi C profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_faculDenitC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Facult_denitrifierC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Facultative denitrifier C biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_fermentorN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermentorC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Fermentor C biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_acetometgC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Acetic_methanogenC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aceticlastic methanogen C biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_aeroN2fixC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_N2fixerC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic N2 fixer C biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_anaeN2FixC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Anaerobic_N2fixerC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Anaerobic N2 fixer C biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_NH3OxiBactC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ammonia_OxidizerBactC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Ammonia oxidize bacteria C biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_NO2OxiBactC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nitrie_OxidizerBactC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Nitrite oxidize bacteria C profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_CH4AeroOxiC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_methanotrophC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic methanotroph C biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_H2MethogenC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Hygrogen_methanogenC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Hydrogenotrophic methanogen C biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_RCH4ProdHydrog_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Prod_hydro_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved hydrogenotrophic CH4 production rate',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_RCH4ProdAcetcl_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Prod_aceto_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved acetoclastic CH4 production rate',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_RCH4Oxi_aero_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Oxi_Aero_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved aerobic CH4 oxidation rate',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_RFerment_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Fermentation_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved fermentation rate',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_nh3oxi_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='NH3Oxid_vr',units='gN/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved NH3 oxidation rate',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_n2oprod_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='N2OProd_vr',units='gN/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total N2O production rate',ptr_col=data2d_ptr)      

!------
  data2d_ptr =>  this%h2D_AeroHrBactN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrBacterN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic bacteria N profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_AeroHrFungN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrFungiN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic fungi N profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_faculDenitN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Facult_denitrifierN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Facultative denitrifier N biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_fermentorN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermentorN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Fermentor N biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_acetometgN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Acetic_methanogenN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aceticlastic methanogen N biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_aeroN2fixN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_N2fixerN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic N2 fixer N biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_anaeN2FixN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Anaerobic_N2fixerN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Anaerobic N2 fixer N biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_NH3OxiBactN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ammonia_OxidizerBactN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Ammonia oxidize bacteria N biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_NO2OxiBactN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nitrie_OxidizerBactN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Nitrite oxidize bacteria N profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_CH4AeroOxiN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_methanotrophN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic methanotroph N biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_H2MethogenN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Hygrogen_methanogenN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Hydrogenotrophic methanogen N biomass profile',ptr_col=data2d_ptr)      

!------
  data2d_ptr =>  this%h2D_AeroHrBactP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrBacterP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic P biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_AeroHrFungP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrFungiP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic fungi P profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_faculDenitP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Facult_denitrifierP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Facultative denitrifier P biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_fermentorP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermentorP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Fermentor P biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_acetometgP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Acetic_methanogenP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aceticlastic methanogen P biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_aeroN2fixP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_N2fixerP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic N2 fixer P biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_anaeN2FixP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Anaerobic_N2fixerP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Anaerobic N2 fixer P biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_NH3OxiBactP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ammonia_OxidizerBactP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Ammonia oxidize bacteria P biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_NO2OxiBactP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nitrie_OxidizerBactP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Nitrite oxidize bacteria P profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_CH4AeroOxiP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_methanotrophP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic methanotroph P biomass profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_H2MethogenP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Hygrogen_methanogenP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Hydrogenotrophic methanogen P biomass profile',ptr_col=data2d_ptr)      
!---

  data2d_ptr =>  this%h2D_MicroBiomeE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='MicroBiomE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Total micorobial elemental biomass in litter',ptr_col=data2d_ptr)      


  data2d_ptr =>  this%h2D_AeroHrBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_HetrBacterE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_AeroHrFungE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_HetrFungiE_litr',units='gP/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic fungi elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_faculDenitE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Facult_denitrifierE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Facultative denitrifier elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_fermentorE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='FermentorE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Fermentor elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_acetometgE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Acetic_methanogenE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aceticlastic methanogen elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_aeroN2fixE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_N2fixerE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic N2 fixer elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_anaeN2FixE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Anaerobic_N2fixerE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Anaerobic N2 fixer elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_NH3OxiBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Ammonia_OxidizerBactE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Ammonia oxidize bacteria elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_NO2OxiBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Nitrie_OxidizerBactE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Nitrite oxidize bacteria elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_CH4AeroOxiE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_methanotrophE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic methanotroph elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_H2MethogenE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Hygrogen_methanogenE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Hydrogenotrophic methanogen elemental biomass in litter',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_HeatFlow_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HeatFlow_vr',units='MJ m-3 hr-1',type2d='levsoi',avgflag='A',&
    long_name='soil heat flow profile',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_HeatUptk_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HeatUptk_vr',units='MJ m-3 hr-1',type2d='levsoi',avgflag='A',&
    long_name='soil heat flow by plant water uptake',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_VSPore_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='VSPore_vr',units='m3 pore/m3 soil',type2d='levsoi',avgflag='A',&
    long_name='Volumetric soil porosity',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_VSM_vr(beg_col:end_col,1:JZ)        !ThetaH2OZ_vr(1:JZ,NY,NX)
  call hist_addfld2d(fname='rWatFLP_vr',units='m3 H2O/m3 soil pore',type2d='levsoi',avgflag='A',&
    long_name='Fraction of soil porosity filled by water',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_VSICE_vr(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='rIceFLP_vr',units='m3 ice/m3 soil pore',type2d='levsoi',avgflag='A',&
    long_name='fraction of soil porosity filled by ice',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_PSI_vr(beg_col:end_col,1:JZ)         
  call hist_addfld2d(fname='PSI_vr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='soil matric pressure+osmotic pressure',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_PsiO_vr(beg_col:end_col,1:JZ)         
  call hist_addfld2d(fname='PsiO_vr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='soil osmotic pressure',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_RootH2OUP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootH2OUptake_vr',units='mmH2O/hr',type2d='levsoi',avgflag='A',&
    long_name='soil water taken up by root',ptr_col=data2d_ptr)      
  
  data2d_ptr => this%h2D_cNH4t_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='cNH4t_vr',units='gN/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='soil NH4x concentration',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_cNO3t_vr(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='cNO3t_vr',units='gN/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='Soil NO3+NO2 concentration',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_cPO4_vr(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='cPO4_vr',units='gP/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='soil dissolved PO4 concentration',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_cEXCH_P_vr(beg_col:end_col,1:JZ)     
  call hist_addfld2d(fname='cEXCH_P_vr',units='gP/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='total exchangeable soil PO4 concentration',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_TEMP_vr(beg_col:end_col,1:JZ)  
  call hist_addfld2d(fname='TMAX_SOIL_vr',units='oC',type2d='levsoi',avgflag='X',&
    long_name='Soil maximum temperature profile',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_TEMP_vr(beg_col:end_col,1:JZ)  
  call hist_addfld2d(fname='TMIN_SOIL_vr',units='oC',type2d='levsoi',avgflag='M',&
    long_name='Soil minimum temperature profile',ptr_col=data2d_ptr)      

  data1d_ptr => this%h1D_TEMP_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='TMAX_LITR',units='oC',avgflag='X',&
    long_name='Litter maximum temperature',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_TEMP_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='TMIN_LITR',units='oC',avgflag='M',&
    long_name='Litter minimum temperature',ptr_col=data1d_ptr)      

  data2d_ptr => this%h2D_ECND_vr(beg_col:end_col,1:JZ)     
  call hist_addfld2d(fname='ECND_vr',units='dS m-1',type2d='levsoi',avgflag='A',&
    long_name='electrical conductivity',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_PSI_RT_vr_ptc(beg_ptc:end_ptc,1:JZ)  
  call hist_addfld2d(fname='PSI_RT_pvr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='root total water potential of each pft',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_Root1stStrutC_ptc(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='ROOTC_1st_pvr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Primary root structural biomass C',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_Root1stStrutN_ptc(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='ROOTN_1st_pvr',units='gN/m2',type2d='levsoi',avgflag='A',&
    long_name='Primary root structural biomass N',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_Root1stStrutP_ptc(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='ROOTP_1st_pvr',units='gP/m2',type2d='levsoi',avgflag='A',&
    long_name='Primary root structural biomass P',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_Root2ndStrutC_ptc(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='ROOTC_2nd_pvr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Secondary root structural biomass C',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_Root2ndStrutN_ptc(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='ROOTN_2nd_pvr',units='gN/m2',type2d='levsoi',avgflag='A',&
    long_name='Secondary root structural biomass N',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_Root2ndStrutP_ptc(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='ROOTP_2nd_pvr',units='gP/m2',type2d='levsoi',avgflag='A',&
    long_name='Secondary root structural biomass P',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_prtUP_NH4_vr_ptc(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='prtUP_NH4_vr',units='gN/m3/hr',type2d='levsoi',avgflag='A',&
    long_name='root uptake of NH4',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_prtUP_NO3_vr_ptc(beg_ptc:end_ptc,1:JZ)      
  call hist_addfld2d(fname='prtUP_NO3_vr',units='gN/m3/hr',type2d='levsoi',&
    avgflag='A',long_name='root uptake of NO3',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_prtUP_PO4_vr_ptc(beg_ptc:end_ptc,1:JZ)     
  call hist_addfld2d(fname='prtUP_PO4_vr',units='gP/m3/hr',type2d='levsoi',avgflag='A',&
    long_name='root uptake of PO4',ptr_patch=data2d_ptr)      

  data2d_ptr => this%h2D_DNS_RT_vr_ptc(beg_ptc:end_ptc,1:JZ)       
  call hist_addfld2d(fname='DNS_RT_vr',units='m/m3',type2d='levsoi',avgflag='A',&
    long_name='root layer length density',ptr_patch=data2d_ptr)      

  do nbr=1,MaxNumBranches
    data2d_ptr => this%h3D_PARTS_ptc(beg_ptc:end_ptc,1:NumOfPlantMorphUnits,nbr)
    write(fieldname,'(I2.2)')nbr
    call hist_addfld2d(fname='C_PARTS_brch_'//trim(fieldname),units='none',&
      type2d='pmorphunits',avgflag='A',&
      long_name='C allocation to different morph unit in branch '//trim(fieldname),ptr_patch=data2d_ptr)      
  enddo
  end subroutine init_hist_data

!----------------------------------------------------------------------
  subroutine hist_update(this,I,J,bounds)
  implicit none
  class(histdata_type) :: this
  integer, intent(in) :: I,J
  type(bounds_type), intent(in) :: bounds
  integer :: ncol,nptc
  integer :: L,NZ,NY,NX,KN,NB,NR
  real(r8) :: micBE(1:NumPlantChemElms)
  real(r8) :: DOM(idom_beg:idom_end)
  real(r8),parameter :: secs1hour=3600._r8
  real(r8),parameter :: MJ2W=1.e6_r8/secs1hour
  real(r8),parameter :: m2mm=1000._r8
  real(r8),parameter :: million=1.e6_r8
  real(r8),parameter :: gC1hour2umol1sec=1.e6_r8/(12._r8*3600._r8)
  real(r8),parameter :: gO1hour2umol1sec=1.e6_r8/(32._r8*3600._r8)
  character(len=15) :: grow_stage_str(11)=(/'Planting      ','Emergence     ','Floral_init   ', &
                                            'Jointing      ','Elongation    ','Heading       ', &
                                            'Anthesis      ','Seed_fill     ','See_no_set    ', &
                                            'Seed_mass_set ','End_seed_fill '/)
  real(r8) :: DVOLL                                            

  DO NX=bounds%NHW,bounds%NHE   
    DO NY=bounds%NVN,bounds%NVS
      ncol=get_col(NY,NX)
      this%h1D_tFIRE_CO2_col(ncol)        =  CO2byFire_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tFIRE_CH4_col(ncol)        =  CH4byFire_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_cNH4_LITR_col(ncol)        =  safe_adb(trc_solml_vr(ids_NH4,0,NY,NX)+&
        natomw*trcx_solml_vr(idx_NH4,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      this%h1D_cNO3_LITR_col(ncol)        =  safe_adb(trc_solml_vr(ids_NO3,0,NY,NX)+&
        trc_solml_vr(ids_NO2,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      
!      if(this%h1D_cNO3_LITR_col(ncol)>1.e7*VLSoilMicPMass_vr(0,NY,NX))then                  
!        write(112,*)I+J/24.,NY,NX,trc_solml_vr(idg_NH3,0,NY,NX),trc_solml_vr(ids_NO3,0,NY,NX),trc_solml_vr(ids_NO2,0,NY,NX),&
!          VLSoilMicPMass_vr(0,NY,NX)*million
!      endif
!      if(this%h1D_cNO3_LITR_col(ncol)>1.e8*VLSoilMicPMass_vr(0,NY,NX))then                  
!        stop
!      endif  
      this%h1D_ECO_HVST_N_col(ncol)   = EcoHavstElmnt_CumYr_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_NET_N_MIN_col(ncol)    = -NetNH4Mineralize_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tLITR_P_col(ncol)      = tLitrOM_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_HUMUS_C_col(ncol)      = tHumOM_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_HUMUS_N_col(ncol)      = tHumOM_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_HUMUS_P_col(ncol)      = tHumOM_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_AMENDED_P_col(ncol)    = FerPFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tLITRf_C_FLX_col(ncol) = LiterfalOrgM_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tLITRf_N_FLX_col(ncol) = LiterfalOrgM_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tLITRf_P_FLX_col(ncol) = LiterfalOrgM_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tEXCH_PO4_col(ncol)        = tHxPO4_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUR_DOP_FLX_col(ncol)      = HydroSufDOPFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUB_DOP_FLX_col(ncol)      = HydroSubsDOPFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUR_DIP_FLX_col(ncol)      = HydroSufDIPFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUB_DIP_FLX_col(ncol)      = HydroSubsDIPFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)  
      this%h1D_HeatFlx2Grnd_col(ncol)     = HeatFlx2Grnd_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_Qinfl2soi_col(ncol)        = m2mm*Qinflx2Soil_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_Qdrain_col(ncol)           = m2mm*QDrain_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUR_DON_FLX_col(ncol)      = HydroSufDONFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUB_DON_FLX_col(ncol)      = HydroSubsDONFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tSALT_DISCHG_FLX_col(ncol) = HydroIonFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUR_DIN_FLX_col(ncol)      = HydroSufDINFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUB_DIN_FLX_col(ncol)      = HydroSubsDINFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUR_DOC_FLX_col(ncol)      = HydroSufDOCFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUB_DOC_FLX_col(ncol)      = HydroSubsDOCFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUR_DIC_FLX_col(ncol)      = HydroSufDICFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUB_DIC_FLX_col(ncol)      = HydroSubsDICFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SUR_DIP_FLX_col(ncol)      = HydroSufDIPFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tPRECIP_P_col(ncol)        = tXPO4_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tMICRO_P_col(ncol)         = tMicBiome_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_PO4_FIRE_col(ncol)         = PO4byFire_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_cPO4_LITR_col(ncol)        = safe_adb(trc_solml_vr(ids_H2PO4,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      this%h1D_cEXCH_P_LITR_col(ncol)     =  patomw*safe_adb(trcx_solml_vr(idx_HPO4,0,NY,NX)+&
        trcx_solml_vr(idx_H2PO4,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      this%h1D_ECO_HVST_P_col(ncol) = EcoHavstElmnt_CumYr_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_NET_P_MIN_col(ncol)  = -NetPO4Mineralize_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_PSI_SURF_col(ncol)   = PSISoilMatricP_vr(0,NY,NX)
      this%h1D_SURF_ELEV_col(ncol)  = -CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX)+DLYR(3,0,NY,NX)
      this%h1D_tLITR_N_col(ncol)    = tLitrOM_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_AMENDED_N_col(ncol)  = FertNFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tNH4X_col(ncol)      = tNH4_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tNO3_col(ncol)       = tNO3_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tRAD_col(ncol)       = TRAD(NY,NX)
      if(this%h1D_tNH4X_col(ncol)<0._r8)then
      write(*,*)'negative',this%h1D_tNH4X_col(ncol),this%h1D_tNO3_col(ncol)
      stop
      endif      
      this%h1D_tMICRO_N_col(ncol)         = tMicBiome_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_TEMP_LITR_col(ncol)        = TCS(0,NY,NX)
      if(VcumSnowWE_col(NY,NX)<=ZEROS(NY,NX))then
        this%h1D_TEMP_SNOW_col(ncol)   = spval
      else
        this%h1D_TEMP_SNOW_col(ncol)   = TCSnow_snvr(1,NY,NX)
      endif
      this%h1D_tLITR_C_col(ncol) = tLitrOM_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)

      this%h1D_AMENDED_C_col(ncol)        = AmendCFlx_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tMICRO_C_col(ncol)         = tMicBiome_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tSoilOrgC_col(ncol)        = tSoilOrgM_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tSoilOrgN_col(ncol)        = tSoilOrgM_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tSoilOrgP_col(ncol)        = tSoilOrgM_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_OMC_LITR_col(ncol)         = SoilOrgM_vr(ielmc,0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_OMN_LITR_col(ncol)         = SoilOrgM_vr(ielmn,0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_OMP_LITR_col(ncol)         = SoilOrgM_vr(ielmp,0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_ATM_CO2_col(ncol)          = CO2E_col(NY,NX)
      this%h1D_ATM_CH4_col(ncol)          = CH4E_col(NY,NX)
      this%h1D_NBP_col(ncol)              = Eco_NBP_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_ECO_HVST_C_col(ncol)       = EcoHavstElmnt_CumYr_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_ECO_LAI_col(ncol)          = CanopyLeafArea_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_Eco_GPP_CumYr_col(ncol)    = Eco_GPP_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_ECO_RA_col(ncol)           = Eco_AutoR_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_Eco_NPP_CumYr_col(ncol)    = Eco_NPP_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_Eco_HR_CumYr_col(ncol)     = Eco_HR_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tDIC_col(ncol)             = DIC_mass_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tSTANDING_DEAD_C_col(ncol) = StandingDeadStrutElms_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tSTANDING_DEAD_N_col(ncol) = StandingDeadStrutElms_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tSTANDING_DEAD_P_col(ncol) = StandingDeadStrutElms_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tPRECN_col(ncol)           = m2mm*QRain_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_ECO_ET_col(ncol)           = m2mm*QEvap_CumYr_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_N2O_LITR_col(ncol)         = trc_solcl_vr(idg_N2O,0,NY,NX)
      this%h1D_NH3_LITR_col(ncol)         = trc_solcl_vr(idg_NH3,0,NY,NX)
      this%h1D_SOL_RADN_col(ncol)         = RadSWSolarBeam_col(NY,NX)*units%get_SecondsPerDay()
      this%h1D_AIR_TEMP_col(ncol)         = TCA_col(NY,NX)
      this%h1D_HUM_col(ncol)              = VPK_col(NY,NX)
      this%h1D_PATM_col(ncol)             = PBOT_col(NY,NX)
      this%h1D_WIND_col(ncol)             = WindSpeedAtm_col(NY,NX)/secs1hour
      this%h1D_PREC_col(ncol)             = (RainFalPrec(NY,NX)+SnoFalPrec_col(NY,NX))*m2mm/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SOIL_RN_col(ncol)          = HeatByRadiation_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SOIL_LE_col(ncol)          = HeatEvapAir2Surf_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SOIL_H_col(ncol)           = HeatSensAir2Surf_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SOIL_G_col(ncol)           = -(HeatNet2Surf_col(NY,NX)-HeatSensVapAir2Surf_col(NY,NX))*MJ2W/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_ECO_RN_col(ncol)           = Eco_NetRad_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_ECO_LE_col(ncol)           = Eco_Heat_Latent_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_Eco_HeatSen_col(ncol)      = Eco_Heat_Sens_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_ECO_Heat2G_col(ncol)       = Eco_Heat_Grnd_col(NY,NX)*MJ2W/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_O2_LITR_col(ncol)          = trc_solcl_vr(idg_O2,0,NY,NX)
      this%h1D_SOIL_CO2_FLX_col(ncol)     = SurfGasFlx_col(idg_CO2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*gC1hour2umol1sec
      this%h1D_ECO_CO2_FLX_col(ncol)      = Eco_NEE_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*gC1hour2umol1sec
      this%h1D_CH4_FLX_col(ncol)          = SurfGasFlx_col(idg_CH4,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*gC1hour2umol1sec
      this%h1D_O2_FLX_col(ncol)           = SurfGasFlx_col(idg_O2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*gO1hour2umol1sec
      this%h1D_CH4_EBU_flx_col(ncol)      = trcg_ebu_flx_col(idg_CH4,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*gO1hour2umol1sec
      this%h1D_CH4_PLTROOT_flx_col(ncol)  = trcg_pltroot_flx_col(idg_CH4,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*gO1hour2umol1sec
      this%h1D_CO2_LITR_col(ncol)         = trc_solcl_vr(idg_CO2,0,NY,NX)
      this%h1D_EVAPN_col(ncol)            = VapXAir2GSurf_col(NY,NX)*m2mm/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_CANET_col(ncol)            = QvET_col(NY,NX)*m2mm/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_RUNOFF_FLX_col(ncol)       = -QRunSurf_col(NY,NX)*m2mm/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SEDIMENT_FLX_col(ncol)     = SedmErossLoss_CumYr_col(NY,NX)*m2mm/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tSWC_col(ncol)             = WatMass_col(NY,NX)*m2mm/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_tHeat_col(ncol)            = HeatStore_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_QDISCHG_FLX_col(ncol)      = QDischar_col(NY,NX)*m2mm/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_HeatDISCHG_FLX_col(ncol)   = HeatDischar_col(NY,NX)*m2mm/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_SNOWPACK_col(ncol)         = AZMAX1((VcumSnowWE_col(NY,NX))*m2mm/AREA(3,NU(NY,NX),NY,NX))
      this%h1D_SURF_WTR_col(ncol)         = ThetaH2OZ_vr(0,NY,NX)
      this%h1D_SURF_ICE_col(ncol)         = ThetaICEZ_vr(0,NY,NX)
      this%h1D_ACTV_LYR_col(ncol)         = -(ActiveLayDepth(NY,NX)-CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX))
      this%h1D_WTR_TBL_col(ncol)          = -(DepthInternalWTBL(NY,NX)-CumDepz2LayerBot_vr(NU(NY,NX)-1,NY,NX))
      this%h1D_sN2O_FLX_col(ncol)         = SurfGasFlx_col(idg_N2O,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_sN2G_FLX_col(ncol)         = SurfGasFlx_col(idg_N2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_sNH3_FLX_col(ncol)         = SurfGasFlx_col(idg_NH3,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_PAR_col(ncol)              = RadPARSolarBeam_col(NY,NX)
      this%h1D_VHeatCap_litr_col(ncol)    = VHeatCapacity_vr(0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)

      call sumMicBiomLayL(0,NY,NX,micBE)      
      this%h2D_MicroBiomeE_litr_col(ncol,1:NumPlantChemElms) =micBE/AREA(3,NU(NY,NX),NY,NX)  

      call SumMicbGroup(0,NY,NX,micpar%mid_Aerob_HeteroBacter,MicbE)
      this%h2D_AeroHrBactE_litr_col(ncol,1:NumPlantChemElms) =micBE/AREA(3,NU(NY,NX),NY,NX)     !aerobic heterotropic bacteria

      call SumMicbGroup(0,NY,NX,micpar%mid_Aerob_Fungi,MicbE)
      this%h2D_AeroHrFungE_litr_col(ncol,1:NumPlantChemElms) = micBE/AREA(3,NU(NY,NX),NY,NX)   !aerobic heterotropic fungi

      call SumMicbGroup(0,NY,NX,micpar%mid_Facult_DenitBacter,MicbE)
      this%h2D_faculDenitE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA(3,NU(NY,NX),NY,NX)  !facultative denitrifier

      call SumMicbGroup(0,NY,NX,micpar%mid_fermentor,MicbE)
      this%h2D_fermentorE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA(3,NU(NY,NX),NY,NX)  !fermentor

      call SumMicbGroup(0,NY,NX,micpar%mid_AcetoMethanogArchea,MicbE)
      this%h2D_acetometgE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA(3,NU(NY,NX),NY,NX)  !acetogenic methanogen

      call SumMicbGroup(0,NY,NX,micpar%mid_aerob_N2Fixer,MicbE)
      this%h2D_aeroN2fixE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA(3,NU(NY,NX),NY,NX)  !aerobic N2 fixer

      call SumMicbGroup(0,NY,NX,micpar%mid_Anaerob_N2Fixer,MicbE)
      this%h2D_anaeN2FixE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA(3,NU(NY,NX),NY,NX)  !anaerobic N2 fixer

      call SumMicbGroup(0,NY,NX,micpar%mid_AmmoniaOxidBacter,MicbE,isauto=.true.)     
      this%h2D_NH3OxiBactE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA(3,NU(NY,NX),NY,NX)

      call SumMicbGroup(0,NY,NX,micpar%mid_NitriteOxidBacter,MicbE,isauto=.true.)     
      this%h2D_NO2OxiBactE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA(3,NU(NY,NX),NY,NX)

      call SumMicbGroup(0,NY,NX,micpar%mid_AerobicMethanotrofBacter,MicbE,isauto=.true.)     
      this%h2D_CH4AeroOxiE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA(3,NU(NY,NX),NY,NX)

      call SumMicbGroup(0,NY,NX,micpar%mid_H2GenoMethanogArchea,MicbE,isauto=.true.)     
      this%h2D_H2MethogenE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA(3,NU(NY,NX),NY,NX)

      call sumDOML(0,NY,NX,DOM)
      this%h1D_DOC_LITR_col(ncol)=DOM(idom_doc)
      this%h1D_DON_LITR_col(ncol)=DOM(idom_don)
      this%h1D_DOP_LITR_col(ncol)=DOM(idom_dop)
      this%h1D_acetate_LITR_col(ncol)=DOM(idom_acetate)

      this%h1D_RCH4ProdHydrog_litr_col(ncol) = RCH4ProdHydrog_vr(0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_RCH4ProdAcetcl_litr_col(ncol) = RCH4ProdAcetcl_vr(0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_RCH4Oxi_aero_litr_col(ncol)=RCH4Oxi_aero_vr(0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_RFermen_litr_col(ncol)=RFerment_vr(0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_nh3oxi_litr_col(ncol)=RNH3oxi_vr(0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%h1D_n2oprod_litr_col(ncol)= (RN2ODeniProd_vr(0,NY,NX)+RN2ONitProd_vr(0,NY,NX) &
                               +RN2OChemoProd_vr(0,NY,NX)-RN2ORedux_vr(0,NY,NX))/AREA(3,NU(NY,NX),NY,NX)


      DO L=1,JZ        
        DVOLL=DLYR(3,L,NY,NX)*AREA(3,NU(NY,NX),NY,NX)

        if(DVOLL<=0._r8)cycle
        this%h2D_RootMassC_vr(ncol,L)     = RootMassElm_vr(ielmc,L,NY,NX)/DVOLL
        this%h2D_RootMassN_vr(ncol,L)     = RootMassElm_vr(ielmn,L,NY,NX)/DVOLL
        this%h2D_RootMassP_vr(ncol,L)     = RootMassElm_vr(ielmp,L,NY,NX)/DVOLL                

        call sumDOML(L,NY,NX,DOM)
        this%h2D_DOC_vr(ncol,L)           = DOM(idom_doc)/DVOLL
        this%h2D_DON_vr(ncol,L)           = DOM(idom_don)/DVOLL
        this%h2D_DOP_vr(ncol,L)           = DOM(idom_dop)/DVOLL
        this%h2D_BotDEPZ_vr(ncol,L)       = CumDepz2LayerBot_vr(L,NY,NX)
        this%h2D_acetate_vr(ncol,L)       = DOM(idom_acetate)/DVOLL
        this%h2D_litrC_vr(ncol,L)         = litrOM_vr(ielmc,L,NY,NX)/DVOLL
        this%h2D_litrN_vr(ncol,L)         = litrOM_vr(ielmn,L,NY,NX)/DVOLL
        this%h2D_litrP_vr(ncol,L)         = litrOM_vr(ielmp,L,NY,NX)/DVOLL
        this%h2D_tSOC_vr(ncol,L)          = SoilOrgM_vr(ielmc,L,NY,NX)/DVOLL
        this%h2D_tSON_vr(ncol,L)          = SoilOrgM_vr(ielmn,L,NY,NX)/DVOLL
        this%h2D_tSOP_vr(ncol,L)          = SoilOrgM_vr(ielmp,L,NY,NX)/DVOLL
        this%h2D_VHeatCap_vr(ncol,L)      = VHeatCapacity_vr(L,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h2D_CO2_vr(ncol,L)           = trc_solcl_vr(idg_CO2,L,NY,NX)
        this%h2D_CH4_vr(ncol,L)           = trc_solcl_vr(idg_CH4,L,NY,NX)
        this%h2D_O2_vr(ncol,L)            = trc_solcl_vr(idg_O2,L,NY,NX)
        this%h2D_N2O_vr(ncol,L)           = trc_solcl_vr(idg_N2O,L,NY,NX)
        this%h2D_NH3_vr(ncol,L)           = trc_solcl_vr(idg_NH3,L,NY,NX)
        this%h2D_TEMP_vr(ncol,L)          = TCS(L,NY,NX)
        this%h2D_HeatFlow_vr(ncol,L)      = THeatFlow2Soil_vr(L,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h2D_HeatUptk_vr(ncol,L)      = THeatRootUptake_vr(L,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h2D_VSPore_vr(ncol,L)        = POROS_vr(L,NY,NX)
        this%h2D_VSM_vr    (ncol,L)       = ThetaH2OZ_vr(L,NY,NX)
        this%h2D_VSICE_vr  (ncol,L)       = ThetaICEZ_vr(L,NY,NX)
        this%h2D_PSI_vr(ncol,L)           = PSISoilMatricP_vr(L,NY,NX)+PSISoilOsmotic_vr(L,NY,NX)
        this%h2D_PsiO_vr(ncol,L)          = PSISoilOsmotic_vr(L,NY,NX)
        this%h2D_RootH2OUP_vr(ncol,L)     = TPlantRootH2OUptake_vr(L,NY,NX)
        this%h2D_cNH4t_vr(ncol,L)=  safe_adb(trc_solml_vr(ids_NH4,L,NY,NX)+trc_solml_vr(ids_NH4B,L,NY,NX) &
                                               +natomw*(trcx_solml_vr(idx_NH4,L,NY,NX)+trcx_solml_vr(idx_NH4B,L,NY,NX)),&
                                               VLSoilMicPMass_vr(L,NY,NX))
  
        this%h2D_cNO3t_vr(ncol,L)= safe_adb(trc_solml_vr(ids_NO3,L,NY,NX)+trc_solml_vr(ids_NO3B,L,NY,NX) &
                                               +trc_solml_vr(ids_NO2,L,NY,NX)+trc_solml_vr(ids_NO2B,L,NY,NX),&
                                               VLSoilMicPMass_vr(L,NY,NX))

        this%h2D_cPO4_vr(ncol,L) = safe_adb(trc_solml_vr(ids_H1PO4,L,NY,NX)+trc_solml_vr(ids_H1PO4B,L,NY,NX) &
                                               +trc_solml_vr(ids_H2PO4,L,NY,NX)+trc_solml_vr(ids_H2PO4B,L,NY,NX),&
                                               VLWatMicP_vr(L,NY,NX))
        this%h2D_cEXCH_P_vr(ncol,L)= patomw*safe_adb(trcx_solml_vr(idx_HPO4,L,NY,NX)+trcx_solml_vr(idx_H2PO4,L,NY,NX) &
                                               +trcx_solml_vr(idx_HPO4B,L,NY,NX)+trcx_solml_vr(idx_H2PO4B,L,NY,NX),&
                                               VLSoilMicPMass_vr(L,NY,NX))
        this%h2D_ECND_vr(ncol,L)     = ECND_vr(L,NY,NX)
        !aerobic heterotropic bacteria
        call SumMicbGroup(L,NY,NX,micpar%mid_Aerob_HeteroBacter,MicbE)
        this%h2D_AeroHrBactC_vr(ncol,L) = MicbE(ielmc)/DVOLL   
        this%h2D_AeroHrBactN_vr(ncol,L) = MicbE(ielmn)/DVOLL   
        this%h2D_AeroHrBactP_vr(ncol,L) = MicbE(ielmp)/DVOLL   

        !facultative denitrifier
        call SumMicbGroup(L,NY,NX,micpar%mid_Facult_DenitBacter,MicbE)
        this%h2D_faculDenitC_vr(ncol,L) = micBE(ielmc)/DVOLL  
        this%h2D_faculDenitN_vr(ncol,L) = micBE(ielmn)/DVOLL 
        this%h2D_faculDenitP_vr(ncol,L) = micBE(ielmp)/DVOLL 

        !aerobic heterotropic fungi
        call SumMicbGroup(L,NY,NX,micpar%mid_Aerob_Fungi,MicbE)
        this%h2D_AeroHrFungC_vr(ncol,L) = micBE(ielmc)/DVOLL  
        this%h2D_AeroHrFungN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_AeroHrFungP_vr(ncol,L) = micBE(ielmp)/DVOLL

        !fermentor
        call SumMicbGroup(L,NY,NX,micpar%mid_fermentor,MicbE)
        this%h2D_fermentorC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_fermentorN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_fermentorP_vr(ncol,L) = micBE(ielmp)/DVOLL

        !acetogenic methanogen
        call SumMicbGroup(L,NY,NX,micpar%mid_AcetoMethanogArchea,MicbE)
        this%h2D_acetometgC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_acetometgN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_acetometgP_vr(ncol,L) = micBE(ielmp)/DVOLL

        !aerobic N2 fixer
        call SumMicbGroup(L,NY,NX,micpar%mid_aerob_N2Fixer,MicbE)
        this%h2D_aeroN2fixC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_aeroN2fixN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_aeroN2fixP_vr(ncol,L) = micBE(ielmp)/DVOLL

        !anaerobic N2 fixer
        call SumMicbGroup(L,NY,NX,micpar%mid_Anaerob_N2Fixer,MicbE)
        this%h2D_anaeN2FixC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_anaeN2FixN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_anaeN2FixP_vr(ncol,L) = micBE(ielmp)/DVOLL

        call SumMicbGroup(L,NY,NX,micpar%mid_AmmoniaOxidBacter,MicbE,isauto=.true.)
        this%h2D_NH3OxiBactC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_NH3OxiBactN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_NH3OxiBactP_vr(ncol,L) = micBE(ielmp)/DVOLL

        call SumMicbGroup(L,NY,NX,micpar%mid_NitriteOxidBacter,MicbE,isauto=.true.)
        this%h2D_NO2OxiBactC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_NO2OxiBactN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_NO2OxiBactP_vr(ncol,L) = micBE(ielmp)/DVOLL

        call SumMicbGroup(L,NY,NX,micpar%mid_AerobicMethanotrofBacter,MicbE,isauto=.true.)
        this%h2D_CH4AeroOxiC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_CH4AeroOxiN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_CH4AeroOxiP_vr(ncol,L) = micBE(ielmp)/DVOLL

        call SumMicbGroup(L,NY,NX,micpar%mid_H2GenoMethanogArchea,MicbE,isauto=.true.)
        this%h2D_H2MethogenC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_H2MethogenN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_H2MethogenP_vr(ncol,L) = micBE(ielmp)/DVOLL

        this%h2D_RCH4ProdAcetcl_vr(ncol,L) = RCH4ProdAcetcl_vr(L,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h2D_RCH4ProdHydrog_vr(ncol,L) = RCH4ProdHydrog_vr(L,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h2D_RCH4Oxi_aero_vr(ncol,L)   = RCH4Oxi_aero_vr(L,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h2D_RFerment_vr(ncol,L)       = RFerment_vr(L,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h2D_nh3oxi_vr(ncol,L)         = RNH3oxi_vr(L,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h2D_n2oprod_vr(ncol,L)        = (RN2ODeniProd_vr(L,NY,NX)+RN2ONitProd_vr(L,NY,NX) &
                               +RN2OChemoProd_vr(L,NY,NX)-RN2ORedux_vr(L,NY,NX))/AREA(3,NU(NY,NX),NY,NX)

      ENDDO

      DO NZ=1,NP0(NY,NX)
        nptc=get_pft(NZ,NY,NX)
        this%h1D_ROOT_NONSTC_ptc(nptc)  = RootMycoNonstElms_pft(ielmc,ipltroot,NZ,NY,NX)
        this%h1D_ROOT_NONSTN_ptc(nptc)  = RootMycoNonstElms_pft(ielmn,ipltroot,NZ,NY,NX)
        this%h1D_ROOT_NONSTP_ptc(nptc)  = RootMycoNonstElms_pft(ielmp,ipltroot,NZ,NY,NX)
        this%h1D_SHOOT_NONSTC_ptc(nptc) = CanopyNonstElms_pft(ielmc,NZ,NY,NX)
        this%h1D_SHOOT_NONSTN_ptc(nptc) = CanopyNonstElms_pft(ielmn,NZ,NY,NX)
        this%h1D_SHOOT_NONSTP_ptc(nptc) = CanopyNonstElms_pft(ielmp,NZ,NY,NX)
        if(CO2FixLL_pft(NZ,NY,NX)/=spval)then
          if(CO2FixLL_pft(NZ,NY,NX)>0._r8)this%h1D_CFIX_lmtf_pft(nptc)  = CO2FixCL_pft(NZ,NY,NX)/CO2FixLL_pft(NZ,NY,NX)
        endif
        this%h1D_MIN_LWP_ptc(nptc)      = PSICanPDailyMin(NZ,NY,NX)
        this%h1D_LEAF_PC_ptc(nptc)       = safe_adb(LeafStrutElms_pft(ielmp,NZ,NY,NX)+CanopyNonstElms_pft(ielmp,NZ,NY,NX), &
                                                 LeafStrutElms_pft(ielmc,NZ,NY,NX)+CanopyNonstElms_pft(ielmc,NZ,NY,NX))
        this%h1D_CAN_RN_ptc(nptc)        = MJ2W*RadNet2Canopy_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_CAN_LE_ptc(nptc)        = MJ2W*EvapTransHeat_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_CAN_H_ptc(nptc)         = MJ2W*HeatXAir2PCan_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_CAN_G_ptc(nptc)         = MJ2W*HeatStorCanopy_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_CAN_TEMPC_ptc(nptc)     = TCelciusCanopy_pft(NZ,NY,NX)
        this%h1D_CAN_TEMPFN_ptc(nptc)    = fTCanopyGroth_pft(NZ,NY,NX)
        this%h1D_CAN_CO2_FLX_ptc(nptc)   = CO2NetFix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*gC1hour2umol1sec
        this%h1D_CAN_GPP_ptc(nptc)       = GrossCO2Fix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_CAN_RA_ptc(nptc)        = CanopyRespC_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_cTNC_ptc(nptc)          = CanopyNonstElmConc_pft(ielmc,NZ,NY,NX)
        this%h1D_cTNN_ptc(nptc)          = CanopyNonstElmConc_pft(ielmn,NZ,NY,NX)
        this%h1D_cTNP_ptc(nptc)          = CanopyNonstElmConc_pft(ielmp,NZ,NY,NX)
        this%h1D_STOML_RSC_CO2_ptc(nptc) = CanPStomaResistH2O_pft(NZ,NY,NX)*1.56_r8*secs1hour
        this%h1D_BLYR_RSC_CO2_ptc(nptc)  = CanopyBndlResist_pft(NZ,NY,NX)*1.34_r8*secs1hour
        this%h1D_CAN_CO2_ptc(nptc)       = CanopyGasCO2_pft(NZ,NY,NX)
        this%h1D_LAI_ptc(nptc)           = LeafStalkArea_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_PSI_CAN_ptc(nptc)       = PSICanopy_pft(NZ,NY,NX)
        this%h1D_TURG_CAN_ptc(nptc)      = PSICanopyTurg_pft(NZ,NY,NX)
        this%h1D_STOM_RSC_H2O_ptc(nptc)  = CanPStomaResistH2O_pft(NZ,NY,NX)*secs1hour
        this%h1D_BLYR_RSC_H2O_ptc(nptc)  = CanopyBndlResist_pft(NZ,NY,NX)*secs1hour
        this%h1D_TRANSPN_ptc(nptc)       = Transpiration_pft(NZ,NY,NX)*m2mm/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_NH4_UPTK_FLX_ptc(nptc)  = RootNH4Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_NO3_UPTK_FLX_ptc(nptc)  = RootNO3Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)        
        this%h1D_N2_FIXN_FLX_ptc(nptc)   = RootN2Fix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_cNH3_FLX_ptc(nptc)      = NH3Dep2Can_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_PO4_UPTK_FLX_ptc(nptc)  = RootH2PO4Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_frcPARabs_ptc(nptc)     = FracPARads2Canopy_pft(NZ,NY,NX)
        this%h1D_PAR_CAN_ptc(nptc)       = RadPARbyCanopy_pft(NZ,NY,NX)   !umol /m2/s        
        this%h1D_SHOOTST_C_ptc(nptc)     = ShootStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_SHOOT_C_ptc(nptc)       = ShootElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_Plant_C_ptc(nptc)       = (ShootElms_pft(ielmc,NZ,NY,NX) &
          +RootElms_pft(ielmc,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_LEAF_C_ptc(nptc)        = LeafStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_Petole_C_ptc(nptc)      = PetoleStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_STALK_C_ptc(nptc)       = StalkStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_RESERVE_C_ptc(nptc)     = StalkRsrvElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_HUSK_C_ptc(nptc)        = (HuskStrutElms_pft(ielmc,NZ,NY,NX) &
          +EarStrutElms_pft(ielmc,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_GRAIN_C_ptc(nptc)       = GrainStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_ROOT_C_ptc(nptc)        = RootElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_ROOTST_C_ptc(nptc)      = RootStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_ROOTST_N_ptc(nptc)      = RootStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_ROOTST_P_ptc(nptc)      = RootStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)        
        this%h1D_NODULE_C_ptc(nptc)      = NodulStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_STORED_C_ptc(nptc)      = SeasonalNonstElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_GRAIN_NO_ptc(nptc)      = CanopySeedNum_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_LAIb_ptc(nptc)          = CanopyLeafArea_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_EXUD_CumYr_C_FLX_ptc(nptc)       = PlantExudElm_CumYr_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_LITRf_C_FLX_ptc(nptc)      = LitrfalStrutElms_CumYr_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_SURF_LITRf_C_FLX_ptc(nptc) = SurfLitrfalStrutElms_CumYr_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_AUTO_RESP_FLX_ptc(nptc)    = GrossResp_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_HVST_C_FLX_ptc(nptc)       = EcoHavstElmnt_CumYr_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_HVST_N_FLX_ptc(nptc)       = EcoHavstElmnt_CumYr_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_HVST_P_FLX_ptc(nptc)       = EcoHavstElmnt_CumYr_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_PLANT_BALANCE_C_ptc(nptc)  = ElmBalanceCum_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_PLANT_BALANCE_N_ptc(nptc)  = ElmBalanceCum_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_PLANT_BALANCE_P_ptc(nptc)  = ElmBalanceCum_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_STANDING_DEAD_C_ptc(nptc)  = StandDeadStrutElms_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_FIREp_CO2_FLX_ptc(nptc)    = CO2ByFire_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_FIREp_CH4_FLX_ptc(nptc)    = CH4ByFire_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_NPP_ptc(nptc)              = NetPrimProduct_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_CAN_HT_ptc(nptc)           = CanopyHeight_pft(NZ,NY,NX)
        this%h1D_POPN_ptc(nptc)             = PlantPopulation_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_tTRANSPN_ptc(nptc)         = -ETCanopy_CumYr_pft(NZ,NY,NX)*m2mm/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_WTR_STRESS_ptc(nptc)       = HoursTooLowPsiCan_pft(NZ,NY,NX)
                
        this%h1D_OXY_STRESS_ptc(nptc)   = PlantO2Stress_pft(NZ,NY,NX)
        this%h1D_SHOOTST_N_ptc(nptc)    = ShootStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_SHOOT_N_ptc(nptc)      = ShootElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)        
        this%h1D_Plant_N_ptc(nptc)      = (ShootElms_pft(ielmn,NZ,NY,NX)&
          +RootElms_pft(ielmn,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_LEAF_N_ptc(nptc)    = LeafStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_Petole_N_ptc(nptc)  = PetoleStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_STALK_N_ptc(nptc)   = StalkStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_RESERVE_N_ptc(nptc) = StalkRsrvElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_HUSK_N_ptc(nptc)    = (HuskStrutElms_pft(ielmn,NZ,NY,NX) &
          +EarStrutElms_pft(ielmn,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_GRAIN_N_ptc(nptc)          = GrainStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_ROOT_N_ptc(nptc)           = RootElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_NODULE_N_ptc(nptc)         = NodulStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_STORED_N_ptc(nptc)         = SeasonalNonstElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_EXUD_N_FLX_ptc(nptc)       = PlantExudElm_CumYr_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_Uptk_N_Flx_ptc(nptc)       = RootUptk_N_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_Uptk_P_Flx_ptc(nptc)       = RootUptk_P_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_LITRf_N_FLX_ptc(nptc)      = LitrfalStrutElms_CumYr_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_TL_N_FIXED_FLX_ptc(nptc)   = PlantN2Fix_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_NH3can_FLX_ptc(nptc)       = NH3Emis_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_STANDING_DEAD_N_ptc(nptc)  = StandDeadStrutElms_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_FIREp_N_FLX_ptc(nptc)      = NH3byFire_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_SURF_LITRf_N_FLX_ptc(nptc) = SurfLitrfalStrutElms_CumYr_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_SHOOTST_P_ptc(nptc)        = ShootStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_SHOOT_P_ptc(nptc)          = ShootElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_Plant_P_ptc(nptc)          = (ShootElms_pft(ielmp,NZ,NY,NX) &
          +RootElms_pft(ielmp,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)        
        this%h1D_stomatal_stress_ptc(nptc) = StomatalStress_pft(NZ,NY,NX)
        this%h1D_CANDew_ptc(nptc)          = QdewCanopy_CumYr_pft(NZ,NY,NX)*m2mm/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_LEAF_P_ptc(nptc)          = LeafStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_Petole_P_ptc(nptc)        = PetoleStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_STALK_P_ptc(nptc)         = StalkStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_RESERVE_P_ptc(nptc)       = StalkRsrvElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_HUSK_P_ptc(nptc)          = (HuskStrutElms_pft(ielmp,NZ,NY,NX) &
          +EarStrutElms_pft(ielmp,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_GRAIN_P_ptc(nptc)          = GrainStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_ROOT_P_ptc(nptc)           = RootElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_NODULE_P_ptc(nptc)         = NodulStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_STORED_P_ptc(nptc)         = SeasonalNonstElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_EXUD_P_FLX_ptc(nptc)       = PlantExudElm_CumYr_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_LITRf_P_FLX_ptc(nptc)      = LitrfalStrutElms_CumYr_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_STANDING_DEAD_P_ptc(nptc)  = StandDeadStrutElms_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_FIREp_P_FLX_ptc(nptc)      = PO4byFire_CumYr_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_SURF_LITRf_P_FLX_ptc(nptc) = SurfLitrfalStrutElms_CumYr_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%h1D_BRANCH_NO_ptc(nptc)        = NumOfBranches_pft(NZ,NY,NX)
        this%h1D_MainBranchNO_pft(nptc)     = MainBranchNum_pft(NZ,NY,NX)
        this%h1D_LEAF_NC_ptc(nptc)      = safe_adb(LeafStrutElms_pft(ielmn,NZ,NY,NX)+CanopyNonstElms_pft(ielmn,NZ,NY,NX),&
                                                 LeafStrutElms_pft(ielmc,NZ,NY,NX)+CanopyNonstElms_pft(ielmc,NZ,NY,NX))
        if(MainBranchNum_pft(NZ,NY,NX)> 0)then
          DO KN=NumGrowthStages,0,-1
            IF(KN==0)THEN
              this%h1D_Growth_Stage_ptc(nptc)   =KN
            ELSE
              if(iPlantCalendar_brch(KN,MainBranchNum_pft(NZ,NY,NX),NZ,NY,NX)>0)then
                this%h1D_Growth_Stage_ptc(nptc) =KN
                exit    
              endif  
            ENDIF  
          ENDDO

          DO NB=1,MainBranchNum_pft(NZ,NY,NX)
            this%h2D_LEAF_NODE_NO_ptc(nptc,NB) = NumOfLeaves_brch(NB,NZ,NY,NX)
            this%h2D_RUB_ACTVN_ptc(nptc,NB)  = RubiscoActivity_brch(NB,NZ,NY,NX)
            this%h3D_PARTS_ptc(nptc,1:NumOfPlantMorphUnits,NB)=PARTS_brch(1:NumOfPlantMorphUnits,NB,NZ,NY,NX)
          ENDDO
        endif
        DO L=1,JZ
          this%h2D_PSI_RT_vr_ptc(nptc,L)    = PSIRoot_pvr(ipltroot,L,NZ,NY,NX)
          this%h2D_prtUP_NH4_vr_ptc(nptc,L) = (sum(RootNutUptake_pvr(ids_NH4,:,L,NZ,NY,NX))+&
            sum(RootNutUptake_pvr(ids_NH4B,:,L,NZ,NY,NX)))/AREA(3,L,NY,NX)
          this%h2D_prtUP_NO3_vr_ptc(nptc,L)  = (sum(RootNutUptake_pvr(ids_NO3,:,L,NZ,NY,NX))+&
            sum(RootNutUptake_pvr(ids_NO3B,:,L,NZ,NY,NX)))/AREA(3,L,NY,NX)
          this%h2D_prtUP_PO4_vr_ptc(nptc,L)  = (sum(RootNutUptake_pvr(ids_H2PO4,:,L,NZ,NY,NX))+&
            sum(RootNutUptake_pvr(ids_H2PO4B,:,L,NZ,NY,NX)))/AREA(3,L,NY,NX)
          this%h2D_DNS_RT_vr_ptc(nptc,L) = RootLenDensPerPlant_pvr(ipltroot,L,NZ,NY,NX)* &
            PlantPopulation_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)

          this%h2D_Root1stStrutC_ptc(nptc,L)=0._r8
          this%h2D_Root1stStrutN_ptc(nptc,L)=0._r8
          this%h2D_Root1stStrutP_ptc(nptc,L)=0._r8
          this%h2D_Root2ndStrutC_ptc(nptc,L)=0._r8
          this%h2D_Root2ndStrutN_ptc(nptc,L)=0._r8
          this%h2D_Root2ndStrutP_ptc(nptc,L)=0._r8

          DO NR=1,NumRootAxes_pft(NZ,NY,NX)
            this%h2D_Root1stStrutC_ptc(nptc,L)= this%h2D_Root1stStrutC_ptc(nptc,L) + &
              RootMyco1stStrutElms_rpvr(ielmc,ipltroot,L,NR,NZ,NY,NX)
            this%h2D_Root1stStrutN_ptc(nptc,L)= this%h2D_Root1stStrutN_ptc(nptc,L) + &
              RootMyco1stStrutElms_rpvr(ielmn,ipltroot,L,NR,NZ,NY,NX)
            this%h2D_Root1stStrutP_ptc(nptc,L)= this%h2D_Root1stStrutP_ptc(nptc,L) + &
              RootMyco1stStrutElms_rpvr(ielmp,ipltroot,L,NR,NZ,NY,NX)
            this%h2D_Root2ndStrutC_ptc(nptc,L)=this%h2D_Root2ndStrutC_ptc(nptc,L) + &
             RootMyco2ndStrutElms_rpvr(ielmc,ipltroot,L,NR,NZ,NY,NX)            
            this%h2D_Root2ndStrutN_ptc(nptc,L)=this%h2D_Root2ndStrutN_ptc(nptc,L) + &
             RootMyco2ndStrutElms_rpvr(ielmn,ipltroot,L,NR,NZ,NY,NX)            
            this%h2D_Root2ndStrutP_ptc(nptc,L)=this%h2D_Root2ndStrutP_ptc(nptc,L) + &
             RootMyco2ndStrutElms_rpvr(ielmp,ipltroot,L,NR,NZ,NY,NX)            
          ENDDO

        ENDDO
      ENDDO 
    ENDDO 
  ENDDO    
  end subroutine hist_update
  
end module HistDataType
