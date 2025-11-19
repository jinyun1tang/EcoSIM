module HistDataType
!
! this module is an intermediate step to support ascii output
! when output is done with netcdf, no id is needed.
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use data_const_mod,   only: spval  => DAT_CONST_SPVAL, ispval => DAT_CONST_ISPVAL
  use SoilBGCNLayMod,   only: SumMicbGroup, sumDOML, sumMicBiomLayL,SumSolidOML
  use UnitMod,          only: units
  use MiniMathMod,      only: safe_adb, AZMAX1,AZERO
  use EcoSiMParDataMod, only: pltpar, micpar
  use DebugToolMod     
  use EcoSIMCtrlMod
  use MicrobialDataType
  use CanopyRadDataType
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
  use SurfLitterDataType
  use SoilBGCDataType
  use AqueChemDatatype
  use SurfSoilDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  type, public :: histdata_type
  real(r8),pointer   :: h1D_cumFIRE_CO2_col(:)    
  real(r8),pointer   :: h1D_cumFIRE_CH4_col(:)    
  real(r8),pointer   :: h1D_cNH4_LITR_col(:)    
  real(r8),pointer   :: h1D_cNO3_LITR_col(:)    
  real(r8),pointer   :: h1D_ECO_HVST_N_col(:)   
  real(r8),pointer   :: h1D_NET_N_MIN_col(:)    
  real(r8),pointer   :: h1D_tLITR_P_col(:)   
  real(r8),pointer   :: h1D_HUMUS_C_col(:)   
  real(r8),pointer   :: h1D_HUMUS_N_col(:)   
  real(r8),pointer   :: h1D_HUMUS_P_col(:)   
  real(r8),pointer   :: h1D_AMENDED_P_col(:)    
  real(r8),pointer   :: h1D_tLITRf_P_FLX_col(:) 
  real(r8),pointer   :: h1D_tEXCH_PO4_col(:)    
  real(r8),pointer   :: h1D_SUR_DOP_FLX_col(:)  
  real(r8),pointer   :: h1D_SUB_DOP_FLX_col(:)  
  real(r8),pointer   :: h1D_SUR_DIP_FLX_col(:)  
  real(r8),pointer   :: h1D_SUB_DIP_FLX_col(:)  
  real(r8),pointer   :: h1D_HeatFlx2Grnd_col(:)  
  real(r8),pointer   :: h1D_RadSW_Grnd_col(:)  
  real(r8),pointer   :: h1D_Qinfl2soi_col(:)   
  real(r8),pointer   :: h1D_QTRANSP_col(:)
  real(r8),pointer   :: h1D_Qdrain_col(:)            
  real(r8),pointer   :: h1D_tPREC_P_col(:)     
  real(r8),pointer   :: h1D_tMICRO_P_col(:)    
  real(r8),pointer   :: h1D_tSoilOrgC_col(:)
  real(r8),pointer   :: h1D_tSoilOrgN_col(:)
  real(r8),pointer   :: h1D_tSoilOrgP_col(:)
  real(r8),pointer   :: h1D_PO4_FIRE_col(:)    
  real(r8),pointer   :: h1D_cPO4_LITR_col(:)    
  real(r8),pointer   :: h1D_cEXCH_P_LITR_col(:) 
  real(r8),pointer   :: h1D_ECO_HVST_P_col(:)   
  real(r8),pointer   :: h1D_NET_P_MIN_col(:)    
  real(r8),pointer   :: h1D_tSALT_DISCHG_FLX_col(:) 
  real(r8),pointer   :: h1D_PSI_SURF_col(:)     
  real(r8),pointer   :: h1D_SURF_ELEV_col(:)    
  real(r8),pointer   :: h1D_tLITR_N_col(:)      
  real(r8),pointer   :: h2D_RootAR_vr(:,:)      
  real(r8),pointer   :: h1D_RootAR_col(:)  
  real(r8),pointer   :: h1D_RootAR_ptc(:)
  real(r8),pointer   :: h1D_RootLenPerPlant_ptc(:)
  real(r8),pointer   :: h1D_RootCO2Relez_col(:)     
  real(r8),pointer   :: h2D_BotDEPZ_vr(:,:)
  real(r8),pointer   :: h1D_tRAD_col(:)
  real(r8),pointer   :: h1D_AMENDED_N_col(:)    
  real(r8),pointer   :: h1D_tLITRf_N_FLX_col(:) 
  real(r8),pointer   :: h1D_tLITRf_C_FLX_col(:) 
  real(r8),pointer   :: h1D_tNH4X_col(:)        
  real(r8),pointer   :: h1D_tNO3_col(:)         
  real(r8),pointer   :: h1D_SUR_DON_FLX_col(:)  
  real(r8),pointer   :: h1D_SUB_DON_FLX_col(:)  
  real(r8),pointer   :: h1D_SUR_DIN_FLX_col(:)  
  real(r8),pointer   :: h1D_SUB_DIN_FLX_col(:)  
  real(r8),pointer   :: h1D_tMICRO_N_col(:)     
  real(r8),pointer   :: h1D_TEMP_LITR_col(:)    
  real(r8),pointer   :: h1D_TEMP_surf_col(:)
  real(r8),pointer   :: h1D_TEMP_SNOW_col(:)    
  real(r8),pointer   :: h1D_FracBySnow_col(:)   
  real(r8),pointer   :: h1D_FracByLitr_col(:)   
  real(r8),pointer   :: h1D_tLITR_C_col(:)      
  real(r8),pointer   :: h1D_AMENDED_C_col(:)    
  real(r8),pointer   :: h1D_tMICRO_C_col(:)     
  real(r8),pointer   :: h1D_OMC_LITR_col(:)     
  real(r8),pointer   :: h1D_OMN_LITR_col(:)     
  real(r8),pointer   :: h1D_OMP_LITR_col(:)     
  real(r8),pointer   :: h1D_SUR_DOC_FLX_col(:)  
  real(r8),pointer   :: h1D_SUB_DOC_FLX_col(:)  
  real(r8),pointer   :: h1D_SUR_DIC_FLX_col(:)  
  real(r8),pointer   :: h1D_SUB_DIC_FLX_col(:)  
  real(r8),pointer   :: h1D_ATM_CO2_col(:)      
  real(r8),pointer   :: h1D_ATM_CH4_col(:)      
  real(r8),pointer   :: h1D_NBP_col(:)          
  real(r8),pointer   :: h1D_CanSWRad_col(:)
  real(r8),pointer   :: h1D_ECO_HVST_C_col(:)  
  real(r8),pointer   :: h1D_ECO_LAI_col(:)
  real(r8),pointer   :: h1D_ECO_SAI_col(:)
  real(r8),pointer   :: h1D_Eco_GPP_CumYr_col(:)   
  real(r8),pointer   :: h1D_ECO_RA_col(:)        
  real(r8),pointer   :: h1D_Eco_NPP_CumYr_col(:) 
  real(r8),pointer   :: h1D_Eco_HR_CumYr_col(:)  
  real(r8),pointer   :: h1D_Eco_HR_CO2_col(:)   
  real(r8),pointer   :: h1D_Eco_HR_CO2_litr_col(:)   
  real(r8),pointer   :: h2D_Eco_HR_CO2_vr(:,:)   
  real(r8),pointer   :: h2D_Gchem_CO2_prod_vr(:,:)
  real(r8),pointer   :: h1D_tDIC_col(:)       
  real(r8),pointer   :: h1D_tSTANDING_DEAD_C_col(:)  
  real(r8),pointer   :: h1D_tSTANDING_DEAD_N_col(:)  
  real(r8),pointer   :: h1D_tSTANDING_DEAD_P_col(:)  
  real(r8),pointer   :: h1D_tPRECIP_col(:)        
  real(r8),pointer   :: h1D_ECO_ET_col(:)         
  real(r8),pointer   :: h1D_trcg_Ar_cumerr_col(:)
  real(r8),pointer   :: h1D_trcg_CO2_cumerr_col(:)
  real(r8),pointer   :: h1D_trcg_CH4_cumerr_col(:)
  real(r8),pointer   :: h1D_trcg_O2_cumerr_col(:)
  real(r8),pointer   :: h1D_trcg_N2_cumerr_col(:)
  real(r8),pointer   :: h1D_trcg_NH3_cumerr_col(:)
  real(r8),pointer   :: h1D_trcg_H2_cumerr_col(:)
  real(r8),pointer   :: h1D_ECO_RADSW_col(:)
  real(r8),pointer   :: h1d_CAN_NEE_col(:)
  real(r8),pointer   :: h1D_N2O_LITR_col(:)     
  real(r8),pointer   :: h1D_NH3_LITR_col(:)     
  real(r8),pointer   :: h1D_SOL_RADN_col(:)     
  real(r8),pointer   :: h1D_AIR_TEMP_col(:)     
  real(r8),pointer   :: h1D_FreeNFix_col(:)
  real(r8),pointer   :: h1D_PATM_col(:)         
  real(r8),pointer   :: h1D_HUM_col(:)          
  real(r8),pointer   :: h1D_WIND_col(:)         
  real(r8),pointer   :: h1D_Snofall_col(:)
  real(r8),pointer   :: h1D_PREC_col(:)   
  real(r8),pointer   :: h1D_LWSky_col(:)
  real(r8),pointer   :: h1D_SOIL_RN_col(:)      
  real(r8),pointer   :: h1D_SOIL_LE_col(:)      
  real(r8),pointer   :: h1D_SOIL_H_col(:)       
  real(r8),pointer   :: h1D_SOIL_G_col(:)       
  real(r8),pointer   :: h1D_ECO_RN_col(:)       
  real(r8),pointer   :: h1D_ECO_LE_col(:)       
  real(r8),pointer   :: h1D_Eco_HeatSen_col(:)  
  real(r8),pointer   :: h1D_ECO_Heat2G_col(:)   
  real(r8),pointer   :: h1D_O2_LITR_col(:)      
  real(r8),pointer   :: h1D_MIN_LWP_ptc(:)      
  real(r8),pointer   :: h1D_AR_SEMIS_FLX_col(:)   
  real(r8),pointer   :: h1D_CO2_SEMIS_FLX_col(:)  
  real(r8),pointer   :: h1D_ECO_CO2_FLX_col(:)  
  real(r8),pointer   :: h1D_CH4_SEMIS_FLX_col(:)
  real(r8),pointer   :: h1D_CH4_EBU_flx_col(:)
  real(r8),pointer   :: h1D_Ar_EBU_flx_col(:)
  real(r8),pointer   :: h1D_CO2_TPR_err_col(:)
  real(r8),pointer   :: h1D_CO2_Drain_flx_col(:)
  real(r8),pointer   :: h1D_CO2_hydloss_flx_col(:)
  real(r8),pointer   :: h1D_Ar_TPR_err_col(:)
  real(r8),pointer   :: h1D_AR_PLTROOT_flx_col(:)  
  real(r8),pointer   :: h1D_CH4_PLTROOT_flx_col(:)
  real(r8),pointer   :: h1D_CO2_PLTROOT_flx_col(:)
  real(r8),pointer   :: h1D_O2_PLTROOT_flx_col(:)
  real(r8),pointer   :: h1D_CO2_DIF_flx_col(:)
  real(r8),pointer   :: h1D_O2_DIF_flx_col(:)  
  real(r8),pointer   :: h1D_Ar_DIF_flx_col(:)
  real(r8),pointer   :: h1D_CH4_DIF_flx_col(:)
  real(r8),pointer   :: h1D_NH3_DIF_flx_col(:)
  real(r8),pointer   :: h1D_Ar_soilMass_col(:)
  real(r8),pointer   :: h1D_O2_SEMIS_FLX_col(:)  
  real(r8),pointer   :: h1D_CO2_LITR_col(:)      
  real(r8),pointer   :: h1D_EVAPN_col(:)         
  real(r8),pointer   :: h1D_CANET_col(:)         
  real(r8),pointer   :: h1D_RUNOFF_FLX_col(:)    
  real(r8),pointer   :: h1D_SEDIMENT_FLX_col(:)  
  real(r8),pointer   :: h1D_tSWC_col(:)       
  real(r8),pointer   :: h1D_tHeat_col(:) 
  real(r8),pointer   :: h1D_QDISCHG_FLX_col(:)  
  real(r8),pointer   :: h1D_HeatDISCHG_FLX_col(:)
  real(r8),pointer   :: h1D_SNOWPACK_col(:)    
  real(r8),pointer   :: h1D_SNOWDENS_col(:) 
  real(r8),pointer   :: h1D_SURF_WTR_col(:)     
  real(r8),pointer   :: h1D_ThetaW_litr_col(:)  
  real(r8),pointer   :: h1D_ThetaI_litr_col(:)    
  real(r8),pointer   :: h1D_SURF_ICE_col(:)     
  real(r8),pointer   :: h1D_ACTV_LYR_col(:)     
  real(r8),pointer   :: h1D_WTR_TBL_col(:)      
  real(r8),pointer   :: h1D_CO2_WetDep_FLX_col(:)
  real(r8),pointer   :: h1D_RootN_Fix_col(:)
  real(r8),pointer   :: h1D_AR_WetDep_FLX_col(:)
  real(r8),pointer   :: h1D_RootXO2_flx_col(:)
  real(r8),pointer   :: h1D_N2O_SEMIS_FLX_col(:)     
  real(r8),pointer   :: h1D_N2_SEMIS_FLX_col(:)      
  real(r8),pointer   :: h1D_NH3_SEMIS_FLX_col(:)     
  real(r8),pointer   :: h1D_H2_SEMIS_FLX_col(:)
  real(r8),pointer   :: h1D_frcPARabs_ptc(:)  
  real(r8),pointer   :: h1D_PAR_CAN_ptc(:)    
  real(r8),pointer   :: h1D_PAR_col(:)        
  real(r8),pointer   :: h1d_fPAR_col(:)
  real(r8),pointer   :: h1D_Plant_C_ptc(:)    
  real(r8),pointer   :: h1D_Plant_N_ptc(:)    
  real(r8),pointer   :: h1D_Plant_P_ptc(:)    
  real(r8),pointer   :: h1D_stomatal_stress_ptc(:)
  real(r8),pointer   :: h1D_CANDew_ptc(:)
  real(r8),pointer   :: h1D_VHeatCap_litr_col(:)
  real(r8),pointer   :: h1D_LEAF_PC_ptc(:)    
  real(r8),pointer   :: h3D_SOC_Cps_vr(:,:,:)  
  real(r8),pointer   :: h2D_tSOC_vr(:,:)      
  real(r8),pointer   :: h2D_POM_C_vr(:,:)
  real(r8),pointer   :: h2D_MAOM_C_vr(:,:)
  real(r8),pointer   :: h2D_microbC_vr(:,:)
  real(r8),pointer   :: h2D_microbN_vr(:,:)
  real(r8),pointer   :: h2D_microbP_vr(:,:)
  real(r8),pointer   :: h2D_AeroBact_PrimS_lim_vr(:,:)
  real(r8),pointer   :: h2D_AeroFung_PrimS_lim_vr(:,:)
  real(r8),pointer   :: h2D_tSOCL_vr(:,:)
  real(r8),pointer   :: h2D_fTRootGro_pvr(:,:) !
  real(r8),pointer   :: h2D_fRootGrowPSISense_pvr(:,:)  !
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
  real(r8),pointer   :: h1D_tDOC_soil_col(:)
  real(r8),pointer   :: h1D_tDON_soil_col(:)
  real(r8),pointer   :: h1D_tDOP_soil_col(:)
  real(r8),pointer   :: h1D_tAcetate_soil_col(:)
  real(r8),pointer   :: h1D_DOC_litr_col(:)
  real(r8),pointer   :: h1D_DON_litr_col(:)
  real(r8),pointer   :: h1D_DOP_litr_col(:)
  real(r8),pointer   :: h1D_acetate_litr_col(:)
  real(r8),pointer   :: h1D_CAN_RN_ptc(:)     
  real(r8),pointer   :: h1D_CAN_LE_ptc(:)     
  real(r8),pointer   :: h1D_CAN_H_ptc(:)      
  real(r8),pointer   :: h1D_CAN_G_ptc(:)      
  real(r8),pointer   :: h1D_CAN_TEMPC_ptc(:)  
  real(r8),pointer   :: h1D_CAN_TEMPFN_ptc(:)  
  real(r8),pointer   :: h1D_CAN_CO2_FLX_ptc(:) 
  real(r8),pointer   :: h1D_CAN_GPP_ptc(:)     
  real(r8),pointer   :: h1D_CAN_cumGPP_ptc(:)       
  real(r8),pointer   :: h1D_dCAN_GPP_CLIM_ptc(:)
  real(r8),pointer   :: h1D_dCAN_GPP_eLIM_ptc(:)  
  real(r8),pointer   :: h1D_CAN_RA_ptc(:)      
  real(r8),pointer   :: h1D_CAN_GROWTH_ptc(:)
  real(r8),pointer   :: h1D_cTNC_ptc(:)        
  real(r8),pointer   :: h1D_cTNN_ptc(:)       
  real(r8),pointer   :: h1D_cTNP_ptc(:)    
  real(r8),pointer   :: h1D_CanNonstBConc_ptc(:)      
  real(r8),pointer   :: h1D_STOML_RSC_CO2_ptc(:)
  real(r8),pointer   :: h1D_STOML_Min_RSC_CO2_ptc(:)
  real(r8),pointer   :: h1D_Km_CO2_carboxy_ptc(:)
  real(r8),pointer   :: h1D_Ci_mesophyll_ptc(:)
  real(r8),pointer   :: h1D_BLYR_RSC_CO2_ptc(:) 
  real(r8),pointer   :: h1D_CAN_CO2_ptc(:)      
  real(r8),pointer   :: h1D_O2L_ptc(:)
  real(r8),pointer   :: h1D_LAI_ptc(:)          
  real(r8),pointer   :: h1D_PSI_CAN_ptc(:)      
  real(r8),pointer   :: h1D_TURG_CAN_ptc(:)     
  real(r8),pointer   :: h1D_STOML_RSC_H2O_ptc(:) 
  real(r8),pointer   :: h1D_BLYR_RSC_H2O_ptc(:) 
  real(r8),pointer   :: h1D_TRANSPN_ptc(:)      
  real(r8),pointer   :: h1D_NH4_UPTK_FLX_ptc(:)  
  real(r8),pointer   :: h1D_NO3_UPTK_FLX_ptc(:)  
  real(r8),pointer   :: h1D_N2_FIXN_FLX_ptc(:)   
  real(r8),pointer   :: h1D_cNH3_FLX_ptc(:)      
  real(r8),pointer   :: h1D_PO4_UPTK_FLX_ptc(:)   
  real(r8),pointer   :: h1D_TC_Groth_ptc(:)
  real(r8),pointer   :: h2D_SoilRest4RootGroth_vr(:,:)
  real(r8),pointer   :: h2D_RootNonstBConc_pvr(:,:)   
  real(r8),pointer   :: h2D_Root1stStrutC_pvr(:,:) 
  real(r8),pointer   :: h2D_Root1stStrutN_pvr(:,:)
  real(r8),pointer   :: h2D_Root1stStrutP_pvr(:,:)
  real(r8),pointer   :: h2D_Root2ndStrutC_pvr(:,:)
  real(r8),pointer   :: h2D_Root2ndStrutN_pvr(:,:)
  real(r8),pointer   :: h2D_Root2ndStrutP_pvr(:,:)
  real(r8),pointer   :: h2D_QDrainloss_vr(:,:)  
  real(r8),pointer   :: h1D_SHOOT_C_ptc(:)      
  real(r8),pointer   :: h1D_RCanMaintDef_CO2_pft(:)           
  real(r8),pointer   :: h1D_LEAF_C_ptc(:)       
  real(r8),pointer   :: h1D_Petole_C_ptc(:)     
  real(r8),pointer   :: h1D_STALK_C_ptc(:)      
  real(r8),pointer   :: h1D_RESERVE_C_ptc(:)    
  real(r8),pointer   :: h1D_HUSK_C_ptc(:)       
  real(r8),pointer   :: h1D_GRAIN_C_ptc(:)      
  real(r8),pointer   :: h1D_ROOT_NONSTC_ptc(:)
  real(r8),pointer   :: h1D_ROOT_NONSTN_ptc(:)
  real(r8),pointer   :: h1D_ROOT_NONSTP_ptc(:)
  real(r8),pointer   :: h1D_ROOT_C_ptc(:)       
  real(r8),pointer   :: h1D_ROOTST_C_ptc(:)     
  real(r8),pointer   :: h1D_ROOTST_N_ptc(:) 
  real(r8),pointer   :: h1D_ROOTST_P_ptc(:) 
  real(r8),pointer   :: h1D_ShootNodule_C_ptc(:)     
  real(r8),pointer   :: h1D_ShootNodule_N_ptc(:)     
  real(r8),pointer   :: h1D_ShootNodule_P_ptc(:)             
  real(r8),pointer   :: h1D_RootNodule_C_ptc(:)     
  real(r8),pointer   :: h1D_STORED_C_ptc(:)     
  real(r8),pointer   :: h1D_GRAIN_NO_ptc(:)     
  real(r8),pointer   :: h1D_CondGasXSurf_col(:)
  real(r8),pointer   :: h1D_LAIb_ptc(:)         
  real(r8),pointer   :: h1D_EXUD_CumYr_C_FLX_ptc(:)   
  real(r8),pointer   :: h1D_LITRf_C_FLX_ptc(:)     
  real(r8),pointer   :: h1D_LITRf_P_FLX_ptc(:)     
  real(r8),pointer   :: h1D_SURF_LITRf_C_FLX_ptc(:)
  real(r8),pointer   :: h1D_AUTO_RESP_FLX_ptc(:)   
  real(r8),pointer   :: h1D_HVST_C_FLX_ptc(:)      
  real(r8),pointer   :: h1D_PLANT_BALANCE_C_ptc(:) 
  real(r8),pointer   :: h1D_STANDING_DEAD_C_ptc(:) 
  real(r8),pointer   :: h1D_FIREp_CO2_FLX_ptc(:)   
  real(r8),pointer   :: h1D_FIREp_CH4_FLX_ptc(:)   
  real(r8),pointer   :: h1D_cNPP_ptc(:)        
  real(r8),pointer   :: h1D_CAN_HT_ptc(:)     
  real(r8),pointer   :: h1D_POPN_ptc(:)       
  real(r8),pointer   :: h1D_tTRANSPN_ptc(:)   
  real(r8),pointer   :: h1D_WTR_STRESS_ptc(:) 
  real(r8),pointer   :: h1D_OXY_STRESS_ptc(:) 
  real(r8),pointer   :: h1D_SHOOT_N_ptc(:)    
  real(r8),pointer   :: h1D_LEAF_N_ptc(:)  
  real(r8),pointer   :: h1D_fCNLFW_ptc(:)  
  real(r8),pointer   :: h1D_fCPLFW_ptc(:)    
  real(r8),pointer   :: h1D_LEAFN2LAI_ptc(:)
  real(r8),pointer   :: h1D_Petole_N_ptc(:)   
  real(r8),pointer   :: h1D_STALK_N_ptc(:)    
  real(r8),pointer   :: h1D_RESERVE_N_ptc(:)  
  real(r8),pointer   :: h1D_HUSK_N_ptc(:)     
  real(r8),pointer   :: h1D_GRAIN_N_ptc(:)    
  real(r8),pointer   :: h1D_ROOT_N_ptc(:)     
  real(r8),pointer   :: h1D_RootNodule_N_ptc(:)   
  real(r8),pointer   :: h1D_STORED_N_ptc(:)   
  real(r8),pointer   :: h1D_EXUD_N_FLX_ptc(:)     
  real(r8),pointer   :: h1D_Uptk_N_Flx_ptc(:)
  real(r8),pointer   :: h1D_Uptk_P_Flx_ptc(:)
  real(r8),pointer   :: h1D_LITRf_N_FLX_ptc(:)    
  real(r8),pointer   :: h1D_cum_N_FIXED_ptc(:) 
  real(r8),pointer   :: h1D_TreeRingRadius_ptc(:)
  real(r8),pointer   :: h1D_HVST_N_FLX_ptc(:)     
  real(r8),pointer   :: h1D_NH3can_FLX_ptc(:)   
  real(r8),pointer   :: h1D_PLANT_BALANCE_N_ptc(:)  
  real(r8),pointer   :: h1D_STANDING_DEAD_N_ptc(:)  
  real(r8),pointer   :: h1D_FIREp_N_FLX_ptc(:)      
  real(r8),pointer   :: h1D_VcMaxRubisco_ptc(:)
  real(r8),pointer   :: h1D_LeafProteinCperm2_ptc(:)
  real(r8),pointer   :: h1D_PARSunlit_ptc(:)
  real(r8),pointer   :: h1D_PARSunsha_ptc(:)
  real(r8),pointer   :: h1D_CH2OSunlit_ptc(:)
  real(r8),pointer   :: h1D_CH2OSunsha_ptc(:)  
  real(r8),pointer   :: h1D_VoMaxRubisco_ptc(:)  
  real(r8),pointer   :: h1D_VcMaxPEP_ptc(:)  
  real(r8),pointer   :: h1D_JMaxPhoto_ptc(:)
  real(r8),pointer   :: h1D_TFN_Carboxy_ptc(:)
  real(r8),pointer   :: h1D_SLA_ptc(:)
  real(r8),pointer   :: h1D_TFN_Oxygen_ptc(:)
  real(r8),pointer   :: h1D_TFN_eTranspt_ptc(:)  
  real(r8),pointer   :: h1D_LeafAreaSunlit_ptc(:)
  real(r8),pointer   :: h1D_fClump_ptc(:)
  real(r8),pointer   :: h1D_SURF_LITRf_N_FLX_ptc(:) 
  real(r8),pointer   :: h1D_SHOOT_P_ptc(:)     
  real(r8),pointer   :: h1D_LEAF_P_ptc(:)      
  real(r8),pointer   :: h1D_Petole_P_ptc(:)    
  real(r8),pointer   :: h1D_STALK_P_ptc(:)     
  real(r8),pointer   :: h1D_RESERVE_P_ptc(:)  
  real(r8),pointer   :: h1D_HUSK_P_ptc(:)     
  real(r8),pointer   :: h1D_GRAIN_P_ptc(:)    
  real(r8),pointer   :: h1D_ROOT_P_ptc(:)     
  real(r8),pointer   :: h1D_RootNodule_P_ptc(:)   
  real(r8),pointer   :: h1D_STORED_P_ptc(:)    
  real(r8),pointer   :: h1D_EXUD_P_FLX_ptc(:)  
  real(r8),pointer   :: h1D_LITTERf_P_ptc(:)   
  real(r8),pointer   :: h1D_HVST_P_FLX_ptc(:)    
  real(r8),pointer   :: h1D_PLANT_BALANCE_P_ptc(:)   
  real(r8),pointer   :: h1D_STANDING_DEAD_P_ptc(:)   
  real(r8),pointer   :: h1D_FIREp_P_FLX_ptc(:)       
  real(r8),pointer   :: h1D_SURF_LITRf_P_FLX_ptc(:) 
  real(r8),pointer   :: h1D_ShootRootXferC_ptc(:)
  real(r8),pointer   :: h1D_ShootRootXferN_ptc(:)
  real(r8),pointer   :: h1D_ShootRootXferP_ptc(:)    
  real(r8),pointer   :: h2D_RootMaintDef_CO2_pvr(:,:)
  real(r8),pointer   :: h2D_RootSurfAreaPP_pvr(:,:)
  real(r8),pointer   :: h2D_ROOTNLim_rpvr(:,:)
  real(r8),pointer   :: h2D_ROOTPLim_rpvr(:,:)  
  real(r8),pointer   :: h1D_RootMaintDef_CO2_pft(:)
  real(r8),pointer   :: h1D_BRANCH_NO_ptc(:)    
  real(r8),pointer   :: h1D_SHOOT_NONSTC_ptc(:) 
  real(r8),pointer   :: h1D_SHOOT_NONSTN_ptc(:)  
  real(r8),pointer   :: h1D_SHOOT_NONSTP_ptc(:) 
  real(r8),pointer   :: h1D_LeafC3ChlCperm2LA_ptc(:) 
  real(r8),pointer   :: h1D_LeafC4ChlCperm2LA_ptc(:) 
  real(r8),pointer   :: h1D_LeafRubiscoCperm2LA_ptc(:)
  real(r8),pointer   :: h1D_LeafPEPCperm2LA_ptc(:)    
!  real(r8),pointer   :: h1D_CFIX_lmtf_ptc(:)
  real(r8),pointer   :: h1D_MainBranchNO_ptc(:)
  real(r8),pointer   :: h1D_Ar_mass_col(:)       
  real(r8),pointer   :: h1D_CO2_mass_col(:)
  real(r8),pointer   :: h1D_Gchem_CO2_prod_col(:)
  real(r8),pointer   :: h1D_LEAF_NC_ptc(:)      
  real(r8),pointer   :: h1D_Growth_Stage_ptc(:) 
  real(r8),pointer   :: h2D_LEAF_NODE_NO_ptc(:,:) 
  real(r8),pointer   :: h1D_RUB_ACTVN_ptc(:)    
  real(r8),pointer   :: h1D_CanopyNLim_ptc(:)
  real(r8),pointer   :: h1D_CanopyPLim_ptc(:)
  real(r8),pointer   :: h3D_PARTS_ptc(:,:,:)      
  real(r8),pointer   :: h2D_Gas_Pressure_vr(:,:)
  real(r8),pointer   :: h2D_N2O_Gas_ppmv_vr(:,:)  
  real(r8),pointer   :: h2D_CO2_Gas_ppmv_vr(:,:)
  real(r8),pointer   :: h2D_CH4_Gas_ppmv_vr(:,:)
  real(r8),pointer   :: h2D_H2_Gas_ppmv_vr(:,:)
  real(r8),pointer   :: h2D_N2_Gas_ppmv_vr(:,:)  
  real(r8),pointer   :: h2D_Ar_Gas_ppmv_vr(:,:)
  real(r8),pointer   :: h2D_O2_Gas_ppmv_vr(:,:)
  real(r8),pointer   :: h2D_NH3_Gas_ppmv_vr(:,:)
  real(r8),pointer   :: h2D_AeroHrBactC_vr(:,:)   
  real(r8),pointer   :: h2D_AeroHrFungC_vr(:,:)   
  real(r8),pointer   :: h2D_faculDenitC_vr(:,:)  
  real(r8),pointer   :: h2D_fermentorC_vr(:,:)  
  real(r8),pointer   :: h2D_acetometgC_vr(:,:)  
  real(r8),pointer   :: h2D_fermentor_frac_vr(:,:)
  real(r8),pointer   :: h2D_acetometh_frac_vr(:,:)  
  real(r8),pointer   :: h2D_hydrogMeth_frac_vr(:,:)
  real(r8),pointer   :: h2D_aeroN2fixC_vr(:,:)  
  real(r8),pointer   :: h2D_anaeN2FixC_vr(:,:)  
  real(r8),pointer   :: h2D_NH3OxiBactC_vr(:,:)
  real(r8),pointer   :: h2D_NO2OxiBactC_vr(:,:)
  real(r8),pointer   :: h2D_CH4AeroOxiC_vr(:,:)
  real(r8),pointer   :: h2D_H2MethogenC_vr(:,:)
  real(r8),pointer   :: h2D_tOMActCDens_vr(:,:)
  real(r8),pointer   :: h1D_tOMActCDens_litr_col(:)  
  real(r8),pointer   :: h2D_TSolidOMActC_vr(:,:)
  real(r8),pointer   :: h2D_TSolidOMActCDens_vr(:,:)
  real(r8),pointer   :: h2D_RCH4ProdHydrog_vr(:,:)
  real(r8),pointer   :: h2D_RCH4ProdAcetcl_vr(:,:)
  real(r8),pointer   :: h2D_RCH4Oxi_aero_vr(:,:)
  real(r8),pointer   :: h1D_RCH4Oxi_aero_col(:)
  real(r8),pointer   :: h2D_RCH4Oxi_anmo_vr(:,:)
  real(r8),pointer   :: h1D_RCH4Oxi_anmo_col(:)
  real(r8),pointer   :: h2D_RFerment_vr(:,:)
  real(r8),pointer   :: h2D_nh3oxi_vr(:,:)
  real(r8),pointer   :: h2D_n2oprod_vr(:,:)
  real(r8),pointer   :: h2D_RootMassC_vr(:,:)
  real(r8),pointer   :: h2D_RootMassC_pvr(:,:)
  real(r8),pointer   :: h2D_RootRadialKond2H2O_pvr(:,:)
  real(r8),pointer   :: h2D_RootAxialKond2H2O_pvr(:,:)
  real(r8),pointer   :: h2D_VmaxNH4Root_pvr(:,:)
  real(r8),pointer   :: h2D_VmaxNO3Root_pvr(:,:)
  real(r8),pointer   :: h2D_RootMassN_vr(:,:)
  real(r8),pointer   :: h2D_RootMassP_vr(:,:)
  real(r8),pointer   :: h1D_TSolidOMActC_litr_col(:)  
  real(r8),pointer   :: h1D_TSolidOMActCDens_litr_col(:)
  real(r8),pointer   :: h1D_RCH4ProdHydrog_litr_col(:)
  real(r8),pointer   :: h1D_RCH4ProdAcetcl_litr_col(:)  
  real(r8),pointer   :: h1D_RCH4Oxi_aero_litr_col(:)
  real(r8),pointer   :: h1D_RCH4Oxi_anmo_litr_col(:)  
  real(r8),pointer   :: h1D_RFermen_litr_col(:)
  real(r8),pointer   :: h1D_NH3oxi_litr_col(:)
  real(r8),pointer   :: h1D_N2oprod_litr_col(:)
  real(r8),pointer   :: h1D_RDECOMPC_SOM_litr_col(:)
  real(r8),pointer   :: h1D_MicrobAct_litr_col(:)
  real(r8),pointer   :: h1D_RDECOMPC_BReSOM_litr_col(:)
  real(r8),pointer   :: h1D_RDECOMPC_SorpSOM_litr_col(:)
  real(r8),pointer   :: h1D_tRespGrossHete_litr_col(:)
  real(r8),pointer   :: h1D_tRespGrossHeteUlm_litr_col(:)

  real(r8),pointer   :: h2D_AeroHrBactN_vr(:,:)  
  real(r8),pointer   :: h2D_AeroHrFungN_vr(:,:)  
  real(r8),pointer   :: h2D_faculDenitN_vr(:,:)  
  real(r8),pointer   :: h2D_fermentorN_vr(:,:) 
  real(r8),pointer   :: h2D_acetometgN_vr(:,:) 
  real(r8),pointer   :: h2D_aeroN2fixN_vr(:,:) 
  real(r8),pointer   :: h2D_anaeN2FixN_vr(:,:) 
  real(r8),pointer   :: h2D_NH3OxiBactN_vr(:,:)
  real(r8),pointer   :: h2D_NO2OxiBactN_vr(:,:)
  real(r8),pointer   :: h2D_CH4AeroOxiN_vr(:,:)
  real(r8),pointer   :: h2D_H2MethogenN_vr(:,:)
  real(r8),pointer   :: h2D_tRespGrossHeter_vr(:,:)
  real(r8),pointer   :: h2D_tRespGrossHeterUlm_vr(:,:)
  real(r8),pointer   :: h2D_FermOXYI_vr(:,:)
  real(r8),pointer   :: h2D_RDECOMPC_SOM_vr(:,:)
  real(r8),pointer   :: h2D_RDECOMPC_BReSOM_vr(:,:)
  real(r8),pointer   :: h2D_RDECOMPC_SorpSOM_vr(:,:)
  real(r8),pointer   :: h2D_MicrobAct_vr(:,:)
  real(r8),pointer   :: h3D_SOMHydrylScalCps_vr(:,:,:)
  real(r8),pointer   :: h3D_MicrobActCps_vr(:,:,:)
  real(r8),pointer   :: h3D_HydrolCSOMCps_vr(:,:,:)
  real(r8),pointer   :: h3D_Aerobic_Bacteria_vr(:,:,:)
  real(r8),pointer   :: h3D_Aerobic_Fungi_vr(:,:,:)
  real(r8),pointer   :: h3D_Facult_Denitf_vr(:,:,:)
  real(r8),pointer   :: h2D_AeroHrBactP_vr(:,:)   
  real(r8),pointer   :: h2D_AeroHrFungP_vr(:,:)   
  real(r8),pointer   :: h2D_faculDenitP_vr(:,:)  
  real(r8),pointer   :: h2D_fermentorP_vr(:,:)  
  real(r8),pointer   :: h2D_acetometgP_vr(:,:)  
  real(r8),pointer   :: h2D_aeroN2fixP_vr(:,:)  
  real(r8),pointer   :: h2D_anaeN2FixP_vr(:,:)  
  real(r8),pointer   :: h2D_NH3OxiBactP_vr(:,:)
  real(r8),pointer   :: h2D_NO2OxiBactP_vr(:,:)
  real(r8),pointer   :: h2D_CH4AeroOxiP_vr(:,:)
  real(r8),pointer   :: h2D_H2MethogenP_vr(:,:)
  real(r8),pointer   :: h2D_CanopyLAIZ_plyr(:,:)
  real(r8),pointer   :: h2D_MicroBiomeE_litr_col(:,:) 
  real(r8),pointer   :: h2D_AeroHrBactE_litr_col(:,:) 
  real(r8),pointer   :: h2D_AeroHrFungE_litr_col(:,:) 
  real(r8),pointer   :: h2D_faculDenitE_litr_col(:,:) 
  real(r8),pointer   :: h2D_fermentorE_litr_col(:,:) 
  real(r8),pointer   :: h2D_acetometgE_litr_col(:,:) 
  real(r8),pointer   :: h2D_aeroN2fixE_litr_col(:,:) 
  real(r8),pointer   :: h2D_anaeN2FixE_litr_col(:,:)  
  real(r8),pointer   :: h2D_NH3OxiBactE_litr_col(:,:)
  real(r8),pointer   :: h2D_NO2OxiBactE_litr_col(:,:)
  real(r8),pointer   :: h2D_CH4AeroOxiE_litr_col(:,:)
  real(r8),pointer   :: h2D_H2MethogenE_litr_col(:,:)
  real(r8),pointer   :: h2D_O2_rootconduct_pvr(:,:)
  real(r8),pointer   :: h2D_CO2_rootconduct_pvr(:,:)
  real(r8),pointer   :: h2D_ProteinNperm2LeafArea_pnd(:,:) !by node
  real(r8),pointer   :: h2D_RNITRIF_vr(:,:)
  real(r8),pointer   :: h2D_Root_CO2_vr(:,:)
  real(r8),pointer   :: h2D_Aqua_CO2_vr(:,:)     
  real(r8),pointer   :: h2D_Aqua_CH4_vr(:,:)     
  real(r8),pointer   :: h2D_Aqua_O2_vr(:,:)      
  real(r8),pointer   :: h2D_Aqua_N2_vr(:,:)        
  real(r8),pointer   :: h2D_Aqua_H2_vr(:,:)        
  real(r8),pointer   :: h2D_Aqua_Ar_vr(:,:)        
  real(r8),pointer   :: h2D_Aqua_N2O_vr(:,:)     
  real(r8),pointer   :: h2D_Aqua_NH3_vr(:,:)     
  real(r8),pointer   :: h2D_RootAR2soil_vr(:,:)  
  real(r8),pointer   :: h2D_RootAR2Root_vr(:,:)
  real(r8),pointer   :: h2D_TEMP_vr(:,:)      
  real(r8),pointer   :: h2D_decomp_OStress_vr(:,:) 
  real(r8),pointer   :: h2D_RO2Decomp_vr(:,:)
  real(r8),pointer   :: h2D_Decomp_temp_FN_vr(:,:)
  real(r8),pointer   :: h1D_Decomp_temp_FN_litr_col(:)
  real(r8),pointer   :: h1D_FracLitMix_litr_col(:)
  real(r8),pointer   :: h2D_FracLitMix_vr(:,:)
  real(r8),pointer   :: h2D_Decomp_Moist_FN_vr(:,:)
  real(r8),pointer   :: h1D_Decomp_Moist_FN_litr_col(:)
  real(r8),pointer   :: h1D_decomp_OStress_litr_col(:)
  real(r8),pointer   :: h1D_RO2Decomp_litr_col(:)  
  real(r8),pointer   :: h2D_HeatUptk_vr(:,:)   
  real(r8),pointer   :: h2D_HeatFlow_vr(:,:)   
  real(r8),pointer   :: h2D_VSPore_vr(:,:)
  real(r8),pointer   :: h2D_FLO_MICP_vr(:,:)
  real(r8),pointer   :: h2D_FLO_MACP_vr(:,:)
  real(r8),pointer   :: h2D_rVSM_vr(:,:)    
  real(r8),pointer   :: h2D_rVSICE_vr(:,:)   
  real(r8),pointer   :: h2D_PSI_vr(:,:)        
  real(r8),pointer   :: h2D_PsiO_vr(:,:)
  real(r8),pointer   :: h2D_cNH4t_vr(:,:)                                                     
  real(r8),pointer   :: h2D_RootH2OUP_vr(:,:)  
  real(r8),pointer   :: h2D_cNO3t_vr(:,:)                                                     
  real(r8),pointer   :: h2D_cPO4_vr(:,:)       
  real(r8),pointer   :: h2D_cEXCH_P_vr(:,:)    
  real(r8),pointer   :: h2D_microb_N2fix_vr(:,:)
  real(r8),pointer   :: h2D_ElectricConductivity_vr(:,:) 
  real(r8),pointer   :: h2D_HydCondSoil_vr(:,:)
  real(r8),pointer   :: h2D_PSI_RT_pvr(:,:)     
  real(r8),pointer   :: h2D_RootH2OUptkStress_pvr(:,:)
  real(r8),pointer   :: h2D_RootH2OUptk_pvr(:,:)
  real(r8),pointer   :: h2D_ROOT_OSTRESS_pvr(:,:)
  real(r8),pointer   :: h2D_prtUP_NH4_pvr(:,:)                                                                     
  real(r8),pointer   :: h2D_prtUP_NO3_pvr(:,:)  
  real(r8),pointer   :: h2D_prtUP_PO4_pvr(:,:)                                                  
  real(r8),pointer   :: h2D_DNS_RT_pvr(:,:)   
  real(r8),pointer   :: h2D_RootNutupk_fClim_pvr(:,:)
  real(r8),pointer   :: h2D_RootNutupk_fNlim_pvr(:,:)
  real(r8),pointer   :: h2D_RootNutupk_fPlim_pvr(:,:)
  real(r8),pointer   :: h2D_RootNutupk_fProtC_pvr(:,:)
  real(r8),pointer   :: h2D_Root1stAxesNumL_pvr(:,:)
  real(r8),pointer   :: h2D_Root2ndAxesNumL_pvr(:,:)
  REAL(R8),pointer   :: h2D_RootKond2H2O_pvr(:,:)
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
  integer :: nbr,jj
  character(len=32) :: fieldname

  beg_col=1;end_col=bounds%ncols
  beg_ptc=1;end_ptc=bounds%npfts
  allocate(this%h1D_cumFIRE_CO2_col(beg_col:end_col)) ;this%h1D_cumFIRE_CO2_col(:)=spval    
  allocate(this%h1D_cumFIRE_CH4_col(beg_col:end_col)) ;this%h1D_cumFIRE_CH4_col(:)=spval    
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
  allocate(this%h1D_RadSW_Grnd_col(beg_col:end_col)); this%h1D_RadSW_Grnd_col(:)=spval
  allocate(this%h1D_CanSWRad_col(beg_col:end_col)); this%h1D_CanSWRad_col(:)=spval
  allocate(this%h1D_Qinfl2soi_col(beg_col:end_col))     ;this%h1D_Qinfl2soi_col(:)=spval
  allocate(this%h1D_Qdrain_col(beg_col:end_col))       ; this%h1D_Qdrain_col(:)=spval
  allocate(this%h2D_QDrainloss_vr(beg_col:end_col,1:JZ)); this%h2D_QDrainloss_vr=spval  
  allocate(this%h2D_FermOXYI_vr(beg_col:end_col,1:JZ)); this%h2D_FermOXYI_vr=spval
  allocate(this%h1D_tSALT_DISCHG_FLX_col(beg_col:end_col)) ;this%h1D_tSALT_DISCHG_FLX_col(:)=spval
  allocate(this%h1D_SUR_DON_FLX_col(beg_col:end_col))    ;this%h1D_SUR_DON_FLX_col(:)=spval
  allocate(this%h1D_SUB_DON_FLX_col(beg_col:end_col))    ;this%h1D_SUB_DON_FLX_col(:)=spval
  allocate(this%h1D_SUR_DIN_FLX_col(beg_col:end_col))    ;this%h1D_SUR_DIN_FLX_col(:)=spval
  allocate(this%h1D_SUB_DIN_FLX_col(beg_col:end_col))    ;this%h1D_SUB_DIN_FLX_col(:)=spval 
  allocate(this%h1D_SUR_DOC_FLX_col(beg_col:end_col))    ;this%h1D_SUR_DOC_FLX_col(:)=spval
  allocate(this%h1D_SUB_DOC_FLX_col(beg_col:end_col))    ;this%h1D_SUB_DOC_FLX_col(:)=spval 
  allocate(this%h1D_SUR_DIC_FLX_col(beg_col:end_col))    ;this%h1D_SUR_DIC_FLX_col(:)=spval 
  allocate(this%h1D_SUB_DIC_FLX_col(beg_col:end_col))    ;this%h1D_SUB_DIC_FLX_col(:)=spval 

  allocate(this%h1D_tPREC_P_col(beg_col:end_col))      ;this%h1D_tPREC_P_col(:)=spval 
  allocate(this%h1D_tMICRO_P_col(beg_col:end_col))       ;this%h1D_tMICRO_P_col(:)=spval
  allocate(this%h1D_PO4_FIRE_col(beg_col:end_col))       ;this%h1D_PO4_FIRE_col(:)=spval
  allocate(this%h1D_cPO4_LITR_col(beg_col:end_col))      ;this%h1D_cPO4_LITR_col(:)=spval
  allocate(this%h1D_cEXCH_P_LITR_col(beg_col:end_col))   ;this%h1D_cEXCH_P_LITR_col(:)=spval
  allocate(this%h1D_NET_P_MIN_col(beg_col:end_col))      ;this%h1D_NET_P_MIN_col(:)=spval
  allocate(this%h1D_N2_SEMIS_FLX_col(beg_col:end_col))   ;this%h1D_N2_SEMIS_FLX_col(:)=spval
  allocate(this%h1D_NH3_SEMIS_FLX_col(beg_col:end_col))  ;this%h1D_NH3_SEMIS_FLX_col(:)=spval
  allocate(this%h1D_H2_SEMIS_FLX_col(beg_col:end_col))   ;this%h1D_H2_SEMIS_FLX_col(:)=spval
  allocate(this%h1D_PSI_SURF_col(beg_col:end_col))       ;this%h1D_PSI_SURF_col(:)=spval
  allocate(this%h1D_SURF_ELEV_col(beg_col:end_col))      ;this%h1D_SURF_ELEV_col(:)=spval
  allocate(this%h1D_tLITR_C_col(beg_col:end_col));this%h1D_tLITR_C_col(:)=spval
  allocate(this%h1D_tLITR_N_col(beg_col:end_col));this%h1D_tLITR_N_col(:)=spval
  allocate(this%h1D_RootAR_col(beg_col:end_col)); this%h1D_RootAR_col(:)=spval
  allocate(this%h1D_RootCO2Relez_col(beg_col:end_col)); this%h1D_RootCO2Relez_col(:)=spval
  allocate(this%h1D_tRAD_col(beg_col:end_col)); this%h1D_tRAD_col(:)=spval
  allocate(this%h2D_BotDEPZ_vr(beg_col:end_col,1:JZ)); this%h2D_BotDEPZ_vr(:,:)=spval
  allocate(this%h1D_AMENDED_N_col(beg_col:end_col))       ;this%h1D_AMENDED_N_col(:)=spval
  allocate(this%h1D_tNH4X_col(beg_col:end_col))           ;this%h1D_tNH4X_col(:)=spval
  allocate(this%h1D_tNO3_col(beg_col:end_col))            ;this%h1D_tNO3_col(:)=spval
  allocate(this%h1D_tMICRO_N_col(beg_col:end_col))        ;this%h1D_tMICRO_N_col(:)=spval
  allocate(this%h1D_TEMP_LITR_col(beg_col:end_col))       ;this%h1D_TEMP_LITR_col(:)=spval
  allocate(this%h1D_TEMP_surf_col(beg_col:end_col))       ;this%h1D_TEMP_surf_col(:)=spval
  allocate(this%h1D_TEMP_SNOW_col(beg_col:end_col))       ;this%h1D_TEMP_SNOW_col(:)=spval
  allocate(this%h1D_FracBySnow_col(beg_col:end_col))      ;this%h1D_FracBySnow_col(:)=spval
  allocate(this%h1D_FracByLitr_col(beg_col:end_col))      ;this%h1D_FracByLitr_col(:)=spval
  allocate(this%h1D_AMENDED_C_col(beg_col:end_col))       ;this%h1D_AMENDED_C_col(:)=spval
  allocate(this%h1D_tMICRO_C_col(beg_col:end_col))        ;this%h1D_tMICRO_C_col(:)=spval
  allocate(this%h1D_tSoilOrgC_col(beg_col:end_col))       ;this%h1D_tSoilOrgC_col(:)=spval
  allocate(this%h1D_tSoilOrgN_col(beg_col:end_col))       ;this%h1D_tSoilOrgN_col(:)=spval
  allocate(this%h1D_tSoilOrgP_col(beg_col:end_col))       ;this%h1D_tSoilOrgP_col(:)=spval
  allocate(this%h1D_OMC_LITR_col(beg_col:end_col))        ;this%h1D_OMC_LITR_col(:)=spval
  allocate(this%h1D_OMN_LITR_col(beg_col:end_col))        ;this%h1D_OMN_LITR_col(:)=spval
  allocate(this%h1D_OMP_LITR_col(beg_col:end_col))        ;this%h1D_OMP_LITR_col(:)=spval
  allocate(this%h1D_Ar_mass_col(beg_col:end_col))         ;this%h1D_Ar_mass_col(:) = spval
  allocate(this%h1D_CO2_mass_col(beg_col:end_col))         ;this%h1D_CO2_mass_col(:) = spval
  allocate(this%h1D_Gchem_CO2_prod_col(beg_col:end_col)) ; this%h1D_Gchem_CO2_prod_col(:)=spval
  allocate(this%h1D_Ar_soilMass_col(beg_col:end_col))     ;this%h1D_Ar_soilMass_col(:)=spval
  allocate(this%h1D_ATM_CO2_col(beg_col:end_col))         ;this%h1D_ATM_CO2_col(:)=spval
  allocate(this%h1D_ATM_CH4_col(beg_col:end_col))         ;this%h1D_ATM_CH4_col(:)=spval
  allocate(this%h1D_NBP_col(beg_col:end_col))             ;this%h1D_NBP_col(:)=spval
  allocate(this%h1D_CO2_WetDep_FLX_col(beg_col:end_col)) ; this%h1D_CO2_WetDep_FLX_col(:)=spval
  allocate(this%h1D_RootN_Fix_col(beg_col:end_col)); this%h1D_RootN_Fix_col(:)=spval
  allocate(this%h1D_Ar_WetDep_FLX_col(beg_col:end_col)) ; this%h1D_Ar_WetDep_FLX_col(:)=spval
  allocate(this%h1D_RootXO2_flx_col(beg_col:end_col)); this%h1D_RootXO2_flx_col(:)=spval
  allocate(this%h1D_RUNOFF_FLX_col(beg_col:end_col)) ; this%h1D_RUNOFF_FLX_col(:)=spval
  allocate(this%h1D_SEDIMENT_FLX_col(beg_col:end_col))    ;this%h1D_SEDIMENT_FLX_col(:)=spval
  allocate(this%h1D_QDISCHG_FLX_col(beg_col:end_col))      ;this%h1D_QDISCHG_FLX_col(:)=spval
  allocate(this%h1D_HeatDISCHG_FLX_col(beg_col:end_col))  ; this%h1D_HeatDISCHG_FLX_col(:)=spval
  allocate(this%h1D_ECO_LAI_col(beg_col:end_col))         ;this%h1D_ECO_LAI_col(:)=spval
  allocate(this%h1D_ECO_SAI_col(beg_col:end_col))         ;this%h1D_ECO_SAI_col(:)=spval
  allocate(this%h1D_Eco_GPP_CumYr_col(beg_col:end_col))         ;this%h1D_Eco_GPP_CumYr_col(:)=spval
  allocate(this%h1D_ECO_RA_col(beg_col:end_col))          ;this%h1D_ECO_RA_col(:)=spval
  allocate(this%h1D_Eco_NPP_CumYr_col(beg_col:end_col))         ;this%h1D_Eco_NPP_CumYr_col(:)=spval
  allocate(this%h1D_Eco_HR_CumYr_col(beg_col:end_col))          ;this%h1D_Eco_HR_CumYr_col(:)=spval
  allocate(this%h1D_Eco_HR_CO2_col(beg_col:end_col));  this%h1D_Eco_HR_CO2_col(:)=spval
  allocate(this%h1D_Eco_HR_CO2_litr_col(beg_col:end_col));  this%h1D_Eco_HR_CO2_litr_col(:)=spval
!  allocate(this%h1D_Eco_HR_CH4_col(beg_col:end_col));  this%h1D_Eco_HR_CH4_col(:)=spval
  allocate(this%h1D_tDIC_col(beg_col:end_col))            ;this%h1D_tDIC_col=spval
  allocate(this%h1D_tSTANDING_DEAD_C_col(beg_col:end_col));this%h1D_tSTANDING_DEAD_C_col=spval 
  allocate(this%h1D_tSTANDING_DEAD_N_col(beg_col:end_col));this%h1D_tSTANDING_DEAD_N_col=spval  
  allocate(this%h1D_tSTANDING_DEAD_P_col(beg_col:end_col));this%h1D_tSTANDING_DEAD_P_col=spval  
  allocate(this%h1D_tPRECIP_col(beg_col:end_col))          ;this%h1D_tPRECIP_col(:)=spval
  allocate(this%h1D_ECO_ET_col(beg_col:end_col))              ;this%h1D_ECO_ET_col(:)=spval
  allocate(this%h1D_trcg_NH3_cumerr_col(beg_col:end_col)); this%h1D_trcg_NH3_cumerr_col(:)=spval
  allocate(this%h1D_trcg_H2_cumerr_col(beg_col:end_col)); this%h1D_trcg_H2_cumerr_col(:)=spval
  allocate(this%h1D_trcg_Ar_cumerr_col(beg_col:end_col)); this%h1D_trcg_Ar_cumerr_col(:)=spval
  allocate(this%h1D_trcg_N2_cumerr_col(beg_col:end_col)); this%h1D_trcg_N2_cumerr_col(:)=spval
  allocate(this%h1D_trcg_O2_cumerr_col(beg_col:end_col)); this%h1D_trcg_O2_cumerr_col(:)=spval
  allocate(this%h1D_trcg_CO2_cumerr_col(beg_col:end_col)); this%h1D_trcg_CO2_cumerr_col(:)=spval
  allocate(this%h1D_trcg_CH4_cumerr_col(beg_col:end_col)); this%h1D_trcg_CH4_cumerr_col(:)=spval
  allocate(this%h1d_CAN_NEE_col(beg_col:end_col))        ; this%h1d_CAN_NEE_col(:)=spval
  allocate(this%h1D_ECO_RADSW_col(beg_col:end_col))       ; this%h1D_ECO_RADSW_col(:)=spval
  allocate(this%h1D_N2O_LITR_col(beg_col:end_col))        ;this%h1D_N2O_LITR_col(:)=spval
  allocate(this%h1D_NH3_LITR_col(beg_col:end_col))        ;this%h1D_NH3_LITR_col(:)=spval
  allocate(this%h1D_SOL_RADN_col(beg_col:end_col))        ;this%h1D_SOL_RADN_col(:)=spval
  allocate(this%h1D_AIR_TEMP_col(beg_col:end_col))        ;this%h1D_AIR_TEMP_col(:)=spval
  allocate(this%h1D_FreeNFix_col(beg_col:end_col))       ;this%h1D_FreeNFix_col(:)=spval
  allocate(this%h1D_PATM_col(beg_col:end_col))            ;this%h1D_PATM_col(:)=spval
  allocate(this%h1D_HUM_col(beg_col:end_col))             ;this%h1D_HUM_col(:)=spval
  allocate(this%h1D_WIND_col(beg_col:end_col))            ;this%h1D_WIND_col(:)=spval
  allocate(this%h1D_PREC_col(beg_col:end_col))            ;this%h1D_PREC_col(:)=spval
  allocate(this%h1D_Snofall_col(beg_col:end_col))         ;this%h1D_Snofall_col(:)=spval
  allocate(this%h1D_SOIL_RN_col(beg_col:end_col))         ;this%h1D_SOIL_RN_col(:)=spval
  allocate(this%h1D_LWSky_col(beg_col:end_col))           ;this%h1D_LWSky_col(:)=spval
  allocate(this%h1D_SOIL_LE_col(beg_col:end_col))         ;this%h1D_SOIL_LE_col(:)=spval
  allocate(this%h1D_SOIL_H_col(beg_col:end_col))          ;this%h1D_SOIL_H_col(:)=spval
  allocate(this%h1D_SOIL_G_col(beg_col:end_col))          ;this%h1D_SOIL_G_col(:)=spval
  allocate(this%h1D_ECO_RN_col(beg_col:end_col))          ;this%h1D_ECO_RN_col(:)=spval
  allocate(this%h1D_ECO_LE_col(beg_col:end_col))          ;this%h1D_ECO_LE_col(:)=spval
  allocate(this%h1D_Eco_HeatSen_col(beg_col:end_col))        ;this%h1D_Eco_HeatSen_col(:)=spval
  allocate(this%h1D_ECO_Heat2G_col(beg_col:end_col))           ;this%h1D_ECO_Heat2G_col(:)=spval
  allocate(this%h1D_O2_LITR_col(beg_col:end_col))         ;this%h1D_O2_LITR_col(:)=spval
  allocate(this%h1D_MIN_LWP_ptc(beg_ptc:end_ptc))         ;this%h1D_MIN_LWP_ptc(:)=spval
  allocate(this%h1D_SLA_ptc(beg_ptc:end_ptc)); this%h1D_SLA_ptc(:)=spval
  allocate(this%h1D_CO2_SEMIS_FLX_col(beg_col:end_col))    ;this%h1D_CO2_SEMIS_FLX_col(:)=spval
  allocate(this%h1D_AR_SEMIS_FLX_col(beg_col:end_col))     ;this%h1D_AR_SEMIS_FLX_col(:)=spval
  allocate(this%h1D_ECO_CO2_FLX_col(beg_col:end_col))     ;this%h1D_ECO_CO2_FLX_col(:)=spval
  allocate(this%h1D_CH4_SEMIS_FLX_col(beg_col:end_col))         ;this%h1D_CH4_SEMIS_FLX_col(:)=spval
  allocate(this%h1D_CH4_EBU_flx_col(beg_col:end_col))     ;this%h1D_CH4_EBU_flx_col(:)=spval
  allocate(this%h1D_Ar_EBU_flx_col(beg_col:end_col))      ;this%h1D_Ar_EBU_flx_col(:)=spval
  allocate(this%h1D_CO2_TPR_err_col(beg_col:end_col))    ; this%h1D_CO2_TPR_err_col(:)=spval
  allocate(this%h1D_Ar_TPR_err_col(beg_col:end_col))     ;this%h1D_Ar_TPR_err_col(:)=spval
  allocate(this%h1D_CO2_Drain_flx_col(beg_col:end_col))  ; this%h1D_CO2_Drain_flx_col(:)=spval
  allocate(this%h1D_CO2_hydloss_flx_col(beg_col:end_col));  this%h1D_CO2_hydloss_flx_col(:)=spval
  allocate(this%h1D_CH4_PLTROOT_flx_col(beg_col:end_col)) ;this%h1D_CH4_PLTROOT_flx_col(:)=spval
  allocate(this%h1D_AR_PLTROOT_flx_col(beg_col:end_col)); this%h1D_AR_PLTROOT_flx_col(:)=spval
  allocate(this%h1D_CO2_PLTROOT_flx_col(beg_col:end_col)) ;this%h1D_CO2_PLTROOT_flx_col(:)=spval
  allocate(this%h1D_O2_PLTROOT_flx_col(beg_col:end_col)) ;this%h1D_O2_PLTROOT_flx_col(:)=spval
  allocate(this%h1D_CO2_DIF_flx_col(beg_col:end_col)); this%h1D_CO2_DIF_flx_col(:)=spval
  allocate(this%h1D_O2_DIF_flx_col(beg_col:end_col)); this%h1D_O2_DIF_flx_col(:)=spval
  allocate(this%h1D_CH4_DIF_flx_col(beg_col:end_col)); this%h1D_CH4_DIF_flx_col(:)=spval  
  allocate(this%h1D_NH3_DIF_flx_col(beg_col:end_col)); this%h1D_NH3_DIF_flx_col(:)=spval
  allocate(this%h1D_Ar_DIF_flx_col(beg_col:end_col)); this%h1D_Ar_DIF_flx_col(:)=spval
  allocate(this%h1D_O2_SEMIS_FLX_col(beg_col:end_col))          ;this%h1D_O2_SEMIS_FLX_col(:)=spval
  allocate(this%h1D_CO2_LITR_col(beg_col:end_col))        ;this%h1D_CO2_LITR_col(:)=spval
  allocate(this%h1D_EVAPN_col(beg_col:end_col))           ;this%h1D_EVAPN_col(:)=spval
  allocate(this%h1D_CondGasXSurf_col(beg_col:end_col)); this%h1D_CondGasXSurf_col(:)=spval
  allocate(this%h1D_CANET_col(beg_col:end_col))           ;this%h1D_CANET_col(:)=spval
  allocate(this%h1D_PAR_col(beg_col:end_col))             ;this%h1D_PAR_col(:)=spval
  allocate(this%h1d_fPAR_col(beg_col:end_col));   this%h1d_fPAR_col(:)=spval
  allocate(this%h1D_tSWC_col(beg_col:end_col))            ;this%h1D_tSWC_col(:)=spval
  allocate(this%h1D_tHeat_col(beg_col:end_col))           ;this%h1D_tHeat_col(:)=spval
  allocate(this%h1D_SNOWPACK_col(beg_col:end_col))        ;this%h1D_SNOWPACK_col(:)=spval
  allocate(this%h1D_SNOWDENS_col(beg_col:end_col))        ;this%h1D_SNOWDENS_col(:)=spval
  allocate(this%h1D_SURF_WTR_col(beg_col:end_col))        ;this%h1D_SURF_WTR_col(:)=spval
  allocate(this%h1D_ThetaW_litr_col(beg_col:end_col))    ;this%h1D_ThetaW_litr_col(:)=spval
  allocate(this%h1D_ThetaI_litr_col(beg_col:end_col))    ;this%h1D_ThetaI_litr_col(:)=spval
  allocate(this%h1D_SURF_ICE_col(beg_col:end_col))        ;this%h1D_SURF_ICE_col(:)=spval
  allocate(this%h1D_ACTV_LYR_col(beg_col:end_col))        ;this%h1D_ACTV_LYR_col(:)=spval
  allocate(this%h1D_WTR_TBL_col(beg_col:end_col))         ;this%h1D_WTR_TBL_col(:)=spval
  allocate(this%h1D_N2O_SEMIS_FLX_col(beg_col:end_col))        ;this%h1D_N2O_SEMIS_FLX_col(:)=spval
  allocate(this%h1D_CAN_G_ptc(beg_ptc:end_ptc))           ;this%h1D_CAN_G_ptc(:)=spval
  allocate(this%h1D_CAN_TEMPC_ptc(beg_ptc:end_ptc))        ;this%h1D_CAN_TEMPC_ptc(:)=spval
  allocate(this%h1D_CAN_TEMPFN_ptc(beg_ptc:end_ptc))      ;this%h1D_CAN_TEMPFN_ptc(:)=spval
  allocate(this%h1D_CAN_CO2_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_CAN_CO2_FLX_ptc(:)=spval
  allocate(this%h1D_CAN_GPP_ptc(beg_ptc:end_ptc))         ;this%h1D_CAN_GPP_ptc(:)=spval
  allocate(this%h1D_CAN_cumGPP_ptc(beg_ptc:end_ptc))         ;this%h1D_CAN_cumGPP_ptc(:)=spval  
  allocate(this%h1D_dCAN_GPP_CLIM_ptc(beg_ptc:end_ptc))    ;this%h1D_dCAN_GPP_CLIM_ptc(:)=spval
  allocate(this%h1D_dCAN_GPP_eLIM_ptc(beg_ptc:end_ptc))    ;this%h1D_dCAN_GPP_eLIM_ptc(:)=spval
  allocate(this%h1D_CAN_RA_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_RA_ptc(:)=spval
  allocate(this%h1D_LEAF_PC_ptc(beg_ptc:end_ptc))         ;this%h1D_LEAF_PC_ptc(:)=spval
  allocate(this%h1D_CAN_RN_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_RN_ptc(:)=spval
  allocate(this%h1D_CAN_LE_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_LE_ptc(:)=spval
  allocate(this%h1D_CAN_H_ptc(beg_ptc:end_ptc))           ;this%h1D_CAN_H_ptc(:)=spval  
  allocate(this%h1D_CAN_GROWTH_ptc(beg_ptc:end_ptc))      ;this%h1D_CAN_GROWTH_ptc(:)=spval
  allocate(this%h1D_cTNC_ptc(beg_ptc:end_ptc))            ;this%h1D_cTNC_ptc(:)=spval
  allocate(this%h1D_cTNN_ptc(beg_ptc:end_ptc))            ;this%h1D_cTNN_ptc(:)=spval
  allocate(this%h1D_cTNP_ptc(beg_ptc:end_ptc))            ;this%h1D_cTNP_ptc(:)=spval
  allocate(this%h1D_CanNonstBConc_ptc(beg_ptc:end_ptc)) ;this%h1D_CanNonstBConc_ptc(:)=spval
  allocate(this%h1D_STOML_RSC_CO2_ptc(beg_ptc:end_ptc))   ;this%h1D_STOML_RSC_CO2_ptc(:)=spval
  allocate(this%h1D_STOML_Min_RSC_CO2_ptc(beg_ptc:end_ptc));this%h1D_STOML_Min_RSC_CO2_ptc(:)=spval
  allocate(this%h1D_Km_CO2_carboxy_ptc(beg_ptc:end_ptc)); this%h1D_Km_CO2_carboxy_ptc(:)=spval
  allocate(this%h1D_Ci_mesophyll_ptc(beg_ptc:end_ptc)); this%h1D_Ci_mesophyll_ptc(:)=spval
  allocate(this%h1D_BLYR_RSC_CO2_ptc(beg_ptc:end_ptc))    ;this%h1D_BLYR_RSC_CO2_ptc(:)=spval
  allocate(this%h1D_CAN_CO2_ptc(beg_ptc:end_ptc))         ;this%h1D_CAN_CO2_ptc(:)=spval
  allocate(this%h1D_O2L_ptc(beg_ptc:end_ptc)); this%h1D_O2L_ptc(:)=spval
  allocate(this%h1D_LAI_ptc(beg_ptc:end_ptc))             ;this%h1D_LAI_ptc(:)=spval
  allocate(this%h1D_PSI_CAN_ptc(beg_ptc:end_ptc))         ;this%h1D_PSI_CAN_ptc(:)=spval
  allocate(this%h1D_RootAR_ptc(beg_ptc:end_ptc))          ;this%h1D_RootAR_ptc(:)=spval
  allocate(this%h1D_RootLenPerPlant_ptc(beg_ptc:end_ptc))         ;this%h1D_RootLenPerPlant_ptc(:)=spval
  allocate(this%h1D_TURG_CAN_ptc(beg_ptc:end_ptc))        ;this%h1D_TURG_CAN_ptc(:)=spval
  allocate(this%h1D_STOML_RSC_H2O_ptc(beg_ptc:end_ptc))    ;this%h1D_STOML_RSC_H2O_ptc(:)=spval
  allocate(this%h1D_BLYR_RSC_H2O_ptc(beg_ptc:end_ptc))    ;this%h1D_BLYR_RSC_H2O_ptc(:)=spval
  allocate(this%h1D_TRANSPN_ptc(beg_ptc:end_ptc))         ;this%h1D_TRANSPN_ptc(:)=spval
  allocate(this%h1D_QTRANSP_col(beg_col:end_col)); this%h1D_QTRANSP_col(:)=spval
  allocate(this%h1D_NH4_UPTK_FLX_ptc(beg_ptc:end_ptc))    ;this%h1D_NH4_UPTK_FLX_ptc(:)=spval
  allocate(this%h1D_NO3_UPTK_FLX_ptc(beg_ptc:end_ptc))    ;this%h1D_NO3_UPTK_FLX_ptc(:)=spval
  allocate(this%h1D_N2_FIXN_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_N2_FIXN_FLX_ptc(:)=spval
  allocate(this%h1D_cNH3_FLX_ptc(beg_ptc:end_ptc))        ;this%h1D_cNH3_FLX_ptc(:)=spval
  allocate(this%h1D_TC_Groth_ptc(beg_ptc:end_ptc))        ;this%h1D_TC_Groth_ptc(:)=spval
  allocate(this%h1D_PO4_UPTK_FLX_ptc(beg_ptc:end_ptc))    ;this%h1D_PO4_UPTK_FLX_ptc(:)=spval
  allocate(this%h1D_SHOOT_C_ptc(beg_ptc:end_ptc))         ;this%h1D_SHOOT_C_ptc(:)=spval

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
  allocate(this%h1D_ShootNodule_C_ptc(beg_ptc:end_ptc))    ;this%h1D_ShootNodule_C_ptc(:)=spval
  allocate(this%h1D_ShootNodule_N_ptc(beg_ptc:end_ptc))    ;this%h1D_ShootNodule_N_ptc(:)=spval
  allocate(this%h1D_ShootNodule_P_ptc(beg_ptc:end_ptc))    ;this%h1D_ShootNodule_P_ptc(:)=spval      
  allocate(this%h1D_RootNodule_C_ptc(beg_ptc:end_ptc))    ;this%h1D_RootNodule_C_ptc(:)=spval
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
  allocate(this%h1D_cNPP_ptc(beg_ptc:end_ptc))             ;this%h1D_cNPP_ptc(:)=spval
  allocate(this%h1D_CAN_HT_ptc(beg_ptc:end_ptc))          ;this%h1D_CAN_HT_ptc(:)=spval
  allocate(this%h1D_POPN_ptc(beg_ptc:end_ptc))            ;this%h1D_POPN_ptc(:)=spval
  allocate(this%h1D_tTRANSPN_ptc(beg_ptc:end_ptc))        ;this%h1D_tTRANSPN_ptc(:)=spval
  allocate(this%h1D_WTR_STRESS_ptc(beg_ptc:end_ptc))      ;this%h1D_WTR_STRESS_ptc(:)=spval
  allocate(this%h1D_OXY_STRESS_ptc(beg_ptc:end_ptc))      ;this%h1D_OXY_STRESS_ptc(:)=spval
  allocate(this%h1D_VcMaxRubisco_ptc(beg_ptc:end_ptc))    ;this%h1D_VcMaxRubisco_ptc(:)=spval
  allocate(this%h1D_LeafProteinCperm2_ptc(beg_ptc:end_ptc));this%h1D_LeafProteinCperm2_ptc(:)=spval
  allocate(this%h1D_VoMaxRubisco_ptc(beg_ptc:end_ptc))    ;this%h1D_VoMaxRubisco_ptc(:)=spval
  allocate(this%h1D_VcMaxPEP_ptc(beg_ptc:end_ptc))        ;this%h1D_VcMaxPEP_ptc(:)=spval
  allocate(this%h1D_JMaxPhoto_ptc(beg_ptc:end_ptc))       ;this%h1D_JMaxPhoto_ptc(:)=spval
  allocate(this%h1D_TFN_Carboxy_ptc(beg_ptc:end_ptc))     ;this%h1D_TFN_Carboxy_ptc(:)=spval
  allocate(this%h1D_TFN_Oxygen_ptc(beg_ptc:end_ptc))      ;this%h1D_TFN_Oxygen_ptc(:)=spval
  allocate(this%h1D_TFN_eTranspt_ptc(beg_ptc:end_ptc))    ;this%h1D_TFN_eTranspt_ptc(:)=spval

  allocate(this%h1D_PARSunlit_ptc(beg_ptc:end_ptc))  ;this%h1D_PARSunlit_ptc(:)=spval
  allocate(this%h1D_PARSunsha_ptc(beg_ptc:end_ptc))  ;this%h1D_PARSunsha_ptc(:)=spval
  allocate(this%h1D_CH2OSunlit_ptc(beg_ptc:end_ptc)) ;this%h1D_CH2OSunlit_ptc(:)=spval
  allocate(this%h1D_CH2OSunsha_ptc(beg_ptc:end_ptc)); this%h1D_CH2OSunsha_ptc(:)=spval
  allocate(this%h1D_fClump_ptc(beg_ptc:end_ptc)); this%h1D_fClump_ptc(:)=spval
  allocate(this%h1D_LeafAreaSunlit_ptc(beg_ptc:end_ptc))  ;this%h1D_LeafAreaSunlit_ptc(:)=spval  
  allocate(this%h1D_SHOOT_N_ptc(beg_ptc:end_ptc))         ;this%h1D_SHOOT_N_ptc(:)=spval
  allocate(this%h1D_Plant_N_ptc(beg_ptc:end_ptc))         ;this%h1D_Plant_N_ptc(:)=spval
  allocate(this%h1D_LEAF_N_ptc(beg_ptc:end_ptc))          ;this%h1D_LEAF_N_ptc(:)=spval
  allocate(this%h1D_fCNLFW_ptc(beg_ptc:end_ptc)) ; this%h1D_fCNLFW_ptc(:)=spval
  allocate(this%h1D_fCPLFW_ptc(beg_ptc:end_ptc)); this%h1D_fCPLFW_ptc(:)=spval
  allocate(this%h1D_LEAFN2LAI_ptc(beg_ptc:end_ptc))      ;this%h1D_LEAFN2LAI_ptc(:)=spval
  allocate(this%h1D_Petole_N_ptc(beg_ptc:end_ptc))       ;this%h1D_Petole_N_ptc(:)=spval
  allocate(this%h1D_STALK_N_ptc(beg_ptc:end_ptc))         ;this%h1D_STALK_N_ptc(:)=spval
  allocate(this%h1D_RESERVE_N_ptc(beg_ptc:end_ptc))       ;this%h1D_RESERVE_N_ptc(:)=spval
  allocate(this%h1D_HUSK_N_ptc(beg_ptc:end_ptc))          ;this%h1D_HUSK_N_ptc(:)=spval
  allocate(this%h1D_GRAIN_N_ptc(beg_ptc:end_ptc))         ;this%h1D_GRAIN_N_ptc(:)=spval
  allocate(this%h1D_ROOT_N_ptc(beg_ptc:end_ptc))          ;this%h1D_ROOT_N_ptc(:)=spval
  allocate(this%h1D_RootNodule_N_ptc(beg_ptc:end_ptc))        ;this%h1D_RootNodule_N_ptc(:)=spval
  allocate(this%h1D_STORED_N_ptc(beg_ptc:end_ptc))        ;this%h1D_STORED_N_ptc(:)=spval
  allocate(this%h1D_EXUD_N_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_EXUD_N_FLX_ptc(:)=spval
  allocate(this%h1D_Uptk_N_Flx_ptc(beg_ptc:end_ptc))      ;this%h1D_Uptk_N_Flx_ptc(:)=spval
  allocate(this%h1D_Uptk_P_Flx_ptc(beg_ptc:end_ptc))      ;this%h1D_Uptk_P_Flx_ptc(:)=spval
  allocate(this%h1D_LITRf_N_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_LITRf_N_FLX_ptc(:)=spval
  allocate(this%h1D_cum_N_FIXED_ptc(beg_ptc:end_ptc))  ;this%h1D_cum_N_FIXED_ptc(:)=spval
  allocate(this%h1D_TreeRingRadius_ptc(beg_ptc:end_ptc));this%h1D_TreeRingRadius_ptc(:)=spval
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
  allocate(this%h1D_RootNodule_P_ptc(beg_ptc:end_ptc))        ;this%h1D_RootNodule_P_ptc(:)=spval
  allocate(this%h1D_STORED_P_ptc(beg_ptc:end_ptc))        ;this%h1D_STORED_P_ptc(:)=spval
  allocate(this%h1D_EXUD_P_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_EXUD_P_FLX_ptc(:)=spval
  allocate(this%h1D_LITTERf_P_ptc(beg_ptc:end_ptc))       ;this%h1D_LITTERf_P_ptc(:)=spval
  allocate(this%h1D_HVST_P_FLX_ptc(beg_ptc:end_ptc))      ;this%h1D_HVST_P_FLX_ptc(:)=spval
  allocate(this%h1D_STANDING_DEAD_P_ptc(beg_ptc:end_ptc)) ;this%h1D_STANDING_DEAD_P_ptc(:)=spval
  allocate(this%h1D_FIREp_P_FLX_ptc(beg_ptc:end_ptc))     ;this%h1D_FIREp_P_FLX_ptc(:)=spval
  allocate(this%h1D_SURF_LITRf_P_FLX_ptc(beg_ptc:end_ptc));this%h1D_SURF_LITRf_P_FLX_ptc(:)=spval
  allocate(this%h1D_ShootRootXferC_ptc(beg_ptc:end_ptc)); this%h1D_ShootRootXferC_ptc(:)=spval
  allocate(this%h1D_ShootRootXferN_ptc(beg_ptc:end_ptc)); this%h1D_ShootRootXferN_ptc(:)=spval
  allocate(this%h1D_ShootRootXferP_ptc(beg_ptc:end_ptc)); this%h1D_ShootRootXferP_ptc(:)=spval    
  allocate(this%h1D_BRANCH_NO_ptc(beg_ptc:end_ptc))       ;this%h1D_BRANCH_NO_ptc(:)=spval
  allocate(this%h1D_Growth_Stage_ptc(beg_ptc:end_ptc))    ;this%h1D_Growth_Stage_ptc(:)=spval
  allocate(this%h1D_LEAF_NC_ptc(beg_ptc:end_ptc));        ;this%h1D_LEAF_NC_ptc(:)=spval
  allocate(this%h1D_SHOOT_NONSTC_ptc(beg_ptc:end_ptc))     ;this%h1D_SHOOT_NONSTC_ptc(:)=spval
  allocate(this%h1D_SHOOT_NONSTN_ptc(beg_ptc:end_ptc))     ;this%h1D_SHOOT_NONSTN_ptc(:)=spval
  allocate(this%h1D_SHOOT_NONSTP_ptc(beg_ptc:end_ptc))     ;this%h1D_SHOOT_NONSTP_ptc(:)=spval

  allocate(this%h1D_LeafC3ChlCperm2LA_ptc(beg_ptc:end_ptc));this%h1D_LeafC3ChlCperm2LA_ptc=spval
  allocate(this%h1D_LeafC4ChlCperm2LA_ptc(beg_ptc:end_ptc));this%h1D_LeafC4ChlCperm2LA_ptc=spval 
  allocate(this%h1D_LeafRubiscoCperm2LA_ptc(beg_ptc:end_ptc));this%h1D_LeafRubiscoCperm2LA_ptc=spval
  allocate(this%h1D_LeafPEPCperm2LA_ptc(beg_ptc:end_ptc));this%h1D_LeafPEPCperm2LA_ptc=spval    
!  allocate(this%h1D_CFIX_lmtf_ptc(beg_ptc:end_ptc))        ;this%h1D_CFIX_lmtf_ptc(:)=spval
  allocate(this%h1D_MainBranchNO_ptc(beg_ptc:end_ptc))     ;this%h1D_MainBranchNO_ptc(:)=ispval
  allocate(this%h1D_RCanMaintDef_CO2_pft(beg_ptc:end_ptc)) ;this%h1D_RCanMaintDef_CO2_pft(:)=spval
  allocate(this%h1D_RootMaintDef_CO2_pft(beg_ptc:end_ptc)) ;this%h1D_RootMaintDef_CO2_pft(:)=spval
  allocate(this%h1D_ROOT_NONSTC_ptc(beg_ptc:end_ptc))     ;this%h1D_ROOT_NONSTC_ptc(:)=spval
  allocate(this%h1D_ROOT_NONSTN_ptc(beg_ptc:end_ptc))     ;this%h1D_ROOT_NONSTN_ptc(:)=spval
  allocate(this%h1D_ROOT_NONSTP_ptc(beg_ptc:end_ptc))     ;this%h1D_ROOT_NONSTP_ptc(:)=spval
  allocate(this%h2D_tSOC_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSOC_vr(:,:)=spval
  allocate(this%h2D_POM_C_vr(beg_col:end_col,1:JZ)) ; this%h2D_POM_C_vr(:,:)=spval
  allocate(this%h2D_MAOM_C_vr(beg_col:end_col,1:JZ)); this%h2D_MAOM_C_vr(:,:)=spval
  allocate(this%h2D_microbC_vr(beg_col:end_col,1:JZ))  ;this%h2D_microbC_vr(:,:)=spval
  allocate(this%h2D_microbN_vr(beg_col:end_col,1:JZ))  ;this%h2D_microbN_vr(:,:)=spval
  allocate(this%h2D_microbP_vr(beg_col:end_col,1:JZ))  ;this%h2D_microbP_vr(:,:)=spval    
  allocate(this%h2D_AeroBact_PrimS_lim_vr(beg_col:end_col,1:JZ)); this%h2D_AeroBact_PrimS_lim_vr(:,:)=spval
  allocate(this%h2D_AeroFung_PrimS_lim_vr(beg_col:end_col,1:JZ)); this%h2D_AeroFung_PrimS_lim_vr(:,:)=spval  
  allocate(this%h2D_tSOCL_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSOCL_vr(:,:)=spval
  allocate(this%h2D_NO3_vr(beg_col:end_col,1:JZ)); this%h2D_NO3_vr(:,:)=spval
  allocate(this%h2D_NH4_vr(beg_col:end_col,1:JZ)); this%h2D_NH4_vr(:,:)=spval
  allocate(this%h2D_tSON_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSON_vr(:,:)=spval
  allocate(this%h2D_tSOP_vr(beg_col:end_col,1:JZ))    ;this%h2D_tSOP_vr(:,:)=spval
  allocate(this%h2D_litrC_vr(beg_col:end_col,1:JZ))   ;this%h2D_litrC_vr(:,:)=spval
  allocate(this%h2D_litrN_vr(beg_col:end_col,1:JZ))   ;this%h2D_litrN_vr(:,:)=spval
  allocate(this%h2D_litrP_vr(beg_col:end_col,1:JZ))   ;this%h2D_litrP_vr(:,:)=spval
  allocate(this%h2D_VHeatCap_vr(beg_col:end_col,1:JZ));this%h2D_VHeatCap_vr(:,:)=spval
  allocate(this%h2D_LEAF_NODE_NO_ptc(beg_ptc:end_ptc,1:MaxNumBranches));this%h2D_LEAF_NODE_NO_ptc(:,:)=spval
  allocate(this%h1D_RUB_ACTVN_ptc(beg_ptc:end_ptc));  this%h1D_RUB_ACTVN_ptc(:)=spval
  allocate(this%h1D_CanopyNLim_ptc(beg_ptc:end_ptc)); this%h1D_CanopyNLim_ptc(:)=spval
  allocate(this%h1D_CanopyPLim_ptc(beg_ptc:end_ptc)); this%h1D_CanopyPLim_ptc(:)=spval  
  allocate(this%h2D_RNITRIF_vr(beg_col:end_col,1:JZ))    ;this%h2D_RNITRIF_vr(:,:)=spval
  allocate(this%h2D_Aqua_N2_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_N2_vr(:,:)=spval  
  allocate(this%h2D_Aqua_H2_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_H2_vr(:,:)=spval  
  allocate(this%h2D_Aqua_Ar_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_Ar_vr(:,:)=spval  
  allocate(this%h2D_Aqua_CO2_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_CO2_vr(:,:)=spval
  allocate(this%h2D_Root_CO2_vr(beg_col:end_col,1:JZ))  ; this%h2D_Root_CO2_vr(:,:)=spval
  allocate(this%h2D_RootNutupk_fClim_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootNutupk_fClim_pvr=spval
  allocate(this%h2D_RootNutupk_fNlim_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootNutupk_fNlim_pvr=spval
  allocate(this%h2D_RootNutupk_fPlim_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootNutupk_fPlim_pvr=spval
  allocate(this%h2D_RootNutupk_fProtC_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootNutupk_fProtC_pvr=spval
  allocate(this%h2D_O2_rootconduct_pvr(beg_ptc:end_ptc,1:JZ))     ;this%h2D_O2_rootconduct_pvr(:,:)=spval
  allocate(this%h2D_CO2_rootconduct_pvr(beg_ptc:end_ptc,1:JZ))     ;this%h2D_CO2_rootconduct_pvr(:,:)=spval
  allocate(this%h2D_ProteinNperm2LeafArea_pnd(beg_ptc:end_ptc,1:MaxNodesPerBranch)); this%h2D_ProteinNperm2LeafArea_pnd(:,:)=spval
  allocate(this%h2D_Aqua_CH4_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_CH4_vr(:,:)=spval
  allocate(this%h2D_Aqua_O2_vr(beg_col:end_col,1:JZ))         ;this%h2D_Aqua_O2_vr(:,:)=spval
  allocate(this%h2D_Aqua_N2O_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_N2O_vr(:,:)=spval
  allocate(this%h2D_Aqua_NH3_vr(beg_col:end_col,1:JZ))        ;this%h2D_Aqua_NH3_vr(:,:)=spval
  allocate(this%h2D_TEMP_vr(beg_col:end_col,1:JZ))       ;this%h2D_TEMP_vr(:,:)=spval
  allocate(this%h1D_decomp_OStress_litr_col(beg_col:end_col)); this%h1D_decomp_OStress_litr_col(:)=spval
  allocate(this%h1D_Decomp_temp_FN_litr_col(beg_col:end_col)); this%h1D_Decomp_temp_FN_litr_col(:)=spval
  allocate(this%h1D_FracLitMix_litr_col(beg_col:end_col)); this%h1D_FracLitMix_litr_col(:)=spval
  allocate(this%h1D_Decomp_Moist_FN_litr_col(beg_col:end_col)); this%h1D_Decomp_Moist_FN_litr_col(:)=spval
  allocate(this%h1D_RO2Decomp_litr_col(beg_col:end_col)); this%h1D_RO2Decomp_litr_col(:)=spval
  allocate(this%h2D_decomp_OStress_vr(beg_col:end_col,1:JZ)); this%h2D_decomp_OStress_vr(:,:)=spval
  allocate(this%h2D_RO2Decomp_vr(beg_col:end_col,1:JZ)); this%h2D_RO2Decomp_vr(:,:)=spval
  allocate(this%h2D_Decomp_temp_FN_vr(beg_col:end_col,1:JZ)); this%h2D_Decomp_temp_FN_vr(:,:)=spval
  allocate(this%h2D_FracLitMix_vr(beg_col:end_col,1:JZ)); this%h2D_FracLitMix_vr(:,:)=spval
  allocate(this%h2D_Decomp_Moist_FN_vr(beg_col:end_col,1:JZ)); this%h2D_Decomp_Moist_FN_vr(:,:)=spval
  allocate(this%h2D_RootMassC_vr(beg_col:end_col,1:JZ)); this%h2D_RootMassC_vr(:,:)=spval
  allocate(this%h2D_RootMassN_vr(beg_col:end_col,1:JZ)); this%h2D_RootMassN_vr(:,:)=spval
  allocate(this%h2D_RootMassP_vr(beg_col:end_col,1:JZ)); this%h2D_RootMassP_vr(:,:)=spval
  allocate(this%h2D_RootMassC_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootMassC_pvr(:,:)=spval
  allocate(this%h2D_RootRadialKond2H2O_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootRadialKond2H2O_pvr(:,:)=spval
  allocate(this%h2D_RootAxialKond2H2O_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootAxialKond2H2O_pvr(:,:)=spval
  allocate(this%h2D_VmaxNH4Root_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_VmaxNH4Root_pvr(:,:)=spval
  allocate(this%h2D_VmaxNO3Root_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_VmaxNO3Root_pvr(:,:)=spval
  allocate(this%h2D_DOC_vr(beg_col:end_col,1:JZ)); this%h2D_DOC_vr(:,:)=spval
  allocate(this%h2D_DON_vr(beg_col:end_col,1:JZ)); this%h2D_DON_vr(:,:)=spval
  allocate(this%h2D_DOP_vr(beg_col:end_col,1:JZ)); this%h2D_DOP_vr(:,:)=spval
  allocate(this%h2D_SoilRest4RootGroth_vr(beg_col:end_col,1:JZ)); this%h2D_SoilRest4RootGroth_vr(:,:)=spval
  allocate(this%h2D_acetate_vr(beg_col:end_col,1:JZ)); this%h2D_acetate_vr(:,:)=spval
  allocate(this%h1D_tDOC_soil_col(beg_col:end_col)); this%h1D_tDOC_soil_col(:)=spval
  allocate(this%h1D_tDON_soil_col(beg_col:end_col)); this%h1D_tDON_soil_col(:)=spval
  allocate(this%h1D_tDOP_soil_col(beg_col:end_col)); this%h1D_tDOP_soil_col(:)=spval
  allocate(this%h1D_tAcetate_soil_col(beg_col:end_col)); this%h1D_tAcetate_soil_col(:)=spval
  allocate(this%h1D_DOC_litr_col(beg_col:end_col));  this%h1D_DOC_litr_col(:)=spval
  allocate(this%h1D_DON_litr_col(beg_col:end_col)); this%h1D_DON_litr_col(:)=spval
  allocate(this%h1D_DOP_litr_col(beg_col:end_col)); this%h1D_DOP_litr_col(:)=spval
  allocate(this%h1D_acetate_litr_col(beg_col:end_col)); this%h1D_acetate_litr_col(:)=spval  
  allocate(this%h1D_VHeatCap_litr_col(beg_col:end_col)); this%h1D_VHeatCap_litr_col(:)=spval
  allocate(this%h2D_Gas_Pressure_vr(beg_col:end_col,1:JZ)); this%h2D_Gas_Pressure_vr(:,:)=spval
  allocate(this%h2D_CO2_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_CO2_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_CH4_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_CH4_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_H2_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_H2_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_Ar_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_Ar_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_N2_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_N2_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_N2O_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_N2O_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_NH3_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_NH3_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_O2_Gas_ppmv_vr(beg_col:end_col,1:JZ)); this%h2D_O2_Gas_ppmv_vr(:,:)=spval
  allocate(this%h2D_AeroHrBactC_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrBactC_vr(:,:)=spval
  allocate(this%h2D_AeroHrFungC_vr(beg_col:end_col,1:JZ)); this%h2D_AeroHrFungC_vr(:,:)=spval
  allocate(this%h2D_faculDenitC_vr(beg_col:end_col,1:JZ)); this%h2D_faculDenitC_vr(:,:)=spval
  allocate(this%h2D_fermentorC_vr(beg_col:end_col,1:JZ));  this%h2D_fermentorC_vr(:,:)=spval

  allocate(this%h2D_fermentor_frac_vr(beg_col:end_col,1:JZ)); this%h2D_fermentor_frac_vr(:,:)=spval
  allocate(this%h2D_acetometh_frac_vr(beg_col:end_col,1:JZ)); this%h2D_acetometh_frac_vr(:,:)=spval  
  allocate(this%h2D_hydrogMeth_frac_vr(beg_col:end_col,1:JZ)); this%h2D_hydrogMeth_frac_vr(:,:)=spval
  allocate(this%h2D_HydCondSoil_vr(beg_col:end_col,1:JZ)); this%h2D_HydCondSoil_vr=spval
  allocate(this%h2D_acetometgC_vr(beg_col:end_col,1:JZ));  this%h2D_acetometgC_vr(:,:)=spval
  allocate(this%h2D_aeroN2fixC_vr(beg_col:end_col,1:JZ));  this%h2D_aeroN2fixC_vr(:,:)=spval
  allocate(this%h2D_anaeN2FixC_vr(beg_col:end_col,1:JZ));  this%h2D_anaeN2FixC_vr(:,:)=spval
  allocate(this%h2D_NH3OxiBactC_vr(beg_col:end_col,1:JZ)); this%h2D_NH3OxiBactC_vr(:,:)=spval
  allocate(this%h2D_NO2OxiBactC_vr(beg_col:end_col,1:JZ)); this%h2D_NO2OxiBactC_vr(:,:)=spval
  allocate(this%h2D_CH4AeroOxiC_vr(beg_col:end_col,1:JZ)); this%h2D_CH4AeroOxiC_vr(:,:)=spval
  allocate(this%h2D_H2MethogenC_vr(beg_col:end_col,1:JZ)); this%h2D_H2MethogenC_vr(:,:)=spval
  allocate(this%h2D_tOMActCDens_vr(beg_col:end_col,1:JZ)); this%h2D_tOMActCDens_vr(:,:)=spval
  allocate(this%h2D_TSolidOMActC_vr(beg_col:end_col,1:JZ));this%h2D_TSolidOMActC_vr(:,:)=spval
  allocate(this%h1D_tOMActCDens_litr_col(beg_col:end_col)); this%h1D_tOMActCDens_litr_col(:)=spval
  allocate(this%h1D_TSolidOMActC_litr_col(beg_col:end_col)); this%h1D_TSolidOMActC_litr_col(:)=spval
  allocate(this%h2D_TSolidOMActCDens_vr(beg_col:end_col,1:JZ));this%h2D_TSolidOMActCDens_vr(:,:)=spval
  allocate(this%h1D_TSolidOMActCDens_litr_col(beg_col:end_col)); this%h1D_TSolidOMActCDens_litr_col(:)=spval
  allocate(this%h2D_RCH4ProdHydrog_vr(beg_col:end_col,1:JZ));   this%h2D_RCH4ProdHydrog_vr(:,:)=spval
  allocate(this%h2D_RCH4ProdAcetcl_vr(beg_col:end_col,1:JZ));   this%h2D_RCH4ProdAcetcl_vr(:,:)=spval
  allocate(this%h2D_RCH4Oxi_aero_vr(beg_col:end_col,1:JZ)); this%h2D_RCH4Oxi_aero_vr(:,:)=spval
  allocate(this%h2D_RCH4Oxi_anmo_vr(beg_col:end_col,1:JZ)); this%h2D_RCH4Oxi_anmo_vr(:,:)=spval  
  allocate(this%h2D_RFerment_vr(beg_col:end_col,1:JZ)); this%h2D_RFerment_vr(:,:)=spval
  allocate(this%h2D_nh3oxi_vr(beg_col:end_col,1:JZ));  this%h2D_nh3oxi_vr(:,:)=spval
  allocate(this%h2D_n2oprod_vr(beg_col:end_col,1:JZ));  this%h2D_n2oprod_vr(:,:)=spval
  allocate(this%h2D_RootAR_vr(beg_col:end_col,1:JZ)); this%h2D_RootAR_vr(:,:)=spval
  allocate(this%h2D_RootAR2soil_vr(beg_col:end_col,1:JZ)); this%h2D_RootAR2soil_vr(:,:)=spval
  allocate(this%h2D_RootAR2Root_vr(beg_col:end_col,1:JZ)); this%h2D_RootAR2Root_vr(:,:)=spval
  allocate(this%h1D_RCH4ProdHydrog_litr_col(beg_col:end_col));  this%h1D_RCH4ProdHydrog_litr_col(:)=spval
  allocate(this%h1D_RCH4ProdAcetcl_litr_col(beg_col:end_col));  this%h1D_RCH4ProdAcetcl_litr_col(:)=spval
  allocate(this%h1D_RCH4Oxi_aero_litr_col(beg_col:end_col)); this%h1D_RCH4Oxi_aero_litr_col(:)=spval
  allocate(this%h1D_RCH4Oxi_anmo_litr_col(beg_col:end_col)); this%h1D_RCH4Oxi_anmo_litr_col(:)=spval  
  allocate(this%h1D_RCH4Oxi_aero_col(beg_col:end_col));this%h1D_RCH4Oxi_aero_col(:)=spval
  allocate(this%h1D_RCH4Oxi_anmo_col(beg_col:end_col));this%h1D_RCH4Oxi_anmo_col(:)=spval  
  allocate(this%h1D_RFermen_litr_col(beg_col:end_col));  this%h1D_RFermen_litr_col(:)=spval
  allocate(this%h1D_nh3oxi_litr_col(beg_col:end_col)); this%h1D_nh3oxi_litr_col(:)=spval
  allocate(this%h1D_n2oprod_litr_col(beg_col:end_col));  this%h1D_n2oprod_litr_col(:)=spval
  allocate(this%h2D_Gchem_CO2_prod_vr(beg_col:end_col,1:JZ)); this%h2D_Gchem_CO2_prod_vr(:,:)=spval
  allocate(this%h2D_Eco_HR_CO2_vr(beg_col:end_col,1:JZ)); this%h2D_Eco_HR_CO2_vr(:,:)=spval
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
  allocate(this%h2D_tRespGrossHeterUlm_vr(beg_col:end_col,1:JZ)); this%h2D_tRespGrossHeterUlm_vr(:,:)=spval
  allocate(this%h2D_tRespGrossHeter_vr(beg_col:end_col,1:JZ)); this%h2D_tRespGrossHeter_vr(:,:)=spval

  allocate(this%h2D_RDECOMPC_SOM_vr(beg_col:end_col,1:JZ)); this%h2D_RDECOMPC_SOM_vr(:,:)=spval
  allocate(this%h2D_RDECOMPC_BReSOM_vr(beg_col:end_col,1:JZ));this%h2D_RDECOMPC_BReSOM_vr(:,:)=spval
  allocate(this%h2D_RDECOMPC_SorpSOM_vr(beg_col:end_col,1:JZ));this%h2D_RDECOMPC_SorpSOM_vr(:,:)=spval
  allocate(this%h2D_MicrobAct_vr(beg_col:end_col,1:JZ)); this%h2D_MicrobAct_vr(:,:)=spval

  allocate(this%h3D_MicrobActCps_vr(beg_col:end_col,1:JZ,jcplx)); this%h3D_MicrobActCps_vr(:,:,:)=spval

  allocate(this%h3D_SOMHydrylScalCps_vr(beg_col:end_col,1:JZ,jcplx));this%h3D_SOMHydrylScalCps_vr(:,:,:)=spval

  allocate(this%h3D_SOC_Cps_vr(beg_col:end_col,1:JZ,jcplx)); this%h3D_SOC_Cps_vr(:,:,:)=spval

  allocate(this%h3D_HydrolCSOMCps_vr(beg_col:end_col,1:JZ,jcplx)); this%h3D_HydrolCSOMCps_vr(:,:,:)=spval

  allocate(this%h1D_RDECOMPC_SOM_litr_col(beg_col:end_col)); this%h1D_RDECOMPC_SOM_litr_col(:)=spval
  allocate(this%h1D_RDECOMPC_BReSOM_litr_col(beg_col:end_col));this%h1D_RDECOMPC_BReSOM_litr_col(:)=spval
  allocate(this%h1D_RDECOMPC_SorpSOM_litr_col(beg_col:end_col));this%h1D_RDECOMPC_SorpSOM_litr_col(:)=spval
  allocate(this%h1D_MicrobAct_litr_col(beg_col:end_col)); this%h1D_MicrobAct_litr_col(:)=spval
  allocate(this%h1D_tRespGrossHeteUlm_litr_col(beg_col:end_col)); this%h1D_tRespGrossHeteUlm_litr_col=spval
  allocate(this%h1D_tRespGrossHete_litr_col(beg_col:end_col)); this%h1D_tRespGrossHete_litr_col=spval

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
  allocate(this%h2D_rVSM_vr    (beg_col:end_col,1:JZ))     ;this%h2D_rVSM_vr    (:,:)=spval
  allocate(this%h2D_VSPore_vr (beg_col:end_col,1:JZ))     ;this%h2D_VSPore_vr (:,:)=spval
  allocate(this%h2D_FLO_MICP_vr(beg_col:end_col,1:JZ)) ;this%h2D_FLO_MICP_vr(:,:)=spval
  allocate(this%h2D_FLO_MACP_vr(beg_col:end_col,1:JZ)) ;this%h2D_FLO_MACP_vr(:,:)=spval  
  allocate(this%h2D_rVSICE_vr    (beg_col:end_col,1:JZ))       ;this%h2D_rVSICE_vr    (:,:)=spval
  allocate(this%h2D_PSI_vr(beg_col:end_col,1:JZ))        ;this%h2D_PSI_vr(:,:)=spval
  allocate(this%h2D_PsiO_vr(beg_col:end_col,1:JZ)) ; this%h2D_PsiO_vr(:,:)=spval
  allocate(this%h2D_RootH2OUP_vr(beg_col:end_col,1:JZ))  ;this%h2D_RootH2OUP_vr(:,:)=spval
  allocate(this%h2D_cNH4t_vr(beg_col:end_col,1:JZ))      ;this%h2D_cNH4t_vr(:,:)=spval
  allocate(this%h2D_cNO3t_vr(beg_col:end_col,1:JZ))      ;this%h2D_cNO3t_vr(:,:)=spval
                                                             
  allocate(this%h2D_cPO4_vr(beg_col:end_col,1:JZ))       ;this%h2D_cPO4_vr(:,:)=spval
  allocate(this%h2D_cEXCH_P_vr(beg_col:end_col,1:JZ))    ;this%h2D_cEXCH_P_vr(:,:)=spval
  allocate(this%h2D_microb_N2fix_vr(beg_col:end_col,1:JZ)); this%h2D_microb_N2fix_vr(:,:)=spval
  allocate(this%h2D_ElectricConductivity_vr(beg_col:end_col,1:JZ))       ;this%h2D_ElectricConductivity_vr(:,:)=spval
  allocate(this%h2D_PSI_RT_pvr(beg_ptc:end_ptc,1:JZ))     ;this%h2D_PSI_RT_pvr(:,:)=spval
  allocate(this%h2D_RootH2OUptkStress_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootH2OUptkStress_pvr(:,:)=spval
  allocate(this%h2D_RootH2OUptk_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_RootH2OUptk_pvr(:,:)=spval
  allocate(this%h2D_RootMaintDef_CO2_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootMaintDef_CO2_pvr(:,:)=spval
  allocate(this%h2D_RootSurfAreaPP_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootSurfAreaPP_pvr(:,:)=spval
  allocate(this%h2D_ROOTNLim_rpvr(beg_ptc:end_ptc,1:JZ)); this%h2D_ROOTNLim_rpvr(:,:)=spval
  allocate(this%h2D_ROOTPLim_rpvr(beg_ptc:end_ptc,1:JZ)); this%h2D_ROOTPLim_rpvr(:,:)=spval
  allocate(this%h2D_ROOT_OSTRESS_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_ROOT_OSTRESS_pvr(:,:)=spval
  allocate(this%h2D_prtUP_NH4_pvr(beg_ptc:end_ptc,1:JZ))  ;this%h2D_prtUP_NH4_pvr(:,:)=spval                                                              
  allocate(this%h2D_prtUP_NO3_pvr(beg_ptc:end_ptc,1:JZ))  ;this%h2D_prtUP_NO3_pvr(:,:)=spval                                                              
  allocate(this%h2D_prtUP_PO4_pvr(beg_ptc:end_ptc,1:JZ))  ;this%h2D_prtUP_PO4_pvr(:,:)=spval                                                              
  allocate(this%h2D_DNS_RT_pvr(beg_ptc:end_ptc,1:JZ))     ;this%h2D_DNS_RT_pvr(:,:)=spval
  allocate(this%h2D_RootNonstBConc_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootNonstBConc_pvr(:,:)=spval   
  allocate(this%h2D_Root1stStrutC_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root1stStrutC_pvr=spval
  allocate(this%h2D_Root1stStrutN_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root1stStrutN_pvr=spval
  allocate(this%h2D_Root1stStrutP_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root1stStrutP_pvr=spval
  allocate(this%h2D_Root2ndStrutC_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root2ndStrutC_pvr=spval
  allocate(this%h2D_Root2ndStrutN_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root2ndStrutN_pvr=spval
  allocate(this%h2D_Root2ndStrutP_pvr(beg_ptc:end_ptc,1:JZ)) ;this%h2D_Root2ndStrutP_pvr=spval
  allocate(this%h2D_Root2ndAxesNumL_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_Root2ndAxesNumL_pvr=spval
  allocate(this%h2D_RootKond2H2O_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_RootKond2H2O_pvr=spval
  allocate(this%h2D_Root1stAxesNumL_pvr(beg_ptc:end_ptc,1:JZ));this%h2D_Root1stAxesNumL_pvr=spval
  allocate(this%h2D_fTRootGro_pvr(beg_ptc:end_ptc,1:JZ)) ; this%h2D_fTRootGro_pvr=spval
  allocate(this%h2D_fRootGrowPSISense_pvr(beg_ptc:end_ptc,1:JZ)); this%h2D_fRootGrowPSISense_pvr=spval
  allocate(this%h3D_PARTS_ptc(beg_ptc:end_ptc,1:NumOfPlantMorphUnits,1:MaxNumBranches));this%h3D_PARTS_ptc(:,:,:)=spval
  allocate(this%h2D_CanopyLAIZ_plyr(beg_ptc:end_ptc,1:NumCanopyLayers)); this%h2D_CanopyLAIZ_plyr(:,:)=spval
  !-----------------------------------------------------------------------
  ! initialize history fields 
  !--------------------------------------------------------------------
  data1d_ptr => this%h1D_cumFIRE_CO2_col(beg_col:end_col) 
  call hist_addfld1d(fname='cumFIRE_CO2_col',units='gC m-2',avgflag='I',&
    long_name='cumulative CO2 flux from fire (<0 into atmosphere)',ptr_col=data1d_ptr,default='inactive')        

  data1d_ptr => this%h1D_cumFIRE_CH4_col(beg_col:end_col)  
  call hist_addfld1d(fname='cumFIRE_CH4_col',units='gC d-2',avgflag='I', &
    long_name='cumulative CH4 flux from fire (<0 into atmosphere)', ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_cNH4_LITR_col(beg_col:end_col) 
  call hist_addfld1d(fname='cNH4_LITR_col',units='gN NH4/g litter',avgflag='A', &
    long_name='NH4 concentration in litter',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_cNO3_LITR_col(beg_col:end_col)        
  call hist_addfld1d(fname='cNO3_LITR_col',units='gN NO3/g litter',avgflag='A',&
    long_name='NO3 concentration in litter',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_ECO_HVST_C_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_HVST_C_col',units='gC/m2',avgflag='A',&
    long_name='Harvested C',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_HVST_N_col(beg_col:end_col)     
  call hist_addfld1d(fname='ECO_HVST_N_col',units='gN/m2',avgflag='A',&
    long_name='Harvested N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_HVST_P_col(beg_col:end_col)   
  call hist_addfld1d(fname='ECO_HVST_P_col',units='gP/m2',avgflag='A',&
    long_name='Harvested P',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_NET_N_MIN_col(beg_col:end_col)  
  call hist_addfld1d(fname='NET_N_MIN_col',units='gN/m2',avgflag='I',&
    long_name='Cumulative net microbial NH4 mineralization (<0 immobilization)',&
    ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_tLITR_C_col(beg_col:end_col)   
  call hist_addfld1d(fname='tLITR_C_col',units='gC/m2',avgflag='A',&
    long_name='Column integrated total (above+belowground) litter C',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tRAD_col(beg_col:end_col)
  call hist_addfld1d(fname='RADN_col',units='MJ/m2/hr',avgflag='A',&
    long_name='Total incoming solar radiation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tLITR_N_col(beg_col:end_col)      
  call hist_addfld1d(fname='tLITR_N_col',units='gN/m2',avgflag='A',&
    long_name='Column integrated total (above+belowground) litter N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_RootAR_col(beg_col:end_col)
  call hist_addfld1d(fname='Root_AR_col',units='gC/m2/h',avgflag='A',&
    long_name='Column integrated root autotrophic respiration',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_RootCO2Relez_col(beg_col:end_col)
  call hist_addfld1d(fname='Root_CO2Relez_col',units='gC/m2/h',avgflag='A',&
    long_name='Column integrated root CO2 flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tLITR_P_col(beg_col:end_col)  
  call hist_addfld1d(fname='tLITR_P_col',units='gP/m2',avgflag='A',&
    long_name='Column integrated total (above+belowground) litter P',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_HUMUS_C_col(beg_col:end_col) 
  call hist_addfld1d(fname='HUMUS_C_col',units='gC/m2',avgflag='A',&
    long_name='colum integrated humus C',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HUMUS_N_col(beg_col:end_col)         
  call hist_addfld1d(fname='HUMUS_N_col',units='gN/m2',avgflag='A',&
    long_name='colum integrated humus N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HUMUS_P_col(beg_col:end_col)       
  call hist_addfld1d(fname='HUMUS_P_col',units='gP/m2',avgflag='A',&
    long_name='colum integrated humus P',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_AMENDED_C_col(beg_col:end_col)   
  call hist_addfld1d(fname='AMENDED_C_col',units='gC/m2',avgflag='A',&
    long_name='Column-integrated total organic fertilizer C amendment',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_AMENDED_N_col(beg_col:end_col)   
  call hist_addfld1d(fname='AMENDED_N_col',units='gN/m2',avgflag='A',&
    long_name='Column-integrated total organic fertilizer N amendment',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_AMENDED_P_col(beg_col:end_col)        
  call hist_addfld1d(fname='AMENDED_P_col',units='gP/m2',avgflag='A',&
    long_name='Column-integrated total organic fertilizer P amendment',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_tLITRf_C_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='tLITRf_C_col',units='gC/m2/hr',avgflag='A',&
    long_name='Column-integrated total (above+belowground) litrFall C flux',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tLITRf_N_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='tLITRf_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='Column-integrated total (above+belowground) litrFall N',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_tLITRf_P_FLX_col(beg_col:end_col) 
  call hist_addfld1d(fname='tLITRf_P',units='gP/m2/hr',avgflag='A',&
    long_name='Column-integrated total (above+belowground) litrFall P',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_tEXCH_PO4_col(beg_col:end_col)      
  call hist_addfld1d(fname='tEXCH_PO4_col',units='gP/m2',avgflag='A',&
    long_name='Column-integrated total bioavailable (exchangeable) mineral P: H2PO4+HPO4',&
    ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_SUR_DOC_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUR_DOC_FLX_col',units='gC/m2/hr',avgflag='A',&
    long_name='Column-integrated surface DOC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUR_DON_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUR_DON_FLX_col',units='gN/m2',avgflag='A',&
    long_name='Column-integrated surface DON flux',ptr_col=data1d_ptr)            

  data1d_ptr => this%h1D_SUR_DOP_FLX_col(beg_col:end_col) 
  call hist_addfld1d(fname='SUR_DOP_FLX',units='gP/m2',avgflag='A',&
    long_name='Column-integrated surface DOP flux',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_SUB_DOC_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUB_DOC_FLX_col',units='gC/m2/hr',avgflag='A',&
    long_name='total subsurface DOC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUB_DON_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUB_DON_FLX_col',units='gN/m2/hr',avgflag='A',&
    long_name='total subsurface DON flux',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_SUB_DOP_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUB_DOP_FLX_col',units='gP/m2/hr',avgflag='A',&
    long_name='total subsurface DOP flux',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_SUR_DIC_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUR_DIC_FLX_col',units='gC/m2/hr',avgflag='A',&
    long_name='total surface DIC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUR_DIN_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUR_DIN_FLX_col',units='gN/m2',avgflag='I',&
    long_name='cumulative total surface DIN flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUR_DIP_FLX_col(beg_col:end_col)     
  call hist_addfld1d(fname='SUR_DIP_FLX_col',units='gP/m2',avgflag='I',&
    long_name='Cumulative surface DIP flux',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_SUB_DIC_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUB_DIC_FLX_col',units='gC/m2/hr',avgflag='A',&
    long_name='total subsurface DIC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUB_DIN_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUB_DIN_FLX_col',units='gN/m2/hr',avgflag='A',&
    long_name='landscape total subsurface DIN flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SUB_DIP_FLX_col(beg_col:end_col)     
  call hist_addfld1d(fname='SUB_DIP_FLX_col',units='gP/m2/hr',avgflag='A',&
    long_name='total subsurface DIP flux',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_HeatFlx2Grnd_col(beg_col:end_col)
  call hist_addfld1d(fname='HeatFlx2Grnd_col',units='MJ/m2/hr',avgflag='A',&
    long_name='Heat flux into the ground',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_RadSW_Grnd_col(beg_col:end_col)
  call hist_addfld1d(fname='RadSW_Grnd_col',units='W/m2',avgflag='A',&
    long_name='Shortwave Radiation onto the ground',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CanSWRad_col(beg_col:end_col)
  call hist_addfld1d(fname='RadSW_Canopy_col',units='W/m2',avgflag='A',&
    long_name='Shortwave Radiation onto the grid canopy',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Qinfl2soi_col(beg_col:end_col)
  call hist_addfld1d(fname='Qinfl2soi_col',units='mm H2O/hr',avgflag='A',&
    long_name='Water flux into the ground',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_QTRANSP_col(beg_col:end_col)
  call hist_addfld1d(fname='QTransp_col',units='mm H2O/hr',avgflag='A',&
    long_name='Soil water loss through transpiration (<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Qdrain_col(beg_col:end_col)
  call hist_addfld1d(fname='Qdrain_col',units='mm H2O/hr',avgflag='A',&
    long_name='Drainage water flux out (>0) of the soil column',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Ar_mass_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_mass_col',units='g/m2',avgflag='A',&
    long_name='total Ar mass of the soil column, include that in snow and roots',&
    ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_CO2_mass_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_mass_col',units='g/m2',avgflag='A',&
    long_name='total CO2 mass of the soil column, include that in snow and roots',&
    ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_Gchem_CO2_prod_col(beg_col:end_col)
  call hist_addfld1d(fname='Gchem_CO2_prod_col',units='gC/m2',avgflag='A',&
    long_name='Column integrated CO2 production rate from geochemistry',&
    ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_Ar_soilMass_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_soil_mass_col',units='g/m2',avgflag='A',&
    long_name='total Ar mass of the soil column, excluding that in snow',&
    ptr_col=data1d_ptr,default='inactive')            

  IF(salt_model)THEN
    data1d_ptr => this%h1D_tSALT_DISCHG_FLX_col(beg_col:end_col) 
    call hist_addfld1d(fname='tSALT_DISCHG_FLX_col',units='mol/m2/hr',avgflag='A',&
      long_name='total subsurface ion flux',ptr_col=data1d_ptr)      
  endif

  data1d_ptr => this%h1D_tPREC_P_col(beg_col:end_col)  
  call hist_addfld1d(fname='tPREC_P_col',units='gP/m2',avgflag='A',&
    long_name='column integrated total soil precipited P',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_tSoilOrgC_col(beg_col:end_col)    
  call hist_addfld1d(fname='tSoilOrgC_col',units='gC/m2',avgflag='A', &
    long_name='Column-integrated total soil organic C',ptr_col=data1d_ptr,default='inactive')

  data1d_ptr => this%h1D_tSoilOrgN_col(beg_col:end_col)    
  call hist_addfld1d(fname='tSoilOrgN_col',units='gN/m2',avgflag='A', &
    long_name='Column-integrated total soil organic N',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_tSoilOrgP_col(beg_col:end_col)    
  call hist_addfld1d(fname='tSoilOrgP_col',units='gP/m2',avgflag='A', &
    long_name='Column-integrated total soil organic P',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_tMICRO_C_col(beg_col:end_col)    
  call hist_addfld1d(fname='tMICROB_C_col',units='gC/m2',avgflag='A', &
    long_name='Column-integrated micriobial C (include surface litter)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tMICRO_N_col(beg_col:end_col)    
  call hist_addfld1d(fname='tMICROB_N_col',units='gN/m2',avgflag='A', &
    long_name='Column-integrated micriobial N (include surface litter)',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_tMICRO_P_col(beg_col:end_col)    
  call hist_addfld1d(fname='tMICROB_P_col',units='gP/m2',avgflag='A', &
    long_name='Column-integrated micriobial P (include surface litter)',ptr_col=data1d_ptr, &
    default='inactive')            

  data1d_ptr => this%h1D_PO4_FIRE_col(beg_col:end_col)  
  call hist_addfld1d(fname='PO4_FIRE_col',units='gP/m2',avgflag='I',&
    long_name='Cumulative PO4 flux from fire',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_cPO4_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='cPO4_LITR_col',units='gP/g litr',avgflag='A',&
    long_name='PO4 concentration in litter',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_cEXCH_P_LITR_col(beg_col:end_col)     
  call hist_addfld1d(fname='cEXCH_P_LITR_col',units='gP/g litr',avgflag='A',&
    long_name='concentration of exchangeable inorganic P in litterr',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_NET_P_MIN_col(beg_col:end_col)    
  call hist_addfld1d(fname='NET_P_MIN_col',units='gP/m2',avgflag='I',&
    long_name='Cumulative net microbial P mineralization (<0 immobilization)',&
    ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_HUM_col(beg_col:end_col)       
  call hist_addfld1d(fname='HMAX_AIR_col',units='kPa',avgflag='X',&
    long_name='daily maximum vapor pressure',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_HUM_col(beg_col:end_col)    
  call hist_addfld1d(fname='HMIN_AIR_col',units='kPa',avgflag='M',&
    long_name='daily maximum vapor pressure',ptr_col=data1d_ptr,default='inactive')      
    
  data1d_ptr => this%h1D_PSI_SURF_col(beg_col:end_col)    
  call hist_addfld1d(fname='PSI_LITR_col',units='MPa',avgflag='A',&
    long_name='Litter layer micropore matric water potential',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_SURF_ELEV_col(beg_col:end_col)  
  call hist_addfld1d(fname='SURF_ELEV_col',units='m',avgflag='A',&
    long_name='Surface elevation, including litter layer',ptr_col=data1d_ptr,default='inactive')            
    
  data1d_ptr => this%h1D_tNH4X_col(beg_col:end_col)     
  call hist_addfld1d(fname='tNH4_col',units='gN/m2',avgflag='A', &
    long_name='Column-integrated NH4+NH3',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tNO3_col(beg_col:end_col)          
  call hist_addfld1d(fname='tNO3_col',units='gN/m2',avgflag='A',&
    long_name='Column integrated NO3+NO2 content',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_TEMP_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='TEMP_LITR_col',units='oC',avgflag='A',&
    long_name='Litter layer temperature',ptr_col=data1d_ptr)            

  data1d_ptr => this%h1D_TEMP_surf_col(beg_col:end_col)    
  call hist_addfld1d(fname='TEMP_SURF_col',units='oC',avgflag='A',&
    long_name='Ground surface temperature',ptr_col=data1d_ptr,default='inactive')           

  data1d_ptr => this%h1D_TEMP_SNOW_col(beg_col:end_col)    
  call hist_addfld1d(fname='TEMP_SNOW_col',units='oC',avgflag='A',&
    long_name='First snow layer temperature',ptr_col=data1d_ptr,default='inactive')            
    
  data1d_ptr => this%h1D_FracBySnow_col(beg_col:end_col)    
  call hist_addfld1d(fname='Frac_Snow_Ground_col',units='none',avgflag='A',&
    long_name='Fraction of ground covered by snow',ptr_col=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_FracByLitr_col(beg_col:end_col)    
  call hist_addfld1d(fname='Frac_Litr_Ground_col',units='none',avgflag='A',&
    long_name='Fraction of ground covered by litter',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_OMC_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='Surf_LitrC_col',units='gC/m2',avgflag='A',&
    long_name='Total surface litter C, including microbes',ptr_col=data1d_ptr)            

  data1d_ptr => this%h1D_OMN_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='Surf_LitrN_col',units='gN/m2',avgflag='A',&
    long_name='Total surface litter N, including microbes',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_OMP_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='Surf_LitrP_col',units='gP/m2',avgflag='A',&
    long_name='Total surface litter P, including microbes',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_ATM_CO2_col(beg_col:end_col)   
  call hist_addfld1d(fname='ATM_CO2_col',units='umol/mol',avgflag='A',&
    long_name='Atmospheric CO2 concentration',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ATM_CH4_col(beg_col:end_col)   
  call hist_addfld1d(fname='ATM_CH4_col',units='umol/mol',avgflag='A',&
    long_name='Atmospheric CH4 concentration',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_NBP_col(beg_col:end_col)   
  call hist_addfld1d(fname='NBP_col',units='gC/m2',avgflag='I',&
    long_name='Cumulative net biosphere productivity (<0 into atmosphere)',&
    ptr_col=data1d_ptr,default='inactive')      
    
  data1d_ptr => this%h1D_ECO_LAI_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_LAI_col',units='m2/m2',avgflag='A',&
    long_name='Ecosystem LAI',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_SAI_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_SAI_col',units='m2/m2',avgflag='A',&
    long_name='Ecosystem stem Area index for all live branches',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Eco_GPP_CumYr_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_GPP_col',units='gC/m2',avgflag='I',&
    long_name='cumulative ecosystem GPP',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_RA_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_RA_col',units='gC/m2',avgflag='I',&
    long_name='cumulative ecosystem autotrophic respiration',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Eco_NPP_CumYr_col(beg_col:end_col)      
  call hist_addfld1d(fname='ECO_NPP_col',units='gC/m2',avgflag='I',&
    long_name='cumulative ecosystem NPP',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Eco_HR_CumYr_col(beg_col:end_col)      
  call hist_addfld1d(fname='ECO_RH_col',units='gC/m2',avgflag='I',&
    long_name='Cumulative ecosystem heterotrophic respiration (<0 into atmosphere)',&
    ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Eco_HR_CO2_col(beg_col:end_col)
  call hist_addfld1d(fname='ECO_HR_CO2_col',units='gC/m2/hr',avgflag='A',&
    long_name='Ecosystem heterotrophic respiration as CO2 (<0 into atmosphere)',&
    ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_Eco_HR_CO2_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='HR_CO2_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Heterotrophic respiration as CO2 in litter (<0 into atmosphere)',&
    ptr_col=data1d_ptr)      

!  data1d_ptr => this%h1D_Eco_HR_CH4_col(beg_col:end_col)
!  call hist_addfld1d(fname='ECO_RH_CH4',units='gC/m2/hr',avgflag='A',&
!    long_name='Ecosystem heterotrophic respiration as CH4',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tDIC_col(beg_col:end_col)       
  call hist_addfld1d(fname='tDIC_col',units='gC/m2',avgflag='A',&
    long_name='column integrated total soil DIC: CO2+CH4',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_tSTANDING_DEAD_C_col(beg_col:end_col)     
  call hist_addfld1d(fname='tSTANDING_DEAD_C_col',units='gC/m2',avgflag='A',&
    long_name='total standing dead C',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_tSTANDING_DEAD_N_col(beg_col:end_col)     
  call hist_addfld1d(fname='tSTANDING_DEAD_N_col',units='gN/m2',avgflag='A',&
    long_name='total standing dead N',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tSTANDING_DEAD_P_col(beg_col:end_col)     
  call hist_addfld1d(fname='tSTANDING_DEAD_P_col',units='gP/m2',avgflag='A',&
    long_name='total standing dead P',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_tPRECIP_col(beg_col:end_col)      
  call hist_addfld1d(fname='tPRECIP_col',units='mm/m2',avgflag='I',&
    long_name='cumulative precipitation, including irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_ET_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_ET_col',units='mm H2O/m2',avgflag='I',&
    long_name='cumulative total evapotranspiration',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_trcg_Ar_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_cumerr_col',units='gAr/m2',avgflag='I',&
    long_name='cumulative mass error for Ar',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_trcg_O2_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='O2_cumerr_col',units='gO/m2',avgflag='I',&
    long_name='cumulative mass error for O2',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_trcg_N2_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='N2_cumerr_col',units='gN/m2',avgflag='I',&
    long_name='cumulative mass error for N2',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_trcg_NH3_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='NH3_cumerr_col',units='gN/m2',avgflag='I',&
    long_name='cumulative mass error for NH3',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_trcg_H2_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='H2_cumerr_col',units='gH/m2',avgflag='I',&
    long_name='cumulative mass error for H2',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_trcg_CO2_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_cumerr_col',units='gC/m2',avgflag='I',&
    long_name='cumulative mass error for CO2',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_trcg_CH4_cumerr_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4_cumerr_col',units='gC/m2',avgflag='I',&
    long_name='cumulative mass error for CH4',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1d_CAN_NEE_col(beg_col:end_col)
  call hist_addfld1d(fname='CAN_NEE_col',units='umol C/m2/s',avgflag='A',&
    long_name='Canopy net CO2 exchange (<0 into atmosphere)',ptr_col=data1d_ptr)            

  data1d_ptr => this%h1D_ECO_RADSW_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_RADSW_col',units='W/m2',avgflag='A',&
    long_name='Shortwave radiation absorbed by the ecosystem',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_N2O_LITR_col(beg_col:end_col)      
  call hist_addfld1d(fname='N2O_LITR_col',units='g/m3',avgflag='A',&
    long_name='N2O solute concentration in soil micropores',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_NH3_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='NH3_LITR_col',units='g/m3',avgflag='A',&
    long_name='NH3 solute concentration in soil micropores',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_SOL_RADN_col(beg_col:end_col)      
  call hist_addfld1d(fname='SOL_RADN_col',units='W/m2',avgflag='A',&
    long_name='Incoming shortwave radiation on the ecosystem',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_AIR_TEMP_col(beg_col:end_col)      
  call hist_addfld1d(fname='AIR_TEMP_col',units='oC',avgflag='A',&
    long_name='air temperature',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_FreeNFix_col(beg_col:end_col)      
  call hist_addfld1d(fname='FreeNFix_col',units='gN/hr/m2',avgflag='A',&
    long_name='N fixation by free-living microbes',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_PATM_col(beg_col:end_col)      
  call hist_addfld1d(fname='PATM_col',units='kPa',avgflag='A',&
    long_name='atmospheric pressure',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HUM_col(beg_col:end_col)    
  call hist_addfld1d(fname='HUM_col',units='kPa',avgflag='A',&
    long_name='vapor pressure',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_WIND_col(beg_col:end_col)   
  call hist_addfld1d(fname='WIND_col',units='m/s',avgflag='A',&
    long_name='wind speed',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_PREC_col(beg_col:end_col)   
  call hist_addfld1d(fname='PREC_col',units='mm H2O/m2/hr',avgflag='A',&
    long_name='Total precipitation, excluding irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Snofall_col(beg_col:end_col)   
  call hist_addfld1d(fname='SNOFAL_col',units='mm H2O/m2/hr',avgflag='A',&
    long_name='Precipitation as snowfall',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SOIL_RN_col(beg_col:end_col)   
  call hist_addfld1d(fname='SOIL_RN_col',units='W/m2',avgflag='A',&
    long_name='total net radiation at ground surface (incoming short/long wave - outgoing short/long wave at soil/snow/litter)',&
    ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_LWSky_col(beg_col:end_col)      
  call hist_addfld1d(fname='LW_Sky_col',units='W/m2',avgflag='A',&
    long_name='Incoming sky long wave radiation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SOIL_LE_col(beg_col:end_col)      
  call hist_addfld1d(fname='SOIL_LE_col',units='W/m2',avgflag='A',&
    long_name='Latent heat flux into ground surface (exclude plant canopy)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SOIL_H_col(beg_col:end_col)      
  call hist_addfld1d(fname='SOIL_H_col',units='W/m2',avgflag='A',&
    long_name='Sensible heat flux into ground surface (exclude plant canopy)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SOIL_G_col(beg_col:end_col)  
  call hist_addfld1d(fname='SOIL_G_col',units='W/m2',avgflag='A',&
    long_name='total heat flux out of ground surface (>0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_RN_col(beg_col:end_col)     
  call hist_addfld1d(fname='ECO_Radnet_col',units='W/m2',avgflag='A',&
    long_name='Ecosystem net radiation (>0 into ecosystem, short+sky_long - plant_long-surf_long)',&
    ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_ECO_LE_col(beg_col:end_col)        
  call hist_addfld1d(fname='ECO_LE_col',units='W/m2',avgflag='A',&
    long_name='Ecosystem latent heat flux (>0 into surface)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Eco_HeatSen_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_HeatS_col',units='W/m2',avgflag='A',&
    long_name='Ecosystem sensible heat flux (>0 into surface)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_Heat2G_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_Heat2G_col',units='W/m2',avgflag='A',&
    long_name='Heat flux to warm the ecosystem (<0 into atmosphere),' &
    //' including canopy, snow, litter and exposed soil',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_O2_LITR_col(beg_col:end_col)      
  call hist_addfld1d(fname='O2w_conc_LITR_col',units='g/m3',avgflag='A',&
    long_name='O2 solute concentration in litter layer',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_MIN_LWP_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='MIN_LWP_pft',units='MPa',avgflag='A',&
    long_name='minimum daily canopy water potential',ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_SLA_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='SLA_pft',units='cm2 leaf (gC leaf)-1',avgflag='A',&
    long_name='Specific leaf area',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_CO2_SEMIS_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_SEMIS_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='Surface CO2 flux (< 0 into atmosphere), '// &
    'excluding wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_AR_SEMIS_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_SEMIS_FLX_col',units='umol Ar/m2/s',avgflag='A',&
    long_name='soil Ar flux (< 0 into atmosphere), '// &
    'excluding wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ECO_CO2_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_NEE_CO2_col',units='umol C/m2/s',avgflag='A',&
    long_name='ecosystem net CO2 exchange (<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CH4_SEMIS_FLX_col(beg_col:end_col)     
  call hist_addfld1d(fname='CH4_SEMIS_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='Surface CH4 flux (<0 into atmosphere), '// &
    'excluding wet deposition from rainfall and surface irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CH4_EBU_flx_col(beg_col:end_col)     
  call hist_addfld1d(fname='CH4_EBU_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='soil CH4 ebullition flux (<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Ar_EBU_flx_col(beg_col:end_col)     
  call hist_addfld1d(fname='Ar_EBU_FLX_col',units='umol Ar/m2/s',avgflag='A',&
    long_name='soil Ar ebullition flux (<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CO2_TPR_err_col(beg_col:end_col)
  call hist_addfld1d(fname='CumCO2_Transpt_Residual_col',units='gC/m2',avgflag='I',&
    long_name='Cumulative difference between soil CO2 production and surface CO2 flux',&
    ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CO2_Drain_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_DRAINLOSS_col',units='gC/m2/hr',avgflag='A',&
    long_name='CO2 loss flux through subsurface drainage',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CO2_hydloss_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_Cum_Hyd_Loss_col',units='gC/m2',avgflag='I',&
    long_name='Cumulative hydrological CO2 loss flux, including subsurface drainage',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_Ar_TPR_err_col(beg_col:end_col)
  call hist_addfld1d(fname='CumAr_Transpt_Residual_col',units='g/m2',avgflag='I',&
    long_name='Cumulative difference between soil Ar production and surface Ar flux',ptr_col=data1d_ptr,&
    default='inactive')      

  data1d_ptr => this%h1D_CH4_PLTROOT_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4_PLTROOT_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='soil CH4 flux through plants(<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_AR_PLTROOT_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_PLTROOT_FLX_col',units='umol Ar/m2/s',avgflag='A',&
    long_name='soil AR flux through plants(<0 into atmosphere)',ptr_col=data1d_ptr, &
    default='inactive')            

  data1d_ptr => this%h1D_CO2_PLTROOT_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_PLTROOT_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='soil CO2 flux through plants(<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_O2_PLTROOT_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='O2_PLTROOT_FLX_col',units='umol O2/m2/s',avgflag='A',&
    long_name='soil O2 flux through plants(<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CO2_DIF_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_DIF_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='soil CO2 flux through advection+diffusion (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_O2_DIF_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='O2_DIF_FLX_col',units='umol O2/m2/s',avgflag='A',&
    long_name='soil O2 flux through advection+diffusion (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CH4_DIF_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='CH4_DIF_FLX_col',units='umol C/m2/s',avgflag='A',&
    long_name='soil CH4 flux through advection+diffusion (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_NH3_DIF_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='NH3_DIF_FLX_col',units='umol N/m2/s',avgflag='A',&
    long_name='soil NH3 flux through advection+diffusion (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_Ar_DIF_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_DIF_FLX_col',units='umol Ar/m2/s',avgflag='A',&
    long_name='soil Ar flux through advection+diffusion (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_O2_SEMIS_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='O2_SEMIS_FLX_col',units='umol O2/m2/s',avgflag='A',&
    long_name='Surface O2 flux (<0 into atmosphere)',ptr_col=data1d_ptr)            

  data1d_ptr => this%h1D_CO2_LITR_col(beg_col:end_col)      
  call hist_addfld1d(fname='CO2_LITR_col',units='gC/m3',avgflag='A',&
    long_name='CO2 solute concentration in litter',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_EVAPN_col(beg_col:end_col)      
  call hist_addfld1d(fname='EVAPN_col',units='mm H2O/m2/hr',avgflag='A',&
    long_name='Column-integrated ground surface evaporation(<0 into atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CondGasXSurf_col(beg_col:end_col)      
  call hist_addfld1d(fname='GasXSurfConduct_col',units='m/hr',avgflag='A',&
    long_name='Conductance for soil-air gas exchange',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CANET_col(beg_col:end_col)
  call hist_addfld1d(fname='CANET_col',units='mm H2O/m2/hr',avgflag='A',&
    long_name='Column-integrated canopy evapotranspiration(<0 int atmosphere)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tSWC_col(beg_col:end_col)  
  call hist_addfld1d(fname='tSWC_col',units='mmH2O/m2',avgflag='A', &
    long_name='column integrated water content (include snow)',ptr_col=data1d_ptr)        

  data1d_ptr => this%h1D_tHeat_col(beg_col:end_col)
  call hist_addfld1d(fname='tSHeat_col',units='MJ/m2',avgflag='A', &
    long_name='column integrated heat content (include snow)',ptr_col=data1d_ptr,&
    default='inactive')              

  data1d_ptr => this%h1D_SNOWPACK_col(beg_col:end_col)      
  call hist_addfld1d(fname='SNOWPACK_col',units='mmH2O/m2',&
    avgflag='A',long_name='total water equivalent snow',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SNOWDENS_col(beg_col:end_col)      
  call hist_addfld1d(fname='SNOWDENS_col',units='kg/m3',&
    avgflag='A',long_name='snow density',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SURF_WTR_col(beg_col:end_col)   
  call hist_addfld1d(fname='SURF_WTR_col',units='m3/m3',avgflag='A',&
    long_name='Volumetric water content in surface litter layer',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ThetaW_litr_col(beg_col:end_col)   
  call hist_addfld1d(fname='ThetaW_litr_col',units='none',avgflag='A',&
    long_name='Relative saturation of water content in surface litter layer [0-1]',&
    ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_ThetaI_litr_col(beg_col:end_col)   
  call hist_addfld1d(fname='ThetaI_litr_col',units='none',avgflag='A',&
    long_name='Relative volume of ice content in surface litter layer [0-1]',&
    ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_SURF_ICE_col(beg_col:end_col)    
  call hist_addfld1d(fname='SURF_ICE_col',units='m3/m3',avgflag='A',&
    long_name='Volumetric ice content in surface litter layer',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_ACTV_LYR_col(beg_col:end_col)    
  call hist_addfld1d(fname='ACTV_LYR_col',units='m',avgflag='A',&
    long_name='active layer depth',ptr_col=data1d_ptr,default='inactive')              

  data1d_ptr => this%h1D_WTR_TBL_col(beg_col:end_col)      
  call hist_addfld1d(fname='WTR_TBL_col',units='m',avgflag='A',&
    long_name='internal water table depth (<0 below soil surface)',&
    ptr_col=data1d_ptr)            

  data1d_ptr => this%h1D_N2O_SEMIS_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='N2O_SEMIS_FLX_col',units='gN/m2/hr',&
    avgflag='A',long_name='Surface N2O flux (<0 into atmosphere), '// &
    'including wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_PAR_col(beg_col:end_col)
  call hist_addfld1d(fname='PAR_col',units='umol m-2 s-1',avgflag='A',&
    long_name='Direct plus diffusive incoming photosynthetic photon flux density',&
    ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1d_fPAR_col(beg_col:end_col)
  call hist_addfld1d(fname='fPAR_col',units='-',avgflag='P',&
    long_name='Fraction of absorbed PAR by canopy',&
    ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_CO2_WetDep_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='CO2_WetDep_FLX_col',units='gC/m2/hr',&
    avgflag='A',long_name='Wet deposition CO2 flux to soil, '// &
    'from rainfall and irrigation (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_RootN_Fix_col(beg_col:end_col)
  call hist_addfld1d(fname='Root_N_FIX_col',units='gN/m2/hr',&
    avgflag='A',long_name='Root N2 fixation',ptr_col=data1d_ptr)            

  data1d_ptr => this%h1D_AR_WetDep_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='Ar_WetDep_FLX_col',units='gAr/m2/hr',&
    avgflag='A',long_name='Wet deposition Ar flux to soil, '// &
    'from rainfall and irrigation (<0 into atmosphere)',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_RootXO2_flx_col(beg_col:end_col)
  call hist_addfld1d(fname='RootO2_X_Flx_col',units='gO2/m2/hr',&
    avgflag='A',long_name='O2 consumption rates in roots',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_N2_SEMIS_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='N2_SEMIS_FLX_col',units='gN/m2/hr',&
    avgflag='A',long_name='Surface N2 flux (<0 into atmosphere), '// &
    'including wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_NH3_SEMIS_FLX_col(beg_col:end_col)       
  call hist_addfld1d(fname='NH3_SEMIS_FLX_col',units='gN/m2/hr',avgflag='A',&
    long_name='Surface NH3 flux (<0 into atmosphere), '// &
    'including wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)            

  data1d_ptr => this%h1D_H2_SEMIS_FLX_col(beg_col:end_col)       
  call hist_addfld1d(fname='H2_SEMIS_FLX_col',units='gH/m2/hr',avgflag='A',&
    long_name='Surface H2 flux (<0 into atmosphere), '// &
    'including wet deposition from rainfall and irrigation',ptr_col=data1d_ptr)        

  data1d_ptr => this%h1D_VHeatCap_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='vHeatCap_litr_col',units='MJ/m3/K',avgflag='A',&
    long_name='surface litter heat capacity',ptr_col=data1d_ptr, &
    default='inactive')            

  data1d_ptr => this%h1D_RUNOFF_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='RUNOFF_FLX_col',units='mmH2O/m2/hr',avgflag='A',&
    long_name='Surface runoff from surface water',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_SEDIMENT_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='SEDIMENT_FLX_col',units='kg/m2/hr',avgflag='A',&
    long_name='total sediment subsurface flux',ptr_col=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_QDISCHG_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='QDischarge_FLX_col',units='mmH2O/m2/hr',avgflag='A',&
    long_name='grid water lateral discharge with respect external water table (>0 out of grid)',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_HeatDISCHG_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='HeatDischarge_FLX_col',units='MJ/m2/hr',avgflag='A',&
    long_name='Column-integrated heat flux through discharge',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_LEAF_PC_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='LEAF_rPC_pft',units='gP/gC',avgflag='I',&
    long_name='Mass based leaf PC ratio',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CAN_RN_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CAN_RN_pft',units='W/m2',avgflag='A',&
    long_name='Canopy net radiation',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CAN_LE_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CAN_LE_pft',units='W/m2',avgflag='A',&
    long_name='Canopy latent heat flux (<0 to ATM)',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CAN_H_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='CAN_H_pft',units='W/m2',avgflag='A',&
    long_name='Canopy sensible heat flux',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CAN_G_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='CAN_G_pft',units='W/m2',avgflag='A',&
    long_name='Canopy storage heat flux',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CAN_TEMPC_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_TEMPC_pft',units='oC',avgflag='A',&
    long_name='Canopy temperature',ptr_patch=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_CAN_TEMPFN_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CANGRO_TEMP_FN_pft',units='none',avgflag='A',&
    long_name='Canopy temperature growth function/stress',ptr_patch=data1d_ptr,&
    default='inactive')              

  data1d_ptr => this%h1D_CAN_CO2_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_CO2_FLX_pft',units='umol C/m2/s',avgflag='A',&
    long_name='Canopy net CO2 exchange',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CAN_GPP_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_GPP_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant canopy gross CO2 fixation',ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_CAN_cumGPP_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_cumGPP_pft',units='gC/m2',avgflag='I',&
    long_name='Plant canopy cumulative gross CO2 fixation',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_dCAN_GPP_CLIM_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='dCAN_GPP_CLIM_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant canopy CO2-limited gross CO2 fixation minus acutal value (>0 light-limitation)',ptr_patch=data1d_ptr,&
    default='inactive')      

  data1d_ptr => this%h1D_dCAN_GPP_eLIM_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='dCAN_GPP_eLIM_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant canopy light-limited gross CO2 fixation minus actual value (>0 C-limitation)',ptr_patch=data1d_ptr,&
    default='inactive')      

  data1d_ptr => this%h1D_CAN_RA_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_RA_pft',units='gC/m2/hr',avgflag='A',&
    long_name='total aboveground autotrophic respiration',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CAN_GROWTH_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_GROWTH_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Canopy structural growth rate',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_cTNC_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='cAbvNonstC_pft',units='gC/gC',avgflag='A',&
    long_name='Canopy nonstructural C concentration',ptr_patch=data1d_ptr,&
    default='inactive')      

  data1d_ptr => this%h1D_cTNN_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='cAbvNonstN_pft',units='gN/gC',avgflag='A',&
    long_name='Canopy nonstructural N concentration',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_cTNP_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='cAbvNonstP_pft',units='gP/gC',avgflag='A',&
    long_name='Canopy nonstructural P concentration',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CanNonstBConc_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='CanNonstBConc_pft',units='g',avgflag='A',&
    long_name='Canopy nonstructural biomass concentration',ptr_patch=data1d_ptr,&
    default='inactive')          

  data1d_ptr => this%h1D_STOML_RSC_CO2_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='STOML_RSC_CO2_pft',units='s/m',avgflag='A',&
    long_name='Canopy stomatal resistance for CO2',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STOML_Min_RSC_CO2_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='STOML_MinRSC_CO2_pft',units='s/m',avgflag='A',&
    long_name='Canopy minimal stomatal resistance for CO2',ptr_patch=data1d_ptr,&
    default='inactive')      

  data1d_ptr => this%h1D_Km_CO2_carboxy_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Km_CO2_carboxy_pft',units='uM',avgflag='A',&
    long_name='MM parameter for CO2 carboxylation by Rubisco',ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_Ci_mesophyll_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Ci_mesophyll_pft',units='uM',avgflag='A',&
    long_name='Intracellular CO2 concentration for photosynthesis',ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_BLYR_RSC_CO2_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='BLYR_RSC_CO2_pft',units='s/m',avgflag='A',&
    long_name='Canopy boundary layer resistance for CO2',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CAN_CO2_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CAN_CO2_pft',units='umol/mol',avgflag='A',&
    long_name='Canopy gaesous CO2 concentration',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_O2L_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='Leaf_O2_pft',units='umol/mol',avgflag='A',&
    long_name='Leaf aqueous O2 concentration',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_LAI_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='LAIstk_pft',units='m2/m2',avgflag='A',&
    long_name='whole plant leaf area, including stalk',ptr_patch=data1d_ptr)            

  data1d_ptr => this%h1D_PSI_CAN_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='PSI_CAN_pft',units='MPa',avgflag='A',&
    long_name='Canopy total water potential',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_RootAR_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='RootAR_pft',units='gC/m2/h',avgflag='A',&
    long_name='Root autotrophic respiraiton',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_RootLenPerPlant_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='RootLen_pft',units='m plant-1',avgflag='A',&
    long_name='Root length per pft (excluding root hair)',ptr_patch=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_TURG_CAN_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='TURG_CAN_pft',units='MPa',avgflag='A',&
    long_name='Canopy turgor water potential',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_STOML_RSC_H2O_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='STOM_RSC_H2O_pft',units='s/m',avgflag='A',&
    long_name='Canopy stomatal resistance for H2O',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_BLYR_RSC_H2O_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='BLYR_RSC_H2O_pft',units='s/m',avgflag='A',&
    long_name='Canopy boundary layer resistance for H2O',ptr_patch=data1d_ptr,&
    default='inactive')            
  
  data1d_ptr => this%h1D_TRANSPN_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='QvTransp_pft',units='mmH2O/m2/h',avgflag='A',&
    long_name='Canopy transpiration (<0 into atmosphere)',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_NH4_UPTK_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='UPTK_NH4_FLX_pft',units='gN/m2/hr',&
    avgflag='A',long_name='total root uptake of NH4',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_NO3_UPTK_FLX_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='UPTK_NO3_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='total root uptake of NO3',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_N2_FIXN_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='N2_FIX_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='total root N2 fixation',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_cNH3_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='CAN_NH3_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='*canopy NH3 flux',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_TC_Groth_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TC_Groth_pft',units='oC',avgflag='A',&
    long_name='Plant growth temperature',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_PO4_UPTK_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='UPTK_PO4_FLX_pft',units='gP/m2/hr',avgflag='A',&
    long_name='total root uptake of PO4',ptr_patch=data1d_ptr,&
    default='inactive')                  

  data1d_ptr => this%h1D_SHOOT_C_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='SHOOT_C_pft',units='gC/m2',avgflag='A',&
    long_name='Live plant shoot C',ptr_patch=data1d_ptr,&
    default='inactive')                  

  data1d_ptr => this%h1D_frcPARabs_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='frcPARabs_pft',units='none',avgflag='A',&
    long_name='fraction of PAR absorbed by plant',ptr_patch=data1d_ptr,&
    default='inactive')                  

  data1d_ptr => this%h1D_PAR_CAN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Canopy_PAR_pft',units='umol m-2 s-1',avgflag='A',&
    long_name='PAR absorbed by plant canopy',ptr_patch=data1d_ptr,&
    default='inactive')                  

  data1d_ptr => this%h1D_Plant_C_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='Plant_C_pft',units='gC/m2',avgflag='A',&
    long_name='Plant C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_LEAF_C_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='LEAF_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy leaf C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_Petole_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='PETIOLE_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy sheath C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_STALK_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='STALK_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy stalk C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_RESERVE_C_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='RESERVE_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy reserve C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_HUSK_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='HUSK_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy husk C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_GRAIN_C_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='GRAIN_C_pft',units='gC/m2',avgflag='A',&
    long_name='Canopy grain C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_ROOT_C_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='Root_C_pft',units='gC/m2',avgflag='A',&
    long_name='Plant root C, exluding nodule',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_ROOTST_C_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootST_C_pft',units='gC/m2',avgflag='A',&
    long_name='Plant root structural C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_ROOTST_N_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootST_N_pft',units='gN/m2',avgflag='A',&
    long_name='Plant root structural N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_ROOTST_P_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='RootST_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant root structural P',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_RootNodule_C_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='RootNodule_C_pft',units='gC/m2',avgflag='A',&
    long_name='Root total nodule C',ptr_patch=data1d_ptr,default='inactive')     

  data1d_ptr => this%h1D_ShootNodule_C_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='ShootNodule_C_pft',units='gC/m2',avgflag='A',&
    long_name='Shoot total nodule C',ptr_patch=data1d_ptr,default='inactive')    

  data1d_ptr => this%h1D_ShootNodule_N_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='ShootNodule_N_pft',units='gN/m2',avgflag='A',&
    long_name='Shoot total nodule N',ptr_patch=data1d_ptr,default='inactive')    

  data1d_ptr => this%h1D_ShootNodule_P_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='ShootNodule_P_pft',units='gP/m2',avgflag='A',&
    long_name='Shoot total nodule P',ptr_patch=data1d_ptr,default='inactive')    

  data1d_ptr => this%h1D_STORED_C_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='SSTORED_C_pft',units='gC/m2',avgflag='A',&
    long_name='Plant seasonal storage of nonstructural C',ptr_patch=data1d_ptr) 

  data1d_ptr => this%h1D_ROOT_NONSTC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root_NONSTC_pft',units='gC/m2',avgflag='A',&
    long_name='Plant root storage of nonstructural C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_ROOT_NONSTN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root_NONSTN_pft',units='gN/m2',avgflag='A',&
    long_name='Plant root storage of nonstructural N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_ROOT_NONSTP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Root_NONSTP_pft',units='gP/m2',avgflag='A',&
    long_name='Plant root storage of nonstructural P',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_SHOOT_NONSTC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_NONSTC_pft',units='gC/m2',avgflag='A',&
    long_name='Plant leaf storage of nonstructural C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_SHOOT_NONSTN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_NONSTN_pft',units='gN/m2',avgflag='A',&
    long_name='Plant leaf storage of nonstructural N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_SHOOT_NONSTP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='SHOOT_NONSTP_pft',units='gP/m2',avgflag='A',&
    long_name='Plant leaf storage of nonstructural P',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_LeafC3ChlCperm2LA_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafC3ChlC_pft',units='mgC chl m-2 leaf area',avgflag='A',&
    long_name='Cholorophyll carbon in mesophyll for C3/bundle sheath for C4 plants',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_LeafC4ChlCperm2LA_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafC4ChlC_pft',units='mgC chl m-2 leaf area',avgflag='A',&
    long_name='chlorophyll carbon in mesophyll for C4 plants',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_LeafRubiscoCperm2LA_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafRubisco_pft',units='gC rubisco m-2 leaf area',avgflag='A',&
    long_name='Rubisco carbon in mesophyll for C3/bundle sheath for C4 plants',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_LeafPEPCperm2LA_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafPEPC_pft',units='gC PEP m-2 leaf area',avgflag='A',&
    long_name='PEP carbon in mesophyll for C4 plants',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_GRAIN_NO_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='GRAIN_NO_pft',units='1/m2',avgflag='A',&
    long_name='Canopy grain number',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_LAIb_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='LAI_xstk_pft',units='m2/m2',avgflag='A',&
    long_name='whole plant leaf area, exclude stalk',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_EXUD_CumYr_C_FLX_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='EXUD_CumYr_C_FLX_pft',units='gC/m2',avgflag='I',&
    long_name='Cumulative root organic C uptake (<0 exudation into soil)',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_LITRf_C_FLX_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='LITRf_C_pft',units='gC/m2/hr',avgflag='A',&
    long_name='total plant LitrFall C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_SURF_LITRf_C_FLX_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='SURF_LITRf_C_FLX_pft',units='gC/m2',avgflag='I',&
    long_name='Cumulative plant LitrFall C to the soil surface',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_AUTO_RESP_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='AUTO_RESP_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Whole plant autotrophic respiration',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_HVST_C_FLX_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='HVST_C_FLX_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant C harvest',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_PLANT_BALANCE_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='Plant_BALANCE_C_pft',units='gC/m2',avgflag='A',&
    long_name='Cumulative plant C conservation error',ptr_patch=data1d_ptr)                  

  data1d_ptr => this%h1D_STANDING_DEAD_C_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='STANDING_DEAD_C_pft',units='gC/m2',avgflag='A',&
    long_name='pft Standing dead C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_FIREp_CO2_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='FIREp_CO2_FLX_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant CO2 from fire',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_FIREp_CH4_FLX_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='FIREp_CH4_FLX_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Plant CH4 emission from fire',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_cNPP_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='cNPP_pft',units='gC/m2',avgflag='I',&
    long_name='Plant cumulative net primary productivity',ptr_patch=data1d_ptr)                  

  data1d_ptr => this%h1D_CAN_HT_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_HT_pft',units='m',avgflag='A',&
    long_name='Canopy height',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_POPN_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='POPN_pft',units='1/m2',avgflag='A',&
    long_name='Plant population',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_tTRANSPN_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='tTRANSPN_pft',units='mmH2O/m2',avgflag='I',&
    long_name='Cumulative canopy evapotranspiration (>0 into atmosphere)',ptr_patch=data1d_ptr,&
    default='inactive')                  

  data1d_ptr => this%h1D_WTR_STRESS_ptc(beg_ptc:end_ptc)    !HoursTooLowPsiCan_pft(NZ,NY,NX)
  call hist_addfld1d(fname='WTR_STRESS_pft',units='hr',avgflag='A',&
    long_name='Canopy plant water stress indicator: number of ' &
    //'hours PSICanopy_pft(< PSILY)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_OXY_STRESS_ptc(beg_ptc:end_ptc)    !OSTR(NZ,NY,NX)
  call hist_addfld1d(fname='OXY_STRESS_pft',units='none',avgflag='A',&
    long_name='Plant root O2 stress indicator [0->1 weaker stress]',ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_VcMaxRubisco_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='VcMax25C_RUBISCO_pft',units='umol CO2 s-1 m-2 leaf area',avgflag='A',&
    long_name='Maximum carboxylation rate by Rubisco at 25oC',&
    ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_LeafProteinCperm2_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafProteinC_pft',units='gC protein m-2 leaf area',avgflag='A',&
    long_name='Protein carbon mass per unit of leaf area',&
    ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_VoMaxRubisco_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='VoMax25C_RUBISCO_pft',units='umol O2 s-1 m-2 leaf area',avgflag='A',&
    long_name='Maximum oxygenation rate by Rubisco at 25oC',&
    ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_VcMaxPEP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='VcMax25C_PEP_pft',units='umol CO2 s-1 m-2 leaf area',avgflag='A',&
    long_name='Maximum carboxylation rate by PEP at 25oC',&
    ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_JMaxPhoto_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='JMax25C_photo_pft',units='umol e- s-1 m-2 leaf area',avgflag='A',&
    long_name='Maximum electron transport rate at 25oC',&
    ptr_patch=data1d_ptr)      

  data1d_ptr => this%h1D_TFN_Carboxy_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TFN_Carboxy_pft',units='none',avgflag='A',&
    long_name='Temperature response of carboyxlation in photosynthesis',&
    ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_TFN_Oxygen_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TFN_Oxygen_pft',units='none',avgflag='A',&
    long_name='Temperature response of oxygenation in photosynthesis',&
    ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_TFN_eTranspt_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='TFN_eTranspt_pft',units='none',avgflag='A',&
    long_name='Temperature response of electron transport in photosynthesis',&
    ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_PARSunlit_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='PARSunlit_pft',units='umol m-2 s-1',avgflag='A',&
    long_name='PAR absorbed by sunlit leaves',&
    ptr_patch=data1d_ptr,default='inactive')      
  
  data1d_ptr => this%h1D_PARSunsha_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='PARSunsha_pft',units='umol m-2 s-1',avgflag='A',&
    long_name='PAR absorbed by sun-shaded leaves',&
    ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_CH2OSunlit_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CH2OSunlit_pft',units='gC m-2 h-1',avgflag='A',&
    long_name='Photosynthesis by sunlit leaves',&
    ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_CH2OSunsha_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CH2OSunsha_pft',units='gC m-2 h-1',avgflag='A',&
    long_name='Photosynthesis by sun-shaded leaves',&
    ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_LeafAreaSunlit_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='LeafArea_sunlit_pft',units='m2 m-2',avgflag='A',&
    long_name='Irridiance-lit leaf area',ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_fClump_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='fClump_pft',units='none',avgflag='A',&
    long_name='Clumping factor of leaf area',ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_SHOOT_N_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='SHOOT_N_pft',units='gN/m2',avgflag='A',&
    long_name='Live plant shoot N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_Plant_N_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='Plant_N_pft',units='gN/m2',avgflag='A',&
    long_name='Plant N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_LEAF_N_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='LEAF_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy leaf N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_fCNLFW_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='fLEAF_CN_pft',units='gN/gC',avgflag='A',&
    long_name='Canopy new leaf C:N mass ratio',ptr_patch=data1d_ptr) 

  data1d_ptr => this%h1D_fCPLFW_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='fLEAF_CP_pft',units='gP/gC',avgflag='A',&
    long_name='Canopy new leaf C:P mass ratio',ptr_patch=data1d_ptr) 

  data1d_ptr => this%h1D_LEAFN2LAI_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='LEAFN2LAI_pft',units='gN m-2 LA',avgflag='A',&
    long_name='Canopy leaf N per m2 leaf area',ptr_patch=data1d_ptr)

  data1d_ptr => this%h1D_Petole_N_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='Petiole_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy sheath N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_STALK_N_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='STALK_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy stalk N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_RESERVE_N_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='RESERVE_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy reserve N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_HUSK_N_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='HUSK_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy husk N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_GRAIN_N_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='GRAIN_N_pft',units='gN/m2',avgflag='A',&
    long_name='Canopy grain C',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_ROOT_N_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='Root_N_pft',units='gN/m2',avgflag='A',&
    long_name='Root nitrogen',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_RootNodule_N_ptc(beg_ptc:end_ptc)        
  call hist_addfld1d(fname='RootNodule_N_pft',units='gN/m2',avgflag='A',&
    long_name='Root total nodule N',ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_STORED_N_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='SSTORED_N_pft',units='gN/m2',avgflag='A',&
    long_name='Plant seasonal storage of nonstructural N',ptr_patch=data1d_ptr,&
    default='inactive')                  

  data1d_ptr => this%h1D_EXUD_N_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='EXUD_CumYr_N_FLX_pft',units='gN/m2',avgflag='I',&
    long_name='Cumulative Root organic N uptake (<0 exudation into soil)',&
    ptr_patch=data1d_ptr,default='inactive')                  

  data1d_ptr => this%h1D_Uptk_N_Flx_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Uptk_N_CumYr_FLX_pft',units='gN/m2',avgflag='I',&
    long_name='Cumulative Root N uptake (including <0 exudation to soil)',ptr_patch=data1d_ptr,&
    default='inactive')      

  data1d_ptr => this%h1D_Uptk_P_Flx_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Uptk_P_CumYr_FLX_pft',units='gP/m2',avgflag='I',&
    long_name='Cumulative Root P uptake (including <0 exudation to soil)',ptr_patch=data1d_ptr,&
    default='inactive')      

  data1d_ptr => this%h1D_LITRf_N_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='LITRf_N_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='total plant LitrFall N',ptr_patch=data1d_ptr,&
    default='inactive')      

  data1d_ptr => this%h1D_cum_N_FIXED_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='cumN_FIXED_pft',units='gN/m2',avgflag='I',&
    long_name='cumulative plant N2 fixation',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_TreeRingRadius_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='TreeRingRadius_pft',units='m',avgflag='I',&
    long_name='Mean tree ring radius',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_HVST_N_FLX_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='HVST_N_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='Plant N harvest',ptr_patch=data1d_ptr,&
    default='inactive')      

  data1d_ptr => this%h1D_NH3can_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='NH3can_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='total canopy NH3 flux',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_PLANT_BALANCE_N_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='Plant_BALANCE_N_pft',units='gC/m2',avgflag='A',&
    long_name='Cumulative plant N conservation error',ptr_patch=data1d_ptr)            

  data1d_ptr => this%h1D_STANDING_DEAD_N_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='STANDING_DEAD_N_pft',units='gN/m2',avgflag='A',&
    long_name='pft standing dead N',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_FIREp_N_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='FIREp_N_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='Plant N emission from fire',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_SURF_LITRf_N_FLX_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='SURF_LITRf_N_FLX_pft',units='gN/m2/hr',avgflag='A',&
    long_name='total surface LitrFall N',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_SHOOT_P_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='SHOOT_P_pft',units='gP/m2',avgflag='A',&
    long_name='Live plant shoot P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_Plant_P_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='Plant_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_stomatal_stress_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='STOMATAL_STRESS_pft',units='none',avgflag='A',&
    long_name='stomatal stress from root turogr [0->1 increasing stress]',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_CANDew_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='Canopy_DEW_pft',units='mm H2O/m2',avgflag='I',&
    long_name='Cumulative canopy dew deposition',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_LEAF_P_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='LEAF_P_pft',units='gP/m2',avgflag='A',&
    long_name='Canopy leaf P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_Petole_P_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='Petiole_P_pft',units='gP/m2',avgflag='A',&
    long_name='Canopy sheath P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_STALK_P_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='STALK_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant stalk P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_RESERVE_P_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='RESERVE_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant reserve P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_HUSK_P_ptc(beg_ptc:end_ptc)        
  call hist_addfld1d(fname='HUSK_P_pft',units='gP/m2',avgflag='A',&
    long_name='Husk P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_GRAIN_P_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='GRAIN_P_pft',units='gP/m2',avgflag='A',&
    long_name='Canopy grain P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_ROOT_P_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='Root_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant root P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_RootNodule_P_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='RootNodule_P_pft',units='gP/m2',avgflag='A',&
    long_name='Root total nodule P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_STORED_P_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='SSTORED_P_pft',units='gP/m2',avgflag='A',&
    long_name='Plant seasonal storage of nonstructural P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_EXUD_P_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='EXUD_CumYr_P_FLX_pft',units='gP/m2',avgflag='I',&
    long_name='Cumulative root organic P uptake (<0 exudation into soil)',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_RDECOMPC_SOM_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RDecompC_SOM_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Hydrolysis of SOM C in litter layer',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_MicrobAct_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='MicrobAct_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Respiration-based micoribal activity for hydrolysis in litter layer',ptr_col=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_RDECOMPC_BReSOM_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RDecompC_BReSOM_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Hydrolysis of microbial residual OM C in litter layer',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_RDECOMPC_SorpSOM_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RDecompC_SorpSOM_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Hydrolysis of adsorbed OM C in litter layer',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_tRespGrossHeteUlm_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='tRespGrossHeteUlm_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Oxygen unlimited gross heterotrophic respiraiton in litter layer',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_tRespGrossHete_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='tRespGrossHete_litr_col',units='gC/m2/hr',avgflag='A',&
    long_name='Oxygen-limited gross heterotrophic respiraiton in litter layer',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_LITRf_P_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='LITRf_P_FLX_pft',units='gP/m2/hr',avgflag='A',&
    long_name='total plant LitrFall P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_HVST_P_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='HVST_P_FLX_pft',units='gP/m2/hr',avgflag='A',&
    long_name='Plant P harvest',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_PLANT_BALANCE_P_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='Plant_BALANCE_P_pft',units='gP/m2',avgflag='A',&
    long_name='Cumulative plant P conservation error',ptr_patch=data1d_ptr)            

  data1d_ptr => this%h1D_STANDING_DEAD_P_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='STANDING_DEAD_P_pft',units='gP/m2',avgflag='A',&
    long_name='pft Standing dead P',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_FIREp_P_FLX_ptc(beg_ptc:end_ptc)             
  call hist_addfld1d(fname='FIREp_P_FLX_pft',units='gP/m2/hr',avgflag='A',&
    long_name='Plant PO4 emission from fire',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_SURF_LITRf_P_FLX_ptc(beg_ptc:end_ptc)         !SurfLitrfallElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
  call hist_addfld1d(fname='SURF_LITRf_P_FLX_pft',units='gP/m2/hr',avgflag='A',&
    long_name='Plant LitrFall P to the soil surface',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_ShootRootXferC_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ShootRoot_XFER_C_pft',units='gC/m2/hr',avgflag='A',&
    long_name='Shoot C transfered to root',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_ShootRootXferN_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ShootRoot_XFER_N_pft',units='gN/m2/hr',avgflag='A',&
    long_name='Shoot N transfered to root',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_ShootRootXferP_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='ShootRoot_XFER_P_pft',units='gP/m2/hr',avgflag='A',&
    long_name='Shoot P transfered to root',ptr_patch=data1d_ptr,&
    default='inactive')            

  data1d_ptr => this%h1D_BRANCH_NO_ptc(beg_ptc:end_ptc)            !NumOfBranches_pft(NZ,NY,NX)
  call hist_addfld1d(fname='BRANCH_NO_pft',units='none',avgflag='I',&
    long_name='Plant branch number',ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_Growth_Stage_ptc(beg_ptc:end_ptc)  !plant development stage, integer, 0-10, planting, emergence, floral_init, jointing, 
                                                               !elongation, heading, anthesis, seed_fill, see_no_set, seed_mass_set, end_seed_fill
  call hist_addfld1d(fname='Growth_Stage_pft',units='none',avgflag='I',&
    long_name='Plant development stage, integer, 0-planting, 1-emergence, 2-floral_init, 3-jointing,'// &
    '4-elongation, 5-heading, 6-anthesis, 7-seed_fill, 8-see_no_set, 9-seed_mass_set, 10-end_seed_fill',&
    ptr_patch=data1d_ptr,default='inactive')            

  data1d_ptr => this%h1D_LEAF_NC_ptc(beg_ptc:end_ptc)            
  call hist_addfld1d(fname='LEAF_rNC_pft',units='gN/gC',avgflag='A',&
    long_name='Mass based plant leaf NC ratio',ptr_patch=data1d_ptr,default='inactive')       

    data1d_ptr => this%h1D_MainBranchNO_ptc(beg_ptc:end_ptc)            
  call hist_addfld1d(fname='MainBranchNO_pft',units='-',avgflag='A',&
    long_name='Main branch number',ptr_patch=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_RCanMaintDef_CO2_pft(beg_ptc:end_ptc)            
  call hist_addfld1d(fname='RCanMaintDef_CO2_pft',units='gC m-2 h-1',avgflag='A',&
    long_name='Canopy maintenance respiraiton deficit as CO2 (<0 deficit)',ptr_patch=data1d_ptr,default='inactive')       

  data1d_ptr => this%h1D_RootMaintDef_CO2_pft(beg_ptc:end_ptc)            
  call hist_addfld1d(fname='RootMaintDef_CO2_pft',units='gC m-2 h-1',avgflag='A',&
    long_name='Root maintenance respiraiton deficit as CO2 (<0 deficit)',ptr_patch=data1d_ptr,default='inactive')       

  data2d_ptr => this%h2D_QDrainloss_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='QDrainloss_vr',units='mm H2O/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved water drainage',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_FermOXYI_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermOXYI_vr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved oxygen inhibitor of fermentation [0->1, weaker inhibition]',&
    ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_litrC_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='litrC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Column-level Vertically resolved litter C',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_litrN_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='litrN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved litter N',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_litrP_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='litrP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved litter P',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_tSOC_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tSOC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total soil organic C (everything organic)',&
    ptr_col=data2d_ptr)       

  data2d_ptr => this%h2D_POM_C_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='POM_C_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved particulate organic C (resulting from the humification process)',&
    ptr_col=data2d_ptr)       

  data2d_ptr => this%h2D_MAOM_C_vr(beg_col:end_col,1:JZ)      
  call hist_addfld2d(fname='MAOM_C_vr',units='gC (kg soil)-1',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved mineral associated organic C (resulting from sorption process)',&
    ptr_col=data2d_ptr)       

  data2d_ptr => this%h2D_microbC_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tMicrobeC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total live microbial C',&
    ptr_col=data2d_ptr)       

  data2d_ptr => this%h2D_microbN_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tMicrobeN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total live microbial N',&
    ptr_col=data2d_ptr)       

  data2d_ptr => this%h2D_microbP_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tMicrobeP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total live microbial P',&
    ptr_col=data2d_ptr)       

  data2d_ptr => this%h2D_AeroBact_PrimS_lim_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='AeroBact_PrimS_lim_vr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved primary substrate limitation for aerobic heterotrophic bacteria',&
    ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_AeroFung_PrimS_lim_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='AeroFung_PrimS_lim_vr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved primary substrate limitation for aerobic heterotrophic fungi',&
    ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_tSOCL_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tSOCL_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved total soil organic C (everything organic)',ptr_col=data2d_ptr,&
    default='inactive')       

  data2d_ptr => this%h2D_tSON_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tSON_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total soil organic N (everything organic)',ptr_col=data2d_ptr,&
    default='inactive')       

  data2d_ptr => this%h2D_BotDEPZ_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='BOTDEPZ_vr',units='m',type2d='levsoi',avgflag='A',&
    long_name='Bottom depth of soil layer',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_tSOP_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tSOP_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total soil organic P (everything organic)',&
    ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_NO3_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='NO3_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved dissolved NO3 concentration',&
    ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_NH4_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='NH4_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved dissolved NH4 concentration',&
    ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_VHeatCap_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='VHeatCap_vr',units='MJ/m3/K',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved Volumetric heat capacity',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_LEAF_NODE_NO_ptc(beg_ptc:end_ptc,1:MaxNumBranches)        !NumOfLeaves_brch(MainBranchNum_pft(NZ,NY,NX),NZ,NY,NX), leaf NO
  call hist_addfld2d(fname='LEAF_NODE_NO_pft',units='none',type2d='nbranches',avgflag='I',&
    long_name='Leaf number',ptr_patch=data2d_ptr,default='inactive')       

  data1d_ptr => this%h1D_RUB_ACTVN_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='RUB_ACTVN_pft',units='none',avgflag='A',&
    long_name='mean rubisco activity for CO2 fixation across branches, 0-1',ptr_patch=data1d_ptr,default='inactive')       

  data1d_ptr => this%h1D_CanopyNLim_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CanopyNLim_pft',units='none',avgflag='A',&
    long_name='mean canopy nitrogen limitation across branches, 0->1 weaker limitation',ptr_patch=data1d_ptr,default='inactive')       

  data1d_ptr => this%h1D_CanopyPLim_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CanopyPLim_pft',units='none',avgflag='A',&
    long_name='mean canopy phosphorus limitation across branches, 0->1 weaker limitation',ptr_patch=data1d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RNITRIF_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RNITRIF_vr',units='gN/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Nitrification rate from NH3 oxidation in each soil layer',ptr_col=data2d_ptr,&
    default='inactive')       

  data2d_ptr => this%h2D_Root_CO2_vr(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_CO2,1:JZ,NY,NX)
  call hist_addfld2d(fname='Root_CO2_mass_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved CO2 mass in roots',ptr_col=data2d_ptr,default='inactive')      

  data2d_ptr => this%h2D_O2_rootconduct_pvr(beg_ptc:end_ptc,1:JZ)    
  call hist_addfld2d(fname='O2_root_conductance_pvr',units='1/h',type2d='levsoi',avgflag='A',&
    long_name='Root conductance for O2 gaseous in soil layer',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_CO2_rootconduct_pvr(beg_ptc:end_ptc,1:JZ)    
  call hist_addfld2d(fname='CO2_root_conductance_pvr',units='1/h',type2d='levsoi',avgflag='A',&
    long_name='Root conductance for CO2 gaseous in soil layer',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_ProteinNperm2LeafArea_pnd(beg_ptc:end_ptc,1:JZ)    
  call hist_addfld2d(fname='ProteinNperm2LeafArea_pnd',units='gC/(m2 leaf area)',type2d='node',avgflag='A',&
    long_name='Areal leaf protein N concentration by node',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Aqua_CO2_vr(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_CO2,1:JZ,NY,NX)
  call hist_addfld2d(fname='CO2w_conc_vr',units='gC/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous CO2 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Aqua_CH4_vr(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='CH4w_conc_vr',units='gC/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous CH4 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Aqua_O2_vr(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='O2w_conc_vr',units='g/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous O2 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Aqua_N2O_vr(beg_col:end_col,1:JZ)    
  call hist_addfld2d(fname='N2Ow_conc_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous N2O concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Aqua_NH3_vr(beg_col:end_col,1:JZ)    
  call hist_addfld2d(fname='NH3w_conc_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous NH3 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Aqua_H2_vr(beg_col:end_col,1:JZ)    
  call hist_addfld2d(fname='H2w_conc_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous H2 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Aqua_Ar_vr(beg_col:end_col,1:JZ)    
  call hist_addfld2d(fname='Arw_conc_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous Ar concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Aqua_N2_vr(beg_col:end_col,1:JZ)    
  call hist_addfld2d(fname='N2w_conc_vr',units='gN/m3 water',type2d='levsoi',avgflag='A',&
    long_name='Aqueous N2 concentration in soil micropore water',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_TEMP_vr(beg_col:end_col,1:JZ)         !TCS_vr(1:JZ,NY,NX)
  call hist_addfld2d(fname='TEMP_vr',units='oC',type2d='levsoi',avgflag='A',&
    long_name='soil temperature profile',ptr_col=data2d_ptr)      

  data2d_ptr => this%h2D_decomp_OStress_vr(beg_col:end_col,1:JZ)         !
  call hist_addfld2d(fname='Decomp_OStress_vr',units='none',type2d='levsoi',avgflag='A',&
    long_name='decomposition oxygen stress [0->1 weaker]',ptr_col=data2d_ptr,default='inactive')           

  data2d_ptr => this%h2D_RO2Decomp_vr(beg_col:end_col,1:JZ)    
  call hist_addfld2d(fname='RO2Decomp_flx_vr',units='gO2/m2/h',type2d='levsoi',avgflag='A',&
    long_name='Decomposition O2 uptake in soil layers',ptr_col=data2d_ptr,default='inactive')                 

  data2d_ptr => this%h2D_Decomp_temp_FN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Decomp_TEMP_FN_vr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Temeprature dependence of microbial decomposition in soil layers',&
    ptr_col=data2d_ptr,default='inactive')                 

  data2d_ptr => this%h2D_FracLitMix_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FracLitMix_vr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Fraction of litter to mixed with the next layer (>0 downward mixing)',ptr_col=data2d_ptr,&
    default='inactive')                 

  data2d_ptr => this%h2D_Decomp_Moist_FN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Decomp_Moist_FN_vr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Moisture dependence of microbial decomposition in soil layers',ptr_col=data2d_ptr,&
    default='inactive') 
!-----

  data2d_ptr =>  this%h2D_RootMassC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Root C density profile',ptr_col=data2d_ptr)       

  data2d_ptr =>  this%h2D_RootMassC_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootC_pvr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Root C density profile of different pft',ptr_patch=data2d_ptr)       

  data2d_ptr =>  this%h2D_RootRadialKond2H2O_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='KH2ORadial_pvr',units='m H2O h-1 MPa-1',type2d='levsoi',avgflag='A',&
    long_name='Radial root conductance for water uptake',ptr_patch=data2d_ptr)       

  data2d_ptr =>  this%h2D_RootAxialKond2H2O_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='KH2OAxial_pvr',units='m3 H2O h-1 MPa-1',type2d='levsoi',avgflag='A',&
    long_name='Axial root conductance for water uptake',ptr_patch=data2d_ptr)       

  data2d_ptr =>  this%h2D_VmaxNH4Root_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='VmaxNH4Root_pvr',units='umolN h-1 (gC root)-1',type2d='levsoi',avgflag='A',&
    long_name='Maximum NH4 uptake rate for given pft',ptr_patch=data2d_ptr)       

  data2d_ptr =>  this%h2D_VmaxNO3Root_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='VmaxNO3Root_pvr',units='umolN h-1 (gC root)-1',type2d='levsoi',avgflag='A',&
    long_name='Maximum NO3 uptake rate for given pft',ptr_patch=data2d_ptr)       

  data2d_ptr =>  this%h2D_RootMassN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Root N density profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_RootMassP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Root P density profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_DOC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='DOC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='DOC profile',ptr_col=data2d_ptr)      

  data2d_ptr =>  this%h2D_DON_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='DON_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='DON profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_DOP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='DOP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='DOP profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_SoilRest4RootGroth_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='SoilRest4RootGroth_vr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='Soil resistance for root penetration',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_acetate_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='acetate_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Acetate profile',ptr_col=data2d_ptr,default='inactive')      

  data1d_ptr => this%h1D_tDOC_soil_col(beg_col:end_col)    
  call hist_addfld1d(fname='DOC_soil_col',units='gC/m2',avgflag='A',&
    long_name='DOC in soil',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tDON_soil_col(beg_col:end_col)    
  call hist_addfld1d(fname='DON_soil_col',units='gN/m2',avgflag='A',&
    long_name='DON in soil',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_RCH4Oxi_aero_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RCH4_AMOX_litr_col',units='gC/m2/h',avgflag='A',&
    long_name='Aerobic CH4 oxidation in litter layer',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_RCH4Oxi_anmo_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RCH4_ANMO_litr_col',units='gC/m2/h',avgflag='A',&
    long_name='Anaerobic CH4 oxidation in litter layer',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_RCH4Oxi_aero_col(beg_col:end_col)
  call hist_addfld1d(fname='RCH4Oxi_aero_col',units='gC/m2/h',avgflag='A',&
    long_name='Aerobic CH4 oxidation integrated over all layers',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_RCH4Oxi_anmo_col(beg_col:end_col)
  call hist_addfld1d(fname='RCH4Oxi_anmo_col',units='gC/m2/h',avgflag='A',&
    long_name='Anaerobic CH4 oxidation integrated over all layers',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_RCH4ProdHydrog_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RCH4ProdHg_litr_col',units='gC/m2/h',avgflag='A',&
    long_name='Hydrogenotrophic CH4 produciton in litter layer',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_RCH4ProdAcetcl_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='RCH4ProdAcet_litr_col',units='gC/m2/h',avgflag='A',&
    long_name='Acetoclastic CH4 produciton in litter layer',ptr_col=data1d_ptr,default='inactive')      

  data1d_ptr => this%h1D_tDOP_soil_col(beg_col:end_col)    
  call hist_addfld1d(fname='DOP_soil_col',units='gP/m2',avgflag='A',&
    long_name='DOP in soil',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_tAcetate_soil_col(beg_col:end_col)    
  call hist_addfld1d(fname='Acetate_soil_col',units='gC/m2',avgflag='A',&
    long_name='Acetate in soil',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_DOC_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='DOC_litr_col',units='gC/m3',avgflag='A',&
    long_name='DOC in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%h1D_DON_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='DON_litr_col',units='gN/m3',avgflag='A',&
    long_name='DON in litter',ptr_col=data1d_ptr,default='inactive')       

  data1d_ptr => this%h1D_DOP_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='DOP_litr_col',units='gP/m3',avgflag='A',&
    long_name='DOP in litter',ptr_col=data1d_ptr,default='inactive')       

  data1d_ptr => this%h1D_acetate_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='Acetate_litr_col',units='gC/m3',avgflag='A',&
    long_name='Acetate in litter',ptr_col=data1d_ptr,default='inactive')      

!------
  data2d_ptr =>  this%h2D_AeroHrBactC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrBacterC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic bacteria C profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_AeroHrFungC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrFungiC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic fungi C profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_faculDenitC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Facult_denitrifierC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Facultative denitrifier C biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_fermentorC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermentorC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Fermentor C biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_fermentor_frac_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Fermentor_frac_vr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Fraction of microbial C in fermentor',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_acetometh_frac_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='AcetoMethgen_frac_vr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Fraction of microbial C in acetoclastic methanogen',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_hydrogMeth_frac_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HydrogMethgen_frac_vr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Fraction of microbial C in hydrogenotrohpic methanogen',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_acetometgC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Acetic_methanogenC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aceticlastic methanogen C biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_aeroN2fixC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_N2fixerC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic N2 fixer C biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_Gas_Pressure_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='GAS_PRESSURE_vr',units='Pa',type2d='levsoi',avgflag='A',&
    long_name='Soil gas pressure profile',ptr_col=data2d_ptr)       

  data2d_ptr =>  this%h2D_CO2_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CO2_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous CO2 profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_CH4_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous CH4 profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_H2_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='H2_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous H2 profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_Ar_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ar_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous Ar profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_O2_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='O2_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous O2 profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_N2_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='N2_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous N2 profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_N2O_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='N2O_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous N2O profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_NH3_Gas_ppmv_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='NH3_gas_ppmv_vr',units='ppmv',type2d='levsoi',avgflag='A',&
    long_name='Equivalent soil gaseous NH3 profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_anaeN2FixC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Anaerobic_N2fixerC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Anaerobic N2 fixer C biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_NH3OxiBactC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ammonia_OxidizerBactC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Ammonia oxidize bacteria C biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_NO2OxiBactC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nitrie_OxidizerBactC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Nitrite oxidize bacteria C profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_CH4AeroOxiC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_methanotrophC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic methanotroph C biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_H2MethogenC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Hygrogen_methanogenC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Hydrogenotrophic methanogen C biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_TSolidOMActC_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='SoildOM_Act_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
    long_name='Active solid organic C in soil layer',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_TSolidOMActCDens_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='SoildOM_Act_Dens_vr',units='gC/gC',type2d='levsoi',avgflag='A',&
    long_name='Active solid organic C Density in soil layer',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_tOMActCDens_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tOMActC_Dens_vr',units='gC/gC',type2d='levsoi',avgflag='A',&
    long_name='Active heterotrophic microbial C Density in soil layer',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_RCH4ProdHydrog_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Prod_hydro_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved hydrogenotrophic CH4 production rate',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_RCH4ProdAcetcl_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Prod_aceto_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved acetoclastic CH4 production rate',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_RCH4Oxi_aero_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Oxi_Aero_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved aerobic CH4 oxidation rate',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_RCH4Oxi_anmo_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='CH4Oxi_ANMO_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved anaerobic CH4 oxidation rate',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_RFerment_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Fermentation_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved fermentation rate',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_RootAR_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootAR_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved root respiration rate',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_RootAR2soil_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootAR2Soil_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved root respiratory CO2 rate to soil',ptr_col=data2d_ptr,default='inactive')      

  data2d_ptr =>  this%h2D_RootAR2Root_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootAR2Root_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved root respiratory CO2 rate to roots',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_nh3oxi_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='NH3Oxid_vr',units='gN/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved NH3 oxidation rate',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_n2oprod_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='N2OProd_vr',units='gN/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved total N2O production rate',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_Eco_HR_CO2_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HR_CO2_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved heterotrophic respiration rate',ptr_col=data2d_ptr,default='inactive')      

  data2d_ptr =>  this%h2D_Gchem_CO2_prod_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Gchem_CO2_Prod_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Vertically resolved geochemical CO2 production rate',ptr_col=data2d_ptr,default='inactive')       
!------
  data2d_ptr =>  this%h2D_AeroHrBactN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrBacterN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic bacteria N profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_AeroHrFungN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrFungiN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic fungi N profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_faculDenitN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Facult_denitrifierN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Facultative denitrifier N biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_fermentorN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermentorN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Fermentor N biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_acetometgN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Acetic_methanogenN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aceticlastic methanogen N biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_aeroN2fixN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_N2fixerN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic N2 fixer N biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_anaeN2FixN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Anaerobic_N2fixerN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Anaerobic N2 fixer N biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_NH3OxiBactN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ammonia_OxidizerBactN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Ammonia oxidize bacteria N biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_NO2OxiBactN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nitrie_OxidizerBactN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Nitrite oxidize bacteria N profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_CH4AeroOxiN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_methanotrophN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic methanotroph N biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_tRespGrossHeterUlm_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tRespGrossHeterUlm_vr',units='gC/m2/h',type2d='levsoi',avgflag='A',&
    long_name='Total oxygen-unlimited gross respiraiton by heterotrophs',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_tRespGrossHeter_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='tRespGrossHeter_vr',units='gC/m2/h',type2d='levsoi',avgflag='A',&
    long_name='Total oxygen gross respiraiton by heterotrophs',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_H2MethogenN_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Hygrogen_methanogenN_vr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Hydrogenotrophic methanogen N biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_MicrobAct_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='MicrobAct_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved respiration-based microbial acitivity for hydrolysis',&
    ptr_col=data2d_ptr,default='inactive')       

  DO jj=1,jcplx
    data2d_ptr =>  this%h3D_HydrolCSOMCps_vr(beg_col:end_col,1:JZ,jj)
    call hist_addfld2d(fname='HydrolCSOM_'//trim(micpar%cplxname(jj))//'_cplx_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
      long_name='Layer resolved carbon hydrolysis in '//trim(micpar%cplxname(jj))//' complex',&
      ptr_col=data2d_ptr,default='inactive')       

    data2d_ptr =>  this%h3D_SOC_Cps_vr(beg_col:end_col,1:JZ,jj)
    call hist_addfld2d(fname='SOC_'//trim(micpar%cplxname(jj))//'_cplx_vr',units='gC/m2',type2d='levsoi',avgflag='A',&
      long_name='Layer resolved carbon in '//trim(micpar%cplxname(jj))//' complex',&
      ptr_col=data2d_ptr,default='inactive')       

    data2d_ptr =>  this%h3D_SOMHydrylScalCps_vr(beg_col:end_col,1:JZ,jj)
    call hist_addfld2d(fname='SOMHydrlScal_'//trim(micpar%cplxname(jj))//'_cplx_vr',units='none',type2d='levsoi',avgflag='A',&
      long_name='Layer resolved SOM hydrolysis scalar in '//trim(micpar%cplxname(jj))//' complex',&
      ptr_col=data2d_ptr,default='inactive')       

    data2d_ptr =>  this%h3D_MicrobActCps_vr(beg_col:end_col,1:JZ,jj)
    call hist_addfld2d(fname='MicrobAct_'//trim(micpar%cplxname(jj))//'_cplx_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
      long_name='Layer resolved respiration-based microbial acitivity for hydrolysis in '//trim(micpar%cplxname(jj))//' complex',&
      ptr_col=data2d_ptr,default='inactive')       
  ENDDO

  DO jj=1,micpar%FG_guilds_heter(micpar%mid_HeterAerobBacter)

  enddo

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_Aerob_Fungi)   
  ENDDO

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_Facult_DenitBacter) 
  ENDDO

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_HeterAerobN2Fixer)  

  ENDDO

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_HeterAnaerobN2Fixer) 
  ENDDO

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_fermentor)
  ENDDO

  DO JJ=1,micpar%FG_guilds_heter(micpar%mid_HeterAcetoCH4GenArchea) 

  ENDDO

  DO JJ=1,micpar%FG_guilds_autor(micpar%mid_AutoH2GenoCH4GenArchea)  

  ENDDO

  DO JJ=1,micpar%FG_guilds_autor(micpar%mid_AutoAmmoniaOxidBacter)  

  ENDDO

  DO JJ=1,micpar%FG_guilds_autor(micpar%mid_AutoNitriteOxidBacter)  

  ENDDO

  DO JJ=1,micpar%FG_guilds_autor(micpar%mid_AutoAeroCH4OxiBacter)
  
  ENDDO

  data2d_ptr =>  this%h2D_RDECOMPC_SOM_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RDecompC_SOM_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved Hydrolysis of solid OM C',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_RDECOMPC_BReSOM_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RDecompC_BReSOM_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved Hydrolysis of microbial residual OM C',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_RDECOMPC_SorpSOM_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RDecompC_SorpSOM_vr',units='gC/m2/hr',type2d='levsoi',avgflag='A',&
    long_name='Layer resolved Hydrolysis of adsorbed OM C',ptr_col=data2d_ptr,default='inactive')       

!------
  data2d_ptr =>  this%h2D_AeroHrBactP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrBacterP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic P biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_AeroHrFungP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_HetrFungiP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic fungi P profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_faculDenitP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Facult_denitrifierP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Facultative denitrifier P biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_fermentorP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='FermentorP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Fermentor P biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_acetometgP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Acetic_methanogenP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aceticlastic methanogen P biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_aeroN2fixP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_N2fixerP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic N2 fixer P biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_anaeN2FixP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Anaerobic_N2fixerP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Anaerobic N2 fixer P biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_NH3OxiBactP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Ammonia_OxidizerBactP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Ammonia oxidize bacteria P biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_NO2OxiBactP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Nitrie_OxidizerBactP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Nitrite oxidize bacteria P profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_CH4AeroOxiP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Aerobic_methanotrophP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Aerobic methanotroph P biomass profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_H2MethogenP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='Hygrogen_methanogenP_vr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Hydrogenotrophic methanogen P biomass profile',ptr_col=data2d_ptr,default='inactive')       
!---

  data2d_ptr =>  this%h2D_MicroBiomeE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='MicroBiomE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Total micorobial elemental biomass in litter',ptr_col=data2d_ptr)

  data2d_ptr =>  this%h2D_AeroHrBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_HetrBacterE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic heterotrophic bacterial elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_AeroHrFungE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_HetrFungiE_litr',units='gP/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic heterotrophic fungi elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_faculDenitE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Facult_denitrifierE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Facultative denitrifier elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_fermentorE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='FermentorE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Fermentor elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_acetometgE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Acetic_methanogenE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aceticlastic methanogen elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_aeroN2fixE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_N2fixerE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic N2 fixer elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_anaeN2FixE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Anaerobic_N2fixerE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Anaerobic N2 fixer elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_NH3OxiBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Ammonia_OxidizerBactE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Ammonia oxidize bacteria elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_NO2OxiBactE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Nitrie_OxidizerBactE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Nitrite oxidize bacteria elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_CH4AeroOxiE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Aerobic_methanotrophE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Aerobic methanotroph elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr =>  this%h2D_H2MethogenE_litr_col(beg_col:end_col,1:NumPlantChemElms)
  call hist_addfld2d(fname='Hygrogen_methanogenE_litr',units='g/m2',type2d='elements',avgflag='A',&
    long_name='Hydrogenotrophic methanogen elemental biomass in litter',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_HeatFlow_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HeatFlow_vr',units='MJ m-3 hr-1',type2d='levsoi',avgflag='A',&
    long_name='soil heat flow profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_HeatUptk_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='HeatUptk_vr',units='MJ m-3 hr-1',type2d='levsoi',avgflag='A',&
    long_name='soil heat flow by plant water uptake (<0 into roots)',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_VSPore_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='VSPore_vr',units='m3 pore/m3 soil',type2d='levsoi',avgflag='A',&
    long_name='Volumetric soil porosity',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_rVSM_vr(beg_col:end_col,1:JZ)        !ThetaH2OZ_vr(1:JZ,NY,NX)
  call hist_addfld2d(fname='rWatFLP_vr',units='m3 H2O/m3 soil pore',type2d='levsoi',avgflag='A',&
    long_name='Fraction of soil porosity filled by water (relative saturation)',ptr_col=data2d_ptr)       

  data2d_ptr => this%h2D_FLO_MICP_vr(beg_col:end_col,1:JZ)      
  call hist_addfld2d(fname='MicPFlo_vr',units='mm H2O h-1',type2d='levsoi',avgflag='A',&
    long_name='Micropore water flow (>0) into soil layer',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_FLO_MACP_vr(beg_col:end_col,1:JZ)        !ThetaH2OZ_vr(1:JZ,NY,NX)
  call hist_addfld2d(fname='MacPFlo_vr',units='mm H2O h-1',type2d='levsoi',avgflag='A',&
    long_name='Macropore water flow (>0) into soil',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_rVSICE_vr(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='rIceFLP_vr',units='m3 ice/m3 soil pore',type2d='levsoi',avgflag='A',&
    long_name='fraction of soil porosity filled by ice',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_PSI_vr(beg_col:end_col,1:JZ)         
  call hist_addfld2d(fname='PSI_vr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='soil matric pressure+osmotic pressure',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_PsiO_vr(beg_col:end_col,1:JZ)         
  call hist_addfld2d(fname='PsiO_vr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='soil osmotic pressure',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RootH2OUP_vr(beg_col:end_col,1:JZ)
  call hist_addfld2d(fname='RootH2OUptake_vr',units='mmH2O/hr',type2d='levsoi',avgflag='A',&
    long_name='soil water taken up by root (<0 into roots)',ptr_col=data2d_ptr,default='inactive')       
  
  data2d_ptr => this%h2D_cNH4t_vr(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='cNH4t_vr',units='gN/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='soil NH4x concentration',ptr_col=data2d_ptr,default='inactive')      

  data2d_ptr => this%h2D_cNO3t_vr(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='cNO3t_vr',units='gN/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='Soil NO3+NO2 concentration',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_cPO4_vr(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='cPO4_vr',units='gP/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='soil dissolved PO4 concentration',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_cEXCH_P_vr(beg_col:end_col,1:JZ)     
  call hist_addfld2d(fname='cEXCH_P_vr',units='gP/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='total exchangeable soil PO4 concentration',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_microb_N2fix_vr(beg_col:end_col,1:JZ)     
  call hist_addfld2d(fname='Free_N2Fix_vr',units='ugN/m3 h-1',type2d='levsoi',avgflag='A',&
    long_name='Free dizotrophic N2 fixation in soil',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_TEMP_vr(beg_col:end_col,1:JZ)  
  call hist_addfld2d(fname='TMAX_SOIL_vr',units='oC',type2d='levsoi',avgflag='X',&
    long_name='Soil maximum temperature profile',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_TEMP_vr(beg_col:end_col,1:JZ)  
  call hist_addfld2d(fname='TMIN_SOIL_vr',units='oC',type2d='levsoi',avgflag='M',&
    long_name='Soil minimum temperature profile',ptr_col=data2d_ptr,default='inactive')       

  data1d_ptr => this%h1D_decomp_ostress_litr_col(beg_col:end_col)    
  call hist_addfld1d(fname='Decomp_OStress_LITR',units='none',avgflag='A',&
    long_name='Decomposition O2 stress in litter layer [0->1: weaker]',ptr_col=data1d_ptr,&
    default='inactive')       

  data1d_ptr => this%h1D_Decomp_temp_FN_litr_col(beg_col:end_col)    
  call hist_addfld1d(fname='Decomp_TEMP_FN_LITR',units='none',avgflag='A',&
    long_name='Decomposition temperature sensitivity in litter layer',ptr_col=data1d_ptr,&
    default='inactive')       

  data1d_ptr => this%h1D_FracLitMix_litr_col(beg_col:end_col)    
  call hist_addfld1d(fname='FracLitMix_LITR',units='none',avgflag='A',&
    long_name='Fraction of surface litter layer to be mixed downward',ptr_col=data1d_ptr,&
    default='inactive')       

  data1d_ptr => this%h1D_Decomp_Moist_FN_litr_col(beg_col:end_col)
  call hist_addfld1d(fname='Decomp_Moist_FN_LITR',units='none',avgflag='A',&
    long_name='Decomposition moisture sensitivity in litter layer',ptr_col=data1d_ptr,&
    default='inactive')       

  data1d_ptr => this%h1D_RO2Decomp_litr_col(beg_col:end_col)    
  call hist_addfld1d(fname='RO2Decomp_flx_LITR',units='gO2/m2/h',avgflag='A',&
    long_name='Decomposition O2 uptake in litter layer ',ptr_col=data1d_ptr,default='inactive')       

  data1d_ptr => this%h1D_TSolidOMActC_litr_col(beg_col:end_col)    
  call hist_addfld1d(fname='SoildOM_Act_litr',units='gC/m2',avgflag='M',&
    long_name='Active solid OM in litter',ptr_col=data1d_ptr,default='inactive')       

  data1d_ptr => this%h1D_tOMActCDens_litr_col(beg_col:end_col)    
  call hist_addfld1d(fname='tOMActC_Dens_litr',units='gC/gC',avgflag='M',&
    long_name='Active heterotrophic microbial C density in litter',ptr_col=data1d_ptr,&
    default='inactive')       

  data1d_ptr => this%h1D_TSolidOMActCDens_litr_col(beg_col:end_col)    
  call hist_addfld1d(fname='SoildOM_Act_Dens_litr',units='gC/gC',avgflag='M',&
    long_name='Active solid OM density in litter',ptr_col=data1d_ptr,default='inactive')       

  data2d_ptr => this%h2D_ElectricConductivity_vr(beg_col:end_col,1:JZ)     
  call hist_addfld2d(fname='ElectricConductivity_vr',units='dS m-1',type2d='levsoi',avgflag='A',&
    long_name='electrical conductivity',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_HydCondSoil_vr(beg_col:end_col,1:JZ)     
  call hist_addfld2d(fname='HydCondSoil_vr',units='m MPa-1 h-1',type2d='levsoi',avgflag='A',&
    long_name='Vertical hydraulic conductivity',ptr_col=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_PSI_RT_pvr(beg_ptc:end_ptc,1:JZ)  
  call hist_addfld2d(fname='PSI_RT_pvr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='Root total water potential of each pft',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RootH2OUptkStress_pvr(beg_ptc:end_ptc,1:JZ)  
  call hist_addfld2d(fname='RootH2OUptkStress_pvr',units='m3 m-2 h-1',type2d='levsoi',avgflag='A',&
    long_name='Rate indicated root water uptake stress of each pft (>0 hydraulic stress)',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RootH2OUptk_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootH2OUptk_pvr',units='mm H2O h-1',type2d='levsoi',avgflag='A',&
    long_name='Plant root water uptake from soil (>0 release water to soil)',ptr_patch=data2d_ptr)

  data2d_ptr => this%h2D_RootMaintDef_CO2_pvr(beg_ptc:end_ptc,1:JZ)  
  call hist_addfld2d(fname='RootMaintDef_CO2_pvr',units='g CO2 m-2 h-1',type2d='levsoi',avgflag='A',&
    long_name='Plant root maintenance deficit each pft (<0 deficit)',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RootSurfAreaPP_pvr(beg_ptc:end_ptc,1:JZ)  
  call hist_addfld2d(fname='RootSurfAreaPP_pvr',units='m2 surface',type2d='levsoi',avgflag='A',&
    long_name='Root surface area per plant (for nutrient uptake)',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_ROOTNLim_rpvr(beg_ptc:end_ptc,1:JZ)  
  call hist_addfld2d(fname='RootNlim_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Plant root nitrogen limitation for each pft',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_ROOTPLim_rpvr(beg_ptc:end_ptc,1:JZ)  
  call hist_addfld2d(fname='RootPlim_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Plant root phosphorus limitation for each pft',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_ROOT_OSTRESS_pvr(beg_ptc:end_ptc,1:JZ)  
  call hist_addfld2d(fname='Root_OXYSTRESS_pvr',units='None',type2d='levsoi',avgflag='A',&
    long_name='Root Oxygen stress profile [0->1 weaker stress]',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Root1stStrutC_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='RootC_1st_pvr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Primary root structural biomass C density',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RootNonstBConc_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='RootNonstBConc_pvr',units='g/gC',type2d='levsoi',avgflag='A',&
    long_name='Primary root nonstructural biomass density',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Root1stStrutN_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='RootN_1st_pvr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Primary root structural biomass N density',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Root1stStrutP_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='RootP_1st_pvr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Primary root structural biomass P density',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Root2ndStrutC_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='RootC_2nd_pvr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='Secondary root structural biomass C density',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Root2ndStrutN_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='RootN_2nd_pvr',units='gN/m3',type2d='levsoi',avgflag='A',&
    long_name='Secondary root structural biomass N density',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Root2ndStrutP_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='RootP_2nd_pvr',units='gP/m3',type2d='levsoi',avgflag='A',&
    long_name='Secondary root structural biomass P density',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Root1stAxesNumL_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='Root1st_AxesNumL_pvr',units='1/d2',type2d='levsoi',avgflag='A',&
    long_name='Primary root axes number in soil layer',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_Root2ndAxesNumL_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='Root2nd_AxesNumL_pvr',units='1/d2',type2d='levsoi',avgflag='A',&
    long_name='Secondary root axes number in soil layer',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RootKond2H2O_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='RootKond2H2O_pvr',units='x1.e7 m s-1 MPa-1',type2d='levsoi',avgflag='A',&
    long_name='Total root conductance to water uptake in soil layer',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_prtUP_NH4_pvr(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='prtUP_NH4_pvr',units='gN/m3/hr',type2d='levsoi',avgflag='A',&
    long_name='Root uptake of NH4',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_prtUP_NO3_pvr(beg_ptc:end_ptc,1:JZ)      
  call hist_addfld2d(fname='prtUP_NO3_pvr',units='gN/m3/hr',type2d='levsoi',&
    avgflag='A',long_name='Root uptake of NO3',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_prtUP_PO4_pvr(beg_ptc:end_ptc,1:JZ)     
  call hist_addfld2d(fname='prtUP_PO4_pvr',units='gP/m3/hr',type2d='levsoi',avgflag='A',&
    long_name='Root uptake of PO4',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_DNS_RT_pvr(beg_ptc:end_ptc,1:JZ)       
  call hist_addfld2d(fname='RootLDS_pvr',units='cm/cm3',type2d='levsoi',avgflag='A',&
    long_name='Root length density (including root hair)',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RootNutupk_fClim_pvr(beg_ptc:end_ptc,1:JZ)       
  call hist_addfld2d(fname='RootNutUptk_fClim_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='C-availability for root nutrient uptake limitation, 0->1 stronger limitation',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RootNutupk_fNlim_pvr(beg_ptc:end_ptc,1:JZ)       
  call hist_addfld2d(fname='RootNutUptk_fNlim_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='N-limitation for root nutrient uptake limitation, 0->1 stronger limitation',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RootNutupk_fPlim_pvr(beg_ptc:end_ptc,1:JZ)       
  call hist_addfld2d(fname='RootNutUptk_fPlim_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='P-limitation for root nutrient uptake limitation, 0->1 stronger limitation',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_RootNutupk_fProtC_pvr(beg_ptc:end_ptc,1:JZ)       
  call hist_addfld2d(fname='RootNutUptk_fProtC_pvr',units='-',type2d='levsoi',avgflag='A',&
    long_name='Protein-limitation for root nutrient uptake capacity, 0->1 weaker limitation',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_fTRootGro_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootGRO_TEMP_FN_pvr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Root growth temperature dependence function',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_CanopyLAIZ_plyr(beg_ptc:end_ptc,1:NumCanopyLayers)
  call hist_addfld2d(fname='CanopyLAIZ_plyr',units='m2/m2',type2d='levcan',avgflag='A',&
    long_name='Vertically distributed leaf area',ptr_patch=data2d_ptr,default='inactive')       

  data2d_ptr => this%h2D_fRootGrowPSISense_pvr(beg_ptc:end_ptc,1:JZ)
  call hist_addfld2d(fname='RootGRO_PSI_FN_pvr',units='none',type2d='levsoi',avgflag='A',&
    long_name='Root growth moisture dependence function',ptr_patch=data2d_ptr,default='inactive')       
  ![terminate]
  do nbr=1,MaxNumBranches
    data2d_ptr => this%h3D_PARTS_ptc(beg_ptc:end_ptc,1:NumOfPlantMorphUnits,nbr)
    write(fieldname,'(I2.2)')nbr
    call hist_addfld2d(fname='C_PARTS_brch_'//trim(fieldname),units='none',&
      type2d='pmorphunits',avgflag='A',&
      long_name='C allocation to different morph unit in branch '//trim(fieldname),ptr_patch=data2d_ptr,default='inactive')      
  enddo
  end subroutine init_hist_data

!----------------------------------------------------------------------
  subroutine hist_update(this,I,J,bounds)
  use TracerPropMod, only : GramPerHr2umolPerSec
  implicit none
  class(histdata_type) :: this
  integer, intent(in) :: I,J
  type(bounds_type), intent(in) :: bounds
  integer :: ncol,nptc
  integer :: L,NZ,NY,NX,KN,NB,NR,NB1,K
  real(r8) :: micBE(1:NumPlantChemElms)
  real(r8) :: DOM(idom_beg:idom_end)
  real(r8),parameter :: secs1hour=3600._r8
  real(r8),parameter :: MJ2W=1.e6_r8/secs1hour
  real(r8),parameter :: m2mm=1000._r8
  real(r8),parameter :: million=1.e6_r8
  character(len=*), parameter :: subname='hist_update'
  character(len=15) :: grow_stage_str(11)=(/'Planting      ','Emergence     ','Floral_init   ', &
                                            'Jointing      ','Elongation    ','Heading       ', &
                                            'Anthesis      ','Seed_fill     ','See_no_set    ', &
                                            'Seed_mass_set ','End_seed_fill '/)
  real(r8) :: DVOLL,SOMC(jcplx)                                            
  integer :: jj
  call PrintInfo('beg '//subname)
  DO NX=bounds%NHW,bounds%NHE   
    DO NY=bounds%NVN,bounds%NVS
      
      ncol=get_col(NY,NX)
      this%h1D_cumFIRE_CO2_col(ncol)        =  CO2byFire_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_cumFIRE_CH4_col(ncol)        =  CH4byFire_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_cNH4_LITR_col(ncol)        =  safe_adb(trcs_solml_vr(ids_NH4,0,NY,NX)+&
        natomw*trcx_solml_vr(idx_NH4,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      this%h1D_cNO3_LITR_col(ncol)        =  safe_adb(trcs_solml_vr(ids_NO3,0,NY,NX)+&
        trcs_solml_vr(ids_NO2,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      
      this%h1D_ECO_HVST_N_col(ncol)   = EcoHavstElmnt_CumYr_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NET_N_MIN_col(ncol)    = -NetNH4Mineralize_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tLITR_P_col(ncol)      = tLitrOM_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_HUMUS_C_col(ncol)      = tHumOM_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_HUMUS_N_col(ncol)      = tHumOM_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_HUMUS_P_col(ncol)      = tHumOM_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_AMENDED_P_col(ncol)    = FerP_Flx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tLITRf_C_FLX_col(ncol) = LiterfalOrgM_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tLITRf_N_FLX_col(ncol) = LiterfalOrgM_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tLITRf_P_FLX_col(ncol) = LiterfalOrgM_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tEXCH_PO4_col(ncol)        = tHxPO4_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DOP_FLX_col(ncol)      = HydroSufDOPFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DOP_FLX_col(ncol)      = HydroSubsDOPFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DIP_FLX_col(ncol)      = HydroSufDIPFlx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DIP_FLX_col(ncol)      = HydroSubsDIPFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)  
      this%h1D_HeatFlx2Grnd_col(ncol)     = HeatFlx2Grnd_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

      this%h1D_CanSWRad_col(ncol)         = MJ2W*RadSW_Canopy_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RadSW_Grnd_col(ncol)       = MJ2W*RadSWGrnd_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Qinfl2soi_col(ncol)        = m2mm*Qinflx2Soil_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Qdrain_col(ncol)           = m2mm*QDrain_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

      this%h1D_SUR_DON_FLX_col(ncol)      = HydroSufDONFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DON_FLX_col(ncol)      = HydroSubsDONFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSALT_DISCHG_FLX_col(ncol) = HydroIonFlx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DIN_FLX_col(ncol)      = HydroSufDINFlx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DIN_FLX_col(ncol)      = HydroSubsDINFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DOC_FLX_col(ncol)      = HydroSufDOCFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DOC_FLX_col(ncol)      = HydroSubsDOCFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DIC_FLX_col(ncol)      = HydroSufDICFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUB_DIC_FLX_col(ncol)      = HydroSubsDICFlx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SUR_DIP_FLX_col(ncol)      = HydroSufDIPFlx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tPREC_P_col(ncol)          = tXPO4_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tMICRO_P_col(ncol)         = tMicBiome_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_PO4_FIRE_col(ncol)         = PO4byFire_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_cPO4_LITR_col(ncol)        = safe_adb(trcs_solml_vr(ids_H2PO4,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      this%h1D_cEXCH_P_LITR_col(ncol)     = patomw*safe_adb(trcx_solml_vr(idx_HPO4,0,NY,NX)+&
        trcx_solml_vr(idx_H2PO4,0,NY,NX),VLSoilMicPMass_vr(0,NY,NX)*million)
      this%h1D_ECO_HVST_P_col(ncol) = EcoHavstElmnt_CumYr_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NET_P_MIN_col(ncol)  = -NetPO4Mineralize_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_PSI_SURF_col(ncol)   = PSISoilMatricP_vr(0,NY,NX)
      this%h1D_SURF_ELEV_col(ncol)  = -CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)+DLYR_3D(3,0,NY,NX)
      this%h1D_tLITR_N_col(ncol)    = tLitrOM_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_AMENDED_N_col(ncol)  = FertN_Flx_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tNH4X_col(ncol)      = tNH4_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tNO3_col(ncol)       = tNO3_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tRAD_col(ncol)       = TRAD_col(NY,NX)
      if(this%h1D_tNH4X_col(ncol)<0._r8)then
        write(*,*)'negative tNH4X',this%h1D_tNH4X_col(ncol),this%h1D_tNO3_col(ncol)
        stop
      endif      
      this%h1D_tMICRO_N_col(ncol)         = tMicBiome_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_TEMP_LITR_col(ncol)        = TCS_vr(0,NY,NX)
      this%h1D_TEMP_surf_col(ncol)        = FracSurfByLitR_col(NY,NX)*TCS_vr(0,NY,NX)+(1._r8-FracSurfByLitR_col(NY,NX))*TCS_vr(NU_col(NY,NX),NY,NX)
      if(VcumSnowWE_col(NY,NX)<=ZEROS(NY,NX))then
        this%h1D_TEMP_SNOW_col(ncol)   = spval
      else
        this%h1D_TEMP_SNOW_col(ncol)   = TCSnow_snvr(1,NY,NX)
      endif
      this%h1D_FracBySnow_col(ncol) = FracSurfAsSnow_col(NY,NX)
      this%h1D_FracByLitr_col(ncol) = FracSurfByLitR_col(NY,NX)
      this%h1D_tLITR_C_col(ncol)    = tLitrOM_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      
      this%h1D_AMENDED_C_col(ncol)        = AmendC_CumYr_flx_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tMICRO_C_col(ncol)         = tMicBiome_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSoilOrgC_col(ncol)        = tSoilOrgM_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSoilOrgN_col(ncol)        = tSoilOrgM_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSoilOrgP_col(ncol)        = tSoilOrgM_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_OMC_LITR_col(ncol)         = SoilOrgM_vr(ielmc,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_OMN_LITR_col(ncol)         = SoilOrgM_vr(ielmn,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_OMP_LITR_col(ncol)         = SoilOrgM_vr(ielmp,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ATM_CO2_col(ncol)          = CO2E_col(NY,NX)
      this%h1D_ATM_CH4_col(ncol)          = CH4E_col(NY,NX)
      this%h1D_NBP_col(ncol)              = Eco_NBP_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_HVST_C_col(ncol)       = EcoHavstElmnt_CumYr_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_LAI_col(ncol)          = CanopyLeafArea_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_SAI_col(ncol)          = StemArea_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Eco_GPP_CumYr_col(ncol)    = Eco_GPP_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_RA_col(ncol)           = Eco_AutoR_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Eco_NPP_CumYr_col(ncol)    = Eco_NPP_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Eco_HR_CumYr_col(ncol)     = Eco_HR_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Eco_HR_CO2_col(ncol)       = ECO_HR_CO2_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
!      this%h1D_Eco_HR_CH4_col(ncol)       = ECO_HR_CH4_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tDIC_col(ncol)             = DIC_mass_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSTANDING_DEAD_C_col(ncol) = StandingDeadStrutElms_col(ielmc,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSTANDING_DEAD_N_col(ncol) = StandingDeadStrutElms_col(ielmn,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSTANDING_DEAD_P_col(ncol) = StandingDeadStrutElms_col(ielmp,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tPRECIP_col(ncol)           = m2mm*QRain_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_ET_col(ncol)           = m2mm*QEvap_CumYr_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_Ar_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_CO2_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_CH4_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_CH4,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_O2_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_O2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_N2_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_N2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_NH3_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_NH3,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_trcg_H2_cumerr_col(ncol)   = trcg_mass_cumerr_col(idg_H2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

      this%h1D_ECO_RADSW_col(ncol)        = MJ2W*Eco_RadSW_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_N2O_LITR_col(ncol)         = trc_solcl_vr(idg_N2O,0,NY,NX)
      this%h1D_NH3_LITR_col(ncol)         = trc_solcl_vr(idg_NH3,0,NY,NX)
      this%h1D_SOL_RADN_col(ncol)         = RadSWSolarBeam_col(NY,NX)*MJ2W
      this%h1D_AIR_TEMP_col(ncol)         = TairK_col(NY,NX)-273.15_r8
      this%h1D_HUM_col(ncol)              = VPK_col(NY,NX)
      this%h1D_PATM_col(ncol)             = PBOT_col(NY,NX)
      this%h1D_WIND_col(ncol)             = WindSpeedAtm_col(NY,NX)/secs1hour
      this%h1D_PREC_col(ncol)             = (RainFalPrec_col(NY,NX)+SnoFalPrec_col(NY,NX))*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Snofall_col(ncol)          = SnoFalPrec_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SOIL_RN_col(ncol)          = HeatByRad2Surf_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_LWSky_col(ncol)            = LWRadSky_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SOIL_LE_col(ncol)          = HeatEvapAir2Surf_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SOIL_H_col(ncol)           = HeatSensAir2Surf_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SOIL_G_col(ncol)           = -(HeatNet2Surf_col(NY,NX)-HeatSensVapAir2Surf_col(NY,NX))*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_RN_col(ncol)           = Eco_NetRad_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_LE_col(ncol)           = Eco_Heat_Latent_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Eco_HeatSen_col(ncol)      = Eco_Heat_Sens_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_ECO_Heat2G_col(ncol)       = Eco_Heat_GrndSurf_col(NY,NX)*MJ2W/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_O2_LITR_col(ncol)          = trc_solcl_vr(idg_O2,0,NY,NX)
      this%h1D_CO2_SEMIS_FLX_col(ncol)    = SurfGasEmiss_all_flx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1D_AR_SEMIS_FLX_col(ncol)     = SurfGasEmiss_all_flx_col(idg_AR,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_AR)
      this%h1D_ECO_CO2_FLX_col(ncol)      = Eco_NEE_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1d_CAN_NEE_col(ncol)          = Canopy_NEE_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1D_CH4_SEMIS_FLX_col(ncol)    = SurfGasEmiss_all_flx_col(idg_CH4,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1D_O2_SEMIS_FLX_col(ncol)     = SurfGasEmiss_all_flx_col(idg_O2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_O2)
      this%h1D_CH4_EBU_flx_col(ncol)      = trcg_ebu_flx_col(idg_CH4,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CH4)
      this%h1D_Ar_EBU_flx_col(ncol)       = trcg_ebu_flx_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_Ar)
      this%h1D_AR_PLTROOT_flx_col(ncol)   = trcg_air2root_flx_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_Ar)
      this%h1D_CH4_PLTROOT_flx_col(ncol)  = trcg_air2root_flx_col(idg_CH4,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CH4)
      this%h1D_CO2_PLTROOT_flx_col(ncol)  = trcg_air2root_flx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1D_O2_PLTROOT_flx_col(ncol)   = trcg_air2root_flx_col(idg_O2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_O2)
      this%h1D_CO2_DIF_flx_col(ncol)      = GasDiff2Surf_flx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
      this%h1D_CH4_DIF_flx_col(ncol)      = GasDiff2Surf_flx_col(idg_CH4,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CH4)      
      this%h1D_Ar_DIF_flx_col(ncol)       = GasDiff2Surf_flx_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_Ar)  
      this%h1D_NH3_DIF_flx_col(ncol)      = GasDiff2Surf_flx_col(idg_NH3,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_NH3)    
      this%h1D_O2_DIF_flx_col(ncol)        =GasDiff2Surf_flx_col(idg_O2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_O2)     
      this%h1D_CO2_TPR_err_col(ncol)      = Gas_Prod_TP_cumRes_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Ar_TPR_err_col(ncol)       = Gas_Prod_TP_cumRes_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_Drain_flx_col(ncol)    = trcs_drainage_flx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_hydloss_flx_col(ncol)  = GasHydroLoss_cumflx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_LITR_col(ncol)         = trc_solcl_vr(idg_CO2,0,NY,NX)
      this%h1D_EVAPN_col(ncol)            = VapXAir2GSurf_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CondGasXSurf_col(ncol)     = CondGasXSurf_col(NY,NX)
      this%h1D_CANET_col(ncol)            = QVegET_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RUNOFF_FLX_col(ncol)       = -QRunSurf_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SEDIMENT_FLX_col(ncol)     = SedmErossLoss_CumYr_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tSWC_col(ncol)             = WatMass_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tHeat_col(ncol)            = HeatStore_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_QDISCHG_FLX_col(ncol)      = QDischarg2WTBL_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_HeatDISCHG_FLX_col(ncol)   = HeatDischar_col(NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_SNOWPACK_col(ncol)         = AZMAX1((VcumSnowWE_col(NY,NX))*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX))
      if(this%h1D_SNOWPACK_col(ncol)>0._r8)then
        this%h1D_SNOWDENS_col(ncol)         = this%h1D_SNOWPACK_col(ncol)/SnowDepth_col(NY,NX)
      else
        this%h1D_SNOWDENS_col(ncol)         =0._r8
      endif
      this%h1D_SURF_WTR_col(ncol)         = ThetaH2OZ_vr(0,NY,NX)
      this%h1D_SURF_ICE_col(ncol)         = ThetaICEZ_vr(0,NY,NX)
      this%h1D_ThetaW_litr_col(ncol)      = safe_adb(VLWatMicP_vr(0,NY,NX),VLitR_col(NY,NX))
      this%h1D_ThetaI_litr_col(ncol)      = safe_adb(VLiceMicP_vr(0,NY,NX),VLitR_col(NY,NX))
      this%h1D_ACTV_LYR_col(ncol)         = -(ActiveLayDepZ_col(NY,NX)-CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX))
      this%h1D_WTR_TBL_col(ncol)          = -(DepzIntWTBL_col(NY,NX)-CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX))
      this%h1D_N2O_SEMIS_FLX_col(ncol)         = SurfGasEmiss_all_flx_col(idg_N2O,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_N2_SEMIS_FLX_col(ncol)         = SurfGasEmiss_all_flx_col(idg_N2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NH3_SEMIS_FLX_col(ncol)         = SurfGasEmiss_all_flx_col(idg_NH3,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_H2_SEMIS_FLX_col(ncol)          = SurfGasEmiss_all_flx_col(idg_H2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_PAR_col(ncol)            = RadPARSolarBeam_col(NY,NX)
      this%h1D_VHeatCap_litr_col(ncol)  = VHeatCapacity_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_AR_WetDep_FLX_col(ncol)  = Gas_WetDeposit_flx_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_WetDep_FLX_col(ncol) = Gas_WetDeposit_flx_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RootXO2_flx_col(ncol)    = RUptkRootO2_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RootN_Fix_col(ncol)      = RootN2Fix_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      call sumMicBiomLayL(0,NY,NX,micBE)      
      this%h2D_MicroBiomeE_litr_col(ncol,1:NumPlantChemElms) =micBE/AREA_3D(3,NU_col(NY,NX),NY,NX)  

      call SumMicbGroup(0,NY,NX,micpar%mid_HeterAerobBacter,MicbE)
      this%h2D_AeroHrBactE_litr_col(ncol,1:NumPlantChemElms) =micBE/AREA_3D(3,NU_col(NY,NX),NY,NX)     !aerobic heterotropic bacteria

      call SumMicbGroup(0,NY,NX,micpar%mid_Aerob_Fungi,MicbE)
      this%h2D_AeroHrFungE_litr_col(ncol,1:NumPlantChemElms) = micBE/AREA_3D(3,NU_col(NY,NX),NY,NX)   !aerobic heterotropic fungi

      call SumMicbGroup(0,NY,NX,micpar%mid_Facult_DenitBacter,MicbE)
      this%h2D_faculDenitE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !facultative denitrifier

      call SumMicbGroup(0,NY,NX,micpar%mid_fermentor,MicbE)
      this%h2D_fermentorE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !fermentor

      call SumMicbGroup(0,NY,NX,micpar%mid_HeterAcetoCH4GenArchea,MicbE)
      this%h2D_acetometgE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !acetogenic methanogen

      call SumMicbGroup(0,NY,NX,micpar%mid_HeterAerobN2Fixer,MicbE)
      this%h2D_aeroN2fixE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !aerobic N2 fixer

      call SumMicbGroup(0,NY,NX,micpar%mid_HeterAnaerobN2Fixer,MicbE)
      this%h2D_anaeN2FixE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)  !anaerobic N2 fixer

      call SumMicbGroup(0,NY,NX,micpar%mid_AutoAmmoniaOxidBacter,MicbE,isauto=.true.)     
      this%h2D_NH3OxiBactE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)

      call SumMicbGroup(0,NY,NX,micpar%mid_AutoNitriteOxidBacter,MicbE,isauto=.true.)     
      this%h2D_NO2OxiBactE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)

      call SumMicbGroup(0,NY,NX,micpar%mid_AutoAeroCH4OxiBacter,MicbE,isauto=.true.)     
      this%h2D_CH4AeroOxiE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)

      call SumMicbGroup(0,NY,NX,micpar%mid_AutoH2GenoCH4GenArchea,MicbE,isauto=.true.)     
      this%h2D_H2MethogenE_litr_col(ncol,1:NumPlantChemElms) = MicbE/AREA_3D(3,NU_col(NY,NX),NY,NX)
      
      call sumDOML(0,NY,NX,DOM)
      
      this%h1D_DOC_LITR_col(ncol)     = safe_adb(DOM(idom_doc),VLitR_col(NY,NX))
      this%h1D_DON_LITR_col(ncol)     = safe_adb(DOM(idom_don),VLitR_col(NY,NX))
      this%h1D_DOP_LITR_col(ncol)     = safe_adb(DOM(idom_dop),VLitR_col(NY,NX))
      this%h1D_acetate_LITR_col(ncol) = safe_adb(DOM(idom_acetate),VLitR_col(NY,NX))

      this%h1D_RCH4ProdHydrog_litr_col(ncol) = RCH4ProdHydrog_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RCH4ProdAcetcl_litr_col(ncol) = RCH4ProdAcetcl_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RCH4Oxi_aero_litr_col(ncol)   = RCH4Oxi_aero_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)   
      this%h1D_RCH4Oxi_aero_col(ncol) =    RCH4Oxi_aero_vr(0,NY,NX)
      this%h1D_RCH4Oxi_ANMO_litr_col(ncol)   = RCH4Oxi_anmo_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)   
      this%h1D_RCH4Oxi_ANMO_col(ncol) =    RCH4Oxi_anmo_vr(0,NY,NX)
      this%h1D_RFermen_litr_col(ncol)        = RFerment_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_NH3oxi_litr_col(ncol)         = RNH3oxi_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_N2Oprod_litr_col(ncol)        = (RN2ODeniProd_vr(0,NY,NX)+RN2ONitProd_vr(0,NY,NX) &
                               +RN2OChemoProd_vr(0,NY,NX)-RN2ORedux_vr(0,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
      
      this%h1D_decomp_OStress_litr_col(ncol)   = OxyDecompLimiter_vr(0,NY,NX)
      this%h1D_MicrobAct_litr_col(ncol)        = TMicHeterActivity_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RO2Decomp_litr_col(ncol)        = RO2DecompUptk_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RDECOMPC_SOM_litr_col(ncol)     = tRHydlySOM_vr(ielmc,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RDECOMPC_BReSOM_litr_col(ncol)  = tRHydlyBioReSOM_vr(ielmc,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_RDECOMPC_SorpSOM_litr_col(ncol) = tRHydlySoprtOM_vr(ielmc,0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tRespGrossHete_litr_col(ncol)   = tRespGrossHeter_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_tRespGrossHeteUlm_litr_col(ncol) = tRespGrossHeterUlm_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

      this%h1D_Decomp_temp_FN_litr_col(ncol)   = TempSensDecomp_vr(0,NY,NX)
      this%h1D_Decomp_moist_FN_litr_col(ncol)  = MoistSensDecomp_vr(0,NY,NX)
      this%h1D_FracLitMix_litr_col(ncol)       = FracLitrMix_vr(0,NY,NX)
      this%h1D_Eco_HR_CO2_litr_col(ncol)       = ECO_HR_CO2_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_TSolidOMActC_litr_col(ncol)     = TSolidOMActC_vr(0,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_TSolidOMActCDens_litr_col(ncol) = safe_adb(TSolidOMActC_vr(0,NY,NX),TSolidOMC_vr(0,NY,NX))
      this%h1D_tOMActCDens_litr_col(ncol)      = safe_adb(tOMActC_vr(0,NY,NX),(TSolidOMC_vr(0,NY,NX)+tOMActC_vr(0,NY,NX)))
      this%h1D_Ar_mass_col(ncol)               = trcg_TotalMass_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Ar_soilMass_col(ncol)           = trcg_soilMass_col(idg_Ar,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_CO2_mass_col(ncol)               = trcg_TotalMass_col(idg_CO2,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
      this%h1D_Gchem_CO2_prod_col(ncol)         = sum(TProd_CO2_geochem_soil_vr(1:JZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
      
      this%h1D_tDOC_soil_col(ncol)=0._r8
      this%h1D_tDON_soil_col(ncol)=0._r8
      this%h1D_tDOP_soil_col(ncol)=0._r8            
      this%h1D_tAcetate_soil_col(ncol)=0._r8      
      this%h1D_FreeNFix_col(ncol)=Micb_N2Fixation_vr(0,NY,NX)

      DO L=1,JZ        
        call SumSolidOML(ielmc,L,NY,NX,SOMC)
        DO jj=1,jcplx
          this%h3D_SOC_Cps_vr(ncol,L,jj) = SOMC(jj)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        ENDDO
        this%h2D_Gas_Pressure_vr(ncol,L) = Soil_Gas_pressure_vr(L,NY,NX)
        this%h2D_CO2_Gas_ppmv_vr(ncol,L) = Soil_Gas_Frac_vr(idg_CO2,L,NY,NX)
        this%h2D_CH4_Gas_ppmv_vr(ncol,L) = Soil_Gas_Frac_vr(idg_CH4,L,NY,NX)
        this%h2D_Ar_Gas_ppmv_vr(ncol,L)  = Soil_Gas_Frac_vr(idg_Ar,L,NY,NX)
        this%h2D_O2_Gas_ppmv_vr(ncol,L)  = Soil_Gas_Frac_vr(idg_O2,L,NY,NX)
        this%h2D_H2_Gas_ppmv_vr(ncol,L)  = Soil_Gas_Frac_vr(idg_H2,L,NY,NX)
        this%h2D_N2O_Gas_ppmv_vr(ncol,L) = Soil_Gas_Frac_vr(idg_N2O,L,NY,NX)
        this%h2D_N2_Gas_ppmv_vr(ncol,L)  = Soil_Gas_Frac_vr(idg_N2,L,NY,NX)
        this%h2D_NH3_Gas_ppmv_vr(ncol,L) = Soil_Gas_Frac_vr(idg_NH3,L,NY,NX)
        
        DVOLL=DLYR_3D(3,L,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
        
        this%h2D_Eco_HR_CO2_vr(ncol,L)    = ECO_HR_CO2_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_Gchem_CO2_prod_vr(ncol,L)= TProd_CO2_geochem_soil_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_SoilRest4RootGroth_vr(ncol,L) = SoilResist4RootPentrate_vr(L,NY,NX)
        if(DVOLL<=1.e-8_r8)cycle

        call sumDOML(L,NY,NX,DOM)
        call sumMicBiomLayL(L,NY,NX,micBE,I,J)    
        this%h1D_tDOC_soil_col(ncol)=this%h1D_tDOC_soil_col(ncol)+DOM(idom_doc)
        this%h1D_tDON_soil_col(ncol)=this%h1D_tDON_soil_col(ncol)+DOM(idom_don)
        this%h1D_tDOP_soil_col(ncol)=this%h1D_tDOP_soil_col(ncol)+DOM(idom_dop)
        this%h1D_tAcetate_soil_col(ncol)=this%h1D_tAcetate_soil_col(ncol)+DOM(idom_acetate)
        this%h2D_QDrainloss_vr(ncol,L)=QDrainloss_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_FermOXYI_vr(ncol,L)=FermOXYI_vr(L,NY,NX)
        this%h2D_DOC_vr(ncol,L)             = DOM(idom_doc)/DVOLL
        this%h2D_DON_vr(ncol,L)             = DOM(idom_don)/DVOLL
        this%h2D_DOP_vr(ncol,L)             = DOM(idom_dop)/DVOLL
        this%h2D_BotDEPZ_vr(ncol,L)         = CumDepz2LayBottom_vr(L,NY,NX)
        this%h2D_acetate_vr(ncol,L)         = DOM(idom_acetate)/DVOLL
        this%h2D_litrC_vr(ncol,L)           = litrOM_vr(ielmc,L,NY,NX)/DVOLL
        this%h2D_litrN_vr(ncol,L)           = litrOM_vr(ielmn,L,NY,NX)/DVOLL
        this%h2D_litrP_vr(ncol,L)           = litrOM_vr(ielmp,L,NY,NX)/DVOLL
        this%h2D_tSOC_vr(ncol,L)            = SoilOrgM_vr(ielmc,L,NY,NX)/DVOLL

        this%h2D_POM_C_vr(ncol,L)           = sum(SolidOM_vr(ielmc,:,micpar%k_POM,L,NY,NX))*safe_adb(1.e-3_r8,VLSoilMicPMass_vr(L,NY,NX))
        this%h2D_MAOM_C_vr(ncol,L)          = (sum(SorbedOM_vr(idom_doc,:,L,NY,NX))+sum(SorbedOM_vr(idom_acetate,:,L,NY,NX)))*safe_adb(1.e-3_r8,VLSoilMicPMass_vr(L,NY,NX))
        this%h2D_microbC_vr(ncol,L)         =  micBE(ielmc)/DVOLL
        this%h2D_microbN_vr(ncol,L)         =  micBE(ielmn)/DVOLL
        this%h2D_microbP_vr(ncol,L)         =  micBE(ielmp)/DVOLL
        this%h2D_AeroBact_PrimS_lim_vr(ncol,L)=AeroBact_PrimeS_lim_vr(L,NY,NX)
        this%h2D_AeroFung_PrimS_lim_vr(ncol,L)=AeroFung_PrimeS_lim_vr(L,NY,NX)
        this%h2D_tSOCL_vr(ncol,L)           = SoilOrgM_vr(ielmc,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_tSON_vr(ncol,L)            = SoilOrgM_vr(ielmn,L,NY,NX)/DVOLL
        this%h2D_tSOP_vr(ncol,L)            = SoilOrgM_vr(ielmp,L,NY,NX)/DVOLL
        this%h2D_VHeatCap_vr(ncol,L)        = VHeatCapacity_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_Root_CO2_vr(ncol,L)        = trcg_root_vr(idg_CO2,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_NO3_vr(ncol,L)             = trc_solcl_vr(ids_NO3,L,NY,NX)
        this%h2D_NH4_vr(ncol,L)             = trc_solcl_vr(ids_NH4,L,NY,NX)
        this%h2D_Aqua_CO2_vr(ncol,L)        = trc_solcl_vr(idg_CO2,L,NY,NX)
        this%h2D_Aqua_CH4_vr(ncol,L)        = trc_solcl_vr(idg_CH4,L,NY,NX)
        this%h2D_Aqua_O2_vr(ncol,L)         = trc_solcl_vr(idg_O2,L,NY,NX)
        this%h2D_Aqua_N2O_vr(ncol,L)        = trc_solcl_vr(idg_N2O,L,NY,NX)
        this%h2D_Aqua_NH3_vr(ncol,L)        = trc_solcl_vr(idg_NH3,L,NY,NX)
        this%h2D_Aqua_N2_vr(ncol,L)         = trc_solcl_vr(idg_N2,L,NY,NX)
        this%h2D_Aqua_Ar_vr(ncol,L)         = trc_solcl_vr(idg_Ar,L,NY,NX)
        this%h2D_Aqua_H2_vr(ncol,L)         = trc_solcl_vr(idg_H2,L,NY,NX)

        this%h2D_TEMP_vr(ncol,L)            = TCS_vr(L,NY,NX)
        this%h2D_decomp_OStress_vr(ncol,L)  = OxyDecompLimiter_vr(L,NY,NX)
        this%h2D_RO2Decomp_vr(ncol,L)       = RO2DecompUptk_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_Decomp_temp_FN_vr(ncol,L)  = TempSensDecomp_vr(L,NY,NX)
        this%h2D_FracLitMix_vr(ncol,L)      = FracLitrMix_vr(L,NY,NX)
        this%h2D_Decomp_Moist_FN_vr(ncol,L) = MoistSensDecomp_vr(L,NY,NX)
        this%h2D_HeatFlow_vr(ncol,L)        = THeatFlowCellSoil_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_HeatUptk_vr(ncol,L)        = THeatLossRoot2Soil_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_VSPore_vr(ncol,L)          = POROS_vr(L,NY,NX)
        this%h2D_FLO_MICP_vr(ncol,L)        = m2mm*TWatFlowCellMicP_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_FLO_MACP_vr(ncol,L)        = m2mm*TWatFlowCellMacP_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_rVSM_vr    (ncol,L)        = ThetaH2OZ_vr(L,NY,NX)
        this%h2D_rVSICE_vr  (ncol,L)        = ThetaICEZ_vr(L,NY,NX)
        this%h2D_PSI_vr(ncol,L)             = PSISoilMatricP_vr(L,NY,NX)+PSISoilOsmotic_vr(L,NY,NX)
        this%h2D_PsiO_vr(ncol,L)            = PSISoilOsmotic_vr(L,NY,NX)
        this%h2D_RootH2OUP_vr(ncol,L)       = TWaterPlantRoot2SoilPrev_vr(L,NY,NX)
        this%h2D_cNH4t_vr(ncol,L)           = safe_adb(trcs_solml_vr(ids_NH4,L,NY,NX)+trcs_solml_vr(ids_NH4B,L,NY,NX) &
                                               +natomw*(trcx_solml_vr(idx_NH4,L,NY,NX)+trcx_solml_vr(idx_NH4B,L,NY,NX)),&
                                               VLSoilMicPMass_vr(L,NY,NX))
  
        this%h2D_cNO3t_vr(ncol,L)= safe_adb(trcs_solml_vr(ids_NO3,L,NY,NX)+trcs_solml_vr(ids_NO3B,L,NY,NX) &
                                               +trcs_solml_vr(ids_NO2,L,NY,NX)+trcs_solml_vr(ids_NO2B,L,NY,NX),&
                                               VLSoilMicPMass_vr(L,NY,NX))

        this%h2D_cPO4_vr(ncol,L) = safe_adb(trcs_solml_vr(ids_H1PO4,L,NY,NX)+trcs_solml_vr(ids_H1PO4B,L,NY,NX) &
                                               +trcs_solml_vr(ids_H2PO4,L,NY,NX)+trcs_solml_vr(ids_H2PO4B,L,NY,NX),&
                                               VLWatMicP_vr(L,NY,NX))
        this%h2D_cEXCH_P_vr(ncol,L)= patomw*safe_adb(trcx_solml_vr(idx_HPO4,L,NY,NX)+trcx_solml_vr(idx_H2PO4,L,NY,NX) &
                                               +trcx_solml_vr(idx_HPO4B,L,NY,NX)+trcx_solml_vr(idx_H2PO4B,L,NY,NX),&
                                               VLSoilMicPMass_vr(L,NY,NX))
        this%h2D_microb_N2fix_vr(ncol,L)         = 1.e6_r8*Micb_N2Fixation_vr(L,NY,NX)/DVOLL
        this%h1D_FreeNFix_col(ncol)              = this%h1D_FreeNFix_col(ncol)+ Micb_N2Fixation_vr(L,NY,NX)
        this%h2D_ElectricConductivity_vr(ncol,L) = ElectricConductivity_vr(L,NY,NX)
        
        this%h2D_HydCondSoil_vr(ncol,L) = HydCondSoil_3D(3,L,NY,NX)

        !aerobic heterotropic bacteria
        call SumMicbGroup(L,NY,NX,micpar%mid_HeterAerobBacter,MicbE)
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

          this%h2D_fermentor_frac_vr(ncol,L)=safe_adb(this%h2D_fermentorC_vr(ncol,L),this%h2D_microbC_vr(ncol,L))
        
        !acetogenic methanogen
        call SumMicbGroup(L,NY,NX,micpar%mid_HeterAcetoCH4GenArchea,MicbE)
        this%h2D_acetometgC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_acetometgN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_acetometgP_vr(ncol,L) = micBE(ielmp)/DVOLL
        
        this%h2D_acetometh_frac_vr(ncol,L)=safe_adb(this%h2D_acetometgC_vr(ncol,L),this%h2D_microbC_vr(ncol,L))

        !aerobic N2 fixer
        call SumMicbGroup(L,NY,NX,micpar%mid_HeterAerobN2Fixer,MicbE)
        this%h2D_aeroN2fixC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_aeroN2fixN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_aeroN2fixP_vr(ncol,L) = micBE(ielmp)/DVOLL

        this%h2D_RDECOMPC_SOM_vr(ncol,L)     = tRHydlySOM_vr(ielmc,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RDECOMPC_BReSOM_vr(ncol,L)  = tRHydlyBioReSOM_vr(ielmc,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RDECOMPC_SorpSOM_vr(ncol,L) = tRHydlySoprtOM_vr(ielmc,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_MicrobAct_vr(ncol,L)        = TMicHeterActivity_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

        DO jj=1,jcplx
          this%h3D_MicrobActCps_vr(ncol,L,jj) = ROQC4HeterMicActCmpK_vr(jj,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

          this%h3D_SOMHydrylScalCps_vr(ncol,L,jj) = AZMAX1(RHydrolysisScalCmpK_vr(jj,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX))
        
          this%h3D_HydrolCSOMCps_vr(ncol,L,jj) = RHydlySOCK_vr(jj,L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        ENDDO
        !anaerobic N2 fixer
        call SumMicbGroup(L,NY,NX,micpar%mid_HeterAnaerobN2Fixer,MicbE)
        this%h2D_anaeN2FixC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_anaeN2FixN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_anaeN2FixP_vr(ncol,L) = micBE(ielmp)/DVOLL

        call SumMicbGroup(L,NY,NX,micpar%mid_AutoAmmoniaOxidBacter,MicbE,isauto=.true.)
        this%h2D_NH3OxiBactC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_NH3OxiBactN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_NH3OxiBactP_vr(ncol,L) = micBE(ielmp)/DVOLL

        call SumMicbGroup(L,NY,NX,micpar%mid_AutoNitriteOxidBacter,MicbE,isauto=.true.)
        this%h2D_NO2OxiBactC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_NO2OxiBactN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_NO2OxiBactP_vr(ncol,L) = micBE(ielmp)/DVOLL

        call SumMicbGroup(L,NY,NX,micpar%mid_AutoAeroCH4OxiBacter,MicbE,isauto=.true.)
        this%h2D_CH4AeroOxiC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_CH4AeroOxiN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_CH4AeroOxiP_vr(ncol,L) = micBE(ielmp)/DVOLL
        
        this%h2D_tRespGrossHeterUlm_vr(ncol,L)= tRespGrossHeterUlm_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_tRespGrossHeter_vr(ncol,L)= tRespGrossHeter_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

        call SumMicbGroup(L,NY,NX,micpar%mid_AutoH2GenoCH4GenArchea,MicbE,isauto=.true.)
        this%h2D_H2MethogenC_vr(ncol,L) = micBE(ielmc)/DVOLL
        this%h2D_H2MethogenN_vr(ncol,L) = micBE(ielmn)/DVOLL
        this%h2D_H2MethogenP_vr(ncol,L) = micBE(ielmp)/DVOLL        
        
        this%h2D_hydrogMeth_frac_vr(ncol,L)=safe_adb(this%h2D_H2MethogenC_vr(ncol,L),this%h2D_microbC_vr(ncol,L))

        this%h2D_TSolidOMActC_vr(ncol,L)     = TSolidOMActC_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_TSolidOMActCDens_vr(ncol,L) = safe_adb(TSolidOMActC_vr(L,NY,NX),TSolidOMC_vr(L,NY,NX))
        this%h2D_tOMActCDens_vr(ncol,L)      = safe_adb(tOMActC_vr(L,NY,NX),(tOMActC_vr(L,NY,NX)+TSolidOMC_vr(L,NY,NX)))
        this%h2D_RCH4ProdAcetcl_vr(ncol,L)   = RCH4ProdAcetcl_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RCH4ProdHydrog_vr(ncol,L)   = RCH4ProdHydrog_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)        

        this%h2D_RCH4Oxi_aero_vr(ncol,L)     = RCH4Oxi_aero_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RCH4Oxi_aero_col(ncol) = this%h1D_RCH4Oxi_aero_col(ncol) + RCH4Oxi_aero_vr(L,NY,NX)
        this%h2D_RCH4Oxi_anmo_vr(ncol,L)     = RCH4Oxi_anmo_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RCH4Oxi_anmo_col(ncol) = this%h1D_RCH4Oxi_anmo_col(ncol) + RCH4Oxi_anmo_vr(L,NY,NX)

        this%h2D_RFerment_vr(ncol,L)         = RFerment_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_nh3oxi_vr(ncol,L)           = RNH3oxi_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_n2oprod_vr(ncol,L)          = (RN2ODeniProd_vr(L,NY,NX)+RN2ONitProd_vr(L,NY,NX) &
                               +RN2OChemoProd_vr(L,NY,NX)-RN2ORedux_vr(L,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RootAR_vr(ncol,L) = -RootCO2Autor_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RootAR2soil_vr(ncol,L)=-RootCO2Ar2Soil_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h2D_RootAR2Root_vr(ncol,L)=-RootCO2Ar2Root_vr(L,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        if(plant_model)then
          this%h2D_RootMassC_vr(ncol,L)     = RootMycoMassElm_vr(ielmc,ipltroot,L,NY,NX)/DVOLL
          this%h2D_RootMassN_vr(ncol,L)     = RootMycoMassElm_vr(ielmn,ipltroot,L,NY,NX)/DVOLL
          this%h2D_RootMassP_vr(ncol,L)     = RootMycoMassElm_vr(ielmp,ipltroot,L,NY,NX)/DVOLL                
        endif
      ENDDO
      this%h1D_RCH4Oxi_aero_col(ncol) = this%h1D_RCH4Oxi_aero_col(ncol)/AREA_3D(3,NU_col(NY,NX),NY,NX)      
      this%h1D_RCH4Oxi_anmo_col(ncol) = this%h1D_RCH4Oxi_anmo_col(ncol)/AREA_3D(3,NU_col(NY,NX),NY,NX)      

      this%h1D_RootAR_col(ncol)  = -AZERO(RootCO2Autor_col(NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)      
      this%h1D_RootCO2Relez_col(ncol)=RootCO2Emis2Root_col(NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)      
      this%h1d_fPAR_col(ncol) = 0._r8
      DO NZ=1,NP0_col(NY,NX)
        nptc=get_pft(NZ,NY,NX)
        this%h1D_ROOT_NONSTC_ptc(nptc)  = RootMycoNonstElms_pft(ielmc,ipltroot,NZ,NY,NX)
        this%h1D_ROOT_NONSTN_ptc(nptc)  = RootMycoNonstElms_pft(ielmn,ipltroot,NZ,NY,NX)
        this%h1D_ROOT_NONSTP_ptc(nptc)  = RootMycoNonstElms_pft(ielmp,ipltroot,NZ,NY,NX)
        this%h1D_SHOOT_NONSTC_ptc(nptc) = CanopyNonstElms_pft(ielmc,NZ,NY,NX)
        this%h1D_SHOOT_NONSTN_ptc(nptc) = CanopyNonstElms_pft(ielmn,NZ,NY,NX)
        this%h1D_SHOOT_NONSTP_ptc(nptc) = CanopyNonstElms_pft(ielmp,NZ,NY,NX)

        if(CO2FixLL_pft(NZ,NY,NX)/=spval .and. CO2FixCL_pft(NZ,NY,NX)/=spval)then
          this%h1D_dCAN_GPP_CLIM_ptc(nptc) = AZERO(CO2FixCL_pft(NZ,NY,NX)-GrossCO2Fix_pft(NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h1D_dCAN_GPP_eLIM_ptc(nptc) = AZERO(CO2FixLL_pft(NZ,NY,NX)-GrossCO2Fix_pft(NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        endif

        this%h1D_LeafC3ChlCperm2LA_ptc(nptc)   = LeafC3ChlCperm2LA_pft(NZ,NY,NX)*1.e3_r8
        this%h1D_LeafC4ChlCperm2LA_ptc(nptc)   = LeafC4ChlCperm2LA_pft(NZ,NY,NX)*1.e3_r8
        this%h1D_LeafRubiscoCperm2LA_ptc(nptc) = LeafRubiscoCperm2LA_pft(NZ,NY,NX)
        this%h1D_LeafPEPCperm2LA_ptc(nptc)     = LeafPEPCperm2LA_pft(NZ,NY,NX)

        this%h1D_MIN_LWP_ptc(nptc)      = PSICanPDailyMin_pft(NZ,NY,NX)
        this%h1D_SLA_ptc(nptc)          = 1.e4_r8*SpecificLeafArea_pft(NZ,NY,NX)        
        this%h1D_LEAF_PC_ptc(nptc)       = safe_adb(LeafStrutElms_pft(ielmp,NZ,NY,NX)+CanopyNonstElms_pft(ielmp,NZ,NY,NX), &
                                                 LeafStrutElms_pft(ielmc,NZ,NY,NX)+CanopyNonstElms_pft(ielmc,NZ,NY,NX))
        this%h1D_CAN_RN_ptc(nptc)        = MJ2W*RadNet2Canopy_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_LE_ptc(nptc)        = MJ2W*EvapTransLHeat_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_H_ptc(nptc)         = MJ2W*HeatXAir2PCan_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_G_ptc(nptc)         = MJ2W*HeatStorCanopy_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_TEMPC_ptc(nptc)     = TdegCCanopy_pft(NZ,NY,NX)
        this%h1D_CAN_TEMPFN_ptc(nptc)    = fTCanopyGroth_pft(NZ,NY,NX)
        this%h1D_CAN_CO2_FLX_ptc(nptc)   = CO2NetFix_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*GramPerHr2umolPerSec(idg_CO2)
        this%h1D_CAN_GPP_ptc(nptc)       = GrossCO2Fix_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_cumGPP_ptc(nptc)       = GrossCO2Fix_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_RA_ptc(nptc)        = CanopyGrosRCO2_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_GROWTH_ptc(nptc)    = canopy_growth_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_cTNC_ptc(nptc)          = CanopyNonstElmConc_pft(ielmc,NZ,NY,NX)
        this%h1D_cTNN_ptc(nptc)          = CanopyNonstElmConc_pft(ielmn,NZ,NY,NX)
        this%h1D_cTNP_ptc(nptc)          = CanopyNonstElmConc_pft(ielmp,NZ,NY,NX)
        this%h1D_CanNonstBconc_ptc(nptc) = sum(CanopyNonstElmConc_pft(1:NumPlantChemElms,NZ,NY,NX))
        this%h1D_STOML_RSC_CO2_ptc(nptc) = CanPStomaResistH2O_pft(NZ,NY,NX)*1.56_r8*secs1hour
        this%h1D_STOML_Min_RSC_CO2_ptc(nptc)=CanopyMinStomaResistH2O_pft(NZ,NY,NX)*1.56_r8*secs1hour
        this%h1D_Km_CO2_carboxy_ptc(nptc)= Km4RubiscoCarboxy_pft(NZ,NY,NX)
        this%h1D_Ci_mesophyll_ptc(nptc)  = LeafIntracellularCO2_pft(NZ,NY,NX)
        this%h1D_BLYR_RSC_CO2_ptc(nptc)  = CanopyBndlResist_pft(NZ,NY,NX)*1.34_r8*secs1hour
        this%h1D_CAN_CO2_ptc(nptc)       = CanopyGasCO2_pft(NZ,NY,NX)
        this%h1D_O2L_ptc(nptc)           = O2L_pft(NZ,NY,NX)
        this%h1D_LAI_ptc(nptc)           = LeafStalkArea_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_PSI_CAN_ptc(nptc)       = PSICanopy_pft(NZ,NY,NX)
        this%h1D_TURG_CAN_ptc(nptc)      = PSICanopyTurg_pft(NZ,NY,NX)
        this%h1D_STOML_RSC_H2O_ptc(nptc)  = CanPStomaResistH2O_pft(NZ,NY,NX)*secs1hour
        this%h1D_BLYR_RSC_H2O_ptc(nptc)  = CanopyBndlResist_pft(NZ,NY,NX)*secs1hour
        this%h1D_TRANSPN_ptc(nptc)       = Transpiration_pft(NZ,NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_QTRANSP_col(ncol)       = this%h1D_QTRANSP_col(ncol)+this%h1D_TRANSPN_ptc(nptc)
        this%h1D_NH4_UPTK_FLX_ptc(nptc)  = RootNH4Uptake_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_NO3_UPTK_FLX_ptc(nptc)  = RootNO3Uptake_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)        
        this%h1D_N2_FIXN_FLX_ptc(nptc)   = RootN2Fix_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_cNH3_FLX_ptc(nptc)      = NH3Dep2Can_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_TC_Groth_ptc(nptc)      = TCGroth_pft(NZ,NY,NX)
        this%h1D_PO4_UPTK_FLX_ptc(nptc)  = RootH2PO4Uptake_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_frcPARabs_ptc(nptc)     = FracPARads2Canopy_pft(NZ,NY,NX)
        this%h1D_PAR_CAN_ptc(nptc)       = RadPARbyCanopy_pft(NZ,NY,NX)   !umol /m2/s      
        this%h1d_fPAR_col(ncol)          = this%h1d_fPAR_col(ncol)+RadPARbyCanopy_pft(NZ,NY,NX)
        this%h1D_SHOOT_C_ptc(nptc)       = ShootElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Plant_C_ptc(nptc)       = (ShootElms_pft(ielmc,NZ,NY,NX) &
          +RootElms_pft(ielmc,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LEAF_C_ptc(nptc)        = LeafStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Petole_C_ptc(nptc)      = PetoleStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STALK_C_ptc(nptc)       = StalkStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RESERVE_C_ptc(nptc)     = StalkRsrvElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HUSK_C_ptc(nptc)        = (HuskStrutElms_pft(ielmc,NZ,NY,NX) &
          +EarStrutElms_pft(ielmc,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_GRAIN_C_ptc(nptc)       = GrainStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOT_C_ptc(nptc)        = RootElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOTST_C_ptc(nptc)      = RootStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOTST_N_ptc(nptc)      = RootStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOTST_P_ptc(nptc)      = RootStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)        
        this%h1D_RootNodule_C_ptc(nptc)  = RootNoduleElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ShootNodule_C_ptc(nptc)  = ShootNoduleElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ShootNodule_N_ptc(nptc)  = ShootNoduleElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ShootNodule_P_ptc(nptc)  = ShootNoduleElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

        this%h1D_STORED_C_ptc(nptc)      = SeasonalNonstElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_GRAIN_NO_ptc(nptc)      = CanopySeedNum_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LAIb_ptc(nptc)          = CanopyLeafArea_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_EXUD_CumYr_C_FLX_ptc(nptc) = PlantExudElm_CumYr_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LITRf_C_FLX_ptc(nptc)      = LitrfalStrutElms_CumYr_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_SURF_LITRf_C_FLX_ptc(nptc) = SurfLitrfalStrutElms_CumYr_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_AUTO_RESP_FLX_ptc(nptc)    = GrossResp_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HVST_C_FLX_ptc(nptc)       = EcoHavstElmnt_CumYr_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HVST_N_FLX_ptc(nptc)       = EcoHavstElmnt_CumYr_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HVST_P_FLX_ptc(nptc)       = EcoHavstElmnt_CumYr_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_PLANT_BALANCE_C_ptc(nptc)  = PlantElmBalCum_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_PLANT_BALANCE_N_ptc(nptc)  = PlantElmBalCum_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_PLANT_BALANCE_P_ptc(nptc)  = PlantElmBalCum_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STANDING_DEAD_C_ptc(nptc)  = StandDeadStrutElms_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_FIREp_CO2_FLX_ptc(nptc)    = CO2ByFire_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_FIREp_CH4_FLX_ptc(nptc)    = CH4ByFire_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_cNPP_ptc(nptc)              = cumNPP_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CAN_HT_ptc(nptc)           = CanopyHeight_pft(NZ,NY,NX)
        this%h1D_POPN_ptc(nptc)             = PlantPopulation_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_tTRANSPN_ptc(nptc)         = -ETCanopy_CumYr_pft(NZ,NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_WTR_STRESS_ptc(nptc)       = HoursTooLowPsiCan_pft(NZ,NY,NX)
        this%h1D_LeafProteinCperm2_ptc(nptc)=LeafProteinCperm2LA_pft(NZ,NY,NX)
        this%h1D_VcMaxRubisco_ptc(nptc) = CanopyVcMaxRubisco25C_pft(NZ,NY,NX)
        this%h1D_VoMaxRubisco_ptc(nptc) = CanopyVoMaxRubisco25C_pft(NZ,NY,NX)
        this%h1D_VcMaxPEP_ptc(nptc)     = CanopyVcMaxPEP25C_pft(NZ,NY,NX)
        this%h1D_JMaxPhoto_ptc(nptc)    = ElectronTransptJmax25C_pft(NZ,NY,NX)
        this%h1D_TFN_Carboxy_ptc(nptc)  = TFN_Carboxy_pft(NZ,NY,NX)
        this%h1D_TFN_Oxygen_ptc(nptc)   = TFN_Oxygen_pft(NZ,NY,NX)
        this%h1D_TFN_eTranspt_ptc(nptc) = TFN_eTranspt_pft(NZ,NY,NX)
        this%h1D_PARSunlit_ptc(nptc)    = PARSunlit_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_PARSunsha_ptc(nptc)    = PARSunsha_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CH2OSunlit_ptc(nptc)   = CH2OSunlit_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_CH2OSunsha_ptc(nptc)   = CH2OSunsha_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)

        this%h1D_fClump_ptc(nptc)  = ClumpFactorNow_pft(NZ,NY,NX)
        this%h1D_LeafAreaSunlit_ptc(nptc)=LeafAreaSunlit_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_OXY_STRESS_ptc(nptc)   = PlantO2Stress_pft(NZ,NY,NX)
        this%h1D_SHOOT_N_ptc(nptc)      = ShootElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)        
        this%h1D_Plant_N_ptc(nptc)      = (ShootElms_pft(ielmn,NZ,NY,NX)&
          +RootElms_pft(ielmn,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_fCNLFW_ptc(nptc) = safe_adb(1._r8,fNCLFW_pft(NZ,NY,NX))  
        this%h1D_fCPLFW_ptc(nptc) = safe_adb(1._r8,fPCLFW_pft(NZ,NY,NX))
        this%h1D_LEAF_N_ptc(nptc)    = LeafStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LEAFN2LAI_ptc(nptc) = safe_adb(LeafStrutElms_pft(ielmn,NZ,NY,NX),CanopyLeafArea_pft(NZ,NY,NX))
        this%h1D_Petole_N_ptc(nptc)  = PetoleStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STALK_N_ptc(nptc)   = StalkStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RESERVE_N_ptc(nptc) = StalkRsrvElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HUSK_N_ptc(nptc)    = (HuskStrutElms_pft(ielmn,NZ,NY,NX) &
          +EarStrutElms_pft(ielmn,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_GRAIN_N_ptc(nptc)          = GrainStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOT_N_ptc(nptc)           = RootElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RootNodule_N_ptc(nptc)     = RootNoduleElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STORED_N_ptc(nptc)         = SeasonalNonstElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_EXUD_N_FLX_ptc(nptc)       = PlantExudElm_CumYr_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Uptk_N_Flx_ptc(nptc)       = RootUptk_N_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Uptk_P_Flx_ptc(nptc)       = RootUptk_P_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LITRf_N_FLX_ptc(nptc)      = LitrfalStrutElms_CumYr_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_cum_N_FIXED_ptc(nptc)      = PlantN2Fix_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_TreeRingRadius_ptc(nptc)   = TreeRingAveRadius_pft(NZ,NY,NX)
        this%h1D_NH3can_FLX_ptc(nptc)       = NH3Emis_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STANDING_DEAD_N_ptc(nptc)  = StandDeadStrutElms_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_FIREp_N_FLX_ptc(nptc)      = NH3byFire_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_SURF_LITRf_N_FLX_ptc(nptc) = SurfLitrfalStrutElms_CumYr_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_SHOOT_P_ptc(nptc)          = ShootElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Plant_P_ptc(nptc)          = (ShootElms_pft(ielmp,NZ,NY,NX) &
          +RootElms_pft(ielmp,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)        
        this%h1D_stomatal_stress_ptc(nptc) = StomatalStress_pft(NZ,NY,NX)
        this%h1D_CANDew_ptc(nptc)          = QdewCanopy_CumYr_pft(NZ,NY,NX)*m2mm/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LEAF_P_ptc(nptc)          = LeafStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_Petole_P_ptc(nptc)        = PetoleStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STALK_P_ptc(nptc)         = StalkStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RESERVE_P_ptc(nptc)       = StalkRsrvElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_HUSK_P_ptc(nptc)          = (HuskStrutElms_pft(ielmp,NZ,NY,NX) &
          +EarStrutElms_pft(ielmp,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_GRAIN_P_ptc(nptc)          = GrainStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ROOT_P_ptc(nptc)           = RootElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_RootNodule_P_ptc(nptc)     = RootNoduleElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STORED_P_ptc(nptc)         = SeasonalNonstElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_EXUD_P_FLX_ptc(nptc)       = PlantExudElm_CumYr_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_LITRf_P_FLX_ptc(nptc)      = LitrfalStrutElms_CumYr_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_STANDING_DEAD_P_ptc(nptc)  = StandDeadStrutElms_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_FIREp_P_FLX_ptc(nptc)      = PO4byFire_CumYr_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_SURF_LITRf_P_FLX_ptc(nptc) = SurfLitrfalStrutElms_CumYr_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ShootRootXferC_ptc(nptc)   = ShootRootXferElm_pft(ielmc,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ShootRootXferN_ptc(nptc)   = ShootRootXferElm_pft(ielmn,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_ShootRootXferP_ptc(nptc)   = ShootRootXferElm_pft(ielmp,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        this%h1D_BRANCH_NO_ptc(nptc)        = NumOfBranches_pft(NZ,NY,NX)
        this%h1D_MainBranchNO_ptc(nptc)     = MainBranchNum_pft(NZ,NY,NX)
        this%h1D_RCanMaintDef_CO2_pft(nptc) = RCanMaintDef_CO2_pft(NZ,NY,NX)
        this%h1D_LEAF_NC_ptc(nptc)      = safe_adb(LeafStrutElms_pft(ielmn,NZ,NY,NX)+CanopyNonstElms_pft(ielmn,NZ,NY,NX),&
                                                 LeafStrutElms_pft(ielmc,NZ,NY,NX)+CanopyNonstElms_pft(ielmc,NZ,NY,NX))
        this%h1D_RootMaintDef_CO2_pft(nptc) = sum(RootMaintDef_CO2_pvr(ipltroot,1:JZ,NZ,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)

        IF(NumOfBranches_pft(NZ,NY,NX)>0)then
          DO K=1,MaxNodesPerBranch
            this%h2D_ProteinNperm2LeafArea_pnd(nptc,K)=0._r8
            DO NB=1,NumOfBranches_pft(NZ,NY,NX)
              this%h2D_ProteinNperm2LeafArea_pnd(nptc,K)=this%h2D_ProteinNperm2LeafArea_pnd(nptc,K)+ProteinCperm2LeafArea_node(K,NB,NZ,NY,NX)
            ENDDO
            !assuming protein C to N mass ratio is 3.3 (median value)
            this%h2D_ProteinNperm2LeafArea_pnd(nptc,K)=this%h2D_ProteinNperm2LeafArea_pnd(nptc,K)/(NumOfBranches_pft(NZ,NY,NX)*3.3_r8)
          ENDDO  
        ENDIF
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
          this%h1D_RUB_ACTVN_ptc(nptc)=0._r8; NB1=0
          this%h1D_CanopyNLim_ptc(nptc)=0._r8
          this%h1D_CanopyPLim_ptc(nptc)=0._r8
          DO NB=1,NumOfBranches_pft(NZ,NY,NX)
            this%h2D_LEAF_NODE_NO_ptc(nptc,NB)              = NumOfLeaves_brch(NB,NZ,NY,NX)
            if(RubiscoActivity_brch(NB,NZ,NY,NX)>0._r8)then
              NB1=NB1+1
              this%h1D_RUB_ACTVN_ptc(nptc)   = this%h1D_RUB_ACTVN_ptc(nptc)+ RubiscoActivity_brch(NB,NZ,NY,NX)
              this%h1D_CanopyNLim_ptc(nptc)  = this%h1D_CanopyNLim_ptc(nptc)+ CanopyNLimFactor_brch(NB,NZ,NY,NX)
              this%h1D_CanopyPLim_ptc(nptc)  = this%h1D_CanopyPLim_ptc(nptc)+ CanopyPLimFactor_brch(NB,NZ,NY,NX)      
            endif
            this%h3D_PARTS_ptc(nptc,1:NumOfPlantMorphUnits,NB) = PARTS_brch(1:NumOfPlantMorphUnits,NB,NZ,NY,NX)            
          ENDDO
          if(NB1>0)then
            this%h1D_RUB_ACTVN_ptc(nptc)=this%h1D_RUB_ACTVN_ptc(nptc)/real(NB1,kind=r8)
            this%h1D_CanopyNLim_ptc(nptc)=this%h1D_CanopyNLim_ptc(nptc)/real(NB1,kind=r8)
            this%h1D_CanopyPLim_ptc(nptc)=this%h1D_CanopyPLim_ptc(nptc)/real(NB1,kind=r8)            
          endif
        endif
        DO L=1,NumCanopyLayers
          this%h2D_CanopyLAIZ_plyr(nptc,L)=CanopyLeafAreaZ_pft(L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        ENDDO  
        this%h1D_RootAR_ptc(nptc)          = 0._r8
        this%h1D_RootLenPerPlant_ptc(nptc) = 0._r8
        DO L=1,JZ
          this%h1D_RootAR_ptc(nptc)=this%h1D_RootAR_ptc(nptc)-RootCO2Autor_pvr(ipltroot,L,NZ,NY,NX)
          DVOLL                                  = DLYR_3D(3,L,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
          if(DVOLL<1.e-8_r8)cycle        
          this%h2D_RootRadialKond2H2O_pvr(nptc,L)=RootRadialKond2H2O_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h2D_RootAXialKond2H2O_pvr(nptc,L) =RootAxialKond2H2O_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h2D_VmaxNH4Root_pvr(nptc,L)       = VmaxNH4Root_pvr(ipltroot,L,NZ,NY,NX)*1.e6/natomw
          this%h2D_VmaxNO3Root_pvr(nptc,L)       = VmaxNO3Root_pvr(ipltroot,L,NZ,NY,NX)*1.e6/natomw
          this%h2D_RootMassC_pvr(nptc,L)         = RootMycoMassElm_pvr(ielmc,ipltroot,L,NZ,NY,NX)/DVOLL
          this%h2D_RootNutupk_fClim_pvr(nptc,L)  = Nutruptk_fClim_rpvr(ipltroot,L,NZ,NY,NX)
          this%h2D_RootNutupk_fNlim_pvr(nptc,L)  = Nutruptk_fNlim_rpvr(ipltroot,L,NZ,NY,NX)
          this%h2D_RootNutupk_fPlim_pvr(nptc,L)  = Nutruptk_fPlim_rpvr(ipltroot,L,NZ,NY,NX)
          this%h2D_RootNutupk_fProtC_pvr(nptc,L) = Nutruptk_fProtC_rpvr(ipltroot,L,NZ,NY,NX)
          this%h2D_O2_rootconduct_pvr(nptc,L)    = RootGasConductance_rpvr(idg_O2,ipltroot,L,NZ,NY,NX)
          this%h2D_CO2_rootconduct_pvr(nptc,L)   = RootGasConductance_rpvr(idg_CO2,ipltroot,L,NZ,NY,NX)
          this%h2D_fTRootGro_pvr(nptc,L)         = fTgrowRootP_vr(L,NZ,NY,NX)
          this%h2D_fRootGrowPSISense_pvr(nptc,L) = fRootGrowPSISense_pvr(ipltroot,L,NZ,NY,NX)

          this%h2D_RootSurfAreaPP_pvr(nptc,L) = RootAreaPerPlant_pvr(ipltroot,L,NZ,NY,NX)  
          this%h2D_ROOT_OSTRESS_pvr(nptc,L) = RAutoRootO2Limter_rpvr(ipltroot,L,NZ,NY,NX)
          this%h2D_PSI_RT_pvr(nptc,L)       = PSIRoot_pvr(ipltroot,L,NZ,NY,NX)
          this%h2D_RootH2OUptkStress_pvr(nptc,L)=RootH2OUptkStress_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h2D_RootH2OUptk_pvr(nptc,L) = 1.e3*RPlantRootH2OUptk_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h2D_RootMaintDef_CO2_pvr(nptc,L)=RootMaintDef_CO2_pvr(ipltroot,L,NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)
          this%h2D_prtUP_NH4_pvr(nptc,L)    = (sum(RootNutUptake_pvr(ids_NH4,:,L,NZ,NY,NX))+&
            sum(RootNutUptake_pvr(ids_NH4B,:,L,NZ,NY,NX)))/AREA_3D(3,L,NY,NX)
          this%h2D_prtUP_NO3_pvr(nptc,L)  = (sum(RootNutUptake_pvr(ids_NO3,:,L,NZ,NY,NX))+&
            sum(RootNutUptake_pvr(ids_NO3B,:,L,NZ,NY,NX)))/AREA_3D(3,L,NY,NX)
          this%h2D_prtUP_PO4_pvr(nptc,L)  = (sum(RootNutUptake_pvr(ids_H2PO4,:,L,NZ,NY,NX))+&
            sum(RootNutUptake_pvr(ids_H2PO4B,:,L,NZ,NY,NX)))/AREA_3D(3,L,NY,NX)
          this%h2D_DNS_RT_pvr(nptc,L) = RootLenDensPerPlant_pvr(ipltroot,L,NZ,NY,NX)* &
            PlantPopulation_pft(NZ,NY,NX)/AREA_3D(3,NU_col(NY,NX),NY,NX)*1.e-4_r8
          this%h1D_RootLenPerPlant_ptc(nptc)=this%h1D_RootLenPerPlant_ptc(nptc)+RootLenPerPlant_pvr(ipltroot,L,NZ,NY,NX)

          this%h2D_ROOTNLim_rpvr(nptc,L) = ROOTNLim_rpvr(ipltroot,L,NZ,NY,NX)
          this%h2D_ROOTPLim_rpvr(nptc,L) = ROOTPLim_rpvr(ipltroot,L,NZ,NY,NX)
          this%h2D_Root1stStrutC_pvr(nptc,L) = 0._r8
          this%h2D_Root1stStrutN_pvr(nptc,L) = 0._r8
          this%h2D_Root1stStrutP_pvr(nptc,L) = 0._r8
          this%h2D_Root2ndStrutC_pvr(nptc,L) = 0._r8
          this%h2D_Root2ndStrutN_pvr(nptc,L) = 0._r8
          this%h2D_Root2ndStrutP_pvr(nptc,L) = 0._r8
          this%h2D_RootNonstBConc_pvr(nptc,L)=sum(RootNonstructElmConc_rpvr(1:NumPlantChemElms,ipltroot,L,NZ,NY,NX))
          this%h2D_Root1stAxesNumL_pvr(nptc,L)= Root1stXNumL_rpvr(ipltroot,L,NZ,NY,NX)
          this%h2D_Root2ndAxesNumL_pvr(nptc,L)= Root2ndXNumL_rpvr(ipltroot,L,NZ,NY,NX)
          this%h2D_RootKond2H2O_pvr(nptc,L)= safe_adb(1._r8,RootResist4H2O_pvr(ipltroot,L,NZ,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX))*1.e7/3600._r8
          DO NR=1,NumPrimeRootAxes_pft(NZ,NY,NX)
            this%h2D_Root1stStrutC_pvr(nptc,L)= this%h2D_Root1stStrutC_pvr(nptc,L) + &
              RootMyco1stStrutElms_rpvr(ielmc,ipltroot,L,NR,NZ,NY,NX)
            this%h2D_Root1stStrutN_pvr(nptc,L)= this%h2D_Root1stStrutN_pvr(nptc,L) + &
              RootMyco1stStrutElms_rpvr(ielmn,ipltroot,L,NR,NZ,NY,NX)
            this%h2D_Root1stStrutP_pvr(nptc,L)= this%h2D_Root1stStrutP_pvr(nptc,L) + &
              RootMyco1stStrutElms_rpvr(ielmp,ipltroot,L,NR,NZ,NY,NX)
            this%h2D_Root2ndStrutC_pvr(nptc,L)=this%h2D_Root2ndStrutC_pvr(nptc,L) + &
             RootMyco2ndStrutElms_rpvr(ielmc,ipltroot,L,NR,NZ,NY,NX)
            this%h2D_Root2ndStrutN_pvr(nptc,L)=this%h2D_Root2ndStrutN_pvr(nptc,L) + &
             RootMyco2ndStrutElms_rpvr(ielmn,ipltroot,L,NR,NZ,NY,NX)
            this%h2D_Root2ndStrutP_pvr(nptc,L)=this%h2D_Root2ndStrutP_pvr(nptc,L) + &
             RootMyco2ndStrutElms_rpvr(ielmp,ipltroot,L,NR,NZ,NY,NX)
          ENDDO
          this%h2D_Root1stStrutC_pvr(nptc,L) = this%h2D_Root1stStrutC_pvr(nptc,L)/DVOLL
          this%h2D_Root1stStrutN_pvr(nptc,L) = this%h2D_Root1stStrutN_pvr(nptc,L)/DVOLL
          this%h2D_Root1stStrutP_pvr(nptc,L) = this%h2D_Root1stStrutP_pvr(nptc,L)/DVOLL
          this%h2D_Root2ndStrutC_pvr(nptc,L) = this%h2D_Root2ndStrutC_pvr(nptc,L)/DVOLL
          this%h2D_Root2ndStrutN_pvr(nptc,L) = this%h2D_Root2ndStrutN_pvr(nptc,L)/DVOLL
          this%h2D_Root2ndStrutP_pvr(nptc,L) = this%h2D_Root2ndStrutP_pvr(nptc,L)/DVOLL
        ENDDO        
        this%h1D_RootAR_ptc(nptc)=this%h1D_RootAR_ptc(nptc)/AREA_3D(3,NU_col(NY,NX),NY,NX)
        
      ENDDO 
      this%h1d_fPAR_col(ncol)=safe_adb(this%h1d_fPAR_col(ncol),RadPARSolarBeam_col(NY,NX))
    ENDDO 
  ENDDO    
  call PrintInfo('end '//subname)
  end subroutine hist_update
  
end module HistDataType
