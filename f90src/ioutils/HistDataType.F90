module HistDataType
!
! this module is an intermediate step to support ascii output
! when output is done with netcdf, no id is needed.
  use data_kind_mod , only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval  => DAT_CONST_SPVAL
  use GridConsts
  use GridMod
  use HistFileMod
  use MiniMathMod, only : safe_adb,AZMAX1
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
  use PlantMngmtDataType
  use EcosimBGCFluxType
  use SoilPropertyDataType
  use SoilBGCDataType
  use AqueChemDatatype
  use SurfSoilDataType
implicit none

  type, public :: histdata_type
  real(r8),pointer   :: histr_1D_tFIRE_CO2_col(:)     !CO2byFire_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tFIRE_CH4_col(:)     !CH4byFire_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_cNH4_LITR_col(:)     !(trc_solml_vr(ids_NH4,0,NY,NX)+14.0*trcx_solml(idx_NH4,0,NY,NX))/SoilMicPMassLayer(0,NY,NX)
  real(r8),pointer   :: histr_1D_cNO3_LITR_col(:)      !(trc_solml_vr(ids_NO3,0,NY,NX)+trc_solml_vr(ids_NO2,0,NY,NX))/SoilMicPMassLayer(0,NY,NX)                            
  real(r8),pointer   :: histr_1D_ECO_HVST_N_col(:)     !EcoHavstElmnt_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_NET_N_MIN_col(:)      !-NetNH4Mineralize_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_SURF_tLITR_P_FLX_col(:)    !URSDP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_HUMUS_C_col(:)        !UORGC(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_HUMUS_N_col(:)        !UORGN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_HUMUS_P_col(:)        !UORGP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_AMENDED_P_col(:)       !FerPFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tLITRf_P_FLX_col(:)   !LiterfalOrgP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tEXCH_PO4_col(:)       !UPO4(NY,NX)/AREA(3,NU(NY,NX),NY,NX), exchangeable 
  real(r8),pointer   :: histr_1D_SUR_DOP_FLX_col(:)    !HydroDOPFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: histr_1D_SUB_DOP_FLX_col(:)    !HDOPD(NY,NX)/TAREA
  real(r8),pointer   :: histr_1D_SUR_DIP_FLX_col(:)    !HydroDIPFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: histr_1D_SUB_DIP_FLX_col(:)    !HDIPD(NY,NX)/TAREA
  
  real(r8),pointer   :: histr_1D_tPRECIP_P_col(:)       !UPP4(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tMICRO_P_col(:)        !TOPT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_PO4_FIRE_col(:)       !PO4byFire_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_cPO4_LITR_col(:)      !trc_solml_vr(ids_H2PO4,0,NY,NX)/SoilMicPMassLayer(0,NY,NX)
  real(r8),pointer   :: histr_1D_cEXCH_P_LITR_col(:)     !31.0*(trcx_solml(idx_HPO4,0,NY,NX)+trcx_solml(idx_H2PO4,0,NY,NX))/SoilMicPMassLayer(0,NY,NX)
  real(r8),pointer   :: histr_1D_ECO_HVST_P_col(:)     !EcoHavstElmnt_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_NET_P_MIN_col(:)      !-NetPO4Mineralize_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_RADN_col(:)           !TRAD(NY,NX)
  real(r8),pointer   :: histr_1D_tSALT_DISCHG_FLX_col(:)    !HydroIonFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: histr_1D_PSI_SURF_col(:)       !PSISM(0,NY,NX)
  real(r8),pointer   :: histr_1D_SURF_ELEV_col(:)      !-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)+DLYR(3,0,NY,NX)
  real(r8),pointer   :: histr_1D_SURF_tLITR_N_FLX_col(:)      !URSDN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)

  real(r8),pointer   :: histr_1D_AMENDED_N_col(:)       !FertNFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tLITRf_N_FLX_col(:)  !LiterfalOrgN_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tLITRf_C_FLX_col(:)  !LiterfalOrgC_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tNH4X_col(:)           !UNH4(NY,NX)/AREA(3,NU(NY,NX),NY,NX), total NH3+NH4 content
  real(r8),pointer   :: histr_1D_tNO3_col(:)           !UNO3(NY,NX)/AREA(3,NU(NY,NX),NY,NX), total NO3+NO2 content
  real(r8),pointer   :: histr_1D_SUR_DON_FLX_col(:)    !HydroDONFlx_col(NY,NX)/TAREA, daily flux
  real(r8),pointer   :: histr_1D_SUB_DON_FLX_col(:)    !HDOND(NY,NX)/TAREA, daily flux
  real(r8),pointer   :: histr_1D_SUR_DIN_FLX_col(:)    !HydroDINFlx_col(NY,NX)/TAREA
  real(r8),pointer   :: histr_1D_SUB_DIN_FLX_col(:)    !HDIND(NY,NX)/TAREA
  real(r8),pointer   :: histr_1D_tMICRO_N_col(:)        !TONT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_TEMP_LITR_col(:)      !TCS(0,NY,NX)
  real(r8),pointer   :: histr_1D_TEMP_SNOW_col(:)      !TCSnow(1,NY,NX)
  real(r8),pointer   :: histr_1D_SURF_tLITR_C_FLX_col(:)      !URSDC(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_AMENDED_C_col(:)      !AmendCFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_CO2_FLX_col(:)        !UCO2G(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tMICRO_C_col(:)        !TOMT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_OMC_LITR_col(:)       !ORGC(0,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total residual C
  real(r8),pointer   :: histr_1D_SUR_DOC_FLX_col(:)    !HDOCQ(NY,NX)/TAREA
  real(r8),pointer   :: histr_1D_SUB_DOC_FLX_col(:)    !HDOCD(NY,NX)/TAREA
  real(r8),pointer   :: histr_1D_SUR_DIC_FLX_col(:)    !HDICQ(NY,NX)/TAREA
  real(r8),pointer   :: histr_1D_SUB_DIC_FLX_col(:)    !HDICD(NY,NX)/TAREA
  real(r8),pointer   :: histr_1D_ATM_CO2_col(:)        !CO2E(NY,NX)
  real(r8),pointer   :: histr_1D_NBP_col(:)            !Eco_NBP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ECO_HVST_C_col(:)     !EcoHavstElmnt_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ECO_LAI_col(:)        !CanopyLeafArea_grd(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ECO_GPP_col(:)        !Eco_GPP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ECO_RA_col(:)         !Eco_AutoR_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ECO_NPP_col(:)        !Eco_NPP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ECO_RH_col(:)         !Eco_HR_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tDIC_col(:)        !DIC_mass_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX), total soil DIC
  real(r8),pointer   :: histr_1D_tSTANDING_DEAD_C_col(:)       !StandingDeadChemElmnt_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tSTANDING_DEAD_N_col(:)       !StandingDeadChemElmnt_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_tSTANDING_DEAD_P_col(:)       !StandingDeadChemElmnt_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)    
  real(r8),pointer   :: histr_1D_tPRECN_col(:)          !1000.0_r8*URAIN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ET_col(:)             !1000.0_r8*UEVAP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_N2O_LITR_col(:)       !trc_solcl_vr(idg_N2O,0,NY,NX)
  real(r8),pointer   :: histr_1D_NH3_LITR_col(:)       !trc_solcl_vr(idg_NH3,0,NY,NX)
  real(r8),pointer   :: histr_1D_SOL_RADN_col(:)       !RAD(NY,NX)*277.8, W m-2
  real(r8),pointer   :: histr_1D_AIR_TEMP_col(:)       !TCA(NY,NX)
  real(r8),pointer   :: histr_1D_HUM_col(:)            !VPK(NY,NX)
  real(r8),pointer   :: histr_1D_WIND_col(:)           !WindSpeedAtm(NY,NX)/3600.0
  real(r8),pointer   :: histr_1D_PREC_col(:)           !(RainFalPrec(NY,NX)+SnoFalPrec(NY,NX))*1000.0/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_SOIL_RN_col(:)        !HeatByRadiation(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX), 1.e6/3600=277.8
  real(r8),pointer   :: histr_1D_SOIL_LE_col(:)        !HeatEvapAir2Surf(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_SOIL_H_col(:)         !HeatSensAir2Surf(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_SOIL_G_col(:)         !-(HeatNet2Surf(NY,NX)-HeatSensVapAir2Surf(NY,NX))*277.8/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ECO_RN_col(:)         !Eco_NetRad_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ECO_LE_col(:)         !Eco_Heat_Latent_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_Eco_Heat_col(:)          !Eco_Heat_Sens_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ECO_G_col(:)          !Eco_Heat_Grnd_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_O2_LITR_col(:)       !trc_solcl_vr(idg_O2,0,NY,NX)
  real(r8),pointer   :: histr_1D_MIN_LWP_ptc(:)       !PSICanPDailyMin(NZ,NY,NX), minimum daily canopy water potential, [MPa]
  real(r8),pointer   :: histr_1D_SOIL_CO2_FLX_col(:)  !SurfGasFlx(idg_CO2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815, umol m-2 s-1, 1.e6/(12*3600)=23.14815
  real(r8),pointer   :: histr_1D_ECO_CO2_FLX_col(:)   !Eco_NEE_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815
  real(r8),pointer   :: histr_1D_CH4_FLX_col(:)       !SurfGasFlx(idg_CH4,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815
  real(r8),pointer   :: histr_1D_O2_FLX_col(:)        !SurfGasFlx(idg_O2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*8.68056,  umol m-2 s-1, 1.e6/(32*3600)=8.68056
  real(r8),pointer   :: histr_1D_CO2_LITR_col(:)      !trc_solcl_vr(idg_CO2,0,NY,NX)
  real(r8),pointer   :: histr_1D_EVAPN_col(:)          !VapXAir2GSurf(NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_RUNOFF_FLX_col(:)         !-WQRH(NY,NX)*1000.0/TAREA, 
  real(r8),pointer   :: histr_1D_SEDIMENT_FLX_col(:)       !USEDOU(NY,NX)*1000.0/TAREA, soil mass 
  real(r8),pointer   :: histr_1D_tSWC_col(:)        !UVLWatMicP(NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX), volumetric soil water content
  real(r8),pointer   :: histr_1D_DISCHG_FLX_col(:)         !FWatDischarge(NY,NX)*1000.0/TAREA
  real(r8),pointer   :: histr_1D_SNOWPACK_col(:)       !AZMAX1((VOLSS(NY,NX)+VcumIceSnow(NY,NX)*DENSICE+VOLWS(NY,NX))*1000.0/AREA(3,NU(NY,NX),NY,NX))
  real(r8),pointer   :: histr_1D_SURF_WTR_col(:)       !THETWZ(0,NY,NX)
  real(r8),pointer   :: histr_1D_SURF_ICE_col(:)       !THETIZ(0,NY,NX)
  real(r8),pointer   :: histr_1D_ACTV_LYR_col(:)       !-(ActiveLayDepth(NY,NX)-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX))
  real(r8),pointer   :: histr_1D_WTR_TBL_col(:)        !-(DepthInternalWTBL(NY,NX)-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX))
  real(r8),pointer   :: histr_1D_sN2O_FLX_col(:)        !SurfGasFlx(idg_N2O,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_sN2G_FLX_col(:)        !SurfGasFlx(idg_N2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_sNH3_FLX_col(:)        !SurfGasFlx(idg_NH3,NY,NX)/AREA(3,NU(NY,NX),NY,NX)

  real(r8),pointer   :: histr_1D_Plant_C_ptc(:)        !whole plant C  
  real(r8),pointer   :: histr_1D_Plant_N_ptc(:)        !whole plant N  
  real(r8),pointer   :: histr_1D_Plant_P_ptc(:)        !whole plant P  
  real(r8),pointer   :: histr_1D_LEAF_PC_ptc(:)       !(LeafChemElmnts_pft(ielmp,NZ,NY,NX)+CanopyNonstructElements_pft(ielmp,NZ,NY,NX))/(LeafChemElmnts_pft(ielmc,NZ,NY,NX)+CanopyNonstructElements_pft(ielmc,NZ,NY,NX)),mass based CP ratio of leaf
  real(r8),pointer   :: histr_2D_tSOC_vr_col(:,:)        !ORGC(1:JZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total soil C
  real(r8),pointer   :: histr_1D_CAN_RN_ptc(:)        !277.8*RadNet2CanP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), W m-2
  real(r8),pointer   :: histr_1D_CAN_LE_ptc(:)        !277.8*EvapTransHeatP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_CAN_H_ptc(:)         !277.8*HeatXAir2PCan(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_CAN_G_ptc(:)         !277.8*HeatStorCanP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_CAN_TEMP_ptc(:)      !TCelciusCanopy_pft(NZ,NY,NX)
  real(r8),pointer   :: histr_1D_TEMP_FN_ptc(:)       !fTgrowCanP(NZ,NY,NX), canopy temperature growth function/stress
  real(r8),pointer   :: histr_1D_CAN_CO2_FLX_ptc(:)   !CO2NetFix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.148, umol m-2 s-1
  real(r8),pointer   :: histr_1D_CAN_GPP_ptc(:)       !GrossCO2Fix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total gross CO2 fixation, gC m-2
  real(r8),pointer   :: histr_1D_CAN_RA_ptc(:)        !CanopyPlusNoduRespC_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total autotrophic respiration
  real(r8),pointer   :: histr_1D_cTNC_ptc(:)          !CanopyNonstructElementConc_pft(ielmc,NZ,NY,NX), canopy nonstructural C concentration, 
  real(r8),pointer   :: histr_1D_cTNN_ptc(:)          !CanopyNonstructElementConc_pft(ielmn,NZ,NY,NX)
  real(r8),pointer   :: histr_1D_cTNP_ptc(:)          !CanopyNonstructElementConc_pft(ielmp,NZ,NY,NX)
  real(r8),pointer   :: histr_1D_STOML_RSC_CO2_ptc(:) !CanPStomaResistH2O_pft(NZ,NY,NX)*1.56*3600.0_r8, s m-1, for CO2
  real(r8),pointer   :: histr_1D_BLYR_RSC_CO2_ptc(:)  !CanopyBndlResist_pft(NZ,NY,NX)*1.34*3600.0_r8, s m-1, for CO2
  real(r8),pointer   :: histr_1D_CAN_CO2_ptc(:)       !CanopyGasCO2_pft(NZ,NY,NX), umol mol-1
  real(r8),pointer   :: histr_1D_LAI_ptc(:)           !CanopyArea_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant leaf area, include stalk
  real(r8),pointer   :: histr_1D_PSI_CAN_ptc(:)       !PSICanopy_pft(NZ,NY,NX), canopy total water potential , MPa
  real(r8),pointer   :: histr_1D_TURG_CAN_ptc(:)      !PSICanopyTurg_pft(NZ,NY,NX), canopy turgor water potential, MPa
  real(r8),pointer   :: histr_1D_STOM_RSC_H2O_ptc(:)  !CanPStomaResistH2O_pft(NZ,NY,NX)*3600.0_r8, s m-1, for H2O
  real(r8),pointer   :: histr_1D_BLYR_RSC_H2O_ptc(:)  !CanopyBndlResist_pft(NZ,NY,NX)*3600.0_r8, s m-1, for H2O
  real(r8),pointer   :: histr_1D_TRANSPN_ptc(:)       !Transpiration_pft(NZ,NY,NX)*1000.0_r8/AREA(3,NU(NY,NX),NY,NX), canopy transpiration mm H2O/m2/h
  real(r8),pointer   :: histr_1D_O2_STRESS_ptc(:)     !OSTR(NZ,NY,NX), plant O2 stress indicator
  real(r8),pointer   :: histr_1D_NH4_UPTK_FLX_ptc(:)      !RootNH4Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_NO3_UPTK_FLX_ptc(:)      !RootNO3Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_N2_FIXN_FLX_ptc(:)       !RootN2Fix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_cNH3_FLX_ptc(:)       !RNH3C(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_PO4_UPTK_FLX_ptc(:)      !RootH2PO4Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_SHOOT_C_ptc(:)       !ShootChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_LEAF_C_ptc(:)        !LeafChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_Petiole_C_ptc(:)        !PetioleChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), canopy sheath element
  real(r8),pointer   :: histr_1D_STALK_C_ptc(:)       !StalkChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_RESERVE_C_ptc(:)     !ReserveChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_HUSK_C_ptc(:)        !(HuskChemElmnts_pft(ielmc,NZ,NY,NX)+EarChemElmnts_pft(ielmc,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_GRAIN_C_ptc(:)       !GrainChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ROOT_C_ptc(:)        !RootChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_NODULE_C_ptc(:)      !NoduleChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), nodule
  real(r8),pointer   :: histr_1D_STORED_C_ptc(:)      !NonstructalChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_GRAIN_NO_ptc(:)      !CanopySeedNumber_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_LAIb_ptc(:)          !CanopyLeafArea_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total branch leaf area
  real(r8),pointer   :: histr_1D_EXUD_C_FLX_ptc(:)        !PlantExudChemElmntCum_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_LITRf_C_FLX_ptc(:)       !LitrfallChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_LITRf_P_FLX_ptc(:)       !LitrfallChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_SURF_LITRf_C_FLX_ptc(:)  !SurfLitrfallChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_AUTO_RESP_FLX_ptc(:)     !GrossResp_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ABV_GRD_RESP_FLX_ptc(:)  !CanopyPlusNoduRespC_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_HVST_C_FLX_ptc(:)        !EcoHavstElmnt_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_PLANT_BALANCE_C_ptc(:)     !ElmntBalanceCum_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_STANDING_DEAD_C_ptc(:)    !StandingDeadChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_FIREp_CO2_FLX_ptc(:)     !CO2ByFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant CO2 from fire
  real(r8),pointer   :: histr_1D_FIREp_CH4_FLX_ptc(:)     !CH4ByFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_NPP_ptc(:)           !NetPrimaryProductvity_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_CAN_HT_ptc(:)        !CanopyHeight_pft(NZ,NY,NX), canopy height, m
  real(r8),pointer   :: histr_1D_POPN_ptc(:)          !PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant population
  real(r8),pointer   :: histr_1D_tTRANSPN_ptc(:)      !-ETCanopy_pft(NZ,NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX), total transpiration
  real(r8),pointer   :: histr_1D_WTR_STRESS_ptc(:)    !HoursCanopyPSITooLow(NZ,NY,NX)
  real(r8),pointer   :: histr_1D_OXY_STRESS_ptc(:)    !OSTR(NZ,NY,NX)
  real(r8),pointer   :: histr_1D_SHOOT_N_ptc(:)       !ShootChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_LEAF_N_ptc(:)        !LeafChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_Petiole_N_ptc(:)        !PetioleChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_STALK_N_ptc(:)       !StalkChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_RESERVE_N_ptc(:)     !ReserveChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_HUSK_N_ptc(:)        !(HuskChemElmnts_pft(ielmn,NZ,NY,NX)+EarChemElmnts_pft(ielmn,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_GRAIN_N_ptc(:)       !GrainChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ROOT_N_ptc(:)        !RootChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_NODULE_N_ptc(:)         !NoduleChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_STORED_N_ptc(:)      !NonstructalChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_EXUD_N_FLX_ptc(:)        !PlantExudChemElmntCum_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_LITRf_N_FLX_ptc(:)       !LitrfallChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total plant litterfall N
  real(r8),pointer   :: histr_1D_TL_N_FIXED_FLX_ptc(:)    !PlantN2FixCum_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total plant N2 fixation
  real(r8),pointer   :: histr_1D_HVST_N_FLX_ptc(:)        !EcoHavstElmnt_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_NH3can_FLX_ptc(:)    !NH3EmiCum_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_PLANT_BALANCE_N_ptc(:)     !ElmntBalanceCum_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_STANDING_DEAD_N_ptc(:)    !StandingDeadChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_FIREp_N_FLX_ptc(:)        !NH3byFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant N emission from fire
  real(r8),pointer   :: histr_1D_SURF_LITRf_N_FLX_ptc(:)   !SurfLitrfallChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), surface litter fall
  real(r8),pointer   :: histr_1D_SHOOT_P_ptc(:)       !ShootChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_LEAF_P_ptc(:)        !LeafChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_Petiole_P_ptc(:)        !PetioleChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_STALK_P_ptc(:)       !StalkChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_RESERVE_P_ptc(:)     !ReserveChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_HUSK_P_ptc(:)        !(HuskChemElmnts_pft(ielmp,NZ,NY,NX)+EarChemElmnts_pft(ielmp,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_GRAIN_P_ptc(:)       !GrainChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_ROOT_P_ptc(:)        !RootChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_NODULE_P_ptc(:)         !NoduleChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_STORED_P_ptc(:)      !NonstructalChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_EXUD_P_FLX_ptc(:)        !PlantExudChemElmntCum_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_LITTERf_P_ptc(:)     !LitrfallChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_HVST_P_FLX_ptc(:)        !EcoHavstElmnt_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_PLANT_BALANCE_P_ptc(:)     !ElmntBalanceCum_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_STANDING_DEAD_P_ptc(:)    !StandingDeadChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_FIREp_P_FLX_ptc(:)        !PO4byFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_SURF_LITRf_P_FLX_ptc(:)  !SurfLitrfallChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  real(r8),pointer   :: histr_1D_BRANCH_NO_ptc(:)     !NumOfBranches_pft(NZ,NY,NX)
  real(r8),pointer   :: histr_1D_LEAF_NC_ptc(:)       !(LeafChemElmnts_pft(ielmn,NZ,NY,NX)+CanopyNonstructElements_pft(ielmn,NZ,NY,NX))/(LeafChemElmnts_pft(ielmc,NZ,NY,NX)+CanopyNonstructElements_pft(ielmc,NZ,NY,NX)),mass based CN ratio of leaf  
  real(r8),pointer   :: histr_1D_Growth_Stage_ptc(:)    !plant development stage, integer, 0-10, planting, emergence, floral_init, jointing, 
                                      !elongation, heading, anthesis, seed_fill, see_no_set, seed_mass_set, end_seed_fill
  real(r8),pointer   :: histr_2D_LEAF_NODE_NO_ptc(:,:)       !NumOfLeaves_brch(NumOfMainBranch_pft(NZ,NY,NX),NZ,NY,NX), leaf NO
  real(r8),pointer   :: histr_2D_RUB_ACTVN_ptc(:,:)     !RubiscoActivity_brpft(NumOfMainBranch_pft(NZ,NY,NX),NZ,NY,NX), branch down-regulation of CO2 fixation
  real(r8),pointer   :: histr_2D_CO2_vr_col(:,:)        !trc_solcl_vr(idg_CO2,1:JZ,NY,NX)
  real(r8),pointer   :: histr_2D_CH4_vr_col(:,:)        !trc_solcl_vr(idg_CH4,1:JZ,NY,NX)
  real(r8),pointer   :: histr_2D_O2_vr_col(:,:)         !trc_solcl_vr(idg_O2,1:JZ,NY,NX)
  real(r8),pointer   :: histr_2D_N2O_vr_col(:,:)         !trc_solcl_vr(idg_N2O,1:JZ,NY,NX)
  real(r8),pointer   :: histr_2D_NH3_vr_col(:,:)         !trc_solcl_vr(idg_NH3,1:JZ,NY,NX)
  real(r8),pointer   :: histr_2D_TEMP_vr_col(:,:)        !TCS(1:JZ,NY,NX)
  real(r8),pointer   :: histr_2D_vWATER_vr_col(:,:)       !THETWZ(1:JZ,NY,NX)
  real(r8),pointer   :: histr_2D_vICE_vr_col(:,:)         !THETIZ(1:JZ,NY,NX)
  real(r8),pointer   :: histr_2D_PSI_vr_col(:,:)         !PSISM(1:JZ,NY,NX)+PSISO(1:JZ,NY,NX)
  real(r8),pointer   :: histr_2D_cNH4t_vr_col(:,:)       !(trc_solml_vr(ids_NH4,1:JZ,NY,NX)+trc_solml_vr(ids_NH4B,1:JZ,NY,NX) &
                                                                  !+14.0*(trcx_solml(idx_NH4,1:JZ,NY,NX)+trcx_solml(idx_NH4B,1:JZ,NY,NX)))/SoilMicPMassLayer(1:JZ,NY,NX)
  real(r8),pointer   :: histr_2D_cNO3t_vr_col(:,:)       !(trc_solml_vr(ids_NO3,1:JZ,NY,NX)+trc_solml_vr(ids_NO3B,1:JZ,NY,NX) &
                                                                  !+trc_solml_vr(ids_NO2,1,NY,NX)+trc_solml_vr(ids_NO2B,1,NY,NX))/SoilMicPMassLayer(1,NY,NX)
  real(r8),pointer   :: histr_2D_cPO4_vr_col(:,:)        !(trc_solml_vr(ids_H1PO4,1:JZ,NY,NX)+trc_solml_vr(ids_H1PO4B,1,NY,NX)+trc_solml_vr(ids_H2PO4,1,NY,NX)+trc_solml_vr(ids_H2PO4B,1,NY,NX))/VLWatMicP(1,NY,NX)
  real(r8),pointer   :: histr_2D_cEXCH_P_vr_col(:,:)     !31.0*(trcx_solml(idx_HPO4,1,NY,NX)+trcx_solml(idx_H2PO4,1,NY,NX)+trcx_solml(idx_HPO4B,1,NY,NX)+trcx_solml(idx_H2PO4B,1,NY,NX))/SoilMicPMassLayer(1,NY,NX)
  real(r8),pointer   :: histr_2D_ECND_vr_col(:,:)         !ECND(1:JZ,NY,NX)

  real(r8),pointer   :: histr_2D_PSI_RT_vr_ptc(:,:)     !PSIRoot_vr(1,1:JZ,NZ,NY,NX), root total water potential , MPa
  real(r8),pointer   :: histr_2D_prtUP_NH4_vr_ptc(:,:)     !(RootNutUptake_pvr(ids_NH4,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NH4,2,1:JZ,NZ,NY,NX) &
                                                                   !+RootNutUptake_pvr(ids_NH4B,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NH4B,2,1:JZ,NZ,NY,NX))/AREA(3,1,NY,NX)
  real(r8),pointer   :: histr_2D_prtUP_NO3_vr_ptc(:,:)     !(RootNutUptake_pvr(ids_NO3,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NO3,2,1:JZ,NZ,NY,NX) &
                                                                   !+RootNutUptake_pvr(ids_NO3B,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NO3B,2,1:JZ,NZ,NY,NX))/AREA(3,1,NY,NX)
  real(r8),pointer   :: histr_2D_prtUP_PO4_vr_ptc(:,:)     !(RootNutUptake_pvr(ids_H2PO4,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_H2PO4,2,1:JZ,NZ,NY,NX) &
                                                                   !+RootNutUptake_pvr(ids_H2PO4B,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_H2PO4B,2,1:JZ,NZ,NY,NX))/AREA(3,1,NY,NX)
  real(r8),pointer   :: histr_2D_DNS_RT_vr_ptc(:,:)     !RootLenthDensPerPopu_pvr(1,1:JZ,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
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
  real(r8), pointer :: data2d_ptr(:,:)

  beg_col=1;end_col=bounds%ncols
  beg_ptc=1;end_ptc=bounds%npfts
  allocate(this%histr_1D_tFIRE_CO2_col(beg_col:end_col))        !CO2byFire_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_tFIRE_CH4_col(beg_col:end_col))        !CH4byFire_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_cNH4_LITR_col(beg_col:end_col))        !(trc_solml_vr(ids_NH4,0,NY,NX)+14.0*trcx_solml(idx_NH4,0,NY,NX))/SoilMicPMassLayer(0,NY,NX)
  allocate(this%histr_1D_cNO3_LITR_col(beg_col:end_col))        !(trc_solml_vr(ids_NO3,0,NY,NX)+trc_solml_vr(ids_NO2,0,NY,NX))/SoilMicPMassLayer(0,NY,NX)                            
  allocate(this%histr_1D_ECO_HVST_C_col(beg_col:end_col))       !EcoHavstElmnt_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)  
  allocate(this%histr_1D_ECO_HVST_N_col(beg_col:end_col))      !EcoHavstElmnt_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ECO_HVST_P_col(beg_col:end_col))      !EcoHavstElmnt_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_NET_N_MIN_col(beg_col:end_col))       !-NetNH4Mineralize_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_SURF_tLITR_P_FLX_col(beg_col:end_col))       !URSDP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_HUMUS_C_col(beg_col:end_col))       !UORGC(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_HUMUS_N_col(beg_col:end_col))         !UORGN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_HUMUS_P_col(beg_col:end_col))         !UORGP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_AMENDED_P_col(beg_col:end_col))        !FerPFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_tLITRf_C_FLX_col(beg_col:end_col))   !LiterfalOrgC_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_tLITRf_N_FLX_col(beg_col:end_col))   !LiterfalOrgN_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_tLITRf_P_FLX_col(beg_col:end_col))   !LiterfalOrgP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_tEXCH_PO4_col(beg_col:end_col))       !UPO4(NY,NX)/AREA(3,NU(NY,NX),NY,NX), exchangeable 
  allocate(this%histr_1D_SUR_DOP_FLX_col(beg_col:end_col))     !HydroDOPFlx_col(NY,NX)/TAREA
  allocate(this%histr_1D_SUB_DOP_FLX_col(beg_col:end_col))     !HDOPD(NY,NX)/TAREA
  allocate(this%histr_1D_SUR_DIP_FLX_col(beg_col:end_col))     !HydroDIPFlx_col(NY,NX)/TAREA
  allocate(this%histr_1D_SUB_DIP_FLX_col(beg_col:end_col))     !HDIPD(NY,NX)/TAREA
  allocate(this%histr_1D_tSALT_DISCHG_FLX_col(beg_col:end_col)) !HydroIonFlx_col(NY,NX)/TAREA
  allocate(this%histr_1D_SUR_DON_FLX_col(beg_col:end_col))     !HydroDONFlx_col(NY,NX)/TAREA, daily flux
  allocate(this%histr_1D_SUB_DON_FLX_col(beg_col:end_col))     !HDOND(NY,NX)/TAREA, daily flux
  allocate(this%histr_1D_SUR_DIN_FLX_col(beg_col:end_col))     !HydroDINFlx_col(NY,NX)/TAREA
  allocate(this%histr_1D_SUB_DIN_FLX_col(beg_col:end_col))     !HDIND(NY,NX)/TAREA
  allocate(this%histr_1D_SUR_DOC_FLX_col(beg_col:end_col))     !HDOCQ(NY,NX)/TAREA
  allocate(this%histr_1D_SUB_DOC_FLX_col(beg_col:end_col))     !HDOCD(NY,NX)/TAREA
  allocate(this%histr_1D_SUR_DIC_FLX_col(beg_col:end_col))     !HDICQ(NY,NX)/TAREA
  allocate(this%histr_1D_SUB_DIC_FLX_col(beg_col:end_col))     !HDICD(NY,NX)/TAREA

  allocate(this%histr_1D_tPRECIP_P_col(beg_col:end_col))        !UPP4(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_tMICRO_P_col(beg_col:end_col))         !TOPT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_PO4_FIRE_col(beg_col:end_col))        !PO4byFire_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_cPO4_LITR_col(beg_col:end_col))       !trc_solml_vr(ids_H2PO4,0,NY,NX)/SoilMicPMassLayer(0,NY,NX)
  allocate(this%histr_1D_cEXCH_P_LITR_col(beg_col:end_col))      !31.0*(trcx_solml(idx_HPO4,0,NY,NX)+trcx_solml(idx_H2PO4,0,NY,NX))/SoilMicPMassLayer(0,NY,NX)
  allocate(this%histr_1D_NET_P_MIN_col(beg_col:end_col))       !-NetPO4Mineralize_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_RADN_col(beg_col:end_col))            !TRAD(NY,NX)

  allocate(this%histr_1D_PSI_SURF_col(beg_col:end_col))        !PSISM(0,NY,NX)
  allocate(this%histr_1D_SURF_ELEV_col(beg_col:end_col))       !-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)+DLYR(3,0,NY,NX)
  allocate(this%histr_1D_SURF_tLITR_C_FLX_col(beg_col:end_col))       !URSDC(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_SURF_tLITR_N_FLX_col(beg_col:end_col))       !URSDN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)

  allocate(this%histr_1D_AMENDED_N_col(beg_col:end_col))        !FertNFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)

  allocate(this%histr_1D_tNH4X_col(beg_col:end_col))          !UNH4(NY,NX)/AREA(3,NU(NY,NX),NY,NX), total NH3+NH4 content
  allocate(this%histr_1D_tNO3_col(beg_col:end_col))          !UNO3(NY,NX)/AREA(3,NU(NY,NX),NY,NX), total NO3+NO2 content
  allocate(this%histr_1D_tMICRO_N_col(beg_col:end_col))       !TONT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_TEMP_LITR_col(beg_col:end_col))     !TCS(0,NY,NX), oC
  allocate(this%histr_1D_TEMP_SNOW_col(beg_col:end_col))     !TCSnow(1,NY,NX), oC

  allocate(this%histr_1D_AMENDED_C_col(beg_col:end_col))     !AmendCFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_CO2_FLX_col(beg_col:end_col))       !UCO2G(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_tMICRO_C_col(beg_col:end_col))       !TOMT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_OMC_LITR_col(beg_col:end_col))      !ORGC(0,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total residual C
  allocate(this%histr_1D_ATM_CO2_col(beg_col:end_col))       !CO2E(NY,NX)
  allocate(this%histr_1D_NBP_col(beg_col:end_col))           !Eco_NBP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)

  allocate(this%histr_1D_ECO_LAI_col(beg_col:end_col))       !CanopyLeafArea_grd(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ECO_GPP_col(beg_col:end_col))       !Eco_GPP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ECO_RA_col(beg_col:end_col))        !Eco_AutoR_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ECO_NPP_col(beg_col:end_col))       !Eco_NPP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ECO_RH_col(beg_col:end_col))        !Eco_HR_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_tDIC_col(beg_col:end_col))       ;  this%histr_1D_tDIC_col=spval
  allocate(this%histr_1D_tSTANDING_DEAD_C_col(beg_col:end_col));  this%histr_1D_tSTANDING_DEAD_C_col=spval 
  allocate(this%histr_1D_tSTANDING_DEAD_N_col(beg_col:end_col));  this%histr_1D_tSTANDING_DEAD_N_col=spval  
  allocate(this%histr_1D_tSTANDING_DEAD_P_col(beg_col:end_col));  this%histr_1D_tSTANDING_DEAD_P_col=spval  
  allocate(this%histr_1D_tPRECN_col(beg_col:end_col))        !1000.0_r8*URAIN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ET_col(beg_col:end_col))            !1000.0_r8*UEVAP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_N2O_LITR_col(beg_col:end_col))      !trc_solcl_vr(idg_N2O,0,NY,NX)
  allocate(this%histr_1D_NH3_LITR_col(beg_col:end_col))      !trc_solcl_vr(idg_NH3,0,NY,NX)
  allocate(this%histr_1D_SOL_RADN_col(beg_col:end_col))      !RAD(NY,NX)*277.8, W m-2
  allocate(this%histr_1D_AIR_TEMP_col(beg_col:end_col))      !TCA(NY,NX), oC
  allocate(this%histr_1D_HUM_col(beg_col:end_col))           !VPK(NY,NX), kPa
  allocate(this%histr_1D_WIND_col(beg_col:end_col))          !WindSpeedAtm(NY,NX)/3600.0, m/s
  allocate(this%histr_1D_PREC_col(beg_col:end_col))          !(RainFalPrec(NY,NX)+SnoFalPrec(NY,NX))*1000.0/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_SOIL_RN_col(beg_col:end_col))       !HeatByRadiation(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX), W m-2
  allocate(this%histr_1D_SOIL_LE_col(beg_col:end_col))       !HeatEvapAir2Surf(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_SOIL_H_col(beg_col:end_col))        !HeatSensAir2Surf(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_SOIL_G_col(beg_col:end_col))        !-(HeatNet2Surf(NY,NX)-HeatSensVapAir2Surf(NY,NX))*277.8/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ECO_RN_col(beg_col:end_col))        !Eco_NetRad_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX), W m-2
  allocate(this%histr_1D_ECO_LE_col(beg_col:end_col))        !Eco_Heat_Latent_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_Eco_Heat_col(beg_col:end_col))         !Eco_Heat_Sens_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ECO_G_col(beg_col:end_col))         !Eco_Heat_Grnd_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_O2_LITR_col(beg_col:end_col))       !trc_solcl_vr(idg_O2,0,NY,NX)
  allocate(this%histr_1D_MIN_LWP_ptc(beg_ptc:end_ptc))       !PSICanPDailyMin(NZ,NY,NX), minimum daily canopy water potential, [MPa]
  allocate(this%histr_1D_SOIL_CO2_FLX_col(beg_col:end_col))  !SurfGasFlx(idg_CO2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815, umol m-2 s-1, 1.e6/(12*3600)=23.14815
  allocate(this%histr_1D_ECO_CO2_FLX_col(beg_col:end_col))   !TCO2NetFix_pft(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815
  allocate(this%histr_1D_CH4_FLX_col(beg_col:end_col))       !SurfGasFlx(idg_CH4,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815
  allocate(this%histr_1D_O2_FLX_col(beg_col:end_col))        !SurfGasFlx(idg_O2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*8.68056,  umol m-2 s-1, 1.e6/(32*3600)=8.68056
  allocate(this%histr_1D_CO2_LITR_col(beg_col:end_col))      !trc_solcl_vr(idg_CO2,0,NY,NX)
  allocate(this%histr_1D_EVAPN_col(beg_col:end_col))         !VapXAir2GSurf(NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX), mm H2O/h/m2

  allocate(this%histr_1D_tSWC_col(beg_col:end_col))       !UVLWatMicP(NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX), volumetric soil water content
  allocate(this%histr_1D_SNOWPACK_col(beg_col:end_col))      !AZMAX1((VcumDrySnoWE(NY,NX)+VcumIceSnow(NY,NX)*DENSICE+VOLWS(NY,NX))*1000.0/AREA(3,NU(NY,NX),NY,NX))
  allocate(this%histr_1D_SURF_WTR_col(beg_col:end_col))      !THETWZ(0,NY,NX)
  allocate(this%histr_1D_SURF_ICE_col(beg_col:end_col))      !THETIZ(0,NY,NX)
  allocate(this%histr_1D_ACTV_LYR_col(beg_col:end_col))      !-(ActiveLayDepth(NY,NX)-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX))
  allocate(this%histr_1D_WTR_TBL_col(beg_col:end_col))       !-(DepthInternalWTBL(NY,NX)-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX))
  allocate(this%histr_1D_sN2O_FLX_col(beg_col:end_col))       !SurfGasFlx(idg_N2O,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_sN2G_FLX_col(beg_col:end_col))       !SurfGasFlx(idg_N2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_sNH3_FLX_col(beg_col:end_col))       !SurfGasFlx(idg_NH3,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_RUNOFF_FLX_col(beg_col:end_col))        !-WQRH(NY,NX)*1000.0/TAREA, 
  allocate(this%histr_1D_SEDIMENT_FLX_col(beg_col:end_col))      !USEDOU(NY,NX)*1000.0/TAREA, soil mass 
  allocate(this%histr_1D_DISCHG_FLX_col(beg_col:end_col))        !FWatDischarge(NY,NX)*1000.0/TAREA
  allocate(this%histr_1D_LEAF_PC_ptc(beg_ptc:end_ptc))       !(LeafChemElmnts_pft(ielmp,NZ,NY,NX)+CanopyNonstructElements_pft(ielmp,NZ,NY,NX))/(LeafChemElmnts_pft(ielmc,NZ,NY,NX)+CanopyNonstructElements_pft(ielmc,NZ,NY,NX)),mass based CP ratio of leaf
  allocate(this%histr_1D_CAN_RN_ptc(beg_ptc:end_ptc))        !277.8*RadNet2CanP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), W m-2
  allocate(this%histr_1D_CAN_LE_ptc(beg_ptc:end_ptc))        !277.8*EvapTransHeatP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_CAN_H_ptc(beg_ptc:end_ptc))         !277.8*HeatXAir2PCan(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_CAN_G_ptc(beg_ptc:end_ptc))         !277.8*HeatStorCanP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_CAN_TEMP_ptc(beg_ptc:end_ptc))     !TCelciusCanopy_pft(NZ,NY,NX)
  allocate(this%histr_1D_TEMP_FN_ptc(beg_ptc:end_ptc))      !fTgrowCanP(NZ,NY,NX), canopy temperature growth function/stress
  allocate(this%histr_1D_CAN_CO2_FLX_ptc(beg_ptc:end_ptc))  !CO2NetFix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.148, umol m-2 s-1
  allocate(this%histr_1D_CAN_GPP_ptc(beg_ptc:end_ptc))      !GrossCO2Fix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total gross CO2 fixation, gC m-2
  allocate(this%histr_1D_CAN_RA_ptc(beg_ptc:end_ptc))       !CanopyPlusNoduRespC_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total autotrophic respiration
  allocate(this%histr_1D_cTNC_ptc(beg_ptc:end_ptc))         !CanopyNonstructElementConc_pft(ielmc,NZ,NY,NX), canopy nonstructural C concentration, 
  allocate(this%histr_1D_cTNN_ptc(beg_ptc:end_ptc))         !CanopyNonstructElementConc_pft(ielmn,NZ,NY,NX)
  allocate(this%histr_1D_cTNP_ptc(beg_ptc:end_ptc))         !CanopyNonstructElementConc_pft(ielmp,NZ,NY,NX)
  allocate(this%histr_1D_STOML_RSC_CO2_ptc(beg_ptc:end_ptc))!CanPStomaResistH2O_pft(NZ,NY,NX)*1.56*3600.0_r8, s m-1, for CO2
  allocate(this%histr_1D_BLYR_RSC_CO2_ptc(beg_ptc:end_ptc)) !CanopyBndlResist_pft(NZ,NY,NX)*1.34*3600.0_r8, s m-1, for CO2
  allocate(this%histr_1D_CAN_CO2_ptc(beg_ptc:end_ptc))      !CanopyGasCO2_pft(NZ,NY,NX), umol mol-1
  allocate(this%histr_1D_LAI_ptc(beg_ptc:end_ptc))          !CanopyArea_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant leaf area, include stalk
  allocate(this%histr_1D_PSI_CAN_ptc(beg_ptc:end_ptc))      !PSICanopy_pft(NZ,NY,NX), canopy total water potential , MPa
  allocate(this%histr_1D_TURG_CAN_ptc(beg_ptc:end_ptc))     !PSICanopyTurg_pft(NZ,NY,NX), canopy turgor water potential, MPa
  allocate(this%histr_1D_STOM_RSC_H2O_ptc(beg_ptc:end_ptc)) !CanPStomaResistH2O_pft(NZ,NY,NX)*3600.0_r8, s m-1, for H2O
  allocate(this%histr_1D_BLYR_RSC_H2O_ptc(beg_ptc:end_ptc)) !CanopyBndlResist_pft(NZ,NY,NX)*3600.0_r8, s m-1, for H2O
  allocate(this%histr_1D_TRANSPN_ptc(beg_ptc:end_ptc))      !Transpiration_pft(NZ,NY,NX)*1000.0_r8/AREA(3,NU(NY,NX),NY,NX), canopy transpiration mm H2O/m2/h
  allocate(this%histr_1D_O2_STRESS_ptc(beg_ptc:end_ptc))    !OSTR(NZ,NY,NX), plant O2 stress indicator
  allocate(this%histr_1D_NH4_UPTK_FLX_ptc(beg_ptc:end_ptc))     !RootNH4Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_NO3_UPTK_FLX_ptc(beg_ptc:end_ptc))     !RootNO3Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_N2_FIXN_FLX_ptc(beg_ptc:end_ptc))      !RootN2Fix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_cNH3_FLX_ptc(beg_ptc:end_ptc))      !RNH3C(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_PO4_UPTK_FLX_ptc(beg_ptc:end_ptc))     !RootH2PO4Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_SHOOT_C_ptc(beg_ptc:end_ptc))      !ShootChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_Plant_C_ptc(beg_ptc:end_ptc))
  allocate(this%histr_1D_LEAF_C_ptc(beg_ptc:end_ptc))       !LeafChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_Petiole_C_ptc(beg_ptc:end_ptc))       !PetioleChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), canopy sheath element
  allocate(this%histr_1D_STALK_C_ptc(beg_ptc:end_ptc))      !StalkChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_RESERVE_C_ptc(beg_ptc:end_ptc))    !ReserveChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_HUSK_C_ptc(beg_ptc:end_ptc))       !(HuskChemElmnts_pft(ielmc,NZ,NY,NX)+EarChemElmnts_pft(ielmc,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_GRAIN_C_ptc(beg_ptc:end_ptc))      !GrainChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ROOT_C_ptc(beg_ptc:end_ptc))       !RootChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_NODULE_C_ptc(beg_ptc:end_ptc))        !NoduleChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), nodule
  allocate(this%histr_1D_STORED_C_ptc(beg_ptc:end_ptc))     !NonstructalChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_GRAIN_NO_ptc(beg_ptc:end_ptc))     !CanopySeedNumber_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_LAIb_ptc(beg_ptc:end_ptc))         !CanopyLeafArea_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total branch leaf area
  allocate(this%histr_1D_EXUD_C_FLX_ptc(beg_ptc:end_ptc))       !PlantExudChemElmntCum_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_LITRf_C_FLX_ptc(beg_ptc:end_ptc));  this%histr_1D_LITRf_C_FLX_ptc=spval
  allocate(this%histr_1D_LITRf_P_FLX_ptc(beg_ptc:end_ptc));  this%histr_1D_LITRf_P_FLX_ptc=spval
  allocate(this%histr_1D_SURF_LITRf_C_FLX_ptc(beg_ptc:end_ptc)) !SurfLitrfallChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_AUTO_RESP_FLX_ptc(beg_ptc:end_ptc))    !GrossResp_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ABV_GRD_RESP_FLX_ptc(beg_ptc:end_ptc))  !CanopyPlusNoduRespC_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_HVST_C_FLX_ptc(beg_ptc:end_ptc))        !EcoHavstElmnt_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_STANDING_DEAD_C_ptc(beg_ptc:end_ptc))    !StandingDeadChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_FIREp_CO2_FLX_ptc(beg_ptc:end_ptc))     !CO2ByFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant CO2 from fire
  allocate(this%histr_1D_FIREp_CH4_FLX_ptc(beg_ptc:end_ptc))     !CH4ByFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_NPP_ptc(beg_ptc:end_ptc))           !NetPrimaryProductvity_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_CAN_HT_ptc(beg_ptc:end_ptc))        !CanopyHeight_pft(NZ,NY,NX), canopy height, m
  allocate(this%histr_1D_POPN_ptc(beg_ptc:end_ptc))          !PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant population
  allocate(this%histr_1D_tTRANSPN_ptc(beg_ptc:end_ptc))      !-ETCanopy_pft(NZ,NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX), total transpiration
  allocate(this%histr_1D_WTR_STRESS_ptc(beg_ptc:end_ptc))    !HoursCanopyPSITooLow(NZ,NY,NX)
  allocate(this%histr_1D_OXY_STRESS_ptc(beg_ptc:end_ptc))    !OSTR(NZ,NY,NX)
  allocate(this%histr_1D_SHOOT_N_ptc(beg_ptc:end_ptc))       !ShootChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_Plant_N_ptc(beg_ptc:end_ptc))  
  allocate(this%histr_1D_LEAF_N_ptc(beg_ptc:end_ptc))        !LeafChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_Petiole_N_ptc(beg_ptc:end_ptc))        !PetioleChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_STALK_N_ptc(beg_ptc:end_ptc))       !StalkChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_RESERVE_N_ptc(beg_ptc:end_ptc))     !ReserveChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_HUSK_N_ptc(beg_ptc:end_ptc))        !(HuskChemElmnts_pft(ielmn,NZ,NY,NX)+EarChemElmnts_pft(ielmn,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_GRAIN_N_ptc(beg_ptc:end_ptc))       !GrainChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ROOT_N_ptc(beg_ptc:end_ptc))        !RootChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_NODULE_N_ptc(beg_ptc:end_ptc))         !NoduleChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_STORED_N_ptc(beg_ptc:end_ptc))      !NonstructalChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_EXUD_N_FLX_ptc(beg_ptc:end_ptc))        !PlantExudChemElmntCum_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_LITRf_N_FLX_ptc(beg_ptc:end_ptc))       !LitrfallChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total plant litterfall N
  allocate(this%histr_1D_TL_N_FIXED_FLX_ptc(beg_ptc:end_ptc))    !PlantN2FixCum_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total plant N2 fixation
  allocate(this%histr_1D_HVST_N_FLX_ptc(beg_ptc:end_ptc))        !EcoHavstElmnt_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_NH3can_FLX_ptc(beg_ptc:end_ptc))    !NH3EmiCum_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_PLANT_BALANCE_C_ptc(beg_ptc:end_ptc)); this%histr_1D_PLANT_BALANCE_C_ptc=spval
  allocate(this%histr_1D_PLANT_BALANCE_N_ptc(beg_ptc:end_ptc)); this%histr_1D_PLANT_BALANCE_N_ptc=spval
  allocate(this%histr_1D_PLANT_BALANCE_P_ptc(beg_ptc:end_ptc)); this%histr_1D_PLANT_BALANCE_P_ptc=spval
  allocate(this%histr_1D_STANDING_DEAD_N_ptc(beg_ptc:end_ptc))    !StandingDeadChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_FIREp_N_FLX_ptc(beg_ptc:end_ptc))        !NH3byFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), plant N emission from fire
  allocate(this%histr_1D_SURF_LITRf_N_FLX_ptc(beg_ptc:end_ptc))   !SurfLitrfallChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), surface litter fall
  allocate(this%histr_1D_SHOOT_P_ptc(beg_ptc:end_ptc))       !ShootChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_Plant_P_ptc(beg_ptc:end_ptc))
  allocate(this%histr_1D_LEAF_P_ptc(beg_ptc:end_ptc))        !LeafChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_Petiole_P_ptc(beg_ptc:end_ptc))        !PetioleChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_STALK_P_ptc(beg_ptc:end_ptc))       !StalkChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_RESERVE_P_ptc(beg_ptc:end_ptc))     !ReserveChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_HUSK_P_ptc(beg_ptc:end_ptc))        !(HuskChemElmnts_pft(ielmp,NZ,NY,NX)+EarChemElmnts_pft(ielmp,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_GRAIN_P_ptc(beg_ptc:end_ptc))       !GrainChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_ROOT_P_ptc(beg_ptc:end_ptc))        !RootChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_NODULE_P_ptc(beg_ptc:end_ptc))         !NoduleChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_STORED_P_ptc(beg_ptc:end_ptc))      !NonstructalChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_EXUD_P_FLX_ptc(beg_ptc:end_ptc))        !PlantExudChemElmntCum_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_LITTERf_P_ptc(beg_ptc:end_ptc))     !LitrfallChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_HVST_P_FLX_ptc(beg_ptc:end_ptc))        !EcoHavstElmnt_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_STANDING_DEAD_P_ptc(beg_ptc:end_ptc))    !StandingDeadChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_FIREp_P_FLX_ptc(beg_ptc:end_ptc))               !PO4byFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_SURF_LITRf_P_FLX_ptc(beg_ptc:end_ptc))         !SurfLitrfallChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  allocate(this%histr_1D_BRANCH_NO_ptc(beg_ptc:end_ptc))            !NumOfBranches_pft(NZ,NY,NX)
  allocate(this%histr_1D_Growth_Stage_ptc(beg_ptc:end_ptc));      this%histr_1D_Growth_Stage_ptc=spval
  allocate(this%histr_1D_LEAF_NC_ptc(beg_ptc:end_ptc))              !(LeafChemElmnts_pft(ielmn,NZ,NY,NX)+CanopyNonstructElements_pft(ielmn,NZ,NY,NX))/(LeafChemElmnts_pft(ielmc,NZ,NY,NX)+CanopyNonstructElements_pft(ielmc,NZ,NY,NX)),mass based CN ratio of leaf  
  allocate(this%histr_2D_tSOC_vr_col(beg_col:end_col,1:JZ))         !ORGC(1:JZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX), total soil C                                          
  allocate(this%histr_2D_LEAF_NODE_NO_ptc(beg_ptc:end_ptc,1:MaxNumBranches))        !NumOfLeaves_brch(NumOfMainBranch_pft(NZ,NY,NX),NZ,NY,NX), leaf NO
  allocate(this%histr_2D_RUB_ACTVN_ptc(beg_ptc:end_ptc,1:MaxNumBranches));      this%histr_2D_RUB_ACTVN_ptc=spval
  allocate(this%histr_2D_CO2_vr_col(beg_col:end_col,1:JZ))          !trc_solcl_vr(idg_CO2,1:JZ,NY,NX)
  allocate(this%histr_2D_CH4_vr_col(beg_col:end_col,1:JZ));    this%histr_2D_CH4_vr_col=spval
  allocate(this%histr_2D_O2_vr_col(beg_col:end_col,1:JZ))           !trc_solcl_vr(idg_O2,1:JZ,NY,NX)
  allocate(this%histr_2D_N2O_vr_col(beg_col:end_col,1:JZ))          !trc_solcl_vr(idg_N2O,1:JZ,NY,NX)
  allocate(this%histr_2D_NH3_vr_col(beg_col:end_col,1:JZ))          !trc_solcl_vr(idg_NH3,1:JZ,NY,NX)
  allocate(this%histr_2D_TEMP_vr_col(beg_col:end_col,1:JZ))         !TCS(1:JZ,NY,NX)
  allocate(this%histr_2D_vWATER_vr_col(beg_col:end_col,1:JZ))        !THETWZ(1:JZ,NY,NX)
  allocate(this%histr_2D_vICE_vr_col(beg_col:end_col,1:JZ))          !THETIZ(1:JZ,NY,NX)
  allocate(this%histr_2D_PSI_vr_col(beg_col:end_col,1:JZ))          !PSISM(1:JZ,NY,NX)+PSISO(1:JZ,NY,NX)
  allocate(this%histr_2D_cNH4t_vr_col(beg_col:end_col,1:JZ))        !(trc_solml_vr(ids_NH4,1:JZ,NY,NX)+trc_solml_vr(ids_NH4B,1:JZ,NY,NX) &
                                                               !+14.0*(trcx_solml(idx_NH4,1:JZ,NY,NX)+trcx_solml(idx_NH4B,1:JZ,NY,NX)))/SoilMicPMassLayer(1:JZ,NY,NX)
  allocate(this%histr_2D_cNO3t_vr_col(beg_col:end_col,1:JZ))        !(trc_solml_vr(ids_NO3,1:JZ,NY,NX)+trc_solml_vr(ids_NO3B,1:JZ,NY,NX) &
                                                               !+trc_solml_vr(ids_NO2,1,NY,NX)+trc_solml_vr(ids_NO2B,1,NY,NX))/SoilMicPMassLayer(1,NY,NX)
  allocate(this%histr_2D_cPO4_vr_col(beg_col:end_col,1:JZ))         !(trc_solml_vr(ids_H1PO4,1:JZ,NY,NX)+trc_solml_vr(ids_H1PO4B,1,NY,NX)+trc_solml_vr(ids_H2PO4,1,NY,NX)+trc_solml_vr(ids_H2PO4B,1,NY,NX))/VLWatMicP(1,NY,NX)
  allocate(this%histr_2D_cEXCH_P_vr_col(beg_col:end_col,1:JZ))      !31.0*(trcx_solml(idx_HPO4,1:JZ,NY,NX)+trcx_solml(idx_H2PO4,1:JZ,NY,NX)+trcx_solml(idx_HPO4B,1:JZ,NY,NX)+trcx_solml(idx_H2PO4B,1,NY,NX))/SoilMicPMassLayer(1,NY,NX)
  allocate(this%histr_2D_ECND_vr_col(beg_col:end_col,1:JZ))         !ECND(1:JZ,NY,NX)
  allocate(this%histr_2D_PSI_RT_vr_ptc(beg_ptc:end_ptc,1:JZ))       !PSIRoot_vr(1,1:JZ,NZ,NY,NX), root total water potential , MPa
  allocate(this%histr_2D_prtUP_NH4_vr_ptc(beg_ptc:end_ptc,1:JZ))       !(RootNutUptake_pvr(ids_NH4,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NH4,2,1:JZ,NZ,NY,NX) &
                                                               !+RootNutUptake_pvr(ids_NH4B,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NH4B,2,1:JZ,NZ,NY,NX))/AREA(3,1,NY,NX)
  allocate(this%histr_2D_prtUP_NO3_vr_ptc(beg_ptc:end_ptc,1:JZ))       !(RootNutUptake_pvr(ids_NO3,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NO3,2,1:JZ,NZ,NY,NX) &
                                                               !+RootNutUptake_pvr(ids_NO3B,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_NO3B,2,1:JZ,NZ,NY,NX))/AREA(3,1,NY,NX)
  allocate(this%histr_2D_prtUP_PO4_vr_ptc(beg_ptc:end_ptc,1:JZ))       !(RootNutUptake_pvr(ids_H2PO4,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_H2PO4,2,1:JZ,NZ,NY,NX) &
                                                               !+RootNutUptake_pvr(ids_H2PO4B,1,1:JZ,NZ,NY,NX)+RootNutUptake_pvr(ids_H2PO4B,2,1:JZ,NZ,NY,NX))/AREA(3,1,NY,NX)
  allocate(this%histr_2D_DNS_RT_vr_ptc(beg_ptc:end_ptc,1:JZ))       !RootLenthDensPerPopu_pvr(1,1:JZ,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)


  !-----------------------------------------------------------------------
  ! initialize history fields 
  !--------------------------------------------------------------------
  data1d_ptr => this%histr_1D_tFIRE_CO2_col(beg_col:end_col) 
  call hist_addfld1d(fname='tFIRE_CO2',units='gC m-2',avgflag='A',&
    long_name='total CO2 flux from fire',ptr_col=data1d_ptr)        

  data1d_ptr => this%histr_1D_tFIRE_CH4_col(beg_col:end_col)  
  call hist_addfld1d(fname='tFIRE_CH4',units='',avgflag='A', &
    long_name='total CH4 flux from fire', ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_cNH4_LITR_col(beg_col:end_col) 
  call hist_addfld1d(fname='cNH4_LITR',units='gN NH4/Mg litter',avgflag='A', &
    long_name='NH4 concentration in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_cNO3_LITR_col(beg_col:end_col)        
  call hist_addfld1d(fname='cNO3_LITR',units='gN NO3/Mg litter',avgflag='A',&
    long_name='NO3 concentration in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_HVST_C_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_HVST_C',units='gC/m2',avgflag='A',&
    long_name='Harvested C',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_HVST_N_col(beg_col:end_col)     
  call hist_addfld1d(fname='ECO_HVST_N',units='gN/m2',avgflag='A',&
    long_name='Harvested N',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_HVST_P_col(beg_col:end_col)   
  call hist_addfld1d(fname='ECO_HVST_P',units='gP/m2',avgflag='A',&
    long_name='Harvested P',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_NET_N_MIN_col(beg_col:end_col)  
  call hist_addfld1d(fname='NET_N_MIN',units='gN/m2/hr',avgflag='A',&
    long_name='total NH4 net mineralization (- means immobilization)',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SURF_tLITR_C_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SURF_tLITR_C_FLX',units='gC/m2',avgflag='A',&
    long_name='column integrated total litter C',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SURF_tLITR_N_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='SURF_tLITR_N_FLX',units='gN/m2',avgflag='A',&
    long_name='column integrated total litter N',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SURF_tLITR_P_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SURF_tLITR_P_FLX',units='gP/m2',avgflag='A',&
    long_name='column integrated total litter P',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_HUMUS_C_col(beg_col:end_col) 
  call hist_addfld1d(fname='HUMUS_C',units='gC/m2',avgflag='A',&
    long_name='colum integrated humus C',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_HUMUS_N_col(beg_col:end_col)         
  call hist_addfld1d(fname='HUMUS_N',units='gN/m2',avgflag='A',&
    long_name='colum integrated humus N',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_HUMUS_P_col(beg_col:end_col)       
  call hist_addfld1d(fname='HUMUS_P',units='gP/m2',avgflag='A',&
    long_name='colum integrated humus P',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_AMENDED_C_col(beg_col:end_col)   
  call hist_addfld1d(fname='AMENDED_C',units='gC/m2',avgflag='A',&
    long_name='total fertilizer C amendment',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_AMENDED_N_col(beg_col:end_col)   
  call hist_addfld1d(fname='AMENDED_N',units='gN/m2',avgflag='A',&
    long_name='total fertilizer N amendment',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_AMENDED_P_col(beg_col:end_col)        
  call hist_addfld1d(fname='AMENDED_P',units='gP/m2',avgflag='A',&
    long_name='total fertilizer P amendment',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tLITRf_C_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='tLITRf_C',units='gC/m2/hr',avgflag='A',&
    long_name='total litterfall C',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tLITRf_N_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='tLITRf_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total litterfall N',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tLITRf_P_FLX_col(beg_col:end_col) 
  call hist_addfld1d(fname='tLITRf_P',units='gP/m2/hr',avgflag='A',&
    long_name='total litterfall P',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tEXCH_PO4_col(beg_col:end_col)      
  call hist_addfld1d(fname='tEXCH_PO4',units='gP/m2',avgflag='A',&
    long_name='total bioavailable (exchangeable) mineral P',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUR_DOC_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUR_DOC_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='total surface DOC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUR_DON_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUR_DON_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total surface DON flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUR_DOP_FLX_col(beg_col:end_col) 
  call hist_addfld1d(fname='SUR_DOP_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total surface DOP flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUB_DOC_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUB_DOC_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='total subsurface DOC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUB_DON_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUB_DON_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total subsurface DON flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUB_DOP_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='SUB_DOP_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total subsurface DOP flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUR_DIC_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUR_DIC_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='total surface DIC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUR_DIN_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUR_DIN_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total surface DIN flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUR_DIP_FLX_col(beg_col:end_col)     
  call hist_addfld1d(fname='SUR_DIP_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total surface DIP flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUB_DIC_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUB_DIC_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='total subsurface DIC flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUB_DIN_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='SUB_DIN_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='landscape total subsurface DIN flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SUB_DIP_FLX_col(beg_col:end_col)     
  call hist_addfld1d(fname='SUB_DIP_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total subsurface DIP flux',ptr_col=data1d_ptr)      

  IF(salt_model)THEN
    data1d_ptr => this%histr_1D_tSALT_DISCHG_FLX_col(beg_col:end_col) 
    call hist_addfld1d(fname='tSALT_DISCHG_FLX',units='mol/m2/hr',avgflag='A',&
      long_name='total subsurface ion flux',ptr_col=data1d_ptr)      
  endif

  data1d_ptr => this%histr_1D_tPRECIP_P_col(beg_col:end_col)  
  call hist_addfld1d(fname='tPRECIP_P',units='gP/m2',avgflag='A',&
    long_name='column integrated total soil precipited P',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tMICRO_C_col(beg_col:end_col)    
  call hist_addfld1d(fname='tMICRO_C',units='gC/m2',avgflag='A', &
    long_name='total micriobial C',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tMICRO_N_col(beg_col:end_col)    
  call hist_addfld1d(fname='tMICRO_N',units='gN/m2',avgflag='A', &
    long_name='total micriobial N',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tMICRO_P_col(beg_col:end_col)    
  call hist_addfld1d(fname='tMICRO_P',units='gP/m2',avgflag='A', &
    long_name='total micriobial P',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_PO4_FIRE_col(beg_col:end_col)  
  call hist_addfld1d(fname='PO4_FIRE',units='gP/m2/hr',avgflag='A',&
    long_name='total PO4 flux from fire',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_cPO4_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='cPO4_LITR',units='gP/Mg litr',avgflag='A',&
    long_name='PO4 concentration in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_cEXCH_P_LITR_col(beg_col:end_col)     
  call hist_addfld1d(fname='cEXCH_P_LITR',units='gP/Mg litr',avgflag='A',&
    long_name='concentration of exchangeable inorganic P in litterr',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_NET_P_MIN_col(beg_col:end_col)    
  call hist_addfld1d(fname='NET_P_MIN',units='gP/m2/hr',avgflag='A',&
    long_name='total inorganic P net mineralization (-ve) or immobilization (+ve)',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_RADN_col(beg_col:end_col)         
  call hist_addfld1d(fname='RADN',units='MJ/day',avgflag='A',&
    long_name='*total daily solar radiation',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_AIR_TEMP_col(beg_col:end_col)    
  call hist_addfld1d(fname='TMAX_AIR',units='oC',avgflag='X',&
    long_name='daily maximum air temperature',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_AIR_TEMP_col(beg_col:end_col)     
  call hist_addfld1d(fname='TMIN_AIR',units='oC',avgflag='M',&
    long_name='daily minimum air temperature',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_HUM_col(beg_col:end_col)       
  call hist_addfld1d(fname='HMAX_AIR',units='kPa',avgflag='X',&
    long_name='daily maximum vapor pressure',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_HUM_col(beg_col:end_col)    
  call hist_addfld1d(fname='HMIN_AIR',units='kPa',avgflag='M',&
    long_name='daily maximum vapor pressure',ptr_col=data1d_ptr)      
    
  data1d_ptr => this%histr_1D_PSI_SURF_col(beg_col:end_col)    
  call hist_addfld1d(fname='PSI_SURF',units='MPa',avgflag='A',&
    long_name='soil micropore matric water potential',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SURF_ELEV_col(beg_col:end_col)  
  call hist_addfld1d(fname='SURF_ELEV',units='m',avgflag='A',&
    long_name='Surface elevation, including litter layer',ptr_col=data1d_ptr)      
    
  data1d_ptr => this%histr_1D_tNH4X_col(beg_col:end_col)     
  call hist_addfld1d(fname='tNH4',units='gN/m2',avgflag='A', &
    long_name='column integrated NH4+NH3',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tNO3_col(beg_col:end_col)          
  call hist_addfld1d(fname='tNO3',units='gN/m2',avgflag='A',&
    long_name='Column integrated NO3+NO2 content',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_TEMP_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='TEMP_LITR',units='oC',avgflag='A',&
    long_name='Litter layer temperature',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_TEMP_SNOW_col(beg_col:end_col)    
  call hist_addfld1d(fname='TEMP_SNOW',units='oC',avgflag='A',&
    long_name='First snow layer temperature',ptr_col=data1d_ptr)      
    
  data1d_ptr => this%histr_1D_CO2_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='CO2_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='total soil CO2 flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_OMC_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='OMC_LITR',units='gC/m2',avgflag='A',&
    long_name='total litter residual C',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ATM_CO2_col(beg_col:end_col)   
  call hist_addfld1d(fname='ATM_CO2',units='umol/mol',avgflag='A',&
    long_name='Atmospheric CO2 concentration',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_NBP_col(beg_col:end_col)   
  call hist_addfld1d(fname='NBP',units='gC/m2/hr',avgflag='A',&
    long_name='Net biosphere productivity',ptr_col=data1d_ptr)      
    
  data1d_ptr => this%histr_1D_ECO_LAI_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_LAI',units='m2/m2',avgflag='A',&
    long_name='ecosystem LAI',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_GPP_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_GPP',units='gC/m2/hr',avgflag='A',&
    long_name='ecosystem GPP',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_RA_col(beg_col:end_col)       
  call hist_addfld1d(fname='ECO_RA',units='gC/m2/hr',avgflag='A',&
    long_name='ecosystem autotrophic respiration',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_NPP_col(beg_col:end_col)      
  call hist_addfld1d(fname='ECO_NPP',units='gC/m2/hr',avgflag='A',&
    long_name='ecosystem NPP',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_RH_col(beg_col:end_col)      
  call hist_addfld1d(fname='ECO_RH',units='gC/m2/hr',avgflag='A',&
    long_name='ecosystem heterotrophic respiration',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tDIC_col(beg_col:end_col)       
  call hist_addfld1d(fname='tDIC',units='gC/m2',avgflag='A',&
    long_name='column integrated total soil DIC: CO2+CH4',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tSTANDING_DEAD_C_col(beg_col:end_col)     
  call hist_addfld1d(fname='tSTANDING_DEAD_C',units='gC/m2',avgflag='A',&
    long_name='total standing dead C',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tSTANDING_DEAD_N_col(beg_col:end_col)     
  call hist_addfld1d(fname='tSTANDING_DEAD_N',units='gN/m2',avgflag='A',&
    long_name='total standing dead N',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tSTANDING_DEAD_P_col(beg_col:end_col)     
  call hist_addfld1d(fname='tSTANDING_DEAD_P',units='gP/m2',avgflag='A',&
    long_name='total standing dead P',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_tPRECN_col(beg_col:end_col)      
  call hist_addfld1d(fname='tPRECN',units='mm/m2',avgflag='A',&
    long_name='total precipitation, including irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ET_col(beg_col:end_col)  
  call hist_addfld1d(fname='ET',units='mm/m2',avgflag='A',&
    long_name='total evapotranspiration',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_N2O_LITR_col(beg_col:end_col)      
  call hist_addfld1d(fname='N2O_LITR',units='g/m3',avgflag='A',&
    long_name='N2O solute concentration in soil micropres',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_NH3_LITR_col(beg_col:end_col)   
  call hist_addfld1d(fname='NH3_LITR',units='g/m3',avgflag='A',&
    long_name='NH3 solute concentration in soil micropres',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SOL_RADN_col(beg_col:end_col)      
  call hist_addfld1d(fname='SOL_RADN',units='W/m2',avgflag='A',&
    long_name='shortwave radiation in solar beam',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_AIR_TEMP_col(beg_col:end_col)      
  call hist_addfld1d(fname='AIR_TEMP',units='oC',avgflag='A',&
    long_name='air temperature',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_HUM_col(beg_col:end_col)    
  call hist_addfld1d(fname='HUM',units='kPa',avgflag='A',&
    long_name='vapor pressure',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_WIND_col(beg_col:end_col)   
  call hist_addfld1d(fname='WIND',units='m/s',avgflag='A',&
    long_name='wind spee',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_PREC_col(beg_col:end_col)   
  call hist_addfld1d(fname='PREC',units='mm H2O/m2',avgflag='A',&
    long_name='Total precipitation, excluding irrigation',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SOIL_RN_col(beg_col:end_col)   
  call hist_addfld1d(fname='SOIL_RN',units='W/m2',avgflag='A',&
    long_name='total net radiation at ground surface',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SOIL_LE_col(beg_col:end_col)      
  call hist_addfld1d(fname='SOIL_LE',units='W/m2',avgflag='A',&
    long_name='total latent heat flux at ground surface',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SOIL_H_col(beg_col:end_col)      
  call hist_addfld1d(fname='SOIL_H',units='W/m2',avgflag='A',&
    long_name='total sensible heat flux at ground surface',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SOIL_G_col(beg_col:end_col)  
  call hist_addfld1d(fname='SOIL_G',units='W/m2',avgflag='A',&
    long_name='*total heat flux into ground surface',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_RN_col(beg_col:end_col)     
  call hist_addfld1d(fname='ECO_RN',units='W/m2',avgflag='A',&
    long_name='ecosystem net radiation',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_LE_col(beg_col:end_col)        
  call hist_addfld1d(fname='ECO_LE',units='W/m2',avgflag='A',&
    long_name='ecosystem latent heat flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_Eco_Heat_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_H',units='W/m2',avgflag='A',&
    long_name='ecosystem sensible heat flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_G_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_G',units='W/m2',avgflag='A',&
    long_name='ecosystem storage heat flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_O2_LITR_col(beg_col:end_col)      
  call hist_addfld1d(fname='O2_LITR',units='g/m3',avgflag='A',&
    long_name='O2 solute concentration in litter layer',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_MIN_LWP_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='MIN_LWP',units='MPa',avgflag='A',&
    long_name='minimum daily canopy water potential',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_SOIL_CO2_FLX_col(beg_col:end_col)
  call hist_addfld1d(fname='SOIL_CO2_FLX',units='umol/m2/s',avgflag='A',&
    long_name='soil CO2 flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ECO_CO2_FLX_col(beg_col:end_col)  
  call hist_addfld1d(fname='ECO_NEE_CO2',units='umol/m2/s',avgflag='A',&
    long_name='ecosystem net CO2 exchange',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_CH4_FLX_col(beg_col:end_col)     
  call hist_addfld1d(fname='CH4_FLX',units='umol/m2/s',avgflag='A',&
    long_name='soil CH4 flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_O2_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='O2_FLX',units='umol/m2/s',avgflag='A',&
    long_name='soil O2 flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_CO2_LITR_col(beg_col:end_col)      
  call hist_addfld1d(fname='CO2_LITR',units='g/m3',avgflag='A',&
    long_name='CO2 solute concentration in litter',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_EVAPN_col(beg_col:end_col)      
  call hist_addfld1d(fname='EVAPN',units='mm H2O/m2/hr',avgflag='A',&
    long_name='total evaporation at ground surface',ptr_col=data1d_ptr)      
    
  data1d_ptr => this%histr_1D_tSWC_col(beg_col:end_col)  
  call hist_addfld1d(fname='tSWC',units='mmH2O/m2',avgflag='A', &
    long_name='column integrated water content',ptr_col=data1d_ptr)        

  data1d_ptr => this%histr_1D_SNOWPACK_col(beg_col:end_col)      
  call hist_addfld1d(fname='SNOWPACK',units='mmH2O/m2',&
    avgflag='A',long_name='total water equivalent snow',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SURF_WTR_col(beg_col:end_col)   
  call hist_addfld1d(fname='SURF_WTR',units='m3/m3',avgflag='A',&
    long_name='Volumetric water content in surface litter layer',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SURF_ICE_col(beg_col:end_col)    
  call hist_addfld1d(fname='SURF_ICE',units='m3/m3',avgflag='A',&
    long_name='Volumetric ice content in surface litter layer',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_ACTV_LYR_col(beg_col:end_col)    
  call hist_addfld1d(fname='ACTV_LYR',units='m',avgflag='A',&
    long_name='active layer depth',ptr_col=data1d_ptr)        

  data1d_ptr => this%histr_1D_WTR_TBL_col(beg_col:end_col)      
  call hist_addfld1d(fname='WTR_TBL',units='m',avgflag='A',&
    long_name='internal water table depth',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_sN2O_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='sN2O_FLX',units='g/m2/hr',avgflag='A',&
    long_name='*soil N2O flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_sN2G_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='sN2G_FLX',units='g/m2/hr',&
    avgflag='A',long_name='soil N2 flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_sNH3_FLX_col(beg_col:end_col)       
  call hist_addfld1d(fname='sNH3_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='soil NH3 flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_RUNOFF_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='RUNOFF_FLX',units='mmH2O/m2/hr',avgflag='A',&
    long_name='landscape runoff from surface water',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_SEDIMENT_FLX_col(beg_col:end_col)      
  call hist_addfld1d(fname='SEDIMENT_FLX',units='kg/m2/hr',avgflag='A',&
    long_name='total sediment subsurface flux',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_DISCHG_FLX_col(beg_col:end_col)   
  call hist_addfld1d(fname='DISCHG_FLX',units='mmH2O/m2/hr',avgflag='A',&
    long_name='landscape water discharge',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_LEAF_PC_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='LEAF_PC',units='gP/gC',avgflag='I',&
    long_name='mass based leaf PC ratio',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_CAN_RN_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CAN_RN',units='W/m2',avgflag='A',&
    long_name='canopy net radiation',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_CAN_LE_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CAN_LE',units='W/m2',avgflag='A',&
    long_name='canopy latent heat flux (<0 to ATM)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_CAN_H_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='CAN_H',units='W/m2',avgflag='A',&
    long_name='canopy sensible heat flux',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_CAN_G_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='CAN_G',units='W/m2',avgflag='A',&
    long_name='canopy storage heat flux',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_CAN_TEMP_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_TEMP',units='oC',avgflag='A',&
    long_name='canopy temperature',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_TEMP_FN_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='TEMP_FN',units='none',avgflag='A',&
    long_name='canopy temperature growth function/stress,[0-1]',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_CAN_CO2_FLX_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='CAN_CO2_FLX',units='umol/m2/hr',avgflag='A',&
    long_name='canopy net CO2 exchange',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_CAN_GPP_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_GPP',units='gC/m2',avgflag='A',&
    long_name='cumulative total gross CO2 fixation',ptr_patch=data1d_ptr)      
  
  data1d_ptr => this%histr_1D_CAN_RA_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_RA',units='gC/m2/hr',avgflag='A',&
    long_name='total autotrophic respiration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_cTNC_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='cTNC',units='gC/gC',avgflag='A',&
    long_name='canopy nonstructural C concentration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_cTNN_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='cTNN',units='gN/gC',avgflag='A',&
    long_name='canopy nonstructural N concentration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_cTNP_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='cTNP',units='gP/gC',avgflag='A',&
    long_name='canopy nonstructural P concentration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STOML_RSC_CO2_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='STOML_RSC_CO2',units='s/m',avgflag='A',&
    long_name='canopy stomatal resistance for CO2',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_BLYR_RSC_CO2_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='BLYR_RSC_CO2',units='s/m',avgflag='A',&
    long_name='canopy boundary layer resistance for CO2',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_CAN_CO2_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='CAN_CO2',units='umol/mol',avgflag='A',&
    long_name='canopy gaesous CO2 concentration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_LAI_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='LAI',units='m2/m2',avgflag='A',&
    long_name='plant leaf area, including stalk',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_PSI_CAN_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='PSI_CAN',units='MPa',avgflag='A',&
    long_name='canopy total water potential',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_TURG_CAN_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='TURG_CAN',units='MPa',avgflag='A',&
    long_name='canopy turgor water potential',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STOM_RSC_H2O_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='STOM_RSC',units='s/m',avgflag='A',&
    long_name='canopy stomatal resistance for H2O',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_BLYR_RSC_H2O_ptc(beg_ptc:end_ptc)
  call hist_addfld1d(fname='BLYR_RSC_H2O',units='s/m',avgflag='A',&
    long_name='canopy boundary layer resistance for H2O',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_TRANSPN_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='TRANSPN',units='mmH2O/m2/h',avgflag='A',&
    long_name='canopy transpiration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_O2_STRESS_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='O2_STRESS',units='none',avgflag='A',&
    long_name='plant O2 stress indicator [0-1]',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_NH4_UPTK_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='NH4_UPTK_FLX',units='gN/m2/hr',&
    avgflag='A',long_name='total root uptake of NH4',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_NO3_UPTK_FLX_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='NO3_UPTK_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total root uptake of NO3',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_N2_FIXN_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='N2_FIXN_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total root N2 fixation',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_cNH3_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='cNH3_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='*canopy NH3 flux',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_PO4_UPTK_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='PO4_UPTK_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total root uptake of PO4',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_SHOOT_C_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='SHOOT_C',units='gC/m2',avgflag='A',&
    long_name='canopy shoot C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_Plant_C_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='PLANT_C',units='gC/m2',avgflag='A',&
    long_name='plant C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_LEAF_C_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='LEAF_C',units='gC/m2',avgflag='A',&
    long_name='canopy leaf C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_Petiole_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='PETIOLE_C',units='gC/m2',avgflag='A',&
    long_name='canopy sheath C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STALK_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='STALK_C',units='gC/m2',avgflag='A',&
    long_name='canopy stalk C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_RESERVE_C_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='RESERVE_C',units='gC/m2',avgflag='A',&
    long_name='canopy reserve C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_HUSK_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='HUSK_C',units='gC/m2',avgflag='A',&
    long_name='canopy husk C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_GRAIN_C_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='GRAIN_C',units='gC/m2',avgflag='A',&
    long_name='canopy grain C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_ROOT_C_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='ROOT_C',units='gC/m2',avgflag='A',&
    long_name='plant root C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_NODULE_C_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='NODULE_C',units='gC/m2',avgflag='A',&
    long_name='root total nodule C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STORED_C_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='STORED_C',units='gC/m2',avgflag='A',&
    long_name='plant stored nonstructural C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_GRAIN_NO_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='GRAIN_NO',units='1/m2',avgflag='A',&
    long_name='canopy grain number',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_LAIb_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='LAIb',units='m2/m2',avgflag='A',&
    long_name=' total plant leaf area, exclude stalk',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_EXUD_C_FLX_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='EXUD_C_FLX',units='gC/m2',avgflag='A',&
    long_name='total net root C uptake (+ve) - exudation',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_LITRf_C_FLX_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='LITRf_C',units='gC/m2/hr',avgflag='A',&
    long_name='total plant litterfall C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_SURF_LITRf_C_FLX_ptc(beg_ptc:end_ptc) 
  call hist_addfld1d(fname='SURF_LITRf_C_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='plant litterfall C to the soil surface',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_AUTO_RESP_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='AUTO_RESP',units='gC/m2',avgflag='A',&
    long_name='cumulative plant autotrophic respiration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_ABV_GRD_RESP_FLX_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='ABV_GRD_RESP',units='gC/m2/hr',&
    avgflag='A',long_name='plant shoot autotrophic respiration',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_HVST_C_FLX_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='HVST_C_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='plant C harvest',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_PLANT_BALANCE_C_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='PLANT_BALANCE_C',units='gC/m2',avgflag='A',&
    long_name='plant C balance',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STANDING_DEAD_C_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='STANDING_DEAD_C',units='gC/m2',avgflag='A',&
    long_name='pft Standing dead C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_FIREp_CO2_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='FIREp_CO2_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='plant CO2 from fire',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_FIREp_CH4_FLX_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='FIREp_CH4_FLX',units='gC/m2/hr',avgflag='A',&
    long_name='plant CH4 emission from fire',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_NPP_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='NPP',units='gC/m2',avgflag='A',&
    long_name='Cumulative net primary productivity',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_CAN_HT_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='CAN_HT',units='m',avgflag='A',&
    long_name='Canopy height',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_POPN_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='POPN',units='1/m2',avgflag='A',&
    long_name='Plant population',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_tTRANSPN_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='tTRANSPN',units='mmH2O/m2/hr',avgflag='A',&
    long_name='Total evapotranspiration (>0 into atmosphere)',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_WTR_STRESS_ptc(beg_ptc:end_ptc)    !HoursCanopyPSITooLow(NZ,NY,NX)
  call hist_addfld1d(fname='WTR_STRESS',units='hr',avgflag='A',&
    long_name='canopy plant water stress indicator: number of hours PSICanopy_pft(< PSILY',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_OXY_STRESS_ptc(beg_ptc:end_ptc)    !OSTR(NZ,NY,NX)
  call hist_addfld1d(fname='OXY_STRESS',units='none',avgflag='A',&
    long_name='plant O2 stress indicator',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_SHOOT_N_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='SHOOT_N',units='gN/m2',avgflag='A',&
    long_name='canopy shoot N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_Plant_N_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='PLANT_N',units='gN/m2',avgflag='A',&
    long_name='plant N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_LEAF_N_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='LEAF_N',units='gN/m2',avgflag='A',&
    long_name='Canopy leaf N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_Petiole_N_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='Petiole_N',units='gN/m2',avgflag='A',&
    long_name='canopy sheath N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STALK_N_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='STALK_N',units='gN/m2',avgflag='A',&
    long_name='Canopy stalk N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_RESERVE_N_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='RESERVE_N',units='gN/m2',avgflag='A',&
    long_name='Canopy reserve N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_HUSK_N_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='HUSK_N',units='gN/m2',avgflag='A',&
    long_name='Canopy husk N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_GRAIN_N_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='GRAIN_N',units='gN/m2',avgflag='A',&
    long_name='Canopy grain C',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_ROOT_N_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='ROOT_N',units='gN/m2',avgflag='A',&
    long_name='gN/m2',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_NODULE_N_ptc(beg_ptc:end_ptc)        
  call hist_addfld1d(fname='NODULE_N',units='gN/m2',avgflag='A',&
    long_name='root total nodule N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STORED_N_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='STORED_N',units='gN/m2',avgflag='A',&
    long_name='plant stored nonstructural N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_EXUD_N_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='EXUD_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total net root N uptake (+ve) - exudation',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_LITRf_N_FLX_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='LITRf_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total plant litterfall N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_TL_N_FIXED_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='TL_N_FIXED_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total plant N2 fixation',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_HVST_N_FLX_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='HVST_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='plant N harvest',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_NH3can_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='NH3can_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total canopy NH3 flux',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_PLANT_BALANCE_N_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='PLANT_BALANCE_N',units='gC/m2',avgflag='A',&
    long_name='plant N balance',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STANDING_DEAD_N_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='STANDING_DEAD_N',units='gN/m2',avgflag='A',&
    long_name='pft standing dead N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_FIREp_N_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='FIREp_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='plant N emission from fire',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_SURF_LITRf_N_FLX_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='SURF_LITRf_N_FLX',units='gN/m2/hr',avgflag='A',&
    long_name='total surface litterfall N',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_SHOOT_P_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='SHOOT_P',units='gP/m2',avgflag='A',&
    long_name='canopy shoot P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_Plant_P_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='PLANT_P',units='gP/m2',avgflag='A',&
    long_name='plant P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_LEAF_P_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='LEAF_P',units='gP/m2',avgflag='A',&
    long_name='Canopy leaf P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_Petiole_P_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='Petiole_P',units='gP/m2',avgflag='A',&
    long_name='canopy sheath P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STALK_P_ptc(beg_ptc:end_ptc)      
  call hist_addfld1d(fname='STALK_P',units='gP/m2',avgflag='A',&
    long_name='Plant stalk P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_RESERVE_P_ptc(beg_ptc:end_ptc)  
  call hist_addfld1d(fname='RESERVE_P',units='gP/m2',avgflag='A',&
    long_name='Plant reserve P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_HUSK_P_ptc(beg_ptc:end_ptc)        
  call hist_addfld1d(fname='HUSK_P',units='gP/m2',avgflag='A',&
    long_name='Husk P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_GRAIN_P_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='GRAIN_P',units='gP/m2',avgflag='A',&
    long_name='canopy grain P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_ROOT_P_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='ROOT_P',units='gP/m2',avgflag='A',&
    long_name='plant root P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_NODULE_P_ptc(beg_ptc:end_ptc)       
  call hist_addfld1d(fname='NODULE_P',units='gP/m2',avgflag='A',&
    long_name='root total nodule P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STORED_P_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='STORED_P',units='gP/m2',avgflag='A',&
    long_name='plant stored nonstructural P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_EXUD_P_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='EXUD_P_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total net root P uptake (+ve) - exudation',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_LITRf_P_FLX_ptc(beg_ptc:end_ptc)     
  call hist_addfld1d(fname='LITRf_P_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='total plant litterfall P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_HVST_P_FLX_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='HVST_P_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='Plant P harvest',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_PLANT_BALANCE_P_ptc(beg_ptc:end_ptc)    
  call hist_addfld1d(fname='PLANT_BALANCE_P',units='gP/m2',avgflag='A',&
    long_name='plant P balance',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_STANDING_DEAD_P_ptc(beg_ptc:end_ptc)   
  call hist_addfld1d(fname='STANDING_DEAD_P',units='gP/m2',avgflag='A',&
    long_name='pft Standing dead P',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_FIREp_P_FLX_ptc(beg_ptc:end_ptc)             
  call hist_addfld1d(fname='FIREp_P_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='plant PO4 emission from fire',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_SURF_LITRf_P_FLX_ptc(beg_ptc:end_ptc)         !SurfLitrfallChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
  call hist_addfld1d(fname='SURF_LITRf_P_FLX',units='gP/m2/hr',avgflag='A',&
    long_name='plant litterfall P to the soil surface',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_BRANCH_NO_ptc(beg_ptc:end_ptc)            !NumOfBranches_pft(NZ,NY,NX)
  call hist_addfld1d(fname='BRANCH_NO',units='none',avgflag='I',&
    long_name='Plant branch number',ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_Growth_Stage_ptc(beg_ptc:end_ptc)  !plant development stage, integer, 0-10, planting, emergence, floral_init, jointing, 
                                                               !elongation, heading, anthesis, seed_fill, see_no_set, seed_mass_set, end_seed_fill
  call hist_addfld1d(fname='Growth_Stage',units='none',avgflag='I',&
    long_name='plant development stage, integer, 0-planting, 1-emergence, 2-floral_init, 3-jointing,'// &
    '4-elongation, 5-heading, 6-anthesis, 7-seed_fill, 8-see_no_set, 9-seed_mass_set, 10-end_seed_fill',&
    ptr_patch=data1d_ptr)      

  data1d_ptr => this%histr_1D_LEAF_NC_ptc(beg_ptc:end_ptc)            
  call hist_addfld1d(fname='LEAF_NC',units='gN/gC',avgflag='A',&
    long_name='mass based plant leaf NC ratio',ptr_patch=data1d_ptr)      

  data2d_ptr => this%histr_2D_tSOC_vr_col(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='tSOC_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='*Vertically resolved total soil organic C',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_LEAF_NODE_NO_ptc(beg_ptc:end_ptc,1:MaxNumBranches)        !NumOfLeaves_brch(NumOfMainBranch_pft(NZ,NY,NX),NZ,NY,NX), leaf NO
  call hist_addfld2d(fname='LEAF_NODE_NO',units='none',type2d='nbranches',avgflag='I',&
    long_name='Leaf number',ptr_patch=data2d_ptr)      

  data2d_ptr => this%histr_2D_RUB_ACTVN_ptc(beg_ptc:end_ptc,1:MaxNumBranches)      !RubiscoActivity_brpft(NumOfMainBranch_pft(NZ,NY,NX),NZ,NY,NX), branch down-regulation of CO2 fixation
  call hist_addfld2d(fname='RUB_ACTVN',units='none',type2d='nbranches',avgflag='A',&
    long_name='branch rubisco activity for CO2 fixation, 0-1',ptr_patch=data2d_ptr)      

  data2d_ptr => this%histr_2D_CO2_vr_col(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_CO2,1:JZ,NY,NX)
  call hist_addfld2d(fname='CO2_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='solute concentration of CO2 in soil micropre',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_CH4_vr_col(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_CH4,1:JZ,NY,NX)
  call hist_addfld2d(fname='CH4_vr',units='gC/m3',type2d='levsoi',avgflag='A',&
    long_name='solute concentration of CH4 in soil micropre',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_O2_vr_col(beg_col:end_col,1:JZ)           !trc_solcl_vr(idg_O2,1:JZ,NY,NX)
  call hist_addfld2d(fname='O2_vr',units='g/m3',type2d='levsoi',avgflag='A',&
    long_name='solute concentration of O2 in soil micropre',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_N2O_vr_col(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_N2O,1:JZ,NY,NX)
  call hist_addfld2d(fname='N2O_vr',units='g/m3',type2d='levsoi',avgflag='A',&
    long_name='solute concentration of N2O in soil micropre',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_NH3_vr_col(beg_col:end_col,1:JZ)          !trc_solcl_vr(idg_NH3,1:JZ,NY,NX)
  call hist_addfld2d(fname='NH3_vr',units='g/m3',type2d='levsoi',avgflag='A',&
    long_name='solute concentration of NH3 in soil micropre',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_TEMP_vr_col(beg_col:end_col,1:JZ)         !TCS(1:JZ,NY,NX)
  call hist_addfld2d(fname='TEMP_vr',units='oC',type2d='levsoi',avgflag='A',&
    long_name='soil temperature profile',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_vWATER_vr_col(beg_col:end_col,1:JZ)        !THETWZ(1:JZ,NY,NX)
  call hist_addfld2d(fname='vWATER_vr',units='m3/m3',type2d='levsoi',avgflag='A',&
    long_name='volumetric soil water content',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_vICE_vr_col(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='vICE_vr',units='m3/m3',type2d='levsoi',avgflag='A',&
    long_name='volumetric soil ice content',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_PSI_vr_col(beg_col:end_col,1:JZ)         
  call hist_addfld2d(fname='PSI_vr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='soil matric pressure+osmotic pressure',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_cNH4t_vr_col(beg_col:end_col,1:JZ)       
  call hist_addfld2d(fname='cNH4t_vr',units='gN/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='soil NH4x concentration',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_cNO3t_vr_col(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='cNO3t_vr',units='gN/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='Soil NO3+NO2 concentration',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_cPO4_vr_col(beg_col:end_col,1:JZ)        
  call hist_addfld2d(fname='cPO4_vr',units='gP/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='soil dissolved PO4 concentration',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_cEXCH_P_vr_col(beg_col:end_col,1:JZ)     
  call hist_addfld2d(fname='cEXCH_P_vr',units='gP/Mg soil',type2d='levsoi',avgflag='A',&
    long_name='total exchangeable soil PO4 concentration',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_TEMP_vr_col(beg_col:end_col,1:JZ)  
  call hist_addfld2d(fname='TMAX_SOIL_vr',units='oC',type2d='levsoi',avgflag='X',&
    long_name='Soil maximum temperature profile',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_TEMP_vr_col(beg_col:end_col,1:JZ)  
  call hist_addfld2d(fname='TMIN_SOIL_vr',units='oC',type2d='levsoi',avgflag='M',&
    long_name='Soil minimum temperature profile',ptr_col=data2d_ptr)      

  data1d_ptr => this%histr_1D_TEMP_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='TMAX_LITR',units='oC',avgflag='X',&
    long_name='Litter maximum temperature',ptr_col=data1d_ptr)      

  data1d_ptr => this%histr_1D_TEMP_LITR_col(beg_col:end_col)    
  call hist_addfld1d(fname='TMIN_LITR',units='oC',avgflag='M',&
    long_name='Litter minimum temperature',ptr_col=data1d_ptr)      

  data2d_ptr => this%histr_2D_ECND_vr_col(beg_col:end_col,1:JZ)     
  call hist_addfld2d(fname='ECND_vr',units='dS m-1',type2d='levsoi',avgflag='A',&
    long_name='electrical conductivity',ptr_col=data2d_ptr)      

  data2d_ptr => this%histr_2D_PSI_RT_vr_ptc(beg_ptc:end_ptc,1:JZ)  
  call hist_addfld2d(fname='PSI_RT_vr',units='MPa',type2d='levsoi',avgflag='A',&
    long_name='root total water potential',ptr_patch=data2d_ptr)      

  data2d_ptr => this%histr_2D_prtUP_NH4_vr_ptc(beg_ptc:end_ptc,1:JZ) 
  call hist_addfld2d(fname='prtUP_NH4_vr',units='gN/m3/hr',type2d='levsoi',avgflag='A',&
    long_name='root uptake of NH4',ptr_patch=data2d_ptr)      

  data2d_ptr => this%histr_2D_prtUP_NO3_vr_ptc(beg_ptc:end_ptc,1:JZ)      
  call hist_addfld2d(fname='prtUP_NO3_vr',units='gN/m3/hr',type2d='levsoi',&
    avgflag='A',long_name='root uptake of NO3',ptr_patch=data2d_ptr)      

  data2d_ptr => this%histr_2D_prtUP_PO4_vr_ptc(beg_ptc:end_ptc,1:JZ)     
  call hist_addfld2d(fname='prtUP_PO4_vr',units='gP/m3/hr',type2d='levsoi',avgflag='A',&
    long_name='root uptake of PO4',ptr_patch=data2d_ptr)      

  data2d_ptr => this%histr_2D_DNS_RT_vr_ptc(beg_ptc:end_ptc,1:JZ)       
  call hist_addfld2d(fname='DNS_RT_vr',units='m/m3',type2d='levsoi',avgflag='A',&
    long_name='root layer length density',ptr_patch=data2d_ptr)      

  end subroutine init_hist_data

!----------------------------------------------------------------------
  subroutine hist_update(this,bounds)
  implicit none
  class(histdata_type) :: this
  type(bounds_type), intent(in) :: bounds
  integer :: ncol,nptc
  integer :: L,NZ,NY,NX,KN,NB
  
  DO NX=bounds%NHW,bounds%NHE   
    DO NY=bounds%NVN,bounds%NVS
      ncol=get_col(NY,NX)
      this%histr_1D_tFIRE_CO2_col(ncol) =  CO2byFire_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tFIRE_CH4_col(ncol) =  CH4byFire_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_cNH4_LITR_col(ncol) =  safe_adb(trc_solml_vr(ids_NH4,0,NY,NX)+natomw*trcx_solml(idx_NH4,0,NY,NX),SoilMicPMassLayer(0,NY,NX))
      this%histr_1D_cNO3_LITR_col(ncol) =  safe_adb(trc_solml_vr(ids_NO3,0,NY,NX)+trc_solml_vr(ids_NO2,0,NY,NX),SoilMicPMassLayer(0,NY,NX))
      this%histr_1D_ECO_HVST_N_col(ncol)=  EcoHavstElmnt_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_NET_N_MIN_col(ncol) = -NetNH4Mineralize_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_SURF_tLITR_P_FLX_col(ncol) =  URSDP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_HUMUS_C_col(ncol)     = UORGC(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_HUMUS_N_col(ncol)     = UORGN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_HUMUS_P_col(ncol)     = UORGP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_AMENDED_P_col(ncol)   = FerPFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tLITRf_C_FLX_col(ncol)= LiterfalOrgC_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tLITRf_N_FLX_col(ncol)= LiterfalOrgN_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tLITRf_P_FLX_col(ncol)= LiterfalOrgP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tEXCH_PO4_col(ncol)   = UPO4(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_SUR_DOP_FLX_col(ncol)= HydroDOPFlx_col(NY,NX)/TAREA
      this%histr_1D_SUB_DOP_FLX_col(ncol)= HDOPD(NY,NX)/TAREA
      this%histr_1D_SUR_DIP_FLX_col(ncol)= HydroDIPFlx_col(NY,NX)/TAREA
      this%histr_1D_SUB_DIP_FLX_col(ncol)= HDIPD(NY,NX)/TAREA  
      this%histr_1D_SUR_DON_FLX_col(ncol) = HydroDONFlx_col(NY,NX)/TAREA
      this%histr_1D_SUB_DON_FLX_col(ncol) = HDOND(NY,NX)/TAREA
      this%histr_1D_tSALT_DISCHG_FLX_col(ncol)= HydroIonFlx_col(NY,NX)/TAREA
      this%histr_1D_SUR_DIN_FLX_col(ncol) = HydroDINFlx_col(NY,NX)/TAREA
      this%histr_1D_SUB_DIN_FLX_col(ncol) = HDIND(NY,NX)/TAREA
      this%histr_1D_SUR_DOC_FLX_col(ncol) = HDOCQ(NY,NX)/TAREA
      this%histr_1D_SUB_DOC_FLX_col(ncol) = HDOCD(NY,NX)/TAREA
      this%histr_1D_SUR_DIC_FLX_col(ncol) = HDICQ(NY,NX)/TAREA
      this%histr_1D_SUB_DIC_FLX_col(ncol) = HDICD(NY,NX)/TAREA
      this%histr_1D_SUR_DIP_FLX_col(ncol)  = HydroDIPFlx_col(NY,NX)/TAREA
      this%histr_1D_tPRECIP_P_col(ncol)    = UPP4(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tMICRO_P_col(ncol)     = TOPT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_PO4_FIRE_col(ncol)    = PO4byFire_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_cPO4_LITR_col(ncol)   = safe_adb(trc_solml_vr(ids_H2PO4,0,NY,NX),SoilMicPMassLayer(0,NY,NX))
      this%histr_1D_cEXCH_P_LITR_col(ncol)=  patomw*safe_adb(trcx_solml(idx_HPO4,0,NY,NX)+trcx_solml(idx_H2PO4,0,NY,NX),SoilMicPMassLayer(0,NY,NX))
      this%histr_1D_ECO_HVST_P_col(ncol)  = EcoHavstElmnt_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_NET_P_MIN_col(ncol)   =  -NetPO4Mineralize_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_RADN_col(ncol)        = TRAD(NY,NX)
      this%histr_1D_PSI_SURF_col(ncol)    = PSISoilMatricP(0,NY,NX)
      this%histr_1D_SURF_ELEV_col(ncol)   = -CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX)+DLYR(3,0,NY,NX)
      this%histr_1D_SURF_tLITR_N_FLX_col(ncol)   = URSDN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_AMENDED_N_col(ncol)    = FertNFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tNH4X_col(ncol)        = UNH4(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tNO3_col(ncol)        = UNO3(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tMICRO_N_col(ncol)     = TONT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_TEMP_LITR_col(ncol)   = TCS(0,NY,NX)
      this%histr_1D_TEMP_SNOW_col(ncol)   = TCSnow(1,NY,NX)
      this%histr_1D_SURF_tLITR_C_FLX_col(ncol)   = URSDC(NY,NX)/AREA(3,NU(NY,NX),NY,NX)

      this%histr_1D_AMENDED_C_col(ncol)   = AmendCFlx_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_CO2_FLX_col(ncol)     = SurfGasFlx(idg_CO2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tMICRO_C_col(ncol)     = TOMT(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_OMC_LITR_col(ncol)    = ORGC(0,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ATM_CO2_col(ncol)     = CO2E(NY,NX)
      this%histr_1D_NBP_col(ncol)         = Eco_NBP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ECO_HVST_C_col(ncol)  = EcoHavstElmnt_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ECO_LAI_col(ncol)     = CanopyLeafArea_grd(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ECO_GPP_col(ncol)     = Eco_GPP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ECO_RA_col(ncol)      = Eco_AutoR_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ECO_NPP_col(ncol)     = Eco_NPP_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ECO_RH_col(ncol)      = Eco_HR_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tDIC_col(ncol)     = DIC_mass_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tSTANDING_DEAD_C_col(ncol)    = StandingDeadChemElmnt_col(ielmc,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tSTANDING_DEAD_N_col(ncol)    = StandingDeadChemElmnt_col(ielmn,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_tSTANDING_DEAD_P_col(ncol)    = StandingDeadChemElmnt_col(ielmp,NY,NX)/AREA(3,NU(NY,NX),NY,NX)            
      this%histr_1D_tPRECN_col(ncol)       = 1000.0_r8*URAIN(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ET_col(ncol)          = 1000.0_r8*UEVAP(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_N2O_LITR_col(ncol)    = trc_solcl_vr(idg_N2O,0,NY,NX)
      this%histr_1D_NH3_LITR_col(ncol)    = trc_solcl_vr(idg_NH3,0,NY,NX)
      this%histr_1D_SOL_RADN_col(ncol)    = RadSWSolarBeam_col(NY,NX)*277.8_r8
      this%histr_1D_AIR_TEMP_col(ncol)    = TCA(NY,NX)
      this%histr_1D_HUM_col(ncol)         = VPK(NY,NX)
      this%histr_1D_WIND_col(ncol)        = WindSpeedAtm(NY,NX)/3600.0_r8
      this%histr_1D_PREC_col(ncol)        = (RainFalPrec(NY,NX)+SnoFalPrec(NY,NX))*1000.0_r8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_SOIL_RN_col(ncol)     = HeatByRadiation(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_SOIL_LE_col(ncol)     = HeatEvapAir2Surf(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_SOIL_H_col(ncol)      = HeatSensAir2Surf(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_SOIL_G_col(ncol)      =-(HeatNet2Surf(NY,NX)-HeatSensVapAir2Surf(NY,NX))*277.8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ECO_RN_col(ncol)      = Eco_NetRad_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ECO_LE_col(ncol)      = Eco_Heat_Latent_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_Eco_Heat_col(ncol)       = Eco_Heat_Sens_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_ECO_G_col(ncol)       = Eco_Heat_Grnd_col(NY,NX)*277.8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_O2_LITR_col(ncol)     = trc_solcl_vr(idg_O2,0,NY,NX)
      this%histr_1D_SOIL_CO2_FLX_col(ncol)= SurfGasFlx(idg_CO2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815_r8
      this%histr_1D_ECO_CO2_FLX_col(ncol) = Eco_NEE_col(NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815_r8
      this%histr_1D_CH4_FLX_col(ncol)     = SurfGasFlx(idg_CH4,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.14815_r8
      this%histr_1D_O2_FLX_col(ncol)      = SurfGasFlx(idg_O2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*8.68056_r8
      this%histr_1D_CO2_LITR_col(ncol)    = trc_solcl_vr(idg_CO2,0,NY,NX)
      this%histr_1D_EVAPN_col(ncol)       = VapXAir2GSurf(NY,NX)*1000.0_r8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_RUNOFF_FLX_col(ncol)      = -WQRH(NY,NX)*1000.0_r8/TAREA 
      this%histr_1D_SEDIMENT_FLX_col(ncol)    = USEDOU(NY,NX)*1000.0_r8/TAREA
      this%histr_1D_tSWC_col(ncol)     = UVLWatMicP(NY,NX)*1000.0_r8/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_DISCHG_FLX_col(ncol)      = FWatDischarge(NY,NX)*1000.0_r8/TAREA
      this%histr_1D_SNOWPACK_col(ncol)    = AZMAX1((VcumDrySnoWE(NY,NX)+VcumIceSnow(NY,NX)*DENSICE+VcumWatSnow(NY,NX))*1000.0_r8/AREA(3,NU(NY,NX),NY,NX))
      this%histr_1D_SURF_WTR_col(ncol)    = THETWZ(0,NY,NX)
      this%histr_1D_SURF_ICE_col(ncol)    = THETIZ(0,NY,NX)
      this%histr_1D_ACTV_LYR_col(ncol)    = -(ActiveLayDepth(NY,NX)-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX))
      this%histr_1D_WTR_TBL_col(ncol)     = -(DepthInternalWTBL(NY,NX)-CumDepth2LayerBottom(NU(NY,NX)-1,NY,NX))
      this%histr_1D_sN2O_FLX_col(ncol)     =  SurfGasFlx(idg_N2O,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_sN2G_FLX_col(ncol)     =  SurfGasFlx(idg_N2,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      this%histr_1D_sNH3_FLX_col(ncol)     =  SurfGasFlx(idg_NH3,NY,NX)/AREA(3,NU(NY,NX),NY,NX)

      DO L=1,JZ
        this%histr_2D_tSOC_vr_col(ncol,L) =  ORGC(L,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_2D_CO2_vr_col(ncol,L)  =  trc_solcl_vr(idg_CO2,L,NY,NX)
        this%histr_2D_CH4_vr_col(ncol,L)  =  trc_solcl_vr(idg_CH4,L,NY,NX)
        this%histr_2D_O2_vr_col(ncol,L)   =  trc_solcl_vr(idg_O2,L,NY,NX)
        this%histr_2D_N2O_vr_col(ncol,L)  =  trc_solcl_vr(idg_N2O,L,NY,NX)
        this%histr_2D_NH3_vr_col(ncol,L)  =  trc_solcl_vr(idg_NH3,L,NY,NX)
        this%histr_2D_TEMP_vr_col(ncol,L) =  TCS(L,NY,NX)
        this%histr_2D_vWATER_vr_col(ncol,L)=  THETWZ(L,NY,NX)
        this%histr_2D_vICE_vr_col(ncol,L)  =  THETIZ(L,NY,NX)
        this%histr_2D_PSI_vr_col(ncol,L)  =  PSISoilMatricP(L,NY,NX)+PSISoilOsmotic(L,NY,NX)     
        this%histr_2D_cNH4t_vr_col(ncol,L)=  safe_adb(trc_solml_vr(ids_NH4,L,NY,NX)+trc_solml_vr(ids_NH4B,L,NY,NX) &
                                               +natomw*(trcx_solml(idx_NH4,L,NY,NX)+trcx_solml(idx_NH4B,L,NY,NX)),SoilMicPMassLayer(L,NY,NX))
        this%histr_2D_cNO3t_vr_col(ncol,L)= safe_adb(trc_solml_vr(ids_NO3,L,NY,NX)+trc_solml_vr(ids_NO3B,L,NY,NX) &
                                               +trc_solml_vr(ids_NO2,L,NY,NX)+trc_solml_vr(ids_NO2B,L,NY,NX),SoilMicPMassLayer(L,NY,NX))
        this%histr_2D_cPO4_vr_col(ncol,L) = safe_adb(trc_solml_vr(ids_H1PO4,L,NY,NX)+trc_solml_vr(ids_H1PO4B,L,NY,NX) &
                                               +trc_solml_vr(ids_H2PO4,L,NY,NX)+trc_solml_vr(ids_H2PO4B,L,NY,NX),VLWatMicP(L,NY,NX))
        this%histr_2D_cEXCH_P_vr_col(ncol,L)= patomw*safe_adb(trcx_solml(idx_HPO4,L,NY,NX)+trcx_solml(idx_H2PO4,L,NY,NX) &
                                               +trcx_solml(idx_HPO4B,L,NY,NX)+trcx_solml(idx_H2PO4B,L,NY,NX),SoilMicPMassLayer(L,NY,NX))
        this%histr_2D_ECND_vr_col(ncol,L)     = ECND(L,NY,NX)
      ENDDO

      DO NZ=1,NP0(NY,NX)
        nptc=get_pft(NZ,NY,NX)
        this%histr_1D_MIN_LWP_ptc(nptc)      = PSICanPDailyMin(NZ,NY,NX)
        this%histr_1D_LEAF_PC_ptc(nptc)      = safe_adb(LeafChemElmnts_pft(ielmp,NZ,NY,NX)+CanopyNonstructElements_pft(ielmp,NZ,NY,NX), &
                                                 LeafChemElmnts_pft(ielmc,NZ,NY,NX)+CanopyNonstructElements_pft(ielmc,NZ,NY,NX))
        this%histr_1D_CAN_RN_ptc(nptc)       = 277.8_r8*RadNet2CanP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_CAN_LE_ptc(nptc)       = 277.8_r8*EvapTransHeatP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_CAN_H_ptc(nptc)        = 277.8_r8*HeatXAir2PCan(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_CAN_G_ptc(nptc)        = 277.8_r8*HeatStorCanP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_CAN_TEMP_ptc(nptc)     = TCelciusCanopy_pft(NZ,NY,NX)
        this%histr_1D_TEMP_FN_ptc(nptc)      = fTgrowCanP(NZ,NY,NX)
        this%histr_1D_CAN_CO2_FLX_ptc(nptc)  = CO2NetFix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)*23.148_r8
        this%histr_1D_CAN_GPP_ptc(nptc)      = GrossCO2Fix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_CAN_RA_ptc(nptc)       = CanopyPlusNoduRespC_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_cTNC_ptc(nptc)         = CanopyNonstructElementConc_pft(ielmc,NZ,NY,NX)
        this%histr_1D_cTNN_ptc(nptc)         = CanopyNonstructElementConc_pft(ielmn,NZ,NY,NX)
        this%histr_1D_cTNP_ptc(nptc)         = CanopyNonstructElementConc_pft(ielmp,NZ,NY,NX)
        this%histr_1D_STOML_RSC_CO2_ptc(nptc)= CanPStomaResistH2O_pft(NZ,NY,NX)*1.56_r8*3600.0_r8
        this%histr_1D_BLYR_RSC_CO2_ptc(nptc) = CanopyBndlResist_pft(NZ,NY,NX)*1.34_r8*3600.0_r8
        this%histr_1D_CAN_CO2_ptc(nptc)      = CanopyGasCO2_pft(NZ,NY,NX)
        this%histr_1D_LAI_ptc(nptc)          = CanopyArea_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_PSI_CAN_ptc(nptc)      = PSICanopy_pft(NZ,NY,NX)
        this%histr_1D_TURG_CAN_ptc(nptc)     = PSICanopyTurg_pft(NZ,NY,NX)
        this%histr_1D_STOM_RSC_H2O_ptc(nptc) = CanPStomaResistH2O_pft(NZ,NY,NX)*3600.0_r8
        this%histr_1D_BLYR_RSC_H2O_ptc(nptc) = CanopyBndlResist_pft(NZ,NY,NX)*3600.0_r8
        this%histr_1D_TRANSPN_ptc(nptc)      = Transpiration_pft(NZ,NY,NX)*1000.0_r8/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_O2_STRESS_ptc(nptc)    = PlantO2Stress(NZ,NY,NX)
        this%histr_1D_NH4_UPTK_FLX_ptc(nptc)     = RootNH4Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_NO3_UPTK_FLX_ptc(nptc)     = RootNO3Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_N2_FIXN_FLX_ptc(nptc)      = RootN2Fix_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_cNH3_FLX_ptc(nptc)      = RNH3C(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_PO4_UPTK_FLX_ptc(nptc)     = RootH2PO4Uptake_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_SHOOT_C_ptc(nptc)      = ShootChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_Plant_C_ptc(nptc)      = (ShootChemElmnts_pft(ielmc,NZ,NY,NX)+RootChemElmnts_pft(ielmc,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_LEAF_C_ptc(nptc)       = LeafChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_Petiole_C_ptc(nptc)       = PetioleChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_STALK_C_ptc(nptc)      = StalkChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_RESERVE_C_ptc(nptc)    = ReserveChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_HUSK_C_ptc(nptc)       = (HuskChemElmnts_pft(ielmc,NZ,NY,NX)+EarChemElmnts_pft(ielmc,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_GRAIN_C_ptc(nptc)      = GrainChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_ROOT_C_ptc(nptc)       = RootChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_NODULE_C_ptc(nptc)        = NoduleChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_STORED_C_ptc(nptc)     = NonstructalChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_GRAIN_NO_ptc(nptc)     = CanopySeedNumber_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_LAIb_ptc(nptc)         = CanopyLeafArea_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_EXUD_C_FLX_ptc(nptc)       = PlantExudChemElmntCum_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_LITRf_C_FLX_ptc(nptc)      = LitrfallChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_SURF_LITRf_C_FLX_ptc(nptc) = SurfLitrfallChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_AUTO_RESP_FLX_ptc(nptc)    = GrossResp_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_ABV_GRD_RESP_FLX_ptc(nptc) = CanopyPlusNoduRespC_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_HVST_C_FLX_ptc(nptc)       = EcoHavstElmnt_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_HVST_N_FLX_ptc(nptc)       = EcoHavstElmnt_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)        
        this%histr_1D_HVST_P_FLX_ptc(nptc)       = EcoHavstElmnt_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)        
        this%histr_1D_PLANT_BALANCE_C_ptc(nptc) = ElmntBalanceCum_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_PLANT_BALANCE_N_ptc(nptc)    = ElmntBalanceCum_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_PLANT_BALANCE_P_ptc(nptc)    = ElmntBalanceCum_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_STANDING_DEAD_C_ptc(nptc)   = StandingDeadChemElmnts_pft(ielmc,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_FIREp_CO2_FLX_ptc(nptc)    = CO2ByFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_FIREp_CH4_FLX_ptc(nptc)    = CH4ByFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_NPP_ptc(nptc)          = NetPrimaryProductvity_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_CAN_HT_ptc(nptc)       = CanopyHeight_pft(NZ,NY,NX)
        this%histr_1D_POPN_ptc(nptc)         = PlantPopulation_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_tTRANSPN_ptc(nptc)     =-ETCanopy_pft(NZ,NY,NX)*1000.0/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_WTR_STRESS_ptc(nptc)   = HoursCanopyPSITooLow(NZ,NY,NX)
        this%histr_1D_OXY_STRESS_ptc(nptc)   = PlantO2Stress(NZ,NY,NX)
        this%histr_1D_SHOOT_N_ptc(nptc)      = ShootChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_Plant_N_ptc(nptc)      = (ShootChemElmnts_pft(ielmn,NZ,NY,NX)+RootChemElmnts_pft(ielmn,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_LEAF_N_ptc(nptc)       = LeafChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_Petiole_N_ptc(nptc)       = PetioleChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_STALK_N_ptc(nptc)      = StalkChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_RESERVE_N_ptc(nptc)    = ReserveChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_HUSK_N_ptc(nptc)       = (HuskChemElmnts_pft(ielmn,NZ,NY,NX)+EarChemElmnts_pft(ielmn,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_GRAIN_N_ptc(nptc)      = GrainChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_ROOT_N_ptc(nptc)       = RootChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_NODULE_N_ptc(nptc)        = NoduleChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_STORED_N_ptc(nptc)     = NonstructalChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_EXUD_N_FLX_ptc(nptc)       = PlantExudChemElmntCum_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_LITRf_N_FLX_ptc(nptc)      = LitrfallChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_TL_N_FIXED_FLX_ptc(nptc)   = PlantN2FixCum_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_NH3can_FLX_ptc(nptc)   = NH3EmiCum_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_STANDING_DEAD_N_ptc(nptc)   = StandingDeadChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_FIREp_N_FLX_ptc(nptc)       = NH3byFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_SURF_LITRf_N_FLX_ptc(nptc)  = SurfLitrfallChemElmnts_pft(ielmn,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_SHOOT_P_ptc(nptc)      = ShootChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_Plant_P_ptc(nptc)      = (ShootChemElmnts_pft(ielmp,NZ,NY,NX)+RootChemElmnts_pft(ielmp,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)        
        this%histr_1D_LEAF_P_ptc(nptc)       = LeafChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_Petiole_P_ptc(nptc)       = PetioleChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_STALK_P_ptc(nptc)      = StalkChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_RESERVE_P_ptc(nptc)    = ReserveChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_HUSK_P_ptc(nptc)       = (HuskChemElmnts_pft(ielmp,NZ,NY,NX)+EarChemElmnts_pft(ielmp,NZ,NY,NX))/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_GRAIN_P_ptc(nptc)      = GrainChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_ROOT_P_ptc(nptc)       = RootChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_NODULE_P_ptc(nptc)        = NoduleChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_STORED_P_ptc(nptc)     = NonstructalChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_EXUD_P_FLX_ptc(nptc)       = PlantExudChemElmntCum_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_LITRf_P_FLX_ptc(nptc)      = LitrfallChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_STANDING_DEAD_P_ptc(nptc)   = StandingDeadChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_FIREp_P_FLX_ptc(nptc)       = PO4byFire_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_SURF_LITRf_P_FLX_ptc(nptc) = SurfLitrfallChemElmnts_pft(ielmp,NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        this%histr_1D_BRANCH_NO_ptc(nptc)    = NumOfBranches_pft(NZ,NY,NX)
        this%histr_1D_LEAF_NC_ptc(nptc)      = safe_adb(LeafChemElmnts_pft(ielmn,NZ,NY,NX)+CanopyNonstructElements_pft(ielmn,NZ,NY,NX),&
                                                 LeafChemElmnts_pft(ielmc,NZ,NY,NX)+CanopyNonstructElements_pft(ielmc,NZ,NY,NX))
        if(NumOfMainBranch_pft(NZ,NY,NX)> 0)then
          DO KN=10,0,-1
            IF(KN==0)THEN
              this%histr_1D_Growth_Stage_ptc(nptc)   =KN
            ELSE
              if(iPlantCalendar_brch(KN,NumOfMainBranch_pft(NZ,NY,NX),NZ,NY,NX)>0)then
                this%histr_1D_Growth_Stage_ptc(nptc) =KN
                exit    
              endif  
            ENDIF  
          ENDDO
          DO NB=1,NumOfMainBranch_pft(NZ,NY,NX)
            this%histr_2D_LEAF_NODE_NO_ptc(nptc,NB) = NumOfLeaves_brch(NB,NZ,NY,NX)
            this%histr_2D_RUB_ACTVN_ptc(nptc,NB)  = RubiscoActivity_brpft(NB,NZ,NY,NX)
          ENDDO
        endif
        DO L=1,JZ
          this%histr_2D_PSI_RT_vr_ptc(nptc,L)  = PSIRoot_vr(ipltroot,L,NZ,NY,NX)
          this%histr_2D_prtUP_NH4_vr_ptc(nptc,L)  = (sum(RootNutUptake_pvr(ids_NH4,:,L,NZ,NY,NX))+sum(RootNutUptake_pvr(ids_NH4B,:,L,NZ,NY,NX)))/AREA(3,L,NY,NX)
          this%histr_2D_prtUP_NO3_vr_ptc(nptc,L)  = (sum(RootNutUptake_pvr(ids_NO3,:,L,NZ,NY,NX))+sum(RootNutUptake_pvr(ids_NO3B,:,L,NZ,NY,NX)))/AREA(3,L,NY,NX)
          this%histr_2D_prtUP_PO4_vr_ptc(nptc,L)  = (sum(RootNutUptake_pvr(ids_H2PO4,:,L,NZ,NY,NX))+sum(RootNutUptake_pvr(ids_H2PO4B,:,L,NZ,NY,NX)))/AREA(3,L,NY,NX)
          this%histr_2D_DNS_RT_vr_ptc(nptc,L)  = RootLenthDensPerPopu_pvr(ipltroot,L,NZ,NY,NX)*PlantPopulation_pft(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
        ENDDO
      ENDDO 
    ENDDO 
  ENDDO    
  end subroutine hist_update
  
end module HistDataType
