module AqueChemDatatype
  !!
  !Description:
  !Data type for aqueous chemistry and transport
  !Note: [d-2] represents per grid area
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMCtrlMod, only : salt_model
  use EcoSIMConfig, only : jcplx=> jcplxc
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  

  real(r8),target,allocatable ::  CAL_vr(:,:,:)                         !soil Al content, [mg Al kg-1]
  real(r8),target,allocatable ::  CFE_vr(:,:,:)                         !soil Fe content, [mg Fe kg-3]
  real(r8),target,allocatable ::  CCA_vr(:,:,:)                         !soil Ca content, [mg Ca kg-3]
  real(r8),target,allocatable ::  CMG_vr(:,:,:)                         !soil Mg content, [mg Mg kg-3]
  real(r8),target,allocatable ::  CNA_vr(:,:,:)                         !soil Na content, [mg Na kg-3]
  real(r8),target,allocatable ::  CKA_vr(:,:,:)                         !soil K content, [mg K kg-3]
  real(r8),target,allocatable ::  CSO4_vr(:,:,:)                        !soil SO4 content, [mg S kg-3]
  real(r8),target,allocatable ::  CCL_vr(:,:,:)                         !soil Cl content, [mg Cl kg-1]
  real(r8),target,allocatable ::  CALOH_vr(:,:,:)                       !soil AlOH3 content, [mg Al kg-1]
  real(r8),target,allocatable ::  CFEOH_vr(:,:,:)                       !soil FeOH3 content, [mg Fe kg-1]
  real(r8),target,allocatable ::  CCACO_vr(:,:,:)                       !soil CaCO3 content, [mg Ca kg-1]
  real(r8),target,allocatable ::  CCASO_vr(:,:,:)                       !soil CaSO4 content, [mg Ca kg-1]
  real(r8),target,allocatable ::  CALPO_vr(:,:,:)                       !soil AlPO4 content, [mg P kg-1]
  real(r8),target,allocatable ::  CFEPO_vr(:,:,:)                       !soil FePO4 content, [mg P kg-1]
  real(r8),target,allocatable ::  CCAPD_vr(:,:,:)                       !soil CaHPO4 content, [mg P kg-1]
  real(r8),target,allocatable ::  CCAPH_vr(:,:,:)                       !soil apatite content, [mg P kg-1]
  real(r8),target,allocatable ::  GKC4_vr(:,:,:)                        !Ca-NH4 Gapon selectivity coefficient, [-]
  real(r8),target,allocatable ::  GKCH_vr(:,:,:)                        !Ca-H Gapon selectivity coefficient, [-]
  real(r8),target,allocatable ::  GKCA_vr(:,:,:)                        !Ca-Al Gapon selectivity coefficient, [-]
  real(r8),target,allocatable ::  GKCM_vr(:,:,:)                        !Ca-Mg Gapon selectivity coefficient, [-]
  real(r8),target,allocatable ::  GKCN_vr(:,:,:)                        !Ca-Na Gapon selectivity coefficient, [-]
  real(r8),target,allocatable ::  GKCK_vr(:,:,:)                        !Ca-K Gapon selectivity coefficient, [-]

  real(r8),target,allocatable :: trcsalt_rain_conc(:,:,:)            !salt tracer concentration in rain [g m-3]
  real(r8),target,allocatable :: trcSalt_solml_vr(:,:,:,:)           !soil aqueous salt content micropre, [mol d-2]
  real(r8),target,allocatable :: trcx_solml_vr(:,:,:,:)              !exchangeable tracers, [mol d-2]
  real(r8),target,allocatable :: trcp_saltpml_vr(:,:,:,:)            !salt precipitate in micropore

  real(r8),target,allocatable ::  ElectricConductivity_vr(:,:,:)     !electrical conductivity , [dS m-1]
  real(r8),target,allocatable ::  SolutesIonStrenth_vr(:,:,:)        !solution ion strength, [mol m-3]
  real(r8),target,allocatable ::  SolutesIonConc_vr(:,:,:)           !solution ion concentratiom, [mol m-3]

  real(r8),target,allocatable :: trcSalt_soHml_vr(:,:,:,:)
  real(r8),target,allocatable :: trcSalt_XFHS_3D(:,:,:,:,:)
  real(r8),target,allocatable :: trcSalt3DFlo2Cell_3D(:,:,:,:,:)
  real(r8),target,allocatable :: trcSaltIonNumber(:)                 !number of ions when the salt is fully dissociated
  real(r8),target,allocatable ::  DOM_Mac2MicPore_flx_vr(:,:,:,:,:)       !total DOC micropore-macropore transfer, [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_Mac2MicPore_flx_vr(:,:,:,:)       !total non-salt solute micropore->macropore transfer, [g d-2 h-1]
  real(r8),target,allocatable ::  trcSalt_XFXS_vr(:,:,:,:)              !total salt micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  trcn_GeoChem_soil_vr(:,:,:,:)        !total solute NH4 transformation non-band, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_sol_NH3_soil_vr(:,:,:)          !total solute NH3 transformation non-band, [mol d-2 h-1]

  real(r8),target,allocatable ::  trcn_RChem_band_soil_vr(:,:,:,:)   !total solute nutrient transformation band, [mol d-2 h-1]
  real(r8),target,allocatable ::  trcSalt_TR_vr(:,:,:,:)                !total salt solute transformation non-band, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_HCO3_col(:,:,:)                 !total solute HCO3 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TProd_CO2_geochem_soil_vr(:,:,:)        !total solute CO2 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TRH2O(:,:,:)                       !total solute H2O transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_FeO3H3_soil_vr(:,:,:)              !total solute FeOH3 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_H_p_sorbed_soil_vr(:,:,:)                       !total adsorbed H transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_Al_sorbed_soil_vr(:,:,:)                       !total adsorbed Al transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_Ca_sorbed_soil_vr(:,:,:)                       !total adsorbed Ca transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_Mg_sorbed_soil_vr(:,:,:)                       !total adsorbed Mg transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_Na_sorbed_soil_vr(:,:,:)                       !total adsorbed Na transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_K_sorbed_soil_vr(:,:,:)                       !total adsorbed K transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_HCO3_sorbed_soil_vr(:,:,:)                       !total adsorbed COOH transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_AlO2H2_sorbed_soil_vr(:,:,:)                      !total adsorbed AlOH2 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_KSO4_soil_soil_vr(:,:,:)                       !total solute KSO4 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_Fe_sorbed_soil_vr(:,:,:)                       !total Fe adsorption
  real(r8),target,allocatable ::  TR_FeO2H2_sorbed_soil_vr(:,:,:)                      !total FeOH2 adsorption
  real(r8),target,allocatable ::  trcx_TRSoilChem_vr(:,:,:,:)                   !total adsorbed OH- transformation non-band, [mol d-2 h-1]
  real(r8),target,allocatable ::  Txchem_CO2_vr(:,:,:)                       !total solute CO2 transformation boundary, [mol d-2 h-1]
  real(r8),target,allocatable ::  TBION_vr(:,:,:)                       !total solute ion transformation boundary, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_gas_NH3_geochem_vr(:,:,:)                       !total gaseous NH3 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  trcp_RChem_soil(:,:,:,:)                   !total precipitated P containing transformation non-band, [mol d-2 h-1]
  real(r8),target,allocatable ::  trcg_AquaAdv_flx_snvr(:,:,:,:)
  real(r8),target,allocatable ::  trcn_AquaAdv_flx_snvr(:,:,:,:)
  real(r8),target,allocatable ::  trcSaltFlo2SnowLay(:,:,:,:)

  private :: InitAllocate
  contains

  subroutine InitAquaChem
  implicit none

  call InitAllocate

  end subroutine InitAquaChem

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(CAL_vr(JZ,JY,JX));      CAL_vr=0._r8
  allocate(CFE_vr(JZ,JY,JX));      CFE_vr=0._r8
  allocate(CCA_vr(JZ,JY,JX));      CCA_vr=0._r8
  allocate(CMG_vr(JZ,JY,JX));      CMG_vr=0._r8
  allocate(CNA_vr(JZ,JY,JX));      CNA_vr=0._r8
  allocate(CKA_vr(JZ,JY,JX));      CKA_vr=0._r8
  allocate(CSO4_vr(JZ,JY,JX));     CSO4_vr=0._r8
  allocate(CCL_vr(JZ,JY,JX));      CCL_vr=0._r8
  allocate(CALOH_vr(JZ,JY,JX));    CALOH_vr=0._r8
  allocate(CFEOH_vr(JZ,JY,JX));    CFEOH_vr=0._r8
  allocate(CCACO_vr(JZ,JY,JX));    CCACO_vr=0._r8
  allocate(CCASO_vr(JZ,JY,JX));    CCASO_vr=0._r8
  allocate(CALPO_vr(JZ,JY,JX));    CALPO_vr=0._r8
  allocate(CFEPO_vr(JZ,JY,JX));    CFEPO_vr=0._r8
  allocate(CCAPD_vr(JZ,JY,JX));    CCAPD_vr=0._r8
  allocate(CCAPH_vr(JZ,JY,JX));    CCAPH_vr=0._r8
  allocate(GKC4_vr(JZ,JY,JX));     GKC4_vr=0._r8
  allocate(GKCH_vr(JZ,JY,JX));     GKCH_vr=0._r8
  allocate(GKCA_vr(JZ,JY,JX));     GKCA_vr=0._r8
  allocate(GKCM_vr(JZ,JY,JX));     GKCM_vr=0._r8
  allocate(GKCN_vr(JZ,JY,JX));     GKCN_vr=0._r8
  allocate(GKCK_vr(JZ,JY,JX));     GKCK_vr=0._r8
  allocate(trcx_solml_vr(idx_beg:idx_end,0:JZ,JY,JX));trcx_solml_vr=0._r8
  allocate(trcp_saltpml_vr(idsp_beg:idsp_end,0:JZ,JY,JX)); trcp_saltpml_vr=0._r8
  allocate(ElectricConductivity_vr(JZ,JY,JX));     ElectricConductivity_vr=0._r8
  allocate(SolutesIonStrenth_vr(JZ,JY,JX));     SolutesIonStrenth_vr=0._r8
  allocate(SolutesIonConc_vr(JZ,JY,JX));     SolutesIonConc_vr=0._r8

  allocate(DOM_Mac2MicPore_flx_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Mac2MicPore_flx_vr=0._r8
  allocate(trcs_Mac2MicPore_flx_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_Mac2MicPore_flx_vr=0._r8

  allocate(trcn_GeoChem_soil_vr(ids_nut_beg:ids_nuts_end,0:JZ,JY,JX)); trcn_GeoChem_soil_vr=0._r8

  allocate(TR_sol_NH3_soil_vr(0:JZ,JY,JX));  TR_sol_NH3_soil_vr=0._r8

  allocate(trcn_RChem_band_soil_vr(ids_nutb_beg:ids_nutb_end,JZ,JY,JX)); trcn_RChem_band_soil_vr=0._r8

  allocate(trcx_TRSoilChem_vr(idx_beg:idx_end,0:JZ,JY,JX));  trcx_TRSoilChem_vr=0._r8

  allocate(TR_HCO3_col(JZ,JY,JX));    TR_HCO3_col=0._r8
  allocate(TProd_CO2_geochem_soil_vr(JZ,JY,JX));    TProd_CO2_geochem_soil_vr=0._r8
  allocate(TRH2O(0:JZ,JY,JX));  TRH2O=0._r8
  allocate(TR_FeO3H3_soil_vr(JZ,JY,JX));    TR_FeO3H3_soil_vr=0._r8
  allocate(TR_H_p_sorbed_soil_vr(JZ,JY,JX));    TR_H_p_sorbed_soil_vr=0._r8
  allocate(TR_Al_sorbed_soil_vr(JZ,JY,JX));    TR_Al_sorbed_soil_vr=0._r8
  allocate(TR_Ca_sorbed_soil_vr(JZ,JY,JX));    TR_Ca_sorbed_soil_vr=0._r8
  allocate(TR_Mg_sorbed_soil_vr(JZ,JY,JX));    TR_Mg_sorbed_soil_vr=0._r8
  allocate(TR_Na_sorbed_soil_vr(JZ,JY,JX));    TR_Na_sorbed_soil_vr=0._r8
  allocate(TR_K_sorbed_soil_vr(JZ,JY,JX));    TR_K_sorbed_soil_vr=0._r8
  allocate(TR_HCO3_sorbed_soil_vr(JZ,JY,JX));    TR_HCO3_sorbed_soil_vr=0._r8
  allocate(TR_AlO2H2_sorbed_soil_vr(JZ,JY,JX));   TR_AlO2H2_sorbed_soil_vr=0._r8
  allocate(TR_KSO4_soil_soil_vr(JZ,JY,JX));    TR_KSO4_soil_soil_vr=0._r8
  allocate(TR_Fe_sorbed_soil_vr(JZ,JY,JX));    TR_Fe_sorbed_soil_vr=0._r8
  allocate(TR_FeO2H2_sorbed_soil_vr(JZ,JY,JX));   TR_FeO2H2_sorbed_soil_vr=0._r8
  allocate(Txchem_CO2_vr(JZ,JY,JX));    Txchem_CO2_vr=0._r8
  allocate(TBION_vr(0:JZ,JY,JX));  TBION_vr=0._r8
  allocate(TR_gas_NH3_geochem_vr(0:JZ,JY,JX));  TR_gas_NH3_geochem_vr=0._r8
  allocate(trcp_RChem_soil(idsp_beg:idsp_end,0:JZ,JY,JX)); trcp_RChem_soil=0._r8
  allocate(trcg_AquaAdv_flx_snvr(idg_beg:idg_NH3,JS,JY,JX)); trcg_AquaAdv_flx_snvr=0._r8
  allocate(trcn_AquaAdv_flx_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX)); trcn_AquaAdv_flx_snvr=0._r8
  if(salt_model)then
    allocate(trcSaltFlo2SnowLay(idsalt_beg:idsalt_end,JS,JY,JX)); trcSaltFlo2SnowLay=0._r8
    allocate(trcSalt3DFlo2Cell_3D(idsalt_beg:idsaltb_end,3,0:JD,JV,JH));trcSalt3DFlo2Cell_3D=0._r8
    allocate(trcSalt_XFHS_3D(idsalt_beg:idsaltb_end,3,JD,JV,JH));trcSalt_XFHS_3D=0._r8
    allocate(trcSalt_solml_vr(idsalt_beg:idsaltb_end,0:JZ,JY,JX));trcSalt_solml_vr=0._r8
    allocate(trcsalt_rain_conc(idsalt_beg:idsalt_end,JY,JX));trcsalt_rain_conc=0._r8
    allocate(trcSaltIonNumber(idsalt_beg:idsaltb_end))
    allocate(trcSalt_soHml_vr(idsalt_beg:idsaltb_end,JZ,JY,JX)); trcSalt_soHml_vr=0._r8
    allocate(trcSalt_XFXS_vr(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_XFXS_vr=0._r8
    allocate(trcSalt_TR_vr(idsalt_beg:idsaltb_end,JZ,JY,JX));    trcSalt_TR_vr=0._r8
  endif

  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructAquaChem
  use abortutils, only : destroy
  implicit none

  call destroy(CAL_vr)
  call destroy(CFE_vr)
  call destroy(CCA_vr)
  call destroy(CMG_vr)
  call destroy(CNA_vr)
  call destroy(CKA_vr)
  call destroy(CSO4_vr)
  call destroy(CCL_vr)
  call destroy(CALOH_vr)
  call destroy(CFEOH_vr)
  call destroy(CCACO_vr)
  call destroy(CCASO_vr)
  call destroy(CALPO_vr)
  call destroy(CFEPO_vr)
  call destroy(CCAPD_vr)
  call destroy(CCAPH_vr)
  call destroy(GKC4_vr)
  call destroy(GKCH_vr)
  call destroy(GKCA_vr)
  call destroy(GKCM_vr)
  call destroy(GKCN_vr)
  call destroy(GKCK_vr)

  call destroy(trcSalt_TR_vr)
  call destroy(trcx_solml_vr)
  call destroy(trcSalt_solml_vr)
  call destroy(trcsalt_rain_conc)
  call destroy(trcSalt_soHml_vr)
  call destroy(trcp_saltpml_vr)
  call destroy(ElectricConductivity_vr)
  call destroy(SolutesIonStrenth_vr)
  call destroy(SolutesIonConc_vr)
  call destroy(trcSalt_XFXS_vr)
  call destroy(trcSalt3DFlo2Cell_3D)
  call destroy(trcSalt_XFHS_3D)
  call destroy(TR_sol_NH3_soil_vr)
  call destroy(trcn_RChem_band_soil_vr)
  call destroy(TR_HCO3_col)
  call destroy(TProd_CO2_geochem_soil_vr)
  call destroy(TRH2O)
  call destroy(TR_FeO3H3_soil_vr)
  call destroy(TR_Al_sorbed_soil_vr)
  call destroy(TR_Ca_sorbed_soil_vr)
  call destroy(TR_Mg_sorbed_soil_vr)
  call destroy(TR_Na_sorbed_soil_vr)
  call destroy(TR_K_sorbed_soil_vr)
  call destroy(TR_HCO3_sorbed_soil_vr)
  call destroy(TR_AlO2H2_sorbed_soil_vr)
  call destroy(TR_KSO4_soil_soil_vr)
  call destroy(TR_Fe_sorbed_soil_vr)
  call destroy(TR_FeO2H2_sorbed_soil_vr)
  call destroy(Txchem_CO2_vr)
  call destroy(TBION_vr)
  call destroy(TR_gas_NH3_geochem_vr)
  call destroy(trcSaltIonNumber)
  call destroy(trcs_Mac2MicPore_flx_vr)
  call destroy(DOM_Mac2MicPore_flx_vr)
  call destroy(trcn_GeoChem_soil_vr)
  end subroutine DestructAquaChem

end module AqueChemDatatype
