module ChemIDMod

implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  integer :: cid_CO2S     !aqueous CO2  micropore	[g d-2]
  integer :: cid_H1PO4_2e_conc    !soil aqueous HPO4 content micropore non-band, [mol m-3]
  integer :: cid_H2PO4_1e_conc    !soil aqueous H2PO4 content micropore non-band, [mol m-3]
  integer :: cid_NH3_aqu_conc     !soil NH3 concentration in non-band soil, [mol m-3]
  integer :: cid_NH4_1p_conc     !soil NH4 concentration in non-band soil, [mol m-3]
  integer :: cid_ZNO3S    !NO3 mass non-band micropore, [g d-2]
  integer :: cid_ZHY      !soil aqueous H content micropore, [mol d-2]
  integer :: cid_ZOH      !soil aqueous OH content micropore, [mol d-2]
  integer :: cid_ZAL      !soil aqueous Al content micropore, [mol d-2]
  integer :: cid_ZFE      !soil aqueous Fe content micropore, [mol d-2]
  integer :: cid_ZCA      !soil aqueous Ca content micropore, [mol d-2]
  integer :: cid_ZMG      !soil aqueous Mg content micropore, [mol d-2]
  integer :: cid_ZNA      !soil aqueous Na content micropore, [mol d-2]
  integer :: cid_ZKA      !soil aqueous K content micropore, [mol d-2]
  integer :: cid_ZSO4     !soil aqueous SO4 content micropore, [mol d-2]
  integer :: cid_ZCL      !soil aqueous Cl content micropore, [mol d-2]
  integer :: cid_ZCO3     !soil aqueous CO3 content micropore, [mol d-2]
  integer :: cid_ZHCO3    !soil aqueous HCO3 content micropore, [mol d-2]
  integer :: cid_ZALOH1   !soil aqueous AlOH content micropore, [mol d-2]
  integer :: cid_ZALOH2   !soil aqueous AlOH2 content micropore, [mol d-2]
  integer :: cid_ZALOH3   !soil aqueous AlOH3 content micropore, [mol d-2]
  integer :: cid_ZALOH4   !soil aqueous AlOH4 content micropore, [mol d-2]
  integer :: cid_ZALS     !soil aqueous AlSO4 content micropore, [mol d-2]
  integer :: cid_ZFEOH1   !soil aqueous FeOH content micropore, [mol d-2]
  integer :: cid_ZFEOH2   !soil aqueous FeOH2 content micropore, [mol d-2]
  integer :: cid_ZFEOH3   !soil aqueous FeOH3 content micropore, [mol d-2]
  integer :: cid_ZFEOH4   !soil aqueous FeOH4 content micropore, [mol d-2]
  integer :: cid_ZFES     !soil aqueous FeSO4 content micropore, [mol d-2]
  integer :: cid_ZCAO     !soil aqueous CaOH2 content micropore, [mol d-2]
  integer :: cid_ZCAC     !soil aqueous CACO3 content micropore, [mol d-2]
  integer :: cid_ZCAH     !soil aqueous CaHCO3 content micropore, [mol d-2]
  integer :: cid_ZCAS     !soil aqueous CaSO4 content micropore, [mol d-2]
  integer :: cid_ZMGO     !soil aqueous MgOH content micropore, [mol d-2]
  integer :: cid_ZMGC     !soil aqueous MgCO3 content micropore, [mol d-2]
  integer :: cid_ZMGH     !soil aqueous MgHCO3 content micropore, [mol d-2]
  integer :: cid_ZMGS     !soil aqueous MgSO4 content micropore, [mol d-2]
  integer :: cid_ZNAC     !soil aqueous NaCO3 content micropore, [mol d-2]
  integer :: cid_ZNAS     !soil aqueous NaSO4 content micropore, [mol d-2]
  integer :: cid_ZKAS     !soil aqueous KSO4 content micropore, [mol d-2]
  integer :: cid_H0PO4    !soil aqueous PO4 content micropore non-band, [mol d-2]
  integer :: cid_H3PO4    !soil aqueous H3PO4 content micropore non-band, [mol d-2]
  integer :: cid_ZFE1P    !soil aqueous FeHPO4 content micropore non-band, [mol d-2]
  integer :: cid_ZFE2P    !soil aqueous FeH2PO4 content micropore non-band, [mol d-2]
  integer :: cid_ZCA0P    !soil aqueous CaPO4 content micropore non-band, [mol d-2]
  integer :: cid_ZCA1P    !soil aqueous CaHPO4 content micropore non-band, [mol d-2]
  integer :: cid_ZCA2P    !soil aqueous CaH4P2O8 content micropore non-band, [mol d-2]
  integer :: cid_ZMG1P    !soil aqueous MgHPO4 content micropore non-band, [mol d-2]
  integer :: cid_Precp_AlPO4_conc   !precipitated AlPO4 non-band, [mol m-3]
  integer :: cid_Precp_CaHPO4_conc   !precipitated CaHPO4 non-band soil, [mol m-3]
  integer :: cid_Precp_Ca5P3O12O3H3_conc   !precipitated Ca5(PO4)3OH hydroxyapatite non-band soil, [mol m-3]
  integer :: cid_Precp_CaH4P2O8_conc   !precipitated Ca(H2PO4)2 non-band soil, [mol m-3]
  integer :: cid_Precp_FePO4_conc   !precipitated FePO4 non-band soil, [mol m-3]
  integer :: cid_PALOH    !precipitated Al(OH)3, [mol d-2]
  integer :: cid_PFEOH    !precipitated Fe(OH)3, [mol d-2]
  integer :: cid_PCACO    !precipitated CaCO3, [mol d-2]
  integer :: cid_PCASO    !precipitated CaSO4, [mol d-2]
  integer :: cid_XHY      !exchangeable H(+) , [mol d-2]
  integer :: cid_XAL      !exchangeable Al, [mol d-2]
  integer :: cid_XFE      !exchangeable Fe, [mol d-2]
  integer :: cid_XCA      !exchangeable Ca, [mol d-2]
  integer :: cid_XMG      !exchangeable Mg, [mol d-2]
  integer :: cid_XNA      !exchangeable Na, [mol d-2]
  integer :: cid_XKA      !exchangeable K, [mol d-2]
  integer :: cid_XHC      !exchangeable COOH , [mol d-2]
  integer :: cid_XROH1_conc    !exchangeable OH  non-band, [mol d-2]
  integer :: cid_XALO2    !exchangeable AlOH2 , [mol d-2]
  integer :: cid_XFEO2    !exchangeable Fe(OH)2, [mol d-2]
  integer :: cid_XNH4_conc     !exchangeable NH4 non-band soil, [mol d-2]
  integer :: cid_XOH_conc    !exchangeable OH- non-band, [mol d-2]
  integer :: cid_XHPO4_conc    !exchangeable HPO4  non-band, [mol m-3]
  integer :: cid_XROH2_conc    !exchangeable OH2  non-band soil, [mol m-3]
  integer :: cid_XH2PO4_conc    !exchangeable H2PO4  non-band soil, [mol m-3]

  integer :: cid_NH4_1p_band_conc     !soil NH4 concentration in band soil, [mol m-3]
  integer :: cid_NH3_aqu_band_conc     !soil NH3 concentration in band soil, [mol m-3]
  integer :: cid_ZNO3B    !NO3 mass band micropore, [g d-2]
  integer :: cid_H1PO4_2e_band_conc    !soil aqueous HPO4 content micropore band, [mol m-3]
  integer :: cid_H2PO4_1e_band_conc    !soil aqueous H2PO4 content micropore  band, [mol m-3]
  integer :: cid_H0POB    !soil aqueous PO4 content micropore band, [mol d-2]
  integer :: cid_H3POB    !soil aqueous H3PO4 content micropore band, [mol d-2]
  integer :: cid_ZFE1PB   !soil aqueous FeHPO4 content micropore band, [mol d-2]
  integer :: cid_ZFE2PB   !soil aqueous FeH2PO4 content micropore band, [mol d-2
  integer :: cid_ZCA0PB   !soil aqueous CaPO4 content micropore band, [mol d-2]
  integer :: cid_ZCA1PB   !soil aqueous CaHPO4 content micropore band, [mol d-2]
  integer :: cid_ZCA2PB   !soil aqueous CaH4P2O8 content micropore band, [mol d-2]
  integer :: cid_ZMG1PB   !soil aqueous MgHPO4 content micropore band, [mol d-2]
  integer :: cid_PrecpB_AlPO4_conc   !precipitated AlPO4 band soil, [mol m-3]
  integer :: cid_PrecpB_CaHPO4_conc   !precipitated CaHPO4 band soil, [mol m-3]
  integer :: cid_PrecpB_Ca5P3O12O3H3_conc   !precipitated Ca5(PO4)3OH hydroxyapatite band soil, [mol m-3]
  integer :: cid_PrecpB_CaH4P2O8_conc   !precipitated CaH4P2O8 band soil, [mol m-3]
  integer :: cid_PrecpB_FePO4_con   !precipitated FePO4 band soil, [mol m-3]
  integer :: cid_XNH4_band_conc     !exchangeable NH4 band soil, [mol d-2]
  integer :: cid_XROH1_band_conc    !exchangeable OH- band, [mol d-2]
  integer :: cid_XHPO4_band_conc    !exchangeable HPO4 concentration band-soil, [mol m-3]
  integer :: cid_XH2PO4_band_conc    !exchangeable H2PO4 concentration band-soil, [mol m-3]
  integer :: cid_XROH_band_conc    !exchangeable OH band-soil, [mol m-3]
  integer :: cid_XROH2_band_conc    !exchangeable OH2 band-soil, [mol m-3]

  integer :: fid_TR_NH4_soil    !total solute NH4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_NH3_soil_vr    !total solute NH3 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_H1PO4_soil    !total solute HPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_H2PO4_soil    !total solute H2PO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_NH4_sorbed_soil    !total adsorbed NH4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_ROH_sorbed_soil    !total adsorbed OH transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_ROH2_sorbed_soil    !total adsorbed OH2 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_RHPO4_sorbed_soil    !total adsorbed HPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_RH2PO4_sorbed_soil    !total adsorbed H2PO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_ROH_sorbed_band_soil    !total adsorbed OH transformation band, [mol d-2 h-1]
  integer :: fid_TR_ROH2_sorbed_band_soil    !total adsorbed OH2 transformation band, [mol d-2 h-1]
  integer :: fid_TR_RHPO4_sorbed_band_soil    !total adsorbed HPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_RH2PO4_sorbed_band_soil    !total adsorbed H2PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_AlPO4_precip_soil   !total precipitated AlPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_FePO4_precip_soil   !total precipitated FePO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_CaHPO4_precip_soil   !total precipitated CaHPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_apatite_precip_soil   !total precipitated CaH4P2O8 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_CaH4P2O8_precip_soil   !total precipitated apatite transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_Al_3p_soil     !total solute Al transformation, [mol d-2 h-1]
  integer :: fid_TR_Fe_3p_soil     !total solute Fe transformation, [mol d-2 h-1]
  integer :: fid_TR_H_p_soil     !total solute H transformation, [mol d-2 h-1]
  integer :: fid_TR_Ca_2p_soil     !total solute Ca transformation, [mol d-2 h-1]
  integer :: fid_TR_Mg_2p_soil     !total solute Mg transformation, [mol d-2 h-1]
  integer :: fid_TR_Na_p_soil     !total solute Na transformation, [mol d-2 h-1]
  integer :: fid_TR_K_1p_soil     !total solute K transformation, [mol d-2 h-1]
  integer :: fid_TR_OH_1e_soil     !total solute OH transformation, [mol d-2 h-1]
  integer :: fid_TR_SO4_2e_soil    !total solute SO4 transformation, [mol d-2 h-1]
  integer :: fid_TR_CO3_2e_soil    !total solute CO3 transformation, [mol d-2 h-1]
  integer :: fid_TR_HCO3    !total solute HCO3 transformation, [mol d-2 h-1]
  integer :: fid_TR_CO2_gchem_soil_vr    !total solute CO2 transformation, [mol d-2 h-1]
  integer :: fid_TR_AlOH_soil    !total solute AlOH transformation, [mol d-2 h-1]
  integer :: fid_TR_AlO2H2_soil    !total solute AlOH2 transformation, [mol d-2 h-1]
  integer :: fid_TR_AlO3H3_soil    !total solute AlOH3 transformation, [mol d-2 h-1]
  integer :: fid_TR_AlO4H4_soil    !total solute AlOH4 transformation, [mol d-2 h-1]
  integer :: fid_TR_AlSO4_soil    !total solute AlSO4 transformation, [mol d-2 h-1]
  integer :: fid_TR_FeOH_soil    !total solute FeOH transformation, [mol d-2 h-1]
  integer :: fid_TR_FeO2H2_soil    !total solute FeOH2 transformation, [mol d-2 h-1]
  integer :: fid_TR_FeO3H3_soil    !total solute FeOH3 transformation, [mol d-2 h-1]
  integer :: fid_TR_FeO4H4_soil    !total solute FeOH4 transformation, [mol d-2 h-1]
  integer :: fid_TR_FeSO4_soil    !total solute FeSO4 transformation, [mol d-2 h-1]
  integer :: fid_TR_CaOH_soil    !total solute CaOH transformation, [mol d-2 h-1]
  integer :: fid_TR_CaCO3_soil    !total solute CaCO3 transformation, [mol d-2 h-1]
  integer :: fid_TR_CaHCO3_soil    !total solute CaHCO3 transformation, [mol d-2 h-1]
  integer :: fid_TR_CaSO4_soil    !total solute CaSO4 transformation, [mol d-2 h-1]
  integer :: fid_TR_MgOH_soil    !total solute MgOH transformation, [mol d-2 h-1]
  integer :: fid_TR_MgCO3_soil    !total solute MgCO3 transformation, [mol d-2 h-1]
  integer :: fid_TR_MgHCO3_soil    !total solute MgHCO3(+) transformation, [mol d-2 h-1]
  integer :: fid_TR_MgSO4_soil    !total solute MgSO4 transformation, [mol d-2 h-1]
  integer :: fid_TR_NaCO3_soil    !total solute NaCO3(-) transformation, [mol d-2 h-1]
  integer :: fid_TR_NaSO4_soil    !total solute NaSO4(-) transformation, [mol d-2 h-1]
  integer :: fid_TR_KSO4_soil    !total solute KSO4(-) transformation, [mol d-2 h-1]
  integer :: fid_TR_PO4_soil    !total solute PO4(---) transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_H3PO4_sorbed_soil    !total solute H3PO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_FeHPO4_soil    !total solute FeHPO4(+) transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_FeH2PO4_soil    !total solute FeH2PO4(++) transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_CaPO4_soil    !total solute CaPO4(-) transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_CaHPO4_soil    !total solute CaHPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_CaH4P2O8_soil    !total solute CaH4P2O8 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_MgHPO4_soil    !total solute MgHPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_H_p_sorbed_soil    !total adsorbed H transformation, [mol d-2 h-1]
  integer :: fid_TR_Al_sorbed_soil    !total adsorbed Al transformation, [mol d-2 h-1]
  integer :: fid_TR_Fe_sorbed_soil    !total Fe adsorption
  integer :: fid_TR_Ca_sorbed_soil    !total adsorbed Ca transformation, [mol d-2 h-1]
  integer :: fid_TR_Mg_sorbed_soil    !total adsorbed Mg transformation, [mol d-2 h-1]
  integer :: fid_TR_Na_sorbed_soil    !total adsorbed Na transformation, [mol d-2 h-1]
  integer :: fid_TR_K_sorbed_soil    !total adsorbed K transformation, [mol d-2 h-1]
  integer :: fid_TR_HCO3_sorbed_soil    !total adsorbed COOH transformation, [mol d-2 h-1]
  integer :: fid_TR_AlO2H2_sorbed_soil   !total adsorbed AlOH2 transformation, [mol d-2 h-1]
  integer :: fid_TR_FeO2H2_sorbed_soil   !total FeOH2 adsorption, [mol d-2 h-1]
  integer :: fid_TR_RO_sorbed_soil    !total adsorbed OH- transformation non-band, [mol d-2 h-1]
  integer :: fid_TR_RO_sorbed_band_soil    !total adsorbed OH- transformation band, [mol d-2 h-1]
  integer :: fid_TR_AlOH3_precip_soil   !total precipitated AlOH3 transformation, [mol d-2 h-1]
  integer :: fid_TR_FeOH3_precip_soil   !total precipitated FeOH3 transformation, [mol d-2 h-1]
  integer :: fid_TR_CaCO3_precip_soil   !total precipitated CaCO3 transformation, [mol d-2 h-1]
  integer :: fid_TR_CaSO4_precip_soil   !total precipitated CaSO4 transformation, [mol d-2 h-1]
  integer :: fid_TRH2O    !total solute H2O transformation, [mol d-2 h-1]
  integer :: fid_TBION    !total solute ion transformation, [mol d-2 h-1]
  integer :: fid_Txchem_CO2    !CO2 net change from all solute equilibria

  integer :: fid_TR_NH4_band_soil    !total solute NH4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_NH3_band_soil    !total solute NH3 transformation band, [mol d-2 h-1]
  integer :: fid_TR_H1PO4_band_soil    !total solute HPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_H2PO4_band_soil    !total solute H2PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_PO4_band_soil    !total solute PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_H3PO4_band_soil    !total solute H3PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_FeHPO4_band_soil    !total solute FeHPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_FeH2PO4_band_soil    !total solute FeH2PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_CaPO4_band_soil    !total solute CaPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_CaHPO4_band_soil    !total solute CaHPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_CaH4P2O8_band_soil    !total solute CaH4P2O8 transformation band, [mol d-2 h-1]
  integer :: fid_TR_MgHPO4_band_soil    !total solute MgHPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_NH4_sorbed_band_soil    !total adsorbed NH4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_AlPO4_precip_band_soil   !total precipitated AlPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_FePO4_precip_band_soil   !total precipitated FePO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_CaHPO4_precip_band_soil   !total precipitated CaHPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TR_apatite_precip_band_soil   !total precipitated CaH4P2O8 transformation band, [mol d-2 h-1]
  integer :: fid_TR_CaH4P2O8_precip_band_soil   !total precipitated apatite transformation band, [mol d-2 h-1]
  public :: getvarlist_nosalt
  contains

!------------------------------------------------------------------
  subroutine getvarlist_nosalt(nvars, varl, varlnml, unitl, vartypes)

  use bhistMod, only : hist_var_str_len,hist_unit_str_len,hist_var_lon_str_len
  use fileUtil, only :  var_flux_type, var_state_type
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(inout) :: varl(:)     !variable name
  character(len=hist_var_lon_str_len), intent(inout) :: varlnml(:)     !variable name
  character(len=hist_unit_str_len),intent(inout) :: unitl(:)
  integer                         ,intent(inout) :: vartypes(:)

  varl(cid_H1PO4_2e_conc)='H1PO4_2e_conc';varlnml(cid_H1PO4_2e_conc)='non-band soil aqueous HPO4 content micropore'
  unitl(cid_H1PO4_2e_conc)='mol m-3';vartypes(cid_H1PO4_2e_conc)=var_state_type

  varl(cid_H1PO4_2e_band_conc)='H1PO4_2e_band_conc';varlnml(cid_H1PO4_2e_band_conc)='band soil aqueous HPO4 content micropore'
  unitl(cid_H1PO4_2e_band_conc)='mol m-3';vartypes(cid_H1PO4_2e_band_conc)=var_state_type

  varl(cid_H2PO4_1e_conc)='H2PO4_1e_conc';varlnml(cid_H2PO4_1e_conc)='non-band soil micropore aqueous H2PO4 content'
  unitl(cid_H2PO4_1e_conc)='mol m-3';vartypes(cid_H2PO4_1e_conc)=var_state_type

  varl(cid_H2PO4_1e_band_conc)='H2PO4_1e_band_conc';varlnml(cid_H2PO4_1e_band_conc)='band soil micropore aqueous H2PO4 content'
  unitl(cid_H2PO4_1e_band_conc)='mol m-3';vartypes(cid_H2PO4_1e_band_conc)=var_state_type

  varl(cid_NH3_aqu_conc) ='NH3_aqu_conc';varlnml(cid_NH3_aqu_conc)='non-band soil NH3 concentration'
  unitl(cid_NH3_aqu_conc)='mol m-3';vartypes(cid_NH3_aqu_conc)=var_state_type

  varl(cid_NH3_aqu_band_conc) ='NH3_aqu_band_conc';varlnml(cid_NH3_aqu_band_conc)='band soil NH3 concentration'
  unitl(cid_NH3_aqu_band_conc)='mol m-3';vartypes(cid_NH3_aqu_band_conc)=var_state_type

  varl(cid_NH4_1p_conc)  ='NH4_1p_conc';varlnml(cid_NH4_1p_conc)='non-band soil NH4 concentration'
  unitl(cid_NH4_1p_conc)='mol m-3';vartypes(cid_NH4_1p_conc)=var_state_type

  varl(cid_NH4_1p_band_conc)  ='NH4_1p_band_conc';varlnml(cid_NH4_1p_band_conc)='band soil NH4 concentration'
  unitl(cid_NH4_1p_band_conc) ='mol m-3';vartypes(cid_NH4_1p_band_conc)=var_state_type

  varl(cid_XNH4_conc)  ='XNH4_conc';varlnml(cid_XNH4_conc)='non-band soil adsorbed NH4 concentration'
  unitl(cid_XNH4_conc) ='mol m-3';vartypes(cid_XNH4_conc)=var_state_type

  varl(cid_XNH4_band_conc)  ='XNH4_band_conc';varlnml(cid_XNH4_band_conc)='band soil adsorbed NH4 concentration'
  unitl(cid_XNH4_band_conc) ='mol m-3';vartypes(cid_XNH4_band_conc)=var_state_type

  varl(cid_XHPO4_band_conc) ='XHPO4_band_conc';varlnml(cid_XHPO4_band_conc)='band soil exchangeable HPO4 concentration'
  unitl(cid_XHPO4_band_conc)='mol m-3';vartypes(cid_XHPO4_band_conc)=var_state_type

  varl(cid_XH2PO4_band_conc) ='XH2PO4_band_conc';varlnml(cid_XH2PO4_band_conc)='band soil exchangeable H2PO4 concentration'
  unitl(cid_XH2PO4_band_conc)='mol m-3';vartypes(cid_XH2PO4_band_conc)=var_state_type

  varl(cid_XROH_band_conc) ='XROH_band_conc';varlnml(cid_XROH_band_conc)='band soil exchangeable site R-OH'
  unitl(cid_XROH_band_conc)='mol m-3';vartypes(cid_XROH_band_conc)=var_state_type

  varl(cid_XHPO4_conc) ='XHPO4_conc';varlnml(cid_XHPO4_conc)='non-band soil concentration of adsorbed HPO4'
  unitl(cid_XHPO4_conc)='mol m-3';vartypes(cid_XHPO4_conc)=var_state_type

  varl(cid_XROH2_band_conc) ='XROH2_band_conc';varlnml(cid_XROH2_band_conc)='band soil exchangeable site R-OH2'
  unitl(cid_XROH2_band_conc)='mol m-3';vartypes(cid_XROH2_band_conc)=var_state_type

  varl(cid_XH2PO4_conc) ='XH2PO4_conc';varlnml(cid_XH2PO4_conc)='non-band soil concentration of adsorbed HPO4'
  unitl(cid_XH2PO4_conc)='mol m-3';vartypes(cid_XH2PO4_conc)=var_state_type

  varl(cid_XROH1_conc) ='XROH1_conc';varlnml(cid_XROH1_conc)='non-band soil concentration of adsorption sites R-OH'
  unitl(cid_XROH1_conc)='mol m-3';vartypes(cid_XROH1_conc)=var_state_type

  varl(cid_XROH2_conc) ='XROH2_conc';varlnml(cid_XROH2_conc)='non-band soil exchangeable site R-OH2'
  unitl(cid_XROH2_conc)='mol m-3';vartypes(cid_XROH2_conc)=var_state_type

  varl(cid_Precp_AlPO4_conc)='Precp_AlPO4_conc';varlnml(cid_Precp_AlPO4_conc)='non-band soil precipitated AlPO4'
  unitl(cid_Precp_AlPO4_conc)='mol m-3';vartypes(cid_Precp_AlPO4_conc)=var_state_type

  varl(cid_PrecpB_AlPO4_conc) ='PrecpB_AlPO4_conc';varlnml(cid_PrecpB_AlPO4_conc)='band soil precipitated AlPO4'
  unitl(cid_PrecpB_AlPO4_conc)='mol m-3';vartypes(cid_PrecpB_AlPO4_conc)=var_state_type

  varl(cid_Precp_CaHPO4_conc) ='Precp_CaHPO4_conc';varlnml(cid_Precp_CaHPO4_conc)='non-band soil precipitated CaHPO4'
  unitl(cid_Precp_CaHPO4_conc)='mol m-3';vartypes(cid_Precp_CaHPO4_conc)=var_state_type

  varl(cid_PrecpB_CaHPO4_conc) ='PrecpB_CaHPO4_conc';varlnml(cid_PrecpB_CaHPO4_conc)='band soil precipitated CaHPO4'
  unitl(cid_PrecpB_CaHPO4_conc)='mol m-3';vartypes(cid_PrecpB_CaHPO4_conc)=var_state_type

  varl(cid_Precp_Ca5P3O12O3H3_conc) ='Precp_Ca5P3O12O3H3_conc';varlnml(cid_Precp_Ca5P3O12O3H3_conc)='non-band soil precipitated Ca5(PO4)3OH hydroxyapatite'
  unitl(cid_Precp_Ca5P3O12O3H3_conc)='mol m-3';vartypes(cid_Precp_Ca5P3O12O3H3_conc)=var_state_type

  varl(cid_PrecpB_Ca5P3O12O3H3_conc) ='PrecpB_Ca5P3O12O3H3_conc';varlnml(cid_PrecpB_Ca5P3O12O3H3_conc)='band soil precipitated Ca5(PO4)3OH hydroxyapatite'
  unitl(cid_PrecpB_Ca5P3O12O3H3_conc)='mol m-3';vartypes(cid_PrecpB_Ca5P3O12O3H3_conc)=var_state_type

  varl(cid_Precp_CaH4P2O8_conc) ='Precp_CaH4P2O8_conc';varlnml(cid_Precp_CaH4P2O8_conc)='non-band soil precipitated Ca(H2PO4)2'
  unitl(cid_Precp_CaH4P2O8_conc)='mol m-3';vartypes(cid_Precp_CaH4P2O8_conc)=var_state_type

  varl(cid_PrecpB_CaH4P2O8_conc) ='PrecpB_CaH4P2O8_conc';varlnml(cid_PrecpB_CaH4P2O8_conc)='band soil precipitated CaH4P2O8'
  unitl(cid_PrecpB_CaH4P2O8_conc)='mol m-3';vartypes(cid_PrecpB_CaH4P2O8_conc)=var_state_type

  varl(cid_Precp_FePO4_conc) ='Precp_FePO4_conc';varlnml(cid_Precp_FePO4_conc)='non-band soil precipitated FePO4'
  unitl(cid_Precp_FePO4_conc)='mol m-3';vartypes(cid_Precp_FePO4_conc)=var_state_type

  varl(cid_PrecpB_FePO4_con) ='PrecpB_FePO4_con';varlnml(cid_PrecpB_FePO4_con)='band soil precipitated FePO4'
  unitl(cid_PrecpB_FePO4_con)='mol m-3';vartypes(cid_PrecpB_FePO4_con)=var_state_type

  varl(fid_TR_NH4_soil) = 'TR_NH4_soil';varlnml(fid_TR_NH4_soil)='non-band soil total solute NH4 transformation'
  unitl(fid_TR_NH4_soil)= 'mol d-2 h-1';vartypes(fid_TR_NH4_soil)=var_flux_type

  varl(fid_TR_NH4_band_soil) = 'TR_NH4_band_soil';varlnml(fid_TR_NH4_band_soil)='band soil total solute NH4 transformation'
  unitl(fid_TR_NH4_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_NH4_band_soil)=var_flux_type

  varl(fid_TR_NH3_soil_vr) = 'TR_NH3_soil_vr';varlnml(fid_TR_NH3_soil_vr)='non-band total solute NH3 transformation'
  unitl(fid_TR_NH3_soil_vr)= 'mol d-2 h-1';vartypes(fid_TR_NH3_soil_vr)=var_flux_type

  varl(fid_TR_NH3_band_soil) = 'TR_NH3_band_soil';varlnml(fid_TR_NH3_band_soil)='band soil total solute NH3 transformation'
  unitl(fid_TR_NH3_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_NH3_band_soil)=var_flux_type

  varl(fid_TR_H1PO4_soil) = 'TR_H1PO4_soil';varlnml(fid_TR_H1PO4_soil)='non-band soil total solute HPO4 transformation'
  unitl(fid_TR_H1PO4_soil) = 'mol d-2 h-1';vartypes(fid_TR_H1PO4_soil)=var_flux_type

  varl(fid_TR_H2PO4_soil) = 'TR_H2PO4_soil';varlnml(fid_TR_H2PO4_soil)='non-band soil total solute H2PO4 transformation'
  unitl(fid_TR_H2PO4_soil)= 'mol d-2 h-1';vartypes(fid_TR_H2PO4_soil)=var_flux_type

  varl(fid_TR_H1PO4_band_soil) = 'TR_H1PO4_band_soil';varlnml(fid_TR_H1PO4_band_soil)='band soil total solute HPO4 transformation'
  unitl(fid_TR_H1PO4_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_H1PO4_band_soil)=var_flux_type

  varl(fid_TR_H2PO4_band_soil) = 'TR_H2PO4_band_soil';varlnml(fid_TR_H2PO4_band_soil)='band soil total solute H2PO4 transformation'
  unitl(fid_TR_H2PO4_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_H2PO4_band_soil)=var_flux_type

  varl(fid_TR_NH4_sorbed_soil) = 'TR_NH4_sorbed_soil';varlnml(fid_TR_NH4_sorbed_soil)='non-band soil total adsorbed NH4 transformation'
  unitl(fid_TR_NH4_sorbed_soil)= 'mol d-2 h-1';vartypes(fid_TR_NH4_sorbed_soil)=var_flux_type

  varl(fid_TR_NH4_sorbed_band_soil) = 'TR_NH4_sorbed_band_soil';varlnml(fid_TR_NH4_sorbed_band_soil)='band soil total adsorbed NH4 transformation'
  unitl(fid_TR_NH4_sorbed_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_NH4_sorbed_band_soil)=var_flux_type

  varl(fid_TR_ROH_sorbed_soil) = 'TR_ROH_sorbed_soil';varlnml(fid_TR_ROH_sorbed_soil)='non-band soil total adsorbed OH transformation'
  unitl(fid_TR_ROH_sorbed_soil)= 'mol d-2 h-1';vartypes(fid_TR_ROH_sorbed_soil)=var_flux_type

  varl(fid_TR_ROH2_sorbed_soil) = 'TR_ROH2_sorbed_soil';varlnml(fid_TR_ROH2_sorbed_soil)='non-band soil total adsorbed OH2 transformation'
  unitl(fid_TR_ROH2_sorbed_soil)= 'mol d-2 h-1';vartypes(fid_TR_ROH2_sorbed_soil)=var_flux_type

  varl(fid_TR_RHPO4_sorbed_soil) = 'TR_RHPO4_sorbed_soil';varlnml(fid_TR_RHPO4_sorbed_soil)='non-band soil total adsorbed HPO4 transformation'
  unitl(fid_TR_RHPO4_sorbed_soil)= 'mol d-2 h-1';vartypes(fid_TR_RHPO4_sorbed_soil)=var_flux_type

  varl(fid_TR_RH2PO4_sorbed_soil) = 'TR_RH2PO4_sorbed_soil';varlnml(fid_TR_RH2PO4_sorbed_soil)='non-band soil total adsorbed H2PO4 transformation'
  unitl(fid_TR_RH2PO4_sorbed_soil)= 'mol d-2 h-1';vartypes(fid_TR_RH2PO4_sorbed_soil)=var_flux_type

  varl(fid_TR_ROH_sorbed_band_soil) = 'TR_ROH_sorbed_band_soil';varlnml(fid_TR_ROH_sorbed_band_soil)='band soil total adsorbed OH transformation'
  unitl(fid_TR_ROH_sorbed_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_ROH_sorbed_band_soil)=var_flux_type

  varl(fid_TR_ROH2_sorbed_band_soil) = 'TR_ROH2_sorbed_band_soil';varlnml(fid_TR_ROH2_sorbed_band_soil)='band soil total adsorbed OH2 transformation'
  unitl(fid_TR_ROH2_sorbed_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_ROH2_sorbed_band_soil)=var_flux_type

  varl(fid_TR_RHPO4_sorbed_band_soil) = 'TR_RHPO4_sorbed_band_soil';varlnml(fid_TR_RHPO4_sorbed_band_soil)='band soil total adsorbed HPO4 transformation'
  unitl(fid_TR_RHPO4_sorbed_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_RHPO4_sorbed_band_soil)=var_flux_type

  varl(fid_TR_RH2PO4_sorbed_band_soil) = 'TR_RH2PO4_sorbed_band_soil';varlnml(fid_TR_RH2PO4_sorbed_band_soil)='band soil total adsorbed H2PO4 transformation'
  unitl(fid_TR_RH2PO4_sorbed_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_RH2PO4_sorbed_band_soil)=var_flux_type

  varl(fid_TR_AlPO4_precip_soil)= 'TR_AlPO4_precip_soil';varlnml(fid_TR_AlPO4_precip_soil)='non-band total precipitated AlPO4 transformation'
  unitl(fid_TR_AlPO4_precip_soil)= 'mol d-2 h-1';vartypes(fid_TR_AlPO4_precip_soil)=var_flux_type

  varl(fid_TR_FePO4_precip_soil)= 'TR_FePO4_precip_soil';varlnml(fid_TR_FePO4_precip_soil)='non-band soil total precipitated FePO4 transformation'
  unitl(fid_TR_FePO4_precip_soil)= 'mol d-2 h-1';vartypes(fid_TR_FePO4_precip_soil)=var_flux_type

  varl(fid_TR_CaHPO4_precip_soil)= 'TR_CaHPO4_precip_soil';varlnml(fid_TR_CaHPO4_precip_soil)='non-band soil total precipitated CaHPO4 transformation'
  unitl(fid_TR_CaHPO4_precip_soil)= 'mol d-2 h-1';vartypes(fid_TR_CaHPO4_precip_soil)=var_flux_type

  varl(fid_TR_apatite_precip_soil)= 'TR_apatite_precip_soil';varlnml(fid_TR_apatite_precip_soil)='non-band soil total precipitated CaH4P2O8 transformation'
  unitl(fid_TR_apatite_precip_soil)= 'mol d-2 h-1';vartypes(fid_TR_apatite_precip_soil)=var_flux_type

  varl(fid_TR_CaH4P2O8_precip_soil)= 'TR_CaH4P2O8_precip_soil';varlnml(fid_TR_CaH4P2O8_precip_soil)='non-band soil total precipitated apatite transformation'
  unitl(fid_TR_CaH4P2O8_precip_soil)= 'mol d-2 h-1';vartypes(fid_TR_CaH4P2O8_precip_soil)=var_flux_type

  varl(fid_TR_AlPO4_precip_band_soil)= 'TR_AlPO4_precip_band_soil';varlnml(fid_TR_AlPO4_precip_band_soil)='band soil total precipitated AlPO4 transformation'
  unitl(fid_TR_AlPO4_precip_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_AlPO4_precip_band_soil)=var_flux_type

  varl(fid_TR_FePO4_precip_band_soil)= 'TR_FePO4_precip_band_soil';varlnml(fid_TR_FePO4_precip_band_soil)='band soil total precipitated FePO4 transformation'
  unitl(fid_TR_FePO4_precip_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_FePO4_precip_band_soil)=var_flux_type

  varl(fid_TR_CaHPO4_precip_band_soil)= 'TR_CaHPO4_precip_band_soil';varlnml(fid_TR_CaHPO4_precip_band_soil)='band soil total precipitated CaHPO4 transformation'
  unitl(fid_TR_CaHPO4_precip_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_CaHPO4_precip_band_soil)=var_flux_type

  varl(fid_TR_apatite_precip_band_soil)= 'TR_apatite_precip_band_soil';varlnml(fid_TR_apatite_precip_band_soil)='band soil total precipitated CaH4P2O8 transformation'
  unitl(fid_TR_apatite_precip_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_apatite_precip_band_soil)=var_flux_type

  varl(fid_TR_CaH4P2O8_precip_band_soil)= 'TR_CaH4P2O8_precip_band_soil';varlnml(fid_TR_CaH4P2O8_precip_band_soil)='band soil total precipitated apatite transformation'
  unitl(fid_TR_CaH4P2O8_precip_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_CaH4P2O8_precip_band_soil)=var_flux_type

  varl(fid_TR_Al_3p_soil)  = 'TR_Al_3p_soil';varlnml(fid_TR_Al_3p_soil)='non-band soil total solute Al transformation'
  unitl(fid_TR_Al_3p_soil) = 'mol d-2 h-1';vartypes(fid_TR_Al_3p_soil)=var_flux_type

  end subroutine getvarlist_nosalt

end module ChemIDMod
