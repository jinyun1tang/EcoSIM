module AquaSaltChemMod
  use  MiniMathMod   , only : addone
  use ModelStatusType, only : model_status_type
  use data_kind_mod  , only : r8 => DAT_KIND_R8
  use ChemIDMod
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__


  public :: Init_geochem_salt
  public :: getvarlist_salt
  public :: initmodel_salt
  public :: RunModel_salt
  contains

! ----------------------------------------------------------------------
  subroutine RunModel_salt(nvars,ystates0l, ystatesfl, err_status)
  use SoluteChemDataType, only : chem_var_type, solute_flx_type
  use SaltChemEquilibriaMod
  implicit none
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  real(r8), intent(out) :: ystatesfl(nvars)
  type(model_status_type), intent(out) :: err_status
! Gapon selectivity coefficient

  type(solute_flx_type) :: solflx
  type(chem_var_type)   :: chemvar

  integer :: jj

  call err_status%reset()

  call SetChemVar(nvars, ystates0l, chemvar)

  call SaltChemEquilibria(chemvar,solflx)

  call RetrieveYstatef(nvars, ystates0l, ystatesfl, chemvar, solflx)


  end subroutine RunModel_salt

! ----------------------------------------------------------------------

  subroutine initmodel_salt(nvars, ystatesfl,err_status)

  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(inout) :: ystatesfl(nvars)
  type(model_status_type), intent(out) :: err_status


  end subroutine initmodel_salt
!------------------------------------------------------------------
  subroutine Init_geochem_salt(nvars)
  implicit none
  integer, intent(out) :: nvars
  integer :: itemp
  itemp=0
  cid_CO2S     =addone(itemp)  !aqueous CO2  micropore	[g d-2]
  cid_H1PO4_2e_conc    =addone(itemp)  !soil aqueous HPO4 content micropore non-band, [mol m-3]
  cid_H2PO4_1e_conc    =addone(itemp)  !soil aqueous H2PO4 content micropore non-band, [mol m-3]
  cid_NH3_aqu_conc     =addone(itemp)  !soil NH3 concentration in non-band soil, [mol m-3]
  cid_NH4_1p_conc     =addone(itemp)  !soil NH4 concentration in non-band soil, [mol m-3]
  cid_ZNO3S    =addone(itemp)  !NO3 mass non-band micropore, [g d-2]
  cid_ZHY      =addone(itemp)  !soil aqueous H content micropore, [mol d-2]
  cid_ZOH      =addone(itemp)  !soil aqueous OH content micropore, [mol d-2]
  cid_ZAL      =addone(itemp)  !soil aqueous Al content micropore, [mol d-2]
  cid_ZFE      =addone(itemp)  !soil aqueous Fe content micropore, [mol d-2]
  cid_ZCA      =addone(itemp)  !soil aqueous Ca content micropore, [mol d-2]
  cid_ZMG      =addone(itemp)  !soil aqueous Mg content micropore, [mol d-2]
  cid_ZNA      =addone(itemp)  !soil aqueous Na content micropore, [mol d-2]
  cid_ZKA      =addone(itemp)  !soil aqueous K content micropore, [mol d-2]
  cid_ZSO4     =addone(itemp)  !soil aqueous SO4 content micropore, [mol d-2]
  cid_ZCL      =addone(itemp)  !soil aqueous Cl content micropore, [mol d-2]
  cid_ZCO3     =addone(itemp)  !soil aqueous CO3 content micropore, [mol d-2]
  cid_ZHCO3    =addone(itemp)  !soil aqueous HCO3 content micropore, [mol d-2]
  cid_ZALOH1   =addone(itemp)  !soil aqueous AlOH content micropore, [mol d-2]
  cid_ZALOH2   =addone(itemp)  !soil aqueous AlOH2 content micropore, [mol d-2]
  cid_ZALOH3   =addone(itemp)  !soil aqueous AlOH3 content micropore, [mol d-2]
  cid_ZALOH4   =addone(itemp)  !soil aqueous AlOH4 content micropore, [mol d-2]
  cid_ZALS     =addone(itemp)  !soil aqueous AlSO4 content micropore, [mol d-2]
  cid_ZFEOH1   =addone(itemp)  !soil aqueous FeOH content micropore, [mol d-2]
  cid_ZFEOH2   =addone(itemp)  !soil aqueous FeOH2 content micropore, [mol d-2]
  cid_ZFEOH3   =addone(itemp)  !soil aqueous FeOH3 content micropore, [mol d-2]
  cid_ZFEOH4   =addone(itemp)  !soil aqueous FeOH4 content micropore, [mol d-2]
  cid_ZFES     =addone(itemp)  !soil aqueous FeSO4 content micropore, [mol d-2]
  cid_ZCAO     =addone(itemp)  !soil aqueous CaOH2 content micropore, [mol d-2]
  cid_ZCAC     =addone(itemp)  !soil aqueous CACO3 content micropore, [mol d-2]
  cid_ZCAH     =addone(itemp)  !soil aqueous CaHCO3 content micropore, [mol d-2]
  cid_ZCAS     =addone(itemp)  !soil aqueous CaSO4 content micropore, [mol d-2]
  cid_ZMGO     =addone(itemp)  !soil aqueous MgOH content micropore, [mol d-2]
  cid_ZMGC     =addone(itemp)  !soil aqueous MgCO3 content micropore, [mol d-2]
  cid_ZMGH     =addone(itemp)  !soil aqueous MgHCO3 content micropore, [mol d-2]
  cid_ZMGS     =addone(itemp)  !soil aqueous MgSO4 content micropore, [mol d-2]
  cid_ZNAC     =addone(itemp)  !soil aqueous NaCO3 content micropore, [mol d-2]
  cid_ZNAS     =addone(itemp)  !soil aqueous NaSO4 content micropore, [mol d-2]
  cid_ZKAS     =addone(itemp)  !soil aqueous KSO4 content micropore, [mol d-2]
  cid_H0PO4    =addone(itemp)  !soil aqueous PO4 content micropore non-band, [mol d-2]
  cid_H3PO4    =addone(itemp)  !soil aqueous H3PO4 content micropore non-band, [mol d-2]
  cid_ZFE1P    =addone(itemp)  !soil aqueous FeHPO4 content micropore non-band, [mol d-2]
  cid_ZFE2P    =addone(itemp)  !soil aqueous FeH2PO4 content micropore non-band, [mol d-2]
  cid_ZCA0P    =addone(itemp)  !soil aqueous CaPO4 content micropore non-band, [mol d-2]
  cid_ZCA1P    =addone(itemp)  !soil aqueous CaHPO4 content micropore non-band, [mol d-2]
  cid_ZCA2P    =addone(itemp)  !soil aqueous CaH4P2O8 content micropore non-band, [mol d-2]
  cid_ZMG1P    =addone(itemp)  !soil aqueous MgHPO4 content micropore non-band, [mol d-2]
  cid_Precp_AlPO4_conc   =addone(itemp)  !precipitated AlPO4 non-band, [mol m-3]
  cid_Precp_CaHPO4_conc   =addone(itemp)  !precipitated CaHPO4 non-band soil, [mol m-3]
  cid_Precp_Ca5P3O12O3H3_conc   =addone(itemp)  !precipitated Ca5(PO4)3OH hydroxyapatite non-band soil, [mol m-3]
  cid_Precp_CaH4P2O8_conc   =addone(itemp)  !precipitated Ca(H2PO4)2 non-band soil, [mol m-3]
  cid_Precp_FePO4_conc   =addone(itemp)  !precipitated FePO4 non-band soil, [mol m-3]
  cid_PALOH    =addone(itemp)  !precipitated Al(OH)3, [mol d-2]
  cid_PFEOH    =addone(itemp)  !precipitated Fe(OH)3, [mol d-2]
  cid_PCACO    =addone(itemp)  !precipitated CaCO3, [mol d-2]
  cid_PCASO    =addone(itemp)  !precipitated CaSO4, [mol d-2]
  cid_XHY      =addone(itemp)  !exchangeable H(+) , [mol d-2]
  cid_XAL      =addone(itemp)  !exchangeable Al, [mol d-2]
  cid_XFE      =addone(itemp)  !exchangeable Fe, [mol d-2]
  cid_XCA      =addone(itemp)  !exchangeable Ca, [mol d-2]
  cid_XMG      =addone(itemp)  !exchangeable Mg, [mol d-2]
  cid_XNA      =addone(itemp)  !exchangeable Na, [mol d-2]
  cid_XKA      =addone(itemp)  !exchangeable K, [mol d-2]
  cid_XHC      =addone(itemp)  !exchangeable COOH , [mol d-2]
  cid_XROH1_conc    =addone(itemp)  !exchangeable OH  non-band, [mol d-2]
  cid_XALO2    =addone(itemp)  !exchangeable AlOH2 , [mol d-2]
  cid_XFEO2    =addone(itemp)  !exchangeable Fe(OH)2, [mol d-2]
  cid_XNH4_conc     =addone(itemp)  !exchangeable NH4 non-band soil, [mol d-2]
  cid_XOH_conc    =addone(itemp)  !exchangeable OH- non-band, [mol d-2]
  cid_XHPO4_conc    =addone(itemp)  !exchangeable HPO4  non-band, [mol m-3]
  cid_XROH2_conc    =addone(itemp)  !exchangeable OH2  non-band soil, [mol m-3]
  cid_XH2PO4_conc    =addone(itemp)  !exchangeable H2PO4  non-band soil, [mol m-3]

  cid_H1PO4_2e_band_conc    =addone(itemp)  !soil aqueous HPO4 content micropore band, [mol m-3]
  cid_H2PO4_1e_band_conc    =addone(itemp)  !soil aqueous H2PO4 content micropore  band, [mol m-3]
  cid_NH3_aqu_band_conc     =addone(itemp)  !soil NH3 concentration in band soil, [mol m-3]
  cid_NH4_1p_band_conc     =addone(itemp)  !soil NH4 concentration in band soil, [mol m-3]
  cid_ZNO3B    =addone(itemp)  !NO3 mass band micropore, [g d-2]
  cid_H0POB    =addone(itemp)  !soil aqueous PO4 content micropore band, [mol d-2]
  cid_H3POB    =addone(itemp)  !soil aqueous H3PO4 content micropore band, [mol d-2]
  cid_ZFE1PB   =addone(itemp)  !soil aqueous FeHPO4 content micropore band, [mol d-2]
  cid_ZFE2PB   =addone(itemp)  !soil aqueous FeH2PO4 content micropore band, [mol d-2
  cid_ZCA0PB   =addone(itemp)  !soil aqueous CaPO4 content micropore band, [mol d-2]
  cid_ZCA1PB   =addone(itemp)  !soil aqueous CaHPO4 content micropore band, [mol d-2]
  cid_ZCA2PB   =addone(itemp)  !soil aqueous CaH4P2O8 content micropore band, [mol d-2]
  cid_ZMG1PB   =addone(itemp)  !soil aqueous MgHPO4 content micropore band, [mol d-2]
  cid_PrecpB_AlPO4_conc   =addone(itemp)  !precipitated AlPO4 band soil, [mol m-3]
  cid_PrecpB_CaHPO4_conc   =addone(itemp)  !precipitated CaHPO4 band soil, [mol m-3]
  cid_PrecpB_Ca5P3O12O3H3_conc   =addone(itemp)  !precipitated Ca5(PO4)3OH hydroxyapatite band soil, [mol m-3]
  cid_PrecpB_CaH4P2O8_conc   =addone(itemp)  !precipitated CaH4P2O8 band soil, [mol m-3]
  cid_PrecpB_FePO4_con   =addone(itemp)  !precipitated FePO4 band soil, [mol m-3]
  cid_XNH4_band_conc     =addone(itemp)  !exchangeable NH4 band soil, [mol d-2]
  cid_XROH1_band_conc    =addone(itemp)  !exchangeable OH- band, [mol d-2]
  cid_XHPO4_band_conc    =addone(itemp)  !exchangeable HPO4 concentration band-soil, [mol m-3]
  cid_XH2PO4_band_conc    =addone(itemp)  !exchangeable H2PO4 concentration band-soil, [mol m-3]
  cid_XROH_band_conc    =addone(itemp)  !exchangeable OH band-soil, [mol m-3]
  cid_XROH2_band_conc    =addone(itemp)  !exchangeable OH2 band-soil, [mol m-3]

  fid_TR_CaCO3_precip_soil   =addone(itemp)
  fid_TR_NaCO3_soil    =addone(itemp)
  fid_TR_MgCO3_soil    =addone(itemp)
  fid_TR_CaCO3_soil    =addone(itemp)
  fid_TR_MgHCO3_soil    =addone(itemp)
  fid_TR_CaHCO3_soil    =addone(itemp)
  fid_TR_HCO3    =addone(itemp)
  fid_TR_HCO3_sorbed_soil    =addone(itemp)
  fid_TR_CaH4P2O8_soil    =addone(itemp)
  fid_TR_FeH2PO4_soil    =addone(itemp)
  fid_TR_MgHPO4_soil    =addone(itemp)
  fid_TR_CaHPO4_soil    =addone(itemp)
  fid_TR_FeHPO4_soil    =addone(itemp)
  fid_TR_CaPO4_soil    =addone(itemp)
  fid_TR_PO4_soil    =addone(itemp)
  fid_TR_FePO4_precip_soil   =addone(itemp)
  fid_TR_AlPO4_precip_soil   =addone(itemp)
  fid_TR_CaH4P2O8_precip_soil   =addone(itemp)
  fid_TR_CaHPO4_precip_soil   =addone(itemp)
  fid_TR_apatite_precip_soil   =addone(itemp)
  fid_TR_FeOH_soil    =addone(itemp)
  fid_TR_AlOH_soil    =addone(itemp)
  fid_TR_FeO2H2_soil    =addone(itemp)
  fid_TR_AlO2H2_soil    =addone(itemp)
  fid_TR_FeO3H3_soil    =addone(itemp)
  fid_TR_AlO3H3_soil    =addone(itemp)
  fid_TR_FeO4H4_soil    =addone(itemp)
  fid_TR_AlO4H4_soil    =addone(itemp)
  fid_TR_MgOH_soil    =addone(itemp)
  fid_TR_CaOH_soil    =addone(itemp)
  fid_TR_FeOH3_precip_soil   =addone(itemp)
  fid_TR_AlOH3_precip_soil   =addone(itemp)
  fid_TR_FeO2H2_sorbed_soil   =addone(itemp)
  fid_TR_Al_3p_soil     =addone(itemp)
  fid_TR_AlSO4_soil   =addone(itemp)
  fid_TR_RHPO4_sorbed_band_soil   =addone(itemp)
  fid_TR_RH2PO4_sorbed_band_soil   =addone(itemp)
  fid_TR_RO_sorbed_band_soil   =addone(itemp)
  fid_TR_ROH_sorbed_band_soil   =addone(itemp)
  fid_TR_ROH2_sorbed_band_soil   =addone(itemp)
  fid_TR_Ca_2p_soil    =addone(itemp)
  fid_TR_CaSO4_soil   =addone(itemp)
  fid_TR_CaSO4_precip_soil  =addone(itemp)
  fid_TR_CO2_gchem_soil_vr   =addone(itemp)
  fid_TR_CO3_2e_soil   =addone(itemp)
  fid_TR_Fe_3p_soil    =addone(itemp)
  fid_TR_FeSO4_soil   =addone(itemp)

  fid_TR_H1PO4_soil   =addone(itemp)
  fid_TR_H2PO4_soil   =addone(itemp)
  fid_TR_H3PO4_sorbed_soil   =addone(itemp)
  fid_TR_H_p_soil    =addone(itemp)
  fid_TR_K_1p_soil    =addone(itemp)
  fid_TR_KSO4_soil   =addone(itemp)
  fid_TR_Mg_2p_soil    =addone(itemp)
  fid_TR_MgSO4_soil   =addone(itemp)
  fid_TR_NH3_soil_vr   =addone(itemp)
  fid_TR_NH4_soil   =addone(itemp)
  fid_TR_Na_p_soil    =addone(itemp)
  fid_TR_NaSO4_soil   =addone(itemp)
  fid_TR_OH_1e_soil    =addone(itemp)
  fid_TR_SO4_2e_soil   =addone(itemp)
  fid_TR_RHPO4_sorbed_soil   =addone(itemp)
  fid_TR_RH2PO4_sorbed_soil   =addone(itemp)
  fid_TR_Al_sorbed_soil   =addone(itemp)
  fid_TR_AlO2H2_sorbed_soil  =addone(itemp)
  fid_TR_Ca_sorbed_soil   =addone(itemp)
  fid_TR_Fe_sorbed_soil   =addone(itemp)
  fid_TR_RO_sorbed_soil   =addone(itemp)
  fid_TR_ROH_sorbed_soil   =addone(itemp)
  fid_TR_ROH2_sorbed_soil   =addone(itemp)
  fid_TR_H_p_sorbed_soil   =addone(itemp)
  fid_TR_K_sorbed_soil   =addone(itemp)
  fid_TR_Mg_sorbed_soil   =addone(itemp)
  fid_TR_NH4_sorbed_soil   =addone(itemp)
  fid_TR_Na_sorbed_soil   =addone(itemp)
  fid_Txchem_CO2   =addone(itemp)
  fid_TBION   =addone(itemp)
  fid_TRH2O   =addone(itemp)

  fid_TR_CaH4P2O8_band_soil    =addone(itemp)
  fid_TR_FeH2PO4_band_soil    =addone(itemp)
  fid_TR_MgHPO4_band_soil    =addone(itemp)
  fid_TR_CaHPO4_band_soil    =addone(itemp)
  fid_TR_FeHPO4_band_soil    =addone(itemp)
  fid_TR_CaPO4_band_soil    =addone(itemp)
  fid_TR_PO4_band_soil    =addone(itemp)
  fid_TR_FePO4_precip_band_soil   =addone(itemp)
  fid_TR_AlPO4_precip_band_soil   =addone(itemp)
  fid_TR_CaH4P2O8_precip_band_soil   =addone(itemp)
  fid_TR_CaHPO4_precip_band_soil   =addone(itemp)
  fid_TR_apatite_precip_band_soil   =addone(itemp)
  fid_TR_H1PO4_band_soil   =addone(itemp)
  fid_TR_H2PO4_band_soil   =addone(itemp)
  fid_TR_H3PO4_band_soil   =addone(itemp)
  fid_TR_NH3_band_soil   =addone(itemp)
  fid_TR_NH4_band_soil   =addone(itemp)
  fid_TR_NH4_sorbed_band_soil   =addone(itemp)
  nvars=itemp
  end subroutine Init_geochem_salt


!------------------------------------------------------------------
  subroutine getvarlist_salt(nvars, varl, varlnml, unitl, vartypes)

  use bhistMod, only : hist_var_str_len,hist_unit_str_len,hist_var_lon_str_len
  use fileUtil, only :  var_flux_type, var_state_type
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(out) :: varl(:)     !variable name
  character(len=hist_var_lon_str_len), intent(out) :: varlnml(:)     !variable name
  character(len=hist_unit_str_len),intent(out) :: unitl(:)
  integer                         ,intent(out) :: vartypes(:)


  varl(cid_CO2S)='CO2S';varlnml(cid_CO2S)='aqueous CO2 concentration micropore';
  unitl(cid_CO2S)='gC d-2';vartypes(cid_CO2S)=var_state_type

  varl(cid_H1PO4_2e_conc)='H1PO4_2e_conc';varlnml(cid_H1PO4_2e_conc)='non-band soil aqueous HPO4 content micropore'
  unitl(cid_H1PO4_2e_conc)='mol m-3';vartypes(cid_H1PO4_2e_conc)=var_state_type

  varl(cid_H2PO4_1e_conc)='H2PO4_1e_conc';varlnml(cid_H2PO4_1e_conc)='non-band soil micropore aqueous H2PO4 content'
  unitl(cid_H2PO4_1e_conc)='mol m-3';vartypes(cid_H2PO4_1e_conc)=var_state_type

  varl(cid_NH3_aqu_conc) ='NH3_aqu_conc';varlnml(cid_NH3_aqu_conc)='non-band soil NH3 concentration'
  unitl(cid_NH3_aqu_conc)='mol m-3';vartypes(cid_NH3_aqu_conc)=var_state_type

  varl(cid_NH4_1p_conc)  ='NH4_1p_conc';varlnml(cid_NH4_1p_conc)='non-band soil NH4 concentration'
  unitl(cid_NH4_1p_conc)='mol m-3';vartypes(cid_NH4_1p_conc)=var_state_type

  varl(cid_ZNO3S) = 'ZNO3S';varlnml(cid_ZNO3S)='non-band soil micropore NO3 mass'
  unitl(cid_ZNO3S)= 'g d-2';vartypes(cid_ZNO3S)=var_state_type

  varl(cid_ZHY) ='ZHY'; varlnml(cid_ZHY) ='soil micropore aqueous H content'
  unitl(cid_ZHY)='mol d-2';vartypes(cid_ZHY)=var_state_type

  varl(cid_ZOH) ='ZOH'; varlnml(cid_ZOH)='soil micropore aqueous OH content micropore'
  unitl(cid_ZOH)='mol d-2';vartypes(cid_ZOH)=var_state_type

  varl(cid_ZAL) ='ZAL'; varlnml(cid_ZAL)= 'soil micropore aqueous Al content micropore'
  unitl(cid_ZAL)='mol d-2';vartypes(cid_ZAL)=var_state_type

  varl(cid_ZFE) ='ZFE'; varlnml(cid_ZFE)='soil micropore aqueous Fe content'
  unitl(cid_ZFE)='mol d-2';vartypes(cid_ZFE)=var_state_type

  varl(cid_ZCA) ='ZCA'; varlnml(cid_ZFE)='soil micropore aqueous Ca content'
  unitl(cid_ZFE)='mol d-2'; vartypes(cid_ZCA)=var_state_type

  varl(cid_ZMG) ='ZMG';varlnml(cid_ZMG)='soil aqueous Mg content micropore'
  unitl(cid_ZMG)='mol d-2';vartypes(cid_ZMG)=var_state_type

  varl(cid_ZNA) ='ZNA';varlnml(cid_ZNA)='soil aqueous Na content micropore'
  unitl(cid_ZNA) ='mol d-2';vartypes(cid_ZNA)=var_state_type

  varl(cid_ZKA) ='ZKA';varlnml(cid_ZKA)='soil aqueous K content micropore'
  unitl(cid_ZKA) ='mol d-2';vartypes(cid_ZKA)=var_state_type

  varl(cid_ZSO4)='ZSO4';varlnml(cid_ZSO4)='soil aqueous micropore SO4 content'
  unitl(cid_ZSO4)='mol d-2';vartypes(cid_ZSO4)=var_state_type

  varl(cid_ZCL) ='ZCL';varlnml(cid_ZCL)= 'soil aqueous micropore Cl content'
  unitl(cid_ZCL)='mol d-2';vartypes(cid_ZCL)=var_state_type

  varl(cid_ZCO3)='ZCO3';varlnml(cid_ZCO3)='soil aqueous CO3 content micropore'
  unitl(cid_ZCO3)='mol d-2';vartypes(cid_ZCO3)=var_state_type

  varl(cid_ZHCO3)='ZHCO3';varlnml(cid_ZHCO3)='soil aqueous HCO3 content micropore'
  unitl(cid_ZHCO3)='mol d-2';vartypes(cid_ZHCO3)=var_state_type

  varl(cid_ZALOH1)='ZALOH1';varlnml(cid_ZALOH1)='soil aqueous AlOH content micropore'
  unitl(cid_ZALOH1)='mol d-2';vartypes(cid_ZALOH1)=var_state_type

  varl(cid_ZALOH2)='ZALOH2';varlnml(cid_ZALOH2)='soil aqueous AlOH2 content micropore'
  unitl(cid_ZALOH2)='mol d-2';vartypes(cid_ZALOH2)=var_state_type

  varl(cid_ZALOH3)='ZALOH3'; varlnml(cid_ZALOH3)='soil aqueous AlOH3 content micropore'
  unitl(cid_ZALOH3)='mol d-2';vartypes(cid_ZALOH3)=var_state_type

  varl(cid_ZALOH4)='ZALOH4'; varlnml(cid_ZALOH4)='soil aqueous AlOH4 content micropore'
  unitl(cid_ZALOH4)='mol d-2';vartypes(cid_ZALOH4)=var_state_type

  varl(cid_ZALS)='ZALS'; varlnml(cid_ZALS)='soil aqueous AlSO4 content micropore'
  unitl(cid_ZALS)='mol d-2'; vartypes(cid_ZALS)=var_state_type

  varl(cid_ZFEOH1)='ZFEOH1';varlnml(cid_ZFEOH1)='soil aqueous FeOH content micropore'
  unitl(cid_ZFEOH1)='mol d-2';vartypes(cid_ZFEOH1)=var_state_type

  varl(cid_ZFEOH2)='ZFEOH2';varlnml(cid_ZFEOH2)='soil aqueous FeOH2 content micropore'
  unitl(cid_ZFEOH2)='mol d-2';vartypes(cid_ZFEOH2)=var_state_type

  varl(cid_ZFEOH3)='ZFEOH3';varlnml(cid_ZFEOH3)='soil aqueous FeOH3 content micropore'
  unitl(cid_ZFEOH3)='mol d-2';vartypes(cid_ZFEOH3)=var_state_type

  varl(cid_ZFEOH4)='ZFEOH4';varlnml(cid_ZFEOH4)='soil aqueous FeOH4 content micropore'
  unitl(cid_ZFEOH4)='mol d-2';vartypes(cid_ZFEOH4)=var_state_type

  varl(cid_ZFES)='ZFES';varlnml(cid_ZFES)='soil aqueous FeSO4 content micropore'
  unitl(cid_ZFES)='mol d-2';vartypes(cid_ZFES)=var_state_type

  varl(cid_ZCAO)='ZCAO';varlnml(cid_ZCAO)='soil aqueous CaOH2 content micropore'
  unitl(cid_ZCAO)='mol d-2';vartypes(cid_ZCAO)=var_state_type

  varl(cid_ZCAC)='ZCAC';varlnml(cid_ZCAC)='soil aqueous CACO3 content micropore'
  unitl(cid_ZCAC)='mol d-2';vartypes(cid_ZCAC)=var_state_type

  varl(cid_ZCAH)='ZCAH';varlnml(cid_ZCAH)='soil aqueous CaHCO3 content micropore'
  unitl(cid_ZCAH)='mol d-2';vartypes(cid_ZCAH)=var_state_type

  varl(cid_ZCAS)='ZCAS';varlnml(cid_ZCAS)='soil aqueous CaSO4 content micropore'
  unitl(cid_ZCAS)='mol d-2';vartypes(cid_ZCAS)=var_state_type

  varl(cid_ZMGO)='ZMGO';varlnml(cid_ZMGO)='soil aqueous MgOH content micropore'
  unitl(cid_ZMGO)='mol d-2';vartypes(cid_ZMGO)=var_state_type

  varl(cid_ZMGC)='ZMGC';varlnml(cid_ZMGC)='soil aqueous MgCO3 content micropore'
  unitl(cid_ZMGC)='mol d-2';vartypes(cid_ZMGC)=var_state_type

  varl(cid_ZMGH)='ZMGH';varlnml(cid_ZMGH)='soil aqueous MgHCO3 content micropore'
  unitl(cid_ZMGH)='mol d-2';vartypes(cid_ZMGH)=var_state_type

  varl(cid_ZMGS)='ZMGS';varlnml(cid_ZMGS)='soil aqueous MgSO4 content micropore'
  unitl(cid_ZMGS)='mol d-2';vartypes(cid_ZMGS)=var_state_type

  varl(cid_ZNAC)='ZNAC';varlnml(cid_ZNAC)='soil aqueous NaCO3 content micropore'
  unitl(cid_ZNAC)='mol d-2';vartypes(cid_ZNAC)=var_state_type

  varl(cid_ZNAS)='ZNAS';varlnml(cid_ZNAS)='soil aqueous NaSO4 content micropore'
  unitl(cid_ZNAS)='mol d-2';vartypes(cid_ZNAS)=var_state_type

  varl(cid_ZKAS)='ZKAS';varlnml(cid_ZKAS)='soil aqueous KSO4 content micropore'
  unitl(cid_ZKAS)='mol d-2';vartypes(cid_ZKAS)=var_state_type

  varl(cid_H0PO4)='H0PO4';varlnml(cid_H0PO4)='soil aqueous PO4 content micropore non-band'
  unitl(cid_H0PO4)='mol d-2';vartypes(cid_H0PO4)=var_state_type

  varl(cid_H3PO4)='H3PO4';varlnml(cid_H3PO4)='soil aqueous H3PO4 content micropore non-band'
  unitl(cid_H3PO4)='mol d-2';vartypes(cid_H3PO4)=var_state_type

  varl(cid_ZFE1P)='ZFE1P';varlnml(cid_ZFE1P)='soil aqueous FeHPO4 content micropore non-band'
  unitl(cid_ZFE1P)='mol d-2';vartypes(cid_ZFE1P)=var_state_type

  varl(cid_ZFE2P)='ZFE2P';varlnml(cid_ZFE2P)='soil aqueous FeH2PO4 content micropore non-band'
  unitl(cid_ZFE2P)='mol d-2';vartypes(cid_ZFE2P)=var_state_type

  varl(cid_ZCA0P)='ZCA0P';varlnml(cid_ZCA0P)='soil aqueous CaPO4 content micropore non-band'
  unitl(cid_ZCA0P)='mol d-2';vartypes(cid_ZCA0P)=var_state_type

  varl(cid_ZCA1P)='ZCA1P'; varlnml(cid_ZCA1P)='soil aqueous CaHPO4 content micropore non-band'
  unitl(cid_ZCA1P)='mol d-2';vartypes(cid_ZCA1P)=var_state_type

  varl(cid_ZCA2P)='ZCA2P'; varlnml(cid_ZCA2P)='soil aqueous CaH4P2O8 content micropore non-band'
  unitl(cid_ZCA2P)='mol d-2';vartypes(cid_ZCA2P)=var_state_type

  varl(cid_ZMG1P)='ZMG1P';varlnml(cid_ZMG1P)='soil aqueous MgHPO4 content micropore non-band'
  unitl(cid_ZMG1P)='mol d-2';vartypes(cid_ZMG1P)=var_state_type

  varl(cid_Precp_AlPO4_conc)='Precp_AlPO4_conc';varlnml(cid_Precp_AlPO4_conc)='non-band soil precipitated AlPO4'
  unitl(cid_Precp_AlPO4_conc)='mol m-3';vartypes(cid_Precp_AlPO4_conc)=var_state_type

  varl(cid_Precp_CaHPO4_conc) ='Precp_CaHPO4_conc';varlnml(cid_Precp_CaHPO4_conc)='non-band soil precipitated CaHPO4'
  unitl(cid_Precp_CaHPO4_conc)='mol m-3';vartypes(cid_Precp_CaHPO4_conc)=var_state_type

  varl(cid_Precp_Ca5P3O12O3H3_conc) ='Precp_Ca5P3O12O3H3_conc';varlnml(cid_Precp_Ca5P3O12O3H3_conc)='non-band soil precipitated Ca5(PO4)3OH hydroxyapatite'
  unitl(cid_Precp_Ca5P3O12O3H3_conc)='mol m-3';vartypes(cid_Precp_Ca5P3O12O3H3_conc)=var_state_type

  varl(cid_Precp_CaH4P2O8_conc) ='Precp_CaH4P2O8_conc';varlnml(cid_Precp_CaH4P2O8_conc)='non-band soil precipitated Ca(H2PO4)2'
  unitl(cid_Precp_CaH4P2O8_conc)='mol m-3';vartypes(cid_Precp_CaH4P2O8_conc)=var_state_type

  varl(cid_Precp_FePO4_conc) ='Precp_FePO4_conc';varlnml(cid_Precp_FePO4_conc)='non-band soil precipitated FePO4'
  unitl(cid_Precp_FePO4_conc)='mol m-3';vartypes(cid_Precp_FePO4_conc)=var_state_type

  varl(cid_PALOH)='PALOH';varlnml(cid_PALOH)='precipitated Al(OH)3'
  unitl(cid_PALOH)='mol d-2';vartypes(cid_PALOH)=var_state_type

  varl(cid_PFEOH)='PFEOH';varlnml(cid_PFEOH)='precipitated Fe(OH)3'
  unitl(cid_PFEOH)='mol d-2';vartypes(cid_PFEOH)=var_state_type

  varl(cid_PCACO)='PCACO';varlnml(cid_PCACO)='precipitated CaCO3'
  unitl(cid_PCACO)='mol d-2';vartypes(cid_PCACO)=var_state_type

  varl(cid_PCASO)='PCASO';varlnml(cid_PCASO)='precipitated CaSO4'
  unitl(cid_PCASO)='mol d-2';;vartypes(cid_PCASO)=var_state_type

  varl(cid_XHY)='XHY';varlnml(cid_XHY)='exchangeable H(+)'
  unitl(cid_XHY)='mol d-2';vartypes(cid_XHY)=var_state_type

  varl(cid_XAL)='XAL';varlnml(cid_XAL)='exchangeable Al'
  unitl(cid_XAL)='mol d-2';vartypes(cid_XAL)=var_state_type

  varl(cid_XFE)='XFE';varlnml(cid_XFE)='exchangeable Fe'
  unitl(cid_XFE)='mol d-2';vartypes(cid_XFE)=var_state_type

  varl(cid_XCA)='XCA';varlnml(cid_XCA)='exchangeable Ca'
  unitl(cid_XCA)='mol d-2';vartypes(cid_XCA)=var_state_type

  varl(cid_XMG)='XMG';varlnml(cid_XMG)='exchangeable Mg'
  unitl(cid_XMG)='mol d-2';vartypes(cid_XMG)=var_state_type

  varl(cid_XNA)='XNA';varlnml(cid_XNA)='exchangeable Na'
  unitl(cid_XNA)='mol d-2';vartypes(cid_XNA)=var_state_type

  varl(cid_XKA)='XKA';varlnml(cid_XKA)='exchangeable K'
  unitl(cid_XKA)='mol d-2';vartypes(cid_XKA)=var_state_type

  varl(cid_XHC)='XHC';varlnml(cid_XHC)='exchangeable COOH'
  unitl(cid_XHC)='mol d-2';vartypes(cid_XHC)=var_state_type

  varl(cid_XROH1_conc)='XROH1_conc';varlnml(cid_XROH1_conc)='exchangeable OH  non-band'
  unitl(cid_XROH1_conc)='mol d-2';vartypes(cid_XROH1_conc)=var_state_type

  varl(cid_XALO2)='XALO2';varlnml(cid_XALO2)='exchangeable Al(OH)2'
  unitl(cid_XALO2)='mol d-2';vartypes(cid_XALO2)=var_state_type

  varl(cid_XFEO2)='XFEO2';varlnml(cid_XFEO2)='exchangeable Fe(OH)2'
  unitl(cid_XFEO2)='mol d-2';vartypes(cid_XFEO2)=var_state_type

  varl(cid_XNH4_conc)  ='XNH4_conc';varlnml(cid_XNH4_conc)='non-band soil adsorbed NH4 concentration'
  unitl(cid_XNH4_conc) ='mol m-3';vartypes(cid_XNH4_conc)=var_state_type

  varl(cid_XOH_conc)='XOH_conc';varlnml(cid_XOH_conc)='exchangeable OH- non-band'
  unitl(cid_XOH_conc)='mol d-2';vartypes(cid_XOH_conc)=var_state_type

  varl(cid_XHPO4_conc) ='XHPO4_conc';varlnml(cid_XHPO4_conc)='non-band soil concentration of adsorbed HPO4'
  unitl(cid_XHPO4_conc)='mol m-3';vartypes(cid_XHPO4_conc)=var_state_type

  varl(cid_XROH2_conc) ='XROH2_conc';varlnml(cid_XROH2_conc)='non-band soil exchangeable site R-OH2'
  unitl(cid_XROH2_conc)='mol m-3';vartypes(cid_XROH2_conc)=var_state_type

  varl(cid_XH2PO4_conc) ='XH2PO4_conc';varlnml(cid_XH2PO4_conc)='non-band soil concentration of adsorbed H2PO4'
  unitl(cid_XH2PO4_conc)='mol m-3';vartypes(cid_XH2PO4_conc)=var_state_type

  varl(cid_H1PO4_2e_band_conc)='H1PO4_2e_band_conc';varlnml(cid_H1PO4_2e_band_conc)='band soil aqueous HPO4 content micropore'
  unitl(cid_H1PO4_2e_band_conc)='mol m-3';vartypes(cid_H1PO4_2e_band_conc)=var_state_type

  varl(cid_H2PO4_1e_band_conc)='H2PO4_1e_band_conc';varlnml(cid_H2PO4_1e_band_conc)='band soil micropore aqueous H2PO4 content'
  unitl(cid_H2PO4_1e_band_conc)='mol m-3';vartypes(cid_H2PO4_1e_band_conc)=var_state_type

  varl(cid_NH3_aqu_band_conc) ='NH3_aqu_band_conc';varlnml(cid_NH3_aqu_band_conc)='band soil NH3 concentration'
  unitl(cid_NH3_aqu_band_conc)='mol m-3';vartypes(cid_NH3_aqu_band_conc)=var_state_type

  varl(cid_NH4_1p_band_conc)  ='NH4_1p_band_conc';varlnml(cid_NH4_1p_band_conc)='band soil NH4 concentration'
  unitl(cid_NH4_1p_band_conc) ='mol m-3';vartypes(cid_NH4_1p_band_conc)=var_state_type

  varl(cid_ZNO3B) ='ZNO3B';varlnml(cid_ZNO3B) ='NO3 mass band micropore'
  unitl(cid_ZNO3B)='g N d-2';vartypes(cid_ZNO3B)=var_state_type

  varl(cid_H0POB)='H0POB';varlnml(cid_H0POB)='band soil aqueous PO4 content micropore'
  unitl(cid_H0POB)='mol d-2'; vartypes(cid_H0POB)=var_state_type

  varl(cid_H3POB)='H3POB';varlnml(cid_H3POB)='band soil aqueous H3PO4 content micropore'
  unitl(cid_H3POB)='mol d-2';vartypes(cid_H3POB)=var_state_type

  varl(cid_ZFE1PB)='ZFE1PB';varlnml(cid_ZFE1PB)='band soil aqueous FeHPO4 content micropore'
  unitl(cid_ZFE1PB)='mol d-2';vartypes(cid_ZFE1PB)=var_state_type

  varl(cid_ZFE2PB)='ZFE2PB';varlnml(cid_ZFE2PB)='band soil aqueous FeH2PO4 content micropore'
  unitl(cid_ZFE2PB)='mol d-2';vartypes(cid_ZFE2PB)=var_state_type

  varl(cid_ZCA0PB) ='ZCA0PB';varlnml(cid_ZCA0PB)='band soil aqueous CaPO4 content micropore'
  unitl(cid_ZCA0PB)='mol d-2';vartypes(cid_ZCA0PB)=var_state_type

  varl(cid_ZCA1PB) ='ZCA1PB';varlnml(cid_ZCA1PB)='band soil aqueous CaHPO4 content micropore'
  unitl(cid_ZCA1PB)='mol d-2';vartypes(cid_ZCA1PB)=var_state_type

  varl(cid_ZCA2PB) ='ZCA2PB';varlnml(cid_ZCA2PB)='band soil aqueous CaH4P2O8 content micropore'
  unitl(cid_ZCA2PB)= 'mol d-2';vartypes(cid_ZCA2PB)=var_state_type

  varl(cid_ZMG1PB)='ZMG1PB';varlnml(cid_ZMG1PB)='band soil aqueous MgHPO4 content micropore'
  unitl(cid_ZMG1PB)='mol d-2';vartypes(cid_ZMG1PB)=var_state_type

  varl(cid_PrecpB_AlPO4_conc) ='PrecpB_AlPO4_conc';varlnml(cid_PrecpB_AlPO4_conc)='band soil precipitated AlPO4'
  unitl(cid_PrecpB_AlPO4_conc)='mol m-3';vartypes(cid_PrecpB_AlPO4_conc)=var_state_type

  varl(cid_PrecpB_CaHPO4_conc) ='PrecpB_CaHPO4_conc';varlnml(cid_PrecpB_CaHPO4_conc)='band soil precipitated CaHPO4'
  unitl(cid_PrecpB_CaHPO4_conc)='mol m-3';vartypes(cid_PrecpB_CaHPO4_conc)=var_state_type

  varl(cid_PrecpB_Ca5P3O12O3H3_conc) ='PrecpB_Ca5P3O12O3H3_conc';varlnml(cid_PrecpB_Ca5P3O12O3H3_conc)='band soil precipitated Ca5(PO4)3OH hydroxyapatite'
  unitl(cid_PrecpB_Ca5P3O12O3H3_conc)='mol m-3';vartypes(cid_PrecpB_Ca5P3O12O3H3_conc)=var_state_type

  varl(cid_PrecpB_CaH4P2O8_conc) ='PrecpB_CaH4P2O8_conc';varlnml(cid_PrecpB_CaH4P2O8_conc)='band soil precipitated CaH4P2O8'
  unitl(cid_PrecpB_CaH4P2O8_conc)='mol m-3';vartypes(cid_PrecpB_CaH4P2O8_conc)=var_state_type

  varl(cid_PrecpB_FePO4_con) ='PrecpB_FePO4_con';varlnml(cid_PrecpB_FePO4_con)='band soil precipitated FePO4'
  unitl(cid_PrecpB_FePO4_con)='mol m-3';vartypes(cid_PrecpB_FePO4_con)=var_state_type

  varl(cid_XNH4_band_conc)  ='XNH4_band_conc';varlnml(cid_XNH4_band_conc)='band soil adsorbed NH4 concentration'
  unitl(cid_XNH4_band_conc) ='mol m-3';vartypes(cid_XNH4_band_conc)=var_state_type

  varl(cid_XROH1_band_conc) ='XROH1_band_conc';varlnml(cid_XROH1_band_conc)='band soil exchangeable R-OH'
  unitl(cid_XROH1_band_conc)='mol d-2';vartypes(cid_XROH1_band_conc)=var_state_type

  varl(cid_XHPO4_band_conc) ='XHPO4_band_conc';varlnml(cid_XHPO4_band_conc)='band soil exchangeable HPO4 concentration'
  unitl(cid_XHPO4_band_conc)='mol m-3';vartypes(cid_XHPO4_band_conc)=var_state_type

  varl(cid_XH2PO4_band_conc) ='XH2PO4_band_conc';varlnml(cid_XH2PO4_band_conc)='band soil exchangeable H2PO4 concentration'
  unitl(cid_XH2PO4_band_conc)='mol m-3';vartypes(cid_XH2PO4_band_conc)=var_state_type

  varl(cid_XROH_band_conc) ='XROH_band_conc';varlnml(cid_XROH_band_conc)='band soil concentration of adsorbed HPO4'
  unitl(cid_XROH_band_conc)='mol m-3';vartypes(cid_XROH_band_conc)=var_state_type

  varl(cid_XROH2_band_conc) ='XROH2_band_conc';varlnml(cid_XROH2_band_conc)='band soil exchangeable site R-OH2'
  unitl(cid_XROH2_band_conc)='mol m-3';vartypes(cid_XROH2_band_conc)=var_state_type

  varl(fid_TR_CaCO3_precip_soil)='TR_CaCO3_precip_soil';varlnml(fid_TR_CaCO3_precip_soil)='total precipitated CaCO3 transformation'
  unitl(fid_TR_CaCO3_precip_soil)='mol d-2 h-1';vartypes(fid_TR_CaCO3_precip_soil)=var_flux_type

  varl(fid_TR_NaCO3_soil)='TR_NaCO3_soil';varlnml(fid_TR_NaCO3_soil)='total solute NaCO3 transformation'
  unitl(fid_TR_NaCO3_soil)='mol d-2 h-1';vartypes(fid_TR_NaCO3_soil)=var_flux_type

  varl(fid_TR_MgCO3_soil)='TR_MgCO3_soil';varlnml(fid_TR_MgCO3_soil)='total solute MgCO3 transformation'
  unitl(fid_TR_MgCO3_soil)='mol d-2 h-1'; vartypes(fid_TR_MgCO3_soil)=var_flux_type

  varl(fid_TR_CaCO3_soil)='TR_CaCO3_soil';varlnml(fid_TR_CaCO3_soil)='total solute CaCO3 transformation'
  unitl(fid_TR_CaCO3_soil)='mol d-2 h-1';vartypes(fid_TR_CaCO3_soil)=var_flux_type

  varl(fid_TR_MgHCO3_soil)='TR_MgHCO3_soil';varlnml(fid_TR_MgHCO3_soil)='total solute MgHCO3 transformation'
  unitl(fid_TR_MgHCO3_soil)='mol d-2 h-1';vartypes(fid_TR_MgHCO3_soil)=var_flux_type

  varl(fid_TR_CaHCO3_soil)='TR_CaHCO3_soil';varlnml(fid_TR_CaHCO3_soil)='total solute CaHCO3 transformation'
  unitl(fid_TR_CaHCO3_soil)='mol d-2 h-1';vartypes(fid_TR_CaHCO3_soil)=var_flux_type

  varl(fid_TR_HCO3)='TR_HCO3';varlnml(fid_TR_HCO3)='total solute HCO3 transformation'
  unitl(fid_TR_HCO3)='mol d-2 h-1';vartypes(fid_TR_HCO3)=var_flux_type

  varl(fid_TR_HCO3_sorbed_soil)='TR_HCO3_sorbed_soil';varlnml(fid_TR_HCO3_sorbed_soil)='total adsorbed COOH transformation'
  unitl(fid_TR_HCO3_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_HCO3_sorbed_soil)=var_flux_type

  varl(fid_TR_CaH4P2O8_soil)='TR_CaH4P2O8_soil';varlnml(fid_TR_CaH4P2O8_soil)='non-band soil total solute CaH4P2O8 transformation'
  unitl(fid_TR_CaH4P2O8_soil)='mol d-2 h-1';vartypes(fid_TR_CaH4P2O8_soil)=var_flux_type

  varl(fid_TR_FeH2PO4_soil)='TR_FeH2PO4_soil';varlnml(fid_TR_FeH2PO4_soil)='non-band soil total solute FeH2PO4 transformation'
  unitl(fid_TR_FeH2PO4_soil)='mol d-2 h-1';vartypes(fid_TR_FeH2PO4_soil)=var_flux_type

  varl(fid_TR_MgHPO4_soil)='TR_MgHPO4_soil';varlnml(fid_TR_MgHPO4_soil)='non-band soil total solute MgHPO4 transformation'
  unitl(fid_TR_MgHPO4_soil)='mol d-2 h-1';vartypes(fid_TR_MgHPO4_soil)=var_flux_type

  varl(fid_TR_CaHPO4_soil)='TR_CaHPO4_soil';varlnml(fid_TR_CaHPO4_soil)='non-band soil total solute CaHPO4 transformation'
  unitl(fid_TR_CaHPO4_soil)='mol d-2 h-1'; vartypes(fid_TR_CaHPO4_soil)=var_flux_type

  varl(fid_TR_FeHPO4_soil)='TR_FeHPO4_soil';varlnml(fid_TR_FeHPO4_soil)='non-band soil total solute FeHPO4 transformation'
  unitl(fid_TR_FeHPO4_soil)='mol d-2 h-1';vartypes(fid_TR_FeHPO4_soil)=var_flux_type

  varl(fid_TR_CaPO4_soil)='TR_CaPO4_soil';varlnml(fid_TR_CaPO4_soil)='non-band soil total solute CaPO4 transformation'
  unitl(fid_TR_CaPO4_soil)='mol d-2 h-1'; vartypes(fid_TR_CaPO4_soil)=var_flux_type

  varl(fid_TR_PO4_soil)='TR_PO4_soil';varlnml(fid_TR_PO4_soil)='non-band soil total solute PO4 transformation'
  unitl(fid_TR_PO4_soil)='mol d-2 h-1';vartypes(fid_TR_PO4_soil)=var_flux_type

  varl(fid_TR_FePO4_precip_soil)= 'TR_FePO4_precip_soil';varlnml(fid_TR_FePO4_precip_soil)='non-band soil total precipitated FePO4 transformation'
  unitl(fid_TR_FePO4_precip_soil)= 'mol d-2 h-1';vartypes(fid_TR_FePO4_precip_soil)=var_flux_type

  varl(fid_TR_AlPO4_precip_soil)= 'TR_AlPO4_precip_soil';varlnml(fid_TR_AlPO4_precip_soil)='non-band total precipitated AlPO4 transformation'
  unitl(fid_TR_AlPO4_precip_soil)= 'mol d-2 h-1';vartypes(fid_TR_AlPO4_precip_soil)=var_flux_type

  varl(fid_TR_CaH4P2O8_precip_soil)= 'TR_CaH4P2O8_precip_soil';varlnml(fid_TR_CaH4P2O8_precip_soil)='non-band soil total precipitated apatite transformation'
  unitl(fid_TR_CaH4P2O8_precip_soil)= 'mol d-2 h-1';vartypes(fid_TR_CaH4P2O8_precip_soil)=var_flux_type

  varl(fid_TR_CaHPO4_precip_soil)= 'TR_CaHPO4_precip_soil';varlnml(fid_TR_CaHPO4_precip_soil)='non-band soil total precipitated CaHPO4 transformation'
  unitl(fid_TR_CaHPO4_precip_soil)= 'mol d-2 h-1';vartypes(fid_TR_CaHPO4_precip_soil)=var_flux_type

  varl(fid_TR_apatite_precip_soil)= 'TR_apatite_precip_soil';varlnml(fid_TR_apatite_precip_soil)='non-band soil total precipitated CaH4P2O8 transformation'
  unitl(fid_TR_apatite_precip_soil)= 'mol d-2 h-1';vartypes(fid_TR_apatite_precip_soil)=var_flux_type

  varl(fid_TR_FeOH_soil)='TR_FeOH_soil';varlnml(fid_TR_FeOH_soil)='total solute FeOH transformation'
  unitl(fid_TR_FeOH_soil)='mol d-2 h-1';    vartypes(fid_TR_FeOH_soil)=var_flux_type

  varl(fid_TR_AlOH_soil)='TR_AlOH_soil';varlnml(fid_TR_AlOH_soil)='total solute AlOH transformation'
  unitl(fid_TR_AlOH_soil)='mol d-2 h-1'; vartypes(fid_TR_AlOH_soil)=var_flux_type

  varl(fid_TR_FeO2H2_soil)='TR_FeO2H2_soil';varlnml(fid_TR_FeO2H2_soil)='total solute FeOH2 transformation'
  unitl(fid_TR_FeO2H2_soil)='mol d-2 h-1';    vartypes(fid_TR_FeO2H2_soil)=var_flux_type

  varl(fid_TR_AlO2H2_soil)='TR_AlO2H2_soil';varlnml(fid_TR_AlO2H2_soil)='total solute AlOH2 transformation'
  unitl(fid_TR_AlO2H2_soil)='mol d-2 h-1'; vartypes(fid_TR_AlO2H2_soil)=var_flux_type

  varl(fid_TR_FeO3H3_soil)='TR_FeO3H3_soil';varlnml(fid_TR_FeO3H3_soil)='total solute FeOH3 transformation'
  unitl(fid_TR_FeO3H3_soil)='mol d-2 h-1';    vartypes(fid_TR_FeO3H3_soil)=var_flux_type

  varl(fid_TR_AlO3H3_soil)='TR_AlO3H3_soil';varlnml(fid_TR_AlO3H3_soil)='total solute AlOH3 transformation'
  unitl(fid_TR_AlO3H3_soil)='mol d-2 h-1'; vartypes(fid_TR_AlO3H3_soil)=var_flux_type

  varl(fid_TR_FeO4H4_soil)='TR_FeO4H4_soil';varlnml(fid_TR_FeO4H4_soil)='total solute FeOH4 transformation'
  unitl(fid_TR_FeO4H4_soil)='mol d-2 h-1';    vartypes(fid_TR_FeO4H4_soil)=var_flux_type

  varl(fid_TR_AlO4H4_soil)='TR_AlO4H4_soil';varlnml(fid_TR_AlO4H4_soil)='total solute AlOH4 transformation'
  unitl(fid_TR_AlO4H4_soil)='mol d-2 h-1'; vartypes(fid_TR_AlO4H4_soil)=var_flux_type


  varl(fid_TR_MgOH_soil)='TR_MgOH_soil';varlnml(fid_TR_MgOH_soil)='total solute MgOH transformation'
  unitl(fid_TR_MgOH_soil)='mol d-2 h-1';vartypes(fid_TR_MgOH_soil)=var_flux_type

  varl(fid_TR_CaOH_soil)='TR_CaOH_soil';varlnml(fid_TR_CaOH_soil)='total solute CaOH transformation'
  unitl(fid_TR_CaOH_soil)='mol d-2 h-1';vartypes(fid_TR_CaOH_soil)=var_flux_type

  varl(fid_TR_FeOH3_precip_soil) ='TR_FeOH3_precip_soil';varlnml(fid_TR_FeOH3_precip_soil)='total precipitated FeOH3 transformation'
  unitl(fid_TR_FeOH3_precip_soil)='mol d-2 h-1'; vartypes(fid_TR_FeOH3_precip_soil)=var_flux_type

  varl(fid_TR_AlOH3_precip_soil)='TR_AlOH3_precip_soil';varlnml(fid_TR_AlOH3_precip_soil)='total precipitated AlOH3 transformation'
  unitl(fid_TR_AlOH3_precip_soil)='mol d-2 h-1'; vartypes(fid_TR_AlOH3_precip_soil)=var_flux_type

  varl(fid_TR_FeO2H2_sorbed_soil)='TR_FeO2H2_sorbed_soil';varlnml(fid_TR_FeO2H2_sorbed_soil)='total FeOH2 adsorption'
  unitl(fid_TR_FeO2H2_sorbed_soil)='mol d-2 h-1'; vartypes(fid_TR_FeO2H2_sorbed_soil)=var_flux_type

  varl(fid_TR_Al_3p_soil)  = 'TR_Al_3p_soil';varlnml(fid_TR_Al_3p_soil)='non-band soil total solute Al transformation'
  unitl(fid_TR_Al_3p_soil) = 'mol d-2 h-1';vartypes(fid_TR_Al_3p_soil)=var_flux_type

  varl(fid_TR_AlSO4_soil)='TR_AlSO4_soil'; varlnml(fid_TR_AlSO4_soil)='total solute AlSO4 transformation'
  unitl(fid_TR_AlSO4_soil)='mol d-2 h-1'; vartypes(fid_TR_AlSO4_soil)=var_flux_type

  varl(fid_TR_RHPO4_sorbed_band_soil)='TR_RHPO4_sorbed_band_soil';varlnml(fid_TR_RHPO4_sorbed_band_soil)='band soil total adsorbed HPO4 transformation'
  unitl(fid_TR_RHPO4_sorbed_band_soil)='mol d-2 h-1';vartypes(fid_TR_RHPO4_sorbed_band_soil)=var_flux_type

  varl(fid_TR_RH2PO4_sorbed_band_soil)='TR_RH2PO4_sorbed_band_soil';varlnml(fid_TR_RH2PO4_sorbed_band_soil)='band soil total adsorbed H2PO4 transformation'
  unitl(fid_TR_RH2PO4_sorbed_band_soil)='mol d-2 h-1';vartypes(fid_TR_RH2PO4_sorbed_band_soil)=var_flux_type

  varl(fid_TR_RO_sorbed_band_soil)='TR_RO_sorbed_band_soil';varlnml(fid_TR_RO_sorbed_band_soil)='band soil total adsorbed OH- transformation'
  unitl(fid_TR_RO_sorbed_band_soil)='mol d-2 h-1'; vartypes(fid_TR_RO_sorbed_band_soil)=var_flux_type

  varl(fid_TR_ROH_sorbed_band_soil)='TR_ROH_sorbed_band_soil';varlnml(fid_TR_ROH_sorbed_band_soil)='band soil total adsorbed OH transformation'
  unitl(fid_TR_ROH_sorbed_band_soil)='mol d-2 h-1'; vartypes(fid_TR_ROH_sorbed_band_soil)=var_flux_type

  varl(fid_TR_ROH2_sorbed_band_soil)='TR_ROH2_sorbed_band_soil';varlnml(fid_TR_ROH2_sorbed_band_soil)='band soil total adsorbed OH2 transformation'
  unitl(fid_TR_ROH2_sorbed_band_soil)='mol d-2 h-1'; vartypes(fid_TR_ROH2_sorbed_band_soil)=var_flux_type

  varl(fid_TR_Ca_2p_soil) ='TR_Ca_2p_soil';varlnml(fid_TR_Ca_2p_soil)='total solute Ca transformation'
  unitl(fid_TR_Ca_2p_soil)='mol d-2 h-1';vartypes(fid_TR_Ca_2p_soil)=var_flux_type

  varl(fid_TR_CaSO4_soil)='TR_CaSO4_soil';varlnml(fid_TR_CaSO4_soil)='total solute CaSO4 transformation'
  unitl(fid_TR_CaSO4_soil)='mol d-2 h-1';vartypes(fid_TR_CaSO4_soil)=var_flux_type

  varl(fid_TR_CaSO4_precip_soil)='TR_CaSO4_precip_soil';varlnml(fid_TR_CaSO4_precip_soil)='total precipitated CaSO4 transformation'
  unitl(fid_TR_CaSO4_precip_soil)='mol d-2 h-1';vartypes(fid_TR_CaSO4_precip_soil)=var_flux_type

  varl(fid_TR_CO2_gchem_soil_vr)='TR_CO2_gchem_soil_vr';varlnml(fid_TR_CO2_gchem_soil_vr)='total solute CO2 transformation'
  unitl(fid_TR_CO2_gchem_soil_vr)='mol d-2 h-1';vartypes(fid_TR_CO2_gchem_soil_vr)=var_flux_type

  varl(fid_TR_CO3_2e_soil)='TR_CO3_2e_soil';varlnml(fid_TR_CO3_2e_soil)='total solute CO3 transformation'
  unitl(fid_TR_CO3_2e_soil)='mol d-2 h-1';vartypes(fid_TR_CO3_2e_soil)=var_flux_type

  varl(fid_TR_Fe_3p_soil) ='TR_Fe_3p_soil';varlnml(fid_TR_Fe_3p_soil)='total solute Fe transformation'
  unitl(fid_TR_Fe_3p_soil)='mol d-2 h-1';vartypes(fid_TR_Fe_3p_soil)= var_flux_type

  varl(fid_TR_FeSO4_soil)='TR_FeSO4_soil';varlnml(fid_TR_FeSO4_soil)='total solute FeSO4 transformation'
  unitl(fid_TR_FeSO4_soil)='mol d-2 h-1';vartypes(fid_TR_FeSO4_soil)=var_flux_type

  varl(fid_TR_H1PO4_soil) = 'TR_H1PO4_soil';varlnml(fid_TR_H1PO4_soil)='non-band soil total solute HPO4 transformation'
  unitl(fid_TR_H1PO4_soil) = 'mol d-2 h-1';vartypes(fid_TR_H1PO4_soil)=var_flux_type

  varl(fid_TR_H2PO4_soil) = 'TR_H2PO4_soil';varlnml(fid_TR_H2PO4_soil)='non-band soil total solute H2PO4 transformation'
  unitl(fid_TR_H2PO4_soil)= 'mol d-2 h-1';vartypes(fid_TR_H2PO4_soil)=var_flux_type

  varl(fid_TR_H3PO4_sorbed_soil) = 'TR_H3PO4_sorbed_soil';varlnml(fid_TR_H3PO4_sorbed_soil)='non-band soil total solute H3PO4 transformation'
  unitl(fid_TR_H3PO4_sorbed_soil)= 'mol d-2 h-1';vartypes(fid_TR_H3PO4_sorbed_soil)=var_flux_type

  varl(fid_TR_H_p_soil)='TR_H_p_soil';varlnml(fid_TR_H_p_soil)='total solute H transformation'
  unitl(fid_TR_H_p_soil)='mol d-2 h-1'; vartypes(fid_TR_H_p_soil)=var_flux_type

  varl(fid_TR_K_1p_soil)='TR_K_1p_soil';varlnml(fid_TR_K_1p_soil)='total solute K transformation'
  unitl(fid_TR_K_1p_soil)='mol d-2 h-1'; vartypes(fid_TR_K_1p_soil)=var_flux_type

  varl(fid_TR_KSO4_soil)='TR_KSO4_soil';varlnml(fid_TR_KSO4_soil)='total solute KSO4 transformation'
  unitl(fid_TR_KSO4_soil)='mol d-2 h-1'; vartypes(fid_TR_KSO4_soil)=var_flux_type

  varl(fid_TR_Mg_2p_soil)='TR_Mg_2p_soil';varlnml(fid_TR_Mg_2p_soil)='total solute Mg transformation'
  unitl(fid_TR_Mg_2p_soil)='mol d-2 h-1'; vartypes(fid_TR_Mg_2p_soil)=var_flux_type

  varl(fid_TR_MgSO4_soil)='TR_MgSO4_soil';varlnml(fid_TR_MgSO4_soil)='total solute MgSO4 transformation'
  unitl(fid_TR_MgSO4_soil)='mol d-2 h-1';vartypes(fid_TR_MgSO4_soil)=var_flux_type

  varl(fid_TR_NH3_soil_vr) = 'TR_NH3_soil_vr';varlnml(fid_TR_NH3_soil_vr)='non-band total solute NH3 transformation'
  unitl(fid_TR_NH3_soil_vr)= 'mol d-2 h-1';vartypes(fid_TR_NH3_soil_vr)=var_flux_type

  varl(fid_TR_NH4_soil) = 'TR_NH4_soil';varlnml(fid_TR_NH4_soil)='non-band soil total solute NH4 transformation'
  unitl(fid_TR_NH4_soil)= 'mol d-2 h-1';vartypes(fid_TR_NH4_soil)=var_flux_type

  varl(fid_TR_Na_p_soil) ='TR_Na_p_soil';varlnml(fid_TR_Na_p_soil)='total solute Na transformation'
  unitl(fid_TR_Na_p_soil)='mol d-2 h-1';vartypes(fid_TR_Na_p_soil)=var_flux_type

  varl(fid_TR_NaSO4_soil)='TR_NaSO4_soil';varlnml(fid_TR_NaSO4_soil)='total solute NaSO4 transformation'
  unitl(fid_TR_NaSO4_soil)='mol d-2 h-1';vartypes(fid_TR_NaSO4_soil)=var_flux_type

  varl(fid_TR_OH_1e_soil)='TR_OH_1e_soil';varlnml(fid_TR_OH_1e_soil)='total solute OH transformation'
  unitl(fid_TR_OH_1e_soil)='mol d-2 h-1';vartypes(fid_TR_OH_1e_soil)=var_flux_type

  varl(fid_TR_SO4_2e_soil)='TR_SO4_2e_soil';varlnml(fid_TR_SO4_2e_soil)='total solute SO4 transformation'
  unitl(fid_TR_SO4_2e_soil)='mol d-2 h-1';vartypes(fid_TR_SO4_2e_soil)=var_flux_type

  varl(fid_TR_RHPO4_sorbed_soil)='TR_RHPO4_sorbed_soil';varlnml(fid_TR_RHPO4_sorbed_soil)='non-band soil total adsorbed HPO4 transformation'
  unitl(fid_TR_RHPO4_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_RHPO4_sorbed_soil)=var_flux_type

  varl(fid_TR_RH2PO4_sorbed_soil)='TR_RH2PO4_sorbed_soil';varlnml(fid_TR_RH2PO4_sorbed_soil)='non-band soil total adsorbed H2PO4 transformation'
  unitl(fid_TR_RH2PO4_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_RH2PO4_sorbed_soil)=var_flux_type

  varl(fid_TR_Al_sorbed_soil)='TR_Al_sorbed_soil';varlnml(fid_TR_Al_sorbed_soil)='total adsorbed Al transformation'
  unitl(fid_TR_Al_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_Al_sorbed_soil)=var_flux_type

  varl(fid_TR_AlO2H2_sorbed_soil)='TR_AlO2H2_sorbed_soil';varlnml(fid_TR_AlO2H2_sorbed_soil)='total adsorbed AlOH2 transformation'
  unitl(fid_TR_AlO2H2_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_AlO2H2_sorbed_soil)=var_flux_type

  varl(fid_TR_Ca_sorbed_soil)='TR_Ca_sorbed_soil';varlnml(fid_TR_Ca_sorbed_soil)='total adsorbed Ca transformation'
  unitl(fid_TR_Ca_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_Ca_sorbed_soil)=var_flux_type

  varl(fid_TR_Fe_sorbed_soil)='TR_Fe_sorbed_soil';varlnml(fid_TR_Fe_sorbed_soil)='total Fe adsorption'
  unitl(fid_TR_Fe_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_Fe_sorbed_soil)=var_flux_type

  varl(fid_TR_RO_sorbed_soil)='TR_RO_sorbed_soil';varlnml(fid_TR_RO_sorbed_soil)='non-band soil total adsorbed OH- transformation'
  unitl(fid_TR_RO_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_RO_sorbed_soil)=var_flux_type

  varl(fid_TR_ROH_sorbed_soil)='TR_ROH_sorbed_soil';varlnml(fid_TR_ROH_sorbed_soil)='non-band total adsorbed OH transformation'
  unitl(fid_TR_ROH_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_ROH_sorbed_soil)=var_flux_type

  varl(fid_TR_ROH2_sorbed_soil)='TR_ROH2_sorbed_soil';varlnml(fid_TR_ROH2_sorbed_soil)='non-band soil total adsorbed OH2 transformation'
  unitl(fid_TR_ROH2_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_ROH2_sorbed_soil)=var_flux_type

  varl(fid_TR_H_p_sorbed_soil)='TR_H_p_sorbed_soil';varlnml(fid_TR_H_p_sorbed_soil)='total adsorbed H transformation'
  unitl(fid_TR_H_p_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_H_p_sorbed_soil)=var_flux_type

  varl(fid_TR_K_sorbed_soil)='TR_K_sorbed_soil';varlnml(fid_TR_K_sorbed_soil)='total adsorbed K transformation'
  unitl(fid_TR_K_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_K_sorbed_soil)=var_flux_type

  varl(fid_TR_Mg_sorbed_soil)='TR_Mg_sorbed_soil';varlnml(fid_TR_Mg_sorbed_soil)='total adsorbed Mg transformation'
  unitl(fid_TR_Mg_sorbed_soil)='mol d-2 h-1'; vartypes(fid_TR_Mg_sorbed_soil)=var_flux_type

  varl(fid_TR_NH4_sorbed_soil) = 'TR_NH4_sorbed_soil';varlnml(fid_TR_NH4_sorbed_soil)='non-band soil total adsorbed NH4 transformation'
  unitl(fid_TR_NH4_sorbed_soil)= 'mol d-2 h-1';vartypes(fid_TR_NH4_sorbed_soil)=var_flux_type

  varl(fid_TR_Na_sorbed_soil)='TR_Na_sorbed_soil';varlnml(fid_TR_Na_sorbed_soil)='total adsorbed Na transformation'
  unitl(fid_TR_Na_sorbed_soil)='mol d-2 h-1';vartypes(fid_TR_Na_sorbed_soil)=var_flux_type

  varl(fid_Txchem_CO2)='Txchem_CO2';varlnml(fid_Txchem_CO2)='total solute CO2 transformation'
  unitl(fid_Txchem_CO2)='mol d-2 h-1';vartypes(fid_Txchem_CO2)=var_flux_type

  varl(fid_TBION)='TBION';varlnml(fid_TBION)='total solute ion transformation'
  unitl(fid_TBION)='mol d-2 h-1';vartypes(fid_TBION)=var_flux_type

  varl(fid_TRH2O)='TRH2O';varlnml(fid_TRH2O)='total solute H2O transformation'
  unitl(fid_TRH2O)='mol d-2 h-1';vartypes(fid_TRH2O)=var_flux_type

  varl(fid_TR_CaH4P2O8_band_soil)='TR_CaH4P2O8_band_soil';varlnml(fid_TR_CaH4P2O8_band_soil)='band soil total solute CaH4P2O8 transformation'
  unitl(fid_TR_CaH4P2O8_band_soil)='mol d-2 h-1';vartypes(fid_TR_CaH4P2O8_band_soil)=var_flux_type

  varl(fid_TR_FeH2PO4_band_soil)='TR_FeH2PO4_band_soil';varlnml(fid_TR_FeH2PO4_band_soil)='band soil total solute FeH2PO4 transformation'
  unitl(fid_TR_FeH2PO4_band_soil)='mol d-2 h-1';vartypes(fid_TR_FeH2PO4_band_soil)=var_flux_type

  varl(fid_TR_MgHPO4_band_soil)='TR_MgHPO4_band_soil';varlnml(fid_TR_MgHPO4_band_soil)='band soil total solute MgHPO4 transformation'
  unitl(fid_TR_MgHPO4_band_soil)='mol d-2 h-1';vartypes(fid_TR_MgHPO4_band_soil)=var_flux_type

  varl(fid_TR_CaHPO4_band_soil)='TR_CaHPO4_band_soil';varlnml(fid_TR_CaHPO4_band_soil)='band soil total solute CaHPO4 transformation'
  unitl(fid_TR_CaHPO4_band_soil)='mol d-2 h-1';vartypes(fid_TR_CaHPO4_band_soil)=var_flux_type

  varl(fid_TR_FeHPO4_band_soil)='TR_FeHPO4_band_soil';varlnml(fid_TR_FeHPO4_band_soil)='band soil total solute FeHPO4 transformation'
  unitl(fid_TR_FeHPO4_band_soil)='mol d-2 h-1';vartypes(fid_TR_FeHPO4_band_soil)=var_flux_type

  varl(fid_TR_CaPO4_band_soil)='TR_CaPO4_band_soil';varlnml(fid_TR_CaPO4_band_soil)='band soil total solute CaPO4 transformation'
  unitl(fid_TR_CaPO4_band_soil)='mol d-2 h-1';vartypes(fid_TR_CaPO4_band_soil)=var_flux_type

  varl(fid_TR_PO4_band_soil)='TR_PO4_band_soil';varlnml(fid_TR_PO4_band_soil)='band soil total solute PO4 transformation'
  unitl(fid_TR_PO4_band_soil)='mol d-2 h-1';vartypes(fid_TR_PO4_band_soil)=var_flux_type

  varl(fid_TR_FePO4_precip_band_soil)= 'TR_FePO4_precip_band_soil';varlnml(fid_TR_FePO4_precip_band_soil)='band soil total precipitated FePO4 transformation'
  unitl(fid_TR_FePO4_precip_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_FePO4_precip_band_soil)=var_flux_type

  varl(fid_TR_AlPO4_precip_band_soil)= 'TR_AlPO4_precip_band_soil';varlnml(fid_TR_AlPO4_precip_band_soil)='band soil total precipitated AlPO4 transformation'
  unitl(fid_TR_AlPO4_precip_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_AlPO4_precip_band_soil)=var_flux_type

  varl(fid_TR_CaH4P2O8_precip_band_soil)= 'TR_CaH4P2O8_precip_band_soil';varlnml(fid_TR_CaH4P2O8_precip_band_soil)='band soil total precipitated apatite transformation'
  unitl(fid_TR_CaH4P2O8_precip_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_CaH4P2O8_precip_band_soil)=var_flux_type

  varl(fid_TR_CaHPO4_precip_band_soil)= 'TR_CaHPO4_precip_band_soil';varlnml(fid_TR_CaHPO4_precip_band_soil)='band soil total precipitated CaHPO4 transformation'
  unitl(fid_TR_CaHPO4_precip_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_CaHPO4_precip_band_soil)=var_flux_type

  varl(fid_TR_apatite_precip_band_soil)= 'TR_apatite_precip_band_soil';varlnml(fid_TR_apatite_precip_band_soil)='band soil total precipitated CaH4P2O8 transformation'
  unitl(fid_TR_apatite_precip_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_apatite_precip_band_soil)=var_flux_type

  varl(fid_TR_H1PO4_band_soil) = 'TR_H1PO4_band_soil';varlnml(fid_TR_H1PO4_band_soil)='band soil total solute HPO4 transformation'
  unitl(fid_TR_H1PO4_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_H1PO4_band_soil)=var_flux_type

  varl(fid_TR_H2PO4_band_soil) = 'TR_H2PO4_band_soil';varlnml(fid_TR_H2PO4_band_soil)='band soil total solute H2PO4 transformation'
  unitl(fid_TR_H2PO4_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_H2PO4_band_soil)=var_flux_type

  varl(fid_TR_H3PO4_band_soil) = 'TR_H3PO4_band_soil';varlnml(fid_TR_H3PO4_band_soil)='band soil total solute H3PO4 transformation'
  unitl(fid_TR_H3PO4_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_H3PO4_band_soil)=var_flux_type

  varl(fid_TR_NH3_band_soil) = 'TR_NH3_band_soil';varlnml(fid_TR_NH3_band_soil)='band soil total solute NH3 transformation'
  unitl(fid_TR_NH3_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_NH3_band_soil)=var_flux_type

  varl(fid_TR_NH4_band_soil) = 'TR_NH4_band_soil';varlnml(fid_TR_NH4_band_soil)='band soil total solute NH4 transformation'
  unitl(fid_TR_NH4_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_NH4_band_soil)=var_flux_type

  varl(fid_TR_NH4_sorbed_band_soil) = 'TR_NH4_sorbed_band_soil';varlnml(fid_TR_NH4_sorbed_band_soil)='band soil total adsorbed NH4 transformation'
  unitl(fid_TR_NH4_sorbed_band_soil)= 'mol d-2 h-1';vartypes(fid_TR_NH4_sorbed_band_soil)=var_flux_type
  end subroutine getvarlist_salt

! ----------------------------------------------------------------------
  subroutine SetChemVar(nvars, ystates0l, chemvar)
!
  use SoluteChemDataType, only : chem_var_type
  implicit none
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  type(chem_var_type), intent(out)   :: chemvar

! Gapon selectivity coefficient
  real(r8), parameter :: GKCA=0.25_r8   !for Ca-Al
  real(r8), parameter :: GKCH=0.25_r8   !for Ca-H
  real(r8), parameter :: GKCM=0.60_r8   !for Ca-Mg
  real(r8), parameter :: GKCK=3.0_r8    !for Ca-K
  real(r8), parameter :: GKCN=0.16_r8   !for Ca-Na
  real(r8), parameter :: GKC4=2.5e-2_r8 !for Ca-NH4

  chemvar%GKC4 = GKC4
  chemvar%GKCA = GKCA
  chemvar%GKCH = GKCH
  chemvar%GKCM = GKCM
  chemvar%GKCK = GKCK
  chemvar%GKCN = GKCN

  end subroutine SetChemVar
! ----------------------------------------------------------------------

  subroutine RetrieveYstatef(nvars, ystates0l, ystatesfl, chemvar, solflx)
!
  use SoluteChemDataType, only : chem_var_type, solute_flx_type
  implicit none
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  type(chem_var_type)  , intent(in) :: chemvar
  type(solute_flx_type), intent(in) :: solflx
  real(r8), intent(out) :: ystatesfl(nvars)

  !ZNH4S=ZNH4S+TR_NH4_soil*Natomw
  ystatesfl(cid_NH4_1p_conc)=ystates0l(cid_NH4_1p_conc)+solflx%TR_NH4_soil/chemvar%VLWatMicPNH
  ystatesfl(fid_TR_NH4_soil)=solflx%TR_NH4_soil

  !ZNH3S=ZNH3S+TR_NH3_soil_vr*Natomw
  ystatesfl(cid_NH3_aqu_conc)=ystates0l(cid_NH3_aqu_conc)+solflx%TR_NH3_soil_vr/chemvar%VLWatMicPNH
  ystatesfl(fid_TR_NH3_soil_vr)=solflx%TR_NH3_soil_vr

  !XN4  =XN4+TR_NH4_sorbed_soil
  ystatesfl(cid_XNH4_conc)=ystates0l(cid_XNH4_conc)+solflx%TR_NH4_sorbed_soil/chemvar%VLWatMicPNH
  ystatesfl(fid_TR_NH4_sorbed_soil)=solflx%TR_NH4_sorbed_soil

  !ZNH4B=ZNH4B+TR_NH4_band_soil*Natomw
  ystatesfl(cid_NH4_1p_band_conc)=ystates0l(cid_NH4_1p_band_conc)+solflx%TR_NH3_band_soil/chemvar%VLWatMicPNB
  ystatesfl(fid_TR_NH3_band_soil)=solflx%TR_NH3_band_soil

  !ZNH3B=ZNH3B+TR_NH3_band_soil*Natomw
  ystatesfl(cid_NH3_aqu_band_conc)=ystates0l(cid_NH3_aqu_band_conc)+solflx%TR_NH3_band_soil/chemvar%VLWatMicPNB
  ystatesfl(fid_TR_NH3_band_soil)=solflx%TR_NH3_band_soil

  !XNB  = XNB+TR_NH4_sorbed_band_soil
  ystatesfl(cid_XNH4_band_conc)=ystates0l(cid_XNH4_band_conc)+solflx%TR_NH4_sorbed_band_soil/chemvar%VLWatMicPNB
  ystatesfl(fid_TR_NH4_sorbed_band_soil)=solflx%TR_NH4_sorbed_band_soil

  !H1PO4=H1PO4+TR_H1PO4_soil*Patomw
  ystatesfl(cid_H1PO4_2e_conc)=ystates0l(cid_H1PO4_2e_conc)+solflx%TR_H1PO4_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_H1PO4_soil)=solflx%TR_H1PO4_soil

  !H2PO4=H2PO4+TR_H2PO4_soil*Patomw
  ystatesfl(cid_H2PO4_1e_conc)=ystates0l(cid_H2PO4_1e_conc)+solflx%TR_H2PO4_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_H2PO4_soil)=solflx%TR_H2PO4_soil

  !XOH1 =XOH1+TR_ROH_sorbed_soil
  ystatesfl(cid_XROH1_conc)=ystates0l(cid_XROH1_conc)+solflx%TR_ROH_sorbed_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_ROH_sorbed_soil)=solflx%TR_ROH_sorbed_soil

  !XOH2 =XOH2+TR_ROH2_sorbed_soil
  ystatesfl(cid_XROH2_conc)=ystates0l(cid_XROH2_conc)+solflx%TR_ROH2_sorbed_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_ROH2_sorbed_soil)=solflx%TR_ROH2_sorbed_soil

  !XH1P =XH1P+TR_RHPO4_sorbed_soil
  ystatesfl(cid_XHPO4_conc)=ystates0l(cid_XHPO4_conc)+solflx%TR_RHPO4_sorbed_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_RHPO4_sorbed_soil)=solflx%TR_RHPO4_sorbed_soil

  !XH2P =XH2P+TR_RH2PO4_sorbed_soil
  ystatesfl(cid_XH2PO4_conc)=ystates0l(cid_XH2PO4_conc)+solflx%TR_RH2PO4_sorbed_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_RH2PO4_sorbed_soil)=solflx%TR_RH2PO4_sorbed_soil

  !PALPO=PALPO+TR_AlPO4_precip_soil
  ystatesfl(cid_Precp_AlPO4_conc)=ystates0l(cid_Precp_AlPO4_conc)+solflx%TR_AlPO4_precip_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_AlPO4_precip_soil)=solflx%TR_AlPO4_precip_soil

  !PFEPO=PFEPO+TR_FePO4_precip_soil
  ystatesfl(cid_Precp_FePO4_conc)=ystates0l(cid_Precp_FePO4_conc)+solflx%TR_FePO4_precip_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_FePO4_precip_soil)=solflx%TR_FePO4_precip_soil

  !PCAPD=PCAPD+TR_CaHPO4_precip_soil
  ystatesfl(cid_Precp_CaHPO4_conc)=ystates0l(cid_Precp_CaHPO4_conc)+solflx%TR_CaHPO4_precip_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_CaHPO4_precip_soil)=solflx%TR_CaHPO4_precip_soil

  !PCAPH=PCAPH+TR_apatite_precip_soil
  ystatesfl(cid_Precp_Ca5P3O12O3H3_conc)=ystates0l(cid_Precp_Ca5P3O12O3H3_conc)+solflx%TR_apatite_precip_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_apatite_precip_soil)=solflx%TR_apatite_precip_soil

  !PCAPM=PCAPM+TR_CaH4P2O8_precip_soil
  ystatesfl(cid_Precp_CaH4P2O8_conc)=ystates0l(cid_Precp_CaH4P2O8_conc)+solflx%TR_CaH4P2O8_precip_soil/chemvar%VLWatMicPPO
  ystatesfl(fid_TR_CaH4P2O8_precip_soil)=solflx%TR_CaH4P2O8_precip_soil

  !H1POB=H1POB+TR_H1PO4_band_soil*Patomw
  ystatesfl(cid_H1PO4_2e_band_conc)=ystates0l(cid_H1PO4_2e_band_conc)+solflx%TR_H1PO4_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_H1PO4_band_soil)=solflx%TR_H1PO4_band_soil

  !H2POB=H2POB+TR_H2PO4_band_soil*Patomw
  ystatesfl(cid_H2PO4_1e_band_conc)=ystates0l(cid_H2PO4_1e_band_conc)+solflx%TR_H2PO4_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_H2PO4_band_soil)=solflx%TR_H2PO4_band_soil

  !XOH1B=XOH1B+TR_ROH_sorbed_band_soil
  ystatesfl(cid_XROH_band_conc)=ystates0l(cid_XROH_band_conc)+solflx%TR_ROH_sorbed_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_ROH_sorbed_band_soil)=solflx%TR_ROH_sorbed_band_soil

  !XOH2B=XOH2B+TR_ROH2_sorbed_band_soil
  ystatesfl(cid_XROH2_band_conc)=ystates0l(cid_XROH2_band_conc)+solflx%TR_ROH2_sorbed_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_ROH2_sorbed_band_soil)=solflx%TR_ROH2_sorbed_band_soil

  !XHPO4_band_conc=XHPO4_band_conc+TR_RHPO4_sorbed_band_soil
  ystatesfl(cid_XHPO4_band_conc)=ystates0l(cid_XHPO4_band_conc)+solflx%TR_RHPO4_sorbed_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_RHPO4_sorbed_band_soil)=solflx%TR_RHPO4_sorbed_band_soil

  !XH2PB=XH2PB+TR_RH2PO4_sorbed_band_soil
  ystatesfl(cid_XH2PO4_band_conc)=ystates0l(cid_XH2PO4_band_conc)+solflx%TR_RH2PO4_sorbed_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_RH2PO4_sorbed_band_soil)=solflx%TR_RH2PO4_sorbed_band_soil

  !PALPB=PALPB+TR_AlPO4_precip_band_soil
  ystatesfl(cid_PrecpB_AlPO4_conc)=ystates0l(cid_PrecpB_AlPO4_conc)+solflx%TR_AlPO4_precip_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_AlPO4_precip_band_soil)=solflx%TR_AlPO4_precip_band_soil

  !PFEPB=PFEPB+TR_FePO4_precip_band_soil
  ystatesfl(cid_PrecpB_FePO4_con)=ystates0l(cid_PrecpB_FePO4_con)+solflx%TR_FePO4_precip_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_FePO4_precip_band_soil)=solflx%TR_FePO4_precip_band_soil

  !PCPDB=PCPDB+TR_CaHPO4_precip_band_soil
  ystatesfl(cid_PrecpB_CaHPO4_conc)=ystates0l(cid_PrecpB_CaHPO4_conc)+solflx%TR_CaHPO4_precip_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_CaHPO4_precip_band_soil)=solflx%TR_CaHPO4_precip_band_soil

  !PCPHB=PCPHB+TR_apatite_precip_band_soil
  ystatesfl(cid_PrecpB_Ca5P3O12O3H3_conc)=ystates0l(cid_PrecpB_Ca5P3O12O3H3_conc)+solflx%TR_apatite_precip_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_apatite_precip_band_soil)=solflx%TR_apatite_precip_band_soil

  !PCPMB=PCPMB+TR_CaH4P2O8_precip_band_soil
  ystatesfl(cid_PrecpB_CaH4P2O8_conc)=ystates0l(cid_PrecpB_CaH4P2O8_conc)+solflx%TR_CaH4P2O8_precip_band_soil/chemvar%VLWatMicPPB
  ystatesfl(fid_TR_CaH4P2O8_precip_band_soil)=solflx%TR_CaH4P2O8_precip_band_soil
  end subroutine RetrieveYstatef

end module AquaSaltChemMod
