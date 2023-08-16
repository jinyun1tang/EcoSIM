module ChemIDMod

implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  integer :: cid_CO2S     !aqueous CO2  micropore	[g d-2]
  integer :: cid_CH1P1    !soil aqueous HPO4 content micropore non-band, [mol m-3]
  integer :: cid_CH2P1    !soil aqueous H2PO4 content micropore non-band, [mol m-3]
  integer :: cid_CN31     !soil NH3 concentration in non-band soil, [mol m-3]
  integer :: cid_CN41     !soil NH4 concentration in non-band soil, [mol m-3]
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
  integer :: cid_ZCA2P    !soil aqueous CaH2PO4 content micropore non-band, [mol d-2]
  integer :: cid_ZMG1P    !soil aqueous MgHPO4 content micropore non-band, [mol d-2]
  integer :: cid_PALPO1   !precipitated AlPO4 non-band, [mol m-3]
  integer :: cid_PCAPD1   !precipitated CaHPO4 non-band soil, [mol m-3]
  integer :: cid_PCAPH1   !precipitated Ca5(PO4)3OH hydroxyapatite non-band soil, [mol m-3]
  integer :: cid_PCAPM1   !precipitated Ca(H2PO4)2 non-band soil, [mol m-3]
  integer :: cid_PFEPO1   !precipitated FePO4 non-band soil, [mol m-3]
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
  integer :: cid_XOH11    !exchangeable OH  non-band, [mol d-2]
  integer :: cid_XALO2    !exchangeable AlOH2 , [mol d-2]
  integer :: cid_XFEO2    !exchangeable Fe(OH)2, [mol d-2]
  integer :: cid_XN41     !exchangeable NH4 non-band soil, [mol d-2]
  integer :: cid_XOH01    !exchangeable OH- non-band, [mol d-2]
  integer :: cid_XH1P1    !exchangeable HPO4  non-band, [mol m-3]
  integer :: cid_XOH21    !exchangeable OH2  non-band soil, [mol m-3]
  integer :: cid_XH2P1    !exchangeable H2PO4  non-band soil, [mol m-3]

  integer :: cid_CN4B     !soil NH4 concentration in band soil, [mol m-3]
  integer :: cid_CN3B     !soil NH3 concentration in band soil, [mol m-3]
  integer :: cid_ZNO3B    !NO3 mass band micropore, [g d-2]
  integer :: cid_CH1PB    !soil aqueous HPO4 content micropore band, [mol m-3]
  integer :: cid_CH2PB    !soil aqueous H2PO4 content micropore  band, [mol m-3]
  integer :: cid_H0POB    !soil aqueous PO4 content micropore band, [mol d-2]
  integer :: cid_H3POB    !soil aqueous H3PO4 content micropore band, [mol d-2]
  integer :: cid_ZFE1PB   !soil aqueous FeHPO4 content micropore band, [mol d-2]
  integer :: cid_ZFE2PB   !soil aqueous FeH2PO4 content micropore band, [mol d-2
  integer :: cid_ZCA0PB   !soil aqueous CaPO4 content micropore band, [mol d-2]
  integer :: cid_ZCA1PB   !soil aqueous CaHPO4 content micropore band, [mol d-2]
  integer :: cid_ZCA2PB   !soil aqueous CaH2PO4 content micropore band, [mol d-2]
  integer :: cid_ZMG1PB   !soil aqueous MgHPO4 content micropore band, [mol d-2]
  integer :: cid_PALPOB   !precipitated AlPO4 band soil, [mol m-3]
  integer :: cid_PCAPDB   !precipitated CaHPO4 band soil, [mol m-3]
  integer :: cid_PCAPHB   !precipitated Ca5(PO4)3OH hydroxyapatite band soil, [mol m-3]
  integer :: cid_PCAPMB   !precipitated CaH2PO4 band soil, [mol m-3]
  integer :: cid_PFEPOB   !precipitated FePO4 band soil, [mol m-3]
  integer :: cid_XN4B     !exchangeable NH4 band soil, [mol d-2]
  integer :: cid_XH01B    !exchangeable OH- band, [mol d-2]
  integer :: cid_X1P1B    !exchangeable HPO4 concentration band-soil, [mol m-3]
  integer :: cid_X2P1B    !exchangeable H2PO4 concentration band-soil, [mol m-3]
  integer :: cid_XH11B    !exchangeable OH band-soil, [mol m-3]
  integer :: cid_XH21B    !exchangeable OH2 band-soil, [mol m-3]

  integer :: fid_TRN4S    !total solute NH4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRN3S    !total solute NH3 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRH1P    !total solute HPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRH2P    !total solute H2PO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRXN4    !total adsorbed NH4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRXH1    !total adsorbed OH transformation non-band, [mol d-2 h-1]
  integer :: fid_TRXH2    !total adsorbed OH2 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRX1P    !total adsorbed HPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRX2P    !total adsorbed H2PO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRBH1    !total adsorbed OH transformation band, [mol d-2 h-1]
  integer :: fid_TRBH2    !total adsorbed OH2 transformation band, [mol d-2 h-1]
  integer :: fid_TRB1P    !total adsorbed HPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRB2P    !total adsorbed H2PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRALPO   !total precipitated AlPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRFEPO   !total precipitated FePO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRCAPD   !total precipitated CaHPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRCAPH   !total precipitated CaH2PO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRCAPM   !total precipitated apatite transformation non-band, [mol d-2 h-1]
  integer :: fid_TRAL     !total solute Al transformation, [mol d-2 h-1]
  integer :: fid_TRFE     !total solute Fe transformation, [mol d-2 h-1]
  integer :: fid_TRHY     !total solute H transformation, [mol d-2 h-1]
  integer :: fid_TRCA     !total solute Ca transformation, [mol d-2 h-1]
  integer :: fid_TRMG     !total solute Mg transformation, [mol d-2 h-1]
  integer :: fid_TRNA     !total solute Na transformation, [mol d-2 h-1]
  integer :: fid_TRKA     !total solute K transformation, [mol d-2 h-1]
  integer :: fid_TROH     !total solute OH transformation, [mol d-2 h-1]
  integer :: fid_TRSO4    !total solute SO4 transformation, [mol d-2 h-1]
  integer :: fid_TRCO3    !total solute CO3 transformation, [mol d-2 h-1]
  integer :: fid_TRHCO    !total solute HCO3 transformation, [mol d-2 h-1]
  integer :: fid_TRCO2    !total solute CO2 transformation, [mol d-2 h-1]
  integer :: fid_TRAL1    !total solute AlOH transformation, [mol d-2 h-1]
  integer :: fid_TRAL2    !total solute AlOH2 transformation, [mol d-2 h-1]
  integer :: fid_TRAL3    !total solute AlOH3 transformation, [mol d-2 h-1]
  integer :: fid_TRAL4    !total solute AlOH4 transformation, [mol d-2 h-1]
  integer :: fid_TRALS    !total solute AlSO4 transformation, [mol d-2 h-1]
  integer :: fid_TRFE1    !total solute FeOH transformation, [mol d-2 h-1]
  integer :: fid_TRFE2    !total solute FeOH2 transformation, [mol d-2 h-1]
  integer :: fid_TRFE3    !total solute FeOH3 transformation, [mol d-2 h-1]
  integer :: fid_TRFE4    !total solute FeOH4 transformation, [mol d-2 h-1]
  integer :: fid_TRFES    !total solute FeSO4 transformation, [mol d-2 h-1]
  integer :: fid_TRCAO    !total solute CaOH transformation, [mol d-2 h-1]
  integer :: fid_TRCAC    !total solute CaCO3 transformation, [mol d-2 h-1]
  integer :: fid_TRCAH    !total solute CaHCO3 transformation, [mol d-2 h-1]
  integer :: fid_TRCAS    !total solute CaSO4 transformation, [mol d-2 h-1]
  integer :: fid_TRMGO    !total solute MgOH transformation, [mol d-2 h-1]
  integer :: fid_TRMGC    !total solute MgCO3 transformation, [mol d-2 h-1]
  integer :: fid_TRMGH    !total solute MgHCO3(+) transformation, [mol d-2 h-1]
  integer :: fid_TRMGS    !total solute MgSO4 transformation, [mol d-2 h-1]
  integer :: fid_TRNAC    !total solute NaCO3(-) transformation, [mol d-2 h-1]
  integer :: fid_TRNAS    !total solute NaSO4(-) transformation, [mol d-2 h-1]
  integer :: fid_TRKAS    !total solute KSO4(-) transformation, [mol d-2 h-1]
  integer :: fid_TRH0P    !total solute PO4(---) transformation non-band, [mol d-2 h-1]
  integer :: fid_TRH3P    !total solute H3PO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRF1P    !total solute FeHPO4(+) transformation non-band, [mol d-2 h-1]
  integer :: fid_TRF2P    !total solute FeH2PO4(++) transformation non-band, [mol d-2 h-1]
  integer :: fid_TRC0P    !total solute CaPO4(-) transformation non-band, [mol d-2 h-1]
  integer :: fid_TRC1P    !total solute CaHPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRC2P    !total solute CaH2PO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRM1P    !total solute MgHPO4 transformation non-band, [mol d-2 h-1]
  integer :: fid_TRXHY    !total adsorbed H transformation, [mol d-2 h-1]
  integer :: fid_TRXAL    !total adsorbed Al transformation, [mol d-2 h-1]
  integer :: fid_TRXFE    !total Fe adsorption
  integer :: fid_TRXCA    !total adsorbed Ca transformation, [mol d-2 h-1]
  integer :: fid_TRXMG    !total adsorbed Mg transformation, [mol d-2 h-1]
  integer :: fid_TRXNA    !total adsorbed Na transformation, [mol d-2 h-1]
  integer :: fid_TRXKA    !total adsorbed K transformation, [mol d-2 h-1]
  integer :: fid_TRXHC    !total adsorbed COOH transformation, [mol d-2 h-1]
  integer :: fid_TRXAL2   !total adsorbed AlOH2 transformation, [mol d-2 h-1]
  integer :: fid_TRXFE2   !total FeOH2 adsorption, [mol d-2 h-1]
  integer :: fid_TRXH0    !total adsorbed OH- transformation non-band, [mol d-2 h-1]
  integer :: fid_TRBH0    !total adsorbed OH- transformation band, [mol d-2 h-1]
  integer :: fid_TRALOH   !total precipitated AlOH3 transformation, [mol d-2 h-1]
  integer :: fid_TRFEOH   !total precipitated FeOH3 transformation, [mol d-2 h-1]
  integer :: fid_TRCACO   !total precipitated CaCO3 transformation, [mol d-2 h-1]
  integer :: fid_TRCASO   !total precipitated CaSO4 transformation, [mol d-2 h-1]
  integer :: fid_TRH2O    !total solute H2O transformation, [mol d-2 h-1]
  integer :: fid_TBION    !total solute ion transformation, [mol d-2 h-1]
  integer :: fid_TBCO2    !CO2 net change from all solute equilibria

  integer :: fid_TRN4B    !total solute NH4 transformation band, [mol d-2 h-1]
  integer :: fid_TRN3B    !total solute NH3 transformation band, [mol d-2 h-1]
  integer :: fid_TRH1B    !total solute HPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRH2B    !total solute H2PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRH0B    !total solute PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRH3B    !total solute H3PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRF1B    !total solute FeHPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRF2B    !total solute FeH2PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRC0B    !total solute CaPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRC1B    !total solute CaHPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRC2B    !total solute CaH2PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRM1B    !total solute MgHPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRXNB    !total adsorbed NH4 transformation band, [mol d-2 h-1]
  integer :: fid_TRALPB   !total precipitated AlPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRFEPB   !total precipitated FePO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRCPDB   !total precipitated CaHPO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRCPHB   !total precipitated CaH2PO4 transformation band, [mol d-2 h-1]
  integer :: fid_TRCPMB   !total precipitated apatite transformation band, [mol d-2 h-1]
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

  varl(cid_CH1P1)='CH1P1';varlnml(cid_CH1P1)='non-band soil aqueous HPO4 content micropore'
  unitl(cid_CH1P1)='mol m-3';vartypes(cid_CH1P1)=var_state_type

  varl(cid_CH1PB)='CH1PB';varlnml(cid_CH1PB)='band soil aqueous HPO4 content micropore'
  unitl(cid_CH1PB)='mol m-3';vartypes(cid_CH1PB)=var_state_type

  varl(cid_CH2P1)='CH2P1';varlnml(cid_CH2P1)='non-band soil micropore aqueous H2PO4 content'
  unitl(cid_CH2P1)='mol m-3';vartypes(cid_CH2P1)=var_state_type

  varl(cid_CH2PB)='CH2PB';varlnml(cid_CH2PB)='band soil micropore aqueous H2PO4 content'
  unitl(cid_CH2PB)='mol m-3';vartypes(cid_CH2PB)=var_state_type

  varl(cid_CN31) ='CN31';varlnml(cid_CN31)='non-band soil NH3 concentration'
  unitl(cid_CN31)='mol m-3';vartypes(cid_CN31)=var_state_type

  varl(cid_CN3B) ='CN3B';varlnml(cid_CN3B)='band soil NH3 concentration'
  unitl(cid_CN3B)='mol m-3';vartypes(cid_CN3B)=var_state_type

  varl(cid_CN41)  ='CN41';varlnml(cid_CN41)='non-band soil NH4 concentration'
  unitl(cid_CN41)='mol m-3';vartypes(cid_CN41)=var_state_type

  varl(cid_CN4B)  ='CN4B';varlnml(cid_CN4B)='band soil NH4 concentration'
  unitl(cid_CN4B) ='mol m-3';vartypes(cid_CN4B)=var_state_type

  varl(cid_XN41)  ='XN41';varlnml(cid_XN41)='non-band soil adsorbed NH4 concentration'
  unitl(cid_XN41) ='mol m-3';vartypes(cid_XN41)=var_state_type

  varl(cid_XN4B)  ='XN4B';varlnml(cid_XN4B)='band soil adsorbed NH4 concentration'
  unitl(cid_XN4B) ='mol m-3';vartypes(cid_XN4B)=var_state_type

  varl(cid_X1P1B) ='X1P1B';varlnml(cid_X1P1B)='band soil exchangeable HPO4 concentration'
  unitl(cid_X1P1B)='mol m-3';vartypes(cid_X1P1B)=var_state_type

  varl(cid_X2P1B) ='X2P1B';varlnml(cid_X2P1B)='band soil exchangeable H2PO4 concentration'
  unitl(cid_X2P1B)='mol m-3';vartypes(cid_X2P1B)=var_state_type

  varl(cid_XH11B) ='XH11B';varlnml(cid_XH11B)='band soil exchangeable site R-OH'
  unitl(cid_XH11B)='mol m-3';vartypes(cid_XH11B)=var_state_type

  varl(cid_XH1P1) ='XH1P1';varlnml(cid_XH1P1)='non-band soil concentration of adsorbed HPO4'
  unitl(cid_XH1P1)='mol m-3';vartypes(cid_XH1P1)=var_state_type

  varl(cid_XH21B) ='XH21B';varlnml(cid_XH21B)='band soil exchangeable site R-OH2'
  unitl(cid_XH21B)='mol m-3';vartypes(cid_XH21B)=var_state_type

  varl(cid_XH2P1) ='XH2P1';varlnml(cid_XH2P1)='non-band soil concentration of adsorbed HPO4'
  unitl(cid_XH2P1)='mol m-3';vartypes(cid_XH2P1)=var_state_type

  varl(cid_XOH11) ='XOH11';varlnml(cid_XOH11)='non-band soil concentration of adsorption sites R-OH'
  unitl(cid_XOH11)='mol m-3';vartypes(cid_XOH11)=var_state_type

  varl(cid_XOH21) ='XOH21';varlnml(cid_XOH21)='non-band soil exchangeable site R-OH2'
  unitl(cid_XOH21)='mol m-3';vartypes(cid_XOH21)=var_state_type

  varl(cid_PALPO1)='PALPO1';varlnml(cid_PALPO1)='non-band soil precipitated AlPO4'
  unitl(cid_PALPO1)='mol m-3';vartypes(cid_PALPO1)=var_state_type

  varl(cid_PALPOB) ='PALPOB';varlnml(cid_PALPOB)='band soil precipitated AlPO4'
  unitl(cid_PALPOB)='mol m-3';vartypes(cid_PALPOB)=var_state_type

  varl(cid_PCAPD1) ='PCAPD1';varlnml(cid_PCAPD1)='non-band soil precipitated CaHPO4'
  unitl(cid_PCAPD1)='mol m-3';vartypes(cid_PCAPD1)=var_state_type

  varl(cid_PCAPDB) ='PCAPDB';varlnml(cid_PCAPDB)='band soil precipitated CaHPO4'
  unitl(cid_PCAPDB)='mol m-3';vartypes(cid_PCAPDB)=var_state_type

  varl(cid_PCAPH1) ='PCAPH1';varlnml(cid_PCAPH1)='non-band soil precipitated Ca5(PO4)3OH hydroxyapatite'
  unitl(cid_PCAPH1)='mol m-3';vartypes(cid_PCAPH1)=var_state_type

  varl(cid_PCAPHB) ='PCAPHB';varlnml(cid_PCAPHB)='band soil precipitated Ca5(PO4)3OH hydroxyapatite'
  unitl(cid_PCAPHB)='mol m-3';vartypes(cid_PCAPHB)=var_state_type

  varl(cid_PCAPM1) ='PCAPM1';varlnml(cid_PCAPM1)='non-band soil precipitated Ca(H2PO4)2'
  unitl(cid_PCAPM1)='mol m-3';vartypes(cid_PCAPM1)=var_state_type

  varl(cid_PCAPMB) ='PCAPMB';varlnml(cid_PCAPMB)='band soil precipitated CaH2PO4'
  unitl(cid_PCAPMB)='mol m-3';vartypes(cid_PCAPMB)=var_state_type

  varl(cid_PFEPO1) ='PFEPO1';varlnml(cid_PFEPO1)='non-band soil precipitated FePO4'
  unitl(cid_PFEPO1)='mol m-3';vartypes(cid_PFEPO1)=var_state_type

  varl(cid_PFEPOB) ='PFEPOB';varlnml(cid_PFEPOB)='band soil precipitated FePO4'
  unitl(cid_PFEPOB)='mol m-3';vartypes(cid_PFEPOB)=var_state_type

  varl(fid_TRN4S) = 'TRN4S';varlnml(fid_TRN4S)='non-band soil total solute NH4 transformation'
  unitl(fid_TRN4S)= 'mol d-2 h-1';vartypes(fid_TRN4S)=var_flux_type

  varl(fid_TRN4B) = 'TRN4B';varlnml(fid_TRN4B)='band soil total solute NH4 transformation'
  unitl(fid_TRN4B)= 'mol d-2 h-1';vartypes(fid_TRN4B)=var_flux_type

  varl(fid_TRN3S) = 'TRN3S';varlnml(fid_TRN3S)='non-band total solute NH3 transformation'
  unitl(fid_TRN3S)= 'mol d-2 h-1';vartypes(fid_TRN3S)=var_flux_type

  varl(fid_TRN3B) = 'TRN3B';varlnml(fid_TRN3B)='band soil total solute NH3 transformation'
  unitl(fid_TRN3B)= 'mol d-2 h-1';vartypes(fid_TRN3B)=var_flux_type

  varl(fid_TRH1P) = 'TRH1P';varlnml(fid_TRH1P)='non-band soil total solute HPO4 transformation'
  unitl(fid_TRH1P) = 'mol d-2 h-1';vartypes(fid_TRH1P)=var_flux_type

  varl(fid_TRH2P) = 'TRH2P';varlnml(fid_TRH2P)='non-band soil total solute H2PO4 transformation'
  unitl(fid_TRH2P)= 'mol d-2 h-1';vartypes(fid_TRH2P)=var_flux_type

  varl(fid_TRH1B) = 'TRH1B';varlnml(fid_TRH1B)='band soil total solute HPO4 transformation'
  unitl(fid_TRH1B)= 'mol d-2 h-1';vartypes(fid_TRH1B)=var_flux_type

  varl(fid_TRH2B) = 'TRH2B';varlnml(fid_TRH2B)='band soil total solute H2PO4 transformation'
  unitl(fid_TRH2B)= 'mol d-2 h-1';vartypes(fid_TRH2B)=var_flux_type

  varl(fid_TRXN4) = 'TRXN4';varlnml(fid_TRXN4)='non-band soil total adsorbed NH4 transformation'
  unitl(fid_TRXN4)= 'mol d-2 h-1';vartypes(fid_TRXN4)=var_flux_type

  varl(fid_TRXNB) = 'TRXNB';varlnml(fid_TRXNB)='band soil total adsorbed NH4 transformation'
  unitl(fid_TRXNB)= 'mol d-2 h-1';vartypes(fid_TRXNB)=var_flux_type

  varl(fid_TRXH1) = 'TRXH1';varlnml(fid_TRXH1)='non-band soil total adsorbed OH transformation'
  unitl(fid_TRXH1)= 'mol d-2 h-1';vartypes(fid_TRXH1)=var_flux_type

  varl(fid_TRXH2) = 'TRXH2';varlnml(fid_TRXH2)='non-band soil total adsorbed OH2 transformation'
  unitl(fid_TRXH2)= 'mol d-2 h-1';vartypes(fid_TRXH2)=var_flux_type

  varl(fid_TRX1P) = 'TRX1P';varlnml(fid_TRX1P)='non-band soil total adsorbed HPO4 transformation'
  unitl(fid_TRX1P)= 'mol d-2 h-1';vartypes(fid_TRX1P)=var_flux_type

  varl(fid_TRX2P) = 'TRX2P';varlnml(fid_TRX2P)='non-band soil total adsorbed H2PO4 transformation'
  unitl(fid_TRX2P)= 'mol d-2 h-1';vartypes(fid_TRX2P)=var_flux_type

  varl(fid_TRBH1) = 'TRBH1';varlnml(fid_TRBH1)='band soil total adsorbed OH transformation'
  unitl(fid_TRBH1)= 'mol d-2 h-1';vartypes(fid_TRBH1)=var_flux_type

  varl(fid_TRBH2) = 'TRBH2';varlnml(fid_TRBH2)='band soil total adsorbed OH2 transformation'
  unitl(fid_TRBH2)= 'mol d-2 h-1';vartypes(fid_TRBH2)=var_flux_type

  varl(fid_TRB1P) = 'TRB1P';varlnml(fid_TRB1P)='band soil total adsorbed HPO4 transformation'
  unitl(fid_TRB1P)= 'mol d-2 h-1';vartypes(fid_TRB1P)=var_flux_type

  varl(fid_TRB2P) = 'TRB2P';varlnml(fid_TRB2P)='band soil total adsorbed H2PO4 transformation'
  unitl(fid_TRB2P)= 'mol d-2 h-1';vartypes(fid_TRB2P)=var_flux_type

  varl(fid_TRALPO)= 'TRALPO';varlnml(fid_TRALPO)='non-band total precipitated AlPO4 transformation'
  unitl(fid_TRALPO)= 'mol d-2 h-1';vartypes(fid_TRALPO)=var_flux_type

  varl(fid_TRFEPO)= 'TRFEPO';varlnml(fid_TRFEPO)='non-band soil total precipitated FePO4 transformation'
  unitl(fid_TRFEPO)= 'mol d-2 h-1';vartypes(fid_TRFEPO)=var_flux_type

  varl(fid_TRCAPD)= 'TRCAPD';varlnml(fid_TRCAPD)='non-band soil total precipitated CaHPO4 transformation'
  unitl(fid_TRCAPD)= 'mol d-2 h-1';vartypes(fid_TRCAPD)=var_flux_type

  varl(fid_TRCAPH)= 'TRCAPH';varlnml(fid_TRCAPH)='non-band soil total precipitated CaH2PO4 transformation'
  unitl(fid_TRCAPH)= 'mol d-2 h-1';vartypes(fid_TRCAPH)=var_flux_type

  varl(fid_TRCAPM)= 'TRCAPM';varlnml(fid_TRCAPM)='non-band soil total precipitated apatite transformation'
  unitl(fid_TRCAPM)= 'mol d-2 h-1';vartypes(fid_TRCAPM)=var_flux_type

  varl(fid_TRALPB)= 'TRALPB';varlnml(fid_TRALPB)='band soil total precipitated AlPO4 transformation'
  unitl(fid_TRALPB)= 'mol d-2 h-1';vartypes(fid_TRALPB)=var_flux_type

  varl(fid_TRFEPB)= 'TRFEPB';varlnml(fid_TRFEPB)='band soil total precipitated FePO4 transformation'
  unitl(fid_TRFEPB)= 'mol d-2 h-1';vartypes(fid_TRFEPB)=var_flux_type

  varl(fid_TRCPDB)= 'TRCPDB';varlnml(fid_TRCPDB)='band soil total precipitated CaHPO4 transformation'
  unitl(fid_TRCPDB)= 'mol d-2 h-1';vartypes(fid_TRCPDB)=var_flux_type

  varl(fid_TRCPHB)= 'TRCPHB';varlnml(fid_TRCPHB)='band soil total precipitated CaH2PO4 transformation'
  unitl(fid_TRCPHB)= 'mol d-2 h-1';vartypes(fid_TRCPHB)=var_flux_type

  varl(fid_TRCPMB)= 'TRCPMB';varlnml(fid_TRCPMB)='band soil total precipitated apatite transformation'
  unitl(fid_TRCPMB)= 'mol d-2 h-1';vartypes(fid_TRCPMB)=var_flux_type

  varl(fid_TRAL)  = 'TRAL';varlnml(fid_TRAL)='non-band soil total solute Al transformation'
  unitl(fid_TRAL) = 'mol d-2 h-1';vartypes(fid_TRAL)=var_flux_type

  end subroutine getvarlist_nosalt

end module ChemIDMod
