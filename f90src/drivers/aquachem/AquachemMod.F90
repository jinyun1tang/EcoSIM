module AquachemMod
  use data_kind_mod  , only : r8 => SHR_KIND_R8
  use MiscUtilMod    , only : addone
  use ModelStatusType, only : model_status_type
implicit none
  public
  character(len=*), private, parameter :: mod_filename = __FILE__

  integer, private :: cid_CO2S     !aqueous CO2  micropore	[g d-2]
  integer, private :: cid_CH1P1    !soil aqueous HPO4 content micropore non-band, [mol m-3]
  integer, private :: cid_CH2P1    !soil aqueous H2PO4 content micropore non-band, [mol m-3]
  integer, private :: cid_CN31     !soil NH3 concentration in non-band soil, [mol m-3]
  integer, private :: cid_CN41     !soil NH4 concentration in non-band soil, [mol m-3]
  integer, private :: cid_ZNO3S    !NO3 mass non-band micropore, [g d-2]
  integer, private :: cid_ZHY      !soil aqueous H content micropore, [mol d-2]
  integer, private :: cid_ZOH      !soil aqueous OH content micropore, [mol d-2]
  integer, private :: cid_ZAL      !soil aqueous Al content micropore, [mol d-2]
  integer, private :: cid_ZFE      !soil aqueous Fe content micropore, [mol d-2]
  integer, private :: cid_ZCA      !soil aqueous Ca content micropore, [mol d-2]
  integer, private :: cid_ZMG      !soil aqueous Mg content micropore, [mol d-2]
  integer, private :: cid_ZNA      !soil aqueous Na content micropore, [mol d-2]
  integer, private :: cid_ZKA      !soil aqueous K content micropore, [mol d-2]
  integer, private :: cid_ZSO4     !soil aqueous SO4 content micropore, [mol d-2]
  integer, private :: cid_ZCL      !soil aqueous Cl content micropore, [mol d-2]
  integer, private :: cid_ZCO3     !soil aqueous CO3 content micropore, [mol d-2]
  integer, private :: cid_ZHCO3    !soil aqueous HCO3 content micropore, [mol d-2]
  integer, private :: cid_ZALOH1   !soil aqueous AlOH content micropore, [mol d-2]
  integer, private :: cid_ZALOH2   !soil aqueous AlOH2 content micropore, [mol d-2]
  integer, private :: cid_ZALOH3   !soil aqueous AlOH3 content micropore, [mol d-2]
  integer, private :: cid_ZALOH4   !soil aqueous AlOH4 content micropore, [mol d-2]
  integer, private :: cid_ZALS     !soil aqueous AlSO4 content micropore, [mol d-2]
  integer, private :: cid_ZFEOH1   !soil aqueous FeOH content micropore, [mol d-2]
  integer, private :: cid_ZFEOH2   !soil aqueous FeOH2 content micropore, [mol d-2]
  integer, private :: cid_ZFEOH3   !soil aqueous FeOH3 content micropore, [mol d-2]
  integer, private :: cid_ZFEOH4   !soil aqueous FeOH4 content micropore, [mol d-2]
  integer, private :: cid_ZFES     !soil aqueous FeSO4 content micropore, [mol d-2]
  integer, private :: cid_ZCAO     !soil aqueous CaOH2 content micropore, [mol d-2]
  integer, private :: cid_ZCAC     !soil aqueous CACO3 content micropore, [mol d-2]
  integer, private :: cid_ZCAH     !soil aqueous CaHCO3 content micropore, [mol d-2]
  integer, private :: cid_ZCAS     !soil aqueous CaSO4 content micropore, [mol d-2]
  integer, private :: cid_ZMGO     !soil aqueous MgOH content micropore, [mol d-2]
  integer, private :: cid_ZMGC     !soil aqueous MgCO3 content micropore, [mol d-2]
  integer, private :: cid_ZMGH     !soil aqueous MgHCO3 content micropore, [mol d-2]
  integer, private :: cid_ZMGS     !soil aqueous MgSO4 content micropore, [mol d-2]
  integer, private :: cid_ZNAC     !soil aqueous NaCO3 content micropore, [mol d-2]
  integer, private :: cid_ZNAS     !soil aqueous NaSO4 content micropore, [mol d-2]
  integer, private :: cid_ZKAS     !soil aqueous KSO4 content micropore, [mol d-2]
  integer, private :: cid_H0PO4    !soil aqueous PO4 content micropore non-band, [mol d-2]
  integer, private :: cid_H3PO4    !soil aqueous H3PO4 content micropore non-band, [mol d-2]
  integer, private :: cid_ZFE1P    !soil aqueous FeHPO4 content micropore non-band, [mol d-2]
  integer, private :: cid_ZFE2P    !soil aqueous FeH2PO4 content micropore non-band, [mol d-2]
  integer, private :: cid_ZCA0P    !soil aqueous CaPO4 content micropore non-band, [mol d-2]
  integer, private :: cid_ZCA1P    !soil aqueous CaHPO4 content micropore non-band, [mol d-2]
  integer, private :: cid_ZCA2P    !soil aqueous CaH2PO4 content micropore non-band, [mol d-2]
  integer, private :: cid_ZMG1P    !soil aqueous MgHPO4 content micropore non-band, [mol d-2]
  integer, private :: cid_PALPO1   !precipitated AlPO4 non-band, [mol m-3]
  integer, private :: cid_PCAPD1   !precipitated CaHPO4 non-band soil, [mol m-3]
  integer, private :: cid_PCAPH1   !precipitated Ca5(PO4)3OH hydroxyapatite non-band soil, [mol m-3]
  integer, private :: cid_PCAPM1   !precipitated Ca(H2PO4)2 non-band soil, [mol m-3]
  integer, private :: cid_PFEPO1   !precipitated FePO4 non-band soil, [mol m-3]
  integer, private :: cid_PALOH    !precipitated Al(OH)3, [mol d-2]
  integer, private :: cid_PFEOH    !precipitated Fe(OH)3, [mol d-2]
  integer, private :: cid_PCACO    !precipitated CaCO3, [mol d-2]
  integer, private :: cid_PCASO    !precipitated CaSO4, [mol d-2]
  integer, private :: cid_XHY      !exchangeable H(+) , [mol d-2]
  integer, private :: cid_XAL      !exchangeable Al, [mol d-2]
  integer, private :: cid_XFE      !exchangeable Fe, [mol d-2]
  integer, private :: cid_XCA      !exchangeable Ca, [mol d-2]
  integer, private :: cid_XMG      !exchangeable Mg, [mol d-2]
  integer, private :: cid_XNA      !exchangeable Na, [mol d-2]
  integer, private :: cid_XKA      !exchangeable K, [mol d-2]
  integer, private :: cid_XHC      !exchangeable COOH , [mol d-2]
  integer, private :: cid_XOH11    !exchangeable OH  non-band, [mol d-2]
  integer, private :: cid_XALO2    !exchangeable AlOH2 , [mol d-2]
  integer, private :: cid_XFEO2    !exchangeable Fe(OH)2, [mol d-2]
  integer, private :: cid_XN41     !exchangeable NH4 non-band soil, [mol d-2]
  integer, private :: cid_XOH01    !exchangeable OH- non-band, [mol d-2]
  integer, private :: cid_XH1P1    !exchangeable HPO4  non-band, [mol m-3]
  integer, private :: cid_XOH21    !exchangeable OH2  non-band soil, [mol m-3]
  integer, private :: cid_XH2P1    !exchangeable H2PO4  non-band soil, [mol m-3]

  integer, private :: cid_CN4B     !soil NH4 concentration in band soil, [mol m-3]
  integer, private :: cid_CN3B     !soil NH3 concentration in band soil, [mol m-3]
  integer, private :: cid_ZNO3B    !NO3 mass band micropore, [g d-2]
  integer, private :: cid_CH1PB    !soil aqueous HPO4 content micropore band, [mol m-3]
  integer, private :: cid_CH2PB    !soil aqueous H2PO4 content micropore  band, [mol m-3]
  integer, private :: cid_H0POB    !soil aqueous PO4 content micropore band, [mol d-2]
  integer, private :: cid_H3POB    !soil aqueous H3PO4 content micropore band, [mol d-2]
  integer, private :: cid_ZFE1PB   !soil aqueous FeHPO4 content micropore band, [mol d-2]
  integer, private :: cid_ZFE2PB   !soil aqueous FeH2PO4 content micropore band, [mol d-2
  integer, private :: cid_ZCA0PB   !soil aqueous CaPO4 content micropore band, [mol d-2]
  integer, private :: cid_ZCA1PB   !soil aqueous CaHPO4 content micropore band, [mol d-2]
  integer, private :: cid_ZCA2PB   !soil aqueous CaH2PO4 content micropore band, [mol d-2]
  integer, private :: cid_ZMG1PB   !soil aqueous MgHPO4 content micropore band, [mol d-2]
  integer, private :: cid_PALPOB   !precipitated AlPO4 band soil, [mol m-3]
  integer, private :: cid_PCAPDB   !precipitated CaHPO4 band soil, [mol m-3]
  integer, private :: cid_PCAPHB   !precipitated Ca5(PO4)3OH hydroxyapatite band soil, [mol m-3]
  integer, private :: cid_PCAPMB   !precipitated CaH2PO4 band soil, [mol m-3]
  integer, private :: cid_PFEPOB   !precipitated FePO4 band soil, [mol m-3]
  integer, private :: cid_XN4B     !exchangeable NH4 band soil, [mol d-2]
  integer, private :: cid_XH01B    !exchangeable OH- band, [mol d-2]
  integer, private :: cid_X1P1B    !exchangeable HPO4 concentration band-soil, [mol m-3]
  integer, private :: cid_X2P1B    !exchangeable H2PO4 concentration band-soil, [mol m-3]
  integer, private :: cid_XH11B    !exchangeable OH band-soil, [mol m-3]
  integer, private :: cid_XH21B    !exchangeable OH2 band-soil, [mol m-3]

  integer, private :: fid_TRN4S    !total solute NH4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRN3S    !total solute NH3 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRH1P    !total solute HPO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRH2P    !total solute H2PO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRXN4    !total adsorbed NH4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRXH1    !total adsorbed OH transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRXH2    !total adsorbed OH2 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRX1P    !total adsorbed HPO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRX2P    !total adsorbed H2PO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRBH1    !total adsorbed OH transformation band, [mol d-2 h-1]
  integer, private :: fid_TRBH2    !total adsorbed OH2 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRB1P    !total adsorbed HPO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRB2P    !total adsorbed H2PO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRALPO   !total precipitated AlPO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRFEPO   !total precipitated FePO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRCAPD   !total precipitated CaHPO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRCAPH   !total precipitated CaH2PO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRCAPM   !total precipitated apatite transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRAL     !total solute Al transformation, [mol d-2 h-1]
  integer, private :: fid_TRFE     !total solute Fe transformation, [mol d-2 h-1]
  integer, private :: fid_TRHY     !total solute H transformation, [mol d-2 h-1]
  integer, private :: fid_TRCA     !total solute Ca transformation, [mol d-2 h-1]
  integer, private :: fid_TRMG     !total solute Mg transformation, [mol d-2 h-1]
  integer, private :: fid_TRNA     !total solute Na transformation, [mol d-2 h-1]
  integer, private :: fid_TRKA     !total solute K transformation, [mol d-2 h-1]
  integer, private :: fid_TROH     !total solute OH transformation, [mol d-2 h-1]
  integer, private :: fid_TRSO4    !total solute SO4 transformation, [mol d-2 h-1]
  integer, private :: fid_TRCO3    !total solute CO3 transformation, [mol d-2 h-1]
  integer, private :: fid_TRHCO    !total solute HCO3 transformation, [mol d-2 h-1]
  integer, private :: fid_TRCO2    !total solute CO2 transformation, [mol d-2 h-1]
  integer, private :: fid_TRAL1    !total solute AlOH transformation, [mol d-2 h-1]
  integer, private :: fid_TRAL2    !total solute AlOH2 transformation, [mol d-2 h-1]
  integer, private :: fid_TRAL3    !total solute AlOH3 transformation, [mol d-2 h-1]
  integer, private :: fid_TRAL4    !total solute AlOH4 transformation, [mol d-2 h-1]
  integer, private :: fid_TRALS    !total solute AlSO4 transformation, [mol d-2 h-1]
  integer, private :: fid_TRFE1    !total solute FeOH transformation, [mol d-2 h-1]
  integer, private :: fid_TRFE2    !total solute FeOH2 transformation, [mol d-2 h-1]
  integer, private :: fid_TRFE3    !total solute FeOH3 transformation, [mol d-2 h-1]
  integer, private :: fid_TRFE4    !total solute FeOH4 transformation, [mol d-2 h-1]
  integer, private :: fid_TRFES    !total solute FeSO4 transformation, [mol d-2 h-1]
  integer, private :: fid_TRCAO    !total solute CaOH transformation, [mol d-2 h-1]
  integer, private :: fid_TRCAC    !total solute CaCO3 transformation, [mol d-2 h-1]
  integer, private :: fid_TRCAH    !total solute CaHCO3 transformation, [mol d-2 h-1]
  integer, private :: fid_TRCAS    !total solute CaSO4 transformation, [mol d-2 h-1]
  integer, private :: fid_TRMGO    !total solute MgOH transformation, [mol d-2 h-1]
  integer, private :: fid_TRMGC    !total solute MgCO3 transformation, [mol d-2 h-1]
  integer, private :: fid_TRMGH    !total solute MgHCO3(+) transformation, [mol d-2 h-1]
  integer, private :: fid_TRMGS    !total solute MgSO4 transformation, [mol d-2 h-1]
  integer, private :: fid_TRNAC    !total solute NaCO3(-) transformation, [mol d-2 h-1]
  integer, private :: fid_TRNAS    !total solute NaSO4(-) transformation, [mol d-2 h-1]
  integer, private :: fid_TRKAS    !total solute KSO4(-) transformation, [mol d-2 h-1]
  integer, private :: fid_TRH0P    !total solute PO4(---) transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRH3P    !total solute H3PO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRF1P    !total solute FeHPO4(+) transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRF2P    !total solute FeH2PO4(++) transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRC0P    !total solute CaPO4(-) transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRC1P    !total solute CaHPO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRC2P    !total solute CaH2PO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRM1P    !total solute MgHPO4 transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRXHY    !total adsorbed H transformation, [mol d-2 h-1]
  integer, private :: fid_TRXAL    !total adsorbed Al transformation, [mol d-2 h-1]
  integer, private :: fid_TRXFE    !total Fe adsorption
  integer, private :: fid_TRXCA    !total adsorbed Ca transformation, [mol d-2 h-1]
  integer, private :: fid_TRXMG    !total adsorbed Mg transformation, [mol d-2 h-1]
  integer, private :: fid_TRXNA    !total adsorbed Na transformation, [mol d-2 h-1]
  integer, private :: fid_TRXKA    !total adsorbed K transformation, [mol d-2 h-1]
  integer, private :: fid_TRXHC    !total adsorbed COOH transformation, [mol d-2 h-1]
  integer, private :: fid_TRXAL2   !total adsorbed AlOH2 transformation, [mol d-2 h-1]
  integer, private :: fid_TRXFE2   !total FeOH2 adsorption, [mol d-2 h-1]
  integer, private :: fid_TRXH0    !total adsorbed OH- transformation non-band, [mol d-2 h-1]
  integer, private :: fid_TRBH0    !total adsorbed OH- transformation band, [mol d-2 h-1]
  integer, private :: fid_TRALOH   !total precipitated AlOH3 transformation, [mol d-2 h-1]
  integer, private :: fid_TRFEOH   !total precipitated FeOH3 transformation, [mol d-2 h-1]
  integer, private :: fid_TRCACO   !total precipitated CaCO3 transformation, [mol d-2 h-1]
  integer, private :: fid_TRCASO   !total precipitated CaSO4 transformation, [mol d-2 h-1]
  integer, private :: fid_TRH2O    !total solute H2O transformation, [mol d-2 h-1]
  integer, private :: fid_TBION    !total solute ion transformation, [mol d-2 h-1]
  integer, private :: fid_TBCO2    !CO2 net change from all solute equilibria

  integer, private :: fid_TRN4B    !total solute NH4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRN3B    !total solute NH3 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRH1B    !total solute HPO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRH2B    !total solute H2PO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRH0B    !total solute PO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRH3B    !total solute H3PO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRF1B    !total solute FeHPO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRF2B    !total solute FeH2PO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRC0B    !total solute CaPO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRC1B    !total solute CaHPO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRC2B    !total solute CaH2PO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRM1B    !total solute MgHPO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRXNB    !total adsorbed NH4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRALPB   !total precipitated AlPO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRFEPB   !total precipitated FePO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRCPDB   !total precipitated CaHPO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRCPHB   !total precipitated CaH2PO4 transformation band, [mol d-2 h-1]
  integer, private :: fid_TRCPMB   !total precipitated apatite transformation band, [mol d-2 h-1]

  private :: Init_Salt_geochem, Init_NoSalt_geochem
  private :: getvarlist_salt, getvarlist_nosalt
  private :: initmodel_salt, initmodel_nosalt
contains
  function getvarllen(salton)result(nvars)
  implicit none
  logical, intent(in) :: salton
  integer :: nvars

  if (salton)then
    call Init_Salt_geochem(nvars)
  else
    call Init_NoSalt_geochem(nvars)
  endif
  end function getvarllen


! ----------------------------------------------------------------------

  subroutine Init_NoSalt_geochem(nvars)

  implicit none
  integer, intent(out) :: nvars
  integer :: itemp
  itemp=0
  cid_ZMG    =addone(itemp)
  cid_ZNA    =addone(itemp)
  cid_ZKA    =addone(itemp)
  cid_CO2S   =addone(itemp)
  cid_CH1P1  =addone(itemp)
  cid_CH1PB  =addone(itemp)
  cid_CH2P1  =addone(itemp)
  cid_CH2PB  =addone(itemp)
  cid_CN31   =addone(itemp)
  cid_CN3B   =addone(itemp)
  cid_CN41   =addone(itemp)
  cid_CN4B   =addone(itemp)
  cid_XN41   =addone(itemp)
  cid_XN4B   =addone(itemp)
  cid_X1P1B  =addone(itemp)
  cid_X2P1B  =addone(itemp)
  cid_XH11B  =addone(itemp)
  cid_XH1P1  =addone(itemp)
  cid_XH21B  =addone(itemp)
  cid_XH2P1  =addone(itemp)
  cid_XOH11  =addone(itemp)
  cid_XOH21  =addone(itemp)
  cid_PALPO1 =addone(itemp)
  cid_PALPOB =addone(itemp)
  cid_PCAPD1 =addone(itemp)
  cid_PCAPDB =addone(itemp)
  cid_PCAPH1 =addone(itemp)
  cid_PCAPHB =addone(itemp)
  cid_PCAPM1 =addone(itemp)
  cid_PCAPMB =addone(itemp)
  cid_PFEPO1 =addone(itemp)
  cid_PFEPOB =addone(itemp)

  fid_TRN4S = addone(itemp)
  fid_TRN4B = addone(itemp)
  fid_TRN3S = addone(itemp)
  fid_TRN3B = addone(itemp)
  fid_TRH1P = addone(itemp)
  fid_TRH2P = addone(itemp)
  fid_TRH1B = addone(itemp)
  fid_TRH2B = addone(itemp)
  fid_TRXN4 = addone(itemp)
  fid_TRXNB = addone(itemp)
  fid_TRXH1 = addone(itemp)
  fid_TRXH2 = addone(itemp)
  fid_TRX1P = addone(itemp)
  fid_TRX2P = addone(itemp)
  fid_TRBH1 = addone(itemp)
  fid_TRBH2 = addone(itemp)
  fid_TRB1P = addone(itemp)
  fid_TRB2P = addone(itemp)
  fid_TRALPO= addone(itemp)
  fid_TRFEPO= addone(itemp)
  fid_TRCAPD= addone(itemp)
  fid_TRCAPH= addone(itemp)
  fid_TRCAPM= addone(itemp)
  fid_TRALPB= addone(itemp)
  fid_TRFEPB= addone(itemp)
  fid_TRCPDB= addone(itemp)
  fid_TRCPHB= addone(itemp)
  fid_TRCPMB= addone(itemp)
  fid_TRAL  = addone(itemp)
  nvars = itemp
  end subroutine Init_NoSalt_geochem

!------------------------------------------------------------------
  subroutine Init_Salt_geochem(nvars)
  implicit none
  integer, intent(out) :: nvars
  integer :: itemp
  itemp=0
  cid_CO2S     =addone(itemp)  !aqueous CO2  micropore	[g d-2]
  cid_CH1P1    =addone(itemp)  !soil aqueous HPO4 content micropore non-band, [mol m-3]
  cid_CH2P1    =addone(itemp)  !soil aqueous H2PO4 content micropore non-band, [mol m-3]
  cid_CN31     =addone(itemp)  !soil NH3 concentration in non-band soil, [mol m-3]
  cid_CN41     =addone(itemp)  !soil NH4 concentration in non-band soil, [mol m-3]
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
  cid_ZCA2P    =addone(itemp)  !soil aqueous CaH2PO4 content micropore non-band, [mol d-2]
  cid_ZMG1P    =addone(itemp)  !soil aqueous MgHPO4 content micropore non-band, [mol d-2]
  cid_PALPO1   =addone(itemp)  !precipitated AlPO4 non-band, [mol m-3]
  cid_PCAPD1   =addone(itemp)  !precipitated CaHPO4 non-band soil, [mol m-3]
  cid_PCAPH1   =addone(itemp)  !precipitated Ca5(PO4)3OH hydroxyapatite non-band soil, [mol m-3]
  cid_PCAPM1   =addone(itemp)  !precipitated Ca(H2PO4)2 non-band soil, [mol m-3]
  cid_PFEPO1   =addone(itemp)  !precipitated FePO4 non-band soil, [mol m-3]
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
  cid_XOH11    =addone(itemp)  !exchangeable OH  non-band, [mol d-2]
  cid_XALO2    =addone(itemp)  !exchangeable AlOH2 , [mol d-2]
  cid_XFEO2    =addone(itemp)  !exchangeable Fe(OH)2, [mol d-2]
  cid_XN41     =addone(itemp)  !exchangeable NH4 non-band soil, [mol d-2]
  cid_XOH01    =addone(itemp)  !exchangeable OH- non-band, [mol d-2]
  cid_XH1P1    =addone(itemp)  !exchangeable HPO4  non-band, [mol m-3]
  cid_XOH21    =addone(itemp)  !exchangeable OH2  non-band soil, [mol m-3]
  cid_XH2P1    =addone(itemp)  !exchangeable H2PO4  non-band soil, [mol m-3]

  cid_CH1PB    =addone(itemp)  !soil aqueous HPO4 content micropore band, [mol m-3]
  cid_CH2PB    =addone(itemp)  !soil aqueous H2PO4 content micropore  band, [mol m-3]
  cid_CN3B     =addone(itemp)  !soil NH3 concentration in band soil, [mol m-3]
  cid_CN4B     =addone(itemp)  !soil NH4 concentration in band soil, [mol m-3]
  cid_ZNO3B    =addone(itemp)  !NO3 mass band micropore, [g d-2]
  cid_H0POB    =addone(itemp)  !soil aqueous PO4 content micropore band, [mol d-2]
  cid_H3POB    =addone(itemp)  !soil aqueous H3PO4 content micropore band, [mol d-2]
  cid_ZFE1PB   =addone(itemp)  !soil aqueous FeHPO4 content micropore band, [mol d-2]
  cid_ZFE2PB   =addone(itemp)  !soil aqueous FeH2PO4 content micropore band, [mol d-2
  cid_ZCA0PB   =addone(itemp)  !soil aqueous CaPO4 content micropore band, [mol d-2]
  cid_ZCA1PB   =addone(itemp)  !soil aqueous CaHPO4 content micropore band, [mol d-2]
  cid_ZCA2PB   =addone(itemp)  !soil aqueous CaH2PO4 content micropore band, [mol d-2]
  cid_ZMG1PB   =addone(itemp)  !soil aqueous MgHPO4 content micropore band, [mol d-2]
  cid_PALPOB   =addone(itemp)  !precipitated AlPO4 band soil, [mol m-3]
  cid_PCAPDB   =addone(itemp)  !precipitated CaHPO4 band soil, [mol m-3]
  cid_PCAPHB   =addone(itemp)  !precipitated Ca5(PO4)3OH hydroxyapatite band soil, [mol m-3]
  cid_PCAPMB   =addone(itemp)  !precipitated CaH2PO4 band soil, [mol m-3]
  cid_PFEPOB   =addone(itemp)  !precipitated FePO4 band soil, [mol m-3]
  cid_XN4B     =addone(itemp)  !exchangeable NH4 band soil, [mol d-2]
  cid_XH01B    =addone(itemp)  !exchangeable OH- band, [mol d-2]
  cid_X1P1B    =addone(itemp)  !exchangeable HPO4 concentration band-soil, [mol m-3]
  cid_X2P1B    =addone(itemp)  !exchangeable H2PO4 concentration band-soil, [mol m-3]
  cid_XH11B    =addone(itemp)  !exchangeable OH band-soil, [mol m-3]
  cid_XH21B    =addone(itemp)  !exchangeable OH2 band-soil, [mol m-3]

  fid_TRCACO   =addone(itemp)
  fid_TRNAC    =addone(itemp)
  fid_TRMGC    =addone(itemp)
  fid_TRCAC    =addone(itemp)
  fid_TRMGH    =addone(itemp)
  fid_TRCAH    =addone(itemp)
  fid_TRHCO    =addone(itemp)
  fid_TRXHC    =addone(itemp)
  fid_TRC2P    =addone(itemp)
  fid_TRF2P    =addone(itemp)
  fid_TRM1P    =addone(itemp)
  fid_TRC1P    =addone(itemp)
  fid_TRF1P    =addone(itemp)
  fid_TRC0P    =addone(itemp)
  fid_TRH0P    =addone(itemp)
  fid_TRFEPO   =addone(itemp)
  fid_TRALPO   =addone(itemp)
  fid_TRCAPM   =addone(itemp)
  fid_TRCAPD   =addone(itemp)
  fid_TRCAPH   =addone(itemp)
  fid_TRFE1    =addone(itemp)
  fid_TRAL1    =addone(itemp)
  fid_TRFE2    =addone(itemp)
  fid_TRAL2    =addone(itemp)
  fid_TRFE3    =addone(itemp)
  fid_TRAL3    =addone(itemp)
  fid_TRFE4    =addone(itemp)
  fid_TRAL4    =addone(itemp)
  fid_TRMGO    =addone(itemp)
  fid_TRCAO    =addone(itemp)
  fid_TRFEOH   =addone(itemp)
  fid_TRALOH   =addone(itemp)
  fid_TRXFE2   =addone(itemp)
  fid_TRAL     =addone(itemp)
  fid_TRALS   =addone(itemp)
  fid_TRB1P   =addone(itemp)
  fid_TRB2P   =addone(itemp)
  fid_TRBH0   =addone(itemp)
  fid_TRBH1   =addone(itemp)
  fid_TRBH2   =addone(itemp)
  fid_TRCA    =addone(itemp)
  fid_TRCAS   =addone(itemp)
  fid_TRCASO  =addone(itemp)
  fid_TRCO2   =addone(itemp)
  fid_TRCO3   =addone(itemp)
  fid_TRFE    =addone(itemp)
  fid_TRFES   =addone(itemp)

  fid_TRH1P   =addone(itemp)
  fid_TRH2P   =addone(itemp)
  fid_TRH3P   =addone(itemp)
  fid_TRHY    =addone(itemp)
  fid_TRKA    =addone(itemp)
  fid_TRKAS   =addone(itemp)
  fid_TRMG    =addone(itemp)
  fid_TRMGS   =addone(itemp)
  fid_TRN3S   =addone(itemp)
  fid_TRN4S   =addone(itemp)
  fid_TRNA    =addone(itemp)
  fid_TRNAS   =addone(itemp)
  fid_TROH    =addone(itemp)
  fid_TRSO4   =addone(itemp)
  fid_TRX1P   =addone(itemp)
  fid_TRX2P   =addone(itemp)
  fid_TRXAL   =addone(itemp)
  fid_TRXAL2  =addone(itemp)
  fid_TRXCA   =addone(itemp)
  fid_TRXFE   =addone(itemp)
  fid_TRXH0   =addone(itemp)
  fid_TRXH1   =addone(itemp)
  fid_TRXH2   =addone(itemp)
  fid_TRXHY   =addone(itemp)
  fid_TRXKA   =addone(itemp)
  fid_TRXMG   =addone(itemp)
  fid_TRXN4   =addone(itemp)
  fid_TRXNA   =addone(itemp)
  fid_TBCO2   =addone(itemp)
  fid_TBION   =addone(itemp)
  fid_TRH2O   =addone(itemp)

  fid_TRC2B    =addone(itemp)
  fid_TRF2B    =addone(itemp)
  fid_TRM1B    =addone(itemp)
  fid_TRC1B    =addone(itemp)
  fid_TRF1B    =addone(itemp)
  fid_TRC0B    =addone(itemp)
  fid_TRH0B    =addone(itemp)
  fid_TRFEPB   =addone(itemp)
  fid_TRALPB   =addone(itemp)
  fid_TRCPMB   =addone(itemp)
  fid_TRCPDB   =addone(itemp)
  fid_TRCPHB   =addone(itemp)
  fid_TRH1B   =addone(itemp)
  fid_TRH2B   =addone(itemp)
  fid_TRH3B   =addone(itemp)
  fid_TRN3B   =addone(itemp)
  fid_TRN4B   =addone(itemp)
  fid_TRXNB   =addone(itemp)
  nvars=itemp
  end subroutine Init_Salt_geochem


!------------------------------------------------------------------
  subroutine getvarlist(nvars, varl, varlnml,unitl, vartypes, salton)

  use histMod, only : hist_var_str_len,hist_unit_str_len,hist_var_lon_str_len
  use fileUtil, only :  var_flux_type, var_state_type
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(out) :: varl(:)     !variable name
  character(len=hist_var_lon_str_len), intent(out) :: varlnml(:)     !variable name
  character(len=hist_unit_str_len),intent(out) :: unitl(:)
  integer                         ,intent(out) :: vartypes(:)
  logical, intent(in) :: salton

  if(salton)then
    call getvarlist_salt(nvars, varl, varlnml, unitl, vartypes)
  else
    call getvarlist_nosalt(nvars, varl, varlnml, unitl, vartypes)
  endif
  end subroutine getvarlist

!------------------------------------------------------------------
  subroutine getvarlist_salt(nvars, varl, varlnml, unitl, vartypes)

  use histMod, only : hist_var_str_len,hist_unit_str_len,hist_var_lon_str_len
  use fileUtil, only :  var_flux_type, var_state_type
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(out) :: varl(:)     !variable name
  character(len=hist_var_lon_str_len), intent(out) :: varlnml(:)     !variable name
  character(len=hist_unit_str_len),intent(out) :: unitl(:)
  integer                         ,intent(out) :: vartypes(:)


  varl(cid_CO2S)='CO2S';varlnml(cid_CO2S)='aqueous CO2 concentration micropore';
  unitl(cid_CO2S)='gC d-2';vartypes(cid_CO2S)=var_state_type

  varl(cid_CH1P1)='CH1P1';varlnml(cid_CH1P1)='non-band soil aqueous HPO4 content micropore'
  unitl(cid_CH1P1)='mol m-3';vartypes(cid_CH1P1)=var_state_type

  varl(cid_CH2P1)='CH2P1';varlnml(cid_CH2P1)='non-band soil micropore aqueous H2PO4 content'
  unitl(cid_CH2P1)='mol m-3';vartypes(cid_CH2P1)=var_state_type

  varl(cid_CN31) ='CN31';varlnml(cid_CN31)='non-band soil NH3 concentration'
  unitl(cid_CN31)='mol m-3';vartypes(cid_CN31)=var_state_type

  varl(cid_CN41)  ='CN41';varlnml(cid_CN41)='non-band soil NH4 concentration'
  unitl(cid_CN41)='mol m-3';vartypes(cid_CN41)=var_state_type

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

  varl(cid_ZCA2P)='ZCA2P'; varlnml(cid_ZCA2P)='soil aqueous CaH2PO4 content micropore non-band'
  unitl(cid_ZCA2P)='mol d-2';vartypes(cid_ZCA2P)=var_state_type

  varl(cid_ZMG1P)='ZMG1P';varlnml(cid_ZMG1P)='soil aqueous MgHPO4 content micropore non-band'
  unitl(cid_ZMG1P)='mol d-2';vartypes(cid_ZMG1P)=var_state_type

  varl(cid_PALPO1)='PALPO1';varlnml(cid_PALPO1)='non-band soil precipitated AlPO4'
  unitl(cid_PALPO1)='mol m-3';vartypes(cid_PALPO1)=var_state_type

  varl(cid_PCAPD1) ='PCAPD1';varlnml(cid_PCAPD1)='non-band soil precipitated CaHPO4'
  unitl(cid_PCAPD1)='mol m-3';vartypes(cid_PCAPD1)=var_state_type

  varl(cid_PCAPH1) ='PCAPH1';varlnml(cid_PCAPH1)='non-band soil precipitated Ca5(PO4)3OH hydroxyapatite'
  unitl(cid_PCAPH1)='mol m-3';vartypes(cid_PCAPH1)=var_state_type

  varl(cid_PCAPM1) ='PCAPM1';varlnml(cid_PCAPM1)='non-band soil precipitated Ca(H2PO4)2'
  unitl(cid_PCAPM1)='mol m-3';vartypes(cid_PCAPM1)=var_state_type

  varl(cid_PFEPO1) ='PFEPO1';varlnml(cid_PFEPO1)='non-band soil precipitated FePO4'
  unitl(cid_PFEPO1)='mol m-3';vartypes(cid_PFEPO1)=var_state_type

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

  varl(cid_XOH11)='XOH11';varlnml(cid_XOH11)='exchangeable OH  non-band'
  unitl(cid_XOH11)='mol d-2';vartypes(cid_XOH11)=var_state_type

  varl(cid_XALO2)='XALO2';varlnml(cid_XALO2)='exchangeable AlOH2'
  unitl(cid_XALO2)='mol d-2';vartypes(cid_XALO2)=var_state_type

  varl(cid_XFEO2)='XFEO2';varlnml(cid_XFEO2)='exchangeable Fe(OH)2'
  unitl(cid_XFEO2)='mol d-2';vartypes(cid_XFEO2)=var_state_type

  varl(cid_XN41)  ='XN41';varlnml(cid_XN41)='non-band soil adsorbed NH4 concentration'
  unitl(cid_XN41) ='mol m-3';vartypes(cid_XN41)=var_state_type

  varl(cid_XOH01)='XOH01';varlnml(cid_XOH01)='exchangeable OH- non-band'
  unitl(cid_XOH01)='mol d-2';vartypes(cid_XOH01)=var_state_type

  varl(cid_XH1P1) ='XH1P1';varlnml(cid_XH1P1)='non-band soil concentration of adsorbed HPO4'
  unitl(cid_XH1P1)='mol m-3';vartypes(cid_XH1P1)=var_state_type

  varl(cid_XOH21) ='XOH21';varlnml(cid_XOH21)='non-band soil exchangeable site R-OH2'
  unitl(cid_XOH21)='mol m-3';vartypes(cid_XOH21)=var_state_type

  varl(cid_XH2P1) ='XH2P1';varlnml(cid_XH2P1)='non-band soil concentration of adsorbed HPO4'
  unitl(cid_XH2P1)='mol m-3';vartypes(cid_XH2P1)=var_state_type

  varl(cid_CH1PB)='CH1PB';varlnml(cid_CH1PB)='band soil aqueous HPO4 content micropore'
  unitl(cid_CH1PB)='mol m-3';vartypes(cid_CH1PB)=var_state_type

  varl(cid_CH2PB)='CH2PB';varlnml(cid_CH2PB)='band soil micropore aqueous H2PO4 content'
  unitl(cid_CH2PB)='mol m-3';vartypes(cid_CH2PB)=var_state_type

  varl(cid_CN3B) ='CN3B';varlnml(cid_CN3B)='band soil NH3 concentration'
  unitl(cid_CN3B)='mol m-3';vartypes(cid_CN3B)=var_state_type

  varl(cid_CN4B)  ='CN4B';varlnml(cid_CN4B)='band soil NH4 concentration'
  unitl(cid_CN4B) ='mol m-3';vartypes(cid_CN4B)=var_state_type

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

  varl(cid_ZCA2PB) ='ZCA2PB';varlnml(cid_ZCA2PB)='band soil aqueous CaH2PO4 content micropore'
  unitl(cid_ZCA2PB)= 'mol d-2';vartypes(cid_ZCA2PB)=var_state_type

  varl(cid_ZMG1PB)='ZMG1PB';varlnml(cid_ZMG1PB)='band soil aqueous MgHPO4 content micropore'
  unitl(cid_ZMG1PB)='mol d-2';vartypes(cid_ZMG1PB)=var_state_type

  varl(cid_PALPOB) ='PALPOB';varlnml(cid_PALPOB)='band soil precipitated AlPO4'
  unitl(cid_PALPOB)='mol m-3';vartypes(cid_PALPOB)=var_state_type

  varl(cid_PCAPDB) ='PCAPDB';varlnml(cid_PCAPDB)='band soil precipitated CaHPO4'
  unitl(cid_PCAPDB)='mol m-3';vartypes(cid_PCAPDB)=var_state_type

  varl(cid_PCAPHB) ='PCAPHB';varlnml(cid_PCAPHB)='band soil precipitated Ca5(PO4)3OH hydroxyapatite'
  unitl(cid_PCAPHB)='mol m-3';vartypes(cid_PCAPHB)=var_state_type

  varl(cid_PCAPMB) ='PCAPMB';varlnml(cid_PCAPMB)='band soil precipitated CaH2PO4'
  unitl(cid_PCAPMB)='mol m-3';vartypes(cid_PCAPMB)=var_state_type

  varl(cid_PFEPOB) ='PFEPOB';varlnml(cid_PFEPOB)='band soil precipitated FePO4'
  unitl(cid_PFEPOB)='mol m-3';vartypes(cid_PFEPOB)=var_state_type

  varl(cid_XN4B)  ='XN4B';varlnml(cid_XN4B)='band soil adsorbed NH4 concentration'
  unitl(cid_XN4B) ='mol m-3';vartypes(cid_XN4B)=var_state_type

  varl(cid_XH01B) ='XH01B';varlnml(cid_XH01B)='band soil exchangeable R-OH'
  unitl(cid_XH01B)='mol d-2';vartypes(cid_XH01B)=var_state_type

  varl(cid_X1P1B) ='X1P1B';varlnml(cid_X1P1B)='band soil exchangeable HPO4 concentration'
  unitl(cid_X1P1B)='mol m-3';vartypes(cid_X1P1B)=var_state_type

  varl(cid_X2P1B) ='X2P1B';varlnml(cid_X2P1B)='band soil exchangeable H2PO4 concentration'
  unitl(cid_X2P1B)='mol m-3';vartypes(cid_X2P1B)=var_state_type

  varl(cid_XH11B) ='XH11B';varlnml(cid_XH11B)='band soil concentration of adsorbed HPO4'
  unitl(cid_XH11B)='mol m-3';vartypes(cid_XH11B)=var_state_type

  varl(cid_XH21B) ='XH21B';varlnml(cid_XH21B)='band soil exchangeable site R-OH2'
  unitl(cid_XH21B)='mol m-3';vartypes(cid_XH21B)=var_state_type

  varl(fid_TRCACO)='TRCACO';varlnml(fid_TRCACO)='total precipitated CaCO3 transformation'
  unitl(fid_TRCACO)='mol d-2 h-1';vartypes(fid_TRCACO)=var_flux_type

  varl(fid_TRNAC)='TRNAC';varlnml(fid_TRNAC)='total solute NaCO3 transformation'
  unitl(fid_TRNAC)='mol d-2 h-1';vartypes(fid_TRNAC)=var_flux_type

  varl(fid_TRMGC)='TRMGC';varlnml(fid_TRMGC)='total solute MgCO3 transformation'
  unitl(fid_TRMGC)='mol d-2 h-1'; vartypes(fid_TRMGC)=var_flux_type

  varl(fid_TRCAC)='TRCAC';varlnml(fid_TRCAC)='total solute CaCO3 transformation'
  unitl(fid_TRCAC)='mol d-2 h-1';vartypes(fid_TRCAC)=var_flux_type

  varl(fid_TRMGH)='TRMGH';varlnml(fid_TRMGH)='total solute MgHCO3 transformation'
  unitl(fid_TRMGH)='mol d-2 h-1';vartypes(fid_TRMGH)=var_flux_type

  varl(fid_TRCAH)='TRCAH';varlnml(fid_TRCAH)='total solute CaHCO3 transformation'
  unitl(fid_TRCAH)='mol d-2 h-1';vartypes(fid_TRCAH)=var_flux_type

  varl(fid_TRHCO)='TRHCO';varlnml(fid_TRHCO)='total solute HCO3 transformation'
  unitl(fid_TRHCO)='mol d-2 h-1';vartypes(fid_TRHCO)=var_flux_type

  varl(fid_TRXHC)='TRXHC';varlnml(fid_TRXHC)='total adsorbed COOH transformation'
  unitl(fid_TRXHC)='mol d-2 h-1';vartypes(fid_TRXHC)=var_flux_type

  varl(fid_TRC2P)='TRC2P';varlnml(fid_TRC2P)='non-band soil total solute CaH2PO4 transformation'
  unitl(fid_TRC2P)='mol d-2 h-1';vartypes(fid_TRC2P)=var_flux_type

  varl(fid_TRF2P)='TRF2P';varlnml(fid_TRF2P)='non-band soil total solute FeH2PO4 transformation'
  unitl(fid_TRF2P)='mol d-2 h-1';vartypes(fid_TRF2P)=var_flux_type

  varl(fid_TRM1P)='TRM1P';varlnml(fid_TRM1P)='non-band soil total solute MgHPO4 transformation'
  unitl(fid_TRM1P)='mol d-2 h-1';vartypes(fid_TRM1P)=var_flux_type

  varl(fid_TRC1P)='TRC1P';varlnml(fid_TRC1P)='non-band soil total solute CaHPO4 transformation'
  unitl(fid_TRC1P)='mol d-2 h-1'; vartypes(fid_TRC1P)=var_flux_type

  varl(fid_TRF1P)='TRF1P';varlnml(fid_TRF1P)='non-band soil total solute FeHPO4 transformation'
  unitl(fid_TRF1P)='mol d-2 h-1';vartypes(fid_TRF1P)=var_flux_type

  varl(fid_TRC0P)='TRC0P';varlnml(fid_TRC0P)='non-band soil total solute CaPO4 transformation'
  unitl(fid_TRC0P)='mol d-2 h-1'; vartypes(fid_TRC0P)=var_flux_type

  varl(fid_TRH0P)='TRH0P';varlnml(fid_TRH0P)='non-band soil total solute PO4 transformation'
  unitl(fid_TRH0P)='mol d-2 h-1';vartypes(fid_TRH0P)=var_flux_type

  varl(fid_TRFEPO)= 'TRFEPO';varlnml(fid_TRFEPO)='non-band soil total precipitated FePO4 transformation'
  unitl(fid_TRFEPO)= 'mol d-2 h-1';vartypes(fid_TRFEPO)=var_flux_type

  varl(fid_TRALPO)= 'TRALPO';varlnml(fid_TRALPO)='non-band total precipitated AlPO4 transformation'
  unitl(fid_TRALPO)= 'mol d-2 h-1';vartypes(fid_TRALPO)=var_flux_type

  varl(fid_TRCAPM)= 'TRCAPM';varlnml(fid_TRCAPM)='non-band soil total precipitated apatite transformation'
  unitl(fid_TRCAPM)= 'mol d-2 h-1';vartypes(fid_TRCAPM)=var_flux_type

  varl(fid_TRCAPD)= 'TRCAPD';varlnml(fid_TRCAPD)='non-band soil total precipitated CaHPO4 transformation'
  unitl(fid_TRCAPD)= 'mol d-2 h-1';vartypes(fid_TRCAPD)=var_flux_type

  varl(fid_TRCAPH)= 'TRCAPH';varlnml(fid_TRCAPH)='non-band soil total precipitated CaH2PO4 transformation'
  unitl(fid_TRCAPH)= 'mol d-2 h-1';vartypes(fid_TRCAPH)=var_flux_type

  varl(fid_TRFE1)='TRFE1';varlnml(fid_TRFE1)='total solute FeOH transformation'
  unitl(fid_TRFE1)='mol d-2 h-1';    vartypes(fid_TRFE1)=var_flux_type

  varl(fid_TRAL1)='TRAL1';varlnml(fid_TRAL1)='total solute AlOH transformation'
  unitl(fid_TRAL1)='mol d-2 h-1'; vartypes(fid_TRAL1)=var_flux_type

  varl(fid_TRFE2)='TRFE2';varlnml(fid_TRFE2)='total solute FeOH2 transformation'
  unitl(fid_TRFE2)='mol d-2 h-1';    vartypes(fid_TRFE2)=var_flux_type

  varl(fid_TRAL2)='TRAL2';varlnml(fid_TRAL2)='total solute AlOH2 transformation'
  unitl(fid_TRAL2)='mol d-2 h-1'; vartypes(fid_TRAL2)=var_flux_type

  varl(fid_TRFE3)='TRFE3';varlnml(fid_TRFE3)='total solute FeOH3 transformation'
  unitl(fid_TRFE3)='mol d-2 h-1';    vartypes(fid_TRFE3)=var_flux_type

  varl(fid_TRAL3)='TRAL3';varlnml(fid_TRAL3)='total solute AlOH3 transformation'
  unitl(fid_TRAL3)='mol d-2 h-1'; vartypes(fid_TRAL3)=var_flux_type

  varl(fid_TRFE4)='TRFE4';varlnml(fid_TRFE4)='total solute FeOH4 transformation'
  unitl(fid_TRFE4)='mol d-2 h-1';    vartypes(fid_TRFE4)=var_flux_type

  varl(fid_TRAL4)='TRAL4';varlnml(fid_TRAL4)='total solute AlOH4 transformation'
  unitl(fid_TRAL4)='mol d-2 h-1'; vartypes(fid_TRAL4)=var_flux_type


  varl(fid_TRMGO)='TRMGO';varlnml(fid_TRMGO)='total solute MgOH transformation'
  unitl(fid_TRMGO)='mol d-2 h-1';vartypes(fid_TRMGO)=var_flux_type

  varl(fid_TRCAO)='TRCAO';varlnml(fid_TRCAO)='total solute CaOH transformation'
  unitl(fid_TRCAO)='mol d-2 h-1';vartypes(fid_TRCAO)=var_flux_type

  varl(fid_TRFEOH) ='TRFEOH';varlnml(fid_TRFEOH)='total precipitated FeOH3 transformation'
  unitl(fid_TRFEOH)='mol d-2 h-1'; vartypes(fid_TRFEOH)=var_flux_type

  varl(fid_TRALOH)='TRALOH';varlnml(fid_TRALOH)='total precipitated AlOH3 transformation'
  unitl(fid_TRALOH)='mol d-2 h-1'; vartypes(fid_TRALOH)=var_flux_type

  varl(fid_TRXFE2)='TRXFE2';varlnml(fid_TRXFE2)='total FeOH2 adsorption'
  unitl(fid_TRXFE2)='mol d-2 h-1'; vartypes(fid_TRXFE2)=var_flux_type

  varl(fid_TRAL)  = 'TRAL';varlnml(fid_TRAL)='non-band soil total solute Al transformation'
  unitl(fid_TRAL) = 'mol d-2 h-1';vartypes(fid_TRAL)=var_flux_type

  varl(fid_TRALS)='TRALS'; varlnml(fid_TRALS)='total solute AlSO4 transformation'
  unitl(fid_TRALS)='mol d-2 h-1'; vartypes(fid_TRALS)=var_flux_type

  varl(fid_TRB1P)='TRB1P';varlnml(fid_TRB1P)='band soil total adsorbed HPO4 transformation'
  unitl(fid_TRB1P)='mol d-2 h-1';vartypes(fid_TRB1P)=var_flux_type

  varl(fid_TRB2P)='TRB2P';varlnml(fid_TRB2P)='band soil total adsorbed H2PO4 transformation'
  unitl(fid_TRB2P)='mol d-2 h-1';vartypes(fid_TRB2P)=var_flux_type

  varl(fid_TRBH0)='TRBH0';varlnml(fid_TRBH0)='band soil total adsorbed OH- transformation'
  unitl(fid_TRBH0)='mol d-2 h-1'; vartypes(fid_TRBH0)=var_flux_type

  varl(fid_TRBH1)='TRBH1';varlnml(fid_TRBH1)='band soil total adsorbed OH transformation'
  unitl(fid_TRBH1)='mol d-2 h-1'; vartypes(fid_TRBH1)=var_flux_type

  varl(fid_TRBH2)='TRBH2';varlnml(fid_TRBH2)='band soil total adsorbed OH2 transformation'
  unitl(fid_TRBH2)='mol d-2 h-1'; vartypes(fid_TRBH2)=var_flux_type

  varl(fid_TRCA) ='TRCA';varlnml(fid_TRCA)='total solute Ca transformation'
  unitl(fid_TRCA)='mol d-2 h-1';vartypes(fid_TRCA)=var_flux_type

  varl(fid_TRCAS)='TRCAS';varlnml(fid_TRCAS)='total solute CaSO4 transformation'
  unitl(fid_TRCAS)='mol d-2 h-1';vartypes(fid_TRCAS)=var_flux_type

  varl(fid_TRCASO)='TRCASO';varlnml(fid_TRCASO)='total precipitated CaSO4 transformation'
  unitl(fid_TRCASO)='mol d-2 h-1';vartypes(fid_TRCASO)=var_flux_type

  varl(fid_TRCO2)='TRCO2';varlnml(fid_TRCO2)='total solute CO2 transformation'
  unitl(fid_TRCO2)='mol d-2 h-1';vartypes(fid_TRCO2)=var_flux_type

  varl(fid_TRCO3)='TRCO3';varlnml(fid_TRCO3)='total solute CO3 transformation'
  unitl(fid_TRCO3)='mol d-2 h-1';vartypes(fid_TRCO3)=var_flux_type

  varl(fid_TRFE) ='TRFE';varlnml(fid_TRFE)='total solute Fe transformation'
  unitl(fid_TRFE)='mol d-2 h-1';vartypes(fid_TRFE)= var_flux_type

  varl(fid_TRFES)='TRFES';varlnml(fid_TRFES)='total solute FeSO4 transformation'
  unitl(fid_TRFES)='mol d-2 h-1';vartypes(fid_TRFES)=var_flux_type

  varl(fid_TRH1P) = 'TRH1P';varlnml(fid_TRH1P)='non-band soil total solute HPO4 transformation'
  unitl(fid_TRH1P) = 'mol d-2 h-1';vartypes(fid_TRH1P)=var_flux_type

  varl(fid_TRH2P) = 'TRH2P';varlnml(fid_TRH2P)='non-band soil total solute H2PO4 transformation'
  unitl(fid_TRH2P)= 'mol d-2 h-1';vartypes(fid_TRH2P)=var_flux_type

  varl(fid_TRH3P) = 'TRH3P';varlnml(fid_TRH3P)='non-band soil total solute H3PO4 transformation'
  unitl(fid_TRH3P)= 'mol d-2 h-1';vartypes(fid_TRH3P)=var_flux_type

  varl(fid_TRHY)='TRHY';varlnml(fid_TRHY)='total solute H transformation'
  unitl(fid_TRHY)='mol d-2 h-1'; vartypes(fid_TRHY)=var_flux_type

  varl(fid_TRKA)='TRKA';varlnml(fid_TRKA)='total solute K transformation'
  unitl(fid_TRKA)='mol d-2 h-1'; vartypes(fid_TRKA)=var_flux_type

  varl(fid_TRKAS)='TRKAS';varlnml(fid_TRKAS)='total solute KSO4 transformation'
  unitl(fid_TRKAS)='mol d-2 h-1'; vartypes(fid_TRKAS)=var_flux_type

  varl(fid_TRMG)='TRMG';varlnml(fid_TRMG)='total solute Mg transformation'
  unitl(fid_TRMG)='mol d-2 h-1'; vartypes(fid_TRMG)=var_flux_type

  varl(fid_TRMGS)='TRMGS';varlnml(fid_TRMGS)='total solute MgSO4 transformation'
  unitl(fid_TRMGS)='mol d-2 h-1';vartypes(fid_TRMGS)=var_flux_type

  varl(fid_TRN3S) = 'TRN3S';varlnml(fid_TRN3S)='non-band total solute NH3 transformation'
  unitl(fid_TRN3S)= 'mol d-2 h-1';vartypes(fid_TRN3S)=var_flux_type

  varl(fid_TRN4S) = 'TRN4S';varlnml(fid_TRN4S)='non-band soil total solute NH4 transformation'
  unitl(fid_TRN4S)= 'mol d-2 h-1';vartypes(fid_TRN4S)=var_flux_type

  varl(fid_TRNA) ='TRNA';varlnml(fid_TRNA)='total solute Na transformation'
  unitl(fid_TRNA)='mol d-2 h-1';vartypes(fid_TRNA)=var_flux_type

  varl(fid_TRNAS)='TRNAS';varlnml(fid_TRNAS)='total solute NaSO4 transformation'
  unitl(fid_TRNAS)='mol d-2 h-1';vartypes(fid_TRNAS)=var_flux_type

  varl(fid_TROH)='TROH';varlnml(fid_TROH)='total solute OH transformation'
  unitl(fid_TROH)='mol d-2 h-1';vartypes(fid_TROH)=var_flux_type

  varl(fid_TRSO4)='TRSO4';varlnml(fid_TRSO4)='total solute SO4 transformation'
  unitl(fid_TRSO4)='mol d-2 h-1';vartypes(fid_TRSO4)=var_flux_type

  varl(fid_TRX1P)='TRX1P';varlnml(fid_TRX1P)='non-band soil total adsorbed HPO4 transformation'
  unitl(fid_TRX1P)='mol d-2 h-1';vartypes(fid_TRX1P)=var_flux_type

  varl(fid_TRX2P)='TRX2P';varlnml(fid_TRX2P)='non-band soil total adsorbed H2PO4 transformation'
  unitl(fid_TRX2P)='mol d-2 h-1';vartypes(fid_TRX2P)=var_flux_type

  varl(fid_TRXAL)='TRXAL';varlnml(fid_TRXAL)='total adsorbed Al transformation'
  unitl(fid_TRXAL)='mol d-2 h-1';vartypes(fid_TRXAL)=var_flux_type

  varl(fid_TRXAL2)='TRXAL2';varlnml(fid_TRXAL2)='total adsorbed AlOH2 transformation'
  unitl(fid_TRXAL2)='mol d-2 h-1';vartypes(fid_TRXAL2)=var_flux_type

  varl(fid_TRXCA)='TRXCA';varlnml(fid_TRXCA)='total adsorbed Ca transformation'
  unitl(fid_TRXCA)='mol d-2 h-1';vartypes(fid_TRXCA)=var_flux_type

  varl(fid_TRXFE)='TRXFE';varlnml(fid_TRXFE)='total Fe adsorption'
  unitl(fid_TRXFE)='mol d-2 h-1';vartypes(fid_TRXFE)=var_flux_type

  varl(fid_TRXH0)='TRXH0';varlnml(fid_TRXH0)='non-band soil total adsorbed OH- transformation'
  unitl(fid_TRXH0)='mol d-2 h-1';vartypes(fid_TRXH0)=var_flux_type

  varl(fid_TRXH1)='TRXH1';varlnml(fid_TRXH1)='non-band total adsorbed OH transformation'
  unitl(fid_TRXH1)='mol d-2 h-1';vartypes(fid_TRXH1)=var_flux_type

  varl(fid_TRXH2)='TRXH2';varlnml(fid_TRXH2)='non-band soil total adsorbed OH2 transformation'
  unitl(fid_TRXH2)='mol d-2 h-1';vartypes(fid_TRXH2)=var_flux_type

  varl(fid_TRXHY)='TRXHY';varlnml(fid_TRXHY)='total adsorbed H transformation'
  unitl(fid_TRXHY)='mol d-2 h-1';vartypes(fid_TRXHY)=var_flux_type

  varl(fid_TRXKA)='TRXKA';varlnml(fid_TRXKA)='total adsorbed K transformation'
  unitl(fid_TRXKA)='mol d-2 h-1';vartypes(fid_TRXKA)=var_flux_type

  varl(fid_TRXMG)='TRXMG';varlnml(fid_TRXMG)='total adsorbed Mg transformation'
  unitl(fid_TRXMG)='mol d-2 h-1'; vartypes(fid_TRXMG)=var_flux_type

  varl(fid_TRXN4) = 'TRXN4';varlnml(fid_TRXN4)='non-band soil total adsorbed NH4 transformation'
  unitl(fid_TRXN4)= 'mol d-2 h-1';vartypes(fid_TRXN4)=var_flux_type

  varl(fid_TRXNA)='TRXNA';varlnml(fid_TRXNA)='total adsorbed Na transformation'
  unitl(fid_TRXNA)='mol d-2 h-1';vartypes(fid_TRXNA)=var_flux_type

  varl(fid_TBCO2)='TBCO2';varlnml(fid_TBCO2)='total solute CO2 transformation'
  unitl(fid_TBCO2)='mol d-2 h-1';vartypes(fid_TBCO2)=var_flux_type

  varl(fid_TBION)='TBION';varlnml(fid_TBION)='total solute ion transformation'
  unitl(fid_TBION)='mol d-2 h-1';vartypes(fid_TBION)=var_flux_type

  varl(fid_TRH2O)='TRH2O';varlnml(fid_TRH2O)='total solute H2O transformation'
  unitl(fid_TRH2O)='mol d-2 h-1';vartypes(fid_TRH2O)=var_flux_type

  varl(fid_TRC2B)='TRC2B';varlnml(fid_TRC2B)='band soil total solute CaH2PO4 transformation'
  unitl(fid_TRC2B)='mol d-2 h-1';vartypes(fid_TRC2B)=var_flux_type

  varl(fid_TRF2B)='TRF2B';varlnml(fid_TRF2B)='band soil total solute FeH2PO4 transformation'
  unitl(fid_TRF2B)='mol d-2 h-1';vartypes(fid_TRF2B)=var_flux_type

  varl(fid_TRM1B)='TRM1B';varlnml(fid_TRM1B)='band soil total solute MgHPO4 transformation'
  unitl(fid_TRM1B)='mol d-2 h-1';vartypes(fid_TRM1B)=var_flux_type

  varl(fid_TRC1B)='TRC1B';varlnml(fid_TRC1B)='band soil total solute CaHPO4 transformation'
  unitl(fid_TRC1B)='mol d-2 h-1';vartypes(fid_TRC1B)=var_flux_type

  varl(fid_TRF1B)='TRF1B';varlnml(fid_TRF1B)='band soil total solute FeHPO4 transformation'
  unitl(fid_TRF1B)='mol d-2 h-1';vartypes(fid_TRF1B)=var_flux_type

  varl(fid_TRC0B)='TRC0B';varlnml(fid_TRC0B)='band soil total solute CaPO4 transformation'
  unitl(fid_TRC0B)='mol d-2 h-1';vartypes(fid_TRC0B)=var_flux_type

  varl(fid_TRH0B)='TRH0B';varlnml(fid_TRH0B)='band soil total solute PO4 transformation'
  unitl(fid_TRH0B)='mol d-2 h-1';vartypes(fid_TRH0B)=var_flux_type

  varl(fid_TRFEPB)= 'TRFEPB';varlnml(fid_TRFEPB)='band soil total precipitated FePO4 transformation'
  unitl(fid_TRFEPB)= 'mol d-2 h-1';vartypes(fid_TRFEPB)=var_flux_type

  varl(fid_TRALPB)= 'TRALPB';varlnml(fid_TRALPB)='band soil total precipitated AlPO4 transformation'
  unitl(fid_TRALPB)= 'mol d-2 h-1';vartypes(fid_TRALPB)=var_flux_type

  varl(fid_TRCPMB)= 'TRCPMB';varlnml(fid_TRCPMB)='band soil total precipitated apatite transformation'
  unitl(fid_TRCPMB)= 'mol d-2 h-1';vartypes(fid_TRCPMB)=var_flux_type

  varl(fid_TRCPDB)= 'TRCPDB';varlnml(fid_TRCPDB)='band soil total precipitated CaHPO4 transformation'
  unitl(fid_TRCPDB)= 'mol d-2 h-1';vartypes(fid_TRCPDB)=var_flux_type

  varl(fid_TRCPHB)= 'TRCPHB';varlnml(fid_TRCPHB)='band soil total precipitated CaH2PO4 transformation'
  unitl(fid_TRCPHB)= 'mol d-2 h-1';vartypes(fid_TRCPHB)=var_flux_type

  varl(fid_TRH1B) = 'TRH1B';varlnml(fid_TRH1B)='band soil total solute HPO4 transformation'
  unitl(fid_TRH1B)= 'mol d-2 h-1';vartypes(fid_TRH1B)=var_flux_type

  varl(fid_TRH2B) = 'TRH2B';varlnml(fid_TRH2B)='band soil total solute H2PO4 transformation'
  unitl(fid_TRH2B)= 'mol d-2 h-1';vartypes(fid_TRH2B)=var_flux_type

  varl(fid_TRH3B) = 'TRH3B';varlnml(fid_TRH3B)='band soil total solute H3PO4 transformation'
  unitl(fid_TRH3B)= 'mol d-2 h-1';vartypes(fid_TRH3B)=var_flux_type

  varl(fid_TRN3B) = 'TRN3B';varlnml(fid_TRN3B)='band soil total solute NH3 transformation'
  unitl(fid_TRN3B)= 'mol d-2 h-1';vartypes(fid_TRN3B)=var_flux_type

  varl(fid_TRN4B) = 'TRN4B';varlnml(fid_TRN4B)='band soil total solute NH4 transformation'
  unitl(fid_TRN4B)= 'mol d-2 h-1';vartypes(fid_TRN4B)=var_flux_type

  varl(fid_TRXNB) = 'TRXNB';varlnml(fid_TRXNB)='band soil total adsorbed NH4 transformation'
  unitl(fid_TRXNB)= 'mol d-2 h-1';vartypes(fid_TRXNB)=var_flux_type
  end subroutine getvarlist_salt

!------------------------------------------------------------------
  subroutine getvarlist_nosalt(nvars, varl, varlnml, unitl, vartypes)

  use histMod, only : hist_var_str_len,hist_unit_str_len,hist_var_lon_str_len
  use fileUtil, only :  var_flux_type, var_state_type
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(out) :: varl(:)     !variable name
  character(len=hist_var_lon_str_len), intent(out) :: varlnml(:)     !variable name
  character(len=hist_unit_str_len),intent(out) :: unitl(:)
  integer                         ,intent(out) :: vartypes(:)

  varl(cid_ZMG) ='ZMG';varlnml(cid_ZMG)='soil aqueous Mg content micropore'
  unitl(cid_ZMG)='mol d-2';vartypes(cid_ZMG)=var_state_type

  varl(cid_ZNA) ='ZNA';varlnml(cid_ZNA)='soil aqueous Na content micropore'
  unitl(cid_ZNA) ='mol d-2';vartypes(cid_ZNA)=var_state_type

  varl(cid_ZKA) ='ZKA';varlnml(cid_ZKA)='soil aqueous K content micropore'
  unitl(cid_ZKA) ='mol d-2';vartypes(cid_ZKA)=var_state_type

  varl(cid_CO2S)='CO2S';varlnml(cid_CO2S)='aqueous CO2 concentration micropore';
  unitl(cid_CO2S)='gC d-2';vartypes(cid_CO2S)=var_state_type

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

  varl(cid_XH11B) ='XH11B';varlnml(cid_XH11B)='band soil concentration of adsorbed HPO4'
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
! ----------------------------------------------------------------------

  subroutine initmodel(nvars, ystatesfl, salton, err_status)

  implicit none
  integer, intent(in) :: nvars
  logical, intent(in) :: salton
  real(r8), intent(inout) :: ystatesfl(nvars)
  type(model_status_type), intent(out) :: err_status

  call err_status%reset()
  if (salton)then
    call initmodel_salt(nvars, ystatesfl,err_status)
  else
    call initmodel_nosalt(nvars, ystatesfl,err_status)
  endif

  end subroutine initmodel


! ----------------------------------------------------------------------

  subroutine initmodel_salt(nvars, ystatesfl,err_status)

  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(inout) :: ystatesfl(nvars)
  type(model_status_type), intent(out) :: err_status


  end subroutine initmodel_salt

! ----------------------------------------------------------------------

  subroutine initmodel_nosalt(nvars, ystatesfl,err_status)
!
! DESCRIPTION:
!
  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(inout) :: ystatesfl(nvars)
  type(model_status_type), intent(out) :: err_status


  end subroutine initmodel_nosalt
end module AquachemMod
