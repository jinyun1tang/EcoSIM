module AquachemMod
  use data_kind_mod     , only : r8 => SHR_KIND_R8
  use MiscUtilMod, only : addone
implicit none
  public
  character(len=*), private, parameter :: mod_filename = __FILE__

  integer :: cid_CO2S     !aqueous CO2  micropore	[g d-2]
  integer :: cid_CH1P1    !soil aqueous HPO4 content micropore non-band, [mol m-3]
  integer :: cid_CH1PB    !soil aqueous HPO4 content micropore band, [mol m-3]
  integer :: cid_CH2P1    !soil aqueous H2PO4 content micropore non-band, [mol m-3]
  integer :: cid_CH2PB    !soil aqueous H2PO4 content micropore  band, [mol m-3]
  integer :: cid_CN31     !soil NH3 concentration in non-band soil, [mol m-3]
  integer :: cid_CN3B     !soil NH3 concentration in band soil, [mol m-3]
  integer :: cid_CN41     !soil NH4 concentration in non-band soil, [mol m-3]
  integer :: cid_CN4B     !soil NH4 concentration in band soil, [mol m-3]
  integer :: cid_ZNO3S    !NO3 mass non-band micropore, [g d-2]
  integer :: cid_ZNO3B    !NO3 mass band micropore, [g d-2]
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
  integer :: cid_H0POB    !soil aqueous PO4 content micropore band, [mol d-2]
  integer :: cid_H3POB    !soil aqueous H3PO4 content micropore band, [mol d-2]
  integer :: cid_ZFE1PB   !soil aqueous FeHPO4 content micropore band, [mol d-2]
  integer :: cid_ZFE2PB   !soil aqueous FeH2PO4 content micropore band, [mol d-2
  integer :: cid_ZCA0PB   !soil aqueous CaPO4 content micropore band, [mol d-2]
  integer :: cid_ZCA1PB   !soil aqueous CaHPO4 content micropore band, [mol d-2]
  integer :: cid_ZCA2PB   !soil aqueous CaH2PO4 content micropore band, [mol d-2]
  integer :: cid_ZMG1PB   !soil aqueous MgHPO4 content micropore band, [mol d-2]
  integer :: cid_PALPO1   !precipitated AlPO4 non-band, [mol m-3]
  integer :: cid_PALPOB   !precipitated AlPO4 band soil, [mol m-3]
  integer :: cid_PCAPD1   !precipitated CaHPO4 non-band soil, [mol m-3]
  integer :: cid_PCAPDB   !precipitated CaHPO4 band soil, [mol m-3]
  integer :: cid_PCAPH1   !precipitated Ca5(PO4)3OH hydroxyapatite non-band soil, [mol m-3]
  integer :: cid_PCAPHB   !precipitated Ca5(PO4)3OH hydroxyapatite band soil, [mol m-3]
  integer :: cid_PCAPM1   !precipitated Ca(H2PO4)2 non-band soil, [mol m-3]
  integer :: cid_PCAPMB   !precipitated CaH2PO4 band soil, [mol m-3]
  integer :: cid_PFEPO1   !precipitated FePO4 non-band soil, [mol m-3]
  integer :: cid_PFEPOB   !precipitated FePO4 band soil, [mol m-3]
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
  integer :: cid_XN4B     !exchangeable NH4 band soil, [mol d-2]
  integer :: cid_XH01B    !exchangeable OH- band, [mol d-2]
  integer :: cid_XOH01    !exchangeable OH- non-band, [mol d-2]
  integer :: cid_X1P1B    !exchangeable HPO4 concentration band-soil, [mol m-3]
  integer :: cid_X2P1B    !exchangeable H2PO4 concentration band-soil, [mol m-3]
  integer :: cid_XH11B    !exchangeable OH band-soil, [mol m-3]
  integer :: cid_XH1P1    !exchangeable HPO4  non-band, [mol m-3]
  integer :: cid_XH21B    !exchangeable OH2 band-soil, [mol m-3]
  integer :: cid_XOH21    !exchangeable OH2  non-band soil, [mol m-3]
  integer :: cid_XH2P1    !exchangeable H2PO4  non-band soil, [mol m-3]
  integer :: nvars
contains

  suborutine Init_NoSalt_geochem

  implicit none
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
  nvars = itemp
  end suborutine Init_NoSalt_geochem

  subroutine Init_Salt_geochem
  implicit none

  itemp=0
  cid_CO2S     =addone(itemp)  !aqueous CO2  micropore	[g d-2]
  cid_CH1P1    =addone(itemp)  !soil aqueous HPO4 content micropore non-band, [mol m-3]
  cid_CH1PB    =addone(itemp)  !soil aqueous HPO4 content micropore band, [mol m-3]
  cid_CH2P1    =addone(itemp)  !soil aqueous H2PO4 content micropore non-band, [mol m-3]
  cid_CH2PB    =addone(itemp)  !soil aqueous H2PO4 content micropore  band, [mol m-3]
  cid_CN31     =addone(itemp)  !soil NH3 concentration in non-band soil, [mol m-3]
  cid_CN3B     !soil NH3 concentration in band soil, [mol m-3]
  cid_CN41     !soil NH4 concentration in non-band soil, [mol m-3]
  cid_CN4B     !soil NH4 concentration in band soil, [mol m-3]
  cid_ZNO3S    !NO3 mass non-band micropore, [g d-2]
  cid_ZNO3B    !NO3 mass band micropore, [g d-2]
  cid_ZHY      !soil aqueous H content micropore, [mol d-2]
  cid_ZOH      !soil aqueous OH content micropore, [mol d-2]
  cid_ZAL      !soil aqueous Al content micropore, [mol d-2]
  cid_ZFE      !soil aqueous Fe content micropore, [mol d-2]
  cid_ZCA      !soil aqueous Ca content micropore, [mol d-2]
  cid_ZMG      !soil aqueous Mg content micropore, [mol d-2]
  cid_ZNA      !soil aqueous Na content micropore, [mol d-2]
  cid_ZKA      !soil aqueous K content micropore, [mol d-2]
  cid_ZSO4     !soil aqueous SO4 content micropore, [mol d-2]
  cid_ZCL      !soil aqueous Cl content micropore, [mol d-2]
  cid_ZCO3     !soil aqueous CO3 content micropore, [mol d-2]
  cid_ZHCO3    !soil aqueous HCO3 content micropore, [mol d-2]
  cid_ZALOH1   !soil aqueous AlOH content micropore, [mol d-2]
  cid_ZALOH2   !soil aqueous AlOH2 content micropore, [mol d-2]
  cid_ZALOH3   !soil aqueous AlOH3 content micropore, [mol d-2]
  cid_ZALOH4   !soil aqueous AlOH4 content micropore, [mol d-2]
  cid_ZALS     !soil aqueous AlSO4 content micropore, [mol d-2]
  cid_ZFEOH1   !soil aqueous FeOH content micropore, [mol d-2]
  cid_ZFEOH2   !soil aqueous FeOH2 content micropore, [mol d-2]
  cid_ZFEOH3   !soil aqueous FeOH3 content micropore, [mol d-2]
  cid_ZFEOH4   !soil aqueous FeOH4 content micropore, [mol d-2]
  cid_ZFES     !soil aqueous FeSO4 content micropore, [mol d-2]
  cid_ZCAO     !soil aqueous CaOH2 content micropore, [mol d-2]
  cid_ZCAC     !soil aqueous CACO3 content micropore, [mol d-2]
  cid_ZCAH     !soil aqueous CaHCO3 content micropore, [mol d-2]
  cid_ZCAS     !soil aqueous CaSO4 content micropore, [mol d-2]
  cid_ZMGO     !soil aqueous MgOH content micropore, [mol d-2]
  cid_ZMGC     !soil aqueous MgCO3 content micropore, [mol d-2]
  cid_ZMGH     !soil aqueous MgHCO3 content micropore, [mol d-2]
  cid_ZMGS     !soil aqueous MgSO4 content micropore, [mol d-2]
  cid_ZNAC     !soil aqueous NaCO3 content micropore, [mol d-2]
  cid_ZNAS     !soil aqueous NaSO4 content micropore, [mol d-2]
  cid_ZKAS     !soil aqueous KSO4 content micropore, [mol d-2]
  cid_H0PO4    !soil aqueous PO4 content micropore non-band, [mol d-2]
  cid_H3PO4    !soil aqueous H3PO4 content micropore non-band, [mol d-2]
  cid_ZFE1P    !soil aqueous FeHPO4 content micropore non-band, [mol d-2]
  cid_ZFE2P    !soil aqueous FeH2PO4 content micropore non-band, [mol d-2]
  cid_ZCA0P    !soil aqueous CaPO4 content micropore non-band, [mol d-2]
  cid_ZCA1P    !soil aqueous CaHPO4 content micropore non-band, [mol d-2]
  cid_ZCA2P    !soil aqueous CaH2PO4 content micropore non-band, [mol d-2]
  cid_ZMG1P    !soil aqueous MgHPO4 content micropore non-band, [mol d-2]
  cid_H0POB    !soil aqueous PO4 content micropore band, [mol d-2]
  cid_H3POB    !soil aqueous H3PO4 content micropore band, [mol d-2]
  cid_ZFE1PB   !soil aqueous FeHPO4 content micropore band, [mol d-2]
  cid_ZFE2PB   !soil aqueous FeH2PO4 content micropore band, [mol d-2
  cid_ZCA0PB   !soil aqueous CaPO4 content micropore band, [mol d-2]
  cid_ZCA1PB   !soil aqueous CaHPO4 content micropore band, [mol d-2]
  cid_ZCA2PB   !soil aqueous CaH2PO4 content micropore band, [mol d-2]
  cid_ZMG1PB   !soil aqueous MgHPO4 content micropore band, [mol d-2]
  cid_PALPO1   !precipitated AlPO4 non-band, [mol m-3]
  cid_PALPOB   !precipitated AlPO4 band soil, [mol m-3]
  cid_PCAPD1   !precipitated CaHPO4 non-band soil, [mol m-3]
  cid_PCAPDB   !precipitated CaHPO4 band soil, [mol m-3]
  cid_PCAPH1   !precipitated Ca5(PO4)3OH hydroxyapatite non-band soil, [mol m-3]
  cid_PCAPHB   !precipitated Ca5(PO4)3OH hydroxyapatite band soil, [mol m-3]
  cid_PCAPM1   !precipitated Ca(H2PO4)2 non-band soil, [mol m-3]
  cid_PCAPMB   !precipitated CaH2PO4 band soil, [mol m-3]
  cid_PFEPO1   !precipitated FePO4 non-band soil, [mol m-3]
  cid_PFEPOB   !precipitated FePO4 band soil, [mol m-3]
  cid_PALOH    !precipitated Al(OH)3, [mol d-2]
  cid_PFEOH    !precipitated Fe(OH)3, [mol d-2]
  cid_PCACO    !precipitated CaCO3, [mol d-2]
  cid_PCASO    !precipitated CaSO4, [mol d-2]
  cid_XHY      !exchangeable H(+) , [mol d-2]
  cid_XAL      !exchangeable Al, [mol d-2]
  cid_XFE      !exchangeable Fe, [mol d-2]
  cid_XCA      !exchangeable Ca, [mol d-2]
  cid_XMG      !exchangeable Mg, [mol d-2]
  cid_XNA      !exchangeable Na, [mol d-2]
  cid_XKA      !exchangeable K, [mol d-2]
  cid_XHC      !exchangeable COOH , [mol d-2]
  cid_XOH11    !exchangeable OH  non-band, [mol d-2]
  cid_XALO2    !exchangeable AlOH2 , [mol d-2]
  cid_XFEO2    !exchangeable Fe(OH)2, [mol d-2]
  cid_XN41     !exchangeable NH4 non-band soil, [mol d-2]
  cid_XN4B     !exchangeable NH4 band soil, [mol d-2]
  cid_XH01B    !exchangeable OH- band, [mol d-2]
  cid_XOH01    !exchangeable OH- non-band, [mol d-2]
  cid_X1P1B    !exchangeable HPO4 concentration band-soil, [mol m-3]
  cid_X2P1B    !exchangeable H2PO4 concentration band-soil, [mol m-3]
  cid_XH11B    !exchangeable OH band-soil, [mol m-3]
  cid_XH1P1    !exchangeable HPO4  non-band, [mol m-3]
  cid_XH21B    !exchangeable OH2 band-soil, [mol m-3]
  cid_XOH21    !exchangeable OH2  non-band soil, [mol m-3]
  cid_XH2P1    !exchangeable H2PO4  non-band soil, [mol m-3]

  end subroutine Init_Salt_geochem



end module AquachemMod
