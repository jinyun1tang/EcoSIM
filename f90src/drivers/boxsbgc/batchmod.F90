module batchmod

  use abortutils, only : endrun
  use MiscUtilMod, only : addone
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use ModelStatusType, only : model_status_type
  use MicBGCPars, only : micpar
  use MicForcTypeMod, only : micforctype
  use MicFLuxTypeMod, only : micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use fileUtil
implicit none
  private
!  public :: BatchModelConfig
  public :: getvarllen, getvarlist,initmodel
  logical :: Litlayer

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

  integer :: id_ZNH4B     !NH4 band micropore, [g d-2]
  integer :: id_ZNH4S     !NH4 non-band micropore, [g d-2]
  integer :: id_ZNO3B     !NO3 band micropore, [g d-2]
  integer :: id_ZNO3S     !NO3 non-band micropore, [g d-2]
  integer :: id_H1POB     !soil aqueous HPO4 content micropore band, [g d-2]
  integer :: id_H1PO4     !soil aqueous HPO4 content micropore non-band, [g d-2]
  integer :: id_ZNO2B     !NO2(-)  band micropore, [g d-2]
  integer :: id_ZNO2S     !NO2(-)  non-band micropore, [g d-2]
  integer :: id_H2POB     !H2PO4 band micropore, [g d-2]
  integer :: id_H2PO4     !H2PO4 non-band micropore, [g d-2]
  integer :: id_CCO2S     !aqueous CO2 concentration micropore	[g m-3]
  integer :: id_CNO2S     !NO2(-) concentration non-band micropore	[g m-3]
  integer :: id_CNO2B     !NO2 concentration band micropore	[g m-3]
  integer :: id_CZ2OS     !aqueous N2O concentration micropore	[g m-3]
  integer :: id_Z2OS      !aqueous N2O mass in micropore, [g d-2]
  integer :: id_COXYS     !aqueous O2 concentration micropore	[g m-3]
  integer :: id_OXYS      !aqueous O2 mass in micropore	[g m-2]
  integer :: id_COXYG     !gaseous O2 concentration	[g m-3]
  integer :: id_CZ2GS     !gaseous N2 concentration [g m-3]
  integer :: id_CH2GS     !gaseous H2 concentration [g m-3]
  integer :: id_H2GS      !gaseous H2 mass [g m-3]
  integer :: id_CCH4G     !gaseous CH4 concentration in micropore [g m-3]
  integer :: id_CH4S      !aqueous CH4 mass in  micropore	[g d-2]
  integer :: id_ZNFN0     !initial nitrification inhibition activity
  integer :: id_ZNFNI     !current nitrification inhibition activity
  integer :: id_oqc_b,id_oqc_e     !dissolved mass organic C micropore	[gC d-2]
  integer :: id_oqn_b,id_oqn_e     !dissolved mass organic N micropore	[gN d-2]
  integer :: id_oqp_b,id_oqp_e     !dissolved mass organic P micropore	[gP d-2]
  integer :: id_oqa_b,id_oqa_e     !dissolved mass acetate micropore [gC d-2]
  integer :: id_ohc_b,id_ohc_e     !adsorbed mass soil C	[gC d-2]
  integer :: id_ohn_b,id_ohn_e     !adsorbed mass soil N	[gN d-2]
  integer :: id_ohp_b,id_ohp_e     !adsorbed mass soil P	[gP d-2]
  integer :: id_oha_b,id_oha_e     !adsorbed mass soil acetate	[gC d-2]
  integer :: id_osc_b,id_osc_e     !humus mass soil C	[gC d-2]
  integer :: id_osa_b,id_osa_e     !humus mass soil acetate	[gC d-2]
  integer :: id_osn_b,id_osn_e     !humus mass soil N	[gN d-2]
  integer :: id_osp_b,id_osp_e     !humus mass soil P	[gP d-2]
  integer :: id_orc_b,id_orc_e     !microbial mass residue C [gC d-2]
  integer :: id_orn_b,id_orn_e     !microbial mass residue N [gN d-2]
  integer :: id_orp_b,id_orp_e     !microbial mass residue P [gP d-2]
  integer :: id_omc_b,id_omc_e     !microbial biomass C	[gC d-2]
  integer :: id_omn_b,id_omn_e     !microbial biomass N	[gN d-2]
  integer :: id_omp_b,id_omp_e     !microbial biomass P	[gP d-2]
  integer :: id_omcff_b,id_omcff_e   !autotrophic microbial biomass C	[gC d-2]
  integer :: id_omnff_b,id_omnff_e   !autotrophic microbial biomass N	[gN d-2]
  integer :: id_ompff_b,id_ompff_e   !autotrophic microbial biomass C	[gC d-2]
  integer :: fid_ROXYY               !total root + microbial O2 uptake, [g d-2 h-1]
  integer :: fid_RNH4Y               !total root + microbial NH4(+) uptake non-band, [gN d-2 h-1]
  integer :: fid_RNO3Y               !total root + microbial NO3(-) uptake non-band, [gN d-2 h-1]
  integer :: fid_RNO2Y               !total root + microbial NO2(-) uptake non-band, [gN d-2 h-1]
  integer :: fid_RN2OY               !total root + microbial N2O uptake, [g d-2 h-1]
  integer :: fid_RPO4Y               !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  integer :: fid_RP14Y               !HPO4 demand in non-band by all microbial,root,myco populations [gP d-2 h-1]
  integer :: fid_RNHBY               !total root + microbial NH4 uptake band, [gN d-2 h-1]
  integer :: fid_RN3BY               !total root + microbial NO3(-) uptake band, [gN d-2 h-1]
  integer :: fid_RN2BY               !total root + microbial NO2(-) uptake band, [gN d-2 h-1]
  integer :: fid_RPOBY               !total root + microbial PO4 uptake band, [gP d-2 h-1]
  integer :: fid_RP1BY               !HPO4 demand in band by all microbial,root,myco populations [gP d-2 h-1]
  integer :: fid_ROQCY_b,fid_ROQCY_e !total root + microbial DOC uptake, [gC d-2 h-1]
  integer :: fid_ROQAY_b,fid_ROQAY_e !total root + microbial acetate uptake, [gC d-2 h-1]

contains

  function getvarllen()result(nvars)
  implicit none
  integer :: nvars

  call Initboxbgc(nvars)

  end function getvarllen
! ----------------------------------------------------------------------
  subroutine initmodel(nvars, ystates0l, err_status)
!
! set initial conditions for the boxsbgc

  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(inout) :: ystates0l(nvars)
  type(model_status_type), intent(out) :: err_status

  call err_status%reset()


  end subroutine initmodel

! ----------------------------------------------------------------------

!  subroutine BatchModelConfig(nvars,ystatesf0l,forc,micfor,micstt,micflx,err_status)
!
! DESCRIPTION:
! configure the batch mode of the soil bgc
!  use MicStateTraitTypeMod, only : micsttype
!  use MicForcTypeMod, only : micforctype
!  use MicBGCPars, only : micpar
!  implicit none
!  integer, intent(in) :: nvars
!  real(r8), intent(in) :: ystatesf0l(nvars)
!  type(forc_type), intent(in) :: forc
!  type(micforctype), intent(inout) :: micfor
!  type(micsttype), intent(inout) :: micstt
!  type(micfluxtype), intent(inout) :: micflx
!  type(model_status_type), intent(out) :: err_status

!  real(r8), parameter :: ZERO=1.0E-15_r8
!  real(r8), parameter :: ZEROS2=1.0E-08_r8
!  integer :: kk

!  call err_status%reset()

!  associate(                      &
!    nlbiomcp => micpar%nlbiomcp , &
!    ndbiomcp => micpar%ndbiomcp , &
!    NFGs    => micpar%NFGs      , &
!    jcplx1  => micpar%jcplx1    , &
!    JG      => micpar%jguilds     &
!  )
!  micfor%ZERO  =ZERO
!  micfor%ZEROS2=ZERO2
!  micfor%ZEROS =ZERO

!  micfor%CCH4E =CCH4E
!  micfor%COXYE =COXYE
!  micfor%COXQ  =COXQ
!  micfor%COXR  =COXR
!  micfor%FLQRI =FLQRI
!  micfor%FLQRQ =FLQRQ
!  micfor%OFFSET=OFFSET
!  micfor%VOLR  =VOLR
!  micfor%VOLWRX=VOLWRX
!  micfor%VOLY  =VOLY
!  micfor%THETY =THETY
!  micfor%POROS =POROS
!  micfor%FC    =FC
!  micfor%TKS   =TKS
!  micfor%THETW =THETW
!  micfor%PH    =PH
!  micfor%BKVL  =BKVL
!  micfor%VOLX  =VOLX
!  micfor%TFND  =TFND
!  micfor%VLNOB =VLNOB
!  micfor%VLNO3 =VLNO3
!  micfor%VLNH4 =VLNH4
!  micfor%VLNHB =VLNHB
!  micfor%VLPO4 =VLPO4
!  micfor%VLPOB =VLPOB
!  micfor%PSISM =PSISM
!  micfor%OLSGL =OLSGL
!  micfor%ORGC  =ORGC
!  micfor%RNO2Y =ystatesf0l(fid_RNO2Y)
!  micfor%RN2OY =ystatesf0l(fid_RN2OY)
!  micfor%RN2BY =ystatesf0l(fid_RN2BY)
!  micfor%ROXYY =ystatesf0l(fid_ROXYY)
!  micfor%ROXYF =ystatesf0l(fid_ROXYF)
!  micfor%RNHBY =ystatesf0l(fid_RNHBY)
!  micfor%RN3BY =ystatesf0l(fid_RN3BY)
!  micfor%RPOBY =ystatesf0l(fid_RPOBY)
!  micfor%RP1BY =ystatesf0l(fid_RP1BY)
!  micfor%ROQCY(0:jcplx1)=ROQCY(0:jcplx1)
!  micfor%ROQAY(0:jcplx1)=ROQAY(0:jcplx1)
!  micfor%RCH4L =RCH4L
!  micfor%ROXYL = ROXYL
!  micfor%CFOMC =CFOMC(1:ndbiomcp)
!  micfor%litrm=(L==0)
!  micfor%Lsurf=.True.
!  if(micfor%litrm)then
!  !  micstt%ZNH4TU=AMAX1(0.0,ZNH4S)+AMAX1(0.0,ZNH4B)
!  !  micstt%ZNO3TU=AMAX1(0.0,ZNO3S)+AMAX1(0.0,ZNO3B)
!  !  micstt%H1P4TU=AMAX1(0.0,H1PO4)+AMAX1(0.0,H1POB)
!  !  micstt%H2P4TU=AMAX1(0.0,H2PO4)+AMAX1(0.0,H2POB)
!  !  micstt%CNH4BU=CNH4B
!  !  micstt%CNH4SU=CNH4S
!  !  micstt%CH2P4U=CH2P4
!  !  micstt%CH2P4BU=CH2P4B
!  !  micstt%CH1P4U=CH1P4
!  !  micstt%CH1P4BU=CH1P4B
!  !  micstt%CNO3SU=CNO3S
!  !  micstt%CNO3BU=CNO3B
!  !  micstt%OSC13U=OSC(1,3)
!  !  micstt%OSN13U=OSN(1,3)
!  !  micstt%OSP13U=OSP(1,3)
!  !  micstt%OSC14U=OSC(1,4)
!  !  micstt%OSN14U=OSN(1,4)
!  !  micstt%OSP14U=OSP(1,4)
!  !  micstt%OSC24U=OSC(2,4)
!  !  micstt%OSN24U=OSN(2,4)
!  !  micstt%OSP24U=OSP(2,4)
  !  micfor%RNH4YU =RNH4Y
  !  micfor%RNO3YU =RNO3Y
  !  micfor%RPO4YU =RPO4Y
  !  micfor%RP14YU =RP14Y
  !  micfor%VOLWU  =VOLW
  !  micfor%CFOMCU=CFOMC(1:ndbiomcp)
!  else
  !  micfor%AEC=AEC
!  !  micstt%OXYG=OXYG
!  endif
!  micstt%CNH4B =CNH4B
!  micstt%CNH4S =CNH4S
!  micstt%CNO3S =CNO3S
!  micstt%CNO3B =CNO3B
!  micstt%CH2P4 =CH2P4
!  micstt%CH2P4B=CH2P4B
!  micstt%CH1P4 =CH1P4
!  micstt%CH1P4B=CH1P4B
!  micfor%RNH4Y =RNH4Y
!  micfor%RNO3Y =RNO3Y
!  micfor%RPO4Y =RPO4Y
!  micfor%RP14Y =RP14Y
!  micfor%VOLW  =VOLW

!  if(micfor%Lsurf)then
  !  micfor%BKVL0=BKVL
!  endif
!  micfor%DFGS(1:NPH)=DFGS(1:NPH)
!  micfor%FILM(1:NPH)=FILM(1:NPH)
!  micfor%THETPM(1:NPH)=THETPM(1:NPH)
!  micfor%VOLWM(1:NPH)=VOLWM(1:NPH)
!  micfor%TORT(1:NPH)=TORT(1:NPH)
!  micfor%VOLPM(1:NPH)=VOLPM(1:NPH)

!  micstt%EPOC=EPOC
!  micstt%EHUM=EHUM
!  micstt%ZNH4B=ZNH4B
!  micstt%ZNH4S=ZNH4S
!  micstt%ZNO3B=ZNO3B
!  micstt%ZNO3S=ZNO3S
!  micstt%H1POB=H1POB
!  micstt%H1PO4=H1PO4
!  micstt%ZNO2B=ZNO2B
!  micstt%ZNO2S=ZNO2S
!  micstt%H2POB=H2POB
!  micstt%H2PO4=H2PO4
!  micstt%CCO2S=CCO2S
!  micstt%CNO2S=CNO2S
!  micstt%CNO2B=CNO2B
!  micstt%CZ2OS=CZ2OS
!  micstt%Z2OS=Z2OS
!  micstt%COXYS=COXYS
!  micstt%OXYS=OXYS
!  micstt%SOXYL=SOXYL
!  micstt%COXYG=COXYG
!  micstt%CZ2GS=CZ2GS
!  micstt%CH2GS=CH2GS
!  micstt%H2GS=H2GS
!  micstt%CCH4G=CCH4G
!  micstt%CH4S=CH4S
!  micstt%SCH4L=SCH4L
!  micstt%ZNFN0=ZNFN0
!  micstt%ZNFNI=ZNFNI

!  micstt%FOSRH(0:jcplx1)=FOSRH(0:jcplx1)
!  micstt%OQC(0:jcplx1)=OQC(0:jcplx1)
!  micstt%OQN(0:jcplx1)=OQN(0:jcplx1)
!  micstt%OQP(0:jcplx1)=OQP(0:jcplx1)
!  micstt%OQA(0:jcplx1)=OQA(0:jcplx1)
!  micstt%OHC(0:jcplx1)=OHC(0:jcplx1)
!  micstt%OHN(0:jcplx1)=OHN(0:jcplx1)
!  micstt%OHP(0:jcplx1)=OHP(0:jcplx1)
!  micstt%OHA(0:jcplx1)=OHA(0:jcplx1)

!  micstt%OSC(1:jsken,0:jcplx1)=OSC(1:jsken,0:jcplx1)
!  micstt%OSA(1:jsken,0:jcplx1)=OSA(1:jsken,0:jcplx1)
!  micstt%OSN(1:jsken,0:jcplx1)=OSN(1:jsken,0:jcplx1)
!  micstt%OSP(1:jsken,0:jcplx1)=OSP(1:jsken,0:jcplx1)
!  micstt%ORC(1:ndbiomcp,0:jcplx1)=ORC(1:ndbiomcp,0:jcplx1)
!  micstt%ORN(1:ndbiomcp,0:jcplx1)=ORN(1:ndbiomcp,0:jcplx1)
!  micstt%ORP(1:ndbiomcp,0:jcplx1)=ORP(1:ndbiomcp,0:jcplx1)
!  micstt%CNOSC(1:jsken,0:jcplx1)=CNOSC(1:jsken,0:jcplx1)
!  micstt%CPOSC(1:jsken,0:jcplx1)=CPOSC(1:jsken,0:jcplx1)
!  micstt%OMC(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)=OMC(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)
!  micstt%OMN(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)=OMN(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)
!  micstt%OMP(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)=OMP(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1)
!  micstt%OMCff(1:nlbiomcp,1:JG,1:NFGs)=OMCff(1:nlbiomcp,1:JG,1:NFGs)
!  micstt%OMNff(1:nlbiomcp,1:JG,1:NFGs)=OMNff(1:nlbiomcp,1:JG,1:NFGs)
!  micstt%OMPff(1:nlbiomcp,1:JG,1:NFGs)=OMPff(1:nlbiomcp,1:JG,1:NFGs)

!  micflx%RVMXC=RVMXC
!  micflx%RVMBC=RVMBC
!  micflx%RINHOff(1:JG,1:NFGs)=RINHOff(1:JG,1:NFGs)
!  micflx%RINHBff(1:JG,1:NFGs)=RINHBff(1:JG,1:NFGs)
!  micflx%RINOOff(1:JG,1:NFGs)=RINOOff(1:JG,1:NFGs)
!  micflx%RINOBff(1:JG,1:NFGs)=RINOBff(1:JG,1:NFGs)
!  micflx%RIPOOff(1:JG,1:NFGs)=RIPOOff(1:JG,1:NFGs)
!  micflx%RIPBOff(1:JG,1:NFGs)=RIPBOff(1:JG,1:NFGs)
!  micflx%RIPO1ff(1:JG,1:NFGs)=RIPO1ff(1:JG,1:NFGs)
!  micflx%RIPB1ff(1:JG,1:NFGs)=RIPB1ff(1:JG,1:NFGs)
!  micflx%ROXYSff(1:JG,1:NFGs)=ROXYSff(1:JG,1:NFGs)

!  micflx%RINHO(1:JG,1:NFGs,0:JCPLX1)=RINHO(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RINHB(1:JG,1:NFGs,0:JCPLX1)=RINHB(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RINOO(1:JG,1:NFGs,0:JCPLX1)=RINOO(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RINOB(1:JG,1:NFGs,0:JCPLX1)=RINOB(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RIPOO(1:JG,1:NFGs,0:JCPLX1)=RIPOO(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RIPBO(1:JG,1:NFGs,0:JCPLX1)=RIPBO(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RIPO1(1:JG,1:NFGs,0:JCPLX1)=RIPO1(1:JG,1:NFGs,0:JCPLX1)
!  micflx%RIPB1(1:JG,1:NFGs,0:JCPLX1)=RIPB1(1:JG,1:NFGs,0:JCPLX1)
!  micflx%ROXYS(1:JG,1:NFGs,0:JCPLX1)=ROXYS(1:JG,1:NFGs,0:JCPLX1)
!  end associate
!  end subroutine BatchModelConfig

! ----------------------------------------------------------------------
!  subroutine ReadForc(forc)
!
! DESCRIPTION
! read forcing data

!  use MicForcTypeMod, only : micforctype
!  implicit none
!  type(forc_type), intent(in) :: forc


!  end subroutine ReadForc
! ----------------------------------------------------------------------


  subroutine Initboxbgc(nvars)
!
! DESCRIPTION
! Initialize the id of relevant model variables &
! obtain total number of variables
  implicit none

  integer, intent(out) :: nvars
  integer :: itemp

  associate(                        &
    jcplx    => micpar%jcplx      , &
    jcplx1   => micpar%jcplx1     , &
    jsken    => micpar%jsken      , &
    NFGs     => micpar%NFGs       , &
    JG       => micpar%jguilds    , &
    ndbiomcp => micpar%ndbiomcp   , &
    nlbiomcp => micpar%nlbiomcp     &
  )
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

  id_oqc_b=addone(itemp);id_oqc_e=id_oqc_b+jcplx1;itemp=id_oqc_e
  id_oqn_b=addone(itemp);id_oqn_e=id_oqn_b+jcplx1;itemp=id_oqn_e
  id_oqp_b=addone(itemp);id_oqp_e=id_oqp_b+jcplx1;itemp=id_oqp_e
  id_oqa_b=addone(itemp);id_oqa_e=id_oqa_b+jcplx1;itemp=id_oqa_e
  id_ohc_b=addone(itemp);id_ohc_e=id_ohc_b+jcplx1;itemp=id_ohc_e
  id_ohn_b=addone(itemp);id_ohn_e=id_ohn_b+jcplx1;itemp=id_ohn_e
  id_ohp_b=addone(itemp);id_ohp_e=id_ohp_b+jcplx1;itemp=id_ohp_e
  id_oha_b=addone(itemp);id_oha_e=id_oha_b+jcplx1;itemp=id_oha_e
  id_osc_b=addone(itemp);id_osc_e=id_osc_b+jsken*jcplx-1;itemp=id_osc_e
  id_osa_b=addone(itemp);id_osa_e=id_osa_b+jsken*jcplx-1;itemp=id_osa_e
  id_osn_b=addone(itemp);id_osn_e=id_osn_b+jsken*jcplx-1;itemp=id_osn_e
  id_osp_b=addone(itemp);id_osp_e=id_osp_b+jsken*jcplx-1;itemp=id_osp_e
  id_orc_b=addone(itemp);id_orc_e=id_orc_b+ndbiomcp*jcplx-1;itemp=id_orc_e
  id_orn_b=addone(itemp);id_orn_e=id_orn_b+ndbiomcp*jcplx-1;itemp=id_orn_e
  id_orp_b=addone(itemp);id_orp_e=id_orp_b+ndbiomcp*jcplx-1;itemp=id_orp_e
  id_omc_b=addone(itemp);id_omc_e=id_omc_b+nlbiomcp*JG*NFGs*jcplx-1;itemp=id_omc_e
  id_omn_b=addone(itemp);id_omn_e=id_omn_b+nlbiomcp*JG*NFGs*jcplx-1;itemp=id_omn_e
  id_omp_b=addone(itemp);id_omp_e=id_omp_b+nlbiomcp*JG*NFGs*jcplx-1;itemp=id_omp_e
  id_omcff_b=addone(itemp);id_omcff_e=id_omcff_b+nlbiomcp*JG*NFGs-1;itemp=id_omcff_e
  id_omnff_b=addone(itemp);id_omnff_e=id_omnff_b+nlbiomcp*JG*NFGs-1;itemp=id_omnff_e
  id_ompff_b=addone(itemp);id_ompff_e=id_ompff_b+nlbiomcp*JG*NFGs-1;itemp=id_ompff_e

  fid_ROXYY=addone(itemp)
  fid_RNH4Y=addone(itemp)
  fid_RNO3Y=addone(itemp)
  fid_RNO2Y=addone(itemp)
  fid_RN2OY=addone(itemp)
  fid_RPO4Y=addone(itemp)
  fid_RP14Y=addone(itemp)
  fid_RNHBY=addone(itemp)
  fid_RN3BY=addone(itemp)
  fid_RN2BY=addone(itemp)
  fid_RPOBY=addone(itemp)
  fid_RP1BY=addone(itemp)
  fid_ROQCY_b=addone(itemp);fid_ROQCY_e=fid_ROQCY_b+jcplx1;itemp=fid_ROQCY_e
  fid_ROQAY_b=addone(itemp);fid_ROQAY_e=fid_ROQAY_b+jcplx1;itemp=fid_ROQAY_e

  id_ZNH4B=addone(itemp)
  id_ZNH4S=addone(itemp)
  id_ZNO3B=addone(itemp)
  id_ZNO3S=addone(itemp)
  id_H1POB=addone(itemp)
  id_H1PO4=addone(itemp)
  id_ZNO2B=addone(itemp)
  id_ZNO2S=addone(itemp)
  id_H2POB=addone(itemp)
  id_H2PO4=addone(itemp)
  id_CCO2S=addone(itemp)
  id_CNO2S=addone(itemp)
  id_CNO2B=addone(itemp)
  id_CZ2OS=addone(itemp)
  id_Z2OS =addone(itemp)
  id_COXYS=addone(itemp)
  id_OXYS =addone(itemp)
  id_COXYG=addone(itemp)
  id_CZ2GS=addone(itemp)
  id_CH2GS=addone(itemp)
  id_H2GS =addone(itemp)
  id_CCH4G=addone(itemp)
  id_CH4S =addone(itemp)
  id_ZNFN0=addone(itemp)
  id_ZNFNI=addone(itemp)

  nvars=itemp
  end associate
  end subroutine Initboxbgc
! ----------------------------------------------------------------------

  subroutine UpdateStateVars(micfor,micstt,micflx,nvars,ystatesfl)


  implicit none
  type(micforctype), intent(in) :: micfor
  type(micfluxtype), intent(in) :: micflx
  type(micsttype)  , intent(in) :: micstt
  integer , intent(in) :: nvars
  real(r8), intent(inout) :: ystatesfl(nvars)

  real(r8) :: DC,DN,DP,OC,ON,OP
  real(r8) :: ORGC,ORGN,ORGR
  integer :: K,N,NGL,M
  associate(                        &
    jcplx1    => micpar%jcplx1    , &
    jcplx     => micpar%jcplx     , &
    JG        => micpar%jguilds   , &
    NFGs      => micpar%NFGs      , &
    jsken     => micpar%jsken     , &
    nlbiomcp => micpar%nlbiomcp   , &
    ndbiomcp => micpar%ndbiomcp   , &
    is_litter => micpar%is_litter , &
    VOLW  => micfor%VOLW            &
  )
!atmospheric gaseous CO2,CH4,O2,NH3,N2,N2O,H2

  ystatesfl(id_ZNH4B)=micstt%ZNH4B
  ystatesfl(id_ZNH4S)=micstt%ZNH4S
  ystatesfl(id_ZNO3B)=micstt%ZNO3B
  ystatesfl(id_ZNO3S)=micstt%ZNO3S
  ystatesfl(id_H1POB)=micstt%H1POB
  ystatesfl(id_H1PO4)=micstt%H1PO4
  ystatesfl(id_ZNO2B)=micstt%ZNO2B
  ystatesfl(id_ZNO2S)=micstt%ZNO2S
  ystatesfl(id_H2POB)=micstt%H2POB
  ystatesfl(id_H2PO4)=micstt%H2PO4
  ystatesfl(id_CCO2S)=micstt%CCO2S !-RCO2O/VOLW
  ystatesfl(id_CNO2S)=micstt%ZNO2S/(VOLW*micfor%VLNO3)
  ystatesfl(id_CNO2B)=micstt%ZNO2B/(VOLW*micfor%VLNOB)
  ystatesfl(id_CZ2OS)=micstt%Z2OS/VOLW
  ystatesfl(id_Z2OS) =micstt%Z2OS
  ystatesfl(id_COXYS)=micstt%OXYS/VOLW
  ystatesfl(id_OXYS) =micstt%OXYS  !-RUPOXO
  ystatesfl(id_COXYG)=micstt%COXYG
  ystatesfl(id_CZ2GS)=micstt%CZ2GS
  ystatesfl(id_CH2GS)=micstt%CH2GS
  ystatesfl(id_H2GS) =micstt%H2GS
  ystatesfl(id_CCH4G)=micstt%CCH4G
  ystatesfl(id_CH4S) =micstt%CH4S  !-RCH4O
  ystatesfl(id_ZNFN0)=micstt%ZNFN0
  ystatesfl(id_ZNFNI)=micstt%ZNFNI

! the following variables are updated in the microbial model
  ystatesfl(id_oqc_b:id_oqc_e)=micstt%OQC(0:jcplx1)
  ystatesfl(id_oqn_b:id_oqn_e)=micstt%OQN(0:jcplx1)
  ystatesfl(id_oqp_b:id_oqp_e)=micstt%OQP(0:jcplx1)
  ystatesfl(id_oqa_b:id_oqa_e)=micstt%OQA(0:jcplx1)
  ystatesfl(id_ohc_b:id_ohc_e)=micstt%OHC(0:jcplx1)
  ystatesfl(id_ohn_b:id_ohn_e)=micstt%OHN(0:jcplx1)
  ystatesfl(id_ohp_b:id_ohp_e)=micstt%OHP(0:jcplx1)
  ystatesfl(id_oha_b:id_oha_e)=micstt%OHA(0:jcplx1)
  ystatesfl(id_osc_b:id_osc_e)=reshape(micstt%OSC(1:jsken,0:jcplx1), &
    (/jsken*jcplx/))
  ystatesfl(id_osa_b:id_osa_e)=reshape(micstt%OSA(1:jsken,0:jcplx1), &
    (/jsken*jcplx/))
  ystatesfl(id_osn_b:id_osn_e)=reshape(micstt%OSN(1:jsken,0:jcplx1), &
    (/jsken*jcplx/))
  ystatesfl(id_osp_b:id_osp_e)=reshape(micstt%OSP(1:jsken,0:jcplx1), &
    (/jsken*jcplx/))
  ystatesfl(id_orc_b:id_orc_e)=reshape(micstt%ORC(1:ndbiomcp,0:jcplx1),&
    (/ndbiomcp*jcplx/))
  ystatesfl(id_orn_b:id_orn_e)=reshape(micstt%ORN(1:ndbiomcp,0:jcplx1),&
    (/ndbiomcp*jcplx/))
  ystatesfl(id_orp_b:id_orp_e)=reshape(micstt%ORP(1:ndbiomcp,0:jcplx1),&
    (/ndbiomcp*jcplx/))
  ystatesfl(id_omc_b:id_omc_e)=reshape(micstt%OMC(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1),&
    (/nlbiomcp*JG*NFGs*jcplx/))
  ystatesfl(id_omn_b:id_omn_e)=reshape(micstt%OMN(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1),&
    (/nlbiomcp*JG*NFGs*jcplx/))
  ystatesfl(id_omp_b:id_omp_e)=reshape(micstt%OMP(1:nlbiomcp,1:JG,1:NFGs,0:jcplx1),&
    (/nlbiomcp*JG*NFGs*jcplx/))
  ystatesfl(id_omcff_b:id_omcff_e)=reshape(micstt%OMCff(1:nlbiomcp,1:JG,1:NFGs),&
    (/nlbiomcp*JG*NFGs/))
  ystatesfl(id_omnff_b:id_omnff_e)=reshape(micstt%OMNff(1:nlbiomcp,1:JG,1:NFGs),&
    (/nlbiomcp*JG*NFGs/))
  ystatesfl(id_ompff_b:id_ompff_e)=reshape(micstt%OMPff(1:nlbiomcp,1:JG,1:NFGs),&
    (/nlbiomcp*JG*NFGs/))

! summarize diagnostic fluxes
  DO K=0,jcplx1
    IF(.not.micfor%litrm.or.(micpar%is_litter(K)))THEN
      DO N=1,NFGs
        DO NGL=1,JG
          ystatesfl(fid_ROXYY)=ystatesfl(fid_ROXYY)!+ROXYS(NGL,N,K)
          ystatesfl(fid_RNH4Y)=ystatesfl(fid_RNH4Y)!+RVMX4(NGL,N,K)+RINHO(NGL,N,K)
          ystatesfl(fid_RNO3Y)=ystatesfl(fid_RNO3Y)!+RVMX3(NGL,N,K)+RINOO(NGL,N,K)
          ystatesfl(fid_RNO2Y)=ystatesfl(fid_RNO2Y)!+RVMX2(NGL,N,K)
          ystatesfl(fid_RN2OY)=ystatesfl(fid_RN2OY)!+RVMX1(NGL,N,K)
          ystatesfl(fid_RPO4Y)=ystatesfl(fid_RPO4Y)!+RIPOO(NGL,N,K)
          ystatesfl(fid_RP14Y)=ystatesfl(fid_RP14Y)!+RIPO1(NGL,N,K)
          ystatesfl(fid_RNHBY)=ystatesfl(fid_RNHBY)!+RVMB4(NGL,N,K)+RINHB(NGL,N,K)
          ystatesfl(fid_RN3BY)=ystatesfl(fid_RN3BY)!+RVMB3(NGL,N,K)+RINOB(NGL,N,K)
          ystatesfl(fid_RN2BY)=ystatesfl(fid_RN2BY)!+RVMB2(NGL,N,K)
          ystatesfl(fid_RPOBY)=ystatesfl(fid_RPOBY)!+RIPBO(NGL,N,K)
          ystatesfl(fid_RP1BY)=ystatesfl(fid_RP1BY)!+RIPB1(NGL,N,K)
          ystatesfl(fid_ROQCY_b+K)=ystatesfl(fid_ROQCY_b+K)!+ROQCS(NGL,N,K)
          ystatesfl(fid_ROQAY_b+K)=ystatesfl(fid_ROQAY_b+K)!+ROQAS(NGL,N,K)
        enddo
      ENDDO
    ENDIF
  ENDDO

  DO  N=1,NFGs
    DO NGL=1,JG
      ystatesfl(fid_ROXYY)=ystatesfl(fid_ROXYY) !+ROXYSff(NGL,N)
      ystatesfl(fid_RNH4Y)=ystatesfl(fid_RNH4Y) !+RVMX4ff(NGL,N)+RINHOff(NGL,N)
      ystatesfl(fid_RNO3Y)=ystatesfl(fid_RNO3Y) !+RVMX3ff(NGL,N)+RINOOff(NGL,N)
      ystatesfl(fid_RNO2Y)=ystatesfl(fid_RNO2Y) !+RVMX2ff(NGL,N)
      ystatesfl(fid_RN2OY)=ystatesfl(fid_RN2OY) !+RVMX1ff(NGL,N)
      ystatesfl(fid_RPO4Y)=ystatesfl(fid_RPO4Y) !+RIPOOff(NGL,N)
      ystatesfl(fid_RP14Y)=ystatesfl(fid_RP14Y) !+RIPO1ff(NGL,N)
      ystatesfl(fid_RNHBY)=ystatesfl(fid_RNHBY) !+RVMB4ff(NGL,N)+RINHBff(NGL,N)
      ystatesfl(fid_RN3BY)=ystatesfl(fid_RN3BY) !+RVMB3ff(NGL,N)+RINOBff(NGL,N)
      ystatesfl(fid_RN2BY)=ystatesfl(fid_RN2BY) !+RVMB2ff(NGL,N)
      ystatesfl(fid_RPOBY)=ystatesfl(fid_RPOBY) !+RIPBOff(NGL,N)
      ystatesfl(fid_RP1BY)=ystatesfl(fid_RP1BY) !+RIPB1ff(NGL,N)
    enddo
  ENDDO

  ystatesfl(fid_RNO2Y)=ystatesfl(fid_RNO2Y) !+RVMXC
  ystatesfl(fid_RN2BY)=ystatesfl(fid_RN2BY) !+RVMBC

  DC=0.0_r8
  DN=0.0_r8
  DP=0.0_r8
  OC=0.0_r8
  ON=0.0_r8
  OP=0.0_r8

  DO K=0,jcplx1
    IF(is_litter(K))THEN
      DO N=1,NFGs
        DO  M=1,nlbiomcp
          DO NGL=1,JG
!            DC=DC+OMC(M,NGL,N,K)
!            DN=DN+OMN(M,NGL,N,K)
!            DP=DP+OMP(M,NGL,N,K)
          ENDDO
        enddo
      ENDDO
    ELSE
      DO N=1,NFGs
        DO  M=1,nlbiomcp
          DO NGL=1,JG
!            OC=OC+OMC(M,NGL,N,K)
!            ON=ON+OMN(M,NGL,N,K)
!            OP=OP+OMP(M,NGL,N,K)
          enddo
        enddo
      ENDDO
    ENDIF
  ENDDO
! abstract complex
  DO  N=1,NFGs
    DO  M=1,nlbiomcp
      DO NGL=1,JG
!        OC=OC+OMCff(M,NGL,N)
!        ON=ON+OMNff(M,NGL,N)
!        OP=OP+OMPff(M,NGL,N)
      enddo
    enddo
  ENDDO
! microbial residue
  DO K=0,jcplx1
    IF(is_litter(K))THEN
      DO M=1,ndbiomcp
!        DC=DC+ORC(M,K)
!        DN=DN+ORN(M,K)
!        DP=DP+ORP(M,K)
      ENDDO
!      DC=DC+OQC(K)+OQCH(K)+OHC(K)+OQA(K)+OQAH(K)+OHA(K)
!      DN=DN+OQN(K)+OQNH(K)+OHN(K)
!      DP=DP+OQP(K)+OQPH(K)+OHP(K)
      DO M=1,jsken
!        DC=DC+OSC(M,K)
!        DN=DN+OSN(M,K)
!        DP=DP+OSP(M,K)
      ENDDO
    ELSE
      DO M=1,ndbiomcp
!        OC=OC+ORC(M,K)
!        ON=ON+ORN(M,K)
!        OP=OP+ORP(M,K)
      ENDDO
!      OC=OC+OQC(K)+OQCH(K)+OHC(K)+OQA(K)+OQAH(K)+OHA(K)
!      ON=ON+OQN(K)+OQNH(K)+OHN(K)
!      OP=OP+OQP(K)+OQPH(K)+OHP(K)
      DO M=1,jsken
!        OC=OC+OSC(M,K)
!        ON=ON+OSN(M,K)
!        OP=OP+OSP(M,K)
      ENDDO
    ENDIF
  ENDDO
! DC is for litter complex, and OC is for POM and humus complex
  ORGC=DC+OC
  ORGN=DN+ON
  ORGR=DC

  end associate
  end subroutine UpdateStateVars

! ----------------------------------------------------------------------

  subroutine getvarlist(nvars, varl, varlnml, unitl, vartypes)

  use histMod, only : hist_var_str_len,hist_unit_str_len,hist_var_lon_str_len
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(out) :: varl(nvars)            !variable name
  character(len=hist_var_lon_str_len), intent(out) :: varlnml(nvars)     !variable name
  character(len=hist_unit_str_len),intent(out) :: unitl(nvars)           !variable unit
  integer                         ,intent(out) :: vartypes(nvars)        !variable type, flx or state

  integer :: iknen,icplx
  integer :: jj,ll,k,m,n,ngl

  associate(                        &
    jcplx1    => micpar%jcplx1    , &
    JG        => micpar%jguilds   , &
    jsken     => micpar%jsken     , &
    NFGs      => micpar%NFGs      , &
    nlbiomcp  => micpar%nlbiomcp  , &
    ndbiomcp  => micpar%ndbiomcp    &
  )


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

  varl(id_ZNH4B)='ZNH4B';varlnml(id_ZNH4B)='band soil micropore NH4(+) mass'
  unitl(id_ZNH4B)='gN d-2';vartypes(id_ZNH4B)=var_state_type

  varl(id_ZNH4S)='ZNH4S';varlnml(id_ZNH4S)='non-band soil micropore NH4(+) mass'
  unitl(id_ZNH4S)='gN d-2';vartypes(id_ZNH4S)=var_state_type

  varl(id_ZNO3B)='ZNO3B';varlnml(id_ZNO3B)='band soil micropore NO3(-) mass'
  unitl(id_ZNO3B)='gN d-2';vartypes(id_ZNO3B)=var_state_type

  varl(id_ZNO3S)='ZNO3S';varlnml(id_ZNO3S)='non-band soil micropore NO3(-) mass'
  unitl(id_ZNO3S)='gN d-2';vartypes(id_ZNO3S)=var_state_type

  varl(id_H1POB)='H1POB';varlnml(id_H1POB)='band soil micropore aqueous HPO4 content'
  unitl(id_H1POB)='gP m-2';vartypes(id_H1POB)=var_state_type

  varl(id_H1PO4)='H1PO4';varlnml(id_H1PO4)='non-band soil micropore aqueous HPO4(--) content';
  unitl(id_H1PO4)='gP m-2';vartypes(id_H1PO4)=var_state_type

  varl(id_ZNO2B)='ZNO2B';varlnml(id_ZNO2B)='band soil micropore NO2(-) mass'
  unitl(id_ZNO2B)='gN m-2';vartypes(id_ZNO2B)=var_state_type

  varl(id_ZNO2S)='ZNO2S';varlnml(id_ZNO2S)='non-band soil micropore NO2(-) mass'
  unitl(id_ZNO2S)='gN m-2';vartypes(id_ZNO2S)=var_state_type

  varl(id_H2POB)='H2POB';varlnml(id_H2POB)='band soil micropore H2PO4 mass'
  unitl(id_H2POB)='gP m-2';vartypes(id_H2POB)=var_state_type

  varl(id_H2PO4)='H2PO4';varlnml(id_H2PO4)='non-band soil micropore H2PO4 mass'
  unitl(id_H2PO4)='gP m-2';vartypes(id_H2PO4)=var_state_type

  varl(id_CCO2S)='CCO2S';varlnml(id_CCO2S)='soil micropore aqueous CO2 concentration'
  unitl(id_CCO2S)='gC m-3';vartypes(id_CCO2S)=var_state_type

  varl(id_CNO2S)='CNO2S';varlnml(id_CNO2S)='non-band soil micropore NO2 concentration'
  unitl(id_CNO2S)='gN m-3';vartypes(id_CNO2S)=var_state_type

  varl(id_CNO2B)='CNO2B';varlnml(id_CNO2B)='band soil micropore NO2 concentration'
  unitl(id_CNO2B)='gN m-3';vartypes(id_CNO2B)=var_state_type

  varl(id_CZ2OS)='CZ2OS';varlnml(id_CZ2OS)='soil micropore aqueous N2O concentration'
  unitl(id_CZ2OS)='gN m-3';vartypes(id_CZ2OS)=var_state_type

  varl(id_Z2OS) ='Z2OS';varlnml(id_Z2OS)='soil micropore aqueous N2O mass'
  unitl(id_Z2OS)='gN d-2';vartypes(id_Z2OS)=var_state_type

  varl(id_COXYS)='COXYS';varlnml(id_COXYS)='soil micropore aqueous O2 concentration'
  unitl(id_COXYS)='g m-3';vartypes(id_COXYS)=var_state_type

  varl(id_OXYS) ='OXYS';varlnml(id_OXYS)='soil micropore aqueous O2 mass'
  unitl(id_OXYS)='g d-2';vartypes(id_OXYS)=var_state_type

  varl(id_COXYG)='COXYG';varlnml(id_COXYG)='soil micropore gaseous O2 concentration'
  unitl(id_COXYG)='g m-3';vartypes(id_COXYG)=var_state_type

  varl(id_CZ2GS)='CZ2GS';varlnml(id_CZ2GS)='soil micropore aqueous N2 concentration'
  unitl(id_CZ2GS)='gN m-3';vartypes(id_CZ2GS)=var_state_type

  varl(id_CH2GS)='CH2GS';varlnml(id_CH2GS)='soil micropore aqueous H2 concentration'
  unitl(id_CH2GS)='g m-3';vartypes(id_CH2GS)=var_state_type

  varl(id_H2GS) ='H2GS';varlnml(id_H2GS)='soil micropore aqueous H2 mass'
  unitl(id_H2GS)='g d-2';vartypes(id_H2GS)=var_state_type

  varl(id_CCH4G)='CCH4G';varlnml(id_CCH4G)='soil micropore gaseous CH4 concentration'
  unitl(id_CCH4G)='gC m-3';vartypes(id_CCH4G)=var_state_type

  varl(id_CH4S) ='CH4S';varlnml(id_CH4S)='soil micropore aqueous CH4 mass'
  unitl(id_CH4S)='gC d-2';vartypes(id_CH4S)=var_state_type

  varl(id_ZNFN0)='ZNFN0';varlnml(id_ZNFN0)='initial nitrification inhibition activity'
  unitl(id_ZNFN0)='none';vartypes(id_ZNFN0)=var_state_type

  varl(id_ZNFNI)='ZNFNI';varlnml(id_ZNFNI)='current nitrification inhibition activity'
  unitl(id_ZNFNI)='none';vartypes(id_ZNFNI)=var_state_type

  do jj=id_oqc_b,id_oqc_e
    write(varl(jj),'(A,I1)')'OQC',jj-id_oqc_b
    varlnml(jj)='micropore dissolved organic C mass in complex '//trim(micpar%cplxname(jj-id_oqc_b))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=id_oqn_b,id_oqn_e
    write(varl(jj),'(A,I1)')'OQN',jj-id_oqn_b
    varlnml(jj)='micropore dissolved N mass in complex '//trim(micpar%cplxname(jj-id_oqn_b))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=id_oqp_b,id_oqp_e
    write(varl(jj),'(A,I1)')'OQP',jj-id_oqn_b
    varlnml(jj)='micropore dissolved N mass in complex '//trim(micpar%cplxname(jj-id_oqp_b))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=id_oqa_b,id_oqa_e
    write(varl(jj),'(A,I1)')'OQA',jj-id_oqa_b
    varlnml(jj)='micropore dissolved acetate mass in complex '//trim(micpar%cplxname(jj-id_oqa_b))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=id_ohc_b,id_ohc_e
    write(varl(jj),'(A,I1)')'OHC',jj-id_ohc_b
    varlnml(jj)='adsorbed soil C mass in complex'//trim(micpar%cplxname(jj-id_ohc_b))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=id_ohn_b,id_ohn_e
    write(varl(jj),'(A,I1)')'OHN',jj-id_ohn_b
    varlnml(jj)='adsorbed soil N mass in complex'//trim(micpar%cplxname(jj-id_ohn_b))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=id_ohp_b,id_ohp_e
    write(varl(jj),'(A,I1)')'OHP',jj-id_ohp_b
    varlnml(jj)='adsorbed soil P mass in complex'//trim(micpar%cplxname(jj-id_ohp_b))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=id_oha_b,id_oha_e
    write(varl(jj),'(A,I1)')'OHA',jj-id_ohp_b
    varlnml(jj)='adsorbed soil acetate mass in complex'//trim(micpar%cplxname(jj-id_oha_b))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=id_osc_b,id_osc_e
    iknen=jj-id_osc_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=mod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSC',iknen+1,icplx
    varlnml(jj)='humus soil C as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=id_osn_b,id_osn_e
    iknen=jj-id_osn_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=mod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSN',iknen+1,icplx
    varlnml(jj)='humus soil N as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=id_osp_b,id_osp_e
    iknen=jj-id_osp_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=mod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSP',iknen+1,icplx
    varlnml(jj)='humus soil P as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=id_osa_b,id_osa_e
    iknen=jj-id_osc_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=mod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSA',iknen+1,icplx
    varlnml(jj)='colonized humus soil C as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=id_orc_b,id_orc_e
    iknen=jj-id_orc_b
    icplx=floor((iknen-1.e-3_r8)/ndbiomcp)
    iknen=mod(iknen,ndbiomcp)
    write(varl(jj),'(A,I1,I1)')'ORC',iknen+1,icplx
    varlnml(jj)='microbial residue C as '//trim(micpar%micresb(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=id_orn_b,id_orn_e
    iknen=jj-id_orc_b
    icplx=floor((iknen-1.e-3_r8)/ndbiomcp)
    iknen=mod(iknen,ndbiomcp)
    write(varl(jj),'(A,I1,I1)')'ORN',iknen+1,icplx
    varlnml(jj)='microbial residue N as '//trim(micpar%micresb(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=id_orp_b,id_orp_e
    iknen=jj-id_orc_b
    icplx=floor((iknen-1.e-3_r8)/ndbiomcp)
    iknen=mod(iknen,ndbiomcp)
    write(varl(jj),'(A,I1,I1)')'ORP',iknen+1,icplx
    varlnml(jj)='microbial residue P as '//trim(micpar%micresb(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo


  jj=0
  DO k=0,jcplx1
  DO N=1,NFGs
  DO NGL=1,JG
  DO M=1,nlbiomcp
    ll=id_omc_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMC'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass C in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))//' in complex '//trim(micpar%cplxname(k))
    unitl(ll)='gC d-2'
    vartypes(ll)=var_state_type

    ll=id_omn_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMN'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass N in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))//' in complex '//trim(micpar%cplxname(k))
    unitl(ll)='gN d-2'
    vartypes(ll)=var_state_type

    ll=id_omp_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMP'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass P in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))//' in complex '//trim(micpar%cplxname(k))
    unitl(ll)='gP d-2'
    vartypes(ll)=var_state_type
    jj=jj+1
  enddo
  ENDDO
  ENDDO
  ENDDO

  jj=0
  DO N=1,NFGs
  DO NGL=1,JG
  DO M=1,nlbiomcp
    ll=id_omcff_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMC'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass C in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))
    unitl(ll)='gC d-2'
    vartypes(ll)=var_state_type

    ll=id_omnff_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMN'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass N in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))
    unitl(ll)='gN d-2'
    vartypes(ll)=var_state_type

    ll=id_ompff_b+jj
    write(varl(ll),'(A,I2.2,A)')'OMP'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass P in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))
    unitl(ll)='gP d-2'
    vartypes(ll)=var_state_type

    jj=jj+1
  enddo
  ENDDO
  ENDDO

  varl(fid_ROXYY)='ROXYY';varlnml(fid_ROXYY)='total root + microbial O2 uptake potential'
  unitl(fid_ROXYY)='g d-2 h-1'; vartypes(fid_ROXYY)=var_flux_type

  varl(fid_RNH4Y)='RNH4Y';varlnml(fid_RNH4Y)='total root + microbial NH4 uptake potential non-band soil'
  unitl(fid_RNH4Y)='gN d-2 h-1'; vartypes(fid_RNH4Y)=var_flux_type

  varl(fid_RNO3Y)='RNO3Y';varlnml(fid_RNO3Y)='total root + microbial NO3 uptake potential non-band soil'
  unitl(fid_RNO3Y)='gN d-2 h-1'; vartypes(fid_RNO3Y)=var_flux_type

  varl(fid_RNO2Y)='RNO2Y';varlnml(fid_RNO2Y)='total root + microbial NO2 uptake potential non-band soil'
  unitl(fid_RNO2Y)='gN d-2 h-1'; vartypes(fid_RNO2Y)=var_flux_type

  varl(fid_RN2OY)='RN2OY';varlnml(fid_RN2OY)='total root + microbial N2O uptake potential';
  unitl(fid_RN2OY)='gN d-2 h-1'; vartypes(fid_RN2OY)=var_flux_type

  varl(fid_RPO4Y)='RPO4Y';varlnml(fid_RPO4Y)='total root + microbial PO4 uptake potential non-band soil'
  unitl(fid_RPO4Y)='gP d-2 h-1'; vartypes(fid_RPO4Y)=var_flux_type

  varl(fid_RP14Y)='RP14Y';varlnml(fid_RP14Y)='total root + microbial HPO4 uptake non-band soil'
  unitl(fid_RP14Y)='gP d-2 h-1'; vartypes(fid_RP14Y)=var_flux_type

  varl(fid_RNHBY)='RNHBY';varlnml(fid_RNHBY)='total root + microbial NH4 uptake potential band soil'
  unitl(fid_RNHBY)='gN d-2 h-1'; vartypes(fid_RNHBY)=var_flux_type

  varl(fid_RN3BY)='RN3BY';varlnml(fid_RN3BY)='total root + microbial NO3 uptake potential band soil'
  unitl(fid_RN3BY)='gN d-2 h-1'; vartypes(fid_RN3BY)=var_flux_type

  varl(fid_RN2BY)='RN2BY';varlnml(fid_RN2BY)='total root + microbial NO2 uptake potential band soil'
  unitl(fid_RN2BY)='gN d-2 h-1'; vartypes(fid_RN2BY)=var_flux_type

  varl(fid_RPOBY)='RPOBY';varlnml(fid_RPOBY)='total root + microbial PO4 uptake potential band soil'
  unitl(fid_RPOBY)='gP d-2 h-1'; vartypes(fid_RPOBY)=var_flux_type

  varl(fid_RP1BY)='RP1BY';varlnml(fid_RP1BY)='total root + microbial HPO4 uptake potential band soil';
  unitl(fid_RP1BY)='gP d-2 h-1'; vartypes(fid_RP1BY)=var_flux_type

  do jj =fid_ROQCY_b,fid_ROQCY_e
    write(varl(jj),'(A,I2.2)')'ROQCY',jj-fid_ROQCY_b
    varlnml(jj)='total root + microbial DOC uptake in complex ' &
      //micpar%cplxname(jj-fid_ROQCY_b)
    vartypes(jj)=var_flux_type
    unitl(jj)='gC d-2 h-1'
  enddo
  do jj =fid_ROQAY_b,fid_ROQAY_e
    write(varl(jj),'(A,I2.2)')'ROQAY',jj-fid_ROQAY_b
    varlnml(jj)='total root + microbial acetate uptake in complex ' &
      //micpar%cplxname(jj-fid_ROQCY_b)
    vartypes(jj)=var_flux_type
    unitl(jj)='gC d-2 h-1'
  enddo
  end associate
  end subroutine getvarlist
! ----------------------------------------------------------------------

end module batchmod
