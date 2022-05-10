module AqueChemDatatype
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable ::  CAL(:,:,:)                         !soil Al content, [mg kg-1]
  real(r8),allocatable ::  CFE(:,:,:)                         !soil Fe content, [mg kg-1]
  real(r8),allocatable ::  CCA(:,:,:)                         !soil Ca content, [mg kg-1]
  real(r8),allocatable ::  CMG(:,:,:)                         !soil Mg content, [mg kg-1]
  real(r8),allocatable ::  CNA(:,:,:)                         !soil Na content, [mg kg-1]
  real(r8),allocatable ::  CKA(:,:,:)                         !soil K content, [mg kg-1]
  real(r8),allocatable ::  CSO4(:,:,:)                        !soil SO4 content, [mg kg-1]
  real(r8),allocatable ::  CCL(:,:,:)                         !soil Cl content, [mg kg-1]
  real(r8),allocatable ::  CALOH(:,:,:)                       !soil AlOH3 content, [mg kg-1]
  real(r8),allocatable ::  CFEOH(:,:,:)                       !soil FeOH3 content, [mg kg-1]
  real(r8),allocatable ::  CCACO(:,:,:)                       !soil CaCO3 content, [mg kg-1]
  real(r8),allocatable ::  CCASO(:,:,:)                       !soil CaSO4 content, [mg kg-1]
  real(r8),allocatable ::  CALPO(:,:,:)                       !soil AlPO4 content, [mg kg-1]
  real(r8),allocatable ::  CFEPO(:,:,:)                       !soil FePO4 content, [mg kg-1]
  real(r8),allocatable ::  CCAPD(:,:,:)                       !soil CaHPO4 content, [mg kg-1]
  real(r8),allocatable ::  CCAPH(:,:,:)                       !soil apatite content, [mg kg-1]
  real(r8),allocatable ::  GKC4(:,:,:)                        !Ca-NH4 Gapon selectivity coefficient, [-]
  real(r8),allocatable ::  GKCH(:,:,:)                        !Ca-H Gapon selectivity coefficient, [-]
  real(r8),allocatable ::  GKCA(:,:,:)                        !Ca-Al Gapon selectivity coefficient, [-]
  real(r8),allocatable ::  GKCM(:,:,:)                        !Ca-Mg Gapon selectivity coefficient, [-]
  real(r8),allocatable ::  GKCN(:,:,:)                        !Ca-Na Gapon selectivity coefficient, [-]
  real(r8),allocatable ::  GKCK(:,:,:)                        !Ca-K Gapon selectivity coefficient, [-]
  real(r8),allocatable ::  ZCO3(:,:,:)                        !soil aqueous CO3 content micropore, [mol d-2]
  real(r8),allocatable ::  ZHCO3(:,:,:)                       !soil aqueous HCO3 content micropore, [mol d-2]
  real(r8),allocatable ::  XN4(:,:,:)                         !exchangeable NH4 non-band, [mol d-2]
  real(r8),allocatable ::  XNB(:,:,:)                         !exchangeable NH4 band, [mol d-2]
  real(r8),allocatable ::  ZAL(:,:,:)                         !soil aqueous Al content micropore, [mol d-2]
  real(r8),allocatable ::  ZFE(:,:,:)                         !soil aqueous Fe content micropore, [mol d-2]
  real(r8),allocatable ::  ZHY(:,:,:)                         !soil aqueous H content micropore, [mol d-2]
  real(r8),allocatable ::  ZCA(:,:,:)                         !soil aqueous Ca content micropore, [mol d-2]
  real(r8),allocatable ::  ZMG(:,:,:)                         !soil aqueous Mg content micropore, [mol d-2]
  real(r8),allocatable ::  ZNA(:,:,:)                         !soil aqueous Na content micropore, [mol d-2]
  real(r8),allocatable ::  ZKA(:,:,:)                         !soil aqueous K content micropore, [mol d-2]
  real(r8),allocatable ::  ZOH(:,:,:)                         !soil aqueous OH content micropore, [mol d-2]
  real(r8),allocatable ::  ZSO4(:,:,:)                        !soil aqueous SO4 content micropore, [mol d-2]
  real(r8),allocatable ::  ZCL(:,:,:)                         !soil aqueous Cl content micropore, [mol d-2]
  real(r8),allocatable ::  ZALOH1(:,:,:)                      !soil aqueous AlOH content micropore, [mol d-2]
  real(r8),allocatable ::  ZALOH2(:,:,:)                      !soil aqueous AlOH2 content micropore, [mol d-2]
  real(r8),allocatable ::  ZALOH3(:,:,:)                      !soil aqueous AlOH3 content micropore, [mol d-2]
  real(r8),allocatable ::  ZALOH4(:,:,:)                      !soil aqueous AlOH4 content micropore, [mol d-2]
  real(r8),allocatable ::  ZALS(:,:,:)                        !soil aqueous AlSO4 content micropore, [mol d-2]
  real(r8),allocatable ::  ZFEOH1(:,:,:)                      !soil aqueous FeOH content micropore, [mol d-2]
  real(r8),allocatable ::  ZFEOH2(:,:,:)                      !soil aqueous FeOH2 content micropore, [mol d-2]
  real(r8),allocatable ::  ZFEOH3(:,:,:)                      !soil aqueous FeOH3 content micropore, [mol d-2]
  real(r8),allocatable ::  ZFEOH4(:,:,:)                      !soil aqueous FeOH4 content micropore, [mol d-2]
  real(r8),allocatable ::  ZFES(:,:,:)                        !soil aqueous FeSO4 content micropore, [mol d-2]
  real(r8),allocatable ::  ZCAO(:,:,:)                        !soil aqueous CaOH2 content micropore, [mol d-2]
  real(r8),allocatable ::  ZCAC(:,:,:)                        !soil aqueous CACO3 content micropore, [mol d-2]
  real(r8),allocatable ::  ZCAH(:,:,:)                        !soil aqueous CaHCO3 content micropore, [mol d-2]
  real(r8),allocatable ::  ZCAS(:,:,:)                        !soil aqueous CaSO4 content micropore, [mol d-2]
  real(r8),allocatable ::  ZMGO(:,:,:)                        !soil aqueous MgOH content micropore, [mol d-2]
  real(r8),allocatable ::  ZMGC(:,:,:)                        !soil aqueous MgCO3 content micropore, [mol d-2]
  real(r8),allocatable ::  ZMGH(:,:,:)                        !soil aqueous MgHCO3 content micropore, [mol d-2]
  real(r8),allocatable ::  ZMGS(:,:,:)                        !soil aqueous MgSO4 content micropore, [mol d-2]
  real(r8),allocatable ::  ZNAC(:,:,:)                        !soil aqueous NaCO3 content micropore, [mol d-2]
  real(r8),allocatable ::  ZNAS(:,:,:)                        !soil aqueous NaSO4 content micropore, [mol d-2]
  real(r8),allocatable ::  ZKAS(:,:,:)                        !soil aqueous KSO4 content micropore, [mol d-2]
  real(r8),allocatable ::  H0PO4(:,:,:)                       !soil aqueous PO4 content micropore non-band, [mol d-2]
  real(r8),allocatable ::  H3PO4(:,:,:)                       !soil aqueous H3PO4 content micropore non-band, [mol d-2]
  real(r8),allocatable ::  ZFE1P(:,:,:)                       !soil aqueous FeHPO4 content micropore non-band, [mol d-2]
  real(r8),allocatable ::  ZFE2P(:,:,:)                       !soil aqueous FeH2PO4 content micropore non-band, [mol d-2]
  real(r8),allocatable ::  ZCA0P(:,:,:)                       !soil aqueous CaPO4 content micropore non-band, [mol d-2]
  real(r8),allocatable ::  ZCA1P(:,:,:)                       !soil aqueous CaHPO4 content micropore non-band, [mol d-2]
  real(r8),allocatable ::  ZCA2P(:,:,:)                       !soil aqueous CaH2PO4 content micropore non-band, [mol d-2]
  real(r8),allocatable ::  ZMG1P(:,:,:)                       !soil aqueous MgHPO4 content micropore non-band, [mol d-2]
  real(r8),allocatable ::  H0POB(:,:,:)                       !soil aqueous PO4 content micropore band, [mol d-2]
  real(r8),allocatable ::  H3POB(:,:,:)                       !soil aqueous H3PO4 content micropore band, [mol d-2]
  real(r8),allocatable ::  ZFE1PB(:,:,:)                      !soil aqueous FeHPO4 content micropore band, [mol d-2]
  real(r8),allocatable ::  ZFE2PB(:,:,:)                      !soil aqueous FeH2PO4 content micropore band, [mol d-2]
  real(r8),allocatable ::  ZCA0PB(:,:,:)                      !soil aqueous CaPO4 content micropore band, [mol d-2]
  real(r8),allocatable ::  ZCA1PB(:,:,:)                      !soil aqueous CaHPO4 content micropore band, [mol d-2]
  real(r8),allocatable ::  ZCA2PB(:,:,:)                      !soil aqueous CaH2PO4 content micropore band, [mol d-2]
  real(r8),allocatable ::  ZMG1PB(:,:,:)                      !soil aqueous MgHPO4 content micropore band, [mol d-2]
  real(r8),allocatable ::  XHY(:,:,:)                         !exchangeable H , [mol d-2]
  real(r8),allocatable ::  XAL(:,:,:)                         !exchangeable Al, [mol d-2]
  real(r8),allocatable ::  XCA(:,:,:)                         !exchangeable Ca, [mol d-2]
  real(r8),allocatable ::  XMG(:,:,:)                         !exchangeable Mg , [mol d-2]
  real(r8),allocatable ::  XNA(:,:,:)                         !exchangeable Na, [mol d-2]
  real(r8),allocatable ::  XKA(:,:,:)                         !exchangeable K, [mol d-2]
  real(r8),allocatable ::  XHC(:,:,:)                         !exchangeable COOH , [mol d-2]
  real(r8),allocatable ::  XALO2(:,:,:)                       !exchangeable AlOH2 , [mol d-2]
  real(r8),allocatable ::  XOH0(:,:,:)                        !exchangeable OH- non-band, [mol d-2]
  real(r8),allocatable ::  XOH1(:,:,:)                        !exchangeable OH  non-band, [mol d-2]
  real(r8),allocatable ::  XOH2(:,:,:)                        !exchangeable OH2  non-band, [mol d-2]
  real(r8),allocatable ::  XH1P(:,:,:)                        !exchangeable HPO4  non-band, [mol d-2]
  real(r8),allocatable ::  XH2P(:,:,:)                        !exchangeable H2PO4  non-band, [mol d-2]
  real(r8),allocatable ::  XOH0B(:,:,:)                       !exchangeable OH- band, [mol d-2]
  real(r8),allocatable ::  XFE(:,:,:)                         !exchangeable Fe, [mol d-2]
  real(r8),allocatable ::  XFEO2(:,:,:)                       !exchangeable Fe(OH)2, [mol d-2]
  real(r8),allocatable ::  XOH1B(:,:,:)                       !exchangeable OH  band, [mol d-2]
  real(r8),allocatable ::  XOH2B(:,:,:)                       !exchangeable OH2  band, [mol d-2]
  real(r8),allocatable ::  XH1PB(:,:,:)                       !exchangeable HPO4  band, [mol d-2]
  real(r8),allocatable ::  XH2PB(:,:,:)                       !exchangeable H2PO4  band, [mol d-2]
  real(r8),allocatable ::  PCAPD(:,:,:)                       !precipitated CaHPO4 non-band, [mol d-2]
  real(r8),allocatable ::  PCAPH(:,:,:)                       !precipitated hydroxyapatite non-band, [mol d-2]
  real(r8),allocatable ::  PALOH(:,:,:)                       !precipitated AlOH3, [mol d-2]
  real(r8),allocatable ::  PFEOH(:,:,:)                       !precipitated FeOH3, [mol d-2]
  real(r8),allocatable ::  PCACO(:,:,:)                       !precipitated CaCO3, [mol d-2]
  real(r8),allocatable ::  PCASO(:,:,:)                       !precipitated CaSO4, [mol d-2]
  real(r8),allocatable ::  PALPO(:,:,:)                       !precipitated AlPO4 non-band, [mol d-2]
  real(r8),allocatable ::  PFEPO(:,:,:)                       !precipitated FePO4 non-band, [mol d-2]
  real(r8),allocatable ::  PCAPM(:,:,:)                       !precipitated CaH2PO4 non-band, [mol d-2]
  real(r8),allocatable ::  PALPB(:,:,:)                       !precipitated AlPO4 band, [mol d-2]
  real(r8),allocatable ::  PFEPB(:,:,:)                       !precipitated FePO4 band, [mol d-2]
  real(r8),allocatable ::  PCPDB(:,:,:)                       !precipitated CaHPO4 band, [mol d-2]
  real(r8),allocatable ::  PCPMB(:,:,:)                       !precipitated CaH2PO4 , [mol d-2]
  real(r8),allocatable ::  ECND(:,:,:)                        !electrical conductivity , [dS m-1]
  real(r8),allocatable ::  CSTR(:,:,:)                        !solution ion strength, [mol m-3]
  real(r8),allocatable ::  CION(:,:,:)                        !solution ion concentratiom, [mol m-3]
  real(r8),allocatable ::  XCEC(:,:,:)                        !cation exchange capacity, [mol d-2]
  real(r8),allocatable ::  XAEC(:,:,:)                        !anion exchange capacity, [mol d-2]
  real(r8),allocatable ::  ZALH(:,:,:)                        !aqueous Al content macropore, [mol d-2]
  real(r8),allocatable ::  ZFEH(:,:,:)                        !aqueous Fe content macropore, [mol d-2]
  real(r8),allocatable ::  ZHYH(:,:,:)                        !aqueous H content macropore, [mol d-2]
  real(r8),allocatable ::  ZCCH(:,:,:)                        !aqueous CO3 content macropore, [mol d-2]
  real(r8),allocatable ::  ZMAH(:,:,:)                        !aqueous Mg content macropore, [mol d-2]
  real(r8),allocatable ::  ZNAH(:,:,:)                        !aqueous Na content macropore, [mol d-2]
  real(r8),allocatable ::  ZKAH(:,:,:)                        !aqueous K content macropore, [mol d-2]
  real(r8),allocatable ::  ZOHH(:,:,:)                        !aqueous OH content macropore, [mol d-2]
  real(r8),allocatable ::  ZSO4H(:,:,:)                       !aqueous SO4 content macropore, [mol d-2]
  real(r8),allocatable ::  ZCLH(:,:,:)                        !aqueous Cl content macropore, [mol d-2]
  real(r8),allocatable ::  ZCO3H(:,:,:)                       !aqueous CO3 content macropore, [mol d-2]
  real(r8),allocatable ::  ZHCO3H(:,:,:)                      !aqueous HCO3 content macropore, [mol d-2]
  real(r8),allocatable ::  ZALO1H(:,:,:)                      !aqueous AlOH content macropore, [mol d-2]
  real(r8),allocatable ::  ZALO2H(:,:,:)                      !aqueous AlOH2 content macropore, [mol d-2]
  real(r8),allocatable ::  ZALO3H(:,:,:)                      !aqueous AlOH3 content macropore, [mol d-2]
  real(r8),allocatable ::  ZALO4H(:,:,:)                      !aqueous AlOH4 content macropore, [mol d-2]
  real(r8),allocatable ::  ZALSH(:,:,:)                       !aqueous AlSO4 content macropore, [mol d-2]
  real(r8),allocatable ::  ZFEO1H(:,:,:)                      !aqueous FeOH content macropore, [mol d-2]
  real(r8),allocatable ::  ZFEO2H(:,:,:)                      !aqueous FeOH2 content macropore, [mol d-2]
  real(r8),allocatable ::  ZFEO3H(:,:,:)                      !aqueous FeOH3 content macropore, [mol d-2]
  real(r8),allocatable ::  ZFEO4H(:,:,:)                      !aqueous FeOH4 content macropore, [mol d-2]
  real(r8),allocatable ::  ZFESH(:,:,:)                       !aqueous FeSO4 content macropore, [mol d-2]
  real(r8),allocatable ::  ZCAOH(:,:,:)                       !aqueous CaOH content macropore, [mol d-2]
  real(r8),allocatable ::  ZCACH(:,:,:)                       !aqueous CaCO3 content macropore, [mol d-2]
  real(r8),allocatable ::  ZCAHH(:,:,:)                       !aqueous CaHCO3 content macropore, [mol d-2]
  real(r8),allocatable ::  ZCASH(:,:,:)                       !aqueous CaSO4 content macropore, [mol d-2]
  real(r8),allocatable ::  PCPHB(:,:,:)                       !precipitated hydroxyapatite band, [mol d-2]
  real(r8),allocatable ::  ZMGOH(:,:,:)                       !aqueous MgOH content macropore, [mol d-2]
  real(r8),allocatable ::  ZMGCH(:,:,:)                       !aqueous MgCO3 content macropore, [mol d-2]
  real(r8),allocatable ::  ZMGHH(:,:,:)                       !aqueous MgHCO3 content macropore, [mol d-2]
  real(r8),allocatable ::  ZMGSH(:,:,:)                       !aqueous MgSO4 content macropore, [mol d-2]
  real(r8),allocatable ::  ZNACH(:,:,:)                       !aqueous NaCO3 content macropore, [mol d-2]
  real(r8),allocatable ::  ZNASH(:,:,:)                       !aqueous NaSO4 content macropore, [mol d-2]
  real(r8),allocatable ::  ZKASH(:,:,:)                       !aqueous KSO4 content macropore, [mol d-2]
  real(r8),allocatable ::  H0PO4H(:,:,:)                      !soil aqueous PO4 content non-band macropore , [mol d-2]
  real(r8),allocatable ::  H3PO4H(:,:,:)                      !soil aqueous H3PO4 content non-band macropore, [mol d-2]
  real(r8),allocatable ::  ZFE1PH(:,:,:)                      !soil aqueous FeHPO4 content non-band macropore, [mol d-2]
  real(r8),allocatable ::  ZFE2PH(:,:,:)                      !soil aqueous FeH2PO4 content non-band macropore, [mol d-2]
  real(r8),allocatable ::  ZCA0PH(:,:,:)                      !soil aqueous CaPO4 content non-band macropore, [mol d-2]
  real(r8),allocatable ::  ZCA1PH(:,:,:)                      !soil aqueous CaHPO4 content non-band macropore, [mol d-2]
  real(r8),allocatable ::  ZCA2PH(:,:,:)                      !soil aqueous CaH2PO4 content non-band macropore, [mol d-2]
  real(r8),allocatable ::  ZMG1PH(:,:,:)                      !soil aqueous MgHPO4 content non-band macropore, [mol d-2]
  real(r8),allocatable ::  H0POBH(:,:,:)                      !soil aqueous PO4 content band macropore , [mol d-2]
  real(r8),allocatable ::  H3POBH(:,:,:)                      !soil aqueous H3PO4 content band macropore, [mol d-2]
  real(r8),allocatable ::  ZFE1BH(:,:,:)                      !soil aqueous FeHPO4 content band macropore, [mol d-2]
  real(r8),allocatable ::  ZFE2BH(:,:,:)                      !soil aqueous FeH2PO4 content band macropore, [mol d-2]
  real(r8),allocatable ::  ZCA0BH(:,:,:)                      !soil aqueous CaPO4 content band macropore, [mol d-2]
  real(r8),allocatable ::  ZCA1BH(:,:,:)                      !soil aqueous CaHPO4 content band macropore, [mol d-2]
  real(r8),allocatable ::  ZCA2BH(:,:,:)                      !soil aqueous CaH2PO4 content band macropore, [mol d-2]
  real(r8),allocatable ::  ZMG1BH(:,:,:)                      !soil aqueous MgHPO4 content band macropore, [mol d-2]
  real(r8),allocatable ::  XALFLS(:,:,:,:)                    !total Al in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFEFLS(:,:,:,:)                    !total Fe in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XHYFLS(:,:,:,:)                    !total H in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCAFLS(:,:,:,:)                    !total Ca in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XMGFLS(:,:,:,:)                    !total Mg in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XNAFLS(:,:,:,:)                    !total Na in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XKAFLS(:,:,:,:)                    !total K in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XOHFLS(:,:,:,:)                    !total OH in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XSOFLS(:,:,:,:)                    !total SO4 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCLFLS(:,:,:,:)                    !total Cl in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XC3FLS(:,:,:,:)                    !total CO3 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XHCFLS(:,:,:,:)                    !total HCO3 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XAL1FS(:,:,:,:)                    !total AlOH in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XAL2FS(:,:,:,:)                    !total AlOH2 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XAL3FS(:,:,:,:)                    !total AlOH3 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XAL4FS(:,:,:,:)                    !total AlOH4 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XALSFS(:,:,:,:)                    !total AlSO4 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFE1FS(:,:,:,:)                    !total FeOH in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFE2FS(:,:,:,:)                    !total FeOH2 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFE3FS(:,:,:,:)                    !total FeOH3 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFE4FS(:,:,:,:)                    !total FeOH4 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFESFS(:,:,:,:)                    !total FeSO4 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCAOFS(:,:,:,:)                    !total CaOH in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCACFS(:,:,:,:)                    !total CaCO3 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCAHFS(:,:,:,:)                    !total CaHCO3 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCASFS(:,:,:,:)                    !total CaSO4 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XMGOFS(:,:,:,:)                    !total MgOH in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XMGCFS(:,:,:,:)                    !total MgCO3 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XMGHFS(:,:,:,:)                    !total MgHCO3 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XMGSFS(:,:,:,:)                    !total MgSO4 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XNACFS(:,:,:,:)                    !total NaCO3 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XNASFS(:,:,:,:)                    !total NaSO4 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XKASFS(:,:,:,:)                    !total KSO4 in micropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XH0PFS(:,:,:,:)                    !total PO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XH3PFS(:,:,:,:)                    !total H3PO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XF1PFS(:,:,:,:)                    !total FeHPO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XF2PFS(:,:,:,:)                    !total FeH2PO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XC0PFS(:,:,:,:)                    !total CaPO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XC1PFS(:,:,:,:)                    !total CaHPO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XC2PFS(:,:,:,:)                    !total CaH2PO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XM1PFS(:,:,:,:)                    !total MgHPO4 in micropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XH0BFB(:,:,:,:)                    !total PO4 in micropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XH3BFB(:,:,:,:)                    !total H3PO4 in micropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XF1BFB(:,:,:,:)                    !total FeHPO4 in micropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XF2BFB(:,:,:,:)                    !total FeH2PO4 in micropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XC0BFB(:,:,:,:)                    !total CaPO4 in micropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XC1BFB(:,:,:,:)                    !total CaHPO4 in micropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XC2BFB(:,:,:,:)                    !total CaH2PO4 in micropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XM1BFB(:,:,:,:)                    !total MgHPO4 in micropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XALFHS(:,:,:,:)                    !total Al in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFEFHS(:,:,:,:)                    !total Fe in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XHYFHS(:,:,:,:)                    !total H in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCAFHS(:,:,:,:)                    !total Ca in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XMGFHS(:,:,:,:)                    !total Mg in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XNAFHS(:,:,:,:)                    !total Na in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XKAFHS(:,:,:,:)                    !total K in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XOHFHS(:,:,:,:)                    !total OH in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XSOFHS(:,:,:,:)                    !total SO4 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCLFHS(:,:,:,:)                    !total Cl in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XC3FHS(:,:,:,:)                    !total CO3 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XHCFHS(:,:,:,:)                    !total HCO3 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XAL1HS(:,:,:,:)                    !total AlOH in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XAL2HS(:,:,:,:)                    !total AlOH2 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XAL3HS(:,:,:,:)                    !total AlOH3 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XAL4HS(:,:,:,:)                    !total AlOH4 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XALSHS(:,:,:,:)                    !total AlSO4 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFE1HS(:,:,:,:)                    !total FeOH in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFE2HS(:,:,:,:)                    !total FeOH2 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFE3HS(:,:,:,:)                    !total FeOH3 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFE4HS(:,:,:,:)                    !total FeOH4 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XFESHS(:,:,:,:)                    !total FeSO4 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCAOHS(:,:,:,:)                    !total CaOH in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCACHS(:,:,:,:)                    !total CaCO3 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCAHHS(:,:,:,:)                    !total CaHCO3 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XCASHS(:,:,:,:)                    !total CaSO4 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XMGOHS(:,:,:,:)                    !total MgOH in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XMGCHS(:,:,:,:)                    !total MgCO3 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XMGHHS(:,:,:,:)                    !total MgHCO3 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XMGSHS(:,:,:,:)                    !total MgSO4 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XNACHS(:,:,:,:)                    !total NaCO3 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XNASHS(:,:,:,:)                    !total NaSO4 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XKASHS(:,:,:,:)                    !total KSO4 in macropore water flux, [mol d-2 h-1]
  real(r8),allocatable ::  XH0PHS(:,:,:,:)                    !total PO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XH3PHS(:,:,:,:)                    !total H3PO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XF1PHS(:,:,:,:)                    !total FeHPO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XF2PHS(:,:,:,:)                    !total FeH2PO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XC0PHS(:,:,:,:)                    !total CaPO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XC1PHS(:,:,:,:)                    !total CaHPO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XC2PHS(:,:,:,:)                    !total CaH2PO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XM1PHS(:,:,:,:)                    !total MgHPO4 in macropore water flux non-band, [mol d-2 h-1]
  real(r8),allocatable ::  XH0BHB(:,:,:,:)                    !total PO4 in macropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XH3BHB(:,:,:,:)                    !total H3PO4 in macropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XF1BHB(:,:,:,:)                    !total FeHPO4 in macropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XF2BHB(:,:,:,:)                    !total FeH2PO4 in macropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XC0BHB(:,:,:,:)                    !total CaPO4 in macropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XC1BHB(:,:,:,:)                    !total CaHPO4 in macropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XC2BHB(:,:,:,:)                    !total CaH2PO4 in macropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XM1BHB(:,:,:,:)                    !total MgHPO4 in macropore water flux band, [mol d-2 h-1]
  real(r8),allocatable ::  XOCFXS(:,:,:,:)                    !total DOC micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XONFXS(:,:,:,:)                    !total DON micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XOPFXS(:,:,:,:)                    !total DOP micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XOAFXS(:,:,:,:)                    !total acetate micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XCOFXS(:,:,:)                      !total CO2 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XCHFXS(:,:,:)                      !total CH4 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XOXFXS(:,:,:)                      !total O2 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XHGFXS(:,:,:)                      !total H2 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XNGFXS(:,:,:)                      !total N2 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XN2FXS(:,:,:)                      !total N2O micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XN4FXW(:,:,:)                      !total NH4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XN3FXW(:,:,:)                      !total NH3 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XNOFXW(:,:,:)                      !total NO3 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XH2PXS(:,:,:)                      !total H2PO4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XN4FXB(:,:,:)                      !total NH4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XN3FXB(:,:,:)                      !total NH3 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XNOFXB(:,:,:)                      !total NO3 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XH2BXB(:,:,:)                      !total H2PO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XNXFXS(:,:,:)                      !total NO2 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XH1PXS(:,:,:)                      !total HPO4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XH1BXB(:,:,:)                      !total HPO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XH0PXS(:,:,:)                      !total PO4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XH3PXS(:,:,:)                      !total H3PO4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XALFXS(:,:,:)                      !total Al micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XFEFXS(:,:,:)                      !total Fe micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XHYFXS(:,:,:)                      !total H micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XCAFXS(:,:,:)                      !total Ca micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XMGFXS(:,:,:)                      !total Mg micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XNAFXS(:,:,:)                      !total Na micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XKAFXS(:,:,:)                      !total K micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XOHFXS(:,:,:)                      !total OH micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XSOFXS(:,:,:)                      !total SO4 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XCLFXS(:,:,:)                      !total Cl micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XC3FXS(:,:,:)                      !total CO3 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XHCFXS(:,:,:)                      !total HCO3 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XAL1XS(:,:,:)                      !total AlOH micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XAL2XS(:,:,:)                      !total AlOH2 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XAL3XS(:,:,:)                      !total AlOH3 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XAL4XS(:,:,:)                      !total AlOH4 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XALSXS(:,:,:)                      !total AlSO4 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XFE1XS(:,:,:)                      !total FeOH micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XFE2XS(:,:,:)                      !total FeOH2 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XFE3XS(:,:,:)                      !total FeOH3 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XFE4XS(:,:,:)                      !total FeOH4 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XFESXS(:,:,:)                      !total FeSO4 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XCAOXS(:,:,:)                      !total CaOH micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XCACXS(:,:,:)                      !total CaCO3 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XCAHXS(:,:,:)                      !total CaHCO3 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XCASXS(:,:,:)                      !total CaSO4 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XMGOXS(:,:,:)                      !total MgOH micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XMGCXS(:,:,:)                      !total MgCO3 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XMGHXS(:,:,:)                      !total MgHCO3 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XMGSXS(:,:,:)                      !total MgSO4 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XNACXS(:,:,:)                      !total NaCO3 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XNASXS(:,:,:)                      !total NaSO4 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XKASXS(:,:,:)                      !total KSO4 micropore-macropore transfer, [g d-2 h-1]
  real(r8),allocatable ::  XF1PXS(:,:,:)                      !total FeHPO4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XF2PXS(:,:,:)                      !total FeH2PO4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XC0PXS(:,:,:)                      !total CaPO4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XC1PXS(:,:,:)                      !total CaHPO4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XC2PXS(:,:,:)                      !total CaH2PO4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XM1PXS(:,:,:)                      !total MgHPO4 micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),allocatable ::  XH0BXB(:,:,:)                      !total PO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XH3BXB(:,:,:)                      !total H3PO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XF1BXB(:,:,:)                      !total FeHPO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XF2BXB(:,:,:)                      !total FeH2PO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XC0BXB(:,:,:)                      !total CaPO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XC1BXB(:,:,:)                      !total CaHPO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XC2BXB(:,:,:)                      !total CaH2PO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XM1BXB(:,:,:)                      !total MgHPO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  XNXFXB(:,:,:)                      !total MgHPO4 micropore-macropore transfer band, [g d-2 h-1]
  real(r8),allocatable ::  TRN4S(:,:,:)                       !total solute NH4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRN3S(:,:,:)                       !total solute NH3 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRN4B(:,:,:)                       !total solute NH4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRNO3(:,:,:)                       !total solute NO3 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRN3B(:,:,:)                       !total solute NH3 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRNOB(:,:,:)                       !total solute NO3 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRH0P(:,:,:)                       !total solute PO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRH1P(:,:,:)                       !total solute HPO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRH2P(:,:,:)                       !total solute H2PO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRH3P(:,:,:)                       !total solute H3PO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRH0B(:,:,:)                       !total solute PO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRH1B(:,:,:)                       !total solute HPO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRH2B(:,:,:)                       !total solute H2PO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRH3B(:,:,:)                       !total solute H3PO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRAL(:,:,:)                        !total solute Al transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRFE(:,:,:)                        !total solute Fe transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRHY(:,:,:)                        !total solute H transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRCA(:,:,:)                        !total solute Ca transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRMG(:,:,:)                        !total solute Mg transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRNA(:,:,:)                        !total solute Na transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRKA(:,:,:)                        !total solute K transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TROH(:,:,:)                        !total solute OH transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRSO4(:,:,:)                       !total solute SO4 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRCO3(:,:,:)                       !total solute CO3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRHCO(:,:,:)                       !total solute HCO3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRCO2(:,:,:)                       !total solute CO2 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRH2O(:,:,:)                       !total solute H2O transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRAL1(:,:,:)                       !total solute AlOH transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRAL2(:,:,:)                       !total solute AlOH2 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRAL3(:,:,:)                       !total solute AlOH3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRAL4(:,:,:)                       !total solute AlOH4 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRALS(:,:,:)                       !total solute AlSO4 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRFE1(:,:,:)                       !total solute FeOH transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRFE2(:,:,:)                       !total solute FeOH2 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRFE3(:,:,:)                       !total solute FeOH3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRFE4(:,:,:)                       !total solute FeOH4 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRFES(:,:,:)                       !total solute FeSO4 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRCAO(:,:,:)                       !total solute CaOH transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRCAC(:,:,:)                       !total solute CaCO3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRCAH(:,:,:)                       !total solute CaHCO3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRCAS(:,:,:)                       !total solute CaSO4 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRMGO(:,:,:)                       !total solute MgOH transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRMGC(:,:,:)                       !total solute MgCO3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRMGH(:,:,:)                       !total solute MgHCO3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRMGS(:,:,:)                       !total solute MgSO4 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRNAC(:,:,:)                       !total solute NaCO3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRNAS(:,:,:)                       !total solute NaSO4 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRF1P(:,:,:)                       !total solute FeHPO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRF2P(:,:,:)                       !total solute FeH2PO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRC0P(:,:,:)                       !total solute CaPO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRC1P(:,:,:)                       !total solute CaHPO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRC2P(:,:,:)                       !total solute CaH2PO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRM1P(:,:,:)                       !total solute MgHPO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRF1B(:,:,:)                       !total solute FeHPO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRF2B(:,:,:)                       !total solute FeH2PO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRC0B(:,:,:)                       !total solute CaPO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRC1B(:,:,:)                       !total solute CaHPO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRC2B(:,:,:)                       !total solute CaH2PO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRM1B(:,:,:)                       !total solute MgHPO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRXN4(:,:,:)                       !total adsorbed NH4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRXNB(:,:,:)                       !total adsorbed NH4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRXHY(:,:,:)                       !total adsorbed H transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRXAL(:,:,:)                       !total adsorbed Al transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRXCA(:,:,:)                       !total adsorbed Ca transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRXMG(:,:,:)                       !total adsorbed Mg transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRXNA(:,:,:)                       !total adsorbed Na transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRXKA(:,:,:)                       !total adsorbed K transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRXHC(:,:,:)                       !total adsorbed COOH transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRXAL2(:,:,:)                      !total adsorbed AlOH2 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRKAS(:,:,:)                       !total solute KSO4 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRXFE(:,:,:)                       !total Fe adsorption
  real(r8),allocatable ::  TRXFE2(:,:,:)                      !total FeOH2 adsorption
  real(r8),allocatable ::  TRXH0(:,:,:)                       !total adsorbed OH- transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRXH1(:,:,:)                       !total adsorbed OH transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRXH2(:,:,:)                       !total adsorbed OH2 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRX1P(:,:,:)                       !total adsorbed HPO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRX2P(:,:,:)                       !total adsorbed H2PO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRBH0(:,:,:)                       !total adsorbed OH- transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRBH1(:,:,:)                       !total adsorbed OH transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRBH2(:,:,:)                       !total adsorbed OH2 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRB1P(:,:,:)                       !total adsorbed HPO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRB2P(:,:,:)                       !total adsorbed H2PO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRCPMB(:,:,:)                      !total precipitated apatite transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TBCO2(:,:,:)                       !total solute CO2 transformation boundary, [mol d-2 h-1]
  real(r8),allocatable ::  TBION(:,:,:)                       !total solute ion transformation boundary, [mol d-2 h-1]
  real(r8),allocatable ::  TRNO2(:,:,:)                       !total solute NO2 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRN2B(:,:,:)                       !total solute NO2 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRN3G(:,:,:)                       !total gaseous NH3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRALOH(:,:,:)                      !total precipitated AlOH3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRFEOH(:,:,:)                      !total precipitated FeOH3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRCACO(:,:,:)                      !total precipitated CaCO3 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRCASO(:,:,:)                      !total precipitated CaSO4 transformation, [mol d-2 h-1]
  real(r8),allocatable ::  TRALPO(:,:,:)                      !total precipitated AlPO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRFEPO(:,:,:)                      !total precipitated FePO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRCAPD(:,:,:)                      !total precipitated CaHPO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRCAPH(:,:,:)                      !total precipitated CaH2PO4 transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRCAPM(:,:,:)                      !total precipitated apatite transformation non-band, [mol d-2 h-1]
  real(r8),allocatable ::  TRALPB(:,:,:)                      !total precipitated AlPO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRFEPB(:,:,:)                      !total precipitated FePO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRCPDB(:,:,:)                      !total precipitated CaHPO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  TRCPHB(:,:,:)                      !total precipitated CaH2PO4 transformation band, [mol d-2 h-1]
  real(r8),allocatable ::  XCOBLS(:,:,:)                      !wet deposition of CO2, [g d-2 h-1]
  real(r8),allocatable ::  XCHBLS(:,:,:)                      !wet deposition of CH4, [g d-2 h-1]
  real(r8),allocatable ::  XOXBLS(:,:,:)                      !wet deposition of O2, [g d-2 h-1]
  real(r8),allocatable ::  XNGBLS(:,:,:)                      !wet deposition of N2, [g d-2 h-1]
  real(r8),allocatable ::  XN2BLS(:,:,:)                      !wet deposition of N2O, [g d-2 h-1]
  real(r8),allocatable ::  XN4BLW(:,:,:)                      !wet deposition of NH4, [g d-2 h-1]
  real(r8),allocatable ::  XN3BLW(:,:,:)                      !wet deposition of NH3, [g d-2 h-1]
  real(r8),allocatable ::  XNOBLW(:,:,:)                      !wet deposition of NO3, [g d-2 h-1]
  real(r8),allocatable ::  XH1PBS(:,:,:)                      !wet deposition of HPO4, [g d-2 h-1]
  real(r8),allocatable ::  XH2PBS(:,:,:)                      !wet deposition of H2PO4, [g d-2 h-1]
  real(r8),allocatable ::  XALBLS(:,:,:)                      !wet deposition of Al, [g d-2 h-1]
  real(r8),allocatable ::  XFEBLS(:,:,:)                      !wet deposition of Fe, [g d-2 h-1]
  real(r8),allocatable ::  XHYBLS(:,:,:)                      !wet deposition of H, [g d-2 h-1]
  real(r8),allocatable ::  XCABLS(:,:,:)                      !wet deposition of Ca, [g d-2 h-1]
  real(r8),allocatable ::  XMGBLS(:,:,:)                      !wet deposition of Mg, [g d-2 h-1]
  real(r8),allocatable ::  XNABLS(:,:,:)                      !wet deposition of Na, [g d-2 h-1]
  real(r8),allocatable ::  XKABLS(:,:,:)                      !wet deposition of K, [g d-2 h-1]
  real(r8),allocatable ::  XOHBLS(:,:,:)                      !wet deposition of OH, [g d-2 h-1]
  real(r8),allocatable ::  XSOBLS(:,:,:)                      !wet deposition of SO4, [g d-2 h-1]
  real(r8),allocatable ::  XCLBLS(:,:,:)                      !wet deposition of Cl, [g d-2 h-1]
  real(r8),allocatable ::  XC3BLS(:,:,:)                      !wet deposition of CO3, [g d-2 h-1]
  real(r8),allocatable ::  XHCBLS(:,:,:)                      !wet deposition of HCO3, [g d-2 h-1]
  real(r8),allocatable ::  XAL1BS(:,:,:)                      !wet deposition of AlOH, [g d-2 h-1]
  real(r8),allocatable ::  XAL2BS(:,:,:)                      !wet deposition of AlOH2, [g d-2 h-1]
  real(r8),allocatable ::  XAL3BS(:,:,:)                      !wet deposition of AlOH3, [g d-2 h-1]
  real(r8),allocatable ::  XAL4BS(:,:,:)                      !wet deposition of AlOH4, [g d-2 h-1]
  real(r8),allocatable ::  XALSBS(:,:,:)                      !wet deposition of AlSO4, [g d-2 h-1]
  real(r8),allocatable ::  XFE1BS(:,:,:)                      !wet deposition of FeOH, [g d-2 h-1]
  real(r8),allocatable ::  XFE2BS(:,:,:)                      !wet deposition of FeOH2, [g d-2 h-1]
  real(r8),allocatable ::  XFE3BS(:,:,:)                      !wet deposition of FeOH3, [g d-2 h-1]
  real(r8),allocatable ::  XFE4BS(:,:,:)                      !wet deposition of FeOH4, [g d-2 h-1]
  real(r8),allocatable ::  XFESBS(:,:,:)                      !wet deposition of FeSO4, [g d-2 h-1]
  real(r8),allocatable ::  XCAOBS(:,:,:)                      !wet deposition of CaOH, [g d-2 h-1]
  real(r8),allocatable ::  XCACBS(:,:,:)                      !wet deposition of CaCO3, [g d-2 h-1]
  real(r8),allocatable ::  XCAHBS(:,:,:)                      !wet deposition of CaHCO3, [g d-2 h-1]
  real(r8),allocatable ::  XCASBS(:,:,:)                      !wet deposition of CaSO4, [g d-2 h-1]
  real(r8),allocatable ::  XMGOBS(:,:,:)                      !wet deposition of MgOH, [g d-2 h-1]
  real(r8),allocatable ::  XMGCBS(:,:,:)                      !wet deposition of MgCO3, [g d-2 h-1]
  real(r8),allocatable ::  XHGBLS(:,:,:)                      !wet deposition of H2, [g d-2 h-1]
  real(r8),allocatable ::  XMGHBS(:,:,:)                      !wet deposition of MgHCO3, [g d-2 h-1]
  real(r8),allocatable ::  XMGSBS(:,:,:)                      !wet deposition of MgSO4, [g d-2 h-1]
  real(r8),allocatable ::  XNACBS(:,:,:)                      !wet deposition of NaCO3, [g d-2 h-1]
  real(r8),allocatable ::  XNASBS(:,:,:)                      !wet deposition of NaSO4, [g d-2 h-1]
  real(r8),allocatable ::  XKASBS(:,:,:)                      !wet deposition of KSO4, [g d-2 h-1]
  real(r8),allocatable ::  XH0PBS(:,:,:)                      !wet deposition of PO4, [g d-2 h-1]
  real(r8),allocatable ::  XH3PBS(:,:,:)                      !wet deposition of H3PO4, [g d-2 h-1]
  real(r8),allocatable ::  XF1PBS(:,:,:)                      !wet deposition of FeHPO4, [g d-2 h-1]
  real(r8),allocatable ::  XF2PBS(:,:,:)                      !wet deposition of FeH2PO4, [g d-2 h-1]
  real(r8),allocatable ::  XC0PBS(:,:,:)                      !wet deposition of CaPO4, [g d-2 h-1]
  real(r8),allocatable ::  XC1PBS(:,:,:)                      !wet deposition of CHPO4, [g d-2 h-1]
  real(r8),allocatable ::  XC2PBS(:,:,:)                      !wet deposition of CaH2PO4, [g d-2 h-1]
  real(r8),allocatable ::  XM1PBS(:,:,:)                      !wet deposition of MgHPO4, [g d-2 h-1]
  private :: InitAllocate
  contains

  subroutine InitAquaChem
  implicit none

  call InitAllocate

  end subroutine InitAquaChem

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(CAL(JZ,JY,JX));      CAL=0._r8
  allocate(CFE(JZ,JY,JX));      CFE=0._r8
  allocate(CCA(JZ,JY,JX));      CCA=0._r8
  allocate(CMG(JZ,JY,JX));      CMG=0._r8
  allocate(CNA(JZ,JY,JX));      CNA=0._r8
  allocate(CKA(JZ,JY,JX));      CKA=0._r8
  allocate(CSO4(JZ,JY,JX));     CSO4=0._r8
  allocate(CCL(JZ,JY,JX));      CCL=0._r8
  allocate(CALOH(JZ,JY,JX));    CALOH=0._r8
  allocate(CFEOH(JZ,JY,JX));    CFEOH=0._r8
  allocate(CCACO(JZ,JY,JX));    CCACO=0._r8
  allocate(CCASO(JZ,JY,JX));    CCASO=0._r8
  allocate(CALPO(JZ,JY,JX));    CALPO=0._r8
  allocate(CFEPO(JZ,JY,JX));    CFEPO=0._r8
  allocate(CCAPD(JZ,JY,JX));    CCAPD=0._r8
  allocate(CCAPH(JZ,JY,JX));    CCAPH=0._r8
  allocate(GKC4(JZ,JY,JX));     GKC4=0._r8
  allocate(GKCH(JZ,JY,JX));     GKCH=0._r8
  allocate(GKCA(JZ,JY,JX));     GKCA=0._r8
  allocate(GKCM(JZ,JY,JX));     GKCM=0._r8
  allocate(GKCN(JZ,JY,JX));     GKCN=0._r8
  allocate(GKCK(JZ,JY,JX));     GKCK=0._r8
  allocate(ZCO3(0:JZ,JY,JX));   ZCO3=0._r8
  allocate(ZHCO3(0:JZ,JY,JX));  ZHCO3=0._r8
  allocate(XN4(0:JZ,JY,JX));    XN4=0._r8
  allocate(XNB(0:JZ,JY,JX));    XNB=0._r8
  allocate(ZAL(0:JZ,JY,JX));    ZAL=0._r8
  allocate(ZFE(0:JZ,JY,JX));    ZFE=0._r8
  allocate(ZHY(0:JZ,JY,JX));    ZHY=0._r8
  allocate(ZCA(0:JZ,JY,JX));    ZCA=0._r8
  allocate(ZMG(0:JZ,JY,JX));    ZMG=0._r8
  allocate(ZNA(0:JZ,JY,JX));    ZNA=0._r8
  allocate(ZKA(0:JZ,JY,JX));    ZKA=0._r8
  allocate(ZOH(0:JZ,JY,JX));    ZOH=0._r8
  allocate(ZSO4(0:JZ,JY,JX));   ZSO4=0._r8
  allocate(ZCL(0:JZ,JY,JX));    ZCL=0._r8
  allocate(ZALOH1(0:JZ,JY,JX)); ZALOH1=0._r8
  allocate(ZALOH2(0:JZ,JY,JX)); ZALOH2=0._r8
  allocate(ZALOH3(0:JZ,JY,JX)); ZALOH3=0._r8
  allocate(ZALOH4(0:JZ,JY,JX)); ZALOH4=0._r8
  allocate(ZALS(0:JZ,JY,JX));   ZALS=0._r8
  allocate(ZFEOH1(0:JZ,JY,JX)); ZFEOH1=0._r8
  allocate(ZFEOH2(0:JZ,JY,JX)); ZFEOH2=0._r8
  allocate(ZFEOH3(0:JZ,JY,JX)); ZFEOH3=0._r8
  allocate(ZFEOH4(0:JZ,JY,JX)); ZFEOH4=0._r8
  allocate(ZFES(0:JZ,JY,JX));   ZFES=0._r8
  allocate(ZCAO(0:JZ,JY,JX));   ZCAO=0._r8
  allocate(ZCAC(0:JZ,JY,JX));   ZCAC=0._r8
  allocate(ZCAH(0:JZ,JY,JX));   ZCAH=0._r8
  allocate(ZCAS(0:JZ,JY,JX));   ZCAS=0._r8
  allocate(ZMGO(0:JZ,JY,JX));   ZMGO=0._r8
  allocate(ZMGC(0:JZ,JY,JX));   ZMGC=0._r8
  allocate(ZMGH(0:JZ,JY,JX));   ZMGH=0._r8
  allocate(ZMGS(0:JZ,JY,JX));   ZMGS=0._r8
  allocate(ZNAC(0:JZ,JY,JX));   ZNAC=0._r8
  allocate(ZNAS(0:JZ,JY,JX));   ZNAS=0._r8
  allocate(ZKAS(0:JZ,JY,JX));   ZKAS=0._r8
  allocate(H0PO4(0:JZ,JY,JX));  H0PO4=0._r8
  allocate(H3PO4(0:JZ,JY,JX));  H3PO4=0._r8
  allocate(ZFE1P(0:JZ,JY,JX));  ZFE1P=0._r8
  allocate(ZFE2P(0:JZ,JY,JX));  ZFE2P=0._r8
  allocate(ZCA0P(0:JZ,JY,JX));  ZCA0P=0._r8
  allocate(ZCA1P(0:JZ,JY,JX));  ZCA1P=0._r8
  allocate(ZCA2P(0:JZ,JY,JX));  ZCA2P=0._r8
  allocate(ZMG1P(0:JZ,JY,JX));  ZMG1P=0._r8
  allocate(H0POB(JZ,JY,JX));    H0POB=0._r8
  allocate(H3POB(JZ,JY,JX));    H3POB=0._r8
  allocate(ZFE1PB(JZ,JY,JX));   ZFE1PB=0._r8
  allocate(ZFE2PB(JZ,JY,JX));   ZFE2PB=0._r8
  allocate(ZCA0PB(JZ,JY,JX));   ZCA0PB=0._r8
  allocate(ZCA1PB(JZ,JY,JX));   ZCA1PB=0._r8
  allocate(ZCA2PB(JZ,JY,JX));   ZCA2PB=0._r8
  allocate(ZMG1PB(JZ,JY,JX));   ZMG1PB=0._r8
  allocate(XHY(JZ,JY,JX));      XHY=0._r8
  allocate(XAL(JZ,JY,JX));      XAL=0._r8
  allocate(XCA(JZ,JY,JX));      XCA=0._r8
  allocate(XMG(JZ,JY,JX));      XMG=0._r8
  allocate(XNA(JZ,JY,JX));      XNA=0._r8
  allocate(XKA(JZ,JY,JX));      XKA=0._r8
  allocate(XHC(JZ,JY,JX));      XHC=0._r8
  allocate(XALO2(JZ,JY,JX));    XALO2=0._r8
  allocate(XOH0(0:JZ,JY,JX));   XOH0=0._r8
  allocate(XOH1(0:JZ,JY,JX));   XOH1=0._r8
  allocate(XOH2(0:JZ,JY,JX));   XOH2=0._r8
  allocate(XH1P(0:JZ,JY,JX));   XH1P=0._r8
  allocate(XH2P(0:JZ,JY,JX));   XH2P=0._r8
  allocate(XOH0B(0:JZ,JY,JX));  XOH0B=0._r8
  allocate(XFE(JZ,JY,JX));      XFE=0._r8
  allocate(XFEO2(JZ,JY,JX));    XFEO2=0._r8
  allocate(XOH1B(0:JZ,JY,JX));  XOH1B=0._r8
  allocate(XOH2B(0:JZ,JY,JX));  XOH2B=0._r8
  allocate(XH1PB(0:JZ,JY,JX));  XH1PB=0._r8
  allocate(XH2PB(0:JZ,JY,JX));  XH2PB=0._r8
  allocate(PCAPD(0:JZ,JY,JX));  PCAPD=0._r8
  allocate(PCAPH(0:JZ,JY,JX));  PCAPH=0._r8
  allocate(PALOH(JZ,JY,JX));    PALOH=0._r8
  allocate(PFEOH(JZ,JY,JX));    PFEOH=0._r8
  allocate(PCACO(JZ,JY,JX));    PCACO=0._r8
  allocate(PCASO(JZ,JY,JX));    PCASO=0._r8
  allocate(PALPO(0:JZ,JY,JX));  PALPO=0._r8
  allocate(PFEPO(0:JZ,JY,JX));  PFEPO=0._r8
  allocate(PCAPM(0:JZ,JY,JX));  PCAPM=0._r8
  allocate(PALPB(0:JZ,JY,JX));  PALPB=0._r8
  allocate(PFEPB(0:JZ,JY,JX));  PFEPB=0._r8
  allocate(PCPDB(0:JZ,JY,JX));  PCPDB=0._r8
  allocate(PCPMB(0:JZ,JY,JX));  PCPMB=0._r8
  allocate(ECND(JZ,JY,JX));     ECND=0._r8
  allocate(CSTR(JZ,JY,JX));     CSTR=0._r8
  allocate(CION(JZ,JY,JX));     CION=0._r8
  allocate(XCEC(JZ,JY,JX));     XCEC=0._r8
  allocate(XAEC(JZ,JY,JX));     XAEC=0._r8
  allocate(ZALH(JZ,JY,JX));     ZALH=0._r8
  allocate(ZFEH(JZ,JY,JX));     ZFEH=0._r8
  allocate(ZHYH(JZ,JY,JX));     ZHYH=0._r8
  allocate(ZCCH(JZ,JY,JX));     ZCCH=0._r8
  allocate(ZMAH(JZ,JY,JX));     ZMAH=0._r8
  allocate(ZNAH(JZ,JY,JX));     ZNAH=0._r8
  allocate(ZKAH(JZ,JY,JX));     ZKAH=0._r8
  allocate(ZOHH(JZ,JY,JX));     ZOHH=0._r8
  allocate(ZSO4H(JZ,JY,JX));    ZSO4H=0._r8
  allocate(ZCLH(JZ,JY,JX));     ZCLH=0._r8
  allocate(ZCO3H(JZ,JY,JX));    ZCO3H=0._r8
  allocate(ZHCO3H(JZ,JY,JX));   ZHCO3H=0._r8
  allocate(ZALO1H(JZ,JY,JX));   ZALO1H=0._r8
  allocate(ZALO2H(JZ,JY,JX));   ZALO2H=0._r8
  allocate(ZALO3H(JZ,JY,JX));   ZALO3H=0._r8
  allocate(ZALO4H(JZ,JY,JX));   ZALO4H=0._r8
  allocate(ZALSH(JZ,JY,JX));    ZALSH=0._r8
  allocate(ZFEO1H(JZ,JY,JX));   ZFEO1H=0._r8
  allocate(ZFEO2H(JZ,JY,JX));   ZFEO2H=0._r8
  allocate(ZFEO3H(JZ,JY,JX));   ZFEO3H=0._r8
  allocate(ZFEO4H(JZ,JY,JX));   ZFEO4H=0._r8
  allocate(ZFESH(JZ,JY,JX));    ZFESH=0._r8
  allocate(ZCAOH(JZ,JY,JX));    ZCAOH=0._r8
  allocate(ZCACH(JZ,JY,JX));    ZCACH=0._r8
  allocate(ZCAHH(JZ,JY,JX));    ZCAHH=0._r8
  allocate(ZCASH(JZ,JY,JX));    ZCASH=0._r8
  allocate(PCPHB(0:JZ,JY,JX));  PCPHB=0._r8
  allocate(ZMGOH(JZ,JY,JX));    ZMGOH=0._r8
  allocate(ZMGCH(JZ,JY,JX));    ZMGCH=0._r8
  allocate(ZMGHH(JZ,JY,JX));    ZMGHH=0._r8
  allocate(ZMGSH(JZ,JY,JX));    ZMGSH=0._r8
  allocate(ZNACH(JZ,JY,JX));    ZNACH=0._r8
  allocate(ZNASH(JZ,JY,JX));    ZNASH=0._r8
  allocate(ZKASH(JZ,JY,JX));    ZKASH=0._r8
  allocate(H0PO4H(JZ,JY,JX));   H0PO4H=0._r8
  allocate(H3PO4H(JZ,JY,JX));   H3PO4H=0._r8
  allocate(ZFE1PH(JZ,JY,JX));   ZFE1PH=0._r8
  allocate(ZFE2PH(JZ,JY,JX));   ZFE2PH=0._r8
  allocate(ZCA0PH(JZ,JY,JX));   ZCA0PH=0._r8
  allocate(ZCA1PH(JZ,JY,JX));   ZCA1PH=0._r8
  allocate(ZCA2PH(JZ,JY,JX));   ZCA2PH=0._r8
  allocate(ZMG1PH(JZ,JY,JX));   ZMG1PH=0._r8
  allocate(H0POBH(JZ,JY,JX));   H0POBH=0._r8
  allocate(H3POBH(JZ,JY,JX));   H3POBH=0._r8
  allocate(ZFE1BH(JZ,JY,JX));   ZFE1BH=0._r8
  allocate(ZFE2BH(JZ,JY,JX));   ZFE2BH=0._r8
  allocate(ZCA0BH(JZ,JY,JX));   ZCA0BH=0._r8
  allocate(ZCA1BH(JZ,JY,JX));   ZCA1BH=0._r8
  allocate(ZCA2BH(JZ,JY,JX));   ZCA2BH=0._r8
  allocate(ZMG1BH(JZ,JY,JX));   ZMG1BH=0._r8
  allocate(XALFLS(3,0:JD,JV,JH));XALFLS=0._r8
  allocate(XFEFLS(3,0:JD,JV,JH));XFEFLS=0._r8
  allocate(XHYFLS(3,0:JD,JV,JH));XHYFLS=0._r8
  allocate(XCAFLS(3,0:JD,JV,JH));XCAFLS=0._r8
  allocate(XMGFLS(3,0:JD,JV,JH));XMGFLS=0._r8
  allocate(XNAFLS(3,0:JD,JV,JH));XNAFLS=0._r8
  allocate(XKAFLS(3,0:JD,JV,JH));XKAFLS=0._r8
  allocate(XOHFLS(3,0:JD,JV,JH));XOHFLS=0._r8
  allocate(XSOFLS(3,0:JD,JV,JH));XSOFLS=0._r8
  allocate(XCLFLS(3,0:JD,JV,JH));XCLFLS=0._r8
  allocate(XC3FLS(3,0:JD,JV,JH));XC3FLS=0._r8
  allocate(XHCFLS(3,0:JD,JV,JH));XHCFLS=0._r8
  allocate(XAL1FS(3,0:JD,JV,JH));XAL1FS=0._r8
  allocate(XAL2FS(3,0:JD,JV,JH));XAL2FS=0._r8
  allocate(XAL3FS(3,0:JD,JV,JH));XAL3FS=0._r8
  allocate(XAL4FS(3,0:JD,JV,JH));XAL4FS=0._r8
  allocate(XALSFS(3,0:JD,JV,JH));XALSFS=0._r8
  allocate(XFE1FS(3,0:JD,JV,JH));XFE1FS=0._r8
  allocate(XFE2FS(3,0:JD,JV,JH));XFE2FS=0._r8
  allocate(XFE3FS(3,0:JD,JV,JH));XFE3FS=0._r8
  allocate(XFE4FS(3,0:JD,JV,JH));XFE4FS=0._r8
  allocate(XFESFS(3,0:JD,JV,JH));XFESFS=0._r8
  allocate(XCAOFS(3,0:JD,JV,JH));XCAOFS=0._r8
  allocate(XCACFS(3,0:JD,JV,JH));XCACFS=0._r8
  allocate(XCAHFS(3,0:JD,JV,JH));XCAHFS=0._r8
  allocate(XCASFS(3,0:JD,JV,JH));XCASFS=0._r8
  allocate(XMGOFS(3,0:JD,JV,JH));XMGOFS=0._r8
  allocate(XMGCFS(3,0:JD,JV,JH));XMGCFS=0._r8
  allocate(XMGHFS(3,0:JD,JV,JH));XMGHFS=0._r8
  allocate(XMGSFS(3,0:JD,JV,JH));XMGSFS=0._r8
  allocate(XNACFS(3,0:JD,JV,JH));XNACFS=0._r8
  allocate(XNASFS(3,0:JD,JV,JH));XNASFS=0._r8
  allocate(XKASFS(3,0:JD,JV,JH));XKASFS=0._r8
  allocate(XH0PFS(3,0:JD,JV,JH));XH0PFS=0._r8
  allocate(XH3PFS(3,0:JD,JV,JH));XH3PFS=0._r8
  allocate(XF1PFS(3,0:JD,JV,JH));XF1PFS=0._r8
  allocate(XF2PFS(3,0:JD,JV,JH));XF2PFS=0._r8
  allocate(XC0PFS(3,0:JD,JV,JH));XC0PFS=0._r8
  allocate(XC1PFS(3,0:JD,JV,JH));XC1PFS=0._r8
  allocate(XC2PFS(3,0:JD,JV,JH));XC2PFS=0._r8
  allocate(XM1PFS(3,0:JD,JV,JH));XM1PFS=0._r8
  allocate(XH0BFB(3,0:JD,JV,JH));XH0BFB=0._r8
  allocate(XH3BFB(3,0:JD,JV,JH));XH3BFB=0._r8
  allocate(XF1BFB(3,0:JD,JV,JH));XF1BFB=0._r8
  allocate(XF2BFB(3,0:JD,JV,JH));XF2BFB=0._r8
  allocate(XC0BFB(3,0:JD,JV,JH));XC0BFB=0._r8
  allocate(XC1BFB(3,0:JD,JV,JH));XC1BFB=0._r8
  allocate(XC2BFB(3,0:JD,JV,JH));XC2BFB=0._r8
  allocate(XM1BFB(3,0:JD,JV,JH));XM1BFB=0._r8
  allocate(XALFHS(3,JD,JV,JH)); XALFHS=0._r8
  allocate(XFEFHS(3,JD,JV,JH)); XFEFHS=0._r8
  allocate(XHYFHS(3,JD,JV,JH)); XHYFHS=0._r8
  allocate(XCAFHS(3,JD,JV,JH)); XCAFHS=0._r8
  allocate(XMGFHS(3,JD,JV,JH)); XMGFHS=0._r8
  allocate(XNAFHS(3,JD,JV,JH)); XNAFHS=0._r8
  allocate(XKAFHS(3,JD,JV,JH)); XKAFHS=0._r8
  allocate(XOHFHS(3,JD,JV,JH)); XOHFHS=0._r8
  allocate(XSOFHS(3,JD,JV,JH)); XSOFHS=0._r8
  allocate(XCLFHS(3,JD,JV,JH)); XCLFHS=0._r8
  allocate(XC3FHS(3,JD,JV,JH)); XC3FHS=0._r8
  allocate(XHCFHS(3,JD,JV,JH)); XHCFHS=0._r8
  allocate(XAL1HS(3,JD,JV,JH)); XAL1HS=0._r8
  allocate(XAL2HS(3,JD,JV,JH)); XAL2HS=0._r8
  allocate(XAL3HS(3,JD,JV,JH)); XAL3HS=0._r8
  allocate(XAL4HS(3,JD,JV,JH)); XAL4HS=0._r8
  allocate(XALSHS(3,JD,JV,JH)); XALSHS=0._r8
  allocate(XFE1HS(3,JD,JV,JH)); XFE1HS=0._r8
  allocate(XFE2HS(3,JD,JV,JH)); XFE2HS=0._r8
  allocate(XFE3HS(3,JD,JV,JH)); XFE3HS=0._r8
  allocate(XFE4HS(3,JD,JV,JH)); XFE4HS=0._r8
  allocate(XFESHS(3,JD,JV,JH)); XFESHS=0._r8
  allocate(XCAOHS(3,JD,JV,JH)); XCAOHS=0._r8
  allocate(XCACHS(3,JD,JV,JH)); XCACHS=0._r8
  allocate(XCAHHS(3,JD,JV,JH)); XCAHHS=0._r8
  allocate(XCASHS(3,JD,JV,JH)); XCASHS=0._r8
  allocate(XMGOHS(3,JD,JV,JH)); XMGOHS=0._r8
  allocate(XMGCHS(3,JD,JV,JH)); XMGCHS=0._r8
  allocate(XMGHHS(3,JD,JV,JH)); XMGHHS=0._r8
  allocate(XMGSHS(3,JD,JV,JH)); XMGSHS=0._r8
  allocate(XNACHS(3,JD,JV,JH)); XNACHS=0._r8
  allocate(XNASHS(3,JD,JV,JH)); XNASHS=0._r8
  allocate(XKASHS(3,JD,JV,JH)); XKASHS=0._r8
  allocate(XH0PHS(3,JD,JV,JH)); XH0PHS=0._r8
  allocate(XH3PHS(3,JD,JV,JH)); XH3PHS=0._r8
  allocate(XF1PHS(3,JD,JV,JH)); XF1PHS=0._r8
  allocate(XF2PHS(3,JD,JV,JH)); XF2PHS=0._r8
  allocate(XC0PHS(3,JD,JV,JH)); XC0PHS=0._r8
  allocate(XC1PHS(3,JD,JV,JH)); XC1PHS=0._r8
  allocate(XC2PHS(3,JD,JV,JH)); XC2PHS=0._r8
  allocate(XM1PHS(3,JD,JV,JH)); XM1PHS=0._r8
  allocate(XH0BHB(3,JD,JV,JH)); XH0BHB=0._r8
  allocate(XH3BHB(3,JD,JV,JH)); XH3BHB=0._r8
  allocate(XF1BHB(3,JD,JV,JH)); XF1BHB=0._r8
  allocate(XF2BHB(3,JD,JV,JH)); XF2BHB=0._r8
  allocate(XC0BHB(3,JD,JV,JH)); XC0BHB=0._r8
  allocate(XC1BHB(3,JD,JV,JH)); XC1BHB=0._r8
  allocate(XC2BHB(3,JD,JV,JH)); XC2BHB=0._r8
  allocate(XM1BHB(3,JD,JV,JH)); XM1BHB=0._r8
  allocate(XOCFXS(0:jcplx1,JZ,JY,JX));XOCFXS=0._r8
  allocate(XONFXS(0:jcplx1,JZ,JY,JX));XONFXS=0._r8
  allocate(XOPFXS(0:jcplx1,JZ,JY,JX));XOPFXS=0._r8
  allocate(XOAFXS(0:jcplx1,JZ,JY,JX));XOAFXS=0._r8
  allocate(XCOFXS(JZ,JY,JX));   XCOFXS=0._r8
  allocate(XCHFXS(JZ,JY,JX));   XCHFXS=0._r8
  allocate(XOXFXS(JZ,JY,JX));   XOXFXS=0._r8
  allocate(XHGFXS(JZ,JY,JX));   XHGFXS=0._r8
  allocate(XNGFXS(JZ,JY,JX));   XNGFXS=0._r8
  allocate(XN2FXS(JZ,JY,JX));   XN2FXS=0._r8
  allocate(XN4FXW(JZ,JY,JX));   XN4FXW=0._r8
  allocate(XN3FXW(JZ,JY,JX));   XN3FXW=0._r8
  allocate(XNOFXW(JZ,JY,JX));   XNOFXW=0._r8
  allocate(XH2PXS(JZ,JY,JX));   XH2PXS=0._r8
  allocate(XN4FXB(JZ,JY,JX));   XN4FXB=0._r8
  allocate(XN3FXB(JZ,JY,JX));   XN3FXB=0._r8
  allocate(XNOFXB(JZ,JY,JX));   XNOFXB=0._r8
  allocate(XH2BXB(JZ,JY,JX));   XH2BXB=0._r8
  allocate(XNXFXS(JZ,JY,JX));   XNXFXS=0._r8
  allocate(XH1PXS(JZ,JY,JX));   XH1PXS=0._r8
  allocate(XH1BXB(JZ,JY,JX));   XH1BXB=0._r8
  allocate(XH0PXS(JZ,JY,JX));   XH0PXS=0._r8
  allocate(XH3PXS(JZ,JY,JX));   XH3PXS=0._r8
  allocate(XALFXS(JZ,JY,JX));   XALFXS=0._r8
  allocate(XFEFXS(JZ,JY,JX));   XFEFXS=0._r8
  allocate(XHYFXS(JZ,JY,JX));   XHYFXS=0._r8
  allocate(XCAFXS(JZ,JY,JX));   XCAFXS=0._r8
  allocate(XMGFXS(JZ,JY,JX));   XMGFXS=0._r8
  allocate(XNAFXS(JZ,JY,JX));   XNAFXS=0._r8
  allocate(XKAFXS(JZ,JY,JX));   XKAFXS=0._r8
  allocate(XOHFXS(JZ,JY,JX));   XOHFXS=0._r8
  allocate(XSOFXS(JZ,JY,JX));   XSOFXS=0._r8
  allocate(XCLFXS(JZ,JY,JX));   XCLFXS=0._r8
  allocate(XC3FXS(JZ,JY,JX));   XC3FXS=0._r8
  allocate(XHCFXS(JZ,JY,JX));   XHCFXS=0._r8
  allocate(XAL1XS(JZ,JY,JX));   XAL1XS=0._r8
  allocate(XAL2XS(JZ,JY,JX));   XAL2XS=0._r8
  allocate(XAL3XS(JZ,JY,JX));   XAL3XS=0._r8
  allocate(XAL4XS(JZ,JY,JX));   XAL4XS=0._r8
  allocate(XALSXS(JZ,JY,JX));   XALSXS=0._r8
  allocate(XFE1XS(JZ,JY,JX));   XFE1XS=0._r8
  allocate(XFE2XS(JZ,JY,JX));   XFE2XS=0._r8
  allocate(XFE3XS(JZ,JY,JX));   XFE3XS=0._r8
  allocate(XFE4XS(JZ,JY,JX));   XFE4XS=0._r8
  allocate(XFESXS(JZ,JY,JX));   XFESXS=0._r8
  allocate(XCAOXS(JZ,JY,JX));   XCAOXS=0._r8
  allocate(XCACXS(JZ,JY,JX));   XCACXS=0._r8
  allocate(XCAHXS(JZ,JY,JX));   XCAHXS=0._r8
  allocate(XCASXS(JZ,JY,JX));   XCASXS=0._r8
  allocate(XMGOXS(JZ,JY,JX));   XMGOXS=0._r8
  allocate(XMGCXS(JZ,JY,JX));   XMGCXS=0._r8
  allocate(XMGHXS(JZ,JY,JX));   XMGHXS=0._r8
  allocate(XMGSXS(JZ,JY,JX));   XMGSXS=0._r8
  allocate(XNACXS(JZ,JY,JX));   XNACXS=0._r8
  allocate(XNASXS(JZ,JY,JX));   XNASXS=0._r8
  allocate(XKASXS(JZ,JY,JX));   XKASXS=0._r8
  allocate(XF1PXS(JZ,JY,JX));   XF1PXS=0._r8
  allocate(XF2PXS(JZ,JY,JX));   XF2PXS=0._r8
  allocate(XC0PXS(JZ,JY,JX));   XC0PXS=0._r8
  allocate(XC1PXS(JZ,JY,JX));   XC1PXS=0._r8
  allocate(XC2PXS(JZ,JY,JX));   XC2PXS=0._r8
  allocate(XM1PXS(JZ,JY,JX));   XM1PXS=0._r8
  allocate(XH0BXB(JZ,JY,JX));   XH0BXB=0._r8
  allocate(XH3BXB(JZ,JY,JX));   XH3BXB=0._r8
  allocate(XF1BXB(JZ,JY,JX));   XF1BXB=0._r8
  allocate(XF2BXB(JZ,JY,JX));   XF2BXB=0._r8
  allocate(XC0BXB(JZ,JY,JX));   XC0BXB=0._r8
  allocate(XC1BXB(JZ,JY,JX));   XC1BXB=0._r8
  allocate(XC2BXB(JZ,JY,JX));   XC2BXB=0._r8
  allocate(XM1BXB(JZ,JY,JX));   XM1BXB=0._r8
  allocate(XNXFXB(JZ,JY,JX));   XNXFXB=0._r8
  allocate(TRN4S(0:JZ,JY,JX));  TRN4S=0._r8
  allocate(TRN3S(0:JZ,JY,JX));  TRN3S=0._r8
  allocate(TRN4B(JZ,JY,JX));    TRN4B=0._r8
  allocate(TRNO3(0:JZ,JY,JX));  TRNO3=0._r8
  allocate(TRN3B(JZ,JY,JX));    TRN3B=0._r8
  allocate(TRNOB(JZ,JY,JX));    TRNOB=0._r8
  allocate(TRH0P(JZ,JY,JX));    TRH0P=0._r8
  allocate(TRH1P(0:JZ,JY,JX));  TRH1P=0._r8
  allocate(TRH2P(0:JZ,JY,JX));  TRH2P=0._r8
  allocate(TRH3P(JZ,JY,JX));    TRH3P=0._r8
  allocate(TRH0B(JZ,JY,JX));    TRH0B=0._r8
  allocate(TRH1B(JZ,JY,JX));    TRH1B=0._r8
  allocate(TRH2B(JZ,JY,JX));    TRH2B=0._r8
  allocate(TRH3B(JZ,JY,JX));    TRH3B=0._r8
  allocate(TRAL(JZ,JY,JX));     TRAL=0._r8
  allocate(TRFE(JZ,JY,JX));     TRFE=0._r8
  allocate(TRHY(JZ,JY,JX));     TRHY=0._r8
  allocate(TRCA(JZ,JY,JX));     TRCA=0._r8
  allocate(TRMG(JZ,JY,JX));     TRMG=0._r8
  allocate(TRNA(JZ,JY,JX));     TRNA=0._r8
  allocate(TRKA(JZ,JY,JX));     TRKA=0._r8
  allocate(TROH(JZ,JY,JX));     TROH=0._r8
  allocate(TRSO4(JZ,JY,JX));    TRSO4=0._r8
  allocate(TRCO3(JZ,JY,JX));    TRCO3=0._r8
  allocate(TRHCO(JZ,JY,JX));    TRHCO=0._r8
  allocate(TRCO2(JZ,JY,JX));    TRCO2=0._r8
  allocate(TRH2O(0:JZ,JY,JX));  TRH2O=0._r8
  allocate(TRAL1(JZ,JY,JX));    TRAL1=0._r8
  allocate(TRAL2(JZ,JY,JX));    TRAL2=0._r8
  allocate(TRAL3(JZ,JY,JX));    TRAL3=0._r8
  allocate(TRAL4(JZ,JY,JX));    TRAL4=0._r8
  allocate(TRALS(JZ,JY,JX));    TRALS=0._r8
  allocate(TRFE1(JZ,JY,JX));    TRFE1=0._r8
  allocate(TRFE2(JZ,JY,JX));    TRFE2=0._r8
  allocate(TRFE3(JZ,JY,JX));    TRFE3=0._r8
  allocate(TRFE4(JZ,JY,JX));    TRFE4=0._r8
  allocate(TRFES(JZ,JY,JX));    TRFES=0._r8
  allocate(TRCAO(JZ,JY,JX));    TRCAO=0._r8
  allocate(TRCAC(JZ,JY,JX));    TRCAC=0._r8
  allocate(TRCAH(JZ,JY,JX));    TRCAH=0._r8
  allocate(TRCAS(JZ,JY,JX));    TRCAS=0._r8
  allocate(TRMGO(JZ,JY,JX));    TRMGO=0._r8
  allocate(TRMGC(JZ,JY,JX));    TRMGC=0._r8
  allocate(TRMGH(JZ,JY,JX));    TRMGH=0._r8
  allocate(TRMGS(JZ,JY,JX));    TRMGS=0._r8
  allocate(TRNAC(JZ,JY,JX));    TRNAC=0._r8
  allocate(TRNAS(JZ,JY,JX));    TRNAS=0._r8
  allocate(TRF1P(JZ,JY,JX));    TRF1P=0._r8
  allocate(TRF2P(JZ,JY,JX));    TRF2P=0._r8
  allocate(TRC0P(JZ,JY,JX));    TRC0P=0._r8
  allocate(TRC1P(JZ,JY,JX));    TRC1P=0._r8
  allocate(TRC2P(JZ,JY,JX));    TRC2P=0._r8
  allocate(TRM1P(JZ,JY,JX));    TRM1P=0._r8
  allocate(TRF1B(JZ,JY,JX));    TRF1B=0._r8
  allocate(TRF2B(JZ,JY,JX));    TRF2B=0._r8
  allocate(TRC0B(JZ,JY,JX));    TRC0B=0._r8
  allocate(TRC1B(JZ,JY,JX));    TRC1B=0._r8
  allocate(TRC2B(JZ,JY,JX));    TRC2B=0._r8
  allocate(TRM1B(JZ,JY,JX));    TRM1B=0._r8
  allocate(TRXN4(0:JZ,JY,JX));  TRXN4=0._r8
  allocate(TRXNB(JZ,JY,JX));    TRXNB=0._r8
  allocate(TRXHY(JZ,JY,JX));    TRXHY=0._r8
  allocate(TRXAL(JZ,JY,JX));    TRXAL=0._r8
  allocate(TRXCA(JZ,JY,JX));    TRXCA=0._r8
  allocate(TRXMG(JZ,JY,JX));    TRXMG=0._r8
  allocate(TRXNA(JZ,JY,JX));    TRXNA=0._r8
  allocate(TRXKA(JZ,JY,JX));    TRXKA=0._r8
  allocate(TRXHC(JZ,JY,JX));    TRXHC=0._r8
  allocate(TRXAL2(JZ,JY,JX));   TRXAL2=0._r8
  allocate(TRKAS(JZ,JY,JX));    TRKAS=0._r8
  allocate(TRXFE(JZ,JY,JX));    TRXFE=0._r8
  allocate(TRXFE2(JZ,JY,JX));   TRXFE2=0._r8
  allocate(TRXH0(0:JZ,JY,JX));  TRXH0=0._r8
  allocate(TRXH1(0:JZ,JY,JX));  TRXH1=0._r8
  allocate(TRXH2(0:JZ,JY,JX));  TRXH2=0._r8
  allocate(TRX1P(0:JZ,JY,JX));  TRX1P=0._r8
  allocate(TRX2P(0:JZ,JY,JX));  TRX2P=0._r8
  allocate(TRBH0(JZ,JY,JX));    TRBH0=0._r8
  allocate(TRBH1(JZ,JY,JX));    TRBH1=0._r8
  allocate(TRBH2(JZ,JY,JX));    TRBH2=0._r8
  allocate(TRB1P(JZ,JY,JX));    TRB1P=0._r8
  allocate(TRB2P(JZ,JY,JX));    TRB2P=0._r8
  allocate(TRCPMB(JZ,JY,JX));   TRCPMB=0._r8
  allocate(TBCO2(JZ,JY,JX));    TBCO2=0._r8
  allocate(TBION(0:JZ,JY,JX));  TBION=0._r8
  allocate(TRNO2(0:JZ,JY,JX));  TRNO2=0._r8
  allocate(TRN2B(JZ,JY,JX));    TRN2B=0._r8
  allocate(TRN3G(0:JZ,JY,JX));  TRN3G=0._r8
  allocate(TRALOH(JZ,JY,JX));   TRALOH=0._r8
  allocate(TRFEOH(JZ,JY,JX));   TRFEOH=0._r8
  allocate(TRCACO(JZ,JY,JX));   TRCACO=0._r8
  allocate(TRCASO(JZ,JY,JX));   TRCASO=0._r8
  allocate(TRALPO(0:JZ,JY,JX)); TRALPO=0._r8
  allocate(TRFEPO(0:JZ,JY,JX)); TRFEPO=0._r8
  allocate(TRCAPD(0:JZ,JY,JX)); TRCAPD=0._r8
  allocate(TRCAPH(0:JZ,JY,JX)); TRCAPH=0._r8
  allocate(TRCAPM(0:JZ,JY,JX)); TRCAPM=0._r8
  allocate(TRALPB(JZ,JY,JX));   TRALPB=0._r8
  allocate(TRFEPB(JZ,JY,JX));   TRFEPB=0._r8
  allocate(TRCPDB(JZ,JY,JX));   TRCPDB=0._r8
  allocate(TRCPHB(JZ,JY,JX));   TRCPHB=0._r8
  allocate(XCOBLS(JS,JY,JX));   XCOBLS=0._r8
  allocate(XCHBLS(JS,JY,JX));   XCHBLS=0._r8
  allocate(XOXBLS(JS,JY,JX));   XOXBLS=0._r8
  allocate(XNGBLS(JS,JY,JX));   XNGBLS=0._r8
  allocate(XN2BLS(JS,JY,JX));   XN2BLS=0._r8
  allocate(XN4BLW(JS,JY,JX));   XN4BLW=0._r8
  allocate(XN3BLW(JS,JY,JX));   XN3BLW=0._r8
  allocate(XNOBLW(JS,JY,JX));   XNOBLW=0._r8
  allocate(XH1PBS(JS,JY,JX));   XH1PBS=0._r8
  allocate(XH2PBS(JS,JY,JX));   XH2PBS=0._r8
  allocate(XALBLS(JS,JY,JX));   XALBLS=0._r8
  allocate(XFEBLS(JS,JY,JX));   XFEBLS=0._r8
  allocate(XHYBLS(JS,JY,JX));   XHYBLS=0._r8
  allocate(XCABLS(JS,JY,JX));   XCABLS=0._r8
  allocate(XMGBLS(JS,JY,JX));   XMGBLS=0._r8
  allocate(XNABLS(JS,JY,JX));   XNABLS=0._r8
  allocate(XKABLS(JS,JY,JX));   XKABLS=0._r8
  allocate(XOHBLS(JS,JY,JX));   XOHBLS=0._r8
  allocate(XSOBLS(JS,JY,JX));   XSOBLS=0._r8
  allocate(XCLBLS(JS,JY,JX));   XCLBLS=0._r8
  allocate(XC3BLS(JS,JY,JX));   XC3BLS=0._r8
  allocate(XHCBLS(JS,JY,JX));   XHCBLS=0._r8
  allocate(XAL1BS(JS,JY,JX));   XAL1BS=0._r8
  allocate(XAL2BS(JS,JY,JX));   XAL2BS=0._r8
  allocate(XAL3BS(JS,JY,JX));   XAL3BS=0._r8
  allocate(XAL4BS(JS,JY,JX));   XAL4BS=0._r8
  allocate(XALSBS(JS,JY,JX));   XALSBS=0._r8
  allocate(XFE1BS(JS,JY,JX));   XFE1BS=0._r8
  allocate(XFE2BS(JS,JY,JX));   XFE2BS=0._r8
  allocate(XFE3BS(JS,JY,JX));   XFE3BS=0._r8
  allocate(XFE4BS(JS,JY,JX));   XFE4BS=0._r8
  allocate(XFESBS(JS,JY,JX));   XFESBS=0._r8
  allocate(XCAOBS(JS,JY,JX));   XCAOBS=0._r8
  allocate(XCACBS(JS,JY,JX));   XCACBS=0._r8
  allocate(XCAHBS(JS,JY,JX));   XCAHBS=0._r8
  allocate(XCASBS(JS,JY,JX));   XCASBS=0._r8
  allocate(XMGOBS(JS,JY,JX));   XMGOBS=0._r8
  allocate(XMGCBS(JS,JY,JX));   XMGCBS=0._r8
  allocate(XHGBLS(JS,JY,JX));   XHGBLS=0._r8
  allocate(XMGHBS(JS,JY,JX));   XMGHBS=0._r8
  allocate(XMGSBS(JS,JY,JX));   XMGSBS=0._r8
  allocate(XNACBS(JS,JY,JX));   XNACBS=0._r8
  allocate(XNASBS(JS,JY,JX));   XNASBS=0._r8
  allocate(XKASBS(JS,JY,JX));   XKASBS=0._r8
  allocate(XH0PBS(JS,JY,JX));   XH0PBS=0._r8
  allocate(XH3PBS(JS,JY,JX));   XH3PBS=0._r8
  allocate(XF1PBS(JS,JY,JX));   XF1PBS=0._r8
  allocate(XF2PBS(JS,JY,JX));   XF2PBS=0._r8
  allocate(XC0PBS(JS,JY,JX));   XC0PBS=0._r8
  allocate(XC1PBS(JS,JY,JX));   XC1PBS=0._r8
  allocate(XC2PBS(JS,JY,JX));   XC2PBS=0._r8
  allocate(XM1PBS(JS,JY,JX));   XM1PBS=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructAquaChem
  use abortutils, only : destroy
  implicit none

  call destroy(CAL)
  call destroy(CFE)
  call destroy(CCA)
  call destroy(CMG)
  call destroy(CNA)
  call destroy(CKA)
  call destroy(CSO4)
  call destroy(CCL)
  call destroy(CALOH)
  call destroy(CFEOH)
  call destroy(CCACO)
  call destroy(CCASO)
  call destroy(CALPO)
  call destroy(CFEPO)
  call destroy(CCAPD)
  call destroy(CCAPH)
  call destroy(GKC4)
  call destroy(GKCH)
  call destroy(GKCA)
  call destroy(GKCM)
  call destroy(GKCN)
  call destroy(GKCK)
  call destroy(ZCO3)
  call destroy(ZHCO3)
  call destroy(XN4)
  call destroy(XNB)
  call destroy(ZAL)
  call destroy(ZFE)
  call destroy(ZHY)
  call destroy(ZCA)
  call destroy(ZMG)
  call destroy(ZNA)
  call destroy(ZKA)
  call destroy(ZOH)
  call destroy(ZSO4)
  call destroy(ZCL)
  call destroy(ZALOH1)
  call destroy(ZALOH2)
  call destroy(ZALOH3)
  call destroy(ZALOH4)
  call destroy(ZALS)
  call destroy(ZFEOH1)
  call destroy(ZFEOH2)
  call destroy(ZFEOH3)
  call destroy(ZFEOH4)
  call destroy(ZFES)
  call destroy(ZCAO)
  call destroy(ZCAC)
  call destroy(ZCAH)
  call destroy(ZCAS)
  call destroy(ZMGO)
  call destroy(ZMGC)
  call destroy(ZMGH)
  call destroy(ZMGS)
  call destroy(ZNAC)
  call destroy(ZNAS)
  call destroy(ZKAS)
  call destroy(H0PO4)
  call destroy(H3PO4)
  call destroy(ZFE1P)
  call destroy(ZFE2P)
  call destroy(ZCA0P)
  call destroy(ZCA1P)
  call destroy(ZCA2P)
  call destroy(ZMG1P)
  call destroy(H0POB)
  call destroy(H3POB)
  call destroy(ZFE1PB)
  call destroy(ZFE2PB)
  call destroy(ZCA0PB)
  call destroy(ZCA1PB)
  call destroy(ZCA2PB)
  call destroy(ZMG1PB)
  call destroy(XHY)
  call destroy(XAL)
  call destroy(XCA)
  call destroy(XMG)
  call destroy(XNA)
  call destroy(XKA)
  call destroy(XHC)
  call destroy(XALO2)
  call destroy(XOH0)
  call destroy(XOH1)
  call destroy(XOH2)
  call destroy(XH1P)
  call destroy(XH2P)
  call destroy(XOH0B)
  call destroy(XFE)
  call destroy(XFEO2)
  call destroy(XOH1B)
  call destroy(XOH2B)
  call destroy(XH1PB)
  call destroy(XH2PB)
  call destroy(PCAPD)
  call destroy(PCAPH)
  call destroy(PALOH)
  call destroy(PFEOH)
  call destroy(PCACO)
  call destroy(PCASO)
  call destroy(PALPO)
  call destroy(PFEPO)
  call destroy(PCAPM)
  call destroy(PALPB)
  call destroy(PFEPB)
  call destroy(PCPDB)
  call destroy(PCPMB)
  call destroy(ECND)
  call destroy(CSTR)
  call destroy(CION)
  call destroy(XCEC)
  call destroy(XAEC)
  call destroy(ZALH)
  call destroy(ZFEH)
  call destroy(ZHYH)
  call destroy(ZCCH)
  call destroy(ZMAH)
  call destroy(ZNAH)
  call destroy(ZKAH)
  call destroy(ZOHH)
  call destroy(ZSO4H)
  call destroy(ZCLH)
  call destroy(ZCO3H)
  call destroy(ZHCO3H)
  call destroy(ZALO1H)
  call destroy(ZALO2H)
  call destroy(ZALO3H)
  call destroy(ZALO4H)
  call destroy(ZALSH)
  call destroy(ZFEO1H)
  call destroy(ZFEO2H)
  call destroy(ZFEO3H)
  call destroy(ZFEO4H)
  call destroy(ZFESH)
  call destroy(ZCAOH)
  call destroy(ZCACH)
  call destroy(ZCAHH)
  call destroy(ZCASH)
  call destroy(PCPHB)
  call destroy(ZMGOH)
  call destroy(ZMGCH)
  call destroy(ZMGHH)
  call destroy(ZMGSH)
  call destroy(ZNACH)
  call destroy(ZNASH)
  call destroy(ZKASH)
  call destroy(H0PO4H)
  call destroy(H3PO4H)
  call destroy(ZFE1PH)
  call destroy(ZFE2PH)
  call destroy(ZCA0PH)
  call destroy(ZCA1PH)
  call destroy(ZCA2PH)
  call destroy(ZMG1PH)
  call destroy(H0POBH)
  call destroy(H3POBH)
  call destroy(ZFE1BH)
  call destroy(ZFE2BH)
  call destroy(ZCA0BH)
  call destroy(ZCA1BH)
  call destroy(ZCA2BH)
  call destroy(ZMG1BH)
  call destroy(XALFLS)
  call destroy(XFEFLS)
  call destroy(XHYFLS)
  call destroy(XCAFLS)
  call destroy(XMGFLS)
  call destroy(XNAFLS)
  call destroy(XKAFLS)
  call destroy(XOHFLS)
  call destroy(XSOFLS)
  call destroy(XCLFLS)
  call destroy(XC3FLS)
  call destroy(XHCFLS)
  call destroy(XAL1FS)
  call destroy(XAL2FS)
  call destroy(XAL3FS)
  call destroy(XAL4FS)
  call destroy(XALSFS)
  call destroy(XFE1FS)
  call destroy(XFE2FS)
  call destroy(XFE3FS)
  call destroy(XFE4FS)
  call destroy(XFESFS)
  call destroy(XCAOFS)
  call destroy(XCACFS)
  call destroy(XCAHFS)
  call destroy(XCASFS)
  call destroy(XMGOFS)
  call destroy(XMGCFS)
  call destroy(XMGHFS)
  call destroy(XMGSFS)
  call destroy(XNACFS)
  call destroy(XNASFS)
  call destroy(XKASFS)
  call destroy(XH0PFS)
  call destroy(XH3PFS)
  call destroy(XF1PFS)
  call destroy(XF2PFS)
  call destroy(XC0PFS)
  call destroy(XC1PFS)
  call destroy(XC2PFS)
  call destroy(XM1PFS)
  call destroy(XH0BFB)
  call destroy(XH3BFB)
  call destroy(XF1BFB)
  call destroy(XF2BFB)
  call destroy(XC0BFB)
  call destroy(XC1BFB)
  call destroy(XC2BFB)
  call destroy(XM1BFB)
  call destroy(XALFHS)
  call destroy(XFEFHS)
  call destroy(XHYFHS)
  call destroy(XCAFHS)
  call destroy(XMGFHS)
  call destroy(XNAFHS)
  call destroy(XKAFHS)
  call destroy(XOHFHS)
  call destroy(XSOFHS)
  call destroy(XCLFHS)
  call destroy(XC3FHS)
  call destroy(XHCFHS)
  call destroy(XAL1HS)
  call destroy(XAL2HS)
  call destroy(XAL3HS)
  call destroy(XAL4HS)
  call destroy(XALSHS)
  call destroy(XFE1HS)
  call destroy(XFE2HS)
  call destroy(XFE3HS)
  call destroy(XFE4HS)
  call destroy(XFESHS)
  call destroy(XCAOHS)
  call destroy(XCACHS)
  call destroy(XCAHHS)
  call destroy(XCASHS)
  call destroy(XMGOHS)
  call destroy(XMGCHS)
  call destroy(XMGHHS)
  call destroy(XMGSHS)
  call destroy(XNACHS)
  call destroy(XNASHS)
  call destroy(XKASHS)
  call destroy(XH0PHS)
  call destroy(XH3PHS)
  call destroy(XF1PHS)
  call destroy(XF2PHS)
  call destroy(XC0PHS)
  call destroy(XC1PHS)
  call destroy(XC2PHS)
  call destroy(XM1PHS)
  call destroy(XH0BHB)
  call destroy(XH3BHB)
  call destroy(XF1BHB)
  call destroy(XF2BHB)
  call destroy(XC0BHB)
  call destroy(XC1BHB)
  call destroy(XC2BHB)
  call destroy(XM1BHB)
  call destroy(XOCFXS)
  call destroy(XONFXS)
  call destroy(XOPFXS)
  call destroy(XOAFXS)
  call destroy(XCOFXS)
  call destroy(XCHFXS)
  call destroy(XOXFXS)
  call destroy(XHGFXS)
  call destroy(XNGFXS)
  call destroy(XN2FXS)
  call destroy(XN4FXW)
  call destroy(XN3FXW)
  call destroy(XNOFXW)
  call destroy(XH2PXS)
  call destroy(XN4FXB)
  call destroy(XN3FXB)
  call destroy(XNOFXB)
  call destroy(XH2BXB)
  call destroy(XNXFXS)
  call destroy(XH1PXS)
  call destroy(XH1BXB)
  call destroy(XH0PXS)
  call destroy(XH3PXS)
  call destroy(XALFXS)
  call destroy(XFEFXS)
  call destroy(XHYFXS)
  call destroy(XCAFXS)
  call destroy(XMGFXS)
  call destroy(XNAFXS)
  call destroy(XKAFXS)
  call destroy(XOHFXS)
  call destroy(XSOFXS)
  call destroy(XCLFXS)
  call destroy(XC3FXS)
  call destroy(XHCFXS)
  call destroy(XAL1XS)
  call destroy(XAL2XS)
  call destroy(XAL3XS)
  call destroy(XAL4XS)
  call destroy(XALSXS)
  call destroy(XFE1XS)
  call destroy(XFE2XS)
  call destroy(XFE3XS)
  call destroy(XFE4XS)
  call destroy(XFESXS)
  call destroy(XCAOXS)
  call destroy(XCACXS)
  call destroy(XCAHXS)
  call destroy(XCASXS)
  call destroy(XMGOXS)
  call destroy(XMGCXS)
  call destroy(XMGHXS)
  call destroy(XMGSXS)
  call destroy(XNACXS)
  call destroy(XNASXS)
  call destroy(XKASXS)
  call destroy(XF1PXS)
  call destroy(XF2PXS)
  call destroy(XC0PXS)
  call destroy(XC1PXS)
  call destroy(XC2PXS)
  call destroy(XM1PXS)
  call destroy(XH0BXB)
  call destroy(XH3BXB)
  call destroy(XF1BXB)
  call destroy(XF2BXB)
  call destroy(XC0BXB)
  call destroy(XC1BXB)
  call destroy(XC2BXB)
  call destroy(XM1BXB)
  call destroy(XNXFXB)
  call destroy(TRN4S)
  call destroy(TRN3S)
  call destroy(TRN4B)
  call destroy(TRNO3)
  call destroy(TRN3B)
  call destroy(TRNOB)
  call destroy(TRH0P)
  call destroy(TRH1P)
  call destroy(TRH2P)
  call destroy(TRH3P)
  call destroy(TRH0B)
  call destroy(TRH1B)
  call destroy(TRH2B)
  call destroy(TRH3B)
  call destroy(TRAL)
  call destroy(TRFE)
  call destroy(TRHY)
  call destroy(TRCA)
  call destroy(TRMG)
  call destroy(TRNA)
  call destroy(TRKA)
  call destroy(TROH)
  call destroy(TRSO4)
  call destroy(TRCO3)
  call destroy(TRHCO)
  call destroy(TRCO2)
  call destroy(TRH2O)
  call destroy(TRAL1)
  call destroy(TRAL2)
  call destroy(TRAL3)
  call destroy(TRAL4)
  call destroy(TRALS)
  call destroy(TRFE1)
  call destroy(TRFE2)
  call destroy(TRFE3)
  call destroy(TRFE4)
  call destroy(TRFES)
  call destroy(TRCAO)
  call destroy(TRCAC)
  call destroy(TRCAH)
  call destroy(TRCAS)
  call destroy(TRMGO)
  call destroy(TRMGC)
  call destroy(TRMGH)
  call destroy(TRMGS)
  call destroy(TRNAC)
  call destroy(TRNAS)
  call destroy(TRF1P)
  call destroy(TRF2P)
  call destroy(TRC0P)
  call destroy(TRC1P)
  call destroy(TRC2P)
  call destroy(TRM1P)
  call destroy(TRF1B)
  call destroy(TRF2B)
  call destroy(TRC0B)
  call destroy(TRC1B)
  call destroy(TRC2B)
  call destroy(TRM1B)
  call destroy(TRXN4)
  call destroy(TRXNB)
  call destroy(TRXHY)
  call destroy(TRXAL)
  call destroy(TRXCA)
  call destroy(TRXMG)
  call destroy(TRXNA)
  call destroy(TRXKA)
  call destroy(TRXHC)
  call destroy(TRXAL2)
  call destroy(TRKAS)
  call destroy(TRXFE)
  call destroy(TRXFE2)
  call destroy(TRXH0)
  call destroy(TRXH1)
  call destroy(TRXH2)
  call destroy(TRX1P)
  call destroy(TRX2P)
  call destroy(TRBH0)
  call destroy(TRBH1)
  call destroy(TRBH2)
  call destroy(TRB1P)
  call destroy(TRB2P)
  call destroy(TRCPMB)
  call destroy(TBCO2)
  call destroy(TBION)
  call destroy(TRNO2)
  call destroy(TRN2B)
  call destroy(TRN3G)
  call destroy(TRALOH)
  call destroy(TRFEOH)
  call destroy(TRCACO)
  call destroy(TRCASO)
  call destroy(TRALPO)
  call destroy(TRFEPO)
  call destroy(TRCAPD)
  call destroy(TRCAPH)
  call destroy(TRCAPM)
  call destroy(TRALPB)
  call destroy(TRFEPB)
  call destroy(TRCPDB)
  call destroy(TRCPHB)
  call destroy(XCOBLS)
  call destroy(XCHBLS)
  call destroy(XOXBLS)
  call destroy(XNGBLS)
  call destroy(XN2BLS)
  call destroy(XN4BLW)
  call destroy(XN3BLW)
  call destroy(XNOBLW)
  call destroy(XH1PBS)
  call destroy(XH2PBS)
  call destroy(XALBLS)
  call destroy(XFEBLS)
  call destroy(XHYBLS)
  call destroy(XCABLS)
  call destroy(XMGBLS)
  call destroy(XNABLS)
  call destroy(XKABLS)
  call destroy(XOHBLS)
  call destroy(XSOBLS)
  call destroy(XCLBLS)
  call destroy(XC3BLS)
  call destroy(XHCBLS)
  call destroy(XAL1BS)
  call destroy(XAL2BS)
  call destroy(XAL3BS)
  call destroy(XAL4BS)
  call destroy(XALSBS)
  call destroy(XFE1BS)
  call destroy(XFE2BS)
  call destroy(XFE3BS)
  call destroy(XFE4BS)
  call destroy(XFESBS)
  call destroy(XCAOBS)
  call destroy(XCACBS)
  call destroy(XCAHBS)
  call destroy(XCASBS)
  call destroy(XMGOBS)
  call destroy(XMGCBS)
  call destroy(XHGBLS)
  call destroy(XMGHBS)
  call destroy(XMGSBS)
  call destroy(XNACBS)
  call destroy(XNASBS)
  call destroy(XKASBS)
  call destroy(XH0PBS)
  call destroy(XH3PBS)
  call destroy(XF1PBS)
  call destroy(XF2PBS)
  call destroy(XC0PBS)
  call destroy(XC1PBS)
  call destroy(XC2PBS)
  call destroy(XM1PBS)

  end subroutine DestructAquaChem

end module AqueChemDatatype
