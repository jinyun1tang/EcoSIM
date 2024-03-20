module AqueChemDatatype
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

  real(r8),target,allocatable ::  CAL(:,:,:)                         !soil Al content, [mg Al kg-1]
  real(r8),target,allocatable ::  CFE(:,:,:)                         !soil Fe content, [mg Fe kg-3]
  real(r8),target,allocatable ::  CCA(:,:,:)                         !soil Ca content, [mg Ca kg-3]
  real(r8),target,allocatable ::  CMG(:,:,:)                         !soil Mg content, [mg Mg kg-3]
  real(r8),target,allocatable ::  CNA(:,:,:)                         !soil Na content, [mg Na kg-3]
  real(r8),target,allocatable ::  CKA(:,:,:)                         !soil K content, [mg K kg-3]
  real(r8),target,allocatable ::  CSO4(:,:,:)                        !soil SO4 content, [mg S kg-3]
  real(r8),target,allocatable ::  CCL(:,:,:)                         !soil Cl content, [mg Cl kg-1]
  real(r8),target,allocatable ::  CALOH(:,:,:)                       !soil AlOH3 content, [mg Al kg-1]
  real(r8),target,allocatable ::  CFEOH(:,:,:)                       !soil FeOH3 content, [mg Fe kg-1]
  real(r8),target,allocatable ::  CCACO(:,:,:)                       !soil CaCO3 content, [mg Ca kg-1]
  real(r8),target,allocatable ::  CCASO(:,:,:)                       !soil CaSO4 content, [mg Ca kg-1]
  real(r8),target,allocatable ::  CALPO(:,:,:)                       !soil AlPO4 content, [mg P kg-1]
  real(r8),target,allocatable ::  CFEPO(:,:,:)                       !soil FePO4 content, [mg P kg-1]
  real(r8),target,allocatable ::  CCAPD(:,:,:)                       !soil CaHPO4 content, [mg P kg-1]
  real(r8),target,allocatable ::  CCAPH(:,:,:)                       !soil apatite content, [mg P kg-1]
  real(r8),target,allocatable ::  GKC4(:,:,:)                        !Ca-NH4 Gapon selectivity coefficient, [-]
  real(r8),target,allocatable ::  GKCH(:,:,:)                        !Ca-H Gapon selectivity coefficient, [-]
  real(r8),target,allocatable ::  GKCA(:,:,:)                        !Ca-Al Gapon selectivity coefficient, [-]
  real(r8),target,allocatable ::  GKCM(:,:,:)                        !Ca-Mg Gapon selectivity coefficient, [-]
  real(r8),target,allocatable ::  GKCN(:,:,:)                        !Ca-Na Gapon selectivity coefficient, [-]
  real(r8),target,allocatable ::  GKCK(:,:,:)                        !Ca-K Gapon selectivity coefficient, [-]

  real(r8),target,allocatable :: trcsalt_rain_conc(:,:,:)            !salt tracer concentration in rain [g m-3]
  real(r8),target,allocatable :: trcSalt_solml(:,:,:,:)                !soil aqueous salt content micropre, [mol d-2]
  real(r8),target,allocatable :: trcx_solml(:,:,:,:)                 !exchangeable tracers
  real(r8),target,allocatable :: trcp_salml(:,:,:,:)                !salt precipitate in micropore

  real(r8),target,allocatable ::  ECND(:,:,:)                        !electrical conductivity , [dS m-1]
  real(r8),target,allocatable ::  CSTR(:,:,:)                        !solution ion strength, [mol m-3]
  real(r8),target,allocatable ::  CION(:,:,:)                        !solution ion concentratiom, [mol m-3]

  real(r8),target,allocatable :: trcSalt_soHml(:,:,:,:)
  real(r8),target,allocatable :: trcSalt_XFHS(:,:,:,:,:)
  real(r8),target,allocatable :: trcSalt3DFlo2Cell(:,:,:,:,:)
  real(r8),target,allocatable :: trcSaltIonNumber(:)                 !number of ions when the salt is fully dissociated
  real(r8),target,allocatable ::  DOM_PoreTranspFlx(:,:,:,:,:)                    !total DOC micropore-macropore transfer, [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_PoreTranspFlx_vr(:,:,:,:)                 !total non-salt solute micropore-macropore transfer, [g d-2 h-1]
  real(r8),target,allocatable ::  trcSalt_XFXS(:,:,:,:)                !total salt micropore-macropore transfer non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  trcn_RChem_soil_vr(:,:,:,:)                       !total solute NH4 transformation non-band, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_NH3_soil_vr(:,:,:)                       !total solute NH3 transformation non-band, [mol d-2 h-1]

  real(r8),target,allocatable ::  trcn_RChem_band_soil_vr(:,:,:,:)     !total solute nutrient transformation band, [mol d-2 h-1]
  real(r8),target,allocatable ::  trcSalt_TR(:,:,:,:)                  !total salt solute transformation non-band, [mol d-2 h-1]
  real(r8),target,allocatable ::  TRHCO(:,:,:)                       !total solute HCO3 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_CO2_aqu_soil_vr(:,:,:)                       !total solute CO2 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TRH2O(:,:,:)                       !total solute H2O transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_FeO3H3_soil(:,:,:)                       !total solute FeOH3 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_H_p_sorbed_soil(:,:,:)                       !total adsorbed H transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_Al_sorbed_soil(:,:,:)                       !total adsorbed Al transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_Ca_sorbed_soil(:,:,:)                       !total adsorbed Ca transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_Mg_sorbed_soil(:,:,:)                       !total adsorbed Mg transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_Na_sorbed_soil(:,:,:)                       !total adsorbed Na transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_K_sorbed_soil(:,:,:)                       !total adsorbed K transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_HCO3_sorbed_soil(:,:,:)                       !total adsorbed COOH transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_AlO2H2_sorbed_soil(:,:,:)                      !total adsorbed AlOH2 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_KSO4_soil(:,:,:)                       !total solute KSO4 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  TR_Fe_sorbed_soil(:,:,:)                       !total Fe adsorption
  real(r8),target,allocatable ::  TR_FeO2H2_sorbed_soil(:,:,:)                      !total FeOH2 adsorption
  real(r8),target,allocatable ::  trcx_TRSoilChem_vr(:,:,:,:)                   !total adsorbed OH- transformation non-band, [mol d-2 h-1]
  real(r8),target,allocatable ::  TBCO2(:,:,:)                       !total solute CO2 transformation boundary, [mol d-2 h-1]
  real(r8),target,allocatable ::  TBION(:,:,:)                       !total solute ion transformation boundary, [mol d-2 h-1]
  real(r8),target,allocatable ::  TRN3G(:,:,:)                       !total gaseous NH3 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  trcp_RChem_soil(:,:,:,:)                   !total precipitated P containing transformation non-band, [mol d-2 h-1]
  real(r8),target,allocatable ::  trcg_XBLS(:,:,:,:)
  real(r8),target,allocatable ::  trcn_XBLS(:,:,:,:)
  real(r8),target,allocatable ::  trcSaltFlo2SnowLay(:,:,:,:)
  real(r8),target,allocatable ::  XCOBLS(:,:,:)                      !wet deposition of CO2, [g d-2 h-1]
  real(r8),target,allocatable ::  XCHBLS(:,:,:)                      !wet deposition of CH4, [g d-2 h-1]
  real(r8),target,allocatable ::  XOXBLS(:,:,:)                      !wet deposition of O2, [g d-2 h-1]
  real(r8),target,allocatable ::  XNGBLS(:,:,:)                      !wet deposition of N2, [g d-2 h-1]
  real(r8),target,allocatable ::  XN2BLS(:,:,:)                      !wet deposition of N2O, [g d-2 h-1]
  real(r8),target,allocatable ::  XN4BLW(:,:,:)                      !wet deposition of NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  XN3BLW(:,:,:)                      !wet deposition of NH3, [g d-2 h-1]
  real(r8),target,allocatable ::  XNOBLW(:,:,:)                      !wet deposition of NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  XH1PBS(:,:,:)                      !wet deposition of HPO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XH2PBS(:,:,:)                      !wet deposition of H2PO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XALBLS(:,:,:)                      !wet deposition of Al, [g d-2 h-1]
  real(r8),target,allocatable ::  XFEBLS(:,:,:)                      !wet deposition of Fe, [g d-2 h-1]
  real(r8),target,allocatable ::  XHYBLS(:,:,:)                      !wet deposition of H, [g d-2 h-1]
  real(r8),target,allocatable ::  XCABLS(:,:,:)                      !wet deposition of Ca, [g d-2 h-1]
  real(r8),target,allocatable ::  XMGBLS(:,:,:)                      !wet deposition of Mg, [g d-2 h-1]
  real(r8),target,allocatable ::  XNABLS(:,:,:)                      !wet deposition of Na, [g d-2 h-1]
  real(r8),target,allocatable ::  XKABLS(:,:,:)                      !wet deposition of K, [g d-2 h-1]
  real(r8),target,allocatable ::  XOHBLS(:,:,:)                      !wet deposition of OH, [g d-2 h-1]
  real(r8),target,allocatable ::  XSOBLS(:,:,:)                      !wet deposition of SO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XCLBLS(:,:,:)                      !wet deposition of Cl, [g d-2 h-1]
  real(r8),target,allocatable ::  XC3BLS(:,:,:)                      !wet deposition of CO3, [g d-2 h-1]
  real(r8),target,allocatable ::  XHCBLS(:,:,:)                      !wet deposition of HCO3, [g d-2 h-1]
  real(r8),target,allocatable ::  XAL1BS(:,:,:)                      !wet deposition of AlOH, [g d-2 h-1]
  real(r8),target,allocatable ::  XAL2BS(:,:,:)                      !wet deposition of AlOH2, [g d-2 h-1]
  real(r8),target,allocatable ::  XAL3BS(:,:,:)                      !wet deposition of AlOH3, [g d-2 h-1]
  real(r8),target,allocatable ::  XAL4BS(:,:,:)                      !wet deposition of AlOH4, [g d-2 h-1]
  real(r8),target,allocatable ::  XALSBS(:,:,:)                      !wet deposition of AlSO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XFE1BS(:,:,:)                      !wet deposition of FeOH, [g d-2 h-1]
  real(r8),target,allocatable ::  XFE2BS(:,:,:)                      !wet deposition of FeOH2, [g d-2 h-1]
  real(r8),target,allocatable ::  XFE3BS(:,:,:)                      !wet deposition of FeOH3, [g d-2 h-1]
  real(r8),target,allocatable ::  XFE4BS(:,:,:)                      !wet deposition of FeOH4, [g d-2 h-1]
  real(r8),target,allocatable ::  XFESBS(:,:,:)                      !wet deposition of FeSO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XCAOBS(:,:,:)                      !wet deposition of CaOH, [g d-2 h-1]
  real(r8),target,allocatable ::  XCACBS(:,:,:)                      !wet deposition of CaCO3, [g d-2 h-1]
  real(r8),target,allocatable ::  XCAHBS(:,:,:)                      !wet deposition of CaHCO3, [g d-2 h-1]
  real(r8),target,allocatable ::  XCASBS(:,:,:)                      !wet deposition of CaSO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XMGOBS(:,:,:)                      !wet deposition of MgOH, [g d-2 h-1]
  real(r8),target,allocatable ::  XMGCBS(:,:,:)                      !wet deposition of MgCO3, [g d-2 h-1]
  real(r8),target,allocatable ::  XHGBLS(:,:,:)                      !wet deposition of H2, [g d-2 h-1]
  real(r8),target,allocatable ::  XMGHBS(:,:,:)                      !wet deposition of MgHCO3, [g d-2 h-1]
  real(r8),target,allocatable ::  XMGSBS(:,:,:)                      !wet deposition of MgSO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XNACBS(:,:,:)                      !wet deposition of NaCO3, [g d-2 h-1]
  real(r8),target,allocatable ::  XNASBS(:,:,:)                      !wet deposition of NaSO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XKASBS(:,:,:)                      !wet deposition of KSO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XH0PBS(:,:,:)                      !wet deposition of PO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XH3PBS(:,:,:)                      !wet deposition of H3PO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XF1PBS(:,:,:)                      !wet deposition of FeHPO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XF2PBS(:,:,:)                      !wet deposition of FeH2PO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XC0PBS(:,:,:)                      !wet deposition of CaPO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XC1PBS(:,:,:)                      !wet deposition of CHPO4, [g d-2 h-1]
  real(r8),target,allocatable ::  XC2PBS(:,:,:)                      !wet deposition of CaH4P2O8, [g d-2 h-1]
  real(r8),target,allocatable ::  XM1PBS(:,:,:)                      !wet deposition of MgHPO4, [g d-2 h-1]
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

  allocate(trcSalt_solml(idsalt_beg:idsaltb_end,0:JZ,JY,JX));trcSalt_solml=0._r8
  allocate(trcsalt_rain_conc(idsalt_beg:idsalt_end,JY,JX));trcsalt_rain_conc=0._r8
  allocate(trcx_solml(idx_beg:idx_end,0:JZ,JY,JX));trcx_solml=0._r8
  allocate(trcp_salml(idsp_beg:idsp_end,0:JZ,JY,JX)); trcp_salml=0._r8
  allocate(trcSaltIonNumber(idsalt_beg:idsaltb_end))
  allocate(ECND(JZ,JY,JX));     ECND=0._r8
  allocate(CSTR(JZ,JY,JX));     CSTR=0._r8
  allocate(CION(JZ,JY,JX));     CION=0._r8

  allocate(trcSalt_soHml(idsalt_beg:idsaltb_end,JZ,JY,JX)); trcSalt_soHml=0._r8
  allocate(DOM_PoreTranspFlx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_PoreTranspFlx=0._r8
  allocate(trcs_PoreTranspFlx_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_PoreTranspFlx_vr=0._r8
  allocate(trcSalt_XFXS(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_XFXS=0._r8
  allocate(trcn_RChem_soil_vr(ids_nut_beg:ids_nuts_end,0:JZ,JY,JX)); trcn_RChem_soil_vr=0._r8

  allocate(TR_NH3_soil_vr(0:JZ,JY,JX));  TR_NH3_soil_vr=0._r8

  allocate(trcn_RChem_band_soil_vr(ids_nutb_beg:ids_nutb_end,JZ,JY,JX)); trcn_RChem_band_soil_vr=0._r8

  allocate(trcSalt_TR(idsalt_beg:idsaltb_end,JZ,JY,JX));    trcSalt_TR=0._r8
  allocate(trcx_TRSoilChem_vr(idx_beg:idx_end,0:JZ,JY,JX));  trcx_TRSoilChem_vr=0._r8

  allocate(TRHCO(JZ,JY,JX));    TRHCO=0._r8
  allocate(TR_CO2_aqu_soil_vr(JZ,JY,JX));    TR_CO2_aqu_soil_vr=0._r8
  allocate(TRH2O(0:JZ,JY,JX));  TRH2O=0._r8
  allocate(TR_FeO3H3_soil(JZ,JY,JX));    TR_FeO3H3_soil=0._r8
  allocate(TR_H_p_sorbed_soil(JZ,JY,JX));    TR_H_p_sorbed_soil=0._r8
  allocate(TR_Al_sorbed_soil(JZ,JY,JX));    TR_Al_sorbed_soil=0._r8
  allocate(TR_Ca_sorbed_soil(JZ,JY,JX));    TR_Ca_sorbed_soil=0._r8
  allocate(TR_Mg_sorbed_soil(JZ,JY,JX));    TR_Mg_sorbed_soil=0._r8
  allocate(TR_Na_sorbed_soil(JZ,JY,JX));    TR_Na_sorbed_soil=0._r8
  allocate(TR_K_sorbed_soil(JZ,JY,JX));    TR_K_sorbed_soil=0._r8
  allocate(TR_HCO3_sorbed_soil(JZ,JY,JX));    TR_HCO3_sorbed_soil=0._r8
  allocate(TR_AlO2H2_sorbed_soil(JZ,JY,JX));   TR_AlO2H2_sorbed_soil=0._r8
  allocate(TR_KSO4_soil(JZ,JY,JX));    TR_KSO4_soil=0._r8
  allocate(TR_Fe_sorbed_soil(JZ,JY,JX));    TR_Fe_sorbed_soil=0._r8
  allocate(TR_FeO2H2_sorbed_soil(JZ,JY,JX));   TR_FeO2H2_sorbed_soil=0._r8
  allocate(TBCO2(JZ,JY,JX));    TBCO2=0._r8
  allocate(TBION(0:JZ,JY,JX));  TBION=0._r8
  allocate(TRN3G(0:JZ,JY,JX));  TRN3G=0._r8
  allocate(trcp_RChem_soil(idsp_beg:idsp_end,0:JZ,JY,JX)); trcp_RChem_soil=0._r8
  allocate(trcg_XBLS(idg_beg:idg_end-1,JS,JY,JX)); trcg_XBLS=0._r8
  allocate(trcn_XBLS(ids_nut_beg:ids_nuts_end,JS,JY,JX)); trcn_XBLS=0._r8
  if(salt_model)then
    allocate(trcSaltFlo2SnowLay(idsalt_beg:idsalt_end,JS,JY,JX)); trcSaltFlo2SnowLay=0._r8
    allocate(trcSalt3DFlo2Cell(idsalt_beg:idsaltb_end,3,0:JD,JV,JH));trcSalt3DFlo2Cell=0._r8
    allocate(trcSalt_XFHS(idsalt_beg:idsaltb_end,3,JD,JV,JH));trcSalt_XFHS=0._r8
  endif
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

  call destroy(trcSalt_TR)
  call destroy(trcSalt_TR)
  call destroy(trcx_solml)
  call destroy(trcSalt_solml)
  call destroy(trcsalt_rain_conc)
  call destroy(trcSalt_soHml)
  call destroy(trcp_salml)
  call destroy(ECND)
  call destroy(CSTR)
  call destroy(CION)
  call destroy(trcSalt_XFXS)
  call destroy(trcSalt3DFlo2Cell)
  call destroy(trcSalt_XFHS)
  call destroy(TR_NH3_soil_vr)
  call destroy(trcn_RChem_band_soil_vr)
  call destroy(TRHCO)
  call destroy(TR_CO2_aqu_soil_vr)
  call destroy(TRH2O)
  call destroy(TR_FeO3H3_soil)
  call destroy(TR_Al_sorbed_soil)
  call destroy(TR_Ca_sorbed_soil)
  call destroy(TR_Mg_sorbed_soil)
  call destroy(TR_Na_sorbed_soil)
  call destroy(TR_K_sorbed_soil)
  call destroy(TR_HCO3_sorbed_soil)
  call destroy(TR_AlO2H2_sorbed_soil)
  call destroy(TR_KSO4_soil)
  call destroy(TR_Fe_sorbed_soil)
  call destroy(TR_FeO2H2_sorbed_soil)
  call destroy(TBCO2)
  call destroy(TBION)
  call destroy(TRN3G)
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
  call destroy(trcSaltIonNumber)
  call destroy(trcs_PoreTranspFlx_vr)
  call destroy(DOM_PoreTranspFlx)
  call destroy(trcn_RChem_soil_vr)
  end subroutine DestructAquaChem

end module AqueChemDatatype
