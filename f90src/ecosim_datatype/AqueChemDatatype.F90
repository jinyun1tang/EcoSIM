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
  real(r8),target,allocatable :: trcSalt_solml_vr(:,:,:,:)           !soil aqueous salt content micropre, [mol d-2]
  real(r8),target,allocatable :: trcx_solml_vr(:,:,:,:)              !exchangeable tracers
  real(r8),target,allocatable :: trcp_saltpml_vr(:,:,:,:)            !salt precipitate in micropore

  real(r8),target,allocatable ::  ECND_vr(:,:,:)                     !electrical conductivity , [dS m-1]
  real(r8),target,allocatable ::  CSTR(:,:,:)                        !solution ion strength, [mol m-3]
  real(r8),target,allocatable ::  CION(:,:,:)                        !solution ion concentratiom, [mol m-3]

  real(r8),target,allocatable :: trcSalt_soHml_vr(:,:,:,:)
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
  real(r8),target,allocatable ::  TR_HCO3_col(:,:,:)                       !total solute HCO3 transformation, [mol d-2 h-1]
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
  real(r8),target,allocatable ::  TR_Fe_sorbed_soil_vr(:,:,:)                       !total Fe adsorption
  real(r8),target,allocatable ::  TR_FeO2H2_sorbed_soil_vr(:,:,:)                      !total FeOH2 adsorption
  real(r8),target,allocatable ::  trcx_TRSoilChem_vr(:,:,:,:)                   !total adsorbed OH- transformation non-band, [mol d-2 h-1]
  real(r8),target,allocatable ::  Txchem_CO2_vr(:,:,:)                       !total solute CO2 transformation boundary, [mol d-2 h-1]
  real(r8),target,allocatable ::  TBION(:,:,:)                       !total solute ion transformation boundary, [mol d-2 h-1]
  real(r8),target,allocatable ::  TRN3G(:,:,:)                       !total gaseous NH3 transformation, [mol d-2 h-1]
  real(r8),target,allocatable ::  trcp_RChem_soil(:,:,:,:)                   !total precipitated P containing transformation non-band, [mol d-2 h-1]
  real(r8),target,allocatable ::  trcg_Xbndl_flx(:,:,:,:)
  real(r8),target,allocatable ::  trcn_Xbndl_flx(:,:,:,:)
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
  allocate(trcx_solml_vr(idx_beg:idx_end,0:JZ,JY,JX));trcx_solml_vr=0._r8
  allocate(trcp_saltpml_vr(idsp_beg:idsp_end,0:JZ,JY,JX)); trcp_saltpml_vr=0._r8
  allocate(ECND_vr(JZ,JY,JX));     ECND_vr=0._r8
  allocate(CSTR(JZ,JY,JX));     CSTR=0._r8
  allocate(CION(JZ,JY,JX));     CION=0._r8

  allocate(DOM_PoreTranspFlx(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_PoreTranspFlx=0._r8
  allocate(trcs_PoreTranspFlx_vr(ids_beg:ids_end,JZ,JY,JX));   trcs_PoreTranspFlx_vr=0._r8

  allocate(trcn_RChem_soil_vr(ids_nut_beg:ids_nuts_end,0:JZ,JY,JX)); trcn_RChem_soil_vr=0._r8

  allocate(TR_NH3_soil_vr(0:JZ,JY,JX));  TR_NH3_soil_vr=0._r8

  allocate(trcn_RChem_band_soil_vr(ids_nutb_beg:ids_nutb_end,JZ,JY,JX)); trcn_RChem_band_soil_vr=0._r8

  allocate(trcx_TRSoilChem_vr(idx_beg:idx_end,0:JZ,JY,JX));  trcx_TRSoilChem_vr=0._r8

  allocate(TR_HCO3_col(JZ,JY,JX));    TR_HCO3_col=0._r8
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
  allocate(TR_Fe_sorbed_soil_vr(JZ,JY,JX));    TR_Fe_sorbed_soil_vr=0._r8
  allocate(TR_FeO2H2_sorbed_soil_vr(JZ,JY,JX));   TR_FeO2H2_sorbed_soil_vr=0._r8
  allocate(Txchem_CO2_vr(JZ,JY,JX));    Txchem_CO2_vr=0._r8
  allocate(TBION(0:JZ,JY,JX));  TBION=0._r8
  allocate(TRN3G(0:JZ,JY,JX));  TRN3G=0._r8
  allocate(trcp_RChem_soil(idsp_beg:idsp_end,0:JZ,JY,JX)); trcp_RChem_soil=0._r8
  allocate(trcg_Xbndl_flx(idg_beg:idg_end-1,JS,JY,JX)); trcg_Xbndl_flx=0._r8
  allocate(trcn_Xbndl_flx(ids_nut_beg:ids_nuts_end,JS,JY,JX)); trcn_Xbndl_flx=0._r8
  if(salt_model)then
    allocate(trcSaltFlo2SnowLay(idsalt_beg:idsalt_end,JS,JY,JX)); trcSaltFlo2SnowLay=0._r8
    allocate(trcSalt3DFlo2Cell(idsalt_beg:idsaltb_end,3,0:JD,JV,JH));trcSalt3DFlo2Cell=0._r8
    allocate(trcSalt_XFHS(idsalt_beg:idsaltb_end,3,JD,JV,JH));trcSalt_XFHS=0._r8
    allocate(trcSalt_solml_vr(idsalt_beg:idsaltb_end,0:JZ,JY,JX));trcSalt_solml_vr=0._r8
    allocate(trcsalt_rain_conc(idsalt_beg:idsalt_end,JY,JX));trcsalt_rain_conc=0._r8
    allocate(trcSaltIonNumber(idsalt_beg:idsaltb_end))
    allocate(trcSalt_soHml_vr(idsalt_beg:idsaltb_end,JZ,JY,JX)); trcSalt_soHml_vr=0._r8
    allocate(trcSalt_XFXS(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_XFXS=0._r8
    allocate(trcSalt_TR(idsalt_beg:idsaltb_end,JZ,JY,JX));    trcSalt_TR=0._r8
  endif

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
  call destroy(trcx_solml_vr)
  call destroy(trcSalt_solml_vr)
  call destroy(trcsalt_rain_conc)
  call destroy(trcSalt_soHml_vr)
  call destroy(trcp_saltpml_vr)
  call destroy(ECND_vr)
  call destroy(CSTR)
  call destroy(CION)
  call destroy(trcSalt_XFXS)
  call destroy(trcSalt3DFlo2Cell)
  call destroy(trcSalt_XFHS)
  call destroy(TR_NH3_soil_vr)
  call destroy(trcn_RChem_band_soil_vr)
  call destroy(TR_HCO3_col)
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
  call destroy(TR_Fe_sorbed_soil_vr)
  call destroy(TR_FeO2H2_sorbed_soil_vr)
  call destroy(Txchem_CO2_vr)
  call destroy(TBION)
  call destroy(TRN3G)


  call destroy(trcSaltIonNumber)
  call destroy(trcs_PoreTranspFlx_vr)
  call destroy(DOM_PoreTranspFlx)
  call destroy(trcn_RChem_soil_vr)
  end subroutine DestructAquaChem

end module AqueChemDatatype
