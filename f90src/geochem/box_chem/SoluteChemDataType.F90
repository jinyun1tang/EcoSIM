module SoluteChemDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8

  implicit none
  public
  CHARACTER(LEN=*), private, PARAMETER :: MOD_FILENAME = &
  __FILE__

  type, public :: solutedtype
    real(r8) :: H2CO3_aqu_conc
    real(r8) :: CH4_aqu_conc
    real(r8) :: O2_aqu_conc
    real(r8) :: N2_aqu_conc
    real(r8) :: N2O_aqu_conc
    real(r8) :: NH4_1p_conc
    real(r8) :: NH3_aqu_conc
    real(r8) :: Al_3p_conc
    real(r8) :: Fe_3p_conc
    real(r8) :: H_1p_conc
    real(r8) :: Ca_2p_conc
    real(r8) :: Mg_2p_conc
    real(r8) :: Na_1p_conc
    real(r8) :: K_1p_conc
    real(r8) :: OH_1e_conc
    real(r8) :: SO4_2e_conc
    real(r8) :: Cl_e_conc
    real(r8) :: CO3_2e_conc
    real(r8) :: HCO3_e_conc
    real(r8) :: AlOH_2p_conc
    real(r8) :: AlO2H2_1p_conc
    real(r8) :: AlO3H3_conc
    real(r8) :: AlO4H4_1e_conc
    real(r8) :: AlSO4_1p_conc
    real(r8) :: FeOH_2p_conc
    real(r8) :: FeO2H2_p_conc
    real(r8) :: FeO3H3_conc
    real(r8) :: FeO4H4_1e_conc
    real(r8) :: FeSO4_1p_conc
    real(r8) :: CaO2H2_conc
    real(r8) :: CaCO3_conc
    real(r8) :: CaHCO3_1p_conc
    real(r8) :: CaSO4_conc
    real(r8) :: MgOH_1p_conc
    real(r8) :: MgCO3_conc
    real(r8) :: MgHCO3_1p_conc
    real(r8) :: MgSO4_conc
    real(r8) :: NaCO3_1e_conc
    real(r8) :: NaSO4_1e_conc
    real(r8) :: KSO4_1e_conc
    real(r8) :: H0PO4_3e_conc
    real(r8) :: H1PO4_2e_conc
    real(r8) :: H2PO4_1e_conc
    real(r8) :: H3PO4_conc
    real(r8) :: FeHPO4_p_conc
    real(r8) :: FeH2PO4_2p_conc
    real(r8) :: CaPO4_1e_con
    real(r8) :: CaHPO4_conc
    real(r8) :: CaH4P2O8_1p_conc
    real(r8) :: MgHPO4_conc
    real(r8) :: CSTR1
    real(r8) :: CCO2M
    real(r8) :: CCH4M
    real(r8) :: COXYM
    real(r8) :: CZ2GM
    real(r8) :: CZ2OM
    real(r8) :: CN4Z
    real(r8) :: CNOZ
    real(r8) :: CNAZ
    real(r8) :: CKAZ
    real(r8) :: CSOZ
    real(r8) :: CCLZ
    real(r8) :: CNOX
    real(r8) :: CCASOX
    real(r8) :: CN4X
    real(r8) :: CPOZ
    real(r8) :: CALZ
    real(r8) :: CFEZ
    real(r8) :: CCAZ
    real(r8) :: CMGZ
    real(r8) :: CALX
    real(r8) :: CFEX
    real(r8) :: CaX_conc
    real(r8) :: MgX_conc
    real(r8) :: CNAX
    real(r8) :: CKAX
    real(r8) :: CSOX
    real(r8) :: CCLX
    real(r8) :: CALPOX
    real(r8) :: CFEPOX
    real(r8) :: CCAPDX
    real(r8) :: CCAPHX
    real(r8) :: CALOHX
    real(r8) :: CFEOHX
    real(r8) :: CCACOX
    real(r8) :: XNH4_conc
    real(r8) :: XHY1
    real(r8) :: XAl_conc
    real(r8) :: XFe_conc
    real(r8) :: XCa_conc
    real(r8) :: Precp_Ca5P3O12O3H3_conc
    real(r8) :: XMg_conc
    real(r8) :: XNa_conc
    real(r8) :: XK_conc
    real(r8) :: XHC1
    real(r8) :: XAlO2H2_conc
    real(r8) :: XFeO2H2_conc
    real(r8) :: XOH_conc
    real(r8) :: XROH1_conc
    real(r8) :: XROH2_conc
    real(r8) :: XHPO4_conc
    real(r8) :: XH2PO4_conc
    real(r8) :: Precp_AlO3H3_conc
    real(r8) :: Precp_FeO3H3_conc
    real(r8) :: Precp_CaCO3_conc
    real(r8) :: Precp_CaSO4_conc
    real(r8) :: Precp_AlPO4_conc
    real(r8) :: Precp_FePO4_conc
    real(r8) :: Precp_CaHPO4_conc
    real(r8) :: FH2O
    real(r8) :: ATCA
    real(r8) :: XAEC
    real(r8) :: CEC
    real(r8) :: ORGC
    real(r8) :: VLPO4
    real(r8) :: XCEC
    real(r8) :: GKC4
    real(r8) :: GKCA
    real(r8) :: GKCH
    real(r8) :: GKCK
    real(r8) :: GKCN
    real(r8) :: GKCM
    real(r8) :: ZEROS
    real(r8) :: VLWatMicP
  end type solutedtype

  type, public :: solute_flx_type
    real(r8) :: TR_NH4_soil
    real(r8) :: TR_NH4_band_soil
    real(r8) :: TR_NH3_soil_vr
    real(r8) :: TR_NH3_band_soil
    real(r8) :: TR_H1PO4_soil
    real(r8) :: TR_H2PO4_soil
    real(r8) :: TR_H1PO4_band_soil
    real(r8) :: TR_H2PO4_band_soil
    real(r8) :: TR_NH4_sorbed_soil
    real(r8) :: TR_NH4_sorbed_band_soil
    real(r8) :: TR_ROH_sorbed_soil
    real(r8) :: TR_ROH2_sorbed_soil
    real(r8) :: TR_RHPO4_sorbed_soil
    real(r8) :: TR_RH2PO4_sorbed_soil
    real(r8) :: TR_ROH_sorbed_band_soil
    real(r8) :: TR_ROH2_sorbed_band_soil
    real(r8) :: TR_RHPO4_sorbed_band_soil
    real(r8) :: TR_RH2PO4_sorbed_band_soil
    real(r8) :: TR_AlPO4_precip_soil
    real(r8) :: TR_FePO4_precip_soil
    real(r8) :: TR_CaHPO4_precip_soil
    real(r8) :: TR_apatite_precip_soil
    real(r8) :: TR_CaH4P2O8_precip_soil
    real(r8) :: TR_AlPO4_precip_band_soil
    real(r8) :: TR_FePO4_precip_band_soil
    real(r8) :: TR_CaHPO4_precip_band_soil
    real(r8) :: TR_apatite_precip_band_soil
    real(r8) :: TR_CaH4P2O8_precip_band_soil
    real(r8) :: TR_Al_3p_soil
    real(r8) :: TR_Fe_3p_soil
    real(r8) :: TR_H_p_soil
    real(r8) :: TR_Ca_2p_soil
    real(r8) :: TR_Mg_2p_soil
    real(r8) :: TR_Na_p_soil
    real(r8) :: TR_K_1p_soil
    real(r8) :: TR_OH_1e_soil
    real(r8) :: TR_SO4_2e_soil
    real(r8) :: TR_CO3_2e_soil
    real(r8) :: TR_HCO3
    real(r8) :: TR_CO2_aqu_soil_vr
    real(r8) :: TR_AlOH_soil
    real(r8) :: TR_AlO2H2_soil
    real(r8) :: TR_AlO3H3_soil
    real(r8) :: TR_AlO4H4_soil
    real(r8) :: TR_AlSO4_soil
    real(r8) :: TR_FeOH_soil
    real(r8) :: TR_FeO2H2_soil
    real(r8) :: TR_FeO3H3_soil
    real(r8) :: TR_FeO4H4_soil
    real(r8) :: TR_FeSO4_soil
    real(r8) :: TR_CaOH_soil
    real(r8) :: TR_CaCO3_soil
    real(r8) :: TR_CaHCO3_soil
    real(r8) :: TR_CaSO4_soil
    real(r8) :: TR_MgOH_soil
    real(r8) :: TR_MgCO3_soil
    real(r8) :: TR_MgHCO3_soil
    real(r8) :: TR_MgSO4_soil
    real(r8) :: TR_NaCO3_soil
    real(r8) :: TR_NaSO4_soil
    real(r8) :: TR_KSO4_soil
    real(r8) :: TR_PO4_soil
    real(r8) :: TR_H3PO4_sorbed_soil
    real(r8) :: TR_FeHPO4_soil
    real(r8) :: TR_FeH2PO4_soil
    real(r8) :: TR_CaPO4_soil
    real(r8) :: TR_CaHPO4_soil
    real(r8) :: TR_CaH4P2O8_soil
    real(r8) :: TR_MgHPO4_soil
    real(r8) :: TR_PO4_band_soil
    real(r8) :: TR_H3PO4_band_soil
    real(r8) :: TR_FeHPO4_band_soil
    real(r8) :: TR_FeH2PO4_band_soil
    real(r8) :: TR_CaPO4_band_soil
    real(r8) :: TR_CaHPO4_band_soil
    real(r8) :: TR_CaH4P2O8_band_soil
    real(r8) :: TR_MgHPO4_band_soil
    real(r8) :: TR_H_p_sorbed_soil
    real(r8) :: TR_Al_sorbed_soil
    real(r8) :: TR_Fe_sorbed_soil
    real(r8) :: TR_Ca_sorbed_soil
    real(r8) :: TR_Mg_sorbed_soil
    real(r8) :: TR_Na_sorbed_soil
    real(r8) :: TR_K_sorbed_soil
    real(r8) :: TR_HCO3_sorbed_soil
    real(r8) :: TR_AlO2H2_sorbed_soil
    real(r8) :: TR_FeO2H2_sorbed_soil
    real(r8) :: TR_RO_sorbed_soil
    real(r8) :: TR_RO_sorbed_band_soil
    real(r8) :: TR_AlOH3_precip_soil
    real(r8) :: TR_FeOH3_precip_soil
    real(r8) :: TR_CaCO3_precip_soil
    real(r8) :: TR_CaSO4_precip_soil
    real(r8) :: TRH2O
    real(r8) :: TBION
    real(r8) :: Txchem_CO2
  contains
    procedure, public :: SetZero
  end type solute_flx_type

  type, public :: chem_var_type
  real(r8) :: H1PO4_2e_conc
  real(r8) :: H1PO4_2e_band_conc
  real(r8) :: H2PO4_1e_conc
  real(r8) :: H2PO4_1e_band_conc
  real(r8) :: NH3_aqu_conc
  real(r8) :: NH3_aqu_band_conc
  real(r8) :: NH4_1p_conc
  real(r8) :: NH4_1p_band_conc
  real(r8) :: Precp_AlPO4_conc
  real(r8) :: PrecpB_AlPO4_conc
  real(r8) :: Precp_CaHPO4_conc
  real(r8) :: PrecpB_CaHPO4_conc
  real(r8) :: Precp_Ca5P3O12O3H3_conc
  real(r8) :: PrecpB_Ca5P3O12O3H3_conc
  real(r8) :: Precp_CaH4P2O8_conc
  real(r8) :: PrecpB_CaH4P2O8_con
  real(r8) :: Precp_FePO4_conc
  real(r8) :: PrecpB_FePO4_con
  real(r8) :: SoilMicPMassLayerX
  real(r8) :: VLWatMicPNB
  real(r8) :: VLWatMicPNH
  real(r8) :: VLWatMicPPB
  real(r8) :: VLWatMicPPO
  real(r8) :: XCEC
  real(r8) :: PH
  real(r8) :: CAL
  real(r8) :: CFE
  real(r8) :: VLWatMicPM
  real(r8) :: ZMG
  real(r8) :: ZNA
  real(r8) :: ZKA
  real(r8) :: CCO2S
  real(r8) :: CCA
  real(r8) :: SoilMicPMassLayer
  real(r8) :: XAEC
  real(r8) :: GKC4
  real(r8) :: GKCA
  real(r8) :: GKCH
  real(r8) :: GKCM
  real(r8) :: GKCK
  real(r8) :: GKCN
  real(r8) :: ZNO3S
  real(r8) :: ZNO3B
  real(r8) :: XZHYS
  real(r8) :: ZHY
  real(r8) :: ZOH
  real(r8) :: ZAL
  real(r8) :: ZFE
  real(r8) :: ZCA
  real(r8) :: ZSO4
  real(r8) :: ZCL
  real(r8) :: ZCO3
  real(r8) :: ZHCO3
  real(r8) :: CO2S
  real(r8) :: ZALOH1
  real(r8) :: ZALOH2
  real(r8) :: ZALOH3
  real(r8) :: ZALOH4
  real(r8) :: ZALS
  real(r8) :: ZFEOH1
  real(r8) :: ZFEOH2
  real(r8) :: ZFEOH3
  real(r8) :: ZFEOH4
  real(r8) :: ZFES
  real(r8) :: ZCAO
  real(r8) :: ZCAC
  real(r8) :: ZCAH
  real(r8) :: ZCAS
  real(r8) :: ZMGO
  real(r8) :: ZMGC
  real(r8) :: ZMGH
  real(r8) :: ZMGS
  real(r8) :: ZNAC
  real(r8) :: ZNAS
  real(r8) :: ZKAS
  real(r8) :: H0PO4
  real(r8) :: H3PO4
  real(r8) :: ZFE1P
  real(r8) :: ZFE2P
  real(r8) :: ZCA0P
  real(r8) :: ZCA1P
  real(r8) :: ZCA2P
  real(r8) :: ZMG1P
  real(r8) :: H0POB
  real(r8) :: H3POB
  real(r8) :: ZFE1PB
  real(r8) :: ZFE2PB
  real(r8) :: ZCA0PB
  real(r8) :: ZCA1PB
  real(r8) :: ZCA2PB
  real(r8) :: ZMG1PB
  real(r8) :: XHY
  real(r8) :: XAL
  real(r8) :: XFE
  real(r8) :: XCA
  real(r8) :: XMG
  real(r8) :: XNA
  real(r8) :: XKA
  real(r8) :: XHC
  real(r8) :: XALO2
  real(r8) :: XFEO2
  real(r8) :: ORGC
  real(r8) :: PALOH
  real(r8) :: PFEOH
  real(r8) :: PCACO
  real(r8) :: PCASO
  real(r8) :: VLPOB
  real(r8) :: VLPO4
  real(r8) :: VLNOB
  real(r8) :: VLNO3
  real(r8) :: VLNH4
  real(r8) :: VLNHB
  real(r8) :: XHPO4_band_conc
  real(r8) :: XH2PO4_band_conc
  real(r8) :: XROH_band_conc
  real(r8) :: XHPO4_conc
  real(r8) :: XROH2_band_conc
  real(r8) :: XH2PO4_conc
  real(r8) :: XNH4_conc
  real(r8) :: XNH4_band_conc
  real(r8) :: XROH1_conc
  real(r8) :: XROH2_conc
  real(r8) :: XROH1_band_conc
  real(r8) :: XOH_conc
  real(r8) :: VLWatMicPNZ
  real(r8) :: VLWatMicPNO
  real(r8) :: BKVLNH
  real(r8) :: BKVLNB
  real(r8) :: BKVLPO
  real(r8) :: BKVLPB
  real(r8) :: BKVLNO
  real(r8) :: BKVLNZ
  end type chem_var_type
contains

  subroutine SetZero(solflx)
  implicit none
  class(solute_flx_type)  :: solflx

  solflx%TR_NH4_soil = 0._r8
  solflx%TR_NH4_band_soil = 0._r8
  solflx%TR_NH3_soil_vr = 0._r8
  solflx%TR_NH3_band_soil = 0._r8
  solflx%TR_H1PO4_soil = 0._r8
  solflx%TR_H2PO4_soil = 0._r8
  solflx%TR_H1PO4_band_soil = 0._r8
  solflx%TR_H2PO4_band_soil = 0._r8
  solflx%TR_NH4_sorbed_soil = 0._r8
  solflx%TR_NH4_sorbed_band_soil = 0._r8
  solflx%TR_ROH_sorbed_soil = 0._r8
  solflx%TR_ROH2_sorbed_soil = 0._r8
  solflx%TR_RHPO4_sorbed_soil = 0._r8
  solflx%TR_RH2PO4_sorbed_soil = 0._r8
  solflx%TR_ROH_sorbed_band_soil = 0._r8
  solflx%TR_ROH2_sorbed_band_soil = 0._r8
  solflx%TR_RHPO4_sorbed_band_soil = 0._r8
  solflx%TR_RH2PO4_sorbed_band_soil = 0._r8
  solflx%TR_AlPO4_precip_soil= 0._r8
  solflx%TR_FePO4_precip_soil= 0._r8
  solflx%TR_CaHPO4_precip_soil= 0._r8
  solflx%TR_apatite_precip_soil= 0._r8
  solflx%TR_CaH4P2O8_precip_soil= 0._r8
  solflx%TR_AlPO4_precip_band_soil= 0._r8
  solflx%TR_FePO4_precip_band_soil= 0._r8
  solflx%TR_CaHPO4_precip_band_soil= 0._r8
  solflx%TR_apatite_precip_band_soil= 0._r8
  solflx%TR_CaH4P2O8_precip_band_soil= 0._r8
  solflx%TR_Al_3p_soil  = 0._r8
  solflx%TR_Fe_3p_soil  = 0._r8
  solflx%TR_H_p_soil  = 0._r8
  solflx%TR_Ca_2p_soil  = 0._r8
  solflx%TR_Mg_2p_soil  = 0._r8
  solflx%TR_Na_p_soil  = 0._r8
  solflx%TR_K_1p_soil  = 0._r8
  solflx%TR_OH_1e_soil  = 0._r8
  solflx%TR_SO4_2e_soil = 0._r8
  solflx%TR_CO3_2e_soil = 0._r8
  solflx%TR_HCO3 = 0._r8
  solflx%TR_CO2_aqu_soil_vr = 0._r8
  solflx%TR_AlOH_soil = 0._r8
  solflx%TR_AlO2H2_soil = 0._r8
  solflx%TR_AlO3H3_soil = 0._r8
  solflx%TR_AlO4H4_soil = 0._r8
  solflx%TR_AlSO4_soil = 0._r8
  solflx%TR_FeOH_soil = 0._r8
  solflx%TR_FeO2H2_soil = 0._r8
  solflx%TR_FeO3H3_soil = 0._r8
  solflx%TR_FeO4H4_soil = 0._r8
  solflx%TR_FeSO4_soil = 0._r8
  solflx%TR_CaOH_soil = 0._r8
  solflx%TR_CaCO3_soil = 0._r8
  solflx%TR_CaHCO3_soil = 0._r8
  solflx%TR_CaSO4_soil = 0._r8
  solflx%TR_MgOH_soil = 0._r8
  solflx%TR_MgCO3_soil = 0._r8
  solflx%TR_MgHCO3_soil = 0._r8
  solflx%TR_MgSO4_soil = 0._r8
  solflx%TR_NaCO3_soil = 0._r8
  solflx%TR_NaSO4_soil = 0._r8
  solflx%TR_KSO4_soil = 0._r8
  solflx%TR_PO4_soil = 0._r8
  solflx%TR_H3PO4_sorbed_soil = 0._r8
  solflx%TR_FeHPO4_soil = 0._r8
  solflx%TR_FeH2PO4_soil = 0._r8
  solflx%TR_CaPO4_soil = 0._r8
  solflx%TR_CaHPO4_soil = 0._r8
  solflx%TR_CaH4P2O8_soil = 0._r8
  solflx%TR_MgHPO4_soil = 0._r8
  solflx%TR_PO4_band_soil = 0._r8
  solflx%TR_H3PO4_band_soil = 0._r8
  solflx%TR_FeHPO4_band_soil = 0._r8
  solflx%TR_FeH2PO4_band_soil = 0._r8
  solflx%TR_CaPO4_band_soil = 0._r8
  solflx%TR_CaHPO4_band_soil = 0._r8
  solflx%TR_CaH4P2O8_band_soil = 0._r8
  solflx%TR_MgHPO4_band_soil = 0._r8
  solflx%TR_H_p_sorbed_soil = 0._r8
  solflx%TR_Al_sorbed_soil = 0._r8
  solflx%TR_Fe_sorbed_soil = 0._r8
  solflx%TR_Ca_sorbed_soil = 0._r8
  solflx%TR_Mg_sorbed_soil = 0._r8
  solflx%TR_Na_sorbed_soil = 0._r8
  solflx%TR_K_sorbed_soil = 0._r8
  solflx%TR_HCO3_sorbed_soil = 0._r8
  solflx%TR_AlO2H2_sorbed_soil= 0._r8
  solflx%TR_FeO2H2_sorbed_soil= 0._r8
  solflx%TR_RO_sorbed_soil = 0._r8
  solflx%TR_RO_sorbed_band_soil = 0._r8
  solflx%TR_AlOH3_precip_soil= 0._r8
  solflx%TR_FeOH3_precip_soil= 0._r8
  solflx%TR_CaCO3_precip_soil= 0._r8
  solflx%TR_CaSO4_precip_soil= 0._r8
  solflx%TRH2O = 0._r8
  solflx%TBION = 0._r8
  solflx%Txchem_CO2 = 0._r8

  end subroutine SetZero
end module SoluteChemDataType
