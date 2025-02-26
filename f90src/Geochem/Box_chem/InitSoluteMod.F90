module InitSoluteMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : isclose,AZMAX1
  use SoluteChemDataType, only : solutedtype
  use AqueChemDatatype, only : trcSaltIonNumber
  use ChemTracerParsMod
  use SoluteParMod
  use TracerPropMod
  use TracerIDMod
  use EcoSiMParDataMod, only : micpar
  use EcoSIMCtrlMod, only : salt_model
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  real(r8), parameter :: COOH1=2.5E-02_r8
  real(r8), parameter :: TAD=5.0E-02_r8
  real(r8), parameter :: CALMX=10.0_r8
  real(r8), parameter :: CFEMX=10.0_r8
  real(r8), parameter :: ZERO=1.0E-15_r8

  real(r8) :: A1,A2,CEC_conc,XCOOH_conc,SPOH2,SPOH1,NO3_1e_aqua_mole_conc
  real(r8), pointer :: H2CO3_aqua_mole_conc
  real(r8), pointer :: CH4_aqua_mole_conc
  real(r8), pointer :: O2_aqua_mole_conc
  real(r8), pointer :: N2_aqua_mole_conc
  real(r8), pointer :: N2O_aqua_mole_conc
  real(r8), pointer :: NH4_1p_aqua_mole_conc
  real(r8), pointer :: NH3_aqua_mole_conc
  real(r8), pointer :: Al_3p_aqua_mole_conc
  real(r8), pointer :: Fe_3p_aqua_mole_conc
  real(r8), pointer :: H_1p_aqua_mole_conc
  real(r8), pointer :: Ca_2p_aqua_mole_conc
  real(r8), pointer :: Mg_2p_aqua_mole_conc
  real(r8), pointer :: Na_1p_aqua_mole_conc
  real(r8), pointer :: K_1p_aqua_mole_conc
  real(r8), pointer :: OH_1e_aqua_mole_conc
  real(r8), pointer :: SO4_2e_aqua_mole_conc
  real(r8), pointer :: Cl_e_conc
  real(r8), pointer :: CO3_2e_aqua_mole_conc
  real(r8), pointer :: HCO3_e_conc
  real(r8), pointer :: AlOH_2p_aqua_mole_conc
  real(r8), pointer :: AlO2H2_1p_aqua_mole_conc
  real(r8), pointer :: AlO3H3_conc
  real(r8), pointer :: AlO4H4_1e_aqua_mole_conc
  real(r8), pointer :: AlSO4_1p_aqua_mole_conc
  real(r8), pointer :: FeOH_2p_aqua_mole_conc
  real(r8), pointer :: FeO2H2_p_conc
  real(r8), pointer :: FeO3H3_conc
  real(r8), pointer :: FeO4H4_1e_aqua_mole_conc
  real(r8), pointer :: FeSO4_1p_aqua_mole_conc
  real(r8), pointer :: CaO2H2_conc
  real(r8), pointer :: CaCO3_conc
  real(r8), pointer :: CaHCO3_1p_aqua_mole_conc
  real(r8), pointer :: CaSO4_conc
  real(r8), pointer :: MgOH_1p_aqua_mole_conc
  real(r8), pointer :: MgCO3_conc
  real(r8), pointer :: MgHCO3_1p_aqua_mole_conc
  real(r8), pointer :: MgSO4_conc
  real(r8), pointer :: NaCO3_1e_aqua_mole_conc
  real(r8), pointer :: NaSO4_1e_aqua_mole_conc
  real(r8), pointer :: KSO4_1e_aqua_mole_conc
  real(r8), pointer :: H0PO4_3e_conc
  real(r8), pointer :: H1PO4_2e_aqua_mole_conc
  real(r8), pointer :: H2PO4_1e_aqua_mole_conc
  real(r8), pointer :: H3PO4_conc
  real(r8), pointer :: FeHPO4_p_conc
  real(r8), pointer :: FeH2PO4_2p_aqua_mole_conc
  real(r8), pointer :: CaPO4_1e_con
  real(r8), pointer :: CaHPO4_conc
  real(r8), pointer :: CaH4P2O8_1p_aqua_mole_conc
  real(r8), pointer :: MgHPO4_conc
  real(r8), pointer :: CSTR1
  real(r8), pointer :: CCO2M
  real(r8), pointer :: CCH4M
  real(r8), pointer :: COXYM
  real(r8), pointer :: CZ2GM
  real(r8), pointer :: CZ2OM
  real(r8), pointer :: CN4Z
  real(r8), pointer :: CNOZ
  real(r8), pointer :: CNAZ
  real(r8), pointer :: CKAZ
  real(r8), pointer :: CSOZ
  real(r8), pointer :: CCLZ
  real(r8), pointer :: CNOX
  real(r8), pointer :: CCASOX
  real(r8), pointer :: CN4X
  real(r8), pointer :: CPOZ
  real(r8), pointer :: Al_mole_conc
  real(r8), pointer :: CFEZ
  real(r8), pointer :: CCAZ
  real(r8), pointer :: CMGZ
  real(r8), pointer :: CALX
  real(r8), pointer :: CFEX
  real(r8), pointer :: CaX_conc
  real(r8), pointer :: MgX_conc
  real(r8), pointer :: CNAX
  real(r8), pointer :: CKAX
  real(r8), pointer :: CSOX
  real(r8), pointer :: CCLX
  real(r8), pointer :: CALPOX
  real(r8), pointer :: CFEPOX
  real(r8), pointer :: CCAPDX
  real(r8), pointer :: CCAPHX
  real(r8), pointer :: CALOHX
  real(r8), pointer :: CFEOHX
  real(r8), pointer :: CCACOX
  real(r8), pointer :: XNH4_mole_conc
  real(r8), pointer :: XHY1
  real(r8), pointer :: XAl_conc
  real(r8), pointer :: XFe_conc
  real(r8), pointer :: XCa_conc
  real(r8), pointer :: Precp_Ca5P3O12O3H3_conc
  real(r8), pointer :: XMg_conc
  real(r8), pointer :: XNa_conc
  real(r8), pointer :: XK_conc
  real(r8), pointer :: XHC1
  real(r8), pointer :: XAlO2H2_conc
  real(r8), pointer :: XFeO2H2_conc
  real(r8), pointer :: XOH_conc
  real(r8), pointer :: XROH1_conc
  real(r8), pointer :: XROH2_conc
  real(r8), pointer :: XHPO4_conc
  real(r8), pointer :: XH2PO4_conc
  real(r8), pointer :: Precp_AlO3H3_conc
  real(r8), pointer :: Precp_FeO3H3_conc
  real(r8), pointer :: Precp_CaCO3_conc
  real(r8), pointer :: Precp_CaSO4_conc
  real(r8), pointer :: Precp_AlPO4_conc
  real(r8), pointer :: Precp_FePO4_conc
  real(r8), pointer :: Precp_CaHPO4_conc
  real(r8), pointer :: FH2O
  real(r8), pointer :: ATCA
  real(r8), pointer :: XAEC
  real(r8), pointer :: CEC
  real(r8), pointer :: ORGC
  real(r8), pointer :: VLPO4
  real(r8), pointer :: XCEC
  real(r8), pointer :: GKC4
  real(r8), pointer :: GKCA
  real(r8), pointer :: GKCH
  real(r8), pointer :: GKCK
  real(r8), pointer :: GKCN
  real(r8), pointer :: GKCM
  real(r8), pointer :: ZEROS
  real(r8), pointer :: VLWatMicP

  public :: InitSoluteModel
  public :: InitSoluteProperty
  contains

  subroutine InitSoluteProperty
  implicit none

  trcSaltIonNumber(idsalt_Al)        = 1._r8
  trcSaltIonNumber(idsalt_Fe)        = 1._r8
  trcSaltIonNumber(idsalt_Hp)        = 1._r8
  trcSaltIonNumber(idsalt_Ca)        = 1._r8
  trcSaltIonNumber(idsalt_Mg)        = 1._r8
  trcSaltIonNumber(idsalt_Na)        = 1._r8
  trcSaltIonNumber(idsalt_K)         = 1._r8
  trcSaltIonNumber(idsalt_OH)        = 1._r8
  trcSaltIonNumber(idsalt_SO4)       = 1._r8
  trcSaltIonNumber(idsalt_Cl)        = 1._r8
  trcSaltIonNumber(idsalt_CO3)       = 1._r8
  trcSaltIonNumber(idsalt_H0PO4)     = 1._r8
  trcSaltIonNumber(idsalt_H0PO4B)    = 1._r8
  trcSaltIonNumber(idsalt_HCO3)      = 2._r8
  trcSaltIonNumber(idsalt_AlOH)      = 2._r8
  trcSaltIonNumber(idsalt_AlSO4)     = 2._r8
  trcSaltIonNumber(idsalt_FeOH)      = 2._r8
  trcSaltIonNumber(idsalt_FeSO4)     = 2._r8
  trcSaltIonNumber(idsalt_CaOH)      = 2._r8
  trcSaltIonNumber(idsalt_CaCO3)     = 2._r8
  trcSaltIonNumber(idsalt_CaSO4)     = 2._r8
  trcSaltIonNumber(idsalt_MgOH2)     = 3._r8
  trcSaltIonNumber(idsalt_MgCO3)     = 2._r8
  trcSaltIonNumber(idsalt_MgSO4)     = 2._r8
  trcSaltIonNumber(idsalt_NaCO3)     = 2._r8
  trcSaltIonNumber(idsalt_NaSO4)     = 2._r8
  trcSaltIonNumber(idsalt_KSO4)      = 2._r8
  trcSaltIonNumber(idsalt_CaPO4)     = 2._r8
  trcSaltIonNumber(idsalt_CaPO4B)    = 2._r8
  trcSaltIonNumber(idsalt_AlOH2)     = 3._r8
  trcSaltIonNumber(idsalt_FeOH2)     = 3._r8
  trcSaltIonNumber(idsalt_CaHCO3)    = 3._r8
  trcSaltIonNumber(idsalt_MgHCO3)    = 3._r8
  trcSaltIonNumber(idsalt_FeHPO4)    = 3._r8
  trcSaltIonNumber(idsalt_CaHPO4)    = 3._r8
  trcSaltIonNumber(idsalt_MgHPO4)    = 3._r8
  trcSaltIonNumber(idsalt_FeHPO4B)   = 3._r8
  trcSaltIonNumber(idsalt_CaHPO4B)   = 3._r8
  trcSaltIonNumber(idsalt_MgHPO4B)   = 3._r8
  trcSaltIonNumber(idsalt_AlOH3)     = 4._r8
  trcSaltIonNumber(idsalt_FeOH3)     = 4._r8
  trcSaltIonNumber(idsalt_H3PO4)     = 4._r8
  trcSaltIonNumber(idsalt_FeH2PO4)   = 4._r8
  trcSaltIonNumber(idsalt_CaH4P2O8)  = 4._r8
  trcSaltIonNumber(idsalt_H3PO4B)    = 4._r8
  trcSaltIonNumber(idsalt_FeH2PO4B)  = 4._r8
  trcSaltIonNumber(idsalt_CaH4P2O8B) = 4._r8
  trcSaltIonNumber(idsalt_AlOH4)     = 5._r8
  trcSaltIonNumber(idsalt_FeOH4)     = 5._r8

  end subroutine InitSoluteProperty
!------------------------------------------------------------------------------------------

  SUBROUTINE InitSoluteModel(K,BulkSoilMass,solutevar)
!
!     THIS SUBROUTINE INITIALIZES ALL SOIL CHEMISTRY VARIABLES
!

  implicit none
  integer, intent(in) :: K
  real(r8), intent(in) :: BulkSoilMass
  type(solutedtype), target, intent(inout) :: solutevar

  integer, parameter :: MRXN1=1000
  integer :: M
!     begin_execution
  H2CO3_aqua_mole_conc       => solutevar%H2CO3_aqua_mole_conc
  CH4_aqua_mole_conc         => solutevar%CH4_aqua_mole_conc
  O2_aqua_mole_conc          => solutevar%O2_aqua_mole_conc
  N2_aqua_mole_conc          => solutevar%N2_aqua_mole_conc
  N2O_aqua_mole_conc         => solutevar%N2O_aqua_mole_conc
  NH4_1p_aqua_mole_conc      => solutevar%NH4_1p_aqua_mole_conc
  NH3_aqua_mole_conc         => solutevar%NH3_aqua_mole_conc
  Al_3p_aqua_mole_conc       => solutevar%Al_3p_aqua_mole_conc
  Fe_3p_aqua_mole_conc       => solutevar%Fe_3p_aqua_mole_conc
  H_1p_aqua_mole_conc        => solutevar%H_1p_aqua_mole_conc
  Ca_2p_aqua_mole_conc       => solutevar%Ca_2p_aqua_mole_conc
  Mg_2p_aqua_mole_conc       => solutevar%Mg_2p_aqua_mole_conc
  Na_1p_aqua_mole_conc       => solutevar%Na_1p_aqua_mole_conc
  K_1p_aqua_mole_conc        => solutevar%K_1p_aqua_mole_conc
  OH_1e_aqua_mole_conc       => solutevar%OH_1e_aqua_mole_conc
  SO4_2e_aqua_mole_conc      => solutevar%SO4_2e_aqua_mole_conc
  Cl_e_conc                  => solutevar%Cl_e_conc
  CO3_2e_aqua_mole_conc      => solutevar%CO3_2e_aqua_mole_conc
  HCO3_e_conc                => solutevar%HCO3_e_conc
  AlOH_2p_aqua_mole_conc     => solutevar%AlOH_2p_aqua_mole_conc
  AlO2H2_1p_aqua_mole_conc   => solutevar%AlO2H2_1p_aqua_mole_conc
  AlO3H3_conc                => solutevar%AlO3H3_conc
  AlO4H4_1e_aqua_mole_conc   => solutevar%AlO4H4_1e_aqua_mole_conc
  AlSO4_1p_aqua_mole_conc    => solutevar%AlSO4_1p_aqua_mole_conc
  FeOH_2p_aqua_mole_conc     => solutevar%FeOH_2p_aqua_mole_conc
  FeO2H2_p_conc              => solutevar%FeO2H2_p_conc
  FeO3H3_conc                => solutevar%FeO3H3_conc
  FeO4H4_1e_aqua_mole_conc   => solutevar%FeO4H4_1e_aqua_mole_conc
  FeSO4_1p_aqua_mole_conc    => solutevar%FeSO4_1p_aqua_mole_conc
  CaO2H2_conc                => solutevar%CaO2H2_conc
  CaCO3_conc                 => solutevar%CaCO3_conc
  CaHCO3_1p_aqua_mole_conc   => solutevar%CaHCO3_1p_aqua_mole_conc
  CaSO4_conc                 => solutevar%CaSO4_conc
  MgOH_1p_aqua_mole_conc     => solutevar%MgOH_1p_aqua_mole_conc
  MgCO3_conc                 => solutevar%MgCO3_conc
  MgHCO3_1p_aqua_mole_conc   => solutevar%MgHCO3_1p_aqua_mole_conc
  MgSO4_conc                 => solutevar%MgSO4_conc
  NaCO3_1e_aqua_mole_conc    => solutevar%NaCO3_1e_aqua_mole_conc
  NaSO4_1e_aqua_mole_conc    => solutevar%NaSO4_1e_aqua_mole_conc
  KSO4_1e_aqua_mole_conc     => solutevar%KSO4_1e_aqua_mole_conc
  H0PO4_3e_conc              => solutevar%H0PO4_3e_conc
  H1PO4_2e_aqua_mole_conc    => solutevar%H1PO4_2e_aqua_mole_conc
  H2PO4_1e_aqua_mole_conc    => solutevar%H2PO4_1e_aqua_mole_conc
  H3PO4_conc                 => solutevar%H3PO4_conc
  FeHPO4_p_conc              => solutevar%FeHPO4_p_conc
  FeH2PO4_2p_aqua_mole_conc  => solutevar%FeH2PO4_2p_aqua_mole_conc
  CaPO4_1e_con               => solutevar%CaPO4_1e_con
  CaHPO4_conc                => solutevar%CaHPO4_conc
  CaH4P2O8_1p_aqua_mole_conc => solutevar%CaH4P2O8_1p_aqua_mole_conc
  MgHPO4_conc                => solutevar%MgHPO4_conc
  CSTR1                      => solutevar%CSTR1
  CCO2M                      => solutevar%CCO2M
  CCH4M                      => solutevar%CCH4M
  COXYM                      => solutevar%COXYM
  CZ2GM                      => solutevar%CZ2GM
  CZ2OM                      => solutevar%CZ2OM
  CN4Z                       => solutevar%CN4Z
  CNOZ                       => solutevar%CNOZ
  CNAZ                       => solutevar%CNAZ
  CKAZ                       => solutevar%CKAZ
  CSOZ                       => solutevar%CSOZ
  CCLZ                       => solutevar%CCLZ
  CNOX                       => solutevar%CNOX
  CCASOX                     => solutevar%CCASOX
  CN4X                       => solutevar%CN4X
  CPOZ                       => solutevar%CPOZ
  Al_mole_conc                       => solutevar%Al_mole_conc
  CFEZ                       => solutevar%CFEZ
  CCAZ                       => solutevar%CCAZ
  CMGZ                       => solutevar%CMGZ
  CALX                       => solutevar%CALX
  CFEX                       => solutevar%CFEX
  CaX_conc                   => solutevar%CaX_conc
  MgX_conc                   => solutevar%MgX_conc
  CNAX                       => solutevar%CNAX
  CKAX                       => solutevar%CKAX
  CSOX                       => solutevar%CSOX
  CCLX                       => solutevar%CCLX
  CALPOX                     => solutevar%CALPOX
  CFEPOX                     => solutevar%CFEPOX
  CCAPDX                     => solutevar%CCAPDX
  CCAPHX                     => solutevar%CCAPHX
  CALOHX                     => solutevar%CALOHX
  CFEOHX                     => solutevar%CFEOHX
  CCACOX                     => solutevar%CCACOX
  XNH4_mole_conc                  => solutevar%XNH4_mole_conc
  XHY1                       => solutevar%XHY1
  XAl_conc                   => solutevar%XAl_conc
  XFe_conc                   => solutevar%XFe_conc
  XCa_conc                   => solutevar%XCa_conc
  Precp_Ca5P3O12O3H3_conc    => solutevar%Precp_Ca5P3O12O3H3_conc
  XMg_conc                   => solutevar%XMg_conc
  XNa_conc                   => solutevar%XNa_conc
  XK_conc                    => solutevar%XK_conc
  XHC1                       => solutevar%XHC1
  XAlO2H2_conc               => solutevar%XAlO2H2_conc
  XFeO2H2_conc               => solutevar%XFeO2H2_conc
  XOH_conc                   => solutevar%XOH_conc
  XROH1_conc                 => solutevar%XROH1_conc
  XROH2_conc                 => solutevar%XROH2_conc
  XHPO4_conc                 => solutevar%XHPO4_conc
  XH2PO4_conc                => solutevar%XH2PO4_conc
  Precp_AlO3H3_conc          => solutevar%Precp_AlO3H3_conc
  Precp_FeO3H3_conc          => solutevar%Precp_FeO3H3_conc
  Precp_CaCO3_conc           => solutevar%Precp_CaCO3_conc
  Precp_CaSO4_conc           => solutevar%Precp_CaSO4_conc
  Precp_AlPO4_conc           => solutevar%Precp_AlPO4_conc
  Precp_FePO4_conc           => solutevar%Precp_FePO4_conc
  Precp_CaHPO4_conc          => solutevar%Precp_CaHPO4_conc
  FH2O                       => solutevar%FH2O
  ATCA                       => solutevar%ATCA
  XAEC                       => solutevar%XAEC
  CEC                        => solutevar%CEC
  ORGC                       => solutevar%ORGC
  VLPO4                      => solutevar%VLPO4
  XCEC                       => solutevar%XCEC
  GKC4                       => solutevar%GKC4
  GKCA                       => solutevar%GKCA
  GKCH                       => solutevar%GKCH
  GKCK                       => solutevar%GKCK
  GKCN                       => solutevar%GKCN
  GKCM                       => solutevar%GKCM
  ZEROS                      => solutevar%ZEROS
  VLWatMicP                  => solutevar%VLWatMicP

  call InitEquilibria(K,BulkSoilMass)
!
!     CONVERGE TOWARDS ALL SOLUBILITY EQUILIBRIA
!     IF SALT OPTION IS SELECTED
!
  IF(salt_model)THEN
    DO   M=1,MRXN1
      call SolubilityEquilibiriaSalt(K,M,BulkSoilMass)
    ENDDO
!
!     CONVERGE TOWARDS ALL SOLUBILITY EQUILIBRIA
!     IF SALT OPTION IS NOT SELECTED
!
  ELSE
    DO  M=1,MRXN1
      call SolubilityEquilibriaNoSalt(K,M,BulkSoilMass)
    ENDDO
  ENDIF
  end SUBROUTINE InitSoluteModel
!------------------------------------------------------------------------------------------

  subroutine InitEquilibria(K,BulkSoilMass)

  implicit none
  integer, intent(in) :: K
  real(r8), intent(in) :: BulkSoilMass
  integer :: MM
  real(r8) :: XNAQ,XPT,XN4Q,XMGQ,XKAQ,XHYQ,XFEQ
  real(r8) :: XHP,XCECQ,XTLQ,XALQ,XCAX,SPH2P
  real(r8) :: XCAQ,XOH,SPH1P,FX,FXH2,FXP1,FSTR2
  real(r8) :: FHPA,FXH1,FXP2,FHP1,FHP3,FHCO
  real(r8) :: FCO3,CX1,CX2,CSTRZ,FHP2,CCO2Y
  real(r8) :: CION2,A3,anion_1e_aqua_mole_conc,anion_2e_aqua_mole_conc,anion_3e_conc,cation_1p_aqua_mole_conc,cation_2p_aqua_mole_conc,cation_3p_aqua_mole_conc
  real(r8) :: CCO2X,CCO2Z,CN,CSTR2,FHP0,FXH0,R,Z
!
!     INITIALIZE SOLUTE EQUILIBRIA
!
  cation_3p_aqua_mole_conc=AZMAX1(Al_mole_conc)+AZMAX1(CFEZ)
  anion_3e_conc=0._r8
  cation_2p_aqua_mole_conc=AZMAX1(CCAZ)+CMGZ
  anion_2e_aqua_mole_conc=CSOZ
  cation_1p_aqua_mole_conc=CN4Z+CNAZ+CKAZ
  anion_1e_aqua_mole_conc=CNOZ+CCLZ+CPOZ
  CX2=0._r8
  CX1=0._r8
  CN=0._r8
!
!     INITIALIZE ION STRENGTH AND ACTIVITIES
!
  CION2=AZMAX1(cation_3p_aqua_mole_conc+anion_3e_conc+cation_2p_aqua_mole_conc+anion_2e_aqua_mole_conc+cation_1p_aqua_mole_conc+anion_1e_aqua_mole_conc+CN)
  CSTR1=0.5E-03*(9.0*(cation_3p_aqua_mole_conc+anion_3e_conc)+4.0*(cation_2p_aqua_mole_conc+anion_2e_aqua_mole_conc)+cation_1p_aqua_mole_conc+anion_1e_aqua_mole_conc)
  CSTRZ=0.5E-03*(9.0*(cation_3p_aqua_mole_conc+anion_3e_conc)+4.0*(cation_2p_aqua_mole_conc+CX2)+cation_1p_aqua_mole_conc+CX1)
  CSTR2=SQRT(CSTR1)
  FSTR2=CSTR2/(1.0_r8+CSTR2)
  FH2O=5.56E+04/(5.56E+04+CION2)

  IF(salt_model)THEN
    A1=AMIN1(1.0,10.0**(-0.509*1.0*FSTR2+0.20*CSTR2))
    A2=AMIN1(1.0,10.0**(-0.509*4.0*FSTR2+0.20*CSTR2))
    A3=AMIN1(1.0,10.0**(-0.509*9.0*FSTR2+0.20*CSTR2))
  ELSE
    A1=1.0_r8
    A2=1.0_r8
    A3=1.0_r8
  ENDIF
!
!     INITIALIZE GASES
!

  CCO2X          = CCO2M*gas_solubility(idg_CO2,ATCA)/(EXP(GasSechenovConst(idg_CO2)*CSTRZ))*FH2O
  CCO2Y          = LOG(CCO2X)
  CCO2Z          = ABS(CCO2Y)
  H2CO3_aqua_mole_conc = CCO2X
  FCO3           = DPCO3*A0/(H_1p_aqua_mole_conc**2._r8*A2)
  FHCO           = DPCO2*A0/(H_1p_aqua_mole_conc*A1)
  Z              = GasSechenovConst(idg_CO2)*(2.0E-03_r8*FCO3+0.5E-03_r8*FHCO)
  DO MM   = 1, 25
    R=(LOG(H2CO3_aqua_mole_conc)+Z*H2CO3_aqua_mole_conc-CCO2Y)/CCO2Z
    IF(R.LT.1.0E-03_r8)exit
    H2CO3_aqua_mole_conc=H2CO3_aqua_mole_conc/SQRT(1.0_r8+R)
  ENDDO

  CH4_aqua_mole_conc = CCH4M*gas_solubility(idg_CH4,ATCA)/(EXP(GasSechenovConst(idg_CH4)*CSTR1))*FH2O
  O2_aqua_mole_conc  = COXYM*gas_solubility(idg_O2,ATCA)/(EXP(GasSechenovConst(idg_O2)*CSTR1))*FH2O
  N2_aqua_mole_conc  = CZ2GM*gas_solubility(idg_N2,ATCA)/(EXP(GasSechenovConst(idg_N2)*CSTR1))*FH2O
  N2O_aqua_mole_conc = CZ2OM*gas_solubility(idg_N2O,ATCA)/(EXP(GasSechenovConst(idg_N2O)*CSTR1))*FH2O

  CO3_2e_aqua_mole_conc=H2CO3_aqua_mole_conc*DPCO3*A0/(H_1p_aqua_mole_conc**2._r8*A2)
  HCO3_e_conc=H2CO3_aqua_mole_conc*DPCO2*A0/(H_1p_aqua_mole_conc*A1)
  NO3_1e_aqua_mole_conc=CNOZ
!
!     INITIALIZE ION PAIR EQUILIBRIA
!
  IF(K.NE.micpar%k_POM)THEN
    NH4_1p_aqua_mole_conc=CN4Z/(1.0_r8+DPN4*A1/(H_1p_aqua_mole_conc*A0))
    NH3_aqua_mole_conc=NH4_1p_aqua_mole_conc*DPN4*A1/(H_1p_aqua_mole_conc*A0)
  ELSE
    NH4_1p_aqua_mole_conc=ZERO
    NH3_aqua_mole_conc=ZERO
  ENDIF
  IF(Al_mole_conc.LT.0.0_r8)THEN
    Al_3p_aqua_mole_conc=AMIN1(CALMX,SPALO/(OH_1e_aqua_mole_conc**3*A3))
  ELSE
    Al_3p_aqua_mole_conc=AMIN1(Al_mole_conc,SPALO/(OH_1e_aqua_mole_conc**3*A3))
  ENDIF
  IF(CFEZ.LT.0.0_r8)THEN
    Fe_3p_aqua_mole_conc=AMIN1(CFEMX,SPFEO/(OH_1e_aqua_mole_conc**3*A3))
  ELSE
    Fe_3p_aqua_mole_conc=AMIN1(CFEZ,SPFEO/(OH_1e_aqua_mole_conc**3*A3))
  ENDIF
  IF(CCAZ.LT.0.0_r8)THEN
    Ca_2p_aqua_mole_conc=AMIN1(CCAMX,SPCAC/(CO3_2e_aqua_mole_conc*A2**2))
  ELSE
    Ca_2p_aqua_mole_conc=AMIN1(CCAZ,SPCAC/(CO3_2e_aqua_mole_conc*A2**2))
  ENDIF
  Mg_2p_aqua_mole_conc   = CMGZ
  Na_1p_aqua_mole_conc             = CNAZ
  K_1p_aqua_mole_conc              = CKAZ
  SO4_2e_aqua_mole_conc            = CSOZ
  Cl_e_conc              = CCLZ
  AlOH_2p_aqua_mole_conc = Al_3p_aqua_mole_conc*OH_1e_aqua_mole_conc*A3/(DPAL1*A2)
  AlO2H2_1p_aqua_mole_conc         = Al_3p_aqua_mole_conc*OH_1e_aqua_mole_conc**2*A3/(DPAL1*DPAL2*A1)
  AlO3H3_conc            = Al_3p_aqua_mole_conc*OH_1e_aqua_mole_conc**3*A3/(DPAL1*DPAL2*DPAL3*A0)
  AlO4H4_1e_aqua_mole_conc         = Al_3p_aqua_mole_conc*OH_1e_aqua_mole_conc**4*A3/(DPAL1*DPAL2*DPAL3*DPAL4*A1)
  AlSO4_1p_aqua_mole_conc          = 0._r8
  FeOH_2p_aqua_mole_conc = Fe_3p_aqua_mole_conc*OH_1e_aqua_mole_conc*A3/(DPFE1*A2)
  FeO2H2_p_conc          = Fe_3p_aqua_mole_conc*OH_1e_aqua_mole_conc**2*A3/(DPFE1*DPFE2*A1)
  FeO3H3_conc            = Fe_3p_aqua_mole_conc*OH_1e_aqua_mole_conc**3*A3/(DPFE1*DPFE2*DPFE3*A0)
  FeO4H4_1e_aqua_mole_conc         = Fe_3p_aqua_mole_conc*OH_1e_aqua_mole_conc**4*A3/(DPFE1*DPFE2*DPFE3*DPFE4*A1)
  FeSO4_1p_aqua_mole_conc          = 0._r8
  CaO2H2_conc            = Ca_2p_aqua_mole_conc*OH_1e_aqua_mole_conc*A2/(DPCAO*A1)
  CaCO3_conc=Ca_2p_aqua_mole_conc*CO3_2e_aqua_mole_conc*A2**2/(DPCAC*A0)
  CaHCO3_1p_aqua_mole_conc=Ca_2p_aqua_mole_conc*HCO3_e_conc*A2/DPCAH
  CaSO4_conc=0._r8
  MgOH_1p_aqua_mole_conc=Mg_2p_aqua_mole_conc*OH_1e_aqua_mole_conc*A2/(DPMGO*A1)
  MgCO3_conc=Mg_2p_aqua_mole_conc*CO3_2e_aqua_mole_conc*A2**2/(DPMGC*A0)
  MgHCO3_1p_aqua_mole_conc=Mg_2p_aqua_mole_conc*HCO3_e_conc*A2/DPMGH
  MgSO4_conc=0._r8
  NaCO3_1e_aqua_mole_conc=Na_1p_aqua_mole_conc*CO3_2e_aqua_mole_conc*A2/DPNAC
  NaSO4_1e_aqua_mole_conc=0._r8
  KSO4_1e_aqua_mole_conc=0._r8
  FeHPO4_p_conc=0._r8
  FeH2PO4_2p_aqua_mole_conc=0._r8
  CaPO4_1e_con=0._r8
  CaHPO4_conc=0._r8
  CaH4P2O8_1p_aqua_mole_conc=0._r8
  MgHPO4_conc=0._r8
!
!     INITIALIZE PHOSPHORUS EQUILIBRIA AMONG SOLUBLE, ADSORBED
!     AND PRECIPITATED FORMS
! NOT POM complex
  IF(K.NE.micpar%k_POM)THEN
    H3PO4_conc=CPOZ/(1.0_r8+DPH3P*A0/(H_1p_aqua_mole_conc*A1)+DPH3P*DPH2P*A0 &
      /(H_1p_aqua_mole_conc**2*A2)+DPH3P*DPH2P*DPH1P*A0/(H_1p_aqua_mole_conc**3*A3))
    H2PO4_1e_aqua_mole_conc=H3PO4_conc*DPH3P*A0/(H_1p_aqua_mole_conc*A1)
    H1PO4_2e_aqua_mole_conc=H3PO4_conc*DPH3P*DPH2P*A0/(H_1p_aqua_mole_conc**2*A2)
    H0PO4_3e_conc=H3PO4_conc*DPH3P*DPH2P*DPH1P*A0/(H_1p_aqua_mole_conc**3*A3)
! POM complex
  ELSE
    XHP=CPOZ
    XOH=XAEC/BulkSoilMass
    FHP3=1.0_r8/(1.0_r8+DPH3P*A0/(H_1p_aqua_mole_conc*A1)+DPH3P*DPH2P*A0 &
      /(H_1p_aqua_mole_conc**2*A2)+DPH3P*DPH2P*DPH1P*A0/(H_1p_aqua_mole_conc**3*A3))
    FHP2=FHP3*DPH3P*A0/(H_1p_aqua_mole_conc*A1)
    FHP1=FHP3*DPH3P*DPH2P*A0/(H_1p_aqua_mole_conc**2*A2)
    FHP0=FHP3*DPH3P*DPH2P*DPH1P*A0/(H_1p_aqua_mole_conc**3*A3)
    SPOH2=SXOH2/A1
    SPOH1=SXOH1/A1
    SPH2P=SXH2P*DPH2O/A1
    SPH1P=SXH1P*DPH2O*A1/A2
    FXH2=1.0_r8/(1.0_r8+SPOH2/H_1p_aqua_mole_conc+SPOH2*SPOH1/H_1p_aqua_mole_conc**2)
    FXH1=FXH2*SPOH2/H_1p_aqua_mole_conc
    FXH0=FXH1*SPOH1/H_1p_aqua_mole_conc
    FXP2=1.0_r8/(1.0_r8+SXH2P*DPH2P/(SXH1P*H_1p_aqua_mole_conc))
    FXP1=FXP2*SXH2P*DPH2P/(SXH1P*H_1p_aqua_mole_conc)
    FHPA=FHP2*A1
    XPT=(XOH+XHP+SXH2P*FXP2*OH_1e_aqua_mole_conc/(FXH1*FHPA)-SQRT(XOH**2*FXH1**2 &
      *FHPA**2-2.0*XOH*FXH1**2*XHP*FHPA**2+FXH1**2*XHP**2*FHPA**2 &
      +2.0*XOH*FXH1*FHPA*SXH2P*FXP2*OH_1e_aqua_mole_conc+2.0*FXH1*XHP*FHPA*SXH2P &
      *FXP2*OH_1e_aqua_mole_conc+SXH2P**2*FXP2**2*OH_1e_aqua_mole_conc**2)/(FXH1*FHPA))/2.0_r8
    XROH2_conc=(XOH-XPT)*FXH2
    XROH1_conc=(XOH-XPT)*FXH1
    XOH_conc=(XOH-XPT)*FXH0
    XHPO4_conc=XPT*FXP1
    XH2PO4_conc=XPT*FXP2
    H3PO4_conc=(XHP-XPT)*FHP3
    H2PO4_1e_aqua_mole_conc=(XHP-XPT)*FHP2
    H1PO4_2e_aqua_mole_conc=(XHP-XPT)*FHP1
    H0PO4_3e_conc=(XHP-XPT)*FHP0
!
!     INITIALIZE CATION EQILIBRIA BETWEEN SOLUBLE
!     AND EXCHANGEABLE FORMS
!
    XCECQ=AMAX1(CN4X,CEC)
    XN4Q=CN4X
    XHYQ=0._r8
    XALQ=0._r8
    XFEQ=0._r8
    XCAQ=0._r8
    XMGQ=0._r8
    XNAQ=0._r8
    XKAQ=0._r8
    XHC1=0._r8
    XAlO2H2_conc=0._r8
    XFeO2H2_conc=0._r8
    XCOOH_conc=AZMAX1(COOH1*ORGC)
  ENDIF
  cation_3p_aqua_mole_conc=Al_3p_aqua_mole_conc+Fe_3p_aqua_mole_conc
  anion_3e_conc=H0PO4_3e_conc
  cation_2p_aqua_mole_conc=Ca_2p_aqua_mole_conc+Mg_2p_aqua_mole_conc+AlOH_2p_aqua_mole_conc+FeOH_2p_aqua_mole_conc+FeH2PO4_2p_aqua_mole_conc
  anion_2e_aqua_mole_conc=SO4_2e_aqua_mole_conc+CO3_2e_aqua_mole_conc+H1PO4_2e_aqua_mole_conc
  cation_1p_aqua_mole_conc=NH4_1p_aqua_mole_conc+H_1p_aqua_mole_conc+Na_1p_aqua_mole_conc+K_1p_aqua_mole_conc+AlO2H2_1p_aqua_mole_conc+FeO2H2_p_conc+AlSO4_1p_aqua_mole_conc &
    +FeSO4_1p_aqua_mole_conc+CaO2H2_conc+CaHCO3_1p_aqua_mole_conc+MgOH_1p_aqua_mole_conc+MgHCO3_1p_aqua_mole_conc+FeHPO4_p_conc+CaH4P2O8_1p_aqua_mole_conc
  anion_1e_aqua_mole_conc=NO3_1e_aqua_mole_conc+OH_1e_aqua_mole_conc+HCO3_e_conc+Cl_e_conc+AlO4H4_1e_aqua_mole_conc+FeO4H4_1e_aqua_mole_conc+NaCO3_1e_aqua_mole_conc &
    +NaSO4_1e_aqua_mole_conc+KSO4_1e_aqua_mole_conc+H2PO4_1e_aqua_mole_conc+CaPO4_1e_con
  CN=H2CO3_aqua_mole_conc+CH4_aqua_mole_conc+O2_aqua_mole_conc+N2_aqua_mole_conc+N2O_aqua_mole_conc+NH3_aqua_mole_conc &
    +AlO3H3_conc+FeO3H3_conc+CaCO3_conc+CaSO4_conc &
    +MgCO3_conc+MgSO4_conc+H3PO4_conc+CaHPO4_conc+MgHPO4_conc
  CX2=anion_2e_aqua_mole_conc-CO3_2e_aqua_mole_conc
  CX1=anion_1e_aqua_mole_conc-HCO3_e_conc
!
!     INITIALIZE EQUILIBRIA BETWEEN SOLUBLE AND PRECIPITATED FORMS
! POM complex
  IF(K.EQ.micpar%k_POM)THEN
    Precp_AlO3H3_conc=CALOHX
    Precp_FeO3H3_conc=CFEOHX
    Precp_CaCO3_conc=CCACOX
    Precp_CaSO4_conc=CCASOX
    Precp_AlPO4_conc=CALPOX*VLPO4
    Precp_FePO4_conc=CFEPOX*VLPO4
    Precp_CaHPO4_conc=CCAPDX*VLPO4
    Precp_Ca5P3O12O3H3_conc=CCAPHX*VLPO4
    CEC_conc=AMAX1(ZERO,XCEC/BulkSoilMass)
    CALX=Al_3p_aqua_mole_conc**0.333_r8
    CFEX=Fe_3p_aqua_mole_conc**0.333_r8
    CaX_conc=Ca_2p_aqua_mole_conc**0.500
    MgX_conc=Mg_2p_aqua_mole_conc**0.500
    XCAX=CEC_conc/(1.0_r8+GKC4*NH4_1p_aqua_mole_conc/CaX_conc &
      +GKCH*H_1p_aqua_mole_conc/CaX_conc+GKCA*CALX/CaX_conc &
      +GKCA*CFEX/CaX_conc+GKCA*Mg_2p_aqua_mole_conc/CaX_conc &
      +GKCN*Na_1p_aqua_mole_conc/CaX_conc+GKCK*K_1p_aqua_mole_conc/CaX_conc)
    XN4Q=CN4X
    XHYQ=XCAX*GKCH
    XALQ=XCAX*GKCA
    XFEQ=XCAX*GKCA
    XCAQ=XCAX*CaX_conc
    XMGQ=XCAX*GKCM
    XNAQ=XCAX*GKCN
    XKAQ=XCAX*GKCK
    XTLQ=XN4Q+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ
    IF(XTLQ.GT.ZERO)THEN
      FX=CEC_conc/XTLQ
    ELSE
      FX=0._r8
    ENDIF
    XNH4_mole_conc=CN4X
    XHY1=FX*XHYQ
    XAl_conc=FX*XALQ/3.0_r8
    XFe_conc=FX*XFEQ/3.0_r8
    XCa_conc=FX*XCAQ/2.0_r8
    XMg_conc=FX*XMGQ/2.0_r8
    XNa_conc=FX*XNAQ
    XK_conc=FX*XKAQ
  ENDIF
  end subroutine InitEquilibria
!------------------------------------------------------------------------------------------

  subroutine SolubilityEquilibiriaSalt(K,M,BulkSoilMass)

  implicit none
  integer, intent(in) :: K,M
  real(r8), intent(in) :: BulkSoilMass  !Mg/d-2, soil mass whole layer
  real(r8) :: XN4Q,XMGQ,XKAQ,XHYQ,XFEQ,XTLQ
  real(r8) :: XALQ,XCAX,SPH2P,XCAQ,SPH1P,FX
  real(r8) :: RN4S,RN3S,RNA,RNAC,RNH4
  real(r8) :: RPALOX,H2PO4_1e_AlPO4_dissol_flx,RPCACX,H2PO4_1e_CaHPO4_dissol_flx
  real(r8) :: H2PO4_e_to_HPO4_2e_flx,RHP1,RHP2,H2PO4_1e_apatite_dissol_flx
  real(r8) :: H2PO4_1e_FePO4_dissol_flx,H1PO4_to_XHPO4_ROH_flx
  real(r8) :: FSTR2,RCAO,FHCO,FHP2
  real(r8) :: H2PO4_1e_to_XH2PO4_ROH2_flx,RXN4,RF1P,FCO3,SP
  real(r8) :: XCOO,XNAQ,RHP0,RHP3,RKAS,RM1P
  real(r8) :: CX1,CX2,CSTRZ,CCO2Y,CION2
  real(r8) :: RFE1,RFE2,RFE3,RFE4
  real(r8) :: RFEO1,RFEO2,RFEO3,RFEO4
  real(r8) :: RHA4P1,RHA4P2,RHAL1
  real(r8) :: RHALO1,RHALO2,RHALO3,RHALO4
  real(r8) :: H2PO4_to_XH2PO4_ROH_flx,S0,S1,VLWatMicPBK,RCA,RCAC,RCAH
  real(r8) :: RCAS,RF2P,RFE,RFES,RH1P,RH3P
  real(r8) :: RHA0P1,RHA0P2,RHA1P1,RHA1P2,RHA2P1
  real(r8) :: RHA2P2,RHA3P1,RHCAC3,RHCACH,RHA3P2
  real(r8) :: RHCACO,RHCAH1,RHCAH2,RHF0P1,RHCAD2
  real(r8) :: RHF1P2,RHF2P2,RHF3P1,RHF0P2,RHF3P2
  real(r8) :: RHF1P1,RHF2P1,RHF4P1,RHF4P2,RHFE1
  real(r8) :: RHFEO1,RHFEO2,RHFEO3,RHFEO4,RKA,RMG
  real(r8) :: RMGC,RMGH,RMGO,RMGS,RNAS,RPCAD1,RPCASO
  real(r8) :: RPFEOX,RSO4,RXCA,RXFE,RXHY,RXKA,RXMG
  real(r8) :: RXAL,RXNA,RXOH1,RXOH2,SPX,R,Z
  real(r8) :: A3,Al_3p_activity,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity
  real(r8) :: AlSO4_1p_activity,AALX,CaPO4_1e_activity,CaHPO4_activity,CaH4P2O8_1p_activity,Ca_2p_activity,CaCO3_activity
  real(r8) :: CaHCO3_1p_activity,CaSO4_activity,CaO2H2_activity,ACAX,H2CO3_activity,CO3_2e_activity
  real(r8) :: FeHPO4_p_activity,FeH2PO4_2p_activity,Fe_3p_activity,FeOH_2p_activity
  real(r8) :: FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity
  real(r8) :: FeSO4_1p_activity,AFEX,H0PO4_3e_activity,H1PO4_2e_activity,H2PO4_1e_activity,H3PO4_activity,HCO3_e_activity
  real(r8) :: H_1p_activity,K_1p_activity,KSO4_1e_activity,MgHPO4_activity,Mg_2p_activity,MgCO3_activity,MgOH_1p_activity
  real(r8) :: AMGX,MgHCO3_1p_activity,MgSO4_activity,NH3_activity,NH4_1p_activity,Na_1p_activity,NaCO3_1e_activity
  real(r8) :: OH_1e_activity,SO4_2e_activity,NaSO4_1e_activity,anion_1e_aqua_mole_conc,anion_2e_aqua_mole_conc
  real(r8) :: anion_3e_conc,cation_1p_aqua_mole_conc,cation_2p_aqua_mole_conc
  real(r8) :: cation_3p_aqua_mole_conc,CCO2X,CCO2Z,CN,CSTR2,DP,P1,P2,P3,PX
  real(r8) :: PY,R1,RAL,RAL1,RAL2,RAL3,RAL4,RALO1
  real(r8) :: RALO2,RALO3,RALO4,RALS,RC0P,RC1P,RC2P
  integer :: NR1,NP2,NP3
  integer :: MM

! begin_execution
  H2CO3_aqua_mole_conc       = AMAX1(ZERO,H2CO3_aqua_mole_conc)
  CO3_2e_aqua_mole_conc      = H2CO3_aqua_mole_conc*DPCO3*A0/(H_1p_aqua_mole_conc**2*A2)
  HCO3_e_conc                = H2CO3_aqua_mole_conc*DPCO2*A0/(H_1p_aqua_mole_conc*A1)
  NH4_1p_aqua_mole_conc      = AMAX1(ZERO,NH4_1p_aqua_mole_conc)
  NH3_aqua_mole_conc         = AMAX1(ZERO,NH3_aqua_mole_conc)
  Al_3p_aqua_mole_conc       = AMAX1(ZERO,Al_3p_aqua_mole_conc)
  Fe_3p_aqua_mole_conc       = AMAX1(ZERO,Fe_3p_aqua_mole_conc)
  Ca_2p_aqua_mole_conc       = AMAX1(ZERO,Ca_2p_aqua_mole_conc)
  Ca_2p_aqua_mole_conc       = AMIN1(Ca_2p_aqua_mole_conc,SPCAC/(CO3_2e_aqua_mole_conc*A2**2))
  Mg_2p_aqua_mole_conc       = AMAX1(ZERO,Mg_2p_aqua_mole_conc)
  Na_1p_aqua_mole_conc       = AMAX1(ZERO,Na_1p_aqua_mole_conc)
  K_1p_aqua_mole_conc        = AMAX1(ZERO,K_1p_aqua_mole_conc)
  SO4_2e_aqua_mole_conc      = AMAX1(ZERO,SO4_2e_aqua_mole_conc)
  AlOH_2p_aqua_mole_conc     = AMAX1(ZERO,AlOH_2p_aqua_mole_conc)
  AlO2H2_1p_aqua_mole_conc   = AMAX1(ZERO,AlO2H2_1p_aqua_mole_conc)
  AlO3H3_conc                = AMAX1(ZERO,AlO3H3_conc)
  AlO4H4_1e_aqua_mole_conc   = AMAX1(ZERO,AlO4H4_1e_aqua_mole_conc)
  AlSO4_1p_aqua_mole_conc    = AMAX1(ZERO,AlSO4_1p_aqua_mole_conc)
  FeOH_2p_aqua_mole_conc     = AMAX1(ZERO,FeOH_2p_aqua_mole_conc)
  FeO2H2_p_conc              = AMAX1(ZERO,FeO2H2_p_conc)
  FeO3H3_conc                = AMAX1(ZERO,FeO3H3_conc)
  FeO4H4_1e_aqua_mole_conc   = AMAX1(ZERO,FeO4H4_1e_aqua_mole_conc)
  FeSO4_1p_aqua_mole_conc    = AMAX1(ZERO,FeSO4_1p_aqua_mole_conc)
  CaO2H2_conc                = AMAX1(ZERO,CaO2H2_conc)
  CaCO3_conc                 = AMAX1(ZERO,CaCO3_conc)
  CaHCO3_1p_aqua_mole_conc   = AMAX1(ZERO,CaHCO3_1p_aqua_mole_conc)
  CaSO4_conc                 = AMAX1(ZERO,CaSO4_conc)
  MgOH_1p_aqua_mole_conc     = AMAX1(ZERO,MgOH_1p_aqua_mole_conc)
  MgCO3_conc                 = AMAX1(ZERO,MgCO3_conc)
  MgHCO3_1p_aqua_mole_conc   = AMAX1(ZERO,MgHCO3_1p_aqua_mole_conc)
  MgSO4_conc                 = AMAX1(ZERO,MgSO4_conc)
  NaCO3_1e_aqua_mole_conc    = AMAX1(ZERO,NaCO3_1e_aqua_mole_conc)
  NaSO4_1e_aqua_mole_conc    = AMAX1(ZERO,NaSO4_1e_aqua_mole_conc)
  KSO4_1e_aqua_mole_conc     = AMAX1(ZERO,KSO4_1e_aqua_mole_conc)
  H0PO4_3e_conc              = AMAX1(ZERO,H0PO4_3e_conc)
  H1PO4_2e_aqua_mole_conc    = AMAX1(ZERO,H1PO4_2e_aqua_mole_conc)
  H2PO4_1e_aqua_mole_conc    = AMAX1(ZERO,H2PO4_1e_aqua_mole_conc)
  H3PO4_conc                 = AMAX1(ZERO,H3PO4_conc)
  FeHPO4_p_conc              = AMAX1(ZERO,FeHPO4_p_conc)
  FeH2PO4_2p_aqua_mole_conc  = AMAX1(ZERO,FeH2PO4_2p_aqua_mole_conc)
  CaPO4_1e_con               = AMAX1(ZERO,CaPO4_1e_con)
  CaHPO4_conc                = AMAX1(ZERO,CaHPO4_conc)
  CaH4P2O8_1p_aqua_mole_conc = AMAX1(ZERO,CaH4P2O8_1p_aqua_mole_conc)
  MgHPO4_conc                = AMAX1(ZERO,MgHPO4_conc)
!
!     ION ACTIVITY COEFFICIENTS
!
  cation_3p_aqua_mole_conc=Al_3p_aqua_mole_conc+Fe_3p_aqua_mole_conc
  anion_3e_conc=H0PO4_3e_conc
  cation_2p_aqua_mole_conc=Ca_2p_aqua_mole_conc+Mg_2p_aqua_mole_conc+AlOH_2p_aqua_mole_conc+FeOH_2p_aqua_mole_conc+FeH2PO4_2p_aqua_mole_conc
  anion_2e_aqua_mole_conc=SO4_2e_aqua_mole_conc+CO3_2e_aqua_mole_conc+H1PO4_2e_aqua_mole_conc
  cation_1p_aqua_mole_conc=NH4_1p_aqua_mole_conc+H_1p_aqua_mole_conc+Na_1p_aqua_mole_conc+K_1p_aqua_mole_conc+AlO2H2_1p_aqua_mole_conc &
    +FeO2H2_p_conc+AlSO4_1p_aqua_mole_conc+FeSO4_1p_aqua_mole_conc+CaO2H2_conc &
    +CaHCO3_1p_aqua_mole_conc+MgOH_1p_aqua_mole_conc+MgHCO3_1p_aqua_mole_conc+FeHPO4_p_conc+CaH4P2O8_1p_aqua_mole_conc
  anion_1e_aqua_mole_conc=NO3_1e_aqua_mole_conc+OH_1e_aqua_mole_conc+HCO3_e_conc+Cl_e_conc+AlO4H4_1e_aqua_mole_conc &
    +FeO4H4_1e_aqua_mole_conc+NaCO3_1e_aqua_mole_conc+NaSO4_1e_aqua_mole_conc+KSO4_1e_aqua_mole_conc &
    +H2PO4_1e_aqua_mole_conc+CaPO4_1e_con
  CN=H2CO3_aqua_mole_conc+CH4_aqua_mole_conc+O2_aqua_mole_conc+N2_aqua_mole_conc+N2O_aqua_mole_conc+NH3_aqua_mole_conc &
    +AlO3H3_conc+FeO3H3_conc+CaCO3_conc+CaSO4_conc &
    +MgCO3_conc+MgSO4_conc+H3PO4_conc+CaHPO4_conc+MgHPO4_conc
  CX2=anion_2e_aqua_mole_conc-CO3_2e_aqua_mole_conc
  CX1=anion_1e_aqua_mole_conc-HCO3_e_conc
  CION2=AZMAX1(cation_3p_aqua_mole_conc+anion_3e_conc+cation_2p_aqua_mole_conc+anion_2e_aqua_mole_conc+cation_1p_aqua_mole_conc+anion_1e_aqua_mole_conc+CN)
  CSTR1=0.5E-03_r8*(9.0_r8*(cation_3p_aqua_mole_conc+anion_3e_conc)+4.0*(cation_2p_aqua_mole_conc+anion_2e_aqua_mole_conc)+cation_1p_aqua_mole_conc+anion_1e_aqua_mole_conc)
  CSTRZ=0.5E-03_r8*(9.0_r8*(cation_3p_aqua_mole_conc+anion_3e_conc)+4.0*(cation_2p_aqua_mole_conc+CX2)+cation_1p_aqua_mole_conc+CX1)
  CSTR2=SQRT(CSTR1)
  FSTR2=CSTR2/(1.0_r8+CSTR2)
  FH2O=5.56E+04_r8/(5.56E+04_r8+CION2)
  A1=AMIN1(1.0_r8,10.0_r8**(-0.509_r8*1.0_r8*FSTR2+0.20_r8*CSTR2))
  A2=AMIN1(1.0_r8,10.0_r8**(-0.509_r8*4.0_r8*FSTR2+0.20_r8*CSTR2))
  A3=AMIN1(1.0_r8,10.0_r8**(-0.509_r8*9.0_r8*FSTR2+0.20_r8*CSTR2))
!
!     PRECIPITATION-DISSOLUTION EQUILIBRIA
!
  H_1p_activity=H_1p_aqua_mole_conc*A1
  OH_1e_activity=OH_1e_aqua_mole_conc*A1
  Al_3p_activity=Al_3p_aqua_mole_conc*A3
  AlOH_2p_activity=AlOH_2p_aqua_mole_conc*A2
  AlO2H2_1p_activity=AlO2H2_1p_aqua_mole_conc*A1
  AlO3H3_activity=AlO3H3_conc*A0
  AlO4H4_1e_activity=AlO4H4_1e_aqua_mole_conc*A1
  Fe_3p_activity=Fe_3p_aqua_mole_conc*A3
  FeOH_2p_activity=FeOH_2p_aqua_mole_conc*A2
  FeO2H2_p_activity=FeO2H2_p_conc*A1
  FeO3H3_activity=FeO3H3_conc*A0
  FeO4H4_1e_activity=FeO4H4_1e_aqua_mole_conc*A1
  Ca_2p_activity=Ca_2p_aqua_mole_conc*A2
  CO3_2e_activity=CO3_2e_aqua_mole_conc*A2
  HCO3_e_activity=HCO3_e_conc*A1
  H2CO3_activity=H2CO3_aqua_mole_conc*A0
  SO4_2e_activity=SO4_2e_aqua_mole_conc*A2
  H0PO4_3e_activity=H0PO4_3e_conc*A3
  H1PO4_2e_activity=H1PO4_2e_aqua_mole_conc*A2
  H2PO4_1e_activity=H2PO4_1e_aqua_mole_conc*A1
  H3PO4_activity=H3PO4_conc*A0
  FeHPO4_p_activity=FeHPO4_p_conc*A2
  FeH2PO4_2p_activity=FeH2PO4_2p_aqua_mole_conc*A2
  CaPO4_1e_activity=CaPO4_1e_con*A1
  CaHPO4_activity=CaHPO4_conc*A0
  CaH4P2O8_1p_activity=CaH4P2O8_1p_aqua_mole_conc*A1
  MgHPO4_activity=MgHPO4_conc*A0
  NH4_1p_activity=NH4_1p_aqua_mole_conc*A1
  NH3_activity=NH3_aqua_mole_conc*A0
  Mg_2p_activity=Mg_2p_aqua_mole_conc*A2
  Na_1p_activity=Na_1p_aqua_mole_conc*A1
  K_1p_activity=K_1p_aqua_mole_conc*A1
  AALX=Al_3p_activity**0.333_r8
  AFEX=Fe_3p_activity**0.333_r8
  ACAX=Ca_2p_activity**0.500_r8
  AMGX=Mg_2p_activity**0.500_r8
  AlSO4_1p_activity=AlSO4_1p_aqua_mole_conc*A1
  FeSO4_1p_activity=FeSO4_1p_aqua_mole_conc*A1
  CaO2H2_activity=CaO2H2_conc*A0
  CaCO3_activity=CaCO3_conc*A0
  CaSO4_activity=CaSO4_conc*A0
  CaHCO3_1p_activity=CaHCO3_1p_aqua_mole_conc*A1
  MgOH_1p_activity=MgOH_1p_aqua_mole_conc*A1
  MgCO3_activity=MgCO3_conc*A0
  MgHCO3_1p_activity=MgHCO3_1p_aqua_mole_conc*A1
  MgSO4_activity=MgSO4_conc*A0
  NaCO3_1e_activity=NaCO3_1e_aqua_mole_conc*A1
  NaSO4_1e_activity=NaSO4_1e_aqua_mole_conc*A1
  KSO4_1e_activity=KSO4_1e_aqua_mole_conc*A1
!
!     ALUMINUM HYDROXIDE (GIBBSITE)
!
  IF(K.EQ.micpar%k_POM)THEN
    PX=AMAX1(Al_3p_activity,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity)
    IF(isclose(PX,Al_3p_activity))THEN
      R1=H_1p_activity
      P1=Al_3p_activity
      P2=OH_1e_activity
      NR1=3
      NP2=0
      SP=SHALO
    ELSEIF(isclose(PX,AlOH_2p_activity))THEN
      R1=H_1p_activity
      P1=AlOH_2p_activity
      P2=OH_1e_activity
      NR1=2
      NP2=0
      SP=SHAL1
    ELSEIF(isclose(PX,AlO2H2_1p_activity))THEN
      R1=H_1p_activity
      P1=AlO2H2_1p_activity
      P2=OH_1e_activity
      NR1=1
      NP2=0
      SP=SHAL2
    ELSEIF(isclose(PX,AlO3H3_activity))THEN
      R1=H_1p_activity
      P1=AlO3H3_activity
      P2=OH_1e_activity
      NR1=0
      NP2=0
      SP=SPAL3
    ELSEIF(isclose(PX,AlO4H4_1e_activity))THEN
      R1=OH_1e_activity
      P1=AlO4H4_1e_activity
      P2=H_1p_activity
      NR1=0
      NP2=1
      SP=SHAL4
    ENDIF
    RHAL1=0._r8
    RHALO1=0._r8
    RHALO2=0._r8
    RHALO3=0._r8
    RHALO4=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1/P2**NP2
    RPALOX=AMAX1(-Precp_AlO3H3_conc,TPD*(P1-SPX))
    IF(isclose(PX,Al_3p_activity))THEN
      RHAL1=RPALOX
    ELSEIF(isclose(PX,AlOH_2p_activity))THEN
      RHALO1=RPALOX
    ELSEIF(isclose(PX,AlO2H2_1p_activity))THEN
      RHALO2=RPALOX
    ELSEIF(isclose(PX,AlO3H3_activity))THEN
      RHALO3=RPALOX
    ELSEIF(isclose(PX,AlO4H4_1e_activity))THEN
      RHALO4=RPALOX
    ENDIF
!
!     IRON HYDROXIDE
!
    PX=AMAX1(Fe_3p_activity,FeOH_2p_activity,FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity)
    IF(isclose(PX,Fe_3p_activity))THEN
      R1=H_1p_activity
      P1=Fe_3p_activity
      P2=OH_1e_activity
      NR1=3
      NP2=0
      SP=SHFEO
    ELSEIF(isclose(PX,FeOH_2p_activity))THEN
      R1=H_1p_activity
      P1=FeOH_2p_activity
      P2=OH_1e_activity
      NR1=2
      NP2=0
      SP=SHFE1
    ELSEIF(isclose(PX,FeO2H2_p_activity))THEN
      R1=H_1p_activity
      P1=FeO2H2_p_activity
      P2=OH_1e_activity
      NR1=1
      NP2=0
      SP=SHFE2
    ELSEIF(isclose(PX,FeO3H3_activity))THEN
      R1=H_1p_activity
      P1=FeO3H3_activity
      P2=OH_1e_activity
      NR1=0
      NP2=0
      SP=SPFE3
    ELSEIF(isclose(PX,FeO4H4_1e_activity))THEN
      R1=OH_1e_activity
      P1=FeO4H4_1e_activity
      P2=H_1p_activity
      NR1=0
      NP2=1
      SP=SHFE4
    ENDIF
    RHFE1=0._r8
    RHFEO1=0._r8
    RHFEO2=0._r8
    RHFEO3=0._r8
    RHFEO4=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1/P2**NP2
    RPFEOX=AMAX1(-Precp_FeO3H3_conc,TPD*(P1-SPX))
    IF(isclose(PX,Fe_3p_activity))THEN
      RHFE1=RPFEOX
    ELSEIF(isclose(PX,FeOH_2p_activity))THEN
      RHFEO1=RPFEOX
    ELSEIF(isclose(PX,FeO2H2_p_activity))THEN
      RHFEO2=RPFEOX
    ELSEIF(isclose(PX,FeO3H3_activity))THEN
      RHFEO3=RPFEOX
    ELSEIF(isclose(PX,FeO4H4_1e_activity))THEN
      RHFEO4=RPFEOX
    ENDIF
!
!     CALCITE
!
    PX=AMAX1(CO3_2e_activity,HCO3_e_activity,H2CO3_activity)
    R1=H_1p_activity
    P1=Ca_2p_activity
    IF(isclose(PX,CO3_2e_activity))THEN
      P2=CO3_2e_activity
      NR1=0
      SP=SPCAC
    ELSEIF(isclose(PX,HCO3_e_activity))THEN
      P2=HCO3_e_activity
      NR1=1
      SP=SHCAC1
    ELSEIF(isclose(PX,H2CO3_activity))THEN
      P2=H2CO3_activity
      NR1=2
      SP=SHCAC2
    ENDIF
    RHCAC3=0._r8
    RHCACH=0._r8
    RHCACO=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1
    S0=P1+P2
    S1=AZMAX1(S0**2_r8-4.0_r8*(P1*P2-SPX))
    RPCACX=AMAX1(-Precp_CaCO3_conc,TPD*(S0-SQRT(S1)))
    IF(isclose(PX,CO3_2e_activity))THEN
      RHCAC3=RPCACX
    ELSEIF(isclose(PX,HCO3_e_activity))THEN
      RHCACH=RPCACX
    ELSEIF(isclose(PX,H2CO3_activity))THEN
      RHCACO=RPCACX
    ENDIF
!
!     GYPSUM
!
    P1=Ca_2p_activity
    P2=SO4_2e_activity
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SPCAS
    S0=P1+P2
    S1=AZMAX1(S0**2-4.0*(P1*P2-SPX))
    RPCASO=AMAX1(-Precp_CaSO4_conc,TPD*(S0-SQRT(S1)))
!
!     PHOSPHORUS PRECIPITATION-DISSOLUTION IN NON-BAND SOIL ZONE
!
!
!     ALUMINUM PHOSPHATE (VARISCITE)
!
    PX=AMAX1(Al_3p_activity,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity)
    PY=AMAX1(H1PO4_2e_activity,H2PO4_1e_activity)
    R1=H_1p_activity
    P3=H_1p_activity
    IF(isclose(PY,H1PO4_2e_activity))THEN
      P2=H1PO4_2e_activity
      IF(isclose(PX,Al_3p_activity))THEN
        P1=Al_3p_activity
        NR1=1
        NP3=0
        SP=SHA0P1
      ELSEIF(isclose(PX,AlOH_2p_activity))THEN
        P1=AlOH_2p_activity
        NR1=0
        NP3=0
        SP=SPA1P1
      ELSEIF(isclose(PX,AlO2H2_1p_activity))THEN
        P1=AlO2H2_1p_activity
        NR1=0
        NP3=1
        SP=SHA2P1
      ELSEIF(isclose(PX,AlO3H3_activity))THEN
        P1=AlO3H3_activity
        NR1=0
        NP3=2
        SP=SHA3P1
      ELSEIF(isclose(PX,AlO4H4_1e_activity))THEN
        P1=AlO4H4_1e_activity
        NR1=0
        NP3=3
        SP=SHA4P1
      ENDIF
    ELSE
      P2=H2PO4_1e_activity
      IF(isclose(PX,Al_3p_activity))THEN
        P1=Al_3p_activity
        NR1=2
        NP3=0
        SP=SHA0P2
      ELSEIF(isclose(PX,AlOH_2p_activity))THEN
        P1=AlOH_2p_activity
        NR1=1
        NP3=0
        SP=SHA1P2
      ELSEIF(isclose(PX,AlO2H2_1p_activity))THEN
        P1=AlO2H2_1p_activity
        NR1=0
        NP3=0
        SP=SPA2P2
      ELSEIF(isclose(PX,AlO3H3_activity))THEN
        P1=AlO3H3_activity
        NR1=0
        NP3=1
        SP=SHA3P2
      ELSEIF(isclose(PX,AlO4H4_1e_activity))THEN
        P1=AlO4H4_1e_activity
        NR1=0
        NP3=2
        SP=SHA4P2
      ENDIF
    ENDIF
    RHA0P1=0._r8
    RHA1P1=0._r8
    RHA2P1=0._r8
    RHA3P1=0._r8
    RHA4P1=0._r8
    RHA0P2=0._r8
    RHA1P2=0._r8
    RHA2P2=0._r8
    RHA3P2=0._r8
    RHA4P2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    P3=AMAX1(ZERO,P3)
    SPX=SP*R1**NR1/P3**NP3
    S0=P1+P2
    S1=AZMAX1(S0**2-4.0*(P1*P2-SPX))
    H2PO4_1e_AlPO4_dissol_flx=AMAX1(-Precp_AlPO4_conc,TPD*(S0-SQRT(S1)))
    IF(isclose(PY,H1PO4_2e_activity))THEN
      IF(isclose(PX,Al_3p_activity))THEN
        RHA0P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AlOH_2p_activity))THEN
        RHA1P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AlO2H2_1p_activity))THEN
        RHA2P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AlO3H3_activity))THEN
        RHA3P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AlO4H4_1e_activity))THEN
        RHA4P1=H2PO4_1e_AlPO4_dissol_flx
      ENDIF
    ELSE
      IF(isclose(PX,Al_3p_activity))THEN
        RHA0P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AlOH_2p_activity))THEN
        RHA1P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AlO2H2_1p_activity))THEN
        RHA2P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AlO3H3_activity))THEN
        RHA3P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AlO4H4_1e_activity))THEN
        RHA4P2=H2PO4_1e_AlPO4_dissol_flx
      ENDIF
    ENDIF
!
!     IRON PHOSPHATE (STRENGITE)
!
    PX=AMAX1(Fe_3p_activity,FeOH_2p_activity,FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity)
    PY=AMAX1(H1PO4_2e_activity,H2PO4_1e_activity)
    R1=H_1p_activity
    P3=H_1p_activity
    IF(isclose(PY,H1PO4_2e_activity))THEN
      P2=H1PO4_2e_activity
      IF(isclose(PX,Fe_3p_activity))THEN
        P1=Fe_3p_activity
        NR1=1
        NP3=0
        SP=SHF0P1
      ELSEIF(isclose(PX,FeOH_2p_activity))THEN
        P1=FeOH_2p_activity
        NR1=0
        NP3=0
        SP=SPF1P1
      ELSEIF(isclose(PX,FeO2H2_p_activity))THEN
        P1=FeO2H2_p_activity
        NR1=0
        NP3=1
        SP=SHF2P1
      ELSEIF(isclose(PX,FeO3H3_activity))THEN
        P1=FeO3H3_activity
        NR1=0
        NP3=2
        SP=SHF3P1
      ELSEIF(isclose(PX,FeO4H4_1e_activity))THEN
        P1=FeO4H4_1e_activity
        NR1=0
        NP3=3
        SP=SHF4P1
      ENDIF
    ELSE
      P2=H2PO4_1e_activity
      IF(isclose(PX,Fe_3p_activity))THEN
        P1=Fe_3p_activity
        NR1=2
        NP3=0
        SP=SHF0P2
      ELSEIF(isclose(PX,FeOH_2p_activity))THEN
        P1=FeOH_2p_activity
        NR1=1
        NP3=0
        SP=SHF1P2
      ELSEIF(isclose(PX,FeO2H2_p_activity))THEN
        P1=FeO2H2_p_activity
        NR1=0
        NP3=0
        SP=SPF2P2
      ELSEIF(isclose(PX,FeO3H3_activity))THEN
        P1=FeO3H3_activity
        NR1=0
        NP3=1
        SP=SHF3P2
      ELSEIF(isclose(PX,FeO4H4_1e_activity))THEN
        P1=FeO4H4_1e_activity
        NR1=0
        NP3=2
        SP=SHF4P2
      ENDIF
    ENDIF
    RHF0P1=0._r8
    RHF1P1=0._r8
    RHF2P1=0._r8
    RHF3P1=0._r8
    RHF4P1=0._r8
    RHF0P2=0._r8
    RHF1P2=0._r8
    RHF2P2=0._r8
    RHF3P2=0._r8
    RHF4P2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    P3=AMAX1(ZERO,P3)
    SPX=SP*R1**NR1/P3**NP3
    S0=P1+P2
    S1=AZMAX1(S0**2_r8-4.0_r8*(P1*P2-SPX))
    H2PO4_1e_FePO4_dissol_flx=AMAX1(-Precp_FePO4_conc,TPD*(S0-SQRT(S1)))
    IF(isclose(PY,H1PO4_2e_activity))THEN
      IF(isclose(PX,Fe_3p_activity))THEN
        RHF0P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,FeOH_2p_activity))THEN
        RHF1P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,FeO2H2_p_activity))THEN
        RHF2P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,FeO3H3_activity))THEN
        RHF3P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,FeO4H4_1e_activity))THEN
        RHF4P1=H2PO4_1e_FePO4_dissol_flx
      ENDIF
    ELSE
      IF(isclose(PX,Fe_3p_activity))THEN
        RHF0P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,FeOH_2p_activity))THEN
        RHF1P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,FeO2H2_p_activity))THEN
        RHF2P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,FeO3H3_activity))THEN
        RHF3P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,FeO4H4_1e_activity))THEN
        RHF4P2=H2PO4_1e_FePO4_dissol_flx
      ENDIF
    ENDIF
!
!     DICALCIUM PHOSPHATE
!
    PX=AMAX1(H1PO4_2e_activity,H2PO4_1e_activity)
    R1=H_1p_activity
    P1=Ca_2p_activity
    IF(isclose(PX,H1PO4_2e_activity))THEN
      P2=H1PO4_2e_activity
      NR1=0
      SP=SPCAD
    ELSEIF(isclose(PX,H2PO4_1e_activity))THEN
      P2=H2PO4_1e_activity
      NR1=1
      SP=SHCAD2
    ENDIF
    RPCAD1=0._r8
    RHCAD2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=SP*R1**NR1
    S0=P1+P2
    S1=AZMAX1(S0**2_r8-4.0_r8*(P1*P2-SPX))
    H2PO4_1e_CaHPO4_dissol_flx=AMAX1(-Precp_CaHPO4_conc,TPD*(S0-SQRT(S1)))
    IF(isclose(PX,H1PO4_2e_activity))THEN
      RPCAD1=H2PO4_1e_CaHPO4_dissol_flx
    ELSEIF(isclose(PX,H2PO4_1e_activity))THEN
      RHCAD2=H2PO4_1e_CaHPO4_dissol_flx
    ENDIF
!
!     HYDROXYAPATITE
!
    PX=AMAX1(H1PO4_2e_activity,H2PO4_1e_activity)
    R1=H_1p_activity
    P1=Ca_2p_activity
    IF(isclose(PX,H1PO4_2e_activity))THEN
      P2=H1PO4_2e_activity
      NR1=4
      SP=SHCAH1
    ELSEIF(isclose(PX,H2PO4_1e_activity))THEN
      P2=H2PO4_1e_activity
      NR1=7
      SP=SHCAH2
    ENDIF
    RHCAH1=0._r8
    RHCAH2=0._r8
    R1=AMAX1(ZERO,R1)
    P1=AMAX1(ZERO,P1)
    P2=AMAX1(ZERO,P2)
    SPX=(SP*R1**NR1/P1**5_r8)**0.333_r8
    H2PO4_1e_apatite_dissol_flx=AMAX1(-Precp_Ca5P3O12O3H3_conc,TPD*(P2-SPX))
    IF(isclose(PX,H1PO4_2e_activity))THEN
      RHCAH1=H2PO4_1e_apatite_dissol_flx
    ELSEIF(isclose(PX,H2PO4_1e_activity))THEN
      RHCAH2=H2PO4_1e_apatite_dissol_flx
    ENDIF
    Precp_AlO3H3_conc=Precp_AlO3H3_conc+RPALOX
    Precp_FeO3H3_conc=Precp_FeO3H3_conc+RPFEOX
    Precp_CaCO3_conc=Precp_CaCO3_conc+RPCACX
    Precp_CaSO4_conc=Precp_CaSO4_conc+RPCASO
    Precp_AlPO4_conc=Precp_AlPO4_conc+H2PO4_1e_AlPO4_dissol_flx
    Precp_FePO4_conc=Precp_FePO4_conc+H2PO4_1e_FePO4_dissol_flx
    Precp_CaHPO4_conc=Precp_CaHPO4_conc+H2PO4_1e_CaHPO4_dissol_flx
    Precp_Ca5P3O12O3H3_conc=Precp_Ca5P3O12O3H3_conc+H2PO4_1e_apatite_dissol_flx
!
!     ANION EXCHANGE EQILIBRIA
!
    IF(VLWatMicP.GT.ZEROS)THEN
      VLWatMicPBK=AMIN1(1.0_r8,BulkSoilMass/VLWatMicP)
    ELSE
      VLWatMicPBK=1.0_r8
    ENDIF
    IF(XAEC.GT.ZEROS)THEN
      RXOH2=TAD*(XROH1_conc*H_1p_activity-SXOH2*XROH2_conc)/(XROH1_conc+SPOH2)*VLWatMicPBK
      RXOH1=TAD*(XOH_conc*H_1p_activity-SXOH1*XROH1_conc)/(XOH_conc+SPOH1)*VLWatMicPBK
      SPH2P=SXH2P*DPH2O
      H2PO4_1e_to_XH2PO4_ROH2_flx=TAD*(XROH2_conc*H2PO4_1e_activity-SPH2P*XH2PO4_conc)/(XROH2_conc+SPH2P)*VLWatMicPBK
      H2PO4_to_XH2PO4_ROH_flx=TAD*(XROH1_conc*H2PO4_1e_activity-SXH2P*XH2PO4_conc*OH_1e_activity) &
        /(XROH1_conc+SXH2P*OH_1e_activity)*VLWatMicPBK
!
!     HPO4 EXCHANGE
!
      SPH1P=SXH1P*DPH2O/DPH2P
      H1PO4_to_XHPO4_ROH_flx=TAD*(XROH1_conc*H2PO4_1e_activity-SPH1P*XHPO4_conc)/(XROH1_conc+SPH1P)*VLWatMicPBK
      XOH_conc=XOH_conc-RXOH1
      XROH1_conc=XROH1_conc+RXOH1-RXOH2-H2PO4_to_XH2PO4_ROH_flx-H1PO4_to_XHPO4_ROH_flx
      XROH2_conc=XROH2_conc+RXOH2-H2PO4_1e_to_XH2PO4_ROH2_flx
      XHPO4_conc=XHPO4_conc+H1PO4_to_XHPO4_ROH_flx
      XH2PO4_conc=XH2PO4_conc+H2PO4_1e_to_XH2PO4_ROH2_flx+H2PO4_to_XH2PO4_ROH_flx
    ELSE
      RXOH2=0._r8
      RXOH1=0._r8
      H2PO4_1e_to_XH2PO4_ROH2_flx=0._r8
      H2PO4_to_XH2PO4_ROH_flx=0._r8
      H1PO4_to_XHPO4_ROH_flx=0._r8
    ENDIF
!
!     CATION EXCHANGE
!
    IF(XCEC.GT.ZEROS)THEN
      AALX=Al_3p_activity**0.333_r8
      AFEX=Fe_3p_activity**0.333_r8
      ACAX=Ca_2p_activity**0.500_r8
      AMGX=Mg_2p_activity**0.500_r8
      XCAX=CEC_conc/(1.0_r8+GKC4*NH4_1p_activity/ACAX &
       +GKCH*H_1p_activity/ACAX+GKCA*AALX/ACAX &
       +GKCA*AFEX/ACAX+GKCM*AMGX/ACAX &
       +GKCN*Na_1p_activity/ACAX+GKCK*K_1p_activity/ACAX)
      XN4Q=XCAX*NH4_1p_activity*GKC4
      XHYQ=XCAX*H_1p_activity*GKCH
      XALQ=XCAX*AALX*GKCA
      XFEQ=XCAX*AFEX*GKCA
      XCAQ=XCAX*ACAX
      XMGQ=XCAX*AMGX*GKCM
      XNAQ=XCAX*Na_1p_activity*GKCN
      XKAQ=XCAX*K_1p_activity*GKCK
      XTLQ=XN4Q+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ
      IF(XTLQ.GT.ZERO)THEN
        FX=CEC_conc/XTLQ
      ELSE
        FX=0._r8
      ENDIF
      XN4Q=FX*XN4Q
      XHYQ=FX*XHYQ
      XALQ=FX*XALQ/3.0_r8
      XFEQ=FX*XFEQ/3.0_r8
      XCAQ=FX*XCAQ/2.0_r8
      XMGQ=FX*XMGQ/2.0_r8
      XNAQ=FX*XNAQ
      XKAQ=FX*XKAQ
      RXN4=TAD*AMIN1((XN4Q-XNH4_mole_conc)*NH4_1p_activity/XN4Q,NH4_1p_aqua_mole_conc)
      RXHY=TAD*AMIN1((XHYQ-XHY1)*H_1p_activity/XHYQ,H_1p_aqua_mole_conc)
      RXAL=TAD*AMIN1((XALQ-XAl_conc)*AALX/XALQ,Al_3p_aqua_mole_conc)
      RXFE=TAD*AMIN1((XFEQ-XFe_conc)*AFEX/XFEQ,Fe_3p_aqua_mole_conc)
      RXCA=TAD*AMIN1((XCAQ-XCa_conc)*ACAX/XCAQ,Ca_2p_aqua_mole_conc)
      RXMG=TAD*AMIN1((XMGQ-XMg_conc)*AMGX/XMGQ,Mg_2p_aqua_mole_conc)
      RXNA=TAD*AMIN1((XNAQ-XNa_conc)*Na_1p_activity/XNAQ,Na_1p_aqua_mole_conc)
      RXKA=TAD*AMIN1((XKAQ-XK_conc)*K_1p_activity/XKAQ,K_1p_aqua_mole_conc)
      XNH4_mole_conc=XNH4_mole_conc+RXN4
      XHY1=XHY1+RXHY
      XAl_conc=XAl_conc+RXAL
      XFe_conc=XFe_conc+RXFE
      XCa_conc=XCa_conc+RXCA
      XMg_conc=XMg_conc+RXMG
      XNa_conc=XNa_conc+RXNA
      XK_conc=XK_conc+RXKA
    ELSE
      RXN4=0._r8
      RXHY=0._r8
      RXAL=0._r8
      RXFE=0._r8
      RXCA=0._r8
      RXMG=0._r8
      RXNA=0._r8
      RXKA=0._r8
    ENDIF
!
!     ORGANIC MATTER
!
    DP=DPCOH*DPALO
    XHC1=H_1p_activity*(XCOOH_conc-XAlO2H2_conc-XFeO2H2_conc)/(H_1p_activity+DPCOH)
    XAlO2H2_conc=AlO2H2_1p_activity*(XCOOH_conc-XHC1)/(AlO2H2_1p_activity+DPALO)
    XFeO2H2_conc=FeO2H2_p_activity*(XCOOH_conc-XHC1)/(FeO2H2_p_activity+DPFEO)
    XCOO=AZMAX1(XCOOH_conc-XHC1-XAlO2H2_conc-XFeO2H2_conc)
  ELSE
    RHAL1=0._r8
    RHALO1=0._r8
    RHALO2=0._r8
    RHALO3=0._r8
    RHALO4=0._r8
    RHFE1=0._r8
    RHFEO1=0._r8
    RHFEO2=0._r8
    RHFEO3=0._r8
    RHFEO4=0._r8
    RHCAC3=0._r8
    RHCACH=0._r8
    RHCACO=0._r8
    RPCACX=0._r8
    RPCASO=0._r8
    H2PO4_1e_CaHPO4_dissol_flx=0._r8
    H2PO4_1e_apatite_dissol_flx=0._r8
    RHA0P1=0._r8
    RHA1P1=0._r8
    RHA2P1=0._r8
    RHA3P1=0._r8
    RHA4P1=0._r8
    RHA0P2=0._r8
    RHA1P2=0._r8
    RHA2P2=0._r8
    RHA3P2=0._r8
    RHA4P2=0._r8
    RHF0P1=0._r8
    RHF1P1=0._r8
    RHF2P1=0._r8
    RHF3P1=0._r8
    RHF4P1=0._r8
    RHF0P2=0._r8
    RHF1P2=0._r8
    RHF2P2=0._r8
    RHF3P2=0._r8
    RHF4P2=0._r8
    RPCAD1=0._r8
    RHCAD2=0._r8
    RHCAH1=0._r8
    RHCAH2=0._r8
    RXOH2=0._r8
    RXOH1=0._r8
    H2PO4_1e_to_XH2PO4_ROH2_flx=0._r8
    H2PO4_to_XH2PO4_ROH_flx=0._r8
    H1PO4_to_XHPO4_ROH_flx=0._r8
    RXN4=0._r8
    RXHY=0._r8
    RXAL=0._r8
    RXFE=0._r8
    RXCA=0._r8
    RXMG=0._r8
    RXNA=0._r8
    RXKA=0._r8
  ENDIF
!
!     ION SPECIATION
!
  S0=H_1p_activity+NH3_activity+DPN4
  S1=S0**2_r8-4.0_r8*(H_1p_activity*NH3_activity-DPN4*NH4_1p_activity)
  RNH4=TSL*(S0-SQRT(S1))
  S0=Al_3p_activity+OH_1e_activity+DPAL1
  S1=S0**2_r8-4.0_r8*(Al_3p_activity*OH_1e_activity-DPAL1*AlOH_2p_activity)
  RALO1=TSL*(S0-SQRT(S1))
  S0=AlOH_2p_activity+OH_1e_activity+DPAL2
  S1=S0**2_r8-4.0_r8*(AlOH_2p_activity*OH_1e_activity-DPAL2*AlO2H2_1p_activity)
  RALO2=TSL*(S0-SQRT(S1))
  S0=AlO2H2_1p_aqua_mole_conc+OH_1e_aqua_mole_conc+DPAL3
  S1=S0**2-4.0*(AlO2H2_1p_activity*OH_1e_activity-DPAL3*AlO3H3_activity)
  RALO3=TSL*(S0-SQRT(S1))
  S0=AlO3H3_activity+OH_1e_activity+DPAL4
  S1=S0**2-4.0*(AlO3H3_activity*OH_1e_activity-DPAL4*AlO4H4_1e_activity)
  RALO4=TSL*(S0-SQRT(S1))
  S0=Al_3p_activity+SO4_2e_activity+DPALS
  S1=S0**2-4.0*(Al_3p_activity*SO4_2e_activity-DPALS*AlSO4_1p_activity)
  RALS=TSL*(S0-SQRT(S1))
  S0=Fe_3p_activity+OH_1e_activity+DPFE1
  S1=S0**2-4.0*(Fe_3p_activity*OH_1e_activity-DPFE1*FeOH_2p_activity)
  RFEO1=TSL*(S0-SQRT(S1))
  S0=FeOH_2p_activity+OH_1e_activity+DPFE2
  S1=S0**2-4.0*(FeOH_2p_activity*OH_1e_activity-DPFE2*FeO2H2_p_activity)
  RFEO2=TSL*(S0-SQRT(S1))
  S0=FeO2H2_p_activity+OH_1e_activity+DPFE3
  S1=S0**2-4.0*(FeO2H2_p_activity*OH_1e_activity-DPFE3*FeO3H3_activity)
  RFEO3=TSL*(S0-SQRT(S1))
  S0=FeO3H3_activity+OH_1e_activity+DPFE4
  S1=S0**2-4.0*(FeO3H3_activity*OH_1e_activity-DPFE4*FeO4H4_1e_activity)
  RFEO4=TSL*(S0-SQRT(S1))
  S0=Fe_3p_activity+SO4_2e_activity+DPFES
  S1=S0**2-4.0*(Fe_3p_activity*SO4_2e_activity-DPFES*FeSO4_1p_activity)
  RFES=TSL*(S0-SQRT(S1))
  S0=Ca_2p_activity+OH_1e_activity+DPCAO
  S1=S0**2-4.0*(Ca_2p_activity*OH_1e_activity-DPCAO*CaO2H2_activity)
  RCAO=TSL*(S0-SQRT(S1))
  S0=Ca_2p_activity+CO3_2e_activity+DPCAC
  S1=S0**2-4.0*(Ca_2p_activity*CO3_2e_activity-DPCAC*CaCO3_activity)
  RCAC=TSL*(S0-SQRT(S1))
  S0=Ca_2p_activity+HCO3_e_activity+DPCAH
  S1=S0**2-4.0*(Ca_2p_activity*HCO3_e_activity-DPCAH*CaHCO3_1p_activity)
  RCAH=TSL*(S0-SQRT(S1))
  S0=Ca_2p_activity+SO4_2e_activity+DPCAS
  S1=S0**2-4.0*(Ca_2p_activity*SO4_2e_activity-DPCAS*CaSO4_activity)
  RCAS=TSL*(S0-SQRT(S1))
  S0=Mg_2p_activity+OH_1e_activity+DPMGO
  S1=S0**2-4.0*(Mg_2p_activity*OH_1e_activity-DPMGO*MgOH_1p_activity)
  RMGO=TSL*(S0-SQRT(S1))
  S0=Mg_2p_activity+CO3_2e_activity+DPMGC
  S1=S0**2-4.0*(Mg_2p_activity*CO3_2e_activity-DPMGC*MgCO3_activity)
  RMGC=TSL*(S0-SQRT(S1))
  S0=Mg_2p_activity+HCO3_e_activity+DPMGH
  S1=S0**2-4.0*(Mg_2p_activity*HCO3_e_activity-DPMGH*MgHCO3_1p_activity)
  RMGH=TSL*(S0-SQRT(S1))
  S0=Mg_2p_activity+SO4_2e_activity+DPMGS
  S1=S0**2-4.0*(Mg_2p_activity*SO4_2e_activity-DPMGS*MgSO4_activity)
  RMGS=TSL*(S0-SQRT(S1))
  S0=Na_1p_activity+CO3_2e_activity+DPNAC
  S1=S0**2-4.0*(Na_1p_activity*CO3_2e_activity-DPNAC*NaCO3_1e_activity)
  RNAC=TSL*(S0-SQRT(S1))
  S0=Na_1p_activity+SO4_2e_activity+DPNAS
  S1=S0**2-4.0*(Na_1p_activity*SO4_2e_activity-DPNAS*NaSO4_1e_activity)
  RNAS=TSL*(S0-SQRT(S1))
  S0=K_1p_activity+SO4_2e_activity+DPKAS
  S1=S0**2-4.0*(K_1p_activity*SO4_2e_activity-DPKAS*KSO4_1e_activity)
  RKAS=TSL*(S0-SQRT(S1))
  S0=H0PO4_3e_activity+H_1p_activity+DPH1P
  S1=S0**2-4.0*(H0PO4_3e_activity*H_1p_activity-DPH1P*H1PO4_2e_activity)
  RH1P=TSL*(S0-SQRT(S1))
  S0=H1PO4_2e_activity+H_1p_activity+DPH2P
  S1=S0**2-4.0*(H1PO4_2e_activity*H_1p_activity-DPH2P*H2PO4_1e_activity)
  H2PO4_e_to_HPO4_2e_flx=TSL*(S0-SQRT(S1))
  S0=H2PO4_1e_activity+H_1p_activity+DPH3P
  S1=S0**2-4.0*(H2PO4_1e_activity*H_1p_activity-DPH3P*H3PO4_activity)
  RH3P=TSL*(S0-SQRT(S1))
  S0=Fe_3p_activity+H1PO4_2e_activity+DPF1P
  S1=S0**2-4.0*(Fe_3p_activity*H1PO4_2e_activity-DPF1P*FeHPO4_p_activity)
  RF1P=TSL*(S0-SQRT(S1))
  S0=Fe_3p_activity+H2PO4_1e_activity+DPF2P
  S1=S0**2-4.0*(Fe_3p_activity*H2PO4_1e_activity-DPF2P*FeH2PO4_2p_activity)
  RF2P=TSL*(S0-SQRT(S1))
  S0=Ca_2p_activity+H0PO4_3e_activity+DPC0P
  S1=S0**2-4.0*(Ca_2p_activity*H0PO4_3e_activity-DPC0P*CaPO4_1e_activity)
  RC0P=TSL*(S0-SQRT(S1))
  S0=Ca_2p_activity+H1PO4_2e_activity+DPC1P
  S1=S0**2-4.0*(Ca_2p_activity*H1PO4_2e_activity-DPC1P*CaHPO4_activity)
  RC1P=TSL*(S0-SQRT(S1))
  S0=Ca_2p_activity+H2PO4_1e_activity+DPC2P
  S1=S0**2-4.0*(Ca_2p_activity*H2PO4_1e_activity-DPC2P*CaH4P2O8_1p_activity)
  RC2P=TSL*(S0-SQRT(S1))
  S0=Mg_2p_activity+H1PO4_2e_activity+DPM1P
  S1=S0**2-4.0*(Mg_2p_activity*H1PO4_2e_activity-DPM1P*MgHPO4_activity)
  RM1P=TSL*(S0-SQRT(S1))
!
!     TOTAL ION FLUX FOR EACH ION SPECIES
!
  RN4S=RNH4-RXN4
  RN3S=-RNH4
  RAL=-RHAL1-RHA0P1-RHA0P2-RALO1-RALS-RXAL
  RFE=-RHFE1-RHF0P1-RHF0P2-RFEO1-RFES-RXFE-RF1P-RF2P
  RCA=-RPCACX-RPCASO-H2PO4_1e_CaHPO4_dissol_flx-5.0*H2PO4_1e_apatite_dissol_flx-RXCA-RCAO-RCAC-RCAH-RCAS-RC0P-RC1P-RC2P
  RMG=-RMGO-RMGC-RMGH-RMGS-RM1P-RXMG
  RNA=-RNAC-RNAS-RXNA
  RKA=-RKAS-RXKA
  RSO4=-RPCASO-RALS-RFES-RCAS-RMGS-RNAS-RKAS
  RAL1=-RHALO1-RHA1P1-RHA1P2+RALO1-RALO2
  RAL2=-RHALO2-RHA2P1-RHA2P2+RALO2-RALO3
  RAL3=-RHALO3-RHA3P1-RHA3P2+RALO3-RALO4
  RAL4=-RHALO4-RHA4P1-RHA4P2+RALO4
  RFE1=-RHFEO1-RHF1P1-RHF1P2+RFEO1-RFEO2
  RFE2=-RHFEO2-RHF2P1-RHF2P2+RFEO2-RFEO3
  RFE3=-RHFEO3-RHF3P1-RHF3P2+RFEO3-RFEO4
  RFE4=-RHFEO4-RHF4P1-RHF4P2+RFEO4
  RHP0=-RH1P-RC0P
  RHP1=-RHA0P1-RHA1P1-RHA2P1-RHA3P1 &
    -RHA4P1-RHF0P1-RHF1P1-RHF2P1 &
    -RHF3P1-RHF4P1-RPCAD1-3.0*RHCAH1 &
    -H1PO4_to_XHPO4_ROH_flx+RH1P-H2PO4_e_to_HPO4_2e_flx-RF1P-RC1P-RM1P
  RHP2=-RHA0P2-RHA1P2-RHA2P2-RHA3P2 &
    -RHA4P2-RHF0P2-RHF1P2-RHF2P2 &
    -RHF3P2-RHF4P2-RHCAD2-3.0*RHCAH2 &
    -H2PO4_1e_to_XH2PO4_ROH2_flx-H2PO4_to_XH2PO4_ROH_flx+H2PO4_e_to_HPO4_2e_flx-RH3P-RF2P-RC2P
  RHP3=RH3P
!
!     ION CONCENTRATIONS
!
  CCO2X=CCO2M*SCO2X/(EXP(GasSechenovConst(idg_CO2)*CSTRZ))*EXP(0.843-0.0281*ATCA)*FH2O
  CCO2Y=LOG(CCO2X)
  CCO2Z=ABS(CCO2Y)
  H2CO3_aqua_mole_conc=CCO2X
  FCO3=DPCO3*A0/(H_1p_activity**2*A2)
  FHCO=DPCO2*A0/(H_1p_activity*A1)
  Z=GasSechenovConst(idg_CO2)*(2.0E-03*FCO3+0.5E-03*FHCO)
  DO  MM=1,25
    R=(LOG(H2CO3_aqua_mole_conc)+Z*H2CO3_aqua_mole_conc-CCO2Y)/CCO2Z
    IF(R.LT.1.0E-03)exit
    H2CO3_aqua_mole_conc=H2CO3_aqua_mole_conc/SQRT(1.0_r8+R)
  ENDDO
  CH4_aqua_mole_conc    = CCH4M*SCH4X/(EXP(GasSechenovConst(idg_CH4)*CSTR1))*EXP(0.597-0.0199*ATCA)*FH2O
  O2_aqua_mole_conc     = COXYM*SOXYX/(EXP(GasSechenovConst(idg_O2)*CSTR1))*EXP(0.516-0.0172*ATCA)*FH2O
  N2_aqua_mole_conc     = CZ2GM*SN2GX/(EXP(GasSechenovConst(idg_N2)*CSTR1))*EXP(0.456-0.0152*ATCA)*FH2O
  N2O_aqua_mole_conc    = CZ2OM*SN2OX/(EXP(GasSechenovConst(idg_N2O)*CSTR1))*EXP(0.897-0.0299*ATCA)*FH2O
  NH4_1p_aqua_mole_conc = NH4_1p_aqua_mole_conc+RN4S
  NH3_aqua_mole_conc    = NH3_aqua_mole_conc+RN3S
  Al_3p_aqua_mole_conc  = Al_3p_aqua_mole_conc+RAL
  Fe_3p_aqua_mole_conc  = Fe_3p_aqua_mole_conc+RFE
  Ca_2p_aqua_mole_conc  = Ca_2p_aqua_mole_conc+RCA
  Mg_2p_aqua_mole_conc  = Mg_2p_aqua_mole_conc+RMG
  Na_1p_aqua_mole_conc  = Na_1p_aqua_mole_conc+RNA
  K_1p_aqua_mole_conc   = K_1p_aqua_mole_conc+RKA
  SO4_2e_aqua_mole_conc    = SO4_2e_aqua_mole_conc+RSO4
  CO3_2e_aqua_mole_conc    = H2CO3_aqua_mole_conc*DPCO3*A0/(H_1p_activity**2*A2)
  HCO3_e_conc              = H2CO3_aqua_mole_conc*DPCO2*A0/(H_1p_activity*A1)
  AlOH_2p_aqua_mole_conc   = AlOH_2p_aqua_mole_conc+RAL1
  AlO2H2_1p_aqua_mole_conc = AlO2H2_1p_aqua_mole_conc+RAL2
  AlO3H3_conc              = AlO3H3_conc+RAL3
  AlO4H4_1e_aqua_mole_conc = AlO4H4_1e_aqua_mole_conc+RAL4
  AlSO4_1p_aqua_mole_conc  = AlSO4_1p_aqua_mole_conc+RALS
  FeOH_2p_aqua_mole_conc   = FeOH_2p_aqua_mole_conc+RFE1
  FeO2H2_p_conc            = FeO2H2_p_conc+RFE2
  FeO3H3_conc              = FeO3H3_conc+RFE3
  FeO4H4_1e_aqua_mole_conc = FeO4H4_1e_aqua_mole_conc+RFE4
  FeSO4_1p_aqua_mole_conc  = FeSO4_1p_aqua_mole_conc+RFES
  CaO2H2_conc              = CaO2H2_conc+RCAO
  CaCO3_conc               = CaCO3_conc+RCAC
  CaHCO3_1p_aqua_mole_conc = CaHCO3_1p_aqua_mole_conc+RCAH
  CaSO4_conc               = CaSO4_conc+RCAS
  MgOH_1p_aqua_mole_conc   = MgOH_1p_aqua_mole_conc+RMGO
  MgCO3_conc               = MgCO3_conc+RMGC
  MgHCO3_1p_aqua_mole_conc   = MgHCO3_1p_aqua_mole_conc+RMGH
  MgSO4_conc                 = MgSO4_conc+RMGS
  NaCO3_1e_aqua_mole_conc    = NaCO3_1e_aqua_mole_conc+RNAC
  NaSO4_1e_aqua_mole_conc    = NaSO4_1e_aqua_mole_conc+RNAS
  KSO4_1e_aqua_mole_conc     = KSO4_1e_aqua_mole_conc+RKAS
  H0PO4_3e_conc              = H0PO4_3e_conc+RHP0
  H1PO4_2e_aqua_mole_conc    = H1PO4_2e_aqua_mole_conc+RHP1
  H2PO4_1e_aqua_mole_conc    = H2PO4_1e_aqua_mole_conc+RHP2
  H3PO4_conc                 = H3PO4_conc+RHP3
  FeHPO4_p_conc              = FeHPO4_p_conc+RF1P
  FeH2PO4_2p_aqua_mole_conc  = FeH2PO4_2p_aqua_mole_conc+RF2P
  CaPO4_1e_con               = CaPO4_1e_con+RC0P
  CaHPO4_conc                = CaHPO4_conc+RC1P
  CaH4P2O8_1p_aqua_mole_conc = CaH4P2O8_1p_aqua_mole_conc+RC2P
  MgHPO4_conc                = MgHPO4_conc+RM1P
!  IF(K.EQ.micpar%k_POM.AND.(M/1)*1.EQ.M)THEN
!     WRITE(*,1112)'A1I',I,NX,NY,L,K,M,A1,A2,A3,FSTR2,CSTR1
!    2,CSTR2,cation_3p_aqua_mole_conc,anion_3e_conc,cation_2p_aqua_mole_conc,anion_2e_aqua_mole_conc,cation_1p_aqua_mole_conc,CA1,VLWatMicP
!     WRITE(*,1112)'ALPO4I',I,NX,NY,L,K,M,Precp_AlPO4_conc,Al_3p_activity
!    2,AlOH_2p_activity,AlO2H2_1p_activity,AlO3H3_activity,AlO4H4_1e_activity
!    2,H0PO4_3e_activity,H1PO4_2e_activity,H2PO4_1e_activity,H_1p_activity,OH_1e_activity,H2PO4_1e_AlPO4_dissol_flx,RHA0P1,RHA1P1,RHA2P1,RHA3P1
!    3,RHA4P1,RHA0P2,RHA1P2,RHA2P2,RHA3P2,RHA4P2,SP,SPX,Al_3p_activity*H0PO4_3e_activity
!    4,SPALP,H0PO4_3e_conc,H1PO4_2e_aqua_mole_conc,H2PO4_1e_aqua_mole_conc,H3PO4_conc,RHP0,RHP1,RHP2,RHP3
!    5,RAL,RHAL1,RHA0P1,RHA0P2,RALO1,RALS,RXAL
!     WRITE(*,1112)'FEPO4I',I,NX,NY,L,K,M,Precp_FePO4_conc,Fe_3p_activity
!    2,FeOH_2p_activity,FeO2H2_p_activity,FeO3H3_activity,FeO4H4_1e_activity
!    2,H0PO4_3e_activity,H1PO4_2e_activity,H2PO4_1e_activity,H_1p_activity,OH_1e_activity,H2PO4_1e_FePO4_dissol_flx,RHF0P1,RHF1P1,RHF2P1,RHF3P1
!    3,RHF4P1,RHF0P2,RHF1P2,RHF2P2,RHF3P2,RHF4P2,SP,SPX,Fe_3p_activity*H0PO4_3e_activity
!    4,SPFEP
!     WRITE(*,1112)'APATITEI',I,NX,NY,L,K,M,Precp_Ca5P3O12O3H3_conc,Ca_2p_activity,XCa_conc
!    2,H0PO4_3e_activity,H1PO4_2e_activity,H2PO4_1e_activity,H_1p_activity,OH_1e_activity,H2PO4_1e_apatite_dissol_flx,RHCAH1,RHCAH2
!    3,SP,SPX,Ca_2p_activity**5*H0PO4_3e_activity**3*OH_1e_activity,SPCAH,SHCAH1,SHCAH2
!    3,H0PO4_3e_conc,H1PO4_2e_aqua_mole_conc,H2PO4_1e_aqua_mole_conc,XOH_conc,XROH1_conc,XROH2_conc,XHPO4_conc,XH2PO4_conc
!    4,RHA0P1,RHA1P1,RHA2P1,RHA3P1
!    2,RHA4P1,RHF0P1,RHF1P1,RHF2P1
!    3,RHF3P1,RHF4P1,RPCAD1,3.0*RHCAH1
!    4,H1PO4_to_XHPO4_ROH_flx,RH1P,H2PO4_e_to_HPO4_2e_flx,RF1P,RC1P,RM1P
!    5,RHA0P2,RHA1P2,RHA2P2,RHA3P2
!    2,RHA4P2,RHF0P2,RHF1P2,RHF2P2
!    3,RHF3P2,RHF4P2,RHCAD2,3.0*RHCAH2
!    4,H2PO4_1e_to_XH2PO4_ROH2_flx,H2PO4_to_XH2PO4_ROH_flx,H2PO4_e_to_HPO4_2e_flx,RH3P,RF2P,RC2P,RH3P
!      ENDIF
  end subroutine SolubilityEquilibiriaSalt
!------------------------------------------------------------------------------------------

  subroutine SolubilityEquilibriaNoSalt(K,M,BulkSoilMass)

  implicit none
  integer, intent(in) :: K, M
  real(r8), intent(in) :: BulkSoilMass
  real(r8) :: H2PO4_1e_AlPO4_eqv,H2PO4_1e_CaHPO4_eqv,XNAQ,XCAQ,SPH1P
  real(r8) :: H2PO4_1e_FePO4_eqv,H2PO4_1e_apatite_eqv,FX,RN4S,RN3S,RNH4
  real(r8) :: H2PO4_1e_AlPO4_dissol_flx,H2PO4_1e_CaHPO4_dissol_flx,H2PO4_e_to_HPO4_2e_flx,RHP1,RHP2
  real(r8) :: H2PO4_1e_apatite_dissol_flx,H2PO4_1e_FePO4_dissol_flx,H1PO4_to_XHPO4_ROH_flx,H2PO4_1e_to_XH2PO4_ROH2_flx
  real(r8) :: RXN4,H2PO4_to_XH2PO4_ROH_flx,SPH2P,VLWatMicPBK,S0,S1
  real(r8) :: XALQ,XCAX,XFEQ,XHYQ,XKAQ,XMGQ,XN4Q
  real(r8) :: XTLQ
! begin_execution
  H2CO3_aqua_mole_conc=AMAX1(ZERO,H2CO3_aqua_mole_conc)
  CO3_2e_aqua_mole_conc=H2CO3_aqua_mole_conc*DPCO3*A0/(H_1p_aqua_mole_conc**2*A2)
  HCO3_e_conc=H2CO3_aqua_mole_conc*DPCO2*A0/(H_1p_aqua_mole_conc*A1)
  NH4_1p_aqua_mole_conc=AMAX1(ZERO,NH4_1p_aqua_mole_conc)
  NH3_aqua_mole_conc=AMAX1(ZERO,NH3_aqua_mole_conc)
  Al_3p_aqua_mole_conc=AMAX1(ZERO,Al_3p_aqua_mole_conc)
  Fe_3p_aqua_mole_conc=AMAX1(ZERO,Fe_3p_aqua_mole_conc)
  Ca_2p_aqua_mole_conc=AMAX1(ZERO,Ca_2p_aqua_mole_conc)
  Ca_2p_aqua_mole_conc=AMIN1(Ca_2p_aqua_mole_conc,SPCAC/(CO3_2e_aqua_mole_conc*A2**2))
  Mg_2p_aqua_mole_conc=AMAX1(ZERO,Mg_2p_aqua_mole_conc)
  Na_1p_aqua_mole_conc=AMAX1(ZERO,Na_1p_aqua_mole_conc)
  K_1p_aqua_mole_conc=AMAX1(ZERO,K_1p_aqua_mole_conc)
  H1PO4_2e_aqua_mole_conc=AMAX1(ZERO,H1PO4_2e_aqua_mole_conc)
  H2PO4_1e_aqua_mole_conc=AMAX1(ZERO,H2PO4_1e_aqua_mole_conc)
!
!     PRECIPITATION-DISSOLUTION FLUXES
!
  IF(K.EQ.micpar%k_POM)THEN
    H2PO4_1e_AlPO4_eqv=SYA0P2/(Al_3p_aqua_mole_conc*OH_1e_aqua_mole_conc**2)
    H2PO4_1e_AlPO4_dissol_flx=AMAX1(-Precp_AlPO4_conc,TPD*(H2PO4_1e_aqua_mole_conc-H2PO4_1e_AlPO4_eqv))
    H2PO4_1e_FePO4_eqv=SYF0P2/(Fe_3p_aqua_mole_conc*OH_1e_aqua_mole_conc**2._r8)
    H2PO4_1e_FePO4_dissol_flx=AMAX1(-Precp_FePO4_conc,TPD*(H2PO4_1e_aqua_mole_conc-H2PO4_1e_FePO4_eqv))
    H2PO4_1e_CaHPO4_eqv=SYCAD2/(Ca_2p_aqua_mole_conc*OH_1e_aqua_mole_conc)
    H2PO4_1e_CaHPO4_dissol_flx=AMAX1(-Precp_CaHPO4_conc,TPD*(H2PO4_1e_aqua_mole_conc-H2PO4_1e_CaHPO4_eqv))
    H2PO4_1e_apatite_eqv=(SYCAH2/(Ca_2p_aqua_mole_conc**5*OH_1e_aqua_mole_conc**7))**0.333
    H2PO4_1e_apatite_dissol_flx=AMAX1(-Precp_Ca5P3O12O3H3_conc,TPD*(H2PO4_1e_aqua_mole_conc-H2PO4_1e_apatite_eqv))
    Precp_AlPO4_conc=Precp_AlPO4_conc+H2PO4_1e_AlPO4_dissol_flx
    Precp_FePO4_conc=Precp_FePO4_conc+H2PO4_1e_FePO4_dissol_flx
    Precp_CaHPO4_conc=Precp_CaHPO4_conc+H2PO4_1e_CaHPO4_dissol_flx
    Precp_Ca5P3O12O3H3_conc=Precp_Ca5P3O12O3H3_conc+H2PO4_1e_apatite_dissol_flx
!
!     ANION EXCHANGE FLUXES
!
    IF(VLWatMicP.GT.ZEROS)THEN
      VLWatMicPBK=AMIN1(1.0,BulkSoilMass/VLWatMicP)
    ELSE
      VLWatMicPBK=1.0
    ENDIF
    IF(XAEC.GT.ZEROS)THEN
      SPH2P=SXH2P*DPH2O
      H2PO4_1e_to_XH2PO4_ROH2_flx=TAD*(XROH2_conc*H2PO4_1e_aqua_mole_conc-SPH2P*XH2PO4_conc)/(XROH2_conc+SPH2P)*VLWatMicPBK
      H2PO4_to_XH2PO4_ROH_flx=TAD*(XROH1_conc*H2PO4_1e_aqua_mole_conc-SXH2P*OH_1e_aqua_mole_conc*XH2PO4_conc) &
        /(XROH1_conc+SXH2P*OH_1e_aqua_mole_conc)*VLWatMicPBK
      SPH1P=SXH1P*DPH2O/DPH2P
      H1PO4_to_XHPO4_ROH_flx=TAD*(XROH1_conc*H1PO4_2e_aqua_mole_conc-SPH1P*XHPO4_conc)/(XROH1_conc+SPH1P)*VLWatMicPBK
      XROH1_conc=XROH1_conc-H2PO4_to_XH2PO4_ROH_flx-H1PO4_to_XHPO4_ROH_flx
      XROH2_conc=XROH2_conc-H2PO4_1e_to_XH2PO4_ROH2_flx
      XHPO4_conc=XHPO4_conc+H1PO4_to_XHPO4_ROH_flx
      XH2PO4_conc=XH2PO4_conc+H2PO4_1e_to_XH2PO4_ROH2_flx+H2PO4_to_XH2PO4_ROH_flx
    ELSE
      H2PO4_1e_to_XH2PO4_ROH2_flx=0._r8
      H2PO4_to_XH2PO4_ROH_flx=0._r8
      H1PO4_to_XHPO4_ROH_flx=0._r8
    ENDIF
!
!     CATION EXCHANGE
!
    IF(XCEC.GT.ZEROS)THEN
      CALX=Al_3p_aqua_mole_conc**0.333
      CFEX=Fe_3p_aqua_mole_conc**0.333
      CaX_conc=Ca_2p_aqua_mole_conc**0.500_r8
      MgX_conc=Mg_2p_aqua_mole_conc**0.500_r8
      XCAX=CEC_conc/(1.0_r8+GKC4*NH4_1p_aqua_mole_conc/CaX_conc &
        +GKCH*H_1p_aqua_mole_conc/CaX_conc+GKCA*CALX/CaX_conc &
        +GKCA*CFEX/CaX_conc+GKCM*MgX_conc/CaX_conc &
        +GKCN*Na_1p_aqua_mole_conc/CaX_conc+GKCK*K_1p_aqua_mole_conc/CaX_conc)
      XN4Q=XCAX*NH4_1p_aqua_mole_conc*GKC4
      XHYQ=XCAX*H_1p_aqua_mole_conc*GKCH
      XALQ=XCAX*CALX*GKCA
      XFEQ=XCAX*CFEX*GKCA
      XCAQ=XCAX*CaX_conc
      XMGQ=XCAX*MgX_conc*GKCM
      XNAQ=XCAX*Na_1p_aqua_mole_conc*GKCN
      XKAQ=XCAX*K_1p_aqua_mole_conc*GKCK
      XTLQ=XN4Q+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ
      IF(XTLQ.GT.ZERO)THEN
        FX=CEC_conc/XTLQ
      ELSE
        FX=0._r8
      ENDIF
      XN4Q=FX*XN4Q
      RXN4=TSL*AMIN1((XN4Q-XNH4_mole_conc)*NH4_1p_aqua_mole_conc/XN4Q,NH4_1p_aqua_mole_conc)
      XNH4_mole_conc=XNH4_mole_conc+RXN4
    ELSE
      RXN4=0._r8
    ENDIF
  ELSE
    H2PO4_1e_AlPO4_dissol_flx=0._r8
    H2PO4_1e_FePO4_dissol_flx=0._r8
    H2PO4_1e_CaHPO4_dissol_flx=0._r8
    H2PO4_1e_apatite_dissol_flx=0._r8
    H2PO4_1e_to_XH2PO4_ROH2_flx=0._r8
    H2PO4_to_XH2PO4_ROH_flx=0._r8
    H1PO4_to_XHPO4_ROH_flx=0._r8
    RXN4=0._r8
  ENDIF
!
!     NH4 <-> NH3
!
  S0=H_1p_aqua_mole_conc+NH3_aqua_mole_conc+DPN4
  S1=AZMAX1(S0**2-4.0*(H_1p_aqua_mole_conc*NH3_aqua_mole_conc-DPN4*NH4_1p_aqua_mole_conc))
  RNH4=TSL*(S0-SQRT(S1))
!
!     H2PO4 <-> HPO4
!
  S0=H1PO4_2e_aqua_mole_conc+H_1p_aqua_mole_conc+DPH2P
  S1=AZMAX1(S0**2-4.0*(H1PO4_2e_aqua_mole_conc*H_1p_aqua_mole_conc-DPH2P*H2PO4_1e_aqua_mole_conc))
  H2PO4_e_to_HPO4_2e_flx=TSL*(S0-SQRT(S1))
!
!     ION FLUXES
!
  RN4S=RNH4-RXN4
  RN3S=-RNH4
  RHP1=-H1PO4_to_XHPO4_ROH_flx-H2PO4_e_to_HPO4_2e_flx
  RHP2=-H2PO4_1e_to_XH2PO4_ROH2_flx-H2PO4_to_XH2PO4_ROH_flx+H2PO4_e_to_HPO4_2e_flx &
    -H2PO4_1e_AlPO4_dissol_flx-H2PO4_1e_FePO4_dissol_flx-H2PO4_1e_CaHPO4_dissol_flx-3.0_r8*H2PO4_1e_apatite_dissol_flx

  NH4_1p_aqua_mole_conc=NH4_1p_aqua_mole_conc+RN4S
  NH3_aqua_mole_conc=NH3_aqua_mole_conc+RN3S
  H1PO4_2e_aqua_mole_conc=H1PO4_2e_aqua_mole_conc+RHP1
  H2PO4_1e_aqua_mole_conc=H2PO4_1e_aqua_mole_conc+RHP2
  end subroutine SolubilityEquilibriaNoSalt

end module InitSoluteMod
