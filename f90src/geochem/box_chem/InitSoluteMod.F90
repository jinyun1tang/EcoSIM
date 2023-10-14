module InitSoluteMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : isclose,AZMAX1
  use SoluteChemDataType, only : solutedtype
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

  real(r8) :: A1,A2,CEC_conc,XCOOH,SPOH2,SPOH1,NO3_1e_conc
  real(r8), pointer :: H2CO3_aqu_conc
  real(r8), pointer :: CH4_aqu_conc
  real(r8), pointer :: O2_aqu_conc
  real(r8), pointer :: N2_aqu_conc
  real(r8), pointer :: N2O_aqu_conc
  real(r8), pointer :: NH4_1p_conc
  real(r8), pointer :: NH3_aqu_conc
  real(r8), pointer :: Al_3p_conc
  real(r8), pointer :: Fe_3p_conc
  real(r8), pointer :: H_1p_conc
  real(r8), pointer :: Ca_2p_conc
  real(r8), pointer :: Mg_2p_conc
  real(r8), pointer :: Na_1p_conc
  real(r8), pointer :: K_1p_conc
  real(r8), pointer :: OH_1e_conc
  real(r8), pointer :: SO4_2e_conc
  real(r8), pointer :: Cl_e_conc
  real(r8), pointer :: CO3_2e_conc
  real(r8), pointer :: HCO3_e_conc
  real(r8), pointer :: AlOH_2p_conc
  real(r8), pointer :: AlOH2_p_conc
  real(r8), pointer :: AlOH3_conc
  real(r8), pointer :: AlOH4_1e_conc
  real(r8), pointer :: AlSO4_1p_conc
  real(r8), pointer :: FeOH_2p_conc
  real(r8), pointer :: FeO2H2_p_conc
  real(r8), pointer :: FeO3H3_conc
  real(r8), pointer :: FeO4H4_1e_conc
  real(r8), pointer :: FeSO4_1p_conc
  real(r8), pointer :: CaO2H2_conc
  real(r8), pointer :: CaCO3_conc
  real(r8), pointer :: CaHCO3_1p_conc
  real(r8), pointer :: CaSO4_conc
  real(r8), pointer :: MgOH_1p_conc
  real(r8), pointer :: MgCO3_conc
  real(r8), pointer :: MgHCO3_1p_conc
  real(r8), pointer :: MgSO4_conc
  real(r8), pointer :: NaCO3_1e_conc
  real(r8), pointer :: NaSO4_1e_conc
  real(r8), pointer :: KSO4_1e_conc
  real(r8), pointer :: H0PO4_conc
  real(r8), pointer :: H1PO4_2e_conc
  real(r8), pointer :: H2PO4_1e_conc
  real(r8), pointer :: H3PO4_conc
  real(r8), pointer :: FeHPO4_conc
  real(r8), pointer :: FeH2PO4_conc
  real(r8), pointer :: CaPO4_1e_con
  real(r8), pointer :: CaHPO4_conc
  real(r8), pointer :: CaH2PO4_1p_conc
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
  real(r8), pointer :: CALZ
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
  real(r8), pointer :: XNH4_conc
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
  contains

  SUBROUTINE InitSoluteModel(K,BulkSoilMass,ISALTG,solutevar)
!
!     THIS SUBROUTINE INITIALIZES ALL SOIL CHEMISTRY VARIABLES
!

  implicit none
  integer, intent(in) :: K,ISALTG
  real(r8), intent(in) :: BulkSoilMass
  type(solutedtype), target, intent(inout) :: solutevar

  integer, parameter :: MRXN1=1000
  integer :: M
!     begin_execution
  H2CO3_aqu_conc     => solutevar%H2CO3_aqu_conc
  CH4_aqu_conc     => solutevar%CH4_aqu_conc
  O2_aqu_conc     => solutevar%O2_aqu_conc
  N2_aqu_conc     => solutevar%N2_aqu_conc
  N2O_aqu_conc     => solutevar%N2O_aqu_conc
  NH4_1p_conc      => solutevar%NH4_1p_conc
  NH3_aqu_conc      => solutevar%NH3_aqu_conc
  Al_3p_conc      => solutevar%Al_3p_conc
  Fe_3p_conc      => solutevar%Fe_3p_conc
  H_1p_conc      => solutevar%H_1p_conc
  Ca_2p_conc      => solutevar%Ca_2p_conc
  Mg_2p_conc      => solutevar%Mg_2p_conc
  Na_1p_conc      => solutevar%Na_1p_conc
  K_1p_conc      => solutevar%K_1p_conc
  OH_1e_conc      => solutevar%OH_1e_conc
  SO4_2e_conc     => solutevar%SO4_2e_conc
  Cl_e_conc      => solutevar%Cl_e_conc
  CO3_2e_conc     => solutevar%CO3_2e_conc
  HCO3_e_conc    => solutevar%HCO3_e_conc
  AlOH_2p_conc     => solutevar%AlOH_2p_conc
  AlOH2_p_conc     => solutevar%AlOH2_p_conc
  AlOH3_conc     => solutevar%AlOH3_conc
  AlOH4_1e_conc     => solutevar%AlOH4_1e_conc
  AlSO4_1p_conc     => solutevar%AlSO4_1p_conc
  FeOH_2p_conc     => solutevar%FeOH_2p_conc
  FeO2H2_p_conc     => solutevar%FeO2H2_p_conc
  FeO3H3_conc     => solutevar%FeO3H3_conc
  FeO4H4_1e_conc     => solutevar%FeO4H4_1e_conc
  FeSO4_1p_conc     => solutevar%FeSO4_1p_conc
  CaO2H2_conc     => solutevar%CaO2H2_conc
  CaCO3_conc     => solutevar%CaCO3_conc
  CaHCO3_1p_conc     => solutevar%CaHCO3_1p_conc
  CaSO4_conc     => solutevar%CaSO4_conc
  MgOH_1p_conc     => solutevar%MgOH_1p_conc
  MgCO3_conc     => solutevar%MgCO3_conc
  MgHCO3_1p_conc     => solutevar%MgHCO3_1p_conc
  MgSO4_conc     => solutevar%MgSO4_conc
  NaCO3_1e_conc     => solutevar%NaCO3_1e_conc
  NaSO4_1e_conc     => solutevar%NaSO4_1e_conc
  KSO4_1e_conc     => solutevar%KSO4_1e_conc
  H0PO4_conc     => solutevar%H0PO4_conc
  H1PO4_2e_conc     => solutevar%H1PO4_2e_conc
  H2PO4_1e_conc     => solutevar%H2PO4_1e_conc
  H3PO4_conc     => solutevar%H3PO4_conc
  FeHPO4_conc     => solutevar%FeHPO4_conc
  FeH2PO4_conc     => solutevar%FeH2PO4_conc
  CaPO4_1e_con     => solutevar%CaPO4_1e_con
  CaHPO4_conc     => solutevar%CaHPO4_conc
  CaH2PO4_1p_conc     => solutevar%CaH2PO4_1p_conc
  MgHPO4_conc     => solutevar%MgHPO4_conc
  CSTR1     => solutevar%CSTR1
  CCO2M     => solutevar%CCO2M
  CCH4M     => solutevar%CCH4M
  COXYM     => solutevar%COXYM
  CZ2GM     => solutevar%CZ2GM
  CZ2OM     => solutevar%CZ2OM
  CN4Z      => solutevar%CN4Z
  CNOZ      => solutevar%CNOZ
  CNAZ      => solutevar%CNAZ
  CKAZ      => solutevar%CKAZ
  CSOZ      => solutevar%CSOZ
  CCLZ      => solutevar%CCLZ
  CNOX      => solutevar%CNOX
  CCASOX    => solutevar%CCASOX
  CN4X      => solutevar%CN4X
  CPOZ      => solutevar%CPOZ
  CALZ      => solutevar%CALZ
  CFEZ      => solutevar%CFEZ
  CCAZ      => solutevar%CCAZ
  CMGZ      => solutevar%CMGZ
  CALX      => solutevar%CALX
  CFEX      => solutevar%CFEX
  CaX_conc      => solutevar%CaX_conc
  MgX_conc      => solutevar%MgX_conc
  CNAX      => solutevar%CNAX
  CKAX      => solutevar%CKAX
  CSOX      => solutevar%CSOX
  CCLX      => solutevar%CCLX
  CALPOX    => solutevar%CALPOX
  CFEPOX    => solutevar%CFEPOX
  CCAPDX    => solutevar%CCAPDX
  CCAPHX    => solutevar%CCAPHX
  CALOHX    => solutevar%CALOHX
  CFEOHX    => solutevar%CFEOHX
  CCACOX    => solutevar%CCACOX
  XNH4_conc      => solutevar%XNH4_conc
  XHY1      => solutevar%XHY1
  XAl_conc      => solutevar%XAl_conc
  XFe_conc      => solutevar%XFe_conc
  XCa_conc      => solutevar%XCa_conc
  Precp_Ca5P3O12O3H3_conc    => solutevar%Precp_Ca5P3O12O3H3_conc
  XMg_conc      => solutevar%XMg_conc
  XNa_conc      => solutevar%XNa_conc
  XK_conc      => solutevar%XK_conc
  XHC1      => solutevar%XHC1
  XAlO2H2_conc    => solutevar%XAlO2H2_conc
  XFeO2H2_conc    => solutevar%XFeO2H2_conc
  XOH_conc     => solutevar%XOH_conc
  XROH1_conc     => solutevar%XROH1_conc
  XROH2_conc     => solutevar%XROH2_conc
  XHPO4_conc     => solutevar%XHPO4_conc
  XH2PO4_conc     => solutevar%XH2PO4_conc
  Precp_AlO3H3_conc    => solutevar%Precp_AlO3H3_conc
  Precp_FeO3H3_conc    => solutevar%Precp_FeO3H3_conc
  Precp_CaCO3_conc    => solutevar%Precp_CaCO3_conc
  Precp_CaSO4_conc    => solutevar%Precp_CaSO4_conc
  Precp_AlPO4_conc    => solutevar%Precp_AlPO4_conc
  Precp_FePO4_conc    => solutevar%Precp_FePO4_conc
  Precp_CaHPO4_conc    => solutevar%Precp_CaHPO4_conc
  FH2O      => solutevar%FH2O
  ATCA      => solutevar%ATCA
  XAEC      => solutevar%XAEC
  CEC       => solutevar%CEC
  ORGC      => solutevar%ORGC
  VLPO4     => solutevar%VLPO4
  XCEC      => solutevar%XCEC
  GKC4      => solutevar%GKC4
  GKCA      => solutevar%GKCA
  GKCH      => solutevar%GKCH
  GKCK      => solutevar%GKCK
  GKCN      => solutevar%GKCN
  GKCM      => solutevar%GKCM
  ZEROS     => solutevar%ZEROS
  VLWatMicP      => solutevar%VLWatMicP

  call InitEquilibria(K,BulkSoilMass,ISALTG)
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

  subroutine InitEquilibria(K,BulkSoilMass,ISALTG)

  implicit none
  integer, intent(in) :: K,ISALTG
  real(r8), intent(in) :: BulkSoilMass
  integer :: MM
  real(r8) :: XNAQ,XPT,XN4Q,XMGQ,XKAQ,XHYQ,XFEQ
  real(r8) :: XHP,XCECQ,XTLQ,XALQ,XCAX,SPH2P
  real(r8) :: XCAQ,XOH,SPH1P,FX,FXH2,FXP1,FSTR2
  real(r8) :: FHPA,FXH1,FXP2,FHP1,FHP3,FHCO
  real(r8) :: FCO3,CX1,CX2,CSTRZ,FHP2,CCO2Y
  real(r8) :: CION2,A3,CA1,CA2,CA3,CC1,CC2,CC3
  real(r8) :: CCO2X,CCO2Z,CN,CSTR2,FHP0,FXH0,R,Z
!
!     INITIALIZE SOLUTE EQUILIBRIA
!
  CC3=AZMAX1(CALZ)+AZMAX1(CFEZ)
  CA3=0._r8
  CC2=AZMAX1(CCAZ)+CMGZ
  CA2=CSOZ
  CC1=CN4Z+CNAZ+CKAZ
  CA1=CNOZ+CCLZ+CPOZ
  CX2=0._r8
  CX1=0._r8
  CN=0._r8
!
!     INITIALIZE ION STRENGTH AND ACTIVITIES
!
  CION2=AZMAX1(CC3+CA3+CC2+CA2+CC1+CA1+CN)
  CSTR1=0.5E-03*(9.0*(CC3+CA3)+4.0*(CC2+CA2)+CC1+CA1)
  CSTRZ=0.5E-03*(9.0*(CC3+CA3)+4.0*(CC2+CX2)+CC1+CX1)
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

  CCO2X=CCO2M*gas_solubility(idg_CO2,ATCA)/(EXP(ACTCG(idg_CO2)*CSTRZ))*FH2O
  CCO2Y=LOG(CCO2X)
  CCO2Z=ABS(CCO2Y)
  H2CO3_aqu_conc=CCO2X
  FCO3=DPCO3*A0/(H_1p_conc**2._r8*A2)
  FHCO=DPCO2*A0/(H_1p_conc*A1)
  Z=ACTCG(idg_CO2)*(2.0E-03_r8*FCO3+0.5E-03_r8*FHCO)
  DO  MM=1,25
    R=(LOG(H2CO3_aqu_conc)+Z*H2CO3_aqu_conc-CCO2Y)/CCO2Z
    IF(R.LT.1.0E-03_r8)exit
    H2CO3_aqu_conc=H2CO3_aqu_conc/SQRT(1.0_r8+R)
  ENDDO
  CH4_aqu_conc=CCH4M*gas_solubility(idg_CH4,ATCA)/(EXP(ACTCG(idg_CH4)*CSTR1))*FH2O
  O2_aqu_conc=COXYM*gas_solubility(idg_O2,ATCA)/(EXP(ACTCG(idg_O2)*CSTR1))*FH2O
  N2_aqu_conc=CZ2GM*gas_solubility(idg_N2,ATCA)/(EXP(ACTCG(idg_N2)*CSTR1))*FH2O
  N2O_aqu_conc=CZ2OM*gas_solubility(idg_N2O,ATCA)/(EXP(ACTCG(idg_N2O)*CSTR1))*FH2O

  CO3_2e_conc=H2CO3_aqu_conc*DPCO3*A0/(H_1p_conc**2._r8*A2)
  HCO3_e_conc=H2CO3_aqu_conc*DPCO2*A0/(H_1p_conc*A1)
  NO3_1e_conc=CNOZ
!
!     INITIALIZE ION PAIR EQUILIBRIA
!
  IF(K.NE.micpar%k_POM)THEN
    NH4_1p_conc=CN4Z/(1.0_r8+DPN4*A1/(H_1p_conc*A0))
    NH3_aqu_conc=NH4_1p_conc*DPN4*A1/(H_1p_conc*A0)
  ELSE
    NH4_1p_conc=ZERO
    NH3_aqu_conc=ZERO
  ENDIF
  IF(CALZ.LT.0.0_r8)THEN
    Al_3p_conc=AMIN1(CALMX,SPALO/(OH_1e_conc**3*A3))
  ELSE
    Al_3p_conc=AMIN1(CALZ,SPALO/(OH_1e_conc**3*A3))
  ENDIF
  IF(CFEZ.LT.0.0_r8)THEN
    Fe_3p_conc=AMIN1(CFEMX,SPFEO/(OH_1e_conc**3*A3))
  ELSE
    Fe_3p_conc=AMIN1(CFEZ,SPFEO/(OH_1e_conc**3*A3))
  ENDIF
  IF(CCAZ.LT.0.0_r8)THEN
    Ca_2p_conc=AMIN1(CCAMX,SPCAC/(CO3_2e_conc*A2**2))
  ELSE
    Ca_2p_conc=AMIN1(CCAZ,SPCAC/(CO3_2e_conc*A2**2))
  ENDIF
  Mg_2p_conc=CMGZ
  Na_1p_conc=CNAZ
  K_1p_conc=CKAZ
  SO4_2e_conc=CSOZ
  Cl_e_conc=CCLZ
  AlOH_2p_conc=Al_3p_conc*OH_1e_conc*A3/(DPAL1*A2)
  AlOH2_p_conc=Al_3p_conc*OH_1e_conc**2*A3/(DPAL1*DPAL2*A1)
  AlOH3_conc=Al_3p_conc*OH_1e_conc**3*A3/(DPAL1*DPAL2*DPAL3*A0)
  AlOH4_1e_conc=Al_3p_conc*OH_1e_conc**4*A3/(DPAL1*DPAL2*DPAL3*DPAL4*A1)
  AlSO4_1p_conc=0._r8
  FeOH_2p_conc=Fe_3p_conc*OH_1e_conc*A3/(DPFE1*A2)
  FeO2H2_p_conc=Fe_3p_conc*OH_1e_conc**2*A3/(DPFE1*DPFE2*A1)
  FeO3H3_conc=Fe_3p_conc*OH_1e_conc**3*A3/(DPFE1*DPFE2*DPFE3*A0)
  FeO4H4_1e_conc=Fe_3p_conc*OH_1e_conc**4*A3/(DPFE1*DPFE2*DPFE3*DPFE4*A1)
  FeSO4_1p_conc=0._r8
  CaO2H2_conc=Ca_2p_conc*OH_1e_conc*A2/(DPCAO*A1)
  CaCO3_conc=Ca_2p_conc*CO3_2e_conc*A2**2/(DPCAC*A0)
  CaHCO3_1p_conc=Ca_2p_conc*HCO3_e_conc*A2/DPCAH
  CaSO4_conc=0._r8
  MgOH_1p_conc=Mg_2p_conc*OH_1e_conc*A2/(DPMGO*A1)
  MgCO3_conc=Mg_2p_conc*CO3_2e_conc*A2**2/(DPMGC*A0)
  MgHCO3_1p_conc=Mg_2p_conc*HCO3_e_conc*A2/DPMGH
  MgSO4_conc=0._r8
  NaCO3_1e_conc=Na_1p_conc*CO3_2e_conc*A2/DPNAC
  NaSO4_1e_conc=0._r8
  KSO4_1e_conc=0._r8
  FeHPO4_conc=0._r8
  FeH2PO4_conc=0._r8
  CaPO4_1e_con=0._r8
  CaHPO4_conc=0._r8
  CaH2PO4_1p_conc=0._r8
  MgHPO4_conc=0._r8
!
!     INITIALIZE PHOSPHORUS EQUILIBRIA AMONG SOLUBLE, ADSORBED
!     AND PRECIPITATED FORMS
! NOT POM complex
  IF(K.NE.micpar%k_POM)THEN
    H3PO4_conc=CPOZ/(1.0_r8+DPH3P*A0/(H_1p_conc*A1)+DPH3P*DPH2P*A0 &
      /(H_1p_conc**2*A2)+DPH3P*DPH2P*DPH1P*A0/(H_1p_conc**3*A3))
    H2PO4_1e_conc=H3PO4_conc*DPH3P*A0/(H_1p_conc*A1)
    H1PO4_2e_conc=H3PO4_conc*DPH3P*DPH2P*A0/(H_1p_conc**2*A2)
    H0PO4_conc=H3PO4_conc*DPH3P*DPH2P*DPH1P*A0/(H_1p_conc**3*A3)
! POM complex
  ELSE
    XHP=CPOZ
    XOH=XAEC/BulkSoilMass
    FHP3=1.0_r8/(1.0_r8+DPH3P*A0/(H_1p_conc*A1)+DPH3P*DPH2P*A0 &
      /(H_1p_conc**2*A2)+DPH3P*DPH2P*DPH1P*A0/(H_1p_conc**3*A3))
    FHP2=FHP3*DPH3P*A0/(H_1p_conc*A1)
    FHP1=FHP3*DPH3P*DPH2P*A0/(H_1p_conc**2*A2)
    FHP0=FHP3*DPH3P*DPH2P*DPH1P*A0/(H_1p_conc**3*A3)
    SPOH2=SXOH2/A1
    SPOH1=SXOH1/A1
    SPH2P=SXH2P*DPH2O/A1
    SPH1P=SXH1P*DPH2O*A1/A2
    FXH2=1.0_r8/(1.0_r8+SPOH2/H_1p_conc+SPOH2*SPOH1/H_1p_conc**2)
    FXH1=FXH2*SPOH2/H_1p_conc
    FXH0=FXH1*SPOH1/H_1p_conc
    FXP2=1.0_r8/(1.0_r8+SXH2P*DPH2P/(SXH1P*H_1p_conc))
    FXP1=FXP2*SXH2P*DPH2P/(SXH1P*H_1p_conc)
    FHPA=FHP2*A1
    XPT=(XOH+XHP+SXH2P*FXP2*OH_1e_conc/(FXH1*FHPA)-SQRT(XOH**2*FXH1**2 &
      *FHPA**2-2.0*XOH*FXH1**2*XHP*FHPA**2+FXH1**2*XHP**2*FHPA**2 &
      +2.0*XOH*FXH1*FHPA*SXH2P*FXP2*OH_1e_conc+2.0*FXH1*XHP*FHPA*SXH2P &
      *FXP2*OH_1e_conc+SXH2P**2*FXP2**2*OH_1e_conc**2)/(FXH1*FHPA))/2.0_r8
    XROH2_conc=(XOH-XPT)*FXH2
    XROH1_conc=(XOH-XPT)*FXH1
    XOH_conc=(XOH-XPT)*FXH0
    XHPO4_conc=XPT*FXP1
    XH2PO4_conc=XPT*FXP2
    H3PO4_conc=(XHP-XPT)*FHP3
    H2PO4_1e_conc=(XHP-XPT)*FHP2
    H1PO4_2e_conc=(XHP-XPT)*FHP1
    H0PO4_conc=(XHP-XPT)*FHP0
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
    XCOOH=AZMAX1(COOH1*ORGC)
  ENDIF
  CC3=Al_3p_conc+Fe_3p_conc
  CA3=H0PO4_conc
  CC2=Ca_2p_conc+Mg_2p_conc+AlOH_2p_conc+FeOH_2p_conc+FeH2PO4_conc
  CA2=SO4_2e_conc+CO3_2e_conc+H1PO4_2e_conc
  CC1=NH4_1p_conc+H_1p_conc+Na_1p_conc+K_1p_conc+AlOH2_p_conc+FeO2H2_p_conc+AlSO4_1p_conc &
    +FeSO4_1p_conc+CaO2H2_conc+CaHCO3_1p_conc+MgOH_1p_conc+MgHCO3_1p_conc+FeHPO4_conc+CaH2PO4_1p_conc
  CA1=NO3_1e_conc+OH_1e_conc+HCO3_e_conc+Cl_e_conc+AlOH4_1e_conc+FeO4H4_1e_conc+NaCO3_1e_conc &
    +NaSO4_1e_conc+KSO4_1e_conc+H2PO4_1e_conc+CaPO4_1e_con
  CN=H2CO3_aqu_conc+CH4_aqu_conc+O2_aqu_conc+N2_aqu_conc+N2O_aqu_conc+NH3_aqu_conc+AlOH3_conc+FeO3H3_conc+CaCO3_conc+CaSO4_conc &
    +MgCO3_conc+MgSO4_conc+H3PO4_conc+CaHPO4_conc+MgHPO4_conc
  CX2=CA2-CO3_2e_conc
  CX1=CA1-HCO3_e_conc
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
    CALX=Al_3p_conc**0.333_r8
    CFEX=Fe_3p_conc**0.333_r8
    CaX_conc=Ca_2p_conc**0.500
    MgX_conc=Mg_2p_conc**0.500
    XCAX=CEC_conc/(1.0_r8+GKC4*NH4_1p_conc/CaX_conc &
      +GKCH*H_1p_conc/CaX_conc+GKCA*CALX/CaX_conc &
      +GKCA*CFEX/CaX_conc+GKCA*Mg_2p_conc/CaX_conc &
      +GKCN*Na_1p_conc/CaX_conc+GKCK*K_1p_conc/CaX_conc)
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
    XNH4_conc=CN4X
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
  real(r8) :: H2PO4_e_to_HPO4_2e_flx,RHP1,RHP2,H2PO4_1e_apatite_dissol_flx,H2PO4_1e_FePO4_dissol_flx,H1PO4_to_XHPO4_ROH_flx
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
  real(r8) :: A3,AAL1,AALO1,AALO2,AALO3,AALO4
  real(r8) :: AALS1,AALX,AC0P1,AC1P1,AC2P1,ACA1,ACAC1
  real(r8) :: ACAH1,ACAS1,ACAO1,ACAX,ACO21,ACO31
  real(r8) :: AF1P1,AF2P1,AFE1,AFEO1,AFEO2,AFEO3,AFEO4
  real(r8) :: AFES1,AFEX,AH0P1,AH1P1,AH2P1,AH3P1,AHCO31
  real(r8) :: AHY1,AKA1,AKAS1,AM1P1,AMG1,AMGC1,AMGO1
  real(r8) :: AMGX,AMGH1,AMGS1,AN31,AN41,ANA1,ANAC1
  real(r8) :: AOH1,ASO41,ANAS1,CA1,CA2,CA3,CC1,CC2
  real(r8) :: CC3,CCO2X,CCO2Z,CN,CSTR2,DP,P1,P2,P3,PX
  real(r8) :: PY,R1,RAL,RAL1,RAL2,RAL3,RAL4,RALO1
  real(r8) :: RALO2,RALO3,RALO4,RALS,RC0P,RC1P,RC2P
  integer :: NR1,NP2,NP3
  integer :: MM

! begin_execution
  H2CO3_aqu_conc=AMAX1(ZERO,H2CO3_aqu_conc)
  CO3_2e_conc=H2CO3_aqu_conc*DPCO3*A0/(H_1p_conc**2*A2)
  HCO3_e_conc=H2CO3_aqu_conc*DPCO2*A0/(H_1p_conc*A1)
  NH4_1p_conc=AMAX1(ZERO,NH4_1p_conc)
  NH3_aqu_conc=AMAX1(ZERO,NH3_aqu_conc)
  Al_3p_conc=AMAX1(ZERO,Al_3p_conc)
  Fe_3p_conc=AMAX1(ZERO,Fe_3p_conc)
  Ca_2p_conc=AMAX1(ZERO,Ca_2p_conc)
  Ca_2p_conc=AMIN1(Ca_2p_conc,SPCAC/(CO3_2e_conc*A2**2))
  Mg_2p_conc=AMAX1(ZERO,Mg_2p_conc)
  Na_1p_conc=AMAX1(ZERO,Na_1p_conc)
  K_1p_conc=AMAX1(ZERO,K_1p_conc)
  SO4_2e_conc=AMAX1(ZERO,SO4_2e_conc)
  AlOH_2p_conc=AMAX1(ZERO,AlOH_2p_conc)
  AlOH2_p_conc=AMAX1(ZERO,AlOH2_p_conc)
  AlOH3_conc=AMAX1(ZERO,AlOH3_conc)
  AlOH4_1e_conc=AMAX1(ZERO,AlOH4_1e_conc)
  AlSO4_1p_conc=AMAX1(ZERO,AlSO4_1p_conc)
  FeOH_2p_conc=AMAX1(ZERO,FeOH_2p_conc)
  FeO2H2_p_conc=AMAX1(ZERO,FeO2H2_p_conc)
  FeO3H3_conc=AMAX1(ZERO,FeO3H3_conc)
  FeO4H4_1e_conc=AMAX1(ZERO,FeO4H4_1e_conc)
  FeSO4_1p_conc=AMAX1(ZERO,FeSO4_1p_conc)
  CaO2H2_conc=AMAX1(ZERO,CaO2H2_conc)
  CaCO3_conc=AMAX1(ZERO,CaCO3_conc)
  CaHCO3_1p_conc=AMAX1(ZERO,CaHCO3_1p_conc)
  CaSO4_conc=AMAX1(ZERO,CaSO4_conc)
  MgOH_1p_conc=AMAX1(ZERO,MgOH_1p_conc)
  MgCO3_conc=AMAX1(ZERO,MgCO3_conc)
  MgHCO3_1p_conc=AMAX1(ZERO,MgHCO3_1p_conc)
  MgSO4_conc=AMAX1(ZERO,MgSO4_conc)
  NaCO3_1e_conc=AMAX1(ZERO,NaCO3_1e_conc)
  NaSO4_1e_conc=AMAX1(ZERO,NaSO4_1e_conc)
  KSO4_1e_conc=AMAX1(ZERO,KSO4_1e_conc)
  H0PO4_conc=AMAX1(ZERO,H0PO4_conc)
  H1PO4_2e_conc=AMAX1(ZERO,H1PO4_2e_conc)
  H2PO4_1e_conc=AMAX1(ZERO,H2PO4_1e_conc)
  H3PO4_conc=AMAX1(ZERO,H3PO4_conc)
  FeHPO4_conc=AMAX1(ZERO,FeHPO4_conc)
  FeH2PO4_conc=AMAX1(ZERO,FeH2PO4_conc)
  CaPO4_1e_con=AMAX1(ZERO,CaPO4_1e_con)
  CaHPO4_conc=AMAX1(ZERO,CaHPO4_conc)
  CaH2PO4_1p_conc=AMAX1(ZERO,CaH2PO4_1p_conc)
  MgHPO4_conc=AMAX1(ZERO,MgHPO4_conc)
!
!     ION ACTIVITY COEFFICIENTS
!
  CC3=Al_3p_conc+Fe_3p_conc
  CA3=H0PO4_conc
  CC2=Ca_2p_conc+Mg_2p_conc+AlOH_2p_conc+FeOH_2p_conc+FeH2PO4_conc
  CA2=SO4_2e_conc+CO3_2e_conc+H1PO4_2e_conc
  CC1=NH4_1p_conc+H_1p_conc+Na_1p_conc+K_1p_conc+AlOH2_p_conc+FeO2H2_p_conc+AlSO4_1p_conc+FeSO4_1p_conc+CaO2H2_conc &
    +CaHCO3_1p_conc+MgOH_1p_conc+MgHCO3_1p_conc+FeHPO4_conc+CaH2PO4_1p_conc
  CA1=NO3_1e_conc+OH_1e_conc+HCO3_e_conc+Cl_e_conc+AlOH4_1e_conc+FeO4H4_1e_conc+NaCO3_1e_conc+NaSO4_1e_conc+KSO4_1e_conc &
    +H2PO4_1e_conc+CaPO4_1e_con
  CN=H2CO3_aqu_conc+CH4_aqu_conc+O2_aqu_conc+N2_aqu_conc+N2O_aqu_conc+NH3_aqu_conc+AlOH3_conc+FeO3H3_conc+CaCO3_conc+CaSO4_conc &
    +MgCO3_conc+MgSO4_conc+H3PO4_conc+CaHPO4_conc+MgHPO4_conc
  CX2=CA2-CO3_2e_conc
  CX1=CA1-HCO3_e_conc
  CION2=AZMAX1(CC3+CA3+CC2+CA2+CC1+CA1+CN)
  CSTR1=0.5E-03*(9.0*(CC3+CA3)+4.0*(CC2+CA2)+CC1+CA1)
  CSTRZ=0.5E-03*(9.0*(CC3+CA3)+4.0*(CC2+CX2)+CC1+CX1)
  CSTR2=SQRT(CSTR1)
  FSTR2=CSTR2/(1.0_r8+CSTR2)
  FH2O=5.56E+04/(5.56E+04+CION2)
  A1=AMIN1(1.0,10.0**(-0.509*1.0*FSTR2+0.20*CSTR2))
  A2=AMIN1(1.0,10.0**(-0.509*4.0*FSTR2+0.20*CSTR2))
  A3=AMIN1(1.0,10.0**(-0.509*9.0*FSTR2+0.20*CSTR2))
!
!     PRECIPITATION-DISSOLUTION EQUILIBRIA
!
  AHY1=H_1p_conc*A1
  AOH1=OH_1e_conc*A1
  AAL1=Al_3p_conc*A3
  AALO1=AlOH_2p_conc*A2
  AALO2=AlOH2_p_conc*A1
  AALO3=AlOH3_conc
  AALO4=AlOH4_1e_conc*A1
  AFE1=Fe_3p_conc*A3
  AFEO1=FeOH_2p_conc*A2
  AFEO2=FeO2H2_p_conc*A1
  AFEO3=FeO3H3_conc
  AFEO4=FeO4H4_1e_conc*A1
  ACA1=Ca_2p_conc*A2
  ACO31=CO3_2e_conc*A2
  AHCO31=HCO3_e_conc*A1
  ACO21=H2CO3_aqu_conc*A0
  ASO41=SO4_2e_conc*A2
  AH0P1=H0PO4_conc*A3
  AH1P1=H1PO4_2e_conc*A2
  AH2P1=H2PO4_1e_conc*A1
  AH3P1=H3PO4_conc*A0
  AF1P1=FeHPO4_conc*A2
  AF2P1=FeH2PO4_conc*A1
  AC0P1=CaPO4_1e_con*A1
  AC1P1=CaHPO4_conc*A0
  AC2P1=CaH2PO4_1p_conc*A1
  AM1P1=MgHPO4_conc*A0
  AN41=NH4_1p_conc*A1
  AN31=NH3_aqu_conc*A0
  AMG1=Mg_2p_conc*A2
  ANA1=Na_1p_conc*A1
  AKA1=K_1p_conc*A1
  AALX=AAL1**0.333
  AFEX=AFE1**0.333
  ACAX=ACA1**0.500
  AMGX=AMG1**0.500
  AALS1=AlSO4_1p_conc*A1
  AFES1=FeSO4_1p_conc*A1
  ACAO1=CaO2H2_conc*A1
  ACAC1=CaCO3_conc*A0
  ACAS1=CaSO4_conc*A0
  ACAH1=CaHCO3_1p_conc*A1
  AMGO1=MgOH_1p_conc*A1
  AMGC1=MgCO3_conc*A0
  AMGH1=MgHCO3_1p_conc*A1
  AMGS1=MgSO4_conc*A0
  ANAC1=NaCO3_1e_conc*A1
  ANAS1=NaSO4_1e_conc*A1
  AKAS1=KSO4_1e_conc*A1
!
!     ALUMINUM HYDROXIDE (GIBBSITE)
!
  IF(K.EQ.micpar%k_POM)THEN
    PX=AMAX1(AAL1,AALO1,AALO2,AALO3,AALO4)
    IF(isclose(PX,AAL1))THEN
      R1=AHY1
      P1=AAL1
      P2=AOH1
      NR1=3
      NP2=0
      SP=SHALO
    ELSEIF(isclose(PX,AALO1))THEN
      R1=AHY1
      P1=AALO1
      P2=AOH1
      NR1=2
      NP2=0
      SP=SHAL1
    ELSEIF(isclose(PX,AALO2))THEN
      R1=AHY1
      P1=AALO2
      P2=AOH1
      NR1=1
      NP2=0
      SP=SHAL2
    ELSEIF(isclose(PX,AALO3))THEN
      R1=AHY1
      P1=AALO3
      P2=AOH1
      NR1=0
      NP2=0
      SP=SPAL3
    ELSEIF(isclose(PX,AALO4))THEN
      R1=AOH1
      P1=AALO4
      P2=AHY1
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
    IF(isclose(PX,AAL1))THEN
      RHAL1=RPALOX
    ELSEIF(isclose(PX,AALO1))THEN
      RHALO1=RPALOX
    ELSEIF(isclose(PX,AALO2))THEN
      RHALO2=RPALOX
    ELSEIF(isclose(PX,AALO3))THEN
      RHALO3=RPALOX
    ELSEIF(isclose(PX,AALO4))THEN
      RHALO4=RPALOX
    ENDIF
!
!     IRON HYDROXIDE
!
    PX=AMAX1(AFE1,AFEO1,AFEO2,AFEO3,AFEO4)
    IF(isclose(PX,AFE1))THEN
      R1=AHY1
      P1=AFE1
      P2=AOH1
      NR1=3
      NP2=0
      SP=SHFEO
    ELSEIF(isclose(PX,AFEO1))THEN
      R1=AHY1
      P1=AFEO1
      P2=AOH1
      NR1=2
      NP2=0
      SP=SHFE1
    ELSEIF(isclose(PX,AFEO2))THEN
      R1=AHY1
      P1=AFEO2
      P2=AOH1
      NR1=1
      NP2=0
      SP=SHFE2
    ELSEIF(isclose(PX,AFEO3))THEN
      R1=AHY1
      P1=AFEO3
      P2=AOH1
      NR1=0
      NP2=0
      SP=SPFE3
    ELSEIF(isclose(PX,AFEO4))THEN
      R1=AOH1
      P1=AFEO4
      P2=AHY1
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
    IF(isclose(PX,AFE1))THEN
      RHFE1=RPFEOX
    ELSEIF(isclose(PX,AFEO1))THEN
      RHFEO1=RPFEOX
    ELSEIF(isclose(PX,AFEO2))THEN
      RHFEO2=RPFEOX
    ELSEIF(isclose(PX,AFEO3))THEN
      RHFEO3=RPFEOX
    ELSEIF(isclose(PX,AFEO4))THEN
      RHFEO4=RPFEOX
    ENDIF
!
!     CALCITE
!
    PX=AMAX1(ACO31,AHCO31,ACO21)
    R1=AHY1
    P1=ACA1
    IF(isclose(PX,ACO31))THEN
      P2=ACO31
      NR1=0
      SP=SPCAC
    ELSEIF(isclose(PX,AHCO31))THEN
      P2=AHCO31
      NR1=1
      SP=SHCAC1
    ELSEIF(isclose(PX,ACO21))THEN
      P2=ACO21
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
    IF(isclose(PX,ACO31))THEN
      RHCAC3=RPCACX
    ELSEIF(isclose(PX,AHCO31))THEN
      RHCACH=RPCACX
    ELSEIF(isclose(PX,ACO21))THEN
      RHCACO=RPCACX
    ENDIF
!
!     GYPSUM
!
    P1=ACA1
    P2=ASO41
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
    PX=AMAX1(AAL1,AALO1,AALO2,AALO3,AALO4)
    PY=AMAX1(AH1P1,AH2P1)
    R1=AHY1
    P3=AHY1
    IF(isclose(PY,AH1P1))THEN
      P2=AH1P1
      IF(isclose(PX,AAL1))THEN
        P1=AAL1
        NR1=1
        NP3=0
        SP=SHA0P1
      ELSEIF(isclose(PX,AALO1))THEN
        P1=AALO1
        NR1=0
        NP3=0
        SP=SPA1P1
      ELSEIF(isclose(PX,AALO2))THEN
        P1=AALO2
        NR1=0
        NP3=1
        SP=SHA2P1
      ELSEIF(isclose(PX,AALO3))THEN
        P1=AALO3
        NR1=0
        NP3=2
        SP=SHA3P1
      ELSEIF(isclose(PX,AALO4))THEN
        P1=AALO4
        NR1=0
        NP3=3
        SP=SHA4P1
      ENDIF
    ELSE
      P2=AH2P1
      IF(isclose(PX,AAL1))THEN
        P1=AAL1
        NR1=2
        NP3=0
        SP=SHA0P2
      ELSEIF(isclose(PX,AALO1))THEN
        P1=AALO1
        NR1=1
        NP3=0
        SP=SHA1P2
      ELSEIF(isclose(PX,AALO2))THEN
        P1=AALO2
        NR1=0
        NP3=0
        SP=SPA2P2
      ELSEIF(isclose(PX,AALO3))THEN
        P1=AALO3
        NR1=0
        NP3=1
        SP=SHA3P2
      ELSEIF(isclose(PX,AALO4))THEN
        P1=AALO4
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
    IF(isclose(PY,AH1P1))THEN
      IF(isclose(PX,AAL1))THEN
        RHA0P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AALO1))THEN
        RHA1P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AALO2))THEN
        RHA2P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AALO3))THEN
        RHA3P1=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AALO4))THEN
        RHA4P1=H2PO4_1e_AlPO4_dissol_flx
      ENDIF
    ELSE
      IF(isclose(PX,AAL1))THEN
        RHA0P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AALO1))THEN
        RHA1P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AALO2))THEN
        RHA2P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AALO3))THEN
        RHA3P2=H2PO4_1e_AlPO4_dissol_flx
      ELSEIF(isclose(PX,AALO4))THEN
        RHA4P2=H2PO4_1e_AlPO4_dissol_flx
      ENDIF
    ENDIF
!
!     IRON PHOSPHATE (STRENGITE)
!
    PX=AMAX1(AFE1,AFEO1,AFEO2,AFEO3,AFEO4)
    PY=AMAX1(AH1P1,AH2P1)
    R1=AHY1
    P3=AHY1
    IF(isclose(PY,AH1P1))THEN
      P2=AH1P1
      IF(isclose(PX,AFE1))THEN
        P1=AFE1
        NR1=1
        NP3=0
        SP=SHF0P1
      ELSEIF(isclose(PX,AFEO1))THEN
        P1=AFEO1
        NR1=0
        NP3=0
        SP=SPF1P1
      ELSEIF(isclose(PX,AFEO2))THEN
        P1=AFEO2
        NR1=0
        NP3=1
        SP=SHF2P1
      ELSEIF(isclose(PX,AFEO3))THEN
        P1=AFEO3
        NR1=0
        NP3=2
        SP=SHF3P1
      ELSEIF(isclose(PX,AFEO4))THEN
        P1=AFEO4
        NR1=0
        NP3=3
        SP=SHF4P1
      ENDIF
    ELSE
      P2=AH2P1
      IF(isclose(PX,AFE1))THEN
        P1=AFE1
        NR1=2
        NP3=0
        SP=SHF0P2
      ELSEIF(isclose(PX,AFEO1))THEN
        P1=AFEO1
        NR1=1
        NP3=0
        SP=SHF1P2
      ELSEIF(isclose(PX,AFEO2))THEN
        P1=AFEO2
        NR1=0
        NP3=0
        SP=SPF2P2
      ELSEIF(isclose(PX,AFEO3))THEN
        P1=AFEO3
        NR1=0
        NP3=1
        SP=SHF3P2
      ELSEIF(isclose(PX,AFEO4))THEN
        P1=AFEO4
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
    IF(isclose(PY,AH1P1))THEN
      IF(isclose(PX,AFE1))THEN
        RHF0P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,AFEO1))THEN
        RHF1P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,AFEO2))THEN
        RHF2P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,AFEO3))THEN
        RHF3P1=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,AFEO4))THEN
        RHF4P1=H2PO4_1e_FePO4_dissol_flx
      ENDIF
    ELSE
      IF(isclose(PX,AFE1))THEN
        RHF0P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,AFEO1))THEN
        RHF1P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,AFEO2))THEN
        RHF2P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,AFEO3))THEN
        RHF3P2=H2PO4_1e_FePO4_dissol_flx
      ELSEIF(isclose(PX,AFEO4))THEN
        RHF4P2=H2PO4_1e_FePO4_dissol_flx
      ENDIF
    ENDIF
!
!     DICALCIUM PHOSPHATE
!
    PX=AMAX1(AH1P1,AH2P1)
    R1=AHY1
    P1=ACA1
    IF(isclose(PX,AH1P1))THEN
      P2=AH1P1
      NR1=0
      SP=SPCAD
    ELSEIF(isclose(PX,AH2P1))THEN
      P2=AH2P1
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
    IF(isclose(PX,AH1P1))THEN
      RPCAD1=H2PO4_1e_CaHPO4_dissol_flx
    ELSEIF(isclose(PX,AH2P1))THEN
      RHCAD2=H2PO4_1e_CaHPO4_dissol_flx
    ENDIF
!
!     HYDROXYAPATITE
!
    PX=AMAX1(AH1P1,AH2P1)
    R1=AHY1
    P1=ACA1
    IF(isclose(PX,AH1P1))THEN
      P2=AH1P1
      NR1=4
      SP=SHCAH1
    ELSEIF(isclose(PX,AH2P1))THEN
      P2=AH2P1
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
    IF(isclose(PX,AH1P1))THEN
      RHCAH1=H2PO4_1e_apatite_dissol_flx
    ELSEIF(isclose(PX,AH2P1))THEN
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
      RXOH2=TAD*(XROH1_conc*AHY1-SXOH2*XROH2_conc)/(XROH1_conc+SPOH2)*VLWatMicPBK
      RXOH1=TAD*(XOH_conc*AHY1-SXOH1*XROH1_conc)/(XOH_conc+SPOH1)*VLWatMicPBK
      SPH2P=SXH2P*DPH2O
      H2PO4_1e_to_XH2PO4_ROH2_flx=TAD*(XROH2_conc*AH2P1-SPH2P*XH2PO4_conc)/(XROH2_conc+SPH2P)*VLWatMicPBK
      H2PO4_to_XH2PO4_ROH_flx=TAD*(XROH1_conc*AH2P1-SXH2P*XH2PO4_conc*AOH1)/(XROH1_conc+SXH2P*AOH1)*VLWatMicPBK
!
!     HPO4 EXCHANGE
!
      SPH1P=SXH1P*DPH2O/DPH2P
      H1PO4_to_XHPO4_ROH_flx=TAD*(XROH1_conc*AH2P1-SPH1P*XHPO4_conc)/(XROH1_conc+SPH1P)*VLWatMicPBK
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
      AALX=AAL1**0.333_r8
      AFEX=AFE1**0.333_r8
      ACAX=ACA1**0.500_r8
      AMGX=AMG1**0.500_r8
      XCAX=CEC_conc/(1.0_r8+GKC4*AN41/ACAX &
       +GKCH*AHY1/ACAX+GKCA*AALX/ACAX &
       +GKCA*AFEX/ACAX+GKCM*AMGX/ACAX &
       +GKCN*ANA1/ACAX+GKCK*AKA1/ACAX)
      XN4Q=XCAX*AN41*GKC4
      XHYQ=XCAX*AHY1*GKCH
      XALQ=XCAX*AALX*GKCA
      XFEQ=XCAX*AFEX*GKCA
      XCAQ=XCAX*ACAX
      XMGQ=XCAX*AMGX*GKCM
      XNAQ=XCAX*ANA1*GKCN
      XKAQ=XCAX*AKA1*GKCK
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
      RXN4=TAD*AMIN1((XN4Q-XNH4_conc)*AN41/XN4Q,NH4_1p_conc)
      RXHY=TAD*AMIN1((XHYQ-XHY1)*AHY1/XHYQ,H_1p_conc)
      RXAL=TAD*AMIN1((XALQ-XAl_conc)*AALX/XALQ,Al_3p_conc)
      RXFE=TAD*AMIN1((XFEQ-XFe_conc)*AFEX/XFEQ,Fe_3p_conc)
      RXCA=TAD*AMIN1((XCAQ-XCa_conc)*ACAX/XCAQ,Ca_2p_conc)
      RXMG=TAD*AMIN1((XMGQ-XMg_conc)*AMGX/XMGQ,Mg_2p_conc)
      RXNA=TAD*AMIN1((XNAQ-XNa_conc)*ANA1/XNAQ,Na_1p_conc)
      RXKA=TAD*AMIN1((XKAQ-XK_conc)*AKA1/XKAQ,K_1p_conc)
      XNH4_conc=XNH4_conc+RXN4
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
    XHC1=AHY1*(XCOOH-XAlO2H2_conc-XFeO2H2_conc)/(AHY1+DPCOH)
    XAlO2H2_conc=AALO2*(XCOOH-XHC1)/(AALO2+DPALO)
    XFeO2H2_conc=AFEO2*(XCOOH-XHC1)/(AFEO2+DPFEO)
    XCOO=AZMAX1(XCOOH-XHC1-XAlO2H2_conc-XFeO2H2_conc)
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
  S0=AHY1+AN31+DPN4
  S1=S0**2_r8-4.0_r8*(AHY1*AN31-DPN4*AN41)
  RNH4=TSL*(S0-SQRT(S1))
  S0=AAL1+AOH1+DPAL1
  S1=S0**2_r8-4.0_r8*(AAL1*AOH1-DPAL1*AALO1)
  RALO1=TSL*(S0-SQRT(S1))
  S0=AALO1+AOH1+DPAL2
  S1=S0**2_r8-4.0_r8*(AALO1*AOH1-DPAL2*AALO2)
  RALO2=TSL*(S0-SQRT(S1))
  S0=AlOH2_p_conc+OH_1e_conc+DPAL3
  S1=S0**2-4.0*(AALO2*AOH1-DPAL3*AALO3)
  RALO3=TSL*(S0-SQRT(S1))
  S0=AALO3+AOH1+DPAL4
  S1=S0**2-4.0*(AALO3*AOH1-DPAL4*AALO4)
  RALO4=TSL*(S0-SQRT(S1))
  S0=AAL1+ASO41+DPALS
  S1=S0**2-4.0*(AAL1*ASO41-DPALS*AALS1)
  RALS=TSL*(S0-SQRT(S1))
  S0=AFE1+AOH1+DPFE1
  S1=S0**2-4.0*(AFE1*AOH1-DPFE1*AFEO1)
  RFEO1=TSL*(S0-SQRT(S1))
  S0=AFEO1+AOH1+DPFE2
  S1=S0**2-4.0*(AFEO1*AOH1-DPFE2*AFEO2)
  RFEO2=TSL*(S0-SQRT(S1))
  S0=AFEO2+AOH1+DPFE3
  S1=S0**2-4.0*(AFEO2*AOH1-DPFE3*AFEO3)
  RFEO3=TSL*(S0-SQRT(S1))
  S0=AFEO3+AOH1+DPFE4
  S1=S0**2-4.0*(AFEO3*AOH1-DPFE4*AFEO4)
  RFEO4=TSL*(S0-SQRT(S1))
  S0=AFE1+ASO41+DPFES
  S1=S0**2-4.0*(AFE1*ASO41-DPFES*AFES1)
  RFES=TSL*(S0-SQRT(S1))
  S0=ACA1+AOH1+DPCAO
  S1=S0**2-4.0*(ACA1*AOH1-DPCAO*ACAO1)
  RCAO=TSL*(S0-SQRT(S1))
  S0=ACA1+ACO31+DPCAC
  S1=S0**2-4.0*(ACA1*ACO31-DPCAC*ACAC1)
  RCAC=TSL*(S0-SQRT(S1))
  S0=ACA1+AHCO31+DPCAH
  S1=S0**2-4.0*(ACA1*AHCO31-DPCAH*ACAH1)
  RCAH=TSL*(S0-SQRT(S1))
  S0=ACA1+ASO41+DPCAS
  S1=S0**2-4.0*(ACA1*ASO41-DPCAS*ACAS1)
  RCAS=TSL*(S0-SQRT(S1))
  S0=AMG1+AOH1+DPMGO
  S1=S0**2-4.0*(AMG1*AOH1-DPMGO*AMGO1)
  RMGO=TSL*(S0-SQRT(S1))
  S0=AMG1+ACO31+DPMGC
  S1=S0**2-4.0*(AMG1*ACO31-DPMGC*AMGC1)
  RMGC=TSL*(S0-SQRT(S1))
  S0=AMG1+AHCO31+DPMGH
  S1=S0**2-4.0*(AMG1*AHCO31-DPMGH*AMGH1)
  RMGH=TSL*(S0-SQRT(S1))
  S0=AMG1+ASO41+DPMGS
  S1=S0**2-4.0*(AMG1*ASO41-DPMGS*AMGS1)
  RMGS=TSL*(S0-SQRT(S1))
  S0=ANA1+ACO31+DPNAC
  S1=S0**2-4.0*(ANA1*ACO31-DPNAC*ANAC1)
  RNAC=TSL*(S0-SQRT(S1))
  S0=ANA1+ASO41+DPNAS
  S1=S0**2-4.0*(ANA1*ASO41-DPNAS*ANAS1)
  RNAS=TSL*(S0-SQRT(S1))
  S0=AKA1+ASO41+DPKAS
  S1=S0**2-4.0*(AKA1*ASO41-DPKAS*AKAS1)
  RKAS=TSL*(S0-SQRT(S1))
  S0=AH0P1+AHY1+DPH1P
  S1=S0**2-4.0*(AH0P1*AHY1-DPH1P*AH1P1)
  RH1P=TSL*(S0-SQRT(S1))
  S0=AH1P1+AHY1+DPH2P
  S1=S0**2-4.0*(AH1P1*AHY1-DPH2P*AH2P1)
  H2PO4_e_to_HPO4_2e_flx=TSL*(S0-SQRT(S1))
  S0=AH2P1+AHY1+DPH3P
  S1=S0**2-4.0*(AH2P1*AHY1-DPH3P*AH3P1)
  RH3P=TSL*(S0-SQRT(S1))
  S0=AFE1+AH1P1+DPF1P
  S1=S0**2-4.0*(AFE1*AH1P1-DPF1P*AF1P1)
  RF1P=TSL*(S0-SQRT(S1))
  S0=AFE1+AH2P1+DPF2P
  S1=S0**2-4.0*(AFE1*AH2P1-DPF2P*AF2P1)
  RF2P=TSL*(S0-SQRT(S1))
  S0=ACA1+AH0P1+DPC0P
  S1=S0**2-4.0*(ACA1*AH0P1-DPC0P*AC0P1)
  RC0P=TSL*(S0-SQRT(S1))
  S0=ACA1+AH1P1+DPC1P
  S1=S0**2-4.0*(ACA1*AH1P1-DPC1P*AC1P1)
  RC1P=TSL*(S0-SQRT(S1))
  S0=ACA1+AH2P1+DPC2P
  S1=S0**2-4.0*(ACA1*AH2P1-DPC2P*AC2P1)
  RC2P=TSL*(S0-SQRT(S1))
  S0=AMG1+AH1P1+DPM1P
  S1=S0**2-4.0*(AMG1*AH1P1-DPM1P*AM1P1)
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
  CCO2X=CCO2M*SCO2X/(EXP(ACTCG(idg_CO2)*CSTRZ))*EXP(0.843-0.0281*ATCA)*FH2O
  CCO2Y=LOG(CCO2X)
  CCO2Z=ABS(CCO2Y)
  H2CO3_aqu_conc=CCO2X
  FCO3=DPCO3*A0/(AHY1**2*A2)
  FHCO=DPCO2*A0/(AHY1*A1)
  Z=ACTCG(idg_CO2)*(2.0E-03*FCO3+0.5E-03*FHCO)
  DO  MM=1,25
    R=(LOG(H2CO3_aqu_conc)+Z*H2CO3_aqu_conc-CCO2Y)/CCO2Z
    IF(R.LT.1.0E-03)exit
    H2CO3_aqu_conc=H2CO3_aqu_conc/SQRT(1.0_r8+R)
  ENDDO
  CH4_aqu_conc=CCH4M*SCH4X/(EXP(ACTCG(idg_CH4)*CSTR1))*EXP(0.597-0.0199*ATCA)*FH2O
  O2_aqu_conc=COXYM*SOXYX/(EXP(ACTCG(idg_O2)*CSTR1))*EXP(0.516-0.0172*ATCA)*FH2O
  N2_aqu_conc=CZ2GM*SN2GX/(EXP(ACTCG(idg_N2)*CSTR1))*EXP(0.456-0.0152*ATCA)*FH2O
  N2O_aqu_conc=CZ2OM*SN2OX/(EXP(ACTCG(idg_N2O)*CSTR1))*EXP(0.897-0.0299*ATCA)*FH2O
  NH4_1p_conc=NH4_1p_conc+RN4S
  NH3_aqu_conc=NH3_aqu_conc+RN3S
  Al_3p_conc=Al_3p_conc+RAL
  Fe_3p_conc=Fe_3p_conc+RFE
  Ca_2p_conc=Ca_2p_conc+RCA
  Mg_2p_conc=Mg_2p_conc+RMG
  Na_1p_conc=Na_1p_conc+RNA
  K_1p_conc=K_1p_conc+RKA
  SO4_2e_conc=SO4_2e_conc+RSO4
  CO3_2e_conc=H2CO3_aqu_conc*DPCO3*A0/(AHY1**2*A2)
  HCO3_e_conc=H2CO3_aqu_conc*DPCO2*A0/(AHY1*A1)
  AlOH_2p_conc=AlOH_2p_conc+RAL1
  AlOH2_p_conc=AlOH2_p_conc+RAL2
  AlOH3_conc=AlOH3_conc+RAL3
  AlOH4_1e_conc=AlOH4_1e_conc+RAL4
  AlSO4_1p_conc=AlSO4_1p_conc+RALS
  FeOH_2p_conc=FeOH_2p_conc+RFE1
  FeO2H2_p_conc=FeO2H2_p_conc+RFE2
  FeO3H3_conc=FeO3H3_conc+RFE3
  FeO4H4_1e_conc=FeO4H4_1e_conc+RFE4
  FeSO4_1p_conc=FeSO4_1p_conc+RFES
  CaO2H2_conc=CaO2H2_conc+RCAO
  CaCO3_conc=CaCO3_conc+RCAC
  CaHCO3_1p_conc=CaHCO3_1p_conc+RCAH
  CaSO4_conc=CaSO4_conc+RCAS
  MgOH_1p_conc=MgOH_1p_conc+RMGO
  MgCO3_conc=MgCO3_conc+RMGC
  MgHCO3_1p_conc=MgHCO3_1p_conc+RMGH
  MgSO4_conc=MgSO4_conc+RMGS
  NaCO3_1e_conc=NaCO3_1e_conc+RNAC
  NaSO4_1e_conc=NaSO4_1e_conc+RNAS
  KSO4_1e_conc=KSO4_1e_conc+RKAS
  H0PO4_conc=H0PO4_conc+RHP0
  H1PO4_2e_conc=H1PO4_2e_conc+RHP1
  H2PO4_1e_conc=H2PO4_1e_conc+RHP2
  H3PO4_conc=H3PO4_conc+RHP3
  FeHPO4_conc=FeHPO4_conc+RF1P
  FeH2PO4_conc=FeH2PO4_conc+RF2P
  CaPO4_1e_con=CaPO4_1e_con+RC0P
  CaHPO4_conc=CaHPO4_conc+RC1P
  CaH2PO4_1p_conc=CaH2PO4_1p_conc+RC2P
  MgHPO4_conc=MgHPO4_conc+RM1P
!  IF(K.EQ.micpar%k_POM.AND.(M/1)*1.EQ.M)THEN
!     WRITE(*,1112)'A1I',I,NX,NY,L,K,M,A1,A2,A3,FSTR2,CSTR1
!    2,CSTR2,CC3,CA3,CC2,CA2,CC1,CA1,VLWatMicP
!     WRITE(*,1112)'ALPO4I',I,NX,NY,L,K,M,Precp_AlPO4_conc,AAL1
!    2,AALO1,AALO2,AALO3,AALO4
!    2,AH0P1,AH1P1,AH2P1,AHY1,AOH1,H2PO4_1e_AlPO4_dissol_flx,RHA0P1,RHA1P1,RHA2P1,RHA3P1
!    3,RHA4P1,RHA0P2,RHA1P2,RHA2P2,RHA3P2,RHA4P2,SP,SPX,AAL1*AH0P1
!    4,SPALP,H0PO4_conc,H1PO4_2e_conc,H2PO4_1e_conc,H3PO4_conc,RHP0,RHP1,RHP2,RHP3
!    5,RAL,RHAL1,RHA0P1,RHA0P2,RALO1,RALS,RXAL
!     WRITE(*,1112)'FEPO4I',I,NX,NY,L,K,M,Precp_FePO4_conc,AFE1
!    2,AFEO1,AFEO2,AFEO3,AFEO4
!    2,AH0P1,AH1P1,AH2P1,AHY1,AOH1,H2PO4_1e_FePO4_dissol_flx,RHF0P1,RHF1P1,RHF2P1,RHF3P1
!    3,RHF4P1,RHF0P2,RHF1P2,RHF2P2,RHF3P2,RHF4P2,SP,SPX,AFE1*AH0P1
!    4,SPFEP
!     WRITE(*,1112)'APATITEI',I,NX,NY,L,K,M,Precp_Ca5P3O12O3H3_conc,ACA1,XCa_conc
!    2,AH0P1,AH1P1,AH2P1,AHY1,AOH1,H2PO4_1e_apatite_dissol_flx,RHCAH1,RHCAH2
!    3,SP,SPX,ACA1**5*AH0P1**3*AOH1,SPCAH,SHCAH1,SHCAH2
!    3,H0PO4_conc,H1PO4_2e_conc,H2PO4_1e_conc,XOH_conc,XROH1_conc,XROH2_conc,XHPO4_conc,XH2PO4_conc
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
  H2CO3_aqu_conc=AMAX1(ZERO,H2CO3_aqu_conc)
  CO3_2e_conc=H2CO3_aqu_conc*DPCO3*A0/(H_1p_conc**2*A2)
  HCO3_e_conc=H2CO3_aqu_conc*DPCO2*A0/(H_1p_conc*A1)
  NH4_1p_conc=AMAX1(ZERO,NH4_1p_conc)
  NH3_aqu_conc=AMAX1(ZERO,NH3_aqu_conc)
  Al_3p_conc=AMAX1(ZERO,Al_3p_conc)
  Fe_3p_conc=AMAX1(ZERO,Fe_3p_conc)
  Ca_2p_conc=AMAX1(ZERO,Ca_2p_conc)
  Ca_2p_conc=AMIN1(Ca_2p_conc,SPCAC/(CO3_2e_conc*A2**2))
  Mg_2p_conc=AMAX1(ZERO,Mg_2p_conc)
  Na_1p_conc=AMAX1(ZERO,Na_1p_conc)
  K_1p_conc=AMAX1(ZERO,K_1p_conc)
  H1PO4_2e_conc=AMAX1(ZERO,H1PO4_2e_conc)
  H2PO4_1e_conc=AMAX1(ZERO,H2PO4_1e_conc)
!
!     PRECIPITATION-DISSOLUTION FLUXES
!
  IF(K.EQ.micpar%k_POM)THEN
    H2PO4_1e_AlPO4_eqv=SYA0P2/(Al_3p_conc*OH_1e_conc**2)
    H2PO4_1e_AlPO4_dissol_flx=AMAX1(-Precp_AlPO4_conc,TPD*(H2PO4_1e_conc-H2PO4_1e_AlPO4_eqv))
    H2PO4_1e_FePO4_eqv=SYF0P2/(Fe_3p_conc*OH_1e_conc**2._r8)
    H2PO4_1e_FePO4_dissol_flx=AMAX1(-Precp_FePO4_conc,TPD*(H2PO4_1e_conc-H2PO4_1e_FePO4_eqv))
    H2PO4_1e_CaHPO4_eqv=SYCAD2/(Ca_2p_conc*OH_1e_conc)
    H2PO4_1e_CaHPO4_dissol_flx=AMAX1(-Precp_CaHPO4_conc,TPD*(H2PO4_1e_conc-H2PO4_1e_CaHPO4_eqv))
    H2PO4_1e_apatite_eqv=(SYCAH2/(Ca_2p_conc**5*OH_1e_conc**7))**0.333
    H2PO4_1e_apatite_dissol_flx=AMAX1(-Precp_Ca5P3O12O3H3_conc,TPD*(H2PO4_1e_conc-H2PO4_1e_apatite_eqv))
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
      H2PO4_1e_to_XH2PO4_ROH2_flx=TAD*(XROH2_conc*H2PO4_1e_conc-SPH2P*XH2PO4_conc)/(XROH2_conc+SPH2P)*VLWatMicPBK
      H2PO4_to_XH2PO4_ROH_flx=TAD*(XROH1_conc*H2PO4_1e_conc-SXH2P*OH_1e_conc*XH2PO4_conc) &
        /(XROH1_conc+SXH2P*OH_1e_conc)*VLWatMicPBK
      SPH1P=SXH1P*DPH2O/DPH2P
      H1PO4_to_XHPO4_ROH_flx=TAD*(XROH1_conc*H1PO4_2e_conc-SPH1P*XHPO4_conc)/(XROH1_conc+SPH1P)*VLWatMicPBK
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
      CALX=Al_3p_conc**0.333
      CFEX=Fe_3p_conc**0.333
      CaX_conc=Ca_2p_conc**0.500_r8
      MgX_conc=Mg_2p_conc**0.500_r8
      XCAX=CEC_conc/(1.0_r8+GKC4*NH4_1p_conc/CaX_conc &
        +GKCH*H_1p_conc/CaX_conc+GKCA*CALX/CaX_conc &
        +GKCA*CFEX/CaX_conc+GKCM*MgX_conc/CaX_conc &
        +GKCN*Na_1p_conc/CaX_conc+GKCK*K_1p_conc/CaX_conc)
      XN4Q=XCAX*NH4_1p_conc*GKC4
      XHYQ=XCAX*H_1p_conc*GKCH
      XALQ=XCAX*CALX*GKCA
      XFEQ=XCAX*CFEX*GKCA
      XCAQ=XCAX*CaX_conc
      XMGQ=XCAX*MgX_conc*GKCM
      XNAQ=XCAX*Na_1p_conc*GKCN
      XKAQ=XCAX*K_1p_conc*GKCK
      XTLQ=XN4Q+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ
      IF(XTLQ.GT.ZERO)THEN
        FX=CEC_conc/XTLQ
      ELSE
        FX=0._r8
      ENDIF
      XN4Q=FX*XN4Q
      RXN4=TSL*AMIN1((XN4Q-XNH4_conc)*NH4_1p_conc/XN4Q,NH4_1p_conc)
      XNH4_conc=XNH4_conc+RXN4
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
  S0=H_1p_conc+NH3_aqu_conc+DPN4
  S1=AZMAX1(S0**2-4.0*(H_1p_conc*NH3_aqu_conc-DPN4*NH4_1p_conc))
  RNH4=TSL*(S0-SQRT(S1))
!
!     H2PO4 <-> HPO4
!
  S0=H1PO4_2e_conc+H_1p_conc+DPH2P
  S1=AZMAX1(S0**2-4.0*(H1PO4_2e_conc*H_1p_conc-DPH2P*H2PO4_1e_conc))
  H2PO4_e_to_HPO4_2e_flx=TSL*(S0-SQRT(S1))
!
!     ION FLUXES
!
  RN4S=RNH4-RXN4
  RN3S=-RNH4
  RHP1=-H1PO4_to_XHPO4_ROH_flx-H2PO4_e_to_HPO4_2e_flx
  RHP2=-H2PO4_1e_to_XH2PO4_ROH2_flx-H2PO4_to_XH2PO4_ROH_flx+H2PO4_e_to_HPO4_2e_flx &
    -H2PO4_1e_AlPO4_dissol_flx-H2PO4_1e_FePO4_dissol_flx-H2PO4_1e_CaHPO4_dissol_flx-3.0_r8*H2PO4_1e_apatite_dissol_flx
  NH4_1p_conc=NH4_1p_conc+RN4S
  NH3_aqu_conc=NH3_aqu_conc+RN3S
  H1PO4_2e_conc=H1PO4_2e_conc+RHP1
  H2PO4_1e_conc=H2PO4_1e_conc+RHP2
  end subroutine SolubilityEquilibriaNoSalt

end module InitSoluteMod
