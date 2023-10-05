module IrrigationDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx => jcplxc
  use TracerIDMod
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8) :: DIRRI(12)                         !change factor for irrigation, [-]

  real(r8) :: CCOU                              !subsurface irrigation  CO2 concentration	[g m-3]
  real(r8) :: CCHU                              !subsurface irrigation  CH4 concentration	[g m-3]
  real(r8) :: COXU                              !subsurface irrigation  O2 concentration	[g m-3]
  real(r8) :: CNNU                              !subsurface irrigation  N2 concentration	[g m-3]
  real(r8) :: CN2U                              !subsurface irrigation  N2O concentration	[g m-3]
  integer ,target,allocatable ::  IIRRA(:,:,:)                       !start and end dates of automated irrigation, [-]
  real(r8),target,allocatable ::  RRIG(:,:,:,:)                      !irrigation application, [mm h-1]
  real(r8),target,allocatable ::  WDPTH(:,:,:)                       !depth of irrigation application, [m]
  real(r8),target,allocatable ::  IrrigSubsurf(:,:)                  !underground irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  IrrigSurface(:,:)                        !surface irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  FIRRA(:,:)                        !fraction of FC - WP below which automatic irrigation applied, [-]
  real(r8),target,allocatable ::  CIRRA(:,:)                        !fraction of FC - WP to which automatic irrigation applied, [-]
  real(r8),target,allocatable ::  DIRRA(:,:,:)                      !depth to which automatic irrigation applied, [m]
  real(r8),target,allocatable ::  TDIRI(:,:,:)                      !accumulated change  for irrigation, [-]
  real(r8),target,allocatable ::  PHQ(:,:,:)                        !surface irrigation  pH, [-]
  real(r8),target,allocatable ::  CN4Q(:,:,:)                       !surface irrigation  NH4 concentration, [g m-3]
  real(r8),target,allocatable ::  CN3Q(:,:,:)                       !surface irrigation  NH3 concentration, [g m-3]
  real(r8),target,allocatable ::  CNOQ(:,:,:)                       !surface irrigation  NO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CPOQ(:,:,:)                       !surface irrigation  H2PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CALQ(:,:,:)                       !surface irrigation  Al concentration, [g m-3]
  real(r8),target,allocatable ::  CFEQ(:,:,:)                       !surface irrigation  Fe concentration, [g m-3]
  real(r8),target,allocatable ::  CHYQ(:,:,:)                       !surface irrigation  H concentration, [g m-3]
  real(r8),target,allocatable ::  CCAQ(:,:,:)                       !surface irrigation  Ca concentration, [g m-3]
  real(r8),target,allocatable ::  CMGQ(:,:,:)                       !surface irrigation  Mg concentration, [g m-3]
  real(r8),target,allocatable ::  CNAQ(:,:,:)                       !surface irrigation  Na concentration, [g m-3]
  real(r8),target,allocatable ::  CKAQ(:,:,:)                       !surface irrigation  K concentration, [g m-3]
  real(r8),target,allocatable ::  COHQ(:,:,:)                       !surface irrigation  OH concentration, [g m-3]
  real(r8),target,allocatable ::  CSOQ(:,:,:)                       !surface irrigation  SO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CCLQ(:,:,:)                       !surface irrigation  Cl concentration, [g m-3]
  real(r8),target,allocatable ::  CC3Q(:,:,:)                       !surface irrigation  CO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CHCQ(:,:,:)                       !surface irrigation  HCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CAL1Q(:,:,:)                      !surface irrigation  AlOH concentration, [g m-3]
  real(r8),target,allocatable ::  CAL2Q(:,:,:)                      !surface irrigation  AlOH2 concentration, [g m-3]
  real(r8),target,allocatable ::  CAL3Q(:,:,:)                      !surface irrigation  AlOH3 concentration, [g m-3]
  real(r8),target,allocatable ::  CAL4Q(:,:,:)                      !surface irrigation  AlOH4 concentration, [g m-3]
  real(r8),target,allocatable ::  CALSQ(:,:,:)                      !surface irrigation  AlSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CFE1Q(:,:,:)                      !surface irrigation  FeOH concentration, [g m-3]
  real(r8),target,allocatable ::  CFE2Q(:,:,:)                      !surface irrigation  FeOH2 concentration, [g m-3]
  real(r8),target,allocatable ::  CFE3Q(:,:,:)                      !surface irrigation  FeOH3 concentration, [g m-3]
  real(r8),target,allocatable ::  CFE4Q(:,:,:)                      !surface irrigation  FeOH4 concentration, [g m-3]
  real(r8),target,allocatable ::  CFESQ(:,:,:)                      !surface irrigation  FeSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CCAOQ(:,:,:)                      !surface irrigation  CaOH concentration, [g m-3]
  real(r8),target,allocatable ::  CCACQ(:,:,:)                      !surface irrigation  CaCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CCAHQ(:,:,:)                      !surface irrigation  CaHCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CCASQ(:,:,:)                      !surface irrigation  CaSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CMGOQ(:,:,:)                      !surface irrigation  MgOH concentration, [g m-3]
  real(r8),target,allocatable ::  CMGCQ(:,:,:)                      !surface irrigation  MgCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CMGHQ(:,:,:)                      !surface irrigation  MgHCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CMGSQ(:,:,:)                      !surface irrigation  MgSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CNACQ(:,:,:)                      !surface irrigation  NaCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CNASQ(:,:,:)                      !surface irrigation  NaSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CKASQ(:,:,:)                      !surface irrigation  K concentration, [g m-3]
  real(r8),target,allocatable ::  CH0PQ(:,:,:)                      !surface irrigation  PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CH1PQ(:,:,:)                      !surface irrigation  HPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CH3PQ(:,:,:)                      !surface irrigation  H3PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CF1PQ(:,:,:)                      !surface irrigation  FeHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CF2PQ(:,:,:)                      !surface irrigation  FeH2PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC0PQ(:,:,:)                      !surface irrigation  CaPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC1PQ(:,:,:)                      !surface irrigation  CaHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC2PQ(:,:,:)                      !surface irrigation  CaH2PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CM1PQ(:,:,:)                      !surface irrigation  MgHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CSTRQ(:,:,:)                      !surface irrigation ion strength, [g m-3]
  real(r8),target,allocatable ::  CSTRR(:,:)                        !surface irrigation ion strength, [g m-3]
  real(r8),target,allocatable ::  CCOQ(:,:)                         !surface irrigation  CO2 concentration, [g m-3]
  real(r8),target,allocatable ::  CCHQ(:,:)                         !surface irrigation  CH4 concentration, [g m-3]
  real(r8),target,allocatable ::  COXQ(:,:)                         !surface irrigation  O2 concentration, [g m-3]
  real(r8),target,allocatable ::  CNNQ(:,:)                         !surface irrigation  N2 concentration, [g m-3]
  real(r8),target,allocatable ::  CN2Q(:,:)                         !surface irrigation  N2O concentration, [g m-3]
  real(r8),target,allocatable ::  COCU(:,:,:,:)                     !subsurface irrigation  DOC concentration, [g m-3]
  real(r8),target,allocatable ::  CONU(:,:,:,:)                     !subsurface irrigation  DON concentration, [g m-3]
  real(r8),target,allocatable ::  COAU(:,:,:,:)                     !subsurface irrigation  acetate concentration, [g m-3]
  real(r8),target,allocatable ::  trcn_irrig(:,:,:,:)                !subsurface irrigation  nutrient concentration, [g m-3]
  real(r8),target,allocatable ::  CNZU(:,:,:)                       !subsurface irrigation  Al concentration, [g m-3]
  real(r8),target,allocatable ::  CALU(:,:,:)                       !subsurface irrigation  Fe concentration, [g m-3]
  real(r8),target,allocatable ::  CFEU(:,:,:)                       !subsurface irrigation  H concentration, [g m-3]
  real(r8),target,allocatable ::  CHYU(:,:,:)                       !subsurface irrigation  Ca concentration, [g m-3]
  real(r8),target,allocatable ::  CCAU(:,:,:)                       !subsurface irrigation  Mg concentration, [g m-3]
  real(r8),target,allocatable ::  CMGU(:,:,:)                       !subsurface irrigation  Na concentration, [g m-3]
  real(r8),target,allocatable ::  CNAU(:,:,:)                       !subsurface irrigation  K concentration, [g m-3]
  real(r8),target,allocatable ::  CKAU(:,:,:)                       !subsurface irrigation  OH concentration, [g m-3]
  real(r8),target,allocatable ::  COHU(:,:,:)                       !subsurface irrigation  SO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CSOU(:,:,:)                       !subsurface irrigation  Cl concentration, [g m-3]
  real(r8),target,allocatable ::  CCLU(:,:,:)                       !subsurface irrigation  CO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CC3U(:,:,:)                       !subsurface irrigation  HCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CHCU(:,:,:)                       !subsurface irrigation  CH4 concentration, [g m-3]
  real(r8),target,allocatable ::  CAL1U(:,:,:)                      !subsurface irrigation  AlOH concentration, [g m-3]
  real(r8),target,allocatable ::  CAL2U(:,:,:)                      !subsurface irrigation  AlOH2 concentration, [g m-3]
  real(r8),target,allocatable ::  CAL3U(:,:,:)                      !subsurface irrigation  AlOH3 concentration, [g m-3]
  real(r8),target,allocatable ::  CAL4U(:,:,:)                      !subsurface irrigation  AlOH4 concentration, [g m-3]
  real(r8),target,allocatable ::  CALSU(:,:,:)                      !subsurface irrigation  AlSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CFE1U(:,:,:)                      !subsurface irrigation  FeOH concentration, [g m-3]
  real(r8),target,allocatable ::  CFE2U(:,:,:)                      !subsurface irrigation  FeOH2 concentration, [g m-3]
  real(r8),target,allocatable ::  CFE3U(:,:,:)                      !subsurface irrigation  FeOH3 concentration, [g m-3]
  real(r8),target,allocatable ::  CFE4U(:,:,:)                      !subsurface irrigation  FeOH4 concentration, [g m-3]
  real(r8),target,allocatable ::  CFESU(:,:,:)                      !subsurface irrigation  FeSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CCAOU(:,:,:)                      !subsurface irrigation  CaOH concentration, [g m-3]
  real(r8),target,allocatable ::  CCACU(:,:,:)                      !subsurface irrigation  CaCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CCAHU(:,:,:)                      !subsurface irrigation  CaHCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CCASU(:,:,:)                      !subsurface irrigation  CaSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CMGOU(:,:,:)                      !subsurface irrigation  MgOH concentration, [g m-3]
  real(r8),target,allocatable ::  CMGCU(:,:,:)                      !subsurface irrigation  MgCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CMGHU(:,:,:)                      !subsurface irrigation  MgHCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CMGSU(:,:,:)                      !subsurface irrigation  MgSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CNACU(:,:,:)                      !subsurface irrigation  NaCO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CNASU(:,:,:)                      !subsurface irrigation  NaSO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CKASU(:,:,:)                      !subsurface irrigation  K concentration, [g m-3]
  real(r8),target,allocatable ::  CH0PU(:,:,:)                      !subsurface irrigation  PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CH3PU(:,:,:)                      !subsurface irrigation  H3PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CF1PU(:,:,:)                      !subsurface irrigation  FeHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CF2PU(:,:,:)                      !subsurface irrigation  FeH2PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC0PU(:,:,:)                      !subsurface irrigation  CaPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC1PU(:,:,:)                      !subsurface irrigation  CaHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC2PU(:,:,:)                      !subsurface irrigation  CaH2PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CM1PU(:,:,:)                      !subsurface irrigation  MgHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  COPU(:,:,:,:)                     !subsurface irrigation  DOP concentration, [g m-3]
  real(r8),target,allocatable ::  FWatIrrigate2MicP(:,:,:)                        !underground irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HeatIrrigation(:,:,:)                      !convective heat of underground irrigation, [MJ d-2 h-1]
  real(r8),target,allocatable ::  trcs_RFLU(:,:,:,:)                     !aqueous non-salt solutes in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  trcSalt_RFLU(:,:,:,:)                     !aqueous PO4 in underground irrigation non-band, [g d-2 h-1]
  private :: InitAllocate
  contains

  subroutine InitIrrigation

  implicit none

  call InitAllocate
  end subroutine InitIrrigation

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none

  allocate(PHQ(366,JY,JX));     PHQ=0._r8
  allocate(CN4Q(366,JY,JX));    CN4Q=0._r8
  allocate(CN3Q(366,JY,JX));    CN3Q=0._r8
  allocate(CNOQ(366,JY,JX));    CNOQ=0._r8
  allocate(CPOQ(366,JY,JX));    CPOQ=0._r8
  allocate(CALQ(366,JY,JX));    CALQ=0._r8
  allocate(CFEQ(366,JY,JX));    CFEQ=0._r8
  allocate(CHYQ(366,JY,JX));    CHYQ=0._r8
  allocate(CCAQ(366,JY,JX));    CCAQ=0._r8
  allocate(CMGQ(366,JY,JX));    CMGQ=0._r8
  allocate(CNAQ(366,JY,JX));    CNAQ=0._r8
  allocate(CKAQ(366,JY,JX));    CKAQ=0._r8
  allocate(COHQ(366,JY,JX));    COHQ=0._r8
  allocate(CSOQ(366,JY,JX));    CSOQ=0._r8
  allocate(CCLQ(366,JY,JX));    CCLQ=0._r8
  allocate(CC3Q(366,JY,JX));    CC3Q=0._r8
  allocate(CHCQ(366,JY,JX));    CHCQ=0._r8
  allocate(CAL1Q(366,JY,JX));   CAL1Q=0._r8
  allocate(CAL2Q(366,JY,JX));   CAL2Q=0._r8
  allocate(CAL3Q(366,JY,JX));   CAL3Q=0._r8
  allocate(CAL4Q(366,JY,JX));   CAL4Q=0._r8
  allocate(CALSQ(366,JY,JX));   CALSQ=0._r8
  allocate(CFE1Q(366,JY,JX));   CFE1Q=0._r8
  allocate(CFE2Q(366,JY,JX));   CFE2Q=0._r8
  allocate(CFE3Q(366,JY,JX));   CFE3Q=0._r8
  allocate(CFE4Q(366,JY,JX));   CFE4Q=0._r8
  allocate(CFESQ(366,JY,JX));   CFESQ=0._r8
  allocate(CCAOQ(366,JY,JX));   CCAOQ=0._r8
  allocate(CCACQ(366,JY,JX));   CCACQ=0._r8
  allocate(CCAHQ(366,JY,JX));   CCAHQ=0._r8
  allocate(CCASQ(366,JY,JX));   CCASQ=0._r8
  allocate(CMGOQ(366,JY,JX));   CMGOQ=0._r8
  allocate(CMGCQ(366,JY,JX));   CMGCQ=0._r8
  allocate(CMGHQ(366,JY,JX));   CMGHQ=0._r8
  allocate(CMGSQ(366,JY,JX));   CMGSQ=0._r8
  allocate(CNACQ(366,JY,JX));   CNACQ=0._r8
  allocate(CNASQ(366,JY,JX));   CNASQ=0._r8
  allocate(CKASQ(366,JY,JX));   CKASQ=0._r8
  allocate(CH0PQ(366,JY,JX));   CH0PQ=0._r8
  allocate(CH1PQ(366,JY,JX));   CH1PQ=0._r8
  allocate(CH3PQ(366,JY,JX));   CH3PQ=0._r8
  allocate(CF1PQ(366,JY,JX));   CF1PQ=0._r8
  allocate(CF2PQ(366,JY,JX));   CF2PQ=0._r8
  allocate(CC0PQ(366,JY,JX));   CC0PQ=0._r8
  allocate(CC1PQ(366,JY,JX));   CC1PQ=0._r8
  allocate(CC2PQ(366,JY,JX));   CC2PQ=0._r8
  allocate(CM1PQ(366,JY,JX));   CM1PQ=0._r8
  allocate(CSTRQ(366,JY,JX));   CSTRQ=0._r8
  allocate(RRIG(24,366,JY,JX)); RRIG=0._r8
  allocate(WDPTH(366,JY,JX));   WDPTH=0._r8

  allocate(IIRRA(4,JY,JX));     IIRRA=0
  allocate(IrrigSubsurf(JY,JX));       IrrigSubsurf=0._r8
  allocate(IrrigSurface(JY,JX));       IrrigSurface=0._r8
  allocate(FIRRA(JY,JX));       FIRRA=0._r8
  allocate(CIRRA(JY,JX));       CIRRA=0._r8
  allocate(DIRRA(2,JY,JX));     DIRRA=0._r8
  allocate(TDIRI(12,JY,JX));    TDIRI=0._r8
  allocate(CSTRR(JY,JX));       CSTRR=0._r8
  allocate(CCOQ(JY,JX));        CCOQ=0._r8
  allocate(CCHQ(JY,JX));        CCHQ=0._r8
  allocate(COXQ(JY,JX));        COXQ=0._r8
  allocate(CNNQ(JY,JX));        CNNQ=0._r8
  allocate(CN2Q(JY,JX));        CN2Q=0._r8
  allocate(COCU(1:jcplx,JZ,JY,JX)); COCU=0._r8
  allocate(CONU(1:jcplx,JZ,JY,JX)); CONU=0._r8
  allocate(COAU(1:jcplx,JZ,JY,JX)); COAU=0._r8
  allocate(trcn_irrig(ids_nuts_beg:ids_nuts_end,JZ,JY,JX)); trcn_irrig=0._r8
  allocate(CNZU(JZ,JY,JX));     CNZU=0._r8
  allocate(CALU(JZ,JY,JX));     CALU=0._r8
  allocate(CFEU(JZ,JY,JX));     CFEU=0._r8
  allocate(CHYU(JZ,JY,JX));     CHYU=0._r8
  allocate(CCAU(JZ,JY,JX));     CCAU=0._r8
  allocate(CMGU(JZ,JY,JX));     CMGU=0._r8
  allocate(CNAU(JZ,JY,JX));     CNAU=0._r8
  allocate(CKAU(JZ,JY,JX));     CKAU=0._r8
  allocate(COHU(JZ,JY,JX));     COHU=0._r8
  allocate(CSOU(JZ,JY,JX));     CSOU=0._r8
  allocate(CCLU(JZ,JY,JX));     CCLU=0._r8
  allocate(CC3U(JZ,JY,JX));     CC3U=0._r8
  allocate(CHCU(JZ,JY,JX));     CHCU=0._r8
  allocate(CAL1U(JZ,JY,JX));    CAL1U=0._r8
  allocate(CAL2U(JZ,JY,JX));    CAL2U=0._r8
  allocate(CAL3U(JZ,JY,JX));    CAL3U=0._r8
  allocate(CAL4U(JZ,JY,JX));    CAL4U=0._r8
  allocate(CALSU(JZ,JY,JX));    CALSU=0._r8
  allocate(CFE1U(JZ,JY,JX));    CFE1U=0._r8
  allocate(CFE2U(JZ,JY,JX));    CFE2U=0._r8
  allocate(CFE3U(JZ,JY,JX));    CFE3U=0._r8
  allocate(CFE4U(JZ,JY,JX));    CFE4U=0._r8
  allocate(CFESU(JZ,JY,JX));    CFESU=0._r8
  allocate(CCAOU(JZ,JY,JX));    CCAOU=0._r8
  allocate(CCACU(JZ,JY,JX));    CCACU=0._r8
  allocate(CCAHU(JZ,JY,JX));    CCAHU=0._r8
  allocate(CCASU(JZ,JY,JX));    CCASU=0._r8
  allocate(CMGOU(JZ,JY,JX));    CMGOU=0._r8
  allocate(CMGCU(JZ,JY,JX));    CMGCU=0._r8
  allocate(CMGHU(JZ,JY,JX));    CMGHU=0._r8
  allocate(CMGSU(JZ,JY,JX));    CMGSU=0._r8
  allocate(CNACU(JZ,JY,JX));    CNACU=0._r8
  allocate(CNASU(JZ,JY,JX));    CNASU=0._r8
  allocate(CKASU(JZ,JY,JX));    CKASU=0._r8
  allocate(CH0PU(JZ,JY,JX));    CH0PU=0._r8
  allocate(CH3PU(JZ,JY,JX));    CH3PU=0._r8
  allocate(CF1PU(JZ,JY,JX));    CF1PU=0._r8
  allocate(CF2PU(JZ,JY,JX));    CF2PU=0._r8
  allocate(CC0PU(JZ,JY,JX));    CC0PU=0._r8
  allocate(CC1PU(JZ,JY,JX));    CC1PU=0._r8
  allocate(CC2PU(JZ,JY,JX));    CC2PU=0._r8
  allocate(CM1PU(JZ,JY,JX));    CM1PU=0._r8
  allocate(COPU(1:jcplx,JZ,JY,JX)); COPU=0._r8
  allocate(FWatIrrigate2MicP(JZ,JY,JX));      FWatIrrigate2MicP=0._r8
  allocate(HeatIrrigation(JZ,JY,JX));    HeatIrrigation=0._r8
  allocate(trcs_RFLU(ids_beg:ids_end,JZ,JY,JX));   trcs_RFLU=0._r8
  allocate(trcSalt_RFLU(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_RFLU=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructIrrigation
  use abortutils, only : destroy
  implicit none

  call destroy(trcn_irrig)
  call destroy(trcs_RFLU)
  call destroy(PHQ)
  call destroy(CN4Q)
  call destroy(CN3Q)
  call destroy(CNOQ)
  call destroy(CPOQ)
  call destroy(CALQ)
  call destroy(CFEQ)
  call destroy(CHYQ)
  call destroy(CCAQ)
  call destroy(CMGQ)
  call destroy(CNAQ)
  call destroy(CKAQ)
  call destroy(COHQ)
  call destroy(CSOQ)
  call destroy(CCLQ)
  call destroy(CC3Q)
  call destroy(CHCQ)
  call destroy(CAL1Q)
  call destroy(CAL2Q)
  call destroy(CAL3Q)
  call destroy(CAL4Q)
  call destroy(CALSQ)
  call destroy(CFE1Q)
  call destroy(CFE2Q)
  call destroy(CFE3Q)
  call destroy(CFE4Q)
  call destroy(CFESQ)
  call destroy(CCAOQ)
  call destroy(CCACQ)
  call destroy(CCAHQ)
  call destroy(CCASQ)
  call destroy(CMGOQ)
  call destroy(CMGCQ)
  call destroy(CMGHQ)
  call destroy(CMGSQ)
  call destroy(CNACQ)
  call destroy(CNASQ)
  call destroy(CKASQ)
  call destroy(CH0PQ)
  call destroy(CH1PQ)
  call destroy(CH3PQ)
  call destroy(CF1PQ)
  call destroy(CF2PQ)
  call destroy(CC0PQ)
  call destroy(CC1PQ)
  call destroy(CC2PQ)
  call destroy(CM1PQ)
  call destroy(CSTRQ)
  call destroy(RRIG)
  call destroy(WDPTH)
  call destroy(IIRRA)
  call destroy(IrrigSubsurf)
  call destroy(IrrigSurface)
  call destroy(FIRRA)
  call destroy(CIRRA)
  call destroy(DIRRA)
  call destroy(TDIRI)
  call destroy(CSTRR)
  call destroy(CCOQ)
  call destroy(CCHQ)
  call destroy(COXQ)
  call destroy(CNNQ)
  call destroy(CN2Q)
  call destroy(COCU)
  call destroy(CONU)
  call destroy(COAU)
  call destroy(CNZU)
  call destroy(CALU)
  call destroy(CFEU)
  call destroy(CHYU)
  call destroy(CCAU)
  call destroy(CMGU)
  call destroy(CNAU)
  call destroy(CKAU)
  call destroy(COHU)
  call destroy(CSOU)
  call destroy(CCLU)
  call destroy(CC3U)
  call destroy(CHCU)
  call destroy(CAL1U)
  call destroy(CAL2U)
  call destroy(CAL3U)
  call destroy(CAL4U)
  call destroy(CALSU)
  call destroy(CFE1U)
  call destroy(CFE2U)
  call destroy(CFE3U)
  call destroy(CFE4U)
  call destroy(CFESU)
  call destroy(CCAOU)
  call destroy(CCACU)
  call destroy(CCAHU)
  call destroy(CCASU)
  call destroy(CMGOU)
  call destroy(CMGCU)
  call destroy(CMGHU)
  call destroy(CMGSU)
  call destroy(CNACU)
  call destroy(CNASU)
  call destroy(CKASU)
  call destroy(CH0PU)
  call destroy(CH3PU)
  call destroy(CF1PU)
  call destroy(CF2PU)
  call destroy(CC0PU)
  call destroy(CC1PU)
  call destroy(CC2PU)
  call destroy(CM1PU)
  call destroy(COPU)
  call destroy(FWatIrrigate2MicP)
  call destroy(HeatIrrigation)

  end subroutine DestructIrrigation

end module IrrigationDataType
