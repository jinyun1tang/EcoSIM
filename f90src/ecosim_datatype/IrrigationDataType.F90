module IrrigationDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx => jcplxc
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: DIRRI(12)                         !change factor for irrigation, [-]

  real(r8) :: CCOU                              !subsurface irrigation  CO2 concentration	[g m-3]
  real(r8) :: CCHU                              !subsurface irrigation  CH4 concentration	[g m-3]
  real(r8) :: COXU                              !subsurface irrigation  O2 concentration	[g m-3]
  real(r8) :: CNNU                              !subsurface irrigation  N2 concentration	[g m-3]
  real(r8) :: CN2U                              !subsurface irrigation  N2O concentration	[g m-3]
  integer,target,allocatable ::  IIRRA(:,:,:)                       !start and end dates of automated irrigation, [-]
  real(r8),target,allocatable ::  RRIG(:,:,:,:)                     !irrigation application, [mm h-1]
  real(r8),target,allocatable ::  WDPTH(:,:,:)                      !depth of irrigation application, [m]
  real(r8),target,allocatable ::  PRECU(:,:)                        !underground irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  PRECI(:,:)                        !surface irrigation, [m3 d-2 h-1]
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
  real(r8),target,allocatable ::  CN4U(:,:,:)                       !subsurface irrigation  NH4 concentration, [g m-3]
  real(r8),target,allocatable ::  CN3U(:,:,:)                       !subsurface irrigation  NH3 concentration, [g m-3]
  real(r8),target,allocatable ::  CNOU(:,:,:)                       !subsurface irrigation  NO3 concentration, [g m-3]
  real(r8),target,allocatable ::  CH2PU(:,:,:)                      !subsurface irrigation H2PO4 concentrations, [g m-3]
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
  real(r8),target,allocatable ::  CH1PU(:,:,:)                      !subsurface irrigation  HPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CH3PU(:,:,:)                      !subsurface irrigation  H3PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CF1PU(:,:,:)                      !subsurface irrigation  FeHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CF2PU(:,:,:)                      !subsurface irrigation  FeH2PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC0PU(:,:,:)                      !subsurface irrigation  CaPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC1PU(:,:,:)                      !subsurface irrigation  CaHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CC2PU(:,:,:)                      !subsurface irrigation  CaH2PO4 concentration, [g m-3]
  real(r8),target,allocatable ::  CM1PU(:,:,:)                      !subsurface irrigation  MgHPO4 concentration, [g m-3]
  real(r8),target,allocatable ::  COPU(:,:,:,:)                     !subsurface irrigation  DOP concentration, [g m-3]
  real(r8),target,allocatable ::  FLU(:,:,:)                        !underground irrigation, [m3 d-2 h-1]
  real(r8),target,allocatable ::  HWFLU(:,:,:)                      !convective heat of underground irrigation, [MJ d-2 h-1]
  real(r8),target,allocatable ::  RCOFLU(:,:,:)                     !aqueous CO2 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RCHFLU(:,:,:)                     !aqueous CH4 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  ROXFLU(:,:,:)                     !aqueous O2 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RNGFLU(:,:,:)                     !aqueousN2 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RN2FLU(:,:,:)                     !aqueous N2O in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RN4FLU(:,:,:)                     !aqueous NH4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RN3FLU(:,:,:)                     !aqueous NH3 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNOFLU(:,:,:)                     !aqueous NO3 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH2PFU(:,:,:)                     !aqueous H2PO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RN4FBU(:,:,:)                     !aqueous NH4 in underground irrigation band, [g d-2 h-1]
  real(r8),target,allocatable ::  RN3FBU(:,:,:)                     !aqueous NH3 in underground irrigation band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNOFBU(:,:,:)                     !aqueous NO3 in underground irrigation band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH2BBU(:,:,:)                     !aqueous H2PO4 in underground irrigation band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH0PFU(:,:,:)                     !aqueous PO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH1PFU(:,:,:)                     !aqueous HPO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH3PFU(:,:,:)                     !aqueous H3PO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH1BBU(:,:,:)                     !aqueous HPO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RALFLU(:,:,:)                     !aqueous Al in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RFEFLU(:,:,:)                     !aqueous Fe in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RHYFLU(:,:,:)                     !aqueous H in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RCAFLU(:,:,:)                     !aqueous Ca in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RMGFLU(:,:,:)                     !aqueous Mg in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RNAFLU(:,:,:)                     !aqueous Na in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RKAFLU(:,:,:)                     !aqueous K in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  ROHFLU(:,:,:)                     !aqueous OH in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RSOFLU(:,:,:)                     !aqueous SO4 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RCLFLU(:,:,:)                     !aqueous Cl in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RC3FLU(:,:,:)                     !aqueous CO3 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RHCFLU(:,:,:)                     !aqueous HCO3 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RAL1FU(:,:,:)                     !aqueous AlOH in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RAL2FU(:,:,:)                     !aqueous AlOH2 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RAL3FU(:,:,:)                     !aqueous AlOH3 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RAL4FU(:,:,:)                     !aqueous AlOH4 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RALSFU(:,:,:)                     !aqueous AlSO4 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RFE1FU(:,:,:)                     !aqueous FeOH in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RFE2FU(:,:,:)                     !aqueous FeOH2 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RFE3FU(:,:,:)                     !aqueous FeOH3 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RFE4FU(:,:,:)                     !aqueous FeOH4 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RFESFU(:,:,:)                     !aqueous FeSO4 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RCAOFU(:,:,:)                     !aqueous CaOH in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RCACFU(:,:,:)                     !aqueous CaCO3 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RCAHFU(:,:,:)                     !aqueous CaHCO3 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RCASFU(:,:,:)                     !aqueous CaSO4 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RMGOFU(:,:,:)                     !aqueous MgOH in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RMGCFU(:,:,:)                     !aqueous MgCO3 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RMGHFU(:,:,:)                     !aqueous MgHCO3 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RMGSFU(:,:,:)                     !aqueous MgSO4 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RNACFU(:,:,:)                     !aqueous NaCO3 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RNASFU(:,:,:)                     !aqueous NaSO4 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RKASFU(:,:,:)                     !aqueous KSO4 in underground irrigation, [g d-2 h-1]
  real(r8),target,allocatable ::  RF1PFU(:,:,:)                     !aqueous FeHPO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RF2PFU(:,:,:)                     !aqueous FeH2PO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RC0PFU(:,:,:)                     !aqueous CaPO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RC1PFU(:,:,:)                     !aqueous CaHPO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RHGFLU(:,:,:)                     !aqueous H2 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RC2PFU(:,:,:)                     !aqueous CaH2PO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RM1PFU(:,:,:)                     !aqueous MgHPO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RF1BBU(:,:,:)                     !aqueous FeHPO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RF2BBU(:,:,:)                     !aqueous FeH2PO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RC0BBU(:,:,:)                     !aqueous CaPO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RC1BBU(:,:,:)                     !aqueous CaHPO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RC2BBU(:,:,:)                     !aqueous CaH2PO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RM1BBU(:,:,:)                     !aqueous MgHPO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH0BBU(:,:,:)                     !aqueous PO4 in underground irrigation non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH3BBU(:,:,:)                     !aqueous H3PO4 in underground irrigation non-band, [g d-2 h-1]


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
  allocate(PRECU(JY,JX));       PRECU=0._r8
  allocate(PRECI(JY,JX));       PRECI=0._r8
  allocate(FIRRA(JY,JX));       FIRRA=0._r8
  allocate(CIRRA(JY,JX));       CIRRA=0._r8
  allocate(DIRRA(2,JY,JX));     DIRRA=0._r8
  allocate(TDIRI(JY,JX,12));    TDIRI=0._r8
  allocate(CSTRR(JY,JX));       CSTRR=0._r8
  allocate(CCOQ(JY,JX));        CCOQ=0._r8
  allocate(CCHQ(JY,JX));        CCHQ=0._r8
  allocate(COXQ(JY,JX));        COXQ=0._r8
  allocate(CNNQ(JY,JX));        CNNQ=0._r8
  allocate(CN2Q(JY,JX));        CN2Q=0._r8
  allocate(COCU(1:jcplx,JZ,JY,JX)); COCU=0._r8
  allocate(CONU(1:jcplx,JZ,JY,JX)); CONU=0._r8
  allocate(COAU(1:jcplx,JZ,JY,JX)); COAU=0._r8
  allocate(CN4U(JZ,JY,JX));     CN4U=0._r8
  allocate(CN3U(JZ,JY,JX));     CN3U=0._r8
  allocate(CNOU(JZ,JY,JX));     CNOU=0._r8
  allocate(CH2PU(JZ,JY,JX));    CH2PU=0._r8
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
  allocate(CH1PU(JZ,JY,JX));    CH1PU=0._r8
  allocate(CH3PU(JZ,JY,JX));    CH3PU=0._r8
  allocate(CF1PU(JZ,JY,JX));    CF1PU=0._r8
  allocate(CF2PU(JZ,JY,JX));    CF2PU=0._r8
  allocate(CC0PU(JZ,JY,JX));    CC0PU=0._r8
  allocate(CC1PU(JZ,JY,JX));    CC1PU=0._r8
  allocate(CC2PU(JZ,JY,JX));    CC2PU=0._r8
  allocate(CM1PU(JZ,JY,JX));    CM1PU=0._r8
  allocate(COPU(1:jcplx,JZ,JY,JX)); COPU=0._r8
  allocate(FLU(JZ,JY,JX));      FLU=0._r8
  allocate(HWFLU(JZ,JY,JX));    HWFLU=0._r8
  allocate(RCOFLU(JZ,JY,JX));   RCOFLU=0._r8
  allocate(RCHFLU(JZ,JY,JX));   RCHFLU=0._r8
  allocate(ROXFLU(JZ,JY,JX));   ROXFLU=0._r8
  allocate(RNGFLU(JZ,JY,JX));   RNGFLU=0._r8
  allocate(RN2FLU(JZ,JY,JX));   RN2FLU=0._r8
  allocate(RN4FLU(JZ,JY,JX));   RN4FLU=0._r8
  allocate(RN3FLU(JZ,JY,JX));   RN3FLU=0._r8
  allocate(RNOFLU(JZ,JY,JX));   RNOFLU=0._r8
  allocate(RH2PFU(JZ,JY,JX));   RH2PFU=0._r8
  allocate(RN4FBU(JZ,JY,JX));   RN4FBU=0._r8
  allocate(RN3FBU(JZ,JY,JX));   RN3FBU=0._r8
  allocate(RNOFBU(JZ,JY,JX));   RNOFBU=0._r8
  allocate(RH2BBU(JZ,JY,JX));   RH2BBU=0._r8
  allocate(RH0PFU(JZ,JY,JX));   RH0PFU=0._r8
  allocate(RH1PFU(JZ,JY,JX));   RH1PFU=0._r8
  allocate(RH3PFU(JZ,JY,JX));   RH3PFU=0._r8
  allocate(RH1BBU(JZ,JY,JX));   RH1BBU=0._r8
  allocate(RALFLU(JZ,JY,JX));   RALFLU=0._r8
  allocate(RFEFLU(JZ,JY,JX));   RFEFLU=0._r8
  allocate(RHYFLU(JZ,JY,JX));   RHYFLU=0._r8
  allocate(RCAFLU(JZ,JY,JX));   RCAFLU=0._r8
  allocate(RMGFLU(JZ,JY,JX));   RMGFLU=0._r8
  allocate(RNAFLU(JZ,JY,JX));   RNAFLU=0._r8
  allocate(RKAFLU(JZ,JY,JX));   RKAFLU=0._r8
  allocate(ROHFLU(JZ,JY,JX));   ROHFLU=0._r8
  allocate(RSOFLU(JZ,JY,JX));   RSOFLU=0._r8
  allocate(RCLFLU(JZ,JY,JX));   RCLFLU=0._r8
  allocate(RC3FLU(JZ,JY,JX));   RC3FLU=0._r8
  allocate(RHCFLU(JZ,JY,JX));   RHCFLU=0._r8
  allocate(RAL1FU(JZ,JY,JX));   RAL1FU=0._r8
  allocate(RAL2FU(JZ,JY,JX));   RAL2FU=0._r8
  allocate(RAL3FU(JZ,JY,JX));   RAL3FU=0._r8
  allocate(RAL4FU(JZ,JY,JX));   RAL4FU=0._r8
  allocate(RALSFU(JZ,JY,JX));   RALSFU=0._r8
  allocate(RFE1FU(JZ,JY,JX));   RFE1FU=0._r8
  allocate(RFE2FU(JZ,JY,JX));   RFE2FU=0._r8
  allocate(RFE3FU(JZ,JY,JX));   RFE3FU=0._r8
  allocate(RFE4FU(JZ,JY,JX));   RFE4FU=0._r8
  allocate(RFESFU(JZ,JY,JX));   RFESFU=0._r8
  allocate(RCAOFU(JZ,JY,JX));   RCAOFU=0._r8
  allocate(RCACFU(JZ,JY,JX));   RCACFU=0._r8
  allocate(RCAHFU(JZ,JY,JX));   RCAHFU=0._r8
  allocate(RCASFU(JZ,JY,JX));   RCASFU=0._r8
  allocate(RMGOFU(JZ,JY,JX));   RMGOFU=0._r8
  allocate(RMGCFU(JZ,JY,JX));   RMGCFU=0._r8
  allocate(RMGHFU(JZ,JY,JX));   RMGHFU=0._r8
  allocate(RMGSFU(JZ,JY,JX));   RMGSFU=0._r8
  allocate(RNACFU(JZ,JY,JX));   RNACFU=0._r8
  allocate(RNASFU(JZ,JY,JX));   RNASFU=0._r8
  allocate(RKASFU(JZ,JY,JX));   RKASFU=0._r8
  allocate(RF1PFU(JZ,JY,JX));   RF1PFU=0._r8
  allocate(RF2PFU(JZ,JY,JX));   RF2PFU=0._r8
  allocate(RC0PFU(JZ,JY,JX));   RC0PFU=0._r8
  allocate(RC1PFU(JZ,JY,JX));   RC1PFU=0._r8
  allocate(RHGFLU(JZ,JY,JX));   RHGFLU=0._r8
  allocate(RC2PFU(JZ,JY,JX));   RC2PFU=0._r8
  allocate(RM1PFU(JZ,JY,JX));   RM1PFU=0._r8
  allocate(RF1BBU(JZ,JY,JX));   RF1BBU=0._r8
  allocate(RF2BBU(JZ,JY,JX));   RF2BBU=0._r8
  allocate(RC0BBU(JZ,JY,JX));   RC0BBU=0._r8
  allocate(RC1BBU(JZ,JY,JX));   RC1BBU=0._r8
  allocate(RC2BBU(JZ,JY,JX));   RC2BBU=0._r8
  allocate(RM1BBU(JZ,JY,JX));   RM1BBU=0._r8
  allocate(RH0BBU(JZ,JY,JX));   RH0BBU=0._r8
  allocate(RH3BBU(JZ,JY,JX));   RH3BBU=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructIrrigation

  implicit none

  if (allocated(PHQ))      deallocate(PHQ)
  if (allocated(CN4Q))     deallocate(CN4Q)
  if (allocated(CN3Q))     deallocate(CN3Q)
  if (allocated(CNOQ))     deallocate(CNOQ)
  if (allocated(CPOQ))     deallocate(CPOQ)
  if (allocated(CALQ))     deallocate(CALQ)
  if (allocated(CFEQ))     deallocate(CFEQ)
  if (allocated(CHYQ))     deallocate(CHYQ)
  if (allocated(CCAQ))     deallocate(CCAQ)
  if (allocated(CMGQ))     deallocate(CMGQ)
  if (allocated(CNAQ))     deallocate(CNAQ)
  if (allocated(CKAQ))     deallocate(CKAQ)
  if (allocated(COHQ))     deallocate(COHQ)
  if (allocated(CSOQ))     deallocate(CSOQ)
  if (allocated(CCLQ))     deallocate(CCLQ)
  if (allocated(CC3Q))     deallocate(CC3Q)
  if (allocated(CHCQ))     deallocate(CHCQ)
  if (allocated(CAL1Q))    deallocate(CAL1Q)
  if (allocated(CAL2Q))    deallocate(CAL2Q)
  if (allocated(CAL3Q))    deallocate(CAL3Q)
  if (allocated(CAL4Q))    deallocate(CAL4Q)
  if (allocated(CALSQ))    deallocate(CALSQ)
  if (allocated(CFE1Q))    deallocate(CFE1Q)
  if (allocated(CFE2Q))    deallocate(CFE2Q)
  if (allocated(CFE3Q))    deallocate(CFE3Q)
  if (allocated(CFE4Q))    deallocate(CFE4Q)
  if (allocated(CFESQ))    deallocate(CFESQ)
  if (allocated(CCAOQ))    deallocate(CCAOQ)
  if (allocated(CCACQ))    deallocate(CCACQ)
  if (allocated(CCAHQ))    deallocate(CCAHQ)
  if (allocated(CCASQ))    deallocate(CCASQ)
  if (allocated(CMGOQ))    deallocate(CMGOQ)
  if (allocated(CMGCQ))    deallocate(CMGCQ)
  if (allocated(CMGHQ))    deallocate(CMGHQ)
  if (allocated(CMGSQ))    deallocate(CMGSQ)
  if (allocated(CNACQ))    deallocate(CNACQ)
  if (allocated(CNASQ))    deallocate(CNASQ)
  if (allocated(CKASQ))    deallocate(CKASQ)
  if (allocated(CH0PQ))    deallocate(CH0PQ)
  if (allocated(CH1PQ))    deallocate(CH1PQ)
  if (allocated(CH3PQ))    deallocate(CH3PQ)
  if (allocated(CF1PQ))    deallocate(CF1PQ)
  if (allocated(CF2PQ))    deallocate(CF2PQ)
  if (allocated(CC0PQ))    deallocate(CC0PQ)
  if (allocated(CC1PQ))    deallocate(CC1PQ)
  if (allocated(CC2PQ))    deallocate(CC2PQ)
  if (allocated(CM1PQ))    deallocate(CM1PQ)
  if (allocated(CSTRQ))    deallocate(CSTRQ)
  if (allocated(RRIG))     deallocate(RRIG)
  if (allocated(WDPTH))    deallocate(WDPTH)

  if (allocated(IIRRA))    deallocate(IIRRA)
  if (allocated(PRECU))    deallocate(PRECU)
  if (allocated(PRECI))    deallocate(PRECI)
  if (allocated(FIRRA))    deallocate(FIRRA)
  if (allocated(CIRRA))    deallocate(CIRRA)
  if (allocated(DIRRA))    deallocate(DIRRA)
  if (allocated(TDIRI))    deallocate(TDIRI)

  if (allocated(CSTRR))    deallocate(CSTRR)
  if (allocated(CCOQ))     deallocate(CCOQ)
  if (allocated(CCHQ))     deallocate(CCHQ)
  if (allocated(COXQ))     deallocate(COXQ)
  if (allocated(CNNQ))     deallocate(CNNQ)
  if (allocated(CN2Q))     deallocate(CN2Q)
  if (allocated(COCU))     deallocate(COCU)
  if (allocated(CONU))     deallocate(CONU)
  if (allocated(COAU))     deallocate(COAU)
  if (allocated(CN4U))     deallocate(CN4U)
  if (allocated(CN3U))     deallocate(CN3U)
  if (allocated(CNOU))     deallocate(CNOU)
  if (allocated(CH2PU))    deallocate(CH2PU)
  if (allocated(CNZU))     deallocate(CNZU)
  if (allocated(CALU))     deallocate(CALU)
  if (allocated(CFEU))     deallocate(CFEU)
  if (allocated(CHYU))     deallocate(CHYU)
  if (allocated(CCAU))     deallocate(CCAU)
  if (allocated(CMGU))     deallocate(CMGU)
  if (allocated(CNAU))     deallocate(CNAU)
  if (allocated(CKAU))     deallocate(CKAU)
  if (allocated(COHU))     deallocate(COHU)
  if (allocated(CSOU))     deallocate(CSOU)
  if (allocated(CCLU))     deallocate(CCLU)
  if (allocated(CC3U))     deallocate(CC3U)
  if (allocated(CHCU))     deallocate(CHCU)
  if (allocated(CAL1U))    deallocate(CAL1U)
  if (allocated(CAL2U))    deallocate(CAL2U)
  if (allocated(CAL3U))    deallocate(CAL3U)
  if (allocated(CAL4U))    deallocate(CAL4U)
  if (allocated(CALSU))    deallocate(CALSU)
  if (allocated(CFE1U))    deallocate(CFE1U)
  if (allocated(CFE2U))    deallocate(CFE2U)
  if (allocated(CFE3U))    deallocate(CFE3U)
  if (allocated(CFE4U))    deallocate(CFE4U)
  if (allocated(CFESU))    deallocate(CFESU)
  if (allocated(CCAOU))    deallocate(CCAOU)
  if (allocated(CCACU))    deallocate(CCACU)
  if (allocated(CCAHU))    deallocate(CCAHU)
  if (allocated(CCASU))    deallocate(CCASU)
  if (allocated(CMGOU))    deallocate(CMGOU)
  if (allocated(CMGCU))    deallocate(CMGCU)
  if (allocated(CMGHU))    deallocate(CMGHU)
  if (allocated(CMGSU))    deallocate(CMGSU)
  if (allocated(CNACU))    deallocate(CNACU)
  if (allocated(CNASU))    deallocate(CNASU)
  if (allocated(CKASU))    deallocate(CKASU)
  if (allocated(CH0PU))    deallocate(CH0PU)
  if (allocated(CH1PU))    deallocate(CH1PU)
  if (allocated(CH3PU))    deallocate(CH3PU)
  if (allocated(CF1PU))    deallocate(CF1PU)
  if (allocated(CF2PU))    deallocate(CF2PU)
  if (allocated(CC0PU))    deallocate(CC0PU)
  if (allocated(CC1PU))    deallocate(CC1PU)
  if (allocated(CC2PU))    deallocate(CC2PU)
  if (allocated(CM1PU))    deallocate(CM1PU)
  if (allocated(COPU))     deallocate(COPU)
  if (allocated(FLU))      deallocate(FLU)
  if (allocated(HWFLU))    deallocate(HWFLU)
  if (allocated(RCOFLU))   deallocate(RCOFLU)
  if (allocated(RCHFLU))   deallocate(RCHFLU)
  if (allocated(ROXFLU))   deallocate(ROXFLU)
  if (allocated(RNGFLU))   deallocate(RNGFLU)
  if (allocated(RN2FLU))   deallocate(RN2FLU)
  if (allocated(RN4FLU))   deallocate(RN4FLU)
  if (allocated(RN3FLU))   deallocate(RN3FLU)
  if (allocated(RNOFLU))   deallocate(RNOFLU)
  if (allocated(RH2PFU))   deallocate(RH2PFU)
  if (allocated(RN4FBU))   deallocate(RN4FBU)
  if (allocated(RN3FBU))   deallocate(RN3FBU)
  if (allocated(RNOFBU))   deallocate(RNOFBU)
  if (allocated(RH2BBU))   deallocate(RH2BBU)
  if (allocated(RH0PFU))   deallocate(RH0PFU)
  if (allocated(RH1PFU))   deallocate(RH1PFU)
  if (allocated(RH3PFU))   deallocate(RH3PFU)
  if (allocated(RH1BBU))   deallocate(RH1BBU)
  if (allocated(RALFLU))   deallocate(RALFLU)
  if (allocated(RFEFLU))   deallocate(RFEFLU)
  if (allocated(RHYFLU))   deallocate(RHYFLU)
  if (allocated(RCAFLU))   deallocate(RCAFLU)
  if (allocated(RMGFLU))   deallocate(RMGFLU)
  if (allocated(RNAFLU))   deallocate(RNAFLU)
  if (allocated(RKAFLU))   deallocate(RKAFLU)
  if (allocated(ROHFLU))   deallocate(ROHFLU)
  if (allocated(RSOFLU))   deallocate(RSOFLU)
  if (allocated(RCLFLU))   deallocate(RCLFLU)
  if (allocated(RC3FLU))   deallocate(RC3FLU)
  if (allocated(RHCFLU))   deallocate(RHCFLU)
  if (allocated(RAL1FU))   deallocate(RAL1FU)
  if (allocated(RAL2FU))   deallocate(RAL2FU)
  if (allocated(RAL3FU))   deallocate(RAL3FU)
  if (allocated(RAL4FU))   deallocate(RAL4FU)
  if (allocated(RALSFU))   deallocate(RALSFU)
  if (allocated(RFE1FU))   deallocate(RFE1FU)
  if (allocated(RFE2FU))   deallocate(RFE2FU)
  if (allocated(RFE3FU))   deallocate(RFE3FU)
  if (allocated(RFE4FU))   deallocate(RFE4FU)
  if (allocated(RFESFU))   deallocate(RFESFU)
  if (allocated(RCAOFU))   deallocate(RCAOFU)
  if (allocated(RCACFU))   deallocate(RCACFU)
  if (allocated(RCAHFU))   deallocate(RCAHFU)
  if (allocated(RCASFU))   deallocate(RCASFU)
  if (allocated(RMGOFU))   deallocate(RMGOFU)
  if (allocated(RMGCFU))   deallocate(RMGCFU)
  if (allocated(RMGHFU))   deallocate(RMGHFU)
  if (allocated(RMGSFU))   deallocate(RMGSFU)
  if (allocated(RNACFU))   deallocate(RNACFU)
  if (allocated(RNASFU))   deallocate(RNASFU)
  if (allocated(RKASFU))   deallocate(RKASFU)
  if (allocated(RF1PFU))   deallocate(RF1PFU)
  if (allocated(RF2PFU))   deallocate(RF2PFU)
  if (allocated(RC0PFU))   deallocate(RC0PFU)
  if (allocated(RC1PFU))   deallocate(RC1PFU)
  if (allocated(RHGFLU))   deallocate(RHGFLU)
  if (allocated(RC2PFU))   deallocate(RC2PFU)
  if (allocated(RM1PFU))   deallocate(RM1PFU)
  if (allocated(RF1BBU))   deallocate(RF1BBU)
  if (allocated(RF2BBU))   deallocate(RF2BBU)
  if (allocated(RC0BBU))   deallocate(RC0BBU)
  if (allocated(RC1BBU))   deallocate(RC1BBU)
  if (allocated(RC2BBU))   deallocate(RC2BBU)
  if (allocated(RM1BBU))   deallocate(RM1BBU)
  if (allocated(RH0BBU))   deallocate(RH0BBU)
  if (allocated(RH3BBU))   deallocate(RH3BBU)

  end subroutine DestructIrrigation

end module IrrigationDataType
