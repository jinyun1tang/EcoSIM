module FertilizerDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable :: ZNH4FA(:,:,:)
  real(r8),allocatable :: ZNH3FA(:,:,:)
  real(r8),allocatable :: ZNHUFA(:,:,:)
  real(r8),allocatable :: ZNO3FA(:,:,:)
  real(r8),allocatable :: ZNH4FB(:,:,:)
  real(r8),allocatable :: ZNH3FB(:,:,:)
  real(r8),allocatable :: ZNHUFB(:,:,:)
  real(r8),allocatable :: ZNO3FB(:,:,:)

  real(r8) :: FERT(20,366,JY,JX)                !fertilizer application, [g m-2]
  real(r8) :: FDPTH(366,JY,JX)                  !depth of fertilizer application, [m]
  real(r8) :: ROWN(JY,JX)                       !row spacing of NH4 fertilizer band, [m]
  real(r8) :: ROWO(JY,JX)                       !row spacing of NO3 fertilizer band, [m]
  real(r8) :: ROWP(JY,JX)                       !row spacing of PO4 fertilizer band, [m]
  real(r8) :: ROWI(366,JY,JX)                   !row spacing of fertilizer band, [m]

  real(r8) :: IIRRA(4,JY,JX)                    !start and end dates of automated irrigation, [-]
  real(r8) :: RRIG(24,366,JY,JX)                !irrigation application, [mm h-1]
  real(r8) :: WDPTH(366,JY,JX)                  !depth of irrigation application, [m]
  real(r8) :: DCORP(366,JY,JX)                  !soil mixing fraction with tillage, [-]
  real(r8) :: PRECU(JY,JX)                      !underground irrigation, [m3 d-2 h-1]
  real(r8) :: PRECI(JY,JX)                      !surface irrigation, [m3 d-2 h-1]
  real(r8) :: FIRRA(JY,JX)                      !fraction of FC - WP below which automatic irrigation applied, [-]
  real(r8) :: CIRRA(JY,JX)                      !fraction of FC - WP to which automatic irrigation applied, [-]
  real(r8) :: DIRRA(2,JY,JX)                    !depth to which automatic irrigation applied, [m]
  real(r8) :: PHQ(366,JY,JX)                    !surface irrigation  pH, [-]
  real(r8) :: CN4Q(366,JY,JX)                   !surface irrigation  NH4 concentration, [g m-3]
  real(r8) :: CN3Q(366,JY,JX)                   !surface irrigation  NH3 concentration, [g m-3]
  real(r8) :: CNOQ(366,JY,JX)                   !surface irrigation  NO3 concentration, [g m-3]
  real(r8) :: CPOQ(366,JY,JX)                   !surface irrigation  H2PO4 concentration, [g m-3]
  real(r8) :: CALQ(366,JY,JX)                   !surface irrigation  Al concentration, [g m-3]
  real(r8) :: CFEQ(366,JY,JX)                   !surface irrigation  Fe concentration, [g m-3]
  real(r8) :: CHYQ(366,JY,JX)                   !surface irrigation  H concentration, [g m-3]
  real(r8) :: CCAQ(366,JY,JX)                   !surface irrigation  Ca concentration, [g m-3]
  real(r8) :: CMGQ(366,JY,JX)                   !surface irrigation  Mg concentration, [g m-3]
  real(r8) :: CNAQ(366,JY,JX)                   !surface irrigation  Na concentration, [g m-3]
  real(r8) :: CKAQ(366,JY,JX)                   !surface irrigation  K concentration, [g m-3]
  real(r8) :: COHQ(366,JY,JX)                   !surface irrigation  OH concentration, [g m-3]
  real(r8) :: CSOQ(366,JY,JX)                   !surface irrigation  SO4 concentration, [g m-3]
  real(r8) :: CCLQ(366,JY,JX)                   !surface irrigation  Cl concentration, [g m-3]
  real(r8) :: CC3Q(366,JY,JX)                   !surface irrigation  CO3 concentration, [g m-3]
  real(r8) :: CHCQ(366,JY,JX)                   !surface irrigation  HCO3 concentration, [g m-3]
  real(r8) :: CAL1Q(366,JY,JX)                  !surface irrigation  AlOH concentration, [g m-3]
  real(r8) :: CAL2Q(366,JY,JX)                  !surface irrigation  AlOH2 concentration, [g m-3]
  real(r8) :: CAL3Q(366,JY,JX)                  !surface irrigation  AlOH3 concentration, [g m-3]
  real(r8) :: CAL4Q(366,JY,JX)                  !surface irrigation  AlOH4 concentration, [g m-3]
  real(r8) :: CALSQ(366,JY,JX)                  !surface irrigation  AlSO4 concentration, [g m-3]
  real(r8) :: CFE1Q(366,JY,JX)                  !surface irrigation  FeOH concentration, [g m-3]
  real(r8) :: CFE2Q(366,JY,JX)                  !surface irrigation  FeOH2 concentration, [g m-3]
  real(r8) :: CFE3Q(366,JY,JX)                  !surface irrigation  FeOH3 concentration, [g m-3]
  real(r8) :: CFE4Q(366,JY,JX)                  !surface irrigation  FeOH4 concentration, [g m-3]
  real(r8) :: CFESQ(366,JY,JX)                  !surface irrigation  FeSO4 concentration, [g m-3]
  real(r8) :: CCAOQ(366,JY,JX)                  !surface irrigation  CaOH concentration, [g m-3]
  real(r8) :: CCACQ(366,JY,JX)                  !surface irrigation  CaCO3 concentration, [g m-3]
  real(r8) :: CCAHQ(366,JY,JX)                  !surface irrigation  CaHCO3 concentration, [g m-3]
  real(r8) :: CCASQ(366,JY,JX)                  !surface irrigation  CaSO4 concentration, [g m-3]
  real(r8) :: CMGOQ(366,JY,JX)                  !surface irrigation  MgOH concentration, [g m-3]
  real(r8) :: CMGCQ(366,JY,JX)                  !surface irrigation  MgCO3 concentration, [g m-3]
  real(r8) :: CMGHQ(366,JY,JX)                  !surface irrigation  MgHCO3 concentration, [g m-3]
  real(r8) :: CMGSQ(366,JY,JX)                  !surface irrigation  MgSO4 concentration, [g m-3]
  real(r8) :: CNACQ(366,JY,JX)                  !surface irrigation  NaCO3 concentration, [g m-3]
  real(r8) :: CNASQ(366,JY,JX)                  !surface irrigation  NaSO4 concentration, [g m-3]
  real(r8) :: CKASQ(366,JY,JX)                  !surface irrigation  K concentration, [g m-3]
  real(r8) :: CH0PQ(366,JY,JX)                  !surface irrigation  PO4 concentration, [g m-3]
  real(r8) :: CH1PQ(366,JY,JX)                  !surface irrigation  HPO4 concentration, [g m-3]
  real(r8) :: CH3PQ(366,JY,JX)                  !surface irrigation  H3PO4 concentration, [g m-3]
  real(r8) :: CF1PQ(366,JY,JX)                  !surface irrigation  FeHPO4 concentration, [g m-3]
  real(r8) :: CF2PQ(366,JY,JX)                  !surface irrigation  FeH2PO4 concentration, [g m-3]
  real(r8) :: CC0PQ(366,JY,JX)                  !surface irrigation  CaPO4 concentration, [g m-3]
  real(r8) :: CC1PQ(366,JY,JX)                  !surface irrigation  CaHPO4 concentration, [g m-3]
  real(r8) :: CC2PQ(366,JY,JX)                  !surface irrigation  CaH2PO4 concentration, [g m-3]
  real(r8) :: CM1PQ(366,JY,JX)                  !surface irrigation  MgHPO4 concentration, [g m-3]
  real(r8) :: CSTRQ(366,JY,JX)                  !surface irrigation ion strength, [g m-3]
  real(r8) :: CSTRR(JY,JX)                      !surface irrigation ion strength, [g m-3]

  real(r8) :: CCOU                              !subsurface irrigation  CO2 concentration	[g m-3]
  real(r8) :: CCHU                              !subsurface irrigation  CH4 concentration	[g m-3]
  real(r8) :: COXU                              !subsurface irrigation  O2 concentration	[g m-3]
  real(r8) :: CNNU                              !subsurface irrigation  N2 concentration	[g m-3]
  real(r8) :: CN2U                              !subsurface irrigation  N2O concentration	[g m-3]
  real(r8) :: CCOQ(JY,JX)                       !surface irrigation  CO2 concentration, [g m-3]
  real(r8) :: CCHQ(JY,JX)                       !surface irrigation  CH4 concentration, [g m-3]
  real(r8) :: COXQ(JY,JX)                       !surface irrigation  O2 concentration, [g m-3]
  real(r8) :: CNNQ(JY,JX)                       !surface irrigation  N2 concentration, [g m-3]
  real(r8) :: CN2Q(JY,JX)                       !surface irrigation  N2O concentration, [g m-3]
  real(r8) :: COCU(0:4,JZ,JY,JX)                !subsurface irrigation  DOC concentration, [g m-3]
  real(r8) :: CONU(0:4,JZ,JY,JX)                !subsurface irrigation  DON concentration, [g m-3]
  real(r8) :: COAU(0:4,JZ,JY,JX)                !subsurface irrigation  acetate concentration, [g m-3]
  real(r8) :: CN4U(JZ,JY,JX)                    !subsurface irrigation  NH4 concentration, [g m-3]
  real(r8) :: CN3U(JZ,JY,JX)                    !subsurface irrigation  NH3 concentration, [g m-3]
  real(r8) :: CNOU(JZ,JY,JX)                    !subsurface irrigation  NO3 concentration, [g m-3]
  real(r8) :: CH2PU(JZ,JY,JX)                   !subsurface irrigation H2PO4 concentrations, [g m-3]
  real(r8) :: CNZU(JZ,JY,JX)                    !subsurface irrigation  Al concentration, [g m-3]
  real(r8) :: CALU(JZ,JY,JX)                    !subsurface irrigation  Fe concentration, [g m-3]
  real(r8) :: CFEU(JZ,JY,JX)                    !subsurface irrigation  H concentration, [g m-3]
  real(r8) :: CHYU(JZ,JY,JX)                    !subsurface irrigation  Ca concentration, [g m-3]
  real(r8) :: CCAU(JZ,JY,JX)                    !subsurface irrigation  Mg concentration, [g m-3]
  real(r8) :: CMGU(JZ,JY,JX)                    !subsurface irrigation  Na concentration, [g m-3]
  real(r8) :: CNAU(JZ,JY,JX)                    !subsurface irrigation  K concentration, [g m-3]
  real(r8) :: CKAU(JZ,JY,JX)                    !subsurface irrigation  OH concentration, [g m-3]
  real(r8) :: COHU(JZ,JY,JX)                    !subsurface irrigation  SO4 concentration, [g m-3]
  real(r8) :: CSOU(JZ,JY,JX)                    !subsurface irrigation  Cl concentration, [g m-3]
  real(r8) :: CCLU(JZ,JY,JX)                    !subsurface irrigation  CO3 concentration, [g m-3]
  real(r8) :: CC3U(JZ,JY,JX)                    !subsurface irrigation  HCO3 concentration, [g m-3]
  real(r8) :: CHCU(JZ,JY,JX)                    !subsurface irrigation  CH4 concentration, [g m-3]
  real(r8) :: CAL1U(JZ,JY,JX)                   !subsurface irrigation  AlOH concentration, [g m-3]
  real(r8) :: CAL2U(JZ,JY,JX)                   !subsurface irrigation  AlOH2 concentration, [g m-3]
  real(r8) :: CAL3U(JZ,JY,JX)                   !subsurface irrigation  AlOH3 concentration, [g m-3]
  real(r8) :: CAL4U(JZ,JY,JX)                   !subsurface irrigation  AlOH4 concentration, [g m-3]
  real(r8) :: CALSU(JZ,JY,JX)                   !subsurface irrigation  AlSO4 concentration, [g m-3]
  real(r8) :: CFE1U(JZ,JY,JX)                   !subsurface irrigation  FeOH concentration, [g m-3]
  real(r8) :: CFE2U(JZ,JY,JX)                   !subsurface irrigation  FeOH2 concentration, [g m-3]
  real(r8) :: CFE3U(JZ,JY,JX)                   !subsurface irrigation  FeOH3 concentration, [g m-3]
  real(r8) :: CFE4U(JZ,JY,JX)                   !subsurface irrigation  FeOH4 concentration, [g m-3]
  real(r8) :: CFESU(JZ,JY,JX)                   !subsurface irrigation  FeSO4 concentration, [g m-3]
  real(r8) :: CCAOU(JZ,JY,JX)                   !subsurface irrigation  CaOH concentration, [g m-3]
  real(r8) :: CCACU(JZ,JY,JX)                   !subsurface irrigation  CaCO3 concentration, [g m-3]
  real(r8) :: CCAHU(JZ,JY,JX)                   !subsurface irrigation  CaHCO3 concentration, [g m-3]
  real(r8) :: CCASU(JZ,JY,JX)                   !subsurface irrigation  CaSO4 concentration, [g m-3]
  real(r8) :: CMGOU(JZ,JY,JX)                   !subsurface irrigation  MgOH concentration, [g m-3]
  real(r8) :: CMGCU(JZ,JY,JX)                   !subsurface irrigation  MgCO3 concentration, [g m-3]
  real(r8) :: CMGHU(JZ,JY,JX)                   !subsurface irrigation  MgHCO3 concentration, [g m-3]
  real(r8) :: CMGSU(JZ,JY,JX)                   !subsurface irrigation  MgSO4 concentration, [g m-3]
  real(r8) :: CNACU(JZ,JY,JX)                   !subsurface irrigation  NaCO3 concentration, [g m-3]
  real(r8) :: CNASU(JZ,JY,JX)                   !subsurface irrigation  NaSO4 concentration, [g m-3]
  real(r8) :: CKASU(JZ,JY,JX)                   !subsurface irrigation  K concentration, [g m-3]
  real(r8) :: CH0PU(JZ,JY,JX)                   !subsurface irrigation  PO4 concentration, [g m-3]
  real(r8) :: CH1PU(JZ,JY,JX)                   !subsurface irrigation  HPO4 concentration, [g m-3]
  real(r8) :: CH3PU(JZ,JY,JX)                   !subsurface irrigation  H3PO4 concentration, [g m-3]
  real(r8) :: CF1PU(JZ,JY,JX)                   !subsurface irrigation  FeHPO4 concentration, [g m-3]
  real(r8) :: CF2PU(JZ,JY,JX)                   !subsurface irrigation  FeH2PO4 concentration, [g m-3]
  real(r8) :: CC0PU(JZ,JY,JX)                   !subsurface irrigation  CaPO4 concentration, [g m-3]
  real(r8) :: CC1PU(JZ,JY,JX)                   !subsurface irrigation  CaHPO4 concentration, [g m-3]
  real(r8) :: CC2PU(JZ,JY,JX)                   !subsurface irrigation  CaH2PO4 concentration, [g m-3]
  real(r8) :: CM1PU(JZ,JY,JX)                   !subsurface irrigation  MgHPO4 concentration, [g m-3]
  real(r8) :: COPU(0:4,JZ,JY,JX)                !subsurface irrigation  DOP concentration, [g m-3]
  real(r8) :: TDIRI(JY,JX,12)                   !accumulated change  for irrigation, [-]
  real(r8) :: DIRRI(12)                         !change factor for irrigation, [-]


  private :: InitAllocate

  contains

  subroutine InitFertilizerData

  implicit none

  call InitAllocate

  end subroutine InitFertilizerData

  subroutine InitAllocate

  implicit none



  allocate(ZNH4FA(0:JZ,JY,JX))
  allocate(ZNH3FA(0:JZ,JY,JX))
  allocate(ZNHUFA(0:JZ,JY,JX))
  allocate(ZNO3FA(0:JZ,JY,JX))
  allocate(ZNH4FB(0:JZ,JY,JX))
  allocate(ZNH3FB(0:JZ,JY,JX))
  allocate(ZNHUFB(0:JZ,JY,JX))
  allocate(ZNO3FB(0:JZ,JY,JX))

  end subroutine InitAllocate

  subroutine DestructFertilizerData
  implicit none


  deallocate(ZNH4FA)
  deallocate(ZNH3FA)
  deallocate(ZNHUFA)
  deallocate(ZNO3FA)
  deallocate(ZNH4FB)
  deallocate(ZNH3FB)
  deallocate(ZNHUFB)
  deallocate(ZNO3FB)

  end subroutine DestructFertilizerData

end module FertilizerDataType
