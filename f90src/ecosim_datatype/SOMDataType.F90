module SOMDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  public
  save

  real(r8),allocatable :: RSC(:,:,:,:)
  real(r8),allocatable :: RSN(:,:,:,:)
  real(r8),allocatable :: RSP(:,:,:,:)
  real(r8),allocatable :: CFOSC(:,:,:,:,:)  !fraction of SOC in kinetic components
  real(r8),allocatable :: CNOSC(:,:,:,:,:)  !N:C,ratios of SOC kinetic components
  real(r8),allocatable :: CPOSC(:,:,:,:,:)  !P:C ratios of SOC kinetic components
  real(r8),allocatable :: CNOFC(:,:)        !fractions to allocate N to kinetic components
  real(r8),allocatable :: CPOFC(:,:)        !fractions to allocate P to kinetic components
  real(r8),allocatable :: CNRH(:)           !default N:C ratios in SOC complexes
  real(r8),allocatable :: CPRH(:)           !default P:C ratios in SOC complexes
  real(r8),allocatable :: OSC(:,:,:,:,:)
  real(r8),allocatable :: OSN(:,:,:,:,:)
  real(r8),allocatable :: OSP(:,:,:,:,:)
  real(r8),allocatable :: OHC(:,:,:,:)
  real(r8),allocatable :: OHN(:,:,:,:)
  real(r8),allocatable :: OHP(:,:,:,:)
  real(r8),allocatable :: OHA(:,:,:,:)
  real(r8),allocatable :: ORC(:,:,:,:,:)
  real(r8),allocatable :: ORN(:,:,:,:,:)
  real(r8),allocatable :: ORP(:,:,:,:,:)
  real(r8),allocatable :: OQC(:,:,:,:)
  real(r8),allocatable :: OQN(:,:,:,:)
  real(r8),allocatable :: OQP(:,:,:,:)
  real(r8),allocatable :: OQA(:,:,:,:)
  real(r8),allocatable :: OQCH(:,:,:,:)
  real(r8),allocatable :: OQNH(:,:,:,:)
  real(r8),allocatable :: OQPH(:,:,:,:)
  real(r8),allocatable :: OQAH(:,:,:,:)
  real(r8),allocatable :: ORGC(:,:,:) !total soil organic C [g d-2]
  real(r8),allocatable :: ORGN(:,:,:) !total soil organic N [g d-2]
  real(r8),allocatable :: OMCF(:)  !hetero microbial biomass composition in SOC
  real(r8),allocatable :: OMCA(:)  !autotrophic microbial biomass composition in SOC
  real(r8),allocatable :: RC0(:,:,:)
  real(r8),allocatable :: ORGCX(:,:,:)
  real(r8),allocatable :: OSA(:,:,:,:,:) !colonized humus C in each complex [g d-2]
  real(r8),allocatable :: ORGR(:,:,:)  !total particulate organic C [g d-2]
  real(r8),allocatable :: CORGC(:,:,:) !soil organic C content [g kg-1]
  real(r8),allocatable :: CORGN(:,:,:) !soil organic N content [mg kg-1]
  real(r8),allocatable :: CORGP(:,:,:) !soil organic P content  [mg kg-1]
  real(r8),allocatable :: CORGR(:,:,:)  !soil particulate C content [g kg-1]
  real(r8),allocatable :: CFOMC(:,:,:,:)
  real(r8) :: EPOC(0:JZ,JY,JX)                  !partitioning coefficient between POC and litter, []
  real(r8) :: EHUM(0:JZ,JY,JX)                  !partitioning coefficient between humus and microbial residue, []
  real(r8) :: COQC(0:4,0:JZ,JY,JX)              !DOC concentration, [g m-3]
  real(r8) :: COQA(0:4,0:JZ,JY,JX)              !acetate concentration, [g m-3]
  real(r8) :: FOSRH(0:4,0:JZ,JY,JX)             !fraction of total organic C in complex, [-]
  real(r8) :: HCO2G(JY,JX)                      !soil CO2 flux, [g d-2 h-1]
  real(r8) :: HCH4G(JY,JX)                      !soil CH4 flux, [g d-2 h-1]
  real(r8) :: HOXYG(JY,JX)                      !soil O2 flux, [g d-2 h-1]
  real(r8) :: HN2OG(JY,JX)                      !soil N2O flux, [g d-2 h-1]
  real(r8) :: HNH3G(JY,JX)                      !soil NH3 flux, [g d-2 h-1]
!  real(r8) :: URAIQ(JY,JX)
  real(r8) :: TOMT(JY,JX)                       !total micriobial C, [g d-2]
  real(r8) :: TONT(JY,JX)                       !total micriobial N, [g d-2]
  real(r8) :: TOPT(JY,JX)                       !total micriobial P, [g d-2]
  real(r8) :: URSDC(JY,JX)                      !total litter C, [g d-2]
  real(r8) :: UORGC(JY,JX)                      !total humus C, [g d-2]
  real(r8) :: UORGF(JY,JX)                      !total C amendment, [g d-2]
  real(r8) :: UXCSN(JY,JX)                      !total litterfall C, [g d-2]
  real(r8) :: UCO2S(JY,JX)                      !total soil DIC, [g d-2]
  real(r8) :: UDOCQ(JY,JX)                      !total surface DOC flux, [g d-2]
  real(r8) :: UDOCD(JY,JX)                      !total subsurface DOC flux, [g d-2]
  real(r8) :: TNBP(JY,JX)                       !total NBP, [g d-2]
  real(r8) :: URSDN(JY,JX)                      !total litter N, [g d-2]
  real(r8) :: UORGN(JY,JX)                      !total humus N, [g d-2]
  real(r8) :: UFERTN(JY,JX)                     !total fertilizer N amendment, [g d-2]
  real(r8) :: UXZSN(JY,JX)                      !total litterfall N, [g d-2]
  real(r8) :: UNH4(JY,JX)                       !total soil NH4 + NH3 content, [g d-2]
  real(r8) :: UNO3(JY,JX)                       !total soil NO3 + NO2 content, [g d-2]
  real(r8) :: UPO4(JY,JX)                       !total soil PO4 content, [g d-2]
  real(r8) :: UDONQ(JY,JX)                      !total surface DON flux, [g d-2]
  real(r8) :: UDOND(JY,JX)                      !total subsurface DON flux, [g d-2]
  real(r8) :: UDOPQ(JY,JX)                      !total surface DOP flux, [g d-2]
  real(r8) :: UDOPD(JY,JX)                      !total subsurface DOP flux, [g d-2]
  real(r8) :: UPP4(JY,JX)                       !total soil precipited P, [g d-2]
  real(r8) :: UN2GS(JY,JX)                      !total N2 fixation, [g d-2]
  real(r8) :: URSDP(JY,JX)                      !total litter P, [g d-2]
  real(r8) :: UORGP(JY,JX)                      !total humus P, [g d-2]
  real(r8) :: UFERTP(JY,JX)                     !total fertilizer P amendment, [g d-2]
  real(r8) :: UH2GG(JY,JX)                      !total H2 flux, []
  real(r8) :: UXPSN(JY,JX)                      !total litterfall P, [g d-2]
  real(r8) :: HN2GG(JY,JX)                      !soil N2 flux, [g d-2 h-1]
  real(r8) :: UN2GG(JY,JX)                      !total soil N2 flux, [g d-2]
  real(r8) :: UION(JY,JX)                       !total soil ion content, [mol d-2]
  real(r8) :: UIONOU(JY,JX)                     !total subsurface ion flux, [mol d-2]
  real(r8) :: UCO2G(JY,JX)                      !total soil CO2 flux, [g d-2]
  real(r8) :: UCH4G(JY,JX)                      !total soil CH4 flux, [g d-2]
  real(r8) :: UOXYG(JY,JX)                      !total soil O2 flux, [g d-2]
  real(r8) :: UNH3G(JY,JX)                      !total soil NH3 flux, [g d-2]
  real(r8) :: UN2OG(JY,JX)                      !total soil N2O flux, [g d-2]
  real(r8) :: ZDRAIN(JY,JX)                     !total N drainage below root zone, [g d-2]
  real(r8) :: PDRAIN(JY,JX)                     !total P drainage below root zone, [g d-2]
  real(r8) :: UCOP(JY,JX)                       !total soil autotrophic respiration, [g d-2]
  real(r8) :: USEDOU(JY,JX)                     !total sediment subsurface flux, [Mg d-2]
  real(r8) :: UDICQ(JY,JX)                      !total surface DIC flux, [g d-2]
  real(r8) :: UDICD(JY,JX)                      !total subsurface DIC flux, [g d-2]
  real(r8) :: UDINQ(JY,JX)                      !total surface DIN flux, [g d-2]
  real(r8) :: UDIND(JY,JX)                      !total subsurface DIN flux, [g d-2]
  real(r8) :: UDIPQ(JY,JX)                      !total surface DIP flux, [g d-2]
  real(r8) :: UDIPD(JY,JX)                      !total subsurface DIP flux, [g d-2]
  real(r8) :: WTSTGT(JY,JX)                     !total standing dead C, [g d-2]
  private :: InitAllocate
  contains

  subroutine InitSOMData

  implicit none

  call InitAllocate

  CNRH=(/3.33E-02_r8,3.33E-02_r8,3.33E-02_r8,5.00E-02_r8,12.50E-02_r8/)
  CPRH=(/3.33E-03_r8,3.33E-03_r8,3.33E-03_r8,5.00E-03_r8,12.50E-03_r8/)
  OMCF=(/0.20_r8,0.20_r8,0.30_r8,0.20_r8,0.050_r8,0.025_r8,0.025_r8/)
  OMCA=(/0.06_r8,0.02_r8,0.01_r8,0.0_r8,0.01_r8,0.0_r8,0.0_r8/)
  end subroutine InitSOMData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none



  allocate(RSC(0:2,0:JZ,JY,JX))
  allocate(RSN(0:2,0:JZ,JY,JX))
  allocate(RSP(0:2,0:JZ,JY,JX))
  allocate(CFOSC(4,0:4,0:JZ,JY,JX))
  allocate(CNOSC(4,0:4,0:JZ,JY,JX))
  allocate(CPOSC(4,0:4,0:JZ,JY,JX))
  allocate(CNOFC(4,0:2))
  allocate(CPOFC(4,0:2))
  allocate(CNRH(0:4))
  allocate(CPRH(0:4))
  allocate(OSC(4,0:4,0:JZ,JY,JX))
  allocate(OSN(4,0:4,0:JZ,JY,JX))
  allocate(OSP(4,0:4,0:JZ,JY,JX))
  allocate(OHC(0:4,0:JZ,JY,JX))
  allocate(OHN(0:4,0:JZ,JY,JX))
  allocate(OHP(0:4,0:JZ,JY,JX))
  allocate(OHA(0:4,0:JZ,JY,JX))
  allocate(ORC(2,0:4,0:JZ,JY,JX))
  allocate(ORN(2,0:4,0:JZ,JY,JX))
  allocate(ORP(2,0:4,0:JZ,JY,JX))
  allocate(OQC(0:4,0:JZ,JY,JX))
  allocate(OQN(0:4,0:JZ,JY,JX))
  allocate(OQP(0:4,0:JZ,JY,JX))
  allocate(OQA(0:4,0:JZ,JY,JX))
  allocate(OQCH(0:4,0:JZ,JY,JX))
  allocate(OQNH(0:4,0:JZ,JY,JX))
  allocate(OQPH(0:4,0:JZ,JY,JX))
  allocate(OQAH(0:4,0:JZ,JY,JX))
  allocate(ORGC(0:JZ,JY,JX))
  allocate(ORGN(0:JZ,JY,JX))
  allocate(OMCF(7))
  allocate(OMCA(7))
!  allocate(IXTYP(2,JY,JX))
  allocate(RC0(0:5,JY,JX))
  allocate(ORGCX(0:JZ,JY,JX))
  allocate(OSA(4,0:4,0:JZ,JY,JX))
  allocate(ORGR(0:JZ,JY,JX))
  allocate(CORGC(0:JZ,JY,JX))
  allocate(CORGN(JZ,JY,JX))
  allocate(CORGP(JZ,JY,JX))
  allocate(CORGR(JZ,JY,JX))
  allocate(CFOMC(2,JZ,JY,JX))
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------

  subroutine DestructSOMData

  implicit none

  deallocate(RSC)
  deallocate(RSN)
  deallocate(RSP)
  deallocate(CFOSC)
  deallocate(CNOSC)
  deallocate(CPOSC)
  deallocate(CNOFC)
  deallocate(CPOFC)
  deallocate(CNRH)
  deallocate(CPRH)
  deallocate(OSC)
  deallocate(OSN)
  deallocate(OSP)
  deallocate(OHC)
  deallocate(OHN)
  deallocate(OHP)
  deallocate(OHA)
  deallocate(ORC)
  deallocate(ORN)
  deallocate(ORP)
  deallocate(OQC)
  deallocate(OQN)
  deallocate(OQP)
  deallocate(OQA)
  deallocate(OQCH)
  deallocate(OQNH)
  deallocate(OQPH)
  deallocate(OQAH)
  deallocate(ORGC)
  deallocate(ORGN)
  deallocate(OMCF)
  deallocate(OMCA)
!  deallocate(IXTYP)
  deallocate(RC0)
  deallocate(ORGCX)
  deallocate(OSA)
  deallocate(ORGR)
  deallocate(CORGC)
  deallocate(CORGN)
  deallocate(CORGP)
  deallocate(CORGR)
  deallocate(CFOMC)

  end subroutine DestructSOMData
end module SOMDataType
