module SOMDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
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
  real(r8),allocatable :: ORGC(:,:,:)
  real(r8),allocatable :: ORGN(:,:,:)
  integer,allocatable :: IXTYP(:,:,:)
  real(r8),allocatable :: OMCI(:,:)
  real(r8),allocatable :: OMCF(:)  !hetero microbial biomass composition in SOC
  real(r8),allocatable :: OMCA(:)  !autotrophic microbial biomass composition in SOC
  real(r8),allocatable :: RC0(:,:,:)
  real(r8),allocatable :: ORGCX(:,:,:)
  real(r8),allocatable :: OSA(:,:,:,:,:)
  real(r8),allocatable :: ORGR(:,:,:)
  real(r8),allocatable :: CORGC(:,:,:)
  real(r8),allocatable :: CORGN(:,:,:)
  real(r8),allocatable :: CORGP(:,:,:)
  real(r8),allocatable :: CORGR(:,:,:)
  real(r8),allocatable :: CFOMC(:,:,:,:)

  private :: InitAllocate
  contains

  subroutine InitSOMData

  implicit none

  call InitAllocate

  CNRH=(/3.33E-02_r8,3.33E-02_r8,3.33E-02_r8,5.00E-02_r8,12.50E-02_r8/)
  CPRH=(/3.33E-03_r8,3.33E-03_r8,3.33E-03_r8,5.00E-03_r8,12.50E-03_r8/)
  OMCF=(/0.20,0.20,0.30,0.20,0.050,0.025,0.025/)
  OMCA=(/0.06,0.02,0.01,0.0,0.01,0.0,0.0/)
  end subroutine InitSOMData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none

  include "parameters.h"

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
  allocate(OMCI(3,0:4))
  allocate(OMCF(7))
  allocate(OMCA(7))
  allocate(IXTYP(2,JY,JX))
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

end module SOMDataType
