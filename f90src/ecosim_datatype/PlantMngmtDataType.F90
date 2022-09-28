module PlantMngmtDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable ::  THIN(:,:,:,:)                      !thinning of plant population, [-]
  real(r8),allocatable ::  EHVST(:,:,:,:,:,:)                 !harvest efficiency, [-]
  real(r8),allocatable ::  HVST(:,:,:,:)                      !harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
  integer ,allocatable ::  IHVST(:,:,:,:)                      !type of harvest, [-]
  integer ,allocatable ::  JHVST(:,:,:,:)                      !flag for stand replacing disturbance, [-]

  integer ,allocatable ::  IYR0(:,:,:)                         !year of planting, [-]
  integer ,allocatable ::  IYRH(:,:,:)                         !year of harvest, [-]
  integer ,allocatable ::  IDAY0(:,:,:)                        !day of planting, [-]
  integer ,allocatable ::  IDAYH(:,:,:)                        !day of harvest, [-]
  integer ,allocatable ::  IDTH(:,:,:)                         !flag for species death, [-]
  integer ,allocatable ::  IYRX(:,:,:)                         !alternate year of planting, [-]
  integer ,allocatable ::  IDAYX(:,:,:)                        !alternate day of planting, [-]
  integer ,allocatable ::  IYRY(:,:,:)                         !alternate year of harvest, [-]
  integer ,allocatable ::  IDAYY(:,:,:)                        !alternate day of harvest, [-]
  real(r8),allocatable ::  UCO2F(:,:)                         !total CO2 flux from fire, [g d-2]
  real(r8),allocatable ::  UCH4F(:,:)                         !total CH4 flux from fire, [g d-2]
  real(r8),allocatable ::  UOXYF(:,:)                         !total O2 flux from fire, [g d-2]
  real(r8),allocatable ::  UNH3F(:,:)                         !total NH3 flux from fire, [g d-2]
  real(r8),allocatable ::  UN2OF(:,:)                         !total N2O flux from fire, [g d-2]
  real(r8),allocatable ::  UPO4F(:,:)                         !total PO4 flux from fire, [g d-2]

  private :: InitAllocate
  contains

  subroutine InitPlantMngmtData
  implicit none

  call InitAllocate

  end subroutine InitPlantMngmtData

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(THIN(05,366,JY,JX)); THIN=0._r8
  allocate(EHVST(2,4,05,366,JY,JX));EHVST=0._r8
  allocate(HVST(05,366,JY,JX)); HVST=0._r8
  allocate(IHVST(05,366,JY,JX));IHVST=0
  allocate(JHVST(05,366,JY,JX));JHVST=0

  allocate(IYR0(JP,JY,JX));     IYR0=0
  allocate(IYRH(JP,JY,JX));     IYRH=0
  allocate(IDAY0(JP,JY,JX));    IDAY0=0
  allocate(IDAYH(JP,JY,JX));    IDAYH=0
  allocate(IDTH(JP,JY,JX));     IDTH=0
  allocate(IYRX(JP,JY,JX));     IYRX=0
  allocate(IDAYX(JP,JY,JX));    IDAYX=0
  allocate(IYRY(JP,JY,JX));     IYRY=0
  allocate(IDAYY(JP,JY,JX));    IDAYY=0
  allocate(UCO2F(JY,JX));       UCO2F=0._r8
  allocate(UCH4F(JY,JX));       UCH4F=0._r8
  allocate(UOXYF(JY,JX));       UOXYF=0._r8
  allocate(UNH3F(JY,JX));       UNH3F=0._r8
  allocate(UN2OF(JY,JX));       UN2OF=0._r8
  allocate(UPO4F(JY,JX));       UPO4F=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructPlantMngmtData
  if (allocated(THIN))     deallocate(THIN)
  if (allocated(EHVST))    deallocate(EHVST)
  if (allocated(HVST))     deallocate(HVST)
  if (allocated(IHVST))    deallocate(IHVST)
  if (allocated(JHVST))    deallocate(JHVST)

  if (allocated(IYR0))     deallocate(IYR0)
  if (allocated(IYRH))     deallocate(IYRH)
  if (allocated(IDAY0))    deallocate(IDAY0)
  if (allocated(IDAYH))    deallocate(IDAYH)
  if (allocated(IDTH))     deallocate(IDTH)
  if (allocated(IYRX))     deallocate(IYRX)
  if (allocated(IDAYX))    deallocate(IDAYX)
  if (allocated(IYRY))     deallocate(IYRY)
  if (allocated(IDAYY))    deallocate(IDAYY)
  if (allocated(UCO2F))    deallocate(UCO2F)
  if (allocated(UCH4F))    deallocate(UCH4F)
  if (allocated(UOXYF))    deallocate(UOXYF)
  if (allocated(UNH3F))    deallocate(UNH3F)
  if (allocated(UN2OF))    deallocate(UN2OF)
  if (allocated(UPO4F))    deallocate(UPO4F)
  end subroutine DestructPlantMngmtData

end module PlantMngmtDataType
