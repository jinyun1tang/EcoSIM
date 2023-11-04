module PlantMngmtDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  THIN_pft(:,:,:,:)                      !thinning of plant population, [-]
  real(r8),target,allocatable ::  EHVST(:,:,:,:,:,:)                 !harvest efficiency, [-]
  real(r8),target,allocatable ::  HVST(:,:,:,:)                      !harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
  integer ,target,allocatable ::  IHVST(:,:,:,:)                      !type of harvest, [-]
  integer ,target,allocatable ::  JHVST(:,:,:,:)                      !flag for stand replacing disturbance, [-]

  integer ,target,allocatable ::  iYearPlanting(:,:,:)                         !year of planting, [-]
  integer ,target,allocatable ::  iYearPlantHarvest(:,:,:)                         !year of harvest, [-]
  integer ,target,allocatable ::  iDayPlanting(:,:,:)                        !day of planting, [-]
  integer ,target,allocatable ::  iDayPlantHarvest(:,:,:)                        !day of harvest, [-]
  integer ,target,allocatable ::  iPlantState(:,:,:)                         !flag for species death, [-]
  integer ,target,allocatable ::  IYRX(:,:,:)                         !alternate year of planting, [-]
  integer ,target,allocatable ::  IDAYX(:,:,:)                        !alternate day of planting, [-]
  integer ,target,allocatable ::  IYRY(:,:,:)                         !alternate year of harvest, [-]
  integer ,target,allocatable ::  IDAYY(:,:,:)                        !alternate day of harvest, [-]
  real(r8),target,allocatable ::  UCO2F(:,:)                         !total CO2 flux from fire, [g d-2]
  real(r8),target,allocatable ::  UCH4F(:,:)                         !total CH4 flux from fire, [g d-2]
  real(r8),target,allocatable ::  UOXYF(:,:)                         !total O2 flux from fire, [g d-2]
  real(r8),target,allocatable ::  UNH3F(:,:)                         !total NH3 flux from fire, [g d-2]
  real(r8),target,allocatable ::  UN2OF(:,:)                         !total N2O flux from fire, [g d-2]
  real(r8),target,allocatable ::  UPO4F(:,:)                         !total PO4 flux from fire, [g d-2]

  private :: InitAllocate
  contains

  subroutine InitPlantMngmtData
  implicit none

  call InitAllocate

  end subroutine InitPlantMngmtData

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(THIN_pft(05,366,JY,JX)); THIN_pft=0._r8
  allocate(EHVST(2,4,05,366,JY,JX));EHVST=0._r8
  allocate(HVST(05,366,JY,JX)); HVST=0._r8
  allocate(IHVST(05,366,JY,JX));IHVST=0
  allocate(JHVST(05,366,JY,JX));JHVST=0

  allocate(iYearPlanting(JP,JY,JX));     iYearPlanting=0
  allocate(iYearPlantHarvest(JP,JY,JX));     iYearPlantHarvest=0
  allocate(iDayPlanting(JP,JY,JX));    iDayPlanting=0
  allocate(iDayPlantHarvest(JP,JY,JX));    iDayPlantHarvest=0
  allocate(iPlantState(JP,JY,JX));     iPlantState=0
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
  use abortutils, only : destroy
  implicit none

  call destroy(THIN_pft)
  call destroy(EHVST)
  call destroy(HVST)
  call destroy(IHVST)
  call destroy(JHVST)

  call destroy(iYearPlanting)
  call destroy(iYearPlantHarvest)
  call destroy(iDayPlanting)
  call destroy(iDayPlantHarvest)
  call destroy(iPlantState)
  call destroy(IYRX)
  call destroy(IDAYX)
  call destroy(IYRY)
  call destroy(IDAYY)
  call destroy(UCO2F)
  call destroy(UCH4F)
  call destroy(UOXYF)
  call destroy(UNH3F)
  call destroy(UN2OF)
  call destroy(UPO4F)
  end subroutine DestructPlantMngmtData

end module PlantMngmtDataType
