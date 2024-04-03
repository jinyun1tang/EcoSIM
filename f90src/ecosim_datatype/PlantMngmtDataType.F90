module PlantMngmtDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  logical, target,allocatable ::  flag_pft_active(:,:,:)
  real(r8),target,allocatable ::  THIN_pft(:,:,:,:)                      !thinning of plant population, [-]
  real(r8),target,allocatable ::  EHVST(:,:,:,:,:,:)                 !harvest efficiency, [-]
  real(r8),target,allocatable ::  HVST(:,:,:,:)                      !harvest cutting height (+ve) or fractional LAI removal (-ve), [m or -]
  integer ,target,allocatable ::  iHarvstType_pft(:,:,:,:)                      !type of harvest, [-]
  integer ,target,allocatable ::  jHarvst_pft(:,:,:,:)                      !flag for stand replacing disturbance, [-]

  integer ,target,allocatable ::  iYearPlanting_pft(:,:,:)                         !year of planting, [-]
  integer ,target,allocatable ::  iYearPlantHarvest_pft(:,:,:)                         !year of harvest, [-]
  integer ,target,allocatable ::  iDayPlanting_pft(:,:,:)                        !day of planting, [-]
  integer ,target,allocatable ::  iDayPlantHarvest_pft(:,:,:)                        !day of harvest, [-]
  integer ,target,allocatable ::  iPlantState_pft(:,:,:)                         !flag for species death, [-]
  integer ,target,allocatable ::  IYRX(:,:,:)                         !alternate year of planting, [-]
  integer ,target,allocatable ::  IDAYX(:,:,:)                        !alternate day of planting, [-]
  integer ,target,allocatable ::  IYRY(:,:,:)                         !alternate year of harvest, [-]
  integer ,target,allocatable ::  IDAYY(:,:,:)                        !alternate day of harvest, [-]
  real(r8),target,allocatable ::  CO2byFire_col(:,:)                         !total CO2 flux from fire, [g d-2]
  real(r8),target,allocatable ::  CH4byFire_col(:,:)                         !total CH4 flux from fire, [g d-2]
  real(r8),target,allocatable ::  O2byFire_col(:,:)                         !total O2 flux from fire, [g d-2]
  real(r8),target,allocatable ::  NH3byFire_col(:,:)                         !total NH3 flux from fire, [g d-2]
  real(r8),target,allocatable ::  N2ObyFire_col(:,:)                         !total N2O flux from fire, [g d-2]
  real(r8),target,allocatable ::  PO4byFire_col(:,:)                         !total PO4 flux from fire, [g d-2]

  private :: InitAllocate
  contains

  subroutine InitPlantMngmtData
  implicit none

  call InitAllocate

  end subroutine InitPlantMngmtData

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(flag_pft_active(JP,JY,JX));  flag_pft_active=.false.
  allocate(THIN_pft(JP,366,JY,JX)); THIN_pft=0._r8
  allocate(EHVST(2,4,JP,366,JY,JX));EHVST=0._r8
  allocate(HVST(JP,366,JY,JX)); HVST=0._r8
  allocate(iHarvstType_pft(JP,366,JY,JX));iHarvstType_pft=0
  allocate(jHarvst_pft(JP,366,JY,JX));jHarvst_pft=0

  allocate(iYearPlanting_pft(JP,JY,JX));     iYearPlanting_pft=0
  allocate(iYearPlantHarvest_pft(JP,JY,JX));     iYearPlantHarvest_pft=0
  allocate(iDayPlanting_pft(JP,JY,JX));    iDayPlanting_pft=0
  allocate(iDayPlantHarvest_pft(JP,JY,JX));    iDayPlantHarvest_pft=0
  allocate(iPlantState_pft(JP,JY,JX));     iPlantState_pft=0
  allocate(IYRX(JP,JY,JX));     IYRX=0
  allocate(IDAYX(JP,JY,JX));    IDAYX=0
  allocate(IYRY(JP,JY,JX));     IYRY=0
  allocate(IDAYY(JP,JY,JX));    IDAYY=0
  allocate(CO2byFire_col(JY,JX));       CO2byFire_col=0._r8
  allocate(CH4byFire_col(JY,JX));       CH4byFire_col=0._r8
  allocate(O2byFire_col(JY,JX));       O2byFire_col=0._r8
  allocate(NH3byFire_col(JY,JX));       NH3byFire_col=0._r8
  allocate(N2ObyFire_col(JY,JX));       N2ObyFire_col=0._r8
  allocate(PO4byFire_col(JY,JX));       PO4byFire_col=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructPlantMngmtData
  use abortutils, only : destroy
  implicit none

  call destroy(flag_pft_active)
  call destroy(THIN_pft)
  call destroy(EHVST)
  call destroy(HVST)
  call destroy(iHarvstType_pft)
  call destroy(jHarvst_pft)

  call destroy(iYearPlanting_pft)
  call destroy(iYearPlantHarvest_pft)
  call destroy(iDayPlanting_pft)
  call destroy(iDayPlantHarvest_pft)
  call destroy(iPlantState_pft)
  call destroy(IYRX)
  call destroy(IDAYX)
  call destroy(IYRY)
  call destroy(IDAYY)
  call destroy(CO2byFire_col)
  call destroy(CH4byFire_col)
  call destroy(O2byFire_col)
  call destroy(NH3byFire_col)
  call destroy(N2ObyFire_col)
  call destroy(PO4byFire_col)
  end subroutine DestructPlantMngmtData

end module PlantMngmtDataType
