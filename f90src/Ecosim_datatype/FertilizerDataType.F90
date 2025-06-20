module FertilizerDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable :: FertN_mole_soil_vr(:,:,:,:)      !fertilizer in soil from broadcast, [mol d-2]
  real(r8),target,allocatable :: FertN_mole_Band_vr(:,:,:,:)      !fertilizer in band  from side-dressing, [mol d-2]

  real(r8),target,allocatable :: DepzCorp_col(:,:,:)         !soil mixing fraction with tillage, [-]
  real(r8),target,allocatable :: FERT(:,:,:,:)               !fertilizer application, [g m-2]
  real(r8),target,allocatable :: FDPTH(:,:,:)                !depth of fertilizer application, [m]
  real(r8),target,allocatable :: ROWI(:,:,:)                 !row spacing of fertilizer band from side-dressing, [m]

  real(r8),target,allocatable :: ROWSpaceNH4_col(:,:)        !row spacing of NH4 fertilizer band from side-dressing, [m]
  real(r8),target,allocatable :: ROWSpaceNO3_col(:,:)        !row spacing of NO3 fertilizer band from side-dressing, [m]
  real(r8),target,allocatable :: ROWSpacePO4_col(:,:)        !row spacing of PO4 fertilizer band from side-dressing, [m]

  private :: InitAllocate

  contains

  subroutine InitFertilizerData

  implicit none

  call InitAllocate

  end subroutine InitFertilizerData

  subroutine InitAllocate

  implicit none

  allocate(ROWSpaceNH4_col(JY,JX)) ;ROWSpaceNH4_col=0._r8                      !row spacing of NH4 fertilizer band, [m]
  allocate(ROWSpaceNO3_col(JY,JX)) ;ROWSpaceNO3_col=0._r8                      !row spacing of NO3 fertilizer band, [m]
  allocate(ROWSpacePO4_col(JY,JX)) ;ROWSpacePO4_col=0._r8                      !row spacing of PO4 fertilizer band, [m]

  allocate(DepzCorp_col(366,JY,JX)) ;DepzCorp_col=0._r8                 !soil mixing fraction with tillage, [-]
  allocate(FERT(20,366,JY,JX)); FERT=0._r8                !fertilizer application, [g m-2]
  allocate(FDPTH(366,JY,JX)); FDPTH=0._r8                  !depth of fertilizer application, [m]
  allocate(ROWI(366,JY,JX)) ;ROWI=0._r8                  !row spacing of fertilizer band, [m]

  allocate(FertN_mole_soil_vr(ifertn_beg:ifertn_end,0:JZ,JY,JX)); FertN_mole_soil_vr=0._r8
  allocate(FertN_mole_Band_vr(ifertnb_beg:ifertnb_end,1:JZ,JY,JX)); FertN_mole_Band_vr=0._r8

  end subroutine InitAllocate
!-------------------------------------------------------
  subroutine DestructFertilizerData
  use abortutils, only : destroy

  implicit none

  call destroy(FertN_mole_soil_vr)
  call destroy(FertN_mole_Band_vr)

  call destroy(ROWSpaceNH4_col)                       !row spacing of NH4 fertilizer band, [m]
  call destroy(ROWSpaceNO3_col)                       !row spacing of NO3 fertilizer band, [m]
  call destroy(ROWSpacePO4_col)                       !row spacing of PO4 fertilizer band, [m]

  call destroy(DepzCorp_col)                  !soil mixing fraction with tillage, [-]
  call destroy(FERT)                !fertilizer application, [g m-2]
  call destroy(FDPTH)                  !depth of fertilizer application, [m]
  call destroy(ROWI)                   !row spacing of fertilizer band, [m]

  end subroutine DestructFertilizerData

end module FertilizerDataType
