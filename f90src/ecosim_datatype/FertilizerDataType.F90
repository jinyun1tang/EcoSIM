module FertilizerDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable :: FertN_soil_vr(:,:,:,:)
  real(r8),target,allocatable :: FertN_Band_vr(:,:,:,:)

  real(r8),target,allocatable :: DCORP(:,:,:)                  !soil mixing fraction with tillage, [-]
  real(r8),target,allocatable :: FERT(:,:,:,:)                !fertilizer application, [g m-2]
  real(r8),target,allocatable :: FDPTH(:,:,:)                  !depth of fertilizer application, [m]
  real(r8),target,allocatable :: ROWI(:,:,:)                   !row spacing of fertilizer band, [m]

  real(r8),target,allocatable :: ROWN(:,:)                       !row spacing of NH4 fertilizer band, [m]
  real(r8),target,allocatable :: ROWO(:,:)                       !row spacing of NO3 fertilizer band, [m]
  real(r8),target,allocatable :: ROWP(:,:)                       !row spacing of PO4 fertilizer band, [m]

  private :: InitAllocate

  contains

  subroutine InitFertilizerData

  implicit none

  call InitAllocate

  end subroutine InitFertilizerData

  subroutine InitAllocate

  implicit none

  allocate(ROWN(JY,JX))                       !row spacing of NH4 fertilizer band, [m]
  allocate(ROWO(JY,JX))                       !row spacing of NO3 fertilizer band, [m]
  allocate(ROWP(JY,JX))                       !row spacing of PO4 fertilizer band, [m]

  allocate(DCORP(366,JY,JX))                  !soil mixing fraction with tillage, [-]
  allocate(FERT(20,366,JY,JX))                !fertilizer application, [g m-2]
  allocate(FDPTH(366,JY,JX))                  !depth of fertilizer application, [m]
  allocate(ROWI(366,JY,JX))                   !row spacing of fertilizer band, [m]

  allocate(FertN_soil_vr(ifertn_beg:ifertn_end,0:JZ,JY,JX)); FertN_soil_vr=0._r8
  allocate(FertN_Band_vr(ifertnb_beg:ifertnb_end,1:JZ,JY,JX)); FertN_Band_vr=0._r8

  end subroutine InitAllocate
!-------------------------------------------------------
  subroutine DestructFertilizerData
  use abortutils, only : destroy

  implicit none


  call destroy(FertN_soil_vr)
  call destroy(FertN_Band_vr)

  call destroy(ROWN)                       !row spacing of NH4 fertilizer band, [m]
  call destroy(ROWO)                       !row spacing of NO3 fertilizer band, [m]
  call destroy(ROWP)                       !row spacing of PO4 fertilizer band, [m]

  call destroy(DCORP)                  !soil mixing fraction with tillage, [-]
  call destroy(FERT)                !fertilizer application, [g m-2]
  call destroy(FDPTH)                  !depth of fertilizer application, [m]
  call destroy(ROWI)                   !row spacing of fertilizer band, [m]

  end subroutine DestructFertilizerData

end module FertilizerDataType
