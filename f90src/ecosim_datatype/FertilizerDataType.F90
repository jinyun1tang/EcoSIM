module FertilizerDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
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

  real(r8),allocatable :: DCORP(:,:,:)                  !soil mixing fraction with tillage, [-]
  real(r8),allocatable :: FERT(:,:,:,:)                !fertilizer application, [g m-2]
  real(r8),allocatable :: FDPTH(:,:,:)                  !depth of fertilizer application, [m]
  real(r8),allocatable :: ROWI(:,:,:)                   !row spacing of fertilizer band, [m]

  real(r8), allocatable :: ROWN(:,:)                       !row spacing of NH4 fertilizer band, [m]
  real(r8), allocatable :: ROWO(:,:)                       !row spacing of NO3 fertilizer band, [m]
  real(r8), allocatable :: ROWP(:,:)                       !row spacing of PO4 fertilizer band, [m]

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
  use abortutils, only : destroy

  implicit none


  call destroy(ZNH4FA)
  call destroy(ZNH3FA)
  call destroy(ZNHUFA)
  call destroy(ZNO3FA)
  call destroy(ZNH4FB)
  call destroy(ZNH3FB)
  call destroy(ZNHUFB)
  call destroy(ZNO3FB)

  call destroy(ROWN)                       !row spacing of NH4 fertilizer band, [m]
  call destroy(ROWO)                       !row spacing of NO3 fertilizer band, [m]
  call destroy(ROWP)                       !row spacing of PO4 fertilizer band, [m]

  call destroy(DCORP)                  !soil mixing fraction with tillage, [-]
  call destroy(FERT)                !fertilizer application, [g m-2]
  call destroy(FDPTH)                  !depth of fertilizer application, [m]
  call destroy(ROWI)                   !row spacing of fertilizer band, [m]

  end subroutine DestructFertilizerData

end module FertilizerDataType
