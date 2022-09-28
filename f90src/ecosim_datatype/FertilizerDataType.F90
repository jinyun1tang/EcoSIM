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

  real(r8) :: DCORP(366,JY,JX)                  !soil mixing fraction with tillage, [-]
  real(r8) :: FERT(20,366,JY,JX)                !fertilizer application, [g m-2]
  real(r8) :: FDPTH(366,JY,JX)                  !depth of fertilizer application, [m]
  real(r8) :: ROWI(366,JY,JX)                   !row spacing of fertilizer band, [m]

  real(r8) :: ROWN(JY,JX)                       !row spacing of NH4 fertilizer band, [m]
  real(r8) :: ROWO(JY,JX)                       !row spacing of NO3 fertilizer band, [m]
  real(r8) :: ROWP(JY,JX)                       !row spacing of PO4 fertilizer band, [m]

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
