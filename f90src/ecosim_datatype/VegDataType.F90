module VegDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  implicit none
  public

  real(r8),allocatable :: ZSIN(:)
  real(r8),allocatable :: ZCOS(:)

  private :: InitAllocate

  contains

  subroutine  InitVegData

  implicit none
  include "parameters.h"
  real(r8) :: da
  real(r8) :: aa
  integer :: N
! JLI: number of leaf inclination groups, 90 deg into JLI groups

  call InitAllocate

  da=PICON2/real(JLI,r8)         !the angle section width

  DO N = 1, JLI
    aa=real(N-0.5,r8)*da
    ZSIN(N)=sin(aa)
    ZCOS(N)=cos(aa)
  ENDDO
  end subroutine InitVegData

!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none
  include "parameters.h"
  !JLS: number of leaf sectors divide the whole horizontal circle

  allocate(ZSIN(JLI))
  allocate(ZCOS(JLI))

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructVegData
  implicit none
  deallocate(ZSIN)
  deallocate(ZCOS)

  end subroutine DestructVegData
end module VegDataType
