module VegDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  implicit none
  public

  real(r8),allocatable :: ZSIN(:)
  real(r8),allocatable :: ZCOS(:)
  real(r8),allocatable :: OMEGA(:,:,:)
  real(r8),allocatable :: OMEGX(:,:,:)
  integer,allocatable :: IALBY(:,:,:)

  real(r8) :: TYSIN

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
  TYSIN = 0._r8
  end subroutine InitVegData

!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none
  include "parameters.h"
  !JLS: number of leaf sectors divide the whole horizontal circle

  allocate(ZSIN(JLI))
  allocate(ZCOS(JLI))
  allocate(OMEGA(JSA,JLI,JLA))
  allocate(OMEGX(JSA,JLI,JLA))
  allocate(IALBY(JSA,JLI,JLA))
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructVegData
  implicit none
  deallocate(ZSIN)
  deallocate(ZCOS)
  deallocate(OMEGA)
  deallocate(OMEGX)

  deallocate(IALBY)
  end subroutine DestructVegData
end module VegDataType
