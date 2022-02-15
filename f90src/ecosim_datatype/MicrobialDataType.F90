module MicrobialDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
implicit none
  public
  save


  integer :: JG

  real(r8),allocatable :: OMC(:,:,:,:,:,:,:)
  real(r8),allocatable :: OMN(:,:,:,:,:,:,:)
  real(r8),allocatable :: OMP(:,:,:,:,:,:,:)
  real(r8),allocatable :: ROXYS(:,:,:,:,:,:)
  real(r8),allocatable :: ROQCS(:,:,:,:,:,:)
  real(r8),allocatable :: ROQAS(:,:,:,:,:,:)
  real(r8),allocatable :: CNOMC(:,:,:,:)
  real(r8),allocatable :: CPOMC(:,:,:,:)
  real(r8),allocatable :: RINHO(:,:,:,:,:,:)
  real(r8),allocatable :: RINOO(:,:,:,:,:,:)
  real(r8),allocatable :: RIPOO(:,:,:,:,:,:)
  real(r8),allocatable :: RINHOR(:,:,:,:,:)
  real(r8),allocatable :: RIPOOR(:,:,:,:,:)
  real(r8),allocatable :: RINOOR(:,:,:,:,:)
  real(r8),allocatable :: RVMX4(:,:,:,:,:,:)
  real(r8),allocatable :: RVMX3(:,:,:,:,:,:)
  real(r8),allocatable :: RVMX2(:,:,:,:,:,:)
  real(r8),allocatable :: RVMB4(:,:,:,:,:,:)
  real(r8),allocatable :: RVMB3(:,:,:,:,:,:)
  real(r8),allocatable :: RVMB2(:,:,:,:,:,:)
  real(r8),allocatable :: RVMX1(:,:,:,:,:,:)
  real(r8),allocatable :: RINHB(:,:,:,:,:,:)
  real(r8),allocatable :: RINOB(:,:,:,:,:,:)
  real(r8),allocatable :: RIPBO(:,:,:,:,:,:)
  real(r8),allocatable :: RIPO1(:,:,:,:,:,:)
  real(r8),allocatable :: RIPB1(:,:,:,:,:,:)
  real(r8),allocatable :: RIPO1R(:,:,:,:,:)
  real(r8),allocatable :: OMCER(:,:,:,:,:,:,:)
  real(r8),allocatable :: OMNER(:,:,:,:,:,:,:)
  real(r8),allocatable :: OMPER(:,:,:,:,:,:,:)
  real(r8),allocatable :: OMCI(:,:)
  real(r8) :: FL(2)
  private :: InitAllocate

  contains

  subroutine InitMicrobialData(nguilds)

  implicit none
  integer, intent(in) :: nguilds

  include "parameters.h"

  JG=nguilds
  call InitAllocate()

  end subroutine InitMicrobialData

  subroutine InitAllocate

  implicit none
  include "parameters.h"

  allocate(OMC(3,JG,7,0:5,0:JZ,JY,JX))
  allocate(OMN(3,JG,7,0:5,0:JZ,JY,JX))
  allocate(OMP(3,JG,7,0:5,0:JZ,JY,JX))
  allocate(ROXYS(JG,7,0:5,0:JZ,JY,JX))
  allocate(ROQCS(JG,7,0:4,0:JZ,JY,JX))
  allocate(ROQAS(JG,7,0:4,0:JZ,JY,JX))
  allocate(CNOMC(3,JG,7,0:5))
  allocate(CPOMC(3,JG,7,0:5))
  allocate(RINHO(JG,7,0:5,0:JZ,JY,JX))
  allocate(RINOO(JG,7,0:5,0:JZ,JY,JX))
  allocate(RIPOO(JG,7,0:5,0:JZ,JY,JX))
  allocate(RINHOR(JG,7,0:5,JY,JX))
  allocate(RIPOOR(JG,7,0:5,JY,JX))
  allocate(RINOOR(JG,7,0:5,JY,JX))
  allocate(RVMX4(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMX3(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMX2(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMB4(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMB3(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMB2(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMX1(JG,7,0:5,0:JZ,JY,JX))
  allocate(RINHB(JG,7,0:5,0:JZ,JY,JX))
  allocate(RINOB(JG,7,0:5,0:JZ,JY,JX))
  allocate(RIPBO(JG,7,0:5,0:JZ,JY,JX))
  allocate(RIPO1(JG,7,0:5,0:JZ,JY,JX))
  allocate(RIPB1(JG,7,0:5,0:JZ,JY,JX))
  allocate(RIPO1R(JG,7,0:5,JY,JX))
  allocate(OMCER(3*JG,7,0:5,2,2,JV,JH))
  allocate(OMNER(3*JG,7,0:5,2,2,JV,JH))
  allocate(OMPER(3*JG,7,0:5,2,2,JV,JH))
  allocate(OMCI(3*JG,0:4))

  end subroutine InitAllocate
!----------------------------------------------------------------------------------------------

  subroutine DestructMicrobialData
  implicit none

  deallocate(OMC)
  deallocate(OMN)
  deallocate(OMP)
  deallocate(ROXYS)
  deallocate(ROQCS)
  deallocate(ROQAS)
  deallocate(CNOMC)
  deallocate(CPOMC)
  deallocate(RINHO)
  deallocate(RINOO)
  deallocate(RIPOO)
  deallocate(RINHOR)
  deallocate(RIPOOR)
  deallocate(RINOOR)
  deallocate(RVMX4)
  deallocate(RVMX3)
  deallocate(RVMX2)
  deallocate(RVMB4)
  deallocate(RVMB3)
  deallocate(RVMB2)
  deallocate(RVMX1)
  deallocate(RINHB)
  deallocate(RINOB)
  deallocate(RIPBO)
  deallocate(RIPO1)
  deallocate(RIPB1)
  deallocate(RIPO1R)
  deallocate(OMCER)
  deallocate(OMNER)
  deallocate(OMPER)
  deallocate(OMCI)

  end subroutine DestructMicrobialData

end module MicrobialDataType
