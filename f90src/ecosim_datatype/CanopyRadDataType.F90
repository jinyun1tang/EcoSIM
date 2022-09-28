module CanopyRadDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use GridConsts
  implicit none
  public
  character(len=*),private, parameter :: mod_filename = __FILE__
  real(r8),allocatable :: ZSIN(:)                        !sine of leaf angle	-
  real(r8),allocatable :: ZCOS(:)                        !cosine of leaf angle	-
  real(r8),allocatable :: OMEGA(:,:,:)                   !sine of indirect sky radiation on leaf surface
  real(r8),allocatable :: OMEGX(:,:,:)                   !sine of indirect sky radiation on leaf surface/sine of indirect sky radiation
  integer,allocatable :: IALBY(:,:,:)                    !flag for calculating backscattering of radiation in canopy
  real(r8), allocatable :: PARDIF(:,:,:,:,:,:)           !diffuse incoming PAR, [umol m-2 s-1]
  real(r8), allocatable :: PAR(:,:,:,:,:,:)              !direct incoming PAR, [umol m-2 s-1]
  real(r8), allocatable :: CLASS(:,:,:,:)                !fractionction of leaves in different angle classes, [-]
  real(r8), allocatable :: SURF(:,:,:,:,:,:,:)           !leaf surface area, [m2 d-2]
  real(r8), allocatable :: SURFX(:,:,:,:,:,:,:)          !leaf irradiated surface area, [m2 d-2]
  real(r8), allocatable :: SURFB(:,:,:,:,:,:)            !stem surface area, [m2 d-2]

  real(r8) :: TYSIN
  real(r8) :: dangle
  private :: InitAllocate

  contains

  subroutine  InitCanopyRad

  implicit none

  real(r8) :: da
  real(r8) :: aa
  integer :: N
! JLI: number of leaf inclination groups, 90 deg into JLI groups

  call InitAllocate

  dangle=PICON2/real(JLI,r8)         !the angle section width

  DO N = 1, JLI
    aa=real(N-0.5,r8)*dangle
    ZSIN(N)=sin(aa)
    ZCOS(N)=cos(aa)
  ENDDO
  TYSIN = 0._r8
  end subroutine InitCanopyRad

!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none

  !JLS: number of leaf sectors divide the whole horizontal circle

  allocate(ZSIN(JLI))
  allocate(ZCOS(JLI))
  allocate(OMEGA(JSA,JLI,JLA))
  allocate(OMEGX(JSA,JLI,JLA))
  allocate(IALBY(JSA,JLI,JLA))
  allocate(CLASS(JLI,JP,JY,JX))
  allocate(SURF(JLI,JC,JNODS,JC,JP,JY,JX))
  allocate(SURFX(JLI,JC,JNODS,JC,JP,JY,JX))
  allocate(PAR(JLI,JSA,JC,JP,JY,JX))
  allocate(PARDIF(JLI,JSA,JC,JP,JY,JX))
  allocate(SURFB(JLI,JC,JC,JP,JY,JX))

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructCanopyRad
  implicit none
  if(allocated(ZSIN))  deallocate(ZSIN)
  if(allocated(ZCOS))  deallocate(ZCOS)
  if(allocated(OMEGA)) deallocate(OMEGA)
  if(allocated(OMEGX)) deallocate(OMEGX)
  if(allocated(IALBY)) deallocate(IALBY)
  if(allocated(CLASS)) deallocate(CLASS)
  if(allocated(SURF))  deallocate(SURF)
  if(allocated(SURFX)) deallocate(SURFX)
  if(allocated(PAR))   deallocate(PAR)
  if(allocated(PARDIF))deallocate(PARDIF)
  if(allocated(SURFB)) deallocate(SURFB)
  end subroutine DestructCanopyRad
end module CanopyRadDataType
