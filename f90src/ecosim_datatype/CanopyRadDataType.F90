module CanopyRadDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GridConsts
  implicit none
  public
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  real(r8),target,allocatable :: ZSIN(:)                        !sine of leaf angle	-
  real(r8),target,allocatable :: ZCOS(:)                        !cosine of leaf angle	-
  real(r8),target,allocatable :: OMEGA(:,:,:)                   !sine of indirect sky radiation on leaf surface
  real(r8),target,allocatable :: OMEGX(:,:,:)                   !sine of indirect sky radiation on leaf surface/sine of indirect sky radiation
  integer,target,allocatable :: IALBY(:,:,:)                    !flag for calculating backscattering of radiation in canopy
  real(r8),target,allocatable :: PARDIF(:,:,:,:,:,:)           !diffuse incoming PAR, [umol m-2 s-1]
  real(r8),target,allocatable :: PAR(:,:,:,:,:,:)              !direct incoming PAR, [umol m-2 s-1]
  real(r8),target,allocatable :: CLASS(:,:,:,:)                !fractionction of leaves in different angle classes, [-]
  real(r8),target,allocatable :: LeafA_lyrnodbrchpft(:,:,:,:,:,:,:)           !leaf surface area, [m2 d-2]
  real(r8),target,allocatable :: LeafAUnshaded_seclyrnodbrpft(:,:,:,:,:,:,:)          !leaf irradiated surface area, [m2 d-2]
  real(r8),target,allocatable :: StemA_lyrnodbrchpft(:,:,:,:,:,:)            !stem surface area, [m2 d-2]

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

  dangle=PICON2h/real(JLI,r8)         !the angle section width

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
  integer :: ncols
  ncols = bounds%ncols

  allocate(ZSIN(JLI))
  allocate(ZCOS(JLI))
  allocate(OMEGA(JSA,JLI,JLA))
  allocate(OMEGX(JSA,JLI,JLA))
  allocate(IALBY(JSA,JLI,JLA))
  allocate(CLASS(JLI,JP,JY,JX))
  allocate(LeafA_lyrnodbrchpft(JLI,JC,JNODS,JBR,JP,JY,JX))
  allocate(LeafAUnshaded_seclyrnodbrpft(JLI,JC,JNODS,JBR,JP,JY,JX))
  allocate(PAR(JLI,JSA,JC,JP,JY,JX))
  allocate(PARDIF(JLI,JSA,JC,JP,JY,JX))
  allocate(StemA_lyrnodbrchpft(JLI,JC,JBR,JP,JY,JX))

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructCanopyRad
  use abortutils, only : destroy
  implicit none
  call destroy(ZSIN)
  call destroy(ZCOS)
  call destroy(OMEGA)
  call destroy(OMEGX)
  call destroy(IALBY)
  call destroy(CLASS)
  call destroy(LeafA_lyrnodbrchpft)
  call destroy(LeafAUnshaded_seclyrnodbrpft)
  call destroy(PAR)
  call destroy(PARDIF)
  call destroy(StemA_lyrnodbrchpft)
  end subroutine DestructCanopyRad
end module CanopyRadDataType
