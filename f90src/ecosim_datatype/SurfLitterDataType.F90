module SurfLitterDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx=> jcplxc
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: BKRS(0:2)                         !surface litter bulk density	[Mg m-3]
  real(r8) ,allocatable ::  PARR(:,:)                        !surface litter boundary layer conductance, [m t-1]
  integer  ,allocatable ::  IXTYP(:,:,:)                     !surface litter type:1=plant,2=manure
  real(r8) ,allocatable ::  XCORP(:,:)                       !factor for surface litter incorporation and soil mixing
  real(r8) ,allocatable ::  FLWRM(:,:,:)                     !water transfer between soil surface and surface litter, [g d-2 t-1]
  real(r8) ,allocatable ::  FLQRM(:,:,:)                     !meltwater flux into surface litter, [m3 d-2 h-1]
  real(r8) ,allocatable ::  CVRD(:,:)                        !fraction of soil surface covered by surface litter, [-]
  real(r8) ,allocatable ::  HFLWR(:,:)                       !net heat transfer to surface litter, [MJ d-2 t-1]
  real(r8) ,allocatable ::  VOLR(:,:)                        !surface litter volume, [m3 d-2]
  real(r8) ,allocatable ::  VHCPRX(:,:)                      !surface litter heat capacity from previous time step, [MJ d-2 K-1]
  real(r8) ,allocatable ::  VOLWRX(:,:)                      !surface litter water holding capacity, [m3 d-2]
  real(r8) ,allocatable ::  FLWR(:,:)                        !net water transfer to surface litter, [MJ d-2 t-1]
  real(r8) ,allocatable ::  THAWR(:,:)                       !freeze (-ve) - thaw (+ve) in surface litter, [m3 d-2 h-1]
  real(r8) ,allocatable ::  HTHAWR(:,:)                      !latent heat of freeze (-ve) - thaw (+ve) in surface litter, [MJ d-2 h-1]
  real(r8) ,allocatable ::  FLQRQ(:,:)                       !precipitation flux into surface litter, [m3 d-2 h-1]
  real(r8) ,allocatable ::  FLQRI(:,:)                       !irrigation flux into surface litter, [m3 d-2 h-1]
  real(r8) ,allocatable ::  POROS0(:,:)                      !litter porosity
  real(r8) ,allocatable ::  RC0(:,:,:)                       !surface litter in each complex	g d-2
  real(r8),allocatable :: RC0ff(:,:)

  private :: InitAllocate
  contains
!----------------------------------------------------------------------

  subroutine InitSurfLitter

  implicit none

  call InitAllocate
  end subroutine InitSurfLitter
!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(PARR(JY,JX));         PARR=0._r8
  allocate(IXTYP(2,JY,JX));      IXTYP=0
  allocate(XCORP(JY,JX));        XCORP=0._r8
  allocate(FLWRM(60,JY,JX));     FLWRM=0._r8
  allocate(FLQRM(60,JY,JX));     FLQRM=0._r8
  allocate(CVRD(JY,JX));         CVRD=0._r8
  allocate(HFLWR(JY,JX));        HFLWR=0._r8
  allocate(VOLR(JY,JX));         VOLR=0._r8
  allocate(VHCPRX(JY,JX));       VHCPRX=0._r8
  allocate(VOLWRX(JY,JX));       VOLWRX=0._r8
  allocate(FLWR(JY,JX));         FLWR=0._r8
  allocate(THAWR(JY,JX));        THAWR=0._r8
  allocate(HTHAWR(JY,JX));       HTHAWR=0._r8
  allocate(FLQRQ(JY,JX));        FLQRQ=0._r8
  allocate(FLQRI(JY,JX));        FLQRI=0._r8
  allocate(POROS0(JY,JX));       POROS0=0._r8
  allocate(RC0(0:jcplx,JY,JX));      Rc0=0._r8
  allocate(RC0ff(JY,JX)); RC0ff=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSurfLitter
  if (allocated(PARR))     deallocate(PARR)
  if (allocated(IXTYP))    deallocate(IXTYP)
  if (allocated(XCORP))    deallocate(XCORP)
  if (allocated(FLWRM))    deallocate(FLWRM)
  if (allocated(FLQRM))    deallocate(FLQRM)
  if (allocated(CVRD))     deallocate(CVRD)
  if (allocated(HFLWR))    deallocate(HFLWR)
  if (allocated(VOLR))     deallocate(VOLR)
  if (allocated(VHCPRX))   deallocate(VHCPRX)
  if (allocated(VOLWRX))   deallocate(VOLWRX)
  if (allocated(FLWR))     deallocate(FLWR)
  if (allocated(THAWR))    deallocate(THAWR)
  if (allocated(HTHAWR))   deallocate(HTHAWR)
  if (allocated(FLQRQ))    deallocate(FLQRQ)
  if (allocated(FLQRI))    deallocate(FLQRI)
  if (allocated(POROS0))   deallocate(POROS0)
  if (allocated(RC0))      deallocate(RC0)
  if (allocated(RC0ff))     deallocate(RC0ff)

  end subroutine DestructSurfLitter

end module SurfLitterDataType
