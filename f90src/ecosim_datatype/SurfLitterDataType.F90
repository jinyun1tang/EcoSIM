module SurfLitterDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx=> jcplxc, jcplx1=>jcplx1c
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) ,target,allocatable ::   BKRS(:)                              !surface litter bulk density	[Mg m-3]
  real(r8) ,target,allocatable ::   PARR(:,:)                        !surface litter boundary layer conductance, [m t-1]
  integer  ,target,allocatable ::   IXTYP(:,:,:)                     !surface litter type:1=plant,2=manure
  real(r8) ,target,allocatable ::   XCORP(:,:)                       !factor for surface litter incorporation and soil mixing
  real(r8) ,target,allocatable ::   FLWRM(:,:,:)                     !water transfer between soil surface and surface litter, [g d-2 t-1]
  real(r8) ,target,allocatable ::   FLQRM(:,:,:)                     !meltwater flux into surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   CVRD(:,:)                        !fraction of soil surface covered by surface litter, [-]
  real(r8) ,target,allocatable ::   HFLWR(:,:)                       !net heat transfer to surface litter, [MJ d-2 t-1]
  real(r8) ,target,allocatable ::   VLitR(:,:)                        !surface litter volume, [m3 d-2]
  real(r8) ,target,allocatable ::   VHCPRX(:,:)                      !surface litter heat capacity from previous time step, [MJ d-2 K-1]
  real(r8) ,target,allocatable ::   VWatLitrX(:,:)                      !surface litter water holding capacity, [m3 d-2]
  real(r8) ,target,allocatable ::   FLWR(:,:)                        !net water transfer to surface litter, [MJ d-2 t-1]
  real(r8) ,target,allocatable ::   THAWR(:,:)                       !freeze (-ve) - thaw (+ve) in surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   HTHAWR(:,:)                      !latent heat of freeze (-ve) - thaw (+ve) in surface litter, [MJ d-2 h-1]
  real(r8) ,target,allocatable ::   FLQRQ(:,:)                       !precipitation flux into surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   FLQRI(:,:)                       !irrigation flux into surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   POROS0(:,:)                      !litter porosity
  real(r8) ,target,allocatable ::   RC0(:,:,:)                       !surface litter in each complex	g d-2
  real(r8),target,allocatable ::  RC0ff(:,:)

  private :: InitAllocate
  contains
!----------------------------------------------------------------------

  subroutine InitSurfLitter(n_litrsfk)

  implicit none
  integer, intent(in) :: n_litrsfk

  call InitAllocate(n_litrsfk)
  end subroutine InitSurfLitter
!----------------------------------------------------------------------

  subroutine InitAllocate(n_litrsfk)

  implicit none
  integer, intent(in) :: n_litrsfk

  allocate(BKRS(1:n_litrsfk));    BKRS =0._r8
  allocate(PARR(JY,JX));         PARR=0._r8
  allocate(IXTYP(2,JY,JX));      IXTYP=0
  allocate(XCORP(JY,JX));        XCORP=0._r8
  allocate(FLWRM(60,JY,JX));     FLWRM=0._r8
  allocate(FLQRM(60,JY,JX));     FLQRM=0._r8
  allocate(CVRD(JY,JX));         CVRD=0._r8
  allocate(HFLWR(JY,JX));        HFLWR=0._r8
  allocate(VLitR(JY,JX));         VLitR=0._r8
  allocate(VHCPRX(JY,JX));       VHCPRX=0._r8
  allocate(VWatLitrX(JY,JX));       VWatLitrX=0._r8
  allocate(FLWR(JY,JX));         FLWR=0._r8
  allocate(THAWR(JY,JX));        THAWR=0._r8
  allocate(HTHAWR(JY,JX));       HTHAWR=0._r8
  allocate(FLQRQ(JY,JX));        FLQRQ=0._r8
  allocate(FLQRI(JY,JX));        FLQRI=0._r8
  allocate(POROS0(JY,JX));       POROS0=0._r8
  allocate(RC0(1:jcplx,JY,JX));      Rc0=0._r8
  allocate(RC0ff(JY,JX)); RC0ff=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSurfLitter
  use abortutils, only : destroy
  implicit none

  call destroy(BKRS)
  call destroy(PARR)
  call destroy(IXTYP)
  call destroy(XCORP)
  call destroy(FLWRM)
  call destroy(FLQRM)
  call destroy(CVRD)
  call destroy(HFLWR)
  call destroy(VLitR)
  call destroy(VHCPRX)
  call destroy(VWatLitrX)
  call destroy(FLWR)
  call destroy(THAWR)
  call destroy(HTHAWR)
  call destroy(FLQRQ)
  call destroy(FLQRI)
  call destroy(POROS0)
  call destroy(RC0)
  call destroy(RC0ff)

  end subroutine DestructSurfLitter

end module SurfLitterDataType
