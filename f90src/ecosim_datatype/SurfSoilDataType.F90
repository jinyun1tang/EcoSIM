module SurfSoilDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__
  real(r8),target,allocatable ::  FSNW(:,:)                          !fraction of snow cover
  real(r8),target,allocatable ::  FSNX(:,:)                          !fraction of snow-free cover
  real(r8),target,allocatable ::  THRMG(:,:)                         !longwave radiation emitted from ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HEATI(:,:)                         !total net radiation at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HEATE(:,:)                         !total latent heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HEATS(:,:)                         !total sensible heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HEATV(:,:)                         !total convective heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  HEATH(:,:)                         !total ground heat flux at ground surface, [MJ d-2 t-1]
  real(r8),target,allocatable ::  TEVAPG(:,:)                        !total evaporation at ground surface, [m3 d-2 t-1]
  real(r8),target,allocatable ::  BARE(:,:)                          !fraction of exposed soil surface, [-]
  real(r8),target,allocatable ::  VOLWG(:,:)                         !surface water storage capacity, [m3 d-2]
  real(r8),target,allocatable ::  VOLWD(:,:)                         !soil surface water retention capacity
  real(r8),target,allocatable ::  VHCPNX(:,:)                        !minimum heat capacities
  real(r8),target,allocatable ::  PARG(:,:,:)                        !soil surface boundary layer conductance, [m t-1]
  real(r8),target,allocatable ::  FLQGQ(:,:)                         !precipitation flux into soil surface , [m3 d-2 h-1]
  real(r8),target,allocatable ::  FLQGI(:,:)                         !irrifation flux into soil surface , [m3 d-2 h-1]
  real(r8),target,allocatable ::  LakeSurfFlow(:,:)                         !lake surface water flux
  real(r8),target,allocatable ::  FLWXNU(:,:)                        !lake surface water flux
  real(r8),target,allocatable ::  FLWHNU(:,:)                        !lake surface water flux
  real(r8),target,allocatable ::  LakeSurfHeatFlux(:,:)              !lake surface heat flux, outgoing positive
!----------------------------------------------------------------------

contains
  subroutine InitSurfSoilData

  implicit none
  allocate(FSNW(JY,JX));        FSNW=0._r8
  allocate(FSNX(JY,JX));        FSNX=0._r8
  allocate(THRMG(JY,JX));       THRMG=0._r8
  allocate(HEATI(JY,JX));       HEATI=0._r8
  allocate(HEATE(JY,JX));       HEATE=0._r8
  allocate(HEATS(JY,JX));       HEATS=0._r8
  allocate(HEATV(JY,JX));       HEATV=0._r8
  allocate(HEATH(JY,JX));       HEATH=0._r8
  allocate(TEVAPG(JY,JX));      TEVAPG=0._r8
  allocate(BARE(JY,JX));        BARE=0._r8
  allocate(VOLWG(JY,JX));       VOLWG=0._r8
  allocate(VOLWD(JY,JX));       VOLWD=0._r8
  allocate(VHCPNX(JY,JX));      VHCPNX=0._r8
  allocate(PARG(60,JY,JX));     PARG=0._r8
  allocate(FLQGQ(JY,JX));       FLQGQ=0._r8
  allocate(FLQGI(JY,JX));       FLQGI=0._r8
  allocate(LakeSurfFlow(JY,JX));       LakeSurfFlow=0._r8
  allocate(FLWXNU(JY,JX));      FLWXNU=0._r8
  allocate(FLWHNU(JY,JX));      FLWHNU=0._r8
  allocate(LakeSurfHeatFlux(JY,JX));      LakeSurfHeatFlux=0._r8
  end subroutine InitSurfSoilData

!----------------------------------------------------------------------
  subroutine DestructSurfSoilData
  use abortutils, only : destroy
  implicit none
  call destroy(FSNW)
  call destroy(FSNX)
  call destroy(THRMG)
  call destroy(HEATI)
  call destroy(HEATE)
  call destroy(HEATS)
  call destroy(HEATV)
  call destroy(HEATH)
  call destroy(TEVAPG)
  call destroy(BARE)
  call destroy(VOLWG)
  call destroy(VOLWD)
  call destroy(VHCPNX)
  call destroy(PARG)
  call destroy(FLQGQ)
  call destroy(FLQGI)
  call destroy(LakeSurfFlow)
  call destroy(FLWXNU)
  call destroy(FLWHNU)
  call destroy(LakeSurfHeatFlux)
  end subroutine DestructSurfSoilData

end module SurfSoilDataType
