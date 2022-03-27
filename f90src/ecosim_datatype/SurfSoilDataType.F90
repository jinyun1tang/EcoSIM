module SurfSoilDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: FSNW(JY,JX)                       !fraction of snow cover
  real(r8) :: FSNX(JY,JX)                       !fraction of snow-free cover
  real(r8) :: THRMG(JY,JX)                      !longwave radiation emitted from ground surface, [MJ d-2 t-1]
  real(r8) :: HEATI(JY,JX)                      !total net radiation at ground surface, [MJ d-2 t-1]
  real(r8) :: HEATE(JY,JX)                      !total latent heat flux at ground surface, [MJ d-2 t-1]
  real(r8) :: HEATS(JY,JX)                      !total sensible heat flux at ground surface, [MJ d-2 t-1]
  real(r8) :: HEATV(JY,JX)                      !total convective heat flux at ground surface, [MJ d-2 t-1]
  real(r8) :: HEATH(JY,JX)                      !total ground heat flux at ground surface, [MJ d-2 t-1]
  real(r8) :: TEVAPG(JY,JX)                     !total evaporation at ground surface, [m3 d-2 t-1]
  real(r8) :: BARE(JY,JX)                       !fraction of exposed soil surface, [-]
  real(r8) :: VOLWG(JY,JX)                      !surface water storage capacity, [m3 d-2]
  real(r8) :: VOLWD(JY,JX)                      !soil surface water retention capacity
  real(r8) :: VHCPNX(JY,JX)                     !minimum heat capacities
  real(r8) :: PARG(60,JY,JX)                    !soil surface boundary layer conductance, [m t-1]

  real(r8) :: FLQGQ(JY,JX)                      !precipitation flux into soil surface , [m3 d-2 h-1]
  real(r8) :: FLQGI(JY,JX)                      !irrifation flux into soil surface , [m3 d-2 h-1]
  real(r8) :: FLWNU(JY,JX)                      !lake surface water flux
  real(r8) :: FLWXNU(JY,JX)                     !lake surface water flux
  real(r8) :: FLWHNU(JY,JX)                     !lake surface water flux
  real(r8) :: HFLWNU(JY,JX)                     !lake surface heat flux

end module SurfSoilDataType
