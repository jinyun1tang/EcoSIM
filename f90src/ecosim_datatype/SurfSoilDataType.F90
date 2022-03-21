module SurfSoilDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  public
  save
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


  real(r8) :: QR(2,2,JV,JH)                     !soil surface runoff water, [m3 d-2 h-1]
  real(r8) :: HQR(2,2,JV,JH)                    !soil surface runoff heat, [MJ d-2 h-1]
  real(r8) :: XCOQRS(2,2,JV,JH)                 !surface runoff CO2 flux, [g d-2 h-1]
  real(r8) :: XCHQRS(2,2,JV,JH)                 !surface brunoff CH4 flux, [g d-2 h-1]
  real(r8) :: XOXQRS(2,2,JV,JH)                 !surface runoff O2 flux, [g d-2 h-1]
  real(r8) :: XNGQRS(2,2,JV,JH)                 !surface runoff N2 flux, [g d-2 h-1]
  real(r8) :: XN2QRS(2,2,JV,JH)                 !surface runoff N2O flux, [g d-2 h-1]
  real(r8) :: XHGQRS(2,2,JV,JH)                 !surface runoff H2 flux, [g d-2 h-1]
  real(r8) :: XN4QRW(2,2,JV,JH)                 !surface runoff NH4 flux non-band, [g d-2 h-1]
  real(r8) :: XN3QRW(2,2,JV,JH)                 !surface runoff NH3 flux non-band, [g d-2 h-1]
  real(r8) :: XNOQRW(2,2,JV,JH)                 !surface runoff NO3 flux non-band, [g d-2 h-1]
  real(r8) :: XNXQRS(2,2,JV,JH)                 !surface runoff NO2 flux, [g d-2 h-1]
  real(r8) :: XP4QRW(2,2,JV,JH)                 !surface runoff PO4 flux, [g d-2 h-1]
  real(r8) :: XP1QRW(2,2,JV,JH)
  real(r8) :: XOCQRS(0:4,2,2,JV,JH)             !surface runoff DOC flux, [g d-2 h-1]
  real(r8) :: XONQRS(0:4,2,2,JV,JH)             !surface runoff DON flux, [g d-2 h-1]
  real(r8) :: XOPQRS(0:4,2,2,JV,JH)             !surface runoff DOP flux, [g d-2 h-1]
  real(r8) :: XOAQRS(0:4,2,2,JV,JH)             !surface runoff acetate flux, [g d-2 h-1]

  real(r8) :: FLQGQ(JY,JX)                      !precipitation flux into soil surface , [m3 d-2 h-1]
  real(r8) :: FLQGI(JY,JX)                      !irrifation flux into soil surface , [m3 d-2 h-1]
  real(r8) :: FLWNU(JY,JX)                      !lake surface water flux
  real(r8) :: FLWXNU(JY,JX)                     !lake surface water flux
  real(r8) :: FLWHNU(JY,JX)                     !lake surface water flux
  real(r8) :: HFLWNU(JY,JX)                     !lake surface heat flux

end module SurfSoilDataType
