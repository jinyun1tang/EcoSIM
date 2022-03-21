module SurfLitterDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: PARR(JY,JX)                       !surface litter boundary layer conductance, [m t-1]
  real(r8) :: PARG(60,JY,JX)                    !soil surface boundary layer conductance, [m t-1]
  real(r8) :: BKRS(0:2)                         !surface litter bulk density	[Mg m-3]
  integer  :: IXTYP(2,JY,JX)                    !surface litter type:1=plant,2=manure
  real(r8) :: XCORP(JY,JX)                      !factor for surface litter incorporation and soil mixing
  real(r8) :: FLWRM(60,JY,JX)                   !water transfer between soil surface and surface litter, [g d-2 t-1]
  real(r8) :: FLQRM(60,JY,JX)                   !meltwater flux into surface litter, [m3 d-2 h-1]
  real(r8) :: CVRD(JY,JX)                       !fraction of soil surface covered by surface litter, [-]
  real(r8) :: HFLWR(JY,JX)                      !net heat transfer to surface litter, [MJ d-2 t-1]
  real(r8) :: VOLR(JY,JX)                       !surface litter volume, [m3 d-2]
  real(r8) :: VHCPRX(JY,JX)                     !surface litter heat capacity from previous time step, [MJ d-2 K-1]
  real(r8) :: VOLWRX(JY,JX)                     !surface litter water holding capacity, [m3 d-2]
  real(r8) :: FLWR(JY,JX)                       !net water transfer to surface litter, [MJ d-2 t-1]
  real(r8) :: THAWR(JY,JX)                      !freeze (-ve) - thaw (+ve) in surface litter, [m3 d-2 h-1]
  real(r8) :: HTHAWR(JY,JX)                     !latent heat of freeze (-ve) - thaw (+ve) in surface litter, [MJ d-2 h-1]
  real(r8) :: FLQRQ(JY,JX)                      !precipitation flux into surface litter, [m3 d-2 h-1]
  real(r8) :: FLQRI(JY,JX)                      !irrigation flux into surface litter, [m3 d-2 h-1]
end module SurfLitterDataType
