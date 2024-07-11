module EcoSIMSolverPar

! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  public

  real(r8) :: dts_wat        !1/hour, time step for water transport calculation
  real(r8) :: dts_sno        !1/hour, time step for snow process
  real(r8) :: XNPB
  real(r8) :: dt_watvap      !1/hour, time step for water vapor fluxes
  real(r8) :: XNPD
  real(r8) :: dts_HeatWatTP   !1/hour, time step for heat, water, solute transport
  real(r8) :: dt_GasCyc
  real(r8) :: dts_gas         !1/hour, time step for gas flx update
  real(r8) :: XNPR
  real(r8) :: XNPS
  integer :: NPH             !number of model cycles per hour for heat and water fluxes
  integer :: NPT             !number of model cycles per hour for gas fluxes within each NPH
  integer :: NPG             !number of model cycles per hour for gas fluxes
  integer :: NPR             !number of cycles per time step used to calculate hat and water transfer in surface litter
  integer :: NPS             !number of cycles for solving snowpack heat and water fluxes

end module EcoSIMSolverPar
