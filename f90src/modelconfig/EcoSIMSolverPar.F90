module EcoSIMSolverPar

! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__
  public

  real(r8) :: XNPX
  real(r8) :: XNPA
  real(r8) :: XNPB
  real(r8) :: XNPC
  real(r8) :: XNPD
  real(r8) :: XNPY
  real(r8) :: XNPZ
  real(r8) :: XNPQ
  real(r8) :: XNPH
  real(r8) :: XNPT
  real(r8) :: XNPG
  real(r8) :: XNPR
  real(r8) :: XNPS
  integer :: NPH             !number of model cycles per hour for heat and water fluxes
  integer :: NPT             !number of model cycles per hour for gas fluxes within each NPH
  integer :: NPG             !number of model cycles per hour for gas fluxes
  integer :: NPR             !number of cycles per time step used to calculate hat and water transfer in surface litter
  integer :: NPS             !number of cycles for solving snowpack heat and water fluxes

end module EcoSIMSolverPar
