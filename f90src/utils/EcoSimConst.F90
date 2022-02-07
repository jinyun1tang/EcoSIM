module EcosimConst
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  real(r8), parameter :: Cpw=4.19_r8
  real(r8), parameter :: cpi=1.9274_r8
  real(r8), parameter :: cpo=2.496E-06_r8
  real(r8), parameter :: cps=2.095_r8
  real(r8), parameter :: TFice=273.15_r8
  real(r8), parameter :: TC2K=273.15_r8
  real(r8), parameter :: Tref=273.15_r8
end module EcosimConst
