module EcosimConst
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  real(r8), parameter :: Cpw=4.19_r8           !heat capacity for water
  real(r8), parameter :: cpi=1.9274_r8         !heat capacity for ice
  real(r8), parameter :: cpo=2.496E-06_r8      !heat capacity for organic matter
  real(r8), parameter :: cps=2.095_r8          !heat capacity for fresh snow
  real(r8), parameter :: TFice=273.15_r8       !frozen temperature [K]
  real(r8), parameter :: TC2K=273.15_r8        !temperature for converting celcius to Kelvin
  real(r8), parameter :: Tref=273.15_r8        !reference temperature for atmospheric variables [K]
  real(r8), parameter :: VHCPWMin=2.095E-04_r8 !minimum heat capacities for solving
                                               !snowpack layered water and heat fluxes, [MJ/K]
  real(r8), parameter :: VHCPRMin=4.190E-05_r8 !minimum heat capacities for solving
                                               !surface litter water and heat fluxes, [MJ/K]
  real(r8), parameter :: VHCPNMin=4.190E-03_r8 !minimum heat capacities for solving
                                               !soil water and heat fluxes, [MJ/K]
  real(r8), parameter :: PICON=3.14159265358979323846_r8
  real(r8), parameter :: PICON2=PICON*0.5_r8
  real(r8), PARAMETER :: PSIPS=-0.5E-03_r8
  real(r8), parameter :: RDN=57.29577951_r8

end module EcosimConst
