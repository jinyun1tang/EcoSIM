module EcosimConst
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  character(len=*),private, parameter :: mod_filename = __FILE__
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
  real(r8), parameter :: VAP=2465.0_r8         !latent heat of vaporization of water, kJ/kg
  real(r8), parameter :: VAPS=2834.0_r8        !sublimation of water, kJ/kg
  real(r8), parameter :: TCNDG=8.1E-03_r8      !assumed thermal conductivity below lower soil boundary,[MJ m-1 h-1 K-1]
  real(r8), parameter :: RGAS=8.3143_r8        !universal gas constant, J/mole.oK
  real(r8) :: FCI=0.05_r8                      !field capacity of ice
  real(r8) :: WPI=0.025_r8                     !wilting point of ice
  real(r8) :: POROQ=0.66_r8                    !soil porosity ^ 2/3
  real(r8) :: FORGC=0.1E+06_r8                 !minimum SOC for combustion 	[g Mg-1]
  real(r8) :: FVLWB=1.0_r8                     !maximum soil water content for combustion	[m3 m-3]
  real(r8) :: FCH4F=0.01_r8                    !fraction of combusted C released as CH4
  real(r8) :: PSIHY=-2500.0_r8                 !hygroscopic water potential [MPa]
  real(r8) :: OXKM=0.080_r8                    !Km for heterotrophic O2 uptake
  real(r8) :: THETX=1.0E-03_r8                 !minimum air-filled porosity for gas transfer	[m3 m-3]
  real(r8), parameter :: THETPI=0.00_r8        !air content of ice	-
  real(r8), parameter :: DENSI=0.92_r8-THETPI  !ice density
  real(r8), parameter :: ZW=0.01_r8            !snowpack surface roughness [m]
  real(r8), parameter :: Catomw=12._r8         !C 12 molar mass,[g/mol]
  real(r8), parameter :: Natomw=14._r8         !N 14 molar mass [g/mol]
  real(r8), parameter :: Patomw=31._r8         !P 31 molar mass [g/mol]
  integer :: ICOR(12)=(/1,-1,0,0,1,1,2,3,3,4,4,5/)
  real(r8), parameter :: TWILGT=0.06976_r8
  real(r8), parameter :: MWC2Soil=1.82E-06_r8  !convert organic carbon (g C/d2) to soil mass (Mg/d2)
end module EcosimConst
