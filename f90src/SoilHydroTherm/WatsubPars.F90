module WatsubPars

! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  public
  save
  character(len=*),private, parameter :: mod_filename = __FILE__
!
! EMMS,EMMW,EMMR=emissivities of surface soil, snow and litter
! RACX,RARX=minimum boundary layer resistances of canopy,litter (h m-1)
! RZ=minimum resistance to evaporation of surface water (h m-1)
! RAM=minimum boundary layer resistance (h m-1)
! DPTHSX=minimum snowpack depth for full cover (m)
! Z1S,Z2SW,Z2SD,Z3SX=parameters for air-water gas transfers in soil
! Z1R,Z2RW,Z2RD,Z3RX=parameters for air-water gas transfers in litter
!
!
! Parameters for calculating convective effects on heat transfer
! in porous media (air and water)
! VISCW,VISCA=water,air viscosity (Mg m-1 s)
! CP[w,i,s]=heat capacity of water, ice, snow, kJ/kg/K
! CPo = heat capacity of organic matter, kJ/g/K
!
! FVOLAH=parameter for clay effect on macropore volume
! DTHETW=difference between saturation and effective saturation
! HCNDRR=saturated hydraulic conductivity of surface litter
! FENGYP=rate constant for restoring surface Ksat
!

  real(r8) :: EMMS
  real(r8) :: EMMW
  real(r8) :: EMMR
  real(r8) :: RACX
  real(r8) :: RARX
  real(r8) :: RZ
  real(r8) :: RAM
  real(r8) :: DPTHSX
  real(r8) :: Z1S
  real(r8) :: Z2SW
  real(r8) :: Z2SD
  real(r8) :: Z3SX
  real(r8) :: Z1R
  real(r8) :: Z2RW
  real(r8) :: Z2RD
  real(r8) :: Z3R
  real(r8) :: VISCW
  real(r8) :: VISCA
  real(r8) :: DIFFW
  real(r8) :: DIFFA
  real(r8) :: EXPNW
  real(r8) :: EXPNA
  real(r8) :: GRAV
  real(r8) :: RYLXW
  real(r8) :: RYLXA
  real(r8) :: PRNTW
  real(r8) :: PRNTA
  real(r8) :: DNUSW
  real(r8) :: DNUSA
  real(r8) :: TRBW
  real(r8) :: TRBA
  real(r8) :: FVOLAH
  real(r8) :: DTHETW
  real(r8) :: HCNDRR
  real(r8) :: FENGYP

  contains
  subroutine InitWatsubPars
  implicit none

  EMMS=0.97_r8
  EMMW=0.97_r8
  EMMR=0.97_r8
  RACX=0.0139_r8
  RARX=0.0139_r8
  RZ=0.0139_r8
  RAM=1.39E-03_r8   !this value differs from that in Hour1Mod.F90
  DPTHSX=0.075_r8
  Z1S=0.010_r8
  Z2SW=12.0_r8
  Z2SD=12.0_r8
  Z3SX=0.50_r8
  Z1R=0.01_r8
  Z2RW=12.0_r8
  Z2RD=12.0_r8
  Z3R=0.50_r8
  VISCW=1.0E-06_r8
  VISCA=2.0E-08_r8
  DIFFW=1.45E-07_r8
  DIFFA=2.01E-05_r8
  EXPNW=2.07E-04_r8
  EXPNA=3.66E-03_r8
  GRAV=9.8_r8
  RYLXW=GRAV*EXPNW/(VISCW*DIFFW)
  RYLXA=GRAV*EXPNA/(VISCA*DIFFA)
  PRNTW=VISCW/DIFFW
  PRNTA=VISCA/DIFFA
  DNUSW=(1.0_r8+(0.492_r8/PRNTW)**0.5625_r8)**0.4444_r8
  DNUSA=(1.0_r8+(0.492_r8/PRNTA)**0.5625_r8)**0.4444_r8
  TRBW=0.375_r8
  TRBA=0.000_r8
  FVOLAH=0.0_r8
  DTHETW=1.0E-06_r8
  HCNDRR=25.0_r8
  FENGYP=1.0E-03_r8
  end subroutine InitWatsubPars
end module WatsubPars
