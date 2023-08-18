module EcoSiMParDataMod
 use GrosubPars, only : plant_bgc_par_type
 use MicBGCPars, only : MicParType
implicit none

  type(plant_bgc_par_type), target, public :: pltpar
  type(MicParType), target, public :: micpar
end module EcoSiMParDataMod
