module InitEcoSIM
!
!! DESCRIPTION
! initialize the data structure for EcoSIM
  use EcoSIMCtrlMod
  use EcoSiMParDataMod    , only : micpar,pltpar
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: InitModules
  contains

  subroutine InitModules(NOMicrobeGuilds)

  use PlantDisturbsMod, only: InitPlantDisturbance
  use UptakesMod,       only: InitUptake
  use WatsubMod,        only: initWatsub
  use SoilBGCNLayMod,   only: InitNitro
  use RedistMod,        only: InitRedist
  use ErosionMod,       only: InitErosion
  use Hour1Mod,         only: InitHour1
  use HistDataType,     only: hist_ecosim
  use InitAllocMod,     only: InitAlloc
  use UnitMod,          only: units
  use TranspNoSaltMod,  only: InitTranspNoSalt
  use GridConsts
  implicit  none

  integer                 , intent(in) :: NOMicrobeGuilds   !number of microbial guilds per group


  call InitAlloc(NOMicrobeGuilds)

  call units%Initailize()

  call InitPlantDisturbance

  call InitUptake

  call initNitro

  call InitRedist

  call InitErosion

  call InitHour1(micpar%NumOfLitrCmplxs)

  call InitTranspNoSalt

  call hist_ecosim%Init(bounds)

  end subroutine InitModules


end module InitEcoSIM
