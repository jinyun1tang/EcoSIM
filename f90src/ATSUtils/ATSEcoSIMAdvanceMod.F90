module ATSEcoSIMAdvanceMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilWaterDataType
  use SharedDataMod
  use GridDataType
  use SOMDataType
  USE SoilPhysDataType
  use LandSurfDataType
  use ClimForcDataType
  use SoilPropertyDataType
implicit none
  character(len=*), private, parameter :: mod_filename=&
  __FILE__
  public :: RunEcoSIMSurfaceBalance
  contains

  subroutine RunEcoSIMSurfaceBalance(NYS)
  use EcoSimConst
  use GridMod           , only : SetMeshATS
  use SurfPhysMod       , only : RunSurfacePhysModel, StageSurfacePhysModel
  implicit none
  integer :: NY,NX,L,NHW,NHE,NVN,NVS
  integer, intent(in) :: NYS
  real(r8) :: YSIN(JSA),YCOS(JSA),YAZI(JSA)

  NHW=1;NHE=1;NVN=1;NVS=NYS

  write(*,*) "In Run Surface Balance"

  call SetMeshATS(NHW,NVN,NHE,NVS)

  NX=1

  !What are I and J are these a loop?
  write(*,*) "Running StageSurfacePhysModel"
  call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)

  write(*,*) "Done; Running RunSurfacePhysModel"

  call RunSurfacePhysModel(M,NHE,NHW,NVS,NVN,ResistanceLitRLay,&
      KSatReductByRainKineticEnergyS,TopLayWatVol,HeatFlux2Ground,Qinfl2MicP)

  write(*,*) "Finished Subroutine RunEcoSIMSurfaceBalance"
  end subroutine RunEcoSIMSurfaceBalance

end module ATSEcoSIMAdvanceMod
