module ATSEcoSIMAdvanceMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilWaterDataType
  use SharedDataMod
  use GridDataType
  use GridConsts
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
  integer :: NY,NX,L,NHW,NHE,NVN,NVS, I, J, M
  integer, intent(in) :: NYS
  real(r8) :: YSIN(JSA),YCOS(JSA),YAZI(JSA)
  !real(r8), dimension(:,:), allocatable :: ResistanceLitRLay
  !real(r8), dimension(:,:), allocatable :: KSatReductByRainKineticEnergyS
  !real(r8), dimension(:,:), allocatable :: HeatFlux2Ground
  !real(r8), dimension(:,:), allocatable :: Qinfl2MicP
  !real(r8), dimension(:,:), allocatable :: TopLayWatVol
  
  !real(r8), dimension(:,:),intent(inout) :: ResistanceLitRLay
  !REAL(R8), dimension(:,:),INTENT(OUT) :: KSatReductByRainKineticEnergy
  !real(r8), dimension(:,:),intent(out) :: HeatFlux2Ground 
  !real(r8),dimension(:,:),intent(inout) :: TopLayWatVol     
  real(r8) :: ResistanceLitRLay(JY,JX)
 
  !real(r8), dimension(:), allocatable :: HeatFlux2Ground

  NHW=1;NHE=1;NVN=1;NVS=NYS

  write(*,*) "In Run Surface Balance"

  call SetMeshATS(NHW,NVN,NHE,NVS)
  
  write(*,*) "finish set mesh"

  NX=1
  !write(*,*) "HeatFlux2Ground"
  !HeatFlux2Ground=a_PORO
  !write(*,*) "KSatREduct"
  !KSatReductByRainKineticEnergyS=a_PORO
  !write(*,*) "Qinfl2MicP"
  !Qinfl2MicP=a_PORO
  !write(*,*) "ResistanceLitRLay"
  !ResistanceLitRLay=a_PORO
  !write(*,*) "TopLayWatVol" 
  !TopLayWatVol=a_PORO

  !What are I and J are these a loop?
  write(*,*) "Running StageSurfacePhysModel"
  call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)

  write(*,*) "Done; Running RunSurfacePhysModel"

  !call RunSurfacePhysModel(M,NHE,NHW,NVS,NVN,ResistanceLitRLay,&
  !    KSatReductByRainKineticEnergyS,HeatFlux2Ground,TopLayWatVol)

  write(*,*) "Finished Subroutine RunEcoSIMSurfaceBalance"
  end subroutine RunEcoSIMSurfaceBalance

end module ATSEcoSIMAdvanceMod
