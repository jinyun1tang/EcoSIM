module PlantMod
  use data_kind_mod,   only: r8 => DAT_KIND_R8
  use grosubsMod,      only: GrowPlants
  use PlantPhenolMod,  only: hfuncs
  use UptakesMod,      only: RootUptakes
  use PlantDisturbMod, only: PrepLandscapeGrazing
  use EcoSimSumDataType
  use EcoSIMCtrlMod, only : lverb
  use PlantAPIData  
  use PlantAPI
  use ExtractsMod
  use PlantMgmtDataType, only : NP
  use PlantBalMod
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: PlantModel
  public :: PlantCanopyRadsModel

  contains

  subroutine PlantModel(I,J,NHW,NHE,NVN,NVS)


  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: t1
  integer :: NY,NX,NZ

333   FORMAT(A8)

  call PrepLandscapeGrazing(I,J,NHW,NHE,NVN,NVS)

  plt_site%PlantElemntStoreLandscape(:)=PlantElemntStoreLandscape(:)
  DO NX=NHW,NHE
    DO NY=NVN,NVS

      call  PlantAPISend(I,J,NY,NX)

      call EnterPlantBalance(I,J,NP(NY,NX))

!   UPDATE PLANT PHENOLOGY IN 'HFUNC'
!     zero out plant hourly fluxes
      call ZeroGrosub()  

      if(lverb)WRITE(*,333)'HFUNC'
      !phenological update, determine living/active branches
!     DO NZ=1,NP(NY,NX)
!       call SumPlantBiom(I,J,NZ,'bfHFUNCS')
!     ENDDO
      
      CALL HFUNCs(I,J)

!      DO NZ=1,NP(NY,NX)
!        call SumPlantBiom(I,J,NZ,'bfUPTAKES')
!      ENDDO
      !predict uptake fluxes of nutrients and O2
      if(lverb)write(*,*)'uptake'
      CALL ROOTUPTAKES(I,J)
!      DO NZ=1,NP(NY,NX)
!        call SumPlantBiom(I,J,NZ,'bfGROSUBS')
!      ENDDO
      !do growth of active branches
      if(lverb)write(*,*)'grosub'
      CALL GROWPLANTS(I,J)
      if(lverb)write(*,*)'EXTRACT'
      !aggregate varaibles
      CALL EXTRACTs(I,J)

      call ExitPlantBalance(I,J,NP(NY,NX))

      call PlantAPIRecv(I,J,NY,NX)
    ENDDO
  ENDDO
  PlantElemntStoreLandscape(:)=plt_site%PlantElemntStoreLandscape(:)


  end subroutine PlantModel
!------------------------------------------------------------------------------------------

  subroutine PlantCanopyRadsModel(I,J,NY,NX,DPTH0)
  use CanopyCondsMod
  use CanopyDataType, only : CanopyLeafArea_lpft
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(in) :: DPTH0

  call PlantAPICanMSend(NY,NX)

  call CanopyConditionModel(I,J,DPTH0)

  call PlantAPICanMRecv(NY,NX)

  end subroutine PlantCanopyRadsModel

end module PlantMod
