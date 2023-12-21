module PlantMod
  use data_kind_mod   , only : r8 => DAT_KIND_R8
  use grosubsMod      , only : grosubs
  use HfuncsMod       , only : hfuncs
  use UptakesMod      , only : uptakes
  use PlantDisturbMod , only : PrepLandscapeGrazing
  use EcoSimSumDataType
  use PlantAPIData  
  use PlantAPI
  use ExtractsMod
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
  integer :: NY,NX

333   FORMAT(A8)

  call PrepLandscapeGrazing(I,J,NHW,NHE,NVN,NVS)

  plt_site%PlantElemntStoreLandscape(:)=PlantElemntStoreLandscape(:)
  DO NX=NHW,NHE
    DO NY=NVN,NVS
!
      call  PlantAPISend(I,J,NY,NX)
!   UPDATE PLANT PHENOLOGY IN 'HFUNC'
!
    !if(lverb)WRITE(*,333)'HFUNC'

      CALL HFUNCs(I,J)

      CALL UPTAKES(I,J)

      CALL GROSUBs(I,J)

      CALL EXTRACTs(I,J)

      call PlantAPIRecv(I,J,NY,NX)
    ENDDO
  ENDDO
  PlantElemntStoreLandscape(:)=plt_site%PlantElemntStoreLandscape(:)


  end subroutine PlantModel
!------------------------------------------------------------------------------------------

  subroutine PlantCanopyRadsModel(I,J,NY,NX,DPTH0)
  use CanopyCondsMod
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(in) :: DPTH0

  call PlantAPICanMSend(NY,NX)

  call CanopyConditionModel(I,J,DPTH0)

  call PlantAPICanMRecv(NY,NX)

  end subroutine PlantCanopyRadsModel

end module PlantMod
