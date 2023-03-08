module PlantInfoMod
!
! DESCRIPTION
! code to read plant information
  use EcoSIMCtrlDataType, only : IYRR
  use EcoSIMCtrlMod, only : lverb
implicit none

  save
  character(len=*),private, parameter :: mod_filename = __FILE__

  public :: ReadPlantInfo
  contains

  subroutine ReadPlantInfo(yearc,yeari,NE,NEX,NHW,NHE,NVN,NVS)

  use readqmod     , only : readq
  use routqMod     , only : routq

  implicit none
  integer, intent(in) :: yearc, yeari
  integer, intent(in) :: NE   !
  integer, intent(in) :: NEX  !
  integer, intent(in) :: NHW,NHE,NVN,NVS     !simulation domain specification

! RECOVER PLANT SPECIES DISTRIBUTION IN 'ROUTQ'
!
  if(lverb)WRITE(*,333)'ROUTQ'
  CALL ROUTQ(yearc,yeari,NE,NEX,NHW,NHE,NVN,NVS)
!
!   READ INPUT DATA FOR PLANT SPECIES AND MANAGEMENT IN 'READQ'
!   AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTP'
!
  if(lverb)WRITE(*,333)'READQ'
  CALL READQ(yearc,yeari,NE,NEX,NHW,NHE,NVN,NVS)

333   FORMAT(A8)

  end subroutine ReadPlantInfo

end module PlantInfoMod
