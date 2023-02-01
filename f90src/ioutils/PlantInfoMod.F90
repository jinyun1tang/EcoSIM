module PlantInfoMod
!
! DESCRIPTION
! code to read plant information
  use EcoSIMCtrlDataType, only : lverb,IYRR
implicit none

  save
  character(len=*),private, parameter :: mod_filename = __FILE__

  public :: ReadPlantInfo
  contains

  subroutine ReadPlantInfo(NA,ND,NT,NE,NTX,NEX,NF,NFX,NTZ,NTZX,NHW,NHE,NVN,NVS)

  use readqmod     , only : readq
  use routqMod     , only : routq

  implicit none
  integer, intent(in) :: NT   !
  integer, intent(in) :: NE   !
  integer, intent(in) :: NFX  !
  integer, intent(in) :: NTZX  !
  integer, intent(in) :: NTX  !
  integer, intent(in) :: NEX  !
  integer, intent(in) :: NF,NTZ
  integer, intent(in) :: NHW,NHE,NVN,NVS     !simulation domain specification
  integer, intent(in) :: NA(1:NEX),ND(1:NEX)

! RECOVER PLANT SPECIES DISTRIBUTION IN 'ROUTQ'
!
  if(lverb)WRITE(*,333)'ROUTQ'
  print*,'ROUTQ',IYRR
  CALL ROUTQ(NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)
!
!   READ INPUT DATA FOR PLANT SPECIES AND MANAGEMENT IN 'READQ'
!   AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTP'
!
  if(lverb)WRITE(*,333)'READQ'
  CALL READQ(NA,ND,NT,NE,NTX,NEX,NF,NFX,NTZ,NTZX,NHW,NHE,NVN,NVS)

333   FORMAT(A8)

  end subroutine ReadPlantInfo

end module PlantInfoMod
