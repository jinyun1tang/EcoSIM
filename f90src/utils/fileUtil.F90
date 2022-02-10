module fileUtil
!!
! subroutines for file open with error check
  use abortutils, only : endrun
  implicit none

  public :: open_safe
  contains

  subroutine open_safe(lun,prefix,fname,status,location,lineno,lverb)

  implicit none
  integer, intent(in) :: lun
  character(len=*), intent(in) :: prefix
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: status
  character(len=*), intent(in) :: location
  integer, intent(in) :: lineno
  logical, optional, intent(in) :: lverb
  character(len=2560) :: pathfile
  integer :: ierr

  call getfilename(prefix,fname,pathfile)
  if(present(lverb))then
    if(lverb)write(*,*)'open file ',trim(pathfile)
  endif
  OPEN(lun,FILE=pathfile,STATUS=status,iostat=ierr)
  if(ierr/=0)then
    call endrun(msg="error in "//location//" while reading file " &
      //TRIM(pathfile), line=lineno)
  endif
  end subroutine open_safe
end module fileUtil
