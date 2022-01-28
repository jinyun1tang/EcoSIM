module fileUtil
!!
! subroutines for file open with error check
  use abortutils, only : endrun
  implicit none

  public :: open_safe
  contains

  subroutine open_safe(lun,prefix,fname,status,location,lineno)

  implicit none
  integer, intent(in) :: lun
  character(len=*), intent(in) :: prefix
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: status
  character(len=*), intent(in) :: location
  integer, intent(in) :: lineno

  character(len=2560) :: pathfile
  integer :: ierr

  call getfilename(prefix,fname,pathfile)
  OPEN(lun,FILE=pathfile,STATUS=status,iostat=ierr)
  if(ierr/=0)then
    call endrun(msg="error in "//location//" while reading file " &
      //TRIM(pathfile), line=lineno)
  endif
  end subroutine open_safe
end module fileUtil
