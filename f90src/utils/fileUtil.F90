module fileUtil
!!
! subroutines for file open with error check
  use abortutils, only : endrun
  implicit none

  public :: open_safe
  public :: check_read
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

!------------------------------------------------------------------------------------------

  subroutine check_read(ios,nelm,lineno,modfile)
  implicit none
  integer, intent(in) :: ios
  integer, intent(in) :: nelm
  integer, intent(in) :: lineno
  character(len=*), intent(in) :: modfile

  if(ios>0)then
    write(*,*)'input error'
    call endrun('read error in '//trim(modfile),lineno)
  elseif(ios<0)then
    write(*,*)'number of inputs is less than ',nelm
    call endrun('read error in '//trim(modfile),lineno)
  endif
  end subroutine check_read
end module fileUtil
