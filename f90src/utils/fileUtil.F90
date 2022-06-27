module fileUtil
!!
! subroutines for file open with error check
  use abortutils, only : endrun
  implicit none

  public :: open_safe
  public :: check_read
  public :: remove_filename_extension
  integer, parameter :: ecosim_filename_length=128
  integer, parameter :: stdout=6
  integer, parameter :: iulog=6
  integer, parameter :: ecosim_string_length_long=256
  integer , public,  parameter :: var_flux_type =1
  integer , public,  parameter :: var_state_type=2
  logical, save :: continue_run = .false.
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

!------------------------------------------------------------------------------------------


  function remove_filename_extension(filename) result(basename)
    !
    ! Remove the extension from a file name to get a base filename.
    !
    ! We start at the end of the filename and assume that the extension
    ! is marked by a period.
    !
    implicit none

    character(len=ecosim_filename_length), intent(in) :: filename
    character(len=ecosim_filename_length) :: basename
    integer :: ext_index

    ext_index = scan(filename, '.', .true.)
    if (ext_index == 0) then
       ! no period marking an extension...
       ext_index = len(trim(filename)) + 1

    end if
    basename = filename(1:ext_index-1)
  end function remove_filename_extension

end module fileUtil
