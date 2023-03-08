module fileUtil
!!
! subroutines for file open with error check
  use abortutils, only : endrun
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__

  public :: open_safe
  public :: check_read
  public :: remove_filename_extension
  public :: file_exists
  integer, parameter :: ecosim_filename_length=128
  integer, parameter :: stdout=6
  integer, parameter :: iulog=6
  integer, parameter :: error_errmsg_len=256
  integer, parameter :: ecosim_string_length_long=256
  integer , public, parameter :: var_flux_type =1
  integer , public, parameter :: var_state_type=2
  integer , public, parameter :: ecosim_namelist_buffer_size = 4096
  logical, save :: continue_run = .false.
  contains
!------------------------------------------------------------------------------------------

  function file_exists(filename) result(res)
!
!! DESCRIPTION
! check existence of a file
  implicit none
  character(len=*),intent(in) :: filename
  logical                     :: res

  ! Check if the file exists
  inquire( file=trim(filename), exist=res )
  end function
!------------------------------------------------------------------------------------------

  subroutine open_safe(lun,prefix,fname,status,location,lineno,lverb)
!
! Description
! safely open a file for IO
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

  call getfilenamef(prefix,fname,pathfile)

  if(trim(status)/='UNKNOWN' .and. .not.file_exists(pathfile))then
    call endrun(msg='Fail to find file '//trim(pathfile)//' in ' &
      //mod_filename,line=__LINE__)
  endif

  if(present(lverb))then
    if(lverb)write(*,*)'open file ',trim(pathfile)
  endif
  OPEN(lun,FILE=trim(pathfile),STATUS=status,iostat=ierr)
  if(ierr/=0)then
    call endrun(msg="error in "//location//" while reading file " &
      //TRIM(pathfile), line=lineno)
  endif
  end subroutine open_safe

!------------------------------------------------------------------------------------------

  subroutine check_read(ios,nelm,ifile,lineno,modfile)
  implicit none
  integer, intent(in) :: ios
  integer, intent(in) :: nelm
  character(len=*),intent(in) :: ifile
  integer, intent(in) :: lineno
  character(len=*), intent(in) :: modfile

  if(ios>0)then
    write(*,*)'input error'
    call endrun('read error when reading '//trim(ifile)//' in'//trim(modfile),lineno)
  elseif(ios<0)then
    write(*,*)'number of inputs is less than ',nelm
    call endrun('read error when reading '//trim(ifile)//' in'//trim(modfile),lineno)
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

!------------------------------------------------------------------------------------------

  subroutine namelist_to_buffer(namelist_filename, namelist_buffer)
  !DESCRIPTION
  !read in namelist
  !USES
  implicit none
  character(len=*)                         , intent(in)  :: namelist_filename
  character(len=ecosim_namelist_buffer_size) , intent(out) :: namelist_buffer

  character(len=*), parameter                            :: subname = 'namelist_to_buffer'
  character(len=ecosim_string_length_long)                 :: ioerror_msg
  integer :: nml_unit, nml_error

  nml_unit = 16

  ! read the namelist file into a buffer.
  open(unit=nml_unit, file=trim(namelist_filename), action='read', access='stream', &
       form='unformatted', iostat=nml_error)
  if (nml_error == 0) then
     read(unit=nml_unit, iostat=nml_error, iomsg=ioerror_msg) namelist_buffer

     ! we should always reach the EOF to capture the entire file...
     if (.not. is_iostat_end(nml_error)) then
        write(stdout, '(a, a, i8)') subname, &
             ": IO ERROR reading namelist file into buffer: ", nml_error
        write(stdout, '(a)') ioerror_msg
        call abort()
     else
        write(stdout, '(a, a, a)') "Read '", trim(namelist_filename), "' until EOF."
     end if

     write(stdout, '(a, a, i7, a)') subname, ": Read buffer of ", &
          len_trim(namelist_buffer), " characters."
     write(stdout, '(a)') "------------------------------"
     write(stdout, '(a)') "  If it looks like part of the namelist is missing, "
     write(stdout, '(a)') "  compare the number of characters read to the actual "
     write(stdout, '(a,a,a)') "  size of your file ($ wc -c ", trim(namelist_filename), ") and increase "
     write(stdout, '(a)') "  the buffer size if necessary."
     write(stdout, '(a)') "------------------------------"
     write(stdout, '(a)') trim(namelist_buffer)
     write(stdout, '(a)') "------------------------------"
  else
     write(stdout, '(a, a, i8, a, a)') subname, ": IO ERROR ", nml_error, &
          " opening namelist file : ", trim(namelist_filename)
     call abort()
  end if
  close(nml_unit)
  end subroutine namelist_to_buffer

!------------------------------------------------------------------------------------------
  subroutine getfilenamef(s1,s2,s3)
  implicit none
  character(len=*), intent(in) :: s1
  character(len=*), intent(in) :: s2
  character(len=*), intent(out) :: s3

  integer :: k,k1,l3
  
  l3=len(s3)
  k=1
  do while(.true.)
    if(s1(k:k)/=' '.and. ichar(s1(k:k))/=0 .and. k<=l3)then
       s3(k:k)=s1(k:k)
    else
       exit
    endif
    k=k+1
  enddo
  if(k>l3)then
    call endrun("Not sufficient memory allocated for string s3 in " &
      //trim(mod_filename),__LINE__)
  endif
  k1=1
  do while(.true.)
    if(s2(k1:k1)/=' '.and. ichar(s2(k1:k1))/=0 .and. k<=l3)then
       s3(k:k)=s2(k1:k1)
    else
       exit
    endif
    k=k+1
    k1=k1+1
  enddo
  s3(k:)=' '
  end subroutine getfilenamef
  
!------------------------------------------------------------------------------------------

  logical function isprtablesymb(c)
  implicit none
  character(len=1), intent(in) :: c
  
  isprtablesymb=ICHAR(c) >= ICHAR(' ') .AND. ICHAR(c) <= ICHAR('~')
  end function isprtablesymb
end module fileUtil
