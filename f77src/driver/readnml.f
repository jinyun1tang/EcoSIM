      subroutine readnml(nmlfile,runfile, case_name,
     2prefix,do_rgres,LYRG,lverb)

      !Description
      !read control namelist
      implicit none
      character(len=*), intent(in) :: nmlfile
      character(len=80), intent(out) :: runfile
      character(len=36)    , intent(out) :: case_name
      character(len=80)    , intent(out) :: prefix
      logical              , intent(out) :: do_rgres
      integer              , intent(out) :: LYRG
      logical              , intent(out) :: lverb
      integer, parameter :: stdout = 6

      logical :: do_regression_test
      integer :: num_of_simdays
      logical :: lverbose
      namelist / ecosys/case_name, prefix, runfile, do_regression_test,
     2num_of_simdays,lverbose

      !local variables
      character(len=256) :: ioerror_msg
      integer :: rc, fu
      integer :: nml_error

      num_of_simdays=-1
      do_regression_test=.false.
      lverbose=.false.
      inquire (file=nmlfile, iostat=rc)
      if (rc /= 0) then
        write (stdout, '(3a)') 'Error: input file ', trim(nmlfile), 
     &' does not exist.'
        call endrun('stopped in readnml', 25)
      end if

      open (action='read', file=nmlfile, iostat=rc, newunit=fu)
      if (rc /= 0) then
        write (stdout, '(2a)') 'Error openning input file "',
     &trim(nmlfile)
        call endrun('stopped in readnml', 32)
      end if

      read(unit=fu, nml=ecosys, iostat=nml_error, iomsg=ioerror_msg)
      if (nml_error /= 0) then
         write(stdout,'(a)')"ERROR reading ecosys namelist "
         call endrun('stopped in readnml', 38)
      end if

      close(fu)
      if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout,ecosys)
       write(stdout, *)
       write(stdout, *) '--------------------'
      endif     
      do_rgres=do_regression_test 
      LYRG=num_of_simdays
      lverb=lverbose
      end subroutine readnml 
