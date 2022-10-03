 subroutine readnamelist(nmlfile,runfile, case_name, &
  prefix,do_rgres,LYRG,lverb,nmicbguilds)
!!
! Description:
! read control namelist
  use abortutils, only : endrun
  use GridConsts, only : JX,JY,JZ,JH,JV,JD
  use ForcWriterMod, only : bgc_forc_conf,do_bgcforc_write
  implicit none
      character(len=*), intent(in) :: nmlfile
      character(len=80), intent(out) :: runfile
      character(len=36)    , intent(out) :: case_name
      character(len=80)    , intent(out) :: prefix
      logical              , intent(out) :: do_rgres
      integer              , intent(out) :: LYRG
      logical              , intent(out) :: lverb
      integer              , intent(out) :: nmicbguilds
  character(len=*), parameter :: mod_filename = __FILE__
      integer, parameter :: stdout = 6

      logical :: do_regression_test
      integer :: num_of_simdays
      logical :: lverbose
      integer :: num_microbial_guilds
      integer :: do_doy,do_year,do_layer
      character(len=64) :: bgc_fname
      namelist /ecosys/case_name, prefix, runfile, do_regression_test, &
      num_of_simdays,lverbose,num_microbial_guilds,JX,JY,JZ

      logical :: laddband
      namelist /bbgcforc/do_bgcforc_write,do_year,do_doy,laddband,do_layer,&
        bgc_fname

      !local variables
      character(len=256) :: ioerror_msg
      integer :: rc, fu
      integer :: nml_error

      JX=2
      JY=2
      JZ=14
      num_of_simdays=-1
      do_year=-1
      do_doy=0
      do_layer=1
      laddband=.false.
      do_regression_test=.false.
      lverbose=.false.
      num_microbial_guilds=1
      do_bgcforc_write=.false.
      bgc_fname='bbforc.nc'
      inquire (file=nmlfile, iostat=rc)
      if (rc /= 0) then
        write (stdout, '(3a)') 'Error: input file ', trim(nmlfile), &
      ' does not exist.'
        call endrun('stopped in readnml', 25)
      end if

      open (action='read', file=nmlfile, iostat=rc, newunit=fu)
      if (rc /= 0) then
        write (stdout, '(2a)') 'Error openning input file "', &
      trim(nmlfile)
        call endrun('stopped in readnml', 32)
      end if

      read(unit=fu, nml=ecosys, iostat=nml_error, iomsg=ioerror_msg)
      if (nml_error /= 0) then
         write(stdout,'(a)')"ERROR reading ecosys namelist "
         call endrun('stopped in readnml', 38)
      end if

      read(unit=fu, nml=bbgcforc, iostat=nml_error, iomsg=ioerror_msg)
      if (nml_error /= 0) then
         write(stdout,'(a)')"ERROR reading bbgcforc namelist "
         call endrun('stopped in readnml', 38)
      end if

  close(fu)
  if (.true.) then
    write(stdout, *)
    write(stdout, *) '--------------------'
    write(stdout,ecosys)
    write(stdout, *)
    write(stdout, *) '--------------------'
    write(stdout,bbgcforc)
    write(stdout, *)
    write(stdout, *) '--------------------'

  endif
  if(do_bgcforc_write)then
    bgc_forc_conf%doy =do_doy
    bgc_forc_conf%year=do_year
    bgc_forc_conf%laddband=laddband
    bgc_forc_conf%Layer=do_layer
    bgc_forc_conf%bgc_fname=bgc_fname
  endif
  do_rgres=do_regression_test
  LYRG=num_of_simdays
  lverb=lverbose
  nmicbguilds=num_microbial_guilds

  JH=JX+1
  JV=JY+1
  JD=JZ+1

end subroutine readnamelist
