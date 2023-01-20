 subroutine readnamelist(nmlfile,runfile, case_name, &
  prefix,do_rgres,LYRG,lverb,nmicbguilds)
!!
! Description:
! read control namelist
  use abortutils   , only : endrun
  use EcoSIMConfig , only : transport_on,column_mode, do_instequil
  use ForcWriterMod, only : bgc_forc_conf,do_bgcforc_write
  use fileUtil     , only : iulog
  use EcoSIMCtrlMod, only : salt_model, pft_file_in,grid_file_in

  implicit none
  character(len=*), parameter :: mod_filename = __FILE__
  character(len=*), intent(in) :: nmlfile
  character(len=80), intent(out) :: runfile
  character(len=36)    , intent(out) :: case_name
  character(len=80)    , intent(out) :: prefix
  logical              , intent(out) :: do_rgres
  integer              , intent(out) :: LYRG
  logical              , intent(out) :: lverb
  integer              , intent(out) :: nmicbguilds

  logical :: do_regression_test
  integer :: num_of_simdays
  logical :: lverbose
  integer :: num_microbial_guilds
  integer :: do_doy,do_year,do_layer
  character(len=64) :: bgc_fname

  namelist /ecosys/case_name, prefix, runfile, do_regression_test, &
  num_of_simdays,lverbose,num_microbial_guilds,transport_on,column_mode,&
  do_instequil,salt_model, pft_file_in,grid_file_in


  logical :: laddband
  namelist /bbgcforc/do_bgcforc_write,do_year,do_doy,laddband,do_layer,&
    bgc_fname

  !local variables
  character(len=256) :: ioerror_msg
  integer :: rc, fu
  integer :: nml_error

  num_of_simdays=-1
  do_year=-1
  do_doy=0
  do_layer=1
  salt_model=.false.
  laddband=.false.
  do_regression_test=.false.
  lverbose=.false.
  num_microbial_guilds=1
  do_bgcforc_write=.false.
  bgc_fname='bbforc.nc'
  do_instequil=.false.
  pft_file_in=''
  grid_file_in=''
  inquire (file=nmlfile, iostat=rc)
  if (rc /= 0) then
    write (iulog, '(3a)') 'Error: input file ', trim(nmlfile), &
  ' does not exist.'
    call endrun('stopped in readnml', 25)
  end if

  open (action='read', file=nmlfile, iostat=rc, newunit=fu)
  if (rc /= 0) then
    write (iulog, '(2a)') 'Error openning input file "', &
  trim(nmlfile)
    call endrun('stopped in readnml', 32)
  end if

  read(unit=fu, nml=ecosys, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(iulog,'(a)')"ERROR reading ecosys namelist "
     call endrun('stopped in readnml', 38)
  end if

  read(unit=fu, nml=bbgcforc, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(iulog,'(a)')"ERROR reading bbgcforc namelist "
     call endrun('stopped in readnml', 38)
  end if

  close(fu)
  if (.true.) then
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,ecosys)
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,bbgcforc)
    write(iulog, *)
    write(iulog, *) '--------------------'

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

  !below is a temporary setup

end subroutine readnamelist
