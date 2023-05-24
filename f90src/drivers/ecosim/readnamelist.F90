 subroutine readnamelist(nml_buffer, case_name, prefix,LYRG,nmicbguilds)
!!
! Description:
! read control namelist
  use abortutils     , only : endrun
  use EcoSIMConfig   , only : transport_on,column_mode, do_instequil,ref_date
  use EcoSIMConfig   , only : finidat,nrevsn,brnch_retain_casename,start_date
  use ForcWriterMod  , only : bgc_forc_conf,do_bgcforc_write
  use fileUtil       , only : iulog,ecosim_namelist_buffer_size
  use EcoSIMHistMod  , only : DATAC
  use EcoSIMCtrlMod
  use HistFileMod
  use RestartMod     , only : rest_frq,rest_opt
  implicit none
  character(len=*), parameter :: mod_filename = __FILE__
  character(len=ecosim_namelist_buffer_size), intent(in) :: nml_buffer
  character(len=36)    , intent(out) :: case_name
  character(len=80)    , intent(out) :: prefix
  integer              , intent(out) :: LYRG
  integer              , intent(out) :: nmicbguilds

  logical :: do_regression_test
  integer :: num_of_simdays
  logical :: lverbose
  integer :: num_microbial_guilds
  integer :: do_doy,do_year,do_layer
  character(len=64) :: bgc_fname

  namelist /ecosim/case_name, prefix, do_regression_test, &
    num_of_simdays,lverbose,num_microbial_guilds,transport_on,column_mode,&
    do_instequil,salt_model, pft_file_in,grid_file_in,pft_mgmt_in, clm_factor_in,&
    clm_file_in,soil_mgmt_in,forc_periods,&
    NPXS,NPYS,JOUTS,continue_run,visual_out,restart_out,&
    finidat,nrevsn,brnch_retain_casename,plant_model,micb_model,&
    soichem_model,atm_ghg_in

  namelist /ecosim/hist_nhtfrq,hist_mfilt,hist_fincl1,hist_fincl2,hist_yrclose, &
    do_budgets,ref_date,start_date

  logical :: laddband
  namelist /bbgcforc/do_bgcforc_write,do_year,do_doy,laddband,do_layer,&
    bgc_fname

  !local variables
  character(len=256) :: ioerror_msg
  integer :: rc, fu
  integer :: nml_error
  integer :: year0

  continue_run=.false.
  NPXS=30   !number of cycles per hour for water,heat,solute flux calcns
  NPYS=20   !number of cycles per NPX for gas flux calcns
  JOUTS=1   !frequency on hourly scale

  visual_out =.false.
  restart_out=.false.
  do_budgets =.false.
  plant_model=.true.
  soichem_model=.true.
  micb_model=.true.
  ref_date  = '18000101000000'   !place holder for future
  start_date= '18000101000000'   !start date of the simulation, differ from the forcing date
  finidat=' '
  nrevsn = ' '

  brnch_retain_casename=.false.
  hist_yrclose=.false.

  forc_periods=(/1980,1980,1,1981,1988,2,1989,2008,1/)

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

  clm_factor_in=''
  pft_file_in=''
  grid_file_in=''
  pft_mgmt_in=''
  clm_file_in=''
  soil_mgmt_in=''
  atm_ghg_in=''
  inquire (file=nmlfile, iostat=rc)
  if (rc /= 0) then
    write (iulog, '(3a)') 'Error: input file ', trim(nmlfile), &
  ' does not exist.'
    call endrun('stopped in '//trim(mod_filename), __LINE__)
  end if

  open (action='read', file=nmlfile, iostat=rc, newunit=fu)
  if (rc /= 0) then
    write (iulog, '(2a)') 'Error openning input file "', &
    trim(nmlfile)
    call endrun('stopped in '//trim(mod_filename), __LINE__)
  end if

  read(nml_buffer, nml=ecosim, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(iulog,'(a)')"ERROR reading ecosim namelist ",nml_error,ioerror_msg
     call endrun('stopped in '//trim(mod_filename), __LINE__)
  end if

  read(nml_buffer, nml=bbgcforc, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(iulog,'(a)')"ERROR reading bbgcforc namelist "
     call endrun('stopped in '//trim(mod_filename), __LINE__)
  end if

  if (.true.) then
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,ecosim)
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,bbgcforc)
    write(iulog, *)
    write(iulog, *) '--------------------'

  endif
  read(start_date,'(I4)')year0
  call etimer%Init(nml_buffer,year0=year0)
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
  if(.not. soichem_model)then
    salt_model=.false.
  endif
end subroutine readnamelist
