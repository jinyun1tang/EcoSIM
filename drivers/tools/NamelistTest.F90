PROGRAM NamelistTest

  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod     , only : etimer
  use ecosim_time_mod   , only : getdow
  use abortutils        , only : endrun
  implicit none
  character(len=*), parameter :: mod_filename = &
  __FILE__

  logical :: continue_run
  character(len=16) :: hist_config(10)
  character(len=8)  :: sim_yyyymmdd
  integer :: forc_year0
  integer :: forc_periods(5)
  integer :: NPXS(5),NPYS(5),J
  namelist /ecosys_ctrl/hist_config,sim_yyyymmdd,forc_year0,forc_periods,&
    NPXS,NPYS,continue_run

  !local variables
  character(len=256) :: ioerror_msg
  character(len=14) :: ymdhs
  character(len=14) :: ymdhs0
  integer :: rc, fu
  integer :: nml_error
  character(len=*), parameter :: nmlfile='example_nl'
  integer :: iulog
  integer :: nn,nyr,year0,year,year1
  logical :: lskip_loop

  iulog=6
  lskip_loop=.true.

  continue_run=.false.
  NPXS=30   !number of cycles per hour for water,heat,solute flux calcns
  NPYS=20   !number of cycles per NPX for gas flux calcns

  hist_config='NO'
  sim_yyyymmdd='18000101'
  forc_year0=1980
  forc_periods=(/1981,1988,2,1989,2008/)

  inquire (file=nmlfile, iostat=rc)
  if (rc /= 0) then
    write (iulog, '(3a)') 'Error: input file ', trim(nmlfile), &
  ' does not exist.'
    call endrun('stopped in readnml ', __LINE__)
  end if

  open (action='read', file=nmlfile, iostat=rc, newunit=fu)
  if (rc /= 0) then
    write (iulog, '(2a)') 'Error openning input file "', &
  trim(nmlfile)
    call endrun('stopped in readnml ', __LINE__)
  end if

  read(unit=fu, nml=ecosys_ctrl, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     print*,nml_error,ioerror_msg
     write(iulog,'(a)')"ERROR reading ecosys_ctrl namelist "
     call endrun('stopped in '//mod_filename, __LINE__)
  end if

  close(fu)
  if (.true.) then
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,ecosys_ctrl)
    write(iulog, *)
    write(iulog, *) '--------------------'
  endif


  read(sim_yyyymmdd,'(I4)')year

  year1=0
  if(continue_run)then
    print*,'read restart file'
    ymdhs0='18820114000000'
    read(ymdhs0,'(I4)')year1
    call etimer%Init(year0=year1)
  else
    ymdhs0='00000000000000'
    ymdhs0(1:8)=sim_yyyymmdd
    call etimer%Init(year0=year)
  endif
  call etimer%setClock(dtime=3600._r8,nelapstep=0)

  year0=forc_year0
  read(sim_yyyymmdd,'(I4)')year
  print*,'year',year1,year
  call etimer%get_ymdhs(ymdhs)
  print*,ymdhs,' ',ymdhs0
  do while(year1<=year)
    call etimer%get_ymdhs(ymdhs)
    if(ymdhs==ymdhs0)lskip_loop=.false.
    call etimer%update_time_stamp()
    if (etimer%its_a_new_year())exit
  enddo

  call etimer%get_ymdhs(ymdhs)
  if(ymdhs==ymdhs0)lskip_loop=.false.
  if(.not.lskip_loop)print*,year,nyr,ymdhs

  DO nn=1,forc_periods(3)
    do nyr=forc_periods(1),forc_periods(2)
      year=year+1

      do while(year1<=year)
        call etimer%get_ymdhs(ymdhs)
        if(ymdhs==ymdhs0)lskip_loop=.false.
        call etimer%update_time_stamp()
        if (etimer%its_a_new_year())exit
      enddo
      call etimer%get_ymdhs(ymdhs)
      if(ymdhs==ymdhs0)lskip_loop=.false.
      if(.not.lskip_loop)print*,year,nyr,ymdhs
    enddo
  enddo

  do nyr=forc_periods(4),forc_periods(5)
    year=year+1
    do while(year1<=year)
      call etimer%get_ymdhs(ymdhs)
      if(ymdhs==ymdhs0)lskip_loop=.false.
      call etimer%update_time_stamp()
      if (etimer%its_a_new_year())exit
    enddo
    call etimer%get_ymdhs(ymdhs)
    if(ymdhs==ymdhs0)lskip_loop=.false.
    if(.not.lskip_loop)print*,year,nyr,ymdhs
  enddo
  print*,lskip_loop
end program NamelistTest
