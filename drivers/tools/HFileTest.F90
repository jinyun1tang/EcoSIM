program HFileTest
  use HistFileMod
  use data_kind_mod , only : r8 => DAT_KIND_R8
  use data_const_mod, only : pi=> DAT_CONST_PI
  use GridConsts    , only : bounds,JP
  use fileUtil      , only : iulog
  use abortutils        , only : endrun
  use EcoSIMCtrlMod     , only : etimer
  use EcoSIMConfig
  implicit none
  character(len=*), parameter :: prognm=&
  __FILE__
  real(r8), pointer :: ptr_1d(:)
  real(r8), pointer :: TSOI(:)
  real(r8), pointer :: VSM(:)
  integer :: jh,jd
  character(len=36):: nmlfile

  character(len=256) :: ioerror_msg
  integer :: rc, fu
  integer :: nml_error
  character(len=16) :: ymdhs
  integer :: doy
  logical :: hist_yrclose
  logical :: rstwr, nlend, lnyr
  namelist /hist_nml/&
    hist_nhtfrq,hist_mfilt,hist_fincl1,hist_fincl2,hist_yrclose


  CALL GETARG(1,nmlfile)

  hist_yrclose=.false.
  inquire (file=nmlfile, iostat=rc)
  if (rc /= 0) then
    write (iulog, '(3a)') 'Error: input file ', trim(nmlfile), ' does not exist.'
    call endrun('stopped in '//trim(prognm), __LINE__)
  end if

  open (action='read', file=nmlfile, iostat=rc, newunit=fu)
  if (rc /= 0) then
    write (iulog, '(2a)') 'Error openning input file "', trim(nmlfile)
    call endrun('stopped in '//trim(prognm), __LINE__)
  end if

  read(unit=fu, nml=hist_nml, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     print*,ioerror_msg
     write(iulog,'(a)')"ERROR reading namelist file ",trim(nmlfile)
     call endrun('stopped in '//trim(prognm), __LINE__)
  end if

  close(fu)
  if (.true.) then
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,hist_nml)

  endif

  bounds%ngrid=1

  bounds%ntopou=1

  bounds%NHW =1
  bounds%NVN =1
  bounds%NHE =1
  bounds%NVS =1

  bounds%begg=1;bounds%endg=bounds%ngrid
  bounds%begt=1;bounds%endt=bounds%ntopou
  bounds%ncols=2
  bounds%npfts=bounds%ncols*JP
  bounds%begc=1;bounds%endc=bounds%ncols
  bounds%begp=1;bounds%endp=bounds%npfts

  allocate(TSOI(2))
  allocate(VSM(2))
  ptr_1d=>TSOI
  call hist_addfld1d(fname='TSOI',units='oC',avgflag='A',long_name='soil temperature',ptr_col=ptr_1d)
  ptr_1d=>TSOI
  call hist_addfld1d(fname='TMXSOI',units='oC',avgflag='X',long_name='maximum soil temperature',ptr_col=ptr_1d)
  ptr_1d=>TSOI
  call hist_addfld1d(fname='TMNSOI',units='oC',avgflag='M',long_name='minimum soil temperature',ptr_col=ptr_1d)

  ptr_1d=>VSM
  call hist_addfld1d(fname='VSM', units='m3/m3',avgflag='A',long_name='soil moisture',ptr_col=ptr_1d)
  ptr_1d=>VSM
  call hist_addfld1d(fname='VMXSM', units='m3/m3',avgflag='X',long_name='maximum soil moisture',ptr_col=ptr_1d)
  ptr_1d=>VSM
  call hist_addfld1d(fname='VMNSM', units='m3/m3',avgflag='M',long_name='minimum soil moisture',ptr_col=ptr_1d)

  case_name='htest'

  call hist_htapes_build()

  call etimer%Init(year0=1804,nyears=12)

  call etimer%setClock(dtime=3600._r8,nelapstep=0)

  DO while(.true.)
    do while(.true.)
      call etimer%get_ymdhs(ymdhs)
      doy=etimer%get_curr_doy()
      DO Jh=1,24
        TSOI(1)=1._r8+cos(2._r8*pi*(jh-12._r8)/24._r8)+sin(2._r8*pi*(doy-182.5_r8)/365._r8)*3._r8
        TSOI(2)=2._r8+cos(2._r8*pi*(jh-12._r8)/24._r8)*0.2_r8+sin(2._r8*pi*(doy-182.5_r8)/365._r8)*3._r8
        VSM(1) =0.2_r8+sin(2._r8*pi*jh/24._r8)*0.1_r8+sin(2._r8*pi*(doy-182.5_r8)/365._r8)*0.1_r8
        VSM(2) =0.3_r8+sin(2._r8*pi*jh/24._r8)*0.2_r8+sin(2._r8*pi*(doy-182.5_r8)/365._r8)*0.1_r8
        call hist_update_hbuf(bounds)
        call etimer%update_time_stamp()
        nlend=etimer%its_time_to_exit()
        rstwr=.false.
        lnyr=etimer%its_a_new_year().and.hist_yrclose
        call hist_htapes_wrapup( rstwr, nlend, bounds, lnyr )
      ENDDO

      if(etimer%its_a_new_month())print*,ymdhs,etimer%get_curr_day()
      if(etimer%its_a_new_year())exit
    enddo
    if(nlend)exit
  enddo

end program HFileTest
