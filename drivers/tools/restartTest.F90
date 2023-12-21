program restartTest
  use HistFileMod
  use RestartMod
  use abortutils        , only : endrun
  use fileUtil          , only : iulog
  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod     , only : etimer
  use GridConsts        , only : JZ,JS,MaxNumBranches,bounds,JP
implicit none
  character(len=*), parameter :: prognm=&
  __FILE__

  character(len=36):: nmlfile

  character(len=256) :: ioerror_msg
  integer :: rc, fu
  integer :: nml_error
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
  print*,prognm

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

  call etimer%Init(year0=1804,nyears=12)

  call etimer%setClock(dtime=3600._r8,nelapstep=0)

  call restFile(flag='write')
end program restartTest
