program PlantManagementReader
  !
  !! Description
  ! test the plant management reading code
  use PlantMgmtDataType
  use PlantInfoMod
  use EcoSIMCtrlMod, only: pft_mgmt_in, lverb
  use abortutils,    only: endrun
  USE fileUtil,      ONLY: iulog, ecosim_namelist_buffer_size, namelist_to_buffer

  implicit none
  character(len=*), parameter :: mod_filename = &
  __FILE__

  character(len=ecosim_namelist_buffer_size) :: nml_buffer
  character(len=40) :: progname
  character(len=256):: nmlfile
  integer :: num_args
  integer :: yearc  !current model year
  integer :: yeari  !current forcing year 
  integer :: NHW,NHE,NVN,NVS
  integer :: nyears   !number of simulation years 
  integer :: rc, fu
  integer :: nml_error
  integer :: year0,yearf,nyr,ndyr
  integer :: NPS
  character(len=256) :: ioerror_msg

  namelist /pft_mgmnt/year0,yearf,NHW,NHE,NVN,NVS,pft_mgmt_in,NPS,lverb,nyears

  CALL GETARG(0,progname)

  num_args = iargc()

  if (num_args<1)then
    call endrun(msg='Please use code as: '//trim(progname)//' infile')
  endif
  call GETARG(1,nmlfile)

  call namelist_to_buffer(nmlfile,nml_buffer)

  open (action='read', file=nmlfile, iostat=rc, newunit=fu)
  if (rc /= 0) then
    write (iulog, '(2a)') 'Error openning input file "', &
  trim(nmlfile)
    call endrun('stopped in readnml ', __LINE__)
  end if
  NHW=1
  NHE=1
  NVN=1
  NVS=1
  call InitPlantMngmtData

  NP(NVN:NVS,NHW:NHE)=NPS
  read(nml_buffer, nml=pft_mgmnt, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(iulog,'(a)')"ERROR reading pft management namelist ",nml_error,ioerror_msg
     call endrun('stopped in '//trim(mod_filename), __LINE__)
  end if

  if (.true.) then
    write(iulog, *)
    write(iulog, *) '--------------------'
    write(iulog,pft_mgmnt)
    write(iulog, *)
    write(iulog, *) '--------------------'
  endif

  do nyr=1,nyears
    yearc=year0+nyr-1
    ndyr=mod(nyr,yearf-year0+1)
    if(ndyr==0)then
      yeari=yearf
    else
      yeari=year0+ndyr
    endif
    write(*,*)'yeari',yeari
    call ReadPlantManagementNC(yearc,yeari,NHW,NHE,NVN,NVS)
  enddo
  call DestructPlantMngmtData
end program PlantManagementReader
