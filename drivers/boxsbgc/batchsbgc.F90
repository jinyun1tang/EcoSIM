program main
!
! Single layer model
  use batchmod
  use abortutils, only : endrun
  use fileUtil
  use data_kind_mod, only : r8 => DAT_KIND_R8
implicit none
  character(len=*), parameter :: mod_filename = &
  __FILE__

  integer :: arg_count
  character(len=ecosim_filename_length)      :: namelist_filename
  character(len=ecosim_namelist_buffer_size) :: namelist_buffer
  character(len=ecosim_filename_length)      :: base_filename

  arg_count = command_argument_count()
  if (arg_count /= 1) then
     write(stdout, '(a, i3)') 'ERROR: must pass exactly one arguement one the command line, received ', arg_count
     call usage()
     call endrun()
  end if

  call get_command_argument(1, namelist_filename)

  write(stdout, '(a, a)') 'Reading namelist filename : ', trim(namelist_filename)
  namelist_buffer = ''

  call namelist_to_buffer(namelist_filename, namelist_buffer)

  call RunModel(namelist_buffer)

end program main


! ----------------------------------------------------------------------
  subroutine usage()
  !DESCRIPTION
  !display something
  use fileUtil, only : stdout
  implicit none
  write(stdout, *) 'boxsbgc.x- standalone driver for ecosim 1-layer soilbgc library.'
  write(stdout, *) 'usage: boxsbgc.x namelist_filename'
  end subroutine usage
! ----------------------------------------------------------------------

  subroutine RunModel(namelist_buffer)
  !DESCRIPTION
  !display something
  use data_kind_mod       , only : r8 => DAT_KIND_R8
  use MicFLuxTypeMod      , only : micfluxtype
  use ModelStatusType     , only : model_status_type
  use MicStateTraitTypeMod, only : micsttype
  use MicForcTypeMod      , only : micforctype
  use ecosim_Time_Mod     , only : ecosim_time_type
  use ecosim_log_mod      , only : errMsg => shr_log_errMsg
  use batchmod
  use ForcTypeMod         , only : forc_type,ReadForc,UpdateForc
  use abortutils          , only : endrun

  use bhistMod
  use fileUtil
  use EcoSIMSolverPar
  implicit none
  character(len=*), intent(in) :: namelist_buffer
!
! local argument
  type(forc_type)   :: forc
  type(micforctype) :: micfor
  type(micsttype) :: micstt
  type(micfluxtype) :: micflx

  character(len=hist_var_str_len) , allocatable :: varl(:)
  character(len=hist_var_lon_str_len) , allocatable :: varlnml(:)
  character(len=hist_unit_str_len), allocatable :: unitl(:)
  character(len=hist_freq_str_len), allocatable :: freql(:)
  integer                         , allocatable :: vartypes(:)
  real(r8), allocatable :: ystatesf(:,:)
  real(r8), allocatable :: ystates0l(:)
  real(r8), allocatable :: ystatesfl(:)
  integer :: nvars,ncols,jj
  real(r8) :: dtime
  type(ecosim_time_type) :: timer
  type(model_status_type) :: err_status
  type(histf_type) :: hist
  character(len=256) :: ioerror_msg
  character(len=256):: gname
  character(len=14) :: yymmddhhss
  integer :: rc, fu
  integer :: nml_error
  character(len=hist_freq_str_len) :: hist_freq
  character(len=32) :: model_name
  character(len=32) :: case_id
  character(len=32) :: prefix
  character(len=128):: forc_file
  logical :: salton
  real(r8) :: CO2E                       !initial atmospheric CO2 concentration, [umol mol-1],ppmv
  real(r8) :: OXYE                       !atmospheric O2 concentration, [umol mol-1],ppmv
  real(r8) :: Z2OE                       !atmospheric N2O concentration, [umol mol-1],ppmv
  real(r8) :: Z2GE                       !atmospheric N2 concentration, [umol mol-1],ppmv
  real(r8) :: ZNH3E                      !atmospheric NH3 concentration, [umol mol-1],ppmv
  real(r8) :: CH4E                       !atmospheric CH4 concentration, [umol mol-1],ppmv
  real(r8) :: H2GE                       !atmospheric H2 concentration, [umol mol-1],ppmv
  integer :: forctype                    ! 0: (transient), 1: T const, 2: water const, 3: T and water const
  logical :: disvolonly                  !dissolution/volatilization only
  namelist /driver_nml/model_name,case_id,hist_freq,salton,forc_file,&
    CO2E,OXYE,Z2OE,Z2GE,ZNH3E,CH4E,H2GE,forctype,disvolonly

  character(len=*), parameter :: mod_filename=&
  __FILE__

  hist_freq='day'
  model_name='boxsbgc'
  case_id='exp0'
  salton=.false.
  forc_file='bbforc.nc'

  disvolonly=.false.
  OXYE=2.1E+05_r8
  Z2GE=7.8E+05_r8
  CO2E=370.0_r8
  CH4E=1.8_r8
  Z2OE= 0.3_r8
  ZNH3E=0.005_r8
  H2GE=1.e-3_r8
  forctype=0

  if ( .true. )then
     ioerror_msg=''
     read(namelist_buffer, nml=driver_nml, iostat=nml_error, iomsg=ioerror_msg)
     if (nml_error /= 0) then
        write(*,*)ioerror_msg
        call endrun(msg="ERROR reading driver_nml namelist "//errmsg(mod_filename, __LINE__))
     end if
  end if
  if (.true.) then
    write(stdout, *)
    write(stdout, *) '--------------------'
    write(stdout,driver_nml)
    write(stdout, *)
    write(stdout, *) '--------------------'
  endif

! set up solver

  NPH=1;NPT=15;NPG=NPT*NPH;dts_gas=1.0_r8/NPG;dts_HeatWatTP=1._r8/NPH;dt_GasCyc=1.0_r8/NPT

  ncols=1
!!============================================================
! customized model varlist creation
  nvars=getvarllen()
  allocate(varl(nvars)); allocate(varlnml(nvars))
  allocate(unitl(nvars)); allocate(freql(nvars)); allocate(vartypes(nvars))
  call getvarlist(nvars,varl, varlnml, unitl, vartypes)
!============================================================

  freql(:) = hist_freq
  allocate(ystatesf(1:ncols,1:nvars));ystatesf(:,:)=0._r8
  allocate(ystates0l(1:nvars));       ystates0l(:) = 0._r8
  allocate(ystatesfl(1:nvars));       ystatesfl(:) = 0._r8


  forc%OXYE =OXYE
  forc%Z2GE =Z2GE
  forc%CO2E =CO2E
  forc%CH4E =CH4E
  forc%Z2OE =Z2OE
  forc%ZNH3E=ZNH3E
  forc%H2GE =H2GE
  forc%disvolonly=disvolonly

  call ReadFORC(forc,forc_file)

  call micfor%Init()
  call micstt%Init()
  call micflx%Init()

! customized model initialization
  call initmodel(nvars, ystates0l, forc, err_status)


  if(err_status%check_status())then
    call endrun(msg=err_status%print_msg())
  endif

  !initialize timer
  call timer%Init(namelist_buffer=namelist_buffer)

  dtime=timer%get_step_size()

  if(len(trim(case_id))==0)then
    write(gname,'(A)')'boxsbgc'//'.'//trim(model_name)
  else
    write(gname,'(A)')'boxsbgc'//'.'//trim(case_id)//'.'//trim(model_name)
  endif

  call hist%init(ncols,varl, varlnml, unitl, vartypes, freql, gname, dtime)


  DO

    call timer%update_time_stamp()

!   setup forcing, e.g., add litter/om/nutrients, set up temperature/moisture
!
    call UpdateFORC(forc,forctype)

    call BatchModelConfig(nvars,ystates0l,forc,micfor,micstt,micflx,err_status)

    if(err_status%check_status())then
      call endrun(msg=err_status%print_msg())
    endif

    call RunMicBGC(nvars, ystates0l, ystatesfl, forc,micfor,micstt,micflx, err_status)

    call timer%update_time_stamp()

    do jj = 1, nvars
      ystatesf(1,jj)=ystatesfl(jj)
      ystates0l(jj) =ystatesfl(jj)
      forc%ORGC=micfor%ORGC
    enddo

    call hist%hist_wrap(ystatesf, timer)

    if(timer%its_a_new_year())then
      write(iulog,*)'year ',timer%get_curr_year()
    endif
    if(timer%its_time_to_exit())exit

  enddo
  call timer%get_ymdhs(yymmddhhss)

  call hist%histrst('boxsbgc.x', 'write', yymmddhhss)

  call micfor%destroy()
  call micflx%destroy()
  call micstt%destroy()
  end subroutine runModel
