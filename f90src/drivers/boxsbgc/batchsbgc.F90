program main
!
! Single layer model
  use batchmod
  use abortutils, only : endrun
  use fileUtil, only : ecosim_filename_length
implicit none
  character(len=*), parameter :: mod_filename = __FILE__

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

  call readNML(namelist_filename)

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

  subroutine RunModel()
  !DESCRIPTION
  !display something
  use MicFLuxTypeMod, only : micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use MicForcTypeMod, only : micforctype
  use ecosim_Time_Mod, only : ecosim_time_type
  use histMod
  implicit none

  type(micforctype) :: micfor
  type(micsttype) :: micstt
  type(micfluxtype) :: micflx
  type(ecosim_time_type) :: timer
  real(r8) :: dtime
  integer  :: jj,nvars

  character(len=hist_var_str_len) , allocatable :: varl(:)
  character(len=hist_unit_str_len), allocatable :: unitl(:)
  character(len=hist_freq_str_len), allocatable :: freql(:)
  integer                         , allocatable :: vartypes(:)
  real(r8), allocatable :: ystatesf(:,:)
  real(r8), allocatable :: ystates0l(:)
  real(r8), allocatable :: ystatesfl(:)

  character(len=64) :: case_name
  logical :: is_surflit = .false.  !logical switch for litter decomposition
  character(len=hist_freq_str_len) :: hist_freq

  namelist / model_driver / case_name,is_surflit,hist_freq

  case_name='boxbgc'
  is_surflit=.false.
  hist_freq='day'
  if ( .true. )then
     ioerror_msg=''
     read(namelist_buffer, nml=model_driver, iostat=nml_error, iomsg=ioerror_msg)
     if (nml_error /= 0) then
        call endrun(msg="ERROR reading forcing_inparm namelist "//errmsg(mod_filename, __LINE__))
     end if
  end if

  call micfor%Init()
  call micstt%Init()
  call micflx%Init()

  call Initboxbgc(nvars)

  allocate(varl(nvars)); allocate(unitl(nvars)); allocate(freql(nvars)); allocate(vartypes(nvars))

  call getvarlist(nvars, varl, unitl, vartypes)

  freql(:) = hist_freq
  allocate(ystatesf(1,1:nvars));ystatesf(1,:)=0._r8
  allocate(ystates0l(1:nvars));ystates0l(:) = 0._r8
  allocate(ystatesfl(1:nvars));ystatesfl(:) = 0._r8


  !initialize timer
  call timer%Init(namelist_buffer=namelist_buffer)

  dtime=timer%get_step_size()

  if(len(trim(case_id))==0)then
    write(gname,'(A)')'jarmodel'//'.'//trim(jarmodel_name)
  else
    write(gname,'(A)')'jarmodel'//'.'//trim(case_id)//'.'//trim(jarmodel_name)
  endif

  call hist%init(1,varl,unitl, vartypes, freql, gname, dtime)

  call ReadForc(forc)

  call Initboxbgc()

  DO
    call BatchModelConfig(forc,micfor,micstt,micflx, ystatesf0l)

!   computes the fluxes
    call SoilBGCOneLayer(I,J,micfor,micstt,micflx)

    call UpdateStateVars(micstt,micflx,ystatesfl)
!
    call timer%update_time_stamp()

    do jj = 1, nvars
      ystatesf(1,jj)=ystatesfl(jj)
    enddo

    call hist%hist_wrap(ystatesf, timer)
    if(timer%its_a_new_year())then
      write(iulog,*)'year ',timer%get_cur_year()
    enddo
    if(timer%its_time_to_exit())exit
  ENDDO
  call micfor%destroy()
  call micflx%destroy()
  call micstt%destroy()
  end subroutine runModel
