module ecosim_Time_Mod
!
! DESCRIPTION
! the module contains subroutine to march the time

  use data_kind_mod  , only : r8 => shr_kind_r8
  use data_kind_mod  , only : i8 => shr_kind_i8
  use fileUtil, only : stdout, ecosim_string_length_long

  implicit none
  private
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8), parameter :: secpday=86400._r8
  real(r8), parameter :: secpyear=86400._r8*365._r8   !seconds in a normal year
  integer , parameter :: daz(12)=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)  
  integer , parameter :: cdaz(12)=(/31,59,90,120,151,181,212,243,273,304,334,365/)
  type, public:: ecosim_time_type
     ! NOTE(bja, 201603) all real variables have units of seconds!
     real(r8) :: delta_time
     real(r8) :: stop_time
     real(r8) :: time
     real(r8) :: time0
     real(r8) :: timef
     real(r8) :: toy                 !time of year
     real(r8) :: restart_dtime
     integer  :: tstep               !number of time steps
     integer  :: nelapstep
     integer  :: dow, dom, doy       !day of week, day month, day of year
     integer  :: moy                 !mon of year
     integer  :: cyears              !cumulative years
     integer  :: cdays               !cumulative days
     integer  :: year0   !beginning year AD
     real(r8) :: tod     !time of day
     integer  :: hist_freq   !negative number, steps, positive number, 1: day, 30:mon, 365:year
     integer  :: stop_opt
     integer  :: leap_yr     
   contains
     procedure, public :: Init
     procedure, public :: its_time_to_write_restart
     procedure, public :: its_time_to_exit
     procedure, public :: update_time_stamp
     procedure, public :: set_nstep
     procedure, public :: get_nstep
     procedure, public :: set_time_offset
     procedure, public :: get_days_per_year
     procedure, public :: get_step_size
     procedure, public :: get_cur_time
     procedure, public :: get_cur_timef

     procedure, public :: proc_initstep
     procedure, public :: print_cur_time
     procedure, public :: its_time_to_histflush

     procedure, public :: setClock
     procedure, public :: its_a_new_hour
     procedure, public :: its_a_new_day
     procedure, public :: its_a_new_week
     procedure, public :: its_a_new_month
     procedure, public :: its_a_new_year
     procedure, public :: get_ymdhs
     procedure, public :: get_cur_year
     procedure, public :: get_cur_yearAD
     procedure, public :: get_cur_day
     procedure, public :: get_cur_mon     
     procedure, public :: is_first_step
     procedure, public :: print_model_time_stamp
     procedure, private:: proc_nextstep
     procedure, private:: ReadNamelist
  end type ecosim_time_type
  public :: getdow

contains

  subroutine setClock(this, dtime, nelapstep)

  implicit none
  class(ecosim_time_type), intent(inout) :: this
  real(r8), intent(in) :: dtime
  integer, intent(in) :: nelapstep

  this%delta_time = dtime
  this%nelapstep = nelapstep

  if(this%year0>0)then
    if(isLeap(this%year0))then
      this%leap_yr=1
    else
      this%leap_yr=0
    endif  
  endif
  end subroutine setClock
  !-------------------------------------------------------------------------------
  function get_cur_time(this)result(ans)

  implicit none
  class(ecosim_time_type), intent(in) :: this
  real(r8) :: ans

  ans = this%time
  end function get_cur_time
  !-------------------------------------------------------------------------------
  function get_cur_timef(this)result(ans)

  implicit none
  class(ecosim_time_type), intent(in) :: this
  real(r8) :: ans

  ans = this%time+this%time0
  end function get_cur_timef
  !-------------------------------------------------------------------------------
  subroutine Init(this, namelist_buffer, masterproc, year0, nyears)

    implicit none

    class(ecosim_time_type), intent(inout) :: this    
    character(len=*), optional, intent(in) :: namelist_buffer
    logical, optional, intent(in) :: masterproc
    integer, optional, intent(in) :: year0    !beginning year, used to record in the AD format
    integer, optional, intent(in) :: nyears   !run the timer for nyears
    integer :: N

    this%tstep = 0
    this%time0 = 0._r8
    this%time  = 0._r8
    this%tod   = 0._r8
    this%toy   = 0._r8
    this%cyears = 0
    this%year0 = 0
    if(present(year0))this%year0=year0
    this%cdays  = 0
    this%dow    = 1
    this%dom    = 1
    this%doy    = 1
    this%moy    = 1
    this%hist_freq=-1
    this%nelapstep=0
    this%stop_opt=3
    this%leap_yr=0
    
    if(present(namelist_buffer))then
      if(present(masterproc))then
        call this%ReadNamelist(namelist_buffer, masterproc)
      else
        call this%ReadNamelist(namelist_buffer)
      endif
    else 
      if(present(nyears))then
        this%stop_time= nyears * 365._r8
        if(this%year0 > 0)then
          do n=0,nyears-1
            if(isLeap(this%year0+n))this%stop_time=this%stop_time+1._r8
          enddo
          if(isLeap(this%year0))this%leap_yr=1
        endif
        this%stop_time=this%stop_time * 86400._r8
      endif   
    endif
  end subroutine Init

  ! ----------------------------------------------------------------------

  subroutine ReadNamelist(this, namelist_buffer, masterproc)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use abortutils    , only : endrun
    use ecosim_log_mod   , only : errMsg => shr_log_errMsg
    use fileUtil , only : stdout, ecosim_string_length_long

    implicit none
    ! !ARGUMENTS:
    class(ecosim_time_type), intent(inout) :: this
    character(len=*), intent(in) :: namelist_buffer
    logical, optional, intent(in) :: masterproc
    !
    ! !LOCAL VARIABLES:
    integer :: nml_error
    character(len=*), parameter :: subname = 'ecosim_time%ReadNamelist'
    character(len=ecosim_string_length_long) :: ioerror_msg
    real(r8) :: delta_time         !model time step
    real(r8) :: stop_time          !when to stop
    integer  :: stop_n
    integer  :: year0
    integer  :: hist_freq
    character(len=8)  :: stop_option        !1: step, 2:day, 3:year
    real(r8) :: restart_dtime      !when to write restart file
    logical :: masterproc_loc
    integer :: jj
    !-----------------------------------------------------------------------

    namelist / ecosim_time / delta_time, stop_n, stop_option, restart_dtime, hist_freq, year0

    ! FIXME(bja, 201603) Only reading time variables in seconds!
    ! Should enable other values with unit coversions.

    ! FIXME(bja, 201603) assign some defaults, should eventually remove
    ! when all input files are updated.
    masterproc_loc=.true.
    if(present(masterproc))masterproc_loc=masterproc
    delta_time = 1800._r8                !half hourly time step
    stop_n=2       !by default 2 cycle
    stop_option='nyears'  !by default years
    this%stop_opt=3
    hist_freq=-1   !write every time step
    restart_dtime = -1._r8
    year0=0
    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( index(trim(namelist_buffer),'ecosim_time')/=0 )then
       ioerror_msg=''
       read(namelist_buffer, nml=ecosim_time, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          call endrun(msg="ERROR reading secosim_driver namelist "//errmsg(mod_filename, __LINE__))
       end if
    end if

    if (masterproc_loc) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' ecosim time :'
       write(stdout, *)
       write(stdout, *) ' ecosim_time namelist settings :'
       write(stdout, *)
       write(stdout, ecosim_time)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif
    this%hist_freq = hist_freq
    this%delta_time = delta_time
    this%stop_time = delta_time*stop_n
    if(year0>0)this%year0=year0
    select case (trim(stop_option))
    case ('nsteps')
      this%stop_time = stop_n * this%delta_time
      this%stop_opt=1
    case ('ndays')
      !day
      this%stop_time= stop_n * 86400._r8
      this%stop_opt=2
    case ('nyears')
      !year
      this%stop_time= stop_n * 365._r8
      this%stop_opt=3
      if(this%year0>0)then
        do jj = 1, stop_n
          if(isLeap(jj+this%year0-1))this%stop_time=this%stop_time+1._r8          
        enddo
      endif
      this%stop_time=this%stop_time*86400._r8
    case default
      call endrun(msg="ERROR setting up stop_option "//errmsg(mod_filename, __LINE__))
    end select
    if(restart_dtime < 0._r8)then
      this%restart_dtime  = this%stop_time
    else
      this%restart_dtime = restart_dtime
    endif
  end subroutine ReadNamelist

  !-------------------------------------------------------------------------------
    function its_time_to_write_restart(this)result(ans)
     !
     ! DESCRIPTION
     ! decide if to write restart file
     !

     implicit none
     class(ecosim_time_type), intent(in) :: this
     logical :: ans

     character(len=80) :: subname = 'its_time_to_write_restart'

     ans = (mod(this%time,this%restart_dtime) == 0)
     end function its_time_to_write_restart

  !-------------------------------------------------------------------------------
  function its_time_to_exit(this) result(ans)
    !
    ! DESCRIPTION
    ! decide if to exit the loop
    !

    implicit none
    class(ecosim_time_type), intent(in) :: this
    logical :: ans

    character(len=80) :: subname = 'its_time_to_exit'

    ans = (this%time >= this%stop_time)

  end function its_time_to_exit

  !-------------------------------------------------------------------------------
  subroutine update_time_stamp(this)
    !
    ! DESCRIPTION
    !

    implicit none

    class(ecosim_time_type), intent(inout) :: this

    character(len=80) :: subname='update_time_stamp'

    this%time = this%time + this%delta_time
    this%toy  = this%toy + this%delta_time

    this%tstep = this%tstep + 1
    !
    ! reset the clock every year, and assuming the time step
    ! size is always
    if(mod(this%toy, secpyear+86400._r8*this%leap_yr) == 0) then
       this%tstep = 0
    end if

    !update time of the day
    this%tod=this%tod+this%delta_time

    !update varaibles when it is to start a new day
    if(this%its_a_new_day())then
      !it reaches 24:00, reset to 00:00
      this%tod=0._r8
      this%dom=this%dom+1
      this%dow = mod(this%dow + 1, 7)
      if(this%dow==0)this%dow=7
      this%doy = this%doy + 1
      this%cdays= this%cdays + 1
    endif

    if(this%its_a_new_month())then
      this%moy=this%moy+1
      if(this%moy==13)then
        this%moy=1
        this%doy=mod(this%doy,365+this%leap_yr)
        if(this%doy==0)this%doy=1
      endif
      this%dom=1
    endif
    if(this%its_a_new_year())then
      this%cyears=this%cyears+1
      if(this%year0>0 .and. isLeap(this%cyears+this%year0))then
        this%leap_yr=1
      else
        this%leap_yr=0
      endif        
    endif  
    call this%proc_nextstep()
  end subroutine update_time_stamp

  !-------------------------------------------------------------------------------
  real(r8) function get_step_size(this) result(rc)

    ! Return the step size in seconds.

    implicit none

    class(ecosim_time_type), intent(in) :: this

    rc = this%delta_time
    return

  end function get_step_size

  !-------------------------------------------------------------------------------
  integer function get_nstep(this)

    ! Return the timestep number.
    implicit none

    class(ecosim_time_type), intent(in) :: this

    character(len=*), parameter :: sub = 'betr::get_nstep'

    get_nstep = this%nelapstep

  end function get_nstep


  !-------------------------------------------------------------------------------
  subroutine set_nstep(this, nstep)

    ! Return the timestep number.
    implicit none
    class(ecosim_time_type), intent(inout) :: this

    character(len=*), parameter :: sub = 'betr::get_nstep'
    integer, intent(in) :: nstep

    this%nelapstep = nstep

    if(this%its_a_new_year())then
      this%tstep = 0
    endif

  end subroutine set_nstep


  !-------------------------------------------------------------------------------
  subroutine set_time_offset(this, nstep, tstep, continue_run)

    ! Return the timestep number.
    implicit none
    class(ecosim_time_type), intent(inout) :: this

    character(len=*), parameter :: sub = 'betr::get_nstep'
    integer, intent(in) :: nstep
    integer, intent(in) :: tstep
    logical, intent(in) :: continue_run
    real(r8) :: cdtime
    integer :: nn,cdays

    this%time0  = nstep*this%delta_time
    if(continue_run)then
      this%nelapstep = nstep
    else
      this%nelapstep = 0
    endif
    this%tstep=tstep
    if(this%its_a_new_year())then
      this%tstep = 0
    endif
    cdtime=this%nelapstep*this%delta_time
    this%cdays=int(cdtime/secpday)
    if(cdtime>0._r8)then
      if(this%year0>0)then        
        if(this%cdays>365)then
          nn=0
          do while(.true.)
            if(isLeap(this%year0+nn))then
              if(cdays>366)then
                nn=nn+1
                cdays=cdays-366
              else
                exit
              endif
            else
              nn=nn+1
              cdays=cdays-365
              if(cdays<365)exit
            endif
          enddo
          this%cyears=nn
          this%doy=cdays
          this%dow=getdow(this%year0+this%cyears,this%doy)
          do nn=1,12
            if(this%doy<=cdaz(nn))then
              this%moy=nn
              exit
            endif
          enddo
        endif
      else
        this%cyears=int(cdtime/secpyear)
        
        this%doy=mod(this%cdays,365)
        this%dow=mod(this%cdays,7)
        do nn=1,12
          if(this%doy<=cdaz(nn))then
            this%moy=nn
            exit
          endif
        enddo
        this%tod=mod(cdtime,secpday)        
        this%toy=mod(cdtime,secpyear)
      endif
    endif
  end subroutine set_time_offset

  !-------------------------------------------------------------------------------
  subroutine proc_nextstep(this)

    implicit none

    class(ecosim_time_type), intent(inout) :: this

    this%nelapstep = this%nelapstep + 1
  end subroutine proc_nextstep

  !-------------------------------------------------------------------------------
  subroutine proc_initstep(this)

    implicit none
    class(ecosim_time_type), intent(inout) :: this

    this%nelapstep = 0
  end subroutine proc_initstep

  !-------------------------------------------------------------------------------
  integer function get_days_per_year(this, offset)

    implicit none

    class(ecosim_time_type), intent(in) :: this
    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.

    ! Positive for future times, negative
    ! for previous times.

    ! remove unused dummy arg compiler warning
    if (offset > 0) continue
    if (this%tstep > 0) continue

    !hardwire at the moment
    get_days_per_year = 365
  end function get_days_per_year

  !-------------------------------------------------------------------------------

  subroutine print_cur_time(this)
  implicit none
  class(ecosim_time_type), intent(in) :: this

  print*,'time=',this%time
  print*,'tod =',this%tod
  print*,'dow =',this%dow
  print*,'dom =',this%dom
  print*,'doy =',this%doy
  print*,'moy =',this%moy
  print*,'cyears=',this%cyears

  end subroutine print_cur_time


  !-------------------------------------------------------------------------------
  subroutine get_ymdhs(this, ymdhs)

  implicit none
  class(ecosim_time_type), intent(in) :: this
  character(len=*), intent(out) :: ymdhs

  write(ymdhs,'(I4.4,I2.2,I2.2,I2.2,I4.4)')this%cyears+this%year0,this%moy,this%dom,&
     int(this%tod/3600.0),int(mod(this%tod,3600._r8))
  end subroutine get_ymdhs
  !-------------------------------------------------------------------------------
  function its_a_new_hour(this)result(yesno)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno


  yesno= abs(mod(this%tod, 3600._r8))<1.e-3_r8

  end function its_a_new_hour

  !-------------------------------------------------------------------------------

  function its_a_new_day(this)result(yesno)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno


  yesno= abs(mod(this%tod, 86400._r8))<1.e-3_r8

  end function its_a_new_day

  !-------------------------------------------------------------------------------

  function its_a_new_week(this)result(yesno)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno


  yesno= (this%dow == 0 .and. this%tod<1.e-3_r8)

  end function its_a_new_week

  !-------------------------------------------------------------------------------

  function its_a_new_month(this)result(yesno)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno

  if(this%moy/=2)then
    !not Feb
    yesno = ((this%dom == daz(this%moy)+1 .or. this%dom==1) .and. this%tod < 1.e-3_r8)
  else
    !Feb
    yesno = ((this%dom == daz(this%moy)+this%leap_yr+1 .or. this%dom==1) .and. this%tod < 1.e-3_r8)
  endif
  end function its_a_new_month

  !-------------------------------------------------------------------------------
  function its_a_new_year(this)result(yesno)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno

  yesno = (this%moy == 1 .and. this%dom == 1 .and. this%tod < 1.e-3_r8)

  end function its_a_new_year
  !-------------------------------------------------------------------------------

  function its_time_to_histflush(this)result(yesno)

  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: yesno

  if(this%hist_freq<0)then
    !by time steps
    yesno = (mod(this%nelapstep,this%hist_freq)==0)
  elseif(this%hist_freq==0)then
    !no history file output until the last time step
    yesno=this%its_time_to_exit()
  elseif(this%hist_freq==1)then
    !by day
    yesno=this%its_a_new_day()
  elseif(this%hist_freq==30)then
    !by month
    yesno=this%its_a_new_month()
  elseif(this%hist_freq==365)then
    !by year
    yesno=this%its_a_new_year()
  elseif(this%hist_freq==9999)then
    yesno=.false.
  endif
  end function its_time_to_histflush
  !-------------------------------------------------------------------------------
  function get_cur_year(this)result(ans)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer :: ans
  ans = this%cyears
  end function get_cur_year
  !-------------------------------------------------------------------------------
  function get_cur_yearAD(this)result(ans)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer :: ans
  ans = this%cyears+this%year0
  end function get_cur_yearAD
  !-------------------------------------------------------------------------------
  function get_cur_mon(this)result(ans)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer :: ans

  ans = this%moy  
  return
  end function get_cur_mon
  
  !-------------------------------------------------------------------------------
  function get_cur_day(this)result(ans)
  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer :: ans
  ans = this%cdays
  end function get_cur_day
  !-------------------------------------------------------------------------------
  function is_first_step(this)result(ans)
  !
  !is it the beginning of the first time step?
  implicit none
  class(ecosim_time_type), intent(in) :: this
  logical :: ans
  ans = (this%nelapstep==0)
  end function is_first_step

  !-------------------------------------------------------------------------------
  subroutine print_model_time_stamp(this, iulog)

  implicit none
  class(ecosim_time_type), intent(in) :: this
  integer, intent(in) :: iulog

  select case(this%stop_opt)
  case (1)
    write(iulog,*)'step', this%get_nstep()
  case (2)
    if (this%its_a_new_day()) then
      write(iulog,*)'day', this%get_cur_day()
    end if
  case (3)
    if (this%its_a_new_year()) then
      write(iulog,*)'year', this%get_cur_year()
    end if
  end select

  end subroutine print_model_time_stamp

!------------------------------------------------------------------------------------------

  function isLeap(year)result(ans)
!
! Description
! Determine if it is a leap year

  implicit none
  integer, intent(in) :: year
  logical :: ans

  ans =(mod(year,400)== 0) .or. (mod(year,4)==0 .and. mod(year,100)/=0)
  end function isLeap  
!------------------------------------------------------------------------------------------

  function getdow(year,doy)result(ans)
  !
  !get day of week
  implicit none
  integer, intent(in) :: year    !inquired year
  integer, intent(in) :: doy     !inquired day of year

  integer :: ans,nd,nn

  nd=0
  if (year<1900)then  
    do nn=1900,year,-1
      if(isLeap(nn))then
        nd=mod(-366+nd,7)
      else
        nd=mod(-365+nd,7)
      endif
    enddo
    ans=mod(nd+doy+1,7)
    if(ans==0)ans=7
    if(ans<0)ans=ans+7
  elseif (year==1900)then
    ans=mod(doy,7)
    if(ans==0)ans=7
  else  
    do nn=1900,year-1
      if(isLeap(nn))then
        nd=mod(366+nd,7)
      else
        nd=mod(365+nd,7)
      endif
    enddo
    ans=mod(nd+doy,7)
    if(ans==0)ans=7
  endif
  end function getdow
end module ecosim_Time_Mod
