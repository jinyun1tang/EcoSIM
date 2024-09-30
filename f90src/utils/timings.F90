module timings
!!
! Description:
! code to do runtime timing
!
  use data_kind_mod, only : r8 => DAT_KIND_R8

  implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename =&
   __FILE__  
  integer,parameter :: varlen=30
  integer, parameter :: maxprocs=200
  integer :: timer_number
  integer :: loop_number
  integer :: lun
  logical :: init_flag
  real(r8) :: timer_array(maxprocs)
  character(len=varlen) :: name_array(maxprocs)
  public :: init_timer
  public :: start_timer
  public :: end_timer
  public :: end_timer_loop
  contains

    subroutine init_timer(outdir)
!!
!   Description:
!
    implicit none
    character(len=*), intent(in) :: outdir
    init_flag=.false.
    timer_number=0
    loop_number=0
    lun=3001
    timer_array(:)=0

    call system('mkdir -p '//trim(outdir)//'/timing' )
    OPEN(UNIT=lun, FILE=trim(outdir)//"/timing/time.txt",STATUS='UNKNOWN')

    end subroutine init_timer
! ----------------------------------------------------------------------

    subroutine start_timer(t1)
!!
!   Description:
!   This is a simple subroutine to start a timer
    implicit none
    real(r8), intent(out) :: t1

!    call timec(t1)

    call timer(.true., t1)

    END subroutine start_timer
! ----------------------------------------------------------------------

    subroutine end_timer(timer_name,t1)
!!
!   Description:
    implicit none
    real(r8), intent(in) :: t1
    character(len=*), intent(in) :: timer_name

    real(r8) :: elapsed_time
    real(r8) :: t2
!   calculate the time between start_timer and end_timer
!   and add it to the timer_array
!    call timec(t2)
    call timer(.false., t2)
    elapsed_time = (t2 - t1)

    timer_array(timer_number+1) = elapsed_time
    name_array(timer_number+1) = trim(timer_name)
    timer_number = timer_number+1

    END subroutine end_timer
! ----------------------------------------------------------------------

    subroutine end_timer_loop()
!!
!   Description:
    implicit none

    character(len=4) :: fmt1
    character(len=30) :: fmt

    integer :: j

    if (loop_number.EQ.0)then
      write(fmt1, '(I3)')timer_number
      fmt='(A20,'//trim(fmt1)//'(A,A20))'
      write(lun,fmt)'steps',(',',name_array(j),j=1,timer_number)
      init_flag =.true.
    end if

    write(fmt1,'(I3)')timer_number
    fmt='(I20,'//trim(fmt1)//'(A,E30.10))'
    write(lun,fmt)loop_number,(',',timer_array(j),j=1,timer_number)

    timer_array=0.
    loop_number = loop_number+1
    timer_number = 0

    end subroutine end_timer_loop

! ----------------------------------------------------------------------

  subroutine timer(start, time)
    implicit none
    logical, intent(in) :: start
    real(r8), intent(out) :: time
    integer :: count, count_rate, count_max
    save count, count_rate, count_max
    integer :: new_count

    if (start) then
      call system_clock(count, count_rate, count_max)
      time=0._r8
    else
      call system_clock(new_count, count_rate, count_max)
      time = real(new_count - count,kind=r8) / real(count_rate,kind=r8)
    end if
  end subroutine timer

END module timings
