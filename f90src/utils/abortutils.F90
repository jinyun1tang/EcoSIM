module abortutils

  !-----------------------------------------------------------------------
  ! !MODULE: abortutils
  !
  ! !DESCRIPTION:
  ! Abort the model for abnormal termination
  !-----------------------------------------------------------------------
  use data_kind_mod, only : r8 => SHR_KIND_R8

  implicit none
  private
  save

  public :: endrun
  interface endrun
     module procedure endrun_vanilla
     module procedure endrun_line
     module procedure endrun_globalindex
  end interface

  public :: padl, padr
  public :: print_info
  interface print_info
    module procedure print_info_arr
    module procedure print_info_msg
  end interface
  integer  :: iulog = 6        ! "stdout" log file unit number, default is 6

CONTAINS
  !-----------------------------------------------------------------------
  subroutine endrun_vanilla(msg)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    !
    !
    ! !ARGUMENTS:

    implicit none
    character(len=*), intent(in), optional :: msg    ! string to be printed
    !-----------------------------------------------------------------------

    if (present (msg)) then
       write(iulog,*)'ENDRUN:', msg
    else
       write(iulog,*)'ENDRUN: called without a message string'
    end if

    stop

  end subroutine endrun_vanilla

  !-----------------------------------------------------------------------
  subroutine endrun_line(msg,line)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    !
    !
    ! !ARGUMENTS:

    implicit none
    character(len=*), intent(in) :: msg    ! string to be printed
    integer, intent(in) :: line
    !-----------------------------------------------------------------------

    write(iulog,*)'ENDRUN:', msg, 'at line ',line
    stop

  end subroutine endrun_line

  !-----------------------------------------------------------------------
  subroutine endrun_globalindex(decomp_index, elmlevel, msg)

    !-----------------------------------------------------------------------
    ! Description:
    ! Abort the model for abnormal termination
    !

    !
    ! Arguments:
    implicit none
    integer          , intent(in)           :: decomp_index
    character(len=*) , intent(in)           :: elmlevel
    character(len=*) , intent(in), optional :: msg    ! string to be printed
    !
    ! Local Variables:
    integer :: igrc, ilun, icol
    !-----------------------------------------------------------------------

    write(6,*)'calling getglobalwrite with decomp_index= ',decomp_index,' and elmlevel= ',trim(elmlevel)


    if (present (msg)) then
       write(iulog,*)'ENDRUN:', msg
    else
       write(iulog,*)'ENDRUN: called without a message string'
    end if

    stop

  end subroutine endrun_globalindex

  !-----------------------------------------------------------------------

  function padl(str,width,symbol) result(padded_str)
    !-----------------------------------------------------------------------
    ! Description:
    ! pad a string from left with symbol to width
    !
    !
    ! !ARGUMENTS:

    implicit none
    character(len=*), intent(in) :: str
    integer, intent(in) :: width
    character(len=1), intent(in), optional :: symbol
    character(len=max(len(str),width)) :: padded_str

    integer :: ns, jj, kk, kk1
    character(len=1) :: symboll

    if(present(symbol))then
      symboll=symbol
    else
      symboll=' '
    endif
    ns = len(str)
    if (ns >= width)then
      padded_str=str
    else
      do jj = width, 1, -1
        kk = width-jj
        if (kk <= ns-1) then
          padded_str(jj:jj) = str(ns-kk:ns-kk)
        else
          padded_str(jj:jj) = symboll
        endif
      enddo
    endif
    return
  end function padl
  !-----------------------------------------------------------------------

  ! right - padding
  function padr(str,width,symbol) result(padded_str)
    !-----------------------------------------------------------------------
    ! Description:
    ! pad a string from right with symbol to width
    !
    !
    ! !ARGUMENTS:

    implicit none
    character(len=*), intent(in) :: str
    integer, intent(in) :: width
    character(len=1), intent(in), optional :: symbol
    character(len=max(len(str),width)) :: padded_str

    integer :: ns, jj
    character(len=1) :: symboll

    if(present(symbol))then
      symboll=symbol
    else
      symboll=' '
    endif
    ns = len(str)
    if (ns >= width)then
      padded_str=str
    else
      do jj = 1, width
        if (jj <= ns) then
          padded_str(jj:jj) = str(jj:jj)
        else
          padded_str(jj:jj) = symboll
        endif
      enddo
    endif
    return

  end function padr

  !-----------------------------------------------------------------------

  subroutine print_info_arr(msg,strarr,valarr)

    !-----------------------------------------------------------------------
    ! Description:
    ! print information
    !
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: msg
    character(len=*), dimension(:), intent(in) :: strarr(:)
    real(r8), dimension(:), intent(in) :: valarr

    integer :: ns1
    integer :: jj

    ns1=size(valarr)
    write(*,*)trim(msg)
    do jj = 1 , ns1
      write(*,*)trim(strarr(jj))//'=',valarr(jj)
    enddo
  end subroutine print_info_arr

  !-----------------------------------------------------------------------
  subroutine print_info_msg(msg,lineno)

  implicit none
  character(len=*), intent(in) :: msg
  integer, intent(in) :: lineno

  write(*,*)''
  write(*,*)trim(msg)
  write(*,*)'at line',lineno
  end subroutine print_info_msg
end module abortutils
