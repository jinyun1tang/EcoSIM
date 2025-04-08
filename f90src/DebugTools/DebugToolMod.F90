module DebugToolMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use fileUtil,      only: iulog
  use EcoSIMCtrlMod, only: lverb
implicit none
  private
  character(len=*), parameter :: mod_filename =&
   __FILE__  

  public :: PrintInfo
  public :: DebugPrint

  interface DebugPrint
    module procedure DebugPrint_real_arr
    module procedure DebugPrint_real_sp
    module procedure DebugPrint_int
  end interface DebugPrint

contains

  !-----------------------------------------------------------------------
    subroutine PrintInfo(message)
    implicit none
    character(len=*), intent(in) :: message

    if(lverb)write(iulog,*)message
    end subroutine PrintInfo

  !-----------------------------------------------------------------------

   subroutine DebugPrint_real_arr(vnames,vals)
   implicit none
   character(len=32), intent(in) :: vnames(:)
   real(r8), intent(in) :: vals(:)
   integer :: nvars
   character(len=256) :: fmt
   nvars=size(vals)
   
   if(.not.lverb)return

   write(fmt,'(A,I2)')'(A,A,',nvars   
   fmt=trim(fmt)//'E14.6)'   
   write(*,fmt)vnames,'=',vals

   end subroutine DebugPrint_real_arr

  !-----------------------------------------------------------------------

   subroutine DebugPrint_real_sp(vname,val)
   implicit none
   character(len=*), intent(in) :: vname
   real(r8), intent(in) :: val
   
   if(.not.lverb)return
   write(iulog,*)vname,'=',val
   end subroutine DebugPrint_real_sp
  !-----------------------------------------------------------------------

   subroutine DebugPrint_int(vname,val)
   implicit none
   character(len=*), intent(in) :: vname
   integer, intent(in) :: val
   
   if(.not.lverb)return
   write(iulog,*)vname,'=',val
   end subroutine DebugPrint_int

end module DebugToolMod