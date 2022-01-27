module minimathmod
!!
! Description:
! Some small subroutines/function to do safe math.

  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none

  private
  public :: safe_adb
  real(r8), parameter :: tiny_val=1.e-20_r8
  contains

   function safe_adb(a,b)result(ans)
   !!
   ! Description:
   ! damp division by zero to zero
   implicit none
   real(r8), intent(in) :: a,b
   real(r8) :: ans

   ans = a*b/(b*b+tiny_val)
   
   return
   end function safe_adb

end module minimathmod
