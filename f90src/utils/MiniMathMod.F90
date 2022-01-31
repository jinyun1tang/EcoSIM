module minimathmod
!!
! Description:
! Some small subroutines/function to do safe math.

  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none

  private
  public :: safe_adb
  public :: p_adb
  public :: test_aeqb
  public :: test_aneb
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
!------------------------------------------------------------------------------------------

   function p_adb(a,b)result(ans)
   !!
   ! Description:
   ! ans=max(0.,a/b)
   implicit none
   real(r8), intent(in) :: a,b
   real(r8) :: ans

   ans=AMAX1(0._r8,a/b)
   return
   end function p_adb
!------------------------------------------------------------------------------------------
   function test_aeqb(a,b)result(ans)
   !!
   ! Description:
   ! a and b are equal within difference tiny_val
   ! ans = abs(a-b)<tiny_val
   implicit none
   real(r8), intent(in) :: a,b
   logical :: ans

   ans = abs(a-b)<tiny_val
   return
   end function test_aeqb

!------------------------------------------------------------------------------------------
   function test_aneb(a,b)result(ans)
   !!
   ! Description:
   ! a and b are not equal within difference tiny_val
   ! ans = abs(a-b)>=tiny_val
   implicit none
   real(r8), intent(in) :: a,b
   logical :: ans

   ans = abs(a-b)>=tiny_val
   return
   end function test_aneb

end module minimathmod
