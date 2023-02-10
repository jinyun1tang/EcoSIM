module minimathmod
!!
! Description:
! Some small subroutines/function to do safe math.

  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  character(len=*),private, parameter :: mod_filename = __FILE__
  private
  public :: safe_adb
  public :: p_adb
  public :: test_aeqb     !a equals b test, with precision tiny_val
  public :: test_aneb     !a not equals b test
  public :: vapsat, vapsat0
  public :: isLeap
  public :: AZMAX1,AZMIN1

  interface AZMAX1
    module procedure AZMAX1_s
    module procedure AZMAX1_d
  end interface AZMAX1

  interface AZMIN1
    module procedure AZMIN1_s
    module procedure AZMIN1_d
  end interface AZMIN1

  public :: addone
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


!------------------------------------------------------------------------------------------

  function vapsat(tempK)result(ans)
  !
  ! Description
  ! compute saturated vapor pressure, based on temperature tempK (in K)
  implicit none
  real(r8), intent(in) :: tempK

  real(r8) :: ans  !(10^3 kg/m3)
  ans=2.173E-03_r8/tempK*0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/tempK))
  end function vapsat

!------------------------------------------------------------------------------------------

  function vapsat0(tempK)result(ans)
  !
  ! Description
  ! compute saturated vapor pressure, based on temperature tempK (in K)
  implicit none
  real(r8), intent(in) :: tempK

  real(r8) :: ans  !(kPa)
  ans=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/tempK))
  end function vapsat0

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

  function AZMAX1_s(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  ans=AMAX1(0.0_r8,val)

  end function AZMAX1_s
!------------------------------------------------------------------------------------------

  function AZMAX1_d(val1,val2)result(ans)
  implicit none
  real(r8), intent(in) :: val1,val2

  real(r8) :: ans

  ans=AMAX1(0.0_r8,val1,val2)

  end function AZMAX1_d

!------------------------------------------------------------------------------------------

  function AZMIN1_s(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  ans=AMIN1(0.0_r8,val)

  end function AZMIN1_s


!------------------------------------------------------------------------------------------

  function AZMIN1_d(val1,val2)result(ans)
  implicit none
  real(r8), intent(in) :: val1,val2

  real(r8) :: ans

  ans=AMIN1(0.0_r8,val1,val2)

  end function AZMIN1_d

! ----------------------------------------------------------------------

  function addone(itemp)result(ans)
!
!  DESCRIPTION
! increase itemp by one
  implicit none
  integer, intent(inout) :: itemp

  integer :: ans

  itemp=itemp+1
  ans=itemp
  end function addone

end module minimathmod
