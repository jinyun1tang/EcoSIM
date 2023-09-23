module minimathmod
!!
! Description:
! Some small subroutines/function to do safe math.

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSimConst
  implicit none
  character(len=*),private, parameter :: mod_filename = __FILE__
  private
  public :: safe_adb
  public :: p_adb
  public :: isclose         !test if two values a and b are close in magnitude
  public :: vapsat, vapsat0
  public :: isLeap
  public :: isnan
  public :: AZMAX1,AZMIN1,AZMAX1t
  public :: GetMolAirPerm3
  interface AZMAX1
    module procedure AZMAX1_s
    module procedure AZMAX1_d
  end interface AZMAX1

  interface AZMIN1
    module procedure AZMIN1_s
    module procedure AZMIN1_d
  end interface AZMIN1

  public :: addone
  public :: RichardsonNumber
  real(r8), parameter :: tiny_val=1.e-20_r8
  contains

   pure function isnan(a)result(ans)
   implicit none
   real(r8), intent(in) :: a
   logical :: ans

   ans=(a/=a)
   return
   end function isnan
!------------------------------------------------------------------------------------------

   pure function safe_adb(a,b)result(ans)
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

   pure function p_adb(a,b)result(ans)
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

  pure function vapsat(tempK)result(ans)
  !
  ! Description
  ! compute saturated vapor pressure, based on temperature tempK (in K)
  implicit none
  real(r8), intent(in) :: tempK

  real(r8) :: ans  !ton/m3, i.e. (ans*10^3=kg/m3) in terms vapor concentration, 2.173~18/8.314
  ans=2.173E-03_r8/tempK*0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/tempK))
  end function vapsat

!------------------------------------------------------------------------------------------

  pure function vapsat0(tempK)result(ans)
  !
  ! Description
  ! compute saturated vapor pressure, based on temperature tempK (in K)
  implicit none
  real(r8), intent(in) :: tempK

  real(r8) :: ans  !(kPa)
  ans=0.61_r8*EXP(5360.0_r8*(3.661E-03_r8-1.0_r8/tempK))
  end function vapsat0

!------------------------------------------------------------------------------------------

  pure function isLeap(year)result(ans)
!
! Description
! Determine if it is a leap year

  implicit none
  integer, intent(in) :: year
  logical :: ans

  ans =(mod(year,400)== 0) .or. (mod(year,4)==0 .and. mod(year,100)/=0)
  end function isLeap
!------------------------------------------------------------------------------------------

  pure function AZMAX1t(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  ans=AMAX1(val,tiny_val)

  end function AZMAX1t
!------------------------------------------------------------------------------------------

  pure function AZMAX1_s(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  ans=AMAX1(0.0_r8,val)

  end function AZMAX1_s  
!------------------------------------------------------------------------------------------

  pure function AZMAX1_d(val1,val2)result(ans)
  implicit none
  real(r8), intent(in) :: val1,val2

  real(r8) :: ans

  ans=AMAX1(0.0_r8,val1,val2)

  end function AZMAX1_d

!------------------------------------------------------------------------------------------

  pure function AZMIN1_s(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  ans=AMIN1(0.0_r8,val)

  end function AZMIN1_s


!------------------------------------------------------------------------------------------

  pure function AZMIN1_d(val1,val2)result(ans)
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


! ----------------------------------------------------------------------
  pure function isclose(a,b)result(ans)
  !DESCRIPTION
  !determine if a is close to b in magnitude by relative magnitude tiny_val

  implicit none
  real(r8), intent(in) :: a,b
  real(r8) :: c,ac,bc
  logical :: ans
  
  c=max(abs(a),(b))  
  if (c==0._r8) then
    ans=.True.
    return 
  endif

  ac=a/c;bc=b/c  
  ans=abs((ac-bc)/(ac+bc))<tiny_val
  end function isclose

! ----------------------------------------------------------------------

  pure function RichardsonNumber(RIB,TK1,TK2)result(ans)
  implicit none
  real(r8), intent(in) :: RIB  !isothermal RI
  real(r8), intent(in) :: TK1, TK2

  real(r8) :: ans

  ans = AMAX1(-0.3_r8,AMIN1(0.075_r8,RIB*(TK1-TK2)))
  
  end function RichardsonNumber
! ----------------------------------------------------------------------
  pure function GetMolAirPerm3(TKair,Patm_Pa)result(ans)
  implicit none
  real(r8), intent(in) :: TKair
  real(r8), optional, intent(in) :: Patm_Pa  !atmospheric pressure in Pascal

  real(r8) :: ans

  if(present(Patm_Pa))then
    ans = Patm_Pa/(TKair*RGAS)
  else
    ans= 1.01325E5_r8/(TKair*RGAS)
  endif

  end function GetMolAirPerm3


end module minimathmod
