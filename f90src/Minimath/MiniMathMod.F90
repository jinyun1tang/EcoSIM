module minimathmod
!!
! Description:
! Some small subroutines/function to do safe math.

  use data_kind_mod, only: r8 => DAT_KIND_R8
  use fileUtil,      only: iulog
  use EcoSimConst

  implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: safe_adb
  public :: p_adb
  public :: isclose         !test if two values a and b are close in magnitude
  public :: vapsat, vapsat0
  public :: isLeap
  public :: isnan
  public :: AZMAX1,AZMIN1,AZMAX1t,AZMAX1d,AZMIN1d
  public :: GetMolAirPerm3
  public :: fSiLU
  public :: fixnegmass
  public :: fixEXConsumpFlux
  public :: yearday,isletter
  public :: dssign
  public :: flux_mass_limiter
  public :: AZERO
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
  real(r8), parameter :: tiny_val=1.e-14_r8

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

   if(abs(b)<tiny_val)then
     ans=tiny_val
   else
     ans = a/b
   endif

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

  function dssign(snow)result(ans)
  implicit none
  real(r8), intent(in) :: snow
  real(r8) :: ans

  if(snow<=1.e-8_r8)then
    ans=0._r8
  else
    ans=1._r8
  endif
  end function dssign

!------------------------------------------------------------------------------------------

  pure function vapsat(tempK)result(ans)
  !
  ! Description
  ! compute saturated vapor pressure, based on temperature tempK (in K)
  ! ep=\mu R T= (m/M)/V0 RT, V0=1m3, m=ep*V0*M/(RT), [kPa]*[m3]*[g/mol]/[Pa m3]~kg 
  implicit none
  real(r8), intent(in) :: tempK

  real(r8) :: ans  !ton, i.e. (ans*10^3=kg/m3) in terms vapor concentration, 2.173~18/8.314
  
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

  pure function iisLeap(year)result(ans)
!
! Description
! Determine if it is a leap year

  implicit none
  integer, intent(in) :: year
  integer :: ans

  if (isLeap(year))then
    ans=1
  else
    ans=0
  endif
  end function iisLeap

!------------------------------------------------------------------------------------------

  pure function AZMAX1t(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  ans=AMAX1(val,tiny_val)

  end function AZMAX1t
!------------------------------------------------------------------------------------------

  pure function AZMAX1d(val,tiny_val2)result(ans)
  !
  implicit none
  real(r8), intent(in) :: val
  real(r8), intent(in) :: tiny_val2 !positive tiny threshold value

  real(r8) :: ans

  if(val>tiny_val2)then
    ans=val 
  else
    ans=0._r8
  endif
  end function AZMAX1d

!------------------------------------------------------------------------------------------

  pure function AZMIN1d(val,tiny_val2)result(ans)
  implicit none
  real(r8), intent(in) :: val
  real(r8), intent(in) :: tiny_val2 !positive tiny threshold value

  real(r8) :: ans

  if(val<-tiny_val2)then
    ans=val 
  else
    ans=0._r8
  endif
  end function AZMIN1d

!------------------------------------------------------------------------------------------
  pure function AZERO(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  if(abs(val)>=tiny_val)then
    ans=val
  else  
    ans=0._r8
  endif

  end function AZERO
!------------------------------------------------------------------------------------------

  pure function AZMAX1_s(val)result(ans)
  implicit none
  real(r8), intent(in) :: val

  real(r8) :: ans

  if(val>=tiny_val)then
    ans=val
  else  
    ans=0._r8
  endif

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
  
  c=max(abs(a),abs(b))  
  if (c<tiny_val) then
    ans=.True.
    return     
  endif

  if(abs(a+b)<1.e-20_r8)then
    ans=.false.
  else  
    ac=a/c;bc=b/c    
    ans=abs((ac-bc)/(ac+bc))<tiny_val
  endif
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
    ans = Patm_Pa/(TKair*RGASC)
  else
    ans= 1.01325E5_r8/(TKair*RGASC)
  endif

  end function GetMolAirPerm3
! ----------------------------------------------------------------------
  function fSiLU(x,b)result(ans)
  !Sigmoid linear unit function
  implicit none
  real(r8), intent(in) :: x
  real(r8),optional, intent(in) :: b
  real(r8) :: b_loc
  
  real(r8) :: ans
  
  if(present(b))then
    b_loc=b
  else
    b_loc=1._r8
  endif
  
  ans=x/(1._r8+exp(-b_loc*x))
  
  end function fSiLU
! ----------------------------------------------------------------------
  function fixnegmass(val,refcon)result(ans)
  implicit none
  real(r8), intent(in) :: val
  real(r8), optional, intent(in) :: refcon
  real(r8) :: ans

  ans=val

  if(present(refcon))then
    if(val<0._r8 .and. val>AMIN1(-refcon*1.e-3_r8,-1.e-5))ans=0._r8
  else
    if(val<0._r8 .and. val >-1.e-5_r8)ans=0._r8
  endif

  end function fixnegmass
! ----------------------------------------------------------------------

  subroutine fixEXConsumpFlux(mass,consum_flux,dsgn)
  implicit none
  real(r8), intent(inout) :: mass
  real(r8), intent(inout) :: consum_flux
  integer, optional, intent(in) :: dsgn
  integer :: dsgnl

  dsgnl=1
  if(present(dsgn))dsgnl=dsgn


  !cut off too small mass
  if(abs(mass)<1.e-12_r8)mass=0._r8    
  !return for zero flux
  if(isclose(consum_flux,0._r8))return

  !udpate as mass=mass-flux
  if(dsgnl>0)then  
    if(mass<consum_flux)then
      consum_flux = mass
      mass        = 0._r8
    else
      mass=mass-consum_flux  
    endif
  !update as mass=mass+flux
  else  
    if(mass<-consum_flux)then  
      consum_flux = -mass
      mass        = 0._r8
    else
      mass=mass+consum_flux
    endif
  endif
  end subroutine fixEXConsumpFlux

! ----------------------------------------------------------------------

  pure function yearday(year,month,day)result(doy)
  implicit none
  integer, intent(in) :: year
  integer, intent(in) :: month
  integer, intent(in) :: day

  integer, parameter :: daz(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  integer :: doy
  integer :: jj

  doy=0
  do jj = 1, month-1
    doy=doy+daz(jj)
  enddo
  if(month>2)then
    doy=doy+iisleap(year);
  endif
  doy=doy+day
  end function yearday
! ----------------------------------------------------------------------

  pure function isletter(c)result(ans)
  implicit none
  character(len=1), intent(in) :: c
  logical :: ans

  ans=(c>='a' .and. c<='z') .or. (c>='A' .and. c<='Z')

  end function isletter

! ----------------------------------------------------------------------

  function flux_mass_limiter(flux,massa,massb)result(ans)
  !
  !limit flux by massa and massb 
  !assuming
  !massa=massa+flux
  !massb=massb-flux
  implicit none
  real(r8), intent(in) :: flux
  real(r8), intent(in) :: massa
  real(r8), intent(in) :: massb
  
  real(r8) :: ans
  if(flux>0._r8)then
    ans=AMIN1(flux,massb)-tiny_val
  else 
    ans=-AMIN1(-flux,massa)+tiny_val
  endif
  end function flux_mass_limiter
end module minimathmod
