module Hour1Pars
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  public
  save
!     BKDSX=maximm soil bulk density
!     THETPW=minimum air-filled porosity for saturation (m3 m-3)
  real(r8) :: BKDSX
  real(r8) :: FORGW
  real(r8) :: DTHETW  !difference between saturation and effective saturation
  real(r8) :: THETPW
  real(r8) :: THETWP
  real(r8) :: XVOLWC(0:3),THETRX(0:2)
!
!     XVOLWC=foliar water retention capacity (m3 m-2)
!     THETRX=litter water retention capacity (m3 g C-1)
!

  contains

  subroutine initHour1Pars
  implicit none

  BKDSX=1.89_r8
  FORGW=0.25E+06_r8 !threshold for  C concentration in organic soil 	g Mg-1
  DTHETW=1.0E-06_r8
  THETPW=0.01_r8
  THETWP=1.0_r8-THETPW

  XVOLWC=real((/5.0E-04,2.5E-04,2.5E-04,2.5E-04/),r8)
  THETRX=real((/4.0E-06,8.0E-06,8.0E-06/),r8)

  end subroutine initHour1Pars

end module Hour1Pars
