module Hour1Pars
  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  public
  save
  character(len=*),private, parameter :: mod_filename = __FILE__


  real(r8) :: BKDSX  !     BKDSX=maximm soil bulk density
  real(r8) :: THETPW !     THETPW=minimum air-filled porosity for saturation (m3 m-3)
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

  THETPW=0.01_r8
  THETWP=1.0_r8-THETPW

  XVOLWC=real((/5.0E-04,2.5E-04,2.5E-04,2.5E-04/),r8)
  THETRX=real((/4.0E-06,8.0E-06,8.0E-06/),r8)

  end subroutine initHour1Pars

end module Hour1Pars
