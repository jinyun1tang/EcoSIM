module TracerPropMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use ChemTracerParsMod
  use abortutils, only : endrun
implicit none
  private
  
  character(len=*), parameter :: mod_filename = __FILE__

  integer, public, parameter :: id_co2g=1
  integer, public, parameter :: id_ch4g=2
  integer, public, parameter :: id_o2g =3
  integer, public, parameter :: id_nh3g=4
  integer, public, parameter :: id_n2g =5
  integer, public, parameter :: id_n2og=6
  integer, public, parameter :: id_h2g =7
  public :: gas_solubility
  contains


  function gas_solubility(gid,tempC)result(coef)
!
! DESCRIPTION
! compute gas solubility

  implicit none
  integer :: gid  !gas id
  real(r8), intent(in) :: tempC
  real(r8) :: coef     !defined based on molar concentration
  select case (gid)
  case (id_co2g)
    coef=SCO2X*EXP(0.843_R8-0.0281_R8*tempC)
  case (id_ch4g)
    coef=SCH4X*EXP(0.597_r8-0.0199_r8*tempC)
  case (id_o2g)
    coef=SOXYX*EXP(0.516_r8-0.0172_r8*tempC)
  case (id_nh3g)
    coef=SNH3X*EXP(0.513_r8-0.0171_r8*tempC)
  case (id_n2og)
    coef=SN2OX*EXP(0.897_r8-0.0299_r8*tempC)
  case (id_n2g)
    coef=SN2GX*EXP(0.456_r8-0.0152_r8*tempC)
  case (id_h2g)
    coef=SH2GX*EXP(0.597_r8-0.0199_r8*tempC)
  case default
    call endrun("tracer not defined in "//mod_filename,__LINE__)
  end select
  end function gas_solubility

end module TracerPropMod
