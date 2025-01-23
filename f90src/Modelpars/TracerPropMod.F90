module TracerPropMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use ChemTracerParsMod
  use abortutils, only : endrun
  use TracerIDMod
implicit none
  private

  character(len=*), parameter :: mod_filename = &
  __FILE__


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
  case (idg_CO2)
    coef=SCO2X*EXP(0.843_R8-0.0281_R8*tempC)
  case (idg_CH4)
    coef=SCH4X*EXP(0.597_r8-0.0199_r8*tempC)
  case (idg_O2)
    coef=SOXYX*EXP(0.516_r8-0.0172_r8*tempC)
  case (idg_NH3)
    coef=SNH3X*EXP(0.513_r8-0.0171_r8*tempC)
  case (idg_N2O)
    coef=SN2OX*EXP(0.897_r8-0.0299_r8*tempC)
  case (idg_N2)
    coef=SN2GX*EXP(0.456_r8-0.0152_r8*tempC)
  case (idg_H2)
    coef=SH2GX*EXP(0.597_r8-0.0199_r8*tempC)
  case default
    call endrun("tracer not defined in "//mod_filename,__LINE__)
  end select
  end function gas_solubility

end module TracerPropMod
