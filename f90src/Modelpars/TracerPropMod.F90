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
  public :: GramPerHr2umolPerSec
  public :: GasSechenovConst
  public :: MolecularWeight
  contains

!------------------------------------------------------------------------------------------

  function gas_solubility(gid,tempC)result(coef)
!
! DESCRIPTION
! compute gas solubility

  implicit none
  integer :: gid  !gas id
  real(r8), intent(in) :: tempC   !temperature in celcius degree
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
  case (idg_AR)
    coef=SARX*EXP(-1500._r8*(1./(tempC+273.15)-1./298.15_r8))  
  case default
    call endrun("tracer not defined in "//mod_filename,__LINE__)
  end select
  end function gas_solubility

!------------------------------------------------------------------------------------------

  function GramPerHr2umolPerSec(gid)result(coef)
  implicit none
  integer, intent(in) :: gid

  real(r8) :: coef
  real(r8),parameter :: gH1hour2umol1sec=1.e6_r8/3600._r8

  select case (gid)
  case (idg_CO2,idg_CH4)
    coef=gH1hour2umol1sec/12._r8
  case (idg_O2)
    coef=gH1hour2umol1sec/32._r8
  case (idg_NH3)
    coef=gH1hour2umol1sec/14._r8
  case (idg_N2O,idg_N2)
    coef=gH1hour2umol1sec/28._r8
  case (idg_H2)
    coef=gH1hour2umol1sec/2._r8
  case (idg_AR)
    coef=gH1hour2umol1sec/39.95
  case default
    call endrun("tracer not defined in "//mod_filename,__LINE__)
  end select
  end function GramPerHr2umolPerSec
!------------------------------------------------------------------------------------------

  function GasSechenovConst(gid)result(ans)
  !
  !The temperature dependence of the Sechenov Constants is not considered here
  implicit none
  integer, intent(in) :: gid

  real(r8) :: ans

  select case(gid)
  case (idg_CO2,idg_CH4)
    ans = 0.14_r8
  case (idg_O2)
    ans = 0.31_r8
  case (idg_N2, idg_N2O)
    ans = 0.23_r8
  case  (idg_NH3)
    ans =0.07_r8
  case (idg_H2)
    ans =0.14_r8
  case (idg_AR)
    ans =0.06_r8 
  case default
    call endrun("tracer not defined in "//mod_filename,__LINE__)
  end select            
  end function GasSechenovConst

!------------------------------------------------------------------------------------------

  function MolecularWeight(gid)result(ans)  

  implicit none
  integer, intent(in) :: gid

  real(r8) :: ans

  select case(gid)
  case (idg_CO2,idg_CH4)  !C-based
    ans = 12._r8
  case (idg_O2)           !O2-based
    ans = 32._r8
  case (idg_N2, idg_N2O)  !N-based
    ans = 28._r8
  case  (idg_NH3)         !N-based
    ans =14._r8
  case (idg_H2)  
    ans =2._r8
  case (idg_AR)
    ans =39.95_r8
  case default
    call endrun("tracer not defined in "//mod_filename,__LINE__)
  end select            
  end function MolecularWeight
end module TracerPropMod
