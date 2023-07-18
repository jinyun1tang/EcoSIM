module ATSCPLMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SharedDataMod
  use ATSEcoSIMInitMod
  implicit none

  public 
  character(len=*), private, parameter :: mod_filename=__FILE__


contains
!------------------------------------------------------------------------------------------

  subroutine ATS2EcoSIMData(filter_col,data_1d,var_1d,data_2d,var_2d,data_3d,var_3d)
  !send data from ATS to ecosim
  implicit none
  integer, intent(in) :: filter_col(:)
  real(r8), optional, intent(in) :: data_1d(:)              !1:nvar
  character(len=*), optional, intent(in) :: var_1d(:)       !1:nvar
  real(r8), optional,intent(in) :: data_2d(:,:)             !1:nvar,1:ncol, column specific scalar
  character(len=*), optional,intent(in) :: var_2d(:)        !1:nvar
  real(r8), optional, intent(in) :: data_3d(:,:,:)          !1:jz, 1:nvar,1:ncol, 1D vector column specific
  character(len=*), optional, intent(in) :: var_3d(:)       !
  character(len=*), parameter :: subname=trim(mod_filename)//'::ATS2EcoSIMData'
  integer :: ncol,nvar
  integer :: j1,j2,j3

  ncol=size(filter_col)
  
  if (present(data_1d) .and. present(var_1d))then
  !domain specific scalar
    nvar=size(var_1d)
    do j1=1,nvar
      select case(var_1d(j1))
      case ('ATM_N2')  !ppmv
        atm_N2=data_1d(j1)
      case ('ATM_O2')  !ppmv    
        atm_o2=data_1d(j1)
      case ('ATM_CO2') !ppmv    
        atm_co2=data_1d(j1)
      case ('ATM_CH4') !ppmv    
        atm_ch4=data_1d(j1)
      case ('ATM_N2O') !ppmv    
        atm_n2o=data_1d(j1)
      case ('ATM_H2')  !ppmv    
        atm_h2=data_1d(j1)
      case ('ATM_NH3') !ppmv    
        atm_NH3=data_1d(j1)
      end select
    enddo
  endif

  if (present(data_2d) .and. present(var_2d))then  
  !columun specific scalar 
    nvar=size(var_2d)
    do j2=1,ncol
      do j1=1,nvar
        select case(var_2d(j1))
        case ('TAIRC')     !air temperature, oC
          tairc(j2)=data_2d(j1,j2)
        case ('PREC')      !precipitation, mm H2O/hr
          prec(j2)=data_2d(j1,j2)
        case ('WINDH')     !horizontal wind speed,   m/s
          uwind(j2)=data_2d(j1,j2)
        case ('DWPTH')     !atmospheric vapor pressure, kPa
          vpa(j2)=data_2d(j1,j2)
        case ('SRADH')     !Incident solar radiation, W/m2 
          srad(j2)=data_2d(j1,j2) 
        end select
      enddo
    enddo
  endif

  if (present(data_2d) .and. present(var_2d))then  
  !1D vertical vector, 
    nvar=size(var_2d)
    do j3=1,ncol
      do j2=1,nvar
        select case (var_2d(j2))
        case ('CSAND')  !g/kg soil
          a_csand(1:JZSOI,j3)=data_3d(1:JZSOI,j2,j3)
        case ('CSILT')
          a_CSILT(1:JZSOI,j3)=data_3d(1:JZSOI,j2,j3)
        end select
      enddo
    enddo
  endif

  end subroutine ATS2EcoSIMData
!------------------------------------------------------------------------------------------

  subroutine EcoSIM2ATSData()
  !!grab data from ecosim and return it to ATS
  implicit none
  character(len=*), parameter :: subname=trim(mod_filename)//'::EcoSIM2ATSData'

  
  end subroutine EcoSIM2ATSData

!------------------------------------------------------------------------------------------

  subroutine Run_EcoSIM_one_step()
  !advance ecosim one time step
  implicit none
  character(len=*), parameter :: subname=trim(mod_filename)//'::Run_EcoSIM_one_step'


  !copy data from compuler to EcoSIM

  !run surface energy balance
  call SurfaceEBalance()

  !copy data back to coupler
  end subroutine Run_EcoSIM_one_step
!------------------------------------------------------------------------------------------

  subroutine Init_EcoSIM(jz,js,ncol)
  !initialize ecosim
  implicit none
  character(len=*), parameter :: subname=trim(mod_filename)//'::Init_EcoSIM'
  integer, intent(in) :: jz   !number of soil layers
  integer, intent(in) :: js   !number of snow layers
  integer, intent(in) :: ncol !number of column


  call InitSharedData(JZ,NCOL)

  call Init_EcoSIM_Soil()
  end subroutine Init_EcoSIM
!------------------------------------------------------------------------------------------

  subroutine SurfaceEBalance()

  implicit none

  end subroutine SurfaceEBalance
  


end module ATSCPLMod
