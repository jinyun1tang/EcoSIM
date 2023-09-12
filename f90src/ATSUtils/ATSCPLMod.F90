module ATSCPLMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SharedDataMod
  use ATSEcoSIMInitMod
  use BGCContainers_module
  implicit none

  public
  character(len=*), private, parameter :: mod_filename=&
  __FILE__
!  integer :: JZSOI !number of soil layers
!  integer :: JSNO  !number of snow layers

! temporary data holder in ecosim
  real(r8), allocatable :: sw_rad(:)
  real(r8), allocatable :: lw_rad(:)
  real(r8), allocatable :: air_temp(:)
  real(r8), allocatable :: p_vap(:)
  real(r8), allocatable :: wind_speed(:)
  real(r8), allocatable :: precipitation_rain(:)

  !ATS variables
  real(r8), allocatable :: PORO(:,:) !porosity
  real(r8), allocatable :: L_DENS(:,:) !liquid density
  real(r8), allocatable :: WC(:,:) !Soil water content
  real(r8), allocatable :: WC_OLD(:,:) !saving the wc for testing
  real(r8), allocatable :: LSAT(:,:) !liquid saturation
  real(r8), allocatable :: RELPERM(:,:) !relative_permeability
  real(r8), allocatable :: HCOND(:,:) !hydraulic conductivity
  real(r8), allocatable :: TEMP(:,:) !temperature
  real(r8), allocatable :: FIELD_CAPACITY(:,:)
  real(r8), allocatable :: WILTING_POINT(:,:)

contains
!------------------------------------------------------------------------------------------

  subroutine ATS2EcoSIMData(ncol, state, props, sizes)
  !send data from ATS to ecosim
  implicit none

  ! BGC coupler variables
  type (BGCState), intent(in) :: state
  type (BGCProperties), intent(in) :: props
  type (BGCSizes), intent(out) :: sizes

  ! Ecosim variables
  real(r8), pointer :: data(:)
  real(r8), pointer :: data2D(:,:)
  integer :: ncol, nvar, size_col, num_cols
  integer :: j1,j2,j3

  write(*,*) "In the driver...."

  write(*,*) "Setting sizes"
  call SetBGCSizes(sizes)

  !ncol=size(filter_col)

  !1D vertical vector,
  write(*,*) "computing column size"

  size_col = sizes%ncells_per_col_
  num_cols = props%shortwave_radiation%size

  call c_f_pointer(props%shortwave_radiation%data, data, (/num_cols/))
  sw_rad = data(:)

  call c_f_pointer(props%longwave_radiation%data, data, (/num_cols/))
  lw_rad = data(:)

  call c_f_pointer(props%air_temperature%data, data, (/num_cols/))
  air_temp = data(:)

  call c_f_pointer(props%vapor_pressure_air%data, data, (/num_cols/))
  p_vap = data(:)

  call c_f_pointer(props%wind_speed%data, data, (/num_cols/))
  wind_speed = data(:)

  call c_f_pointer(props%precipitation%data, data, (/num_cols/))
  precipitation_rain = data(:)

  atm_n2 = props%atm_n2
  atm_o2 = props%atm_o2
  atm_co2 = props%atm_co2
  atm_ch4 = props%atm_ch4
  atm_n2o = props%atm_n2o
  atm_h2 = props%atm_h2
  atm_nh3 = props%atm_nh3

  call c_f_pointer(state%porosity%data, data2D, [(/size_col/),(/num_cols/)])
  PORO=data2D(:,:)

  call c_f_pointer(state%water_content%data, data2D, [(/size_col/),(/num_cols/)])
  WC=data2D(:,:)

  call c_f_pointer(props%liquid_saturation%data, data2D, [(/size_col/),(/num_cols/)])
  LSAT=data2D(:,:)

  call c_f_pointer(props%relative_permeability%data, data2D, [(/size_col/),(/num_cols/)])
  RELPERM=data2D(:,:)

  call c_f_pointer(state%hydraulic_conductivity%data, data2D, [(/size_col/),(/num_cols/)])
  HCOND=data2D(:,:)

  call c_f_pointer(state%temperature%data, data2D, [(/size_col/),(/num_cols/)])
  TEMP=data2D(:,:)

  write(*,*) "Data Transfer Finished"
  end subroutine ATS2EcoSIMData
!------------------------------------------------------------------------------------------

  subroutine EcoSIM2ATSData(ncol, state, sizes)
  !!grab data from ecosim and return it to ATS
  implicit none
  !character(len=*), parameter :: subname=trim(mod_filename)//'::EcoSIM2ATSData'

  type (BGCState), intent(in) :: state
  type (BGCSizes), intent(out) :: sizes

  ! Ecosim variables
  real(r8), pointer :: data(:)
  real(r8), pointer :: data2D(:,:)
  integer :: ncol, nvar, size_col, size_procs
  integer :: j1,j2,j3

  write(*,*) "Copying back"
  call SetBGCSizes(sizes)

  size_col = sizes%ncells_per_col_
  size_procs = state%porosity%cols

  write(*,*) "column size: ", size_col, " columns on this process: ", size_procs

  !seems like we call the pointer as normal,
  !then just reverse the data
  call c_f_pointer(state%water_content%data, data2D, [(/size_col/),(/size_procs/)])
  data2D(:,:)=WC

  call c_f_pointer(state%hydraulic_conductivity%data, data2D, [(/size_col/),(/size_procs/)])
  data2D(:,:)=HCOND

  call c_f_pointer(state%temperature%data, data2D, [(/size_col/),(/size_procs/)])
  data2D(:,:)=TEMP

  write(*,*) "finished copying back in driver"

  end subroutine EcoSIM2ATSData

!------------------------------------------------------------------------------------------

  subroutine Run_EcoSIM_one_step()
  !advance ecosim one time step
  implicit none
  !character(len=*), parameter :: subname=trim(mod_filename)//'::Run_EcoSIM_one_step'


  !copy data from compuler to EcoSIM

  !run surface energy balance
  call SurfaceEBalance()

  !copy data back to coupler
  end subroutine Run_EcoSIM_one_step
!------------------------------------------------------------------------------------------

  subroutine Init_EcoSIM(sizes)
  !initialize ecosim
  implicit none
  !character(len=*), parameter :: subname=trim(mod_filename)//'::Init_EcoSIM'
  type (BGCSizes), intent(in) :: sizes
  integer :: size_col, num_cols

  size_col = sizes%ncells_per_col_
  num_cols = sizes%num_columns

  call InitSharedData(size_col,num_cols)

  call Init_EcoSIM_Soil(size_col)
  end subroutine Init_EcoSIM
!------------------------------------------------------------------------------------------

  subroutine SurfaceEBalance()

  implicit none

  end subroutine SurfaceEBalance

  subroutine SetBGCSizes(sizes)

    use BGCContainers_module, only : BGCSizes

    implicit none

    type (BGCSizes), intent(out) :: sizes

    sizes%num_components = 1
    sizes%ncells_per_col_ = 100
    sizes%num_columns = 25

  end subroutine SetBGCSizes

!-----------------------------------------------------------------------------------------
end module ATSCPLMod
