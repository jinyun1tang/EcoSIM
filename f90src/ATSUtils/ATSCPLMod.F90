module ATSCPLMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SharedDataMod
  use ATSEcoSIMInitMod
  use ATSEcoSIMAdvanceMod
  use BGCContainers_module
  implicit none

  public
  character(len=*), private, parameter :: mod_filename=&
  __FILE__

contains
!------------------------------------------------------------------------------------------

  ! Function to check for NaN in an array
  !function is_nan(x) result(mask)
  !  real(r8), intent(in) :: x(:)
  !  logical, dimension(size(x)) :: mask
  !  integer :: i

    !allocate(mask(size(x)))
  !  do i = 1, size(x)
  !    mask(i) = (x(i) /= x(i))  ! NaN is the only value that is not equal to itself
  !  end do
  !end function is_nan


  subroutine ATS2EcoSIMData(ncol, state, props, sizes)
  implicit none

  ! BGC coupler variables
  type (BGCState), intent(in) :: state
  type (BGCProperties), intent(in) :: props
  type (BGCSizes), intent(out) :: sizes

  ! Ecosim variables
  real(r8), pointer :: data(:), nested_ptr(:)
  real(r8), pointer :: data2D(:,:)
  real(r8), pointer :: temp, temp_deref, double_deref
  real(r8) :: real_deref
  type(c_ptr), pointer :: cptr_temp, ptr1
  real(r8), target :: target_val
  real(r8), pointer :: ptr(:,:)
  real(r8), pointer :: temp_array(:,:)
  integer :: ncol, nvar, size_col, num_cols, size_col_pad
  integer :: j1,j2,j3,i,j
  integer :: test_rows, test_columns
  real(r8) :: temp_eq, double_eq
  type(c_ptr) :: data_ptr

  !call SetBGCSizes(sizes)

  size_col = sizes%ncells_per_col_
  num_cols = props%shortwave_radiation%size

  size_col_pad = size_col+30

  allocate(temp_array(size_col, num_cols))

  data_ptr = state%temperature%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_TEMP=data2D(:,:)

  data_ptr = props%depth%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_CumDepz2LayBottom_vr = data2D(:,:)

  data_ptr = props%dz%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_dz = data2D(:,:)

  data_ptr = props%volume%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_Volume = data2D(:,:)

  data_ptr = props%volume%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_AREA3 = data2D(:,:)

  data_ptr = state%water_content%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_WC = data2D(:,:)

  data_ptr = props%volume%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_AreaZ = data2D(:,:)

  do i = 1, size_col
     a_AreaZ(i,1) = a_Volume(i,1)/a_dz(i,1)
  end do

  data_ptr = state%temperature%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_TEMP = data2D(:,:)

  data_ptr = state%bulk_density%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_BKDSI = data2D(:,:)

  data_ptr = state%liquid_density%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_LDENS = data2D(:,:)  

  data_ptr = state%matric_pressure%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_MATP = data2D(:,:)

  data_ptr = state%porosity%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_PORO = data2D(:,:)

  data_ptr = props%liquid_saturation%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_LSAT = data2D(:,:)

  data_ptr = props%relative_permeability%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_RELPERM = data2D(:,:)

  data_ptr = state%hydraulic_conductivity%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_HCOND = data2D(:,:)

  data_ptr = props%rooting_depth_fraction%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_FC = data2D(:,:)

  data_ptr = state%subsurface_water_source%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_SSWS = data2D(:,:)

  data_ptr = state%subsurface_energy_source%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_SSES = data2D(:,:)

  call c_f_pointer(props%shortwave_radiation%data, data, (/num_cols/))
  swrad = data(:)

  call c_f_pointer(props%longwave_radiation%data, data, (/num_cols/))
  sunrad = data(:)

  call c_f_pointer(props%air_temperature%data, data, (/num_cols/))
  tairc = data(:)

  call c_f_pointer(props%vapor_pressure_air%data, data, (/num_cols/))
  vpair = data(:)

  call c_f_pointer(props%wind_speed%data, data, (/num_cols/))
  uwind = data(:)

  !call c_f_pointer(props%precipitation%data, data, (/num_cols/))
  !p_rain = data(:)

  !call c_f_pointer(props%precipitation_snow%data, data, (/num_cols/))
  !p_snow = data(:)

  call c_f_pointer(props%aspect%data, data, (/num_cols/))
  a_ASP = data(:)

  call c_f_pointer(props%LAI%data, data, (/num_cols/))
  a_LAI = data(:)  

  call c_f_pointer(props%SAI%data, data, (/num_cols/))
  a_SAI = data(:)

  call c_f_pointer(props%vegetation_type%data, data, (/num_cols/))
  a_VEG = data(:)

  !a_MATP(:,:) = -6.9
  !do i = 1, size_col
  !  a_MATP(i, 1) = 100.0
  !end do

  atm_n2 = props%atm_n2
  atm_o2 = props%atm_o2
  atm_co2 = props%atm_co2
  atm_ch4 = props%atm_ch4
  atm_n2o = props%atm_n2o
  atm_h2 = props%atm_h2
  atm_nh3 = props%atm_nh3
  heat_capacity = props%heat_capacity
  pressure_at_field_capacity = props%field_capacity
  pressure_at_wilting_point = props%wilting_point
  p_bool = props%p_bool

  if(p_bool)THEN
    call c_f_pointer(props%precipitation%data, data, (/num_cols/))
    p_total = data(:)
  else
    call c_f_pointer(props%precipitation%data, data, (/num_cols/))
    p_rain = data(:)

    call c_f_pointer(props%precipitation_snow%data, data, (/num_cols/))
    p_snow = data(:)
  endif

  call c_f_pointer(state%surface_water_source%data, data, (/num_cols/))
  surf_w_source = data(:)

  call c_f_pointer(state%surface_energy_source%data, data, (/num_cols/))
  surf_e_source = data(:)

  call c_f_pointer(state%snow_depth%data, data, (/num_cols/))
  surf_snow_depth = data(:)

  end subroutine ATS2EcoSIMData
!------------------------------------------------------------------------------------------

  subroutine EcoSIM2ATSData(ncol, state, sizes)
  implicit none
  type (BGCState), intent(in) :: state
  type (BGCSizes), intent(out) :: sizes

  real(r8), pointer :: data(:)
  real(r8), pointer :: data2D(:,:)
  integer :: num_cols, ncol, nvar, size_col, size_procs
  integer :: j1,j2,j3, i

  !call SetBGCSizes(sizes)

  size_col = sizes%ncells_per_col_
  size_procs = state%porosity%cols

  call c_f_pointer(state%subsurface_water_source%data, data2D, [(/size_col/),(/num_cols/)])
  data2D(:,:)=a_SSWS

  call c_f_pointer(state%subsurface_energy_source%data, data2D, [(/size_col/),(/num_cols/)])
  data2D(:,:)=a_SSES

  call c_f_pointer(state%surface_water_source%data, data, (/num_cols/))
  data(:) = surf_w_source

  call c_f_pointer(state%surface_energy_source%data, data, (/num_cols/))
  data(:) = surf_e_source

  !write(*,*) "surf_e_source (ATSCPL): ", surf_e_source

  call c_f_pointer(state%snow_depth%data, data, (/num_cols/))
  data(:) = surf_snow_depth

  end subroutine EcoSIM2ATSData

!------------------------------------------------------------------------------------------

  subroutine Run_EcoSIM_one_step(sizes)
  implicit none

  type (BGCSizes), intent(in) :: sizes

  !copy data from compuler to EcoSIM

  !run surface energy balance
  call SurfaceEBalance(sizes)

  !copy data back to coupler
  end subroutine Run_EcoSIM_one_step
!------------------------------------------------------------------------------------------

  subroutine Init_EcoSIM(sizes)
  !initialize ecosim
  implicit none

  type (BGCSizes), intent(in) :: sizes
  integer :: size_col, num_cols

  size_col = sizes%ncells_per_col_
  num_cols = sizes%num_columns

  call InitSharedData(size_col,num_cols)

  call Init_EcoSIM_Soil(num_cols)
  end subroutine Init_EcoSIM
!------------------------------------------------------------------------------------------

  subroutine SurfaceEBalance(sizes)
  implicit none

  integer :: K, vec_size
  type (BGCSizes), intent(in) :: sizes
  integer :: size_col, num_cols

  size_col = sizes%ncells_per_col_
  num_cols = sizes%num_columns

  !vec_size = size(surf_e_source)
  call RunEcoSIMSurfaceBalance(num_cols)

  end subroutine SurfaceEBalance

!------------------------------------------------------------------------------------------

  subroutine SetBGCSizes(sizes)

    use BGCContainers_module, only : BGCSizes

    implicit none

    type (BGCSizes), intent(out) :: sizes

    sizes%num_components = 1
    sizes%ncells_per_col_ = 100
    sizes%num_columns = 1

    write(*,*) "(SetBGCSizes f): N_cells, N_cols ", sizes%ncells_per_col_, sizes%num_columns
  end subroutine SetBGCSizes

!-----------------------------------------------------------------------------------------
end module ATSCPLMod
