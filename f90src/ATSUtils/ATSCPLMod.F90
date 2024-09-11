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
  integer :: ncol, nvar, size_col, num_cols
  integer :: j1,j2,j3,i,j
  real(r8) :: temp_eq, double_eq

  call SetBGCSizes(sizes)

  size_col = sizes%ncells_per_col_
  num_cols = props%shortwave_radiation%size

  allocate(temp_array(size_col, num_cols))

  write(*,*) "size_col: ", size_col
  write(*,*) "num_cols: ", num_cols
  
  call c_f_pointer(props%depth%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_CumDepz2LayerBot_vr=data2D(:,:)

  call c_f_pointer(props%dz%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_dz=data2D(:,:)

  call c_f_pointer(props%volume%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_Volume=data2D(:,:)

  call c_f_pointer(props%volume%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_AREA3=data2D(:,:) 

  call c_f_pointer(state%water_content%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_WC=data2D(:,:)

  call c_f_pointer(props%volume%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_AreaZ=data2D(:,:)

  do i = 1, size_col
     a_AreaZ(i,1) = a_Volume(i,1)/a_dz(i,1)
  end do

  call c_f_pointer(state%temperature%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_TEMP=data2D(:,:)

  call c_f_pointer(state%bulk_density%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_BKDSI=data2D(:,:)

  call c_f_pointer(state%matric_pressure%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_MATP=data2D(:,:)

  do i = 1, size_col
     a_MATP(i, 1) = 100.0 
  end do

  call c_f_pointer(state%porosity%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_PORO=data2D(:,:)

  call c_f_pointer(props%liquid_saturation%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_LSAT=data2D(:,:)

  call c_f_pointer(props%relative_permeability%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_RELPERM=data2D(:,:)

  call c_f_pointer(state%hydraulic_conductivity%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_HCOND=data2D(:,:)

  call c_f_pointer(props%rooting_depth_fraction%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_FC=data2D(:,:)

  call c_f_pointer(state%subsurface_water_source%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_SSWS=data2D(:,:)

  call c_f_pointer(state%subsurface_energy_source%data, cptr_temp)
  call c_f_pointer(cptr_temp,data2D,[(/size_col/),(/num_cols/)])
  a_SSES=data2D(:,:)

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

  call c_f_pointer(props%precipitation%data, data, (/num_cols/))
  p_rain = data(:)

  !call c_f_pointer(props%precipitation_snow%data, data, (/num_cols/))
  !p_snow = data(:)

  call c_f_pointer(props%aspect%data, data, (/num_cols/))
  a_ASP = data(:)

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

  call c_f_pointer(state%surface_water_source%data, data, (/num_cols/))
  surf_w_source = data(:)

  call c_f_pointer(state%surface_energy_source%data, data, (/num_cols/))
  surf_e_source = data(:)
  
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

  call SetBGCSizes(sizes)

  size_col = sizes%ncells_per_col_
  size_procs = state%porosity%cols

  !call c_f_pointer(state%bulk_density%data, data2D, [(/size_col/),(/size_procs/)])
  !data2D(:,:)=a_BKDSI

  !call c_f_pointer(state%water_content%data, data2D, [(/size_col/),(/size_procs/)])
  !data2D(:,:)=a_WC

  !call c_f_pointer(state%hydraulic_conductivity%data, data2D, [(/size_col/),(/size_procs/)])
  !data2D(:,:)=a_HCOND

  !call c_f_pointer(state%temperature%data, data2D, [(/size_col/),(/size_procs/)])
  !data2D(:,:)=a_TEMP

  !call c_f_pointer(state%subsurface_water_source%data, data2D, [(/size_col/),(/size_procs/)])
  !data2D(:,:)=a_SSWS

  !call c_f_pointer(state%subsurface_energy_source%data, data2D, [(/size_col/),(/size_procs/)])
  !data2D(:,:)=a_SSES

  call c_f_pointer(state%surface_water_source%data, data, (/num_cols/))
  data(:) = surf_w_source

  call c_f_pointer(state%surface_energy_source%data, data, (/num_cols/))
  data(:) = surf_e_source

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

  vec_size = size(surf_e_source)
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

  end subroutine SetBGCSizes

!-----------------------------------------------------------------------------------------
end module ATSCPLMod
