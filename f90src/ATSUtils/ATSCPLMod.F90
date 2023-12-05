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

  call SetBGCSizes(sizes)

  !ncol=size(filter_col)

  !1D vertical vector

  size_col = sizes%ncells_per_col_
  num_cols = props%shortwave_radiation%size

  !Pass the data over:
  !Variables that need to be passed still:
  !PSISM1, soil matric pressure, MPa
  !X - VLWatMicP1, top layer volumetric liquid water content, m3/d2, d2 means grid area,
  !VLiceMicP1,  top layer volumetric ice content, (better pass in ice mass and then convert)
  !SoilMicPMassLayer, total topsoil layer mass, 10^6 g (or Mg),
  !X - TKSoi1, topsoil temperature, K
  !X - SoiBulkDensity, topsoil layer bulk density, Mg m-3
  !VLHeatCapacity, topsoil layer heat capacity, MJ m-3 K-1
  !X - SoilFracAsMicP, Fraction of soil as micropore, this will be set to 1

  call c_f_pointer(state%water_content%data, data2D, [(/size_col/),(/num_cols/)])
  a_WC=data2D(:,:)

  call c_f_pointer(state%temperature%data, data2D, [(/size_col/),(/num_cols/)])
  a_TEMP=data2D(:,:)

  call c_f_pointer(state%bulk_density%data, data2D, [(/size_col/),(/num_cols/)])
  a_BKDSI=data2D(:,:)

  call c_f_pointer(state%matric_pressure%data, data2D, [(/size_col/),(/num_cols/)])
  a_MATP=data2D(:,:)

  call c_f_pointer(state%porosity%data, data2D, [(/size_col/),(/num_cols/)])
  a_PORO=data2D(:,:)

  call c_f_pointer(props%liquid_saturation%data, data2D, [(/size_col/),(/num_cols/)])
  a_LSAT=data2D(:,:)

  call c_f_pointer(props%relative_permeability%data, data2D, [(/size_col/),(/num_cols/)])
  a_RELPERM=data2D(:,:)

  call c_f_pointer(state%hydraulic_conductivity%data, data2D, [(/size_col/),(/num_cols/)])
  a_HCOND=data2D(:,:)

  call c_f_pointer(props%depth_c%data, data2D, [(/size_col/),(/num_cols/)])
  a_CumDepth2LayerBottom=data2D(:,:)

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
  prec = data(:)

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

  ! Move write statements to the end
  !write(*, *) "Shortwave Radiation (SRAD):", srad
  !write(*, *) "Longwave Radiation (SUNRAD):", sunrad
  !write(*, *) "Air Temperature (TAIRC):", tairc
  !write(*, *) "Vapor Pressure in Air (VPAIR):", vpair
  !write(*, *) "Wind Speed (UWIND):", uwind
  !write(*, *) "Precipitation (PREC):", prec
  !write(*, *) "Atmospheric N2 (ATM_N2):", atm_n2
  !write(*, *) "Atmospheric O2 (ATM_O2):", atm_o2
  !write(*, *) "Atmospheric CO2 (ATM_CO2):", atm_co2
  !write(*, *) "Atmospheric CH4 (ATM_CH4):", atm_ch4
  !write(*, *) "Atmospheric N2O (ATM_N2O):", atm_n2o
  !write(*, *) "Atmospheric H2 (ATM_H2):", atm_h2
  !write(*, *) "Atmospheric NH3 (ATM_NH3):", atm_nh3
  !write(*, *) "Heat Capacity (HEAT_CAPACITY):", heat_capacity
  !write(*, *) "Pressure at Field Capacity (PRESSURE_AT_FIELD_CAPACITY):", pressure_at_field_capacity
  !write(*, *) "Pressure at Wilting Point (PRESSURE_AT_WILTING_POINT):", pressure_at_wilting_point

  !call c_f_pointer(props%plant_wilting_factor%data, data2D, [(/size_col/),(/num_cols/)])
  !a_WP=data2D(:,:)

  call c_f_pointer(props%rooting_depth_fraction%data, data2D, [(/size_col/),(/num_cols/)])
  a_FC=data2D(:,:)

  call c_f_pointer(state%bulk_density%data, data2D, [(/size_col/),(/num_cols/)])
  a_BKDSI=data2D(:,:)

  call c_f_pointer(state%subsurface_water_source%data, data2D, [(/size_col/),(/num_cols/)])
  a_SSWS=data2D(:,:)

  call c_f_pointer(state%subsurface_energy_source%data, data2D, [(/size_col/),(/num_cols/)])
  a_SSES=data2D(:,:)

  !call c_f_pointer(state%surface_water_source%data, data, (/num_cols/))
  !surf_w_source = data(:)

  call c_f_pointer(state%surface_energy_source%data, data, (/num_cols/))
  surf_e_source = data(:)

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
  integer :: num_cols, ncol, nvar, size_col, size_procs
  integer :: j1,j2,j3, i

  write(*,*) "Copying back"
  call SetBGCSizes(sizes)

  size_col = sizes%ncells_per_col_
  size_procs = state%porosity%cols

  !Pass the data over:
  !Variables that need to be passed still:
  !PSISM1, soil matric pressure, MPa
  !VLWatMicP1, top layer volumetric liquid water content, m3/d2, d2 means grid area,
  !VLiceMicP1,  top layer volumetric ice content, (better pass in ice mass and then convert)
  !SoilMicPMassLayer, total topsoil layer mass, 10^6 g (or Mg),
  !TKSoi1, topsoil temperature, K
  !SoiBulkDensity, topsoil layer bulk density, Mg m-3
  !VLHeatCapacity, topsoil layer heat capacity, MJ m-3 K-1
  !SoilFracAsMicP, Fraction of soil as micropore, this will be set to 1


  call c_f_pointer(state%bulk_density%data, data2D, [(/size_col/),(/size_procs/)])
  data2D(:,:)=a_BKDSI

  call c_f_pointer(state%water_content%data, data2D, [(/size_col/),(/size_procs/)])
  data2D(:,:)=a_WC

  call c_f_pointer(state%hydraulic_conductivity%data, data2D, [(/size_col/),(/size_procs/)])
  data2D(:,:)=a_HCOND

  call c_f_pointer(state%temperature%data, data2D, [(/size_col/),(/size_procs/)])
  data2D(:,:)=a_TEMP

  call c_f_pointer(state%subsurface_water_source%data, data2D, [(/size_col/),(/size_procs/)])
  data2D(:,:)=a_SSWS

  call c_f_pointer(state%subsurface_energy_source%data, data2D, [(/size_col/),(/size_procs/)])
  data2D(:,:)=a_SSES

  !call c_f_pointer(state%surface_water_source%data, data, (/num_cols/))
  !data(:) = surf_w_source

  call c_f_pointer(state%surface_energy_source%data, data, (/num_cols/))
  data(:) = surf_e_source

  write(*,*) "finished copying back in driver"

  end subroutine EcoSIM2ATSData

!------------------------------------------------------------------------------------------

  subroutine Run_EcoSIM_one_step(sizes)
  !advance ecosim one time step
  implicit none
  !character(len=*), parameter :: subname=trim(mod_filename)//'::Run_EcoSIM_one_step'
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
  !character(len=*), parameter :: subname=trim(mod_filename)//'::Init_EcoSIM'
  type (BGCSizes), intent(in) :: sizes
  integer :: size_col, num_cols

  size_col = sizes%ncells_per_col_
  num_cols = sizes%num_columns

  call InitSharedData(size_col,num_cols)

  call Init_EcoSIM_Soil(size_col)
  end subroutine Init_EcoSIM
!------------------------------------------------------------------------------------------

  subroutine SurfaceEBalance(sizes)
  implicit none

  !For this we need to either specifically take energy away from the top layer
  !or add it to the surface_energy-source variable
  integer :: K, vec_size
  type (BGCSizes), intent(in) :: sizes
  integer :: size_col, num_cols

  size_col = sizes%ncells_per_col_
  num_cols = sizes%num_columns

  vec_size = size(surf_e_source)
  write(*,*) "Surface energy before running RunEcoSIMSurfaceBalance"

  do K=1,vec_size
    write(*,*) "surface energy (", K, ") = ", surf_e_source(K)
  !  write(*,*) "surface water (", K, ") = ", surf_w_source(K)
  end do

  !needs number of columns
  call RunEcoSIMSurfaceBalance(num_cols)

  write(*,*) "Surface energy After running RunEcoSIMSurfaceBalance"

  do K=1,vec_size
    write(*,*) "surface energy (", K, ") = ", surf_e_source(K)
  !  write(*,*) "surface water (", K, ") = ", surf_w_source(K)
  end do

  write(*,*) "leaving surface energy balance"

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
