program EcoATSTest
  use ATSCPLMod
  implicit none

  ! Declare variables
  type (BGCState) :: state
  type (BGCProperties) :: props
  type (BGCSizes) :: sizes

  ! Initialize variables with simple values (sample values)
  integer, parameter :: ncol = 5    ! Number of columns
  integer, parameter :: size_col = 100 ! Number of cells per column
  integer, parameter :: num_cols = ncol / size_col ! Number of columns

  ! Allocate and initialize arrays with sample values
  allocate(state%plant_wilting_factor%data(size_col, num_cols))
  allocate(state%rooting_depth_fraction%data(size_col, num_cols))
  allocate(state%bulk_density%data(size_col, num_cols))
  allocate(props%shortwave_radiation%data(num_cols))
  allocate(props%longwave_radiation%data(num_cols))
  allocate(props%air_temperature%data(num_cols))
  allocate(props%vapor_pressure_air%data(num_cols))
  allocate(props%wind_speed%data(num_cols))
  allocate(props%precipitation%data(num_cols))

  ! Set sample values for the arrays
  state%plant_wilting_factor%data = 0.5_r8  ! Sample values for plant_wilting_factor
  state%rooting_depth_fraction%data = 0.2_r8 ! Sample values for rooting_depth_fraction
  state%bulk_density%data = 1.2_r8            ! Sample values for bulk_density
  props%shortwave_radiation%data = 300.0_r8    ! Sample values for shortwave_radiation
  props%longwave_radiation%data = 100.0_r8     ! Sample values for longwave_radiation
  props%air_temperature%data = 25.0_r8         ! Sample values for air_temperature
  props%vapor_pressure_air%data = 18.0_r8      ! Sample values for vapor_pressure_air
  props%wind_speed%data = 5.0_r8              ! Sample values for wind_speed
  props%precipitation%data = 0.0_r8           ! Sample values for precipitation

  ! Initialize sizes structure
  sizes%num_components = 1
  sizes%ncells_per_col_ = size_col
  sizes%num_columns = num_cols

  ! Call initialization subroutine from ATSCPLMod module
  call Init_EcoSIM(sizes)

  ! Call ATS2EcoSIMData to send data from ATS to EcoSIM
  call ATS2EcoSIMData(ncol, state, props, sizes)

  ! Call Run_EcoSIM_one_step to advance EcoSIM
  call Run_EcoSIM_one_step()

  ! Call EcoSIM2ATSData to get data back from EcoSIM to ATS
  call EcoSIM2ATSData(ncol, state, sizes)

  ! Clean up (deallocate arrays)
  deallocate(state%plant_wilting_factor%data)
  deallocate(state%rooting_depth_fraction%data)
  deallocate(state%bulk_density%data)
  deallocate(props%shortwave_radiation%data)
  deallocate(props%longwave_radiation%data)
  deallocate(props%air_temperature%data)
  deallocate(props%vapor_pressure_air%data)
  deallocate(props%wind_speed%data)
  deallocate(props%precipitation%data)

end program EcoATSTest
