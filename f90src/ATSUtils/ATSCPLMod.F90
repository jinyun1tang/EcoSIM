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

  call c_f_pointer(props%shortwave_radiation%data, data, (/num_cols/))
  srad = data(:)

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

!!** Currently used by coupler **!!

  atm_n2 = props%atm_n2
  atm_o2 = props%atm_o2
  atm_co2 = props%atm_co2
  atm_ch4 = props%atm_ch4
  atm_n2o = props%atm_n2o
  atm_h2 = props%atm_h2
  atm_nh3 = props%atm_nh3

  call c_f_pointer(props%plant_wilting_factor%data, data2D, [(/size_col/),(/num_cols/)])
  a_WP=data2D(:,:)

  call c_f_pointer(props%rooting_depth_fraction%data, data2D, [(/size_col/),(/num_cols/)])
  a_FC=data2D(:,:)

  call c_f_pointer(state%bulk_density%data, data2D, [(/size_col/),(/num_cols/)])
  a_BKDSI=data2D(:,:)

  call c_f_pointer(state%subsurface_water_source%data, data2D, [(/size_col/),(/num_cols/)])
  a_SSWS=data2D(:,:)

  call c_f_pointer(state%subsurface_energy_source%data, data2D, [(/size_col/),(/num_cols/)])
  a_SSES=data2D(:,:)

  call c_f_pointer(state%surface_water_source%data, data, (/num_cols/))
  surf_w_source = data(:)

  call c_f_pointer(state%surface_energy_source%data, data, (/num_cols/))
  surf_e_source = data(:)

!!***********************************!!

  call c_f_pointer(state%porosity%data, data2D, [(/size_col/),(/num_cols/)])
  a_PORO=data2D(:,:)

  call c_f_pointer(state%water_content%data, data2D, [(/size_col/),(/num_cols/)])
  a_WC=data2D(:,:)

  call c_f_pointer(props%liquid_saturation%data, data2D, [(/size_col/),(/num_cols/)])
  a_LSAT=data2D(:,:)

  call c_f_pointer(props%relative_permeability%data, data2D, [(/size_col/),(/num_cols/)])
  a_RELPERM=data2D(:,:)

  call c_f_pointer(state%hydraulic_conductivity%data, data2D, [(/size_col/),(/num_cols/)])
  a_HCOND=data2D(:,:)

  call c_f_pointer(state%temperature%data, data2D, [(/size_col/),(/num_cols/)])
  a_TEMP=data2D(:,:)

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
  integer :: j1,j2,j3

  write(*,*) "Copying back"
  call SetBGCSizes(sizes)

  size_col = sizes%ncells_per_col_
  size_procs = state%porosity%cols

  !seems like we call the pointer as normal,
  !then just reverse the data

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

  call c_f_pointer(state%surface_water_source%data, data, (/num_cols/))
  data(:) = surf_w_source

  call c_f_pointer(state%surface_energy_source%data, data, (/num_cols/))
  data(:) = surf_e_source

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

  !For this we need to either specifically take energy away from the top layer
  !or add it to the surface_energy-source variable
  integer :: K, vec_size
  type (BGCSizes), intent(in) :: sizes
  integer :: size_col, num_cols

  size_col = sizes%ncells_per_col_
  num_cols = sizes%num_columns

  vec_size = size(surf_e_source)
  write(*,*) "in surface energy balance"

  call RunEcoSIMSurfaceBalance(size_col)

  !call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistanceLitRLay)
  !call RunSurfacePhysModel(M,NHE,NHW,NVS,NVN,ResistanceLitRLay,&
  !    KSatReductByRainKineticEnergyS,TopLayWatVol,HeatFlux2Ground,Qinfl2MicP)
  !do K=1,vec_size
  !  write(*,*) "surface energy (", K, ") = ", surf_e_source(K)
  !  write(*,*) "surface water (", K, ") = ", surf_w_source(K)
  !
  !  write(*,*) "changing surface vars"
  !  surf_e_source(K) = 8.08
  !  surf_w_source(K) = 8.08
  !end do

  write(*,*) "leaving surface energy balance"

  end subroutine SurfaceEBalance

!------------------------------------------------------------------------------------------

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
