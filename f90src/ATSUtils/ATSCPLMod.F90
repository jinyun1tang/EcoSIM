module ATSCPLMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SharedDataMod
  use ATSEcoSIMInitMod
  implicit none

  public
  character(len=*), private, parameter :: &
    mod_filename=__FILE__
  integer :: JZSOI   !number of soil layers
  integer :: JSNO    !number of snow layers

! temporary data holder in ecosim
  real(r8) :: atm_n2, atm_o2,atm_co2,atm_ch4,atm_n2o,atm_h2,atm_nh3
  real(r8) :: sw_rad, lw_rad, air_temp, p_vap, wind_speed, precipitation_rain
  real(r8), allocatable :: csand(:,:)
  real(r8), allocatable :: CSILT(:,:)
  real(r8), allocatable :: tairc(:)
  real(r8), allocatable :: uwind(:)
  real(r8), allocatable :: prec(:)
  real(r8), allocatable :: srad(:)
  real(r8), allocatable :: vpa(:)

  !ATS variables
  real(r8), allocatable :: PORO(:) !porosity
  real(r8), allocatable :: L_DENS(:) !liquid density
  real(r8), allocatable :: WC(:) !Soil water content
  real(r8), allocatable :: WC_OLD(:) !saving the wc for testing
  real(r8), allocatable :: L_SAT(:) !liquid saturation
  real(r8), allocatable :: REL_PERM(:) !relative_permeability
  real(r8), allocatable :: H_COND(:) !hydraulic conductivity
  real(r8), allocatable :: TEMP(:) !temperature

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
  integer :: ncol, nvar, size_col, size_procs
  integer :: j1,j2,j3

  write(*,*) "In the driver...."


  write(*,*) "Setting sizes"
  call SetBGCSizes(sizes)
  !ncol=size(filter_col)

  !1D vertical vector,
  !variables that take on a different value in each cell
  !Bulk of data will go here
  !nvar=size(var_2d)

  !do j1=1,nvar
  !case ('CSAND')  !g/kg soil
  !  csand(1:JZSOI,ncol)=data_3d(1:JZSOI,j1)
  !case ('CSILT')
  !  CSILT(1:JZSOI,ncol)=data_3d(1:JZSOI,j1)
  !Variables related to flow:
  write(*,*) "computing column size"

  size_col = sizes%ncells_per_col_
  size_procs = props%shortwave_radiation%size

  write(*,*) "Column size is: ", size_col
  write(*,*) "surface properties"

  call c_f_pointer(props%shortwave_radiation%data, data, (/size_procs/))
  sw_rad = data(:)

  call c_f_pointer(props%longwave_radiation%data, data, (/size_procs/))
  lw_rad = data(:)

  call c_f_pointer(props%air_temperature%data, data, (/size_procs/))
  air_temp = data(:)

  call c_f_pointer(props%vapor_pressure_air%data, data, (/size_procs/))
  p_vap = data(:)

  call c_f_pointer(props%wind_speed%data, data, (/size_procs/))
  wind_speed = data(:)

  call c_f_pointer(props%precipitation%data, data, (/size_procs/))
  precipitation_rain = data(:)

  write(*,*) "writing atm abundances"
  atm_n2 = props%atm_n2
  atm_o2 = props%atm_o2
  atm_co2 = props%atm_co2
  atm_ch4 = props%atm_ch4
  atm_n2o = props%atm_n2o
  atm_h2 = props%atm_h2
  atm_nh3 = props%atm_nh3

  write(*,*) "looping over datasets starting with porosity"
  call c_f_pointer(state%porosity%data, data, (/size_col/))
  PORO = data(:)

  write(*,*) "finished copying poro"
  !do j3 = 1, size_col
  !  PORO(1:JZSOI)=data(j3)
  !enddo

  write(*,*) "Porosity finished, continuing"
  call c_f_pointer(state%liquid_density%data, data, [(/size_col/),(/size_procs/)])
  data(:,:)=L_DENS

  call c_f_pointer(state%water_content%data, data, [(/size_col/),(/size_procs/)])
  data(:,:)=WC

  call c_f_pointer(props%liquid_saturation%data, data, [(/size_col/),(/size_procs/)])
  data(:,:)=L_SAT

  call c_f_pointer(props%relative_permeability%data, data, [(/size_col/),(/size_procs/)])
  data(:,:)=REL_PERM

  call c_f_pointer(state%hydraulic_conductivity%data, data, [(/size_col/),(/size_procs/)])
  data(:,:)=H_COND

  call c_f_pointer(state%temperature%data, data, [(/size_col/),(/size_procs/)])
  data(:,:)=TEMP

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

  WC_OLD = WC

  do j3 = 1, size_col
    WC(j3) = 2.0*WC(j3)
    !write(*,*) "Old value: ", WC_OLD(j3), " New value: ", WC(j3)
  enddo

  !seems like we call the pointer as normal,
  !then just reverse the data
  call c_f_pointer(state%liquid_density%data, data2D, [(/size_col/),(/size_procs/)])
  data(:,:)=L_DENS

  call c_f_pointer(state%water_content%data, data2D, [(/size_col/),(/size_procs/)])
  data(:,:)=WC

  call c_f_pointer(state%hydraulic_conductivity%data, data2D, [(/size_col/),(/size_procs/)])
  data(:,:)=H_COND

  call c_f_pointer(state%temperature%data, data2D, [(/size_col/),(/size_procs/)])
  data(:,:)=TEMP

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

  subroutine Init_EcoSIM(jz,js,ncol)
  !initialize ecosim
  implicit none
  !character(len=*), parameter :: subname=trim(mod_filename)//'::Init_EcoSIM'
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

  subroutine SetBGCSizes(sizes)

    use BGCContainers_module, only : BGCSizes

    implicit none

    type (BGCSizes), intent(out) :: sizes

    sizes%num_components = 1
    sizes%ncells_per_col_ = 100

  end subroutine SetBGCSizes

!-----------------------------------------------------------------------------------------
end module ATSCPLMod
