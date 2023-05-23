module ATSCPLMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SharedDataMod
  use ATSEcoSIMInitMod
  implicit none

  public
<<<<<<< HEAD
  character(len=*), private, parameter :: &
    mod_filename=__FILE__
=======
  character(len=*), private, parameter :: mod_filename=__FILE__
>>>>>>> d038c8f (Initial reworking, addition of water variables, and alquimia data funcs)
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
<<<<<<< HEAD
  real(r8), allocatable :: PORO(:) !porosity
  real(r8), allocatable :: L_DENS(:,:) !liquid density
  real(r8), allocatable :: WC(:,:) !Soil water content
  real(r8), allocatable :: WC_OLD(:,:) !saving the wc for testing
=======
  real(r8), allocatable :: PORO(:,:) !porosity
  real(r8), allocatable :: L_DENS(:,:) !liquid density
  real(r8), allocatable :: WC(:,:) !Soil water content
>>>>>>> d038c8f (Initial reworking, addition of water variables, and alquimia data funcs)
  real(r8), allocatable :: L_SAT(:,:) !liquid saturation
  real(r8), allocatable :: REL_PERM(:,:) !relative_permeability
  real(r8), allocatable :: H_COND(:,:) !hydraulic conductivity
  real(r8), allocatable :: TEMP(:,:) !temperature

contains
!------------------------------------------------------------------------------------------

<<<<<<< HEAD
  subroutine ATS2EcoSIMData(ncol, state, props, sizes)
=======
  subroutine ATS2EcoSIMData(ncol,state, aux_data, prop)
>>>>>>> d038c8f (Initial reworking, addition of water variables, and alquimia data funcs)
  !send data from ATS to ecosim
  implicit none

  ! BGC coupler variables
  type (BGCState), intent(in) :: state
<<<<<<< HEAD
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
  call c_f_pointer(state%porosity%data, data2D, [(/size_col/),(/size_procs/)])
  PORO=data(:,:)

  write(*,*) "finished copying poro"
  !do j3 = 1, size_col
  !  PORO(1:JZSOI)=data(j3)
  !enddo

  write(*,*) "Porosity finished, continuing"
  call c_f_pointer(state%liquid_density%data, data2D, [(/size_col/),(/size_procs/)])
  L_DENS=data2D(:,:)

  call c_f_pointer(state%water_content%data, data2D, [(/size_col/),(/size_procs/)])
  WC=data2D(:,:)
=======
  type (BGCAuxiliaryData), intent(in) :: aux_data
  type (BGCProperties), intent(in) :: prop


  ! Ecosim variables
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

  !ncol=size(filter_col)

  if (present(data_1d) .and. present(var_1d) .and. ncol .EQ. 0)then
  !domain specific scalar
  !Variables with only a single value over the domain
  !loops over j1 (number of vars)
  !Only do this on the first column (if ncol=0)
  !I think we basically have to hardcode this with foreknowledge of what is
  !going to be in the Alquimia-like dictionary
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
  !Variables that are single valued along a column
  !j1 - number of variables
  ! changed the loop over j2 because we call this column by column
    nvar=size(var_2d)
    do j1=1,nvar
      select case(var_2d(j1))
      case ('TAIRC')     !air temperature, oC
        tairc(ncol)=data_2d(j1,ncol)
      case ('PREC')      !precipitation, mm H2O/hr
        prec(ncol)=data_2d(j1,ncol)
      case ('WINDH')     !horizontal wind speed,   m/s
        uwind(ncol)=data_2d(j1,ncol)
      case ('DWPTH')     !atmospheric vapor pressure, kPa
        vpa(ncol)=data_2d(j1,ncol)
      case ('SRADH')     !Incident solar radiation, W/m2
        srad(ncol)=data_2d(j1,ncol)
      end select
    enddo
  endif

  if (present(data_3d) .and. present(var_3d))then
  !1D vertical vector,
  !variables that take on a different value in each cell
  !Bulk of data will go here
    nvar=size(var_2d)
    do j1=1,nvar
      select case (var_2d(j1))
      case ('CSAND')  !g/kg soil
        csand(1:JZSOI,ncol)=data_3d(1:JZSOI,j1)
      case ('CSILT')
        CSILT(1:JZSOI,ncol)=data_3d(1:JZSOI,j1)
      !Variables related to flow:
      case ('PORO')
        call c_f_pointer(state%porosity%data, data, (/size_col/))
        do j3 = 1, size_col
          PORO(1:JZSOI,ncol)=data
        enddo
      case('L_DENS')
        call c_f_pointer(state%liquid_density%data, data, (/size_col/))
        do j3 = 1,size_col
          L_DENS(1:JZSOI,ncol)=data
        enddo
      case('WC')
        call c_f_pointer(state%water_content%data, data, (/size_col/))
        do j3 = 1,size_col
          WC(1:JZSOI,ncol)=data
        enddo
      case('L_SAT')
        call c_f_pointer(props%liquid_saturation%data, data, (/size_col/))
        do j3 = 1,size_col
          L_SAT(1:JZSOI,ncol)=data
        enddo
      case('REL_PERM')
        call c_f_pointer(props%relative_permeability%data, data, (/size_col/))
        do j3 = 1,size_col
          REL_PERM(1:JZSOI,ncol)=data
        enddo
      case('H_COND')
        call c_f_pointer(state%hydraulic_conductivity%data, data, (/size_col/))
        do j3 = 1,size_col
          H_COND(1:JZSOI,ncol)=data
        enddo
      case('TEMP')
        call c_f_pointer(state%temperature%data, data, (/size_col/))
        do j3 = 1,size_col
          TEMP(1:JZSOI,ncol)=data
        enddo

      end select
    enddo
  endif
>>>>>>> d038c8f (Initial reworking, addition of water variables, and alquimia data funcs)

  call c_f_pointer(props%liquid_saturation%data, data2D, [(/size_col/),(/size_procs/)])
  L_SAT=data2D(:,:)

  call c_f_pointer(props%relative_permeability%data, data2D, [(/size_col/),(/size_procs/)])
  REL_PERM=data2D(:,:)

  call c_f_pointer(state%hydraulic_conductivity%data, data2D, [(/size_col/),(/size_procs/)])
  H_COND=data2D(:,:)

  call c_f_pointer(state%temperature%data, data2D, [(/size_col/),(/size_procs/)])
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

  WC_OLD = WC

<<<<<<< HEAD
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
=======

>>>>>>> d038c8f (Initial reworking, addition of water variables, and alquimia data funcs)
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
