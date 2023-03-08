module EcoSIMCtrlDataType

! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__
  public
  real(r8) :: ZERO
  real(r8) :: ZERO2
  real(r8),target,allocatable :: ZEROS(:,:)
  real(r8),target,allocatable :: ZEROP(:,:,:)
  real(r8),target,allocatable :: ZEROQ(:,:,:)
  real(r8),target,allocatable :: ZEROL(:,:,:)
  real(r8),target,allocatable :: ZEROS2(:,:)

  integer :: NPX             !number of cycles per hour for water, heat, and solute flux calcns
  integer :: NPY             !number of cycles per NPX for gas flux calculations
  integer :: IGO             !flag for first scenario
  integer :: IDAYR           !day of recovery from earlier run
  integer :: IYRC            !current year
  integer :: IYRR            !year of recovery from earlier run
  integer :: NYR             !flag for new year
  integer :: ITERM           !end date for reading weather data
  integer :: IFIN            !end date for reading weather data
  integer :: JOUT            !frequency of hourly ouput
  integer :: IOUT            !frequency of daily ouput
  integer :: KOUT            !frequency of storage ouput
  integer :: IOLD            !last day of previous scenario
  integer :: ILAST           !last day of previous scenario
  integer :: IRUN            !start date of current scenario
  integer :: IBEGIN          !start date of model run
  integer :: ISTART          !start date of model run
  integer :: IEND            !end date of  current scenario
  integer :: LYRX            !last year
  integer :: LYRC            !number of days in current year
  integer :: LYRG            !num_of_simdays, defined for regression tests

  integer :: iYear_cur        !current year
  integer :: iyear_pre        !previous year
  integer :: iyear_rest       !restart year
  contains

  subroutine InitEcoSIMCtrlData

  implicit none

  allocate(ZEROS(JY,JX)); ZEROS(:,:) = 0._r8
  allocate(ZEROP(JP,JY,JX)); ZEROP(:,:,:)=0._r8
  allocate(ZEROQ(JP,JY,JX)); ZEROQ(:,:,:)=0._r8
  allocate(ZEROL(JP,JY,JX)); ZEROL(:,:,:)=0._r8
  allocate(ZEROS2(JY,JX));  ZEROS2(:,:)=0._r8

  end subroutine InitEcoSIMCtrlData

  subroutine DestructEcoSIMCtrlData
  use abortutils, only : destroy

  implicit none

  call destroy(ZEROS)
  call destroy(ZEROP)
  call destroy(ZEROQ)
  call destroy(ZEROL)
  call destroy(ZEROS2)

  end subroutine DestructEcoSIMCtrlData

end module EcoSIMCtrlDataType
