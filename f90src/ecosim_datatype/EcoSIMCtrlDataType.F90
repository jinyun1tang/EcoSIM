module EcoSIMCtrlDataType

! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__
  public
  real(r8) :: ZERO
  real(r8) :: ZERO2
  real(r8) :: ZEROS(JY,JX)
  real(r8) :: ZEROP(JP,JY,JX)
  real(r8) :: ZEROQ(JP,JY,JX)
  real(r8) :: ZEROL(JP,JY,JX)
  real(r8) :: ZEROS2(JY,JX)

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

  logical :: lverb           !logical switch for verbose output
  logical :: do_rgres        !logical switch for regression tests
end module EcoSIMCtrlDataType
