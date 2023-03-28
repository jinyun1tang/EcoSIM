module EcoSIMConfig
  
implicit none
  character(len=*),private, parameter :: mod_filename = __FILE__
  logical :: is_restart_run=.false.
  logical :: is_first_year=.false.
  logical :: transport_on=.true.
  logical :: column_mode=.false.
  logical :: do_instequil=.false.
  integer :: IFLGW                 !flag for raising Z0G with vegn

  integer, parameter :: ndbiomcpc = 2 !# of microbial residue components
  integer, parameter :: nlbiomcpc = 3 !# of living biomass components
  integer, parameter :: jskenc    = 4 !# of kinetic components of the substrates
  integer, parameter :: jcplxc    = 5 !# of microbe-substrate complexes
  integer, parameter :: jcplx1c   = jcplxc-1
  integer, parameter :: NFGsc     = 7 !# of microbial functional groups in each complex
  character(len=36)  :: case_name
  character(len=256) :: finidat =' '

  character(len=256), public :: hostname = ' '
  character(len=256), public :: username = ' '          ! Current user
  ! description of this source
  character(len=256), public :: source   = "ECOSIM"
  ! version of program
  character(len=256), public :: version  = " "
 
  integer, public :: nsrest  
  integer, public, parameter :: nsrBranch   = 2        ! Branch from restart files 
  integer, public, parameter :: nsrStartup  = 0        ! Startup from initial conditions
  integer, public, parameter :: nsrContinue = 1        ! Continue from restart files
  public :: is_restart,   &
            set_sim_type, &
            cold_run
contains

  logical function is_restart( )
  !
  ! Determine if restart run
  implicit none
  
  if (nsrest == nsrContinue) then
    is_restart = .true.
  else
    is_restart = .false.
  end if
  end function is_restart
!-----------------------------------------------------------------------

  subroutine set_sim_type()
  use EcoSIMCtrlMod, only : continue_run
  implicit none
  !determine simulation type

  if(continue_run)then
    nsrest=nsrContinue
  else
    if(finidat ==' ')then
      nsrest=nsrStartup
    else
      nsrest=nsrBranch
    endif
  endif
  end subroutine set_sim_type
!-----------------------------------------------------------------------

logical function cold_run()
implicit none

cold_run=(nsrest==nsrStartup)

end function cold_run
end module EcoSIMConfig
