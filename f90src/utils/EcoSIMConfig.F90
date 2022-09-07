module EcoSIMConfig

implicit none
  character(len=*),private, parameter :: mod_filename = __FILE__
  logical :: is_restart_run=.false.
  logical :: is_first_year=.false.
  integer :: IFLGW                 !flag for raising Z0G with vegn
end module EcoSIMConfig
